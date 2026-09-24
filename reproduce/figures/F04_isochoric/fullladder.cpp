// F04 reference: every stage of C, N, O, Si, S integrated along the T(t), n_e(t), x_H0(t)
// and x_H+(t) of an isochoric run (isochoric.cpp output), with the same rate files
// (Badnell RR+DR, Voronov CI, Kingdon & Ferland CT) the network uses. Each element's
// ladder takes backward-Euler (tridiagonal) steps, 400 per output interval, from CIE
// at the run's first temperature. So the comparison isolates the one approximation of
// the network: the stages above the top tracked one lumped with a CIE split.
// usage: fullladder isochoric_n1.txt N_H > fullladder_n1.txt
// Output per row: t, T, then per element (C, Si, O, S, N) the fractions of every stage.
// Build: OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2
//   -I$HOME/Projects/tigris-gow17/src -I/opt/homebrew/opt/hdf5-mpi/include fullladder.cpp
//   $HOME/Projects/tigris-gow17/src/photchem/network/gow17_thermo_table.cpp -o fullladder
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "photchem/network/gow17_cie_tables.hpp"
#include "photchem/rates/coll_ion_rate.hpp"
#include "photchem/rates/ct_rate.hpp"
#include "photchem/rates/recomb_rate.hpp"

int main(int argc, char **argv) {
  const std::string home = getenv("HOME");
  const std::string rd = home + "/Projects/tigris-gow17/inputs/tables/rates/";
  const std::string cd = home + "/Projects/tigris-gow17/inputs/tables/chianti_v11/";
  const Real nH = std::atof(argv[2]);
  Rates::RecombRate rec;
  rec.Load(rd + "badnell_rr_2023.dat", rd + "badnell_dr_C_2023.dat", rd + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci;
  ci.Load(rd + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct;
  ct.LoadRecomb(rd + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(rd + "kingdon_ferland96_ct_ion.dat");

  // the isochoric run: t, T, x_e, then species; H2 at column 9, H+ at column 10
  std::vector<std::vector<Real>> run;
  std::ifstream f(argv[1]);
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::vector<Real> r;
    Real v;
    while (ss >> v) r.push_back(v);
    run.push_back(r);
  }
  struct El { const char *name; int z; };
  const El els[5] = {{"C", 6}, {"Si", 14}, {"O", 8}, {"S", 16}, {"N", 7}};
  std::vector<std::vector<Real>> x(5);
  for (int k = 0; k < 5; ++k) {
    std::vector<Real> lq;
    gow17::CIEElement(cd, els[k].name, els[k].z).Stages(run[0][1], &x[k], &lq);
  }
  auto state = [&](const std::vector<Real> &r, Real *T, Real *ne, Real *xh0, Real *xhp) {
    *T = r[1];
    *ne = nH*r[2];
    *xhp = r[10];
    *xh0 = std::fmax(1.0 - r[10] - 2.0*r[9], 0.0);
  };
  std::printf("# full-ladder reference along %s, n_H = %g; t T then per element C Si O S N "
              "fractions of stages 0..Z\n", argv[1], nH);
  auto print = [&](Real t, Real T) {
    std::printf("%.6e %.6e", t, T);
    for (int k = 0; k < 5; ++k) for (Real v : x[k]) std::printf(" %.6e", v);
    std::printf("\n");
  };
  print(run[0][0], run[0][1]);
  const int nsub = 400;
  for (std::size_t i = 1; i < run.size(); ++i) {
    const Real t0 = run[i - 1][0], t1 = run[i][0];
    const Real h = (t1 - t0)/nsub*3.0856776e18/1.0e5;  // code time -> s
    for (int s = 0; s < nsub; ++s) {
      // state at the end of the substep, linear in time between the two rows
      const Real w = (s + 1.0)/nsub;
      std::vector<Real> r(run[i].size());
      for (std::size_t c = 0; c < r.size(); ++c) r[c] = (1 - w)*run[i - 1][c] + w*run[i][c];
      Real T, ne, xh0, xhp;
      state(r, &T, &ne, &xh0, &xhp);
      for (int k = 0; k < 5; ++k) {
        const int z = els[k].z;
        std::vector<Real> up(z + 1, 0.0), dn(z + 1, 0.0);  // q -> q+1, q -> q-1 [s^-1]
        for (int q = 0; q <= z; ++q) {
          const int nel = z - q;
          if (q < z) {
            up[q] = ci.Rate(z, nel, T)*ne;
            if (ct.HasIon(z, nel)) up[q] += ct.CtIon(z, nel, T)*nH*xhp;
          }
          if (q > 0) {
            dn[q] = rec.Rate(z, nel, T)*ne;
            if (ct.HasRec(z, nel)) dn[q] += ct.CtRec(z, nel, T)*nH*xh0;
          }
        }
        // (1 + h (up_q + dn_q)) x_q - h up_{q-1} x_{q-1} - h dn_{q+1} x_{q+1} = x_q^old
        std::vector<Real> a(z + 1), b(z + 1), c(z + 1), d = x[k];
        for (int q = 0; q <= z; ++q) {
          a[q] = (q > 0) ? -h*up[q - 1] : 0.0;
          b[q] = 1.0 + h*(up[q] + dn[q]);
          c[q] = (q < z) ? -h*dn[q + 1] : 0.0;
        }
        for (int q = 1; q <= z; ++q) {
          const Real m = a[q]/b[q - 1];
          b[q] -= m*c[q - 1];
          d[q] -= m*d[q - 1];
        }
        x[k][z] = d[z]/b[z];
        for (int q = z - 1; q >= 0; --q) x[k][q] = (d[q] - c[q]*x[k][q + 1])/b[q];
        Real sum = 0.0;
        for (Real v : x[k]) sum += v;
        for (Real &v : x[k]) v /= sum;
      }
    }
    print(t1, run[i][1]);
  }
}
