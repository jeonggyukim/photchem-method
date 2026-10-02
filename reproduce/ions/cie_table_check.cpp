// CIEPoolTable on the ThermoTable rows: (1) read through Locate at the CHIANTI T, against
// pyathena's cie_high_pool_<El>.txt; (2) linear-in-logT interpolation error, from the
// error of 0.1-dex interpolation at the skipped 0.05-dex points (the 0.05-dex error
// is about a quarter of it); (3) tracked stages at CIE plus the pool reproduce the
// element's full CIE cooling; (4) return-rate columns at 1e4 K and 1e5 K.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "photchem/network/gow17_cie_tables.hpp"
#include "../paths.hpp"
using namespace gow17;
int main() {
  const std::string tdir = TigrisDir() + "inputs/tables/";
  const std::string cdir = tdir + "chianti_v11/", rdir = tdir + "rates/";
  Rates::RecombRate rec; rec.Load(rdir + "badnell_rr_2023.dat", rdir + "badnell_dr_C_2023.dat",
                                  rdir + "badnell_dr_E_2023.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(rdir + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(rdir + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid; BuildThermoTable(grid);
  const std::vector<CIEElementSpec> els = {{"O", 8, 2, 0.0}, {"S", 16, 2, 0.0},
                                           {"C", 6, 1, 0.0}, {"N", 7, 1, 0.0},
                                           {"Ne", 10, -1, 1.0}, {"Fe", 26, -1, 2.0}};
  CIEPoolTable tab; BuildCIEPoolTable(grid, cdir, els, rec, ct, &tab);

  // (1)
  for (int p = 0; p < 4; ++p) {
    std::ifstream f(PyathenaDir() + "data/chemistry/cie_high_pool_" +
                    els[p].name + ".txt");
    std::string line; double worst[3] = {0, 0, 0}; int nmatch = 0;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line); double lt, xr, qr, lr; ss >> lt >> xr >> qr >> lr;
      const ThermoTable::Slot sl = grid.Locate(std::pow(10.0, lt));
      ++nmatch;
      const double got[3] = {tab.At(sl, p, CIEPoolTable::IX), tab.At(sl, p, CIEPoolTable::IQ),
                             tab.At(sl, p, CIEPoolTable::IX)*tab.At(sl, p, CIEPoolTable::ILAMBDA)};
      const double ref[3] = {xr, qr, lr};
      for (int c = 0; c < 3; ++c)
        if (xr > 1e-6) worst[c] = std::fmax(worst[c], std::fabs(got[c]/ref[c] - 1));
    }
    std::printf("(1) %s q_max=%d: %d CHIANTI points, max rel diff where x_high > 1e-6: x_high %.1e q_mean %.1e "
                "x_high*Lambda %.1e\n", els[p].name.c_str(), els[p].q_max, nmatch,
                worst[0], worst[1], worst[2]);
  }

  // (2) and (3)
  for (const CIEElementSpec &e : els) {
    CIEElement el(cdir, e.name, e.z);
    double worst_l = 0, worst_x = 0, worst_sum = 0; double t_l = 0;
    for (int i = 1; i + 1 < el.npoint(); i += 2) {
      Real x0, q0, l0, x1, q1, l1, xm, qm, lm;
      const int qm_ = (e.q_max < 0) ? -1 : e.q_max;
      el.Pool(i - 1, -1, &x0, &q0, &l0); el.Pool(i + 1, -1, &x1, &q1, &l1);
      el.Pool(i, -1, &xm, &qm, &lm);
      const double le = std::fabs(0.5*(l0 + l1)/lm - 1);
      if (le > worst_l) { worst_l = le; t_l = el.log_t(i); }
      if (qm_ >= 0) {
        el.Pool(i - 1, qm_, &x0, &q0, &l0); el.Pool(i + 1, qm_, &x1, &q1, &l1);
        el.Pool(i, qm_, &xm, &qm, &lm);
        worst_x = std::fmax(worst_x, std::fabs(0.5*(x0 + x1) - xm));
      }
    }
    // (3) at every thermo row from 1e3 K up
    std::vector<Real> x, lam;
    int p = 0; for (const CIEElementSpec &e2 : els) { if (&e2 == &e) break; if (e2.q_max >= 0) ++p; }
    for (int i = 0; i < ThermoTable::n_T && e.q_max >= 0; ++i) {
      const Real T = grid.data(i, ThermoTable::ITEMP);
      if (T < 1e3) continue;
      el.Stages(T, &x, &lam);
      double full = 0, low = 0;
      for (int q = 0; q <= e.z; ++q) { full += x[q]*lam[q]; if (q <= e.q_max) low += x[q]*lam[q]; }
      const int c0 = p*CIEPoolTable::NCOL;
      const double pool = tab.data(i, c0 + CIEPoolTable::IX)*tab.data(i, c0 + CIEPoolTable::ILAMBDA);
      if (full > 1e-35) worst_sum = std::fmax(worst_sum, std::fabs((low + pool)/full - 1));
    }
    std::printf("(2) %-2s full-CIE cooling, 0.1-dex interp error max %.1e (logT %.2f); "
                "pool fraction abs error max %.1e\n", e.name.c_str(), worst_l, t_l, worst_x);
    if (e.q_max >= 0) std::printf("(3) %-2s tracked + pool vs full CIE cooling: max rel diff %.1e\n",
                                  e.name.c_str(), worst_sum);
  }

  // (4)
  for (Real T : {1e4, 1e5, 1e6}) {
    const ThermoTable::Slot s = grid.Locate(T);
    for (int p = 0; p < 4; ++p)
      std::printf("(4) T %.0e %s: x_pool %.3e q_mean %.2f Lambda/ion %.3e alpha_ret %.3e "
                  "kCT_ret %.3e\n", T, els[p].name.c_str(), tab.At(s, p, CIEPoolTable::IX),
                  tab.At(s, p, CIEPoolTable::IQ), tab.At(s, p, CIEPoolTable::ILAMBDA),
                  tab.At(s, p, CIEPoolTable::IREC), tab.At(s, p, CIEPoolTable::ICT_REC));
    std::printf("(4) T %.0e untracked Ne + 2 Fe: %.3e\n", T, tab.Untracked(s));
  }
  return 0;
}
