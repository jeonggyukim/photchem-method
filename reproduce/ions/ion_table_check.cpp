// Ion-rate table (gow17_ion_table.hpp) against the exact evaluators of photchem/rates
// for O0-O++ and S0-S++, at 4000 log-uniform temperatures in 10 K - 1e9 K.
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include "photchem/network/gow17_ion_table.hpp"
using namespace gow17;
int main() {
  const std::string d = std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/rates/";
  Rates::RecombRate rec; rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat",
                                  d + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci; ci.Load(d + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid; BuildThermoTable(grid);
  const std::vector<IonStage> st = {{8, 0}, {8, 1}, {8, 2}, {16, 0}, {16, 1}, {16, 2}};
  IonRateTable tab; BuildIonRateTable(grid, st, rec, ci, ct, &tab);
  const char *cname[4] = {"rec", "coll.ion", "CT rec", "CT ion"};
  std::mt19937 g(3); std::uniform_real_distribution<double> u(1.0, 9.0);
  for (int k = 0; k < (int)st.size(); ++k) {
    for (int c = 0; c < IonRateTable::NCOL; ++c) {
      std::vector<double> e; double emax = 0, Tmax = 0; int nz = 0;
      for (int t = 0; t < 4000; ++t) {
        const double T = std::pow(10.0, u(g));
        const int z = st[k].z, q = st[k].q, n = z - q;
        double ex = 0;
        if (c == 0) ex = (q >= 1) ? rec.Rate(z, n, T) : 0.0;
        if (c == 1) ex = (n >= 1) ? ci.Rate(z, n, T) : 0.0;
        if (c == 2) ex = (q >= 1) ? ct.CtRec(z, n, T) : 0.0;
        if (c == 3) ex = (n >= 1) ? ct.CtIon(z, n, T) : 0.0;
        const double tb = tab.At(grid.Locate(T), k, c);
        if (ex < 1e-30) { if (tb > 1e-30) ++nz; continue; }
        const double r = std::fabs(tb/ex - 1); e.push_back(r);
        if (r > emax) { emax = r; Tmax = T; }
      }
      if (e.empty()) continue;
      std::sort(e.begin(), e.end());
      std::printf("Z=%2d q=%d %-8s n=%4zu median %.1e  99th %.1e  max %.1e at T=%.3g K%s\n",
                  st[k].z, st[k].q, cname[c], e.size(), e[e.size()/2], e[e.size()*99/100],
                  emax, Tmax, nz ? "  (nonzero where exact is 0!)" : "");
    }
  }
  return 0;
}
