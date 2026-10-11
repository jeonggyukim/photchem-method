// Per-stage interpolation of CHIANTI x_q and L_q from 0.1-dex spacing onto the skipped
// 0.05-dex points: linear vs logarithmic (floor 1e-300, fractions renormalized to sum 1).
// Errors weighted where they matter: full-element cooling relative error where the
// cooling is above 1e-3 of its peak over T >= 1e4 K... reported per T band.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "photchem/network/gow17_cie_tables.hpp"
#include "../paths.hpp"
using namespace gow17;
int main() {
  const std::string cdir = TigrisDir() + "inputs/tables/chianti_v11/";
  struct E { const char *n; int z; int qmax; };
  for (E e : {E{"H", 1, 0}, E{"He", 2, 0}, E{"C", 6, 1}, E{"N", 7, 1}, E{"O", 8, 2},
              E{"Ne", 10, 1}, E{"Si", 14, 1}, E{"S", 16, 2}, E{"Fe", 26, 1}}) {
    CIEElement el(cdir, e.n, e.z);
    std::vector<Real> xa, la, xb, lb, xm, lm;
    double peak = 0;
    for (int i = 0; i < el.npoint(); ++i) {
      el.Stages(std::pow(10.0, el.log_t(i)), &xm, &lm);
      double s = 0; for (int q = 0; q <= e.z; ++q) s += xm[q]*lm[q]; peak = std::max(peak, s);
    }
    double wl[2] = {0, 0}, wp[2] = {0, 0}; double tl[2] = {0, 0};
    for (int i = 1; i + 1 < el.npoint(); i += 2) {
      el.Stages(std::pow(10.0, el.log_t(i - 1)), &xa, &la);
      el.Stages(std::pow(10.0, el.log_t(i + 1)), &xb, &lb);
      el.Stages(std::pow(10.0, el.log_t(i)), &xm, &lm);
      double truth = 0, tpool = 0;
      for (int q = 0; q <= e.z; ++q) { truth += xm[q]*lm[q]; if (q > e.qmax) tpool += xm[q]; }
      for (int s = 0; s < 2; ++s) {
        std::vector<Real> x(e.z + 1), l(e.z + 1); double sx = 0;
        for (int q = 0; q <= e.z; ++q) {
          if (s == 0) { x[q] = 0.5*(xa[q] + xb[q]); l[q] = 0.5*(la[q] + lb[q]); }
          else {
            x[q] = std::sqrt(std::max(xa[q], 1e-300)*std::max(xb[q], 1e-300));
            l[q] = std::sqrt(std::max(la[q], 1e-300)*std::max(lb[q], 1e-300));
          }
          sx += x[q];
        }
        double c = 0, pool = 0;
        for (int q = 0; q <= e.z; ++q) { if (s == 1) x[q] /= sx; c += x[q]*l[q]; if (q > e.qmax) pool += x[q]; }
        if (truth > 1e-3*peak) {
          const double r = std::fabs(c/truth - 1);
          if (r > wl[s]) { wl[s] = r; tl[s] = el.log_t(i); }
        }
        wp[s] = std::max(wp[s], std::fabs(pool - tpool));
      }
    }
    std::printf("%-2s cooling rel err (where > 1e-3 peak): linear %.1e (logT %.2f) log %.1e (logT %.2f); "
                "pool frac abs err: linear %.1e log %.1e\n", e.n, wl[0], tl[0], wl[1], tl[1], wp[0], wp[1]);
  }
  return 0;
}
