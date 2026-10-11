// CIEElement pool sums at the CHIANTI grid points against pyathena's
// data/chemistry/cie_high_pool_<El>.txt (build_cie_high_pool.py).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "photchem/network/gow17_cie_tables.hpp"
#include "../paths.hpp"
int main() {
  const std::string dir = TigrisDir() + "inputs/tables/chianti_v11/";
  struct E { const char *name; int z; int q_max; };
  for (E e : {E{"C", 6, 1}, E{"N", 7, 1}, E{"O", 8, 2}, E{"S", 16, 2}}) {
    gow17::CIEElement el(dir, e.name, e.z);
    std::ifstream f(PyathenaDir() + "data/chemistry/cie_high_pool_" +
                    e.name + ".txt");
    std::string line; int i = 0; double worst[3] = {0, 0, 0};
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line); double lt, xr, qr, lr; ss >> lt >> xr >> qr >> lr;
      Real x, q, l; el.Pool(i, e.q_max, &x, &q, &l);
      if (std::fabs(el.log_t(i) - lt) > 1e-6) { std::printf("grid mismatch\n"); return 1; }
      const double ref[3] = {xr, qr, lr}, got[3] = {x, q, l};
      for (int c = 0; c < 3; ++c)
        if (ref[c] > 1e-30) worst[c] = std::fmax(worst[c], std::fabs(got[c]/ref[c] - 1));
      ++i;
    }
    std::printf("%s q_max=%d: %d points, max rel diff x_high %.1e q_mean %.1e Lambda_high %.1e\n",
                e.name, e.q_max, i, worst[0], worst[1], worst[2]);
  }
  return 0;
}
