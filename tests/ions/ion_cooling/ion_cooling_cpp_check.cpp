// C++ IonCoolingTable (resampled onto the ThermoTable rows) at (T, n_e) points
// read from stdin, one "ion T ne" per line; prints Lambda for each.
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "photchem/network/gow17_ion_cooling.hpp"
#include "../../paths.hpp"
using namespace gow17;
int main() {
  ThermoTable grid; BuildThermoTable(grid);
  const std::vector<std::string> names = {"o_2", "o_3", "s_2", "s_3", "s_4",
                                          "n_2", "n_3", "ne_2", "ne_3"};
  IonCoolingTable tab;
  BuildIonCoolingTable(grid, TigrisDir() + "inputs/tables/chianti_v11/ion_cooling/", names,
                       &tab);
  std::string ion; double T, ne;
  while (std::cin >> ion >> T >> ne) {
    int k = 0; while (names[k] != ion) ++k;
    std::printf("%s %.6e %.6e %.6e\n", ion.c_str(), T, ne, tab.Lambda(grid.Locate(T), k, ne));
  }
}
