#include <cstdio>
#include "photchem/network/gow17_thermo_table.hpp"
using namespace gow17;
int main() {
  ThermoTable g; BuildThermoTable(g);
  for (int i : {0, 750, 775, 1000, 1249, 1250, 1251, 1500, 2259})
    std::printf("row %d T %.6e expect %.6e\n", i, g.data(i, ThermoTable::ITEMP),
                std::pow(10.0, ThermoTable::logT_min + i*ThermoTable::dlogT));
  auto s = g.Locate(1e6); std::printf("Locate(1e6): i0 %d w %.3f\n", s.i0, s.w);
}
