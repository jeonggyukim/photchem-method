// Build the GOW17 coefficient table in AthenaK (-DATHENAK) and in the Tigris copy;
// print every entry. The two outputs must be identical.
#include <cstdio>
#ifdef ATHENAK
#include "athena.hpp"
#include "chemistry/thermo/thermo_table.hpp"
using chemistry::ThermoTable; using chemistry::BuildThermoTable;
#else
#include "photchem/network/gow17_thermo_table.hpp"
using gow17::ThermoTable; using gow17::BuildThermoTable;
#endif
int main(int argc, char **argv) {
#ifdef ATHENAK
  Kokkos::initialize(argc, argv);
  {
#else
  (void)argc; (void)argv;
#endif
  ThermoTable tab;
  BuildThermoTable(tab);
  std::printf("%.17g %.17g\n", tab.nqt_t_min, tab.nqt_idt);
  for (int i = 0; i < ThermoTable::n_T; ++i) {
    for (int c = 0; c < ThermoTable::NCOEF; ++c) std::printf("%.17g ", tab.data(i, c));
    std::printf("\n");
  }
  const Real Ts[4] = {7.0, 63.0, 812.0, 9.5e3};
  for (Real T : Ts) {
    auto s = tab.Locate(T);
    std::printf("%d %d %.17g %.17g\n", s.i0, s.i1, s.w, tab.At(s, ThermoTable::ICII_k10e));
  }
#ifdef ATHENAK
  }
  Kokkos::finalize();
#endif
  return 0;
}
