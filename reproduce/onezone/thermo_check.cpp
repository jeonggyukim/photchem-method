// Same calls into AthenaK chemistry::Thermo (-DATHENAK) and the Tigris copy
// gow17::Thermo; the two outputs must be identical.
#include <cstdio>
#ifdef ATHENAK
#include "athena.hpp"
#include "chemistry/thermo/thermo.hpp"
using chemistry::Thermo;
#else
#include "photchem/network/gow17_thermo.hpp"
using gow17::Thermo;
#endif
int main(int argc, char **argv) {
#ifdef ATHENAK
  Kokkos::initialize(argc, argv);
#endif
  const Real Ts[6] = {5.0, 20.0, 100.0, 1.0e3, 5.0e3, 1.0e4};
  const Real nHs[3] = {1.0, 1.0e2, 1.0e4};
  for (Real nH : nHs) for (Real T : Ts) {
    const Real nHI = 0.4*nH, nH2 = 0.3*nH, ne = 1e-4*nH;
    std::printf("%g %g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", nH, T,
                Thermo::CvCold(0.3, 0.1, 1e-4, 5.0/3.0),
                Thermo::HeatingPE(1.2, 1.0, T, ne),
                Thermo::CoolingCII(1e-4, nHI, nH2, ne, T),
                Thermo::CoolingOI(3e-4, nHI, nH2, ne, T),
                Thermo::CoolingLya(0.4, ne, T),
                Thermo::CoolingCOR(1e-4, nHI, nH2, ne, T, 1e15),
                Thermo::CoolingRec(1.0, T, ne, 1.2));
  }
#ifdef ATHENAK
  Kokkos::finalize();
#endif
  return 0;
}
