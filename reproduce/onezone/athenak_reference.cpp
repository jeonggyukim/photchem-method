// Reference values from AthenaK's own GOW17Network (Kokkos Serial) for the Tigris
// port: ghosts, creation and destruction per species, Edot, on a grid of states.
// Build: see build_athenak_reference.sh beside this file.
#include <cstdio>
#include <limits>
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "chemistry/network/gow17.hpp"

using namespace chemistry;

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  {
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
    s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = inf;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1.0; s.length_cgs = 1.0; s.multi_d = false; s.three_d = false;

    DvceArray5D<Real> w0("w0", 1, 5, 1, 1, 3);
    DualArray1D<RegionSize> sizes("sizes", 1);
    sizes.h_view(0).dx1 = sizes.h_view(0).dx2 = sizes.h_view(0).dx3 = 1.0;
    sizes.template modify<HostMemSpace>(); sizes.template sync<DevMemSpace>();
    DvceArray1D<Real> ir("ir", 8);
    const Real rad[8] = {1.0, 0.5, 0.2, 0.8, 0.1, 0.9, 1.2, 2.0e-16};
    for (int n = 0; n < 8; ++n) ir(n) = rad[n];

    const Real nHs[3] = {1.0, 1.0e2, 1.0e4};
    const Real Ts[4] = {15.0, 100.0, 1.0e3, 8.0e3};
    // He+, OHx, CHx, CO, C+, HCO+, H2, H+, H3+, H2+, O+, Si+
    const Real ys[2][12] = {
      {1e-7, 1e-9, 1e-9, 1e-7, 1.5e-4, 1e-11, 1e-3, 2e-4, 1e-10, 1e-12, 1e-8, 1.5e-6},
      {1e-9, 1e-7, 1e-8, 1.3e-4, 5e-6, 1e-9, 0.49, 1e-7, 1e-8, 1e-12, 1e-11, 1e-7}};
    std::printf("# nH T_in iy | T ghosts(Si C O He e H) | C[12] | D[12] | Edot\n");
    for (Real nH : nHs) for (Real T : Ts) for (int iy = 0; iy < 2; ++iy) {
      for (int i = 0; i < 3; ++i) w0(0, IDN, 0, 0, i) = nH;
      GOW17Network net(s, 0, 0, 0, 1, w0, sizes, ir, 1.0, 1.0, 5.0/3.0, 1.0, 1.0, 1.0);
      Kokkos::View<Real[13], HostMemSpace> y("y");
      for (int n = 0; n < 12; ++n) y(n) = ys[iy][n];
      // energy density from the target temperature, as AthenaK's Temperature inverts
      Real xe = 0.0; for (int n : {0, 4, 5, 7, 8, 9, 10, 11}) xe += ys[iy][n];
      y(12) = T*Thermo::CvCold(ys[iy][6], 0.1, xe, 5.0/3.0)*nH;
      auto g = net.SetupNextStep(y);
      Real Tout = net.Temperature(y, g);
      auto cd = net.CDRates(y, g);
      Real edot = net.Edot(y, g);
      std::printf("%.17g %.17g %d | %.17g %.17g %.17g %.17g %.17g %.17g %.17g |",
                  nH, T, iy, Tout, g.Si, g.C, g.O, g.He, g.e, g.H);
      for (int n = 0; n < 12; ++n) std::printf(" %.17g", cd.creation[n]);
      std::printf(" |");
      for (int n = 0; n < 12; ++n) std::printf(" %.17g", cd.destruction[n]);
      std::printf(" | %.17g\n", edot);
    }
  }
  Kokkos::finalize();
  return 0;
}
