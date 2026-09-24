// Same states and calls as athenak_reference.cpp, on the Tigris copy of the network
// (gow17::GOW17Network). The two outputs must agree to round-off.
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;

int main() {
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = inf;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1.0; s.length_cgs = 1.0; s.multi_d = false; s.three_d = false;

  Real rad[8] = {1.0, 0.5, 0.2, 0.8, 0.1, 0.9, 1.2, 2.0e-16};
  const Real nHs[3] = {1.0, 1.0e2, 1.0e4};
  const Real Ts[4] = {15.0, 100.0, 1.0e3, 8.0e3};
  const Real ys[2][12] = {
    {1e-7, 1e-9, 1e-9, 1e-7, 1.5e-4, 1e-11, 1e-3, 2e-4, 1e-10, 1e-12, 1e-8, 1.5e-6},
    {1e-9, 1e-7, 1e-8, 1.3e-4, 5e-6, 1e-9, 0.49, 1e-7, 1e-8, 1e-12, 1e-11, 1e-7}};
  std::printf("# nH T_in iy | T ghosts(Si C O He e H) | C[12] | D[12] | Edot\n");
  for (Real nH : nHs) for (Real T : Ts) for (int iy = 0; iy < 2; ++iy) {
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), 5.0/3.0, 1.0, 1.0);
    Real ybuf[13];
    View1D<Real> y(ybuf, 13);
    for (int n = 0; n < 12; ++n) y(n) = ys[iy][n];
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
  return 0;
}
