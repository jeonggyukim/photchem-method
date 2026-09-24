// Size of the ionized-hydrogen cooling added from NCR: Edot with and without it.
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
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = false;
  s.velocity_cgs = 1.0; s.length_cgs = 1.0; s.multi_d = false; s.three_d = false;
  Real rad[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 2e-16};
  // end state of the uniform run and an early, atomic state
  const Real ys[2][12] = {
    {4.47e-7, 1.43e-5, 6.1e-9, 1.335e-4, 3.5e-7, 6.05e-8, 0.4499, 2.94e-6, 1.26e-6,
     1.97e-9, 8.7e-11, 6.63e-7},
    {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6}};
  const Real Ts[2] = {45.0, 154.0};
  for (int iy = 0; iy < 2; ++iy) {
    for (int on = 0; on < 2; ++on) {
      GOW17Network net(s, 109.21, 0.0, View1D<Real>(rad, 8), 5.0/3.0, 1.0, 1.0);
      if (on) net.SetIonizedCooling(true, 1.0);
      Real xe = 0.0; for (int n : {0, 4, 5, 7, 8, 9, 10, 11}) xe += ys[iy][n];
      for (int n = 0; n < 12; ++n) net.y(n) = ys[iy][n];
      net.y(12) = Ts[iy]*Thermo::CvCold(ys[iy][6], 0.1, xe, 5.0/3.0)*109.21;
      auto g = net.SetupNextStep(net.y);
      std::printf("state %d T %.1f ionized_cooling %d: Edot %.6e\n", iy, Ts[iy], on,
                  net.Edot(net.y, g));
    }
  }
  return 0;
}
