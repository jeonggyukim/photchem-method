// Grain-assisted recombination in hot gas: D(H+), D(C+) with the hot table set
// against without, fully ionized H at n_H = 1, chi = 1.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;
int main() {
  HotCIETable hot; hot.Load(std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt");
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10, nH = 1.0;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  Real rad[8] = {1, 1, 1, 1, 1, 1, 1, 2e-16};
  std::printf("%8s %12s %12s %12s\n", "T", "grainH+ ratio", "grainC+ ratio", "w1 expected");
  for (Real T : {1.5e4, 2.2e4, 2.75e4, 3.2e4, 5e4, 1e6}) {
    Real d[2][2], ee[2];
    for (int hotp = 0; hotp < 2; ++hotp) {
     Real dz[2][2];
     for (int zdust = 0; zdust < 2; ++zdust) {
      GOW17Settings sz = s; sz.zd = zdust ? 1.0 : 0.0;
      GOW17Network net(sz, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
      if (hotp) net.SetHotCooling(&hot, 2.0e4, 3.5e4, 1.0);
      for (int n = 0; n < 13; ++n) net.y(n) = 0.0;
      net.y(IH_plus) = 1.0; net.y(IC_plus) = 1.6e-4; net.y(IHe_plus) = 0.1;
      auto g0 = net.SetupNextStep(net.y);
      net.y(12) = T*Thermo::CvCold(0.0, 0.1, g0.e, gamma)*nH/edens;
      auto g = net.SetupNextStep(net.y);
      // same x_e in both, so the only difference is the grain term
      const auto r = net.CDRates(net.y, g);
      dz[zdust][0] = r.destruction(IH_plus); dz[zdust][1] = r.destruction(IC_plus); ee[hotp] = g.e;
     }
     d[hotp][0] = dz[1][0] - dz[0][0]; d[hotp][1] = dz[1][1] - dz[0][1];
    }
    const Real w2 = T >= 3.5e4 ? 1.0 : 1.0/(1.0 + std::exp(-10.0*(T - 2.75e4)/1.5e4));
    std::printf("%8.2e %12.4f %12.4f %12.4f   x_e %.4f/%.4f\n", T, d[1][0]/d[0][0], d[1][1]/d[0][1], T >= 2e4 ? 1 - w2 : 1.0, ee[1], ee[0]);
  }
}
