// Hot-gas electrons: E set for temperature T with NCR's x_e = x_H+ + xe_He(T) +
// xe_metal(T); the network must read back that T and that x_e (fully ionized H,
// He+ and the tracked metal ions at zero, so the table supplies all of theirs).
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
  std::printf("%8s %10s %10s %10s %10s\n", "T", "T_back", "x_e", "x_e NCR", "rel");
  for (Real T : {1.5e4, 2.5e4, 3.0e4, 5e4, 1e5, 1e6, 1e7, 1e8}) {
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
    net.SetHotCooling(&hot, 2.0e4, 3.5e4, 1.0);
    for (int n = 0; n < 13; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0;
    const Real xe_ncr = 1.0 + hot.At(HotCIETable::IXE_HE, T) + hot.At(HotCIETable::IXE_METAL, T);
    net.y(12) = T*Thermo::CvCold(0.0, 0.1, xe_ncr, gamma)*nH/edens;
    auto g = net.SetupNextStep(net.y);
    std::printf("%8.2e %10.4e %10.5f %10.5f %10.2e\n", T, net.Temperature(net.y, g), g.e,
                xe_ncr, g.e/xe_ncr - 1);
  }
}
