// GOW17 cooling of fully ionized gas (x_H+ = 1, He+ = 0.1, C+ O+ Si+ at their
// totals) at n_H = 1 cm^-3 from 1e4 to 1e8 K, per n_H^2, against the NCR CIE
// table total Lambda_H + Lambda_He + Lambda_metal (tigress_coolftn_ncr.txt).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;
int main() {
  std::ifstream f(std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt");
  std::string l; std::vector<std::vector<double>> tab;
  while (std::getline(f, l)) { if (l[0] == '#') continue; std::istringstream s(l); std::vector<double> r; double v; while (s >> v) r.push_back(v); if (r.size() == 10) tab.push_back(r); }
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
  HotCIETable hot; hot.Load(std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt");
  const bool use_hot = std::getenv("HOT") != nullptr;
  Real rad[8] = {1, 1, 1, 1, 1, 1, 1, 2e-16};
  std::printf("%8s %12s %12s %12s %8s\n", "logT", "GOW17 C/nH^2", "NCR CIE", "ratio", "heat/cool");
  for (const auto &r : tab) {
    if (std::fmod(r[0] + 1e-9, 0.5) > 0.02 && std::fabs(r[0] - 4.4) > 0.01) continue;
    const Real T = std::pow(10.0, r[0]);
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
    net.SetIonizedCooling(true, 1.0);
    if (use_hot) net.SetHotCooling(&hot, 2.0e4, 3.5e4, 1.0);
    for (int n = 0; n < 13; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0; net.y(IHe_plus) = 0.1; net.y(IC_plus) = 1.6e-4; net.y(IO_plus) = 3.2e-4;
    net.y(ISi_plus) = 1.7e-6;
    auto g0 = net.SetupNextStep(net.y);
    net.y(12) = T*Thermo::CvCold(0.0, 0.1, g0.e, gamma)*nH/edens;
    auto g = net.SetupNextStep(net.y);
    const Real Tn = net.Temperature(net.y, g);
    const Real c = net.CoolingTerm(net.y, g, Tn), h = net.HeatingTerm(net.y, g, Tn);
    const Real ncr = (r[1] + r[2] + r[3]);  // per n_H^2? see header
    std::printf("%8.2f %12.3e %12.3e %12.3e %8.2e  (T %.3g)\n", r[0], c/nH, ncr, (c/nH)/ncr, h/c, Tn);
  }
}
