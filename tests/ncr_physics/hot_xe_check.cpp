// Hot-gas x_e with tracked O/S ions: species electrons + CIE extras should equal
// H+ + the table's He and metal electrons (no double counting of tracked ions).
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "../paths.hpp"
using namespace gow17;
int main() {
  ThermoTable grid; BuildThermoTable(grid);
  HotCIETable hot; hot.Load(TigrisDir() + "inputs/tables/tigress_coolftn_ncr.txt");
  GOW17Settings s; s.use_thermo_table = true; s.thermo_table = grid; s.isothermal = false;
  s.isothermal_temperature = 0; s.zd = 0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1; const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf; s.temperature_min_cooling = 1;
  s.temperature_max_cooling_nm = 1e9; s.Leff_CO_max = 3e20; s.H2_rovib_cooling = true;
  s.is_kgrH2_const = false; s.velocity_cgs = 1; s.length_cgs = 1; s.multi_d = false; s.three_d = false;
  Real rad[GOW17Network::n_freq] = {};
  const double T = 5e5, nH = 1.0;
  GOW17Network net(s, nH, 0.0, View1D<Real>(rad, GOW17Network::n_freq), 5.0/3.0, 1.0, 1.0);
  net.SetHotCooling(&hot, 2e4, 3.5e4, 1.0, false);
  for (int n = 0; n < GOW17Network::neqs; ++n) net.y(n) = 0.0;
  net.y(IH_plus) = 1.0; net.y(IHe_plus) = 0.1; net.y(IC_plus) = 1.6e-4; net.y(IO_2plus) = 3.2e-4;
  net.y(ISi_plus) = 1.7e-6; net.y(IS_3plus) = 1.45e-5;
  const double xe0 = 1.0 + 0.1 + 1.6e-4 + 2*3.2e-4 + 1.7e-6 + 3*1.45e-5;
  net.y(GOW17Network::neqs - 1) = T*Thermo::CvCold(0.0, 0.1, xe0, 5.0/3.0)*nH;
  const auto g = net.SetupNextStep(net.y);
  std::printf("x_e network %.6f | expected H+ + table(He) + table(metal) = %.6f\n", g.e,
              1.0 + hot.At(HotCIETable::IXE_HE, T) + hot.At(HotCIETable::IXE_METAL, T));
}
