// The pdr_slab GOW17 cell with the largest heating/cooling imbalance, integrated
// alone in its own fields: how far a capped postproc iteration takes it, and
// whether H/C reaches 1.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
int main() {
  std::ifstream f("worst_cell.txt");
  Real rho, press, y[12], chi_ci, chi_lw, chi_h2, chi_pe, xi_cr;
  f >> rho >> press; for (Real &v : y) f >> v; f >> chi_ci >> chi_lw >> chi_h2 >> chi_pe >> xi_cr;
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = rho/1.4;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  SemiImplicitSettings si{0.1, 1e-12, 1000000000, 0, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  Real rad[8] = {chi_ci, chi_lw, chi_lw, chi_lw, chi_h2, chi_lw, chi_pe, xi_cr};
  GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  for (int n = 0; n < 12; ++n) net.y(n) = y[n];
  net.y(12) = press/(gamma - 1.0);
  Real t = 0.0;
  for (Real dt : {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}) {
    SemiImplicit<GOW17Network> solver(si, net, 0.0, dt);
    solver.SolveODE(); t += dt;
    auto g = net.SetupNextStep(net.y);
    const Real T = net.Temperature(net.y, g);
    std::printf("t %8.1e code (+%8d substeps, mean %.2e yr)  T %.3f  H/C %.4f  x_H2 %.4e x_CO %.3e\n",
                t, solver.n_substeps, dt*0.978e6/solver.n_substeps, T,
                net.HeatingTerm(net.y, g, T)/net.CoolingTerm(net.y, g, T), net.y(6), net.y(3));
  }
}
