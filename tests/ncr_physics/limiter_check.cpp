// At the equilibrium state of equil_one_cell.cpp: the timescales the substep
// controller uses (y/|C - D y| for H2 and CO, E/|Edot|), in code time, and the
// substep count of SolveODE for steps of 1e-3 to 1e1 code time.
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
int main() {
  const Real time_cgs = 3.0856776e18/1.0e5, edens_cgs = 1.6738234e-24*1.0e10;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = false;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  SemiImplicitSettings si;
  si.semi_implicit_cfl = 0.1; si.semi_implicit_yfloor = 1.0e-12;
  si.semi_implicit_n_substep_max = 100000; si.semi_implicit_n_substep_fixed = 0;
  si.semi_implicit_stiff_threshold = 1.0; si.semi_implicit_n_iter = 1;
  si.semi_implicit_energy = true; si.semi_implicit_renormalize = true;
  si.semi_implicit_gauss_seidel = true; si.semi_implicit_hep_first = true;
  si.semi_implicit_exact_map = false; si.semi_implicit_exact_block = false;
  si.semi_implicit_co_block = true; si.semi_implicit_exact_ghosts = false;
  si.semi_implicit_h2_first = false; si.semi_implicit_adaptive = true;
  si.semi_implicit_nbad_max = 3; si.semi_implicit_refresh_rates = false;
  si.semi_implicit_table_deriv = false;
  Real rad[8] = {1, 1, 1, 1, 0, 1, 1, 2e-16};
  const Real nH = 100.0, gamma = 5.0/3.0;
  Real y[13] = {0};
  y[7] = 0.01; y[4] = 1.6e-4; y[11] = 1.7e-6;
  y[12] = 20.0*Thermo::CvCold(0.0, 0.1, 0.01 + 1.6e-4 + 1.7e-6, gamma)*nH/edens_cgs;
  {  // bring to the equilibrium state first
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
    net.SetIonizedCooling(true, 1.0);
    for (int n = 0; n < 13; ++n) net.y(n) = y[n];
    SemiImplicit<GOW17Network> solver(si, net, 0.0, 1.0e5);
    solver.SolveODE();
    for (int n = 0; n < 13; ++n) y[n] = net.y(n);
  }
  GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
  net.SetIonizedCooling(true, 1.0);
  for (int n = 0; n < 13; ++n) net.y(n) = y[n];
  auto g = net.SetupNextStep(net.y);
  auto cd = net.CDRates(net.y, g);
  const char *nm[12] = {"He+","OHx","CHx","CO","C+","HCO+","H2","H+","H3+","H2+","O+","Si+"};
  std::printf("T %.2f K; timescales y/|C-Dy| [code time]:\n", net.Temperature(net.y, g));
  for (int n = 0; n < 12; ++n) {
    const Real f = cd.creation[n] - cd.destruction[n]*y[n];
    std::printf("  %-4s y %.3e  C %.3e  D %.3e  y/|f| %.3e  1/D %.3e\n", nm[n], y[n],
                cd.creation[n], cd.destruction[n], y[n]/std::fabs(f), 1.0/cd.destruction[n]);
  }
  const Real edot = net.Edot(net.y, g);
  std::printf("  E %.4e  Edot %.4e  E/|Edot| %.3e\n", y[12], edot, y[12]/std::fabs(edot));
  for (Real dt : {1e-3, 1e-2, 1e-1, 1.0, 10.0}) {
    GOW17Network n2(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
    n2.SetIonizedCooling(true, 1.0);
    for (int n = 0; n < 13; ++n) n2.y(n) = y[n];
    SemiImplicit<GOW17Network> solver(si, n2, 0.0, dt);
    solver.SolveODE();
    std::printf("dt %.0e code time: %d substeps\n", dt, solver.n_substeps);
  }
  return 0;
}
