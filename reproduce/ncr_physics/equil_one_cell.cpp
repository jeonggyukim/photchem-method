// One cell of the hii problem generator's GOW17 equilibration: n_H = 100 cm^-3,
// 20 K, H+ 0.01, C+ and Si+ at their totals, chi0 = 1 on every channel but H2
// (h2_diss_bg_flag = false), CR 2e-16 s^-1, advanced by dt = 1e5 code time in one
// SolveODE call as AdvanceCells does. Tigris ism units. Reports substeps, wall time
// and the end state; then the same after splitting dt into pieces.
#include <chrono>
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
  Real rad[8] = {1, 1, 1, 1, 0, 1, 1, 2e-16};  // H2 channel 0, G_PE 1, CR
  const Real nH = 100.0, gamma = 5.0/3.0;
  for (int nsplit : {1, 100}) {
    Real y[13] = {0};
    y[7] = 0.01; y[4] = 1.6e-4; y[11] = 1.7e-6;  // H+, C+, Si+
    const Real xe = y[7] + y[4] + y[11];
    y[12] = 20.0*Thermo::CvCold(0.0, 0.1, xe, gamma)*nH/edens_cgs;
    long nsub = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int p = 0; p < nsplit; ++p) {
      GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
      net.SetIonizedCooling(true, 1.0);
      for (int n = 0; n < 13; ++n) net.y(n) = y[n];
      SemiImplicit<GOW17Network> solver(si, net, 0.0, 1.0e5/nsplit);
      solver.SolveODE();
      nsub += solver.n_substeps;
      for (int n = 0; n < 13; ++n) y[n] = net.y(n);
    }
    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
    for (int n = 0; n < 13; ++n) net.y(n) = y[n];
    auto g = net.SetupNextStep(net.y);
    std::printf("dt split in %3d: substeps %ld, %.2f ms, T %.2f K, x_H2 %.4g, x_H+ %.3g, "
                "x_CO %.3g, x_C+ %.3g, x_e %.3g\n", nsplit, nsub, ms,
                net.Temperature(net.y, g), y[6], y[7], y[3], y[4], g.e);
  }
  return 0;
}
