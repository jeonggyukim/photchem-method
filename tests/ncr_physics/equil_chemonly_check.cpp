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
  auto state = [&](const Real *yy, Real *T, Real *heat, Real *cool) {
    GOW17Network n2(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
    n2.SetIonizedCooling(true, 1.0);
    for (int n = 0; n < 13; ++n) n2.y(n) = yy[n];
    auto g = n2.SetupNextStep(n2.y);
    *T = n2.Temperature(n2.y, g);
    *heat = n2.HeatingTerm(n2.y, g, *T);
    *cool = n2.CoolingTerm(n2.y, g, *T);
  };
  // (a) successive adaptive substeps from the stuck state: one call of dt = 1e-3
  // code time takes 3 substeps; trace 8 such calls
  Real ya[13]; for (int n = 0; n < 13; ++n) ya[n] = y[n];
  std::printf("(a) trace, calls of dt = 1e-3 code time:\n");
  for (int c = 0; c < 8; ++c) {
    Real T, h, cl; state(ya, &T, &h, &cl);
    std::printf("  call %d: T %.4f K  heat %.4e cool %.4e  x_H+ %.4e x_CO %.4e\n", c, T, h, cl,
                ya[7], ya[3]);
    GOW17Network n2(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
    n2.SetIonizedCooling(true, 1.0);
    for (int n = 0; n < 13; ++n) n2.y(n) = ya[n];
    SemiImplicit<GOW17Network> solver(si, n2, 0.0, 1e-3);
    solver.SolveODE();
    for (int n = 0; n < 13; ++n) ya[n] = n2.y(n);
  }
  // (b) NCR-style: repeat one backward-Euler step of the full dt (n_substep_fixed = 1),
  // from the pgen's initial state, until x_H, x_H+ change < 1e-6 and |H-C|/(H+C) < 1e-6
  SemiImplicitSettings s1 = si; s1.semi_implicit_n_substep_fixed = 1; s1.semi_implicit_energy = false;
  for (Real dt : {1e5, 1e2, 1.0}) {
    Real yb[13] = {0};
    yb[7] = 0.01; yb[4] = 1.6e-4; yb[11] = 1.7e-6;
    yb[12] = 20.0*Thermo::CvCold(0.0, 0.1, 0.01 + 1.6e-4 + 1.7e-6, gamma)*nH/edens_cgs;
    int it = 0; Real fh = 1, fhii = 1, fc = 1, T = 0, h = 0, cl = 0;
    Real xh_prev = 1.0, xhii_prev = yb[7];
    for (it = 1; it <= 2000; ++it) {
      GOW17Network n2(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
      n2.SetIonizedCooling(true, 1.0);
      for (int n = 0; n < 13; ++n) n2.y(n) = yb[n];
      SemiImplicit<GOW17Network> solver(s1, n2, 0.0, dt);
      solver.SolveODE();
      for (int n = 0; n < 13; ++n) yb[n] = n2.y(n);
      auto g = n2.SetupNextStep(n2.y);
      state(yb, &T, &h, &cl);
      fh = std::fabs(g.H/xh_prev - 1); fhii = std::fabs(yb[7]/xhii_prev - 1);
      fc = std::fabs(h - cl)/(h + cl);
      xh_prev = g.H; xhii_prev = yb[7];
      if (fh < 1e-6 && fhii < 1e-6 && fc < 2.0) break;
    }
    std::printf("(b, chem_only at 20 K) dt %.0e: %s after %d iterations: T %.3f K x_H2 %.4f x_H+ %.3e x_CO %.3e "
                "|H-C|/(H+C) %.1e\n", dt, it <= 2000 ? "converged" : "NOT converged",
                std::min(it, 2000), T, yb[6], yb[7], yb[3], fc);
  }
  return 0;
}
