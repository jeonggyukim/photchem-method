// Energy-species coupling within a substep, case C of F05 (CO formation: cosmic rays
// only, n_H 1e3, from atomic H and C+ at 100 K, 2.9 Myr), Gauss-Seidel + CO/HCO+ group.
// Modes: 0 default (energy and species both read the start-of-substep state: Jacobi
// between them), 1 semi_implicit_refresh_rates (species read rates at the updated T:
// Gauss-Seidel), 2 semi_implicit_n_iter = 2, 3 both. A run with N equal substeps is N
// one-substep calls over t_end/N, printing the state after each. The reference (N = 0)
// is mode 0 with 200000 substeps, printed at 64 times; N = 100 is mode 2 with 200000.
// Rows: mode N k t_Myr T x_H2 x_CO x_C+ x_HCO+ x_e.
// Build as solver_sweep.cpp (../figures/F05_solver).
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;

void Trajectory(unsigned nout, unsigned nsub_per, unsigned label, int mode) {
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = 1.0e3, t_end = 3.0;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  s.riccati_hp = false;
  SemiImplicitSettings si{0.1, 1e-12, 100000000, nsub_per, 1.0, 1, true, true, true, true,
                          false, false, true, false, false, true, 3, false, false};
  si.semi_implicit_gauss_seidel = true;
  si.semi_implicit_co_block = true;
  si.semi_implicit_refresh_rates = (mode == 1 || mode == 3);
  si.semi_implicit_n_iter = (mode >= 2) ? 2 : 1;
  Real rad[G::n_freq];
  for (int n = 0; n < G::n_freq - 1; ++n) rad[n] = 0.0;
  rad[G::n_freq - 1] = 2.0e-16;
  G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  const Real ev = 1.602176634e-12;
  net.SetLyC(0.0, 0.0, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
  net.y(IC_plus) = 1.6e-4;
  net.y(ISi_plus) = 1.7e-6;
  net.y(IH_plus) = 1.0e-4;
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = 100.0*Thermo::CvCold(0.0, 0.1, g0.e, gamma)*nH/edens;
  const Real myr = time_cgs/3.15576e13, dt = t_end/nout;
  for (unsigned k = 0; k <= nout; ++k) {
    if (k > 0) {
      SemiImplicit<G> solver(si, net, (k - 1)*dt, dt);
      solver.SolveODE();
    }
    auto g = net.SetupNextStep(net.y);
    std::printf("%d %u %u %.4e %.5e %.5e %.5e %.5e %.5e %.5e\n", mode, label, k, k*dt*myr,
                net.Temperature(net.y, g), net.y(gow17::IH2), net.y(ICO), net.y(IC_plus),
                net.y(IHCO_plus), g.e);
  }
}

int main() {
  std::printf("# mode N k t_Myr T x_H2 x_CO x_C+ x_HCO+ x_e; N = 0 is the reference\n");
  Trajectory(64, 3125, 0, 0);
  Trajectory(64, 3125, 100, 2);
  for (int mode = 0; mode < 4; ++mode)
    for (unsigned n : {8u, 16u, 32u, 64u, 256u, 1024u, 4096u}) Trajectory(n, 1, n, mode);
}
