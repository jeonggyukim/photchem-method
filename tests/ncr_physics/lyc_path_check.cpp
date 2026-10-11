// A cell exposed to Lyman-continuum photoionization, integrated to steady state by
// the Jacobi (CDRates) and the Gauss-Seidel update paths. Both must reach the
// photoionization balance x_HI xi = alpha_B(T) n_H x_e x_H+.
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;

int main() {
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5;
  const Real edens = 1.6738234e-24*1.0e10, nH = 32.0, xi_hi = 1e-9, xi_h2 = 1.5e-9;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = false;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  SemiImplicitSettings si{0.1, 1e-12, 100000, 0, 1.0, 1, true, true, false, true, false,
                          false, true, false, false, true, 3, false, false};
  for (int gs = 0; gs <= 1; ++gs) {
    si.semi_implicit_gauss_seidel = (gs == 1);
    Real rad[8] = {1, 1, 1, 1, 0, 1, 1, 2e-16};
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
    net.SetLyC(xi_hi, xi_h2, 3.45*1.602176634e-12, 4.42*1.602176634e-12);
    net.SetIonizedCooling(true, 1.0);
    for (int n = 0; n < 13; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 0.01; net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6;
    net.y(12) = 100.0*Thermo::CvCold(0.0, 0.1, 0.01, gamma)*nH/edens;
    int nsub = 0;
    for (int step = 0; step < 20; ++step) {
      SemiImplicit<GOW17Network> solver(si, net, 0.0, 0.05);
      solver.SolveODE(); nsub += solver.n_substeps;
      if (step % 4 == 3) { auto gg = net.SetupNextStep(net.y); std::printf("  t %.2f T %.1f x_HI %.4e\n", 0.05*(step+1), net.Temperature(net.y, gg), gg.H); }
    }
    auto g = net.SetupNextStep(net.y);
    const Real T = net.Temperature(net.y, g);
    const Real alphaB = 2.59e-13*std::pow(T*1e-4, -0.7);
    std::printf("%s: %d substeps, T %.0f K, x_H+ %.6f x_HI %.3e x_H2 %.2e x_e %.4f; "
                "x_HI xi/(alpha_B n_H x_e x_H+) = %.3f\n",
                gs ? "Gauss-Seidel" : "Jacobi", nsub, T, net.y(IH_plus), g.H,
                net.y(gow17::IH2), g.e, g.H*xi_hi/(alphaB*nH*g.e*net.y(IH_plus)));
  }
  return 0;
}
