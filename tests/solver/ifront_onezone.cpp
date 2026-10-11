// One-zone ionization (A: xi = 1e-9 s^-1 on atomic gas) and recombination (B: source
// off in ionized gas), n_H = 100 cm^-3. Backward-Euler vs Riccati H+ with N fixed
// substeps against backward Euler with 200000 substeps; x_H+, x_e and T at the end.
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
struct Out { Real xhp, xe, T; };
Out Run(int scenario, bool riccati, unsigned nsub, Real t_end) {
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10, nH = 100.0;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  s.riccati_hp = riccati;
  SemiImplicitSettings si{0.1, 1e-12, 100000000, nsub, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  Real rad[8] = {1, 1, 1, 1, 0, 1, 1, 2e-16};
  GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  const Real ev = 1.602176634e-12;
  net.SetLyC(scenario == 0 ? 1e-9 : 0.0, scenario == 0 ? 1.5e-9 : 0.0, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
  net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6;
  const Real xhp0 = scenario == 0 ? 1e-4 : 0.99, T0 = scenario == 0 ? 100.0 : 8000.0;
  net.y(IH_plus) = xhp0;
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = T0*Thermo::CvCold(0.0, 0.1, g0.e, gamma)*nH/edens;
  SemiImplicit<GOW17Network> solver(si, net, 0.0, t_end);
  solver.SolveODE();
  auto g = net.SetupNextStep(net.y);
  return {net.y(IH_plus), g.e, net.Temperature(net.y, g)};
}
int main() {
  const char *name[2] = {"A ionization (xi 1e-9)", "B recombination (xi 0)"};
  for (int sc = 0; sc < 2; ++sc) {
    for (Real t : {1e-3, 1e-2}) {  // code time: 978 yr, 9780 yr
      const Out ref = Run(sc, false, 200000, t);
      std::printf("%s, t = %.0f yr: reference x_H+ %.6f x_e %.6f T %.1f K\n", name[sc], t*0.978e6,
                  ref.xhp, ref.xe, ref.T);
      std::printf("  %6s | %28s | %28s\n", "N", "BE: dx_H+/x  dT/T", "Riccati: dx_H+/x  dT/T");
      for (unsigned n : {1u, 4u, 16u, 64u}) {
        const Out be = Run(sc, false, n, t), ri = Run(sc, true, n, t);
        std::printf("  %6u | %13.2e %13.2e | %13.2e %13.2e\n", n, be.xhp/ref.xhp - 1, be.T/ref.T - 1,
                    ri.xhp/ref.xhp - 1, ri.T/ref.T - 1);
      }
    }
  }
}
