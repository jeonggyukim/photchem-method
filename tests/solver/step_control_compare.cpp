// Substeps and accuracy: current limiter (cfl 0.1) vs change-based control, one
// zone, against backward Euler with 200000 fixed substeps. Errors: T, and the
// largest species error measured relative to max(y_ref, 1e-3 x element scale).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
struct Case { const char *name; Real nH, chi, xi, T0, xhp0, xh2_0, t_end; bool b1; };
struct Out { Real T; Real y[NSPEC]; int nsub; };
Out Run(const Case &c, int mode, bool riccati) {  // mode 0 ref, 1 limiter, 2 change
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
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
  SemiImplicitSettings si{0.1, 1e-12, 100000000, mode == 0 ? 200000u : 0u, 1.0, 1, true, true,
                          true, true, false, false, true, false, false, true, 3, false, false};
  si.semi_implicit_change_control = (mode == 2);
  if (const char *v = std::getenv("FLOOR")) si.semi_implicit_floor_frac = std::atof(v);
  if (const char *v = std::getenv("TARGET")) { si.semi_implicit_change_target = std::atof(v); si.semi_implicit_change_reject = 2*std::atof(v); }
  Real rad[8] = {c.chi, c.chi, c.chi, c.chi, c.b1 ? c.chi : 0.0, c.chi, c.chi, 2e-16};
  GOW17Network net(s, c.nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  const Real ev = 1.602176634e-12;
  net.SetLyC(c.xi, 1.5*c.xi, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC; ++n) net.y(n) = c.b1 ? 1e-6 : 0.0;
  if (!c.b1) { net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6; net.y(IH_plus) = c.xhp0; net.y(gow17::IH2) = c.xh2_0; }
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = c.T0*Thermo::CvCold(net.y(gow17::IH2), 0.1, g0.e, gamma)*c.nH/edens;
  Out o{0, {}, 0};
  const int nstep = c.b1 ? 28 : 1;  // B1: 28 hydro steps; else one
  for (int k = 0; k < nstep; ++k) {
    SemiImplicit<GOW17Network> solver(si, net, 0.0, c.t_end/nstep);
    solver.SolveODE(); o.nsub += solver.n_substeps;
  }
  auto g = net.SetupNextStep(net.y);
  o.T = net.Temperature(net.y, g);
  for (int n = 0; n < NSPEC; ++n) o.y[n] = net.y(n);
  return o;
}
int main() {
  const Case cases[] = {
    {"ionization  n100 xi1e-9 1e3yr", 100, 1, 1e-9, 100, 1e-4, 0, 1e-3, false},
    {"ionization  n100 xi1e-9 1e4yr", 100, 1, 1e-9, 100, 1e-4, 0, 1e-2, false},
    {"recombine   n100 xi0    1e4yr", 100, 1, 0, 8000, 0.99, 0, 1e-2, false},
    {"unshielded  n100 chi1   1e5yr", 100, 1, 0, 100, 1e-3, 1e-4, 0.1, false},
    {"B1 shielded 28 steps to 1e4",  109.21, 1e-6, 0, 0, 0, 0, 1e4, true},
  };
  // B1 energy: cs = 1 code, e = nH/(gamma-1) in code units, T0 set below via E
  std::printf("%-30s %8s | %8s %9s %9s | %8s %9s %9s | %8s %9s %9s\n", "case", "ref",
              "lim nsub", "dT/T", "dy max", "chg nsub", "dT/T", "dy max", "chg+Ric", "dT/T", "dy max");
  for (const Case &c : cases) {
    Case cc = c;
    if (c.b1) cc.T0 = 45.0;  // start warm; the run relaxes
    const Out ref = Run(cc, 0, false);
    std::printf("%-30s %8.1f |", c.name, ref.T);
    for (int m : {1, 2, 3}) {
      const Out o = Run(cc, m == 3 ? 2 : m, m == 3);
      Real dy = 0;
      for (int n = 0; n < NSPEC; ++n) {
        GOW17Settings dummy{}; (void)dummy;
        const Real scale = std::max(ref.y[n], 1e-3*1.6e-4);
        dy = std::max(dy, std::fabs(o.y[n] - ref.y[n])/scale);
      }
      std::printf(" %8d %9.2e %9.2e |", o.nsub, o.T/ref.T - 1, dy);
    }
    std::printf("\n");
  }
}
