// Gauss-Seidel coupling factor of every pair of integrated species along the reference
// trajectories of the F05 cases (A ionization, B recombination, C CO formation; setup as
// ../figures/F05_solver/solver_sweep.cpp). At 16 times of each reference run (default
// update, 12500 substeps per sixteenth of t_end) the rate coefficients and ghost species
// are frozen, D_s is the destruction rate and K_sr = dC_s/dx_r the rate at which member r
// forms member s (finite difference of CDRates). For substep h = t_end/N,
// theta_sr = h K_sr/(1 + h D_s) and rho = theta_sr theta_rs, the factor by which one
// Gauss-Seidel pass over the 2x2 pair reduces the error of the lagged partner.
// Rows: case N pair max_t(rho) t_Myr_at_max theta_sr theta_rs  (pairs with max rho > 1e-4).
// Build as solver_sweep.cpp.
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
constexpr int NS = NSPEC;

struct Sample { Real t; Real D[NS]; Real K[NS][NS]; };

std::vector<Sample> Reference(int scenario, Real& t_end_out, Real& myr) {
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = (scenario == 2) ? 1.0e3 : 100.0;
  const Real t_end = (scenario == 2) ? 3.0 : 1.0e-3;
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
  const unsigned nout = 16;
  SemiImplicitSettings si{0.1, 1e-12, 100000000, 200000/nout, 1.0, 1, true, true, true, true,
                          false, false, true, false, false, true, 3, false, false};
  si.semi_implicit_gauss_seidel = true;
  si.semi_implicit_co_block = true;
  Real rad[G::n_freq];
  const Real chi = (scenario == 2) ? 0.0 : 1.0;
  for (int n = 0; n < G::n_freq - 1; ++n) rad[n] = chi;
  rad[G::n_freq - 1] = 2.0e-16;
  G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  const Real ev = 1.602176634e-12;
  net.SetLyC(scenario == 0 ? 1e-9 : 0.0, scenario == 0 ? 1.5e-9 : 0.0, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
  net.y(IC_plus) = 1.6e-4;
  net.y(ISi_plus) = 1.7e-6;
  net.y(IH_plus) = (scenario == 1) ? 0.99 : 1.0e-4;
  const Real T0 = (scenario == 1) ? 8000.0 : 100.0;
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = T0*Thermo::CvCold(0.0, 0.1, g0.e, gamma)*nH/edens;
  const Real dt = t_end/nout;
  std::vector<Sample> out;
  for (unsigned k = 1; k <= nout; ++k) {
    SemiImplicit<G> solver(si, net, (k - 1)*dt, dt);
    solver.SolveODE();
    Sample smp;
    smp.t = k*dt;
    // net.y is a View: copy the values, not the handle, before perturbing.
    Real yb[NSPEC + 1], yq[NSPEC + 1];
    for (int i = 0; i <= NSPEC; ++i) yb[i] = net.y(i);
    View1D<Real> y(yb, NSPEC + 1), yp(yq, NSPEC + 1);
    const auto g = net.SetupNextStep(y);
    const auto base = net.CDRates(y, g);
    for (int i = 0; i < NS; ++i) smp.D[i] = base.destruction(i);
    for (int r = 0; r < NS; ++r) {
      for (int i = 0; i <= NSPEC; ++i) yq[i] = yb[i];
      const Real dx = 1e-6*std::fabs(yb[r]) + 1e-12;
      yq[r] += dx;
      const auto pert = net.CDRates(yp, g);
      for (int s2 = 0; s2 < NS; ++s2)
        smp.K[s2][r] = (s2 == r) ? 0.0 : (pert.creation(s2) - base.creation(s2))/dx;
    }
    out.push_back(smp);
  }
  t_end_out = t_end;
  myr = time_cgs/3.15576e13;
  return out;
}

int main() {
  std::printf("# case N pair max_rho t_Myr theta_sr theta_rs\n");
  const char* cname[3] = {"A", "B", "C"};
  for (int sc = 0; sc < 3; ++sc) {
    Real t_end, myr;
    const auto samples = Reference(sc, t_end, myr);
    for (unsigned n : {4u, 16u, 64u, 256u}) {
      const Real h = t_end/n;
      for (int s = 0; s < NS; ++s) for (int r = s + 1; r < NS; ++r) {
        Real best = 0, tb = 0, a = 0, b = 0;
        for (const auto& p : samples) {
          const Real th_sr = h*p.K[s][r]/(1 + h*p.D[s]);
          const Real th_rs = h*p.K[r][s]/(1 + h*p.D[r]);
          const Real rho = th_sr*th_rs;
          if (rho > best) { best = rho; tb = p.t*myr; a = th_sr; b = th_rs; }
        }
        if (best > 1e-4)
          std::printf("%s %u %s/%s %.3e %.3e %.3e %.3e\n", cname[sc], n,
                      std::string(G::species_names[s]).c_str(),
                      std::string(G::species_names[r]).c_str(), best, tb, a, b);
      }
    }
  }
}
