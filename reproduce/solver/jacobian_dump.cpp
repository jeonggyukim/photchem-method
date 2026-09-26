// Species Jacobian along the reference run of the F05 CO-formation case (cosmic rays
// only, n_H 1e3, from atomic H and C+ at 100 K, 2.9 Myr; setup as
// ../figures/F05_solver/solver_sweep.cpp), for the block Gauss-Seidel error matrix of
// Appendix E. At 16 times of the reference run (default update, 12500 substeps per
// sixteenth of t_end) the rate coefficients and ghost species are frozen and
// J_ij = d(C_i - D_i x_i)/dx_j by finite difference of CDRates.
// Rows: k t_code i j J_ij (code time units), with i, j in network enum order; rows
// with j = -1 give the abundance x_i at that time.
// Build as ../figures/F05_solver/solver_sweep.cpp.
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
constexpr int NS = NSPEC;

struct Sample { Real t; Real x[NS]; Real J[NS][NS]; };

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
    for (int i = 0; i < NS; ++i) smp.x[i] = yb[i];
    View1D<Real> y(yb, NSPEC + 1), yp(yq, NSPEC + 1);
    const auto g = net.SetupNextStep(y);
    const auto base = net.CDRates(y, g);
    Real f0[NS];
    for (int i = 0; i < NS; ++i) f0[i] = base.creation(i) - yb[i]*base.destruction(i);
    for (int r = 0; r < NS; ++r) {
      for (int i = 0; i <= NSPEC; ++i) yq[i] = yb[i];
      const Real dx = 1e-6*std::fabs(yb[r]) + 1e-12;
      yq[r] += dx;
      const auto pert = net.CDRates(yp, g);
      for (int s2 = 0; s2 < NS; ++s2)
        smp.J[s2][r] = (pert.creation(s2) - yq[s2]*pert.destruction(s2) - f0[s2])/dx;
    }
    out.push_back(smp);
  }
  t_end_out = t_end;
  myr = time_cgs/3.15576e13;
  return out;
}

int main() {
  std::printf("# k t_code i j J_ij; species order:");
  for (int i = 0; i < NS; ++i) std::printf(" %s", std::string(G::species_names[i]).c_str());
  std::printf("\n");
  Real t_end, myr;
  const auto samples = Reference(2, t_end, myr);
  std::printf("# t_end_code %.6e myr_per_code %.6e\n", t_end, myr);
  for (size_t k = 0; k < samples.size(); ++k)
  {
    for (int i = 0; i < NS; ++i)
      std::printf("%zu %.6e %d -1 %.8e\n", k, samples[k].t, i, samples[k].x[i]);
    for (int i = 0; i < NS; ++i)
      for (int j = 0; j < NS; ++j)
        std::printf("%zu %.6e %d %d %.8e\n", k, samples[k].t, i, j, samples[k].J[i][j]);
  }
}
