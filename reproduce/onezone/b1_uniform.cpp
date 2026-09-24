// Test B1: AthenaK's GOW17_uniform problem in one zone with the Tigris copy.
// n_H = 109.21 cm^-3, all species 1e-6 at t = 0, e = n_H cs^2/(gamma-1) code units
// (cs = 1), G0 = 1e-6 in every photo slot and G_PE, CR 2e-16 s^-1, evolved to
// t = 1e4 code units (1 code unit = 3.0857e13 s) in 28 equal steps, as the AthenaK
// test's 28 cycles. Prints T and the 12 species at the end (Gauss-Seidel update).
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;

int main() {
  const Real length_cgs = 3.0856775809623245e18, mass_cgs = 6.882615081124725e31,
             time_cgs = 3.08567758e13;
  const Real dens_cgs = mass_cgs/(length_cgs*length_cgs*length_cgs);
  const Real vel_cgs = length_cgs/time_cgs;
  const Real edens_cgs = dens_cgs*vel_cgs*vel_cgs;
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1.0e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = false;
  s.velocity_cgs = vel_cgs; s.length_cgs = length_cgs; s.multi_d = true; s.three_d = true;
  SemiImplicitSettings si;
  si.semi_implicit_cfl = 0.1; si.semi_implicit_yfloor = 1.0e-12;
  si.semi_implicit_n_substep_max = 100000; si.semi_implicit_n_substep_fixed = 0;
  si.semi_implicit_stiff_threshold = 1.0; si.semi_implicit_n_iter = 1;
  si.semi_implicit_energy = true; si.semi_implicit_renormalize = true;
  si.semi_implicit_hep_first = true; si.semi_implicit_exact_map = false;
  si.semi_implicit_exact_block = false; si.semi_implicit_co_block = true;
  si.semi_implicit_exact_ghosts = false; si.semi_implicit_h2_first = false;
  si.semi_implicit_adaptive = true; si.semi_implicit_nbad_max = 3;
  si.semi_implicit_refresh_rates = false; si.semi_implicit_table_deriv = false;
  Real rad[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 2e-16};
  const Real nH = 109.21, tlim = 1.0e4, gamma = 5.0/3.0;
  const int nstep = 28;
  si.semi_implicit_gauss_seidel = true;
  for (int gs = 1; gs < 2; ++gs) {
    Real y[13];
    for (int n = 0; n < 12; ++n) y[n] = 1.0e-6;
    y[12] = nH*1.0/(gamma - 1.0);  // code units, cs = 1
    long nsub = 0;
    for (int step = 0; step < nstep; ++step) {
      GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens_cgs);
      for (int n = 0; n < 13; ++n) net.y(n) = y[n];
      SemiImplicit<GOW17Network> solver(si, net, 0.0, tlim/nstep);
      solver.SolveODE();
      nsub += solver.n_substeps;
      for (int n = 0; n < 13; ++n) y[n] = net.y(n);
      if (step == nstep - 1) {
        auto g = net.SetupNextStep(net.y);
        std::printf("%d %.10g %.10g", gs, net.Temperature(net.y, g), y[12]);
        for (int n = 0; n < 12; ++n) std::printf(" %.10g", y[n]);
        std::printf(" %ld\n", nsub);
      }
    }
  }
  return 0;
}
