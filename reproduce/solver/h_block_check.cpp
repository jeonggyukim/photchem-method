// Does solving H+ and H2 together (semi_implicit_h_block) help? Error of the
// semi-implicit solver (Gauss-Seidel + CO/HCO+ block, backward Euler) against the number
// of fixed substeps, one zone, GOW17 core network, default sweep (H+ at its place, H2
// last) against the H+/H2 2x2 block. Cases: A ionization of atomic gas, B recombination,
// C CO formation, D photodissociation of molecular gas (the cases of h2_order_check),
// and E an ionization front into molecular gas: 2 x_H2 = 1 at n_H 100, chi = 1, LyC
// photoionization of H (1e-9 s^-1) and H2 (1.5e-9 s^-1), from 100 K, 978 yr.
// Reference: the default at 200000 substeps; the block at 200000 is printed against it.
// Output per row: case N, then per variant |dT/T| |dx_H+/x_H+| |dx_H2/x_H2| |dx_CO/x_CO|.
// Build: configure tigris-gow17 with -gow17 (defs.hpp only), then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$HOME/Projects/tigris-gow17/src
//   -I/opt/homebrew/opt/hdf5-mpi/include h_block_check.cpp -o h_block_check
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
struct Out { Real T, hp, h2, co; };

Out Run(int scenario, bool h_block, unsigned nsub) {
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
  SemiImplicitSettings si{0.1, 1e-12, 100000000, nsub, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  si.semi_implicit_gauss_seidel = true;
  si.semi_implicit_co_block = true;
  si.semi_implicit_h_block = h_block;
  Real rad[G::n_freq];
  const Real chi = (scenario == 2) ? 0.0 : 1.0;
  for (int n = 0; n < G::n_freq - 1; ++n) rad[n] = chi;
  rad[G::n_freq - 1] = 2.0e-16;
  G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
  net.SetIonizedCooling(true, 1.0);
  const Real ev = 1.602176634e-12;
  const bool lyc = (scenario == 0 || scenario == 4);
  net.SetLyC(lyc ? 1e-9 : 0.0, lyc ? 1.5e-9 : 0.0, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
  net.y(IC_plus) = 1.6e-4;
  net.y(ISi_plus) = 1.7e-6;
  net.y(IH_plus) = (scenario == 1) ? 0.99 : 1.0e-4;
  if (scenario == 3 || scenario == 4) net.y(gow17::IH2) = 0.5*(1.0 - 1.0e-4);
  const Real T0 = (scenario == 1) ? 8000.0 : 100.0;
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = T0*Thermo::CvCold(net.y(gow17::IH2), 0.1, g0.e, gamma)*nH/edens;
  SemiImplicit<G> solver(si, net, 0.0, t_end);
  solver.SolveODE();
  auto g = net.SetupNextStep(net.y);
  return {net.Temperature(net.y, g), net.y(IH_plus), net.y(gow17::IH2), net.y(ICO)};
}

int main() {
  std::printf("# case N then per variant (default, H+/H2 block) |dT/T| |dx_H+/x_H+| "
              "|dx_H2/x_H2| |dx_CO/x_CO|; reference: default, 200000 substeps\n");
  auto rel = [](Real a, Real b) { return std::fabs(a/b - 1); };
  for (int sc = 0; sc < 5; ++sc) {
    const Out ref = Run(sc, false, 200000);
    const Out chk = Run(sc, true, 200000);
    std::printf("# case %d reference T %.5e H+ %.5e H2 %.5e CO %.5e; block at 200000: "
                "dT %.1e dH+ %.1e dH2 %.1e dCO %.1e\n", sc, ref.T, ref.hp, ref.h2, ref.co,
                rel(chk.T, ref.T), rel(chk.hp, ref.hp), rel(chk.h2, ref.h2), rel(chk.co, ref.co));
    for (unsigned n : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u, 256u, 1024u, 4096u, 16384u}) {
      std::printf("%d %u", sc, n);
      for (bool blk : {false, true}) {
        const Out o = Run(sc, blk, n);
        std::printf(" %.4e %.4e %.4e %.4e", rel(o.T, ref.T), rel(o.hp, ref.hp),
                    rel(o.h2, ref.h2), rel(o.co, ref.co));
      }
      std::printf("\n");
    }
  }
}
