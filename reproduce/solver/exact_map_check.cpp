// Does the exponential map help? Error of the semi-implicit solver against the number of
// fixed substeps, one zone, GOW17 core network, the cases of F05: A ionization (xi_H
// 1e-9 s^-1, 3.45 eV per ionization, n_H 100, from atomic 100 K, 978 yr); B recombination
// (from x_H+ 0.99 at 8000 K, same time); C CO formation (cosmic rays only, n_H 1e3, from
// atomic H and C+ at 100 K, 2.9 Myr). Variants, all Gauss-Seidel with the CO/HCO+ block:
// backward Euler (default), exact_map (exponential map for the scalar species),
// exact_block (2x2 matrix exponential for CO/HCO+), and both. Reference: the default at
// 200000 substeps; the exact_map run at 200000 is printed against it as a check.
// Output per row: case N, then per variant |dT/T| |dx/x| (x = x_H+ for A and B, x_CO for C)
// and the CPU time of the run [s].
// Build: configure tigris-gow17 with -gow17 (defs.hpp only), then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$TIGRIS_DIR/src
//   -I/opt/homebrew/opt/hdf5-mpi/include exact_map_check.cpp -o exact_map_check
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
struct Out { Real T, x, cpu; };
enum Variant { DEFAULT, EXACT_MAP, EXACT_BLOCK, BOTH, NVARIANT };

Out Run(int scenario, Variant v, unsigned nsub) {
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
  si.semi_implicit_exact_map = (v == EXACT_MAP || v == BOTH);
  si.semi_implicit_exact_block = (v == EXACT_BLOCK || v == BOTH);
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
  const auto t0 = std::chrono::steady_clock::now();
  SemiImplicit<G> solver(si, net, 0.0, t_end);
  solver.SolveODE();
  const Real cpu = std::chrono::duration<Real>(std::chrono::steady_clock::now() - t0).count();
  auto g = net.SetupNextStep(net.y);
  return {net.Temperature(net.y, g), (scenario == 2) ? net.y(ICO) : net.y(IH_plus), cpu};
}

int main() {
  std::printf("# case N then per variant (default, exact_map, exact_block, both) "
              "|dT/T| |dx/x| cpu[s]; reference: default, 200000 substeps\n");
  for (int sc = 0; sc < 3; ++sc) {
    const Out ref = Run(sc, DEFAULT, 200000);
    const Out chk = Run(sc, BOTH, 200000);
    std::printf("# case %d reference T %.6e K x %.6e; both-exact at 200000: dT/T %.2e "
                "dx/x %.2e; cpu default %.3f s, both %.3f s\n", sc, ref.T, ref.x,
                chk.T/ref.T - 1, chk.x/ref.x - 1, ref.cpu, chk.cpu);
    for (unsigned n : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u, 256u, 512u, 1024u, 2048u, 4096u,
                       16384u}) {
      std::printf("%d %u", sc, n);
      for (int v = 0; v < NVARIANT; ++v) {
        const Out o = Run(sc, static_cast<Variant>(v), n);
        std::printf(" %.6e %.6e %.3e", std::fabs(o.T/ref.T - 1), std::fabs(o.x/ref.x - 1),
                    o.cpu);
      }
      std::printf("\n");
    }
  }
}
