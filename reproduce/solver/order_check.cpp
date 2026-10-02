// Which Gauss-Seidel update order? Error of the semi-implicit solver (Gauss-Seidel +
// CO/HCO+ block, backward Euler) against the number of fixed substeps, one zone, GOW17
// core network, for the orders set by semi_implicit_order:
//   0 default          He+ Si+ SN H2+ H3+ H+ O+ C+ CHx OHx CO H2 (built-in sequence)
//   1 default, listed  the same as a list; must agree with 0 to round-off
//   2 low lag          H2 Si+ H2+ H3+ OHx CO He+ C+ O+ CHx H+   (order_search.py, lag)
//   3 rho+lag          H2 H2+ H3+ CHx Si+ OHx CO He+ O+ H+ C+   (order_search.py, rho+lag)
// Cases as h2_order_check: A ionization, B recombination, C CO formation, D
// photodissociation. Reference: order 0 at 200000 substeps. Output per row: case N,
// then per order |dT/T|; order_check <order> ... compares the given orders with order 0
// instead of 1-3. Per order |dy/y| max_s|dx_s/x_s| (y = x_H+ for A and B, x_CO for C, x_H2
// for D; the max over species whose reference abundance exceeds 1e-8).
// Build: tigris-gow17 configured with -gow17 (defs.hpp only), then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$TIGRIS_DIR/src
//   -I/opt/homebrew/opt/hdf5-mpi/include order_check.cpp -o order_check
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
struct Out { Real T; std::vector<Real> x; };

std::vector<std::string> ORDERS = {"", "He+,Si+,SN,H2+,H3+,H+,O+,C+,CHx,OHx,CO,H2",
                                   "H2,Si+,H2+,H3+,OHx,CO,He+,C+,O+,CHx,H+",
                                   "H2,H2+,H3+,CHx,Si+,OHx,CO,He+,O+,H+,C+"};
int NORDER = 4;

Out Run(int scenario, int iorder, unsigned nsub) {
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
  std::string order = ORDERS[iorder];
  if (order.rfind("b:", 0) == 0) {
    si.semi_implicit_h_block = true;
    order = order.substr(2);
  }
  si.semi_implicit_n_order = SemiImplicit<G>::ParseOrder(order, si.semi_implicit_order);
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
  if (scenario == 3) net.y(gow17::IH2) = 0.5*(1.0 - 1.0e-4);
  const Real T0 = (scenario == 1) ? 8000.0 : 100.0;
  auto g0 = net.SetupNextStep(net.y);
  net.y(NSPEC) = T0*Thermo::CvCold(net.y(gow17::IH2), 0.1, g0.e, gamma)*nH/edens;
  SemiImplicit<G> solver(si, net, 0.0, t_end);
  solver.SolveODE();
  auto g = net.SetupNextStep(net.y);
  Out o{net.Temperature(net.y, g), std::vector<Real>(NSPEC)};
  for (int n = 0; n < NSPEC; ++n) o.x[n] = net.y(n);
  return o;
}

// Largest relative error over the species whose reference abundance exceeds 1e-8, and
// the species that carries it
Real MaxError(const Out& o, const Out& ref, int* which) {
  Real emax = 0.0;
  *which = 0;
  for (int s = 0; s < NSPEC; ++s) {
    const Real e = std::fabs(o.x[s]/ref.x[s] - 1);
    if (ref.x[s] > 1e-8 && e > emax) { emax = e; *which = s; }
  }
  return emax;
}

// Local search over unit orders by moving one unit to another position, scored by
// the mean over N = 16, 64, 256 of log10 max(|dT/T|, max_s |dx_s/x_s|) in one case.
void Search(int sc) {
  const unsigned NS[] = {16u, 64u, 256u};
  const Out ref = Run(sc, 0, 200000);
  auto score = [&](const std::vector<std::string>& u) {
    std::string s;
    for (const auto& n : u) s += (s.empty() ? "" : ",") + n;
    ORDERS.push_back(s);
    Real sum = 0.0;
    for (unsigned n : NS) {
      const Out o = Run(sc, static_cast<int>(ORDERS.size()) - 1, n);
      int w;
      sum += std::log10(std::fmax(std::fabs(o.T/ref.T - 1), MaxError(o, ref, &w)));
    }
    ORDERS.pop_back();
    return sum/3.0;
  };
  std::vector<std::string> best = {"He+", "Si+", "H2+", "H3+", "H+", "O+", "C+", "CHx",
                                   "OHx", "CO", "H2"};
  Real best_s = score(best);
  std::printf("# case %d default score %.3f\n", sc, best_s);
  for (bool improved = true; improved;) {
    improved = false;
    for (std::size_t a = 0; a < best.size(); ++a) {
      for (std::size_t b = 0; b < best.size(); ++b) {
        if (a == b) continue;
        auto cand = best;
        const std::string u = cand[a];
        cand.erase(cand.begin() + a);
        cand.insert(cand.begin() + b, u);
        const Real s = score(cand);
        if (s < best_s - 1e-3) { best = cand; best_s = s; improved = true; }
      }
    }
  }
  std::string s;
  for (const auto& n : best) s += (s.empty() ? "" : ",") + n;
  std::printf("# case %d best score %.3f  %s\n", sc, best_s, s.c_str());
}

int main(int argc, char** argv) {
  // order_check search <case>: local search for the best order in one case
  if (argc == 3 && std::string(argv[1]) == "search") {
    Search(std::atoi(argv[2]));
    return 0;
  }
  // extra arguments: orders to compare against order 0 in place of 1-3; a "b:"
  // prefix also sets semi_implicit_h_block
  if (argc > 1) {
    ORDERS.resize(1);
    for (int a = 1; a < argc; ++a) ORDERS.push_back(argv[a]);
    NORDER = static_cast<int>(ORDERS.size());
  }
  std::printf("# case N then per order (0 default, 1 default listed, 2 low lag, 3 rho+lag) "
              "|dT/T| |dy/y| max_s|dx/x|; reference: order 0, 200000 substeps\n");
  for (int sc = 0; sc < 4; ++sc) {
    const int iy = (sc == 2) ? ICO : (sc == 3) ? gow17::IH2 : IH_plus;
    const Out ref = Run(sc, 0, 200000);
    std::printf("# case %d reference T %.6e K y %.6e;", sc, ref.T, ref.x[iy]);
    for (int io = 1; io < NORDER; ++io) {
      const Out chk = Run(sc, io, 200000);
      std::printf(" order %d at 200000: dT/T %.2e dy/y %.2e", io, chk.T/ref.T - 1,
                  chk.x[iy]/ref.x[iy] - 1);
    }
    std::printf("\n");
    for (unsigned n : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u, 256u, 512u, 1024u, 4096u,
                       16384u}) {
      std::printf("%d %u", sc, n);
      for (int io = 0; io < NORDER; ++io) {
        const Out o = Run(sc, io, n);
        int w;
        const Real emax = MaxError(o, ref, &w);
        std::printf(" %.6e %.6e %.6e %s", std::fabs(o.T/ref.T - 1),
                    std::fabs(o.x[iy]/ref.x[iy] - 1), emax, G::species_names[w].data());
      }
      std::printf("\n");
    }
  }
}
