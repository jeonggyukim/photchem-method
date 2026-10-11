// M4 check: an isothermal cell with no radiation, O2,S3,N2 ion set, integrated to
// steady state, against the CHIANTI v11 CIE fractions (inputs/tables/chianti_v11/
// ioneq_<El>.txt). The top tracked stage has no ionization out, so it is compared
// with CHIANTI's stages at and above it summed. CHIANTI's CIE has no charge
// transfer; the network does (with the H0 left at that T).
// Build: configure tigris-gow17 with -gow17 --gow17_ions=O2,S3,N2 (defs.hpp), then
//   g++-16 -std=c++17 -O2 -I$TIGRIS_DIR/src cie_convergence_check.cpp
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
#include "../paths.hpp"
using namespace gow17;

// CHIANTI CIE fractions of element el at T, linear in log T between rows.
std::vector<Real> Ioneq(const std::string &el, Real T) {
  const std::string fn = TigrisDir() + "inputs/tables/chianti_v11/ioneq_" + el + ".txt";
  std::ifstream f(fn);
  std::string line;
  std::vector<std::vector<Real>> rows;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::vector<Real> r;
    Real v;
    while (ss >> v) r.push_back(v);
    rows.push_back(r);
  }
  const Real lt = std::log10(T);
  for (std::size_t i = 0; i + 1 < rows.size(); ++i) {
    if (rows[i][0] <= lt && lt <= rows[i + 1][0]) {
      const Real w = (lt - rows[i][0])/(rows[i + 1][0] - rows[i][0]);
      std::vector<Real> x(rows[i].size() - 1);
      for (std::size_t q = 0; q < x.size(); ++q) x[q] = (1 - w)*rows[i][q + 1] + w*rows[i + 1][q + 1];
      return x;
    }
  }
  return {};
}

int main() {
  const std::string d = TigrisDir() + "inputs/tables/rates/";
  Rates::RecombRate rec;
  rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat", d + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci;
  ci.Load(d + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct;
  ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid;
  BuildThermoTable(grid);
  std::vector<IonStage> st;
  for (int k = 0; k < GOW17Network::n_ion_stage; ++k) st.push_back(GOW17Network::IonLadderStage(k));
  IonRateTable tab;
  BuildIonRateTable(grid, st, rec, ci, ct, &tab);
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = 1.0, xO = 3.2e-4, xS = 1.45e-5, xN = 7.4e-5;
  std::printf("# T[K] el stage tigris chianti (top stage: chianti summed over stages >= top)\n");
  for (Real T : {2.0e4, 3.0e4, 5.0e4, 7.0e4, 1.0e5}) {
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = true; s.isothermal_temperature = T;
    s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = xO; s.xSi = 1.7e-6; s.xS = xS; s.xN = xN;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
    SemiImplicitSettings si{0.1, 1e-12, 100000000, 0, 1.0, 1, true, true, true, true, false,
                            false, true, false, false, true, 3, false, false};
    Real rad[GOW17Network::n_freq] = {};
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, GOW17Network::n_freq), gamma, time_cgs,
                     edens);
    net.SetIonRates(&tab, &grid);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0; net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6;
    net.y(NSPEC) = 1.0;
    for (int k = 0; k < 40; ++k) {
      SemiImplicit<GOW17Network> sv(si, net, 0.0, 10.0);
      sv.SolveODE();
    }
    auto g = net.SetupNextStep(net.y);
    struct El { const char *name; Real tot; std::vector<Real> x; };
    El els[3] = {{"O", xO, {g.O, net.y(IO_plus), net.y(IO_2plus)}},
                 {"S", xS, {g.S, net.y(IS_plus), net.y(IS_2plus), net.y(IS_3plus)}},
                 {"N", xN, {g.N, net.y(IN_plus), net.y(IN_2plus)}}};
    for (auto &e : els) {
      const auto c = Ioneq(e.name, T);
      const int top = static_cast<int>(e.x.size()) - 1;
      for (int q = 0; q <= top; ++q) {
        Real cq = c[q];
        if (q == top) for (std::size_t k = top + 1; k < c.size(); ++k) cq += c[k];
        std::printf("%.0e %s %d %.4e %.4e\n", T, e.name, q, e.x[q]/e.tot, cq);
      }
    }
    std::printf("%.0e H+ - %.4e -\n", T, net.y(IH_plus));
  }
}
