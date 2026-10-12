// He collisional ionization equilibrium: an isothermal cell with no radiation and the
// O3,S3,N3 ion set (so He++ is tracked), integrated to steady state, against the
// CHIANTI v11 CIE fractions (inputs/tables/chianti_v11/ioneq_He.txt). He0 is the
// remainder of the He total.
// Build: configure tigris-gow17 with -gow17 --photchem_ions=O3,S3,N3 (defs.hpp), then
//   mpicxx -std=c++17 -O2 -I$TIGRIS_DIR/src -I<hdf5 prefix>/include he_cie_check.cpp \
//     $TIGRIS_DIR/src/photchem/network/gow17_thermo_table.cpp -o he_cie_check
#include <cmath>
#include <cstdio>
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
  std::ifstream f(TigrisDir() + "inputs/tables/chianti_v11/ioneq_" + el + ".txt");
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
      for (std::size_t q = 0; q < x.size(); ++q) {
        x[q] = (1 - w)*rows[i][q + 1] + w*rows[i + 1][q + 1];
      }
      return x;
    }
  }
  return {};
}

int main() {
  const std::string d = TigrisDir() + "inputs/tables/rates/";
  Rates::RecombRate rec;
  rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat",
           d + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci;
  ci.Load(d + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct;
  ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid;
  BuildThermoTable(grid);
  std::vector<IonStage> st;
  for (int k = 0; k < GOW17Network::n_ion_stage; ++k) {
    st.push_back(GOW17Network::IonLadderStage(k));
  }
  IonRateTable tab;
  BuildIonRateTable(grid, st, rec, ci, ct, &tab);
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = 1.0, xHe = 0.1;
  std::printf("# He CIE, n_H = 1 cm^-3: fractions of He from Tigris (case A) and CHIANTI v11\n");
  std::printf("# T[K] He0 He+ He++ | chianti He0 He+ He++\n");
  for (int i = 0; i <= 24; ++i) {
    const Real T = std::pow(10.0, 4.2 + 0.1*i);
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = true; s.isothermal_temperature = T;
    s.zd = 1.0; s.xHe = xHe; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
    s.xS = 1.45e-5; s.xN = 7.4e-5;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
    SemiImplicitSettings si{0.1, 1e-12, 100000000, 0, 1.0, 1, true, true, true, true, false,
                            false, true, false, false, true, 3, false, false};
    Real rad[GOW17Network::n_freq] = {};
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, GOW17Network::n_freq), gamma,
                     time_cgs, edens);
    net.SetIonRates(&tab, &grid);
    net.SetRecombinationOTS(false);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0; net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6;
    net.y(IHe_plus) = 0.5*xHe;
    net.y(NSPEC) = 1.0;
    for (int k = 0; k < 40; ++k) {
      SemiImplicit<GOW17Network> sv(si, net, 0.0, 10.0);
      sv.SolveODE();
    }
    auto g = net.SetupNextStep(net.y);
    const auto c = Ioneq("He", T);
    std::printf("%.4e %.4e %.4e %.4e | %.4e %.4e %.4e\n", T, g.He/xHe,
                net.y(IHe_plus)/xHe, net.y(GOW17Network::IHE_2plus)/xHe, c[0], c[1], c[2]);
  }
}
