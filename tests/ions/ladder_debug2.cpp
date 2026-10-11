// O and S ladders (--gow17_ions=O2,S2) in a collisionally ionized cell at fixed T,
// no radiation, integrated to steady state: each pair must satisfy
// rate(q -> q+1) x_q = rate(q+1 -> q) x_q+1 with the ladder's own rates.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
#include "../paths.hpp"
using namespace gow17;
int main() {
  const std::string d = TigrisDir() + "inputs/tables/rates/";
  Rates::RecombRate rec; rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat", d + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci; ci.Load(d + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat"); ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid; BuildThermoTable(grid);
  std::vector<IonStage> st; for (int k = 0; k < GOW17Network::n_ion_stage; ++k) st.push_back(GOW17Network::IonLadderStage(k));
  IonRateTable tab; BuildIonRateTable(grid, st, rec, ci, ct, &tab);
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10, nH = 1.0;
  for (Real T : {2.0e4, 3.0e4}) {
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = true; s.isothermal_temperature = T;
    s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6; s.xS = 1.45e-5;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
    SemiImplicitSettings si{0.1, 1e-12, 100000000, 0, 1.0, 1, true, true, true, true, false,
                            false, true, false, false, true, 3, false, false};
    Real rad[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
    net.SetIonRates(&tab, &grid);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0; net.y(IC_plus) = 1.6e-4; net.y(ISi_plus) = 1.7e-6;
    net.y(NSPEC) = 1.0;
    for (int k = 0; k < 40; ++k) { SemiImplicit<GOW17Network> sv(si, net, 0.0, 10.0); sv.SolveODE(); }
    auto g = net.SetupNextStep(net.y);
    const Real ne = nH*g.e, xh = g.H, xhp = net.y(IH_plus);
    auto K = [&](int k, int c) { return tab.At(grid.Locate(T), k, c); };
    using G = GOW17Network;
    const Real xS0 = g.S, xS1 = net.y(IS_plus), xS2 = net.y(IS_2plus), xO1 = net.y(IO_plus), xO2 = net.y(IO_2plus);
    const Real up_s0 = (K(G::kS0, IonRateTable::ICI)*ne + K(G::kS0, IonRateTable::ICT_ION)*nH*xhp)*xS0;
    const Real dn_s1 = (K(G::kS1, IonRateTable::IREC)*ne + K(G::kS1, IonRateTable::ICT_REC)*nH*xh)*xS1;
    const Real up_s1 = (K(G::kS1, IonRateTable::ICI)*ne + K(G::kS1, IonRateTable::ICT_ION)*nH*xhp)*xS1;
    const Real dn_s2 = (K(G::kS2, IonRateTable::IREC)*ne + K(G::kS2, IonRateTable::ICT_REC)*nH*xh)*xS2;
    const Real up_o1 = (K(G::kO1, IonRateTable::ICI)*ne + K(G::kO1, IonRateTable::ICT_ION)*nH*xhp)*xO1;
    const Real dn_o2 = (K(G::kO2, IonRateTable::IREC)*ne + K(G::kO2, IonRateTable::ICT_REC)*nH*xh)*xO2;
    {
      // one Gauss-Seidel sweep with a huge step: O++ lands on its C/D
      Real yold[NSPEC + 1]; for (int n = 0; n < NSPEC + 1; ++n) yold[n] = net.y(n);
      auto g2 = net.SetupNextStep(net.y);
      net.OrderedGaussSeidelUpdate(net.y, yold, g2, 1e12, true, false, false, true, false, false);
      const Real cd = (K(G::kO1, IonRateTable::ICI)*ne + K(G::kO1, IonRateTable::ICT_ION)*nH*net.y(IH_plus))*net.y(IO_plus)
                      /(K(G::kO2, IonRateTable::IREC)*ne + K(G::kO2, IonRateTable::ICT_REC)*nH*g2.H);
      std::printf("  huge-step sweep: O++ %.6e, table C/D with the swept O+ %.6e; kO1 CI %.3e CTion %.3e kO2 rec %.3e CTrec %.3e\n",
                  net.y(IO_2plus), cd, K(G::kO1, IonRateTable::ICI), K(G::kO1, IonRateTable::ICT_ION),
                  K(G::kO2, IonRateTable::IREC), K(G::kO2, IonRateTable::ICT_REC));
      for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = yold[n];
    }
    std::printf("T %.0e: x_H+ %.4f  S0/S+/S++ %.3e %.3e %.3e (sum/xS %.6f)  O0/O+/O++ %.3e %.3e %.3e\n",
                T, xhp, xS0/1.45e-5, xS1/1.45e-5, xS2/1.45e-5, (xS0 + xS1 + xS2)/1.45e-5,
                g.O/3.2e-4, xO1/3.2e-4, xO2/3.2e-4);
    std::printf("  balance residuals: S0<->S+ %.2e  S+<->S++ %.2e  O+<->O++ %.2e\n",
                up_s0/dn_s1 - 1, up_s1/dn_s2 - 1, up_o1/dn_o2 - 1);
  }
}
