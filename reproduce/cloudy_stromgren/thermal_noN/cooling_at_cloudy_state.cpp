// Tigris GOW17 cooling [erg cm^-3 s^-1] at Cloudy zone states read from stdin:
// "r T ne xHp xHep xOp xO2p xSp xS2p xS3p Ctot_cloudy" per line (fractions of each
// element; C taken as C+). Ion set O2,S3, per-ion cooling on, nebular lump off.
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;
int main() {
  const std::string d = "/Users/jgkim/Projects/tigris-gow17/inputs/tables/";
  ThermoTable grid; BuildThermoTable(grid);
  Rates::RecombRate rec; rec.Load(d + "rates/badnell_rr_2023.dat",
      d + "rates/badnell_dr_C_2023.dat", d + "rates/badnell_dr_E_2023.dat");
  Rates::CollIonRate ci; ci.Load(d + "rates/voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(d + "rates/kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "rates/kingdon_ferland96_ct_ion.dat");
  std::vector<IonStage> st;
  for (int k = 0; k < GOW17Network::n_ion_stage; ++k) st.push_back(GOW17Network::IonLadderStage(k));
  IonRateTable itab; BuildIonRateTable(grid, st, rec, ci, ct, &itab);
  std::vector<std::string> nm;
  for (int k = 0; k < GOW17Network::n_ion_cool; ++k) nm.push_back(GOW17Network::IonCoolingName(k));
  IonCoolingTable ctab; BuildIonCoolingTable(grid, d + "chianti_v11/ion_cooling/", nm, &ctab);
  GOW17Settings s;
  s.use_thermo_table = true; s.thermo_table = grid; s.isothermal = false;
  s.isothermal_temperature = 0.0; s.zd = 0.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4;
  s.xSi = 1.7e-6; s.xS = 1.45e-5; s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = false;
  s.velocity_cgs = 1.0; s.length_cgs = 1.0; s.multi_d = false; s.three_d = false;
  Real rad[GOW17Network::n_freq] = {};
  double r, T, ne, xHp, xHep, xOp, xO2p, xSp, xS2p, xS3p, Cc;
  const double nH = 100.0;
  std::printf("# r[pc] T ne  C_tigris C_cloudy ratio | C_ion(O,S) C_Hion C_other [erg cm^-3 s^-1]\n");
  while (std::cin >> r >> T >> ne >> xHp >> xHep >> xOp >> xO2p >> xSp >> xS2p >> xS3p >> Cc) {
    Real out[3];
    for (int pass = 0; pass < 3; ++pass) {
      GOW17Network net(s, nH, 0.0, View1D<Real>(rad, GOW17Network::n_freq), 5.0/3.0, 1.0, 1.0);
      if (pass >= 1) net.SetIonizedCooling(false, 1.0);
      if (pass >= 2) { net.SetIonRates(&itab, &grid); net.SetIonCooling(&ctab); }
      for (int n = 0; n < GOW17Network::neqs; ++n) net.y(n) = 0.0;
      net.y(IHe_plus) = 0.1*xHep; net.y(IH_plus) = xHp; net.y(IC_plus) = 1.6e-4;
      net.y(IO_plus) = 3.2e-4*xOp; net.y(IO_2plus) = 3.2e-4*xO2p; net.y(ISi_plus) = 1.7e-6;
      net.y(IS_plus) = 1.45e-5*xSp; net.y(IS_2plus) = 1.45e-5*xS2p; net.y(IS_3plus) = 1.45e-5*xS3p;
      const Real xe = xHp + 0.1*xHep + 1.6e-4 + 3.2e-4*(xOp + 2*xO2p) + 1.7e-6
                      + 1.45e-5*(xSp + 2*xS2p + 3*xS3p);
      net.y(GOW17Network::neqs - 1) = T*Thermo::CvCold(0.0, 0.1, xe, 5.0/3.0)*nH;
      const auto g = net.SetupNextStep(net.y);
      out[pass] = net.CoolingTerm(net.y, g, T) * nH;
    }
    std::printf("%.3f %.0f %.1f  %.3e %.3e %.3f | %.3e %.3e %.3e\n", r, T, ne, out[2], Cc,
                out[2]/Cc, out[2]-out[1], out[1]-out[0], out[0]);
  }
}
