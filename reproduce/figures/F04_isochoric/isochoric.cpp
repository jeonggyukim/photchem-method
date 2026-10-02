// F04 (E1): isochoric cooling of one cell from CIE at 1e7 K with the GOW17 network and
// --photchem_ions=O3,S3,N3, no radiation, hot-gas cooling as in the SNR runs (CIE
// hydrogen, returned C, O, Si in the untracked column). Writes t [code], T [K], x_e and
// every species per H at each output, until T < 1e4 K or t = t_max.
// usage: isochoric N_H [T0] > out.txt
// Build: configure tigris-gow17 with -mpi -gow17 --photchem_ions=O3,S3,N3, then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$TIGRIS_DIR/src
//   -I/opt/homebrew/opt/hdf5-mpi/include isochoric.cpp
//   $TIGRIS_DIR/src/photchem/network/gow17_thermo_table.cpp -o isochoric
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
#include "../../paths.hpp"
using namespace gow17;
using G = GOW17Network;

int main(int argc, char **argv) {
  const Real nH = (argc > 1) ? std::atof(argv[1]) : 1.0;
  const Real T0 = (argc > 2) ? std::atof(argv[2]) : 1.0e7;
  const std::string rd = TigrisDir() + "inputs/tables/rates/";
  const std::string cd = TigrisDir() + "inputs/tables/chianti_v11/";
  const std::string hot = TigrisDir() + "inputs/tables/tigress_coolftn_ncr.txt";
  Rates::RecombRate rec;
  rec.Load(rd + "badnell_rr_2023.dat", rd + "badnell_dr_C_2023.dat", rd + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci;
  ci.Load(rd + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct;
  ct.LoadRecomb(rd + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(rd + "kingdon_ferland96_ct_ion.dat");
  ThermoTable grid;
  BuildThermoTable(grid);
  std::vector<IonStage> st;
  for (int k = 0; k < G::n_ion_stage; ++k) st.push_back(G::IonLadderStage(k));
  IonRateTable itab;
  BuildIonRateTable(grid, st, rec, ci, ct, &itab);
  const Real xC = 1.6e-4, xSi = 1.7e-6, xO = 3.2e-4, xS = 1.45e-5, xN = 7.4e-5;
  // the higher-ion table as photchem_gow17.cpp builds it
  std::vector<CIEElementSpec> els;
  for (int k = 0; k < G::n_high; ++k) {
    const char *nm;
    int z, q;
    G::HighIonElement(k, &nm, &z, &q);
    els.push_back({nm, z, q, 0.0});
  }
  els.push_back({"Ne", 10, -1, 8.51e-5});
  els.push_back({"Mg", 12, -1, 3.98e-5});
  els.push_back({"Ar", 18, -1, 2.51e-6});
  els.push_back({"Ca", 20, -1, 2.19e-6});
  els.push_back({"Fe", 26, -1, 3.16e-5});
  els.push_back({"C", 6, -1, 2.69e-4 - xC});
  els.push_back({"O", 8, -1, 4.90e-4 - xO});
  els.push_back({"Si", 14, -1, 3.24e-5 - xSi});
  HighIonTable htab;
  BuildHighIonTable(grid, cd, els, rec, ct, &htab);
  std::vector<std::string> cnames;
  for (int k = 0; k < G::n_ion_cool; ++k) cnames.push_back(G::IonCoolingName(k));
  IonCoolingTable ctab;
  BuildIonCoolingTable(grid, cd + "ion_cooling/", cnames, &ctab);
  HotCIETable htcie;
  htcie.Load(hot);

  GOW17Settings s;
  s.use_thermo_table = true; s.thermo_table = grid; s.isothermal = false;
  s.zd = 1.0; s.xHe = 0.1; s.xC = xC; s.xO = xO; s.xSi = xSi; s.xS = xS; s.xN = xN;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  // the default solver settings of photchem_gow17 (Gauss-Seidel, adaptive limiters)
  SemiImplicitSettings si{0.1, 1e-12, 2000, 0, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  Real rad[G::n_freq] = {};
  G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
  net.SetIonRates(&itab, &grid);
  net.SetHighIons(&htab);
  net.SetIonCooling(&ctab);
  net.SetIonizedCooling(false, 1.0);
  net.SetHotCooling(&htcie, 2.0e4, 3.5e4, 1.0, true);

  // CIE at T0: H+ and He+ from CHIANTI, the ladders' tracked stages and higher ions
  for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
  std::vector<Real> xq, lq;
  CIEElement(cd, "H", 1).Stages(T0, &xq, &lq);
  net.y(IH_plus) = xq[1];
  CIEElement(cd, "He", 2).Stages(T0, &xq, &lq);
  net.y(IHe_plus) = 0.1*xq[1];
  const Real xel[5] = {xC, xSi, xO, xS, xN};
  const std::vector<std::vector<int>> tracked = {
      {IC_plus}, {ISi_plus}, {IO_plus, IO_2plus}, {IS_plus, IS_2plus}, {IN_plus, IN_2plus}};
  for (int k = 0; k < G::n_high; ++k) {
    const char *nm;
    int z, qmax;
    G::HighIonElement(k, &nm, &z, &qmax);
    CIEElement(cd, nm, z).Stages(T0, &xq, &lq);
    int q = 1;
    for (int sp : tracked[k]) net.y(sp) = xel[k]*xq[q++];
    Real xs = 0.0;
    for (int qq = qmax + 1; qq <= z; ++qq) xs += xq[qq];
    net.y(G::HighIonSpecies(k)) = xel[k]*xs;
  }
  // energy for T0: T is linear in E at fixed composition
  net.y(NSPEC) = 1.0;
  for (int it = 0; it < 4; ++it) {
    auto g = net.SetupNextStep(net.y);
    net.y(NSPEC) *= T0/net.Temperature(net.y, g);
  }

  std::printf("# isochoric cooling, n_H = %g cm^-3, T0 = %g K; t [code time = 0.978 Myr] T x_e",
              nH, T0);
  for (int n = 0; n < NSPEC; ++n) std::printf(" %s", std::string(G::species_names[n]).c_str());
  std::printf("\n");
  // output every 0.02 dex in T, integrating in chunks of a tenth of the cooling time
  Real t = 0.0, last_logT = 99.0;
  const Real t_max = 1.0e4;
  while (t < t_max) {
    auto g = net.SetupNextStep(net.y);
    const Real T = net.Temperature(net.y, g);
    if (std::fabs(std::log10(T) - last_logT) >= 0.02 || t == 0.0) {
      std::printf("%.6e %.6e %.6e", t, T, g.e);
      for (int n = 0; n < NSPEC; ++n) std::printf(" %.6e", net.y(n));
      std::printf("\n");
      last_logT = std::log10(T);
    }
    if (T < 1.0e4) break;
    const Real edot = net.Edot(net.y, g);
    const Real tcool = (edot != 0.0) ? std::fabs(net.y(NSPEC)/edot) : 1.0;
    const Real dt = std::fmin(0.02*tcool, t_max - t);
    SemiImplicit<G> sv(si, net, 0.0, dt);
    sv.SolveODE();
    t += dt;
  }
}
