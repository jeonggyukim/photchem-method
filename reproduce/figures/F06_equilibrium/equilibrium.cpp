// F06 (E8): thermal and chemical equilibrium of one cell with the GOW17 network and
// --photchem_ions=O3,S3,N3 against n_H, in an unshielded field chi = 1 (every photo-rate,
// H2 dissociation and the photoelectric field; S I photoionization 6e-10 s^-1 chi) and
// xi_cr = 2e-16 s^-1, as the uniform-field branch of PhotochemistryGOW17::SetCellRadiation
// sets them, and with the hot-gas cooling of the SNR runs. Each density starts at
// T0, atomic H, C+ and Si+, and is integrated with growing steps until the state
// changes by less than 1e-8 (relative) over a step, or t = 1e4 code.
// usage: equilibrium [T0] > out.txt
// Output per row: n_H [cm^-3], T [K], P/k [K cm^-3], t_end [code], x_e, then every
// species per H.
// Build: configure tigris-gow17 with -mpi -gow17 --photchem_ions=O3,S3,N3, then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$HOME/Projects/tigris-gow17/src
//   -I/opt/homebrew/opt/hdf5-mpi/include equilibrium.cpp
//   $HOME/Projects/tigris-gow17/src/photchem/network/gow17_thermo_table.cpp -o equilibrium
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
#include "photchem/sed_average.hpp"
using namespace gow17;
using G = GOW17Network;

int main(int argc, char **argv) {
  const Real T0 = (argc > 1) ? std::atof(argv[1]) : 1.0e4;
  const Real chi = 1.0, xi_cr = 2.0e-16;
  const std::string home = getenv("HOME");
  const std::string rd = home + "/Projects/tigris-gow17/inputs/tables/rates/";
  const std::string cd = home + "/Projects/tigris-gow17/inputs/tables/chianti_v11/";
  const std::string hot = home + "/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt";
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
  SemiImplicitSettings si{0.1, 1e-12, 2000, 0, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;

  std::printf("# equilibrium, chi = %g, xi_cr = %g s^-1, T0 = %g K; n_H T P/k t_end x_e",
              chi, xi_cr, T0);
  for (int n = 0; n < NSPEC; ++n) std::printf(" %s", std::string(G::species_names[n]).c_str());
  std::printf("\n");
  for (int in = 0; in <= 60; ++in) {
    const Real nH = std::pow(10.0, -2.0 + 0.1*in);
    Real rad[G::n_freq];
    for (int n = 0; n < G::n_ph; ++n) rad[n] = chi;
    rad[G::n_freq - 2] = chi;
    rad[G::n_freq - 1] = xi_cr;
    G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
    net.SetLyC(0.0, 0.0, 0.0, 0.0);
    net.SetHePhotoionization(0.0, 0.0);
    Real xi[NPHOTOION_GOW17] = {};
    xi[ISPI_SI] = 6.0e-10*chi;
    net.SetIonPhotoRates(xi);
    net.SetIonRates(&itab, &grid);
    net.SetHighIons(&htab);
    net.SetIonCooling(&ctab);
    net.SetIonizedCooling(false, 1.0);
    net.SetHotCooling(&htcie, 2.0e4, 3.5e4, 1.0, true);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0e-2;
    net.y(IC_plus) = xC;
    net.y(ISi_plus) = xSi;
    net.y(NSPEC) = 1.0;
    for (int it = 0; it < 4; ++it) {
      auto g = net.SetupNextStep(net.y);
      net.y(NSPEC) *= T0/net.Temperature(net.y, g);
    }
    Real t = 0.0, dt = 1.0e-6;
    const Real t_max = 1.0e4;
    while (t < t_max) {
      Real yold[NSPEC + 1];
      for (int n = 0; n <= NSPEC; ++n) yold[n] = net.y(n);
      SemiImplicit<G> sv(si, net, 0.0, dt);
      sv.SolveODE();
      t += dt;
      Real dmax = 0.0;
      for (int n = 0; n <= NSPEC; ++n) {
        const Real scale = std::fmax(std::fabs(yold[n]), 1.0e-10);
        dmax = std::fmax(dmax, std::fabs(net.y(n) - yold[n])/scale);
      }
      if (t > 1.0 && dmax < 1.0e-8) break;
      dt = std::fmin(dt*1.5, t_max - t);
    }
    auto g = net.SetupNextStep(net.y);
    const Real T = net.Temperature(net.y, g);
    const Real pok = (1.0 + 0.1 + g.e - net.y(gow17::IH2))*nH*T;
    std::printf("%.6e %.6e %.6e %.6e %.6e", nH, T, pok, t, g.e);
    for (int n = 0; n < NSPEC; ++n) std::printf(" %.6e", net.y(n));
    std::printf("\n");
  }
}
