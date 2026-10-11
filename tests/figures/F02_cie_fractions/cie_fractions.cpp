// F02: stage fractions of C, Si, O, S, N at fixed T with no radiation, from the GOW17
// network with --photchem_ions=O3,S3,N3 run to steady state, against CHIANTI v11 CIE.
// Writes one row per log T (4.0 to 7.0, 0.05 dex):
//   logT, then per element: tigris neutral, tracked stages..., X_high,
//                           chianti neutral, tracked stages..., summed above the top,
//   then tigris and chianti metal electrons per H.
// Build: configure tigris-gow17 with -mpi -gow17 --photchem_ions=O3,S3,N3, then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$TIGRIS_DIR/src
//   -I/opt/homebrew/opt/hdf5-mpi/include cie_fractions.cpp
//   $TIGRIS_DIR/src/photchem/network/gow17_thermo_table.cpp -o cie_fractions
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

int main() {
  const std::string rd = TigrisDir() + "inputs/tables/rates/";
  const std::string cd = TigrisDir() + "inputs/tables/chianti_v11/";
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
  std::vector<CIEElementSpec> els;
  for (int k = 0; k < G::n_high; ++k) {
    const char *nm;
    int z, q;
    G::HighIonElement(k, &nm, &z, &q);
    els.push_back({nm, z, q, 0.0});
  }
  HighIonTable htab;
  BuildHighIonTable(grid, cd, els, rec, ct, &htab);
  const Real xC = 1.6e-4, xSi = 1.7e-6, xO = 3.2e-4, xS = 1.45e-5, xN = 7.4e-5;
  const Real xel[5] = {xC, xSi, xO, xS, xN};
  const std::vector<std::vector<int>> tracked = {
      {IC_plus}, {ISi_plus}, {IO_plus, IO_2plus}, {IS_plus, IS_2plus}, {IN_plus, IN_2plus}};
  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real nH = 1.0;
  std::printf("# logT | per element C Si O S N: tigris x0 x1..top x_high | chianti x0 x1..top "
              "x_sum>top | xe_metal tigris chianti. Fractions of the element.\n");
  for (int it = 0; it <= 60; ++it) {
    const Real lT = 4.0 + 0.05*it, T = std::pow(10.0, lT);
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = true; s.isothermal_temperature = T;
    s.zd = 1.0; s.xHe = 0.1; s.xC = xC; s.xO = xO; s.xSi = xSi; s.xS = xS; s.xN = xN;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
    SemiImplicitSettings si{0.1, 1e-12, 100000000, 0, 1.0, 1, true, true, true, true, false,
                            false, true, false, false, true, 3, false, false};
    Real rad[G::n_freq] = {};
    G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
    net.SetIonRates(&itab, &grid);
    net.SetHighIons(&htab);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    net.y(IH_plus) = 1.0; net.y(IC_plus) = xC; net.y(ISi_plus) = xSi;
    net.y(NSPEC) = 1.0;
    for (int k = 0; k < 60; ++k) {
      SemiImplicit<G> sv(si, net, 0.0, 10.0);
      sv.SolveODE();
    }
    auto g = net.SetupNextStep(net.y);
    const Real neutral[5] = {g.C, g.Si, g.O, g.S, g.N};
    const ThermoTable::Slot sl = grid.Locate(T);
    Real xe_t = 0.0, xe_c = 0.0;
    std::printf("%.2f", lT);
    for (int k = 0; k < G::n_high; ++k) {
      const char *nm;
      int z, qmax;
      G::HighIonElement(k, &nm, &z, &qmax);
      CIEElement cie(cd, nm, z);
      std::vector<Real> xq, lq;
      cie.Stages(T, &xq, &lq);
      const Real xh = net.y(G::HighIonSpecies(k));
      std::printf(" %.5e", neutral[k]/xel[k]);
      int q = 1;
      for (int sp : tracked[k]) {
        std::printf(" %.5e", net.y(sp)/xel[k]);
        xe_t += q*net.y(sp);
        ++q;
      }
      std::printf(" %.5e", xh/xel[k]);
      xe_t += htab.At(sl, k, HighIonTable::IQ)*xh;
      Real xs = 0.0;
      for (int qq = qmax + 1; qq <= z; ++qq) xs += xq[qq];
      for (int qq = 0; qq <= qmax; ++qq) std::printf(" %.5e", xq[qq]);
      std::printf(" %.5e", xs);
      for (int qq = 1; qq <= z; ++qq) xe_c += qq*xq[qq]*xel[k];
    }
    std::printf(" %.5e %.5e\n", xe_t, xe_c);
  }
}
