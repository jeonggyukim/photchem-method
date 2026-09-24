// Lambda(T) of the GOW17 + multi-ion network at CIE, from the network's own
// CoolingTerm, with its components, against NCR's hot-gas table (Gnat & Ferland 2012,
// Asplund 2009) and a CHIANTI all-element solar total.
// Every tracked species and higher-ion total is set to its CHIANTI CIE value at T; H+
// and He+ from CHIANTI ioneq_H/He; the electrons are the network's own ghost
// (species charges, the higher ions' mean charge and the hot-gas extra).
// Output: one row per T, Lambda_N = cooling per H / n_H [erg cm^3 s^-1]
// (cooling rate per volume / n_H^2).
// Build (against a source tree whose defs.hpp has GOW17_ENABLED 1, PHOTCHEM_IONS 1):
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I<src> -I/opt/homebrew/opt/hdf5-mpi/include
//     lambda_curve.cpp <src>/photchem/network/gow17_thermo_table.cpp -o lambda_curve
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;
using G = GOW17Network;

int main(int argc, char **argv) {
  const std::string home = getenv("HOME");
  const std::string rd = home + "/Projects/tigris-gow17/inputs/tables/rates/";
  const std::string cd = home + "/Projects/tigris-gow17/inputs/tables/chianti_v11/";
  const std::string hot = home + "/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt";
  const Real nH = (argc > 1) ? std::atof(argv[1]) : 1.0;
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
  // The higher-ion table as photchem_gow17.cpp builds it
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
  CIEElement cH(cd, "H", 1), cHe(cd, "He", 2);

  // CHIANTI all-element reference at solar (Asplund 2009) abundances
  struct Ab { const char *nm; int z; Real x; };
  const std::vector<Ab> solar = {{"H", 1, 1.0}, {"He", 2, 0.0851}, {"C", 6, 2.69e-4},
      {"N", 7, 6.76e-5}, {"O", 8, 4.90e-4}, {"Ne", 10, 8.51e-5}, {"Mg", 12, 3.98e-5},
      {"Si", 14, 3.24e-5}, {"S", 16, 1.32e-5}, {"Ar", 18, 2.51e-6}, {"Ca", 20, 2.19e-6},
      {"Fe", 26, 3.16e-5}};

  const Real gamma = 5.0/3.0, time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1e10;
  const Real xel[5] = {xC, xSi, xO, xS, xN};
  std::printf("# n_H = %g cm^-3. Lambda_N = cooling per H / n_H [erg cm^3 s^-1]\n", nH);
  std::printf("# logT x_e tigris_total H He tracked_ions higher_ions untracked_and_returned "
              "gow17_fits_rest sum_of_parts ncr_total ncr_H ncr_He ncr_metal chianti_solar\n");
  for (int i = 0; i <= 80; ++i) {
    const Real lT = 4.0 + 0.05*i, T = std::pow(10.0, lT);
    GOW17Settings s;
    s.use_thermo_table = false; s.isothermal = true; s.isothermal_temperature = T;
    s.zd = 1.0; s.xHe = 0.1; s.xC = xC; s.xO = xO; s.xSi = xSi; s.xS = xS; s.xN = xN;
    s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
    const Real inf = std::numeric_limits<Real>::infinity();
    s.temperature_max_rates = inf; s.temperature_max_heating = inf;
    s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
    s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
    s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
    Real rad[G::n_freq] = {};
    G net(s, nH, 0.0, View1D<Real>(rad, G::n_freq), gamma, time_cgs, edens);
    net.SetIonRates(&itab, &grid);
    net.SetHighIons(&htab);
    net.SetIonCooling(&ctab);
    net.SetIonizedCooling(false, 1.0);
    net.SetHotCooling(&htcie, 2.0e4, 3.5e4, 1.0, true);
    for (int n = 0; n < NSPEC + 1; ++n) net.y(n) = 0.0;
    std::vector<Real> xq, lq;
    cH.Stages(T, &xq, &lq);
    net.y(IH_plus) = xq[1];
    cHe.Stages(T, &xq, &lq);
    net.y(IHe_plus) = 0.1*xq[1];  // He++ is not a species; its electrons come from the hot extra
    struct Tr { int sp, q; };
    const std::vector<std::vector<Tr>> tracked = {
        {{IC_plus, 1}}, {{ISi_plus, 1}}, {{IO_plus, 1}, {IO_2plus, 2}},
        {{IS_plus, 1}, {IS_2plus, 2}}, {{IN_plus, 1}, {IN_2plus, 2}}};
    for (int k = 0; k < G::n_high; ++k) {
      const char *nm;
      int z, qmax;
      G::HighIonElement(k, &nm, &z, &qmax);
      CIEElement ce(cd, nm, z);
      ce.Stages(T, &xq, &lq);
      for (const Tr &t : tracked[k]) net.y(t.sp) = xel[k]*xq[t.q];
      Real xs = 0.0;
      for (int q = qmax + 1; q <= z; ++q) xs += xq[q];
      net.y(G::HighIonSpecies(k)) = xel[k]*xs;
    }
    net.y(NSPEC) = 1.0;
    auto g = net.SetupNextStep(net.y);
    const Real ne = nH*g.e;
    const Real total = net.CoolingTerm(net.y, g, T)/nH;
    // components, with the network's formulas
    const ThermoTable::Slot sl = grid.Locate(T);
    Real lam[G::n_ion_cool];
    ctab.LambdaAll(sl, ne, lam);
    Real trk = 0.0, hi = 0.0;
    for (int k = 0; k < G::n_ion_cool; ++k) trk += net.y(G::IonCoolingSpecies(k))*lam[k];
    for (int k = 0; k < G::n_high; ++k) {
      hi += net.y(G::HighIonSpecies(k))*htab.At(sl, k, HighIonTable::ILAMBDA);
    }
    trk *= ne/nH;
    hi *= ne/nH;
    const Real w2 = (T <= 2.0e4) ? 0.0 : (T >= 3.5e4 ? 1.0 :
        1.0/(1.0 + std::exp(-10.0*(T - 2.75e4)/(3.5e4 - 2.0e4))));
    const Real unt = w2*ne*htab.Untracked(sl)/nH;
    const Real lH = w2*htcie.At(HotCIETable::ILAMBDA_H, T);
    const Real lHe = htcie.At(HotCIETable::ILAMBDA_HE, T);
    const Real rest = total - trk - hi - unt - lH - lHe;
    const Real ncr_m = htcie.At(HotCIETable::ILAMBDA_METAL, T);
    // CHIANTI solar: x_e from CHIANTI too
    Real xe_c = 0.0, lam_c = 0.0;
    for (const Ab &a : solar) {
      CIEElement ce(cd, a.nm, a.z);
      ce.Stages(T, &xq, &lq);
      for (int q = 0; q <= a.z; ++q) {
        xe_c += a.x*q*xq[q];
        lam_c += a.x*xq[q]*lq[q];
      }
    }
    std::printf("%.2f %.4f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                lT, g.e, total, lH, lHe, trk, hi, unt, rest, lH + lHe + trk + hi + unt,
                htcie.At(HotCIETable::ILAMBDA_H, T) + lHe + ncr_m,
                htcie.At(HotCIETable::ILAMBDA_H, T), lHe, ncr_m, xe_c*lam_c);
  }
}
