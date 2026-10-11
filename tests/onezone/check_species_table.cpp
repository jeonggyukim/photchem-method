// Table-generated GOW17 ghosts and x_e vs AthenaK gow17.hpp:1802-1818 (transcribed).
// Build: c++ -std=c++17 -O2 -I $TIGRIS_DIR/src check_species_table.cpp
#include <cmath>
#include <cstdio>
#include <random>
#include "photchem/network/gow17_core.hpp"
namespace g = gow17;
using g::species_rows; using g::NSPEC; using g::IHe_plus; using g::IOHx; using g::ICHx; using g::ICO;
using g::IC_plus; using g::IHCO_plus; using g::IH_plus; using g::IH3_plus; using g::IH2_plus;
using g::IO_plus; using g::ISi_plus;
using photchem_net::ElementBudget; using photchem_net::ChargeSum;
int main() {
  const Real xHe = 0.1, xC = 1.6e-4, xO = 3.2e-4, xSi = 1.7e-6;
  std::mt19937 g(1); std::uniform_real_distribution<double> u(-12.0, 0.0);
  double worst = 0.0, worst_ghost = 0.0; int nclamp = 0;
  for (int t = 0; t < 100000; ++t) {
    Real y[NSPEC];
    for (int n = 0; n < NSPEC; ++n) y[n] = std::pow(10.0, u(g));
    y[g::IH2] *= 0.5; y[IHe_plus] *= xHe; y[ISi_plus] *= xSi;
    for (int n : {ICHx, IC_plus}) y[n] *= xC;
    for (int n : {IOHx, IO_plus}) y[n] *= xO;
    y[ICO] *= xC; y[IHCO_plus] *= xC;
    // AthenaK reference
    Real Si = std::fmax(xSi - y[ISi_plus], 0.0);
    Real C = std::fmax(xC - y[IHCO_plus] - y[ICHx] - y[ICO] - y[IC_plus], 0.0);
    Real O = std::fmax(xO - y[IHCO_plus] - y[IOHx] - y[ICO] - y[IO_plus], 0.0);
    Real He = std::fmax(xHe - y[IHe_plus], 0.0);
    Real e = y[IHe_plus] + y[IC_plus] + y[IHCO_plus] + y[IH3_plus] + y[IH2_plus]
           + y[IH_plus] + y[IO_plus] + y[ISi_plus];
    Real H = std::fmax(1.0 - (y[IOHx] + y[ICHx] + y[IHCO_plus] + 3.0*y[IH3_plus]
           + 2.0*y[IH2_plus] + y[IH_plus] + 2.0*y[g::IH2]), 0.0);
    if (H == 0.0 || C == 0.0 || O == 0.0) ++nclamp;
    const Real ref[6] = {H, He, C, O, Si, e};
    const Real tab[6] = {ElementBudget(species_rows, y, photchem_net::ELEM_H, 1.0),
                         ElementBudget(species_rows, y, photchem_net::ELEM_HE, xHe),
                         ElementBudget(species_rows, y, photchem_net::ELEM_C, xC),
                         ElementBudget(species_rows, y, photchem_net::ELEM_O, xO),
                         ElementBudget(species_rows, y, photchem_net::ELEM_SI, xSi),
                         ChargeSum(species_rows, y)};
    // Budgets: difference relative to the element total (the ghost itself can be a
    // near-cancellation of that total). x_e: relative to x_e.
    const Real scale[6] = {1.0, xHe, xC, xO, xSi, e};
    for (int k = 0; k < 6; ++k) {
      worst = std::fmax(worst, std::fabs(tab[k] - ref[k])/scale[k]);
      worst_ghost = std::fmax(worst_ghost, ref[k] > 0.0 ?
          std::fabs(tab[k] - ref[k])/ref[k] : std::fabs(tab[k]));
    }
  }
  std::printf("100000 random states, %d with a clamped ghost\n"
              "worst difference / element total (x_e: / x_e): %.3e\n"
              "worst difference / ghost value:               %.3e\n", nclamp, worst,
              worst_ghost);
  return worst < 1e-12 ? 0 : 1;
}
