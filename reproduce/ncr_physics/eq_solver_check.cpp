// GOW17 equilibrium by NCR's method (ncr_solver.hpp SolveToEquilibrium): repeat
// [temperature step, then one Gauss-Seidel backward-Euler species step] with a long
// step, halving the step when T would change by more than 20 per cent, until H, H+,
// H2 and CO change by less than tol and heating = cooling to tol. chem_only holds T
// fixed. Checked on (1) the hii cell against a long time integration, (2) the B1
// state against Athena++, (3) chem_only at the temperature of (1).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include "photchem/network/gow17_network.hpp"
using namespace gow17;

struct EqResult { int niter; int nhalve; Real balance; bool converged; };

EqResult SolveToEquilibrium(GOW17Network &net, const Real dt, const int niter_max,
                            const bool chem_only, const Real T_fixed, const Real nH,
                            const Real xHe, const Real gamma, const Real edens_cgs,
                            const Real tol) {
  constexpr int iie = GOW17Network::neqs - 1;
  constexpr int iH2 = 6, iHp = 7, iCO = 3;
  EqResult res{0, 0, 1.0, false};
  Real y_old[GOW17Network::neqs];
  for (int it = 1; it <= niter_max; ++it) {
    res.niter = it;
    auto g = net.SetupNextStep(net.y);
    const Real xH_prev = g.H, xHp_prev = net.y(iHp), xH2_prev = net.y(iH2),
               xCO_prev = net.y(iCO);
    Real h = dt;
    if (!chem_only) {
      // Semi-implicit energy step, Kim+23 eq. 59; rejected and halved beyond 20%.
      const Real E = net.y(iie);
      const Real edot = net.Edot(net.y, g);
      const Real dE = net.EnergyPerturbationScale()*std::abs(E);
      net.y(iie) = E + dE;
      const Real dedot_de = (net.Edot(net.y, g) - edot)/(net.y(iie) - E);
      net.y(iie) = E;
      for (;;) {
        const Real E_new = E + edot*h/(1.0 - h*dedot_de);
        if (std::isfinite(E_new) && E_new > 0.8*E && E_new < 1.2*E) {
          net.y(iie) = E_new;
          break;
        }
        h *= 0.5;
        ++res.nhalve;
      }
      g = net.SetupNextStep(net.y);
    }
    for (int n = 0; n < GOW17Network::neqs; ++n) y_old[n] = net.y(n);
    net.OrderedGaussSeidelUpdate(net.y, y_old, g, h, true, false, false, true, false,
                                 false);
    net.RenormalizeElements(net.y);
    g = net.SetupNextStep(net.y);
    if (chem_only) {
      net.y(iie) = T_fixed*Thermo::CvCold(net.y(iH2), xHe, g.e, gamma)*nH/edens_cgs;
      g = net.SetupNextStep(net.y);
    }
    const Real T = net.Temperature(net.y, g);
    const Real heat = net.HeatingTerm(net.y, g, T), cool = net.CoolingTerm(net.y, g, T);
    res.balance = std::abs(heat - cool)/(heat + cool);
    const Real f = std::max({std::abs(g.H/xH_prev - 1), std::abs(net.y(iHp)/xHp_prev - 1),
                             std::abs(net.y(iH2)/xH2_prev - 1),
                             std::abs(net.y(iCO)/xCO_prev - 1)});
    if (f < tol && (chem_only || res.balance < tol)) { res.converged = true; break; }
  }
  return res;
}

GOW17Settings Settings(bool kgr_const) {
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = kgr_const;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  return s;
}

void Report(const char *label, GOW17Network &net, const EqResult &r) {
  auto g = net.SetupNextStep(net.y);
  std::printf("%s: %s, %d iterations, %d halvings, |H-C|/(H+C) %.1e\n  T %.4f K "
              "x_H2 %.6f x_H+ %.4e x_CO %.4e x_C+ %.4e x_e %.4e x_HCO+ %.4e x_Si+ %.4e\n",
              label, r.converged ? "converged" : "NOT converged", r.niter, r.nhalve,
              r.balance, net.Temperature(net.y, g), net.y(6), net.y(7), net.y(3),
              net.y(4), g.e, net.y(5), net.y(11));
}

int main() {
  const Real gamma = 5.0/3.0, tol = 1e-6;
  // (1) hii cell, Tigris ism units
  {
    const Real time_cgs = 3.0856776e18/1.0e5, edens = 1.6738234e-24*1.0e10, nH = 100.0;
    GOW17Settings s = Settings(false);
    Real rad[8] = {1, 1, 1, 1, 0, 1, 1, 2e-16};
    for (Real dt : {1e5, 1e2}) {
      GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
      net.SetIonizedCooling(true, 1.0);
      for (int n = 0; n < 13; ++n) net.y(n) = 0.0;
      net.y(7) = 0.01; net.y(4) = 1.6e-4; net.y(11) = 1.7e-6;
      net.y(12) = 20.0*Thermo::CvCold(0.0, 0.1, 0.01 + 1.6e-4 + 1.7e-6, gamma)*nH/edens;
      EqResult r = SolveToEquilibrium(net, dt, 2000, false, 0.0, nH, 0.1, gamma, edens, tol);
      char lab[64]; std::snprintf(lab, sizeof(lab), "(1) hii cell, dt %.0e", dt);
      Report(lab, net, r);
      if (dt == 1e5) {
        // (3) chem_only at this T, from the same initial species
        auto g = net.SetupNextStep(net.y);
        const Real T1 = net.Temperature(net.y, g);
        GOW17Network n3(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
        n3.SetIonizedCooling(true, 1.0);
        for (int n = 0; n < 13; ++n) n3.y(n) = 0.0;
        n3.y(7) = 0.01; n3.y(4) = 1.6e-4; n3.y(11) = 1.7e-6;
        n3.y(12) = T1*Thermo::CvCold(0.0, 0.1, 0.01 + 1.6e-4 + 1.7e-6, gamma)*nH/edens;
        EqResult r3 = SolveToEquilibrium(n3, dt, 2000, true, T1, nH, 0.1, gamma, edens, tol);
        Report("(3) chem_only at the T of (1)", n3, r3);
      }
    }
    std::printf("  long time integration (earlier): T 94.04 K x_H2 0.449 x_H+ 9.87e-06 "
                "x_CO 8.26e-09 (noisy, heating/cooling within 0.3-1.6%%)\n");
  }
  // (2) B1 state, AthenaK units, against the Athena++ reference
  {
    const Real length = 3.0856775809623245e18, mass = 6.882615081124725e31,
               time_cgs = 3.08567758e13;
    const Real dens = mass/(length*length*length), vel = length/time_cgs;
    const Real edens = dens*vel*vel, nH = 109.21;
    GOW17Settings s = Settings(false);
    s.velocity_cgs = vel; s.length_cgs = length;
    Real rad[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 2e-16};
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), gamma, time_cgs, edens);
    for (int n = 0; n < 12; ++n) net.y(n) = 1e-6;
    net.y(12) = nH/(gamma - 1.0);
    EqResult r = SolveToEquilibrium(net, 1e5, 2000, false, 0.0, nH, 0.1, gamma, edens, tol);
    Report("(2) B1 state", net, r);
    const Real e_ref = 28.234391212463375;
    std::printf("  Athena++: e %.4f (this %.4f, %+.2f%%) x_H2 0.449894 x_H+ 2.943e-06 "
                "x_CO 1.3346e-04 x_C+ 3.498e-07 x_HCO+ 6.033e-08 x_Si+ 6.626e-07\n",
                e_ref, net.y(12), 100*(net.y(12)/e_ref - 1));
  }
  return 0;
}
