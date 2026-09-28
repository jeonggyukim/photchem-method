// The gas-dust balance (DustTemperature::EquilibriumWithGas) against
// the dust balance H + A (T - T_d) = C Z_d (T_d^6 - T_bg^6) solved by bisection.
//   c++ -std=c++17 -O2 -I<tigris>/src td_gas_check.cpp -o td_gas_check && ./td_gas_check
#include <cmath>
#include <cstdio>

#include "radiative_properties/dust_temperature.hpp"

int main() {
  namespace D = DustTemperature;
  const double c = D::EmissionCoefficient(), bg = D::kTempCMB;
  double worst = 0.0;
  std::printf("# H[erg/s/H]  n_H     T     Td_rad   Td_newton  Td_exact  rel\n");
  for (double h : {1e-27, 5.73e-24}) {
    for (double nh : {1e3, 1e5, 1e6, 1e7, 1e8}) {
      for (double t : {6.0, 10.0, 30.0, 100.0}) {
        const double a = D::AlphaGasDustHM89(t);
        const double tn = D::EquilibriumWithGas(h, 1.0, nh, t, a);
        const double aa = a*nh*std::sqrt(t);
        double lo = 1.0, hi = 300.0;
        for (int it = 0; it < 200; ++it) {
          const double m = 0.5*(lo + hi);
          const double f = c*(std::pow(m, 6) - std::pow(bg, 6)) - h - aa*(t - m);
          (f > 0 ? hi : lo) = m;
        }
        const double te = 0.5*(lo + hi), r = tn/te - 1.0;
        worst = std::fmax(worst, std::fabs(r));
        std::printf("%9.2e %8.0e %6.1f %8.3f %9.3f %9.3f %+.4f\n", h, nh, t,
                    D::Equilibrium(h, 1.0), tn, te, r);
      }
    }
  }
  std::printf("# largest |rel|: %.4f\n", worst);
  return 0;
}
