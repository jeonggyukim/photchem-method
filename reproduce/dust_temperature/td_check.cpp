// Closed-form T_d (radiative_properties/dust_temperature.hpp) against the equilibrium
// temperature of the Zubko+04 BARE-GR-S effective grain (effective_grain.hpp), which
// inverts 4 pi int C_abs B_lam(T) dlam.  Background off (T_bg -> 0) for the comparison.
//   c++ -std=c++17 -O2 -I<tigris>/src td_check.cpp -o td_check && ./td_check <model file>
#include <cmath>
#include <cstdio>
#include <string>

#include "radiative_properties/dust_temperature.hpp"
#include "radiative_properties/effective_grain.hpp"

int main(int argc, char **argv) {
  Emission::EffectiveGrain grain;
  grain.Load(argv[1]);
  double worst = 0.0;
  std::printf("# H [erg/s/H]   T_grain [K]   T_closed [K]   rel diff\n");
  for (double t = 3.0; t <= 30.01; t *= 1.25) {
    // absorbed power for which the grain sits at t: invert by bisection on the table
    double lo = 1e-40, hi = 1e-10;
    for (int it = 0; it < 200; ++it) {
      const double mid = std::sqrt(lo*hi);
      (grain.EquilibriumTemperature(mid*1e-7) < t ? lo : hi) = mid;
    }
    const double h = std::sqrt(lo*hi);
    const double tc = DustTemperature::Equilibrium(h, 1.0, 0.0);
    const double d = tc/t - 1.0;
    worst = std::fmax(worst, std::fabs(d));
    std::printf("%12.4e   %8.3f   %8.3f   %+.4f\n", h, t, tc, d);
  }
  std::printf("# largest |rel diff| over 3-30 K: %.4f\n", worst);
  std::printf("# full ISRF, 5.73e-24 erg/s/H: grain %.2f K, closed %.2f K (with CMB %.2f K)\n",
              grain.EquilibriumTemperature(5.73e-24*1e-7),
              DustTemperature::Equilibrium(5.73e-24, 1.0, 0.0),
              DustTemperature::Equilibrium(5.73e-24, 1.0));
  return 0;
}
