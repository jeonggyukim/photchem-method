// ShieldingFactorCO (O(1) index, ln table) against the Athena++ fShield_CO_V09
// algorithm (linear search, logs per call) on the same table.
#include <cmath>
#include <cstdio>
#include <random>
#include "radiative_properties/shielding_co.hpp"
using namespace shielding_co;
int Idx(int len, const Real *x, Real v) {
  int i = 0;
  if (v < x[0]) return 0;
  if (v > x[len-1]) return len - 2;
  for (i = 0; v > x[i]; i++) {}
  return i - 1;
}
Real Lin(Real x0, Real x1, Real y0, Real y1, Real x) { return y0 + (y1 - y0)/(x1 - x0)*(x - x0); }
Real Ref(Real NCO, Real NH2) {
  const Real ns = 1e10; Real lco, lh2;
  if (NCO < ns && NH2 < ns) return 1.0;
  if (NCO < ns) { lco = 10.; lh2 = std::log10(NH2); }
  else if (NH2 < ns) { lh2 = 10.; lco = std::log10(NCO); }
  else { lco = std::log10(NCO); lh2 = std::log10(NH2); }
  int i0 = Idx(len_nco, log_nco, lco), j0 = Idx(len_nh2, log_nh2, lh2);
  Real f1 = Lin(log_nco[i0], log_nco[i0+1], std::log(theta[j0][i0]), std::log(theta[j0][i0+1]), lco);
  Real f2 = Lin(log_nco[i0], log_nco[i0+1], std::log(theta[j0+1][i0]), std::log(theta[j0+1][i0+1]), lco);
  return std::exp(Lin(log_nh2[j0], log_nh2[j0+1], f1, f2, lh2));
}
int main() {
  std::mt19937_64 g(1); std::uniform_real_distribution<Real> u(8.0, 24.0);
  Real worst = 0; int n = 0;
  auto test = [&](Real a, Real b) {
    const Real r = Ref(a, b), t = ShieldingFactorCO(a, b);
    worst = std::fmax(worst, std::fabs(t/r - 1)); ++n;
  };
  for (int k = 0; k < 1000000; ++k) test(std::pow(10.0, u(g)), std::pow(10.0, u(g)));
  for (int i = 1; i < len_nco; ++i) for (int j = 1; j < len_nh2; ++j)
    test(std::pow(10.0, log_nco[i]), std::pow(10.0, log_nh2[j]));
  std::printf("%d evaluations, max relative difference %.2e\n", n, worst);
  std::printf("examples: f(1e15,1e21) %.4f  f(1e17,1e21) %.4e  f(1e12,1e20) %.4f\n",
              ShieldingFactorCO(1e15, 1e21), ShieldingFactorCO(1e17, 1e21), ShieldingFactorCO(1e12, 1e20));
}
