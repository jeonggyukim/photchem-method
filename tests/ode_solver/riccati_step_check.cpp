// GOW17Network::RiccatiStep_ against RK4 with fine steps on dx/dt = A + Bx + Cx^2.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include "photchem/network/gow17_network.hpp"
using gow17::GOW17Network;
double RK4(double x, double A, double B, double C, double h) {
  const double lam = std::sqrt(B*B - 4*A*C) + std::fabs(B) + 2*std::fabs(C)*(std::fabs(x) + 1);
  const int n = std::max(20000, static_cast<int>(40*lam*h));
  const double dt = h/n;
  auto f = [&](double y) { return A + B*y + C*y*y; };
  for (int i = 0; i < n; ++i) {
    double k1 = f(x), k2 = f(x + 0.5*dt*k1), k3 = f(x + 0.5*dt*k2), k4 = f(x + dt*k3);
    x += dt/6*(k1 + 2*k2 + 2*k3 + k4);
  }
  return x;
}
int main() {
  std::mt19937_64 g(7); std::uniform_real_distribution<double> u(0, 1);
  double worst = 0; int n = 0;
  for (int k = 0; k < 3000; ++k) {
    double A = std::pow(10.0, -6 + 6*u(g)), B = -10 + 20*u(g), C = -std::pow(10.0, -3 + 4*u(g));
    if (k % 10 == 0) C = 0.0;
    if (k % 10 == 1) A = 0.0;
    const double x0 = u(g), h = std::pow(10.0, -3 + 4*u(g));
    const double r = RK4(x0, A, B, C, h), t = GOW17Network::RiccatiStep_(x0, A, B, C, h);
    if (C == 0.0 && B > 0 && B*h > 30) continue;  // unbounded growth, not a GOW17 case
    worst = std::max(worst, std::fabs(t - r)/std::max(std::fabs(r), 1e-12)); ++n;
  }
  std::printf("%d cases, max relative difference vs RK4 %.2e\n", n, worst);
}
