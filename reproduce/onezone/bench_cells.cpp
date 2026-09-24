// One-zone GOW17 evolution with the semi-implicit solver, AthenaK (-DATHENAK) and the
// Tigris copy; 200 steps of 1e4 yr from an atomic state, Gauss-Seidel on and off.
// The two outputs must be identical.
#include <cstdio>
#include <limits>
#ifdef ATHENAK
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "chemistry/network/gow17.hpp"
#include "ode_solvers/semi_implicit.hpp"
using namespace chemistry; using ode_solvers::SemiImplicit;
using ode_solvers::SemiImplicitSettings;
#else
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
#endif

#include <chrono>
#include <random>
int main(int argc, char **argv) {
#ifdef ATHENAK
  Kokkos::initialize(argc, argv);
  {
#else
  (void)argc; (void)argv;
#endif
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = inf;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1.0; s.length_cgs = 1.0; s.multi_d = false; s.three_d = false;

  SemiImplicitSettings si;
  si.semi_implicit_cfl = 0.1; si.semi_implicit_yfloor = 1.0e-12;
  si.semi_implicit_n_substep_max = 100000; si.semi_implicit_n_substep_fixed = 0;
  si.semi_implicit_stiff_threshold = 1.0; si.semi_implicit_n_iter = 1;
  si.semi_implicit_energy = true; si.semi_implicit_renormalize = true;
  si.semi_implicit_hep_first = true; si.semi_implicit_exact_map = false;
  si.semi_implicit_exact_block = false; si.semi_implicit_co_block = true;
  si.semi_implicit_exact_ghosts = false; si.semi_implicit_h2_first = false;
  si.semi_implicit_adaptive = true; si.semi_implicit_nbad_max = 3;
  si.semi_implicit_refresh_rates = false; si.semi_implicit_table_deriv = false;

  const Real rad0[8] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0e-16};  // chi = 1
#ifdef ATHENAK
  DvceArray5D<Real> w0("w0", 1, 5, 1, 1, 3);
  DualArray1D<RegionSize> sizes("sizes", 1);
  sizes.h_view(0).dx1 = 1.0;
  sizes.template modify<HostMemSpace>(); sizes.template sync<DevMemSpace>();
  DvceArray1D<Real> ir("ir", 8);
  for (int n = 0; n < 8; ++n) ir(n) = rad0[n];
#else
  Real rad[8]; for (int n = 0; n < 8; ++n) rad[n] = rad0[n];
#endif
  const Real dt = 1.0e4*3.15576e7;
  si.semi_implicit_gauss_seidel = true;
  const int ncell = 20000;
  std::mt19937 g(7); std::uniform_real_distribution<double> u(0.0, 1.0);
  double checksum = 0.0; long nsub = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int c = 0; c < ncell; ++c) {
    const Real nH = std::pow(10.0, 4.0*u(g));
    const Real T0 = std::pow(10.0, 1.3 + 2.6*u(g));
    const Real fH2 = u(g);
    Real y[12] = {1e-8, 1e-9, 1e-9, 1.6e-4*fH2, 1.6e-4*(1.0 - fH2), 1e-11, 0.5*fH2,
                  1e-5, 1e-10, 1e-13, 1e-9, 1.7e-6};
#ifdef ATHENAK
    for (int i = 0; i < 3; ++i) w0(0, IDN, 0, 0, i) = nH;
    GOW17Network net(s, 0, 0, 0, 1, w0, sizes, ir, 1.0, 1.0, 5.0/3.0, 1.0, 1.0, 1.0);
#else
    GOW17Network net(s, nH, 0.0, View1D<Real>(rad, 8), 5.0/3.0, 1.0, 1.0);
#endif
    for (int n = 0; n < 12; ++n) net.y(n) = y[n];
    Real xe = 0.0; for (int n : {0, 4, 5, 7, 8, 9, 10, 11}) xe += y[n];
    net.y(12) = T0*Thermo::CvCold(y[6], 0.1, xe, 5.0/3.0)*nH;
    SemiImplicit<GOW17Network> solver(si, net, 0.0, dt);
    solver.SolveODE();
    nsub += solver.n_substeps;
    for (int n = 0; n < 13; ++n) checksum += net.y(n);
  }
  auto t1 = std::chrono::steady_clock::now();
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
  std::printf("cells %d  substeps/cell %.2f  ns/cell %.0f  ns/substep %.0f  checksum %.17g\n",
              ncell, double(nsub)/ncell, ns/ncell, ns/nsub, checksum);
#ifdef ATHENAK
  }
  Kokkos::finalize();
#endif
  return 0;
}
