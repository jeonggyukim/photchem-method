// Is the 200000-substep semi-implicit run a good reference? The F05 cases (A ionization,
// B recombination, n_H 100, 978 yr; C CO formation, n_H 1e3, cosmic rays only, 2.9 Myr),
// one zone, GOW17 core, integrated by CVODE BDF (SUNDIALS 7, dense direct solver, CVODE's
// own difference-quotient Jacobian) on the network's full right-hand side
// (evaluate_function: species and internal energy), at rtol 1e-8 and 1e-10, against the
// semi-implicit default (Gauss-Seidel + CO/HCO+ block) with 200000 equal substeps.
// Output per case: T and x (x_H+ for A and B, x_CO for C) of each, and the relative
// differences.
// Build: configure tigris-gow17 with -gow17 (defs.hpp only), then
//   OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -I$HOME/Projects/tigris-gow17/src
//   -I/opt/homebrew/opt/hdf5-mpi/include -I/opt/homebrew/opt/sundials/include
//   bdf_reference_check.cpp -L/opt/homebrew/opt/sundials/lib -lsundials_cvode
//   -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense
//   -lsundials_core -o bdf_reference_check
#include <cmath>
#include <cstdio>
#include <limits>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "photchem/network/gow17_network.hpp"
#include "photchem/network/gow17_semi_implicit.hpp"
using namespace gow17;
using G = GOW17Network;
struct Out { Real T, x; };
constexpr Real GAMMA = 5.0/3.0, TIME_CGS = 3.0856776e18/1.0e5, EDENS = 1.6738234e-24*1e10;

// The network with the case's field and initial state; rad must outlive it.
void Setup(int scenario, G *net) {
  const Real ev = 1.602176634e-12;
  const Real nH = (scenario == 2) ? 1.0e3 : 100.0;
  net->SetIonizedCooling(true, 1.0);
  net->SetLyC(scenario == 0 ? 1e-9 : 0.0, scenario == 0 ? 1.5e-9 : 0.0, 3.45*ev, 4.42*ev);
  for (int n = 0; n < NSPEC + 1; ++n) net->y(n) = 0.0;
  net->y(IC_plus) = 1.6e-4;
  net->y(ISi_plus) = 1.7e-6;
  net->y(IH_plus) = (scenario == 1) ? 0.99 : 1.0e-4;
  const Real T0 = (scenario == 1) ? 8000.0 : 100.0;
  auto g0 = net->SetupNextStep(net->y);
  net->y(NSPEC) = T0*Thermo::CvCold(0.0, 0.1, g0.e, GAMMA)*nH/EDENS;
}

GOW17Settings Settings() {
  GOW17Settings s;
  s.use_thermo_table = false; s.isothermal = false; s.isothermal_temperature = 0.0;
  s.zd = 1.0; s.xHe = 0.1; s.xC = 1.6e-4; s.xO = 3.2e-4; s.xSi = 1.7e-6;
  s.jacobian_hoist = false; s.temperature_min_rates = 1.0;
  const Real inf = std::numeric_limits<Real>::infinity();
  s.temperature_max_rates = inf; s.temperature_max_heating = inf;
  s.temperature_min_cooling = 1.0; s.temperature_max_cooling_nm = 1e9;
  s.Leff_CO_max = 3.0e20; s.H2_rovib_cooling = true; s.is_kgrH2_const = true;
  s.velocity_cgs = 1e5; s.length_cgs = 3.0856776e18; s.multi_d = true; s.three_d = true;
  return s;
}

Out Result(int scenario, G *net) {
  auto g = net->SetupNextStep(net->y);
  return {net->Temperature(net->y, g), (scenario == 2) ? net->y(ICO) : net->y(IH_plus)};
}

int Rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *user) {
  const G *net = static_cast<const G*>(user);
  View1D<Real> yv(N_VGetArrayPointer(y), G::neqs);
  RegisterArray<Real, G::neqs> f;
  net->evaluate_function(t, 0.0, yv, f);
  Real *fd = N_VGetArrayPointer(ydot);
  for (int i = 0; i < G::neqs; ++i) fd[i] = f(i);
  return 0;
}

Out RunBDF(int scenario, Real rtol) {
  const Real nH = (scenario == 2) ? 1.0e3 : 100.0;
  const Real t_end = (scenario == 2) ? 3.0 : 1.0e-3;
  Real rad[G::n_freq];
  const Real chi = (scenario == 2) ? 0.0 : 1.0;
  for (int n = 0; n < G::n_freq - 1; ++n) rad[n] = chi;
  rad[G::n_freq - 1] = 2.0e-16;
  G net(Settings(), nH, 0.0, View1D<Real>(rad, G::n_freq), GAMMA, TIME_CGS, EDENS);
  Setup(scenario, &net);
  SUNContext ctx;
  SUNContext_Create(SUN_COMM_NULL, &ctx);
  N_Vector y = N_VNew_Serial(G::neqs, ctx), atol = N_VNew_Serial(G::neqs, ctx);
  for (int i = 0; i < G::neqs; ++i) {
    NV_Ith_S(y, i) = net.y(i);
    // species per H: 1e-20, far below any abundance compared; energy: relative only
    NV_Ith_S(atol, i) = (i < NSPEC) ? 1.0e-20 : 1.0e-12*net.y(NSPEC);
  }
  void *mem = CVodeCreate(CV_BDF, ctx);
  CVodeInit(mem, Rhs, 0.0, y);
  CVodeSVtolerances(mem, rtol, atol);
  CVodeSetUserData(mem, &net);
  CVodeSetMaxNumSteps(mem, 100000000);
  SUNMatrix A = SUNDenseMatrix(G::neqs, G::neqs, ctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, ctx);
  CVodeSetLinearSolver(mem, LS, A);
  sunrealtype t = 0.0;
  const int flag = CVode(mem, t_end, y, &t, CV_NORMAL);
  long nst = 0;
  CVodeGetNumSteps(mem, &nst);
  for (int i = 0; i < G::neqs; ++i) net.y(i) = NV_Ith_S(y, i);
  std::printf("#   BDF rtol %.0e: flag %d, %ld steps\n", rtol, flag, nst);
  const Out o = Result(scenario, &net);
  CVodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(atol);
  SUNContext_Free(&ctx);
  return o;
}

Out RunSemiImplicit(int scenario, unsigned nsub) {
  const Real nH = (scenario == 2) ? 1.0e3 : 100.0;
  const Real t_end = (scenario == 2) ? 3.0 : 1.0e-3;
  Real rad[G::n_freq];
  const Real chi = (scenario == 2) ? 0.0 : 1.0;
  for (int n = 0; n < G::n_freq - 1; ++n) rad[n] = chi;
  rad[G::n_freq - 1] = 2.0e-16;
  G net(Settings(), nH, 0.0, View1D<Real>(rad, G::n_freq), GAMMA, TIME_CGS, EDENS);
  Setup(scenario, &net);
  SemiImplicitSettings si{0.1, 1e-12, 100000000, nsub, 1.0, 1, true, true, true, true, false,
                          false, true, false, false, true, 3, false, false};
  SemiImplicit<G> solver(si, net, 0.0, t_end);
  solver.SolveODE();
  return Result(scenario, &net);
}

int main() {
  const char *name[3] = {"A ionization", "B recombination", "C CO formation"};
  for (int sc = 0; sc < 3; ++sc) {
    std::printf("# case %d %s\n", sc, name[sc]);
    const Out si = RunSemiImplicit(sc, 200000);
    const Out b8 = RunBDF(sc, 1.0e-8), b10 = RunBDF(sc, 1.0e-10);
    std::printf("  semi-implicit 200000: T %.8e x %.8e\n", si.T, si.x);
    std::printf("  BDF rtol 1e-8       : T %.8e x %.8e\n", b8.T, b8.x);
    std::printf("  BDF rtol 1e-10      : T %.8e x %.8e\n", b10.T, b10.x);
    std::printf("  BDF 1e-8 vs 1e-10   : dT/T %+.2e dx/x %+.2e\n", b8.T/b10.T - 1, b8.x/b10.x - 1);
    std::printf("  semi-implicit vs BDF 1e-10: dT/T %+.2e dx/x %+.2e\n", si.T/b10.T - 1,
                si.x/b10.x - 1);
  }
}
