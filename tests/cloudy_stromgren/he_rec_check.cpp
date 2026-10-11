// He+ -> He recombination: Badnell total (case A) vs GOW17 case B, 5000-20000 K.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "photchem/network/gow17_ion_table.hpp"
#include "../paths.hpp"
int main() {
  const std::string d = TigrisDir() + "inputs/tables/rates/";
  Rates::RecombRate rec; rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat",
                                  d + "badnell_dr_E_2023.dat");
  for (double T : {5000.0, 8000.0, 1.0e4, 2.0e4}) {
    const double L = std::log10(T);
    const double aB = 1.0e-11/std::sqrt(T)*(11.19 + (-1.676 + (-0.2852 + 0.04433*L)*L)*L);
    const double aA = rec.Rate(2, 1, T);
    std::printf("T=%6.0f  alpha_A(Badnell) %.3e  alpha_B(GOW17) %.3e  alpha_1 %.3e\n",
                T, aA, aB, aA - aB);
  }
}
