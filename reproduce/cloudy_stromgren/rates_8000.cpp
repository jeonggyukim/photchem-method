// S and O ladder rate coefficients at 8000 K from the Tigris ion-rate evaluators.
#include <cstdio>
#include <cstdlib>
#include <string>
#include "photchem/network/gow17_ion_table.hpp"
int main() {
  const std::string d = std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/rates/";
  Rates::RecombRate rec; rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat",
                                  d + "badnell_dr_E_2023.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  const double T = 8000.0;
  const int zs[2] = {8, 16};
  for (int z : zs)
    for (int q = 1; q <= 3; ++q) {
      const int n = z - q;
      std::printf("Z=%2d q=%d  alpha(RR+DR) %.3e  CT_rec(+H) %.3e  | q-1=%d CT_ion(+H+) %.3e\n",
                  z, q, rec.Rate(z, n, T), ct.CtRec(z, n, T), q - 1, ct.CtIon(z, n + 1, T));
    }
}
