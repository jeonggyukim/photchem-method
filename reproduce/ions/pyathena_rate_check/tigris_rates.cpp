// Raw Tigris ion rate coefficients [cm^3 s^-1] for the stages compared with pyathena.
// Columns keyed by the lower stage q of element Z:
//   rec   = RecombRate::Rate(Z, Z-q-1)        X^(q+1) + e -> X^q (RR+DR)
//   ci    = CollIonRate::Rate(Z, Z-q)         X^q + e -> X^(q+1)
//   ctrec = ChargeTransferRate::CtRec(Z, Z-q-1)  X^(q+1) + H0 -> X^q + H+
//   ction = ChargeTransferRate::CtIon(Z, Z-q)    X^q + H+ -> X^(q+1) + H0
#include <cstdio>
#include <cstdlib>
#include <string>
#include "photchem/rates/recomb_rate.hpp"
#include "photchem/rates/coll_ion_rate.hpp"
#include "photchem/rates/ct_rate.hpp"
int main() {
  const std::string d = std::string(getenv("HOME")) + "/Projects/tigris-gow17/inputs/tables/rates/";
  Rates::RecombRate rec; rec.Load(d + "badnell_rr_2023.dat", d + "badnell_dr_C_2023.dat",
                                  d + "badnell_dr_E_2023.dat");
  Rates::CollIonRate ci; ci.Load(d + "voronov97_coll_ion.dat");
  Rates::ChargeTransferRate ct; ct.LoadRecomb(d + "kingdon_ferland96_ct_rec.dat");
  ct.LoadIon(d + "kingdon_ferland96_ct_ion.dat");
  struct S { const char *name; int z, q; };
  const S st[] = {{"He0", 2, 0}, {"He1", 2, 1}, {"N0", 7, 0}, {"N1", 7, 1}, {"N2", 7, 2},
                  {"O0", 8, 0}, {"O1", 8, 1}, {"O2", 8, 2}, {"S0", 16, 0}, {"S1", 16, 1},
                  {"S2", 16, 2}, {"S3", 16, 3}};
  const double temps[] = {1e3, 3e3, 5e3, 8e3, 1e4, 2e4, 5e4, 1e5, 1e6};
  std::printf("ion,Z,q,T,rec,ci,ctrec,ction,hasDR\n");
  for (const S &s : st) {
    const int nrec = s.z - s.q - 1, nion = s.z - s.q;
    for (double T : temps) {
      std::printf("%s,%d,%d,%.6e,%.10e,%.10e,%.10e,%.10e,%d\n", s.name, s.z, s.q, T,
                  rec.Rate(s.z, nrec, T), ci.Rate(s.z, nion, T), ct.CtRec(s.z, nrec, T),
                  ct.CtIon(s.z, nion, T), (int)(nrec > 0 && rec.HasDR(s.z, nrec)));
    }
  }
  return 0;
}
