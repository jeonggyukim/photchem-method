// Band-averaged photoionization cross sections of the GOW17 photoionizing species
// over the three NCR bands, SB99 Z=0.014 2 Myr spectrum.
#include <cstdio>
#include <cstdlib>
#include <string>
#include "photchem/rates/photoion_band.hpp"
#include "../paths.hpp"
int main() {
  const std::string root = TigrisDir();
  Radiation::Spectrum sed; sed.LoadTable(root + "inputs/tables/sed/sb99_Z014_GenevaV00_2Myr.txt");
  PhotoIon::VernerXsec xs; xs.Load(root + "inputs/tables/rates/verner96_photx.dat");
  const Real edge[4] = {sed.LambdaMin(), 911.6, 1108.0, 2066.4};
  std::printf("%-5s %8s %8s | %12s %12s %12s | dhnu LyC [eV]\n", "ion", "E_th", "lam_th", "sigma LyC", "sigma LW", "sigma PE");
  for (int is = 0; is < NPHOTOION_GOW17; ++is) {
    std::printf("%-5s %8.3f %8.1f |", kPhotoIonName[is], SEDAvg::ThresholdEnergy(is, xs), SEDAvg::ThresholdWavelength(is, xs));
    for (int b = 0; b < 3; ++b) std::printf(" %12.4e", SEDAvg::BandAverageSigmaPhotoIon(sed, edge[b], edge[b+1], is, xs));
    std::printf(" | %.3f\n", SEDAvg::BandAverageDhnuPhotoIon(sed, edge[0], edge[1], is, xs));
  }
}
