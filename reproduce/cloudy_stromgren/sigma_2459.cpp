#include <cstdio>
#include "radiative_properties/photoion_xsec.hpp"
#include "../paths.hpp"
int main() {
  PhotoIon::VernerXsec x;
  x.Load(TigrisDir() + "inputs/tables/rates/verner96_photx.dat");
  const double e = x.Ethreshold(2, 2);
  std::printf("E_th(He I) %.4f eV  sigma_H %.4e  sigma_He %.4e  ratio %.3f\n",
              e, x.Sigma(1, 1, e*1.0001), x.Sigma(2, 2, e*1.0001),
              x.Sigma(2, 2, e*1.0001)/x.Sigma(1, 1, e*1.0001));
}
