"""Convert the SB99 photon spectrum (Angstrom, photons/s/Angstrom) into a
Cloudy 'table SED' file (energy in Ryd, F_nu in arbitrary units), and check
that photon-number ratios above 24.59 eV and 35.12 eV to those above 13.6 eV
are preserved."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402

SRC = paths.tigris("inputs/tables/sed/sb99_Z014_GenevaV00_2Myr.txt")
OUT = "sb99_Z014_GenevaV00_2Myr.sed"
HC_EV_A = 12398.419843    # h c [eV Angstrom]
RYD_EV = 13.605693123     # Rydberg (infinite mass) [eV]
THRESH_EV = {"H": 13.598434, "He": 24.587389, "O+": 35.121125}

lam, xi = np.loadtxt(SRC, comments="#", unpack=True)
order = np.argsort(lam)
lam, xi = lam[order], xi[order]
assert np.all(np.diff(lam) > 0)

e_ryd = HC_EV_A / lam / RYD_EV
# F_nu = L_lambda lambda^2 / c and L_lambda = Xi h c / lambda, so F_nu ∝ Xi lambda
fnu = xi * lam
fnu = fnu / fnu.max()
fnu = np.maximum(fnu, 1e-30)

with open(OUT, "w") as f:
    f.write("# Starburst99 Z=0.014 Geneva v00, age 2 Myr, from %s\n" % SRC.split("/")[-1])
    f.write("# energy [Ryd]   F_nu [arbitrary]\n")
    for e, fn in zip(e_ryd[::-1], fnu[::-1]):
        f.write("%.8e  %.8e\n" % (e, fn))
    f.write("*****\n")

def photons_above_orig(e_ev):
    """Integrate Xi d(lambda) for lambda < hc/E with log-log interpolation."""
    lmax = HC_EV_A / e_ev
    g = np.geomspace(lam[0], lmax, 200001)
    y = np.exp(np.interp(np.log(g), np.log(lam), np.log(xi)))
    return np.trapezoid(y, g)

def photons_above_sed(e_ev):
    """Read the SED file back; integrate F_nu / (h nu) d nu above E."""
    d = np.loadtxt(OUT, comments=("#", "*"))
    e, fn = d[:, 0], d[:, 1]
    emin = e_ev / RYD_EV
    g = np.geomspace(emin, e[-1], 200001)
    y = np.exp(np.interp(np.log(g), np.log(e), np.log(fn)))
    return np.trapezoid(y / g, g)

ro = {k: photons_above_orig(v) for k, v in THRESH_EV.items()}
rs = {k: photons_above_sed(v) for k, v in THRESH_EV.items()}
for k in ("He", "O+"):
    a, b = ro[k] / ro["H"], rs[k] / rs["H"]
    print("Q(>%s)/Q(>H): original %.5e  sed %.5e  rel.diff %.2e" % (k, a, b, b / a - 1))
