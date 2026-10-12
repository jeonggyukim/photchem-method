"""Write a blackbody photon spectrum for both codes: bb_T<T>.txt for Tigris
(<photchem_gow17>/sed_file: wavelength [Angstrom], photons s^-1 Angstrom^-1, any
normalization) and bb_T<T>.sed for Cloudy ('table SED': energy [Ryd], F_nu, any
normalization). Both sample the same Planck function on the same wavelengths, from
40 Angstrom (310 eV; at 1e5 K the photon spectrum there is 1e-15 of its peak) to
1.6e6 Angstrom. Prints the photon fractions above the H, He I, O II and He II edges.
Usage: python make_bb_sed.py T [OUTDIR]   (T as written in the file names, e.g. 5e4)"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TSTR = sys.argv[1] if len(sys.argv) > 1 else '1e5'
T = float(TSTR)
TAG = 'bb_T' + TSTR
if len(sys.argv) > 2:
    HERE = os.path.join(HERE, sys.argv[2])
HC_EV_A = 12398.419843    # h c [eV Angstrom]
RYD_EV = 13.605693123     # Rydberg [eV]
KT_EV = 8.617333262e-5*T  # k_B T [eV]
EDGES = {'H': 13.598434, 'He I': 24.587389, 'O II': 35.121125, 'He II': 54.417763}

lam = np.geomspace(40.0, 1.6e6, 4000)
x = HC_EV_A/lam/KT_EV
# photons per unit wavelength of a Planck spectrum: B_lambda/(h c/lambda) ~ lambda^-4
xi = lam**-4/np.expm1(x)
xi /= xi.max()

with open(os.path.join(HERE, TAG + '.txt'), 'w') as f:
    f.write('# Photon spectrum of a blackbody at T = %.4g K, for band averaging.\n' % T)
    f.write('# Written by make_bb_sed.py; the same samples as %s.sed.\n#\n' % TAG)
    f.write('# lambda [Angstrom]  photons s^-1 Angstrom^-1 (arbitrary normalization)\n')
    for lv, v in zip(lam, np.maximum(xi, 1e-300)):
        f.write('%14.4f  %.8e\n' % (lv, v))

# F_nu = L_lambda lambda^2 / c and L_lambda = Xi h c / lambda, so F_nu ~ Xi lambda
fnu = xi*lam
fnu = np.maximum(fnu/fnu.max(), 1e-30)
e_ryd = HC_EV_A/lam/RYD_EV
with open(os.path.join(HERE, TAG + '.sed'), 'w') as f:
    f.write('# blackbody T = %.4g K, written by make_bb_sed.py\n' % T)
    f.write('# energy [Ryd]   F_nu [arbitrary]\n')
    for e, v in zip(e_ryd[::-1], fnu[::-1]):
        f.write('%.8e  %.8e\n' % (e, v))
    f.write('*****\n')

e_ev = HC_EV_A/lam
q = {k: np.trapezoid(np.where(e_ev >= v, xi, 0.0)[::-1], lam[::-1])*-1
     for k, v in EDGES.items()}
for k in EDGES:
    print('Q(>%s)/Q(H) = %.4e' % (k, q[k]/q['H']))
print('mean photon energy above 13.6 eV: %.2f eV'
      % (np.trapezoid(np.where(e_ev >= 13.6, xi*e_ev, 0)[::-1], lam[::-1])
         / np.trapezoid(np.where(e_ev >= 13.6, xi, 0)[::-1], lam[::-1])))
