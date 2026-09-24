"""Ionizing photon fractions per band from a -gow17 --gow17_bands=5 run's
athinput.runtime, against the ratios Cloudy's SED conversion reproduces.

usage: python check_photon_frac.py athinput.runtime
"""
import re
import sys

CLOUDY = {'Q(>35.12 eV)/Q(H)': 0.0109256, 'Q(>24.59 eV)/Q(H)': 0.141267}
txt = open(sys.argv[1]).read()
frac = {int(m.group(1)): float(m.group(2))
        for m in re.finditer(r'^photon_frac\[(\d)\]\s*=\s*([-+.\deE]+)', txt, re.M)}
tig = {'Q(>35.12 eV)/Q(H)': frac[0], 'Q(>24.59 eV)/Q(H)': frac[0] + frac[1]}
for k in CLOUDY:
    print('%-20s Tigris %.6f  Cloudy %.6f  rel %+.2e'
          % (k, tig[k], CLOUDY[k], tig[k]/CLOUDY[k] - 1))
print('sum of fractions %.12f' % sum(frac.values()))
