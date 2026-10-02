"""Ion fractions in the radiative SNR (GOW17 + O2,S3,N2) against CHIANTI CIE at the
cell temperature, mass-weighted in bins of log T, for one snapshot.

usage: python shell_ions_vs_cie.py SNAPSHOT.athdf [OUT.png]
The top tracked stage is compared with CHIANTI's stages at and above it summed,
since the network has no ion above it (the CIE pool return is not wired).
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
athena_read = paths.athena_read()

# ism units: density unit m_H/cm^3, rho = 1.4 n_H; pressure unit m_H (km/s)^2
MH, KB, MU_H = 1.6735575e-24, 1.380649e-16, 1.4
X_HE = 0.1
TOT = {'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}
STAGES = {'O': ['O+', 'O++'], 'S': ['S+', 'S++', 'S3+'], 'N': ['N+', 'N++']}
CHIANTI = paths.tigris('inputs/tables/chianti_v11', 'ioneq_%s.txt')

d = athena_read.athdf(sys.argv[1])
names = CORE + IONS
# Feedback's three scalars (rmetal, rSN, rret) come first; photochemistry starts at 3.
i0 = 3
keys = ['rHI'] + ['r%d' % n for n in range(i0 + 1, i0 + len(names))]
x = {nm: d[k].astype(np.float64).ravel() for nm, k in zip(names, keys)}
nH = d['rho'].astype(np.float64).ravel()/MU_H
xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
T = d['press'].ravel()*MH*1e10/KB/(nH*(1.0 - x['H2'] + X_HE + xe))
print('%s t = %.4f' % (sys.argv[1].split('/')[-1], d['Time']))

edges = np.arange(4.0, 7.01, 0.25)
lt = np.log10(np.maximum(T, 1.0))
for el, st in STAGES.items():
    tab = np.loadtxt(CHIANTI % el)
    xs = {q: x[q]/TOT[el] for q in st}
    x0 = 1.0 - sum(xs.values())
    print('%s: log T bin, mass frac of %s, tigris [0, %s] / CHIANTI [0, %s (top: summed)]'
          % (el, el, ', '.join(st), ', '.join(st)))
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (lt >= lo) & (lt < hi)
        if m.sum() < 5:
            continue
        w = nH[m]
        tig = [np.sum(x0[m]*w)/w.sum()] + [np.sum(xs[q][m]*w)/w.sum() for q in st]
        # CHIANTI at each cell's T, then mass-weighted
        c = np.array([np.interp(lt[m], tab[:, 0], tab[:, 1 + k]) for k in range(tab.shape[1] - 1)])
        top = len(st)
        cie = [np.sum(c[k]*w)/w.sum() for k in range(top)] + [np.sum(c[top:].sum(0)*w)/w.sum()]
        print('  %.2f-%.2f n=%6d  ' % (lo, hi, m.sum())
              + ' '.join('%.3f/%.3f' % (a, b) for a, b in zip(tig, cie)))
