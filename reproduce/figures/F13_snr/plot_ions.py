"""F13b: ion fractions against temperature in the radiative SNR with GOW17 + O3,S3,N3,
before and after shell formation, against CHIANTI v11 CIE at the same T. Per element, the
tracked stages and X_high (all stages above the top tracked one; CIE summed the same way).
Points: n_H-weighted mean over cells in 0.1 dex bins of T (bins with >= 5 cells).
Writes ../../../figures/F13_snr_ions.{pdf,png}.

Inputs: run_ions/snr.out2.00006 (t = 0.030) and .00009 (t = 0.045) of run_series.sh."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F13_snr_ions')
RUN = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion/M6_rad_snr/'
                         'F13_series/run_ions')
CHIANTI = '/Users/jgkim/Projects/tigris-gow17/inputs/tables/chianti_v11/ioneq_%s.txt'
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
# element, gas-phase total per H, tracked ions, higher ions, stage labels
ELEMENTS = [('C', 1.6e-4, ['C+'], 'C_high', ['C I', 'C II', r'C$^{\geq 2+}$']),
            ('Si', 1.7e-6, ['Si+'], 'Si_high', ['Si I', 'Si II', r'Si$^{\geq 2+}$']),
            ('O', 3.2e-4, ['O+', 'O++'], 'O_high', ['O I', 'O II', 'O III', r'O$^{\geq 3+}$']),
            ('S', 1.45e-5, ['S+', 'S++'], 'S_high', ['S I', 'S II', 'S III', r'S$^{\geq 3+}$']),
            ('N', 7.4e-5, ['N+', 'N++'], 'N_high', ['N I', 'N II', 'N III', r'N$^{\geq 3+}$'])]
MH, KB = 1.6735575e-24, 1.380649e-16


def binned(fn):
    d = athena_read.athdf(fn)
    keys = ['rHI'] + ['r%d' % n for n in range(4, 3 + len(NAMES))]
    x = {nm: d[k].astype(np.float64).ravel() for nm, k in zip(NAMES, keys)}
    nh = d['rho'].astype(np.float64).ravel()/1.4
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    T = d['press'].ravel()*MH*1e10/KB/(nh*(1.0 - x['H2'] + 0.1 + xe))
    lt = np.log10(np.maximum(T, 1.0))
    edges = np.arange(4.0, 7.51, 0.1)
    mid = 0.5*(edges[1:] + edges[:-1])
    out = {}
    for el, xtot, tracked, high, _ in ELEMENTS:
        ions = [x[q]/xtot for q in tracked + [high]]
        fr = [1.0 - sum(ions)] + ions
        prof = np.full((len(fr), len(mid)), np.nan)
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            m = (lt >= lo) & (lt < hi)
            if m.sum() >= 5:
                w = nh[m]
                prof[:, i] = [np.sum(f[m]*w)/w.sum() for f in fr]
        out[el] = prof
    return mid, out, d['Time']


snaps = [os.path.join(RUN, 'snr.out2.%05d.athdf' % i) for i in (6, 9)]
fig, axs = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
for row, fn in enumerate(snaps):
    mid, out, t = binned(fn)
    for ax, (el, _, tracked, _, labels) in zip(axs[row], ELEMENTS):
        tab = np.loadtxt(CHIANTI % el)
        top = len(tracked)
        for q, lab in enumerate(labels):
            cie = tab[:, 1 + q] if q <= top else tab[:, 1 + q:].sum(1)
            ax.plot(tab[:, 0], cie, color='C%d' % q, lw=1.2, alpha=0.5)
            ax.plot(mid, out[el][q], color='C%d' % q, marker='o', ms=3, lw=0.8, label=lab)
        ax.set(xlim=(4, 7.5), ylim=(1e-3, 1.5), yscale='log',
               title='%s, t = %.3f code (%s shell formation)'
               % (el, t, 'before' if t < 0.0385 else 'after'))
        ax.grid(alpha=0.25)
        if row == 0:
            ax.legend(fontsize=8, loc='lower left')
for ax in axs[-1]:
    ax.set_xlabel('log T [K]')
for ax in axs[:, 0]:
    ax.set_ylabel('fraction of the element')
fig.suptitle('Radiative SNR, GOW17 + O3,S3,N3: points, $n_{\\rm H}$-weighted means in T '
             'bins; lines, CHIANTI v11 CIE (stages above the top tracked one summed)',
             fontsize=12)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
