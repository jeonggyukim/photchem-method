"""F13b: ion fractions against temperature in the radiative SNR with GOW17 + O3,S3,N3,
before and after shell formation, against CHIANTI v11 CIE at the same T. Per element, the
tracked stages and X_high (all stages above the top tracked one; CIE summed the same way).
Points: n_H-weighted mean over cells in 0.1 dex bins of T (bins with >= 5 cells).
Writes ../../../figures/F13_snr_ions.{pdf,png}.

Inputs: run_ions/snr.out2.00006 (t = 0.030) and .00009 (t = 0.045) of run_series.sh.
The binned fractions and the CIE curves come from reduced_ions.txt when it exists;
`python plot_ions.py --from-runs` reads the run ($PHOTCHEM_RUNS/M6_rad_snr/F13_series)
and the CHIANTI tables ($TIGRIS_DIR) and rewrites it."""
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
import reduced  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F13_snr_ions')
REDUCED = os.path.join(HERE, 'reduced_ions.txt')
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


if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    run = paths.runs('M6_rad_snr', 'F13_series', 'run_ions')
    data = {}
    for row, i in enumerate((6, 9)):
        data['mid'], out, data['t_%d' % row] = binned(
            os.path.join(run, 'snr.out2.%05d.athdf' % i))
        for el in out:
            data['prof_%s_%d' % (el, row)] = out[el]
    for el, _, tracked, _, labels in ELEMENTS:
        tab = np.loadtxt(paths.tigris('inputs/tables/chianti_v11', 'ioneq_%s.txt' % el))
        top = len(tracked)
        data['cie_logT_' + el] = tab[:, 0]
        data['cie_' + el] = [tab[:, 1 + q] if q <= top else tab[:, 1 + q:].sum(1)
                             for q in range(len(labels))]
    reduced.save(REDUCED, 'F13b: n_H-weighted stage fractions in log T bins (mid) of '
                 'M6_rad_snr/F13_series/run_ions\nsnr.out2.00006 (row 0) and .00009 '
                 '(row 1), one row per stage, and the CHIANTI v11 CIE\nfractions '
                 '(cie_logT, cie) of inputs/tables/chianti_v11/ioneq_<el>.txt; t [code]',
                 data)
else:
    data = reduced.load(REDUCED)

fig, axs = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
for row in range(2):
    mid, t = data['mid'], float(data['t_%d' % row])
    for ax, (el, _, tracked, _, labels) in zip(axs[row], ELEMENTS):
        for q, lab in enumerate(labels):
            ax.plot(data['cie_logT_' + el], data['cie_' + el][q], color='C%d' % q, lw=1.2,
                    alpha=0.5)
            ax.plot(mid, data['prof_%s_%d' % (el, row)][q], color='C%d' % q, marker='o', ms=3,
                    lw=0.8, label=lab)
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
