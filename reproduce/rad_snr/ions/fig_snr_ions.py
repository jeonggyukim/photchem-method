"""Radiative SNR with GOW17 + O2,S3,N2: ion fractions against temperature, before
and after shell formation, with CHIANTI CIE at the same T; and the hot-gas mass
and radial momentum histories of the two SNR variants with t_sf and the NCR
reference values marked.

usage: python fig_snr_ions.py OUT.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = '/Users/jgkim/Documents/tigris-photchem-gow17-multi-ion/M6_rad_snr/'
MH, KB, MU_H, X_HE = 1.6735575e-24, 1.380649e-16, 1.4, 0.1
TOT = {'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}
STAGES = {'O': ['O+', 'O++'], 'S': ['S+', 'S++', 'S3+'], 'N': ['N+', 'N++']}
LABEL = {'O': ['O I', 'O II', 'O III'], 'S': ['S I', 'S II', 'S III', 'S IV'],
         'N': ['N I', 'N II', 'N III']}
CHIANTI = '/Users/jgkim/Projects/tigris-gow17/inputs/tables/chianti_v11/ioneq_%s.txt'
REF = {'M_hot': 1483.49, 'pr': 192798., 't_sf': 0.038473}


def binned(fn):
    d = athena_read.athdf(fn)
    names = CORE + IONS
    keys = ['rHI'] + ['r%d' % n for n in range(4, 3 + len(names))]
    x = {nm: d[k].astype(np.float64).ravel() for nm, k in zip(names, keys)}
    nH = d['rho'].astype(np.float64).ravel()/MU_H
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
    T = d['press'].ravel()*MH*1e10/KB/(nH*(1.0 - x['H2'] + X_HE + xe))
    lt = np.log10(np.maximum(T, 1.0))
    edges = np.arange(4.0, 7.01, 0.1)
    mid = 0.5*(edges[1:] + edges[:-1])
    out = {}
    for el, st in STAGES.items():
        fr = [1.0 - sum(x[q] for q in st)/TOT[el]] + [x[q]/TOT[el] for q in st]
        prof = np.full((len(fr), len(mid)), np.nan)
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            m = (lt >= lo) & (lt < hi)
            if m.sum() >= 5:
                w = nH[m]
                prof[:, i] = [np.sum(f[m]*w)/w.sum() for f in fr]
        out[el] = prof
    return mid, out, d['Time']


fig = plt.figure(figsize=(17, 10))
gs = fig.add_gridspec(3, 3)
snaps = [(HERE + 'ions/run_hdf5/snr_ions.out2.00006.athdf', '--'),
         (HERE + 'ions/run_hdf5/snr_ions.out2.00010.athdf', '-')]
for j, el in enumerate(STAGES):
    tab = np.loadtxt(CHIANTI % el)
    top = len(STAGES[el])
    for row in range(2):
        ax = fig.add_subplot(gs[row, j])
        fn, ls = snaps[row]
        mid, out, t = binned(fn)
        for q in range(top + 1):
            cie = tab[:, 1 + q] if q < top else tab[:, 1 + top:].sum(1)
            ax.plot(tab[:, 0], cie, color='C%d' % q, lw=1, alpha=0.5)
            ax.plot(mid, out[el][q], color='C%d' % q, ls='-', marker='o', ms=3, alpha=0.8,
                    label=LABEL[el][q] + ('+' if q == top else ''))
        ax.set(xlim=(4, 7), ylim=(1e-3, 1.5), yscale='log', xlabel='log T [K]',
               title='%s, t = %.3f (%s shell formation)'
               % (el, t, 'before' if t < REF['t_sf'] else 'after'))
        ax.legend(fontsize=8, loc='lower left', title='points: GOW17; lines: CHIANTI CIE',
                  title_fontsize=7)

runs = [(HERE + 'ions/run/snr_ions.hst', 'SNR, Subcell'),
        (HERE + 'ions_particle/run/snr_ions_SB99.hst', 'SNR particle, Classic')]
for k, (q, lab) in enumerate([('M_hot', r'$M_{\rm hot}$ [M$_\odot$]'),
                              ('pr', r'$p_r$ [M$_\odot$ km s$^{-1}$]')]):
    ax = fig.add_subplot(gs[2, k])
    for i, (fn, name) in enumerate(runs):
        h = athena_read.hst(fn)
        par = athena_read.athinput(fn.rsplit('/', 1)[0] + '/athinput.runtime')
        to_msun = float(par['units']['mass_cgs'])/1.98841e33
        if q in h:
            ax.plot(h['time'], h[q]*to_msun, color='C%d' % i, alpha=0.7, label=name)
    ax.axvline(REF['t_sf'], color='k', ls=':', lw=1, label=r'NCR $t_{\rm sf}$')
    ax.axhline(REF[q], color='k', ls='--', lw=1, label='NCR value at $t_{\\rm sf}$')
    ax.set(xlabel='t [code time = 0.978 Myr]', title=lab)
    ax.legend(fontsize=8)
ax = fig.add_subplot(gs[2, 2])
ax.axis('off')
ax.text(0, 1, 'Shell formation vs NCR reference (tol 5%)\n\n'
        '              Subcell    particle\n'
        't_sf          +0.0%      +0.0%\n'
        'M_sf          -2.2%      -1.8%\n'
        'p_sf          -0.0%      -0.0%\n'
        'f_ej          +2.2%      +1.8%\n\n'
        'Top tracked stage (O III, S IV, N III) holds\n'
        'every higher stage: the CIE pool is not wired.\n'
        'Above 3.5e4 K the metal cooling is the CIE table.',
        va='top', family='monospace', fontsize=10)
fig.suptitle('Radiative SNR, GOW17 + O2,S3,N2: mass-weighted ion fractions in 0.1-dex '
             'T bins vs CHIANTI CIE at the same T', fontsize=13)
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=100)
print(sys.argv[1])
