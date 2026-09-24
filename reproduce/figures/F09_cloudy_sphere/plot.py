"""F09: static H II region, Tigris against Cloudy 25. Radial T, n_e and the stage fractions
of H, He, O, S and N from the Tigris post-processed sphere (GOW17 + O3,S3,N3, 7 bands, 64^3)
against Cloudy 25 with its default diffuse field (solid) and `diffuse OTS` (dashed). Tigris
X_high (stages above the top tracked one) is compared with Cloudy's X^3+. Writes
../../../figures/F09_cloudy_sphere.{pdf,png}.

Inputs: run64_b7 of ../F10_band_count/run_series.sh; Cloudy radial_profiles.txt of
reproduce/cloudy_stromgren and its thermal_with_N_ots."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.join(HERE, '..', '..')
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F09_cloudy_sphere')
TIGRIS = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion/T6_multi_ion/'
                            'F10_bands/run64_b7')
CLOUDY = [('default', os.path.join(REPRO, 'cloudy_stromgren', 'radial_profiles.txt'), '-'),
          ('diffuse OTS', os.path.join(REPRO, 'cloudy_stromgren', 'thermal_with_N_ots',
                                       'radial_profiles.txt'), '--')]
CCOL = {'T': 1, 'ne': 2, 'H+': 3, 'He+': 4, 'He++': 5, 'O+': 7, 'O++': 8, 'O_high': 9,
        'S+': 11, 'S++': 12, 'S_high': 13, 'N+': 14, 'N++': 15}
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
# element totals per H of the run (Cloudy input uses the same)
XTOT = {'H': 1.0, 'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
MUH, KB = 2.34335276e-24, 1.380649e-16


def tigris_profile(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    nh = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    # rho is in m_H cm^-3 (1.4 n_H); the T formula takes rho with mu_H, which cancels
    f = {'T': d['press']*MUH*1e10/KB/(nh*(1.0 - x['H2'] + 0.1 + xe)), 'ne': xe*nh/1.4}
    for nm in NAMES:
        el = nm.split('+')[0].split('_')[0]
        if el in XTOT:
            f[nm] = x[nm]/XTOT[el]
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
    ri = np.digitize(R.ravel(), rb)
    w = nh.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    prof = {k: np.bincount(ri, v.ravel()*w, len(rb) + 1)[1:len(rb)]/np.maximum(den, 1e-300)
            for k, v in f.items()}
    return 0.5*(rb[1:] + rb[:-1]), prof


def radius_half(r, xhp):
    i = np.nonzero(xhp < 0.5)[0][0]
    return np.interp(0.5, [xhp[i], xhp[i - 1]], [r[i], r[i - 1]])


r, tg = tigris_profile(TIGRIS)
cl = [(lab, np.loadtxt(fn), ls) for lab, fn, ls in CLOUDY]
PANELS = [('T [K]', ['T'], 'lin'), (r'$n_e$ [cm$^{-3}$]', ['ne'], 'lin'),
          ('H, He', ['H+', 'He+'], 'log'), ('O', ['O+', 'O++', 'O_high'], 'log'),
          ('S', ['S+', 'S++', 'S_high'], 'log'), ('N', ['N+', 'N++'], 'log')]
LABEL = {'H+': 'H$^+$', 'He+': 'He$^+$', 'O+': 'O$^+$', 'O++': 'O$^{2+}$',
         'O_high': r'O$^{\geq 3+}$ (Cloudy O$^{3+}$)', 'S+': 'S$^+$', 'S++': 'S$^{2+}$',
         'S_high': r'S$^{\geq 3+}$ (Cloudy S$^{3+}$)', 'N+': 'N$^+$', 'N++': 'N$^{2+}$'}
fig, axs = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
for ax, (ylab, keys, scale) in zip(axs.ravel(), PANELS):
    for q, k in enumerate(keys):
        col = 'k' if len(keys) == 1 else 'C%d' % q
        for lab, c, ls in cl:
            ax.plot(c[:, 0], c[:, CCOL[k]], color=col, ls=ls, lw=1.8, alpha=0.6)
        ax.plot(r, tg[k], 'o', color=col, ms=4, label=LABEL.get(k))
    if scale == 'log':
        ax.set_yscale('log')
        ax.set_ylim(1e-3, 1.5)
        ax.set_ylabel('fraction of ' + ylab)
        ax.legend(fontsize=8, loc='lower left')
    else:
        ax.set_ylabel(ylab)
    ax.set_xlim(0, 3.5)
    ax.grid(alpha=0.25)
axs[0, 0].set_ylim(4000, 11000)
axs[0, 0].plot([], [], 'k-', label='Cloudy 25, default')
axs[0, 0].plot([], [], 'k--', label='Cloudy 25, diffuse OTS')
axs[0, 0].plot([], [], 'ko', ms=4, label='Tigris GOW17 + O3,S3,N3, 7 bands')
axs[0, 0].legend(fontsize=8, loc='upper left')
for ax in axs[-1]:
    ax.set_xlabel('r [pc]')
rs = [radius_half(r, tg['H+'])] + [radius_half(c[:, 0], c[:, 3]) for _, c, _ in cl]
fig.suptitle(r'Static H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, '
             r'SB99 2 Myr. $R_s$: Tigris %.3f pc, Cloudy %.3f pc (default), %.3f pc (OTS)'
             % tuple(rs), fontsize=11)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
# ionization-weighted means inside 2.5 pc, Tigris vs Cloudy default
m = r < 2.5
for k in ['T', 'ne', 'He+', 'O+', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++']:
    c = cl[0][1]
    mc = c[:, 0] < 2.5
    wt, wc = r[m]**2, c[mc, 0]**2*np.gradient(c[mc, 0])
    vt = np.sum(tg[k][m]*wt)/np.sum(wt)
    vc = np.sum(c[mc, CCOL[k]]*wc)/np.sum(wc)
    print('%-7s Tigris %.4g Cloudy %.4g ratio %.3f' % (k, vt, vc, vt/vc))
print(OUT)
