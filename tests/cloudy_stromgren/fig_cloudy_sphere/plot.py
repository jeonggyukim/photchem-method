"""Static H II region, Tigris against Cloudy 25. Radial T, n_e and the stage fractions
of H, He, O, S and N from two Tigris post-processed spheres (GOW17 + O3,S3,N3, 7 bands,
64^3): recombination photons carried by the diffuse ray tracer (filled), against Cloudy's
default diffuse field (solid), and absorbed on the spot (open), against Cloudy's
`diffuse OTS` (dashed). Tigris X_high (stages above the top tracked one) is compared with
Cloudy's X^3+. Writes ../../../figures/cloudy_sphere.png.

Inputs: ots_false and ots_true of ../fig_interior_temperature/run_recomb.sh; Cloudy
radial_profiles.txt of tests/cloudy_stromgren and its thermal_with_N_ots. The Tigris
profiles come from reduced_<run>.txt when they exist; `python plot.py --from-runs` reads
the runs ($PHOTCHEM_RUNS) and rewrites them."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from photchem_method import paths  # noqa: E402
from photchem_method import reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.join(HERE, '..', '..')
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'cloudy_sphere')
os.makedirs(os.path.dirname(OUT), exist_ok=True)
REDUCED = os.path.join(HERE, 'reduced_%s.txt')
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
    rb = np.arange(0, d['x1f'][-1] + 1e-9, 0.25*dx)
    ri = np.digitize(R.ravel(), rb)
    w = nh.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    prof = {k: np.bincount(ri, v.ravel()*w, len(rb) + 1)[1:len(rb)]/np.maximum(den, 1e-300)
            for k, v in f.items()}
    # shells a quarter cell wide; near the centre some hold no cell centre
    full = den > 0
    return 0.5*(rb[1:] + rb[:-1])[full], {k: v[full] for k, v in prof.items()}


def radius_half(r, xhp):
    i = np.nonzero(xhp < 0.5)[0][0]
    return np.interp(0.5, [xhp[i], xhp[i - 1]], [r[i], r[i - 1]])


# (run, legend label, matching Cloudy model, marker face)
RUNS = [('ots_false', r'Tigris ($N = 64^3$), no OTS', 'default', None),
        ('ots_true', r'Tigris ($N = 64^3$), OTS', 'diffuse OTS', 'white')]
TG = {}
for run, _, _, _ in RUNS:
    fn = REDUCED % run
    if '--from-runs' in sys.argv or not os.path.exists(fn):
        athena_read = paths.athena_read()
        src = paths.runs('T6_multi_ion', 'recomb_ots', run)
        r, tg = tigris_profile(src)
        reduced.save(fn, 'density-weighted radial profiles of the Tigris run '
                     'T6_multi_ion/recomb_ots/%s\nr [pc], T [K], ne [cm^-3], ions as '
                     'fractions of the element' % run, dict(r=r, **tg), sources=[src])
    else:
        tg = reduced.load(fn)
        r = tg.pop('r')
    TG[run] = (r, tg)
cl = [(lab, np.loadtxt(fn), ls) for lab, fn, ls in CLOUDY]
PANELS = [(r'$T$ [K]', ['T'], 'lin'), (r'$n_e$ [cm$^{-3}$]', ['ne'], 'lin'),
          ('H, He', ['H+', 'He+'], 'log'), ('O', ['O+', 'O++', 'O_high'], 'log'),
          ('S', ['S+', 'S++', 'S_high'], 'log'), ('N', ['N+', 'N++'], 'log')]
LABEL = {'H+': 'H$^+$', 'He+': 'He$^+$', 'O+': 'O$^+$', 'O++': 'O$^{2+}$',
         'O_high': r'O$^{\geq 3+}$', 'S+': 'S$^+$', 'S++': 'S$^{2+}$',
         'S_high': r'S$^{\geq 3+}$', 'N+': 'N$^+$', 'N++': 'N$^{2+}$'}
# one colour per charge, the same in every element panel (Okabe-Ito)
STAGE = {'+': '#0072B2', '++': '#D55E00', '_high': '#009E73'}
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'axes.linewidth': 0.7, 'xtick.direction': 'in',
                     'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'xtick.minor.visible': True, 'ytick.minor.visible': True,
                     'legend.frameon': False, 'legend.fontsize': 8})


def stage_colour(k):
    if k in ('T', 'ne'):
        return 'k'
    if k == 'He+':
        return STAGE['++']
    return STAGE['_high'] if k.endswith('_high') else STAGE['++' if k.endswith('++')
                                                             else '+']


fig, axs = plt.subplots(3, 2, figsize=(7.1, 7.6), sharex=True)
for ax, (ylab, keys, scale) in zip(axs.ravel(), PANELS):
    for k in keys:
        col = stage_colour(k)
        for lab, c, ls in cl:
            ax.plot(c[:, 0], c[:, CCOL[k]], color=col, ls=ls, lw=1.1)
        for i, (run, _, _, mfc) in enumerate(RUNS):
            r, tg = TG[run]
            ax.plot(r, tg[k], 'o', color=col, ms=3.0, mfc=mfc or col, mew=0.8,
                    label=LABEL.get(k) if i == 0 else None)
    if scale == 'log':
        ax.set_yscale('log')
        ax.set_ylim(1e-3, 8)
        ax.set_ylabel('fraction of ' + ylab)
        ax.legend(loc='upper left', ncol=3, borderaxespad=0.25, handletextpad=0.1, columnspacing=1.0)
    else:
        ax.set_ylabel(ylab)
    ax.set_xlim(0, 3.3)
for ax, t in zip(axs.ravel(), 'abcdef'):
    ax.text(0.97, 0.94, '(%s)' % t, transform=ax.transAxes, ha='right', va='top')
axs[0, 0].set_ylim(4000, 11000)
axs[0, 1].set_ylim(-5, 130)
L2D = matplotlib.lines.Line2D
handles = [L2D([], [], color='k', ls='-', lw=1.1),
           L2D([], [], color='k', ls='none', marker='o', ms=3.0, mew=0.8),
           L2D([], [], color='k', ls='--', lw=1.1),
           L2D([], [], color='k', ls='none', marker='o', ms=3.0, mfc='white', mew=0.8)]
axs[0, 0].legend(handles, ['Cloudy, no OTS', RUNS[0][1], 'Cloudy, OTS',
                           RUNS[1][1]], loc='upper left', handlelength=1.8)
for ax in axs[-1]:
    ax.set_xlabel(r'$r$ [pc]')
fig.tight_layout(h_pad=0.4, w_pad=1.0)
for ext in figure_formats():
    fig.savefig(OUT + '.' + ext, dpi=300)
cmodel = {lab: c for lab, c, _ in cl}
for lab, c in cmodel.items():
    print('R_s [pc] Cloudy %-11s %.3f' % (lab, radius_half(c[:, 0], c[:, 3])))
# ionization-weighted means inside 2.5 pc, each Tigris run against its Cloudy model
for run, rlab, clab, _ in RUNS:
    r, tg = TG[run]
    print('R_s [pc] %-26s %.3f' % (rlab, radius_half(r, tg['H+'])))
    c = cmodel[clab]
    m, mc = r < 2.5, c[:, 0] < 2.5
    wt, wc = r[m]**2, c[mc, 0]**2*np.gradient(c[mc, 0])
    for k in ['T', 'ne', 'He+', 'O+', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++']:
        vt = np.sum(tg[k][m]*wt)/np.sum(wt)
        vc = np.sum(c[mc, CCOL[k]]*wc)/np.sum(wc)
        print('  %-7s Tigris %.4g Cloudy %s %.4g ratio %.3f' % (k, vt, clab, vc, vt/vc))
print(OUT)
