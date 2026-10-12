"""Static H II region (Q(H) = 1e49 photons s^-1, uniform dust-free gas of n_H = 100
cm^-3), Tigris (GOW17 + O3,S3,N3 with He++, 10 bands, 64^3): slices through the plane
of cell centres nearest the source (the source sits on a cell corner at the origin,
the plane is at z = 0.094 pc), laid out as ../fig_ion_stages. First column: T for
recombination photons carried by the diffuse ray tracer (no OTS) and absorbed on the
spot (OTS), n_e and the heating per H, Gamma; then H, He, O, S, N, one row per stage
(neutral, singly, doubly and higher ionized). Everything but the first T is OTS; the
no-OTS ion fractions differ by at most 45% inside the front (see ../fig_ion_stages).
Writes ../../../figures/slices_<spectrum>.png.

Reads ../fig_ion_stages/run_tigris.sh output in
$PHOTCHEM_RUNS/T6_multi_ion/hii_static_bands10/<spectrum>_{ots,no_ots} through
reduced_<spectrum>.txt (|x|, |y| < 3.75 pc); `python plot.py SPECTRUM --from-runs`
rewrites it.
Usage: python plot.py [SPECTRUM] [--from-runs]   (default: bb5e4)"""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from photchem_method import paths  # noqa: E402
from photchem_method import reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

ARGS = [a for a in sys.argv[1:] if a != '--from-runs']
SPEC = ARGS[0] if ARGS else 'bb5e4'
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'slices_' + SPEC)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
REDUCED = os.path.join(HERE, 'reduced_%s.txt' % SPEC)
HALF = 3.75   # [pc], half width of the region shown
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++',
         'N_high', 'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
XTOT = {'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
IONS = {'He': ['He+', 'He++'], 'O': ['O+', 'O++', 'O_high'],
        'S': ['S+', 'S++', 'S_high'], 'N': ['N+', 'N++', 'N_high']}
# columns top to bottom; None leaves a panel empty
LAYOUT = [['T_no_ots', 'T_ots', 'ne', 'heat'], ['H0', 'H+', None, None],
          ['He0', 'He+', 'He++', None]] + \
         [[el + '0', el + '+', el + '++', el + '_high'] for el in 'OSN']
LAB = {'T_no_ots': r'$T$, no OTS', 'T_ots': r'$T$, OTS', 'ne': r'$n_e$ [cm$^{-3}$]',
       'heat': r'$\Gamma$ [$10^{-22}$ erg s$^{-1}$]', 'H0': r'H$^0$/H',
       'H+': r'H$^+$/H', 'He0': r'He$^0$/He', 'He+': r'He$^+$/He',
       'He++': r'He$^{2+}$/He'}
for el in 'OSN':
    LAB.update({el + '0': r'%s$^0$/%s' % (el, el), el + '+': r'%s$^+$/%s' % (el, el),
                el + '++': r'%s$^{2+}$/%s' % (el, el),
                el + '_high': r'%s$^{\geq3+}$/%s' % (el, el)})
QTY = [q for col in LAYOUT for q in col if q]


def tigris_slice(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    k = int(np.where(d['x3v'] > 0, d['x3v'], np.inf).argmin())  # first plane above 0
    i = np.nonzero(np.abs(d['x1v']) < HALF)[0]
    j = np.nonzero(np.abs(d['x2v']) < HALF)[0]
    cut = (k, slice(j[0], j[-1] + 1), slice(i[0], i[-1] + 1))
    x = {nm: d[key][cut].astype(np.float64)
         for nm, key in zip(NAMES, ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))])}
    s = {'T': d['temp'][cut].astype(np.float64),
         'ne': 100.0*sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES),
         'heat': 1e22*d['heat_rate'][cut].astype(np.float64),
         'H0': np.maximum(1.0 - x['H+'] - 2*x['H2'], 0.0), 'H+': x['H+']}
    for el, ions in IONS.items():
        for nm in ions:
            s[nm] = x[nm]/XTOT[el]
        s[el + '0'] = np.maximum(1.0 - sum(s[nm] for nm in ions), 0.0)
    x1f = d['x1f'][i[0]:i[-1] + 2]
    x2f = d['x2f'][j[0]:j[-1] + 2]
    return x1f, x2f, float(d['x3v'][k]), s


if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    arrays, sources = {}, []
    for mode in ('ots', 'no_ots'):
        tdir = paths.runs('T6_multi_ion', 'hii_static_bands10', '%s_%s' % (SPEC, mode))
        x1f, x2f, z, s = tigris_slice(tdir)
        arrays['T_' + mode] = s['T']
        if mode == 'ots':
            arrays.update({q: s[q] for q in QTY if not q.startswith('T_')})
            arrays['x1f'], arrays['x2f'], arrays['z'] = x1f, x2f, np.array([z])
        sources.append(tdir)
    reduced.save(REDUCED, 'slices of the static H II region, spectrum %s, at z [pc]\n'
                 'x1f, x2f: cell faces [pc]; T_<mode> [K]; ne [cm^-3]; heat per H '
                 '[1e-22 erg s^-1]; ions as fractions of the element (OTS); arrays '
                 '[x2, x1]' % SPEC, arrays, fmts={k: '%.6g' for k in arrays},
                 sources=sources)
data = reduced.load(REDUCED)

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
NCOL, NROW = len(LAYOUT), 4
fig, allax = plt.subplots(NROW, NCOL + 1, figsize=(10.0, 6.9), sharex=True,
                          sharey=True,
                          gridspec_kw={'width_ratios': [1, 0.30] + [1]*(NCOL - 1)})
for a in allax[:, 1]:
    a.set_visible(False)
axs = allax[:, [0] + list(range(2, NCOL + 1))]
NORM = {'T': Normalize(5e3, 1.1e4), 'ne': Normalize(0, 125), 'heat': Normalize(0, 6),
        'ion': LogNorm(1e-5, 1)}
CMAP = {'T': 'magma', 'ne': 'cividis', 'heat': 'inferno',
        'ion': plt.get_cmap('viridis').copy()}
# a fraction of exactly 0 (a neutral that the ladder emptied) takes the lowest colour
CMAP['ion'].set_bad(CMAP['ion'](0.0))
ims = {}
for col, acol in zip(LAYOUT, axs.T):
    for q, a in zip(col, acol):
        if q is None:
            a.set_visible(False)
            continue
        kind = 'T' if q.startswith('T_') else q if q in ('ne', 'heat') else 'ion'
        ims[kind] = a.pcolormesh(data['x1f'], data['x2f'], data[q], norm=NORM[kind],
                                 cmap=CMAP[kind], shading='flat', rasterized=True)
        a.set_aspect('equal')
        a.set_title(LAB[q], fontsize=9, pad=2)
        a.set_xticks([-3, 0, 3])
        a.set_yticks([-3, 0, 3])
for acol in axs.T:
    vis = [a for a in acol if a.get_visible()]
    vis[-1].tick_params(labelbottom=True)
    vis[-1].set_xlabel(r'$x$ [pc]', labelpad=1)
for arow in axs:
    arow[0].set_ylabel(r'$y$ [pc]', labelpad=1)
    for a in arow[1:]:
        if a.get_visible():
            a.tick_params(labelleft=True)
            break
fig.subplots_adjust(left=0.05, right=0.92, bottom=0.07, top=0.93, wspace=0.06,
                    hspace=0.22)
# colour bars: the ion fractions right of the grid; T (one bar for both rows), n_e and
# Gamma beside the first column, in the spacer
cax_ion = fig.add_axes([0.935, 0.07, 0.012, 0.86])
fig.colorbar(ims['ion'], cax=cax_ion, label='fraction of the element')
for rows, kind in [((0, 1), 'T'), ((2,), 'ne'), ((3,), 'heat')]:
    top, bot = axs[rows[0], 0].get_position(), axs[rows[-1], 0].get_position()
    cax = fig.add_axes([top.x1 + 0.006, bot.y0, 0.008, top.y1 - bot.y0])
    fig.colorbar(ims[kind], cax=cax).ax.tick_params(labelsize=6)
fig.suptitle(r'%s, slice at $z = %.3f$ pc, Tigris ($N = 64^3$)'
             % ('5e4 K blackbody' if SPEC == 'bb5e4' else 'SB99 2 Myr',
                float(data['z'][0])), fontsize=10)
for ext in figure_formats():
    fig.savefig(OUT + '.' + ext, dpi=250)
print(OUT)
