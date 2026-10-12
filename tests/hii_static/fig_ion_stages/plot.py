"""Static H II region, Tigris against Cloudy 25, every ionization stage of H, He, O, S
and N. One spectrum per figure: sb99 (Starburst99 at 2 Myr) or bb5e4 (a 5e4 K
blackbody); Q(H) = 1e49 photons s^-1, uniform dust-free gas of n_H = 100 cm^-3. Tigris
is GOW17 + O3,S3,N3 with He++, 10 bands, 64^3, 20 point rays per cell, recombination
photons on the spot (OTS) or carried by the diffuse ray tracer (no OTS); Cloudy is
`diffuse OTS` and its default transfer of the diffuse field. Columns: T and n_e, then
H, He, O, S, N; rows: neutral, singly, doubly and higher ionized. Tigris X_high (every
stage above X^2+) is compared with Cloudy's stages 3+ and above summed. Writes
../../../figures/ion_stages_<spectrum>.png.

Inputs: run_tigris.sh output in $PHOTCHEM_RUNS/T6_multi_ion/hii_static_bands10
(<spectrum>_ots, <spectrum>_no_ots) and ../blackbody/run_cloudy_bb.sh output in
$PHOTCHEM_RUNS/T6_multi_ion/cloudy_<spectrum> (ots, no_ots). The profiles come from
reduced_<spectrum>.txt when it exists; `python plot.py SPECTRUM --from-runs` reads the
runs and rewrites it.
Usage: python plot.py [SPECTRUM] [--from-runs]   (default: sb99)"""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from photchem_method import paths  # noqa: E402
from photchem_method import reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

ARGS = [a for a in sys.argv[1:] if a != '--from-runs']
SPEC = ARGS[0] if ARGS else 'sb99'
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'ion_stages_' + SPEC)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
REDUCED = os.path.join(HERE, 'reduced_%s.txt' % SPEC)
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high',
         'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
# element totals per H of the runs (the Cloudy inputs use the same)
XTOT = {'H': 1.0, 'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
MUH, KB, PC = 2.34335276e-24, 1.380649e-16, 3.0856776e18
# columns top to bottom: neutral, singly, doubly, higher ions; None leaves a panel empty
LAYOUT = [['T', 'ne', 'heat', None], ['H0', 'H+', None, None],
          ['He0', 'He+', 'He++', None]] + \
         [[el + '0', el + '+', el + '++', el + '_high'] for el in 'OSN']
# heat: total heating; cool: total cooling, drawn in the heat panel
QTY = [q for col in LAYOUT for q in col if q] + ['cool']
LAB = {'T': r'$T$ [K]', 'ne': r'$n_e$ [cm$^{-3}$]',
       'heat': r'$\Gamma$, $n_{\rm H}\Lambda$ [$10^{-22}$ erg s$^{-1}$]',
       'H0': r'H$^0$/H', 'H+': r'H$^+$/H',
       'He0': r'He$^0$/He', 'He+': r'He$^+$/He', 'He++': r'He$^{2+}$/He'}
for el in 'OSN':
    LAB.update({el + '0': r'%s$^0$/%s' % (el, el), el + '+': r'%s$^+$/%s' % (el, el),
                el + '++': r'%s$^{2+}$/%s' % (el, el),
                el + '_high': r'%s$^{\geq3+}$/%s' % (el, el)})
MODES = [('no_ots', 'no OTS'), ('ots', 'OTS')]
# panels drawn for both modes; inside the front the no-OTS abundances and n_e lie within
# 16% of the OTS ones (H0 45%, He++ 27%, N0 26%), so the other panels show OTS only
BOTH_MODES = ('T', 'heat')


def tigris(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    nh = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    # rho is in m_H cm^-3 (1.4 n_H); the T formula takes rho with mu_H, which cancels
    f = {'T': d['press']*MUH*1e10/KB/(nh*(1.0 - x['H2'] + 0.1 + xe)), 'ne': xe*nh/1.4,
         'H0': np.maximum(1.0 - x['H+'] - 2*x['H2'], 0.0), 'H+': x['H+'],
         'heat': d['heat_rate'], 'cool': d['cool_rate']}
    for el, ions in [('He', ['He+', 'He++']), ('O', ['O+', 'O++', 'O_high']),
                     ('S', ['S+', 'S++', 'S_high']), ('N', ['N+', 'N++', 'N_high'])]:
        for nm in ions:
            f[nm] = x[nm]/XTOT[el]
        f[el + '0'] = np.maximum(1.0 - sum(f[nm] for nm in ions), 0.0)
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    # density-weighted shells a quarter cell wide; near the centre some hold no cell
    rb = np.arange(0, d['x1f'][-1] + 1e-9, 0.25*dx)
    ri = np.digitize(R.ravel(), rb)
    w = nh.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    full = den > 0
    prof = {k: (np.bincount(ri, v.ravel()*w, len(rb) + 1)[1:len(rb)]
                / np.maximum(den, 1e-300))[full] for k, v in f.items()}
    return 0.5*(rb[1:] + rb[:-1])[full], prof


def cloudy(mdir):
    c = np.loadtxt(os.path.join(mdir, 'radial_profiles.txt'))
    rad = np.loadtxt(os.path.join(mdir, 'hii.rad'))
    r_rad = rad[:, 1]/PC
    p = {'T': c[:, 1], 'ne': c[:, 2]}
    # save cooling: depth, T, Htot, Ctot [erg cm^-3 s^-1]; per H with n_H = 100 cm^-3
    hc = np.loadtxt(os.path.join(mdir, 'hii.cool'), usecols=(2, 3))
    p['heat'] = np.interp(c[:, 0], r_rad, hc[:, 0])/100.0
    p['cool'] = np.interp(c[:, 0], r_rad, hc[:, 1])/100.0
    # save element columns: depth, X0, X+, X2+, ...; Cloudy appends extra columns after
    # the bare nucleus, so the stages stop at column 1 + Z
    for el, z in [('H', 1), ('He', 2), ('O', 8), ('S', 16), ('N', 7)]:
        a = np.loadtxt(os.path.join(mdir, 'hii.ele_' + el))
        st = lambda i0, i1: np.interp(c[:, 0], r_rad, a[:, i0:i1].sum(1))  # noqa: E731
        p[el + '0'], p[el + '+'] = st(1, 2), st(2, 3)
        if z >= 2:
            p[el + '++'] = st(3, 4)
        if z >= 3:
            p[el + '_high'] = st(4, 2 + z)
    return c[:, 0], p


def radius_half(r, x):
    i = np.nonzero(x < 0.5)[0][0]
    return np.interp(0.5, [x[i], x[i - 1]], [r[i], r[i - 1]])


if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    arrays, sources = {}, []
    for mode, _ in MODES:
        tdir = paths.runs('T6_multi_ion', 'hii_static_bands10', '%s_%s' % (SPEC, mode))
        cdir = paths.runs('T6_multi_ion', 'cloudy_' + SPEC, mode)
        for code, (r, p) in [('tigris', tigris(tdir)), ('cloudy', cloudy(cdir))]:
            arrays['%s_%s_r' % (code, mode)] = r
            arrays.update({'%s_%s_%s' % (code, mode, q): p[q] for q in QTY})
        sources += [tdir, cdir]
    reduced.save(REDUCED, 'radial profiles of the static H II region, spectrum %s: Tigris '
                 'density-weighted shells, Cloudy zones\n<code>_<mode>_<quantity>; r [pc], '
                 'T [K], ne [cm^-3], ions as fractions of the element, heat (Gamma) and '
                 'cool (n_H Lambda) per H [erg s^-1]' % SPEC, arrays,
                 sources=sources)
data = reduced.load(REDUCED)

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'xtick.major.pad': 2, 'ytick.major.pad': 1.5,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
                     'legend.fontsize': 9})
COL = {'no_ots': '#D55E00', 'ots': '#0072B2'}
MARK = {'no_ots': 'o', 'ots': 's'}
NCOL, NROW = len(LAYOUT), 4
# a narrow empty column after the first leaves room for the abundance tick labels
fig, allax = plt.subplots(NROW, NCOL + 1, figsize=(9.5, 6.8), sharex=True,
                          gridspec_kw={'width_ratios': [1, 0.22] + [1]*(NCOL - 1)})
for a in allax[:, 1]:
    a.set_visible(False)
axs = allax[:, [0] + list(range(2, NCOL + 1))]
ax = {}
for col, acol in zip(LAYOUT, axs.T):
    for q, a in zip(col, acol):
        if q is None:
            a.set_visible(False)
        else:
            ax[q] = a
rmax = 0.0
for i, (mode, label) in enumerate(MODES):
    rc, rt = data['cloudy_%s_r' % mode], data['tigris_%s_r' % mode]
    pc = {q: data['cloudy_%s_%s' % (mode, q)] for q in QTY}
    pt = {q: data['tigris_%s_%s' % (mode, q)] for q in QTY}
    rmax = max(rmax, radius_half(rc, pc['H+']))
    print('%s: R_s [pc] Tigris %.3f Cloudy %.3f' % (label, radius_half(rt, pt['H+']),
                                                     radius_half(rc, pc['H+'])))
    for q in ('heat', 'cool'):
        pc[q], pt[q] = 1e22*pc[q], 1e22*pt[q]
    for q in QTY:
        # OTS and no OTS differ in T and the heating; elsewhere OTS alone is shown
        if q == 'cool' or (mode == 'no_ots' and q not in BOTH_MODES):
            continue
        ax[q].plot(rc, pc[q], color=COL[mode], lw=1.8, alpha=0.5,
                   ls='--' if mode == 'ots' else '-', zorder=1)
        ax[q].plot(rt, pt[q], MARK[mode], color=COL[mode], ms=2.6,
                   mfc='white' if mode == 'ots' else COL[mode], mew=0.7, ls='none',
                   markevery=(i, 2), alpha=0.75, zorder=3)
    # Tigris cooling on its heating: a cell in balance puts the + on the marker
    ax['heat'].plot(rt, pt['cool'], '+', color=COL[mode], ms=4.5, mew=0.8, ls='none',
                    markevery=(i, 2), zorder=4)
for q in QTY:
    if q == 'cool':
        continue
    a = ax[q]
    if q in ('T', 'ne', 'heat'):
        a.set_ylabel(LAB[q])
        a.yaxis.labelpad = 1
    else:
        a.set_title(LAB[q], fontsize=9.5, pad=2)
        a.set_yscale('log')
        a.set_ylim(1e-5, 1e1)
        a.set_yticks([1e-4, 1e-2, 1])
    if q == 'T':
        a.set_ylim(5000, 12000)
    if q == 'ne':
        a.set_ylim(0, 1.25*a.get_ylim()[1])
    if q == 'heat':
        a.set_ylim(0, 1.3*a.get_ylim()[1])
        key = [Line2D([], [], color='0.3', marker='o', ls='none', ms=3),
               Line2D([], [], color='0.3', marker='+', ls='none', ms=5, mew=0.9)]
        a.legend(key, [r'$\Gamma$', r'$n_{\rm H}\Lambda$'], loc='upper right', ncol=2,
                 fontsize=8, handletextpad=0.2, columnspacing=0.8, borderaxespad=0.3)
for a in axs.ravel():
    a.set_xlim(0, 1.15*rmax)
# the abundance columns share one y range: tick labels only on the first visible one
for arow in axs:
    first = True
    for a in arow[1:]:
        if a.get_visible():
            if not first:
                a.tick_params(labelleft=False)
            first = False
# x labels on the lowest visible panel of each column
for acol in axs.T:
    low = [a for a in acol if a.get_visible()][-1]
    low.tick_params(labelbottom=True)
    low.set_xlabel(r'$r$ [pc]')
h = [Line2D([], [], color=COL['no_ots'], lw=1.8, alpha=0.5),
     Line2D([], [], color=COL['no_ots'], marker='o', ls='none', ms=3),
     Line2D([], [], color=COL['ots'], lw=1.8, alpha=0.5, ls='--'),
     Line2D([], [], color=COL['ots'], marker='s', ls='none', ms=3, mfc='white')]
fig.legend(h, ['Cloudy, no OTS', r'Tigris ($N=64^3$), no OTS', 'Cloudy, OTS',
               r'Tigris ($N=64^3$), OTS'], loc='upper center', ncol=4,
           bbox_to_anchor=(0.5, 1.0))
fig.subplots_adjust(left=0.065, right=0.995, bottom=0.07, top=0.92, wspace=0.08,
                    hspace=0.28)
for ext in figure_formats():
    fig.savefig(OUT + '.' + ext, dpi=250)
print(OUT)
