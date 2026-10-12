"""z = 0 slice movie of the expanding H II region of run_tigris.sh (GOW17 + O3,S3,N3
and He++, 8 bands, no dust, OTS, 128^3), laid out as ../../hii_static/fig_slices:
first column n_H, T, P/k and |v|, then H, He, O, S, N, one row per stage (neutral,
singly, doubly and higher ionized). Every panel keeps one colour range over all
frames.

usage: python slice_movie.py RUNDIR OUTDIR [--nproc N] [--frame I]
Writes OUTDIR/frame_NNNNN.png, one per HII.out2 dump, then OUTDIR/slices.mp4 (ffmpeg).
The frames are independent and are drawn by N processes (default 8). --frame I draws
only dump I, without the movie, to check the layout first.
"""
import argparse
import glob
import multiprocessing
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm, Normalize  # noqa: E402
from photchem_method import paths  # noqa: E402

athena_read = paths.athena_read()
# ism units: density unit m_H/cm^3, rho = 1.4 n_H; pressure unit m_H (km/s)^2;
# time unit 0.978 Myr
MH, KB, MU_H, X_HE = 1.6735575e-24, 1.380649e-16, 1.4, 0.1
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++',
         'N_high', 'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
XTOT = {'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
IONS = {'He': ['He+', 'He++'], 'O': ['O+', 'O++', 'O_high'],
        'S': ['S+', 'S++', 'S_high'], 'N': ['N+', 'N++', 'N_high']}


def layout(molecular):
    """Columns top to bottom; --molecular adds H2 and the carbon carriers."""
    return [['nH', 'T', 'Pk', 'v']] + \
        ([['H2', 'CO', 'C+', 'C0']] if molecular else []) + \
        [['H0', 'H+', None, None], ['He0', 'He+', 'He++', None]] + \
        [[el + '0', el + '+', el + '++', el + '_high'] for el in 'OSN']


LAB = {'nH': r'$n_{\rm H}$ [cm$^{-3}$]', 'T': r'$T$ [K]', 'Pk': r'$P/k$ [K cm$^{-3}$]',
       'v': r'$|v|$ [km s$^{-1}$]', 'H0': r'H$^0$/H', 'H+': r'H$^+$/H',
       'He0': r'He$^0$/He', 'He+': r'He$^+$/He', 'He++': r'He$^{2+}$/He',
       'H2': r'2H$_2$/H', 'CO': 'CO/C', 'C+': r'C$^+$/C', 'C0': r'C$^0$/C'}
X_C = 1.6e-4
for el in 'OSN':
    LAB.update({el + '0': r'%s$^0$/%s' % (el, el), el + '+': r'%s$^+$/%s' % (el, el),
                el + '++': r'%s$^{2+}$/%s' % (el, el),
                el + '_high': r'%s$^{\geq3+}$/%s' % (el, el)})
NORM = {'nH': LogNorm(10, 1e3), 'T': Normalize(5e3, 1.1e4), 'Pk': LogNorm(1e3, 3e6),
        'v': Normalize(0, 15), 'ion': LogNorm(1e-5, 1)}
CMAP = {'nH': 'cividis', 'T': 'inferno', 'Pk': 'magma', 'v': 'plasma',
        'ion': plt.get_cmap('viridis').copy()}
CMAP['ion'].set_bad(CMAP['ion'](0.0))   # a fraction of exactly 0 takes the lowest colour
# --temp-log: T on a log scale, 10-2e4 K, in pyathena's temperature colour map (RdYlBu_r
# with its middle moved to 1e4 K, as pyathena/fields/fields.py cmap_shift)
T_LOG = (10.0, 2.0e4)


def shifted_cmap(base, mid):
    """`base` resampled so that its middle sits at `mid` (0-1) instead of 0.5."""
    s = np.linspace(0.0, 1.0, 256)
    src = np.where(s < mid, 0.5*s/mid, 0.5 + 0.5*(s - mid)/(1.0 - mid))
    return ListedColormap(plt.get_cmap(base)(src), name='T_pyathena')


CMAP_T_LOG = shifted_cmap('RdYlBu_r', np.log10(1e4/T_LOG[0])/np.log10(T_LOG[1]/T_LOG[0]))


def fields(fn):
    d = athena_read.athdf(fn)
    k = int(np.where(d['x3v'] > 0, d['x3v'], np.inf).argmin())  # first plane above 0
    x = {nm: d[key][k].astype(np.float64)
         for nm, key in zip(NAMES, ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))])}
    f = {'nH': d['rho'][k]/MU_H}
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    f['Pk'] = d['press'][k]*MH*1e10/KB
    f['T'] = f['Pk']/(f['nH']*(1.0 - x['H2'] + X_HE + xe))
    f['v'] = np.sqrt(d['vel1'][k]**2 + d['vel2'][k]**2 + d['vel3'][k]**2)
    f['H0'] = np.maximum(1.0 - x['H+'] - 2*x['H2'], 0.0)
    f['H+'] = x['H+']
    f['H2'] = 2*x['H2']
    f['CO'], f['C+'] = x['CO']/X_C, x['C+']/X_C
    f['C0'] = np.maximum(1.0 - (x['CO'] + x['C+'] + x['CHx'] + x['HCO+']
                                + x['C_high'])/X_C, 0.0)
    for el, ions in IONS.items():
        for nm in ions:
            f[nm] = x[nm]/XTOT[el]
        f[el + '0'] = np.maximum(1.0 - sum(f[nm] for nm in ions), 0.0)
    return f, d['x1f'], d['x2f'], float(d['Time'])


plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})


def frame(job):
    """Draw dump fn as OUTDIR/frame_<n>.png."""
    n, fn, outdir, molecular, title, scale, temp_log = job
    # the n_H and P/k ranges follow the ambient density (--density-scale)
    norm = dict(NORM, nH=LogNorm(10*scale, 1e3*scale), Pk=LogNorm(1e3*scale, 3e6*scale))
    cmap = dict(CMAP)
    if temp_log:
        norm['T'], cmap['T'] = LogNorm(*T_LOG), CMAP_T_LOG
    LAYOUT = layout(molecular)
    NCOL = len(LAYOUT)
    f, x1f, x2f, t = fields(fn)
    # panels square: the figure widens by one panel per column
    fig, allax = plt.subplots(4, NCOL + 1, figsize=(10.0*(NCOL + 0.45)/6.45, 6.1),
                              sharex=True,
                              sharey=True,
                              gridspec_kw={'width_ratios': [1, 0.45] + [1]*(NCOL - 1),
                                           'wspace': 0.0, 'hspace': 0.0})
    for a in allax[:, 1]:
        a.set_visible(False)
    axs = allax[:, [0] + list(range(2, NCOL + 1))]
    ims = {}
    for col, acol in zip(LAYOUT, axs.T):
        for q, a in zip(col, acol):
            if q is None:
                a.set_visible(False)
                continue
            kind = q if q in NORM else 'ion'
            ims[kind] = a.pcolormesh(x1f, x2f, f[q], norm=norm[kind], cmap=cmap[kind],
                                     shading='flat', rasterized=True)
            a.set_aspect('equal')
            a.text(0.05, 0.94, LAB[q], transform=a.transAxes, ha='left', va='top',
                   color='white', fontsize=9,
                   path_effects=[pe.withStroke(linewidth=2.0, foreground='black')])
            tick = np.round(2.0*x1f[-1]/3.0)   # two thirds of the half width
            a.set_xticks([-tick, 0, tick])
            a.set_yticks([-tick, 0, tick])
    for acol in axs.T:
        vis = [a for a in acol if a.get_visible()]
        vis[-1].tick_params(labelbottom=True)
        vis[-1].set_xlabel(r'$x$ [pc]', labelpad=1)
    for arow in axs:
        arow[0].set_ylabel(r'$y$ [pc]', labelpad=1)
    fig.subplots_adjust(left=0.045, right=0.92, bottom=0.075, top=0.955)
    fig.colorbar(ims['ion'], cax=fig.add_axes((0.93, 0.075, 0.012, 0.88)),
                 label='fraction of the element')
    for row, q in enumerate(LAYOUT[0]):
        p = axs[row, 0].get_position()
        cax = fig.add_axes((p.x1 + 0.006, p.y0, 0.008, p.height))
        fig.colorbar(ims[q], cax=cax).ax.tick_params(labelsize=6)
    fig.text(0.5, 0.985, title + r', $t = %.3f$ Myr' % (0.978*t), ha='center',
             va='top', fontsize=10)
    fig.savefig('%s/frame_%05d.png' % (outdir, n), dpi=150)
    plt.close(fig)
    print(os.path.basename(fn), 't = %.3f' % t, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('outdir')
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--frame', type=int, default=None)
    ap.add_argument('--molecular', action='store_true')
    ap.add_argument('--temp-log', action='store_true',
                    help='T on a log scale (10-2e4 K) in the pyathena colour map')
    ap.add_argument('--density-scale', type=float, default=1.0,
                    help='ambient n_H / 100 cm^-3, which scales the n_H and P/k ranges')
    ap.add_argument('--title', default=r'expanding H II region, SB99 2 Myr, '
                    r'$Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, no dust; '
                    r'Tigris ($N = 128^3$)')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    files = sorted(glob.glob(args.rundir + '/HII.out2.*.athdf'))
    jobs = [(n, fn, args.outdir, args.molecular, args.title, args.density_scale,
             args.temp_log)
            for n, fn in enumerate(files)]
    if args.frame is not None:
        frame(jobs[args.frame])
    else:
        with multiprocessing.Pool(args.nproc) as pool:
            pool.map(frame, jobs)
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '8', '-i',
                        args.outdir + '/frame_%05d.png', '-vf',
                        'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-pix_fmt', 'yuv420p',
                        '-vcodec', 'libx264', args.outdir + '/slices.mp4'], check=True)
        print(args.outdir + '/slices.mp4')
