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
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
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
LAYOUT = [['nH', 'T', 'Pk', 'v'], ['H0', 'H+', None, None],
          ['He0', 'He+', 'He++', None]] + \
         [[el + '0', el + '+', el + '++', el + '_high'] for el in 'OSN']
LAB = {'nH': r'$n_{\rm H}$ [cm$^{-3}$]', 'T': r'$T$ [K]', 'Pk': r'$P/k$ [K cm$^{-3}$]',
       'v': r'$|v|$ [km s$^{-1}$]', 'H0': r'H$^0$/H', 'H+': r'H$^+$/H',
       'He0': r'He$^0$/He', 'He+': r'He$^+$/He', 'He++': r'He$^{2+}$/He'}
for el in 'OSN':
    LAB.update({el + '0': r'%s$^0$/%s' % (el, el), el + '+': r'%s$^+$/%s' % (el, el),
                el + '++': r'%s$^{2+}$/%s' % (el, el),
                el + '_high': r'%s$^{\geq3+}$/%s' % (el, el)})
NORM = {'nH': LogNorm(10, 1e3), 'T': Normalize(5e3, 1.1e4), 'Pk': LogNorm(1e3, 3e6),
        'v': Normalize(0, 15), 'ion': LogNorm(1e-5, 1)}
CMAP = {'nH': 'cividis', 'T': 'inferno', 'Pk': 'magma', 'v': 'plasma',
        'ion': plt.get_cmap('viridis').copy()}
CMAP['ion'].set_bad(CMAP['ion'](0.0))   # a fraction of exactly 0 takes the lowest colour


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
    for el, ions in IONS.items():
        for nm in ions:
            f[nm] = x[nm]/XTOT[el]
        f[el + '0'] = np.maximum(1.0 - sum(f[nm] for nm in ions), 0.0)
    return f, d['x1f'], d['x2f'], float(d['Time'])


plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
NCOL = len(LAYOUT)


def frame(job):
    """Draw dump fn as OUTDIR/frame_<n>.png."""
    n, fn, outdir = job
    f, x1f, x2f, t = fields(fn)
    fig, allax = plt.subplots(4, NCOL + 1, figsize=(10.0, 6.1), sharex=True,
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
            ims[kind] = a.pcolormesh(x1f, x2f, f[q], norm=NORM[kind], cmap=CMAP[kind],
                                     shading='flat', rasterized=True)
            a.set_aspect('equal')
            a.text(0.05, 0.94, LAB[q], transform=a.transAxes, ha='left', va='top',
                   color='white', fontsize=9,
                   path_effects=[pe.withStroke(linewidth=2.0, foreground='black')])
            a.set_xticks([-8, 0, 8])
            a.set_yticks([-8, 0, 8])
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
    fig.text(0.5, 0.985, r'expanding H II region, SB99 2 Myr, $Q = 10^{49}$ s$^{-1}$, '
             r'$n_{\rm H} = 100$ cm$^{-3}$, no dust; Tigris ($N = 128^3$), '
             r'$t = %.3f$ Myr' % (0.978*t), ha='center', va='top', fontsize=10)
    fig.savefig('%s/frame_%05d.png' % (outdir, n), dpi=150)
    plt.close(fig)
    print(os.path.basename(fn), 't = %.3f' % t, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('outdir')
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--frame', type=int, default=None)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    files = sorted(glob.glob(args.rundir + '/HII.out2.*.athdf'))
    jobs = [(n, fn, args.outdir) for n, fn in enumerate(files)]
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
