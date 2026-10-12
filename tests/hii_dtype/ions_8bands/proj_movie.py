"""Projection movie of an hii_dtype run (GOW17 + O3,S3,N3, run_tigris.sh or the grammar
jobs): along x3, the density-weighted temperature, the columns of hydrogen nuclei, H2
and CO, and the emission measure of the ionized gas, int n_e^2 dl, beside the mass in
H+, H0, H2 and CO over the whole run (HII.hst) with a line at the time of the frame.
Every panel keeps one colour range over all frames.

usage: python proj_movie.py RUNDIR OUTDIR [--nproc N] [--frame I] [--title TEXT]
Writes OUTDIR/proj_NNNNN.png, one per HII.out2 dump, then OUTDIR/proj.mp4 when ffmpeg
is found. --frame I draws only dump I.
"""
import argparse
import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import slice_movie as sm  # noqa: E402

PC = 3.0857e18
PANELS = [('T', r'$\langle T\rangle_\rho$ [K]', LogNorm(10, 1e4), 'RdYlBu_r'),
          ('NH', r'$N_{\rm H}$ [cm$^{-2}$]', LogNorm(5e21, 1e23), 'cividis'),
          ('NH2', r'$N_{\rm H_2}$ [cm$^{-2}$]', LogNorm(1e21, 5e22), 'viridis'),
          ('NCO', r'$N_{\rm CO}$ [cm$^{-2}$]', LogNorm(1e16, 1e19), 'magma'),
          ('EM', r'EM [cm$^{-6}$ pc]', LogNorm(1e3, 1e7), 'inferno')]
# ism units: the code density is m_H cm^-3 and the code length a pc
CODE_MASS_MSUN = 1.6735575e-24*PC**3/1.989e33
MYR = 0.978


def mass_history(rundir):
    """Mass [Msun] in H+, H0, H2 (hydrogen nuclei times mu_H) and CO (28 m_H per
    molecule) against time [Myr], from the passive-scalar sums of HII.hst, which are
    integrals of rho x over the box."""
    fn = os.path.join(rundir, 'HII.hst')
    hdr = [line for line in open(fn) if line.startswith('#')][-1]
    names = [m.split('=')[1] for m in re.findall(r'\[\d+\]=\S+', hdr)]
    a = np.atleast_2d(np.loadtxt(fn))
    col = {nm: a[:, names.index('%d-scalar' % n)] for n, nm in enumerate(sm.NAMES)}
    total = a[:, names.index('mass')]
    m = CODE_MASS_MSUN
    hp, h2 = col['H+']*m, 2.0*col['H2']*m
    return {'t': a[:, 0]*MYR, 'H+': hp, 'H2': h2, 'H0': total*m - hp - h2,
            'CO': col['CO']*m*28.0/sm.MU_H}


def project(fn):
    d = sm.athena_read.athdf(fn)
    x = {nm: d[key] for nm, key in
         zip(sm.NAMES, ['rHI'] + ['r%d' % n for n in range(1, len(sm.NAMES))])}
    nh = d['rho'].astype(np.float64)/sm.MU_H
    xe = sum(sm.CHARGE.get(nm, 0)*x[nm] for nm in sm.NAMES)
    pk = d['press']*sm.MH*1e10/sm.KB
    temp = pk/(nh*(1.0 - x['H2'] + sm.X_HE + xe))
    dz = float(d['x3f'][1] - d['x3f'][0])
    f = {'NH': nh.sum(axis=0)*dz*PC,
         'NH2': (nh*x['H2']).sum(axis=0)*dz*PC,
         'NCO': (nh*x['CO']).sum(axis=0)*dz*PC,
         'EM': ((nh*xe)**2).sum(axis=0)*dz,
         'T': (nh*temp).sum(axis=0)/nh.sum(axis=0)}
    return f, d['x1f'], d['x2f'], float(d['Time'])


def frame(job):
    n, fn, outdir, title, hist = job
    f, x1f, x2f, t = project(fn)
    fig, ax = plt.subplots(1, len(PANELS) + 1, figsize=(3.0*(len(PANELS) + 1.2), 3.25),
                           gridspec_kw={'wspace': 0.0,
                                        'width_ratios': [1]*len(PANELS) + [1.3]},
                           constrained_layout=False)
    fig.subplots_adjust(left=0.02, right=0.99, bottom=0.2, top=0.88)
    h = ax[-1]
    for q, col in [('H+', 'C3'), ('H0', 'C0'), ('H2', 'C2'), ('CO', 'C4')]:
        h.semilogy(hist['t'], hist[q], color=col, lw=1.0, label=q)
    h.axvline(t*MYR, color='k', lw=0.6)
    h.set_xlim(0, hist['t'][-1])
    lo = min(np.min(hist[q][hist[q] > 0]) for q in ['H+', 'H0', 'H2', 'CO'])
    hi = max(np.max(hist[q]) for q in ['H+', 'H0', 'H2', 'CO'])
    h.set_ylim(0.5*lo, 2.0*hi)
    h.set_xlabel('t [Myr]', fontsize=8)
    h.yaxis.tick_right()
    h.yaxis.set_label_position('right')
    h.set_ylabel(r'mass [M$_\odot$]', fontsize=8)
    h.tick_params(labelsize=7)
    h.legend(fontsize=7, loc='lower right', frameon=False)
    h.set_box_aspect(1)
    for a, (q, lab, norm, cmap) in zip(ax[:-1], PANELS):
        cm = plt.get_cmap(cmap).copy()
        cm.set_bad(cm(0.0))
        z = np.ma.masked_invalid(np.where(f[q] > 0, f[q], np.nan))
        im = a.pcolormesh(x1f, x2f, z, norm=norm, cmap=cm, shading='flat')
        a.set_aspect('equal')
        a.set_xticks([])
        a.set_yticks([])
        a.text(0.04, 0.95, lab, transform=a.transAxes, va='top', color='w', fontsize=9,
               path_effects=[pe.withStroke(linewidth=2, foreground='k')])
        cb = fig.colorbar(im, ax=a, orientation='horizontal', fraction=0.05, pad=0.03)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle('%s, $t = %.3f$ Myr, projected along $z$' % (title, t*0.978),
                 fontsize=10)
    out = os.path.join(outdir, 'proj_%05d.png' % n)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(os.path.basename(fn), 't = %.3f' % t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('outdir')
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--frame', type=int, default=None)
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fns = sorted(glob.glob(os.path.join(a.rundir, 'HII.out2.*.athdf')))
    hist = mass_history(a.rundir)
    jobs = [(n, fn, a.outdir, a.title, hist) for n, fn in enumerate(fns)]
    if a.frame is not None:
        frame(jobs[a.frame])
        return
    with multiprocessing.Pool(a.nproc) as p:
        p.map(frame, jobs)
    if shutil.which('ffmpeg'):
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '8', '-i',
                        os.path.join(a.outdir, 'proj_%05d.png'), '-pix_fmt', 'yuv420p',
                        '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-vcodec', 'libx264',
                        os.path.join(a.outdir, 'proj.mp4')], check=True)


if __name__ == '__main__':
    main()
