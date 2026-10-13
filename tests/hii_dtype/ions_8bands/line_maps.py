"""Optically thin line maps of an hii_dtype run (GOW17 + O3,S3,N3): H-alpha, [O III]
5007 and [S II] 6716+6731 surface brightness, int j dl along z, and their composite in
the narrowband palette of Hubble images, [S II] red, H-alpha green, [O III] blue. No
dust extinction or scattering; this is the numpy reference the renderer is checked
against.

H-alpha is the case B fit of Tigris emission.hpp (Draine 2011 eq. 14.8); the forbidden
lines come from the 5-level populations of pyathena's IonCoolant (CHIANTI 11), with
4 pi j = n_ion f_u A_ul h nu_ul. The ism-unit dump holds rho in m_H cm^-3.

usage: python line_maps.py RUNDIR OUTDIR [--frame I ...] [--nproc N] [--title TEXT]
Writes OUTDIR/lines_NNNNN.png per HII.out2 dump (all, or those of --frame) and, for the
whole run, OUTDIR/lines.mp4 when ffmpeg is found. Each line keeps one stretch over all
frames drawn.
needs pyathena (the Dropbox checkout) on PYTHONPATH.
"""
import argparse
import glob
import multiprocessing
import os
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from pyathena.chemistry.coolants.base import IonCoolant  # noqa: E402
import slice_movie as sm  # noqa: E402

PC = 3.0857e18
E_HALPHA = 6.62607015e-27*2.99792458e10/6562.80e-8
# label, ion file, abundance key, (upper, lower) level pairs summed
FORBIDDEN = {'OIII': ('o_3.txt', 'O++', [(3, 2)]),
             'SII': ('s_2.txt', 'S+', [(1, 0), (2, 0)])}
LINES = ['SII', 'Halpha', 'OIII']           # red, green, blue
LABEL = {'Halpha': r'H$\alpha$', 'OIII': r'[O III] 5007',
         'SII': r'[S II] 6716+6731'}
ION = {}


def ion(fname):
    if fname not in ION:
        ION[fname] = IonCoolant(fname, fname[:-4])
    return ION[fname]


def surface_brightness(fn):
    """int j dl along x3 [erg s^-1 cm^-2 sr^-1] of every line."""
    d = sm.athena_read.athdf(fn)
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(sm.NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(sm.NAMES, keys)}
    nh = d['rho'].astype(np.float64)/sm.MU_H
    xe = sum(sm.CHARGE.get(nm, 0)*x[nm] for nm in sm.NAMES)
    T = d['press'].astype(np.float64)*sm.MH*1e10/sm.KB/(nh*(1.0 - x['H2'] + sm.X_HE + xe))
    ne = xe*nh
    dz = float(d['x3f'][1] - d['x3f'][0])*PC
    m = (ne > 1e-3*nh) & (T > 0)
    out = {}
    j = np.zeros_like(nh)
    t4 = T[m]*1e-4
    j[m] = (E_HALPHA*1.17e-13*t4**(-0.942 - 0.031*np.log(t4))*ne[m]*x['H+'][m]*nh[m])
    out['Halpha'] = j.sum(axis=0)*dz/(4*np.pi)
    for line, (fname, key, pairs) in FORBIDDEN.items():
        c = ion(fname)
        tab = c.table()
        f = c.populations(T[m], ne[m], x['H+'][m]*nh[m])
        j = np.zeros_like(nh)
        E, A = tab['E_erg'], tab['A']
        j[m] = x[key][m]*nh[m]*sum(f[u]*A[u, lo]*(E[u] - E[lo]) for u, lo in pairs)
        out[line] = j.sum(axis=0)*dz/(4*np.pi)
    return out, d['x1f'], d['x2f'], float(d['Time'])


def stretch(img, lo, hi):
    """Logarithmic stretch of img onto 0-1 between lo and hi."""
    return np.clip(np.log10(np.maximum(img, lo)/lo)/np.log10(hi/lo), 0.0, 1.0)


def frame(job):
    n, fn, outdir, title, rng = job
    sb, x1f, x2f, t = surface_brightness(fn)
    fig, ax = plt.subplots(1, 4, figsize=(13.0, 3.6), constrained_layout=True)
    for a, line in zip(ax[:3], ['Halpha', 'OIII', 'SII']):
        lo, hi = rng[line]
        im = a.pcolormesh(x1f, x2f, np.maximum(sb[line], lo), norm=LogNorm(lo, hi),
                          cmap='magma', shading='flat')
        a.text(0.04, 0.95, LABEL[line], transform=a.transAxes, va='top', color='w',
               fontsize=10, path_effects=[pe.withStroke(linewidth=2, foreground='k')])
        cb = fig.colorbar(im, ax=a, orientation='horizontal', fraction=0.05, pad=0.02)
        cb.set_label(r'$\int j\,dl$ [erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$]', fontsize=8)
        cb.ax.tick_params(labelsize=7)
    rgb = np.dstack([stretch(sb[line], *rng[line]) for line in LINES])
    ax[3].imshow(rgb, origin='lower', extent=[x1f[0], x1f[-1], x2f[0], x2f[-1]],
                 interpolation='nearest')
    ax[3].text(0.04, 0.95, '[S II] H$\\alpha$ [O III]', transform=ax[3].transAxes,
               va='top', color='w', fontsize=10)
    for a in ax:
        a.set_aspect('equal')
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle('%s, $t = %.3f$ Myr, optically thin, along $z$' % (title, t*0.978),
                 fontsize=11)
    out = os.path.join(outdir, 'lines_%05d.png' % n)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(os.path.basename(fn), 't = %.3f' % t)


def ranges(fns):
    """One stretch per line over the frames drawn: from 1e-3.5 of the brightest pixel
    to the brightest."""
    hi = {line: 0.0 for line in LINES}
    for fn in fns:
        sb = surface_brightness(fn)[0]
        for line in LINES:
            hi[line] = max(hi[line], float(sb[line].max()))
    return {line: (hi[line]*10**-3.5, hi[line]) for line in LINES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('outdir')
    ap.add_argument('--frame', type=int, nargs='*', default=None)
    ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fns = sorted(glob.glob(os.path.join(a.rundir, 'HII.out2.*.athdf')))
    pick = list(range(len(fns))) if a.frame is None else a.frame
    rng = ranges([fns[n] for n in pick])
    print('stretch [erg s^-1 cm^-2 sr^-1]:',
          {k: '%.3g-%.3g' % v for k, v in rng.items()})
    jobs = [(n, fns[n], a.outdir, a.title, rng) for n in pick]
    if a.frame is not None:
        for job in jobs:
            frame(job)
        return
    with multiprocessing.Pool(a.nproc) as p:
        p.map(frame, jobs)
    if shutil.which('ffmpeg'):
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '8', '-i',
                        os.path.join(a.outdir, 'lines_%05d.png'), '-pix_fmt', 'yuv420p',
                        '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-vcodec', 'libx264',
                        os.path.join(a.outdir, 'lines.mp4')], check=True)


if __name__ == '__main__':
    main()
