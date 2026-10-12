"""Which cells set the chemistry timestep dt_user of an hii_dtype run (GOW17 or NCR,
built after the per-cell timescale outputs t_cool, t_Hp, t_H2, t_CO were added).

dt_user = cfl_photchem * min over cells and kinds of the start-of-step timescales.
For one dump this lists the cells with the shortest timescales, with their state, and
draws the z slice through the shortest one: n_H, T, the shortest timescale of each
cell and the kind that sets it.

usage: python timestep_cells.py RUNDIR DUMP [--ncell N] [--out PNG]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm  # noqa: E402
import slice_movie as sm  # noqa: E402

KINDS = ['t_cool', 't_Hp', 't_H2', 't_CO']
LABELS = ['cooling', r'H$^+$', r'H$_2$', 'CO']
COLORS = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd']


def load(rundir, dump):
    fn = os.path.join(rundir, 'HII.out2.%05d.athdf' % dump)
    d = sm.athena_read.athdf(fn)
    x = {nm: d[key] for nm, key in
         zip(sm.NAMES, ['rHI'] + ['r%d' % n for n in range(1, len(sm.NAMES))])}
    nh = d['rho']/sm.MU_H
    xe = sum(sm.CHARGE.get(nm, 0)*x[nm] for nm in sm.NAMES)
    pk = d['press']*sm.MH*1e10/sm.KB
    f = {'nH': nh, 'T': pk/(nh*(1.0 - x['H2'] + sm.X_HE + xe)), 'xHp': x['H+'],
         'xH2': x['H2'], 'xCO': x['CO']}
    ts = np.stack([d[k] for k in KINDS])
    return d, f, ts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('dump', type=int)
    ap.add_argument('--ncell', type=int, default=15)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    d, f, ts = load(a.rundir, a.dump)
    tmin = ts.min(axis=0)
    kind = ts.argmin(axis=0)
    x3, x2, x1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    r = np.sqrt(x1**2 + x2**2 + x3**2)
    print('t = %.4f [code]; shortest timescale per kind [code time]:' % d['Time'])
    for n, k in enumerate(KINDS):
        print('  %-7s %.3e' % (k, ts[n].min()))
    print('cells setting the shortest timescale, by kind: '
          + ', '.join('%s %d' % (k, (kind == n).sum()) for n, k in enumerate(KINDS)))
    order = np.argsort(tmin, axis=None)[:a.ncell]
    print('\n%-4s %-9s %-8s %-6s %-9s %-9s %-7s %-9s %-9s' %
          ('rank', 't_min', 'kind', 'r[pc]', 'nH', 'T[K]', 'xH+', 'xH2', 'xCO/xC'))
    for m, flat in enumerate(order):
        c = np.unravel_index(flat, tmin.shape)
        print('%-4d %-9.3e %-8s %-6.1f %-9.3g %-9.3g %-7.3f %-9.2e %-9.2e' %
              (m, tmin[c], KINDS[kind[c]], r[c], f['nH'][c], f['T'][c], f['xHp'][c],
               f['xH2'][c], f['xCO'][c]/sm.X_C))
    k0 = np.unravel_index(order[0], tmin.shape)[0]
    ext = [d['x1f'][0], d['x1f'][-1], d['x2f'][0], d['x2f'][-1]]
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.8), constrained_layout=True)
    panels = [(f['nH'][k0], r'$n_{\rm H}$ [cm$^{-3}$]', LogNorm(0.3, 1e3), 'cividis'),
              (f['T'][k0], r'$T$ [K]', LogNorm(10, 2e4), 'RdYlBu_r'),
              (tmin[k0], r'shortest timescale [code]',
               LogNorm(tmin.min(), tmin.min()*1e3), 'magma'),
              (kind[k0], 'kind of the shortest', BoundaryNorm(np.arange(-0.5, 4), 4),
               ListedColormap(COLORS))]
    for p, (z, lab, norm, cmap) in zip(ax, panels):
        im = p.imshow(z, origin='lower', extent=ext, norm=norm, cmap=cmap)
        p.set_title(lab, fontsize=9)
        p.set_xlabel('x [pc]')
        cb = fig.colorbar(im, ax=p, shrink=0.85)
        if p is ax[3]:
            cb.set_ticks(range(4))
            cb.set_ticklabels(LABELS)
    sel = np.argsort(tmin[k0], axis=None)[:a.ncell]
    jj, ii = np.unravel_index(sel, tmin[k0].shape)
    for p in ax[:3]:
        p.plot(d['x1v'][ii], d['x2v'][jj], 'o', mfc='none', mec='w', ms=6, mew=1.0)
    ax[0].set_ylabel('y [pc]')
    fig.suptitle('t = %.3f code, slice z = %.1f pc through the cell with the '
                 'shortest timescale (open circles: %d shortest in the slice)'
                 % (d['Time'], d['x3v'][k0], a.ncell), fontsize=9)
    out = a.out or os.path.join(a.rundir, 'timestep_cells_%05d.png' % a.dump)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == '__main__':
    main()
