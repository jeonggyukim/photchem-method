"""F15: slices and profiles of the Bonnor-Ebert sphere in an isotropic FUV field, one
figure per network (NCR, GOW17 + ions), same layout as
ai-notes/docs-code/tigris-rayt/scripts/postproc_background_slices.py.

    python plot.py [RUNS]    (default ~/Documents/photchem-postproc/runs/be_net)

RUNS/ncr and RUNS/gow17 hold the tab slices out3 (prim) and out4 (uov) of run.sh.
Writes ../../../figures/F15_be_sphere_{ncr,gow17}.{pdf,png}.

GOW17 writes its species unnamed (r0 ... r21, the first labelled rHI by the NCR output
code, though it is He+); they are mapped from the species order below.  x(H) and x(e)
are derived: x(H) from the hydrogen budget, x(e) from the charges of the ions, with each
element's higher stages counted at their lowest charge (negligible without ionizing light).
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', '..', 'figures')
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/Documents/photchem-postproc/runs/be_net')
X_SRC = 0.046875
XC_TOT = 1.6e-4
XI_CR0 = 2.0e-16
R_CLOUD = 2.001
CLOUD = ('A Bonnor-Ebert sphere of central density 1000 H cm$^{-3}$ cut off at 2.00 pc '
         '($\\xi_{\\rm max} = 8.4$) in gas of 0.05 H cm$^{-3}$')
DUST = {'ncr': ', dust temperature solved and coupled to the gas',
        'gow17': ', no dust temperature (GOW17), CO self-shielding on'}

GOW17 = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}


def read_plane(pattern):
    data, names = [], None
    for fname in sorted(glob.glob(pattern)):
        with open(fname) as f:
            f.readline()
            names = f.readline().strip('# \n').split()
        data.append(np.atleast_2d(np.loadtxt(fname)))
    data = np.vstack(data)
    x = np.unique(data[:, names.index('x1v')])
    y = np.unique(data[:, names.index('x2v')])
    i = np.searchsorted(x, data[:, names.index('x1v')])
    j = np.searchsorted(y, data[:, names.index('x2v')])
    out = {}
    for n, name in enumerate(names):
        out[name] = np.full((y.size, x.size), np.nan)
        out[name][j, i] = data[:, n]
    return x, y, out


def gow17_species(p):
    """Name the GOW17 scalars and add the derived x(H) and x(e) under NCR's names."""
    cols = ['rHI'] + ['r%d' % k for k in range(1, len(GOW17))]
    s = {sp: p[c] for sp, c in zip(GOW17, cols)}
    q = dict(p)
    q['rH2'] = s['H2']
    q['rHI'] = (1.0 - 2.0*s['H2'] - s['H+'] - 2.0*s['H2+'] - 3.0*s['H3+'] - s['CHx']
                - s['OHx'] - s['HCO+'])
    q['rEL'] = sum(z*s[sp] for sp, z in CHARGE.items())
    return q


def decade_limits(lo, hi):
    return 10.0**np.floor(np.log10(lo)), 10.0**np.ceil(np.log10(hi))


def nice_limits(lo, hi):
    if hi <= lo:
        return lo - 1.0, hi + 1.0
    raw = (hi - lo)/5.0
    mag = 10.0**np.floor(np.log10(raw))
    step = mag*min(s for s in (1.0, 2.0, 5.0, 10.0) if s*mag >= raw)
    return step*np.floor(lo/step + 1e-9), step*np.ceil(hi/step - 1e-9)


COLS = (
    (r'$n_{\rm H}$ [cm$^{-3}$]', lambda p, u: p['rho']/1.4, 'Greys', 'log'),
    (r'$x_{\rm e}$', lambda p, u: p['rEL'], 'viridis', 'log'),
    ('2 x(H2)', lambda p, u: 2.0*p['rH2'], 'viridis', (0.0, 1.0)),
    ('x(C+)/x(C,tot)', lambda p, u: u['xCII']/XC_TOT, 'viridis', (0.0, 1.0)),
    ('x(C)/x(C,tot)', lambda p, u: u['xCI']/XC_TOT, 'viridis', (0.0, 1.0)),
    ('x(CO)/x(C,tot)', lambda p, u: u['xCO']/XC_TOT, 'viridis', (0.0, 1.0)),
    ('T [K]', lambda p, u: u['temp'], 'magma', 'cloud'),
    (r'$T_{\rm d}$ [K]', lambda p, u: u['temp_dust'], 'magma', 'cloud'),
    (r'$\chi_{\rm PE}$ [Draine]', lambda p, u: u['chi_PE'], 'inferno', 'log'),
    (r'$\chi_{\rm H_2}$ [Draine]', lambda p, u: u['chi_H2'], 'inferno', 'log'),
    (r'$\chi_{\rm C\,I}$ [Draine]', lambda p, u: u['chi_CI'], 'inferno', 'log'),
    # the field that photodissociates CO, with dust and CO self-shielding (GOW17 only;
    # runs from before the output existed have no chi_CO column)
    (r'$\chi_{\rm CO}$ [Draine]', lambda p, u: u.get('chi_CO', 0.0*u['temp']), 'inferno',
     'log'),
)
CUTS = (
    (r'$n_{\rm H}$ [cm$^{-3}$]', 'log', None,
     ((r'$n_{\rm H}$', lambda p, u: p['rho']/1.4),)),
    ('temperature [K]', 'log', 'cloud_tight',
     (('gas', lambda p, u: u['temp']), ('dust', lambda p, u: u['temp_dust']))),
    ('hydrogen and electrons', 'log', 'cloud',
     ((r'$x_{\rm e}$', lambda p, u: p['rEL']), ('x(H)', lambda p, u: p['rHI']),
      ('2 x(H2)', lambda p, u: 2.0*p['rH2']))),
    ('carbon, fraction of x(C,tot)', 'log', 'cloud',
     (('C+', lambda p, u: u['xCII']/XC_TOT), ('C', lambda p, u: u['xCI']/XC_TOT),
      ('CO', lambda p, u: u['xCO']/XC_TOT, False))),
    ('radiation field [Draine]', 'log', 'cloud',
     ((r'$\chi_{\rm PE}$', lambda p, u: u['chi_PE']),
      (r'$\chi_{\rm LW}$', lambda p, u: u['chi_LW']),
      (r'$\chi_{\rm H_2}$', lambda p, u: u['chi_H2']),
      (r'$\chi_{\rm C\,I}$', lambda p, u: u['chi_CI']),
      (r'$\chi_{\rm CO}$', lambda p, u: u.get('chi_CO', 0.0*u['temp'])))),
    (r'$\xi_{\rm cr}/\xi_{\rm cr,0}$', 'linear', 'line',
     ((r'$\xi_{\rm cr}/\xi_{\rm cr,0}$', lambda p, u: u['xi_CR']/XI_CR0),)),
)


def figure(net, label):
    run_dir = os.path.join(RUNS, net)
    dumps = sorted(set(f.split('.')[-2] for f in glob.glob(
        os.path.join(run_dir, 'cloud.block*.out4.*.tab'))))
    runlog = os.path.join(run_dir, 'run.out')
    if len(dumps) < 2 or 'cpu time used' not in open(runlog).read():
        print('skipped', net, '(no finished run in %s)' % run_dir)
        return
    x, y, p = read_plane(os.path.join(run_dir, 'cloud.block*.out3.%s.tab' % dumps[-1]))
    _, _, u = read_plane(os.path.join(run_dir, 'cloud.block*.out4.%s.tab' % dumps[-1]))
    if net == 'gow17':
        p = gow17_species(p)
    ncol = (len(COLS) + 1)//2
    fig, grid = plt.subplots(3, ncol, figsize=(3.3*ncol, 3.1*3), squeeze=False)
    axes = grid.reshape(1, 3*ncol)
    dx = x[1] - x[0]
    ext = (x[0] - 0.5*dx - X_SRC, x[-1] + 0.5*dx - X_SRC,
           y[0] - 0.5*dx - X_SRC, y[-1] + 0.5*dx - X_SRC)
    jcut = int(np.argmin(np.abs(y - X_SRC)))
    xx, yy = np.meshgrid(x - X_SRC, y - X_SRC)
    inside = np.hypot(xx, yy) < R_CLOUD - dx
    for c, (lab, func, cmap, scale) in enumerate(COLS):
        ax = axes[0, c]
        q = func(p, u)
        if not np.any(q):
            # a quantity this network does not compute (T_dust under GOW17)
            ax.text(0.5, 0.5, 'not computed\nby %s' % label.split()[0], ha='center',
                    va='center', transform=ax.transAxes, fontsize=10, color='0.4')
            ax.set_xticks([]); ax.set_yticks([]); ax.set_title(lab, fontsize=10)
            continue
        if isinstance(scale, tuple):
            vmin, vmax = scale[:2]
            log = False
        elif scale == 'log' and np.nanmax(q) > 30.0*np.nanmin(q):
            vmin, vmax = decade_limits(np.nanmin(q), np.nanmax(q))
            log = True
        else:
            sel = inside if scale == 'cloud' else np.isfinite(q)
            vmin, vmax = nice_limits(np.nanmin(q[sel]), np.nanmax(q[sel]))
            log = False
        kw = {'norm': LogNorm(vmin, vmax)} if log else {'vmin': vmin, 'vmax': vmax}
        im = ax.imshow(q, origin='lower', extent=ext, cmap=cmap, **kw)
        if c == 0:
            ax.axhline(y[jcut] - X_SRC, c='C1', lw=0.8, ls='--')
        ax.set_xticks([-2, 0, 2]); ax.set_yticks([-2, 0, 2]); ax.tick_params(labelsize=7)
        ax.set_title(lab, fontsize=10)
        if c % ncol == 0:
            ax.set_ylabel(('isotropic FUV, 1 Draine\n\n' if c == 0 else '') + 'x2 [pc]',
                          fontsize=9)
        if c >= ncol:
            ax.set_xlabel('x1 [pc]', fontsize=9)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6, which='both')
    for c, (lab, yscale, ylim, curves) in enumerate(CUTS):
        ax = axes[0, 2*ncol + c]
        lo, hi = np.inf, -np.inf
        sel = inside[jcut] if ylim in ('cloud', 'cloud_tight') else np.ones(x.size, bool)
        for curve in curves:
            q = curve[1](p, u)[jcut]
            if not np.any(q):
                continue
            ax.plot(x - X_SRC, q, lw=1.2, label=curve[0])
            if curve[2:] != (False,):
                lo = min(lo, np.nanmin(q[sel])); hi = max(hi, np.nanmax(q[sel]))
        ax.set_yscale(yscale)
        if ylim == 'cloud_tight':
            ax.set_ylim(lo/1.25, 1.25*hi)
        elif ylim in ('cloud', 'line') and yscale == 'log':
            lo, hi = decade_limits(lo, hi); ax.set_ylim(lo, 2.0*hi)
        elif ylim in ('cloud', 'line'):
            lo, hi = nice_limits(lo, hi)
            ax.set_ylim(lo - 0.04*(hi - lo), hi + 0.04*(hi - lo))
        ax.set_xlim(ext[0], ext[1]); ax.set_xticks([-2, 0, 2])
        ax.axvspan(-R_CLOUD, R_CLOUD, color='0.9', zorder=0)
        ax.tick_params(labelsize=7); ax.set_title(lab, fontsize=10)
        ax.set_xlabel('x1 [pc]', fontsize=9); ax.grid(alpha=0.3)
        if len(curves) > 1:
            ax.legend(fontsize=7, loc='best', framealpha=0.8)
    text = open(runlog).read()
    ndir = re.search(r'on (\d+) directions', text)
    ndir = ndir.group(1) if ndir else '48'
    m = re.search(r'converged after (\d+) iterations', text)
    ended = ('%s post-processed to convergence (%s iterations)' % (label, m.group(1)) if m
             else '%s post-processed, NOT converged' % label)
    fig.suptitle('%s%s, plane through its middle\nisotropic background through the diffuse '
                 'solver (HEALPix, %s directions), %s' % (CLOUD, DUST[net], ndir, ended),
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    for ext_ in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG, 'F15_be_sphere_%s.%s' % (net, ext_)), dpi=150)
    plt.close(fig)
    print('wrote', os.path.join(FIG, 'F15_be_sphere_%s.png' % net))


if __name__ == '__main__':
    figure('ncr', 'NCR')
    figure('gow17', 'GOW17 + ions (O3,S3,N3)')
