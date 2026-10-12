"""The Wolfire, Hollenbach & McKee (2010) cloud (athinput.cloud_w10, run_w10.sh): the
radial structure against the visual extinction from the surface, and the CO-dark
molecular gas fraction. Writes ../../../figures/pdr_sphere_w10.png.

Their definitions (sec. 3): R_H2 is where 2 n(H2) = n(HI); R_CO where the CO column
from the surface reaches 2e16 cm^-2 (tau_CO = 1); A_V = N_H/1.9e21 cm^-2 at Z' = 1;
f_DG = [M(R_H2) - M(R_CO)]/M(R_H2) with M(R) the gas mass inside R, and their eq. 10,
f_DG = 1 - exp(-4 Delta A_V/A_V_bar), A_V_bar = 7.9 for N = 1.5e22 cm^-2. Profiles
are averages over shells one cell wide about the cloud centre; columns are integrated
inward along the radius of the averaged profile.

Beside the sphere, the plane-parallel slab W10 solve (make_slab_w10.py, run_slab_w10.sh)
at the mean density (f = 1) and at their clump density (f = 4), lit by one beam that
stands for theirs at 60 degrees: the normal A_V and N(CO) are half the beam's, and f_DG
takes the enclosed masses of the sphere shells at the normal depths of R_H2 and R_CO.
The slabs of 512 and 1024 cells are printed as a resolution check, not drawn.

The profiles come from reduced_w10.txt; `python plot.py --from-runs` reads the last
dumps of $PHOTCHEM_RUNS/T8_pdr_sphere_w10/{<mode>,slab_f<F>_n<NX>} and rewrites it.
Usage: python plot.py [--from-runs]
"""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from photchem_method import paths, reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'pdr_sphere_w10')
REDUCED = os.path.join(HERE, 'reduced_w10.txt')
PC, MH, MSUN, MU_H = 3.0857e18, 1.6735575e-24, 1.989e33, 1.4
AV_PER_NH = 1.0/1.9e21          # [mag cm^2] at Z' = 1, their N_0
N_CO_TAU1 = 2.0e16              # [cm^-2], their tau_CO = 1 column
X_C = 1.6e-4
CENTRE = 0.15625                # [pc], athinput.cloud_w10
R_OUT = 17.5                    # [pc], make_density_w10.py
AV_BAR = 5.26*1.5               # their eq. 3 for N = 1.5e22 cm^-2
R_C, N_COLUMN = 1.0, 1.5e22     # [pc], [cm^-2], make_density_w10.py and make_slab_w10.py
MODES = ['gow17']
# make_slab_w10.py factor and cells along the beam (run_slab_w10.sh); the first two
# are drawn, the others check the resolution
SLABS = [(1, 4096), (4, 4096), (1, 512), (4, 1024)]
# GOW17 species without ion ladders; the first is written as rHI
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+']


def profile(run):
    athena_read = paths.athena_read()
    fn = sorted(glob.glob(os.path.join(run, 'cloud.out2.*.athdf')))[-1]
    d = athena_read.athdf(fn)
    x = d['x1v'].astype(np.float64) - CENTRE
    z, y, xx = np.meshgrid(x, x, x, indexing='ij')
    r = np.sqrt(xx**2 + y**2 + z**2)
    dx = float(x[1] - x[0])
    edges = np.arange(0.0, R_OUT + dx, dx)
    ib = np.digitize(r.ravel(), edges) - 1
    keep = (ib >= 0) & (ib < len(edges) - 1)
    cnt = np.bincount(ib[keep], minlength=len(edges) - 1)

    def avg(q):
        s = np.bincount(ib[keep], weights=q.ravel()[keep].astype(np.float64),
                        minlength=len(edges) - 1)
        return s/np.maximum(cnt, 1)
    nh = d['rho']/MU_H
    sp = {nm: d['rHI' if n == 0 else 'r%d' % n] for n, nm in enumerate(NAMES)}
    x_h2, x_hp = sp['H2'], sp['H+']
    x_hi = np.clip(1.0 - 2.0*x_h2 - x_hp, 0.0, 1.0)
    out = {'r': 0.5*(edges[1:] + edges[:-1]), 'nH': avg(nh),
           'nH2': avg(nh*x_h2), 'nHI': avg(nh*x_hi), 'nCO': avg(nh*sp['CO']),
           'nCp': avg(nh*sp['C+']), 'T': avg(d['temp'])}
    return out, fn


def slab_profile(run):
    """Cells along the beam of a slab run: x [pc] and the same quantities as profile."""
    athena_read = paths.athena_read()
    fn = sorted(glob.glob(os.path.join(run, 'pdr_gow17.out2.*.athdf')))[-1]
    d = athena_read.athdf(fn)
    nh = d['rho'][0, 0].astype(np.float64)/MU_H
    sp = {nm: d['rHI' if n == 0 else 'r%d' % n][0, 0].astype(np.float64)
          for n, nm in enumerate(NAMES)}
    x_hi = np.clip(1.0 - 2.0*sp['H2'] - sp['H+'], 0.0, 1.0)
    out = {'x': d['x1v'].astype(np.float64), 'nH': nh, 'nH2': nh*sp['H2'],
           'nHI': nh*x_hi, 'nCO': nh*sp['CO'], 'nCp': nh*sp['C+'],
           'T': d['temp'][0, 0].astype(np.float64)}
    return out, fn


def from_runs():
    base = os.path.join(paths.env('PHOTCHEM_RUNS'), 'T8_pdr_sphere_w10')
    arrays, sources = {}, []
    for mode in MODES:
        p, fn = profile(os.path.join(base, mode))
        sources.append(fn)
        for k, v in p.items():
            arrays['%s_%s' % (mode, k)] = v
    for f, nx in SLABS:
        p, fn = slab_profile(os.path.join(base, 'slab_f%d_n%d' % (f, nx)))
        sources.append(fn)
        for k, v in p.items():
            arrays['slab_f%d_n%d_%s' % (f, nx, k)] = v
    reduced.save(REDUCED, 'shell averages about the cloud centre: r [pc], number '
                 'densities [cm^-3], T [K]; slab_f<F>_n<NX>: cells along the beam, '
                 'x [pc]', arrays, sources=sources)


def sphere_mass_inside(r):
    """Gas mass [Msun] inside radius r [pc] of the 1/r profile of make_density_w10.py."""
    rr = np.linspace(0.0, R_OUT, 4001)
    n = N_COLUMN/(2.0*np.sqrt(rr**2 + R_C**2)*PC)
    dm = 4*np.pi*(rr*PC)**2*n*MU_H*MH/MSUN
    m = np.concatenate([[0.0], np.cumsum(0.5*(dm[1:] + dm[:-1])*np.diff(rr)*PC)])
    return np.interp(r, rr, m)


def analyse_slab(f, x, nh, nh2, nhi, nco):
    """The W10 quantities of a slab: the beam crosses twice the normal column, so the
    normal A_V and N(CO) are half the beam's, and a layer at normal depth z = f x/2
    stands for the sphere shell at r_out - z, whose enclosed mass gives f_DG."""
    dx = (x[1] - x[0])*PC
    av = (np.cumsum(nh*dx) - 0.5*nh*dx)*AV_PER_NH/2.0
    nco_col = (np.cumsum(nco*dx) - 0.5*nco*dx)/2.0
    i_h2 = np.where(2*nh2 >= nhi)[0].min()
    i_co = np.where(nco_col >= N_CO_TAU1)[0].min()
    z_h2, z_co = f*x[i_h2]/2.0, f*x[i_co]/2.0
    m_h2, m_co = sphere_mass_inside(R_OUT - z_h2), sphere_mass_inside(R_OUT - z_co)
    a_h2, a_co = av[i_h2], av[i_co]
    return dict(av=av, av_h2=a_h2, av_co=a_co, f_mass=(m_h2 - m_co)/m_h2,
                f_eq10=1 - np.exp(-4*(a_co - a_h2)/AV_BAR))


def analyse(r, nh, nh2, nhi, nco):
    dr = (r[1] - r[0])*PC
    av = np.cumsum((nh*dr)[::-1])[::-1]*AV_PER_NH          # from the surface inward
    nco_col = np.cumsum((nco*dr)[::-1])[::-1]
    shell_mass = 4*np.pi*(r*PC)**2*dr*nh*MU_H*MH/MSUN
    m_in = np.cumsum(shell_mass)                           # gas mass inside r
    r_h2 = r[np.where(2*nh2 >= nhi)[0].max()]             # outermost molecular shell
    r_co = r[np.where(nco_col >= N_CO_TAU1)[0].max()]
    a_h2, a_co = np.interp(r_h2, r, av), np.interp(r_co, r, av)
    m_h2, m_co = np.interp(r_h2, r, m_in), np.interp(r_co, r, m_in)
    return dict(av=av, r_h2=r_h2, r_co=r_co, av_h2=a_h2, av_co=a_co,
                f_mass=(m_h2 - m_co)/m_h2, f_eq10=1 - np.exp(-4*(a_co - a_h2)/AV_BAR),
                m_co=m_co, m_h2=m_h2)


plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'xtick.major.pad': 2, 'ytick.major.pad': 1.5,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
                     'legend.fontsize': 8})
# Wolfire et al. (2010): the 1e5 Msun, G0' = 10 model (sec. 5.1) and the range of
# their grid (figs. 7 and 9)
W10 = {'av_h2': 0.23, 'av_co': 0.95, 'f_dg': 0.30}
W10_GRID = {'dav': (0.6, 0.8), 'f_dg': (0.25, 0.33)}
COL = {'H2': '#009E73', 'HI': '#0072B2', 'CO': '#CC79A7', 'C+': '#D55E00'}
# line style, marker and label of each geometry
GEOM = {'sphere': ('-', 'o', 'sphere'),
        'slab1': ('--', '^', r'slab, $f = 1$'),
        'slab4': (':', 'D', r'slab, $f = 4$')}


def main():
    if '--from-runs' in sys.argv:
        from_runs()
    rd = reduced.load(REDUCED)
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.0), constrained_layout=True)
    ax = axs.ravel()
    keys = ['nH', 'nH2', 'nHI', 'nCO', 'nCp', 'T']

    def fractions(p, av, g, ls, label):
        for q, y, lab in [('H2', 2*g['nH2']/g['nH'], r'2H$_2$/H'),
                          ('HI', g['nHI']/g['nH'], r'H$^0$/H'),
                          ('CO', g['nCO']/(X_C*g['nH']), 'CO/C'),
                          ('C+', g['nCp']/(X_C*g['nH']), r'C$^+$/C')]:
            p.plot(av, y, color=COL[q], lw=1.1, ls=ls, label=lab if label else None)

    def report(name, a):
        print('%-20s A_V(R_H2) = %.2f, A_V(R_CO) = %.2f, Delta A_V = %.2f, '
              'f_DG = %.2f (masses), %.2f (eq. 10)'
              % (name, a['av_h2'], a['av_co'], a['av_co'] - a['av_h2'],
                 a['f_mass'], a['f_eq10']))

    points = []
    for mode in MODES:
        g = {k: rd['%s_%s' % (mode, k)] for k in ['r'] + keys}
        a = analyse(g['r'], g['nH'], g['nH2'], g['nHI'], g['nCO'])
        report('sphere %s' % mode, a)
        fractions(ax[0], a['av'], g, '-', True)
        ax[2].plot(a['av'], g['T'], color='k', lw=1.1, ls=GEOM['sphere'][0],
                   label=GEOM['sphere'][2])
        for key, col in [('av_h2', COL['H2']), ('av_co', COL['CO'])]:
            ax[0].axvline(a[key], color=col, ls='-', lw=0.8)
        points.append(('sphere', a))
    for f, nx in SLABS:
        name = 'slab_f%d_n%d' % (f, nx)
        g = {k: rd['%s_%s' % (name, k)] for k in ['x'] + keys}
        a = analyse_slab(f, g['x'], g['nH'], g['nH2'], g['nHI'], g['nCO'])
        report(name, a)
        if (f, nx) not in SLABS[:2]:
            continue
        geom = 'slab%d' % f
        ls = GEOM[geom][0]
        fractions(ax[1], a['av'], g, ls, False)
        ax[2].plot(a['av'], g['T'], color='k', lw=1.1, ls=ls, label=GEOM[geom][2])
        for key, col in [('av_h2', COL['H2']), ('av_co', COL['CO'])]:
            ax[1].axvline(a[key], color=col, ls=ls, lw=0.8)
        points.append((geom, a))
    for p in ax[:2]:
        for key, col in [('av_h2', COL['H2']), ('av_co', COL['CO'])]:
            p.plot(W10[key], 1.0, 'v', color=col, ms=6, clip_on=False,
                   transform=p.get_xaxis_transform())
    for geom, a in points:
        dav = a['av_co'] - a['av_h2']
        mk = GEOM[geom][1]
        ax[3].plot(dav, a['f_mass'], mk, color='k', ms=5, mfc='none', mew=1.0)
        ax[3].plot(dav, a['f_eq10'], mk, color='k', ms=5, mew=1.0)
    x = np.linspace(0, 2.0, 100)
    ax[3].plot(x, 1 - np.exp(-4*x/AV_BAR), color='0.5', lw=0.8,
               label=r'eq. (10), $\bar A_V = 7.9$')
    (x0, x1), (y0, y1) = W10_GRID['dav'], W10_GRID['f_dg']
    ax[3].fill_between([x0, x1], y0, y1, color='0.85', lw=0, label='W10 grid')
    ax[3].plot(W10['av_co'] - W10['av_h2'], W10['f_dg'], 's', color='0.45', ms=5,
               label=r'W10, $10^5\,M_\odot$')
    for geom in ['sphere', 'slab1', 'slab4']:
        ax[3].plot([], [], GEOM[geom][1], color='k', ms=5, mew=1.0, label=GEOM[geom][2])
    ax[3].text(0.97, 0.03, 'open: masses; filled: eq. (10)', transform=ax[3].transAxes,
               ha='right', fontsize=7)
    for p in ax[:2]:
        p.set_yscale('log')
        p.set_ylim(1e-4, 2)
        p.set_ylabel('fraction of the element')
    ax[0].legend(loc='lower right', ncol=2, handlelength=1.4, columnspacing=0.8)
    ax[0].set_title('sphere, 128$^3$', loc='right', fontsize=8, pad=3)
    ax[1].set_title('slab: $f = 1$ dashed, $f = 4$ dotted', loc='right', fontsize=8,
                    pad=3)
    ax[1].plot([], [], 'v', color='0.3', ms=6,
               label=r'$R_{\rm H_2}$, $R_{\rm CO}$ of W10')
    ax[1].plot([], [], color='0.3', lw=0.8, label=r'$R_{\rm H_2}$, $R_{\rm CO}$, Tigris')
    ax[1].legend(loc='lower left', fontsize=7, handlelength=1.8)
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$T$ [K]')
    ax[2].legend(loc='lower left', fontsize=7, handlelength=2.2)
    for p in ax[:3]:
        p.set_xscale('log')
        p.set_xlim(1e-2, 15)
        p.set_xlabel(r'$A_V$ from the surface')
    ax[3].set_xlim(0, 2.0)
    ax[3].set_ylim(0, 0.8)
    ax[3].set_xlabel(r'$\Delta A_V = A_V(R_{\rm CO}) - A_V(R_{\rm H_2})$')
    ax[3].set_ylabel(r'$f_{\rm DG}$')
    ax[3].legend(loc='upper left', fontsize=7, ncol=2, columnspacing=0.8)
    for p, lab in zip(ax, ['(a)', '(b)', '(c)', '(d)']):
        p.set_title(lab, loc='left', fontsize=9, pad=3)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for ext in figure_formats():
        fig.savefig(OUT + '.' + ext, dpi=250)
    print(OUT)


if __name__ == '__main__':
    main()
