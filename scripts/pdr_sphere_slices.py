"""The post-processed sphere: slices through its middle and profiles along x1.

    python pdr_sphere_slices.py [npz ...]

Reads the reduced mid-plane slices that reproduce/pdr_sphere/reduce.py writes
(default: the NCR and GOW17 ones in ../data/pdr_sphere, beside this repository) and writes figures/pdr_sphere_<mode>.{pdf,png}.
A uniform sphere of 1000 H cm^-3, 2 pc in radius, in gas of 0.05 H cm^-3, lit from
every side by one Draine field through the diffuse solver.
"""
import ast
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REDUCED = os.path.join(HERE, '..', '..', 'data', 'pdr_sphere')
FIGDIR = os.path.join(HERE, '..', 'figures')
XC_TOT = 1.6e-4       # total gas-phase carbon per H
XI_CR0 = 2.0e-16      # unattenuated cosmic-ray ionization rate [s^-1]
RADIUS = 2.0          # [pc]
LABEL = {'ncr': 'NCR', 'gow17': 'GOW17'}

# slice panels: title, quantity, colour map, norm
SLICES = (
    (r'$n_{\rm H}$ [cm$^{-3}$]', lambda d: d['nH'], 'Greys', None),  # set from the data
    (r'$x_{\rm e}$', lambda d: d['x_e'], 'viridis', LogNorm(1e-7, 1e-1)),
    (r'$2x_{\rm H_2}$', lambda d: 2.0*d['x_h2'], 'viridis', Normalize(0, 1)),
    (r'$x_{\rm C^+}/x_{\rm C,tot}$', lambda d: d['xCII']/XC_TOT, 'viridis',
     Normalize(0, 1)),
    (r'$x_{\rm C}/x_{\rm C,tot}$', lambda d: d['xCI']/XC_TOT, 'viridis', Normalize(0, 1)),
    (r'$x_{\rm CO}/x_{\rm C,tot}$', lambda d: d['xCO']/XC_TOT, 'viridis',
     Normalize(0, 1)),
    # cloud only: the ambient gas (~3000 K) saturates
    (r'$T$ [K]', lambda d: d['temp'], 'magma', Normalize(5, 35)),
    (r'$\chi_{\rm PE}$', lambda d: d['chi_PE'], 'inferno', LogNorm(1e-3, 1)),
    (r'$\chi_{\rm LW}$', lambda d: d['chi_LW'], 'inferno', LogNorm(1e-6, 1)),
    (r'$\chi_{\rm H_2}$', lambda d: d['chi_H2'], 'inferno', LogNorm(1e-12, 1)),
    (r'$\chi_{\rm C\,I}$', lambda d: d['chi_CI'], 'inferno', LogNorm(1e-9, 1)),
    (r'$\xi_{\rm cr}/\xi_{\rm cr,0}$', lambda d: d['xi_CR']/XI_CR0, 'cividis',
     Normalize(0, 1)),
)
# profile panels: title, y scale, y limits, curves (label, quantity)
PROFILES = (
    (r'$n_{\rm H}$ [cm$^{-3}$]', 'log', None, ((None, lambda d: d['nH']),)),
    (r'$T$ [K]', 'linear', (0, 40), (('gas', lambda d: d['temp']),
                                    ('dust', lambda d: d['temp_dust']))),
    ('hydrogen and electrons', 'log', (1e-7, 2),
     ((r'$x_{\rm e}$', lambda d: d['x_e']), (r'$x_{\rm H}$', lambda d: d['x_h']),
      (r'$2x_{\rm H_2}$', lambda d: 2.0*d['x_h2']))),
    ('carbon, fraction of the total', 'log', (1e-4, 2),
     ((r'C$^+$', lambda d: d['xCII']/XC_TOT), ('C', lambda d: d['xCI']/XC_TOT),
      ('CO', lambda d: d['xCO']/XC_TOT))),
    ('radiation field [Draine]', 'log', (1e-12, 2),
     ((r'$\chi_{\rm PE}$', lambda d: d['chi_PE']), (r'$\chi_{\rm LW}$', lambda d: d['chi_LW']),
      (r'$\chi_{\rm H_2}$', lambda d: d['chi_H2']), (r'$\chi_{\rm C\,I}$', lambda d: d['chi_CI']),
      (r'$\chi_{\rm CO}$', lambda d: d['chi_CO']))),
    (r'$\xi_{\rm cr}/\xi_{\rm cr,0}$', 'linear', (0, 1.05),
     ((None, lambda d: d['xi_CR']/XI_CR0),)),
)


def suptitle_text(d, mode, nh0, radius):
    """Three lines: the cloud; the radiation; the chemistry and the dust."""
    s = ast.literal_eval(str(d['setup'])) if 'setup' in d else {}
    bands = ['LyC', 'LW 912-1108 Å', 'PE 1108-2066 Å']
    lo = 0.20664
    for e in s.get('edges_um', []):
        bands.append('%.2g-%.2g µm' % (lo, e))
        lo = e
    field = {'draine78': 'Draine (1978) UV + Mathis et al. (1983) stars',
             'mmp83': 'Mathis et al. (1983)',
             'mmp83_draine11': 'Mathis et al. (1983) with Draine (2011) dilutions'}[
                 s.get('isrf', 'draine78')]
    rad = 'isotropic %s x %g' % (field, s.get('chi0', 1))
    if s.get('chi_edge'):
        rad += r' ($\chi_{\rm LW}$ = %.2f, $\chi_{\rm PE}$ = %.2f)' % tuple(s['chi_edge'])
    rad += '; bands: %s; %d directions' % (', '.join(bands), s.get('ndir', -1))
    if ((mode == 'ncr' and s.get('cool_dust') == 'true')
            or (mode == 'gow17' and s.get('dust_cooling') == 'true')):
        heat = ('absorbed power over the transported bands'
                if mode == 'gow17' or s.get('dust_heating') == 'absorbed'
                else 'FUV bands + attenuated interstellar floor (legacy)')
        alpha = ('Hollenbach & McKee (1989)' if float(s.get('alpha_gd', -1)) <= 0
                 else s['alpha_gd'])
        dust = ('dust: T_d from %s, $\\sigma_{10}$ = %s cm$^2$ H$^{-1}$, '
                'gas-dust %s, background %s K'
                % (heat, s.get('sigma10'), alpha, s.get('temp_bg')))
    else:
        dust = 'no gas-dust coupling'
    return ('Uniform sphere of %.0e H cm$^{-3}$, radius %.2g pc, in gas of 0.05 H cm$^{-3}$; '
            'plane through the centre\n%s\n%s post-processed to convergence '
            '(%d iterations); %s'
            % (nh0, radius, rad, LABEL[mode], int(d['niter']), dust))


def plot(path):
    d = dict(np.load(path))
    mode = str(d['mode'])
    tag = os.path.basename(path)[:-4]
    nh0 = float(np.nanmax(d['nH']))
    radius = float(np.sum(d['nH'][int(np.argmin(np.abs(d['x2v']))), :] > 0.5*nh0)
                   * (d['x1v'][1] - d['x1v'][0]))/2.0
    x, y = d['x1v'], d['x2v']
    dx = x[1] - x[0]
    ext = (x[0] - dx/2, x[-1] + dx/2, y[0] - dx/2, y[-1] + dx/2)
    j = int(np.argmin(np.abs(y)))
    fig, axes = plt.subplots(3, 6, figsize=(20, 10))
    fig.subplots_adjust(left=0.05, right=0.97, top=0.85, bottom=0.06, wspace=0.45,
                        hspace=0.35)
    for n, (title, f, cmap, norm) in enumerate(SLICES):
        ax = axes[n // 6, n % 6]
        if norm is None:
            norm = LogNorm(1e-2, 3.0*nh0)
        im = ax.imshow(f(d), origin='lower', extent=ext, cmap=cmap,
                       norm=norm, interpolation='nearest')
        fig.colorbar(im, ax=ax, shrink=0.85)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r'$x_1$ [pc]', fontsize=9)
        if n % 6 == 0:
            ax.set_ylabel(r'$x_2$ [pc]', fontsize=9)
            if n == 0:
                ax.axhline(y[j], color='C1', lw=0.8, ls='--')
    for c, (title, yscale, ylim, curves) in enumerate(PROFILES):
        ax = axes[2, c]
        ax.axvspan(-radius, radius, color='0.9', zorder=0)
        for label, f in curves:
            if not (np.nan_to_num(f(d)) > 0).any():  # a quantity this mode does not have
                continue
            ax.plot(x, f(d)[j], label=label, lw=1.4)
        ax.set_yscale(yscale)
        ax.set_ylim(*(ylim if ylim is not None else (1e-2, 3.0*nh0)))
        ax.set_xlim(ext[0], ext[1])
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r'$x_1$ [pc]', fontsize=9)
        ax.grid(alpha=0.3)
        if len(curves) > 1:
            ax.legend(fontsize=8, loc='lower center', framealpha=0.8)
    fig.suptitle(suptitle_text(d, mode, nh0, radius), fontsize=12, linespacing=1.5)
    os.makedirs(FIGDIR, exist_ok=True)
    for ext_ in ('pdf', 'png'):
        out = os.path.join(FIGDIR, '%s.%s' % (tag, ext_))
        fig.savefig(out, dpi=150)
    plt.close(fig)
    print('wrote', out)


if __name__ == '__main__':
    paths = sys.argv[1:] or [os.path.join(REDUCED, 'pdr_sphere_%s.npz' % m)
                             for m in ('ncr', 'gow17')]
    for p in paths:
        if os.path.exists(p):
            plot(p)
