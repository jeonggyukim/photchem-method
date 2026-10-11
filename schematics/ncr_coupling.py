"""Schematic of the TIGRESS-NCR coupling, for a talk.  Gas dynamics is on top.
Radiation transfer by adaptive ray tracing from point sources (no diffuse field) is at
lower left.  Photochemistry and thermodynamics is at lower right, with a cosmic ray
ionizing H2.  Cosmic rays enter only as a prescribed ionization rate; there is no
transport module.

    python ncr_coupling.py [--no-cr-icon]

--no-cr-icon draws H2 without the cosmic ray that ionizes it (ncr_coupling_nocricon).

Writes ncr_coupling.png and .pdf to ../figures. Drawing helpers (boxes, arrows,
molecules, the adaptive ray tree) come from ncr_cr_coupling.py beside this script.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    allow_abbrev=False)
_parser.add_argument('--no-cr-icon', action='store_true',
                     help='draw H2 without the cosmic ray that ionizes it')
CR_ICON = not _parser.parse_args().no_cr_icon
sys.argv = sys.argv[:1]   # ncr_cr_coupling reads its style from argv at import
from ncr_cr_coupling import (COL_CHEM, COL_MHD, COL_RAD, POINT_SOURCES, TXT,  # noqa: E402
                             arrow, atom, box, electron, make_cloud, photon)
from rays_adaptive import draw_ray_tree  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
COL_CRION = '#0072B2'
W, H = 14.0, 9.0


def draw_chem(ax, x, y, k=1.0):
    """Top row: a LyC photon ionizing H, an FUV photon ejecting a photoelectron from a
    grain. Bottom row: a cosmic ray ionizing H2 (an electron leaves), then C+ and CO.
    k scales the drawing about (x, y)."""
    yt, yb = y + 0.2*k, y - 0.4*k
    photon(ax, x, yt, 'LyC', '#8E44AD', k)
    atom(ax, x, yt, 0.13, 'white', r'H$^+$', '0.2', k)
    electron(ax, x, yt, k)
    photon(ax, x + 1.25*k, yt, 'FUV', '#D55E00', k)
    atom(ax, x + 1.25*k, yt, 0.15, '#A0785A', 'gr', k=k)
    electron(ax, x + 1.25*k, yt, k)
    # the cosmic ray passes through H2 and survives: H2 + CR -> H2+ + e- + CR
    if CR_ICON:
        ax.plot([x - 0.5*k, x + 0.55*k], [yb + 0.07*k, yb - 0.04*k], color=COL_CRION,
                lw=2.0*k, zorder=1)
        atom(ax, x + 0.63*k, yb - 0.05*k, 0.085, COL_CRION, k=k)
        ax.text(x + 0.63*k, yb - 0.05*k, 'CR', ha='center', va='center',
                fontsize=6*k, color='white', fontweight='bold', zorder=4)
    ax.plot([x - 0.04*k, x + 0.24*k], [yb, yb], color='0.25', lw=2.5*k, zorder=2)
    atom(ax, x - 0.04*k, yb, 0.13, 'white', 'H', '0.2', k)
    atom(ax, x + 0.24*k, yb, 0.13, 'white', 'H', '0.2', k)
    if CR_ICON:
        ax.text(x + 0.1*k, yb + 0.17*k, '+', ha='center', va='center',
                fontsize=10*k, color=TXT, fontweight='bold')
        electron(ax, x + 0.24*k, yb, k)
    atom(ax, x + 1.1*k, yb, 0.15, '0.3', r'C$^+$', k=k)
    ax.plot([x + 1.5*k, x + 1.74*k], [yb, yb], color='0.25', lw=2.5*k, zorder=2)
    atom(ax, x + 1.5*k, yb, 0.14, '0.3', 'C', k=k)
    atom(ax, x + 1.74*k, yb, 0.14, '#D55E00', 'O', k=k)


def main():
    fig, ax = plt.subplots(figsize=(W, H))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect('equal')
    ax.axis('off')

    # gas dynamics, top centre
    gx0, gx1, gy0, gy1 = 3.6, 10.4, 6.55, 8.85
    gxc = 0.5*(gx0 + gx1)
    box(ax, gx0, gy0, gx1, gy1, COL_MHD, '')
    ax.text(gxc, gy1 - 0.12, 'Gas Dynamics', ha='center', va='top', fontsize=23,
            fontweight='bold', color='0.1', linespacing=1.1)
    ax.text(gxc, 7.85, r'$\rho,\ \mathbf{v},\ \mathbf{B}$', ha='center', va='center',
            fontsize=32, color='0.1')
    ax.text(gxc, 7.08, 'gravity · turbulence · galactic shear\n'
            'star formation · supernovae', ha='center', va='center', fontsize=22,
            color=TXT, linespacing=1.3)

    # radiation transfer, lower left: adaptive rays from point sources only
    rx0, rx1, ry0, ry1 = 0.25, 6.25, 0.3, 5.1
    box(ax, rx0, ry0, rx1, ry1, COL_RAD, '')
    ax.text(0.5*(rx0 + rx1), ry1 - 0.12, 'Radiation Transfer', ha='center', va='top', fontsize=23,
            fontweight='bold', color='0.1', linespacing=1.1)
    fig.canvas.draw()
    to_fig = ax.transData + fig.transFigure.inverted()
    (fx0, fy0), (fx1, fy1) = to_fig.transform([(0.3, 1.05), (3.55, 4.75)])
    inset = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0], projection='3d')
    inset.set_facecolor('none')
    draw_ray_tree(inset, lw=0.6, zoom=1.05, ms=0.35, sources=POINT_SOURCES,
                  cloud=make_cloud(16.0), kappa=0.25)
    ax.text(1.95, 1.0, 'adaptive rays from\nstar clusters',
            ha='center', va='top', fontsize=18, color=TXT, linespacing=1.25)
    ax.text(4.9, 3.55, r'$\mathbf{\hat{n}}\cdot\nabla I_\nu = -\chi_\nu I_\nu$',
            ha='center', va='center', fontsize=26, color='0.1')
    ax.text(4.9, 2.6, r'$\mathcal{E}_\nu = \dfrac{1}{c}\oint I_\nu\,d\Omega$',
            ha='center', va='center', fontsize=26, color='0.1')
    ax.text(4.9, 1.65, r'$\mathbf{F}_\nu = \oint I_\nu\,\mathbf{\hat{n}}\,d\Omega$',
            ha='center', va='center', fontsize=26, color='0.1')

    # photochemistry and thermodynamics, lower right
    cx0, cx1, cy0, cy1 = 7.75, 13.75, 0.3, 5.1
    box(ax, cx0, cy0, cx1, cy1, COL_CHEM, '')
    ax.text(0.5*(cx0 + cx1), cy1 - 0.12, 'Photochemistry and\nThermodynamics', ha='center', va='top', fontsize=23,
            fontweight='bold', color='0.1', linespacing=1.1)
    draw_chem(ax, 9.85, 3.3, k=1.55)
    ax.text(9.5, 1.75, r'$\dfrac{dx_i}{dt} = C_i - D_i\,x_i$', ha='center', va='center',
            fontsize=21, color='0.1')
    ax.text(12.3, 1.75, r'$\dfrac{de}{dt} = n\Gamma - n^2\Lambda$', ha='center',
            va='center', fontsize=21, color='0.1')
    ax.text(10.75, 0.82, 'CR ionization ' r'$\xi_{\rm cr}$' ' prescribed:\n'
            r'$\propto\Sigma_{\rm SFR}/\Sigma_{\rm gas}$, attenuated by $N_{\rm eff}$',
            ha='center', va='center', fontsize=20, color=COL_CRION, linespacing=1.3)

    # radiation <-> chemistry, across the gap between the lower boxes
    arrow(ax, (rx1 + 0.03, 3.3), (cx0 - 0.03, 3.3), COL_RAD)
    arrow(ax, (cx0 - 0.03, 2.1), (rx1 + 0.03, 2.1), COL_CHEM)
    xm = 0.5*(rx1 + cx0)
    ax.text(xm, 3.45, r'$\mathcal{E}_{\rm LyC}$,' '\n' r'$\mathcal{E}_{\rm FUV}$,' '\n'
            r'$\tilde{\mathcal{E}}_{\rm H_2}$,' '\n' r'$\tilde{\mathcal{E}}_{\rm C}$',
            ha='center', va='bottom', fontsize=21,
            color=COL_RAD, linespacing=1.2)
    ax.text(xm, 1.95, r'$\chi_\nu$,' '\n' r'$N_{\rm H_2}$, $N_{\rm C}$,' '\n' r'$N_{\rm CO}$',
            ha='center', va='top', fontsize=21, color=COL_CHEM, linespacing=1.2)

    # the lower boxes -> gas dynamics: radiation force and net heating (the gas state the
    # other modules read is left implicit)
    for xu, col, lab, ha in ((4.15, COL_RAD, r'$\mathbf{f}_{\rm rad}$', 'left'),
                             (9.85, COL_CHEM, r'$n\Gamma - n^2\Lambda$', 'right')):
        arrow(ax, (xu, ry1 + 0.03), (xu, gy0 - 0.03), col)
        ax.text(xu + (0.15 if ha == 'left' else -0.15), 5.85, lab, ha=ha, va='center',
                fontsize=24, color=col)

    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.normpath(os.path.join(
            OUT, 'ncr_coupling' + ('' if CR_ICON else '_nocricon') + '.' + ext))
        fig.savefig(path, dpi=200)
        print('wrote', path)


if __name__ == '__main__':
    main()
