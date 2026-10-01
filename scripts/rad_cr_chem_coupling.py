"""Schematic of how radiation, cosmic rays and the photochemistry of the gas couple.

    python rad_cr_chem_coupling.py [3d2|3d|2d] [talk|paper] [--mhd]

The first argument picks the drawing of the ray tracers: both in 3D (3d2, default),
point sources in 3D only (3d), or flat (2d). The second picks the labels: symbols in
large type for slides (talk, default), or symbols with a few words each for a paper.

--mhd adds a gas-dynamics box for a coupled simulation, with the forces and heating
the other modules exert on the gas.

Writes rad_cr_chem_coupling[_3d2|_3d]_{talk,paper}[_mhd] .png and .pdf to ../figures.
The 3D ray panels come from adaptive_rays_3d.py and diffuse_rays_3d.py beside this
script.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

from adaptive_rays_3d import draw_ray_tree
from diffuse_rays_3d import draw_diffuse_rays

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
_ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
MODE = _ARGS[0] if len(_ARGS) > 0 else '3d2'
STYLE = _ARGS[1] if len(_ARGS) > 1 else 'talk'
FLAT = MODE == '2d'
PAPER = STYLE == 'paper'
# --mhd: a coupled simulation rather than post-processing; a gas-dynamics box between
# the others takes the radiation force, the cosmic-ray force and heating, and the
# net heating of the gas, and supplies rho, v, B to all three.
SIM = '--mhd' in sys.argv
NAME = ('rad_cr_chem_coupling' + {'2d': '', '3d': '_3d', '3d2': '_3d2'}[MODE]
        + '_' + STYLE + ('_mhd' if SIM else ''))

# Labels per style: (text, font size).
LABELS = {
    'talk': {
        'point': ('point sources', 15),
        'diffuse': ('diffuse field', 15),
        'rad_to_chem': (r'$J_{\rm LyC}$, $\tilde{J}_{\rm LW}$, $J_{\rm PE}$', 19),
        'chem_to_rad': (r'$\chi_\nu$, $N_{\rm H_2}$, $N_{\rm C}$, $N_{\rm CO}$', 19),
        'cr_to_chem': (r'$\xi_{\rm cr}\propto e_{\rm c}$, $\Gamma_{\rm cr}$', 19),
        'chem_to_cr': (r'$x_{\rm i}\ \rightarrow\ \sigma_{\parallel},\ v_{\rm A,i}$', 19),
        'species': ('', 10),
        'rad_to_cr': (r'$\mathcal{E}_{\rm rad}$', 18),
        'rad_to_cr_note': ('inverse\nCompton\n(CR e$^-$)', 14),
        'rad_to_mhd': (r'$\mathbf{f}_{\rm rad}$', 19),
        'cr_to_mhd': (r'$\mathbf{G}$', 19),
        'chem_to_mhd': (r'$n\Gamma - n^2\Lambda$', 19),
        'mhd_body': (r'$\rho,\ \mathbf{v},\ \mathbf{B}$', 17),
    },
    'paper': {
        'point': ('point sources: adaptive rays\nthat split as they spread', 10.5),
        'diffuse': ('diffuse field: parallel rays,\nexternal background, scattering', 10.5),
        'rad_to_chem': ('photoionization,\nphotodissociation,\nPE heating\n'
                        r'$J_{\rm LyC}$, $\tilde{J}_{\rm LW}$, $J_{\rm PE}$', 14.5),
        'chem_to_rad': ('opacity ' r'$\chi_\nu$ (dust, H I, H$_2$)' '\n'
                        r'shielding $N_{\rm H_2}$, $N_{\rm C}$, $N_{\rm CO}$', 14.5),
        'cr_to_chem': ('ionization, heating\n'
                       r'$\xi_{\rm cr}\propto e_{\rm c}$, $\Gamma_{\rm cr}$', 14.5),
        'chem_to_cr': ('ion fraction $\\rightarrow$\nion–neutral damping\n'
                       r'$x_{\rm i}\ \rightarrow\ \sigma_{\parallel},\ v_{\rm A,i}$', 14.5),
        'rad_to_cr': (r'$\mathcal{E}_{\rm rad}$', 15.5),
        'rad_to_cr_note': ('inverse-\nCompton\nlosses\nof CR e$^-$', 11.5),
        'rad_to_mhd': ('force\n' r'$\mathbf{f}_{\rm rad}$', 14.5),
        'cr_to_mhd': ('force,\nheating\n'
                      r'$\mathbf{G}$, $\mathbf{v}_{\rm s}\!\cdot\!\mathbf{G}$', 14.5),
        'chem_to_mhd': ('net\nheating\n' r'$n\Gamma - n^2\Lambda$', 14.5),
        'mhd_body': (r'$\rho,\ \mathbf{v},\ \mathbf{B}$', 15),
        'species': ('NCR or GOW17 network (H, C, O, ... and O, S, N ion ladders); '
                    r'gas temperature $T$', 10),
    },
}[STYLE]

COL_RAD = '#E69F00'
COL_CR = '#0072B2'
COL_CHEM = '#009E73'
COL_MHD = '0.35'
TXT = '0.2'


def box(ax, x0, y0, x1, y1, color, title):
    for fc, ec, alpha in ((color, 'none', 0.10), ('none', color, 1.0)):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle='round,pad=0.02,rounding_size=0.12',
                                    fc=fc, ec=ec, lw=1.8, alpha=alpha))
    ax.text(0.5*(x0 + x1), y1 - 0.12, title, ha='center', va='top', fontsize=16,
            fontweight='bold', color='0.1')


def arrow(ax, p, q, color):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>,head_length=8,head_width=4.5',
                                 color=color, lw=2.0, shrinkA=0, shrinkB=0))


def draw_point_source(ax, x, y):
    """A star whose photon packets split in two as they spread."""
    rng_split, rng_end = 0.42, 0.78
    for a in np.deg2rad(np.arange(0, 360, 45) + 22.5):
        xs, ys = x + rng_split*np.cos(a), y + rng_split*np.sin(a)
        ax.plot([x, xs], [y, ys], color=COL_RAD, lw=1.3, zorder=2)
        for da in (-0.2, 0.2):
            ax.plot([xs, x + rng_end*np.cos(a + da)], [ys, y + rng_end*np.sin(a + da)],
                    color=COL_RAD, lw=0.9, alpha=0.8, zorder=2)
        ax.plot(xs, ys, 'o', ms=2.5, color=COL_RAD, zorder=3)
    ax.plot(x, y, marker='*', ms=24, color='#FFD23F', mec='#C77800', mew=1.0, zorder=4)


def draw_diffuse(ax, x, y, w, h):
    """Parallel rays in three fixed directions crossing a cloud."""
    clip = Rectangle((x - w/2, y - h/2), w, h, fc='none', ec='0.6', lw=0.8, zorder=1)
    ax.add_patch(clip)
    ax.add_patch(Ellipse((x + 0.05, y - 0.03), 0.8*w, 0.55*h, angle=20, fc='0.72',
                         ec='none', alpha=0.8, zorder=1))
    ax.add_patch(Ellipse((x - 0.12, y + 0.05), 0.42*w, 0.3*h, angle=-15, fc='0.55',
                         ec='none', alpha=0.8, zorder=1))
    for ang, col in ((0.0, COL_RAD), (60.0, '#CC79A7'), (120.0, '#56B4E9')):
        a = np.deg2rad(ang)
        d = np.array([np.cos(a), np.sin(a)])
        nrm = np.array([-d[1], d[0]])
        for s in np.linspace(-1.0, 1.0, 7):
            p = np.array([x, y]) + s*0.75*max(w, h)*nrm
            line, = ax.plot([p[0] - 2*d[0], p[0] + 2*d[0]], [p[1] - 2*d[1], p[1] + 2*d[1]],
                            color=col, lw=0.8, alpha=0.9, zorder=2)
            line.set_clip_path(clip)


def draw_cosmic_ray(ax, x0, x1, y):
    """A cosmic ray gyrating along a magnetic field line."""
    s = np.linspace(x0, x1, 400)
    ax.plot(s, y + 0.06*np.sin(2*np.pi*(s - x0)/(x1 - x0)), color='0.35', lw=1.4,
            zorder=2)
    ax.annotate('', xy=(x1 + 0.15, y), xytext=(x1 - 0.05, y),
                arrowprops=dict(arrowstyle='-|>', color='0.35', lw=1.4))
    ax.text(x1 + 0.22, y + 0.02, r'$\mathbf{B}$', fontsize=13, va='center', color='0.3')
    t = np.linspace(0, 1, 800)
    phi = 2*np.pi*7*t
    xs = x0 + 0.25 + (x1 - x0 - 0.7)*t + 0.13*np.cos(phi)
    ys = y + 0.06*np.sin(2*np.pi*(xs - x0)/(x1 - x0)) + 0.32*np.sin(phi)
    front = np.cos(phi) > 0
    ax.plot(np.where(front, np.nan, xs), np.where(front, np.nan, ys), color=COL_CR,
            lw=1.0, alpha=0.45, zorder=1)
    ax.plot(np.where(front, xs, np.nan), np.where(front, ys, np.nan), color=COL_CR,
            lw=1.8, zorder=3)
    ax.plot(xs[-1], ys[-1], 'o', ms=6, color=COL_CR, zorder=4)


def atom(ax, x, y, r, fc, label=None, lc='white'):
    ax.add_patch(Circle((x, y), r, fc=fc, ec='0.25', lw=0.8, zorder=3))
    if label:
        ax.text(x, y, label, ha='center', va='center', fontsize=8, color=lc,
                fontweight='bold', zorder=4)


def draw_molecules(ax, x, y):
    """H2 being dissociated by a far-ultraviolet photon, C+, CO and an electron."""
    ax.plot([x - 0.14, x + 0.14], [y, y], color='0.25', lw=2.5, zorder=2)
    atom(ax, x - 0.14, y, 0.13, 'white', 'H', '0.2')
    atom(ax, x + 0.14, y, 0.13, 'white', 'H', '0.2')
    ax.text(x, y - 0.3, r'H$_2$', ha='center', va='top', fontsize=9.5, color=TXT)
    s = np.linspace(0, 1, 200)
    ax.plot(x - 0.75 + 0.5*s, y + 0.45 - 0.3*s + 0.04*np.sin(2*np.pi*6*s), color='#8E44AD',
            lw=1.2)
    ax.annotate('', xy=(x - 0.2, y + 0.12), xytext=(x - 0.27, y + 0.16),
                arrowprops=dict(arrowstyle='-|>', color='#8E44AD', lw=1.2))
    ax.text(x - 0.78, y + 0.5, r'$h\nu$', fontsize=9.5, color='#8E44AD', va='bottom')

    atom(ax, x + 0.8, y, 0.15, '0.3', r'C$^+$')
    ax.text(x + 0.8, y - 0.3, r'C$^+$', ha='center', va='top', fontsize=9.5, color=TXT)
    ax.plot([x + 1.38, x + 1.62], [y, y], color='0.25', lw=2.5, zorder=2)
    atom(ax, x + 1.38, y, 0.14, '0.3', 'C')
    atom(ax, x + 1.62, y, 0.14, '#D55E00', 'O')
    ax.text(x + 1.5, y - 0.3, 'CO', ha='center', va='top', fontsize=9.5, color=TXT)
    atom(ax, x + 1.12, y + 0.5, 0.06, COL_CR)
    ax.text(x + 1.22, y + 0.5, r'e$^-$', fontsize=9.5, va='center', color=TXT)


def main():
    top = 8.9   # upper edge of the two upper boxes
    cy = -1.8 if SIM else 0.0   # the chemistry box moves down to make room for MHD
    bot = 0.6 + cy
    fig, ax = plt.subplots(figsize=(12.0, 12.0*(top + 0.25 - bot)/12.0))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, 12.0)
    ax.set_ylim(bot, top + 0.25)
    ax.set_aspect('equal')
    ax.axis('off')

    # radiation transfer
    box(ax, 0.25, 4.4, 5.5, top, COL_RAD, 'Radiation transfer')
    if FLAT:
        draw_point_source(ax, 1.55, 6.9)
    else:
        fig.canvas.draw()
        to_fig = ax.transData + fig.transFigure.inverted()
        (fx0, fy0), (fx1, fy1) = to_fig.transform([(0.3, 5.75), (2.85, 8.45)])
        inset = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0], projection='3d')
        inset.set_facecolor('none')
        draw_ray_tree(inset, lw=0.5, zoom=1.05, ms=0.8)
    s, fs = LABELS['point']
    ax.text(1.55, 5.75, s, ha='center', va='top', fontsize=fs, color=TXT, linespacing=1.25)
    if MODE == '3d2':
        (fx0, fy0), (fx1, fy1) = to_fig.transform([(2.95, 5.75), (5.45, 8.45)])
        inset2 = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0], projection='3d')
        inset2.set_facecolor('none')
        draw_diffuse_rays(inset2, lw=0.45, zoom=1.05)
    else:
        draw_diffuse(ax, 4.28, 6.9, 1.5, 1.25)
    s, fs = LABELS['diffuse']
    ax.text(4.28, 5.75, s, ha='center', va='top', fontsize=fs, color=TXT, linespacing=1.25)
    ax.text(2.92, 4.62, r'$\hat{\mathbf{n}}\cdot\nabla I_\nu = -\chi_\nu I_\nu + \eta_\nu$'
            r'$,\qquad J_\nu = \frac{1}{4\pi}\oint I_\nu\,d\Omega$', ha='center',
            va='bottom', fontsize=15, color='0.1')

    # cosmic-ray transport: the two-moment equations of Armillotta et al. (2021,
    # Eqs. 1, 2 and 4), with the isotropic pressure P_c = e_c/3 that turns the bracket
    # F_c - v.(P_c + e_c I) into F_c - (4/3) v e_c; then the scattering coefficient of
    # Armillotta et al. (2022, Eqs. 16 and 17), the smaller of the nonlinear-Landau and
    # the ion-neutral values, with their dependence on n_i, n_n and T
    # (c_s^{-1/2} ~ T^{-1/4}; 1/v_A,i ~ n_i^{1/2}).
    box(ax, 6.5, 4.4, 11.75, top, COL_CR, 'Cosmic-ray transport')
    draw_cosmic_ray(ax, 6.95, 10.9, 7.95)

    # radiation -> cosmic rays: the inverse-Compton losses of cosmic-ray electrons need
    # the radiation energy density (Linzer et al. 2025, App. C.2), so this link exists
    # only when the electrons are followed; dashed for that reason.
    ax.add_patch(FancyArrowPatch((5.52, 6.4), (6.48, 6.4),
                                 arrowstyle='-|>,head_length=8,head_width=4.5',
                                 color=COL_RAD, lw=2.0, ls=(0, (4, 2)), shrinkA=0,
                                 shrinkB=0))
    s, fs = LABELS['rad_to_cr']
    ax.text(6.0, 6.55, s, ha='center', va='bottom', fontsize=fs, color=COL_RAD)
    s, fs = LABELS['rad_to_cr_note']
    ax.text(6.0, 6.25, s, ha='center', va='top', fontsize=fs, color=COL_RAD,
            linespacing=1.15)
    cr_eq = [
        r'$\dfrac{\partial e_{\rm c}}{\partial t} + \nabla\cdot\mathbf{F}_{\rm c} = '
        r'-(\mathbf{v} + \mathbf{v}_{\rm s})\cdot\mathbf{G}'
        r' + Q - \Lambda_{\rm coll}n_{\rm H}e_{\rm c}$',
        r'$\dfrac{1}{v_{\rm m}^2}\dfrac{\partial \mathbf{F}_{\rm c}}{\partial t} + '
        r'\nabla P_{\rm c} = -\mathbf{G} - \Lambda_{\rm coll}n_{\rm H}'
        r'\dfrac{\mathbf{F}_{\rm c}}{v_{\rm p}^2}$',
        r'$\mathbf{G} \equiv \sigma_{\rm tot}\cdot'
        r'[\mathbf{F}_{\rm c} - \frac{4}{3}\mathbf{v}e_{\rm c}]$',
        r'$\sigma_{\rm tot}^{-1} = \sigma_{\parallel}^{-1} + '
        r'\dfrac{v_{\rm A,i}\,(P_{\rm c} + e_{\rm c})}'
        r'{|\hat{\mathbf{B}}\cdot\nabla P_{\rm c}|}$',
        r'$\sigma_{\parallel,\rm NLL} \propto |\hat{\mathbf{B}}\cdot\nabla P_{\rm c}|^{1/2}'
        r'\,T^{-1/4}\,n_{\rm i}^{-1/4}$',
        r'$\sigma_{\parallel,\rm IN} \propto |\hat{\mathbf{B}}\cdot\nabla P_{\rm c}|'
        r'\,n_{\rm i}^{-1/2}\,n_{\rm n}^{-1}$',
    ]
    for k, s in enumerate(cr_eq):
        ax.text(9.125, 7.2 - 0.52*k, s, ha='center', va='center', fontsize=13,
                color='0.1')

    # photochemistry
    box(ax, 2.9, 0.75 + cy, 9.1, 2.55 + cy, COL_CHEM, 'Photochemistry and thermal balance')
    dy = (0.17 if PAPER else 0.0) + cy   # room for the species line below
    draw_molecules(ax, 3.75, 1.45 + dy)
    ax.text(5.6, 1.6 + dy, r'$\dfrac{dx_i}{dt} = C_i - D_i\,x_i$', ha='left',
            va='center', fontsize=15, color='0.1')
    ax.text(7.35, 1.6 + dy, r'$\dfrac{de}{dt} = n\Gamma - n^2\Lambda$', ha='left',
            va='center', fontsize=15, color='0.1')
    s, fs = LABELS['species']
    if s:
        ax.text(6.0, 0.95 + cy, s, ha='center', va='center', fontsize=fs, color=TXT)

    # radiation <-> chemistry; the two labels sit at different heights so that the
    # one beside each inner arrow does not run into the other. With the MHD box
    # between the inner arrows, both labels go to the outer side.
    arrow(ax, (4.65, 4.38), (4.65, 2.57 + cy), COL_RAD)
    arrow(ax, (4.3, 2.57 + cy), (4.3, 4.38), COL_CHEM)
    s, fs = LABELS['rad_to_chem']
    if SIM:
        ax.text(4.1, 3.6, s, ha='right', va='center', fontsize=fs, color=COL_RAD,
                linespacing=1.25)
    else:
        ax.text(4.85, 3.78 if PAPER else 3.85, s, ha='left', va='center', fontsize=fs,
                color=COL_RAD, linespacing=1.25)
    s, fs = LABELS['chem_to_rad']
    ax.text(4.1, 1.7 if SIM else 3.5, s, ha='right', va='center', fontsize=fs,
            color=COL_CHEM, linespacing=1.25)

    # cosmic rays <-> chemistry
    arrow(ax, (7.35, 4.38), (7.35, 2.57 + cy), COL_CR)
    arrow(ax, (7.7, 2.57 + cy), (7.7, 4.38), COL_CHEM)
    s, fs = LABELS['cr_to_chem']
    if SIM:
        ax.text(7.9, 3.6, s, ha='left', va='center', fontsize=fs, color=COL_CR,
                linespacing=1.25)
    else:
        ax.text(7.15, 2.88 if PAPER else 3.0, s, ha='right', va='center', fontsize=fs,
                color=COL_CR, linespacing=1.25)
    s, fs = LABELS['chem_to_cr']
    ax.text(7.9, 1.7 if SIM else 3.5, s, ha='left', va='center', fontsize=fs,
            color=COL_CHEM, linespacing=1.25)

    # gas dynamics: radiation force, cosmic-ray force and heating, and the net heating
    # from the chemistry act on the gas; rho, v, B go back to every module.
    if SIM:
        box(ax, 4.9, 1.75, 7.1, 3.25, COL_MHD, 'MHD')
        s, fs = LABELS['mhd_body']
        ax.text(6.0, 2.35, s, ha='center', va='center', fontsize=fs, color='0.1',
                linespacing=1.2)
        arrow(ax, (5.2, 4.38), (5.2, 3.27), COL_RAD)
        arrow(ax, (6.8, 4.38), (6.8, 3.27), COL_CR)
        arrow(ax, (6.0, 2.57 + cy), (6.0, 1.73), COL_CHEM)
        s, fs = LABELS['rad_to_mhd']
        ax.text(5.32, 3.85, s, ha='left', va='center', fontsize=fs, color=COL_RAD,
                linespacing=1.15)
        s, fs = LABELS['cr_to_mhd']
        ax.text(6.68, 3.85, s, ha='right', va='center', fontsize=fs, color=COL_CR,
                linespacing=1.15)
        s, fs = LABELS['chem_to_mhd']
        ax.text(6.12, 1.25, s, ha='left', va='center', fontsize=fs, color=COL_CHEM,
                linespacing=1.15)

    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.normpath(os.path.join(OUT, NAME + '.' + ext))
        fig.savefig(path, dpi=200)
        print('wrote', path)


if __name__ == '__main__':
    main()
