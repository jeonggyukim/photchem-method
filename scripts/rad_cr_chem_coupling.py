"""Schematic of how radiation, cosmic rays and the photochemistry of the gas couple.

    python rad_cr_chem_coupling.py [3d2|3d|2d] [talk|paper] [--mhd] [--cr-left]

The first argument picks the drawing of the ray tracers: both in 3D (3d2, default),
point sources in 3D only (3d), or flat (2d). The second picks the labels: symbols in
large type for slides (talk, default), or symbols with a few words each for a paper.

--mhd adds a gas-dynamics box for a coupled simulation, with the forces and heating
the other modules exert on the gas. --cr-left swaps the radiation and cosmic-ray boxes.

Writes rad_cr_chem_coupling[_3d2|_3d]_{talk,paper}[_mhd][_crleft] .png and .pdf to
../figures.
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
from clumpy_cloud import Cloud
from diffuse_rays_3d import draw_diffuse_rays

# Two point sources on different sides of the cloud, each with its adaptive ray tree
# aimed at it: (position, aim, narrow) on a 16^3 grid; the second follows one child
# at its first two splits, so its cone is a quarter as wide.
POINT_SOURCES = [((5.0, 6.0, 5.0), (0.8, 0.5, 0.35), 1),
                 ((14.0, 3.0, 4.5), (-0.47, 0.81, 0.40), 2)]

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
# --cr-left: cosmic-ray transport in the upper left and radiation in the upper right.
CRL = '--cr-left' in sys.argv
# --crpic=A|B: the picture at the top of the cosmic-ray box. Default, a cosmic ray
# gyrating along B; A, a cosmic ray scattered by the waves it drives; B, the
# self-confinement loop as a diagram.
CRPIC = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--crpic=')), '')
NAME = ('rad_cr_chem_coupling' + {'2d': '', '3d': '_3d', '3d2': '_3d2'}[MODE]
        + '_' + STYLE + ('_mhd' if SIM else '') + ('_crleft' if CRL else '')
        + ('_cr' + CRPIC if CRPIC else ''))

# Labels per style: (text, font size).
LABELS = {
    'talk': {
        'point': ('point sources', 15),
        'diffuse': ('diffuse field', 15),
        'rad_to_chem': (r'$\mathcal{E}_{\rm LyC}$, $\mathcal{E}_{\rm FUV}$, '
                        r'$\tilde{\mathcal{E}}_{\rm H_2}$, $\cdots$', 19),
        'chem_to_rad': (r'$\chi_\nu$, $N_{\rm H_2}$, $N_{\rm C}$, $N_{\rm CO}$', 19),
        'cr_to_chem': (r'$\xi_{\rm cr}\propto e_{\rm c}$, $\Gamma_{\rm cr}$', 19),
        'chem_to_cr': (r'$x_{\rm i},\ x_{\rm n},\ T\ \rightarrow\ \sigma_{\parallel},\ '
                       r'v_{\rm A,i}$', 19),
        'species': ('', 10),
        'rad_to_cr': (r'$\mathcal{E}_{\rm rad}$', 18),
        'rad_to_cr_note': ('inverse\nCompton\n(CR e$^-$)', 14),
        'rad_to_mhd': (r'$\mathbf{f}_{\rm rad}$', 19),
        'cr_to_mhd': (r'$\mathbf{G}$', 19),
        'chem_to_mhd': (r'$n\Gamma - n^2\Lambda$', 19),
        'mhd_body': (r'$\rho,\ \mathbf{v},\ \mathbf{B}$', 17),
        'mhd_note': ('star formation\ngravity\ngalactic shear\nstellar feedback', 15),
    },
    'paper': {
        'point': ('point sources: adaptive rays\nthat split as they spread', 10.5),
        'diffuse': ('diffuse field: parallel rays,\nexternal background, scattering', 10.5),
        'rad_to_chem': ('photoionization,\nphotodissociation,\nPE heating\n'
                        r'$\mathcal{E}_{\rm LyC}$, $\mathcal{E}_{\rm FUV}$, '
                        r'$\tilde{\mathcal{E}}_{\rm H_2}$, $\cdots$', 14.5),
        'chem_to_rad': ('opacity ' r'$\chi_\nu$ (dust, HI, H$_2$)' '\n'
                        r'shielding $N_{\rm H_2}$, $N_{\rm C}$, $N_{\rm CO}$', 14.5),
        'cr_to_chem': ('ionization, heating\n'
                       r'$\xi_{\rm cr}\propto e_{\rm c}$, $\Gamma_{\rm cr}$', 14.5),
        'chem_to_cr': ('ion and neutral fractions, $T$\n$\\rightarrow$ wave damping\n'
                       r'$x_{\rm i},\ x_{\rm n},\ T\ \rightarrow\ \sigma_{\parallel},\ '
                       r'v_{\rm A,i}$', 14.5),
        'rad_to_cr': (r'$\mathcal{E}_{\rm rad}$', 15.5),
        'rad_to_cr_note': ('inverse-\nCompton\nlosses\nof CR e$^-$', 11.5),
        'rad_to_mhd': ('force\n' r'$\mathbf{f}_{\rm rad}$', 14.5),
        'cr_to_mhd': ('force,\nheating\n'
                      r'$\mathbf{G}$, $\mathbf{v}_{\rm s}\!\cdot\!\mathbf{G}$', 14.5),
        'chem_to_mhd': ('net\nheating\n' r'$n\Gamma - n^2\Lambda$', 14.5),
        'mhd_body': (r'$\rho,\ \mathbf{v},\ \mathbf{B}$', 15),
        'mhd_note': ('star formation\ngravity\ngalactic shear\nstellar feedback', 14),
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


def rippled_field(ax, xa, xb, y, packets, amp=0.075, width=0.32):
    """A field line from xa to xb carrying Alfven-wave packets, given as (centre,
    wavelength) pairs; returns the line's height as a function of x."""
    def yf(x):
        return y + amp*sum(np.exp(-((x - c)/width)**2)*np.sin(2*np.pi*(x - c)/lam)
                           for c, lam in packets)
    s = np.linspace(xa, xb, 800)
    ax.plot(s, yf(s), color='0.35', lw=1.4, zorder=2)
    ax.annotate('', xy=(xb + 0.15, y), xytext=(xb - 0.05, y),
                arrowprops=dict(arrowstyle='-|>', color='0.35', lw=1.4))
    ax.text(xb + 0.22, y + 0.02, r'$\mathbf{B}$', fontsize=13, va='center', color='0.3')
    return yf


def scattered_helix(ax, xa, xb, y, xk, r0=0.3, pitch=(0.42, 0.2)):
    """A cosmic ray gyrating from xa, with its pitch (advance per gyration) changing
    from pitch[0] to pitch[1] at the wave packet at xk: pitch-angle scattering.  The
    orbit is drawn front and back, as in draw_cosmic_ray."""
    xs, ys, x, phi = [], [], xa, 0.0
    while x < xb:
        p = pitch[0] if x < xk else pitch[1]
        dphi = 2*np.pi/80
        x += p*dphi/(2*np.pi)
        phi += dphi
        xs.append(x + 0.1*np.cos(phi))
        ys.append(y + r0*np.sin(phi))
    xs, ys = np.array(xs), np.array(ys)
    front = np.cos(np.linspace(0, phi, len(xs))) > 0
    ax.plot(np.where(front, np.nan, xs), np.where(front, np.nan, ys), color=COL_CR,
            lw=1.0, alpha=0.45, zorder=1)
    ax.plot(np.where(front, xs, np.nan), np.where(front, ys, np.nan), color=COL_CR,
            lw=1.8, zorder=3)
    ax.plot(xs[-1], ys[-1], 'o', ms=6, color=COL_CR, zorder=4)


# Helix pitch (advance along B per gyration) before and after the scattering; each
# wave packet's wavelength equals the pitch of the orbit it resonates with, one
# wavelength per gyration, so the shorter wave pairs with the smaller pitch angle
# cosine (k = 1/(mu r_L)).
PITCH = (0.42, 0.2)


def resonant_waves(ax, x0, y, r0):
    """Field line, two wave packets and the scattered helix shared by A and C."""
    xk, xk2 = x0 + 1.9, x0 + 3.4
    rippled_field(ax, x0 + 0.45, x0 + 4.35, y, ((xk, PITCH[0]), (xk2, PITCH[1])))
    scattered_helix(ax, x0 + 0.5, x0 + 4.1, y, xk, r0=r0, pitch=PITCH)


def draw_cr_scatter(ax, x0, y):
    """A: a cosmic ray streaming along B, scattered in pitch angle by the resonant
    Alfven waves it drives: first a long wave, then, at smaller pitch-angle cosine,
    a shorter one."""
    resonant_waves(ax, x0, y, 0.3)


def draw_cr_loop(ax, x0, y):
    """B: the self-confinement loop.  Streaming down the CR pressure gradient drives
    waves, the waves scatter the cosmic rays, and scattering limits the streaming;
    damping by the gas removes wave energy."""
    nodes = {'grad': (x0 + 0.75, y, r'$\nabla P_{\rm c}$'),
             'wave': (x0 + 2.6, y, r'waves $\delta\mathbf{B}$'),
             'scat': (x0 + 4.45, y, r'$\sigma_\parallel$'),
             'damp': (x0 + 2.6, y - 0.5, r'damping: IN ($n_{\rm n}$), NLL ($T$)')}
    for xn, yn, s in nodes.values():
        ax.text(xn, yn, s, ha='center', va='center', fontsize=13, color='0.1', zorder=4,
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=COL_CR, lw=1.2))

    def link(a, b, text, rad=0.0, dy=0.13, col=COL_CR):
        (xa_, ya_, _), (xb_, yb_, _) = nodes[a], nodes[b]
        ax.add_patch(FancyArrowPatch((xa_, ya_), (xb_, yb_), shrinkA=26, shrinkB=26,
                                     arrowstyle='-|>,head_length=6,head_width=3.5',
                                     connectionstyle='arc3,rad=%g' % rad, color=col,
                                     lw=1.6, zorder=3))
        ax.text(0.5*(xa_ + xb_), 0.5*(ya_ + yb_) + dy, text, ha='center',
                va='bottom' if dy > 0 else 'top', fontsize=10.5, color=col)
    link('grad', 'wave', r'streaming, $v_{\rm D}>v_{\rm A,i}$', dy=-0.12)
    link('wave', 'scat', 'scattering', dy=-0.12)
    ax.add_patch(FancyArrowPatch((x0 + 4.45, y + 0.2), (x0 + 0.75, y + 0.2),
                                 shrinkA=2, shrinkB=2,
                                 arrowstyle='-|>,head_length=6,head_width=3.5',
                                 connectionstyle='arc3,rad=0.08', color='0.45', lw=1.3,
                                 ls=(0, (3, 2)), zorder=3))
    ax.text(x0 + 2.6, y + 0.35, r'limits $v_{\rm D}\rightarrow v_{\rm A,i}$', ha='center',
            va='center', fontsize=10.5, color='0.4', zorder=4,
            bbox=dict(boxstyle='square,pad=0.1', fc='#E5EFF7', ec='none'))
    ax.add_patch(FancyArrowPatch((x0 + 2.6, y - 0.36), (x0 + 2.6, y - 0.16),
                                 arrowstyle='-[,widthB=0.6,lengthB=0.2', color='#D55E00',
                                 lw=1.6, zorder=3))


def atom(ax, x, y, r, fc, label=None, lc='white', k=1.0):
    ax.add_patch(Circle((x, y), k*r, fc=fc, ec='0.25', lw=0.8, zorder=3))
    if label:
        ax.text(x, y, label, ha='center', va='center', fontsize=8*k, color=lc,
                fontweight='bold', zorder=4)


def photon(ax, x, y, label, color, k=1.0):
    """A wavy photon arriving at (x, y) from the upper left, labelled at its tail;
    k scales the drawing."""
    s = np.linspace(0, 1, 200)
    ax.plot(x + k*(-0.5 + 0.36*s), y + k*(0.22 - 0.16*s + 0.035*np.sin(2*np.pi*5*s)),
            color=color, lw=1.2*k)
    ax.annotate('', xy=(x - 0.1*k, y + 0.04*k), xytext=(x - 0.16*k, y + 0.07*k),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=1.2*k))
    ax.text(x - 0.5*k, y + 0.25*k, label, fontsize=11*k, color=color, ha='center',
            va='bottom')


def electron(ax, x, y, k=1.0):
    """An electron leaving (x, y) toward the upper right."""
    ax.annotate('', xy=(x + 0.27*k, y + 0.2*k), xytext=(x + 0.1*k, y + 0.07*k),
                arrowprops=dict(arrowstyle='-|>', color=COL_CR, lw=1.0*k))
    atom(ax, x + 0.31*k, y + 0.23*k, 0.05, COL_CR, k=k)
    ax.text(x + 0.39*k, y + 0.25*k, r'e$^-$', fontsize=9*k, va='center', color=TXT)


def draw_molecules(ax, x, y, k=1.0):
    """Top row: a Lyman-continuum photon ionizing H, and a far-ultraviolet photon
    ejecting a photoelectron from a grain.  Bottom row: H2, C+ and CO.  k scales
    the drawing about (x, y)."""
    yt, yb = y + 0.2*k, y - 0.33*k
    photon(ax, x, yt, 'LyC', '#8E44AD', k)
    atom(ax, x, yt, 0.13, 'white', r'H$^+$', '0.2', k)
    electron(ax, x, yt, k)
    photon(ax, x + 1.25*k, yt, 'FUV', '#D55E00', k)
    atom(ax, x + 1.25*k, yt, 0.15, '#A0785A', 'gr', k=k)
    electron(ax, x + 1.25*k, yt, k)

    ax.plot([x - 0.04*k, x + 0.24*k], [yb, yb], color='0.25', lw=2.5*k, zorder=2)
    atom(ax, x - 0.04*k, yb, 0.13, 'white', 'H', '0.2', k)
    atom(ax, x + 0.24*k, yb, 0.13, 'white', 'H', '0.2', k)
    atom(ax, x + 0.75*k, yb, 0.15, '0.3', r'C$^+$', k=k)
    ax.plot([x + 1.18*k, x + 1.42*k], [yb, yb], color='0.25', lw=2.5*k, zorder=2)
    atom(ax, x + 1.18*k, yb, 0.14, '0.3', 'C', k=k)
    atom(ax, x + 1.42*k, yb, 0.14, '#D55E00', 'O', k=k)


def main():
    top = 8.9   # upper edge of the upper boxes
    # With --mhd the gas-dynamics box is a third column between radiation and cosmic
    # rays, so the figure widens to about 16:9 instead of growing taller.
    xc = 3.7 if SIM else 0.0    # shift of the cosmic-ray box
    xm = 1.85 if SIM else 0.0   # shift of the chemistry box
    width = 12.0 + xc
    # The chemistry box spans cb..ct. The talk labels are one line each, so the gap
    # between it and the upper boxes is narrower there; the paper labels run to four
    # lines and keep the full gap.
    ct = 2.55 if PAPER else 3.1
    cb = 0.6 if PAPER else 1.1
    bot = cb - 0.15

    def gy(y):
        """Map a height in the gap as laid out for ct = 2.55 to the actual gap."""
        return 0.5*(4.4 + ct) + (y - 3.475)*(4.4 - ct)/1.85
    fig, ax = plt.subplots(figsize=(width, top + 0.25 - bot))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, width)
    ax.set_ylim(bot, top + 0.25)
    ax.set_aspect('equal')
    ax.axis('off')

    # --cr-left mirrors the arrows and their labels about the vertical centre line and
    # swaps the two upper boxes; the chemistry and gas boxes are centred already.
    def mx(x):
        return width - x if CRL else x

    def hflip(h):
        return {'left': 'right', 'right': 'left'}.get(h, h) if CRL else h

    orad = width - 5.75 if CRL else 0.0   # shift of the radiation box
    ocr = 0.25 - 6.5 if CRL else xc       # shift of the cosmic-ray box

    # radiation transfer
    box(ax, 0.25 + orad, 4.4, 5.5 + orad, top, COL_RAD, 'Radiation Transfer')
    if FLAT:
        draw_point_source(ax, 1.55 + orad, 6.9)
    else:
        fig.canvas.draw()
        to_fig = ax.transData + fig.transFigure.inverted()
        (fx0, fy0), (fx1, fy1) = to_fig.transform([(0.3 + orad, 5.95), (2.85 + orad, 8.6)])
        inset = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0], projection='3d')
        inset.set_facecolor('none')
        draw_ray_tree(inset, lw=0.5, zoom=1.05, ms=0.3, sources=POINT_SOURCES,
                      cloud=Cloud(16.0, (10.0, 10.0, 8.0), 3.8), kappa=0.25)
    s, fs = LABELS['point']
    ax.text(1.55 + orad, 5.95, s, ha='center', va='top', fontsize=fs, color=TXT,
            linespacing=1.25)
    if MODE == '3d2':
        (fx0, fy0), (fx1, fy1) = to_fig.transform([(2.95 + orad, 5.95),
                                                   (5.45 + orad, 8.6)])
        inset2 = fig.add_axes([fx0, fy0, fx1 - fx0, fy1 - fy0], projection='3d')
        inset2.set_facecolor('none')
        draw_diffuse_rays(inset2, lw=0.45, zoom=1.05, nbundle=4, spacing=1.6,
                          cloud=Cloud(8.0, (4.0, 4.0, 4.0), 2.6, seed=7))
    else:
        draw_diffuse(ax, 4.28 + orad, 6.9, 1.5, 1.25)
    s, fs = LABELS['diffuse']
    ax.text(4.28 + orad, 5.95, s, ha='center', va='top', fontsize=fs, color=TXT,
            linespacing=1.25)
    # the transfer equation, then the energy density (inverse-Compton losses) and the
    # flux (radiation force) that the other boxes take from it
    ax.text(2.92 + orad, 5.25, r'$\hat{\mathbf{n}}\cdot\nabla I_\nu = -\chi_\nu I_\nu'
            r' + \eta_\nu$', ha='center', va='center', fontsize=15, color='0.1')
    ax.text(2.92 + orad, 4.8, r'$\mathcal{E}_{\rm rad} = \dfrac{1}{c}\int\!\oint I_\nu'
            r'\,d\Omega\,d\nu,\quad \mathbf{F}_{\rm rad} = \int\!\oint I_\nu\,'
            r'\hat{\mathbf{n}}\,d\Omega\,d\nu$', ha='center', va='center', fontsize=14,
            color='0.1')

    # cosmic-ray transport: the two-moment equations of Armillotta et al. (2021,
    # Eqs. 1, 2 and 4), with the isotropic pressure P_c = e_c/3 that turns the bracket
    # F_c - v.(P_c + e_c I) into F_c - (4/3) v e_c; then the scattering coefficient of
    # Armillotta et al. (2022, Eqs. 16 and 17), the smaller of the nonlinear-Landau and
    # the ion-neutral values, with their dependence on n_i, n_n and T
    # (c_s^{-1/2} ~ T^{-1/4}; 1/v_A,i ~ n_i^{1/2}).
    box(ax, 6.5 + ocr, 4.4, 11.75 + ocr, top, COL_CR, 'Cosmic-ray Transport')
    if CRPIC == 'A':
        draw_cr_scatter(ax, 6.5 + ocr, 7.95)
    elif CRPIC == 'B':
        draw_cr_loop(ax, 6.5 + ocr, 7.95)
    else:
        draw_cosmic_ray(ax, 6.95 + ocr, 10.9 + ocr, 7.95)

    # radiation -> cosmic rays: the inverse-Compton losses of cosmic-ray electrons need
    # the radiation energy density (Linzer et al. 2025, App. C.2), so this link exists
    # only when the electrons are followed; dashed for that reason. With --mhd it
    # passes over the gas-dynamics box.
    y_ic = 8.15 if SIM else 6.4
    x_ic = 0.5*(5.5 + 6.5 + xc)
    ax.add_patch(FancyArrowPatch((mx(5.52), y_ic), (mx(6.48 + xc), y_ic),
                                 arrowstyle='-|>,head_length=8,head_width=4.5',
                                 color=COL_RAD, lw=2.0, ls=(0, (4, 2)), shrinkA=0,
                                 shrinkB=0))
    s, fs = LABELS['rad_to_cr']
    ax.text(x_ic, y_ic + 0.15, s, ha='center', va='bottom', fontsize=fs, color=COL_RAD)
    s, fs = LABELS['rad_to_cr_note']
    if SIM:
        s = s.replace('-\n', '-').replace('\n', ' ')
    ax.text(x_ic, y_ic - 0.15, s, ha='center', va='top', fontsize=fs, color=COL_RAD,
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
        r'\,T^{-1/4}\,n_{\rm i}^{-1/4},\ \ '
        r'\sigma_{\parallel,\rm IN} \propto |\hat{\mathbf{B}}\cdot\nabla P_{\rm c}|'
        r'\,n_{\rm i}^{-1/2}\,n_{\rm n}^{-1}$',
    ]
    # picture B reaches lower, so the equations start lower and sit closer
    y_eq, dy_eq = (6.95, 0.53) if CRPIC == 'B' else (7.2, 0.58)
    for k, s in enumerate(cr_eq):
        ax.text(9.125 + ocr, y_eq - dy_eq*k, s, ha='center', va='center', fontsize=13,
                color='0.1')

    # photochemistry
    box(ax, 2.6 + xm, cb, 9.4 + xm, ct, COL_CHEM, 'Photochemistry and Thermodynamics')
    # molecules and equations sit 1.0-1.2 below the top; the paper's species line
    # takes the bottom 0.3
    ym = ct - (1.0 if PAPER else 1.2)
    draw_molecules(ax, (3.35 if PAPER else 3.45) + xm, ym - 0.15, k=1.1 if PAPER else 1.3)
    # the two equations, and under them what every rate, heating and cooling term
    # depends on (Kim et al. 2023, Eqs. 12-14): density, temperature, abundances,
    # gas metallicity, dust abundance, radiation and cosmic rays
    ye = ym + 0.2
    ax.text(5.85 + xm, ye, r'$\dfrac{dx_i}{dt} = C_i - D_i\,x_i$', ha='left',
            va='center', fontsize=15, color='0.1')
    ax.text(7.6 + xm, ye, r'$\dfrac{de}{dt} = n\Gamma - n^2\Lambda$', ha='left',
            va='center', fontsize=15, color='0.1')
    ax.text(7.55 + xm, ye - 0.6, r'$C_i,\ D_i,\ \Gamma,\ \Lambda\,(n,\,T,\,x_s,\,'
            r"Z'_{\rm g},\,Z'_{\rm d},\,\mathcal{E},\,\xi_{\rm cr})$", ha='center',
            va='center', fontsize=14, color='0.1')
    s, fs = LABELS['species']
    if s:
        ax.text(6.0 + xm, cb + 0.15, s, ha='center', va='center', fontsize=fs, color=TXT)

    # radiation <-> chemistry: each label sits on the side of its own arrow, the
    # downward one inside and the upward one outside, at different heights so that
    # neither runs into the other.
    xr = 5.3 if SIM else 4.65   # radiation -> chemistry; chemistry -> radiation 0.35 left
    arrow(ax, (mx(xr), 4.38), (mx(xr), ct + 0.02), COL_RAD)
    arrow(ax, (mx(xr - 0.35), ct + 0.02), (mx(xr - 0.35), 4.38), COL_CHEM)
    s, fs = LABELS['rad_to_chem']
    ax.text(mx(xr + 0.2), gy(3.78 if PAPER else 3.85), s, ha=hflip('left'),
            va='center', fontsize=fs, color=COL_RAD, linespacing=1.25)
    s, fs = LABELS['chem_to_rad']
    ax.text(mx(xr - 0.55), gy(3.5), s, ha=hflip('right'), va='center', fontsize=fs,
            color=COL_CHEM, linespacing=1.25)

    # cosmic rays <-> chemistry; with --mhd the cosmic-ray label goes higher, above
    # the label of the chemistry -> gas arrow.
    xq = 10.4 if SIM else 7.35   # cosmic rays -> chemistry; chemistry -> CR 0.35 right
    arrow(ax, (mx(xq), 4.38), (mx(xq), ct + 0.02), COL_CR)
    arrow(ax, (mx(xq + 0.35), ct + 0.02), (mx(xq + 0.35), 4.38), COL_CHEM)
    s, fs = LABELS['cr_to_chem']
    ax.text(mx(xq - 0.2), gy(3.85 if SIM else (2.88 if PAPER else 3.0)), s,
            ha=hflip('right'), va='center', fontsize=fs, color=COL_CR, linespacing=1.25)
    s, fs = LABELS['chem_to_cr']
    ax.text(mx(xq + 0.55), gy(3.5), s, ha=hflip('left'), va='center', fontsize=fs,
            color=COL_CHEM, linespacing=1.25)

    # gas dynamics: radiation force, cosmic-ray force and heating, and the net heating
    # from the chemistry act on the gas, which moves, compresses and carries B.
    if SIM:
        x0, x1, y0, y1 = 6.65, 9.05, 4.75, 7.15
        xg = 0.5*(x0 + x1)
        box(ax, x0, y0, x1, y1, COL_MHD, 'Gas Dynamics')
        s, fs = LABELS['mhd_body']
        ax.text(xg, 6.45, s, ha='center', va='center', fontsize=fs, color='0.1')
        s, fs = LABELS['mhd_note']
        ax.text(xg, 5.45, s, ha='center', va='center', fontsize=fs, color=TXT,
                linespacing=1.3)
        y_f = 5.75
        arrow(ax, (mx(5.52), y_f), (mx(x0 - 0.02), y_f), COL_RAD)
        arrow(ax, (mx(6.48 + xc), y_f), (mx(x1 + 0.02), y_f), COL_CR)
        arrow(ax, (xg, ct + 0.02), (xg, y0 - 0.02), COL_CHEM)
        s, fs = LABELS['rad_to_mhd']
        ax.text(mx(0.5*(5.5 + x0)), y_f + 0.12, s, ha='center', va='bottom',
                fontsize=fs, color=COL_RAD, linespacing=1.15)
        s, fs = LABELS['cr_to_mhd']
        ax.text(mx(0.5*(x1 + 6.5 + xc)), y_f + 0.12, s, ha='center', va='bottom',
                fontsize=fs, color=COL_CR, linespacing=1.15)
        s, fs = LABELS['chem_to_mhd']
        ax.text(mx(xg + 0.15), gy(3.0), s, ha=hflip('left'), va='center',
                fontsize=fs,
                color=COL_CHEM, linespacing=1.15)

    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.normpath(os.path.join(OUT, NAME + '.' + ext))
        fig.savefig(path, dpi=200)
        print('wrote', path)


if __name__ == '__main__':
    main()
