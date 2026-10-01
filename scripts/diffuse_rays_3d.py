"""The diffuse field solver drawn in 3D: for each of a few fixed directions of the angle
set, a bundle of parallel rays crosses the grid and a cloud in it; each ray fades with
the optical depth it has crossed.  The solver launches one ray per face cell through
every entry face; the bundle is a sample of those.

    python diffuse_rays_3d.py

Writes diffuse_rays_3d.png and .pdf to
~/Dropbox/Research/Rayt-Method/figures/photchem-postproc and copies them to
~/Documents/photchem-postproc.  draw_diffuse_rays() is also used by
rad_cr_chem_coupling.py.
"""
import os
import shutil

import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

OUT = os.path.expanduser('~/Dropbox/Research/Rayt-Method/figures/photchem-postproc')
COPY = os.path.expanduser('~/Documents/photchem-postproc')

COLORS = ('#E69F00', '#CC79A7', '#56B4E9')


def ray_paths(d, n, nbundle=3, spacing=2.0):
    """Entry and exit points of a bundle of nbundle^2 parallel rays of direction d,
    spacing cells apart, aimed through the centre of an n^3 grid of unit cells: one
    sample of the rays the solver launches, one per face cell, through every entry
    face."""
    centre = np.full(3, 0.5*n)
    e1 = np.cross(d, [0.0, 0.0, 1.0] if abs(d[2]) < 0.9 else [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    offs = spacing*(np.arange(nbundle) - 0.5*(nbundle - 1))
    paths = []
    for a in offs:
        for b in offs:
            q = centre + a*e1 + b*e2
            t_in = max(((0.0 if c > 0 else n) - s)/c for c, s in zip(d, q)
                       if abs(c) > 1e-12)
            t_out = min(((n if c > 0 else 0.0) - s)/c for c, s in zip(d, q)
                        if abs(c) > 1e-12)
            if t_out > t_in:
                paths.append((q + t_in*d, q + t_out*d))
    return paths


def optical_depth(p, q, centre, radius, kappa, nseg=40):
    """Points along p->q and the optical depth accumulated at each through a uniform
    sphere of opacity kappa per unit length."""
    s = np.linspace(0.0, 1.0, nseg + 1)
    pts = p[None, :] + s[:, None]*(q - p)[None, :]
    inside = np.linalg.norm(pts - centre, axis=1) < radius
    ds = np.linalg.norm(q - p)/nseg
    tau = np.concatenate(([0.0], np.cumsum(kappa*ds*0.5*(inside[1:] + inside[:-1]))))
    return pts, tau


def draw_diffuse_rays(ax, n=8, nbundle=3, pixels=(0, 5, 11), lw=0.7, kappa=0.55,
                      zoom=1.0, cloud=None, spacing=2.0):
    """With `cloud` (a clumpy_cloud.Cloud on the n^3 box) the rays cross that density
    field instead of a uniform sphere."""
    centre, radius = np.full(3, 0.5*n), 0.3*n
    for a in (0, n):
        for b in (0, n):
            ax.plot([0, n], [a, a], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [0, n], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [b, b], [0, n], color='0.55', lw=0.6)
    if cloud is None:
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        ax.plot_surface(centre[0] + radius*np.cos(u)*np.sin(v),
                        centre[1] + radius*np.sin(u)*np.sin(v),
                        centre[2] + radius*np.cos(v), color='0.45', alpha=0.25,
                        linewidth=0, shade=False)
    else:
        cloud.draw(ax)
    dirs = np.array(hp.pix2vec(1, list(pixels), nest=True)).T
    for d, col in zip(dirs, COLORS):
        rgba = np.array(to_rgba(col))
        for p, q in ray_paths(d, n, nbundle, spacing):
            if cloud is None:
                pts, tau = optical_depth(p, q, centre, radius, kappa)
            else:
                pts, tau = cloud.tau(p, q, kappa)
            for k in range(len(pts) - 1):
                c = rgba.copy()
                c[3] = 0.15 + 0.85*np.exp(-tau[k])
                ax.plot(pts[k:k + 2, 0], pts[k:k + 2, 1], pts[k:k + 2, 2], color=c,
                        lw=lw)
        tip = centre + 0.62*n*d
        ax.quiver(*(tip - 0.25*n*d), *(0.25*n*d), color=col, lw=1.6,
                  arrow_length_ratio=0.35)
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_zlim(0, n)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.view_init(elev=20, azim=-60)
    ax.set_axis_off()


def main():
    fig = plt.figure(figsize=(5.2, 5.0))
    ax = fig.add_subplot(projection='3d')
    draw_diffuse_rays(ax, nbundle=4, pixels=(0, 11), lw=1.5)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(COPY, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUT, 'diffuse_rays_3d.' + ext)
        fig.savefig(path, dpi=200)
        shutil.copy(path, COPY)
        print('wrote', path)


if __name__ == '__main__':
    main()
