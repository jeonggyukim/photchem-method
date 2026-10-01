"""Adaptive ray tracing from a point source, after Abel & Wandelt (2002, Fig. 2): one
HEALPix base ray splits into four child rays whenever the solid angle it carries would
cover more than 1/f of a cell's face, A(l) r^2 > dx^2/f with A(l) = 4 pi/(12 4^l).

    python adaptive_rays_3d.py

Writes adaptive_rays_3d.png and .pdf to
~/Dropbox/Research/Rayt-Method/figures/photchem-postproc, with a legend-free
adaptive_rays_3d_panel for Figure 1 of the method paper, and copies them to
~/Documents/photchem-postproc.  draw_ray_tree() is also used by rad_cr_chem_coupling.py.
"""
import os
import shutil

import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Line3DCollection

OUT = os.path.expanduser('~/Dropbox/Research/Rayt-Method/figures/photchem-postproc')
COPY = os.path.expanduser('~/Documents/photchem-postproc')


def ray_tree(ncell=16, f=1.0, src=(5.0, 6.0, 5.0), narrow=1, aim=(0.8, 0.5, 0.35)):
    """Segments (start, end, level) and split points of the rays grown from the base
    pixel of level 0 that points closest to `aim`, for a source at `src` in a grid of
    ncell^3 cells of unit width.  A ray splits wherever the criterion says, until it
    leaves the grid.  To keep the drawing readable, at the first `narrow` splits only
    the child nearest `aim` is followed, so the rays drawn fill 1/(12 4^narrow) of the
    sky."""
    src = np.asarray(src, dtype=float)
    diag = np.asarray(aim, dtype=float)/np.linalg.norm(aim)
    vec0 = np.array(hp.pix2vec(1, np.arange(12), nest=True)).T
    base = int(np.argmax(vec0 @ diag))
    segs, splits = [], []
    stack = [(0, base, 0.0)]
    while stack:
        lev, ipix, r0 = stack.pop()
        v = np.array(hp.pix2vec(2**lev, ipix, nest=True))
        r_exit = min(((ncell if c > 0 else 0.0) - s)/c for c, s in zip(v, src)
                     if abs(c) > 1e-12)
        r_split = 1.0/np.sqrt(f*4.0*np.pi/(12*4**lev))
        r1 = min(r_split, r_exit)
        if r1 > r0:
            segs.append((src + r0*v, src + r1*v, lev))
        if r_split < r_exit:
            splits.append((src + r_split*v, lev))
            children = [4*ipix + k for k in range(4)]
            if lev < narrow:
                cv = np.array(hp.pix2vec(2**(lev + 1), children, nest=True)).T
                children = [children[int(np.argmax(cv @ diag))]]
            stack.extend((lev + 1, c, r_split) for c in children)
    return segs, splits


def draw_ray_tree(ax, ncell=16, f=1.0, lw=0.7, grid=True, cmap='plasma', zoom=1.0,
                  ms=1.5, narrow=1, sources=None, cloud=None, kappa=0.5):
    """Draw the grid's edges and the ray tree on a 3D axes.  `sources` is a list of
    (position, aim) pairs, one tree each; with `cloud` (a clumpy_cloud.Cloud on the
    ncell^3 box) the cloud is drawn and each ray fades with the optical depth it has
    crossed since leaving its source."""
    if sources is not None or cloud is not None:
        return draw_ray_forest(ax, ncell, f, lw, grid, cmap, zoom, ms, narrow,
                               sources, cloud, kappa)
    segs, splits = ray_tree(ncell, f, narrow=narrow)
    lmax = max(s[2] for s in segs)
    src = segs[0][0] if segs[0][2] == 0 else np.full(3, 0.5*ncell)
    cm = plt.get_cmap(cmap)
    col = [cm(0.1 + 0.75*lev/lmax) for lev in range(lmax + 1)]
    n = ncell
    for a in (0, n):
        for b in (0, n):
            ax.plot([0, n], [a, a], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [0, n], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [b, b], [0, n], color='0.55', lw=0.6)
    if grid:
        g = np.arange(2, n, 2)
        for c in g:
            ax.plot([c, c], [n, n], [0, n], color='0.85', lw=0.4)
            ax.plot([0, n], [n, n], [c, c], color='0.85', lw=0.4)
            ax.plot([0, 0], [c, c], [0, n], color='0.85', lw=0.4)
            ax.plot([0, 0], [0, n], [c, c], color='0.85', lw=0.4)
            ax.plot([c, c], [0, n], [0, 0], color='0.85', lw=0.4)
            ax.plot([0, n], [c, c], [0, 0], color='0.85', lw=0.4)
    for p, q, lev in segs:
        ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], color=col[lev], lw=lw)
    for p, lev in splits:
        ax.scatter(*p, color=col[lev], s=ms, depthshade=False)
    ax.scatter(*src, marker='*', s=160, color='#FFD23F', edgecolor='#C77800',
               depthshade=False, zorder=10)
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_zlim(0, n)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.view_init(elev=20, azim=-60)
    ax.set_axis_off()
    return col


def draw_box(ax, n, grid):
    for a in (0, n):
        for b in (0, n):
            ax.plot([0, n], [a, a], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [0, n], [b, b], color='0.55', lw=0.6)
            ax.plot([a, a], [b, b], [0, n], color='0.55', lw=0.6)
    if grid:
        for c in np.arange(2, n, 2):
            ax.plot([c, c], [n, n], [0, n], color='0.85', lw=0.4)
            ax.plot([0, n], [n, n], [c, c], color='0.85', lw=0.4)
            ax.plot([0, 0], [c, c], [0, n], color='0.85', lw=0.4)
            ax.plot([0, 0], [0, n], [c, c], color='0.85', lw=0.4)
            ax.plot([c, c], [0, n], [0, 0], color='0.85', lw=0.4)
            ax.plot([0, n], [c, c], [0, 0], color='0.85', lw=0.4)


def draw_ray_forest(ax, ncell, f, lw, grid, cmap, zoom, ms, narrow, sources, cloud,
                    kappa, nsub=10):
    """Several ray trees, optionally through a cloud; see draw_ray_tree."""
    if sources is None:
        sources = [((5.0, 6.0, 5.0), (0.8, 0.5, 0.35))]
    n = ncell
    draw_box(ax, n, grid)
    if cloud is not None:
        cloud.draw(ax)
    # a source may carry its own `narrow` as a third element
    trees = [(np.asarray(sa[0], dtype=float),)
             + ray_tree(n, f, sa[0], sa[2] if len(sa) > 2 else narrow, sa[1])
             for sa in sources]
    lmax = max(sg[2] for _, segs, _ in trees for sg in segs)
    cm = plt.get_cmap(cmap)
    col = [np.array(cm(0.1 + 0.75*lev/lmax)) for lev in range(lmax + 1)]
    lines, colors = [], []
    for src, segs, splits in trees:
        for p, q, lev in segs:
            if cloud is None:
                lines.append([p, q])
                colors.append(col[lev])
                continue
            tau0 = cloud.tau(src, p, kappa, nseg=20)[1][-1]
            pts, tau = cloud.tau(p, q, kappa, nseg=nsub)
            for k in range(nsub):
                c = col[lev].copy()
                c[3] = 0.12 + 0.88*np.exp(-(tau0 + tau[k]))
                lines.append([pts[k], pts[k + 1]])
                colors.append(c)
        for p, lev in splits:
            ax.scatter(*p, color=col[lev], s=ms, depthshade=False)
    ax.add_collection3d(Line3DCollection(lines, colors=colors, linewidths=lw))
    for src, _, _ in trees:
        ax.scatter(*src, marker='*', s=160, color='#FFD23F', edgecolor='#C77800',
                   depthshade=False, zorder=10)
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_zlim(0, n)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.view_init(elev=20, azim=-60)
    ax.set_axis_off()
    return col


def main():
    fig = plt.figure(figsize=(5.2, 5.0))
    ax = fig.add_subplot(projection='3d')
    col = draw_ray_tree(ax)
    segs, _ = ray_tree()
    for lev, c in enumerate(col):
        n = sum(1 for s in segs if s[2] == lev)
        ax.plot([], [], color=c, lw=2, label='level %d: %d ray%s drawn' % (
            lev, n, '' if n == 1 else 's'))
    ax.legend(loc='upper left', fontsize=8, frameon=False)
    fig.tight_layout(pad=0.1)
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(COPY, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUT, 'adaptive_rays_3d.' + ext)
        fig.savefig(path, dpi=200)
        shutil.copy(path, COPY)
        print('wrote', path)
    # Figure 1 of the method paper shrinks this to a 34 mm panel, where the
    # legend would be 2 pt type.
    ax.get_legend().remove()
    for ext in ('png', 'pdf'):
        path = os.path.join(OUT, 'adaptive_rays_3d_panel.' + ext)
        fig.savefig(path, dpi=200, bbox_inches='tight', pad_inches=0.02)
        shutil.copy(path, COPY)
        print('wrote', path)


if __name__ == '__main__':
    main()
