"""The adaptive ray-tracing panel of ncr_coupling.py on its own, for a slide: rays from
two star clusters split as they spread through a clumpy cloud.

    python rays_adaptive_talk.py

Writes rays_adaptive_talk.png to ../figures. Uses draw_ray_tree from
rays_adaptive.py and the same sources and cloud as ncr_cr_coupling.py.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Bbox
from photchem_method.paths import figure_formats  # noqa: E402

sys.argv = sys.argv[:1]   # ncr_cr_coupling reads its style from argv at import
from ncr_cr_coupling import POINT_SOURCES, make_cloud  # noqa: E402
from rays_adaptive import draw_ray_tree  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def main():
    fig = plt.figure(figsize=(7.0, 7.0))
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.set_facecolor('none')
    cloud = make_cloud(16.0)
    # denser and larger cloud points than the shared default, for a slide
    base_draw = cloud.draw
    cloud.draw = lambda a: base_draw(a, npts=14000, alpha=0.07, ms=9.0)
    draw_ray_tree(ax, lw=0.9, zoom=1.05, ms=0.6, sources=POINT_SOURCES, cloud=cloud,
                  kappa=0.25)
    os.makedirs(OUT, exist_ok=True)
    # A 3D axes fills its whole rectangle, so bbox_inches='tight' keeps the blank space
    # around the cube; crop to the drawn pixels instead and reuse that box for the PDF.
    dpi = 250
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    ys, xs = np.nonzero((img < 245).any(axis=2))
    pad = 0.02*dpi
    w, h = fig.get_size_inches()*fig.dpi
    x0, x1 = (xs.min() - pad)/fig.dpi, (xs.max() + pad)/fig.dpi
    y0, y1 = (h - ys.max() - pad)/fig.dpi, (h - ys.min() + pad)/fig.dpi
    box = Bbox([[max(x0, 0), max(y0, 0)], [min(x1, w/fig.dpi), min(y1, h/fig.dpi)]])
    for ext in figure_formats():
        path = os.path.normpath(os.path.join(OUT, 'rays_adaptive_talk.' + ext))
        fig.savefig(path, dpi=dpi, bbox_inches=box)
        print('wrote', path)


if __name__ == '__main__':
    main()
