"""H-alpha, [O III] 5007 and [S II] 6716+6731 images from the volume renderer
(run_lines.sh), without dust (top row) and with dust extinction (bottom row), each
with the SHO composite of ../line_maps.py ([S II] red, H-alpha green, [O III] blue).
Every line keeps one log stretch, from the brightest pixel of the dust-free images
over all frames down 3.5 decades, so the two rows compare directly.

usage: python line_frames.py LINESDIR OUTDIR [--title TEXT]
  LINESDIR: run_lines.sh output (thin/ and dust/ with lines.out1.NNNNN.h5)
Writes OUTDIR/volren_lines_NNNNN.png per frame.
"""
import argparse
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import line_maps as lm  # noqa: E402
from photchem_method import paths  # noqa: E402
sys.path.insert(0, os.path.join(paths.env('TIGRIS_DIR'), 'vis', 'volren'))
from volren_reader import read_volren  # noqa: E402

BANDS = {'Halpha': [0], 'OIII': [1], 'SII': [2, 3]}   # emit[] order in run_lines.sh


def lines(fn):
    fr = read_volren(fn)[0]
    img = fr['intensity']
    return {k: sum(img[:, :, b] for b in bs).astype(np.float64)
            for k, bs in BANDS.items()}, float(fr['sim_time'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('linesdir')
    ap.add_argument('outdir')
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    thin = sorted(glob.glob(os.path.join(a.linesdir, 'thin', 'lines.out1.*.h5')))
    dust = sorted(glob.glob(os.path.join(a.linesdir, 'dust', 'lines.out1.*.h5')))
    hi = {k: max(float(lines(fn)[0][k].max()) for fn in thin) for k in BANDS}
    rng = {k: (hi[k]*10**-3.5, hi[k]) for k in BANDS}
    for ft, fd in zip(thin, dust):
        fig, ax = plt.subplots(2, 4, figsize=(13.0, 6.9), constrained_layout=True)
        for row, fn, lab in [(0, ft, 'no dust'), (1, fd, r'dust, $Z_d = 1$')]:
            sb, t = lines(fn)
            for col, k in enumerate(['Halpha', 'OIII', 'SII']):
                im = ax[row, col].imshow(np.maximum(sb[k], rng[k][0]), origin='lower',
                                         norm=LogNorm(*rng[k]), cmap='magma',
                                         interpolation='nearest')
                ax[row, col].text(0.04, 0.95, lm.LABEL[k] + ', ' + lab, color='w',
                                  transform=ax[row, col].transAxes, va='top',
                                  fontsize=9, path_effects=[
                                      pe.withStroke(linewidth=2, foreground='k')])
                if row == 1:
                    cb = fig.colorbar(im, ax=ax[:, col], orientation='horizontal',
                                      fraction=0.04, pad=0.02)
                    cb.set_label(r'$I$ [erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$]', fontsize=8)
                    cb.ax.tick_params(labelsize=7)
            ax[row, 3].imshow(lm.compose(sb, rng, False), origin='lower',
                              interpolation='nearest')
            ax[row, 3].text(0.04, 0.95, r'[S II], H$\alpha$, [O III], ' + lab,
                            transform=ax[row, 3].transAxes, va='top', color='w',
                            fontsize=9)
        for x in ax.ravel():
            x.set_xticks([])
            x.set_yticks([])
        fig.suptitle('%s, $t = %.3f$ Myr, volume renderer, along $z$'
                     % (a.title, t*0.978), fontsize=11)
        num = os.path.basename(ft).split('.')[-2]
        out = os.path.join(a.outdir, 'volren_lines_%s.png' % num)
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print(out)


if __name__ == '__main__':
    main()
