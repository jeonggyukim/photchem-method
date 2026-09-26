"""Figures for the Gauss-Seidel order tests, one zone, GOW17 core: error against the
number of equal substeps, measured against the default solver at 200000 substeps.
  h2_order_check.txt -> ../../figures/h2_order_check.{pdf,png}: H2 last (default) vs first
  h_block_check.txt  -> ../../figures/h_block_check.{pdf,png}: default vs H+/H2 2x2 block"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'figures')
CASE_TITLE = [r'A: ionization of atomic gas, $n_{\rm H} = 100$, 978 yr',
              r'B: recombination, $n_{\rm H} = 100$, 978 yr',
              r'C: CO formation, $n_{\rm H} = 10^3$, 2.9 Myr',
              r'D: photodissociation of H$_2$, $n_{\rm H} = 100$, $\chi = 1$, 978 yr',
              r'E: ionization front into H$_2$, $n_{\rm H} = 100$, 978 yr']


def panel(ax, n, series, ylabel):
    for (lab, col, ls, mk), y in series:
        ax.loglog(n, np.maximum(y, 1e-9), color=col, ls=ls, marker=mk, ms=4, label=lab)
    ax.axhspan(0.01, 0.1, color='C2', alpha=0.08)
    ax.set_ylim(1e-7, 1e3)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25, which='both')


# H2 order: columns case N, then per variant |dT/T| |dy/y|
d = np.loadtxt(os.path.join(HERE, 'h2_order_check.txt'))
ylab = [r'$y_{\rm H^+}$', r'$y_{\rm H^+}$', r'$y_{\rm CO}$', r'$y_{\rm H_2}$']
styles = [('H$_2$ last (default)', 'k', '-', 'o'), ('H$_2$ first', 'C3', '--', 's')]
fig, axs = plt.subplots(2, 4, figsize=(19, 8), sharex=True, layout='constrained')
for c in range(4):
    m = d[:, 0] == c
    n = d[m, 1]
    panel(axs[0, c], n, [(s, d[m, 2 + 2*v]) for v, s in enumerate(styles)], r'$|\Delta T/T|$')
    panel(axs[1, c], n, [(s, d[m, 3 + 2*v]) for v, s in enumerate(styles)],
          r'$|\Delta y/y|$, $y$ = ' + ylab[c])
    axs[0, c].set_title(CASE_TITLE[c], fontsize=10)
    axs[1, c].set_xlabel('substeps (fixed, equal)')
axs[0, 0].legend(fontsize=9)
fig.suptitle('Where H$_2$ sits in the Gauss-Seidel update: error against 200000 substeps '
             '(green band: 1-10%)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG, 'h2_order_check.' + ext), dpi=300)

# H+/H2 block: columns case N, then per variant |dT/T| |dH+| |dH2| |dCO|
d = np.loadtxt(os.path.join(HERE, 'h_block_check.txt'))
styles = [('default (H$^+$ in the sweep, H$_2$ last)', 'k', '-', 'o'),
          ('H$^+$/H$_2$ 2$\\times$2 block', 'C0', '--', 's')]
rows = [(0, r'$|\Delta T/T|$'), (1, r'$|\Delta y/y|$, $y = y_{\rm H^+}$'),
        (2, r'$|\Delta y/y|$, $y = y_{\rm H_2}$')]
fig, axs = plt.subplots(3, 5, figsize=(23, 11), sharex=True, layout='constrained')
for c in range(5):
    m = d[:, 0] == c
    n = d[m, 1]
    for r, (q, lab) in enumerate(rows):
        panel(axs[r, c], n, [(s, d[m, 2 + 4*v + q]) for v, s in enumerate(styles)], lab)
    axs[0, c].set_title(CASE_TITLE[c], fontsize=10)
    axs[-1, c].set_xlabel('substeps (fixed, equal)')
axs[0, 0].legend(fontsize=9)
fig.suptitle('H$^+$ and H$_2$ solved together or in sequence: error against 200000 '
             'substeps (green band: 1-10%)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG, 'h_block_check.' + ext), dpi=300)
print('wrote', os.path.join(FIG, 'h2_order_check.png'), os.path.join(FIG, 'h_block_check.png'))
