"""Schematic of Jacobi, Gauss-Seidel and block Gauss-Seidel on one sparse matrix.
-> ../../figures/block_gs_schematic.{pdf,png}. Each non-zero entry a_ij (i != j) is the
coupling of unknown i to unknown j; dot size marks its strength. Colour says what
value of x_j equation i reads in one iteration: solved together (grey block), the new
value (green) or the old value (orange, the lag). x_4 and x_5 (1-based) are
strongly coupled both ways, as CO and HCO+ are."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'figures')

n = 6
S = np.zeros((n, n))
for (i, j, s) in [(1, 0, 0.3), (2, 1, 0.3), (0, 2, 0.15), (3, 2, 0.3), (4, 1, 0.15),
                  (5, 4, 0.3), (2, 5, 0.15), (5, 0, 0.15),
                  (3, 4, 1.0), (4, 3, 1.0)]:
    S[i, j] = s

GREY, GREEN, ORANGE = '#d9d9d9', '#2ca02c', '#ff7f0e'


def panel(ax, title, blocks, jacobi):
    blk = np.zeros(n, int)
    for b, members in enumerate(blocks):
        blk[members] = b
    for members in blocks:
        lo, hi = min(members), max(members)
        ax.add_patch(Rectangle((lo - 0.5, lo - 0.5), hi - lo + 1, hi - lo + 1,
                               facecolor=GREY, edgecolor='k', lw=1.2, zorder=1))
    for i in range(n):
        ax.plot(i, i, 's', color='k', ms=7, zorder=3)
        for j in range(n):
            if i == j or S[i, j] == 0:
                continue
            if blk[i] == blk[j]:
                col = 'k'
            elif jacobi or blk[j] > blk[i]:
                col = ORANGE
            else:
                col = GREEN
            ax.plot(j, i, 'o', color=col, ms=6 + 16 * S[i, j], zorder=3,
                    mec='k' if col != 'k' else 'k', mew=0.6)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(n), [f'$x_{k+1}$' for k in range(n)])
    ax.set_yticks(range(n), [f'eq. {k+1}' for k in range(n)])
    ax.xaxis.tick_top()
    ax.set_aspect('equal')
    ax.grid(color='0.85', lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=11, pad=28)


singles = [[k] for k in range(n)]
with_block = [[0], [1], [2], [3, 4], [5]]
fig, axs = plt.subplots(1, 3, figsize=(13, 5.0), layout='constrained')
panel(axs[0], 'Jacobi\n(every equation reads old values)', singles, True)
panel(axs[1], 'Gauss–Seidel\n(new values from equations above)', singles, False)
panel(axs[2], 'Block Gauss–Seidel\n($x_4$, $x_5$ solved together)', with_block, False)
handles = [plt.Line2D([], [], marker='s', ls='', color='k', ms=7, label='own term (diagonal)'),
           plt.Line2D([], [], marker='o', ls='', color=GREEN, mec='k', ms=9,
                      label='reads the new value'),
           plt.Line2D([], [], marker='o', ls='', color=ORANGE, mec='k', ms=9,
                      label='reads the old value (lag)'),
           plt.Line2D([], [], marker='o', ls='', color='k', ms=9,
                      label='inside a block (solved together)'),
           Rectangle((0, 0), 1, 1, facecolor=GREY, edgecolor='k', label='solved exactly')]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=10, frameon=False,
           bbox_to_anchor=(0.5, -0.06))
fig.suptitle('Which values each equation reads in one iteration '
             '(dot size: coupling strength; $x_4$ and $x_5$ strongly coupled both ways)',
             fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG, 'block_gs_schematic.' + ext), dpi=200, bbox_inches='tight')
print('wrote', os.path.join(FIG, 'block_gs_schematic.png'))
