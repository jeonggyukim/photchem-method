"""The block Gauss-Seidel schematic on the real GOW17 network, molecular gas.
jacobian_dump.txt (CO-formation test at 2.93 Myr, 64 substeps) ->
../../figures/block_gs_gow17.{pdf,png} (sensitivity) and
../../figures/block_gs_gow17_lag.{pdf,png} (lag carried).

A = I - hJ in update order.  Each off-diagonal dot is the coupling of species i to
species j, sized by |A_ij| x_j / (A_ii x_i): the fractional change of the new x_i per
fractional error in the x_j it reads (unit-free; the entries of one Jacobi iteration's
fractional error matrix).  Colour: what value of x_j the update of i reads in one
iteration -- solved together inside a block (black), the new value (green) or the
start-of-substep value (orange, the lag).  Couplings below 1e-3 are not drawn.

The second figure sizes each dot instead by the error it carries in one substep: the
coupling times delta_j, the fractional change of the species read over one substep,
|dx_j/dt| h / x_j, with dx_j/dt from the last two samples of the reference run."""
import io
import contextlib
import os
import runpy

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'figures')
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(HERE, 'plot_gs_error_matrix.py'))
J, X, perm, ORDER, LABEL, t_myr = g['J'], g['X'], g['perm'], g['ORDER'], g['LABEL'], g['t_myr']
ns = len(ORDER)
N = 64
h = 3.0/N


def couplings(kk):
    """Fractional couplings C at sample kk, the fractional change delta of each species
    over one substep (from samples kk-1 and kk), and the lag carried E = C delta."""
    A = np.eye(ns) - h*J[kk][np.ix_(perm, perm)]
    x = np.maximum(X[kk][perm], 1e-30)
    C = np.abs(A)*x[None, :]/(np.abs(np.diag(A))[:, None]*x[:, None])
    np.fill_diagonal(C, 0.0)
    x_prev = np.maximum(X[kk - 1][perm], 1e-30)
    dt_samples = (t_myr[kk] - t_myr[kk - 1])/g['myr']
    delta = np.abs(x - x_prev)/dt_samples*h/x
    return C, delta, C*delta[None, :]


k = len(t_myr) - 1
C, delta, E = couplings(k)
CUT = 1e-3
# the sample where CO changes fastest: the lag carried is largest while CO forms
ico_ = ORDER.index('CO')
k_act = 1 + int(np.argmax([couplings(kk)[1][ico_] for kk in range(1, len(t_myr))]))
C_act, delta_act, E_act = couplings(k_act)

GREY, GREEN, ORANGE = '#d9d9d9', '#2ca02c', '#ff7f0e'
ico = ORDER.index('CO')


def size(c, lo=-3.0, hi=0.0):
    """Marker size, logarithmic between 10^lo (small) and 10^hi (large)."""
    return 4 + 22*np.clip((np.log10(c) - lo)/(hi - lo), 0, 1)


def panel(ax, title, blocks, jacobi, M, cut, lo, hi):
    blk = np.zeros(ns, int)
    for b, members in enumerate(blocks):
        blk[members] = b
    for members in blocks:
        b0, b1 = min(members), max(members)
        ax.add_patch(Rectangle((b0 - 0.5, b0 - 0.5), b1 - b0 + 1, b1 - b0 + 1,
                               facecolor=GREY, edgecolor='k', lw=1.0, zorder=1))
    lag = 0.0
    for i in range(ns):
        ax.plot(i, i, 's', color='k', ms=6, zorder=3)
        for j in range(ns):
            if i == j or M[i, j] < cut:
                continue
            if blk[i] == blk[j]:
                col = 'k'
            elif jacobi or blk[j] > blk[i]:
                col = ORANGE
                lag = max(lag, M[i, j])
            else:
                col = GREEN
            ax.plot(j, i, 'o', color=col, ms=size(M[i, j], lo, hi), mec='k', mew=0.5,
                    zorder=3)
    ax.set_xlim(-0.5, ns - 0.5)
    ax.set_ylim(ns - 0.5, -0.5)
    ax.set_xticks(range(ns), [LABEL[s] for s in ORDER], rotation=90, fontsize=9)
    ax.set_yticks(range(ns), [LABEL[s] for s in ORDER], fontsize=9)
    ax.xaxis.tick_top()
    ax.set_aspect('equal')
    ax.grid(color='0.9', lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=11, pad=48)
    return lag


singles = [[i] for i in range(ns)]
grouped = [[i] for i in range(ico)] + [[ico, ico + 1]] + [[i] for i in range(ico + 2, ns)]


def make(M, cut, lo, hi, name, what, note, kk):
    fig, axs = plt.subplots(1, 3, figsize=(16, 6.8), layout='constrained')
    panel(axs[0], 'Jacobi (index-order option)\nevery species reads old values', singles,
          True, M, cut, lo, hi)
    panel(axs[1], 'Gauss–Seidel, no group\nnew values from species above', singles,
          False, M, cut, lo, hi)
    panel(axs[2], 'Gauss–Seidel + CO/HCO$^+$ group (default)\nCO and HCO$^+$ solved '
          'together', grouped, False, M, cut, lo, hi)
    handles = [plt.Line2D([], [], marker='s', ls='', color='k', ms=6, label='own term'),
               plt.Line2D([], [], marker='o', ls='', color=GREEN, mec='k', ms=9,
                          label='reads the new value'),
               plt.Line2D([], [], marker='o', ls='', color=ORANGE, mec='k', ms=9,
                          label='reads the old value (lag)'),
               plt.Line2D([], [], marker='o', ls='', color='k', ms=9,
                          label='solved together in a block')]
    for e in np.arange(lo, hi + 0.5):
        handles.append(plt.Line2D([], [], marker='o', ls='', color='0.6', mec='k', mew=0.5,
                                  ms=size(10.0**e, lo, hi), label='%g' % 10.0**e))
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle('Which value each species reads in one iteration, GOW17 in molecular gas '
                 r'(CO formation, $n_{\rm H} = 10^3\,{\rm cm^{-3}}$, no FUV, $t$ = %.2f Myr, '
                 '%d substeps)\nrow: species being updated;  column: species it reads;  '
                 'dot size: %s\n%s' % (t_myr[kk], N, what, note), fontsize=10)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG, name + '.' + ext), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', os.path.join(FIG, name + '.png'))


make(C, CUT, -3.0, 0.0, 'block_gs_gow17',
     r'fractional coupling $|A_{ij}|x_j/(A_{ii}x_i)$ (sensitivity)',
     'A large orange dot is harmless if the species read barely changes in a substep '
     '(H$_2$, CO) or if the coupling runs one way only; see block_gs_gow17_lag for the '
     'error each dot carries.', k)
hi = np.ceil(np.log10(E_act.max()))
make(E_act, 10.0**(hi - 5), hi - 5, hi, 'block_gs_gow17_lag',
     r'lag carried: coupling $\times$ fractional change of $x_j$ in one substep',
     'At the time CO forms fastest; changes from the reference run; dots below '
     '$10^{%d}$ not drawn.' % (hi - 5), k_act)
i_hco = ORDER.index('HCO+')
print('C[CO,HCO+] = %.3g, C[HCO+,CO] = %.3g' % (C[ico, i_hco], C[i_hco, ico]))
big = sorted(((C[i, j], ORDER[i], ORDER[j]) for i in range(ns) for j in range(ns)
              if j > i and C[i, j] >= 0.1), reverse=True)
print('couplings above the diagonal >= 0.1 (lag in Gauss-Seidel):', big[:8])

print('lag figure at t = %.2f Myr; fractional change per substep:' % t_myr[k_act],
      ', '.join('%s %.1e' % (ORDER[i], delta_act[i]) for i in np.argsort(-delta_act)[:6]))
bigE = sorted(((E_act[i, j], ORDER[i], ORDER[j]) for i in range(ns) for j in range(ns)
               if j > i), reverse=True)
print('largest lag carried above the diagonal:', [(round(v, 8), a, b) for v, a, b in bigE[:6]])
