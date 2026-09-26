"""Exponential map against backward Euler, one zone, GOW17 core, from exact_map_check.txt:
|dT/T| (top) and |dx/x| (bottom) against the number of equal substeps for the three F05
cases, with backward Euler (default), exact_map, exact_block and both. Reference: the default
with 200000 substeps. Writes ../../figures/exact_map_check.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'figures', 'exact_map_check')
d = np.loadtxt(os.path.join(HERE, 'exact_map_check.txt'))
VARIANTS = [('backward Euler (default)', 'k', '-', 'o'),
            ('exact_map', 'C0', '--', 's'),
            ('exact_block', 'C1', ':', '^'),
            ('exact_map + exact_block', 'C3', '-.', 'x')]
CASES = [r'A: ionization, $n_{\rm H} = 100$, 978 yr', r'B: recombination, $n_{\rm H} = 100$, 978 yr',
         r'C: CO formation, $n_{\rm H} = 10^3$, 2.9 Myr']
XLAB = [r'$x_{\rm H^+}$', r'$x_{\rm H^+}$', r'$x_{\rm CO}$']
fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True, layout='constrained')
for c in range(3):
    m = d[:, 0] == c
    n = d[m, 1]
    for v, (lab, col, ls, mk) in enumerate(VARIANTS):
        axs[0, c].loglog(n, np.maximum(d[m, 2 + 3*v], 1e-9), color=col, ls=ls, marker=mk,
                         ms=4, label=lab)
        axs[1, c].loglog(n, np.maximum(d[m, 3 + 3*v], 1e-9), color=col, ls=ls, marker=mk, ms=4)
    for ax in axs[:, c]:
        ax.axhspan(0.01, 0.1, color='C2', alpha=0.08)
        ax.set_ylim(1e-7, 20)
        ax.grid(alpha=0.25, which='both')
    axs[0, c].set_title(CASES[c], fontsize=11)
    axs[1, c].set_xlabel('substeps (fixed, equal)')
    axs[1, c].set_ylabel(r'$|\Delta x/x|$, x = ' + XLAB[c])
    axs[0, c].set_ylabel(r'$|\Delta T/T|$')
axs[0, 0].legend(fontsize=8, loc='lower left')
fig.suptitle('Exponential map vs backward Euler (Gauss-Seidel + CO/HCO$^+$ block), one zone, '
             'GOW17 core: error against 200000 substeps (green band: 1-10%)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
