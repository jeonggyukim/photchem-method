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
# Gauss-Seidel: C_i, D_i take the species earlier in the order at n+1 and the rest at n,
# x~ = (x_{<i}^{n+1}, x_{>=i}^n); T and the ghost species stay at t^n
BE = (r'$x_i^{n+1} = \frac{x_i^n + \mathcal{C}_i(\tilde{\mathbf{x}})\,h}'
      r'{1 + \mathcal{D}_i(\tilde{\mathbf{x}})\,h}$')
EXP = (r'$x_i^{n+1} = \frac{\mathcal{C}_i}{\mathcal{D}_i} + \left(x_i^n - '
       r'\frac{\mathcal{C}_i}{\mathcal{D}_i}\right)e^{-\mathcal{D}_i h}$')
BLOCK_BE = r'$\mathbf{x}^{n+1} = (\mathbf{I} + \mathbf{M}h)^{-1}(\mathbf{x}^n + \mathbf{s}h)$'
# written as the analogue of the scalar map; the code uses e^{-Mh} x^n + h phi1(-Mh) s,
# which needs no M^-1 (det M = D_HCO+ (D_CO - c) vanishes when no carbon leaves the pair)
BLOCK_EXP = (r'$\mathbf{x}^{n+1} = \mathbf{M}^{-1}\mathbf{s} + e^{-\mathbf{M}h}'
             r'\left(\mathbf{x}^n - \mathbf{M}^{-1}\mathbf{s}\right)$')
# (label, colour, line style, marker): scalar species update; CO/HCO+ pair update
VARIANTS = [('species: ' + BE + '\nCO/HCO$^+$: ' + BLOCK_BE + '  (default)', 'k', '-', 'o'),
            ('species: ' + EXP + '\nCO/HCO$^+$: ' + BLOCK_BE, 'C0', '--', 's'),
            ('species: ' + BE + '\nCO/HCO$^+$: ' + BLOCK_EXP, 'C1', ':', '^'),
            ('species: ' + EXP + '\nCO/HCO$^+$: ' + BLOCK_EXP, 'C3', '-.', 'x')]
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
fig.legend(*axs[0, 0].get_legend_handles_labels(), loc='outside lower center', ncol=2,
           fontsize=10,
           title=r'$\tilde{\mathbf{x}} = (x_{<i}^{n+1}, x_{\geq i}^{n})$ (Gauss-Seidel order; '
           r'$T$ and ghost species at $t^n$);  '
           r'$\mathbf{x} = (x_{\rm CO}, x_{\rm HCO^+})$, '
           r'$\mathbf{M} = [[\mathcal{D}_{\rm CO}, -b], [-c, \mathcal{D}_{\rm HCO^+}]]$, '
           r'$\mathbf{s} = (a, f)$', title_fontsize=10)
fig.suptitle('Exponential map vs backward Euler (Gauss-Seidel + CO/HCO$^+$ block), one zone, '
             'GOW17 core: error against 200000 substeps (green band: 1-10%)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
