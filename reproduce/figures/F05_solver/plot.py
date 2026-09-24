"""F05: error of the semi-implicit solver against the number of fixed substeps, one zone,
GOW17 core, from solver_sweep.txt (solver_sweep.cpp). Columns: A ionization, B
recombination (n_H 100, 978 yr), C CO formation (n_H 1e3, cosmic rays only, 2.9 Myr); rows:
|dT/T| and |dx/x| (x = x_H+ for A and B, x_CO for C) against the default with 200000
substeps. Writes ../../../figures/F05_solver.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F05_solver')
d = np.loadtxt(os.path.join(HERE, 'solver_sweep.txt'))
VARIANTS = [('Gauss-Seidel + CO/HCO$^+$ block (default)', 'k', '-', 'o'),
            ('Gauss-Seidel, no block', 'C1', '--', 's'),
            ('Jacobi', 'C3', ':', '^'),
            ('default + Riccati H$^+$', 'C0', '-.', 'x')]
CASES = [r'A: ionization, $n_{\rm H} = 100$, 978 yr', r'B: recombination, $n_{\rm H} = 100$, 978 yr',
         r'C: CO formation, $n_{\rm H} = 10^3$, 2.9 Myr']
XLAB = [r'$x_{\rm H^+}$', r'$x_{\rm H^+}$', r'$x_{\rm CO}$']
fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True, layout='constrained')
for c in range(3):
    m = d[:, 0] == c
    n = d[m, 1]
    for v, (lab, col, ls, mk) in enumerate(VARIANTS):
        axs[0, c].loglog(n, np.maximum(d[m, 2 + 2*v], 1e-9), color=col, ls=ls, marker=mk,
                         ms=4, label=lab)
        axs[1, c].loglog(n, np.maximum(d[m, 3 + 2*v], 1e-9), color=col, ls=ls, marker=mk, ms=4)
    for ax in axs[:, c]:
        ax.loglog(n, 0.1*n[0]/n, color='0.7', lw=0.8)
        ax.axhspan(0.01, 0.1, color='C2', alpha=0.08)
        ax.set_ylim(1e-7, 20)
        ax.grid(alpha=0.25, which='both')
    axs[0, c].set_title(CASES[c], fontsize=11)
    axs[1, c].set_xlabel('substeps (fixed, equal)')
    axs[1, c].set_ylabel(r'$|\Delta x/x|$, x = ' + XLAB[c])
    axs[0, c].set_ylabel(r'$|\Delta T/T|$')
axs[0, 0].legend(fontsize=8, loc='lower left')
fig.suptitle('Semi-implicit solver, one zone, GOW17 core: error against 200000 substeps '
             '(grey: first order; green band: 1-10% budget)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
