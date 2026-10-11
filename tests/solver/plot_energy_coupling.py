"""Energy-species coupling within a substep, CO-formation case of F05.
energy_coupling.txt -> ../../figures/energy_coupling.{pdf,png}: T and x_CO against time for
the default (energy and species both read the start-of-substep state), refresh_rates (species
read rates at the updated T) and n_iter = 2, at 16, 32 and 64 equal substeps, against the
reference (default, 200000 substeps)."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'figures')
d = np.loadtxt(os.path.join(HERE, 'energy_coupling.txt'))
mode, N, t, T, xco = d[:, 0], d[:, 1], d[:, 3], d[:, 4], d[:, 6]
ref = (N == 0)
styles = {0: ('default: energy and species from start-of-substep state', 'C0', 'o', '-'),
          1: ('refresh_rates: species use rates at the updated $T$', 'C3', 's', '--'),
          2: ('n_iter = 2: species update repeated', 'C2', '^', ':')}
cols = [16, 32, 64]
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True, layout='constrained')
for c, n in enumerate(cols):
    for r, (y, lab) in enumerate([(T, r'$T$ [K]'), (xco / 1.6e-4, r'$x_{\rm CO}/x_{\rm C,tot}$')]):
        ax = axs[r, c]
        ax.plot(t[ref], y[ref], color='k', lw=2.5, alpha=0.4,
                label=r'reference ($2\times10^5$ substeps)')
        for m, (lab_m, col, mk, ls) in styles.items():
            sel = (mode == m) & (N == n)
            ax.plot(t[sel], y[sel], color=col, marker=mk, ms=4, ls=ls, lw=1.2, label=lab_m)
        if r == 0:
            ax.set_yscale('log')
            ax.set_title(f'{n} equal substeps', fontsize=11)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.25, which='both')
    axs[1, c].set_xlabel('t [Myr]')
axs[0, 0].legend(fontsize=8, loc='upper right')
fig.suptitle(r'CO formation, $n_{\rm H} = 10^3\,{\rm cm^{-3}}$, cosmic rays only: '
             'coupling of the energy and species updates within a substep', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG, 'energy_coupling.' + ext), dpi=200)
print('wrote', os.path.join(FIG, 'energy_coupling.png'))
