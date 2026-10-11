"""F02: collisional-equilibrium stage fractions of C, Si, O, S, N from the GOW17 network
with --photchem_ions=O3,S3,N3 (fixed T, no radiation, steady state) against CHIANTI v11,
and the metal electrons per H. Reads cie_fractions.txt; writes
../../../figures/F02_cie_fractions.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F02_cie_fractions')
d = np.loadtxt(os.path.join(HERE, 'cie_fractions.txt'))
T = 10**d[:, 0]
# element, spectroscopic labels of neutral .. top tracked, then the higher ions
ELEMENTS = [('C', ['C I', 'C II'], r'C$^{\geq 2+}$'),
            ('Si', ['Si I', 'Si II'], r'Si$^{\geq 2+}$'),
            ('O', ['O I', 'O II', 'O III'], r'O$^{\geq 3+}$'),
            ('S', ['S I', 'S II', 'S III'], r'S$^{\geq 3+}$'),
            ('N', ['N I', 'N II', 'N III'], r'N$^{\geq 3+}$')]
fig, axs = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
axs = axs.ravel()
col = 1
for ax, (el, stages, high) in zip(axs, ELEMENTS):
    n = len(stages) + 1
    tig = d[:, col:col + n]
    chi = d[:, col + n:col + 2*n]
    col += 2*n
    for q, lab in enumerate(stages + [high]):
        ax.loglog(T, np.maximum(chi[:, q], 1e-30), color='C%d' % q, lw=2.2, alpha=0.45)
        ax.loglog(T, np.maximum(tig[:, q], 1e-30), color='C%d' % q, ls='--', lw=1.4,
                  label=lab)
    ax.set(ylim=(1e-3, 1.5), title=el)
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(alpha=0.25, which='both')
ax = axs[5]
ax.loglog(T, d[:, col + 1], color='0.5', lw=2.2, alpha=0.6, label='CHIANTI v11 CIE')
ax.loglog(T, d[:, col], 'k--', lw=1.4, label='GOW17 + ions (tracked + higher ions)')
ax.set(title=r'metal electrons per H (C, Si, O, S, N)', ylim=(1e-7, 1e-2))
ax.legend(fontsize=8, loc='lower right')
ax.grid(alpha=0.25, which='both')
for ax in axs[3:]:
    ax.set_xlabel('T [K]')
for ax in (axs[0], axs[3]):
    ax.set_ylabel('fraction of the element')
fig.suptitle('Fixed T, no radiation, steady state: thick lines CHIANTI v11 CIE (stages above '
             'the top tracked one summed), dashed GOW17 + ions', fontsize=11)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
