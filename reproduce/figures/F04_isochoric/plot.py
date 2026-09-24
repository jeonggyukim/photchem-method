"""F04: isochoric cooling from CIE at 1e7 K, n_H = 1 cm^-3: stage fractions against T for
C, Si, O, S, N from the GOW17 network with --photchem_ions=O3,S3,N3 (dashed) and from every
stage integrated along the same T(t), n_e(t), x_H(t) with the same rates (solid; stages
above the top tracked one summed). The last panel shows the full ladder's O4+..O7+ and the
temperature history. Reads isochoric_n1.txt and fullladder_n1.txt; writes
../../../figures/F04_isochoric.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F04_isochoric')
tg = np.loadtxt(os.path.join(HERE, 'isochoric_n1.txt'))
fl = np.loadtxt(os.path.join(HERE, 'fullladder_n1.txt'))
T = tg[:, 1]
# element, total per H, tigris columns of the tracked ions then X_high, full-ladder
# first column and Z, stage labels (neutral, tracked..., higher ions)
ELEMENTS = [('C', 1.6e-4, [7], 15, 2, 6, ['C I', 'C II', r'C$^{\geq 2+}$']),
            ('Si', 1.7e-6, [14], 16, 9, 14, ['Si I', 'Si II', r'Si$^{\geq 2+}$']),
            ('O', 3.2e-4, [13, 17], 18, 24, 8, ['O I', 'O II', 'O III', r'O$^{\geq 3+}$']),
            ('S', 1.45e-5, [19, 20], 21, 33, 16, ['S I', 'S II', 'S III', r'S$^{\geq 3+}$']),
            ('N', 7.4e-5, [22, 23], 24, 50, 7, ['N I', 'N II', 'N III', r'N$^{\geq 3+}$'])]
fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
axs = axs.ravel()
for ax, (el, xtot, tcols, hcol, f0, z, labels) in zip(axs, ELEMENTS):
    ntr = len(tcols)
    tig = [tg[:, c]/xtot for c in tcols] + [tg[:, hcol]/xtot]
    tig = [np.clip(1.0 - np.sum(tig, axis=0), 0.0, None)] + tig
    ref = [fl[:, f0 + q] for q in range(ntr + 1)] + [fl[:, f0 + ntr + 1:f0 + z + 1].sum(1)]
    for q, lab in enumerate(labels):
        ax.semilogx(T, ref[q], color='C%d' % q, lw=2.2, alpha=0.45)
        ax.semilogx(T, tig[q], color='C%d' % q, ls='--', lw=1.4, label=lab)
    ax.set(title=el, xlim=(1e7, 9e3), ylim=(0, 1.02), xscale='log')
    ax.legend(fontsize=8, loc='center right')
    ax.grid(alpha=0.25, which='both')
ax = axs[5]
for q in range(4, 8):
    ax.semilogx(T, fl[:, 24 + q], color='C%d' % (q - 4), lw=1.8, label='O$^{%d+}$ (full ladder)' % q)
ax.set(title='full ladder: the higher ions of O', xlim=(1e7, 9e3), ylim=(0, 1.02))
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.25, which='both')
ax2 = ax.twinx()
ax2.loglog(T, tg[:, 0], 'k:', lw=1.2)
ax2.set_ylabel('t [code time = 0.978 Myr] (dotted)', fontsize=8)
for ax in axs[3:]:
    ax.set_xlabel('T [K] (decreasing: time runs to the right)')
for ax in (axs[0], axs[3]):
    ax.set_ylabel('fraction of the element')
fig.suptitle(r'Isochoric cooling from CIE at $10^7$ K, $n_{\rm H} = 1$ cm$^{-3}$: dashed GOW17 + '
             'ions, thick every stage along the same T(t) (stages above the top summed)',
             fontsize=11)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
