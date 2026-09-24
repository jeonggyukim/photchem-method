"""F14: cost per cycle, split by task, for the thermal/chemical modes on hii_dtype (64^3,
4 ranks, to t = 1 code, no HDF5) and rad_snr (64^3, 8 ranks, to t = 0.05 code), from
cost_table.txt; cycles and total per-cycle cost x cycles written above each bar.
Writes ../../../figures/F14_cost.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F14_cost')
rows = [ln.split() for ln in open(os.path.join(HERE, 'cost_table.txt'))
        if not ln.startswith('#')]
PARTS = [('photochemistry', 4, 'C0'), ('hydro', 5, 'C7'), ('operator split', 6, 'C2'),
         ('ray tracing', 7, 'C1'), ('other (output)', 8, '0.85')]
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), layout='constrained')
for ax, prob, title in [(axs[0], 'hii_dtype', 'hii_dtype, 4 ranks, to t = 1 code, no HDF5'),
                        (axs[1], 'rad_snr', 'rad_snr, 8 ranks, to t = 0.05 code '
                                            '(GOW17 runs write HDF5)')]:
    sel = [r for r in rows if r[0] == prob]
    x = np.arange(len(sel))
    bottom = np.zeros(len(sel))
    for lab, col, c in PARTS:
        v = np.array([float(r[col]) for r in sel])
        ax.bar(x, v, bottom=bottom, color=c, label=lab, width=0.6)
        bottom += v
    for i, r in enumerate(sel):
        ax.text(i, bottom[i]*1.02, '%s cycles\n%.1f s in loop' % (r[2], float(r[3])*int(r[2])),
                ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x, [r[1].replace('_', ' ') for r in sel], fontsize=9)
    ax.set(ylabel='rank-0 time per cycle [s]', title=title, ylim=(0, bottom.max()*1.3))
    ax.grid(alpha=0.25, axis='y')
axs[0].legend(fontsize=8, loc='upper left')
fig.suptitle(r'Cost per cycle by task, $64^3$', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
