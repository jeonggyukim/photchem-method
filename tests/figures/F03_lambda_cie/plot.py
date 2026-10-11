"""Lambda(T) of GOW17 + multi-ion at CIE (lambda_curve.cpp output at n_H = 1), components,
against NCR's hot-gas table and a CHIANTI all-element solar total. n_H = 0.01 differs from
n_H = 1 by at most 0.9% (at 1e4 K), so one density is shown. Reads lambda_curve_n1.txt
beside this script; writes ../../../figures/F03_lambda_cie.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F03_lambda_cie')
a = np.loadtxt(os.path.join(HERE, 'lambda_curve_n1.txt'))
lT = a[:, 0]; T = 10**lT
fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 8.5), sharex=True,
                              gridspec_kw={'height_ratios': [2.4, 1]})
parts = [(3, 'H (CIE)', 'C0'), (4, 'He (CIE)', 'C1'), (5, 'tracked ions (O II-III, S II-III, N II-III)', 'C2'),
         (6, 'higher ions (C, Si, O, S, N above top stage)', 'C3'),
         (7, 'untracked (Ne, Mg, Ar, Ca, Fe) + returned C, O, Si', 'C4')]
for c, lab, col in parts:
    ax.loglog(T, np.maximum(a[:, c], 1e-30), color=col, lw=1.3, alpha=0.8, label=lab)
ax.loglog(T, a[:, 2], 'k-', lw=2.4, label=r'Tigris total ($n_{\rm H}=1$; $n_{\rm H}=0.01$ within 1%)')
ax.loglog(T, a[:, 10], color='0.45', ls='--', lw=2, label='NCR (Gnat & Ferland 2012, Asplund 2009)')
ax.loglog(T, a[:, 14], color='C5', ls='-.', lw=1.6, label='CHIANTI v11 all elements, solar')
ax.set(ylim=(1e-24, 2e-21), ylabel=r'$\Lambda_N = \mathcal{L}/n_{\rm H}^2$ [erg cm$^3$ s$^{-1}$]',
       title='CIE cooling: GOW17 + multi-ion (gas-phase C, O, Si + returned grain metals)')
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.25, which='both')
m = lT >= 4.4
axr.semilogx(T[m], a[m, 2]/a[m, 10], 'k-', lw=2, label='Tigris / NCR')
axr.semilogx(T[m], a[m, 2]/a[m, 14], color='C5', lw=1.6, label='Tigris / CHIANTI solar')
axr.axhline(1, color='0.5', lw=0.8)
axr.set(ylim=(0.5, 1.3), xlabel='T [K]', ylabel='ratio')
axr.legend(fontsize=8)
axr.grid(alpha=0.25, which='both')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
