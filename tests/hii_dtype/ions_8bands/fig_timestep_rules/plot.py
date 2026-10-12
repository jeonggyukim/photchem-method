"""dt, the domain-minimum timescales and the ionized mass against time for runs of one
setup under different timestep rules (grammar diag.sbatch: cfl10, cfl100, nolimit).

usage: python plot_timestep_rules.py DIR OUT.png
DIR holds one subdirectory per case with its HII.hst.
"""
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

CASES = [('cfl10', r'$c_{\rm pc} = 10$', '-'), ('cfl100', r'$c_{\rm pc} = 100$', '--'),
         ('nolimit', 'no limit', ':')]
KINDS = [('tmin_cool', 'cooling', 'C3'), ('tmin_Hp', r'H$^+$', 'C0'),
         ('tmin_H2', r'H$_2$', 'C2'), ('tmin_CO', 'CO', 'C4')]
MYR = 0.978  # code time unit [Myr]


def history(fn):
    hdr = [line for line in open(fn) if line.startswith('#')][-1]
    names = [m.split('=')[1] for m in re.findall(r'\[\d+\]=\S+', hdr)]
    a = np.atleast_2d(np.loadtxt(fn))
    return {n: a[:, i] for i, n in enumerate(names)}


def main():
    d, out = sys.argv[1], sys.argv[2]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
    for case, lab, ls in CASES:
        h = history('%s/%s/HII.hst' % (d, case))
        t = h['time']*MYR
        ax[0].semilogy(t[1:], h['dt'][1:], ls, color='k', lw=1, label=lab)
        if case != 'nolimit':
            for k, klab, col in KINDS:
                ax[1].semilogy(t[1:], h[k][1:], ls, color=col, lw=0.9,
                               label=klab if case == 'cfl10' else None)
        ax[2].plot(t, h['mass_ion'], ls, color='k', lw=1, marker='o' if len(t) < 5
                   else None, ms=3, label=lab)
    ax[0].set_ylabel('dt [code]')
    ax[1].set_ylabel('shortest timescale over the domain [code]')
    ax[2].set_ylabel('ionized mass [code]')
    for a in ax:
        a.set_xlabel('t [Myr]')
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8, title=r'solid $c_{\rm pc} = 10$, dashed 100',
                 title_fontsize=7, ncol=2)
    ax[2].legend(fontsize=8)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == '__main__':
    main()
