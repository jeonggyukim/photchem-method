"""dt, the shortest timescales over the domain and the ionized mass against time for
one setup under different timestep rules: the n_H = 1e3 cm^-3 dusty GOW17 H II
region from a molecular start at 128^3. Two sets:
  v1  grammar/diag.sbatch to t = 0.02: cfl_photchem = 10 and 100 with the cooling and
      H+ limits, and no chemistry limit on dt (flag_update_dt_main = false, before
      the first step was limited); writes ../../../../figures/timestep_rules.png
  v2  grammar/diag2.sbatch to t = 0.05: the cooling, H+, H2 and CO limits, the H+
      limit alone (dt_limit = Hp), and the chemistry limiting the first step only;
      writes ../../../../figures/timestep_rules_v2.png

The curves come from reduced_n128.txt (v1) or reduced_v2_n128.txt (v2);
`python plot.py SET --from-runs DIR` reads DIR/<case>/HII.hst (DIR = the diag/n128 or
diag/v2_n128 directory of the grammar job, copied here) and rewrites it.
Usage: python plot.py [v1|v2] [--from-runs DIR]   (default v1)
"""
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from photchem_method import reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SET = 'v2' if 'v2' in sys.argv[1:] else 'v1'
OUT = os.path.join(HERE, '..', '..', '..', '..', 'figures',
                   'timestep_rules' + ('_v2' if SET == 'v2' else ''))
REDUCED = os.path.join(HERE, 'reduced_v2_n128.txt' if SET == 'v2' else 'reduced_n128.txt')
# (case, label, line style); the first two also draw their timescales
CASES = {'v1': [('cfl10', 'cfl_photchem = 10', '-'),
                ('cfl100', 'cfl_photchem = 100', '--'), ('nolimit', 'no limit', ':')],
         'v2': [('all', 'cooling, H+, H2, CO', '-'), ('hponly', 'H+ only', '--'),
                ('nolimit', 'first step only', ':')]}[SET]
KINDS = [('tmin_cool', 'cooling', 'C3'), ('tmin_Hp', r'H$^+$', 'C0'),
         ('tmin_H2', r'H$_2$', 'C2'), ('tmin_CO', 'CO', 'C4')]
FIELDS = ['time', 'dt', 'mass_ion'] + [k for k, _, _ in KINDS]
MYR = 0.978  # code time unit [Myr]


def history(fn):
    hdr = [line for line in open(fn) if line.startswith('#')][-1]
    names = [m.split('=')[1] for m in re.findall(r'\[\d+\]=\S+', hdr)]
    a = np.atleast_2d(np.loadtxt(fn))
    return {n: a[:, i] for i, n in enumerate(names)}


def from_runs(d):
    arrays = {}
    for case, _, _ in CASES:
        h = history(os.path.join(d, case, 'HII.hst'))
        for f in FIELDS:
            arrays['%s_%s' % (case, f)] = h[f]
    reduced.save(REDUCED, 'history columns of each case: <case>_<column>', arrays,
                 sources=[os.path.join(d, c) for c, _, _ in CASES])


def main():
    if '--from-runs' in sys.argv:
        from_runs(sys.argv[sys.argv.index('--from-runs') + 1])
    r = reduced.load(REDUCED)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
    for case, lab, ls in CASES:
        t = r[case + '_time']*MYR
        mk = 'o' if len(t) < 5 else None
        ax[0].semilogy(t, r[case + '_dt'], ls, color='k', lw=1, marker=mk, ms=3,
                       label=lab)
        if case != 'nolimit':
            for k, klab, col in KINDS:
                ax[1].semilogy(t[1:], r['%s_%s' % (case, k)][1:], ls, color=col, lw=0.9,
                               label=klab if case == CASES[0][0] else None)
        ax[2].plot(t, r[case + '_mass_ion'], ls, color='k', lw=1, marker=mk, ms=3,
                   label=lab)
    ax[0].set_ylabel('dt [code]')
    ax[1].set_ylabel('shortest timescale over the domain [code]')
    ax[2].set_ylabel('ionized mass [code]')
    for a in ax:
        a.set_xlabel('t [Myr]')
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8, title='solid: %s, dashed: %s' % (CASES[0][1], CASES[1][1]),
                 title_fontsize=7, ncol=2)
    ax[2].legend(fontsize=8)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for ext in figure_formats():
        fig.savefig(OUT + '.' + ext, dpi=110)
    print(OUT)


if __name__ == '__main__':
    main()
