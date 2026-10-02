"""F11: where the interior temperature of the static H II region comes from. T(r) and
x(H+)(r) of the Tigris post-processed sphere (GOW17 + O3,S3,N3, 7 bands) with the
recombination photons of H and He absorbed on the spot (recombination_ots = true: H
case B, He on-the-spot split) and carried by the diffuse tracer (false: H and He case A,
case A recombination cooling), against Cloudy 25 with `diffuse OTS` and with its default
transfer of the diffuse field; each Tigris run is paired with the Cloudy run of the same
treatment, and the bottom panel gives their T ratio. Writes
../../../figures/F11_interior_temperature.{pdf,png}.

Inputs: run_recomb.sh output in WORKDIR (ots_true, ots_false); Cloudy radial_profiles.txt
of reproduce/cloudy_stromgren (default) and reproduce/cloudy_stromgren/thermal_with_N_ots.
The Tigris profiles come from reduced<SUFFIX>.txt when it exists;
`python plot.py [SUFFIX] --from-runs` reads WORKDIR
($PHOTCHEM_RUNS/T6_multi_ion/recomb_ots) and rewrites it."""
import glob
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
import reduced  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.join(HERE, '..', '..')
ARGS = [a for a in sys.argv[1:] if a != '--from-runs']
# optional suffix of the run directories and the figure, e.g. _r20 for the runs with
# rayt_point/rays_per_cell = 20 and tau_max = 30
SUFFIX = ARGS[0] if ARGS else ''
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F11_interior_temperature' + SUFFIX)
REDUCED = os.path.join(HERE, 'reduced%s.txt' % SUFFIX)
CLOUDY_DEFAULT = os.path.join(REPRO, 'cloudy_stromgren', 'radial_profiles.txt')
CLOUDY_OTS = os.path.join(REPRO, 'cloudy_stromgren', 'thermal_with_N_ots',
                          'radial_profiles.txt')
# (Tigris run, its label, Cloudy file, Cloudy label, colour)
PAIRS = [('ots_true', 'Tigris, recombination_ots = true (on the spot)', CLOUDY_OTS,
          'Cloudy 25, diffuse OTS', 'C0'),
         ('ots_false', 'Tigris, recombination_ots = false (diffuse tracer)', CLOUDY_DEFAULT,
          'Cloudy 25, default', 'C3')]
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
# rho is in m_H cm^-3; with mu_H in the T formula the factor 1.4 cancels
MUH, KB = 2.34335276e-24, 1.380649e-16


def tigris_profile(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    rho = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    T = d['press']*MUH*1e10/KB/(rho*(1.0 - x['H2'] + 0.1 + xe))
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    # bins a quarter of a cell wide; the cell centres fall at discrete radii, so some
    # bins hold no cell and are dropped
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, 0.25*dx)
    ri = np.digitize(R.ravel(), rb)
    w = rho.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    full = den > 0.0

    def radial(f):
        return (np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)])[full]/den[full]

    return (0.5*(rb[1:] + rb[:-1]))[full], radial(T), radial(x['H+'])


def radius_half(r, xhp):
    i = np.nonzero(xhp < 0.5)[0][0]
    return np.interp(0.5, [xhp[i], xhp[i - 1]], [r[i], r[i - 1]])


if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    data = {}
    for run, *_ in PAIRS:
        data['r_' + run], data['T_' + run], data['xhp_' + run] = tigris_profile(
            paths.runs('T6_multi_ion', 'recomb_ots', run + SUFFIX))
    reduced.save(REDUCED, 'F11: density-weighted radial profiles of the Tigris runs '
                 'T6_multi_ion/recomb_ots/ots_{true,false}%s\nr [pc], T [K], xhp = x_H+'
                 % SUFFIX, data)
else:
    data = reduced.load(REDUCED)

fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7, 9.5), sharex=True,
                                 gridspec_kw={'height_ratios': [2, 1, 1]})
for run, lab, cfile, clab, col in PAIRS:
    r, T, xhp = data['r_' + run], data['T_' + run], data['xhp_' + run]
    c = np.loadtxt(cfile)
    a1.plot(c[:, 0], c[:, 1], color=col, lw=2, alpha=0.6,
            label='%s: $R_s$ = %.3f pc' % (clab, radius_half(c[:, 0], c[:, 3])))
    a1.plot(r, T, 'o', color=col, ms=2.5, label='%s: $R_s$ = %.3f pc' % (lab, radius_half(r, xhp)))
    a2.semilogy(c[:, 0], c[:, 3], color=col, lw=2, alpha=0.6)
    a2.semilogy(r, xhp, 'o', color=col, ms=2.5)
    # inside the ionization front only: across it T drops by 10^3 over 0.05 pc
    m = (r > 0.1) & (r < 2.8)
    ratio = T[m]/np.interp(r[m], c[:, 0], c[:, 1])
    a3.plot(r[m], ratio, 'o-', color=col, ms=2, lw=0.8, label='%s / %s' % (lab.split(' (')[0], clab))
    print('%-9s T at 0.3, 1.0, 2.0, 2.7 pc:' % run, np.interp([0.3, 1.0, 2.0, 2.7], r, T).round(0),
          ' ratio to %s at the same radii:' % clab,
          (np.interp([0.3, 1.0, 2.0, 2.7], r, T)
           / np.interp([0.3, 1.0, 2.0, 2.7], c[:, 0], c[:, 1])).round(3))
a1.set(ylabel='T [K]', ylim=(4000, 11000),
       title=r'Static H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, SB99 2 Myr')
a1.legend(fontsize=7.5, loc='upper left')
a2.set(ylabel=r'$x_{\rm H^+}$', ylim=(1e-3, 1.5))
a3.axhline(1, color='0.5', lw=0.8)
a3.set(ylabel='T ratio', xlabel='r [pc]', ylim=(0.8, 1.15), xlim=(0, 3.5))
a3.legend(fontsize=7.5, loc='lower left')
for a in (a1, a2, a3):
    a.grid(alpha=0.25)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
