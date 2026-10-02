"""F10: how many point-source bands the static H II region needs. T(r) and x(H+)(r) of
the Tigris post-processed sphere (GOW17 + O3,S3,N3, 64^3) with 3, 5, 6 and 7 bands against
Cloudy 25 (default and `diffuse OTS`); R_s (x_H+ = 0.5) and wall time per run in the
legend. Writes ../../../figures/F10_band_count.{pdf,png}.

Inputs: run_series.sh output in WORKDIR (run64_b3 .. run64_b7, each with time.out);
Cloudy radial_profiles.txt of reproduce/cloudy_stromgren and its thermal_with_N_ots.
The Tigris profiles and wall times come from reduced.txt when it exists;
`python plot.py --from-runs` reads WORKDIR ($PHOTCHEM_RUNS/T6_multi_ion/F10_bands) and
rewrites it."""
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
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F10_band_count')
REDUCED = os.path.join(HERE, 'reduced.txt')
CLOUDY = {'Cloudy 25, default': os.path.join(REPRO, 'cloudy_stromgren', 'radial_profiles.txt'),
          'Cloudy 25, diffuse OTS': os.path.join(REPRO, 'cloudy_stromgren',
                                                 'thermal_with_N_ots', 'radial_profiles.txt')}
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
MUH, KB = 2.34335276e-24, 1.380649e-16


def tigris_profile(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    nh = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    T = d['press']*MUH*1e10/KB/(nh*(1.0 - x['H2'] + 0.1 + xe))
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
    ri = np.digitize(R.ravel(), rb)
    w = nh.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]

    def radial(f):
        return np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)]/np.maximum(den, 1e-300)

    return 0.5*(rb[1:] + rb[:-1]), radial(T), radial(x['H+'])


def radius_half(r, xhp):
    i = np.nonzero(xhp < 0.5)[0][0]
    return np.interp(0.5, [xhp[i], xhp[i - 1]], [r[i], r[i - 1]])


NBANDS = [3, 5, 6, 7]
if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    data = {}
    for nb in NBANDS:
        run = paths.runs('T6_multi_ion', 'F10_bands', 'run64_b%d' % nb)
        data['r_b%d' % nb], data['T_b%d' % nb], data['xhp_b%d' % nb] = tigris_profile(run)
        data['wall_b%d' % nb] = float(open(os.path.join(run, 'time.out')).read().split()[1])
    reduced.save(REDUCED, 'F10: density-weighted radial profiles of the Tigris runs '
                 'T6_multi_ion/F10_bands/run64_b<n>\nand their wall time\nr [pc], T [K], '
                 'xhp = x_H+, wall [s]', data)
else:
    data = reduced.load(REDUCED)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 7.5), sharex=True,
                             gridspec_kw={'height_ratios': [2, 1]})
for (lab, fn), col in zip(CLOUDY.items(), ['0.35', '0.65']):
    c = np.loadtxt(fn)
    rs = radius_half(c[:, 0], c[:, 3])
    a1.plot(c[:, 0], c[:, 1], color=col, lw=2.5, label='%s: $R_s$ = %.3f pc' % (lab, rs))
    a2.semilogy(c[:, 0], c[:, 3], color=col, lw=2.5)
# 6 bands add only the He II 54.4 eV edge, beyond which SB99 at 2 Myr emits almost no
# photons, so 5 and 6 bands coincide: 6 is drawn dashed with open markers
for nb, col, ls, mfc in zip(NBANDS, ['C3', 'C1', 'C2', 'C0'], ['-', '-', '--', '-'],
                            [None, None, 'none', None]):
    r, T, xhp = data['r_b%d' % nb], data['T_b%d' % nb], data['xhp_b%d' % nb]
    wall = float(data['wall_b%d' % nb])
    lab = 'Tigris %d bands: $R_s$ = %.3f pc, %.1f s wall' % (nb, radius_half(r, xhp), wall)
    print(lab, ' T(0.3, 1.0, 2.0, 2.7 pc) =', np.interp([0.3, 1.0, 2.0, 2.7], r, T).round(0))
    kw = dict(color=col, ls=ls, marker='o', mfc=mfc, ms=5 if mfc else 3.5, lw=1)
    a1.plot(r, T, label=lab, **kw)
    a2.semilogy(r, xhp, **kw)
a1.set(ylabel='T [K]', ylim=(4000, 11000),
       title=r'Static H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, SB99 2 Myr')
a1.legend(fontsize=8, loc='upper left')
a2.set(ylabel=r'$x_{\rm H^+}$', ylim=(1e-3, 1.5), xlabel='r [pc]', xlim=(0, 3.5))
for a in (a1, a2):
    a.grid(alpha=0.25)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
