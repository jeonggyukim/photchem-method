"""F11: where the interior temperature of the static H II region comes from. T(r) and
x(H+)(r) of the Tigris post-processed sphere (GOW17 + O3,S3,N3, 7 bands, case B H and
on-the-spot He) against Cloudy 25 with its default treatment of the diffuse field and with
`diffuse OTS`, and the ratio of Tigris to each. Writes
../../../figures/F11_interior_temperature.{pdf,png}.

Inputs: Tigris dump TIGRIS (postproc7/run64_high, 64^3, the last out2); Cloudy
radial_profiles.txt of reproduce/cloudy_stromgren (default) and
reproduce/cloudy_stromgren/thermal_with_N_ots."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.join(HERE, '..', '..')
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F11_interior_temperature')
TIGRIS = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion/T6_multi_ion/'
                            'postproc7/run64_high')
CLOUDY = {'Cloudy 25, default': os.path.join(REPRO, 'cloudy_stromgren', 'radial_profiles.txt'),
          'Cloudy 25, diffuse OTS': os.path.join(REPRO, 'cloudy_stromgren',
                                                 'thermal_with_N_ots', 'radial_profiles.txt')}
# species order of the dump: GOW17 core then the ion ladders with higher ions
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
# ism units in this run: density unit m_H/cm^3 with n_H = rho (photchem_postproc sets
# rho per H), pressure unit m_H (km/s)^2 per mu_H, as compare_cloudy.py
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


r, T, xhp = tigris_profile(TIGRIS)
fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7, 9), sharex=True,
                                 gridspec_kw={'height_ratios': [2, 1, 1]})
a1.plot(r, T, 'ko', ms=4, label='Tigris GOW17 + ions (case B H, on-the-spot He)')
a2.semilogy(r, xhp, 'ko', ms=4)
for (lab, fn), col in zip(CLOUDY.items(), ['C5', 'C0']):
    c = np.loadtxt(fn)
    a1.plot(c[:, 0], c[:, 1], color=col, lw=2, label=lab)
    a2.semilogy(c[:, 0], c[:, 3], color=col, lw=2)
    # inside the ionization front only: across it T drops by 10^3 over 0.05 pc
    m = (r > 0.1) & (r < 2.8)
    a3.plot(r[m], T[m]/np.interp(r[m], c[:, 0], c[:, 1]), 'o-', color=col, ms=3,
            label='Tigris / ' + lab)
a1.set(ylabel='T [K]', ylim=(4000, 11000),
       title=r'Static H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, SB99 2 Myr')
a1.legend(fontsize=8, loc='upper left')
a2.set(ylabel=r'$x_{\rm H^+}$', ylim=(1e-3, 1.5))
a3.axhline(1, color='0.5', lw=0.8)
a3.set(ylabel='T ratio', xlabel='r [pc]', ylim=(0.7, 1.2), xlim=(0, 3.5))
a3.legend(fontsize=8, loc='lower left')
for a in (a1, a2, a3):
    a.grid(alpha=0.25)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
