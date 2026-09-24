"""Radial profiles of hii_dtype with GOW17 + O2,S3,N2: n_H, T, x_H+ and the ion
fractions per element, at several snapshots.

usage: python ion_profiles.py RUNDIR OUT.png [snapshot numbers, default 10 25 50]
Density-weighted spherical averages about the source at the origin.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

rundir, out = sys.argv[1], sys.argv[2]
snaps = [int(s) for s in sys.argv[3:]] or [10, 25, 50]
# ism units: density unit m_H/cm^3 and rho = 1.4 n_H; pressure unit m_H (km/s)^2
MH, KB, MU_H = 1.6735575e-24, 1.380649e-16, 1.4
P_TO_PK = MH*1e10/KB
X_HE, X_O, X_S, X_N = 0.1, 3.2e-4, 1.45e-5, 7.4e-5
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}
PANELS = [('H', {'H+': 1.0}), ('He', {'He+': X_HE}),
          ('O', {'O+': X_O, 'O++': X_O}), ('S', {'S+': X_S, 'S++': X_S, 'S3+': X_S}),
          ('N', {'N+': X_N, 'N++': X_N})]


def profiles(fn):
    d = athena_read.athdf(fn)
    nscal = sum(1 for k in d if k == 'rHI' or (k.startswith('r') and k[1:].isdigit()))
    names = CORE + IONS[:nscal - len(CORE)]
    keys = ['rHI'] + ['r%d' % n for n in range(1, nscal)]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(names, keys)}
    nH = d['rho'].astype(np.float64)/MU_H
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
    T = d['press']*P_TO_PK/(nH*(1.0 - x['H2'] + X_HE + xe))
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
    ri = np.digitize(R.ravel(), rb)
    w = nH.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    cnt = np.bincount(ri, None, len(rb) + 1)[1:len(rb)]

    def radial(f):
        num = np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)]
        return np.where(den > 0, num/np.maximum(den, 1e-300), np.nan)

    prof = {nm: radial(v) for nm, v in x.items()}
    prof['T'] = radial(T)
    prof['nH'] = np.where(cnt > 0, den/np.maximum(cnt, 1), np.nan)
    return 0.5*(rb[1:] + rb[:-1]), prof, d['Time']


fig, axs = plt.subplots(2, 4, figsize=(19, 8.5))
axs = axs.ravel()
for i, s in enumerate(snaps):
    r, p, t = profiles('%s/HII.out2.%05d.athdf' % (rundir, s))
    ls = ['-', '--', ':'][i % 3]
    lab = 't = %.2f' % t
    axs[0].semilogy(r, p['nH'], 'k', ls=ls, label=lab)
    axs[1].plot(r, p['T'], 'k', ls=ls, label=lab)
    for n, (el, ions) in enumerate(PANELS):
        ax = axs[2 + n]
        for j, (q, tot) in enumerate(ions.items()):
            ax.semilogy(r, p[q]/tot, color='C%d' % j, ls=ls, alpha=0.7,
                        label=q if i == 0 else None)
axs[0].set(title=r'$n_{\rm H}$ [cm$^{-3}$]', xlabel='r [pc]', xlim=(0, 16))
axs[1].set(title='T [K]', xlabel='r [pc]', xlim=(0, 16), ylim=(0, 1.2e4))
for n, (el, _) in enumerate(PANELS):
    axs[2 + n].set(title=el + ' ion fractions', xlabel='r [pc]', xlim=(0, 16),
                   ylim=(1e-4, 1.5))
for ax in axs[:7]:
    ax.legend(fontsize=8)
axs[7].axis('off')
axs[7].text(0, 1, 'line style = snapshot time (code units, 0.978 Myr)\n'
            'solid, dashed, dotted in order of time', va='top', fontsize=10)
fig.suptitle(r'hii_dtype, GOW17 + O2,S3,N2, 7 bands, $64^3$: density-weighted '
             'spherical averages', fontsize=13)
fig.tight_layout()
fig.savefig(out, dpi=110)
print(out)
