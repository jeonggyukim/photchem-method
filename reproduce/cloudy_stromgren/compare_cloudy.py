"""Tigris static sphere (hydro fixed) against the Cloudy Stromgren reference.

usage: [CLOUDY_PROFILES=radial_profiles.txt] python compare_cloudy.py OUT.png \
           LABEL1:RUNDIR1 [LABEL2:RUNDIR2 ...]

Radial profiles (density-weighted spherical averages of the last snapshot) of
x_H+, x_He+, O+/O, O++/O, S+/S, S++/S, S3+/S and T, and the H+-zone averages
printed beside Cloudy's. A run whose input sets GOW17_temperature_fixed is drawn
at that T, since the fixed-T path leaves the pressure unchanged.
"""
import glob
import os
import re
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read()

CLOUDY = (os.environ.get('CLOUDY_PROFILES')
          or paths.runs('T6_multi_ion', 'cloudy_stromgren', 'radial_profiles.txt'))
MUH, KB = 2.34335276e-24, 1.380649e-16
P_TO_PK = MUH*1e10/KB
X_HE, X_O, X_S, X_N = 0.1, 3.2e-4, 1.45e-5, 7.4e-5
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}
TOTAL = {'H+': 1.0, 'He+': X_HE, 'O+': X_O, 'O++': X_O, 'S+': X_S, 'S++': X_S, 'S3+': X_S,
         'N+': X_N, 'N++': X_N}
QTY = ['H+', 'He+', 'O+', 'O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
# Cloudy radial_profiles.txt columns
CCOL = {'H+': 3, 'He+': 4, 'O+': 7, 'O++': 8, 'S+': 11, 'S++': 12, 'S3+': 13,
        'N+': 14, 'N++': 15}


def fixed_temperature(rundir):
    fn = rundir + '/athinput.runtime'
    if not os.path.exists(fn):
        return None
    m = re.search(r'^GOW17_temperature_fixed\s*=\s*([-+.\deE]+)', open(fn).read(), re.M)
    return float(m.group(1)) if m and float(m.group(1)) > 0 else None


def tigris_profiles(rundir):
    fn = sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1]
    d = athena_read.athdf(fn)
    nscal = sum(1 for k in d if k == 'rHI' or (k.startswith('r') and k[1:].isdigit()))
    names = CORE + IONS[:nscal - len(CORE)]
    keys = ['rHI'] + ['r%d' % n for n in range(1, nscal)]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(names, keys)}
    nH = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
    tfix = fixed_temperature(rundir)
    T = (np.full_like(nH, tfix) if tfix else
         d['press']*P_TO_PK/(nH*(1.0 - x['H2'] + X_HE + xe)))
    frac = {q: x[q]/TOTAL[q] for q in QTY if q in x}
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
    ri = np.digitize(R.ravel(), rb)
    w = nH.ravel()

    def radial(f):
        num = np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)]
        den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
        return np.where(den > 0, num/np.maximum(den, 1e-300), np.nan)

    r = 0.5*(rb[1:] + rb[:-1])
    prof = {k: radial(v) for k, v in frac.items()}
    prof['T'] = radial(T)
    # H+ zone: cells with x_H+ > 0.5, volume averages as Cloudy's table
    m = x['H+'] > 0.5
    avg = {k: np.mean(v[m]) for k, v in frac.items()}
    avg['T'] = np.mean(T[m])
    rs = (3*m.sum()*dx**3/(4*np.pi))**(1/3)
    return r, prof, avg, rs


c = np.loadtxt(CLOUDY)
cr, cT = c[:, 0], c[:, 1]
cfr = {q: c[:, CCOL[q]] for q in QTY}
mz = cfr['H+'] > 0.5
vol = np.gradient(cr)*cr**2
cavg = {k: np.sum(v[mz]*vol[mz])/np.sum(vol[mz]) for k, v in cfr.items()}
cavg['T'] = np.sum(cT[mz]*vol[mz])/np.sum(vol[mz])
crs = cr[mz].max()

out = sys.argv[1]
runs = [a.split(':', 1) for a in sys.argv[2:]]
fig, axs = plt.subplots(3, 4, figsize=(21, 13))
axs = axs.ravel()
for n, q in enumerate(QTY):
    axs[n].semilogy(cr, cfr[q], 'k-', lw=2, label='Cloudy 25')
    axs[n].set(title=q + (' / element' if q != 'H+' else ''), xlabel='r [pc]',
               ylim=(1e-4, 1.5), xlim=(0, 4))
iT = len(QTY)
axs[iT].plot(cr, cT, 'k-', lw=2, label='Cloudy 25')
axs[iT].set(title='T [K]', xlabel='r [pc]', xlim=(0, 4), ylim=(0, 1.5e4))

res = []
for irun, (lab, rd) in enumerate(runs):
    col = 'C%d' % irun
    r, prof, avg, rs = tigris_profiles(rd)
    res.append((avg, rs))
    for n, q in enumerate(QTY):
        if q in prof:
            axs[n].semilogy(r, prof[q], '-o', ms=3, color=col, alpha=0.6, label=lab)
    axs[iT].plot(r, prof['T'], '-o', ms=3, color=col, alpha=0.6, label=lab)

lines = ['%-9s %8s' % ('', 'Cloudy') + ''.join('%12s' % lab[:12] for lab, _ in runs)]
for q in QTY + ['T']:
    row = '%-9s %8.3g' % (q, cavg[q])
    for avg, _ in res:
        row += '%12.3g' % avg[q] if q in avg else '%12s' % '-'
    lines.append(row)
lines.append('%-9s %8.3f' % ('R_s [pc]', crs) + ''.join('%12.3f' % rs for _, rs in res))
print('\n'.join(lines))
axs[iT + 1].axis('off')
axs[iT + 1].text(0, 1, 'H+-zone volume averages (x_H+ > 0.5)\n\n' + '\n'.join(lines),
                 va='top', family='monospace', fontsize=9)
for ax in axs[:iT + 1]:
    ax.legend(fontsize=8)
fig.suptitle('Static Stromgren sphere, Q = 1e49 s^-1, n_H = 100 cm^-3, SB99 2 Myr: '
             'Tigris GOW17 vs Cloudy', fontsize=13)
fig.tight_layout()
fig.savefig(out, dpi=90)
print(out)
