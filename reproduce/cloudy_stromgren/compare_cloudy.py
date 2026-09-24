"""Tigris static sphere (hydro fixed) against the Cloudy Stromgren reference.

usage: python compare_cloudy.py OUT.png LABEL1:RUNDIR1 [LABEL2:RUNDIR2 ...]

Radial profiles (density-weighted spherical averages of the last snapshot) of
x_H+, x_He+, O+/O, O++/O, S+/S, S++/S and T, and the H+-zone averages printed
beside Cloudy's.
"""
import glob
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

import os
CLOUDY = os.environ.get('CLOUDY_PROFILES',
                        '/Users/jgkim/Documents/tigris-photchem-gow17-multi-ion/'
                        'T6_multi_ion/cloudy_stromgren/radial_profiles.txt')
MUH, KB = 2.34335276e-24, 1.380649e-16
P_TO_PK = MUH*1e10/KB
X_HE, X_O, X_S = 0.1, 3.2e-4, 1.45e-5
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2}
QTY = ['H+', 'He+', 'O+', 'O++', 'S+', 'S++']


def tigris_profiles(rundir):
    fn = sorted(glob.glob(rundir + '/HII.out2.*.athdf'))[-1]
    d = athena_read.athdf(fn)
    nscal = sum(1 for k in d if k == 'rHI' or (k.startswith('r') and k[1:].isdigit()))
    names = CORE + (IONS if nscal > len(CORE) else [])
    keys = ['rHI'] + ['r%d' % n for n in range(1, nscal)]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(names, keys)}
    nH = d['rho'].astype(np.float64)
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
    T = d['press']*P_TO_PK/(nH*(1.0 - x['H2'] + X_HE + xe))
    frac = {'H+': x['H+'], 'He+': x['He+']/X_HE, 'O+': x['O+']/X_O}
    if 'O++' in x:
        frac.update({'O++': x['O++']/X_O, 'S+': x['S+']/X_S, 'S++': x['S++']/X_S})
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
    return r, prof, avg, rs, d['Time']


c = np.loadtxt(CLOUDY)
cr, cT = c[:, 0], c[:, 1]
cfr = {'H+': c[:, 3], 'He+': c[:, 4], 'O+': c[:, 7], 'O++': c[:, 8], 'S+': c[:, 11],
       'S++': c[:, 12]}
mz = cfr['H+'] > 0.5
vol = np.gradient(cr)*cr**2
cavg = {k: np.sum(v[mz]*vol[mz])/np.sum(vol[mz]) for k, v in cfr.items()}
cavg['T'] = np.sum(cT[mz]*vol[mz])/np.sum(vol[mz])
crs = cr[mz].max()

out = sys.argv[1]
runs = [a.split(':', 1) for a in sys.argv[2:]]
fig, axs = plt.subplots(2, 4, figsize=(20, 9))
axs = axs.ravel()
for n, q in enumerate(QTY):
    axs[n].semilogy(cr, cfr[q], 'k-', lw=2, label='Cloudy 25')
    axs[n].set(title=q + (' / element' if q != 'H+' else ''), xlabel='r [pc]',
               ylim=(1e-4, 1.5), xlim=(0, 4))
axs[6].plot(cr, cT, 'k-', lw=2, label='Cloudy 25')
axs[6].set(title='T [K]', xlabel='r [pc]', xlim=(0, 4), ylim=(0, 1.5e4))

print('%-10s %8s' % ('', 'Cloudy') + ''.join('%10s' % lab for lab, _ in runs))
res = []
for lab, rd in runs:
    r, prof, avg, rs, t = tigris_profiles(rd)
    res.append((avg, rs))
    for n, q in enumerate(QTY):
        if q in prof:
            axs[n].semilogy(r, prof[q], '-o', ms=3, label=lab)
    axs[6].plot(r, prof['T'], '-o', ms=3, label=lab)
for q in QTY + ['T']:
    row = '%-10s %8.3g' % (q, cavg[q])
    for avg, _ in res:
        row += '%10.3g' % avg[q] if q in avg else '%10s' % '-'
    print(row)
print('%-10s %8.3f' % ('R_s [pc]', crs) + ''.join('%10.3f' % rs for _, rs in res))
axs[7].axis('off')
axs[7].text(0, 1, 'H+-zone volume averages\n(x_H+ > 0.5)\n\n' + '\n'.join(
    '%-4s Cloudy %.3g | ' % (q, cavg[q]) + ' | '.join('%s %.3g' % (lab, a[q])
                                                    for (lab, _), (a, _) in zip(runs, res)
                                                    if q in a)
    for q in QTY + ['T']), va='top', family='monospace', fontsize=9)
for ax in axs[:7]:
    ax.legend(fontsize=8)
fig.suptitle('Static Stromgren sphere, Q = 1e49 s^-1, n_H = 100 cm^-3, SB99 2 Myr: '
             'Tigris GOW17 vs Cloudy', fontsize=13)
fig.tight_layout()
fig.savefig(out, dpi=90)
print(out)
