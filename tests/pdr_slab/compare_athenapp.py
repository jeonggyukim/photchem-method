"""GOW17 pdr_slab in Tigris against Athena++ chem_pdr_static.vtk (GOW17, CVODE,
six-ray with one lit ray, n_H = 100 cm^-3, chi = 1 at the face): species and
temperature against column."""
import glob, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read('ATHENAPP_PDR_DIR')
run, out = sys.argv[1], sys.argv[2]
pc, nh0, kB, e_unit = 3.0856776e18, 100.0, 1.380649e-16, 1.6738234e-24*1e10
names = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
charge = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1}

x, _, _, d = athena_read.vtk(paths.athenapp_pdr('tst/regression/data/chem_pdr_static.vtk'))
xc = 0.5*(x[1:] + x[:-1])
ap = {n: d['r' + n].ravel() for n in names}
ap['e'] = sum(charge.get(n, 0)*ap[n] for n in names)
# Athena++ chemistry units: [rho] = 1.4 m_H cm^-3 (src/units/units.cpp)
ap['T'] = d['press'].ravel()*1.4*e_unit/kB/(nh0*(1.1 + ap['e'] - ap['H2']))
ap['N'] = xc*pc*nh0

rows, head = [], None
for f in sorted(glob.glob(f'{run}.block*.out3.00001.tab')):
    for l in open(f):
        if l.startswith('# i'): head = l[1:].split()
        elif not l.startswith('#'): rows.append([float(v) for v in l.split()])
a = np.array(rows); a = a[np.argsort(a[:, 1])]
col = {h: a[:, n] for n, h in enumerate(head)}
tg = {n: col['rHI' if k == 0 else f'r{k}'] for k, n in enumerate(names)}
tg['e'] = sum(charge.get(n, 0)*tg[n] for n in names)
tg['T'] = col['press']*e_unit/kB/(nh0*(1.1 + tg['e'] - tg['H2']))
tg['N'] = col['x1v']*pc*nh0

def at(prof, N):
    return np.exp(np.interp(np.log(N), np.log(prof['N']), np.log(np.maximum(prof[q], 1e-30))))
print(f"{'N [cm^-2]':>10s} " + ' '.join(f'{q:>16s}' for q in ('T', 'H2', 'C+', 'CO', 'e')))
for N in (1e19, 1e20, 3e20, 1e21, 2e21, 3e21, 5e21, 1e22):
    s = f'{N:10.1e} '
    for q in ('T', 'H2', 'C+', 'CO', 'e'):
        s += f' {at(tg, N):8.3g}/{at(ap, N):<7.3g}'
    print(s)
print('(Tigris / Athena++)')
def front(p, y, level):
    i = np.where(y >= level)[0][0]; return p['N'][i]
for lab, q, lev in (('2x_H2 = 0.5', 'H2', 0.25), ('x_CO = xC/2', 'CO', 0.8e-4)):
    ft, fa = front(tg, tg[q], lev), front(ap, ap[q], lev)
    print(f'{lab}: Tigris {ft:.3e}, Athena++ {fa:.3e}, relative {ft/fa - 1:+.3f}')

fig, ax = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True)
for q, c in (('H2', 'C0'), ('C+', 'C1'), ('CO', 'C2'), ('e', 'C3'), ('H+', 'C4')):
    ax[0].loglog(tg['N'], tg[q], c=c, label=q)
    ax[0].loglog(ap['N'], ap[q], c=c, ls='--')
ax[0].set_ylim(1e-9, 1); ax[0].set_ylabel('x per H')
ax[0].legend(title='solid Tigris GOW17, dashed Athena++', fontsize=8, ncol=2)
ax[1].semilogx(tg['N'], tg['T'], 'k-', label='Tigris GOW17 (semi-implicit, postproc)')
ax[1].semilogx(ap['N'], ap['T'], 'k--', label='Athena++ GOW17 (CVODE)')
ax[1].set_ylim(0, 100); ax[1].set_ylabel('T [K]'); ax[1].set_xlabel(r'$N_{\rm H}$ [cm$^{-2}$]')
ax[1].legend(fontsize=8); ax[1].set_xlim(1e18, 1.3e22)
fig.suptitle(r'PDR slab, $n_{\rm H} = 100\,{\rm cm^{-3}}$, $\chi = 1$ at the face')
fig.tight_layout(); fig.savefig(out, dpi=150)
