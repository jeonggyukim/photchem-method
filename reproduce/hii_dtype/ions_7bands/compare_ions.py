"""hii_dtype with GOW17 + ions (O2,S3,N2), 7 bands (0-4 ionizing), against the
stored NCR reference and the GOW17 run without ions.

usage: python compare_ions.py RUNDIR OUT.png [M2_RUNDIR]
Prints f_ion over the ionizing bands, the per-band photon budget residual
|Ltot - Labs - Lesc|/Ltot, r_sh vs NCR, and H+-zone T and ion fractions at the
last snapshot.
"""
import glob
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
athena_read = paths.athena_read()

rundir, out = sys.argv[1], sys.argv[2]
m2 = sys.argv[3] if len(sys.argv) > 3 else None
ref = paths.tigris('tst/regression/data/ref_rayt_solutions/hii_dtype_ncr.txt')
NBAND, NION = 7, 5
MUH, KB = 2.34335276e-24, 1.380649e-16
P_TO_PK = MUH*1e10/KB
X_HE, X_O, X_S, X_N = 0.1, 3.2e-4, 1.45e-5, 7.4e-5
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}
TOTAL = {'He+': X_HE, 'O+': X_O, 'O++': X_O, 'S+': X_S, 'S++': X_S, 'S3+': X_S,
         'N+': X_N, 'N++': X_N}
trapz = getattr(np, 'trapezoid', None) or np.trapz


def shell(hst):
    h = athena_read.hst(hst)
    m = h['sh_mass'] > 0
    return h, h['time'][m], h['sh_mass_r'][m]/h['sh_mass'][m]


h, t, r = shell(rundir + '/HII.hst')
tt = h['time']
mL = h['Ltot0'] > 0
Ltot = sum(h['Ltot%d' % b] for b in range(NION))
Lgas = sum(h['Labs%d' % b] - h['Ldust%d' % b] for b in range(NION))
f_ion = trapz(Lgas[mL], tt[mL])/trapz(Ltot[mL], tt[mL])
print(f'f_ion (gas share of ionizing luminosity, bands 0-{NION - 1}) {f_ion:.3f}')
for b in range(NBAND):
    lt = h['Ltot%d' % b][mL]
    if np.all(lt == 0):
        print(f'band {b}: no source luminosity')
        continue
    res = np.abs(lt - h['Labs%d' % b][mL] - h['Lesc%d' % b][mL])/lt
    print(f'band {b}: budget residual max {res.max():.2e}; dust share '
          f'{trapz(h["Ldust%d" % b][mL], tt[mL])/trapz(lt, tt[mL]):.3f}')

dd = np.genfromtxt(ref, dtype=None, names=True, encoding=None)
tr, rr = dd['time'], dd['r_sh']
ok = (t >= tr[0]) & (t <= tr[-1])
diff = (r[ok] - np.interp(t[ok], tr, rr))/np.interp(t[ok], tr, rr)
print(f'r_sh at t = {t[-1]:.3f}: {r[-1]:.3f} pc; NCR reference {np.interp(t[-1], tr, rr):.3f} pc')
print(f'(r - r_NCR)/r_NCR: median {np.median(diff):+.4f} min {diff.min():+.4f} '
      f'max {diff.max():+.4f}')

fn = sorted(glob.glob(rundir + '/HII.out2.*.athdf'))[-1]
d = athena_read.athdf(fn)
nscal = sum(1 for k in d if k == 'rHI' or (k.startswith('r') and k[1:].isdigit()))
names = CORE + IONS[:nscal - len(CORE)]
keys = ['rHI'] + ['r%d' % n for n in range(1, nscal)]
x = {nm: d[k].astype(np.float64) for nm, k in zip(names, keys)}
nH = d['rho'].astype(np.float64)
xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
T = d['press']*P_TO_PK/(nH*(1.0 - x['H2'] + X_HE + xe))
mi = x['H+'] > 0.5
print(f'{fn.split("/")[-1]} t = {d["Time"]:.3f}: H+ zone {mi.sum()} cells, '
      f'T mass-weighted {np.sum(T[mi]*nH[mi])/np.sum(nH[mi]):.0f} K, '
      f'median {np.median(T[mi]):.0f} K')
print('H+-zone mass-weighted fractions: ' + ', '.join(
    f'{q} {np.sum(x[q][mi]*nH[mi])/np.sum(nH[mi])/TOTAL[q]:.3f}'
    for q in TOTAL if q in x))

fig, (ax, axd) = plt.subplots(2, 1, figsize=(6, 6.5), sharex=True,
                              gridspec_kw={'height_ratios': [2.2, 1]})
ax.plot(tr, rr, 'k--', lw=1.5, label='NCR (stored reference)')
if m2:
    _, t2, r2 = shell(m2 + '/HII.hst')
    ax.plot(t2, r2, 'C0-', lw=1.5, alpha=0.6, label='GOW17, 3 bands, no ions')
    ok2 = (t2 >= tr[0]) & (t2 <= tr[-1])
    axd.plot(t2[ok2], (r2[ok2] - np.interp(t2[ok2], tr, rr))/np.interp(t2[ok2], tr, rr),
             'C0-', alpha=0.6)
ax.plot(t, r, 'C3-', lw=2, alpha=0.6, label='GOW17 + O2,S3,N2, 7 bands')
ax.set_ylabel(r'$r_{\rm sh}$ [pc]')
ax.legend(fontsize=9, loc='lower right')
ax.set_title(r'hii_dtype, $64^3$, $Q_i = 10^{49}\,{\rm s^{-1}}$, $n_{\rm H} = 100\,{\rm cm^{-3}}$')
axd.plot(t[ok], diff, 'C3-', alpha=0.6)
axd.axhline(0, color='k', lw=0.6)
axd.set_ylabel(r'$(r - r_{\rm NCR})/r_{\rm NCR}$')
axd.set_xlabel('t [code time = 0.978 Myr]')
fig.tight_layout()
fig.savefig(out, dpi=150)
print(out)
