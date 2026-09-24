"""Shell radius of hii_dtype with GOW17 against the stored NCR reference and the
Spitzer / Hosokawa-Inutsuka curves, as tst/.../rayt_point/hii_dtype.py draws them."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

hst, out = sys.argv[1], sys.argv[2]
ref = '/Users/jgkim/Projects/tigris-gow17/tst/regression/data/ref_rayt_solutions/hii_dtype_ncr.txt'
muH, pc, kB, Qi, nH0 = 1.4*1.6738234e-24, 3.08567758e18, 1.38065e-16, 1e49, 1e2
tunit = pc/1e5

def analytic(t, tgas, f_ion):
    cion = np.sqrt(2.1*kB*tgas/muH)
    alphaB = 2.59e-13*(tgas*1e-4)**-0.7
    R0 = (3.0*f_ion*Qi/(4.0*np.pi*alphaB*nH0**2))**(1/3)
    t = t*tunit
    return (R0*(1 + 7/4*cion*t/R0)**(4/7)/pc,
            R0*(1 + 7/4*np.sqrt(4/3)*cion*t/R0)**(4/7)/pc)

h = athena_read.hst(hst)
m = h['sh_mass'] > 0
t, r = h['time'][m], h['sh_mass_r'][m]/h['sh_mass'][m]
trapz = getattr(np, 'trapezoid', None) or np.trapz
mL = h['Ltot0'] > 0
f_ion = trapz((h['Labs0'] - h['Ldust0'])[mL], h['time'][mL])/trapz(h['Ltot0'][mL], h['time'][mL])
dd = np.genfromtxt(ref, dtype=None, names=True, encoding=None)
tr, rr = dd['time'], dd['r_sh']
ok = (t >= tr[0]) & (t <= tr[-1])
diff = (r[ok] - np.interp(t[ok], tr, rr))/np.interp(t[ok], tr, rr)
print(f'f_ion (gas-absorbed share of ionizing luminosity) {f_ion:.3f}')
print(f'r_sh at t = {t[-1]:.3f}: {r[-1]:.3f} pc; NCR reference {np.interp(t[-1], tr, rr):.3f} pc')
print(f'(r - r_NCR)/r_NCR over t in [{t[ok][0]:.3f}, {t[ok][-1]:.3f}]: '
      f'median {np.median(diff):+.4f} min {diff.min():+.4f} max {diff.max():+.4f}')
for k in ('Lesc0', 'Ldust0', 'Labs0'):
    if k in h.dtype.names if hasattr(h, 'dtype') else k in h:
        print(k, 'share', trapz(h[k][mL], h['time'][mL])/trapz(h['Ltot0'][mL], h['time'][mL]))

fig, (ax, axd) = plt.subplots(2, 1, figsize=(6, 6.5), sharex=True,
                              gridspec_kw={'height_ratios': [2.2, 1]})
ax.plot(t, r, 'C3-', lw=2, label='GOW17')
ax.plot(tr, rr, 'k--', lw=1.5, label='NCR (stored reference)')
s78, h06 = analytic(t, 7200.0, f_ion)
ax.plot(t, s78, 'C0:', label=f'Spitzer 78 (T = 7200 K, f_ion = {f_ion:.2f})')
ax.plot(t, h06, 'C2:', label=f'Hosokawa 06 (T = 7200 K, f_ion = {f_ion:.2f})')
ax.set_ylabel(r'$r_{\rm sh}$ [pc]'); ax.legend(fontsize=9, loc='lower right')
ax.set_title(r'hii_dtype, $64^3$, $Q_i = 10^{49}\,{\rm s^{-1}}$, $n_{\rm H} = 100\,{\rm cm^{-3}}$')
axd.plot(t[ok], diff, 'C3-'); axd.axhline(0, color='k', lw=0.6)
axd.set_ylabel(r'$(r - r_{\rm NCR})/r_{\rm NCR}$'); axd.set_xlabel('t [code time = 0.978 Myr]')
fig.tight_layout(); fig.savefig(out, dpi=150)
