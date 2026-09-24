"""F12: D-type H II region, Q = 1e49 s^-1 into n_H = 100 cm^-3, 64^3. Top: shell radius
against time for Simple (T fixed at 8000 K) and NCR (stored regression references) and
GOW17 + O3,S3,N3 (this run), with the Spitzer (1978) and Hosokawa & Inutsuka (2006)
solutions at 8000 K; mass-weighted T of the ionized gas (x_H+ > 0.9) against time; radial
profiles at t = 1 code. Bottom: z = 0 slices of n_H at four times, one colour range.
Writes ../../../figures/F12_hii_dtype.{pdf,png}.

Inputs: run_f12.sh output in WORKDIR/run; tigris-gow17 tst/regression/data/ref_rayt_solutions."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F12_hii_dtype')
RUN = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion/M5_hii_dtype_ions/'
                         'F12_run/run')
REF = '/Users/jgkim/Projects/tigris-gow17/tst/regression/data/ref_rayt_solutions/'
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2, 'S_high': 3,
          'N+': 1, 'N++': 2, 'N_high': 3}
XTOT = {'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
MUH, KB, PC, KMS = 1.4*1.6738234e-24, 1.38065e-16, 3.08567758e18, 1e5
TUNIT = PC/KMS


def state(fn):
    d = athena_read.athdf(fn)
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    rho = d['rho'].astype(np.float64)  # m_H cm^-3, so n_H = rho/1.4
    nh = rho/1.4
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    T = d['press']*MUH*1e10/KB/(rho*(1.0 - x['H2'] + 0.1 + xe))
    return d, x, nh, T


def radius_shell(t, tion, qi=1e49, n=100.0):
    cion = np.sqrt(2.1*KB*tion/MUH)
    alpha_b = 2.59e-13*(tion*1e-4)**-0.7
    r0 = (3.0*qi/(4.0*np.pi*alpha_b*n**2))**(1/3)
    s78 = r0*(1.0 + 7/4*cion*t/r0)**(4/7)
    h06 = r0*(1.0 + 7/4*np.sqrt(4/3)*cion*t/r0)**(4/7)
    return s78/PC, h06/PC


fig = plt.figure(figsize=(17, 9.5), layout='constrained')
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05])
ax_r, ax_t = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2])
ax_p = fig.add_subplot(gs[0, 3])

# shell radius
h = athena_read.hst(os.path.join(RUN, 'hii_dtype_gow17_ions.hst'))
m = h['sh_mass'] > 0
t_g, r_g = h['time'][m], (h['sh_mass_r']/h['sh_mass'])[m]
tt = np.linspace(0, 1.0, 200)
s78, h06 = radius_shell(tt*TUNIT, 8000.0)
ax_r.plot(tt, s78, 'k:', lw=1.2, label='Spitzer (1978), 8000 K')
ax_r.plot(tt, h06, 'k--', lw=1.2, label='Hosokawa & Inutsuka (2006), 8000 K')
for fn, lab, col in [('hii_dtype_simple.txt', 'Simple (T = 8000 K)', 'C7'),
                     ('hii_dtype_ncr.txt', 'NCR', 'C1')]:
    dd = np.genfromtxt(REF + fn, names=True)
    ax_r.plot(dd['time'], dd['r_sh'], color=col, lw=2, label=lab)
ax_r.plot(t_g, r_g, color='C0', lw=2, label='GOW17 + O3,S3,N3')
ax_r.set(xlabel='t [code = 0.978 Myr]', ylabel=r'$r_{\rm sh}$ [pc]', xlim=(0, 1),
         ylim=(0, 10))
ax_r.legend(fontsize=9, loc='lower right')
ax_r.grid(alpha=0.25)

# ionized-gas temperature
dumps = sorted(glob.glob(os.path.join(RUN, '*.out2.*.athdf')))
times, tions = [], []
for fn in dumps[1:]:
    d, x, nh, T = state(fn)
    sel = x['H+'] > 0.9
    if sel.any():
        times.append(float(d['Time']))
        tions.append(np.sum(T[sel]*nh[sel])/np.sum(nh[sel]))
ax_t.plot(times, tions, 'C0o-', ms=3)
ax_t.set(xlabel='t [code]', ylabel=r'mass-weighted $T(x_{\rm H^+} > 0.9)$ [K]',
         ylim=(6000, 10000), xlim=(0, 1))
ax_t.grid(alpha=0.25)
print('T_ion at t = 1: %.0f K; r_sh(t = 1) %.3f pc' % (tions[-1], r_g[-1]))

# radial profiles at t = 1
d, x, nh, T = state(dumps[-1])
X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
R = np.sqrt(X1**2 + X2**2 + X3**2)
dx = d['x1f'][1] - d['x1f'][0]
rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
ri = np.digitize(R.ravel(), rb)
w = nh.ravel()
den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
rc = 0.5*(rb[1:] + rb[:-1])


def radial(f):
    return np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)]/np.maximum(den, 1e-300)


for k, el, col, lab in [('H+', None, 'k', 'H$^+$'), ('O+', 'O', 'C1', 'O$^+$'),
                        ('O++', 'O', 'C2', 'O$^{2+}$'), ('S+', 'S', 'C4', 'S$^+$'),
                        ('S++', 'S', 'C5', 'S$^{2+}$'), ('N++', 'N', 'C6', 'N$^{2+}$')]:
    ax_p.semilogy(rc, radial(x[k]/(XTOT[el] if el else 1.0)), color=col, lw=1.5, label=lab)
ax_p.set(xlabel='r [pc]', ylabel='fraction of the element', ylim=(1e-3, 1.5), xlim=(0, 12),
         title='t = 1 code')
ax_p.legend(fontsize=8, loc='lower left', ncol=2)
axT = ax_p.twinx()
axT.plot(rc, radial(T), 'r:', lw=1.5)
axT.set_ylabel('T [K] (dotted)', color='r')
axT.set_yscale('log')
axT.set_ylim(30, 3e4)

# slices
for j, tsel in enumerate([0.1, 0.3, 0.6, 1.0]):
    fn = min(dumps, key=lambda f: abs(athena_read.athdf(f, quantities=['rho'])['Time']
                                      - tsel))
    d, x, nh, T = state(fn)
    k = T.shape[0]//2
    ax = fig.add_subplot(gs[1, j])
    im = ax.pcolormesh(d['x1f'], d['x2f'], nh[k], norm=LogNorm(10, 1e3), cmap='viridis',
                       shading='flat')
    ax.contour(d['x1v'], d['x2v'], x['H+'][k], levels=[0.5], colors='c', linewidths=0.8)
    ax.set_aspect('equal')
    ax.set(title=r'$n_{\rm H}$, z = 0, t = %.2f code' % d['Time'], xlabel='x [pc]')
    if j == 0:
        ax.set_ylabel('y [pc]')
fig.colorbar(im, ax=fig.axes[-4:], shrink=0.9, label=r'$n_{\rm H}$ [cm$^{-3}$] (cyan: $x_{\rm H^+} = 0.5$)')
fig.suptitle(r'D-type H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, '
             r'$64^3$ (0.56 pc cells)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
