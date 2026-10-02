"""F12: D-type H II region, Q = 1e49 s^-1 into n_H = 100 cm^-3, 64^3. Top: shell radius
against time for Simple (T fixed at 8000 K) and NCR (stored regression references) and
GOW17 + O3,S3,N3 (this run), with the Spitzer (1978) and Hosokawa & Inutsuka (2006)
solutions at 8000 K; mass-weighted T of the ionized gas (x_H+ > 0.9) against time; radial
profiles at t = 1 code. Bottom: z = 0 slices of n_H at four times, one colour range.
Writes ../../../figures/F12_hii_dtype.{pdf,png}.

Inputs: run_f12.sh output in WORKDIR/run; tigris-gow17 tst/regression/data/ref_rayt_solutions.
The plotted arrays come from reduced.txt when it exists; `python plot.py --from-runs`
reads the run ($PHOTCHEM_RUNS/M5_hii_dtype_ions/F12_run/run) and the references
($TIGRIS_DIR) and rewrites it."""
import glob
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
import reduced  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F12_hii_dtype')
REDUCED = os.path.join(HERE, 'reduced.txt')
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


PROFILES = [('H+', None, 'k', 'H$^+$'), ('O+', 'O', 'C1', 'O$^+$'),
            ('O++', 'O', 'C2', 'O$^{2+}$'), ('S+', 'S', 'C4', 'S$^+$'),
            ('S++', 'S', 'C5', 'S$^{2+}$'), ('N++', 'N', 'C6', 'N$^{2+}$')]
TSLICE = [0.1, 0.3, 0.6, 1.0]
if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    run = paths.runs('M5_hii_dtype_ions', 'F12_run', 'run')
    data = {}
    # shell radius
    h = athena_read.hst(os.path.join(run, 'hii_dtype_gow17_ions.hst'))
    m = h['sh_mass'] > 0
    data['t_g'], data['r_g'] = h['time'][m], (h['sh_mass_r']/h['sh_mass'])[m]
    for name in ('simple', 'ncr'):
        dd = np.genfromtxt(paths.tigris('tst/regression/data/ref_rayt_solutions',
                                        'hii_dtype_%s.txt' % name), names=True)
        data['time_' + name], data['r_sh_' + name] = dd['time'], dd['r_sh']
    # ionized-gas temperature
    dumps = sorted(glob.glob(os.path.join(run, '*.out2.*.athdf')))
    times, tions = [], []
    for fn in dumps[1:]:
        d, x, nh, T = state(fn)
        sel = x['H+'] > 0.9
        if sel.any():
            times.append(float(d['Time']))
            tions.append(np.sum(T[sel]*nh[sel])/np.sum(nh[sel]))
    data['times'], data['tions'] = times, tions
    # radial profiles at t = 1
    d, x, nh, T = state(dumps[-1])
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, dx)
    ri = np.digitize(R.ravel(), rb)
    w = nh.ravel()
    den = np.bincount(ri, w, len(rb) + 1)[1:len(rb)]
    data['rc'] = 0.5*(rb[1:] + rb[:-1])

    def radial(f):
        return (np.bincount(ri, f.ravel()*w, len(rb) + 1)[1:len(rb)]
                / np.maximum(den, 1e-300))

    for k, el, _, _ in PROFILES:
        data['prof_' + k] = radial(x[k]/(XTOT[el] if el else 1.0))
    data['prof_T'] = radial(T)
    # slices
    for j, tsel in enumerate(TSLICE):
        fn = min(dumps, key=lambda f: abs(athena_read.athdf(f, quantities=['rho'])['Time']
                                          - tsel))
        d, x, nh, T = state(fn)
        k = T.shape[0]//2
        data['time_%d' % j] = d['Time']
        data['nh_%d' % j], data['xhp_%d' % j] = nh[k], x['H+'][k]
    for c in ('x1f', 'x2f', 'x1v', 'x2v'):
        data[c] = d[c]
    reduced.save(REDUCED, 'F12: shell radius (t_g, r_g) and ionized-gas T (times, tions) '
                 'of M5_hii_dtype_ions/F12_run/run,\nthe Simple and NCR references of '
                 'tst/regression/data/ref_rayt_solutions,\nradial profiles at t = 1 '
                 '(rc, prof_*) and z = 0 slices of n_H and x_H+ at time_0 .. time_3\n'
                 't [code], r [pc], T [K], n_H [cm^-3]', data)
else:
    data = reduced.load(REDUCED)

fig = plt.figure(figsize=(17, 9.5), layout='constrained')
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05])
ax_r, ax_t = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2])
ax_p = fig.add_subplot(gs[0, 3])

# shell radius
tt = np.linspace(0, 1.0, 200)
s78, h06 = radius_shell(tt*TUNIT, 8000.0)
ax_r.plot(tt, s78, 'k:', lw=1.2, label='Spitzer (1978), 8000 K')
ax_r.plot(tt, h06, 'k--', lw=1.2, label='Hosokawa & Inutsuka (2006), 8000 K')
for name, lab, col in [('simple', 'Simple (T = 8000 K)', 'C7'), ('ncr', 'NCR', 'C1')]:
    ax_r.plot(data['time_' + name], data['r_sh_' + name], color=col, lw=2, label=lab)
ax_r.plot(data['t_g'], data['r_g'], color='C0', lw=2, label='GOW17 + O3,S3,N3')
ax_r.set(xlabel='t [code = 0.978 Myr]', ylabel=r'$r_{\rm sh}$ [pc]', xlim=(0, 1),
         ylim=(0, 10))
ax_r.legend(fontsize=9, loc='lower right')
ax_r.grid(alpha=0.25)

# ionized-gas temperature
ax_t.plot(data['times'], data['tions'], 'C0o-', ms=3)
ax_t.set(xlabel='t [code]', ylabel=r'mass-weighted $T(x_{\rm H^+} > 0.9)$ [K]',
         ylim=(6000, 10000), xlim=(0, 1))
ax_t.grid(alpha=0.25)
print('T_ion at t = 1: %.0f K; r_sh(t = 1) %.3f pc'
      % (data['tions'][-1], data['r_g'][-1]))

# radial profiles at t = 1
for k, el, col, lab in PROFILES:
    ax_p.semilogy(data['rc'], data['prof_' + k], color=col, lw=1.5, label=lab)
ax_p.set(xlabel='r [pc]', ylabel='fraction of the element', ylim=(1e-3, 1.5), xlim=(0, 12),
         title='t = 1 code')
ax_p.legend(fontsize=8, loc='lower left', ncol=2)
axT = ax_p.twinx()
axT.plot(data['rc'], data['prof_T'], 'r:', lw=1.5)
axT.set_ylabel('T [K] (dotted)', color='r')
axT.set_yscale('log')
axT.set_ylim(30, 3e4)

# slices
for j in range(len(TSLICE)):
    ax = fig.add_subplot(gs[1, j])
    im = ax.pcolormesh(data['x1f'], data['x2f'], data['nh_%d' % j], norm=LogNorm(10, 1e3),
                       cmap='viridis', shading='flat')
    ax.contour(data['x1v'], data['x2v'], data['xhp_%d' % j], levels=[0.5], colors='c',
               linewidths=0.8)
    ax.set_aspect('equal')
    ax.set(title=r'$n_{\rm H}$, z = 0, t = %.2f code' % data['time_%d' % j],
           xlabel='x [pc]')
    if j == 0:
        ax.set_ylabel('y [pc]')
fig.colorbar(im, ax=fig.axes[-4:], shrink=0.9, label=r'$n_{\rm H}$ [cm$^{-3}$] (cyan: $x_{\rm H^+} = 0.5$)')
fig.suptitle(r'D-type H II region, $Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, '
             r'$64^3$ (0.56 pc cells)', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
