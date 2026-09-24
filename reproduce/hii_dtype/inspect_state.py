"""Radial profile of T, x_H+, x_H2 from the per-block vtk dumps of a hii run."""
import glob, sys
import numpy as np
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read
athena_read.check_nan_flag = False
run, idx = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(f'{run}.block*.out2.{idx}.vtk'))
r_all, d_all, T_all, xp_all, xh2_all = [], [], [], [], []
names = None
for f in files:
    xf, yf, zf, dat = athena_read.vtk(f)
    if names is None:
        names = list(dat.keys()); print('fields', names)
    xc, yc, zc = [0.5*(a[1:] + a[:-1]) for a in (xf, yf, zf)]
    Z, Y, X = np.meshgrid(zc, yc, xc, indexing='ij')
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho, P = dat['rho'], dat['press']
    s = {k: dat[k]/1.0 for k in names if k.startswith('r')}
    r_all.append(r.ravel()); d_all.append(rho.ravel())
    # GOW17 scalar order: He+ OHx CHx CO C+ HCO+ H2 H+ H3+ H2+ O+ Si+
    sc = [dat['rHI']] + [dat[f'r{n}'] for n in range(1, 12)]  # r0 is written as rHI
    xH2, xHp = sc[6], sc[7]
    xe = sc[0] + sc[4] + sc[5] + sc[7] + sc[8] + sc[9] + sc[10] + sc[11]
    nH = rho/1.4
    T = P*1.6738234e-14/1.38065e-16/(nH*(1.1 + xe - xH2))
    T_all.append(T.ravel()); xp_all.append(xHp.ravel()); xh2_all.append(xH2.ravel())
r, d, T, xp, xh2 = map(np.concatenate, (r_all, d_all, T_all, xp_all, xh2_all))
bins = np.arange(0, 12, 1.0)
print(' r[pc]   n_H     T[K]    x_H+     x_H2   (bin means)')
for a, b in zip(bins[:-1], bins[1:]):
    m = (r >= a) & (r < b)
    if m.any():
        print(f'{a:4.0f}-{b:<3.0f} {np.mean(d[m])/1.4:7.2f} {np.mean(T[m]):8.1f} {np.mean(xp[m]):8.4f} {np.mean(xh2[m]):8.4f}')
print('max T', T.max(), 'at r', r[np.argmax(T)], 'n_H', d[np.argmax(T)]/1.4)
