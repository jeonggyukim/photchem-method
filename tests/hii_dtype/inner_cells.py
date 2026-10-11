import glob, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read()
athena_read.check_nan_flag = False
run, idx = sys.argv[1], sys.argv[2]
rows = []
for f in sorted(glob.glob(f'{run}.block*.out2.{idx}.vtk')):
    xf, yf, zf, d = athena_read.vtk(f)
    xc, yc, zc = [0.5*(a[1:] + a[:-1]) for a in (xf, yf, zf)]
    Z, Y, X = np.meshgrid(zc, yc, xc, indexing='ij')
    r = np.sqrt(X**2 + Y**2 + Z**2)
    sc = [d['rHI']] + [d[f'r{n}'] for n in range(1, 12)]
    xe = sc[0] + sc[4] + sc[5] + sc[7] + sc[8] + sc[9] + sc[10] + sc[11]
    nH = d['rho']/1.4
    T = d['press']*1.6738234e-14/1.38065e-16/(nH*(1.1 + xe - sc[6]))
    xHI = 1 - sc[7] - 2*sc[6] - 3*sc[8] - 2*sc[9] - sc[1] - sc[2] - sc[5]
    for m in zip(*[a.ravel() for a in (r, nH, T, sc[7], xHI, sc[6], sc[9], sc[8], sc[0], xe, d['Er_rayt0'])]):
        rows.append(m)
rows = np.array(rows); rows = rows[np.argsort(rows[:, 0])]
print('   r     n_H      T      x_H+     x_HI     x_H2     x_H2+    x_H3+    x_He+    x_e      Er0')
for row in rows[:12]:
    print(' '.join(f'{v:8.3g}' for v in row))
m = rows[:, 3] > 0.01
print('cells with x_H+ > 0.01:', m.sum(), 'max r', rows[m, 0].max())
