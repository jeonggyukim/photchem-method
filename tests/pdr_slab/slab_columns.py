"""H I -> H2 and C+ -> C transition columns of the GOW17 pdr_slab run, measured as
tst/regression/scripts/tests/photchem/pdr_slab.py measures them, against the
Gong et al. (2017) code values the NCR test encodes."""
import glob, sys
import numpy as np
run, idx = sys.argv[1], sys.argv[2]
def load(out):
    rows, head = [], None
    for f in sorted(glob.glob(f'{run}.block*.{out}.{idx}.tab')):
        for l in open(f):
            if l.startswith('# i'): head = l[1:].split()
            elif not l.startswith('#'): rows.append([float(v) for v in l.split()])
    a = np.array(rows); a = a[np.argsort(a[:, 1])]
    return {h: a[:, n] for n, h in enumerate(head)}
prim, uov = load('out3'), load('out4')
pc, nh0 = 3.0856776e18, 100.0
dx = (prim['x1v'][1] - prim['x1v'][0])*pc
col = (np.arange(len(prim['x1v'])) + 0.5)*nh0*dx
xh2 = prim['r6']                                   # GOW17 H2 slot
def crossing(col, a, b):
    d = a - b; i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0][0]
    f = d[i]/(d[i] - d[i + 1])
    return np.exp(np.log(col[i]) + f*(np.log(col[i + 1]) - np.log(col[i])))
col_h2 = col[np.where(2.0*xh2 >= 0.5)[0][0]]
col_c = crossing(col, uov['xCII'], uov['xCI'])
for name, val, ref in (('2 x(H2) = 0.5', col_h2, 2.213e20), ('x(C+) = x(C)', col_c, 1.844e21)):
    print(f'{name}: {val:.4e} cm^-2, GOW17 code {ref:.4e}, relative {val/ref - 1:+.3f}')
print('max |heat/cool - 1|:', f"{np.max(np.abs(uov['heat_rate']/uov['cool_rate'] - 1)):.2e}")
for c in (1e21, 3e21, 1e22):
    i = np.argmin(abs(col - c)); print(f'N {col[i]:.2e}: T {uov["temp"][i]:.2f} K')
