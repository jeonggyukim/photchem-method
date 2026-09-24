"""Species and T of two Tigris GOW17 slab runs against Athena++ chem_pdr_static.vtk,
as Tigris/Athena++ ratios at a few columns."""
import glob, sys
import numpy as np
sys.path.insert(0, '/Users/jgkim/Projects/athena-pp-pdr1d/vis/python')
import athena_read
names = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
ch = {'He+', 'C+', 'HCO+', 'H+', 'H3+', 'H2+', 'O+', 'Si+'}
pc, nh0, kB, eu = 3.0856776e18, 100.0, 1.380649e-16, 1.6738234e-24*1e10
def prof(get, x, press, fac):
    p = {n: get(n) for n in names}
    p['e'] = sum(p[n] for n in ch)
    p['T'] = press*fac*eu/kB/(nh0*(1.1 + p['e'] - p['H2'])); p['N'] = x*pc*nh0
    return p
x, _, _, d = athena_read.vtk('/Users/jgkim/Projects/athena-pp-pdr1d/tst/regression/data/chem_pdr_static.vtk')
ap = prof(lambda n: d['r' + n].ravel(), 0.5*(x[1:] + x[:-1]), d['press'].ravel(), 1.4)
def tig(run):
    rows, head = [], None
    for f in sorted(glob.glob(f'{run}.block*.out3.00001.tab')):
        for l in open(f):
            if l.startswith('# i'): head = l[1:].split()
            elif not l.startswith('#'): rows.append([float(v) for v in l.split()])
    a = np.array(rows); a = a[np.argsort(a[:, 1])]; c = {h: a[:, k] for k, h in enumerate(head)}
    return prof(lambda n: c['rHI' if names.index(n) == 0 else f'r{names.index(n)}'], c['x1v'], c['press'], 1.0)
runs = [(lab, tig(r)) for lab, r in zip(sys.argv[1::2], sys.argv[2::2])]
at = lambda p, q, N: np.exp(np.interp(np.log(N), np.log(p['N']), np.log(np.maximum(p[q], 1e-30))))
print(f"{'q':>5s} {'N':>7s} " + ' '.join(f'{lab:>9s}' for lab, _ in runs) + '   (Tigris/Athena++)')
for q in ('T', 'CO', 'CHx', 'OHx', 'HCO+', 'C+', 'e', 'H2'):
    for N in (1e21, 3e21, 1e22):
        print(f'{q:>5s} {N:7.0e} ' + ' '.join(f'{at(p, q, N)/at(ap, q, N):9.3f}' for _, p in runs))
