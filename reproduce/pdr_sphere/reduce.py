"""Mid-plane slices of the post-processed sphere, for scripts/pdr_sphere_slices.py.

    python reduce.py <run_dir> <out.npz> [athena_read directory]

Write out.npz to ../data/pdr_sphere/pdr_sphere_<mode>.npz, outside the repository.

Reads the last HDF5 dump of a run of run.sh (prim and uov, cloud.out2.*.athdf) and
keeps the plane x3 = 0 through the centre of the sphere: n_H, x_e, x_H, x_H2, the
carbon fractions, T, T_d (NaN where the mode does not solve it), the fields and xi_cr, plus the iteration count from run.log.
The GOW17 output names its first scalar rHI although it holds He+; its species
follow gow17_core.hpp: He+ OHx CHx CO C+ HCO+ H2 H+ H3+ H2+ O+ Si+.
"""
import glob
import re
import sys

import numpy as np

run, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, sys.argv[3] if len(sys.argv) > 3 else
                '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read  # noqa: E402

dump = sorted(glob.glob(run + '/cloud.out2.*.athdf'))[-1]
d = athena_read.athdf(dump)
k = int(np.argmin(np.abs(d['x3v'])))
gow17 = 'r6' in d
p = {n: d[n][k] for n in d if isinstance(d[n], np.ndarray) and d[n].ndim == 3}
if gow17:
    x_e = sum(p['r%d' % i] for i in (4, 5, 7, 8, 9, 10, 11)) + p['rHI']
    x_h2 = p['r6']
    x_h = 1.0 - 2.0*p['r6'] - p['r7'] - 2.0*p['r8'] - 1.5*p['r9']
else:
    x_e, x_h2, x_h = p['rEL'], p['rH2'], p['rHI']
text = open(run + '/run.log').read()
m = re.search(r'converged after (\d+) iterations', text)
np.savez_compressed(
    out, mode='gow17' if gow17 else 'ncr', niter=int(m.group(1)) if m else -1,
    x1v=d['x1v'], x2v=d['x2v'], x3=d['x3v'][k],
    nH=p['rho']/1.4, x_e=x_e, x_h=x_h, x_h2=x_h2,
    xCII=p['xCII'], xCI=p['xCI'], xCO=p['xCO'], temp=p['temp'],
    chi_PE=p['chi_PE'], chi_LW=p['chi_LW'], chi_H2=p['chi_H2'], chi_CI=p['chi_CI'],
    chi_CO=p['chi_CO'] if 'chi_CO' in p else np.full_like(p['temp'], np.nan),
    xi_CR=p['xi_CR'],
    temp_dust=p['temp_dust'] if 'temp_dust' in p and np.any(p['temp_dust'] > 0)
    else np.full_like(p['temp'], np.nan))
print('wrote', out, 'from', dump)
