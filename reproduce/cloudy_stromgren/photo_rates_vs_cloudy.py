"""Photoionization rate per ion in the Tigris static sphere (5 bands) against
Cloudy's total ionization rate per ion (save ionization rates), at fixed radii."""
import sys, re, glob, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read()
run, cdir = sys.argv[1], sys.argv[2]
txt = open(run + '/athinput.runtime').read()
blk = txt[txt.index('<photchem_gow17>'):]
blk = blk[:blk.find('\n<', 5)] if blk.find('\n<', 5) > 0 else blk
def key(k):
    return float(re.search(r'^%s\s*=\s*([-+.\deE]+)' % re.escape(k), blk, re.M).group(1))
L, Tu, Eu = 3.0856776e18, 3.0856776e13, 4.91769147364387e31*1e10
d = athena_read.athdf(sorted(glob.glob(run + '/HII.out2.*.athdf'))[-1])
X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
R = np.sqrt(X1**2 + X2**2 + X3**2)
nion = sum(1 for g in range(8) if 'photon_frac[%d]' % g in blk) or 1
def xi(sp):
    # Er_rayt is the energy density E; the rate needs c E, c = 2.998e5 km/s in code units
    return sum(2.99792458e5*d['Er_rayt%d' % g].astype(float)*(key('sigma_pi_%s[%d]' % (sp, g))/L**2)
               /(key('hnu[%d]' % g)*1.602176634e-12/Eu)/Tu for g in range(nion))
def load(f):
    rows = [l.split() for l in open(f) if not l.startswith('#')]
    n = min(len(r) for r in rows); return np.array([[float(x) for x in r[:n]] for r in rows])
cS, cO = load(cdir + '/hii.ionr_S'), load(cdir + '/hii.ionr_O')
cr = (10**16.5 + cS[:, 0])/3.0857e18
print('ionization rate per ion [s^-1]: Tigris photo (%d ionizing bands) / Cloudy total' % nion)
for r0 in (0.3, 1.0, 1.5, 2.0, 2.5):
    m = (R > r0 - 0.1) & (R < r0 + 0.1); i = np.argmin(abs(cr - r0)); out = 'r=%.1f pc' % r0
    for lab, sp, A, st in (('S+', 'SII', cS, 1), ('S++', 'SIII', cS, 2), ('O+', 'OII', cO, 1)):
        t = np.mean(xi(sp)[m]); c = A[i, 3 + 4*st + 1]
        out += '  %s %.2e / %.2e = %.2f' % (lab, t, c, t/c)
    print(out)
