"""Heat per H photoionization [eV]: Cloudy (H I heating / (n_H0 Gamma_H)) against
Tigris (band-rate-weighted dhnu from the local band fields), per radius.
usage: python heat_per_ionization.py CLOUDY_DIR TIGRIS_RUNDIR"""
import sys, re, glob, numpy as np
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python'); import athena_read
C, run = sys.argv[1].rstrip('/') + '/', sys.argv[2]
eV, pc, Rin = 1.602176634e-12, 3.0857e18, 10**16.5
H = [l.split() for l in open(C + 'hii.ionr_H') if l.startswith('hion')]
gam = np.array([float(r[2]) for r in H]); hiihi = np.array([float(r[6]) for r in H])
heat = [l.rstrip('\n').split('\t') for l in open(C + 'hii.heat') if not l.startswith('#')]
depth = np.array([float(r[0]) for r in heat]); Htot = np.array([float(r[2]) for r in heat])
def frac(r, lab):
    for k in range(5, len(r) - 1, 2):
        if r[k].strip() == lab: return float(r[k + 1])
    return 0.0
fH = np.array([frac(r, 'H  1') for r in heat])
n = min(len(gam), len(depth)); r_pc = (Rin + depth[:n])/pc
dE_c = fH[:n]*Htot[:n]/(100/(1 + hiihi[:n])*gam[:n])/eV
txt = open(run + '/athinput.runtime').read(); blk = txt[txt.index('<photchem_gow17>'):]
key = lambda k: float(re.search(r'^%s\s*=\s*([-+.\deE]+)' % re.escape(k), blk, re.M).group(1))
nion = sum(1 for g in range(8) if re.search(r'^photon_frac\[%d\]' % g, blk, re.M)) - 2
d = athena_read.athdf(sorted(glob.glob(run + '/HII.out2.*.athdf'))[-1])
X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij'); R = np.sqrt(X1**2 + X2**2 + X3**2)
w = [d['Er_rayt%d' % g].astype(float)*key('sigma_pi_HI[%d]' % g)/key('hnu[%d]' % g) for g in range(nion)]
num = sum(wg*key('dhnu_pi_HI[%d]' % g) for g, wg in enumerate(w)); den = sum(w)
print(' r[pc]  heat per H ionization [eV]: Cloudy  Tigris')
for r0 in (0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8):
    i = np.argmin(abs(r_pc - r0)); m = (R > r0 - 0.1) & (R < r0 + 0.1)
    print('%5.1f   %6.2f  %6.2f' % (r0, dE_c[i], np.sum(num[m])/np.sum(den[m])))
