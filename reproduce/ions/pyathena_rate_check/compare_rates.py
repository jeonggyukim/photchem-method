"""Compare pyathena (PhotChem) ion rate coefficients with Tigris (tigris_rates.csv)."""
import os, sys, types
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.expanduser('~/Dropbox/Projects/pyathena'))
from pyathena.microphysics.rec_rate import RecRate
from pyathena.microphysics.ci_rate import CollIonRate
from pyathena.microphysics.ct_rate import ChargeTransferRate
from pyathena.microphysics.photchem import PhotChem

here = os.path.dirname(os.path.abspath(__file__))
tg = pd.read_csv(os.path.join(here, 'tigris_rates.csv'))
rc, ci, ct = RecRate(), CollIonRate(), ChargeTransferRate()
pc = types.SimpleNamespace(ct=ct)


def ct_photchem(kind, Z, N, T):
    # PhotChem.evolve_one_species skips CT for He; metals go through _ct_rate_safe.
    if Z == 2:
        return 0.0
    return float(PhotChem._ct_rate_safe(pc, kind, Z, N, T))


rows = []
for _, r in tg.iterrows():
    Z, q, T = int(r.Z), int(r.q), float(r['T'])
    nrec, nion = Z - q - 1, Z - q
    py = dict(rec=float(rc.get_rec_rate(Z, nrec, T)),
              ci=float(ci.get_ci_rate(Z, nion, T)),
              ctrec=ct_photchem('rec', Z, nrec, T),
              ction=ct_photchem('ion', Z, nion, T))
    raw = dict(ctrec=float(ct.get_ct_rec_rate(Z, nrec, T)),
               ction=float(ct.get_ct_ion_rate(Z, nion, T)))
    for p in ('rec', 'ci', 'ctrec', 'ction'):
        rows.append(dict(ion=r.ion, T=T, proc=p, tigris=r[p], pyathena=py[p],
                         py_raw=raw.get(p, py[p])))
df = pd.DataFrame(rows)
thr = 1e-15
big = np.maximum(df.tigris, df.pyathena)
df['rel'] = np.where((df.tigris > 0) & (df.pyathena > 0), df.pyathena/df.tigris - 1, np.nan)
df['zero_mismatch'] = ((df.tigris == 0) != (df.pyathena == 0)) & (big > thr)

out = []
out.append('# Tigris vs pyathena ion rate coefficients [cm^3 s^-1]')
out.append('# rel = pyathena/tigris - 1; stage q = lower stage (rec, ctrec: X^(q+1)->X^q; ci, ction: X^q->X^(q+1))')
out.append('# pyathena CT = values PhotChem uses (_ct_rate_safe; 0 for He); py_raw = ChargeTransferRate.get_ct_*_rate')
out.append('')
out.append('## Summary: max |rel| over T where both rates > 1e-15 cm^3 s^-1')
out.append(f"{'ion':4s} {'proc':6s} {'n_T':>3s} {'max|rel|':>9s} {'at T[K]':>8s}  zero-mismatch T [K] (one code 0, other > 1e-15)")
for (ion, p), g in df.groupby(['ion', 'proc'], sort=False):
    ok = g[(g.tigris > thr) & (g.pyathena > thr)]
    zm = g[g.zero_mismatch]
    zs = ','.join(f'{t:.0e}' for t in zm['T']) if len(zm) else '-'
    if len(ok):
        i = ok.rel.abs().idxmax()
        out.append(f"{ion:4s} {p:6s} {len(ok):3d} {abs(ok.rel[i]):9.2e} {ok['T'][i]:8.0e}  {zs}")
    else:
        out.append(f"{ion:4s} {p:6s} {0:3d} {'-':>9s} {'-':>8s}  {zs}")
out.append('')
out.append('## Full table')
out.append(f"{'ion':4s} {'proc':6s} {'T[K]':>8s} {'tigris':>11s} {'pyathena':>11s} {'py_raw':>11s} {'rel':>10s}")
for _, r in df.iterrows():
    rel = f'{r.rel:+10.3e}' if np.isfinite(r.rel) else f"{'n/a':>10s}"
    out.append(f"{r.ion:4s} {r.proc:6s} {r['T']:8.0e} {r.tigris:11.4e} {r.pyathena:11.4e} {r.py_raw:11.4e} {rel}")
txt = '\n'.join(out) + '\n'
open(os.path.join(here, 'rate_comparison.txt'), 'w').write(txt)
print(txt.split('## Full table')[0])
