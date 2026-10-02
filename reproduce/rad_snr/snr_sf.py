"""Shell formation of the radiative SNR run, measured as
tst/regression/scripts/tests/feedback/rad_snr.py measures it, against the NCR
reference values stored there (tolerance 5%)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read()
hst, runtime = sys.argv[1], sys.argv[2]
h = athena_read.hst(hst); par = athena_read.athinput(runtime); u = athena_read.units(par)
isf = h['M_hot'].argmax()
vals = {'t_sf': h['time'][isf], 'M_sf': h['M_hot'][isf]/u['solar_mass'],
        'p_sf': h['pr'][isf]/(u['solar_mass']*u['km_s'])}
if 'Mej_hot' in h: vals['f_ej'] = (h['Mej_hot']/h['M_hot'])[isf]
ref = {'t_sf': 0.038473, 'M_sf': 1483.49, 'p_sf': 192798., 'f_ej': 0.00673852}
for k, v in vals.items():
    print(f'{k}: {v:.6g}  NCR ref {ref[k]:.6g}  relative {v/ref[k] - 1:+.3f}')
m = h['mass'][isf]
base = int(par['scalar_indices']['IHI'])  # first photchem scalar
x = lambda n: h[f'{base + n}-scalar'][isf]/m
# GOW17 order: He+ OHx CHx CO C+ HCO+ H2 H+ H3+ H2+ O+ Si+; x_e without hot-gas electrons
xhp, xh2 = x(7), x(6)
xe = sum(x(n) for n in (0, 4, 5, 7, 8, 9, 10, 11))
print(f'mass-weighted at t_sf: x_H+ {xhp:.4f}, x_HI {1 - xhp - 2*xh2:.4f}, x_e(species) {xe:.4f}; '
      'NCR refs x_HI 0.8159, x_e 0.2128')
