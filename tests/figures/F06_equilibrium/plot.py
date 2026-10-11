"""F06: thermal and chemical equilibrium in an unshielded field chi = 1, xi_cr = 2e-16 s^-1,
z = 1: P/k, T, x_e and x_H2 against n_H. GOW17 + O3,S3,N3 (one-zone harness,
equilibrium.cpp), NCR in Tigris (photchem_equil, run_ncr.sh) and the NCR python reference
of the regression test (Kim J.-G. et al. 2023). Writes
../../../figures/F06_equilibrium.{pdf,png}.

The NCR curves come from reduced.txt when it exists; `python plot.py --from-runs` reads
the NCR run ($PHOTCHEM_RUNS) and the reference ($TIGRIS_DIR) and rewrites it."""
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402
import reduced  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F06_equilibrium')
REDUCED = os.path.join(HERE, 'reduced.txt')
KB = 1.380649e-16

# GOW17 harness: n_H T P/k t_end x_e, then species (H2 is column 5 + 6)
g = np.loadtxt(os.path.join(HERE, 'equilibrium_T1e4.txt'))
with open(os.path.join(HERE, 'equilibrium_T1e4.txt')) as f:
    names = f.readline().split(';')[1].split()
gow = {'nH': g[:, 0], 'T': g[:, 1], 'pok': g[:, 2], 'xe': g[:, 4],
       'xH2': g[:, names.index('H2')], 'xCO': g[:, names.index('CO')]}

if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    ncr_run = paths.runs('T3_onezone', 'F06_ncr')
    # NCR in Tigris: rho in m_H cm^-3 (1.4 n_H)
    d = athena_read.tab(os.path.join(ncr_run, 'equil_ncr.block0.out1.00001.tab'))
    u = athena_read.athinput(os.path.join(ncr_run, 'athinput.runtime'))['units']
    p_cgs = u['mass_cgs']/u['time_cgs']**2/u['length_cgs']
    nh = d['rho']*u['mass_cgs']/u['length_cgs']**3/u['mean_mass_per_hydrogen']
    r = np.genfromtxt(paths.tigris('tst/regression/data/ref_photchem_solutions',
                                   'noshld_zg0.0_zd0.0_chi0.0_xi0.00.txt'), names=True)
    data = {'ncr_nH': nh, 'ncr_pok': d['press']*p_cgs/KB, 'ncr_xe': d['rEL'],
            'ncr_xH2': d['rH2'], 'ref_nH': r['nH'], 'ref_pok': r['pok'], 'ref_xe': r['xe'],
            'ref_xH2': r['xH2']}
    reduced.save(REDUCED, 'F06: NCR in Tigris (T3_onezone/F06_ncr) and the NCR python '
                 'reference\n(tst/regression/data/ref_photchem_solutions/'
                 'noshld_zg0.0_zd0.0_chi0.0_xi0.00.txt)\nn_H [cm^-3], P/k [K cm^-3]',
                 data)
else:
    data = reduced.load(REDUCED)
ncr = {k: data['ncr_' + k] for k in ('nH', 'pok', 'xe', 'xH2')}
ncr['T'] = ncr['pok']/((1.1 + ncr['xe'] - ncr['xH2'])*ncr['nH'])
ref = {k: data['ref_' + k] for k in ('nH', 'pok', 'xe', 'xH2')}
ref['T'] = ref['pok']/((1.1 + ref['xe'] - ref['xH2'])*ref['nH'])

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), layout='constrained')
for dd, lab, kw in [(ref, 'NCR python reference (Kim+ 2023)', dict(color='0.6', lw=4)),
                    (ncr, 'NCR, Tigris', dict(color='C1', lw=1.8)),
                    (gow, 'GOW17 + O3,S3,N3, Tigris', dict(color='C0', lw=1.8, marker='o',
                                                           ms=3))]:
    axs[0].loglog(dd['nH'], dd['pok'], label=lab, **kw)
    axs[1].loglog(dd['nH'], dd['T'], **kw)
    axs[2].loglog(dd['nH'], dd['xe'], **kw)
    axs[2].loglog(dd['nH'], dd['xH2'], ls='--', **{k: v for k, v in kw.items() if k != 'ls'})
axs[0].set(ylabel=r'$P/k$ [K cm$^{-3}$]', ylim=(1e2, 1e5))
axs[1].set(ylabel='T [K]', ylim=(10, 2e4))
axs[2].set(ylabel=r'$x_e$ (solid), $x_{\rm H_2}$ (dashed)', ylim=(1e-5, 1.5))
axs[0].legend(fontsize=9, loc='upper left')
for ax in axs:
    ax.set_xlabel(r'$n_{\rm H}$ [cm$^{-3}$]')
    ax.set_xlim(1e-2, 1e4)
    ax.grid(alpha=0.25, which='both')
fig.suptitle(r'Equilibrium in an unshielded field, $\chi = 1$, $\xi_{\rm CR} = 2\times10^{-16}$ '
             r's$^{-1}$, $Z = 1$', fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
for key in ['T', 'pok', 'xe']:
    for n0 in [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]:
        vg = np.interp(np.log10(n0), np.log10(gow['nH']), gow[key])
        vn = np.interp(np.log10(n0), np.log10(ncr['nH']), ncr[key])
        vr = np.interp(np.log10(n0), np.log10(ref['nH']), ref[key])
        print('%-4s n_H %-7g GOW17 %.4g  NCR %.4g  ref %.4g  GOW17/NCR %.3f'
              % (key, n0, vg, vn, vr, vg/vn))
print(OUT)
