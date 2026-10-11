"""F13: radiative SNR in a uniform n_H = 1 medium with five thermal treatments: classic
cooling, NCR, GOW17 core, GOW17 + O3,S3,N3, and the same without the returned grain C, O, Si.
Hot-gas mass, radial momentum and total cooling rate (-dE_tot/dt; periodic box, one SN at
t = 0) against time; shell formation (maximum of M_hot) marked. Writes
../../../figures/F13_snr.{pdf,png}.

Inputs: run_series.sh output in WORKDIR (run_<name>/snr.hst, athinput.runtime, time.out).
The histories and wall times come from reduced.txt when it exists;
`python plot.py --from-runs` reads WORKDIR ($PHOTCHEM_RUNS/M6_rad_snr/F13_series) and
rewrites it."""
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
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F13_snr')
REDUCED = os.path.join(HERE, 'reduced.txt')
RUNS = [('classic', 'classic cooling', 'C7', ':'),
        ('ncr', 'NCR', 'k', '-'),
        ('core', 'GOW17 core', 'C2', '--'),
        ('ions', 'GOW17 + O3,S3,N3', 'C0', '-'),
        ('noreturn', 'GOW17 + O3,S3,N3, no returned C, O, Si', 'C3', '-.')]

if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    data = {}
    for name, *_ in RUNS:
        run = paths.runs('M6_rad_snr', 'F13_series', 'run_' + name)
        h = athena_read.hst(os.path.join(run, 'snr.hst'))
        u = athena_read.units(athena_read.athinput(os.path.join(run, 'athinput.runtime')))
        data['t_' + name] = h['time']/u['Myr']
        data['M_hot_' + name] = h['M_hot']/u['solar_mass']
        data['pr_' + name] = h['pr']/(u['solar_mass']*u['km_s'])
        data['E_tot_' + name] = h['tot-E']/u['erg']
        data['wall_' + name] = [float(ln.split()[1]) for ln in
                                open(os.path.join(run, 'time.out'))
                                if ln.startswith('real')][0]
    reduced.save(REDUCED, 'F13: histories of the runs M6_rad_snr/F13_series/run_<name> and '
                 'their wall time\nt [Myr], M_hot [M_sun], pr [M_sun km/s], E_tot [erg], '
                 'wall [s]', data)
else:
    data = reduced.load(REDUCED)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ref = {}
for name, lab, col, ls in RUNS:
    t, mhot, pr = data['t_' + name], data['M_hot_' + name], data['pr_' + name]
    isf = mhot.argmax()
    msf, psf = mhot[isf], pr[isf]
    if name == 'ncr':
        ref = {'t': t[isf], 'M': msf, 'p': psf}
    wall = float(data['wall_' + name])
    # cooling rate from the total energy, on a uniform grid of 200 times to damp the
    # per-cycle noise of the history
    tg = np.linspace(t[0], t[-1], 200)
    eg = np.interp(tg, t, data['E_tot_' + name])
    lcool = -np.gradient(eg, tg*1e6*3.15576e7)
    label = '%s (%.1f s wall)' % (lab, wall)
    axs[0].plot(t, mhot, color=col, ls=ls, lw=1.6, label=label)
    axs[1].plot(t, pr, color=col, ls=ls, lw=1.6)
    axs[2].semilogy(tg[1:], lcool[1:], color=col, ls=ls, lw=1.6)
    for ax in axs:
        ax.axvline(t[isf], color=col, ls=ls, lw=0.8, alpha=0.6)
    print('%-9s t_sf %.4f Myr  M_sf %.1f M_sun  p_sf %.4g M_sun km/s  wall %.1f s'
          % (name, t[isf], msf, psf, wall), end='')
    if ref:
        print('  vs NCR: %+.3f %+.3f %+.3f' % (t[isf]/ref['t'] - 1, msf/ref['M'] - 1,
                                              psf/ref['p'] - 1))
    else:
        print()
axs[0].set(ylabel=r'$M_{\rm hot}$ [$M_\odot$]')
axs[1].set(ylabel=r'$p_r$ [$M_\odot$ km s$^{-1}$]')
axs[2].set(ylabel=r'$-dE_{\rm tot}/dt$ [erg s$^{-1}$]', ylim=(1e36, None))
axs[1].legend(*axs[0].get_legend_handles_labels(), fontsize=8, loc='lower right')
for ax in axs:
    ax.set_xlabel('t [Myr]')
    ax.grid(alpha=0.25)
fig.suptitle(r'Radiative SNR, $10^{51}$ erg in $n_{\rm H} = 1$ cm$^{-3}$, $64^3$: '
             'shell formation (thin vertical lines) at the maximum of $M_{\\rm hot}$',
             fontsize=11)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT + '.' + ext, dpi=300)
print(OUT)
