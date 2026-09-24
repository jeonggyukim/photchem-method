"""F13: radiative SNR in a uniform n_H = 1 medium with five thermal treatments: classic
cooling, NCR, GOW17 core, GOW17 + O3,S3,N3, and the same without the returned grain C, O, Si.
Hot-gas mass, radial momentum and total cooling rate (-dE_tot/dt; periodic box, one SN at
t = 0) against time; shell formation (maximum of M_hot) marked. Writes
../../../figures/F13_snr.{pdf,png}.

Inputs: run_series.sh output in WORKDIR (run_<name>/snr.hst, athinput.runtime, time.out)."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'F13_snr')
WORKDIR = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion/M6_rad_snr/'
                             'F13_series')
RUNS = [('classic', 'classic cooling', 'C7', ':'),
        ('ncr', 'NCR', 'k', '-'),
        ('core', 'GOW17 core', 'C2', '--'),
        ('ions', 'GOW17 + O3,S3,N3', 'C0', '-'),
        ('noreturn', 'GOW17 + O3,S3,N3, no returned C, O, Si', 'C3', '-.')]

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ref = {}
for name, lab, col, ls in RUNS:
    run = os.path.join(WORKDIR, 'run_' + name)
    h = athena_read.hst(os.path.join(run, 'snr.hst'))
    u = athena_read.units(athena_read.athinput(os.path.join(run, 'athinput.runtime')))
    t = h['time']/u['Myr']
    isf = h['M_hot'].argmax()
    msf, psf = h['M_hot'][isf]/u['solar_mass'], h['pr'][isf]/(u['solar_mass']*u['km_s'])
    if name == 'ncr':
        ref = {'t': t[isf], 'M': msf, 'p': psf}
    wall = [float(ln.split()[1]) for ln in open(os.path.join(run, 'time.out'))
            if ln.startswith('real')][0]
    # cooling rate from the total energy, on a uniform grid of 200 times to damp the
    # per-cycle noise of the history
    tg = np.linspace(t[0], t[-1], 200)
    eg = np.interp(tg, t, h['tot-E']/u['erg'])
    lcool = -np.gradient(eg, tg*1e6*3.15576e7)
    label = '%s (%.1f s wall)' % (lab, wall)
    axs[0].plot(t, h['M_hot']/u['solar_mass'], color=col, ls=ls, lw=1.6, label=label)
    axs[1].plot(t, h['pr']/(u['solar_mass']*u['km_s']), color=col, ls=ls, lw=1.6)
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
