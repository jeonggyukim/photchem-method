"""The turbulent cloud M5n20 (athinput.hii_cloud_m5n20_64, GOW17, and its NCR twin
athinput.hii_cloud_m5n20_64_ncr), 64^3: gas mass in ionized, atomic and molecular
hydrogen and the timestep against time, for the two chemical networks. Writes
../../../../figures/ncr_gow17_m5n20.png.

The curves come from reduced_m5n20.txt; `python plot.py --from-runs` reads
$PHOTCHEM_RUNS/M5_hii_dtype_ions/cloud_m5n20_64{,_ncr}/run and rewrites it. The GOW17
run was restarted at t = 3 (code) with the timescale history; its history rows before
the restart come from part1 of the run, which has four columns fewer.
Usage: python plot.py [--from-runs]
"""
import glob
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from photchem_method import paths, reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', '..', 'figures', 'ncr_gow17_m5n20')
REDUCED = os.path.join(HERE, 'reduced_m5n20.txt')
MYR = 0.978                  # code time [Myr]
PC_CM, MH, MSUN = 3.0857e18, 1.6735575e-24, 1.989e33
MU_H = 1.4                   # gas mass per H [m_H]
GOW17_NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+']
T0_SRC = 2.58                # source on [code]


def phase_masses(fn, ncr):
    athena_read = paths.athena_read()
    d = athena_read.athdf(fn)
    dv = (float(np.diff(d['x1f'])[0])*PC_CM)**3
    nh = d['rho'].astype(np.float64)/MU_H
    if ncr:
        x_h2 = d['rH2']
        x_hp = np.clip(1.0 - d['rHI'] - 2.0*x_h2, 0.0, 1.0)
    else:
        x_h2 = d['r%d' % GOW17_NAMES.index('H2')]
        x_hp = d['r%d' % GOW17_NAMES.index('H+')]
    x_hi = np.clip(1.0 - x_hp - 2.0*x_h2, 0.0, 1.0)
    m = MU_H*MH*dv/MSUN
    return (float(d['Time']), m*np.sum(nh*x_hp), m*np.sum(nh*x_hi),
            m*np.sum(nh*2.0*x_h2))


def history(fn):
    hdr = [line for line in open(fn) if line.startswith('#')][-1]
    names = [m.split('=')[1] for m in re.findall(r'\[\d+\]=\S+', hdr)]
    rows = [line.split() for line in open(fn) if not line.startswith('#')]
    t = np.array([float(r[0]) for r in rows])
    dt = np.array([float(r[names.index('dt')]) for r in rows])
    return t, dt


def from_runs():
    base = os.path.join(paths.env('PHOTCHEM_RUNS'), 'M5_hii_dtype_ions')
    arrays, sources = {}, []
    for key, sub, ncr in [('gow17', 'cloud_m5n20_64', False),
                          ('ncr', 'cloud_m5n20_64_ncr', True)]:
        run = os.path.join(base, sub, 'run')
        sources.append(run)
        rows = [phase_masses(fn, ncr)
                for fn in sorted(glob.glob(os.path.join(run, 'HII.out2.*.athdf')))]
        a = np.array(rows)
        for n, name in enumerate(['time', 'M_Hp', 'M_HI', 'M_H2']):
            arrays['%s_%s' % (key, name)] = a[:, n]
        t, dt = history(os.path.join(run, 'HII.hst'))
        keep = np.r_[True, np.diff(t) > 0]   # the restart repeats t = 3 to 3.37
        arrays[key + '_hst_time'], arrays[key + '_hst_dt'] = t[keep], dt[keep]
    reduced.save(REDUCED, 'gas masses [Msun] of each hydrogen phase from the dumps; '
                 'dt [code] from the history', arrays, sources=sources)


def main():
    if '--from-runs' in sys.argv:
        from_runs()
    r = reduced.load(REDUCED)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    for key, lab, ls in [('gow17', 'GOW17', '-'), ('ncr', 'NCR', '--')]:
        t = (r[key + '_time'] - T0_SRC)*MYR
        for q, qlab, col in [('M_Hp', r'H$^+$', 'C3'), ('M_HI', r'H$^0$', 'C0'),
                             ('M_H2', r'H$_2$', 'C2')]:
            ax[0].plot(t, r['%s_%s' % (key, q)], ls, color=col, lw=1.2,
                       label='%s %s' % (lab, qlab))
        th = (r[key + '_hst_time'] - T0_SRC)*MYR
        ax[1].semilogy(th, r[key + '_hst_dt'], ls, color='k', lw=0.8, label=lab)
    ax[0].set_ylabel(r'gas mass [M$_\odot$]')
    ax[0].set_yscale('log')
    ax[1].set_ylabel('dt [code]')
    for a in ax:
        a.set_xlabel(r'$t - t_{\rm src}$ [Myr]')
        a.axvline(0.0, color='0.6', lw=0.6)
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].legend(fontsize=8)
    fig.suptitle(r'M5n20, $64^3$ (1.82 pc cells), GOW17 and NCR', fontsize=10)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for ext in figure_formats():
        fig.savefig(OUT + '.' + ext, dpi=120)
    print(OUT)


if __name__ == '__main__':
    main()
