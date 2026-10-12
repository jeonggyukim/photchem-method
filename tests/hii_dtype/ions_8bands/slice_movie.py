"""z = 0 slice movie of the expanding H II region of run_tigris.sh (GOW17 + O3,S3,N3
and He++, 8 bands, OTS), laid out as ../../hii_static/fig_slices: first column n_H, T,
P/k and |v|; with --radiation log U (ionizing photons per H), chi_PE and the shielded
chi_H2 and chi_CO; with --molecular H2, CO, C+ and C0; then H, He, O, S, N, one row per
stage (neutral, singly, doubly and higher ionized). Every panel keeps one colour range
over all frames.

usage: python slice_movie.py RUNDIR OUTDIR [--nproc N] [--frame I] [--molecular]
       [--radiation] [--temp-log] [--density-scale S] [--title TEXT]
Writes OUTDIR/frame_NNNNN.png, one per HII.out2 dump, then OUTDIR/slices.mp4 (ffmpeg).
The frames are independent and are drawn by N processes (default 8). --frame I draws
only dump I, without the movie, to check the layout first.
"""
import argparse
import functools
import glob
import re
import multiprocessing
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import FuncNorm, LinearSegmentedColormap  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from photchem_method import paths  # noqa: E402

athena_read = paths.athena_read()
# ism units: density unit m_H/cm^3, rho = 1.4 n_H; pressure unit m_H (km/s)^2;
# time unit 0.978 Myr
MH, KB, MU_H, X_HE = 1.6735575e-24, 1.380649e-16, 1.4, 0.1
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++',
         'N_high', 'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
XTOT = {'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
IONS = {'He': ['He+', 'He++'], 'O': ['O+', 'O++', 'O_high'],
        'S': ['S+', 'S++', 'S_high'], 'N': ['N+', 'N++', 'N_high']}


def layout(molecular, radiation):
    """Columns top to bottom, as (columns with their own colour bars, columns that share
    the fraction colour bar). --radiation adds the band fields, --molecular H2 and the
    carbon carriers."""
    own = [['nH', 'T', 'Pk', 'v']] + \
        ([['logU', 'chiPE', 'chiH2', 'chiCO']] if radiation else [])
    shared = ([['H2', 'CO', 'C+', 'C0']] if molecular else []) + \
        [['H0', 'H+', None, None], ['He0', 'He+', 'He++', None]] + \
        [[el + '0', el + '+', el + '++', el + '_high'] for el in 'OSN']
    return own, shared


# ism units: energy density unit m_H (km/s)^2 / cm^3 [erg cm^-3]; the point-source bands
# of the 8-band set are 0-5 ionizing (0-2 above the He I edge, 24.59 eV), 6 LW, 7 PE;
# chi in Draine units, as the chemical network normalizes them
E_CODE = 1.6738234e-14
U_DRAINE_LW, U_DRAINE_PE = 1.335e-14, 7.613e-14
# speed of light [km/s]: the diffuse groups' mean intensity Jdiff gives u = 4 pi J/c
C_CODE = 2.99792458e5
EV = 1.602176634e-12


@functools.lru_cache(maxsize=None)
def band_energies(rundir):
    """Mean photon energy of each band [eV], hnu[b] of the run's athinput.runtime."""
    hnu = {}
    with open(os.path.join(rundir, 'athinput.runtime')) as fin:
        for line in fin:
            m = re.match(r'\s*hnu\[(\d+)\]\s*=\s*([0-9.eE+-]+)', line)
            if m:
                hnu[int(m.group(1))] = float(m.group(2))
    return [hnu[b] for b in sorted(hnu)]


LAB = {'nH': r'$n_{\rm H}$ [cm$^{-3}$]', 'T': r'$T$ [K]', 'Pk': r'$P/k$ [K cm$^{-3}$]',
       'v': r'$|v|$ [km s$^{-1}$]', 'H0': r'H$^0$/H', 'H+': r'H$^+$/H',
       'He0': r'He$^0$/He', 'He+': r'He$^+$/He', 'He++': r'He$^{2+}$/He',
       'H2': r'2H$_2$/H', 'CO': 'CO/C', 'C+': r'C$^+$/C', 'C0': r'C$^0$/C',
       'logU': r'$\log U = \log(n_\gamma/n_{\rm H})$', 'chiPE': r'$\chi_{\rm PE}$',
       'chiH2': r'$\chi_{\rm H_2}$', 'chiCO': r'$\chi_{\rm CO}$'}
X_C = 1.6e-4
for el in 'OSN':
    LAB.update({el + '0': r'%s$^0$/%s' % (el, el), el + '+': r'%s$^+$/%s' % (el, el),
                el + '++': r'%s$^{2+}$/%s' % (el, el),
                el + '_high': r'%s$^{\geq3+}$/%s' % (el, el)})
NORM = {'nH': LogNorm(10, 1e3), 'T': Normalize(5e3, 1.1e4), 'Pk': LogNorm(1e3, 3e6),
        'v': Normalize(0, 15), 'ion': LogNorm(1e-5, 1), 'logU': Normalize(-6.0, 0.0),
        'chiPE': LogNorm(1e-2, 1e5), 'chiH2': LogNorm(1e-6, 1e3),
        'chiCO': LogNorm(1e-6, 1e3)}
CMAP = {'nH': 'cividis', 'T': 'inferno', 'Pk': 'magma', 'v': 'plasma',
        'ion': plt.get_cmap('viridis').copy(), 'logU': 'magma', 'chiPE': 'BuPu_r',
        'chiH2': 'BuPu_r', 'chiCO': 'BuPu_r'}
CMAP['ion'].set_bad(CMAP['ion'](0.0))   # a fraction of exactly 0 takes the lowest colour
# --temp-log: T in a colour map with one hue family per phase, on a two-segment scale so
# that both the neutral and the ionized gas show their structure: 10-5000 K on a log
# scale over the lower T_SPLIT of the bar, 5000-11000 K linear over the rest.
#   10-100 K      molecular, cold      navy to blue
#   100-1000 K    warm PDR and shell   blue to teal to pale green
#   1000-5000 K   warm neutral         pale green to pale yellow
#   5000-11000 K  ionized              gold to orange to red to dark red
T_LOG = (10.0, 5.0e3, 1.1e4)
T_SPLIT = 0.6


def _t_forward(t):
    t = np.clip(np.asarray(t, dtype=float), T_LOG[0], T_LOG[2])
    lo = T_SPLIT*np.log10(t/T_LOG[0])/np.log10(T_LOG[1]/T_LOG[0])
    hi = T_SPLIT + (1.0 - T_SPLIT)*(t - T_LOG[1])/(T_LOG[2] - T_LOG[1])
    return np.where(t < T_LOG[1], lo, hi)


def _t_inverse(s):
    s = np.asarray(s, dtype=float)
    lo = T_LOG[0]*(T_LOG[1]/T_LOG[0])**(s/T_SPLIT)
    hi = T_LOG[1] + (s - T_SPLIT)/(1.0 - T_SPLIT)*(T_LOG[2] - T_LOG[1])
    return np.where(s < T_SPLIT, lo, hi)


def _at(t):
    return float(_t_forward(t))


CMAP_T_LOG = LinearSegmentedColormap.from_list('T_phases', [
    (0.0, '#0b1d51'), (_at(100.0), '#2b8cbe'), (_at(3.2e2), '#41b6c4'),
    (_at(1e3), '#a1dab4'), (_at(4.5e3), '#f4f9c4'), (_at(5e3), '#fed976'),
    (_at(7.5e3), '#fd8d3c'), (_at(9.5e3), '#d7301f'), (1.0, '#67000d')])
T_TICKS = [10, 100, 1000, 5000, 7000, 9000]


def fields(fn):
    d = athena_read.athdf(fn)
    k = int(np.where(d['x3v'] > 0, d['x3v'], np.inf).argmin())  # first plane above 0
    x = {nm: d[key][k].astype(np.float64)
         for nm, key in zip(NAMES, ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))])}
    f = {'nH': d['rho'][k]/MU_H}
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    f['Pk'] = d['press'][k]*MH*1e10/KB
    f['T'] = f['Pk']/(f['nH']*(1.0 - x['H2'] + X_HE + xe))
    f['v'] = np.sqrt(d['vel1'][k]**2 + d['vel2'][k]**2 + d['vel3'][k]**2)
    f['H0'] = np.maximum(1.0 - x['H+'] - 2*x['H2'], 0.0)
    f['H+'] = x['H+']
    f['H2'] = 2*x['H2']
    f['CO'], f['C+'] = x['CO']/X_C, x['C+']/X_C
    f['C0'] = np.maximum(1.0 - (x['CO'] + x['C+'] + x['CHx'] + x['HCO+']
                                + x['C_high'])/X_C, 0.0)
    for el, ions in IONS.items():
        for nm in ions:
            f[nm] = x[nm]/XTOT[el]
        f[el + '0'] = np.maximum(1.0 - sum(f[nm] for nm in ions), 0.0)
    if 'chi_H2' in d:
        u = [d['Er_rayt%d' % b][k].astype(np.float64)*E_CODE for b in range(8)]
        # ionizing photons per H: the band energy densities over their mean photon
        # energies, which the band averaging wrote back into athinput.runtime
        hnu = band_energies(os.path.dirname(os.path.abspath(fn)))
        n_gamma = sum(u[b]/(hnu[b]*EV) for b in range(6))
        with np.errstate(divide='ignore'):
            f['logU'] = np.log10(n_gamma/f['nH'])
        u_pe = u[7]
        if 'Jdiff7' in d:
            u_pe = u_pe + 4.0*np.pi*d['Jdiff7'][k].astype(np.float64)/C_CODE*E_CODE
        f['chiPE'] = u_pe/U_DRAINE_PE
        f['chiH2'], f['chiCO'] = d['chi_H2'][k], d['chi_CO'][k]
    return f, d['x1f'], d['x2f'], float(d['Time'])


plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})


def frame(job):
    """Draw dump fn as OUTDIR/frame_<n>.png."""
    n, fn, outdir, molecular, radiation, title, scale, temp_log = job
    # the n_H and P/k ranges follow the ambient density (--density-scale)
    norm = dict(NORM, nH=LogNorm(10*scale, 1e3*scale), Pk=LogNorm(1e3*scale, 3e6*scale))
    cmap = dict(CMAP)
    if temp_log:
        norm['T'] = FuncNorm((_t_forward, _t_inverse), vmin=T_LOG[0], vmax=T_LOG[2])
        cmap['T'] = CMAP_T_LOG
    for c in cmap:
        cmap[c] = plt.get_cmap(cmap[c]).copy()
        cmap[c].set_bad(cmap[c](0.0))  # 0 on a log scale, or no light, takes the lowest
    own, shared = layout(molecular, radiation)
    LAYOUT = own + shared
    f, x1f, x2f, t = fields(fn)
    # grid: each column with its own colour bars is followed by a spacer that holds them
    ratios = [1, 0.45]*len(own) + [1]*len(shared)
    cols = [2*i for i in range(len(own))] + list(range(2*len(own), len(ratios)))
    # panels square: the figure widens with the columns (10 in for 6.45 units)
    fig, allax = plt.subplots(4, len(ratios), figsize=(10.0*sum(ratios)/6.45, 6.1),
                              sharex=True, sharey=True,
                              gridspec_kw={'width_ratios': ratios,
                                           'wspace': 0.0, 'hspace': 0.0})
    for i in range(len(ratios)):
        if i not in cols:
            for a in allax[:, i]:
                a.set_visible(False)
    axs = allax[:, cols]
    ims = {}
    for col, acol in zip(LAYOUT, axs.T):
        for q, a in zip(col, acol):
            if q is None:
                a.set_visible(False)
                continue
            kind = q if q in NORM else 'ion'
            ims[q if kind != 'ion' else 'ion'] = a.pcolormesh(
                x1f, x2f, f[q], norm=norm[kind], cmap=cmap[kind], shading='flat',
                rasterized=True)
            a.set_aspect('equal')
            a.text(0.05, 0.94, LAB[q], transform=a.transAxes, ha='left', va='top',
                   color='white', fontsize=9,
                   path_effects=[pe.withStroke(linewidth=2.0, foreground='black')])
            tick = np.round(2.0*x1f[-1]/3.0)   # two thirds of the half width
            a.set_xticks([-tick, 0, tick])
            a.set_yticks([-tick, 0, tick])
    for acol in axs.T:
        vis = [a for a in acol if a.get_visible()]
        vis[-1].tick_params(labelbottom=True)
        vis[-1].set_xlabel(r'$x$ [pc]', labelpad=1)
    for arow in axs:
        arow[0].set_ylabel(r'$y$ [pc]', labelpad=1)
    fig.subplots_adjust(left=0.045, right=0.92, bottom=0.075, top=0.955)
    fig.colorbar(ims['ion'], cax=fig.add_axes((0.93, 0.075, 0.012, 0.88)),
                 label='fraction of the element')
    for ic, col in enumerate(own):
        for row, q in enumerate(col):
            p = axs[row, ic].get_position()
            cax = fig.add_axes((p.x1 + 0.004, p.y0 + 0.004, 0.007, p.height - 0.008))
            cb = fig.colorbar(ims[q], cax=cax)
            if q == 'T' and temp_log:
                cb.set_ticks(T_TICKS)
                cb.set_ticklabels(['%g' % tk for tk in T_TICKS])
            cb.ax.tick_params(labelsize=6)
    fig.text(0.5, 0.985, title + r', $t = %.3f$ Myr' % (0.978*t), ha='center',
             va='top', fontsize=10)
    fig.savefig('%s/frame_%05d.png' % (outdir, n), dpi=150)
    plt.close(fig)
    print(os.path.basename(fn), 't = %.3f' % t, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('outdir')
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--frame', type=int, default=None)
    ap.add_argument('--molecular', action='store_true')
    ap.add_argument('--radiation', action='store_true')
    ap.add_argument('--temp-log', action='store_true',
                    help='T in the phase colour map (10-5000 K log, 5000-11000 K linear)')
    ap.add_argument('--density-scale', type=float, default=1.0,
                    help='ambient n_H / 100 cm^-3, which scales the n_H and P/k ranges')
    ap.add_argument('--title', default=r'expanding H II region, SB99 2 Myr, '
                    r'$Q = 10^{49}$ s$^{-1}$, $n_{\rm H} = 100$ cm$^{-3}$, no dust; '
                    r'Tigris ($N = 128^3$)')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    files = sorted(glob.glob(args.rundir + '/HII.out2.*.athdf'))
    jobs = [(n, fn, args.outdir, args.molecular, args.radiation, args.title,
             args.density_scale, args.temp_log)
            for n, fn in enumerate(files)]
    if args.frame is not None:
        frame(jobs[args.frame])
    else:
        with multiprocessing.Pool(args.nproc) as pool:
            pool.map(frame, jobs)
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '8', '-i',
                        args.outdir + '/frame_%05d.png', '-vf',
                        'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-pix_fmt', 'yuv420p',
                        '-vcodec', 'libx264', args.outdir + '/slices.mp4'], check=True)
        print(args.outdir + '/slices.mp4')
