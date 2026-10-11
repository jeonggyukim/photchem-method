"""Diagnostic figure of one hii_dtype run with GOW17 (and ion ladders if built).

usage: python hii_dtype_diag.py RUNDIR [SNAP] [OUT.png]
       python hii_dtype_diag.py RUNDIR movie [OUT.mp4]

Every slice uses one colour range for all frames, so a movie's frames compare.

Row 1: history -- shell radius (with the NCR reference), timestep, masses, photon
budget per band, radial momentum, ionization-front radius; a line marks the frame.
Row 2: midplane slices of n_H, T, |v|, P/k and the band energy densities.
Rows 3-4: midplane slices of every species, as a fraction of its element.
Row 5: density-weighted spherical averages of T, n_H, the ion fractions and v_r.
"""
import glob
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
athena_read = paths.athena_read()

# ISM units: rho = n_H [cm^-3], pressure unit mu_H cm^-3 (km/s)^2, time unit pc/(km/s)
MUH, KB = 2.34335276e-24, 1.380649e-16
P_TO_PK = MUH*1e10/KB
TUNIT_MYR = 3.0856776e13/3.15576e13
X_HE, X_C, X_O, X_S, X_SI = 0.1, 1.6e-4, 3.2e-4, 1.45e-5, 1.7e-6
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2}
REF = paths.tigris('tst/regression/data/ref_rayt_solutions/hii_dtype_ncr.txt')
BANDS = ['LyC', 'LW', 'PE']


def global_er_max(files):
    """log10 of the largest band energy density over all frames, per band."""
    er_max = {}
    for b in range(3):
        k = 'Er_rayt%d' % b
        vals = [np.max(athena_read.athdf(f, quantities=[k])[k]) for f in files]
        er_max[k] = np.log10(max(max(vals), 1e-300))
    return er_max


def make_frame(fn, h, er_max, out):
    d = athena_read.athdf(fn)
    nscal = sum(1 for k in d if k == 'rHI' or (k.startswith('r') and k[1:].isdigit()))
    names = CORE + (IONS if nscal > len(CORE) else [])
    keys = ['rHI'] + ['r%d' % n for n in range(1, nscal)]
    x = {nm: d[k] for nm, k in zip(names, keys)}

    nH = d['rho']
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in names)
    T = d['press']*P_TO_PK/(nH*(1.0 - x['H2'] + X_HE + xe))
    vel = np.sqrt(d['vel1']**2 + d['vel2']**2 + d['vel3']**2)
    x_hi = np.clip(1.0 - x['H+'] - 2*x['H2'] - 3*x['H3+'] - 2*x['H2+'], 0, 1)

    kz = nH.shape[0]//2
    X, Y = np.meshgrid(d['x1f'], d['x2f'])
    tcode = d['Time']
    t = h['time']

    fig = plt.figure(figsize=(28, 19))
    gs = fig.add_gridspec(5, 7, hspace=0.42, wspace=0.38, top=0.95)

    def hist_axis(pos, title, ylabel, log=False):
        ax = fig.add_subplot(pos)
        ax.axvline(tcode, color='0.6', lw=0.8)
        ax.set(xlabel='t [code]', ylabel=ylabel, title=title, xlim=(0, t.max()*1.02))
        if log:
            ax.set_yscale('log')
        return ax

    # --- Row 1: history --------------------------------------------------------------
    ax = hist_axis(gs[0, 0], 'shell radius', 'r_sh [pc]')
    m = h['sh_mass'] > 0
    ax.plot(t[m], h['sh_mass_r'][m]/h['sh_mass'][m], 'k-', label='GOW17')
    dd = np.genfromtxt(REF, dtype=None, names=True, encoding=None)
    ax.plot(dd[dd.dtype.names[0]], dd[dd.dtype.names[1]], 'r--', label='NCR ref')
    ax.legend(fontsize=8)

    ax = hist_axis(gs[0, 1], 'timestep', 'dt [yr]', log=True)
    ax.plot(t, h['dt']*TUNIT_MYR*1e6, 'k-')

    ax = hist_axis(gs[0, 2], 'mass', 'mass [code]', log=True)
    ax.plot(t, h['mass_ion'], label='M_ion')
    ax.plot(t[m], h['sh_mass'][m], label='M_sh')
    ax.legend(fontsize=8)

    ax = hist_axis(gs[0, 3], 'photon budget per band', 'L / L_tot')
    for b, c in zip(range(3), 'bgr'):
        L = h['Ltot%d' % b]
        ok = L > 0
        if not ok.any():
            continue
        ax.plot(t[ok], (h['Labs%d' % b] - h['Ldust%d' % b])[ok]/L[ok], c + '-',
                label='gas ' + BANDS[b])
        ax.plot(t[ok], h['Ldust%d' % b][ok]/L[ok], c + '--', label='dust ' + BANDS[b])
        ax.plot(t[ok], h['Lesc%d' % b][ok]/L[ok], c + ':', label='esc ' + BANDS[b])
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, ncol=3)

    ax = hist_axis(gs[0, 4], 'radial momentum', 'p_r [code]')
    ax.plot(t, h['pr_ion'], label='ionized')
    ax.plot(t, h['pr_neu'], label='neutral')
    ax.legend(fontsize=8)

    ax = hist_axis(gs[0, 5], 'ionization-front radius', '<r>_IF [pc]')
    ax.plot(t, h['IF_vol_r']/np.maximum(h['IF_vol'], 1e-30), 'k-')

    ax = fig.add_subplot(gs[0, 6])
    ax.axis('off')
    ax.text(0, 1, '%s\nt = %.3f code\n  = %.3f Myr\ncells %s\nspecies %d (%s)'
            % (os.path.basename(fn), tcode, tcode*TUNIT_MYR, str(nH.shape), nscal,
               'core + O2,S2' if nscal > len(CORE) else 'core'),
            va='top', fontsize=12, family='monospace')

    def slice_panel(pos, f, title, vmin, vmax, cmap='viridis', log=True):
        ax = fig.add_subplot(pos)
        z = f[kz].astype(np.float64)
        if log:
            z = np.log10(np.maximum(z, 1e-300))
            title = 'log ' + title
        im = ax.pcolormesh(X, Y, z, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=10)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

    # --- Row 2: gas and radiation -------------------------------------------------------
    slice_panel(gs[1, 0], nH, 'n_H [cm^-3]', 0, 3.5)
    slice_panel(gs[1, 1], T, 'T [K]', 1, 4.3, cmap='inferno')
    slice_panel(gs[1, 2], vel, '|v| [km/s]', 0, 15, cmap='magma', log=False)
    slice_panel(gs[1, 3], d['press']*P_TO_PK, 'P/k [K cm^-3]', 3.5, 7, cmap='cividis')
    for b in range(3):
        k = 'Er_rayt%d' % b
        title = 'cE %s [code]' % BANDS[b] + ('  (zero in all frames)' if er_max[k] < -299
                                               else '')
        vmax = max(er_max[k], -290.0)
        slice_panel(gs[1, 4 + b], d[k], title, vmax - 10, vmax, cmap='plasma')

    # --- Rows 3-4: species -------------------------------------------------------------
    spec = [('H I', x_hi, 1.0), ('H+', x['H+'], 1.0), ('H2', 2*x['H2'], 1.0),
            ('He+', x['He+'], X_HE), ('C+', x['C+'], X_C), ('CO', x['CO'], X_C),
            ('e', xe, 1.0), ('O+', x['O+'], X_O), ('Si+', x['Si+'], X_SI),
            ('OHx', x['OHx'], X_O), ('CHx', x['CHx'], X_C)]
    if 'O++' in x:
        spec += [('O++', x['O++'], X_O), ('S+', x['S+'], X_S), ('S++', x['S++'], X_S)]
    for n, (nm, f, tot) in enumerate(spec):
        title = 'x_e' if nm == 'e' else '%s / element' % nm
        slice_panel(gs[2 + n//7, n % 7], f/tot, title, -4, 0)

    # --- Row 5: radial profiles ---------------------------------------------------------
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    rb = np.linspace(0, d['x1f'][-1], 49)
    ri = np.digitize(R.ravel(), rb)

    def radial(f, w):
        num = np.bincount(ri, (f*w).ravel(), len(rb) + 1)[1:len(rb)]
        den = np.bincount(ri, w.ravel(), len(rb) + 1)[1:len(rb)]
        prof = np.where(den > 0, num/np.maximum(den, 1e-300), np.nan)
        return 0.5*(rb[1:] + rb[:-1]), prof

    def prof_axis(pos, title, ylabel=None):
        ax = fig.add_subplot(pos)
        ax.set(xlabel='r [pc]', title=title, xlim=(0, rb[-1]))
        if ylabel:
            ax.set_ylabel(ylabel)
        return ax

    ones = np.ones_like(nH)
    ax = prof_axis(gs[4, 0], 'temperature', '<T>_n [K]')
    ax.semilogy(*radial(T, nH), 'k-')
    ax.set_ylim(10, 2e4)
    ax = prof_axis(gs[4, 1], 'density', '<n_H> [cm^-3]')
    ax.semilogy(*radial(nH, ones), 'k-')
    ax.set_ylim(1, 3e3)
    groups = [('H, He', ['H I', 'H+', 'H2', 'He+']), ('C, Si', ['C+', 'CO', 'Si+']),
              ('O', ['O+', 'O++']), ('S', ['S+', 'S++'])]
    lookup = {s[0]: s for s in spec}
    for n, (title, members) in enumerate(groups):
        present = [lookup[nm] for nm in members if nm in lookup]
        if not present:
            continue
        ax = prof_axis(gs[4, 2 + n], title + ' fractions')
        for nm, f, tot in present:
            ax.semilogy(*radial(f/tot, nH), label=nm)
        ax.set_ylim(1e-5, 1.5)
        ax.legend(fontsize=8)
    ax = prof_axis(gs[4, 6], 'radial velocity', '<v_r>_n [km/s]')
    ax.plot(*radial((d['vel1']*X1 + d['vel2']*X2 + d['vel3']*X3)/np.maximum(R, 1e-30),
                    nH), 'k-')
    ax.set_ylim(-2, 15)

    fig.suptitle('hii_dtype GOW17   t = %.3f Myr' % (tcode*TUNIT_MYR), fontsize=15, y=0.98)
    fig.savefig(out, dpi=70, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    rundir = sys.argv[1].rstrip('/')
    files = sorted(glob.glob(rundir + '/HII.out2.*.athdf'))
    h = athena_read.hst(rundir + '/HII.hst')
    er_max = global_er_max(files)
    if len(sys.argv) > 2 and sys.argv[2] == 'movie':
        out = sys.argv[3] if len(sys.argv) > 3 else rundir + '/hii_dtype_diag.mp4'
        fdir = rundir + '/diag_frames'
        os.makedirs(fdir, exist_ok=True)
        for n, fn in enumerate(files):
            make_frame(fn, h, er_max, '%s/frame%04d.png' % (fdir, n))
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '4',
                        '-i', fdir + '/frame%04d.png', '-vf',
                        'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-pix_fmt', 'yuv420p',
                        out], check=True)
    else:
        snap = int(sys.argv[2]) if len(sys.argv) > 2 else -1
        out = sys.argv[3] if len(sys.argv) > 3 else rundir + '/hii_dtype_diag.png'
        make_frame(files[snap], h, er_max, out)
    print(out)
