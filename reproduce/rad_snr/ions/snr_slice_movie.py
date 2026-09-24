"""z = 0 slice movie of the radiative SNR with GOW17 + multi-ion (C, Si, O, S, N higher ions).

usage: python snr_slice_movie.py RUNDIR OUTDIR
Writes OUTDIR/frame_NNNNN.png, one per snr_ions.out2 dump, then OUTDIR/snr_slices.mp4.
Every panel keeps one colour range over all frames. The text panel carries t,
the cycle, the hydro step, and the hot-gas mass and radial momentum from the .hst.

Scalar mapping (hdf5 names): rmetal, rSN, rret are feedback; the 22 photochem
species follow in rHI, r4, ..., r24 in this order:
  rHI=He+, r4=OHx, r5=CHx, r6=CO, r7=C+, r8=HCO+, r9=H2, r10=H+, r11=H3+, r12=H2+,
  r13=O+, r14=Si+, r15=C_high, r16=Si_high, r17=O++, r18=O_high, r19=S+, r20=S++,
  r21=S_high, r22=N+, r23=N++, r24=N_high.
Verified at t = 0 (ambient: rHI 1.4e-3 <= 0.1, r7 = 1.6e-4 = C total, r9 = 3e-8,
r14 = 1.7e-6 = Si total, r19 = 1.45e-5 = S total) and at t = 0.05, where the
per-cell maxima of the element sums C, O, S, N, Si divided by their totals are
1.0000001 for every element.
"""
import glob
import os
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
sys.path.insert(0, '/Users/jgkim/Projects/tigris-gow17/vis/python')
import athena_read

rundir, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
# ism units: density unit m_H/cm^3, rho = 1.4 n_H; pressure unit m_H (km/s)^2
MH, KB, MU_H, X_HE, MSUN = 1.6735575e-24, 1.380649e-16, 1.4, 0.1, 1.98841e33
TOT = {'He': 0.1, 'C': 1.6e-4, 'Si': 1.7e-6, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++', 'N_high']
KEYS = ['rHI'] + ['r%d' % n for n in range(4, 4 + len(NAMES) - 1)]
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3}

units = athena_read.athinput(rundir + '/athinput.runtime')['units']
M_UNIT = float(units['mass_cgs'])/MSUN


def frac(title, sp, el):
    return (title, lambda f: f[sp]/TOT[el], LogNorm(1e-3, 1.0), 'viridis')


PANELS = [
    (r'$n_{\rm H}$ [cm$^{-3}$]', lambda f: f['nH'], LogNorm(3e-2, 10), 'cividis'),
    ('T [K]', lambda f: f['T'], LogNorm(1e2, 1e8), 'inferno'),
    (r'$P/k$ [K cm$^{-3}$]', lambda f: f['Pk'], LogNorm(1e3, 1e8), 'magma'),
    (r'$|v|$ [km s$^{-1}$]', lambda f: f['v'], Normalize(0, 500), 'plasma'),
    (r'$x_{\rm H^+}$', lambda f: f['H+'], LogNorm(1e-3, 1.0), 'viridis'),
    (r'$2x_{\rm H_2}$', lambda f: 2*f['H2'], LogNorm(1e-3, 1.0), 'viridis'),
    frac(r'He$^+$/He', 'He+', 'He'),
    frac(r'C$^+$/C', 'C+', 'C'),
    frac(r'C$_{\rm high}$/C', 'C_high', 'C'),
    frac(r'Si$^+$/Si', 'Si+', 'Si'),
    frac(r'Si$_{\rm high}$/Si', 'Si_high', 'Si'),
    frac(r'O$^+$/O', 'O+', 'O'),
    frac(r'O$^{2+}$/O', 'O++', 'O'),
    frac(r'O$_{\rm high}$/O', 'O_high', 'O'),
    frac(r'S$^+$/S', 'S+', 'S'),
    frac(r'S$^{2+}$/S', 'S++', 'S'),
    frac(r'S$_{\rm high}$/S', 'S_high', 'S'),
    frac(r'N$^+$/N', 'N+', 'N'),
    frac(r'N$^{2+}$/N', 'N++', 'N'),
    frac(r'N$_{\rm high}$/N', 'N_high', 'N'),
]

h = athena_read.hst(rundir + '/snr_ions.hst')


def fields(fn):
    d = athena_read.athdf(fn)
    k = d['x3v'].size//2
    f = {nm: d[kk][k].astype(np.float64) for nm, kk in zip(NAMES, KEYS)}
    f['nH'] = d['rho'][k]/MU_H
    xe = sum(CHARGE.get(nm, 0)*f[nm] for nm in NAMES)
    f['Pk'] = d['press'][k]*MH*1e10/KB
    f['T'] = f['Pk']/(f['nH']*(1.0 - f['H2'] + X_HE + xe))
    f['v'] = np.sqrt(d['vel1'][k]**2 + d['vel2'][k]**2 + d['vel3'][k]**2)
    ext = [d['x1f'][0], d['x1f'][-1], d['x2f'][0], d['x2f'][-1]]
    return f, ext, d['Time'], d['NumCycles']


files = sorted(glob.glob(rundir + '/snr_ions.out2.*.athdf'))
for n, fn in enumerate(files):
    f, ext, t, ncyc = fields(fn)
    fig, axs = plt.subplots(3, 7, figsize=(28, 12))
    axs = axs.ravel()
    for ax, (title, fun, norm, cmap) in zip(axs, PANELS):
        im = ax.imshow(np.maximum(fun(f), 1e-30), origin='lower', extent=ext, norm=norm,
                       cmap=cmap, interpolation='nearest')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        ax.set_title(title, fontsize=12)
        ax.set_xticks([-30, 0, 30])
        ax.set_yticks([-30, 0, 30])
    i = np.argmin(np.abs(h['time'] - t))
    axs[-1].axis('off')
    axs[-1].text(0.0, 1.0,
                 'radiative SNR\nGOW17 + multi-ion, 64$^3$\n'
                 'z = 0 slice, x and y in pc\n\n'
                 't = %.4f code\n  = %.4f Myr\ncycle %d\ndt = %.2e code\n\n'
                 'M_hot = %.3g Msun\np_r = %.3g Msun km/s'
                 % (t, t*0.978, ncyc, h['dt'][i], h['M_hot'][i]*M_UNIT, h['pr'][i]*M_UNIT),
                 va='top', family='monospace', fontsize=13, transform=axs[-1].transAxes)
    fig.tight_layout()
    fig.savefig('%s/frame_%05d.png' % (outdir, n), dpi=70)
    plt.close(fig)
    print(fn.split('/')[-1], 't = %.4f' % t, flush=True)

subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '8', '-i',
                outdir + '/frame_%05d.png', '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                '-pix_fmt', 'yuv420p', '-vcodec', 'libx264', outdir + '/snr_slices.mp4'],
               check=True)
print(outdir + '/snr_slices.mp4')
