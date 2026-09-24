"""z = 0 slice movie of hii_dtype with GOW17 + O2,S3,N2, 7 bands.

usage: python slice_movie.py RUNDIR OUTDIR
Writes OUTDIR/frame_NNNNN.png, one per HII.out2 dump, then OUTDIR/slices.mp4.
Every panel keeps one colour range over all frames. The text panel carries t,
the cycle, the hydro step at that time, r_sh, and the ionizing-photon budget.
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
MH, KB, MU_H, X_HE = 1.6735575e-24, 1.380649e-16, 1.4, 0.1
TOT = {'He': 0.1, 'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5, 'C': 1.6e-4}
CORE = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+']
IONS = ['O++', 'S+', 'S++', 'S3+', 'N+', 'N++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1,
          'Si+': 1, 'O++': 2, 'S+': 1, 'S++': 2, 'S3+': 3, 'N+': 1, 'N++': 2}

# (title, function of the fields, norm, cmap)
FRAC = dict(norm=LogNorm(1e-3, 1.0), cmap='viridis')
PANELS = [
    (r'$n_{\rm H}$ [cm$^{-3}$]', lambda f: f['nH'], LogNorm(5, 1e3), 'cividis'),
    ('T [K]', lambda f: f['T'], LogNorm(50, 2e4), 'inferno'),
    (r'$P/k$ [K cm$^{-3}$]', lambda f: f['Pk'], LogNorm(1e4, 3e6), 'magma'),
    (r'$|v|$ [km s$^{-1}$]', lambda f: f['v'], Normalize(0, 15), 'plasma'),
    (r'$x_{\rm H^+}$', lambda f: f['H+'], FRAC['norm'], FRAC['cmap']),
    (r'$2x_{\rm H_2}$', lambda f: 2*f['H2'], FRAC['norm'], FRAC['cmap']),
    (r'He$^+$/He', lambda f: f['He+']/TOT['He'], FRAC['norm'], FRAC['cmap']),
    (r'C$^+$/C', lambda f: f['C+']/TOT['C'], FRAC['norm'], FRAC['cmap']),
    (r'O$^+$/O', lambda f: f['O+']/TOT['O'], FRAC['norm'], FRAC['cmap']),
    (r'O$^{2+}$/O', lambda f: f['O++']/TOT['O'], FRAC['norm'], FRAC['cmap']),
    (r'N$^+$/N', lambda f: f['N+']/TOT['N'], FRAC['norm'], FRAC['cmap']),
    (r'N$^{2+}$/N', lambda f: f['N++']/TOT['N'], FRAC['norm'], FRAC['cmap']),
    (r'S$^+$/S', lambda f: f['S+']/TOT['S'], FRAC['norm'], FRAC['cmap']),
    (r'S$^{2+}$/S', lambda f: f['S++']/TOT['S'], FRAC['norm'], FRAC['cmap']),
    (r'S$^{3+}$/S', lambda f: f['S3+']/TOT['S'], FRAC['norm'], FRAC['cmap']),
]

h = athena_read.hst(rundir + '/HII.hst')
m_sh = h['sh_mass'] > 0
t_sh, r_sh = h['time'][m_sh], h['sh_mass_r'][m_sh]/h['sh_mass'][m_sh]
NION = 5
l_tot = sum(h['Ltot%d' % b] for b in range(NION))
l_gas = sum(h['Labs%d' % b] - h['Ldust%d' % b] for b in range(NION))
l_dust = sum(h['Ldust%d' % b] for b in range(NION))
l_esc = sum(h['Lesc%d' % b] for b in range(NION))


def fields(fn):
    d = athena_read.athdf(fn)
    k = d['x3v'].size//2
    names = CORE + IONS
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(names))]
    f = {nm: d[kk][k].astype(np.float64) for nm, kk in zip(names, keys)}
    f['nH'] = d['rho'][k]/MU_H
    xe = sum(CHARGE.get(nm, 0)*f[nm] for nm in names)
    f['Pk'] = d['press'][k]*MH*1e10/KB
    f['T'] = f['Pk']/(f['nH']*(1.0 - f['H2'] + X_HE + xe))
    f['v'] = np.sqrt(d['vel1'][k]**2 + d['vel2'][k]**2 + d['vel3'][k]**2)
    ext = [d['x1f'][0], d['x1f'][-1], d['x2f'][0], d['x2f'][-1]]
    return f, ext, d['Time'], d['NumCycles']


files = sorted(glob.glob(rundir + '/HII.out2.*.athdf'))
for n, fn in enumerate(files):
    f, ext, t, ncyc = fields(fn)
    fig, axs = plt.subplots(4, 4, figsize=(17, 16))
    axs = axs.ravel()
    for ax, (title, fun, norm, cmap) in zip(axs, PANELS):
        im = ax.imshow(np.maximum(fun(f), 1e-30), origin='lower', extent=ext, norm=norm,
                       cmap=cmap, interpolation='nearest')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        ax.set_title(title, fontsize=11)
        ax.set_xticks([-15, 0, 15])
        ax.set_yticks([-15, 0, 15])
    i = np.argmin(np.abs(h['time'] - t))
    rs = np.interp(t, t_sh, r_sh) if t >= t_sh[0] else 0.0
    lt = l_tot[i] if l_tot[i] > 0 else np.nan
    axs[-1].axis('off')
    axs[-1].text(0.0, 1.0,
                 'hii_dtype, GOW17 + O2,S3,N2\n7 bands (5 ionizing), 64$^3$\n'
                 'z = 0 slice, x and y in pc\n\n'
                 't = %.3f code (%.3f Myr)\ncycle %d\ndt = %.2e code\n'
                 'r_sh = %.2f pc\n\nionizing photons at this step:\n'
                 '  gas %.3f  dust %.3f  esc %.3f'
                 % (t, t*0.978, ncyc, h['dt'][i], rs, l_gas[i]/lt, l_dust[i]/lt,
                    l_esc[i]/lt),
                 va='top', family='monospace', fontsize=11, transform=axs[-1].transAxes)
    fig.tight_layout()
    fig.savefig('%s/frame_%05d.png' % (outdir, n), dpi=70)
    plt.close(fig)
    print(fn.split('/')[-1], 't = %.3f' % t, flush=True)

subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '6', '-i',
                outdir + '/frame_%05d.png', '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                '-pix_fmt', 'yuv420p', '-vcodec', 'libx264', outdir + '/slices.mp4'],
               check=True)
print(outdir + '/slices.mp4')
