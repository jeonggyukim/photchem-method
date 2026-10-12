"""Static H II region (5e4 K blackbody, Q(H) = 1e49 photons s^-1, uniform dust-free gas
of n_H = 100 cm^-3): heating and cooling by process, Tigris against Cloudy 25.

Top row: at two cells, r = 1.14 and 2.44 pc (0.39 and 0.83 R_s), the heating and
cooling rates per n_H^2 as functions of T with the cell's abundances and radiation held
fixed (the probes of ../fig_ion_stages/athinput.stromgren_stages). The crossing of total heating and cooling is the
equilibrium temperature these abundances predict; vertical lines mark the temperature
of the Tigris cell and of Cloudy at the same radius. Middle row: radial profiles of
the heating and cooling by agent. Bottom row: the mean heat per H photoionization and
the share of H photoionizations by recombination photons, for recombination photons
carried by the diffuse ray tracer (no OTS) and absorbed on the spot (OTS).

Tigris: GOW17 + O3,S3,N3 with He++, 10 bands, 64^3, 20 point rays per cell; Cloudy:
`diffuse OTS` and its default transfer of the diffuse field. Reads
$PHOTCHEM_RUNS/T6_multi_ion/hii_static_bands10/bb5e4_{ots,no_ots}
(../fig_ion_stages/run_tigris.sh) and
$PHOTCHEM_RUNS/T6_multi_ion/cloudy_bb5e4/{ots,no_ots}; the profiles come from
reduced_bb5e4.txt when it exists, and `python plot.py --from-runs` rewrites it.
Writes ../../../figures/heat_cool_bb5e4.png."""
import glob
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from photchem_method import paths  # noqa: E402
from photchem_method import reduced  # noqa: E402
from photchem_method.paths import figure_formats  # noqa: E402

SPEC = 'bb5e4'
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '..', 'figures', 'heat_cool_' + SPEC)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
REDUCED = os.path.join(HERE, 'reduced_%s.txt' % SPEC)
MODES = [('no_ots', 'no OTS'), ('ots', 'OTS')]
NH, PC, EV = 100.0, 3.0856776e18, 1.602176634e-12
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++',
         'N_high', 'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
# cooling agents: Tigris slots, Cloudy `save cooling each` columns
COOL = [('H', ['Lya', 'HII_rr', 'HI_ci'], ['H']),
        ('free-free', ['HII_ff', 'He2_ff'], ['FF_H', 'FF_M']),
        ('He', ['He2_rr', 'HeI_ci', 'HeII_ci'], ['He']),
        ('O', ['OI', 'O+', 'O++', 'O_high'], ['O']),
        ('N', ['N+', 'N++', 'N_high'], ['N']),
        ('S', ['S+', 'S++', 'S_high'], ['S']),
        ('C', ['CII', 'CI', 'C_high'], ['C'])]
# heating agents: Tigris slots, Cloudy `save heating` labels
HEAT = [('HI', ['HI', 'He_ots'], ['H  1']), ('HeI', ['HeI'], ['He 1']),
        ('HeII', ['HeII'], ['He 2']), ('metals', [], None)]
HEAT_LABEL = {'HI': 'H I', 'HeI': 'He I', 'HeII': 'He II', 'metals': 'metals (Cloudy)'}
CCOL = {'H': '#0072B2', 'free-free': '#56B4E9', 'He': '#009E73', 'O': '#D55E00',
        'N': '#CC79A7', 'S': '#E69F00', 'C': '#999999', 'HI': '#0072B2',
        'HeI': '#009E73', 'HeII': '#000000', 'metals': '#D55E00'}


def probe_file(path):
    with open(path) as f:
        f.readline()
        pos = f.readline().split()
        names = f.readline().split()[1:]
    a = np.loadtxt(path)
    x = np.array([float(v) for v in pos[3:6]])
    return float(np.sqrt((x**2).sum())), {n: a[:, i] for i, n in enumerate(names)}


def alpha_hg(T, case):
    """H+ recombination coefficient [cm^3 s^-1], Hui & Gnedin (1997)."""
    lam = 2.0*157807.0/T
    if case == 'A':
        return 1.269e-13*lam**1.503/(1.0 + (lam/0.522)**0.470)**1.923
    return 2.753e-14*lam**1.500/(1.0 + (lam/2.740)**0.407)**2.242


def tigris(rundir):
    d = athena_read.athdf(sorted(glob.glob(rundir + '/*.out2.*.athdf'))[-1])
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64) for nm, k in zip(NAMES, keys)}
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    xhi = np.maximum(1.0 - x['H+'] - 2*x['H2'], 0.0)
    T = d['temp'].astype(np.float64)
    f = {'T': T}
    for name, slots, _ in COOL:
        f['cool_' + name] = sum(d['cool_' + s] for s in slots)/NH
    for name, slots, _ in HEAT:
        if slots:
            f['heat_' + name] = sum(d['heat_' + s] for s in slots)/NH
    f['cool_total'] = d['cool_rate']/NH
    f['heat_total'] = d['heat_rate']/NH
    # H photoionizations per volume and their heat: stellar, recombination photons
    # (diffuse ray tracer), on the spot (case A minus case B recombinations, no heat)
    f['ion_star'] = NH*xhi*d['xi_HI_star']
    f['ion_rec'] = NH*xhi*d['xi_HI_rec']
    f['ion_ots'] = (NH*xe)*(NH*x['H+'])*(alpha_hg(T, 'A') - alpha_hg(T, 'B')) \
        if 'ots' == os.path.basename(rundir).split('_', 1)[1] else 0.0*T
    f['q_rec'] = NH*d['heat_HI_rec']
    f['q_star'] = NH*d['heat_HI'] - f['q_rec']
    X3, X2, X1 = np.meshgrid(d['x3v'], d['x2v'], d['x1v'], indexing='ij')
    R = np.sqrt(X1**2 + X2**2 + X3**2)
    dx = d['x1f'][1] - d['x1f'][0]
    rb = np.arange(0, d['x1f'][-1] + 1e-9, 0.5*dx)
    ri = np.digitize(R.ravel(), rb)
    cnt = np.bincount(ri, minlength=len(rb) + 1)[1:len(rb)]
    full = cnt > 0
    prof = {k: (np.bincount(ri, v.ravel(), len(rb) + 1)[1:len(rb)]
                / np.maximum(cnt, 1))[full] for k, v in f.items()}
    return 0.5*(rb[1:] + rb[:-1])[full], prof


def cloudy(mdir):
    r = np.loadtxt(os.path.join(mdir, 'hii.rad'))[:, 1]/PC
    p = {}
    with open(os.path.join(mdir, 'hii.cool_each')) as f:
        cols = [c.strip() for c in f.readline().lstrip('#').split('\t')]
    ce = np.loadtxt(os.path.join(mdir, 'hii.cool_each'))
    for name, _, ccols in COOL:
        p['cool_' + name] = sum(ce[:, cols.index(c)] for c in ccols)/NH**2
    p['cool_total'] = ce[:, 2]/NH**2
    p['T'] = ce[:, 1]
    tot, frac = [], {name: [] for name, _, _ in HEAT}
    with open(os.path.join(mdir, 'hii.heat')) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tok = line.rstrip('\n').split('\t')
            tot.append(float(tok[2]))
            pairs = {tok[i].strip(): float(tok[i + 1]) for i in range(5, len(tok) - 1, 2)}
            for name, _, labs in HEAT:
                if labs is None:
                    continue
                frac[name].append(sum(pairs.get(lab.strip(), 0.0) for lab in labs))
            metals = sum(v for k, v in pairs.items()
                         if k.split()[0] not in ('H', 'He', 'Comp'))
            frac['metals'].append(metals)
    tot = np.array(tot)
    p['heat_total'] = tot/NH**2
    for name in frac:
        p['heat_' + name] = np.array(frac[name])*tot/NH**2
    return r, p


if '--from-runs' in sys.argv or not os.path.exists(REDUCED):
    athena_read = paths.athena_read()
    arrays, sources = {}, []
    for mode, _ in MODES:
        tdir = paths.runs('T6_multi_ion', 'hii_static_bands10', '%s_%s' % (SPEC, mode))
        cdir = paths.runs('T6_multi_ion', 'cloudy_' + SPEC, mode)
        for code, (r, p) in [('tigris', tigris(tdir)), ('cloudy', cloudy(cdir))]:
            arrays['%s_%s_r' % (code, mode)] = r
            arrays.update({'%s_%s_%s' % (code, mode, q): v for q, v in p.items()})
        for n in range(2):
            rp, pr = probe_file(os.path.join(tdir, 'cloud.probe%d.txt' % n))
            arrays['probe%d_%s_r' % (n, mode)] = np.array([rp])
            arrays.update({'probe%d_%s_%s' % (n, mode, k): v for k, v in pr.items()})
        sources += [tdir, cdir]
    reduced.save(REDUCED, 'heating and cooling of the static H II region, spectrum %s\n'
                 '<code>_<mode>_<quantity>: radial profiles, r [pc], T [K], heat_/cool_ '
                 'per n_H^2 [erg cm^3 s^-1], ion_ per volume [cm^-3 s^-1], q_ [erg cm^-3 '
                 's^-1]\nprobe<n>_<mode>_<column>: the probe files, per H [erg s^-1]'
                 % SPEC, arrays, sources=sources)
data = reduced.load(REDUCED)

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
                     'font.size': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
                     'legend.fontsize': 7.5})
fig, ax = plt.subplots(3, 2, figsize=(7.2, 8.6))


def at(r, rr, v):
    return float(np.interp(r, rr, v))


# top row: the T sweeps
for n, a in enumerate(ax[0]):
    g = {k[len('probe%d_no_ots_' % n):]: v for k, v in data.items()
         if k.startswith('probe%d_no_ots_' % n)}
    go = {k[len('probe%d_ots_' % n):]: v for k, v in data.items()
          if k.startswith('probe%d_ots_' % n)}
    rp = float(g['r'][0])
    T = g['T'][1:]
    for name, slots, _ in COOL:
        y = sum(g['cool_' + s][1:] for s in slots)/NH
        a.loglog(T, y, color=CCOL[name], lw=1.0, alpha=0.85, label=name)
    a.loglog(T, g['cool'][1:]/NH, 'k-', lw=2.0, label=r'total $\Lambda$')
    a.loglog(T, g['heat'][1:]/NH, 'k--', lw=2.0, label=r'total $\Gamma$')
    a.loglog(go['T'][1:], go['cool'][1:]/NH, '-', color='0.55', lw=1.2)
    a.loglog(go['T'][1:], go['heat'][1:]/NH, '--', color='0.55', lw=1.2)
    a.axvline(g['T'][0], color='k', lw=0.8, ls=':')
    a.axvline(go['T'][0], color='0.55', lw=0.8, ls=':')
    tc = at(rp, data['cloudy_no_ots_r'], data['cloudy_no_ots_T'])
    a.axvline(tc, color='#D55E00', lw=1.2, alpha=0.6)
    a.plot(tc, at(rp, data['cloudy_no_ots_r'], data['cloudy_no_ots_cool_total']), 'o',
           color='#D55E00', ms=4, alpha=0.7)
    a.set(xlim=(3e3, 3e4), ylim=(1e-27, 3e-23), xlabel=r'$T$ [K]',
          title=r'$r = %.2f$ pc' % rp)
    if n == 0:
        a.set_ylabel(r'rate$/n_{\rm H}^2$ [erg cm$^3$ s$^{-1}$]')
        a.legend(loc='lower right', ncol=2, fontsize=6.5)
    print('probe %d r = %.2f pc: T no OTS %.0f, OTS %.0f, Cloudy no OTS %.0f K'
          % (n, rp, g['T'][0], go['T'][0], tc))

# middle row: radial heating and cooling by agent, no OTS
rt, rc = data['tigris_no_ots_r'], data['cloudy_no_ots_r']
rmax = 3.4
for name, _, _ in HEAT:
    if 'tigris_no_ots_heat_' + name in data:
        ax[1, 0].semilogy(rt, data['tigris_no_ots_heat_' + name], color=CCOL[name],
                          lw=1.4, label=HEAT_LABEL[name])
    ax[1, 0].semilogy(rc, data['cloudy_no_ots_heat_' + name], color=CCOL[name], lw=1.4,
                      ls='--', alpha=0.6,
                      label=HEAT_LABEL[name] if name == 'metals' else None)
for name, _, _ in COOL:
    ax[1, 1].semilogy(rt, data['tigris_no_ots_cool_' + name], color=CCOL[name], lw=1.4,
                      label=name)
    ax[1, 1].semilogy(rc, data['cloudy_no_ots_cool_' + name], color=CCOL[name], lw=1.4,
                      ls='--', alpha=0.6)
for a, t in zip(ax[1], ['heating', 'cooling']):
    a.set(xlim=(0, rmax), ylim=(1e-26, 1e-23), xlabel=r'$r$ [pc]', title=t)
    for n in range(2):
        a.axvline(float(data['probe%d_no_ots_r' % n][0]), color='k', lw=0.6, ls=':')
ax[1, 0].set_ylabel(r'rate$/n_{\rm H}^2$ [erg cm$^3$ s$^{-1}$]')
ax[1, 0].legend(loc='lower left', fontsize=7)
ax[1, 1].legend(loc='lower left', ncol=2, fontsize=7)

# bottom row: heat per H photoionization and the recombination-photon share
COL = {'no_ots': '#D55E00', 'ots': '#0072B2'}
for mode, label in MODES:
    r = data['tigris_%s_r' % mode]
    s, rec, ots = (data['tigris_%s_ion_%s' % (mode, k)] for k in ('star', 'rec', 'ots'))
    qs, qr = data['tigris_%s_q_star' % mode], data['tigris_%s_q_rec' % mode]
    ok = s > 0
    ax[2, 0].plot(r[ok], (qs + qr)[ok]/(s + rec + ots)[ok]/EV, color=COL[mode], lw=1.6,
                  label='all photons, ' + label)
    if mode == 'no_ots':
        ax[2, 0].plot(r[ok], qs[ok]/s[ok]/EV, color=COL[mode], lw=1.0, ls='--',
                      label='stellar photons, ' + label)
        okr = rec > 0
        ax[2, 0].plot(r[okr], qr[okr]/rec[okr]/EV, color=COL[mode], lw=1.0, ls=':',
                      label='recombination photons, ' + label)
    ax[2, 1].plot(r[ok], ((rec + ots)/(s + rec + ots))[ok], color=COL[mode], lw=1.6,
                  label=label)
ax[2, 0].set(xlim=(0, rmax), ylim=(0, 6), xlabel=r'$r$ [pc]',
             ylabel=r'heat per H photoionization [eV]')
ax[2, 0].legend(loc='upper left', fontsize=7)
ax[2, 1].set(xlim=(0, rmax), ylim=(0, 1), xlabel=r'$r$ [pc]',
             ylabel='H photoionizations by\nrecombination photons')
ax[2, 1].legend(loc='upper left', fontsize=7)
h = [Line2D([], [], color='k', lw=1.4), Line2D([], [], color='k', lw=1.4, ls='--',
                                               alpha=0.6)]
fig.legend(h, [r'Tigris ($N=64^3$), no OTS', 'Cloudy, no OTS'], loc='upper center',
           ncol=2, bbox_to_anchor=(0.5, 1.0))
fig.subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.94, wspace=0.25,
                    hspace=0.38)
for ext in figure_formats():
    fig.savefig(OUT + '.' + ext, dpi=250)
