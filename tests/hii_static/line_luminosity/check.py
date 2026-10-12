"""Line luminosities of a converged hii_static run against Cloudy: H-alpha from the case B
fit of Tigris emission.hpp (Draine 2011 eq. 14.8), and the [O III] and [S III] forbidden
lines from the 5-level populations of pyathena's IonCoolant (CHIANTI 11), two ways: an
explicit population solve in every cell, and an eps(T, n_e) table built by the same solver
and interpolated bilinearly in (log T, log n_e) as the C++ reader will.

    4 pi j = n_ion f_u A_ul h nu_ul = n_e n_ion eps(T, n_e)

The explicit solve takes n_p = n_H+ of the cell; the table takes n_p = 0.85 n_e, as the
ion-cooling tables do.

usage: python check.py [RUNDIR [CLOUDY_OUT]]
needs pyathena (the Dropbox checkout) on PYTHONPATH.
"""
import os
import re
import sys
import numpy as np
from pyathena.chemistry.coolants.base import IonCoolant
from photchem_method import paths

athena_read = paths.athena_read()
# ism units: the code density is m_H cm^-3, so n_H = rho/MU_H
MH, KB, MU_H, X_HE = 1.6735575e-24, 1.380649e-16, 1.4, 0.1
NAMES = ['He+', 'OHx', 'CHx', 'CO', 'C+', 'HCO+', 'H2', 'H+', 'H3+', 'H2+', 'O+', 'Si+',
         'C_high', 'Si_high', 'O++', 'O_high', 'S+', 'S++', 'S_high', 'N+', 'N++',
         'N_high', 'He++']
CHARGE = {'He+': 1, 'C+': 1, 'HCO+': 1, 'H+': 1, 'H3+': 1, 'H2+': 1, 'O+': 1, 'Si+': 1,
          'C_high': 2, 'Si_high': 2, 'O++': 2, 'O_high': 3, 'S+': 1, 'S++': 2,
          'S_high': 3, 'N+': 1, 'N++': 2, 'N_high': 3, 'He++': 2}
PC = 3.0857e18
HC = 6.62607015e-27*2.99792458e10
E_HALPHA = HC/6562.80e-8
# (label, ion file, abundance key, upper, lower, Cloudy label, Cloudy wavelength)
LINES = [('[O III] 5007', 'o_3.txt', 'O++', 3, 2, 'O  3', '5006.84A'),
         ('[O III] 4959', 'o_3.txt', 'O++', 3, 1, 'O  3', '4958.91A'),
         ('[O III] 4363', 'o_3.txt', 'O++', 4, 3, 'O  3', '4363.21A'),
         ('[O III] 88um', 'o_3.txt', 'O++', 1, 0, 'O  3', '88.3323m'),
         ('[O III] 52um', 'o_3.txt', 'O++', 2, 1, 'O  3', '51.8004m'),
         ('[S III] 9531', 's_3.txt', 'S++', 3, 2, 'S  3', '9530.62A'),
         ('[S III] 9069', 's_3.txt', 'S++', 3, 1, 'S  3', '9068.62A'),
         ('[S III] 6312', 's_3.txt', 'S++', 4, 3, 'S  3', '6312.06A'),
         ('[S III] 18.7um', 's_3.txt', 'S++', 2, 1, 'S  3', '18.7078m'),
         ('[S III] 33.5um', 's_3.txt', 'S++', 1, 0, 'S  3', '33.4704m')]
LOGT = np.arange(3.0, 5.0 + 1e-9, 0.02)
LOGNE = np.arange(-2.0, 7.0 + 1e-9, 0.1)
NP_OVER_NE = 0.85


def alpha_eff_halpha(T):
    t4 = T*1e-4
    return 1.17e-13*t4**(-0.942 - 0.031*np.log(t4))


def cloudy_lum(fn, label, wl):
    """log10 L [erg/s] of one line from the emission-line block of a Cloudy .out."""
    m = re.search(re.escape(label) + r'\s+' + re.escape(wl) + r'\s+([-\d.]+)',
                  open(fn).read())
    return float(m.group(1)) if m else np.nan


def cells(rundir):
    fn = sorted(f for f in os.listdir(rundir)
                if re.match(r'.*\.out2\.\d+\.athdf$', f))[-1]
    d = athena_read.athdf(os.path.join(rundir, fn))
    keys = ['rHI'] + ['r%d' % n for n in range(1, len(NAMES))]
    x = {nm: d[k].astype(np.float64).ravel() for nm, k in zip(NAMES, keys)}
    nH = d['rho'].astype(np.float64).ravel()/MU_H
    xe = sum(CHARGE.get(nm, 0)*x[nm] for nm in NAMES)
    T = (d['press'].astype(np.float64).ravel()*MH*1e10/KB
         / (nH*(1.0 - x['H2'] + X_HE + xe)))
    dV = np.prod([float(d['x%df' % i][1] - d['x%df' % i][0]) for i in (1, 2, 3)])*PC**3
    return fn, nH, T, xe*nH, x, dV


def line_power(tab, f, u, lo):
    """f_u A_ul h nu_ul [erg/s per ion]."""
    return f[u]*tab['A'][u, lo]*(tab['E_erg'][u] - tab['E_erg'][lo])


def build_table(ion, u, lo):
    """log10 eps(T, n_e) [erg cm^3/s] on (LOGT, LOGNE)."""
    tt, nn = np.meshgrid(10**LOGT, 10**LOGNE, indexing='ij')
    f = ion.populations(tt.ravel(), nn.ravel(), NP_OVER_NE*nn.ravel())
    eps = line_power(ion.table(), f, u, lo)/nn.ravel()
    return np.log10(eps).reshape(tt.shape)


def interp(logeps, T, ne):
    """Bilinear in (log T, log n_e), clamped to the grid edges."""
    def index(v, grid):
        s = (np.clip(v, grid[0], grid[-1]) - grid[0])/(grid[1] - grid[0])
        i = np.minimum(s.astype(int), len(grid) - 2)
        return i, s - i
    i, a = index(np.log10(T), LOGT)
    j, b = index(np.log10(ne), LOGNE)
    v = ((1 - a)*(1 - b)*logeps[i, j] + a*(1 - b)*logeps[i + 1, j]
         + (1 - a)*b*logeps[i, j + 1] + a*b*logeps[i + 1, j + 1])
    return 10**v


def main():
    rundir = (sys.argv[1] if len(sys.argv) > 1 else
              paths.runs('T6_multi_ion', 'hii_static_bands10', 'sb99_ots'))
    cloudy = (sys.argv[2] if len(sys.argv) > 2 else
              paths.runs('T6_multi_ion', 'cloudy_sb99', 'ots', 'hii.out'))
    fn, nH, T, ne, x, dV = cells(rundir)
    m = (ne > 0) & (T > 0) & (x['H+'] > 1e-3)
    nH, T, ne = nH[m], T[m], ne[m]
    x = {k: v[m] for k, v in x.items()}
    print('%s/%s: %d cells with x_H+ > 1e-3, T %.0f-%.0f K, n_e %.3g-%.3g cm^-3'
          % (rundir, fn, m.sum(), T.min(), T.max(), ne.min(), ne.max()))
    print('Cloudy: %s\n' % cloudy)
    rows = []
    L_ha = np.sum(E_HALPHA*alpha_eff_halpha(T)*ne*x['H+']*nH)*dV
    rows.append(('H-alpha (case B)', np.log10(L_ha), np.nan,
                 cloudy_lum(cloudy, 'H  1', '6562.80A'), np.nan, np.nan, np.nan, np.nan))
    lum = {'H-alpha (case B)': L_ha}
    ions = {}
    for lab, ifile, key, u, lo, clab, cwl in LINES:
        if ifile not in ions:
            ion = IonCoolant(ifile, ifile[:-4])
            ions[ifile] = (ion, ion.populations(T, ne, x['H+']*nH),
                           ion.populations(T, ne, NP_OVER_NE*ne))
        ion, f, f_np = ions[ifile]
        nion = x[key]*nH
        j_exp = nion*line_power(ion.table(), f, u, lo)
        j_np = nion*line_power(ion.table(), f_np, u, lo)
        j_tab = ne*nion*interp(build_table(ion, u, lo), T, ne)
        w = j_exp > 1e-6*j_exp.max()
        on_grid = w & (T >= 10**LOGT[0])
        err_interp = np.abs(j_tab[on_grid]/j_np[on_grid] - 1.0).max()
        err_np = np.abs(j_np[w]/j_exp[w] - 1.0).max()
        off_grid = j_exp[T < 10**LOGT[0]].sum()/j_exp.sum()
        L_exp, L_tab = j_exp.sum()*dV, j_tab.sum()*dV
        lum[lab] = L_exp
        rows.append((lab, np.log10(L_exp), np.log10(L_tab), cloudy_lum(cloudy, clab, cwl),
                     L_tab/L_exp - 1.0, err_interp, err_np, off_grid))
    tmin = 10**LOGT[0]
    print('interp max: largest cell error of the table against the solve at the same')
    print('  n_p, over cells with T >= %.0f K.' % tmin)
    print('n_p max: largest cell change from n_p = n_H+ to n_p = %.2f n_e.' % NP_OVER_NE)
    print('L(T<1e3K): fraction of the line from cells below the grid, where the table')
    print('  holds its edge value.\n')
    print('%-17s %9s %9s %9s %11s %12s %11s %9s %11s'
          % ('line', 'logL exp', 'logL tab', 'Cloudy', 'exp/Cloudy', 'tab/exp - 1',
             'interp max', 'n_p max', 'L(T<1e3K)'))
    for lab, le, lt, lc, dt, ei, en, og in rows:
        print('%-17s %9.3f %9.3f %9.3f %11.3f %12.2e %11.2e %9.1e %11.1e'
              % (lab, le, lt, lc, 10**(le - lc), dt, ei, en, og))

    def cl(lab):
        r = [r for r in rows if r[0] == lab][0]
        return 10**r[3]
    print('\nratio              Tigris     Cloudy')
    for a, b in [('[O III] 4363', '[O III] 5007'), ('[S III] 6312', '[S III] 9069'),
                 ('[S III] 9531', '[S III] 9069'), ('[O III] 5007', 'H-alpha (case B)'),
                 ('[S III] 9069', 'H-alpha (case B)')]:
        print('%-12s/%-5s %9.4f %10.4f' % (a, b.split()[-1], lum[a]/lum[b], cl(a)/cl(b)))


if __name__ == '__main__':
    main()
