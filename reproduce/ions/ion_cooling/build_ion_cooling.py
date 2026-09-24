"""Build per-ion collisional line-cooling tables Lambda_q(T, n_e) from CHIANTI.

Line cooling rate per unit volume by ion q = n_e * n_q * Lambda_q(T, n_e).

Level populations come from ChiantiPy's full-level solver (all levels in the
CHIANTI .elvlc file) with electron collisions (.scups) and, where CHIANTI has
them (.psplups), proton collisions at n_p = P2E * n_e (with the rate
correction in fix_proton_rates).  Recombination and ionization into/out of
the level system are switched off (CHIANTI 11 has no .rrlvl files for these
ions, so ChiantiPy's default already leaves them out), so Lambda_q is
bound-bound emission driven by collisional excitation only.

Lambda_q = sum_lines n_u A_ul h nu_ul / n_e, with n_u the fractional population
of the upper level (sum over the ion's levels = 1), all lines in .wgfa
including those with only theoretical wavelengths.

Run:
    export XUVTOP=$HOME/Dropbox/Projects/CHIANTI_db
    python build_ion_cooling.py
"""
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings('ignore')

IONS = ['o_2', 'o_3', 's_2', 's_3', 's_4', 'n_2', 'n_3', 'ne_2', 'ne_3']
LOGT = np.round(np.arange(3.0, 9.0 + 1e-9, 0.05), 4)
LOGNE = np.round(np.arange(-2.0, 6.0 + 1e-9, 0.5), 4)
# n_p/n_e for fully ionized H with He singly to doubly ionized (0.91 to 0.83).
P2E = 0.85
LOG_FLOOR = -300.0

PLANCK = 6.62607015e-27
CLIGHT = 2.99792458e10


def chianti_version():
    with open(os.path.join(os.environ['XUVTOP'], 'VERSION')) as f:
        return f.read().strip()


def spectroscopic(name):
    el, st = name.split('_')
    roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']
    return el.capitalize() + ' ' + roman[int(st) - 1]


def line_loss(ion):
    """Lambda per ion per electron [erg cm^3 s^-1] from ion.Population."""
    pop = np.atleast_2d(ion.Population['population'])
    pop = pop / pop.sum(axis=1, keepdims=True)
    wvl = np.abs(np.asarray(ion.Wgfa['wvl'], float))
    ok = wvl > 0.0
    l2 = np.asarray(ion.Wgfa['lvl2'], int)[ok] - 1
    a = np.asarray(ion.Wgfa['avalue'], float)[ok]
    hnu = PLANCK * CLIGHT / (wvl[ok] * 1e-8)
    power = pop[:, l2] @ (a * hnu)
    return power / np.asarray(ion.EDensity, float)


def fix_proton_rates(ion):
    """CHIANTI .psplups tabulates the proton excitation rate coefficient
    [cm^3 s^-1] directly.  ChiantiPy 0.16 Ion.upsilonDescale(prot=1) treats
    that number as an effective collision strength and multiplies it by
    8.63e-6 / (g sqrt(T)), which makes proton rates ~1e8 times too small.
    Replace exRate/dexRate with the tabulated rate and detailed balance."""
    orig = ion.upsilonDescale

    def patched(prot=False):
        orig(prot=prot)
        if not prot or ion.PUpsilon is None:
            return
        rate = np.maximum(np.asarray(ion.PUpsilon['upsilon'], float), 0.0)
        T = np.asarray(ion.Temperature, float)
        lvl = list(ion.Elvlc['lvl'])
        ecm = np.where(np.asarray(ion.Elvlc['ecm']) >= 0,
                       ion.Elvlc['ecm'], ion.Elvlc['ecmth'])
        mult = np.asarray(ion.Elvlc['mult'], float)
        i1 = np.array([lvl.index(l) for l in ion.Psplups['lvl1']])
        i2 = np.array([lvl.index(l) for l in ion.Psplups['lvl2']])
        dE_K = 1.4387769 * np.abs(ecm[i2] - ecm[i1])
        x = dE_K[:, None] / T[None, :]
        # For dE >> kT the spline fit is outside its useful range and the
        # proton rate is negligible anyway.
        rate = np.where(x < 50.0, rate, 0.0)
        ion.PUpsilon['exRate'] = rate
        ion.PUpsilon['dexRate'] = (rate * (mult[i1] / mult[i2])[:, None]
                                   * np.exp(np.minimum(x, 50.0)))

    ion.upsilonDescale = patched


def compute(name, logT, logne, p2e=P2E, recomb=False):
    """Lambda_q on the flattened (T, n_e) pairs; returns array (len(logT),)."""
    import ChiantiPy.core as ch
    T = 10.0 ** np.asarray(logT, float)
    ne = 10.0 ** np.asarray(logne, float)
    ion = ch.ion(name, temperature=T, eDensity=ne)
    ion.PDensity = p2e * ne
    fix_proton_rates(ion)
    if not recomb:
        ion.Nrrlvl = 0
        ion.Nauto = 0
    ion.populate(popCorrect=recomb)
    return line_loss(ion), ion


def write_table(path, name, version, lam2d, nlvls, elapsed):
    with np.errstate(divide='ignore'):
        loglam = np.where(lam2d > 0, np.log10(np.maximum(lam2d, 1e-300)),
                          LOG_FLOOR)
    with open(path, 'w') as f:
        f.write(f"# ion: {spectroscopic(name)} (CHIANTI name {name})\n")
        f.write(f"# source: CHIANTI {version} via ChiantiPy 0.16.0 "
                f"ion.populate, {nlvls} levels\n")
        f.write("# quantity: log10 Lambda_q, Lambda_q in erg cm^3 s^-1; "
                "line cooling per volume = n_e n_q Lambda_q\n")
        f.write("# process: collisional excitation -> bound-bound emission "
                "only; e- collisions + p collisions (n_p = "
                f"{P2E} n_e, where CHIANTI has .psplups); "
                "recombination/ionization into levels off\n")
        f.write(f"# rows: log10 T [K], {len(LOGT)} points, "
                f"{LOGT[0]:.2f} to {LOGT[-1]:.2f} step 0.05\n")
        f.write(f"# cols: log10 n_e [cm^-3], {len(LOGNE)} points, "
                f"{LOGNE[0]:.2f} to {LOGNE[-1]:.2f} step 0.5\n")
        f.write(f"# floor: {LOG_FLOOR:g} where Lambda underflows\n")
        f.write(f"# build time: {elapsed:.1f} s\n")
        f.write("# logT_grid: " + " ".join(f"{x:.2f}" for x in LOGT) + "\n")
        f.write("# logne_grid: " + " ".join(f"{x:.2f}" for x in LOGNE)
                + "\n")
        for i in range(len(LOGT)):
            f.write(" ".join(f"{v:10.5f}" for v in loglam[i]) + "\n")


def main():
    if not os.environ.get('XUVTOP'):
        sys.exit("set XUVTOP to the CHIANTI database directory")
    out = os.path.dirname(os.path.abspath(__file__))
    version = chianti_version()
    ions = sys.argv[1:] or IONS
    TT, NN = np.meshgrid(LOGT, LOGNE, indexing='ij')
    t_all = time.time()
    for name in ions:
        t0 = time.time()
        lam, ion = compute(name, TT.ravel(), NN.ravel())
        lam2d = lam.reshape(TT.shape)
        dt = time.time() - t0
        write_table(os.path.join(out, f"lambda_{name}.txt"), name, version,
                    lam2d, ion.Nlvls, dt)
        print(f"{name:5s} {ion.Nlvls:4d} levels  {dt:7.1f} s  "
              f"Lambda(1e4 K, 1 cm^-3) = "
              f"{lam2d[np.argmin(abs(LOGT-4)), np.argmin(abs(LOGNE))]:.3e}",
              flush=True)
    print(f"total {time.time() - t_all:.1f} s for {len(ions)} ions, "
          f"{TT.size} (T, n_e) points each")


if __name__ == '__main__':
    main()
