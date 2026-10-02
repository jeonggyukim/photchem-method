"""Checks of the lambda_<ion>.txt tables.

1. Low-density limit and the n_e at which Lambda falls to half of it.
2. Cooling per ion vs Cloudy 25 (constant T = 8000 K HII region).
3. Lambda at n_e = 1 cm^-3 vs pyathena's 1D CHIANTI tables and 5-level
   coolants.
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np

from lambda_table import LambdaTable
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths  # noqa: E402

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
CLOUDY = os.path.join(HERE, 'cloudy_rerun')
IONS = ['o_2', 'o_3', 's_2', 's_3', 's_4', 'n_2', 'n_3', 'ne_2', 'ne_3']


def half_density(tab, T):
    """n_e where Lambda = 0.5 * Lambda(n_e = 1e-2), log-interpolated."""
    y = np.array([tab.log_lambda(T, 10.0**x) for x in tab.logne])
    d = y - y[0] + np.log10(2.0)
    k = np.where(d < 0)[0]
    if k.size == 0:
        return np.inf
    k = k[0]
    return 10.0 ** np.interp(0.0, [d[k], d[k - 1]],
                             [tab.logne[k], tab.logne[k - 1]])


def check_density():
    print("1. low-density limit and half-drop density")
    print(f"{'ion':5s} {'T[K]':>7s} {'Lam(1e-2)':>10s} {'L(1)/L(1e-2)':>12s} "
          f"{'L(1e2)/L(1e-2)':>14s} {'n_e,1/2 [cm^-3]':>16s}")
    for ion in IONS:
        tab = LambdaTable(ion)
        for T in (8.0e3, 1.0e4, 1.0e5):
            L0 = tab(T, 1e-2)
            print(f"{ion:5s} {T:7.0f} {L0:10.3e} {tab(T, 1.0)/L0:12.4f} "
                  f"{tab(T, 1e2)/L0:14.4f} {half_density(tab, T):16.3e}")


def load_cloudy():
    ovr = np.loadtxt(os.path.join(CLOUDY, 'hii.ovr'))
    rad = np.loadtxt(os.path.join(CLOUDY, 'hii.rad'))
    ele = {e: np.loadtxt(os.path.join(CLOUDY, f'hii.ele_{e}'))
           for e in ('H', 'O', 'S', 'N')}
    cool = []
    with open(os.path.join(CLOUDY, 'hii.cool')) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n ').split('\t')
            ctot = float(p[3])
            agents = {p[i].rsplit(' ', 1)[0].strip(): float(p[i + 1])
                      for i in range(5, len(p) - 1, 2)}
            cool.append((ctot, agents))
    return ovr, rad, ele, cool


def check_cloudy():
    print("\n2. Tigris-table cooling / Cloudy cooling, T = 8000 K, "
          "n_H = 100 cm^-3")
    ovr, rad, ele, cool = load_cloudy()
    abund = {'O': 3.2e-4, 'S': 1.45e-5, 'N': 7.4e-5}
    ions = [('O', 2, 'o_2'), ('O', 3, 'o_3'), ('S', 2, 's_2'),
            ('S', 3, 's_3'), ('S', 4, 's_4'), ('N', 2, 'n_2'),
            ('N', 3, 'n_3')]
    tabs = {name: LambdaTable(name) for _, _, name in ions}
    r = rad[:, 1]
    xH = ele['H'][:, 2]
    rs = r[np.where(xH < 0.5)[0][0]]
    zones = [np.argmin(abs(r - f * rs)) for f in (0.3, 0.5, 0.7, 0.9, 0.98)]
    hdr = f"{'ion':6s}" + "".join(f"  z{z:<4d}" for z in zones)
    print("zones (index: r/R_s, n_e, x_H+):")
    for z in zones:
        print(f"  z{z}: r/R_s={r[z]/rs:.3f} n_e={ovr[z, 4]:.2f} "
              f"x_H+={xH[z]:.4f} Ctot={cool[z][0]:.3e}")
    rows = {}
    for el, stage, name in ions:
        label = f"{el:2s}{stage:2d}"
        ratios, ours, theirs, fq = [], [], [], []
        for z in zones:
            ne, nH = ovr[z, 4], ovr[z, 3]
            xq = ele[el][z, stage]
            c_ours = ne * xq * abund[el] * nH * tabs[name](8000.0, ne)
            frac = cool[z][1].get(label, 0.0)
            c_cl = frac * cool[z][0]
            ours.append(c_ours)
            theirs.append(c_cl)
            fq.append(frac)
            ratios.append(c_ours / c_cl if c_cl > 0 else np.nan)
        rows[name] = (ratios, ours, theirs, fq)
    print("\nratio Tigris-table / Cloudy")
    print(hdr)
    for name, (ratios, *_ ) in rows.items():
        print(f"{name:6s}" + "".join(f"  {v:6.3f}" for v in ratios))
    print("\nCloudy fraction of total cooling by that ion")
    print(hdr)
    for name, (_, _, _, fq) in rows.items():
        print(f"{name:6s}" + "".join(f"  {v:6.4f}" for v in fq))
    print("\nTigris-table cooling [erg cm^-3 s^-1]")
    print(hdr)
    for name, (_, ours, _, _) in rows.items():
        print(f"{name:6s}" + "".join(f"  {v:8.2e}" for v in ours))


def read_pyathena_bb(element):
    path = paths.pyathena('data', 'microphysics', 'chianti_v11', f'cool_BB_{element}.txt')
    d = np.loadtxt(path, comments='#')
    return d[:, 0], d[:, 1:]


def check_pyathena():
    print("\n3. Lambda at n_e = 1 cm^-3 vs pyathena [erg cm^3 s^-1]")
    import importlib
    el_name = {'o': 'O', 's': 'S', 'n': 'N', 'ne': 'Ne'}
    print(f"{'ion':5s} {'T[K]':>7s} {'this':>10s} {'pa_1D':>10s} "
          f"{'pa_5lev':>10s} {'this/1D':>8s} {'this/5lev':>9s}")
    for ion in IONS:
        tab = LambdaTable(ion)
        el, st = ion.split('_')
        logT, bb = read_pyathena_bb(el_name[el])
        try:
            mod = importlib.import_module(
                f'pyathena.chemistry.coolants.{ion}')
        except Exception:
            mod = None
        for T in (5.0e3, 8.0e3, 1.0e4, 2.0e4):
            this = tab(T, 1.0)
            one_d = 10 ** np.interp(np.log10(T), logT,
                                    np.log10(bb[:, int(st) - 1] + 1e-300))
            five = np.nan
            if mod is not None:
                Ta = np.array([T])
                ne = np.array([1.0])
                five = float(np.asarray(
                    mod.cooling(Ta, ne, np.array([1.0]), n_p=0.85 * ne))[0])
            print(f"{ion:5s} {T:7.0f} {this:10.3e} {one_d:10.3e} "
                  f"{five:10.3e} {this/one_d:8.3f} {this/five:9.3f}")


if __name__ == '__main__':
    check_density()
    check_cloudy()
    check_pyathena()
