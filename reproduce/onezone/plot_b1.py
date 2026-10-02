"""Test B1 figure: Tigris GOW17 one-zone end state vs Athena++ (reference values in
AthenaK test_GOW17_uniform_gpu.py) and AthenaK kokkos_BDF (last record of
gow17_time_series.npz). Input: b1_uniform.txt from b1_uniform.cpp."""
import sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402
names = ["He+","OHx","CHx","CO","C+","HCO+","H2","H+","H3+","H2+","O+","Si+"]
ref = dict(eint=28.234391212463375, **dict(zip(names, [4.4769834062208247e-07,
    1.4335562809719704e-05, 6.079905467970548e-09, 0.00013345545448828489,
    3.4978197049895243e-07, 6.033158683749207e-08, 0.44989413022994995,
    2.9430591439449927e-06, 1.2544340961540001e-06, 1.969339358254274e-09,
    8.677133317425145e-11, 6.625605806220847e-07])))
ts = np.load(paths.athenak("tst/test_suite/chemistry/data/gow17_time_series.npz"))["time_series"][-1]
bdf = dict(eint=float(ts["eint"]), **{n: float(ts[f"s_{i:02d}_chem_{n}"]) for i, n in enumerate(names)})
row = np.atleast_2d(np.loadtxt("b1_uniform.txt"))[0]
keys = ["eint"] + names
tig = dict(zip(keys, [row[2]] + list(row[3:15])))
fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
xs = np.arange(len(keys))
for off, (lab, d, c) in enumerate([("AthenaK kokkos_BDF (last record of its time series)", bdf, "0.5"),
                                   ("Tigris semi-implicit (Gauss-Seidel)", tig, "C3")]):
    rel = np.array([d[k] for k in keys])/np.array([ref[k] for k in keys]) - 1
    ax.bar(xs + (off - 0.5)*0.38, 100*rel, width=0.38, color=c, label=lab)
ax.axhline(0, color="k", lw=0.8)
for y in (-10, -1, 1, 10):
    ax.axhline(y, color="k", lw=0.5, ls=":")
ax.set_yscale("symlog", linthresh=1)
ax.set_xticks(xs); ax.set_xticklabels(["e_int"] + names)
ax.set_ylabel("difference from Athena++ [%]")
ax.set_title("GOW17 one zone, n$_H$ = 109.21 cm$^{-3}$, t = 9.8 Gyr: "
             "end state vs Athena++ (AthenaK test reference)")
ax.legend(loc="upper left", frameon=False)
fig.savefig("fig_b1_uniform_endstate.png", dpi=150)
for k in keys:
    print(f"{k:5s} BDF {100*(bdf[k]/ref[k]-1):+7.2f}%  Tigris {100*(tig[k]/ref[k]-1):+7.2f}%")
