"""Build radial_profiles.txt and H+-zone summary from the Cloudy save files."""
import numpy as np
from collections import defaultdict

PC = 3.0856776e18
R0 = 10**16.5

def load(fn):
    return np.loadtxt(fn, comments="#")

ovr = load("hii.ovr")
rad = load("hii.rad")
H, He, N, O, S = (load("hii.ele_" + e) for e in ("H", "He", "N", "O", "S"))
depth = ovr[:, 0]
assert np.allclose(rad[:, 2], depth, rtol=1e-4)
r = rad[:, 1]
dr = rad[:, 3]
T, ne = ovr[:, 1], ovr[:, 4]

cols = [r / PC, T, ne, H[:, 2], He[:, 2], He[:, 3],
        O[:, 1], O[:, 2], O[:, 3], O[:, 4],
        S[:, 1], S[:, 2], S[:, 3], S[:, 4], N[:, 2], N[:, 3]]
names = ["r_pc", "T_K", "n_e_cm-3", "x_H+", "x_He+", "x_He++",
         "x_O0", "x_O+", "x_O++", "x_O3+", "x_S0", "x_S+", "x_S++", "x_S3+",
         "x_N+", "x_N++"]
tab = np.column_stack(cols)
hdr = ("Cloudy 25 static HII region: Q(H)=1e49 s^-1, SB99 Z=0.014 2 Myr SED, n_H=100 cm^-3, "
       "R_in=10^16.5 cm\nr is the zone-centre distance from the source; x_X is the fraction "
       "of element X in that stage\n" + "  ".join(names))
np.savetxt("radial_profiles.txt", tab, fmt="%.5e", header=hdr)

xH = H[:, 2]
i = np.where(xH < 0.5)[0][0]
Rs = np.interp(0.5, [xH[i], xH[i - 1]], [r[i], r[i - 1]])
print("zones %d, outer edge r = %.3f pc" % (len(r), (r[-1] + dr[-1] / 2) / PC))
print("Stromgren radius (x_H+=0.5): %.4f pc = %.4e cm" % (Rs / PC, Rs))
print("zone thickness at front: dr = %.3e pc (%.2e of R_s)" % (dr[i] / PC, dr[i] / Rs))
print("width 0.9>x_H+>0.1: %.4f pc" % ((np.interp(0.1, xH[::-1], r[::-1]) - np.interp(0.9, xH[::-1], r[::-1])) / PC))

dV = 4 * np.pi * r**2 * dr
inz = r < Rs
wv = dV * inz
wi = dV * xH
lab = ["O+", "O++", "O3+", "S+", "S++", "S3+", "He+", "He++", "N+", "N++"]
dat = [O[:, 2], O[:, 3], O[:, 4], S[:, 2], S[:, 3], S[:, 4], He[:, 2], He[:, 3], N[:, 2], N[:, 3]]
print("\n%-6s %12s %12s" % ("ion", "<x>_V(r<Rs)", "<x>_{H+ mass}"))
for l, d in zip(lab, dat):
    print("%-6s %12.4f %12.4f" % (l, np.sum(d * wv) / wv.sum(), np.sum(d * wi) / wi.sum()))
print("T: volume avg (r<Rs) %.0f K; H+-weighted %.0f K; n_e n_H+ weighted %.0f K" % (
    np.sum(T * wv) / wv.sum(), np.sum(T * wi) / wi.sum(),
    np.sum(T * ne * xH * dV) / np.sum(ne * xH * dV)))

# top coolants inside r<Rs, from 'save cooling' (per-zone agent list)
agents = defaultdict(float)
ctot_all = 0.0
with open("hii.cool") as f:
    for k, line in enumerate(l for l in f if not l.startswith("#")):
        p = line.rstrip("\n ").split("\t")
        if not inz[k]:
            continue
        ctot = float(p[3]) * dV[k]
        ctot_all += ctot
        for j in range(5, len(p) - 1, 2):
            agents[" ".join(p[j].split()[:-1])] += float(p[j + 1]) * ctot
print("\nTop coolants in r<Rs (share of total volume-integrated cooling):")
for a, c in sorted(agents.items(), key=lambda x: -x[1])[:8]:
    print("  %-12s %.3f" % (a, c / ctot_all))
print("  (sum of listed agents covers %.3f of total)" % (sum(agents.values()) / ctot_all))

ce = np.loadtxt("hii.cool_each", comments="#")
with open("hii.cool_each") as f:
    hn = f.readline().lstrip("#").rstrip("\n").split("\t")
tot = np.sum(ce[:, 2] * wv)
sh = sorted(((np.sum(ce[:, j] * wv) / tot, hn[j].strip()) for j in range(3, ce.shape[1])), reverse=True)
print("\nBy element/process (save cooling each), r<Rs:")
for s, n in sh[:7]:
    print("  %-8s %.3f" % (n, s))
