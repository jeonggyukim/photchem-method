#!/usr/bin/env python
"""Figure 1: the GOW17 network and its Tigris extension, in four versions.

    core        GOW17 core only (the AthenaK diagram of plot_gow17_network.py)
    greycore    full Tigris network, additions in colour, GOW17 core in grey
    2panel      (a) core + ionized-gas additions, (b) ion ladders
    colour      full Tigris network in full colour, additions highlighted

and a contact sheet of the four. The core reactions and their layout are copied
from athenak-chem/scripts/plot_gow17_network.py. Every Tigris arrow cites the
code in tigris-rayt/src/photchem/network/gow17_network.hpp (branch
rayt-photchem-updates, configure -gow17 --photchem_ions=O3,S3,N3: O, S and N
tracked to O III, S III, N III, plus one higher-ion species per element for
C, Si, O, S and N; gow17_core.hpp:25-31).

Usage:
    fig_network.py [--outdir ../figures]
"""
import argparse
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                     "pdf.fonttype": 42})

WIDTH_IN = 7.0
FS_SPECIES, FS_LABEL, FS_LEGEND, FS_NOTE = 8.5, 7.0, 7.0, 7.0

# Okabe & Ito (2008) colours.
VERMILLION, PURPLE, BLUE, ORANGE = "#D55E00", "#CC79A7", "#0072B2", "#E69F00"
GREEN, SKY, BLACK, DARK = "#009E73", "#56B4E9", "#1A1A1A", "#555555"
GREY, GREY_TEXT, GREY_FILL = "#BDBDBD", "#9E9E9E", "#F4F4F4"
HIGHLIGHT = "#F7E3A1"
PARTNER = "#6B6B6B"

STYLE = {
    "cr": dict(color=VERMILLION, ls="-", lw=1.1),
    "crphoto": dict(color=VERMILLION, ls=(0, (3.5, 1.8)), lw=1.1),
    "uv": dict(color=PURPLE, ls=(0, (3.5, 1.8)), lw=1.1),
    "euv": dict(color=PURPLE, ls="-", lw=1.4),
    "grain": dict(color=BLUE, ls="-", lw=1.1),
    "rec": dict(color=DARK, ls="-", lw=0.7),
    "chem": dict(color=BLACK, ls="-", lw=0.9),
    "coll": dict(color=ORANGE, ls=(0, (4, 1.5, 1, 1.5)), lw=1.1),
    "ct": dict(color=GREEN, ls="-", lw=1.1),
    "ots": dict(color=SKY, ls=(0, (1.2, 1.2)), lw=1.5),
}
LEGEND = {
    "cr": "cosmic-ray ionization",
    "crphoto": "CR-induced UV photoreaction",
    "uv": "FUV photoreaction",
    "euv": "ionizing-band photoionization (eV)",
    "grain": "grain-assisted",
    "rec": "recombination with e",
    "chem": "ion-molecule / neutral-neutral",
    "coll": "collisional (e; H, H$_2$)",
    "ct": "charge transfer (H$^+$ up, H down)",
    "ots": "He recombination photons on H",
}

TEX = {"He+": "He$^+$", "H2": "H$_2$", "H+": "H$^+$", "H2+": "H$_2^+$",
       "H3+": "H$_3^+$", "C+": "C$^+$", "CHx": "CH$_x$", "HCO+": "HCO$^+$",
       "OHx": "OH$_x$", "O+": "O$^+$", "Si+": "Si$^+$", "O2+": "O$^{2+}$",
       "S+": "S$^+$", "S2+": "S$^{2+}$", "N": "N", "N+": "N$^+$", "N2+": "N$^{2+}$",
       "C_high": r"C$^{\geq 2+}$", "Si_high": r"Si$^{\geq 2+}$",
       "O_high": r"O$^{\geq 3+}$", "S_high": r"S$^{\geq 3+}$",
       "N_high": r"N$^{\geq 3+}$"}
GHOSTS = {"H", "He", "C", "O", "Si", "S", "N"}
MOLECULES = {"H2", "CHx", "OHx", "CO"}
LIMITERS = {"H2", "CO"}
# gow17_core.hpp:25-31 (Species enum under PHOTCHEM_IONS); S and N are ghosts
# (GhostSpecies, gow17_network.hpp:318-319).
NEW_SPECIES = {"S", "S+", "S2+", "S_high", "O2+", "O_high", "N", "N+", "N2+",
               "N_high", "C_high", "Si_high"}
# Elements that a ladder could be written for but that have no species or rates.
OPTIONAL = ["Ne", "Mg", "Fe"]
# Block Gauss-Seidel order (OrderedGaussSeidelUpdate, default placement of H2 last;
# CO and HCO+ are one 2x2 block, numbered in that order).
GS_CORE = ["He+", "Si+", "H2+", "H3+", "H+", "O+", "C+", "CHx", "OHx", "CO",
           "HCO+", "H2"]
# gow17_network.hpp:1834-1849 with the defaults semi_implicit_hep_first = true and
# h2_first = false: step_si (Si+, Si_high), step_sn (S+, S++, S_high, N+, N++,
# N_high), step_h2p, step_h3p, step_hp, step_o (O+, O++, O_high), step_c (C+,
# C_high), step_chx, step_ohx, step_co, step_h2.
GS_TIGRIS = ["He+", "Si+", "Si_high", "S+", "S2+", "S_high", "N+", "N2+", "N_high",
             "H2+", "H3+", "H+", "O+", "O2+", "O_high", "C+", "C_high", "CHx", "OHx",
             "CO", "HCO+", "H2"]

POS_CORE = {
    "He": (7.4, 8.0), "He+": (7.4, 5.5),
    "H": (1.0, 8.0), "H2": (4.2, 8.0), "H+": (1.0, 3.0),
    "H2+": (4.2, 5.5), "H3+": (4.2, 3.0),
    "C": (10.6, 8.0), "C+": (11.0, 5.5), "CO": (17.0, 8.0),
    "CHx": (13.8, 8.0), "HCO+": (17.0, 5.0),
    "O": (13.8, 0.6), "OHx": (13.8, 3.0), "O+": (10.6, 0.6),
    "Si": (-1.8, 8.0), "Si+": (-1.8, 5.5),
}
POS_FULL = dict(POS_CORE, **{
    "S": (-5.0, 8.0), "S+": (-5.0, 5.5), "S2+": (-5.0, 3.0), "S_high": (-5.0, 0.6),
    "O2+": (7.4, 0.6), "O_high": (4.2, 0.6),
    "N": (20.2, 8.0), "N+": (20.2, 5.5), "N2+": (20.2, 3.0), "N_high": (20.2, 0.6),
    "Si_high": (-1.8, 3.0), "C_high": (11.0, 3.0),
})
POS_LADDERS = {
    "S": (0.0, 2.6), "S+": (3.2, 2.6), "S2+": (6.4, 2.6), "S_high": (9.6, 2.6),
    "O": (0.0, 0.0), "O+": (3.2, 0.0), "O2+": (6.4, 0.0), "O_high": (9.6, 0.0),
    "N": (0.0, -2.6), "N+": (3.2, -2.6), "N2+": (6.4, -2.6), "N_high": (9.6, -2.6),
    "C+": (12.9, 2.6), "C_high": (16.1, 2.6),
    "Si+": (12.9, 0.0), "Si_high": (16.1, 0.0),
}
BOX_W, BOX_H = 1.5, 0.8


def E(src, dst, kind, label="", arc=0.0, t=0.5, new=False, emph=False, off=None,
      lab_t=0.5, lab_off=0.42):
    return dict(src=src, dst=dst, kind=kind, label=label, arc=arc, t=t, new=new,
                emph=emph, off=off, lab_t=lab_t, lab_off=lab_off)


# The 49 GOW17 reactions, one arrow each to the main product; copied from
# plot_gow17_network.py EDGES. The two H+ <-> O+ arrows are (2body 27) H+ + O and
# (2body 28) O+ + H; the Tigris layouts redraw them as O <-> O+ charge transfer.
CORE_EDGES = [
    E("H2", "H2+", "cr", "CR", 0.0), E("He", "He+", "cr", "CR", 0.25),
    E("H", "H+", "cr", "CR", 0.25), E("C", "C+", "cr", "CR", 0.25),
    E("CO", "C", "crphoto", "CR", 0.35), E("CO", "HCO+", "cr", "CR", 0.3),
    E("Si", "Si+", "cr", "CR", 0.3),
    E("C", "C+", "uv", "", -0.25), E("CHx", "C", "uv", "", 0.3),
    E("CO", "C", "uv", "", -0.25), E("OHx", "O", "uv", "", 0.3),
    E("H2", "H", "uv", "", 0.35), E("Si", "Si+", "uv", "", -0.3),
    E("H", "H2", "grain", "H", 0.4), E("H+", "H", "grain", "e", -0.55),
    E("C+", "C", "grain", "e", -0.6), E("He+", "He", "grain", "e", 0.3),
    E("Si+", "Si", "grain", "e", 0.55),
    E("He+", "He", "rec", "e", 0.0), E("H3+", "H2", "rec", "e", -0.45),
    E("H3+", "H", "rec", "e", 0.2), E("C+", "C", "rec", "e", 0.0),
    E("HCO+", "CO", "rec", "e", 0.3), E("H+", "H", "rec", "e", 0.0),
    E("Si+", "Si", "rec", "e", 0.0),
    E("H3+", "CHx", "chem", "C", 0.0), E("H3+", "OHx", "chem", "O", -0.25),
    E("CO", "HCO+", "chem", "H$_3^+$", -0.3),
    E("He+", "H+", "chem", "H$_2$", 0.12, 0.25),
    E("CO", "C+", "chem", "He$^+$", -0.3), E("C+", "CHx", "chem", "H$_2$", 0.0),
    E("C+", "HCO+", "chem", "OH$_x$", 0.0, 0.3),
    E("OHx", "HCO+", "chem", "C$^+$", 0.2, 0.72),
    E("CHx", "CO", "chem", "O", 0.25), E("OHx", "CO", "chem", "C", 0.25),
    E("H2+", "H3+", "chem", "H$_2$", 0.0), E("He+", "H2+", "chem", "H$_2$", 0.0),
    E("CHx", "C", "chem", "H", -0.3), E("OHx", "O", "chem", "O", -0.3),
    E("OHx", "O+", "chem", "He$^+$", 0.2),
    E("H2+", "H+", "chem", "H", 0.0),
    E("H+", "O+", "chem", "O", 0.12),        # (2body 27)
    E("O+", "H+", "chem", "H", 0.12),        # (2body 28)
    E("O+", "OHx", "chem", "H$_2$", 0.2), E("O+", "O", "chem", "H$_2$", 0.0),
    E("C+", "C", "chem", "H$_2$, e", 0.5), E("H3+", "H2", "chem", "O, e", 0.45, 0.82),
    E("H2", "H", "coll", "H, H$_2$", 0.0), E("H", "H+", "coll", "e", -0.25),
]
PARTNERS = [("CO", "C+", "He+"), ("OHx", "O+", "He+"), ("CO", "HCO+", "H3+")]


def core_edges_tigris():
    """Core arrows without (2body 27, 28); the O ladder redraws them as O <-> O+."""
    return [e for e in CORE_EDGES
            if not ((e["src"], e["dst"]) in (("H+", "O+"), ("O+", "H+")))]


# Tigris, ionized gas. gow17_network.hpp: SetLyC xi_ph_hi_ (H+ step "a") and
# xi_ph_h2_ (H2+ creation, H2 destruction); SetIonPhotoRates xi_ph_he_
# (HePlusStep_); HeOTS_ s_h in the H+ creation, k_he in the He+ destruction.
ION_GAS = [
    E("H", "H+", "euv", "13.6", -0.62, new=True),
    E("H2", "H2+", "euv", "15.4", -0.3, new=True, lab_t=0.72, lab_off=0.75),
    E("He", "He+", "euv", "24.6", 0.55, new=True),
    E("He+", "H+", "ots", "", -0.18, new=True),
]
# C I and Si I gain collisional ionization (ICI) and charge transfer with H+
# (ICT_ION) in step_c (gow17_network.hpp:1683-1684) and step_si (1437-1438).
NEUTRAL_CI = [
    E("C", "C+", "coll", "", 0.1, new=True), E("C", "C+", "ct", "", -0.1, new=True),
    E("Si", "Si+", "coll", "", 0.12, new=True), E("Si", "Si+", "ct", "", -0.12, new=True),
]


def U(kind, label="", new=True, emph=False):
    return (kind, label, new, emph)


STD_UP = lambda eth: [U("euv", eth), U("coll"), U("ct")]  # noqa: E731
STD_DOWN = [U("ct"), U("rec")]


def ladder_step(lo, hi, ups, downs, sign=1.0):
    """Parallel arrows lo -> hi (ups, outermost first) and hi -> lo (downs).

    Offsets are perpendicular to the pair, in units of the box half-extent,
    0.36 apart (0.3 for six); sign puts ionization on one side and
    recombination on the other.
    """
    arrows = [(lo, hi, u) for u in ups] + [(hi, lo, d) for d in downs]
    n = len(arrows)
    gap = 0.36 if n <= 5 else 0.3
    out = []
    for i, (src, dst, (kind, label, new, emph)) in enumerate(arrows):
        off = (i - 0.5 * (n - 1)) * gap * sign
        out.append(E(src, dst, kind, label, new=new, emph=emph, off=off))
    return out


def ladder_edges(sign):
    """The ion ladders of --photchem_ions=O3,S3,N3, gow17_network.hpp.

    Standard step X^q -> X^{q+1}: ionizing-band photoionization (threshold in eV,
    verner96_photx.dat), collisional ionization ion_k_[ICI], charge transfer with
    H+ ion_k_[ICT_ION]; back by charge transfer with H ion_k_[ICT_REC] and
    radiative recombination ion_k_[IREC]. Into X_high: HighIonization (CI and
    CT_ION of the top stage, :3120-3124) plus the top stage's photoionization;
    back: HighReturn (IREC, ICT_REC of the higher-ion table, :3127-3132).

    O -> O+ (step_o, :1620-1645): GOW17 (2body 27, 28) charge transfer, plus
    ICI of O I and O I by cosmic-ray-induced UV (kcr_o_, :2476); back by IREC of
    O+. No O I photoionization: SetIonPhotoRates has no O I rate.
    O+ -> O++ (xi_ph_oii_, 35.1 eV), O++ -> O_high (xi_ph_oiii_, 54.9 eV; :1647-1662).
    S -> S+ (step_sn, :1460-1471): xi_ph_si_ (10.4 eV, all bands and the ISRF,
    photchem_gow17.cpp:764), kcr_s_ (:2475), ICI, ICT_ION. S+ -> S++ (xi_ph_sii_,
    23.3 eV), S++ -> S_high (xi_ph_siii_, 34.8 eV; :1477-1488).
    N (:1491-1520): xi_ph_ni_ (14.5 eV), xi_ph_nii_ (29.6 eV), xi_ph_niii_
    (47.4 eV, into N_high).
    C+ -> C_high (step_c, :1685-1695; xi_ph_cii_, 24.4 eV) and Si+ -> Si_high
    (step_si, :1439-1450; xi_ph_si_ii_, 16.3 eV).
    """
    o_up = [U("crphoto"), U("coll"), U("ct", new=False, emph=True)]
    o_down = [U("ct", new=False, emph=True), U("rec")]
    s_up = [U("uv", "10.4"), U("crphoto"), U("coll"), U("ct")]
    return (ladder_step("O", "O+", o_up, o_down, sign["O"])
            + ladder_step("O+", "O2+", STD_UP("35.1"), STD_DOWN, sign["O"])
            + ladder_step("O2+", "O_high", STD_UP("54.9"), STD_DOWN, sign["O"])
            + ladder_step("S", "S+", s_up, STD_DOWN, sign["S"])
            + ladder_step("S+", "S2+", STD_UP("23.3"), STD_DOWN, sign["S"])
            + ladder_step("S2+", "S_high", STD_UP("34.8"), STD_DOWN, sign["S"])
            + ladder_step("N", "N+", STD_UP("14.5"), STD_DOWN, sign["N"])
            + ladder_step("N+", "N2+", STD_UP("29.6"), STD_DOWN, sign["N"])
            + ladder_step("N2+", "N_high", STD_UP("47.4"), STD_DOWN, sign["N"])
            + ladder_step("C+", "C_high", STD_UP("24.4"), STD_DOWN, sign["C"])
            + ladder_step("Si+", "Si_high", STD_UP("16.3"), STD_DOWN, sign["Si"]))


# ---------------------------------------------------------------------------
# drawing


def is_new(name, version):
    return name in NEW_SPECIES and version != "core"


# Per-figure box sizes: the O row gets taller boxes for its five parallel arrows.
SIZE_OVERRIDE = {}


def box_size(name):
    if name in SIZE_OVERRIDE:
        return SIZE_OVERRIDE[name]
    return (BOX_W, BOX_H)


def straight_ends(p, q, ws, wd, off):
    """End points of a straight arrow between axis-aligned neighbours, shifted
    perpendicular to the pair by off times the smaller half-extent."""
    if abs(q[0] - p[0]) < 1e-9:
        half = 0.5 * min(ws[0], wd[0]) - 0.08
        x = p[0] + off * half
        sg = 1.0 if q[1] > p[1] else -1.0
        return (x, p[1] + sg * ws[1] / 2), (x, q[1] - sg * wd[1] / 2)
    half = 0.5 * min(ws[1], wd[1]) - 0.06
    y = p[1] + off * half
    sg = 1.0 if q[0] > p[0] else -1.0
    return (p[0] + sg * ws[0] / 2, y), (q[0] - sg * wd[0] / 2, y)


def box_edge(p, q, w=BOX_W, h=BOX_H):
    dx, dy = q[0] - p[0], q[1] - p[1]
    sx = (w / 2) / abs(dx) if dx else float("inf")
    sy = (h / 2) / abs(dy) if dy else float("inf")
    s = min(sx, sy)
    return (p[0] + s * dx, p[1] + s * dy)


def on_arc(a, b, arc, t):
    cx = (a[0] + b[0]) / 2 + arc * (b[1] - a[1])
    cy = (a[1] + b[1]) / 2 - arc * (b[0] - a[0])
    s = 1.0 - t
    return (s * s * a[0] + 2 * s * t * cx + t * t * b[0],
            s * s * a[1] + 2 * s * t * cy + t * t * b[1])


def bulge_normal(a, b, arc):
    dx, dy = b[0] - a[0], b[1] - a[1]
    n = (dx * dx + dy * dy) ** 0.5
    sgn = 1.0 if arc >= 0 else -1.0
    return (sgn * dy / n, -sgn * dx / n)


def draw_arrow(ax, a, b, arc, st, shrink=1.5, under=None, zorder=1.0):
    common = dict(connectionstyle=f"arc3,rad={arc}", shrinkA=shrink, shrinkB=shrink)
    if under:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", lw=st["lw"] + 2.4,
                                     color=under, capstyle="round",
                                     zorder=zorder - 0.3, **common))
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", lw=st["lw"], ls=st["ls"],
                                 color=st["color"], alpha=st.get("alpha", 1.0),
                                 zorder=zorder, **common))
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>,head_length=4.2,head_width=2.1",
                                 lw=0, color=st["color"], alpha=st.get("alpha", 1.0),
                                 zorder=zorder, **common))


def grey_style(st):
    return dict(color=GREY, ls=st["ls"], lw=max(0.5, 0.7 * st["lw"]), alpha=0.9)


def text_halo(ax, x, y, s, color, fs=FS_LABEL, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, fontsize=fs, color=color, zorder=3,
            path_effects=[patheffects.withStroke(linewidth=2.5, foreground="white")],
            **kw)


def draw_edges(ax, edges, pos, mode):
    """mode: 'plain', 'grey' (core grey, new in colour) or 'highlight'."""
    mids = {}
    for e in edges:
        p, q = pos[e["src"]], pos[e["dst"]]
        if e["off"] is not None:
            a, b = straight_ends(p, q, box_size(e["src"]), box_size(e["dst"]), e["off"])
        else:
            a = box_edge(p, q, *box_size(e["src"]))
            b = box_edge(q, p, *box_size(e["dst"]))
        st = dict(STYLE[e["kind"]])
        if e["emph"]:
            st["lw"] = 2.3
        greyed = mode == "grey" and not (e["new"] or e["emph"])
        if greyed:
            st = grey_style(st)
        under = HIGHLIGHT if (mode == "highlight" and (e["new"] or e["emph"])) else None
        draw_arrow(ax, a, b, e["arc"], st, under=under,
                   zorder=0.8 if greyed else 1.6, shrink=0.5 if e["off"] is not None else 1.5)
        mids[(e["src"], e["dst"], e["kind"])] = on_arc(a, b, e["arc"], 0.5)
        if not e["label"]:
            continue
        if e["kind"] in ("chem", "coll", "ct"):
            x, y = on_arc(a, b, e["arc"], e["t"])
            text_halo(ax, x, y, e["label"], GREY_TEXT if greyed else st["color"])
        elif e["kind"] in ("euv", "uv") and e["off"] is not None:
            x, y = on_arc(a, b, 0.0, 0.5)
            if abs(a[0] - b[0]) < 1e-9:
                sg = -1.0 if e["off"] < 0 else 1.0
                ax.text(x + 0.12 * sg, y, e["label"], color=PURPLE, fontsize=FS_LABEL,
                        fontweight="bold", ha="right" if sg < 0 else "left",
                        va="center")
            else:
                sg = -1.0 if e["off"] < 0 else 1.0
                ax.text(x, y + 0.1 * sg, e["label"], color=PURPLE, fontsize=FS_LABEL,
                        fontweight="bold", ha="center",
                        va="top" if sg < 0 else "bottom")
        elif e["kind"] == "euv":
            x, y = on_arc(a, b, e["arc"], e["lab_t"])
            nx, ny = bulge_normal(a, b, e["arc"])
            text_halo(ax, x + e["lab_off"] * nx, y + e["lab_off"] * ny, e["label"],
                      PURPLE, fontweight="bold")
    return mids


def draw_partners(ax, mids, pos, grey):
    for src, dst, partner in PARTNERS:
        m = mids.get((src, dst, "chem"))
        if m is None:
            continue
        start = box_edge(pos[partner], m)
        c = GREY if grey else PARTNER
        ax.plot([start[0], m[0]], [start[1], m[1]], ls=(0, (1, 2)), lw=0.9,
                color=c, zorder=0.7)
        ax.plot(*m, "o", ms=2.5, color=c, zorder=0.9)


def draw_box(ax, name, x, y, mode, gs_order, label=None):
    w, h = box_size(name)
    ghost = name in GHOSTS
    new = name in NEW_SPECIES
    greyed = mode == "grey" and not new
    if greyed:
        fc = "white" if ghost else GREY_FILL
        ec, tc = GREY, GREY_TEXT
    else:
        fc = "white" if ghost else ("#E3F2DC" if name in MOLECULES else "#DCEBF7")
        ec = "#7A7A7A" if ghost else BLACK
        tc = "#5A5A5A" if ghost else "black"
    ls = "--" if ghost else "-"
    lw = 2.2 if name in LIMITERS else 0.9
    if mode == "highlight" and new:
        ax.add_patch(FancyBboxPatch((x - w / 2 - 0.13, y - h / 2 - 0.13),
                                    w + 0.26, h + 0.26,
                                    boxstyle="round,pad=0.02,rounding_size=0.18",
                                    fc=HIGHLIGHT, ec="none", zorder=1.4))
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=ec, ls=ls, lw=lw, zorder=2))
    ax.text(x, y, label or TEX.get(name, name), ha="center", va="center",
            fontsize=FS_SPECIES - (0.5 if "_high" in name else 0),
            style="italic" if ghost else "normal", color=tc, zorder=4)
    if gs_order and name in gs_order:
        ax.text(x - w / 2 + 0.03, y + h / 2 - 0.01, str(gs_order.index(name) + 1),
                fontsize=5.5 if len(gs_order) < 10 else 4.8, ha="center",
                va="center", color="white", zorder=5,
                bbox=dict(boxstyle="circle,pad=0.15",
                          fc=GREY if greyed else BLACK, ec="none"))


def draw_note_box(ax, x, y, w, h, text, mode, fs=FS_NOTE, highlight=False, **kw):
    if highlight and mode == "highlight":
        ax.add_patch(FancyBboxPatch((x - w / 2 - 0.13, y - h / 2 - 0.13),
                                    w + 0.26, h + 0.26,
                                    boxstyle="round,pad=0.02,rounding_size=0.18",
                                    fc=HIGHLIGHT, ec="none", zorder=1.4))
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#FAFAFA", ec=DARK, ls=(0, (1, 1.5)), lw=0.9,
                                zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4,
            linespacing=1.25, **kw)


def draw_optional(ax, x, y, w, mode, fs=FS_NOTE - 0.5):
    """Dotted grey placeholder ladders: other elements can be added as O, S, N.

    configure.py accepts one set, --photchem_ions=O3,S3,N3 (photchem_ion_sets);
    another element or stage count needs its species rows, rates and update
    step in the code, so these are drawn without reactions.
    """
    h = 0.42
    bw = 0.78
    rows = [(el, [el, el + "$^+$", "…"]) for el in OPTIONAL]
    top = y + 0.55
    ax.text(x, top + 0.2, "other elements:\ncan be added", ha="center", va="center",
            fontsize=fs, color=GREY_TEXT, style="italic", linespacing=1.15)
    x0 = x - (bw + 0.3)
    for i, (el, names) in enumerate(rows):
        yy = top - 0.6 - i * 0.55
        for j, n in enumerate(names):
            xx = x0 + j * (bw + 0.3)
            if n != "…":
                ax.add_patch(FancyBboxPatch((xx - bw / 2, yy - h / 2), bw, h,
                                            boxstyle="round,pad=0.02,rounding_size=0.08",
                                            fc="white", ec=GREY, ls=(0, (1, 1.5)), lw=0.8,
                                            zorder=2))
            ax.text(xx, yy, n, ha="center", va="center", fontsize=fs, color=GREY_TEXT,
                    zorder=4)
            if j:
                ax.plot([xx - bw / 2 - 0.27, xx - bw / 2 - 0.03], [yy, yy],
                        ls=(0, (1, 1.5)), lw=0.7, color=GREY, zorder=1)


# ---------------------------------------------------------------------------
# legend


def legend_axes(fig, rect):
    lax = fig.add_axes(rect)
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    return lax


def draw_legend(fig, lax, kinds, boxes, mode, ncol=3, nbcol=4, row=0.155):
    """Arrow kinds in ncol columns, box kinds in nbcol columns below; inches."""
    bbox = lax.get_position()
    w_in = bbox.width * fig.get_figwidth()
    h_in = bbox.height * fig.get_figheight()
    lax.set_xlim(0, w_in)
    lax.set_ylim(0, h_in)
    rows = -(-len(kinds) // ncol)
    col_w = w_in / ncol
    for i, kind in enumerate(kinds):
        c, r = divmod(i, rows)
        x, y = c * col_w, h_in - row * (r + 0.5)
        draw_arrow(lax, (x + 0.02, y), (x + 0.36, y), 0.0, dict(STYLE[kind]), shrink=0)
        lax.text(x + 0.42, y, LEGEND[kind], fontsize=FS_LEGEND, va="center")
    y0 = h_in - row * (rows + 0.35)
    bcol_w = w_in / nbcol
    brows = -(-len(boxes) // nbcol)
    for i, (kind, text) in enumerate(boxes):
        c, r = divmod(i, brows)
        x, y = c * bcol_w, y0 - row * (r + 0.5)
        if kind == "gs":
            lax.text(x + 0.17, y, "1", fontsize=5.5, ha="center", va="center",
                     color="white",
                     bbox=dict(boxstyle="circle,pad=0.15", fc=BLACK, ec="none"))
        elif kind == "partner":
            lax.plot([x + 0.02, x + 0.34], [y, y], ls=(0, (1, 2)), lw=0.9, color=PARTNER)
            lax.plot(x + 0.34, y, "o", ms=2.5, color=PARTNER)
        else:
            fc, ec, ls, lw = {
                "ion": ("#DCEBF7", BLACK, "-", 0.9),
                "mol": ("#E3F2DC", BLACK, "-", 0.9),
                "ghost": ("white", "#7A7A7A", "--", 0.9),
                "limit": ("white", BLACK, "-", 2.2),
                "grey": (GREY_FILL, GREY, "-", 0.9),
                "optional": ("white", GREY, (0, (1, 1.5)), 0.9),
                "new": ("#DCEBF7", BLACK, "-", 0.9),
            }[kind]
            if kind == "new":
                lax.add_patch(FancyBboxPatch((x - 0.01, y - 0.07), 0.36, 0.14,
                                             boxstyle="round,pad=0.01,rounding_size=0.04",
                                             fc=HIGHLIGHT, ec="none"))
            lax.add_patch(FancyBboxPatch((x + 0.04, y - 0.042), 0.26, 0.084,
                                         boxstyle="round,pad=0.01,rounding_size=0.03",
                                         fc=fc, ec=ec, ls=ls, lw=lw))
        lax.text(x + 0.42, y, text, fontsize=FS_LEGEND, va="center")


def legend_height(kinds, boxes, ncol=3, nbcol=4, row=0.155):
    return row * (-(-len(kinds) // ncol) + -(-len(boxes) // nbcol) + 0.5) + 0.06


CORE_KINDS = ["cr", "crphoto", "uv", "grain", "rec", "chem", "coll"]
ALL_KINDS = ["cr", "crphoto", "uv", "euv", "grain", "rec", "chem", "coll", "ct", "ots"]


def boxes_legend(version):
    b = [("ion", "ion"), ("mol", "molecule"), ("ghost", "ghost (conservation)"),
         ("limit", "limits the substep"), ("gs", "Gauss–Seidel order"),
         ("partner", "partner reactant")]
    if version == "greycore":
        b.insert(0, ("grey", "GOW17 core (grey)"))
    if version in ("colour", "2panel"):
        b.insert(0, ("new", "added in Tigris"))
    if version != "core":
        b += [("optional", "can be added")]
    return b


# ---------------------------------------------------------------------------
# figures

XE_CORE = (r"$x_e = x({\rm He^+}) + x({\rm C^+}) + x({\rm HCO^+}) + x({\rm H_3^+})"
           r" + x({\rm H_2^+}) + x({\rm H^+}) + x({\rm O^+}) + x({\rm Si^+})$"
           "   (charge neutrality; electrons are a ghost species)")
# species_table.hpp ChargeSum over species_rows (gow17_core.hpp:36-63): X_high
# counts the charge of its lowest stage; above temp_hot0 HigherIonExtraCharge adds
# the rest of its CIE mean charge (gow17_network.hpp:2355-2367, 3135-3150).
XE_TIGRIS = (r"$x_e = \sum_i q_i\,x_i$ over the 18 ions, X$^{\geq q+}$ at charge $q$"
             r"$\;+\;\sum_{\rm X}[\bar q_{\rm X}(T) - q]\,x({\rm X}^{\geq q+})"
             r"\;+\;w(T)\,\Delta x_{e,\rm CIE}$ in hot gas")
HIGH_TEXT = (r"X$^{\geq q+}$: every stage above the top tracked one, one species; "
             "its split among stages, recombination and cooling follow CIE at $T$.\n"
             "The ladders are a configure choice (--photchem_ions=O3,S3,N3 here: O, S, N "
             "to O$^{2+}$, S$^{2+}$, N$^{2+}$); other elements need rates and code.")
HOT_TEXT = ("Hot gas, $T > 2\\times10^4$ K, with weight $w(T)$ rising to 1 at "
            "$3.5\\times10^4$ K: He and metal cooling from the CIE table;\n"
            "$x_e$ adds the CIE electrons of He$^{2+}$ and of metals beyond the tracked ions; "
            "heating and grain-assisted recombination are scaled by $1-w$.")


def canvas(xlim, ylim, legend_h, foot_h, top_pad=0.04):
    xr, yr = xlim[1] - xlim[0], ylim[1] - ylim[0]
    diag_h = WIDTH_IN * yr / xr
    H = top_pad + diag_h + foot_h + legend_h
    fig = plt.figure(figsize=(WIDTH_IN, H))
    ax = fig.add_axes([0, (legend_h + foot_h) / H, 1, diag_h / H])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")
    fax = legend_axes(fig, [0.0, legend_h / H, 1, foot_h / H])
    lax = legend_axes(fig, [0.03, 0.0, 0.97, legend_h / H])
    return fig, ax, fax, lax


def footer(fax):
    fax.text(0.5, 0.88, XE_TIGRIS, ha="center", va="center", fontsize=FS_NOTE)
    fax.text(0.5, 0.56, HIGH_TEXT, ha="center", va="center", fontsize=FS_NOTE,
             linespacing=1.3)
    fax.text(0.5, 0.17, HOT_TEXT, ha="center", va="center", fontsize=FS_NOTE,
             linespacing=1.3)


def fig_core():
    SIZE_OVERRIDE.clear()
    boxes = boxes_legend("core")
    lh = legend_height(CORE_KINDS, boxes)
    fig, ax, fax, lax = canvas((-2.75, 17.85), (0.05, 9.25), lh, 0.3)
    mids = draw_edges(ax, CORE_EDGES, POS_CORE, "plain")
    draw_partners(ax, mids, POS_CORE, grey=False)
    for n, (x, y) in POS_CORE.items():
        draw_box(ax, n, x, y, "plain", GS_CORE)
    fax.text(0.5, 0.5, XE_CORE, ha="center", va="center", fontsize=FS_NOTE)
    draw_legend(fig, lax, CORE_KINDS, boxes, "plain")
    return fig


SIGN_FULL = {"O": -1.0, "S": 1.0, "N": -1.0, "C": 1.0, "Si": 1.0}
FOOT_H = 1.05
# The core O+ -> O (H2) arrow bows under the O <-> O+ ladder arrows.
O_H2_ARC = 1.0


def full_core_edges():
    return [dict(e, arc=O_H2_ARC) if (e["src"], e["dst"], e["kind"]) == ("O+", "O", "chem")
            else e for e in core_edges_tigris()]


def draw_full(ax, mode):
    edges = full_core_edges() + ION_GAS + NEUTRAL_CI + ladder_edges(SIGN_FULL)
    mids = draw_edges(ax, edges, POS_FULL, mode)
    draw_partners(ax, mids, POS_FULL, grey=(mode == "grey"))
    for n, (x, y) in POS_FULL.items():
        draw_box(ax, n, x, y, mode, GS_TIGRIS)
    draw_optional(ax, 17.2, 1.2, 3.4, mode)
    text_halo(ax, 10.9, -0.75, "O$\\leftrightarrow$O$^+$ charge transfer is in GOW17"
              " (2-body 27, 28); near-resonant", GREEN, fs=FS_NOTE)


def fig_full(mode):
    SIZE_OVERRIDE.clear()
    SIZE_OVERRIDE.update({n: (BOX_W, 1.25) for n in ("O", "O+", "O2+", "O_high")})
    version = "greycore" if mode == "grey" else "colour"
    boxes = boxes_legend(version)
    lh = legend_height(ALL_KINDS, boxes)
    fig, ax, fax, lax = canvas((-6.9, 22.5), (-1.0, 9.25), lh, FOOT_H)
    draw_full(ax, mode)
    footer(fax)
    draw_legend(fig, lax, ALL_KINDS, boxes, mode)
    return fig


def fig_2panel():
    SIZE_OVERRIDE.clear()
    xa, ya = (-2.75, 17.85), (0.05, 9.25)
    xb, yb = (-2.75, 17.85), (-4.3, 3.55)
    ha = WIDTH_IN * (ya[1] - ya[0]) / (xa[1] - xa[0])
    hb = WIDTH_IN * (yb[1] - yb[0]) / (xb[1] - xb[0])
    boxes = boxes_legend("2panel")
    legend_h, foot_h, gap = legend_height(ALL_KINDS, boxes), FOOT_H, 0.1
    H = ha + gap + hb + foot_h + legend_h + 0.04
    fig = plt.figure(figsize=(WIDTH_IN, H))
    axa = fig.add_axes([0, (legend_h + foot_h + hb + gap) / H, 1, ha / H])
    axb = fig.add_axes([0, (legend_h + foot_h) / H, 1, hb / H])
    for ax, xl, yl in ((axa, xa, ya), (axb, xb, yb)):
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.axis("off")
    # (a) core with the ionized-gas additions; O <-> O+ is drawn in (b)
    mids = draw_edges(axa, core_edges_tigris() + ION_GAS + NEUTRAL_CI, POS_CORE,
                      "highlight")
    draw_partners(axa, mids, POS_CORE, grey=False)
    for n, (x, y) in POS_CORE.items():
        draw_box(axa, n, x, y, "highlight", GS_TIGRIS)
    axa.text(-2.6, 9.15, "(a)", fontsize=9, fontweight="bold", va="top")
    axa.text(6.6, 0.6, "O ladder and O$\\leftrightarrow$O$^+$\ncharge transfer: see (b)",
             color=DARK, fontsize=FS_NOTE, style="italic", ha="center", va="center")
    # (b) ion ladders, left to right, ionization above each pair
    SIZE_OVERRIDE.update({n: (BOX_W, 1.25) for n in POS_LADDERS})
    pos = POS_LADDERS
    sign = {k: -1.0 for k in SIGN_FULL}
    draw_edges(axb, ladder_edges(sign), pos, "highlight")
    for n, (x, y) in pos.items():
        draw_box(axb, n, x, y, "highlight", GS_TIGRIS)
    axb.text(-2.6, 3.5, "(b)", fontsize=9, fontweight="bold", va="top")
    axb.text(-0.75, -4.0, "O, O$^+$, C$^+$ and Si$^+$ are the boxes of (a); "
             "O$\\leftrightarrow$O$^+$ charge transfer is in GOW17 (near-resonant)",
             fontsize=FS_NOTE, color=DARK, ha="left", va="center", style="italic")
    draw_optional(axb, 14.5, -2.45, 3.4, "highlight")
    fax = legend_axes(fig, [0.0, legend_h / H, 1, foot_h / H])
    footer(fax)
    lax = legend_axes(fig, [0.03, 0.0, 0.97, legend_h / H])
    draw_legend(fig, lax, ALL_KINDS, boxes, "highlight")
    return fig


def contact_sheet(paths, out):
    import matplotlib.image as mpimg
    imgs = [mpimg.imread(p) for p in paths]
    fig, axs = plt.subplots(1, len(imgs), figsize=(4.2 * len(imgs), 3.6))
    titles = ["V1  core", "V2  Tigris, core greyed", "V3  two panels",
              "V4  Tigris, full colour"]
    for ax, im, t in zip(axs, imgs, titles):
        ax.imshow(im)
        ax.set_title(t, fontsize=12)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = pathlib.Path(__file__).resolve().parent
    parser.add_argument("--outdir", default=str(here.parent / "figures"))
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    made = []
    for name, make in [("fig_network_core", fig_core),
                       ("fig_network_tigris_greycore", lambda: fig_full("grey")),
                       ("fig_network_2panel", fig_2panel),
                       ("fig_network_1panel_colour", lambda: fig_full("highlight"))]:
        fig = make()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{name}.{ext}", dpi=300)
        plt.close(fig)
        made.append(outdir / f"{name}.png")
        print("wrote", outdir / name, "{pdf,png}")
    contact_sheet(made, outdir / "fig_network_versions.png")
    print("wrote", outdir / "fig_network_versions.png")


if __name__ == "__main__":
    main()
