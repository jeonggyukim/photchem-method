#!/usr/bin/env python
"""Talk versions of the GOW17 core chemical network (network_gow17.py, version core).

    network_gow17_talk_all  all reactions by kind; CR (vermillion) and FUV (purple) thick
    network_gow17_talk_ncr  what the TIGRESS-NCR chemical network keeps from GOW17

Sized for a TV screen: 13 in wide, species 20 pt, legend and note >= 18 pt.

NCR content from Kim, J.-G. et al. 2023, ApJS 264, 10, Sections 3.1-3.3 and
Table 2, cross-checked with tigris-ncr/src/photchem/ncr_rates.hpp (HIIRates,
H2Rates, GetxCII, GetxeMol, CalculateFreeElectron, CalculateCOAbundance):
    H2   time dependent, Eq. (15): grain formation, CR (1.65 zeta_cr,H2),
         FUV photodissociation, collisional dissociation. H2+ is not tracked;
         the CR arrow is redrawn as H2 -> H.
    H+   time dependent, Eq. (20): CR, collisional and LyC ionization of H;
         radiative and grain-assisted recombination.
    C+   equilibrium, Eqs. (22)-(24): CR (3.85 xi_cr) and CR-induced
         (520 * 2 x_H2 xi_cr) ionization, FUV photoionization, grain-assisted
         and radiative recombination, C+ + H2 (the sum of GOW17 2-body 5 and 22).
    O+   x_O+ / x_O,tot = x_H+ (Section 3.2.2), the near-resonant charge
         exchange O + H+ <-> O+ + H (GOW17 2-body 27, 28).
    CO   fit, Eqs. (25)-(26); no explicit reactions.
    H, C, O by conservation; x_e = x_H+ + x_C+ + x_O+ + x_MH+, Eq. (27).

Usage:
    network_gow17_talk.py [--outdir ../figures]
"""
import argparse
import pathlib

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

import network_gow17 as fn  # noqa: E402

WIDTH_IN = 13.0
FS_SPECIES, FS_LABEL, FS_LEGEND, FS_NOTE = 20, 15, 20, 18
FS_REACTANT = 15
# Labels left off: the collisional H2 -> H partners (H, H2) crowd the H-H2 gap.
NO_REACTANT = {("H2", "H", "coll")}
# Reactant labels moved off neighbouring arrows (data units).
REACTANT_SHIFT = {("H2", "H", "coll"): (0.0, -0.33), ("C+", "C", "chem"): (0.45, 0.0),
                  ("OHx", "CO", "chem"): (0.0, 0.55)}
BOX = (1.75, 0.95)
XLIM, YLIM = (-2.9, 18.1), (0.0, 9.6)

CR_RED = fn.VERMILLION
CR_LW = 1.1 * 2.5
LIGHT, LIGHT_TEXT = "#CFCFCF", "#B5B5B5"
DROP, DROP_TEXT = "#E6E6E6", "#C4C4C4"
KEPT = "#262626"
CR_KINDS = ("cr", "crphoto")

NCR_SPECIES = {"H", "H2", "H+", "C", "C+", "CO", "O", "O+"}
NCR_EDGES = {
    ("H", "H+", "cr"), ("C", "C+", "cr"), ("C", "C+", "crphoto"),
    ("C", "C+", "uv"), ("H2", "H", "uv"),
    ("H", "H2", "grain"), ("H+", "H", "grain"), ("C+", "C", "grain"),
    ("C+", "C", "rec"), ("H+", "H", "rec"),
    ("C+", "C", "chem"),
    ("H+", "O+", "ct"), ("O+", "H+", "ct"),
    ("H2", "H", "coll"), ("H", "H+", "coll"),
}
NCR_H2_CR = fn.E("H2", "H", "cr", "CR", -0.62)
# Photoionization by ionizing (LyC) radiation, in NCR only (Kim J.-G. et al. 2023,
# Eqs. 19-21): H -> H+, and H2 -> H2+ taken to give 2H at once, drawn as H2 -> H.
LYC = "#5E3C99"
NCR_PHOT = [fn.E("H", "H+", "lyc", "", -0.5), fn.E("H2", "H", "lyc", "", -0.32)]

NCR_NOTE = (r"NCR: H$_2$, H$^+$ time-dependent;  C$^+$ in equilibrium;  CO from a fit;"
            r"  $x({\rm O^+}) = x({\rm H^+})\,x_{\rm O,tot}$"
            "\n"
            r"H, C, O by conservation;  "
            r"$x_e = x({\rm H^+}) + x({\rm C^+}) + x({\rm O^+}) + x({\rm MH^+})$,"
            r"  MH$^+$ = molecular ions")


# The core O+ -> H+ arrow runs under the H3+ box; bow it below, beside H+ -> O+.
ARC_FIX = {("H+", "O+", "chem"): 0.08, ("O+", "H+", "chem"): -0.2}


# GOW17 cosmic-ray rates scaled by 2 x(H2) are CR-induced UV photoreactions
# (AthenaK gow17.hpp:1924-1937): C -> C+ (520, beside the direct 3.85), CO -> C + O
# (92) and Si -> Si+ (8400). network_gow17 draws only CO as CR-induced UV; here Si is
# redrawn dashed and C gets a second, dashed arrow for its induced part.
CR_INDUCED = [fn.E("C", "C+", "crphoto", "", 0.55)]


def core_edges():
    out = []
    for e in fn.CORE_EDGES:
        e = dict(e, arc=ARC_FIX.get((e["src"], e["dst"], e["kind"]), e["arc"]))
        if (e["src"], e["dst"], e["kind"]) == ("Si", "Si+", "cr"):
            e["kind"] = "crphoto"
        # (2body 27, 28) H+ + O <-> O+ + H: charge transfer, which ties x(O+) to x(H+)
        if (e["src"], e["dst"]) in (("H+", "O+"), ("O+", "H+")) and e["kind"] == "chem":
            e["kind"] = "ct"
        out.append(e)
    return out + [dict(e) for e in CR_INDUCED]


def arrow(ax, a, b, arc, color, lw, ls="-", head=(7.0, 3.5), zorder=1.0, shrink=2.0):
    common = dict(connectionstyle=f"arc3,rad={arc}", shrinkA=shrink, shrinkB=shrink)
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", lw=lw, ls=ls, color=color,
                                 capstyle="butt", zorder=zorder, **common))
    ax.add_patch(FancyArrowPatch(
        a, b, arrowstyle=f"-|>,head_length={head[0]},head_width={head[1]}",
        lw=0, color=color, zorder=zorder, **common))


def cr_style(kind):
    ls = "-" if kind == "cr" else (0, (3.0, 1.4))
    return dict(color=CR_RED, lw=CR_LW, ls=ls, head=(11.0, 5.0), zorder=2.0)


def halo(ax, x, y, s, color, fs, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, fontsize=fs, color=color, zorder=5,
            path_effects=[patheffects.withStroke(linewidth=4, foreground="white")], **kw)


def draw_edges(ax, edges, style_of, label_of):
    for e in edges:
        p, q = fn.POS_CORE[e["src"]], fn.POS_CORE[e["dst"]]
        a = fn.box_edge(p, q, *BOX)
        b = fn.box_edge(q, p, *BOX)
        arrow(ax, a, b, e["arc"], **style_of(e))
        lab = label_of(e)
        if lab:
            x, y = fn.on_arc(a, b, e["arc"], lab["t"])
            if e["arc"]:
                nx, ny = fn.bulge_normal(a, b, e["arc"])
            else:
                nx, ny = 0.0, 0.0
            halo(ax, x + lab["off"] * nx, y + lab["off"] * ny, lab["text"],
                 CR_RED, FS_LABEL, fontweight="bold")


def draw_box(ax, name, style):
    x, y = fn.POS_CORE[name]
    w, h = BOX
    ghost = name in fn.GHOSTS
    if style == "drop":
        fc, ec, tc, lw = "white", DROP, DROP_TEXT, 1.2
    else:
        fc = "white" if ghost else ("#E3F2DC" if name in fn.MOLECULES else "#DCEBF7")
        ec = "#7A7A7A" if ghost else fn.BLACK
        tc = "#5A5A5A" if ghost else "black"
        lw = 1.4
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.14",
                                fc=fc, ec=ec, ls="--" if ghost else "-", lw=lw,
                                zorder=3))
    ax.text(x, y, fn.TEX.get(name, name), ha="center", va="center",
            fontsize=FS_SPECIES, style="italic" if ghost else "normal", color=tc,
            zorder=4)


def legend(lax, entries, gap=0.45, arrow_w=0.6):
    """One centred row of arrow entries, (style dict, text); axes in inches."""
    fig = lax.figure
    bb = lax.get_position()
    w_in, h_in = bb.width * fig.get_figwidth(), bb.height * fig.get_figheight()
    lax.set_xlim(0, w_in)
    lax.set_ylim(0, h_in)
    y = 0.5 * h_in
    r = fig.canvas.get_renderer()
    texts, widths = [], []
    for st, text in entries:
        t = lax.text(0, y, text, fontsize=FS_LEGEND, va="center")
        texts.append(t)
        widths.append(t.get_window_extent(r).width / fig.dpi)
    total = sum(arrow_w + 0.15 + w for w in widths) + gap * (len(entries) - 1)
    x = 0.5 * (w_in - total)
    for (st, _), t, w in zip(entries, texts, widths):
        arrow(lax, (x, y), (x + arrow_w, y), 0.0, shrink=0, **st)
        t.set_x(x + arrow_w + 0.15)
        x += arrow_w + 0.15 + w + gap


def canvas(legend_h, note_h):
    diag_h = WIDTH_IN * (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])
    H = diag_h + note_h + legend_h
    fig = plt.figure(figsize=(WIDTH_IN, H))
    ax = fig.add_axes([0, (legend_h + note_h) / H, 1, diag_h / H])
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    nax = fn.legend_axes(fig, [0, legend_h / H, 1, note_h / H]) if note_h else None
    lax = fn.legend_axes(fig, [0, 0, 1, legend_h / H])
    return fig, ax, nax, lax


CR_LABEL = {
    ("H2", "H2+"): dict(text="CR", t=0.5, off=0.0),
    ("He", "He+"): dict(text="CR", t=0.5, off=0.42),
    ("H", "H+"): dict(text="CR", t=0.5, off=0.42),
    ("C", "C+"): dict(text="CR", t=0.5, off=0.42),
    ("CO", "C"): dict(text="CR", t=0.5, off=0.4),
    ("CO", "HCO+"): dict(text="CR", t=0.5, off=0.42),
    ("Si", "Si+"): dict(text="CR", t=0.5, off=0.42),
}


def reactant_label(ax, e, color, extra=None):
    """Partner reactant on an ion-molecule or collisional arrow, as in network_gow17."""
    if e["kind"] not in ("chem", "coll", "ct") or not e["label"]:
        return
    if (e["src"], e["dst"], e["kind"]) in NO_REACTANT:
        return
    p, q = fn.POS_CORE[e["src"]], fn.POS_CORE[e["dst"]]
    a, b = fn.box_edge(p, q, *BOX), fn.box_edge(q, p, *BOX)
    x, y = fn.on_arc(a, b, e["arc"], e["t"])
    key = (e["src"], e["dst"], e["kind"])
    dx, dy = (extra or {}).get(key, REACTANT_SHIFT.get(key, (0.0, 0.0)))
    ax.text(x + dx, y + dy, e["label"], fontsize=FS_REACTANT, color=color, ha="center",
            va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=1.0))


def fig_cr():
    fig, ax, _, lax = canvas(0.75, 0.0)
    grey = dict(color=LIGHT, lw=1.0, head=(6.0, 3.0), zorder=1.0)
    edges = sorted(core_edges(), key=lambda e: e["kind"] in CR_KINDS)

    def style_of(e):
        return cr_style(e["kind"]) if e["kind"] in CR_KINDS else grey

    def label_of(e):
        reactant_label(ax, e, "#8A8A8A")
        return None

    draw_edges(ax, edges, style_of, label_of)
    for n in fn.POS_CORE:
        draw_box(ax, n, "plain")
    legend(lax, [(cr_style("cr"), "cosmic-ray ionization"),
                 (cr_style("crphoto"), "CR-induced UV photoreaction"),
                 (dict(color=LIGHT, lw=1.6, head=(8.0, 4.0)), "other reactions")])
    return fig


FUV_PURPLE = fn.PURPLE


def fuv_style():
    return dict(color=FUV_PURPLE, lw=CR_LW, ls=(0, (3.0, 1.4)), head=(11.0, 5.0),
                zorder=1.9)


def thin_style(kind):
    st = fn.STYLE[kind]
    return dict(color=st["color"], lw=1.3, ls=st["ls"], head=(7.0, 3.5), zorder=1.0)


FUV_LABEL = {
    ("H2", "H"): dict(text="FUV", t=0.5, off=0.0),
    ("C", "C+"): dict(text="FUV", t=0.5, off=-0.42),
    ("CHx", "C"): dict(text="FUV", t=0.5, off=0.0),
    ("CO", "C"): dict(text="FUV", t=0.6, off=0.0),
    ("OHx", "O"): dict(text="FUV", t=0.5, off=-0.42),
    ("Si", "Si+"): dict(text="FUV", t=0.5, off=-0.42),
}


def fig_all():
    """All GOW17 core reactions by kind; cosmic-ray (vermillion) and FUV (purple)
    reactions thick, the rest thin in their own colours."""
    # same canvas as fig_ncr (its note strip left blank) so the boxes do not move
    # between the two slides
    fig, ax, _, lax = canvas(1.25, 0.95)
    rank = {"cr": 3, "crphoto": 3, "uv": 2}
    edges = sorted(core_edges(), key=lambda e: rank.get(e["kind"], 0))

    def style_of(e):
        if e["kind"] in CR_KINDS:
            return cr_style(e["kind"])
        if e["kind"] == "uv":
            return fuv_style()
        return thin_style(e["kind"])

    def label_of(e):
        reactant_label(ax, e, fn.STYLE[e["kind"]]["color"])
        return None

    draw_edges(ax, edges, style_of, label_of)
    for n in fn.POS_CORE:
        draw_box(ax, n, "plain")
    fig.delaxes(lax)
    H = fig.get_figheight()
    top = fn.legend_axes(fig, [0, 0.62 / H, 1, 0.62 / H])
    bot = fn.legend_axes(fig, [0, 0.0, 1, 0.62 / H])
    legend(top, [(cr_style("cr"), "cosmic-ray ionization"),
                 (cr_style("crphoto"), "CR-induced UV"),
                 (fuv_style(), "FUV photoreaction")])
    legend(bot, [(thin_style("grain"), "grain"),
                 (thin_style("rec"), "recombination"),
                 (thin_style("chem"), "ion-molecule"),
                 (thin_style("coll"), "collisional"),
                 (thin_style("ct"), "charge transfer")], gap=0.25)
    return fig


def lyc_style():
    return dict(color=LYC, lw=CR_LW, ls="-", head=(11.0, 5.0), zorder=1.95)


def fig_ncr():
    fig, ax, nax, lax = canvas(1.25, 0.95)
    kept = dict(color=KEPT, lw=1.6, head=(8.0, 4.0), zorder=1.5)
    drop = dict(color=DROP, lw=1.0, head=(6.0, 3.0), zorder=0.8)
    edges = [e for e in core_edges() if (e["src"], e["dst"]) != ("H2", "H2+")]
    edges.append(NCR_H2_CR)
    edges.extend(NCR_PHOT)

    def is_kept(e):
        return (e is NCR_H2_CR or e in NCR_PHOT
                or (e["src"], e["dst"], e["kind"]) in NCR_EDGES)

    edges.sort(key=lambda e: (is_kept(e), e["kind"] in CR_KINDS))

    def style_of(e):
        """Kept arrows in the styles of fig_all, dropped ones faint grey."""
        if not is_kept(e):
            return drop
        if e["kind"] in CR_KINDS:
            return cr_style(e["kind"])
        if e["kind"] == "uv":
            return fuv_style()
        if e["kind"] == "lyc":
            return lyc_style()
        return thin_style(e["kind"])

    def label_of(e):
        if e is not NCR_H2_CR and e["kind"] != "lyc":
            # the NCR H2 -> H CR arrow bows below, under the collisional label
            col = fn.STYLE[e["kind"]]["color"] if is_kept(e) else "#B0B0B0"
            reactant_label(ax, e, col, extra={("H2", "H", "coll"): (0.0, -0.75)})
        return None

    draw_edges(ax, edges, style_of, label_of)
    for n in fn.POS_CORE:
        draw_box(ax, n, "plain" if n in NCR_SPECIES else "drop")
    nax.text(0.5, 0.5, NCR_NOTE, ha="center", va="center", fontsize=FS_NOTE,
             linespacing=1.5)
    fig.delaxes(lax)
    H = fig.get_figheight()
    top = fn.legend_axes(fig, [0, 0.62 / H, 1, 0.62 / H])
    bot = fn.legend_axes(fig, [0, 0.0, 1, 0.62 / H])
    legend(top, [(cr_style("cr"), "CR ionization"),
                 (cr_style("crphoto"), "CR-induced UV"),
                 (lyc_style(), "photoionization"),
                 (fuv_style(), "FUV")])
    legend(bot, [(thin_style("grain"), "grain"),
                 (thin_style("rec"), "recombination"),
                 (thin_style("chem"), "ion-molecule"),
                 (thin_style("coll"), "collisional"),
                 (thin_style("ct"), "charge transfer")], gap=0.25)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = pathlib.Path(__file__).resolve().parent
    parser.add_argument("--outdir", default=str(here.parent / "figures"))
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, make in [("network_gow17_talk_all", fig_all),
                       ("network_gow17_talk_ncr", fig_ncr)]:
        fig = make()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{name}.{ext}", dpi=200)
        plt.close(fig)
        print("wrote", outdir / name, "{pdf,png}")


if __name__ == "__main__":
    main()
