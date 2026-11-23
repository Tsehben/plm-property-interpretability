#!/usr/bin/env python3
"""
Three figures for the final localisation results package.

Saved outputs only: every number is read from the committed result tables, none
is recomputed or re-estimated here. No model is run, no sample selected, P0 never
read.

Light and dark variants are rendered as separate, deliberately stepped palettes
(not an automatic flip), PNG at 200 dpi plus SVG.
"""

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "figures")

THEME = {
    "light": dict(surface="#fcfcfb", page="#f9f9f7", ink="#0b0b0b", ink2="#52514e",
                  muted="#898781", grid="#e1e0d9", axis="#c3c2b7",
                  s1="#2a78d6", s2="#eb6834", s3="#1baf7a", on_fill="#ffffff"),
    "dark": dict(surface="#1a1a19", page="#0d0d0d", ink="#f0efec", ink2="#c3c2b7",
                 muted="#898781", grid="#2c2c2a", axis="#383835",
                 s1="#3987e5", s2="#d95926", s3="#199e70", on_fill="#ffffff"),
}


def style(t):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"],
        "savefig.facecolor": t["surface"],
        "text.color": t["ink"], "axes.labelcolor": t["ink2"],
        "xtick.color": t["muted"], "ytick.color": t["muted"],
        "axes.edgecolor": t["axis"], "grid.color": t["grid"],
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 9.5,
    })


def finish(fig, name, t, mode):
    os.makedirs(OUT, exist_ok=True)
    for ext, dpi in (("png", 200), ("svg", None)):
        fig.savefig(os.path.join(OUT, f"{name}_{mode}.{ext}"),
                    dpi=dpi, bbox_inches="tight", facecolor=t["surface"])
    plt.close(fig)


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ============================================================ FIGURE 1
def figure1(mode):
    t = THEME[mode]
    style(t)
    rows = [r for r in read_tsv(os.path.join(
        ROOT, "candidate_discovery", "v2_redundancy", "signal_annotation_cells.tsv"))
        if r["feature_id"] == "1446"]
    want = ["Extracellular, annotated_present", "Cytoplasm, annotated_present",
            "Extracellular, not_annotated", "Cytoplasm, not_annotated"]
    d = {r["cell"]: r for r in rows}
    labels = ["Extracellular\nsignal annotated", "Cytoplasm\nsignal annotated",
              "Extracellular\nnot annotated", "Cytoplasm\nnot annotated"]
    colors = [t["s1"], t["s1"], t["s2"], t["s2"]]
    y = [3, 2, 1, 0]

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(11.4, 4.5), gridspec_kw=dict(width_ratios=[2.15, 1], wspace=0.42))

    for yy, cell, c in zip(y, want, colors):
        r = d[cell]
        p25, p75 = float(r["p25"]), float(r["p75"])
        med, mean = float(r["median"]), float(r["mean"])
        ax.plot([p25, p75], [yy, yy], color=c, lw=9, solid_capstyle="round",
                zorder=2, alpha=0.9)
        ax.plot([med, med], [yy - 0.21, yy + 0.21], color=t["surface"], lw=4, zorder=4)
        ax.plot([med, med], [yy - 0.19, yy + 0.19], color=t["ink"], lw=2, zorder=5)
        ax.plot([mean], [yy], marker="D", ms=7.5, mfc=t["surface"], mec=c,
                mew=2.2, zorder=6, clip_on=False)
        ax.text(max(p75, mean) + 0.22, yy, f"median {med:.2f}   mean {mean:.2f}",
                va="center", ha="left", fontsize=9, color=t["ink2"])

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.5, color=t["ink"])
    ax.set_xlim(-0.18, 8.6)
    ax.set_ylim(-0.62, 3.62)
    ax.set_xlabel("Feature 1446 — mean activation over the first 30 residues", color=t["ink2"])
    ax.xaxis.grid(True, lw=0.8); ax.set_axisbelow(True)
    ax.spines["left"].set_color(t["axis"])
    ax.set_title("Activation distribution", fontsize=10, color=t["ink2"], loc="left", pad=9)
    ax.legend(handles=[
        Patch(facecolor=t["s1"], label="Signal peptide annotated"),
        Patch(facecolor=t["s2"], label="Not annotated"),
        Line2D([0], [0], color=t["muted"], lw=8, alpha=0.45, label="bar = interquartile range"),
        Line2D([0], [0], color=t["ink"], lw=2, label="line = median"),
        Line2D([0], [0], color=t["muted"], marker="D", ms=7.5, mfc=t["surface"], mew=2.2,
               lw=0, label="open diamond = mean")],
        loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False,
        ncol=3, fontsize=8.8, labelcolor=t["ink"], columnspacing=1.4, handlelength=1.8)

    for yy, cell, c in zip(y, want, colors):
        r = d[cell]
        pv = 100 * float(r["prevalence_active"])
        bx.barh(yy, pv, height=0.52, color=c, zorder=3)
        bx.text(pv - 2.5 if pv > 22 else pv + 2.5, yy, f"{pv:.1f}%",
                va="center", ha="right" if pv > 22 else "left", fontsize=9.5,
                fontweight="bold", color=t["on_fill"] if pv > 22 else t["ink"], zorder=4)
        bx.text(102, yy, f"n = {int(r['n']):,}", va="center", ha="left",
                fontsize=9, color=t["ink2"])
    bx.set_yticks(y); bx.set_yticklabels([])
    bx.set_xlim(0, 100); bx.set_ylim(-0.62, 3.62)
    bx.set_xlabel("Proteins with any activity in the first 30 residues (%)", color=t["ink2"])
    bx.xaxis.grid(True, lw=0.8); bx.set_axisbelow(True)
    bx.spines["left"].set_color(t["axis"])
    bx.set_title("Prevalence and group size", fontsize=10, color=t["ink2"],
                 loc="left", pad=9)

    fig.suptitle("Feature 1446 is strongly associated with signal-peptide annotation, "
                 "including in the label-discordant groups",
                 x=0.085, y=1.055, ha="left", fontsize=13.5, fontweight="bold", color=t["ink"])
    fig.text(0.085, 0.99,
             "Discovery sample: 5,190 P1+P2 proteins. The 59 cytoplasmic proteins that carry a "
             "signal annotation activate 1446 about as strongly as the 1,252 extracellular ones.",
             ha="left", fontsize=9.5, color=t["ink2"])
    fig.text(0.085, -0.32,
             "“Not annotated” means no signal-peptide feature was recorded for the frozen "
             "sequence — it is NOT confirmed biological absence. Among the 59 label-discordant\n"
             "cytoplasmic proteins that do carry an annotation, 4 have experimental evidence "
             "(all 4 active, mean 4.72). This is an association measured on discovery data; it "
             "is not exclusive annotation tracking and identifies no mechanism.",
             ha="left", fontsize=8.8, color=t["muted"])
    finish(fig, "fig1_discovery_1446_annotation", t, mode)


# ============================================================ FIGURE 2
def figure2(mode):
    t = THEME[mode]
    style(t)
    rows = read_tsv(os.path.join(ROOT, "evaluation_p4", "v1", "primary_results.tsv"))
    order = [("2564", "0.0", "Feature 2564 · ablation (×0)"),
             ("2564", "2.0", "Feature 2564 · amplification (×2)"),
             ("3900", "0.0", "Feature 3900 · ablation (×0)"),
             ("3900", "2.0", "Feature 3900 · amplification (×2)")]
    d = {(r["candidate"], r["multiplier"]): r for r in rows}
    S = 1e5  # plot in units of 1e-5 nats

    fig, ax = plt.subplots(figsize=(11.4, 4.3))
    y = [3, 2, 1, 0]
    for yy, (c, m, lab) in zip(y, order):
        r = d[(c, m)]
        est = float(r["estimate"]) * S
        lo95, hi95 = float(r["ci95_lo"]) * S, float(r["ci95_hi"]) * S
        lo98, hi98 = float(r["ci9875_lo"]) * S, float(r["ci9875_hi"]) * S
        ax.plot([lo98, hi98], [yy, yy], color=t["s1"], lw=2, solid_capstyle="butt",
                zorder=2, alpha=0.95)
        for xend in (lo98, hi98):
            ax.plot([xend, xend], [yy - 0.11, yy + 0.11], color=t["s1"], lw=2, zorder=2)
        ax.plot([lo95, hi95], [yy, yy], color=t["s1"], lw=8, solid_capstyle="round",
                zorder=3, alpha=0.95)
        ax.plot([est], [yy], marker="o", ms=9, mfc=t["surface"], mec=t["s1"],
                mew=2.4, zorder=5)
        excl = r["ci9875_excludes_zero"] == "1"
        ax.text(hi98 + 1.2, yy + 0.17,
                f"{est:+.2f}   95% [{lo95:+.2f}, {hi95:+.2f}]", va="center",
                ha="left", fontsize=9, color=t["ink2"])
        ax.text(hi98 + 1.2, yy - 0.19,
                ("excludes zero at both levels" if excl else "interval includes zero"),
                va="center", ha="left", fontsize=8.8,
                fontweight="bold" if excl else "normal",
                color=t["ink"] if excl else t["muted"])
        ax.text(-9.5, yy, f"n = {r['n_proteins']} proteins\n{r['n_paired_positions']} paired cells",
                va="center", ha="right", fontsize=8.6, color=t["muted"])

    ax.axvline(0, color=t["ink"], lw=1.4, zorder=1, alpha=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([o[2] for o in order], fontsize=10, color=t["ink"])
    ax.set_xlim(-11, 86)
    ax.set_ylim(-0.6, 3.6)
    ax.set_xlabel("Paired KL difference  (candidate − matched controls), units of 10⁻⁵ nats",
                  color=t["ink2"])
    ax.xaxis.grid(True, lw=0.8); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=96)

    ax.legend(handles=[
        Line2D([0], [0], color=t["s1"], lw=8, label="95% BCa"),
        Line2D([0], [0], color=t["s1"], lw=2, label="98.75% BCa (Bonferroni, 4 contrasts)"),
        Line2D([0], [0], color=t["s1"], marker="o", ms=9, mfc=t["surface"], mew=2.4,
               lw=0, label="point estimate")],
        loc="upper right", frameon=False, fontsize=9, labelcolor=t["ink"],
        bbox_to_anchor=(1.0, 1.02))

    fig.suptitle("One of four internal intervention contrasts excludes zero — on the narrowest coverage in the design",
                 x=0.0, y=1.10, ha="left", fontsize=13.5, fontweight="bold", color=t["ink"])
    fig.text(0.0, 1.035,
             "P4 evaluation, 400 proteins (200 Extracellular / 200 Cytoplasm by design). "
             "Outcome is the model’s own masked-token distribution, not localisation.",
             ha="left", fontsize=9.5, color=t["ink2"])
    fig.text(0.0, -0.165,
             "Coverage of the 2564-ablation result: 61 of 400 proteins — 18.4% of the 332 "
             "candidate-active cells, the narrowest in the design. 55 of those 61 proteins are\n"
             "Extracellular and 50 of the 61 paired cells sit at residue 20. The loss is "
             "structural: matching an ablation’s norm needs a large negative control latent\n"
             "change, which removes exactly the cells where the candidate’s latent was largest. "
             "Median protein difference +0.23×10⁻⁵, with 37 of 61 proteins positive.",
             ha="left", fontsize=8.8, color=t["muted"])
    finish(fig, "fig2_p4_internal_interventions", t, mode)


# ============================================================ FIGURE 3
def figure3(mode):
    t = THEME[mode]
    style(t)
    iv = {r["contrast"]: r for r in read_tsv(os.path.join(
        ROOT, "sequence_variant_pilot", "v1_retrospective", "retrospective_intervals.tsv"))}
    ti = {r["contrast"]: r for r in read_tsv(os.path.join(
        ROOT, "sequence_variant_pilot", "v1_retrospective", "tie_decomposition.tsv"))}
    order = [("B-A", "1446 ×2  vs  unmodified sampling"),
             ("C-A", "1446 ×0  vs  unmodified sampling"),
             ("B-D", "1446 ×2  vs  random direction, matched norm"),
             ("C-E", "1446 ×0  vs  random direction, matched norm")]
    y = [3, 2, 1, 0]

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(12.4, 4.4), gridspec_kw=dict(width_ratios=[1.0, 1.0], wspace=0.30))

    for yy, (k, lab) in zip(y, order):
        r = iv[k]
        est, lo, hi = (float(r["estimate_mean"]), float(r["ci95_lo"]), float(r["ci95_hi"]))
        med = float(r["median_of_paired_differences"])
        ax.plot([lo, hi], [yy, yy], color=t["s1"], lw=8, solid_capstyle="round",
                zorder=3, alpha=0.95)
        ax.plot([est], [yy], marker="o", ms=9, mfc=t["surface"], mec=t["s1"],
                mew=2.4, zorder=5)
        ax.plot([med], [yy], marker="|", ms=16, color=t["ink"], mew=2.4, zorder=6)
        ax.text(hi + 0.035, yy, f"mean {est:+.3f}\n[{lo:+.3f}, {hi:+.3f}]",
                va="center", ha="left", fontsize=8.6, color=t["ink2"])
    ax.axvline(0, color=t["ink"], lw=1.4, zorder=1, alpha=0.55)
    ax.set_yticks(y); ax.set_yticklabels([o[1] for o in order], fontsize=9.5, color=t["ink"])
    ax.set_xlim(-0.40, 0.80); ax.set_ylim(-0.6, 3.6)
    ax.set_xlabel("Paired difference in WoLF PSORT margin  (extr − cyto, predictor score units)",
                  color=t["ink2"])
    ax.xaxis.grid(True, lw=0.8); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0, pad=4)
    ax.legend(handles=[
        Line2D([0], [0], color=t["s1"], lw=8, label="95% interval (exploratory)"),
        Line2D([0], [0], color=t["s1"], marker="o", ms=9, mfc=t["surface"], mew=2.4,
               lw=0, label="mean"),
        Line2D([0], [0], color=t["ink"], marker="|", ms=14, mew=2.4, lw=0,
               label="median of paired differences (all = 0)")],
        loc="upper center", bbox_to_anchor=(0.5, -0.20), frameon=False,
        fontsize=8.8, labelcolor=t["ink"])
    ax.set_title("All four contrasts: n = 40 proteins", fontsize=10, color=t["ink2"],
                 loc="left", pad=9)

    cats = [("n_tie_identical_sequences", t["s1"], "Identical sequence in both arms"),
            ("n_tie_distinct_sequences_same_score", t["s2"],
             "Different sequences, identical score"),
            ("n_not_tied", t["s3"], "Different sequences, different scores")]
    for yy, (k, lab) in zip(y, order):
        r = ti[k]
        left = 0.0
        for fld, col, _ in cats:
            v = int(r[fld])
            bx.barh(yy, v, left=left, height=0.52, color=col, zorder=3,
                    edgecolor=t["surface"], linewidth=2)
            if v >= 14:
                bx.text(left + v / 2, yy, str(v), va="center", ha="center",
                        fontsize=9.5, fontweight="bold", color=t["on_fill"], zorder=5)
            else:
                bx.plot([left + v / 2, left + v / 2], [yy + 0.27, yy + 0.44],
                        color=col, lw=1.2, zorder=4)
                bx.text(left + v / 2, yy + 0.47, str(v), va="bottom", ha="center",
                        fontsize=9, fontweight="bold", color=t["ink"], zorder=5)
            left += v
    bx.set_yticks(y); bx.set_yticklabels([])
    bx.set_xlim(0, 160); bx.set_ylim(-0.6, 3.6)
    bx.set_xlabel("Paired cells (40 proteins × 4 replicates = 160 per contrast)",
                  color=t["ink2"])
    bx.xaxis.grid(True, lw=0.8); bx.set_axisbelow(True)
    bx.spines["left"].set_visible(False); bx.tick_params(axis="y", length=0)
    bx.legend(handles=[Patch(facecolor=c, label=l) for _, c, l in cats],
              loc="upper center", bbox_to_anchor=(0.5, -0.20), frameon=False,
              ncol=1, fontsize=8.8, labelcolor=t["ink"], handlelength=1.6)
    bx.set_title("What the zero medians are made of", fontsize=10, color=t["ink2"],
                 loc="left", pad=9)

    fig.suptitle("The sequence pilot mostly did not change the sequence — so the zero medians are not an evaluator verdict",
                 x=0.0, y=1.10, ha="left", fontsize=13.5, fontweight="bold", color=t["ink"])
    fig.text(0.0, 1.035,
             "P3 pilot, 40 Metazoa proteins (20/20 by design), feature 1446, 4 replicates, "
             "10 fixed positions, shared uniforms across arms.",
             ha="left", fontsize=9.5, color=t["ink2"])
    fig.text(0.0, -0.46,
             "Intervals are RETROSPECTIVE and EXPLORATORY: not pre-registered, no error control, "
             "no multiplicity adjustment. An interval containing zero is not equivalence.\n"
             "Identical sequences across paired arms do NOT imply identical token distributions — "
             "the arms’ distributions differed at every edited position; the shared uniform simply\n"
             "fell in the same inverse-CDF bin. When the arms did diverge they differed by a median "
             "of 1 residue, and 24–33% of those divergent pairs still scored identically.",
             ha="left", fontsize=8.8, color=t["muted"])
    finish(fig, "fig3_sequence_pilot", t, mode)


if __name__ == "__main__":
    for mode in ("light", "dark"):
        figure1(mode)
        figure2(mode)
        figure3(mode)
    print("wrote:", sorted(os.listdir(OUT)))
