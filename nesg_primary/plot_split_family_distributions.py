"""Plot existing split-local mixed-family tables; never reopen test records.

Run with the pinned nesg-validation environment. Each family contributes once
to distribution panels. A multidomain protein may occur in several families;
unique-protein summaries come from the previously audited union counts.
"""

from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nesg-family-plot-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "results/split_mixed_label_families"
OUTPUT = ROOT / "figures/split_family_distributions"
NAVY, SOLUBLE, INSOLUBLE = "#20334E", "#7A65B5", "#D49A68"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    paths = [INPUT / "split_summary.tsv"] + [
        INPUT / f"{split}_all_GA_mixed_families.tsv"
        for split in ("train", "validation", "reserved_test")
    ]
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in paths}
    summary = pd.read_csv(paths[0], sep="\t")
    frames = {p.name.split("_all_GA")[0]: pd.read_csv(p, sep="\t") for p in paths[1:]}
    size_max = max(int(t.total_unique_proteins.max()) for t in frames.values())
    # Common, exhaustive integer-size bins, fixed independently of label ratio.
    edges = [1.5, 2.5, 3.5, 5.5, 10.5, 20.5, 40.5, max(80.5, size_max + 0.5)]
    bin_names = ["2", "3", "4–5", "6–10", "11–20", "21–40", "41–80"]
    assert size_max <= 80, "Update the final size-bin label if cohort sizes change."
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.color": NAVY,
            "axes.labelcolor": NAVY,
            "axes.edgecolor": "#D6DAE1",
            "xtick.color": NAVY,
            "ytick.color": NAVY,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
        }
    )
    stats, bins, memberships = [], [], []
    for split, title in [
        ("train", "Training"),
        ("validation", "Validation"),
        ("reserved_test", "Reserved test"),
    ]:
        table = (
            frames[split]
            .sort_values(
                ["total_unique_proteins", "pfam_accession"], ascending=[False, True]
            )
            .reset_index(drop=True)
        )
        row = summary[
            (summary.split == split) & (summary.annotation_policy == "all_GA")
        ].iloc[0]
        s = table.unique_soluble_proteins.to_numpy()
        i = table.unique_insoluble_proteins.to_numpy()
        n = table.total_unique_proteins.to_numpy()
        fraction = s / n
        assert table.pfam_accession.is_unique and len(table) == row.qualifying_families
        assert (s >= 1).all() and (i >= 1).all() and np.array_equal(s + i, n)
        assert np.isclose(fraction + i / n, 1).all()
        ranks = np.arange(1, len(table) + 1)
        table["split"] = split
        table["family_rank"] = ranks
        table["soluble_percent"] = 100 * fraction
        table["insoluble_percent"] = 100 * i / n
        table["soluble_to_insoluble_ratio"] = s / i
        memberships.append(table)
        counts, _ = np.histogram(n, edges)
        assert counts.sum() == len(table)
        for name, count in zip(bin_names, counts):
            bins.append(
                dict(split=split, proteins_per_family=name, families=int(count))
            )

        fig = plt.figure(figsize=(15, 11.5))
        grid = fig.add_gridspec(
            3,
            2,
            height_ratios=[1.1, 1.35, 1],
            left=0.075,
            right=0.97,
            bottom=0.18,
            top=0.86,
            hspace=0.56,
            wspace=0.24,
        )
        ax_size = fig.add_subplot(grid[0, 0])
        ax_ratio = fig.add_subplot(grid[0, 1])
        ax_counts = fig.add_subplot(grid[1, :])
        ax_percent = fig.add_subplot(grid[2, :], sharex=ax_counts)
        fig.suptitle(
            f"{title} · {len(table):,} mixed-label Pfam families",
            x=0.075,
            y=0.975,
            ha="left",
            fontsize=21,
            fontweight="bold",
        )
        fig.text(
            0.075,
            0.937,
            f"Whole split: {row.split_unique_proteins:,} unique proteins  |  "
            f"Soluble {row.split_soluble:,} ({100 * row.split_soluble / row.split_unique_proteins:.1f}%)  ·  "
            f"Insoluble {row.split_insoluble:,} ({100 * row.split_insoluble / row.split_unique_proteins:.1f}%)",
            fontsize=12,
        )
        fig.text(
            0.075,
            0.909,
            f"In these mixed families (counted once across families): "
            f"{row.unique_proteins_in_qualifying_families:,} proteins  |  "
            f"{row.soluble_in_qualifying_families:,} soluble : "
            f"{row.insoluble_in_qualifying_families:,} insoluble",
            fontsize=11,
        )

        bars = ax_size.bar(bin_names, counts, color="#9BACC0", width=0.73)
        ax_size.bar_label(bars, padding=3, fontsize=10)
        ax_size.set(
            title="How large are the families?",
            xlabel="Unique proteins per family",
            ylabel="Number of families",
            ylim=(0, max(counts) * 1.23),
        )
        ax_size.yaxis.set_major_locator(MaxNLocator(integer=True))
        ratio_counts, _, _ = ax_ratio.hist(
            100 * fraction, bins=np.arange(0, 101, 10), color=SOLUBLE, edgecolor="white"
        )
        assert int(ratio_counts.sum()) == len(table)
        ax_ratio.axvline(50, color=NAVY, ls="--", lw=1)
        ax_ratio.set(
            title="How soluble-dominated is each family?",
            xlabel="Soluble proteins within a family (%)",
            ylabel="Number of families",
            xlim=(0, 100),
        )
        ax_ratio.yaxis.set_major_locator(MaxNLocator(integer=True))

        ax_counts.bar(ranks, s, color=SOLUBLE, width=0.88, label="Soluble")
        ax_counts.bar(
            ranks, i, bottom=s, color=INSOLUBLE, width=0.88, label="Insoluble"
        )
        ax_counts.set(
            title="Every mixed-label family · ranked by total size",
            ylabel="Unique proteins",
            ylim=(0, size_max * 1.12),
        )
        ax_counts.tick_params(axis="x", labelbottom=False)
        ax_counts.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax_counts.legend(
            handles=[
                Patch(color=SOLUBLE, label="Soluble"),
                Patch(color=INSOLUBLE, label="Insoluble"),
            ],
            loc="upper right",
            frameon=False,
            ncol=2,
        )
        ax_percent.bar(ranks, 100 * fraction, color=SOLUBLE, width=0.88)
        ax_percent.bar(
            ranks, 100 * i / n, bottom=100 * fraction, color=INSOLUBLE, width=0.88
        )
        ax_percent.axhline(50, color=NAVY, ls="--", lw=0.8)
        ax_percent.set(
            title="Soluble : insoluble proportions · same families, same order",
            ylabel="Within-family proportion",
            xlabel="Family rank (one bar per family)",
            ylim=(0, 100),
            xlim=(0.3, len(table) + 0.7),
        )
        ax_percent.yaxis.set_major_formatter(PercentFormatter(100))
        ticks = (
            [1]
            + [
                k
                for k in range(20, len(table), 20)
                if len(table) - k >= max(4, 0.025 * len(table))
            ]
            + [len(table)]
        )
        ax_percent.set_xticks(ticks)

        top_lines = ["Largest families · exact soluble : insoluble counts"]
        for rank, r in enumerate(table.head(3).itertuples(), 1):
            top_lines.append(
                f"{rank}. {r.pfam_accession}  {r.biological_family_name}  —  "
                f"{r.unique_soluble_proteins} : {r.unique_insoluble_proteins}  "
                f"(total {r.total_unique_proteins})"
            )
        fig.text(
            0.075, 0.132, "\n".join(top_lines), fontsize=10, va="top", linespacing=1.5
        )
        footer = (
            "All GA-passing Pfam 37.0 assignments; ≥1 protein per class within this split. "
            "No training-support or ≥3/≥5 requirement.\n"
            "Each protein counts once per family; multidomain proteins can appear in several families. "
            "Competing overlaps retained; families are not full-protein architectures."
        )
        if split == "reserved_test":
            footer += "\nHeld-out family–label composition has been inspected; do not use it to revise datasets, model selection or hypotheses."
        fig.text(
            0.075,
            0.037 if split != "reserved_test" else 0.047,
            footer,
            fontsize=8.5,
            va="top",
            color="#566578",
        )
        for ext in ("png", "svg"):
            fig.savefig(OUTPUT / f"{split}_family_distribution.{ext}", dpi=200)
        plt.close(fig)
        stats.append(
            dict(
                split=split,
                mixed_families=len(table),
                min_family_size=int(n.min()),
                median_family_size=float(np.median(n)),
                max_family_size=int(n.max()),
                soluble_majority_families=int((s > i).sum()),
                equally_represented_families=int((s == i).sum()),
                insoluble_majority_families=int((s < i).sum()),
                median_family_soluble_percent=float(np.median(100 * fraction)),
                split_soluble=int(row.split_soluble),
                split_insoluble=int(row.split_insoluble),
                split_soluble_percent=100
                * row.split_soluble
                / row.split_unique_proteins,
                mixed_family_unique_soluble=int(row.soluble_in_qualifying_families),
                mixed_family_unique_insoluble=int(row.insoluble_in_qualifying_families),
                mixed_family_unique_soluble_percent=100
                * row.soluble_in_qualifying_families
                / row.unique_proteins_in_qualifying_families,
            )
        )
    result = pd.DataFrame(stats)
    result.to_csv(OUTPUT / "distribution_summary.tsv", sep="\t", index=False)
    pd.DataFrame(bins).to_csv(OUTPUT / "family_size_bins.tsv", sep="\t", index=False)
    pd.concat(memberships).to_csv(
        OUTPUT / "ranked_family_counts_and_ratios.tsv", sep="\t", index=False
    )
    assert hashes == {str(p.relative_to(ROOT)): sha256(p) for p in paths}
    (OUTPUT / "provenance.json").write_text(
        json.dumps(
            dict(
                input_sha256=hashes,
                script_sha256=sha256(Path(__file__)),
                matplotlib=matplotlib.__version__,
                pandas=pd.__version__,
                numpy=np.__version__,
                raw_test_records_opened=False,
                existing_test_label_report_used=True,
                definition="Split-local >=1 unique soluble and >=1 unique insoluble protein; all GA assignments",
                outputs_sha256={
                    p.name: sha256(p)
                    for p in sorted(OUTPUT.iterdir())
                    if p.suffix in (".png", ".svg", ".tsv")
                },
            ),
            indent=2,
        )
        + "\n"
    )
    print(result.to_string(index=False))
    print(
        "\nFamily-size distribution:\n"
        + pd.DataFrame(bins)
        .pivot(index="proteins_per_family", columns="split", values="families")
        .reindex(bin_names)
        .to_string()
    )


if __name__ == "__main__":
    main()
