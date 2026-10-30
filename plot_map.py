"""Static overview of the map: small multiples by topic family + a year panel.

Usage: python plot_map.py  ->  data/map_overview.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.colors import LinearSegmentedColormap

DATA = Path(__file__).parent / "data"
SURFACE, INK, INK2, MUTED, BG_PTS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
HIGHLIGHT = "#2a78d6"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

NAMES = {int(k): v for k, v in
         json.loads((Path(__file__).parent / "labels.json").read_text(encoding="utf-8"))["topics"].items()}
# Hand grouping of the 30 topics (ids as in labels.json) for the small multiples.
FAMILIES = {
    "Core patient safety": [27, 11, 17, 22, 5],
    "Infection & antimicrobials": [13, 25, 7, 0],
    "Drug & treatment harms": [10, 4, 14, 28],
    "Procedures & specialties": [21, 23, 16, 2, 12, 26],
    "Populations & settings": [19, 29, 15, 1, 18, 3, 9],
    "Health systems & QI": [8, 6, 20, 24],
}


def halo(t):
    t.set_path_effects([pe.withStroke(linewidth=3, foreground=SURFACE)])


def setup(ax, title, lim):
    ax.set_facecolor(SURFACE)
    ax.set_xlim(lim[0]); ax.set_ylim(lim[1])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, loc="left", fontsize=13, color=INK, fontweight="bold")


def main():
    m = pd.read_parquet(DATA / "map.parquet")
    years = pd.read_parquet(DATA / "corpus.parquet", columns=["pmid", "year"])
    m = m.merge(years, on="pmid", how="left")
    topics = {t["id"]: t for t in json.loads((DATA / "topics.json").read_text())["topics"]}
    x, y = m["x"].values, m["y"].values
    pad = 0.5
    lim = ((np.percentile(x, 0.2) - pad, np.percentile(x, 99.8) + pad),
           (np.percentile(y, 0.2) - pad, np.percentile(y, 99.8) + pad))
    rng = np.random.default_rng(0)
    order = rng.permutation(len(m))  # random draw order so no topic paints on top

    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig, axes = plt.subplots(2, 4, figsize=(28, 14.5), facecolor=SURFACE)
    axes = axes.ravel()

    def place_labels(ax, texts):
        # Labels sit at topic centroids; nudge them apart (leader lines show the anchor).
        adjust_text(texts, ax=ax, expand=(1.15, 1.3), only_move="xy",
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))

    # Panel 1: everything, labelled.
    ax = axes[0]
    setup(ax, f"All {len(m):,} papers, 30 topics", lim)
    ax.scatter(x[order], y[order], s=0.15, c=MUTED, alpha=0.25, lw=0, rasterized=True)
    texts = []
    for tid, t in topics.items():
        texts.append(ax.text(t["x"], t["y"], NAMES[tid].replace("\n", " "), fontsize=8,
                             color=INK, ha="center", va="center"))
        halo(texts[-1])
    place_labels(ax, texts)

    # Panels 2-7: one family highlighted over the grey corpus.
    for ax, (fam, ids) in zip(axes[1:7], FAMILIES.items()):
        sel = m["topic"].isin(ids).values
        n = int(sel.sum())
        setup(ax, f"{fam}  ({n:,} papers)", lim)
        ax.scatter(x[~sel], y[~sel], s=0.15, c=BG_PTS, alpha=0.18, lw=0, rasterized=True)
        ax.scatter(x[sel], y[sel], s=0.25, c=HIGHLIGHT, alpha=0.35, lw=0, rasterized=True)
        texts = []
        for tid in ids:
            t = topics[tid]
            texts.append(ax.text(t["x"], t["y"], f"{NAMES[tid].replace(chr(10), ' ')}\n{t['size']:,}",
                                 fontsize=9.5, color=INK, ha="center", va="center", linespacing=1.1))
            halo(texts[-1])
        place_labels(ax, texts)

    # Panel 8: where the recent literature sits (share of papers from 2020 on, per hex cell).
    ax = axes[7]
    setup(ax, "Share of papers published 2020 or later", lim)
    recent = (m["year"].fillna(0).values >= 2020).astype(float)
    cmap = LinearSegmentedColormap.from_list("blues", BLUES)
    hb = ax.hexbin(x, y, C=recent, reduce_C_function=np.mean, gridsize=90, mincnt=8,
                   cmap=cmap, vmin=0, vmax=0.8, linewidths=0, extent=(*lim[0], *lim[1]))
    cb = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.01, ticks=[0, 0.2, 0.4, 0.6, 0.8])
    cb.ax.set_yticklabels(["0%", "20%", "40%", "60%", "≥80%"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK2, labelsize=9, length=0)
    overall = recent.mean()
    ax.text(0.0, -0.02, f"Corpus-wide: {overall:.0%} of papers are from 2020 or later. "
            "Cells with < 8 papers hidden.", transform=ax.transAxes, fontsize=9, color=INK2, va="top")

    fig.suptitle("Patient safety & quality-of-care literature, PubMed 2000–2026  ·  MedCPT embeddings, UMAP, k-means (30)",
                 x=0.01, ha="left", fontsize=16, color=INK)
    fig.text(0.01, 0.005, "Nearby points are semantically similar papers. Topic names are provisional "
             "(from distinctive title terms and MeSH headings). Label number = papers in topic.",
             fontsize=10, color=INK2)
    fig.tight_layout(rect=(0, 0.015, 1, 0.97))
    out = DATA / "map_overview.png"
    fig.savefig(out, dpi=110, facecolor=SURFACE)
    print(out)


if __name__ == "__main__":
    main()
