"""AI / machine-learning papers in the corpus by year (docs/ai_trend.png).

A paper counts as AI / ML if its title or abstract mentions machine learning, deep learning,
artificial intelligence, natural language processing, neural networks or (large) language
models, or it has the MeSH heading Artificial Intelligence / Machine Learning (exploded terms
such as Deep Learning, Natural Language Processing and Large Language Models included).
Trials are papers with the publication type Randomized Controlled Trial; recent years are an
underestimate because MEDLINE indexing lags.

Usage: python plot_ai_trend.py  ->  data/ai_trend.png
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

DATA = Path(__file__).parent / "data"
SURFACE, INK, INK2, GRID, GREY, BLUE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#898781", "#2a78d6"
Y0, Y1 = 2005, 2025
AI_TEXT = re.compile(r"machine[- ]learning|deep[- ]learning|artificial intelligence|natural language processing|"
                     r"neural network|large language model|language model|\bLLMs?\b|\bChatGPT\b|\bGPT-4", re.I)
AI_MESH = {"Artificial Intelligence", "Machine Learning", "Deep Learning", "Natural Language Processing",
           "Neural Networks, Computer", "Supervised Machine Learning", "Unsupervised Machine Learning",
           "Large Language Models", "Generative Artificial Intelligence"}


def main():
    df = pd.read_parquet(DATA / "corpus.parquet", columns=["pmid", "title", "abstract", "year", "mesh_terms",
                                                          "pub_types"])
    df = df[df["pmid"].isin(pd.read_parquet(DATA / "map.parquet", columns=["pmid"])["pmid"])]  # papers on the map
    df = df[df["year"].between(Y0, Y1)]
    text = df["title"].fillna("") + " " + df["abstract"].fillna("")
    df["ai"] = text.str.contains(AI_TEXT) | df["mesh_terms"].map(lambda ms: bool(AI_MESH & set(ms)))
    df["rct"] = df["pub_types"].map(lambda ps: "Randomized Controlled Trial" in set(ps))
    by = df.groupby("year").agg(n=("pmid", "size"), ai=("ai", "sum"))
    by["rct"] = df[df["ai"]].groupby("year")["rct"].sum()
    by = by.fillna(0)
    by["share"] = 100 * by["ai"] / by["n"]
    n_ai, n_rct = int(by["ai"].sum()), int(by["rct"].sum())

    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 5.6), facecolor=SURFACE)
    for ax in (a1, a2):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#c3c2b7")
        ax.tick_params(colors=INK2, labelsize=12, length=0)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xticks(range(Y0, Y1 + 1, 5))
    a1.plot(by.index, by["share"], color=BLUE, lw=3, marker="o", ms=5)
    a1.set_ylim(0, by["share"].max() * 1.15)
    a1.set_ylabel("AI / ML papers, % of all papers", color=INK2, fontsize=12)
    a1.set_title("Share of papers about AI or machine learning", loc="left", fontsize=14, color=INK)
    a1.annotate(f"{by.loc[Y1, 'share']:.1f}% in {Y1}", (Y1, by.loc[Y1, "share"]), xytext=(-12, -4),
                textcoords="offset points", ha="right", va="center", fontsize=13, color=INK)
    a2.bar(by.index, by["ai"], color=GREY, width=0.72, label="All AI / ML papers")
    a2.bar(by.index, by["rct"], color=BLUE, width=0.72, label="Of which randomized trials (est.)")
    a2.set_ylabel("Papers per year", color=INK2, fontsize=12)
    a2.set_title(f"Few are trials: about {100 * n_rct / n_ai:.1f}% of {n_ai:,}", loc="left", fontsize=14,
                 color=INK)
    a2.legend(frameon=False, fontsize=12, loc="upper left")
    fig.suptitle(f"AI in patient-safety research, {Y0}–{Y1} (analysis corpus)", x=0.01, ha="left",
                 fontsize=16, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = DATA / "ai_trend.png"
    fig.savefig(out, dpi=120, facecolor=SURFACE)
    print(out, f"{n_ai:,} AI/ML papers, {n_rct} RCTs, {by.loc[Y1, 'share']:.1f}% in {Y1}")


if __name__ == "__main__":
    main()
