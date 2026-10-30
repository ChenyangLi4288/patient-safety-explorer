"""Build the static site payload in site/ (adapted from flecomet/neurips-explorer).

Inputs:  data/map.parquet, data/corpus.parquet, data/topics.json, data/emb_medcpt.npy,
         labels.json (hand-written topic / subtopic names)
Outputs: site/data.json     coordinates, titles, clusters, filters (loads first)
         site/details-people.json  authors, first affiliation        (loaded afterwards)
         site/details-topics.json  DOI, MeSH, neighbours             (loaded afterwards)
         site/index.html    templates/index.html with config values filled in

Abstracts are deliberately not published (publisher copyright); the panel links to PubMed.

Usage: python build_site.py
"""
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).parent
DATA, SITE = ROOT / "data", ROOT / "site"
TEMPLATE = ROOT / "templates" / "index.html"
N_NEIGHBORS = 5
MAX_AUTHORS = 8

CONFIG = {
    "NAME": "Patient Safety Literature",
    "YEAR": "2026",
    "REPO_URL": "https://github.com/ChenyangLi4288/patient-safety-explorer",
    "SOURCE_NAME": "PubMed",
    "SOURCE_URL": "https://pubmed.ncbi.nlm.nih.gov/",
    "STORAGE_PREFIX": "ps-explorer",
    "CSV_NAME": "patient-safety-saved.csv",
}

PERIODS = [(2000, 2004), (2005, 2009), (2010, 2014), (2015, 2019), (2020, 2022), (2023, 2024),
           (2025, 2100)]
# First match wins, so the more specific design comes first.
STUDY_TYPES = [
    ("Systematic review / meta-analysis", {"Systematic Review", "Meta-Analysis"}),
    ("Trial", {"Randomized Controlled Trial", "Clinical Trial", "Controlled Clinical Trial",
               "Pragmatic Clinical Trial", "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
               "Clinical Trial, Phase II", "Clinical Trial, Phase I", "Equivalence Trial"}),
    ("Review", {"Review", "Scoping Review"}),
    ("Guideline / consensus", {"Practice Guideline", "Guideline", "Consensus Development Conference"}),
    ("Observational / evaluation", {"Observational Study", "Multicenter Study", "Evaluation Study",
                                    "Comparative Study", "Validation Study"}),
    ("Case report", {"Case Reports"}),
    ("Commentary / editorial", {"Editorial", "Comment", "Letter", "News", "Interview"}),
    ("Book chapter", {"Book Chapter"}),
]


def period(year) -> str:
    if year is None or np.isnan(year):
        return "Unknown"
    for a, b in PERIODS:
        if a <= year <= b:
            return f"{a}–{str(b)[2:]}" if b < 2100 else f"{a}–"
    return "Unknown"


def study_type(pub_types) -> str:
    pts = set(pub_types)
    for name, keys in STUDY_TYPES:
        if pts & keys:
            return name
    return "Other article"


EMAIL = re.compile(r"(electronic address:\s*)?\S+@\S+", re.I)


def clean_aff(aff: str) -> str:
    # Drop author e-mail addresses that PubMed appends to affiliations.
    aff = EMAIL.sub("", aff or "").strip(" .,;")
    return aff[:100]


def intern(values, order=None):
    """Replace repeated strings by indices into a table, to keep the payload small."""
    table = order or sorted(set(values))
    index = {v: i for i, v in enumerate(table)}
    return table, [index[v] for v in values]


def neighbours(emb: np.ndarray, k: int) -> list[list[int]]:
    e = torch.from_numpy(emb).cuda()
    out = []
    for i in range(0, len(e), 4096):
        s = e[i:i + 4096] @ e.T
        s[torch.arange(s.shape[0]), torch.arange(i, i + s.shape[0])] = -2  # not yourself
        out.extend(s.topk(k, dim=1).indices.cpu().tolist())
    return out


def render_index(template: str) -> str:
    for key, value in CONFIG.items():
        template = template.replace(f"@@{key}@@", value)
    left = sorted(set(re.findall(r"@@\w+@@", template)))
    assert not left, f"unfilled placeholders in {TEMPLATE}: {left}"
    return template


def main():
    m = pd.read_parquet(DATA / "map.parquet")
    corpus = pd.read_parquet(DATA / "corpus.parquet")
    order = pd.Series(np.arange(len(corpus)), index=corpus["pmid"])
    rows = order.loc[m["pmid"]].values  # map rows -> corpus/embedding rows
    df = corpus.iloc[rows].reset_index(drop=True)
    assert (df["pmid"].values == m["pmid"].values).all()

    topics = json.loads((DATA / "topics.json").read_text(encoding="utf-8"))
    labels = json.loads((ROOT / "labels.json").read_text(encoding="utf-8"))
    x, y = m["x"].round(2).tolist(), m["y"].round(2).tolist()
    sub, top = m["subtopic"].astype(int).values, m["topic"].astype(int).values

    def centre(mask):
        # Median, not mean: robust to the stray points UMAP leaves around a cluster.
        return round(float(np.median(m["x"].values[mask])), 3), round(float(np.median(m["y"].values[mask])), 3)

    clusters, groups = {}, {}
    for s in topics["subtopics"]:
        mask = sub == s["id"]
        cx, cy = centre(mask)
        clusters[str(s["id"])] = {"name": labels["subtopics"][str(s["id"])], "n": int(mask.sum()),
                                  "cx": cx, "cy": cy, "parent": s["parent"]}
    for t in topics["topics"]:
        mask = top == t["id"]
        cx, cy = centre(mask)
        groups[str(t["id"])] = {"name": labels["topics"][str(t["id"])], "n": int(mask.sum()),
                                "cx": cx, "cy": cy}
    # Colour and browse subtopics by the topic most of their papers fall in.
    parent = [clusters[str(s)]["parent"] for s in sub]

    periods, period_idx = intern([period(v) for v in df["year"].astype(float)],
                                 order=[period(a) for a, _ in PERIODS] + ["Unknown"])
    types, type_idx = intern([study_type(p) for p in df["pub_types"]],
                             order=[n for n, _ in STUDY_TYPES] + ["Other article"])
    journals, journal_idx = intern(df["journal_abbr"].fillna("").tolist())
    mesh_table, _ = intern([t for ms in df["mesh_major"] for t in ms])
    mesh_index = {t: i for i, t in enumerate(mesh_table)}

    print("neighbours…", flush=True)
    emb = np.load(DATA / "emb_medcpt.npy")[rows]
    nn = neighbours(emb, N_NEIGHBORS)

    core = {
        "meta": {"name": CONFIG["NAME"], "n": len(df)},
        "clusters": clusters,
        "topics": groups,
        "periods": periods,
        "types": types,
        "journals": journals,
        "id": df["pmid"].tolist(),
        "x": x, "y": y,
        "c": sub.tolist(),
        "t": parent,
        "per": period_idx,
        "typ": type_idx,
        "jr": journal_idx,
        "yr": df["year"].fillna(0).astype(int).tolist(),
        "title": df["title"].str.replace(r"\s+", " ", regex=True).tolist(),
    }

    def authors(a):
        names = [x["name"] for x in a[:MAX_AUTHORS]]
        return ", ".join(names) + (", et al." if len(a) > MAX_AUTHORS else "")

    # Author lists live in the raw records (corpus.parquet keeps only the first author).
    raw = {}
    for p in sorted((DATA / "raw").glob("*.jsonl")):
        with p.open(encoding="utf-8") as f:
            for line in f:
                r = json.loads(line)
                raw.setdefault(r["pmid"], r["authors"])
    au = [raw.get(p, []) for p in df["pmid"]]
    details = {
        "authors": [authors(a) for a in au],
        "aff": [clean_aff(a[0]["aff"]) if a else "" for a in au],
        "doi": df["doi"].fillna("").tolist(),
        "mesh_table": mesh_table,
        "mesh": [[mesh_index[t] for t in ms] for ms in df["mesh_major"]],
        # Only the paper panel needs neighbours, so they stay out of the first payload.
        "nn": nn,
    }

    # Split details in two: GitHub warns on files over 50 MB.
    people = {k: details.pop(k) for k in ("authors", "aff")}
    SITE.mkdir(exist_ok=True)
    (SITE / "details.json").unlink(missing_ok=True)
    for name, payload in (("data.json", core), ("details-people.json", people),
                          ("details-topics.json", details)):
        path = SITE / name
        with path.open("w", encoding="utf-8") as f:
            json.dump(payload, f, separators=(",", ":"), ensure_ascii=False)
        print(f"wrote {path}: {os.path.getsize(path) / 1e6:.1f} MB")
    (SITE / "index.html").write_text(render_index(TEMPLATE.read_text(encoding="utf-8")), encoding="utf-8")
    print(f"wrote {SITE / 'index.html'}")


if __name__ == "__main__":
    main()
