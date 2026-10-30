"""Write a compact per-cluster summary (terms, MeSH, sample titles) for hand/LLM labelling.

Usage: python dump_for_labels.py  ->  data/label_input.md
"""
import json
from pathlib import Path

import pandas as pd

DATA = Path(__file__).parent / "data"
N_TITLES = 7


def main():
    m = pd.read_parquet(DATA / "map.parquet")
    c = pd.read_parquet(DATA / "corpus.parquet", columns=["pmid", "title"])
    m = m.merge(c, on="pmid")
    topics = json.loads((DATA / "topics.json").read_text(encoding="utf-8"))
    lines = []
    for t in sorted(topics["topics"], key=lambda t: t["id"]):
        lines.append(f"\n# T{t['id']} ({t['size']:,}) terms: {', '.join(t['terms'])} | "
                     f"mesh: {'; '.join(t['mesh'])}")
        for s in sorted((s for s in topics["subtopics"] if s["parent"] == t["id"]),
                        key=lambda s: -s["size"]):
            titles = m.loc[m["subtopic"] == s["id"], "title"]
            titles = titles.sample(min(N_TITLES, len(titles)), random_state=0)
            lines.append(f"## S{s['id']} ({s['size']:,}, med yr {s['median_year']:.0f}) "
                         f"terms: {', '.join(s['terms'])} | mesh: {'; '.join(s['mesh'][:4])}")
            lines.extend(f"- {x[:95]}" for x in titles)
    out = DATA / "label_input.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(out, len(lines), "lines")


if __name__ == "__main__":
    main()
