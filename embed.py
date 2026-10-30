"""Merge data/raw/*.jsonl into data/corpus.parquet and embed title+abstract with MedCPT.

Outputs:
  data/corpus.parquet   one row per PMID (deduplicated), metadata only
  data/emb_medcpt.npy   float16 [N, 768], L2-normalised, row order = corpus.parquet

Usage: python embed.py [--limit N] [--batch 64]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

DATA = Path(__file__).parent / "data"
MODEL = "ncbi/MedCPT-Article-Encoder"
MAX_LEN = 512


def load_corpus() -> pd.DataFrame:
    rows = []
    for p in sorted((DATA / "raw").glob("*.jsonl")):
        with p.open(encoding="utf-8") as f:
            rows.extend(json.loads(line) for line in f)
    df = pd.DataFrame(rows).drop_duplicates("pmid").reset_index(drop=True)
    df = df[df["abstract"].str.len() > 0].reset_index(drop=True)
    df["mesh_terms"] = df["mesh"].map(lambda ms: [m["d"] for m in ms])
    df["mesh_major"] = df["mesh"].map(lambda ms: [m["d"] for m in ms if m["major"]])
    df["first_author"] = df["authors"].map(lambda a: a[0]["name"] if a else "")
    df["first_aff"] = df["authors"].map(lambda a: a[0]["aff"] if a else "")
    df["n_authors"] = df["authors"].map(len)
    return df


@torch.inference_mode()
def embed(titles: list[str], abstracts: list[str], batch: int) -> np.ndarray:
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModel.from_pretrained(MODEL, dtype=torch.float16).cuda().eval()
    # Length-sorted batches waste far less padding.
    order = np.argsort([len(a) for a in abstracts])[::-1]
    out = np.zeros((len(titles), model.config.hidden_size), dtype=np.float16)
    for i in range(0, len(order), batch):
        idx = order[i:i + batch]
        enc = tok([[titles[j], abstracts[j]] for j in idx], truncation=True, padding=True,
                  max_length=MAX_LEN, return_tensors="pt").to("cuda")
        cls = model(**enc).last_hidden_state[:, 0, :].float()
        out[idx] = torch.nn.functional.normalize(cls, dim=-1).cpu().numpy().astype(np.float16)
        if (i // batch) % 200 == 0:
            print(f"  {i + len(idx):,}/{len(order):,}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="embed only the first N (smoke test)")
    ap.add_argument("--batch", type=int, default=64)
    args = ap.parse_args()
    df = load_corpus()
    print(f"corpus: {len(df):,} papers", flush=True)
    if args.limit:
        df = df.iloc[:args.limit]
    emb = embed(df["title"].tolist(), df["abstract"].tolist(), args.batch)
    suffix = f"_first{args.limit}" if args.limit else ""
    df.drop(columns=["mesh", "authors"]).to_parquet(DATA / f"corpus{suffix}.parquet")
    np.save(DATA / f"emb_medcpt{suffix}.npy", emb)
    print(f"saved {emb.shape}", flush=True)


if __name__ == "__main__":
    main()
