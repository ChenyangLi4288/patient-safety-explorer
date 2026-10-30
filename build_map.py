"""2D map + two-level topics from MedCPT embeddings.

Inputs:  data/corpus.parquet, data/emb_medcpt.npy  (from embed.py)
Outputs: data/map.parquet   pmid, x, y, topic, subtopic
         data/topics.json   per topic/subtopic: size, top title terms, top MeSH, centroid

Usage: python build_map.py [--k1 30] [--k2 200] [--suffix _first2000] [--drop drop_pmids.csv]
"""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import umap
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS, CountVectorizer

DATA = Path(__file__).parent / "data"
# Words that appear in nearly every cluster of this corpus and say nothing about it.
DOMAIN_STOP = {"patient", "patients", "care", "health", "healthcare", "quality", "safety", "study",
               "results", "methods", "hospital", "clinical", "medical", "using", "use", "based",
               "analysis", "review", "evaluation", "effect", "effects", "outcomes", "associated"}


def ctfidf_terms(texts: pd.Series, labels: np.ndarray, top: int = 6) -> dict[int, list[str]]:
    """Class-based TF-IDF (BERTopic style): one document per cluster."""
    docs = texts.groupby(labels).apply(" ".join)
    vec = CountVectorizer(ngram_range=(1, 2), min_df=5, max_features=200_000,
                          stop_words=list(ENGLISH_STOP_WORDS | DOMAIN_STOP))
    tf = vec.fit_transform(docs.values).astype(np.float64)
    tf = tf.multiply(1 / np.maximum(tf.sum(axis=1), 1)).tocsr()
    avg = np.asarray(tf.sum(axis=0)).ravel() / tf.shape[0]
    idf = np.log1p(1 / np.maximum(avg, 1e-12))
    score = tf.multiply(idf).tocsr()
    vocab = vec.get_feature_names_out()
    out = {}
    for row, c in enumerate(docs.index):
        r = score.getrow(row).toarray().ravel()
        out[int(c)] = [vocab[i] for i in r.argsort()[::-1][:top]]
    return out


def mesh_lift(mesh: pd.Series, labels: np.ndarray, top: int = 5, min_n: int = 10) -> dict[int, list[str]]:
    """Major MeSH headings most over-represented in each cluster relative to the corpus."""
    glob = Counter(t for ms in mesh for t in set(ms))
    n = len(mesh)
    out = {}
    for c in np.unique(labels):
        sub = mesh[labels == c]
        cnt = Counter(t for ms in sub for t in set(ms))
        scored = [(k / len(sub)) / (glob[t] / n) * np.log1p(k)
                  for t, k in cnt.items() if k >= min_n]
        terms = [t for t, k in cnt.items() if k >= min_n]
        out[int(c)] = [terms[i] for i in np.argsort(scored)[::-1][:top]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k1", type=int, default=30)
    ap.add_argument("--k2", type=int, default=200)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--drop", default="", help="CSV with a pmid column of papers to exclude")
    args = ap.parse_args()

    df = pd.read_parquet(DATA / f"corpus{args.suffix}.parquet")
    emb = np.load(DATA / f"emb_medcpt{args.suffix}.npy").astype(np.float32)
    if args.drop:
        # Off-topic papers found on a previous map (e.g. the air/water-quality cluster).
        drop = set(pd.read_csv(args.drop, dtype=str)["pmid"])
        keep = ~df["pmid"].isin(drop).values
        df, emb = df[keep].reset_index(drop=True), emb[keep]
        print(f"dropped {(~keep).sum():,} papers listed in {args.drop}", flush=True)
    print(f"{len(df):,} papers", flush=True)

    # Centre first: MedCPT vectors share a large common direction (mean cos ~0.68).
    emb -= emb.mean(axis=0)
    z = PCA(n_components=50, random_state=0).fit_transform(emb)
    z /= np.linalg.norm(z, axis=1, keepdims=True)
    print("pca done", flush=True)

    k1, k2 = min(args.k1, len(df) // 20), min(args.k2, len(df) // 10)
    topic = KMeans(k1, n_init=4, random_state=0).fit_predict(z)
    sub = KMeans(k2, n_init=2, random_state=0).fit_predict(z)
    # Each subtopic belongs to the topic most of its papers fall in.
    parent = pd.Series(topic).groupby(sub).agg(lambda s: s.mode()[0]).to_dict()
    print("kmeans done", flush=True)

    xy = umap.UMAP(n_neighbors=30, min_dist=0.1, metric="cosine", random_state=0,
                   low_memory=True, verbose=True).fit_transform(z)
    print("umap done", flush=True)

    text = df["title"].fillna("") + " " + df["keywords"].map(" ".join)
    topics = {"topics": [], "subtopics": []}
    for level, labels, key in ((1, topic, "topics"), (2, sub, "subtopics")):
        terms, mesh = ctfidf_terms(text, labels), mesh_lift(df["mesh_major"], labels)
        for c in sorted(terms):
            m = labels == c
            entry = {"id": c, "size": int(m.sum()), "terms": terms[c], "mesh": mesh[c],
                     "x": float(xy[m, 0].mean()), "y": float(xy[m, 1].mean()),
                     "median_year": float(np.nanmedian(df.loc[m, "year"].astype(float)))}
            if level == 2:
                entry["parent"] = int(parent[c])
            topics[key].append(entry)

    pd.DataFrame({"pmid": df["pmid"], "x": xy[:, 0].astype(np.float32),
                  "y": xy[:, 1].astype(np.float32), "topic": topic.astype(np.int16),
                  "subtopic": sub.astype(np.int16)}).to_parquet(DATA / f"map{args.suffix}.parquet")
    (DATA / f"topics{args.suffix}.json").write_text(json.dumps(topics, indent=1), encoding="utf-8")
    for t in sorted(topics["topics"], key=lambda t: -t["size"]):
        print(f"{t['id']:3d} {t['size']:7,d}  {', '.join(t['terms'][:4]):60s} | {'; '.join(t['mesh'][:2])}")


if __name__ == "__main__":
    main()
