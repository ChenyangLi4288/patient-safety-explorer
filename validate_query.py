"""Check the corpus query against hand-labelled PubMed records.

Recall
  docs/query_validation/calibration_labels.csv  David W. Bates and Jason S. Adelman papers
      (2000+, with abstract), each labelled as patient-safety research (1) or not (0). These
      were used to tune the v2 query.
  docs/query_validation/holdout_labels.csv  Hardeep Singh (Houston), Tejal K. Gandhi and
      Gordon D. Schiff papers, not in the calibration set. Labelled after v2 was drafted and
      used only to check it. The diagnostic-safety text arm was then added after looking at the
      misses, so hold-out recall is reported both with and without that component.

Precision
  docs/query_validation/precision_samples.csv  random samples of query or component hits,
      one row per paper, labelled from title + major MeSH:
        Y  patient-safety research: the paper is mainly about harm caused by health care, or
           the risk of it, or how to prevent or measure that harm. This covers medical,
           medication and diagnostic errors; adverse events and ADEs of care (incidence,
           preventability, detection, trigger tools); iatrogenic injury, HAIs, falls and
           pressure injuries studied as safety outcomes, and their prevention (infection
           control, hand hygiene); CPOE, CDS, alerts, barcoding and smart pumps assessed for
           safety; wrong-patient and wrong-site errors; incident reporting, RCA and FMEA;
           safety culture and climate, speaking up, second victims; handoffs and communication
           failures; disclosure; malpractice-claim analyses of errors; safety indicators and
           measurement; deterioration and rescue (rapid response teams); medication
           reconciliation and inappropriate prescribing as medication safety; follow-up of test
           results.
        B  borderline: adjacent to safety but not mainly about it, e.g. pharmacovigilance
           signals for one drug, HAI epidemiology without a prevention angle, general quality
           measurement of complications, defensive medicine.
        N  not patient-safety research: general quality of care (access, guideline adherence,
           satisfaction, cost) without a harm focus; efficacy trials that report adverse events;
           test accuracy and inter-rater reliability; bench microbiology; occupational and
           environmental safety; clinical case reports.
      Precision is reported strictly (Y) and leniently (Y + B), with 95% Wilson intervals.
  docs/query_validation/sample_queries.json  the query each sample was drawn from, and the
      v2 component (or whole query) it is used to estimate.

Usage
  python validate_query.py                      recall + precision report, writes summary.json
  python validate_query.py --refresh-authors    also re-pull the author records and list
                                                new, unlabelled PMIDs
  python validate_query.py --draw NAME --source QUERY|QUERY_V1|<v2 component> [--n 50]
                                                draw a new random sample (empty labels) and
                                                append it to precision_samples.csv for labelling

Uses only public NCBI E-utilities (tool=ps_explorer, under 3 requests/s, no API key).
"""
import argparse
import calendar
import csv
import json
import math
import random
from datetime import date
from pathlib import Path

from count_corpus import (COMPONENTS_V2, HAS_ABSTRACT, QUERY_TERM, QUERY_V1_TERM, SINCE_2000, count)
from fetch_corpus import _request, fetch_records

DIR = Path(__file__).parent / "docs" / "query_validation"
CAL, HOLD = DIR / "calibration_labels.csv", DIR / "holdout_labels.csv"
SAMPLES, META, SUMMARY = DIR / "precision_samples.csv", DIR / "sample_queries.json", DIR / "summary.json"
LIMITS = f"{HAS_ABSTRACT} AND {SINCE_2000}"
QUERIES = {"QUERY": QUERY_TERM, "QUERY_V1": QUERY_V1_TERM, **COMPONENTS_V2}
POST_HOC = "tiab_diagnostic_safety"  # added after inspecting hold-out misses
# Author searches behind the label files (2000+, with abstract). Adelman JS also matches other
# people (e.g. James S. Adelman, Memphis/Warwick); only Jason S. Adelman's papers are labelled.
AUTHORS = {
    CAL: ['Bates DW[au]', 'Adelman JS[au] AND ("Adelman, Jason"[FAU] OR Montefiore[ad] OR Columbia[ad] OR '
          'NewYork-Presbyterian[ad] OR Einstein[ad])'],
    HOLD: ['"Singh, Hardeep"[FAU] AND (Houston[ad] OR Baylor[ad] OR DeBakey[ad])', '"Gandhi, Tejal K"[FAU]',
           '"Schiff, Gordon D"[FAU]'],
}


def esearch(term: str, retmax: int = 9999) -> tuple[int, list[str]]:
    for _ in range(5):
        raw = _request("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
        try:
            js = json.loads(raw)["esearchresult"]
            return int(js["count"]), js.get("idlist", [])
        except (ValueError, KeyError):
            continue
    raise RuntimeError("esearch failed repeatedly")


def members(query: str, pmids, chunk: int = 200) -> set[str]:
    """PMIDs (of those given) that the query retrieves, with the corpus limits applied."""
    pmids, out = sorted(set(pmids)), set()
    for i in range(0, len(pmids), chunk):
        ids = " OR ".join(f"{p}[uid]" for p in pmids[i:i + chunk])
        out |= set(esearch(f"({query}) AND {LIMITS} AND ({ids})")[1])
    return out


def all_ids(query: str, start: int = 2000, end: int = date.today().year) -> list[str]:
    """Every PMID the query retrieves (per year, per month when a year passes 9,999)."""
    out = []
    for y in range(start, end + 1):
        n, ids = esearch(f'({query}) AND {HAS_ABSTRACT} AND ("{y}/01/01"[dp] : "{y}/12/31"[dp])')
        if n <= 9999:
            out += ids
            continue
        for m in range(1, 13):
            last = calendar.monthrange(y, m)[1]
            out += esearch(f'({query}) AND {HAS_ABSTRACT} AND ("{y}/{m:02d}/01"[dp] : "{y}/{m:02d}/{last}"[dp])')[1]
    return sorted(set(out))


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def fmt(k: int, n: int) -> str:
    p, lo, hi = wilson(k, n)
    return f"{k}/{n} = {p:.0%} [{lo:.0%}-{hi:.0%}]"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- recall
def recall_report(path: Path, refresh: bool) -> dict:
    rows = read_csv(path)
    pos = {r["pmid"] for r in rows if r["safety"] == "1"}
    neg = {r["pmid"] for r in rows if r["safety"] == "0"}
    print(f"\n== recall: {path.name} ({len(rows)} papers, {len(pos)} patient-safety)")
    if refresh:
        found = set()
        for q in AUTHORS[path]:
            found |= set(esearch(f"({q}) AND {LIMITS}")[1])
        if path == HOLD:  # papers co-authored with Bates or Adelman are in the calibration set
            found -= {r["pmid"] for r in read_csv(CAL)}
        new = sorted(found - pos - neg)
        print(f"  author search now returns {len(found)} PMIDs; {len(new)} not labelled"
              + (f" (e.g. {', '.join(new[:10])})" if new else ""))
    hits = {k: members(q, pos | neg) for k, q in QUERIES.items()}
    others = lambda k: set().union(*(hits[j] for j in COMPONENTS_V2 if j != k))  # noqa: E731
    out = {"n": len(rows), "n_safety": len(pos)}
    for k in ("QUERY_V1", "QUERY"):
        out[k] = {"recall": len(hits[k] & pos) / len(pos), "safety_hits": len(hits[k] & pos),
                  "non_safety_hits": len(hits[k] & neg)}
        print(f"  {k:9s} recall {fmt(len(hits[k] & pos), len(pos))}; also retrieves "
              f"{len(hits[k] & neg)}/{len(neg)} of the non-safety papers")
    pre = hits["QUERY"] - (hits[POST_HOC] - others(POST_HOC))
    out["QUERY_without_" + POST_HOC] = len(pre & pos) / len(pos)
    print(f"  QUERY without {POST_HOC}: recall {fmt(len(pre & pos), len(pos))}")
    print(f"  {'component':30s} {'safety hits':>11s} {'unique':>7s} {'non-safety':>10s}")
    out["components"] = {}
    for k in COMPONENTS_V2:
        u = hits[k] - others(k)
        out["components"][k] = {"safety_hits": len(hits[k] & pos), "unique_safety": len(u & pos),
                                "non_safety_hits": len(hits[k] & neg)}
        print(f"  {k:30s} {len(hits[k] & pos):>11d} {len(u & pos):>7d} {len(hits[k] & neg):>10d}")
    missed = sorted(pos - hits["QUERY"])
    titles = {r["pmid"]: r["title"] for r in rows}
    out["missed"] = [{"pmid": p, "title": titles[p]} for p in missed]
    print(f"  missed by QUERY ({len(missed)}):")
    for p in missed:
        print(f"    {p} {titles[p][:100]}")
    return out


# ---------------------------------------------------------------- precision
def precision_report() -> dict:
    rows, meta = read_csv(SAMPLES), json.loads(META.read_text(encoding="utf-8"))
    by = {}
    for r in rows:
        by.setdefault(r["sample"], []).append(r)
    print("\n== precision samples (strict Y; lenient Y+B; 95% Wilson interval)")
    out = {}
    for name, rs in by.items():
        rs = [r for r in rs if r["label"] in ("Y", "B", "N")]
        if not rs:
            print(f"  {name}: no labels yet")
            continue
        target = meta.get(name, {}).get("target")
        cur = rs
        if target in QUERIES:  # keep the papers that the target query retrieves today
            m = members(QUERIES[target], [r["pmid"] for r in rs])
            cur = [r for r in rs if r["pmid"] in m]
        y = sum(r["label"] == "Y" for r in cur)
        yb = sum(r["label"] in ("Y", "B") for r in cur)
        out[name] = {"target": target, "n": len(cur), "Y": y, "YB": yb,
                     "strict": wilson(y, len(cur)), "lenient": wilson(yb, len(cur))}
        tgt = f" -> {target} ({len(cur)}/{len(rs)} still match)" if target in QUERIES else ""
        print(f"  {name:42s} strict {fmt(y, len(cur)):24s} lenient {fmt(yb, len(cur))}{tgt}")
    # Whole-query estimates
    pool = [r for n in ("v2_draft_corpus", "v2_topup") for r in by.get(n, [])]
    m = members(QUERY_TERM, [r["pmid"] for r in pool])
    cur = [r for r in pool if r["pmid"] in m]
    out["QUERY"] = {"n": len(cur), "strict": wilson(sum(r["label"] == "Y" for r in cur), len(cur)),
                    "lenient": wilson(sum(r["label"] in ("Y", "B") for r in cur), len(cur))}
    out["QUERY_V1"] = {k: out["v1_corpus"][k] for k in ("n", "strict", "lenient")} if "v1_corpus" in out else {}
    return out


def sizes() -> dict:
    v1, v2 = f"({QUERY_V1_TERM}) AND {LIMITS}", f"({QUERY_TERM}) AND {LIMITS}"
    s = {"QUERY_V1": count(v1), "QUERY": count(v2), "both": count(f"{v1} AND ({QUERY_TERM})")}
    s["added"], s["removed"] = s["QUERY"] - s["both"], s["QUERY_V1"] - s["both"]
    return s


# ---------------------------------------------------------------- drawing new samples
def draw(name: str, source: str, n: int, seed: str) -> None:
    rows, meta = read_csv(SAMPLES), json.loads(META.read_text(encoding="utf-8"))
    if name in meta:
        raise SystemExit(f"sample {name!r} exists already")
    ids = all_ids(QUERIES[source])
    pick = random.Random(seed).sample(ids, min(n, len(ids)))
    recs, _ = fetch_records(pick)
    by = {r["pmid"]: r for r in recs}
    for p in pick:
        r = by.get(p)
        if r:
            rows.append({"sample": name, "pmid": p, "year": r["year"], "journal": r["journal_abbr"],
                         "title": " ".join(r["title"].split()),
                         "mesh_major": "; ".join(m["d"] for m in r["mesh"] if m["major"]), "label": ""})
    with SAMPLES.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    meta[name] = {"drawn_from": QUERIES[source], "target": source, "n": len(pick), "note": f"seed {seed!r}",
                  "drawn": date.today().isoformat()}
    META.write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(f"drew {len(pick)} of {len(ids):,} {source} papers into {SAMPLES.name}; fill in the label column")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh-authors", action="store_true")
    ap.add_argument("--draw", metavar="NAME")
    ap.add_argument("--source", default="QUERY", choices=list(QUERIES))
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--seed", default="ps_explorer")
    args = ap.parse_args()
    if args.draw:
        return draw(args.draw, args.source, args.n, args.seed)
    summary = {"date": date.today().isoformat(), "sizes": sizes()}
    s = summary["sizes"]
    print(f"== corpus size (abstract, 2000+): v1 {s['QUERY_V1']:,}  v2 {s['QUERY']:,}  "
          f"both {s['both']:,}  added {s['added']:,}  removed {s['removed']:,}")
    summary["calibration"] = recall_report(CAL, args.refresh_authors)
    summary["holdout"] = recall_report(HOLD, args.refresh_authors)
    summary["precision"] = precision_report()
    p1, p2 = summary["precision"].get("QUERY_V1", {}), summary["precision"]["QUERY"]
    print("\n== summary")
    for k, p, n in (("v1", p1, s["QUERY_V1"]), ("v2", p2, s["QUERY"])):
        if p:
            st, le = p["strict"], p["lenient"]
            print(f"  {k}: precision strict {st[0]:.0%} [{st[1]:.0%}-{st[2]:.0%}], lenient {le[0]:.0%} "
                  f"[{le[1]:.0%}-{le[2]:.0%}] (n={p['n']}); ~{st[0] * n:,.0f} "
                  f"[{st[1] * n:,.0f}-{st[2] * n:,.0f}] patient-safety papers of {n:,}")
    SUMMARY.write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(f"wrote {SUMMARY}")


if __name__ == "__main__":
    main()
