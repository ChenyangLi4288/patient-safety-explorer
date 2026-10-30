"""Download the patient-safety + quality-of-care PubMed corpus (has abstract, 2000+).

PubMed esearch returns at most 10k PMIDs per query, so the query is split by
publication month. Each month is written to data/raw/YYYY-MM.jsonl; months that
already exist are skipped, so the script is resumable.

Usage: python fetch_corpus.py [--start 2000] [--end 2026]
"""
import argparse
import calendar
import http.client
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from count_corpus import CORE_Q, HAS_ABSTRACT, QUALITY_Q

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERY = f"({CORE_Q} OR {QUALITY_Q}) AND {HAS_ABSTRACT}"
OUT = Path(__file__).parent / "data" / "raw"
BATCH = 400
MIN_INTERVAL = 0.4  # < 3 req/s without an API key (0.35 drew some 429s)

_last = 0.0


def _request(endpoint: str, params: dict, retries: int = 6) -> bytes:
    global _last
    params = {**params, "tool": "ps_explorer"}
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        wait = MIN_INTERVAL - (time.time() - _last)
        if wait > 0:
            time.sleep(wait)
        _last = time.time()
        try:
            with urllib.request.urlopen(EUTILS + endpoint, data=data, timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
            delay = 2 ** attempt
            print(f"  {endpoint} failed ({e}); retry in {delay}s", flush=True)
            time.sleep(delay)
    raise RuntimeError(f"{endpoint} failed after {retries} attempts")


def search_pmids(term: str, retries: int = 6) -> list[str]:
    # esearch sometimes answers with an error payload (no idlist); retry those.
    for attempt in range(retries):
        raw = _request("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": 10000,
                                        "retmode": "json"})
        try:
            js = json.loads(raw)["esearchresult"]
            ids, total = js["idlist"], int(js["count"])
            break
        except (ValueError, KeyError) as e:
            print(f"  bad esearch response ({e!r}): {raw[:300]!r}; retry", flush=True)
            time.sleep(2 ** attempt)
    else:
        raise RuntimeError("esearch failed repeatedly")
    if total > len(ids):
        raise RuntimeError(f"query returned {total} > {len(ids)} ids; split it further")
    return ids


def _text(el) -> str:
    return "".join(el.itertext()).strip() if el is not None else ""


def parse_article(art) -> dict:
    mc = art.find("MedlineCitation")
    a = mc.find("Article")
    abstract = []
    for t in a.findall("Abstract/AbstractText"):
        label, txt = t.get("Label"), _text(t)
        abstract.append(f"{label}: {txt}" if label else txt)
    year = _text(a.find("Journal/JournalIssue/PubDate/Year")) or \
        _text(a.find("Journal/JournalIssue/PubDate/MedlineDate"))[:4]
    authors = []
    for au in a.findall("AuthorList/Author"):
        name = " ".join(x for x in (_text(au.find("ForeName")), _text(au.find("LastName"))) if x) \
            or _text(au.find("CollectiveName"))
        affs = [_text(x) for x in au.findall("AffiliationInfo/Affiliation")]
        authors.append({"name": name, "aff": affs[0] if affs else ""})
    ids = {i.get("IdType"): _text(i) for i in art.findall("PubmedData/ArticleIdList/ArticleId")}
    return {
        "pmid": _text(mc.find("PMID")),
        "title": _text(a.find("ArticleTitle")),
        "abstract": "\n".join(abstract),
        "journal": _text(a.find("Journal/Title")),
        "journal_abbr": _text(a.find("Journal/ISOAbbreviation")),
        "year": int(year) if year.isdigit() else None,
        "authors": authors,
        "mesh": [{"d": _text(m.find("DescriptorName")),
                  "major": m.find("DescriptorName").get("MajorTopicYN") == "Y"
                  or any(q.get("MajorTopicYN") == "Y" for q in m.findall("QualifierName")),
                  "q": [_text(q) for q in m.findall("QualifierName")]}
                 for m in mc.findall("MeshHeadingList/MeshHeading")],
        "keywords": [_text(k) for k in mc.findall("KeywordList/Keyword")],
        "pub_types": [_text(p) for p in a.findall("PublicationTypeList/PublicationType")],
        "language": _text(a.find("Language")),
        "country": _text(mc.find("MedlineJournalInfo/Country")),
        "doi": ids.get("doi", ""),
        "pmcid": ids.get("pmc", ""),
    }


def parse_book(art) -> dict:
    # NCBI Bookshelf chapters (e.g. AHRQ "Advances in Patient Safety"); no MeSH.
    bd = art.find("BookDocument")
    book = bd.find("Book")
    abstract = []
    for t in bd.findall("Abstract/AbstractText"):
        label, txt = t.get("Label"), _text(t)
        abstract.append(f"{label}: {txt}" if label else txt)
    year = _text(book.find("PubDate/Year"))
    authors = []
    for au in bd.findall("AuthorList[@Type='authors']/Author"):
        name = " ".join(x for x in (_text(au.find("ForeName")), _text(au.find("LastName"))) if x) \
            or _text(au.find("CollectiveName"))
        affs = [_text(x) for x in au.findall("AffiliationInfo/Affiliation")]
        authors.append({"name": name, "aff": affs[0] if affs else ""})
    ids = {i.get("IdType"): _text(i) for i in bd.findall("ArticleIdList/ArticleId")}
    return {
        "pmid": _text(bd.find("PMID")),
        "title": _text(bd.find("ArticleTitle")) or _text(book.find("BookTitle")),
        "abstract": "\n".join(abstract),
        "journal": _text(book.find("BookTitle")),
        "journal_abbr": _text(book.find("CollectionTitle")) or _text(book.find("BookTitle")),
        "year": int(year) if year.isdigit() else None,
        "authors": authors,
        "mesh": [],
        "keywords": [_text(k) for k in bd.findall("KeywordList/Keyword")],
        "pub_types": ["Book Chapter"],
        "language": _text(bd.find("Language")),
        "country": "",
        "doi": ids.get("doi", ""),
        "pmcid": ids.get("bookaccession", ""),
    }


def _efetch(pmids: list[str], retries: int = 5) -> list[dict]:
    # NCBI occasionally returns truncated XML; retry until it parses.
    for attempt in range(retries):
        xml = _request("efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
        try:
            root = ET.fromstring(xml)
        except ET.ParseError as e:
            print(f"  bad XML ({e}); retry", flush=True)
            time.sleep(2 ** attempt)
            continue
        return [parse_article(a) for a in root.findall("PubmedArticle")] + \
            [parse_book(a) for a in root.findall("PubmedBookArticle")]
    raise RuntimeError("efetch returned unparseable XML repeatedly")


def fetch_records(pmids: list[str]) -> tuple[list[dict], list[str]]:
    """Return (records, pmids still missing). Batches can come back silently short,
    so missing ids are refetched in smaller batches."""
    out = []
    for i in range(0, len(pmids), BATCH):
        out.extend(_efetch(pmids[i:i + BATCH]))
    for size in (50, 1):
        missing = sorted(set(pmids) - {r["pmid"] for r in out})
        for i in range(0, len(missing), size):
            out.extend(_efetch(missing[i:i + size]))
    missing = sorted(set(pmids) - {r["pmid"] for r in out})
    if len(missing) > 20:  # single-id refetches are slow; flag rather than hide it
        print(f"  warning: {len(missing)} ids still missing", flush=True)
    return out, missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=2000)
    ap.add_argument("--end", type=int, default=2026)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    for y in range(args.start, args.end + 1):
        for m in range(1, 13):
            path = OUT / f"{y}-{m:02d}.jsonl"
            if path.exists() and path.with_suffix(".missing").exists():
                total += sum(1 for _ in path.open(encoding="utf-8"))
                continue
            last = calendar.monthrange(y, m)[1]
            term = f'{QUERY} AND ("{y}/{m:02d}/01"[dp] : "{y}/{m:02d}/{last}"[dp])'
            pmids = search_pmids(term)
            recs, missing = fetch_records(pmids) if pmids else ([], [])
            path.with_suffix(".missing").write_text("\n".join(missing), encoding="utf-8")
            tmp = path.with_suffix(".tmp")
            with tmp.open("w", encoding="utf-8") as f:
                for r in recs:
                    f.write(json.dumps(r, ensure_ascii=False) + "\n")
            tmp.replace(path)
            total += len(recs)
            print(f"{y}-{m:02d}: {len(pmids):5d} ids, {len(recs):5d} records, "
                  f"{len(missing)} missing  (total {total:,})",
                  flush=True)
    print(f"done: {total:,} records", flush=True)


if __name__ == "__main__":
    sys.exit(main())
