"""Count PubMed hits for candidate patient-safety corpus queries (no download)."""
import json
import time
import urllib.parse
import urllib.request

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

COMPONENTS = {
    "mesh_medical_errors": '"Medical Errors"[MeSH]',
    "mesh_patient_safety": '"Patient Safety"[MeSH]',
    "mesh_adr": '"Drug-Related Side Effects and Adverse Reactions"[MeSH]',
    "mesh_iatrogenic": '"Iatrogenic Disease"[MeSH]',
    "mesh_patient_harm": '"Patient Harm"[MeSH]',
    "mesh_safety_mgmt": '"Safety Management"[MeSH]',
    "mesh_risk_mgmt": '"Risk Management"[MeSH]',
    "mesh_quality_of_care": '"Quality of Health Care"[MeSH]',
    "mesh_quality_improvement": '"Quality Improvement"[MeSH]',
    "tiab_adverse_event": '"adverse event*"[tiab] AND (hospital*[tiab] OR inpatient*[tiab] OR "health care"[tiab] OR healthcare[tiab])',
    "tiab_medical_error": '("medical error*"[tiab] OR "medication error*"[tiab] OR "diagnostic error*"[tiab] OR "surgical error*"[tiab])',
    "tiab_patient_safety": '"patient safety"[tiab]',
    "tiab_never_event": '("never event*"[tiab] OR "sentinel event*"[tiab] OR "near miss*"[tiab] OR "incident report*"[tiab] OR "trigger tool"[tiab])',
    "tiab_quality_of_care": '("quality of care"[tiab] OR "quality improvement"[tiab] OR "healthcare quality"[tiab] OR "health care quality"[tiab])',
}

CORE = ["mesh_medical_errors", "mesh_patient_safety", "mesh_iatrogenic", "mesh_patient_harm",
        "tiab_medical_error", "tiab_patient_safety", "tiab_never_event", "tiab_adverse_event"]
QUALITY = ["mesh_quality_improvement", "tiab_quality_of_care"]
CORE_Q = "(" + " OR ".join(COMPONENTS[k] for k in CORE) + ")"
QUALITY_Q = "(" + " OR ".join(COMPONENTS[k] for k in QUALITY) + ")"
HAS_ABSTRACT = "hasabstract[text]"
SINCE_2000 = '("2000/01/01"[dp] : "3000"[dp])'

COMBOS = {
    "CORE (safety/error/AE)": CORE_Q,
    "CORE + quality (narrow QI)": f"{CORE_Q} OR {QUALITY_Q}",
    "CORE + quality, abstract, 2000+": f"({CORE_Q} OR {QUALITY_Q}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "CORE, abstract, 2000+": f"{CORE_Q} AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "CORE, abstract, 2015+": f'{CORE_Q} AND {HAS_ABSTRACT} AND ("2015/01/01"[dp] : "3000"[dp])',
    "CORE x ML/AI, abstract": f'{CORE_Q} AND {HAS_ABSTRACT} AND ("machine learning"[tiab] OR "deep learning"[tiab] OR "artificial intelligence"[tiab] OR "natural language processing"[tiab] OR "large language model*"[tiab] OR "Machine Learning"[MeSH])',
}


def count(term: str) -> int:
    url = ESEARCH + "?" + urllib.parse.urlencode(
        {"db": "pubmed", "term": term, "rettype": "count", "retmode": "json", "tool": "ps_explorer"})
    with urllib.request.urlopen(url, timeout=60) as r:
        n = int(json.load(r)["esearchresult"]["count"])
    time.sleep(0.4)  # stay under 3 req/s without an API key
    return n


def main():
    print("== components ==")
    for k, q in COMPONENTS.items():
        print(f"{k:28s} {count(q):>10,d}")
    print("\n== combined ==")
    for k, q in COMBOS.items():
        print(f"{k:34s} {count(q):>10,d}")
    print("\n== CORE, abstract, by 5-year bin ==")
    for y0 in range(1990, 2030, 5):
        q = f'{CORE_Q} AND {HAS_ABSTRACT} AND ("{y0}/01/01"[dp] : "{y0 + 4}/12/31"[dp])'
        print(f"{y0}-{y0 + 4}  {count(q):>10,d}")


if __name__ == "__main__":
    main()
