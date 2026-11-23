"""PubMed query definitions for the patient-safety corpus, plus hit counts (no download).

Two versions live here:

* v1 (QUERY_V1_TERM): the query behind the first map (255,023 papers, October 2026). Kept
  verbatim for reproducibility. Note that the `tiab_adverse_event` component has no
  parentheses of its own, and PubMed evaluates AND/OR left to right. So v1 actually ran as
  `(all other CORE terms OR "adverse event*") AND (hospital* OR inpatient* OR "health care"
  OR healthcare)`, OR the broad quality-of-care terms.
* v2 (QUERY_TERM): the current query, built and checked with validate_query.py. Each
  component is wrapped in its own parentheses.

Usage: python count_corpus.py [--v1]
"""
import argparse
import json
import time
import urllib.parse
import urllib.request

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
HAS_ABSTRACT = "hasabstract[text]"
SINCE_2000 = '("2000/01/01"[dp] : "3000"[dp])'


def _or(parts) -> str:
    return "(" + " OR ".join(f"({p})" for p in parts) + ")"


# ---------------------------------------------------------------- v1 (first map)
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
# Joined without per-component parentheses, exactly as v1 ran (see the module docstring).
CORE_Q_V1 = "(" + " OR ".join(COMPONENTS[k] for k in CORE) + ")"
QUALITY_Q_V1 = "(" + " OR ".join(COMPONENTS[k] for k in QUALITY) + ")"
QUERY_V1_TERM = f"{CORE_Q_V1} OR {QUALITY_Q_V1}"

# ---------------------------------------------------------------- v2 (current)
# Building blocks
CARE = ('(hospital*[tiab] OR inpatient*[tiab] OR hospitali*[tiab] OR "nursing home*"[tiab] OR '
        '"long-term care"[tiab] OR "intensive care"[tiab] OR ICU[tiab] OR ward*[tiab] OR '
        '"acute care"[tiab] OR "health care"[tiab] OR healthcare[tiab])')
HEALTHCARE_CTX = ('(patient*[tiab] OR hospital*[tiab] OR healthcare[tiab] OR "health care"[tiab] OR '
                  'clinical[tiab] OR nurs*[tiab] OR medical[tiab])')
# Bench microbiology, molecular typing and veterinary papers that share HAI MeSH.
NOT_MICRO = ('("Bacterial Typing Techniques"[MeSH] OR "Genome, Bacterial"[MeSH] OR '
             '"Whole Genome Sequencing"[MeSH] OR "Molecular Epidemiology"[MeSH] OR '
             '"Microbial Sensitivity Tests"[MeSH] OR "Bacterial Proteins"[MeSH] OR "Virulence"[MeSH] OR '
             '"Genes, Bacterial"[MeSH] OR "Biofilms"[Majr] OR "Biosensing Techniques"[MeSH] OR '
             '"Anti-Bacterial Agents/pharmacology"[MeSH] OR "Hospitals, Animal"[MeSH] OR '
             '"Veterinarians"[MeSH] OR (Animals[MeSH] NOT Humans[MeSH]))')
HAI_TIAB = ('("healthcare-associated infection*"[tiab] OR "health care-associated infection*"[tiab] OR '
            '"hospital-acquired infection*"[tiab] OR "nosocomial infection*"[tiab] OR '
            '"central line-associated"[tiab] OR CLABSI[tiab] OR CAUTI[tiab] OR '
            '"catheter-associated urinary"[tiab])')
HIT = ('("Medical Order Entry Systems"[MeSH] OR "Decision Support Systems, Clinical"[MeSH] OR '
       '"Electronic Prescribing"[MeSH] OR "Drug Therapy, Computer-Assisted"[MeSH] OR '
       '"Medication Systems, Hospital"[MeSH] OR "Clinical Alarms"[MeSH] OR "Infusion Pumps"[MeSH] OR '
       '"order entry"[tiab] OR CPOE[tiab] OR "decision support"[tiab] OR "electronic prescribing"[tiab] OR '
       'e-prescribing[tiab] OR "electronic health record*"[tiab] OR "health information technology"[tiab] OR '
       '"smart pump*"[tiab] OR "infusion pump*"[tiab] OR barcode*[tiab] OR "bar code*"[tiab] OR '
       '"bar-code*"[tiab])')
HIT_SAFETY = ('("medication error*"[tiab] OR "prescribing error*"[tiab] OR error*[ti] OR '
              '"patient safety"[tiab] OR "medication safety"[tiab] OR "adverse drug"[tiab] OR harm[tiab] OR '
              'alert*[tiab] OR override*[tiab] OR "Medication Errors"[MeSH] OR "Drug Interactions"[MeSH] OR '
              '"Drug Hypersensitivity"[MeSH] OR "Drug-Related Side Effects and Adverse Reactions"[MeSH])')
QI_FOCUS = ('("Quality Improvement"[Majr] OR "quality improvement"[ti] OR "quality of care"[ti] OR '
            '"healthcare quality"[ti] OR "health care quality"[ti] OR "improvement project"[ti] OR '
            '"improvement initiative"[ti])')
SAFETY_TERMS = ('("patient safety"[tiab] OR "safety culture"[tiab] OR "medical error*"[tiab] OR '
                '"medication error*"[tiab] OR error*[ti] OR "adverse event*"[tiab] OR "adverse drug"[tiab] OR '
                'iatrogenic[tiab] OR "preventable adverse"[tiab] OR "preventable harm"[tiab] OR '
                '"patient harm"[tiab] OR "avoidable harm"[tiab] OR "hospital-acquired"[tiab] OR '
                '"healthcare-associated"[tiab] OR "health care-associated"[tiab] OR "near miss*"[tiab] OR '
                '"never event*"[tiab] OR "Patient Safety"[MeSH] OR "Medical Errors"[MeSH:noexp] OR '
                '"Medication Errors"[MeSH:noexp] OR "Iatrogenic Disease"[MeSH])')

COMPONENTS_V2 = {
    # -- MeSH
    # "Medical Errors"[MeSH] explodes into Diagnostic Errors -> False Positive/Negative Reactions
    # and Observer Variation (test accuracy, inter-rater reliability): ~48% of its 2000+ hits.
    # So the error descriptors are listed individually here.
    "mesh_errors": ('"Medical Errors"[MeSH:noexp] OR "Medication Errors"[MeSH:noexp] OR '
                    '"Medication Reconciliation"[MeSH] OR "Near Miss, Healthcare"[MeSH] OR '
                    '("Inappropriate Prescribing"[MeSH] AND ("Potentially Inappropriate Medication List"[MeSH] OR '
                    '"potentially inappropriate"[tiab] OR older[tiab] OR elderly[tiab] OR adverse[tiab] OR '
                    'harm*[tiab] OR "Drug-Related Side Effects and Adverse Reactions"[MeSH]))'),
    "mesh_diagnostic_errors": ('("Diagnostic Errors"[MeSH:noexp] OR "Missed Diagnosis"[MeSH]) AND '
                               '("Patient Safety"[MeSH] OR "Medical Errors"[MeSH:noexp] OR "patient safety"[tiab] OR '
                               '"diagnostic error*"[tiab] OR "diagnostic safety"[tiab] OR malpractice[tiab] OR '
                               '"missed opportunit*"[tiab])'),
    "mesh_patient_safety": '"Patient Safety"[MeSH] OR "Patient Harm"[MeSH]',
    "mesh_iatrogenic": '"Iatrogenic Disease"[Majr:noexp] NOT "Case Reports"[pt]',
    "mesh_safety_systems": (f'("Safety Management"[MeSH] AND (patient*[tiab] OR {CARE})) OR '
                            '"Root Cause Analysis"[MeSH] OR ("Truth Disclosure"[MeSH] AND (error*[tiab] OR '
                            '"adverse event*"[tiab] OR harm[tiab] OR mistake*[tiab] OR "Medical Errors"[MeSH]))'),
    "mesh_handoff_deterioration": ('"Patient Handoff"[MeSH] OR "Clinical Deterioration"[MeSH] OR '
                                   '"Hospital Rapid Response Team"[MeSH]'),
    "pressure_injury_falls": (f'("Pressure Ulcer"[MeSH] AND {CARE}) OR ("Accidental Falls"[MeSH] AND '
                              '(inpatient*[tiab] OR "in-hospital"[tiab] OR hospitalized[tiab] OR hospitalised[tiab] OR '
                              'ward*[tiab] OR "acute care"[tiab] OR "patient safety"[tiab])) OR "inpatient fall*"[tiab] OR '
                              '"in-hospital fall*"[tiab] OR "hospital fall*"[tiab] OR "patient fall*"[tiab] OR '
                              '"hospital-acquired pressure"[tiab] OR "unplanned extubation*"[tiab] OR '
                              '("fall prevention"[tiab] AND (hospital*[tiab] OR inpatient*[tiab] OR nurs*[tiab]))'),
    "hai_prevention": (f'(("Cross Infection/prevention and control"[MeSH] OR ("Hand Hygiene"[MeSH] AND {CARE}) OR '
                       f'({HAI_TIAB} AND (prevent*[tiab] OR reduc*[tiab] OR bundle*[tiab] OR compliance[tiab] OR '
                       '"quality improvement"[tiab] OR "infection control"[tiab] OR "infection prevention"[tiab] OR '
                       'surveillance[tiab])) OR (("Cross Infection"[MeSH] OR "Catheter-Related Infections"[MeSH]) AND '
                       '("Infection Control"[MeSH] OR "Equipment Contamination"[MeSH] OR "Disinfection"[MeSH] OR '
                       f'"Hand Disinfection"[MeSH]))) NOT {NOT_MICRO}'),
    "hit_safety": (f'({HIT} AND {HIT_SAFETY}) OR "EHR safety"[tiab] OR "health IT safety"[tiab] OR '
                   '"health IT-related"[tiab] OR "HIT-related"[tiab] OR "e-iatrogenesis"[tiab] OR "SAFER guide*"[tiab]'),
    # -- text
    "tiab_medical_error": COMPONENTS["tiab_medical_error"],
    "tiab_patient_safety": ('"patient safety"[tiab] AND ("patient safety"[ti] OR error*[tiab] OR harm[tiab] OR '
                            '"adverse event*"[tiab] OR incident*[tiab] OR culture[tiab] OR climate[tiab] OR '
                            '"near miss*"[tiab] OR preventable[tiab] OR iatrogenic[tiab])'),
    "tiab_incidents": ('"never event*"[tiab] OR "sentinel event*"[tiab] OR "incident report*"[tiab] OR '
                       '"trigger tool*"[tiab] OR (("near miss*"[tiab] AND (error*[tiab] OR medication*[tiab] OR '
                       'incident*[tiab] OR report*[tiab] OR "patient safety"[tiab])) NOT "maternal near miss*"[tiab])'),
    "tiab_ade": ('"adverse drug event*"[tiab] OR "preventable adverse"[tiab] OR "medication safety"[tiab] OR '
                 '"medication-related harm"[tiab] OR ("drug-related problem*"[tiab] AND (pharmac*[tiab] OR '
                 'medication*[tiab])) OR "medication-related problem*"[tiab] OR "medication-related hospital*"[tiab] OR '
                 '"drug-related hospital*"[tiab] OR "drug-related admission*"[tiab] OR "medication-related admission*"[tiab]'),
    "tiab_error_types": ('"prescribing error*"[tiab] OR "dispensing error*"[tiab] OR "administration error*"[tiab] OR '
                         '"dosing error*"[tiab] OR "drug error*"[tiab] OR "transfusion error*"[tiab] OR '
                         '"preanalytical error*"[tiab] OR "pre-analytical error*"[tiab] OR "wrong-patient"[tiab] OR '
                         '"wrong patient"[tiab] OR "wrong-site"[tiab] OR "wrong site"[tiab] OR "wrong-side"[tiab] OR '
                         '"retract-and-reorder"[tiab] OR "sound-alike"[tiab] OR "high-alert medication*"[tiab] OR '
                         '"retained surgical"[tiab]'),
    "tiab_alerts": ('"alert fatigue"[tiab] OR "alarm fatigue"[tiab] OR "alert override*"[tiab] OR "override rate*"[tiab] OR '
                    '"overridden alert*"[tiab] OR "interruptive alert*"[tiab] OR "medication alert*"[tiab] OR '
                    '"drug alert*"[tiab] OR "drug-drug interaction alert*"[tiab] OR "interaction alert*"[tiab] OR '
                    '"allergy alert*"[tiab] OR "alarm safety"[tiab]'),
    "tiab_culture": ('"second victim*"[tiab] OR "just culture"[tiab] OR (("safety culture"[tiab] OR '
                     '"safety climate"[tiab] OR "culture of safety"[tiab] OR "root cause analys*"[tiab] OR '
                     '"failure mode and effect*"[tiab] OR "failure modes and effect*"[tiab] OR '
                     '"high reliability organi*"[tiab] OR "safety-II"[tiab] OR "speaking up"[tiab] OR '
                     f'"safety huddle*"[tiab] OR "safety checklist*"[tiab]) AND {HEALTHCARE_CTX})'),
    "tiab_measurement": ('"patient safety indicator*"[tiab] OR "global trigger tool"[tiab] OR "preventable harm"[tiab] OR '
                         '"avoidable harm"[tiab] OR "patient harm"[tiab] OR "iatrogenic harm"[tiab] OR '
                         '"safety incident*"[tiab] OR "adverse incident*"[tiab] OR "unsafe care"[tiab] OR '
                         '"diagnostic safety"[tiab] OR ("safety event*"[tiab] AND (patient*[tiab] OR hospital*[tiab])) OR '
                         '"critical test result*"[tiab] OR "critical result*"[tiab] OR "test results pending"[tiab] OR '
                         '"pending test result*"[tiab] OR "abnormal test result*"[tiab] OR "test result follow-up"[tiab] OR '
                         '"test result management"[tiab]'),
    "tiab_diagnostic_safety": ('("missed opportunit*"[tiab] AND ("diagnostic error*"[tiab] OR "missed diagnos*"[tiab] OR '
                               '"delayed diagnos*"[tiab] OR "diagnostic delay*"[tiab] OR "timely diagnos*"[tiab] OR '
                               '"earlier diagnos*"[tiab] OR "diagnostic evaluation"[tiab] OR "diagnostic opportunit*"[tiab])) OR '
                               '"missed diagnostic opportunit*"[tiab] OR "delays in follow-up"[tiab] OR '
                               '"delayed follow-up of"[tiab] OR "follow-up of abnormal"[tiab] OR '
                               '"lack of timely follow-up"[tiab] OR "missed test result*"[tiab] OR '
                               '"test results management"[tiab] OR (("closing the loop"[tiab] OR "loop closure"[tiab]) AND '
                               '("test result*"[tiab] OR referral*[tiab])) OR (trigger*[tiab] AND ("diagnostic evaluation"[tiab] OR '
                               '"missed opportunit*"[tiab] OR "delays in diagnos*"[tiab] OR "diagnostic error*"[tiab]))'),
    "tiab_rapid_response": ('"rapid response team*"[tiab] OR "rapid response system*"[tiab] OR '
                            '"medical emergency team*"[tiab] OR "failure to rescue"[tiab]'),
    "tiab_handoff": ('(handoff*[tiab] OR hand-off*[tiab] OR handover*[tiab]) AND (communication[tiab] OR '
                     'safety[tiab] OR error*[tiab] OR shift*[tiab])'),
    "malpractice_claims": ('"closed claim*"[tiab] OR ("Malpractice"[MeSH] AND ("Medical Errors"[MeSH] OR '
                           '"Patient Safety"[MeSH] OR error*[tiab] OR "adverse event*"[tiab] OR "patient safety"[tiab]))'),
    # -- quality arm: QI as the main subject AND an explicit safety / harm term
    "quality_improvement_safety": f'{QI_FOCUS} AND {SAFETY_TERMS}',
}
CORE_V2 = [k for k in COMPONENTS_V2 if k != "quality_improvement_safety"]
QUALITY_V2 = ["quality_improvement_safety"]
CORE_Q = _or(COMPONENTS_V2[k] for k in CORE_V2)
QUALITY_Q = _or(COMPONENTS_V2[k] for k in QUALITY_V2)
QUERY_TERM = f"{CORE_Q} OR {QUALITY_Q}"

# v1 components that are broad on their own, kept for comparison
V1_REFERENCE = {
    "v1 mesh_medical_errors (explodes)": COMPONENTS["mesh_medical_errors"],
    "v1 tiab_adverse_event": f'({COMPONENTS["tiab_adverse_event"]})',
    "v1 mesh_quality_improvement": COMPONENTS["mesh_quality_improvement"],
    "v1 tiab_quality_of_care": COMPONENTS["tiab_quality_of_care"],
}
COMBOS = {
    "v1 query as run, abstract, 2000+": f"({QUERY_V1_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 CORE, abstract, 2000+": f"{CORE_Q} AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 quality arm only (not in CORE)": f"({QUALITY_Q} NOT {CORE_Q}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 query, abstract, 2000+": f"({QUERY_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 AND v1": f"({QUERY_TERM}) AND ({QUERY_V1_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 NOT v1 (added)": f"({QUERY_TERM}) NOT ({QUERY_V1_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v1 NOT v2 (removed)": f"({QUERY_V1_TERM}) NOT ({QUERY_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000}",
    "v2 x ML/AI, abstract": f'({QUERY_TERM}) AND {HAS_ABSTRACT} AND {SINCE_2000} AND ("machine learning"[tiab] OR "deep learning"[tiab] OR "artificial intelligence"[tiab] OR "natural language processing"[tiab] OR "large language model*"[tiab] OR "Machine Learning"[MeSH])',
}


def count(term: str) -> int:
    data = urllib.parse.urlencode({"db": "pubmed", "term": term, "rettype": "count", "retmode": "json",
                                   "tool": "ps_explorer"}).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(ESEARCH, data=data, timeout=120) as r:  # POST: long terms
                n = int(json.load(r)["esearchresult"]["count"])
            time.sleep(0.4)  # stay under 3 req/s without an API key
            return n
        except (OSError, ValueError, KeyError):
            time.sleep(2 ** attempt)
    raise RuntimeError("esearch failed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1", action="store_true", help="also count the v1 components")
    args = ap.parse_args()
    lim = f" AND {HAS_ABSTRACT} AND {SINCE_2000}"
    rest = lambda k: _or(COMPONENTS_V2[j] for j in COMPONENTS_V2 if j != k)  # noqa: E731
    print(f"== v2 components (abstract, 2000+)        {'hits':>8s} {'unique':>8s}")
    for k, q in COMPONENTS_V2.items():
        print(f"{k:40s} {count(f'({q}){lim}'):>8,d} {count(f'(({q}) NOT {rest(k)}){lim}'):>8,d}")
    if args.v1:
        print("\n== v1 components (abstract, 2000+)")
        for k, q in COMPONENTS.items():
            print(f"{k:40s} {count(f'({q}){lim}'):>8,d}")
    print("\n== v1 broad pieces (abstract, 2000+)")
    for k, q in V1_REFERENCE.items():
        print(f"{k:40s} {count(f'({q}){lim}'):>8,d}")
    print("\n== combined")
    for k, q in COMBOS.items():
        print(f"{k:40s} {count(q):>8,d}")
    print("\n== v2 query, abstract, by 5-year bin")
    for y0 in range(2000, 2030, 5):
        q = f'({QUERY_TERM}) AND {HAS_ABSTRACT} AND ("{y0}/01/01"[dp] : "{y0 + 4}/12/31"[dp])'
        print(f"{y0}-{y0 + 4}  {count(q):>10,d}")


if __name__ == "__main__":
    main()
