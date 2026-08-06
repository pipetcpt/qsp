#!/usr/bin/env python3
"""Write mr_references.md from the PubMed-resolved metadata.

Titles, journals, years, authors and PMIDs all come from NCBI esummary, so
nothing in the bibliography is written from memory.  A handful of queries
returned a paper other than the one intended; those are re-labelled here with
an honest note rather than silently mis-cited, and a few junk hits are dropped.
"""
import json
from collections import OrderedDict

SECTIONS = OrderedDict([
 ("1", ("The two pivotal trials, and the proportionality argument they created",
        "The model is calibrated on the COAPT control arm alone; every other trial "
        "number in this list is an out-of-sample prediction target.")),
 ("2", ("Quantifying the lesion: EROA, PISA, and the geometry PISA assumes",
        "Denominator 5. The severity threshold everyone uses was calibrated on a "
        "round orifice and is then applied to a crescentic one.")),
 ("3", ("Acute versus chronic, and left atrial compliance",
        "Denominator 1. The same regurgitant volume floods one lung and is silent "
        "in another, and the difference is a property of the atrium.")),
 ("4", ("Ejection fraction, hidden contractility, and the operative thresholds",
        "Denominator 4. EF is computed against a stroke volume that contains the "
        "leak, which is why the surgical trigger sits at 60% and not 50%.")),
 ("5", ("Secondary MR: tethering, annular geometry, and the dynamic orifice",
        "Where the functional orifice comes from, and the negative results that "
        "constrain any mechanism -- annular dilation alone is often not enough.")),
 ("6", ("Leaflet plasticity: the valve is not a passive gasket",
        "The mitral leaflets grow under tethering, but only partially, which is "
        "why the SPEED of ventricular dilation matters and not only its size.")),
 ("7", ("Remodelling, growth laws, fibrosis and contractile loss",
        "The growth laws the model implements, and the reason the hypertrophy of "
        "mitral regurgitation does not protect the ventricle.")),
 ("8", ("Pulmonary vascular and right heart consequences",
        "The second barrier: fix the valve late enough and the lung has become "
        "the disease.")),
 ("9", ("Medical therapy",
        "Effect sizes the drug block is scaled against, and the vasodilator "
        "argument that follows from the leak and the aorta being in parallel.")),
 ("10", ("Devices, surgery, recurrence and guidelines",
         "Which arm of the feedback loop each procedure actually cuts.")),
 ("11", ("Atrial fibrillation and atrial myopathy",
         "The second arm of the vortex, and the reason a long-compensated patient "
         "decompensates with an unchanged valve.")),
 ("12", ("Modelling methodology",
         "Closed-loop lumped-parameter circulation, elastance mechanics, "
         "stress-driven growth, and the simulation engine.")),
])

# queries whose top hit was NOT the intended paper -> re-labelled honestly
RELABEL = {
 "9011684":  "Volume-overload adaptation measured in a conscious large-animal model (aortic regurgitation); the ventricular response to a regurgitant load.",
 "24050987": "Giant left atrial v waves and their pulmonary consequences -- denominator 1 seen directly in a patient.",
 "3053381":  "Echocardiographic detection of mitral valve prolapse and the geometry of the annulus.",
 "33983831": "Inflammation-to-fibrosis signalling in the myocardium; the pathway behind the model's collagen state.",
 "36017548": "Pulmonary hypertension definitions, including the post-capillary and combined categories the model's PVR state maps onto.",
 "32749493": "Contemporary HFrEF pharmacotherapy, including rate control, as a reference for the drug block.",
 "40122615": "Right-sided consequences of degenerative mitral regurgitation in asymptomatic patients.",
 "29044932": "Right heart dysfunction and RV-PA coupling in left-sided heart failure.",
 "19593121": "Three-dimensional echocardiography of the regurgitant orifice, which is how its non-circular shape became visible.",
 "25166465": "Three- versus two-dimensional assessment of functional MR: the measurement disagreement, quantified.",
 "37984204": "Multiscale finite-element ventricular mechanics with baroreflex control; a higher-resolution counterpart to the lumped model used here.",
 "27040451": "Adding an annuloplasty ring to revascularisation in moderate ischaemic MR.",
 "40799133": "Long-term outcomes of early surgery versus watchful waiting in asymptomatic severe MR.",
 "487530":   "The regurgitant orifice area changes DURING ejection -- the dynamic orifice, measured in 1979.",
}
DROP = {"27885969", "37494779"}   # conference-abstract dump; single case report


def load(fn):
    try:
        return json.load(open(fn))
    except Exception:
        return {"rows": [], "meta": {}}


def main():
    rows, meta = [], {}
    for fn in ("refs_raw.json", "refs_raw2.json", "refs_raw3.json"):
        d = load(fn)
        rows += d["rows"]
        meta.update(d["meta"])

    seen, byec = set(), {k: [] for k in SECTIONS}
    for sec, q, note, pmid in rows:
        if not pmid or pmid in seen or pmid in DROP or pmid not in meta:
            continue
        seen.add(pmid)
        byec.setdefault(sec, []).append((pmid, RELABEL.get(pmid, note)))

    total = sum(len(v) for v in byec.values())
    out = []
    out.append("# Mitral Regurgitation — QSP model references\n")
    out.append(
      "**%d references, every one resolved against PubMed programmatically.**\n" % total)
    out.append(
      "The bibliography was not written from memory. Each entry began as a query\n"
      "against the NCBI E-utilities API; the author, title, journal, year and PMID\n"
      "below are the values NCBI returned. Where a query's top hit turned out to be\n"
      "a paper other than the one intended, the entry has been re-labelled to\n"
      "describe what it actually is rather than silently mis-cited, and two junk\n"
      "hits (a conference-abstract compilation and a single case report) were\n"
      "dropped. Regenerate or re-verify with the scripts noted at the foot of this\n"
      "file.\n")
    out.append(
      "The organising idea of the model is that a regurgitant volume means nothing\n"
      "until it is divided by something, and that the disease is really a family of\n"
      "diseases indexed by the denominator. The sections below follow the five\n"
      "denominators.\n")
    out.append("| | denominator | what it decides | section |")
    out.append("|---|---|---|---|")
    out.append("| 1 | operating LA compliance `C_op` | congestion; acute vs chronic | §3 |")
    out.append("| 2 | total stroke volume | regurgitant fraction, the grade | §2 |")
    out.append("| 3 | LV end-diastolic volume | valve or ventricle; COAPT vs MITRA-FR | §1 |")
    out.append("| 4 | the afterload the leak removes | hidden contractility, the EF 60% rule | §4 |")
    out.append("| 5 | `k_PISA` | the measurement's own bias | §2 |")
    out.append("")
    out.append("---\n")

    for sec, (title, blurb) in SECTIONS.items():
        items = byec.get(sec, [])
        out.append("## %s. %s\n" % (sec, title))
        out.append("*%s*\n" % blurb)
        for pmid, note in items:
            m = meta[pmid]
            au = m.get("authors", [])
            names = [a["name"] for a in au if a.get("authtype") == "Author"]
            if len(names) > 3:
                alist = ", ".join(names[:3]) + ", et al."
            else:
                alist = ", ".join(names) if names else "—"
            yr = (m.get("pubdate") or "")[:4]
            ttl = (m.get("title") or "").rstrip(".")
            src = m.get("source") or ""
            alist = alist.rstrip(".")
            out.append("- %s. **%s.** *%s* %s. "
                       "[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)  \n  %s"
                       % (alist, ttl, src, yr, pmid, pmid, note))
        out.append("")

    out.append("---\n")
    out.append("## Reproducing this bibliography\n")
    out.append("Each entry was resolved with `esearch` (relevance-sorted, top hit) and\n"
               "described with `esummary` against `eutils.ncbi.nlm.nih.gov`. The harvest\n"
               "and rendering scripts are `refs.py`, `refs2.py`, `refs3.py` and\n"
               "`mkrefs.py`; the raw NCBI payloads are kept in `refs_raw*.json`. Any PMID\n"
               "in this file can be checked directly:\n")
    out.append("```bash\ncurl -s \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi\"\\\n"
               "     \"?db=pubmed&retmode=json&id=30280640\" | python3 -m json.tool\n```\n")
    out.append("## A note on what is calibrated and what is predicted\n")
    out.append("Only two numbers in the whole model are fitted to outcome data, and both\n"
               "come from the **COAPT control arm**: the baseline hazard for heart-failure\n"
               "hospitalisation and the baseline hazard for death. Every other trial\n"
               "quantity referenced above — the COAPT device arm, both MITRA-FR arms, and\n"
               "the baseline echocardiographic and haemodynamic values the virtual patients\n"
               "reproduce — is a prediction, and the ones the model gets wrong are reported\n"
               "as misses in `README.md` rather than absorbed by refitting.\n")

    txt = "\n".join(out) + "\n"
    open("/home/user/qsp/mitral-regurgitation/mr_references.md", "w").write(txt)
    print("wrote %d references" % total)
    for sec in SECTIONS:
        print("  §%-3s %2d" % (sec, len(byec.get(sec, []))))


if __name__ == "__main__":
    main()
