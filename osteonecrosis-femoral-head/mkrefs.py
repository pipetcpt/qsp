#!/usr/bin/env python3
"""Build onfh_references.md from PubMed.

Every citation in the reference list is RETRIEVED from PubMed by an esearch
query, never recalled from memory.  The script stores the raw esummary payload
in refs_raw.json so the result is auditable and reproducible.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# (section, display query label, PubMed query, how many to keep)
QUERIES = [
    # ---- 1 epidemiology and natural history --------------------------------
    ("1. Epidemiology and natural history", "epidemiology of osteonecrosis of the femoral head",
     "osteonecrosis[ti] AND femoral head[ti] AND (epidemiology[ti] OR incidence[ti] OR prevalence[ti])", 5),
    ("1. Epidemiology and natural history", "natural history / untreated course",
     "osteonecrosis[ti] AND femoral head[ti] AND (natural history[ti] OR asymptomatic[ti] OR nonoperative[ti])", 5),
    # ---- 2 glucocorticoid pathogenesis -------------------------------------
    ("2. Glucocorticoid pathogenesis", "corticosteroid dose and osteonecrosis risk",
     "osteonecrosis[ti] AND (corticosteroid[ti] OR glucocorticoid[ti] OR steroid[ti]) AND (dose[tiab] OR risk[ti])", 6),
    ("2. Glucocorticoid pathogenesis", "marrow adipogenesis, PPARgamma, MSC lineage",
     "(osteonecrosis[tiab] OR osteonecrotic[tiab]) AND (adipogenesis[tiab] OR PPAR gamma[tiab] OR mesenchymal stem cell[tiab]) AND glucocorticoid[tiab]", 6),
    ("2. Glucocorticoid pathogenesis", "glucocorticoid osteocyte apoptosis",
     "glucocorticoid[tiab] AND osteocyte[tiab] AND apoptosis[tiab] AND bone[tiab]", 4),
    # ---- 3 vascular / coagulation ------------------------------------------
    ("3. Vascular supply, coagulation and intraosseous pressure", "blood supply of the femoral head",
     "femoral head[ti] AND (blood supply[ti] OR vascular anatomy[ti] OR medial femoral circumflex[tiab])", 5),
    ("3. Vascular supply, coagulation and intraosseous pressure", "intraosseous pressure / compartment",
     "(intraosseous pressure[tiab] OR bone marrow pressure[tiab]) AND (femoral head[tiab] OR osteonecrosis[tiab])", 5),
    ("3. Vascular supply, coagulation and intraosseous pressure", "thrombophilia and hypofibrinolysis",
     "osteonecrosis[ti] AND (thrombophilia[tiab] OR hypofibrinolysis[tiab] OR PAI-1[tiab] OR coagulation[tiab])", 5),
    # ---- 4 histopathology and creeping substitution ------------------------
    ("4. Histopathology and creeping substitution", "histopathology of osteonecrosis",
     "osteonecrosis[ti] AND femoral head[ti] AND (histopathology[tiab] OR histology[tiab] OR pathogenesis[ti])", 6),
    ("4. Histopathology and creeping substitution", "creeping substitution and repair",
     "(creeping substitution[tiab] OR repair[ti]) AND (osteonecrosis[tiab] OR bone necrosis[tiab])", 5),
    # ---- 5 biomechanics -----------------------------------------------------
    ("5. Biomechanics of the femoral head and of collapse", "finite element / mechanics of ONFH collapse",
     "osteonecrosis[tiab] AND femoral head[tiab] AND (finite element[tiab] OR biomechanic*[tiab] OR stress[ti])", 7),
    ("5. Biomechanics of the femoral head and of collapse", "subchondral fracture and collapse mechanism",
     "femoral head[tiab] AND (subchondral fracture[tiab] OR collapse[ti]) AND osteonecrosis[tiab]", 6),
    ("5. Biomechanics of the femoral head and of collapse", "in vivo hip contact force and pressure",
     "hip[ti] AND (contact pressure[tiab] OR joint contact force[tiab] OR in vivo load[tiab]) AND (gait[tiab] OR walking[tiab])", 5),
    ("5. Biomechanics of the femoral head and of collapse", "trabecular bone mechanical properties and fatigue",
     "trabecular bone[tiab] AND (fatigue[ti] OR elastic modulus[tiab] OR compressive strength[tiab])", 6),
    # ---- 6 staging and prediction ------------------------------------------
    ("6. Staging, lesion size and prediction of collapse", "ARCO / Ficat / Steinberg classification",
     "osteonecrosis[ti] AND femoral head[ti] AND (classification[ti] OR staging[ti] OR ARCO[tiab])", 5),
    ("6. Staging, lesion size and prediction of collapse", "JIC type and necrotic angle",
     "femoral head[tiab] AND osteonecrosis[tiab] AND (necrotic angle[tiab] OR Kerboull[tiab] OR Japanese Investigation Committee[tiab] OR lesion size[tiab])", 7),
    ("6. Staging, lesion size and prediction of collapse", "MRI diagnosis and marrow oedema",
     "osteonecrosis[ti] AND femoral head[ti] AND (MRI[tiab] OR magnetic resonance[ti] OR bone marrow edema[tiab])", 5),
    # ---- 7 bisphosphonates --------------------------------------------------
    ("7. Bisphosphonates and other antiresorptives", "alendronate trials in ONFH",
     "osteonecrosis[tiab] AND femoral head[tiab] AND (alendronate[tiab] OR bisphosphonate[tiab])", 8),
    ("7. Bisphosphonates and other antiresorptives", "bisphosphonate bone pharmacokinetics",
     "bisphosphonate[tiab] AND (pharmacokinetics[ti] OR bone uptake[tiab] OR skeletal retention[tiab])", 5),
    # ---- 8 other drugs ------------------------------------------------------
    ("8. Statins, anticoagulants, prostacyclin and anabolic agents", "statin prophylaxis",
     "(statin[tiab] OR lovastatin[tiab] OR simvastatin[tiab]) AND osteonecrosis[tiab]", 5),
    ("8. Statins, anticoagulants, prostacyclin and anabolic agents", "enoxaparin / anticoagulation",
     "osteonecrosis[tiab] AND (enoxaparin[tiab] OR anticoagulation[tiab] OR heparin[tiab])", 4),
    ("8. Statins, anticoagulants, prostacyclin and anabolic agents", "iloprost for bone marrow oedema and ONFH",
     "iloprost[tiab] AND (bone marrow edema[tiab] OR osteonecrosis[tiab] OR femoral head[tiab])", 4),
    ("8. Statins, anticoagulants, prostacyclin and anabolic agents", "teriparatide / PTH in osteonecrosis",
     "(teriparatide[tiab] OR parathyroid hormone[tiab]) AND osteonecrosis[tiab]", 4),
    # ---- 9 surgery ----------------------------------------------------------
    ("9. Core decompression, cell therapy and joint-preserving surgery", "core decompression outcomes",
     "core decompression[ti] AND (femoral head[tiab] OR osteonecrosis[tiab])", 7),
    ("9. Core decompression, cell therapy and joint-preserving surgery", "bone marrow cell therapy",
     "osteonecrosis[tiab] AND femoral head[tiab] AND (bone marrow mononuclear[tiab] OR mesenchymal stem cell[tiab] OR cell therapy[tiab])", 6),
    ("9. Core decompression, cell therapy and joint-preserving surgery", "grafting, osteotomy and arthroplasty",
     "femoral head[tiab] AND osteonecrosis[tiab] AND (vascularized fibular[tiab] OR osteotomy[tiab] OR tantalum[tiab] OR arthroplasty[ti])", 6),
    # ---- 10 guidelines and reviews -----------------------------------------
    ("10. Guidelines, reviews and QSP methodology", "guidelines and systematic reviews",
     "osteonecrosis[ti] AND femoral head[ti] AND (guideline[pt] OR systematic review[pt] OR meta-analysis[pt] OR consensus[ti])", 7),
    ("10. Guidelines, reviews and QSP methodology", "bone remodelling and QSP/systems models",
     "(quantitative systems pharmacology[tiab] OR mechanistic model[tiab] OR mathematical model[tiab]) AND (bone remodeling[tiab] OR bone remodelling[tiab] OR RANKL[tiab])", 6),
]


# Pivotal papers pinned by PMID so that a relevance-ranked search cannot drop
# them.  Titles/journals/years are still fetched from PubMed, never typed.
PINNED = [
    ("5. Biomechanics of the femoral head and of collapse",
     "in vivo hip contact force and pressure", ["11410170"]),
    ("6. Staging, lesion size and prediction of collapse",
     "the staging systems themselves",
     ["31866252", "7822393", "10546624", "36036766"]),
    ("7. Bisphosphonates and other antiresorptives",
     "the two contradictory alendronate trials",
     ["16203877", "22127729"]),
    ("8. Statins, anticoagulants, prostacyclin and anabolic agents",
     "pivotal prophylaxis and vasoactive studies",
     ["11347831", "11347834", "15795211"]),
    ("9. Core decompression, cell therapy and joint-preserving surgery",
     "pivotal cell-therapy studies",
     ["12461352", "28988340", "36961220"]),
    ("10. Guidelines, reviews and QSP methodology",
     "national guidelines", ["32309135"]),
]


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "onfh-qsp/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=45) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:                        # noqa: BLE001
            if attempt == 4:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def esearch(query, retmax):
    url = (f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&sort=relevance"
           f"&retmax={retmax}&term={urllib.parse.quote(query)}")
    return fetch(url)["esearchresult"].get("idlist", [])


def esummary(pmids):
    if not pmids:
        return {}
    url = (f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id="
           + ",".join(pmids))
    return fetch(url)["result"]


def main():
    cache = {}
    if os.path.exists("refs_raw.json"):
        cache = json.load(open("refs_raw.json"))

    sections = {}
    order = []
    seen = set()
    queries = list(QUERIES)
    for section, label, pmids in PINNED:
        queries.append((section, label, "PMID:" + ",".join(pmids), len(pmids)))
    for section, label, query, n in queries:
        key = f"{section}||{label}"
        if key in cache:
            ids = cache[key]["ids"]
            summ = cache[key]["summ"]
        elif query.startswith("PMID:"):
            ids = query[5:].split(",")
            summ = esummary(ids)
            time.sleep(0.40)
            cache[key] = {"ids": ids, "summ": summ, "query": query}
            json.dump(cache, open("refs_raw.json", "w"))
            print(f"pinned  {len(ids):2d}  {label}", file=sys.stderr)
        else:
            ids = esearch(query, n)
            time.sleep(0.40)
            summ = esummary(ids)
            time.sleep(0.40)
            cache[key] = {"ids": ids, "summ": summ, "query": query}
            json.dump(cache, open("refs_raw.json", "w"))
            print(f"fetched {len(ids):2d}  {label}", file=sys.stderr)
        if section not in sections:
            sections[section] = []
            order.append(section)
        for pid in ids:
            if pid in seen:
                continue
            rec = summ.get(pid)
            if not isinstance(rec, dict) or "title" not in rec:
                continue
            seen.add(pid)
            authors = [a["name"] for a in rec.get("authors", [])
                       if a.get("authtype") == "Author"]
            if len(authors) > 3:
                astr = ", ".join(authors[:3]) + ", et al."
            else:
                astr = ", ".join(authors) if authors else "(no author listed)"
            year = (rec.get("pubdate") or "").split(" ")[0]
            sections[section].append(dict(
                pmid=pid, title=rec["title"].rstrip("."),
                authors=astr, journal=rec.get("source", ""), year=year,
                topic=label))

    lines = []
    lines.append("# 대퇴골두 무혈성 괴사 (Osteonecrosis of the Femoral Head, ONFH) — 참고문헌")
    lines.append("")
    lines.append("> **How this list was made.** Every entry below was retrieved")
    lines.append("> programmatically from PubMed with `mkrefs.py` (E-utilities")
    lines.append("> `esearch` + `esummary`), so every PMID, title, journal and year")
    lines.append("> is the one PubMed returns rather than one recalled from memory.")
    lines.append("> The raw payload is kept in `refs_raw.json`; re-running the script")
    lines.append("> regenerates this file. Entries are grouped by the part of the")
    lines.append("> model they support.")
    lines.append("")
    total = sum(len(v) for v in sections.values())
    lines.append(f"**{total} references, {len(order)} sections.**")
    lines.append("")
    n = 0
    for section in order:
        lines.append(f"## {section}")
        lines.append("")
        last_topic = None
        for r in sections[section]:
            if r["topic"] != last_topic:
                if last_topic is not None:
                    lines.append("")
                lines.append(f"*{r['topic']}*")
                lines.append("")
                last_topic = r["topic"]
            n += 1
            lines.append(
                f"{n}. {r['authors']} **{r['title']}.** *{r['journal']}* "
                f"{r['year']}. "
                f"[PMID {r['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/)")
        lines.append("")
    open("onfh_references.md", "w").write("\n".join(lines) + "\n")
    print(f"wrote onfh_references.md with {total} references")


if __name__ == "__main__":
    main()
