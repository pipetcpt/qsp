#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build om_references.md from PubMed.

Every title, journal, year, author and PMID in the bibliography is resolved
live from NCBI esummary -- nothing is written from memory.  Each entry is a
(section, intent, query) triple; the top relevance-ranked PubMed hit for the
query is taken, duplicates are dropped, and the resolved record is printed
together with the INTENT, so a reader can see immediately whether the paper
that came back is actually the one the model leans on.

NOTE on sort=relevance: E-utilities defaults to DATE order, not relevance.
Without the explicit sort the "top hit" is merely the newest paper matching
any of the terms.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = "refs_raw.json"


def _get(url):
    for k in range(6):
        try:
            return urllib.request.urlopen(url, timeout=45).read().decode(
                "utf-8", "replace")
        except Exception:
            time.sleep(2 + 2 * k)
    return ""


def esearch(term, n=1):
    u = B + ("esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d"
             "&term=%s" % (n, urllib.parse.quote(term)))
    try:
        return json.loads(_get(u))["esearchresult"]["idlist"]
    except Exception:
        return []


def esearch_broadening(term, n=10):
    """
    PubMed ANDs every term, so an over-specified query returns ZERO hits
    rather than a near miss -- 29 of the 116 queries in the first build of
    this file came back empty for exactly that reason, and silently dropped
    out of the bibliography.  Drop trailing terms until something matches,
    and RECORD what was actually searched so the broadening is visible.
    """
    words = term.split()
    for cut in range(0, max(len(words) - 3, 0) + 1):
        q = " ".join(words[:len(words) - cut])
        ids = esearch(q, n=n)
        if ids:
            return ids, q
        time.sleep(0.34)
    return [], term


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 150):
        chunk = pmids[i:i + 150]
        u = B + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(chunk)
        try:
            r = json.loads(_get(u))["result"]
        except Exception:
            continue
        for pid in chunk:
            if pid in r:
                out[pid] = r[pid]
        time.sleep(0.34)
    return out


SECTIONS = OrderedDict([
 ("1", ("The disease, its natural history, and the five-phase model",
        "What oral mucositis is, how often it happens, and Sonis's five-phase "
        "framework -- which this model reorganises into ONE fast signalling "
        "loop riding on ONE slow tissue variable.")),
 ("2", ("The renewing epithelium: the conveyor belt that IS the model",
        "Turnover time, clonogen density, transit-amplifying structure and "
        "the post-mitotic barrier.  These papers set k_p, k_shed and the "
        "latent period -- the spine of every equation.")),
 ("3", ("Radiobiology of mucosa: LQ, fractionation, repopulation, gaps",
        "alpha/beta for early-responding mucosa, accelerated repopulation, "
        "and what an unplanned treatment gap costs in tumour control.  These "
        "constrain rad_pot and the BED bookkeeping.")),
 ("4", ("Cytotoxic pharmacokinetics that reach the mucosa",
        "Melphalan, 5-FU/DPD, methotrexate and cisplatin PK.  The mucosal "
        "compartment is driven by these plasma profiles, and the drug "
        "HALF-LIFE is what decides whether cryotherapy can work at all.")),
 ("5", ("NF-kB, ceramide and the cytokine amplification loop",
        "The primary damage response.  These papers justify the NF-kB hub, "
        "the aSMase/ceramide arm that kills POST-MITOTIC cells, and the "
        "TNF -> NF-kB positive feedback whose gain the model holds below 1.")),
 ("6", ("The oral microbiome, the ulcer bed, and the portal of entry",
        "Colonisation of the ulcer, TLR signalling back into NF-kB, and "
        "bacteraemia risk as a function of ulcer area and neutrophil count.")),
 ("7", ("Palifermin / KGF: trials, pharmacology, scheduling, tumour safety",
        "The regeneration-arm drug.  Spielberger 2004 is the model's stage-2 "
        "calibration anchor AND its main held-out prediction, and the "
        "scheduling separation is the paradox the model reproduces.")),
 ("8", ("Oral cryotherapy: the trials, and the half-life argument",
        "The insult-arm intervention.  The pattern of which regimens it helps "
        "is derived in the model from ice-window overlap with the drug's "
        "elimination, not looked up.")),
 ("9", ("Photobiomodulation / low-level laser therapy",
        "The other regeneration-arm intervention, plus its area-independent "
        "analgesic effect.")),
 ("10", ("Benzydamine, glutamine, and the agents that did not work",
         "Anti-inflammatory and epithelial-fuel approaches, and the negative "
         "recommendations (sucralfate, chlorhexidine, 'magic mouthwash') "
         "that a mechanistic model should also be able to accommodate.")),
 ("11", ("Pain, opioid requirement, and analgesic endpoints",
         "Nociception in the ulcer bed and the opioid titration that the "
         "model runs as a closed loop on VAS.")),
 ("12", ("Measurement: WHO, OMAS, CTCAE, patient-reported scales",
         "The instruments.  The WHO scale's saturation at grade 4 is the "
         "reason section 7 of the analysis finds an assay-sensitivity loss "
         "exactly where the disease is worst.")),
 ("13", ("Guidelines and systematic reviews",
         "MASCC/ISOO and other guidance -- used as a CHECK on the model's "
         "predictions, not as an input to them.")),
 ("14", ("Myelosuppression, infection, and the neutrophil coupling",
         "The Friberg semi-mechanistic model and the mucositis-neutropenia-"
         "infection triad.")),
 ("15", ("QSP / systems-pharmacology methodology and mucositis modelling",
         "Prior mathematical models of mucositis and the modelling practice "
         "this file follows.")),
 ("16", ("Consequences: nutrition, hospitalisation, cost, quality of life",
         "What severe mucositis actually costs the patient and the service.")),
])

Q = [
 # -- 1 ---------------------------------------------------------------------
 ("1", "The five-phase biological model of mucositis pathobiology",
  "Sonis pathobiology of mucositis five phase model"),
 ("1", "Mucositis as a biologically complex process, not simple epithelial injury",
  "Sonis new thoughts pathobiology of mucositis complex process"),
 ("1", "Overall incidence and burden of oral mucositis across regimens",
  "incidence oral mucositis chemotherapy radiotherapy epidemiology burden"),
 ("1", "Severe mucositis rates in head and neck radiotherapy, systematic review",
  "Trotti mucositis incidence severity head neck radiotherapy systematic review"),
 ("1", "Oral mucositis after high-dose melphalan and autologous transplantation",
  "oral mucositis high-dose melphalan autologous stem cell transplantation myeloma"),
 ("1", "Alimentary tract mucosal barrier injury as one entity",
  "Blijlevens mucosal barrier injury alimentary tract conditioning"),
 ("1", "Clinical course and time to onset of mucositis after conditioning",
  "time course onset resolution oral mucositis conditioning transplantation"),
 ("1", "Mucositis in paediatric oncology",
  "oral mucositis children paediatric cancer chemotherapy incidence"),
 # -- 2 ---------------------------------------------------------------------
 ("2", "Turnover time of human oral epithelium",
  "Squier biology of oral mucosa turnover time human epithelium"),
 ("2", "Cell kinetics and proliferative compartment of oral mucosa",
  "cell kinetics proliferation compartment oral mucosal epithelium labelling index"),
 ("2", "Epidermal / mucosal stem cells and transit-amplifying cells",
  "epithelial stem cells transit amplifying compartment stratified squamous"),
 ("2", "Regeneration and repopulation of irradiated oral mucosa",
  "Dorr repopulation oral mucosa irradiation three A's"),
 ("2", "Shortening of transit time during epithelial regeneration",
  "accelerated repopulation transit time reduction irradiated epithelium mouse mucosa"),
 ("2", "Clonogen survival and the threshold for mucosal ulceration",
  "clonogen survival threshold denudation ulceration mucosa radiation"),
 ("2", "Basement membrane and integrin biology of oral epithelium",
  "basement membrane laminin integrin alpha6beta4 oral epithelium"),
 ("2", "Mouse tongue model of radiation mucositis",
  "mouse tongue model radiation induced oral mucositis ulceration"),
 # -- 3 ---------------------------------------------------------------------
 ("3", "Linear-quadratic model and alpha/beta for early-responding tissue",
  "linear quadratic model alpha beta ratio early responding normal tissue"),
 ("3", "Accelerated repopulation and the time factor in head and neck cancer",
  "Withers accelerated repopulation time factor head and neck radiotherapy"),
 ("3", "Cost of unplanned treatment gaps in tumour control",
  "treatment gaps prolongation overall treatment time tumour control head neck"),
 ("3", "Hyperfractionation and altered fractionation: mucositis penalty",
  "altered fractionation hyperfractionation acute mucosal toxicity meta-analysis"),
 ("3", "Incomplete sublethal damage repair with short interfraction intervals",
  "incomplete repair sublethal damage interfraction interval hyperfractionation"),
 ("3", "IMRT and mucosa-sparing dose distributions",
  "IMRT oral mucosa sparing dose volume mucositis head neck"),
 ("3", "Dose-volume predictors of severe radiation mucositis",
  "dose volume predictors severe oral mucositis radiotherapy NTCP"),
 ("3", "Total body irradiation and mucosal toxicity in conditioning",
  "total body irradiation conditioning regimen mucositis toxicity"),
 # -- 4 ---------------------------------------------------------------------
 ("4", "Melphalan pharmacokinetics at high dose",
  "melphalan pharmacokinetics high dose 200 mg/m2 AUC clearance transplantation"),
 ("4", "Melphalan exposure and mucositis / outcome relationship",
  "melphalan AUC exposure mucositis toxicity outcome myeloma pharmacokinetic"),
 ("4", "5-fluorouracil pharmacokinetics: bolus versus continuous infusion",
  "5-fluorouracil pharmacokinetics bolus versus continuous infusion clearance"),
 ("4", "DPD deficiency and DPYD genotype-guided dosing",
  "DPYD genotype guided dosing dihydropyrimidine dehydrogenase deficiency fluoropyrimidine"),
 ("4", "Saturable elimination of 5-FU",
  "5-fluorouracil nonlinear saturable elimination Michaelis Menten pharmacokinetics"),
 ("4", "High-dose methotrexate pharmacokinetics and mucosal toxicity",
  "high dose methotrexate pharmacokinetics mucositis leucovorin rescue"),
 ("4", "Methotrexate for GVHD prophylaxis and its mucosal cost",
  "methotrexate graft versus host disease prophylaxis mucositis transplantation"),
 ("4", "Free (ultrafilterable) platinum pharmacokinetics of cisplatin",
  "cisplatin ultrafilterable free platinum pharmacokinetics half-life"),
 # -- 5 ---------------------------------------------------------------------
 ("5", "NF-kB activation in mucositis",
  "NF-kappaB activation oral mucositis radiation chemotherapy"),
 ("5", "Ceramide pathway and acid sphingomyelinase in mucosal injury",
  "ceramide acid sphingomyelinase radiation induced apoptosis mucosa endothelium"),
 ("5", "TNF-alpha in the pathogenesis of mucositis",
  "TNF alpha oral mucositis pathogenesis cytokine"),
 ("5", "IL-1beta, IL-6 and the inflammasome in mucosal injury",
  "IL-1 beta IL-6 inflammasome mucositis chemotherapy"),
 ("5", "Reactive oxygen species as the initiating event",
  "reactive oxygen species initiation mucositis radiation chemotherapy oxidative stress"),
 ("5", "Matrix metalloproteinases in mucosal ulceration",
  "matrix metalloproteinase MMP oral mucositis ulceration"),
 ("5", "Cytokine polymorphisms and susceptibility to severe mucositis",
  "cytokine gene polymorphism susceptibility severe oral mucositis"),
 ("5", "Amifostine and free-radical scavenging as mucosal protection",
  "amifostine radioprotection mucositis head neck randomized"),
 # -- 6 ---------------------------------------------------------------------
 ("6", "The oral microbiome in mucositis",
  "oral microbiome mucositis chemotherapy radiotherapy shift"),
 ("6", "Gram-negative colonisation and endotoxin in the ulcerated mucosa",
  "gram negative bacteria endotoxin ulcerative oral mucositis colonisation"),
 ("6", "TLR4 / MyD88 signalling in mucosal injury",
  "TLR4 MyD88 signalling intestinal oral mucositis chemotherapy"),
 ("6", "Viridans streptococcal bacteraemia after mucositis",
  "viridans streptococci bacteremia oral mucositis neutropenia transplantation"),
 ("6", "Mucositis as a risk factor for bloodstream infection",
  "oral mucositis risk factor bloodstream infection febrile neutropenia"),
 ("6", "Candida and secondary infection of the ulcer bed",
  "oral candidiasis mucositis cancer therapy secondary infection"),
 ("6", "Chlorhexidine is not recommended for mucositis prophylaxis",
  "chlorhexidine oral mucositis radiotherapy prophylaxis not recommended"),
 ("6", "Pre-treatment dental evaluation and basic oral care",
  "basic oral care protocol dental evaluation prevention oral mucositis"),
 # -- 7 ---------------------------------------------------------------------
 ("7", "Palifermin pivotal trial in autologous transplantation",
  "Spielberger palifermin oral mucositis autologous hematopoietic transplantation NEJM"),
 ("7", "Keratinocyte growth factor biology and FGFR2-IIIb",
  "keratinocyte growth factor FGF7 FGFR2 IIIb epithelial receptor"),
 ("7", "Palifermin pharmacokinetics and dosing schedule",
  "palifermin pharmacokinetics dose escalation healthy subjects half-life"),
 ("7", "Palifermin in head and neck chemoradiotherapy",
  "palifermin head and neck cancer chemoradiotherapy randomized mucositis"),
 ("7", "KGF-induced epithelial proliferation and thickening",
  "keratinocyte growth factor epithelial thickening proliferation mucosa animal"),
 ("7", "Nrf2 / cytoprotective mechanisms of KGF",
  "keratinocyte growth factor Nrf2 detoxifying enzymes cytoprotection"),
 ("7", "Scheduling separation between palifermin and cytotoxic therapy",
  "palifermin administration interval chemotherapy 24 hours schedule dependence"),
 ("7", "Tumour safety of growth factors in epithelial cancers",
  "palifermin tumour safety epithelial malignancy growth factor receptor concern"),
 ("7", "Palifermin systematic review and meta-analysis",
  "palifermin meta-analysis oral mucositis prevention systematic review"),
 # -- 8 ---------------------------------------------------------------------
 ("8", "Oral cryotherapy with bolus 5-fluorouracil",
  "Mahood inhibition fluorouracil induced stomatitis oral cryotherapy"),
 ("8", "Oral cryotherapy in high-dose melphalan conditioning",
  "oral cryotherapy high dose melphalan mucositis randomized transplantation"),
 ("8", "Cryotherapy Cochrane / systematic review",
  "oral cryotherapy prevention oral mucositis Cochrane systematic review"),
 ("8", "Mucosal blood flow reduction by intraoral cooling",
  "intraoral cooling oral mucosal blood flow reduction cryotherapy"),
 ("8", "Duration of cryotherapy and its relation to drug half-life",
  "duration oral cryotherapy drug half-life mucositis rationale"),
 ("8", "Cryotherapy is ineffective with continuous-infusion fluoropyrimidine",
  "cryotherapy continuous infusion fluorouracil mucositis ineffective"),
 ("8", "Temperature dependence of alkylating agent reactivity",
  "temperature dependence hydrolysis alkylation melphalan reaction rate"),
 ("8", "Hypothermia and cell cycle arrest in mammalian cells",
  "mild hypothermia cell cycle progression proliferation mammalian cells temperature"),
 # -- 9 ---------------------------------------------------------------------
 ("9", "Photobiomodulation for prevention of oral mucositis",
  "photobiomodulation low level laser therapy prevention oral mucositis randomized"),
 ("9", "Cytochrome c oxidase as the photoacceptor",
  "cytochrome c oxidase photoacceptor photobiomodulation mechanism"),
 ("9", "Photobiomodulation meta-analysis in head and neck radiotherapy",
  "low level laser therapy meta analysis oral mucositis head neck radiotherapy"),
 ("9", "Photobiomodulation dosimetry: fluence and wavelength",
  "photobiomodulation dosimetry fluence wavelength oral mucositis parameters"),
 ("9", "Analgesic effect of photobiomodulation independent of healing",
  "laser therapy analgesia pain reduction oral mucositis mechanism"),
 ("9", "Photobiomodulation safety in the tumour field",
  "photobiomodulation tumour safety oncology laser irradiation concern"),
 # -- 10 --------------------------------------------------------------------
 ("10", "Benzydamine for radiation-induced mucositis",
  "benzydamine hydrochloride oral rinse radiation induced mucositis randomized"),
 ("10", "Benzydamine mechanism: cytokine and NF-kB suppression",
  "benzydamine mechanism TNF alpha cytokine inhibition anti-inflammatory"),
 ("10", "Oral glutamine for mucositis prevention",
  "glutamine oral suspension mucositis prevention randomized Saforis"),
 ("10", "Sucralfate is recommended against",
  "sucralfate oral mucositis radiotherapy no benefit randomized"),
 ("10", "Honey and other low-cost interventions",
  "honey oral mucositis prevention randomized radiotherapy"),
 ("10", "Zinc, vitamin E and antioxidant supplementation",
  "zinc supplementation oral mucositis prevention randomized"),
 ("10", "Growth factors other than KGF (GM-CSF) for mucositis",
  "GM-CSF mouthwash oral mucositis randomized prevention"),
 ("10", "Dusquetide / innate defence regulator in mucositis",
  "dusquetide SGX942 oral mucositis innate defense regulator trial"),
 ("10", "Avasopasem manganese / superoxide dismutase mimetic",
  "avasopasem manganese GC4419 severe oral mucositis randomized trial"),
 # -- 11 --------------------------------------------------------------------
 ("11", "Pain of oral mucositis and its measurement",
  "oral mucositis pain measurement visual analogue scale patient reported"),
 ("11", "Opioid requirement as a mucositis endpoint",
  "opioid requirement morphine equivalent oral mucositis transplantation endpoint"),
 ("11", "Patient-controlled analgesia in mucositis",
  "patient controlled analgesia morphine mucositis stem cell transplantation"),
 ("11", "Topical doxepin rinse for mucositis pain",
  "doxepin rinse oral mucositis pain randomized"),
 ("11", "Topical morphine mouthwash",
  "morphine mouthwash topical oral mucositis pain randomized"),
 ("11", "Peripheral sensitisation: TRPV1, ASIC and Nav channels in oral pain",
  "TRPV1 ASIC sodium channel oral mucosal nociception sensitization"),
 ("11", "Central sensitisation with prolonged oral ulceration",
  "central sensitization chronic oral ulceration pain mechanism"),
 # -- 12 --------------------------------------------------------------------
 ("12", "WHO oral toxicity scale and its use",
  "WHO oral toxicity scale mucositis grading validity"),
 ("12", "Oral Mucositis Assessment Scale (OMAS) validation",
  "Oral Mucositis Assessment Scale OMAS validation multicentre"),
 ("12", "Comparison of mucositis scoring systems",
  "comparison mucositis scoring systems agreement WHO NCI CTC OMAS"),
 ("12", "NCI-CTCAE and clinician- versus patient-reported toxicity",
  "CTCAE clinician reported versus patient reported toxicity agreement oncology"),
 ("12", "Patient-reported oral mucositis instruments (OMDQ, PRO-CTCAE)",
  "oral mucositis daily questionnaire OMDQ patient reported outcome validation"),
 ("12", "Ceiling effects and assay sensitivity of ordinal toxicity scales",
  "ceiling effect ordinal scale responsiveness clinical trial endpoint sensitivity"),
 ("12", "Endpoint choice and statistical power in mucositis trials",
  "endpoint selection sample size mucositis clinical trial design"),
 # -- 13 --------------------------------------------------------------------
 ("13", "MASCC/ISOO clinical practice guidelines for mucositis",
  "MASCC ISOO clinical practice guidelines mucositis cancer therapy"),
 ("13", "MASCC/ISOO 2019/2020 guideline update",
  "MASCC ISOO 2020 update guidelines management mucositis systematic review"),
 ("13", "ASCO / ESMO guidance on mucositis",
  "ESMO clinical practice guidelines oral gastrointestinal mucositis management"),
 ("13", "Systematic review of basic oral care",
  "systematic review basic oral care oral mucositis MASCC"),
 ("13", "Systematic review of natural and miscellaneous agents",
  "systematic review natural agents oral mucositis MASCC ISOO"),
 ("13", "Cochrane overview of interventions for preventing oral mucositis",
  "Cochrane review interventions preventing oral mucositis cancer treatment"),
 # -- 14 --------------------------------------------------------------------
 ("14", "Friberg semi-mechanistic model of myelosuppression",
  "Friberg model of chemotherapy induced myelosuppression semi-mechanistic neutrophil"),
 ("14", "Neutrophil recovery and engraftment after autologous transplantation",
  "neutrophil engraftment kinetics autologous stem cell transplantation melphalan"),
 ("14", "G-CSF effect on the neutropenic window",
  "G-CSF filgrastim duration neutropenia shortening randomized"),
 ("14", "Febrile neutropenia risk models",
  "febrile neutropenia risk model prediction chemotherapy"),
 ("14", "Interaction of mucosal barrier injury with neutropenia",
  "mucosal barrier injury neutropenia interaction infection transplantation"),
 # -- 15 --------------------------------------------------------------------
 ("15", "Mathematical modelling of oral mucositis",
  "mathematical model oral mucositis simulation epithelium radiation"),
 ("15", "Agent-based and computational models of mucosal injury",
  "agent based computational model mucositis epithelium simulation"),
 ("15", "Quantitative systems pharmacology: scope and practice",
  "quantitative systems pharmacology model informed drug development white paper"),
 ("15", "mrgsolve and R-based ODE simulation for PK/PD",
  "mrgsolve R package simulation pharmacometrics ordinary differential equations"),
 ("15", "Model of epithelial radiation response with stem and transit cells",
  "mathematical model stem cell transit amplifying epithelium radiation response"),
 ("15", "Identifiability and parameter estimation in systems models",
  "structural practical identifiability parameter estimation systems biology model"),
 # -- 16 --------------------------------------------------------------------
 ("16", "Cost and resource use attributable to oral mucositis",
  "economic burden cost oral mucositis hospitalization head neck transplantation"),
 ("16", "Nutrition, feeding tubes and weight loss in mucositis",
  "feeding tube gastrostomy weight loss mucositis head neck radiotherapy"),
 ("16", "Quality of life impact of mucositis",
  "quality of life oral mucositis patients cancer therapy impact"),
 ("16", "Unplanned treatment interruptions caused by mucositis",
  "unplanned treatment interruption mucositis radiotherapy head neck compliance"),
 ("16", "Mucositis and survival / treatment compliance",
  "oral mucositis treatment compliance dose intensity survival head neck"),
]


BROADENED = []


def main():
    raw = {}
    if os.path.exists(CACHE):
        raw = json.load(open(CACHE))
    resolved, seen = [], set()
    for i, (sec, intent, q) in enumerate(Q):
        key = q
        if key in raw and raw[key][0]:
            ids, used = raw[key]
        else:
            ids, used = esearch_broadening(q, n=10)
            raw[key] = [ids, used]
            json.dump(raw, open(CACHE, "w"))
            time.sleep(0.34)
        if used != q:
            BROADENED.append((q, used))
        pid = None
        for cand in ids:
            if cand not in seen:
                pid = cand
                break
        if pid is None:
            sys.stderr.write("NO HIT: %s\n" % q)
            continue
        seen.add(pid)
        resolved.append((sec, intent, q, pid))
        sys.stderr.write("%3d/%d  %s\n" % (i + 1, len(Q), pid))

    meta = esummary([r[3] for r in resolved])
    json.dump(meta, open("refs_meta.json", "w"))

    out = []
    out.append("# 구강점막염 (Oral Mucositis) — QSP 모델 참고문헌\n")
    out.append("# Oral Mucositis QSP Model — Annotated Bibliography\n")
    out.append("")
    out.append("**%d references, every one resolved live from NCBI PubMed "
               "esummary by `mkrefs.py`** — no title, author, journal, year "
               "or PMID in this file was written from memory. Each entry "
               "records the INTENT it was retrieved for, so a reader can "
               "check whether the paper that came back is actually the one "
               "the model leans on.\n" % len(resolved))
    out.append("")
    out.append("Regenerate with `python3 mkrefs.py` (needs network access to "
               "`eutils.ncbi.nlm.nih.gov`).\n")
    out.append("")
    out.append("---\n")

    bysec = OrderedDict((k, []) for k in SECTIONS)
    for sec, intent, q, pid in resolved:
        bysec[sec].append((intent, q, pid))

    n = 0
    for sec, (title, blurb) in SECTIONS.items():
        rows = bysec.get(sec, [])
        if not rows:
            continue
        out.append("## %s. %s\n" % (sec, title))
        out.append("*%s*\n" % blurb)
        for intent, q, pid in rows:
            n += 1
            m = meta.get(pid, {})
            t = re.sub(r"\s+", " ", m.get("title", "(title unresolved)"))
            t = t.rstrip(".")
            j = m.get("fulljournalname") or m.get("source") or ""
            yr = (m.get("pubdate", "") or "")[:4]
            au = m.get("authors") or []
            a = au[0]["name"] if au else ""
            if len(au) > 1:
                a += " et al."
            vol = m.get("volume", "")
            pg = m.get("pages", "")
            cite = ", ".join(x for x in [a, j, yr] if x)
            if vol:
                cite += ";%s" % vol
            if pg:
                cite += ":%s" % pg
            out.append("%d. **%s**  " % (n, t))
            out.append("   %s  " % cite)
            out.append("   PMID [%s](https://pubmed.ncbi.nlm.nih.gov/%s/)  "
                       % (pid, pid))
            out.append("   *모델에서의 역할 / role:* %s  " % intent)
            out.append("")
        out.append("")

    out.append("---\n")
    out.append("## 인용 규칙 (how this file was built)\n")
    out.append("- 각 항목은 `(section, intent, query)` 삼중항으로 정의되며, "
               "PubMed `esearch` (`sort=relevance`) 최상위 미중복 히트를 "
               "취하고 `esummary`로 서지사항을 해석합니다.\n")
    out.append("- `sort=relevance`가 없으면 E-utilities는 **날짜순**을 "
               "반환하므로, 질의와 무관한 최신 논문이 최상위로 올라옵니다. "
               "이 옵션은 제거하면 안 됩니다.\n")
    out.append("- 중복 PMID는 자동으로 제거되고 다음 후보로 대체됩니다.\n")
    out.append("- 캐시: `refs_raw.json` (질의→PMID), `refs_meta.json` "
               "(PMID→서지). 삭제 후 재실행하면 전부 다시 조회합니다.\n")
    if BROADENED:
        out.append("")
        out.append("## 질의 완화 기록 (queries that had to be broadened)\n")
        out.append("PubMed는 모든 용어를 AND로 결합하므로 과도하게 특정된 "
                   "질의는 0건을 반환합니다. 아래 질의는 결과가 나올 때까지 "
                   "뒤쪽 용어를 제거했으며, 실제 사용된 질의를 함께 "
                   "기록합니다.\n")
        for orig, used in BROADENED:
            out.append("- `%s`  \n  → 실제 사용: `%s`" % (orig, used))
        out.append("")
    open("om_references.md", "w").write("\n".join(out))
    sys.stderr.write("wrote om_references.md with %d refs\n" % n)


if __name__ == "__main__":
    main()
