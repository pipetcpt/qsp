#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build hdfn_references.md from PubMed.

Every title, journal, year, author list and PMID in the bibliography is
resolved live from NCBI esearch + esummary -- none of it is written from
memory.  Each entry is a (section, intent, query) triple: the top
RELEVANCE-ranked hit is taken and printed together with the INTENT, so that a
retrieved paper which does not match its intent is visible at a glance instead
of hidden.

E-utilities defaults to DATE order, not relevance.  sort=relevance is therefore
mandatory: without it the "top hit" is simply the newest paper matching any
term.
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
            return urllib.request.urlopen(url, timeout=60).read().decode(
                "utf-8", "replace")
        except Exception:
            time.sleep(1.5 + 2 * k)
    return ""


def esearch(term, n=1):
    u = B + ("esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d"
             "&term=%s" % (n, urllib.parse.quote(term)))
    try:
        return json.loads(_get(u))["esearchresult"]["idlist"]
    except Exception:
        return []


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 120):
        chunk = pmids[i:i + 120]
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


# (section, intent, query) -- the intent says what the MODEL takes from it
Q = [
 # ---------------------------------------------------------------- 1
 ("1. The disease, and why it is a transport problem before it is an immune one",
  "The modern clinical picture and management of HDFN.",
  "hemolytic disease of the fetus and newborn review management"),
 ("1", "That the causal chain is maternal alloantibody -> transplacental "
       "transfer -> fetal red cell destruction -> anaemia, which is the chain "
       "the model implements.",
  "red cell alloimmunization pregnancy pathophysiology fetal anemia"),
 ("1", "Epidemiology and residual burden in the prophylaxis era.",
  "epidemiology severe hemolytic disease fetus newborn incidence prophylaxis era"),
 ("1", "That non-D antibodies now account for a large share of severe disease.",
  "non-RhD red cell antibodies severe hemolytic disease fetus newborn"),
 ("1", "Guideline framing of surveillance thresholds and referral.",
  "management alloimmunization pregnancy guideline society maternal fetal medicine"),
 # ---------------------------------------------------------------- 2
 ("2. The antigen: RHD genetics, antigen density, and who becomes immunised",
  "RHD gene structure and the molecular basis of the D-negative phenotype.",
  "RHD gene molecular basis D negative phenotype"),
 ("2", "D antigen site density per red cell by genotype -- the site_D "
       "parameter that scales opsonisation.",
  "number of D antigen sites red cell"),
 ("2", "Sensitisation risk per pregnancy without prophylaxis (the 16% "
       "calibration target).",
  "risk of Rh immunization pregnancy without prophylaxis Bowman"),
 ("2", "That ABO incompatibility protects against D sensitisation, because "
       "incompatible fetal cells are cleared before they can prime.",
  "ABO incompatibility protection against Rh immunization"),
 ("2", "HLA class II restriction of the anti-D response and the 'non-responder' "
       "phenotype.",
  "HLA class II anti-D immunization responder"),
 ("2", "Fetomaternal haemorrhage volume distribution -- the log-normal the "
       "prophylaxis model integrates over.",
  "fetomaternal hemorrhage volume distribution flow cytometry delivery"),
 # ---------------------------------------------------------------- 3
 ("3. Anti-D prophylaxis: a stoichiometric race, and where the residual risk is",
  "Practical dosing, the 20 ug per mL of fetal red cells rule, and the 72-hour "
  "window.",
  "prevention fetomaternal rhesus D alloimmunization practical aspects dose"),
 ("3", "Antenatal prophylaxis at 28 weeks and the residual sensitisation rate "
       "of 0.1-0.4% that the model predicts rather than fits.",
  "routine antenatal anti-D prophylaxis residual immunization"),
 ("3", "Postpartum-only prophylaxis leaving ~1.6% -- the second calibration "
       "target.",
  "anti-D immunoglobulin postpartum prevention Rh immunization trial"),
 ("3", "Kleihauer-Betke and flow cytometry quantification of FMH to size the "
       "dose.",
  "Kleihauer Betke flow cytometry quantification fetomaternal hemorrhage dose"),
 ("3", "Cell-free fetal DNA RHD genotyping to target prophylaxis.",
  "cell free fetal DNA RHD genotyping targeted antenatal prophylaxis"),
 ("3", "Mechanism of anti-D-mediated suppression: clearance versus immune "
       "deviation (the dev parameter).",
  "mechanism of action anti-D immunoglobulin immune suppression antibody mediated"),
 ("3", "Timing: efficacy as a function of delay after the sensitising event.",
  "anti-D immunoglobulin 72 hours sensitizing event"),
 # ---------------------------------------------------------------- 4
 ("4. Measuring the antibody: titre, quantitation, subclass, function",
  "Quantitation in IU/mL and the <4 / 4-15 / >15 risk bands used as the "
  "destruction-potency anchor.",
  "anti-D quantitation IU/mL prediction severity hemolytic disease"),
 ("4", "IgG subclass composition and why IgG1 and IgG3 matter differently -- "
       "the pot3 parameter.",
  "IgG subclass anti-D IgG1 IgG3 severity hemolytic disease newborn"),
 ("4", "ADCC and monocyte monolayer assays outperforming serology.",
  "antibody dependent cellular cytotoxicity monocyte monolayer assay predict HDN"),
 ("4", "Critical titre and its limitations as a trigger for surveillance.",
  "critical antibody titer alloimmunized pregnancy"),
 ("4", "Titre rise after intrauterine transfusion -- the boost the model gives "
       "every procedure.",
  "antibody formation after intrauterine transfusion mother"),
 # ---------------------------------------------------------------- 5
 ("5. The placental conveyor: FcRn, gestational age, and subclass selectivity",
  "The fetal:maternal IgG ratio through gestation -- the two numbers the "
  "conveyor is calibrated to (0.075 at 19.5 wk, 1.25 at 39 wk).",
  "evolution of maternofetal transport of immunoglobulins during human pregnancy"),
 ("5", "FcRn as the transporter, and the mechanism of transcytosis across the "
       "syncytiotrophoblast.",
  "FcRn neonatal Fc receptor placental transfer IgG transcytosis mechanism"),
 ("5", "Gestational-age dependence of transfer, i.e. the exponential capacity "
       "term.",
  "gestational age dependence placental transfer IgG exponential increase"),
 ("5", "FcRn-mediated IgG recycling and the 21-day half-life it produces.",
  "FcRn IgG homeostasis half life recycling catabolism"),
 ("5", "Saturability of placental transfer by bulk IgG -- the basis of the "
       "IVIG competition mechanism.",
  "saturation placental IgG transfer competition high dose immunoglobulin"),
 ("5", "Maternal plasma volume expansion in pregnancy, which alone accounts "
       "for much of the fall in maternal IgG concentration.",
  "plasma volume expansion normal pregnancy magnitude"),
 # ---------------------------------------------------------------- 6
 ("6. Fetal erythropoiesis, and what a transfusion does to it",
  "Reference range for fetal haemoglobin by gestational age.",
  "fetal hemoglobin concentration reference range gestation cordocentesis normal"),
 ("6", "Fetal erythropoietin response to anaemia.",
  "erythropoietin concentration in fetal blood anaemia cordocentesis"),
 ("6", "Extramedullary and hepatic erythropoiesis in severe fetal anaemia.",
  "extramedullary hepatic erythropoiesis fetal anemia erythroblastosis"),
 ("6", "Reticulocytosis and nucleated red cells as markers of the "
       "compensatory response.",
  "nucleated red blood cells reticulocyte count fetal anemia alloimmunization"),
 ("6", "Suppression of fetal erythropoiesis by intrauterine transfusion.",
  "suppression fetal erythropoiesis after intrauterine transfusion reticulocyte"),
 ("6", "Fetal red cell lifespan and its difference from the adult.",
  "fetal erythrocyte survival lifespan"),
 # ---------------------------------------------------------------- 7
 ("7. MCA-PSV: the threshold this model derives rather than fits",
  "The original demonstration: 100% sensitivity for moderate/severe anaemia, "
  "12% false positives, and the 1.5 MoM threshold.",
  "noninvasive diagnosis fetal anemia Doppler middle cerebral artery alloimmunization"),
 ("7", "Reference ranges for MCA-PSV and the exp(2.31+0.046*GA) median.",
  "middle cerebral artery peak systolic velocity reference range gestational age"),
 ("7", "Physiological basis: increased cardiac output and reduced viscosity in "
       "fetal anaemia.",
  "hemodynamic mechanism increased middle cerebral artery velocity fetal anemia viscosity"),
 ("7", "Performance of MCA-PSV after the first transfusion, when it degrades -- "
       "why repeat transfusions are scheduled on the calendar.",
  "Performance of middle cerebral artery doppler for prediction of recurrent fetal anemia"),
 ("7", "Systematic review / meta-analysis of MCA-PSV diagnostic accuracy.",
  "meta-analysis accuracy middle cerebral artery peak systolic velocity fetal anemia"),
 ("7", "Fetal cerebral oxygen delivery and its defence, which is the "
       "assumption behind do2_alpha = 1.",
  "fetal cerebral blood flow anemia oxygen delivery"),
 # ---------------------------------------------------------------- 8
 ("8. Amniotic fluid bilirubin: the test the model shows is wrong in Kell disease",
  "The Liley curve and amniotic fluid dOD450.",
  "amniotic fluid optical density 450 Rh disease Liley"),
 ("8", "Extension of the curve into the second trimester (Queenan).",
  "Queenan curve amniotic fluid deviation 450 alloimmunization"),
 ("8", "That dOD450 measures bilirubin, i.e. haem released, not haemoglobin -- "
       "the reason it fails in Kell disease.",
  "Kell alloimmunization amniotic fluid bilirubin unreliable"),
 ("8", "Comparison of amniocentesis and Doppler for the same decision.",
  "comparison amniocentesis Doppler middle cerebral artery management alloimmunization"),
 # ---------------------------------------------------------------- 9
 ("9. Intrauterine transfusion: mass balance, intervals, decline, risk",
  "The measured haemoglobin decline between the first and second transfusion "
  "(0.40 g/dL/day, SD 0.25) -- the model's headline prediction target.",
  "prediction rate of decline fetal hemoglobin between first second transfusion"),
 ("9", "Complication and loss rates per procedure and per fetus.",
  "complications intrauterine intravascular blood transfusion lessons learned 1678"),
 ("9", "Volume calculation and target haematocrit for intravascular transfusion.",
  "intrauterine transfusion volume hematocrit formula"),
 ("9", "Fetoplacental blood volume estimation as a function of fetal weight.",
  "fetoplacental blood volume estimation fetal weight gestational age"),
 ("9", "Survival and outcome of hydropic versus non-hydropic fetuses treated "
       "by IUT.",
  "outcome intrauterine transfusion hydrops survival red cell alloimmunization"),
 ("9", "Rate of haemoglobin fall in hydropic versus non-hydropic fetuses.",
  "effect of fetal hydrops rate of fall of hemoglobin after intravascular transfusion"),
 ("9", "Intrahepatic versus cord insertion route and technique refinements.",
  "intrahepatic vein intrauterine transfusion"),
 ("9", "Timing of delivery after a course of intrauterine transfusions.",
  "late versus early intrauterine blood transfusion fetal anemia neonatal outcome"),
 # ---------------------------------------------------------------- 10
 ("10. Hydrops fetalis: the Starling balance and the haemoglobin threshold",
  "Haemoglobin deficit at which hydrops appears.",
  "fetal hydrops hemoglobin deficit alloimmunization"),
 ("10", "Umbilical venous pressure in anaemic and hydropic fetuses -- rising, "
        "then falling in the extreme, which the model reproduces.",
  "umbilical venous pressure normal growth retarded anemic fetuses"),
 ("10", "Hypoalbuminaemia and colloid osmotic pressure in the hydropic fetus.",
  "fetal plasma protein colloid osmotic pressure hydrops"),
 ("10", "Fetal cardiovascular response to anaemia: cardiac output and its "
        "reserve.",
  "fetal cardiac output response to anemia cardiovascular hemodynamics"),
 ("10", "Atrial natriuretic factor and volume overload in the anaemic fetus.",
  "atrial natriuretic factor anemic hydropic fetuses concentration"),
 ("10", "Lymphatic return and interstitial pressure as the safety factor "
        "against oedema.",
  "lymph flow interstitial pressure safety factor against edema"),
 # ---------------------------------------------------------------- 11
 ("11. FcRn blockade: nipocalimab, the UNITY trial, and the class",
  "The UNITY phase 2 result: 7/13 (54%) live birth >= 32 wk without IUT, no "
  "hydrops.",
  "nipocalimab early-onset severe hemolytic disease of the fetus and newborn"),
 ("11", "Nipocalimab pharmacology: FcRn blockade, IgG lowering, and "
        "pharmacokinetics.",
  "nipocalimab pharmacokinetics pharmacodynamics IgG"),
 ("11", "FcRn antagonists as a class and their use in IgG-mediated disease.",
  "FcRn antagonist therapeutic class IgG mediated autoimmune disease review"),
 ("11", "Efgartigimod pharmacology as the comparator FcRn antagonist.",
  "efgartigimod FcRn antagonist IgG reduction pharmacodynamics"),
 ("11", "Whether an FcRn blocker crosses the placenta, and the consequences "
        "for neonatal IgG.",
  "FcRn antagonist pregnancy placental transfer"),
 ("11", "Safety considerations of lowering maternal and neonatal IgG.",
  "hypogammaglobulinemia infection risk FcRn antagonist safety"),
 ("11", "The 2024-2025 clinical development landscape for HDFN.",
  "emerging treatment hemolytic disease fetus newborn"),
 # ---------------------------------------------------------------- 12
 ("12. IVIG and plasmapheresis: mechanism, and the evidence they do not prevent IUT",
  "Combined plasmapheresis and IVIG: all nine pregnancies still required IUT "
  "(median 4).",
  "combined plasmapheresis intravenous immune globulin severe maternal red cell alloimmunization"),
 ("12", "IVIG mechanisms in antibody-mediated cytopenias: FcRn competition and "
        "FcgammaR blockade.",
  "mechanisms of action intravenous immunoglobulin autoimmune"),
 ("12", "Systematic review of maternal IVIG in HDFN.",
  "systematic review maternal intravenous immunoglobulin hemolytic disease fetus"),
 ("12", "Plasmapheresis kinetics and antibody rebound.",
  "plasma exchange immunoglobulin removal kinetics"),
 ("12", "Neonatal IVIG for isoimmune haemolytic jaundice.",
  "intravenous immunoglobulin neonatal isoimmune hemolytic jaundice exchange transfusion"),
 # ---------------------------------------------------------------- 13
 ("13. Other antibodies: Kell suppresses erythropoiesis, ABO barely does anything",
  "Anti-Kell causes anaemia by suppressing erythropoiesis rather than by "
  "haemolysis -- the kell_kill parameter.",
  "anti-Kell antibodies inhibit erythroid progenitor growth suppression erythropoiesis"),
 ("13", "Clinical severity and management of Kell alloimmunisation.",
  "Kell alloimmunization pregnancy severity management anti-K"),
 ("13", "Anti-c and anti-E as causes of significant HDFN.",
  "anti-c anti-E alloimmunization severity hemolytic disease fetus newborn"),
 ("13", "Why ABO HDFN is a neonatal jaundice and not a fetal anaemia.",
  "ABO hemolytic disease of the newborn pathogenesis"),
 ("13", "Rare antibodies and antigen-negative donor sourcing.",
  "rare blood group antibodies pregnancy antigen negative donor transfusion"),
 # ---------------------------------------------------------------- 14
 ("14. The newborn: bilirubin, phototherapy, exchange transfusion",
  "Bilirubin production rate in the newborn and the 34 mg per g of haemoglobin "
  "stoichiometry.",
  "bilirubin production rate newborn carbon monoxide heme catabolism mg per kg"),
 ("14", "UGT1A1 ontogeny -- the clock that makes neonatal haemolysis a "
        "bilirubin problem.",
  "UGT1A1 ontogeny developmental expression glucuronidation bilirubin infant"),
 ("14", "Bilirubin distribution volume and albumin binding.",
  "bilirubin albumin binding reserve capacity distribution volume neonate"),
 ("14", "Phototherapy mechanism, dose-response and irradiance.",
  "phototherapy mechanism irradiance dose response neonatal hyperbilirubinemia"),
 ("14", "Exchange transfusion efficiency and post-exchange rebound.",
  "exchange transfusion neonatal hyperbilirubinemia efficacy"),
 ("14", "Thresholds for phototherapy and exchange in isoimmune haemolysis.",
  "guideline management hyperbilirubinemia newborn 35 weeks thresholds exchange"),
 ("14", "Directed antiglobulin test and the diagnosis of neonatal HDFN.",
  "direct antiglobulin test neonatal hemolytic disease diagnosis"),
 # ---------------------------------------------------------------- 15
 ("15. Late anaemia: the nadir that arrives six weeks after delivery",
  "Late hyporegenerative anaemia after intrauterine transfusion and the need "
  "for top-up transfusions.",
  "late anemia intrauterine transfusion erythropoiesis suppression"),
 ("15", "Erythropoietin for late anaemia of HDFN.",
  "erythropoietin treatment late anemia hemolytic disease newborn transfusion"),
 ("15", "Persistence of maternal antibody in the newborn and duration of "
        "haemolysis.",
  "persistent hemolysis infant maternal antibody months"),
 ("15", "Iron overload after multiple intrauterine and neonatal transfusions.",
  "iron overload ferritin after intrauterine transfusions neonate"),
 ("15", "Neonatal complications after IUT: neutropenia, thrombocytopenia, "
        "cholestasis.",
  "neonatal morbidity after intrauterine transfusion"),
 # ---------------------------------------------------------------- 16
 ("16. Long-term outcome",
  "LOTUS: 4.8% neurodevelopmental impairment, and severe hydrops as the "
  "dominant preoperative risk factor (OR 11.2).",
  "long-term neurodevelopmental outcome intrauterine transfusion LOTUS"),
 ("16", "Design of the long-term follow-up cohort.",
  "long term follow up after intrauterine transfusions LOTUS study design"),
 ("16", "Bilirubin neurotoxicity and kernicterus spectrum disorder.",
  "bilirubin induced neurologic dysfunction kernicterus spectrum disorder"),
 ("16", "Maternal red cell alloimmunisation as a lifelong transfusion problem.",
  "maternal alloimmunization additional antibodies after intrauterine transfusion"),
 # ---------------------------------------------------------------- 17
 ("17. Methods: QSP, mrgsolve, and model-informed drug development",
  "mrgsolve as the ODE engine.",
  "mrgsolve simulation pharmacometric models R"),
 ("17", "QSP model credibility and validation practice.",
  "credibility assessment of quantitative systems pharmacology models"),
 ("17", "Model-informed drug development in pregnancy, where trials are small.",
  "model informed drug development pediatric rare disease"),
 ("17", "Physiologically based modelling of maternal-fetal drug and antibody "
        "transfer.",
  "physiologically based pharmacokinetic model maternal fetal placental transfer antibody"),
 ("17", "Virtual population methods for small single-arm trials.",
  "virtual patient population clinical trial simulation"),
]


def main():
    cache = {}
    if os.path.exists(CACHE):
        try:
            cache = json.load(open(CACHE))
        except Exception:
            cache = {}
    resolved = OrderedDict()
    sections = OrderedDict()
    cur = None
    for row in Q:
        sec, intent, query = row
        if len(sec) > 3:
            cur = sec
        if cur not in sections:
            sections[cur] = []
        if query in cache:
            ids = cache[query]
        else:
            ids = esearch(query, 3)
            cache[query] = ids
            time.sleep(0.34)
        sections[cur].append((intent, query, ids))
    json.dump(cache, open(CACHE, "w"), indent=1)

    allids = []
    for sec, rows in sections.items():
        for intent, query, ids in rows:
            allids.extend(ids[:3])
    meta = esummary(sorted(set(allids)))
    json.dump(meta, open("refs_meta.json", "w"))

    used = set()
    n = 0
    out = []
    out.append("# Hemolytic Disease of the Fetus and Newborn — QSP model references\n")
    out.append(
        "Every entry below was resolved live from NCBI PubMed (`esearch` + "
        "`esummary`) by [`mkrefs.py`](mkrefs.py); titles, journals, years, "
        "authors and PMIDs are machine-transcribed rather than recalled. Each "
        "entry carries the **intent** — what the model actually takes from it — "
        "so a retrieved paper that does not match its intent is visible at a "
        "glance rather than hidden.\n")
    out.append(
        "> **Eight numbers in this model were fitted; everything else was "
        "predicted.** The fitted eight are: the fetal:maternal IgG ratio at "
        "19.5 and 39 weeks (Malek 1996); the gestational age of the first "
        "intrauterine transfusion at a maternal anti-D of 15 IU/mL (26 weeks, "
        "the cohort mean of Nishie 2012); the population sensitisation risk "
        "without prophylaxis (16%) and with postpartum-only prophylaxis (1.6%); "
        "the peak and timing of physiological jaundice in a **healthy** term "
        "newborn (8 mg/dL at ~4 days); and the progenitor expansion at which "
        "extramedullary erythropoiesis is recruited, placed so that overt "
        "ascites appears at 5–6 g/dL. Sections 7, 9, 11 and 13–16 are, with "
        "those exceptions, out-of-sample prediction targets.\n")

    for sec, rows in sections.items():
        out.append("\n## %s\n" % sec)
        for intent, query, ids in rows:
            pid = None
            for cand in ids:
                if cand not in used and cand in meta:
                    pid = cand
                    break
            if pid is None:
                for cand in ids:
                    if cand in meta:
                        pid = cand
                        break
            if pid is None:
                continue
            used.add(pid)
            m = meta[pid]
            au = m.get("authors", [])
            names = [a.get("name", "") for a in au]
            if len(names) == 0:
                who = ""
            elif len(names) == 1:
                who = names[0]
            elif len(names) <= 3:
                who = ", ".join(names)
            else:
                who = names[0] + " et al."
            title = re.sub(r"\s+", " ", m.get("title", "")).strip().rstrip(".")
            jour = m.get("source", "")
            yr = (m.get("pubdate", "") or "")[:4]
            n += 1
            out.append("%d. **%s**  \n   %s %s. *%s* %s. "
                       "[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)"
                       % (n, intent, who, title, jour, yr, pid, pid))
    out.append("\n---\n")
    out.append("%d references, resolved %s from PubMed. Queries and the raw "
               "esearch/esummary payloads are in "
               "[`refs_raw.json`](refs_raw.json) and "
               "[`refs_meta.json`](refs_meta.json) so that every entry can be "
               "re-derived." % (n, time.strftime("%Y-%m-%d")))
    open("hdfn_references.md", "w").write("\n".join(out) + "\n")
    print("wrote hdfn_references.md with %d references" % n)


if __name__ == "__main__":
    main()
