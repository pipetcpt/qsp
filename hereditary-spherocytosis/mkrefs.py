#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build hsph_references.md from PubMed.

Every title, journal, year, author and PMID below is resolved LIVE from NCBI
esearch + esummary.  Nothing is written from memory.  Each entry is a
(section, intent, query) triple; the top RELEVANCE-ranked hit is taken and
printed together with the INTENT, so a retrieved paper that does not match
its intent is visible rather than hidden.

E-utilities defaults to DATE order.  sort=relevance is therefore mandatory.
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


Q = [
 # ------------------------------------------------------------------- 1
 ("1. 질환 개관 · What hereditary spherocytosis is",
  "Modern clinical picture, severity classification and management.",
  "hereditary spherocytosis diagnosis and management guideline"),
 ("1", "Epidemiology and the 1/2000 northern European prevalence the model "
       "assumes for its virtual population.",
  "hereditary spherocytosis prevalence epidemiology northern European"),
 ("1", "That the causal chain is vertical-linkage defect -> membrane loss -> "
       "spherocyte -> splenic destruction, which is the chain the model "
       "implements.",
  "hereditary spherocytosis pathophysiology red cell membrane review"),
 ("1", "Severity grading by haemoglobin, reticulocytes and bilirubin -- the "
       "mild/moderate/severe bands the model predicts into.",
  "hereditary spherocytosis clinical severity classification trait mild moderate severe"),
 ("1", "Paediatric natural history and how the phenotype changes with age.",
  "hereditary spherocytosis children natural history"),
 ("1", "British Committee for Standards in Haematology guideline.",
  "guidelines for the diagnosis and management of hereditary spherocytosis"),
 # ------------------------------------------------------------------- 2
 ("2. 유전자 · The genes, and what fraction of HS each explains",
  "ANK1 mutations as the commonest cause and the de novo rate.",
  "ANK1 ankyrin gene mutations hereditary spherocytosis"),
 ("2", "SPTB / beta-spectrin defects.",
  "beta spectrin SPTB mutation hereditary spherocytosis"),
 ("2", "SPTA1 / alpha-spectrin and the recessive severe phenotype.",
  "alpha spectrin SPTA1 recessive severe hereditary spherocytosis"),
 ("2", "SLC4A1 / band 3 deficiency -- the genotype the model predicts "
       "responds differently to splenectomy.",
  "band 3 SLC4A1 deficiency hereditary spherocytosis mutation"),
 ("2", "EPB42 / protein 4.2 deficiency and its Japanese enrichment.",
  "protein 4.2 deficiency hereditary spherocytosis"),
 ("2", "Next-generation sequencing panels and the genotype distribution they "
       "reveal.",
  "next generation sequencing panel hereditary spherocytosis genotype spectrum"),
 ("2", "Spectrin content as a quantitative severity marker -- the basis of "
       "the model's fdef parameter.",
  "spectrin content red cell membrane correlates severity hereditary spherocytosis"),
 ("2", "Genotype-phenotype correlation across the membrane protein defects.",
  "genotype phenotype correlation hereditary spherocytosis membrane protein defect"),
 # ------------------------------------------------------------------- 3
 ("3. 막 구조 · The vertical linkage and why losing it loses area",
  "Ankyrin-1 as the bridge between beta-spectrin and band 3.",
  "ankyrin binds band 3 spectrin erythrocyte membrane"),
 ("3", "Band 3 copy number per red cell -- the model's B0 parameter.",
  "band 3 copy number erythrocyte membrane"),
 ("3", "The complete membrane skeleton architecture: junctional complex, "
       "protein 4.1R, actin, adducin.",
  "erythrocyte membrane skeleton junctional complex protein 4.1"),
 ("3", "Protein 4.2 stabilising the ankyrin-band 3 interaction.",
  "protein 4.2 stabilizes ankyrin band 3 interaction erythrocyte"),
 ("3", "The Rh/RhAG macrocomplex as a second vertical anchor.",
  "Rh RhAG complex band 3 macrocomplex erythrocyte"),
 ("3", "Membrane area incompressibility -- why V_sph is a hard ceiling and "
       "not a soft one.",
  "erythrocyte membrane area expansion modulus"),
 ("3", "Microvesiculation as the mechanism of surface area loss.",
  "membrane loss vesiculation hereditary spherocytosis"),
 ("3", "Band 3-containing vesicles released by spectrin/ankyrin-deficient "
       "cells -- the vesicle composition parameter.",
  "band 3 containing vesicles released hereditary spherocytosis membrane loss"),
 # ------------------------------------------------------------------- 4
 ("4. 기하학 · The geometry: area, volume, and the minimum cylindrical diameter",
  "The original measurement of red cell size and shape distributions and the "
  "minimum cylindrical diameter -- the model's central formula.",
  "Canham Burton distribution of size and shape in populations of normal human red cells"),
 ("4", "Surface area and volume of the normal human red cell.",
  "human erythrocyte surface area volume measurement micropipette"),
 ("4", "Surface-to-volume ratio as the determinant of deformability and of "
       "osmotic fragility.",
  "surface area to volume ratio erythrocyte deformability osmotic fragility"),
 ("4", "Loss of surface area AND volume with normal red cell age -- the "
       "Waugh numbers the normal arm is calibrated to.",
  "Waugh rheologic properties senescent erythrocytes loss of surface area and volume"),
 ("4", "Osmotically inactive (haemoglobin) volume of the red cell.",
  "osmotically inactive volume erythrocyte"),
 ("4", "Ektacytometry and the osmotic gradient deformability profile "
       "(O_min, EI_max, O_hyper).",
  "osmotic gradient ektacytometry deformability profile red cell membrane disorders"),
 ("4", "Direct measurement of red cell deformability in hereditary "
       "spherocytosis.",
  "red cell deformability hereditary spherocytosis ektacytometry measurement"),
 # ------------------------------------------------------------------- 5
 ("5. 비장 · The splenic filter: slits, cords, and transit",
  "Anatomy of the splenic red pulp, the cords of Billroth and the sinus wall.",
  "spleen red pulp cords structure function"),
 ("5", "The interendothelial slit and its width -- the physical filter.",
  "interendothelial slit splenic sinus red cell"),
 ("5", "Biophysics of a red cell squeezing through a spleen-like slit.",
  "red blood cell flow narrow slit spleen model"),
 ("5", "Splenic transit times and the fast/slow (closed/open) circulation "
       "split.",
  "splenic transit time red cells"),
 ("5", "The size of the human splenic red cell pool.",
  "splenic red cell pool sequestration"),
 ("5", "Splenic blood flow as a fraction of cardiac output.",
  "splenic blood flow fraction cardiac output human measurement"),
 ("5", "Mechanical filtration and the physical basis of splenic red cell "
       "retention.",
  "spleen mechanical filtration retention rigid red blood cells physical"),
 ("5", "Splenic conditioning: red cells inside the spleen are the most "
       "spherocytic in the body.",
  "splenic conditioning spherocytes"),
 ("5", "Metabolic environment of the cords: hypoxia, acidosis, glucose "
       "deprivation.",
  "splenic red pulp hypoxia acidosis metabolic stress erythrocyte"),
 ("5", "Pitting function of the spleen (Howell-Jolly bodies).",
  "pitting function spleen Howell-Jolly bodies red cell inclusions"),
 # ------------------------------------------------------------------- 6
 ("6. 청소 기전 · Two clearance mechanisms in one organ",
  "Splenectomy prolongs survival differently in spectrin/ankyrin- versus band "
  "3-deficient HS -- the model's central genotype prediction.",
  "splenectomy prolongs in vivo survival erythrocytes spectrin ankyrin band 3 deficient hereditary spherocytosis"),
 ("6", "Naturally occurring anti-band-3 antibodies and the requirement for "
       "band 3 clustering.",
  "naturally occurring anti band 3 antibodies clustering senescent cell removal"),
 ("6", "Band 3 clustering as the senescent-cell signal.",
  "band 3 clustering senescent cell antigen erythrocyte removal mechanism"),
 ("6", "Red cell IgG quantification by flow cytometry.",
  "quantification of red cell bound IgG flow cytometry molecules per cell"),
 ("6", "CD47-SIRPalpha as the 'don't eat me' signal on red cells.",
  "CD47 SIRP alpha marker of self red blood cell clearance"),
 ("6", "Phosphatidylserine exposure and eryptosis.",
  "phosphatidylserine exposure eryptosis red cell clearance macrophage"),
 ("6", "Red pulp macrophage biology and erythrophagocytosis.",
  "red pulp macrophage erythrophagocytosis SPIC iron recycling"),
 ("6", "Kupffer cell clearance of opsonised red cells -- the arm that "
       "survives splenectomy.",
  "hepatic sequestration IgG coated erythrocytes"),
 ("6", "Normal red cell lifespan and how it is measured.",
  "human red blood cell lifespan biotin label"),
 ("6", "Red cell survival in hereditary spherocytosis before and after "
       "splenectomy.",
  "red cell survival hereditary spherocytosis splenectomy"),
 # ------------------------------------------------------------------- 7
 ("7. 진단 검사 · Diagnostics, and what each one actually measures",
  "EMA (eosin-5-maleimide) binding test performance and cut-off.",
  "eosin 5 maleimide binding test hereditary spherocytosis flow cytometry sensitivity"),
 ("7", "That EMA binds band 3 and the Rh proteins, i.e. reads membrane "
       "content rather than shape.",
  "eosin 5 maleimide binds band 3 Rh protein erythrocyte membrane"),
 ("7", "Osmotic fragility testing, its normal range and its limitations.",
  "osmotic fragility test hereditary spherocytosis sensitivity limitations"),
 ("7", "Incubated osmotic fragility and why ATP depletion exaggerates the "
       "defect.",
  "incubated osmotic fragility hereditary spherocytosis"),
 ("7", "The acidified glycerol lysis time test.",
  "acidified glycerol lysis time test hereditary spherocytosis"),
 ("7", "Cryohaemolysis test.",
  "cryohemolysis test hereditary spherocytosis diagnosis"),
 ("7", "MCHC and RDW as a CBC-only screen.",
  "MCHC RDW screening hereditary spherocytosis red cell distribution width"),
 ("7", "Comparison of screening tests: EMA vs osmotic fragility vs "
       "ektacytometry.",
  "comparison screening tests hereditary spherocytosis EMA osmotic fragility ektacytometry"),
 ("7", "Direct antiglobulin test negativity as the discriminator from "
       "autoimmune haemolytic anaemia.",
  "spherocytes direct antiglobulin test autoimmune hemolytic anemia"),
 # ------------------------------------------------------------------- 8
 ("8. 용혈 지표 · Haemolysis read-outs",
  "Bilirubin production per gram of haemoglobin catabolised (the 34 mg "
  "constant).",
  "bilirubin production rate hemoglobin catabolism"),
 ("8", "Haem oxygenase-1 and the haem degradation pathway.",
  "heme oxygenase 1 biliverdin carbon monoxide iron degradation pathway"),
 ("8", "End-tidal carbon monoxide as a direct measure of haemolytic rate.",
  "end tidal carbon monoxide hemolysis"),
 ("8", "Haptoglobin consumption and its interpretation in extravascular "
       "haemolysis.",
  "haptoglobin extravascular hemolysis interpretation low levels"),
 ("8", "LDH in extravascular versus intravascular haemolysis.",
  "lactate dehydrogenase intravascular versus extravascular hemolysis"),
 ("8", "UGT1A1 kinetics and bilirubin conjugation capacity.",
  "UGT1A1 bilirubin glucuronidation kinetics"),
 # ------------------------------------------------------------------- 9
 ("9. 담석과 길버트 증후군 · Pigment gallstones and UGT1A1",
  "Co-inherited Gilbert syndrome and the risk of gallstones in HS.",
  "UGT1A1 polymorphism gallstones hereditary spherocytosis"),
 ("9", "Prevalence and age distribution of gallstones in HS.",
  "gallstones hereditary spherocytosis prevalence age cholelithiasis"),
 ("9", "Black pigment stone formation: calcium bilirubinate chemistry.",
  "black pigment gallstone formation calcium bilirubinate mechanism"),
 ("9", "Biliary beta-glucuronidase and deconjugation in the bile.",
  "beta glucuronidase bile deconjugation unconjugated bilirubin pigment stones"),
 ("9", "Whether to perform concurrent cholecystectomy at splenectomy.",
  "concomitant cholecystectomy splenectomy hereditary spherocytosis gallstones management"),
 ("9", "UGT1A1*28 allele frequency and its effect on serum bilirubin.",
  "UGT1A1 promoter polymorphism Gilbert syndrome bilirubin"),
 # ------------------------------------------------------------------ 10
 ("10. 비장절제 · Splenectomy: what it does and does not fix",
  "Haematological response to splenectomy in HS.",
  "splenectomy hereditary spherocytosis hematologic response hemoglobin reticulocyte outcome"),
 ("10", "That spherocytes and an abnormal osmotic fragility persist after "
        "splenectomy -- the model's 'cures the anaemia, not the cell' result.",
  "osmotic fragility after splenectomy hereditary spherocytosis"),
 ("10", "Indications for and timing of splenectomy in children.",
  "indications timing splenectomy children hereditary spherocytosis age"),
 ("10", "Laparoscopic versus open splenectomy outcomes.",
  "laparoscopic versus open splenectomy hereditary spherocytosis outcomes"),
 ("10", "Accessory spleens and splenectomy failure.",
  "accessory spleen recurrent hemolysis splenectomy"),
 # ------------------------------------------------------------------ 11
 ("11. 부분 비장절제 · Partial splenectomy and remnant regrowth",
  "Partial (subtotal) splenectomy: haematological benefit with retained "
  "splenic function.",
  "partial splenectomy hereditary spherocytosis outcome subtotal splenic function"),
 ("11", "Long-term follow-up of partial splenectomy and regrowth of the "
        "remnant.",
  "long term follow up partial splenectomy splenic remnant regrowth children"),
 ("11", "Splenic regrowth kinetics after subtotal splenectomy.",
  "splenic regrowth partial splenectomy"),
 ("11", "Comparison of total versus partial splenectomy.",
  "total versus partial splenectomy hereditary spherocytosis comparison"),
 # ------------------------------------------------------------------ 12
 ("12. 비장절제 후 위험 · The price of removing the spleen",
  "Overwhelming post-splenectomy infection: risk and organisms.",
  "overwhelming post splenectomy infection risk incidence encapsulated organisms"),
 ("12", "Vaccination and antibiotic prophylaxis after splenectomy.",
  "vaccination antibiotic prophylaxis asplenia splenectomy recommendations"),
 ("12", "Venous and portal thrombosis after splenectomy.",
  "portal vein thrombosis after splenectomy incidence risk factors"),
 ("12", "Pulmonary hypertension and arterial disease after splenectomy for "
        "haemolytic anaemia.",
  "pulmonary hypertension after splenectomy hemolytic anemia thrombotic risk"),
 ("12", "Long-term cardiovascular and mortality outcomes of splenectomy.",
  "long term mortality vascular outcomes after splenectomy cohort study"),
 # ------------------------------------------------------------------ 13
 ("13. 위기 · Aplastic, megaloblastic and haemolytic crises",
  "Parvovirus B19 aplastic crisis in hereditary spherocytosis.",
  "parvovirus B19 aplastic crisis hereditary spherocytosis"),
 ("13", "The P antigen (globoside) as the parvovirus B19 receptor on "
        "erythroid progenitors.",
  "globoside P antigen receptor parvovirus B19"),
 ("13", "Duration of erythroid arrest in B19 infection -- the 7-10 day "
        "window the model uses.",
  "parvovirus B19 erythroid aplasia reticulocytopenia"),
 ("13", "Folate deficiency and megaloblastic crisis in chronic haemolysis.",
  "folate supplementation chronic hemolytic anemia"),
 ("13", "Haemolytic crisis with intercurrent infection.",
  "hemolytic crisis hereditary spherocytosis children"),
 # ------------------------------------------------------------------ 14
 ("14. 신생아기 · Neonatal hereditary spherocytosis",
  "Neonatal jaundice and anaemia as the presentation of HS.",
  "neonatal jaundice hereditary spherocytosis presentation newborn"),
 ("14", "Late anaemia at 3-8 weeks and the physiological EPO nadir.",
  "anemia infants hereditary spherocytosis erythropoietin"),
 ("14", "Erythropoietin treatment for the transfusion-requiring infant.",
  "recombinant erythropoietin treatment infants hereditary spherocytosis anemia"),
 ("14", "Kernicterus risk and exchange transfusion thresholds.",
  "kernicterus risk exchange transfusion threshold hemolytic disease newborn bilirubin"),
 ("14", "UGT1A1 maturation in the newborn.",
  "UGT1A1 ontogeny newborn bilirubin conjugation"),
 # ------------------------------------------------------------------ 15
 ("15. 수혈과 철 · Transfusion and iron",
  "Transfusion practice in severe HS.",
  "transfusion dependence severe hereditary spherocytosis management"),
 ("15", "Iron loading in non-transfused chronic haemolysis.",
  "iron overload non transfused hemolytic anemia hepcidin erythroferrone"),
 ("15", "Erythroferrone suppression of hepcidin in expanded erythropoiesis.",
  "erythroferrone suppresses hepcidin expanded erythropoiesis"),
 ("15", "HFE co-inheritance and iron loading in HS.",
  "HFE iron overload hereditary spherocytosis"),
 ("15", "Liver iron concentration measurement by MRI R2*.",
  "liver iron concentration MRI R2 validation"),
 ("15", "Deferasirox pharmacology and efficacy.",
  "deferasirox iron chelation efficacy liver iron concentration trial"),
 # ------------------------------------------------------------------ 16
 ("16. 미타피밧 · PKR activation as a non-surgical option",
  "Mitapivat in hereditary spherocytosis.",
  "mitapivat hereditary spherocytosis trial"),
 ("16", "Mitapivat mechanism: allosteric activation of pyruvate kinase R.",
  "mitapivat AG-348 allosteric activator pyruvate kinase R mechanism"),
 ("16", "Mitapivat pharmacokinetics and dosing.",
  "mitapivat pharmacokinetics dose healthy volunteers AG-348"),
 ("16", "Mitapivat in pyruvate kinase deficiency (ACTIVATE) as the "
        "reference efficacy dataset.",
  "mitapivat pyruvate kinase deficiency ACTIVATE trial hemoglobin response"),
 ("16", "Mitapivat in thalassaemia and sickle cell disease -- evidence that "
        "PKR activation helps membranes it was not designed for.",
  "mitapivat thalassemia sickle cell disease trial hemoglobin"),
 ("16", "Red cell ATP and 2,3-DPG responses to PKR activation.",
  "pyruvate kinase activation ATP 2,3-DPG red blood cell response"),
 ("16", "ATP dependence of red cell membrane maintenance and shape.",
  "ATP depletion red cell shape change echinocyte membrane vesiculation"),
 # ------------------------------------------------------------------ 17
 ("17. 조혈 반응 · Erythropoiesis and the ceiling on compensation",
  "The inverse log-linear relationship between haemoglobin and "
  "erythropoietin.",
  "relationship hemoglobin serum erythropoietin log linear anemia"),
 ("17", "Maximal erythropoietic expansion in chronic haemolysis.",
  "erythropoietic expansion chronic hemolytic anemia marrow"),
 ("17", "Reticulocyte shift and the maturation-time correction.",
  "reticulocyte production index maturation time"),
 ("17", "Extramedullary haematopoiesis in chronic haemolytic anaemia.",
  "extramedullary hematopoiesis chronic hemolytic anemia paraspinal masses"),
 ("17", "Erythroferrone, EPO and the erythroid regulator of iron.",
  "erythropoietin erythroferrone erythroid regulator iron absorption"),
 # ------------------------------------------------------------------ 18
 ("18. 모델링 방법 · QSP methods",
  "mrgsolve as an ODE/PK-PD simulation engine in R.",
  "mrgsolve simulation ordinary differential equation pharmacometrics R"),
 ("18", "Quantitative systems pharmacology in model-informed drug "
        "development.",
  "quantitative systems pharmacology model informed drug development review"),
 ("18", "Age-structured population models of red cell survival.",
  "mathematical model erythrocyte survival age structured population kinetics"),
 ("18", "Mechanistic modelling of erythropoiesis and anaemia therapy.",
  "model erythropoiesis erythropoietin pharmacodynamic simulation"),
 ("18", "Good practice for QSP model credibility and verification.",
  "quantitative systems pharmacology model credibility"),
]


def main():
    if os.path.exists(CACHE) and "--refresh" not in sys.argv:
        raw = json.load(open(CACHE))
    else:
        raw = {}
    todo = [q for (_, _, q) in Q if q not in raw]
    print("resolving %d queries (%d cached)" % (len(todo), len(Q) - len(todo)))
    for i, q in enumerate(todo):
        ids = esearch(q, 1)
        raw[q] = ids[0] if ids else ""
        if (i + 1) % 10 == 0:
            print("  %d/%d" % (i + 1, len(todo)))
            json.dump(raw, open(CACHE, "w"))
        time.sleep(0.36)
    json.dump(raw, open(CACHE, "w"))

    pmids = sorted({v for v in raw.values() if v})
    meta = esummary(pmids)
    json.dump(meta, open("refs_meta.json", "w"))

    # ---------------- render
    secs = OrderedDict()
    cur = None
    for sec, intent, q in Q:
        if not sec.isdigit() and " " in sec:
            cur = sec
            secs.setdefault(cur, [])
        secs.setdefault(cur, []).append((intent, q))

    out = []
    out.append("# 유전성 구상적혈구증 (Hereditary Spherocytosis) — 참고문헌\n")
    out.append("> **모든 항목은 NCBI E-utilities(`esearch` + `esummary`, "
               "`sort=relevance`)로 실시간 조회되었습니다.** 제목·저널·연도·"
               "저자·PMID 중 어느 것도 기억에서 작성되지 않았습니다. 각 항목의 "
               "**의도(intent)** 는 모델이 그 문헌에서 무엇을 가져왔는지를 "
               "밝힌 것이며, 검색 결과가 의도와 맞지 않는 경우가 한눈에 보이도록 "
               "함께 표기합니다.\n")
    out.append("> Rebuild with `python3 mkrefs.py --refresh`.\n")
    n_tot = 0
    body = []
    for sec, items in secs.items():
        body.append("\n## %s\n" % sec)
        for intent, q in items:
            pid = raw.get(q, "")
            m = meta.get(pid, {})
            if not m:
                body.append("- *(no PubMed hit)* — intent: %s  \n  `query: %s`\n"
                            % (intent, q))
                continue
            n_tot += 1
            au = m.get("authors", [])
            first = au[0]["name"] if au else "—"
            alist = first + (" et al." if len(au) > 2 else
                             (" & " + au[1]["name"] if len(au) == 2 else ""))
            title = re.sub(r"\s+", " ", m.get("title", "")).strip().rstrip(".")
            jr = m.get("source", "")
            yr = (m.get("pubdate", "") or "")[:4]
            body.append(
                "- **%s** — %s. *%s* %s. "
                "[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)  \n"
                "  <sub>모델에서의 역할: %s</sub>\n"
                % (title, alist, jr, yr, pid, pid, intent))
    out.append("총 **%d편**의 문헌이 %d개 섹션에 정리되어 있습니다.\n"
               % (n_tot, len(secs)))
    out.extend(body)
    open("hsph_references.md", "w").write("".join(out))
    print("wrote hsph_references.md with %d resolved references" % n_tot)


if __name__ == "__main__":
    main()
