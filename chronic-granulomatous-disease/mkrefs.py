#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build cgd_references.md from PubMed.

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
 ("1. 질환 개관 · What chronic granulomatous disease is",
  "Modern clinical picture, genetics and management of CGD as a whole.",
  "chronic granulomatous disease review diagnosis management"),
 ("1", "The 368-patient US national registry: the organism spectrum and the "
       "site distribution the model's Section B is scored against.",
  "chronic granulomatous disease national registry 368 patients United States"),
 ("1", "The European registry series, for infection incidence per patient-year.",
  "van den Berg chronic granulomatous disease the European experience"),
 ("1", "Long-term survival and cause of death in a large modern US cohort.",
  "Marciano common severe infections chronic granulomatous disease"),
 ("1", "Incidence and birth prevalence of CGD.",
  "chronic granulomatous disease birth incidence United States registry"),
 ("1", "Discovery that the defect is an absent respiratory burst, not a missing enzyme.",
  "quantitative study of the dysfunction of leukocytes chronic granulomatous disease"),
 # ------------------------------------------------------------------- 2
 ("2. 유전자와 NADPH 산화효소 · The genes and the oxidase",
  "CYBB / gp91phox as the X-linked cause.",
  "CYBB gp91phox mutation X-linked chronic granulomatous disease"),
 ("2", "NCF1 / p47phox and the GT deletion pseudogene recombination.",
  "NCF1 p47phox GT deletion pseudogene chronic granulomatous disease"),
 ("2", "NCF2 / p67phox deficiency.",
  "NCF2 p67phox deficiency chronic granulomatous disease mutation"),
 ("2", "CYBA / p22phox deficiency.",
  "CYBA p22phox deficiency chronic granulomatous disease"),
 ("2", "NCF4 / p40phox deficiency and its colitis-predominant phenotype.",
  "NCF4 p40phox deficiency chronic granulomatous disease colitis"),
 ("2", "CYBC1 / EROS as the chaperone required for gp91phox-p22phox stability.",
  "CYBC1 EROS chaperone gp91phox chronic granulomatous disease"),
 ("2", "RAC2 and the small GTPase requirement for oxidase activation.",
  "Rac2 GTPase NADPH oxidase activation neutrophil"),
 ("2", "Structure of the assembled NOX2 complex.",
  "structure NOX2 NADPH oxidase complex cryo-EM p22phox"),
 ("2", "Electron transfer through FAD and the two haems of gp91phox.",
  "flavocytochrome b558 heme FAD electron transfer phagocyte oxidase"),
 ("2", "p47phox phosphorylation as the activation switch.",
  "phosphorylation of p47phox sites activation phagocyte NADPH oxidase"),
 ("2", "The mutation database and genotype distribution across CGD genes.",
  "hematologically important mutations X-linked chronic granulomatous disease update"),
 # ------------------------------------------------------------------- 3
 ("3. 잔여 활성이 결과를 정한다 · Residual oxidase activity predicts survival",
  "THE central clinical paper: residual superoxide production, not genotype, "
  "predicts survival. The model's whole phi axis is this result.",
  "Kuhns residual NADPH oxidase and survival in chronic granulomatous disease"),
 ("3", "The DHR-123 flow assay and how CGD is actually diagnosed.",
  "dihydrorhodamine 123 flow cytometry diagnosis chronic granulomatous disease"),
 ("3", "That the DHR histogram distinguishes X-linked carriers (bimodal) from "
       "hypomorphs (shifted) -- the model's Section D claim.",
  "dihydrorhodamine carrier detection X-linked chronic granulomatous disease bimodal"),
 ("3", "X-linked carriers: lyonisation, the fraction of oxidase-normal "
       "neutrophils, and when carriers get disease.",
  "X-linked carriers chronic granulomatous disease lyonization clinical manifestations"),
 ("3", "Female carriers with skewed X-inactivation presenting as CGD.",
  "extreme skewing X inactivation female carrier chronic granulomatous disease infection"),
 ("3", "Nitroblue tetrazolium and ferricytochrome c reduction as quantitative "
       "superoxide assays.",
  "nitroblue tetrazolium ferricytochrome c reduction superoxide neutrophil assay"),
 ("3", "Whole-cell superoxide production rate of stimulated human neutrophils, "
       "the independent cross-check on the model's R_ox.",
  "quantitation superoxide anion release human neutrophils phorbol myristate acetate"),
 # ------------------------------------------------------------------- 4
 ("4. 식포 화학 · The phagosome as a chemical reactor",
  "THE kernel's source: a kinetic model of superoxide and myeloperoxidase "
  "reactions in the neutrophil phagosome, with the ~1e8 HOCl per bacterium.",
  "Winterbourn modeling the reactions of superoxide and myeloperoxidase in the neutrophil phagosome"),
 ("4", "Reconciling the chemistry and biology of neutrophil oxidant production.",
  "Winterbourn Kettle redox reactions and microbial killing in the neutrophil phagosome"),
 ("4", "Myeloperoxidase concentration and the azurophil granule content.",
  "myeloperoxidase content of human neutrophil azurophil granules"),
 ("4", "MPO Compound I formation rate constant with H2O2.",
  "Furtmuller reaction of myeloperoxidase compound I hydrogen peroxide rate constant"),
 ("4", "MPO Compound I + chloride -> HOCl rate constant.",
  "myeloperoxidase compound I chloride oxidation rate constant hypochlorous acid"),
 ("4", "Superoxide reactions with MPO: Compound III formation and the "
       "SOD-like cycle the model relies on for its H2O2 supply.",
  "Kettle superoxide myeloperoxidase compound III interaction neutrophil"),
 ("4", "MPO Compound II and its reduction by superoxide.",
  "superoxide converts myeloperoxidase compound II to the ferric enzyme"),
 ("4", "Spontaneous dismutation rate of superoxide and its pH dependence.",
  "reactivity of HO2/O2- radicals in aqueous solution"),
 ("4", "Phagosomal chloride concentration and its depletion.",
  "chloride concentration in neutrophil phagosomes measured"),
 ("4", "Rate constants for HOCl with methionine, cysteine and amines.",
  "rate constants hypochlorous acid reaction methionine cysteine amino acids"),
 ("4", "Chlorination of bacterial proteins as the measurable footprint of HOCl.",
  "chlorination bacterial proteins neutrophil phagosome 3-chlorotyrosine"),
 ("4", "How much of the phagosomal oxygen consumption ends as HOCl.",
  "myeloperoxidase chlorination activity stimulated neutrophils taurine chloramine"),
 # ------------------------------------------------------------------- 5
 ("5. 식포 pH와 이온 · pH, potassium and the granule matrix",
  "Segal's measurement that the normal phagosome alkalinises and the CGD "
  "phagosome does not -- one side of a live controversy.",
  "Segal kinetics of fusion of the cytoplasmic granules with phagocytic vacuoles pH"),
 ("5", "The opposing measurement: CGD phagosomes acidify normally.",
  "Jankowski Grinstein determinants of the phagosomal pH in neutrophils"),
 ("5", "Reeves: K+ flux, not oxidants, releases granule proteases -- the "
       "model's non-oxidative alkaline arm.",
  "Reeves Segal killing activity of neutrophils proteases K+ flux Nature 2002"),
 ("5", "The critique of the protease hypothesis.",
  "role of oxidants versus proteases in neutrophil microbial killing controversy"),
 ("5", "HVCN1 / Hv1 proton channel and charge compensation for the oxidase.",
  "Hv1 HVCN1 proton channel charge compensation NADPH oxidase neutrophil"),
 ("5", "The V-ATPase in the phagosome.",
  "vacuolar ATPase phagosome acidification neutrophil macrophage"),
 ("5", "Neutrophil elastase and cathepsin G in bacterial killing.",
  "neutrophil elastase cathepsin G deficient mice impaired killing bacteria"),
 ("5", "Alpha-defensins HNP1-4: acid pH optimum and antimicrobial spectrum.",
  "human neutrophil alpha defensins HNP antimicrobial activity pH dependence"),
 ("5", "Bactericidal/permeability-increasing protein and lactoferrin.",
  "bactericidal permeability increasing protein lactoferrin neutrophil granule antimicrobial"),
 # ------------------------------------------------------------------- 6
 ("6. 카탈라제 도그마 · The catalase question",
  "THE key experiment: catalase-negative S. aureus is still fully virulent in "
  "CGD mice, which the catalase hypothesis forbids.",
  "Messina catalase negative Staphylococcus aureus virulent chronic granulomatous disease mice"),
 ("6", "Catalase-deficient organisms and CGD: the hypothesis stated and tested.",
  "why catalase positive organisms chronic granulomatous disease reconsidered"),
 ("6", "Aspergillus catalase mutants are not attenuated in CGD models either.",
  "Aspergillus fumigatus catalases virulence pulmonary aspergillosis"),
 ("6", "Bacterial catalase turnover number and kinetics.",
  "steady state kinetics of catalase hydrogen peroxide decomposition"),
 ("6", "Streptococcus pneumoniae hydrogen peroxide production by pyruvate "
       "oxidase -- the model's self-arming index anchor.",
  "Pericone factors contributing to hydrogen peroxide production Streptococcus pneumoniae"),
 ("6", "Hydrogen peroxide production by lactobacilli and other lactic acid bacteria.",
  "hydrogen peroxide production lactobacilli lactic acid bacteria quantitative"),
 ("6", "Bacterial peroxiredoxin AhpC and OxyR as the real H2O2 defence.",
  "AhpC alkyl hydroperoxide reductase OxyR primary scavenger hydrogen peroxide Escherichia coli"),
 ("6", "Superoxide dismutase in bacteria and why superoxide does not cross membranes.",
  "superoxide dismutase bacteria membrane impermeant superoxide anion"),
 ("6", "Granulibacter bethesdensis, a CGD-specific pathogen.",
  "Granulibacter bethesdensis chronic granulomatous disease lymphadenitis"),
 ("6", "Burkholderia cepacia complex infection in CGD.",
  "Burkholderia cepacia infection chronic granulomatous disease"),
 ("6", "Nocardia infection in CGD.",
  "Nocardia infection chronic granulomatous disease"),
 # ------------------------------------------------------------------- 7
 ("7. 아스페르길루스 · Aspergillus, the leading cause of death",
  "Invasive aspergillosis in CGD: presentation, species and outcome.",
  "invasive aspergillosis chronic granulomatous disease outcome species"),
 ("7", "Aspergillus nidulans as a CGD-specific and especially lethal pathogen.",
  "Aspergillus nidulans chronic granulomatous disease invasive infection"),
 ("7", "Neutrophils damage hyphae; macrophages handle conidia -- the model's "
       "central antifungal asymmetry.",
  "neutrophils hyphal damage Aspergillus fumigatus polymorphonuclear leukocytes"),
 ("7", "NADPH oxidase is required for antifungal defence independently of "
       "microbial killing per se.",
  "p47phox knockout mice invasive aspergillosis NADPH oxidase host defense"),
 ("7", "NETs and their contribution to control of Aspergillus.",
  "neutrophil extracellular traps Aspergillus fumigatus hyphae control"),
 ("7", "Galactomannan and beta-D-glucan for diagnosis in CGD.",
  "serum galactomannan assay invasive aspergillosis diagnostic accuracy"),
 ("7", "Mulch pneumonitis: the massive-conidial-inoculum syndrome.",
  "acute pneumonitis after mulch exposure chronic granulomatous disease Aspergillus"),
 # ------------------------------------------------------------------- 8
 ("8. 염증 · Why CGD inflames without an organism",
  "NADPH oxidase restrains the NLRP3 inflammasome: CGD cells overproduce IL-1beta.",
  "reactive oxygen species independent activation of the IL-1beta inflammasome in cells from patients with chronic granulomatous disease"),
 ("8", "Anakinra reverses granulomatous colitis in CGD -- the model's "
       "IL-1-blockade prediction.",
  "de Luca IL-1 receptor blockade restores autophagy and reduces inflammation chronic granulomatous disease"),
 ("8", "Defective efferocytosis of CGD neutrophils and failed resolution.",
  "impaired efferocytosis apoptotic neutrophils chronic granulomatous disease resolution inflammation"),
 ("8", "Phosphatidylserine externalisation is an oxidation-dependent event.",
  "phosphatidylserine externalization oxidation apoptotic cell recognition NADPH oxidase"),
 ("8", "LC3-associated phagocytosis requires NOX2.",
  "LC3 associated phagocytosis NOX2 required Rubicon"),
 ("8", "Nrf2 and the redox brake on inflammation in CGD.",
  "Nrf2 activation chronic granulomatous disease inflammation oxidative"),
 ("8", "IL-17 / Th17 skewing in CGD.",
  "increased Th17 responses in patients with chronic granulomatous disease"),
 ("8", "Granulomatous colitis in CGD: prevalence, endoscopy and histology.",
  "gastrointestinal involvement chronic granulomatous disease colitis patients"),
 ("8", "Gastric outlet and genitourinary obstruction from CGD granulomas.",
  "gastric outlet obstruction bladder granuloma chronic granulomatous disease"),
 ("8", "Corticosteroids for CGD granulomatous complications.",
  "corticosteroids treatment of obstructive gastrointestinal granulomata chronic granulomatous disease"),
 ("8", "The tryptophan/kynurenine hypothesis: the original claim.",
  "Romani defective tryptophan catabolism underlies inflammation chronic granulomatous disease"),
 ("8", "And its failure to replicate.",
  "De Ravin tryptophan catabolism kynurenine chronic granulomatous disease normal indoleamine dioxygenase"),
 # ------------------------------------------------------------------- 9
 ("9. 예방요법 · Prophylaxis: the three trials the model predicts",
  "Trimethoprim-sulfamethoxazole prophylaxis reduces bacterial infection.",
  "trimethoprim-sulfamethoxazole prophylaxis in the management of chronic granulomatous disease"),
 ("9", "Itraconazole prophylaxis prevents fungal infection: the crossover trial.",
  "Gallin itraconazole to prevent fungal infections in chronic granulomatous disease"),
 ("9", "Interferon gamma-1b reduces serious infection by two thirds.",
  "controlled trial of interferon gamma to prevent infection in chronic granulomatous disease"),
 ("9", "Long-term follow-up and the durability of interferon gamma prophylaxis.",
  "long term interferon gamma therapy chronic granulomatous disease follow up"),
 ("9", "The unresolved question of interferon gamma's mechanism in CGD.",
  "mechanism of action interferon gamma chronic granulomatous disease superoxide restoration"),
 ("9", "Trimethoprim-sulfamethoxazole pharmacokinetics.",
  "trimethoprim sulfamethoxazole pharmacokinetics children clearance volume distribution"),
 ("9", "Itraconazole pharmacokinetics and the hydroxy metabolite.",
  "itraconazole hydroxyitraconazole pharmacokinetics steady state trough"),
 ("9", "Voriconazole nonlinear pharmacokinetics and CYP2C19.",
  "voriconazole nonlinear pharmacokinetics CYP2C19 saturable metabolism"),
 ("9", "Azole hepatotoxicity and drug interactions.",
  "azole antifungal hepatotoxicity drug interaction CYP3A4"),
 ("9", "Interferon gamma-1b pharmacokinetics after subcutaneous dosing.",
  "interferon gamma-1b pharmacokinetics subcutaneous administration"),
 ("9", "Anakinra pharmacokinetics.",
  "anakinra pharmacokinetics subcutaneous interleukin-1 receptor antagonist"),
 # ------------------------------------------------------------------- 10
 ("10. 이식과 유전자치료 · Putting the electron flux back",
  "Reduced-intensity conditioning HSCT in CGD: the prospective multicentre study.",
  "Gungor reduced-intensity conditioning haematopoietic cell transplantation chronic granulomatous disease prospective"),
 ("10", "Myeloid chimerism thresholds and resolution of CGD manifestations.",
  "myeloid chimerism after transplantation chronic granulomatous disease outcome"),
 ("10", "Mixed chimerism and the fraction of oxidase-normal neutrophils needed.",
  "stable mixed chimerism sufficient to cure chronic granulomatous disease"),
 ("10", "Long-term outcome of HSCT for CGD in large registries.",
  "outcome allogeneic hematopoietic stem cell transplantation chronic granulomatous disease registry"),
 ("10", "Lentiviral gene therapy for X-linked CGD: the multicentre trial.",
  "Kohn lentiviral gene therapy for X-linked chronic granulomatous disease"),
 ("10", "The earlier gamma-retroviral trials and the MECOM/EVI1 clonal expansions.",
  "correction of X-linked chronic granulomatous disease by gene therapy clonal expansion MDS1 EVI1"),
 ("10", "Gene editing approaches to CGD.",
  "CRISPR gene editing correction chronic granulomatous disease hematopoietic stem cells"),
 ("10", "Granulocyte transfusion in refractory CGD infection.",
  "granulocyte transfusion refractory infection chronic granulomatous disease"),
 # ------------------------------------------------------------------- 11
 ("11. 호중구 동태 · Granulopoiesis and neutrophil kinetics",
  "Neutrophil kinetics in man: marrow pools, blood half-life and transit times.",
  "Dancey neutrophil kinetics in man marrow pool blood transit"),
 ("11", "G-CSF regulation of granulopoiesis and receptor-mediated clearance.",
  "target mediated disposition granulocyte colony stimulating factor neutrophil model"),
 ("11", "Quantitative models of granulopoiesis for PK/PD use.",
  "semi-mechanistic model myelosuppression neutrophil Friberg"),
 # ------------------------------------------------------------------- 12
 ("12. QSP 방법론 · Modelling methodology",
  "mrgsolve as the ODE engine this model targets.",
  "mrgsolve simulation pharmacometrics R package"),
 ("12", "Quantitative systems pharmacology: scope and good practice.",
  "quantitative systems pharmacology models current status future directions"),
 ("12", "Virtual population generation for QSP models.",
  "virtual population generation quantitative systems pharmacology plausible patients"),
 ("12", "Critical inoculum and threshold behaviour in host-pathogen models.",
  "mathematical model bacterial infection innate immune response bistability clearance"),
 ("12", "Within-host models of neutrophil-bacteria interaction.",
  "within host mathematical model neutrophil bacteria dynamics phagocytosis"),
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
    out.append("# 만성 육아종병 (Chronic Granulomatous Disease) — 참고문헌\n")
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
    open("cgd_references.md", "w").write("".join(out))
    print("wrote cgd_references.md with %d resolved references" % n_tot)


if __name__ == "__main__":
    main()
