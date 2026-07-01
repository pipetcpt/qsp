# Mucopolysaccharidosis Type I (MPS I, 뮤코다당증 제1형) — QSP Disease Model

> **Mucopolysaccharidosis type I (MPS I)** is a lysosomal storage disease
> caused by **loss-of-function mutations in IDUA** (4p16.3), which encodes
> **alpha-L-iduronidase**. Enzyme deficiency blocks the stepwise catabolism
> of **dermatan sulfate** and **heparan sulfate**, causing progressive
> multi-organ glycosaminoglycan (GAG) accumulation across a clinical
> spectrum from severe **Hurler syndrome** (CNS-involved, infantile onset)
> through intermediate **Hurler-Scheie** to attenuated **Scheie syndrome**
> (no CNS involvement, adult presentation). Manifestations include
> dysostosis multiplex, cardiac valve thickening/cardiomyopathy,
> upper-airway obstruction/OSA, restrictive lung disease, hepatosplenomegaly,
> corneal clouding, and — in severe Hurler syndrome — progressive
> neurodegeneration behind an intact blood-brain barrier.

The model captures: **IDUA loss-of-function** → alpha-L-iduronidase
deficiency (<1% normal activity) → lysosomal **dermatan/heparan sulfate**
accumulation → secondary lysosomal/autophagy dysfunction and
TLR4/NF-κB-mediated inflammation → organ-specific GAG burden driving
**dysostosis multiplex**, **valvular/myocardial disease**, **upper-airway
obstruction**, **hepatosplenomegaly**, **corneal clouding**, and (severe
Hurler) **CNS neurodegeneration**. **Laronidase** (recombinant human
alpha-L-iduronidase, IV ERT) restores enzyme activity in well-perfused
visceral tissue via mannose-6-phosphate-receptor-mediated uptake but cannot
cross the blood-brain barrier and penetrates avascular cartilage/cornea
poorly. **Hematopoietic stem cell transplantation (HSCT)** achieves donor
engraftment that cross-corrects visceral tissue AND — uniquely — repopulates
CNS microglia with enzyme-competent donor monocytes, the only clinically
validated route to halt neurodegeneration, provided transplantation occurs
early. Investigational modalities modeled: lentiviral HSC gene therapy
(OTL-203), AAV9 CNS-directed gene therapy, and oral substrate-reduction
therapy (genistein).

---

## Deliverables

| File | Purpose |
|------|---------|
| `mps1_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **11 clusters, 133 nodes** |
| `mps1_mrgsolve_model.R` | **21-ODE** mrgsolve model (8 PK/intervention + 13 disease/clinical) with 10 scenarios |
| `mps1_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `mps1_references.md` | **30** PubMed citations grouped by section |

---

## Mechanistic Map — Cluster Index

1. **Genetics & etiology** — IDUA 4p16.3 loss-of-function, autosomal recessive, null/null vs. missense genotypes, genotype-phenotype correlation, newborn screening
2. **GAG substrate & catabolic pathway** — dermatan/heparan sulfate, normal stepwise exo-glycosidase degradation, M6P tagging/M6PR trafficking
3. **Lysosomal storage & secondary pathology** — lysosomal swelling, autophagy dysfunction, secondary ganglioside accumulation, TLR4/NF-κB inflammation
4. **Skeletal — dysostosis multiplex** — growth-plate chondrocyte storage, kyphoscoliosis, joint contractures, odontoid hypoplasia/cervical stenosis, short stature
5. **Cardiovascular manifestations** — valve leaflet GAG infiltration, myocardial infiltration/cardiomyopathy, coronary intimal narrowing, heart failure
6. **Respiratory / upper airway** — adenotonsillar hypertrophy, tracheobronchial narrowing, OSA, restrictive lung physiology
7. **CNS / neurological** — CNS GAG accumulation behind the BBB, hydrocephalus, neurodegeneration, developmental quotient, cervical myelopathy
8. **Ophthalmologic / hepatosplenic / other** — corneal clouding, hepatosplenomegaly, hearing loss, coarse facial features
9. **Laronidase ERT PK/PD** — plasma vs. tissue-retained enzyme, M6PR-mediated uptake, ADA formation, organ-specific penetration (good visceral, poor cartilage, none CNS)
10. **HSCT / gene therapy / investigational SRT** — myeloablative conditioning, donor chimerism, cross-correction, CNS microglial replacement, lentiviral/AAV9 gene therapy, genistein
11. **Clinical endpoints / biomarkers** — urinary GAG, 6MWT, FVC, joint ROM, liver/spleen volume, echocardiographic valve thickness, corneal score, DQ, height Z-score, survival

---

## mrgsolve Model

### ODE Compartments (21)
**PK / intervention (8):** LARO_CENT, LARO_TISSUE (laronidase plasma + slowly-decaying tissue-retained active enzyme pool); ADA (anti-laronidase antibody); GEN_GUT, GEN_CENT (genistein); CHIMERISM, CNS_ENZ_ACCESS (HSCT engraftment / CNS microglial replacement); LVGT_SYS_ACCESS (lentiviral HSC gene therapy)

**Disease / PD / clinical (13):** GAG_SYS, GAG_CNS, GAG_CART (organ-specific GAG burden pools); UGAG, LIVSPLEEN, VALVE, FVC, AHI, JOINTROM, CORNEA, DQ, HEIGHTZ, HAZARD (+ derived SURVIVAL)

### Treatment Scenarios (10)
1. **Untreated natural history** — severe Hurler phenotype reference
2. **Laronidase ERT monotherapy** — attenuated Hurler-Scheie/Scheie, 0.58 mg/kg IV weekly
3. **HSCT alone, early transplant (9 mo)** — severe Hurler, window-of-opportunity favorable
4. **HSCT alone, delayed transplant (30 mo)** — severe Hurler, accrued CNS injury before engraftment
5. **ERT bridging (~15 wk) then HSCT** — standard-of-care combination, ERT stopped at engraftment
6. **ERT, high anti-drug antibody titer** — immunogenicity sensitivity
7. **ERT, poor adherence (60% of infusions)** — real-world adherence sensitivity
8. **ERT + genistein (investigational SRT)** — substrate-reduction adjunct
9. **ERT + investigational AAV9 CNS gene therapy** — CNS-targeted adjunct without HSCT
10. **HSCT + long-term low-dose ERT** — post-transplant maintenance for residual visceral/joint disease

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| Wraith 2004 J Pediatr (PMID 15126990, pivotal RCT) | uGAG, liver volume, FVC, AHI, shoulder flexion | week-26 ERT treatment effect |
| Clarke 2009 Pediatrics (PMID 19117887, 3.5-yr extension) | sustained biomarker/organ response | long-term ERT plateau |
| Sifuentes 2007 Mol Genet Metab (PMID 17011223, 6-yr follow-up) | uGAG -76%, growth, sleep apnea | long-horizon ERT ceiling |
| Peters 1998 Blood / Boelens 2013 Blood (PMID 9516162 / 23493783) | donor engraftment, chimerism | HSCT engraftment kinetics |
| Aldenhoven 2015 Blood / Eisengart 2018 Genet Med (PMID 25624320 / 29517765) | developmental outcome by transplant age | window-of-opportunity DQ relationship |
| Gentner 2021 NEJM (PMID 34788506) | lentiviral HSC gene therapy | illustrative investigational parameters |

---

## Shiny App — 8 Tabs

1. **Patient & Overview** — phenotype/regimen sidebar + mechanistic-map schematic
2. **Drug PK** — laronidase plasma vs. tissue-retained enzyme; ADA/genistein exposure
3. **Enzyme access / GAG burden** — systemic/CNS/cartilage enzyme-access indices and GAG pools
4. **Clinical endpoints** — uGAG, DQ, liver/spleen index, survival; full trajectory panel
5. **Scenario comparison** — all 10 regimens overlaid (uGAG) + endpoint table
6. **Biomarkers** — cardiac valve/LV-mass index, corneal clouding, height Z-score
7. **Safety / Survival** — cumulative mortality hazard and survival probability
8. **References** — key trial citations

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg mps1_qsp_model.dot -o mps1_qsp_model.svg
dot -Tpng -Gdpi=150 mps1_qsp_model.dot -o mps1_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("mps1_mrgsolve_model.R")           # builds `mps1_mod` + `scenarios`
res <- run_scenario("5_ERT_Bridging_then_HSCT_SevereHurler", scenarios[["5_ERT_Bridging_then_HSCT_SevereHurler"]])
plot(res$time/24/365, res$DQ, type = "l")

# Launch the dashboard
shiny::runApp("mps1_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| 뮤코다당증 제1형 | Mucopolysaccharidosis type I (MPS I) |
| 알파-L-이두로니다아제 | Alpha-L-iduronidase (IDUA) |
| 글리코사미노글리칸 | Glycosaminoglycan (GAG) |
| 더마탄황산 / 헤파란황산 | Dermatan sulfate / Heparan sulfate |
| 헐러 증후군 (중증형) | Hurler syndrome (severe phenotype) |
| 샤이에 증후군 (경증형) | Scheie syndrome (attenuated phenotype) |
| 골격이형성(다발성뼈형성이상) | Dysostosis multiplex |
| 조혈모세포이식 | Hematopoietic stem cell transplantation (HSCT) |
| 효소대체요법 | Enzyme replacement therapy (ERT) |
| 발달지수 | Developmental quotient (DQ) |

---

*Built by Claude Code Routine on 2026-07-01 as part of the QSP Disease Model
Library. See root [README.md](../README.md) for the full model gallery.*
