# X-Linked Hypophosphatemia (XLH, 저인산혈증성 구루병) — QSP Disease Model

> **X-linked hypophosphatemia (XLH)** is the most common heritable form of
> rickets/osteomalacia, caused by **loss-of-function mutations in PHEX**
> (Xp22.11). PHEX loss leads osteocytes to chronically **overproduce FGF23**,
> which signals through the renal **FGFR1c/αKlotho** receptor complex to
> internalize the **NPT2a/NPT2c** sodium-phosphate cotransporters and
> suppress **1α-hydroxylase** — producing renal phosphate wasting and
> inappropriately low calcitriol despite hypophosphatemia. The chronic
> phosphate deficit impairs growth-plate and osteoid mineralization,
> causing pediatric rickets, adult osteomalacia, leg bowing, enthesopathy,
> myopathy, and spontaneous dental abscesses.

The model captures: **PHEX loss-of-function** → osteocyte **FGF23**
overproduction (unopposed ASARM-peptide/MEPE dysregulation) → circulating
FGF23 engages **FGFR1c-αKlotho** in the renal proximal tubule → **NPT2a/c
internalization** → **TmP/GFR ↓↓** → renal phosphate wasting + suppressed
**CYP27B1** (1α-hydroxylase) → chronic **hypophosphatemia** and
inappropriately low **1,25(OH)2D** → defective mineralization → **rickets**
(pediatric) / **osteomalacia** (adult), leg bowing, enthesopathy, early OA,
myopathy, dental abscesses. **Burosumab** is a fully human anti-FGF23 IgG1
monoclonal antibody that neutralizes circulating FGF23 directly, restoring
both the renal phosphate-reabsorption and calcitriol-synthesis arms;
**conventional therapy** (oral phosphate salts + active vitamin D) bypasses
the FGF23 axis via direct substrate replacement, correcting biochemistry
only transiently and carrying nephrocalcinosis/hyperparathyroidism risk
with chronic use.

---

## Deliverables

| File | Purpose |
|------|---------|
| `xlh_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **15 clusters, 121 nodes** |
| `xlh_mrgsolve_model.R` | **19-ODE** mrgsolve model (7 PK + 12 disease/clinical) with 10 scenarios |
| `xlh_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `xlh_references.md` | **33** PubMed citations grouped by section |

---

## Mechanistic Map — Cluster Index

1. **Genetics & etiology** — PHEX Xp22.11 loss-of-function, X-linked dominant, ~20-30% de novo, ADHR (FGF23 GOF)/ARHR1 (DMP1)/ARHR2 (ENPP1)/TIO paraneoplastic spectrum
2. **PHEX-FGF23 osteocyte axis** — MEPE/ASARM peptide, PC2/furin proprotein convertase, intact vs. cleaved FGF23, local mineralization feedback disruption
3. **FGF23-FGFR1c/αKlotho signaling** — ternary receptor complex, FRS2-RAS-RAF-MEK-ERK, EGR1, burosumab's direct FGF23-neutralization target
4. **Renal phosphate handling** — NPT2a (SLC34A1)/NPT2c (SLC34A3) internalization, TmP/GFR, urinary phosphate wasting
5. **Vitamin D / calcitriol metabolism** — CYP27B1 suppression, CYP24A1 induction, inappropriately low 1,25(OH)2D, reduced intestinal Ca/Pi absorption
6. **PTH-calcium axis** — secondary → tertiary hyperparathyroidism, PTH-FGF23 amplifying feedback, parathyroidectomy
7. **Bone mineralization pathology** — osteoid accumulation, growth-plate widening, rickets/osteomalacia, craniosynostosis, Chiari I, pseudofractures
8. **Skeletal / orthopedic** — genu varum/valgum, short stature, enthesopathy, early-onset OA, spinal stenosis, corrective osteotomy
9. **Dental manifestations** — interglobular dentin defects, enlarged pulp chambers, spontaneous abscesses, periodontal disease
10. **Muscular / functional** — proximal myopathy, gait abnormality, chronic pain, fatigue, QoL, employment impact
11. **Drug PK — burosumab** — SC 2-compartment mAb PK, FcRn recycling (~19-day t½), immunogenicity
12. **Drug PK — conventional therapy** — oral phosphate salts (transient spike), oral calcitriol/alfacalcidol
13. **Drug PD** — FGF23 neutralization → NPT2/TmP-GFR/calcitriol rescue (burosumab, sustained) vs. direct substrate boost (conventional, transient) and overcorrection risk
14. **Clinical endpoints** — Rickets Severity Score (RSS), RGI-C, AGV, height Z-score, 6MWT, WOMAC, bone-specific ALP (BSAP)
15. **Safety / adverse events** — injection-site reactions, restless legs (burosumab); nephrocalcinosis, hypercalciuria, hyperphosphatemia overcorrection (conventional therapy)

---

## mrgsolve Model

### ODE Compartments (19)
**PK (7):** BURO_DEPOT, BURO_CENT, BURO_PERIPH (burosumab); PHOSORAL_GUT, PHOSORAL_SIG (oral phosphate); CALC_GUT, CALC_CENT (oral calcitriol)

**Disease / PD / clinical (12):** NPT2, TMPGFR, PHOS, CALCITRIOL, PTH, BSAP,
RSS, HEIGHTZ_XLH, SIXMWT, WOMAC, UCACR, NEPHROCALC (+ derived FGF23_NEUT, AGV_CALC_XLH)

### Treatment Scenarios (10)
1. **Untreated** — natural history reference (PHOS ≈ 2.2 mg/dL, RSS ≈ 5.5)
2. **Conventional therapy (ped)** — oral phosphate + calcitriol, standard dosing
3. **Burosumab 0.8 mg/kg SC Q2W** — approved pediatric starting dose
4. **Burosumab 2.0 mg/kg SC Q2W** — maximum approved pediatric dose
5. **Burosumab 1.0 mg/kg SC Q4W** — approved adult regimen (AXLES1)
6. **Switch: conventional → burosumab** — CL303-style crossover (modeled as post-switch burosumab phase)
7. **Conventional therapy, poor GI adherence (60%)** — real-world tolerability sensitivity
8. **Burosumab, supratherapeutic overcorrection** — hyperphosphatemia risk exploration
9. **Conventional therapy, long-term** — chronic tertiary hyperparathyroidism / nephrocalcinosis risk accrual
10. **Burosumab adult, 60% adherence** — injection-adherence sensitivity

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| Carpenter 2018 NEJM (PMID 29791829, phase 2) | Burosumab dose-finding, TmP/GFR & serum Pi rescue | Q2W/Q4W dose-response |
| Imel 2019 Lancet (PMID 31104833, CL303 phase 3) | RSS/RGI-C & growth vs. conventional therapy | burosumab superiority |
| Insogna 2018 JBMR (PMID 29947083, AXLES1 adult) | 6MWT/WOMAC, week-24 primary analysis | adult functional/pain gains |
| Carpenter 2014 JCI (PMID 24569459, KRN23 FIH) | Burosumab PK (t½ ~19 d) | 2-cpt mAb PK parameters |
| Imel 2024 JBMR (PMID 39151033) | Nephrocalcinosis prevalence, conventional vs. burosumab | 22% (ped)/38% (adult) baseline risk |
| Haffner 2019 Nat Rev Nephrol (PMID 31068690) | Diagnosis/monitoring guideline | biochemical monitoring targets |

---

## Shiny App — 8 Tabs

1. **Patient & Overview** — covariate sidebar + mechanistic-map schematic
2. **Drug PK** — burosumab and oral phosphate/calcitriol concentration-time
3. **Pathway PD** — FGF23 neutralization vs. NPT2a/c rescue; TmP/GFR & serum phosphate
4. **Clinical endpoints** — RSS, height Z-score (ped), 6MWT/WOMAC (adult)
5. **Scenario comparison** — all 10 regimens overlaid + endpoint table
6. **Biomarkers** — 1,25(OH)2D, PTH, bone-specific ALP (BSAP)
7. **Safety** — urine Ca/Cr ratio and nephrocalcinosis risk index
8. **References** — key trial citations

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg xlh_qsp_model.dot -o xlh_qsp_model.svg
dot -Tpng -Gdpi=150 xlh_qsp_model.dot -o xlh_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("xlh_mrgsolve_model.R")           # builds `xlh_mod` + `scenarios`
res <- run_scenario("3_Burosumab_Ped_0p8mgkg_Q2W", scenarios[["3_Burosumab_Ped_0p8mgkg_Q2W"]])
plot(res$time/24, res$RSS, type = "l")

# Launch the dashboard
shiny::runApp("xlh_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| X-연관 저인산혈증 | X-linked hypophosphatemia (XLH) |
| 인산조절유전자(내페스티다아제 상동) | PHEX (phosphate-regulating endopeptidase homolog, X-linked) |
| 섬유아세포성장인자23 | Fibroblast growth factor 23 (FGF23) |
| 나트륨-인산 공동수송체 | Sodium-phosphate cotransporter (NPT2a/NPT2c) |
| 신세관 최대재흡수/사구체여과율 | Tubular maximum phosphate reabsorption / GFR (TmP/GFR) |
| 구루병(소아) / 골연화증(성인) | Rickets (pediatric) / Osteomalacia (adult) |
| 구루병 중증도 점수 | Rickets Severity Score (RSS) |
| 신석회화증 | Nephrocalcinosis |
| 이차성/삼차성 부갑상선기능항진증 | Secondary/tertiary hyperparathyroidism |

---

*Built by Claude Code Routine on 2026-07-01 as part of the QSP Disease Model
Library. See root [README.md](../README.md) for the full model gallery.*
