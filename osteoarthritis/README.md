# Osteoarthritis (OA) QSP Model

[![Disease](https://img.shields.io/badge/Disease-Osteoarthritis-orange)]()
[![Category](https://img.shields.io/badge/Category-Musculoskeletal%20%7C%20Chronic-blue)]()
[![Model](https://img.shields.io/badge/mrgsolve-20%20ODE%20compartments-green)]()
[![Shiny](https://img.shields.io/badge/Shiny-8%20tabs-purple)]()
[![References](https://img.shields.io/badge/References-60%20PubMed-red)]()

---

## Disease Overview

**Osteoarthritis (OA)** is the most prevalent musculoskeletal disease worldwide, affecting ~500 million people globally (Vos et al., *Lancet*, 2012). It is a whole-joint disease characterised by:

- Progressive **articular cartilage degradation** driven by an imbalance between anabolic (collagen II / aggrecan synthesis) and catabolic (MMP-13, ADAMTS-5) processes
- **Subchondral bone remodelling** including sclerosis, osteophyte formation, and bone marrow lesions
- **Synovitis** with macrophage M1 activation, IL-1β/TNF-α production, and fibroblast-like synoviocyte (FLS) amplification
- **Chronic pain** mediated by PGE2, NGF/TrkA, substance P, TRPV1, and central sensitisation
- Cellular senescence, oxidative stress, and adipokine (leptin)-driven joint inflammation

Despite being the leading cause of disability in adults, **no disease-modifying OA drug (DMOAD) is currently approved**. Current treatment is symptom-focused (NSAIDs, IA corticosteroids, IA hyaluronic acid) with total knee replacement (TKR) as the endpoint for severe disease.

---

## Mechanistic Map

[![OA QSP Mechanistic Map](oa_qsp_model.png)](oa_qsp_model.svg)

> Click on the image to open the full SVG. ← navigates to individual cluster details.

### Cluster Summary (12 clusters, 130+ nodes)

| Cluster | Key Components |
|---------|---------------|
| 1. Risk Factors | Age, BMI/adipokines, female sex (estrogen loss), joint injury, GDF5/ALDH1A2 SNPs, KL grade |
| 2. Chondrocyte Biology | NF-κB, MAPK-p38/ERK, SOX-9, RUNX2, apoptosis, autophagy, mTORC1, SIRT1, p53/p21 |
| 3. Cartilage ECM | Collagen II, aggrecan, fibronectin-f, lubricin (PRG4), COMP, JSW, cartilage volume |
| 4. Inflammatory Mediators | IL-1β, TNF-α, IL-6, IL-8, IL-17A, HMGB1, S100A8/A9, TGF-β1, PGE2, LTB4, NO |
| 5. Catabolic Proteases | MMP-1/3/9/13/14, ADAM-17, ADAMTS-4/5, cathepsin K/B/D, HTRA1 |
| 6. Anabolic & Protective | IGF-1, BMP-7, FGF-18 (sprifermin target), FGFR3, TIMP-1/2/3, IL-10, IL-4 |
| 7. Subchondral Bone | RANKL/OPG/RANK, osteoclasts, osteoblasts, sclerosis, osteophytes, Wnt/DKK1/sclerostin |
| 8. Synovitis | FLS, M1/M2 macrophages, mast cells, neutrophils, Th17, Treg, effusion |
| 9. Pain Signalling | COX-2/PGE2/EP1-4, NGF/TrkA, SubstP, CGRP, TRPV1, ASIC3, DRG, spinal cord |
| 10. Drug PK/PD | NSAIDs (celecoxib), IA-CS (triamcinolone), IA-HA, sprifermin, tanezumab, lorecivivint |
| 11. Biomarkers | uCTX-II, serum COMP, serum HA, CS-846, hs-CRP, MMP-3 |
| 12. Clinical Endpoints | VAS pain, KOOS, WOMAC, JSW, TKR probability, OARSI grade |

---

## mrgsolve ODE Model

**File:** [`oa_mrgsolve_model.R`](oa_mrgsolve_model.R)

### Model Architecture

| Module | Compartments | Key Equations |
|--------|-------------|---------------|
| NSAID PK | `A_NSAID_gut`, `A_NSAID_plasma`, `A_NSAID_joint` | 2-compartment oral; Kp_joint = 0.4; IC50_COX2 = 0.042 μg/mL |
| IA Corticosteroid | `A_IACS_joint`, `A_IACS_plasma` | Slow joint absorption k = 0.008/h; GR-mediated NF-κB inhibition |
| IA Hyaluronic Acid | `A_HA_joint` | t½ ~17h; viscosupplementation → pain↓ + IL-1β↓ |
| Sprifermin (FGF-18) | `A_Sprif_joint` | t½ ~16min; FGFR3 → ColII synthesis ↑ (FORWARD RCT) |
| Tanezumab | `A_Tanz_depot`, `A_Tanz_plasma` | Anti-NGF mAb; SC F = 73%, t½ = 23d |
| Inflammatory | `IL1b`, `TNFa` | Synovitis-driven; NF-κB feedback |
| Catabolic Enzymes | `MMP13`, `ADAM5` | IL-1β/TNF-α Emax induction; TIMP-mediated inhibition |
| Cartilage ECM | `ColII`, `Aggrecan` | Synthesis (chondrocyte-dependent) vs. MMP/ADAMTS-mediated degradation |
| Chondrocytes | `Chondro` | Proliferation vs. IL-1β-driven apoptosis |
| Synovitis | `Synovitis` | ECM-loss feedback + IL-1β amplification loop |
| Bone | `OC_act`, `OB_act` | RANKL/OPG balance; IL-1β → RANKL ↑ |
| JSW | `JSW` | Degradation-driven loss; sprifermin attenuates |
| Pain | `PGE2_jt`, `VASPain` | PGE2 + synovitis + structural components |
| Biomarkers | `uCTXII`, `COMP_s` | Released from ColII degradation / cartilage damage |

### Treatment Scenarios (8 scenarios over 2 years)

| # | Scenario | Key Mechanism | Clinical Reference |
|---|----------|--------------|-------------------|
| 1 | Natural history | Disease progression without treatment | Historical controls |
| 2 | Celecoxib 200mg BID | COX-2 inhibition → PGE2 ↓ | CONDOR trial (Goldstein 2010) |
| 3 | IA Triamcinolone 40mg Q3mo | GR → NF-κB inhibition → IL-1β ↓ | Raynauld 2003 Arthritis Rheum |
| 4 | IA Hyaluronic Acid 16mg×3 | Viscosupplementation | Bellamy 2006 Cochrane |
| 5 | Sprifermin 30μg Q12w | FGF-18/FGFR3 → ColII synthesis ↑ | FORWARD RCT (Hochberg 2019 JAMA) |
| 6 | Tanezumab 2.5mg SC Q8w | Anti-NGF → TrkA block → pain ↓ | Lane 2010 NEJM; Schnitzer 2019 |
| 7 | Celecoxib + IA Triamcinolone | Combined symptomatic | Expert consensus |
| 8 | Sprifermin + Celecoxib | Structural + symptomatic | Phase-III rationale |

---

## Shiny Dashboard

**File:** [`oa_shiny_app.R`](oa_shiny_app.R)

### Tab Structure (8 tabs)

| Tab | Content |
|-----|---------|
| 1. Patient Profile | KL grade, age selection → disease state table + overview |
| 2. Drug PK | NSAID plasma/joint Cmax, COX-2 inhibition %, IA-CS, HA, sprifermin, tanezumab PK |
| 3. Inflammatory Markers | IL-1β, TNF-α, synovitis score, MMP-13, ADAMTS-5, PGE2 time-courses |
| 4. Cartilage & Bone | Collagen II, aggrecan, chondrocyte pool, JSW (with TKR threshold), OC/OB |
| 5. Pain & Function | VAS pain, estimated KOOS, pain component breakdown (PGE2 / synovitis / structural) |
| 6. Scenario Comparison | 6 pre-defined scenarios; VAS, JSW, KOOS plots + 1-year outcome table |
| 7. Biomarkers | uCTX-II, serum COMP with clinical reference range table |
| 8. Sensitivity Analysis | Tornado plots; parameter sweep (±30%) for any outcome at 1 year |

### How to Run

```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","ggplot2",
                   "tidyr","plotly","DT","scales","purrr"))
shiny::runApp("oa_shiny_app.R")
```

---

## References

Full annotated bibliography in [`oa_references.md`](oa_references.md) — **60 PubMed citations** organised across 12 sections:

1. Disease Overview & Pathophysiology
2. Chondrocyte Biology & ECM
3. Inflammatory Cytokines & Signaling
4. Matrix Metalloproteinases & ADAMTS
5. Subchondral Bone Remodelling
6. Pain Mechanisms
7. Drug Pharmacology (NSAIDs, Corticosteroids, HA)
8. DMOADs & Emerging Therapies
9. Biomarkers of Cartilage Degradation
10. Clinical Outcomes & Trial Design
11. Epidemiology & Risk Factors
12. QSP / Computational Modelling

---

## File List

| File | Description |
|------|-------------|
| [`oa_qsp_model.dot`](oa_qsp_model.dot) | Graphviz source — 130+ nodes, 12 clusters |
| [`oa_qsp_model.svg`](oa_qsp_model.svg) | Scalable mechanistic map (vector) |
| [`oa_qsp_model.png`](oa_qsp_model.png) | Raster mechanistic map (150 dpi) |
| [`oa_mrgsolve_model.R`](oa_mrgsolve_model.R) | mrgsolve ODE model — 20 compartments, 8 scenarios |
| [`oa_shiny_app.R`](oa_shiny_app.R) | Shiny dashboard — 8 interactive tabs |
| [`oa_references.md`](oa_references.md) | 60 PubMed references |
| [`README.md`](README.md) | This file |

---

*Model built: 2026-06-23 | QSP Disease Model Library*
