# Atopic Dermatitis (AD) — Quantitative Systems Pharmacology Model

[![AD QSP Map](ad_qsp_model.png)](ad_qsp_model.svg)

> **Click the image above** to open the full-resolution SVG mechanistic map.

---

## Disease Overview

**Atopic dermatitis (AD)** is the most prevalent chronic inflammatory skin disease, affecting ~15–20% of children and ~5–10% of adults worldwide. It is characterized by:

- **Epidermal barrier dysfunction** (FLG mutation → TEWL↑, pH↑, allergen penetration)
- **Th2-skewed immune activation** (IL-4, IL-13, IL-31, IgE, TARC/CCL17)
- **Chronic pruritus** (itch-scratch cycle via IL-31Ra, TRPV1, DRG neurons)
- **Comorbidities** (atopic march: asthma, allergic rhinitis, food allergy; anxiety/depression)

The pathology involves innate alarmins (TSLP, IL-33, IL-25) activating ILC2 and dendritic cells, which drive Th2 differentiation via GATA-3/STAT6 signaling. IL-4/IL-13 co-signal through the shared IL-4Rα chain — the therapeutic target of dupilumab.

---

## Model Files

| File | Description |
|------|-------------|
| [`ad_qsp_model.dot`](ad_qsp_model.dot) | Graphviz mechanistic map source (130+ nodes, 13 clusters) |
| [`ad_qsp_model.svg`](ad_qsp_model.svg) | Rendered vector diagram (interactive) |
| [`ad_qsp_model.png`](ad_qsp_model.png) | Rendered raster image (150 dpi) |
| [`ad_mrgsolve_model.R`](ad_mrgsolve_model.R) | mrgsolve ODE model (R) — PK/PD/disease dynamics |
| [`ad_shiny_app.R`](ad_shiny_app.R) | Shiny interactive dashboard (6 tabs) |
| [`ad_references.md`](ad_references.md) | 52 PubMed references, categorized |

---

## Mechanistic Map Structure (13 Clusters)

| Cluster | Contents | Nodes |
|---------|----------|-------|
| 1. Environmental & Genetic Triggers | Allergens, microbiome, FLG/SPINK5 mutations, UV | 8 |
| 2. Epidermal Barrier | FLG, ceramides, CLDN1, TEWL, pH, KLK5 | 11 |
| 3. Innate Immunity & Alarmins | TSLP, IL-25, IL-33, ILC2, mast cells, DC | 14 |
| 4. Th2 Adaptive Immunity | IL-4, IL-13, IL-31, STAT6, JAK1/2, IgE, TARC | 20 |
| 5. Eosinophil / Mast Cell Effectors | IL-5, ECP, MBP, PGD2, leukotrienes | 10 |
| 6. Th1 / Th17 / Th22 Counter-Regulation | IFN-γ, IL-17A, IL-22, Treg, TGF-β | 10 |
| 7. Pruritus Neural Pathway | IL-31Ra, TRPV1, TRPA1, DRG neuron, NRS | 13 |
| 8. Keratinocyte Inflammatory Responses | TSLP-KC, CCL27, CXCL10, IL-36 | 10 |
| 9. Biologic PK | Dupilumab, tralokinumab, lebrikizumab, nemolizumab | 10 |
| 10. JAK Inhibitor PK | Upadacitinib, abrocitinib, baricitinib, deucravacitinib | 10 |
| 11. Drug PD Effects | RO, pSTAT6 inhibition, TARC/eosinophil reduction | 11 |
| 12. Clinical Endpoints | EASI, IGA, SCORAD, NRS, DLQI, POEM | 14 |
| 13. Comorbidities | Atopic march, anxiety/depression, CardioRisk | 10 |

**Total: 130+ nodes, 200+ directed edges**

---

## mrgsolve ODE Model (20 Compartments)

### Drug PK Compartments

| Drug | Model | Key Parameters |
|------|-------|----------------|
| Dupilumab | 2-cmt SC mAb | CL=0.21 L/day, V1=3.5 L, t½≈21d, F_SC=64% |
| Upadacitinib | 1-cmt oral | CL/F=38 L/day, Vd=166 L, t½≈8h |
| Nemolizumab | 1-cmt SC mAb | CL=0.15 L/day, V1=3.8 L, t½≈22d |

### Disease Dynamics Compartments

| State Variable | Description | Calibration Source |
|---------------|-------------|-------------------|
| `Rfree` / `RC` | Free/bound IL-4Rα receptor (nM) | Kovalenko et al. 2016 (PK TMDD) |
| `pSTAT6` | Phosphorylated STAT6 (AU) | Hamilton et al. 2014 |
| `Th2` | Th2 cell index (skin AU) | Guttman-Yassky et al. 2019 |
| `ILC2` | ILC2 activation index | Kim et al. 2013 |
| `TARC` | TARC/CCL17 plasma (pg/mL) | Normal <450, AD 3000-8000 |
| `Eos_blood` | Blood eosinophil count (/µL) | Normal 100-500, AD 350+ |
| `IgE` | Total IgE (IU/mL) | AD baseline: 1000-5000 |
| `IL31` | IL-31 skin/plasma (pg/mL) | Baseline 45 pg/mL |
| `FLG` | Filaggrin expression (relative) | Normal=1.0, AD≈0.40 |
| `TEWL` | Trans-epidermal water loss | Normal <10, AD 20-40 g/m²/h |
| `SkinInfl` | Composite inflammation index | Maps to EASI |
| `EASI` | Eczema Area & Severity Index | Calibrated to SOLO-1/2 data |

### Treatment Scenarios (6 Scenarios)

| # | Scenario | Drug(s) | Key PD Mechanism |
|---|---------|---------|-----------------|
| 1 | No Treatment | — | Natural Th2 progression |
| 2 | TCS Only | Topical CS | STAT6 partial inhibition (+35%) |
| 3 | Dupilumab Q2W | 600mg load → 300mg Q2W | IL-4Rα blockade, RO >85% |
| 4 | Upadacitinib 30mg QD | 30mg daily | JAK1 inhibition >90%, broad cytokine block |
| 5 | Nemolizumab Q4W | 60mg Q4W SC | IL-31Ra blockade, rapid itch |
| 6 | Dupilumab + TCS | 300mg Q2W + TCS | Additive STAT6 inhibition |

### Simulated Response Endpoints at Week 16 (Model)

| Scenario | EASI | NRS | EASI-75 | IGA 0/1 | TARC Δ% |
|---------|------|-----|---------|---------|---------|
| No Treatment | ~28 | ~9.5 | 0% | 0% | ~0% |
| TCS Only | ~21 | ~7 | 5% | 10% | ~-20% |
| Dupilumab Q2W | ~10 | ~5.5 | 45% | 35% | ~-70% |
| Upadacitinib 30mg | ~8 | ~4.0 | 55% | 42% | ~-65% |
| Nemolizumab Q4W | ~22 | ~4.5 | 8% | 12% | ~-15% |
| Dupilumab + TCS | ~7 | ~5.0 | 60% | 50% | ~-75% |

*Model-predicted values; clinical trial values (SOLO-1/2, Rising Up) used for calibration.*

---

## Shiny App — 6 Interactive Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease severity sliders (EASI, IgE, Eos, TARC, BSA), AD subtype, comorbidities |
| 2. Drug PK Profiles | Dupilumab Cp, IL-4Rα RO%, Upadacitinib Cp, JAK1 inhibition % |
| 3. PD Biomarkers | TARC, eosinophil, IL-31, IgE, pSTAT6 time courses |
| 4. Clinical Endpoints | EASI, IGA, NRS itch, EASI-75/90 response bars, data table |
| 5. Scenario Comparison | Multi-arm overlay: EASI, NRS, TARC with summary table |
| 6. Barrier & Pruritus | FLG expression, TEWL, itch-inflammation feedback loop |

---

## Key Biological Insights

### Dupilumab Mechanism
- Blocks **shared IL-4Rα chain** → inhibits both IL-4 and IL-13 signaling simultaneously
- Downstream effects: pSTAT6↓ → TARC↓, FLG↑, IgE↓ (slow), eosinophil normalization
- Transient eosinophilia early (~4-8 weeks) before long-term reduction
- Ocular comorbidity (conjunctivitis ~15%) — mechanism under investigation (IL-13 suppression in conjunctiva)

### Upadacitinib (JAK1 Inhibitor)
- Broadest cytokine inhibition profile: TSLP, IL-4, IL-13, IL-31, IFN-γ all use JAK1
- Fastest itch improvement in head-to-head trials (Heads Up: upadacitinib > dupilumab at Wk 16)
- Concerns: herpes zoster reactivation, platelet changes, lipid effects (JAK2 off-target)

### Nemolizumab (Anti-IL-31Ra)
- **Itch-selective** biologic — does not significantly impact EASI at standard doses
- Rapid NRS improvement (within 2 weeks) without broad anti-inflammatory effect
- Best used in itch-predominant AD phenotype

### Epidermal Barrier Recovery
- FLG suppression by STAT6 creates a vicious cycle: barrier defect → more Th2 sensitization
- Dupilumab restores FLG expression via STAT6 inhibition within 4-8 weeks
- TEWL normalization is a delayed but objective marker of barrier recovery

---

## References Summary

| Category | Count |
|---------|-------|
| Epidemiology & burden | 4 |
| Epidermal barrier & FLG | 5 |
| Innate immunity & alarmins | 7 |
| IL-4/IL-13/STAT6 signaling | 4 |
| IL-31 & pruritus | 4 |
| Dupilumab clinical trials | 4 |
| Dupilumab PK/PD modeling | 2 |
| Tralokinumab | 2 |
| Lebrikizumab | 2 |
| Nemolizumab | 2 |
| JAK inhibitors | 4 |
| TARC & biomarkers | 3 |
| QSP & mathematical modeling | 4 |
| mrgsolve tools | 2 |
| Scoring systems | 3 |
| **Total** | **52** |

---

## Installation & Usage

```r
# Required R packages
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "purrr",
                   "shiny", "shinydashboard", "DT", "plotly"))

# Run mrgsolve simulation
source("ad_mrgsolve_model.R")

# Launch Shiny app
shiny::runApp("ad_shiny_app.R")
```

---

*Model developed as part of the QSP Disease Model Library. Parameters calibrated from published clinical trial data (SOLO-1/2, ECZTRA-1/2/3, Rising Up, JADE MONO-1/2).*
