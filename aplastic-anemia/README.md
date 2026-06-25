# Aplastic Anemia (AA) — QSP Model

> **재생불량성 빈혈** | Immune-mediated destruction of hematopoietic stem cells

[![DOT map](aa_qsp_model.png)](aa_qsp_model.svg)

---

## Overview

Aplastic Anemia (AA) is a life-threatening bone marrow failure syndrome caused by immune-mediated destruction of hematopoietic stem cells (HSCs). Autoreactive CD8+ cytotoxic T cells, driven by IFN-γ and TNF-α, eliminate HSCs via Fas/FasL and perforin/granzyme pathways, resulting in pancytopenia and hypocellular bone marrow. This QSP model integrates:

- T-cell immune dynamics (Teff/Treg balance, IFN-γ signaling)
- HSC compartment (destruction kinetics, EPO/TPO-mediated recovery)
- Multi-lineage hematopoiesis (erythroid, myeloid, megakaryocytic)
- Drug PK/PD: ATG (hATG/rATG), Cyclosporine (CsA), Eltrombopag (EPAG), Danazol
- Clinical endpoints: Hgb, ANC, PLT, BM cellularity
- PNH clone expansion and clonal evolution risk

---

## Disease Severity Classification (Camitta Criteria)

| Category | ANC           | Criteria (any 2/3)              |
|----------|---------------|---------------------------------|
| **VSAA** | < 0.2 ×10⁹/L | ANC, PLT, or ARC below threshold |
| **SAA**  | < 0.5 ×10⁹/L | PLT < 20 K/μL or Hgb < 8 g/dL  |
| **nSAA** | 0.5–1.5 ×10⁹/L | Pancytopenia without SAA criteria |

---

## Key Mechanistic Pathways

| Pathway                          | Key Components                              |
|----------------------------------|---------------------------------------------|
| Immune trigger & APC             | HLA-DR, pDC, MHC-I/II, cross-reactive Ag   |
| T cell activation & expansion    | CD4+ Th1, CD8+ CTL, TCR clonal expansion    |
| Cytokine storm                   | IFN-γ (key), TNF-α, IL-2, CXCL9/10         |
| HSC apoptosis pathways           | Fas/FasL, Perforin/Granzyme B, p53, ROS     |
| BM microenvironment              | MSC, Endothelial cells, SCF, CXCL12, Ang-1  |
| Regulatory immune balance        | Treg (FoxP3+) deficiency in AA              |
| PNH clone dynamics               | GPI-anchor deficiency → immune escape       |
| Erythroid recovery               | EPO feedback, CFU-E → Retic → RBC           |
| Megakaryocyte recovery           | TPO/c-Mpl, MK ploidy, platelet output       |

---

## Model Files

| File | Description |
|------|-------------|
| [`aa_qsp_model.dot`](aa_qsp_model.dot) | Graphviz mechanistic map (100+ nodes, 13 clusters) |
| [`aa_qsp_model.svg`](aa_qsp_model.svg) | SVG vector graphic |
| [`aa_qsp_model.png`](aa_qsp_model.png) | PNG thumbnail (150 dpi) |
| [`aa_mrgsolve_model.R`](aa_mrgsolve_model.R) | mrgsolve ODE model (20 compartments, 5 treatment scenarios) |
| [`aa_shiny_app.R`](aa_shiny_app.R) | Interactive Shiny dashboard (6 tabs) |
| [`aa_references.md`](aa_references.md) | 41 curated PubMed references |

---

## mrgsolve Model Specifications

### Compartments (20 ODEs)

| # | Compartment | Description |
|---|-------------|-------------|
| 1 | ATG_C | ATG central plasma (mg/L) |
| 2 | ATG_P | ATG peripheral tissue |
| 3 | CsA_C | Cyclosporine blood concentration (ng/mL) |
| 4 | EPAG_C | Eltrombopag plasma (μg/mL) |
| 5 | Danazol_C | Danazol plasma (ng/mL) |
| 6 | Teff | Autoreactive effector T cells (×10⁶/kg) |
| 7 | Treg | Regulatory T cells (×10⁶/kg) |
| 8 | HSC | Hematopoietic stem cell pool (% normal) |
| 9 | CFU_E | Erythroid progenitor pool |
| 10 | Retic | Reticulocyte pool |
| 11 | RBC | Circulating RBC → Hemoglobin |
| 12 | CFU_G | Granulocyte progenitor pool |
| 13 | ANC_pool | Circulating neutrophil pool → ANC |
| 14 | MK | Megakaryocyte pool |
| 15 | PLT_pool | Platelet pool |
| 16 | BM_score | BM cellularity score (0–1) |
| 17 | IFNg_c | IFN-γ concentration (pg/mL) |
| 18 | TNFa_c | TNF-α concentration (pg/mL) |
| 19 | IL2_c | IL-2 concentration (pg/mL) |
| 20 | PNH_clone | PNH clone fraction (0–1) |

### Treatment Scenarios

| Scenario | Treatment | Model Clinical Ref |
|----------|-----------|--------------------|
| 1 | No Treatment | Natural history — severe pancytopenia |
| 2 | hATG + CsA | Standard IST (Scheinberg 2011 NEJM) |
| 3 | hATG + CsA + EPAG | Triple IST (Townsley 2017 NEJM) |
| 4 | rATG + CsA + EPAG | NIH protocol (Peffault 2022 NEJM) |
| 5 | Allogeneic HSCT | MSD conditioning + engraftment |

---

## Shiny Dashboard Tabs

1. **Patient Profile** — Initial severity, value boxes (Hgb/ANC/PLT/BM), disease overview
2. **Drug PK** — ATG, CsA, EPAG concentration–time profiles; PK parameter table
3. **Hematopoiesis** — HSC pool, BM cellularity, erythroid/myeloid lineage, T-cell dynamics
4. **Clinical Endpoints** — Hgb, ANC, PLT, ARC trajectories; response classification
5. **Scenario Comparison** — Multi-scenario overlay plots; Day-180 CR/PR/NR rates
6. **Biomarkers & Clones** — IFN-γ/TNF-α/IL-2 dynamics; PNH clone expansion; MDS risk

---

## Key Parameter Calibration Notes

| Parameter | Value | Source |
|-----------|-------|--------|
| hATG t½ | ~7 h | Zingman 1990; Fetterly 2019 |
| rATG CL | ~50% of hATG CL | Scheinberg 2011 PK sub-study |
| CsA target trough | 150–250 ng/mL | Killick 2016 BSH Guidelines |
| EPAG EC50 (c-Mpl) | ~60 μg/mL | van der Straaten 2021 PopPK |
| EPAG t½ | ~21 h | FDA label; Olnes 2012 |
| IFN-γ (active AA) | 30–200 pg/mL | Young 2018; rodgers 2020 |
| BM cellularity (SAA) | <25% | Camitta criteria |
| CR rate: hATG+CsA | ~50% | Scheinberg 2011 |
| CR rate: +EPAG | ~68–74% | Townsley 2017; Peffault 2022 |
| PNH clone prevalence | 50–60% at Dx | Dezern 2014 |

---

## Running the Model

```r
# 1. Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "purrr"))

# 2. Run mrgsolve scenarios
source("aa_mrgsolve_model.R")
results <- run_all_scenarios()
plots   <- plot_results(results)
print(plots$Hgb)

# 3. Launch Shiny dashboard
shiny::runApp("aa_shiny_app.R")
```

```bash
# Graphviz rendering
dot -Tsvg aa_qsp_model.dot -o aa_qsp_model.svg
dot -Tpng -Gdpi=150 aa_qsp_model.dot -o aa_qsp_model.png
```

---

## References (Key)

- Young NS. *N Engl J Med.* 2018;379:1643–1656. (Pathophysiology)
- Scheinberg P et al. *N Engl J Med.* 2011;365:430–438. (hATG vs rATG)
- Townsley DM et al. *N Engl J Med.* 2017;376:1540–1550. (EPAG + IST)
- Peffault de Latour R et al. *N Engl J Med.* 2022;386:11–23. (rATG+CsA+EPAG)
- Olnes MJ et al. *N Engl J Med.* 2012;367:11–19. (EPAG monotherapy)

Full list: [`aa_references.md`](aa_references.md) (41 entries)

---

*Model created 2026-06-25 by Claude Code Routine (CCR)*
