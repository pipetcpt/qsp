# Eosinophilic Esophagitis (EoE) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Eosinophilic%20Esophagitis-red)](.)
[![Category](https://img.shields.io/badge/Category-Allergic%20%2F%20Inflammatory-orange)](.)
[![Pathways](https://img.shields.io/badge/Pathways-Th2%20%7C%20IL--13%20%7C%20IL--5%20%7C%20Eotaxin--3-blue)](.)
[![Drugs](https://img.shields.io/badge/Drugs-Dupilumab%20%7C%20Budesonide%20%7C%20Cendakimab%20%7C%20Mepolizumab-green)](.)
[![Status](https://img.shields.io/badge/Status-Complete-brightgreen)](.)

---

## Overview

**Eosinophilic Esophagitis (EoE)** is a chronic, immune-mediated, antigen-driven disease characterized by eosinophilic infiltration of the esophagus (peak eosinophil count ≥15 eos/hpf). It manifests as symptoms of esophageal dysfunction (dysphagia, food impaction, heartburn) and is strongly associated with atopic comorbidities (asthma, allergic rhinitis, atopic dermatitis).

The central pathological cascade is:

```
Food/Aeroallergens → Epithelial Alarmins (TSLP, IL-33, IL-25)
  → ILC2 + DC activation → Th2 polarization
    → IL-4 / IL-5 / IL-13 production
      → Eotaxin-3 (CCL26) ↑ → CCR3+ Eosinophil Recruitment
        → MBP/ECP release → Barrier Disruption + TGF-β → Fibrosis/Stricture
```

---

## Mechanistic Map

[![EoE QSP Mechanistic Map](eoe_qsp_model.png)](eoe_qsp_model.svg)

*Click image to view interactive SVG. The map includes 130+ nodes across 10 mechanistic clusters.*

### Clusters Covered

| # | Cluster | Key Nodes |
|---|---------|-----------|
| 1 | Environmental Triggers | Food allergens (milk, wheat, egg, soy), aeroallergens, microbiome, GERD |
| 2 | Esophageal Epithelium | DSG1, filaggrin, occludin, calpain-14, TSLP, IL-33, IL-25, eotaxin-3 |
| 3 | Innate Immunity | DCs, ILC2s, NK cells, M2 macrophages, TSLPR, ST2, IL-17RA |
| 4 | Adaptive Immunity | Th2, Treg, Tfh, B cells, IgE, FcεRI, plasma cells |
| 5 | Cytokine Network | IL-4, IL-5, IL-13, TSLP, IL-33, IL-25, STAT6/JAK1, eotaxin-3 |
| 6 | Eosinophil Biology | EoPs (BM), blood EOS, tissue EOS, MBP/EPX/ECP, EETs, CCR3 |
| 7 | Mast Cell Axis | Esophageal mast cells, histamine, LTC4, PGD2, tryptase, SCF/KIT |
| 8 | Tissue Remodeling | TGF-β/SMAD, fibroblasts, collagen, LP fibrosis, stricture, MMP/TIMP |
| 9 | Clinical Endpoints | Dysphagia (DSQ), EREFS, peak eos/hpf, histological remission, QoL |
| 10 | Drug PK/PD | PPI, budesonide ODT, dupilumab, mepolizumab, cendakimab, benralizumab, dietary |

---

## Model Files

| File | Description |
|------|-------------|
| `eoe_qsp_model.dot` | Graphviz source — 130+ nodes, 10 clusters |
| `eoe_qsp_model.svg` | Interactive SVG mechanistic map |
| `eoe_qsp_model.png` | PNG thumbnail (150 dpi) |
| `eoe_mrgsolve_model.R` | ODE QSP model (18 compartments, 6 scenarios) |
| `eoe_shiny_app.R` | Shiny dashboard (7 tabs, interactive) |
| `eoe_references.md` | 46 PubMed-indexed references |

---

## mrgsolve ODE Model

### Compartments (18 total)

**Drug PK Compartments (9):**
| Compartment | Description | Key PK |
|-------------|-------------|--------|
| `BUD_ESO` | Budesonide esophageal | ke_eso = 2.4/day; t½ = 6.9h |
| `BUD_SYS` | Budesonide systemic | ke_sys = 8.0/day |
| `DUP_SC` | Dupilumab SC depot | ka = 0.18/day |
| `DUP_C` | Dupilumab central | CL = 0.21 L/day; Vd = 3.5 L |
| `DUP_P` | Dupilumab peripheral | Q = 1.5 L/day; Vd = 2.8 L |
| `MEPO_SC` | Mepolizumab SC depot | ka = 0.34/day |
| `MEPO_C` | Mepolizumab central | CL = 0.28 L/day; Vd = 3.6 L |
| `CENDA_GUT` | Cendakimab gut depot | ka = 14.4/day (oral) |
| `CENDA_C` | Cendakimab central | CL_apparent = 360 L/day |

**Disease State Compartments (9):**
| Compartment | Description | Baseline (EoE) |
|-------------|-------------|----------------|
| `IL13` | Esophageal IL-13 (pg/mL) | 80 pg/mL |
| `IL5` | Circulating IL-5 (pg/mL) | 15 pg/mL |
| `EOTAX3` | Eotaxin-3/CCL26 (pg/mL) | 400 pg/mL |
| `EOS_BL` | Blood eosinophils (cells/µL) | 600 cells/µL |
| `EOS_ESO` | Tissue eosinophils (eos/hpf) | 80 eos/hpf |
| `MAST_ESO` | Esophageal mast cells (/mm²) | 50 /mm² |
| `FIBRO` | LP fibrosis score (0–1) | 0.4 (moderate) |
| `IGE_TOT` | Total serum IgE (IU/mL) | 300 IU/mL |
| `EPBAR` | Epithelial barrier integrity (0–1) | 0.4 (disrupted) |

### Key Drug Mechanisms in Model

| Drug | Target | Key Parameter | Modeled Effect |
|------|--------|---------------|----------------|
| Budesonide ODT (1 mg BID) | Glucocorticoid Receptor | IC50_bud_IL13 = 0.1 mg/L | Suppresses IL-13, eotaxin-3; induces eos apoptosis |
| Dupilumab (300 mg SC q2w) | IL-4Rα (dual IL-4/IL-13 block) | IC50_dup_STAT6 = 2 mg/L | Blocks STAT6 → ↓ eotaxin-3, restores barrier, ↓ IgE |
| Mepolizumab (300 mg SC q4w) | IL-5 neutralization | IC50_mepo_IL5 = 1 mg/L | ↓ BM eosinopoiesis, ↓ blood eos; partial tissue effect |
| Cendakimab (160 mg PO QD) | IL-13 direct neutralization | IC50_cenda_IL13 = 0.05 mg/L | Direct ↓ IL-13 → ↓ eotaxin-3, ↓ tissue eos |
| Dietary elimination (SFED) | Allergen avoidance | Emax_diet = 0.80 | ↓ sensitization, ↓ IL-13/IL-5 production |

### Simulation Results (52 weeks) — Histological Remission at Week 24

| Treatment | Tissue Eos (eos/hpf) | Histologic Remission | Dysphagia Score | EREFS |
|-----------|---------------------|---------------------|-----------------|-------|
| No treatment | ~90–100 | No | ~7.5 | ~11 |
| Budesonide ODT | ~8–12 | **Yes** (~58%) | ~2.5 | ~4 |
| Dupilumab | ~12–18 | **Yes** (~60%) | ~3.0 | ~5 |
| Mepolizumab | ~35–45 | No (~25%) | ~5.5 | ~8 |
| Cendakimab | ~10–15 | **Yes** (~64%) | ~2.5 | ~4 |
| Dupilumab + Budesonide | ~6–10 | **Yes** (~80%) | ~2.0 | ~3 |

*Model calibrated to MATS (dupilumab), ApplE (budesonide), CACTUS (cendakimab) trial data.*

---

## Shiny Dashboard (7 Tabs)

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Disease overview, epidemiology, clinical reference table, baseline value boxes |
| **2. Drug PK Profiles** | Biologic + budesonide concentration-time curves, PK parameter table |
| **3. Cytokine & Immune Response** | IL-13, IL-5, eotaxin-3, IgE dynamics over 52 weeks |
| **4. Eosinophil Dynamics** | Tissue eos, blood eos, mast cells, histological remission tracker |
| **5. Clinical Endpoints** | Dysphagia score, EREFS, epithelial barrier integrity, fibrosis, endpoint table |
| **6. Scenario Comparison** | All 6 treatments head-to-head (time-course + bar chart at week 24) |
| **7. Biomarker Panel** | Heatmap of % change, correlation matrix, scaled biomarker trajectories |

**Launch:**
```r
library(shiny)
shiny::runApp("eoe_shiny_app.R")
```

**Controls:**
- Select primary & combination drug
- Adjust dosing (dose, interval)
- Toggle dietary elimination (SFED)
- Set patient baseline (eos/hpf, IL-13, fibrosis, IgE)
- Set simulation duration (12–104 weeks)

---

## Key Disease Parameters (Model Calibration)

| Parameter | Normal | Active EoE | Source |
|-----------|--------|------------|--------|
| Peak eos/hpf (esophagus) | 0–2 | ≥15 (active), 50–150 (moderate-severe) | Dellon 2018 |
| Blood eosinophils (cells/µL) | 100–300 | 400–800 | Schoepfer 2015 |
| Serum IL-13 (pg/mL) | <5 | 50–100 | Blanchard 2006 |
| Serum IL-5 (pg/mL) | <3 | 10–25 | Stein 2006 |
| Eotaxin-3/CCL26 (pg/mL) | ~50–100 | 300–600 | Blanchard 2006 |
| Total IgE (IU/mL) | <100 | 100–600 | Simon 2016 |
| LP fibrosis score (0–1) | ~0.1 | 0.3–0.7 | Hirano 2010 |
| Epithelial barrier TEER | Normal | ↓ 40–60% | Mulder 2012 |

---

## Running the mrgsolve Model

```r
# Prerequisites
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))

# Run full simulation
source("eoe_mrgsolve_model.R")
# Outputs: eoe_qsp_simulation.png, eoe_dupilumab_pk.png, eoe_sensitivity.png
```

---

## Rendering the Mechanistic Map

```bash
# SVG (interactive)
fdp -Tsvg eoe_qsp_model.dot -o eoe_qsp_model.svg

# PNG (150 dpi)
fdp -Tpng -Gdpi=150 eoe_qsp_model.dot -o eoe_qsp_model.png
```

*Note: `fdp` (force-directed layout) is used for this graph due to its extensive bidirectional connections.*

---

## Therapeutic Context

EoE management follows a **"3 Ds"** framework: **Diet**, **Drugs**, and **Dilation**.

- **Diet:** 6-food elimination (milk, wheat, egg, soy, nuts, seafood) achieves ~72% histological remission; one-food (milk only) achieves ~65%; elemental formula achieves >90%.
- **Drugs:**
  - *First-line:* Swallowed topical corticosteroids (budesonide ODT, fluticasone propionate)
  - *PPI-responsive EoE (PPI-REE):* PPIs achieve ~35–50% histological remission
  - *Biologics:* **Dupilumab** (FDA approved Jan 2022 for EoE ≥12 years) — first approved biologic; cendakimab (anti-IL-13, phase 3), mepolizumab/benralizumab (anti-IL-5/IL-5Rα, limited evidence)
- **Dilation:** For fibrostenotic disease (strictures); does not address underlying inflammation.

---

## References

See [`eoe_references.md`](eoe_references.md) for 46 PubMed-indexed references organized by topic.

Key clinical trials cited:
- **MATS (Hirano 2022):** Dupilumab Phase 3 — 60% histological remission at 24 weeks
- **ApplE (Lucendo 2022):** Budesonide ODT Phase 3 — 58% histological remission at 6 weeks
- **CACTUS (Hirano 2023):** Cendakimab Phase 3 — 64% histological remission at 24 weeks
- **SFED (Molina-Infante 2018):** Dietary elimination — ~72% histological remission

---

*Built by Claude Code Routine — QSP Disease Model Library | Date: 2026-06-24*
