# Cytokine Release Syndrome (CRS) — QSP Model

[![Model](https://img.shields.io/badge/QSP-CRS-red)](.) [![Framework](https://img.shields.io/badge/mrgsolve-ODE-blue)](.) [![Shiny](https://img.shields.io/badge/Shiny-Dashboard-orange)](.)

---

## Disease Overview

**Cytokine Release Syndrome (CRS)** is a life-threatening hyperinflammatory syndrome triggered by immune effector cell therapies, characterized by:

- **Massive T cell and macrophage activation** following CAR-T cell infusion or bispecific antibody administration
- **IL-6 cytokine storm** — the central driver of fever, hypotension, and hypoxia
- **Macrophage amplification loop** — IFN-γ and GM-CSF activate macrophages → IL-6 ↑↑ → further endothelial damage
- **Multi-organ involvement** — cardiovascular (Takotsubo cardiomyopathy), hepatic, renal, coagulopathy
- **ICANS** (Immune Effector Cell-Associated Neurotoxicity Syndrome) — BBB disruption, cerebral edema

CRS occurs in **50-90% of CD19 CAR-T patients** (grade ≥ 2 in 20-30%) and **virtually all bispecific antibody patients** (grade 1-2 most common). Mortality from severe CRS/ICANS has declined from >10% to <3% with modern management.

---

## Subtypes & Triggers

| Type | Typical Trigger | Onset | CRS Grade |
|------|----------------|-------|-----------|
| **CD19 CAR-T (CD28 costim.)** | ALL, DLBCL | Day 1-5 | Grade 2-4 |
| **CD19 CAR-T (4-1BB costim.)** | ALL, DLBCL | Day 2-7 | Grade 1-3 (milder) |
| **BCMA CAR-T** | Multiple Myeloma | Day 1-7 | Grade 1-3 |
| **Bispecific Ab (CD3×CD19)** | Blinatumomab | Hours | Grade 1-2 (cyclical) |
| **BCMA Bispecific Ab** | MM | Day 1-3 | Grade 1-3 |

---

## Mechanistic Map

[![CRS QSP Mechanistic Map](crs_qsp_model.png)](crs_qsp_model.svg)

> Click image to open high-resolution SVG. Contains **10 clusters, 100+ nodes** covering:
> immunotherapy trigger, T cell activation/expansion, T cell–derived cytokines,
> macrophage/myeloid activation, IL-6/JAK-STAT signaling network,
> endothelial/vascular response, ASTCT CRS grading, ICANS neurotoxicity,
> drug PK/PD (tocilizumab, siltuximab, DEX, ruxolitinib, anakinra), and clinical outcomes.

---

## Model Files

| File | Description |
|------|-------------|
| [`crs_qsp_model.dot`](crs_qsp_model.dot) | Graphviz source (10 clusters, 100+ nodes) |
| [`crs_qsp_model.svg`](crs_qsp_model.svg) | High-resolution vector image |
| [`crs_qsp_model.png`](crs_qsp_model.png) | 150 dpi raster image |
| [`crs_mrgsolve_model.R`](crs_mrgsolve_model.R) | mrgsolve ODE model + 5 treatment scenarios |
| [`crs_shiny_app.R`](crs_shiny_app.R) | 6-tab interactive Shiny dashboard |
| [`crs_references.md`](crs_references.md) | 55 PubMed-cited references |

---

## mrgsolve ODE Model

**22 ODE compartments** — CAR-T dynamics, cytokine cascade, macrophage activation, IL-6 signaling, endothelial/vascular response, CRS severity, organ damage + drug PK:

| Compartment Group | Variables |
|-------------------|-----------|
| CAR-T dynamics | CAR-T activated, exhausted, tumor burden |
| T cell cytokines | IFN-γ, IL-2, GM-CSF, TNF-α (T cell) |
| Macrophage cascade | MAC_ACT, IL-6, IL-1β, TNF-α (macrophage) |
| IL-6 downstream | pSTAT3, CRP, Ferritin |
| Endothelial/vascular | ENDO_ACT, VASC_PERM |
| CRS / organ damage | CRS severity (0-4), organ damage (0-1) |
| Drug PK | Tocilizumab (2-cmt), Siltuximab (1-cmt), Dexamethasone (2-cmt), Ruxolitinib (2-cmt), Anakinra (1-cmt) |

**5 Treatment Scenarios:**

| Scenario | Regimen | Indication |
|----------|---------|-----------|
| 1 | **Untreated** | Natural CRS progression (monitoring only) |
| 2 | **Tocilizumab** (8 mg/kg IV) | Grade 2 CRS (IL-6R blockade, 1-2 doses) |
| 3 | **Tocilizumab + Dexamethasone** | Grade 3 CRS (combined anti-cytokine) |
| 4 | **Ruxolitinib + Dexamethasone** | Refractory/Grade 4 CRS (JAK1/2 inhibition) |
| 5 | **Anakinra + Dexamethasone** | MAS-type/NLRP3-dominant / pediatric pattern |

---

## Shiny Dashboard (6 Tabs)

1. **Patient Profile** — CAR-T therapy type selector, tumor burden/dose sliders, treatment checkboxes, CRS grade trajectory, pathophysiology summary
2. **Drug PK** — Concentration-time profiles for all 5 drugs, PK parameter table, mechanism-of-action summary
3. **Cytokine Storm** — IFN-γ & IL-6 (key drivers), IL-1β & TNF-α, macrophage activation, CAR-T expansion, IL-6 downstream (STAT3, endothelial)
4. **Clinical Endpoints** — CRS severity (0-4), survival probability, organ damage, temperature/fibrinogen, ICANS score
5. **Scenario Comparison** — Side-by-side for all 5 treatment scenarios: IFN-γ, IL-6, CRS severity, survival, summary table
6. **Biomarkers & CRS Grading** — CRS biomarker index, ferritin (log scale), CRP, ICANS trajectory, ASTCT criteria tracker

**Run the Shiny app:**
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "tidyr", "plotly"))
shiny::runApp("crs_shiny_app.R")
```

---

## CRS Pathophysiology

```
CAR-T Infusion / Bispecific Ab Administration (Day 0)
         |
         v
Antigen Recognition → CAR-T Activation → IL-2 (autocrine), IFN-γ↑↑, GM-CSF↑
         |                                                     |
         v                                                     v
T Cell Expansion (peak Day 7-14)                    IFN-γ + GM-CSF →
(CD28: faster/severe)                               Macrophage M1 Activation
(4-1BB: slower/milder)                              (JAK1-STAT1 pathway)
                                                              |
                                                              v
                                                    IL-6↑↑↑ (main CRS driver)
                                                    IL-1β↑ (NLRP3 inflammasome)
                                                    TNF-α↑
                                                              |
                    ┌─────────────────────────────────────────┤
                    │                                         │
                    v                                         v
         IL-6 → Endothelial Activation              IL-1β + IL-6 → Fever
         → NO↑ → Vasodilation → HYPOTENSION         (hypothalamus PGE2)
         → VEGF/ANG2 → Vascular Leak → HYPOXIA
                    │
                    v
         BBB disruption → CNS macrophage activation
         → ICANS (encephalopathy, seizures)
                    │
                    v
         ASTCT CRS Grade 1→2→3→4
         (Fever → Hypotension/Hypoxia → Vasopressors → Vent/ECMO)

TREATMENT:
Tocilizumab → IL-6R blockade (gp130 + sIL-6R trans-signaling)  [Grade ≥2 CRS]
Siltuximab  → Direct IL-6 neutralization                        [Alternative]
DEX         → NF-κB inhibition → ↓IL-6, TNF, IFN-γ + BBB       [All ICANS grades]
Ruxolitinib → JAK1/2 inhibition → ↓pSTAT3, ↓pSTAT1             [Refractory]
Anakinra    → IL-1R blockade → ↓NLRP3 pathway                   [MAS-type/pediatric]
```

---

## Drug Target Summary

| Drug | Target | Key Effect in CRS |
|------|--------|------------------|
| **Tocilizumab** | IL-6R (gp130) | First-line CRS; blocks both cis- and trans-IL-6 signaling; FDA-approved |
| **Siltuximab** | Free IL-6 | Neutralizes IL-6 directly; does not reverse already-formed IL-6/sIL-6R complex |
| **Dexamethasone** | GR/NF-κB | Essential for ICANS (stabilizes BBB); suppresses multiple cytokines; used with TOCI |
| **Ruxolitinib** | JAK1/2 | Inhibits downstream IFN-γ/IL-6 signaling; breaks macrophage amplification feedback |
| **Anakinra** | IL-1R | Blocks NLRP3 inflammasome downstream; preferred in MAS overlap / pediatric CRS |

---

## ASTCT CRS Grading (2019 Consensus)

| Grade | Fever | Hypotension | Hypoxia | Treatment |
|-------|-------|-------------|---------|-----------|
| **1** | ≥38°C | None | None | Supportive |
| **2** | ≥38°C | IV fluids | Low-flow O₂ | **Tocilizumab** ± DEX |
| **3** | ≥38°C | 1 vasopressor | High-flow O₂ | **Toci + DEX** |
| **4** | ≥38°C | Multiple vasopressors | Mechanical ventilation / ECMO | **Toci + DEX + Ruxolitinib/Anakinra** |

Note: Any one of fever, hypotension, or hypoxia determines grade.

---

## Predictive Biomarkers for CRS Severity

| Biomarker | Timing | Significance |
|-----------|--------|-------------|
| **IFN-γ** ≥75 pg/mL Day 3 | Early | High specificity for severe CRS |
| **IL-6** >50 pg/mL | Day 2-5 | Best correlate; drives fever and hypotension |
| **CRP** >100 mg/L | Day 3-7 | Accessible clinical biomarker |
| **Ferritin** >2000 ng/mL | Day 5-10 | Macrophage activation; HLH overlap if >10,000 |
| **Fibrinogen** <2 g/L | Day 5+ | DIC risk; high CRS severity |
| **Tumor burden** (pre-infusion) | Pre-CAR-T | Higher → more antigen → more CRS |

---

## References

See [`crs_references.md`](crs_references.md) for 55 PubMed-cited references covering:
- Pivotal CAR-T trials (ZUMA-1, JULIET, TRANSCEND, KymRAH Phase 3)
- Bispecific antibody studies (blinatumomab, teclistamab, glofitamab)
- CRS grading systems (Lee 2014, ASTCT 2019)
- Macrophage biology (Giavridis 2018 Nat Med, Norelli 2018 Nat Med)
- Cytokine biomarkers (Teachey 2016 Cancer Discov)
- Tocilizumab PK/PD
- Ruxolitinib + anakinra for refractory CRS/HLH

---

*Model created: 2026-06-26 | Part of the QSP Disease Model Library*
