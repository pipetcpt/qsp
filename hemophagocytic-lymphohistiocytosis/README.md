# Hemophagocytic Lymphohistiocytosis (HLH) — QSP Model

[![Model](https://img.shields.io/badge/QSP-HLH-red)](.) [![Framework](https://img.shields.io/badge/mrgsolve-ODE-blue)](.) [![Shiny](https://img.shields.io/badge/Shiny-Dashboard-orange)](.)

---

## Disease Overview

**Hemophagocytic Lymphohistiocytosis (HLH)** is a life-threatening hyperinflammatory syndrome driven by uncontrolled immune activation, characterized by:

- **Cytokine storm** led by IFN-γ, IL-18, TNF-α, and IL-6
- **Macrophage activation** and pathological hemophagocytosis (engulfment of blood cells by macrophages in bone marrow, liver, spleen)
- **Multi-organ damage** — liver failure, coagulopathy/DIC, cytopenias, CNS involvement
- **NK cell / CTL cytotoxicity defect** (primary HLH) or trigger-induced immune dysregulation (secondary HLH/MAS)

Without treatment, **mortality exceeds 50% within 5 weeks**.

---

## Subtypes

| Type | Basis | Trigger | Typical Treatment |
|------|-------|---------|-------------------|
| **Primary / Familial HLH (FHL)** | Genetic (PRF1, UNC13D, STX11, STXBP2, RAB27A…) | Often absent or minor | HLH-2004 → HSCT (curative) |
| **Secondary HLH** | Acquired; intact genetics | Infection (EBV, CMV), malignancy | Treat trigger + DEX ± ETOP ± CsA |
| **MAS (Macrophage Activation Syndrome)** | sJIA/AOSD, SLE | Autoimmune disease flare | High-dose steroids ± anakinra/tocilizumab |

---

## Mechanistic Map

[![HLH QSP Mechanistic Map](hlh_qsp_model.png)](hlh_qsp_model.svg)

> Click image to open high-resolution SVG. Contains **10 clusters, 100+ nodes** covering:
> genetic background, triggering events, NK/CTL dysfunction, APC–T cell crosstalk,
> cytokine storm network (IFN-γ, IL-6, TNF-α, IL-18, IL-10, IL-12),
> macrophage activation/hemophagocytosis, multi-organ damage,
> clinical biomarkers/HScore, drug PK/PD, and clinical outcomes.

---

## Model Files

| File | Description |
|------|-------------|
| [`hlh_qsp_model.dot`](hlh_qsp_model.dot) | Graphviz source (10 clusters, 100+ nodes) |
| [`hlh_qsp_model.svg`](hlh_qsp_model.svg) | High-resolution vector image |
| [`hlh_qsp_model.png`](hlh_qsp_model.png) | 150 dpi raster image |
| [`hlh_mrgsolve_model.R`](hlh_mrgsolve_model.R) | mrgsolve ODE model + 5 treatment scenarios |
| [`hlh_shiny_app.R`](hlh_shiny_app.R) | 6-tab interactive Shiny dashboard |
| [`hlh_references.md`](hlh_references.md) | 52 PubMed-cited references |

---

## mrgsolve ODE Model

**20 ODE compartments** — disease state, cytokines, macrophage cascade, organ damage + drug PK:

| Compartment Group | Variables |
|-------------------|-----------|
| Immune cells | NK cells, CTL, activated T cells, APCs |
| Cytokines (pg/mL) | IFN-γ, IL-6, TNF-α, IL-18, IL-10, IL-12 |
| Macrophage cascade | Macrophage activation, hemophagocytosis index, ferritin |
| Organ damage | Bone marrow suppression, liver damage, coagulopathy |
| Drug PK | DEX (2-cmt oral), etoposide (2-cmt IV), CsA (2-cmt oral), emapalumab (mAb TMDD), anakinra (SC 1-cmt), ruxolitinib (2-cmt oral) |

**5 Treatment Scenarios:**

| Scenario | Regimen | Indication |
|----------|---------|-----------|
| 1 | **Untreated** | Natural disease progression control |
| 2 | **HLH-2004** (DEX + Etoposide + CsA) | Standard first-line for primary/secondary HLH |
| 3 | **Emapalumab + DEX** | Primary/refractory HLH (FDA-approved 2018) |
| 4 | **Anakinra + DEX** | MAS / sJIA-HLH (IL-1 pathway dominant) |
| 5 | **Ruxolitinib + DEX** | Refractory/relapsed HLH (salvage) |

---

## Shiny Dashboard (6 Tabs)

1. **Patient Profile** — HLH subtype selector, value boxes (IFN-γ, ferritin, survival %), disease overview plot
2. **Drug PK** — Concentration-time profiles for all 6 drugs, PK parameter table, mechanism-of-action table
3. **Cytokine Storm** — Pro-inflammatory cytokines (IFN-γ, TNF-α, IL-18, IL-12), macrophage activation, hemophagocytosis index
4. **Clinical Endpoints** — Survival probability, organ damage (liver, coagulation), bone marrow suppression
5. **Scenario Comparison** — Side-by-side comparison of all 5 treatment scenarios for IFN-γ, ferritin, survival, HScore
6. **Biomarkers & HScore** — HScore trajectory, ferritin (log scale), sCD25/NK activity, fibrinogen/triglycerides, HLH-2004 criteria tracker

**Run the Shiny app:**
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "tidyr", "plotly"))
shiny::runApp("hlh_shiny_app.R")
```

---

## Key Pathophysiology Captured

```
Genetic defect (PRF1/UNC13D...)    External trigger (EBV/CMV/sJIA)
          ↓                                    ↓
    NK/CTL dysfunction ──────────────→ APCs not cleared
          ↓                                    ↓
    Persistent APC-T cell crosstalk → T cell hyperactivation
          ↓                                    ↓
    IFN-γ ↑↑↑ ──────────────────→ Macrophage activation
    IL-18 ↑  (synergy)                         ↓
    IL-12 ↑  (amplification)        Hemophagocytosis
          ↓                         Ferritin ↑↑↑ (>10,000 ng/mL)
    Cytokine storm                  sCD25 ↑↑↑ (>2,400 U/mL)
    TNF-α, IL-6, IL-1β ↑                       ↓
          ↓                         BM suppression → cytopenias
    Multi-organ damage              Coagulopathy (DIC)
    (Liver, Kidney, CNS, Lung)                  ↓
          ↓                              ↓ DEATH (untreated)
    HScore ≥169 → Diagnosis
```

---

## Drug Target Summary

| Drug | Target | Key Effect in HLH |
|------|--------|------------------|
| Dexamethasone | GR/NF-κB | ↓ IFN-γ, TNF-α, macrophage activation |
| Etoposide | Topoisomerase II | Depletes activated T cells (reduces cytokine source) |
| Cyclosporine A | Calcineurin/NFAT | ↓ IL-2, IFN-γ; blocks T cell proliferation |
| **Emapalumab** | IFN-γ (direct neutralization) | First targeted therapy; FDA-approved for primary HLH |
| Anakinra | IL-1 receptor | Blocks IL-1β signal; most effective in MAS/NLRP3-driven |
| Ruxolitinib | JAK1/2 | Inhibits downstream signaling of IFN-γ, IL-6 |

---

## Diagnostic Criteria

**HLH-2004 Criteria** (≥5 of 8 = diagnosis):
1. Fever (>38.5°C)
2. Splenomegaly
3. Cytopenias (≥2 cell lines: ANC <1×10⁹/L, Hgb <90 g/L, Plt <100×10⁹/L)
4. Hypertriglyceridemia ≥3 mmol/L or hypofibrinogenemia <1.5 g/L
5. Hemophagocytosis in BM/spleen/LN
6. Low/absent NK cell activity
7. Ferritin ≥500 ng/mL (>10,000 has ~93% specificity)
8. sCD25 ≥2,400 U/mL

**HScore ≥169** → 93% sensitivity, 86% specificity for HLH.

---

## References

See [`hlh_references.md`](hlh_references.md) for 52 PubMed-cited references covering:
- Primary HLH genetics (perforin, MUNC13-4, Rab27a, etc.)
- IFN-γ biology and cytokine storm mechanisms
- MAS diagnostic criteria (Ravelli 2016)
- HLH-2004 and HLH-94 treatment outcomes
- Emapalumab pivotal trial (Locatelli 2020, NEJM)
- Novel therapies (anakinra, ruxolitinib)
- HScore validation (Fardet 2014)
- HSCT outcomes
- Pharmacokinetics/pharmacometrics

---

*Model created: 2026-06-26 | Part of the QSP Disease Model Library*
