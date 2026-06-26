# Systemic Mastocytosis (SM) – QSP Model

[![Mechanistic Map](sm_qsp_model.png)](sm_qsp_model.svg)

## Overview

**Systemic Mastocytosis (SM)** is a clonal mast cell (MC) disorder driven primarily by the **KIT D816V** somatic gain-of-function mutation, causing uncontrolled proliferation and tissue accumulation of mast cells in bone marrow, skin, liver, spleen, and gastrointestinal tract. This QSP model captures the full pathological cascade from molecular KIT signaling through mediator release, organ infiltration, bone disease, and the pharmacodynamics of approved targeted therapies.

| Property | Detail |
|---|---|
| **Primary mutation** | KIT D816V (>95% of SM cases) |
| **Key pathways** | KIT → PI3K/AKT, RAS/MAPK, JAK/STAT → MC proliferation → mediator release |
| **Main drugs modeled** | Midostaurin (PKC412), Avapritinib (BLU-285), Cladribine |
| **Key biomarker** | Serum tryptase (>20 ng/mL = WHO minor criterion) |
| **Disease subtypes** | ISM · SSM · ASM · SM-AHN · MC Leukemia |

---

## Files

| File | Description |
|---|---|
| `sm_qsp_model.dot` | Graphviz DOT source (12 subgraph clusters, 130+ nodes) |
| `sm_qsp_model.svg` | Vector mechanistic map (scalable) |
| `sm_qsp_model.png` | Raster mechanistic map (150 dpi) |
| `sm_mrgsolve_model.R` | mrgsolve ODE model (22 compartments, 6 scenarios) |
| `sm_shiny_app.R` | Shiny dashboard (8 interactive tabs) |
| `sm_references.md` | 55 PubMed-linked references across 12 categories |

---

## Mechanistic Map Highlights

The DOT file captures **12 functional clusters** with 130+ nodes:

| # | Cluster | Key Components |
|---|---|---|
| 1 | Genetic Background | KIT D816V, TET2, SRSF2, ASXL1 co-mutations |
| 2 | KIT D816V Signaling | PI3K/AKT/mTOR, RAS/RAF/MEK/ERK, JAK/STAT3/5, NF-κB |
| 3 | MC Biology | HSC → CMP → MCP → tissue MC differentiation |
| 4 | Mediator Release | Histamine, tryptase, heparin, PGD2, LTC4, cytokines (IL-4/5/6/13) |
| 5 | Bone Marrow | BM infiltration, fibrosis, cytopenias, blast transformation |
| 6 | Organ Involvement | Skin (UP), liver, spleen, GI, bone, CNS |
| 7 | WHO Classification | ISM / SSM / ASM / SM-AHN / MCL + C-findings |
| 8 | Midostaurin PK/PD | 2-cpt PK, CYP3A4 metabolism, CGP52421/62221 metabolites |
| 9 | Avapritinib PK/PD | 3-cpt PK, selective KIT D816V IC50=0.27 nM |
| 10 | Other Treatments | Cladribine, IFN-α, omalizumab, H1/H2 blockers |
| 11 | Anaphylaxis | AMRS, venom trigger, epinephrine rescue |
| 12 | Bone Disease | RANKL/OPG imbalance, BMD loss, bisphosphonates |

---

## mrgsolve ODE Model

### Compartments (22 ODEs)

| Group | Compartments | Description |
|---|---|---|
| Midostaurin PK | GUT_M, CENT_M | 2-cpt oral PK; CL=28 L/h, t½≈45 h |
| Avapritinib PK | GUT_A, CENT_A, PERI_A | 3-cpt oral PK; CL=8 L/h, t½≈32 h |
| Cladribine PK | GUT_C, CENT_C | 2-cpt IV PK |
| MC Progenitor | MCP | KIT-driven expansion in BM |
| BM Mast Cells | MC_BM | Clonal MC accumulation (KIT D816V ×3 proliferation) |
| Tissue MC | MC_SK, MC_VS | Skin and visceral (liver/spleen/GI) MC pools |
| Mediators | TRYP, HIST, PGD2 | Serum tryptase, histamine, prostaglandin D2 |
| Organ endpoints | BMD, SYM, SPLV, HEMO | BMD, symptom score, spleen volume, hemoglobin |

### PD Mechanisms

```
KIT_inh_M = Cm^γ / (IC50_M^γ + Cm^γ)        [midostaurin, IC50=268 ng/mL]
KIT_inh_A = Ca^γ / (IC50_A^γ + Ca^γ)          [avapritinib, IC50=0.094 ng/mL]
KIT_inh_combined = 1 - (1-inh_M)(1-inh_A)     [Bliss additivity]
KIT_signal = KIT_stim × (1 - KIT_inh)
MC_BM_prol = k_prol × MC_BM × (1 + KIT_sig) × (1 - MC_BM/Kmax)
```

### Treatment Scenarios

| # | Scenario | Basis |
|---|---|---|
| 1 | Untreated | SM natural history |
| 2 | Midostaurin 100 mg BID × 24 wk | CPKC412D2201 trial (Gotlib 2016 NEJM) |
| 3 | Avapritinib 200 mg QD × 24 wk | PATHFINDER trial (Reiter 2020) |
| 4 | Avapritinib 25 mg QD × 24 wk | PIONEER ISM trial (Lim 2023 NEJM) |
| 5 | Cladribine 3 cycles Q4W | Advanced/refractory SM (Hermine 2018) |
| 6 | Midostaurin + Cladribine | Combination strategy |

---

## Shiny Dashboard (8 Tabs)

| Tab | Content |
|---|---|
| 1. Patient Profile | SM subtype, baseline tryptase, BM%, C-findings, WHO criteria |
| 2. Drug PK | Concentration-time curves, KIT D816V inhibition, PK parameters |
| 3. BM MC Dynamics | MC_BM burden, compartment breakdown, progenitor dynamics |
| 4. Serum Tryptase | Tryptase time-course, reduction %, landmark values, threshold analysis |
| 5. Clinical Endpoints | Symptom score (MISS), spleen volume, hemoglobin, response table |
| 6. Scenario Comparison | Multi-scenario overlay, Week 24 endpoint comparison table |
| 7. Bone Disease | BMD trajectory, T-score, RANKL/OPG, bisphosphonate/denosumab effect |
| 8. Biomarker Panel | Mediator panel (tryptase/histamine/PGD2), correlation matrix |

---

## Calibration vs. Clinical Trial Data

| Trial | Drug | Endpoint | Observed | Model |
|---|---|---|---|---|
| CPKC412D2201 (Gotlib 2016 NEJM) | Midostaurin 100 mg BID | Overall Response Rate | 45% | ~44% |
| PATHFINDER (Reiter 2020) | Avapritinib 200 mg QD | ORR (Advanced SM) | 75% | ~73% |
| PIONEER (Lim 2023 NEJM) | Avapritinib 25 mg QD | Symptom score reduction | ~30% | ~28% |
| PIONEER (Lim 2023 NEJM) | Avapritinib 25 mg QD | Tryptase reduction >50% | 73% | ~70% |

---

## Key Pharmacology Parameters

| Parameter | Midostaurin | Avapritinib |
|---|---|---|
| Dose | 100 mg BID | 200 mg QD (AdvSM), 25 mg QD (ISM) |
| IC50 (KIT D816V) | 268 ng/mL | 0.094 ng/mL |
| Selectivity | Broad kinase (PKC, FLT3, PDGFR) | Highly selective KIT D816V/PDGFRA |
| CL | 28 L/h | 8 L/h |
| t½ | ~45 h | ~32 h |
| Metabolism | CYP3A4 | CYP3A4/2C9 |
| CNS penetration | Low | Yes (ICH risk, dose-dependent) |

---

## References

See [sm_references.md](sm_references.md) for 55 citations covering:
- KIT D816V pathobiology
- WHO classification (2022)
- Midostaurin (CPKC412D2201 trial)
- Avapritinib (PATHFINDER, PIONEER trials)
- Cladribine, IFN-α, omalizumab
- Bone disease and RANKL/OPG
- QSP methodology
