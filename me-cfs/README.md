# ME/CFS QSP Model

**Myalgic Encephalomyelitis / Chronic Fatigue Syndrome**

[![Mechanistic Map](mecfs_qsp_model.png)](mecfs_qsp_model.svg)

---

## Disease Overview

Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) is a complex, debilitating multi-system disease affecting an estimated 17–24 million people in the USA and 0.3–2.8% of the global population. It is characterized by:

- **Post-Exertional Malaise (PEM)**: The pathological hallmark — worsening of symptoms 24-72h after physical or cognitive exertion
- **Profound fatigue**: Not relieved by rest, unrelated to activity level
- **Cognitive impairment** ("brain fog"): Memory deficits, processing speed reduction
- **Autonomic dysfunction**: POTS, orthostatic intolerance, reduced HRV
- **Unrefreshing sleep**: Altered REM/Stage 3 architecture
- **Pain**: Central sensitization, myalgia, headache

The disease frequently follows viral infections (EBV, HHV-6, SARS-CoV-2, enteroviruses), establishing it as a post-infectious syndrome with persistent immune activation and metabolic dysregulation.

---

## Pathophysiology

### Core Mechanistic Pillars

| System | Mechanism | Key Finding |
|--------|-----------|-------------|
| Immune | NK cell dysfunction | 30-50% reduced cytotoxicity |
| Immune | CD8+ T cell exhaustion | PD-1/LAG-3/TIM-3 upregulation |
| Immune | Autoantibodies | β2-Adrenergic receptor + M1-Muscarinic autoAb in ~60% |
| Mitochondrial | PDH inhibition | PDK1↑ (IFN-γ driven) → Warburg-like shift |
| Mitochondrial | ATP deficit | Complex I/IV activity ↓, ETC impairment |
| Neurological | Neuroinflammation | TSPO PET (¹¹C-PK11195) signal ↑ in brain |
| Neurological | Kynurenine pathway | IDO1↑ → Quinolinic acid ↑ → Hippocampal damage |
| Autonomic | POTS | β2AR autoAb → dysautonomia |
| Autonomic | Hypocortisolism | HPA axis suppression, GR resistance |
| MCAS | Mast cell activation | Tryptase ↑, histamine-mediated neuroinflammation |

### Critical Mechanistic Cascade: IFN-γ → PDK1 → PDH → ATP Deficit

```
Viral persistence / Immune activation
    ↓
IFN-γ ↑ (Th1/NK)
    ↓
PDK1 ↑ (Pyruvate Dehydrogenase Kinase 1)
    ↓
PDH ↓ (Pyruvate Dehydrogenase — phosphorylated/inhibited)
    ↓
Pyruvate → Lactate (instead of Acetyl-CoA → TCA cycle)
    ↓
↓ NADH → ↓ Complex I → ↓ ATP synthesis
    ↓
ATP Deficit → Fatigue + PEM (anaerobic threshold ↓)
```

This "Warburg-like" metabolic shift in immune/muscle cells was confirmed by:
- Fluge et al. (2017) JCI Insight: metabolic profiling
- Naviaux et al. (2016) PNAS: CDe metabolomics
- Tomas et al. (2017) PNAS: metabolite profiling

---

## QSP Model Files

| File | Description |
|------|-------------|
| [`mecfs_qsp_model.dot`](mecfs_qsp_model.dot) | Graphviz mechanistic map (220+ nodes, 12 subgraphs) |
| [`mecfs_qsp_model.svg`](mecfs_qsp_model.svg) | High-resolution vector image |
| [`mecfs_qsp_model.png`](mecfs_qsp_model.png) | Rasterized image (150 dpi) |
| [`mecfs_mrgsolve_model.R`](mecfs_mrgsolve_model.R) | mrgsolve ODE model (25 compartments, 6 scenarios) |
| [`mecfs_shiny_app.R`](mecfs_shiny_app.R) | Interactive Shiny dashboard (8 tabs) |
| [`mecfs_references.md`](mecfs_references.md) | 62 PubMed-cited references |

---

## Mechanistic Map Structure (12 Subgraphs)

1. **Triggering Events & Viral Persistence** — EBV, HHV-6, SARS-CoV-2, enteroviruses, latent reservoir
2. **Innate Immune Activation** — TLR3/4/7/9, cGAS-STING, IRF3, NLRP3 inflammasome, NK cells, macrophages
3. **Adaptive Immune Dysregulation** — T cell exhaustion, B cell dysfunction, autoantibodies (β2AR, M1R)
4. **Cytokine Network** — IL-6, TNF-α, IFN-γ, TGF-β, CXCL10, cytokine dysregulation
5. **Mast Cell Activation (MCAS)** — MC degranulation, tryptase, histamine, leukotrienes, prostaglandins
6. **HPA Axis & Neuroendocrine** — CRH, ACTH, hypocortisolism, GR resistance, circadian disruption
7. **Autonomic Nervous System** — POTS, orthostatic intolerance, HRV reduction, RAAS, small fiber neuropathy
8. **Mitochondrial Dysfunction** — PDH/PDK1, Complex I-IV, ROS, ATP deficit, NAD+ pool, PGC-1α
9. **Neuroinflammation & CNS** — Microglial activation, TSPO, BBB disruption, central sensitization, kynurenine pathway
10. **Post-Exertional Malaise (PEM)** — Anaerobic threshold ↓, Day-2 CPET, immune flare, metabolic collapse
11. **Clinical Endpoints & Biomarkers** — FSS, SF-36, CPET VO2max, NK cytotoxicity, morning cortisol, HRV RMSSD
12. **Drug PK/PD** — LDN, Pyridostigmine, Rituximab, BC007, NAD+ precursors, Rintatolimod, Aripiprazole, IVIG

---

## mrgsolve ODE Model

### State Variables (25 Compartments)

| Category | Variables |
|----------|-----------|
| Immune | V (viral), IFN (Type I), NK (cells), Tex (exhaustion), AutoAb |
| Cytokines | IL6, TNFa, IFNg, NLRP3state |
| MCAS | MC_act, Histamine |
| HPA Axis | CRH, Cortisol |
| Autonomic | NE_plasma, HRV_index |
| Mitochondrial | PDH_act, ATP_state, ROS_state |
| CNS | Neuro_inf, Cog_func |
| PEM/Fatigue | PEM_sens, Fatigue |
| Drug PK | LDN_cp, Pyr_cp, Rit_cp, NADpool |

### Treatment Scenarios

| # | Scenario | Drug(s) | Mechanism |
|---|----------|---------|-----------|
| 1 | No Treatment | — | Natural history progression |
| 2 | LDN Monotherapy | LDN 4.5 mg/day | TLR4 inhibition → ↓Neuroinflammation |
| 3 | Pyridostigmine | 30 mg TID | AChE inhibition → ANS restoration → POTS↓ |
| 4 | Rituximab | 1000 mg IV ×2 | B cell depletion → AutoAb reduction |
| 5 | NAD+ Precursors | NMN/NR 500 mg/day | NAD+ restoration → Sirtuin → PGC-1α → Mito↑ |
| 6 | Combination | LDN + NAD+ + Pyridostigmine | Multi-target synergy |

---

## Shiny App Features (8 Tabs)

1. **Patient Profile** — Disease overview, radar chart, symptom burden timeline, pathophysiology table
2. **PK Profiles** — LDN, Pyridostigmine, Rituximab, NAD+ pool kinetics
3. **Immune & Cytokines** — IL-6, TNF-α, IFN-γ, NK cells, T cell exhaustion, autoantibodies, MCAS
4. **Mitochondria & Energy** — PDH activity, ATP state, ROS, metabolic cascade visualization
5. **Neuroinflammation/CNS** — NI index, cognitive function, HPA axis, ANS, kynurenine pathway
6. **Clinical Endpoints** — FSS fatigue score, PEM sensitivity, HRV, cortisol, outcome table at 6 months
7. **Scenario Comparison** — Multi-treatment comparison with fatigue, ATP, NI, PEM outcomes
8. **Biomarker Dashboard** — Bar chart vs. healthy reference, drug target table, response heatmap

---

## Key Pharmacological Targets

| Drug | Target | Evidence | Status |
|------|--------|----------|--------|
| **LDN** (1.5-4.5 mg/d) | TLR4 / Microglial | Multiple open-label RCTs | Phase II |
| **Pyridostigmine** (30-60 mg TID) | AChE / Vagal tone | POTS/ME clinical trials | Phase III |
| **Rituximab** (1g IV ×2) | CD20+ B cells / AutoAb | Phase III RCT (negative) | Completed |
| **BC007** | β2AR autoantibodies | Open-label pilot | Phase II |
| **NMN/NR** (250-500 mg/d) | NAD+ / Sirtuins | Phase I/II | Active |
| **Rintatolimod** (Ampligen) | TLR3 / Antiviral | FDA fast-track | Phase III |
| **Aripiprazole** (0.5-2 mg/d) | D2 partial agonist | Open-label | Phase II |
| **Cromolyn** | Mast cell stabilizer | Off-label | Clinical use |
| **IVIG** | Fc receptors / Immune modulation | Several RCTs | Phase II/III |

---

## References

See [`mecfs_references.md`](mecfs_references.md) for 62 PubMed-cited references organized by:
- Foundational & Epidemiology (6 refs)
- Viral Triggers & Pathogen Persistence (5 refs)
- Immune Dysregulation (6 refs)
- Autoantibodies & Autonomic Dysfunction (5 refs)
- Mitochondrial Dysfunction & Energy Metabolism (7 refs)
- Neuroinflammation & CNS (6 refs)
- HPA Axis & Neuroendocrine (3 refs)
- Post-Exertional Malaise (4 refs)
- MCAS & Mast Cells (3 refs)
- Pharmacological Treatments (9 refs)
- QSP & Systems Biology (7 refs)

---

*Date: 2026-06-28 | QSP Library — ME/CFS Model*
