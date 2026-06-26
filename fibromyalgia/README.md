# Fibromyalgia (FM) — Quantitative Systems Pharmacology Model

[![Mechanistic Map](fm_qsp_model.png)](fm_qsp_model.svg)

## Overview

Fibromyalgia (FM) is a chronic widespread pain disorder affecting 2–5% of the population, predominantly women. Its hallmark is **central sensitization** — amplified pain processing in the central nervous system — combined with deficient descending pain inhibition, neuroendocrine dysregulation, non-restorative sleep, and autonomic imbalance. This QSP model integrates all of these mechanisms with the pharmacokinetics and pharmacodynamics of four approved/commonly used pharmacotherapies.

---

## Disease Pathophysiology Modeled

| Subsystem | Key Components |
|-----------|----------------|
| **Peripheral Sensitization** | NGF↑, PGE2, TRPV1, NaV1.7/1.8, DRG nociceptors (Aδ/C-fiber) |
| **Spinal Dorsal Horn** | WDR neurons, NMDA receptors, Substance P, wind-up, spinal LTP |
| **Central Sensitization** | Long-term potentiation (LTP), BDNF/TrkB, microglia, IL-1β, KCC2 disinhibition |
| **Descending Modulation** | PAG-RVM axis, noradrenergic (LC), serotonergic (raphe), DPMS, DNIC |
| **Supraspinal Brain** | ACC, insula, PFC, amygdala, default mode network hyperactivity |
| **HPA Axis** | CRH-ACTH-cortisol loop, negative feedback, blunted diurnal cortisol, GH/IGF-1 deficiency |
| **Autonomic Nervous System** | Sympathetic hyperactivity, reduced HRV, orthostatic dysregulation |
| **Sleep** | Adenosine pressure, SWS depth, alpha-delta intrusion, non-restorative sleep |
| **Neuroinflammation** | Microglia activation, IL-1β, TNF-α, NLRP3 inflammasome, mast cells |
| **Drug PK/PD** | Duloxetine, Pregabalin, Milnacipran, Amitriptyline |

---

## Model Files

| File | Description |
|------|-------------|
| `fm_qsp_model.dot` | Graphviz source — 100+ nodes, 10 subgraph clusters |
| `fm_qsp_model.svg` | Scalable mechanistic map (vector) |
| `fm_qsp_model.png` | Mechanistic map (150 dpi raster) |
| `fm_mrgsolve_model.R` | mrgsolve ODE model (30 compartments, 6 clinical scenarios) |
| `fm_shiny_app.R` | Shiny interactive dashboard (6 tabs) |
| `fm_references.md` | 60 PubMed-linked references |

---

## ODE Model Structure (mrgsolve, 30 compartments)

```
Drug PK (10 CMT)
  DUL_gut → DUL_cent ⇌ DUL_peri  [2-cpt, CL/F 54 L/h, Vd 1640 L]
  PRE_gut → PRE_cent               [1-cpt, CLr 6.8 L/h, renal]
  MIL_gut → MIL_cent               [1-cpt, CL 50 L/h]
  TCA_gut → TCA_cent               [1-cpt, CL 40 L/h, Vd 1500 L]

Peripheral (3 CMT)
  NGF, PGE2, DRG_act

Spinal/Central (7 CMT)
  SP_csf, NMDA_state, WindUp, LTP_cs, NE_syn, SHT_syn, DPMS

HPA Axis (3 CMT)
  CRH → ACTH → CORT (negative feedback loop)

ANS + Sleep (4 CMT)
  SNS_tone, SWS_depth, Adenosine, DPMS

Neuroinflammation (2 CMT)
  MG_act (microglia), IL1b_sp

Clinical Outcomes (4 CMT)
  Pain_score (NRS 0-10), FIQ_score (FIQR 0-100),
  Fatigue_VAS (0-100), Depression_score (PHQ-9 0-27)
```

---

## Drug PK Parameters

| Drug | Mechanism | Vd (L) | CL (L/h) | t½ (h) | F (%) |
|------|-----------|--------|----------|--------|-------|
| Duloxetine | SERT+NET inhibition | 1640 | 54 | ~12 | 50 |
| Pregabalin | α2δ-1 Ca-channel block | 42 | 6.8 (renal) | ~6 | 90 |
| Milnacipran | SERT+NET inhibition | 300 | 50 | ~8 | 85 |
| Amitriptyline | SERT/NET/H1/muscarinic | 1500 | 40 | ~25 | 48 |

---

## Treatment Scenarios Simulated

| # | Scenario | Drug(s) | Dose |
|---|----------|---------|------|
| 1 | Untreated FM baseline | — | — |
| 2 | Duloxetine monotherapy | DUL | 60 mg QD |
| 3 | Pregabalin monotherapy | PRE | 150 mg BID |
| 4 | Milnacipran monotherapy | MIL | 50 mg BID |
| 5 | DUL + PRE combination | DUL + PRE | 60 mg QD + 150 mg BID |
| 6 | Low-dose amitriptyline (sleep) | TCA | 25 mg QHS |

---

## Shiny App Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease overview, diagnostic criteria, QSP structure, baseline KPIs |
| 2. Drug PK | Plasma concentration-time profiles, PK parameter table, target occupancy |
| 3. PD Key Markers | Spinal LTP, DPMS, NE/5-HT pools, wind-up & NMDA state |
| 4. Clinical Endpoints | Pain NRS, FIQR, Fatigue VAS, PHQ-9; week-12 outcome table |
| 5. Scenario Comparison | All 6 scenarios overlaid; responder analysis (≥30%/≥50%) |
| 6. Biomarkers | CSF Substance P, microglia, cortisol, SWS depth, SNS tone, IL-1β |

---

## Key Findings from Model

1. **Central sensitization (LTP index)** is the primary driver of chronic pain in FM — all effective drugs reduce it via complementary mechanisms.
2. **DUL + PRE combination** shows additive effects via distinct targets (SERT/NET vs α2δ-1).
3. **Amitriptyline** uniquely improves slow-wave sleep (via H1 blockade) reducing the sleep-pain vicious cycle.
4. **Descending inhibition deficit (DPMS)** explains why standard analgesics (NSAIDs, opioids) are ineffective.
5. **HPA axis blunting** (hypocortisolism) removes anti-inflammatory brake → perpetuates neuroinflammation.
6. **CSF Substance P** tracks central sensitization and predicts treatment response.

---

## Running the Model

```r
# Install dependencies
install.packages(c("mrgsolve","ggplot2","dplyr","tidyr","gridExtra"))

# Run mrgsolve simulation
source("fm_mrgsolve_model.R")

# Launch interactive Shiny app
shiny::runApp("fm_shiny_app.R")
```

---

## References

See [`fm_references.md`](fm_references.md) for 60 annotated PubMed references covering:
- Pathophysiology & central sensitization
- Neuroendocrinology & autonomic dysfunction
- Sleep disorders
- Neuroimaging & biomarkers
- Drug PK/PD for duloxetine, pregabalin, milnacipran, amitriptyline
- Clinical trials & outcome measures

---

## Citation

*Generated by Claude Code Routine (CCR) | 2026-06-23*  
*Fibromyalgia QSP Model v1.0 — for research and educational purposes only*
