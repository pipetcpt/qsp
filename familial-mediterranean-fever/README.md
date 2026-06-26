# Familial Mediterranean Fever (FMF) — QSP Model

[![Mechanistic Map](fmf_qsp_model.png)](fmf_qsp_model.svg)

---

## Disease Overview

**Familial Mediterranean Fever (FMF)** is the most common hereditary autoinflammatory disease, caused by gain-of-function mutations in the *MEFV* gene encoding **PYRIN** (marenostrin). It is characterized by episodic self-limited attacks of fever and serositis (peritonitis, pleuritis, synovitis), and carries a major long-term risk of **AA amyloidosis** leading to renal failure.

| Feature | Details |
|---------|---------|
| Gene | *MEFV* (chromosome 16p13.3) |
| Protein | PYRIN / Marenostrin |
| Inheritance | Autosomal recessive (AR); some AD-like phenotypes |
| Prevalence | 1:200–1:1000 in Mediterranean populations |
| Key mutations | M694V, M680I, V726A, E148Q, M694I |
| Onset | Usually < 20 years |
| Episodes | Fever 38–40°C + serositis, 12–72h duration |
| Major complication | AA amyloidosis → CKD/ESRD |

---

## Molecular Pathophysiology

```
MEFV mutation (M694V/M680I)
    │
    ▼
PYRIN protein: impaired RhoA/PKN1-PKN2–mediated phosphorylation
    │         → cannot be inactivated → constitutively active
    ▼
ASC speck assembly (Ser208/Ser242 phosphorylation blocked)
    │
    ▼
Caspase-1 activation
    │
    ├──► Pro-IL-1β cleavage → mature IL-1β (17 kDa) secretion
    ├──► Pro-IL-18 cleavage → IL-18 secretion
    └──► Gasdermin D (GSDMD) cleavage → pyroptosis pores → sterile inflammation
    
IL-1β / IL-18 ──► fever, acute phase response (SAA↑, CRP↑), serositis, neutrophil influx
SAA excess ──► AA amyloid fibril deposition in kidney, liver, spleen
AA amyloidosis ──► nephrotic syndrome → CKD → ESRD
```

---

## Model Architecture

### Mechanistic Map (`fmf_qsp_model.dot/.svg/.png`)

**10 subgraph clusters:**
1. **Genetic Basis** — MEFV alleles, mutation severity scores
2. **PYRIN Inflammasome** — RhoA-GTP, phospho-PYRIN, ASC speck, Caspase-1, GSDMD
3. **Innate Immunity** — neutrophil recruitment, NF-κB, NLRP3 crosstalk
4. **Cytokines & Acute Phase** — IL-1β, IL-18, IL-6, TNF-α, SAA, CRP
5. **Acute Attack Manifestations** — peritonitis, pleuritis, synovitis, fever, erythema
6. **Long-term Complications** — AA amyloidosis, renal amyloidosis, CKD, ESRD
7. **Biomarkers & Clinical Endpoints** — AIDAI, PGA, SAA, CRP, eGFR
8. **Colchicine PK** — oral absorption, 2-compartment, leukocyte accumulation
9. **IL-1 Inhibitor PK** — Anakinra SC, Canakinumab 2-compartment SC
10. **Drug PD** — colchicine tubulin/PYRIN inhibition, anakinra IL-1R blockade, canakinumab IL-1β neutralization

---

### mrgsolve ODE Model (`fmf_mrgsolve_model.R`)

**22 ODE compartments:**

| Module | Compartments |
|--------|-------------|
| Colchicine PK | GUT → CENTRAL → PERIPHERAL + LEUKOCYTE (4) |
| Anakinra PK | SC depot → CENTRAL (2) |
| Canakinumab PK | SC depot → CENTRAL → PERIPHERAL (3) |
| PYRIN inflammasome | RhoA, Pyrin_phospho, ASC_speck, Caspase1 (4) |
| Cytokines | IL-1β_pro, IL-1β_mature, IL-18, SAA, CRP (5) |
| Neutrophils | Circulating ANC, tissue neutrophils (2) |
| Attack dynamics | Trigger signal, severity score (2) |
| Amyloidosis | AA deposits, eGFR (2) |

**5 Treatment scenarios:**
1. No treatment (untreated M694V phenotype)
2. Colchicine 0.5 mg BID (standard first-line)
3. Colchicine 1.0 mg QD (alternative dosing)
4. Anakinra 100 mg SC QD (IL-1Ra, colchicine-resistant)
5. Canakinumab 150 mg SC Q8W (anti-IL-1β mAb, CLUSTER trial)

**Key PK parameters (colchicine):**
- F = 45%; Vd ~5 L/kg; t½ ~30h; leukocyte:plasma ratio ~10:1
- Mechanism: tubulin polymerization inhibition → neutrophil immobility; PYRIN/NLRP3 inhibition

---

### Shiny Dashboard (`fmf_shiny_app.R`)

**8 interactive tabs:**

| Tab | Content |
|-----|---------|
| 1. Patient Profile | MEFV genotype, labs, drug selection, disease classification |
| 2. Drug PK | Colchicine plasma/leukocyte, Anakinra, Canakinumab PK curves |
| 3. Inflammasome Dynamics | PYRIN, ASC, Caspase-1, IL-1β under drug effect |
| 4. Attack Simulation | Attack frequency/severity time series, statistics |
| 5. Clinical Endpoints | SAA, CRP, AIDAI, eGFR, outcome boxes |
| 6. Scenario Comparison | Multi-arm comparison table + bar chart |
| 7. Amyloidosis Risk | Long-term SAA, AA deposits, eGFR projection (≤20 years) |
| 8. Sensitivity Analysis | One-at-a-time parameter sensitivity |

---

## Key Clinical Trials Informing Parameters

| Trial | Drug | N | Key Finding |
|-------|------|---|-------------|
| Zemer et al. 1986 NEJM | Colchicine | 470 | 75% attack reduction; amyloidosis prevention |
| De Benedetti et al. 2018 NEJM (CLUSTER) | Canakinumab | 63 (FMF) | Attack-free 61% vs 6% placebo; SAA normalization |
| Georgin-Lavialle et al. 2020 | Anakinra | 67 | Rapid attack resolution; SAA reduction |
| Varan et al. 2019 | IL-1 inhibitors | 72 | Amyloidosis regression with IL-1 blockade |

---

## Drug Mechanisms

| Drug | Target | Mechanism | Route | Dosing |
|------|--------|-----------|-------|--------|
| **Colchicine** | Tubulin / PYRIN | Blocks neutrophil microtubule polymerization; reduces PYRIN activity | Oral | 0.5–1 mg QD/BID |
| **Anakinra** | IL-1 receptor (IL-1R1) | Competitive antagonist; blocks both IL-1α and IL-1β | SC daily | 100 mg/day |
| **Canakinumab** | IL-1β (free) | Neutralizing anti-IL-1β monoclonal antibody; t½ ~26d | SC Q8W | 150 mg |
| **Rilonacept** | IL-1α + IL-1β | Dimeric IL-1 trap (ILRAP-Fc) | SC weekly | 160 mg |

---

## Files

| File | Description |
|------|-------------|
| `fmf_qsp_model.dot` | Graphviz source (fdp layout, 100+ nodes, 10 clusters) |
| `fmf_qsp_model.svg` | Vector mechanistic map |
| `fmf_qsp_model.png` | Raster thumbnail (150 dpi) |
| `fmf_mrgsolve_model.R` | mrgsolve ODE model (22 compartments, 5 scenarios) |
| `fmf_shiny_app.R` | Shiny dashboard (8 tabs) |
| `fmf_references.md` | 50 curated PubMed references |

---

## Quick Start

```r
# Install dependencies (first time)
install.packages(c("mrgsolve", "dplyr", "tidyr", "ggplot2", "patchwork",
                   "shiny", "DT"))

# Run mrgsolve simulation
source("fmf_mrgsolve_model.R")
print(summary_metrics)
combined_plot

# Launch Shiny app
shiny::runApp("fmf_shiny_app.R")
```

---

## References

See [`fmf_references.md`](fmf_references.md) — 50 citations organized by topic:
genetics · PYRIN inflammasome · IL-1β biology · colchicine · IL-1 inhibitors · amyloidosis · classification criteria · QSP modeling

---

*Added 2026-06-25 | Claude Code Routine (CCR) | QSP Model Library*
