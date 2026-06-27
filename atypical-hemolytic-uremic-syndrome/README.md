# Atypical Hemolytic Uremic Syndrome (aHUS) QSP Model

> **Complement Alternative Pathway Dysregulation → Thrombotic Microangiopathy → Multi-Organ Damage**

[![Disease](https://img.shields.io/badge/disease-aHUS-red)](https://www.ncbi.nlm.nih.gov/books/NBK544267/)
[![Framework](https://img.shields.io/badge/mrgsolve-18%20ODEs-blue)](https://mrgsolve.org)
[![Shiny](https://img.shields.io/badge/Shiny-8%20tabs-green)](https://shiny.posit.co)
[![Drugs](https://img.shields.io/badge/drugs-eculizumab%20%7C%20ravulizumab%20%7C%20iptacopan%20%7C%20danicopan%20%7C%20avacopan-purple)](https://www.accessdata.fda.gov/drugsatfda_docs/label/2011/125166s219lbl.pdf)

---

## Disease Overview

**Atypical Hemolytic Uremic Syndrome (aHUS)** is a rare, life-threatening **thrombotic microangiopathy (TMA)** driven by chronic, uncontrolled activation of the **complement alternative pathway (AP)**. Unlike typical HUS (caused by Shiga toxin-producing *E. coli*), aHUS is caused by:

- **Germline mutations** in complement regulatory genes: *CFH* (~25%), *CD46* (~12%), *CFI* (~8%), *C3* (~8%), *CFB* (~2%), *THBD* (~3%)
- **Anti-CFH autoantibodies** (~10% of cases)
- **Incomplete penetrance** (~50% in *CFH* carriers) — environmental triggers required

### Pathophysiological Cascade

```
AP Dysregulation (CFH/CFI mutation)
        ↓
C3 Tickover → C3b accumulation (uncontrolled)
        ↓
C5 Convertase (C3b₂Bb) → C5 cleavage
        ↓
C5a (anaphylatoxin) + C5b → MAC (C5b-9)
        ↓
Endothelial injury → TMA
        ↓
Thrombocytopenia + MAHA + Acute Kidney Injury
```

### Epidemiology
- Incidence: ~0.5–2 per million/year
- Median age at onset: 2.5 years (children) / 33 years (adults)
- ESRD risk: 50% within first year if untreated
- Mortality: 25% in first acute episode
- Relapse rate: 40–80% with CFH mutations

---

## QSP Model Architecture

### Mechanistic Map

| File | Description |
|------|-------------|
| [`ahus_qsp_model.dot`](ahus_qsp_model.dot) | Graphviz source (13 clusters, 130+ nodes) |
| [`ahus_qsp_model.svg`](ahus_qsp_model.svg) | Scalable vector diagram |
| [`ahus_qsp_model.png`](ahus_qsp_model.png) | High-resolution raster (150 dpi) |

**Clusters included:**
1. Genetic & Acquired Risk Factors
2. Complement Alternative Pathway (AP)
3. Classical & Lectin Pathways (CP/LP)
4. Complement Regulatory Proteins (CFH, CFI, CD46, CD59...)
5. Terminal Complement (C5–C9 → MAC)
6. Endothelial Injury & TMA
7. Platelet & Hematological Compartment
8. Renal Pathophysiology
9. Extrarenal Manifestations
10. Inflammatory & Coagulation Amplification
11. Drug Mechanisms (PK/PD)
12. Clinical Endpoints & Biomarkers
13. Clinical Monitoring & Management

### ODE Model — 18 State Variables

| # | Compartment | Variable | Unit | Notes |
|---|-------------|----------|------|-------|
| 1 | Drug central | `Drug_C` | mg | Eculizumab PK (IV) |
| 2 | Drug peripheral | `Drug_P` | mg | 2-comp PK |
| 3 | Drug-C5 complex | `Drug_C5` | mg | TMDD complex |
| 4 | C3 pool | `C3pool` | nM | Plasma C3 (~8000 nM) |
| 5 | AP C3b | `C3b_AP` | nM | Active C3b accumulation |
| 6 | Free C5 | `C5free` | nM | Key drug target (~395 nM) |
| 7 | C5 convertase | `C5conv` | AU | C3b₂Bb activity |
| 8 | MAC flux | `MACflux` | AU | C5b-9 formation rate |
| 9 | Endothelial injury | `Endo_inj` | 0–10 | TMA driver score |
| 10 | Platelet count | `PLT` | ×10⁹/L | Thrombocytopenia marker |
| 11 | Hemoglobin | `Hgb` | g/dL | MAHA marker |
| 12 | LDH | `LDH` | U/L | Hemolysis marker |
| 13 | Haptoglobin | `Hpg` | g/L | Free Hgb scavenger |
| 14 | eGFR | `GFR` | mL/min/1.73m² | Renal function |
| 15 | Schistocytes | `Schist` | % | Blood smear |
| 16 | CRP | `CRP` | mg/L | Inflammatory marker |
| 17 | sC5b-9 | `sC5b9` | ng/mL | Complement activation biomarker |
| 18 | CH50 | `CH50pct` | % | Functional complement activity |

### Drug PK Parameters (Eculizumab, 2-Compartment)

| Parameter | Value | Unit | Reference |
|-----------|-------|------|-----------|
| CL | 0.31 | L/day | Menne 2015 |
| V1 | 5.30 | L | Menne 2015 |
| Q | 0.54 | L/day | Menne 2015 |
| V2 | 4.10 | L | Menne 2015 |
| t½ | ~11 | days | Calculated |
| C5 Kd | ~10 | pM | Thomas 1996 |
| MW | 148 | kDa | Product info |

### Treatment Scenarios Simulated

| Scenario | Treatment | Dosing | Route |
|----------|-----------|--------|-------|
| 1 | Natural history | — | — |
| 2 | **Eculizumab** (standard) | 900 mg/wk×4 → 1200 mg/q2w | IV |
| 3 | **Ravulizumab** | 2400/3000 mg → 3300 mg/q8w | IV |
| 4 | **Eculizumab + Danicopan** (fDi) | Standard + 150 mg TID oral add-on | IV+PO |
| 5 | **Iptacopan** (oral fBi) | 200 mg BID | PO |

---

## Key Model Equations

### Alternative Pathway Activation
```
dC3b/dt = kAP × C3/C3_ss + kAP_amp × (C3b/C3_ss) × C3/C3_ss 
          - kC3b_deg × CFH_factor × CFI_factor × C3b
```

### TMDD for Eculizumab-C5 Binding
```
dC5free/dt = kC5syn - kC5deg × C5free - kC5conv × C5conv × C5free
             - (kon × Cc_nM × C5free - koff × [Drug-C5])
```

### Endothelial Injury Dynamics
```
dEI/dt = kEI_MAC × MACflux - kEI_rep × EI
```

### GFR Recovery (Complement-dependent)
```
dGFR/dt = kGFR_rep × (GFR_ss - GFR) × (1 - EI/EI_max)
          - kGFR_loss × GFR × (EI/6)
```

---

## Shiny App — 8 Interactive Tabs

| Tab | Content |
|-----|---------|
| **Overview** | Disease overview, mechanistic cascade, KPI value boxes |
| **Patient Profile** | Complement regulatory radar, presentation values, cascade dynamics |
| **PK Analysis** | Eculizumab concentration-time, TMDD (C5 free vs bound), blockade % |
| **Complement PD** | C3 pool, C3b AP, C5 convertase, MAC flux |
| **TMA & Hematology** | Platelets, Hgb, LDH, haptoglobin, schistocytes |
| **Renal Outcomes** | eGFR trajectory, endothelial injury, dialysis risk, CKD staging |
| **Scenario Compare** | All 5 treatment arms overlaid, efficacy summary table |
| **Biomarkers** | sC5b-9, CH50, CRP, biomarker correlation dashboard |

---

## Running the Model

### Requirements

```r
# Install required packages
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork",
                   "shiny","shinydashboard","plotly","DT"))
```

### mrgsolve Simulation

```r
source("ahus_mrgsolve_model.R")
# → runs 5 treatment scenarios × 365 days
# → generates multi-panel comparison plot
# → prints efficacy summary table at days 28, 84, 180, 365
```

### Shiny Dashboard

```r
shiny::runApp("ahus_shiny_app.R")
# → Interactive 8-tab dashboard
# → Parameter sliders: CFH/CFI/CD46 function, AP activation rate, kEI_MAC
# → Scenario selector: 5 treatment arms
# → All plots in plotly (interactive, zoom, download)
```

---

## Key Clinical Predictions

### With Eculizumab (Standard Dosing):
- **Day 7**: CH50 → <10% (complement blockade confirmed)
- **Day 28**: PLT >150 ×10⁹/L in ~80% of patients
- **Day 84**: eGFR improvement of 25–35 mL/min
- **Day 365**: Dialysis-free rate ~85% vs ~30% untreated

### Natural History (untreated):
- Progressive TMA → ESRD within 3–6 months
- 25% mortality in first acute episode
- Persistent dialysis requirement in >50%

### Iptacopan Advantage:
- Proximal complement inhibition → also addresses C3-mediated hemolysis
- Oral bioavailability eliminates IV infusion burden
- Covers C3 fragment deposition (missed by anti-C5 agents)

---

## Complement Biomarker Reference Ranges

| Biomarker | Normal | Active aHUS | On Eculizumab |
|-----------|--------|-------------|---------------|
| CH50 | 60–144 U/mL (75–100%) | ↓↓ (consumed) | <10% |
| sC5b-9 | <244 ng/mL | >2000 ng/mL | <244 ng/mL |
| C3 | 1.0–1.8 g/L | ↓ (consumed) | Normalizes |
| C4 | 0.2–0.4 g/L | Normal (AP) | Normal |
| PLT | 150–400 ×10⁹/L | <150 | Recovers |
| LDH | <ULN (200 U/L) | >500 U/L | Normalizes |
| Haptoglobin | 0.3–2.0 g/L | <0.1 g/L | Normalizes |

---

## Files in This Directory

```
atypical-hemolytic-uremic-syndrome/
├── ahus_qsp_model.dot       ← Graphviz mechanistic map source
├── ahus_qsp_model.svg       ← Scalable diagram (13 clusters, 130+ nodes)
├── ahus_qsp_model.png       ← High-res PNG thumbnail (150 dpi)
├── ahus_mrgsolve_model.R    ← 18-ODE mrgsolve model, 5 scenarios
├── ahus_shiny_app.R         ← Interactive Shiny dashboard (8 tabs)
├── ahus_references.md       ← 55 PubMed references (13 sections)
└── README.md                ← This file
```

---

## References (Selected)

1. Fakhouri F et al. Lancet 2017;390:1847–1860. [PMID: 28499565](https://pubmed.ncbi.nlm.nih.gov/28499565/)
2. Legendre CM et al. NEJM 2013;368:2169–2181. [PMID: 23738544](https://pubmed.ncbi.nlm.nih.gov/23738544/)
3. Jokiranta TS. Blood 2017;129:2847–2856. [PMID: 28416507](https://pubmed.ncbi.nlm.nih.gov/28416507/)
4. Rother RP et al. Nat Biotechnol 2007;25:1256–1264. [PMID: 17989688](https://pubmed.ncbi.nlm.nih.gov/17989688/)
5. Goodship TH et al. (KDIGO). Kidney Int 2017;91:539–551. [PMID: 28131400](https://pubmed.ncbi.nlm.nih.gov/28131400/)

*See [ahus_references.md](ahus_references.md) for all 55 references.*

---

*Built by Claude Code Routine | 2026-06-27 | QSP Disease Model Library*
