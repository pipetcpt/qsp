# Paroxysmal Nocturnal Hemoglobinuria (PNH) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-PNH-red)](.)
[![Nodes](https://img.shields.io/badge/Map_Nodes-130%2B-blue)](pnh_qsp_model.dot)
[![ODE](https://img.shields.io/badge/ODE_Compartments-24-green)](pnh_mrgsolve_model.R)
[![Scenarios](https://img.shields.io/badge/Scenarios-6-orange)](pnh_mrgsolve_model.R)
[![Tabs](https://img.shields.io/badge/Shiny_Tabs-8-purple)](pnh_shiny_app.R)
[![Refs](https://img.shields.io/badge/References-35-teal)](pnh_references.md)

---

## Mechanistic Map

[![PNH QSP Model](pnh_qsp_model.png)](pnh_qsp_model.svg)

*Click image to view interactive SVG. 130+ nodes, 13 clusters covering GPI anchor biosynthesis, complement alternative pathway, terminal complement/MAC, hemolysis, NO biology, thrombosis, drug PK/PD, and clinical endpoints.*

---

## Disease Summary

**Paroxysmal Nocturnal Hemoglobinuria (PNH)** is a rare acquired clonal disorder of hematopoietic stem cells caused by somatic mutation of the *PIGA* gene (X-linked). PIGA encodes a critical enzyme in the first step of glycosylphosphatidylinositol (GPI) anchor biosynthesis, and its loss results in complete absence of all GPI-anchored proteins from affected cells.

### Core Pathophysiology

```
PIGA mutation (somatic)
       ↓
GPI anchor deficiency on blood cells (PNH clone: type I, II, III)
       ↓
Loss of CD55 (DAF) and CD59 (MIRL) from RBC/Plt surfaces
       ↓
   ┌──────────────────────────────────────────────────┐
   │                                                  │
No CD55                                          No CD59
(DAF absent)                               (MIRL absent)
   ↓                                               ↓
C3bBb convertase                           C9 polymerizes
not decayed                                    freely
   ↓                                               ↓
Massive C3b                               MAC (C5b-9)
deposition on PNH RBC                     pore formation
   ↓                                               ↓
Extravascular hemolysis                  Intravascular
(C3b/iC3b → spleen/liver)               hemolysis (IVH)
   ↓                                               ↓
   └──────────────────────┬───────────────────────┘
                          ↓
              Free plasma hemoglobin
              ↓               ↓
         NO scavenging    Haptoglobin depletion
              ↓
         NO deficiency
         ↓           ↓
    Thrombosis   Smooth muscle dystonias
    (DVT, Budd-  (esophageal spasm,
    Chiari, CVST)  erectile dysfunction)
```

### Key Clinical Features
| Feature | Mechanism | Prevalence |
|---------|-----------|------------|
| Chronic hemolytic anemia | MAC + EVH | ~100% |
| Hemoglobinuria | IVH → free Hgb in urine | Classic but not universal |
| Fatigue | Anemia + NO depletion | ~90% |
| Thromboembolism | NO↓, platelet activation | 30–40% lifetime |
| Budd-Chiari syndrome | Hepatic vein thrombosis | 7–10% |
| Aplastic anemia overlap | Shared autoimmune HSPC destruction | 30–40% AA have PNH clone |

---

## Model Architecture

### ODE Compartments (24 total)

| Category | Compartments | Description |
|----------|-------------|-------------|
| **Hematopoiesis** | PNH_Ret, NL_Ret, PNH_RBC, NL_RBC | EPO-driven production, MAC/C3b lysis |
| **Complement** | C3, C3b, C5, MAC | Alternative pathway, amplification, MAC |
| **Hemolysis outputs** | fHgb, Haptoglobin, LDH, NO_rel | Free Hgb, scavengers, biomarkers |
| **Eculizumab PK** | ECU_C, ECU_P, C5_ECU | 2-compartment IV model |
| **Ravulizumab PK** | RAV_C, RAV_P, C5_RAV | 2-compartment IV (longer t½) |
| **Iptacopan PK** | IPC_gut, IPC_plasma | 1-compartment oral (Factor B) |
| **Danicopan PK** | DAN_gut, DAN_plasma | 1-compartment oral (Factor D) |

### Key Parameters Calibrated Against Clinical Data

| Parameter | Value | Source |
|-----------|-------|--------|
| Eculizumab CL | 0.31 L/day | Hills 2019 PopPK |
| Eculizumab t½ | ~11 days | TRIUMPH PK substudy |
| Ravulizumab CL | 0.069 L/day | Lee 2019 PopPK |
| Ravulizumab t½ | ~49 days | ALXN1210-301 |
| Iptacopan F | 69% | APPLY-PNH Phase 3 |
| Iptacopan IC50 (FB) | 0.05 μg/mL | Schubart 2019 |
| Normal RBC lifespan | 120 days | Standard hematology |
| PNH RBC lifespan | ~12 days | Complement-shortened |

---

## Treatment Scenarios

| Scenario | Drugs | Mechanism | Key Trial | TI Rate |
|----------|-------|-----------|-----------|---------|
| **S0: Untreated** | None | Natural history | Brodsky 2014 | 0% |
| **S1: Eculizumab** | ECU 900mg q2w IV | C5 → blocks IVH | TRIUMPH NEJM 2006 | 49% |
| **S2: Ravulizumab** | RAV 3300mg q8w IV | C5 (longer t½) → same | ALXN1210-301 Blood 2019 | 73.6% |
| **S3: Iptacopan** | IPC 200mg BID PO | Factor B → blocks IVH + EVH | APPLY NEJM 2024 | 51.1% |
| **S4: ECU + Danicopan** | ECU + DAN 150mg TID | C5 + FD → reduces EVH | GALAXY trial | IPC > ECU for Hgb |
| **S5: Iptacopan (high clone)** | IPC 200mg BID | Factor B (f_PNH=0.85) | APPLY high-clone subset | — |

### Why Iptacopan Outperforms Eculizumab for TI?
```
Eculizumab blocks C5:
  ✅ Prevents MAC → eliminates IVH
  ❌ C3b still deposits on PNH RBCs
  ❌ C3b-mediated EVH continues (spleen/liver)
  → Residual anemia persists despite LDH normalization

Iptacopan blocks Factor B:
  ✅ Blocks C3 convertase amplification
  ✅ Prevents C3b deposition → eliminates EVH
  ✅ Prevents MAC → eliminates IVH
  → Complete suppression of both hemolysis pathways
  → Superior Hgb improvement and TI rates
```

---

## File Index

| File | Description |
|------|-------------|
| [`pnh_qsp_model.dot`](pnh_qsp_model.dot) | Graphviz mechanistic map source (130+ nodes, 13 clusters) |
| [`pnh_qsp_model.svg`](pnh_qsp_model.svg) | Interactive SVG |
| [`pnh_qsp_model.png`](pnh_qsp_model.png) | Static PNG (150 dpi) |
| [`pnh_mrgsolve_model.R`](pnh_mrgsolve_model.R) | ODE model: 24 compartments, 6 scenarios, full PK/PD |
| [`pnh_shiny_app.R`](pnh_shiny_app.R) | Shiny dashboard: 8 tabs |
| [`pnh_references.md`](pnh_references.md) | 35 PubMed references (12 sections) |

---

## Shiny Dashboard Tabs (8 tabs)

| Tab | Contents |
|-----|----------|
| 1. Patient Profile | ValueBoxes, RBC dynamics, lab table |
| 2. Drug PK | C5i trough levels, oral drug levels, PK summary table |
| 3. Complement | C3/C5 dynamics, C3b/MAC evolution, pathway overview |
| 4. Hemolysis Markers | LDH, free Hgb, haptoglobin, NO levels |
| 5. Clinical Endpoints | Hemoglobin, TI threshold, FACIT-fatigue, trial comparison table |
| 6. Scenario Comparison | 4-way scenario comparison: Hgb, LDH, C3b, FACIT |
| 7. Biomarkers | Clone size, free C5 monitoring, haptoglobin, reference ranges |
| 8. About | Model description, architecture, references |

---

## Quick Start

```r
# Run ODE simulations
source("pnh_mrgsolve_model.R")

# Launch Shiny app
shiny::runApp("pnh_shiny_app.R")

# Render mechanistic map
# dot -Tsvg pnh_qsp_model.dot -o pnh_qsp_model.svg
# dot -Tpng -Gdpi=150 pnh_qsp_model.dot -o pnh_qsp_model.png
```

**Required R packages:**
```r
install.packages(c("mrgsolve","shiny","shinydashboard","dplyr","ggplot2","tidyr","DT"))
```

---

## Key Insights from the Model

1. **C5 inhibitors normalize LDH but not always Hgb:** Eculizumab/ravulizumab stop intravascular hemolysis (MAC-mediated) but cannot prevent extravascular hemolysis driven by C3b opsonization — explaining why ~50% of eculizumab-treated patients remain transfusion-dependent.

2. **Proximal complement inhibition is mechanistically superior:** Iptacopan (Factor B) blocks before C3 cleavage, preventing both C3b deposition (EVH) and MAC formation (IVH), enabling transfusion independence in 51.1% vs. 0% for eculizumab (APPLY-PNH, NEJM 2024).

3. **NO depletion drives thrombotic and vasomotor complications:** Free hemoglobin scavenges nitric oxide, explaining smooth muscle dystonias (esophageal spasm, erectile dysfunction) and elevated thrombosis risk independent of the underlying coagulation cascade.

4. **PNH clone persists on all treatments:** Current complement inhibitors suppress hemolysis but do not eradicate the PNH clone. Only allogeneic HSCT is curative. This underlies the need for lifelong therapy.

5. **Meningococcal infection risk:** C5 inhibition abrogates MAC-dependent killing of Neisseria meningitidis, increasing infection risk ~1000-fold — mandatory meningococcal vaccination and prophylaxis are required.

---

*Generated 2026-06-24 | QSP Disease Model Library | Claude Code Routine (CCR)*
