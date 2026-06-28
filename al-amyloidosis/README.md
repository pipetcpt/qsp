# AL (Immunoglobulin Light Chain) Amyloidosis — QSP Model

[![Mechanistic Map](al_qsp_model.png)](al_qsp_model.svg)

## Disease Overview

**AL amyloidosis** (immunoglobulin light chain amyloidosis) is a systemic disease caused by a clonal plasma cell dyscrasia in which misfolded monoclonal free light chains (FLC) polymerize into amyloid fibrils and deposit in organs—primarily the heart, kidneys, liver, and peripheral nerves—causing progressive multi-organ failure.

| Feature | Details |
|---------|---------|
| **Incidence** | ~10 per million/year (~4,000 new US cases/year) |
| **Median Age at Diagnosis** | ~63 years |
| **Light Chain Type** | λ (lambda) more common than κ (kappa) in AL |
| **Most Affected Organs** | Heart (>75%), Kidney (>50%), Liver (~20%), Nerves (~15%) |
| **Most Common Cause of Death** | Cardiac arrhythmia / sudden cardiac death |
| **Key Prognostic Factor** | Mayo 2012 Stage (NT-proBNP + hs-TnT + dFLC + eGFR) |

## Mechanistic Map Clusters

| Cluster | Description |
|---------|-------------|
| ① Plasma Cell Biology | HSC → MGUS → SMM → clonal PC expansion, CD38 expression |
| ② Light Chain Production | Monoclonal FLC (λ/κ), dFLC, misfolding, UPR, oligomers |
| ③ Amyloid Fibril Formation | Nucleation, elongation, seeding (autocatalytic), SAP/GAG binding |
| ④ Cardiac Amyloidosis | LV thickening, diastolic dysfunction, NT-proBNP, hs-TnT, NYHA HF |
| ⑤ Renal Amyloidosis | Glomerular amyloid, proteinuria, eGFR decline, nephrotic syndrome |
| ⑥ Hepatic & Neuropathic | Hepatomegaly, ALP, autonomic/peripheral neuropathy, macroglossia |
| ⑦ BM Microenvironment | Stromal IL-6, APRIL/BAFF, NF-κB, IRF4/BLIMP-1, MCL-1 |
| ⑧ Daratumumab PK | TMDD 2-CMT model, CD38 receptor dynamics, Q1W→Q2W→Q4W schedule |
| ⑨ CyBorD PK | Bortezomib 1-CMT + proteasome binding; CY active metabolite; DEX |
| ⑩ Pharmacodynamics | ADCC/CDC/ADCP, PC killing, hematologic + organ response |
| ⑪ Staging & Biomarkers | Mayo 2012 (4-stage), RF thresholds, ECOG PS |
| ⑫ Immune Effectors | NK cells, CD8 T cells, Treg, PD-1/PD-L1, NK fratricide by Dara |
| ⑬ Clinical Endpoints | OS, EFS, PFS, TTR, cardiac death, QoL, hospitalization |
| ⑭ Safety Monitor | IRR, cytopenias, neuropathy, DVT/PE, hyperglycemia, REMS |

## Drug Summary

| Drug | Class | Target | Key Dose |
|------|-------|--------|---------|
| Daratumumab | Anti-CD38 mAb | CD38 on PCs | 16 mg/kg IV QW×8→Q2W×16→Q4W |
| Bortezomib | Proteasome inhibitor | 20S proteasome | 1.3 mg/m² SC D1,8,15,22 |
| Cyclophosphamide | Alkylating agent | DNA (via 4-OH-CY) | 300 mg/m² PO D1,8,15,22 |
| Dexamethasone | Glucocorticoid | GR-mediated apoptosis | 20-40 mg PO D1,8,15,22 |
| Melphalan | Alkylating agent | DNA cross-linking | 0.15-0.25 mg/kg PO D1-4 (maintenance) |

## Model Specifications

### mrgsolve ODE Model (`al_mrgsolve_model.R`)
- **20 ODE compartments** across 2 modules:
  - **PK (12 compartments):** Daratumumab C1, C2, RC (TMDD complex), CD38 receptor; BTZ plasma + proteasome complex; CY gut depot + active metabolite; DEX gut depot + plasma; MEL gut + plasma
  - **Disease PD (8 compartments):** Plasma cell pool, dFLC, cardiac amyloid, renal amyloid, NT-proBNP, hs-TnT, eGFR, NK cells
- **7 treatment scenarios:**
  1. Untreated (natural history)
  2. Daratumumab monotherapy
  3. CyBorD (3-drug)
  4. Dara-CyBorD (ANDROMEDA regimen)
  5. VCD (bortezomib + CY + high-dose DEX)
  6. Melphalan + Dexamethasone (MDex)
  7. Dara-CyBorD in CYP2C19 poor metabolizer

### Shiny App (`al_shiny_app.R`)
8 interactive tabs:
1. **Patient Profile** — Mayo 2012 staging, disease overview, drug mechanisms
2. **PK Profiles** — Daratumumab (TMDD), bortezomib-proteasome complex, CY, DEX
3. **Hematologic PD** — Plasma cell pool, dFLC, NK cells, response classification
4. **Organ Biomarkers** — NT-proBNP, hs-TnT, eGFR, 24h proteinuria
5. **Clinical Endpoints** — Cardiac/renal amyloid burden, organ response flags
6. **Scenario Comparison** — Multi-arm dFLC/NT-proBNP/eGFR comparison + Day 180 table
7. **Biomarker Staging** — Mayo 2012 stage evolution, risk factor status gauges
8. **Safety Monitor** — IRR, NK depletion, DDI table, REMS schedule

### References (`al_references.md`)
50 references across 10 sections: Landmark Trials · Pathophysiology · Cardiac Amyloidosis · Biomarkers/Staging · Daratumumab PK/PD · Bortezomib · Renal Amyloidosis · QSP Modeling · Guidelines · Drug Combination

## Calibration

| Clinical Trial | Regimen | CR Rate | Cardiac Organ Response | Reference |
|---------------|---------|---------|----------------------|-----------|
| **ANDROMEDA** | Dara+CyBorD | **53.3%** | ~42% | NEJM 2021 |
| **ANDROMEDA** | CyBorD alone | 18.1% | ~22% | NEJM 2021 |
| MDex | Melphalan+Dex | ~33% | ~26% | Blood 2004 |
| CyBorD | Bortezomib+CY+DEX | ~29% | ~18% | Blood 2012 |

## Response Criteria

| Category | Hematologic | Cardiac Organ | Renal Organ |
|----------|-------------|---------------|-------------|
| **Response** | dFLC ≥50% reduction | NT-proBNP ≥30% + ≥300 pg/mL reduction | Proteinuria ≥30% reduction |
| **VGPR** | dFLC <40 mg/L or >90% reduction | — | — |
| **CR** | dFLC <40 mg/L | — | — |
| **sCR** | CR + normal FLC ratio + negative BMPC | — | — |

## How to Run

```r
# Run mrgsolve model
source("al_mrgsolve_model.R")

# Launch Shiny app
shiny::runApp("al_shiny_app.R")

# Render mechanistic map
# dot -Tsvg al_qsp_model.dot -o al_qsp_model.svg
# dot -Tpng -Gdpi=150 -Gsplines=line al_qsp_model.dot -o al_qsp_model.png
```

## Files

| File | Description |
|------|-------------|
| `al_qsp_model.dot` | Graphviz mechanistic map source (14 clusters, 100+ nodes) |
| `al_qsp_model.svg` | Vector mechanistic map (scalable) |
| `al_qsp_model.png` | Raster mechanistic map (150 dpi) |
| `al_mrgsolve_model.R` | mrgsolve ODE model (20 compartments, 7 scenarios) |
| `al_shiny_app.R` | Interactive Shiny dashboard (8 tabs) |
| `al_references.md` | 50 PubMed-linked references (10 sections) |
