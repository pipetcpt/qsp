# Polycythemia Vera (PV) — QSP Model

> **Directory:** `polycythemia-vera/` | **Abbreviation:** PV | **Date:** 2026-06-24  
> **Disease Category:** Chronic Myeloproliferative Neoplasm (MPN)

[![PV QSP 기계론적 지도](pv_qsp_model.png)](pv_qsp_model.svg)

---

## Disease Overview

**Polycythemia vera (PV)** is a clonal myeloproliferative neoplasm (MPN) defined by
pathological erythrocytosis driven by somatic gain-of-function mutations in the **JAK2 tyrosine
kinase gene** — predominantly **JAK2 V617F** (>95% of cases) or JAK2 exon 12 mutations.
The constitutively active JAK2 → STAT5 → STAT3 axis drives EPO-independent erythropoiesis,
resulting in elevated red cell mass, thrombocytosis, and leukocytosis.

| Parameter | Value |
|-----------|-------|
| Prevalence | ~44/100,000 (Western countries) |
| Incidence  | 0.4–2.8/100,000/year |
| Median age | ~60 years |
| Male:Female | 1.2:1 |
| JAK2 V617F | >95% of cases |
| 15-yr OS | ~60–70% |
| MF transformation | 10–25% at 15yr |
| AML transformation | ~5–10% lifetime |

---

## Core Pathophysiology (11 Mechanistic Clusters)

| Cluster | Key Mechanism | Central Components |
|---------|---------------|-------------------|
| **1. Genetic Basis** | JAK2 V617F / exon12 → constitutive kinase activation; co-mutations (TET2, ASXL1, IDH) | JAK2V617F, Clone_HSC, VF_allele |
| **2. JAK-STAT Signaling** | Persistent p-STAT5 → anti-apoptosis (BCL-2, MCL-1), PI3K/AKT/mTOR, CDK4/6 cell cycle | JAK2_kinase, STAT5_act, PI3K, BCL2 |
| **3. Erythropoiesis** | EPO-independent BFU-E → CFU-E → RBC expansion; endogenous EPO suppressed | BFU_E_mut, RBC, HCT, EEC |
| **4. Myeloid / Megakaryocyte** | Thrombopoietin-independent megakaryopoiesis; platelet activation/TXA2; NETs | PLT, Megakaryocyte, NET, TXA2 |
| **5. BM Microenvironment** | Abnormal MK release TGF-β1, PDGF → reticulin fibrosis → collagen deposition | TGFb1, Reticulin, MF_prog, OPN |
| **6. Spleen / EMH** | CD34+ mobilization → extramedullary hematopoiesis → splenomegaly | Spleen_vol, EMH, Portal_HT |
| **7. Vascular / Thrombosis** | ↑ Viscosity (HCT >45%) → endothelial activation → DVT/stroke/MI | Blood_viscosity, Thromb_risk, DVT |
| **8. Inflammatory Cytokines** | IL-6, IL-1β, TNF-α, histamine (basophils), IL-13 → symptoms | MPN_SAF_TSS, Pruritus, Fatigue |
| **9. Disease Progression** | Allele burden ↑ + additional mutations (ASXL1/IDH) → Post-PV MF → AML | PV_to_MF, MF_to_AML |
| **10. Clinical Endpoints** | HCT <45%, SVR35, TSS50, allele burden, CHR, mol. remission | HCT_target, SVR, CHR, QoL_out |
| **11. Drug PK/PD** | 8 drugs: Ruxolitinib · HU · Ropeg-IFN · Aspirin · Phlebotomy · Fedratinib · Anagrelide · Busulfan | RUX_cent, HU_cent, IFN_cent |

---

## mrgsolve ODE Model (22 Compartments)

| Module | Compartments | Key Dynamics |
|--------|-------------|--------------|
| Clonal HSC | `mut_clone`, `wt_clone`, `allele_burden` | Logistic competition; IFN-mediated clone elimination |
| JAK-STAT | `STAT5` | Allele-burden driven activation; JAK inhibitor suppression |
| Erythroid Progenitors | `BFU_E_mut`, `BFU_E_wt` | EPO-independent (mutant) vs EPO-dependent (wt) |
| Mature Blood Cells | `RBC`, `HCT`, `PLT`, `WBC` | Progenitor-driven production; HU/phlebotomy reduction |
| Organ / Symptoms | `Spleen_vol`, `MPN_SAF`, `MF_score`, `Thromb_hazard` | EMH-driven spleen growth; fibrosis TGF-β axis |
| Ruxolitinib PK | `RUX_gut`, `RUX_cent`, `RUX_periph` | 2-compartment oral; ka=2.16/h, CL=17.7L/h, t½~3h |
| HU PK | `HU_cent` | 1-compartment; EC50=8.5mg/L, t½~3.5h |
| Ropeg-IFN PK | `IFN_sc`, `IFN_cent` | SC depot; ka=0.15/day, t½~80–130h |
| Aspirin PK | `ASP_cent` | COX-1 irreversible inhibition |
| Fedratinib PK | `FED_cent` | 1-compartment; CL=4.2L/h, Vd=212L |

---

## Treatment Scenarios (6 Arms, 5-Year Simulation)

| Scenario | Regimen | Mechanism | Key Calibration |
|----------|---------|-----------|-----------------|
| **S0** | Natural history | None | Baseline PV trajectory |
| **S1** | Phlebotomy (8×/yr) + Aspirin 100mg/day | RBC removal + COX-1 inhibition | ECLAP trial (aspirin) |
| **S2** | Hydroxyurea 1500mg/day + Phlebotomy + Aspirin | Ribonucleotide reductase inhibition | PVSG cytoreductive trial |
| **S3** | Ruxolitinib 10mg BID + Phlebotomy + Aspirin | JAK1/2 inhibition (IC50~3nM) | **RESPONSE (NEJM 2015)** |
| **S4** | Ruxolitinib 20mg BID + Aspirin | High-dose JAK1/2 inhibition | MAJIC-PV (Lancet Haematol 2017) |
| **S5** | Ropeginterferon-α2b 100mcg q2w + Phlebotomy + Aspirin | Clone-selective elimination | **PROUD-PV / CONTINUATION-PV** |
| **S6** | Fedratinib 400mg/day + Phlebotomy + Aspirin | JAK2/FLT3 inhibition | JAKARTA-2 study |

---

## Predicted Response at Week 32 (RESPONSE Trial Design)

| Scenario | HCT (%) | PLT (×10⁹/L) | Spleen (mL) | MPN-SAF TSS | HCT Control | SVR35 | TSS50 |
|----------|---------|--------------|-------------|-------------|-------------|-------|-------|
| S0: Natural History | 58–62 | 600–700 | 950–1100 | 25–30 | No | No | No |
| S1: Phlebotomy+ASA | 44–48 | 550–620 | 900–1000 | 22–26 | Partial | No | No |
| S2: HU+Phlebotomy+ASA | 43–46 | 380–450 | 800–900 | 15–20 | Yes (26%) | No | No |
| S3: Ruxolitinib 10mg BID | **<45** | **<400** | **<553** | **<10** | **Yes (60%)** | **Yes (38%)** | **Yes (55%)** |
| S4: Ruxolitinib 20mg BID | **<44** | **<380** | **<520** | **<8** | **Yes (68%)** | **Yes (45%)** | **Yes (62%)** |
| S5: Ropeg-IFN 100mcg q2w | 44–47 | 400–500 | 780–850 | 14–18 | Partial | No | Partial |
| S6: Fedratinib 400mg | **<45** | **<400** | **<540** | **<12** | **Yes (55%)** | **Yes (35%)** | **Yes (48%)** |

*Week 32 values — model-predicted (calibrated to RESPONSE, MAJIC-PV, PROUD-PV clinical trials)*

---

## Key Pharmacodynamic Relationships

### Ruxolitinib Dose-JAK2 Inhibition
- IC50 (JAK2) ≈ 450 ng/mL; Emax = 92%; Hill coefficient = 1.2
- At 10mg BID: Cmax ~1200 ng/mL → ~75% JAK2 inhibition
- At 20mg BID: Cmax ~2400 ng/mL → ~87% JAK2 inhibition
- Sustained inhibition drives SVR35 and TSS50 responses

### Ropeginterferon vs JAK Inhibitors
| Feature | Ruxolitinib | Ropeg-IFN |
|---------|-------------|-----------|
| JAK2 allele burden | Stable (no change) | **↓ Progressive** |
| Molecular remission | Rare | **41% at 5yr** |
| HCT control | Rapid (weeks) | Slower (months) |
| Spleen reduction (SVR35) | 38–45% | 10–15% |
| MPN-SAF TSS reduction | ~55% | ~30% |
| Clone elimination | No | **Yes** |

---

## ELN 2022 Risk Stratification & Treatment Algorithm

```
PV Diagnosis (JAK2V617F + erythrocytosis)
         │
         ├── Low Risk (age <60 + no thrombosis)
         │       └── Phlebotomy (HCT <45%) + Aspirin 100mg/day
         │
         └── High Risk (age ≥60 OR prior thrombosis)
                 ├── HU 500–2000mg/day (first-line)
                 │       ↓ Resistance/Intolerance
                 ├── Ruxolitinib 10–25mg BID (RESPONSE)
                 ├── Ropeg-IFN α-2b 100–250mcg q2w → qmo (PROUD-PV)
                 └── Fedratinib 400mg/day (HU+Rux failed)
```

---

## File Summary

| File | Description | Specification |
|------|-------------|---------------|
| [`pv_qsp_model.dot`](pv_qsp_model.dot) | Graphviz mechanistic map | **146 nodes, 11 clusters** (fdp layout) |
| [`pv_qsp_model.svg`](pv_qsp_model.svg) | Vector SVG rendering | Full interactive pathway map |
| [`pv_qsp_model.png`](pv_qsp_model.png) | Raster PNG rendering | 150 dpi for publications |
| [`pv_mrgsolve_model.R`](pv_mrgsolve_model.R) | mrgsolve ODE model | **22-compartment ODE, 6 treatment scenarios** |
| [`pv_shiny_app.R`](pv_shiny_app.R) | Interactive Shiny dashboard | **8 tabs**: Profile · PK · Hematology · Spleen/Symptoms · Progression · Scenarios · Biomarkers · About |
| [`pv_references.md`](pv_references.md) | Curated literature | **55 PubMed references** (12 sections) |

---

## Quick Start

```r
# Install required packages
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))

# Run the ODE simulation
source("pv_mrgsolve_model.R")

# Launch the Shiny dashboard
library(shiny)
shiny::runApp("pv_shiny_app.R")
```

---

## Clinical Trial Calibration Summary

| Trial | Design | Drug | Primary Endpoint | Outcome |
|-------|--------|------|-----------------|---------|
| **RESPONSE** (NEJM 2015) | Phase 3 RCT | Ruxolitinib vs BAT | SVR35 + HCT control at Wk32 | **38% vs 1%** |
| **RESPONSE-2** (Lancet Oncol 2017) | Phase 3b | Ruxolitinib vs BAT (no splenomegaly) | HCT control Wk28 | **62% vs 19%** |
| **MAJIC-PV** (Lancet Haematol 2017) | Phase 2 RCT | Ruxolitinib vs HU | CHR at 1yr | **41% vs 26%** |
| **CYTOREDUCE** (NEJM 2013) | Phase 3 RCT | HCT <45% vs <50% | Thrombosis + death | **2.7% vs 9.8%** (HR 0.25) |
| **PROUD-PV** (Lancet Haematol 2020) | Phase 3 RCT | Ropeg-IFN vs HU | Non-inferiority (CHR) | Non-inferior + allele ↓ |
| **ECLAP** (NEJM 2004) | Phase 3 RCT | Aspirin 100mg vs placebo | Thromboembolic events | **HR 0.41** (59% risk ↓) |

---

*Model developed for educational QSP research. Clinical decisions should follow current ELN and NCCN guidelines.*
