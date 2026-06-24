# Sepsis / Systemic Inflammatory Response Syndrome — QSP Model
# 패혈증 정량적 시스템 약리학 모델

> **Directory:** `sepsis/` | **Abbreviation:** SEP | **Date:** 2026-06-24

[![Sepsis QSP Mechanistic Map](sep_qsp_model.png)](sep_qsp_model.svg)

---

## Disease Overview (질환 개요)

**Sepsis** is a life-threatening organ dysfunction caused by a dysregulated host response to infection. Defined by the Sepsis-3 consensus (Singer et al. JAMA 2016), it represents one of the most complex and lethal syndromes in critical care medicine, accounting for >11 million deaths annually worldwide.

**Septic shock** is a subset of sepsis with profound circulatory and cellular/metabolic abnormalities associated with a >40% hospital mortality.

| Parameter | Value |
|-----------|-------|
| Global incidence | ~49 million cases/year (WHO 2020) |
| Global mortality | ~11 million deaths/year |
| ICU mortality (shock) | 40–60% |
| Primary mediators | TNF-α, IL-6, IL-1β, NO (vasoplegia), thrombin (DIC) |
| Diagnosis | SOFA score ≥2 + suspected infection |
| Septic shock criteria | MAP <65 mmHg despite fluids + lactate >2 mmol/L |

---

## Pathophysiology Summary (병태생리)

| Phase | Key Mechanisms |
|-------|----------------|
| **1. Infection & Recognition** | Gram⁻ LPS/Gram⁺ LTA → TLR4/TLR2 → MyD88 → IRAK1/4 → TRAF6 → IKK → NF-κB |
| **2. Cytokine Storm** | NF-κB → TNF-α, IL-1β, IL-6, IL-8 (pro-inflammatory); IL-10, TGF-β (anti-inflammatory) |
| **3. Complement Activation** | C3a/C5a → amplify neutrophil/macrophage activation |
| **4. Coagulation/DIC** | TF expression → thrombin → fibrin; PAI-1↑ → fibrinolysis↓ → thrombocytopenia |
| **5. Vasoplegia** | iNOS↑ → NO↑ → sGC/cGMP → vascular smooth muscle relaxation → MAP↓ |
| **6. Endothelial Dysfunction** | Glycocalyx shedding → capillary leak → tissue edema → oxygen debt |
| **7. HPA Axis** | CRH → ACTH → cortisol (stress response); CIRCI (relative adrenal insufficiency) in 30–40% |
| **8. Organ Failure (MODS)** | Cardiac depression, AKI (tubular injury), ARDS (alveolar flooding), hepatic dysfunction |

---

## Mechanistic Map (기계론적 지도)

**File:** [`sep_qsp_model.dot`](sep_qsp_model.dot) → [`sep_qsp_model.svg`](sep_qsp_model.svg)

**11 Clusters, 124+ nodes, 211 edges:**

| Cluster | Key Components |
|---------|----------------|
| ① Infection & PAMPs | Bacteria, LPS, LTA, Peptidoglycan, PAMP, DAMP, HMGB1, Bacteremia |
| ② Innate Immune Recognition | TLR4/2/9, MD2, MyD88, IRAK1/4, TRAF6, IKK, IκB, NF-κB activation cascade |
| ③ Cytokine Network | TNF-α, IL-1β, IL-6, IL-8, IL-10, IL-12, IL-17, IFN-γ, G-CSF, MCP-1 |
| ④ Coagulation & DIC | TF→FVIIa→FXa→Thrombin→Fibrin→D-dimer; PAI-1/tPA/Plasmin; Protein C, TFPI |
| ⑤ Complement System | Classical/lectin/alternative pathways; C3a, C3b, C5a, MAC |
| ⑥ Vascular Dysfunction | iNOS→NO→vasodilation→SVR↓→MAP↓; glycocalyx shedding→capillary leak→edema |
| ⑦ HPA Axis | CRH→ACTH→Cortisol→GR→anti-inflammatory genes; CIRCI; ADH/vasopressin |
| ⑧ Organ Failure | Cardiac depression, AKI, hepatic dysfunction, ARDS, gut barrier failure, MODS |
| ⑨ Drug PK | PipTazo 2-compartment, vancomycin, norepinephrine, hydrocortisone, vasopressin |
| ⑩ Clinical Endpoints | SOFA, qSOFA, lactate, PCT, CRP, platelets, INR, creatinine, P/F ratio |

**Edge types:**
- 🟢 Activation (green, 162 edges)
- 🔴 Inhibition (red, 33 edges)
- 🔵 Drug effect (blue dashed, 16 edges)

---

## mrgsolve ODE Model (수리 모델)

**File:** [`sep_mrgsolve_model.R`](sep_mrgsolve_model.R)

### Compartments (20 ODEs)

| Module | Compartment | Key Dynamics |
|--------|------------|--------------|
| Infection | `B` | Bacterial load (CFU/mL); logistic growth kb=0.9/h; neutrophil + antibiotic killing |
| Innate immunity | `N` | Neutrophils (cells/µL); TNF-driven mobilisation; infection-site margination |
| | `M` | Activated macrophages; bacteria-driven; IL-10-suppressed |
| Cytokines | `TNF` | TNF-α (pg/mL); macrophage-driven; t½≈2h; cortisol-inhibited |
| | `IL6` | IL-6 (pg/mL); macrophage+TNF-driven; t½≈6h |
| | `IL10` | IL-10 (pg/mL); IL-6-dependent anti-inflammatory feedback |
| | `IL1b` | IL-1β (pg/mL); macrophage-driven; t½≈4h |
| Coagulation | `Th` | Thrombin (nM); cytokine+TF-driven |
| | `F` | Fibrin (µg/mL); DIC indicator |
| | `Plt` | Platelets (×10³/µL); thrombin-mediated consumption |
| Vascular | `NO` | Nitric oxide (µM); iNOS from macrophages/cytokines |
| Organ | `D_tissue` | Tissue damage index (0–1); accumulates from bacteria+cytokines+coagulation |
| | `Lac` | Lactate (mmol/L); hypoperfusion marker; MAP-dependent |
| | `MAP` | Mean arterial pressure (mmHg); NO/damage-driven decrease; vasopressor-restored |
| | `Cr` | Serum creatinine (mg/dL); GFR reduces with tissue damage |
| Antibiotic PK | `AB_C` | Piperacillin central (µg/mL); CL=15 L/h, Vc=10 L |
| | `AB_P` | Piperacillin peripheral; k12=0.5/h, k21=0.28/h |
| Vasopressor | `NE_C` | Norepinephrine plasma; CL=150 L/h, t½≈2 min |
| Steroid PK | `HC_C` | Hydrocortisone central (µg/mL); CL=15 L/h |
| | `Cort` | Total cortisol (µg/dL); endogenous stress response + exogenous HC |

### Drug PK/PD Parameters

| Drug | Route | Dose | t½ | PD Mechanism |
|------|-------|------|----|--------------|
| Piperacillin/Tazobactam | IV | 4.5 g q6h | 1.0 h | %fT>MIC kills bacteria (Emax=0.95, MIC=16 µg/mL) |
| Norepinephrine | IV infusion | 0.1–0.5 µg/kg/min | ~2 min | α1-agonist → SVR↑ → MAP↑ (Emax=30 mmHg, EC50=0.15) |
| Hydrocortisone | IV continuous | 200 mg/day | 1.5 h | GR activation → cytokine suppression (Imax=65%, IC50=5 µg/dL) |
| Vasopressin | IV continuous | 0.03 units/min | 10–20 min | V1a → vasoconstriction → MAP↑ (Emax=15 mmHg) |

### 6 Treatment Scenarios

| Scenario | Treatment | Clinical Trial Basis | Key Outcome |
|----------|-----------|---------------------|-------------|
| S1: Untreated Sepsis | None | Historical cohort | Progressive MAP↓, MOF, >80% mortality model estimate |
| S2: Early Antibiotics (1h) | PipTazo 4.5g q6h from 1h | Kumar et al. Crit Care Med 2006 | Bacterial clearance; 7% mortality reduction per hour |
| S3: Antibiotics + NE | PipTazo + NE 0.2 µg/kg/min | De Backer et al. NEJM 2010 | MAP restored to ≥65; lactate clearance |
| S4: Full Bundle | AB + NE + HC 200mg/day | ADRENAL trial NEJM 2018 | Faster shock reversal; cytokine suppression |
| S5: Delayed Antibiotics | PipTazo starting at 6h | Kumar et al. 2006 (survival cliff) | 6h delay → 14% higher mortality estimate |
| S6: Refractory Shock | AB + NE 0.5 + VP 0.03 units/min + HC | VASST NEJM 2008 | Vasopressin spares NE; steroid reduces pressor dependency |

---

## Shiny App (대시보드)

**File:** [`sep_shiny_app.R`](sep_shiny_app.R)

**6 Interactive Tabs:**

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Age/weight, infection type (E. coli/S. aureus/K. pneumoniae/P. aeruginosa/fungal), comorbidities, initial inoculum, sepsis severity; predicted 28-day mortality |
| 2. PK | Antibiotic concentration-time (48h), AUC/MIC, %fT>MIC for β-lactams, Cmax/Ctrough info boxes |
| 3. Cytokine Dynamics | TNF-α, IL-6, IL-1β (pro-inflammatory) vs IL-10, HMGB1 (anti-inflammatory); tocilizumab toggle |
| 4. Hemodynamics & Organs | MAP with 65 mmHg target, lactate clearance %, SOFA components, creatinine (AKI staging), DIC markers |
| 5. Treatment Scenarios | 6 preset scenarios comparison: bacterial load, MAP, SOFA, 28-day mortality table |
| 6. Biomarkers | 7-panel dashboard: lactate, PCT, CRP, WBC, platelets, INR, bilirubin; SOFA timeline; survival curve |

---

## References (참고문헌)

**File:** [`sep_references.md`](sep_references.md) — **55 PubMed citations in 10 sections**

Key references:
- Rivers et al. NEJM 2001 — Early Goal-Directed Therapy (EGDT)
- Singer et al. JAMA 2016 — Sepsis-3 definitions
- Evans et al. Intensive Care Med 2021 — Surviving Sepsis Campaign 2021
- Russell et al. NEJM 2008 (VASST) — vasopressin
- Venkatesh et al. NEJM 2018 (ADRENAL) — hydrocortisone
- Kumar et al. Crit Care Med 2006 — antibiotic timing and mortality
- Chow et al. Shock 2005 — mathematical model of acute inflammation

---

## Model Deliverables

| Component | File | Specification |
|-----------|------|---------------|
| 🗺️ Mechanistic Map | [`sep_qsp_model.dot`](sep_qsp_model.dot) | **124+ nodes, 11 clusters**, 211 edges |
| 🖼️ SVG | [`sep_qsp_model.svg`](sep_qsp_model.svg) | Vector graphic, scalable |
| 🖼️ PNG | [`sep_qsp_model.png`](sep_qsp_model.png) | 150 dpi raster |
| ⚙️ mrgsolve ODE | [`sep_mrgsolve_model.R`](sep_mrgsolve_model.R) | **20-compartment ODE**, **6 treatment scenarios** |
| 📊 Shiny App | [`sep_shiny_app.R`](sep_shiny_app.R) | **6 tabs** (patient profile·PK·cytokines·hemodynamics·scenarios·biomarkers) |
| 📚 References | [`sep_references.md`](sep_references.md) | **55 PubMed citations** (10 sections) |

---

## How to Run (실행 방법)

```bash
# Render mechanistic map
dot -Tsvg sep_qsp_model.dot -o sep_qsp_model.svg
dot -Tpng -Gdpi=150 sep_qsp_model.dot -o sep_qsp_model.png
```

```r
# Run mrgsolve simulation
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))
source("sep_mrgsolve_model.R")

# Launch Shiny app
install.packages(c("shiny", "shinydashboard", "deSolve", "ggplot2"))
shiny::runApp("sep_shiny_app.R")
```
