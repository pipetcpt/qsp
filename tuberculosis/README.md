# Tuberculosis (결핵) QSP Model

> **Quantitative Systems Pharmacology model for pulmonary tuberculosis caused by *Mycobacterium tuberculosis* (Mtb).**
> Covers bacterial subpopulation dynamics (AR/SR/NR), innate & adaptive immunity, granuloma formation,
> and 4-drug RIPE pharmacokinetics/pharmacodynamics including MDR-TB (Bedaquiline-based) scenarios.

---

## 📋 Disease Overview

| Feature | Details |
|---------|---------|
| **Disease** | Pulmonary Tuberculosis |
| **Pathogen** | *Mycobacterium tuberculosis* (aerobic, slow-growing, acid-fast bacillus) |
| **ICD-10** | A15 (respiratory TB) |
| **Global Burden** | ~10.6 million new cases/year; ~1.6 million deaths (WHO 2023) |
| **Latent TB** | ~1/3 of global population infected (LTBI); 5–10% lifetime reactivation risk |
| **MDR-TB** | ~450,000 cases/year; RIF + INH dual resistance |
| **Key Immune Players** | Alveolar macrophages, Th1 CD4⁺ T cells, IFN-γ, TNF-α, granuloma |
| **Standard Treatment** | 2HRZE / 4HR (6 months total); MDR: BDQ-Pa-LZD ± PZA |

---

## 🔬 Mechanistic Map

[![TB QSP Mechanistic Map](tb_qsp_model.png)](tb_qsp_model.svg)

> *Click the image to open the full interactive SVG. The map contains 115+ nodes across 9 mechanistic clusters.*

### Clusters

| # | Cluster | Key Components |
|---|---------|----------------|
| 1 | **Mtb Bacterial States** | AR_free, AR_intra, SR_intra, NR_persist, DR variants |
| 2 | **Innate Immunity** | AM resting/infected/activated, MDM, DC, NK, neutrophil |
| 3 | **Adaptive Immunity** | Th1, Th17, Treg, CTL (CD8+), Memory T cells, MAIT, γδT |
| 4 | **Cytokine Network** | IFN-γ, TNF-α, IL-12, IL-10, IL-6, IL-17A, TGF-β, CXCL10 |
| 5 | **Granuloma Dynamics** | Early → organized → caseous; cavity formation; calcification |
| 6 | **Drug PK** | RIF, INH (NAT2 polymorphism), PZA, EMB, BDQ plasma/lung |
| 7 | **Drug PD** | Emax models per drug/bacterial state; MIC; resistance (rpoB, katG) |
| 8 | **Clinical Endpoints** | Culture conversion, smear, FEV₁, treatment success/failure, relapse |
| 9 | **Host Risk Factors** | HIV, diabetes, malnutrition, BCG, LTBI, anti-TNF therapy |

---

## ⚙️ mrgsolve ODE Model

**File:** `tb_mrgsolve_model.R`

### State Variables (17 ODEs)

| Group | Variables | Description |
|-------|-----------|-------------|
| **Drug PK** | RIF_gut, RIF_c, INH_c, PZA_c, EMB_c, BDQ_c | 1-cmpt PK; RIF auto-induction; INH NAT2 |
| **Bacteria** | AR, SR, NR | Logistic growth; drug Emax killing |
| **Immune** | UM, IM, AM | Macrophage states (resting → infected → activated) |
| **Adaptive** | Th1, IFNg, TNFa, IL10 | Cytokine-driven Th1 feedback loop |

### Drug PK Parameters (Calibrated)

| Drug | Cmax | t½ | Key PK feature |
|------|------|-----|----------------|
| Rifampicin | 8–24 mg/L | 3–5 h | CYP3A4 auto-induction; lung penetration 25% |
| Isoniazid | 3–5 mg/L | 1–4 h | NAT2 acetylation polymorphism (slow/fast) |
| Pyrazinamide | 20–60 mg/L | 9–10 h | Active at acid pH (≤5.5); renal clearance |
| Ethambutol | 2–5 mg/L | 7 h | Renal clearance; bacteriostatic |
| Bedaquiline | 0.4–1.5 mg/L | 5.5 months | High Vd; M2 metabolite; ATP synthase inhibition |

### 6 Treatment Scenarios

| # | Scenario | Key Parameters |
|---|----------|----------------|
| 1 | **Natural History** | No drugs; Mtb grows to plateau (host immunity limits) |
| 2 | **Standard RIPE** | 2HRZE / 4HR; full adherence; culture conversion ~8–12 weeks |
| 3 | **Poor Adherence** | 30% missed doses; incomplete sterilization; relapse risk |
| 4 | **MDR-TB → BDQ** | RIF+INH resistant (rpoB, katG); BDQ + PZA regimen |
| 5 | **HIV Co-infected** | phi_Th1 = 0.2 (CD4<200); impaired bacterial clearance |
| 6 | **Diabetic Host** | phi_Mact = 0.4 (HbA1c 9%); macrophage dysfunction |

---

## 📊 Shiny App

**File:** `tb_shiny_app.R`

### 7 Tabs

| Tab | Content |
|-----|---------|
| ① Patient Profile | Disease overview, pathophysiology, treatment regimens, mechanistic map thumbnail |
| ② Drug PK | RIPE plasma concentration profiles; NAT2 polymorphism; Cmax vs. targets |
| ③ Bacterial Dynamics | AR/SR/NR subpopulation trajectories; total burden; drug kill rates |
| ④ Immune Response | Macrophage (UM/IM/AM), Th1, cytokine (IFN-γ, TNF-α, IL-10) dynamics |
| ⑤ Clinical Endpoints | Culture conversion, smear positivity, treatment outcome summary table |
| ⑥ Scenario Comparison | All 6 scenarios overlaid; bacteria + immune panels |
| ⑦ PD / Biomarkers | Emax concentration-effect curves; target attainment; biomarker table at selected day |

### Running the App

```r
# Install dependencies
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr",
                   "ggplot2","tidyr","DT"))

# Launch
shiny::runApp("tuberculosis/tb_shiny_app.R")
```

---

## 📚 References

**File:** `tb_references.md` — 46 curated PubMed references covering:
- Epidemiology (WHO 2023 Global TB Report)
- Bacterial biology & pathophysiology (Mtb survival, dormancy)
- Immunology (granuloma, Th1, cytokines)
- QSP / mathematical modeling (Wigginton & Kirschner 2001; Pienaar et al.)
- Drug PK (Wilkins, Peloquin, Smythe — RIF autoinduction)
- Drug PD (Gumbo, Jayaram, Srivastava — Emax models)
- Bedaquiline / MDR-TB (Diacon, Svensson, Conradie)
- HIV–TB co-infection; diabetes–TB interaction
- Clinical trials (STREAM, REMoxTB, OFLOTUB)

---

## 🗂️ Files

```
tuberculosis/
├── tb_qsp_model.dot          # Graphviz mechanistic map source (115+ nodes, 9 clusters)
├── tb_qsp_model.svg          # Vector graphic (interactive)
├── tb_qsp_model.png          # Rasterized map (150 dpi)
├── tb_mrgsolve_model.R       # mrgsolve ODE model + 6 scenarios + plots
├── tb_shiny_app.R            # Shiny dashboard (7 tabs)
├── tb_references.md          # 46 PubMed references
└── README.md                 # This file
```

---

## 🔑 Key Scientific Insights

1. **Bacterial heterogeneity** drives treatment duration: AR killed by INH/RIF in weeks, but SR and NR require PZA/BDQ and months of treatment for sterilization.
2. **Granuloma** is both a host defense (containing Mtb) and a Mtb survival niche (hypoxia → NR state); TNF-α blockade (biologics) destabilizes it.
3. **NAT2 polymorphism** causes 4-fold variability in INH exposure: slow acetylators achieve higher Cmax (better efficacy) but higher hepatotoxicity risk.
4. **RIF auto-induction** reduces its own plasma concentration by ~30–60% after 2 weeks; this must be accounted for in PK models.
5. **HIV co-infection** (CD4<200) impairs Th1 generation, reducing IFN-γ-driven macrophage activation → slower bacterial clearance, higher relapse risk.
6. **MDR-TB** management with BDQ-Pa-LZD (ZeNix/TB-PRACTECAL trials) achieves ~89% success, compared to ~57% for historical injectable-based regimens.

---

*Generated by Claude Code Routine (CCR) — 2026-06-27*
