# Hereditary Angioedema (HAE) — QSP Model
# 유전성 혈관부종 정량적 시스템 약리학 모델

> **Directory:** `hereditary-angioedema/` | **Abbreviation:** HAE | **Date:** 2026-06-24

[![HAE QSP Mechanistic Map](hae_qsp_model.png)](hae_qsp_model.svg)

---

## Disease Overview (질환 개요)

**Hereditary Angioedema (HAE)** is a rare, life-threatening autosomal dominant disorder caused by deficiency (Type I) or dysfunction (Type II) of C1-inhibitor (C1-INH), or by gain-of-function mutations in coagulation factor XII (Type III). The result is uncontrolled activation of the **kallikrein-kinin system (KKS)**, leading to excess **bradykinin (BK)** generation, which drives vascular permeability increases and episodic angioedema.

| Parameter | Value |
|-----------|-------|
| Prevalence | 1:50,000 people |
| Inheritance | Autosomal dominant |
| Gene | *SERPING1* (C1-INH), *F12* (FXII) |
| Primary mediator | Bradykinin (B2R-dependent) |
| Diagnostic | C4 ↓, C1-INH antigen ↓ (Type I), C1-INH function ↓ |
| Untreated mortality | ~30% (laryngeal asphyxiation, historical) |

---

## HAE Types (분류)

| Type | C1-INH Antigen | C1-INH Function | C4 | FXII mutation | Predominance |
|------|----------------|-----------------|-----|---------------|--------------|
| **Type I** | Low (<50%) | Low | Low | No | ~85% |
| **Type II** | Normal or ↑ | Low | Low | No | ~15% |
| **Type III** | Normal | Normal | Normal | Yes (Thr328Lys) | Rare, female |

---

## Mechanistic Map (기계론적 지도)

**File:** [`hae_qsp_model.dot`](hae_qsp_model.dot) → [`hae_qsp_model.svg`](hae_qsp_model.svg)

**12 Clusters, 120+ nodes:**

| Cluster | Key Components |
|---------|----------------|
| ① Genetic Basis | *SERPING1* gene, C1-INH protein, HAE type mutations, *FXII* Thr328Lys |
| ② Contact Activation | FXII → FXIIa (surface autoactivation), prekallikrein, HMWK, amplification loop |
| ③ Kallikrein-Kinin System | Plasma kallikrein, HMWK cleavage, bradykinin (BK), degradation kininases |
| ④ C1-INH Biology | SERPIN mechanism, C1-INH:FXIIa complex, C1-INH:kallikrein, complement C4 depletion |
| ⑤ BK Receptor Signaling | B2R (constitutive), B1R (inducible), Gq/IP3/Ca²⁺/eNOS/NO pathway |
| ⑥ Vascular Effects | EC tight junctions, VE-cadherin, plasma extravasation → subcutaneous/laryngeal edema |
| ⑦ Attack Pathogenesis | Triggers (emotional, surgical, ACEi, estrogen), attack timeline, prodrome, mortality |
| ⑧ Acute Treatment | Icatibant (B2R Ki=0.47nM), C1-INH IV, Ecallantide, recombinant C1-INH |
| ⑨ Prophylaxis | Berotralstat (IC50=3.7nM), Lanadelumab (KD<100pM), C1-INH SC, danazol |
| ⑩ Inflammatory Amplification | IL-1β/B1R upregulation, coagulation crosstalk, vWF/P-selectin |
| ⑪ Clinical Endpoints | Attack frequency, ACE test, AE-QoL, laryngeal risk, biomarkers |
| ⑫ Diagnostic Workup | C4 → C1-INH antigen → C1-INH function → C1q → genetic testing |

---

## mrgsolve ODE Model (수리 모델)

**File:** [`hae_mrgsolve_model.R`](hae_mrgsolve_model.R)

### Compartments (20 total)

| Module | Compartments | Key Dynamics |
|--------|-------------|--------------|
| Icatibant PK | `A_ICA_depot`, `A_ICA_C`, `A_ICA_P` | SC depot → 2-compartment; ka=0.74/h, CL=15.5 L/h, Vc=29 L |
| C1-INH IV PK | `A_C1INH_IV` | 1-compartment; CL=0.051 L/h, Vd=3.3 L, t½=45h |
| C1-INH SC PK | `A_C1INH_SC`, `A_C1INH_SC_C` | SC depot; F=43%, ka=0.025/h |
| Berotralstat PK | `A_BER_gut`, `A_BER_C` | Oral absorption F=57%; Vd=268 L, t½=93h |
| Lanadelumab PK | `A_LAN_depot`, `A_LAN_C`, `A_LAN_P` | SC → 2-compartment; F=61%, t½=17d |
| C1-INH Biology | `C1INH_free` | Liver synthesis – kininase consumption – degradation |
| Contact Activation | `FXII_act` | FXIIa kinetics; C1-INH inhibition; trigger forcing |
| KKS | `Kallikrein_act` | FXIIa-driven; inhibited by C1-INH, berotralstat, lanadelumab |
| Bradykinin | `BK_plasma` | Kallikrein-driven synthesis; ACE degradation |
| B2R Dynamics | `B2R_free`, `B2R_bound` | BK binding (kon/koff); icatibant competitive block; internalization |
| Vascular Permeability | `VP` | B2R_bound → NO/PGI2-driven; Emax model with Hill=1.8 |
| Swelling | `SW_score` | VP threshold-driven edema formation; resolution kinetics |

### Drug PK Parameters

| Drug | Route | Dose | t½ | Key PD |
|------|-------|------|----|--------|
| Icatibant | SC | 30 mg | 1.3 h | Ki(B2R) = 0.47 nM |
| C1-INH IV (Berinert) | IV | 20 IU/kg | 45 h | C1-INH replacement → restored inhibition |
| C1-INH SC (Haegarda) | SC | 60 IU/kg 2×/wk | 45 h | Steady-state C1-INH ↑ >40% |
| Berotralstat | PO | 150 mg QD | 93 h | IC50(Kal) = 3.7 nM, Emax = 92% |
| Lanadelumab | SC | 300 mg Q2W | 17 d | KD(prekallikrein) < 100 pM, Emax = 93% |

### Treatment Scenarios

| Scenario | Drug | Key Outcome | Clinical Trial |
|----------|------|-------------|----------------|
| S1: Untreated HAE | — | BK ↑10×, SW max, laryngeal risk | — |
| S2: Icatibant SC | Icatibant 30 mg | B2R blocked ~98%; resolution 2–4h | FAST-1/3 (NEJM 2010/2011) |
| S3: C1-INH IV | Berinert 20 IU/kg | C1-INH restored; BK ↓↓ | Cicardi 2012 NEJM |
| S4: Berotralstat QD | 150 mg oral | 44% attack reduction; Kal IC50 steady-state | BELO 2020 |
| S5: Lanadelumab Q2W | 300 mg SC | 87% attack reduction; 32% attack-free | HELP OLE 2020 |
| S6: C1-INH SC | Haegarda 60 IU/kg | 95% attack reduction; C4 normalized | CONFIDENT 2017 |

---

## Shiny Dashboard (인터랙티브 대시보드)

**File:** [`hae_shiny_app.R`](hae_shiny_app.R)

**6 Tabs:**

| Tab | Content |
|-----|---------|
| ① Patient Profile | HAE type, C1-INH level, diagnostic criteria, treatment overview, attack timeline |
| ② PK Profiles | Icatibant / Berotralstat / Lanadelumab / C1-INH plasma concentration-time curves |
| ③ KKS Biology | FXIIa activation, kallikrein activity, bradykinin dynamics, B2R occupancy |
| ④ Clinical Endpoints | VP index, swelling score, C1-INH %, max BK fold, attack duration |
| ⑤ Scenario Comparison | Side-by-side comparison of all 6 treatment scenarios (BK, swelling) |
| ⑥ Biomarkers | C4 proxy, C1-INH %, kallikrein suppression, biomarker reference table |

**Run:**
```r
library(shiny)
shiny::runApp("hereditary-angioedema/hae_shiny_app.R")
```

---

## References (참고문헌)

**File:** [`hae_references.md`](hae_references.md) — **58 citations** across 12 sections

Key references:
- Zuraw BL (2008) N Engl J Med 359:1027 — Disease overview
- Maurer M (2018) N Engl J Med 378:1141 — Lanadelumab HELP trial
- Cicardi M (2010) N Engl J Med 363:532 — Icatibant FAST-1 trial
- Farkas H (2020) Allergy 75:1683 — Berotralstat BELO trial
- Craig T (2017) JACI Pract 5:1538 — C1-INH SC CONFIDENT trial

---

## How to Run (실행 방법)

```bash
# 1. Render mechanistic map
fdp -Tsvg hae_qsp_model.dot -o hae_qsp_model.svg
fdp -Tpng -Gdpi=150 hae_qsp_model.dot -o hae_qsp_model.png

# 2. Run mrgsolve ODE model
Rscript hereditary-angioedema/hae_mrgsolve_model.R

# 3. Launch Shiny dashboard
R -e "shiny::runApp('hereditary-angioedema/hae_shiny_app.R')"
```

---

*Model built by Claude Code Routine | Date: 2026-06-24*
*Based on: WAO/EAACI HAE Guidelines 2020; HELP/BELO/CONFIDENT clinical trials*
