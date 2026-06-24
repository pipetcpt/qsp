# Alzheimer's Disease (AD) QSP Model

## 질환 개요 (Disease Overview)

### 한국어

알츠하이머병(Alzheimer's Disease, AD)은 가장 흔한 신경퇴행성 질환으로, 전 세계적으로 약 5,500만 명, 한국에서는 약 90만 명이 이 질환을 앓고 있습니다. AD는 서서히 진행하는 인지기능 저하를 특징으로 하며, 기억력 손상에서 시작하여 언어, 판단, 일상생활 능력의 전반적인 저하로 이어집니다.

**아밀로이드 캐스케이드 가설 (Amyloid Cascade Hypothesis)**

AD의 주요 발병기전으로 널리 받아들여지는 아밀로이드 캐스케이드 가설에 따르면, 아밀로이드 전구 단백질(APP)의 비정상적 처리로 인해 생성된 아밀로이드 베타(Aβ42) 펩타이드가 뇌 조직에 축적되어 노인반(senile plaque)을 형성합니다. 이러한 Aβ42 축적은 신경염증 반응을 유발하고, 시냅스 기능 장애 및 신경세포 사멸로 이어집니다.

**타우 병리 (Tau Pathology)**

타우 단백질의 과인산화(hyperphosphorylation)는 신경섬유다발(neurofibrillary tangles, NFT)의 형성으로 이어집니다. 타우 병리는 내후각 피질에서 시작하여 해마, 신피질로 전파되며(Braak 병기), 인지기능 저하와 밀접하게 연관됩니다.

**신경염증 (Neuroinflammation)**

Aβ 축적에 반응하여 미세아교세포(microglia)와 별아교세포(astrocyte)가 활성화됩니다. 초기에는 Aβ 제거를 위한 보호적 반응이지만, 만성 활성화 상태에서는 친염증성 사이토카인(IL-1β, TNF-α, IL-6) 분비가 증가하여 신경세포 손상을 악화시킵니다.

**콜린성 시스템 손상 (Cholinergic System Impairment)**

마이네르트 기저핵(nucleus basalis of Meynert)의 콜린성 신경세포 소실은 대뇌 피질과 해마의 아세틸콜린(ACh) 농도 감소를 초래합니다. 이는 인지기능, 특히 기억력과 주의력 저하에 직접적으로 기여하며, 도네페질 등 아세틸콜린에스테라제(AChE) 억제제 치료의 근거가 됩니다.

### English

Alzheimer's disease (AD) is the most common neurodegenerative disorder and the leading cause of dementia worldwide. It accounts for 60–80% of all dementia cases, affecting approximately 55 million people globally and approximately 900,000 people in South Korea.

**Pathological Hallmarks**

- **Amyloid beta plaques**: Extracellular deposits of aggregated Aβ42 peptides, formed by sequential cleavage of APP by β-secretase (BACE1) and γ-secretase
- **Neurofibrillary tangles (NFTs)**: Intracellular accumulations of hyperphosphorylated tau protein, causing microtubule destabilization and neuronal death
- **Neuroinflammation**: Chronic activation of microglia and astrocytes, contributing to synaptic loss and neurodegeneration
- **Vascular changes**: Cerebral amyloid angiopathy (CAA), impaired blood-brain barrier (BBB), reduced cerebral blood flow

**Genetic Factors**

- **APOE4**: Major genetic risk factor for late-onset AD; increases Aβ accumulation and reduces clearance (1 allele: 3× risk; 2 alleles: 12× risk)
- **PSEN1/PSEN2**: Presenilin mutations cause early-onset familial AD by altering γ-secretase activity, increasing Aβ42/Aβ40 ratio
- **APP mutations**: Rare mutations near secretase cleavage sites cause familial AD or protective effects (A673T)

**Disease Stages**

| Stage | Description | Cognitive Status | Biomarker Changes |
|-------|-------------|-----------------|-------------------|
| Preclinical | No symptoms | Normal | Aβ PET+, CSF Aβ42↓ |
| MCI due to AD | Subtle changes | Objective impairment, independent ADL | Tau PET+, MRI atrophy |
| Mild AD | Memory/cognitive decline | Impaired ADL | All biomarkers positive |
| Moderate AD | Significant impairment | Dependent for complex ADL | Significant neurodegeneration |
| Severe AD | Profound impairment | Fully dependent | Widespread cortical atrophy |

---

## Model Architecture

The AD QSP model integrates five mechanistic components representing the key pathophysiological pathways and pharmacological interventions:

### 1. APP Processing & Amyloid Cascade

This component models the sequential processing of amyloid precursor protein (APP) through both amyloidogenic and non-amyloidogenic pathways:

- APP production from neurons and metabolism via α-secretase (non-amyloidogenic) or BACE1/γ-secretase (amyloidogenic)
- Aβ40 and Aβ42 monomer production, oligomerization, and fibrillization
- Amyloid plaque formation, seeding, and growth kinetics
- Enzymatic clearance via neprilysin (NEP) and insulin-degrading enzyme (IDE)
- BBB-mediated transport (LRP1 efflux, RAGE influx)
- APOE4-dependent modulation of Aβ clearance efficiency

### 2. Tau Pathology Pathway

Models the progression of tau pathology from normal function to NFT formation:

- Normal tau synthesis and microtubule binding
- GSK3β and CDK5-mediated phosphorylation
- Hyperphosphorylated tau accumulation and oligomer formation
- NFT assembly kinetics (Braak-stage-based spread approximation)
- Tau clearance via autophagy-lysosomal pathway and proteasome
- Cross-talk with Aβ pathology (Aβ-induced GSK3β activation)

### 3. Neuroinflammation (Microglia/Astrocyte)

Captures the dual role of neuroinflammation in AD:

- Resting microglia activation by Aβ oligomers and DAMPs
- M1 (pro-inflammatory) vs. M2 (anti-inflammatory/phagocytic) microglial states
- Pro-inflammatory cytokine production (IL-1β, TNF-α, IL-6, IL-18)
- Astrocyte reactivity and GFAP expression
- Complement system activation (C1q, C3, CR3)
- TREM2-mediated microglial phagocytosis of Aβ
- NF-κB signaling pathway

### 4. Cholinergic/Synaptic Neurotransmission

Models the cholinergic hypothesis and synaptic function:

- ACh synthesis (ChAT activity), release, and hydrolysis (AChEI target)
- Muscarinic (M1, M2) and nicotinic (α7) receptor signaling
- Glutamatergic transmission: AMPA and NMDA receptor-mediated LTP/LTD
- NMDA receptor overstimulation and excitotoxicity (memantine target)
- Synaptic density dynamics (formation vs. Aβ/tau-mediated loss)
- BDNF-TrkB neuroprotective signaling
- Cognitive function scoring linked to synaptic integrity

### 5. Drug PK/PD

Pharmacokinetic and pharmacodynamic models for four approved/emerging therapies:

**Donepezil (AChE Inhibitor)**
- 2-compartment oral PK model (ka, Vd, CL, brain Kp)
- Reversible AChE inhibition (IC50-based Emax model)
- ACh concentration-dependent cognitive improvement

**Memantine (NMDA Receptor Antagonist)**
- 3-compartment oral PK model
- Voltage-dependent, uncompetitive NMDA receptor block
- Neuroprotection via glutamate excitotoxicity attenuation

**Lecanemab (Anti-Aβ Protofibril mAb)**
- 2-compartment IV PK model with target-mediated drug disposition (TMDD)
- BBB penetration via IgG transcytosis (Kp ~ 0.001)
- Dose-dependent Aβ protofibrils and plaque removal
- ARIA (amyloid-related imaging abnormalities) risk modeling

**Donanemab (Anti-N3pE Aβ mAb)**
- Similar TMDD framework targeting N-terminally truncated Aβ (N3pE)
- Potent plaque clearance kinetics (amyloid clearance in ~76% at 76 weeks)
- Dose-suspension upon amyloid clearance (TRAILBLAZER protocol)

---

## Files in This Directory

| File | Description |
|------|-------------|
| `ad_qsp_model.dot` | Graphviz DOT source for mechanistic map (100+ nodes, 10+ clusters) |
| `ad_qsp_model.svg` | Vector format mechanistic map |
| `ad_qsp_model.png` | PNG format mechanistic map (150 dpi) |
| `ad_mrgsolve_model.R` | mrgsolve ODE model with 20+ compartments, 7 treatment scenarios |
| `ad_shiny_app.R` | Interactive Shiny dashboard (6 tabs) |
| `ad_references.md` | 40+ PubMed references organized by topic |
| `README.md` | This file |

---

## Key Model Parameters

| Parameter | Symbol | Value | Unit | Source |
|-----------|--------|-------|------|--------|
| Donepezil ka (absorption) | ka_don | 1.2 | h⁻¹ | Tiseo et al. 1998 |
| Donepezil Vd | Vd_don | 0.7 | L/kg | Tiseo et al. 1998 |
| Donepezil CL | CL_don | 0.13 | L/h/kg | Tiseo et al. 1998 |
| Donepezil brain Kp | Kp_don | 0.6 | — | Estimated from tissue distribution |
| Donepezil AChE IC50 | IC50_ache | 0.006 | μg/mL | Sugimoto et al. 1995 |
| Memantine Vd | Vd_mem | 9 | L/kg | Wesemann et al. 1983 |
| Memantine CL | CL_mem | 0.13 | L/h/kg | Wesemann et al. 1983 |
| Memantine NMDAR IC50 | IC50_nmda | 1.0 | μg/mL | Parsons et al. 2007 |
| Lecanemab Vd (central) | Vd_lec | 0.05 | L/kg | Swanson et al. 2021 |
| Lecanemab t½ (terminal) | t_half_lec | 270 | h (~11 days) | Swanson et al. 2021 |
| Lecanemab CNS Kp | Kp_lec | 0.001 | — | Estimated (IgG BBB penetration) |
| Abeta42 production rate | k_prod | 0.15 | nM/h | Bateman et al. 2006 |
| Abeta42 clearance (NEP/IDE) | k_clear | 0.3 | h⁻¹ | Leissring et al. 2003 |
| Tau phosphorylation rate | k_phos | 0.05 | h⁻¹ | Estimated |
| NFT formation rate | k_agg | 0.01 | h⁻¹ | Estimated |
| Microglial activation rate | k_micro | 0.02 | h⁻¹ | Estimated |
| Synaptic damage rate | k_syn_dam | 0.03 | h⁻¹ | Estimated |
| APOE4 effect on clearance | APOE4_mult | 1.3 | fold | Liu et al. 2013 |

---

## Treatment Scenarios

| Scenario | Drug(s) | Dose | Regimen | Clinical Trial Basis |
|----------|---------|------|---------|---------------------|
| No treatment | — | — | — | Natural history |
| Donepezil | AChE inhibitor | 10 mg | QD oral | Rogers et al. 1996 |
| Memantine | NMDAR antagonist | 20 mg | QD oral | Reisberg et al. 2003 |
| Combination | Don + Mem | 10 + 20 mg | QD oral | Tariot et al. 2004 |
| Lecanemab | Anti-Aβ mAb | 10 mg/kg | Q2W IV | CLARITY-AD (van Dyck 2023) |
| Donanemab | Anti-Aβ mAb (N3pE) | 1500 mg | Q4W IV | TRAILBLAZER-ALZ-2 (Sims 2023) |
| Triple therapy | Don + Mem + Lec | Standard | Combined | Theoretical combination |

---

## Clinical Trial Context

### CLARITY-AD (Lecanemab)

- **Citation**: Van Dyck CH et al. *N Engl J Med* 2023;388:9-21 (PMID: 36449413)
- **Design**: Phase 3, randomized, double-blind, placebo-controlled
- **Population**: 1,795 patients with early AD (MCI or mild dementia with confirmed Aβ pathology)
- **Intervention**: Lecanemab 10 mg/kg IV every 2 weeks for 18 months
- **Primary endpoint**: Change from baseline in CDR-SB at 18 months
- **Result**: **27% slowing of decline** vs placebo (CDR-SB change: 1.21 vs 1.66; difference -0.45, p < 0.001)
- **Amyloid reduction**: -55.48 centiloids vs +3.64 for placebo (p < 0.001)
- **Secondary endpoints**: ADAS-Cog14, ADCOMS, ADCS-MCI-ADL all significantly improved
- **ARIA**: ARIA-E in 12.6%, ARIA-H in 17.3% of lecanemab group

### TRAILBLAZER-ALZ 2 (Donanemab)

- **Citation**: Sims JR et al. *JAMA* 2023;330(6):512-527 (PMID: 37459141)
- **Design**: Phase 3, randomized, double-blind, placebo-controlled
- **Population**: 1,736 patients with early symptomatic AD (low/medium tau stratum primary analysis)
- **Intervention**: Donanemab 1,500 mg IV every 4 weeks (transitioned to placebo after amyloid clearance)
- **Primary endpoint**: Change from baseline in iADRS at 76 weeks
- **Result**: **35% slowing of clinical progression** in combined population; 40% in low/medium tau subgroup
- **Amyloid clearance**: Complete clearance (< 24.1 centiloids) in 76% of patients by 76 weeks
- **iADRS result**: -6.02 vs -9.27 for placebo (difference 3.25, p < 0.001)
- **ARIA**: ARIA-E in 24.0%, ARIA-H in 31.4% of donanemab group

### Supporting Symptomatic Trial Data

| Trial | Drug | N | Duration | Primary Result |
|-------|------|---|----------|----------------|
| Rogers et al. 1996 | Donepezil 5/10 mg | 468 | 24 weeks | ADAS-Cog: -2.49 (10 mg) vs +0.52 placebo |
| Reisberg et al. 2003 | Memantine 20 mg | 252 | 28 weeks | CIBIC-plus significant improvement |
| Tariot et al. 2004 | Donepezil + Memantine | 404 | 24 weeks | SIB: +0.9 vs -2.5 placebo (p = 0.002) |

---

## Running the Models

### Mechanistic Map (DOT to PNG)

```bash
# Generate SVG
dot -Tsvg ad_qsp_model.dot -o ad_qsp_model.svg

# Generate PNG at 150 dpi
dot -Tpng -Gdpi=150 ad_qsp_model.dot -o ad_qsp_model.png
```

### mrgsolve Model

```r
library(mrgsolve)
source("ad_mrgsolve_model.R")
# Runs automatically when sourced, generates simulation plots
# for all 7 treatment scenarios over 104 weeks
```

### Shiny App

```r
library(shiny)
runApp("ad_shiny_app.R")
# Opens interactive dashboard at localhost with 6 tabs:
# 1. Patient Profile
# 2. Drug PK
# 3. PD Biomarkers
# 4. Clinical Endpoints
# 5. Scenario Comparison
# 6. Biomarker Panel
```

---

## Shiny Dashboard Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Age, APOE4 status, baseline cognitive scores, disease stage selection |
| 2. Drug PK | Plasma and CSF concentration-time profiles for all drugs |
| 3. PD Biomarkers | Aβ42, tau, phospho-tau, neuroinflammation markers over time |
| 4. Clinical Endpoints | CDR-SB, ADAS-Cog, MMSE trajectories with clinical trial benchmarks |
| 5. Scenario Comparison | Side-by-side comparison of all 7 treatment scenarios |
| 6. Biomarker Panel | PET amyloid centiloids, CSF biomarker ratios, synaptic density |

---

## Model Limitations & Assumptions

1. **Simplified disease pathology**: Three main pathways modeled rather than full mechanistic detail; spatial spreading (Braak staging) not explicitly implemented
2. **Donepezil bioavailability**: Assumed 100% (actual approximately 100% based on clinical studies; well-supported assumption)
3. **Lecanemab CNS penetration**: Modeled as passive diffusion with fixed Kp = 0.001; actual mechanism involves FcRn-mediated transcytosis and is concentration-dependent
4. **Tau propagation**: Prion-like tau spreading not modeled spatially; Braak stage progression approximated as a monotonic progression variable
5. **Individual variability**: Modeled as fixed effects only; no population PK/PD variability or Monte Carlo simulation implemented in base model
6. **Neurodegeneration irreversibility**: Neuronal loss modeled as a slow, irreversible process; no neuroregeneration or plasticity mechanisms included
7. **ARIA mechanistic model**: ARIA risk modeled empirically based on Aβ efflux rate; the true immunological cascade mediating ARIA is not fully captured
8. **Drug interactions**: Donepezil-memantine combination modeled as additive; potential pharmacokinetic interactions (CYP2D6, CYP3A4) not explicitly simulated
9. **Comorbidities**: Vascular contributions to cognitive impairment (VCID), cerebrovascular disease not included
10. **Biomarker calibration**: PET amyloid centiloid scale conversion is approximate; CSF biomarker thresholds based on published reference ranges

---

## References

Key references — see `ad_references.md` for the complete annotated list of 40+ citations:

### Disease Pathology & Mechanisms
- Jack CR Jr et al. NIA-AA Research Framework 2018 (PMID: 29653606)
- Hardy J, Selkoe DJ. The amyloid hypothesis of AD. *Science* 2002 (PMID: 12130773)
- Braak H, Braak E. Neuropathological stageing of Alzheimer-related changes. *Acta Neuropathol* 1991 (PMID: 1759558)

### Clinical Trials — Symptomatic Therapies
- Rogers SL et al. Donepezil Phase III trial. *Neurology* 1996 (PMID: 8649628)
- Reisberg B et al. Memantine in moderate-to-severe AD. *NEJM* 2003 (PMID: 12672860)
- Tariot PN et al. Memantine + donepezil combination. *JAMA* 2004 (PMID: 15010441)

### Clinical Trials — Disease-Modifying Therapies
- van Dyck CH et al. CLARITY-AD Trial 2023 (PMID: 36449413)
- Sims JR et al. TRAILBLAZER-ALZ 2 2023 (PMID: 37459141)

### PK/PD & QSP Modeling
- Tiseo PJ et al. Donepezil PK. *Br J Clin Pharmacol* 1998 (PMID: 9723830)
- Swanson CJ et al. Lecanemab PK/PD. *Alzheimers Res Ther* 2021 (PMID: 34118953)
- Geerts H et al. QSP model for AD. *CPT Pharmacometrics Syst Pharmacol* 2013 (PMID: 23475765)
- Bateman RJ et al. Human amyloid beta kinetics. *Nat Med* 2006 (PMID: 16906156)

---

*Model developed as part of the QSP Disease Model Library. Date: 2026-06-20.*
