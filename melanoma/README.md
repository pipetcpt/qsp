# Melanoma — QSP Model

> **악성 흑색종 (Melanoma)**  
> BRAF V600E/K 돌연변이 → 구성적 MEK-ERK 활성 → 세포 증식 · 생존
> → 표적치료(BRAFi/MEKi) + 면역관문억제제(PD-1/CTLA-4) 병용

---

## Disease Overview

Cutaneous melanoma is the most lethal skin cancer, arising from neoplastic transformation of melanocytes. Approximately **50% of cutaneous melanomas** harbor BRAF V600E/K driver mutations, making this one of the best-characterized oncogene-driven tumors and an exemplary target for QSP modeling combining targeted therapy and immunotherapy.

| Parameter | Heterozygous BRAF V600E | NRAS-mutant | WT / NF1-loss |
|-----------|------------------------|-------------|----------------|
| Frequency | ~45-50% | ~20% | ~15-25% |
| Key driver | BRAF→MEK→ERK | NRAS→CRAF/MEK→ERK | NF1 loss→RAS↑→MEK→ERK |
| Targeted Rx | BRAFi + MEKi | MEKi (modest) | None |
| ICI response | ORR ~58% (nivo+ipi) | ORR ~50% | ORR ~50-60% (high TMB) |
| 5-yr OS | ~28-52% (depending on Rx) | ~20-35% | ~25-40% |

---

## Mechanistic Map

[![Melanoma QSP Mechanistic Map](melanoma_qsp_model.png)](melanoma_qsp_model.svg)

> Click the image to view the full interactive SVG map.

**Map Statistics:**
- **12 clusters** (pathway modules)
- **159+ nodes** (biological entities, drugs, enzymes, checkpoints)
- **Pathways covered:**
  - BRAF inhibitor PK (vemurafenib/dabrafenib: CYP3A4, OATP1B1)
  - MEK inhibitor PK (trametinib/cobimetinib: CYP3A4, UGT1A1)
  - Immune checkpoint inhibitor PK (pembrolizumab/nivolumab/ipilimumab: TMDD, FcRn)
  - MAPK signaling (RAS→BRAF V600E/K→MEK1/2→ERK1/2→nucleus)
  - Transcription factors & cell cycle (MITF, MYC, CCND1, CDK4/6, RB1)
  - PI3K/AKT/mTOR pathway (PTEN loss, resistance bypass)
  - Intrinsic apoptosis (BCL-2 family, BIM, caspase cascade)
  - Antigen presentation & T cell priming (MHC-I/II, CD28, B7, DC)
  - T cell biology & checkpoint regulation (PD-1, CTLA-4, LAG-3, TIM-3, TIGIT)
  - Tumor biology & microenvironment (VEGF, TGF-β, MDSC, TAM, Treg)
  - Acquired resistance (BRAF amplification, NRAS secondary mutation, MEK1/2 mutation)
  - Clinical endpoints (ORR, PFS, OS, LDH, S100B, ctDNA, TMB)

---

## Files

| File | Description |
|------|-------------|
| `melanoma_qsp_model.dot` | Graphviz source (159+ nodes, 12 clusters) |
| `melanoma_qsp_model.svg` | Vector mechanistic map (interactive) |
| `melanoma_qsp_model.png` | Raster map (150 dpi) |
| `melanoma_mrgsolve_model.R` | ODE-based PK/PD model (16 compartments, 6 scenarios) |
| `melanoma_shiny_app.R` | Shiny dashboard (6 tabs, interactive) |
| `melanoma_references.md` | 63 curated PubMed references |

---

## mrgsolve ODE Model

### Compartments (16 ODEs)

**PK Compartments:**
```
GUT_BRAF     BRAF inhibitor GI absorption (oral)
CENT_BRAF    BRAF inhibitor central plasma (µg/mL)
GUT_MEK      MEK inhibitor GI absorption (oral)
CENT_MEK     MEK inhibitor central plasma (ng/mL)
CENT_ICI     Immune checkpoint inhibitor central (µg/mL)
PERI_ICI     ICI peripheral compartment
```

**PD Compartments:**
```
ERK_act      ERK relative activity (1=uninhibited, 0=fully blocked)
RESIST       Acquired resistance factor (0→1 over months under BRAFi)
TUMOR        Tumor burden (normalized: 1=baseline SLD)
CD8_TIL      CD8+ tumor-infiltrating lymphocytes (relative)
PD1_RO       PD-1 receptor occupancy (0–1)
CTLA4_RO     CTLA-4 receptor occupancy (0–1)
Treg_frac    Regulatory T cell fraction in TME
IFNg_TME     IFN-gamma in tumor microenvironment (relative)
LDH_ser      Serum LDH (U/L) ← tumor necrosis surrogate
S100B_ser    S100B protein (µg/L) ← melanocyte burden marker
```

### Key PD Mechanisms

```
BRAFi → BRAF V600E inhibition (Emax, IC50~0.28 µg/mL, Hill=1.8)
      → ↓MEK → ↓ERK → ↑BIM → BAX/BAK → apoptosis → tumor shrinkage
      → acquired resistance emerges (kR_on * drug_pressure * (1-R))

MEKi → MEK1/2 allosteric inhibition (Emax, IC50~12 ng/mL)
     → ↓ERK → cooperative with BRAFi (Loewe combination model)
     → blocks CRAF paradox activation (monotherapy BRAFi paradox)

BRAFi + MEKi → synergistic ERK suppression → delayed resistance
              → COMBI-d: mPFS 9.3mo vs 5.3mo (Vmfnb mono)

PD-1 blockade → PD-1 receptor occupancy → ↑CD8 TIL recruitment
              → CD8 TIL kill tumor (kd_immune * CD8_eff)
              → IFN-γ → PD-L1↑ (adaptive resistance, kPDL1_IFNg)

CTLA-4 blockade → ↓Treg fraction → ↑CD28 costimulation
                → combined with PD-1 blockade: additive/synergistic

Acquired resistance: BRAF amplification, secondary NRAS Q61 mut,
                     MEK1/2 gain-of-function, COT1 overexpression
                     → ERK reactivation → tumor regrowth
```

### Treatment Scenarios

| # | Scenario | Population | Key Clinical Data |
|---|----------|-----------|------------------|
| 1 | No treatment (untreated) | BRAF V600E met | Median OS ~8mo historically |
| 2 | Vemurafenib 960mg BID | BRAF V600E | BRIM-3: ORR 48%, mPFS 5.3mo |
| 3 | Dabrafenib 150mg BID + Trametinib 2mg QD | BRAF V600E | COMBI-d: ORR 67%, mPFS 9.3mo, 5yr OS 28% |
| 4 | Pembrolizumab 200mg q3w | All subtypes | KEYNOTE-006: ORR 33%, mPFS 5.6mo, 3yr OS 50% |
| 5 | Nivolumab 1mg/kg + Ipilimumab 3mg/kg q3w×4 | All subtypes | CheckMate 067: ORR 58%, mPFS 11.5mo, 5yr OS 52% |
| 6 | Dabrafenib+Trametinib → Pembrolizumab (wk 24) | BRAF V600E | Sequential strategy, KEYNOTE-022 |

---

## Shiny App Features (6 Tabs)

### Tab 1: Patient Profile
- BRAF mutation status selector (V600E / V600K / NRAS / NF1 / WT)
- Clinical stage (IIIB–IVC) and ECOG performance status
- Baseline LDH gauge (prognostic indicator)
- TMB and PD-L1 expression inputs
- Risk stratification summary table

### Tab 2: Drug PK
- BRAF inhibitor plasma concentration (vemurafenib/dabrafenib)
- MEK inhibitor plasma concentration (trametinib/cobimetinib)
- ICI plasma concentration profile (q3w peaks/troughs)
- ERK inhibition profile over time
- PD-1 and CTLA-4 receptor occupancy

### Tab 3: Tumor Response
- Tumor burden (SLD %) over time with CR/PR/SD/PD thresholds
- Best %change gauge at Week 12
- Waterfall plot
- Acquired resistance emergence curve

### Tab 4: Immune Dynamics
- CD8+ TIL dynamics (recruitment boosted by ICI)
- Treg fraction in TME (depleted by CTLA-4 blockade)
- PD-1 receptor occupancy trajectory
- CTLA-4 receptor occupancy trajectory
- IFN-gamma in TME (adaptive PD-L1 induction)
- Immune summary table (weeks 0, 4, 12, 24, 52)

### Tab 5: Scenario Comparison
- Up to 6 treatment arms with selectable checkboxes
- Tumor burden comparison over time
- Best %change waterfall (Week 24)
- Response summary table (tumor %, CD8 TIL, resistance)

### Tab 6: Biomarkers
- Serum LDH trajectory (with ULN reference line)
- S100B dynamics (melanoma burden marker)
- TMB impact on immunotherapy response (sensitivity analysis, 5 TMB levels)
- PD-L1 expression impact on ICI response (5 PD-L1% levels)
- Comprehensive biomarker time-course table

---

## Key References

| Category | Key Papers |
|----------|-----------|
| BRAF mutation | Davies 2002 (*Nature*) — discovery of BRAF V600E |
| BRAFi mono | Chapman 2011 (*NEJM*, BRIM-3) — vemurafenib ORR 48% |
| BRAFi+MEKi | Long 2014 (*NEJM*, COMBI-d) — dab+tram mPFS 9.3mo |
| Encorafenib+binimetinib | Dummer 2018 (*Lancet Oncol*, COLUMBUS) |
| PD-1 mono | Robert 2015 (*NEJM*, KEYNOTE-006) — pembrolizumab |
| CTLA-4 | Hodi 2010 (*NEJM*, MDX010-20) — ipilimumab |
| Dual checkpoint | Larkin 2015 (*NEJM*, CheckMate 067) — nivo+ipi |
| 5-year OS | Wolchok 2019 (*NEJM*) — nivo+ipi 5yr OS 52% |
| Resistance | Nazarian 2010 (*Nature*) — RTK/NRAS resistance |
| TME | Tumeh 2014 (*Nature*) — CD8+ TIL predicts PD-1 response |
| Biomarkers | Ascierto 2017 (*Eur J Cancer*) — LDH prognostic |

Full reference list: **[melanoma_references.md](melanoma_references.md)** (63 PubMed citations)

---

## Clinical Significance

Melanoma is a paradigmatic QSP target because:

1. **Dual therapeutic strategy**: targeted therapy (BRAFi/MEKi) addresses oncogene addiction; immunotherapy exploits immune surveillance
2. **Resistance biology**: acquired BRAFi resistance is mechanistically well-characterized (BRAF amplification, secondary NRAS mutation, MEK1/2 mutation, COT1 bypass) — ideal for ODE modeling
3. **Immune biomarkers**: CD8+ TIL density, TMB, and PD-L1 expression are quantifiable ICI predictors
4. **Sequential opportunity**: BRAFi/MEKi can reshape the TME (reduce Treg, increase MHC-I) before ICI — captured in Scenario 6
5. **Rich clinical trial data**: BRIM-3, COMBI-d/v, CheckMate 067, KEYNOTE-006 provide robust calibration anchors

The QSP model captures the critical interplay:
- **ERK inhibition by BRAFi**: direct tumor kill but also paradox CRAF activation in WT cells
- **MEKi synergy with BRAFi**: blocks CRAF paradox, delays resistance, additional ERK suppression
- **IFN-γ feedback**: CD8 TIL → IFN-γ → PD-L1 upregulation = adaptive resistance to ICI
- **CTLA-4 and Treg**: ipilimumab depletes intra-tumoral Treg → removes CD8 brake

---

*Model created: 2026-06-25 | QSP Disease Library*
