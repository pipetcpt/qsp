# Familial Hypercholesterolemia (FH) — QSP Model

> **가족성 고콜레스테롤혈증 (Familial Hypercholesterolemia)**  
> LDLR/APOB/PCSK9 돌연변이로 인한 LDL 수용체 기능 저하 → LDL-C 극적 상승 → 조기 심혈관 질환

---

## Disease Overview

Familial Hypercholesterolemia (FH) is a monogenic autosomal dominant disorder of cholesterol metabolism caused by pathogenic variants in **LDLR** (∼85%), **APOB** (∼5–10%), or **PCSK9** (∼1–2%) genes. It is the most common single-gene disorder predisposing to premature atherosclerotic cardiovascular disease (ASCVD).

| Parameter | Heterozygous FH | Homozygous FH |
|-----------|-----------------|---------------|
| Prevalence | ~1/250–500 | ~1/300,000 |
| Baseline LDL-C | 190–400 mg/dL | 400–1,000 mg/dL |
| LDLR Function | ~50% | 0–25% |
| Premature MI | Males <55 yr, Females <65 yr | Often <20 yr |
| Xanthomas | Tendon xanthomas | Extensive xanthomas |

---

## Mechanistic Map

[![FH QSP Mechanistic Map](fh_qsp_model.png)](fh_qsp_model.svg)

> Click the image to view the full interactive SVG map.

**Map Statistics:**
- **11 clusters** (pathway modules)
- **157+ nodes** (biological entities, drugs, enzymes)
- **Pathways covered:**
  - Statin PK (CYP3A4/2C9, OATP1B1, UGT1A3)
  - PCSK9 inhibitor & inclisiran PK
  - Other lipid-lowering drug PK (ezetimibe, bempedoic acid, lomitapide)
  - Mevalonate/cholesterol biosynthesis (22 steps)
  - LDLR–PCSK9 biology (SREBP2, SCAP, INSIG, clathrin endocytosis)
  - Plasma lipoprotein metabolism (VLDL→IDL→LDL, HDL reverse transport)
  - Intestinal cholesterol absorption (NPC1L1, ABCG5/G8, FXR-CYP7A1)
  - LDLR genetic variant classes (Class I–V mutations)
  - Vascular atherogenesis (foam cells, NLRP3, plaque rupture)
  - Clinical endpoints (LDL-C, ESC goals, CVD risk)
  - Statin pleiotropic effects (eNOS, Rho, MMP, CRP)

---

## Files

| File | Description |
|------|-------------|
| `fh_qsp_model.dot` | Graphviz source (157+ nodes, 11 clusters) |
| `fh_qsp_model.svg` | Vector mechanistic map (interactive) |
| `fh_qsp_model.png` | Raster map (150 dpi) |
| `fh_mrgsolve_model.R` | ODE-based PK/PD model (18 compartments, 6 scenarios) |
| `fh_shiny_app.R` | Shiny dashboard (6 tabs, interactive) |
| `fh_references.md` | 57 curated PubMed references |

---

## mrgsolve ODE Model

### Compartments (18 ODEs)

**PK Compartments:**
```
GUT_S       Statin GI absorption
CENT_S      Statin central plasma
LIV_S       Statin hepatic (active site)
SC_PCSK9I   Evolocumab subcutaneous depot
CENT_PCSK9I Evolocumab central plasma
PERI_PCSK9I Evolocumab peripheral compartment
COMP_PK9    PCSK9i–PCSK9 complex
GUT_EZE     Ezetimibe GI
CENT_EZE    Ezetimibe plasma
```

**PD Compartments:**
```
HMGCR_rel   HMGCR relative activity (statin target)
LDLR_rel    LDLR surface pool (genotype-scaled + PCSK9 regulated)
PCSK9_pl    Plasma PCSK9 (ng/mL)
VLDL_C      VLDL cholesterol (mg/dL)
IDL_C       IDL cholesterol (mg/dL)
LDL_C       LDL cholesterol (mg/dL) ← primary endpoint
HDL_C       HDL cholesterol (mg/dL)
TG_C        Triglycerides (mg/dL)
```

### Key PD Mechanisms

```
Statin → HMGCR inhibition (Emax model) → ↓intracellular cholesterol
       → SCAP/SREBP2 activation → ↑LDLR transcription
       → ↑LDLR surface → ↑LDL clearance → ↓plasma LDL-C
       → ↑PCSK9 transcription (SREBP2 feedback) [attenuates LDLR↑]

PCSK9i → neutralizes secreted PCSK9
       → prevents PCSK9-mediated LDLR degradation
       → ↑LDLR recycling → sustained ↑LDLR → ↓LDL-C

Inclisiran → hepatic siRNA → RISC → PCSK9 mRNA cleavage
          → ↓PCSK9 synthesis → same downstream effect as PCSK9i

Ezetimibe → NPC1L1 inhibition → ↓intestinal cholesterol absorption
          → ↓VLDL substrate → modest ↓LDL-C (~18%)

Bempedoic acid → ACLY inhibition → ↓acetyl-CoA → ↓HMGCR substrate
              → ↑LDLR (SREBP2) → ↓LDL-C (~18%)

Lomitapide → MTP inhibition → ↓VLDL assembly → ↓LDL-C (HomFH)
```

### Treatment Scenarios

| # | Scenario | Population | Expected ↓LDL-C |
|---|----------|-----------|----------------|
| 1 | No treatment (baseline) | HetFH | 0% |
| 2 | Rosuvastatin 40 mg/d | HetFH | ~45–55% |
| 3 | Rosuvastatin 40 mg + Ezetimibe 10 mg | HetFH | ~60–65% |
| 4 | Evolocumab 420 mg q4w (mono) | HetFH | ~55–60% |
| 5 | Rosuvastatin 40 mg + Evolocumab 420 mg q4w | HetFH | ~70–78% |
| 6 | Lomitapide + Statin + Evolocumab | HomFH | ~50–60% |

---

## Shiny App Features (6 Tabs)

### Tab 1: Patient Profile
- FH genotype selector (HetFH / HomFH / FDB / Normal)
- Baseline lipid panel input
- LDL-C goal attainment gauge
- 10-year CVD risk estimation

### Tab 2: Drug PK
- Statin hepatic concentration profile
- Evolocumab plasma concentration (q4w peaks/troughs)
- Ezetimibe plasma levels
- HMGCR relative activity over time

### Tab 3: Lipid Response
- LDL-C time course with ESC goal lines (55 / 70 mg/dL)
- HDL-C & TG trajectory
- Full lipid panel comparison: baseline vs. Week 52

### Tab 4: Biomarkers
- Total vs. free plasma PCSK9
- Hepatic LDLR surface expression (%)
- Non-HDL-C trajectory
- Biomarker summary table (DataTables)

### Tab 5: Scenario Comparison
- 8-arm scenario comparison (all treatment combos)
- % LDL-C reduction waterfall chart
- ESC goal attainment table (55/70 mg/dL)

### Tab 6: Genetic Risk Profile
- LDLR mutation class (I–V) → LDL-C scatter
- PCSK9 variant (GoF/WT/LoF) → LDL-C impact
- CVD risk reduction by treatment duration (CTT meta-analysis)

---

## Key References

| Category | Key Papers |
|----------|-----------|
| LDLR biology | Brown & Goldstein 1986 (Nobel lecture) |
| PCSK9 discovery | Seidah 2003; Abifadel 2003 |
| Statin meta-analysis | CTT Collaboration 2010 (*Lancet*) |
| FOURIER trial | Sabatine 2017 (*NEJM*) — evolocumab |
| ODYSSEY OUTCOMES | Schwartz 2018 (*NEJM*) — alirocumab |
| ORION inclisiran | Wright 2020; Raal 2020 (*NEJM*) |
| IMPROVE-IT | Cannon 2015 (*NEJM*) — ezetimibe |
| CLEAR Outcomes | Nissen 2023 (*NEJM*) — bempedoic acid |
| HomFH lomitapide | Cuchel 2013 (*Lancet*) |
| ESC guidelines | Mach 2019 (*Eur Heart J*) |

Full reference list: **[fh_references.md](fh_references.md)** (57 PubMed citations)

---

## Clinical Significance

FH represents a critical QSP modeling target because:

1. **Clear genetic mechanism** → receptor activity directly determines LDL-C
2. **Multiple drug classes** with complementary mechanisms (HMGCR, PCSK9, NPC1L1, ACLY, MTP)
3. **Quantifiable biomarkers** (LDL-C, PCSK9, ApoB)
4. **Lifetime cumulative exposure** drives CVD risk (Mendelian randomization evidence)
5. **Gene therapy emerging** (AAV-LDLR, CRISPR-PCSK9 base editing)

The QSP model captures the interplay between:
- **Statin's dual effect**: ↓LDL-C (via LDLR↑) but also ↑PCSK9 (via SREBP2)
- **PCSK9i synergy with statin**: compensates for the PCSK9 feedback
- **Genotype-driven ceiling effects**: HomFH shows blunted statin response

---

*Model created: 2026-06-25 | QSP Disease Library*
