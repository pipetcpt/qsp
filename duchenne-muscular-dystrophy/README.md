# Duchenne Muscular Dystrophy (DMD) QSP Model

> **Directory**: `duchenne-muscular-dystrophy/` | **Date**: 2026-06-23 | **Category**: Neuromuscular Genetic Disease

[![DMD QSP Map](dmd_qsp_model.png)](dmd_qsp_model.svg)

---

## Disease Overview

**Duchenne Muscular Dystrophy (DMD)** is the most common lethal X-linked recessive neuromuscular disease, affecting approximately 1 in 3,500 live male births (~300,000 patients worldwide). It results from mutations in the *DMD* gene (Xp21.2, the largest human gene at 2.4 Mb), leading to complete absence or severe deficiency of dystrophin (427 kDa).

### Pathophysiology Summary

The absence of dystrophin disrupts the Dystrophin-Associated Protein Complex (DAPC), which mechanically links the intracellular F-actin cytoskeleton to the extracellular matrix laminin-211. This causes sarcolemmal fragility → calcium influx → calpain/ROS-mediated necrosis → inflammation (NF-κB/TNF-α pathway) → fibrosis (TGF-β1/SMAD2/3) → progressive loss of muscle function, cardiomyopathy, and respiratory failure.

### Mutation Distribution
| Type | Frequency | Exon-Skip Eligibility |
|------|-----------|----------------------|
| Deletion/Duplication | ~72% | Many skippable |
| Nonsense/PTC | ~13% | Ataluren (PTC124) |
| Point/Other | ~15% | Variable |
| **Exon-51 skippable** | ~13% | Eteplirsen (EXONDYS 51) |
| **Exon-45 skippable** | ~8% | Casimersen (AMONDYS 45) |
| **Exon-53 skippable** | ~8% | Golodirsen (VYONDYS 53) |

### Natural History (Untreated)
| Milestone | Age |
|-----------|-----|
| Gross motor delay, Gower sign | 2–3 yr |
| Calf pseudohypertrophy, waddling gait | 3–6 yr |
| 6MWD peak (~380–400m) | 8–9 yr |
| Loss of ambulation | Median 12–13 yr |
| Respiratory support (NIV) | ~16–20 yr |
| DCM onset (LVEF <55%) | ~10–15 yr |
| Median survival (untreated) | ~26 yr |
| Survival with modern care | >40 yr |

---

## Core Pathophysiological Pathways

| Cluster | Key Molecules | Therapeutic Targets |
|---------|--------------|---------------------|
| **Genetic basis** | DMD exon deletions, nonsense PTC, reading frame disruption | Exon-skipping (ASOs), ataluren, gene therapy |
| **DAPC complex** | α/β-dystroglycan, sarcoglycans, syntrophins, nNOS, laminin-211 | Utrophin upregulation, DAPC restoration |
| **Sarcolemmal Ca²⁺** | TRPC1/6, TRPV2, calpain-1/2/3, mPTP, ROS, ER stress | Ca²⁺ channel blockers, calpain inhibitors |
| **Inflammation** | NF-κB, TNF-α, IL-1β/6, M1/M2 macrophages, NLRP3, complement | Deflazacort/prednisone/vamorolone (GR), anti-TNF |
| **Muscle regeneration** | Satellite cells (Pax7), MyoD, myogenin, IGF-1/mTORC1, Notch/Wnt | SC activation, myostatin inhibition (anti-GDF8) |
| **Fibrosis** | TGF-β1, SMAD2/3, FAPs, myofibroblasts, collagen I/III, MMP-9 | Givinostat (HDAC), anti-TGF-β, losartan |
| **Cardiac (DCM)** | Cp dystrophin deficiency → LVIDd↑, LVEF↓, LGE, NT-proBNP↑ | ACEi (enalapril), ARBs, beta-blockers, eplerenone |
| **Respiratory** | Diaphragm weakness, FVC↓ ~2-4%/yr, PCF↓ | NIV, cough assist, FVC monitoring |
| **Drug PK** | Eteplirsen/casimersen (PMO-ASO IV), DFZ/prednisone (oral), Elevidys (AAV IV) | Dosing optimization |
| **Clinical endpoints** | NSAA (0-34), 6MWD, TLA, FVC%, LVEF, CK, dystrophin IHC | Treatment response biomarkers |

---

## QSP Model Outputs (4 Deliverables)

| Deliverable | File | Description |
|-------------|------|-------------|
| 🗺️ Mechanistic Map | [`dmd_qsp_model.dot`](dmd_qsp_model.dot) / [`.svg`](dmd_qsp_model.svg) / [`.png`](dmd_qsp_model.png) | 160+ nodes, 10 subgraph clusters, 200+ edges covering all DMD pathways |
| ⚙️ mrgsolve ODE | [`dmd_mrgsolve_model.R`](dmd_mrgsolve_model.R) | 18-compartment ODE: 6 drug PK + 12 disease PD; 7 treatment scenarios; calibrated to ESSENCE/EMBARK/VISION-DMD trials |
| 📊 Shiny App | [`dmd_shiny_app.R`](dmd_shiny_app.R) | 7-tab interactive dashboard: Patient Profile · Drug PK · Dystrophin · Motor Function · Clinical Endpoints · Scenario Comparison · Biomarkers |
| 📚 References | [`dmd_references.md`](dmd_references.md) | 61 PubMed citations across 13 sections (epidemiology, genetics, mechanisms, treatments, QSP modeling) |

---

## Model Specifications

### mrgsolve ODE Compartments (18 total)

**Drug PK (6):**
| Compartment | Description |
|-------------|-------------|
| `Ete_C1` | Eteplirsen central (mg) |
| `Ete_C2` | Eteplirsen peripheral (mg) |
| `Ete_Muscle` | Eteplirsen muscle uptake (pmol/g ×1000) |
| `DFZ_Gut` | Deflazacort GI absorption (mg) |
| `DFZ_Plasma` | Deflazacort plasma (mg) |
| `Active_DFZ` | 21-desacetyl-DFZ active metabolite (mg) |

**Disease PD (12):**
| Compartment | Description |
|-------------|-------------|
| `Dystrophin` | Dystrophin level (% of normal, 0-100) |
| `Fiber_H` | Healthy muscle fibers (au) |
| `Fiber_N` | Necrotic fibers (au) |
| `Fiber_R` | Regenerating fibers (au) |
| `Inflam` | Inflammation index (0-100) |
| `Fibrosis` | Fibrosis score (0-1) |
| `SC_Pool` | Satellite cell pool capacity (0-1) |
| `CK_serum` | Serum creatine kinase (U/L) |
| `FVC_pct` | FVC% predicted |
| `LVEF` | Left ventricular ejection fraction (%) |
| `NSAA` | North Star Ambulatory Assessment (0-34) |
| `SixMWD` | 6-minute walk distance (meters) |

### Treatment Scenarios (7)

| # | Scenario | Key Parameters | Clinical Trial |
|---|----------|----------------|----------------|
| 1 | Natural History | No drug | Historical controls |
| 2 | Deflazacort | 0.9 mg/kg/day PO | ESSENCE (Griggs 2016) |
| 3 | Prednisone | 0.75 mg/kg/day PO | CINRG (McDonald 2018) |
| 4 | Eteplirsen + DFZ | 30 mg/kg/wk IV + DFZ | Mendell 2016 Ann Neurol |
| 5 | Casimersen + DFZ | 30 mg/kg/wk IV + DFZ | ESSENCE-DMD45 |
| 6 | Elevidys (Gene Therapy) | 1×10¹⁴ vg/kg IV, single dose | EMBARK Phase 3 (2023) |
| 7 | Vamorolone | 6 mg/kg/day PO | VISION-DMD (Servais 2022) |

---

## Clinical Trial Calibration Data

| Treatment | Trial | Key Result | PMID |
|-----------|-------|-----------|------|
| Deflazacort vs prednisone | ESSENCE (Griggs 2016, Neurology) | DFZ: +2yr ambulation, +23m 6MWD vs pred | [27784946](https://pubmed.ncbi.nlm.nih.gov/27784946/) |
| Deflazacort long-term | McDonald 2018, J Pediatr | +2yr ambulation vs placebo; lean mass ↑ | [29425585](https://pubmed.ncbi.nlm.nih.gov/29425585/) |
| Eteplirsen 30 mg/kg/wk | Mendell 2016, Ann Neurol | Dystrophin 0.28–0.93% vs baseline; 6MWD decline slowed | [26920136](https://pubmed.ncbi.nlm.nih.gov/26920136/) |
| Delandistrogene (Elevidys) | EMBARK Phase 3 (2023) | Dystrophin 28.1% vs 1.7% (placebo) at 52wk; NSAA improvement | FDA approval 2023 |
| Vamorolone 6 mg/kg/d | VISION-DMD (Campbell 2023, Ann Neurol) | +1.88 m/yr 6MWD vs baseline; ↓ bone SE vs pred | [36738106](https://pubmed.ncbi.nlm.nih.gov/36738106/) |
| Givinostat (ongoing) | EPIDYS Phase 3 | Muscle fat fraction ↓ (MRI), Tanner 1–3 boys | NCT02851797 |

---

## Running the Model

```bash
# Render mechanistic map (requires Graphviz)
dot -Tsvg dmd_qsp_model.dot -o dmd_qsp_model.svg
dot -Tpng -Gdpi=150 dmd_qsp_model.dot -o dmd_qsp_model.png
```

```r
# Run mrgsolve simulation
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("dmd_mrgsolve_model.R")

# Interactive Shiny dashboard
install.packages(c("shiny","plotly","bslib","DT"))
shiny::runApp("dmd_shiny_app.R")
```

---

## Key Biological Insights

1. **Dystrophin threshold effect**: Even 3–15% of normal dystrophin can confer significant functional benefit (Becker-like phenotype), making exon-skipping and gene therapy clinically meaningful despite incomplete restoration.

2. **Necrosis-regeneration cycle exhaustion**: Satellite cells undergo ~50 divisions before telomere shortening impairs regenerative capacity. Progressive SC pool depletion accelerates loss of ambulation.

3. **Fibrosis as gatekeeper**: TGF-β1-driven endomysial fibrosis is the primary irreversible driver of functional decline. Anti-fibrotic strategies (givinostat) may be critical adjuncts to dystrophin-restoring therapies.

4. **Cardiac surveillance**: DCM typically lags skeletal muscle disease by 5–7 years. Prophylactic ACEi at age 10yr (even with preserved LVEF) delays cardiac deterioration (TREAT-NMD guidelines).

5. **Dissociated GR modulation (vamorolone)**: By retaining anti-inflammatory transrepression (↓NF-κB) while reducing transactivation of metabolic genes, vamorolone achieves similar motor benefits with better bone/growth safety vs prednisone.

---

*Model generated by Claude Code Routine (CCR) | 2026-06-23*
