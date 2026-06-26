# Chronic Myeloid Leukemia (CML) — QSP Model

> **First Oncology Model** in the QSP Disease Library — 2026-06-20

[![CML QSP Model Map](cml_qsp_model.png)](cml_qsp_model.svg)

---

## Disease Overview

**Chronic Myeloid Leukemia (CML)** is a clonal myeloproliferative neoplasm defined by the **Philadelphia chromosome** — a t(9;22)(q34;q11) translocation generating the *BCR-ABL1* fusion oncogene. The resulting **p210 BCR-ABL** oncoprotein is a constitutively active tyrosine kinase that drives uncontrolled proliferation and survival of myeloid progenitors.

### Epidemiology
| Parameter | Value |
|-----------|-------|
| Incidence | ~1-2 per 100,000 / year |
| Prevalence | ~70,000 in US (rising due to TKI success) |
| Median age at diagnosis | 57 years |
| Male:Female ratio | 1.4:1 |
| 10-year OS (imatinib era) | ~85% |
| Phase distribution | CP 90% · AP 5% · BC 5% |

### CML Phase Criteria (EWI 2022)
| Phase | Blast % | Key Features |
|-------|---------|--------------|
| Chronic (CP) | <10% | Splenomegaly, WBC ↑, Ph+ |
| Accelerated (AP) | 10–19% | Basophils ≥20%, Plt <100×10⁹/L |
| Blast Crisis (BC) | ≥20% | Lymphoid (30%) or Myeloid (70%) |

---

## Pathophysiology

```
t(9;22) Translocation
    ↓
BCR-ABL1 Fusion Gene → p210 Protein (tyrosine kinase)
    ↓
Constitutive Signaling: RAS-MAPK · PI3K-AKT-mTOR · JAK-STAT · MYC
    ↓
Leukemic Stem Cell (LSC) Self-Renewal ↑ · Apoptosis ↓ · Quiescence ↑
    ↓
Clonal Expansion: LSC → LPC → Leukemic Blasts → WBC ↑↑
    ↓
Clinical Manifestations: Splenomegaly · Fatigue · Elevated WBC
```

---

## Key Signaling Pathways

| Pathway | Effectors | Oncogenic Effect |
|---------|-----------|-----------------|
| RAS-MAPK | GRB2 → SOS1 → RAS → RAF → MEK → ERK | Proliferation, MYC activation |
| PI3K-AKT-mTOR | PI3K → PIP3 → AKT → mTORC1/2 | Survival, protein synthesis |
| JAK-STAT | JAK2 → STAT5 | BCL-2/BCL-xL anti-apoptosis |
| MYC/Cyclin | MYC → CDK4/6 → Cyclin D1 | Cell cycle G1/S progression |

---

## Leukemic Hierarchy (Michor 2005 Framework)

```
Quiescent LSC (x0q) ←→ Cycling LSC (x0p)
                              ↓ differentiation
                    Leukemic Progenitor 1 (x1)
                              ↓
                    Leukemic Progenitor 2 (x2)
                              ↓
                    Leukemic WBC/Blasts (x3)   ← BCR-ABL IS% numerator
```

**TKI sensitivity by compartment:**
- x0q (quiescent LSC): TKI-REFRACTORY (primary resistance mechanism)
- x0p (cycling LSC): TKI-sensitive
- x1, x2 (LPC): partially sensitive
- x3 (WBC): highly sensitive → rapid clearance

---

## TKI Drug Parameters

| TKI | Generation | Standard Dose | IC50 (BCR-ABL) | Bioavailability | t½ | T315I |
|-----|-----------|--------------|----------------|-----------------|-----|-------|
| Imatinib | 1st | 400 mg/day | 0.25 µM | 98% | 18h | Resistant |
| Dasatinib | 2nd | 100 mg/day | 0.025 µM | 14-34% | 5-6h | Resistant |
| Nilotinib | 2nd | 300 mg BID | 0.020 µM | 30-40% | 17h | Resistant |
| Bosutinib | 2nd | 400 mg/day | 0.10 µM | ~34% | 22-27h | Resistant |
| Ponatinib | 3rd | 15-45 mg/day | 0.37 nM | ~54% | 24h | **Sensitive** |
| Asciminib | 3rd STAMP | 40 mg BID | 0.5 nM | ~72% | 8.7h | 200mg BID |

**OCT1 pharmacology (imatinib-specific):**
- OCT1 (SLC22A1) mediates ~80% of intracellular imatinib uptake
- Low OCT1 activity → reduced IC → inferior response
- Dasatinib, nilotinib: OCT1-independent (passive diffusion)

---

## Clinical Response Milestones (ELN 2020)

| Timepoint | Optimal Response | Warning | Failure |
|-----------|-----------------|---------|---------|
| 3 months | BCR-ABL ≤10% | BCR-ABL >10% | BCR-ABL >10% + Ph >95% |
| 6 months | BCR-ABL ≤10% | — | BCR-ABL >10%, Ph >35% |
| 12 months | BCR-ABL ≤1% (CCyR) | BCR-ABL 1-10% | BCR-ABL >1%, Ph >0% |
| Any time | MMR (≤0.1%) | — | Loss of CCyR/CHR |

### Molecular Response Definitions
- **CHR**: Complete Hematologic Response — WBC <10×10⁹/L, no immature cells
- **CCyR**: Complete Cytogenetic Response — 0% Ph+ metaphases
- **MMR/MR3**: BCR-ABL IS ≤0.1%
- **MR4**: BCR-ABL IS ≤0.01%
- **MR4.5**: BCR-ABL IS ≤0.0032% (TFR prerequisite)

---

## Treatment-Free Remission (TFR)

**Prerequisites (ELN 2020):**
- CP CML, no history of AP/BC
- Frontline TKI treatment ≥ 5 years
- Sustained MR4.5 for ≥ 2 years (or MR4 ≥ 3 years)
- BCR-ABL IS quantifiable assay (sensitivity ≥MR4.5)

**TFR Success Rates (Clinical Trials):**
| Trial | TKI | TFR Rate @ 12mo |
|-------|-----|----------------|
| STIM1 (Mahon 2010) | Imatinib | 38% |
| EURO-SKI (Saussele 2018) | Mixed | 43% |
| ENESTfreedom (Hochhaus 2017) | Nilotinib | 51% |
| DADI (Imagawa 2015) | Dasatinib | 49% |

**Molecular relapse after TFR:** ~50-60%; >90% respond to TKI re-initiation

---

## Resistance Mechanisms

| Mechanism | Frequency | Implications |
|-----------|-----------|--------------|
| BCR-ABL kinase domain mutations | ~50% | T315I (gatekeeper): all 1st/2nd-gen TKIs |
| T315I | ~15% of mutations | Ponatinib or asciminib 200mg BID |
| BCR-ABL amplification/overexpression | ~10% | Dose escalation, switch TKI |
| LSC quiescence (primary resistance) | All patients | Combination strategies |
| OCT1 downregulation | Variable (imatinib) | Switch to dasatinib/nilotinib |
| ABCB1/MDR1 overexpression | ~10% | Switch TKI |
| Bypass pathways (LYN, FGF2) | Rare | Src inhibitor (dasatinib) |

---

## Model Files

| File | Description |
|------|-------------|
| [`cml_qsp_model.dot`](cml_qsp_model.dot) | Graphviz mechanistic map (100+ nodes, 12 clusters) |
| [`cml_qsp_model.svg`](cml_qsp_model.svg) | SVG vector graphic (interactive) |
| [`cml_qsp_model.png`](cml_qsp_model.png) | PNG thumbnail (150 dpi) |
| [`cml_mrgsolve_model.R`](cml_mrgsolve_model.R) | mrgsolve ODE model (22 compartments, 7 scenarios) |
| [`cml_shiny_app.R`](cml_shiny_app.R) | Shiny interactive dashboard (6 tabs) |
| [`cml_references.md`](cml_references.md) | 60 PubMed-linked references |

---

## ODE Model Structure (mrgsolve)

**22 Compartments:**

| Compartment | Symbol | Description |
|-------------|--------|-------------|
| y0 | HSC | Normal hematopoietic stem cells |
| y1 | CMP | Common myeloid progenitor |
| y2 | GMP | Granulocyte-macrophage progenitor |
| y3 | Neutrophil | Mature normal WBC |
| x0q | LSC_q | Quiescent CML stem cells (TKI-refractory) |
| x0p | LSC_p | Cycling CML stem cells (TKI-sensitive) |
| x1 | LPC1 | Leukemic progenitor committed |
| x2 | LPC2 | Leukemic progenitor differentiated |
| x3 | WBC_leuk | Mature leukemic cells/blasts |
| xs | LSC_sens | TKI-sensitive LSC clone (BCR-ABL WT) |
| xsp | LPC_sens | Sensitive progenitors |
| xr | LSC_res | T315I-resistant LSC clone |
| xrp | LPC_res | Resistant progenitors |
| Gut_imt | — | Imatinib GI absorption |
| Cp_imt | — | Imatinib plasma (µM) |
| Cic_imt | — | Imatinib intracellular (OCT1-dependent, µM) |
| Cp_das | — | Dasatinib plasma (µM) |
| Cic_das | — | Dasatinib intracellular (µM) |
| Gut_nil | — | Nilotinib GI absorption |
| Cp_nil | — | Nilotinib plasma (µM) |
| Cic_nil | — | Nilotinib intracellular (µM) |
| Cp_asc | — | Asciminib plasma (µM) |

---

## Clinical Trial Calibration

| Trial | TKI | Key Endpoint | Observed | Model (12mo) |
|-------|-----|--------------|----------|--------------|
| IRIS | Imatinib 400mg | MMR @ 12mo | 28% | ~25-32% |
| IRIS | Imatinib 400mg | MMR @ 18mo | 39% | ~35-42% |
| ENESTnd | Nilotinib 300mg BID | MMR @ 12mo | 44% | ~40-48% |
| DASISION | Dasatinib 100mg | MMR @ 12mo | 46% | ~42-50% |
| ASCEMBL | Asciminib 40mg BID | MMR @ 24mo | 25% | ~22-28% |

---

## Shiny App Features

| Tab | Contents |
|-----|---------|
| 환자 프로파일 | Sokal/EUTOS score calculator, ELN treatment recommendation |
| 약동학 (PK) | Plasma & intracellular TKI concentration profiles |
| 분자반응 | BCR-ABL IS% trajectory, ELN milestone checker |
| 세포역학 | Normal vs leukemic cell population dynamics |
| 시나리오 비교 | Head-to-head TKI comparison (MMR/MR4/MR4.5 rates) |
| 내성/TFR 평가 | T315I clone kinetics, TFR simulation, mutation table |

---

## Key References

1. Michor F, et al. (2005). **Dynamics of chronic myeloid leukaemia.** *Nature* 435:1267. [PMID: 15988530]
2. O'Brien SG, et al. (2003). **IRIS trial.** *NEJM* 348:994. [PMID: 12637609]
3. Saglio G, et al. (2010). **ENESTnd trial.** *NEJM* 362:2251. [PMID: 20525993]
4. Kantarjian H, et al. (2010). **DASISION trial.** *NEJM* 362:2260. [PMID: 20525994]
5. Réa D, et al. (2021). **ASCEMBL trial.** *Blood* 138:2031. [PMID: 34407542]
6. Wylie AA, et al. (2017). **Asciminib mechanism.** *Nature* 543:733. [PMID: 28329763]
7. Mahon FX, et al. (2010). **STIM trial (TFR).** *Lancet Oncol* 11:1029. [PMID: 20965785]
8. Hochhaus A, et al. (2020). **ELN 2020 Guidelines.** *Leukemia* 34:966. [PMID: 32127639]

See [`cml_references.md`](cml_references.md) for complete 60-reference bibliography.
