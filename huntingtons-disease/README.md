# Huntington's Disease — QSP Model

[![Mechanistic Map](hd_qsp_model.png)](hd_qsp_model.svg)

## Overview

| Item | Detail |
|------|--------|
| **Disease** | Huntington's Disease (HD) |
| **Category** | Neurodegenerative / Genetic |
| **OMIM** | [143100](https://omim.org/entry/143100) |
| **Gene** | *HTT* (Chromosome 4p16.3) |
| **Mechanism** | CAG repeat expansion → polyQ mHTT → aggregation → striatal MSN degeneration |
| **Key Pathways** | mHTT aggregation · BDNF deficit · Excitotoxicity · Mito dysfunction · Neuroinflammation |
| **Treatments Modeled** | TBZ · DTBZ · VBZ (VMAT2 inhibitors) · Tominersen (ASO) · Branaplam · PTC518 |

---

## Disease Mechanism

Huntington's disease is caused by an autosomal dominant CAG trinucleotide repeat expansion (≥36 repeats) in exon 1 of the *HTT* gene, encoding an expanded polyglutamine (polyQ) tract in the huntingtin protein (mHTT). Key pathophysiological events:

1. **mHTT Production**: CAG length-dependent increase in mHTT protein levels and misfolding propensity
2. **Aggregation Cascade**: Monomer → oligomers (toxic) → protofibrils → intranuclear inclusions
3. **Transcriptional Dysregulation**: mHTT sequesters REST/NRSF → BDNF gene repression; CBP sequestration → HDAC overactivity
4. **BDNF Deficit**: Reduced cortical BDNF production + impaired anterograde transport (HAP1/mHTT) → loss of TrkB-mediated striatal MSN survival
5. **Excitotoxicity**: Sensitized extrasynaptic NMDARs → Ca²⁺ overload → calpain/nNOS activation → mitochondrial failure
6. **Striatal Degeneration**: D2-MSN (indirect pathway) degenerate first → loss of GPe inhibition → **chorea** (early); D1-MSN later → **bradykinesia/dystonia** (late)
7. **Neuroinflammation**: mHTT activates microglia TLR4/NF-κB → IL-1β/TNF-α → worsened excitotoxicity

---

## Files

| File | Description |
|------|-------------|
| `hd_qsp_model.dot` | Graphviz mechanistic map (14 clusters, 110+ nodes) |
| `hd_qsp_model.svg` | Rendered SVG (vector, scalable) |
| `hd_qsp_model.png` | Rendered PNG (150 dpi) |
| `hd_mrgsolve_model.R` | mrgsolve ODE model (20 compartments, 7 scenarios) |
| `hd_shiny_app.R` | Interactive Shiny dashboard (6 tabs) |
| `hd_references.md` | 53 curated PubMed references |

---

## Mechanistic Map — Cluster Summary

| Cluster | Components | Key Biology |
|---------|-----------|-------------|
| 1. Genetic Basis | HTT gene, CAG repeat, mHTT mRNA, polyQ, mHTT protein | CAG-length dependent mHTT production |
| 2. mHTT Aggregation | Monomer, oligomers, fibrils, inclusions, chaperones, ubiquitin | Nucleation-elongation aggregation cascade |
| 3. Protein Quality Control | 26S Proteasome, UPS impairment, autophagy (LC3/Beclin-1/p62), lysosome, TFEB | mHTT clearance mechanisms |
| 4. Transcriptional Dysregulation | REST/NRSF, CREB, Sp1, CBP/p300, HDACs, PGC-1α, BDNF gene | mHTT-driven gene suppression |
| 5. BDNF-TrkB Signaling | Cortical BDNF, HAP1, axonal transport, TrkB, PI3K/Akt, MAPK/ERK | Neurotrophic deficit |
| 6. Mitochondrial Dysfunction | Complex I/II/III, ROS, oxidative stress, ΔΨm, cytochrome c, PGC-1α | Bioenergetic failure |
| 7. Excitotoxicity | NMDAR/AMPAR, eNMDAR, Ca²⁺ influx, nNOS, NO, calpain, mGluR5/IP3R | Glutamate-mediated toxicity |
| 8. Striatal Circuitry | D1/D2-MSN, VMAT2, dopamine, GABA direct/indirect, enkephalin, substance P | Basal ganglia circuit dysfunction |
| 9. Neuroinflammation | Microglia M1/M2, TLR4, NF-κB, IL-1β/TNF-α/IL-6, NLRP3, complement | Neuroinflammatory amplification |
| 10. Apoptosis | BAX/BAK, BCL-2, caspase-3/6, apoptosome, p53, PUMA, AIF | Neuronal death pathways |
| 11. Drug PK | TBZ, DTBZ, VBZ, Tominersen (IT), Branaplam, Riluzole | Multi-drug pharmacokinetics |
| 12. Symptomatic PD | VMAT2 inhibition, DA depletion, chorea reduction | VMAT2-targeted therapy |
| 13. Disease-Modifying PD | ASO (tominersen), splicing (branaplam/PTC518), HDAC inhibitors, CSF biomarkers | HTT-lowering therapies |
| 14. Clinical Endpoints | UHDRS-TMS, TFC, cUHDRS, CAP score, caudate atrophy, CSF NfL/mHTT | Measurable disease burden |

---

## mrgsolve ODE Model

**20 compartments** across PK and PD modules:

### PK Compartments (10)
| # | Compartment | Description |
|---|-------------|-------------|
| 1 | TBZ_gut | Tetrabenazine GI absorption |
| 2 | TBZ_plasma | TBZ central plasma |
| 3 | HTBZ_brain | α/β-HTBZ active metabolites (brain) |
| 4 | DTBZ_plasma | Deutetrabenazine plasma |
| 5 | DTBZ_brain | DTBZ brain |
| 6 | VBZ_plasma | Valbenazine plasma |
| 7 | VBZ_brain | VBZ brain |
| 8 | tominersen_CSF | Tominersen CSF (intrathecal) |
| 9 | riluzole_plasma | Riluzole plasma |
| 10 | riluzole_brain | Riluzole brain |

### PD Compartments (10)
| # | Compartment | Description |
|---|-------------|-------------|
| 11 | mHTT_mRNA | mHTT mRNA (CAG-length dependent production) |
| 12 | mHTT_prot | mHTT soluble protein (nM) |
| 13 | mHTT_oligo | mHTT oligomers — toxic species (nM) |
| 14 | BDNF_cmt | BDNF level (ng/mL CSF proxy) |
| 15 | dopamine_cmt | Synaptic dopamine (nmol) |
| 16 | MSN_surv | MSN survival (% of baseline) |
| 17 | oxidative_idx | Oxidative stress index (ROS proxy) |
| 18 | neuroinflam_idx | Neuroinflammation index (IL-1β proxy) |
| 19 | UHDRS_TMS | UHDRS Total Motor Score |
| 20 | TFC_cmt | Total Functional Capacity (13→0) |

---

## Treatment Scenarios

| # | Scenario | Drug | Dose & Route | Mechanism |
|---|----------|------|-------------|-----------|
| 1 | Natural History | — | None | Disease progression only |
| 2 | Tetrabenazine | TBZ | 25 mg/d TID (oral) | VMAT2 inhibition → DA depletion → chorea↓ |
| 3 | Deutetrabenazine | DTBZ | 30 mg/d BID (oral) | VMAT2 inhib (d-KIE → less peak Cp) |
| 4 | Valbenazine | VBZ | 80 mg QD (oral) | VMAT2 inhib (long T½, QD dosing) |
| 5 | Tominersen | ASO | 120 mg Q8W (IT) | RNase H1 → mHTT mRNA↓ ~74% |
| 6 | Branaplam | LMI070 | 50 mg QW (oral) | Exon 49 skipping → mHTT NMD ↓~50% |
| 7 | Combination | DTBZ + Tominersen | Both | Symptomatic + disease-modifying |

---

## Clinical Trial Calibration

| Trial | Drug | Key Finding | Reference |
|-------|------|------------|-----------|
| TETRA-HD (2008, NEJM) | Tetrabenazine | Chorea ↓5.0 TMS units vs placebo | Frank et al., NEJM 2008 |
| FIRST-HD (2016, JAMA) | Deutetrabenazine | Chorea ↓2.5 units (TMS); longer tolerability | HSG HART, JAMA 2016 |
| KINECT-HD (2023, NEJM) | Valbenazine | TMS ↓3.2 (placebo-adjusted); QD once-daily | Videnovic, NEJM 2023 |
| GENERATION-HD1 (2022, NEJM) | Tominersen | CSF mHTT ↓74%; motor progression worsened (high dose arm) | Tabrizi, NEJM 2022 |
| NCT04000594 | Branaplam | HTT mRNA ↓~50%; ongoing Phase 2 | Schulte, Ann Neurol 2023 |
| TRACK-HD (2011, Lancet Neurol) | — | TFC −0.7/yr; caudate 2–4 mL/yr atrophy | Tabrizi, Lancet Neurol 2011 |

---

## Shiny App — Tab Structure

| Tab | Content |
|-----|---------|
| 1. Patient Profile | CAG repeat, age, CAP score, HD staging, 10-yr natural history projection |
| 2. Pharmacokinetics | Drug selection (TBZ/DTBZ/VBZ), plasma & brain PK curves, PK parameter table |
| 3. PD Key Indices | mHTT cascade, oligomers, BDNF, MSN survival, dopamine, ROS, neuroinflammation |
| 4. Clinical Endpoints | UHDRS-TMS, TFC, cUHDRS trajectories with treatment options |
| 5. Scenario Comparison | 7-arm comparison (TMS, TFC, MSN, summary table) |
| 6. Biomarkers | CSF NfL, CSF mHTT, oligomer fraction, target engagement, biomarker-clinical correlation |

---

## Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| CAG onset threshold | ≥36 repeats | Langbehn 2004 |
| mHTT mRNA T½ | ~15 hr | Model calibration |
| UPS saturation Km | 120 nM | Hanser 2021 |
| BDNF healthy (CSF) | 6.5 ng/mL | Zuccato 2009 |
| MSN at manifest onset | ~85% surviving | Vonsattel 1985 |
| TFC decline rate | −0.7 units/yr | ENROLL-HD 2019 |
| TMS progression | +2.5 units/yr | TRACK-HD |
| HTBZ VMAT2 IC50 | 0.012 μM | Frank 2008 |
| Tominersen Emax | 74% mHTT↓ | GENERATION-HD1 |
| CAG-TMS correlation (r) | 0.60 | Snell 1993 |

---

## References

See [`hd_references.md`](hd_references.md) for 53 curated PubMed references organized by:
Genetic Basis · mHTT Aggregation · Transcription · Excitotoxicity · Mitochondria · Striatal Pathology · Neuroinflammation · Clinical Scales · Symptomatic Treatment · HTT-Lowering Therapies · Biomarkers · QSP Modeling · Natural History Studies · Neuroprotective Strategies · Guidelines

---

*Generated 2026-06-23 | QSP Disease Model Library*
