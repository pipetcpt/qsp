# Carcinoid Syndrome (Functional Midgut NET) — QSP Model

Comprehensive Quantitative Systems Pharmacology (QSP) library entry for
**functional, well-differentiated midgut neuroendocrine tumor (NET) with
carcinoid syndrome (CS)** — covering hepatic-metastasis-driven 5-HT
spillover, SSTR2/5 pharmacology, TPH1 inhibition (telotristat),
¹⁷⁷Lu-DOTATATE PRRT, mTORC1 (everolimus), VEGFR-TKI, IFN-α, and CHD
(valvular fibrosis) progression.

## Files

| File | Purpose |
| ---- | ------- |
| [`carcsyn_qsp_model.dot`](carcsyn_qsp_model.dot) | Mechanistic map source (130 nodes, 14 clusters) |
| [`carcsyn_qsp_model.svg`](carcsyn_qsp_model.svg) | Vector render |
| [`carcsyn_qsp_model.png`](carcsyn_qsp_model.png) | 150-dpi PNG render |
| [`carcsyn_mrgsolve_model.R`](carcsyn_mrgsolve_model.R) | 27-state mrgsolve ODE model + 12 scenarios |
| [`carcsyn_shiny_app.R`](carcsyn_shiny_app.R) | 8-tab Shiny dashboard |
| [`carcsyn_references.md`](carcsyn_references.md) | 70+ PubMed-linked references |

## Mechanistic map (14 clusters)

1. Enterochromaffin / midgut NET cell — Trp → TPH1 → 5-HTP → AADC → 5-HT
2. NET driver biology — Ki-67, CDKN1B, MEN1, DAXX/ATRX, PI3K/AKT/mTOR
3. Secreted products — 5-HT, CgA, NKA, SP, bradykinin, histamine, PGE2
4. Serotonin clearance — MAO-A / ALDH → urinary 5-HIAA; SERT → platelets
5. Flushing & vasomotor — H1, NK1, B2-bradykinin, NO → vasodilation
6. Secretory diarrhea — 5-HT3/4 → Cl⁻/H₂O secretion + motility
7. Carcinoid Heart Disease — 5-HT2B → TGF-β1 → valve fibrosis → TR/PS
8. Bronchospasm — 5-HT2A, H1 airway tone
9. Carcinoid crisis — anesthesia/PRRT trigger; octreotide prophylaxis
10. SSTR2/5 axis — Gαi → ↓cAMP → ↓exocytosis + antiproliferation
11. Drug PK/PD — octreotide LAR, lanreotide, pasireotide, telotristat,
    PRRT, IFN-α, everolimus, VEGFR-TKI, HAE/TACE/Y90, ondansetron
12. Clinical endpoints — BM/day, flushing/day, urinary 5-HIAA, PFS, OS,
    Hassan score, EORTC GI-NET21 QoL
13. Drug safety — gallstones, hyperglycemia, steatorrhea, cytopenias,
    pneumonitis, depression
14. Patient covariates — age, weight, renal/hepatic function, grade,
    primary site, SSTR PET Krenning, dietary triggers

## mrgsolve model — 27 ODE compartments

* PK chains: octreotide LAR (3-cmt), lanreotide, pasireotide, telotristat
  (parent → active LP-778902), everolimus, VEGFR-TKI lump, IFN-α, PRRT
* Tumor compartment with Gompertz/logistic growth modulated by mTORi,
  VEGFi, IFNi and PRRT-mediated kill
* Serotonin pathway: TPH1 fractional activity → 5-HTP → tumor 5-HT →
  plasma 5-HT → platelet 5-HT, urinary 5-HIAA
* Symptom indirect-response: BM/day (Emax on 5-HT), flushing/day
* CHD axis: 5-HT2B-driven TGF-β1 → valve collagen → NT-proBNP

### Pre-built scenarios

| ID | Regimen | Source |
| -- | ------- | ------ |
| S1 | Natural history (12 mo) | — |
| S2 | Octreotide LAR 30 mg q28d × 12 | PROMID 2009 |
| S3 | Lanreotide autogel 120 mg SC q28d × 12 | CLARINET 2014 |
| S4 | Octreotide LAR + Telotristat 250 mg t.i.d. | TELESTAR/TELECAST |
| S5 | Pasireotide LAR 60 mg q28d (refractory) | Wolin 2015 |
| S6 | Everolimus 10 mg/day | RADIANT-4 |
| S7 | VEGFR-TKI 37.5 mg/day (sunitinib/surufatinib lump) | Raymond 2011 |
| S8 | ¹⁷⁷Lu-DOTATATE 7.4 GBq q8w × 4 | NETTER-1/2 |
| S9 | IFN-α-2b 5 MU SC TIW | Faiss 2003 |
| S10 | Hepatic artery embolization + octreotide | NCCN |
| S11 | Carcinoid crisis prophylaxis (IV octreotide) | Massimino 2013 |
| S12 | Quad therapy (octreotide+telotristat+PRRT+everolimus) | composite |

## Shiny dashboard — 8 tabs

1. Patient / tumor profile (Ki-67, Krenning, hepatic load, site)
2. Drug PK & SSTR2/5 occupancy
3. Tryptophan → 5-HTP → 5-HT → 5-HIAA pathway
4. Symptoms — bowel movements/day, flushing/day
5. Tumor burden & growth modulation
6. Carcinoid heart disease — valve score, NT-proBNP
7. Scenario comparator
8. Biomarker panel (5-HT, urinary 5-HIAA, platelet 5-HT)

## How to render

```bash
dot -Tsvg carcsyn_qsp_model.dot -o carcsyn_qsp_model.svg
dot -Tpng -Gdpi=150 carcsyn_qsp_model.dot -o carcsyn_qsp_model.png
```

```r
source("carcsyn_mrgsolve_model.R")
shiny::runApp("carcsyn_shiny_app.R")
```

## Notes

* Parameter priors were anchored to published trial-level summaries
  (TELESTAR, CLARINET, PROMID, NETTER-1/2, RADIANT-4); IIV/RUV ranges
  reflect typical NET-NLME literature.
* The 5-HT2B → TGF-β → valve collagen module is exploratory and
  intentionally slow-evolving (decade time-scale) per Møller/Hassan data.
* Carcinoid crisis is encoded as an event-driven IV octreotide bolus
  that pre-suppresses release prior to a procedure (S11).
