# Acute Lymphoblastic Leukemia (ALL) — QSP Model

Comprehensive Quantitative Systems Pharmacology (QSP) library entry for
**precursor B-/T-cell acute lymphoblastic leukemia** — covering
cytogenetic leukemogenesis (BCR-ABL1, ETV6-RUNX1, KMT2A-r, hyperdiploidy,
Ph-like, NOTCH1), bone-marrow niche adhesion, glucocorticoid apoptosis
signaling, and nine therapeutic agents spanning conventional chemotherapy,
tyrosine-kinase inhibition (Ph+ ALL), CD19×CD3 bispecific T-cell engagement
(blinatumomab), anti-CD22 antibody-drug conjugate (inotuzumab ozogamicin),
and CD19 CAR-T cell therapy (tisagenlecleucel).

## Files

| File | Purpose |
| ---- | ------- |
| [`all_qsp_model.dot`](all_qsp_model.dot) | Mechanistic map source (137 nodes, 15 clusters) |
| [`all_qsp_model.svg`](all_qsp_model.svg) | Vector render |
| [`all_qsp_model.png`](all_qsp_model.png) | 150-dpi PNG render |
| [`all_mrgsolve_model.R`](all_mrgsolve_model.R) | 27-compartment mrgsolve ODE model + 10 scenarios |
| [`all_shiny_app.R`](all_shiny_app.R) | 8-tab Shiny dashboard |
| [`all_references.md`](all_references.md) | 40+ PubMed-linked references |

## Mechanistic map (15 clusters)

1. Lymphoid progenitor biology & leukemogenic drivers — BCR-ABL1,
   ETV6-RUNX1, KMT2A-r, hyperdiploidy/hypodiploidy, Ph-like, TCF3-PBX1,
   NOTCH1/TLX1/TLX3 (T-ALL)
2. Bone-marrow microenvironment & niche adhesion — CXCL12/CXCR4,
   VCAM-1/VLA-4, hypoxia, niche-mediated chemoresistance
3. Survival signal transduction — PI3K/AKT/mTOR, JAK/STAT5, RAS/MAPK
4. Apoptosis resistance & glucocorticoid (NR3C1) signaling — BCL-2/MCL-1,
   BIM/NUR77, steroid resistance
5. Surface immunophenotype — CD19, CD20, CD22, CD3 anchor, antigen-loss/
   lineage-switch immune escape
6. Vincristine — microtubule inhibitor PK/PD, peripheral neuropathy
7. Dexamethasone/prednisone — glucocorticoid PK/PD, hyperglycemia,
   osteonecrosis
8. PEG-asparaginase — enzymatic asparagine depletion, ASNS resistance,
   pancreatitis/thrombosis/hepatotoxicity
9. Methotrexate (HD-MTX + IT) — antifolate PK/PD, CSF penetration,
   leucovorin rescue, CNS sanctuary
10. 6-Mercaptopurine (maintenance) — TPMT/NUDT15 pharmacogenomic TGN
    formation
11. TKI (dasatinib/ponatinib) — BCR-ABL1 inhibition, T315I resistance,
    Ph+ ALL
12. Blinatumomab — CD19×CD3 BiTE immunologic synapse, polyclonal T-cell
    activation
13. Inotuzumab ozogamicin & CD19 CAR-T — ADC payload release, antigen-
    driven CAR-T expansion/contraction/persistence
14. CRS & ICANS — IL-6/IFN-γ cytokine release, BBB permeability,
    tocilizumab/steroid rescue
15. Myelosuppression & clinical endpoints — Friberg ANC/platelet chain,
    MRD, CR, relapse, EFS/OS

## mrgsolve model — 27 ODE compartments

* PK chains: vincristine (2-cmt), dexamethasone (1-cmt oral), PEG-
  asparaginase (enzyme-activity depot), methotrexate (3-cmt incl. CSF),
  6-MP → TGN (pharmacogenomic-scaled), dasatinib/ponatinib (1-cmt oral),
  inotuzumab ozogamicin (1-cmt), blinatumomab (1-cmt continuous infusion)
* CD19 CAR-T cell kinetics: antigen-driven expansion/contraction (blood +
  marrow/tissue trafficking), IL-6 cytokine-release coupling
* Leukemic blast burden: logistic growth minus combined multi-agent
  log-kill (Hill/Emax per agent); CNS sanctuary-site compartment killed
  only by CSF-penetrant MTX (± partial dasatinib)
* MRD (log10 scale), CR/MRD-negativity flags
* Friberg semi-mechanistic myelosuppression: 4-compartment ANC transit
  chain + simplified platelet compartment

### Pre-built scenarios

| ID | Regimen | Source |
| -- | ------- | ------ |
| 1 | Untreated (natural history) | — |
| 2 | Pediatric SR B-ALL induction (VCR+DEX+PEG-ASP+IT-MTX) | COG/DFCI-style |
| 3 | Pediatric HR B-ALL + HD-MTX consolidation | COG-style |
| 4 | Adult Ph-negative B-ALL (hyper-CVAD-like) | Larson CALGB 8811-style |
| 5 | Ph+ ALL + Dasatinib + steroid | GIMEMA LAL1509-style |
| 6 | Ph+ ALL (T315I) + Ponatinib | hyper-CVAD+ponatinib-style |
| 7 | R/R B-ALL + Blinatumomab | TOWER |
| 8 | R/R B-ALL + Inotuzumab ozogamicin | INO-VATE |
| 9 | R/R B-ALL + CD19 CAR-T | ELIANA (tisagenlecleucel) |
| 10 | Maintenance 6-MP, TPMT poor- vs normal-metabolizer | CPIC-style |

## Shiny dashboard — 8 tabs

1. Patient / disease profile
2. Drug PK (cytotoxic/targeted + immunotherapy exposure)
3. Disease biology (BM blast %, MRD, CNS sanctuary burden)
4. Hematologic toxicity (Friberg ANC, platelets)
5. Immunotherapy & CRS (IL-6, CD19 CAR-T cell kinetics)
6. Clinical endpoints (CR, MRD-negativity, relapse risk)
7. Scenario comparison (10 regimens)
8. Biomarkers / safety (plasma asparagine, TGN)

## How to render

```bash
dot -Tsvg all_qsp_model.dot -o all_qsp_model.svg
dot -Tpng -Gdpi=150 all_qsp_model.dot -o all_qsp_model.png
```

```r
source("all_mrgsolve_model.R")
shiny::runApp("all_shiny_app.R")
```

## Notes

* Time unit throughout is **days**; clinical PK parameters typically
  reported in L/h were converted to L/day so CL/V ratios stay
  dimensionally consistent with the day-based dosing/leukemic-kinetics
  axis.
* CAR-T cell kinetics follow the antigen-driven expansion/contraction
  structure described by Stein AM et al. (CPT:PSP 2019) for
  tisagenlecleucel; CRS (IL-6) is coupled to both CAR-T tissue
  activation and blinatumomab-driven T-cell engagement.
* Growth-rate (`k_grow`) and baseline blast burden (`BM_blast_init`) are
  exposed as parameters so relapsed/refractory or MRD-level (low-burden)
  starting conditions can be explored without editing compartment code.
* This is an illustrative, mechanism-forward teaching/exploration model,
  not a clinically validated dosing tool.
