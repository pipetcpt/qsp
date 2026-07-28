# 위장관 기질종양 (GIST) — QSP 모델
### Gastrointestinal Stromal Tumor — Quantitative Systems Pharmacology Model

> A **polyclonal** QSP model of GIST: six tumor sub-populations under a
> drug × clone effective-IC50 matrix, five tyrosine kinase inhibitors with
> mechanistic oral PK, and the imaging / ctDNA / toxicity endpoints that
> actually drive clinical decisions in this disease.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`gist_qsp_model.dot`](gist_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`gist_qsp_model.svg`](gist_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`gist_qsp_model.png`](gist_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`gist_mrgsolve_model.R`](gist_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`gist_shiny_app.R`](gist_shiny_app.R) |
| 📚 References (46, all PMIDs verified) | [`gist_references.md`](gist_references.md) |

---

## 1. Disease in one paragraph

Gastrointestinal stromal tumor is the most common mesenchymal neoplasm of the
GI tract, arising from **interstitial cells of Cajal** (the gut's pacemaker
cells) or their precursors under an **ETV1** lineage-survival program. In
roughly 85% of cases a single mutually exclusive gain-of-function mutation
drives the disease: **KIT exon 11** (juxtamembrane, ~65-70%, loss of the
autoinhibitory brake), **KIT exon 9** (extracellular dimerization domain,
~8-10%), or **PDGFRA exon 18 D842V** (~5-6%). The remainder are
**SDH-deficient / quadruple wild-type** tumors in which succinate accumulation
inhibits prolyl hydroxylases, stabilizes HIF-1α and drives VEGF — a
KIT-independent, pseudohypoxic biology with an entirely different treatment
logic. GIST was the first solid tumor to be transformed by a targeted kinase
inhibitor, and remains the cleanest human model of **drug-driven clonal
evolution**: imatinib works in almost everyone and stops working in almost
everyone at 20-24 months, not because the drug stops reaching the target, but
because it selects secondary kinase-domain mutants — **polyclonally**, with
different mutations in different metastases of the same patient.

## 2. Why this model has six tumor compartments

Three clinical facts make a single-compartment "tumor size vs drug
concentration" model unable to represent GIST:

1. **Resistance is target mutation, not drug failure.** Secondary mutations in
   the KIT kinase domain arise under therapy and restore signaling at
   achievable drug concentrations.
2. **Resistance is polyclonal, and the two families respond to opposite
   drugs.** ATP-binding-pocket mutants (exon 13 V654A, exon 14 T670I
   gatekeeper) are suppressed by **sunitinib** but poorly by regorafenib;
   activation-loop mutants (exon 17 D816/D820/N822K, Y823D) are the reverse.
3. **Therefore the winner of a head-to-head trial depends on the clonal
   composition of the enrolled population** — which is exactly what the
   INTRIGUE ctDNA analysis found.

So the model carries: the **primary-driver sensitive clone**, a **quiescent
drug-tolerant persister pool**, an **exon 13/14 ATP-pocket clone**, an **exon
17/18 activation-loop clone**, a **primary-refractory clone** (PDGFRA D842V or
SDH-deficient biology, depending on genotype), and a **KIT-independent bypass
clone** — each with its own row in a 5 × 4 drug × clone effective-IC50 matrix.
Mutant supply is proportional to the **division flux of the sensitive clone**,
so a large, actively cycling sensitive population manufactures its own
resistance. That single coupling is what makes dose, trough and genotype
affect *progression timing* and not merely *response depth*.

## 3. Mechanistic map — 14 clusters, 170 nodes

1. **Cell of origin & lineage dependency** — ICC/ICC-precursor, SCF-KIT, ANO1/DOG1, ETV1 master transcription factor and its ERK-dependent protein stabilization, microscopic "tumorlets"
2. **Primary driver mutations** — KIT exon 11/9/13/17, PDGFRA exon 12/14/18 D842V, SDH-deficient, NF1, BRAF V600E, RAS, NTRK/FGFR fusion, quadruple wild-type, Carney triad, germline familial GIST
3. **Receptor activation & conformational states** — juxtamembrane autoinhibition, ligand-independent dimerization, DFG-in/DFG-out, switch pocket, activation loop, ATP pocket, gatekeeper T670, KIT amplification
4. **Downstream signaling** — PI3K/AKT/mTOR/S6K-4EBP1, RAS/RAF/MEK/ERK, JAK-STAT3, SRC, PLCγ, FOXO3a, GLUT1
5. **Cell cycle, apoptosis, autophagy & persistence** — MYC, cyclin D-CDK4/6-Rb, CDKN2A deletion, BIM/MCL-1/BCL-xL, caspases, LC3-dependent autophagy, drug-tolerant persisters, myxoid degeneration
6. **SDH-deficient pseudohypoxia axis** — complex II loss, succinate, PHD inhibition, HIF-1α, VEGF-A, IGF1R, global DNA hypermethylation, SDHB IHC loss
7. **Microenvironment, angiogenesis & immunity** — VEGF-A/VEGFR2, pericyte PDGFR-β, TIE2/ANGPT2, FGF2 and AXL/MET bypass, TAMs, IDO/Treg, CD8 TILs, PD-1/PD-L1, interstitial pressure
8. **Clonal evolution & secondary resistance** — V654A, T670I, D816H/V, D820Y, N822K, Y823D, polyclonality, intra-lesional heterogeneity, "nodule-in-a-mass", clonal sweep on withdrawal
9. **Imatinib PK/PD** — hOCT1 influx, ABCB1/ABCG2 efflux, CYP3A4 → CGP74588, AAG binding, the 1,100 ng/mL trough threshold, TDM, time-dependent clearance, DDI, adherence
10. **Sunitinib & regorafenib** — SU12662 and M-2/M-5 metabolites, the exon 13/14 vs exon 17/18 activity split, the anti-angiogenic component
11. **Ripretinib, avapritinib & genotype-directed agents** — switch-control dual-site binding, DP-5439, type I active-conformation binding, bezuclastinib, larotrectinib, dabrafenib+trametinib, cytoreductive surgery, adjuvant therapy
12. **Tumor burden, imaging & clinical endpoints** — RECIST, Choi density, FDG-PET SUVmax, clone-resolved ctDNA VAF, molecular-before-radiologic progression, AFIP risk, tumor rupture, GI bleeding
13. **Organ-system toxicity** — periorbital edema, myelosuppression, hepatotoxicity, hypothyroidism, hypertension, hand-foot skin reaction, alopecia, avapritinib cognitive effects, dose reduction → exposure loss → resistance
14. **Management decision nodes** — mandatory genotyping, neoadjuvant therapy, R0 resection, continue-vs-interrupt, line sequencing, the INTRIGUE decision, imatinib rechallenge, trial referral

## 4. mrgsolve model — 34 ODE compartments

| Group | Compartments |
|---|---|
| Imatinib PK | `IMA_DEPOT`, `IMA_CENT`, `IMA_PER` |
| Sunitinib / regorafenib / ripretinib / avapritinib PK | `SUN_*`, `REG_*`, `RIP_*`, `AVA_*` (2 each) |
| Protein binding | `AAG` (acute-phase, tumor-burden driven) |
| Tumor sub-populations | `TS`, `TQ`, `TA`, `TL`, `TD`, `TI` |
| Signal transduction | `PKIT`, `AKTP`, `ERKP`, `ETV1`, `BIM` |
| Vascular support | `ANG` |
| Imaging | `SUV`, `HU`, `MYX` |
| Molecular | `CTDS`, `CTDR` (clone-resolved ctDNA) |
| Toxicity | `EDEMA`, `ANC`, `TSH`, `SBP`, `HFS` |

Five genotypes (`GENO` 1-5) reconfigure the IC50 matrix, the baseline clonal
composition, the vascular dependence and the apoptotic competence. Thirteen
prebuilt scenarios are documented at the bottom of the model file.

## 5. What the model reproduces — and what it doesn't

These are **simulation outputs**, produced by numerically integrating the model
as written and compared against the cited trials. They are not targets.

| Observation | Model output | Source |
|---|---|---|
| First-line imatinib, KIT exon 11 | best response −56%, radiologic PD at **23.9 months** | B2222 / EORTC 62005 (~20-24 mo) |
| FDG-PET metabolic response | SUV **−63% by day 7** | Van den Abbeele 2008 (24-48 h) |
| Choi vs RECIST discordance | Choi PR at **day 25**, RECIST PR at **day 101** — a 76-day lead | Choi 2007 |
| Molecular before radiologic progression | resistant ctDNA detectable **day 531**, PD **day 727** (~6.5 mo lead) | Namløs 2018 |
| KIT exon 9 dose effect | 400 mg −22.9% vs 800 mg −54.4% best response | MetaGIST 2010 |
| Imatinib exposure–response | trough 373 → 560 → 746 ng/mL gives monotonically longer control | Demetri 2009 |
| Interruption at 1 year | −51% → −27% by day 545, above baseline by day 730, never rejoins the continuous curve | BFR14 |
| Adjuvant imatinib × 3 y | recurrence 1.8 y → 4.7 y ("delay, not cure") | SSGXVIII |
| PDGFRA D842V | PD at 2.9 mo on imatinib; **−86.7%** on avapritinib | NAVIGATOR |
| SDH-deficient | imatinib inert (PD 2.9 mo); sunitinib slows growth (PD 6.8 mo) **purely via the anti-angiogenic arm** | Janeway 2011 |
| BIM deletion polymorphism | same exposure and inhibition, best response −47.4% vs −56.0% | Ng 2012 |
| INTRIGUE ctDNA subgroups | exon 13/14-only → **sunitinib ahead**; exon 17/18-only → **ripretinib ahead by far more** | Nat Med 2024 |

**The one that failed.** The model does **not** reproduce the INTRIGUE primary
endpoint. In a mixed second-line population it puts ripretinib clearly ahead of
sunitinib, whereas the trial was negative. The reason is structural: with
sunitinib near-inert against activation-loop mutants, *any* population
containing an appreciable exon-17/18 clone favours ripretinib, because
sunitinib's failure there is catastrophic while ripretinib's disadvantage on
exon 13/14 is only moderate. For the arms to tie, one of three things must hold
and none is encoded here — ATP-pocket-only resistance is far commoner than the
~1:1 split assumed; sunitinib's anti-angiogenic contribution exceeds
`WVASC = 0.45`; or ripretinib's real exposure relative to its IC50s is lower
than `FT_RIP` implies. Each is a falsifiable quantitative hypothesis, and
distinguishing them is a better use of the model than tuning until the curves
overlap.

**A known artifact worth knowing.** Nadir-referenced RECIST progression
misbehaves when arms differ in *response depth*: a shallower responder has a
higher nadir, so the "≥20% above nadir" rule fires **later** even though the
patient is doing worse. This shows up in the exon-9 400 mg and BIM-deletion
arms. When comparing arms of unequal depth, use time-to-return-to-baseline
(helper provided in the model file) or read depth and progression separately.

## 6. Shiny dashboard — 8 tabs

1. **Patient profile** — genotype, burden, clearance, baseline clonal mix, and the live IC50 matrix in force for that patient
2. **Drug exposure (PK)** — all five agents, imatinib trough against the 1,100 ng/mL line, and the AAG/free-drug feedback
3. **Target engagement** — pathway activity against the cytostatic/cytocidal threshold; AKT, ERK, ETV1, BIM
4. **Clonal dynamics** — the six sub-populations on a log scale plus stacked composition over time
5. **Imaging endpoints** — RECIST, Choi density and SUV on one clock
6. **ctDNA & molecular lead time** — resistant-clone VAF against radiologic progression
7. **Scenario comparison** — eight prebuilt head-to-heads including the INTRIGUE subgroup question
8. **Toxicity** — edema, ANC, TSH, blood pressure, hand-foot

## 7. Quick start

```r
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread_cache("gist_mrgsolve_model.R")

# First-line imatinib 400 mg/d in KIT exon 11 metastatic disease, 4 years
out <- mod |>
  param(GENO = 1, TB0 = 500) |>
  mrgsim(ev(amt = 400, cmt = "IMA_DEPOT", ii = 1, addl = 1459),
         end = 1460, delta = 1) |>
  as_tibble()

ggplot(out, aes(time)) +
  geom_line(aes(y = RECIST_pct_change, colour = "RECIST diameter")) +
  geom_line(aes(y = Choi_density_pct,  colour = "Choi CT density")) +
  geom_line(aes(y = SUV_pct_change,    colour = "FDG-PET SUVmax")) +
  labs(x = "Day", y = "% change from baseline")

# The interactive dashboard
# shiny::runApp("gist_shiny_app.R")
```

Rendering the map:

```bash
dot -Tsvg gist_qsp_model.dot -o gist_qsp_model.svg
dot -Tpng -Gdpi=150 gist_qsp_model.dot -o gist_qsp_model.png
```

## 8. Limitations

- Tumor burden is a single well-mixed volume per clone. Real GIST resistance is
  **spatially** organised ("nodule-in-a-mass", different mutations in different
  metastases); the model captures the *composition* consequence of that but not
  its geometry, and so cannot represent focal resection of a resistant nodule
  except as an instantaneous reduction in one compartment.
- `FT_*` and the IC50 matrix are jointly calibrated; only their ratio is
  identifiable from clinical response data. Individual values must not be
  quoted as measured unbound fractions or biochemical potencies.
- Mutant-supply coefficients are lumped phenomenological parameters, not
  per-division point-mutation probabilities.
- No hepatic, renal or gastric-emptying covariates on PK; no explicit CYP3A4
  genotype, only a lumped clearance multiplier.
- Immune terms (IDO, CD8, PD-1) are present in the mechanistic map but are not
  dynamically coupled in the ODE model.
- Not fitted to individual patient data. Parameters are anchored to published
  trial-level summaries.

## ⚠️ 면책 조항 (Disclaimer)

교육 및 연구 목적의 정성적·반정량적 QSP 모델입니다. 공개 문헌과 임상시험
데이터를 바탕으로 구성되었으나 독립적으로 검증·인증되지 않았으며, **실제 임상
의사결정, 처방, 또는 규제 제출에 직접 사용해서는 안 됩니다.**

For research, education and hypothesis generation only. Not validated for
clinical decision-making, prescribing, or regulatory submission.
