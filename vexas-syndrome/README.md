# VEXAS Syndrome — QSP Model

> **VEXAS** = **V**acuoles, **E**1 enzyme, **X**-linked,
> **A**utoinflammatory, **S**omatic — first defined by Beck *et al.*
> (NEJM 2020). A myeloid-restricted somatic mutation in **UBA1** (codon
> Met41) abolishes the cytoplasmic UBA1b isoform, leading to defective
> cytoplasmic ubiquitination, ER stress, NF-κB/NLRP3/IFN
> hyperactivation, and a treatment-refractory hematoinflammatory
> phenotype with **macrocytic anemia, fevers, neutrophilic dermatosis,
> chondritis, pulmonary infiltrates, vasculitis and VTE**.

| File | Description |
|------|-------------|
| `vexas_qsp_model.dot` | Mechanistic map (≈110 nodes, 11 functional clusters) |
| `vexas_qsp_model.svg` | Vector render of the map |
| `vexas_qsp_model.png` | 150-dpi raster preview |
| `vexas_mrgsolve_model.R` | mrgsolve QSP model — 23 ODE compartments, 7 therapy scenarios |
| `vexas_shiny_app.R` | Shiny dashboard — 8 tabs |
| `vexas_references.md` | 70 PubMed-indexed references |

## Mechanistic snapshot

![VEXAS QSP mechanistic map](vexas_qsp_model.png)

**Eleven functional clusters in the map:**

1. **Genetics & HSC mosaicism** — UBA1a/UBA1b isoforms; Met41 alleles; VAF.
2. **Ubiquitin proteostasis** — E1/E2/E3 cascade; proteasome/autophagy clearance.
3. **ER stress & UPR** — PERK–eIF2α–ATF4–CHOP; IRE1α-XBP1s; ATF6.
4. **Mitochondrial stress & ROS** — mtDNA leak, cGAS–STING, oxidized cardiolipin.
5. **Innate inflammation** — NF-κB, NLRP3/AIM2, IL-1β/IL-6/TNF/IL-8/IFN-α.
6. **Bone marrow & cytopenias** — myeloid vacuolation, dyserythropoiesis, MDS overlap.
7. **Multisystem clinical phenotype** — fever, skin, chondritis, lung, vasculitis, VTE.
8. **Diagnostics** — BM vacuoles, UBA1 NGS/ddPCR, IFN score, severity score.
9. **Pharmacology** — prednisone, tocilizumab, anakinra, canakinumab, ruxolitinib, azacitidine, HSCT.
10. **Toxicity** — steroid GIO/DM/infection, IL-6 blockade diverticulitis, transfusion iron overload.
11. **Outcomes** — prednisone-eq dose, Hb, PLT, CRP, fever days, OS, VAF kinetics.

## mrgsolve model (23 ODE compartments)

| Compartment | Pathophysiology |
|-------------|-----------------|
| `VAF`              | UBA1 mutant fraction (0–1) — logistic clonal expansion |
| `MISF`, `ROS`, `ERST` | Proteostasis/oxidative/ER stress feedbacks |
| `IL1B`, `IL6`, `TNFa`, `IFNa`, `CXCL8`, `CCL2` | Core cytokines/chemokines |
| `CRP`, `FER`       | Acute-phase reactants (IL-6 / IL-1β driven) |
| `HB`, `PLT`, `ANC` | Inflammation-suppressed hematopoiesis |
| `FEV`, `SKIN`, `VTE` | Clinical activity indices |
| `HPA`              | Endogenous cortisol output (steroid suppression) |
| `AGUT_PRED, CC_PRED` | Prednisone PK |
| `ADEP_TOC, CC_TOC`   | Tocilizumab SC PK (with IL-6R Kd) |
| `ADEP_ANA, CC_ANA`   | Anakinra SC PK (IL-1R Kd) |
| `ADEP_CAN, CC_CAN`   | Canakinumab SC PK |
| `AGUT_RUX, CC_RUX`   | Ruxolitinib oral PK + JAK Kd |
| `ADEP_AZA, CC_AZA`   | Azacitidine SC PK |

**Seven therapy scenarios**: untreated · prednisone taper · tocilizumab + low-dose GC · anakinra · ruxolitinib + GC · azacitidine · allogeneic HSCT.

## Shiny app

`vexas_shiny_app.R` exposes 8 tabs:

1. **Patient profile** — editable demographics & baselines.
2. **UBA1 clone & VAF** — clonal expansion / shrinkage with HSCT or azacitidine.
3. **Cytokine storm** — IL-6 / IL-1β / TNF, IFN-α / CXCL8 / CCL2 with free-vs-bound resolution.
4. **Hematology** — Hb, PLT, ANC trajectories.
5. **Clinical activity** — fever, skin/chondritis, VTE risk indices.
6. **Drug PK** — prednisone, JAKi, biologics, and receptor occupancy.
7. **Scenario compare** — head-to-head efficacy plots.
8. **Biomarker panel** — CRP, ferritin, composite VEXAS activity score, week-12 endpoint table.

## Calibration anchors (see `vexas_references.md`)

* VAF: Beck 2020 (NEJM) — median ~75% monocyte VAF, M41T allele.
* IL-6 80–250 pg/mL in active disease; CRP 80–200 mg/L; ferritin 1–10 k ng/mL.
* Hb 8–10 g/dL macrocytic; PLT 80–150 ×10⁹/L; ANC variable.
* Tocilizumab SC: t½ 11–13 d (Frenzel 2022); Anakinra t½ 4–6 h; Ruxolitinib t½ 3 h; Azacitidine t½ ~0.7 h SC.
* Survival: ~50–60% at 5 y untreated; improved with HSCT (Hadjadj 2024).

## Build / run

```bash
# Render the mechanistic map
dot -Tsvg vexas_qsp_model.dot -o vexas_qsp_model.svg
dot -Tpng -Gdpi=150 vexas_qsp_model.dot -o vexas_qsp_model.png

# Run the mrgsolve model
Rscript -e 'source("vexas_mrgsolve_model.R")'

# Launch the Shiny app
Rscript -e 'shiny::runApp("vexas_shiny_app.R", launch.browser = TRUE)'
```

*Model status: research prototype; parameters informed by 2020-2024 cohort
publications and biologic PK labels. Not validated for clinical decision making.*
