# Sarcopenia — Age-Related Skeletal Muscle Wasting (QSP)

> Quantitative Systems Pharmacology model of primary (age-related) sarcopenia
> and its principal pharmacological/non-pharmacological interventions.

![Mechanistic map](sarc_qsp_model.png)
[🗺️ Full SVG map](sarc_qsp_model.svg) ·
[⚙️ mrgsolve model](sarc_mrgsolve_model.R) ·
[📚 References](sarc_references.md) ·
[📊 Shiny app](sarc_shiny_app.R)

---

## 1. Disease overview

Sarcopenia is a progressive, generalized skeletal-muscle disorder characterized
by loss of muscle mass *and* strength/function, leading to falls, fractures,
hospitalization and mortality (EWGSOP2 / AWGS2019 / SDOC / FNIH). Prevalence
rises from ~10 % at age 60 to >30 % at age ≥80; sex-specific cutoffs are
applied for grip strength (M < 27, F < 16 kg) and ALMI
(M < 7.0, F < 6.0 kg/m²), with severity flagged by gait speed < 0.8 m/s.

### Primary mechanisms encoded in this model
1. **Aging drivers** — chronological aging, cellular senescence (SASP),
   inflammaging (IL-6, TNFα, CRP), mitochondrial dysfunction, hormonal
   decline (T↓, IGF-1↓, GH↓, DHEA↓), motor-neuron loss, malnutrition,
   physical inactivity, oxidative stress, gut dysbiosis, epigenetic drift.
2. **Anabolic axis** — leucine/EAA → Sestrin2/Ragulator → mTORC1
   (S6K1, 4E-BP1) → MPS; IGF-1/insulin → IR/IGF-1R → IRS-1 → PI3K → Akt;
   Akt ⊣ TSC1/2 → Rheb-GTP → mTORC1.
3. **Catabolic axis** — Akt ⊣ FoxO; FoxO → Atrogin-1 (MAFbx) + MuRF1 →
   K48-Ub → 26S proteasomal proteolysis; FoxO → LC3, Beclin-1 →
   macroautophagy; calpain/caspase-3 cleavage; AMPK → ULK1 → autophagy.
4. **Myostatin–Activin–ActRII** — promyostatin → mature MSTN; MSTN /
   Activin A/B / GDF-11 → ActRIIB(>A) → ALK4/5 → SMAD2/3-SMAD4 →
   atrogene transcription, anabolic resistance; follistatin as endogenous
   antagonist.
5. **Inflammaging & cytokines** — IL-6 → JAK/STAT3 → SOCS3 ⊣ IRS-1
   (anabolic resistance); TNFα → IκBα phospho → NF-κB → atrogenes;
   IL-15 anti-sarcopenic; macrophage M1, COX-2/PGE2.
6. **Mitochondrial energetics** — PGC-1α → NRF1/2 → TFAM → mtDNA;
   OXPHOS / ATP; PINK1/Parkin mitophagy; SIRT1↔NAD⁺; UCP3, Drp1/Mfn.
7. **Satellite-cell pool** — Pax7⁺ → MyoD → Myf5 → myogenin → myotube;
   senescent (p16⁺) niche, Notch/Wnt, GDF-11-mediated aging.
8. **NMJ / fiber-type** — α-MN loss, agrin/MuSK/LRP4, denervation →
   preferential Type IIX atrophy, reinnervation → fiber-type grouping,
   myosteatosis.

---

## 2. Mechanistic map

* **12 clusters · 130+ nodes · 220+ edges** (DOT source).
* SVG and PNG (150 dpi) provided.
* Solid arrows = activation; T-arrows = inhibition; dotted = clinical gates.

## 3. mrgsolve model

* **21 ODE compartments**:
  Bimagrumab (2-cpt PK), Apitegromab (2-cpt PK), Testosterone (oil depot + 1-cpt),
  Vitamin D3 + 25(OH)D pool, Leucine (gut + central), Anamorelin (gut + central),
  dynamic IGF-1, total myostatin, IL-6, GDF-15, ALM, grip, gait, SPPB,
  cumulative falls, cumulative frailty index.
* **Anabolic drive** — multiplicative Emax for leucine, IGF-1, testosterone,
  bimagrumab, apitegromab, exercise, vitamin D; divided by anabolic-resistance
  + SOCS3-IL-6 attenuator.
* **Catabolic drive** — Emax for free myostatin, IL-6/TNFα, disuse;
  attenuated by ActRII-mAb sequestration (bimagrumab).
* **Clinical translation** — grip / gait / SPPB driven by ALM with vitamin-D
  and testosterone modifiers; fall hazard ∝ (gait_threshold / gait)².
* **Calibration anchors**: Janssen 2000 (−0.46 %/yr ALM), BMR-101 bimagrumab
  (+7.1 % ALM @ 24 wk), Bhasin 1996 (testosterone +6.1 % LBM @ 10 wk),
  Bischoff-Ferrari 2009 (Vitamin D fall RR 0.81).

### Six treatment scenarios (driver in script):
| # | Arm |
|---|-----|
| 1 | No treatment |
| 2 | Bimagrumab 30 mg/kg IV q4w × 1 yr |
| 3 | Apitegromab 20 mg/kg IV q4w × 1 yr |
| 4 | Testosterone 100 mg/wk IM × 1 yr |
| 5 | Vit D 2000 IU/d + EAA + RT 3×/wk |
| 6 | Combo (Bimagrumab + RT 3×/wk) |

## 4. Shiny dashboard (8 tabs)

1. Patient profile & EWGSOP2 status (sex/age cutoffs)
2. Drug PK (Bimagrumab, Apitegromab, Testosterone, Vit D, EAA)
3. Anabolic vs catabolic drive, endocrine/cytokine trajectories
4. Muscle mass & function (ALM, grip, gait, SPPB)
5. Clinical endpoints (ALMI vs cutoff, falls, frailty, dx table)
6. Six-arm scenario comparison + ALM/gait at 1 yr
7. Biomarker panel (IGF-1, MSTN, IL-6, GDF-15, 25(OH)D)
8. Curated references (rendered Markdown)

## 5. References

See [`sarc_references.md`](sarc_references.md) — 72 curated PubMed entries
across 18 thematic groups (definition, aging, mTORC1, UPS, myostatin,
mitochondria, satellite cells, inflammaging, anabolic resistance, drug
trials, exercise, NMJ, sarcopenic obesity, biomarkers, outcomes,
methodology).

## 6. Limitations & caveats

* Aging covariates are **lumped** — model does not resolve individual
  sarcomere proteins or single-fiber dynamics.
* Bimagrumab failed primary endpoints in BMR-101 community-dwelling
  elderly (Rooks 2020); secondary outcomes (lean mass, walking) improved.
  EC50 / Emax values used are illustrative, not regulatory anchors.
* SARMs and testosterone carry cardiovascular and prostate-safety signals
  not modeled here.
* Anabolic resistance and inflammaging are coupled but with simplified
  feedback; SOCS3 module is empirical.

---

*Generated as part of the daily QSP library (CLAUDE.md routine).*
