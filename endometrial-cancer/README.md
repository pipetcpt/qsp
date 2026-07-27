# Endometrial Carcinoma (EC) — QSP Model

> Integrated Quantitative Systems Pharmacology model of endometrial carcinoma,
> built around a single organising idea: the **TCGA molecular classes are not
> four diseases but four points on one continuous axis** — tumour mutational
> burden — and essentially every therapeutic difference between them can be
> derived from where they sit on that axis rather than fitted separately.
> Couples the unopposed-oestrogen endocrine engine (adipose aromatase → ERα,
> opposed by PR-B) and the PI3K/PTEN oncogenic engine to a TMB → neoantigen →
> T-cell-priming cascade, with mechanistic PK/PD for carboplatin, paclitaxel,
> anti-PD-1 (dostarlimab / pembrolizumab), lenvatinib, progestins (megestrol
> and the levonorgestrel IUS), letrozole, everolimus and trastuzumab.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`emc_qsp_model.dot`](emc_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`emc_qsp_model.svg`](emc_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`emc_qsp_model.png`](emc_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`emc_mrgsolve_model.R`](emc_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`emc_shiny_app.R`](emc_shiny_app.R) |
| 📚 References (75)        | [`emc_references.md`](emc_references.md) |

---

## 1. Disease in one paragraph

Endometrial carcinoma is the most common gynaecologic malignancy in
high-income countries and one of the very few cancers whose incidence and
mortality are both **rising**, tracking the obesity epidemic. Two engines
drive it. The **endocrine engine**: adipose tissue aromatase (CYP19A1)
converts adrenal and ovarian androstenedione to estrone, hyperinsulinaemia
suppresses SHBG and so raises the *free* estradiol fraction, and in the
absence of luteal progesterone — anovulation, PCOS, nulliparity, the
postmenopausal state — the endometrium experiences **unopposed oestrogen**,
driving ERα transcription, mitosis, hyperplasia, endometrial intraepithelial
neoplasia and finally endometrioid carcinoma. Key and Pike's classic
dose-response argument is the quantitative form of this: endometrial mitotic
rate scales with unopposed oestrogen exposure, and cancer risk scales with
accumulated mitoses. The **genomic engine**: PI3K pathway lesions
(PTEN ~77%, PIK3CA ~53%, PIK3R1 ~33%) are near-universal in endometrioid
tumours, joined by ARID1A, CTNNB1, KRAS and FGFR2; and superimposed on all of
it is DNA-repair status — POLE exonuclease-domain mutation, mismatch-repair
deficiency (Lynch germline or sporadic MLH1 promoter hypermethylation), or
neither — with TP53 abnormality marking the copy-number-high, serous-like,
immune-excluded class that carries most of the mortality.

## 2. Why this model is built the way it is

Endometrial cancer is the solid tumour where **molecular class, not histology
and not stage, reorders the entire therapeutic decision tree**. Most QSP
treatments of a four-class disease would build four parameterisations. This
one does not. Instead it commits to a single equation:

```
IMMUNOGENICITY = TMB^h / (TMB50^h + TMB^h),   TMB50 = 13 mut/Mb,  h = 1.4
```

parameterised from the **mutation-rate literature** (TCGA, Church 2013,
Howitt 2015, Chan 2019) and then *checked against* the immunotherapy trials
rather than fitted to them. Feeding in the four class-default burdens gives:

| Class | TMB (mut/Mb) | Immunogenicity | Prevalence | 5-y PFS |
|---|---|---|---|---|
| POLEmut (ultramutated) | 140 | **0.97** | 7-12% | ~95-100% |
| MMRd / MSI-H | 25 (×1.35 indel uplift) | **0.79** | 25-30% | ~75% |
| NSMP / copy-number-low | 2.9 | **0.11** | 30-40% | ~80% |
| p53abn / copy-number-high | 2.3 | **0.08** | 15-25% | ~50% |

Four consequences follow, and they are the reason to build it this way.

**(a) The biomarker works because the cut sits on the cliff.** TMB50 = 13
mut/Mb lands almost exactly where the clinical dMMR/pMMR dichotomy is drawn.
A binary test placed on the steepest part of a smooth sigmoid will look
extraordinarily predictive — RUBY reports HR 0.28 in dMMR against ~0.76 in
pMMR — while also guaranteeing that the patients it gets *wrong* are precisely
those near the boundary: MSI-low, MMR-heterogeneous, subclonal MLH1
methylation. The dashboard's **TMB continuum** tab holds the drug, dose,
patient and exclusion penalty fixed and moves only the mutational burden;
the result is a smooth curve, not a step.

**(b) Because the curve saturates, POLEmut is an argument for *less*
therapy.** At TMB 140 the immunogenicity is 0.97 — barely distinguishable
from a hypothetical TMB-60 tumour. POLEmut's spectacular prognosis is
therefore not a dose-response case for intensification; it is the
quantitative case for **de-escalation** (the RAINBO POLE-BLUE hypothesis).
In the adjuvant/micrometastatic setting (`TUM0 = 8`, three-year horizon) the
model puts POLEmut at **4.8 mm under observation versus 4.6 mm with adjuvant
chemotherapy** — the chemotherapy makes the fall faster, not deeper. The same
comparison in p53abn moves time-to-progression from **18 days to 146 days**,
mirroring the PORTEC-3 molecular analysis.

**(c) At the floor of the curve, the pharmacology has to differ in kind, not
in dose.** pMMR tumours cannot be rescued by more checkpoint blockade. They
are rescued by changing the microenvironment, and so **lenvatinib in this
model has no tumour-kill term at all**. Its contribution runs through
`LEN_MDSC` (myeloid-suppressor depletion), `LEN_NORM` (vascular normalisation
improving T-cell delivery, peaking at *intermediate* VEGFR occupancy) and the
plain perfusion-starvation effect of losing microvessel density. Knocking out
the first two in an NSMP tumour takes the pembrolizumab+lenvatinib response
from **−60% to −26%** and the survival surrogate from **0.64 to 0.43** — so
about two thirds of the combination benefit is immunological reprogramming
and one third is perfusion starvation. That split is the model's falsifiable
claim about the KEYNOTE-775 mechanism, and it is checkable in one line:

```r
mod |> param(MOLCLASS = 3, LEN_MDSC = 0, LEN_NORM = 0) |> mrgsim(pembro + lenva, end = 730)
```

**(d) A third class axis runs the opposite way from the first two.**
Immunogenicity and chemosensitivity are *not* the same variable, and
collapsing them would make the model unable to say which one drives the
prognostic gap. Losing the p53 G1 checkpoint lets damaged cells enter S and M
phase, so the class that responds worst to checkpoint blockade (p53abn,
`CHS_P53 = 1.55`) responds best to platinum and taxane. Aggressiveness is a
third, separate axis again (`AGG_P53 = 1.35`). p53abn is therefore the most
aggressive, the least immunogenic, and the most chemosensitive class
simultaneously — which is exactly why it needs its own treatment strategy
rather than more of anyone else's.

A fifth, smaller commitment worth flagging: **platinum immunogenic cell
death** (`ICD_CHEMO`) raises neoantigen release during chemotherapy. This is
what lets the model reproduce the otherwise awkward pair of facts that
dostarlimab **monotherapy** does almost nothing in pMMR (GARNET, ORR 15.4%)
while dostarlimab **plus chemotherapy** does something real (GY018 pMMR
HR 0.54). Same drug, same tumour, different partner.

## 3. Mechanistic map — 14 clusters, ~175 nodes

1. **Host risk factors / metabolic-endocrine milieu** — obesity, adipose mass, insulin resistance, free IGF-1, SHBG suppression, adipokine shift, PCOS, nulliparity, tamoxifen exposure, unopposed ERT, Lynch and Cowden syndromes
2. **Steroid hormone biosynthesis & HPO axis** — GnRH/LH/FSH, ovarian and adrenal steroidogenesis, aromatase CYP19A1, estrone → 17β-HSD1 → estradiol, luteal progesterone, CYP1B1 catechol-oestrogen genotoxicity
3. **Endometrial ER / PR nuclear receptor signalling** — ERα/ERβ/GPER1, coactivators, ERE transcription, PR-A/PR-B isoforms, PR down-regulation, HAND2, FOXO1 decidualisation, hyperplasia → EIN
4. **Somatic driver landscape & mutational processes** — PTEN, PIK3CA, PIK3R1, ARID1A, CTNNB1, KRAS, FGFR2, TP53, PPP2R1A, ERBB2, CCNE1, MYC, POLE exonuclease domain, MMR genes, MLH1 hypermethylation, MSI
5. **TCGA / ProMisE molecular classification** — the four classes, the surrogate classifier, histology, grade, ER/PR IHC, L1CAM, FIGO 2023 molecularly integrated staging
6. **Growth-factor signal transduction** — PI3K/PIP3/AKT/mTORC1/mTORC2, S6K1-IRS1 feedback, FOXO3a, RAS-RAF-MEK-ERK, Wnt/β-catenin, HER2, FGFR2, and ligand-independent ERα phosphorylation
7. **Cell cycle, apoptosis, genomic instability** — cyclin D1-CDK4/6, RB, E2F, cyclin E-CDK2, p53 checkpoint, replication stress, aneuploidy, subclonal heterogeneity
8. **TMB → neoantigen → antigen presentation** — the continuous axis, frameshift versus SNV neoantigens, proteasome/TAP, MHC-I, B2M and JAK1/2 escape, dendritic-cell cross-presentation
9. **Adaptive anti-tumour immunity & checkpoints** — CD8 priming, IFN-γ, granzyme/perforin, PD-1/PD-L1 including IFN-γ-driven adaptive resistance, CTLA-4, LAG-3/TIM-3, exhaustion, memory, tertiary lymphoid structures
10. **Immunosuppressive microenvironment** — Treg, MDSC, M2 TAM, TGF-β, IL-10, IDO1, CD39/CD73-adenosine, CAF, hypoxia/HIF-1α, immune exclusion, obesity-driven T-cell dysfunction
11. **Angiogenesis & vascular remodelling** — VEGF-A/VEGFR-2, FGF escape, PDGFR-β, Ang2-Tie2, microvessel density, the vascular-normalisation window, interstitial pressure
12. **Invasion, LVSI & metastasis** — EMT, MMP-2/9, deep myometrial invasion, substantial LVSI, sentinel and para-aortic nodes, peritoneal and distant spread
13. **Drug PK/PD** — carboplatin, paclitaxel, dostarlimab, pembrolizumab, lenvatinib, megestrol, LNG-IUS, letrozole, everolimus, trastuzumab, selinexor, metformin, tamoxifen alternation, weight loss
14. **Clinical endpoints, biomarkers & toxicity** — ORR/PFS/OS, CA-125, ctDNA, HE4, abnormal uterine bleeding, pathological CR and fertility preservation, neutropenia, neuropathy, hypertension, proteinuria, irAEs, VTE

![EC mechanistic map](emc_qsp_model.png)

## 4. mrgsolve model — 35 ODE compartments

| Block | Compartments |
|---|---|
| PK (13) | `CARB_C`; `PTX_C`/`PTX_P`; `IO_C`/`IO_P`; `LEN_G`/`LEN_C`; `PRG_G`/`PRG_C`; `LTZ_G`/`LTZ_C`; `EVE_G`/`EVE_C` |
| Endocrine | `E2` (free estradiol), `ER_ACT` (ERα transcriptional activity), `PR_EXP` (PR-B expression) |
| Oncogenic signalling | `PI3K_A` |
| Tumour | `TUM_S` (sensitive), `TUM_R` (resistant clone) |
| Immunity | `NEO`, `TEFF`, `TEXH`, `TREG`, `MDSC` |
| Vasculature | `VEGF`, `MVD` |
| Biomarkers | `CA125`, `CTDNA` |
| Myelosuppression | `PROL`, `TR1`, `TR2`, `TR3`, `CIRC` (Friberg) |
| Toxicity / outcome | `SBP`, `HAZ` |

Three mechanistic details worth pointing at, because they are where the model
earns its keep rather than merely bookkeeping:

- **PR-B is a state, not a parameter.** Sustained progestin occupancy
  down-regulates the very receptor the drug needs (`K_PRDOWN`), which is the
  mechanism of the plateau and escape seen at 6–12 months of continuous
  progestin in fertility-sparing therapy. Tamoxifen alternation re-induces it
  (`K_TAM_PR`, the GOG-153 rationale). The levonorgestrel IUS is modelled as a
  route-dependent local exposure multiplier (`IUS_FACTOR`), which is the whole
  point of the device: high endometrial concentration from a low systemic one.
- **ERα has two inputs.** Ligand (estradiol) *and* ligand-independent
  activation by mTOR/ERK phosphorylation of Ser167/118 (`W_LIGIND`). This is
  why letrozole alone underperforms in a PI3K-driven tumour and why everolimus
  restores its activity — no added efficacy parameter is involved. Worth being
  precise about what this produces: in an NSMP tumour neither agent drives
  frank regression, but two-year burden falls from **227 mm** (letrozole
  alone) to **118 mm** (the combination), with everolimus alone at 139 mm.
  That is stabilisation with markedly slowed growth rather than response —
  which, for a regimen whose real ORR is 32% (two thirds of patients do *not*
  respond), is arguably the right answer for a deterministic typical patient.
- **Trastuzumab acts twice, and only in HER2-amplified tumours.** It removes
  HER2 pathway drive *and* recruits NK-cell ADCC (`K_HER2KILL`). Both terms
  are gated on `HER2_AMP`, so setting `HER2_TX = 1` in a HER2-normal tumour
  changes the trajectory by exactly nothing (−49.9% either way) while in the
  amplified serous tumour it deepens response from −44% to −91%. That
  on/off behaviour is the entire content of Fader 2018.
- **Lenvatinib's hypertension and its efficacy come from the same term.**
  `K_SBP_LEN` is driven by the identical VEGFR-2 occupancy that produces
  vascular normalisation, so the model structurally cannot separate benefit
  from toxicity without separating VEGFR-2 from the other kinases — which is
  an honest representation of why the 20 → 14 → 10 mg reduction ladder costs
  efficacy.

### Eight prebuilt scenarios

| # | Scenario | Anchor |
|---|---|---|
| 1 | Untreated natural history (any class) | — |
| 2 | Carboplatin–paclitaxel × 6 | GOG-0209 (Miller 2020) |
| 3 | Dostarlimab + chemo, dMMR/MSI-H | RUBY (Mirza 2023), HR 0.28 |
| 4 | Dostarlimab/pembrolizumab + chemo, pMMR | NRG-GY018 (Eskander 2023), HR 0.54 |
| 5 | Pembrolizumab + lenvatinib, pMMR | KEYNOTE-775 (Makker 2022), HR 0.60 |
| 6 | Fertility-sparing LNG-IUS, grade 1 ER+/PR+ | feMMe (Janda 2021), Westin 2021 |
| 7 | Letrozole + everolimus | Slomovitz 2015, ORR 32% |
| 8 | Megestrol ± tamoxifen alternation | GOG-119 / GOG-153 |

Plus a ninth, non-therapeutic experiment — `EMC_tmb_sweep()` — which holds the
regimen fixed and sweeps TMB across two orders of magnitude. That one is the
model's actual thesis.

```r
library(mrgsolve)
mod <- mread_cache("emc_mrgsolve_model.R")

# Scenario 3 — RUBY-like dMMR arm
chemo <- ev(amt = 625, cmt = "CARB_C", ii = 21, addl = 5) +
         ev(amt = 300, cmt = "PTX_C",  ii = 21, addl = 5)
dost  <- ev(amt = 500,  cmt = "IO_C", ii = 21, addl = 3) +
         ev(amt = 1000, cmt = "IO_C", time = 84, ii = 42, addl = 15)

mod |> param(MOLCLASS = 2) |> mrgsim(chemo + dost, end = 730, delta = 1) |> plot()

# The same regimen in a pMMR tumour — change one number
mod |> param(MOLCLASS = 3) |> mrgsim(chemo + dost, end = 730, delta = 1) |> plot()
```

### Verified behaviour against the calibration targets

Every number below was produced by running the committed model file, not
asserted. Advanced-disease rows use `TUM0 = 60`, BMI 34, grade 2, 730-day
horizon; adjuvant rows use `TUM0 = 8`, 1095-day horizon. "Best" is best
percentage change from baseline; "S" is the survival surrogate at the horizon.

| Target | Trial expectation | Model output | Met? |
|---|---|---|---|
| Untreated class ordering | POLE < MMRd < NSMP < p53abn burden | 115 / 158 / 240 / 310 mm | ✅ |
| IO monotherapy, dMMR vs pMMR | GARNET ORR 45.5% vs 15.4% | best −61% vs 0% | ✅ |
| IO + chemo, dMMR | RUBY HR 0.28 | S 0.82 vs 0.09 chemo alone | ✅ |
| IO + chemo, pMMR | GY018 HR 0.54 (real but smaller) | best −49% vs −33%; S 0.028 vs 0.002 | ✅ |
| Lenvatinib + pembro, pMMR | KEYNOTE-775 HR 0.60 | best −60%, S 0.64 vs 0.012 pembro alone | ✅ |
| Lenvatinib single-agent | weak (ORR ~14%) | no regression; growth slowed only | ✅ |
| POLEmut adjuvant de-escalation | PORTEC-3: excellent regardless | 4.8 mm obs vs 4.6 mm chemo | ✅ |
| p53abn adjuvant chemo benefit | PORTEC-3: concentrated here | TTP 18 → 146 d, best 0% → −62% | ✅ |
| Trastuzumab, HER2+ vs HER2-normal | Fader 2018: benefit only if amplified | −91% vs −44%; **identical** if HER2-normal | ✅ |
| Progestin resistance clock | plateau/escape at 6–12 months | PR-B 1.00 → 0.61 → 0.59 over 360 d | ✅ |
| Tamoxifen alternation | GOG-153: re-induces PR | PR-B 0.50 → 0.72 at 360 d; best −42% → −55% | ✅ |
| LNG-IUS fertility-sparing | feMMe pCR ~57–70% | best −45%, S 0.71 | ✅ |
| B2M loss uncouples TMB | escape despite high TMB | MMRd best −90% → −59%, S 0.82 → 0.07 | ✅ |
| MMRd adjuvant chemo benefit | PORTEC-3: **not** significant | TTP 238 d → not reached | ❌ over-predicted |

The last row is a genuine miss and is discussed in section 7.

## 5. Shiny dashboard — 9 tabs

1. **Patient & molecular profile** — the TMB → immunogenicity curve with all four classes plotted on it and the current patient ringed
2. **Drug PK** — chemotherapy, anti-PD-1 with receptor occupancy, lenvatinib with VEGFR-2 inhibition, endocrine agents
3. **Endocrine & PI3K PD** — free estradiol, ERα activity, the PR-B resistance clock, PI3K/AKT/mTOR
4. **Tumour dynamics** — RECIST sum against the ±20%/−30% thresholds, sensitive/resistant clone decomposition, angiogenesis
5. **Immune compartment** — neoantigen load, effector versus exhausted CD8, Treg/MDSC, PD-1 blockade
6. **Clinical endpoints** — survival surrogate, CA-125 and ctDNA, Friberg neutropenia, lenvatinib blood pressure
7. **Scenario comparison** — overlay any subset of the eight regimens, or hold the regimen fixed and overlay the four molecular classes
8. **TMB continuum** — the sweep described above, run live
9. **References**

```r
shiny::runApp("emc_shiny_app.R")   # requires emc_mrgsolve_model.R alongside it
```

## 6. Rendering the map

```bash
dot -Tsvg emc_qsp_model.dot -o emc_qsp_model.svg
dot -Tpng -Gdpi=150 emc_qsp_model.dot -o emc_qsp_model.png
```

## 7. What this model deliberately does *not* do

- It does not model **surgery, sentinel node mapping, or radiotherapy** as
  interventions, though the anatomic spread pathway that motivates them is in
  the map. The pharmacology was the target.
- It has **no spatial structure**. Immune exclusion — the defining feature of
  the p53abn microenvironment — is represented as a scalar penalty (`EXCL_P53`)
  rather than as a geometry. That is a real simplification: a model of
  stromal T-cell exclusion that cannot represent distance is describing the
  consequence, not the cause.
- The **survival output is a hazard surrogate**, not a fitted time-to-event
  model. It is monotone in tumour burden by construction, so it can compare
  arms but should not be read as a calibrated survival prediction.
- **Class prevalences are not simulated.** Each run is one patient of one
  declared class; there is no virtual-population sampling over the 7/28/35/20
  split, and no inter-individual variability at all. Adding that is the
  obvious next step and would let the model produce trial-level rather than
  patient-level predictions — an ORR instead of a single trajectory.
- **It over-predicts adjuvant chemotherapy benefit in MMRd.** Stated here
  rather than buried: in the adjuvant setting the model moves MMRd from
  time-to-progression 238 days under observation to not-reached with
  chemotherapy, whereas PORTEC-3 found no significant chemotherapy benefit in
  that subgroup. The cause is structural — debulking hands the residual
  tumour to an immune compartment that is competent in MMRd, so chemotherapy
  and immunity compound. Reproducing the trial would require
  chemotherapy-induced lymphodepletion, which this model does not have. Every
  other calibration target in section 4 is met; this one is not.

## 8. Disclaimer

This is a qualitative-to-semi-quantitative QSP model built for education,
research and hypothesis generation from public literature and published trial
results. It has not been independently validated or qualified. Parameters are
illustrative approximations; fitting and verification against real patient
data would be required for any applied use. **It must not be used for
clinical decision-making, prescribing, or regulatory submission.**
