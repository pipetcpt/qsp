# Phenylketonuria (페닐케톤뇨증, PKU) — QSP Model

> Integrated Quantitative Systems Pharmacology model of phenylketonuria,
> linking dietary phenylalanine (Phe) intake to saturable hepatic
> phenylalanine-hydroxylase (PAH)/tetrahydrobiopterin (BH4) enzymatic
> clearance, saturable blood-brain-barrier LAT1 large-neutral-amino-acid
> transport competition (Phe vs. Tyr/Trp), downstream cerebral dopamine/
> serotonin synthesis, myelination, and critical-period-weighted cumulative
> neurocognitive injury — across the classic/moderate/mild PKU and mild
> hyperphenylalaninemia (MHP) genotype-severity spectrum. The model spans
> dietary/medical-formula management and the modern pharmacology stack
> (sapropterin dihydrochloride BH4-cofactor therapy, pegvaliase PEG-PAL
> enzyme-substitution therapy including anti-drug-antibody immunogenicity
> and immune tolerization) plus maternal-PKU fetal-risk reporting.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`pku_qsp_model.dot`](pku_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`pku_qsp_model.svg`](pku_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`pku_qsp_model.png`](pku_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`pku_mrgsolve_model.R`](pku_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`pku_shiny_app.R`](pku_shiny_app.R) |
| 📚 References             | [`pku_references.md`](pku_references.md) |

---

## 1. Disease in one paragraph

Phenylketonuria is an autosomal-recessive inborn error of amino-acid
metabolism caused by loss-of-function mutations in the *PAH* gene
(12q23.2), which encodes hepatic phenylalanine hydroxylase — the enzyme
that, with the tetrahydrobiopterin (BH4) cofactor, converts phenylalanine
(Phe) to tyrosine (Tyr). Residual PAH activity defines a phenotypic
spectrum from classic PKU (<1% activity) through moderate and mild PKU to
benign mild hyperphenylalaninemia (MHP, >10-35% activity), detected by
mandatory newborn dried-blood-spot screening (tandem MS/MS). Untreated,
plasma Phe accumulates to levels that saturate and, via the shared LAT1
(SLC7A5) transporter, competitively exclude tyrosine and tryptophan from
crossing the blood-brain barrier, starving the brain of the precursors for
dopamine, norepinephrine and serotonin synthesis and impairing myelination
and protein synthesis — producing, if untreated in infancy, severe and
irreversible intellectual disability, seizures and microcephaly (the
namesake musty odor and eczema stem from the minor phenylpyruvate/
phenylacetate transamination pathway). Lifelong dietary phenylalanine
restriction with a Phe-free medical amino-acid formula (target plasma Phe
120-360 umol/L, ACMG/NIH guideline) begun at birth largely prevents this
injury, though even well-controlled patients show subtle, prefrontal-
dopamine-dependent executive-function vulnerability (Diamond 1997) because
the developing/adult brain remains sensitive to LAT1 competition at any
elevated Phe. Two disease-modifying drugs supplement diet: **sapropterin
dihydrochloride** (Kuvan), an oral synthetic BH4 cofactor/pharmacological
chaperone effective only in BH4-responsive (typically missense/misfolding)
genotypes, and **pegvaliase** (Palynziq), a subcutaneously injected
PEGylated phenylalanine-ammonia-lyase enzyme that metabolizes Phe via an
alternative, PAH-independent route (to trans-cinnamic acid), enabling
substantial diet liberalization in severe classic PKU but complicated by
anti-drug-antibody immunogenicity that typically attenuates early efficacy
before immune tolerization develops over months. Maternal PKU with
periconceptional Phe >900-1200 umol/L risks fetal congenital heart disease
and, with sustained third-trimester elevation, microcephaly and
neurocognitive deficit (maternal PKU syndrome), making strict
pre-conception dietary control essential.

## 2. Mechanistic clusters (16 in the DOT map, 122 nodes)

1. Genetics, inheritance & phenotype spectrum (*PAH* gene, genotype classes, classic/moderate/mild PKU, MHP, BH4-responsive vs. non-responsive alleles)
2. Newborn screening & diagnosis (dried-blood-spot MS/MS, confirmatory plasma amino acids, BH4 loading test, molecular testing, cofactor-deficiency screen)
3. Dietary Phe intake & GI absorption (natural protein, Phe-free medical formula, glycomacropeptide food, tyrosine supplementation)
4. Hepatic PAH enzyme & BH4 catalytic cycle (PAH misfolding/chaperone rescue, cofactor recycling via PCD/DHPR)
5. De novo BH4 biosynthesis (GTPCH, PTPS, sepiapterin reductase)
6. Alternative minor Phe catabolism/transamination (phenylpyruvate, phenylacetate "mousy odor", urinary phenylketones)
7. Plasma amino-acid pool & BBB LAT1 competition (Phe/Tyr/Trp/BCAA, competitive inhibition kinetics)
8. Cerebral neurotoxicity & monoamine synthesis (tyrosine/tryptophan hydroxylase, dopamine/serotonin/norepinephrine, myelin, oxidative stress)
9. Peripheral tyrosine-deficiency downstream effects (hypopigmentation, catecholamine deficiency, growth)
10. Maternal PKU & fetal teratogenicity (placental Phe transport, congenital heart disease, microcephaly)
11. Clinical & neurocognitive endpoints (IQ trajectory, executive function, ADHD-like symptoms, seizures, QoL)
12. Drug PK/PD — sapropterin dihydrochloride (BH4 cofactor/chaperone)
13. Drug PK/PD — pegvaliase (PEG-PAL enzyme substitution, ADA immunogenicity)
14. Investigational & emerging therapies (SYNB1618 gut probiotic, mRNA-PAH, AAV gene therapy, CRISPR base-editing)
15. Non-pharmacologic management & monitoring (lifelong diet, home DBS monitoring, multidisciplinary clinic)
16. Adverse effects & safety monitoring (drug AEs, anaphylaxis REMS, micronutrient deficiency, bone health)

## 3. mrgsolve model (18 ODE compartments)

* **Drug PK (2 agents + immunogenicity, 6 compartments)** — sapropterin
  (oral gut depot + plasma), pegvaliase (SC depot + plasma), anti-drug-
  antibody titer, immune-tolerization fraction.
* **Amino-acid mass-balance network (6 compartments)** — plasma Phe/Tyr/Trp
  via Michaelis-Menten hepatic PAH clearance + pegvaliase alternative
  clearance + linear minor-pathway/renal clearance; brain (BBB) Phe/Tyr/Trp
  via saturable, competitively-inhibited LAT1 transport.
* **Downstream PD & clinical endpoints (6 compartments)** — dopamine- and
  serotonin-synthesis capacity indices, myelin/white-matter integrity,
  cumulative IQ deficit (critical-developmental-period-weighted, largely
  irreversible), executive-function deficit (reversible, tracks current
  control), growth Z-score deficit.
* Genotype residual PAH activity (classic <1%, moderate 1-5%, mild 5-10%,
  MHP >10-35%) and BH4-responsiveness are patient-level covariates;
  Michaelis-Menten and BBB-transport parameters are calibrated so that
  simulated steady-state plasma Phe reproduces clinically reported ranges
  for untreated classic PKU (~1500 umol/L), diet-controlled classic PKU
  (~330 umol/L, within the ACMG band), untreated MHP (~200 umol/L,
  benign), and pegvaliase-maintenance diet-liberalized classic PKU
  (~500 umol/L) — see §7 for validation detail.

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated classic PKU (natural history)                    | Historical untreated-cohort Phe ranges (Hillert 2020 AJHG epidemiology) |
| 2 | Diet-only, well-controlled classic PKU                     | Vockley 2014 ACMG guideline (target 120-360 umol/L) |
| 3 | Diet + sapropterin, BH4-responsive moderate PKU            | Levy 2007 Lancet PKU-001; Trefz 2009 PKU-003 |
| 4 | Sapropterin trial, classic PKU non-responder (null genotype) | Confirms zero drug effect without BH4-responsive allele |
| 5 | Pegvaliase induction→maintenance, severe classic PKU (diet liberalized) | Thomas 2018 PRISM-1/2; Longo 2018 PAL-003 extension |
| 6 | Pegvaliase with heightened immunogenicity (ADA rise → tolerization) | Gupta 2018 pooled immunogenicity analysis |
| 7 | Untreated mild hyperphenylalaninemia (MHP, benign)         | MHP natural history, no treatment indicated |
| 8 | Maternal PKU pregnancy — poor pre-conception control (high fetal risk) | Rouse 2000 J Pediatr / Koch 2000 Maternal PKU Collaborative Study |
| 9 | Maternal PKU pregnancy — optimized pre-conception control (low fetal risk) | Koch 2000; Lenke & Levy 1980 NEJM |
| 10 | Adolescent diet discontinuation after good childhood control ("off-diet") | Waisbren 2007 meta-analysis; Channon 2004 |

## 4. Shiny dashboard (8 tabs)

1. **Patient / genotype** — adjustable weight/starting age/genotype/scenario.
2. **Drug PK** — sapropterin and pegvaliase plasma concentrations.
3. **Amino-acid PD** — plasma and brain Phe/Tyr/Trp, with the ACMG target band overlay.
4. **Neurocognitive endpoints** — dopamine/serotonin indices, myelin integrity, cumulative IQ deficit, executive-function deficit.
5. **Scenario comparison** — runs all 9 single-stage scenarios with the chosen profile.
6. **Immunogenicity / growth** — ADA titer, tolerization fraction, growth Z-score deficit.
7. **Maternal PKU & safety** — gestational Phe trajectory vs. CHD/microcephaly risk thresholds.
8. **References** — key citations and link to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg pku_qsp_model.dot -o pku_qsp_model.svg
dot -Tpng -Gdpi=150 pku_qsp_model.dot -o pku_qsp_model.png
```

```r
# 2) Simulate scenarios in R (validated with mrgsolve 1.x / R 4.3)
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("pku_mrgsolve_model.R")

# Diet-only, well-controlled classic PKU (infant, 8 kg)
out <- mod %>% param(GENOTYPE=0, WT=8, DIET_PHE_MGKGD=12, TYR_DIET_MGKGD=70) %>%
  mrgsim(end=8760*8, delta=24)
plot(out, c("Plasma_Phe_umolL","Brain_Phe_idx","Dopamine_idx","IQ_deficit_pts"))

# Pegvaliase maintenance, severe classic PKU, diet liberalized (adult, 70 kg)
out2 <- mod %>% param(GENOTYPE=0, WT=70, DIET_PHE_MGKGD=45, TYR_DIET_MGKGD=70) %>%
  ev(amt=20, cmt="SC_PEG", ii=24, addl=366) %>% mrgsim(end=8760, delta=24)
plot(out2, c("Plasma_Phe_umolL","conc_pegvaliase","ADA_titer","Tolerization_frac"))

# 3) Launch the dashboard
shiny::runApp("pku_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Plasma Phe target | ACMG/NIH lifelong guideline | 120-360 umol/L (Vockley 2014 Genet Med) |
| Sapropterin responder | ≥30% plasma-Phe reduction on 8.5-20 mg/kg/d | Levy 2007 Lancet PKU-001; Trefz 2009 PKU-003 |
| Pegvaliase maintenance | 20 mg/d s.c., diet liberalized | Substantial Phe reduction vs. baseline (Thomas 2018 PRISM-1) |
| Pegvaliase immunogenicity | Early ADA rise vs. long-term tolerization | Efficacy improves over ~3-6 months (Longo 2018 PAL-003 extension; Gupta 2018 EBioMedicine) |
| Maternal PKU fetal risk | Periconceptional/gestational Phe control | CHD risk above ~900-1200 umol/L; microcephaly with sustained >600 umol/L (Rouse 2000 J Pediatr; Koch 2000 Mol Genet Metab) |
| IQ outcome | Early-and-continuously-treated vs. late/discontinued diet | Early treatment preserves IQ near population mean (Waisbren 2007 meta-analysis) |
| Executive function | Even well-controlled PKU vs. healthy controls | Persistent prefrontal-dopamine-dependent deficits (Diamond 1997) |

## 7. Model validation (actually executed with mrgsolve 1.x / R 4.3.3)

Unlike a purely illustrative sketch, this model was built, compiled and
simulated end-to-end in this environment (R + mrgsolve + RcppArmadillo
installed from source) to confirm internal consistency before delivery.
Representative validated 8-year (or 1-year/40-week where noted)
steady-state outputs:

| Scenario | Plasma Phe (umol/L) | Notes |
|---|---|---|
| 1. Untreated classic PKU | ~1502 | Matches historical untreated adult ranges (900-2400) |
| 2. Diet-only, controlled classic PKU | ~329 | Within ACMG band |
| 3. Diet + sapropterin (BH4-responsive) | ~304 | vs. ~603 without BH4 at the same diet — clear responder benefit |
| 4. Sapropterin, non-responder | ~329 | Identical to diet-only #2 — correctly shows zero drug benefit |
| 5. Pegvaliase maintenance (diet liberalized) | ~503 | Consistent with PRISM-2 diet-liberalized outcomes |
| 6. Pegvaliase + heightened immunogenicity | ~736 | ADA neutralization measurably attenuates efficacy |
| 7. Untreated MHP | ~203 | Within ACMG band — benign, confirms no treatment needed |
| 8. Maternal PKU, poor control | ~1345 | Exceeds CHD-risk threshold (900) |
| 9. Maternal PKU, optimized control | ~270 | Within ACMG band — low fetal risk |
| 10. Off-diet adolescent (after good childhood control) | Phe reverts to ~1502; IQ deficit +0.2 pts (small, critical period mostly over); executive-function deficit index rises 0→17 | Matches literature: diet discontinuation after childhood causes limited further IQ loss but pronounced executive-function/attention vulnerability |

## 8. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* All PK/PD and BBB-transport parameters are illustrative, literature-
  informed surrogates (absolute "brain Phe/Tyr/Trp index" values are not
  validated against a specific in vivo MRS scale); inter-individual
  variability needs `omega()` blocks for population simulations.
* The hepatic PAH Michaelis-Menten Vmax is scaled linearly per kg body
  weight (liver/enzyme mass approximation) — a simplification of true
  allometric scaling from infancy to adulthood.
* Sapropterin and pegvaliase apparent PK (F, CL, V) are illustrative
  surrogates tuned to produce plausible Tmax/half-life and dose-response,
  not fitted to published population-PK models.
* The anti-drug-antibody/immune-tolerization sub-model is a simplified
  two-state (titer + tolerization-fraction) representation of a genuinely
  complex, patient-variable immunogenicity process.
* Multi-stage scenarios (e.g. #10) require carrying state forward via the
  model's `INIT_*` parameters rather than mrgsolve's `init()`, because
  `$MAIN` sets initial conditions from those parameters on every run — see
  the R-side scenario helpers at the bottom of `pku_mrgsolve_model.R`.

## 9. License

Inherits the repository [LICENSE](../LICENSE).
