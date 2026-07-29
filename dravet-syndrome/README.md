# Dravet Syndrome (SCN1A Developmental & Epileptic Encephalopathy) — QSP Model

> A quantitative systems pharmacology model built around a question that the
> Dravet literature has never cleanly answered: **stiripentol and cannabidiol
> both raise norclobazam, and both are suspected of working through that
> interaction rather than on their own merits. How much of each drug's effect
> is really the drug?**
>
> The model answers it without fitting anything to the answer, because the size
> of the norclobazam route is pinned by a *third*, independent observation —
> and that observation says the route is already saturated before either drug
> is added. The same structure then produces, unprompted, the sign flip that
> makes sodium-channel blockers anticonvulsant in a normal brain and
> catastrophic in Dravet.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT, 211 nodes / 22 clusters) | [`dravet_qsp_model.dot`](dravet_qsp_model.dot) |
| 🖼️ Map (SVG) | [`dravet_qsp_model.svg`](dravet_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`dravet_qsp_model.png`](dravet_qsp_model.png) |
| ⚙️ mrgsolve ODE model (44 states, 21 scenarios, 17 diagnostics) | [`dravet_mrgsolve_model.R`](dravet_mrgsolve_model.R) |
| 🐍 Dependency-free reference implementation | [`dravet_reference_impl.py`](dravet_reference_impl.py) |
| 📊 Shiny dashboard (12 tabs) | [`dravet_shiny_app.R`](dravet_shiny_app.R) |
| 📚 References (154, every PMID resolved through NCBI E-utilities) | [`dravet_references.md`](dravet_references.md) |

**Every number on this page comes from one of two runs, and both are
reproducible:**

```bash
python3 dravet_reference_impl.py     # pure-stdlib RK4 twin, no dependencies
Rscript  dravet_mrgsolve_model.R     # real mrgsolve 2.0.1 -> 17/17 diagnostics pass
```

The two implementations are independent — fixed-step RK4 in Python against
LSODA in mrgsolve's compiled C++ — and they agree. The sodium-channel-blocker
aggravation comes out at **+1400.2%** in both.

---

## 1. The disease in one paragraph

Dravet syndrome is a developmental and epileptic encephalopathy caused, in
about 85% of cases, by a *de novo* loss-of-function variant in **SCN1A**, the
gene for the Nav1.1 sodium channel α-subunit. Development is normal until
seizures begin in the first year of life, usually as prolonged hemiclonic or
generalised convulsions triggered by fever. Multiple seizure types follow,
status epilepticus is common, and development plateaus and then regresses
between ages one and four. Mortality is roughly 10–20% by age 20, about half
of it SUDEP. The pharmacology is unusually sharply defined: valproate and
clobazam are the background, stiripentol, cannabidiol and fenfluramine are the
approved add-ons, and **sodium-channel blockers are contraindicated because
they make the disease dramatically worse** — the opposite of their effect in
almost every other epilepsy.

## 2. The four structural commitments

Everything below follows from four choices. Each is stated so it can be
attacked.

### Commitment 1 — the lesion is in inhibition, not excitation

Nav1.1 is expressed preferentially in **GABAergic interneurons**, above all in
parvalbumin-positive fast-spiking cells (Yu 2006, PMID 16921370). The
excitatory compartment runs on Nav1.2/Nav1.6 and is left intact. In the model
this is a single parameter, `ALLELE`: 0.5 in Dravet, 1.0 in a healthy control.
Nothing else distinguishes the two hosts, and no drug-effect calculation
anywhere in the code asks which host it is in.

### Commitment 2 — interneuron firing is a steep threshold function of sodium reserve

Sustaining 100–200 Hz non-adapting firing costs a fast-spiking interneuron far
more available sodium conductance than a pyramidal cell needs for its 5–20 Hz
output. So the model carries **two** sodium-reserve terms with **two different
steepnesses**: interneuron firing capacity is a Hill function with n = 6 and a
threshold at 0.45, pyramidal excitability a Hill function with n = 1.5.

The healthy operating point (reserve 1.0) sits far above the threshold; the
Dravet operating point (0.5) sits just above it. That is the whole mechanism,
and section 6 below is what it produces.

| | Healthy control | Dravet |
|---|---|---|
| Interneuron firing capacity | **0.9918** | **0.6530** |
| Pyramidal excitability | 0.8285 | 0.8285 |
| Baseline seizures/month on valproate + clobazam | 6.15 | 15.00 |

Dravet retains **65.8%** of normal inhibitory capacity — impaired but
functioning, which is what makes it fragile rather than simply broken.

### Commitment 3 — clobazam is two drugs, and the second one is the important one

Clobazam is demethylated by CYP3A4 to **norclobazam**, which is cleared by
polymorphic CYP2C19 and accumulates to 5–20× the parent. Both are GABA-A
positive allosteric modulators. The model reproduces this without being told:

| Regimen | Clobazam (mg/L) | Norclobazam (mg/L) | N-CLB : CLB |
|---|---|---|---|
| Clobazam 0.5 mg/kg/day alone | 0.200 | 1.600 | **8.0** |
| CYP2C19 intermediate metaboliser | 0.200 | 2.91 | 14.6 |
| CYP2C19 poor metaboliser | 0.200 | 7.70 | **38.5** |

**Consequence: every CYP2C19 inhibitor is an indirect GABAergic drug.**
Stiripentol and cannabidiol both are, so each has a *pharmacokinetic* route and
a *direct* route into the same clinical endpoint. The two are kept as separately
switchable code paths (`PK_ROUTE`, `PD_ROUTE`) because separating them is the
point of the model.

One structural detail matters here and was got wrong in an earlier draft.
Stiripentol's own GABA-A action is given a **separate saturable term** from the
benzodiazepine site, because it binds an α3-preferring site that is not the
benzodiazepine site. Folding the two into one occupancy term makes any
GABAergic add-on structurally incapable of helping a patient already taking
clobazam — which is empirically false, and stiripentol is the proof.

### Commitment 4 — efficacy and sedation share the norclobazam node

Somnolence is driven by benzodiazepine-site occupancy, i.e. by norclobazam.
So a drug that works *through* the interaction pays for its efficacy in
sedation, and a drug that works elsewhere does not. Section 7 prices that.

## 3. Calibration — three fitted anchors, two withheld arms

All comparisons are **placebo-adjusted**, because the placebo response in these
trials is not comparable: +19.2% in Study 1, +26.9% in GWPCARE2, but **−7% in
STICLO**, a 1990s trial of 41 children with a one-month baseline. Comparing raw
reductions across them would compare placebo arms, not drugs.

| Arm | Model | Observed | Source |
|---|---|---|---|
| **[FIT]** Stiripentol 50 mg/kg/day | 75.9% | 76.0 | STICLO: −69% vs +7% |
| **[FIT]** Fenfluramine 0.7 mg/kg/day | 55.7% | 55.7 | Study 1: 74.9 vs 19.2 |
| **[FIT]** Cannabidiol 20 mg/kg/day | 22.2% | 22.2 | GWPCARE1+2 average |
| **[PRED]** Fenfluramine 0.2 mg/kg/day | **24.1%** | 23.1 | Study 1: 42.3 vs 19.2 |
| **[PRED]** Cannabidiol 10 mg/kg/day | **20.4%** | 21.8 | GWPCARE2: 48.7 vs 26.9 |

The two withheld arms are dose-response predictions with no freedom left in
them, and they are informative in opposite directions. Fenfluramine turns out
to be **linear** in exposure over the therapeutic range, so the 0.2 mg/kg arm
follows arithmetically from the 0.7 mg/kg fit. Cannabidiol is the opposite: its
direct route is **already saturated at 10 mg/kg/day**, which is why GWPCARE2
found 10 and 20 mg/kg/day nearly indistinguishable.

The measured drug-drug interactions are reproduced too, and these are what
`KI_*` were fitted to:

| | Model | Observed |
|---|---|---|
| Norclobazam × on cannabidiol 20 mg/kg/day | 5.72 | 6.0 (Geffrey 2015) |
| Clobazam parent × on cannabidiol | 1.61 | 1.6 (Geffrey 2015) |
| Norclobazam × on stiripentol | 2.48 | ~2.5 (Jullien 2015) |

### Study 1504 reproduced by enrichment, not by assertion

Study 1504 gave fenfluramine 0.4 mg/kg/day on stiripentol-inclusive regimens
and required ≥6 convulsive seizures in a 6-week baseline **while already taking
stiripentol**. It therefore recruited stiripentol-*refractory* patients.
Simulating a good stiripentol responder and adding fenfluramine is a different
experiment, so the model screens a virtual cohort the way the trial screened a
real one:

| | Model | Observed |
|---|---|---|
| Eligible after screening on stiripentol | 48% (60/124) | — |
| Median residual seizures/month on stiripentol | 8.5 | — |
| Median further reduction on fenfluramine 0.4 | **52.1%** | 49.0 (placebo-adjusted) |
| ≥50% responders | **55%** | 54% |

Nothing in that section was tuned to it. The fenfluramine parameters came from
Study 1 in an unselected population, the stiripentol parameters from STICLO,
and even the **0.4 mg/kg dose cap** follows from the CYP inhibition constants:
stiripentol inhibits the enzymes that clear fenfluramine, so 0.4 mg/kg with
stiripentol delivers 46.1 ng/mL against 54.9 ng/mL for 0.7 mg/kg alone — the
cap recovers **83.8%** of the uncapped exposure.

## 4. The central result — the interaction explains almost none of it

The size of the norclobazam route is set by an observation from a cohort that
received **neither** stiripentol **nor** cannabidiol. Hashi 2015 (PMID
25323806) studied 50 patients on low-dose add-on clobazam and found median
norclobazam of **1103 ng/mL** in those with ≥90% seizure reduction, against
341 and 570 ng/mL in the less well controlled groups. So the benzodiazepine-site
exposure–response is already near-saturated around 1.1 mg/L.

A standard clobazam dose in this model reaches **1.60 mg/L**. The route is past
its useful range *before* any interacting drug arrives. Forcing norclobazam
upward and changing nothing else:

| Norclobazam | Concentration | Seizures/month | vs baseline |
|---|---|---|---|
| ×1.0 | 1600 ng/mL | 15.00 | — |
| ×2.5 (the stiripentol interaction) | 3993 ng/mL | 14.76 | **−1.6%** |
| ×6.0 (the cannabidiol interaction) | 8988 ng/mL | 14.77 | **−1.5%** |
| ×10 | 13009 ng/mL | 14.81 | −1.3% |

A six-fold rise in norclobazam — the entire measured magnitude of the
cannabidiol interaction — buys essentially nothing. The curve even **turns
back up** above about ×4, and that is not a numerical artefact: sustained
occupancy recruits receptor tolerance, so past a point extra norclobazam is
self-defeating as well as merely useless.

Deleting each route in turn:

| Drug | Both routes | PK route only | Direct route only | **PK share** |
|---|---|---|---|---|
| Stiripentol 50 mg/kg/day | 75.9% | 2.7% | 74.6% | **3.6%** |
| Cannabidiol 20 mg/kg/day | 22.2% | 4.5% | 16.5% | **20.2%** |

**Stiripentol's effect is essentially all its own.** The suspicion that STICLO
measured a norclobazam boost does not survive contact with the arithmetic: the
same model that reproduces the cannabidiol interaction at ×6 says a ×2.5 rise
cannot deliver 76 percentage points. Cannabidiol is the drug with the larger
share of borrowed effect — a fifth of it — which is the reverse of the usual
suspicion.

### The consequence that was already tested, and held

If the PK route carried stiripentol's effect, response would track CYP2C19
genotype, since that genotype sets how much 2C19 is left to inhibit. The model
says the opposite — stiripentol response should be nearly
**genotype-independent**:

| CYP2C19 | Stiripentol | Cannabidiol 20 | Baseline N-CLB |
|---|---|---|---|
| Normal | 75.9% | 22.2% | 1.60 |
| Intermediate | 75.6% | 20.3% | 2.91 |
| **Poor** | **75.4%** | 19.1% | 7.70 |

Kouga 2015 (PMID 24819914) looked. In 11 Japanese Dravet patients, stiripentol
response showed no significant relationship to CYP2C19 genotype; 3 of 8 had no
benefit *despite* a rise in norclobazam, and one benefited while norclobazam
**fell**. The authors concluded this "suggests a significant anti-convulsant
action of STP." That is a confirmation the model was not fitted to.

The clobazam-free stratum says the same thing from the other direction:

| Drug | On clobazam | Clobazam-free | Retained |
|---|---|---|---|
| Stiripentol 50 | 75.9% | 73.0% | **96.1%** |
| Cannabidiol 20 | 22.2% | 13.4% | **60.4%** |

## 5. The sodium-channel-blocker paradox is an output

Same drug, same dose, same equations. Only `ALLELE` differs.

| Host | Baseline MCSF | On lamotrigine 5 mg/kg/day | Change | Interneuron capacity |
|---|---|---|---|---|
| Healthy control | 6.15 | 5.63 | **−8.4%** | 0.8071 |
| Dravet | 15.00 | 225.04 | **+1400.2%** | 0.0614 |

The sign flips because a 30% use-dependent block leaves a healthy interneuron
at reserve 0.70 — still far above the 0.45 threshold — but takes a Dravet
interneuron to 0.35, below it, and the Hill-6 term collapses. Meanwhile the
genuine anticonvulsant action (use-dependent block of high-frequency
propagation) is unchanged in both hosts; in the healthy brain it is the only
thing that happens, which is why the drug works there.

Dose-response, and the interaction that makes it worse without changing the
prescription:

| Lamotrigine mg/kg/day | Css (mg/L) | MCSF | vs baseline |
|---|---|---|---|
| 0 | 0.00 | 15.35 | +2.3% |
| 1 | 1.40 | 28.49 | +89.9% |
| 2 | 2.80 | 65.47 | +336.5% |
| 5 | 7.00 | 225.04 | +1400.2% |

Valproate inhibits UGT1A4, so the *same* lamotrigine dose reaches 7.00 mg/L
with valproate against 3.65 mg/L without — and the aggravation goes from
+520.6% to +1400.2%. The prescription looks identical.

Febrile susceptibility falls out of the same threshold, with a
mutant-selective thermal term:

| Host | Peak temperature | Interneuron capacity nadir | Peak rate | Ratio |
|---|---|---|---|---|
| Healthy control | 38.52 °C | 0.9910 | 20.8 | **3.37×** |
| Dravet | 38.52 °C | 0.3858 | 105.2 | **6.92×** |

## 6. Zorevunersen and a ceiling made of arithmetic

Antisense poison-exon skipping (the TANGO approach) upregulates SCN1A from the
**intact** allele. It recovers the transcripts that allele wastes on the
non-productive exon-20N isoform and nothing more:

| | Normalised Nav1.1 function |
|---|---|
| Untreated Dravet | 0.500 |
| After 3 intrathecal doses | 0.536 |
| Complete poison-exon skipping (the ceiling) | **0.746** |
| Healthy | 1.000 |

So this class cannot reach a normal phenotype however well it works — that is
a property of allele arithmetic, not of dose. In the model three quarterly
doses give a 24.7% reduction in seizure frequency **at zero sedation cost**,
because it acts at the root node and never touches the benzodiazepine site.

## 7. Therapeutic index — what each route costs

| Regimen | MCSF | vs reference | Somnolence | pp per unit sedation |
|---|---|---|---|---|
| Clobazam 0.5 alone (reference) | 15.35 | — | 0.388 | — |
| Clobazam 1.0 (dose doubled) | 15.03 | +2.3% | 0.477 | 25 |
| + stiripentol 50 | 3.59 | +76.7% | 0.496 | 712 |
| + cannabidiol 20 | 11.66 | +24.2% | 0.555 | 145 |
| + fenfluramine 0.7 | 6.64 | +56.8% | 0.388 | **free** |
| + zorevunersen | 9.17 | +40.4% | 0.388 | **free** |

"Free" means the arm reduced seizures with no added sedation at all, because it
does not touch the benzodiazepine site. Doubling clobazam is the opposite
trade: it buys 2.3 percentage points for a large sedation increment, because
norclobazam is already on the flat part of its curve.

This is directly actionable, and it is what STICLO reported informally when
side-effects resolved on cutting comedication in 12 of 21 patients. Adding
cannabidiol while *cutting* clobazam improves both endpoints at once:

| Clobazam (mg/kg/day) | Reduction | Somnolence | Norclobazam |
|---|---|---|---|
| 0.50 | 22.3% | 0.555 | 9.16 |
| 0.30 | 25.5% | 0.518 | 5.49 |
| 0.20 | **28.7%** | **0.479** | 3.66 |

## 8. Long-horizon outcomes (5 years, febrile illness each quarter)

| Arm | DQ | Cumulative seizures | SE/year | 5-y SUDEP risk | Weight z | ALT |
|---|---|---|---|---|---|---|
| No add-on | 49.2 | 1275 | 12.79 | 17.4% | 0.00 | 25 |
| Stiripentol | **75.3** | 295 | 1.99 | 4.0% | −0.60 | 25 |
| Fenfluramine 0.7 | 71.6 | 522 | 4.01 | 7.1% | −1.10 | 25 |
| Fenfluramine + clobazam cut | 73.8 | 556 | 4.34 | 7.6% | −1.10 | 25 |
| Zorevunersen q6mo | 63.9 | 779 | 5.91 | 10.4% | 0.00 | 25 |
| Cannabidiol 20 | 52.8 | 977 | 9.00 | 13.4% | 0.00 | **116** |
| Lamotrigine (contraindicated) | **22.0** | 14018 | 424 | **95.4%** | 0.00 | 25 |

Developmental loss is driven by seizure burden **and** by sedation, which is
why "fenfluramine + clobazam cut" beats fenfluramine alone on DQ despite
slightly more seizures. The transaminase column locates the cannabidiol
hepatic signal where it is actually seen — in valproate co-medication:

| Regimen | ALT (U/L) | × baseline |
|---|---|---|
| Cannabidiol 20, no valproate | 50.8 | 2.03 |
| Cannabidiol 20 + valproate 30 | **112.7** | **4.51** |
| Valproate 30 alone | 25.0 | 1.00 |

Benzodiazepine tolerance is a slow state, and starting a naive patient on
clobazam shows it: seizures/month fall to 13.71 by day 30, then climb back to
**16.81 by day 360** on an unchanged dose as receptor availability falls from
1.00 to 0.76.

## 9. What was fitted, and what was not

**Fitted — 8 parameters, each to a named published quantity:**
`KI_CBD_2C19` and `KI_7OH_2C19` (norclobazam ×6.0, Geffrey 2015);
`KI_CBD_3A4` (clobazam ×1.6, Geffrey 2015); `KI_STP_2C19` (norclobazam ×2.5,
Jullien 2015); `EMAX_PAM` (clobazam withdrawal doubles seizure frequency);
`EMAX_CBD` (GWPCARE1/2); `EMAX_5HT` (Study 1); `EMAX_STP_SITE` (STICLO).

**Pinned by an independent observation:** `EC50_PAM_NCLB`, from Hashi 2015.
The central conclusion rests on this one number, from a cohort exposed to
neither drug being decomposed.

**Not fitted — these fall out of the structure:** the sodium-channel-blocker
sign flip; the febrile susceptibility ratio; fenfluramine 0.2 mg/kg and
cannabidiol 10 mg/kg; Study 1504 under trial-faithful enrichment; the
zorevunersen ceiling; the therapeutic-index ordering; and the CYP2C19
genotype-independence of stiripentol response.

## 10. Known misses, stated rather than buried

1. **The stiripentol response distribution is too tightly clustered.** The
   model puts 94% of virtual patients over the 50% responder line against an
   observed 71% in STICLO, and reaches nobody at 100% against an observed 9/21
   seizure-free. Both errors point the same way: there is no genuine
   non-responder subpopulation. Kouga 2015 found 3 of 8 patients with no
   stiripentol benefit at all despite rising norclobazam. A mixture model with
   an explicit non-responder fraction would fit both tails; a single
   log-normal on `EMAX_STP_SITE` cannot.
2. **There is no clobazam-free stiripentol trial in existence.** STICLO was run
   entirely on valproate + clobazam. The 96.1% retention figure in section 4 is
   an extrapolation from the model's structure, not a reproduction of data. The
   nearest real evidence is Kouga 2015, and it is eleven patients.
3. **The Hashi 2015 exposure–response comes from mixed epilepsy, not Dravet,**
   and it is retrospective. If the true Dravet norclobazam response saturates
   much higher than 1.1 mg/L, the headroom argument weakens and the PK share
   rises.
4. **The magnitude of the lamotrigine aggravation is not calibrated.** The
   direction and the dose-dependence are structural results; +1400% is not a
   measured quantity, and no trial will ever measure it. Treat the sign as the
   claim and the magnitude as an illustration.
5. **`EMAX_ADEN` is deliberately zero.** The adenosine (ENT1), GPR55 and TRPV1
   nodes are all driven by the same cannabidiol concentration and are not
   separately identifiable from any clinical dataset. Loading part of the
   fitted effect onto an unfitted accessory term would misreport where the
   evidence is, so `EMAX_CBD` carries the whole direct effect and the model
   does not claim to know which downstream node it travels through.
6. **Development, SUDEP and weight are the softest parts of the model.** They
   are calibrated to order-of-magnitude plausibility, not to cohort data, and
   should be read as structure rather than prediction.

## 11. Model structure

44 ODEs:

- **PK (21 states)** — clobazam (2-compartment) → norclobazam (2-compartment);
  stiripentol with Michaelis-Menten elimination, giving the observed
  dose-supraproportionality; cannabidiol → 7-OH-CBD → 7-COOH-CBD;
  fenfluramine → norfenfluramine; valproate; a lamotrigine-like blocker;
  rescue diazepam; intrathecal antisense in CSF and brain.
- **Enzymology (4 states)** — CYP2C19, CYP3A4, CYP1A2/2B6/2D6 and UGT1A4 as
  turnover states, so inhibition has an onset and an offset rather than being
  instantaneous.
- **Target biology (4 states)** — Nav1.1 function from productive SCN1A
  transcript; GABA-A receptor availability (the tolerance substrate);
  serotonergic and adenosine tone.
- **Circuit (3 states)** — kindling, chronic interneuron attrition, core
  temperature.
- **Clinical (12 states)** — infection burden, cumulative seizures, status
  epilepticus, developmental quotient, somnolence, weight, SUDEP hazard, ALT,
  valve index, filtered seizure frequency.

The seizure rate has a physiological ceiling (`LAM_MAX`), with events crowding
into status epilepticus as it saturates. This is not cosmetic: without it the
sodium-channel-blocker arm diverges numerically to 10^19 seizures/month instead
of predicting what actually happens to those patients, which is status.

---

## Reproducing this page

```bash
# 211-node mechanistic map
dot -Tsvg dravet_qsp_model.dot -o dravet_qsp_model.svg
dot -Tpng -Gdpi=150 dravet_qsp_model.dot -o dravet_qsp_model.png

# every number above, no dependencies
python3 dravet_reference_impl.py

# the same model under real mrgsolve, with 17 diagnostics and 21 scenarios
Rscript dravet_mrgsolve_model.R

# interactive dashboard (needs shiny, ggplot2, dplyr, tidyr, DT)
Rscript -e 'shiny::runApp("dravet_shiny_app.R", port = 8080)'
```

## ⚠️ Disclaimer

This is an **educational and research** QSP model. It is built from public
literature and has not been independently validated or qualified for any
regulatory purpose. It must not be used for clinical decisions, prescribing, or
dose selection. In particular, the lamotrigine arm exists to test the model's
structure — sodium-channel blockers are contraindicated in Dravet syndrome, and
nothing here should be read as exploring their use.
