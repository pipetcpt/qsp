# Gastrointestinal Stromal Tumour (GIST) — QSP model

> A 49-compartment quantitative systems pharmacology model of GIST posed **not**
> as "one tumour, one kinase activity, one IC50, one resistance state" but as a
> **set-cover problem over a clone population**, on top of a **quiescent
> reservoir the drug does not kill** and an **imaged mass that is mostly not
> viable tumour**. The model's headline result is that the INTRIGUE genotype
> crossover — sunitinib 15.0 vs ripretinib 4.0 months in one ctDNA subgroup and
> ripretinib 14.2 vs sunitinib 1.5 in the other, with a null in ITT — falls out
> of five published *in-vitro* potency ratios that were never fitted to any
> clinical endpoint.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT, 197 nodes / 21 clusters / 292 edges) | [`gist_qsp_model.dot`](gist_qsp_model.dot) |
| 🖼️ Map (SVG) | [`gist_qsp_model.svg`](gist_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`gist_qsp_model.png`](gist_qsp_model.png) |
| ⚙️ mrgsolve ODE model (49 compartments, 26 scenarios) | [`gist_mrgsolve_model.R`](gist_mrgsolve_model.R) |
| 🐍 Dependency-free Python twin (runs, prints every number, self-checks) | [`gist_python_twin.py`](gist_python_twin.py) |
| 📊 Shiny dashboard (12 tabs) | [`gist_shiny_app.R`](gist_shiny_app.R) |
| 🧪 mrgsolve scenario driver (26 scenarios, endpoint helpers) | [`gist_scenarios.R`](gist_scenarios.R) |
| 📚 References (57 PubMed-verified) | [`gist_references.md`](gist_references.md) |

---

## 1 · The problem the model is built to solve

GIST is the disease where targeted therapy worked first and best, and it is
therefore the disease where the standard cartoon is most entrenched:

> KIT mutation → constitutive kinase → proliferation → imatinib occupies the ATP
> pocket → tumour shrinks → a resistance mutation appears → switch drug.

One tumour, one kinase activity, one IC50, one resistance state. Six things in
the trial record cannot all be true in that cartoon at once:

1. **Two second-line drugs tie overall and invert by genotype.** INTRIGUE:
   ripretinib 8.0 vs sunitinib 8.3 months in ITT, no difference. In the ctDNA
   analysis of the same trial, patients whose plasma carried *only* KIT exon
   11+13/14 mutations got 15.0 months from sunitinib and 4.0 from ripretinib;
   patients with *only* exon 11+17/18 got 14.2 from ripretinib and 1.5 from
   sunitinib, and the two groups were **mutually exclusive**
   ([Heinrich 2024](https://pubmed.ncbi.nlm.nih.gov/38182785/),
   [Bauer 2022](https://pubmed.ncbi.nlm.nih.gov/35947817/)).
2. **Imatinib rechallenge works after imatinib has failed.** RIGHT: PFS 1.8 vs
   0.9 months on placebo, in patients who had progressed on imatinib, sunitinib
   and everything else ([Kang 2013](https://pubmed.ncbi.nlm.nih.gov/24140183/)).
3. **Dose escalation helps exactly one genotype.** 800 mg reduces the relative
   risk of progression by 61% in KIT exon 9 and does nothing elsewhere
   ([Debiec-Rychter 2006](https://pubmed.ncbi.nlm.nih.gov/16624552/),
   [MetaGIST 2010](https://pubmed.ncbi.nlm.nih.gov/20124181/)).
4. **The same node has opposite drug ranks.** PDGFRA D842V is the most
   imatinib-resistant genotype there is and the most avapritinib-sensitive
   (ORR ~88%, [Heinrich 2020](https://pubmed.ncbi.nlm.nih.gov/32615108/)).
5. **PET and CT disagree by weeks.** FDG SUV collapses within 24–48 h while
   tumour size barely moves, and "stable disease" is a common best response in
   patients who live for years
   ([Van den Abbeele 2012](https://pubmed.ncbi.nlm.nih.gov/22381410/)).
6. **Years of deep response do not eradicate the disease.** Interrupting
   imatinib in responders gives PFS 6.1 / 7.0 / 12.0 months (after 1 / 3 / 5
   years of treatment) versus 27.8 / 67.0 / not reached on continuation — and
   the time to regrowth gets *longer* the longer the drug was given
   ([Blay 2024](https://pubmed.ncbi.nlm.nih.gov/39127063/)).

A single-clone model can be fitted to any one of these. It cannot hold 1 and 2
at the same time, because a tumour that has "become imatinib-resistant" cannot
respond to imatinib; and it cannot hold 5 and 6 at the same time, because a
tumour whose mass has fallen 50% ought to be 50% gone.

## 2 · Three structural commitments

### C1 · The tumour is a population, not a size

The state is a **clone vector**, not a burden. Four clones, each with a cycling
and a quiescent compartment:

| clone | identity | who covers it |
|---|---|---|
| 0 | primary genotype only (exon 11 / exon 9 / D842V / SDH-deficient) | imatinib (genotype-dependent), all others partially |
| 1 | KIT exon 13/14 **ATP-binding pocket** — V654A, T670I | **sunitinib**; ripretinib and regorafenib poorly |
| 2 | KIT exon 17/18 **activation loop** — D816, D820, N822K, Y823D | **ripretinib, regorafenib, avapritinib**; sunitinib poorly |
| 3 | KIT-independent bypass (KIT loss, MET/AXL/FGFR, RAS, CDKN2A) | nobody |

Each clone has an EC50 for each of five drugs, and the coverages multiply:

```
A_i = Π_d  1 / (1 + C_d / EC50_{d,i})          SIG_i = 1 − kitdep_i·(1 − A_i)
```

Subclones are generated **per cell division**, so almost all of them exist before
the first dose: the pre-existing fraction is `μ × TURNOVER`, where `TURNOVER` is
the number of divisions per surviving cell accumulated over the pre-diagnostic
history. Resistance is therefore **selected, not induced**, and efficacy is a
**set-cover problem**: progression is the growth of the least-covered clone, so
two drugs can swap rank between genotypes while their averages are identical.

### C2 · Occupancy is fast, killing is slow, and the drug is mostly cytostatic

Three time scales are separated explicitly, and a fourth variable is the imaged
mass:

| variable | time constant | read out by |
|---|---|---|
| KIT pathway output `SIG` | hours | FDG-PET SUV |
| cycling fraction | days | Ki67 |
| quiescent reservoir | months–years | what happens when you stop |
| imaged mass | months | CT / RECIST |

Cells whose KIT signal falls below threshold mostly **stop rather than die**
(`KAMAX` is capped and needs *deep* suppression: Hill 6 on `1 − SIG`). Those that
stop enter a compartment the drug cannot kill, which dies at 0.065%/day. Exit
from quiescence requires mitogenic signal **steeply** (`SIG³`) and is slowed
further by a dormancy-depth state, so the reservoir survives years of therapy and
wakes within weeks of stopping.

And a **synchronous drug-induced die-off overwhelms phagocytic clearance**,
leaving hyalinised/myxoid residue that resorbs at 0.3%/day — whereas the trickle
of physiological apoptosis is cleared efficiently. That single asymmetry is why
observations 5 and 6 are the same fact: the residual mass at nadir is mostly
**not** viable tumour, so re-expanding it to +20% SLD takes months even though
the cells regrow at their untreated rate.

### C3 · Exposure matters where the genotype puts the EC50

Imatinib PK is explicit — CYP3A4 autoinduction, the equipotent metabolite
CGP74588, and AGP binding entering as `EC50_eff = EC50·(AGP/AGP0)^0.9`. At
400 mg the model reaches a steady-state trough of **1205 ng/mL** (parent +
metabolite; the measured value is ~1100, and the exposure–response threshold
Demetri identified was 1110). That sits ~8-fold above the exon 11 EC50 and only
~1.4-fold above the exon 9 EC50 — i.e. on **opposite sides of the proliferation
threshold**. Escalation moves exon 9 across it and does nothing for exon 11.

---

## 3 · What comes out: first line

| scenario | PFS | best SLD |
|---|---|---|
| KIT exon 11, 400 mg | **28.5 mo** | −48.5% |
| KIT exon 11, 800 mg | 29.4 mo (**+3%**) | −47.5% |
| KIT exon 9, 400 mg | **8.7 mo** | 0% |
| KIT exon 9, 800 mg | **33.1 mo** (+279%) | −28.4% |
| KIT exon 9, 400 → 800 only at progression | 35.2 mo total on imatinib | — |
| exon 11, lowest exposure quartile (trough 399 ng/mL) | 32.4 mo | −15.6% |
| rifampicin co-administration (trough 461) | 29.2 mo | −42.2% |
| PDGFRA D842V, imatinib 400 mg | 3.0 mo | 0% |
| PDGFRA D842V, avapritinib 300 mg | 31.3 mo | −49.1% |
| SDH-deficient, imatinib 400 mg | 9.2 mo | 0% |

Observed: exon 11 median PFS ~25 months (S0033); exon 9 raises the relative risk
of progression 171% and 800 mg reduces it 61% *in exon 9 only*, which for a
25-month exon 11 median implies ~9 months for exon 9 at 400 mg and ~23 at 800.
D842V imatinib ORR ~0%, avapritinib ORR 88%.

The escalation asymmetry (**+279% in exon 9, +3% in exon 11**) is not a fitted
subgroup effect. It is threshold geometry: only one genotype's EC50 lies between
the two doses.

Note what the model does **not** say. Escalating an exon 9 patient only at
progression gives 35.2 months of total time on imatinib against 33.1 for starting
at 800 mg — the threshold has to be crossed eventually, and crossing it late
loses nothing. That is what practice already does, and it is consistent with
MetaGIST finding a PFS but no overall-survival advantage for the higher dose.

## 4 · What comes out: the INTRIGUE crossover  [C1]

Second-line PFS, by which secondary-mutation class the tumour is able to carry:

| | sunitinib | ripretinib |
|---|---|---|
| KIT exon 11 + 13/14 only | **6.2 mo** | 2.8 mo |
| KIT exon 11 + 17/18 only | 2.5 mo | **31.3 mo** |
| both classes in one tumour | 3.0 mo | 5.3 mo |

Observed: 15.0 vs 4.0, and 1.5 vs 14.2 months. The model reproduces **the
inversion** — the product of the two sunitinib/ripretinib ratios is 0.18, i.e.
the ranking flips by a factor of 5.5 — and under-states the magnitude in the
ATP-pocket subgroup.

Two things about this are worth stating precisely:

- **Nothing here was fitted to INTRIGUE.** The relevant numbers are the EC50
  columns, taken from [Serrano 2019](https://pubmed.ncbi.nlm.nih.gov/30792533/),
  [Smith 2019](https://pubmed.ncbi.nlm.nih.gov/31085175/) and
  [Heinrich 2008](https://pubmed.ncbi.nlm.nih.gov/18955458/): sunitinib is potent
  against V654A/T670I and weak against the activation loop; ripretinib and
  regorafenib are the mirror image.
- **The model does not try to reproduce the ITT median, and says why.** The two
  ctDNA subgroups were *mutually exclusive*, so the ITT median is the median of a
  mixture of two opposite-signed effects — a sampling property of the enrolled
  population, not a property of any patient. The model's statement about a
  patient who carries **both** classes is that they do badly on either drug
  (3.0 and 5.3 months), which is a different and testable claim.

## 5 · What comes out: later lines, and why rechallenge is not a paradox

| line | model | observed |
|---|---|---|
| 3rd line regorafenib | 2.5 mo | 4.8 (GRID) |
| 4th line ripretinib | 6.7 mo | 6.3 (INVICTUS) |
| 4th line best supportive care | 2.1 mo | 0.9–1.0 |
| 4th line imatinib rechallenge | 2.5 mo | 1.8 (RIGHT) |

Clone composition at the start of the fourth line: **primary mutation only 15.8%,
ATP-pocket 20.7%, activation loop 63.4%, bypass 0.1%**. A tumour that has "failed"
imatinib three lines ago still contains a substantial sensitive population — which
is Kang's own conclusion from RIGHT ("the disease continues to harbour many clones
that are sensitive to kinase inhibitors") stated as a state variable.

The ceiling on the whole sequence is the **KIT-independent clone**, which nothing
covers and which becomes dominant at around the fourth line. That is why each
successive line is shorter, without any per-line efficacy parameter.

## 6 · What comes out: PET vs CT, and ctDNA

| | 48 h | day 7 | 8 weeks |
|---|---|---|---|
| PET SUV | **−61%** | | |
| Ki67 | | −85% | |
| SLD | −0.0% | | −15% |

At one year of imatinib **94.5% of the surviving viable cells are quiescent** and
the **viable fraction of the imaged mass is 36%** — the mass on CT is mostly
hyalinised residue around a small, largely non-cycling population.
Resistance-mutation ctDNA crosses 1% VAF **5.9 months before** RECIST calls
progression (VAF 5.4% at progression) — which is the practical
point of C1, because that window is when the second line could be chosen by
genotype rather than by label order.

## 7 · What comes out: stopping the drug  [C2]

| randomised after | interruption | continuation | viable fraction of mass | dormancy depth |
|---|---|---|---|---|
| 1 year | **4.3 mo** | 16.5 mo | 0.36 | 0.17 |
| 3 years | **3.6 mo** | 13.5 mo | 0.46 | 0.11 |
| 5 years | **10.0 mo** | 47.4 mo | 0.28 | 0.09 |

Observed (BFR14): 6.1 vs 27.8, 7.0 vs 67.0, 12.0 vs not reached.

Two honest points about this table. First, the phenotype of each randomised
cohort is set from **its own continuation arm**, because BFR14 randomised only
patients who were still non-progressing — so the interruption arm is a prediction
and the continuation arm is the control that defines the patient. Second, the
model *requires* an essentially resistance-free tumour to stay progression-free
for 3–5 years, and that is itself a prediction rather than an excuse: only 50 of
the 434 patients enrolled in BFR14 reached the 3-year randomisation and 27 the
5-year one. The model's reading of the observed 12-month interruption PFS in the
5-year cohort is that it is a property of **those patients**, not of the extra
two years of drug.

## 7b · Adjuvant therapy: delay, not cure

| adjuvant duration | recurrence becomes detectable | delay vs 12 months |
|---|---|---|
| 12 months | 31.4 mo | — |
| 36 months | 53.0 mo | +21.5 mo |
| 60 months | 54.5 mo | +23.1 mo |

Observed (SSGXVIII): 5-year RFS 53.0% for 12 months versus 71.4% for 36 months,
converging to 41.8% versus 52.5% at ten years. There is **no cure term anywhere
in the model** — the occult burden is suppressed while the drug is given and
resumes its untreated kinetics afterwards — and that alone produces both the
separation and the convergence. Note also that the 36→60 month step buys almost
nothing here (+1.6 months): by three years the occult burden has already been
driven down to the floor set by the drug-insensitive reservoir, which decays at
0.065%/day, so two more years of imatinib removes only another third of it.

## 8 · A prediction that goes against intuition

Cover the activation-loop clone *before* it is selected, or *after*?

| strategy | per line | total |
|---|---|---|
| imatinib → ripretinib → sunitinib | 28.5 \| 5.3 \| 7.4 mo | **41.2 mo** |
| imatinib 400 + ripretinib 100 upfront → sunitinib | 32.2 \| 4.8 mo | **37.0 mo** |

The combination buys **+3.7 months of first-line PFS and loses 4.2 months
overall**. This is not the intuitive answer and it is not a fitted one: covering
the activation-loop clone early removes the very heterogeneity the later lines
exploit, so the tumour arrives at the uncoverable KIT-independent clone sooner. A
set-cover model predicts that upfront combination is a bad trade *even though it
lengthens the first line* — which is what the field has found the hard way, and
which is a falsifiable statement about any future combination trial.

## 9 · Falsification: one clone whose potency drifts

The fair competitor is not "no resistance". It is resistance as a **continuous
loss of potency under drug pressure** — the standard acquired-resistance
formulation — **refitted to the same first-line anchor** (`KDRIFT` = 2.1×10⁻³/day
by bisection, giving 28.5 months). Then:

| discriminating test | clone vector | single drifting clone | verdict |
|---|---|---|---|
| INTRIGUE genotype crossover | 6.2/2.8 and 2.5/31.3 mo | 2.5 vs 2.8 mo, **the same pair in every subgroup** | **decisive** — with one clone the subgroup parameters are literally inert, so the crossover is not merely unfitted, it is unreachable |
| 4th-line ripretinib (INVICTUS, observed 6.3) | 6.7 mo | 2.1 mo | **decisive** — with one clone there is nothing broad for a broad-spectrum inhibitor to be good at |
| imatinib rechallenge (RIGHT) | 2.5 vs 2.1 mo | 9.0 vs 1.8 mo | **does not discriminate** — a 7-fold potency shift leaves imatinib partly active, so the drifting model also predicts a benefit, and a larger one |
| escalation confined to exon 9 | +279% / +3% | +290% / +27% | **weak** — the asymmetry is mostly C3 (genotype EC50), not C1 |

So of the four tests, two are decisive, one is null and one is weak. The report
says so rather than claiming a clean sweep.

## 10 · What the model gets wrong

- **M1 · Imatinib trough and outcome in exon 11.** Demetri found TTP 11.3 months
  in the lowest trough quartile versus >30 months above it. In this model exon 11
  first-line PFS is almost **exposure-independent** above a trough of ~250 ng/mL
  (32.4 months at a trough of 399), because progression is driven by a clone
  imatinib never covers at any dose; only the *depth* of response tracks exposure
  (−15.6% vs −48.5%). Either the clinical association is not a potency threshold,
  or exon 11 EC50s are more spread than a single number.
- **M2 · The rechallenge effect is too small** (2.5 vs 2.1 months against 1.8 vs
  0.9 observed), and as §9 shows the observation does not discriminate between
  the two structures anyway.
- **M3 · Third-line regorafenib and second-line agents in a both-classes tumour
  are under-predicted**, because the trial populations are mixtures of
  single-class patients while the model's "unselected" patient carries both
  classes at once.
- **M4 · BFR14 magnitudes are 1–3 months short** in the 1- and 3-year arms, and
  the 3- and 5-year cohorts can only be represented as resistance-free tumours.
- **M5 · No overall survival and no spatial structure**, so no focal progression
  and no cytoreductive surgery of a single progressing nodule — which is exactly
  the setting where polyclonality matters most clinically.

## 11 · Parameters: what is fitted and what is not

Eight of ~150 parameters are fitted, and they are named in `FITTED`:

| parameter | fitted to |
|---|---|
| `KPMAX`, `KD0` | untreated volume doubling time (62 days) |
| `TURNOVER` | KIT exon 11 first-line PFS |
| `FNEC`, `KRES` | depth of the RECIST response |
| `KDQ`, `PHI_D` | BFR14 1-year interruption arm |
| `N_MICRO` | 12-month adjuvant arm |

Everything called a prediction above is downstream of those eight plus potency
ratios that were never fitted to a clinical endpoint. The drug × clone EC50
matrix, the genotype EC50s, all five drugs' PK, and every toxicity slope come
from the literature in [`gist_references.md`](gist_references.md).

## 11b · Toxicity, read inside each line of therapy

| regimen | ANC min | ΔMAP | HFSR | oedema | TSH | dose intensity |
|---|---|---|---|---|---|---|
| imatinib 400 mg | 3.35 | 0.0 | 0.05 | **0.99** | 1.00 | 1.00 |
| sunitinib 50 mg 4/2 | 2.41 | **+11.2** | 0.42 | 0.46 | **1.87** | 1.00 |
| sunitinib 37.5 mg daily | 2.75 | +11.6 | 0.37 | 0.46 | 2.00 | 1.00 |
| regorafenib 160 mg 3/1 | 3.28 | +13.2 | **2.07** | 0.03 | 1.35 | **0.58** |

Each row is read inside its own line of therapy, not over the whole run, because
oedema from the imatinib phase would otherwise be attributed to whatever came
next. Sunitinib's free T4 falls to 0.63 of baseline and TSH rises 1.87-fold;
regorafenib is the only agent whose toxicity forces a sustained dose reduction
(to 58% of the intended dose), and that reduction feeds straight back into
coverage. Second-line PFS is 3.0 months for both 50 mg 4/2 and 37.5 mg
continuous — the intermittent schedule lets every clone regrow for two weeks in
six, and the continuous one is gentler but never lets up.

## 12 · Running it

```bash
# the executable reference: prints every number above and asserts them
python3 gist_python_twin.py            # full report + self-checks
python3 gist_python_twin.py --quiet    # self-checks only

# render the map
dot -Tsvg gist_qsp_model.dot -o gist_qsp_model.svg
dot -Tpng -Gdpi=150 gist_qsp_model.dot -o gist_qsp_model.png
```

```r
library(mrgsolve)
mod <- mread("gist_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt = 400, ii = 1, addl = 900, cmt = "A_IM"), end = 900)
plot(out, SLD + SUV + N_TOT + VAF_R ~ time)

shiny::runApp("gist_shiny_app.R")     # 12-tab explorer
```

The Python twin is the numerical reference. It has no dependencies beyond the
standard library, integrates the identical 49-state system with RK4, and its
dosing is driven by an integer day counter so results do not depend on the step
size — a bug that step-aligned dose testing hid until it was found and fixed.

---

*Educational and research model. Not validated for clinical decision-making.*
