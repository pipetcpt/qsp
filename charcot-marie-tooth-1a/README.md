# Charcot-Marie-Tooth disease type 1A (CMT1A) — QSP model

**샤르코-마리-투스병 1A형 · 유전자 용량(gene dosage) 신경병증**

| File | What it is |
|---|---|
| [`cmt1a_qsp_model.dot`](cmt1a_qsp_model.dot) · [`.svg`](cmt1a_qsp_model.svg) · [`.png`](cmt1a_qsp_model.png) | Mechanistic map — 24 modules, 273 labelled nodes |
| [`cmt1a_mrgsolve_model.R`](cmt1a_mrgsolve_model.R) | 59-ODE mrgsolve model, 32 scenarios, 17 diagnostics |
| [`cmt1a_shiny_app.R`](cmt1a_shiny_app.R) | 10-tab interactive dashboard |
| [`cmt1a_references.md`](cmt1a_references.md) | 152 references, every PMID resolved through NCBI E-utilities |

Validated under **mrgsolve 2.0.1**. `source("cmt1a_mrgsolve_model.R"); run_diagnostics()`
reproduces every number below.

---

## Why this disease needed a different kind of model

CMT1A is the most common inherited neuropathy (~1 in 5000), and it has an
aetiology of unusual simplicity: a 1.4-Mb tandem duplication at 17p11.2
delivers **three copies of *PMP22* instead of two**. Not a loss of
function. Not a toxic misfolded species (that is CMT1E, a different
disease caused by point mutations in the same gene). Just the wrong
number of copies of one structural myelin protein.

That forces an architecture most disease models do not have. **The
therapeutic target is a ratio to be restored, not a quantity to be
minimised** — and human genetics supplies both walls of the window at
integer points:

| copies | *PMP22* mRNA | phenotype |
|---|---|---|
| 1 | 0.5× | HNPP — tomacula, recurrent pressure palsies |
| 2 | 1.0× | normal |
| 3 | 1.5× | **CMT1A** |
| 4 | 2.0× | severe infantile / Dejerine-Sottas-like |

Because 1.5 × (1 − 0.667) = 0.5, the far wall of the window sits at
**67% knockdown** and the optimum at **33%**. Neither number is fitted.
They are arithmetic.

---

## Five structural commitments

**S1 — the disposal system is already at ~80% duty cycle in health.**
PMP22 folds inefficiently even in wild-type Schwann cells; most of what
is translated is degraded before it ever reaches myelin. The model
therefore does not treat a 50% rise in synthesis as producing a 50% rise
in anything. Misfolded flux passes a *saturating* ERAD step, and what
escapes accumulates as an aggregate pool. This is what makes 3 copies a
disease and 2 copies not, and it generates the CN=4 runaway without
being asked to.

**S2 — the two arms of the U-curve are different pathways.** CMT1A and
HNPP are not mirror images. Excess dosage acts through the **aggregate
burden** (proteostatic stress → c-Jun → dedifferentiation →
demyelination). Deficient dosage acts through **membrane stoichiometry**
(gm < 1 → tomacula → focal pressure palsy). The consequence, which the
model produces rather than assumes, is that *the far wall is invisible on
the disability scale* — see D04.

**S3 — length-dependence is geometry, not biology.** Transported support
decays as exp(−L/(v·τ)); locally delivered glial lactate does not; demand
does not depend on L. Nothing about the axon to the foot is biologically
different from the axon to the shoulder. It is 1 m long instead of 0.35 m.

**S4 — the motor unit reserve is large, silent and one-way.** Collateral
sprouting lets a surviving unit adopt the fibres of several lost ones.
"Onset" is reserve exhaustion, not disease onset.

**S5 — conduction velocity is a developmental fossil.** Internodal length
is eroded overwhelmingly while the myelination window is open. Nothing
in the model forces NCV to be flat afterwards; it emerges from putting
the erosion term under a developmental weight.

---

## The model

59 ODEs in days. Every scenario starts from the **same solved healthy
newborn nerve** and integrates copy number over 70 years; no adult
phenotype is typed in anywhere.

The healthy reference state is *solved in `$MAIN`*, not specified: the
folding efficiency at CN=2 (φ = 0.200), the membrane PMP22 reference, the
healthy aggregate burden (a fixed point of escape-vs-autophagy) and the
c-Jun/EGR2 fixed point are all computed by iteration at run time. That is
why the drift test comes out exactly zero:

```
D01  healthy CN=2 over 70 years
     max |mRNA − 1|  = 0.000e+00      max |gm − 1| = 0.000e+00
     max |NCV − 52|  = 0.000e+00
```

### The generated adult phenotype (CN=3, age 40)

| quantity | model | published range |
|---|---|---|
| motor NCV | 17.0 m/s | 15–25 m/s, uniform |
| MUNE, leg | 27% of normal | severely reduced |
| ulnar CMAP | 6.5 mV | reduced |
| sural SNAP | absent | usually absent |
| CMTNS-R | 17.0 | mean ~15, range 2–30 |
| ONLS | 1.4 | ~2–3 at trial entry |
| calf fat fraction | 24.6% (+1.05 %/y at 40–60) | ~30%, +1.2 %/y |
| plasma NfL | 1.7× healthy control | ~2× |

---

## Results

### D02 — one integer generates four diseases

| CN | mRNA | gm | aggregates | tomacula | NCV | MUNE | CMTNS | pressure palsies |
|---|---|---|---|---|---|---|---|---|
| 1 (HNPP) | 0.500 | 0.810 | 0.001 | 0.99 | 50.5 | 86% | 2.3 | **0.89/y** |
| 2 (normal) | 1.000 | 1.000 | 0.047 | 0 | 52.0 | 86% | 2.1 | 0 |
| 3 (CMT1A) | 1.500 | 1.085 | 0.409 | 0 | 17.0 | 27% | 17.0 | 0 |
| 4 (homozygous) | 2.000 | 1.133 | **4.898** | 0 | 0.6 | 0% | 33.4 | 0 |

HNPP has essentially **no aggregate pathology** and near-normal
conduction velocity: its lesion is purely stoichiometric. CMT1A has
essentially **no tomacula**: its lesion is purely proteostatic. The two
diseases come out of the same equations at different points and they do
not resemble each other, which is what S2 predicted.

Note also `gm` at CN=3 is only **1.085**, not 1.5. The folding bottleneck
buffers the membrane pool; the excess ends up in aggregates rather than
in myelin. The model's claim is therefore that CMT1A is driven by
mistrafficked PMP22, not by an excess of PMP22 in compact myelin.

### D03 — NCV really is a fossil

```
CMT1A NCV:  age 2  24.5 → age 5  18.9 → age 20  17.15 → age 70  16.78
            = −0.0074 m/s per YEAR
over a 15-month trial NCV moves −0.009 m/s.  Test-retest SD is 2–3 m/s.
over the same window CMTNS moves +0.227/y and calf fat fraction +0.60 %/y
```

A trial powered on nerve conduction velocity in adult CMT1A is measuring
a variable whose time-derivative is roughly 300× below its own
measurement noise. This is not a statement about drug effect. It is a
property of the endpoint.

### D04 — the therapeutic window, and where its far wall actually shows up

Lifelong constant knockdown, outcomes at age 50:

| knockdown | mRNA | gm | aggregates | tomacula | NCV | **CMTNS** | **pressure palsies/y** |
|---|---|---|---|---|---|---|---|
| 0% | 1.500 | 1.085 | 0.409 | 0 | 16.9 | **18.91** | 0 |
| 20% | 1.200 | 1.041 | 0.127 | 0 | 39.2 | 2.81 | 0 |
| **35%** | 0.975 | 0.994 | 0.041 | 0.01 | 52.0 | **2.58** ← min | 0.005 |
| 50% | 0.750 | 0.927 | 0.009 | 0.23 | 51.6 | 2.63 | 0.21 |
| 65% | 0.525 | 0.825 | 0.001 | 0.87 | 50.7 | 2.77 | **0.79** |
| 80% | 0.300 | 0.646 | 0.000 | 2.50 | 48.8 | 3.10 | **2.25** |

Two things fall out.

The **CMTNS optimum is at 35%** against an arithmetic prediction of
33.3%, and at 65% knockdown the mRNA is 0.525 against the 1-copy level
of 0.500. The window's walls land where the gene copy numbers are.

But the far wall **does not appear on the disability scale at all**.
CMTNS at 80% knockdown is 3.10 against 2.58 at the optimum — a
difference no trial would ever see. The cost of overshooting shows up
entirely as a *different toxicity*: 2.25 pressure palsies per year.
A dose-finding study powered on CMTNS would happily over-dose all the
way into iatrogenic HNPP and read the result as flat-topped efficacy.
That is a direct consequence of S2 and it is the single most
actionable output of this model.

### D05 — CMT1A sits 7% below a proteostatic collapse threshold

| mRNA | aggregates | differentiation | NCV | CMTNS |
|---|---|---|---|---|
| 1.40 | 0.283 | 0.618 | 23.6 | 4.6 |
| **1.50 (CMT1A)** | 0.409 | 0.489 | 17.0 | **17.0** |
| 1.55 | 0.494 | 0.412 | 14.1 | 30.8 |
| 1.60 | 0.606 | 0.326 | 11.5 | 33.3 |
| 1.70 | **3.393** | 0.007 | 1.2 | 33.4 |

Autophagic clearance of aggregates has a *maximum* (it is inhibited by
its own substrate), so above a critical escape flux the aggregate pool
has no stable low branch. That threshold lands at mRNA ≈ 1.6–1.7, and
CMT1A operates at 1.5.

This is the model's account of two otherwise unrelated observations:
why the *homozygous* duplication is not simply "twice as bad" but
categorically different, and why the same duplication produces CMTNS
scores from 2 to 30 in different people. A ±5% shift in expression moves
the simulated patient from CMTNS 4.6 to CMTNS 30.8. **Falsifiable
prediction: measured *PMP22* mRNA in CMT1A patients should span a narrow
band around 1.4–1.55 and should correlate very tightly with severity.**

### D06 — four null ascorbic acid trials, explained before efficacy is discussed

| oral dose | steady-state plasma | achievable *PMP22* suppression |
|---|---|---|
| 1000 mg/d | 62.9 µmol/L | 0.45% |
| 1500 mg/d | 66.0 µmol/L | 0.53% |
| 3000 mg/d | 77.5 µmol/L | 0.75% |
| 4000 mg/d | 88.8 µmol/L | 0.91% |

Two independent saturating mechanisms — SVCT1 intestinal absorption
collapsing from ~90% to ~30% bioavailability, and a sharp renal
reabsorption threshold — clamp plasma ascorbate. **A 2.67-fold dose
difference produces a 1.34-fold plasma difference and a 0.38
percentage-point difference in target engagement**, against a modelled
therapeutic requirement of 33%.

The trials that were run (CMT-TRIAAL, CMT-TRAUK, Micallef, Burns,
Verhamme, Lewis) explored essentially none of the dose–response axis.
The doses tested were pharmacodynamically indistinguishable from each
other and from placebo. That conclusion needs no assumption about
whether ascorbate lowers PMP22 at all: raising the model's assumed Emax
to 100% changes the *absolute* numbers but not the ratio, because both
saturations are upstream of the pharmacology.

### D07 — PXT3003: the model under-predicts, and says so

| | modelled 15-month ONLS difference | published |
|---|---|---|
| dose level 1 | −0.094 | not separated from placebo |
| dose level 2 (high) | −0.099 | **−0.37** |

The model achieves 15.1% and 18.8% *PMP22* suppression at the two dose
levels — in line with the preclinical claim — and turns that into an
ONLS effect roughly **3.7× smaller** than PLEO-CMT reported. It also
fails to separate the two dose levels, which the trial did.

Both are reported rather than tuned away. Either the model's
mRNA-to-disability transfer is too shallow at 15 months, or PXT3003 does
something the PMP22 axis alone does not capture, or the trial's effect
size is fragile. The model cannot distinguish these, and the same model
run for 15 years instead of 15 months gives ΔONLS = −1.52 — so its
disagreement is about *timescale*, not about direction.

### D08 — the silent decade

| age | MUNE | mean unit size | dorsiflexion | CMTNS |
|---|---|---|---|---|
| 10 | 92.6% | 1.08× | 99.3% | 2.4 |
| 20 | 65.2% | 1.52× | 95.5% | 9.9 |
| 30 | 41.4% | 2.37× | 88.9% | 14.4 |
| 40 | 27.0% | 3.55× | 79.1% | 17.0 |
| 60 | 13.1% | 5.00× | 48.4% | 20.7 |

By the time dorsiflexion strength has fallen 5%, 38% of motor units are
gone. By the time it has fallen 20%, 73% are gone. An adult trial is
operating on what is left after the reserve has been spent, which is
why D12 finds start age dominates dose.

### D09 — geometry, not biology

At age 40, all four motor classes obey the same equation and differ only
in length:

| class | length | margin | surviving axons |
|---|---|---|---|
| intrinsic foot | 1.00 m | 0.973 | 7.9% |
| tibialis anterior | 0.85 m | 1.111 | 27.0% |
| hand intrinsics | 0.75 m | 1.218 | 45.5% |
| proximal | 0.35 m | 1.812 | 83.7% |

And because length is the variable, **height is a dose**:

```
stature factor 0.88 / 1.00 / 1.15  →  CMTNS at 50 = 15.2 / 18.9 / 25.2
                                       (identical genotype)
```

### D10 — vincristine: one equation, three hosts

Vincristine 2 mg IV weekly × 4 at age 20. The drug lowers axonal
transport velocity by up to 40%; nothing else about it is CMT-specific.

| host | foot margin | surviving axons | ΔONLS |
|---|---|---|---|
| healthy (CN=2) | 2.397 → 1.536 | 0.924 → 0.920 (**−0.5%**) | +0.001 |
| CMT1A (CN=3) | 0.973 → **0.596** | 0.436 → 0.028 (**−93.5%**) | **+0.654** |
| HNPP (CN=1) | 2.225 → 1.437 | 0.923 → 0.918 (−0.6%) | +0.001 |

The catastrophe is not written into the model as a contraindication. The
healthy nerve absorbs the same transport hit because it has 2.4× the
support it needs; the CMT1A nerve is already at 0.97 and has nowhere to
go. HNPP — despite being a *PMP22* disease — is protected, because its
lesion never touched the margin.

### D13 — the endpoint dominates everything else about the trial

24-month trial, 50% slowing, 80% power, α = 0.05:

| endpoint | placebo change | treated change | **N per arm** |
|---|---|---|---|
| calf MRI fat fraction | +1.028% | −0.001% | **38** |
| CMTNS-R | +0.402 | −0.062 | 236 |
| ONLS | +0.168 | −0.008 | 410 |
| plasma NfL | −1.020 | −0.002 | 545 |
| ulnar CMAP | −0.056 mV | −0.004 mV | 4,644 |
| **motor NCV** | −0.014 m/s | −0.000 m/s | **513,512** |

A 38-patient MRI trial and a 513,512-patient NCV trial test the same drug
with the same effect. Note also that plasma NfL *falls* on placebo,
because the rate of axon loss declines once there is less left to lose —
so NfL is at its most informative in young, rapidly progressing patients
and least informative in exactly the adult cohorts that get recruited.

### D12 — when you start beats what you give

| | CMTNS at 65 | MUNE at 65 |
|---|---|---|
| oligonucleotide from age 5 | **8.09** | 66.4% |
| from age 30 | 15.00 | 34.9% |
| from age 50 | 17.11 | 20.5% |
| age 30–40 then stopped | 19.65 | 16.0% |
| untreated | 21.43 | 11.3% |

Ten years of perfect therapy withdrawn at 40 retains about a quarter of
the benefit of continuing. Reserve is not recoverable.

---

## Negative and self-refuting results

These are reported because they were found, not because they help.

**D14 — cluster 4 of the map is refuted by the model built from it.**
The map was drawn with Schwann-cell c-Jun as a second interior optimum:
the repair Schwann cell is bad for myelin but supplies trophic support
to the axon, so there should be an intermediate best. Sweeping the c-Jun
gain over a decade gives a **monotone** response — CMTNS 7.7 → 33.4 as
c-Jun rises, with no minimum. The trophic gain never outweighs the loss
of lactate/MCT1 support from the myelinating cells that dedifferentiated.
Cluster 4 has been relabelled with its own refutation rather than
redrawn.

**D15 — the long tissue half-life is not protective, and it does not
need to be.** The model was built expecting that a slow-turnover
oligonucleotide would matter because it flattens the dosage trajectory
inside a narrow window. It does flatten it — monthly dosing holds mRNA
at 0.88–0.96 while six-monthly dosing swings it 0.63–1.33 — and the
outcome is **identical** (CMTNS at 65: 15.01 vs 14.98). The downstream
integrators average over months regardless. The useful version of this
result is the opposite of the hypothesis: infrequent dosing is safe.

**D07 — PXT3003 under-predicted ~3.7×** and the two dose levels are not
separated. See above.

**The CN=4 phenotype is over-severe.** The model gives a homozygous
duplication carrier NCV 0.6 m/s and MUNE 0%, i.e. functionally no
peripheral nervous system. Reported homozygous *PMP22* duplication
phenotypes are severe and Dejerine-Sottas-like but not that. The
aggregate runaway is real in the equations; its magnitude is not
calibrated against anything, because there is nothing to calibrate it
against.

**The 5% cliff.** Going from 0% to 5% knockdown takes CMTNS at 50 from
18.9 to 7.8. That is a large therapeutic claim resting on the steepness
of the ERAD-escape function (S1) and of the c-Jun/EGR2 switch. It is the
same steepness that explains the clinical variability (D05), so it is not
free to remove — but if the model is over-promising anywhere, it is here.

**The 80% figure at the heart of S1 has no verified citation.** The
qualitative claim (PMP22 misfolds, aggregates, saturates the proteasome)
is well supported; the specific fraction is a model assumption. It is
listed in `cmt1a_references.md` under "claims with no verified source"
and is the first thing to check before using this model.

**Four PMIDs in the first draft were wrong.** They were written from
memory and pointed at unrelated papers. Every citation in this directory
was subsequently resolved through NCBI E-utilities; the four errors are
tabulated at the top of the references file.

---

## What the model is for, and what it is not for

It is a hypothesis-structuring and trial-design tool. The three outputs
worth taking seriously are (1) the shape of the dose–response for a
*PMP22*-lowering agent and the fact that its safety wall is invisible on
the efficacy endpoint, (2) the ranking of trial endpoints by required
sample size, and (3) the pharmacokinetic account of the ascorbate
trials, which requires no efficacy assumption at all.

It is not calibrated against individual patient data, the free
parameters are numerous and listed as such in the model header, and no
part of it should inform clinical decisions.
