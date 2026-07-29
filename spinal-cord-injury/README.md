# Traumatic Spinal Cord Injury (SCI) — QSP Model

> A quantitative systems pharmacology model of traumatic spinal cord injury,
> built around a single question: **why has almost every neuroprotective drug
> that worked in the laboratory failed in the clinic?** The model's answer is not
> "the drugs don't work". It is that the cascade they target is over in hours,
> that the endpoint they are measured on is a steep saturating function of the
> tissue they save, and that the patients easiest to enrol are exactly the
> patients in whom a real tissue benefit is arithmetically undetectable.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`sci_qsp_model.dot`](sci_qsp_model.dot) |
| 🖼️ Map (SVG) | [`sci_qsp_model.svg`](sci_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`sci_qsp_model.png`](sci_qsp_model.png) |
| ⚙️ mrgsolve ODE model | [`sci_mrgsolve_model.R`](sci_mrgsolve_model.R) |
| 📊 Shiny dashboard | [`sci_shiny_app.R`](sci_shiny_app.R) |
| 🐍 Reference implementation (no dependencies) | [`sci_reference_model.py`](sci_reference_model.py) |
| 📄 Reference output (every number below) | [`sci_reference_output.txt`](sci_reference_output.txt) |
| 📚 References | [`sci_references.md`](sci_references.md) |

**Scale:** 238-node / 18-cluster / 378-edge mechanistic map · 43-ODE mrgsolve
model with 160 annotated parameters · 18 prebuilt scenarios · 12-tab Shiny
dashboard · 70 references (PMID-verified where the record could be confirmed,
PubMed search links otherwise, with an explicit list of what is *not* sourced).

Every quantitative claim in this README is produced by
`python3 sci_reference_model.py`, a dependency-free re-implementation of the same
43-state system, and the full run is saved in `sci_reference_output.txt`. No
number here is asserted from memory.

---

## 1. Disease in one paragraph

Traumatic spinal cord injury begins with an **irreversible** mechanical insult —
contusion, persistent compression, axonal shear — that no drug will ever reach,
and continues as a **secondary injury cascade** that is, in principle, entirely
addressable. Residual compression plus cord swelling raise intraspinal pressure
(ISP) while neurogenic shock lowers mean arterial pressure, and the resulting
collapse in spinal cord perfusion pressure (SCPP = MAP − ISP) drives ischemia →
ATP failure → SUR1-TRPM4/AQP4-mediated cytotoxic edema, which raises ISP further:
a **compression–edema–ischemia vicious cycle**. Inside it, glutamate floods the
extracellular space within minutes, Ca²⁺ overload activates calpain and opens the
mitochondrial permeability transition pore, peroxynitrite and Fenton chemistry
peroxidize membrane lipids, neutrophils arrive at ~24 h and macrophages over days
to weeks, and caspase-3 kills neurons and oligodendrocytes for a fortnight after
the accident. What is left is a cystic cavity walled off by a CSPG-rich glial
scar, a demyelinated and thinned rim of spared white matter, and a spinal cord
below the lesion that slowly reorganizes — restoring reflexes that produce
spasticity, detrusor overactivity and autonomic dysreflexia while a closing
critical period allows rehabilitation to recruit whatever descending drive
survived. Function at one year is not a property of the injury alone: it is the
product of what survived, how well it conducts, and how much of it the patient's
circuits learned to use.

## 2. Mechanistic map — 18 clusters, 238 nodes

1. **Primary mechanical injury** (t = 0, irreversible, *not* a drug target): mechanism, impact energy, compression, contusion, shear, neurological level, initial AIS grade, central hemorrhagic necrosis, primary axonal loss, spared subpial rim
2. **Systemic haemodynamics & cord perfusion**: neurogenic shock, sympathectomy, bradycardia, MAP, vasopressor support, the MAP ≥ 85 mmHg protocol, ISP, SCPP, autoregulation loss, SCBF, ischemic penumbra, reperfusion, duraplasty
3. **Cord swelling & ionic dysregulation**: ATP depletion, Na⁺/K⁺-ATPase failure, SUR1–TRPM4, AQP4, cytotoxic and vasogenic edema, BSCB breakdown, MMP-9, tight-junction loss, rostral ascending edema
4. **Glutamatergic excitotoxicity**: vesicular release, EAAT/GLT-1 reversal, NMDA/AMPA-kainate, persistent Na⁺ current, Ca²⁺ overload, calpain, α-II-spectrin breakdown, mPTP, cytochrome-c, nNOS
5. **Oxidative stress**: superoxide, peroxynitrite, hydroxyl radical, free iron, lipid peroxidation, 4-HNE/MDA, protein nitration, ferroptosis, Nrf2/HO-1, SOD/catalase
6. **Neuroinflammation**: DAMPs/HMGB1, TLR4–NF-κB, NLRP3, IL-1β/TNF-α/IL-6, CXCL1/CCL2, neutrophils, monocytes, microglia, the M1/M2 balance, lymphocytes, debris clearance, SCI-induced immune depression
7. **Cell death programs**: caspase-8/9/3, Bax/Bcl-2, necroptosis, autophagy, neuron death, oligodendrocyte apoptosis, Wallerian degeneration, retrograde dieback, cystic cavitation, late syringomyelia
8. **Demyelination & remyelination**: myelin loss, conduction block, juxtaparanodal K⁺ exposure, Nav redistribution, OPC/NG2 proliferation, Schwann-cell remyelination, thin sheaths
9. **Glial scar & regeneration failure**: reactive astrogliosis, fibrotic core, CSPG, PTPσ/LAR, Nogo-A, MAG/OMgp, NgR1/p75/LINGO-1, RhoA–ROCK, growth-cone collapse, PTEN/mTOR, KLF/RAG silencing, perineuronal nets
10. **Spared tracts, plasticity & circuit reorganization**: CST and reticulospinal sparing, propriospinal relays, effective descending drive, collateral sprouting, synaptogenesis, the critical period, BDNF/NT-3, locomotor CPG, cortical remapping, maladaptive plasticity
11. **Methylprednisolone PK/PD** including start time, cumulative exposure and the complication set (pneumonia, GI haemorrhage, hyperglycaemia, myopathy)
12. **Riluzole PK/PD** (CYP1A2, persistent Na⁺-current block, release inhibition, enhanced uptake, ALT elevation)
13. **Minocycline & glibenclamide PK/PD** (M1 blockade, MMP inhibition, caspase inhibition; SUR1–TRPM4 blockade, hypoglycaemia)
14. **Regenerative & neuromodulation interventions**: intrathecal anti-Nogo-A, Rho-ROCK inhibition, chondroitinase, activity-based rehabilitation, body-weight-supported and robotic gait training, epidural stimulation, FES cycling, cell therapy, decompression timing
15. **Symptomatic pharmacology**: baclofen (and its weakness trade-off), pregabalin, antimuscarinics/onabotulinumtoxinA, α-blockade
16. **Secondary systemic complications**: reflex reorganization, spasticity, neuropathic pain, neurogenic bladder and DSD, autonomic dysreflexia, respiratory failure, sublesional osteoporosis, pressure injury, VTE, cardiometabolic decline, heterotopic ossification
17. **Clinical endpoints**: ISNCSCI motor and sensory scores, AIS conversion, SCIM III, WISCI II, 10-m walk, GRASSP, pain NRS, Ashworth, MRI lesion volume and tissue bridges, MEP/SSEP, CSF NfL/GFAP, QoL, mortality/LOS
18. **Quantitative concepts** encoded in the ODEs, drawn as explicit notes on the map: the cascade flux Φ(t), the window integral, the threshold mapping, surrogate dissociation, the two clocks, and recovery-creates-pathology

## 3. The 43 ODE states

| Group | Compartments |
|---|---|
| Drug PK (9) | `MP_CENT`, `MP_PERIPH`, `RIL_DEPOT`, `RIL_CENT`, `MINO_CENT`, `GLY_CENT`, `NOGO_ITH`, `BAC_CENT`, `PGB_CENT` |
| Exposure ledger (1) | `MP_AUC` |
| Haemodynamics (2) | `SHOCK`, `MAP` |
| Acute drivers (2) | `INJ` (necrosis/DAMP release), `ABLOCK` (spinal shock) |
| Cascade (8) | `ISCHEMIA`, `EDEMA`, `GLU`, `CAI`, `ROS`, `NEUT`, `CYTO`, `APOP` |
| Immune cells (3) | `M1`, `M2`, `DEBRIS` |
| Tissue (5) | `NEURON`, `OLIG`, `OLIG_CAP`, `AXON`, `CAVITY` |
| Scar (2) | `GFAP`, `CSPG` |
| Recovery (4) | `CRIT`, `PLAST`, `ATRO`, `MOTOR` |
| Late complications (5) | `REFLEX`, `SPAST`, `NPAIN`, `BLADDER`, `AD_TRIG` |
| Organ systems (2) | `RESP`, `BMD` |

Two structural choices carry most of the model's behaviour:

- **`AXON` only ever falls; `PLAST` is what rises.** Effective descending drive
  is `CONN = AXON × myelin efficiency × conduction factor × (1 + plasticity +
  stimulation gain)`. Anatomy and function are separate variables, which is the
  only way a plasticity-directed therapy can change outcome without touching the
  imaging surrogate — and vice versa.
- **Only *supra-threshold* death signalling destroys tissue** (`APOP_THR`,
  `ISCH_THR`). Without this the system drains axons forever at its own tiny
  chronic steady state; with it, the chronic state is genuinely quiescent and
  late decline has to be put in deliberately.

## 4. Natural history is the calibration target

Standard care in the model = decompression at 24 h, rehabilitation intensity 0.5
from day 14. For a C5 lesion, initial → 1-year ISNCSCI total motor score:

| Baseline AIS | at presentation | day-3 nadir | 1 year (rehab 0.5) | change | 1 year (rehab 1.0) |
|---|---|---|---|---|---|
| A | 11.2 | 11.1 | 12.1 | **+0.9** | 12.5 (+1.3) |
| B | 17.3 | 16.6 | 22.1 | **+4.8** | 24.3 (+7.0) |
| C | 38.9 | 36.4 | 51.3 | **+12.4** | 55.7 (+16.8) |
| D | 67.8 | 63.7 | 79.7 | **+11.9** | 82.8 (+15.0) |

The **day-3 nadir** is not tissue loss — it is reversible conduction block from
edema and ischemia. That is why patients improve from their admission exam at
all, and why the admission exam systematically understates the surviving
anatomy. About 90% of the modelled change is complete by 3–6 months.

## 5. Three results the model exists to produce

### 5.1 Timing dominates dose — and toxicity does not care about timing

The NASCIS methylprednisolone regimen (30 mg/kg bolus + 5.4 mg/kg/h × 23 h) given
at increasing delays, in a C5 AIS C patient:

| Start | Δ ISNCSCI motor at 1 yr | Δ lesion volume | Cumulative MP exposure | Complication index |
|---|---|---|---|---|
| 1 h | **+2.56** | −4.9% | 25.7 au·d | 0.391 |
| 3 h | +2.26 | −4.6% | 25.7 | 0.391 |
| 6 h | +1.88 | −4.2% | 25.7 | 0.391 |
| 9 h | +1.59 | −3.8% | 25.7 | 0.391 |
| 12 h | +1.36 | −3.5% | 25.7 | 0.391 |
| 24 h | **+0.83** | −2.4% | 25.7 | 0.391 |
| 48-h infusion from 3 h | +2.68 | −5.8% | 47.3 au·d | **0.542 (1.39×)** |

Benefit falls **3.1-fold** from 1 h to 24 h. The complication index is
**identical to three decimal places** at every delay, because efficacy is
`∫ E(t)·Φ(t) dt` against a cascade flux that decays, while exposure is
`∫ C(t) dt` against nothing. Doubling the infusion to 48 h buys 1.39× the
complication index for +0.42 motor points — the drug has already spent its
window by the time the second day of infusion begins. This is the entire NASCIS story
reproduced from mechanism rather than from post-hoc subgroups: in the same
patient at AIS A, the benefit is +0.31 at 1 h and +0.08 at 24 h — a real drug
effect that no trial could measure.

### 5.2 The surrogate and the endpoint are not the same variable

`MOTOR` is a steep saturating Hill function of `CONN` (TH50 = 0.22, HILL = 2.2).
The value of **doubling** a patient's spared drive:

| Spared drive | doubled to | ISNCSCI motor | → | gain | detectable vs ~5-point noise? |
|---|---|---|---|---|---|
| 0.03 | 0.06 | 11.1 | 14.9 | **+3.8** | no |
| 0.06 | 0.12 | 14.9 | 28.8 | +13.9 | yes |
| 0.11 | 0.22 | 26.1 | 55.0 | **+28.9** | yes |
| 0.20 | 0.40 | 50.3 | 81.0 | +30.7 | yes |

Run the same drug across grades and the dissociation becomes explicit. Riluzole
(14 d, RISCIS regimen) reduces the lesion index by **8.7–9.4% in every grade** —
essentially grade-independent, as a tissue-level mechanism should be — while its
effect on the endpoint is not:

| Baseline AIS | Δ lesion index | Δ ISNCSCI motor |
|---|---|---|
| A | −8.7% | +0.18 |
| B | −8.8% | +0.90 |
| C | −9.0% | **+1.77** |
| D | −9.4% | +1.20 |

A neuroprotectant can therefore be mechanistically successful and clinically
negative *at the same time*, and the difference between those two outcomes is
**who was enrolled**, not whether the drug works. AIS A patients are the easiest
to recruit and consent, have the least ambiguous baseline, and sit exactly where
the derivative of the endpoint with respect to the mechanism is smallest.

### 5.3 Two clocks — and recovery that creates its own pathology

The neuroprotection clock runs in hours. The plasticity clock runs in weeks and
is shut by the scar (`CRIT` decays at `K_CRIT·(1 + K_CRIT_SCAR·CSPG)`). The
**same** rehabilitation dose, moved 106 days later (C5 AIS C):

| Arm | plasticity pool | 1-year ISNCSCI motor |
|---|---|---|
| no rehabilitation | 0.154 | 43.7 |
| rehab d14–180 | **0.465** | **55.7** |
| rehab d120–286 (identical dose) | 0.176 | 47.5 |

**+8.2 motor points for the same intervention delivered earlier** — a window
effect that has nothing to do with the acute cascade and cannot be recovered by
doing more rehabilitation later.

Meanwhile the sublesional reflex reorganization that underlies functional
recovery is itself the source of the chronic complications. The same bladder
distension applied at different times after a T4 injury:

| Trigger day | reflex maturity | peak systolic BP | surge |
|---|---|---|---|
| 7 | 0.30 | 115 mmHg | +18 |
| 30 | 0.77 | 161 mmHg | +57 |
| 90 | 0.97 | 189 mmHg | **+85** |
| 180 | 0.98 | 191 mmHg | +87 |

A harmless nursing event in week 1 is a hypertensive emergency at 6 months, and
nothing about the lesion changed in between. The same variable drives Ashworth
from 0.0 → 2.7 and pain NRS from 0.2 → 3.6 over the first three months.

## 6. Interventions, ranked as the model ranks them

C5 AIS A, 1 year, versus standard care (lesion index 0.323, motor 12.1):

| Intervention | lesion index | Δ motor | note |
|---|---|---|---|
| Bundle: decompression 8 h + MAP protocol + riluzole | **0.118 (−63%)** | **+2.4** | interventions compose |
| MAP ≥ 85 mmHg × 7 d | 0.155 (−52%) | +1.9 | cheapest, largest single effect |
| Minocycline 7 d | 0.295 (−9%) | +0.2 | wide window, small flux |
| Riluzole 14 d | 0.295 (−9%) | +0.2 | |
| Methylprednisolone ≤ 3 h | 0.306 (−5%) | +0.2 | plus a 0.391 complication index |
| Glibenclamide 3 d | 0.306 (−5%) | +0.2 | edema-directed |
| Early decompression < 12 h | 0.315 (−2%) | +0.1 | large in AIS C: +1.9 (4 h vs 24 h) |
| Anti-Nogo-A + intensive rehab | 0.323 (0%) | +0.8 | changes function, not the surrogate |
| Intensive rehab from d14 | 0.323 (0%) | +0.4 | |
| Baclofen + pregabalin from d30 | 0.323 (0%) | **−2.6** | Ashworth 2.74 → 2.09, NRS 3.58 → 2.92 |
| No decompression at all | 0.993 (+207%) | −2.1 | permanent compression is a different disease |

Two things worth reading twice. **Haemodynamic management, which involves no new
molecule, dominates every drug arm in this model**, because it acts on the SCPP
term that gates the entire cascade for as long as the cord is swollen — that is a
testable prediction, and `SCPP`/`ISP` are captured as outputs so it can be
argued with. And **symptomatic pharmacology is not free**: baclofen lowers
Ashworth through the same GABA-B mechanism that produces weakness, so the
modelled motor score falls by 2.6 points. The trade-off is structural, not a
side effect.

For an incomplete injury, epidural stimulation from day 120 is the largest single
effect anywhere in the model: AIS C 51.3 → **61.9** motor points, with no change
whatsoever in lesion index or axon sparing — it amplifies residual drive rather
than saving tissue.

## 7. Organ-system consequences at one year

| Patient | vital capacity | sublesional BMD | muscle atrophy | motor |
|---|---|---|---|---|
| C4 AIS A, no rehab | 29% predicted | 0.78 | 0.60 | 1.7 |
| C4 AIS A, rehab + FES | 30% | **0.90** | 0.10 | 3.2 |
| C5 AIS A, no rehab | 54% | 0.78 | 0.60 | 11.6 |
| C5 AIS A, rehab + FES | 54% | **0.90** | 0.10 | 12.9 |
| T10 AIS A, rehab + FES | 95% | 0.91 | 0.09 | 53.3 |

Sublesional bone loss (−22% at one year without loading) is driven by absent
weight-bearing, not by the lesion — so it is the one endpoint that FES and
standing improve **without any change in descending drive**. It is also the
reason the model tracks `BMD` at all: an intervention that never moves the motor
score can still move a hard outcome.

## 8. Running it

```r
# mrgsolve
library(mrgsolve)
mod <- mread_cache("sci_mrgsolve_model.R")

# all 18 scenarios (helpers are documented at the bottom of the model file)
df <- SCI_simulate_scenarios(365)

# the two window experiments
SCI_mp_window(ais = 3)            # methylprednisolone start time, AIS C
SCI_decompression_sweep(ais = 3)  # time to surgical decompression, AIS C
SCI_ad_challenge()                # autonomic dysreflexia, T4, by trigger day
```

```r
# Shiny dashboard (12 tabs; model file must be in the same directory)
shiny::runApp("sci_shiny_app.R")
```

```bash
# dependency-free reference implementation — reproduces every number above
python3 sci_reference_model.py          # full report (~3 min)
python3 sci_reference_model.py --quick  # scenarios only
```

Render the map:

```bash
dot -Tsvg sci_qsp_model.dot -o sci_qsp_model.svg
dot -Tpng -Gdpi=150 sci_qsp_model.dot -o sci_qsp_model.png
```

## 9. Scenarios

| # | Scenario | # | Scenario |
|---|---|---|---|
| 1 | Natural history (no decompression, no rehab) | 10 | Glibenclamide 3 d (SUR1–TRPM4) |
| 2 | Standard care (decompression 24 h + rehab) | 11 | MAP ≥ 85 mmHg × 7 d |
| 3 | Early decompression < 12 h | 12 | Bundle: decompression 8 h + MAP + riluzole |
| 4 | Late decompression 72 h | 13 | Anti-Nogo-A + intensive rehab from d14 |
| 5 | Methylprednisolone within 3 h (NASCIS II) | 14 | Intensive rehab EARLY (d14–180) |
| 6 | Methylprednisolone at 9 h (outside window) | 15 | Intensive rehab LATE (d120–286, same dose) |
| 7 | Methylprednisolone 48-h infusion (NASCIS III) | 16 | Symptomatic: baclofen + pregabalin from d30 |
| 8 | Riluzole 14 d (RISCIS) | 17 | Incomplete AIS C, standard care |
| 9 | Minocycline 7 d | 18 | AIS C + epidural stimulation from d120 |

## 10. Shiny dashboard — 12 tabs

Patient & injury profile · Drug PK · Secondary cascade · Neuroinflammation ·
Tissue & imaging surrogates · Neurological recovery · **Therapeutic-window
explorer** · **Surrogate-vs-endpoint explorer** · Secondary complications ·
Scenario comparison · Steroid safety ledger · References.

The two bolded tabs carry the argument: the first sweeps intervention delay and
plots benefit against the complication index on the same axes; the second marks
where the current patient sits on the Hill curve and shows what a doubling of
their spared drive would be worth in ISNCSCI points.

## 11. Limitations

Read [`sci_references.md` § 16](sci_references.md#16-what-is-not-sourced-here-limitations-to-read-before-reusing)
before reusing anything here. In short: most rate constants are shape-matched to
described time courses rather than fitted; `TH50`/`HILL` are a calibrated
modelling choice, so the *direction* of the surrogate–endpoint argument is robust
but the specific point estimates are not; the lesion index is a bounded 0–1 index
with a nominal mL conversion and must not be compared numerically to MRI
volumetry; drug PK is dose-proportional rather than physiological; the
complication index orders regimens by exposure and is not a probability;
`LEVEL_IDX` is a lookup, not a segmental cord, so central-cord and Brown-Séquard
syndromes cannot be represented; and syringomyelia, pressure injury, VTE,
heterotopic ossification and the cardiometabolic trajectory appear on the map but
have no ODEs.

**Research and education only. Not a substitute for clinical judgment, and not
suitable for regulatory use.**
