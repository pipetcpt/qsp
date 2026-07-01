# Ménière's Disease (메니에르병) — QSP Model

> Integrated Quantitative Systems Pharmacology model of Ménière's disease,
> linking endolymphatic sac absorptive dysfunction and stria vascularis
> ion-transport imbalance (Na+/K+-ATPase, NKCC1, aquaporin/vasopressin-V2
> axis) to endolymphatic hydrops (cochlear duct, saccule and utricle
> distension), episodic Reissner's-membrane micro-rupture with perilymph
> K+ intoxication, and progressive cochlear (outer/inner hair cell) and
> vestibular hair-cell injury — culminating in episodic vertigo attacks,
> fluctuating low-to-pan-frequency sensorineural hearing loss, tinnitus,
> aural fullness and (in advanced disease) Tumarkin otolithic drop attacks.
> The model spans the modern pharmacology stack (betahistine, thiazide/K+-
> sparing diuretics, intratympanic dexamethasone, intratympanic gentamicin)
> together with surgical/procedural options and AAO-HNS 1995 diagnostic
> staging.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`med_qsp_model.dot`](med_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`med_qsp_model.svg`](med_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`med_qsp_model.png`](med_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`med_mrgsolve_model.R`](med_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`med_shiny_app.R`](med_shiny_app.R) |
| 📚 References             | [`med_references.md`](med_references.md) |

---

## 1. Disease in one paragraph

Ménière's disease is a chronic inner-ear disorder characterized pathologically
by endolymphatic hydrops — distension of the endolymph-containing cochlear
duct, saccule and (to a lesser degree) utricle — thought to arise from
endolymphatic sac absorptive failure and/or endolymph overproduction by the
stria vascularis, with contributions proposed from vasopressin-V2-receptor
signaling, autoimmune/anti-inner-ear-antibody mechanisms, viral triggers,
genetic susceptibility (DTNA, FAM136A, PRKCB familial loci) and vascular
factors (Sajjadi & Paparella 2008 Lancet). As hydrops progresses, episodic
mechanical/ionic stress causes micro-ruptures of Reissner's membrane, allowing
K+-rich endolymph to leak into perilymph and transiently depolarize/block
cochlear and vestibular hair cells and afferents — the presumed substrate of
the acute rotational vertigo attack (20 minutes to 12 hours, per Bárány
Society/AAO-HNS criteria). Repeated injury cycles progressively damage outer
and inner hair cells (fluctuating, eventually permanent, low-frequency-first
sensorineural hearing loss and tinnitus) and vestibular hair cells (which,
paradoxically, can eventually reduce attack frequency as the labyrinth
"burns out," while occasionally producing sudden vestibular/otolithic Tumarkin
drop attacks). Diagnosis follows the 1995 AAO-HNS Committee on Hearing and
Equilibrium criteria (updated 2015 by the Bárány Society/AAO-HNS/EAONO/
Japan/Korea consortium) and disease is staged I-IV by four-tone pure-tone
average. First-line therapy is dietary/lifestyle modification plus
betahistine and/or diuretics (modest, heterogeneous trial evidence — the
BEMED trial found no superiority of high-dose betahistine over placebo);
second-line is intratympanic corticosteroid; refractory unilateral disease
is treated with intratympanic gentamicin (chemical labyrinthectomy, dose-
titrated to reduce hearing-loss risk) or, rarely, endolymphatic sac surgery,
vestibular neurectomy or surgical labyrinthectomy.

## 2. Mechanistic clusters (18 in the DOT map, 150 labeled nodes)

1. Etiology, genetic & autoimmune susceptibility (DTNA/FAM136A loci, anti-inner-ear antibodies, migraine/vascular theory)
2. Endolymphatic sac dysfunction & fluid homeostasis (AQP2/AQP3, vasopressin-V2-cAMP axis)
3. Stria vascularis ion transport & endolymph composition (Na+/K+-ATPase, NKCC1, KCNQ1/KCNE1, connexin gap junctions)
4. Endolymphatic hydrops formation (cochlear duct/saccule/utricle distension, Reissner's membrane rupture)
5. Cochlear hair cell & auditory neural injury (OHC/IHC viability, spiral ganglion neurons, oxidative stress/apoptosis)
6. Vestibular hair cell & end-organ dysfunction (semicircular canal/otolith function, VOR, nystagmus)
7. Acute vertigo attack pathophysiology (pressure spikes, Tumarkin otolithic crisis, fall risk)
8. Clinical manifestations, diagnostic criteria & staging (AAO-HNS 1995, FLS, burnout phase, bilateral progression)
9. Patient-reported & composite endpoints (DHI, THI, ECoG SP/AP ratio, VEMP)
10. Drug PK/PD — betahistine (H1 agonism/H3 antagonism, active metabolite, vestibular compensation)
11. Drug PK/PD — diuretics (HCTZ/triamterene, acetazolamide/carbonic anhydrase inhibition)
12. Drug PK/PD — intratympanic dexamethasone (round-window permeation, GR genomic signaling)
13. Drug PK/PD — intratympanic gentamicin (selective vestibular ablation, cochleotoxicity risk)
14. Acute symptomatic pharmacotherapy (meclizine, benzodiazepines, antiemetics)
15. Surgical & procedural therapy (endolymphatic sac decompression/shunt, neurectomy, labyrinthectomy, Meniett device)
16. Investigational & emerging therapy (V2-receptor antagonists, intratympanic IGF-1, hair-cell regeneration)
17. Systemic/renal effects & central adaptation (electrolyte handling, vestibular compensation, psychological comorbidity)
18. Adverse effects & safety monitoring (iatrogenic hearing loss, bilateral vestibular loss, surgical complications)

## 3. mrgsolve model (17 ODE compartments)

* **Drug PK (4 agents, 10 compartments)** — betahistine (oral depot + parent +
  active metabolite, 3-cmpt), diuretic (lumped HCTZ/triamterene signal,
  2-cmpt), intratympanic dexamethasone (middle-ear depot + perilymph +
  genomic GR signal, 3-cmpt), intratympanic gentamicin (middle-ear depot +
  perilymph, 2-cmpt).
* **Inner-ear disease network (7 compartments)** — endolymphatic hydrops
  index, perilymph K+-intoxication index (episodic), cochlear outer-hair-cell
  viability, vestibular hair-cell viability, vertigo attack-frequency state,
  central vestibular compensation index, Dizziness Handicap Inventory.
* **Derived outputs (via $TABLE)** — PTA (dB HL), tinnitus severity, ECoG
  SP/AP ratio, Tumarkin drop-attack risk flag, AAO-HNS stage, QoL score.

### 8 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history (moderate MD, stage II) | Sajjadi & Paparella 2008 Lancet natural history |
| 2 | Betahistine standard-dose (16mg TID) | Nauta 2014 meta-analysis |
| 3 | Betahistine high-dose (48mg TID, BEMED regimen) | Adrion 2016 BMJ (BEMED trial, PMID 27075667) |
| 4 | Diuretic (HCTZ 25mg + triamterene 50mg, once daily) | James & Burton 2001 / Thirlwall & Kundu 2006 Cochrane |
| 5 | Intratympanic dexamethasone series (weekly x4) | Phillips & Westerberg 2011 Cochrane |
| 6 | Low-dose intratympanic gentamicin (titration protocol) | Postema 2008; Boleas-Aguirre 2008; Patel 2016 Lancet |
| 7 | High-dose intratympanic gentamicin (fixed, destructive) | Patel 2016 Lancet (high-dose arm) |
| 8 | Combination: betahistine high-dose + diuretic | Common real-world combination therapy |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — adjustable age/weight/baseline stage/duration/laterality.
2. **Drug PK** — betahistine parent+metabolite, diuretic, perilymph dexamethasone/gentamicin.
3. **Inner-ear PD / Hydrops** — hydrops index, K+ intoxication, hair-cell viability, compensation.
4. **Clinical endpoints** — vertigo frequency, PTA, tinnitus severity, DHI/QoL.
5. **Scenario comparison** — runs all 8 scenarios with the chosen profile.
6. **Biomarkers (ECoG)** — SP/AP ratio, AAO-HNS stage trajectory, Tumarkin risk flag.
7. **Safety** — cumulative gentamicin exposure, cochlear-vs-vestibular injury selectivity.
8. **References** — key citations and link to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg med_qsp_model.dot -o med_qsp_model.svg
dot -Tpng -Gdpi=150 med_qsp_model.dot -o med_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("med_mrgsolve_model.R")
out <- mod %>% ev(amt=48, cmt="GUT_BETA", ii=8, addl=999) %>%
  mrgsim(end=8760, delta=6)
plot(out, c("VertigoFreq_permo","PTA_dB","Hydrops_idx","DHI_score"))

# 3) Launch the dashboard
shiny::runApp("med_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Definitive vertigo-day rate | High-dose betahistine vs placebo | No significant superiority (Adrion 2016 BMJ, BEMED, PMID 27075667) |
| Vertigo control | Low-dose vs high-dose intratympanic gentamicin, 18mo | Equivalent control; low-dose trend to less hearing loss (Patel 2016 Lancet) |
| Vertigo/hearing control | Intratympanic methylprednisolone vs gentamicin | Comparable vertigo control, steroid favored for hearing preservation (Patel 2016 Lancet) |
| Diagnostic staging | 4-tone PTA average | AAO-HNS 1995 stage I (≤25dB) - IV (>70dB) |
| Diuretic benefit | vs placebo/no treatment | Insufficient high-quality RCT evidence despite widespread use (Cochrane 2001/2006) |
| Vertigo control | Intratympanic gentamicin titration protocols | Postema 2008; Boleas-Aguirre 2008 |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* All PK/PD parameters are illustrative literature-informed surrogates
  (Ménière's disease has no validated population PK/PD model in the public
  literature for most of these agents); inter-individual variability needs
  `omega()` blocks for population simulations.
* The episodic, stochastic nature of Reissner's-membrane rupture and vertigo
  attacks is approximated deterministically here (a smooth hydrops-driven
  rupture-rate function); a stochastic/hybrid ODE-Poisson-process
  reformulation would better capture attack unpredictability.
* Betahistine's mechanism (H1 agonism/H3 antagonism, cochlear blood flow,
  vestibular compensation) is biologically plausible but the BEMED trial
  found no significant benefit over placebo on hard endpoints at high dose —
  this is deliberately reflected via small effect-size parameters rather
  than an idealized dose-response.
* Diuretic PK is a lumped single-compartment surrogate for the combined
  HCTZ+triamterene (or acetazolamide) regimen; it does not separately model
  each agent's distinct renal transporter pharmacology.
* Tumarkin drop-attack risk is represented as a threshold flag on hydrops
  and vestibular viability, not a mechanistic falls-biomechanics model.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
