# Peripartum Cardiomyopathy (PPCM) — QSP Model

> Integrated Quantitative Systems Pharmacology model of peripartum
> cardiomyopathy, built around the **prolactin-cleavage / anti-angiogenic**
> mechanism: peripartum haemodynamic and hormonal load on a STAT3/MnSOD-deficient
> myocardium → oxidative stress → **cathepsin D activation** → cleavage of
> 23-kDa prolactin into **16-kDa prolactin (vasoinhibin)** → together with
> **placental sFlt-1** a two-hit anti-angiogenic milieu → coronary microvascular
> capillary rarefaction, plus NF-κB-driven endothelial **miR-146a-5p** delivered
> by exosomes to cardiomyocytes → loss of Erbb4/Nras/Notch1 survival signalling
> → myocyte loss, fibrosis, LV dilation and systolic failure. The pharmacology
> layer covers **bromocriptine** (the disease-specific intervention) and the full
> guideline heart-failure stack under the **antepartum fetotoxicity gate** that
> defines prescribing in this disease.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ppcm_qsp_model.dot`](ppcm_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`ppcm_qsp_model.svg`](ppcm_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ppcm_qsp_model.png`](ppcm_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ppcm_mrgsolve_model.R`](ppcm_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ppcm_shiny_app.R`](ppcm_shiny_app.R) |
| 📚 References             | [`ppcm_references.md`](ppcm_references.md) |

---

## 1. Disease in one paragraph

Peripartum cardiomyopathy is new-onset systolic heart failure (LVEF <45%, no
other identifiable cause) arising in the last month of pregnancy or within about
five months of delivery. It is not simply heart failure that happens to occur
around childbirth: it has a specific, largely reversible pathobiology. In late
gestation the maternal circulation carries a 30–50% higher cardiac output and a
40–50% larger plasma volume, prolactin rises roughly ten-fold, and the placenta
releases the VEGF-trapping decoy receptor **sFlt-1**. In a myocardium whose
antioxidant reserve is compromised — reduced cardiomyocyte **STAT3** and hence
reduced **MnSOD**, frequently on a background of a **TTN truncating variant**
that leaves the sarcomere with less reserve — this load generates oxidative
stress that activates lysosomal **cathepsin D**. Cathepsin D cleaves the
abundant 23-kDa prolactin into a 16-kDa fragment, **vasoinhibin**, which is
strongly anti-angiogenic and pro-apoptotic for endothelium. Acting together with
placental sFlt-1 (the *two-hit* model), it causes coronary microvascular
endothelial apoptosis and **capillary rarefaction**, and it drives endothelial
NF-κB to produce **miR-146a-5p**, which is packaged into exosomes, transferred
to cardiomyocytes, and knocks down **Erbb4, Nras and Notch1** — the myocyte's
own survival and metabolic signalling. The result is myocyte loss, replacement
fibrosis, chamber dilation and systolic failure, amplified by RAAS and
sympathetic activation. The mechanism is also the therapeutic target: the
dopamine D2 agonist **bromocriptine** removes the prolactin substrate, and this
is the only disease-specific therapy PPCM has.

## 2. Mechanistic clusters (15 in the DOT map, 195 nodes, 305 edges)

1. **Genetic substrate and predisposing risk factors** — TTN truncating variant (~15%), BAG3, DSP, MYH7/TNNT2, LMNA, FLNC; shared architecture with idiopathic DCM; African ancestry; multiparity; advanced maternal age; twin gestation; preeclampsia (~22%); prolonged β-agonist tocolysis; anaemia; selenium deficiency; delayed care access
2. **Late-pregnancy and peripartum haemodynamic / hormonal load** — plasma-volume expansion, cardiac-output rise, falling systemic vascular resistance, physiological eccentric hypertrophy, relaxin, oestrogen/progesterone, the delivery event, postpartum uterine autotransfusion, lactation demand, suckling-induced prolactin surge
3. **Prolactin axis and cathepsin D cleavage** — hypothalamic dopamine, lactotroph D2 receptor, 23-kDa prolactin, PRLR/JAK2/STAT5 pro-survival signalling, pro- and active cathepsin D, the Ser–Asn cleavage site, 16-kDa prolactin (vasoinhibin)
4. **Oxidative stress · STAT3 · MnSOD axis** — STAT3-driven MnSOD induction, NOX2/NOX4, peroxynitrite, uncoupled eNOS, Nrf2, lysosomal destabilisation, the mitochondrial ROS feed-forward loop, 8-isoprostane/oxLDL
5. **Anti-angiogenic milieu** — placental sFlt-1, soluble endoglin, low PlGF, cardiac PGC-1α, myocardial VEGF production, free VEGF, VEGFR-2 signalling, the sFlt-1/PlGF ratio, the two-hit node
6. **Coronary microvascular endothelium** — endothelial apoptosis, capillary density, NO bioavailability, endothelin-1, coronary flow reserve, microvascular ischaemia, pericyte detachment, permeability
7. **16K-PRL → endothelial exosome → miR-146a cardiomyocyte axis** — endothelial NF-κB, miR-146a-5p, exosome packaging and transfer, serum miR-146a as a PPCM-specific biomarker, Erbb4/Nras/Notch1/Irak1, neuregulin-1 signalling
8. **Cardiomyocyte injury, energetics and structural remodelling** — mitochondrial dysfunction, ATP and PCr:ATP, fatty-acid oxidation, substrate shift, sarcomere dysfunction, Ca²⁺ handling (SERCA2a, RyR2), apoptosis, hypertrophy, autophagy, viable contractile mass, fibroblast activation, TGF-β1/SMAD, collagen, MMP:TIMP, fibrosis, galectin-3
9. **Immune / inflammatory and viral contribution** — IFN-γ, TNF-α, IL-6, sFas/Fas, hs-CRP, macrophage infiltration, anti-cardiac autoantibodies, fetal microchimerism, parvovirus B19/HHV-6, myocarditis overlap
10. **LV dysfunction, chamber remodelling and haemodynamics** — contractility, LVEF, LVEDD/LVEDV, LVESV, stroke volume, cardiac output, filling pressure, functional mitral regurgitation, RV dysfunction (~30%), congestion, global longitudinal strain, wall stress, reverse remodelling
11. **Neurohormonal activation** — renin–angiotensin–aldosterone, AT1 signalling, sympathetic activation, β1 downregulation, vasopressin, sodium/water retention, cardiorenal interaction, ANP/BNP, neprilysin, NT-proBNP, cGMP counter-regulation
12. **Bromocriptine PK/PD (disease-specific, BOARD concept)** — absorption with high first-pass, CYP3A4 metabolism, pituitary D2 agonism, prolactin suppression, lactation suppression, the 1-week and 8-week regimens, prothrombotic risk mandating anticoagulation, adverse effects
13. **Guideline heart-failure pharmacotherapy** — β-blocker, ACEi/ARB with its **antepartum fetotoxicity block**, ARNI, MRA, SGLT2 inhibitor, loop diuretic, antepartum-safe hydralazine + nitrate, ivabradine, digoxin, anticoagulation, levosimendan and IVIG/pentoxifylline (mapped but unproven), the breastfeeding-compatibility filter
14. **Device, mechanical support and obstetric management** — wearable cardioverter defibrillator, ICD deferral pending recovery, IABP/VA-ECMO, LVAD, transplantation, delivery timing, preterm delivery, CRT, contraception and subsequent-pregnancy counselling
15. **Clinical endpoints, complications and outcomes** — NYHA class, 6-minute walk, LVEF recovery ≥50%, persistent dysfunction, LV thrombus, thromboembolism, ventricular arrhythmia/SCD, cardiogenic shock, rehospitalisation, mortality, MACE, subsequent-pregnancy relapse, quality of life, lactation cessation and maternal–infant bonding, peripartum depression

## 3. The mrgsolve model

**30 ODE compartments** — 11 drug PK/exposure states and 19 disease states.

| Group | Compartments |
|---|---|
| Bromocriptine PK | `BRC_G` (oral depot), `BRC_C` (central), `BRC_P` (peripheral) |
| HF agent exposures | `BB`, `RASI` (ACEi/ARB/ARNI), `MRA`, `DIU`, `SG` (SGLT2i), `VD` (hydralazine+nitrate), `AC` (anticoagulant), `DIG` |
| Prolactin axis | `PRL23`, `CTSD`, `PRL16`, `LACT` |
| Oxidative / angiogenic | `ROS`, `SFLT`, `VEGF`, `CAP` |
| miR-146a axis | `MIR`, `SURV` |
| Myocardium | `CM`, `SCAR`, `FIB`, `INF` |
| Systemic | `NEURO`, `VOL`, `LVEDV`, `NTBNP`, `THR` |

### Three design decisions worth knowing

**The cathepsin D threshold is the disease switch.** 16-kDa prolactin is
generated only when cathepsin D activity exceeds its physiological baseline
(`CTSD_TOL`). This is what makes the model clinically honest: a normal lactating
woman has abundant 23-kDa prolactin and *no* disease, so breastfeeding alone is
not pathogenic here — consistent with the absence of convincing evidence that
breastfeeding worsens PPCM outcome. Disease requires the prolactin substrate
**and** an oxidatively stressed myocardium.

**Recovery has a path-dependent ceiling.** Injury above a physiological
turnover threshold accumulates as irreversible `SCAR`, which caps how far
`CM` can recover. Damage prevented early is therefore worth more than damage
treated late — this is the mechanism by which a time-limited intervention
produces a permanent difference in outcome.

**Lactation involution is absorbing.** `LACT` is an autocatalytic supply-and-demand
state; once bromocriptine drives it below `LACT_MIN` it cannot restart. Milk
supply, once involuted, does not spontaneously return.

### Two hard pharmacological gates

Both are enforced in the model rather than left to the user, because they are
the defining prescribing constraints of PPCM:

- **ACEi/ARB/ARNI effect is zero while `t < T_DELIVERY`** (antepartum
  fetotoxicity: fetal renal failure, oligohydramnios). Dosing records may be
  present; the effect simply does not exist until delivery. Hydralazine +
  nitrate is the antepartum-safe substitute.
- **Bromocriptine effect is zero while `t < T_DELIVERY`** — it is a postpartum
  intervention.

### Verified behaviour

The `[MAIN]`/`[ODE]`/`[TABLE]` blocks were compiled verbatim into a standalone
integrator and run before release. From a presenting LVEF of 29.5%
(`SEV = 0.72`), delivered on day −7, patient breastfeeding:

| Arm | d30 | d90 | **d180** | d365 | Scar at 1 y |
|---|---|---|---|---|---|
| Healthy control (no substrate) | 60.4 | 62.8 | **63.6** | 63.9 | 0.02 |
| No therapy | 25.6 | 23.1 | **20.0** | 6.3 | 0.69 |
| Standard HF therapy only | 34.3 | 39.1 | **42.4** | 38.4 | 0.35 |
| + bromocriptine 1 week | 37.9 | 45.2 | **47.6** | 45.2 | 0.25 |
| + bromocriptine 8 weeks (BOARD) | 38.9 | 53.3 | **56.5** | 56.4 | 0.14 |
| BOARD + SGLT2 inhibitor | 39.7 | 54.1 | **57.1** | 57.0 | 0.13 |
| TTN carrier on full BOARD | 25.2 | 39.6 | **43.8** | 42.4 | 0.29 |
| Antepartum preeclamptic (delivery d14) | 22.6 | 40.3 | **48.0** | 48.4 | 0.22 |

Also verified:

- the **healthy control has a stable normal attractor** (LVEF 64.5% → 63.9% at
  one year, LVEDD 4.6 cm, negligible scar) — the model does not simply decay;
- the **untreated arm reaches a genuine non-recovery attractor** (LVEDD 6.2 cm,
  NT-proBNP ~7400 pg/mL, NYHA 3, scar fraction 0.69);
- the **antepartum arm nadirs at 21.7% *before* delivery**, while ACE inhibition
  is gated off and only hydralazine/nitrate is available, then recovers once the
  placenta and its sFlt-1 are gone;
- **anticoagulation leaves LVEF completely unchanged** and cuts peak LV thrombus
  propensity from 0.40 to 0.21 au — a safety intervention cleanly separated from
  an efficacy one, which is exactly why bromocriptine should not be given
  without it.

### Nine prebuilt scenarios

`PPCM_simulate_scenarios()` runs: (1) no therapy, (2) standard HF therapy,
(3) + bromocriptine 1 week + LMWH, (4) BOARD with bromocriptine 8 weeks,
(5) BOARD + SGLT2 inhibitor, (6) TTN-variant carrier on BOARD, (7) antepartum
preeclamptic presentation with ACE inhibition gated off, and (8)/(9) a matched
severe pair with and without anticoagulation that isolates thrombus safety from
LVEF efficacy.

## 4. The Shiny dashboard

Eight tabs: patient profile · drug PK/exposure · prolactin axis · microvascular
and miR-146a · LV function · clinical endpoints · scenario comparison ·
biomarker panel. The sidebar exposes susceptibility (severity, TTN variant,
preeclampsia, STAT3 deficit), delivery timing, the bromocriptine regimen and
each HF agent independently. When the patient is antepartum the app annotates
the delivery day and reports the ACEi/bromocriptine gates as *dosed but gated
off*, so the constraint is visible rather than implicit.

```r
shiny::runApp("ppcm_shiny_app.R")
```

## 5. Known limitations

- **The bromocriptine duration effect is exaggerated.** The model reproduces the
  *direction* of the German trial's 1-week versus 8-week result and offers a
  mechanism for it (coverage of the early high-oxidative-stress window, with
  early damage locked in as scar), but the mean-LVEF gap it produces (47.6%
  versus 56.5% at 6 months) is larger than the trial's (0.49 versus 0.51, even
  though the *recovery-rate* difference there was 16 points). Read the duration
  comparison as a hypothesis, not a prediction.
- **Bromocriptine efficacy itself is contested** — the BRO-HF retrospective
  cohort did not find an independent benefit. The model encodes the mechanistic
  hypothesis, not a settled fact.
- **Drug exposures are dose-proportional arbitrary units**, not validated plasma
  PK; only bromocriptine has a structured depot-plus-two-compartment model.
- **Fibrosis may be over-weighted as a recovery brake.** Late gadolinium
  enhancement is in fact uncommon in PPCM; the reversible-`FIB`/irreversible-`SCAR`
  split respects that, but the balance between the two is not identified from data.
- **No mortality or competing-risk model.** `MACE_hazard_yr` is an illustrative
  algebraic surrogate; no patient leaves the simulation.
- **Deterministic single-patient simulation** — there is no `$OMEGA`/`$SIGMA`
  variability block, so outputs are typical-patient trajectories. Comparing a
  model LVEF value against a published *recovery percentage* requires that caveat.

## 6. References

63 annotated references in [`ppcm_references.md`](ppcm_references.md), grouped
by the role each plays in the model. The mechanistic backbone rests on
Hilfiker-Kleiner 2007 *Cell* (STAT3/cathepsin D/16K-prolactin), Patten 2012
*Nature* (PGC-1α/VEGF and placental sFlt-1 as the second hit) and Halkein 2013
*J Clin Invest* (the miR-146a exosome axis); the quantitative anchors are IPAC
(McNamara 2015 *JACC*), the German bromocriptine trial (Hilfiker-Kleiner 2017
*Eur Heart J*), Sliwa 2010 *Circulation*, the ESC EORP registry (Sliwa 2020
*Eur Heart J*) and Ware 2016 *NEJM*.

---

*Part of the [QSP Disease Model Library](../README.md). Research, education and
hypothesis generation only — not for clinical decision-making or regulatory use.*
