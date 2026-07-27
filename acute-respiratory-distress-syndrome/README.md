# Acute Respiratory Distress Syndrome (ARDS) — QSP Model
# 급성 호흡곤란 증후군 QSP 모델

> Quantitative Systems Pharmacology model of ARDS in which **the treatment is
> part of the disease**. A direct or indirect insult breaks both sides of the
> alveolar–capillary barrier; protein-rich oedema floods and de-recruits alveoli;
> the aerated lung shrinks to a **"baby lung"** of 20–35% of the parenchyma — and
> every subsequent breath the ventilator delivers goes into *that* lung. The same
> tidal volume that is harmless in a normal lung becomes a strain of 1.5–2.0, and
> the resulting **ventilator-induced lung injury feeds back into the inflammation
> as biotrauma**. Closing that loop is the point of this model, and it is why the
> interventions with the largest mortality effects here (tidal volume, prone
> position, ECMO) change nothing about the biology at all.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ards_qsp_model.dot`](ards_qsp_model.dot) |
| 🖼️ Map (SVG, full detail) | [`ards_qsp_model.svg`](ards_qsp_model.svg) |
| 🖼️ Map (PNG)              | [`ards_qsp_model.png`](ards_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ards_mrgsolve_model.R`](ards_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ards_shiny_app.R`](ards_shiny_app.R) |
| 📚 References             | [`ards_references.md`](ards_references.md) |

**Scale:** 221 nodes · 17 clusters · 428 edges in the map · 38 ODE compartments ·
14 prebuilt scenarios (19 arms validated) · 8 dashboard tabs · 120 references.

---

## 1. Disease in one paragraph

ARDS is not a disease of one organ compartment but a **barrier failure with a
mechanical consequence**. Pneumonia, aspiration, contusion (direct, epithelial-
first) or extrapulmonary sepsis, pancreatitis, trauma, TRALI (indirect,
endothelial-first) release PAMPs and DAMPs that drive TLR4/MyD88→NF-κB and the
NLRP3 inflammasome. IL-8 pulls neutrophils into the alveolus, where elastase,
MPO, NETs and oxidants cleave the tight junctions of the epithelium and the
VE-cadherin junctions of the endothelium; Ang-2 released from Weibel–Palade
bodies antagonises the Tie2 quiescence signal and the glycocalyx is shed. The
barrier becomes permeable, so protein-rich fluid floods the alveolus — and it
stays there, because alveolar fluid clearance depends on the very ENaC and
Na⁺/K⁺-ATPase of the AT2 cells that are being destroyed (AFC <3%/h predicts
death). Flooded and surfactant-depleted alveoli collapse; the aerated
parenchyma shrinks; shunt through that non-aerated lung produces refractory
hypoxaemia while capillary microthrombosis, driven by tissue factor and PAI-1,
produces dead space. Then the ventilator arrives: the same set tidal volume now
lands in a fraction of the lung, driving pressure and mechanical power rise for
free, and volu-/atelec-/ergotrauma generate biotrauma that re-enters NF-κB.
After roughly a week the disease forks — αvβ6-activated TGF-β1 and collagen
toward fibroproliferation and ventilator dependence, or efferocytosis by M2
macrophages, pro-resolving lipid mediators and AT2→AT1 transdifferentiation
toward resolution.

한국어 요약: ARDS는 폐포 상피와 폐 모세혈관 내피가 동시에 무너져 단백질이 풍부한
폐부종이 폐포를 채우고, 남은 환기 가능한 폐("baby lung")가 전체의 20–35%로 줄어드는
질환입니다. 여기에 인공호흡기가 같은 일회호흡량을 계속 넣으면 구동압·변형·기계적
파워가 저절로 올라가 인공호흡기 유발 폐손상(VILI)이 발생하고, 그 기계적 에너지가
biotrauma로 다시 염증 회로에 되먹임됩니다. 이 되먹임 고리를 닫는 것이 본 모델의
핵심이며, 실제로 사망률을 가장 크게 낮춘 중재(일회호흡량 감소, 복와위, ECMO)가
모두 "약이 아닌" 이유이기도 합니다.

---

## 2. Mechanistic map — 17 clusters

[![ARDS QSP map](ards_qsp_model.png)](ards_qsp_model.svg)

| # | Cluster | What it carries |
|---|---|---|
| 1 | Precipitating insult | direct vs indirect aetiologies, PAMPs/DAMPs, host susceptibility |
| 2 | Alveolar epithelium | AT1/AT2, tight junctions, ENaC/Na-K-ATPase, AFC, surfactant, SP-D/RAGE |
| 3 | Endothelium & glycocalyx | VE-cadherin, RhoA/Rac1-S1P, Ang-1/Tie2 vs Ang-2, syndecan-1 shedding |
| 4 | Innate immunity | TLR4/MyD88, NLRP3-caspase-1-GSDMD, neutrophil influx, NETs, ROS, M1/M2 |
| 5 | Cytokine network | TNF-α, IL-1β, IL-6, IL-8, IL-10, IL-1RA, JAK-STAT hub |
| 6 | Alveolar coagulopathy | tissue factor, thrombin, fibrin, PAI-1, microthrombi, dead space |
| 7 | Alveolar flooding | Starling forces, Pcap, EVLW, BAL/plasma protein ratio, lymphatics |
| 8 | Respiratory mechanics | baby lung, compliance, **driving pressure**, strain, **mechanical power**, P-SILI |
| 9 | VILI | volu-/baro-/atelec-/ergotrauma → mechanotransduction → **biotrauma** |
| 10 | Gas exchange | shunt, HPV, PaO2/FiO2, dead-space fraction, oxygenation index |
| 11 | Pulmonary circulation & RV | PVR, acute cor pulmonale, RV dilatation, cardiac output |
| 12 | Late-phase fork | TGF-β/collagen fibroproliferation **vs** efferocytosis/SPM/Treg resolution |
| 13 | Organ crosstalk | systemic spillover, AKI, ICU-acquired weakness, secondary infection, SOFA |
| 14 | Subphenotypes | hyper- vs hypo-inflammatory latent classes and treatment-effect heterogeneity |
| 15 | Pharmacotherapy | corticosteroids, NMBA, iNO, tocilizumab, baricitinib, diuretics, failed agents |
| 16 | Non-pharmacological | LTVV, PEEP, recruitment, prone, conservative fluids, VV-ECMO, ABCDEF |
| 17 | Clinical endpoints | Berlin class, Murray score, ventilator-free days, mortality, PICS |

Failed and harmful interventions are drawn in the map in red (sivelestat,
salbutamol, adult surfactant, aggressive recruitment manoeuvres) because a model
that only contains what worked cannot explain why so much of ARDS pharmacology
did not.

---

## 3. The mrgsolve model

`ards_mrgsolve_model.R` — 38 ODE compartments, time in **hours**, default
horizon **672 h (28 days)**. Developed and validated under mrgsolve 2.0.1.

| Block | Compartments |
|---|---|
| PK / effect (8) | `DEX`, `MPRED`, `GRE` (glucocorticoid genomic effect), `CIS`, `NMB`, `VDIL` (inhaled vasodilator effect site), `TCZ`, `FUR` |
| Injury cascade (12) | `PATH`, `DAMP`, `NFKB`, `TNF`, `IL1B`, `IL6`, `IL8`, `IL10`, `NEUTA`, `NETS`, `ROS`, plus `M1`/`M2` |
| Barrier (7) | `EPI`, `AT2P`, `SURF`, `ENDO`, `ANG2`, `PAI1`, `FIBA` |
| Physiology (3) | `EVLW`, `FLUID`, `VILI` |
| Late phase (2) | `TGFB`, `COLL` |
| Systemic (5) | `SOFA`, `HAZ`, `WEAK`, `GLU`, `VFD` |

The algebraic respiratory-physiology layer (`ARDS_PHYSIO()`) is expanded
verbatim inside both `[ODE]` and `[TABLE]`, so the mechanics the solver
integrates against and the mechanics that are reported are the same code by
construction. It derives, in order: lung water + surfactant + PEEP → aerated
fraction → compliance → driving pressure → spontaneous effort (P-SILI) →
effective tidal volume → strain and mechanical power → shunt and dead space →
PaO₂/FiO₂, PaCO₂, pH → PVR and cardiac output → alveolar fluid clearance and
transvascular filtration.

### Self-calibrating baseline

A clinician does not know a patient's endothelial integrity; they know the
PaO₂/FiO₂. `[MAIN]` therefore **inverts the physiology**:

1. severity from the presenting PaO₂/FiO₂ → the latent-class cytokine profile;
2. that profile sets how badly hypoxic pulmonary vasoconstriction is blunted;
3. the shunt that reproduces PF₀ is solved by fixed-point iteration (the
   diffusion penalty depends on lung water, which depends on the answer);
4. the non-aerated lung is split into an atelectatic share — capped by what a
   surfactant pool of zero can physically produce at that PEEP — and a flooded
   remainder, which fixes EVLW;
5. a second fixed point solves the Starling balance backwards for the
   epithelial and endothelial integrities that hold that much lung water, split
   between the two barriers by the direct/indirect injury pattern;
6. cytokine production constants are back-calculated so the presenting profile
   is stationary for the NF-κB activity the DAMP load supports.

Change `PF0` from 250 to 90 and every downstream state moves the way it moves in
a real patient, with no initial condition set by hand. Verified across 120
parameter combinations (PF₀ 80–280 × both latent classes × direct/indirect ×
Vt 4–12 mL/kg): no non-finite states, and the realised t=0 PaO₂/FiO₂ tracks the
requested value throughout.

---

## 4. Validation against the trial literature

All rows: reference patient PF₀ 130, PBW 60 kg, hypoinflammatory, pneumonia
(DIRECT 0.8), 672 h, real mrgsolve runs.

| Scenario | Model 28-d mortality | Model VFD | Trial comparator |
|---|---|---|---|
| 01 conventional Vt 12 mL/kg | **38.3%** | 10.9 | ARMA 39.8% / 10 d |
| 02 LTVV Vt 6 mL/kg | **28.8%** | 13.6 | ARMA 31.0% / 12 d |
| 03 LTVV + PEEP 16 (poorly recruitable) | 32.2% | 12.4 | ALVEOLI/ART — no benefit, harm when unrecruitable |
| 04 LTVV + prone 16 h/day | 22.9% | 16.2 | PROSEVA direction |
| 05 prone + cisatracurium, deep sedation | 23.2% | 15.9 | ACURASYS — small increment on top of prone |
| 05b cisatracurium, light sedation | 28.4% | 13.7 | ROSE — neutral |
| 06 dexamethasone 20→10 mg × 10 d | 23.9% | 16.0 | DEXA-ARDS 21% vs 36%, VFD +4.8 |
| 07 conservative fluid + furosemide | 24.1% | 16.5 | FACTT — VFD +2.5, no mortality effect |
| 08 inhaled NO 20 ppm | 26.8% | 14.4 | Taylor 2004 — oxygenation up, mortality unchanged |
| 09 VV-ECMO | 22.1% | 16.2 | EOLIA direction |
| 10 COVID-ARDS: dexamethasone + tocilizumab | 38.1% | 8.5 | RECOVERY / REMAP-CAP |
| 11 hyperinflammatory + steroid | 41.2% | 7.7 | Calfee HTE — steroid ARR **9.5%** in hyper … |
| 11b hyperinflammatory, no steroid | 50.7% | 5.2 | … |
| 12 hypoinflammatory + steroid | 23.9% | 16.0 | … **vs 4.9%** in hypo |
| 13 methylprednisolone started day 14 | 29.3% | 13.4 | LaSRS — no benefit; weakness 0.84 vs 0.38 |
| 14 full bundle | 17.8% | 19.4 | — |
| severe (PF₀ 90): supine / prone / ECMO | 31.9% / 24.7% / 23.6% | 12.5 / 15.4 / 15.6 | PROSEVA 32.8→16.0%; EOLIA 46→35% |

**Behaviour the model reproduces that it was not fitted to.** Inhaled NO raises
PaO₂/FiO₂ by ~20% at 24 h (168→203) while moving mortality by 2 percentage
points — the shunt term improves, the mechanical exposure does not, and the
hazard is dominated by mechanics, dead space and organ failure. Higher PEEP
*harms* the poorly recruitable patient and helps the recruitable one, because
`RECRUIT` enters the atelectasis term while PEEP enters plateau pressure, RV
afterload and mechanical power regardless. Cisatracurium is worth something
under deep sedation and nothing under light sedation, because its only lung
effect runs through the P-SILI term. Corticosteroid benefit is nowhere
hard-coded: it emerges from NF-κB transrepression, reduced neutrophil
recruitment, increased alveolar fluid clearance and reduced TGF-β — and the
late-start harm emerges from secondary-infection seeding and weakness acting on
a patient whose collagen is already deposited, not from a "late" rule.

**Where the model misses, stated plainly.** The mortality effects of prone
position (model ARR 7.2% vs PROSEVA 16.8%) and of dexamethasone (model ARR 4.9%
vs DEXA-ARDS 15%) are roughly **half** the size reported in those two trials.
Both are single trials with effect sizes at the large end of the ARDS
literature, and reproducing them exactly would require attributing more of the
mortality signal to regional stress homogenisation, and to NF-κB
transrepression, than the physiological evidence supports. Direction, ordering
and mechanism are reproduced; those two magnitudes are not tuned to match and
should not be read as predictions of trial results.

---

## 5. Modelling choices worth arguing about

- **Three partly redundant mechanical exposures.** VILI is driven by strain,
  driving pressure *and* mechanical power, because no one of them explains the
  data alone: ARMA separates on tidal volume, Amato's mediation analysis on
  driving pressure, the observational literature on J/min. A model with only one
  cannot reproduce all three findings.
- **Two feedbacks in the fluid block, without which it is meaningless.**
  Pressure natriuresis (urine output rises with positive balance) and a volume
  gate on the diuretic response. Without them a fixed intake or a standing
  furosemide order integrates for 28 days into a balance no patient ever has,
  and the capillary-pressure term then dominates the lung water. With them the
  liberal arm sits at +4.5 L at day 7 (FACTT liberal: +6.9 L) and the
  conservative arm at −1.8 L (FACTT: ≈ −0.1 L).
- **Saturating injury.** Every injury intensity is Michaelis–Menten in its
  driver. Without saturation the epithelium→DAMP→NF-κB→IL-8→neutrophil→NET→
  epithelium loop has a gain above 1 and the model runs away — which is a real
  property of the biology being modelled, not a numerical artefact, and is worth
  noticing: ARDS resolves only because repair outcompetes a *bounded* injury.
- **Weaning readiness is judged at a weaning setting.** The gate uses the
  driving pressure a 6 mL/kg breath *would* generate against the current
  compliance, not the driving pressure of whatever the patient is on. Otherwise
  a patient left on 12 mL/kg is scored as unweanable forever.
- **On ECMO, decannulation depends on the native lung.** The reported arterial
  PaO₂/FiO₂ includes the circuit's work, so the weaning gate uses the
  ventilator-only oxygenation.
- **IL-6 measured ≠ IL-6 signalling.** Tocilizumab raises the assay while
  lowering the biology; the model reports both, deliberately.

---

## 6. Running it

```r
# ODE model
library(mrgsolve)
mod <- mread("ards_mrgsolve_model.R")          # from this directory

# reference patient, lung-protective ventilation
sim <- mrgsim(param(mod, PF0 = 130, VT_KG = 6, PEEP = 12, RR = 28),
              end = 672, delta = 6)

# add the DEXA-ARDS dexamethasone schedule
dexa <- c(ev(amt = 20, cmt = "DEX", ii = 24, addl = 4, time = 0),
          ev(amt = 10, cmt = "DEX", ii = 24, addl = 4, time = 120))
sim2 <- mrgsim(param(mod, PF0 = 130, VT_KG = 6, PEEP = 12, RR = 28),
               events = dexa, end = 672, delta = 6)
```

```r
# interactive dashboard (8 tabs)
shiny::runApp("ards_shiny_app.R")
```

Dose-unit conventions and all 14 prebuilt scenarios are documented at the bottom
of `ards_mrgsolve_model.R`. Non-pharmacological interventions are *parameters*,
not doses, because that is what they physically are — a change of ventilator
setting or of body position:

| Intervention | Parameters |
|---|---|
| Lung-protective ventilation | `VT_KG = 6, PEEP = 12, RR = 28` |
| Conventional 1990s | `VT_KG = 12, PEEP = 8, RR = 16` |
| Prone (PROSEVA) | `PRONE_H = 16` |
| Deep vs light sedation | `SEDATION = 0.9` vs `0.35` |
| Conservative fluids | `FLUID_IN = 60` + furosemide events |
| VV-ECMO | `ECMO = 1` (forces `VT_ECMO`, `RR_ECMO`) |
| Inhaled NO / epoprostenol | `NO_PPM = 20` |
| Hyperinflammatory class | `HYPER = 1` |

---

## 7. Shiny dashboard — 8 tabs

1. **Patient & phenotype** — the solved presenting state, Berlin class, barrier trajectories
2. **Ventilator & mechanics** — compliance, driving pressure, plateau, strain, mechanical power, the VILI loop
3. **Gas exchange** — PaO₂/FiO₂ (arterial and ventilator-only), shunt, dead space, PaCO₂, pH
4. **Lung water & barrier** — EVLW, permeability, alveolar fluid clearance, fluid balance
5. **Inflammation** — cytokines (log scale), neutrophils/NETs/ROS, M1↔M2, Ang-2, PAI-1
6. **Drug PK/PD** — glucocorticoid effect, NMB depth vs effort, tocilizumab, glucose and weakness
7. **Clinical endpoints** — mortality, ventilator-free days, SOFA, Murray score, right heart
8. **Scenario comparison** — up to ten arms in the *same* virtual patient, with an endpoint table

---

## 8. Limitations

- Ventilator settings are held fixed unless changed by the user; the model does
  not implement an automatic PEEP/FiO₂ table or a weaning protocol, so the
  recovery-phase settings of a real ICU are not simulated.
- Gas exchange uses a calibrated empirical shunt→PaO₂/FiO₂ mapping rather than a
  full oxyhaemoglobin dissociation inversion. It is anchored to Riley shunt
  analyses and reproduces the Berlin thresholds, but it is not a blood-gas
  simulator and should not be used to predict an individual arterial gas.
- Single virtual patient: no population variability, no inter-individual random
  effects, no covariate model. The library's models are deterministic by design.
- The hazard model is a proportional accrual of risk factors calibrated to
  trial-level mortality, not a validated individual prediction score. It must
  not be used for prognostication.
- Fibroproliferation is represented by a single collagen pool; the model does
  not distinguish the radiological patterns or the long-term DLCO trajectory.

> ⚠️ **Educational and research use only.** This model has not been
> independently validated or certified and must not be used for clinical
> decision-making, prescribing, or regulatory submission.
