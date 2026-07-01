# Opioid Use Disorder (OUD, 오피오이드 사용장애) — QSP Model

> Integrated Quantitative Systems Pharmacology model of opioid use disorder,
> linking competitive mu-opioid-receptor (MOR) occupancy from illicit
> fentanyl/heroin, methadone, buprenorphine, naloxone, and naltrexone to the
> mesolimbic reward pathway, chronic tolerance/opponent-process withdrawal
> (locus-coeruleus noradrenergic rebound), craving/relapse risk, dose-
> dependent respiratory depression/overdose, methadone-associated QTc
> prolongation, and treatment-retention outcomes, together with symptomatic
> alpha2-agonist (lofexidine) withdrawal management and harm-reduction
> comorbidity pathways.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`oud_qsp_model.dot`](oud_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`oud_qsp_model.svg`](oud_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`oud_qsp_model.png`](oud_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`oud_mrgsolve_model.R`](oud_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`oud_shiny_app.R`](oud_shiny_app.R) |
| 📚 References             | [`oud_references.md`](oud_references.md) |

---

## 1. Disease in one paragraph

Opioid use disorder is a chronic relapsing brain disease driven by repeated
mu-opioid-receptor (MOR) agonism in the mesolimbic reward pathway: opioid
binding disinhibits VTA dopaminergic neurons (via MOR-mediated suppression
of tonic GABAergic interneurons), triggering nucleus-accumbens dopamine
release, euphoria, and incentive-salience learning that consolidates
drug-cue associations (Volkow 2016 NEJM). With repeated exposure, homeostatic
counter-adaptation (cAMP superactivation, receptor downregulation) produces
**tolerance**, while the locus coeruleus — acutely suppressed by opioid
MOR-Gi signaling — rebounds into noradrenergic hyperactivity when drug
levels fall, driving the autonomic **withdrawal syndrome** (mydriasis,
piloerection, GI cramping, myalgia, anxiety — captured by the Clinical
Opiate Withdrawal Scale, COWS) and, together with dynorphin/kappa-receptor-
mediated anti-reward dysphoria, intense **craving** and relapse risk (Koob &
Volkow 2016 Lancet Psychiatry). The most acute danger is **dose-dependent
respiratory depression**: at brainstem (pre-Bötzinger complex) MOR
saturation, full agonists (illicit fentanyl — now often adulterated with
xylazine — heroin, methadone) can produce fatal hypoventilation, especially
with benzodiazepine co-use; naloxone competitively displaces agonist to
acutely reverse this, though its short half-life creates renarcotization
risk against longer-acting opioids. Evidence-based pharmacotherapy —
**methadone** (full agonist, QTc-prolonging), **buprenorphine** (high-
affinity partial agonist with a **ceiling effect** that caps respiratory
depression, risking precipitated withdrawal if induced too early), and
**naltrexone** (full antagonist, oral or monthly extended-release depot,
which also creates overdose risk if tolerance is lost during treatment
lapses) — stabilizes MOR occupancy, suppresses withdrawal/craving, and
roughly halves mortality versus no opioid-agonist treatment (Sordo 2017
BMJ), while lofexidine offers non-opioid symptomatic withdrawal relief and
counseling/contingency management improves treatment retention.

## 2. Mechanistic clusters (13 in the DOT map, 130 nodes)

1. Etiology & risk factors (genetic OPRM1 variant, chronic pain/iatrogenic exposure, social/adolescent vulnerability)
2. Opioid receptor pharmacology & signaling (MOR/KOR/DOR, Gi/Go, cAMP, GIRK, β-arrestin, analgesia/GI/miosis)
3. Mesolimbic reward pathway (VTA GABA/DA neurons, NAc dopamine release, euphoria, incentive salience, cue memory)
4. Neuroadaptation — tolerance & withdrawal substrate (cAMP superactivation, CREB, locus coeruleus, dynorphin/KOR anti-reward, HPA/CRF)
5. Withdrawal clinical syndrome (autonomic hyperactivity, GI cramping, myalgia, anxiety/insomnia, COWS composite)
6. Respiratory depression & overdose cascade (brainstem MOR, pre-Bötzinger complex, hypoventilation → apnea → death)
7. Drug PK/PD — illicit fentanyl/heroin (full agonist, high potency, xylazine/carfentanil adulteration)
8. Drug PK/PD — methadone (full agonist, long half-life, hERG/QTc, NMDA antagonism)
9. Drug PK/PD — buprenorphine ±naloxone (partial agonist, ceiling effect, precipitated-withdrawal risk)
10. Drug PK/PD — naloxone (acute reversal) & naltrexone (sustained blockade, XR depot)
11. Drug PK/PD — lofexidine/clonidine (α2-agonist, symptomatic withdrawal management)
12. Craving, relapse & behavioral/psychosocial therapy (cue reactivity, stress, contingency management, retention)
13. Comorbidity, harm reduction & biomarkers (IV-use infection risk, benzodiazepine synergy, naloxone distribution, mortality)

## 3. mrgsolve model (23 ODE compartments)

* **Drug PK (6 drug classes, 15 compartments)** — illicit fentanyl (gut/IN +
  central + peripheral), methadone (gut + central), buprenorphine (SL depot
  + XR-SC depot + central), naloxone (IN/IM depot + central), naltrexone
  (oral + XR-IM depot + central), lofexidine (gut + central).
* **MOR occupancy** — computed algebraically each step from a competitive
  multi-ligand Emax equation (relative Ki affinities, ligand-specific
  intrinsic efficacy), with an enforced partial-agonist ceiling for
  buprenorphine-dominated regimens.
* **Disease/PD network (8 compartments)** — MOR tolerance/adaptation index,
  locus-coeruleus noradrenergic tone, withdrawal (COWS 0-48), craving
  (VAS 0-100), respiratory-drive index (with an overdose/apnea threshold),
  methadone-associated QTc, phasic dopamine/euphoria, and treatment-
  retention index.

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated chronic illicit fentanyl use (IN, q4-6h) | Ciccarone 2021 Curr Opin Psychiatry |
| 2 | Fentanyl overdose (IV bolus) + naloxone rescue at 5 min | Boyer 2012 NEJM; Purssell 2021 Sci Rep |
| 3 | Methadone maintenance treatment (100 mg PO QD) | Mattick 2009 Cochrane; Eap 2002 Clin Pharmacokinet |
| 4 | Buprenorphine/naloxone SL maintenance (16 mg PO QD) | Mattick 2014 Cochrane; Walsh 1994 Clin Pharmacol Ther |
| 5 | Precipitated withdrawal (buprenorphine induced too early) | Rosado 2007 Drug Alcohol Depend |
| 6 | Naltrexone XR-IM depot (380 mg, post-detox relapse prevention) | Krupitsky 2011 Lancet; Lee 2018 Lancet (X:BOT) |
| 7 | Untreated abrupt cessation ("cold turkey") + lofexidine | Gowing 2016 Cochrane; Yu 2008 Drug Alcohol Depend |
| 8 | Illicit fentanyl + benzodiazepine co-use (synergistic depression) | Dahan 2010 Anesthesiology |
| 9 | High-dose methadone (140 mg QD) with QTc-risk comorbidity | Krantz 2009 Ann Intern Med; Chou 2014 J Pain |
| 10 | Relapse after buprenorphine discontinuation (tolerance lost) | Bentzley 2015 J Subst Abuse Treat |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — adjustable age/weight/use-history & baseline severity.
2. **Drug PK** — log-scale plasma concentrations for the six tracked drug moieties.
3. **MOR occupancy/tolerance** — total receptor occupancy, tolerance index, phasic dopamine/euphoria.
4. **Withdrawal & craving** — locus-coeruleus tone, COWS score, craving VAS.
5. **Respiratory/overdose risk** — respiratory-drive index vs. the apnea threshold, overdose-flag events.
6. **Clinical endpoints** — methadone QTc and treatment-retention index.
7. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
8. **References** — key citations and link to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg oud_qsp_model.dot -o oud_qsp_model.svg
dot -Tpng -Gdpi=150 oud_qsp_model.dot -o oud_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("oud_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=16, cmt="GUT_BUP", ii=24, addl=14), end=336, delta=0.1)
plot(out, c("MOR_occupancy_total","COWS_score","Respiratory_idx","Craving_VAS"))

# 3) Launch the dashboard
shiny::runApp("oud_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| MOR occupancy / respiratory suppression | Buprenorphine vs. fentanyl/methadone | Partial-agonist ceiling plateaus at ~40-50% max suppression (Walsh 1994; Dahan 2005) |
| Naloxone reversal duration | vs. fentanyl/methadone half-life | Short t1/2 (30-90 min) risks renarcotization against longer-acting agonists |
| Methadone QTc | Dose-dependent hERG block | Risk rises sharply above ~120-150 mg/day (Krantz 2009; Chou 2014) |
| Naltrexone XR blockade | Single 380 mg depot injection | Near-complete MOR blockade sustained ~28 days (Comer 2006; Krupitsky 2011) |
| Treatment retention & mortality | Agonist/antagonist maintenance vs. none | ~50% mortality reduction (Sordo 2017 BMJ; Mattick 2009/2014 Cochrane) |
| Precipitated withdrawal | Early buprenorphine induction | Sharp COWS spike if induced on high residual full-agonist occupancy (Rosado 2007) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support, addiction-medicine guidance, or a substitute
  for validated overdose-risk prediction tools.
* MOR occupancy is modeled via a simplified competitive multi-ligand Emax
  equation with literature-informed relative affinity/efficacy surrogates,
  not a mechanistic receptor-binding kinetic (association/dissociation
  rate) model; real pharmacology includes slower buprenorphine receptor
  dissociation kinetics not explicitly represented here.
* The respiratory-depression Hill function and apnea threshold are
  illustrative surrogates for hypoxic-ventilatory-response suppression,
  not validated overdose-risk predictors — inter-individual variability
  (opioid-naive vs. tolerant, renal/hepatic impairment, genetic CYP
  variability) needs `omega()` blocks for population simulations.
* Xylazine and other non-opioid adulterant effects are represented as a
  simple additive respiratory-suppression term; they are not
  naloxone-reversible in reality, which the model reflects only as an
  unmodified additive offset rather than a full non-opioid sedation
  pathway.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
