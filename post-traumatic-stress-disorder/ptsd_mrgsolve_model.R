## =============================================================================
## Post-Traumatic Stress Disorder (PTSD) -- mrgsolve QSP model
##
## Trauma exposure -> amygdala-driven fear acquisition (noradrenergic/
## glucocorticoid memory consolidation) -> HPA axis dysregulation (enhanced
## GR negative feedback, paradoxically low cortisol) + locus coeruleus (LC)
## noradrenergic hyperarousal -> amygdala hyperreactivity / vmPFC
## hypoactivation (fear-extinction circuit failure) -> DSM-5 symptom clusters
## (intrusion, avoidance, negative cognition/mood, hyperarousal) -> CAPS-5
## composite endpoint, coupled to SSRI (sertraline/paroxetine), prazosin
## (alpha-1 blocker, nightmares), ketamine/esketamine (rapid NMDA-antagonist
## extinction-facilitation), and MDMA-assisted / trauma-focused psychotherapy
## (session-dose driven extinction-learning boost) PK/PD.
##
## Time unit: hours. Disease horizon simulated over months-years
## (end = 24*7*N weeks or 24*365*N years).
##
## Calibration anchors (see ptsd_references.md for full PMID list):
##   - HPA axis "paradoxical" low cortisol / enhanced GR feedback: Yehuda 2001
##     Biol Psychiatry; Yehuda 2009 Prog Brain Res
##   - Noradrenergic/glucocorticoid memory consolidation: Cahill & McGaugh 1996
##     Trends Neurosci; Pitman 2012 Nat Rev Neurosci
##   - Amygdala hyperreactivity / vmPFC hypoactivation fear-extinction circuit:
##     Rauch 2006 Biol Psychiatry; Milad 2009 Biol Psychiatry; Shin & Liberzon
##     2010 Neuropsychopharmacology
##   - Hippocampal volume reduction: Bremner 2003 Am J Psychiatry; Woon 2010
##     Prog Neuropsychopharmacol Biol Psychiatry (meta-analysis)
##   - DSM-5 symptom clusters & CAPS-5: Weathers 2018 Psychol Assess
##   - Sertraline/paroxetine RCTs (FDA-approved): Brady 2000 JAMA (sertraline);
##     Marshall 2001 Am J Psychiatry (paroxetine); Stein 2006 Cochrane review
##   - Prazosin for nightmares/sleep: Raskind 2013 Am J Psychiatry; Raskind
##     2018 NEJM (negative confirmatory trial, heterogeneous effect)
##   - Ketamine rapid PTSD symptom reduction: Feder 2014 JAMA Psychiatry;
##     Feder 2021 Am J Psychiatry (repeated-dose)
##   - MDMA-assisted psychotherapy: Mitchell 2021 Nat Med (phase 3);
##     Mitchell 2023 Nat Med (confirmatory phase 3)
##   - Prolonged Exposure / CPT / EMDR efficacy: Foa 2005 JAMA; Resick 2002
##     J Consult Clin Psychol; Watts 2013 J Clin Psychiatry (meta-analysis)
## =============================================================================

$PROB
# PTSD QSP model (24-compartment PK/PD + fear-circuit + symptom-cluster + biomarker)

$PARAM
// ---- Trauma / risk-resilience severity ----
TRAUMA_SEVERITY = 1.0     // 0.6=single mild trauma, 1.0=reference (moderate combat/assault), 1.5=severe/repeated (childhood abuse + adult retraumatization)
FKBP5_RISK      = 0       // 1 = FKBP5 risk-allele x childhood-adversity interaction present (amplifies HPA/GR dysregulation)
RESILIENCE      = 1.0      // composite resilience index (NPY/BDNF/social support); >1 = protective, <1 = vulnerable
DISSOCIATION    = 0       // 1 = high peritraumatic dissociation subtype (impairs contextual encoding, worse prognosis)

// ---- Sertraline (oral SSRI) PK ----
SERT_KA   = 0.35    // 1/h
SERT_F    = 0.88
SERT_V1   = 1200     // L (large Vd, lipophilic)
SERT_CL   = 90       // L/h
SERT_EC50 = 12       // relative plasma conc for half-max SERT occupancy
SERT_EMAX = 0.85

// ---- Paroxetine (oral SSRI, alternative) PK ----
PAROX_KA   = 0.40
PAROX_F    = 0.80
PAROX_V1   = 20       // L/kg-scaled central (large Vd)
PAROX_CL   = 5.0      // L/h (nonlinear/autoinhibition simplified linear)
PAROX_EC50 = 8
PAROX_EMAX = 0.85

// ---- Prazosin (oral alpha-1 antagonist) PK ----
PRAZ_KA   = 1.2
PRAZ_F    = 0.65
PRAZ_V1   = 55
PRAZ_CL   = 12
PRAZ_EC50 = 4.0      // relative plasma conc for half-max alpha-1 blockade
PRAZ_EMAX = 0.90

// ---- Ketamine/esketamine (IV/IN, rapid-acting) PK ----
KET_KE    = 0.30     // 1/h elimination (t1/2 ~2-3h)
KET_V1    = 200      // L
KET_KEO   = 0.15     // 1/h effect-compartment equilibration (biophase, mTOR/synaptogenesis delay)
KET_EC50  = 1.5
KET_EMAX  = 0.75     // acute extinction-facilitation ceiling per infusion

// ---- MDMA (oral, session-based adjunct to psychotherapy) PK ----
MDMA_KA   = 0.55
MDMA_F    = 0.90
MDMA_V1   = 300
MDMA_CL   = 24       // saturable CYP2D6 metabolism simplified linear
MDMA_EC50 = 6
MDMA_EMAX = 0.80     // acute amygdala-reactivity dampening during processing session

// ---- HPA axis kinetics ----
CORT_BASE     = 12       // ug/dL reference healthy baseline cortisol
CORT_K_IN     = 0.03     // 1/h rate toward target
GR_FEEDBACK   = 1.0      // gain on enhanced negative feedback (>1 = stronger suppression, "paradoxical low cortisol")

// ---- LC-noradrenergic tone kinetics ----
NE_BASE    = 0.15
NE_K_IN    = 0.05
NE_DRIVE   = 1.0         // scaling for CRH-LC positive feedback drive

// ---- Fear circuit kinetics ----
AMYG_K_IN   = 0.020      // 1/h amygdala reactivity index toward target
VMPFC_K_IN  = 0.015      // 1/h vmPFC inhibitory tone toward target
FEAR_K_IN   = 0.0060     // fear-memory trace accrual
EXT_K_IN    = 0.0060     // extinction-memory trace accrual
EXT_DECAY   = 0.0015     // spontaneous extinction-memory decay (renewal risk)

// ---- Symptom-cluster kinetics (0-40 CAPS-subscale-like units each) ----
K_INTRUSION   = 0.010
K_AVOIDANCE   = 0.010
K_NEGCOG      = 0.009
K_HYPERAROUSE = 0.011
SX_MAX        = 40       // ceiling per DSM-5 cluster subscale (sums to CAPS-5 total, 0-160-like scale rescaled to 0-80)

// ---- Sleep / therapy-dose kinetics ----
K_SLEEP      = 0.030
THERAPY_DECAY = 0.0     // cumulative therapy dose does not decay (permanent skill/extinction consolidation record)

$CMT @annotated
SERT_GUT    : Sertraline gut depot (mg)
SERT_CENT   : Sertraline plasma compartment (mg)
PAROX_GUT   : Paroxetine gut depot (mg)
PAROX_CENT  : Paroxetine plasma compartment (mg)
PRAZ_GUT    : Prazosin gut depot (mg)
PRAZ_CENT   : Prazosin plasma compartment (mg)
KET_CENT    : Ketamine/esketamine plasma compartment (mg)
KET_EFFECT  : Ketamine effect compartment (biophase, mTOR/synaptogenesis signal)
MDMA_GUT    : MDMA gut depot (mg)
MDMA_CENT   : MDMA plasma compartment (mg)
CORTISOL    : Serum cortisol (ug/dL)
NE_TONE     : Locus coeruleus/noradrenergic tone index (0-1+)
AMYG_REACT  : Amygdala reactivity index (0-1)
VMPFC_TONE  : vmPFC top-down inhibitory/extinction-capacity tone (0-1)
FEAR_MEM    : Conditioned fear-memory trace strength (0-1)
EXT_MEM     : Extinction-memory trace strength (0-1)
THERAPY_CUM : Cumulative trauma-focused psychotherapy dose (sessions, arbitrary units)
INTRUSION   : Cluster B intrusion/re-experiencing severity (0-40)
AVOIDANCE   : Cluster C avoidance severity (0-40)
NEGCOG      : Cluster D negative cognition/mood severity (0-40)
HYPERAROUSE : Cluster E hyperarousal severity (0-40)
SLEEP_DIST  : Sleep disturbance/nightmare index (0-1)
CAPS5       : Composite CAPS-5-like total severity score (0-160)
FX_WEEKS    : Elapsed simulation time tracker (weeks)

$MAIN
double SERT_Cp  = SERT_CENT/SERT_V1;
double PAROX_Cp = PAROX_CENT/PAROX_V1;
double PRAZ_Cp  = PRAZ_CENT/PRAZ_V1;
double MDMA_Cp  = MDMA_CENT/MDMA_V1;

// ---- SSRI effect (attenuates amygdala reactivity target, delayed via SS7-like onset already captured by plasma accumulation) ----
double sert_effect  = SERT_EMAX*(SERT_Cp/(SERT_EC50+SERT_Cp+1e-9));
double parox_effect = PAROX_EMAX*(PAROX_Cp/(PAROX_EC50+PAROX_Cp+1e-9));
double ssri_effect  = 1.0 - (1.0-sert_effect)*(1.0-parox_effect);
if (ssri_effect > 0.90) ssri_effect = 0.90;

// ---- Prazosin effect (alpha-1 blockade -> sleep/nightmare normalization, NOT direct daytime amygdala effect) ----
double praz_effect = PRAZ_EMAX*(PRAZ_Cp/(PRAZ_EC50+PRAZ_Cp+1e-9));

// ---- Ketamine effect (via effect-compartment biophase, acute extinction-facilitation pulse) ----
double ket_effect = KET_EMAX*(KET_EFFECT/(KET_EC50+KET_EFFECT+1e-9));

// ---- MDMA acute session effect (amygdala dampening during processing window) ----
double mdma_effect = MDMA_EMAX*(MDMA_Cp/(MDMA_EC50+MDMA_Cp+1e-9));

// ---- Composite trauma/resilience drive ----
double genetic_gr_gain = GR_FEEDBACK*(1.0 + 0.5*FKBP5_RISK);
double dissoc_gain     = 1.0 + 0.3*DISSOCIATION;
double vulnerability   = (TRAUMA_SEVERITY*dissoc_gain)/RESILIENCE;

$ODE
// ---------------- Sertraline PK ----------------
dxdt_SERT_GUT  = -SERT_KA*SERT_GUT;
dxdt_SERT_CENT =  SERT_KA*SERT_GUT*SERT_F - (SERT_CL/SERT_V1)*SERT_CENT;

// ---------------- Paroxetine PK ----------------
dxdt_PAROX_GUT  = -PAROX_KA*PAROX_GUT;
dxdt_PAROX_CENT =  PAROX_KA*PAROX_GUT*PAROX_F - (PAROX_CL/PAROX_V1)*PAROX_CENT;

// ---------------- Prazosin PK ----------------
dxdt_PRAZ_GUT  = -PRAZ_KA*PRAZ_GUT;
dxdt_PRAZ_CENT =  PRAZ_KA*PRAZ_GUT*PRAZ_F - (PRAZ_CL/PRAZ_V1)*PRAZ_CENT;

// ---------------- Ketamine/esketamine PK + effect compartment ----------------
dxdt_KET_CENT   = -KET_KE*KET_CENT;
dxdt_KET_EFFECT =  KET_KEO*(KET_CENT/KET_V1 - KET_EFFECT);

// ---------------- MDMA PK ----------------
dxdt_MDMA_GUT  = -MDMA_KA*MDMA_GUT;
dxdt_MDMA_CENT =  MDMA_KA*MDMA_GUT*MDMA_F - (MDMA_CL/MDMA_V1)*MDMA_CENT;

// ---------------- Cortisol (enhanced-feedback HPA axis; paradoxical low/normal baseline) ----------------
double cort_target = CORT_BASE*(1.0 + 0.25*vulnerability) / (1.0 + 0.6*genetic_gr_gain*vulnerability);
dxdt_CORTISOL = CORT_K_IN*(cort_target-CORTISOL);

// ---------------- Locus coeruleus / noradrenergic tone ----------------
double ne_target = NE_BASE + NE_DRIVE*0.35*vulnerability*(1.0 - 0.3*praz_effect);
if (ne_target < 0) ne_target = 0;
dxdt_NE_TONE = NE_K_IN*(ne_target-NE_TONE);

// ---------------- Amygdala reactivity ----------------
double amyg_target = 0.15 + 0.55*vulnerability*(1.0+0.4*NE_TONE)*(1.0-EXT_MEM)*(1.0-ssri_effect)*(1.0-0.5*mdma_effect);
if (amyg_target < 0.05) amyg_target = 0.05;
if (amyg_target > 1.0) amyg_target = 1.0;
dxdt_AMYG_REACT = AMYG_K_IN*(amyg_target-AMYG_REACT);

// ---------------- vmPFC inhibitory/extinction-capacity tone ----------------
double vmpfc_target = 0.75 - 0.45*vulnerability*(1.0-0.7*EXT_MEM);
if (vmpfc_target < 0.05) vmpfc_target = 0.05;
if (vmpfc_target > 1.0) vmpfc_target = 1.0;
dxdt_VMPFC_TONE = VMPFC_K_IN*(vmpfc_target-VMPFC_TONE);

// ---------------- Fear-memory trace ----------------
double fear_drive = FEAR_K_IN*AMYG_REACT*(1.0-EXT_MEM);
dxdt_FEAR_MEM = fear_drive*(1.0-FEAR_MEM) - EXT_DECAY*0.2*FEAR_MEM;

// ---------------- Extinction-memory trace (therapy/MDMA/ketamine-facilitated plasticity) ----------------
double ext_boost = 1.0 + 3.0*ket_effect + 2.0*mdma_effect + 0.02*THERAPY_CUM;
double ext_drive = EXT_K_IN*VMPFC_TONE*ext_boost;
dxdt_EXT_MEM = ext_drive*(1.0-EXT_MEM) - EXT_DECAY*EXT_MEM;

// ---------------- Cumulative psychotherapy dose (incremented via discrete session events) ----------------
dxdt_THERAPY_CUM = -THERAPY_DECAY*THERAPY_CUM;

// ---------------- DSM-5 symptom clusters ----------------
double intrusion_target   = SX_MAX*FEAR_MEM*(1.0-0.5*EXT_MEM);
double avoidance_target   = SX_MAX*0.85*FEAR_MEM*(1.0-0.4*EXT_MEM);
double negcog_target      = SX_MAX*0.8*(0.5*FEAR_MEM + 0.5*(1.0-VMPFC_TONE));
double hyperarouse_target = SX_MAX*NE_TONE/(0.5+NE_TONE);

dxdt_INTRUSION   = K_INTRUSION*(intrusion_target-INTRUSION);
dxdt_AVOIDANCE   = K_AVOIDANCE*(avoidance_target-AVOIDANCE);
dxdt_NEGCOG      = K_NEGCOG*(negcog_target-NEGCOG);
dxdt_HYPERAROUSE = K_HYPERAROUSE*(hyperarouse_target-HYPERAROUSE);

// ---------------- Sleep disturbance (NE surges during REM, normalized by prazosin) ----------------
double sleep_target = (NE_TONE/(0.4+NE_TONE))*(1.0-0.75*praz_effect);
if (sleep_target < 0) sleep_target = 0;
dxdt_SLEEP_DIST = K_SLEEP*(sleep_target-SLEEP_DIST);

// ---------------- Composite CAPS-5-like total ----------------
double caps_target = INTRUSION+AVOIDANCE+NEGCOG+HYPERAROUSE;
dxdt_CAPS5 = 0.15*(caps_target-CAPS5);

dxdt_FX_WEEKS = 1.0/(24.0*7.0);

$INIT
SERT_GUT = 0, SERT_CENT = 0, PAROX_GUT = 0, PAROX_CENT = 0,
PRAZ_GUT = 0, PRAZ_CENT = 0, KET_CENT = 0, KET_EFFECT = 0,
MDMA_GUT = 0, MDMA_CENT = 0,
CORTISOL = 12, NE_TONE = 0.15, AMYG_REACT = 0.15, VMPFC_TONE = 0.75,
FEAR_MEM = 0, EXT_MEM = 0, THERAPY_CUM = 0,
INTRUSION = 2, AVOIDANCE = 2, NEGCOG = 2, HYPERAROUSE = 2,
SLEEP_DIST = 0.1, CAPS5 = 8, FX_WEEKS = 0

$CAPTURE sert_effect parox_effect ssri_effect praz_effect ket_effect mdma_effect
$CAPTURE genetic_gr_gain dissoc_gain vulnerability cort_target ne_target amyg_target
$CAPTURE vmpfc_target fear_drive ext_boost ext_drive intrusion_target avoidance_target
$CAPTURE negcog_target hyperarouse_target sleep_target caps_target

## =============================================================================
## Treatment scenarios (see ptsd_shiny_app.R for interactive dosing UI)
##
## 1. Natural history, moderate trauma (TRAUMA_SEVERITY=1.0, untreated)
## 2. Natural history, severe/repeated trauma + FKBP5 risk + dissociative
##    subtype (TRAUMA_SEVERITY=1.5, FKBP5_RISK=1, DISSOCIATION=1 -> worst
##    prognosis, chronic course)
## 3. Sertraline 100-200 mg/day (FDA-approved SSRI, chronic oral)
## 4. Paroxetine 20-50 mg/day (FDA-approved SSRI, chronic oral)
## 5. Prazosin 1-15 mg qhs titrated (nightmare/sleep add-on to SSRI)
## 6. Prolonged Exposure / CPT / EMDR weekly sessions x12 (THERAPY_CUM
##    incremented per session event, no drug)
## 7. Ketamine 0.5 mg/kg IV x6 infusions over 2 weeks (rapid-acting adjunct)
## 8. MDMA-assisted psychotherapy: 3 dosing sessions (~120 mg + 60 mg
##    supplemental) embedded in 12 preparatory/integration therapy sessions
## 9. High resilience (RESILIENCE=1.5) + brief trauma -> natural remission
##    trajectory without treatment
## 10. Combination: SSRI + weekly trauma-focused psychotherapy +
##     prazosin for residual nightmares (multimodal standard-of-care)
##
## Example mrgsolve event code (sertraline chronic QD + weekly PE sessions):
##   mod <- mread("ptsd_mrgsolve_model") %>% param(TRAUMA_SEVERITY = 1.0)
##   e_sert <- ev(amt = 100, cmt = "SERT_GUT", time = 0, ii = 24, addl = 180)
##   e_pe   <- ev(amt = 1, cmt = "THERAPY_CUM", time = 168, ii = 168, addl = 11)
##   out <- mod %>% ev(e_sert) %>% ev(e_pe) %>% mrgsim(end = 24*7*52, delta = 24)
## =============================================================================
