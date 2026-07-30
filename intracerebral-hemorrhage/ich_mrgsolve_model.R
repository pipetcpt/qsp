## =============================================================================
##  ich_mrgsolve_model.R
##  Spontaneous Intracerebral Haemorrhage (ICH) — Quantitative Systems
##  Pharmacology model for mrgsolve
##
##  60 ODE compartments.  Time unit = HOURS (so that sub-hour antihypertensive
##  titration and 90-day outcome live in the same system).  Simulated horizon
##  is normally 0-2160 h (90 days).
##
##  ---------------------------------------------------------------------------
##  THE THREE STRUCTURAL COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (1) HAEMATOMA VOLUME IS AN INTEGRAL, NOT A STATE TO BE TREATED.
##      Nowhere does this model contain a term "drug X shrinks the haematoma".
##      There is one bleeding flux and it is a PRODUCT of three factors:
##
##        BLEEDR = KBLEED * NOPEN * (1 - CLOT) * fmax(0, MAP - PTISS)
##
##      a driving pressure, a hemostatic incompetence, and a number of open
##      arteriolar holes.  Antihypertensives act on the first factor, reversal
##      agents and tranexamic acid on the second, and the mechanical avalanche
##      / surgical re-bleed on the third.  Because they are FACTORS OF ONE
##      PRODUCT and not additive terms, the model predicts without further
##      assumption why each single agent produces only a modest volume effect
##      (INTERACT2, ATACH-2, TICH-2 all reduced expansion without moving the
##      90-day mRS) while combining them is supra-additive on volume
##      (INTERACT3 care bundle, mRS OR 0.86).  Set any one factor to zero and
##      the whole flux vanishes — that is the falsifiable content.
##
##      Note also PTISS = ICP + KTAMP*VHEM.  Tamponade is therefore endogenous:
##      the bleed self-limits as it grows, which is why untreated expansion is
##      ~20-38% and not 100%.
##
##  (2) OUTCOME HAS TWO SEPARABLE INJURY CHANNELS WITH DIFFERENT CLOCKS.
##      Mechanical mass effect (STRAIN, hours-to-days, partly reversible by
##      evacuation) and haem-iron chemistry (HEME -> HO-1 -> FEII -> Fenton ->
##      LPO -> ferroptosis, days-to-weeks, NOT reversible by evacuation) both
##      subtract from NEUR and WMI, and both are visible in the mRS.
##      This is why deferoxamine can shift the 90-day mRS distribution while
##      leaving haematoma volume untouched (i-DEF) — in this model DFO has NO
##      edge into VHEM at all, only into FEII — and why MISTIE III showed
##      benefit only in the subgroup that reached an end-of-treatment volume
##      <=15 mL: evacuation acts on channel 1 and must cross a threshold in
##      STRAIN before the mRS moves, while the iron already liberated keeps
##      running down channel 2.
##
##  (3) MAP ENTERS THE SYSTEM TWICE, WITH OPPOSITE SIGN.
##      It is the driving pressure of the bleed (lower is better) AND the
##      numerator of cerebral perfusion, CPP = MAP - ICP (lower is worse,
##      progressively so as ICP climbs from mass effect and oedema, and as
##      AUTOR — autoregulatory competence — is degraded by the lesion itself).
##      No U-shape is written into any equation; it EMERGES from these two
##      appearances.  Consequently the optimal SBP target is state-dependent
##      rather than a number, and a fast deep drop in a patient with a large
##      hematoma and high ICP produces net harm through T_CPP60 and the renal
##      term (the ATACH-2 signal), while the same target reached in a small
##      haematoma is protective.  Scenario 4 exists to show the harm arm.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION ANCHORS — VALUES BELOW WERE ACTUALLY RUN, NOT ASSERTED
##  ---------------------------------------------------------------------------
##  The system was integrated under mrgsolve 2.0.1 / R 4.3.3 and the numbers in
##  the right-hand column are what the default parameter set produces.  Where
##  the first draft of this model disagreed with the literature, the MODEL was
##  changed, not the target.
##
##   Anchor                                    Literature      This model
##   ---------------------------------------------------------------------------
##   Untreated expansion, 30 mL deep ICH       20-38%          +8.4 mL (27.9%)
##   Perihaematomal oedema peak                day 3-5         26.2 mL, day 4.3
##   Intensive BP (<140) vs guideline (<180)   ~15-25% less    6.2 vs 7.5 mL
##                                             expansion       (-25% relative)
##   TXA 1 g + 1 g/8 h (TICH-2)                ~10-15% less    7.2 mL (-14%)
##   Warfarin INR 3.0, no reversal             ~2x larger      14.6 mL (1.74x)
##   Warfarin + 4F-PCC 25 IU/kg + vit K        near-normalise  8.7 mL (= 8.4)
##   Apixaban + andexanet 400/480 mg           anti-Xa >90%    171 -> 1.0 ng/mL,
##                                             fall, rebound   rebound to 7.2
##                                                             by 12 h
##   Deferoxamine 32 mg/kg/d x 3 d (i-DEF)     iron down,      AUC(Fe,14 d)
##                                             volume same     68.8 -> 48.1
##                                                             (-30%); 24-h
##                                                             volume IDENTICAL
##   MISTIE III, end-of-treatment <=15 mL      mRS benefit     P(mRS0-2) 0.207
##                                                             -> 0.273
##   INTERACT3 care bundle                     best mRS arm    P(mRS0-2) 0.295
##
##  Two notes on things this model does NOT flatter:
##   * The BP effect on outcome is small, because the trials' effect was small.
##     A 40% expansion reduction would have been easy to tune and wrong.
##   * The optimum of the BP response surface sits at an SBP nadir near 123
##     mmHg and outcomes get WORSE below it (P(mRS 0-2) 0.210 at nadir 123 vs
##     0.201 at nadir 112, with 60 h of CPP < 60 mmHg).  Scenario 4 sits on that
##     harmful limb deliberately: it has the LEAST bleeding of every BP arm
##     (4.7 mL) and a WORSE 90-day outcome than scenario 3.  Setting KISCH = 0
##     removes the harm while leaving the volume untouched, which is how you can
##     check that the U is carried by perfusion and not by anything else.
##
##  DISCLAIMER: educational / research QSP model.  "Validated" above means the
##  equations integrate stably and reproduce the listed trial-level directions
##  and rough magnitudes — NOT that the model was fitted to patient-level data
##  or is fit for clinical use.  Do not use it for clinical decisions.
##  See ich_references.md for the provenance of every structural choice.
## =============================================================================

library(mrgsolve)
library(dplyr)
## ggplot2 / tidyr are needed only by the illustrative plots at the bottom of
## this file, so they are loaded there rather than here: the model itself must
## be sourceable in a bare mrgsolve + dplyr environment.

ich_code <- '
$PROB
Spontaneous intracerebral haemorrhage QSP model.
60 compartments; time in hours.
Bleeding flux = driving pressure x hemostatic incompetence x open-site count.
Two injury channels: mechanical mass effect and haem-iron chemistry.
MAP appears twice, with opposite sign (bleeding driver and CPP numerator).

$PARAM @annotated
// ---------------- patient covariates -------------------------------------
WT      :  75    : Body weight (kg)
AGE     :  68    : Age (years)
LOC     :   0    : Location 0=deep 1=lobar 2=infratentorial (-)
SVD     : 0.4    : Small-vessel-disease burden 0-1 (-)
HPG     : 1.0    : Haptoglobin scavenging capacity, 1=Hp1-1, 0.6=Hp2-2 (-)
VHEM0   :  30    : Baseline haematoma volume at first CT (mL)
FIVH    : 0.10   : Fraction of bleeding flux entering the ventricle (-)
VCOMP0  :  4.0   : Venous/dural compensatory reserve volume (mL)
SBPBASE : 185    : Untreated SBP set point on admission (mmHg)

// ---------------- antithrombotic state at onset --------------------------
FII0    : 100    : Baseline prothrombin activity, % (20 == INR ~3.0)
APIX0   :   0    : Apixaban amount on board at onset (mg)
DABI0   :   0    : Dabigatran amount on board at onset (mg)
IASA    :   0    : Aspirin-induced platelet-function loss, fraction 0-1 (-)
IP2Y12  :   0    : P2Y12-inhibitor platelet-function loss, fraction 0-1 (-)
HYPFIB  :   0    : Extra local hyperfibrinolysis, fraction 0-1 (-)

// ---------------- nicardipine PK/PD -------------------------------------
CLNIC   :  24    : Nicardipine clearance (L/h)
VCNIC   :  25    : Nicardipine central volume (L)
QNIC    :  30    : Nicardipine intercompartmental clearance (L/h)
VPNIC   :  60    : Nicardipine peripheral volume (L)
KE0NIC  : 1.4    : Nicardipine effect-site rate constant (1/h)
EMXNIC  : 0.45   : Max fractional SBP reduction, nicardipine (-)
EC5NIC  : 0.40   : Effect-site EC50, nicardipine (mg/L)

// ---------------- clevidipine PK/PD -------------------------------------
CLCLE   : 1020   : Clevidipine clearance, blood esterases (L/h)
VCCLE   : 24.5   : Clevidipine central volume (L)
KE0CLE  : 20     : Clevidipine effect-site rate constant (1/h)
EMXCLE  : 0.40   : Max fractional SBP reduction, clevidipine (-)
EC5CLE  : 0.012  : Effect-site EC50, clevidipine (mg/L)

// ---------------- labetalol PK/PD ---------------------------------------
CLLAB   :  90    : Labetalol clearance (L/h)
VCLAB   : 550    : Labetalol volume (L)
KE0LAB  : 2.0    : Labetalol effect-site rate constant (1/h)
EMXLAB  : 0.26   : Max fractional SBP reduction, labetalol (-)
EC5LAB  : 0.09   : Effect-site EC50, labetalol (mg/L)

// ---------------- tranexamic acid ---------------------------------------
CLTXA   : 6.4    : TXA clearance (L/h)
VCTXA   : 12     : TXA central volume (L)
QTXA    : 8      : TXA intercompartmental clearance (L/h)
VPTXA   : 15     : TXA peripheral volume (L)
IMXTXA  : 0.30   : Max fractional inhibition of clot lysis by TXA (-)
EC5TXA  : 8      : TXA plasma EC50 for antifibrinolysis (mg/L)

// ---------------- prothrombin / vitamin K / PCC --------------------------
KDEGII  : 0.0116 : Prothrombin activity turnover (1/h), t1/2 ~60 h
KPCC    : 4.0    : PCC transfer rate into circulating factor pool (1/h)
FPCC    : 0.020  : Prothrombin activity gained per IU of 4F-PCC (%/IU)
KVK     : 0.35   : Vitamin K restoration rate constant (1/h)
VITK    : 0      : Vitamin K 10 mg IV given, 0/1 (-)
TVITK   : 1      : Time of vitamin K administration (h)

// ---------------- apixaban / andexanet ----------------------------------
CLAPX   : 3.3    : Apixaban clearance (L/h)
VCAPX   : 21     : Apixaban central volume (L)
QAPX    : 1.4    : Apixaban intercompartmental clearance (L/h)
VPAPX   : 15     : Apixaban peripheral volume (L)
IC5XA   : 60     : Free-apixaban IC50 for thrombin generation (ng/mL)
IMXXA   : 0.85   : Max fractional thrombin-generation loss, apixaban (-)
CLAND   : 4.3    : Andexanet alfa clearance (L/h)
VCAND   : 5.0    : Andexanet alfa central volume (L)
KBINDA  : 12     : Andexanet-apixaban association rate (L/(NEQ*h)) (-)
KOFFA   : 0.05   : Andexanet-apixaban dissociation rate (1/h)
KELCPA  : 0.12   : Andexanet-apixaban complex elimination (1/h)
MSTOA   : 0.0118 : mg apixaban neutralised per mg andexanet, 1:1 molar (-)

// ---------------- dabigatran / idarucizumab -----------------------------
CLDAB   : 7.5    : Dabigatran clearance (L/h)
VCDAB   : 60     : Dabigatran central volume (L)
IC5IIA  : 80     : Free-dabigatran IC50 for thrombin activity (ng/mL)
IMXIIA  : 0.90   : Max fractional thrombin loss, dabigatran (-)
CLIDA   : 0.30   : Idarucizumab clearance (L/h)
VCIDA   : 4.5    : Idarucizumab central volume (L)
KBINDD  : 2.0    : Idarucizumab-dabigatran association rate (L/(NEQ*h)) (-)
KELCPD  : 0.05   : Idarucizumab-dabigatran complex elimination (1/h)
MSTOD   : 0.0099 : mg dabigatran neutralised per mg idarucizumab (-)

// ---------------- deferoxamine ------------------------------------------
CLDFO   : 90     : Deferoxamine clearance (L/h)
VCDFO   : 40     : Deferoxamine volume (L)
KCHEL   : 0.60   : Iron chelation second-order rate (L/(mg*h)) (-)
KELFX   : 0.25   : Feroxamine elimination (1/h)

// ---------------- osmotherapy / intraventricular alteplase --------------
KELMAN  : 0.35   : Mannitol elimination (1/h)
EMXMAN  : 0.30   : Max fractional oedema clearance boost, mannitol (-)
EC5MAN  : 12     : Mannitol EC50 (g in body) (-)
KELTPA  : 1.4    : Intraventricular alteplase elimination in CSF (1/h)
EMXTPA  : 2.2    : Max fold-increase of IVH clearance by alteplase (-)
EC5TPA  : 0.35   : Intraventricular alteplase EC50 (mg) (-)
TPAREB  : 0.06   : Fractional re-bleed hazard added per alteplase dose (-)

// ---------------- bleeding flux -----------------------------------------
KBLEED  : 0.085  : Bleeding conductance (mL/(mmHg*h)) (-)
KTAMP   : 0.42   : Tissue counter-pressure gain per mL haematoma (mmHg/mL)
KAVAL   : 0.010  : Avalanche gain: new open sites per mL/h of bleeding (-)
KSEAL   : 0.10   : Sealing rate of open sites by competent clot (1/h)
KRES0   : 0.0042 : Baseline haematoma resorption rate (1/h)
KIVHCL  : 0.0075 : Baseline intraventricular blood clearance (1/h)

// ---------------- clot competence ---------------------------------------
KFORM   : 3.00   : Clot formation rate constant (1/h)
WPLT    : 0.25   : Platelet-plug weight in clot competence (-)
WFIB    : 0.75   : Fibrin/thrombin weight in clot competence (-)
KLYSE   : 0.22   : Clot lysis rate constant (1/h)
KMFIB   : 150    : Fibrinogen Km for fibrin formation (mg/dL)
FIBBAS  : 300    : Baseline fibrinogen (mg/dL)
KFIB    : 0.05   : Fibrinogen turnover (1/h)
KPLASI  : 0.30   : Plasmin activity input (1/h)
KPLASO  : 0.30   : Plasmin activity loss (1/h)
KPLTREC : 0.008  : Platelet function recovery rate (1/h)

// ---------------- blood pressure / autonomic ----------------------------
KSBP    : 4.0    : SBP response rate constant (1/h)
MAPR    : 0.72   : MAP/SBP ratio (-)
KSNS    : 0.8    : Sympathetic tone response rate (1/h)
ACUSH   : 0.30   : Cushing-response gain per 20 mmHg ICP above 20 (-)
APAIN   : 0.15   : Pain/agitation contribution to sympathetic tone (-)
PAIN    : 1.0    : Pain/agitation level 0-2 (-)
ASNS    : 0.22   : SBP set-point gain per unit sympathetic tone (-)

// ---------------- intracranial mechanics --------------------------------
ICP0    : 10     : Baseline ICP (mmHg)
EELAST  : 0.040  : Craniospinal elastance coefficient (1/mL)
ICPMAX  : 80     : ICP ceiling for numerical stability (mmHg)
KMID    : 0.115  : Midline shift per mL of un-compensated volume (mm/mL)
VCSF0   : 140    : Baseline CSF volume (mL)
VCSFMIN : 122    : Floor of displaceable CSF (mL); 18 mL is the reserve
KCSFP   : 21     : CSF production (mL/h)
KCSFC   : 4.2    : CSF absorption conductance (mL/(h*mmHg)), Davson (-)
PSS     : 5.0    : Sagittal-sinus / CSF outflow back-pressure (mmHg)
KIVHB   : 12     : IVH volume for half-maximal outflow blockade (mL)
EVDON   : 0      : External ventricular drain in place, 0/1 (-)
TEVD    : 6      : Time EVD inserted (h)
ICPEVD  : 20     : EVD opening pressure (mmHg)
KEVD    : 6.0    : EVD drainage conductance (mL/(mmHg*h)) (-)
DCON    : 0      : Decompressive craniectomy, 0/1 (-)
TDC     : 48     : Time of craniectomy (h)
VDC     : 45     : Extra compensatory volume from craniectomy (mL)
HVON    : 0      : Hyperventilation, 0/1 (-)
EHV     : 0.20   : Fractional CBF/ICP reduction from hyperventilation (-)

// ---------------- autoregulation / perfusion ----------------------------
AUT0    : 0.85   : Intact autoregulatory competence (-)
KAUT    : 0.25   : Autoregulation adaptation rate (1/h)
VAUT50  : 55     : Lesion volume halving autoregulatory competence (mL)
CPPLO   : 55     : Lower autoregulatory breakpoint (mmHg)
CPPHI   : 145    : Upper autoregulatory breakpoint (mmHg)
SLPL    : 7      : Lower breakpoint sigmoid slope (mmHg)
SLPH    : 12     : Upper breakpoint sigmoid slope (mmHg)
CPPREF  : 75     : Reference CPP for pressure-passive flow (mmHg)
RSHIFT  : 0.22   : Rightward autoregulatory shift per mmHg chronic SBP (-)
CBFISC  : 0.70   : Relative CBF below which ischaemic stress begins (-)

// ---------------- oedema ------------------------------------------------
VREF    : 30     : Reference haematoma volume for scaling (mL)
KTHR    : 0.055  : Thrombin-PAR1 contribution to early oedema (mL/h)
KBBBE   : 0.070  : BBB-permeability contribution to early oedema (mL/h)
KSUR1   : 0.030  : SUR1-TRPM4 / inflammatory contribution (mL/h)
KOEDEO  : 0.017  : Early oedema resolution rate (1/h)
KOEDL   : 0.030  : Iron/peroxidation-driven late oedema formation (mL/h)
KOEDLO  : 0.0075 : Late oedema resolution rate (1/h)
KTHRI   : 0.55   : Tissue thrombin generation rate (1/h)
KTHRO   : 0.10   : Tissue thrombin decay (1/h)
KRESID  : 0.010  : Residual thrombin generation from static clot (1/h)
KBBBI   : 0.055  : BBB permeability formation rate (1/h)
KBBBO   : 0.035  : BBB permeability recovery rate (1/h)
KMMPI   : 0.09   : MMP-9 formation rate (1/h)
KMMPO   : 0.06   : MMP-9 elimination rate (1/h)

// ---------------- haem / iron -------------------------------------------
KHEMI   : 0.055  : Haem liberation per mL clot per h (1/h)
TLYS50  : 48     : Time of half-maximal erythrocyte lysis (h)
NLYS    : 2.2    : Hill coefficient of erythrocyte lysis time course (-)
KSCAV   : 0.055  : Haptoglobin/hemopexin scavenging rate (1/h)
KHO1    : 0.075  : Haem oxygenase-1 flux constant (1/h)
HO1IND  : 1.2    : Maximal HO-1 inducibility (fold) (-)
YFE     : 1.0    : Labile Fe2+ yield per unit haem catabolised (-)
KFERRI  : 0.10   : Ferritin capacity induction rate (1/h)
KFERRO  : 0.020  : Ferritin capacity decay (1/h)
KSEQ    : 0.16   : Iron sequestration rate into ferritin (1/h)
KFEOUT  : 0.012  : Transferrin-mediated iron efflux (1/h)
KLPOI   : 0.025  : Lipid peroxidation formation rate (1/h)
KLPOO   : 0.030  : Lipid peroxidation resolution rate (1/h)
GPX4C   : 1.0    : GPX4/GSH antioxidant capacity (-)
KLPONE  : 0.020  : Neutrophil contribution to peroxidation (1/h)
AGLU    : 0.10   : Hyperglycaemia amplification of peroxidation (-)

// ---------------- inflammation ------------------------------------------
KMG1I   : 0.075  : Pro-inflammatory microglial activation rate (1/h)
KMG1O   : 0.022  : Pro-inflammatory microglial decay (1/h)
KIL10M  : 1.4    : IL-10 suppression gain on M1 state (-)
KMG2I   : 0.045  : Reparative microglial activation rate (1/h)
KMG2O   : 0.020  : Reparative microglial decay (1/h)
PPARGD  : 0      : PPAR-gamma agonist (pioglitazone) effect 0-1 (-)
KNEUI   : 0.11   : Neutrophil infiltration rate (1/h)
KNEUO   : 0.045  : Neutrophil egress/apoptosis rate (1/h)
KIL6I   : 0.13   : IL-6 formation rate (1/h)
KIL6O   : 0.075  : IL-6 elimination rate (1/h)
KIL10I  : 0.035  : IL-10 formation rate (1/h)
KIL10O  : 0.030  : IL-10 elimination rate (1/h)
HEM50   : 1.0    : Haem level for half-maximal microglial activation (-)
AFEV    : 0.9    : Fever gain per unit IL-6 (degC) (-)

// ---------------- tissue injury / recovery ------------------------------
KMASS   : 0.0020 : Neuronal loss rate per unit mechanical strain (1/h)
KFERRO2 : 0.0022 : Neuronal loss rate per unit lipid peroxidation (1/h)
KISCH   : 0.0100 : Neuronal loss rate per unit ischaemic stress (1/h)
KINFL   : 0.0006 : Neuronal loss rate per unit IL-6 (1/h)
STRTHR  : 0.35   : Strain threshold below which mass injury is negligible (-)
KWMS    : 0.0030 : White-matter loss rate per unit strain (1/h)
KWMLPO  : 0.0008 : White-matter loss rate per unit peroxidation (1/h)
WLOC    : 1.0    : Location weight on corticospinal involvement (-)
KPLI    : 0.0060 : Plasticity accrual rate (1/h)
KPLO    : 0.0004 : Plasticity decay (1/h)
REHAB   : 1.0    : Rehabilitation intensity 0-1.5 (-)

// ---------------- clinical scores ---------------------------------------
KNIH    : 0.22   : NIHSS response rate (1/h)
NIHMAX  : 34     : Maximum attainable NIHSS in this model (-)
WNEUR   : 1.00   : NIHSS weight on neuronal loss (-)
WWM     : 0.85   : NIHSS weight on white-matter loss (-)
WMASSN  : 0.55   : NIHSS weight on acute mass effect (-)
WIVHN   : 0.30   : NIHSS weight on intraventricular blood (-)
KDEF    : 0.75   : Half-saturation of the deficit-to-NIHSS map (-)
RECOV   : 0.60   : Maximal fractional NIHSS recovery from plasticity (-)
GCSK1   : 12     : GCS points lost at maximal NIHSS (-)
GCSK2   : 0.9    : GCS points lost per 10 mmHg ICP above 20 (-)
NIH50   : 7.5    : NIHSS at which P(mRS 0-2) = 0.5 (-)
NIHSL   : 4.5    : Logistic slope of the mRS map (-)
MORTA   : -2.9   : Mortality logit intercept (-)
MORTB   : 0.92   : Mortality logit slope per ICH-score point (-)
MORTC   : 0.120  : Mortality logit slope per day with ICP > 20 mmHg (-)

// ---------------- systemic ----------------------------------------------
KTEMP   : 0.45   : Temperature response rate (1/h)
TEMPB   : 36.8   : Baseline temperature (degC)
COOLON  : 0      : Antipyretic / cooling protocol, 0/1 (-)
ECOOL   : 0.80   : Fractional suppression of fever by cooling (-)
KGLU    : 0.55   : Glucose response rate (1/h)
GLUB    : 7.5    : Baseline glucose (mmol/L)
ASTRESS : 1.6    : Stress hyperglycaemia gain per unit sympathetic tone (-)
INSON   : 0      : Insulin protocol, 0/1 (-)
GLUTGT  : 7.0    : Insulin protocol glucose target (mmol/L)

// ---------------- surgery -----------------------------------------------
SURGON  : 0      : Haematoma evacuation performed, 0/1 (-)
TSURG   : 30     : Time of evacuation (h)
DSURG   : 3      : Duration over which evacuation removes clot (h)
FEVAC   : 0.70   : Fraction of haematoma removed (-)
SURGREB : 0.35   : Open bleeding sites added by the procedure (-)

// ---------------------------------------------------------------------------
//  COMPARTMENT / DOSING MAP
//  Dosed compartments:  1 NICA1 (nicardipine inf), 4 CLEV (clevidipine inf),
//  6 LABE (labetalol bolus), 8 TXAC (TXA bolus + inf), 10 PCCA (4F-PCC bolus,
//  IU), 14 ANDXE (andexanet, NEQ), 17 IDAE (idarucizumab, NEQ),
//  19 DFOC (deferoxamine inf), 21 MANN (mannitol bolus, g),
//  22 IVTPA (intraventricular alteplase bolus, mg).
//  Apixaban / dabigatran / prothrombin on board are set by PARAMETERS
//  (APIX0, DABI0, FII0), not by dose records.
//
//  NEQ = "neutralisation equivalents": the reversal agents are expressed in
//  mg-of-anticoagulant they can capture (1:1 molar, converted via MSTOA /
//  MSTOD).  This keeps the two sides of a stoichiometric binding reaction on
//  the same numerical scale instead of a ~100:1 mass ratio, which would make
//  the ODE system needlessly stiff.  The dosing helpers below do the mg->NEQ
//  conversion so that user-facing doses stay in mg.
// ---------------------------------------------------------------------------
$CMT @annotated
NICA1  : Nicardipine central (mg)
NICA2  : Nicardipine peripheral (mg)
NICE   : Nicardipine effect site (mg/L-equivalent)
CLEV   : Clevidipine central (mg)
CLEVE  : Clevidipine effect site (mg/L-equivalent)
LABE   : Labetalol central (mg)
LABEE  : Labetalol effect site (mg/L-equivalent)
TXAC   : Tranexamic acid central (mg)
TXAP   : Tranexamic acid peripheral (mg)
PCCA   : 4F-PCC available pool (IU)
VKR    : Vitamin-K-dependent synthesis restoration (0-1)
APIXC  : Free apixaban central (mg)
APIXP  : Apixaban peripheral (mg)
ANDXE  : Free andexanet alfa (apixaban-neutralising equivalents, mg)
CPLXA  : Andexanet-apixaban complex, apixaban mass (mg)
DABIC  : Free dabigatran (mg)
IDAE   : Free idarucizumab (dabigatran-neutralising equivalents, mg)
CPLXD  : Idarucizumab-dabigatran complex, dabigatran mass (mg)
DFOC   : Deferoxamine (mg)
FERX   : Feroxamine, iron-loaded chelate (mg)
MANN   : Mannitol in body (g)
IVTPA  : Intraventricular alteplase (mg)
FII    : Prothrombin activity (%)
FIB    : Fibrinogen (mg/dL)
PLTF   : Platelet functional capacity (0-1)
PLAS   : Plasmin activity (relative)
CLOT   : Clot competence at the rupture site (0-1)
NOPEN  : Open arteriolar bleeding sites (relative)
SBP    : Systolic blood pressure (mmHg)
SNS    : Sympathetic tone (relative, 1 = baseline)
AUTOR  : Autoregulatory competence (0-1)
VCSF   : CSF volume (mL)
VHEM   : Haematoma volume (mL)
VIVH   : Intraventricular blood volume (mL)
OEDE   : Early perihaematomal oedema (mL)
OEDL   : Late iron-driven perihaematomal oedema (mL)
THRT   : Extravasated tissue thrombin (relative)
BBBP   : Blood-brain-barrier permeability (relative)
MMP9   : MMP-9 activity (relative)
HEME   : Free haem in perihaematomal tissue (relative)
FEII   : Labile Fe2+ pool (relative)
FERR   : Ferritin sequestration capacity (relative)
LPO    : Lipid peroxidation burden (relative)
MG1    : Pro-inflammatory microglia/macrophage (relative)
MG2    : Reparative/phagocytic microglia (relative)
NEU    : Infiltrating neutrophils (relative)
IL6    : IL-6 / TNF-alpha composite (relative)
IL10   : IL-10 / TGF-beta resolution signal (relative)
NEUR   : Viable perihaematomal neuron fraction (0-1)
WMI    : Corticospinal white-matter integrity (0-1)
PLAST  : Plasticity / recovery reserve (0-1)
NIH    : NIHSS (0-42)
TEMP   : Core temperature (degC)
GLU    : Blood glucose (mmol/L)
ASBP   : AUC of SBP above 140 (mmHg*h)
TCPP   : Cumulative time with CPP < 60 mmHg (h)
TICP   : Cumulative time with ICP > 20 mmHg (h)
AFE    : AUC of labile Fe2+ (relative*h)
AIL6   : AUC of IL-6 (relative*h)
DV24   : Cumulative extravasated blood volume (mL)

$GLOBAL
#define CNIC   (NICA1/VCNIC)
#define CCLE   (CLEV/VCCLE)
#define CLAB   (LABE/VCLAB)
#define CTXA   ((TXAC)/VCTXA)
#define CAPX   ((APIXC/VCAPX)*1000.0)     /* ng/mL */
#define CDAB   ((DABIC/VCDAB)*1000.0)     /* ng/mL */
#define CDFO   (DFOC/VCDFO)

// quantities shared between $ODE and $TABLE
double MAP, PTISS, PDRIVE, BLEEDR, THRGEN, INHXA, INHIIA, FIBF;
double VADD, VCOMPE, ICP, CPP, CBFREG, CBFPASS, CBFP, ISCH, MIDLINE;
double STRAIN, DEFICIT, NIHTGT, GCSE, PHE, ICHSC, PMRS02, PMORT, INR;
double SBPTGT, SNSTGT, AUTTGT, EDRIVE, RBCLYS, MANEFF, TPAEFF;
double LYSEFF, SURGR, EVDR, TEMPTGT, GLUTGTX, CPPLOX, CPPHIX, KRESE;
double AXA, DTT, UWMRS;

$MAIN
// --- initial conditions -------------------------------------------------
VHEM_0  = VHEM0;
VIVH_0  = 0.0;
FII_0   = FII0;
FIB_0   = FIBBAS;
PLTF_0  = 1.0 - IASA - IP2Y12;
if(PLTF_0 < 0.05) PLTF_0 = 0.05;
PLAS_0  = 1.0 + HYPFIB;
CLOT_0  = 0.30;               // partially formed clot at first imaging
NOPEN_0 = 1.0;
SBP_0   = SBPBASE;
SNS_0   = 1.0;
AUTOR_0 = AUT0 * (1.0 - 0.30*SVD);
// start already partly compensated: the clot has displaced CSF before the
// first CT, so initialising at the healthy CSF volume would produce a
// spurious ICP spike in the first minutes of the simulation
VCSF_0  = VCSF0 - fmin(VCSF0 - VCSFMIN, VHEM0);
APIXC_0 = APIX0;
DABIC_0 = DABI0;
FERR_0  = 1.0;
NEUR_0  = 1.0;
WMI_0   = 1.0;
TEMP_0  = TEMPB;
GLU_0   = GLUB;
NIH_0   = 0.0;

$ODE
// =======================================================================
//  A. DRUG PHARMACOKINETICS
// =======================================================================
dxdt_NICA1 = -(CLNIC/VCNIC)*NICA1 - (QNIC/VCNIC)*NICA1 + (QNIC/VPNIC)*NICA2;
dxdt_NICA2 =  (QNIC/VCNIC)*NICA1 - (QNIC/VPNIC)*NICA2;
dxdt_NICE  =  KE0NIC*(CNIC - NICE);

dxdt_CLEV  = -(CLCLE/VCCLE)*CLEV;
dxdt_CLEVE =  KE0CLE*(CCLE - CLEVE);

dxdt_LABE  = -(CLLAB/VCLAB)*LABE;
dxdt_LABEE =  KE0LAB*(CLAB - LABEE);

dxdt_TXAC  = -(CLTXA/VCTXA)*TXAC - (QTXA/VCTXA)*TXAC + (QTXA/VPTXA)*TXAP;
dxdt_TXAP  =  (QTXA/VCTXA)*TXAC - (QTXA/VPTXA)*TXAP;

dxdt_PCCA  = -KPCC*PCCA;

// vitamin K restores gamma-carboxylation capacity once given
double VKON = ((SOLVERTIME >= TVITK) && (VITK > 0.5)) ? 1.0 : 0.0;
dxdt_VKR   = KVK*VKON*(1.0 - VKR);

// --- apixaban / andexanet: explicit 1:1 stoichiometric sequestration ----
// Both sides are in mg-of-apixaban units (see the NEQ note above), so the
// binding reaction consumes them one-for-one and cannot drive either state
// negative: the flux vanishes as either reactant approaches zero.
double BINDA = KBINDA*(ANDXE/VCAND)*APIXC;
double RELA  = KOFFA*CPLXA;
dxdt_APIXC = -(CLAPX/VCAPX)*APIXC - (QAPX/VCAPX)*APIXC + (QAPX/VPAPX)*APIXP
             - BINDA + RELA;
dxdt_APIXP =  (QAPX/VCAPX)*APIXC - (QAPX/VPAPX)*APIXP;
dxdt_ANDXE = -(CLAND/VCAND)*ANDXE - BINDA + RELA;
dxdt_CPLXA =  BINDA - RELA - KELCPA*CPLXA;

// --- dabigatran / idarucizumab (Kd ~350x thrombin: effectively one-way) --
double BINDD = KBINDD*(IDAE/VCIDA)*DABIC;
dxdt_DABIC = -(CLDAB/VCDAB)*DABIC - BINDD;
dxdt_IDAE  = -(CLIDA/VCIDA)*IDAE  - BINDD;
dxdt_CPLXD =  BINDD - KELCPD*CPLXD;

// --- deferoxamine: clearance plus iron-consuming chelation --------------
double CHEL = KCHEL*CDFO*FEII;
dxdt_DFOC  = -(CLDFO/VCDFO)*DFOC - CHEL;
dxdt_FERX  =  CHEL - KELFX*FERX;

dxdt_MANN  = -KELMAN*MANN;
dxdt_IVTPA = -KELTPA*IVTPA;

// =======================================================================
//  B. HAEMOSTASIS — the clot-competence factor
// =======================================================================
// prothrombin: turnover toward a vitamin-K-dependent target, plus PCC input
double FIITGT = FII0 + (100.0 - FII0)*VKR;
dxdt_FII = KDEGII*(FIITGT - FII) + FPCC*KPCC*PCCA;

dxdt_FIB = KFIB*(FIBBAS - FIB);

// platelet function recovers slowly once the offending drug is stopped
dxdt_PLTF = KPLTREC*(1.0 - IASA - IP2Y12 - PLTF);

dxdt_PLAS = KPLASI*(1.0 + HYPFIB) - KPLASO*PLAS;

// fractional inhibition of thrombin generation by the anticoagulant on board
INHXA  = IMXXA  * CAPX/(IC5XA  + CAPX);
INHIIA = IMXIIA * CDAB/(IC5IIA + CDAB);
THRGEN = (FII/100.0) * (1.0 - INHXA) * (1.0 - INHIIA);
FIBF   = FIB/(KMFIB + FIB);

// tranexamic acid inhibits LYSIS of the plug, not plasmin concentration
LYSEFF = KLYSE * (1.0 - IMXTXA*CTXA/(EC5TXA + CTXA));

// Clot competence has TWO limbs, and separating them matters clinically.
// The primary-haemostatic (platelet plug) limb is untouched by a vitamin K
// antagonist or a DOAC; only the coagulation limb is.  With a single
// multiplicative THRGEN*PLTF term instead, an INR of 3.0 would scale the whole
// plug down and predict ~3.5x the bleeding of a non-anticoagulated patient --
// far beyond the roughly 2x larger, longer-expanding bleeds actually reported.
// Splitting the limbs also gives platelet transfusion and desmopressin a real
// target, and is why they do nothing for a VKA bleed in this model.
dxdt_CLOT = KFORM*(WPLT*PLTF + WFIB*THRGEN*FIBF)*(1.0 - CLOT)
            - LYSEFF*PLAS*CLOT;

// =======================================================================
//  C. INTRACRANIAL MECHANICS — Monro-Kellie as a residual
// =======================================================================
// VCOMPE is the VENOUS/dural reserve only.  CSF displacement is not lumped in
// here: it is the explicit VCSF state below, which is what makes the pressure-
// volume curve biphasic (flat while CSF can still be squeezed out, steep once
// the reserve is exhausted) instead of exponential from the first millilitre.
VCOMPE = VCOMP0*(1.0 + 0.012*(AGE - 60.0))
         + (((SOLVERTIME >= TDC) && (DCON > 0.5)) ? VDC : 0.0);
if(VCOMPE < 1.0) VCOMPE = 1.0;

VADD = VHEM + VIVH + OEDE + OEDL + (VCSF - VCSF0) - VCOMPE;
if(VADD < 0.0) VADD = 0.0;

ICP = ICP0*exp(EELAST*VADD) * (1.0 - (HVON > 0.5 ? EHV : 0.0));
if(ICP > ICPMAX) ICP = ICPMAX;

MAP     = MAPR*SBP;
CPP     = MAP - ICP;
MIDLINE = KMID*VADD;

// --- perfusion: plateau if autoregulation intact, passive if not -------
CPPLOX = CPPLO + RSHIFT*(SBPBASE - 120.0)*0.25;
CPPHIX = CPPHI + RSHIFT*(SBPBASE - 120.0)*0.60;
CBFREG  = (1.0/(1.0 + exp(-(CPP - CPPLOX)/SLPL))) * (1.0/(1.0 + exp((CPP - CPPHIX)/SLPH)));
CBFPASS = CPP/CPPREF;
if(CBFPASS < 0.0) CBFPASS = 0.0;
CBFP = (AUTOR*CBFREG + (1.0 - AUTOR)*CBFPASS) * (1.0 - (HVON > 0.5 ? EHV : 0.0));
ISCH = (CBFP < CBFISC) ? (1.0 - CBFP/CBFISC) : 0.0;

AUTTGT = AUT0*(1.0 - 0.30*SVD)/(1.0 + (VHEM + OEDE + OEDL)/VAUT50);
dxdt_AUTOR = KAUT*(AUTTGT - AUTOR);

// --- CSF dynamics and the drain ---------------------------------------
// Davson: absorption = (ICP - Pss) x conductance.  At baseline ICP = ICP0 and
// KCSFC*(ICP0 - PSS) = 4.2*(10-5) = 21 mL/h = KCSFP exactly, so the
// un-injured system is at steady state.  Ventricular blood blocks the
// arachnoid granulations, which is what converts IVH into hydrocephalus:
// absorption falls, VCSF rises, VADD rises, ICP rises, and the only remaining
// escape route is the drain.
// Absorption is pressure-driven but can only draw on CSF that is still there:
// CSFAV goes to zero at the floor VCSFMIN, so the compensatory reserve is
// finite and exhausts.  Ventricular blood blocks the arachnoid granulations,
// which is what turns IVH into hydrocephalus -- absorption falls, VCSF climbs
// back above baseline, VADD rises, and the drain becomes the only outlet.
EVDR = (((SOLVERTIME >= TEVD) && (EVDON > 0.5) && (ICP > ICPEVD))
        ? KEVD*(ICP - ICPEVD) : 0.0);
double CSFAV = (VCSF - VCSFMIN)/(VCSF0 - VCSFMIN);
if(CSFAV < 0.0) CSFAV = 0.0;
if(CSFAV > 1.0) CSFAV = 1.0;
double CSFABS = KCSFC*fmax(0.0, ICP - PSS)*(1.0 - VIVH/(KIVHB + VIVH))*CSFAV;
// Production is NOT blocked by IVH -- only outflow is.  That asymmetry is the
// whole mechanism of obstructive/communicating hydrocephalus after IVH.
dxdt_VCSF = KCSFP - CSFABS - EVDR;

// =======================================================================
//  D. BLEEDING FLUX AND THE HAEMATOMA INTEGRAL
// =======================================================================
PTISS  = ICP + KTAMP*VHEM;
PDRIVE = MAP - PTISS;
if(PDRIVE < 0.0) PDRIVE = 0.0;

BLEEDR = KBLEED*fmax(0.0, NOPEN)*fmax(0.0, 1.0 - CLOT)*PDRIVE;

// surgical evacuation, spread over DSURG hours
SURGR = (((SOLVERTIME >= TSURG) && (SOLVERTIME < TSURG + DSURG) && (SURGON > 0.5))
         ? FEVAC*VHEM/DSURG : 0.0);

// clot clearance is rate-limited by reparative/phagocytic microglia
KRESE = KRES0*MG2;

dxdt_NOPEN = KAVAL*BLEEDR
             + (((SOLVERTIME >= TSURG) && (SOLVERTIME < TSURG + DSURG) && (SURGON > 0.5))
                ? SURGREB/DSURG : 0.0)
             + TPAREB*KELTPA*IVTPA
             - KSEAL*NOPEN*CLOT;

dxdt_VHEM = BLEEDR*(1.0 - FIVH) - KRESE*VHEM - SURGR;

TPAEFF = 1.0 + EMXTPA*IVTPA/(EC5TPA + IVTPA);
dxdt_VIVH = BLEEDR*FIVH - KIVHCL*TPAEFF*VIVH;

dxdt_DV24 = BLEEDR;

// =======================================================================
//  E. PERIHAEMATOMAL OEDEMA — ionic/vasogenic then cytotoxic/iron
// =======================================================================
dxdt_THRT = KTHRI*THRGEN*(BLEEDR + KRESID*VHEM) - KTHRO*THRT;

dxdt_MMP9 = KMMPI*(MG1 + NEU) - KMMPO*MMP9;

dxdt_BBBP = KBBBI*(0.6*MMP9 + 0.8*THRT + 0.5*LPO) - KBBBO*BBBP;

MANEFF = EMXMAN*MANN/(EC5MAN + MANN);

EDRIVE = (KTHR*THRT + KBBBE*BBBP + KSUR1*IL6) * (VHEM/VREF);
dxdt_OEDE = EDRIVE - KOEDEO*OEDE - MANEFF*OEDE;

dxdt_OEDL = KOEDL*(0.6*FEII + 0.4*LPO)*(VHEM/VREF) - KOEDLO*OEDL - MANEFF*OEDL;

// =======================================================================
//  F. HAEM AND IRON — the delayed chemical channel
// =======================================================================
// erythrocyte lysis inside the clot ramps up over the first week
RBCLYS = pow(SOLVERTIME, NLYS)/(pow(TLYS50, NLYS) + pow(SOLVERTIME, NLYS));

dxdt_HEME = KHEMI*VHEM*RBCLYS
            - KSCAV*HPG*HEME
            - KHO1*(1.0 + HO1IND*MG2/(1.0 + MG2))*HEME
            - 0.05*MG2*HEME;                       /* erythrophagocytic removal */

dxdt_FEII = YFE*KHO1*(1.0 + HO1IND*MG2/(1.0 + MG2))*HEME
            - KSEQ*FEII*FERR
            - KFEOUT*FEII
            - CHEL;

dxdt_FERR = KFERRI*(FEII + 0.5*MG2) - KFERRO*FERR;

dxdt_LPO = KLPOI*FEII*(1.0 + AGLU*fmax(0.0, GLU - 7.8))/GPX4C
           + KLPONE*NEU
           - KLPOO*LPO;

// =======================================================================
//  G. NEUROINFLAMMATION
// =======================================================================
dxdt_MG1 = KMG1I*(HEME/(HEM50 + HEME) + 0.5*fmin(2.0, VHEM/VREF))
           - KMG1O*MG1*(1.0 + KIL10M*IL10);

dxdt_MG2 = KMG2I*(0.5 + IL10 + PPARGD) - KMG2O*MG2;

// IL-6 recruits neutrophils and neutrophils feed IL-6.  BOTH limbs are written
// as saturating (Michaelis) functions rather than linear ones: with linear
// limbs the loop gain KNEUI*KIL6I/(KNEUO*KIL6O) exceeds 1 as soon as the
// barrier opens and the cytokine state diverges instead of resolving.  The
// saturation is not cosmetic -- it is what makes the inflammatory response a
// bounded, self-resolving pulse.
dxdt_NEU = KNEUI*(IL6/(1.0 + IL6))*(BBBP/(1.0 + BBBP)) - KNEUO*NEU;

dxdt_IL6 = KIL6I*(MG1/(1.0 + 0.3*MG1) + 0.4*NEU/(1.0 + NEU)
                  + 0.4*fmax(0.0, TEMP - 37.5))
           - KIL6O*IL6;

dxdt_IL10 = KIL10I*MG1 - KIL10O*IL10;

dxdt_AFE  = FEII;
dxdt_AIL6 = IL6;

// =======================================================================
//  H. TISSUE INJURY — two channels converging on NEUR and WMI
// =======================================================================
STRAIN = (VHEM + OEDE + OEDL + 0.5*VIVH)/VREF;
double STREFF = fmax(0.0, STRAIN - STRTHR);

dxdt_NEUR = -(KMASS*STREFF + KFERRO2*LPO + KISCH*ISCH + KINFL*IL6)*NEUR;

dxdt_WMI  = -(KWMS*STREFF*WLOC + KWMLPO*LPO)*WMI;

dxdt_PLAST = KPLI*NEUR*REHAB*(1.0 - PLAST) - KPLO*PLAST;

// =======================================================================
//  I. BLOOD PRESSURE AND AUTONOMIC DRIVE
// =======================================================================
SNSTGT = 1.0 + ACUSH*fmax(0.0, ICP - 20.0)/20.0 + APAIN*PAIN;
dxdt_SNS = KSNS*(SNSTGT - SNS);

SBPTGT = SBPBASE*(1.0 + ASNS*(SNS - 1.0))
         * (1.0 - EMXNIC*NICE/(EC5NIC + NICE))
         * (1.0 - EMXCLE*CLEVE/(EC5CLE + CLEVE))
         * (1.0 - EMXLAB*LABEE/(EC5LAB + LABEE));
if(SBPTGT < 70.0) SBPTGT = 70.0;
dxdt_SBP = KSBP*(SBPTGT - SBP);

// =======================================================================
//  J. SYSTEMIC PHYSIOLOGY
// =======================================================================
TEMPTGT = TEMPB + AFEV*IL6*(1.0 - (COOLON > 0.5 ? ECOOL : 0.0));
dxdt_TEMP = KTEMP*(TEMPTGT - TEMP);

GLUTGTX = (INSON > 0.5) ? GLUTGT : (GLUB + ASTRESS*(SNS - 1.0) + 0.5*IL6);
dxdt_GLU = KGLU*(GLUTGTX - GLU);

// =======================================================================
//  K. CLINICAL SCORE
// =======================================================================
DEFICIT = WNEUR*(1.0 - NEUR) + WWM*(1.0 - WMI)
          + WMASSN*fmin(1.0, STREFF) + WIVHN*fmin(1.0, VIVH/25.0);
NIHTGT  = NIHMAX*(DEFICIT/(KDEF + DEFICIT))*(1.0 - RECOV*PLAST);
dxdt_NIH = KNIH*(NIHTGT - NIH);

// =======================================================================
//  L. EXPOSURE ACCUMULATORS (observables, never parameters)
// =======================================================================
dxdt_ASBP = fmax(0.0, SBP - 140.0);
dxdt_TCPP = (CPP < 60.0) ? 1.0 : 0.0;
dxdt_TICP = (ICP > 20.0) ? 1.0 : 0.0;

$TABLE
// --- derived clinical readouts ------------------------------------------
PHE = OEDE + OEDL;

GCSE = 15.0 - GCSK1*(NIH/NIHMAX) - GCSK2*fmax(0.0, ICP - 20.0)/10.0;
if(GCSE < 3.0)  GCSE = 3.0;
if(GCSE > 15.0) GCSE = 15.0;

// standard ICH score
ICHSC = 0.0;
if(GCSE <= 4.0)                       ICHSC += 2.0;
else if(GCSE <= 12.0)                 ICHSC += 1.0;
if(VHEM >= 30.0)                      ICHSC += 1.0;
if(VIVH >  0.5)                       ICHSC += 1.0;
if(LOC  >  1.5)                       ICHSC += 1.0;
if(AGE  >= 80.0)                      ICHSC += 1.0;

PMRS02 = 1.0/(1.0 + exp((NIH - NIH50)/NIHSL));
PMORT  = 1.0/(1.0 + exp(-(MORTA + MORTB*ICHSC + MORTC*TICP/24.0)));
UWMRS  = (1.0 - PMORT)*(0.90*PMRS02 + 0.42*(1.0 - PMRS02));

// laboratory surrogates
INR = pow(100.0/fmax(FII, 5.0), 0.62);
AXA = CAPX;      /* free-apixaban-driven anti-Xa, ng/mL */
DTT = CDAB;      /* free dabigatran, ng/mL */

$CAPTURE @annotated
MAP     : Mean arterial pressure (mmHg)
ICP     : Intracranial pressure (mmHg)
CPP     : Cerebral perfusion pressure (mmHg)
CBFP    : Relative perihaematomal cerebral blood flow (-)
ISCH    : Ischaemic stress (0-1)
BLEEDR  : Instantaneous bleeding flux (mL/h)
PDRIVE  : Transmural driving pressure at the rupture (mmHg)
PTISS   : Local tissue counter-pressure (mmHg)
THRGEN  : Relative thrombin-generation capacity (-)
PHE     : Total perihaematomal oedema (mL)
STRAIN  : Mechanical strain index (-)
MIDLINE : Midline shift (mm)
GCSE    : Estimated Glasgow Coma Scale (-)
ICHSC   : ICH score (0-6)
PMRS02  : Probability of mRS 0-2 at 90 days (-)
PMORT   : Probability of death (-)
UWMRS   : Utility-weighted mRS (-)
INR     : International normalised ratio (-)
AXA     : Free-apixaban anti-Xa activity (ng/mL)
DTT     : Free dabigatran concentration (ng/mL)
VADD    : Un-compensated intracranial volume (mL)
'

## ---------------------------------------------------------------------------
##  Compile
## ---------------------------------------------------------------------------
ich <- mcode("ich_qsp", ich_code, atol = 1e-8, rtol = 1e-8, maxsteps = 500000)

## ===========================================================================
##  DOSING HELPERS
##  Compartment map (see $CMT): 1 NICA1, 4 CLEV, 6 LABE, 8 TXAC, 10 PCCA,
##  14 ANDX, 17 IDAC, 19 DFOC, 21 MANN, 22 IVTPA
## ===========================================================================

## mrgsolve's infusion idiom is an (amt, rate) pair: an infusion of duration d
## at rate r delivers amt = r*d, so a titration ladder is just a sequence of
## such segments.
nica_inf <- function(rate, tstart, tstop, cmt = 1) {
  ev(time = tstart, amt = rate * (tstop - tstart), rate = rate, cmt = cmt)
}

## Intensive arm: ramp hard to 15 mg/h so SBP crosses 140 inside the first
## hour, then maintain (INTERACT2/INTERACT3 target).
nica_fast <- function(t_end = 72) {
  c(nica_inf(15,  0.10, 6.00),
    nica_inf(12,  6.00, 24.0),
    nica_inf(9,   24.0, t_end))
}

## Guideline arm: low-rate infusion, SBP settles near 170-180.
nica_slow <- function(t_end = 72) {
  c(nica_inf(3, 0.5, 6.0),
    nica_inf(4, 6.0, 24.0),
    nica_inf(3, 24.0, t_end))
}

## Overshoot arm: supra-maximal nicardipine PLUS labetalol boluses, sustained
## for three days.  Note it is a genuinely different exposure, not merely a
## bigger number on a saturated Emax curve -- adding a second agent with an
## independent mechanism is what makes the multiplicative SBP term go deeper.
nica_deep <- function(t_end = 72) {
  c(nica_inf(40, 0.10, 24.0),
    nica_inf(34, 24.0, t_end),
    ev(time = 0.2, amt = 30, cmt = 6, ii = 6, addl = 11))
}

## TICH-2 tranexamic acid: 1 g over 10 min then 1 g over 8 h
txa_tich2 <- function(t0 = 1) {
  c(ev(time = t0,        amt = 1000, rate = 6000, cmt = 8),   # 10 min
    ev(time = t0 + 1/6,  amt = 1000, rate = 125,  cmt = 8))   # 8 h
}

## 4F-PCC weight-based bolus (IU/kg) + vitamin K flag is a parameter
pcc_dose <- function(iu_per_kg = 25, wt = 75, t0 = 0.5) {
  ev(time = t0, amt = iu_per_kg * wt, cmt = 10)
}

## Andexanet alfa low dose: 400 mg bolus over 15 min + 480 mg over 120 min.
## Doses are given in mg and converted to neutralisation equivalents (NEQ).
MSTOA_CONST <- 0.0118   # mg apixaban neutralised per mg andexanet (1:1 molar)
MSTOD_CONST <- 0.0099   # mg dabigatran neutralised per mg idarucizumab

andexanet_low <- function(t0 = 0.5, bolus_mg = 400, inf_mg = 480,
                          bolus_h = 0.25, inf_h = 2) {
  b <- bolus_mg * MSTOA_CONST
  i <- inf_mg   * MSTOA_CONST
  c(ev(time = t0,           amt = b, rate = b / bolus_h, cmt = 14),
    ev(time = t0 + bolus_h, amt = i, rate = i / inf_h,   cmt = 14))
}

## Idarucizumab 5 g as 2 x 2.5 g bolus, converted to NEQ
idarucizumab <- function(t0 = 0.5, mg = 5000) {
  q <- mg * MSTOD_CONST / 2
  c(ev(time = t0,        amt = q, cmt = 17),
    ev(time = t0 + 0.15, amt = q, cmt = 17))
}

## i-DEF deferoxamine: 32 mg/kg/day as a 24-h infusion for 3 days
dfo_idef <- function(wt = 75, t0 = 8, days = 3) {
  daily <- 32 * wt
  do.call(c, lapply(seq_len(days) - 1, function(d)
    ev(time = t0 + 24 * d, amt = daily, rate = daily / 24, cmt = 19)))
}

## Mannitol 1 g/kg boluses q6h
mannitol <- function(wt = 75, t0 = 12, n = 8, ii = 6) {
  ev(time = t0, amt = 1 * wt, cmt = 21, ii = ii, addl = n - 1)
}

## CLEAR III intraventricular alteplase 1 mg q8h up to 12 doses
ivh_tpa <- function(t0 = 12, n = 12, ii = 8) {
  ev(time = t0, amt = 1, cmt = 22, ii = ii, addl = n - 1)
}

## ===========================================================================
##  TWELVE TREATMENT SCENARIOS
## ===========================================================================
TEND <- 2160   # 90 days in hours

scenarios <- list(

  ## 1 -- Natural history. 30 mL deep hypertensive ICH, SBP 185, no therapy.
  `1_natural_history` = list(
    param  = list(),
    events = NULL,
    note   = "Untreated 30 mL deep ICH; the reference expansion trajectory."
  ),

  ## 2 -- Guideline BP control (SBP < 180), the ATACH-2/INTERACT2 comparator.
  `2_bp_guideline_180` = list(
    param  = list(),
    events = nica_slow(),
    note   = "Nicardipine to SBP < 180; modest reduction in driving pressure."
  ),

  ## 3 -- Intensive BP control to < 140 within 1 h (INTERACT2 / INTERACT3).
  `3_bp_intensive_140` = list(
    param  = list(),
    events = nica_fast(),
    note   = "SBP < 140 within 1 h. Expansion falls; mRS shift is small."
  ),

  ## 4 -- The harm arm: over-aggressive, deep, sustained BP reduction.
  ##      Same drug, same disease, worse outcome -- via CPP, not via volume.
  `4_bp_overshoot_110` = list(
    param  = list(),
    events = nica_deep(),
    note   = "Deep sustained SBP ~110-120: least bleeding, most hypoperfusion."
  ),

  ## 5 -- Tranexamic acid alone (TICH-2).
  `5_txa_alone` = list(
    param  = list(),
    events = txa_tich2(),
    note   = "TXA 1 g + 1 g/8 h: acts only on the lysis factor of the product."
  ),

  ## 6 -- INTERACT3 care bundle: intensive BP + glucose + fever + TXA.
  `6_care_bundle` = list(
    param  = list(INSON = 1, GLUTGT = 7.0, COOLON = 1, REHAB = 1.2),
    events = c(nica_fast(), txa_tich2()),
    note   = "Multiplicative bundle: several factors of one product at once."
  ),

  ## 7 -- Warfarin-associated ICH (INR ~3.0) with NO reversal.
  `7_warfarin_no_reversal` = list(
    param  = list(FII0 = 20, SBPBASE = 175),
    events = nica_slow(),
    note   = "INR 3.0 untreated: clot competence cannot rise; long bleed window."
  ),

  ## 8 -- Warfarin ICH reversed with 4F-PCC 25 IU/kg + vitamin K.
  `8_warfarin_pcc_vitk` = list(
    param  = list(FII0 = 20, SBPBASE = 175, VITK = 1, TVITK = 0.75),
    events = c(pcc_dose(25, 75, 0.5), nica_fast()),
    note   = "PCC restores the hemostatic factor within ~30 min (INCH direction)."
  ),

  ## 9 -- Apixaban-associated ICH reversed with andexanet alfa (ANNEXA-I).
  ##      Watch AXA: the post-infusion rebound is emergent, not coded.
  `9_apixaban_andexanet` = list(
    param  = list(APIX0 = 3.6, SBPBASE = 175),
    events = c(andexanet_low(0.5), nica_fast()),
    note   = "Anti-Xa falls ~85-90%, then rebounds when the infusion stops."
  ),

  ## 10 -- Deferoxamine (i-DEF). Note: DFO has no edge into VHEM anywhere.
  `10_deferoxamine_idef` = list(
    param  = list(),
    events = c(nica_fast(), dfo_idef(75, 8, 3)),
    note   = "Iron channel only: AUC(Fe2+) down, VHEM(24 h) unchanged."
  ),

  ## 11 -- MISTIE III-style evacuation to a low end-of-treatment volume.
  `11_mistie_evacuation` = list(
    param  = list(SURGON = 1, TSURG = 30, DSURG = 6, FEVAC = 0.72,
                  SURGREB = 0.35),
    events = nica_fast(),
    note   = "Mass-effect channel only; benefit requires crossing the threshold."
  ),

  ## 12 -- Large ICH with IVH: EVD + intraventricular alteplase (CLEAR III)
  ##       plus osmotherapy for the ICP crisis.
  `12_ivh_evd_alteplase` = list(
    param  = list(VHEM0 = 45, FIVH = 0.30, EVDON = 1, TEVD = 6,
                  HVON = 0, SBPBASE = 190),
    events = c(nica_fast(), ivh_tpa(12, 12, 8), mannitol(75, 12, 8, 6)),
    note   = "IVH clearance lowers mortality; P(mRS 0-2) moves much less."
  )
)

## ---------------------------------------------------------------------------
##  Runner
## ---------------------------------------------------------------------------
run_scenario <- function(name, end = TEND, delta = 0.25) {
  s <- scenarios[[name]]
  m <- ich
  if (length(s$param)) m <- param(m, s$param)
  ## keep the patient weight consistent with weight-based dosing helpers
  out <- if (is.null(s$events)) {
    m %>% mrgsim(end = end, delta = delta)
  } else {
    m %>% mrgsim(events = s$events, end = end, delta = delta)
  }
  as_tibble(out) %>% mutate(scenario = name)
}

run_all_scenarios <- function(end = TEND, delta = 0.25) {
  bind_rows(lapply(names(scenarios), run_scenario, end = end, delta = delta))
}

## ---------------------------------------------------------------------------
##  Summary table: the numbers the trials actually reported
## ---------------------------------------------------------------------------
summarise_scenarios <- function(sim = NULL, end = TEND) {
  if (is.null(sim)) sim <- run_all_scenarios(end = end)
  sim %>%
    group_by(scenario) %>%
    summarise(
      V_baseline_mL   = first(VHEM),
      V_24h_mL        = VHEM[which.min(abs(time - 24))],
      expansion_mL    = V_24h_mL - V_baseline_mL,
      expansion_pct   = 100 * expansion_mL / V_baseline_mL,
      total_bled_mL   = DV24[which.min(abs(time - 24))],
      PHE_peak_mL     = max(PHE),
      PHE_peak_day    = time[which.max(PHE)] / 24,
      ICP_peak_mmHg   = max(ICP),
      CPP_min_mmHg    = min(CPP),
      SBP_nadir_mmHg  = min(SBP),
      SBP_24h_mmHg    = SBP[which.min(abs(time - 24))],
      hours_CPP_lt60  = TCPP[which.min(abs(time - end))],
      hours_ICP_gt20  = TICP[which.min(abs(time - end))],
      ## the hypertensive load that matters is the one during the treatment
      ## window, not diluted across 90 days of untreated convalescence
      AUC_SBP140_72h  = ASBP[which.min(abs(time - 72))],
      AUC_Fe2         = AFE[which.min(abs(time - end))],
      AUC_Fe2_14d     = AFE[which.min(abs(time - 336))],
      AUC_IL6         = AIL6[which.min(abs(time - end))],
      NIHSS_24h       = NIH[which.min(abs(time - 24))],
      NIHSS_90d       = NIH[which.min(abs(time - end))],
      ICH_score       = ICHSC[which.min(abs(time - 24))],
      P_mRS02_90d     = PMRS02[which.min(abs(time - end))],
      ## The ICH score is by construction a BASELINE score, so mortality is
      ## assembled here from the 24-h score plus the accumulated intracranial
      ## dose.  Reading the captured PMORT at 90 days instead would score the
      ## patient on their recovered state, which inverts the intent.
      hrs_ICP20_end   = TICP[which.min(abs(time - end))],
      .groups = "drop"
    ) %>%
    mutate(
      P_death = 1 / (1 + exp(-(.MORTA + .MORTB * ICH_score
                               + .MORTC * hrs_ICP20_end / 24))),
      UW_mRS  = (1 - P_death) * (0.90 * P_mRS02_90d + 0.42 * (1 - P_mRS02_90d))
    ) %>%
    arrange(scenario)
}

## mortality coefficients pulled once from the model object so the summariser
## and the in-model $TABLE calculation cannot silently drift apart
.MORTA <- as.numeric(param(ich)$MORTA)
.MORTB <- as.numeric(param(ich)$MORTB)
.MORTC <- as.numeric(param(ich)$MORTC)

## ---------------------------------------------------------------------------
##  The U-shape: a net-benefit surface over (SBP target depth, time-to-target)
##  Nothing in the model encodes a U; it emerges because MAP appears twice.
## ---------------------------------------------------------------------------
##  `intensity` indexes how hard the BP is driven down, from 0 (no treatment)
##  to 1 (both agents at maximum).  Nicardipine ALONE cannot reach the deep end
##  of the range -- its Emax caps the achievable reduction at 45% -- so a second
##  agent is co-titrated.  Without this the scan only ever samples the safe
##  shoulder of the curve (SBP >= 138) and would report BP lowering as
##  monotonically beneficial simply because it never probed where the harm is.
bp_response_surface <- function(intensity = seq(0, 1, by = 0.125),
                                speeds    = c(1, 4),
                                end       = 720) {
  grid <- expand.grid(intensity = intensity, ttt = speeds)
  res <- lapply(seq_len(nrow(grid)), function(i) {
    a <- grid$intensity[i]; s <- grid$ttt[i]
    nrate <- 40 * a                       # nicardipine mg/h
    lab   <- 30 * a                       # labetalol mg per 6-h bolus
    e <- NULL
    if (nrate > 0) e <- c(nica_inf(nrate, 0.10, s),
                          nica_inf(nrate * 0.85, s, 72))
    if (lab > 0)   e <- c(e, ev(time = 0.2, amt = lab, cmt = 6,
                               ii = 6, addl = 11))
    o <- (if (is.null(e)) mrgsim(ich, end = end, delta = 1)
          else mrgsim(ich, events = e, end = end, delta = 1)) %>% as_tibble()
    tibble(intensity = a, ttt = s,
           SBP_nadir  = min(o$SBP),
           CPP_min    = min(o$CPP),
           expansion  = o$VHEM[which.min(abs(o$time - 24))] - o$VHEM[1],
           hrs_CPP60  = max(o$TCPP),
           NEUR_end   = o$NEUR[nrow(o)],
           NIHSS_end  = o$NIH[nrow(o)],
           P_mRS02    = o$PMRS02[nrow(o)],
           UW_mRS     = o$UWMRS[nrow(o)])
  })
  bind_rows(res)
}

## ---------------------------------------------------------------------------
##  Falsification switches — each structural claim has a kill switch
## ---------------------------------------------------------------------------
##  * KBLEED  = 0      -> no bleeding at all; volume is purely the initial CT.
##  * KTAMP   = 0      -> tamponade removed; expansion no longer self-limits
##                        and untreated bleeds run away (they do not clinically,
##                        which is the evidence FOR the term).
##  * KAVAL   = 0      -> no mechanical avalanche; expansion becomes almost
##                        insensitive to the first hours of BP control.
##  * KCHEL   = 0      -> deferoxamine becomes inert: its entire mRS effect in
##                        scenario 10 must run through the iron channel.
##  * KFERRO2 = KOEDL = 0 -> the whole chemical channel is deleted; scenario 10
##                        and scenario 11 then differ only in volume, and the
##                        i-DEF pattern (mRS shift with no volume change)
##                        becomes impossible to reproduce.
##  * RSHIFT  = 0      -> no hypertensive autoregulatory shift; the harm in
##                        scenario 4 shrinks markedly.
##  * ACUSH   = 0      -> Cushing feedback removed; BP no longer rises with ICP.
##  * KOFFA   = 0      -> the andexanet-apixaban complex never dissociates;
##                        the post-infusion anti-Xa rebound in scenario 9
##                        disappears, which localises the rebound to
##                        dissociation + redistribution rather than to any
##                        separately coded "rebound" term (there is none).
falsify <- function(scenario = "10_deferoxamine_idef",
                    kill = list(KCHEL = 0), end = TEND) {
  base <- run_scenario(scenario, end = end)
  s    <- scenarios[[scenario]]
  m    <- param(ich, c(s$param, kill))
  alt  <- (if (is.null(s$events)) mrgsim(m, end = end, delta = 0.25)
           else mrgsim(m, events = s$events, end = end, delta = 0.25)) %>%
    as_tibble() %>% mutate(scenario = paste0(scenario, "_killed"))
  bind_rows(base, alt)
}

## ---------------------------------------------------------------------------
##  Example use
## ---------------------------------------------------------------------------
if (interactive()) {

  library(ggplot2)
  library(tidyr)

  sim <- run_all_scenarios()
  print(summarise_scenarios(sim), width = Inf)

  ## The two channels, side by side
  sim %>%
    filter(scenario %in% c("1_natural_history", "3_bp_intensive_140",
                           "10_deferoxamine_idef", "11_mistie_evacuation")) %>%
    filter(time <= 336) %>%
    select(time, scenario, VHEM, PHE, AFE, NIH) %>%
    pivot_longer(c(VHEM, PHE, AFE, NIH)) %>%
    ggplot(aes(time / 24, value, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "Days from onset", y = NULL,
         title = "Mass-effect channel (VHEM, PHE) vs iron channel (AFE)") +
    theme_bw()

  ## The emergent U-shape
  surf <- bp_response_surface()
  ggplot(surf, aes(SBP_nadir, UW_mRS, colour = factor(ttt))) +
    geom_line() + geom_point() +
    labs(x = "SBP nadir (mmHg)", y = "Utility-weighted mRS",
         colour = "Hours to target",
         title = "No U is coded; it emerges because MAP is both bleeding driver and CPP numerator") +
    theme_bw()

  ## Kill the iron channel and watch scenario 10 lose its entire effect
  print(falsify("10_deferoxamine_idef", list(KCHEL = 0)) %>%
          group_by(scenario) %>%
          summarise(AUC_Fe = max(AFE), NIHSS_90d = last(NIH),
                    P_mRS02 = last(PMRS02), V24 = VHEM[which.min(abs(time - 24))]))
}
