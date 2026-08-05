## ===========================================================================
##  osa_mrgsolve_model.R
##  Osteosarcoma (OSA) — Quantitative Systems Pharmacology model
##  골육종 정량적 시스템 약리학 모델
##
##  ORGANISING THESIS
##  -----------------
##  Cure in osteosarcoma is a POISSON BET ON LESIONS NOBODY CAN SEE, and it is
##  paid for out of an EXPOSURE BUDGET with three hard organ ceilings.
##
##      P(cure) = exp( -lambda0 * [ 1 - exp( -n0 * e^-K ) ] )
##                 \________/       \_______________________/
##                  how many          did THIS lesion die
##                  lesions           (n0 ~ 1e6 clonogens)
##
##      K = kill integral,  subject to  E_DOX <= heart,
##                                      E_CIS <= cochlea + tubule,
##                                      E_MTX <= tubule
##
##  lambda0 = 1.80 is fixed by ONE historical number: cure after amputation
##  alone is exp(-lambda0) = 0.165, and the reported figure is 0.15-0.20.
##  Everything else follows from that anchor.  To reach the ~0.60 that MAP
##  achieves, K must be 6.6 log10 — and the model delivers 6.61.
##
##  And the methotrexate ceiling is not a constant, because THE DRUG IS
##  CLEARED THROUGH THE ORGAN IT DESTROYS:
##
##      C_tub = CL_ren(KID) * C_plasma / UrineFlow(KID)
##      S(pH) = 0.86 * 10^(0.682*(pH - 5))  mM
##
##  Everything above S crystallises, obstructs, and lowers BOTH CL_ren and
##  UrineFlow — which raises C_tub.  Positive feedback, so a critical urine
##  pH exists.  It is 7.30 for supersaturation, 7.15 for delayed elimination
##  and 6.95 for AKI, i.e. the protocol instruction "alkalinise to pH >= 7.0"
##  sits 0.30 units INSIDE the unstable region.
##
##  SEVEN RESULTS THAT ARE ARITHMETIC RATHER THAN ASSERTED
##  ------------------------------------------------------
##   A. A bifurcation, not a guideline: critical urine pH 7.30 / 7.15 / 6.95.
##   B. One pH unit is worth 4.81x the fluid, or a 4.81-fold dose reduction,
##      because C_tub scales as dose/flow and S scales as 10^(0.682 pH).
##   C. The loop gain rises as the kidney falls: halving the starting GFR
##      multiplies the same course's AUC by 1.59x and C48 by 3.8x, leaving the
##      kidney at 53.7% instead of 62.9%.  Twelve methotrexate and four
##      cisplatin courses are multiplicative, not additive.
##   D. Glucarpidase is triggered by a 48-h concentration but the AUC is
##      front-loaded: it removes 29% at 12 h, 10% at 24 h, 2% at the standard
##      48-h trigger, nothing at 72 h, and none of the tubular injury ever.
##   E. Huvos necrosis and micrometastatic kill are two integrals of the same
##      parameter (r = 0.97) — hence prognostic — but one integrates over the
##      primary and one over the lung, so re-timing therapy moves the marker
##      without moving the patient.  Prognostic AND manipulable.
##   F. MAPIE's null result is the arithmetic of its design: every kill term
##      carries (1 - RES), a poor responder is a patient in whom that factor
##      is small, and the escalation buys 5 IE cycles by giving up one
##      doxorubicin and one cisplatin.  Predicted EFS HR 0.98 in poor
##      responders, 0.86 in good ones.  EURAMOS-1: 0.98 (0.78-1.23).
##   G. Three plateau mechanisms, none of them the drugs.
##      DETECTABILITY: +0.2 log is worth +4 points of cure and needs 2764 per
##        arm; EURAMOS-1 randomised 618.  The same half-log is worth +19.6
##        points at CV 2% and +3.5 at CV 65%.
##      INDIFFERENCE: under-hydration buys +0.4 log10 of kill and pays 12%
##        treatment-related mortality, for survival 0.595 against MAP's 0.596.
##        Only dexrazoxane moves the wall (dox 600 + dexrazoxane: 0.649).
##      WRONG RECIPIENT: the same total extra log-kill is worth +0.116 in the
##        SECOND quartile of achieved exposure, +0.048 spread over everyone,
##        and +0.030 in the worst quartile.  Not the mean, not the tail — the
##        shoulder, at ~85-89% necrosis.
##
##  MODEL SIZE
##  ----------
##    55 ODE compartments · 141 parameters · 20 therapy scenarios ·
##    protocol go/no-go gate layer · virtual-population cure layer
##
##  UNITS
##  -----
##    time            hours (the 4-h infusion and the +24 h rescue matter;
##                    the disease runs 40 weeks = 6720 h)
##    methotrexate    mmol (amount) / mM (concentration); 1 mM = 1000 uM
##    doxorubicin,
##    cisplatin,
##    ifosfamide,
##    etoposide       mg / mg/L
##    tumour          mL (1 mL ~ 1 g ~ 1e9 cells)
##    LMET            natural log of clonogens in one micrometastatic lesion
##    ANC             10^9/L
##    organ states    fraction of baseline (KID, CMV) — 1.0 is intact
##
##  CALIBRATION (every target is a published number; the achieved value is
##  from osa_reference_model.py / osa_reference_output.txt, which integrates
##  these same equations in Python because R is unavailable in the build
##  container)
##  --------------------------------------------------------------------------
##    HDMTX 12 g/m2 monitoring   C24 8.8 uM, C48 0.56, C72 0.042
##                               (protocol thresholds 10 / 1 / 0.1)  [Evans]
##    MTX solubility vs pH       0.39 / 1.55 / 9.04 mg/mL at pH 5 / 6 / 7
##    cure, surgery alone        0.168   (historical 0.15-0.20)
##    cure, MAP                  0.596   (EURAMOS-1 5-y EFS 0.54-0.59)
##    cure, metastatic at dx     0.322   (reported 0.20-0.30)
##    Huvos necrosis on MAP      90.3%   (good responders ~45-50% of patients;
##                                        this run is the RES0 = 0.02 tumour)
##    doxorubicin cumulative     450 mg/m2 (EURAMOS-1 MAP)
##    cisplatin cumulative       480 mg/m2 (EURAMOS-1 MAP)
##    LVEF after MAP             58.5% from 62.0
##    hearing loss after MAP     12.2 dB high-frequency shift
##    ANC nadir on MAP           0.39 x10^9/L (grade 4)
##    courses on schedule        MAP 100%, MAPIE 76%
##                               (EURAMOS-1 completed protocol 76% / 51%)
##    treatment-related death    MAP 1.0%, MAPIE 3.0% (reported ~1% overall)
##
##  ASSUMPTIONS, FLAGGED
##  --------------------
##    * KPPT, KDIS, PPT50, KINJ_PPT and KINJ_TUB have no direct in-vivo
##      measurement.  They are ASSUMED, and constrained only by the
##      requirement that the derived critical urine pH land near the
##      guideline value, which it does (7.30 vs "keep pH >= 7.0").  The
##      SHAPE of the result — a bifurcation with a 4.81x pH/fluid exchange
##      rate — comes from the solubility law and is not fitted.
##    * A single scalar RES stands for cross-resistance across all four
##      cytotoxics.  This is deliberately conservative: it makes escalation
##      look WORSE than a drug-specific resistance model would, and it is
##      the reason the MAPIE prediction lands where the trial did.
##    * LUNG_IMM = 4 (immune kill weighted 4x higher in lung than primary)
##      is ASSUMED, and it is the sole channel through which mifamurtide
##      acts.  The predicted mifamurtide benefit (+1.4 points of survival)
##      is therefore a consequence of that assumption, not evidence for it.
##    * DEXR and MIFA are EFFECT-SITE surrogates, not plasma drug: the
##      turnover constants are the turnover of the effect (iron chelation /
##      TOP2B depletion; macrophage activation), which far outlasts the
##      parent molecule.
##    * TRM is a hazard integral with an infection arm written as the PRODUCT
##      of neutropenia depth and mucositis grade (both are needed for a
##      barrier-breach death) plus renal, cardiac and baseline arms.
##
##  This is an educational / research model.  It is not validated for
##  clinical use, dosing, or regulatory submission.
## ===========================================================================

library(mrgsolve)
suppressMessages(library(dplyr))

osa_code <- '
$PARAM @annotated
// ---------------- patient ------------------------------------------------
BSA      :  1.7   : Body surface area (m2)
MTX_MW   : 454.44 : Methotrexate molecular weight (g/mol)

// ---------------- methotrexate PK (3 compartments, L and L/h) ------------
V1       : 18.0   : Central volume (L)
V2       : 10.0   : Shallow peripheral volume (L)
V3       :  5.0   : Deep / third-space volume - effusion, ascites (L)
Q2       : 18.0   : Central-shallow intercompartmental clearance (L/h)
Q3       :  0.60  : Central-deep intercompartmental clearance (L/h)
CLREN0   :  7.5   : Renal clearance at intact kidney (L/h) - 86% of total
CLNR     :  1.2   : Non-renal clearance - 7-OH-MTX, biliary (L/h)
FU       :  0.55  : Unbound fraction
GFR0     :  7.2   : GFR at intact kidney (L/h) = 120 mL/min

// ---------------- the renal feedback loop --------------------------------
HYDR     :  3.0   : Protocol hydration (L/m2/day) - reference only
UF0      :  0.21  : Urine flow at 3 L/m2/day, 1.7 m2 (L/h)
URINE_PH :  7.5   : Urine pH under bicarbonate
SOL_A    :  0.86  : Methotrexate solubility at pH 5.0 (mM)
SOL_B    :  0.682 : Decades of solubility gained per pH unit
KPPT     :  0.45  : Crystallisation rate constant (1/h per mM excess) ASSUMED
KDIS     :  0.030 : Crystal redissolution rate (1/h) ASSUMED
PPT50    :  2.0   : Crystal burden at half-maximal obstruction (mmol) ASSUMED
KINJ_PPT :  0.020 : Obstructive tubular injury rate (1/h) ASSUMED
KINJ_TUB :  0.0016: Direct concentration-driven tubular injury (1/h) ASSUMED
KTUB     : 12.0   : Half-maximal concentration for direct injury (mM) ASSUMED
KREP_KID :  0.0022: Tubular regeneration rate (1/h)

// ---------------- folate, polyglutamates and rescue ---------------------
KPG      :  0.055 : Tumour polyglutamation rate
MPGMAX   :  8.0   : Maximum tumour polyglutamate pool
KDPG     :  0.010 : Tumour polyglutamate loss rate (1/h)
KPGM     :  0.035 : Marrow / mucosa polyglutamation rate
MPGMAXM  :  6.0   : Maximum marrow polyglutamate pool
KDPGM    :  0.032 : Marrow polyglutamate loss rate (1/h)
KIPG     :  1.10  : Polyglutamate giving half-maximal TS / DHFR block
CL_LV    :  9.0   : Leucovorin clearance (L/h)
V_LV     : 25.0   : Leucovorin volume (L)
KRF      :  0.60  : Reduced folate that doubles the apparent Ki
KRF_IN   :  0.020 : Reduced-folate uptake, tumour
KRF_INM  :  0.055 : Reduced-folate uptake, marrow (2.75x the tumour)
KRF_OUT  :  0.045 : Reduced-folate turnover (1/h)
RF0      :  0.15  : Baseline reduced folate pool
CL_GLU0  : 90.0   : Glucarpidase-mediated clearance at full effect (L/h)
KDGLU    :  0.35  : Glucarpidase effect decay (1/h)

// ---------------- doxorubicin and the heart -----------------------------
CL_DOX   : 45.0   : Doxorubicin clearance (L/h)
V_DOXC   : 25.0   : Doxorubicin central volume (L)
V_DOXP   : 1200   : Doxorubicin deep tissue volume (L)
Q_DOX    : 60.0   : Doxorubicin intercompartmental clearance (L/h)
KMET_DOXOL: 0.22  : Carbonyl reduction to doxorubicinol (1/h)
CL_DOXOL : 30.0   : Doxorubicinol clearance (L/h)
KCM      :  0.0235: Cardiomyocyte loss per (mg/L)-h of doxorubicinol
KCM_REP  :  1.2e-6: Cardiomyocyte repair - essentially irreversible
DEXR_EMAX:  0.72  : Maximum dexrazoxane protection
DEXR_EC50:  3.0   : Dexrazoxane effect-site EC50 (mg/L)
CL_DEXR  :  3.0   : Dexrazoxane EFFECT turnover (L/h) - not plasma drug
V_DEXR   : 40.0   : Dexrazoxane effect volume (L)

// ---------------- cisplatin, tubule and cochlea -------------------------
CL_CIS   : 30.0   : Cisplatin clearance (L/h)
V_CIS    : 18.0   : Cisplatin volume (L)
KADD_T   :  0.030 : Pt-DNA adduct formation, tumour
KREP_ADDT:  0.020 : Adduct repair, tumour (ERCC1/XPF)
KADD_K   :  0.020 : Pt-DNA adduct formation, tubule (OCT2-mediated)
KREP_ADDK:  0.008 : Adduct repair, tubule
KINJ_CIS :  0.0060: Cisplatin tubular injury rate (1/h)
KADD_C   :  0.014 : Pt-DNA adduct formation, cochlea
KREP_ADDC:  0.0016: Adduct repair, cochlea - slower than therapy
KHL      :  0.055 : dB of high-frequency shift per adduct-hour
KMG      :  0.010 : Renal magnesium wasting rate
MG0      :  0.85  : Baseline serum magnesium (mmol/L)

// ---------------- ifosfamide and etoposide -----------------------------
CL_IFO   : 12.0   : Ifosfamide clearance (L/h)
V_IFO    : 40.0   : Ifosfamide volume (L)
FACT_IFO :  0.55  : Fraction 4-hydroxylated by CYP3A4/2B6
CL_ETO   :  1.9   : Etoposide clearance (L/h)
V_ETO    : 18.0   : Etoposide volume (L)
KINJ_IFO :  0.0022: Ifosfamide / chloroacetaldehyde tubular injury (1/h)

// ---------------- tumour growth and kill --------------------------------
PRIM0    : 150.0  : Primary tumour volume at diagnosis (mL)
PRIMMAX  : 1500   : Gompertz carrying capacity (mL)
KG       :  9.63e-4: Primary growth rate (1/h) - volume doubling ~30 d
KMTX     :  0.1828: Maximum methotrexate kill rate at full TS block (1/h)
KDOX     :  0.2766: Maximum doxorubicin kill rate (1/h)
EC50_DOX :  0.35  : Doxorubicin EC50 (mg/L)
KCIS     :  0.1730: Maximum cisplatin kill rate (1/h)
EC50_CIS :  0.55  : Cisplatin adduct EC50
KIFO     :  0.0200: Maximum ifosfamide kill rate (1/h)
EC50_IFO :  6.0   : Ifosfamide EC50 (mg/L)
KETO     :  0.0160: Maximum etoposide kill rate (1/h)
EC50_ETO :  2.5   : Etoposide EC50 (mg/L)
KIMM     :  6.0e-5: Immune kill rate coefficient
LUNG_IMM :  4.0   : Immune kill weighting, lung vs primary  ASSUMED
KNECCLR  :  6.0e-5: Necrotic osteoid clearance - it is RETAINED (1/h)
PEN_MIN  :  0.30  : Minimum drug penetration into the mineralised primary

// ---------------- micrometastatic pool ----------------------------------
N_LESION_CELLS: 1.0e6 : Clonogens in one micrometastatic lesion
LAMBDA0  :  1.80  : Poisson mean occult lesions, localised presentation
LAMBDA0M : 15.0   : Poisson mean, overt metastatic presentation
KGM      :  3.4e-4: Dormant micrometastatic growth rate (1/h) ~85 d doubling
LMAX     : 25.3   : ln(1e11) - carrying capacity of one lung lesion

// ---------------- resistance -------------------------------------------
RES0     :  0.020 : Resistant clone fraction at diagnosis
KMUT     :  1.2e-5: Spontaneous resistance acquisition (1/h)
KSEL     :  0.220 : Selection coefficient under cytotoxic pressure

// ---------------- bone remodelling vicious cycle -----------------------
OCL0     :  1.0   : Baseline osteoclast activity
KOCL     :  0.010 : Osteoclast turnover (1/h)
RKL0     :  1.0   : Baseline RANKL
KRKL_T   :  0.0035: Tumour-driven RANKL production
KDRKL    :  0.020 : RANKL degradation (1/h)
OPG      :  1.0   : Osteoprotegerin (decoy receptor)
KTGFM    :  0.014 : Matrix TGF-beta / IGF-1 liberation per osteoclast-hour
KDTGFM   :  0.030 : Matrix growth-factor clearance (1/h)
TGF_BOOST:  0.45  : Fractional growth boost per unit matrix TGF-beta
KDENO    :  9.0e-4: Denosumab elimination (1/h) - t1/2 ~32 d
V_DENO   :  3.5   : Denosumab volume (L)
DENO_IC50:  1.4   : Denosumab RANKL-neutralising IC50 (mg/L)
KDZOL    :  1.2e-4: Zoledronate bone-bound loss (1/h)
ZOL_IC50 :  0.30  : Zoledronate osteoclast IC50

// ---------------- vasculature, hypoxia, immunity ------------------------
VASC0    :  0.55  : Baseline functional vascular density
KVASC    :  0.004 : Angiogenic drive
KDVASC   :  0.004 : Vessel regression
HYP_K    :  0.6   : Hypoxia gain on vascular deficit
CTL0     :  1.0   : Baseline CD8 T-cell activity
KCTL     :  0.004 : CTL recruitment
KDCTL    :  0.006 : CTL loss
M20      :  1.0   : Baseline M2 / TAM
KM2      :  0.003 : M2 recruitment
KDM2     :  0.004 : M2 loss
CL_MIFA  :  0.12  : Mifamurtide EFFECT turnover (L/h) - not plasma drug
V_MIFA   : 12.0   : Mifamurtide effect volume (L)
MIFA_EMAX:  0.55  : Maximum mifamurtide macrophage-activation effect
MIFA_EC50:  0.25  : Mifamurtide effect EC50 (mg/L)
VEGFI    :  0.0   : Anti-angiogenic drive (regorafenib etc.), 0 = none

// ---------------- myelosuppression (Friberg transit model) -------------
CIRC0    :  4.0   : Baseline ANC (10^9/L)
MTT      : 125.0  : Mean maturation transit time (h)
GAM      :  0.17  : Feedback exponent
SL_MTX   :  5.00  : Methotrexate slope on proliferating marrow
SL_DOX   : 34.0   : Doxorubicin slope (per mg/L)
SL_CIS   :  4.00  : Cisplatin slope (per mg/L)
SL_IFO   :  0.60  : Ifosfamide slope (per mg/L)
SL_ETO   :  4.60  : Etoposide slope (per mg/L)
GCSF     :  0.0   : G-CSF support (fractional feedback boost)

// ---------------- mucositis and biomarkers ------------------------------
KMUC     :  0.500 : Mucositis formation rate
KMUC_REP :  0.020 : Mucosal healing rate (1/h)
ALP0     :  1.0   : Baseline alkaline phosphatase (x ULN)
KALP     :  0.010 : ALP production per unit tumour x osteoclast
KDALP    :  0.012 : ALP clearance (1/h)

// ---------------- treatment-related mortality --------------------------
KTRM_BASE:  1.40e-6: Baseline hazard - line sepsis, thrombosis, surgery
KTRM_INF :  6.90e-4: Infection hazard coefficient (neutropenia x mucositis)
KTRM_REN :  7.00e-4: Renal-failure hazard coefficient
KTRM_CARD:  2.20e-4: Cardiac hazard coefficient
ANC_G4   :  0.35  : ANC scale of the infection hazard (10^9/L)
KID_CRIT :  0.60  : Kidney fraction below which the renal hazard opens

// ---------------- population layer -------------------------------------
POP_CV_K :  0.25  : Between-patient CV of achieved log-kill

$CMT @annotated
// ---- methotrexate PK and the renal feedback loop -----------------------
A1       : Methotrexate, central (mmol)
A2       : Methotrexate, shallow peripheral (mmol)
A3       : Methotrexate, deep / third space (mmol)
PPT      : Intratubular methotrexate crystal burden (mmol)
KID      : Kidney functional mass (fraction of baseline)
// ---- folate and rescue ------------------------------------------------
MPGT     : Methotrexate polyglutamate, tumour
MPGM     : Methotrexate polyglutamate, marrow and mucosa
LV       : Leucovorin, central (arbitrary mass)
RFT      : Reduced folate pool, tumour
RFM      : Reduced folate pool, marrow and mucosa
GLUC     : Glucarpidase effect
// ---- doxorubicin and the heart ----------------------------------------
DOXC     : Doxorubicin, central (mg)
DOXP     : Doxorubicin, deep tissue (mg)
DOXOL    : Doxorubicinol (mg)
CMV      : Cardiomyocyte viability (fraction)
CUMDOX   : Cumulative doxorubicin dose (mg/m2)
DEXR     : Dexrazoxane effect site (mg)
// ---- cisplatin, tubule and cochlea ------------------------------------
CISC     : Cisplatin, central (mg)
ADDT     : Pt-DNA adducts, tumour
ADDK     : Pt-DNA adducts, renal tubule
ADDC     : Pt-DNA adducts, cochlea
HL       : Cumulative high-frequency hearing loss (dB)
MG       : Serum magnesium (mmol/L)
// ---- ifosfamide and etoposide ----------------------------------------
IFOA     : 4-OH-ifosfamide / isophosphoramide mustard (mg)
ETOC     : Etoposide, central (mg)
// ---- tumour -----------------------------------------------------------
PRIM     : Viable primary tumour volume (mL)
NEC      : Necrotic primary tumour volume (mL)
RES      : Resistant clone fraction
LMET     : ln(clonogens in one micrometastatic lesion)
KI_MTX   : Micrometastatic kill attributed to methotrexate (nats)
KI_DOX   : Micrometastatic kill attributed to doxorubicin (nats)
KI_CIS   : Micrometastatic kill attributed to cisplatin (nats)
KI_IE    : Micrometastatic kill attributed to ifosfamide/etoposide (nats)
KI_IMM   : Micrometastatic kill attributed to immune clearance (nats)
// ---- bone remodelling -------------------------------------------------
OCL      : Osteoclast activity
TGFM     : Matrix-liberated TGF-beta / IGF-1 / BMP-2
RKL      : RANKL
DENO     : Denosumab (mg)
ZOL      : Zoledronate, bone-bound (mg)
// ---- vasculature and hypoxia -----------------------------------------
VASC     : Functional vascular density
HYP      : Hypoxic fraction
// ---- immunity ---------------------------------------------------------
CTL      : CD8 cytotoxic T-cell activity
M2       : M2 / tumour-associated macrophage
MIFA     : Mifamurtide effect site (mg)
// ---- myelosuppression -------------------------------------------------
PROL     : Proliferating marrow cells (10^9/L equivalent)
TR1      : Maturation transit 1
TR2      : Maturation transit 2
TR3      : Maturation transit 3
CIRC     : Circulating absolute neutrophil count (10^9/L)
// ---- mucosa, exposure accumulators, biomarkers, mortality ------------
MUC      : Mucositis grade (0-4)
EMTX     : Methotrexate exposure accumulator (mM*h)
EDOX     : Doxorubicin exposure accumulator (mg/L*h)
ECIS     : Cisplatin exposure accumulator (mg/L*h)
EIFO     : Ifosfamide exposure accumulator (mg/L*h)
ALP      : Alkaline phosphatase (x ULN)
TRM      : Cumulative treatment-related mortality hazard

$MAIN
// ---- initial conditions ------------------------------------------------
KID_0  = 1.0;
RFT_0  = RF0;
RFM_0  = RF0;
CMV_0  = 1.0;
MG_0   = MG0;
PRIM_0 = PRIM0;
RES_0  = RES0;
OCL_0  = OCL0;
RKL_0  = RKL0;
TGFM_0 = 0.4;
VASC_0 = VASC0;
HYP_0  = 0.3;
CTL_0  = CTL0;
M2_0   = M20;
PROL_0 = CIRC0;
TR1_0  = CIRC0;
TR2_0  = CIRC0;
TR3_0  = CIRC0;
CIRC_0 = CIRC0;
ALP_0  = ALP0;
LMET_0 = log(N_LESION_CELLS);

$ODE
// =========================================================================
//  1. METHOTREXATE PK AND THE RENAL POSITIVE FEEDBACK LOOP
//
//  This block is the heart of the model.  Read it as one statement:
//  the concentration that damages the tubule is the excretion RATE divided
//  by the urine FLOW, and injury lowers both the rate constant and the flow.
// =========================================================================
double KIDc = KID;
if (KIDc < 0.02) KIDc = 0.02;
if (KIDc > 1.00) KIDc = 1.00;

double C1  = A1 / V1;                       // mM
double C2  = A2 / V2;
double C3  = A3 / V3;

double OBS   = PPT / (PPT + PPT50);         // fractional obstruction
double CLREN = CLREN0 * KIDc * (1.0 - 0.85 * OBS);
double UF    = UF0 * KIDc * (1.0 - 0.90 * OBS);
if (UF < 0.02) UF = 0.02;

double excr = CLREN * C1;                   // mmol/h into the tubular lumen
double CTUB = excr / UF;                    // mM in tubular fluid
double SOL  = SOL_A * pow(10.0, SOL_B * (URINE_PH - 5.0));
double excess = CTUB - SOL;
if (excess < 0.0) excess = 0.0;

double J_ppt = KPPT * excess * UF;                       // mmol/h -> crystal
double J_dis = KDIS * PPT * (SOL > CTUB ? (SOL - CTUB) / SOL : 0.0);
double CL_GLU = CL_GLU0 * GLUC / (GLUC + 1.0);

dxdt_A1  = -(CLNR + CLREN + CL_GLU) * C1
           - Q2 * (C1 - C2) - Q3 * (C1 - C3) - J_ppt;
dxdt_A2  = Q2 * (C1 - C2);
dxdt_A3  = Q3 * (C1 - C3);
dxdt_PPT = J_ppt - J_dis;
dxdt_GLUC = -KDGLU * GLUC;

// ---- the shared state variable: kidney functional mass ------------------
double inj = KINJ_PPT * OBS
           + KINJ_TUB * CTUB / (CTUB + KTUB)
           + KINJ_CIS * ADDK / (ADDK + 1.0)
           + KINJ_IFO * IFOA / (IFOA + 5.0);
dxdt_KID = -inj * KIDc + KREP_KID * (1.0 - KIDc);

// =========================================================================
//  2. LEUCOVORIN RESCUE — ONE CLOCK FOR TUMOUR AND MARROW ALIKE
// =========================================================================
double CLV = LV / V_LV;
dxdt_LV  = -CL_LV * CLV;
dxdt_RFT = KRF_IN  * CLV - KRF_OUT * (RFT - RF0);
dxdt_RFM = KRF_INM * CLV - KRF_OUT * (RFM - RF0);

// penetration penalty applies to the mineralised PRIMARY only
double PEN = PEN_MIN + (1.0 - PEN_MIN) * VASC / (VASC + VASC0);
double C_tum = C1 * PEN;

dxdt_MPGT = KPG  * C_tum * (1.0 - MPGT / MPGMAX)  - KDPG  * MPGT;
dxdt_MPGM = KPGM * C1    * (1.0 - MPGM / MPGMAXM) - KDPGM * MPGM;

// fractional TS / DHFR block; leucovorin RAISES the apparent Ki in BOTH
double KiT = KIPG * (1.0 + RFT / KRF);
double KiM = KIPG * (1.0 + RFM / KRF);
double blockT = MPGT / (MPGT + KiT);
double blockM = MPGM / (MPGM + KiM);
// a fully perfused micrometastasis sees a higher polyglutamate load
double mpgL = MPGT * 1.25;
if (mpgL > MPGMAX) mpgL = MPGMAX;
double blockL = mpgL / (mpgL + KiT);

// =========================================================================
//  3. DOXORUBICIN AND THE CARDIAC CEILING
// =========================================================================
double CD  = DOXC / V_DOXC;
double CDP = DOXP / V_DOXP;
dxdt_DOXC  = -CL_DOX * CD - Q_DOX * (CD - CDP) - KMET_DOXOL * DOXC;
dxdt_DOXP  = Q_DOX * (CD - CDP);
dxdt_DOXOL = KMET_DOXOL * DOXC - CL_DOXOL * DOXOL / V_DOXC;
dxdt_CUMDOX = 0.0;                          // incremented by dosing events

double CDOL = DOXOL / V_DOXC;
double CDX  = DEXR / V_DEXR;
double prot = 1.0 - DEXR_EMAX * CDX / (CDX + DEXR_EC50);
dxdt_DEXR = -CL_DEXR * CDX;
dxdt_CMV  = -KCM * (CDOL + 0.25 * CD) * prot * CMV + KCM_REP * (1.0 - CMV);

// =========================================================================
//  4. CISPLATIN — TUMOUR, TUBULE AND COCHLEA FROM ONE CONCENTRATION
// =========================================================================
double CC = CISC / V_CIS;
dxdt_CISC = -CL_CIS * CC;
dxdt_ADDT = KADD_T * CC * PEN - KREP_ADDT * ADDT;
dxdt_ADDK = KADD_K * CC       - KREP_ADDK * ADDK;
dxdt_ADDC = KADD_C * CC       - KREP_ADDC * ADDC;
dxdt_HL   = KHL * ADDC / (ADDC + 1.0);
dxdt_MG   = -KMG * ADDK / (ADDK + 1.0) + 0.004 * (MG0 - MG);

// =========================================================================
//  5. IFOSFAMIDE AND ETOPOSIDE (the MAPIE escalation)
// =========================================================================
double CI = IFOA / V_IFO;
double CE = ETOC / V_ETO;
dxdt_IFOA = -CL_IFO * CI;
dxdt_ETOC = -CL_ETO * CE;

// =========================================================================
//  6. BONE REMODELLING VICIOUS CYCLE
//     Note where TGFM lands: it multiplies GROWTH, not survival.  That one
//     structural choice is why zoledronate and denosumab move osteolysis and
//     ALP without moving the cure fraction.
// =========================================================================
double CDEN = DENO / V_DENO;
dxdt_DENO = -KDENO * DENO;
dxdt_ZOL  = -KDZOL * ZOL;
double rkl_free  = RKL / (1.0 + CDEN / DENO_IC50);
double ocl_drive = rkl_free / (rkl_free + OPG);
double zol_inh   = 1.0 / (1.0 + ZOL / ZOL_IC50);
dxdt_RKL  = KRKL_T * PRIM / 100.0 - KDRKL * RKL;
dxdt_OCL  = KOCL * (2.0 * ocl_drive * zol_inh - OCL);
dxdt_TGFM = KTGFM * OCL - KDTGFM * TGFM;

// =========================================================================
//  7. VASCULATURE, HYPOXIA, IMMUNITY
// =========================================================================
dxdt_VASC = KVASC * PRIM / (PRIM + 200.0) - KDVASC * VASC * (1.0 + VEGFI);
dxdt_HYP  = 0.01 * (HYP_K * (1.0 - VASC) - HYP);

double CMF = MIFA / V_MIFA;
dxdt_MIFA = -CL_MIFA * CMF;
double mifa_e = MIFA_EMAX * CMF / (CMF + MIFA_EC50);
dxdt_M2  = KM2 * (1.0 + 0.6 * HYP) - KDM2 * M2 * (1.0 + 2.0 * mifa_e);
dxdt_CTL = KCTL * (1.0 + 1.5 * mifa_e) - KDCTL * CTL * (1.0 + 0.8 * M2);

// =========================================================================
//  8. TUMOUR KILL — TWO COMPARTMENTS, TWO DIFFERENT INTEGRALS
//     The primary carries the penetration penalty and the matrix TGF-beta
//     growth boost.  The micrometastases carry neither, and they carry the
//     immune term at LUNG_IMM-fold weight.  Necrosis and cure are therefore
//     integrals of the same parameter over DIFFERENT operators, which is
//     what makes Huvos grade prognostic and manipulable at the same time.
// =========================================================================
double res = RES;
if (res > 0.95) res = 0.95;
double kill_imm = KIMM * CTL / (1.0 + M2);

double kill_mtx = KMTX * blockT * (1.0 - res);
double kill_dox = KDOX * (CD * PEN) / (CD * PEN + EC50_DOX) * (1.0 - res);
double kill_cis = KCIS * ADDT / (ADDT + EC50_CIS) * (1.0 - res);
double kill_ifo = KIFO * (CI * PEN) / (CI * PEN + EC50_IFO) * (1.0 - res);
double kill_eto = KETO * (CE * PEN) / (CE * PEN + EC50_ETO) * (1.0 - res);
double kill_prim = kill_mtx + kill_dox + kill_cis + kill_ifo + kill_eto
                 + kill_imm;

double gomp = log(PRIMMAX / (PRIM > 1e-6 ? PRIM : 1e-6));
if (gomp < 0.0) gomp = 0.0;
double grow = KG * (1.0 + TGF_BOOST * TGFM) * gomp;
dxdt_PRIM = PRIM * (grow - kill_prim);
dxdt_NEC  = PRIM * kill_prim - KNECCLR * NEC;

// ---- micrometastatic pool: LMET = ln(clonogens per lesion) --------------
double km_mtx = KMTX * blockL * (1.0 - res);
double km_dox = KDOX * CD / (CD + EC50_DOX) * (1.0 - res);
double km_cis = KCIS * ADDT / (ADDT + EC50_CIS) * (1.0 - res);
double km_ie  = (KIFO * CI / (CI + EC50_IFO)
               + KETO * CE / (CE + EC50_ETO)) * (1.0 - res);
double km_imm = LUNG_IMM * kill_imm;
double kill_met = km_mtx + km_dox + km_cis + km_ie + km_imm;

double gm = 1.0 - LMET / LMAX;
if (gm < 0.0) gm = 0.0;
dxdt_LMET = (LMET > -25.0) ? (KGM * gm - kill_met) : 0.0;
dxdt_KI_MTX = km_mtx;
dxdt_KI_DOX = km_dox;
dxdt_KI_CIS = km_cis;
dxdt_KI_IE  = km_ie;
dxdt_KI_IMM = km_imm;

// ---- resistance: acquisition plus selection under pressure -------------
dxdt_RES = KMUT * (1.0 - res) + KSEL * res * (1.0 - res) * kill_prim;

// =========================================================================
//  9. EXPOSURE ACCUMULATORS — THE BUDGET, MADE EXPLICIT
// =========================================================================
dxdt_EMTX = (C1 > 0.001) ? C1 : 0.0;
dxdt_EDOX = CD;
dxdt_ECIS = CC;
dxdt_EIFO = CI;

// =========================================================================
// 10. MYELOSUPPRESSION (Friberg) AND MUCOSITIS — THE GOVERNOR'S SENSORS
// =========================================================================
double ktr = 4.0 / MTT;
double edrug = SL_MTX * blockM + SL_DOX * CD + SL_CIS * CC
             + SL_IFO * CI + SL_ETO * CE;
if (edrug > 0.98) edrug = 0.98;
double circ_s = (CIRC > 0.05) ? CIRC : 0.05;
double fb = pow(CIRC0 / circ_s, GAM);
if (fb > 6.0) fb = 6.0;
fb = fb * (1.0 + GCSF);

dxdt_PROL = ktr * PROL * ((1.0 - edrug) * fb - 1.0);
dxdt_TR1  = ktr * (PROL - TR1);
dxdt_TR2  = ktr * (TR1  - TR2);
dxdt_TR3  = ktr * (TR2  - TR3);
dxdt_CIRC = ktr * (TR3  - CIRC);

double muc_head = 1.0 - MUC / 4.0;
if (muc_head < 0.0) muc_head = 0.0;
dxdt_MUC = KMUC * blockM * 4.0 * muc_head - KMUC_REP * MUC;

// =========================================================================
// 11. BIOMARKER AND TREATMENT-RELATED MORTALITY
//     The infection arm is the PRODUCT of neutropenia depth and mucositis
//     grade: a barrier breach with neutrophils is survivable and neutropenia
//     with an intact barrier is survivable, and the two together are not.
// =========================================================================
dxdt_ALP = KALP * (PRIM / 150.0) * (1.0 + OCL) - KDALP * ALP;

double g_neut = exp(-CIRC / ANC_G4);
double g_muc  = MUC / 4.0;
if (g_muc > 1.0) g_muc = 1.0;
double g_ren  = (KID_CRIT - KIDc) / KID_CRIT;
if (g_ren < 0.0) g_ren = 0.0;
double g_card = (0.85 - CMV) / 0.85;
if (g_card < 0.0) g_card = 0.0;
dxdt_TRM = KTRM_BASE + KTRM_INF * g_neut * g_muc
         + KTRM_REN * g_ren + KTRM_CARD * g_card;

$TABLE
// ---- observable concentrations -----------------------------------------
double CMTX_uM  = A1 / V1 * 1000.0;               // uM, the monitored value
double CDOXo    = DOXC / V_DOXC;
double CCISo    = CISC / V_CIS;
double CIFOo    = IFOA / V_IFO;
double CETOo    = ETOC / V_ETO;

// ---- the loop, made observable -----------------------------------------
double KIDo   = (KID < 0.02) ? 0.02 : ((KID > 1.0) ? 1.0 : KID);
double OBSo   = PPT / (PPT + PPT50);
double CLRENo = CLREN0 * KIDo * (1.0 - 0.85 * OBSo);
double UFo    = UF0 * KIDo * (1.0 - 0.90 * OBSo);
if (UFo < 0.02) UFo = 0.02;
double CTUBo  = CLRENo * (A1 / V1) / UFo;         // mM
double SOLo   = SOL_A * pow(10.0, SOL_B * (URINE_PH - 5.0));
double SSRAT  = CTUBo / SOLo;                     // supersaturation ratio
double eGFRo  = 120.0 * KIDo;                     // mL/min/1.73m2 equivalent

// ---- organ ceilings ----------------------------------------------------
double LVEFo  = 62.0 * (0.35 + 0.65 * pow(CMV > 0.0 ? CMV : 0.0, 0.6));
double HLo    = HL;
double MUCo   = MUC;

// ---- disease read-outs -------------------------------------------------
double TOTV   = PRIM + NEC;
double HUVOS  = (TOTV > 0.0) ? 100.0 * NEC / TOTV : 0.0;
double KNATS  = log(N_LESION_CELLS) - LMET;
if (KNATS < 0.0) KNATS = 0.0;
double KLOG10 = KNATS / log(10.0);

// ---- the double exponential -------------------------------------------
//  surviving clonogens in ONE lesion, then the chance that lesion lives,
//  then the chance that NONE of a Poisson(LAMBDA0) number of them lives
double SURVCELL = N_LESION_CELLS * exp(-KNATS);
double PLES     = 1.0 - exp(-SURVCELL);
double PCURE    = exp(-LAMBDA0 * PLES);
double PCUREM   = exp(-LAMBDA0M * PLES);
double PTRM     = 1.0 - exp(-TRM);
double POS      = PCURE * (1.0 - PTRM);           // overall survival

$CAPTURE @annotated
CMTX_uM : Plasma methotrexate (uM) - the monitored value
CDOXo   : Plasma doxorubicin (mg/L)
CCISo   : Plasma cisplatin (mg/L)
CIFOo   : 4-OH-ifosfamide (mg/L)
CETOo   : Plasma etoposide (mg/L)
CTUBo   : Tubular methotrexate concentration (mM)
SOLo    : Methotrexate solubility at this urine pH (mM)
SSRAT   : Supersaturation ratio C_tub / S - crystallises above 1
eGFRo   : eGFR (mL/min/1.73m2)
LVEFo   : Left ventricular ejection fraction (%)
HLo     : High-frequency hearing loss (dB)
MUCo    : Mucositis grade (0-4)
HUVOS   : Necrosis fraction of the primary (%)
KNATS   : Micrometastatic log-kill (nats)
KLOG10  : Micrometastatic log-kill (log10)
SURVCELL: Surviving clonogens per micrometastatic lesion
PLES    : Probability that one lesion survives
PCURE   : Cure probability, localised presentation
PCUREM  : Cure probability, overt metastatic presentation
PTRM    : Treatment-related mortality probability
POS     : Overall survival = cure x (1 - TRM)
'

mod <- mcode("osa", osa_code, soloc = tempdir())

## ===========================================================================
##  DOSING HELPERS — the EURAMOS-1 backbone, in hours
## ===========================================================================
W   <- 168                                   # hours per week
BSA <- 1.7

## MAP (EURAMOS-1 arm A): doxorubicin 75 mg/m2 x6 = 450 mg/m2;
## cisplatin 120 mg/m2 x4 = 480 mg/m2; HDMTX 12 g/m2 x12.  Surgery week 11.
MAP_DOX_W <- c(1, 6, 12, 17, 22, 27)
MAP_CIS_W <- c(1, 6, 12, 22)
MAP_MTX_W <- c(4, 5, 9, 10, 15, 16, 20, 21, 25, 26, 30, 31)
SURGERY_W <- 11

## High-dose methotrexate 12 g/m2 over 4 h, with leucovorin 15 mg/m2 q6h
## from +24 h for 48 h.  The rescue delay is the therapeutic window.
ev_mtx <- function(week, dose_g_m2 = 12, lv_hours = 48, lv_mg_m2 = 15,
                   gluc_at = NA) {
  t0   <- (week - 1) * W
  mmol <- dose_g_m2 * BSA * 1000 / 454.44
  out  <- ev(amt = mmol, cmt = "A1", time = t0, tinf = 4)
  for (i in seq_len(lv_hours %/% 6)) {
    out <- c(out, ev(amt = lv_mg_m2 * BSA, cmt = "LV",
                     time = t0 + 24 + 6 * (i - 1)))
  }
  if (!is.na(gluc_at)) {
    out <- c(out, ev(amt = 50 * BSA, cmt = "GLUC", time = t0 + gluc_at))
  }
  out
}

## Doxorubicin 75 mg/m2 as a 48-h infusion; dexrazoxane 10:1 concurrent.
ev_dox <- function(week, dose = 75, dexrazoxane = FALSE) {
  t0  <- (week - 1) * W
  out <- c(ev(amt = dose * BSA, cmt = "DOXC", time = t0, tinf = 48),
           ev(amt = dose,       cmt = "CUMDOX", time = t0))
  if (dexrazoxane) {
    out <- c(out, ev(amt = 10 * dose * BSA, cmt = "DEXR", time = t0,
                     tinf = 48))
  }
  out
}

## Cisplatin 120 mg/m2 as a 72-h infusion.
ev_cis <- function(week, dose = 120) {
  ev(amt = dose * BSA, cmt = "CISC", time = (week - 1) * W, tinf = 72)
}

## Ifosfamide 2.8 g/m2/day x5 with etoposide 100 mg/m2/day x5 (MAPIE).
ev_ie <- function(week) {
  t0  <- (week - 1) * W
  out <- NULL
  for (d in 0:4) {
    out <- c(out,
             ev(amt = 2800 * BSA * 0.55, cmt = "IFOA", time = t0 + 24 * d,
                tinf = 3),
             ev(amt = 100 * BSA, cmt = "ETOC", time = t0 + 24 * d, tinf = 1))
  }
  out
}

ev_deno <- function(week) ev(amt = 120, cmt = "DENO", time = (week - 1) * W)
ev_zol  <- function(week) ev(amt = 4,   cmt = "ZOL",  time = (week - 1) * W)
ev_mifa <- function(week) {
  t0 <- (week - 1) * W
  c(ev(amt = 2 * BSA, cmt = "MIFA", time = t0),
    ev(amt = 2 * BSA, cmt = "MIFA", time = t0 + 84))
}

## ===========================================================================
##  THE PROTOCOL GO / NO-GO GATES
##
##  This layer is not decoration.  Real protocols do not administer a
##  scheduled course into an unrecovered patient, and that is the mechanism
##  by which ESCALATION SPENDS THE BUDGET INSTEAD OF ADDING TO IT: every
##  extra ifosfamide/etoposide cycle deepens the neutrophil nadir and shaves
##  the kidney, and the courses that get deferred or dropped are courses of
##  the backbone that was doing the killing.  EURAMOS-1 delivered full
##  protocol therapy to 76% of MAP patients and 51% of MAPIE patients; this
##  layer reproduces the ~25-point gap from the gate rules rather than
##  assuming it.
## ===========================================================================
GATES <- list(
  ## kind = c(min ANC 10^9/L, max mucositis grade, min kidney fraction)
  MTX = c(anc = 0.75, muc = 3.0, kid = 0.60),
  DOX = c(anc = 1.00, muc = 3.0, kid = 0.00),
  CIS = c(anc = 1.00, muc = 3.0, kid = 0.60),
  IE  = c(anc = 1.00, muc = 3.0, kid = 0.60)
)
MAX_DEFER_WEEKS   <- 2
DEFER_DOSE_FACTOR <- 0.75    # a deferred course returns at -25% dose
KID_STOP_MTX      <- 0.55    # grade-3 AKI holds methotrexate
KID_RESUME_MTX    <- 0.70    # ... and it resumes on renal recovery
CARDIAC_STOP_LVEF <- 50.0    # no further anthracycline below this

lvef_of <- function(cmv) 62 * (0.35 + 0.65 * max(cmv, 0)^0.6)

## Walk the schedule one week at a time, deciding each course against the
## current state, and integrate week by week.  Returns the simulation plus
## the delivery bookkeeping that Section 3 of the reference output reads.
simulate_osa <- function(reg = list(), tmax_w = 40) {
  p <- reg$p
  m <- if (length(p)) param(mod, p) else mod

  sched <- list()
  push  <- function(w, spec) {
    k <- as.character(w)
    sched[[k]] <<- c(sched[[k]], list(spec))
  }
  for (w in reg$dox_weeks %||% MAP_DOX_W)
    push(w, list(kind = "DOX", dose = reg$dox_dose %||% 75, defer = 0))
  for (w in reg$cis_weeks %||% MAP_CIS_W)
    push(w, list(kind = "CIS", dose = reg$cis_dose %||% 120, defer = 0))
  for (w in reg$mtx_weeks %||% MAP_MTX_W)
    push(w, list(kind = "MTX", dose = reg$mtx_dose %||% 12, defer = 0))
  for (w in reg$ifo_weeks %||% integer(0))
    push(w, list(kind = "IE", dose = 2800, defer = 0))
  for (w in reg$deno_weeks %||% integer(0))
    push(w, list(kind = "DENO", dose = 120, defer = 0))
  for (w in reg$zol_weeks %||% integer(0))
    push(w, list(kind = "ZOL", dose = 4, defer = 0))
  for (w in reg$mifa_weeks %||% integer(0))
    push(w, list(kind = "MIFA", dose = 2, defer = 0))

  init      <- as.numeric(init(m))
  names(init) <- names(init(m))
  y         <- init
  pending   <- list()
  given     <- c(MTX = 0, DOX = 0, CIS = 0, IE = 0)
  on_time   <- given
  omitted   <- given
  planned   <- given
  for (k in names(sched))
    for (s in sched[[k]])
      if (s$kind %in% names(planned)) planned[s$kind] <- planned[s$kind] + 1
  out       <- NULL
  huvos     <- NA_real_
  mtx_hold  <- FALSE
  surgery_w <- reg$surgery_week %||% SURGERY_W

  for (wk in seq_len(tmax_w)) {
    t0 <- (wk - 1) * W

    ## ---- definitive surgery: the primary leaves the model -------------
    if (wk == surgery_w && (reg$surgery %||% TRUE)) {
      tot   <- y["PRIM"] + y["NEC"]
      huvos <- if (tot > 0) 100 * y["NEC"] / tot else 0
      y["PRIM"] <- y["PRIM"] * 0.002        # macroscopically complete
      y["NEC"]  <- 0
    }

    ## ---- gate this week's courses -------------------------------------
    if (y["KID"] < KID_STOP_MTX)      mtx_hold <- TRUE
    else if (y["KID"] > KID_RESUME_MTX) mtx_hold <- FALSE

    due     <- c(sched[[as.character(wk)]], pending)
    pending <- list()
    evs     <- NULL
    for (c0 in due) {
      if (mtx_hold && c0$kind == "MTX") {
        omitted["MTX"] <- omitted["MTX"] + 1
        next
      }
      g  <- GATES[[c0$kind]]
      ok <- TRUE
      if (!is.null(g)) {
        if (y["CIRC"] < g["anc"])                          ok <- FALSE
        if (y["MUC"]  >= g["muc"])                         ok <- FALSE
        if (g["kid"] > 0 && y["KID"] < g["kid"])           ok <- FALSE
        if (c0$kind == "DOX" &&
            lvef_of(y["CMV"]) < CARDIAC_STOP_LVEF)         ok <- FALSE
      }
      if (ok) {
        evs <- c(evs, switch(c0$kind,
          MTX  = ev_mtx(wk, c0$dose, reg$lv_hours %||% 48,
                        reg$lv_dose %||% 15, reg$gluc_at %||% NA),
          DOX  = ev_dox(wk, c0$dose, reg$dexrazoxane %||% FALSE),
          CIS  = ev_cis(wk, c0$dose),
          IE   = ev_ie(wk),
          DENO = ev_deno(wk),
          ZOL  = ev_zol(wk),
          MIFA = ev_mifa(wk)))
        if (c0$kind %in% names(given)) {
          given[c0$kind] <- given[c0$kind] + 1
          if (c0$defer == 0) on_time[c0$kind] <- on_time[c0$kind] + 1
        }
      } else {
        c0$defer <- c0$defer + 1
        c0$dose  <- c0$dose * DEFER_DOSE_FACTOR
        if (c0$defer <= MAX_DEFER_WEEKS) pending <- c(pending, list(c0))
        else if (c0$kind %in% names(omitted))
          omitted[c0$kind] <- omitted[c0$kind] + 1
      }
    }

    sim <- m %>%
      init(y) %>%
      data_set(if (is.null(evs)) ev(amt = 0, cmt = "A1", time = t0) else
                 as_data_set(evs)) %>%
      mrgsim(start = t0, end = t0 + W, delta = 1, hmax = 4,
             atol = 1e-10, rtol = 1e-7, recover = "*") %>%
      as.data.frame()
    y   <- unlist(sim[nrow(sim), names(init)])
    out <- rbind(out, sim)
  }
  for (c0 in pending)
    if (c0$kind %in% names(omitted)) omitted[c0$kind] <- omitted[c0$kind] + 1

  list(sim = out, huvos = huvos, given = given, planned = planned,
       on_time = on_time, omitted = omitted,
       on_schedule = sum(on_time) / max(sum(planned), 1))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ===========================================================================
##  THE POPULATION CURE LAYER
##
##  A single deterministic run gives the cure probability of a patient sitting
##  exactly at the mean parameters, and because P(cure) is a DOUBLE
##  exponential in K that number is violently sensitive near the operating
##  point.  What a trial measures is the cure FRACTION of a population whose
##  K is spread out.  The scenario supplies the MEAN of K; this layer supplies
##  the SPREAD (SLC19A1/ABCC2 genotype, methotrexate clearance CV ~30%,
##  MTHFR, delivered dose intensity after the gates, intrinsic
##  chemosensitivity); only the two together are comparable to a published
##  EFS curve.  None of those sources of spread has ever been randomised.
## ===========================================================================
POP_N <- 8000
set.seed(20260805)
POP_Z   <- rnorm(POP_N)
POP_LAM <- rpois(POP_N, 1.80)

pop_cure <- function(K_nats, cv = 0.25, n0 = 1e6, lam = POP_LAM) {
  K    <- pmax(POP_Z * cv * max(K_nats, 1e-9) + K_nats, 0)
  surv <- n0 * exp(-K)
  ples <- -expm1(-surv)                     # P(one lesion survives)
  mean((1 - ples)^lam)                      # cure needs ALL of them to fail
}

## Under proportional hazards the long-term event-free proportion satisfies
## S_test = S_ref^HR, so HR = ln(cure_test)/ln(cure_ref).  HR < 1 favours test.
efs_hr <- function(cure_test, cure_ref) {
  log(min(max(cure_test, 1e-9), 1 - 1e-9)) /
    log(min(max(cure_ref, 1e-9), 1 - 1e-9))
}

## ===========================================================================
##  THERAPY SCENARIOS (20)
## ===========================================================================
SCENARIOS <- list(
  S01_surgery_only = list(
    label = "Surgery alone (historical control, no chemotherapy)",
    dox_weeks = integer(0), cis_weeks = integer(0), mtx_weeks = integer(0)),

  S02_MAP = list(
    label = "MAP - EURAMOS-1 arm A (reference standard)"),

  S03_MAPIE = list(
    label = "MAP+IE (MAPIE) - poor-responder escalation, EURAMOS-1",
    ## As randomised: doxorubicin 375 mg/m2 (5 cycles), cisplatin 360 mg/m2
    ## (3 cycles), 12 HDMTX courses, 5 ifosfamide/etoposide cycles.  The
    ## escalation BUYS the IE cycles by giving up one doxorubicin and one
    ## cisplatin cycle from a budget already at its walls.
    ifo_weeks = c(14, 18, 23, 28, 32),
    dox_weeks = c(1, 6, 12, 21, 26),
    cis_weeks = c(1, 6, 12)),

  S04_MAP_pH60 = list(
    label = "MAP with urine pH 6.0 (inadequate alkalinisation)",
    p = list(URINE_PH = 6.0)),

  S05_MAP_pH75 = list(
    label = "MAP with urine pH 7.5 + full hydration",
    p = list(URINE_PH = 7.5)),

  S06_MAP_lowhydration = list(
    label = "MAP at half hydration (1.5 L/m2/day), pH 7.5",
    p = list(UF0 = 0.105)),

  S07_MAP_pH66 = list(
    label = "MAP with urine pH 6.6 (control for the rescue arm)",
    p = list(URINE_PH = 6.6)),

  S07b_MAP_glucarpidase = list(
    label = "MAP, pH 6.6, glucarpidase rescue at 48 h",
    p = list(URINE_PH = 6.6), gluc_at = 48),

  S08_MAP_dexrazoxane = list(
    label = "MAP + dexrazoxane 10:1 with every doxorubicin dose",
    dexrazoxane = TRUE),

  S09_MAP_dox600 = list(
    label = "MAP intensified to doxorubicin 600 mg/m2 (unprotected)",
    dox_dose = 100),

  S10_MAP_dox600_dexra = list(
    label = "MAP intensified to 600 mg/m2 WITH dexrazoxane",
    dox_dose = 100, dexrazoxane = TRUE),

  S11_MAP_zoledronate = list(
    label = "MAP + zoledronate 4 mg q4w (OS2006 design)",
    zol_weeks = seq(1, 29, by = 4)),

  S12_MAP_denosumab = list(
    label = "MAP + denosumab 120 mg q4w",
    deno_weeks = seq(1, 29, by = 4)),

  S13_MAP_mifamurtide = list(
    label = "MAP + mifamurtide (INT-0133 schedule)",
    mifa_weeks = 12:33),

  S14_AP_only = list(
    label = "AP only - adult/elderly, methotrexate omitted",
    mtx_weeks = integer(0)),

  S15_MAP_metastatic = list(
    label = "MAP in overt metastatic disease at diagnosis",
    lambda0 = 15.0),

  S16_MAP_resistant = list(
    label = "MAP in an intrinsically resistant tumour (RES0 = 0.25)",
    p = list(RES0 = 0.25)),

  S17_MAP_regorafenib = list(
    label = "MAP + regorafenib maintenance (anti-angiogenic, REGOBONE)",
    p = list(VEGFI = 1.2)),

  S18_MAP_dosedense = list(
    label = "MAP with 25% higher methotrexate (15 g/m2)",
    mtx_dose = 15),

  S19_MAP_pre_op_IE = list(
    label = "MAP + the same IE cycles moved PRE-operatively",
    ifo_weeks = c(2, 3, 7, 8, 12),
    dox_weeks = c(1, 6, 17, 22, 27),
    cis_weeks = c(1, 6, 22))
)

## ===========================================================================
##  DRIVER
## ===========================================================================
run_all <- function() {
  rows <- lapply(names(SCENARIOS), function(nm) {
    reg <- SCENARIOS[[nm]]
    r   <- simulate_osa(reg)
    s   <- r$sim
    fin <- s[nrow(s), ]
    K   <- fin$KNATS
    lam <- if (!is.null(reg$lambda0)) rpois(POP_N, reg$lambda0) else POP_LAM
    cure <- pop_cure(K, lam = lam)
    data.frame(
      scenario   = nm,
      label      = reg$label,
      huvos_pct  = r$huvos,
      log10_kill = fin$KLOG10,
      cure       = cure,
      TRM        = fin$PTRM,
      survival   = cure * (1 - fin$PTRM),
      on_sched   = 100 * r$on_schedule,
      n_MTX      = r$given["MTX"], n_DOX = r$given["DOX"],
      n_CIS      = r$given["CIS"], n_IE  = r$given["IE"],
      eGFR_nadir = 100 * min(s$eGFRo) / 120,
      LVEF       = fin$LVEFo,
      hearing_dB = fin$HLo,
      ANC_nadir  = min(s$CIRC),
      MUC_max    = max(s$MUCo),
      cum_dox    = fin$CUMDOX,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  ref <- out$survival[out$scenario == "S02_MAP"]
  out$EFS_HR <- vapply(out$survival, efs_hr, 1.0, cure_ref = ref)
  out
}

## ---------------------------------------------------------------------------
##  ANALYSIS 1 — the bifurcation.  Bisect on urine pH for the three
##  thresholds: supersaturation, delayed elimination, and AKI.
## ---------------------------------------------------------------------------
mtx_single_course <- function(pH = 7.5, dose_g_m2 = 12, uf_mult = 1,
                              kid0 = 1, gluc_at = NA, hours = 336) {
  m <- param(mod, list(URINE_PH = pH, UF0 = 0.21 * uf_mult))
  y <- as.numeric(init(m)); names(y) <- names(init(m))
  y["KID"] <- kid0
  s <- m %>% init(y) %>%
    data_set(as_data_set(ev_mtx(1, dose_g_m2, gluc_at = gluc_at))) %>%
    mrgsim(end = hours, delta = 0.25, hmax = 1, atol = 1e-11, rtol = 1e-8,
           recover = "*") %>% as.data.frame()
  a <- function(h) approx(s$time, s$CMTX_uM, h)$y
  list(C24 = a(24), C48 = a(48), C72 = a(72), Cmax = max(s$CMTX_uM),
       SS_max = max(s$SSRAT), kid_final = tail(s$KID, 1),
       ppt_max = max(s$PPT),
       auc = sum(diff(s$time) * head(s$CMTX_uM / 1000, -1)),
       delayed = (a(48) > 1) || (a(72) > 0.1),
       aki = tail(s$KID, 1) < 0.75)
}

critical_ph <- function(pred, lo = 5.5, hi = 9.0, iters = 34) {
  if (!pred(mtx_single_course(pH = lo))) return(NA_real_)
  for (i in seq_len(iters)) {
    mid <- (lo + hi) / 2
    if (pred(mtx_single_course(pH = mid))) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

## Expected (from osa_reference_output.txt, same equations in Python):
##   supersaturation onset  pH 7.30
##   delayed elimination    pH 7.15
##   AKI (>25% GFR lost)    pH 6.95
## and the exchange rate log10(f)/SOL_B: one pH unit = 10^0.682 = 4.81x the
## fluid, or a 4.81-fold dose reduction.  Alkalinisation, not the dose, is
## arithmetically the largest lever on this regimen's safety.

## ---------------------------------------------------------------------------
##  ANALYSIS 2 — MAPIE by responder status.  Every kill term carries
##  (1 - RES); a poor responder is a patient in whom that factor is small,
##  so the escalation was tested in exactly the subgroup where extra
##  cytotoxic is worth least.  Expected: HR 0.98 in poor responders,
##  0.86 in good responders.  EURAMOS-1 reported 0.98 (0.78-1.23).
## ---------------------------------------------------------------------------
mapie_by_responder <- function() {
  do.call(rbind, lapply(c(0.02, 0.20), function(r0) {
    pair <- lapply(c("S02_MAP", "S03_MAPIE"), function(nm) {
      reg <- SCENARIOS[[nm]]
      reg$p <- c(reg$p, list(RES0 = r0))
      r <- simulate_osa(reg)
      fin <- r$sim[nrow(r$sim), ]
      list(K = fin$KNATS, cure = pop_cure(fin$KNATS),
           trm = fin$PTRM, sched = 100 * r$on_schedule)
    })
    a <- pair[[1]]; b <- pair[[2]]
    data.frame(RES0 = r0,
               responder = if (r0 < 0.1) "good" else "poor (trial population)",
               MAP_K = a$K / log(10), MAPIE_K = b$K / log(10),
               MAP_cure = a$cure, MAPIE_cure = b$cure,
               MAPIE_sched = b$sched,
               HR = efs_hr(b$cure * (1 - b$trm), a$cure * (1 - a$trm)))
  }))
}

## ---------------------------------------------------------------------------
##  ANALYSIS 3 — the detectability wall.  At the real spread of achieved
##  log-kill, a realistic intensification is worth a couple of points of
##  survival, and the trial that could see it does not exist.
## ---------------------------------------------------------------------------
detectability <- function(Kbar, cv = 0.25) {
  base <- pop_cure(Kbar, cv)
  do.call(rbind, lapply(c(0.10, 0.20, 0.30, 0.50), function(dl) {
    c1   <- pop_cure(Kbar + dl * log(10), cv)
    gain <- c1 - base
    pbar <- (c1 + base) / 2
    data.frame(delta_log10 = dl, cure = c1, gain = gain,
               HR = efs_hr(c1, base),
               n_per_arm = 2 * pbar * (1 - pbar) * (1.96 + 0.842)^2 / gain^2)
  }))
}
## Expected: +0.20 log10 -> +0.037 absolute, HR 0.88, n = 2764 per arm.
## EURAMOS-1 randomised 618 poor responders to MAPIE.  The trial was not
## underpowered by accident; it was underpowered by the shape of the
## dose-response surface it was built on top of.

## ---------------------------------------------------------------------------
##  ANALYSIS 4 — not the mean, and not the tail either: the SHOULDER.
##
##  P(cure) saturates at BOTH ends, so an extra log of kill is worth nothing
##  in patients who were already going to be cured and nothing in patients no
##  achievable dose can rescue.  Spend the same total extra log-kill on
##  different parts of the achieved-exposure distribution and the answers are
##  not close.  Expected (osa_reference_output.txt section 4b-ii):
##
##    everyone, +0.25 log each            cure 0.649   gain +0.048
##    bottom quartile, +1.00 log each     cure 0.632   gain +0.030  <- WORSE
##    SECOND quartile, +1.00 log each     cure 0.717   gain +0.116  <- 2.4x
##    top quartile, +1.00 log each        cure 0.604   gain +0.002
##
##  The marginal value peaks at 6.2-6.6 log10 of achieved kill and is
##  essentially zero in the bottom two and top two deciles.  Since necrosis
##  correlates with achieved log-kill at r = 0.97, the marker CAN identify
##  that band — the near-misses at roughly 85-89% necrosis.  EURAMOS-1 used
##  the marker correctly as a prognostic tool and then spent its escalation
##  budget on the group the marker had identified as LEAST salvageable.
## ---------------------------------------------------------------------------
budget_targeting <- function(Kbar, cv = 0.25, spend_log10 = 0.25) {
  Kp    <- pmax(POP_Z * cv * Kbar + Kbar, 0)
  cure  <- function(K) mean((1 - -expm1(-1e6 * exp(-K)))^POP_LAM)
  total <- spend_log10 * log(10) * POP_N
  q     <- quantile(Kp, c(0, .25, .50, .75, 1))
  masks <- list(
    "everyone (uniform intensification)" = rep(TRUE, POP_N),
    "bottom quartile (worst responders)" = Kp <= q[2],
    "SECOND quartile (the near-misses)"  = Kp > q[2] & Kp <= q[3],
    "third quartile"                     = Kp > q[3] & Kp <= q[4],
    "top quartile (already cured)"       = Kp > q[4],
    "middle half (25-75th percentile)"   = Kp > q[2] & Kp <= q[4])
  base <- cure(Kp)
  do.call(rbind, lapply(names(masks), function(nm) {
    mk  <- masks[[nm]]
    per <- total / sum(mk)
    K2  <- Kp; K2[mk] <- K2[mk] + per
    c1  <- cure(K2)
    data.frame(who = nm, each_log10 = per / log(10), cure = c1,
               gain = c1 - base, HR = efs_hr(c1, base),
               stringsAsFactors = FALSE)
  }))
}

## Where the marginal value lives, decile by decile.
marginal_by_decile <- function(Kbar, cv = 0.25) {
  Kp   <- pmax(POP_Z * cv * Kbar + Kbar, 0)
  cure <- function(K) mean((1 - -expm1(-1e6 * exp(-K)))^POP_LAM)
  base <- cure(Kp)
  eps  <- 0.05 * log(10)
  do.call(rbind, lapply(0:9, function(d) {
    b  <- quantile(Kp, c(d, d + 1) / 10)
    mk <- Kp >= b[1] & Kp <= b[2]
    K2 <- Kp; K2[mk] <- K2[mk] + eps
    data.frame(decile = d + 1, lo_log10 = b[1] / log(10),
               hi_log10 = b[2] / log(10),
               marginal_value = (cure(K2) - base) /
                 (eps / log(10) * sum(mk) / POP_N))
  }))
}

## ===========================================================================
##  RUN
## ===========================================================================
if (identical(environment(), globalenv()) && !interactive()) {
  res <- run_all()
  print(res[, c("scenario", "huvos_pct", "log10_kill", "cure", "TRM",
                "survival", "EFS_HR", "on_sched", "eGFR_nadir", "LVEF",
                "hearing_dB", "ANC_nadir")], digits = 3)
  cat("\n--- critical urine pH -------------------------------------------\n")
  cat("supersaturation :", critical_ph(function(r) r$SS_max > 1), "\n")
  cat("delayed elim.   :", critical_ph(function(r) r$delayed), "\n")
  cat("AKI             :", critical_ph(function(r) r$aki), "\n")
  cat("1 pH unit is worth", round(10^0.682, 2), "x the fluid\n")
  cat("\n--- MAPIE by responder status -----------------------------------\n")
  print(mapie_by_responder(), digits = 3)
  Kbar <- res$log10_kill[res$scenario == "S02_MAP"] * log(10)
  cat("\n--- the detectability wall --------------------------------------\n")
  print(detectability(Kbar), digits = 3)
  cat("\n--- where the budget belongs (the shoulder) ---------------------\n")
  print(budget_targeting(Kbar), digits = 3)
  print(marginal_by_decile(Kbar), digits = 3)
}
