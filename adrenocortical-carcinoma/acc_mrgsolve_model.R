## ===========================================================================
##  acc_mrgsolve_model.R
##  Adrenocortical Carcinoma (ACC) — QSP model
##  부신피질암 정량적 시스템 약리학 모델
##
##  ORGANISING THESIS
##  -----------------
##  Plasma mitotane is not a dose. It is the OVERFLOW of a slowly filling
##  lipid depot. For a constant intake the trough obeys
##
##      Cp(t) = (D*F/CL) * [1 - exp(-t*CL/Vss)]
##               \______/   \____________________/
##                WHERE            WHEN
##
##  with Vss ~ 6,000 L (CV 81.5%, Arshad 2018) — overwhelmingly ADIPOSE.
##  Efficacy is a threshold on Cp (>= 14 mg/L). Neurotoxicity is a SECOND
##  threshold on the SAME Cp (> 20 mg/L). And the same exposure induces
##  CYP3A4 four-fold, which destroys (a) etoposide, (b) hydrocortisone and
##  (c) mitotane's own clearance.
##
##  ONE state variable (mitotane in fat), THREE opposed reads. Six results
##  follow arithmetically rather than by assertion — see README.md and the
##  scenario block at the bottom of this file:
##
##    A. Fat sets WHEN, not WHERE. Doubling adipose mass leaves
##       Css = D*F/CL untouched and roughly doubles tau = Vss/CL.
##    B. The dose-regimen RCT randomised the SMALLER variable: one SD of
##       body composition moves time-to-target more than doubling the dose.
##    C. Hypercortisolism grows fat, so the disease lengthens its own
##       treatment delay (a feedback on the CLOCK, not the set-point).
##    D. EDP-M is weakest where mitotane is strongest (shared state var).
##    E. Total cortisol lies: CBG doubles, so free falls while total reads
##       reassuringly normal — and hydrocortisone CL doubles at the same
##       time. Two independent errors, same direction.
##    F. Cortisol is an immunosuppressant the tumour makes itself, so
##       PD-1 blockade fails in the secreting subset.
##
##  MODEL SIZE
##  ----------
##    52 ODE compartments · 148 parameters · 18 therapy scenarios ·
##    virtual-population layer · closed-form cross-check
##
##  UNITS
##  -----
##    time            days
##    mitotane        mg (amount) / mg/L (concentration)
##    cortisol        ug/dL (clinical units); volumes in dL
##    tumour          mL (1 mL ~ 1 g)
##    etoposide etc.  mg / mg/L
##    ACTH, DHEAS, ALDO, IGF2, immune states, organ states: relative to
##                    healthy baseline = 1 unless noted
##
##  CALIBRATION TARGETS (see acc_references.md for the numbered sources)
##  --------------------------------------------------------------------
##    Mitotane apparent Vss ~ 6,086 L, CV 81.5%; 1-cmt popPK with linear
##      enzyme autoinduction; first TDM advised day 16          [Arshad 2018]
##    Therapeutic window 14-20 mg/L; >20 mg/L neurotoxicity     [Terzolo 2013]
##    Terminal t1/2 18-160 d (median ~50-60 d)                  [Corso 2021]
##    High-dose vs low-dose start: popPK favours high dose, RCT
##      found no significant difference (n small, IIV huge)      [Kerkhofs 2013]
##    FIRM-ACT EDP-M: ORR 23.2%, PFS 5.0 mo; Sz-M ORR 9.2%,
##      PFS 2.1 mo; OS 14.8 vs 12.0 mo (NS)                     [Fassnacht 2012]
##    Mitotane induces CYP3A4 strongly and durably              [Kroiss 2011]
##    Hydrocortisone requirement roughly doubles                 [Chortis 2013]
##    CBG and SHBG roughly double (oestrogenic action)          [Nader 2006]
##    SOAT1 inhibition -> free cholesterol/oxysterol -> ER
##      stress -> apoptosis (the adrenolytic mechanism)          [Sbiera 2015]
##    Linsitinib phase 3 negative despite IGF2 in ~90%          [Fassnacht 2015]
##    Pembrolizumab ORR ~14-23%                                 [Habra 2019, Raj 2020]
##    Glucocorticoid excess depletes tumour-infiltrating
##      lymphocytes and worsens prognosis                        [Landwehr 2020]
##
##  IMPORTANT: this is an educational / research QSP model. Several
##  couplings (notably the etoposide induction factor, which is inferred
##  from the CYP3A4 induction magnitude and etoposide's CYP3A4-metabolised
##  fraction rather than measured directly in a mitotane-etoposide
##  interaction study) are stated assumptions, flagged inline as ASSUMED.
##  Do not use for clinical decisions.
## ===========================================================================

library(mrgsolve)
suppressMessages(library(dplyr))

acc_code <- '
$PARAM @annotated
// ---------------- patient / body composition -----------------------------
WT0      :  80   : Body weight (kg)
FATKG0   :  22   : Baseline adipose mass (kg) - reference BMI ~27
BSA      : 1.85  : Body surface area (m2)
ALBU     : 4.0   : Serum albumin (g/dL)
SECRETOR :  1    : Cortisol-secreting phenotype (1 yes / 0 no)

// ---------------- mitotane PK --------------------------------------------
MITF     : 0.35  : Mitotane bioavailability - fasting
FMEALMAX : 0.55  : Bioavailability with high-fat meal
MEALFAT  :  1    : Taken with fatty meal (1) or fasting (0)
KAMIT    : 4.0   : Mitotane absorption rate (1/day)
V1MIT    : 400   : Mitotane central volume - blood + rapidly-equilibrating lean tissue (L)
CLMIT0   :  48   : Mitotane baseline (uninduced) clearance (L/day)
QADI     : 900   : Mitotane central<->depot intercompartmental clearance (L/day)
VLEANTIS : 1100  : Non-adipose tissue distribution volume (L)
KPFAT    : 203   : Adipose partition volume per kg fat (L/kg)
FMIND    : 0.35  : Fraction of mitotane CL that is autoinducible
KE0ADR   : 0.50  : Adrenal/tumour effect-site equilibration (1/day)
KE0CNS   : 1.00  : CNS effect-site equilibration (1/day)

// ---------------- CYP3A4 induction ---------------------------------------
KENZ     : 0.15  : CYP3A4 activity turnover (1/day) - t1/2 ~4.6 d
EMAXIND  : 4.5   : Maximum fold-induction above baseline
EC50IND  : 8.0   : Mitotane conc for half-maximal induction (mg/L)
INDON    :  1    : Global induction switch (1 on / 0 off)
INDVIC   :  1    : Victim-drug induction switch - gates the etoposide/cortisol clearance effect only, leaving mitotane autoinduction intact

// ---------------- binding proteins ---------------------------------------
BMAXCBG0 :  22   : CBG binding capacity at baseline (ug/dL)
KDCBG    : 0.90  : CBG dissociation constant (ug/dL)
FALB     : 1.60  : Albumin linear binding coefficient
KCBG     : 0.05  : CBG turnover (1/day)
EMAXCBG  : 1.30  : Max fractional rise in CBG (-> ~2.3x)
EC50CBG  : 10.0  : Mitotane conc for half-maximal CBG rise (mg/L)
KSHBG    : 0.05  : SHBG turnover (1/day)
EMAXSHBG : 1.20  : Max fractional rise in SHBG

// ---------------- steroidogenesis / HPA ----------------------------------
VCORT    : 340   : Cortisol distribution volume (dL)
SYNNORM  : 10000 : Normal adrenal cortisol production (ug/day)
SYNTUM   : 60000 : Max autonomous tumour cortisol production (ug/day)
CLCORT0  : 44440 : Cortisol clearance on FREE conc basis (dL/day)
FMHC     : 0.45  : Fraction of cortisol CL that is CYP3A4-inducible
IC50ENZ  : 10.0  : Mitotane conc for 50% steroidogenic enzyme block (mg/L)
HENZ     : 1.50  : Hill coefficient - enzyme block
KACTH    : 6.0   : ACTH turnover (1/day)
IC50FB   : 0.225 : Free cortisol for 50% ACTH suppression (ug/dL)
HFB      : 2.0   : Hill coefficient - ACTH feedback
KADRREG  : 0.002 : Normal adrenal regrowth rate (1/day) - near-permanent loss
KLYS     : 0.035 : Max adrenolytic rate (1/day)
EC50LYS  : 14.0  : Mitotane effect-site conc for half-max adrenolysis (mg/L)
HLYS     : 3.0   : Hill coefficient - adrenolysis (threshold-like)
KSTCAP   : 0.02  : Tumour steroidogenic capacity turnover (1/day)
KGR      : 0.50  : GR occupancy equilibration (1/day)
EC50GR   : 0.60  : Free cortisol for half-max GR occupancy (ug/dL)
KDHEAS   : 0.08  : DHEAS turnover (1/day)
KOHP     : 0.50  : 17-OH-progesterone turnover (1/day)
KALDO    : 0.20  : Aldosterone turnover (1/day)
FGLOM    : 0.35  : Fraction of glomerulosa spared from adrenolysis

// ---------------- hydrocortisone replacement -----------------------------
FHC      : 0.95  : Hydrocortisone oral bioavailability
KAHC     :  30   : Hydrocortisone absorption rate (1/day)

// ---------------- body composition / GC end-organ ------------------------
KFAT     : 0.012 : Adipose mass turnover (1/day) - t1/2 ~58 d
AFAT     : 0.35  : Max fractional fat gain from GC excess
KMUSC    : 0.02  : Muscle mass turnover (1/day)
AMUSC    : 0.35  : Max fractional muscle loss from GC excess
KBMD     : 0.002 : BMD turnover (1/day)
ABMD     : 0.25  : Max fractional BMD loss from GC excess

// ---------------- tumour --------------------------------------------------
KGROW    : 0.0047: Gompertz growth rate constant (1/day)
TVMAX    : 20000 : Carrying capacity (mL)
TUM0     : 120   : Initial sensitive tumour volume (mL)
TUMR0    : 2.0   : Initial resistant clone volume (mL)
KRES     : 2e-4  : Selection flux into the resistant clone (1/day)
FRESK    : 0.12  : Residual chemo sensitivity of the resistant clone
SLPMIT   : 0.010 : Mitotane cytotoxic slope (1/day per unit Hill effect)
EC50MITK : 16.0  : Mitotane effect-site conc for half-max tumour kill (mg/L)
HMITK    : 3.0   : Hill coefficient - mitotane tumour kill
SLPETO   : 0.35  : Etoposide kill slope (1/day per mg/L)
SLPDOX   : 0.55  : Doxorubicin kill slope (1/day per mg/L)
SLPCIS   : 0.60  : Cisplatin kill slope (1/day per mg/L)
SLPSZ    : 0.030 : Streptozotocin kill slope (1/day per mg/L)
SLPIMM   : 0.010 : Immune kill slope (1/day per unit effector activity)
MDR1     : 1.6   : Constitutive ABCB1 efflux factor in ACC
FPGPIND  : 0.25  : Further P-gp induction by mitotane

// ---------------- etoposide ----------------------------------------------
CLETO0   :  67   : Etoposide clearance (L/day at BSA 1.85)
V1ETO    :  10   : Etoposide central volume (L)
QETO     :  60   : Etoposide intercompartmental clearance (L/day)
V2ETO    :  12   : Etoposide peripheral volume (L)
FM3A4ETO : 0.35  : ASSUMED CYP3A4-metabolised fraction of etoposide

// ---------------- doxorubicin --------------------------------------------
CLDOX    : 1300  : Doxorubicin clearance (L/day)
V1DOX    :  25   : Doxorubicin central volume (L)
QDOX     : 900   : Doxorubicin intercompartmental clearance (L/day)
V2DOX    : 1600  : Doxorubicin peripheral volume (L)

// ---------------- cisplatin ----------------------------------------------
CLCIS    : 500   : Free cisplatin clearance (L/day)
V1CIS    :  20   : Free cisplatin volume (L)
KBIND    :  35   : Irreversible protein binding rate (1/day)

// ---------------- streptozotocin ----------------------------------------
CLSZ     : 700   : Streptozotocin clearance (L/day)
V1SZ     :  25   : Streptozotocin volume (L)

// ---------------- pembrolizumab -----------------------------------------
CLPEM    : 0.22  : Pembrolizumab clearance (L/day)
V1PEM    : 7.4   : Pembrolizumab volume (L)
KDPEM    : 0.05  : Apparent PD-1 binding constant (mg/L)
EMAXPD1  : 2.2   : Max fold increase in effector drive from PD-1 blockade

// ---------------- linsitinib / IGF axis ---------------------------------
CLLIN    : 150   : Linsitinib clearance (L/day) - Css ~0.9 mg/L on 150 mg BID
V1LIN    : 400   : Linsitinib volume (L)
KALIN    :  12   : Linsitinib absorption rate (1/day)
FLIN     : 0.45  : Linsitinib bioavailability
IC50R1   : 0.25  : Linsitinib IC50 at IGF1R (mg/L)
IC50RA   : 1.50  : Linsitinib IC50 at tumour IR-A (mitogenic, IGF2-avid)
IC50RB   : 0.35  : Linsitinib IC50 at metabolic IR-B - more sensitive than tumour IR-A, so hyperglycaemia arrives before antitumour effect
KIGF2    : 0.20  : IGF2 turnover (1/day)
KSIG     : 0.80  : IGF/AKT/mTOR signalling turnover (1/day)
WIGF     : 0.50  : Weight of IGF2-IGF1R input to proliferation drive
WINS     : 0.50  : Weight of insulin-IRA input to proliferation drive
APROLIF  : 0.30  : Fractional change in growth rate per unit IGF signalling - ACC growth is NOT IGF-limited, which is half the reason GALACTIC failed

// ---------------- glucose / insulin -------------------------------------
VGLU     : 100   : Glucose distribution volume (dL)
GPROD    : 11000 : Hepatic glucose production (mg/day)
AGLUGR   : 0.55  : Max fractional rise in glucose production from GC excess
GUPT     : 64.6  : Glucose disposal coefficient (dL/day) - set so GLU 92 mg/dL is a fixed point at INS 1
SIINS    : 0.85  : Insulin sensitivity coefficient
VINS     : 50    : Insulin distribution volume (dL)
KINS     :  75   : Beta-cell secretory gain - set so INS = 1 at GLU = 92
GTHR     : 80    : Glucose threshold for insulin secretion (mg/dL)
CLINS    : 900   : Insulin clearance (dL/day)
KBCELL   : 0.030 : Streptozotocin beta-cell toxicity slope (1/day per mg/L)

// ---------------- immune -------------------------------------------------
KTEFF    : 0.10  : Effector T-cell turnover (1/day)
KTREG    : 0.08  : Treg / MDSC turnover (1/day)
GCIC50   : 0.55  : GR occupancy for 50% T-cell suppression
HGC      : 2.0   : Hill coefficient - GC immunosuppression
ATREG    : 0.80  : Max fractional Treg expansion from GC excess
TMBFAC   : 1.0   : Neoantigen / TMB multiplier on effector drive

// ---------------- neutrophils (Friberg) ---------------------------------
CIRC0    : 4.5   : Baseline neutrophil count (10^9/L)
MTTN     : 5.4   : Mean transit time (days)
GAMN     : 0.17  : Feedback exponent
SLPNETO  : 0.90  : Etoposide myelosuppression slope (per mg/L)
SLPNDOX  : 1.40  : Doxorubicin myelosuppression slope (per mg/L)

// ---------------- toxicity ------------------------------------------------
KNTOX    : 0.10  : CNS injury turnover (1/day)
CNSTHR   : 24.0  : CNS injury EC50 (mg/L) - above the clinical 20 mg/L alert level so 20 sits on the rising limb
HNTOX    : 8.0   : Hill coefficient - CNS threshold steepness
KALT     : 0.10  : ALT turnover (1/day)
ALT0     : 25    : Baseline ALT (U/L)
AALT     : 2.2   : Max fold rise in ALT
EC50ALT  : 14.0  : Mitotane conc for half-max hepatic effect (mg/L)
KFT4     : 0.06  : Free T4 turnover (1/day)
FT40     : 1.20  : Baseline free T4 (ng/dL)
AFT4     : 0.45  : Max fractional free T4 suppression
EC50FT4  : 12.0  : Mitotane conc for half-max FT4 suppression (mg/L)
LT4SUP   :  0    : Levothyroxine replacement (0-1 fractional restoration)
KLVREC   : 0.004 : LVEF recovery rate (1/day)
LVEF0    : 62    : Baseline LVEF (%)
KLVDOX   : 8.00  : Doxorubicin LVEF slope (per mg/L per day)
DOXTHR   : 400   : Cumulative anthracycline inflection (mg/m2)
KGFRREC  : 0.010 : eGFR recovery rate (1/day)
GFR0     : 92    : Baseline eGFR (mL/min/1.73m2)
AGFRCIS  : 0.30  : Max fractional eGFR loss from cumulative platinum
CIS50    : 900   : Cumulative bound platinum for half-max eGFR loss (mg)
HAZAI    : 0.030 : Adrenal crisis hazard scale (1/day)
CFADQ    : 0.16  : Free cortisol considered adequate (ug/dL)

$CMT @annotated
MITG   : Mitotane gut depot (mg)
MITC   : Mitotane central (mg)
MITA   : Mitotane ADIPOSE DEPOT (mg)
MITE   : Mitotane adrenal/tumour effect site (mg/L)
MITN   : Mitotane CNS effect site (mg/L)
ENZ    : CYP3A4 relative activity (baseline 1)
CBG    : CBG relative abundance (baseline 1)
SHBG   : SHBG relative abundance (baseline 1)
ETOC   : Etoposide central (mg)
ETOP   : Etoposide peripheral (mg)
DOXC   : Doxorubicin central (mg)
DOXP   : Doxorubicin peripheral (mg)
DOXCUM : Cumulative anthracycline (mg/m2)
CISC   : Free cisplatin (mg)
CISP   : Bound platinum (mg)
SZC    : Streptozotocin (mg)
HCG    : Hydrocortisone gut depot (ug)
PEMC   : Pembrolizumab central (mg)
LING   : Linsitinib gut depot (mg)
LINC   : Linsitinib central (mg)
TUMS   : Sensitive tumour volume (mL)
TUMR   : Resistant clone volume (mL)
STCAP  : Tumour steroidogenic capacity (relative)
ADRN   : Normal adrenal cortex mass (relative)
CORT   : Total plasma cortisol (ug/dL)
ACTH   : Plasma ACTH (relative)
GRO    : Glucocorticoid receptor occupancy (0-1)
DHEAS  : DHEA-S (relative)
OHP17  : 17-OH-progesterone (relative)
ALDO   : Aldosterone (relative)
FATKG  : Adipose mass (kg)
MUSC   : Muscle mass (relative)
BMD    : Bone mineral density (relative)
IGF2   : IGF2 (relative)
IGFSIG : IGF/AKT/mTOR signalling (relative)
GLU    : Plasma glucose (mg/dL)
INS    : Plasma insulin (relative)
TEFF   : Effector CD8 T-cell activity (relative)
TREG   : Treg / MDSC activity (relative)
PROLN  : Neutrophil progenitors (10^9/L)
TR1N   : Neutrophil transit 1 (10^9/L)
TR2N   : Neutrophil transit 2 (10^9/L)
TR3N   : Neutrophil transit 3 (10^9/L)
CIRCN  : Circulating neutrophils (10^9/L)
NTOX   : Cumulative CNS injury score (0-1)
ALT    : ALT (U/L)
FT4    : Free T4 (ng/dL)
LVEF   : LVEF (%)
GFR    : eGFR (mL/min/1.73m2)
AIHAZ  : Cumulative adrenal crisis hazard
AUCMIT : Cumulative mitotane AUC (mg/L*day)
TIW    : Cumulative time in therapeutic window (days)

$MAIN
// ---- initial conditions -------------------------------------------------
ENZ_0    = 1.0;
CBG_0    = 1.0;
SHBG_0   = 1.0;
TUMS_0   = TUM0;
TUMR_0   = TUMR0;
STCAP_0  = 1.0;
ADRN_0   = 1.0;
ACTH_0   = 1.0;
DHEAS_0  = 1.0;
OHP17_0  = 1.0;
ALDO_0   = 1.0;
FATKG_0  = FATKG0;
MUSC_0   = 1.0;
BMD_0    = 1.0;
IGF2_0   = 1.0;
IGFSIG_0 = 1.0;
GLU_0    = 92.0;
INS_0    = 1.0;
TEFF_0   = 1.0;
TREG_0   = 1.0;
PROLN_0  = CIRC0;
TR1N_0   = CIRC0;
TR2N_0   = CIRC0;
TR3N_0   = CIRC0;
CIRCN_0  = CIRC0;
ALT_0    = ALT0;
FT4_0    = FT40;
LVEF_0   = LVEF0;
GFR_0    = GFR0;

// Baseline total cortisol depends on whether the tumour secretes.
// Solved below from the steady state of the CORT equation; a plain
// non-secretor sits near 5 ug/dL.
CORT_0   = 5.0;
GRO_0    = IC50FB / (EC50GR + IC50FB);   // healthy axis must start AT its fixed point

// ---- mitotane bioavailability (meal-dependent) --------------------------
double FMIT = MEALFAT > 0.5 ? FMEALMAX : MITF;

// ---- the DEPOT VOLUME: this is where body composition enters ------------
// Vss = V1 + VLEANTIS + KPFAT * FATKG.  At FATKG 22 kg this is
// 400 + 1100 + 4466 = 5966 L, matching the reported apparent Vss ~6,086 L.
double VDEPOT = VLEANTIS + KPFAT * FATKG;

$ODE
// =======================================================================
//  MITOTANE PK — the depot that sets the clock
// =======================================================================
double CMIT   = MITC / V1MIT;                      // plasma mg/L
double CDEP   = MITA / VDEPOT;                     // depot pseudo-conc mg/L
double CLMITt = CLMIT0 * (1.0 + FMIND * (ENZ - 1.0));   // AUTOINDUCTION

dxdt_MITG = -KAMIT * MITG;
dxdt_MITC =  KAMIT * MITG * FMIT
             - CLMITt * CMIT
             - QADI * (CMIT - CDEP);
dxdt_MITA =  QADI * (CMIT - CDEP);

// effect sites
dxdt_MITE = KE0ADR * (CMIT - MITE);
dxdt_MITN = KE0CNS * (CMIT - MITN);

// cumulative exposure and time-in-window bookkeeping
dxdt_AUCMIT = CMIT;
dxdt_TIW    = (CMIT >= 14.0 && CMIT <= 20.0) ? 1.0 : 0.0;

// =======================================================================
//  READ 2 — CYP3A4 INDUCTION (harm #1): one state, four victims
// =======================================================================
double INDDRIVE = INDON * EMAXIND * CMIT / (EC50IND + CMIT);
// ENZV is the induction level SEEN BY VICTIM DRUGS. Gating it separately
// from ENZ lets scenario 08 remove the etoposide interaction WITHOUT also
// removing the mitotane autoinduction (which would change mitotane exposure
// and confound the comparison).
dxdt_ENZ = KENZ * (1.0 + INDDRIVE - ENZ);

// =======================================================================
//  BINDING PROTEINS and the free-cortisol solve (result E)
// =======================================================================
dxdt_CBG  = KCBG  * (1.0 + EMAXCBG  * CMIT / (EC50CBG + CMIT) - CBG);
dxdt_SHBG = KSHBG * (1.0 + EMAXSHBG * CMIT / (EC50CBG + CMIT) - SHBG);

// Free cortisol: CORT = CF*(1+FALB) + Bmax*CF/(Kd+CF)  ->  quadratic in CF
double BMAXC = BMAXCBG0 * CBG;
double aq    = 1.0 + FALB;
double bq    = aq * KDCBG + BMAXC - CORT;
double disc  = bq * bq + 4.0 * aq * CORT * KDCBG;
if (disc < 0.0) disc = 0.0;
double CFREE = (-bq + sqrt(disc)) / (2.0 * aq);
if (CFREE < 1e-8) CFREE = 1e-8;

// =======================================================================
//  READ 1 — MITOTANE PD: fast reversible enzyme block + slow adrenolysis
// =======================================================================
// (a) FAST, REVERSIBLE: steroidogenic enzyme inhibition
double ENZBLK = 1.0 / (1.0 + pow(MITE / IC50ENZ, HENZ));

// (b) SLOW, IRREVERSIBLE: adrenolysis (SOAT1 -> ER stress -> apoptosis)
double LYSD = pow(MITE, HLYS) / (pow(EC50LYS, HLYS) + pow(MITE, HLYS));

dxdt_ADRN  = KADRREG * (ACTH - ADRN) - KLYS * LYSD * ADRN;
dxdt_STCAP = KSTCAP * (1.0 - STCAP)  - KLYS * LYSD * STCAP;

// =======================================================================
//  CORTISOL: two inputs (adrenal + tumour + oral HC), one INDUCED clearance
// =======================================================================
double TUMTOT  = TUMS + TUMR;
// Saturating burden term, = 1 at the baseline tumour volume and -> 2 as the
// burden grows. Secretory output and IGF2 release per unit volume FALL as the
// tumour enlarges (dedifferentiation, necrosis), so neither may scale linearly.
double TFRAC   = 2.0 * TUMTOT / (TUMTOT + TUM0 + 1e-9);
double SYNADR  = SYNNORM * ADRN * ACTH * ENZBLK;
double SYNTUMf = SYNTUM * SECRETOR * STCAP * TFRAC * ENZBLK;
double ENZV    = 1.0 + INDVIC * (ENZ - 1.0);
double CLCORTt = CLCORT0 * (1.0 - FMHC + FMHC * ENZV);   // INDUCED clearance

dxdt_HCG  = -KAHC * HCG;
dxdt_CORT = (SYNADR + SYNTUMf + KAHC * HCG * FHC - CLCORTt * CFREE) / VCORT;

// ACTH negative feedback (fast)
// Normalised so that ACTH = 1 when free cortisol sits at IC50FB (its normal
// value): a healthy axis must be a FIXED POINT of this equation, not drift.
double ACTHT = 2.0 / (1.0 + pow(CFREE / IC50FB, HFB));
dxdt_ACTH = KACTH * (ACTHT - ACTH);

// GR occupancy -> the effect variable everything downstream reads
dxdt_GRO = KGR * (CFREE / (EC50GR + CFREE) - GRO);
double GRO0 = IC50FB / (EC50GR + IC50FB);       // occupancy at normal free
double GREX = (GRO - GRO0) / GRO0;              // fractional GC excess
if (GREX < -1.0) GREX = -1.0;

// other steroid readouts
dxdt_DHEAS = KDHEAS * (ADRN * ACTH * ENZBLK
                       + SECRETOR * STCAP * TFRAC * ENZBLK - DHEAS);
dxdt_OHP17 = KOHP   * (ADRN * ACTH + SECRETOR * STCAP * TFRAC
                       * (1.0 + 1.8 * (1.0 - ENZBLK)) - OHP17);   // substrate piles up
double GLOMSPARE = FGLOM + (1.0 - FGLOM) * ADRN;
dxdt_ALDO  = KALDO * (GLOMSPARE - ALDO);

// =======================================================================
//  RESULT C — the feedback on the CLOCK: cortisol builds the depot
// =======================================================================
double FATT = FATKG0 * (1.0 + AFAT * (GREX > 0 ? GREX/(1.0+GREX) : 0.0));
dxdt_FATKG = KFAT * (FATT - FATKG);

dxdt_MUSC = KMUSC * (1.0 - AMUSC * (GREX > 0 ? GREX/(1.0+GREX) : 0.0) - MUSC);
dxdt_BMD  = KBMD  * (1.0 - ABMD  * (GREX > 0 ? GREX/(1.0+GREX) : 0.0) - BMD);

// =======================================================================
//  CYTOTOXIC PK
// =======================================================================
double CETO = ETOC / V1ETO;
double CDOX = DOXC / V1DOX;
double CCIS = CISC / V1CIS;
double CSZ  = SZC  / V1SZ;
double CPEM = PEMC / V1PEM;
double CLIN = LINC / V1LIN;

// RESULT D — HARM #1 lands here: induction doubles etoposide clearance
double CLETOt = CLETO0 * (BSA/1.85) * (1.0 - FM3A4ETO + FM3A4ETO * ENZV);

dxdt_ETOC = -CLETOt * CETO - QETO * (CETO - ETOP/V2ETO);
dxdt_ETOP =  QETO * (CETO - ETOP/V2ETO);

dxdt_DOXC = -CLDOX * CDOX - QDOX * (CDOX - DOXP/V2DOX);
dxdt_DOXP =  QDOX * (CDOX - DOXP/V2DOX);
dxdt_DOXCUM = 0.0;                              // incremented by dosing events

double CLCISt = CLCIS * (GFR / GFR0);           // nephrotoxicity feeds back
dxdt_CISC = -CLCISt * CCIS - KBIND * CISC;
dxdt_CISP =  KBIND * CISC;

dxdt_SZC = -CLSZ * CSZ;

dxdt_PEMC = -CLPEM * CPEM;

dxdt_LING = -KALIN * LING;
dxdt_LINC =  KALIN * LING * FLIN - CLLIN * CLIN;

// =======================================================================
//  IGF AXIS — why a target present in 90% of tumours did not deliver
// =======================================================================
double R1FREE = 1.0 / (1.0 + CLIN / IC50R1);    // IGF1R left unblocked
double RAFREE = 1.0 / (1.0 + CLIN / IC50RA);    // tumour IR-A left unblocked (escape)
double RBFREE = 1.0 / (1.0 + CLIN / IC50RB);    // metabolic IR-B left unblocked

dxdt_IGF2 = KIGF2 * (TFRAC - IGF2);

double SIGT = WIGF * IGF2 * R1FREE + WINS * INS * RAFREE;
dxdt_IGFSIG = KSIG * (SIGT / (WIGF + WINS) - IGFSIG);

// glucose / insulin: linsitinib blocks peripheral IR -> glucose up ->
// insulin up -> insulin re-drives the partially-blocked IR-A on the tumour
double GPRODt = GPROD * (1.0 + AGLUGR * (GREX > 0 ? GREX/(1.0+GREX) : 0.0));
double BCELL  = 1.0 - KBCELL * CSZ; if (BCELL < 0.15) BCELL = 0.15;
dxdt_GLU = (GPRODt - GUPT * GLU * (1.0 + SIINS * INS * RBFREE)) / VGLU;
dxdt_INS = (KINS * BCELL * (GLU > GTHR ? GLU - GTHR : 0.0) - CLINS * INS) / VINS;

// =======================================================================
//  READ 3 / RESULT F — IMMUNE COMPARTMENT
//  cortisol is an immunosuppressant the tumour administers to itself
// =======================================================================
// normalised so that IMMSUP = 1 at normal free cortisol: what matters is
// the DEVIATION in glucocorticoid tone, not its absolute level
double IMMRAW = 1.0 / (1.0 + pow(GRO  / GCIC50, HGC));
double IMMREF = 1.0 / (1.0 + pow(GRO0 / GCIC50, HGC));
double IMMSUP = IMMRAW / IMMREF;
double OCCPD1 = CPEM / (KDPEM + CPEM);
double TEFFT  = TMBFAC * (1.0 + EMAXPD1 * OCCPD1) * IMMSUP;

dxdt_TEFF = KTEFF * (TEFFT - TEFF);
dxdt_TREG = KTREG * (1.0 + ATREG * (GREX > 0 ? GREX/(1.0+GREX) : 0.0) - TREG);

double TREGF = 1.0 + 0.6 * (TREG - 1.0);
if (TREGF < 0.2) TREGF = 0.2;
double TEFFNET = TEFF / TREGF;

// =======================================================================
//  TUMOUR — Gompertz growth with five kill terms
// =======================================================================
double PGP    = MDR1 * (1.0 + FPGPIND * (ENZ - 1.0));
double DELIV  = 1.0 / PGP;                       // intratumoural availability

double KMITK  = SLPMIT * pow(MITE, HMITK) / (pow(EC50MITK, HMITK) + pow(MITE, HMITK));
double KCHEM  = (SLPETO * CETO + SLPDOX * CDOX) * DELIV + SLPCIS * CCIS;
double KSZK   = SLPSZ * CSZ;
double KIMMK  = SLPIMM * TEFFNET;

double GRWTH  = KGROW * (1.0 + APROLIF * (IGFSIG - 1.0));
if (GRWTH < 0.0) GRWTH = 0.0;
double LOGB   = log(TVMAX / (TUMTOT > 1e-6 ? TUMTOT : 1e-6));
if (LOGB < 0.0) LOGB = 0.0;

dxdt_TUMS = GRWTH * TUMS * LOGB
            - (KMITK + KCHEM + KSZK + KIMMK) * TUMS
            - KRES * TUMS;
dxdt_TUMR = GRWTH * TUMR * LOGB
            - (FRESK * (KMITK + KCHEM + KSZK) + KIMMK) * TUMR
            + KRES * TUMS;

// =======================================================================
//  MYELOSUPPRESSION (Friberg 4-transit)
// =======================================================================
double KTRN  = 4.0 / MTTN;
double EDRUG = SLPNETO * CETO + SLPNDOX * CDOX;
if (EDRUG > 0.95) EDRUG = 0.95;
double FBK   = pow(CIRC0 / (CIRCN > 0.05 ? CIRCN : 0.05), GAMN);

dxdt_PROLN = KTRN * PROLN * (1.0 - EDRUG) * FBK - KTRN * PROLN;
dxdt_TR1N  = KTRN * (PROLN - TR1N);
dxdt_TR2N  = KTRN * (TR1N  - TR2N);
dxdt_TR3N  = KTRN * (TR2N  - TR3N);
dxdt_CIRCN = KTRN * TR3N - KTRN * CIRCN;

// =======================================================================
//  TOXICITY
// =======================================================================
double NTD = pow(MITN / CNSTHR, HNTOX) / (1.0 + pow(MITN / CNSTHR, HNTOX));
dxdt_NTOX = KNTOX * (NTD - NTOX);

dxdt_ALT = KALT * (ALT0 * (1.0 + AALT * CMIT / (EC50ALT + CMIT)) - ALT);

dxdt_FT4 = KFT4 * (FT40 * (1.0 - AFT4 * CMIT / (EC50FT4 + CMIT))
                   + LT4SUP * FT40 * AFT4 * CMIT / (EC50FT4 + CMIT) - FT4);

double DOXVULN = 1.0 + DOXCUM / DOXTHR;
dxdt_LVEF = KLVREC * (LVEF0 - LVEF) - KLVDOX * CDOX * DOXVULN;

// Nephrotoxicity tracks CUMULATIVE platinum (CISP, which never clears in
// this model and therefore IS the cumulative dose), not the 40-minute free
// cisplatin spike.
dxdt_GFR = KGFRREC * (GFR0 * (1.0 - AGFRCIS * CISP / (CISP + CIS50)) - GFR);

// adrenal crisis hazard accumulates whenever free cortisol is inadequate
dxdt_AIHAZ = HAZAI * exp(-CFREE / CFADQ);

$TABLE
double CMITo   = MITC / V1MIT;
double VDEPOTo = VLEANTIS + KPFAT * FATKG;
double VSSo    = V1MIT + VDEPOTo;
double TAUo    = VSSo / (CLMIT0 * (1.0 + FMIND * (ENZ - 1.0)));
double CSSo    = 0.0;   // filled in by the driver script per regimen

double BMAXCo  = BMAXCBG0 * CBG;
double aqo     = 1.0 + FALB;
double bqo     = aqo * KDCBG + BMAXCo - CORT;
double disco   = bqo*bqo + 4.0*aqo*CORT*KDCBG;
if (disco < 0.0) disco = 0.0;
double CFREEo  = (-bqo + sqrt(disco)) / (2.0*aqo);
double FFRAC   = (CORT > 1e-6) ? CFREEo / CORT : 0.0;

double TUMTOTo = TUMS + TUMR;
double SLDo    = 20.0 * pow(3.0*TUMTOTo/(4.0*3.14159265), 1.0/3.0);  // mm
double RESPo   = 100.0 * (TUMTOTo - (TUM0 + TUMR0)) / (TUM0 + TUMR0);

double GRO0o   = IC50FB / (EC50GR + IC50FB);
double CETOo   = ETOC / V1ETO;
double ENZVo   = 1.0 + INDVIC * (ENZ - 1.0);
double CLETOo  = CLETO0 * (BSA/1.85) * (1.0 - FM3A4ETO + FM3A4ETO * ENZVo);
double INDFOLD = CLETOo / (CLETO0 * (BSA/1.85));
double OCCPD1o = (PEMC/V1PEM) / (KDPEM + PEMC/V1PEM);
double INWIN   = (CMITo >= 14.0 && CMITo <= 20.0) ? 1.0 : 0.0;
double ABOVE20 = (CMITo > 20.0) ? 1.0 : 0.0;
double DEPFRAC = (MITC + MITA > 1e-9) ? MITA/(MITC+MITA) : 0.0;

$CAPTURE @annotated
CMITo   : Plasma mitotane (mg/L)
VSSo    : Apparent steady-state volume (L)
TAUo    : Mitotane time constant Vss/CL (days)
CFREEo  : FREE cortisol (ug/dL) - the biologically active one
FFRAC   : Free cortisol fraction
TUMTOTo : Total tumour volume (mL)
SLDo    : RECIST-like sum of diameters (mm)
RESPo   : Change from baseline tumour burden (%)
GRO0o   : Reference GR occupancy
CETOo   : Etoposide concentration (mg/L)
CLETOo  : Etoposide clearance (L/day)
INDFOLD : Fold increase in etoposide clearance
OCCPD1o : PD-1 receptor occupancy
INWIN   : Currently in 14-20 mg/L window (0/1)
ABOVE20 : Currently above 20 mg/L (0/1)
DEPFRAC : Fraction of body mitotane sitting in the depot
'

mod <- mcode("acc", acc_code, soloc = tempdir())

## ===========================================================================
##  DOSING HELPERS
## ===========================================================================

## Mitotane: divided daily oral dose, given as one daily amount for the
## purposes of the depot model (absorption t1/2 ~4 h makes intra-day
## structure irrelevant on a 50-110 day time constant).
ev_mit <- function(g_per_day, start = 0, dur_days = 365) {
  if (g_per_day <= 0) return(NULL)
  ev(amt = g_per_day * 1000, cmt = "MITG", ii = 1, addl = dur_days - 1,
     time = start)
}

## Mitotane with a loading/taper schedule: named vector of g/day by segment
ev_mit_sched <- function(sched) {
  out <- NULL
  t <- 0
  for (i in seq_len(nrow(sched))) {
    seg <- sched[i, ]
    if (seg$g > 0) {
      e <- ev(amt = seg$g * 1000, cmt = "MITG", ii = 1,
              addl = seg$days - 1, time = t)
      out <- if (is.null(out)) e else c(out, e)
    }
    t <- t + seg$days
  }
  out
}

## Hydrocortisone tid (ug) - the saturable CBG binding makes the PEAKINESS
## of this input matter, so it is dosed three times daily on purpose.
ev_hc <- function(mg_per_day, dur_days = 365, start = 0) {
  if (mg_per_day <= 0) return(NULL)
  per <- mg_per_day * 1000 / 3
  e1 <- ev(amt = per, cmt = "HCG", ii = 1, addl = dur_days - 1, time = start)
  e2 <- ev(amt = per, cmt = "HCG", ii = 1, addl = dur_days - 1, time = start + 0.25)
  e3 <- ev(amt = per, cmt = "HCG", ii = 1, addl = dur_days - 1, time = start + 0.5)
  c(e1, e2, e3)
}

## EDP q28d: etoposide 100 mg/m2 d2-4, doxorubicin 40 mg/m2 d1,
## cisplatin 40 mg/m2 d3-4   (FIRM-ACT arm A schedule)
ev_edp <- function(bsa = 1.85, cycles = 8, start = 0) {
  out <- NULL
  add <- function(e) out <<- if (is.null(out)) e else c(out, e)
  for (k in 0:(cycles - 1)) {
    t0 <- start + k * 28
    add(ev(amt = 40 * bsa,  cmt = "DOXC",   time = t0 + 0))
    add(ev(amt = 40,        cmt = "DOXCUM", time = t0 + 0))  # mg/m2, not mg
    for (d in 1:3) add(ev(amt = 100 * bsa, cmt = "ETOC", time = t0 + d))
    for (d in 2:3) add(ev(amt = 40 * bsa,  cmt = "CISC", time = t0 + d))
  }
  out
}

## Streptozotocin: 1 g d1-5 induction, then 2 g q21d (FIRM-ACT arm B)
ev_sz <- function(cycles = 8, start = 0) {
  out <- ev(amt = 1000, cmt = "SZC", ii = 1, addl = 4, time = start)
  for (k in 1:cycles) out <- c(out, ev(amt = 2000, cmt = "SZC",
                                       time = start + 21 * k))
  out
}

ev_pem <- function(cycles = 12, start = 0) {
  ev(amt = 200, cmt = "PEMC", ii = 21, addl = cycles - 1, time = start)
}

ev_lin <- function(dur_days = 180, start = 0) {
  ## linsitinib 150 mg BID
  e1 <- ev(amt = 150, cmt = "LING", ii = 1, addl = dur_days - 1, time = start)
  e2 <- ev(amt = 150, cmt = "LING", ii = 1, addl = dur_days - 1, time = start + 0.5)
  c(e1, e2)
}

comb <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(ev(amt = 0, cmt = "MITG", time = 0))
  Reduce(function(a, b) c(a, b), parts)
}

run <- function(events, end = 365, delta = 0.25, ...) {
  p <- list(...)
  m <- mod
  if (length(p)) m <- param(m, p)
  out <- mrgsim(m, events = events, end = end, delta = delta,
                atol = 1e-8, rtol = 1e-6, maxsteps = 500000) %>% as_tibble()
  ## dosing records duplicate their time stamp; keep the post-dose row
  out[!duplicated(out$time, fromLast = TRUE), ]
}

## Standard high-dose start then TDM-guided maintenance
sched_highdose <- data.frame(g = c(6, 4, 3), days = c(42, 42, 400))
sched_lowdose  <- data.frame(g = c(1, 2, 3), days = c(14, 14, 456))

## ===========================================================================
##  SCENARIOS (18)
## ===========================================================================
scenarios <- list(

  ## ---- A. the depot: same regimen, three body compositions -------------
  s01 = list(
    label = "01 High-dose start · reference fat (22 kg)",
    ev    = ev_mit_sched(sched_highdose),
    par   = list(FATKG0 = 22),
    end   = 480),

  s02 = list(
    label = "02 High-dose start · LEAN (10 kg fat)",
    ev    = ev_mit_sched(sched_highdose),
    par   = list(FATKG0 = 10),
    end   = 480),

  s03 = list(
    label = "03 High-dose start · OBESE (45 kg fat)",
    ev    = ev_mit_sched(sched_highdose),
    par   = list(FATKG0 = 45),
    end   = 480),

  s04 = list(
    label = "04 Low-dose ramp · reference fat",
    ev    = ev_mit_sched(sched_lowdose),
    par   = list(FATKG0 = 22),
    end   = 480),

  ## ---- B. secreting tumour, hormone control ----------------------------
  s05 = list(
    label = "05 Secreting ACC · mitotane monotherapy",
    ev    = comb(ev_mit_sched(sched_highdose), ev_hc(20, 480)),
    par   = list(SECRETOR = 1),
    end   = 480),

  s06 = list(
    label = "06 Secreting ACC · NO mitotane (natural history)",
    ev    = ev_hc(0),
    par   = list(SECRETOR = 1),
    end   = 480),

  ## ---- C. EDP-M and the induction counterfactual (result D) ------------
  s07 = list(
    label = "07 EDP-M (FIRM-ACT arm A) · induction ON",
    ev    = comb(ev_mit_sched(sched_highdose), ev_edp(cycles = 10),
                 ev_hc(50, 480)),
    par   = list(INDON = 1),
    end   = 400),

  s08 = list(
    label = "08 EDP-M · victim-drug induction OFF (counterfactual)",
    ev    = comb(ev_mit_sched(sched_highdose), ev_edp(cycles = 10),
                 ev_hc(50, 480)),
    par   = list(INDVIC = 0),
    end   = 400),

  s09 = list(
    label = "09 EDP alone (no mitotane)",
    ev    = ev_edp(cycles = 10),
    par   = list(),
    end   = 400),

  s10 = list(
    label = "10 Sz-M (FIRM-ACT arm B)",
    ev    = comb(ev_mit_sched(sched_highdose), ev_sz(cycles = 10),
                 ev_hc(50, 480)),
    par   = list(),
    end   = 400),

  ## ---- D. the replacement trap (result E) -----------------------------
  s11 = list(
    label = "11 Mitotane + hydrocortisone 20 mg/day (standard dose)",
    ev    = comb(ev_mit_sched(sched_highdose), ev_hc(20, 480)),
    par   = list(SECRETOR = 0),
    end   = 400),

  s12 = list(
    label = "12 Mitotane + hydrocortisone 50 mg/day (doubled)",
    ev    = comb(ev_mit_sched(sched_highdose), ev_hc(50, 480)),
    par   = list(SECRETOR = 0),
    end   = 400),

  s13 = list(
    label = "13 Hydrocortisone 20 mg/day, NO mitotane (control)",
    ev    = ev_hc(20, 480),
    par   = list(SECRETOR = 0),
    end   = 400),

  ## ---- E. immunotherapy and the steroid it meets (result F) -----------
  s14 = list(
    label = "14 Pembrolizumab · cortisol-SECRETING ACC",
    ev    = ev_pem(cycles = 16),
    par   = list(SECRETOR = 1),
    end   = 360),

  s15 = list(
    label = "15 Pembrolizumab · NON-secreting ACC",
    ev    = ev_pem(cycles = 16),
    par   = list(SECRETOR = 0),
    end   = 360),

  s16 = list(
    label = "16 Pembrolizumab + mitotane (steroid controlled)",
    ev    = comb(ev_pem(cycles = 16), ev_mit_sched(sched_highdose),
                 ev_hc(50, 400)),
    par   = list(SECRETOR = 1),
    end   = 360),

  ## ---- F. the IGF target that did not deliver -------------------------
  s17 = list(
    label = "17 Linsitinib monotherapy (GALACTIC)",
    ev    = ev_lin(dur_days = 300),
    par   = list(),
    end   = 300),

  ## ---- G. withdrawal: the long tail -----------------------------------
  s18 = list(
    label = "18 Mitotane stopped at day 240 (washout)",
    ev    = ev_mit_sched(data.frame(g = c(6, 4, 3), days = c(42, 42, 156))),
    par   = list(),
    end   = 600)
)

run_scenario <- function(s) {
  args <- c(list(events = s$ev, end = s$end %||% 365), s$par)
  do.call(run, args)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ===========================================================================
##  VIRTUAL POPULATION — RESULT B
##  Does randomising the STARTING DOSE move time-to-target more than
##  the body-composition variability the trial did not stratify on?
## ===========================================================================
vpop_time_to_target <- function(n = 400, seed = 20260805) {
  set.seed(seed)
  ## Arshad 2018: IIV on V = 81.5% CV. Put that variability where it
  ## physically lives - adipose mass and partition - plus CL variability.
  fat <- pmax(4, rlnorm(n, log(22), 0.55))       # kg adipose
  clf <- rlnorm(n, 0, 0.30)                      # clearance multiplier
  out <- lapply(c("high", "low"), function(reg) {
    sched <- if (reg == "high") sched_highdose else sched_lowdose
    ## Return BOTH the censored endpoint (time to target, NA if never
     ## reached) and an always-defined continuous one (Cp at day 90). The
     ## censored endpoint cannot be used to compare the size of the dose
     ## effect against the fat effect, because conditioning on "reached
     ## target" preferentially DROPS the high-fat patients - exactly the
     ## ones the fat effect acts on. Cp90 has no such selection.
    res <- lapply(seq_len(n), function(i) {
      s <- run(ev_mit_sched(sched), end = 400, delta = 1,
               FATKG0 = fat[i], CLMIT0 = 48 * clf[i])
      hit <- s$CMITo >= 14
      data.frame(regimen = reg, fat = fat[i], clmult = clf[i],
                 ttt  = if (any(hit)) s$time[which(hit)[1]] else NA_real_,
                 cp90 = s$CMITo[s$time == 90],
                 hit90 = as.integer(any(s$CMITo[s$time <= 90] >= 14)))
    })
    bind_rows(res)
  })
  bind_rows(out)
}

## ===========================================================================
##  ANALYTIC CHECK — the two-factor decomposition, in closed form
##  Cp(t) = Css * (1 - exp(-t/tau));  Css = D*F/CL;  tau = Vss/CL
## ===========================================================================
analytic_depot <- function(g_per_day, fat_kg, F_ = 0.55, CL = 98.8,
                           V1 = 400, VLEAN = 1100, KP = 203) {
  Vss <- V1 + VLEAN + KP * fat_kg
  Css <- g_per_day * 1000 * F_ / CL
  tau <- Vss / CL
  ttt <- if (Css <= 14) NA_real_ else -tau * log(1 - 14 / Css)
  list(Vss = Vss, Css = Css, tau_days = tau, t_half_days = log(2) * tau,
       time_to_14 = ttt)
}

## ===========================================================================
##  If sourced interactively, print a short self-test.
## ===========================================================================
if (identical(environment(), globalenv()) && !interactive()) {
  message("acc_mrgsolve_model.R loaded: ", length(scenarios), " scenarios, ",
          length(mrgsolve::init(mod)), " compartments.")
}
