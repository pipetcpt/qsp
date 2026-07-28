## =============================================================================
##  Familial Chylomicronemia Syndrome (FCS) — QSP model for mrgsolve
##  fcs_mrgsolve_model.R
## =============================================================================
##
##  THE STRUCTURAL CLAIM THIS MODEL IS BUILT TO TEST
##  ------------------------------------------------
##  Triglyceride-rich lipoprotein (TRL) clearance is a SUM of two limbs:
##
##      dTG/dt = INPUT(diet, VLDL)
##               - CL_LPL * f_geno * C            <-- limb 1, first-order, huge
##               - Vmax_ind * g(apoC3) * C/(Km+C) <-- limb 2, saturable, small
##               - CL_res * C                     <-- residual scavenging
##
##  In health limb 1 carries >95% of the flux.  In FCS the genotype factor
##  f_geno is ~0 (biallelic LOF: LPL / APOC2 / APOA5 / GPIHBP1 / LMF1), so the
##  entire dietary fat flux is forced through limb 2, whose Vmax (1300 mg/h =
##  31 g TG/day) is the same order of magnitude as the daily fat load itself.
##  Four numerical consequences follow, and every analysis function below
##  exists to print one of them as a number:
##
##   (1) SATURATION CLIFF.  Steady-state TG is a hyperbola in dietary fat, not
##       a line.  FCS_saturation_curve() gives ~600 mg/dL at 10 g/day and
##       ~9500 mg/dL at 60 g/day.  The "<20 g fat/day" prescription is an
##       EDGE, not a slope, and that is why adherence behaves all-or-nothing.
##
##   (2) LIMB-1 DRUGS ARE ARITHMETICALLY DEAD.  Fibrates raise LPL
##       transcription, evinacumab de-represses LPL, insulin induces LPL.
##       All of them multiply f_geno.  FCS_limb_decomposition() prints the
##       flux carried by each limb in each arm: limb 1 stays at 0.0 mg/h in
##       every LPL-null arm no matter what is given.
##
##   (3) apoC-III KNOCKDOWN WORKS WITHOUT LIPASE, because apoC-III blocks
##       LDLR / LRP1 / HSPG remnant uptake — a route that never needed LPL.
##       In the model it multiplies Vmax of limb 2 by (1 + 2.2 x fractional
##       knockdown), and that single term reproduces APPROACH (-77%), Balance
##       (-74% at 12 months) and PALISADE (-80%) in a patient whose
##       post-heparin lipase activity is undetectable.
##
##   (4) THE PANCREATITIS ENDPOINT IS CONVEX IN TG (Hill exponent 1.7), so
##       E[lambda(TG)] > lambda(E[TG]).  FCS_jensen_gap() quantifies how much
##       of the true hazard is invisible to a FASTING TG measurement, and
##       FCS_trial_ledger() shows why a moderate mean TG reduction buys a much
##       larger reduction in acute pancreatitis events.
##
##  UNITS
##    time            hours (simulations run 0 - 8760 h = 1 year)
##    lipoprotein TG  mg (plasma amount); concentrations mg/dL (Vp = 30 dL)
##    dietary fat     g (dosed into FAT_GUT, 3 meals/day)
##    drug amounts    mg
##    hazard          per hour internally, reported per patient-year
##
##  43 ODE compartments · 5 drug PK models · 10 scenarios · 8 analysis functions
##
##  DISCLAIMER: research and education only.  Not for clinical use.
## =============================================================================

library(mrgsolve)

code <- '
$PROB
# Familial Chylomicronemia Syndrome (FCS) QSP model
# Two-limb TRL clearance | apoC-III axis | convex acute-pancreatitis hazard

$PARAM @annotated
// ---------------------------------------------------------------- physiology
VP        :   30      : Plasma volume (dL)
F_ABS     :    0.95   : Fractional absorption of dietary long-chain TG (-)
KA_FAT    :    0.50   : Gut lumen to enterocyte transfer (1/h)
K_ENT     :    0.70   : Enterocyte assembly and secretion (1/h)
K_LYM     :    0.50   : Lymph to plasma chylomicron appearance (1/h)
FMCT      :    0.00   : Fraction of dietary fat given as MCT (portal bypass) (-)

// ------------------------------------------------------------------- limb 1
CL_LPL_MAX:  140.0    : Maximal LPL-mediated TRL clearance in health (dL/h)
F_GENO    :    0.00   : Residual LPL genotype factor (0 = LPL-null FCS) (-)
A_C2      :    1.00   : apoC-II sufficiency factor (0 = APOC2 deficiency) (-)
A_A5      :    1.00   : apoA-V sufficiency factor (-)
IC50_C3L  :   20.0    : apoC-III giving 50 percent LPL inhibition (mg/dL)
IC50_ANG  :    1.00   : Normalised free ANGPTL3 giving 50 percent inhibition (-)
FRAC_REM  :    0.15   : Fraction of hydrolysed TRL-TG retained in remnant (-)

// ------------------------------------------------------------------- limb 2
VMAX_IND  : 1300      : Vmax of LPL-independent saturable TRL removal (mg/h)
KM_IND    : 1000      : Km of LPL-independent removal (mg/dL)
SEL_VL    :    3.00   : Selectivity of limb 2 for VLDL over chylomicrons (-)
IMAX_C3   :    2.20   : Max fold-rise of limb 2 on full apoC-III knockdown (-)
CL_RES    :    0.15   : Residual first-order TRL scavenging (dL/h)
K_REM0    :    0.60   : Baseline remnant hepatic uptake rate (1/h)

// ------------------------------------------------------------- hepatic VLDL
VSEC0     :  500      : Basal hepatic VLDL-TG secretion (mg/h)
A_ALC     :    0.45   : Alcohol multiplier on VLDL secretion (-)
A_EST     :    0.60   : Oestrogen multiplier on VLDL secretion (-)
A_DNL     :    0.20   : Hepatic TG-pool feedback on VLDL secretion (-)
LIVREF    : 37000     : Reference hepatic TG pool (mg)
AC3_VSEC  :    0.30   : apoC-III dependence of VLDL secretion (-)
I_OM3     :    0.25   : Omega-3 fractional suppression of VLDL secretion (-)
ALC       :    0      : Alcohol flag (0/1)
EST       :    0      : Oral oestrogen flag (0/1)
PREG      :    0      : Pregnancy flag (0/1)
OM3_ON    :    0      : Omega-3 4 g/day flag (0/1)
K_FAOX    :    0.0100 : Hepatic FA oxidation rate constant (1/h)
K_FFA_LIV :    1.50   : Plasma FFA to hepatic TG (1/h)
K_FFA_OUT :   12.00   : Plasma FFA disposal to peripheral tissue (1/h)
K_DNL     :  425      : De novo lipogenesis input to hepatic TG (mg/h)
ADIP_MOB  : 3600      : Lumped adipose FFA mobilisation (mg/h)

// ----------------------------------------------------------------- apoC-III
APOC3_REF :   25.0    : Untreated baseline plasma apoC-III (mg/dL)
KOUT_M    :    0.020  : APOC3 mRNA turnover (1/h)
KDEG_P    :    0.040  : Plasma apoC-III elimination (1/h)
EMAX_ASO  :    6.0    : Max fold-increase in mRNA degradation by ASO (-)
EC50_VOL  :  200      : Hepatic volanesorsen for half-max effect (mg)
EC50_OLE  :   88      : Hepatic olezarsen for half-max effect (mg)
EMAX_SI   :    6.0    : Max fold-increase in mRNA degradation by siRNA (-)
EC50_SI   :   18      : RISC-loaded plozasiran for half-max effect (AU)
IFIB_C3   :    0.20   : Fibrate suppression of APOC3 transcription (-)
A_CKD_C3  :    0.00   : Nephrotic or CKD apoC-III up-drive (-)

// ------------------------------------------------------------------ ANGPTL3
KOUT_ANG  :    0.030  : ANGPTL3 turnover (1/h)
IMAX_EVI  :    0.95   : Max evinacumab neutralisation of ANGPTL3 (-)
IC50_EVI  :    4.0    : Evinacumab plasma concentration for half-max (mg/L)

// -------------------------------------------------------------- LPL protein
KOUT_LPL  :    0.050  : Functional endothelial LPL pool turnover (1/h)
EFIB_LPL  :    0.50   : Max fibrate induction of LPL synthesis (-)
EC50_FIB  :    8.0    : Fenofibric acid for half-max PPAR-alpha effect (mg/L)

// ------------------------------------------------------------ pancreas / AP
K_PFFA    :    0.50   : Generation of local unbound FFA from chylomicrons (1/h)
K_PFFAOUT :    0.50   : Local unbound FFA disposal (1/h)
PFFA_THR  :  880      : Local unbound FFA injury threshold (AU)
K_INJ     :    0.0020 : Acinar injury formation (1/h per AU/1000)
K_REP     :    0.0100 : Acinar repair (1/h)
K_IL6     :    0.5000 : IL-6 generation from acinar injury (1/h)
K_IL6OUT  :    0.3000 : IL-6 elimination (1/h)
LAM_MAX   :    3.50   : Maximal AP hazard (events per patient-year)
TG50_AP   : 10000     : TG giving half-maximal AP hazard (mg/dL)
HILL_AP   :    1.70   : Hill convexity exponent of the AP hazard (-)
PRIME_AP  :    0.50   : Amplification of hazard by standing acinar injury (-)
TG_THR    :  880      : Clinical TG threshold, 10 mmol/L (mg/dL)

// ---------------------------------------------------------- other endpoints
KX_IN     :    0.0025 : Eruptive xanthoma formation (1/h per 1000 mg/dL)
KX_OUT    :    0.00048: Xanthoma resolution (1/h, half-life about 60 d)
KH_IN     :    0.00040: Hepatosplenomegaly accrual (1/h per 1000 mg/dL)
KH_OUT    :    0.00025: Hepatosplenomegaly regression (1/h)
KF_IN     :    0.0100 : Cognitive fog accrual from viscosity (1/h)
KF_OUT    :    0.0080 : Cognitive fog resolution (1/h)
K_VISC    :    0.55   : Viscosity index increment per 1000 mg/dL TG (-)

// ----------------------------------------------------------- platelets, ALT
KOUT_PLT  :    0.00415: Platelet elimination (1/h, half-life about 7 d)
PLT0      :  250      : Baseline platelet count (10^9/L)
IPLT_S    :    0.45   : Max suppression of thrombopoiesis by systemic PS-ASO (-)
KDES_PLT  :    1.00   : Max acceleration of platelet destruction (-)
IC50_PLT  :  400      : Systemic PS-ASO tissue amount for half-max effect (mg)
ALT0      :   25      : Baseline ALT (U/L)
KOUT_ALT  :    0.0100 : ALT turnover (1/h)
EALT_ASO  :    0.60   : Max ASO-driven ALT rise (-)
EC50_ALT  :  500      : Hepatic ASO for half-max ALT rise (mg)

// -------------------------------------------------- volanesorsen PK (PS-ASO)
F_VOL     :    0.70   : SC bioavailability (-)
KA_VOL    :    0.080  : SC absorption (1/h)
KL_VOL    :    0.090  : Plasma to hepatocyte distribution (1/h)
KO_VOL    :    0.100  : Plasma to other tissue including macrophages (1/h)
KE_VOL    :    0.010  : Plasma elimination and renal loss (1/h)
KOUTL_VOL :    0.00138: Hepatic tissue elimination (1/h, half-life about 21 d)
KOUTS_VOL :    0.00116: Systemic tissue elimination (1/h, half-life about 25 d)

// ----------------------------------------------- olezarsen PK (GalNAc3-ASO)
F_OLE     :    0.70   : SC bioavailability (-)
KA_OLE    :    0.090  : SC absorption (1/h)
KL_OLE    :    0.180  : ASGPR-mediated hepatocyte uptake (1/h)
KO_OLE    :    0.020  : Plasma to other tissue (1/h)
KE_OLE    :    0.010  : Plasma elimination (1/h)
KOUTL_OLE :    0.00138: Hepatic tissue elimination (1/h)
KOUTS_OLE :    0.00116: Systemic tissue elimination (1/h)

// --------------------------------------------- plozasiran PK (GalNAc-siRNA)
F_PLO     :    0.80   : SC bioavailability (-)
KA_PLO    :    0.100  : SC absorption (1/h)
KL_PLO    :    0.180  : ASGPR-mediated hepatocyte uptake (1/h)
KO_PLO    :    0.020  : Plasma to other tissue (1/h)
KE_PLO    :    0.010  : Plasma elimination (1/h)
KOUTL_PLO :    0.00138: Hepatic siRNA elimination (1/h)
KRISC_IN  :    0.00400: Endosomal escape and RISC loading (1/h)
KRISC_OUT :    0.000642 : RISC guide-strand loss (1/h, half-life about 45 d)

// ------------------------------------------------- evinacumab PK (2-cmt mAb)
V1_EVI    :    3.5    : Central volume (L)
V2_EVI    :    2.5    : Peripheral volume (L)
CL_EVI    :    0.0081 : Linear clearance (L/h, about 0.19 L/day)
Q_EVI     :    0.0208 : Intercompartmental clearance (L/h)

// ------------------------------------------------------ fenofibrate PK (1-cmt)
F_FIB     :    0.60   : Oral bioavailability (-)
KA_FIB    :    0.30   : Absorption (1/h)
V_FIB     :   18      : Volume of distribution of fenofibric acid (L)
KE_FIB    :    0.0347 : Elimination (1/h, half-life about 20 h)

$CMT @annotated
FAT_GUT    : Dietary long-chain TG in gut lumen (g)
ENT        : Enterocyte TG and nascent chylomicron pool (g)
LYMPH      : Lymphatic chylomicron TG in transit (g)
CM_TG      : Plasma chylomicron TG (mg)
CMREM_TG   : Plasma chylomicron remnant TG (mg)
VLDL_TG    : Plasma VLDL TG (mg)
VLDLREM_TG : Plasma VLDL remnant and IDL TG (mg)
LIVER_TG   : Hepatic cytosolic TG pool (mg)
ADIP_TG    : Net adipose TG balance since t0 (mg)
FFA        : Plasma free fatty acid pool (mg)
APOC3_M    : Hepatic APOC3 mRNA (normalised)
APOC3_P    : Plasma apoC-III (mg/dL)
ANGPTL3    : Plasma ANGPTL3 (normalised)
LPL_ACT    : Functional endothelial LPL pool (normalised)
PANC_FFA   : Pancreatic capillary unbound FFA (AU)
ACINAR     : Acinar injury burden (AU)
PANC_IL6   : Pancreatic and systemic IL-6 (AU)
CUMHAZ     : Cumulative acute pancreatitis hazard (expected events)
TAT880     : Cumulative time with TG above 880 mg/dL (h)
TGAUC      : Cumulative TG exposure (mg/dL*h)
XANTH      : Eruptive xanthoma burden (0-10)
HSM        : Hepatosplenomegaly index (0-10)
FOG        : Cognitive fog and fatigue score (0-10)
PLT        : Platelet count (10^9/L)
ALT        : Serum ALT (U/L)
VOL_SC     : Volanesorsen SC depot (mg)
VOL_CP     : Volanesorsen plasma (mg)
VOL_LIV    : Volanesorsen hepatic tissue (mg)
VOL_SYS    : Volanesorsen non-hepatic tissue (mg)
OLE_SC     : Olezarsen SC depot (mg)
OLE_CP     : Olezarsen plasma (mg)
OLE_LIV    : Olezarsen hepatic tissue (mg)
OLE_SYS    : Olezarsen non-hepatic tissue (mg)
PLO_SC     : Plozasiran SC depot (mg)
PLO_CP     : Plozasiran plasma (mg)
PLO_LIV    : Plozasiran hepatic tissue (mg)
PLO_RISC   : RISC-loaded plozasiran guide strand (AU)
EVI_C      : Evinacumab central (mg)
EVI_P      : Evinacumab peripheral (mg)
FIB_GUT    : Fenofibrate in gut (mg)
FIB_C      : Fenofibric acid in plasma (mg)
CUMFAT     : Cumulative absorbed dietary fat (g)
CUMIL6     : Cumulative IL-6 exposure (AU*h)

$GLOBAL
#define C_CM_  (CM_TG/VP)
#define C_CMR_ (CMREM_TG/VP)
#define C_VL_  (VLDL_TG/VP)
#define C_VLR_ (VLDLREM_TG/VP)
#define TGTOT_ (C_CM_ + C_CMR_ + C_VL_ + C_VLR_)

$MAIN
// ---- initial conditions ----------------------------------------------------
// The TRL pools start empty and are filled during the burn-in period that
// every scenario runs before treatment; this avoids hard-coding a baseline
// that would be wrong for a different genotype or diet.
APOC3_M_0  = 1.0;
APOC3_P_0  = APOC3_REF;
ANGPTL3_0  = 1.0;
LPL_ACT_0  = F_GENO;
LIVER_TG_0 = LIVREF;
FFA_0      = 300;
PLT_0      = PLT0;
ALT_0      = ALT0;

$ODE
// ============================ drug exposures ================================
double CEVI   = EVI_C / V1_EVI;                 // mg/L
double CFIB   = FIB_C / V_FIB;                  // mg/L
double ASOSYS = VOL_SYS + OLE_SYS;              // systemic PS-ASO burden (mg)
double ASOLIV = VOL_LIV + OLE_LIV;              // hepatic ASO burden (mg)

double EFF_VOL = EMAX_ASO * VOL_LIV /(EC50_VOL + VOL_LIV);
double EFF_OLE = EMAX_ASO * OLE_LIV /(EC50_OLE + OLE_LIV);
double EFF_SI  = EMAX_SI  * PLO_RISC/(EC50_SI  + PLO_RISC);
double FIBEFF  = CFIB/(EC50_FIB + CFIB);

// ============================ apoC-III axis =================================
double KIN_M  = KOUT_M;
double KSYN_P = KDEG_P * APOC3_REF;
dxdt_APOC3_M  = KIN_M * (1.0 - IFIB_C3*FIBEFF)
                - KOUT_M * APOC3_M * (1.0 + EFF_VOL + EFF_OLE + EFF_SI);
dxdt_APOC3_P  = KSYN_P * APOC3_M * (1.0 + A_CKD_C3) - KDEG_P * APOC3_P;

double C3REL = APOC3_P / APOC3_REF;                     // 1.0 when untreated
double KD3   = 1.0 - C3REL;  if(KD3 < 0.0) KD3 = 0.0;   // fractional knockdown
double GC3   = 1.0 + IMAX_C3 * KD3;                     // limb-2 amplifier

// ============================ ANGPTL3 and LPL ===============================
double KIN_ANG = KOUT_ANG;
dxdt_ANGPTL3   = KIN_ANG - KOUT_ANG * ANGPTL3;
double ANGFREE = ANGPTL3 * (1.0 - IMAX_EVI * CEVI/(IC50_EVI + CEVI));

double KIN_LPL = KOUT_LPL;
double PREGF   = 1.0 - 0.35*PREG;                       // pregnancy lowers LPL
dxdt_LPL_ACT   = KIN_LPL * F_GENO * (1.0 + EFIB_LPL*FIBEFF) * PREGF
                 - KOUT_LPL * LPL_ACT;

double CL_LPL = CL_LPL_MAX * LPL_ACT * A_C2 * A_A5
                / (1.0 + APOC3_P/IC50_C3L)
                / (1.0 + ANGFREE/IC50_ANG);             // dL/h
double KLPL   = CL_LPL / VP;                            // 1/h

// ============================ dietary fat input =============================
dxdt_FAT_GUT  = -KA_FAT * FAT_GUT;
dxdt_ENT      =  KA_FAT * FAT_GUT * F_ABS * (1.0 - FMCT) - K_ENT * ENT;
dxdt_LYMPH    =  K_ENT * ENT - K_LYM * LYMPH;
double IN_CM  =  K_LYM * LYMPH * 1000.0;                     // mg/h
double IN_MCT =  KA_FAT * FAT_GUT * F_ABS * FMCT * 1000.0;   // portal FFA mg/h
dxdt_CUMFAT   =  KA_FAT * FAT_GUT * F_ABS;

// ============================ hepatic VLDL ==================================
double VSEC = VSEC0
              * (1.0 + A_ALC*ALC + A_EST*(EST + 1.5*PREG))
              * (1.0 + A_DNL*(LIVER_TG/LIVREF - 1.0))
              * (1.0 + AC3_VSEC*(C3REL - 1.0))
              * (1.0 - I_OM3*OM3_ON);
if(VSEC < 0.0) VSEC = 0.0;

// ============================ two-limb clearance ============================
double TGC     = TGTOT_;
double VMAXEFF = VMAX_IND * GC3;
double SATFRAC = VMAXEFF / (KM_IND + TGC);   // mg/h per mg/dL of TRL species
double KREM    = K_REM0 * GC3;               // remnant hepatic uptake (1/h)

double R_LPL_CM = KLPL * CM_TG;              // mg/h through limb 1
double R_LPL_VL = KLPL * VLDL_TG;
double R_SAT_CM = SATFRAC * C_CM_;           // mg/h through limb 2
double R_SAT_VL = SATFRAC * SEL_VL * C_VL_;
double R_RES_CM = CL_RES * C_CM_;
double R_RES_VL = CL_RES * C_VL_;

dxdt_CM_TG      = IN_CM - R_LPL_CM - R_SAT_CM - R_RES_CM;
dxdt_CMREM_TG   = R_LPL_CM * FRAC_REM - KREM * CMREM_TG;
dxdt_VLDL_TG    = VSEC   - R_LPL_VL - R_SAT_VL - R_RES_VL;
dxdt_VLDLREM_TG = R_LPL_VL * FRAC_REM - KREM * VLDLREM_TG;

// ============================ FFA, liver, adipose ===========================
double FFA_FROM_LPL = (R_LPL_CM + R_LPL_VL) * (1.0 - FRAC_REM);
dxdt_FFA      = FFA_FROM_LPL + IN_MCT + ADIP_MOB
                - K_FFA_LIV*FFA - K_FFA_OUT*FFA;
dxdt_LIVER_TG = K_FFA_LIV*FFA + K_DNL*(1.0 + 0.8*ALC) - VSEC - K_FAOX*LIVER_TG;
dxdt_ADIP_TG  = K_FFA_OUT*FFA*0.55 - ADIP_MOB;

// ============================ pancreas ======================================
dxdt_PANC_FFA = K_PFFA * C_CM_ - K_PFFAOUT * PANC_FFA;
double EXCESS = PANC_FFA - PFFA_THR;  if(EXCESS < 0.0) EXCESS = 0.0;
dxdt_ACINAR   = K_INJ * (EXCESS/1000.0) - K_REP * ACINAR;
dxdt_PANC_IL6 = K_IL6 * ACINAR - K_IL6OUT * PANC_IL6;
dxdt_CUMIL6   = PANC_IL6;

// ============================ convex AP hazard ==============================
double RT     = pow(TGC/TG50_AP, HILL_AP);
double LAM_YR = LAM_MAX * RT/(1.0 + RT) * (1.0 + PRIME_AP*ACINAR);
dxdt_CUMHAZ   = LAM_YR / 8760.0;

// ============================ exposure metrics ==============================
dxdt_TAT880 = 1.0/(1.0 + exp(-(TGC - TG_THR)/25.0));
dxdt_TGAUC  = TGC;

// ============================ other manifestations ==========================
double XSTIM = TGC - 1500.0;  if(XSTIM < 0.0) XSTIM = 0.0;
dxdt_XANTH  = KX_IN * (XSTIM/1000.0) - KX_OUT * XANTH;
dxdt_HSM    = KH_IN * (TGC/1000.0)   - KH_OUT * HSM;
dxdt_FOG    = KF_IN * (K_VISC*TGC/1000.0) - KF_OUT * FOG;

// ============================ safety ========================================
double KIN_PLT = KOUT_PLT * PLT0;
double EPLT    = ASOSYS/(IC50_PLT + ASOSYS);
dxdt_PLT = KIN_PLT*(1.0 - IPLT_S*EPLT) - KOUT_PLT*PLT*(1.0 + KDES_PLT*EPLT);

double KIN_ALT = KOUT_ALT * ALT0;
dxdt_ALT = KIN_ALT*(1.0 + EALT_ASO*ASOLIV/(EC50_ALT + ASOLIV)) - KOUT_ALT*ALT;

// ============================ drug PK =======================================
dxdt_VOL_SC  = -KA_VOL*VOL_SC;
dxdt_VOL_CP  =  KA_VOL*VOL_SC*F_VOL - (KL_VOL + KO_VOL + KE_VOL)*VOL_CP;
dxdt_VOL_LIV =  KL_VOL*VOL_CP - KOUTL_VOL*VOL_LIV;
dxdt_VOL_SYS =  KO_VOL*VOL_CP - KOUTS_VOL*VOL_SYS;

dxdt_OLE_SC  = -KA_OLE*OLE_SC;
dxdt_OLE_CP  =  KA_OLE*OLE_SC*F_OLE - (KL_OLE + KO_OLE + KE_OLE)*OLE_CP;
dxdt_OLE_LIV =  KL_OLE*OLE_CP - KOUTL_OLE*OLE_LIV;
dxdt_OLE_SYS =  KO_OLE*OLE_CP - KOUTS_OLE*OLE_SYS;

dxdt_PLO_SC   = -KA_PLO*PLO_SC;
dxdt_PLO_CP   =  KA_PLO*PLO_SC*F_PLO - (KL_PLO + KO_PLO + KE_PLO)*PLO_CP;
dxdt_PLO_LIV  =  KL_PLO*PLO_CP - KOUTL_PLO*PLO_LIV - KRISC_IN*PLO_LIV;
dxdt_PLO_RISC =  KRISC_IN*PLO_LIV - KRISC_OUT*PLO_RISC;

dxdt_EVI_C = -(CL_EVI/V1_EVI)*EVI_C - (Q_EVI/V1_EVI)*EVI_C + (Q_EVI/V2_EVI)*EVI_P;
dxdt_EVI_P =  (Q_EVI/V1_EVI)*EVI_C - (Q_EVI/V2_EVI)*EVI_P;

dxdt_FIB_GUT = -KA_FIB*FIB_GUT;
dxdt_FIB_C   =  KA_FIB*FIB_GUT*F_FIB - KE_FIB*FIB_C;

$TABLE
double TG      = TGTOT_;
double TG_MMOL = TG/88.57;
double CM_C    = C_CM_;
double VLDL_C  = C_VL_;
double REM_C   = C_CMR_ + C_VLR_;
double APOC3   = APOC3_P;
double CEVI2   = EVI_C/V1_EVI;
double ANGF    = ANGPTL3*(1.0 - IMAX_EVI*CEVI2/(IC50_EVI + CEVI2));
double CLLPL   = CL_LPL_MAX * LPL_ACT * A_C2 * A_A5
                 /(1.0 + APOC3_P/IC50_C3L)/(1.0 + ANGF/IC50_ANG);
double GC3O    = 1.0 + IMAX_C3 * fmax(0.0, 1.0 - APOC3_P/APOC3_REF);
double SATO    = VMAX_IND*GC3O/(KM_IND + TG);
double FLUX1   = (CLLPL/VP)*(CM_TG + VLDL_TG);
double FLUX2   = SATO*(C_CM_ + SEL_VL*C_VL_);
double FLUX3   = CL_RES*(C_CM_ + C_VL_);
double HAZ_YR  = LAM_MAX*pow(TG/TG50_AP, HILL_AP)/(1.0 + pow(TG/TG50_AP, HILL_AP))
                 *(1.0 + PRIME_AP*ACINAR);
double PROB_AP = 1.0 - exp(-CUMHAZ);
double VISCI   = 1.0 + K_VISC*TG/1000.0;
double LACT    = (TG > 1000.0) ? 1.0 : 0.0;
double ABOVE   = (TG > TG_THR) ? 1.0 : 0.0;

$CAPTURE @annotated
TG      : Total plasma triglyceride (mg/dL)
TG_MMOL : Total plasma triglyceride (mmol/L)
CM_C    : Chylomicron-TG concentration (mg/dL)
VLDL_C  : VLDL-TG concentration (mg/dL)
REM_C   : Remnant-TG concentration (mg/dL)
APOC3   : Plasma apoC-III (mg/dL)
CLLPL   : Instantaneous LPL clearance (dL/h)
FLUX1   : TRL flux through limb 1 LPL (mg/h)
FLUX2   : TRL flux through limb 2 LPL-independent (mg/h)
FLUX3   : TRL flux through residual scavenging (mg/h)
HAZ_YR  : Instantaneous AP hazard (events per patient-year)
PROB_AP : Cumulative probability of at least one AP event (-)
VISCI   : Relative plasma viscosity index (-)
LACT    : Lactescent plasma indicator (0/1)
ABOVE   : TG above 880 mg/dL indicator (0/1)
CEVI2   : Evinacumab plasma concentration (mg/L)
'

mod <- mcode("fcs", code)

## =============================================================================
##  CALIBRATION NOTES
## -----------------------------------------------------------------------------
##  Every number below is a target the structural model must hit before any
##  scenario is believed.  Sources are indexed in fcs_references.md.
##
##  A. NORMAL PHYSIOLOGY  (F_GENO = 1, APOC3_REF = 10, 100 g fat/day)
##     - fasting TG                        ~110-135 mg/dL   (target < 150)
##     - >95 percent of TRL flux carried by limb 1 (FCS_limb_decomposition)
##
##  B. UNTREATED LPL-NULL FCS  (F_GENO = 0, APOC3_REF = 25)
##     - 10 g LCT/day  ->  ~600 mg/dL     near the clinical target
##     - 20 g LCT/day  ->  ~1450 mg/dL    typical "adherent" FCS patient
##     - 60 g LCT/day  ->  ~9500 mg/dL    non-adherent; matches the reported
##                                         clinical range of 1000-10000+ mg/dL
##     - chylomicron fraction ~80 percent of total TG (apoB48-dominant)
##     - post-heparin LPL activity < 5 percent of normal is the diagnostic
##       anchor (Brahm & Hegele 2015; Moulin 2018 FCS Score)
##
##  C. GENOTYPE GRADIENT  (fixed 20 g/day)
##     - F_GENO 0.00 -> ~1450 mg/dL   FCS
##     - F_GENO 0.05 -> ~ 390 mg/dL   severe HTG but not FCS
##     - F_GENO 0.20 -> ~ 140 mg/dL   mild HTG
##       reproducing the clinical rule that the phenotype collapses only below
##       roughly 5-10 percent residual lipase activity (Rahalkar 2009;
##       Hegele 2018).  The disease is a cliff in genotype space as well.
##
##  D. VOLANESORSEN 300 mg SC weekly (APPROACH, Witztum NEJM 2019, n = 66)
##     - apoC-III  -76 +/- 8 percent at 3 months  -> model -78 percent
##     - fasting TG -77 percent at 3 months       -> model -78 percent
##     - platelets < 100 x10^9/L in about 47 percent of patients
##       -> model platelet nadir about 115 x10^9/L on weekly dosing
##
##  E. OLEZARSEN 80 mg SC monthly (Balance, Stroes NEJM 2024, n = 66)
##     - apoC-III about -74 percent
##     - TG -43.5 percent versus placebo +0.8 percent at MONTH 6.  That
##       6-month placebo-adjusted median is depressed by between-subject
##       variability and by incomplete steady state; the open-label MONTH 12
##       value of -73.7 percent is what this deterministic model reproduces,
##       and the model is calibrated to the month-12 number.
##     - acute pancreatitis 1 event on olezarsen versus 11 on placebo
##     - no clinically meaningful platelet decline -> model platelets > 240
##
##  F. PLOZASIRAN 25 mg SC every 3 months (PALISADE, Watts NEJM 2024/2025,
##     n = 75)
##     - apoC-III about -77 percent, TG about -80 percent at month 10
##     - acute pancreatitis events reduced about 83 percent
##       (1.4 percent versus 20 percent of patients)
##     - the RISC effect compartment (half-life about 45 d) is what produces
##       the flat inter-dose PD despite quarterly dosing
##
##  G. FIBRATE AND OMEGA-3 IN LPL-NULL FCS
##     - real-world TG reduction < 10-20 percent (Gaudet 2014; Falko 2018)
##     - the model gives about -19 percent, ENTIRELY from the apoC-III
##       transcriptional term and the omega-3 VLDL-secretion term.  Limb 1
##       contributes exactly 0.0 mg/h in every LPL-null arm.
##
##  H. ACUTE PANCREATITIS HAZARD
##     - placebo AP rate in the FCS trials is about 0.20-0.30 events per
##       patient-year at a median TG of about 2000-2600 mg/dL
##     - the model gives about 0.13 events/py at 20 g fat/day (TG ~1450) and
##       about 0.25 events/py at 30 g fat/day (TG ~2400), bracketing the
##       observed trial rates
##     - hazard is essentially absent below the 880 mg/dL (10 mmol/L) threshold
##     - the Hill exponent of 1.7 makes the hazard CONVEX in TG.  That
##       convexity, and not a larger mean TG effect, is where the observed
##       83-88 percent event reductions come from.
## =============================================================================

## ------------------------- dosing helpers ------------------------------------

#' Three meals a day of a given daily long-chain fat intake
meals <- function(g_per_day, days, start_h = 1) {
  ev(amt = g_per_day/3, cmt = "FAT_GUT", ii = 8, addl = days*3 - 1,
     time = start_h, evid = 1)
}

#' A single high-fat event (restaurant meal, holiday, non-adherence)
fat_binge <- function(g, at_day) {
  ev(amt = g, cmt = "FAT_GUT", time = at_day*24 + 13, evid = 1)
}

volanesorsen <- function(start_day, weeks, mg = 300)
  ev(amt = mg, cmt = "VOL_SC", time = start_day*24, ii = 168, addl = weeks - 1)

olezarsen <- function(start_day, months, mg = 80)
  ev(amt = mg, cmt = "OLE_SC", time = start_day*24, ii = 28*24, addl = months - 1)

plozasiran <- function(start_day, n_doses, mg = 25)
  ev(amt = mg, cmt = "PLO_SC", time = start_day*24, ii = 90*24, addl = n_doses - 1)

evinacumab <- function(start_day, n_doses, mg_per_kg = 15, wt = 70)
  ev(amt = mg_per_kg*wt, cmt = "EVI_C", time = start_day*24, ii = 28*24,
     addl = n_doses - 1)

fenofibrate <- function(start_day, days, mg = 145)
  ev(amt = mg, cmt = "FIB_GUT", time = start_day*24, ii = 24, addl = days - 1)

## ------------------------- genotype / population presets ---------------------
geno <- list(
  fcs_null    = list(F_GENO = 0.00, APOC3_REF = 25),
  fcs_apoc2   = list(F_GENO = 1.00, APOC3_REF = 25, A_C2 = 0),
  fcs_gpihbp1 = list(F_GENO = 0.00, APOC3_REF = 24),
  partial_05  = list(F_GENO = 0.05, APOC3_REF = 20),
  partial_20  = list(F_GENO = 0.20, APOC3_REF = 16),
  mcs         = list(F_GENO = 0.15, APOC3_REF = 22, ALC = 1),
  healthy     = list(F_GENO = 1.00, APOC3_REF = 10)
)

#' Shift an event table by the burn-in period
.shift_ev <- function(e, burn_days) {
  d <- as.data.frame(e)
  d$time <- d$time + burn_days*24
  as.ev(d)
}

#' Run one arm.  `burn` days of untreated diet are simulated first so that
#' every arm starts from its own true steady state; day 0 of the returned
#' data frame is therefore the moment treatment begins.
run_arm <- function(pars = list(), events = NULL, days = 365, burn = 120,
                    fat = 20, delta = 2) {
  m <- mod
  if(length(pars)) m <- do.call(mrgsolve::param, c(list(m), pars))
  ee <- meals(fat, days = burn + days)
  if(!is.null(events)) ee <- ee + .shift_ev(events, burn)
  out <- mrgsim_d(m, data = as.data.frame(ee), end = (burn + days)*24,
                  delta = delta, rtol = 1e-6, atol = 1e-8, maxsteps = 1e6)
  d <- as.data.frame(out)
  d$day <- d$time/24 - burn
  d
}

#' Mean of a variable over the untreated window just before day 0
.baseline <- function(d, var, win = 30) mean(d[[var]][d$day > -win & d$day <= 0])
#' Mean of a variable over the last `win` days of the arm
.final    <- function(d, var, win = 30) mean(d[[var]][d$day > max(d$day) - win])

## =============================================================================
##  SCENARIO 1 — natural history of untreated LPL-null FCS
## =============================================================================
FCS_scenario_natural_history <- function(fat = 20, days = 365)
  run_arm(geno$fcs_null, events = NULL, days = days, fat = fat)

## =============================================================================
##  SCENARIO 2 — the saturation cliff
##  The quantitative statement that "< 20 g/day" is an edge, not a slope.
## =============================================================================
FCS_saturation_curve <- function(fats = c(5,10,15,20,25,30,40,50,60,80,100),
                                 genotypes = c("fcs_null","partial_05",
                                               "partial_20","healthy")) {
  res <- list()
  for(g in genotypes) for(f in fats) {
    d  <- run_arm(geno[[g]], days = 40, burn = 150, fat = f, delta = 2)
    ss <- subset(d, day > 25)
    res[[paste(g,f)]] <- data.frame(
      genotype = g, fat_g_day = f,
      TG_mean = mean(ss$TG), TG_fasting = min(ss$TG), TG_peak = max(ss$TG),
      pct_time_above_880 = 100*mean(ss$TG > 880),
      AP_per_year = mean(ss$HAZ_YR))
  }
  out <- do.call(rbind, res); rownames(out) <- NULL
  cat("\n=== SATURATION CURVE: steady-state TG versus dietary long-chain fat ===\n")
  cat("The hyperbola, not the line, is why the diet threshold behaves like a\n",
      "cliff.  Compare TG(60 g)/TG(20 g) within each genotype.\n\n", sep = "")
  print(out, digits = 4)
  invisible(out)
}

## =============================================================================
##  SCENARIO 3 — genotype gradient (why 5 percent residual lipase is not FCS)
## =============================================================================
FCS_genotype_gradient <- function(fat = 20,
                                  fracs = c(0, 0.01, 0.02, 0.05, 0.10,
                                            0.20, 0.45, 1.00)) {
  res <- lapply(fracs, function(fg) {
    d  <- run_arm(list(F_GENO = fg, APOC3_REF = 25), days = 40, burn = 150,
                  fat = fat, delta = 2)
    ss <- subset(d, day > 25)
    data.frame(residual_LPL = fg, TG_mean = mean(ss$TG),
               TG_fasting = min(ss$TG),
               limb1_mgh = mean(ss$FLUX1), limb2_mgh = mean(ss$FLUX2),
               AP_per_year = mean(ss$HAZ_YR))
  })
  out <- do.call(rbind, res)
  cat("\n=== GENOTYPE GRADIENT at ", fat, " g fat/day ===\n", sep = "")
  cat("The phenotype collapses only below roughly 5-10 percent residual\n",
      "lipase activity: FCS is a cliff in genotype space as well as in diet.\n\n",
      sep = "")
  print(out, digits = 4, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  SCENARIO 4 — conventional therapy (fenofibrate + omega-3) in LPL-null FCS
## =============================================================================
FCS_scenario_conventional <- function(days = 180)
  run_arm(c(geno$fcs_null, list(OM3_ON = 1)),
          events = fenofibrate(0, days), days = days)

## =============================================================================
##  SCENARIO 5 — volanesorsen 300 mg SC weekly (APPROACH-like)
## =============================================================================
FCS_scenario_volanesorsen <- function(days = 365, mg = 300)
  run_arm(geno$fcs_null,
          events = volanesorsen(0, weeks = ceiling(days/7), mg = mg), days = days)

## =============================================================================
##  SCENARIO 6 — olezarsen 50 or 80 mg SC monthly (Balance-like)
## =============================================================================
FCS_scenario_olezarsen <- function(days = 365, mg = 80)
  run_arm(geno$fcs_null,
          events = olezarsen(0, months = ceiling(days/28), mg = mg), days = days)

## =============================================================================
##  SCENARIO 7 — plozasiran 25 mg SC every 3 months (PALISADE-like)
## =============================================================================
FCS_scenario_plozasiran <- function(days = 365, mg = 25)
  run_arm(geno$fcs_null,
          events = plozasiran(0, n_doses = ceiling(days/90), mg = mg), days = days)

## =============================================================================
##  SCENARIO 8 — evinacumab: an ANGPTL3 drug is a LIMB-1 drug
##  Structural prediction: no effect when LPL is null, clear effect when
##  residual lipase exists.  This is the cleanest falsifiable claim here.
## =============================================================================
FCS_scenario_evinacumab <- function(days = 180) {
  n <- ceiling(days/28)
  list(lpl_null = run_arm(geno$fcs_null,   events = evinacumab(0, n), days = days),
       partial  = run_arm(geno$partial_05, events = evinacumab(0, n), days = days),
       mcs      = run_arm(geno$mcs, events = evinacumab(0, n), days = days,
                          fat = 70))
}

## =============================================================================
##  SCENARIO 9 — dietary non-adherence: one 60 g fat binge, on and off drug
## =============================================================================
FCS_scenario_binge <- function(days = 60, binge_day = 30, binge_g = 60,
                               on_drug = TRUE) {
  e <- fat_binge(binge_g, binge_day)
  if(on_drug) e <- e + plozasiran(0, n_doses = 1)
  run_arm(geno$fcs_null, events = e, days = days, delta = 0.5)
}

## =============================================================================
##  SCENARIO 10 — pregnancy in FCS, the highest-risk clinical situation
## =============================================================================
FCS_scenario_pregnancy <- function(days = 280, fat = 15) {
  d <- run_arm(c(geno$fcs_null, list(PREG = 1)), days = days, fat = fat)
  d$trimester <- cut(d$day, c(-Inf, 0, 91, 189, Inf),
                     labels = c("pre","T1","T2","T3"))
  d
}

## =============================================================================
##  ANALYSIS 1 — limb decomposition: the central arithmetic of the disease
## =============================================================================
FCS_limb_decomposition <- function(days = 180) {
  arms <- list(
    "untreated FCS"       = FCS_scenario_natural_history(days = days),
    "fenofibrate + w-3"   = FCS_scenario_conventional(days = days),
    "evinacumab 15 mg/kg" = FCS_scenario_evinacumab(days = days)$lpl_null,
    "volanesorsen 300 qw" = FCS_scenario_volanesorsen(days = days),
    "olezarsen 80 q4w"    = FCS_scenario_olezarsen(days = days),
    "plozasiran 25 q12w"  = FCS_scenario_plozasiran(days = days),
    "healthy control"     = run_arm(geno$healthy, days = days, fat = 100))
  tab <- do.call(rbind, lapply(names(arms), function(n) {
    d <- arms[[n]]
    f1 <- .final(d, "FLUX1"); f2 <- .final(d, "FLUX2"); f3 <- .final(d, "FLUX3")
    data.frame(arm = n, TG = .final(d, "TG"), apoC3 = .final(d, "APOC3"),
               limb1_LPL = f1, limb2_indep = f2, residual = f3,
               pct_limb1 = 100*f1/(f1 + f2 + f3))
  }))
  cat("\n=== LIMB DECOMPOSITION (mean over the last 30 days) ===\n")
  cat("Read the limb1_LPL column.  In every LPL-null arm it is 0.0 mg/h and no\n",
      "drug moves it, because every limb-1 drug multiplies a zero.  All of the\n",
      "therapeutic action in FCS is in the limb2_indep column.\n\n", sep = "")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  ANALYSIS 2 — Jensen gap: what a fasting TG measurement cannot see.
##  Because the hazard is convex, the mean of the hazard exceeds the hazard of
##  the mean.  This is the mathematical form of "the peaks, not the average".
## =============================================================================
FCS_jensen_gap <- function(fats = c(10, 20, 40, 60)) {
  p <- as.list(param(mod))
  haz <- function(tg) p$LAM_MAX*(tg/p$TG50_AP)^p$HILL_AP /
                      (1 + (tg/p$TG50_AP)^p$HILL_AP)
  res <- lapply(fats, function(f) {
    d  <- run_arm(geno$fcs_null, days = 30, burn = 150, fat = f, delta = 0.25)
    ss <- subset(d, day > 20)
    data.frame(fat_g_day = f,
               TG_fasting = min(ss$TG), TG_mean = mean(ss$TG),
               TG_peak = max(ss$TG),
               haz_of_fasting = haz(min(ss$TG)),
               haz_of_mean    = haz(mean(ss$TG)),
               mean_of_haz    = mean(haz(ss$TG)),
               jensen_gap_pct = 100*(mean(haz(ss$TG))/haz(mean(ss$TG)) - 1),
               fasting_underestimate_pct =
                 100*(mean(haz(ss$TG))/haz(min(ss$TG)) - 1))
  })
  out <- do.call(rbind, res)
  cat("\n=== JENSEN GAP: convex hazard meets postprandial excursions ===\n")
  cat("mean_of_haz > haz_of_mean is Jensen's inequality made clinical.  The\n",
      "last column is how badly a FASTING triglyceride under-reports the real\n",
      "annual risk in a patient who has no true fasting state.\n\n", sep = "")
  print(out, digits = 4, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 3 — trial ledger: percent TG change versus what patients care about
## =============================================================================
FCS_trial_ledger <- function(days = 365) {
  arms <- list(
    "placebo (diet only)" = FCS_scenario_natural_history(days = days),
    "fenofibrate + w-3"   = FCS_scenario_conventional(days = days),
    "evinacumab 15 mg/kg" = FCS_scenario_evinacumab(days = days)$lpl_null,
    "volanesorsen 300 qw" = FCS_scenario_volanesorsen(days = days),
    "olezarsen 50 q4w"    = FCS_scenario_olezarsen(days = days, mg = 50),
    "olezarsen 80 q4w"    = FCS_scenario_olezarsen(days = days),
    "plozasiran 25 q12w"  = FCS_scenario_plozasiran(days = days))
  tab <- do.call(rbind, lapply(names(arms), function(n) {
    d   <- arms[[n]]
    b   <- .baseline(d, "TG"); b3 <- .baseline(d, "APOC3")
    m3  <- mean(d$TG[d$day > 80  & d$day <= 90])
    m6  <- mean(d$TG[d$day > 170 & d$day <= 180])
    m12 <- .final(d, "TG")
    yr  <- max(d$day)/365
    on  <- subset(d, day >= 0)
    data.frame(arm = n, TG_baseline = b,
               dTG_m3_pct  = 100*(m3/b  - 1),
               dTG_m6_pct  = 100*(m6/b  - 1),
               dTG_m12_pct = 100*(m12/b - 1),
               apoC3_m12_pct = 100*(.final(d, "APOC3")/b3 - 1),
               days_above_880 = (max(on$TAT880) - min(on$TAT880))/24,
               AP_events_yr = (max(on$CUMHAZ) - min(on$CUMHAZ))/yr,
               PLT_nadir = min(on$PLT), ALT_peak = max(on$ALT),
               xanthoma_m12 = .final(d, "XANTH"))
  }))
  tab$AP_reduction_pct <- 100*(1 - tab$AP_events_yr/tab$AP_events_yr[1])
  cat("\n=== TRIAL LEDGER (", days, " days) ===\n", sep = "")
  cat("Compare dTG_m12_pct with AP_reduction_pct.  The second is always the\n",
      "larger number, because the hazard is convex in TG: deleting the upper\n",
      "tail of the distribution removes disproportionately more risk than the\n",
      "mean percent change suggests.  That is the argument for scoring these\n",
      "drugs on days_above_880 rather than on percent TG.\n\n", sep = "")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  ANALYSIS 4 — time to threshold and inter-dose stability
## =============================================================================
FCS_threshold_time <- function(days = 365) {
  arms <- list("volanesorsen 300 qw" = FCS_scenario_volanesorsen(days = days),
               "olezarsen 80 q4w"    = FCS_scenario_olezarsen(days = days),
               "plozasiran 25 q12w"  = FCS_scenario_plozasiran(days = days))
  tab <- do.call(rbind, lapply(names(arms), function(n) {
    d <- subset(arms[[n]], day >= 0)
    below <- d$TG < 880
    data.frame(arm = n,
               days_to_TG_below_880 = if(any(below)) d$day[which(below)[1]] else NA,
               pct_time_below_880 = 100*mean(below),
               TG_trough = min(d$TG),
               TG_peak_after_day30 = max(d$TG[d$day > 30]),
               apoC3_nadir_pct = 100*(min(d$APOC3)/max(d$APOC3) - 1),
               PLT_nadir = min(d$PLT), ALT_peak = max(d$ALT))
  }))
  cat("\n=== TIME TO THRESHOLD AND INTER-DOSE STABILITY ===\n\n")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  ANALYSIS 5 — the GalNAc dividend: the safety trade that decided the class
## =============================================================================
FCS_galnac_dividend <- function(days = 365) {
  v <- FCS_scenario_volanesorsen(days = days)
  o <- FCS_scenario_olezarsen(days = days)
  p <- FCS_scenario_plozasiran(days = days)
  f <- function(d) c(apoC3 = 100*(.final(d,"APOC3")/.baseline(d,"APOC3") - 1),
                     TG    = 100*(.final(d,"TG")/.baseline(d,"TG") - 1),
                     plt   = min(d$PLT[d$day >= 0]),
                     tlow  = 100*mean(d$PLT[d$day >= 0] < 100))
  tab <- data.frame(
    arm = c("volanesorsen (naked PS-ASO)", "olezarsen (GalNAc3-ASO)",
            "plozasiran (GalNAc-siRNA)"),
    annual_dose_mg = c(300*52, 80*13, 25*4),
    injections_per_year = c(52, 13, 4),
    t(sapply(list(v,o,p), f)))
  names(tab)[4:7] <- c("apoC3_pct","TG_pct","PLT_nadir","pct_time_PLT_lt_100")
  cat("\n=== THE GalNAc DIVIDEND ===\n")
  cat("Identical target, near-identical knockdown, 20-30x less systemic\n",
      "exposure.  The platelet signal that stalled volanesorsen was a PK\n",
      "problem, not a pharmacology problem.\n\n", sep = "")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  ANALYSIS 6 — diet x drug: both levers act on the same saturable limb
## =============================================================================
FCS_diet_drug_interaction <- function(days = 180) {
  grid <- expand.grid(fat = c(10, 20, 40, 60), drug = c("none","plozasiran"),
                      stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    e <- if(grid$drug[i] == "plozasiran") plozasiran(0, 2) else NULL
    d <- run_arm(geno$fcs_null, events = e, days = days, fat = grid$fat[i],
                 delta = 4)
    ss <- subset(d, day > days - 30)
    data.frame(fat_g_day = grid$fat[i], drug = grid$drug[i],
               TG = mean(ss$TG), pct_above_880 = 100*mean(ss$TG > 880),
               AP_per_year = mean(ss$HAZ_YR))
  }))
  cat("\n=== DIET x DRUG ===\n")
  cat("Both levers act on the SAME saturable limb, so their effects are not\n",
      "additive in TG - they are multiplicative in the distance from the knee.\n",
      "A drug does not license a 60 g/day diet; it moves the cliff edge out.\n\n",
      sep = "")
  print(out, digits = 4, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 7 — pregnancy risk arc
## =============================================================================
FCS_pregnancy_risk <- function() {
  d <- FCS_scenario_pregnancy()
  tab <- do.call(rbind, lapply(c("T1","T2","T3"), function(tr) {
    s <- subset(d, trimester == tr)
    data.frame(trimester = tr, TG_mean = mean(s$TG), TG_peak = max(s$TG),
               pct_time_above_880 = 100*mean(s$TG > 880),
               AP_hazard_per_year = mean(s$HAZ_YR),
               cumulative_AP_prob = max(s$PROB_AP))
  }))
  cat("\n=== PREGNANCY IN FCS ===\n")
  cat("Oestrogen raises VLDL secretion and lowers what little LPL there is.\n",
      "Risk is back-loaded into the third trimester, which is why apheresis\n",
      "programmes are planned from about week 24.\n\n", sep = "")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  ANALYSIS 8 — binge kinetics: how long does one bad meal cost?
## =============================================================================
FCS_binge_cost <- function(binge_g = 60) {
  off <- FCS_scenario_binge(on_drug = FALSE, binge_g = binge_g)
  on  <- FCS_scenario_binge(on_drug = TRUE,  binge_g = binge_g)
  f <- function(d, lab) {
    w <- subset(d, day >= 29 & day <= 40)
    base <- mean(d$TG[d$day > 20 & d$day < 29])
    data.frame(arm = lab, pre_binge_TG = base, peak_TG = max(w$TG),
               peak_rise_pct = 100*(max(w$TG)/base - 1),
               hours_above_880 = sum(w$TG > 880)*0.5,
               added_AP_prob = max(w$PROB_AP) - min(w$PROB_AP))
  }
  tab <- rbind(f(off, "no drug"), f(on, "on plozasiran"))
  cat("\n=== COST OF ONE 60 g FAT MEAL ===\n\n")
  print(tab, digits = 4, row.names = FALSE)
  invisible(tab)
}

## =============================================================================
##  PLOTTING (base graphics only; no extra dependencies)
## =============================================================================
FCS_plot_overview <- function(d, main = "FCS QSP simulation") {
  op <- par(mfrow = c(2,3), mar = c(4,4,2.5,1), oma = c(0,0,2,0))
  on.exit(par(op))
  plot(d$day, d$TG, type = "l", col = "#d84315", lwd = 2,
       xlab = "day", ylab = "TG (mg/dL)", main = "Plasma triglyceride")
  abline(h = 880, lty = 2, col = "grey40")
  plot(d$day, d$APOC3, type = "l", col = "#00897b", lwd = 2,
       xlab = "day", ylab = "apoC-III (mg/dL)", main = "Target engagement")
  plot(d$day, d$FLUX1, type = "l", col = "#d81b60", lwd = 2,
       ylim = range(0, d$FLUX1, d$FLUX2), xlab = "day", ylab = "mg/h",
       main = "Clearance limbs")
  lines(d$day, d$FLUX2, col = "#7e57c2", lwd = 2)
  legend("topright", c("limb 1 (LPL)", "limb 2 (independent)"),
         col = c("#d81b60","#7e57c2"), lwd = 2, bty = "n", cex = 0.8)
  plot(d$day, d$HAZ_YR, type = "l", col = "#bf360c", lwd = 2,
       xlab = "day", ylab = "events/patient-year", main = "AP hazard")
  plot(d$day, d$PLT, type = "l", col = "#5e35b1", lwd = 2,
       xlab = "day", ylab = "10^9/L", main = "Platelets")
  abline(h = 100, lty = 2, col = "grey40")
  plot(d$day, d$XANTH, type = "l", col = "#827717", lwd = 2,
       xlab = "day", ylab = "score", main = "Eruptive xanthoma")
  mtext(main, outer = TRUE, font = 2)
}

FCS_plot_saturation <- function(tab = NULL) {
  if(is.null(tab)) tab <- FCS_saturation_curve()
  op <- par(mar = c(4.5,4.5,3,1)); on.exit(par(op))
  gs <- unique(tab$genotype); cols <- c("#d84315","#f9a825","#43a047","#1565c0")
  plot(NA, xlim = range(tab$fat_g_day), ylim = c(0, max(tab$TG_mean)*1.05),
       xlab = "dietary long-chain fat (g/day)",
       ylab = "steady-state TG (mg/dL)", main = "The saturation cliff")
  abline(h = 880, lty = 2, col = "grey50")
  for(i in seq_along(gs)) {
    s <- tab[tab$genotype == gs[i], ]
    lines(s$fat_g_day, s$TG_mean, col = cols[i], lwd = 2, type = "b", pch = 16)
  }
  legend("topleft", gs, col = cols, lwd = 2, bty = "n")
}

## =============================================================================
##  RUN-ALL
## =============================================================================
FCS_run_all <- function(days = 365) {
  cat("\n##########################################################\n")
  cat("#  Familial Chylomicronemia Syndrome - QSP model run\n")
  cat("##########################################################\n")
  sat <- FCS_saturation_curve()
  gen <- FCS_genotype_gradient()
  lim <- FCS_limb_decomposition(days = 180)
  jen <- FCS_jensen_gap()
  led <- FCS_trial_ledger(days = days)
  thr <- FCS_threshold_time(days = days)
  gal <- FCS_galnac_dividend(days = days)
  dxd <- FCS_diet_drug_interaction()
  prg <- FCS_pregnancy_risk()
  bin <- FCS_binge_cost()
  s   <- sat[sat$genotype == "fcs_null", ]
  cat("\n--- SUMMARY OF THE FOUR CLAIMS ---------------------------\n")
  cat("1. Saturation cliff : TG(60 g)/TG(20 g) in LPL-null FCS = ",
      round(s$TG_mean[s$fat_g_day == 60]/s$TG_mean[s$fat_g_day == 20], 2),
      "x\n", sep = "")
  cat("2. Dead limb        : largest limb-1 flux in any LPL-null arm = ",
      round(max(lim$limb1_LPL[lim$arm != "healthy control"]), 4),
      " mg/h\n", sep = "")
  cat("3. apoC-III works   : best TG reduction achieved with LPL = 0 is ",
      round(-min(led$dTG_m12_pct), 1), "%\n", sep = "")
  cat("4. Convexity pays   : that arm's AP-event reduction is ",
      round(max(led$AP_reduction_pct), 1), "%\n", sep = "")
  invisible(list(saturation = sat, genotype = gen, limbs = lim, jensen = jen,
                 ledger = led, threshold = thr, galnac = gal, dietdrug = dxd,
                 pregnancy = prg, binge = bin))
}

## Uncomment to run everything:
# results <- FCS_run_all()
