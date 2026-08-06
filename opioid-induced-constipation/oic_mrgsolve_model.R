## =============================================================================
##  oic_mrgsolve_model.R
##  Opioid-Induced Constipation (OIC) — Quantitative Systems Pharmacology model
##  오피오이드 유발 변비 — 정량적 시스템 약리학 모델
##
##  51 ODEs · 5 drug classes · 16 treatment scenarios
##
##  PROVENANCE — READ THIS FIRST
##  ---------------------------------------------------------------------------
##  The container this model was built in has NO R RUNTIME.  This file has
##  therefore NEVER BEEN EXECUTED.  Every equation in it was first implemented,
##  integrated and calibrated in `oic_reference_model.py` (Python/scipy, in this
##  same directory), and this file is a term-for-term transcription of that
##  implementation.  The committed run log of the Python version is
##  `oic_reference_output.txt`, and every number quoted in README.md comes from
##  there — not from this file.
##
##  Consequence for the reader: the SCIENCE here has been checked by
##  integration; the R SYNTAX has only been checked by eye.  If you are the
##  first person to run it, expect to fix typos, and please do not treat a
##  failure to compile as evidence about the model.
##
##  THE ONE STRUCTURAL CLAIM
##  ---------------------------------------------------------------------------
##  A PAMORA's therapeutic window is not a potency.  It is the ratio of two
##  receptor occupancies computed at the same plasma concentration:
##
##       SI = OCC_antagonist(enteric plexus) / OCC_antagonist(CNS)
##
##  The enteric plexus is outside the blood-brain barrier and the CNS is not.
##  Both hold the same receptor competing with the same agonist.  So SI reduces,
##  in the low-occupancy limit, to 1/Kp_uu of the antagonist — corrected by the
##  agonist's own brain partitioning.  "Peripheral selectivity" is not asserted
##  anywhere in this file; it is a division.
##
##  WHAT WAS FITTED (three parameters, one per drug)
##  ---------------------------------------------------------------------------
##    KIGUT_naloxegol   = 2.685   nM  <- KODIAC-04/05 25 mg arm, mean SBM/wk 4.3
##    KIGUT_naldemedine = 0.00821 nM  <- COMPOSE-1/2 0.2 mg arm, mean SBM/wk 5.0
##    KIGUT_mntx        = 0.795   nM  <- Thomas 2008, placebo-subtracted 4-h
##                                       laxation difference of 33 pp
##  Everything else is a literature value or a stated structural assumption.
##  The three fitted values are each 3–40x more potent than the drug's own
##  published MOR binding Ki; the model does NOT explain that spread and says so
##  (see analysis G in the Python run log).  The SELECTIVITY results do not
##  depend on them: KIGUT cancels out of every ratio.
##
##  Units: time h · drug mg (peptides ug) · colonic solids g · water mL ·
##         osmoles mmol · receptor concentrations nM
## =============================================================================

library(mrgsolve)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

oic_code <- '
$PARAM @annotated
// ---------------- opioid (index drug: oxycodone) ----------------
MWOP    : 315.4  : opioid molecular weight (g/mol)
FOP     : 0.75   : oral bioavailability
KAOP    : 1.2    : absorption rate constant (1/h)
CLOP    : 45.0   : clearance (L/h)
V1OP    : 90.0   : central volume (L)
V2OP    : 120.0  : peripheral volume (L)
QOP     : 40.0   : intercompartmental clearance (L/h)
FUOP    : 0.55   : unbound fraction in plasma
KPUUOP  : 3.0    : unbound brain/plasma ratio (oxycodone has ACTIVE BBB influx)
KEOOP   : 0.60   : brain equilibration rate (1/h)
KIOP    : 30.0   : MOR Ki, unbound (nM)
FLUMOP  : 0.25   : fraction of an oral dose reaching the colonic lumen
KLUMOUT : 0.15   : luminal loss rate (1/h)

// ---------------- antagonist (defaults = naloxegol) ----------------
MWANT   : 651.8  : antagonist molecular weight (g/mol)
FANT    : 0.50   : bioavailability
KAANT   : 1.0    : absorption rate constant (1/h)
CLANT   : 22.0   : clearance (L/h)
V1ANT   : 68.0   : central volume (L)
V2ANT   : 100.0  : peripheral volume (L)
QANT    : 15.0   : intercompartmental clearance (L/h)
FUANT   : 0.96   : unbound fraction in plasma
KPUU    : 0.020  : unbound brain/plasma ratio -- THE selectivity parameter
KEOANT  : 0.50   : brain equilibration rate (1/h)
KIANT   : 7.4    : MOR binding Ki, unbound (nM) -- literature, NOT fitted
KIGUT   : 2.6851 : OPERATIONAL enteric potency (nM) -- FITTED (1 of 3)
FLUMANT : 0.35   : fraction of antagonist dose reaching the lumen
ORALANT : 1      : 1 = oral, 0 = subcutaneous (no luminal/pre-systemic route)
FLUMACCX: 1.0    : per-drug multiplier on luminal mucosal access (quaternary=low)

// ---------------- prokinetic (prucalopride, 5-HT4) ----------------
MWPRO   : 367.9  : molecular weight (g/mol)
FPRO    : 0.90   : bioavailability
CLPRO   : 18.0   : clearance (L/h)
V1PRO   : 320.0  : central volume (L)
V2PRO   : 300.0  : peripheral volume (L)
QPRO    : 25.0   : intercompartmental clearance (L/h)
EC50PRO : 8.0    : enteric 5-HT4 EC50 (nM)
EMAXPRO : 0.55   : max fractional rise in enteric cAMP

// ---------------- gut-wall / lumen geometry ----------------
QPORT   : 60.0   : portal blood flow (L/h) -- pre-systemic gut-wall term
VLUM    : 0.35   : effective colonic luminal volume for plexus access (L)
FLUMACC : 0.002  : fraction of luminal conc reaching the myenteric plexus

// ---------------- receptor trafficking (the tolerance asymmetry) ----------
KINRC   : 0.0060 : CNS MOR resensitisation (1/h)
KDESC   : 0.0045 : CNS MOR desensitisation (1/h per unit occupancy)
KINRG   : 0.0060 : enteric MOR resensitisation (1/h)
KDESG   : 0.00060: enteric MOR desensitisation (1/h) -- 7.5x smaller than CNS
KONARR  : 0.012  : beta-arrestin-2 recruitment (1/h)
KOFFARR : 0.010  : beta-arrestin-2 loss (1/h)
BARR    : 1.20   : arrestin amplification of CNS desensitisation

// ---------------- transduction (operational agonism) ----------------
OCC50G  : 0.30   : enteric occupancy giving half-max motility inhibition
HG      : 1.40   : enteric transduction Hill coefficient
OCC50C  : 0.30   : CNS occupancy giving half-max analgesia
HC      : 2.00   : CNS transduction Hill coefficient

// ---------------- enteric signalling ----------------
KSCAMP  : 2.50   : cAMP synthesis rate (1/h)
KDCAMP  : 2.50   : cAMP degradation rate (1/h)
EMAXMOR : 0.82   : max fractional cAMP suppression by enteric MOR
KACH    : 1.50   : ACh release turnover (1/h)
HACH    : 1.60   : cAMP -> ACh exponent
KNO     : 0.30   : NO/VIP tone turnover (1/h)
ANO     : 0.50   : MOR -> inhibitory NO/VIP tone
KTONE   : 1.20   : segmental tone turnover (1/h)
ATONE   : 0.80   : MOR -> non-propulsive segmental tone
NTONE   : 0.85   : tone exponent damping propulsion
ANOVIP  : 0.50   : NO/VIP damping of propagating contractions
ESECMOR : 0.55   : max fractional suppression of SECRETION by enteric MOR

// ---------------- secretagogues ----------------
KCGMP   : 0.50   : cGMP turnover (1/h)
EC50LIN : 55.0   : linaclotide luminal amount for half-max GC-C (ug)
EMAXCGMP: 2.60   : max cGMP rise
KCLC2   : 0.60   : ClC-2 activation turnover (1/h)
EC50LUB : 14.0   : lubiprostone luminal amount for half-max ClC-2 (ug)
EMAXCLC2: 1.00   : max ClC-2 activation
ICLC2ME : 0.35   : methadone ClC-2 block potency (occupancy units)
ECLC2   : 2.30   : ClC-2 -> secretion gain
ECGMP   : 1.05   : cGMP -> CFTR secretion gain

// ---------------- motility ----------------
HAPC0   : 0.28   : baseline high-amplitude propagating contractions (1/h)
KH      : 1.50   : HAPC turnover (1/h)
KMACH   : 0.85   : ACh -> HAPC half-saturation
HDIST   : 0.45   : colonic filling gate exponent
FILLREF : 300.0  : reference total colonic content for the filling gate (g)

// ---------------- colonic content ----------------
SIN     : 1.46   : dry solids entering the caecum (g/h) = 35 g/day
WIN     : 62.5   : water entering the caecum (mL/h) = 1500 mL/day
WBIND   : 1.60   : water bound per g dry solids, UNABSORBABLE (mL/g)
WPEG    : 1.02   : H-bonded water per g PEG 3350 (mL/g)
KABS1   : 2.050  : ascending colon water absorption (1/h)
KABS2   : 0.137  : transverse colon water absorption (1/h)
KABS3   : 0.211  : descending colon water absorption (1/h)
KABS4   : 0.0279 : rectosigmoid water absorption (1/h)
VSEC1   : 3.2    : ascending colon baseline secretion (mL/h)
VSEC2   : 2.2    : transverse colon baseline secretion (mL/h)
VSEC3   : 1.5    : descending colon baseline secretion (mL/h)
VSEC4   : 0.9    : rectosigmoid baseline secretion (mL/h)
KPROP1  : 0.232  : ascending colon propulsion (1/h)
KPROP2  : 0.278  : transverse colon propulsion (1/h)
KPROP3  : 0.510  : descending colon propulsion (1/h)
W50     : 0.680  : stool water fraction at half-max propulsive efficiency
GW      : 3.5    : hydration -> propulsion Hill coefficient
OSMPL   : 0.290  : osmoles per mL isotonic water (mmol/mL)
KFERM   : 0.075  : lactulose fermentation rate (1/h)
FERMAMP : 3.60   : osmole amplification per lactulose molecule fermented
FGAS    : 0.30   : fraction of fermentation osmoles that are gas

// ---------------- defaecation ----------------
KDEF    : 0.415  : maximal urge/defaecation rate (1/h)
VDEF    : 170.0  : rectosigmoid load at half-max urge (g)
HDEF    : 2.60   : urge Hill coefficient
W50D    : 0.550  : anorectal hydration gate midpoint (milder than colonic)
GWD     : 2.50   : anorectal hydration gate Hill coefficient
HHAPC   : 0.50   : HAPC gating exponent on defaecation
AANO    : 3.00   : MOR -> internal anal sphincter tone / blunted RAIR
AVSENS  : 2.50   : MOR-driven rise in the RECTAL URGE THRESHOLD
FEVAC   : 0.85   : rectosigmoid fraction voided per spontaneous event
KRESC   : 0.012  : maximal rescue-laxative-driven BM rate (1/h) ~2/wk
RTHR    : 0.0139 : spontaneous rate below which rescue is taken (1/h) = 1/72 h
FEVACR  : 0.75   : rectosigmoid fraction voided per rescue BM
FEVACR3 : 0.35   : descending colon fraction recruited by rescue
RESCUEON: 1.0    : 1 = trial-protocol rescue laxative allowed, 0 = disabled

// ---------------- symptoms ----------------
KBSFS   : 0.30   : Bristol score turnover (1/h)
WB50    : 0.766  : water fraction giving Bristol 4
HB      : 9.0    : Bristol Hill coefficient
KSTR    : 0.25   : straining turnover (1/h)
W50S    : 0.660  : water fraction at half-max straining relief
HS      : 8.0    : straining Hill coefficient
KDIST   : 0.10   : distension turnover (1/h)
KSYM    : 0.12   : PAC-SYM turnover (1/h)
KQOL    : 0.045  : PAC-QOL turnover (1/h)
DIST0   : 900.0  : colonic content giving distension score 2 (g)

// ---------------- CNS / pain / withdrawal ----------------
PAIN0   : 7.2    : untreated chronic pain NRS
PAINMAX : 6.0    : maximal opioid analgesia (NRS units)
KPAIN   : 0.25   : pain turnover (1/h)
KTOLA   : 0.0016 : counter-adaptation accrual (1/h)
KTOLAR  : 0.0022 : counter-adaptation loss (1/h)
ATOLA   : 1.4    : counter-adaptation weight on pain
APAINGI : 0.45   : abdominal-pain contribution of constipation to NRS
WDMAX   : 34.0   : COWS ceiling used here
WD50    : 0.145  : CNS antagonist occupancy at half-max withdrawal
HWD     : 2.6    : withdrawal Hill coefficient
KWD     : 0.9    : withdrawal turnover (1/h)
DEPEND  : 1.0    : physical dependence multiplier
KNAUS   : 0.30   : nausea turnover (1/h)
ENAUSLUB: 0.9    : lubiprostone systemic nausea
ENAUSDIS: 0.45   : distension-driven nausea
KIMP    : 0.010  : impaction hazard accrual (1/h)

// ---------------- modifiers ----------------
PGPINH  : 1.0    : multiplier on antagonist Kp_uu (P-gp inhibition)
CYP3A4I : 1.0    : divisor on antagonist clearance (CYP3A4 inhibition)
METHADON: 0.0    : 1 = index opioid is methadone (direct ClC-2 block)
FLUIDX  : 1.0    : patient fluid-intake multiplier on WIN
ABSX    : 1.0    : patient colonic water-absorption multiplier

$CMT @annotated
AOP_DEP  : opioid oral depot (mg)
AOP_CEN  : opioid central (mg)
AOP_PER  : opioid peripheral (mg)
COP_BR   : opioid brain unbound concentration (nM)
AOP_LUM  : opioid in colonic lumen (mg)
APAM_DEP : antagonist depot (mg)
APAM_CEN : antagonist central (mg)
APAM_PER : antagonist peripheral (mg)
CPAM_BR  : antagonist brain unbound concentration (nM)
APAM_LUM : antagonist in colonic lumen (mg)
APEG_LUM : PEG 3350 in lumen (mmol)
ALAC_LUM : lactulose in lumen (mmol)
ALUB_LUM : lubiprostone in lumen (ug)
ALIN_LUM : linaclotide in lumen (ug)
APRO_CEN : prucalopride central (mg)
APRO_PER : prucalopride peripheral (mg)
RG_AV    : enteric MOR availability (fraction)
RC_AV    : CNS MOR availability (fraction)
ARR      : beta-arrestin-2 / GRK state
CAMP     : enteric neuron cAMP (relative)
ACH      : ACh release capacity (relative)
NOVIP    : inhibitory NO/VIP tone (relative)
CGMP     : enterocyte cGMP (relative)
CLC2     : ClC-2 activation state
HAPC     : high-amplitude propagating contractions (1/h)
TONE     : segmental muscle tone (relative)
S1       : ascending colon dry solids (g)
S2       : transverse colon dry solids (g)
S3       : descending colon dry solids (g)
S4       : rectosigmoid dry solids (g)
W1       : ascending colon water (mL)
W2       : transverse colon water (mL)
W3       : descending colon water (mL)
W4       : rectosigmoid water (mL)
O1       : ascending colon non-absorbable osmoles (mmol)
O2       : transverse colon non-absorbable osmoles (mmol)
O3       : descending colon non-absorbable osmoles (mmol)
O4       : rectosigmoid non-absorbable osmoles (mmol)
BSFS     : Bristol Stool Form Scale (1-7)
STRAIN   : straining score (0-4)
DIST     : abdominal distension (0-4)
PACSYM   : PAC-SYM symptom score (0-4)
QOL      : PAC-QOL (0-4, higher worse)
PAIN     : pain NRS (0-10)
TOLA     : analgesic counter-adaptation
WD       : withdrawal (COWS)
NAUSEA   : nausea (0-3)
CUM_SBM  : cumulative spontaneous bowel movements
CUM_CSBM : cumulative complete SBMs
CUM_RESC : cumulative rescue-laxative bowel movements
IMPACT   : faecal impaction hazard integral

$GLOBAL
#define HILLF(x, x50, n) ( ((x) <= 0.0) ? 0.0 : pow((x),(n)) / (pow((x),(n)) + pow((x50),(n))) )
#define WFRAC(W, S) ( ((W)+(S) <= 1e-9) ? 0.75 : ((W)/((W)+(S))) )

$MAIN
// initial conditions -- the healthy steady state
RG_AV_0 = 1.0;   RC_AV_0 = 1.0;
CAMP_0  = 1.0;   ACH_0   = 1.0;   NOVIP_0 = 1.0;   CGMP_0 = 1.0;
HAPC_0  = HAPC0; TONE_0  = 1.0;
S1_0 = 11.7;  S2_0 = 10.2;  S3_0 = 5.8;   S4_0 = 27.5;
W1_0 = 46.8;  W2_0 = 36.2;  W3_0 = 18.5;  W4_0 = 82.5;
BSFS_0 = 4.0; PAIN_0 = PAIN0;

$ODE
// ============================ opioid PK ==================================
double Cop_nM = (AOP_CEN / V1OP) / MWOP * 1e6;
double Aop_u  = FUOP * Cop_nM;
double abs_op = KAOP * AOP_DEP;

dxdt_AOP_DEP = -abs_op;
dxdt_AOP_CEN = FOP*abs_op - (CLOP/V1OP)*AOP_CEN - (QOP/V1OP)*AOP_CEN
               + (QOP/V2OP)*AOP_PER;
dxdt_AOP_PER = (QOP/V1OP)*AOP_CEN - (QOP/V2OP)*AOP_PER;
dxdt_COP_BR  = KEOOP * (KPUUOP*Aop_u - COP_BR);
dxdt_AOP_LUM = FLUMOP*abs_op - KLUMOUT*AOP_LUM;

// pre-systemic gut-wall increment = absorption flux / portal blood flow
double gw_op  = FUOP * (FOP*abs_op/QPORT) / MWOP * 1e6;
double lum_op = (AOP_LUM/VLUM) / MWOP * 1e6 * FLUMACC;
double A_gut  = Aop_u + gw_op + lum_op;
double A_cns  = COP_BR;

// ========================== antagonist PK ================================
double CL_ant  = CLANT / CYP3A4I;
double Cpam_nM = (APAM_CEN / V1ANT) / MWANT * 1e6;
double Bp_u    = FUANT * Cpam_nM;
double abs_pam = KAANT * APAM_DEP;

dxdt_APAM_DEP = -abs_pam;
dxdt_APAM_CEN = FANT*abs_pam - (CL_ant/V1ANT)*APAM_CEN - (QANT/V1ANT)*APAM_CEN
                + (QANT/V2ANT)*APAM_PER;
dxdt_APAM_PER = (QANT/V1ANT)*APAM_CEN - (QANT/V2ANT)*APAM_PER;

// P-gp inhibition acts on Kp_uu, NOT on clearance.  That separation is the
// entire point of the drug-interaction analysis.
double kpuu_eff = (KPUU*PGPINH > 1.0) ? 1.0 : KPUU*PGPINH;
dxdt_CPAM_BR  = KEOANT * (kpuu_eff*Bp_u - CPAM_BR);
dxdt_APAM_LUM = FLUMANT*abs_pam*ORALANT - KLUMOUT*APAM_LUM;

double gw_pam  = FUANT * (FANT*abs_pam*ORALANT/QPORT) / MWANT * 1e6;
double lum_pam = (APAM_LUM/VLUM) / MWANT * 1e6 * FLUMACC * FLUMACCX;
double B_gut   = Bp_u + gw_pam + lum_pam;
double B_cns   = CPAM_BR;

// ====================== competitive occupancy ============================
// enteric site uses the FITTED operational potency KIGUT;
// the CNS site uses the LITERATURE binding Ki -- the safety side never
// borrows a fitted number.
double rg_a = A_gut / KIOP;
double rg_b = B_gut / KIGUT;
double den_g = 1.0 + rg_a + rg_b;
double OCCg_ag  = rg_a / den_g;
double OCCg_ant = rg_b / den_g;

double rc_a = A_cns / KIOP;
double rc_b = B_cns / KIANT;
double den_c = 1.0 + rc_a + rc_b;
double OCCc_ag  = rc_a / den_c;
double OCCc_ant = rc_b / den_c;

double SIGg = HILLF(OCCg_ag*RG_AV, OCC50G, HG);
double SIGc = HILLF(OCCc_ag*RC_AV, OCC50C, HC);

// ======================= receptor trafficking ============================
dxdt_RG_AV = KINRG*(1.0 - RG_AV) - KDESG*OCCg_ag*RG_AV;
dxdt_RC_AV = KINRC*(1.0 - RC_AV) - KDESC*OCCc_ag*(1.0 + BARR*ARR)*RC_AV;
dxdt_ARR   = KONARR*OCCc_ag - KOFFARR*ARR;

// ============ prucalopride: Gs on the SAME adenylyl cyclase ==============
double Cpro  = (APRO_CEN/V1PRO) / MWPRO * 1e6;
double E5HT4 = EMAXPRO * Cpro / (Cpro + EC50PRO);
dxdt_APRO_CEN = -(CLPRO/V1PRO)*APRO_CEN - (QPRO/V1PRO)*APRO_CEN
                + (QPRO/V2PRO)*APRO_PER;
dxdt_APRO_PER = (QPRO/V1PRO)*APRO_CEN - (QPRO/V2PRO)*APRO_PER;

// ====================== enteric second messengers ========================
dxdt_CAMP  = KSCAMP*(1.0 - EMAXMOR*SIGg + E5HT4) - KDCAMP*CAMP;
dxdt_ACH   = KACH*(pow(CAMP, HACH) - ACH);
dxdt_NOVIP = KNO*(1.0 + ANO*SIGg - NOVIP);

dxdt_ALUB_LUM = -KLUMOUT*ALUB_LUM;
dxdt_ALIN_LUM = -KLUMOUT*ALIN_LUM;

// methadone blocks ClC-2 DIRECTLY -- a channel block, not a MOR effect.
// This is why lubiprostone underperformed in methadone-treated patients while
// a PAMORA did not.
double clc2_block = (METHADON > 0.5) ? (1.0/(1.0 + OCCg_ag/ICLC2ME)) : 1.0;
double clc2_drive = EMAXCLC2 * HILLF(ALUB_LUM, EC50LUB, 1.0) * clc2_block;
dxdt_CLC2 = KCLC2*(clc2_drive - CLC2);
double cgmp_drive = 1.0 + EMAXCGMP * HILLF(ALIN_LUM, EC50LIN, 1.0);
dxdt_CGMP = KCGMP*(cgmp_drive - CGMP);

// ============================= motility ==================================
double Stot = S1 + S2 + S3 + S4;
double Wtot = W1 + W2 + W3 + W4;
double fill_tot = (Stot + Wtot < 1.0) ? 1.0 : (Stot + Wtot);
double fill_gate = pow(fill_tot/FILLREF, HDIST);
double ach_term  = (ACH/(ACH + KMACH)) / (1.0/(1.0 + KMACH));
double novip_ex  = (NOVIP > 1.0) ? (NOVIP - 1.0) : 0.0;
double hapc_drive = HAPC0 * ach_term * fill_gate / (1.0 + ANOVIP*novip_ex);
dxdt_HAPC = KH*(hapc_drive - HAPC);
dxdt_TONE = KTONE*(1.0 + ATONE*SIGg - TONE);

double hapc_rel = ((HAPC < 1e-6) ? 1e-6 : HAPC) / HAPC0;
double tone_div = pow(((TONE < 1e-3) ? 1e-3 : TONE), NTONE);

// ========================= colonic content ===============================
double w1 = WFRAC(W1, S1);
double w2 = WFRAC(W2, S2);
double w3 = WFRAC(W3, S3);
double w4 = WFRAC(W4, S4);

double kp1 = KPROP1 * hapc_rel * HILLF(w1, W50, GW) / tone_div;
double kp2 = KPROP2 * hapc_rel * HILLF(w2, W50, GW) / tone_div;
double kp3 = KPROP3 * hapc_rel * HILLF(w3, W50, GW) / tone_div;

// enteric MOR suppresses SECRETION as well as motility.  This is the FAST arm
// of acute laxation: without it the model cannot produce methylnaltrexone-like
// laxation within 4 h at ANY enteric potency.
double sec_gain = (1.0 - ESECMOR*SIGg) + ECLC2*CLC2
                  + ECGMP*((CGMP > 1.0) ? (CGMP - 1.0) : 0.0);

dxdt_APEG_LUM = -kp1*APEG_LUM;
double ferm   = KFERM*ALAC_LUM;
dxdt_ALAC_LUM = -ferm - kp1*ALAC_LUM;
double osm_in = ferm*FERMAMP;

// PEG is osmotically active AND holds water by hydrogen bonding; both are
// entered as osmole-equivalents so they share the unabsorbable-water machinery.
double O_extra1 = APEG_LUM*(1.0 + 3.350*WPEG*OSMPL) + ALAC_LUM;

// ============================ defaecation ================================
double load4 = S4 + W4;
double vdef_eff = VDEF * (1.0 + AVSENS*SIGg);
double Rdef = KDEF * HILLF(load4, vdef_eff, HDEF) * HILLF(w4, W50D, GWD)
              * pow(hapc_rel, HHAPC) / (1.0 + AANO*SIGg);
// trial-protocol rescue laxative (permitted after 72 h with no bowel movement).
// Rescue-driven evacuations empty the colon but are NOT counted as SBMs --
// that is the endpoint definition, and the model has to honour it.
double resc_need = 1.0 - HILLF(Rdef, RTHR, 2.5);
double Rresc = KRESC * RESCUEON * resc_need;
double evac4 = FEVAC*Rdef + FEVACR*Rresc;
double evac3 = FEVACR3*Rresc;

double kout1 = kp1;
double kout2 = kp2;
double kout3 = kp3 + evac3;
double kout4 = evac4;

// absorbable water = free water above BOTH the osmotically obliged volume and
// the water bound to the solid phase.  Without the bound-water floor the colon
// can dry a segment to w -> 0, propulsion -> 0, and solids accumulate without
// bound -- in the HEALTHY arm.
double Wosm1 = O1 + O_extra1;  Wosm1 = Wosm1/OSMPL;
double Wosm2 = O2/OSMPL;
double Wosm3 = O3/OSMPL;
double Wosm4 = O4/OSMPL;
double flr1 = (Wosm1 > WBIND*S1) ? Wosm1 : WBIND*S1;
double flr2 = (Wosm2 > WBIND*S2) ? Wosm2 : WBIND*S2;
double flr3 = (Wosm3 > WBIND*S3) ? Wosm3 : WBIND*S3;
double flr4 = (Wosm4 > WBIND*S4) ? Wosm4 : WBIND*S4;
double abs1 = KABS1*ABSX*((W1-flr1 > 0.0) ? (W1-flr1) : 0.0);
double abs2 = KABS2*ABSX*((W2-flr2 > 0.0) ? (W2-flr2) : 0.0);
double abs3 = KABS3*ABSX*((W3-flr3 > 0.0) ? (W3-flr3) : 0.0);
double abs4 = KABS4*ABSX*((W4-flr4 > 0.0) ? (W4-flr4) : 0.0);

dxdt_S1 = SIN            - kout1*S1;
dxdt_S2 = kout1*S1       - kout2*S2;
dxdt_S3 = kout2*S2       - kout3*S3;
dxdt_S4 = kout3*S3       - kout4*S4;

dxdt_W1 = WIN*FLUIDX + VSEC1*sec_gain - abs1 - kout1*W1;
dxdt_W2 = kout1*W1   + VSEC2*sec_gain - abs2 - kout2*W2;
dxdt_W3 = kout2*W2   + VSEC3*sec_gain - abs3 - kout3*W3;
dxdt_W4 = kout3*W3   + VSEC4*sec_gain - abs4 - kout4*W4;

dxdt_O1 = osm_in     - kout1*O1;
dxdt_O2 = kout1*O1   - kout2*O2;
dxdt_O3 = kout2*O2   - kout3*O3;
dxdt_O4 = kout3*O3   - kout4*O4;

dxdt_CUM_SBM  = Rdef;
dxdt_CUM_RESC = Rresc;
double pcomp_r = (hapc_rel > 1.5) ? 1.5 : hapc_rel;
double pcomp   = HILLF(w4, W50, GW) * pcomp_r;
dxdt_CUM_CSBM = Rdef * ((pcomp > 1.0) ? 1.0 : pcomp);

// ============================= symptoms ==================================
double bsfs_t = 1.0 + 6.0*HILLF(w4, WB50, HB);
dxdt_BSFS = KBSFS*(bsfs_t - BSFS);
double strain_t = 4.0*(1.0 - HILLF(w4, W50S, HS));
dxdt_STRAIN = KSTR*(((strain_t > 4.0) ? 4.0 : strain_t) - STRAIN);
double gas = FGAS*(O1 + O2 + O3 + O4);
double dist_t = 4.0*HILLF(Stot + Wtot + 3.0*gas, DIST0, 2.0);
dxdt_DIST = KDIST*(dist_t - DIST);
double sbm_wk = 168.0*Rdef;
double freq_sym = 4.0*(1.0 - HILLF(sbm_wk, 3.5, 2.0));
double sym_t = 0.34*freq_sym + 0.28*STRAIN + 0.22*DIST
               + 0.16*(4.0*(1.0 - HILLF(BSFS, 3.0, 4.0)));
dxdt_PACSYM = KSYM*(sym_t - PACSYM);
dxdt_QOL    = KQOL*(0.85*PACSYM - QOL);

// ================================ CNS ====================================
double pain_t = PAIN0 - PAINMAX*SIGc + ATOLA*TOLA + APAINGI*DIST*0.5;
if (pain_t < 0.0)  pain_t = 0.0;
if (pain_t > 10.0) pain_t = 10.0;
dxdt_PAIN = KPAIN*(pain_t - PAIN);
dxdt_TOLA = KTOLA*OCCc_ag - KTOLAR*TOLA;
double wd_t = WDMAX*DEPEND*HILLF(OCCc_ant, WD50, HWD);
dxdt_WD = KWD*(wd_t - WD);
double naus_t = ENAUSLUB*HILLF(ALUB_LUM, EC50LUB*2.0, 1.0) + ENAUSDIS*DIST*0.5;
dxdt_NAUSEA = KNAUS*(((naus_t > 3.0) ? 3.0 : naus_t) - NAUSEA);

// =============================== safety ==================================
double imp_f = 1.0 - HILLF(sbm_wk, 1.6, 3.0);
double imp_w = 1.0 - HILLF(w4, 0.62, 8.0);
dxdt_IMPACT = KIMP * ((imp_f > 0.0) ? imp_f : 0.0) * ((imp_w > 0.0) ? imp_w : 0.0);

$TABLE
double Cop_nM_o = (AOP_CEN/V1OP)/MWOP*1e6;
double Aop_u_o  = FUOP*Cop_nM_o;
double Cpam_o   = (APAM_CEN/V1ANT)/MWANT*1e6;
double Bp_u_o   = FUANT*Cpam_o;
double lum_op_o = (AOP_LUM/VLUM)/MWOP*1e6*FLUMACC;
double lum_pm_o = (APAM_LUM/VLUM)/MWANT*1e6*FLUMACC*FLUMACCX;
double rga = (Aop_u_o + lum_op_o)/KIOP;
double rgb = (Bp_u_o + lum_pm_o)/KIGUT;
double dg  = 1.0 + rga + rgb;
double rca = COP_BR/KIOP;
double rcb = CPAM_BR/KIANT;
double dc  = 1.0 + rca + rcb;

capture CP_OP    = AOP_CEN/V1OP;
capture CP_PAM   = APAM_CEN/V1ANT;
capture OCCG_AG  = rga/dg;
capture OCCG_ANT = rgb/dg;
capture OCCC_AG  = rca/dc;
capture OCCC_ANT = rcb/dc;
// the selectivity index: enteric blockade bought per unit of CNS blockade
capture SI       = (rcb/dc > 1e-12) ? ((rgb/dg)/(rcb/dc)) : 0.0;
capture W4FRAC   = WFRAC(W4, S4);
capture STOT     = S1+S2+S3+S4;
capture CTT_H    = (S1+S2+S3+S4)/SIN;
'

mod <- mcode("oic", oic_code)

## =============================================================================
##  DRUG LIBRARY — each entry is a parameter overlay, not a separate model.
##  KPUU is the only thing that separates a PAMORA from naloxone.
## =============================================================================
ANTAGONISTS <- list(
  none = list(MWANT=400,  FANT=0,    KAANT=1,   CLANT=20,   V1ANT=80,  V2ANT=80,
              QANT=10, FUANT=1,     KPUU=0,     KIANT=1e9,  KIGUT=1e9,
              FLUMANT=0,    ORALANT=1, FLUMACCX=1.0),
  naloxegol = list(MWANT=651.8, FANT=0.50, KAANT=1.0, CLANT=22.0, V1ANT=68,
              V2ANT=100, QANT=15, FUANT=0.96, KPUU=0.020, KIANT=7.4,
              KIGUT=2.6851, FLUMANT=0.35, ORALANT=1, FLUMACCX=1.0),
  naldemedine = list(MWANT=570.6, FANT=0.55, KAANT=2.5, CLANT=8.5, V1ANT=155,
              V2ANT=90, QANT=12, FUANT=0.065, KPUU=0.012, KIANT=0.34,
              KIGUT=0.00821, FLUMANT=0.30, ORALANT=1, FLUMACCX=1.0),
  methylnaltrexone_sc = list(MWANT=356.5, FANT=0.82, KAANT=3.0, CLANT=17.6,
              V1ANT=84, V2ANT=60, QANT=20, FUANT=0.885, KPUU=0.005, KIANT=28.0,
              KIGUT=0.7945, FLUMANT=0.02, ORALANT=0, FLUMACCX=0.0015),
  methylnaltrexone_po = list(MWANT=356.5, FANT=0.012, KAANT=0.9, CLANT=17.6,
              V1ANT=84, V2ANT=60, QANT=20, FUANT=0.885, KPUU=0.005, KIANT=28.0,
              KIGUT=0.7945, FLUMANT=0.85, ORALANT=1, FLUMACCX=0.0015),
  naloxone_po = list(MWANT=327.4, FANT=0.02, KAANT=1.5, CLANT=90.0, V1ANT=180,
              V2ANT=120, QANT=30, FUANT=0.55, KPUU=1.00, KIANT=1.1,
              KIGUT=1.1, FLUMANT=0.75, ORALANT=1, FLUMACCX=1.0)
)

## dosing helper ---------------------------------------------------------------
oic_events <- function(days = 84,
                       op_dose = 30, op_int = 12,      # oxycodone 60 mg/day
                       pam_dose = 0, pam_int = 24,
                       pro_dose = 0,
                       peg_g = 0, lac_g = 0,
                       lub_ug = 0, lub_int = 12,
                       lin_ug = 0) {
  tmax <- days * 24
  e <- NULL
  add <- function(e, amt, cmt, ii) {
    if (amt <= 0) return(e)
    ev <- ev(amt = amt, cmt = cmt, ii = ii, addl = floor(tmax/ii) - 1, time = 0)
    if (is.null(e)) ev else e + ev
  }
  e <- add(e, op_dose,  "AOP_DEP",  op_int)
  e <- add(e, pam_dose, "APAM_DEP", pam_int)
  e <- add(e, pro_dose * 0.90, "APRO_CEN", 24)          # F ~0.9, fast absorption
  e <- add(e, peg_g / 3.350,   "APEG_LUM", 24)          # g -> mmol (MW 3350)
  e <- add(e, lac_g / 0.3423,  "ALAC_LUM", 24)          # g -> mmol (MW 342.3)
  e <- add(e, lub_ug, "ALUB_LUM", lub_int)
  e <- add(e, lin_ug, "ALIN_LUM", 24)
  e
}

run_oic <- function(drug = "none", days = 84, params = list(),
                    init0 = NULL, ...) {
  p <- ANTAGONISTS[[drug]]
  p[names(params)] <- params
  m <- mod %>% param(p)
  if (!is.null(init0)) m <- m %>% init(init0)   # carry state across segments
  m %>%
    ev(oic_events(days = days, ...)) %>%
    mrgsim(end = days * 24, delta = 1, atol = 1e-9, rtol = 1e-6) %>%
    as_tibble()
}

## final state of a run, as a named list suitable for init0
final_state <- function(out) {
  cmts <- names(mod@init)          # compartment names, in model order
  as.list(tail(out[, cmts], 1))
}

## weekly readout -------------------------------------------------------------
summarise_oic <- function(out, window_h = 168) {
  tail_out <- out %>% filter(time >= max(time) - window_h)
  dt <- max(tail_out$time) - min(tail_out$time)
  last <- tail(out, 1)
  tibble(
    SBM_wk  = (max(tail_out$CUM_SBM)  - min(tail_out$CUM_SBM))  / dt * 168,
    CSBM_wk = (max(tail_out$CUM_CSBM) - min(tail_out$CUM_CSBM)) / dt * 168,
    RESC_wk = (max(tail_out$CUM_RESC) - min(tail_out$CUM_RESC)) / dt * 168,
    BSFS = last$BSFS, STRAIN = last$STRAIN, DIST = last$DIST,
    PACSYM = last$PACSYM, QOL = last$QOL, PAIN = last$PAIN, COWS = last$WD,
    OCCG_AG = last$OCCG_AG, OCCG_ANT = last$OCCG_ANT,
    OCCC_ANT = last$OCCC_ANT, SI = last$SI,
    W4 = last$W4FRAC, CTT_H = last$CTT_H, IMPACT = last$IMPACT
  )
}

## =============================================================================
##  16 TREATMENT SCENARIOS
##  (Values reproduced by the Python reference implementation are quoted in the
##   comment column; see oic_reference_output.txt, analysis I.)
## =============================================================================
SCENARIOS <- list(
  # id, label, drug, argument list                                   # SBM/wk
  list("S01","healthy (no opioid)",              "none", list(op_dose=0)),        # 11.16
  list("S02","OIC untreated (oxycodone 60/d)",   "none", list()),                 #  1.49
  list("S03","+ PEG 3350 17 g/d",                "none", list(peg_g=17)),         #  1.49
  list("S04","+ lactulose 20 g/d",               "none", list(lac_g=20)),         #  6.73
  list("S05","+ lubiprostone 24 ug bid",         "none", list(lub_ug=24)),        #  1.84
  list("S06","+ linaclotide 145 ug qd",          "none", list(lin_ug=145)),       #  1.81
  list("S07","+ prucalopride 2 mg qd",           "none", list(pro_dose=2)),       #  1.92
  list("S08","+ naloxegol 12.5 mg qd",     "naloxegol", list(pam_dose=12.5)),     #  3.11
  list("S09","+ naloxegol 25 mg qd",       "naloxegol", list(pam_dose=25)),       #  4.30
  list("S10","+ naldemedine 0.2 mg qd",  "naldemedine", list(pam_dose=0.2)),      #  5.00
  list("S11","+ methylnaltrexone 12 mg SC qod",
       "methylnaltrexone_sc", list(pam_dose=12, pam_int=48)),                     #  4.22
  list("S12","+ methylnaltrexone 450 mg PO qd",
       "methylnaltrexone_po", list(pam_dose=450)),                                #  6.23
  list("S13","+ naloxone 20 mg PO tid",  "naloxone_po", list(pam_dose=20, pam_int=8)), # 9.75
  list("S14","+ naloxegol 25 + PEG 17 g",  "naloxegol", list(pam_dose=25, peg_g=17)),  # 4.31
  list("S15","naloxegol 25 + strong P-gp inhibitor", "naloxegol",
       list(pam_dose=25)),                                                        #  7.10
  list("S16","methadone 60/d + lubiprostone 24 ug bid", "none", list(lub_ug=24))  #  1.59
)

run_all_scenarios <- function(days = 84) {
  purrr::map_dfr(SCENARIOS, function(s) {
    args <- s[[4]]
    pars <- list()
    if (s[[1]] == "S15") pars <- list(PGPINH = 10.0, CYP3A4I = 3.4)
    if (s[[1]] == "S16") pars <- list(METHADON = 1.0)
    out <- do.call(run_oic, c(list(drug = s[[3]], days = days, params = pars), args))
    summarise_oic(out) %>% mutate(id = s[[1]], label = s[[2]], .before = 1)
  })
}

## =============================================================================
##  ANALYSIS 1 — the selectivity index, computed rather than asserted
## =============================================================================
selectivity_table <- function(days = 84) {
  specs <- list(c("naloxone_po", 20, 8), c("naloxegol", 25, 24),
                c("naldemedine", 0.2, 24), c("methylnaltrexone_sc", 12, 24))
  purrr::map_dfr(specs, function(sp) {
    out <- run_oic(sp[1], days = days,
                   pam_dose = as.numeric(sp[2]), pam_int = as.numeric(sp[3]))
    summarise_oic(out) %>%
      mutate(drug = sp[1], KPUU = ANTAGONISTS[[sp[1]]]$KPUU,
             inv_KPUU = 1 / ANTAGONISTS[[sp[1]]]$KPUU, .before = 1)
  })
}

## =============================================================================
##  ANALYSIS 2 — a P-gp inhibitor rotates the ratio; a CYP3A4 inhibitor does not
##  This is the analysis that says the label's "reduce the dose" advice fixes
##  the wrong variable when the perpetrator drug also inhibits P-gp.
## =============================================================================
ddi_table <- function(days = 84) {
  cases <- list(
    list("no interaction",                    list(), 25),
    list("CYP3A4 inhibitor only (AUC x3.4)",  list(CYP3A4I = 3.4), 25),
    list("CYP3A4 inhibitor, dose 12.5",       list(CYP3A4I = 3.4), 12.5),
    list("P-gp inhibitor only (Kp_uu x10)",   list(PGPINH = 10), 25),
    list("P-gp inhibitor, dose 12.5",         list(PGPINH = 10), 12.5),
    list("dual (verapamil-like)",             list(CYP3A4I = 3.4, PGPINH = 10), 25),
    list("dual, dose reduced to 12.5",        list(CYP3A4I = 3.4, PGPINH = 10), 12.5)
  )
  purrr::map_dfr(cases, function(cs) {
    out <- run_oic("naloxegol", days = days, params = cs[[2]], pam_dose = cs[[3]])
    summarise_oic(out) %>% mutate(case = cs[[1]], dose = cs[[3]], .before = 1)
  })
}

## =============================================================================
##  ANALYSIS 3 — which brake carries the drug response?
##  Disable each enteric MOR brake in turn and re-measure the SAME 0 -> 25 mg
##  naloxegol comparison.  In the reference implementation this REFUTED the
##  hypothesis it was written to test: removing the transit-hydration feedback
##  loop barely changes the fold-response.  The amplification is the product of
##  several independent brakes released together, not a feedback gain.
## =============================================================================
brake_decomposition <- function(days = 84) {
  variants <- list(
    list("(none -- full model)",                    list()),
    list("hydration -> propulsion loop (GW=GWD=0)", list(GW = 0, GWD = 0)),
    list("rectal urge threshold shift (AVSENS=0)",  list(AVSENS = 0)),
    list("anal sphincter tone (AANO=0)",            list(AANO = 0)),
    list("segmental muscle tone (ATONE=0)",         list(ATONE = 0)),
    list("inhibitory NO/VIP tone (ANO=0)",          list(ANO = 0)),
    list("secretion suppression (ESECMOR=0)",       list(ESECMOR = 0)),
    list("cAMP -> ACh -> HAPC drive (EMAXMOR=0)",   list(EMAXMOR = 0))
  )
  purrr::map_dfr(variants, function(v) {
    u <- summarise_oic(run_oic("none", days = days, params = v[[2]]))$SBM_wk
    a <- summarise_oic(run_oic("naloxegol", days = days, params = v[[2]],
                               pam_dose = 25))$SBM_wk
    tibble(brake_disabled = v[[1]], untreated = u, naloxegol25 = a,
           fold = a / pmax(u, 1e-9))
  })
}

## =============================================================================
##  ANALYSIS 4 — tolerance asymmetry under a pain-holding dose controller
##  Equilibrate 12 weeks, then titrate weekly to hold NRS 4.0 for 24 weeks.
##  The reference run shows CNS receptor availability collapsing while enteric
##  availability holds — and the gut ENDPOINT not deteriorating, because enteric
##  transduction is already saturated at an ordinary analgesic dose.
## =============================================================================
tolerance_titration <- function(weeks = 24, target = 4.0, equilibrate_days = 84) {
  ## Equilibrate FIRST, otherwise week 0 reports a transient rather than the
  ## patient's actual starting state, and then carry the state forward across
  ## every weekly segment -- the whole analysis is about a slow receptor
  ## trafficking process, so restarting from the initial conditions each week
  ## would erase exactly the thing being measured.
  dose <- 30
  st <- final_state(run_oic("none", days = equilibrate_days, op_dose = dose))
  rows <- list()
  for (wk in seq_len(weeks)) {
    out <- run_oic("none", days = 7, op_dose = dose, init0 = st)
    st  <- final_state(out)
    s   <- summarise_oic(out, window_h = 168)
    rows[[wk]] <- s %>% mutate(week = wk, oxy_mg_day = dose * 2, .before = 1)
    dose <- if (s$PAIN > target + 0.25) dose * 1.22
            else if (s$PAIN < target - 0.25) dose * 0.94 else dose
    dose <- min(dose, 100)          # cap at 200 mg/day, a real-world ceiling
  }
  bind_rows(rows)
}

## =============================================================================
##  ANALYSIS 5 — osmole arithmetic: PEG 3350 vs lactulose
##  17 g of a 3350 Da polymer is 5.1 mmol and obliges ~17 mL of water.  20 g of
##  a 342 Da disaccharide is 58 mmol AND is cleaved by colonic flora into ~3.6x
##  more osmoles.  The model therefore predicts PEG 3350 is close to inert in
##  OIC — a prediction the author believes is WRONG, and which is discussed as
##  an open discrepancy in README.md rather than tuned away.
## =============================================================================
osmole_arithmetic <- function() {
  tibble(
    agent = c("PEG 3350 17 g", "lactulose 20 g", "lactulose 40 g"),
    g     = c(17, 20, 40),
    mw    = c(3350, 342.3, 342.3),
    ferm  = c(1, 3.6, 3.6)
  ) %>%
    mutate(mmol = g / mw * 1000,
           obliged_mL = mmol / 0.290,
           after_fermentation_mL = obliged_mL * ferm)
}

## =============================================================================
##  PLOTS
## =============================================================================
plot_scenarios <- function(res) {
  res %>%
    mutate(label = factor(label, levels = rev(label))) %>%
    ggplot(aes(SBM_wk, label)) +
    geom_col(fill = "#4f88a8") +
    geom_vline(xintercept = 3, linetype = "dashed", colour = "#c0392b") +
    labs(x = "SBM / week (dashed line = responder threshold)", y = NULL,
         title = "Opioid-induced constipation: 16 scenarios",
         subtitle = "oxycodone 60 mg/day background unless stated") +
    theme_minimal(base_size = 11)
}

plot_selectivity <- function(sel) {
  sel %>%
    ggplot(aes(OCCC_ANT, OCCG_ANT, label = drug)) +
    geom_point(size = 4, colour = "#2e8b57") +
    geom_text(nudge_y = 0.04, size = 3.4) +
    scale_x_log10() +
    labs(x = "CNS antagonist occupancy (log)", y = "enteric antagonist occupancy",
         title = "One receptor, two compartments",
         subtitle = "Distance from the diagonal is the therapeutic window") +
    theme_minimal(base_size = 11)
}

## =============================================================================
##  MAIN
## =============================================================================
if (sys.nframe() == 0) {
  message("Running 16 OIC scenarios ...")
  res <- run_all_scenarios()
  print(res, n = 20)
  message("\nSelectivity index:")
  print(selectivity_table())
  message("\nDrug-drug interaction (P-gp vs CYP3A4):")
  print(ddi_table())
  message("\nBrake decomposition:")
  print(brake_decomposition())
  message("\nOsmole arithmetic:")
  print(osmole_arithmetic())
  print(plot_scenarios(res))
}
