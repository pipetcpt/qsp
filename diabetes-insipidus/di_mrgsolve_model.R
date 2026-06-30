## =============================================================================
##  Diabetes Insipidus (Arginine Vasopressin Disorders, AVP-D / AVP-R)
##  Quantitative Systems Pharmacology Model — mrgsolve specification
## -----------------------------------------------------------------------------
##  Covers:
##    * Endogenous AVP synthesis / release driven by plasma osmolality
##      and effective circulating volume (Robertson 1976; Bichet 2019).
##    * Renal V2-receptor -> Gs -> cAMP -> PKA -> AQP2 trafficking
##      (Nielsen 1995; Fenton 2007); short-term shuttling + long-term TX.
##    * Free-water balance (TBW, plasma Na+, plasma osmolality) and
##      thirst-driven oral intake (Christ-Crain 2019, 2020 NEJM).
##    * Pharmacology: desmopressin (SC, intranasal, oral, sublingual MELT),
##      thiazide (HCTZ paradoxical antidiuresis in AVP-R), amiloride
##      (lithium-induced NDI), indomethacin (PGE2/AQP2), tolvaptan
##      (V2 antagonist comparator for SIADH), lithium (NDI inducer).
##    * Disease phenotypes via parameter switches:
##         AVP_D_sev    1 = complete central DI (AVP-D), 0.5 = partial
##         AVP_R_sev    1 = complete nephrogenic DI (AVP-R)
##         PSYCH_drive  L/d added to baseline drinking (primary polydipsia)
##         GEST         1 enables placental vasopressinase x4
##  Calibration anchors:
##    * AVP plasma osm slope ~0.4 pmol/L per mOsm/kg (Robertson 1976).
##    * DDAVP SC 1 ug -> Uosm 600-800 mOsm/kg sustained 6-8 h (Vavra 1968).
##    * DDAVP IN 10 ug typical adult dose; PO 100-400 ug TID.
##    * Lithium 30 mmol/d steady-state -> ~30% urinary concentrating defect
##      after >10 y (Walker 1986; Bedford 2008 JASN).
##    * Christ-Crain 2019 NEJM hypertonic saline copeptin cut-off 4.9 pmol/L.
##    * 24-h urine output in untreated severe CDI 8-15 L/d; full DDAVP
##      replacement returns to ~1.5-2.5 L/d (Verbalis 2003; Arima 2014).
##  Units:  time hours; volumes L; mass µg (drugs) and pmol (peptides);
##          osmolality mOsm/kg; Na+ mmol/L.
## =============================================================================

library(mrgsolve)

di_model <- '
$PROB
# Diabetes Insipidus (AVP-D / AVP-R) QSP model
#  20-ODE system: AVP secretion, V2R/cAMP/AQP2 signalling, water balance,
#  desmopressin (SC/IN/PO/SL) PK, plus thiazide/amiloride/indomethacin/
#  tolvaptan/lithium handles.

$PARAM @annotated
// ===== Demographics =====
WT     :  70   : Body weight (kg)
SEX    :  1    : 1=male, 0=female  (affects TBW fraction)
AGE    :  40   : Age (yr)
eGFR   :  95   : eGFR (mL/min/1.73 m2)

// ===== Disease switches =====
AVP_D_sev : 1.0 : AVP-D (central DI) severity 0-1, 1=complete
AVP_R_sev : 0.0 : AVP-R (nephrogenic DI) severity 0-1
PSYCH_drive: 0.0 : Extra drinking (L/d) for primary polydipsia
GEST   : 0    : 1 = pregnancy with placental vasopressinase

// ===== Plasma osmolality / sodium (initial) =====
P_OSM0 : 287  : Initial plasma osmolality (mOsm/kg)
P_NA0  : 140  : Initial plasma Na+ (mmol/L)
OSM_TH_AVP   : 280 : AVP osmotic threshold (mOsm/kg)
OSM_TH_THIRST: 290 : Thirst threshold (mOsm/kg)
SLOPE_AVP    : 0.4 : Plasma AVP slope (pmol/L per mOsm/kg)
AVP_MAX_OSM  : 20  : Max osmotically driven AVP (pmol/L)

// ===== Endogenous AVP kinetics =====
KAVP_SYN : 0.5  : Maximal AVP secretion (pmol/L per h) (drive units)
KAVP_DEG : 3.0  : AVP elimination rate (per h)  -> t1/2 ~14 min
AVP_BASAL_DRIVE : 1.0 : Constitutive low-level drive
VPase_GEST : 4.0 : Placental vasopressinase activity (relative)

// ===== Thirst drive (drinking) =====
DRINK_MAX  : 18  : Maximum voluntary intake (L/d)
DRINK_BASE : 1.8 : Basal water intake (L/d, ad lib)
KTH_SLOPE  : 0.6 : Drink (L/d) per (mOsm/kg above thirst threshold)

// ===== V2R / AQP2 signalling =====
EC50_V2  : 1.5  : EC50 of total V2 agonist (pmol/L AVP-equiv)
HILL_V2  : 1.4  : Hill slope V2R
AQP2_BASE: 0.10 : Constitutive AQP2 apical fraction (AVP-R "leak")
KAQP2_ON : 0.8  : Rate constant AQP2 vesicle insertion (per h)
KAQP2_OFF: 0.5  : Rate constant AQP2 internalization (per h)
KAQP2_TX : 0.04 : AQP2 transcription rate (per h, long term)
KAQP2_DEG: 0.06 : AQP2 mRNA degradation (per h)
AQP2_MAX : 1.0  : Maximum apical AQP2 (normalized)

// ===== Renal water handling =====
GFR_LH   : 5.4  : Glomerular filtrate (L/h) ~ 130 mL/min
PCT_REABS: 0.67 : Proximal tubule fractional H2O reabs
LOOP_REABS: 0.15: Loop of Henle fractional reabs
DCT_REABS: 0.05 : DCT reabs (thiazide-sensitive Na-driven)
UOSM_MAX : 1200 : Maximum urinary osmolality (mOsm/kg) full AQP2
UOSM_MIN : 50   : Minimum urinary osmolality (mOsm/kg) zero AQP2
SOLUTE_LOAD : 900 : Daily solute load (mOsm/d)

// ===== Sodium dynamics =====
NA_INTAKE : 150 : Dietary Na+ intake (mmol/d)
KNA_EX    : 0.05: Na+ exchange rate (per h)

// ===== Drug PK -- Desmopressin (multi-route) =====
KA_DDAVP_SC : 1.5   : SC absorption rate (per h)
KA_DDAVP_IN : 1.2   : Intranasal absorption (per h)
KA_DDAVP_PO : 1.0   : Oral absorption rate (per h)
KA_DDAVP_SL : 1.4   : Sublingual lyophilisate absorption (per h)
F_DDAVP_SC  : 1.00  : SC bioavailability
F_DDAVP_IN  : 0.04  : Intranasal F (3-5%)
F_DDAVP_PO  : 0.0016: Oral F (~0.16%)
F_DDAVP_SL  : 0.0025: Sublingual F (~0.25%)
CL_DDAVP    : 7.6   : DDAVP clearance (L/h, 70 kg)
V1_DDAVP    : 12.0  : Central V (L)
V2_DDAVP    : 25.0  : Peripheral V (L)
Q_DDAVP     : 2.5   : Inter-compartmental CL (L/h)
DDAVP_POT   : 12    : V2 potency (DDAVP : AVP molar equiv)

// ===== Drug PK -- Thiazide (HCTZ) =====
KA_HCTZ : 1.3 : HCTZ absorption (per h)
CL_HCTZ : 22  : HCTZ clearance (L/h)
V_HCTZ  : 60  : Volume of distribution (L)
EC50_HCTZ : 50 : HCTZ plasma EC50 (ng/mL) for NCC inhibition
IMAX_HCTZ : 0.75 : Maximal NCC inhibition (paradoxical antidiuresis)

// ===== Amiloride =====
KA_AMI : 1.4 : Amiloride absorption (per h)
CL_AMI : 10  : Clearance (L/h)
V_AMI  : 300 : Volume (L)
EC50_AMI : 20 : Plasma EC50 (ng/mL) for ENaC/Li+ entry block
IMAX_AMI : 0.6 : Max ENaC inhibition

// ===== Indomethacin =====
KA_IND : 2.0 : Indomethacin absorption (per h)
CL_IND : 8   : Clearance (L/h)
V_IND  : 18  : Volume (L)
EC50_IND : 1.5 : EC50 (mg/L) for COX/PGE2 inhibition
IMAX_IND : 0.5 : Max PGE2 suppression -> AQP2 traffic boost

// ===== Tolvaptan (V2 antagonist comparator) =====
KA_TOL : 1.0 : Tolvaptan absorption (per h)
CL_TOL : 9.7 : Clearance (L/h)
V_TOL  : 230 : Volume (L)
KI_TOL : 0.7 : V2R inhibition constant (ng/mL)

// ===== Lithium (NDI inducer) =====
KA_LI : 1.0 : Li+ absorption (per h)
CL_LI : 1.5 : Lithium clearance (L/h)  ~ 0.025 L/min
V_LI  : 50  : Lithium Vd (L)
KLI_NDI_ON  : 1e-4  : Rate of NDI induction (per h per mmol/L*time)
KLI_NDI_OFF : 5e-4  : Reversal once Li+ stopped

// ===== Safety / hyponatremia hazard =====
KHYPONA : 0.001 : Per-hour hazard scaling when Na+ <130

$CMT @annotated
DEPOT_SC : SC desmopressin depot (ug)
CENT_DDAVP : Central desmopressin (ug)
PERI_DDAVP : Peripheral desmopressin (ug)
DEPOT_IN : Intranasal depot (ug)
DEPOT_PO : Oral depot (ug)
DEPOT_SL : Sublingual depot (ug)
HCTZ_DEPOT : HCTZ gut (mg)
HCTZ_CENT  : HCTZ plasma (mg)
AMI_DEPOT  : Amiloride gut (mg)
AMI_CENT   : Amiloride plasma (mg)
IND_DEPOT  : Indomethacin gut (mg)
IND_CENT   : Indomethacin plasma (mg)
TOL_DEPOT  : Tolvaptan gut (mg)
TOL_CENT   : Tolvaptan plasma (mg)
LI_DEPOT   : Lithium gut (mmol)
LI_CENT    : Lithium plasma (mmol)
AVP_E      : Endogenous plasma AVP (pmol/L)
AQP2_M     : AQP2 mRNA (normalized)
AQP2_A     : AQP2 apical (normalized 0-1)
TBW        : Total body water (L)
NA_BODY    : Total exchangeable Na+ (mmol)
NDI_LI     : Lithium-induced NDI severity (0-1)
CUM_URINE  : Cumulative urine output (L)
CUM_HAZ    : Cumulative hyponatremia hazard

$MAIN
// Initial conditions
double TBW0 = WT * (SEX ? 0.60 : 0.50);
TBW_0      = TBW0;
NA_BODY_0  = TBW0 * P_NA0;
AVP_E_0    = 1.2;
AQP2_M_0   = 1.0 * (1.0 - 0.7 * AVP_R_sev);
AQP2_A_0   = 0.5 * (1.0 - 0.85 * AVP_R_sev);
NDI_LI_0   = 0.0;

// Physiologic eGFR scaling
double GFRsc = eGFR / 95.0;

$ODE
// ---------- Plasma osmolality & sodium ----------
double Posm   = NA_BODY / TBW * 2.0;        // 2 * [Na+] approximation
double P_Na   = NA_BODY / TBW;

// ---------- Endogenous AVP ----------
double osmStim = (Posm > OSM_TH_AVP) ? (Posm - OSM_TH_AVP) * SLOPE_AVP : 0.0;
osmStim       = (osmStim > AVP_MAX_OSM) ? AVP_MAX_OSM : osmStim;
double avpProd_cap = (1.0 - AVP_D_sev) * (AVP_BASAL_DRIVE + osmStim);
double VPase = GEST ? VPase_GEST : 1.0;
dxdt_AVP_E   = KAVP_SYN * avpProd_cap - KAVP_DEG * VPase * AVP_E;

// ---------- Desmopressin PK ----------
dxdt_DEPOT_SC = -KA_DDAVP_SC * DEPOT_SC;
dxdt_DEPOT_IN = -KA_DDAVP_IN * DEPOT_IN;
dxdt_DEPOT_PO = -KA_DDAVP_PO * DEPOT_PO;
dxdt_DEPOT_SL = -KA_DDAVP_SL * DEPOT_SL;

double inSC = F_DDAVP_SC * KA_DDAVP_SC * DEPOT_SC;
double inIN = F_DDAVP_IN * KA_DDAVP_IN * DEPOT_IN;
double inPO = F_DDAVP_PO * KA_DDAVP_PO * DEPOT_PO;
double inSL = F_DDAVP_SL * KA_DDAVP_SL * DEPOT_SL;

double k10 = CL_DDAVP / V1_DDAVP;
double k12 = Q_DDAVP / V1_DDAVP;
double k21 = Q_DDAVP / V2_DDAVP;

dxdt_CENT_DDAVP = inSC + inIN + inPO + inSL
                  - (k10 + k12) * CENT_DDAVP
                  + k21 * PERI_DDAVP;
dxdt_PERI_DDAVP = k12 * CENT_DDAVP - k21 * PERI_DDAVP;

double Cddavp_pgmL = (CENT_DDAVP / V1_DDAVP) * 1000.0; // ug/L -> ng/L = pg/mL
// AVP-equivalent stimulation (DDAVP ~12x molar potency at V2)
double Cddavp_eq_pmol = (CENT_DDAVP / V1_DDAVP) * 1e6 / 1069.0 * DDAVP_POT;

// ---------- Other drug PKs ----------
dxdt_HCTZ_DEPOT = -KA_HCTZ * HCTZ_DEPOT;
dxdt_HCTZ_CENT  =  KA_HCTZ * HCTZ_DEPOT - (CL_HCTZ / V_HCTZ) * HCTZ_CENT;
double Chctz = HCTZ_CENT / V_HCTZ * 1000.0; // ng/mL

dxdt_AMI_DEPOT = -KA_AMI * AMI_DEPOT;
dxdt_AMI_CENT  =  KA_AMI * AMI_DEPOT - (CL_AMI / V_AMI) * AMI_CENT;
double Cami = AMI_CENT / V_AMI * 1000.0; // ng/mL

dxdt_IND_DEPOT = -KA_IND * IND_DEPOT;
dxdt_IND_CENT  =  KA_IND * IND_DEPOT - (CL_IND / V_IND) * IND_CENT;
double Cind = IND_CENT / V_IND; // mg/L

dxdt_TOL_DEPOT = -KA_TOL * TOL_DEPOT;
dxdt_TOL_CENT  =  KA_TOL * TOL_DEPOT - (CL_TOL / V_TOL) * TOL_CENT;
double Ctol = TOL_CENT / V_TOL * 1000.0; // ng/mL
double tolV2block = Ctol / (Ctol + KI_TOL);

dxdt_LI_DEPOT  = -KA_LI * LI_DEPOT;
dxdt_LI_CENT   =  KA_LI * LI_DEPOT - (CL_LI / V_LI) * LI_CENT;
double Cli = LI_CENT / V_LI; // mmol/L plasma Li+

// ---------- Lithium-induced NDI severity (slow accumulation) ----------
double liEffect = Cli > 0.6 ? (Cli - 0.6) : 0.0;
dxdt_NDI_LI = KLI_NDI_ON * liEffect - KLI_NDI_OFF * NDI_LI;
double NDI_total = (AVP_R_sev > NDI_LI) ? AVP_R_sev : NDI_LI;

// ---------- V2R / cAMP / AQP2 dynamics ----------
double V2_signal = AVP_E + Cddavp_eq_pmol;
double V2_act_raw = pow(V2_signal, HILL_V2) /
                    (pow(EC50_V2, HILL_V2) + pow(V2_signal, HILL_V2));
// AVP-R blunts maximal AQP2 trafficking; tolvaptan blocks V2R
double V2_act = V2_act_raw * (1.0 - NDI_total) * (1.0 - tolV2block);

// Indomethacin boosts trafficking (reduced PGE2 -> retain AQP2)
double indomBoost = 1.0 + IMAX_IND * (Cind / (Cind + EC50_IND));

// AQP2 mRNA (long term)
dxdt_AQP2_M = KAQP2_TX * (V2_act * indomBoost) * (1.0 - NDI_total) - KAQP2_DEG * (AQP2_M - 0.05);

// AQP2 apical pool — short term shuttling + bounded by mRNA
double aqp2_target = AQP2_BASE + (AQP2_MAX - AQP2_BASE) * V2_act * indomBoost;
aqp2_target = aqp2_target * AQP2_M;
dxdt_AQP2_A = KAQP2_ON * (aqp2_target - AQP2_A) - KAQP2_OFF * AQP2_A * (1.0 - V2_act);

// ---------- Urine output & concentration ----------
double Uosm = UOSM_MIN + (UOSM_MAX - UOSM_MIN) * AQP2_A;
double UrineFlow_Ld = SOLUTE_LOAD / Uosm; // L/d
// HCTZ paradoxical antidiuresis in AVP-R: thiazide -> ECV down -> PCT reabs up
double thzEff = IMAX_HCTZ * (Chctz / (Chctz + EC50_HCTZ));
UrineFlow_Ld = UrineFlow_Ld * (1.0 - 0.45 * thzEff * NDI_total);
// Amiloride blocks ENaC: protects against Li+-induced AQP2 down (works on NDI severity)
double amiEff = IMAX_AMI * (Cami / (Cami + EC50_AMI));
// Amiloride directly attenuates NDI_LI signaling (built into UrineFlow via NDI cleanup)
UrineFlow_Ld = UrineFlow_Ld * (1.0 - 0.30 * amiEff * (NDI_LI > 0));
// Convert to L/h, apply eGFR scaling
double UrineFlow = UrineFlow_Ld / 24.0 * GFRsc;

dxdt_CUM_URINE = UrineFlow;

// ---------- Drinking behavior ----------
double thirstDrive = (Posm > OSM_TH_THIRST) ? (Posm - OSM_TH_THIRST) * KTH_SLOPE : 0.0;
double drink_Ld    = DRINK_BASE + thirstDrive + PSYCH_drive;
if (drink_Ld > DRINK_MAX) drink_Ld = DRINK_MAX;
double drink = drink_Ld / 24.0;

// ---------- Sodium balance ----------
double Na_loss = UrineFlow * (P_Na * 0.4 + 30.0); // rough urinary Na (mmol/L) ~ tubular
double Na_in   = NA_INTAKE / 24.0;
dxdt_NA_BODY   = Na_in - Na_loss;

// ---------- Water balance ----------
dxdt_TBW = drink - UrineFlow - (0.5/24.0); // insensible losses ~0.5 L/d

// ---------- Hyponatremia hazard ----------
double hypoTrig = (P_Na < 130) ? (130 - P_Na) : 0.0;
dxdt_CUM_HAZ   = KHYPONA * hypoTrig * hypoTrig;

$TABLE
double DDAVP_pg = (CENT_DDAVP / V1_DDAVP) * 1e6; // pg/mL
double Uosm_out = UOSM_MIN + (UOSM_MAX - UOSM_MIN) * AQP2_A;
double P_Na_out = NA_BODY / TBW;
double P_Osm_out= 2.0 * P_Na_out;
double UrineLday = SOLUTE_LOAD / Uosm_out;
double thirst   = (P_Osm_out > OSM_TH_THIRST) ? (P_Osm_out - OSM_TH_THIRST) : 0.0;
double V2occ    = (AVP_E + (CENT_DDAVP/V1_DDAVP)*1e6/1069.0*DDAVP_POT) /
                  (EC50_V2 + AVP_E + (CENT_DDAVP/V1_DDAVP)*1e6/1069.0*DDAVP_POT);
double copeptin = (1.0 - AVP_D_sev) * (2.0 + AVP_E);
capture DDAVP_conc_pgmL = DDAVP_pg;
capture Uosm_mosm       = Uosm_out;
capture Urine_Lday      = UrineLday;
capture Plasma_Na       = P_Na_out;
capture Plasma_Osm      = P_Osm_out;
capture Thirst_score    = thirst;
capture V2_occupancy    = V2occ;
capture Copeptin_pmol   = copeptin;
capture AVP_pmolL       = AVP_E;
capture AQP2_apical     = AQP2_A;
capture NDI_total_sev   = (AVP_R_sev > NDI_LI) ? AVP_R_sev : NDI_LI;
'

di_mod <- mcode("diabetes_insipidus_qsp", di_model)

## =============================================================================
## Scenario library
## =============================================================================

run_scenario <- function(scenario = c("untreated_CDI",
                                      "DDAVP_SC_2ug_BID",
                                      "DDAVP_IN_10ug_BID",
                                      "DDAVP_PO_200ug_TID",
                                      "DDAVP_SL_120ug_TID",
                                      "NDI_lithium_HCTZ",
                                      "NDI_lithium_amiloride",
                                      "NDI_indomethacin",
                                      "tolvaptan_SIADH_comparator",
                                      "primary_polydipsia",
                                      "gestational_DDAVP",
                                      "pediatric_DDAVP_SC"),
                         duration_h = 240) {

  scenario <- match.arg(scenario)
  mod <- di_mod

  ev <- ev()
  param_update <- list()

  if (scenario == "untreated_CDI") {
    param_update <- list(AVP_D_sev = 1.0, AVP_R_sev = 0.0)
  } else if (scenario == "DDAVP_SC_2ug_BID") {
    param_update <- list(AVP_D_sev = 1.0)
    ev <- ev(amt = 2, cmt = "DEPOT_SC", ii = 12, addl = 19) # ug q12h
  } else if (scenario == "DDAVP_IN_10ug_BID") {
    param_update <- list(AVP_D_sev = 1.0)
    ev <- ev(amt = 10, cmt = "DEPOT_IN", ii = 12, addl = 19)
  } else if (scenario == "DDAVP_PO_200ug_TID") {
    param_update <- list(AVP_D_sev = 1.0)
    ev <- ev(amt = 200, cmt = "DEPOT_PO", ii = 8, addl = 29)
  } else if (scenario == "DDAVP_SL_120ug_TID") {
    param_update <- list(AVP_D_sev = 1.0)
    ev <- ev(amt = 120, cmt = "DEPOT_SL", ii = 8, addl = 29)
  } else if (scenario == "NDI_lithium_HCTZ") {
    param_update <- list(AVP_D_sev = 0.0, AVP_R_sev = 0.8)
    ev1 <- ev(amt = 25, cmt = "HCTZ_DEPOT", ii = 24, addl = 9) # mg q24h
    ev  <- ev1
  } else if (scenario == "NDI_lithium_amiloride") {
    param_update <- list(AVP_D_sev = 0.0, AVP_R_sev = 0.5)
    ev1 <- ev(amt = 30, cmt = "LI_DEPOT", ii = 12, addl = 19) # 30 mmol BID
    ev2 <- ev(amt = 10, cmt = "AMI_DEPOT", ii = 24, addl = 9, time = 24)
    ev  <- ev1 + ev2
  } else if (scenario == "NDI_indomethacin") {
    param_update <- list(AVP_D_sev = 0.0, AVP_R_sev = 0.7)
    ev <- ev(amt = 50, cmt = "IND_DEPOT", ii = 8, addl = 29) # 50 mg q8h
  } else if (scenario == "tolvaptan_SIADH_comparator") {
    param_update <- list(AVP_D_sev = 0.0, AVP_R_sev = 0.0)
    ev <- ev(amt = 15, cmt = "TOL_DEPOT", ii = 24, addl = 9) # mg QD
  } else if (scenario == "primary_polydipsia") {
    param_update <- list(AVP_D_sev = 0.0, AVP_R_sev = 0.0, PSYCH_drive = 6)
  } else if (scenario == "gestational_DDAVP") {
    param_update <- list(AVP_D_sev = 0.6, GEST = 1)
    ev <- ev(amt = 10, cmt = "DEPOT_IN", ii = 12, addl = 19)
  } else if (scenario == "pediatric_DDAVP_SC") {
    param_update <- list(WT = 20, SEX = 1, AVP_D_sev = 1.0)
    ev <- ev(amt = 0.3, cmt = "DEPOT_SC", ii = 12, addl = 19) # 0.3 ug q12h
  }

  if (length(param_update)) mod <- param(mod, .list = param_update)

  out <- mrgsim(mod, ev = ev, end = duration_h, delta = 0.25)
  list(scenario = scenario, sim = out)
}

## =============================================================================
## Quick-look helper
## =============================================================================
if (FALSE) {
  res <- run_scenario("DDAVP_SC_2ug_BID", duration_h = 96)
  plot(res$sim, DDAVP_conc_pgmL + Uosm_mosm + Plasma_Na + Urine_Lday ~ time,
       main = "Desmopressin SC 2 ug q12h in complete AVP-D")
}
