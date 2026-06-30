# =============================================================================
# Achalasia (식도이완불능증) — QSP Model
# Disease: Primary idiopathic achalasia — Chicago Classification v4.0
# Phenotype: Impaired LES relaxation (IRP > 15 mmHg) + esophageal aperistalsis
# Pathobiology: Autoimmune (HSV-1/HLA-DQA1*0103/DQB1*0603) loss of
#               inhibitory NO/VIP myenteric neurons → LES non-relaxation +
#               body aperistalsis → dysphagia, regurgitation, weight loss
# Drugs modeled (PK + PD): ISDN (NO donor), Nifedipine (Cav1.2 blocker),
#                          Sildenafil (PDE5 inhibitor), Tadalafil (PDE5),
#                          Botulinum toxin A (SNAP-25 cleavage, intra-LES)
# Procedures: PD (pneumatic dilation), LHM (Heller myotomy), POEM (peroral
#             endoscopic myotomy) modeled as instantaneous LES_PRESS step changes
# Scenarios (10): no Tx, ISDN, Nifedipine, Sildenafil, Tadalafil, Botox,
#                 PD, LHM, POEM, Combo (ISDN + Sild)
# Endpoints: Eckardt 0-12, IRP, LES residual pressure, esophageal diameter,
#            TBE 5-min column, QoL, AE flags
# Calibration anchors:
#   - Boeckxstaens 2011 NEJM (European Achalasia Trial PD vs LHM 86 vs 90% at 2y)
#   - Werner 2019 NEJM (POEM vs LHM 83% vs 81% at 2y, Hannover)
#   - Ponds 2019 JAMA (POEM vs PD 92% vs 54% at 2y)
#   - Pasricha 1995 NEJM, 1996 GE (Botox pilot/RCT)
#   - Triadafilopoulos 1991 Dig Dis Sci (Nifedipine SL crossover)
#   - Bortolotti 2000 Gastroenterology (Sildenafil 50 mg)
#   - Eckardt 1992 Gastroenterology (score validation)
#   - Yadlapati 2021 Neurogastroenterol Motil (Chicago v4.0)
#   - Vaezi 2020 ACG guideline; Zaninotto 2018 ISDE
# =============================================================================

library(mrgsolve)

code <- '
$PARAM @annotated
// ---------- Patient covariates -----------------------------------------------
BWT       :  70   : Body weight (kg)
AGE       :  50   : Age (years)
SEX       :  0    : Sex (0=F, 1=M)
SUBTYPE   :  2    : Chicago subtype (1=I, 2=II, 3=III)
DURATION  :  3    : Disease duration (years)
SURG_OK   :  1    : Surgical candidate (1=yes, 0=no)
CYP3A4    :  1    : CYP3A4 phenotype (0.5 PM, 1 EM, 2 UM)

// ---------- Isosorbide dinitrate (ISDN) — sublingual / PO --------------------
ISDN_F        : 1.00  : SL bioavail (F=1 SL; F=0.22 PO)
ISDN_KA       : 1.40  : SL absorption rate (1/h)
ISDN_CL       : 90    : Clearance L/h (high first-pass mononitrate conversion)
ISDN_V        : 90    : Volume of distribution (L) — SL
ISDN_ISMN_CL  : 10    : ISMN clearance (active mononitrate proxy, L/h)
ISDN_EMAX     : 0.40  : Max ΔLES tone reduction (fraction)
ISDN_EC50     : 5.0   : EC50 in plasma (mg/L) effective in LES
ISDN_HEAD_E   : 0.50  : Headache potency (per mg/L)

// ---------- Nifedipine -------------------------------------------------------
NIF_F      : 0.50  : SL/PO bioavailability
NIF_KA     : 1.20  : Absorption rate (1/h)
NIF_CL     : 30    : Clearance L/h
NIF_V      : 100   : Volume of distribution (L)
NIF_EMAX   : 0.35  : Max ΔLES tone reduction (fraction)
NIF_EC50   : 30    : EC50 (ng/mL)
NIF_HYPO_E : 0.40  : Postural hypotension potency

// ---------- Sildenafil -------------------------------------------------------
SIL_F      : 0.41  : Oral bioavailability
SIL_KA     : 1.50  : Absorption rate (1/h)
SIL_CL     : 41    : Clearance L/h
SIL_V      : 105   : Volume of distribution (L)
SIL_EMAX   : 0.40  : Max ΔLES tone reduction
SIL_EC50   : 100   : EC50 (ng/mL)
SIL_HEAD_E : 0.30  : Headache potency
SIL_HYPO_E : 0.20  : Hypotension potency

// ---------- Tadalafil --------------------------------------------------------
TAD_F      : 0.36  : Oral bioavailability
TAD_KA     : 0.80  : Absorption rate (1/h)
TAD_CL     : 2.5   : Clearance L/h (slow → 17 h t½)
TAD_V      : 60    : Volume of distribution (L)
TAD_EMAX   : 0.40  : Max ΔLES tone reduction
TAD_EC50   : 80    : EC50 (ng/mL)

// ---------- Botulinum toxin A (Botox) intra-LES injection --------------------
BOT_DOSE   : 100   : Total Botox dose (U) per session
BOT_KE     : 0.0040 : Effective LES decay rate (1/d) → ~6 mo duration
BOT_EMAX   : 0.55  : Maximal ΔLES tone reduction (fraction)
BOT_EC50   : 25    : Local U equivalent giving half-max (units)

// ---------- Pneumatic dilation / LHM / POEM (procedural steps) ---------------
PROC_PD    : 0     : Pneumatic dilation flag (one-shot at t=0+)
PROC_LHM   : 0     : Laparoscopic Heller myotomy flag
PROC_POEM  : 0     : POEM flag
POEM_SUB3  : 1     : POEM advantage in Type III (longer myotomy)

// Procedural ΔLES (fraction reduction, immediate but slow re-creep)
PD_DELTA   : 0.45  : LES pressure reduction after PD
LHM_DELTA  : 0.65  : LES pressure reduction after LHM
POEM_DELTA : 0.70  : LES pressure reduction after POEM (Type I/II)
POEM_T3_BONUS : 0.10 : Extra ΔLES with long POEM in Type III

// LES recurrence (creep back) rates (per year)
RECREEP_PD   : 0.06 : Per-year LES recovery fraction after PD (high)
RECREEP_LHM  : 0.015 : Per-year recovery after LHM (low)
RECREEP_POEM : 0.012 : Per-year recovery after POEM (low)

// ---------- Disease physiology baseline (untreated achalasia steady-state) ----
nNOS_BASE  : 0.20  : Residual nNOS activity (fraction of normal)
VIP_BASE   : 0.25  : Residual VIP tone
LES_REST   : 40    : Untreated resting LES pressure (mmHg) — supranormal
IRP_BASE   : 25    : Baseline IRP (mmHg)
ESO_BASE   : 3.0   : Baseline esophageal diameter (cm), normal 2 cm
PERIST_BASE: 0.05  : Residual peristaltic integrity (frac of normal)

// Symptom progression (per year)
NEURON_DECAY : 0.05 : Annual fractional further inhibitory neuron loss
DIL_PROG   : 0.20  : Annual diameter growth (cm) if untreated
STASIS_PROD: 200   : Daily food/saliva input (mL/d)
STASIS_KE  : 0.05  : Stasis clearance (1/d) when functional emptying

// Symptom mapping (sigmoid sensitivity)
DYS_K      : 0.10  : LES Δ→dysphagia mapping
REG_K      : 0.40  : Dilatation→regurgitation mapping
CP_K       : 0.15  : LES pressurization→chest pain (Type III heavy)
WT_K       : 0.020 : Dysphagia→weight loss (kg/wk)

// AE thresholds
HA_EC50    : 0.25  : Headache severity threshold from cGMP elevation
HYPO_EC50  : 0.30  : Hypotension threshold from drug effect

// Misc
SCALE_T    : 1.0   : Time scaling

$CMT @annotated
ISDN_GUT  : ISDN sublingual absorption depot (mg)
ISDN_CEN  : ISDN central plasma (mg)
ISMN_CEN  : ISMN active metabolite (mg)
NIF_GUT   : Nifedipine gut depot (mg)
NIF_CEN   : Nifedipine plasma (mg)
SIL_GUT   : Sildenafil gut depot (mg)
SIL_CEN   : Sildenafil plasma (mg)
TAD_GUT   : Tadalafil gut depot (mg)
TAD_CEN   : Tadalafil plasma (mg)
BOTOX_LES : Botox local LES depot (U)
LES_PRESS : LES resting pressure (mmHg)
IRP_C     : Integrated relaxation pressure (mmHg)
ESO_DIL   : Esophageal diameter (cm)
ESO_STAS  : Food/saliva stasis volume (mL)
PERIST    : Peristaltic integrity (fraction)
DYS_S     : Dysphagia (0-3, Eckardt component)
REG_S     : Regurgitation (0-3)
CP_S      : Chest pain (0-3)
WT_S      : Weight-loss component (0-3)
TBE5      : 5-min TBE column height (cm)
QOL       : Achalasia HRQoL (0-100; 100 best)
AE_HA     : Drug headache severity (0-3)
AE_HYPO   : Drug-induced postural hypotension (0-3)

$MAIN
double F_ISDN  = ISDN_F;
double F_NIF   = NIF_F;
double F_SIL   = SIL_F;
double F_TAD   = TAD_F;

// CYP3A4 adjustment for nifedipine/sildenafil/tadalafil
double CL_NIF  = NIF_CL * CYP3A4;
double CL_SIL  = SIL_CL * CYP3A4;
double CL_TAD  = TAD_CL * CYP3A4;

// Initial conditions reflect untreated achalasia steady-state
LES_PRESS_0 = LES_REST;
IRP_C_0     = IRP_BASE;
ESO_DIL_0   = ESO_BASE + DIL_PROG * DURATION;
PERIST_0    = PERIST_BASE * (SUBTYPE == 3 ? 4 : 1); // some preserved in Type III
DYS_S_0     = 2.0;
REG_S_0     = (SUBTYPE == 1 || SUBTYPE == 2 ? 1.5 : 1.0);
CP_S_0      = (SUBTYPE == 3 ? 2.0 : 0.5);
WT_S_0      = 1.5;
TBE5_0      = 8.0;
QOL_0       = 45;

$ODE
// --------- Drug PK -----------------------------------------------------------
dxdt_ISDN_GUT = -ISDN_KA * ISDN_GUT;
dxdt_ISDN_CEN =  ISDN_KA * ISDN_GUT * F_ISDN - (ISDN_CL/ISDN_V) * ISDN_CEN;
dxdt_ISMN_CEN =  0.5 * (ISDN_CL/ISDN_V) * ISDN_CEN - ISDN_ISMN_CL/ISDN_V * ISMN_CEN;

dxdt_NIF_GUT  = -NIF_KA * NIF_GUT;
dxdt_NIF_CEN  =  NIF_KA * NIF_GUT * F_NIF - (CL_NIF/NIF_V) * NIF_CEN;

dxdt_SIL_GUT  = -SIL_KA * SIL_GUT;
dxdt_SIL_CEN  =  SIL_KA * SIL_GUT * F_SIL - (CL_SIL/SIL_V) * SIL_CEN;

dxdt_TAD_GUT  = -TAD_KA * TAD_GUT;
dxdt_TAD_CEN  =  TAD_KA * TAD_GUT * F_TAD - (CL_TAD/TAD_V) * TAD_CEN;

// Botox local LES depot, slow first-order decay
dxdt_BOTOX_LES = -BOT_KE * BOTOX_LES;

// --------- Plasma concentrations ---------------------------------------------
double Cp_ISDN = ISDN_CEN / ISDN_V * 1000.0;   // mg/L (multiplied → ng/mL effectively)
double Cp_ISMN = ISMN_CEN / 200 * 1000.0;
double Cp_NIF  = NIF_CEN / NIF_V * 1000.0;     // ng/mL approx
double Cp_SIL  = SIL_CEN / SIL_V * 1000.0;
double Cp_TAD  = TAD_CEN / TAD_V * 1000.0;

// --------- Drug PD on LES tone (fractional reduction) ------------------------
double E_ISDN = ISDN_EMAX * (Cp_ISDN + 0.6 * Cp_ISMN) / (ISDN_EC50 + Cp_ISDN + 0.6 * Cp_ISMN);
double E_NIF  = NIF_EMAX * Cp_NIF / (NIF_EC50 + Cp_NIF);
double E_SIL  = SIL_EMAX * Cp_SIL / (SIL_EC50 + Cp_SIL);
double E_TAD  = TAD_EMAX * Cp_TAD / (TAD_EC50 + Cp_TAD);
double E_BOT  = BOT_EMAX * BOTOX_LES / (BOT_EC50 + BOTOX_LES);

// Combined Bliss-independent reduction in LES tone (drugs only; procedural step
// applied separately as initial LES_PRESS adjustment via $MAIN dosing logic).
double drug_red = 1.0 - (1.0-E_ISDN) * (1.0-E_NIF) * (1.0-E_SIL) * (1.0-E_TAD) * (1.0-E_BOT);

// --------- LES pressure dynamics ---------------------------------------------
// Untreated LES tone slowly climbs as more inhibitory neurons are lost
double les_baseline = LES_REST + 5.0 * NEURON_DECAY * DURATION;
// Procedural adjustment encoded by EVID-2 events in regimen; here we capture
// slow recurrence/creep back to baseline after a procedure dropped LES_PRESS.
double recreep = 0.0;
if (PROC_PD   > 0.5) recreep = RECREEP_PD;
if (PROC_LHM  > 0.5) recreep = RECREEP_LHM;
if (PROC_POEM > 0.5) recreep = RECREEP_POEM;

double les_target = les_baseline * (1.0 - drug_red);
dxdt_LES_PRESS = 2.0 * (les_target - LES_PRESS)        // pharmacologic shift
               + recreep * (les_baseline - LES_PRESS); // post-proc creep back

// IRP tracks LES tone (slightly damped, target 0.65× LES_PRESS - 1)
double irp_target = 0.65 * LES_PRESS - 1.0;
dxdt_IRP_C = 1.5 * (irp_target - IRP_C);

// --------- Esophageal diameter and stasis ------------------------------------
double dil_rate = (LES_PRESS > 25.0 ? 0.0008 * (LES_PRESS - 25.0) : -0.001);
dxdt_ESO_DIL = dil_rate;                          // cm/day (slow)

double emptying_eff = 1.0 / (1.0 + exp(0.4 * (LES_PRESS - 18.0)));  // sigmoid
dxdt_ESO_STAS = STASIS_PROD - STASIS_KE * (1.0 + 3.0 * emptying_eff) * ESO_STAS;

// --------- Peristaltic integrity (mostly fixed, minimal recovery) -------------
double per_target = PERIST_BASE * (SUBTYPE == 3 ? 4.0 : 1.0);
dxdt_PERIST = 0.001 * (per_target - PERIST);     // essentially fixed

// --------- Symptom dynamics (Eckardt components) ------------------------------
double dys_target = 3.0 / (1.0 + exp(-DYS_K * (LES_PRESS - 15.0)));
double reg_target = 3.0 / (1.0 + exp(-REG_K * (ESO_DIL - 3.5))) +
                    0.5 * (ESO_STAS > 200 ? 1 : 0);
double cp_target  = (SUBTYPE == 3 ? 3.0 : 1.5) /
                    (1.0 + exp(-CP_K * (LES_PRESS - 15.0)));
double wt_target  = 3.0 / (1.0 + exp(-(dys_target - 1.5) * 1.5));

// Cap at 3 (Eckardt scale)
if (reg_target > 3.0) reg_target = 3.0;
if (cp_target  > 3.0) cp_target  = 3.0;

dxdt_DYS_S = 0.5 * (dys_target - DYS_S);
dxdt_REG_S = 0.5 * (reg_target - REG_S);
dxdt_CP_S  = 0.5 * (cp_target  - CP_S);
dxdt_WT_S  = 0.3 * (wt_target  - WT_S);

// TBE 5-min column tracks dilation and emptying
double tbe_target = 0.7 * ESO_DIL * (1.0 + 0.4 * emptying_eff) + 0.005 * ESO_STAS;
dxdt_TBE5 = 0.7 * (tbe_target - TBE5);

// QoL inversely linked to total Eckardt and proc-related GERD penalty
double eck_inst = DYS_S + REG_S + CP_S + WT_S;
double qol_target = 100.0 - 6.0 * eck_inst - 5.0 * (AE_HA + AE_HYPO);
if (qol_target < 0.0) qol_target = 0.0;
dxdt_QOL = 0.4 * (qol_target - QOL);

// AE dynamics — drugs cause AE proportional to plasma exposure
double ha_drive = ISDN_HEAD_E * (Cp_ISDN/100.0) + SIL_HEAD_E * (Cp_SIL/100.0);
double hypo_drive = NIF_HYPO_E * (Cp_NIF/100.0) + SIL_HYPO_E * (Cp_SIL/100.0)
                  + ISDN_HEAD_E * 0.3 * (Cp_ISDN/100.0);
dxdt_AE_HA   = 1.0 * (3.0 * ha_drive   / (HA_EC50   + ha_drive)   - AE_HA);
dxdt_AE_HYPO = 1.0 * (3.0 * hypo_drive / (HYPO_EC50 + hypo_drive) - AE_HYPO);

$TABLE
double ECKARDT = DYS_S + REG_S + CP_S + WT_S;
double CLIN_RESP  = (ECKARDT <= 3.0 ? 1 : 0);     // clinical remission flag
double IRP_OK     = (IRP_C   < 15.0 ? 1 : 0);     // manometric success
double TBE_OK     = (TBE5    < 5.0  ? 1 : 0);     // adequate emptying

$CAPTURE Cp_ISDN Cp_NIF Cp_SIL Cp_TAD ECKARDT CLIN_RESP IRP_OK TBE_OK drug_red emptying_eff
'

mod_ach <- mcode("achalasia_qsp", code)

# =============================================================================
# Procedural step helper — adjusts LES_PRESS instantaneously
# =============================================================================
apply_procedure <- function(mod, proc = c("none","PD","LHM","POEM"), subtype = 2) {
  proc <- match.arg(proc)
  if (proc == "none") return(mod)
  les_now <- as.numeric(init(mod)$LES_PRESS)
  if (is.na(les_now) || les_now == 0) les_now <- mod$LES_REST
  delta <- switch(proc,
                  PD   = mod$PD_DELTA,
                  LHM  = mod$LHM_DELTA,
                  POEM = mod$POEM_DELTA + ifelse(subtype == 3, mod$POEM_T3_BONUS, 0))
  new_les <- les_now * (1 - delta)
  mod %>% init(LES_PRESS = new_les, IRP_C = 0.65 * new_les - 1.0) %>%
    param(setNames(list(1),
                   switch(proc, PD = "PROC_PD",
                                LHM = "PROC_LHM",
                                POEM = "PROC_POEM")))
}

# =============================================================================
# Scenario list (10 treatment arms)
# =============================================================================
scenarios <- list(
  list(name = "S01_NoTx",   regimen = ev(amt = 0, cmt = 1, time = 0)),
  list(name = "S02_ISDN",   regimen = ev(amt = 10, cmt = "ISDN_GUT", ii = 8, addl = 90 * 3, time = 0)),
  list(name = "S03_Nifedipine",
       regimen = ev(amt = 20, cmt = "NIF_GUT",  ii = 8, addl = 90 * 3, time = 0)),
  list(name = "S04_Sildenafil",
       regimen = ev(amt = 50, cmt = "SIL_GUT",  ii = 8, addl = 90 * 3, time = 0)),
  list(name = "S05_Tadalafil",
       regimen = ev(amt = 20, cmt = "TAD_GUT",  ii = 24, addl = 360, time = 0)),
  list(name = "S06_Botox",
       regimen = ev(amt = 100, cmt = "BOTOX_LES",
                    ii = 24 * 180, addl = 3, time = 0)),
  list(name = "S07_PneumaticDilation", proc = "PD"),
  list(name = "S08_HellerMyotomy",     proc = "LHM"),
  list(name = "S09_POEM",              proc = "POEM"),
  list(name = "S10_ISDN_plus_Sildenafil",
       regimen = c(ev(amt = 10, cmt = "ISDN_GUT", ii = 8, addl = 90 * 3),
                   ev(amt = 25, cmt = "SIL_GUT",  ii = 8, addl = 90 * 3)))
)

# =============================================================================
# Simulation runner
# =============================================================================
run_scenario <- function(scn, mod = mod_ach, end_days = 720, subtype = 2) {
  m <- mod %>% param(SUBTYPE = subtype)
  if (!is.null(scn$proc)) {
    m <- apply_procedure(m, scn$proc, subtype)
    out <- mrgsim(m, end = end_days, delta = 1)
  } else {
    out <- mrgsim(m, events = scn$regimen, end = end_days, delta = 1)
  }
  as.data.frame(out) %>% transform(scenario = scn$name, subtype = subtype)
}

if (FALSE) {
  # ------- Example: simulate all 10 scenarios for Type II achalasia ----------
  library(dplyr); library(ggplot2)
  res_all <- do.call(rbind, lapply(scenarios, run_scenario, subtype = 2))

  ggplot(res_all, aes(time / 7, ECKARDT, colour = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 3, lty = 2) +
    labs(x = "Weeks", y = "Eckardt symptom score",
         title = "Achalasia QSP — Eckardt score trajectories (10 scenarios)") +
    theme_minimal()

  ggplot(res_all, aes(time / 30, LES_PRESS, colour = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = c(15, 35), lty = 3) +
    labs(x = "Months", y = "LES resting pressure (mmHg)") +
    theme_minimal()
}

# =============================================================================
# Calibration notes (anchors with literature targets)
# =============================================================================
# 1. Resting LES pressure target ranges:
#    - Healthy 10–35 mmHg; achalasia typically 35–60 mmHg
#    - Post-LHM ~12 mmHg; post-POEM ~10–14 mmHg; post-PD ~15–18 mmHg
#    [Eckardt 1992 GE; Boeckxstaens 2011 NEJM; Werner 2019 NEJM]
#
# 2. IRP target: <15 mmHg success threshold (Chicago v4.0)
#    [Yadlapati 2021 NGM]
#
# 3. Clinical success (Eckardt ≤3) at 2 years (Boeckxstaens 2011 NEJM):
#    - PD: 86% (recurrence rises after 5 y)
#    - LHM: 90%
#    - POEM (Werner 2019): 83% non-inferior to LHM
#    - POEM > PD (Ponds 2019): 92% vs 54% @ 2 y
#
# 4. Botox single session: ~70% clinical response at 1 mo, dropping to
#    32% at 6 mo, 17% at 1 y (Pasricha 1996; Annese 2000)
#    Re-dose mean every 6–9 months.
#
# 5. Nifedipine 10–20 mg SL: ↓LES pressure by 30–40% lasting 30–60 min
#    (Triadafilopoulos 1991 DDS; Bortolotti 1981 GE)
#
# 6. Sildenafil 50 mg PO: ↓LES pressure by 35% lasting ~2 h
#    (Bortolotti 2000 GE; Eherer 2002 Gut)
#
# 7. AE rates:
#    - Nitrate headache 30–60%, hypotension 5–10%
#    - Nifedipine pedal edema/flushing 15–25%
#    - PDE5 headache 15–25%
#    - PD perforation 1–3%
#    - POEM post-procedural reflux 30–50% (PPI required)
#
# 8. Disease progression untreated:
#    - Eckardt slow rise 0.3–0.6 pts/y; esophageal diameter +0.2 cm/y
#    - End-stage sigmoid Ø > 6 cm typically after 15–20 y
#    - SCC risk 16-fold (Sandler 1995)
#
# 9. Subtype response (Pandolfino 2008 GE; Rohof 2013 GE):
#    - Type II: best response across all therapies (~96% LHM/POEM)
#    - Type I: ~80% LHM/POEM
#    - Type III: spastic — POEM with long myotomy preferred (~92% vs LHM 70%)
# =============================================================================
