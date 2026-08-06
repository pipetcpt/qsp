# =============================================================================
#  ts_mrgsolve_model.R
#  THYROID STORM (갑상선 폭풍 / thyrotoxic crisis) — QSP model for mrgsolve
# =============================================================================
#
#  CENTRAL THESIS
#  --------------
#  Thyroid storm is not "more thyrotoxicosis".  Serum T4 and T3 do not separate
#  storm from uncomplicated thyrotoxicosis.  What separates them is whether a
#  fast positive-feedback loop is CLOSED, and the model writes that as a single
#  dimensionless ratio which happens to be a PRODUCT OF FACTORS:
#
#                Q_prod        80 W · M_thy · M_temp · M_sns · M_unc
#     Lambda =  --------  =  -----------------------------------------------
#                Q_loss       (8 + 52·E) · Vol^1.5 · (1 + cool) · (Tc − 33)
#
#  Every drug used in thyroid storm multiplies one or two of those factors,
#  so the clinical treatment hierarchy becomes arithmetic rather than a list:
#
#     M_thy   thionamide (slowly), iodide (release), PTU/glucocorticoid/
#             iopanoic acid (D1,D2), cholestyramine, plasma exchange
#             ... AND beta-blockade, via NEFA → displacement → free T3
#     M_sns   beta-blockade (a SECOND time)
#     M_unc   ASPIRIN raises it — uncoupling of oxidative phosphorylation
#     E       acetaminophen raises it (lowers the febrile setpoint);
#             CNS injury destroys it (ring C)
#     Vol     fluid resuscitation raises it
#     cool    external cooling raises it
#
#  THREE POSITIVE-FEEDBACK RINGS are represented explicitly:
#     ring A   fT3 → beta1 density & lipolytic machinery → NEFA
#                  → displacement from TBG → free fraction → fT3
#     ring B   Tc → Q10 amplification of metabolism → heat production → Tc
#     ring C   Tc → hypothalamic injury → loss of heat-loss effector → Tc
#
#  38 compartments.  Time unit = HOURS throughout.
#
#  This file is a line-by-line port of ts_verify_python.py, which is the
#  arithmetic of record for every number quoted in README.md.  Three deliberate
#  differences are noted where they occur:
#     (i)   plasma exchange is spread over a 2-hour first-order removal here,
#           rather than applied as an instantaneous fractional removal;
#     (ii)  the state clamps are implemented as flux limiters (the correct way
#           inside an ODE solver) rather than as post-step clipping;
#     (iii) chronic baselines are reached by TIME INTEGRATION with continuation
#           in TRAb (ts_thyrotoxic), never by a fixed-point solve -- the fast
#           subsystem is bistable and a fixed point can land on the hot branch.
#
#  Usage
#  -----
#    install.packages("mrgsolve"); install.packages("dplyr")
#    source("ts_mrgsolve_model.R")
#    res <- ts_run_all()          # all 18 scenarios
#    ts_plot(res)                 # quick-look panels
#    ts_table(res)                # the summary table from README section [1-18]
#
# =============================================================================

suppressMessages({
  library(mrgsolve)
  library(dplyr)
})

# =============================================================================
#  MODEL
# =============================================================================
ts_code <- '
$PROB
# Thyroid storm QSP model
- 38 compartments; time in hours
- Storm = loss of a stable operating point of the FAST subsystem

$PARAM @annotated
// ---------------- distribution volumes ----------------
V_T4    :  10.0 : T4 distribution volume (L)
V_T3    :  40.0 : T3 distribution volume (L)
V_rT3   :  30.0 : rT3 distribution volume (L)
V_I     :  17.0 : iodide distribution volume (L)

// ---------------- T4 disposal (per hour) ----------------
kD1_4   : 0.0011125 : D1 outer-ring T4->T3 (1/h)  [0.0267/24]
kD2_4   : 0.00055417: D2 outer-ring T4->T3 (1/h)  [0.0133/24]
kD3_4   : 0.00145833: D3 inner-ring T4->rT3 (1/h) [0.0350/24]
kEHC_4  : 0.00368750: biliary/conjugated T4 -> gut (1/h) [0.0885/24]
kreab   : 0.60      : fraction of biliary T4 reabsorbed
kgut    : 0.25      : gut transit (1/h)
krT3    : 0.1733    : rT3 elimination (1/h), t1/2 4 h
fD1_rT3 : 0.85      : fraction of rT3 clearance that is D1-mediated
kT3     : 0.028704  : T3 elimination (1/h), t1/2 24.1 h

// ---------------- thyroid gland ----------------
ksec0   : 6.25e-4   : fractional secretion of the colloid store (1/h)
fT3sec0 : 0.080     : basal fraction of secretion appearing as T3
fT3secM : 0.170     : stimulated gland T3 fraction
Vorg0   : 0.040     : organification Vmax at Stim = 1 (umol I/h)
Km_org  : 0.200     : intrathyroidal iodide half-saturating organification (umol)
kleak   : 0.012     : iodide efflux from the gland (1/h)
CLthy   : 1.121     : thyroidal iodide clearance at NIS = Stim = 1 (L/h)
CLren_I : 2.100     : renal iodide clearance (L/h)
I_diet  : 0.04790   : dietary iodide (umol/h) = 146 ug/d
nI_step : 1.0       : iodines released per deiodination step
nI_T3   : 3.0       : iodines released when T3/rT3 are cleared
IC_WC   : 2.00      : intrathyroidal iodide for half Wolff-Chaikoff (umol)
nWC     : 2.0       : Wolff-Chaikoff Hill exponent
ImWCrel : 0.78      : max inhibition of hormone RELEASE by iodide
ImWCorg : 0.90      : max inhibition of ORGANIFICATION by iodide
K_NIS   : 1.00      : plasma iodide for half NIS down-regulation (umol/L)
Km_NIS  : 25.0      : NIS Km for iodide -- transport SATURATES (umol/L)
tau_NIS : 84.0      : Wolff-Chaikoff escape time constant (h)

// ---------------- HPT axis ----------------
tau_TSH : 12.0      : TSH time constant (h)
TSH0    : 1.50      : euthyroid TSH (mIU/L)
kTSHfb  : 8.0       : TSH feedback gain
TRAb    : 1.00      : TSH-receptor stimulating drive (1 = normal)

// ---------------- binding and the free fraction (ring A) ----------------
phi4_0  : 2.0e-4    : euthyroid free T4 fraction
phi3_0  : 3.0e-3    : euthyroid free T3 fraction
NEFA0   : 0.40      : euthyroid NEFA (mmol/L)
EmNEFA  : 1.110     : max NEFA displacement
K_NEFA  : 0.862     : NEFA displacement K (mmol/L above baseline)
n_NEFA  : 2.0       : NEFA displacement Hill exponent (STEEP on purpose)
EmSal   : 1.20      : max salicylate displacement
K_Sal   : 20.0      : salicylate displacement K (mg/dL)
tau_TBG : 48.0      : TBG time constant (h)
eTBGpr  : 0.15      : acute-phase fall in TBG per unit precipitant

// ---------------- nuclear receptor signal ----------------
K_TR    : 1.0       : normalised fT3 for half TR occupancy

// ---------------- lipolysis (ring A) ----------------
tauNEFA : 0.25      : NEFA time constant (h)
eL_T3   : 0.50      : T3 effect on lipolytic capacity
eL_beta : 1.17      : beta effect on lipolysis

// ---------------- sympathetic outflow ----------------
// eT3_S = 0 ON PURPOSE: thyroid hormone does NOT raise circulating
// catecholamines.  The sensitisation is at the receptor (Rb), not the
// transmitter.  This is a modelling commitment, not an omission.
tau_SNS : 0.50      : sympathetic time constant (h)
ePrec_S : 0.70      : precipitant effect on sympathetic outflow
eTemp_S : 0.16      : temperature effect on sympathetic outflow (per degC)
eVol_S  : 0.35      : hypovolaemia effect on sympathetic outflow
eT3_S   : 0.00      : thyroid hormone effect on catecholamine LEVEL

// ---------------- beta1 receptor density ----------------
tau_Rb  : 36.0      : beta1 up-regulation time constant (h)
eR      : 1.00      : beta1 density gain per unit TR signal

// ---------------- heart rate ----------------
HR0     : 70.0      : euthyroid heart rate (bpm)
HRmax   : 210.0     : maximal sinus rate (bpm)
Kd_HR   : 1.205     : heart-rate drive K
eHR_b   : 0.55      : beta contribution to rate drive
eHR_T3  : 0.35      : direct SA-node (If) contribution
eHR_T   : 0.100     : temperature contribution (per degC)
tau_HR  : 0.30      : heart-rate time constant (h)

// ---------------- cardiac pump ----------------
SV0     : 1.10      : stroke-volume scale
eSV_T3  : 0.60      : T3 effect on contractility
eSV_b   : 0.50      : beta effect on contractility
W_crit  : 4.50      : cardiac work index above which reserve depletes
kdep    : 0.020     : reserve depletion rate (1/h per unit excess)
krec    : 0.010     : reserve recovery rate (1/h)

// ---------------- thermal (rings B and C) ----------------
BMR0    : 80.0      : euthyroid heat production (W)
C_body  : 245000.0  : body heat capacity (J/K)
eBMR    : 1.00      : TR signal effect on BMR
Q10     : 2.00      : van t Hoff Q10
eQ_sns  : 0.15      : beta contribution to thermogenesis
h_min   : 8.0       : fully vasoconstricted heat-loss coefficient (W/K)
h_span  : 52.0      : effector span of the heat-loss coefficient (W/K)
T_sink  : 33.0      : effective heat sink temperature (degC)
Tset0   : 37.42     : effector midpoint (degC)
w_eff   : 0.35      : effector width (degC)
ePGEset : 1.90      : febrile setpoint rise per unit PGE2 (degC)
ePGEshv : 0.15      : PGE2 contribution to thermogenesis
eCNSeff : 0.60      : loss of heat-loss effector per unit CNSx  (RING C GAIN)
eVoleff : 1.50      : exponent, h scales with Vol^eVoleff
EmUnc   : 0.50      : max salicylate uncoupling of oxidative phosphorylation
K_Unc   : 30.0      : salicylate uncoupling K (mg/dL)
Tc_ceil : 43.0      : numerical ceiling on core temperature (degC)

// ---------------- CNS (ring C) ----------------
kCNSon  : 0.022     : thermal CNS injury rate (per degC per h)
kCNShyp : 0.080     : hypoperfusion CNS injury rate (1/h)
kCNSoff : 0.045     : CNS recovery rate (1/h)
eCNS_T3 : 0.50      : TR signal amplification of thermal CNS injury
T_CNS   : 38.5      : temperature above which thermal CNS injury begins
perfCNS : 0.65      : perfusion below which hypoperfusion injury begins

// ---------------- liver ----------------
Bili0   : 0.80      : euthyroid bilirubin (mg/dL)
kb_out  : 0.030     : bilirubin turnover (1/h)
eB_hyp  : 6.0       : hypoperfusion effect on bilirubin
eB_hf   : 2.5       : congestion effect on bilirubin
eB_T3   : 0.35      : TR signal effect on bilirubin

// ---------------- cortisol ----------------
Cort0   : 400.0     : euthyroid cortisol (nmol/L)
kCLcort : 0.35      : cortisol clearance (1/h)
eCL_T3  : 1.67      : TR signal acceleration of cortisol clearance
eStrCor : 0.90      : stress increment in adrenal output
AR      : 1.00      : adrenal reserve (1 = intact)
nmolmgL : 2759.0    : nmol/L glucocorticoid activity per mg/L hydrocortisone-eq

// ---------------- GI ----------------
tau_GI  : 6.0       : GI dysfunction time constant (h)
eGI_T3  : 0.55      : TR signal effect on GI dysfunction
eGI_pr  : 0.50      : precipitant effect on GI dysfunction

// ---------------- volume ----------------
kdrink  : 0.120     : oral repletion rate (1/h)
h_base  : 20.04     : euthyroid heat-loss coefficient (W/K)
kins_h  : 0.00025   : insensible loss per (W/K) above baseline (1/h)
kins_T  : 0.00150   : insensible loss per degC above 37.5 (1/h)
kGIloss : 0.0055    : GI fluid loss (1/h per unit GIx)
ePO_GI  : 0.35      : GI symptoms blunting oral intake
VolFlr  : 0.55      : floor on volume status

// ---------------- atrial fibrillation ----------------
tau_AF  : 8.0       : AF burden time constant (h)
eAF_T3  : 0.55      : TR signal effect on AF
eAF_b   : 0.35      : beta effect on AF
AFsusc  : 1.00      : individual AF susceptibility

// ---------------- precipitant ----------------
kPrec   : 0.0110    : precipitant decay (1/h)
kPGE    : 0.0140    : PGE2 decay (1/h)
Prec0   : 0.00      : precipitant intensity at t = 0
PGE00   : 0.00      : PGE2 intensity at t = 0

// ---------------- mortality hazard ----------------
// h0_haz is the ONLY parameter in this model fitted to storm data: it is set
// so that the untreated fulminant arm gives 7-day mortality 85%.
h0_haz  : 2.755e-6  : hazard scale (1/h)
aBW     : 0.45      : hazard sensitivity to Burch-Wartofsky
eShock  : 3.00      : hazard multiplier for shock
eBili   : 1.20      : hazard multiplier for jaundice
eHyper  : 2.00      : hazard multiplier above 40.5 degC
eAI     : 0.90      : hazard multiplier for adrenal insufficiency

// ================= DRUG PK =================
ka_ptu  : 2.00 : PTU absorption (1/h)
ke_ptu  : 0.462: PTU elimination (1/h), t1/2 1.5 h
F_ptu   : 0.80 : PTU bioavailability
V_ptu   : 21.0 : PTU volume (L)
ka_mmi  : 1.60 : methimazole absorption (1/h)
ke_mmi  : 0.139: methimazole elimination (1/h), t1/2 5 h
F_mmi   : 0.93 : methimazole bioavailability
V_mmi   : 28.0 : methimazole volume (L)
ka_pro  : 1.40 : propranolol absorption (1/h)
ke_pro  : 0.173: propranolol elimination (1/h), t1/2 4 h
F_pro   : 0.30 : propranolol bioavailability (high first pass)
V_pro   : 280.0: propranolol volume (L)
ke_esm  : 4.62 : esmolol elimination (1/h), t1/2 9 min
ke_gc   : 0.408: hydrocortisone elimination (1/h), t1/2 1.7 h
kchol   : 0.25 : cholestyramine gut clearance (1/h)
ke_iop  : 0.0347: iopanoic acid elimination (1/h), t1/2 20 h
ke_asa  : 0.1155: salicylate elimination (1/h), t1/2 6 h
ke_apap : 0.277: acetaminophen elimination (1/h), t1/2 2.5 h
ESMRATE : 0.0  : esmolol infusion rate (ug/mL/h)

// ================= DRUG PD =================
IC50TPUp: 0.80 : PTU IC50 on TPO (mg/L)
IC50TPMm: 0.060: methimazole IC50 on TPO (mg/L)
Imax_TPO: 0.97 : max organification inhibition by a thionamide
IC50D1p : 2.50 : PTU IC50 on D1 (mg/L)  -- MMI has NO D1 activity
ImD1_ptu: 0.78 : PTU max D1 inhibition
IC50D1gc: 0.35 : glucocorticoid IC50 on D1 (mg/L)
ImD1_gc : 0.40 : glucocorticoid max D1 inhibition
IC50D2gc: 0.35 : glucocorticoid IC50 on D2 (mg/L)
ImD2_gc : 0.45 : glucocorticoid max D2 inhibition
IC50D1io: 1.50 : iopanoic acid IC50 on D1 (mg/L)
ImD1_iop: 0.92 : iopanoic acid max D1 inhibition
IC50D2io: 1.50 : iopanoic acid IC50 on D2 (mg/L)
ImD2_iop: 0.85 : iopanoic acid max D2 inhibition
IC50D1pr: 95.0 : propranolol IC50 on D1 (ng/mL)
ImD1_pro: 0.30 : propranolol max D1 inhibition
IC50bpro: 32.0 : propranolol IC50 on the beta receptor (ng/mL)
Imaxbpro: 0.90 : propranolol max beta blockade
IC50besm: 0.55 : esmolol IC50 on the beta receptor (ug/mL)
Imaxbesm: 0.90 : esmolol max beta blockade
ekchol  : 0.90 : max block of enterohepatic T4 reuptake
IC50chol: 4.00 : cholestyramine IC50 (g in gut)
eapapPGE: 0.65 : max acetaminophen suppression of PGE2
IC50apap: 12.0 : acetaminophen IC50 (mg/L)

// ================= SUPPORTIVE CARE (constant covariates) =================
COOL    : 0.0  : external cooling, fractional increase in h
FLUID   : 0.0  : intravenous fluid rate (fraction of volume status per h)

// ================= PLASMA EXCHANGE =================
TPE_ON  : 0.0  : 1 = perform a plasma exchange
TPE_T   : 12.0 : start time of the exchange (h)
TPE_DUR : 2.0  : duration of the exchange (h)
TPE_FR  : 0.65 : fraction of the PLASMA pool removed
V_plasma: 3.0  : plasma volume (L) -- only this fraction is exchangeable

$CMT @annotated
S     : thyroidal organified hormone store (nmol T4-eq)
T4    : plasma total T4 (nmol/L)
T3    : plasma total T3 (nmol/L)
rT3   : plasma reverse T3 (nmol/L)
Ggut  : enterohepatic T4 in gut lumen (nmol)
Ipl   : plasma inorganic iodide (umol/L)
Ithy  : intrathyroidal non-organified iodide (umol)
NIS   : NIS activity (relative)
TSH   : plasma TSH (mIU/L)
NEFA  : plasma non-esterified fatty acids (mmol/L)   [FAST, ring A]
TBG   : thyroxine-binding globulin (relative)
SNS   : sympathetic outflow (relative)               [FAST]
Rb    : beta1-adrenoceptor density (relative)
HR    : heart rate (bpm)                             [FAST]
CR    : cardiac contractile reserve (1 = full)
Tc    : core temperature (degC)                      [FAST, rings B and C]
Vol   : volume status (1 = euvolemic)                [FAST]
CNSx  : CNS / hypothalamic dysfunction index         [FAST, ring C]
Bili  : total bilirubin (mg/dL)
Cort  : plasma cortisol (nmol/L)
GIx   : GI dysfunction index
AFb   : atrial fibrillation burden
Hz    : cumulative mortality hazard
BWi   : integral of the Burch-Wartofsky score
Prec  : precipitant intensity
PGE   : PGE2-driven setpoint elevation
Aptu  : PTU in gut (mg)
Cptu  : PTU plasma (mg/L)
Ammi  : methimazole in gut (mg)
Cmmi  : methimazole plasma (mg/L)
Apro  : propranolol in gut (mg)
Cpro  : propranolol plasma (ng/mL)
Cesm  : esmolol plasma (ug/mL)
Cgc   : glucocorticoid, hydrocortisone-eq (mg/L)
Achol : cholestyramine in gut (g)
Ciop  : iopanoic acid plasma (mg/L)
Casa  : salicylate plasma (mg/dL)
Capap : acetaminophen plasma (mg/L)

$GLOBAL
#define _CLIP(x, lo, hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))

// fractional inhibition, 0..Imax
double emaxf(double C, double IC50, double Imax) {
  return (C > 0.0) ? (Imax * C / (IC50 + C)) : 0.0;
}

// derived quantities exported to $TABLE
double g_fT3, g_fT4, g_fT3n, g_Fd, g_phi3, g_Sig, g_Bsig, g_occb;
double g_Qprod, g_Qloss, g_h, g_E, g_Mthy, g_Mtemp, g_Msns, g_Munc;
double g_perf, g_shock, g_CHF, g_Widx, g_CO, g_Lam;
double g_inhD1, g_inhD2, g_inhTPO, g_ATPO, g_relinh, g_Synth, g_Sec, g_Stim;
double g_BW;

$MAIN
// Euthyroid initial condition (mass-balance consistent; see README section N)
S_0    = 8000.0;
T4_0   = 100.0;
T3_0   = 1.80;
rT3_0  = 0.2805;
Ggut_0 = 14.75;
Ipl_0  = 0.019840;
Ithy_0 = 0.20;
NIS_0  = 1.0;
TSH_0  = TSH0;
NEFA_0 = NEFA0;
TBG_0  = 1.0;
SNS_0  = 1.0;
Rb_0   = 1.0;
HR_0   = HR0;
CR_0   = 1.0;
Tc_0   = 37.0;
Vol_0  = 1.0;
CNSx_0 = 0.0;
Bili_0 = Bili0;
Cort_0 = Cort0;
GIx_0  = 0.0;
AFb_0  = 0.0;
Prec_0 = Prec0;
PGE_0  = PGE00;

$ODE
// ===========================================================================
// (a) FREE FRACTION -- the second axis.  Fatty-acid displacement is STEEP
//     (Hill 2), i.e. negligible in ordinary states and large only at the NEFA
//     levels reached in critical illness, heparin exposure, and storm.
// ===========================================================================
double dN   = (NEFA > NEFA0) ? (NEFA - NEFA0) : 0.0;
double dNn  = pow(dN, n_NEFA);
double Dn   = EmNEFA * dNn / (pow(K_NEFA, n_NEFA) + dNn);
double Dsal = emaxf(Casa, K_Sal, EmSal);
double TBGe = (TBG < 0.25) ? 0.25 : TBG;
g_Fd   = (1.0 + Dn + Dsal) / TBGe;
g_phi3 = phi3_0 * g_Fd;
g_fT3  = T3 * g_phi3 * 1000.0;                    // pmol/L
g_fT4  = T4 * phi4_0 * g_Fd * 1000.0;             // pmol/L
g_fT3n = g_fT3 / (1.8 * phi3_0 * 1000.0);

// (b) NUCLEAR TR SIGNAL -- saturating; = 1 when normal, -> 2 maximal
g_Sig = (g_fT3n / (K_TR + g_fT3n)) * (K_TR + 1.0);

// (c) BETA SIGNAL
g_occb = 1.0 - (1.0 - emaxf(Cpro, IC50bpro, Imaxbpro)) *
               (1.0 - emaxf(Cesm, IC50besm, Imaxbesm));
g_Bsig = SNS * Rb * (1.0 - g_occb);

// (d) DEIODINASE INHIBITION -- the FAST axis
g_inhD1 = 1.0 - (1.0 - emaxf(Cptu, IC50D1p , ImD1_ptu)) *
                (1.0 - emaxf(Cgc , IC50D1gc, ImD1_gc )) *
                (1.0 - emaxf(Ciop, IC50D1io, ImD1_iop)) *
                (1.0 - emaxf(Cpro, IC50D1pr, ImD1_pro));
g_inhD2 = 1.0 - (1.0 - emaxf(Cgc , IC50D2gc, ImD2_gc )) *
                (1.0 - emaxf(Ciop, IC50D2io, ImD2_iop));

// (e) THIONAMIDE / TPO -- the SLOW axis
g_inhTPO = 1.0 - (1.0 - emaxf(Cptu, IC50TPUp, Imax_TPO)) *
                 (1.0 - emaxf(Cmmi, IC50TPMm, Imax_TPO));

// (f) IODIDE: Wolff-Chaikoff acts on BOTH release and organification
double IthyN = (Ithy > 0.0) ? Ithy : 0.0;
double wcN   = pow(IthyN, nWC);
double wc    = wcN / (pow(IC_WC, nWC) + wcN);
g_ATPO   = (1.0 - g_inhTPO) * (1.0 - ImWCorg * wc);
g_relinh = ImWCrel * wc;

// (g) GLAND
g_Stim = TRAb + 0.6 * (TSH / TSH0 - 1.0);
if (g_Stim < 0.05) g_Stim = 0.05;
// Organification is ENZYME-limited, not first-order in iodide.  Without the
// Michaelis term an iodide load drives synthesis without bound and the gland
// refills straight through a thionamide block, which is wrong.  Vmax rises
// with stimulation because TPO/DUOX are cAMP-regulated.
double J_org = Vorg0 * g_Stim * g_ATPO * IthyN / (Km_org + IthyN);
g_Synth = J_org / 4.0 * 1000.0;                        // nmol T4-eq / h
g_Sec   = ksec0 * g_Stim * S * (1.0 - g_relinh);
double sr = (g_Stim - 1.0) / 5.0;
sr = _CLIP(sr, 0.0, 1.0);
double fT3sec = fT3sec0 + (fT3secM - fT3sec0) * sr;
// NIS is cAMP/TSH-regulated, so the thyroid SHARE of the plasma iodide pool
// rises with stimulation.  Without this the gland is dietary-iodine-limited.
// NIS SATURATES: written as a clearance times a Michaelis factor so that at
// physiological iodide it is exactly CLthy*NIS*Stim, while a pharmacologic
// load cannot drive the intrathyroidal pool without bound.
double CLthyE = CLthy * NIS * g_Stim / (1.0 + Ipl / Km_NIS);

// ===========================================================================
// (h) THERMAL -- Lambda = Q_prod / Q_loss is a PRODUCT OF FACTORS
// ===========================================================================
g_Mthy  = 1.0 + eBMR * (g_Sig - 1.0);           if (g_Mthy < 0.35) g_Mthy = 0.35;
g_Mtemp = pow(Q10, (Tc - 37.0) / 10.0);
g_Msns  = 1.0 + eQ_sns * (g_Bsig - 1.0) + ePGEshv * PGE;
if (g_Msns < 0.40) g_Msns = 0.40;
g_Munc  = 1.0 + emaxf(Casa, K_Unc, EmUnc);
g_Qprod = BMR0 * g_Mthy * g_Mtemp * g_Msns * g_Munc;

double PGEeff = PGE * (1.0 - emaxf(Capap, IC50apap, eapapPGE));
double Tset   = Tset0 + ePGEset * PGEeff;
double xE     = _CLIP(-(Tc - Tset) / w_eff, -40.0, 40.0);
g_E = 1.0 / (1.0 + exp(xE));
g_E *= (1.0 - eCNSeff * CNSx);                  // <== RING C
double VolE = (Vol < 0.40) ? 0.40 : Vol;
g_h     = (h_min + h_span * g_E) * pow(VolE, eVoleff) * (1.0 + COOL);
g_Qloss = g_h * (Tc - T_sink);
g_Lam   = g_Qprod / ((g_Qloss > 1e-9) ? g_Qloss : 1e-9);

// (i) PUMP / PERFUSION
double contract = SV0 * (1.0 + eSV_T3 * (g_Sig - 1.0) + eSV_b * (g_Bsig - 1.0));
if (contract < 0.15) contract = 0.15;
double CRe = (CR > 0.05) ? CR : 0.05;
double SV  = contract * CRe * ((Vol < 1.0) ? Vol : 1.0);
g_CO   = SV * HR / (SV0 * HR0);
double CO_need = g_Qprod / BMR0;
g_perf = g_CO / ((CO_need > 0.2) ? CO_need : 0.2);
if (g_perf > 1.6) g_perf = 1.6;
g_shock = (g_perf < 0.70) ? (1.0 - g_perf / 0.70) : 0.0;
g_Widx  = contract * HR / HR0;
g_CHF   = (CR < 1.0) ? (1.0 - CR) : 0.0;

// ===========================================================================
//  ODEs
// ===========================================================================
// (0) thyroid store -- the 2-month reservoir a thionamide cannot touch
dxdt_S = g_Synth - g_Sec;

// (1) plasma T4
double kD1 = kD1_4 * (1.0 - g_inhD1);
double kD2 = kD2_4 * (1.0 - g_inhD2);
double reab_blk = emaxf(Achol, IC50chol, ekchol);
// A 1-1.5 plasma-volume exchange removes TPE_FR of the PLASMA pool, not of
// the body.  Plasma is 30% of the T4 distribution volume but only 7.5% of
// the T3 volume, so the whole-body removal differs per hormone.  The model has no
// plasma / extravascular split, so only the net whole-body removal can be
// represented (see caveat L7 in README.md).
double tpe_on = (TPE_ON > 0.5 && SOLVERTIME >= TPE_T &&
                 SOLVERTIME < TPE_T + TPE_DUR) ? 1.0 : 0.0;
double tpe     = tpe_on * (-log(1.0 - TPE_FR * V_plasma / V_T4 ) / TPE_DUR);
double tpe_T3  = tpe_on * (-log(1.0 - TPE_FR * V_plasma / V_T3 ) / TPE_DUR);
double tpe_rT3 = tpe_on * (-log(1.0 - TPE_FR * V_plasma / V_rT3) / TPE_DUR);
dxdt_T4 = g_Sec * (1.0 - fT3sec) / V_T4
        + kgut * Ggut * kreab * (1.0 - reab_blk) / V_T4
        - (kD1 + kD2 + kD3_4 + kEHC_4) * T4
        - tpe * T4;

// (2) plasma T3 -- 80% of it is made OUTSIDE the thyroid
dxdt_T3 = (kD1 + kD2) * T4 * V_T4 / V_T3
        + g_Sec * fT3sec / V_T3
        - kT3 * T3
        - tpe_T3 * T3;

// (3) rT3 -- the observable signature of a D1 block
// rT3 is cleared mainly BY D1, so a D1 block raises it.  This is the
// laboratory signature that PTU -- and not methimazole -- is working on
// the fast axis.
double krT3e = krT3 * (1.0 - fD1_rT3 * g_inhD1);
dxdt_rT3 = kD3_4 * T4 * V_T4 / V_rT3 - krT3e * rT3 - tpe_rT3 * rT3;

// (4) enterohepatic T4 in the gut lumen
dxdt_Ggut = kEHC_4 * T4 * V_T4 - kgut * Ggut;

// (5) plasma inorganic iodide, WITH RECYCLING.  Every deiodination step
//     returns one iodine to the plasma pool and clearance of T3/rT3 returns
//     the remaining three.  Without this the gland could not sustain
//     thyrotoxicosis on an unchanged 150 ug/d intake.
double recyc = (nI_step * (kD1 + kD2 + kD3_4) * T4 * V_T4
              + nI_T3 * kT3 * T3 * V_T3
              + nI_T3 * krT3 * rT3 * V_rT3) / 1000.0;
dxdt_Ipl = (I_diet + recyc - CLren_I * Ipl - CLthyE * Ipl + kleak * IthyN) / V_I;

// (6) intrathyroidal iodide
dxdt_Ithy = CLthyE * Ipl - J_org - kleak * IthyN;

// (7) NIS down-regulation = WOLFF-CHAIKOFF ESCAPE
double NIS_t = 1.0 / (1.0 + Ipl / K_NIS);
dxdt_NIS = (NIS_t - NIS) / tau_NIS;
if (NIS >= 1.0 && dxdt_NIS > 0.0) dxdt_NIS = 0.0;

// (8) TSH
double TSH_t = _CLIP(TSH0 * exp(-kTSHfb * (g_Sig - 1.0)), 0.002, 60.0);
dxdt_TSH = (TSH_t - TSH) / tau_TSH;

// (9) NEFA == RING A, first half
double M_lip = (1.0 + eL_T3 * (g_Sig - 1.0)) * (1.0 + eL_beta * (g_Bsig - 1.0));
if (M_lip < 0.25) M_lip = 0.25;
dxdt_NEFA = (NEFA0 * M_lip - NEFA) / tauNEFA;

// (10) TBG
double Pc = (Prec < 1.0) ? Prec : 1.0;
dxdt_TBG = (1.0 - eTBGpr * Pc - TBG) / tau_TBG;

// (11) sympathetic outflow -- driven by the PRECIPITANT, not by T3
double dT_S  = (Tc  > 37.0) ? (Tc - 37.0)  : 0.0;
double dV_S  = (Vol < 1.0 ) ? (1.0 - Vol)  : 0.0;
double SNS_t = 1.0 + ePrec_S * Prec + eTemp_S * dT_S + eVol_S * dV_S
             + eT3_S * (g_Sig - 1.0);
dxdt_SNS = (SNS_t - SNS) / tau_SNS;

// (12) beta1 receptor density -- the T3-driven sensitisation
dxdt_Rb = (1.0 + eR * (g_Sig - 1.0) - Rb) / tau_Rb;

// (13) heart rate
double drive = eHR_b * (g_Bsig - 1.0) + eHR_T3 * (g_Sig - 1.0) + eHR_T * dT_S;
if (drive < -0.55) drive = -0.55;
double HR_t;
if (drive >= 0.0) HR_t = HR0 + (HRmax - HR0) * drive / (Kd_HR + drive);
else              HR_t = HR0 * (1.0 + drive * 0.55);
if (HR_t < 38.0) HR_t = 38.0;
dxdt_HR = (HR_t - HR) / tau_HR;

// (14) cardiac reserve
double over = (g_Widx / W_crit > 1.0) ? (g_Widx / W_crit - 1.0) : 0.0;
dxdt_CR = -kdep * over + ((over <= 0.0) ? (krec * (1.0 - CR)) : 0.0);
if (CR >= 1.0 && dxdt_CR > 0.0) dxdt_CR = 0.0;

// (15) CORE TEMPERATURE == RINGS B and C
dxdt_Tc = (g_Qprod - g_Qloss) * 3600.0 / C_body;
if (Tc >= Tc_ceil && dxdt_Tc > 0.0) dxdt_Tc = 0.0;

// (16) volume
double ins = kins_h * ((g_h > h_base) ? (g_h - h_base) : 0.0)
           + kins_T * ((Tc > 37.5) ? (Tc - 37.5) : 0.0);
double po  = (1.0 - CNSx) * (1.0 - ePO_GI * GIx);
if (po < 0.0) po = 0.0;
dxdt_Vol = kdrink * po * (1.0 - Vol) + FLUID - ins - kGIloss * GIx;
if (Vol <= VolFlr && dxdt_Vol < 0.0) dxdt_Vol = 0.0;
if (Vol >= 1.0   && dxdt_Vol > 0.0) dxdt_Vol = 0.0;

// (17) CNS / hypothalamic injury == RING C
double dT_C = (Tc > T_CNS) ? (Tc - T_CNS) : 0.0;
double dP_C = (g_perf < perfCNS) ? (perfCNS - g_perf) : 0.0;
double cns_t = (kCNSon * dT_C * (1.0 + eCNS_T3 * (g_Sig - 1.0))
              + kCNShyp * dP_C) / kCNSoff;
if (cns_t > 1.0) cns_t = 1.0;
dxdt_CNSx = (cns_t - CNSx) * kCNSoff;
if (CNSx >= 1.0 && dxdt_CNSx > 0.0) dxdt_CNSx = 0.0;

// (18) bilirubin
double dperf_B = (g_perf < 0.85) ? (0.85 - g_perf) : 0.0;
dxdt_Bili = kb_out * (Bili0 * (1.0 + eB_hyp * dperf_B + eB_hf * g_CHF
                              + eB_T3 * (g_Sig - 1.0)) - Bili);

// (19) cortisol -- clearance is accelerated BY thyroid hormone
dxdt_Cort = kCLcort * Cort0 * AR * (1.0 + eStrCor * Prec)
          - kCLcort * (1.0 + eCL_T3 * (g_Sig - 1.0)) * Cort;

// (20) GI dysfunction
double dS_G = (g_Sig > 1.0) ? (g_Sig - 1.0) : 0.0;
double GI_t = eGI_T3 * dS_G + eGI_pr * Prec;
if (GI_t > 1.0) GI_t = 1.0;
dxdt_GIx = (GI_t - GIx) / tau_GI;

// (21) atrial fibrillation burden
double dB_A = (g_Bsig > 1.0) ? (g_Bsig - 1.0) : 0.0;
double AF_t = AFsusc * (eAF_T3 * dS_G * 2.0 + eAF_b * dB_A);
if (AF_t > 1.0) AF_t = 1.0;
dxdt_AFb = (AF_t - AFb) / tau_AF;
if (AFb >= 1.0 && dxdt_AFb > 0.0) dxdt_AFb = 0.0;

// ---------------------------------------------------------------------------
// BURCH-WARTOFSKY POINT SCALE, computed from the simulated signs
// ---------------------------------------------------------------------------
double bwT;
if      (Tc < 37.2) bwT = 0.0;
else if (Tc < 37.8) bwT = 5.0;
else if (Tc < 38.3) bwT = 10.0;
else if (Tc < 38.9) bwT = 15.0;
else if (Tc < 39.4) bwT = 20.0;
else if (Tc < 40.0) bwT = 25.0;
else                bwT = 30.0;
double bwC;
if      (CNSx < 0.12) bwC = 0.0;
else if (CNSx < 0.42) bwC = 10.0;
else if (CNSx < 0.75) bwC = 20.0;
else                  bwC = 30.0;
double bwG;
if      (Bili >= 3.0) bwG = 20.0;
else if (GIx  > 0.25) bwG = 10.0;
else                  bwG = 0.0;
double bwH;
if      (HR < 99.0)  bwH = 0.0;
else if (HR < 110.0) bwH = 5.0;
else if (HR < 120.0) bwH = 10.0;
else if (HR < 130.0) bwH = 15.0;
else if (HR < 140.0) bwH = 20.0;
else                 bwH = 25.0;
double bwF;
if      (g_CHF < 0.08) bwF = 0.0;
else if (g_CHF < 0.20) bwF = 5.0;
else if (g_CHF < 0.38) bwF = 10.0;
else                   bwF = 15.0;
double bwA = (AFb > 0.5) ? 10.0 : 0.0;
double bwP = (Prec > 0.05) ? 10.0 : 0.0;
g_BW = bwT + bwC + bwG + bwH + bwF + bwA + bwP;

// (22) mortality hazard
double bili_x = (Bili > 1.2) ? (Bili - 1.2) : 0.0;
double hyp_x  = (Tc > 40.5) ? (Tc - 40.5) : 0.0;
// Effective glucocorticoid activity = endogenous cortisol PLUS the exogenous
// hydrocortisone-equivalent.  Without this the stress-dose steroid would
// inhibit D1/D2 but never repair the deficit it is actually given for.
double Cort_eff = Cort + nmolmgL * Cgc;
double cneed  = Cort0 * (1.0 + 0.9 * Prec);
double ai_x   = (Cort_eff < cneed) ? (1.0 - Cort_eff / cneed) : 0.0;
dxdt_Hz = h0_haz * exp(aBW * (g_BW - 45.0) / 10.0)
        * (1.0 + eShock * g_shock)
        * (1.0 + eBili * bili_x / 3.0)
        * (1.0 + eHyper * hyp_x)
        * (1.0 + eAI * ai_x);

// (23) integral of the Burch-Wartofsky score
dxdt_BWi = g_BW;

// (24,25) precipitant decay
dxdt_Prec = -kPrec * Prec;
dxdt_PGE  = -kPGE  * PGE;

// ================= DRUG PK =================
dxdt_Aptu  = -ka_ptu * Aptu;
dxdt_Cptu  =  ka_ptu * Aptu * F_ptu / V_ptu - ke_ptu * Cptu;
dxdt_Ammi  = -ka_mmi * Ammi;
dxdt_Cmmi  =  ka_mmi * Ammi * F_mmi / V_mmi - ke_mmi * Cmmi;
dxdt_Apro  = -ka_pro * Apro;
dxdt_Cpro  =  ka_pro * Apro * F_pro * 1000.0 / V_pro - ke_pro * Cpro;
dxdt_Cesm  =  ESMRATE - ke_esm * Cesm;
dxdt_Cgc   = -ke_gc * Cgc;
dxdt_Achol = -kchol * Achol;
dxdt_Ciop  = -ke_iop * Ciop;
dxdt_Casa  = -ke_asa * Casa;
dxdt_Capap = -ke_apap * Capap;

$TABLE
capture fT3    = g_fT3;
capture fT4    = g_fT4;
capture fT3n   = g_fT3n;
capture Fd     = g_Fd;
capture Sig    = g_Sig;
capture Bsig   = g_Bsig;
capture betaOcc= g_occb;
capture Qprod  = g_Qprod;
capture Qloss  = g_Qloss;
capture hcond  = g_h;
capture Eeff   = g_E;
capture M_thy  = g_Mthy;
capture M_temp = g_Mtemp;
capture M_sns  = g_Msns;
capture M_unc  = g_Munc;
capture Lambda = g_Lam;
capture perf   = g_perf;
capture shock  = g_shock;
capture CHF    = g_CHF;
capture Widx   = g_Widx;
capture BWPS   = g_BW;
capture mort   = 1.0 - exp(-Hz);
capture inhD1  = g_inhD1;
capture inhD2  = g_inhD2;
capture inhTPO = g_inhTPO;
capture relinh = g_relinh;
capture Synth  = g_Synth;
capture Secr   = g_Sec;
capture CortEff= Cort + nmolmgL * Cgc;
capture storeD = (g_Sec > 1e-9) ? (S / (g_Sec * 24.0)) : 0.0;
'

ts_mod <- mcode_cache("thyroid_storm", ts_code, atol = 1e-10, rtol = 1e-8,
                      maxsteps = 500000)

# =============================================================================
#  BASELINES
# =============================================================================
#  ts_euthyroid()   -- integrate to the euthyroid steady state
#  ts_thyrotoxic()  -- integrate to a chronic thyrotoxic steady state
#
#  NOTE ON BISTABILITY.  The fast subsystem has, besides the cool operating
#  point, a self-consistent HOT solution (Tc at the ceiling -> CNSx = 1 ->
#  heat-loss effector destroyed -> Tc at the ceiling).  Baselines must
#  therefore be reached by TIME INTEGRATION from a cool initial condition, and
#  thyrotoxic baselines by CONTINUATION from a less thyrotoxic one.  A
#  fixed-point solve would sometimes land on the hot branch and misreport
#  ordinary florid thyrotoxicosis as a storm.
# =============================================================================
ts_settle <- function(mod, hours = 9000, init = NULL, ...) {
  m <- mod
  if (!is.null(init)) m <- update(m, init = init)
  out <- m %>% param(...) %>%
    mrgsim(end = hours, delta = hours, recsort = 3) %>% as.data.frame()
  tail(out, 1)
}

ts_state <- function(row) {
  cmts <- c("S", "T4", "T3", "rT3", "Ggut", "Ipl", "Ithy", "NIS", "TSH",
            "NEFA", "TBG", "SNS", "Rb", "HR", "CR", "Tc", "Vol", "CNSx",
            "Bili", "Cort", "GIx", "AFb", "Hz", "BWi", "Prec", "PGE",
            "Aptu", "Cptu", "Ammi", "Cmmi", "Apro", "Cpro", "Cesm", "Cgc",
            "Achol", "Ciop", "Casa", "Capap")
  as.list(row[1, cmts])
}

ts_euthyroid <- function() ts_state(ts_settle(ts_mod, 9000, TRAb = 1.0))

#' Chronic thyrotoxic baseline, reached by continuation in TRAb.
ts_thyrotoxic <- function(TRAb = 4.19) {
  st <- ts_euthyroid()
  ladder <- unique(c(seq(1.5, min(TRAb, 4.0), by = 0.5), TRAb))
  for (tr in ladder) {
    st <- ts_state(ts_settle(ts_mod, 4000, init = st, TRAb = tr))
  }
  st$Hz <- 0; st$BWi <- 0          # reset the outcome accumulators
  st
}

# =============================================================================
#  REGIMENS
# =============================================================================
#  Doses are given in the units of the receiving compartment, exactly as in
#  ts_verify_python.py: gut compartments in mg, plasma compartments already
#  divided by their volume of distribution.
# =============================================================================
ev_ptu  <- function() c(ev(amt = 600, cmt = "Aptu", time = 0),
                        ev(amt = 250, cmt = "Aptu", time = 4, ii = 4, addl = 41))
ev_mmi  <- function() c(ev(amt = 40, cmt = "Ammi", time = 0),
                        ev(amt = 25, cmt = "Ammi", time = 6, ii = 6, addl = 27))
ev_pro  <- function() ev(amt = 80, cmt = "Apro", time = 0.25, ii = 4, addl = 41)
ev_iod  <- function(t0 = 1) ev(amt = 1970 / 17, cmt = "Ipl", time = t0,
                               ii = 6, addl = 27)   # SSKI 250 mg iodide q6h
ev_hc   <- function() ev(amt = 100 / 35, cmt = "Cgc", time = 0.25,
                         ii = 8, addl = 20)        # hydrocortisone 100 mg q8h
ev_chol <- function() ev(amt = 4, cmt = "Achol", time = 1, ii = 6, addl = 27)
ev_iop  <- function() ev(amt = 1000 / 25, cmt = "Ciop", time = 1, ii = 8, addl = 20)
ev_asa  <- function() ev(amt = 650 / 12 / 10, cmt = "Casa", time = 0.5,
                         ii = 4, addl = 41)        # aspirin 650 mg q4h
ev_apap <- function() ev(amt = 650 / 50, cmt = "Capap", time = 0.5,
                         ii = 6, addl = 27)        # acetaminophen 650 mg q6h

ts_none <- function() ev(amt = 0, cmt = "Aptu", time = 0)

# canonical fulminant precipitant (see README section [Bi])
TS_PREC <- 1.30
TS_PGE  <- 0.6 * TS_PREC

#' 18 scenarios, matching README section [1-18].
ts_scenarios <- function() list(
  `1  Untreated`                        = list(ev = ts_none(), p = list()),
  `2  PTU alone`                        = list(ev = ev_ptu(), p = list()),
  `3  Methimazole alone`                = list(ev = ev_mmi(), p = list()),
  `4  Iodide (1 h after PTU)`           = list(ev = c(ev_ptu(), ev_iod(1)),
                                               p = list()),
  `5  Propranolol alone`                = list(ev = ev_pro(), p = list()),
  `6  Esmolol alone`                    = list(ev = ts_none(),
                                               p = list(ESMRATE = 1.60)),
  `7  Hydrocortisone alone`             = list(ev = ev_hc(), p = list()),
  `8  Cooling + fluids alone`           = list(ev = ts_none(),
                                               p = list(COOL = 1.2,
                                                        FLUID = 0.0045)),
  `9  FULL BUNDLE`                      = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `10 BUNDLE with methimazole`          = list(
      ev = c(ev_mmi(), ev_iod(1), ev_pro(), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `11 BUNDLE, iodide BEFORE thionamide` = list(
      ev = c(ev_ptu(), ev_iod(0), ev_pro(), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `12 BUNDLE minus beta-blockade`       = list(
      ev = c(ev_ptu(), ev_iod(1), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `13 BUNDLE minus cooling/fluids`      = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()), p = list()),
  `14 BUNDLE + aspirin`                 = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc(), ev_asa()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `15 BUNDLE + acetaminophen`           = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc(), ev_apap()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `16 BUNDLE + iopanoic + cholestyramine` = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc(), ev_iop(), ev_chol()),
      p = list(COOL = 1.2, FLUID = 0.0045)),
  `17 BUNDLE + plasma exchange at 12 h` = list(
      ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045, TPE_ON = 1)),
  `18 BUNDLE with esmolol`              = list(
      ev = c(ev_ptu(), ev_iod(1), ev_hc()),
      p = list(COOL = 1.2, FLUID = 0.0045, ESMRATE = 1.60))
)

# =============================================================================
#  RUNNERS
# =============================================================================
ts_run <- function(spec, base = NULL, TRAb = 4.19, hours = 168, delta = 0.25,
                   prec = TS_PREC, pge = TS_PGE) {
  if (is.null(base)) base <- ts_thyrotoxic(TRAb)
  init <- base
  init$Prec <- prec
  init$PGE  <- pge
  m <- ts_mod %>% update(init = init) %>% param(TRAb = TRAb)
  if (length(spec$p)) m <- do.call(param, c(list(m), spec$p))
  m %>% mrgsim(events = spec$ev, end = hours, delta = delta, recsort = 3) %>%
    as.data.frame()
}

ts_run_all <- function(TRAb = 4.19, hours = 168) {
  base <- ts_thyrotoxic(TRAb)
  sc <- ts_scenarios()
  out <- lapply(names(sc), function(nm) {
    d <- ts_run(sc[[nm]], base = base, TRAb = TRAb, hours = hours)
    d$scenario <- nm
    d
  })
  bind_rows(out)
}

#' The summary table reproduced in README section [1-18].
ts_table <- function(res) {
  res %>% group_by(scenario) %>%
    summarise(
      Tmax      = max(Tc),
      Tc_24h    = Tc[which.min(abs(time - 24))],
      BWPS_24h  = BWPS[which.min(abs(time - 24))],
      BWPS_72h  = BWPS[which.min(abs(time - 72))],
      totT3_24h = T3[which.min(abs(time - 24))],
      freeT3_24 = fT3[which.min(abs(time - 24))],
      HR_24h    = HR[which.min(abs(time - 24))],
      runaway   = max(Tc) >= 40,
      mort_pct  = 100 * max(mort),
      .groups = "drop") %>%
    arrange(scenario) %>% as.data.frame()
}

# =============================================================================
#  DIAGNOSTIC SWEEPS -- the controlled comparisons behind the thesis
# =============================================================================

#' [H] Raise the hormone and NOTHING else.  No precipitant is ever present.
#' The prediction: no hormone level alone produces a storm.
ts_sweep_hormone <- function(trabs = c(1, 1.5, 2, 3, 4.19, 6, 9, 14, 22, 35),
                             hours = 72) {
  st <- ts_euthyroid()
  rows <- list()
  for (tr in trabs) {
    st <- ts_state(ts_settle(ts_mod, 4000, init = st, TRAb = tr))
    s2 <- st; s2$Hz <- 0; s2$BWi <- 0
    d <- ts_mod %>% update(init = s2) %>% param(TRAb = tr) %>%
      mrgsim(end = hours, delta = 0.25) %>% as.data.frame()
    rows[[length(rows) + 1]] <- data.frame(
      TRAb = tr, totT4 = tail(d$T4, 1), totT3 = tail(d$T3, 1),
      freeT3 = tail(d$fT3, 1), Fd = tail(d$Fd, 1), Sig = tail(d$Sig, 1),
      store_d = tail(d$storeD, 1), TSH = tail(d$TSH, 1),
      Tc = tail(d$Tc, 1), HR = tail(d$HR, 1), NEFA = tail(d$NEFA, 1),
      BWPS = tail(d$BWPS, 1), Lambda = tail(d$Lambda, 1),
      Tmax72 = max(d$Tc), storm = max(d$Tc) >= 40)
  }
  bind_rows(rows)
}

#' [Bi] Where is the boundary?  Sweep the precipitant at CONSTANT hormone.
ts_sweep_precipitant <- function(precs = c(0, .3, .6, .9, 1.1, 1.2, 1.25,
                                           1.3, 1.5, 2.0),
                                 TRAb = 4.19, hours = 168) {
  base <- ts_thyrotoxic(TRAb)
  bind_rows(lapply(precs, function(pr) {
    d <- ts_run(list(ev = ts_none(), p = list()), base = base, TRAb = TRAb,
                hours = hours, prec = pr, pge = 0.6 * pr)
    data.frame(Prec = pr, Tmax = max(d$Tc),
               Tc_24h = d$Tc[which.min(abs(d$time - 24))],
               HRmax = max(d$HR), BWPS_24h = d$BWPS[which.min(abs(d$time - 24))],
               BWPSmax = max(d$BWPS), Bilimax = max(d$Bili),
               mort_pct = 100 * max(d$mort),
               storm = max(d$Tc) >= 40)
  }))
}

#' [P] The propranolol trap: beta-blockade of a rate-dependent heart.
ts_sweep_reserve <- function(CRs = c(1.0, 0.8, 0.6, 0.45, 0.35),
                             TRAb = 4.19, hours = 168) {
  base <- ts_thyrotoxic(TRAb)
  arms <- list(
    propranolol = list(ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()),
                       p = list(COOL = 1.2, FLUID = 0.0045)),
    esmolol     = list(ev = c(ev_ptu(), ev_iod(1), ev_hc()),
                       p = list(COOL = 1.2, FLUID = 0.0045, ESMRATE = 1.60)))
  bind_rows(lapply(CRs, function(cr) {
    b <- base; b$CR <- cr
    bind_rows(lapply(names(arms), function(a) {
      d <- ts_run(arms[[a]], base = b, TRAb = TRAb, hours = hours)
      data.frame(CR0 = cr, arm = a, CO_nadir = min(d$perf),
                 shock_max = max(d$shock),
                 BWPS_24h = d$BWPS[which.min(abs(d$time - 24))],
                 Bilimax = max(d$Bili), mort_pct = 100 * max(d$mort))
    }))
  }))
}

#' [W] Wolff-Chaikoff escape: iodide is a bridge, not a treatment.
ts_wolff_chaikoff <- function(TRAb = 4.19, days = 21) {
  base <- ts_thyrotoxic(TRAb)
  d <- ts_run(list(ev = ev_iod(1), p = list()), base = base, TRAb = TRAb,
              hours = 24 * days, delta = 1, prec = 0, pge = 0)
  d %>% filter(time %% 24 == 0) %>%
    transmute(day = time / 24, Ipl, NIS, Ithy, relinh, A_TPO = 1 - inhTPO,
              Synth, S, totT4 = T4, totT3 = T3)
}

#' [F] The second axis: a doubled free fraction IS a doubled hormone, and the
#' total-hormone assay is completely blind to it.
ts_free_fraction <- function(TRAb = 4.19,
                             nefas = c(.4, .6, .8, 1.0, 1.5, 2.0, 2.5, 3.0)) {
  base <- ts_thyrotoxic(TRAb)
  bind_rows(lapply(nefas, function(nf) {
    b <- base; b$NEFA <- nf
    d <- ts_mod %>% update(init = b) %>% param(TRAb = TRAb) %>%
      mrgsim(end = 0, delta = 1) %>% as.data.frame()
    data.frame(NEFA = nf, measured_totT3 = d$T3[1], Fd = d$Fd[1],
               freeT3 = d$fT3[1], Sig = d$Sig[1],
               equivalent_totT3 = d$T3[1] * d$Fd[1])
  }))
}

# =============================================================================
#  QUICK-LOOK PLOTS
# =============================================================================
ts_plot <- function(res,
                    keep = c("1  Untreated", "2  PTU alone",
                             "3  Methimazole alone", "5  Propranolol alone",
                             "8  Cooling + fluids alone", "9  FULL BUNDLE")) {
  d <- res[res$scenario %in% keep, ]
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  vars <- list(c("Tc", "core temperature (degC)"),
               c("HR", "heart rate (bpm)"),
               c("fT3", "FREE T3 (pmol/L)"),
               c("T3", "TOTAL T3 (nmol/L)"),
               c("BWPS", "Burch-Wartofsky score"),
               c("mort", "cumulative mortality"))
  cols <- seq_along(keep)
  for (v in vars) {
    plot(NA, xlim = range(d$time), ylim = range(d[[v[1]]]), xlab = "time (h)",
         ylab = v[2], main = v[2])
    for (i in seq_along(keep)) {
      s <- d[d$scenario == keep[i], ]
      lines(s$time, s[[v[1]]], col = cols[i], lwd = 2)
    }
    if (v[1] == "Tc") abline(h = 40, lty = 3)
    if (v[1] == "BWPS") abline(h = 45, lty = 3)
  }
  plot.new()
  legend("center", legend = keep, col = cols, lwd = 2, bty = "n", cex = 0.8)
  par(op)
  invisible(NULL)
}

# =============================================================================
#  ONE-CALL DEMONSTRATION
# =============================================================================
ts_demo <- function() {
  cat("\n== [H] hormone sweep, NO precipitant ",
      "(prediction: nothing storms) ==\n", sep = "")
  print(ts_sweep_hormone(), row.names = FALSE, digits = 4)
  cat("\n== [Bi] precipitant sweep at CONSTANT hormone ",
      "(prediction: a threshold) ==\n", sep = "")
  print(ts_sweep_precipitant(), row.names = FALSE, digits = 4)
  cat("\n== [1-18] treatment scenarios ==\n")
  res <- ts_run_all()
  print(ts_table(res), row.names = FALSE, digits = 4)
  cat("\n== [F] the second axis ==\n")
  print(ts_free_fraction(), row.names = FALSE, digits = 4)
  cat("\n== [W] Wolff-Chaikoff escape ==\n")
  print(ts_wolff_chaikoff(), row.names = FALSE, digits = 4)
  cat("\n== [P] the propranolol trap ==\n")
  print(ts_sweep_reserve(), row.names = FALSE, digits = 4)
  invisible(res)
}

# Not run on source():
#   res <- ts_demo()
#   ts_plot(res)
