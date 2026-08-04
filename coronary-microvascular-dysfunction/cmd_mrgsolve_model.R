## =============================================================================
##  cmd_mrgsolve_model.R
##  Coronary Microvascular Dysfunction (CMD) / ANOCA-INOCA
##  Quantitative Systems Pharmacology model — mrgsolve implementation
## =============================================================================
##
##  THE STRUCTURAL IDEA
##  -------------------
##  Every index of the coronary microcirculation used in the clinic is a RATIO
##  or a DIFFERENCE of two quantities, and this model refuses to merge them:
##
##      CFR   = v_hyper / v_rest          two flows, one number
##      MR    = Pd / v                    resistance, not flow
##      DPTI  = (Pd - LVEDP) * t_dia      supply is bought only in diastole
##      SPTI  = P_sys * t_sys             demand is spent only in systole
##      def_L = MVO2_L - v_L*k*CaO2*Emax  the only quantity a patient can feel
##
##  Four consequences the model is built to EXPOSE rather than assume:
##
##   (1) A CFR of 1.9 is either a raised DENOMINATOR (resting flow too high:
##       the "functional" endotype) or a floored NUMERATOR (hyperaemic flow
##       capped by structure: "structural").  Rahman et al. Circulation
##       2019;140:1805 measured both; the two groups have the same CFR and
##       their absolute maximal flows differ by 57%.
##   (2) Resting tone is not a free parameter — the metabolic controller sets
##       it so that resting supply meets resting demand.  Given the hyperaemic
##       MR (which fixes RMIN) the RESTING MR is therefore predicted.  The
##       prediction lands on the measured value for the control and structural
##       groups.  It fails for the functional group, and the size of the
##       failure (AUTO_OFF) is this model's quantitative definition of what
##       functional CMD is.
##   (3) Heart rate enters the oxygen balance TWICE with the same sign: it
##       raises demand through the tension-time product and it shortens the
##       diastolic window in which the subendocardium is perfused.  43% of the
##       benefit of slowing 68 -> 55 bpm at fixed workload is the time window.
##   (4) An on-target drug can produce a null trial three ways: internal
##       cancellation (zibotentan / PRIZE), population dilution (ranolazine /
##       RWISE), and background contamination plus a timescale mismatch
##       (statin + ACEi / WARRIOR).  All three are reproduced.
##
##  And one prediction that was not designed in: because the functional
##  endotype's MINIMAL resistance is normal, it has almost no subendocardial
##  supply-demand deficit at any workload, so its angina cannot be ischaemic in
##  the supply-demand sense and must be carried by the nociceptive arm
##  (A1-adenosine afferent signalling plus central sensitisation).  That is why
##  anti-ischaemic drugs fail in it, and why aminophylline — an adenosine
##  receptor antagonist — is predicted to help precisely the group in which
##  every vasodilator fails (Elliott, Heart 1997;77:523, PMID 9227295).
##
##  CALIBRATION
##  -----------
##  Seven constants are solved (not hand-tuned) against seven published
##  targets; everything else is fixed a priori.  The solved values below were
##  produced by `cmd_reference_model.py` section 0 and are hard-coded here:
##
##      RMIN0        = 1.5174   control hyperaemic MR   = 2.20 mmHg/(cm/s)
##      W_ENDO       = 0.6283   control hyperaemic endo/epi = 1.02
##      F_ADO_REV    = 0.8760   functional hyperaemic MR = 2.30
##      RMIN_F(str)  = 1.4306   structural hyperaemic MR = 3.60
##      AUTO_OFF(fn) = 0.9165   functional resting MR    = 4.20
##      K_ANG        = 3.3997   untreated functional SAQ = 55
##      NOCI_THRESH  = 3.8575   untreated functional Bruce = 480 s
##
##  Rahman's three functional-group numbers are internally consistent
##  (4.2 / 2.30 = 1.83 = their CFR) because CFR = MR_rest/MR_hyp identically
##  when both are measured at the same aortic pressure, so calibrating two of
##  them fixes the third: the model's functional CFR is not free.
##
##  RELATIONSHIP TO THE PYTHON REFERENCE
##  ------------------------------------
##  There was no R runtime available when this model was built, so every
##  equation was implemented and executed first in dependency-free Python
##  (`cmd_reference_model.py`), which also carries the full bug log — 25 real
##  defects the numerical work exposed, each marked at the offending line in
##  both files.  Two deliberate differences:
##
##    * The arteriolar control loop settles in ~15 s.  Here TONE_E and TONE_P
##      are genuine ODEs with TAU_TONE = 15 s, integrated by LSODA, which is
##      built for that stiffness.  The Python reference solves the identical
##      expressions as a fixed point by bisection.  Section V of
##      `cmd_reference_output.txt` shows the two agree to 1e-6 in every
##      endotype and from either side.  Every index reported (CFR, MR,
##      endo/epi, DPTI, deficit) is by definition a steady state of that loop,
##      so the choice cannot change a reported number.
##      (Python bug B11: a damped fixed-point iteration at this loop gain
##      oscillated between the clamps and silently returned garbage; the
##      continuous-time ODE was never unstable.)
##    * Daily ischaemic burden.  Here the workload is a real diurnal covariate
##      WL(t) and the burden is a 24 h first-order filter (BI, BN) of the
##      instantaneous deficit and afferent drive.  The Python reference
##      evaluates the same activity distribution by quadrature.  The two agree
##      at steady state because tone equilibrates in 15 s, i.e. instantly on
##      the scale of the activity blocks.
##
##  UNITS.  Time h.  Pressure mmHg.  Velocity cm/s.  Resistance mmHg/(cm/s).
##  MVO2 mL O2/min/100 g.  Perfusion mL/min/g.  Concentrations mg/L.
##
##  COMPARTMENTS.  60 ODEs = 29 physiology + 31 PK/metabolite.
##  DRUGS.  14: ranolazine, ivabradine, bisoprolol, nebivolol, amlodipine,
##  diltiazem, nicorandil, zibotentan, ramipril, atorvastatin, fasudil,
##  sildenafil, aminophylline, imipramine — plus adenosine and acetylcholine
##  as diagnostic agents (computed algebraically in $TABLE, not dosed).
##
##  ⚠ Educational / research model.  Not validated for clinical use.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

cmd_code <- '
$PARAM @annotated
// ---- patient descriptors (the endotype lives here) -------------------------
RMIN_F   : 1.000 : Structural floor multiplier on minimal resistance (-)
AUTO_OFF : 0.000 : Metabolic controller offset, the functional lesion (-)
INFL     : 1.000 : Systemic inflammatory drive (-)
GENO     : 0     : rs9349379-G allele carrier (0/1)
ROCK_D   : 1.000 : Rho-kinase drive (-)
KSBP     : 16.0  : Exercise systolic gain, mmHg per workload unit
SPASM    : 0.000 : Spasm susceptibility (0-1)
SENS_BAS : 0.050 : Primary central-sensitisation drive (-)

// ---- systemic haemodynamics ----------------------------------------------
HR0      : 68.0  : Resting heart rate (bpm)
SBP0     : 132.0 : Resting systolic pressure (mmHg)
DBP0     : 76.0  : Resting diastolic pressure (mmHg)
PV       : 8.0   : Coronary venous pressure (mmHg)
KHR_WL   : 30.0  : Heart rate per workload unit (bpm)
KDBP_WL  : 5.0   : Diastolic pressure per workload unit (mmHg)
FR_DIA_P : 0.35  : Mean diastolic aortic P = DBP + FR_DIA_P*(SBP-DBP)
FR_SYS_P : 0.75  : Mean systolic aortic P  = DBP + FR_SYS_P*(SBP-DBP)
TSYS_K   : 0.30  : Systolic duration coefficient, t_sys = TSYS_K*sqrt(60/HR)
PHI_SYS  : 0.45  : Systolic perfusion efficiency, subepicardium (-)
FFR0     : 0.95  : Epicardial conductance, non-obstructive disease (-)

// ---- microvascular resistance --------------------------------------------
RMIN0    : 1.5174 : Minimal territory resistance, CALIBRATED, mmHg/(cm/s)
A_TONE   : 0.55   : r/r_max = 1 - A_TONE*TONE, so R scales as (.)^-4
W_ENDO   : 0.6283 : Relative subendocardial resistance, CALIBRATED (-)
W_EPI    : 1.15   : Relative subepicardial resistance (-)
RCOMP_K  : 0.0318 : Series compression resistance per mmHg of LVEDP
TONE_MAX : 0.95   : Maximal arteriolar activation (-)
SPASM_PD : 0.50   : Distal pressure lost per unit SPASM during an episode (-)
SPASM_SH : 0.40   : Fraction of the bed occluded per unit SPASM (-)
P_SP_DAY : 0.06   : Fraction of the day in spontaneous spasm per unit SPASM

// ---- oxygen transport ----------------------------------------------------
KCONV    : 0.0610 : Perfusion (mL/min/g) per unit velocity (cm/s)
CAO2     : 0.19   : Arterial oxygen content (mL O2 per mL blood)
E_REST   : 0.70   : Autoregulated resting oxygen extraction (-)
E_MAX    : 0.80   : Maximal oxygen extraction (-)

// ---- myocardial oxygen demand --------------------------------------------
MVO2_B0  : 1.50   : Basal MVO2 (mL O2/min/100 g)
K_TTI    : 0.0882 : MVO2 per bpm per unit normalised wall stress
K_CTR    : 1.50   : MVO2 per unit contractility index
K_DIAST  : 0.60   : MVO2 per unit diastolic tension
SIGMA_RF : 132.0  : Systolic pressure normalising wall stress to 1.0 (mmHg)
F_END_DM : 1.15   : Subendocardial demand multiplier (-)
F_EPI_DM : 0.88   : Subepicardial demand multiplier (-)

// ---- autoregulation ------------------------------------------------------
TONE_REF : 0.50    : Controller intercept (inside the adenosine-sensitive term)
G_AUTO   : 6.0     : Metabolic controller gain (-)
K_MYO    : 0.25    : Myogenic gain on (Pd-90)/90
TAU_TONE : 0.004167: Arteriolar tone time constant (h) = 15 s
K_ET_TN  : 0.12    : Tone per unit ET-1 above baseline
K_ROCK_TN: 0.10    : Tone per unit Rho-kinase above baseline
K_NO_TN  : 0.12    : Tone removed per unit NO above baseline
K_A1_TN  : 0.18    : alpha1 tone per unit sympathetic drive
F_ADO_REV: 0.8760  : Constrictor tone reversed by adenosine, CALIBRATED (-)

// ---- endothelial / constrictor biology ------------------------------------
KROS_BAS : 0.35   : Baseline ROS production coefficient
KBH4_OX  : 0.55   : BH4 oxidation per unit excess ROS
KET_SYN  : 1.45   : ET-1 synthesis coefficient
KET_CLR  : 1.0    : ET-1 clearance rate coefficient
A_INFL_E : 0.55   : ET-1 synthesis per unit inflammation
A_GENE_E : 0.30   : ET-1 synthesis in rs9349379-G carriers
B_NO_ET  : 0.45   : NO suppression of ET-1 synthesis
F_ETA    : 0.75   : Share of ET-1 tone mediated by ETA (-)
F_ETB2   : 0.25   : Share of ET-1 tone mediated by ETB2 (-)
TAU_ROS  : 12.0   : ROS pool time constant (h)
TAU_BH4  : 24.0   : BH4 pool time constant (h)
TAU_NO   : 6.0    : NO bioavailability time constant (h)
TAU_ET1  : 12.0   : ET-1 pool time constant (h)
TAU_ROCK : 12.0   : Rho-kinase activity time constant (h)

// ---- structure (slow) ----------------------------------------------------
TAU_ML   : 4320.0 : Arteriolar media/lumen time constant (h) = 180 d
TAU_CAPD : 5760.0 : Capillary density time constant (h) = 240 d
TAU_PVF  : 4800.0 : Perivascular fibrosis time constant (h) = 200 d
TAU_ICF  : 5760.0 : Interstitial fibrosis time constant (h) = 240 d
TAU_LVH  : 4320.0 : LV mass time constant (h) = 180 d
W_ML     : 0.55   : Minimal-resistance weight on media/lumen
W_PVF    : 0.35   : Minimal-resistance weight on perivascular fibrosis

// ---- filling pressure / volume -------------------------------------------
LVEDP0   : 11.0   : Baseline LV end-diastolic pressure (mmHg)
LVEDP_MX : 30.0   : Ceiling on LV end-diastolic pressure (mmHg)
TAU_LVDP : 24.0   : LVEDP time constant (h)
K_VOL_LP : 9.0    : LVEDP per unit fractional plasma volume gain
K_CAD_LP : 3.2    : LVEDP per unit diastolic calcium overload
K_ICF_LP : 6.0    : LVEDP per unit interstitial fibrosis
TAU_VOL  : 72.0   : Plasma volume time constant (h)
K_LVDP_W : 1.60   : Exercise LVEDP rise per workload unit (mmHg)
K_ICF_LW : 2.50   : Amplification of that rise by interstitial fibrosis

// ---- late sodium current / diastolic tension ------------------------------
TAU_NAI  : 12.0   : Intracellular sodium time constant (h)
TAU_CAD  : 12.0   : Diastolic calcium time constant (h)
K_ISCH_N : 0.85   : Late I_Na drive per unit burden (saturating)
K_ISCH_K : 1.50   : Half-saturation of that drive
K_TEN_CD : 0.55   : Diastolic tension per unit calcium overload

// ---- nociception / symptoms ----------------------------------------------
K_NOC_DF : 2.50   : Afferent drive per unit subendocardial deficit
K_NOC_AD : 0.90   : Afferent drive per unit adenosine signal
A_ADO_DF : 4.00   : Adenosine signal per unit deficit (saturating)
A_ADO_OF : 3.00   : Adenosine signal per unit controller offset x excess demand
MVO2_REF : 9.60   : Reference resting MVO2 (mL O2/min/100 g)
TAU_BURD : 24.0   : Daily burden filter time constant (h)
K_ANG    : 3.3997 : Angina episodes/week per unit afferent drive, CALIBRATED
K_ANG_GN : 0.85   : Amplification of angina by central sensitisation
K_ANG_CN : 8.00   : Angina generated by sensitisation alone
TAU_ANG  : 72.0   : Angina rate time constant (h)
TAU_SENS : 504.0  : Central sensitisation time constant (h) = 21 d
K_SENS_ON: 0.35   : Sensitisation gain on angina rate
K_SYMP_K : 12.0   : Half-saturation of sympathetic drive on angina rate
TAU_SYMP : 48.0   : Sympathetic tone time constant (h)
K_SYMP_HR: 0.22   : Fraction of resting heart rate recruited at SYMP = 1
SAQ_MAX  : 48.0   : Maximal SAQ decrement from angina
SAQ_K    : 3.4    : Half-saturating angina rate for the SAQ decrement
SAQ_SENS : 14.0   : SAQ decrement per unit central sensitisation
TAU_SAQ  : 336.0  : SAQ time constant (h) = 14 d
NOCI_THR : 3.8575 : Afferent drive at the anginal threshold, CALIBRATED
MCID_SAQ : 10.0   : Clinically important SAQ difference (U)

// ---- biomarkers / outcome ------------------------------------------------
KTNI_OFF : 0.35   : hs-cTnI equilibration rate (1/h)
TNI_MAX  : 14.0   : hs-cTnI ceiling above baseline (ng/L)
KBNP_ON  : 6.0    : NT-proBNP production per mmHg of LVEDP above baseline
KBNP_OFF : 0.02   : NT-proBNP elimination (1/h)
H0_MORT  : 0.0155 : Baseline mortality hazard (per year)
H0_HOSP  : 0.075  : Baseline angina-hospitalisation hazard (per year)
K_HOSP_A : 0.105  : Hospitalisation hazard per angina episode/week
LN116    : 0.14842: log(1.16), Kelshiker 2022 mortality HR per -0.1 CFR
LN108    : 0.07696: log(1.08), MACE HR per -0.1 CFR

// ---- drug PK -------------------------------------------------------------
KA_RAN : 0.45 : ranolazine absorption (1/h)
KE_RAN : 0.099: ranolazine elimination (1/h)
V_RAN  : 180.0: ranolazine volume (L)
EC_RAN : 1.60 : ranolazine EC50 (mg/L)
KA_IVA : 1.20 : ivabradine absorption (1/h)
KE_IVA : 0.116: ivabradine elimination (1/h)
V_IVA  : 110.0: ivabradine volume (L)
EC_IVA : 0.021: ivabradine EC50 (mg/L)
KF_IVA : 0.42 : S18982 metabolite formation (1/h)
KM_IVA : 0.155: S18982 metabolite elimination (1/h)
KA_BIS : 0.90 : bisoprolol absorption (1/h)
KE_BIS : 0.058: bisoprolol elimination (1/h)
V_BIS  : 230.0: bisoprolol volume (L)
EC_BIS : 0.012: bisoprolol EC50 (mg/L)
KA_NEB : 1.10 : nebivolol absorption (1/h)
KE_NEB : 0.069: nebivolol elimination (1/h)
V_NEB  : 1500.0 : nebivolol volume (L)
EC_NEB : 0.0016: nebivolol EC50 (mg/L)
KA_AML : 0.35 : amlodipine absorption (1/h)
KE_AML : 0.019: amlodipine elimination (1/h)
V_AML  : 1400.0 : amlodipine volume (L)
EC_AML : 0.0055: amlodipine EC50 (mg/L)
KA_DIL : 0.80 : diltiazem absorption (1/h)
KE_DIL : 0.139: diltiazem elimination (1/h)
V_DIL  : 350.0: diltiazem volume (L)
EC_DIL : 0.115: diltiazem EC50 (mg/L)
KA_NIC : 1.60 : nicorandil absorption (1/h)
KE_NIC : 0.578: nicorandil elimination (1/h)
V_NIC  : 70.0 : nicorandil volume (L)
EC_NIC : 0.055: nicorandil EC50 (mg/L)
KA_ZIB : 0.70 : zibotentan absorption (1/h)
KE_ZIB : 0.058: zibotentan elimination (1/h)
V_ZIB  : 95.0 : zibotentan volume (L)
EC_ZIB : 0.048: zibotentan EC50 (mg/L)
KA_RAM : 1.00 : ramipril absorption (1/h)
KE_RAM : 0.050: ramipril elimination (1/h)
V_RAM  : 100.0: ramipril volume (L)
EC_RAM : 0.0021: ramiprilat EC50 (mg/L)
KF_RAM : 0.62 : ramiprilat formation (1/h)
KM_RAM : 0.038: ramiprilat elimination (1/h)
KA_ATO : 0.60 : atorvastatin absorption (1/h)
KE_ATO : 0.050: atorvastatin elimination (1/h)
V_ATO  : 560.0: atorvastatin volume (L)
EC_ATO : 0.0085: atorvastatin EC50 (mg/L)
KA_FAS : 1.40 : fasudil absorption (1/h)
KE_FAS : 0.347: fasudil elimination (1/h)
V_FAS  : 90.0 : fasudil volume (L)
EC_FAS : 0.090: hydroxyfasudil EC50 (mg/L)
KF_FAS : 0.85 : hydroxyfasudil formation (1/h)
KM_FAS : 0.29 : hydroxyfasudil elimination (1/h)
KA_SIL : 1.30 : sildenafil absorption (1/h)
KE_SIL : 0.173: sildenafil elimination (1/h)
V_SIL  : 105.0: sildenafil volume (L)
EC_SIL : 0.135: sildenafil EC50 (mg/L)
KA_AMI : 1.10 : aminophylline absorption (1/h)
KE_AMI : 0.087: aminophylline elimination (1/h)
V_AMI  : 35.0 : aminophylline volume (L)
EC_AMI : 8.00 : aminophylline EC50 (mg/L)
KA_IMI : 0.90 : imipramine absorption (1/h)
KE_IMI : 0.041: imipramine elimination (1/h)
V_IMI  : 1200.0 : imipramine volume (L)
EC_IMI : 0.020: imipramine EC50 (mg/L)

// ---- maximal drug effects ------------------------------------------------
EM_IVA_HR : 0.19  : ivabradine, fraction of resting heart rate removed
EM_BIS_HR : 0.21  : bisoprolol, fraction of resting heart rate removed
EM_NEB_HR : 0.18  : nebivolol, fraction of resting heart rate removed
EM_DIL_HR : 0.09  : diltiazem, fraction of resting heart rate removed
EM_RAN_HR : 0.052 : ranolazine, fraction of resting heart rate removed
EM_BIS_CT : 0.17  : bisoprolol, contractility
EM_NEB_CT : 0.14  : nebivolol, contractility
EM_DIL_CT : 0.11  : diltiazem, contractility
EM_AML_TN : 0.36  : amlodipine, arteriolar tone
EM_DIL_TN : 0.24  : diltiazem, arteriolar tone
EM_NIC_TN : 0.33  : nicorandil, arteriolar tone
EM_SIL_TN : 0.16  : sildenafil, arteriolar tone
EM_AML_MP : 0.075 : amlodipine, mean arterial pressure
EM_DIL_MP : 0.055 : diltiazem, mean arterial pressure
EM_NIC_MP : 0.030 : nicorandil, mean arterial pressure
EM_ZIB_MP : 0.055 : zibotentan, mean arterial pressure
EM_RAM_MP : 0.070 : ramiprilat, mean arterial pressure
EM_SIL_MP : 0.045 : sildenafil, mean arterial pressure
EM_ZIB_ET : 0.86  : zibotentan, ETA receptor blockade
EM_ZIB_VL : 0.115 : zibotentan, fractional plasma volume gain
EM_FAS_RK : 0.72  : hydroxyfasudil, Rho-kinase
EM_ATO_RK : 0.22  : atorvastatin, Rho-kinase
EM_ATO_RS : 0.30  : atorvastatin, oxidative stress
EM_RAM_RS : 0.18  : ramiprilat, oxidative stress
EM_NEB_NO : 0.22  : nebivolol, NO bioavailability (beta3)
EM_NIC_NO : 0.20  : nicorandil, NO bioavailability
EM_RAM_NO : 0.14  : ramiprilat, NO bioavailability
EM_RAN_LN : 0.62  : ranolazine, late sodium current
EM_AMI_A1 : 0.68  : aminophylline, A1 afferent signalling
EM_AMI_A2 : 0.45  : aminophylline, A2A dilator reserve
EM_IMI_SN : 0.45  : imipramine, central sensitisation off-rate
EM_ATO_ML : 0.16  : atorvastatin, inward remodelling drive
EM_RAM_ML : 0.34  : ramiprilat, inward remodelling drive
EM_RAM_PV : 0.25  : ramiprilat, perivascular fibrosis drive
EM_RHB_CD : 0.28  : exercise rehabilitation, capillary density

// ---- non-pharmacological switches / diagnostics --------------------------
REHAB    : 0 : Exercise rehabilitation programme (0/1)
CBT      : 0 : Cognitive behavioural therapy for chest pain (0/1)
TMZ      : 0 : Trimetazidine, oxygen-sparing metabolic shift (0/1)
NOBPDROP : 0 : Counterfactual: suppress the drug blood-pressure fall (0/1)
NOFLUID  : 0 : Counterfactual: suppress drug fluid retention (0/1)
WL_FIX   : 0 : If > 0, hold workload at this value instead of the diurnal
               profile (used for provocation and sensitivity runs)

$CMT @annotated
// ---- endothelium and constrictors (5) ------------------------------------
NO     : NO bioavailability, normalised to healthy (-)
ROS    : Oxidative stress, normalised (-)
BH4    : Tetrahydrobiopterin, normalised (-)
ET1    : Plasma/interstitial endothelin-1, normalised (-)
ROCK   : Rho-kinase activity, normalised (-)
// ---- structure (5) -------------------------------------------------------
ML     : Arteriolar media/lumen excess (-)
CAPD   : Capillary density, normalised (-)
PVF    : Perivascular fibrosis (-)
ICF    : Interstitial fibrosis / ECV excess (-)
LVH    : LV mass index, normalised (-)
// ---- filling pressure and volume (2) -------------------------------------
LVEDP  : LV end-diastolic pressure (mmHg)
VOL    : Plasma volume, normalised (-)
// ---- late sodium current axis (2) ----------------------------------------
NAI    : Intracellular sodium, normalised (-)
CAD    : Diastolic calcium overload, normalised (-)
// ---- symptoms (4) --------------------------------------------------------
SENS   : Central sensitisation (0-1)
SYMP   : Sympathetic tone (0-1)
ANG    : Angina episode rate (per week)
SAQ    : Seattle Angina Questionnaire summary score (0-100)
// ---- biomarkers (2) ------------------------------------------------------
TNI    : High-sensitivity cardiac troponin I (ng/L)
BNP    : NT-proBNP (pg/mL)
// ---- hazards and exposure (5) --------------------------------------------
CHMORT : Cumulative mortality hazard (-)
CHMACE : Cumulative MACE hazard (-)
CHHOSP : Cumulative angina-hospitalisation hazard (-)
AUCDEF : Cumulative subendocardial oxygen deficit (mL O2/100 g)
AUCNOC : Cumulative afferent drive (arbitrary units x h)
// ---- fast arteriolar loop (2) --------------------------------------------
TONE_E : Subendocardial arteriolar activation (0-1)
TONE_P : Subepicardial arteriolar activation (0-1)
// ---- daily burden filters (2) --------------------------------------------
BI     : 24 h filtered subendocardial deficit (mL O2/min/100 g)
BN     : 24 h filtered afferent drive (-)
// ---- pharmacokinetics (31) ----------------------------------------------
A_RAN  : ranolazine absorption site (mg)
C_RAN  : ranolazine plasma (mg/L)
A_IVA  : ivabradine absorption site (mg)
C_IVA  : ivabradine plasma (mg/L)
M_IVA  : S18982 active metabolite (mg/L)
A_BIS  : bisoprolol absorption site (mg)
C_BIS  : bisoprolol plasma (mg/L)
A_NEB  : nebivolol absorption site (mg)
C_NEB  : nebivolol plasma (mg/L)
A_AML  : amlodipine absorption site (mg)
C_AML  : amlodipine plasma (mg/L)
A_DIL  : diltiazem absorption site (mg)
C_DIL  : diltiazem plasma (mg/L)
A_NIC  : nicorandil absorption site (mg)
C_NIC  : nicorandil plasma (mg/L)
A_ZIB  : zibotentan absorption site (mg)
C_ZIB  : zibotentan plasma (mg/L)
A_RAM  : ramipril absorption site (mg)
C_RAM  : ramipril plasma (mg/L)
M_RAM  : ramiprilat active metabolite (mg/L)
A_ATO  : atorvastatin absorption site (mg)
C_ATO  : atorvastatin plasma (mg/L)
A_FAS  : fasudil absorption site (mg)
C_FAS  : fasudil plasma (mg/L)
M_FAS  : hydroxyfasudil active metabolite (mg/L)
A_SIL  : sildenafil absorption site (mg)
C_SIL  : sildenafil plasma (mg/L)
A_AMI  : aminophylline absorption site (mg)
C_AMI  : aminophylline plasma (mg/L)
A_IMI  : imipramine absorption site (mg)
C_IMI  : imipramine plasma (mg/L)

$GLOBAL
#define CMD_POS(a) ((a) > 0.0 ? (a) : 0.0)

// Emax with a Hill coefficient of 1
double cmd_emax(double conc, double ec50, double emax) {
  if (conc <= 0.0) return 0.0;
  return emax * conc / (ec50 + conc);
}

// Systolic duration (s).  Shortens with rate, but sub-proportionally, which is
// exactly why heart rate costs the subendocardium diastolic time.
double cmd_tsys(double hr) {
  double tcyc = 60.0 / hr;
  double ts   = 0.30 * sqrt(60.0 / hr);
  if (ts > 0.62 * tcyc) ts = 0.62 * tcyc;
  return ts;
}

// Diurnal activity profile: 15 h at rest, 6 h light, 2.5 h moderate, 0.5 h
// stair-climbing.  The Python reference integrates the same distribution by
// quadrature.  bug B22: an earlier profile topped out below the workload at
// which a CFR-2 patient becomes ischaemic, so ordinary daily exertion
// generated no burden at all and the structural endotype came out
// asymptomatic.
double cmd_wl(double t) {
  double h = fmod(t, 24.0);
  if (h < 0.0) h += 24.0;
  if (h <  7.0)  return 1.00;   // sleep
  if (h <  8.0)  return 2.40;   // morning activity
  if (h < 11.0)  return 1.55;
  if (h < 12.0)  return 2.40;
  if (h < 16.0)  return 1.00;
  if (h < 19.0)  return 1.55;
  if (h < 19.5)  return 2.40;
  if (h < 20.0)  return 3.10;   // stairs
  return 1.00;
}

// Layer velocity (cm/s).  layer 0 = subendocardium, 1 = subepicardium.
// bug B3: the subendocardium is perfused ONLY in diastole and against LVEDP,
// and carries a SERIES compression resistance that does not scale with tone;
// the subepicardium is perfused in systole too, at reduced efficiency.  With
// the diastolic fraction applied to both layers it cancels out of the endo/epi
// ratio and heart rate becomes irrelevant to transmural distribution.
double cmd_vlayer(int layer, double Rart, double fdia, double Pd, double Pdsys,
                  double pv, double lvedp, double w_endo, double w_epi,
                  double rcomp_k, double phi_sys) {
  if (layer == 0) {
    double Rl = Rart * w_endo + rcomp_k * lvedp;
    return fdia * CMD_POS(Pd - lvedp) / Rl;
  }
  double Rl = Rart * w_epi;
  return (fdia * CMD_POS(Pd - pv)
          + (1.0 - fdia) * phi_sys * CMD_POS(Pdsys - pv)) / Rl;
}

$MAIN
// ---------------------------------------------------------------------------
// Initial conditions.  The endothelial block has a closed-form steady state
// that does not depend on flow, so the patient starts exactly at it and any
// drift below is a drug effect rather than an initialisation transient.
// ---------------------------------------------------------------------------
double ros0 = KROS_BAS * (1.0 + 0.55 * (INFL - 1.0)) * 2.857;
double bh40 = 1.0 / (1.0 + KBH4_OX * CMD_POS(ros0 - 1.0));
double no0  = bh40 / ros0;
double et0  = KET_SYN * (1.0 + A_INFL_E * (INFL - 1.0) + A_GENE_E * GENO)
              / ((1.0 + B_NO_ET * no0) * KET_CLR);

ROS_0   = ros0;
BH4_0   = bh40;
NO_0    = no0;
ET1_0   = et0;
ROCK_0  = ROCK_D;
ML_0    = RMIN_F - 1.0;
CAPD_0  = 1.0 / (1.0 + 0.45 * (RMIN_F - 1.0));
PVF_0   = 0.6 * (RMIN_F - 1.0);
ICF_0   = 0.25 * (INFL - 1.0);
LVH_0   = 1.0 + 0.10 * (INFL - 1.0);
LVEDP_0 = LVEDP0 + K_ICF_LP * 0.25 * (INFL - 1.0);
VOL_0   = 1.0;
NAI_0   = 1.0;
CAD_0   = 1.0;
SENS_0  = SENS_BAS / (SENS_BAS + 1.0);
SYMP_0  = 0.0;
ANG_0   = 0.0;
SAQ_0   = 70.0;
TNI_0   = 3.0;
BNP_0   = 60.0;
TONE_E_0 = 0.45;
TONE_P_0 = 0.45;

$ODE
// ===========================================================================
// 1.  DRUG EFFECTS
// ===========================================================================
double e_iva_hr = cmd_emax(C_IVA + 0.7 * M_IVA, EC_IVA, EM_IVA_HR);
double e_bis_hr = cmd_emax(C_BIS, EC_BIS, EM_BIS_HR);
double e_neb_hr = cmd_emax(C_NEB, EC_NEB, EM_NEB_HR);
double e_dil_hr = cmd_emax(C_DIL, EC_DIL, EM_DIL_HR);
double e_ran_hr = cmd_emax(C_RAN, EC_RAN, EM_RAN_HR);
double f_HR = e_iva_hr + e_bis_hr + e_neb_hr + e_dil_hr + e_ran_hr;

double f_CTR = cmd_emax(C_BIS, EC_BIS, EM_BIS_CT)
             + cmd_emax(C_NEB, EC_NEB, EM_NEB_CT)
             + cmd_emax(C_DIL, EC_DIL, EM_DIL_CT);

double e_ami_a2 = cmd_emax(C_AMI, EC_AMI, EM_AMI_A2);
double f_TONE = cmd_emax(C_AML, EC_AML, EM_AML_TN)
              + cmd_emax(C_DIL, EC_DIL, EM_DIL_TN)
              + cmd_emax(C_NIC, EC_NIC, EM_NIC_TN)
              + cmd_emax(C_SIL, EC_SIL, EM_SIL_TN)
              - 0.20 * e_ami_a2;      // aminophylline also blocks A2A dilatation

double e_zib_eta = cmd_emax(C_ZIB, EC_ZIB, EM_ZIB_ET);
double f_MAP = cmd_emax(C_AML, EC_AML, EM_AML_MP)
             + cmd_emax(C_DIL, EC_DIL, EM_DIL_MP)
             + cmd_emax(C_NIC, EC_NIC, EM_NIC_MP)
             + cmd_emax(C_ZIB, EC_ZIB, EM_ZIB_MP)
             + cmd_emax(M_RAM, EC_RAM, EM_RAM_MP)
             + cmd_emax(C_SIL, EC_SIL, EM_SIL_MP);
if (NOBPDROP > 0.5) f_MAP = 0.0;
if (f_MAP > 0.25)   f_MAP = 0.25;

double d_NO  = cmd_emax(C_NEB, EC_NEB, EM_NEB_NO)
             + cmd_emax(C_NIC, EC_NIC, EM_NIC_NO)
             + cmd_emax(M_RAM, EC_RAM, EM_RAM_NO);
double f_A1  = cmd_emax(C_AMI, EC_AMI, EM_AMI_A1);
double e_fas_rk = cmd_emax(M_FAS, EC_FAS, EM_FAS_RK);
double e_ato_rk = cmd_emax(C_ATO, EC_ATO, EM_ATO_RK);

// anti-spasm cover: calcium antagonists and Rho-kinase blockade
double anti_sp = cmd_emax(C_AML, EC_AML, EM_AML_TN) / EM_AML_TN
               + cmd_emax(C_DIL, EC_DIL, EM_DIL_TN) / EM_DIL_TN
               + e_fas_rk / EM_FAS_RK;
if (anti_sp > 1.0) anti_sp = 1.0;
double d_SPASM = 0.34 * (1.0 - anti_sp);
double f_MVO2  = (TMZ > 0.5) ? 0.055 : 0.0;

// ===========================================================================
// 2.  SYSTEMIC HAEMODYNAMICS AT THE CURRENT WORKLOAD
// ===========================================================================
double wl = (WL_FIX > 0.0) ? WL_FIX : cmd_wl(SOLVERTIME);

double HR = HR0 * (1.0 - f_HR) * (1.0 + K_SYMP_HR * SYMP) + KHR_WL * (wl - 1.0);
if (HR < 38.0)  HR = 38.0;
if (HR > 190.0) HR = 190.0;

double SBP = (SBP0 + KSBP * (wl - 1.0)) * (1.0 - f_MAP)
             * (1.0 + 0.35 * K_SYMP_HR * SYMP);
double DBP  = (DBP0 + KDBP_WL * (wl - 1.0)) * (1.0 - f_MAP);
double Pdia = DBP + FR_DIA_P * (SBP - DBP);
double Psys = DBP + FR_SYS_P * (SBP - DBP);
double tcyc = 60.0 / HR;
double tsys = cmd_tsys(HR);
double fdia = (tcyc - tsys) / tcyc;

// exercise filling pressure (bug B22): the dominant subendocardial insult
double lvedp_w = LVEDP + K_LVDP_W * (1.0 + K_ICF_LW * ICF) * CMD_POS(wl - 1.0);

// ===========================================================================
// 3.  DEMAND
// ===========================================================================
double sigma = SBP / SIGMA_RF;
double ctr   = (1.0 + 0.50 * (wl - 1.0)) * (1.0 - f_CTR);
double ten   = 1.0 + K_TEN_CD * (CAD - 1.0);
double mvo2  = (MVO2_B0 + K_TTI * HR * sigma + K_CTR * ctr + K_DIAST * ten)
               * (1.0 - f_MVO2);
double dem_e = mvo2 * F_END_DM;
double dem_p = mvo2 * F_EPI_DM;
double kO2   = KCONV * 100.0 * CAO2;

// required velocity per layer.  AUTO_OFF is the controller defect that DEFINES
// the functional endotype: the arterioles behave as though demand were
// AUTO_OFF x 100% higher than it is (bug B1 - a working controller cannot
// produce a low resting MR beside a normal minimal MR any other way).
double vreq_e = dem_e / (kO2 * E_REST) * (1.0 + AUTO_OFF);
double vreq_p = dem_p / (kO2 * E_REST) * (1.0 + AUTO_OFF);

// ===========================================================================
// 4.  RESISTANCE AND THE TWO-LAYER FLOW
// ===========================================================================
double capd = (CAPD < 0.35) ? 0.35 : CAPD;
double structf = (1.0 + W_ML * ML + W_PVF * PVF) / capd;

// spontaneous spasm: a fraction of the day in which a series obstruction and a
// bed occlusion appear.  bug B20: written as extra arteriolar TONE the
// controller simply absorbed it and a vasospastic patient could not become
// ischaemic at all.
double sev   = SPASM * (1.0 - anti_sp) * (1.0 + 0.5 * CMD_POS(ROCK - 1.0));
double shut  = SPASM_SH * sev;  if (shut > 0.80) shut = 0.80;
double p_sp  = P_SP_DAY * SPASM * (1.0 - anti_sp);
if (p_sp > 0.5) p_sp = 0.5;

// constrictor tone that the metabolic controller does not generate.  Only
// (1 - F_ADO_REV) of it survives exogenous adenosine (bug B2 and bug B14).
double et_eff = (ET1 - 1.0) * (F_ETA * (1.0 - e_zib_eta) + F_ETB2);
double tp = K_ET_TN * CMD_POS(et_eff)
          + K_ROCK_TN * CMD_POS(ROCK - 1.0)
          + K_A1_TN * SYMP * (1.0 - f_A1)
          - K_NO_TN * (NO - 1.0 + d_NO) - f_TONE;
if (tp < 0.0) tp = 0.0;
if (tp > TONE_MAX) tp = TONE_MAX;

double myo = K_MYO * (Pdia - 90.0) / 90.0;

// --- non-spasm state ---
double Rmin_ns = RMIN0 * structf;
double Pd_ns   = Pdia * FFR0;
double Pds_ns  = Psys * FFR0;
double Rart_e  = Rmin_ns / pow(1.0 - A_TONE * TONE_E, 4.0);
double Rart_p  = Rmin_ns / pow(1.0 - A_TONE * TONE_P, 4.0);
double v_e_ns  = cmd_vlayer(0, Rart_e, fdia, Pd_ns, Pds_ns, PV, lvedp_w,
                            W_ENDO, W_EPI, RCOMP_K, PHI_SYS);
double v_p_ns  = cmd_vlayer(1, Rart_p, fdia, Pd_ns, Pds_ns, PV, lvedp_w,
                            W_ENDO, W_EPI, RCOMP_K, PHI_SYS);

// --- spasm state (same tone, obstructed bed and lost distal pressure) ---
double ffr_sp  = 1.0 - SPASM_PD * sev;  if (ffr_sp < 0.25) ffr_sp = 0.25;
double Rmin_sp = RMIN0 * structf / (1.0 - shut);
double Pd_sp   = Pdia * FFR0 * ffr_sp;
double Pds_sp  = Psys * FFR0 * ffr_sp;
double v_e_sp  = cmd_vlayer(0, Rmin_sp / pow(1.0 - A_TONE * TONE_E, 4.0), fdia,
                            Pd_sp, Pds_sp, PV, lvedp_w, W_ENDO, W_EPI,
                            RCOMP_K, PHI_SYS);
double v_p_sp  = cmd_vlayer(1, Rmin_sp / pow(1.0 - A_TONE * TONE_P, 4.0), fdia,
                            Pd_sp, Pds_sp, PV, lvedp_w, W_ENDO, W_EPI,
                            RCOMP_K, PHI_SYS);

// duty-cycle weighted flows and deficits
double v_e = (1.0 - p_sp) * v_e_ns + p_sp * v_e_sp;
double v_p = (1.0 - p_sp) * v_p_ns + p_sp * v_p_sp;
double def_e_ns = CMD_POS(dem_e - v_e_ns * kO2 * E_MAX);
double def_e_sp = CMD_POS(dem_e - v_e_sp * kO2 * E_MAX);
double def_e = (1.0 - p_sp) * def_e_ns + p_sp * def_e_sp;
double def_p = CMD_POS(dem_p - v_p * kO2 * E_MAX);

// ===========================================================================
// 5.  THE FAST CONTROL LOOP (tau = 15 s; LSODA handles the stiffness)
// ===========================================================================
double tgt_e = TONE_REF + G_AUTO * (v_e_ns / vreq_e - 1.0) + myo + tp;
double tgt_p = TONE_REF + G_AUTO * (v_p_ns / vreq_p - 1.0) + myo + tp;
if (tgt_e < 0.0) tgt_e = 0.0;  if (tgt_e > TONE_MAX) tgt_e = TONE_MAX;
if (tgt_p < 0.0) tgt_p = 0.0;  if (tgt_p > TONE_MAX) tgt_p = TONE_MAX;
dxdt_TONE_E = (tgt_e - TONE_E) / TAU_TONE;
dxdt_TONE_P = (tgt_p - TONE_P) / TAU_TONE;

// ===========================================================================
// 6.  AFFERENT SIGNAL
// ===========================================================================
// bug B13/B16: the afferent drive has TWO sources - a true deficit and the
// adenosine arm that the controller offset drives on its own, the latter
// scaling with EXCESS demand so that it vanishes at rest and rises on
// exertion.  Only the second is blocked by aminophylline.  Without it the
// functional endotype, whose minimal resistance is normal and whose deficit is
// therefore zero at every workload, has no symptoms at all.
double ado = 1.0 + A_ADO_DF * def_e / (1.0 + def_e)
             + A_ADO_OF * AUTO_OFF * CMD_POS(mvo2 / MVO2_REF - 1.0);
double noci = K_NOC_DF * def_e + K_NOC_AD * (ado - 1.0) * (1.0 - f_A1);

dxdt_BI = (def_e - BI) / TAU_BURD;
dxdt_BN = (noci - BN) / TAU_BURD;

// ===========================================================================
// 7.  ENDOTHELIUM AND CONSTRICTORS
// ===========================================================================
// bug B21: these five were originally written with unit rate constants, i.e. a
// 1 h time constant used as a placeholder.  The regulated pools turn over on
// hours to days; stating the time constants explicitly also keeps an explicit
// integrator stable, which the placeholder did not.
double ros_t = KROS_BAS * (1.0 + 0.55 * (INFL - 1.0)) * 2.857
               * (1.0 - cmd_emax(C_ATO, EC_ATO, EM_ATO_RS)
                      - cmd_emax(M_RAM, EC_RAM, EM_RAM_RS)
                      - 0.16 * REHAB);
dxdt_ROS = (ros_t - ROS) / TAU_ROS;
double bh4_t = 1.0 / (1.0 + KBH4_OX * CMD_POS(ROS - 1.0));
dxdt_BH4 = (bh4_t - BH4) / TAU_BH4;
double no_t = BH4 * (1.0 + 0.10 * REHAB) / ((ROS > 1e-6) ? ROS : 1e-6);
dxdt_NO = (no_t - NO) / TAU_NO;

// bug B8: ETA blockade RAISES plasma ET-1, because it also removes
// ETA-mediated clearance.  PRIZE measured exactly that; a model that lowered
// ET-1 would have predicted a win.
double et_syn = KET_SYN * (1.0 + A_INFL_E * (INFL - 1.0) + A_GENE_E * GENO)
                / (1.0 + B_NO_ET * NO);
double et_t   = et_syn / (KET_CLR * (1.0 - 0.42 * e_zib_eta));
dxdt_ET1 = (et_t - ET1) / TAU_ET1;
double rock_t = ROCK_D * (1.0 - e_ato_rk - e_fas_rk);
dxdt_ROCK = (rock_t - ROCK) / TAU_ROCK;

// ===========================================================================
// 8.  STRUCTURE (months to years)
// ===========================================================================
double ml_t = (RMIN_F - 1.0) * (1.0 - cmd_emax(M_RAM, EC_RAM, EM_RAM_ML)
                                    - cmd_emax(C_ATO, EC_ATO, EM_ATO_ML));
dxdt_ML = (ml_t - ML) / TAU_ML;
double capd_t = 1.0 / (1.0 + 0.45 * (RMIN_F - 1.0)) * (1.0 + EM_RHB_CD * REHAB);
dxdt_CAPD = (capd_t - CAPD) / TAU_CAPD;
double pvf_t = 0.6 * (RMIN_F - 1.0) * (1.0 - cmd_emax(M_RAM, EC_RAM, EM_RAM_PV));
dxdt_PVF = (pvf_t - PVF) / TAU_PVF;
dxdt_ICF = (0.25 * (INFL - 1.0) - ICF) / TAU_ICF;
dxdt_LVH = (1.0 + 0.10 * (INFL - 1.0) - LVH) / TAU_LVH;

// ===========================================================================
// 9.  VOLUME, FILLING PRESSURE, LATE SODIUM CURRENT
// ===========================================================================
double vol_t = 1.0 + ((NOFLUID > 0.5) ? 0.0
                      : EM_ZIB_VL * e_zib_eta / EM_ZIB_ET);
dxdt_VOL = (vol_t - VOL) / TAU_VOL;

// bug B19: this loop (ischaemia -> late I_Na -> diastolic Ca -> LVEDP ->
// (Pd - LVEDP) -> ischaemia) had linear gains and no ceiling and ran away.
double lvedp_t = LVEDP0 + K_VOL_LP * (VOL - 1.0) + K_CAD_LP * (CAD - 1.0)
                 + K_ICF_LP * ICF;
if (lvedp_t > LVEDP_MX) lvedp_t = LVEDP_MX;
dxdt_LVEDP = (lvedp_t - LVEDP) / TAU_LVDP;

double late = (1.0 + K_ISCH_N * BI / (K_ISCH_K + BI))
              * (1.0 - cmd_emax(C_RAN, EC_RAN, EM_RAN_LN));
dxdt_NAI = (late - NAI) / TAU_NAI;
dxdt_CAD = (NAI - CAD) / TAU_CAD;

// ===========================================================================
// 10.  SYMPTOMS
// ===========================================================================
// bug B23: central sensitisation needs a primary drive of its own and angina a
// component sensitisation can generate with no afferent input, or the model
// cannot represent the fifth of every ANOCA cohort with normal coronary
// function and real pain.  bug B5: both feedback states are bounded, or the
// angina -> sympathetic -> angina loop diverges for any realistic gain.
double ang_t = K_ANG * BN * (1.0 + K_ANG_GN * SENS) + K_ANG_CN * SENS * SENS;
dxdt_ANG = (ang_t - ANG) / TAU_ANG;

double sens_off = 1.0 + cmd_emax(C_IMI, EC_IMI, EM_IMI_SN)
                  + 0.30 * CBT + 0.18 * REHAB;
dxdt_SENS = (K_SENS_ON * ANG * (1.0 - SENS) + SENS_BAS * (1.0 - SENS)
             - sens_off * SENS) / TAU_SENS;
dxdt_SYMP = (ANG / (K_SYMP_K + ANG) - SYMP) / TAU_SYMP;

double saq_ss = 100.0 - SAQ_MAX * ANG / (SAQ_K + ANG) - SAQ_SENS * SENS;
dxdt_SAQ = (saq_ss - SAQ) / TAU_SAQ;

// ===========================================================================
// 11.  BIOMARKERS AND HAZARDS
// ===========================================================================
dxdt_TNI = KTNI_OFF * (3.0 + TNI_MAX * BI / (2.5 + BI) - TNI);   // bug B6
dxdt_BNP = KBNP_ON * (LVEDP - LVEDP0) - KBNP_OFF * (BNP - 60.0);

// resting CFR for the hazard model, computed from the passive hyperaemic tone
double tp_h  = tp * (1.0 - F_ADO_REV);
double Rh    = Rmin_ns / pow(1.0 - A_TONE * tp_h, 4.0);
double vh_e  = cmd_vlayer(0, Rh, fdia, Pd_ns, Pds_ns, PV, LVEDP, W_ENDO, W_EPI,
                          RCOMP_K, PHI_SYS);
double vh_p  = cmd_vlayer(1, Rh, fdia, Pd_ns, Pds_ns, PV, LVEDP, W_ENDO, W_EPI,
                          RCOMP_K, PHI_SYS);
double vmean = 0.5 * (v_e + v_p);
double cfr_now = (vmean > 1e-9) ? 0.5 * (vh_e + vh_p) / vmean : 0.0;

dxdt_CHMORT = H0_MORT / 8760.0 * exp(LN116 * CMD_POS(2.5 - cfr_now) / 0.1);
dxdt_CHMACE = H0_HOSP / 8760.0 * exp(LN108 * CMD_POS(2.5 - cfr_now) / 0.1);
dxdt_CHHOSP = (H0_HOSP + K_HOSP_A * ANG) / 8760.0;
dxdt_AUCDEF = def_e;
dxdt_AUCNOC = noci;

// ===========================================================================
// 12.  PHARMACOKINETICS
// ===========================================================================
dxdt_A_RAN = -KA_RAN * A_RAN;
dxdt_C_RAN =  KA_RAN * A_RAN / V_RAN - KE_RAN * C_RAN;
dxdt_A_IVA = -KA_IVA * A_IVA;
dxdt_C_IVA =  KA_IVA * A_IVA / V_IVA - KE_IVA * C_IVA;
dxdt_M_IVA =  KF_IVA * C_IVA - KM_IVA * M_IVA;
dxdt_A_BIS = -KA_BIS * A_BIS;
dxdt_C_BIS =  KA_BIS * A_BIS / V_BIS - KE_BIS * C_BIS;
dxdt_A_NEB = -KA_NEB * A_NEB;
dxdt_C_NEB =  KA_NEB * A_NEB / V_NEB - KE_NEB * C_NEB;
dxdt_A_AML = -KA_AML * A_AML;
dxdt_C_AML =  KA_AML * A_AML / V_AML - KE_AML * C_AML;
dxdt_A_DIL = -KA_DIL * A_DIL;
dxdt_C_DIL =  KA_DIL * A_DIL / V_DIL - KE_DIL * C_DIL;
dxdt_A_NIC = -KA_NIC * A_NIC;
dxdt_C_NIC =  KA_NIC * A_NIC / V_NIC - KE_NIC * C_NIC;
dxdt_A_ZIB = -KA_ZIB * A_ZIB;
dxdt_C_ZIB =  KA_ZIB * A_ZIB / V_ZIB - KE_ZIB * C_ZIB;
dxdt_A_RAM = -KA_RAM * A_RAM;
dxdt_C_RAM =  KA_RAM * A_RAM / V_RAM - KE_RAM * C_RAM;
dxdt_M_RAM =  KF_RAM * C_RAM - KM_RAM * M_RAM;
dxdt_A_ATO = -KA_ATO * A_ATO;
dxdt_C_ATO =  KA_ATO * A_ATO / V_ATO - KE_ATO * C_ATO;
dxdt_A_FAS = -KA_FAS * A_FAS;
dxdt_C_FAS =  KA_FAS * A_FAS / V_FAS - KE_FAS * C_FAS;
dxdt_M_FAS =  KF_FAS * C_FAS - KM_FAS * M_FAS;
dxdt_A_SIL = -KA_SIL * A_SIL;
dxdt_C_SIL =  KA_SIL * A_SIL / V_SIL - KE_SIL * C_SIL;
dxdt_A_AMI = -KA_AMI * A_AMI;
dxdt_C_AMI =  KA_AMI * A_AMI / V_AMI - KE_AMI * C_AMI;
dxdt_A_IMI = -KA_IMI * A_IMI;
dxdt_C_IMI =  KA_IMI * A_IMI / V_IMI - KE_IMI * C_IMI;

$TABLE
// ===========================================================================
//  The invasive coronary function test, as it is actually performed.
//  Everything here is evaluated at RESTING haemodynamics, because that is the
//  condition in the catheter laboratory.
// ===========================================================================
double e_zib_eta_t = cmd_emax(C_ZIB, EC_ZIB, EM_ZIB_ET);
double e_ami_a2_t  = cmd_emax(C_AMI, EC_AMI, EM_AMI_A2);
double f_TONE_t = cmd_emax(C_AML, EC_AML, EM_AML_TN)
                + cmd_emax(C_DIL, EC_DIL, EM_DIL_TN)
                + cmd_emax(C_NIC, EC_NIC, EM_NIC_TN)
                + cmd_emax(C_SIL, EC_SIL, EM_SIL_TN)
                - 0.20 * e_ami_a2_t;
double f_MAP_t = cmd_emax(C_AML, EC_AML, EM_AML_MP)
               + cmd_emax(C_DIL, EC_DIL, EM_DIL_MP)
               + cmd_emax(C_NIC, EC_NIC, EM_NIC_MP)
               + cmd_emax(C_ZIB, EC_ZIB, EM_ZIB_MP)
               + cmd_emax(M_RAM, EC_RAM, EM_RAM_MP)
               + cmd_emax(C_SIL, EC_SIL, EM_SIL_MP);
if (NOBPDROP > 0.5) f_MAP_t = 0.0;
if (f_MAP_t > 0.25) f_MAP_t = 0.25;
double f_HR_t = cmd_emax(C_IVA + 0.7 * M_IVA, EC_IVA, EM_IVA_HR)
              + cmd_emax(C_BIS, EC_BIS, EM_BIS_HR)
              + cmd_emax(C_NEB, EC_NEB, EM_NEB_HR)
              + cmd_emax(C_DIL, EC_DIL, EM_DIL_HR)
              + cmd_emax(C_RAN, EC_RAN, EM_RAN_HR);
double f_CTR_t = cmd_emax(C_BIS, EC_BIS, EM_BIS_CT)
               + cmd_emax(C_NEB, EC_NEB, EM_NEB_CT)
               + cmd_emax(C_DIL, EC_DIL, EM_DIL_CT);
double d_NO_t = cmd_emax(C_NEB, EC_NEB, EM_NEB_NO)
              + cmd_emax(C_NIC, EC_NIC, EM_NIC_NO)
              + cmd_emax(M_RAM, EC_RAM, EM_RAM_NO);
double f_A1_t = cmd_emax(C_AMI, EC_AMI, EM_AMI_A1);
double e_fas_t = cmd_emax(M_FAS, EC_FAS, EM_FAS_RK);
double anti_sp_t = cmd_emax(C_AML, EC_AML, EM_AML_TN) / EM_AML_TN
                 + cmd_emax(C_DIL, EC_DIL, EM_DIL_TN) / EM_DIL_TN
                 + e_fas_t / EM_FAS_RK;
if (anti_sp_t > 1.0) anti_sp_t = 1.0;

double HRr  = HR0 * (1.0 - f_HR_t) * (1.0 + K_SYMP_HR * SYMP);
if (HRr < 38.0) HRr = 38.0;
double SBPr = SBP0 * (1.0 - f_MAP_t) * (1.0 + 0.35 * K_SYMP_HR * SYMP);
double DBPr = DBP0 * (1.0 - f_MAP_t);
double Pdiar = DBPr + FR_DIA_P * (SBPr - DBPr);
double Psysr = DBPr + FR_SYS_P * (SBPr - DBPr);
double tcycr = 60.0 / HRr;
double tsysr = cmd_tsys(HRr);
double fdiar = (tcycr - tsysr) / tcycr;
double Pdr   = Pdiar * FFR0;
double Pdsr  = Psysr * FFR0;

double capd_t2 = (CAPD < 0.35) ? 0.35 : CAPD;
double Rmin_t = RMIN0 * (1.0 + W_ML * ML + W_PVF * PVF) / capd_t2;

double sigr = SBPr / SIGMA_RF;
double ctrr = 1.0 * (1.0 - f_CTR_t);
double tenr = 1.0 + K_TEN_CD * (CAD - 1.0);
double mv_r = (MVO2_B0 + K_TTI * HRr * sigr + K_CTR * ctrr + K_DIAST * tenr)
              * ((TMZ > 0.5) ? (1.0 - 0.055) : 1.0);
double kO2t = KCONV * 100.0 * CAO2;

// resting velocity: the loop is at its fixed point, so the tone states are the
// answer and no iteration is needed
double Rar_e = Rmin_t / pow(1.0 - A_TONE * TONE_E, 4.0);
double Rar_p = Rmin_t / pow(1.0 - A_TONE * TONE_P, 4.0);
double vr_e = cmd_vlayer(0, Rar_e, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO, W_EPI,
                         RCOMP_K, PHI_SYS);
double vr_p = cmd_vlayer(1, Rar_p, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO, W_EPI,
                         RCOMP_K, PHI_SYS);
double V_REST = 0.5 * (vr_e + vr_p);

// hyperaemic velocity: adenosine abolishes the metabolic and myogenic terms
// and F_ADO_REV of the constrictor tone
double et_eff_t = (ET1 - 1.0) * (F_ETA * (1.0 - e_zib_eta_t) + F_ETB2);
double tp_t = K_ET_TN * CMD_POS(et_eff_t) + K_ROCK_TN * CMD_POS(ROCK - 1.0)
            + K_A1_TN * SYMP * (1.0 - f_A1_t)
            - K_NO_TN * (NO - 1.0 + d_NO_t) - f_TONE_t;
if (tp_t < 0.0) tp_t = 0.0;
if (tp_t > TONE_MAX) tp_t = TONE_MAX;
double tone_h = tp_t * (1.0 - F_ADO_REV);
double Rar_h  = Rmin_t / pow(1.0 - A_TONE * tone_h, 4.0);
double vh_e_t = cmd_vlayer(0, Rar_h, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO,
                           W_EPI, RCOMP_K, PHI_SYS);
double vh_p_t = cmd_vlayer(1, Rar_h, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO,
                           W_EPI, RCOMP_K, PHI_SYS);
double V_HYP = 0.5 * (vh_e_t + vh_p_t);

// acetylcholine provocation: NO-mediated dilatation in a healthy endothelium,
// direct M3 constriction plus spasm when it has failed
double ach_t = -0.30 * ((NO < 1.0) ? NO : 1.0)
             + 0.55 * SPASM * (1.0 - anti_sp_t) * (1.0 + 0.5 * CMD_POS(ROCK - 1.0))
             + 0.22 * CMD_POS(1.0 - NO);
double tone_a = tp_t + ach_t + TONE_REF;
if (tone_a < 0.0) tone_a = 0.0;
if (tone_a > TONE_MAX) tone_a = TONE_MAX;
double Rar_a = Rmin_t / pow(1.0 - A_TONE * tone_a, 4.0);
double V_ACH = 0.5 * (cmd_vlayer(0, Rar_a, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO,
                                 W_EPI, RCOMP_K, PHI_SYS)
                    + cmd_vlayer(1, Rar_a, fdiar, Pdr, Pdsr, PV, LVEDP, W_ENDO,
                                 W_EPI, RCOMP_K, PHI_SYS));

double CFR    = (V_REST > 1e-9) ? V_HYP / V_REST : 0.0;
double MR_RST = (V_REST > 1e-9) ? Pdr / V_REST : 0.0;
double MR_HYP = (V_HYP  > 1e-9) ? Pdr / V_HYP  : 0.0;
double IMR    = MR_HYP * 8.35;         // Pd x Tmn, scaled so normal reads 18 U
double MRR    = CFR / FFR0;            // Pa_rest = Pa_hyper here, see note
double ACHFR  = (V_REST > 1e-9) ? V_ACH / V_REST : 0.0;
double ENDOEPI = (vh_p_t > 1e-9) ? vh_e_t / vh_p_t : 0.0;
double MBF_RST = V_REST * KCONV;
double MBF_HYP = V_HYP * KCONV;
double DPTI = CMD_POS(Pdr - LVEDP) * tcycr * fdiar;
double SPTI = Psysr * tsysr;
double SEVR = (SPTI > 0.0) ? DPTI / SPTI : 0.0;
double E_RST = (V_REST > 1e-9) ? mv_r / (V_REST * kO2t) : 0.0;
double HR_RST = HRr;
double ANGINA = (BN >= NOCI_THR / (1.0 + 0.55 * SENS)) ? 1.0 : 0.0;

$CAPTURE @annotated
CFR     : Coronary flow reserve (-)
MR_RST  : Resting microvascular resistance (mmHg/(cm/s))
MR_HYP  : Hyperaemic microvascular resistance (mmHg/(cm/s))
IMR     : Index of microcirculatory resistance (U)
MRR     : Microvascular resistance reserve (-)
ACHFR   : Acetylcholine flow reserve (-)
ENDOEPI : Hyperaemic subendocardial/subepicardial perfusion ratio (-)
MBF_RST : Resting myocardial blood flow (mL/min/g)
MBF_HYP : Hyperaemic myocardial blood flow (mL/min/g)
DPTI    : Diastolic pressure-time index (mmHg s)
SPTI    : Systolic tension-time index (mmHg s)
SEVR    : Subendocardial viability ratio (-)
E_RST   : Resting myocardial oxygen extraction (-)
HR_RST  : Resting heart rate (bpm)
ANGINA  : Afferent drive above the anginal threshold (0/1)
'

cmd_mod <- mcode("cmd_qsp", cmd_code)

## =============================================================================
##  ENDOTYPE PRESETS
## =============================================================================
##  These are the archetypes calibrated in `cmd_reference_model.py`.  The
##  functional endotype's structural floor is fixed at 1.000 BY DEFINITION —
##  a normal hyperaemic resistance is what the label means — so the only thing
##  that distinguishes it from a control is the controller offset.
## =============================================================================
cmd_endotypes <- list(
  control     = list(RMIN_F = 1.000, AUTO_OFF = 0.000, INFL = 1.0, GENO = 0,
                     ROCK_D = 1.00, KSBP = 16.0, SPASM = 0.00, SENS_BAS = 0.05),
  functional  = list(RMIN_F = 1.000, AUTO_OFF = 0.9165, INFL = 1.4, GENO = 1,
                     ROCK_D = 1.30, KSBP = 19.0, SPASM = 0.10, SENS_BAS = 0.75),
  structural  = list(RMIN_F = 1.4306, AUTO_OFF = 0.05, INFL = 2.2, GENO = 0,
                     ROCK_D = 1.50, KSBP = 37.0, SPASM = 0.05, SENS_BAS = 0.65),
  vasospastic = list(RMIN_F = 1.080, AUTO_OFF = 0.12, INFL = 1.9, GENO = 1,
                     ROCK_D = 3.10, KSBP = 19.0, SPASM = 0.85, SENS_BAS = 0.85),
  noncardiac  = list(RMIN_F = 1.000, AUTO_OFF = 0.000, INFL = 1.1, GENO = 0,
                     ROCK_D = 1.00, KSBP = 16.0, SPASM = 0.00, SENS_BAS = 1.60)
)

## =============================================================================
##  DOSING HELPERS
## =============================================================================
## Events are built as plain data frames with a NUMERIC cmt so that
## multi-drug regimens combine with rbind() and pass straight to data_set()
## without depending on a c() method for ev or on character-cmt resolution.
cmd_cmt <- function(nm) {
  i <- match(nm, names(init(cmd_mod)@data))
  if (is.na(i)) stop("unknown compartment: ", nm)
  i
}

cmd_ev <- function(cmt, amt, ii, days, start = 0) {
  data.frame(ID = 1, time = start, amt = amt, cmt = cmd_cmt(cmt), ii = ii,
             addl = max(0, floor(days * 24 / ii) - 1), evid = 1)
}

#' Build an event object for a named regimen.
#' Doses are the ones used in the trials the model is calibrated against.
cmd_regimen <- function(name, days = 168) {
  switch(name,
    untreated     = NULL,
    ranolazine    = cmd_ev("A_RAN", 1000, 12, days),
    ivabradine    = cmd_ev("A_IVA",  7.5, 12, days),
    bisoprolol    = cmd_ev("A_BIS",  5.0, 24, days),
    nebivolol     = cmd_ev("A_NEB",  5.0, 24, days),
    amlodipine    = cmd_ev("A_AML", 10.0, 24, days),
    diltiazem     = cmd_ev("A_DIL",  180, 24, days),
    nicorandil    = cmd_ev("A_NIC", 20.0,  8, days),
    zibotentan    = cmd_ev("A_ZIB", 10.0, 24, days),
    fasudil       = cmd_ev("A_FAS", 80.0,  8, days),
    sildenafil    = cmd_ev("A_SIL", 50.0,  8, days),
    aminophylline = cmd_ev("A_AMI",  225, 12, days),
    imipramine    = cmd_ev("A_IMI", 50.0, 24, days),
    imt_warrior   = rbind(cmd_ev("A_RAM", 10.0, 24, days),
                          cmd_ev("A_ATO", 80.0, 24, days)),
    iva_ran       = rbind(cmd_ev("A_IVA",  7.5, 12, days),
                          cmd_ev("A_RAN", 1000, 12, days)),
    cormica_func  = rbind(cmd_ev("A_IVA",  7.5, 12, days),
                          cmd_ev("A_RAM", 10.0, 24, days),
                          cmd_ev("A_ATO", 80.0, 24, days)),
    cormica_struct = rbind(cmd_ev("A_NEB", 5.0, 24, days),
                           cmd_ev("A_RAM", 10.0, 24, days),
                           cmd_ev("A_ATO", 80.0, 24, days)),
    cormica_vaso  = rbind(cmd_ev("A_AML", 10.0, 24, days),
                          cmd_ev("A_ATO", 80.0, 24, days)),
    cormica_noncard = cmd_ev("A_IMI", 50.0, 24, days),
    stop("unknown regimen: ", name))
}

#' Extra (non-dose) parameter switches that some regimens need
cmd_switches <- list(
  trimetazidine   = list(TMZ = 1),
  rehab           = list(REHAB = 1),
  cormica_noncard = list(CBT = 1),
  zibotentan_cf   = list(NOBPDROP = 1, NOFLUID = 1),
  zib_nobp        = list(NOBPDROP = 1),
  zib_nofluid     = list(NOFLUID = 1)
)

#' Simulate one endotype on one regimen.
#'
#' @param endotype one of names(cmd_endotypes)
#' @param regimen  one of the regimens in cmd_regimen(), or "trimetazidine",
#'                 "rehab", "zibotentan_cf", "zib_nobp", "zib_nofluid"
#' @param days     follow-up (days)
#' @param wl_fix   hold workload at this value instead of the diurnal profile
cmd_simulate <- function(endotype = "functional", regimen = "untreated",
                         days = 168, delta = 1, wl_fix = 0, extra = list()) {
  stopifnot(endotype %in% names(cmd_endotypes))
  par <- cmd_endotypes[[endotype]]
  par$WL_FIX <- wl_fix
  base_reg <- regimen
  if (regimen %in% names(cmd_switches)) {
    par <- c(par, cmd_switches[[regimen]])
    base_reg <- switch(regimen,
                       trimetazidine = "untreated",
                       rehab = "untreated",
                       zibotentan_cf = "zibotentan",
                       zib_nobp = "zibotentan",
                       zib_nofluid = "zibotentan",
                       regimen)
  }
  par <- c(par, extra)
  m <- param(cmd_mod, par)
  dose <- cmd_regimen(base_reg, days = days)
  out <- m
  if (!is.null(dose)) out <- data_set(out, dose)
  sim <- out %>%
    mrgsim(end = days * 24, delta = delta, hmax = 0.05,
           atol = 1e-8, rtol = 1e-6) %>%
    as.data.frame()
  sim$day <- sim$time / 24
  sim$endotype <- endotype
  sim$regimen <- regimen
  sim
}

## =============================================================================
##  THE INVASIVE COronARY FUNCTION TEST, ON DEMAND
## =============================================================================
##  CFR, both microvascular resistances, IMR, MRR, the acetylcholine flow
##  reserve and the transmural ratio are all captured continuously, so a
##  "catheter laboratory visit" is just a row of the output.  Endotype
##  classification follows Rahman 2019 exactly.
## =============================================================================
cmd_classify <- function(sim, day = NULL) {
  row <- if (is.null(day)) sim[nrow(sim), ] else
    sim[which.min(abs(sim$day - day)), ]
  lab <- if (row$CFR >= 2.5) {
    if (row$ACHFR < 0.9) "vasospastic (normal CFR, positive ACh)" else
      "no coronary microvascular dysfunction"
  } else if (row$MR_HYP < 2.5) {
    "CMD, functional endotype (normal minimal resistance)"
  } else {
    "CMD, structural endotype (raised minimal resistance)"
  }
  data.frame(day = row$day, CFR = row$CFR, MR_RST = row$MR_RST,
             MR_HYP = row$MR_HYP, IMR = row$IMR, MRR = row$MRR,
             ACHFR = row$ACHFR, ENDOEPI = row$ENDOEPI,
             SAQ = row$SAQ, ANG = row$ANG, endotype_call = lab,
             stringsAsFactors = FALSE)
}

## =============================================================================
##  BRUCE TREADMILL TEST
## =============================================================================
##  Workload is stepped through the Bruce stages and the test stops when the
##  afferent drive crosses the anginal threshold.  This is the endpoint PRIZE
##  used, and section VI of the reference output shows why it was the wrong
##  instrument for that trial.
## =============================================================================
cmd_bruce_wl <- function(t_s) {
  stg <- rbind(c(0, 180, 1.0, 4.6), c(180, 360, 4.6, 7.0),
               c(360, 540, 7.0, 10.2), c(540, 720, 10.2, 12.9),
               c(720, 900, 12.9, 15.0))
  mets <- 15.0
  for (i in seq_len(nrow(stg))) {
    if (t_s < stg[i, 2]) {
      mets <- stg[i, 3] + (t_s - stg[i, 1]) / (stg[i, 2] - stg[i, 1]) *
        (stg[i, 4] - stg[i, 3])
      break
    }
  }
  1.0 + (mets - 1.0) / 4.35
}

#' Bruce duration (s) for a patient in a given state.
#' The state is taken from a completed cmd_simulate() run, so the treadmill
#' test is performed on the treated patient rather than on a fresh one.
cmd_bruce <- function(endotype, regimen, days = 168, grid = 20) {
  base <- cmd_simulate(endotype, regimen, days = days, delta = 24)
  last <- base[nrow(base), ]
  init_state <- as.list(last[names(last) %in% names(init(cmd_mod)@data)])
  par <- cmd_endotypes[[endotype]]
  if (regimen %in% names(cmd_switches)) par <- c(par, cmd_switches[[regimen]])
  thr <- 3.8575 / (1.0 + 0.55 * last$SENS)
  ts <- seq(5, 900, length.out = grid)
  drive <- vapply(ts, function(t_s) {
    m <- param(cmd_mod, c(par, list(WL_FIX = cmd_bruce_wl(t_s))))
    m <- init(m, init_state)
    s <- mrgsim(m, end = 0.25, delta = 0.25, hmax = 0.01) %>% as.data.frame()
    s$BN[nrow(s)]
  }, numeric(1))
  idx <- which(drive >= thr)
  if (!length(idx)) return(900)
  ts[min(idx)]
}

## =============================================================================
##  TWENTY-FOUR SCENARIOS
## =============================================================================
cmd_scenarios <- data.frame(
  id = 1:24,
  label = c("untreated functional", "untreated structural",
            "untreated vasospastic", "non-cardiac chest pain",
            "ivabradine / functional", "ivabradine / structural",
            "bisoprolol / functional", "nebivolol / functional",
            "ranolazine / functional", "ranolazine / structural",
            "ranolazine / control", "amlodipine / vasospastic",
            "diltiazem / vasospastic", "nicorandil / structural",
            "zibotentan / functional", "zibotentan counterfactual / functional",
            "IMT (WARRIOR) / structural", "fasudil / vasospastic",
            "aminophylline / functional", "imipramine / non-cardiac",
            "exercise rehabilitation / structural",
            "ivabradine + ranolazine / functional",
            "sildenafil / structural", "CorMicA-directed / functional"),
  endotype = c("functional", "structural", "vasospastic", "noncardiac",
               "functional", "structural", "functional", "functional",
               "functional", "structural", "control", "vasospastic",
               "vasospastic", "structural", "functional", "functional",
               "structural", "vasospastic", "functional", "noncardiac",
               "structural", "functional", "structural", "functional"),
  regimen = c("untreated", "untreated", "untreated", "untreated",
              "ivabradine", "ivabradine", "bisoprolol", "nebivolol",
              "ranolazine", "ranolazine", "ranolazine", "amlodipine",
              "diltiazem", "nicorandil", "zibotentan", "zibotentan_cf",
              "imt_warrior", "fasudil", "aminophylline", "imipramine",
              "rehab", "iva_ran", "sildenafil", "cormica_func"),
  stringsAsFactors = FALSE)

#' Run all 24 scenarios and return the 24-week end-points.
cmd_run_all <- function(days = 168) {
  do.call(rbind, lapply(seq_len(nrow(cmd_scenarios)), function(i) {
    s <- cmd_scenarios[i, ]
    sim <- cmd_simulate(s$endotype, s$regimen, days = days, delta = 24)
    last <- sim[nrow(sim), ]
    data.frame(id = s$id, label = s$label, endotype = s$endotype,
               regimen = s$regimen, CFR = last$CFR, MR_HYP = last$MR_HYP,
               IMR = last$IMR, ENDOEPI = last$ENDOEPI, burden = last$BI,
               noci = last$BN, ANG = last$ANG, SAQ = last$SAQ,
               LVEDP = last$LVEDP, TNI = last$TNI, BNP = last$BNP,
               stringsAsFactors = FALSE)
  }))
}

## =============================================================================
##  THE FIVE TRIAL REPRODUCTIONS
## =============================================================================

#' CorMicA (PMID 30266608): stratified therapy vs unguided care.
#' Observed SAQ summary +11.7 U at 6 months (95% CI 5.0-18.4).
cmd_cormica <- function(days = 182) {
  mix <- c(functional = 0.30, structural = 0.19,
           vasospastic = 0.29, noncardiac = 0.22)
  strat <- c(functional = "cormica_func", structural = "cormica_struct",
             vasospastic = "cormica_vaso", noncardiac = "cormica_noncard")
  out <- do.call(rbind, lapply(names(mix), function(k) {
    u <- cmd_simulate(k, "bisoprolol", days = days, delta = 24)
    s <- cmd_simulate(k, strat[[k]], days = days, delta = 24)
    data.frame(stratum = k, weight = mix[[k]],
               SAQ_standard = u$SAQ[nrow(u)], SAQ_stratified = s$SAQ[nrow(s)],
               stringsAsFactors = FALSE)
  }))
  out$dSAQ <- out$SAQ_stratified - out$SAQ_standard
  attr(out, "population_dSAQ") <- sum(out$weight * out$dSAQ)
  out
}

#' RWISE (PMID 26614823): ranolazine, null overall, positive at CFR < 2.5.
#' The dilution is arithmetic, not pharmacological.
cmd_rwise <- function(days = 84) {
  mix <- c(control = 0.47, functional = 0.33, structural = 0.20)
  out <- do.call(rbind, lapply(names(mix), function(k) {
    p <- cmd_simulate(k, "untreated", days = days, delta = 24)
    r <- cmd_simulate(k, "ranolazine", days = days, delta = 24)
    data.frame(stratum = k, weight = mix[[k]],
               SAQ_placebo = p$SAQ[nrow(p)], SAQ_ranolazine = r$SAQ[nrow(r)],
               dCFR = r$CFR[nrow(r)] - p$CFR[nrow(p)], stringsAsFactors = FALSE)
  }))
  out$dSAQ <- out$SAQ_ranolazine - out$SAQ_placebo
  w_cmd <- mix[["functional"]] + mix[["structural"]]
  attr(out, "dSAQ_CMD_stratum") <-
    sum(out$weight[out$stratum != "control"] * out$dSAQ[out$stratum != "control"]) / w_cmd
  attr(out, "dSAQ_whole_cohort") <- sum(out$weight * out$dSAQ)
  out
}

#' PRIZE (PMID 39217504): zibotentan, -4.26 s of Bruce time (95% CI -19.6 to
#' +11.1).  The decomposition below shows the microvascular gain is zero, the
#' blood-pressure fall is a small net BENEFIT, and fluid retention accounts for
#' the whole null: the trial's dominant adverse effect and its failure on the
#' primary endpoint are the same event measured twice.
cmd_prize <- function(days = 84) {
  mix <- c(functional = 0.33, structural = 0.20,
           vasospastic = 0.29, noncardiac = 0.18)
  arms <- c("zibotentan", "zib_nobp", "zib_nofluid", "zibotentan_cf")
  do.call(rbind, lapply(names(mix), function(k) {
    b <- cmd_bruce(k, "untreated", days = days)
    data.frame(endotype = k, weight = mix[[k]], arm = c("placebo", arms),
               bruce_s = c(b, vapply(arms, function(a) cmd_bruce(k, a, days = days),
                                     numeric(1))),
               stringsAsFactors = FALSE)
  }))
}

#' WARRIOR (PMID 41932694): intensive medical therapy, MACE HR 1.13 (0.94-1.37);
#' the trial's own contamination-adjusted estimate was 0.74.  Two problems are
#' separable here: background contamination, and a structural target scored on
#' a symptom endpoint.
cmd_warrior <- function(days = 912) {
  u <- cmd_simulate("structural", "untreated", days = days, delta = 168)
  t <- cmd_simulate("structural", "imt_warrior", days = days, delta = 168)
  hr_true <- t$CHHOSP[nrow(t)] / u$CHHOSP[nrow(u)]
  contam <- c(0, 0.2, 0.4, 0.6, 0.75, 0.85, 0.95)
  list(hr_true = hr_true,
       dCFR = (t$CFR[nrow(t)] - t$CFR[1]) - (u$CFR[nrow(u)] - u$CFR[1]),
       dSAQ = t$SAQ[nrow(t)] - u$SAQ[nrow(u)],
       observable = data.frame(contamination = contam,
                               HR_observed = 1 + (hr_true - 1) * (1 - contam)))
}

#' Rahman 2019 (PMID 31707835): the endotype identity.  The same CFR, two
#' different diseases, and the absolute maximal flows differ by 57%.
cmd_endotype_table <- function() {
  do.call(rbind, lapply(names(cmd_endotypes), function(k) {
    s <- cmd_simulate(k, "untreated", days = 1, delta = 1)
    r <- s[nrow(s), ]
    data.frame(endotype = k, CFR = r$CFR, MR_rest = r$MR_RST,
               MR_hyp = r$MR_HYP, IMR = r$IMR,
               MBF_rest = r$MBF_RST, MBF_hyp = r$MBF_HYP,
               endo_epi_hyp = r$ENDOEPI, E_rest = r$E_RST,
               stringsAsFactors = FALSE)
  }))
}

## =============================================================================
##  HEART RATE ENTERS THE OXYGEN BALANCE TWICE
## =============================================================================
##  Sweeping resting heart rate at a FIXED workload separates the two routes.
##  Holding demand at its HR = 68 value isolates the diastolic time window;
##  43% of the fall in subendocardial deficit from 68 -> 55 bpm is time and 57%
##  is demand (reference output section III).
## =============================================================================
cmd_hr_sweep <- function(endotype = "structural", wl = 3.2,
                         hrs = c(50, 55, 60, 68, 75, 85, 95, 110)) {
  do.call(rbind, lapply(hrs, function(hr) {
    par <- c(cmd_endotypes[[endotype]],
             list(WL_FIX = wl, HR0 = hr, KHR_WL = 30.0))
    s <- param(cmd_mod, par) %>%
      mrgsim(end = 48, delta = 1, hmax = 0.05) %>% as.data.frame()
    r <- s[nrow(s), ]
    data.frame(HR = hr, DPTI = r$DPTI, SPTI = r$SPTI, SEVR = r$SEVR,
               deficit = r$BI, noci = r$BN, endo_epi = r$ENDOEPI,
               stringsAsFactors = FALSE)
  }))
}

## =============================================================================
##  WHAT THE HYPERAEMIC-MR CUT-OFF ACTUALLY SEPARATES
## =============================================================================
##  A hyperaemic MR of 2.5 is read as "structural", but two different things
##  raise it: inward remodelling and capillary rarefaction, which fix the
##  minimal radius, and constrictor tone that adenosine does not fully reverse
##  (87.6% of it is reversed here, calibrated).  The first is not reversible on
##  any useful timescale; the second is reversible in weeks.  The label cannot
##  tell them apart — a repeat measurement on a Rho-kinase inhibitor can.
##
##  PREDICTION: repeat the hyperaemic measurement after acute Rho-kinase or ETA
##  blockade before telling a patient their microvascular disease is
##  structural.  It is a one-visit test and it changes the therapeutic class.
## =============================================================================
cmd_reclassify <- function(endotype = "structural", days = 21) {
  a <- cmd_classify(cmd_simulate(endotype, "untreated", days = days, delta = 24))
  b <- cmd_classify(cmd_simulate(endotype, "fasudil", days = days, delta = 24))
  rbind(cbind(arm = "baseline", a), cbind(arm = "Rho-kinase inhibitor", b))
}

## =============================================================================
##  VIRTUAL POPULATION
## =============================================================================
##  Risk factors are sampled first and the structural floor and controller
##  offset are drawn CONDITIONALLY on them; independent draws produce patients
##  with severe remodelling and no risk factor at all, and miss the observed
##  62/38 functional/structural split among CFR < 2.5 (reference bug B10).
## =============================================================================
cmd_population <- function(n = 200, seed = 20260804, days = 84) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    htn <- runif(1) < 0.55
    dm  <- runif(1) < 0.22
    infl <- max(1, 1 + rnorm(1, 0.35, 0.30) + 0.25 * htn + 0.35 * dm)
    rmin_f <- 1 + max(0, rnorm(1, 0.02 + 0.04 * htn + 0.11 * dm, 0.04))
    auto <- max(0, rnorm(1, 0.26 + 0.10 * htn, 0.33))
    par <- list(RMIN_F = rmin_f, AUTO_OFF = auto, INFL = infl,
                GENO = as.numeric(runif(1) < 0.5),
                ROCK_D = 1 + 0.35 * (infl - 1),
                KSBP = if (rmin_f > 1.15) 37 else 19,
                SPASM = max(0, rnorm(1, 0.25, 0.28)),
                SENS_BAS = max(0, rnorm(1, 0.72, 0.45)))
    s <- param(cmd_mod, par) %>%
      mrgsim(end = days * 24, delta = 24, hmax = 0.05) %>% as.data.frame()
    r <- s[nrow(s), ]
    data.frame(id = i, htn = htn, dm = dm, RMIN_F = rmin_f, AUTO_OFF = auto,
               CFR = r$CFR, MR_rest = r$MR_RST, MR_hyp = r$MR_HYP, IMR = r$IMR,
               SAQ = r$SAQ, ANG = r$ANG, stringsAsFactors = FALSE)
  }))
}

## =============================================================================
##  QUICK LOOK
## =============================================================================
if (interactive()) {
  print(cmd_endotype_table())
  print(cmd_run_all())
  print(cmd_hr_sweep())
  cm <- cmd_cormica(); print(cm)
  cat("population dSAQ:", attr(cm, "population_dSAQ"), "U (observed +11.7)\n")
  rw <- cmd_rwise(); print(rw)
  cat("dSAQ CMD stratum:", attr(rw, "dSAQ_CMD_stratum"),
      " whole cohort:", attr(rw, "dSAQ_whole_cohort"), "\n")
  print(cmd_reclassify())

  sim <- cmd_simulate("functional", "cormica_func", days = 168)
  ggplot(sim, aes(day, SAQ)) + geom_line() +
    labs(title = "Endotype-directed therapy, functional CMD",
         x = "day", y = "SAQ summary score") + theme_bw()
}
