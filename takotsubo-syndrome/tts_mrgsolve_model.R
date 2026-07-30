## =============================================================================
##  tts_mrgsolve_model.R
##  Takotsubo syndrome (stress-induced / catecholamine-mediated cardiomyopathy)
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  64 ODE compartments.  Time unit = HOURS.  A 2-minute adrenaline half-life,
##  a 4-day receptor-dephosphorylation half-life and a 1-year recurrence hazard
##  have to live in the same system, so the unit is chosen for the fastest
##  state and the slow ones are integrated over 8760 h.
##
##  ---------------------------------------------------------------------------
##  THE FOUR STRUCTURAL COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (1) "APICAL BALLOONING" IS A THRESHOLD, NOT A LOCATION.
##
##      There are three myocardial segments (APEX / MID / BASE).  They differ
##      in EXACTLY TWO numbers, both of which are measured quantities and
##      neither of which mentions contraction:
##
##          RHO_x   total beta-AR density        1.40 / 1.15 / 1.00
##          FB2_x   beta2 fraction of that pool  0.42 / 0.32 / 0.24
##          (plus INN_x, sympathetic innervation 0.62 / 0.90 / 1.20, which is
##           the reciprocal gradient measured by 123I-MIBG)
##
##      A SINGLE agonist-occupancy-dependent Gs -> Gi switch function is then
##      applied identically to all three:
##
##          SIG_x  = ((occ2_x - THR_SW)_+)^2 / (... + KSW^2) * (1 + WPKA*PKA_x)
##          PHOS_x' = KPH*SIG_x*(1 - PHOS_x) - KDEPH*PHOS_x
##          GI_x   = b2_x * occ2_x * PHOS_x        <- DENSITY enters here
##          GS_x   = b1_x * occ1_x + b2_x * occ2_x * (1 - PHOS_x)
##          AC_x   = (AC_BAS + EMAX_AC*GS_x/(GS_x+KM_GS)) / (1 + GAM_GI*GI_x)
##
##      Nothing in this file says "the apex is akinetic".  The word APEX
##      appears only in the three parameter values above and in the output
##      names.  Apical akinesis, basal hypercontractility, the LVOT gradient
##      and the ballooning index are all READ OUT of one equation evaluated at
##      three receptor densities.
##
##      A consequence worth stating because it is not obvious: at REST the
##      three segments have almost identical cAMP (0.149 / 0.147 / 0.152 in
##      GS units) because the higher apical receptor density is offset by the
##      lower apical innervation.  The segments are indistinguishable until a
##      BLOOD-BORNE beta2 agonist arrives.  That is why this disease needs
##      ADRENALINE and not sympathetic nerve traffic, and it is why the
##      NE-only scenario (S2) must fail to produce ballooning.
##
##  (2) THE BENEFIT OF AN INOTROPE IS A PRODUCT, AND THE SECOND FACTOR CAN BE
##      NEGATIVE.
##
##          benefit = (contractile drive added)
##                    x (fraction of that drive delivered WITHOUT raising the
##                       phosphorylation that creates the Gi coupling)
##
##      Dobutamine, adrenaline, dopamine and milrinone all enter cAMP_x.  cAMP
##      raises PKA, PKA raises SIG, SIG raises PHOS, PHOS raises GI, and GI
##      divides AC.  The second factor is therefore negative once the apex is
##      past threshold - the same feedback term read at a higher density.
##      Levosimendan/OR-1896 enters CASENS, which multiplies myofilament
##      output and never touches cAMP, so its second factor is 1 by
##      construction.  Scenarios 8/9 titrate the two to the SAME added drive
##      (function `matched_drive()`), because comparing them at matched DOSE
##      compares potency and tells you nothing about routing.
##
##      Milrinone exists in this model for one reason: it raises cAMP
##      DOWNSTREAM of the receptor.  If the harm were a receptor-occupancy
##      artefact, milrinone would be safe.  The model says it is not.  That is
##      a falsifiable prediction (`falsification3()`).
##
##  (3) BETA-BLOCKADE HAS A CONDITIONAL SIGN, AND THE CONDITION IS A NUMBER IN
##      THE MODEL, NOT A CLINICAL RULE.
##
##      Lowering occupancy removes drive from the base (harm if the deficit is
##      forward flow) and removes phosphorylation from the apex (benefit if
##      the deficit is the switch).  Which dominates is decided by the LVOT
##      gradient term, whose existence depends on the patient parameter SEPT
##      (basal septal geometry) multiplied by basal hypercontractility.  The
##      registry observation - short-acting beta-blockade helps obstructive
##      TTS, and beta-blockade does NOT reduce recurrence whereas ACE
##      inhibition does - is a consequence of the routing, not a lookup.
##
##  (4) FOUR COUNTERFACTUAL INTEGRATORS, SO NO MECHANISM'S SHARE IS ASSUMED.
##
##      cAMP has a ~1-minute half-life, orders of magnitude faster than every
##      other state, so its quasi-steady-state value is available in closed
##      form.  That makes exact counterfactuals cheap:
##
##          HD_GI    EF if PHOS were 0        (the pertussis-toxin experiment)
##          HD_OED   EF if oedema were 0
##          HD_ATP   EF if energetics were normal
##          HD_LVOT  forward flow if the gradient were 0
##
##      Each is integrated alongside the real trajectory, so the reported
##      share of the EF deficit attributable to the Gi switch is MEASURED.
##      The QSS assumption is not taken on faith: CAMPQ_AP is captured next to
##      the integrated CAMP_AP so the two can be compared in every run.
##
##  ---------------------------------------------------------------------------
##  WHAT IS DELIBERATELY NOT IN HERE
##  ---------------------------------------------------------------------------
##  * Reverse (basal) takotsubo is NOT derived.  The model reproduces it only
##    by INVERTING RHO/FB2 (see `phenotypes$reverse`).  That is an honest
##    statement of a prediction - reverse TTS should have an inverted receptor
##    gradient - not a derivation, and it is flagged as such in the README.
##  * An isolated MID-VENTRICULAR ring cannot be produced by any monotonic
##    apex>mid>base gradient.  The model gives apical or apical+mid patterns.
##    Reported as a limitation, not tuned away.
##  * There is no plaque, no thrombosis-of-a-coronary and no infarct
##    compartment.  Troponin here comes from a cytosolic pool plus a small
##    necrotic term; that IS the claim being tested (small troponin, large
##    dysfunction), so building in a large necrosis term would beg it.
##
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY -
##  not validated for clinical or regulatory use.
## =============================================================================

library(mrgsolve)

tts_code <- r"---(
$PROB
# Takotsubo syndrome QSP model (64 ODEs, time unit = hours)

$PARAM @annotated
// ===================================================== SEGMENT ARCHITECTURE
// The ONLY inter-segment differences in the whole model.
FAP    : 0.28  : Apical mass fraction of LV (-)
FMD    : 0.34  : Mid mass fraction of LV (-)
FBS    : 0.38  : Basal mass fraction of LV (-)
RHO_AP : 1.40  : Total beta-AR density apex relative to base (-)
RHO_MD : 1.15  : Total beta-AR density mid relative to base (-)
RHO_BS : 1.00  : Total beta-AR density base (reference) (-)
FB2_AP : 0.42  : beta2 fraction of beta-AR pool apex (-)
FB2_MD : 0.32  : beta2 fraction of beta-AR pool mid (-)
FB2_BS : 0.24  : beta2 fraction of beta-AR pool base (-)
INN_AP : 0.62  : Sympathetic innervation apex rel. to mean (123I-MIBG) (-)
INN_MD : 0.90  : Sympathetic innervation mid (-)
INN_BS : 1.20  : Sympathetic innervation base (-)
TH_AP  : 0.36  : Wall thickness apex relative to base (-)
TH_MD  : 0.72  : Wall thickness mid (-)
TH_BS  : 1.00  : Wall thickness base (-)

// ================================================== TRIGGER AND SUSCEPTIBILITY
TRIG_T  : 0.0   : Trigger onset time (h)
AMP_TOT : 250   : Total catecholamine secretion amplitude (nmol/L/h at t=TRIG_T)
FRAC_E  : 0.75  : Fraction of the surge routed to ADRENALINE (adrenal) (-)
TAU_SUR : 1.6   : Surge decay time constant (h)
E2      : 15    : Oestradiol (pg/mL) - 15 postmenopausal, 100 premenopausal
E2K     : 30    : Oestradiol concentration halving sympathetic gain (pg/mL)
SEPT    : 0.35  : Basal septal geometry factor - sets whether LVOTO is possible (-)
QTDRUG  : 0.0   : QT-prolonging comedication burden (0-1)
PTX     : 0.0   : Pertussis toxin - experimental Gi blockade (0-1)

// ============================================== CATECHOLAMINE KINETICS (nmol/L)
EPI0     : 0.30  : Basal plasma adrenaline (nmol/L)
NE0      : 1.80  : Basal plasma noradrenaline (nmol/L)
NEI0     : 8.0   : Basal cardiac interstitial noradrenaline (nmol/L)
KE_EPI   : 20.8  : Adrenaline elimination rate (1/h ; t1/2 2 min)
KE_NE    : 16.6  : Noradrenaline elimination rate (1/h ; t1/2 2.5 min)
KNEI_OUT : 30.0  : Interstitial NE washout (1/h)
VMAX_UP1 : 250   : Uptake-1 Vmax (nmol/L/h)
KM_UP1   : 12.0  : Uptake-1 Km (nmol/L)
FNEI     : 10.0  : Interstitial amplification of local NE release (-)
VCAT     : 15.0  : Catecholamine distribution volume (L)

// ================================================= RECEPTOR PHARMACOLOGY
KD_B2  : 3.0   : Adrenaline Kd at beta2-AR (nmol/L)
KD_B1  : 60.0  : Noradrenaline Kd at beta1-AR (nmol/L)
KD_A1  : 120   : Noradrenaline Kd at alpha1-AR (nmol/L)
AN2    : 0.04  : Relative potency of noradrenaline at beta2 (-)
AE1    : 1.00  : Relative potency of adrenaline at beta1 (-)
AE_A1  : 0.25  : Relative potency of adrenaline at alpha1 (-)

// ============================================ THE SWITCH (Gs -> Gi trafficking)
THR_SW : 0.30    : beta2 occupancy threshold for Gi trafficking (-)
KSW    : 0.25    : Half-maximal supra-threshold occupancy (-)
KPH    : 0.20    : Receptor phosphorylation rate constant (1/h)
KDEPH  : 0.00722 : Receptor dephosphorylation rate (1/h ; t1/2 4.0 d)
WPKA   : 1.20    : PKA amplification of the phosphorylation drive (-)
WB2D   : 2.6336    : beta2-density amplification of microdomain drive (-)
KGRK_I : 0.020   : GRK2 upregulation rate (1/h)
KGRK_O : 0.020   : GRK2 turnover (1/h)
WGRK   : 0.60    : GRK2 contribution to phosphorylation drive (-)
KRDN_I : 0.004   : Receptor internalisation rate (1/h)
KRDN_O : 0.010   : Receptor recycling rate (1/h)
RDNMAX : 0.45    : Maximum receptor downregulation (-)

// ==================================================== cAMP / PKA / CALCIUM
AC_BAS   : 30.0  : Agonist-independent adenylyl cyclase flux (50 pct of resting flux)
EMAX_AC  : 432   : Maximal receptor-driven AC flux (cAMP units/h)
KM_GS    : 2.00  : Gs half-saturation of adenylyl cyclase (-)
GAM_GI   : 9.97164  : Gi inhibitory gain on adenylyl cyclase (-)
KDEG_CA  : 60.0  : cAMP degradation rate (1/h ; PDE3+PDE4)
FPDE3    : 0.55  : Fraction of cAMP degradation via PDE3 (-)
KPKA     : 2.00  : PKA activation half-constant on cAMP (-)
KCA_IN   : 2.00  : Calcium loading rate constant (1/h)
ECA_PKA  : 3.00  : PKA drive on calcium loading (-)
KCA_OUT  : 2.00  : Calcium extrusion rate constant (1/h)
KCA_ATP  : 0.20  : ATP dependence of calcium extrusion (-)

// ================================================= ENERGETICS AND INJURY
KATP     : 0.35  : Energetic reserve relaxation rate (1/h)
KO2_CONT : 0.55  : Contractility weight in oxygen demand (-)
KO2_HR   : 0.30  : Heart-rate weight in oxygen demand (-)
KO2_WS   : 0.35  : Wall-stress weight in oxygen demand (-)
CFR      : 2.50  : Coronary flow reserve - supply can exceed resting demand (-)
PAUTO    : 30.0  : Coronary autoregulation pressure intercept (mmHg)
PAUT50   : 20.0  : Coronary autoregulation half-pressure above intercept (mmHg)
KWSI     : 0.14843  : Systolic-stretch drive on membrane permeability (-)
WSTR_ST  : 0.60  : Systolic-stretch drive on stunning (-)
KMPTP    : 0.010 : mPTP-mediated ATP loss (1/h)
CATHR    : 1.35  : Calcium overload threshold for calpain activation (-)
KSTUN    : 0.0010: Stunning generation rate (1/h)
KSTUNR   : 0.00241: Stunning repair rate (1/h ; t1/2 12 d)
WBARR    : 0.35  : beta-arrestin protection of the stunning term (-)
KOED     : 0.0075: Oedema generation rate (1/h)
KOEDR    : 0.0116 : Oedema resolution rate (1/h ; t1/2 2.5 d)
KOEDC    : 0.45  : Oedema penalty on contractility (-)
WINF_OED : 0.30  : IL-6 contribution to oedema generation (-)
KNECR    : 0.0020: Necrosis generation rate (1/h)
NECTHR   : 0.55  : ATP threshold below which necrosis proceeds (-)
KLEAK    : 4.85864   : Cytosolic troponin leak coefficient (ng/mL per unit flux)
KNEC_TNI : 260   : Troponin release per unit necrosis (ng/mL)
KTNI_OUT : 0.35  : Plasma troponin elimination (1/h ; t1/2 2 h)
KIL6_IN  : 0.55  : IL-6 generation rate (pg/mL/h per unit injury)
KIL6_OUT : 0.14  : IL-6 elimination (1/h ; t1/2 5 h)
KCRP_IN  : 0.0126: CRP synthesis per IL-6 (mg/L per pg/mL per h)
KCRP_OUT : 0.021 : CRP elimination (1/h ; t1/2 33 h)

// ================================================ MECHANICS AND GEOMETRY
SH0      : 0.62  : Baseline segmental fractional shortening (-)
KDR      : 1.20  : cAMP half-constant for contractile drive (-)
NDR      : 1.60  : Hill coefficient of contractile drive (-)
KAFT0    : 0.35  : Afterload sensitivity of shortening at unit thickness (-)
KBULGE   : 0.1537  : Passive systolic bulging coefficient (-)
KBUL50   : 0.50  : Contractility deficit giving half-maximal bulging (-)
SHMIN    : -0.25 : Most negative attainable segmental shortening (-)
LVEDV0   : 105   : Baseline LV end-diastolic volume (mL)
KLV      : 0.10  : LVEDV relaxation rate (1/h)
KDILAT   : 0.55  : Apical akinesia drive on LVEDV (-)
KVPRE    : 0.45  : Volume drive on LVEDV (-)
EF0      : 0.62  : Baseline LVEF (-)

// ==================================================== LVOT OBSTRUCTION
GMAXG  : 150   : Maximal attainable LVOT gradient (mmHg)
KH50   : 0.25  : Effective basal hypercontractility at half-maximal gradient (-)
GTHR   : 1.02  : Basal contractility above which a gradient appears (-)
KGF    : 55.0  : Gradient producing half-maximal forward-flow loss (mmHg)
MRF0   : 0.04  : Baseline mitral regurgitant fraction (-)
KMR    : 0.35  : SAM-mediated mitral regurgitation coefficient (-)

// ======================================================= HAEMODYNAMICS
SV0     : 68.0  : Baseline stroke volume (mL)
HR0     : 72.0  : Baseline heart rate (1/min)
KHR     : 6.0   : Heart-rate relaxation rate (1/h)
EHR_B1  : 110   : beta1 occupancy drive on heart rate (1/min per unit occ)
KBARO   : 0.55  : Baroreflex heart-rate gain (1/min per mmHg)
KREFL   : 3300  : Maximal reflex drive on interstitial NE (nmol/L/h) - held FIXED
REFL50  : 0.15  : Cardiac-output deficit at half-maximal reflex drive (-)
REFLDZ  : 0.030 : Cardiac-output dead-zone before the reflex engages (-)
SVR0    : 1200  : Baseline systemic vascular resistance (dyn s/cm5)
KSVR    : 3.0   : SVR relaxation rate (1/h)
EA1_SVR : 0.45  : alpha1 occupancy drive on SVR (-)
EAII    : 0.35  : Angiotensin II drive on SVR (-)
ENO     : 0.30  : NO brake on SVR (-)
VOL0    : 5.00  : Baseline effective circulating volume (L)
KVIN    : 0.25  : Volume input rate (L/h)
KVOUT   : 0.05  : Volume output rate constant (1/h)
EALD    : 0.35  : Aldosterone drive on volume retention (-)
EBNP    : 0.15  : Maximal natriuretic-peptide drive on volume loss (-)
KBNPV   : 10.0  : NT-proBNP fold-rise at half-maximal natriuresis (-)
MAP0    : 80.0  : Reference mean arterial pressure (mmHg)
PCWP0   : 9.0   : Baseline pulmonary capillary wedge pressure (mmHg)
CVP     : 5.0   : Central venous pressure offset (mmHg)

// ======================================================= NEUROHORMONAL
KREN_I  : 0.35  : Renin-angiotensin generation rate (1/h)
KANG_O  : 0.35  : Angiotensin II elimination (1/h)
ANGII0  : 1.00  : Baseline angiotensin II (normalised)
ERENB1  : 2.20  : beta1 occupancy drive on renin (-)
EREN_MAP: 2.50  : Hypotension drive on renin (-)
KALD_I  : 0.10  : Aldosterone generation (1/h)
KALD_O  : 0.10  : Aldosterone elimination (1/h)
KBSYN   : 0.058 : proBNP transcription rate (1/h ; t1/2 12 h)
GBNP    : 7.03797  : NT-proBNP fold-rise per unit supra-normal transcription (-)
KBNP_O  : 0.416 : NT-proBNP elimination (1/h ; t1/2 100 min)
NTBNP0  : 90    : Baseline NT-proBNP (pg/mL)
NBNP    : 2.20  : Wall-stress exponent for proBNP transcription (-)

// ==================================================== ELECTROPHYSIOLOGY
QTC0    : 418   : Baseline QTc (ms)
KQT     : 0.050 : QTc remodelling rate (1/h ; t1/2 35 h - the lag that puts the QTc peak on day 2-3)
AQT_OED : 130   : Oedema drive on QTc (ms per unit oedema)
AQT_BAL : 136.375   : Regional dispersion (ballooning index) drive on QTc (ms)
AQT_CA  : 45    : Calcium overload drive on QTc (ms per unit)
AQT_D   : 40    : Comedication drive on QTc (ms)
AQT_K   : 28    : Hypokalaemia drive on QTc (ms per mmol/L below 4.0)
K0      : 4.10  : Baseline plasma potassium (mmol/L)
KK_IN   : 0.41  : Potassium input (mmol/L/h)
KK_OUT  : 0.10  : Potassium turnover (1/h)
KK_FUR  : 0.055 : Furosemide-driven potassium loss (mmol/L/h per unit effect)
KK_CAT  : 0.030 : Catecholamine-driven potassium shift (mmol/L/h per unit occ)
H0_TDP  : 1.42e-5 : Torsade hazard coefficient (1/h) - set to a 2-5 pct incidence
BQT     : 0.60  : Hazard exponent per 10 ms of QTc above 450 (-)

// ========================================================= THROMBUS
KTF     : 0.0022: Thrombus formation rate (units/h per unit stasis)
KTL     : 0.0035: Thrombus lysis rate (1/h)
KAC_LYS : 2.60  : Anticoagulant amplification of lysis (-)
KHYPERC : 0.45  : Catecholamine hypercoagulability contribution (-)
H0_EMB  : 0.0016: Embolic hazard per unit thrombus per hour (1/h)

// ========================================================= RECURRENCE
H0_REC  : 2.6e-6: Baseline recurrence hazard (1/h ; ~2.3%/yr)
W_PHREC : 3.5   : Residual phosphorylation drive on recurrence (-)
E_ACEREC: 0.55  : ACE-inhibitor reduction of recurrence hazard (-)
E_BBREC : 0.05  : beta-blocker reduction of recurrence hazard (-)

// ============================================ DRUG PK (amounts mg, conc mg/L)
// Dobutamine
V_DOB    : 25.0  : Dobutamine volume of distribution (L)
CL_DOB   : 520   : Dobutamine clearance (L/h ; t1/2 ~2 min)
EC50_DOB : 0.10  : Dobutamine EC50 at beta1 (mg/L)
EMAX_DOB : 0.75  : Dobutamine intrinsic activity at beta1 (partial) (-)
FB2_DOB  : 0.25  : Dobutamine beta2 potency relative to beta1 (-)
RATE_DOB : 0.0   : Dobutamine infusion rate (mg/h)
// Exogenous catecholamines (ug/min)
RATE_EPI : 0.0   : Adrenaline infusion (ug/min)
RATE_NEP : 0.0   : Noradrenaline infusion (ug/min)
KEPI_INF : 21.9  : nmol/L/h of adrenaline per ug/min infused (-)
KNEP_INF : 23.7  : nmol/L/h of noradrenaline per ug/min infused (-)
// Milrinone
V_MIL    : 31.0  : Milrinone volume of distribution (L)
CL_MIL   : 9.1   : Milrinone clearance (L/h)
IC50_MIL : 0.20  : Milrinone IC50 on PDE3 (mg/L)
IMAX_MIL : 0.90  : Milrinone maximal PDE3 inhibition (-)
RATE_MIL : 0.0   : Milrinone infusion rate (mg/h)
// Levosimendan and OR-1896
V_LEV    : 14.0  : Levosimendan central volume (L)
VP_LEV   : 20.0  : Levosimendan peripheral volume (L)
Q_LEV    : 8.0   : Levosimendan intercompartmental clearance (L/h)
CL_LEV   : 13.0  : Levosimendan clearance (L/h ; t1/2 ~1 h)
FM_LEV   : 0.05  : Fraction of levosimendan converted to OR-1896 (-)
V_OR     : 48.0  : OR-1896 volume of distribution (L)
CL_OR    : 0.43  : OR-1896 clearance (L/h ; t1/2 ~77 h)
EC50_LEV : 0.030 : Levosimendan EC50 for Ca sensitisation (mg/L)
EMAX_LEV : 0.55  : Maximal myofilament Ca sensitisation (-)
POT_OR   : 1.20  : OR-1896 potency relative to parent (-)
KATP_LEV : 0.30  : Levosimendan K_ATP vasodilator effect on SVR (-)
RATE_LEV : 0.0   : Levosimendan infusion rate (mg/h)
// Esmolol
V_ESM    : 240   : Esmolol volume of distribution (L)
CL_ESM   : 1200  : Esmolol clearance (L/h ; t1/2 ~8 min)
KI1_ESM  : 0.15  : Esmolol Ki at beta1 (mg/L)
KI2_ESM  : 5.00  : Esmolol Ki at beta2 (mg/L)
RATE_ESM : 0.0   : Esmolol infusion rate (mg/h)
// Metoprolol (oral)
KA_MET   : 1.50  : Metoprolol absorption rate (1/h)
F_MET    : 0.40  : Metoprolol bioavailability (-)
V_MET    : 280   : Metoprolol volume of distribution (L)
CL_MET   : 63.0  : Metoprolol clearance (L/h)
KI1_MET  : 0.060 : Metoprolol Ki at beta1 (mg/L)
KI2_MET  : 1.500 : Metoprolol Ki at beta2 (mg/L)
// Carvedilol (oral)
KA_CAR   : 1.00  : Carvedilol absorption rate (1/h)
F_CAR    : 0.25  : Carvedilol bioavailability (-)
V_CAR    : 115   : Carvedilol volume of distribution (L)
CL_CAR   : 12.0  : Carvedilol clearance (L/h)
KI1_CAR  : 0.040 : Carvedilol Ki at beta1, TOTAL-concentration referenced (mg/L)
KI2_CAR  : 0.080 : Carvedilol Ki at beta2, total-concentration referenced (mg/L)
KIA_CAR  : 0.200 : Carvedilol Ki at alpha1, total-concentration referenced (mg/L)
// Phenylephrine
V_PHE    : 25.0  : Phenylephrine volume of distribution (L)
CL_PHE   : 200   : Phenylephrine clearance (L/h ; t1/2 ~5 min)
POT_PHE  : 45.0  : Phenylephrine alpha1 potency (nmol/L equivalent per mg/L)
RATE_PHE : 0.0   : Phenylephrine infusion rate (mg/h)
// Ramipril / ramiprilat
KA_RAM   : 1.20  : Ramipril absorption rate (1/h)
F_RAM    : 0.55  : Ramipril bioavailability to ramiprilat (-)
V_RAM    : 100   : Ramiprilat volume of distribution (L)
CL_RAM   : 6.00  : Ramiprilat clearance (L/h)
IC50_RAM : 0.0020: Ramiprilat IC50 on ACE (mg/L)
IMAX_RAM : 0.85  : Maximal ACE inhibition (-)
// Furosemide
KA_FUR   : 1.50  : Furosemide absorption rate (1/h)
F_FUR    : 0.50  : Furosemide bioavailability (-)
V_FUR    : 12.0  : Furosemide volume of distribution (L)
CL_FUR   : 9.00  : Furosemide clearance (L/h)
EC50_FUR : 0.50  : Furosemide EC50 for natriuresis (mg/L)
KFUR_V   : 0.55  : Maximal furosemide volume loss (L/h)
// Apixaban
KA_APX   : 1.00  : Apixaban absorption rate (1/h)
F_APX    : 0.50  : Apixaban bioavailability (-)
V_APX    : 21.0  : Apixaban volume of distribution (L)
CL_APX   : 1.30  : Apixaban clearance (L/h ; t1/2 ~11 h)
EC50_APX : 0.10  : Apixaban EC50 for anticoagulation (mg/L)
// Mechanical support and adjuncts (switches / fractions)
IABP     : 0.0   : IABP support (0-1) - lowers afterload
IMPELLA  : 0.0   : Impella support (0-1) - direct forward flow
E_IABP   : 0.22  : IABP afterload reduction (-)
E_IMP_CO : 0.45  : Impella forward-flow contribution (-)
E_IMP_UNL: 0.18  : Impella LV unloading (-)
SEDATE   : 0.0   : Anxiolysis/sedation (0-1) - cuts the afferent trigger
E_SEDATE : 0.45  : Maximal trigger attenuation by sedation (-)
KSUP      : 0.0  : Potassium/magnesium repletion (mmol/L/h)

$INIT @annotated
// --------------------------------------------------------------- catecholamines
EPI    : 0.30  : Plasma adrenaline (nmol/L)
NE     : 1.80  : Plasma noradrenaline (nmol/L)
NEI    : 8.00  : Cardiac interstitial noradrenaline (nmol/L)
// -------------------------------------------------------------------- drug PK
DOB    : 0.0   : Dobutamine amount (mg)
MIL    : 0.0   : Milrinone amount (mg)
LEVC   : 0.0   : Levosimendan central amount (mg)
LEVP   : 0.0   : Levosimendan peripheral amount (mg)
ORM    : 0.0   : OR-1896 amount (mg)
ESM    : 0.0   : Esmolol amount (mg)
METD   : 0.0   : Metoprolol depot (mg)
METC   : 0.0   : Metoprolol central amount (mg)
CARD   : 0.0   : Carvedilol depot (mg)
CARC   : 0.0   : Carvedilol central amount (mg)
PHE    : 0.0   : Phenylephrine amount (mg)
RAMD   : 0.0   : Ramipril depot (mg)
RAMC   : 0.0   : Ramiprilat amount (mg)
FURD   : 0.0   : Furosemide depot (mg)
FURC   : 0.0   : Furosemide central amount (mg)
APXD   : 0.0   : Apixaban depot (mg)
APXC   : 0.0   : Apixaban central amount (mg)
// ------------------------------------------------------------ receptor / cAMP
CAMP_AP: 1.00  : cAMP apex (normalised)
CAMP_MD: 1.00  : cAMP mid (normalised)
CAMP_BS: 1.00  : cAMP base (normalised)
PHOS_AP: 0.0   : Gi-coupled beta2 fraction apex (-)
PHOS_MD: 0.0   : Gi-coupled beta2 fraction mid (-)
PHOS_BS: 0.0   : Gi-coupled beta2 fraction base (-)
GRK2   : 1.00  : GRK2 activity (normalised)
RDN    : 0.0   : beta-AR downregulation fraction (-)
// ------------------------------------------------------------- injury / tissue
CAL_AP : 1.00  : Calcium overload index apex (-)
CAL_MD : 1.00  : Calcium overload index mid (-)
CAL_BS : 1.00  : Calcium overload index base (-)
ATP_AP : 1.00  : Energetic reserve apex (-)
ATP_MD : 1.00  : Energetic reserve mid (-)
ATP_BS : 1.00  : Energetic reserve base (-)
OED_AP : 0.0   : Myocardial oedema apex (-)
OED_MD : 0.0   : Myocardial oedema mid (-)
OED_BS : 0.0   : Myocardial oedema base (-)
STU_AP : 0.0   : Stunning apex (-)
STU_MD : 0.0   : Stunning mid (-)
STU_BS : 0.0   : Stunning base (-)
NECR   : 0.0   : Cumulative necrosis (fraction of LV)
// ------------------------------------------------------------------ biomarkers
TNI    : 0.010 : Plasma hs-cTnI (ng/mL)
IL6    : 2.50  : Plasma IL-6 (pg/mL)
CRP    : 1.50  : Plasma CRP (mg/L)
BSYN   : 1.00  : proBNP transcriptional signal (normalised)
NTBNP  : 90.0  : Plasma NT-proBNP (pg/mL)
// -------------------------------------------------------------- neurohormonal
ANGII  : 1.00  : Angiotensin II (normalised)
ALDO   : 1.00  : Aldosterone (normalised)
// -------------------------------------------------------------- haemodynamics
SVR    : 1200  : Systemic vascular resistance (dyn s/cm5)
VOL    : 5.00  : Effective circulating volume (L)
HR     : 72.0  : Heart rate (1/min)
LVEDV  : 105   : LV end-diastolic volume (mL)
// --------------------------------------------------------- electrophysiology
KPL    : 4.10  : Plasma potassium (mmol/L)
QTC    : 418   : QTc interval (ms)
HTDP   : 0.0   : Cumulative torsade hazard (-)
// ------------------------------------------------------------------- thrombus
THR    : 0.0   : LV apical thrombus (arbitrary mass units)
HEMB   : 0.0   : Cumulative embolic hazard (-)
HREC   : 0.0   : Cumulative recurrence hazard (-)
// --------------------------------------------------- inference / bookkeeping
AUCEPI : 0.0   : Cumulative adrenaline exposure (nmol/L*h)
CHDGI  : 0.0   : Time-integral of the Gi-attributable EF deficit (EF*h)
CHDOED : 0.0   : Time-integral of the oedema-attributable EF deficit (EF*h)
CHDATP : 0.0   : Time-integral of the energetic EF deficit (EF*h)
CHDLVO : 0.0   : Time-integral of the gradient-attributable SV deficit (mL*h)
CEFDEF : 0.0   : Time-integral of the total EF deficit (EF*h)

$GLOBAL
#define POSPART(x) ((x) > 0.0 ? (x) : 0.0)

/* ---------------------------------------------------------------------------
   THE ALGEBRAIC LAYER, DEFINED EXACTLY ONCE.

   Every quantity that both $ODE and $TABLE need is computed by this macro.
   The 2026-07 lesson from this repository was that a $TABLE block which
   recomputes a flux with even slightly different gating reports a trajectory
   that is not the one being integrated.  A macro cannot drift: there is one
   text, expanded in both places.
   --------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
   Shared algebra variables.  mrgsolve promotes every block-local declaration
   to namespace scope, so a macro that DECLARED its variables could only be
   expanded once.  Declaring them here, and letting the macro do assignments
   only, is what allows $ODE and $TABLE to run literally the same text.
   --------------------------------------------------------------------------- */
double SGAIN, NOF;
/* Resting reference values, recomputed per subject in $MAIN so that the
   untreated baseline is EXACT by construction rather than by hand-typed
   initial conditions.  Every normalisation below divides by one of these. */
double GS0_AP, GS0_MD, GS0_BS, CAMP0_AP, CAMP0_MD, CAMP0_BS;
double PKA0_AP, PKA0_MD, PKA0_BS, DRR_AP, DRR_MD, DRR_BS;
double O10_AP, O10_MD, O10_BS, O20_AP, O20_MD, O20_BS, OCC0, SV0R, CO0R;
double AVAIL, C_ESM, C_MET, C_CAR, C_DOB, C_MIL;
double C_LEV, C_OR, C_PHE, C_RAM, C_FUR, C_APX;
double INH1, INH2, DOBEF, NEL_AP, NEL_MD, NEL_BS;
double AG2_AP, AG2_MD, AG2_BS, O2_AP, O2_MD, O2_BS;
double AG1_AP, AG1_MD, AG1_BS, O1_AP, O1_MD, O1_BS;
double B1D_AP, B1D_MD, B1D_BS, B2D_AP, B2D_MD, B2D_BS;
double PKA_AP, PKA_MD, PKA_BS, SPA, SPM, SPB;
double SIG_AP, SIG_MD, SIG_BS, GIF, GS_AP, GS_MD;
double GS_BS, GI_AP, GI_MD, GI_BS, PDE3INH, KDEG;
double AC_AP, AC_MD, AC_BS, CQ_AP, CQ_MD, CQ_BS;
double LEVEQ, CASENS, ATPF_AP, ATPF_MD, ATPF_BS, DR0;
double DR_AP, DR_MD, DR_BS, CONT_AP, CONT_MD, CONT_BS;
double AG_A1, INHA1, OA1, OA1_0, MAPX, RAD;
double WSF_AP, WSF_MD, WSF_BS, WS0_AP, WS0_MD, WS0_BS;
double SHA_AP, SHA_MD, SHA_BS, BUL_AP, BUL_MD, BUL_BS;
double SH_AP, SH_MD, SH_BS, EFX, BALL, HYPB;
double GRADX, FFLOSS, MRF, SVX, COX, MAPY;
double PCWPX, IMRX, PERF_AP, PERF_MD, PERF_BS, DEM_AP;
double PAUTF, STR_AP, STR_MD, STR_BS;
double DEM_MD, DEM_BS, PRM_AP, PRM_MD, PRM_BS, LEAKF;
double ACEINH, FUREF, ANTIC, ACN_AP, ACN_MD, ACN_BS;
double CN_AP, CN_MD, CN_BS, DN_AP, DN_MD, DN_BS;
double CTN_AP, CTN_MD, CTN_BS, EF_NOGI, CO2_AP, CO2_MD;
double CO2_BS, EF_NOED, CO3_AP, CO3_MD, CO3_BS, EF_NATP;
double SV_NOLV, STASIS, HYPERC;

#define TTS_ALGEBRA                                                           \
  /* --- oestradiol sets the sympathetic gain and the NO brake --- */         \
         SGAIN   = (1.0/(1.0 + E2/E2K))*(1.0 - E_SEDATE*SEDATE);              \
         NOF     = E2/(E2 + E2K);                                             \
  /* --- receptor availability after downregulation and antagonists --- */    \
         AVAIL   = 1.0 - RDN;                                                 \
         C_ESM   = ESM  / V_ESM;                                              \
         C_MET   = METC / V_MET;                                              \
         C_CAR   = CARC / V_CAR;                                              \
         C_DOB   = DOB  / V_DOB;                                              \
         C_MIL   = MIL  / V_MIL;                                              \
         C_LEV   = LEVC / V_LEV;                                              \
         C_OR    = ORM  / V_OR;                                               \
         C_PHE   = PHE  / V_PHE;                                              \
         C_RAM   = RAMC / V_RAM;                                              \
         C_FUR   = FURC / V_FUR;                                              \
         C_APX   = APXC / V_APX;                                              \
         INH1    = 1.0 / (1.0 + C_ESM/KI1_ESM + C_MET/KI1_MET + C_CAR/KI1_CAR); \
         INH2    = 1.0 / (1.0 + C_ESM/KI2_ESM + C_MET/KI2_MET + C_CAR/KI2_CAR); \
  /* --- dobutamine acts as an added beta agonist, not as a cAMP source --- */ \
         DOBEF   = EMAX_DOB * C_DOB / (C_DOB + EC50_DOB);                     \
  /* --- local agonist concentrations per segment --- */                      \
         NEL_AP  = NEI * INN_AP;                                              \
         NEL_MD  = NEI * INN_MD;                                              \
         NEL_BS  = NEI * INN_BS;                                              \
         AG2_AP  = EPI + AN2*NEL_AP;                                          \
         AG2_MD  = EPI + AN2*NEL_MD;                                          \
         AG2_BS  = EPI + AN2*NEL_BS;                                          \
         O2_AP   = AVAIL*INH2*(AG2_AP/(AG2_AP + KD_B2) + FB2_DOB*DOBEF);      \
         O2_MD   = AVAIL*INH2*(AG2_MD/(AG2_MD + KD_B2) + FB2_DOB*DOBEF);      \
         O2_BS   = AVAIL*INH2*(AG2_BS/(AG2_BS + KD_B2) + FB2_DOB*DOBEF);      \
         AG1_AP  = AE1*EPI + NEL_AP;                                          \
         AG1_MD  = AE1*EPI + NEL_MD;                                          \
         AG1_BS  = AE1*EPI + NEL_BS;                                          \
         O1_AP   = AVAIL*INH1*(AG1_AP/(AG1_AP + KD_B1) + DOBEF);              \
         O1_MD   = AVAIL*INH1*(AG1_MD/(AG1_MD + KD_B1) + DOBEF);              \
         O1_BS   = AVAIL*INH1*(AG1_BS/(AG1_BS + KD_B1) + DOBEF);              \
  /* --- receptor densities --- */                                            \
         B1D_AP  = RHO_AP*(1.0 - FB2_AP);                                     \
         B1D_MD  = RHO_MD*(1.0 - FB2_MD);                                     \
         B1D_BS  = RHO_BS*(1.0 - FB2_BS);                                     \
         B2D_AP  = RHO_AP*FB2_AP;                                             \
         B2D_MD  = RHO_MD*FB2_MD;                                             \
         B2D_BS  = RHO_BS*FB2_BS;                                             \
  /* --- PKA from the integrated cAMP --- */                                  \
         PKA_AP  = CAMP_AP*CAMP_AP/(CAMP_AP*CAMP_AP + KPKA*KPKA);             \
         PKA_MD  = CAMP_MD*CAMP_MD/(CAMP_MD*CAMP_MD + KPKA*KPKA);             \
         PKA_BS  = CAMP_BS*CAMP_BS/(CAMP_BS*CAMP_BS + KPKA*KPKA);             \
  /* --- THE SWITCH: one function, three densities --- */                     \
         SPA     = POSPART(O2_AP - THR_SW);                                   \
         SPM     = POSPART(O2_MD - THR_SW);                                   \
         SPB     = POSPART(O2_BS - THR_SW);                                   \
         SIG_AP  = SPA*SPA/(SPA*SPA + KSW*KSW)                                \
                   * (1.0 + WPKA*PKA_AP + WGRK*(GRK2 - 1.0))                  \
                   * (1.0 + WB2D*(B2D_AP/B2D_MD - 1.0));                      \
         SIG_MD  = SPM*SPM/(SPM*SPM + KSW*KSW)                                \
                   * (1.0 + WPKA*PKA_MD + WGRK*(GRK2 - 1.0))                  \
                   * (1.0 + WB2D*(B2D_MD/B2D_MD - 1.0));                      \
         SIG_BS  = SPB*SPB/(SPB*SPB + KSW*KSW)                                \
                   * (1.0 + WPKA*PKA_BS + WGRK*(GRK2 - 1.0))                  \
                   * (1.0 + WB2D*(B2D_BS/B2D_MD - 1.0));                      \
  SIG_AP = POSPART(SIG_AP); SIG_MD = POSPART(SIG_MD); SIG_BS = POSPART(SIG_BS); \
  /* --- Gs and Gi coupling (PTX abolishes Gi: the control experiment) --- */ \
         GIF     = 1.0 - PTX;                                                 \
         GS_AP   = B1D_AP*O1_AP + B2D_AP*O2_AP*(1.0 - PHOS_AP*GIF);           \
         GS_MD   = B1D_MD*O1_MD + B2D_MD*O2_MD*(1.0 - PHOS_MD*GIF);           \
         GS_BS   = B1D_BS*O1_BS + B2D_BS*O2_BS*(1.0 - PHOS_BS*GIF);           \
         GI_AP   = B2D_AP*O2_AP*PHOS_AP*GIF;                                  \
         GI_MD   = B2D_MD*O2_MD*PHOS_MD*GIF;                                  \
         GI_BS   = B2D_BS*O2_BS*PHOS_BS*GIF;                                  \
  /* --- adenylyl cyclase and PDE (milrinone enters the PDE3 term) --- */     \
         PDE3INH = IMAX_MIL*C_MIL/(C_MIL + IC50_MIL);                         \
         KDEG    = KDEG_CA*(1.0 - FPDE3*PDE3INH);                             \
         AC_AP   = (AC_BAS + EMAX_AC*GS_AP/(GS_AP + KM_GS))/(1.0 + GAM_GI*GI_AP); \
         AC_MD   = (AC_BAS + EMAX_AC*GS_MD/(GS_MD + KM_GS))/(1.0 + GAM_GI*GI_MD); \
         AC_BS   = (AC_BAS + EMAX_AC*GS_BS/(GS_BS + KM_GS))/(1.0 + GAM_GI*GI_BS); \
  /* --- quasi-steady-state cAMP: exact because t1/2(cAMP) ~ 1 min --- */     \
         CQ_AP   = AC_AP/KDEG;                                                \
         CQ_MD   = AC_MD/KDEG;                                                \
         CQ_BS   = AC_BS/KDEG;                                                \
  /* --- myofilament Ca sensitisation: levosimendan route, no cAMP --- */     \
         LEVEQ   = C_LEV + POT_OR*C_OR;                                       \
         CASENS  = EMAX_LEV*LEVEQ/(LEVEQ + EC50_LEV);                         \
  /* --- energetic and structural gates --- */                                \
         ATPF_AP = (ATP_AP/(ATP_AP + 0.15))*1.15;                             \
         ATPF_MD = (ATP_MD/(ATP_MD + 0.15))*1.15;                             \
         ATPF_BS = (ATP_BS/(ATP_BS + 0.15))*1.15;                             \
         DR0     = pow(1.0, NDR)/(pow(1.0, NDR) + pow(KDR, NDR));             \
         DR_AP   = pow(CAMP_AP,NDR)/(pow(CAMP_AP,NDR) + pow(KDR,NDR));        \
         DR_MD   = pow(CAMP_MD,NDR)/(pow(CAMP_MD,NDR) + pow(KDR,NDR));        \
         DR_BS   = pow(CAMP_BS,NDR)/(pow(CAMP_BS,NDR) + pow(KDR,NDR));        \
         CONT_AP = (DR_AP/DRR_AP)*ATPF_AP*(1.0 - KOEDC*OED_AP)                \
                   *(1.0 - STU_AP)*(1.0 + CASENS);                            \
         CONT_MD = (DR_MD/DRR_MD)*ATPF_MD*(1.0 - KOEDC*OED_MD)                \
                   *(1.0 - STU_MD)*(1.0 + CASENS);                            \
         CONT_BS = (DR_BS/DRR_BS)*ATPF_BS*(1.0 - KOEDC*OED_BS)                \
                   *(1.0 - STU_BS)*(1.0 + CASENS);                            \
  CONT_AP = POSPART(CONT_AP); CONT_MD = POSPART(CONT_MD);                     \
  CONT_BS = POSPART(CONT_BS);                                                 \
  /* --- alpha1 tone and afterload --- */                                     \
         AG_A1   = NE + AE_A1*EPI + POT_PHE*C_PHE;                            \
         INHA1   = 1.0/(1.0 + C_CAR/KIA_CAR);                                 \
         OA1     = INHA1*AG_A1/(AG_A1 + KD_A1);                               \
         OA1_0   = (NE0 + AE_A1*EPI0)/(NE0 + AE_A1*EPI0 + KD_A1);             \
  /* --- systolic pressure proxy and wall stress --- */                       \
         MAPX    = 0.0;                                                       \
         RAD     = pow(LVEDV/LVEDV0, 1.0/3.0);                                \
  /* --- geometry-dependent shortening; the apex is thin, hence afterload     \
         sensitive and passively distensible (Laplace, not an assumption) --- */ \
         WSF_AP  = RAD/TH_AP;                                                 \
         WSF_MD  = RAD/TH_MD;                                                 \
         WSF_BS  = RAD/TH_BS;                                                 \
         WS0_AP  = 1.0/TH_AP;                                                 \
         WS0_MD  = 1.0/TH_MD;                                                 \
         WS0_BS  = 1.0/TH_BS;                                                 \
         SHA_AP  = SH0*CONT_AP/(1.0 + (KAFT0/TH_AP)*(WSF_AP/WS0_AP - 1.0));   \
         SHA_MD  = SH0*CONT_MD/(1.0 + (KAFT0/TH_MD)*(WSF_MD/WS0_MD - 1.0));   \
         SHA_BS  = SH0*CONT_BS/(1.0 + (KAFT0/TH_BS)*(WSF_BS/WS0_BS - 1.0));   \
         BUL_AP  = KBULGE*(1.0/TH_AP)*POSPART(1.0 - CONT_AP)                  \
                   /(POSPART(1.0 - CONT_AP) + KBUL50);                        \
         BUL_MD  = KBULGE*(1.0/TH_MD)*POSPART(1.0 - CONT_MD)                  \
                   /(POSPART(1.0 - CONT_MD) + KBUL50);                        \
         BUL_BS  = KBULGE*(1.0/TH_BS)*POSPART(1.0 - CONT_BS)                  \
                   /(POSPART(1.0 - CONT_BS) + KBUL50);                        \
         SH_AP   = fmax(SHA_AP - BUL_AP, SHMIN);                              \
         SH_MD   = fmax(SHA_MD - BUL_MD, SHMIN);                              \
         SH_BS   = fmax(SHA_BS - BUL_BS, SHMIN);                              \
         EFX     = FAP*SH_AP + FMD*SH_MD + FBS*SH_BS;                         \
         BALL    = SH_BS - SH_AP;                                             \
  /* --- LVOT gradient: patient geometry x basal hypercontractility --- */    \
         HYPB    = POSPART(CONT_BS - GTHR)*pow(VOL0/VOL, 0.75)                \
                   * sqrt(SVR0/SVR) * (1.0 + E_IABP*IABP);                    \
  /* Saturating: a dynamic gradient is bounded by the pressure the basal      \
     segment can generate, not by the square of its hypercontractility. */    \
         GRADX   = GMAXG*SEPT*HYPB*HYPB/(HYPB*HYPB + KH50*KH50);              \
         FFLOSS  = GRADX/(GRADX + KGF);                                       \
         MRF     = MRF0 + KMR*FFLOSS;                                         \
  /* --- forward flow --- */                                                  \
         SVX     = SV0R*(EFX/EF0)*(LVEDV/LVEDV0)*(1.0 - FFLOSS)*(1.0 - MRF)   \
                   + SV0R*E_IMP_CO*IMPELLA;                                   \
  SVX = POSPART(SVX);                                                         \
         COX     = SVX*HR/1000.0;                                             \
         MAPY    = COX*SVR/80.0 + CVP;                                        \
  MAPX = MAPY;                                                                \
         PCWPX   = PCWP0*pow(LVEDV/LVEDV0, 2.2)*(1.0 + 0.8*MRF)               \
                   /(1.0 + 0.8*MRF0);                                         \
  /* --- perfusion: microvascular constriction plus driving pressure --- */   \
         IMRX    = 1.0 + 1.9*POSPART(OA1 - OA1_0)/fmax(OA1_0,1e-6);           \
  /* coronary autoregulation: flow is nearly pressure-independent above about \
     60 mmHg, so a fall in MAP cannot by itself start an energetic spiral */  \
         PAUTF   = (POSPART(MAPX - PAUTO)/(POSPART(MAPX - PAUTO) + PAUT50))   \
                   /((MAP0 - PAUTO)/((MAP0 - PAUTO) + PAUT50));               \
         PERF_AP = PAUTF/IMRX;                                                \
         PERF_MD = PAUTF/IMRX;                                                \
         PERF_BS = PAUTF/IMRX;                                                \
  /* --- oxygen demand per segment --- */                                     \
         DEM_AP  = KO2_CONT*CONT_AP + KO2_HR*(HR/HR0) + KO2_WS*(WSF_AP/WS0_AP); \
         DEM_MD  = KO2_CONT*CONT_MD + KO2_HR*(HR/HR0) + KO2_WS*(WSF_MD/WS0_MD); \
         DEM_BS  = KO2_CONT*CONT_BS + KO2_HR*(HR/HR0) + KO2_WS*(WSF_BS/WS0_BS); \
  /* --- membrane permeability flux driving the cytosolic troponin leak --- */ \
  /* Systolic stretch: 0 in a normally contracting segment, 1 when akinetic,  \
     >1 when dyskinetic.  Injury, oedema and the troponin leak are driven by  \
     it as well as by calcium, which is what makes them APICAL-PREDOMINANT    \
     without any segment-specific injury parameter: the akinesis comes first  \
     and the thin wall then carries the stress. */                            \
         STR_AP  = POSPART(1.0 - SH_AP/SH0);                                  \
         STR_MD  = POSPART(1.0 - SH_MD/SH0);                                  \
         STR_BS  = POSPART(1.0 - SH_BS/SH0);                                  \
         PRM_AP  = POSPART(CAL_AP - 1.0) + POSPART(1.0 - ATP_AP)              \
                   + KWSI*STR_AP*RAD/TH_AP*TH_BS;                             \
         PRM_MD  = POSPART(CAL_MD - 1.0) + POSPART(1.0 - ATP_MD)              \
                   + KWSI*STR_MD*RAD/TH_MD*TH_BS;                             \
         PRM_BS  = POSPART(CAL_BS - 1.0) + POSPART(1.0 - ATP_BS)              \
                   + KWSI*STR_BS*RAD/TH_BS*TH_BS;                             \
         LEAKF   = FAP*PRM_AP + FMD*PRM_MD + FBS*PRM_BS;                      \
  /* --- drug effects on the neurohormonal and volume arms --- */             \
         ACEINH  = IMAX_RAM*C_RAM/(C_RAM + IC50_RAM);                         \
         FUREF   = C_FUR/(C_FUR + EC50_FUR);                                  \
         ANTIC   = C_APX/(C_APX + EC50_APX);                                  \
  /* --- counterfactuals, using the exact QSS cAMP --- */                     \
         ACN_AP  = AC_BAS + EMAX_AC*(B1D_AP*O1_AP + B2D_AP*O2_AP)             \
                   /((B1D_AP*O1_AP + B2D_AP*O2_AP) + KM_GS);                  \
         ACN_MD  = AC_BAS + EMAX_AC*(B1D_MD*O1_MD + B2D_MD*O2_MD)             \
                   /((B1D_MD*O1_MD + B2D_MD*O2_MD) + KM_GS);                  \
         ACN_BS  = AC_BAS + EMAX_AC*(B1D_BS*O1_BS + B2D_BS*O2_BS)             \
                   /((B1D_BS*O1_BS + B2D_BS*O2_BS) + KM_GS);                  \
         CN_AP   = ACN_AP/KDEG;                                               \
         CN_MD   = ACN_MD/KDEG;                                               \
         CN_BS   = ACN_BS/KDEG;                                               \
         DN_AP   = pow(CN_AP,NDR)/(pow(CN_AP,NDR) + pow(KDR,NDR));            \
         DN_MD   = pow(CN_MD,NDR)/(pow(CN_MD,NDR) + pow(KDR,NDR));            \
         DN_BS   = pow(CN_BS,NDR)/(pow(CN_BS,NDR) + pow(KDR,NDR));            \
         CTN_AP  = POSPART((DN_AP/DRR_AP)*ATPF_AP*(1.0-KOEDC*OED_AP)          \
                   *(1.0-STU_AP)*(1.0+CASENS));                               \
         CTN_MD  = POSPART((DN_MD/DRR_MD)*ATPF_MD*(1.0-KOEDC*OED_MD)          \
                   *(1.0-STU_MD)*(1.0+CASENS));                               \
         CTN_BS  = POSPART((DN_BS/DRR_BS)*ATPF_BS*(1.0-KOEDC*OED_BS)          \
                   *(1.0-STU_BS)*(1.0+CASENS));                               \
         EF_NOGI = FAP*(SH0*CTN_AP/(1.0+(KAFT0/TH_AP)*(WSF_AP/WS0_AP-1.0))    \
                        - KBULGE*(1.0/TH_AP)*POSPART(1.0-CTN_AP)              \
                          /(POSPART(1.0-CTN_AP)+KBUL50))                      \
                 + FMD*(SH0*CTN_MD/(1.0+(KAFT0/TH_MD)*(WSF_MD/WS0_MD-1.0))    \
                        - KBULGE*(1.0/TH_MD)*POSPART(1.0-CTN_MD)              \
                          /(POSPART(1.0-CTN_MD)+KBUL50))                      \
                 + FBS*(SH0*CTN_BS/(1.0+(KAFT0/TH_BS)*(WSF_BS/WS0_BS-1.0))    \
                        - KBULGE*(1.0/TH_BS)*POSPART(1.0-CTN_BS)              \
                          /(POSPART(1.0-CTN_BS)+KBUL50));                     \
         CO2_AP  = POSPART((DR_AP/DRR_AP)*ATPF_AP*(1.0-STU_AP)*(1.0+CASENS)); \
         CO2_MD  = POSPART((DR_MD/DRR_MD)*ATPF_MD*(1.0-STU_MD)*(1.0+CASENS)); \
         CO2_BS  = POSPART((DR_BS/DRR_BS)*ATPF_BS*(1.0-STU_BS)*(1.0+CASENS)); \
         EF_NOED = FAP*(SH0*CO2_AP/(1.0+(KAFT0/TH_AP)*(WSF_AP/WS0_AP-1.0))    \
                        - KBULGE*(1.0/TH_AP)*POSPART(1.0-CO2_AP)              \
                          /(POSPART(1.0-CO2_AP)+KBUL50))                      \
                 + FMD*(SH0*CO2_MD/(1.0+(KAFT0/TH_MD)*(WSF_MD/WS0_MD-1.0))    \
                        - KBULGE*(1.0/TH_MD)*POSPART(1.0-CO2_MD)              \
                          /(POSPART(1.0-CO2_MD)+KBUL50))                      \
                 + FBS*(SH0*CO2_BS/(1.0+(KAFT0/TH_BS)*(WSF_BS/WS0_BS-1.0))    \
                        - KBULGE*(1.0/TH_BS)*POSPART(1.0-CO2_BS)              \
                          /(POSPART(1.0-CO2_BS)+KBUL50));                     \
         CO3_AP  = POSPART((DR_AP/DRR_AP)*1.0*(1.0-KOEDC*OED_AP)              \
                   *(1.0-STU_AP)*(1.0+CASENS));                               \
         CO3_MD  = POSPART((DR_MD/DRR_MD)*1.0*(1.0-KOEDC*OED_MD)              \
                   *(1.0-STU_MD)*(1.0+CASENS));                               \
         CO3_BS  = POSPART((DR_BS/DRR_BS)*1.0*(1.0-KOEDC*OED_BS)              \
                   *(1.0-STU_BS)*(1.0+CASENS));                               \
         EF_NATP = FAP*(SH0*CO3_AP/(1.0+(KAFT0/TH_AP)*(WSF_AP/WS0_AP-1.0))    \
                        - KBULGE*(1.0/TH_AP)*POSPART(1.0-CO3_AP)              \
                          /(POSPART(1.0-CO3_AP)+KBUL50))                      \
                 + FMD*(SH0*CO3_MD/(1.0+(KAFT0/TH_MD)*(WSF_MD/WS0_MD-1.0))    \
                        - KBULGE*(1.0/TH_MD)*POSPART(1.0-CO3_MD)              \
                          /(POSPART(1.0-CO3_MD)+KBUL50))                      \
                 + FBS*(SH0*CO3_BS/(1.0+(KAFT0/TH_BS)*(WSF_BS/WS0_BS-1.0))    \
                        - KBULGE*(1.0/TH_BS)*POSPART(1.0-CO3_BS)              \
                          /(POSPART(1.0-CO3_BS)+KBUL50));                     \
         SV_NOLV = SV0R*(EFX/EF0)*(LVEDV/LVEDV0)*(1.0 - MRF0)                 \
                   + SV0R*E_IMP_CO*IMPELLA;                                   \
  /* --- stasis and hypercoagulability --- */                                 \
         STASIS  = POSPART(SH0 - SH_AP)/SH0 * pow(LVEDV/LVEDV0, 1.5);         \
         HYPERC  = 1.0 + KHYPERC*POSPART(EPI/EPI0 - 1.0)/10.0;
$MAIN
// ---------------------------------------------------------------------------
//  RESTING REFERENCE STATE
//
//  Computed here, not typed.  At rest the three segments receive EPI0 plus
//  their own innervation-weighted interstitial noradrenaline; that fixes their
//  Gs coupling, hence their cAMP, PKA and contractile drive.  Every
//  normalisation in the algebra divides by these numbers, so CONT_x = 1,
//  CAL_x = 1, ATP_x = 1, GRK2 = 1 and RDN = 0 hold EXACTLY at t = 0 for any
//  parameter set - including a subject whose receptor gradient has been
//  changed through idata.
//
//  A worth-stating consequence: with RHO_AP/RHO_BS = 1.40 and the reciprocal
//  MIBG innervation gradient INN_AP/INN_BS = 0.52, the resting Gs coupling of
//  the three segments is almost identical.  The apex is NOT a weak segment at
//  rest.  It only becomes one when a blood-borne beta2 agonist arrives.
// ---------------------------------------------------------------------------
O20_AP = (EPI0 + AN2*NEI0*INN_AP)/(EPI0 + AN2*NEI0*INN_AP + KD_B2);
O20_MD = (EPI0 + AN2*NEI0*INN_MD)/(EPI0 + AN2*NEI0*INN_MD + KD_B2);
O20_BS = (EPI0 + AN2*NEI0*INN_BS)/(EPI0 + AN2*NEI0*INN_BS + KD_B2);
O10_AP = (AE1*EPI0 + NEI0*INN_AP)/(AE1*EPI0 + NEI0*INN_AP + KD_B1);
O10_MD = (AE1*EPI0 + NEI0*INN_MD)/(AE1*EPI0 + NEI0*INN_MD + KD_B1);
O10_BS = (AE1*EPI0 + NEI0*INN_BS)/(AE1*EPI0 + NEI0*INN_BS + KD_B1);
OCC0   = FAP*O20_AP + FMD*O20_MD + FBS*O20_BS;
GS0_AP = RHO_AP*(1.0 - FB2_AP)*O10_AP + RHO_AP*FB2_AP*O20_AP;
GS0_MD = RHO_MD*(1.0 - FB2_MD)*O10_MD + RHO_MD*FB2_MD*O20_MD;
GS0_BS = RHO_BS*(1.0 - FB2_BS)*O10_BS + RHO_BS*FB2_BS*O20_BS;
CAMP0_AP = (AC_BAS + EMAX_AC*GS0_AP/(GS0_AP + KM_GS))/KDEG_CA;
CAMP0_MD = (AC_BAS + EMAX_AC*GS0_MD/(GS0_MD + KM_GS))/KDEG_CA;
CAMP0_BS = (AC_BAS + EMAX_AC*GS0_BS/(GS0_BS + KM_GS))/KDEG_CA;
PKA0_AP  = CAMP0_AP*CAMP0_AP/(CAMP0_AP*CAMP0_AP + KPKA*KPKA);
PKA0_MD  = CAMP0_MD*CAMP0_MD/(CAMP0_MD*CAMP0_MD + KPKA*KPKA);
PKA0_BS  = CAMP0_BS*CAMP0_BS/(CAMP0_BS*CAMP0_BS + KPKA*KPKA);
DRR_AP   = pow(CAMP0_AP,NDR)/(pow(CAMP0_AP,NDR) + pow(KDR,NDR));
DRR_MD   = pow(CAMP0_MD,NDR)/(pow(CAMP0_MD,NDR) + pow(KDR,NDR));
DRR_BS   = pow(CAMP0_BS,NDR)/(pow(CAMP0_BS,NDR) + pow(KDR,NDR));
// Stroke volume is back-solved from the reference pressure, so MAP = MAP0 and
// therefore myocardial perfusion / demand = 1 at rest for any SVR0 or HR0.
SV0R = (MAP0 - CVP)*80.0/SVR0*1000.0/HR0/(1.0 - MRF0);
CO0R = SV0R*(1.0 - MRF0)*HR0/1000.0;
CAMP_AP_0 = CAMP0_AP;
CAMP_MD_0 = CAMP0_MD;
CAMP_BS_0 = CAMP0_BS;

$ODE
TTS_ALGEBRA

// ---------------------------------------------------------------- the surge
// The trigger amplitude is scaled by the oestradiol-dependent sympathetic gain.
// Postmenopausal E2 = 15 pg/mL gives a gain 2.9 x that of a premenopausal
// E2 = 100 pg/mL.  This is the only place sex enters the model.
double SUR   = (SOLVERTIME >= TRIG_T)
               ? AMP_TOT*SGAIN*exp(-(SOLVERTIME - TRIG_T)/TAU_SUR) : 0.0;
double SUR_E = FRAC_E*SUR;
double SUR_N = (1.0 - FRAC_E)*SUR;

// ------------------------------------------------- catecholamine kinetics
dxdt_EPI = SUR_E + KEPI_INF*RATE_EPI + KE_EPI*EPI0 - KE_EPI*EPI;
dxdt_NE  = SUR_N + KNEP_INF*RATE_NEP + KE_NE*NE0 - KE_NE*NE;
double PROD_NEI0 = KNEI_OUT*NEI0 + VMAX_UP1*NEI0/(KM_UP1 + NEI0);
// REFLEX limb.  Once forward flow falls, the baroreflex sustains sympathetic
// outflow for days after the primary surge has cleared.  This is what keeps the
// BASE hyperkinetic (and therefore keeps an LVOT gradient open) long after
// plasma adrenaline has normalised - and it is why breaking the loop with a
// short-acting beta-blocker helps an obstructed patient and hurts an
// unobstructed one.  It is a consequence of the flow deficit, not a schedule.
double CODEF = POSPART(1.0 - COX/fmax(CO0R, 0.1) - REFLDZ);
double REFL  = KREFL*CODEF/(CODEF + REFL50);
dxdt_NEI = PROD_NEI0 + FNEI*SUR_N + REFL + 0.35*KNEP_INF*RATE_NEP
           - KNEI_OUT*NEI - VMAX_UP1*NEI/(KM_UP1 + NEI);

// ------------------------------------------------------------- drug PK
dxdt_DOB  = RATE_DOB - (CL_DOB/V_DOB)*DOB;
dxdt_MIL  = RATE_MIL - (CL_MIL/V_MIL)*MIL;
dxdt_LEVC = RATE_LEV - (CL_LEV/V_LEV)*LEVC
            - (Q_LEV/V_LEV)*LEVC + (Q_LEV/VP_LEV)*LEVP;
dxdt_LEVP = (Q_LEV/V_LEV)*LEVC - (Q_LEV/VP_LEV)*LEVP;
dxdt_ORM  = FM_LEV*(CL_LEV/V_LEV)*LEVC - (CL_OR/V_OR)*ORM;
dxdt_ESM  = RATE_ESM - (CL_ESM/V_ESM)*ESM;
dxdt_METD = -KA_MET*METD;
dxdt_METC = F_MET*KA_MET*METD - (CL_MET/V_MET)*METC;
dxdt_CARD = -KA_CAR*CARD;
dxdt_CARC = F_CAR*KA_CAR*CARD - (CL_CAR/V_CAR)*CARC;
dxdt_PHE  = RATE_PHE - (CL_PHE/V_PHE)*PHE;
dxdt_RAMD = -KA_RAM*RAMD;
dxdt_RAMC = F_RAM*KA_RAM*RAMD - (CL_RAM/V_RAM)*RAMC;
dxdt_FURD = -KA_FUR*FURD;
dxdt_FURC = F_FUR*KA_FUR*FURD - (CL_FUR/V_FUR)*FURC;
dxdt_APXD = -KA_APX*APXD;
dxdt_APXC = F_APX*KA_APX*APXD - (CL_APX/V_APX)*APXC;

// -------------------------------------------------------- cAMP per segment
dxdt_CAMP_AP = AC_AP - KDEG*CAMP_AP;
dxdt_CAMP_MD = AC_MD - KDEG*CAMP_MD;
dxdt_CAMP_BS = AC_BS - KDEG*CAMP_BS;

// ---------------------------------------- THE SWITCH: one law, three densities
dxdt_PHOS_AP = KPH*SIG_AP*(1.0 - PHOS_AP) - KDEPH*PHOS_AP;
dxdt_PHOS_MD = KPH*SIG_MD*(1.0 - PHOS_MD) - KDEPH*PHOS_MD;
dxdt_PHOS_BS = KPH*SIG_BS*(1.0 - PHOS_BS) - KDEPH*PHOS_BS;

// ---------------------------------------------- GRK2 and receptor trafficking
double OCCBAR = FAP*O2_AP + FMD*O2_MD + FBS*O2_BS;
dxdt_GRK2 = KGRK_I*(1.0 + 4.0*POSPART(OCCBAR - OCC0)) - KGRK_O*GRK2;
dxdt_RDN  = KRDN_I*POSPART(OCCBAR - OCC0)*(RDNMAX - RDN) - KRDN_O*RDN;

// ----------------------------------------------------------- calcium loading
dxdt_CAL_AP = KCA_IN*(1.0 + ECA_PKA*(PKA_AP - PKA0_AP))
              - KCA_OUT*CAL_AP*(ATP_AP/(ATP_AP + KCA_ATP))/(1.0/(1.0+KCA_ATP));
dxdt_CAL_MD = KCA_IN*(1.0 + ECA_PKA*(PKA_MD - PKA0_MD))
              - KCA_OUT*CAL_MD*(ATP_MD/(ATP_MD + KCA_ATP))/(1.0/(1.0+KCA_ATP));
dxdt_CAL_BS = KCA_IN*(1.0 + ECA_PKA*(PKA_BS - PKA0_BS))
              - KCA_OUT*CAL_BS*(ATP_BS/(ATP_BS + KCA_ATP))/(1.0/(1.0+KCA_ATP));

// --------------------------------------------------------- energetic reserve
double SUPD_AP = (DEM_AP > 0.05) ? PERF_AP/DEM_AP : 20.0;
double SUPD_MD = (DEM_MD > 0.05) ? PERF_MD/DEM_MD : 20.0;
double SUPD_BS = (DEM_BS > 0.05) ? PERF_BS/DEM_BS : 20.0;
double DEM0_AP = KO2_CONT + KO2_HR + KO2_WS;
dxdt_ATP_AP = KATP*(fmin(1.0, CFR*SUPD_AP*DEM0_AP) - ATP_AP)
              - KMPTP*POSPART(CAL_AP - CATHR)*ATP_AP;
dxdt_ATP_MD = KATP*(fmin(1.0, CFR*SUPD_MD*DEM0_AP) - ATP_MD)
              - KMPTP*POSPART(CAL_MD - CATHR)*ATP_MD;
dxdt_ATP_BS = KATP*(fmin(1.0, CFR*SUPD_BS*DEM0_AP) - ATP_BS)
              - KMPTP*POSPART(CAL_BS - CATHR)*ATP_BS;

// ------------------------------------------------------------------- oedema
dxdt_OED_AP = KOED*(PRM_AP + WINF_OED*POSPART(IL6/2.5 - 1.0))*(1.0 - OED_AP)
              - KOEDR*OED_AP;
dxdt_OED_MD = KOED*(PRM_MD + WINF_OED*POSPART(IL6/2.5 - 1.0))*(1.0 - OED_MD)
              - KOEDR*OED_MD;
dxdt_OED_BS = KOED*(PRM_BS + WINF_OED*POSPART(IL6/2.5 - 1.0))*(1.0 - OED_BS)
              - KOEDR*OED_BS;

// ----------------------------------------------------------------- stunning
double BPROT = 1.0 - WBARR*PHOS_AP;
dxdt_STU_AP = KSTUN*(pow(POSPART(CAL_AP - CATHR), 2.0)
                     + WSTR_ST*STR_AP*STR_AP)*BPROT*(1.0 - STU_AP)
              - KSTUNR*STU_AP;
dxdt_STU_MD = KSTUN*(pow(POSPART(CAL_MD - CATHR), 2.0)
                     + WSTR_ST*STR_MD*STR_MD)*(1.0 - WBARR*PHOS_MD)
              *(1.0 - STU_MD) - KSTUNR*STU_MD;
dxdt_STU_BS = KSTUN*(pow(POSPART(CAL_BS - CATHR), 2.0)
                     + WSTR_ST*STR_BS*STR_BS)*(1.0 - WBARR*PHOS_BS)
              *(1.0 - STU_BS) - KSTUNR*STU_BS;

// ------------------------------------------------------- necrosis and troponin
double NECDRV = FAP*POSPART(NECTHR - ATP_AP) + FMD*POSPART(NECTHR - ATP_MD)
                + FBS*POSPART(NECTHR - ATP_BS);
dxdt_NECR = KNECR*NECDRV;
dxdt_TNI  = KLEAK*LEAKF + KNEC_TNI*KNECR*NECDRV
            + KTNI_OUT*0.010 - KTNI_OUT*TNI;

// -------------------------------------------------------------- inflammation
dxdt_IL6 = KIL6_IN*(LEAKF + 40.0*NECDRV) + KIL6_OUT*2.5 - KIL6_OUT*IL6;
dxdt_CRP = KCRP_IN*IL6 - KCRP_OUT*CRP;

// ------------------------------------------------------------- NT-proBNP
// NT-proBNP transcription follows END-DIASTOLIC wall stress (filling pressure
// x radius / thickness).  That is why it can be very high while troponin is
// barely raised: the two markers report different physical quantities.
double WSBAR = RAD*fmax(PCWPX/PCWP0, 0.2);
double STRS  = pow(fmax(WSBAR, 0.2), NBNP);
dxdt_BSYN  = KBSYN*(STRS - BSYN);
dxdt_NTBNP = KBNP_O*(NTBNP0*(1.0 + GBNP*POSPART(BSYN - 1.0)) - NTBNP);

// ------------------------------------------------------------ neurohormonal
double RENDRV = 1.0 + ERENB1*POSPART(O1_BS - O10_BS)
                + EREN_MAP*POSPART(1.0 - MAPX/MAP0);
dxdt_ANGII = KREN_I*RENDRV*(1.0 - ACEINH) - KANG_O*ANGII;
dxdt_ALDO  = KALD_I*ANGII - KALD_O*ALDO;

// ------------------------------------------------------------ haemodynamics
// Oestradiol-dependent NO tone BUFFERS the alpha1 response rather than
// shifting the resting SVR: a constant offset would only move the baseline and
// would contribute nothing to the dynamics.
double SVRT = SVR0*(1.0 + EA1_SVR*(1.0 - ENO*NOF)*(OA1 - OA1_0)/fmax(OA1_0,1e-6)
                    + EAII*(ANGII/ANGII0 - 1.0)
                    - KATP_LEV*CASENS - E_IABP*IABP);
dxdt_SVR = KSVR*(fmin(fmax(SVRT, 350.0), 3000.0) - SVR);
double HRT = HR0*(1.0 + EHR_B1*POSPART(O1_BS - O10_BS)/HR0)
             + KBARO*POSPART(MAP0 - MAPX);
dxdt_HR = KHR*(fmin(fmax(HRT, 40.0), 190.0) - HR);
// The natriuretic response SATURATES.  Left proportional to NT-proBNP, a
// 20-fold rise removes 1.9 L/h and empties the circulation - which is not a
// subtle calibration error but a structural one.
double BNPF = POSPART(NTBNP/NTBNP0 - 1.0);
dxdt_VOL = KVIN*(1.0 + EALD*(ALDO - 1.0))
           - KVOUT*VOL*(1.0 + EBNP*BNPF/(BNPF + KBNPV)) - KFUR_V*FUREF;
double LVEDVT = LVEDV0*(1.0 + KDILAT*POSPART(1.0 - SH_AP/SH0)
                        + KVPRE*(VOL/VOL0 - 1.0) - E_IMP_UNL*IMPELLA);
dxdt_LVEDV = KLV*(fmax(LVEDVT, 55.0) - LVEDV);

// -------------------------------------------------------- electrophysiology
dxdt_KPL = KK_IN + KSUP - KK_OUT*KPL - KK_FUR*FUREF - KK_CAT*POSPART(O2_BS - O20_BS);
double OEDBAR = FAP*OED_AP + FMD*OED_MD + FBS*OED_BS;
double CALBAR = FAP*CAL_AP + FMD*CAL_MD + FBS*CAL_BS;
// Repolarisation follows tissue oedema AND the apex-to-base dispersion of
// contractile state (regional repolarisation heterogeneity).  Both feed a slow
// first-order remodelling state, which is why the QTc peak LAGS the event
// instead of coinciding with the catecholamine surge.
double QTT = QTC0 + AQT_OED*OEDBAR + AQT_BAL*POSPART(BALL)
             + AQT_CA*POSPART(CALBAR - 1.0)
             + AQT_D*QTDRUG + AQT_K*POSPART(4.0 - KPL);
dxdt_QTC  = KQT*(QTT - QTC);
// The torsade hazard must be ZERO at a normal QTc, not merely small: a
// constant term integrated over 90 days gave a 4.9% cumulative probability in
// the no-trigger control arm, which is not a calibration error but a wrong
// functional form.  Subtracting 1 makes the hazard vanish below 450 ms, and the
// coefficient is then set so that the reference arm accumulates the reported
// 2-5% incidence rather than the 17% the original form produced.
dxdt_HTDP = H0_TDP*(exp(BQT*POSPART(QTC - 450.0)/10.0) - 1.0);

// ---------------------------------------------------------------- thrombus
dxdt_THR  = KTF*STASIS*HYPERC*(1.0 - 0.75*ANTIC) - KTL*(1.0 + KAC_LYS*ANTIC)*THR;
dxdt_HEMB = H0_EMB*THR*(1.0 - 0.70*ANTIC);

// --------------------------------------------------------------- recurrence
double PHBAR = FAP*PHOS_AP + FMD*PHOS_MD + FBS*PHOS_BS;
dxdt_HREC = H0_REC*(1.0 + W_PHREC*PHBAR)*(1.0 - E_ACEREC*ACEINH)
            *(1.0 - E_BBREC*(1.0 - INH1));

// ------------------------------------------------- inference / bookkeeping
dxdt_AUCEPI = EPI;
dxdt_CHDGI  = POSPART(EF_NOGI - EFX);
dxdt_CHDOED = POSPART(EF_NOED - EFX);
dxdt_CHDATP = POSPART(EF_NATP - EFX);
dxdt_CHDLVO = POSPART(SV_NOLV - SVX);
dxdt_CEFDEF = POSPART(EF0 - EFX);

$TABLE
TTS_ALGEBRA

double LVEF   = 100.0*EFX;
double SV     = SVX;
double CO     = COX;
double MAP    = MAPX;
double PCWP   = PCWPX;
double GRAD   = GRADX;
double WMSI   = 1.0 + 2.0*POSPART(1.0 - SH_AP/SH0)*FAP/0.28*0.9
                + 1.2*POSPART(1.0 - SH_MD/SH0) + 0.6*POSPART(1.0 - SH_BS/SH0);
double GLS    = -18.5*(EFX/EF0);
double IMR    = IMRX;
double BNPTNI = NTBNP/fmax(TNI, 1e-4);
double PTDP   = 1.0 - exp(-HTDP);
double PEMB   = 1.0 - exp(-HEMB);
double PREC   = 1.0 - exp(-HREC);
double SHR_GI = (CEFDEF > 0.5) ? CHDGI/CEFDEF : 0.0;
double SHR_OE = (CEFDEF > 0.5) ? CHDOED/CEFDEF : 0.0;
double SHR_AT = (CEFDEF > 0.5) ? CHDATP/CEFDEF : 0.0;
double CQERR  = fabs(CQ_AP - CAMP_AP)/fmax(CAMP_AP, 1e-6);
double VARIANT = (BALL > 0.25) ? 1.0 : ((BALL < -0.25) ? -1.0 : 0.0);
double DRIVE  = FAP*CONT_AP + FMD*CONT_MD + FBS*CONT_BS;

$CAPTURE @annotated
LVEF   : LV ejection fraction (%)
SV     : Stroke volume (mL)
CO     : Cardiac output (L/min)
MAP    : Mean arterial pressure (mmHg)
PCWP   : Pulmonary capillary wedge pressure (mmHg)
GRAD   : LVOT pressure gradient (mmHg)
BALL   : Ballooning index = basal minus apical shortening (-)
VARIANT: Morphological variant (+1 apical, 0 none, -1 reverse)
WMSI   : Wall motion score index (-)
GLS    : Global longitudinal strain (%)
SH_AP  : Apical fractional shortening (-)
SH_MD  : Mid fractional shortening (-)
SH_BS  : Basal fractional shortening (-)
CONT_AP: Apical contractility (-)
CONT_BS: Basal contractility (-)
DRIVE  : Mass-weighted contractile drive (-)
O2_AP  : Apical beta2 occupancy (-)
O2_BS  : Basal beta2 occupancy (-)
GI_AP  : Apical Gi coupling (-)
GI_BS  : Basal Gi coupling (-)
GS_AP  : Apical Gs coupling (-)
GS_BS  : Basal Gs coupling (-)
CQ_AP  : Apical quasi-steady-state cAMP (-)
CQERR  : Relative QSS error in apical cAMP (-)
IMR    : Microvascular resistance index (-)
BNPTNI : NT-proBNP to troponin ratio (pg/mL per ng/mL)
FFLOSS : Fractional forward-flow loss from LVOT obstruction (-)
CASENS : Myofilament calcium sensitisation (-)
ACEINH : Fractional ACE inhibition (-)
ANTIC  : Fractional anticoagulation (-)
PTDP   : Cumulative probability of torsade (-)
PEMB   : Cumulative probability of embolism (-)
PREC   : Cumulative probability of recurrence (-)
EF_NOGI: Counterfactual EF with no Gi coupling (-)
EF_NOED: Counterfactual EF with no oedema (-)
EF_NATP: Counterfactual EF with normal energetics (-)
SV_NOLV: Counterfactual stroke volume with no gradient (mL)
SHR_GI : Share of the EF deficit attributable to Gi (-; 0 below 0.5 EF-hours)
SHR_OE : Share attributable to oedema (-)
SHR_AT : Share attributable to energetics (-)
STASIS : Apical stasis index (-)

// Between-subject variability is applied through `idata` in `population()`
// (receptor density, septal geometry, surge amplitude, adrenaline fraction and
// oestradiol are all named parameters, so they do not need ETAs).
$OMEGA 0

$SIGMA 0
)---"

## -----------------------------------------------------------------------------
##  BUILD
## -----------------------------------------------------------------------------
mod <- mcode_cache("tts_qsp", tts_code, atol = 1e-8, rtol = 1e-8,
                   maxsteps = 2000000)

## =============================================================================
##  PHENOTYPES
##
##  A phenotype here is NOT a diagnosis.  It is (i) a receptor field, (ii) a
##  surge amplitude and adrenaline fraction, (iii) an oestradiol level and
##  (iv) a septal geometry.  Everything clinical - LVEF, ballooning index,
##  gradient, troponin - is an OUTPUT of running those four inputs.
## =============================================================================
phenotypes <- list(
  control = list(AMP_TOT = 0,    FRAC_E = 0.75, TAU_SUR = 1.6,  E2 = 15,
                 SEPT = 0.35, label = "No trigger (control)"),
  emo     = list(AMP_TOT = 250, FRAC_E = 0.75, TAU_SUR = 1.6,  E2 = 15,
                 SEPT = 0.35, label = "Emotional trigger, postmenopausal"),
  emo_pre = list(AMP_TOT = 250, FRAC_E = 0.75, TAU_SUR = 1.6,  E2 = 100,
                 SEPT = 0.35, label = "Emotional trigger, premenopausal"),
  phys    = list(AMP_TOT = 300, FRAC_E = 0.30, TAU_SUR = 12.0, E2 = 15,
                 SEPT = 0.35, label = "Physical trigger (NE-dominant)"),
  ## NOTE ON THE SAH ARM.  This model has no representation of death, so any
  ## arm that drives cardiac output below about 1.5 L/min leaves the region
  ## where its haemodynamic algebra means anything (MAP collapses onto the CVP
  ## floor and the integrator fails).  AMP_TOT = 600 with TAU_SUR = 48 h did
  ## exactly that.  The neurologic arm is therefore set just inside that
  ## boundary - severe, survivable - and the boundary itself is reported in the
  ## README as a limit of validity rather than hidden behind a solver tolerance.
  neuro   = list(AMP_TOT = 420, FRAC_E = 0.55, TAU_SUR = 24.0, E2 = 15,
                 SEPT = 0.35, label = "Neurologic trigger (SAH)"),
  pheo    = list(AMP_TOT = 160, FRAC_E = 0.92, TAU_SUR = 72.0, E2 = 15,
                 SEPT = 0.35, label = "Adrenaline-secreting phaeochromocytoma"),
  obstruct= list(AMP_TOT = 250, FRAC_E = 0.75, TAU_SUR = 1.6,  E2 = 15,
                 SEPT = 1.00, label = "Emotional trigger WITH obstructive septal geometry"),
  ne_only = list(AMP_TOT = 250, FRAC_E = 0.00, TAU_SUR = 1.6,  E2 = 15,
                 SEPT = 0.35, label = "Pure noradrenaline surge (falsification test)"),
  reverse = list(AMP_TOT = 250, FRAC_E = 0.75, TAU_SUR = 1.6,  E2 = 15,
                 SEPT = 0.35, RHO_AP = 1.00, RHO_BS = 1.40,
                 FB2_AP = 0.24, FB2_BS = 0.42, INN_AP = 1.20, INN_BS = 0.62,
                 label = "Reverse (basal) TTS - requires an INVERTED gradient")
)

pheno_param <- function(ph) {
  p <- phenotypes[[ph]]
  p$label <- NULL
  p
}

## =============================================================================
##  BASELINE
##
##  Every scenario starts from a 14-day burn-in with NO trigger, so the
##  baseline is derived rather than typed.  A hand-set initial condition would
##  let the untreated trajectory be asserted.
## =============================================================================
baseline <- function(ph = "emo", burn = 336) {
  p  <- pheno_param(ph)
  p$AMP_TOT <- 0
  m  <- param(mod, p)
  o  <- mrgsim_df(zero_re(m), end = burn, delta = burn, recsort = 3)
  last <- o[nrow(o), ]
  st <- as.list(last[intersect(names(last), names(as.list(init(mod))))])
  list(mod = init(param(mod, pheno_param(ph)), st), init = st,
       param = pheno_param(ph), label = phenotypes[[ph]]$label)
}

## =============================================================================
##  EVENT HELPERS  (all times in HOURS from the trigger)
## =============================================================================
## Infusions are driven by RATE_* parameters rather than by dosing records, so
## an arm is applied by SPLITTING the simulation at each parameter change and
## carrying the state across the join.  `run_arm()` does that; oral drugs use
## ordinary mrgsolve events into their depot compartments.
run_arm <- function(b, changes = list(), end = 2160, delta = 0.25) {
  ## changes: list of list(t=, par=list(...)) applied in increasing time order
  if (!length(changes)) {
    return(mrgsim_df(zero_re(b$mod), end = end, delta = delta, recsort = 3))
  }
  ts   <- vapply(changes, function(x) x$t, numeric(1))
  ord  <- order(ts); changes <- changes[ord]; ts <- ts[ord]
  segs <- c(0, ts, end)
  m    <- zero_re(b$mod)
  outs <- list(); st <- b$init
  for (i in seq_len(length(segs) - 1)) {
    if (i > 1 && length(changes[[i - 1]]$par))
      m <- param(m, changes[[i - 1]]$par)
    m2 <- init(m, st)
    o  <- mrgsim_df(m2, start = 0, end = segs[i + 1] - segs[i], delta = delta,
                    recsort = 3)
    o$time <- o$time + segs[i]
    st <- as.list(o[nrow(o), intersect(names(o), names(as.list(init(mod))))])
    outs[[i]] <- if (i == 1) o else o[-1, ]
  }
  do.call(rbind, outs)
}

run_oral <- function(b, events, end = 2160, delta = 0.5, par = list()) {
  m <- zero_re(b$mod)
  if (length(par)) m <- param(m, par)
  mrgsim_df(m, events = events, end = end, delta = delta, recsort = 3)
}

## =============================================================================
##  SCENARIOS
##
##  28 arms.  They are built as MATCHED PAIRS wherever a claim is being tested,
##  because an arm on its own only shows that the model produces a number.
## =============================================================================
scenarios <- list(

  ## --- 1-6 : the disease itself, by trigger -------------------------------
  S01 = list(ph = "control",  ch = list(), label = "S01 No trigger (control)"),
  S02 = list(ph = "emo",      ch = list(), label = "S02 Emotional trigger, untreated"),
  S03 = list(ph = "phys",     ch = list(), label = "S03 Physical trigger (NE-dominant)"),
  S04 = list(ph = "neuro",    ch = list(), label = "S04 Neurologic trigger (SAH)"),
  S05 = list(ph = "pheo",     ch = list(), label = "S05 Adrenaline-secreting phaeo"),
  S06 = list(ph = "ne_only",  ch = list(), label = "S06 Pure NE surge (falsification 1)"),

  ## --- 7 : the mechanistic control ---------------------------------------
  S07 = list(ph = "emo", ch = list(list(t = 0.01, par = list(PTX = 1))),
             label = "S07 Emotional trigger + pertussis toxin (falsification 2)"),

  ## --- 8-11 : inotropes, the routing claim -------------------------------
  S08 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_DOB = 105)),
                       list(t = 52, par = list(RATE_DOB = 0))),
             label = "S08 Dobutamine 5 ug/kg/min, h4-52"),
  S09 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_LEV = 0.42)),
                       list(t = 28, par = list(RATE_LEV = 0))),
             label = "S09 Levosimendan 0.1 ug/kg/min x 24 h from h4"),
  S10 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_MIL = 3.0)),
                       list(t = 52, par = list(RATE_MIL = 0))),
             label = "S10 Milrinone 0.5 ug/kg/min, h4-52 (falsification 3)"),
  S11 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_EPI = 5)),
                       list(t = 28, par = list(RATE_EPI = 0))),
             label = "S11 Adrenaline 5 ug/min, h4-28"),

  ## --- 12-13 : noradrenaline and dopamine-like pressor support -----------
  ## 8 ug/min drove this arm past the limit of validity described in the README
  ## (cardiac output below ~1.5 L/min, where the haemodynamic algebra stops
  ## meaning anything); 5 ug/min sits inside it and still shows the harm.
  S12 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_NEP = 5)),
                       list(t = 52, par = list(RATE_NEP = 0))),
             label = "S12 Noradrenaline 5 ug/min, h4-52"),
  S13 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_PHE = 30)),
                       list(t = 52, par = list(RATE_PHE = 0))),
             label = "S13 Phenylephrine 30 mg/h, h4-52"),

  ## --- 14-17 : beta-blockade, the conditional-sign pair ------------------
  S14 = list(ph = "emo",
             ch = list(list(t = 4, par = list(RATE_ESM = 210)),
                       list(t = 52, par = list(RATE_ESM = 0))),
             label = "S14 Esmolol, NO obstruction (expected harm)"),
  S15 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(RATE_ESM = 210)),
                       list(t = 52, par = list(RATE_ESM = 0))),
             label = "S15 Esmolol, WITH obstruction (expected benefit)"),
  S16 = list(ph = "obstruct", ch = list(),
             label = "S16 Obstructive TTS, untreated (pair for S15)"),
  S17 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(RATE_PHE = 30)),
                       list(t = 52, par = list(RATE_PHE = 0))),
             label = "S17 Phenylephrine for obstruction"),

  ## --- 18-19 : mechanical support, same gradient term --------------------
  S18 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(IABP = 1)),
                       list(t = 72, par = list(IABP = 0))),
             label = "S18 IABP in obstructive TTS (expected harm)"),
  S19 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(IMPELLA = 1)),
                       list(t = 72, par = list(IMPELLA = 0))),
             label = "S19 Impella in obstructive TTS"),

  ## --- 20-21 : dobutamine in the obstructive phenotype -------------------
  S20 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(RATE_DOB = 105)),
                       list(t = 52, par = list(RATE_DOB = 0))),
             label = "S20 Dobutamine in obstructive TTS"),
  S21 = list(ph = "obstruct",
             ch = list(list(t = 4, par = list(RATE_LEV = 0.42)),
                       list(t = 28, par = list(RATE_LEV = 0))),
             label = "S21 Levosimendan in obstructive TTS"),

  ## --- 22-23 : sex / oestradiol pair -------------------------------------
  S22 = list(ph = "emo_pre", ch = list(),
             label = "S22 Identical trigger, premenopausal (pair for S02)"),
  S23 = list(ph = "emo", ch = list(list(t = 0.01, par = list(SEDATE = 1))),
             label = "S23 Anxiolysis at presentation"),

  ## --- 24 : QT-prolonging comedication -----------------------------------
  S24 = list(ph = "emo", ch = list(list(t = 0.01, par = list(QTDRUG = 1))),
             label = "S24 QT-prolonging comedication"),

  ## --- 25 : reverse TTS (inverted gradient, stated as a prediction) ------
  S25 = list(ph = "reverse", ch = list(),
             label = "S25 Inverted receptor gradient (reverse TTS prediction)"),

  ## --- 26-28 : oral chronic arms (built with `oral()` in run_all) --------
  S26 = list(ph = "emo", ch = list(), oral = "ramipril",
             label = "S26 Ramipril 5 mg daily from day 1"),
  S27 = list(ph = "emo", ch = list(), oral = "metoprolol",
             label = "S27 Metoprolol 50 mg bd from day 1"),
  S28 = list(ph = "emo", ch = list(), oral = "apixaban",
             label = "S28 Apixaban 5 mg bd from day 1")
)

oral_events <- list(
  ramipril   = ev(amt = 5,  cmt = "RAMD", ii = 24, addl = 364, time = 24),
  metoprolol = ev(amt = 50, cmt = "METD", ii = 12, addl = 729, time = 24),
  apixaban   = ev(amt = 5,  cmt = "APXD", ii = 12, addl = 179, time = 24),
  furosemide = ev(amt = 40, cmt = "FURD", ii = 12, addl = 13,  time = 12),
  carvedilol = ev(amt = 12.5, cmt = "CARD", ii = 12, addl = 729, time = 24)
)

run_scenario <- function(id, end = 2160, delta = 0.25) {
  s <- scenarios[[id]]
  b <- baseline(s$ph)
  if (!is.null(s$oral)) {
    o <- run_oral(b, oral_events[[s$oral]], end = end, delta = delta)
  } else {
    o <- run_arm(b, s$ch, end = end, delta = delta)
  }
  o$scenario <- id
  o$label    <- s$label
  o
}

run_all <- function(ids = names(scenarios), end = 2160, delta = 0.25) {
  do.call(rbind, lapply(ids, run_scenario, end = end, delta = delta))
}

## =============================================================================
##  SUMMARY
## =============================================================================
summarise_arm <- function(o) {
  d <- function(t) o[which.min(abs(o$time - t)), ]
  data.frame(
    scenario   = o$scenario[1],
    label      = o$label[1],
    LVEF_nadir = round(min(o$LVEF), 1),
    t_nadir_h  = round(o$time[which.min(o$LVEF)], 1),
    LVEF_d1    = round(d(24)$LVEF, 1),
    LVEF_d7    = round(d(168)$LVEF, 1),
    LVEF_d30   = round(d(720)$LVEF, 1),
    LVEF_d90   = round(d(2159)$LVEF, 1),
    BALL_peak  = round(max(o$BALL), 3),
    BALL_min   = round(min(o$BALL), 3),
    SH_AP_min  = round(min(o$SH_AP), 3),
    SH_BS_max  = round(max(o$SH_BS), 3),
    GRAD_peak  = round(max(o$GRAD), 1),
    CO_min     = round(min(o$CO), 2),
    MAP_min    = round(min(o$MAP), 1),
    PCWP_max   = round(max(o$PCWP), 1),
    TNI_peak   = round(max(o$TNI), 3),
    NTBNP_peak = round(max(o$NTBNP), 0),
    BNPTNI     = round(max(o$NTBNP)/max(max(o$TNI), 1e-4), 0),
    QTC_peak   = round(max(o$QTC), 0),
    t_QTC_d    = round(o$time[which.max(o$QTC)]/24, 1),
    PTDP_pct   = round(100*max(o$PTDP), 2),
    THR_peak   = round(max(o$THR), 4),
    PEMB_pct   = round(100*max(o$PEMB), 2),
    PHOSAP_max = round(max(o$GI_AP), 3),
    SHR_GI     = round(tail(o$SHR_GI, 1), 3),
    SHR_OE     = round(tail(o$SHR_OE, 1), 3),
    SHR_AT     = round(tail(o$SHR_AT, 1), 3),
    QSSERR     = signif(max(o$CQERR), 2),
    row.names  = NULL
  )
}

summarise_all <- function(all) {
  do.call(rbind, lapply(split(all, all$scenario), summarise_arm))
}

## =============================================================================
##  CALIBRATION
##
##  Six parameters are fitted by Nelder-Mead to eight published anchors.  The
##  anchors are all from the acute apical phenotype (InterTAK / Templin 2015,
##  Ghadri 2018 consensus, Eitel 2011 CMR, Nef 2007, Madhavan 2009):
##
##      LVEF nadir                       38 %
##      LVEF at day 30                   58 %
##      Ballooning index peak            0.75  (basal minus apical shortening)
##      hs-cTnI peak                     6.0  ng/mL   (~7-10 x URL)
##      NT-proBNP peak                   4200 pg/mL
##      NT-proBNP:troponin ratio         700
##      QTc peak                         488  ms
##      Time to QTc peak                 2.5  d
##
##  Reported objective and per-anchor residuals come from `calibrate()`.
## =============================================================================
## The LVEF anchor is the ADMISSION echocardiogram (24 h), not the instantaneous
## minimum of the trajectory.  Anchoring min(LVEF) makes the objective depend on
## output resolution: at delta = 1 h the reported minimum is 38%, at delta =
## 0.25 h it is 33%, because the model dips transiently during the catecholamine
## peak in the first two hours - a moment no patient is scanned at.  Reporting a
## resolution-dependent quantity as a fitted value would be a silent artefact.
anchors <- data.frame(
  name   = c("LVEF_h24", "LVEF_d30", "BALL_peak", "TNI_peak",
             "NTBNP_peak", "BNPTNI", "QTC_peak", "t_QTC_d"),
  target = c(38, 58, 0.75, 6.0, 4200, 700, 488, 2.5),
  weight = c(3, 2, 2, 1.5, 1, 1, 1.5, 1)
)

## Eight parameters are estimated.  The oedema kinetics (KOED, KOEDR), the
## stunning kinetics (KSTUN, KSTUNR) and the oedema-to-QTc slope (AQT_OED) are
## held FIXED at their literature-derived values: an unconstrained fit trades
## KOED against AQT_OED without limit (their product is what the QTc anchor
## sees), which produces an excellent objective and an uninterpretable model.
## KREFL is also FIXED, and the reason is worth recording.  Left free, the
## optimiser drove it to 70 000 nmol/L/h.  Because the reflex saturates
## (CODEF/(CODEF+REFL50)), any value that large turns a graded baroreflex into a
## near-binary switch that fires whenever cardiac output falls 3%: interstitial
## noradrenaline leaves the physiological range, beta1 occupancy saturates, and
## the SIZE of the catecholamine surge stops mattering.  The visible symptom was
## that the premenopausal arm (a 2.9-fold smaller surge) stopped being milder -
## the oestradiol effect was washed out by a reflex that no longer cared how
## large the primary insult was.  KREFL is therefore pinned at 3300, the value
## that keeps peak interstitial noradrenaline near 60 nmol/L (about 8 x resting),
## and the remaining seven parameters are refitted around it.  A better
## objective bought with a non-physiological reflex is not a better model.
fit_pars <- c(GAM_GI = 9.97164, WB2D = 2.63360, KBULGE = 0.15370,
              KLEAK = 4.85864, GBNP = 7.03797, AQT_BAL = 136.375,
              KWSI = 0.14843)

objective <- function(theta, verbose = FALSE) {
  p <- as.list(exp(theta)); names(p) <- names(fit_pars)
  b <- try(baseline("emo"), silent = TRUE)
  if (inherits(b, "try-error")) return(1e6)
  b$mod <- param(b$mod, p)
  o <- try(mrgsim_df(zero_re(b$mod), end = 1200, delta = 0.25, recsort = 3),
           silent = TRUE)
  if (inherits(o, "try-error") || any(!is.finite(o$LVEF))) return(1e6)
  pred <- c(LVEF_h24   = o$LVEF[which.min(abs(o$time - 24))],
            LVEF_d30   = o$LVEF[which.min(abs(o$time - 720))],
            BALL_peak  = max(o$BALL),
            TNI_peak   = max(o$TNI),
            NTBNP_peak = max(o$NTBNP),
            BNPTNI     = max(o$NTBNP)/max(max(o$TNI), 1e-4),
            QTC_peak   = max(o$QTC),
            t_QTC_d    = o$time[which.max(o$QTC)]/24)
  res <- (pred[anchors$name] - anchors$target)/anchors$target
  if (verbose) print(data.frame(anchors, pred = round(pred[anchors$name], 3),
                                rel = round(res, 3)))
  sum(anchors$weight * res^2)/sum(anchors$weight)
}

calibrate <- function(maxit = 1500, restarts = 3) {
  th <- log(fit_pars); op <- NULL
  for (i in seq_len(restarts)) {
    op <- optim(th, objective, method = "Nelder-Mead",
                control = list(maxit = maxit, reltol = 1e-8))
    th <- op$par
  }
  fitted <- exp(op$par); names(fitted) <- names(fit_pars)
  list(par = fitted, value = op$value, convergence = op$convergence)
}

## =============================================================================
##  INFERENCE 1 - MATCHED CONTRACTILE DRIVE
##
##  Dobutamine and levosimendan are titrated so that the mass-weighted
##  contractile drive added at hour 8 is the SAME.  Only then does the
##  comparison isolate the ROUTE rather than the potency.
## =============================================================================
drive_at <- function(par, t = 8) {
  b <- baseline("emo")
  o <- run_arm(b, list(list(t = 4, par = par), list(t = 52, par = list())),
               end = 24, delta = 0.25)
  o$DRIVE[which.min(abs(o$time - t))]
}

matched_drive <- function(target = NULL) {
  base_drive <- drive_at(list())
  if (is.null(target)) target <- base_drive * 1.18
  f <- function(rate, which) {
    p <- if (which == "dob") list(RATE_DOB = rate) else list(RATE_LEV = rate)
    drive_at(p) - target
  }
  rd <- try(uniroot(f, c(1, 900), which = "dob")$root, silent = TRUE)
  rl <- try(uniroot(f, c(0.005, 6), which = "lev")$root, silent = TRUE)
  out <- list(target_drive = target, base_drive = base_drive,
              rate_dob = if (inherits(rd, "try-error")) NA else rd,
              rate_lev = if (inherits(rl, "try-error")) NA else rl)
  b <- baseline("emo")
  if (!is.na(out$rate_dob)) {
    od <- run_arm(b, list(list(t = 4, par = list(RATE_DOB = out$rate_dob)),
                          list(t = 52, par = list(RATE_DOB = 0))), end = 720)
    out$dob <- summarise_arm(transform(od, scenario = "matched_dob",
                                       label = "matched dobutamine"))
  }
  if (!is.na(out$rate_lev)) {
    ol <- run_arm(b, list(list(t = 4, par = list(RATE_LEV = out$rate_lev)),
                          list(t = 28, par = list(RATE_LEV = 0))), end = 720)
    out$lev <- summarise_arm(transform(ol, scenario = "matched_lev",
                                       label = "matched levosimendan"))
  }
  out
}

## =============================================================================
##  INFERENCE 2 - WHERE DOES THE THRESHOLD SIT?
##
##  The claim is that the apex crosses the Gs->Gi threshold at a LOWER plasma
##  adrenaline than the base.  That is a model OUTPUT and it is measured here
##  by ramping a steady adrenaline infusion and recording the concentration at
##  which each segment's Gi coupling exceeds 20% of its beta2 pool.
## =============================================================================
threshold_scan <- function(rates = c(0.5, 1, 1.5, 2, 3, 4, 6, 8, 12)) {
  b <- baseline("control")
  do.call(rbind, lapply(rates, function(r) {
    o <- run_arm(b, list(list(t = 1, par = list(RATE_EPI = r))),
                 end = 12, delta = 0.25)
    l <- o[nrow(o), ]
    data.frame(rate_ug_min = r, EPI_nM = round(l$EPI, 2),
               occ2_AP = round(l$O2_AP, 3), occ2_BS = round(l$O2_BS, 3),
               GI_AP = round(l$GI_AP, 3), GI_BS = round(l$GI_BS, 3),
               fracGi_AP = round(l$GI_AP/(1.40*0.42*l$O2_AP + 1e-9), 3),
               fracGi_BS = round(l$GI_BS/(1.00*0.24*l$O2_BS + 1e-9), 3),
               LVEF = round(l$LVEF, 1), BALL = round(l$BALL, 3),
               row.names = NULL)
  }))
}

## =============================================================================
##  INFERENCE 3 - THE THREE FALSIFICATION TESTS
## =============================================================================
falsification <- function() {
  s <- function(id) summarise_arm(run_scenario(id, end = 720))
  rbind(
    transform(s("S02"), test = "reference: adrenaline-rich surge"),
    transform(s("S06"), test = "F1: NE-only must NOT balloon"),
    transform(s("S07"), test = "F2: pertussis toxin must abolish it"),
    transform(s("S08"), test = "F3a: dobutamine (receptor-level cAMP)"),
    transform(s("S10"), test = "F3b: milrinone (post-receptor cAMP)")
  )[, c("test", "LVEF_nadir", "BALL_peak", "SH_AP_min", "SH_BS_max",
        "PHOSAP_max", "SHR_GI")]
}

## =============================================================================
##  INFERENCE 4 - THE CONDITIONAL SIGN OF BETA-BLOCKADE
## =============================================================================
bblocker_sign <- function() {
  s <- function(id) summarise_arm(run_scenario(id, end = 720))
  a <- s("S02"); b <- s("S14"); c <- s("S16"); d <- s("S15")
  data.frame(
    setting = c("no obstruction", "no obstruction",
                "obstruction", "obstruction"),
    arm     = c("untreated", "esmolol", "untreated", "esmolol"),
    LVEF_nadir = c(a$LVEF_nadir, b$LVEF_nadir, c$LVEF_nadir, d$LVEF_nadir),
    GRAD_peak  = c(a$GRAD_peak, b$GRAD_peak, c$GRAD_peak, d$GRAD_peak),
    CO_min     = c(a$CO_min, b$CO_min, c$CO_min, d$CO_min),
    MAP_min    = c(a$MAP_min, b$MAP_min, c$MAP_min, d$MAP_min),
    row.names  = NULL
  )
}

## =============================================================================
##  INFERENCE 4b - IS THERE ANY ESMOLOL DOSE THAT HELPS AN OBSTRUCTED PATIENT?
##
##  A single dose cannot answer a question about a CONDITIONAL SIGN: if the
##  chosen dose is past the point where the inotropic cost exceeds the gradient
##  benefit, the arm looks like harm in both phenotypes and the conditionality
##  is invisible.  This scans the dose in both phenotypes and reports the
##  forward flow (which is what the gradient steals) rather than the ejection
##  fraction (which it does not touch).
## =============================================================================
bblocker_dose_scan <- function(rates = c(0, 26, 52, 105, 210, 420),
                               end = 336) {
  do.call(rbind, lapply(c("emo", "obstruct"), function(ph) {
    b <- baseline(ph)
    do.call(rbind, lapply(rates, function(r) {
      ch <- if (r == 0) list() else
        list(list(t = 4, par = list(RATE_ESM = r)),
             list(t = 52, par = list(RATE_ESM = 0)))
      o <- run_arm(b, ch, end = end, delta = 1)
      w <- o$time >= 4 & o$time <= 52
      data.frame(phenotype = ph, esmolol_mg_h = r,
                 GRAD_mean_on_drug = round(mean(o$GRAD[w]), 1),
                 CO_mean_on_drug   = round(mean(o$CO[w]), 3),
                 CO_min            = round(min(o$CO), 3),
                 MAP_min           = round(min(o$MAP), 1),
                 LVEF_nadir        = round(min(o$LVEF), 1),
                 LVEF_d14          = round(o$LVEF[which.min(abs(o$time-336))], 1),
                 row.names = NULL)
    }))
  }))
}

## =============================================================================
##  INFERENCE 5 - RECURRENCE OVER ONE YEAR: ACEi vs BETA-BLOCKER
## =============================================================================
recurrence_1y <- function() {
  b <- baseline("emo")
  arms <- list(
    none       = NULL,
    ramipril   = oral_events$ramipril,
    metoprolol = oral_events$metoprolol,
    carvedilol = oral_events$carvedilol
  )
  do.call(rbind, lapply(names(arms), function(a) {
    o <- if (is.null(arms[[a]]))
           mrgsim_df(zero_re(b$mod), end = 8760, delta = 24, recsort = 3)
         else run_oral(b, arms[[a]], end = 8760, delta = 24)
    data.frame(arm = a,
               PREC_1y_pct = round(100*tail(o$PREC, 1), 2),
               LVEF_d365 = round(tail(o$LVEF, 1), 1),
               ACEINH = round(max(o$ACEINH), 3),
               row.names = NULL)
  }))
}

## =============================================================================
##  INFERENCE 6 - THROMBUS AND ANTICOAGULATION
## =============================================================================
thrombus_arms <- function() {
  b <- baseline("emo")
  o1 <- mrgsim_df(zero_re(b$mod), end = 1440, delta = 1, recsort = 3)
  o2 <- run_oral(b, oral_events$apixaban, end = 1440, delta = 1)
  data.frame(arm = c("no anticoagulation", "apixaban 5 mg bd"),
             THR_peak = round(c(max(o1$THR), max(o2$THR)), 4),
             PEMB_60d_pct = round(100*c(tail(o1$PEMB,1), tail(o2$PEMB,1)), 2),
             STASIS_peak = round(c(max(o1$STASIS), max(o2$STASIS)), 3),
             row.names = NULL)
}

## =============================================================================
##  INFERENCE 7 - HOW MUCH OF THE EF DEFICIT IS THE SWITCH?
## =============================================================================
headroom <- function(ids = c("S02","S03","S04","S05","S08","S09","S16")) {
  do.call(rbind, lapply(ids, function(i) {
    o <- run_scenario(i, end = 720)
    data.frame(scenario = i, label = scenarios[[i]]$label,
               share_Gi   = round(tail(o$SHR_GI, 1), 3),
               share_oed  = round(tail(o$SHR_OE, 1), 3),
               share_atp  = round(tail(o$SHR_AT, 1), 3),
               EFdef_hours = round(tail(o$CEFDEF, 1), 2),
               row.names = NULL)
  }))
}

## =============================================================================
##  POPULATION RUN (for the Shiny "population" tab)
## =============================================================================
population <- function(n = 200, seed = 20260730, verbose = TRUE) {
  set.seed(seed)
  b <- baseline("emo")
  idata <- data.frame(
    ID     = seq_len(n),
    RHO_AP = pmax(1.0, rnorm(n, 1.40, 0.16)),
    SEPT   = pmax(0.05, rlnorm(n, log(0.35), 0.75)),
    AMP_TOT= pmax(40, rlnorm(n, log(250), 0.45)),
    FRAC_E = pmin(0.98, pmax(0.05, rnorm(n, 0.68, 0.18))),
    E2     = pmax(5, rlnorm(n, log(18), 0.5)))
  ## One subject at a time, with try().  Running all IDs in a single mrgsim call
  ## is faster but a single non-integrable subject kills the whole cohort, and
  ## the temptation is then to shrink the sampled range until nothing fails -
  ## which silently removes exactly the severe tail the cohort exists to show.
  ## Failures are counted and REPORTED instead (see the F6 limit of validity in
  ## the README: this model has no death, so the sickest draws leave the domain
  ## where its haemodynamics mean anything).
  rows <- list(); nfail <- 0
  for (i in seq_len(n)) {
    p <- as.list(idata[i, setdiff(names(idata), "ID")])
    o <- try(mrgsim_df(zero_re(param(b$mod, p)), end = 2160, delta = 6,
                       recsort = 3), silent = TRUE)
    if (inherits(o, "try-error") || any(!is.finite(o$LVEF))) { nfail <- nfail + 1; next }
    rows[[length(rows) + 1]] <- data.frame(
      ID = idata$ID[i], RHO_AP = p$RHO_AP, SEPT = p$SEPT, AMP_TOT = p$AMP_TOT,
      FRAC_E = p$FRAC_E, E2 = p$E2,
      LVEF_h24 = o$LVEF[which.min(abs(o$time - 24))],
      LVEF_nadir = min(o$LVEF), GRAD_peak = max(o$GRAD),
      ## The registry definition of LVOTO is a SUSTAINED bedside gradient, not
      ## the transient peak during the catecholamine surge, so the comparable
      ## quantity is the mean over the first two days.
      GRAD_mean_48h = mean(o$GRAD[o$time >= 4 & o$time <= 48]),
      BALL_peak = max(o$BALL), TNI_peak = max(o$TNI),
      QTC_peak = max(o$QTC),
      LVEF_d30 = o$LVEF[which.min(abs(o$time - 720))], row.names = NULL)
  }
  out <- do.call(rbind, rows)
  attr(out, "n_failed") <- nfail
  if (verbose && nfail)
    message(sprintf("population(): %d of %d subjects left the domain of validity and were dropped",
                    nfail, n))
  out
}

## =============================================================================
##  RUN EVERYTHING (guarded, so sourcing the file is cheap)
## =============================================================================
if (identical(environment(), globalenv()) &&
    !is.null(getOption("tts.run.scenarios"))) {
  message("--- calibration ---");      print(objective(log(fit_pars), TRUE))
  message("--- all scenarios ---");    print(summarise_all(run_all()))
  message("--- threshold scan ---");   print(threshold_scan())
  message("--- falsification ---");    print(falsification())
  message("--- beta-blocker sign ---");print(bblocker_sign())
  message("--- matched drive ---");    print(matched_drive())
  message("--- headroom ---");         print(headroom())
  message("--- recurrence 1y ---");    print(recurrence_1y())
  message("--- thrombus ---");         print(thrombus_arms())
}
