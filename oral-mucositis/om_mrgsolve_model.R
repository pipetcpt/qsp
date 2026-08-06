# =============================================================================
#  Oral Mucositis (OM) — QSP model for mrgsolve
#  구강점막염 · 항암화학요법 / 방사선치료 유발 — 정량적 시스템 약리학 모델
# =============================================================================
#
#  THE ONE IDEA
#  ------------
#  Oral mucositis is a BALANCE-SHEET FAILURE OF A RENEWING EPITHELIUM.  The
#  mucosa is a conveyor belt
#
#         S -> P1 -> P2 -> P3 -> D -> shed
#
#  with a clonogenic basal pool S, three transit-amplifying generations, and a
#  POST-MITOTIC barrier D.  An ulcer exists exactly while D < D_crit.
#
#  Every cytotoxic insult acts on the PROLIFERATIVE end and spares the
#  POST-MITOTIC end.  That single asymmetry forces TWO CLOCKS:
#
#      ONSET    <- the INSULT term        (hours to days)
#      DURATION <- the REGENERATION term  (days to weeks)
#
#  and the two terms belong to different drugs — cryotherapy / dose /
#  fractionation on one side, palifermin / photobiomodulation / glutamine on
#  the other.  An agent aimed at the wrong clock cannot buy the other
#  endpoint, however large its effect size.
#
#  50 ODEs.  17 scenarios at the bottom of this file.
#
#  PROVENANCE OF THE NUMBERS
#  -------------------------
#  Eight published numbers were spent on eight parameters (see om_calibrate.py
#  and calib.log).  Everything else in this file is structural.  The values
#  below marked  <<FITTED>>  were produced by that calibration; the rest come
#  from the literature cited in om_references.md.  pot_mtx is a PLACEHOLDER:
#  methotrexate is carried for completeness but constrained by nothing here.
#
#  HONEST LIMITATION
#  -----------------
#  No R toolchain was available in the environment that produced this file.
#  The equations below mirror om_python_reference.py — which WAS executed, and
#  whose output is in om_reference_output.txt — line for line, but this file
#  itself has not been run.  Treat a first `mread()` as a syntax check.
#
#  USAGE
#    library(mrgsolve); library(dplyr)
#    mod <- mread("om_mrgsolve_model.R")
#    out <- mod %>% data_set(scen_hdm()) %>% mrgsim(end = 45, delta = 0.05)
#    plot(out, ULC + WHO + VAS + ANC ~ time)
# =============================================================================

$PROB
# Oral mucositis: renewing-epithelium QSP model with two clocks

$PARAM @annotated
// ---- epithelial conveyor belt ---------------------------------------------
lamS   : 0.0604688   : basal clonogen renewal rate (1/d)               <<FITTED>>
S0     : 1.0        : homeostatic clonogen density (normalised)
aS     : 0.3     : flux S -> P1 per unit S (1/d, normalisation-corrected)
k_p    : 0.62       : transit-amplifying maturation rate (1/d)
k_shed : 0.30       : desquamation of the barrier layer (1/d)
fmax   : 6.0        : ceiling on the regeneration multiplier
greg   : 5.84941   : gain of the barrier-deficit -> proliferation loop <<FITTED>>
Kreg   : 0.35       : half-saturation of that loop (deficit fraction)
k_ab   : 0.55       : abortive-division death of doomed basal cells (1/d)
ndiv   : 0.35       : residual flux a doomed cell still contributes
acc_kp : 0.35       : regenerative shortening of transit time
// ---- the QUIESCENT CLONOGEN RESERVE ----------------------------------------
// WITHOUT THIS THE MODEL CANNOT HEAL, AND THE FAILURE IS SILENT.  dS/dt is
// proportional to S, so S = 0 is an ABSORBING state and any insult that drives
// the cycling pool to numerical zero sterilises the tissue permanently.  The
// first calibration run hit exactly that: severe mucositis lasted 39 days and
// never resolved.  The missing biology is the label-retaining, out-of-cycle
// stem cell, and adding it also puts a second time constant (kact) into the
// REGENERATION arm.
Sq0    : 0.16       : size of the quiescent reserve, relative to S0
q_res  : 0.18       : quiescent kill rate, as a fraction of the cycling rate
kact   : 1.35       : activation of reserve -> cycling pool on deficit (1/d)
krest  : 0.10       : restoration of the reserve once S has recovered (1/d)
wP1    : 1.0        : cytotoxic weight on P1
wP2    : 0.6        : cytotoxic weight on P2
wP3    : 0.3        : cytotoxic weight on P3

// ---- ulceration link -------------------------------------------------------
Dcrit  : 0.34       : barrier fraction at which the epithelium ulcerates
wD     : 0.055      : steepness of the ulceration transition

// ---- inflammation ----------------------------------------------------------
kROSd  : 6.0        : ROS clearance (1/d)
kNFb   : 0.5        : NF-kB turnover (1/d)
gROS   : 1.10       : ROS -> NF-kB gain
gTNF   : 0.42       : TNF -> NF-kB gain (POSITIVE FEEDBACK; loop gain < 1)
gMB    : 0.75       : PAMP -> NF-kB gain
KTNF   : 1.0        : TNF half-saturation
KMB    : 1.0        : microbial half-saturation
KROS   : 1.0        : ROS half-saturation
kTNF   : 2.4        : TNF turnover (1/d)
kIL1   : 2.2        : IL-1beta turnover (1/d)
kIL6   : 1.8        : IL-6 turnover (1/d)
pTNF   : 2.4        : TNF production per unit NF-kB (1/d)
pIL1   : 2.2        : IL-1beta production per unit NF-kB (1/d)
pIL6   : 1.8        : IL-6 production per unit NF-kB (1/d)
kCER   : 1.6        : ceramide turnover (1/d)
pCER   : 1.9        : ceramide production (1/d)
aCER   : 0.22       : ceramide -> apoptosis of the POST-MITOTIC layer
aTNFs  : 0.30       : TNF -> extra basal apoptosis

// ---- microbial colonisation of the ulcer bed --------------------------------
kMBg   : 1.5        : colonisation rate (1/d)
kMBd   : 0.9        : neutrophil-dependent clearance (1/d)
MBmax  : 1.0        : carrying capacity
KANC   : 0.5        : ANC half-saturation (x10^9/L)

// ---- nociception ------------------------------------------------------------
kNS    : 0.45       : sensitisation turnover (1/d)
gNS    : 1.5        : ulcer -> sensitisation gain

// ---- regeneration modifiers -------------------------------------------------
kKGFe  : 0.35       : KGF effect-pool turnover (1/d)
eKGF   : 2.6        : max fold-increase of proliferation at saturating KGF
KpalPD : 8.0        : palifermin EC50 on KGFR (ng/mL)
fcycKGF: 0.85       : fraction of the KGF gain that also enlarges the TARGET
                    : (acts ONLY through the cycle-specific agents: setting
                    : it to 0 barely moves an alkylator regimen and turns
                    : palifermin from marginal to near-curative in bolus 5-FU)
kPBMe  : 0.5        : photobiomodulation effect turnover (1/d)
ePBM   : 0.55       : max photobiomodulation effect
KPBM   : 4.0        : PBM half-saturation (J/cm2 accumulated)
kGLNe  : 0.7        : glutamine effect turnover (1/d)
eGLN   : 0.30       : max glutamine effect
eBZD   : 0.45       : max benzydamine suppression of NF-kB production

// ---- systemic PK (V in mL, CL in mL/d unless noted) -------------------------
V_mel  : 37500      : melphalan central volume (mL)
CL_mel : 520000     : melphalan clearance (mL/d)
Q_mel  : 150000     : melphalan intercompartmental clearance (mL/d)
Vp_mel : 30000      : melphalan peripheral volume (mL)
V_5fu  : 25000      : 5-FU volume (mL)
Vmx5fu : 2.2e6      : 5-FU saturable elimination Vmax (ug/d)
Km5fu  : 2.0        : 5-FU Km (ug/mL)
CLl5fu : 2.0e6      : 5-FU linear clearance (mL/d)
V_mtx  : 20         : methotrexate central volume (L)
CL_mtx : 140        : methotrexate clearance (L/d)
Q1_mtx : 90         : methotrexate Q1 (L/d)
Vp1mtx : 12         : methotrexate V2 (L)
Q2_mtx : 8          : methotrexate Q2 (L/d)
Vp2mtx : 30         : methotrexate V3 (L)
V_cis  : 45000      : cisplatin (free Pt) volume (mL)
CL_cis : 1.83e6     : cisplatin clearance (mL/d)
Q_cis  : 2.0e5      : cisplatin Q (mL/d)
Vp_cis : 60000      : cisplatin peripheral volume (mL)
V_pal  : 120000     : palifermin central volume (mL)
CL_pal : 1.386e6    : palifermin clearance (mL/d)
Q_pal  : 4.0e5      : palifermin Q (mL/d)
Vp_pal : 255000     : palifermin peripheral volume (mL)
V_mor  : 1.0e5      : morphine central volume (mL)
CL_mor : 2.16e6     : morphine clearance (mL/d)
Q_mor  : 3.0e5      : morphine Q (mL/d)
Vp_mor : 1.5e5      : morphine peripheral volume (mL)
ke0mor : 6.0        : morphine effect-site rate (1/d)
EC50mor: 25.0       : morphine effect-site EC50 (ng/mL)
Emaxmor: 6.5        : max morphine analgesia (VAS points)
k_titr : 1.35e7     : opioid titration gain (ng/d per VAS point)
opitgt : 4.0        : VAS target for opioid titration

// ---- mucosal delivery compartment --------------------------------------------
keqmuc : 500.0      : mucosal equilibration rate constant (1/d)
Kp_mel : 0.75       : melphalan tissue:plasma partition
Kp_5fu : 0.90       : 5-FU tissue:plasma partition
Kp_mtx : 0.55       : methotrexate tissue:plasma partition
Kp_cis : 0.45       : cisplatin tissue:plasma partition
fQcryo : 0.22       : mucosal blood flow during ice (fraction of basal)
Q10cry : 2.5        : Q10 for alkylation / cell-cycle rate
dTcryo : 13.0       : mucosal temperature drop during ice (deg C)

// ---- potencies (per unit mucosal concentration per day) ----------------------
pot_mel: 23.2215 : melphalan potency on the clonogen pool           <<FITTED>>
pot_5fu: 1.83921 : 5-FU potency on the clonogen pool                <<FITTED>>
pot_mtx: 0.085      : methotrexate potency
pot_cis: 15.625 : cisplatin potency                                <<FITTED>>
cs_5fu : 1          : 5-FU is S-phase specific (kill scales with cycling)
cs_mtx : 1          : methotrexate is S-phase specific
cs_mel : 0          : melphalan is not
cs_cis : 0          : cisplatin is not
sens   : 1.0        : individual sensitivity multiplier (IIV handle)

// ---- radiation ---------------------------------------------------------------
alpha_m: 0.30       : mucosa alpha (1/Gy)
beta_m : 0.030      : mucosa beta (1/Gy^2)   -> alpha/beta = 10 Gy
alpha_t: 0.30       : tumour alpha (1/Gy)
beta_t : 0.030      : tumour beta (1/Gy^2)
rad_pot: 0.505097 : conversion of mucosal log-kill to clonogen loss   <<FITTED>>
Tk_tum : 28.0       : onset of tumour accelerated repopulation (d)
lamtum : 0.198      : tumour repopulation rate (1/d)  = ln2 / 3.5 d

// ---- myelosuppression (Friberg) ----------------------------------------------
MTT    : 5.2083     : mean transit time (d) = 125 h
gamanc : 0.161      : feedback exponent
Circ0  : 5.0        : baseline ANC (x10^9/L)
ktrscl : 4.0        : number of transit compartments
slpmel : 0.055      : melphalan marrow slope (per ug/mL)
slp5fu : 0.0016     : 5-FU marrow slope
slpcis : 0.030      : cisplatin marrow slope
slpmtx : 0.020      : methotrexate marrow slope

// ---- clinical linkage ---------------------------------------------------------
painmx : 10.0       : maximum VAS
paink  : 8.5        : ulcer-area -> pain gain
infbas : 0.004      : baseline bacteraemia hazard (1/d)
infgn  : 0.55       : ulcer/microbe -> bacteraemia hazard gain

// ---- TIME-VARYING COVARIATES, supplied as data-set columns --------------------
// These are how the schedule enters the model.  RRATE is the instantaneous
// radiotherapy dose rate (Gy/d): a 2 Gy fraction over 10 min is RRATE = 288
// for 0.00694 d.  See the scenario builders at the bottom.
RRATE  : 0          : radiotherapy dose rate (Gy/d)
DPF    : 2.0        : dose per fraction (Gy) — enters the quadratic LQ term
CRYO   : 0          : 1 while ice is in the mouth
PBMR   : 0          : photobiomodulation delivery rate (J/cm2/d)
BZDL   : 0          : benzydamine exposure level (0-1)
GLNL   : 0          : glutamine exposure level (0-1)
OPION  : 1          : 1 = run the opioid titration loop

$CMT @annotated
// ---- epithelium (8) ---------------------------------------------------------
S       : cycling clonogenic basal cells (normalised)
Sq      : quiescent clonogen reserve (makes S = 0 non-absorbing)
Sd      : lethally hit basal cells, still present
P1      : transit-amplifying generation 1
P2      : transit-amplifying generation 2
P3      : transit-amplifying generation 3
D       : post-mitotic barrier (normalised to 1 at homeostasis)
Dh      : reserved (barrier history, unused)
// ---- inflammation (8) --------------------------------------------------------
ROS     : reactive oxygen species
NFkB    : NF-kB activity
TNF     : TNF-alpha
IL1b    : IL-1beta
IL6     : IL-6
CER     : ceramide
MB      : microbial / PAMP load on the ulcer bed
NS      : nociceptor sensitisation state
// ---- regeneration modifiers (3) -----------------------------------------------
KGFe    : KGF pharmacodynamic effect pool
PBMe    : photobiomodulation effect pool
GLNe    : glutamine effect pool
// ---- PK (12) --------------------------------------------------------------------
A_mel_c : melphalan central (ug)
A_mel_p : melphalan peripheral (ug)
A_5fu_c : 5-FU central (ug)
A_5fu_g : reserved
A_mtx_c : methotrexate central (umol)
A_mtx_1 : methotrexate peripheral 1 (umol)
A_mtx_2 : methotrexate peripheral 2 (umol)
A_cis_c : cisplatin central (ug)
A_cis_p : cisplatin peripheral (ug)
A_pal_c : palifermin central (ng)
A_pal_p : palifermin peripheral (ng)
A_mor_c : morphine central (ng)
A_mor_p : morphine peripheral (ng)
A_mor_e : morphine effect site (ng/mL)
// ---- mucosal tissue concentrations (4) ------------------------------------------
Cm_mel  : mucosal melphalan (ug/mL)
Cm_5fu  : mucosal 5-FU (ug/mL)
Cm_mtx  : mucosal methotrexate (umol/L)
Cm_cis  : mucosal cisplatin (ug/mL)
// ---- myelosuppression (5) ---------------------------------------------------------
Prol    : proliferating marrow pool
Tr1     : marrow transit 1
Tr2     : marrow transit 2
Tr3     : marrow transit 3
Circ    : circulating ANC (x10^9/L)
// ---- accumulators (8) -------------------------------------------------------------
cAUCmuc : cumulative mucosal drug exposure
cUlcD   : cumulative ulcer days
cOpiD   : cumulative opioid days
cBEDt   : tumour BED with repopulation (Gy10)
cBEDm   : mucosal BED (Gy10)
cInf    : cumulative bacteraemia hazard
cSev    : cumulative severe-mucositis days
cPain   : cumulative pain (VAS.d)

$GLOBAL
#define HILL(x, K) (((x) > 0 ? (x) : 0) / ((K) + ((x) > 0 ? (x) : 0)))
#define POS(x)     ((x) > 0 ? (x) : 0)

$MAIN
// Homeostatic initial condition.  The belt is at its drug-free attractor and
// the barrier is normalised to 1, which is what makes Dcrit a fraction of a
// real steady state rather than of a guess.  aS is supplied already corrected
// for that normalisation (see mkmrgsolve.py).
double Pss = aS * S0 / k_p;
S_0    = S0;
Sq_0   = Sq0;
P1_0   = Pss;
P2_0   = Pss;
P3_0   = Pss;
D_0    = k_p * Pss / k_shed;
Prol_0 = Circ0;
Tr1_0  = Circ0;
Tr2_0  = Circ0;
Tr3_0  = Circ0;
Circ_0 = Circ0;

$ODE
// ===========================================================================
// PHARMACOKINETICS
// ===========================================================================
double Cmel  = A_mel_c / V_mel;
double Cmelp = A_mel_p / Vp_mel;
dxdt_A_mel_c = -CL_mel * Cmel - Q_mel * (Cmel - Cmelp);
dxdt_A_mel_p =  Q_mel * (Cmel - Cmelp);

double C5fu = A_5fu_c / V_5fu;
dxdt_A_5fu_c = -(Vmx5fu * C5fu / (Km5fu + C5fu) + CLl5fu * C5fu);
dxdt_A_5fu_g = 0.0;

double Cmtx = A_mtx_c / V_mtx;
double Cmt1 = A_mtx_1 / Vp1mtx;
double Cmt2 = A_mtx_2 / Vp2mtx;
dxdt_A_mtx_c = -CL_mtx * Cmtx - Q1_mtx * (Cmtx - Cmt1) - Q2_mtx * (Cmtx - Cmt2);
dxdt_A_mtx_1 =  Q1_mtx * (Cmtx - Cmt1);
dxdt_A_mtx_2 =  Q2_mtx * (Cmtx - Cmt2);

double Ccis  = A_cis_c / V_cis;
double Ccisp = A_cis_p / Vp_cis;
dxdt_A_cis_c = -CL_cis * Ccis - Q_cis * (Ccis - Ccisp);
dxdt_A_cis_p =  Q_cis * (Ccis - Ccisp);

double Cpal  = A_pal_c / V_pal;
double Cpalp = A_pal_p / Vp_pal;
dxdt_A_pal_c = -CL_pal * Cpal - Q_pal * (Cpal - Cpalp);
dxdt_A_pal_p =  Q_pal * (Cpal - Cpalp);

// ===========================================================================
// MUCOSAL DELIVERY -- the ONLY place cryotherapy enters.
// It touches no potency, no repair and no cytokine: every cryotherapy result
// in om_reference_output.txt is a consequence of flow and temperature alone.
// ===========================================================================
double iced  = (CRYO > 0.5) ? 1.0 : 0.0;
double fQ    = iced ? fQcryo : 1.0;
double Qd    = keqmuc * fQ;
double ftemp = iced ? (1.0 / pow(Q10cry, dTcryo / 10.0)) : 1.0;

dxdt_Cm_mel = Qd * (Cmel - Cm_mel / Kp_mel);
dxdt_Cm_5fu = Qd * (C5fu - Cm_5fu / Kp_5fu);
dxdt_Cm_mtx = Qd * (Cmtx - Cm_mtx / Kp_mtx);
dxdt_Cm_cis = Qd * (Ccis - Cm_cis / Kp_cis);

// ===========================================================================
// REGENERATION SIGNAL
// ===========================================================================
double deficit = POS(1.0 - D);
double fdef    = 1.0 + greg * deficit / (Kreg + deficit);

double kgfdrv  = HILL(Cpal, KpalPD);
dxdt_KGFe      = kKGFe * (kgfdrv - KGFe);
double fkgf    = 1.0 + eKGF * KGFe;

dxdt_PBMe      = PBMR - kPBMe * PBMe;
double fpbm    = 1.0 + ePBM * HILL(PBMe, KPBM);

dxdt_GLNe      = GLNL - kGLNe * GLNe;
double fgln    = 1.0 + eGLN * HILL(GLNe, 1.0);

double freg    = fdef * fkgf * fpbm * fgln;
if (freg > fmax) freg = fmax;

// THE CYCLING FRACTION SEEN BY A CYCLE-ACTIVE CYTOTOXIC.
// KGF raises it.  That is the whole palifermin scheduling paradox, and
// fcycKGF is the single parameter that carries it.
double fcyc = (1.0 + fcycKGF * (fkgf - 1.0)) * sqrt(fdef);

// ===========================================================================
// CYTOTOXIC INSULT ON THE PROLIFERATIVE POOL
// ===========================================================================
double km_mel = ftemp * pot_mel * Cm_mel * (cs_mel > 0.5 ? fcyc : 1.0);
double km_5fu = ftemp * pot_5fu * Cm_5fu * (cs_5fu > 0.5 ? fcyc : 1.0);
double km_mtx = ftemp * pot_mtx * Cm_mtx * (cs_mtx > 0.5 ? fcyc : 1.0);
double km_cis = ftemp * pot_cis * Cm_cis * (cs_cis > 0.5 ? fcyc : 1.0);
double k_chem = sens * (km_mel + km_5fu + km_mtx + km_cis);

// linear-quadratic radiation kill, expressed as an instantaneous rate
double k_rad  = sens * ftemp * rad_pot * RRATE * (alpha_m + beta_m * DPF);

double k_tot  = k_chem + k_rad + aTNFs * HILL(TNF, 1.5);

// ===========================================================================
// THE EPITHELIAL CONVEYOR BELT -- the spine
// ===========================================================================
double sdef = POS(1.0 - S / S0);
double act  = kact * Sq * sdef;
dxdt_S  = lamS * freg * S * (1.0 - S / S0) + act - k_tot * S;
dxdt_Sq = -q_res * k_tot * Sq - act + krest * POS(Sq0 - Sq) * (S / S0);
dxdt_Sd = k_tot * S + q_res * k_tot * Sq - k_ab * Sd;

// Regenerating epithelium does not merely divide faster, it also SHORTENS the
// transit.  Without this the healing limb is floored at 3/k_p + 1/k_shed and
// no drug could shorten it -- which would make the "duration is regeneration"
// claim untestable by construction.
double kpe     = k_p * (1.0 + acc_kp * (freg - 1.0));
double flux_in = aS * freg * (S + ndiv * Sd);

dxdt_P1 = flux_in - kpe * P1 - wP1 * k_tot * P1;
dxdt_P2 = kpe * P1 - kpe * P2 - wP2 * k_tot * P2;
dxdt_P3 = kpe * P2 - kpe * P3 - wP3 * k_tot * P3;

// D is POST-MITOTIC: there is NO cytotoxic term here.  That absence is the
// reason the whole disease has a latent period.
dxdt_D  = kpe * P3 - k_shed * D - aCER * HILL(CER, 1.4) * D;
dxdt_Dh = 0.0;

double A_ulc = 1.0 / (1.0 + exp((D - Dcrit) / wD));

// ===========================================================================
// INFLAMMATION
// ===========================================================================
double ros_in = 2.2 * RRATE
              + 1.6 * (Cm_mel * pot_mel + Cm_5fu * pot_5fu
                       + Cm_mtx * pot_mtx + Cm_cis * pot_cis);
dxdt_ROS = ros_in - kROSd * ROS;

double bzd   = eBZD * HILL(BZDL, 0.5);
double drive = gROS * HILL(ROS, KROS)
             + gTNF * HILL(TNF, KTNF)
             + gMB  * HILL(MB,  KMB);
dxdt_NFkB = kNFb * ((1.0 - bzd) * drive - NFkB);

dxdt_TNF  = pTNF * NFkB - kTNF * TNF;
dxdt_IL1b = pIL1 * NFkB - kIL1 * IL1b;
dxdt_IL6  = pIL6 * NFkB - kIL6 * IL6;
dxdt_CER  = pCER * (0.6 * NFkB + 0.4 * HILL(ROS, 0.8)) - kCER * CER;

double anc = Circ > 1e-3 ? Circ : 1e-3;
dxdt_MB = kMBg * A_ulc * (1.0 - MB / MBmax) - kMBd * MB * anc / (KANC + anc);

// ===========================================================================
// NOCICEPTION AND THE OPIOID TITRATION LOOP
// ===========================================================================
double opi = Emaxmor * HILL(A_mor_e, EC50mor);
dxdt_NS    = kNS * (gNS * A_ulc * (1.0 + 0.5 * HILL(IL1b, 1.0)) - NS);

double xp      = paink * A_ulc * (0.4 + 0.6 * NS);
double vas_raw = painmx * xp / (1.0 + xp);
double vas     = vas_raw - opi;
if (vas < 0.0) vas = 0.0;

double Cmor  = A_mor_c / V_mor;
double Cmorp = A_mor_p / Vp_mor;
double need  = (OPION > 0.5 && vas_raw > opitgt)
             ? k_titr * (vas_raw - opitgt) : 0.0;
dxdt_A_mor_c = need - CL_mor * Cmor - Q_mor * (Cmor - Cmorp);
dxdt_A_mor_p = Q_mor * (Cmor - Cmorp);
dxdt_A_mor_e = ke0mor * (Cmor - A_mor_e);

// ===========================================================================
// MYELOSUPPRESSION (Friberg)
// ===========================================================================
double ktr   = ktrscl / MTT;
double edrug = slpmel * Cmel + slp5fu * C5fu + slpcis * Ccis + slpmtx * Cmtx;
if (edrug > 0.98) edrug = 0.98;
double fb = pow(Circ0 / anc, gamanc);
if (fb > 6.0) fb = 6.0;
dxdt_Prol = ktr * Prol * (1.0 - edrug) * fb - ktr * Prol;
dxdt_Tr1  = ktr * (Prol - Tr1);
dxdt_Tr2  = ktr * (Tr1 - Tr2);
dxdt_Tr3  = ktr * (Tr2 - Tr3);
dxdt_Circ = ktr * (Tr3 - Circ);

// ===========================================================================
// ACCUMULATORS
// ===========================================================================
double rt_on = (RRATE > 0.0) ? 1.0 : 0.0;
double rep   = (SOLVERTIME > Tk_tum) ? lamtum * rt_on : 0.0;

dxdt_cAUCmuc = Cm_mel + Cm_5fu + Cm_mtx + Cm_cis;
dxdt_cUlcD   = (A_ulc > 0.5) ? 1.0 : 0.0;
dxdt_cOpiD   = (need > 0.0) ? 1.0 : 0.0;
dxdt_cBEDt   = RRATE * (1.0 + DPF / (alpha_t / beta_t)) - rep / alpha_t;
dxdt_cBEDm   = RRATE * (1.0 + DPF / (alpha_m / beta_m));
dxdt_cInf    = infbas + infgn * A_ulc * MB / (KANC + anc);
dxdt_cSev    = (A_ulc > 0.5) ? 1.0 : 0.0;
dxdt_cPain   = vas;

$TABLE
double ULC = 1.0 / (1.0 + exp((D - Dcrit) / wD));
double ERY = 1.0 / (1.0 + exp((D - 0.72) / 0.06));

// WHO oral toxicity grade.  ORDINAL, and grade 4 is ABSORBING -- which is
// exactly why om_analysis.py section 7 finds an assay-sensitivity loss where
// the disease is worst.
double WHO = 0.0;
if (ERY > 0.40) WHO = 1.0;
if (ULC > 0.08) WHO = 2.0;
if (ULC > 0.35) WHO = 3.0;
if (ULC > 0.72) WHO = 4.0;

// OMAS: continuous, area-based, does NOT saturate
double OMAS = 5.0 * pow(ULC, 0.8);

double OPIEFF  = Emaxmor * HILL(A_mor_e, EC50mor);
double XP      = paink * ULC * (0.4 + 0.6 * NS);
double VAS_RAW = painmx * XP / (1.0 + XP);
double VAS     = VAS_RAW - OPIEFF;
if (VAS < 0.0) VAS = 0.0;

double CPAL  = A_pal_c / V_pal;
double CMEL  = A_mel_c / V_mel;
double C5FU  = A_5fu_c / V_5fu;
double CCIS  = A_cis_c / V_cis;
double CMOR  = A_mor_c / V_mor;
double MEDMG = (OPION > 0.5 && VAS_RAW > opitgt)
             ? k_titr * (VAS_RAW - opitgt) / 1.0e6 : 0.0;   // mg/d IV morphine
double ANC   = Circ;
double BARRIER = D;

$CAPTURE @annotated
ULC     : ulcerated fraction of the oral mucosa (0-1)
WHO     : WHO oral toxicity grade (0-4)
OMAS    : OMAS ulceration score (0-5)
VAS     : pain VAS after analgesia (0-10)
VAS_RAW : pain VAS before analgesia (0-10)
MEDMG   : IV morphine-equivalent requirement (mg/d)
ANC     : absolute neutrophil count (x10^9/L)
BARRIER : post-mitotic barrier mass (1 = intact)
CMEL    : melphalan plasma concentration (ug/mL)
C5FU    : 5-FU plasma concentration (ug/mL)
CCIS    : cisplatin free plasma concentration (ug/mL)
CPAL    : palifermin plasma concentration (ng/mL)
CMOR    : morphine plasma concentration (ng/mL)

# =============================================================================
#  SCENARIOS
# =============================================================================
#  Radiotherapy and the topical/physical interventions enter as TIME-VARYING
#  COVARIATES in the data set, not as dosing records: RRATE is the
#  instantaneous dose rate in Gy/d (a 2 Gy fraction over 10 min is
#  RRATE = 288 for 0.00694 d), CRYO/BZDL/GLNL are 0/1 levels, PBMR is a
#  delivery rate in J/cm2/d.
#
#  IMPORTANT: give the solver the fraction boundaries as records.  A 10-minute
#  window inside a 60-day simulation will otherwise be stepped over, and the
#  delivered dose becomes a property of the integrator rather than of the
#  prescription.  This exact defect cost the Python reference a factor of
#  three in delivered dose before it was caught (om_python_reference.py,
#  `breakpoints`).
#
#  $ Rscript-style helpers -------------------------------------------------
#  library(dplyr); library(mrgsolve)
#
#  # ---- building blocks ---------------------------------------------------
#  cov_frame <- function(times, ...) {
#      d <- tibble::tibble(ID = 1, time = times, evid = 0, cmt = 0, amt = 0,
#                          rate = 0, RRATE = 0, DPF = 2, CRYO = 0, PBMR = 0,
#                          BZDL = 0, GLNL = 0)
#      mods <- list(...)
#      for (nm in names(mods)) d[[nm]] <- mods[[nm]]
#      d
#  }
#
#  # Fractions delivered over 10 minutes each.  per_day > 1 gives genuine
#  # b.i.d./t.i.d. delivery with a bid_gap_h interfraction interval -- do NOT
#  # express hyperfractionation as per_week = 10 with one fraction per day.
#  # That stretches 68 fractions over 68 DAYS instead of 34, and overall
#  # treatment time is the entire point of an altered-fractionation
#  # comparison.  The Python reference had exactly that bug, and it showed up
#  # as a LOWER tumour BED for the hyperfractionated arm.
#  rt_records <- function(nfx = 35, dpf = 2, per_week = 5, per_day = 1,
#                         bid_gap_h = 6, t0 = 0, dur = 10/1440) {
#      day <- t0; given <- 0; starts <- c()
#      while (given < nfx && (day - t0) < 300) {
#          if ((round(day - t0) %% 7) < per_week) {
#              for (j in seq_len(per_day)) {
#                  if (given >= nfx) break
#                  starts <- c(starts, day + (j - 1) * bid_gap_h / 24)
#                  given <- given + 1
#              }
#          }
#          day <- day + 1
#      }
#      bind_rows(lapply(starts, function(s)
#          tibble::tibble(time = c(s, s + dur), RRATE = c(dpf/dur, 0))))
#  }
#
#  # ---- 1. HDM 200 mg/m2, no prophylaxis ---------------------------------
#  scen_hdm <- function(dose_mgm2 = 200, bsa = 1.8) {
#      ev(ID = 1, time = 0, amt = dose_mgm2 * bsa * 1000, cmt = "A_mel_c",
#         rate = dose_mgm2 * bsa * 1000 / (0.5/24)) %>% as_data_frame()
#  }
#
#  # ---- 2. HDM 140 mg/m2 (reduced dose, INSULT arm) ----------------------
#  scen_hdm140 <- function() scen_hdm(140)
#
#  # ---- 3/4. HDM + oral cryotherapy (INSULT arm) -------------------------
#  #   30 min and 6 h windows; CRYO switches 0 -> 1 -> 0
#  scen_hdm_cryo <- function(hours = 6) {
#      dosing <- scen_hdm()
#      cov <- tibble::tibble(ID = 1, time = c(0, hours/24), CRYO = c(1, 0))
#      bind_rows(dosing, cov) %>% arrange(time)
#  }
#
#  # ---- 5. HDM + palifermin, SEPARATED (label schedule) ------------------
#  #   60 ug/kg/d IV on d -3,-2,-1 and d +3,+4,+5
#  scen_hdm_pal <- function(wt = 75) {
#      amt <- 60 * wt * 1000                       # ng
#      bind_rows(scen_hdm(),
#                ev(ID = 1, time = c(0, 1, 2, 6, 7, 8) - 3 + 3,
#                   amt = amt, cmt = "A_pal_c") %>% as_data_frame())
#  }
#
#  # ---- 6. HDM + palifermin, CONCURRENT (the paradox arm) ----------------
#  #   identical total dose, moved so the KGF effect pool is HIGH when the
#  #   alkylator lands.  The model should show the benefit shrink or invert.
#  scen_hdm_pal_conc <- function(wt = 75) {
#      amt <- 60 * wt * 1000
#      bind_rows(scen_hdm(),
#                ev(ID = 1, time = c(-1, -0.5, 0, 3, 4, 5) + 1,
#                   amt = amt, cmt = "A_pal_c") %>% as_data_frame())
#  }
#
#  # ---- 7. HDM + photobiomodulation (REGENERATION arm) -------------------
#  scen_hdm_pbm <- function(days = 0:20, fluence = 6) {
#      cov <- bind_rows(lapply(days, function(d)
#          tibble::tibble(ID = 1, time = c(d, d + 5/1440),
#                         PBMR = c(fluence / (5/1440), 0))))
#      bind_rows(scen_hdm(), cov) %>% arrange(time)
#  }
#
#  # ---- 8. HDM + oral glutamine ------------------------------------------
#  scen_hdm_gln <- function() bind_rows(
#      scen_hdm(), tibble::tibble(ID = 1, time = c(0, 21), GLNL = c(1, 0)))
#
#  # ---- 9. HDM + cryotherapy + palifermin (both clocks) ------------------
#  scen_hdm_both <- function() bind_rows(scen_hdm_cryo(6), scen_hdm_pal()) %>%
#      arrange(time)
#
#  # ---- 10/11. TBI-VP16-Cy conditioning, +/- palifermin ------------------
#  #   fractionated TBI 12 Gy in 8 fx of 1.5 Gy over 4 days, then the
#  #   alkylator block as a melphalan equivalent (see om_calibrate.py stage 2)
#  scen_tbi <- function(cy_equiv = 50.1437, bsa = 1.8, pal = FALSE,
#                       wt = 75) {
#      dur <- 10/1440
#      rt  <- bind_rows(lapply(c(0, 0.35, 1, 1.35, 2, 2.35, 3, 3.35),
#          function(s) tibble::tibble(ID = 1, time = c(s, s + dur),
#                                     RRATE = c(1.5/dur, 0), DPF = 1.5)))
#      amt <- cy_equiv * bsa * 1000
#      alk <- ev(ID = 1, time = c(4, 6), amt = amt/2, cmt = "A_mel_c",
#                rate = (amt/2) / (0.5/24)) %>% as_data_frame()
#      out <- bind_rows(rt, alk)
#      if (pal) out <- bind_rows(out,
#          ev(ID = 1, time = c(0, 1, 2, 11, 12, 13),
#             amt = 60 * wt * 1000, cmt = "A_pal_c") %>% as_data_frame())
#      arrange(out, time)
#  }
#
#  # ---- 12/13. Head-and-neck 70 Gy/35 fx, +/- cisplatin ------------------
#  scen_chemort <- function(cis_mgm2 = 100, bsa = 1.8, nfx = 35, dpf = 2,
#                           per_week = 5, per_day = 1) {
#      rt <- rt_records(nfx, dpf, per_week, per_day) %>%
#            mutate(ID = 1, DPF = dpf)
#      out <- rt
#      if (cis_mgm2 > 0) {
#          amt <- cis_mgm2 * bsa * 1000
#          out <- bind_rows(out,
#              ev(ID = 1, time = c(0, 21, 42), amt = amt, cmt = "A_cis_c",
#                 rate = amt / (2/24)) %>% as_data_frame())
#      }
#      arrange(out, time)
#  }
#
#  # ---- 14. Hyperfractionation 81.6 Gy / 68 fx b.i.d. --------------------
#  scen_hyperfx <- function() scen_chemort(cis_mgm2 = 0, nfx = 68, dpf = 1.2,
#                                          per_week = 5, per_day = 2)
#
#  # ---- 15. Head-and-neck chemoRT + benzydamine --------------------------
#  scen_chemort_bzd <- function() bind_rows(
#      scen_chemort(), tibble::tibble(ID = 1, time = c(0, 60), BZDL = c(1, 0))
#  ) %>% arrange(time)
#
#  # ---- 16. 5-FU bolus 425 mg/m2 d1-5 ------------------------------------
#  scen_5fu_bolus <- function(dose = 425, bsa = 1.8) {
#      amt <- dose * bsa * 1000
#      ev(ID = 1, time = 0:4, amt = amt, cmt = "A_5fu_c",
#         rate = amt / (5/1440)) %>% as_data_frame()
#  }
#
#  # ---- 17. 5-FU 96-h continuous infusion 4000 mg/m2 ---------------------
#  #   Same drug, same total dose, ~200x longer exposure.  This is the arm
#  #   that makes the cryotherapy criterion visible: 30 min of ice covers
#  #   1 - 2^(-T/t_half) of a 12-minute half-life and essentially none of a
#  #   96-hour infusion.
#  scen_5fu_ci <- function(dose = 4000, bsa = 1.8, dur_d = 4) {
#      amt <- dose * bsa * 1000
#      ev(ID = 1, time = 0, amt = amt, cmt = "A_5fu_c",
#         rate = amt / dur_d) %>% as_data_frame()
#  }
#
#  # ---- running one ------------------------------------------------------
#  #  mod <- mread("om_mrgsolve_model.R")
#  #  out <- mod %>% data_set(scen_hdm_cryo(6)) %>%
#  #         mrgsim(end = 45, delta = 0.05, recsort = 3)
#  #  plot(out, ULC + WHO + VAS + ANC ~ time)
#
#  # ---- virtual population ------------------------------------------------
#  #  Inter-individual variability goes on exactly three handles, because the
#  #  structure says those are the three places a patient can differ:
#  #  overall cytotoxic sensitivity, regenerative capacity, barrier threshold.
#  #
#  #  idata <- tibble::tibble(
#  #      ID    = 1:400,
#  #      sens  = exp(rnorm(400, 0, 0.38) - 0.5 * 0.38^2),
#  #      lamS  = 0.0604688 * exp(rnorm(400, 0, 0.30) - 0.5 * 0.30^2),
#  #      Dcrit = 0.34    * exp(rnorm(400, 0, 0.07) - 0.5 * 0.07^2))
#  #  mod %>% idata_set(idata) %>% data_set(scen_hdm()) %>% mrgsim(end = 45)
# =============================================================================
