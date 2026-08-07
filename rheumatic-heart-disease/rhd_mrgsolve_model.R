## =====================================================================
##  Rheumatic heart disease (RHD)
##  Acute rheumatic fever -> rheumatic mitral stenosis
##  A quantitative systems pharmacology model for mrgsolve
## ---------------------------------------------------------------------
##  WHAT THIS MODEL IS FOR
##
##  RHD is managed on two numbers, and neither of them is the number
##  that matters.
##
##  1. Secondary prophylaxis is audited as ADHERENCE (doses received /
##     doses due).  What protects is TIME ABOVE 0.02 ug/mL of penicillin
##     G.  Benzathine penicillin G is a depot with flip-flop kinetics:
##     plasma concentration tracks RELEASE (t1/2 ~9 d), not clearance
##     (t1/2 30 min), and its amplitude is inversely proportional to
##     body size.  One 1.2 MU injection covers ~18 d in a 70 kg adult
##     and ~23 d in a 40 kg child.  Adherence and protection therefore
##     rank patients differently, and CAN RANK THEM IN OPPOSITE ORDERS:
##     an 80%-adherent 55 kg patient on a 3-weekly interval spends 79%
##     of the year protected, a 100%-adherent 95 kg patient on a
##     4-weekly interval spends 60%.
##
##  2. Valve area is destroyed by TWO arms with different owners:
##
##       immune arm      -2 sqrt(pi A) * KFI * VIT
##       autonomous arm  -2 sqrt(pi A) * KFS * (MVG/4) * BRAKE
##
##     The autonomous arm is shear.  By the simplified Bernoulli
##     relation the squared mean transvalvular velocity IS the mean
##     gradient over four, so the shear driver is the gradient itself,
##     and the gradient rises as (flow / area)^2.  The arm is therefore
##     self-accelerating, and BELOW a crossover area (2.75 cm2 here --
##     computed, not assumed) it exceeds the expected immune loss of a
##     completely unprotected patient.  Prophylaxis buys most of its
##     VALVE benefit above that line, which is the latent-RHD window
##     that echocardiographic screening finds and that GOAL randomised.
##     It buys its RECURRENCE benefit everywhere; these are different
##     endpoints and the model separates them.
##
##  3. The Gorlin block is why a patient decompensates without the
##     valve changing:  MVG = (CO / (37.7 * A * DFPmin))^2.  Pregnancy
##     raises CO; atrial fibrillation shortens diastole and removes the
##     atrial contribution.  At a FIXED 1.5 cm2 orifice, pregnancy
##     takes the mean gradient from 6.4 to 23.0 mmHg and left atrial
##     pressure from 12.4 to 29.0 mmHg.  Read backwards, the same
##     equation yields an OPTIMAL heart rate that falls as the valve
##     narrows -- 108 bpm at 2.0 cm2, 66 at 1.0, 55 at 0.8 -- because
##     total diastolic time per minute rises as rate falls while stroke
##     volume is capped by the ventricle.
##
## ---------------------------------------------------------------------
##  STATUS.  All 39 ODEs were independently re-implemented in
##  Python/scipy (rhd_verify_python.py) and run against 47 published or
##  derived anchors; 47/47 pass.  That exercise found and fixed seven
##  real defects, which are listed in README.md rather than quietly
##  corrected.  Parameters are illustrative and calibrated to published
##  central estimates; this is a teaching and hypothesis-generating
##  model, not a clinical tool.
##
##  Units: time = DAYS, concentrations = ug/mL (= mg/L), amounts = mg
##  unless noted, areas = cm2, pressures = mmHg, flow = L/min.
## =====================================================================

library(mrgsolve)

code <- '
$PARAM @annotated
// ---- demographics -------------------------------------------------
WT     :  70    : Body weight (kg)
CRCL   : 100    : Creatinine clearance (mL/min)
AGE    :  25    : Age (years)

// ---- penicillin pharmacokinetics -----------------------------------
VPEN0  :  15    : Apparent V/F of penicillin G per 70 kg (L)
CLPEN0 : 500    : Penicillin G clearance per 70 kg at CrCl 100 (L/day)
KAFBPG :   0.35 : Fast-release rate constant of the BPG depot (1/day)
KASBPG :   0.075: Slow-release rate constant of the BPG depot (1/day)
FFBPG  :   0.30 : Fraction of a BPG dose entering the fast depot (-)
FSBPG  :   0.70 : Fraction of a BPG dose entering the slow depot (-)
KAORAL :  24    : Oral penicillin V absorption rate constant (1/day)
MICP   :   0.02 : Protective plasma threshold (ug/mL)
KEO    :   5.545: Plasma to tonsillar effect-site rate constant (1/day)

// ---- group A streptococcal pharyngitis ------------------------------
KGROW  :   2.0  : GAS maximum growth rate (1/day)
GMAX   :   1e8  : GAS carrying capacity (CFU-equivalents)
KGCL   :   0.35 : Innate clearance of GAS (1/day)
KMIK   :   1.2  : Clearance per unit mucosal immunity (1/day)
KPKILL :   6.0  : Maximum penicillin kill rate (1/day)
GEXT   : 100    : Extinction floor below which GAS cannot regrow (CFU-eq)
HPEN   :   3.0  : Hill coefficient of the penicillin kill effect (-)
EC50P  :   0.02 : Effect-site EC50 for GAS killing (ug/mL)
GM50   :   1e6  : Burden giving half-maximal systemic antigen release
GM50I  :   1e4  : Burden giving half-maximal immune stimulus
KMI    :   0.30 : Mucosal immunity formation rate (1/day)
KMD    :   0.015: Mucosal immunity decay rate (1/day)

// ---- antigen, antibody, memory ---------------------------------------
KAG    :   1.0  : Antigen release rate constant (AU/day)
KAGD   :   0.15 : Antigen elimination rate constant (1/day)
KASO   :  40    : ASO formation per unit antigen (Todd U/AU/day)
KASOD  :   0.0154 : ASO decay rate constant (1/day)
ASOBASE: 100    : Baseline ASO titre (Todd units)
KMEM   :   4e-3 : Memory formation rate (1/AU/day)
KMEMD  :   3.8e-4 : Memory decay rate constant (1/day)
MEMMAX :   1.0  : Memory pool ceiling (AU)
SUSC   :   1.0  : Host susceptibility multiplier (HLA-DR7 etc.)
KXAB   :   4e-3 : Cross-reactive antibody formation (1/AU/day)
BMEM   :   1.5  : Memory amplification of the antibody response (-)
KXABD  :   0.01155 : Cross-reactive antibody decay (1/day), t1/2 60 d
KVIT   :   0.60 : Valvulitis formation per unit antibody (1/day)
BSCAR  :   2.0  : Amplification of valvulitis by existing scar (-)
KVITD  :   0.0231 : Valvulitis resolution rate (1/day), t1/2 30 d

// ---- valve -------------------------------------------------------------
MVAN   :   4.5  : Normal mitral valve area (cm2)
MVAMIN :   0.30 : Floor on mitral valve area (cm2)
KFI    :   3.004e-4 : Immune-arm deposition (cm per AU-day of valvulitis)
KFS    :   9.132e-5 : Shear-arm deposition (cm/day per unit shear)
KHALT  :   1.6  : Rigidity brake on the shear arm (cm2)
RESTEN :   1.0  : Multiplier on the shear arm after balloon valvotomy (-)
KED    :   0.5  : Leaflet oedema formation (1/day)
KEDD   :   0.033: Leaflet oedema resolution (1/day), t1/2 21 d
KEDA   :   0.012: Effective area lost per unit oedema (cm2/AU)
KCA    :   4e-3 : Valve calcium accrual (AU/day)
KMRA   :   0.096: Acute mitral regurgitation formation (grade/AU/day)
KMRAD  :   0.03 : Acute mitral regurgitation resolution (1/day)
KMRC   :   1.2e-3 : Chronic (fixed) mitral regurgitation (grade/AU/day)
KAOI   :   6e-4 : Aortic valve involvement (grade/AU/day)

// ---- haemodynamics ------------------------------------------------------
HRB    :  72    : Baseline heart rate (beats/min)
DHRAF  :  35    : Heart-rate increment at full AF burden (beats/min)
HRMIN  :  45    : Floor on heart rate (beats/min)
HRSNS  :  15    : Maximum sympathetic heart-rate increment (beats/min)
KHRDEM :   0.5  : Fractional heart-rate rise per unit extra CO demand (-)
SVMAX  : 110    : Maximum stroke volume the ventricle can deliver (mL)
GORLIN :  37.7  : Gorlin constant for the mitral valve
AFPEN  :   0.12 : Effective-orifice penalty at full AF burden (-)
CO0    :   5.0  : Resting cardiac output demand (L/min)
DEMF   :   1.0  : Demand multiplier: 1 rest, 1.5 pregnancy, 2 exercise
LVEDP0 :   6.0  : Left ventricular end-diastolic pressure at euvolaemia
KVOLP  :   8.0  : LVEDP rise per unit volume excess (mmHg)
LAPMAX :  32    : Maximum sustainable left atrial pressure (mmHg)
LAPTH  :  18    : LA pressure above which congestion accumulates (mmHg)
KCON   :   0.10 : Congestion formation (1/mmHg/day)
KCOND  :   0.25 : Congestion resolution (1/day)
LAV0   :  55    : Normal left atrial volume (mL)
KLAV   :   0.06 : LA remodelling rate (mL/mmHg/day)
KLAVD  :   0.002: LA reverse remodelling rate (1/day)
LAPB   :  12    : LA pressure above which the atrium remodels (mmHg)
KAFO   :   3e-4 : Atrial fibrillation onset rate constant (1/day)
LAVREF :  60    : Reference LA volume excess for the AF hazard (mL)
NAF    :   2.0  : Exponent of the AF hazard in LA size (-)
KAFR   :   0    : Rhythm-control removal of AF burden (1/day)
PVRB   :   1.4  : Baseline pulmonary vascular resistance (Wood units)
KPVR   :   1.8e-3 : Pulmonary vascular remodelling (WU/mmHg/day)
PVRMAX :  12    : Ceiling on pulmonary vascular resistance (Wood units)
KPVRR  :   4e-3 : Pulmonary vascular reverse remodelling (1/day)
LAPPV  :  20    : LA pressure above which PVR rises (mmHg)
KRV    :   3e-4 : RV deterioration rate (1/mmHg/day)
PAPTH  :  35    : Mean PA pressure above which the RV fails (mmHg)
KRVR   :   2e-3 : RV recovery rate (1/day)
KVOLR  :   0.05 : Return of volume state to euvolaemia (1/day)
EDIUR  :   0    : Diuretic effect on the volume state (1/day)

// ---- rhythm, embolism, anticoagulation ------------------------------------
KEMB   :   1.85e-5 : Embolic hazard scale (1/day)
AF50E  :   0.15 : AF burden giving half-maximal embolic risk (-)
EMAXAC :   0.85 : Maximum anticoagulant risk reduction (-)
INRE50 :   1.55 : INR giving half-maximal risk reduction (-)

// ---- recurrence hazard ------------------------------------------------------
LAMEXP :   3.0  : Sore-throat exposures per year in this setting
PGASP  :   0.25 : P(GAS | sore throat) in this setting
PRHEUM :   0.0444 : P(ARF | untreated GAS in a susceptible host)
BMEM2  :   2.0  : Memory amplification of the recurrence hazard (-)
QARF   :   2.483: Antigen delivered by one full ARF episode (AU)
EPISODIC:  0    : 1 = discrete inoculations only, 0 = expected value

// ---- other drug pharmacokinetics ---------------------------------------------
KAASA  :  24    : Aspirin absorption rate constant (1/day)
VASA   :  12    : Salicylate volume of distribution (L)
CLASA  :  25    : Salicylate clearance (L/day)
EMAXASA:   0.55 : Maximum aspirin suppression of valvulitis (-)
EC50ASA: 120    : Salicylate EC50 (mg/L)
KAPRED :  24    : Prednisolone absorption rate constant (1/day)
VPRED  :  40    : Prednisolone volume of distribution (L)
CLPRED : 222    : Prednisolone clearance (L/day)
EMAXPR :   3.0  : Maximum prednisolone suppression of valvulitis (-)
EC50PR :   0.05 : Prednisolone EC50 (mg/L)
KABB   :  12    : Beta-blocker absorption rate constant (1/day)
VBB    : 250    : Beta-blocker volume of distribution (L)
CLBB   :1000    : Beta-blocker clearance (L/day)
EMAXBB :   0.35 : Maximum fractional heart-rate reduction (-)
EC50BB :   0.04 : Beta-blocker EC50 (mg/L)
KADIG  :  12    : Digoxin absorption rate constant (1/day)
VDIG   : 500    : Digoxin volume of distribution (L)
CLDIG  : 180    : Digoxin clearance (L/day)
FDIG   :   0.7  : Digoxin oral bioavailability (-)
EMAXDIG:   0.18 : Maximum digoxin heart-rate reduction (-)
EC50DIG:   1.0  : Digoxin EC50 (ng/mL)
KAW    :  12    : Warfarin absorption rate constant (1/day)
VW     :  10    : Warfarin volume of distribution (L)
CLW    :   4.2  : Warfarin clearance (L/day)
EMAXW  :   0.92 : Maximum inhibition of clotting factor synthesis (-)
IC50W  :   0.6  : Warfarin IC50 (mg/L)
KDEGP  :   0.35 : Prothrombin complex turnover (1/day)
GAMINR :   0.85 : Exponent relating factor activity to INR (-)

// ---- acute phase ----------------------------------------------------------------
KCRP   :  26.4  : CRP formation per unit valvulitis (mg/L/day)
KCRPD  :   0.35 : CRP elimination (1/day)
CRPB   :   2.0  : Baseline CRP (mg/L)

$CMT @annotated
BPGF   : Benzathine penicillin G, fast depot (mg penicillin-G equivalent)
BPGS   : Benzathine penicillin G, slow depot (mg)
PVA    : Oral penicillin absorption site (mg)
PENC   : Penicillin G, central compartment (mg)
CE     : Penicillin at the tonsillar effect site (ug/mL)
GAS    : Pharyngeal group A streptococcal burden (CFU-equivalents)
MIMM   : Mucosal / type-specific anti-GAS immunity (AU)
AG     : Systemic streptococcal antigen load (AU)
ASO    : Anti-streptolysin O titre (Todd units)
MEM    : Cross-reactive memory B/T pool (AU)
XAB    : Cross-reactive antibody (AU)
VIT    : Valvulitis activity (AU)
EDEM   : Acute leaflet oedema (AU)
MVA    : Mitral valve area (cm2)
CA     : Valve calcium (AU)
MRA    : Acute, partly reversible mitral regurgitation (grade)
MRC    : Chronic, fixed mitral regurgitation (grade)
AOI    : Aortic valve involvement (grade)
LAV    : Left atrial volume (mL)
AFB    : Atrial fibrillation burden (0-1)
PVR    : Pulmonary vascular resistance (Wood units)
RVF    : Right ventricular function index (1 = normal)
CONG   : Pulmonary congestion index (AU)
VOL    : Volume / preload state (1 = euvolaemic)
CRP    : C-reactive protein (mg/L)
ASAA   : Aspirin absorption site (mg)
ASAC   : Salicylate, central (mg)
PREDA  : Prednisolone absorption site (mg)
PREDC  : Prednisolone, central (mg)
BBA    : Beta-blocker absorption site (mg)
BBC    : Beta-blocker, central (mg)
DIGA   : Digoxin absorption site (ug)
DIGC   : Digoxin, central (ug)
WA     : Warfarin absorption site (mg)
WC     : Warfarin, central (mg)
PCF    : Prothrombin complex activity (fraction of normal)
CUMARF : Cumulative expected recurrent ARF episodes (count)
CUMEMB : Cumulative expected embolic events (count)
TPROT  : Cumulative days with plasma penicillin above threshold (days)

$MAIN
  // dose a BPG injection into BOTH depots with a single record by giving
  // two records (cmt 1 and cmt 2) with the same amt; F splits them
  F_BPGF = FFBPG;
  F_BPGS = FSBPG;

  MVA_0  = MVAN;
  ASO_0  = ASOBASE;
  LAV_0  = LAV0;
  PVR_0  = PVRB;
  RVF_0  = 1.0;
  VOL_0  = 1.0;
  PCF_0  = 1.0;
  CRP_0  = CRPB;

$ODE
  // =================================================================
  //  A.  ALGEBRAIC HAEMODYNAMIC BLOCK  (the arithmetic core)
  // =================================================================
  double MVAc  = (MVA > MVAMIN) ? MVA : MVAMIN;

  double CBB   = BBC  / VBB;                       // mg/L
  double CDIG  = DIGC / VDIG;                      // ug/L = ng/mL
  double EBB   = EMAXBB  * CBB  / (CBB  + EC50BB);
  double EDIG  = EMAXDIG * CDIG / (CDIG + EC50DIG);

  // Heart rate.  Note the SIGN of the congestion term: the sympathetic
  // response to pulmonary congestion RAISES the rate, and in mitral
  // stenosis a higher rate shortens diastole and makes the gradient
  // worse.  The reflex is maladaptive here, and the model shows it.
  double HRraw = (HRB * (1.0 + KHRDEM * (DEMF - 1.0))
                  + DHRAF * AFB
                  + HRSNS * CONG / (CONG + 1.0)) * (1.0 - EBB) * (1.0 - EDIG);
  double HR    = (HRraw > HRMIN) ? HRraw : HRMIN;

  double RR    = 60.0 / HR;                        // s per beat
  double SYS   = 0.36 * sqrt(RR);                  // systolic period, s
  double DFPB  = RR - SYS; if (DFPB < 0.05) DFPB = 0.05;
  double DFPM  = HR * DFPB;                        // diastolic s per minute

  // Effective orifice: acute oedema narrows it reversibly, and losing
  // the atrial contribution to filling behaves like a smaller orifice.
  double MVAo  = MVAc - KEDA * EDEM; if (MVAo < MVAMIN) MVAo = MVAMIN;
  double MVAeff= MVAo * (1.0 - AFPEN * AFB);

  double LVEDP = LVEDP0 + KVOLP * (VOL - 1.0);
  double dPmax = LAPMAX - LVEDP; if (dPmax < 1.0) dPmax = 1.0;

  // Three separate limits on cardiac output.  The min() is what makes
  // an optimal heart rate exist: COvalve falls with rate, COsv rises.
  double CODEM  = CO0 * DEMF;
  double COvalve= GORLIN * MVAeff * sqrt(dPmax) * DFPM / 1000.0;
  double COsv   = SVMAX * HR / 1000.0;
  double CO     = CODEM;
  if (COvalve < CO) CO = COvalve;
  if (COsv    < CO) CO = COsv;

  double MVF   = CO * 1000.0 / DFPM;               // mL/s
  double MVG   = pow(MVF / (GORLIN * MVAeff), 2.0);// mmHg, mean gradient
  double LAP   = LVEDP + MVG;
  double PAPm  = LAP + PVR * CO;

  double SXI   = 0.6 * CONG / (CONG + 1.2) + 0.4 * (1.0 - CO / CODEM);
  double NYHA  = 1.0 + 3.0 * SXI; if (NYHA > 4.0) NYHA = 4.0;

  // =================================================================
  //  B.  PENICILLIN
  // =================================================================
  double VPEN  = VPEN0  * (WT / 70.0);
  double CLPEN = CLPEN0 * pow(WT / 70.0, 0.75) * (CRCL / 100.0);
  double K10   = CLPEN / VPEN;

  dxdt_BPGF = -KAFBPG * BPGF;
  dxdt_BPGS = -KASBPG * BPGS;
  dxdt_PVA  = -KAORAL * PVA;
  dxdt_PENC = KAFBPG * BPGF + KASBPG * BPGS + KAORAL * PVA - K10 * PENC;

  double CPEN = PENC / VPEN;                       // ug/mL, PLASMA
  dxdt_CE     = KEO * (CPEN - CE);                 // tonsillar effect site

  // Bacterial kill follows the EFFECT SITE.  The 0.02 ug/mL programme
  // benchmark is a PLASMA number, so TPROT below is computed on plasma.
  double EPEN = pow(CE, HPEN) / (pow(CE, HPEN) + pow(EC50P, HPEN));

  // =================================================================
  //  C.  PHARYNGEAL INFECTION
  // =================================================================
  // The GAS/(GAS+GEXT) factor is an extinction floor.  Without it the
  // ODE never reaches zero, mucosal immunity wanes, and a single
  // pharyngitis relapses forever.
  dxdt_GAS  = KGROW * GAS * (1.0 - GAS / GMAX) * (GAS / (GAS + GEXT))
              - (KGCL + KMIK * MIMM) * GAS
              - KPKILL * EPEN * GAS;

  dxdt_MIMM = KMI * GAS / (GAS + GM50I) - KMD * MIMM;

  double agpres = GAS / (GAS + GM50);   // systemic antigen RELEASE

  // =================================================================
  //  D.  RECURRENCE HAZARD  (a product: no single term abolishes it)
  // =================================================================
  double HAZARF = LAMEXP * PGASP * PRHEUM * SUSC
                  * (1.0 - EPEN) * (1.0 + BMEM2 * MEM);      // per year
  dxdt_CUMARF = HAZARF / 365.0;

  // Expected-value channel: for decade-scale runs the valve sees the
  // EXPECTED antigen rather than sampled episodes.  Set EPISODIC = 1
  // and dose the GAS compartment to simulate discrete infections.
  double AGBG = (EPISODIC > 0.5) ? 0.0 : QARF * HAZARF / 365.0;

  // =================================================================
  //  E.  ANTIGEN, ANTIBODY, MEMORY
  // =================================================================
  dxdt_AG  = KAG * agpres + AGBG - KAGD * AG;
  dxdt_ASO = KASO * AG - KASOD * (ASO - ASOBASE);
  // MEM SATURATES.  Written unbounded, the loop
  // AG -> MEM -> recurrence hazard -> AG has gain 2.9 and diverges.
  dxdt_MEM = KMEM * AG * SUSC * (1.0 - MEM / MEMMAX) - KMEMD * MEM;
  dxdt_XAB = KXAB * AG * SUSC * (1.0 + BMEM * MEM) - KXABD * XAB;

  // =================================================================
  //  F.  ANTI-INFLAMMATORY DRUGS
  // =================================================================
  double CASA  = ASAC  / VASA;
  double CPRED = PREDC / VPRED;
  double EASA  = EMAXASA * CASA  / (CASA  + EC50ASA);
  double EPRED = EMAXPR  * CPRED / (CPRED + EC50PR);
  double ANTIINF = 1.0 + EASA + EPRED;

  dxdt_ASAA  = -KAASA * ASAA;
  dxdt_ASAC  =  KAASA * ASAA  - CLASA  / VASA  * ASAC;
  dxdt_PREDA = -KAPRED * PREDA;
  dxdt_PREDC =  KAPRED * PREDA - CLPRED / VPRED * PREDC;

  // =================================================================
  //  G.  VALVULITIS
  // =================================================================
  double SCARF = 1.0 - MVAc / MVAN; if (SCARF < 0.0) SCARF = 0.0;
  dxdt_VIT  = KVIT * XAB * (1.0 + BSCAR * SCARF) / ANTIINF - KVITD * VIT;
  dxdt_EDEM = KED * VIT - KEDD * EDEM;
  dxdt_CRP  = KCRP * VIT - KCRPD * (CRP - CRPB);

  // =================================================================
  //  H.  VALVE REMODELLING -- the two arms and the geometry
  // =================================================================
  double SH    = MVG / 4.0;      // = (mean transvalvular velocity)^2
  double BRAKE = (MVAc - MVAMIN) / (MVAc - MVAMIN + KHALT);
  double DEPO  = KFI * VIT + KFS * RESTEN * SH * BRAKE;
  double GEOM  = 2.0 * sqrt(M_PI * MVAc);          // dA = -2 sqrt(pi A) dH

  dxdt_MVA = (MVA > MVAMIN) ? (-GEOM * DEPO) : 0.0;
  dxdt_CA  = KCA * (1.0 - MVAc / MVAN) * (1.0 + AGE / 50.0);
  dxdt_MRA = KMRA * VIT - KMRAD * MRA;
  dxdt_MRC = KMRC * VIT;
  dxdt_AOI = KAOI * VIT;

  // =================================================================
  //  I.  CHAMBERS, RHYTHM, LUNG
  // =================================================================
  double dLAP = LAP - LAPB;  if (dLAP < 0.0) dLAP = 0.0;
  dxdt_LAV = KLAV * dLAP - KLAVD * (LAV - LAV0);

  double lavx = LAV - LAV0; if (lavx < 0.0) lavx = 0.0;
  double fAF  = pow(lavx / LAVREF, NAF);
  dxdt_AFB = KAFO * (1.0 - AFB) * fAF - KAFR * AFB;

  double dLPV = LAP - LAPPV; if (dLPV < 0.0) dLPV = 0.0;
  dxdt_PVR = KPVR * dLPV * (1.0 - PVR / PVRMAX) - KPVRR * (PVR - PVRB);

  double dPAP = PAPm - PAPTH; if (dPAP < 0.0) dPAP = 0.0;
  // damage proportional to RVF so the index decays TOWARDS zero, never through
  dxdt_RVF = -KRV * dPAP * RVF + KRVR * (1.0 - RVF);

  double dCON = LAP - LAPTH; if (dCON < 0.0) dCON = 0.0;
  dxdt_CONG = KCON * dCON - KCOND * CONG;
  dxdt_VOL  = KVOLR * (1.0 - VOL) - EDIUR * VOL;

  // =================================================================
  //  J.  ANTICOAGULATION AND EMBOLISM
  // =================================================================
  dxdt_WA = -KAW * WA;
  dxdt_WC =  KAW * WA - CLW / VW * WC;
  double CW  = WC / VW;
  double IW  = EMAXW * CW / (CW + IC50W);
  dxdt_PCF = KDEGP * (1.0 - IW) - KDEGP * PCF;

  double PCFc = (PCF > 1e-3) ? PCF : 1e-3;
  double INR  = pow(1.0 / PCFc, GAMINR);
  double dINR = INR - 1.0; if (dINR < 0.0) dINR = 0.0;
  double EANTI= EMAXAC * dINR / (dINR + (INRE50 - 1.0));

  // AFB**0.7 was the first form here; its derivative is infinite at
  // AFB = 0, which is where every simulation starts.  A saturating Hill
  // has the same shape with a finite Jacobian.
  double EMBAF = AFB / (AFB + AF50E);
  dxdt_CUMEMB = KEMB * (0.15 + EMBAF)
                * pow(LAV / LAV0, 1.2)
                * pow(MVAN / MVAc, 0.5)
                * (1.0 - EANTI);

  // =================================================================
  //  K.  RATE-CONTROL PHARMACOKINETICS
  // =================================================================
  dxdt_BBA  = -KABB * BBA;
  dxdt_BBC  =  KABB * BBA - CLBB / VBB * BBC;
  dxdt_DIGA = -KADIG * DIGA;
  dxdt_DIGC =  KADIG * DIGA - CLDIG / VDIG * DIGC;

  // =================================================================
  //  L.  TIME ABOVE THE PROTECTIVE THRESHOLD  (the quantity that protects)
  // =================================================================
  dxdt_TPROT = pow(CPEN, 6.0) / (pow(CPEN, 6.0) + pow(MICP, 6.0));

$TABLE
  // ---- recompute the algebra for reporting ------------------------------
  double oMVA   = (MVA > MVAMIN) ? MVA : MVAMIN;
  double oCBB   = BBC / VBB;
  double oCDIG  = DIGC / VDIG;
  double oEBB   = EMAXBB  * oCBB  / (oCBB  + EC50BB);
  double oEDIG  = EMAXDIG * oCDIG / (oCDIG + EC50DIG);
  double oHRraw = (HRB * (1.0 + KHRDEM * (DEMF - 1.0)) + DHRAF * AFB
                   + HRSNS * CONG / (CONG + 1.0)) * (1.0 - oEBB) * (1.0 - oEDIG);
  double oHR    = (oHRraw > HRMIN) ? oHRraw : HRMIN;
  double oRR    = 60.0 / oHR;
  double oDFPB  = oRR - 0.36 * sqrt(oRR); if (oDFPB < 0.05) oDFPB = 0.05;
  double oDFPM  = oHR * oDFPB;
  double oMVAo  = oMVA - KEDA * EDEM; if (oMVAo < MVAMIN) oMVAo = MVAMIN;
  double oMVAeff= oMVAo * (1.0 - AFPEN * AFB);
  double oLVEDP = LVEDP0 + KVOLP * (VOL - 1.0);
  double odPmax = LAPMAX - oLVEDP; if (odPmax < 1.0) odPmax = 1.0;
  double oCODEM = CO0 * DEMF;
  double oCOv   = GORLIN * oMVAeff * sqrt(odPmax) * oDFPM / 1000.0;
  double oCOsv  = SVMAX * oHR / 1000.0;
  double oCO    = oCODEM; if (oCOv < oCO) oCO = oCOv; if (oCOsv < oCO) oCO = oCOsv;
  double oMVF   = oCO * 1000.0 / oDFPM;
  double oMVG   = pow(oMVF / (GORLIN * oMVAeff), 2.0);
  double oLAP   = oLVEDP + oMVG;
  double oPAP   = oLAP + PVR * oCO;
  double oSXI   = 0.6 * CONG / (CONG + 1.2) + 0.4 * (1.0 - oCO / oCODEM);
  double oNYHA  = 1.0 + 3.0 * oSXI; if (oNYHA > 4.0) oNYHA = 4.0;

  double oVPEN  = VPEN0 * (WT / 70.0);

  // Wilkins / Abascal echo score, 4-16.  Above 8 predicts a poor balloon result.
  double narrow = (MVAN - oMVA) / (MVAN - 0.8);
  if (narrow < 0.0) narrow = 0.0; if (narrow > 1.0) narrow = 1.0;
  double c8  = CA / (CA + 8.0);
  double c10 = CA / (CA + 10.0);
  double smob = 0.35 * narrow + 0.45 * c8;  if (smob > 1.0) smob = 1.0;
  double sthk = 0.30 * narrow + 0.50 * c8;  if (sthk > 1.0) sthk = 1.0;
  double ssub = 0.30 * narrow + 0.45 * c10; if (ssub > 1.0) ssub = 1.0;
  double oWILK = (1.0 + 3.0 * smob) + (1.0 + 3.0 * sthk)
               + (1.0 + 3.0 * c10)  + (1.0 + 3.0 * ssub);

  double oPCF  = (PCF > 1e-3) ? PCF : 1e-3;

  capture CPEN   = PENC / oVPEN;      // plasma penicillin G, ug/mL
  capture CEFF   = CE;                // effect-site penicillin, ug/mL
  capture PROT   = (CPEN > MICP) ? 1.0 : 0.0;   // above threshold now?
  capture MVAcm2 = oMVA;
  capture MVAEFF = oMVAeff;
  capture MVGRAD = oMVG;
  capture LAPRES = oLAP;
  capture PAPMN  = oPAP;
  capture HRBPM  = oHR;
  capture DFPMIN = oDFPM;
  capture COLMIN = oCO;
  capture NYHACL = oNYHA;
  capture WILKIN = oWILK;
  capture INRV   = pow(1.0 / oPCF, GAMINR);
  capture MRTOT  = MRA + MRC;
  capture CRPMGL = CRP;
  capture ASOT   = ASO;
  capture AFBURD = AFB;
  capture LOGGAS = (GAS > 1.0) ? log10(GAS) : 0.0;
'

mod <- mcode("rhd", code)

## =====================================================================
##  DOSING
##  Event tables are built as plain NM-TRAN data frames (ID, time, amt,
##  cmt, evid) rather than with ev() helpers, so that the compartment
##  numbering below is explicit and nothing depends on argument-matching
##  behaviour that varies between mrgsolve versions.
##
##  Compartment numbers, in $CMT order:
##    1 BPGF   2 BPGS   3 PVA    4 PENC   5 CE     6 GAS    7 MIMM
##    8 AG     9 ASO   10 MEM   11 XAB   12 VIT   13 EDEM  14 MVA
##   15 CA    16 MRA   17 MRC   18 AOI   19 LAV   20 AFB   21 PVR
##   22 RVF   23 CONG  24 VOL   25 CRP   26 ASAA  27 ASAC  28 PREDA
##   29 PREDC 30 BBA   31 BBC   32 DIGA  33 DIGC  34 WA    35 WC
##   36 PCF   37 CUMARF 38 CUMEMB 39 TPROT
## =====================================================================

CMT <- c(BPGF = 1, BPGS = 2, PVA = 3, GAS = 6, MVA = 14,
         ASAA = 26, PREDA = 28, BBA = 30, DIGA = 32, WA = 34)

rec <- function(time, amt, cmt, evid = 1L)
  data.frame(ID = 1L, time = time, amt = amt, cmt = as.integer(cmt),
             evid = as.integer(evid))

## One benzathine penicillin G injection is TWO records at the same time,
## one into each depot; F_BPGF = 0.30 and F_BPGS = 0.70 split the dose.
## 1.2 MU ~ 720 mg penicillin-G equivalent (360 mg if body weight < 27 kg).
## `keep` selects which scheduled doses were actually given, so that
## imperfect adherence is expressed as missing records, not as a fudge.
bpg <- function(interval = 28, n = 13, mg = 720, start = 0, keep = NULL) {
  i <- seq_len(n); if (!is.null(keep)) i <- i[keep]
  tt <- start + (i - 1) * interval
  rbind(rec(tt, mg, CMT[["BPGF"]]), rec(tt, mg, CMT[["BPGS"]]))
}

## Oral penicillin V 500 mg twice daily for 10 days (F ~ 0.5 -> 250 mg)
penv  <- function(start = 2, days = 10)
  rec(seq(start, start + days, by = 0.5), 250, CMT[["PVA"]])

## A pharyngeal inoculum, for the episodic channel (set EPISODIC = 1)
inoc  <- function(time = 0, n = 1e4) rec(time, n, CMT[["GAS"]])

asa   <- function(start, days, mg = 1000)      # aspirin 1 g every 6 h
  rec(seq(start, start + days, by = 0.25), mg, CMT[["ASAA"]])
pred  <- function(start, days, mg = 40)        # prednisolone once daily
  rec(seq(start, start + days, by = 1), mg, CMT[["PREDA"]])
bblok <- function(start, days, mg = 50)        # metoprolol once daily
  rec(seq(start, start + days, by = 1), mg, CMT[["BBA"]])
warf  <- function(start, days, mg = 5)         # warfarin once daily
  rec(seq(start, start + days, by = 1), mg, CMT[["WA"]])

## A balloon valvotomy is an instantaneous ADDITION to the valve-area
## compartment.  The step change in the shear coefficient that follows it
## is applied by running the post-procedure period as a separate call.
pmbv  <- function(time, gain) rec(time, gain, CMT[["MVA"]])

## order records by time, as mrgsolve expects
mk <- function(...) {
  d <- do.call(rbind, list(...))
  d[order(d$time, d$cmt), ]
}

run <- function(m, data = NULL, end = 365, delta = 1) {
  if (is.null(data)) mrgsim(m, end = end, delta = delta)
  else mrgsim(m, data = data, end = end, delta = delta)
}


## =====================================================================
##  SCENARIOS
##  Sixteen runs.  Each is the numerical form of a claim made in
##  README.md, and each is re-derived independently in
##  rhd_verify_python.py, where the printed numbers below come from.
## =====================================================================

## --- 1.  A single 1.2 MU injection, by body weight ---------------------
##  Days above 0.02 ug/mL: 23.4 (40 kg), 18.1 (70 kg), 14.8 (100 kg).
##  NO adult weight covers a 28-day interval on one injection.
s01 <- lapply(c(40, 55, 70, 85, 100), function(w)
  run(param(mod, WT = w), bpg(n = 1), end = 45, delta = 0.1))

## --- 2.  Four-weekly versus three-weekly in the same patient -----------
##  70 kg: 70.0% of calendar time protected at 28 d, 94%+ at 21 d.
s02a <- run(param(mod, WT = 70), bpg(28, 40), end = 1095)
s02b <- run(param(mod, WT = 70), bpg(21, 53), end = 1095)

## --- 3.  THE RANK INVERSION -------------------------------------------
##  80% of doses given, 55 kg, 3-weekly -> 79.0% of time protected
## 100% of doses given, 95 kg, 4-weekly -> 60.0% of time protected
##  The audited number and the protective number disagree on the ordering.
s03a <- run(param(mod, WT = 55), bpg(21, 53, keep = which(seq_len(53) %% 5 != 0)),
            end = 1095)
s03b <- run(param(mod, WT = 95), bpg(28, 40), end = 1095)

## --- 4.  Untreated streptococcal pharyngitis ---------------------------
##  GAS peaks day 5, clears by day 12; ASO peaks day 22 at ~608 Todd units.
s04 <- run(param(mod, EPISODIC = 1), inoc(), end = 60, delta = 0.25)

## --- 5.  The nine-day window ------------------------------------------
##  Treated on day 2: cumulative antigen falls 94%.  Treated on day 8 in
##  this fast-clearing host: only 7%, because the organism has already
##  gone.  In a slow-clearing (carrier) host the day-8 benefit is 76%.
##  The window is therefore a statement about CARRIAGE, not about the drug.
s05a <- run(param(mod, EPISODIC = 1), mk(inoc(), penv(2)), end = 60, delta = 0.25)
s05b <- run(param(mod, EPISODIC = 1), mk(inoc(), penv(8)), end = 60, delta = 0.25)
s05c <- run(param(mod, EPISODIC = 1, KMI = 0.08), mk(inoc(), penv(8)),
            end = 60, delta = 0.25)

## --- 6.  A first attack of acute rheumatic fever, with aspirin --------
##  Peak CRP 62 mg/L, peak acute MR grade 2.3, 0.36 cm2 of valve lost.
s06 <- run(param(mod, EPISODIC = 1), mk(inoc(), asa(19, 42)), end = 400)

## --- 7.  Steroid duration ---------------------------------------------
##  Valve area saved: 5% (2 wk), 14% (6 wk), 27% (12 wk), 42% (26 wk).
##  The antibody half-life is 60 d, so a 6-week course covers only 38%
##  of the antibody-time integral -- which is why the trials are null,
##  and what would have to change for them not to be.
s07 <- lapply(c(0, 2, 6, 12, 26), function(wk)
  run(param(mod, EPISODIC = 1),
      if (wk > 0) mk(inoc(), asa(19, 42), pred(19, wk * 7))
      else        mk(inoc(), asa(19, 42)),
      end = 400))

## --- 8.  Twenty-five years, prophylaxis on and off --------------------
##  Recurrent ARF 0.85 -> 0.05 episodes (-94%), but valve area only
##  3.28 -> 3.69 cm2 (+0.41).  Two different endpoints, two different
##  effect sizes, one drug.
latent <- function(m) init(m, MVA = 3.6, MEM = 1.0)
s08a <- run(latent(param(mod, WT = 55)), NULL,          end = 9125, delta = 5)
s08b <- run(latent(param(mod, WT = 55)), bpg(28, 327),  end = 9125, delta = 5)
s08c <- run(latent(param(mod, WT = 55)), bpg(21, 435),  end = 9125, delta = 5)

## --- 9.  Latent RHD over two years (the GOAL population) --------------
##  This is where the valve benefit lives: above the 2.75 cm2 crossover.
s09a <- run(init(param(mod, WT = 40), MVA = 4.2, MEM = 1.0), NULL, end = 730, delta = 2)
s09b <- run(init(param(mod, WT = 40), MVA = 4.2, MEM = 1.0), bpg(28, 27),
            end = 730, delta = 2)

## --- 10. Pregnancy at a FIXED 1.5 cm2 orifice -------------------------
##  Mean gradient 6.4 -> 23.0 mmHg, LA pressure 12.4 -> 29.0, NYHA 1 -> 2.4.
##  The valve did not change.  The demand did.
s10a <- run(init(mod, MVA = 1.5), NULL, end = 200)
s10b <- run(init(param(mod, DEMF = 1.5), MVA = 1.5), NULL, end = 200)

## --- 11. New-onset AF at the same orifice, and rate control -----------
##  Gradient 6.4 -> 11.9 mmHg on the same valve; metoprolol recovers part.
s11a <- run(init(mod, MVA = 1.5, AFB = 1.0), NULL, end = 200)
s11b <- run(init(mod, MVA = 1.5, AFB = 1.0), bblok(0, 200, 100), end = 200)

## --- 12. At the LA-pressure ceiling, rate control buys OUTPUT ---------
##  Severe MS, permanent AF, exertion: the gradient is pinned at
##  LAPMAX - LVEDP and CANNOT fall, so the diastole beta blockade buys is
##  spent on forward flow.  The measured gradient is the one number that
##  does not move while the patient improves.  Below the ceiling
##  (MVA 1.4, at rest) the same drug lowers the gradient instead.
s12a <- run(init(param(mod, DEMF = 1.6), MVA = 1.0, AFB = 1.0, LAV = 130),
            NULL, end = 30, delta = 0.5)
s12b <- run(init(param(mod, DEMF = 1.6), MVA = 1.0, AFB = 1.0, LAV = 130),
            bblok(0, 30, 200), end = 30, delta = 0.5)
s12c <- run(init(mod, MVA = 1.4, AFB = 1.0, LAV = 110), NULL, end = 30, delta = 0.5)
s12d <- run(init(mod, MVA = 1.4, AFB = 1.0, LAV = 110), bblok(0, 30, 200),
            end = 30, delta = 0.5)

## --- 13. Anticoagulation in rheumatic AF ------------------------------
##  Embolic events 0.116/yr untreated; warfarin at INR ~2.2 cuts 65%.
s13a <- run(init(mod, MVA = 1.0, AFB = 1.0, LAV = 140), NULL, end = 365)
s13b <- run(init(mod, MVA = 1.0, AFB = 1.0, LAV = 140), warf(0, 365, 5), end = 365)

## --- 14. Balloon valvotomy and restenosis -----------------------------
##  0.95 -> 2.10 cm2 immediately, 1.57 at 5 years, 0.81 at 10 years.
##  RESTEN steps to 1.35 after the procedure, so the run is split in two
##  and the second call starts from the post-procedure state.
s14pre  <- run(init(mod, MVA = 0.95, CA = 3.0, LAV = 85), NULL, end = 1, delta = 1)
s14post <- run(init(param(mod, RESTEN = 1.35),
                    MVA = 0.95 + 1.15, CA = 3.0, LAV = 85),
               bpg(28, 131), end = 3649, delta = 5)

## --- 15. Untreated severe mitral stenosis over eight years ------------
##  PVR 1.4 -> 4.7 Wood units, mean PA pressure 39 mmHg, RV index 0.58.
s15 <- run(init(mod, MVA = 0.9, LAV = 120, AFB = 0.8), NULL, end = 2920, delta = 5)

## --- 16. Host susceptibility: the same exposure, two hosts ------------
s16a <- run(init(param(mod, SUSC = 0.4), MVA = 3.6, MEM = 1.0), NULL,
            end = 9125, delta = 10)
s16b <- run(init(param(mod, SUSC = 2.0), MVA = 3.6, MEM = 1.0), NULL,
            end = 9125, delta = 10)


## =====================================================================
##  READING THE OUTPUT
##
##  TPROT / time            fraction of calendar time protected -- compare
##                          this with doses-given, not with each other
##  MVAcm2                  mitral valve area; the ratchet only turns down
##  MVGRAD, LAPRES          mean gradient and LA pressure; these move when
##                          DEMAND moves, at a completely unchanged orifice
##  DFPMIN                  diastolic seconds per minute -- the quantity
##                          rate control actually buys
##  CUMARF                  the endpoint prophylaxis moves a lot
##  MVAcm2 trajectory       the endpoint it moves much less once stenosis
##                          exists; the two are not the same claim
##  WILKIN                  Wilkins score; > 8 argues against a balloon
##  CUMEMB                  expected embolic events
##
##  CALIBRATION AND ITS LIMITS.  Fitted or anchored: BPG depot release
##  (Kaplan 1989), penicillin size scaling (Hand 2019, Neely 2014),
##  progression of established MS (Sagie 1996, Gordon 1992), Gorlin
##  constants (Gorlin 1951, Hakki 1981), gradient/severity grades
##  (Baumgartner 2009), recurrence reduction on BPG (Manyemba 2002),
##  embolic rate in rheumatic AF, warfarin dose-INR, PMBV gain and
##  5-year area (Iung 1999, Palacios 2002), reactive PH in severe MS.
##  PREDICTIONS, not fitted: the adherence/protection rank inversion,
##  the 2.75 cm2 crossover, the area-dependent optimal heart rate, the
##  steroid duration-response, and the claim that at the LA-pressure
##  ceiling rate control raises output while leaving the measured
##  gradient untouched.
##
##  MOST EXPOSED CLAIM.  The shear arm is calibrated to a cohort MEAN
##  progression of 0.09 cm2/yr at MVA 1.5 and then accelerates as the
##  valve narrows.  Sagie 1996 reported the opposite association --
##  progression was SLOWER in tighter valves.  The model's 10-year
##  post-valvotomy area (0.8 cm2) is correspondingly more pessimistic
##  than the ~40%-restenosis literature.  A sensitivity in
##  rhd_verify_python.py shows this is insensitive to how the shear arm
##  depends on area, because the gradient saturates at LAPMAX - LVEDP
##  once the valve is tight; so if it is wrong, the error is in reading
##  a cohort mean as if every valve progressed.
## =====================================================================
