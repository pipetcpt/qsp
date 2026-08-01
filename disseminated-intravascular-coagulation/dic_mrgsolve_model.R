## =====================================================================
##  dic_mrgsolve_model.R
##  Disseminated Intravascular Coagulation (DIC) — QSP / PK-PD model
##  49 ODE compartments · 15 therapeutic scenarios
##
##  파종성 혈관내 응고 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  DIC is ONE consumption process (tissue-factor-driven thrombin
##  generation) read through TWO independently-set clocks:
##
##   CLOCK 1  the SIGN of the acute-phase term on each protein's synthesis
##            fibrinogen : synthesis x (1 + APOSF * APF)      <- POSITIVE
##            AT, PC, PS : synthesis / (1 + ANEG*  * APF)     <- NEGATIVE
##            Sepsis sets APF high (IL-6 storm); APL leaves it near zero.
##            The SAME consumption rate therefore leaves fibrinogen at
##            281 mg/dL in sepsis and 87 mg/dL in APL.  Thrombotic vs
##            haemorrhagic phenotype is a CONSEQUENCE, not an input.
##
##   CLOCK 2  the PAI-1 : t-PA MOLAR ratio
##            free tPA = tPA_nM * (FLOC + (1-FLOC)/(1 + PAI1_nM/KPAI))
##            Sepsis: 42:1 molar excess -> free tPA falls to 0.66x of the
##            healthy value even though total tPA rose 4-fold  (SHUTDOWN).
##            APL   :  3.9:1 plus annexin-A2 cofactor -> 13x sepsis
##                     (HYPERFIBRINOLYSIS).
##            Tranexamic acid enters ONE equation with ONE sign and its
##            clinical benefit changes sign automatically.
##
##  THREE PRODUCTS THAT EXPLAIN THE TRIALS
##  --------------------------------------
##   (1) Heparin is a CATALYST.  Its effect is
##         KATII * (AT + (AMP-1) * ATcat),  ATcat = MM-saturating in AT.
##       At AT 40% the achievable ceiling is 72% of the AT-100% ceiling —
##       a ceiling loss, not a potency loss, so a dose increase cannot
##       recover it.  Conversely raising AT to 220% buys only +26% on the
##       heparin arm but 3.7x on the UNCATALYSED arm — which is why
##       antithrombin concentrate helps patients who are NOT on heparin
##       (KyberSept subgroup) far more than those who are.
##
##   (2) APC generation = KPCA x THR x TM x EPCR x PC — a product of FOUR
##       terms.  In simulated septic DIC TM 0.318 x EPCR 0.559 x PC 0.333
##       = 0.059, so APC per unit thrombin falls 17-fold.  That is the
##       pharmacological logic of giving APC itself (drotrecogin alfa).
##
##   (3) Haemostatic failure is a PRODUCT of independent deficits, so the
##       WORST term dominates: correcting fibrinogen while platelets are
##       30 x10^9/L moves the bleeding index almost not at all.
##
##  Calibration anchors (see dic_references.md for sources) and the
##  achieved values of this parameter set are listed at the bottom of
##  this file under CALIBRATION NOTES.
##
##  Requires: mrgsolve (>= 1.0), dplyr, ggplot2, tidyr
## =====================================================================

library(mrgsolve)
library(dplyr)

dic_code <- '

$PROB
# Disseminated Intravascular Coagulation — QSP model (49 ODEs)

$PARAM @annotated
// ---------------- trigger ----------------
KCLRP   : 0.030   : pathogen clearance rate (1/h)
SRCCTL  : 1.0     : source-control multiplier on pathogen clearance (-)
KHP     : 0.30    : half-saturating pathogen load for the inflammatory hazard (-)
KGB     : 0.010   : APL blast growth rate (1/h)
BMAX    : 120     : APL blast carrying capacity (10^9/L)
KDIFF   : 0.020   : ATRA-driven blast differentiation rate (1/h)
EC50ATRA: 0.05    : ATRA EC50 for differentiation (mg/L)
EMAXATRA: 1.0     : ATRA maximal differentiation effect (-)
TFEXO   : 0.0     : direct tissue-factor input (trauma/obstetric) (1/h)

// ---------------- cytokines ----------------
KDIL6   : 0.10    : IL-6 turnover (1/h)
IL6B    : 5.0     : IL-6 baseline (pg/mL)
KSIL6   : 4000    : IL-6 gain per unit pathogen load (pg/mL)
SIL6    : 1.0     : aetiology-specific IL-6 gain; 1 sepsis, 0.10 APL (-)
KIL6BL  : 6.0     : IL-6 gain per blast unit (pg/mL per 10^9/L)
KDTNF   : 0.25    : TNF turnover (1/h)
TNFB    : 10.0    : TNF baseline (pg/mL)
KSTNF   : 300     : TNF gain per unit pathogen load (pg/mL)
KTNFBL  : 1.0     : TNF gain per blast unit (pg/mL per 10^9/L)
KTNF    : 120     : TNF half-effect increment (pg/mL)
KIL6    : 400     : IL-6 half-effect increment; sets the acute-phase factor (pg/mL)

// ---------------- NETs / histones ----------------
KSNET   : 1.2     : NET release rate (ug/mL/h)
KDNET   : 0.09    : cf-DNA clearance (1/h)
NETB    : 0.30    : cf-DNA baseline (ug/mL)
KHNET   : 0.55    : histone release per unit NET (1/h)
KDHIST  : 0.15    : histone clearance (1/h)
HISTB   : 0.10    : histone baseline (ug/mL)
KHIST   : 8.0     : histone half-effect (ug/mL)

// ---------------- tissue factor ----------------
KDTF    : 0.10    : tissue-factor turnover (1/h)
TFB     : 0.02    : constitutive tissue-factor activity (nM-eq)
KSTF    : 0.55    : monocyte/endothelial TF induction gain (nM-eq/h)
KTFBL   : 0.010   : blast-derived TF gain (nM-eq/h per 10^9/L)

// ---------------- endothelium ----------------
KRTM    : 0.020   : thrombomodulin recovery (1/h)
KSHTM   : 0.040   : thrombomodulin shedding gain (1/h)
KREP    : 0.030   : EPCR recovery (1/h)
KSHEP   : 0.022   : EPCR shedding gain (1/h)
KRGL    : 0.025   : glycocalyx recovery (1/h)
KSHGL   : 0.038   : glycocalyx shedding gain (1/h)

// ---------------- thrombin generation ----------------
KXA     : 9.0     : extrinsic tenase FXa generation gain (nM/h)
RTFPI   : 1.8     : TFPI inhibition coefficient on TF-FVIIa (-)
KTEN    : 0.80    : intrinsic tenase amplification gain (nM/h)
KTHRF   : 1.2     : thrombin half-effect for FXI/FVIII feedback (nM)
PSMV    : 0.40    : microvesicle-derived procoagulant surface weight (-)
KTFS    : 0.25    : TF half-effect for the microvesicle surface term (nM-eq)
PLT0    : 250     : platelet reference count (10^9/L)
KPRO    : 12.0    : prothrombinase gain (1/h per nM)
KDXA    : 3.0     : AT-independent FXa clearance (1/h)
KATXA   : 1.6     : AT-mediated FXa inactivation (1/h)
KATII   : 2.2     : AT-mediated thrombin inactivation (1/h)
KDTHR   : 1.5     : AT-independent thrombin clearance (1/h)
KTMB    : 1.2     : thrombomodulin-bound thrombin removal (1/h)
EMAXHX  : 10.0    : maximal heparin anti-Xa rate enhancement (-)
EMAXHI  : 12.0    : maximal heparin anti-IIa rate enhancement (-)
EC50H   : 0.35    : heparin EC50 (anti-Xa IU/mL)
KMAT    : 0.25    : Michaelis constant for AT on heparin, in AT-fraction units (-)
KHEPN   : 25.0    : histone/PF4 concentration neutralising half the heparin (ug/mL)
KDTI50  : 0.45    : argatroban EC50 (ug/mL)
EMAXDTI : 9.0     : argatroban maximal thrombin-clearance enhancement (-)
KIAPC   : 0.30    : APC concentration halving FVa/FVIIIa cofactor activity (nM)

// ---------------- the common consumption clock ----------------
KC0     : 0.032   : maximal fractional consumption rate (1/h)
KTHRC   : 0.70    : thrombin half-effect for consumption (nM)
WFIB    : 1.00    : fibrinogen consumption weight (-)
WF2     : 0.95    : prothrombin consumption weight (-)
WF5     : 1.35    : factor V consumption weight (-)
WF7     : 0.55    : factor VII consumption weight (-)
WF8     : 1.30    : factor VIII consumption weight (-)
WF10    : 0.95    : factor X consumption weight (-)
WPC     : 0.45    : protein C consumption weight (-)
WPS     : 0.50    : protein S consumption weight (-)
WTFPI   : 1.20    : TFPI consumption weight (-)
KATC    : 0.00060 : AT consumed per unit thrombin-inactivation flux (1/nM)

// ---------------- synthesis / acute phase (CLOCK 1) ----------------
THFIB   : 100     : fibrinogen half-life (h)
FIB0    : 300     : fibrinogen baseline (mg/dL)
APOSF   : 3.0     : maximal IL-6-driven rise in fibrinogen synthesis (-)
THF2    : 65      : prothrombin half-life (h)
THF5    : 15      : factor V half-life (h)
THF7    : 5       : factor VII half-life (h)
THF8    : 10      : factor VIII half-life (h)
APOS8   : 1.5     : acute-phase rise in factor VIII synthesis (-)
THF10   : 40      : factor X half-life (h)
THAT    : 65      : antithrombin half-life (h)
ANEGAT  : 0.20    : acute-phase suppression of AT synthesis (-)
THPC    : 6       : protein C half-life (h)
ANEGPC  : 1.20    : acute-phase suppression of protein C synthesis (-)
THPS    : 42      : free protein S half-life (h)
ANEGPS  : 1.00    : acute-phase suppression of free protein S, via C4BP (-)
THTFPI  : 3       : TFPI half-life (h)
KLEAK   : 0.015   : capillary-leak loss of AT/PC/PS at zero glycocalyx (1/h)
KELAT   : 0.002   : elastase-mediated AT degradation (1/h)
KELPC   : 0.004   : elastase-mediated protein C / TFPI degradation (1/h)

// ---------------- protein C pathway ----------------
KPCA    : 1.10    : protein C activation gain (1/h per nM thrombin)
KDAPC   : 1.73    : APC clearance (1/h)
KAPCV   : 0.020   : residual APC-mediated factor V loss (1/h per nM)
KAPC8   : 0.040   : residual APC-mediated factor VIII loss (1/h per nM)
KAPCPAI : 0.10    : APC-mediated PAI-1 inhibition (1/h per nM)
EMTM    : 1.30    : maximal thrombomodulin restored by rTM (fraction)
EC50TM  : 0.40    : rTM EC50 (ug/mL)
FTAFIS  : 0.05    : soluble-rTM efficiency for TAFI activation vs membrane TM (-)
EHMGB   : 0.45    : rTM lectin-domain neutralisation of histone/HMGB1 shedding (-)

// ---------------- platelets ----------------
KDPLT   : 0.0034657 : platelet senescent turnover (1/h)
ETPO    : 0.80    : TPO-driven compensatory production exponent (-)
MSUP    : 1.0     : marrow output multiplier; 0.20 in APL (-)
KPC0    : 0.045   : maximal thrombin-driven platelet consumption (1/h)
KTHRP   : 0.80    : thrombin half-effect on platelets (nM)
KPFD    : 0.010   : platelet consumption onto deposited fibrin (1/h per AU)
KACT    : 1.60    : platelet activation rate (1/h)
KDEACT  : 0.80    : platelet de-activation rate (1/h)
AHISTP  : 1.20    : histone weight in platelet activation (-)
WTHRP   : 0.35    : thrombin weight in platelet activation (-)

// ---------------- fibrin: the RATE and the STOCK ----------------
KFPA    : 0.85    : fibrinopeptide release gain (1/h per nM)
KSFMD   : 2.5     : soluble fibrin monomer clearance (1/h)
KLYSS   : 1.6     : plasmin-mediated soluble fibrin lysis (1/h per nM)
KDEP    : 0.055   : maximal microvascular fibrin deposition rate (AU/h)
KTHD    : 0.30    : thrombin at half-maximal deposition (nM)
NHD     : 1.5     : Hill coefficient for deposition (-)
KRES    : 0.008   : plasmin-independent (RES) fibrin clearance (1/h)
KLYS    : 0.150   : plasmin-mediated deposited-fibrin lysis (1/h per nM)
ATAFI   : 2.5     : TAFIa potency in blocking plasmin action on fibrin (-)

// ---------------- fibrinolysis (CLOCK 2) ----------------
KDTPA   : 0.35    : t-PA clearance (1/h)
TPAB    : 5.0     : t-PA baseline (ng/mL)
KSTPA   : 5.0     : t-PA release gain (ng/mL/h)
AENDO   : 1.0     : glycocalyx-loss weight in t-PA release (-)
KTPABL  : 0.05    : blast-derived t-PA/u-PA release (ng/mL/h per 10^9/L)
KDPAI   : 0.20    : PAI-1 clearance (1/h)
PAIB    : 20.0    : PAI-1 baseline (ng/mL)
KSPAI   : 140     : PAI-1 induction gain (ng/mL/h)
SPAI    : 1.0     : aetiology-specific PAI-1 gain; 1 sepsis, 0.15 APL (-)
PAIACT  : 0.60    : active (non-latent) fraction of PAI-1 (-)
FLOC    : 0.05    : fibrin-surface-protected fraction of t-PA (-)
KPAI    : 0.20    : PAI-1 concentration halving free t-PA (nM)
KANX    : 3.0     : maximal annexin-A2 amplification of t-PA on blasts (-)
KBL50   : 6.0     : blast burden at half-maximal annexin-A2 effect (10^9/L)
KPG     : 95      : plasminogen activation gain (1/h per nM free t-PA)
KFS     : 0.35    : fibrin surface at half-maximal plasmin generation (AU)
KUPA    : 0.16    : blast u-PA fibrin-INDEPENDENT plasmin generation (nM/h)
KPLGC   : 0.0005  : plasminogen consumed per unit plasmin generated (-)
KDPLN   : 1.2     : plasmin clearance (1/h)
KA2AP   : 5.0     : alpha-2-antiplasmin neutralisation of plasmin (1/h)
CA2     : 0.030   : alpha-2-antiplasmin consumption per unit plasmin (1/h per nM)
KDA2    : 0.030   : alpha-2-antiplasmin resynthesis (1/h)
KFGN    : 0.020   : direct fibrinogenolysis gain, gated on a2AP depletion (1/h per nM)
IC50TXA : 6.0     : tranexamic acid IC50 on plasmin generation (ug/mL)
KTAFI   : 7.63    : TAFI activation gain (1/h per nM thrombin)
KDTAFI  : 4.16    : TAFIa decay, t1/2 ~10 min (1/h)
KDDG    : 20.0    : D-dimer generated per unit fibrin lysis flux (mg/L per AU)
KDDS    : 0.060   : weight of soluble-fibrin lysis in D-dimer generation (-)
KDDD    : 0.086643 : D-dimer clearance, t1/2 8 h (1/h)
DDB     : 0.25    : D-dimer baseline (mg/L FEU)

// ---------------- organ injury ----------------
KIKID   : 0.020   : fibrin-driven renal injury (1/h per AU)
KIKIDH  : 0.0002  : histone-driven renal injury (1/h per ug/mL)
KRKID   : 0.022   : renal recovery (1/h)
KIIK    : 0.012   : inflammation-driven renal injury (1/h)
KILIV   : 0.010   : fibrin-driven hepatic injury (1/h per AU)
KRLIV   : 0.018   : hepatic recovery (1/h)
KIIL    : 0.006   : inflammation-driven hepatic injury (1/h)
KILUN   : 0.019   : fibrin-driven pulmonary injury (1/h per AU)
KRLUN   : 0.020   : pulmonary recovery (1/h)
KIIU    : 0.010   : inflammation-driven pulmonary injury (1/h)
KICNS   : 0.009   : fibrin-driven CNS injury (1/h per AU)
KRCNS   : 0.025   : CNS recovery (1/h)
KIIC    : 0.005   : inflammation-driven CNS injury (1/h)

// ---------------- bleeding ----------------
KBPLT   : 45      : platelet count at half-maximal bleeding hazard (10^9/L)
GBP     : 2.2     : platelet Hill coefficient (-)
KBFIB   : 95      : fibrinogen at half-maximal bleeding hazard (mg/dL)
GBF     : 2.6     : fibrinogen Hill coefficient (-)
KBPT    : 7.0     : PT prolongation at half-maximal bleeding hazard (s)
GBT     : 2.0     : PT Hill coefficient (-)
KBLYS   : 0.45    : plasmin at half-maximal lysis-driven bleeding (nM)
KBHEP   : 6.0     : anticoagulant intensity at half-maximal drug bleeding (-)
FHEP    : 0.50    : maximal drug-attributable bleeding hazard (-)
KBLD    : 0.10    : accumulation rate of transfusion-requiring bleeding (1/h)
KBLR    : 0.05    : resolution rate of the cumulative bleeding burden (1/h)
BIDX0   : 0.10    : bleeding-index offset (-)
HFAIL0  : 0.30    : haemostatic-deficit offset on the log scale (-)

// ---------------- outcome ----------------
HZ0     : 0.0000075 : background hazard (1/h)
HZI     : 0.00370   : inflammation/shock hazard weight, NOT reachable by anticoagulants (1/h)
HZS     : 0.0168217 : organ-failure hazard weight (1/h)
HZB     : 0.0000174 : haemostatic-failure hazard weight (1/h)
HZF     : 0.0       : optional direct microthrombus hazard, off by default (1/h)
KHF     : 1.2       : fibrin stock at half-maximal direct hazard (AU)

// ---------------- PT ----------------
PT0     : 12.0    : reference prothrombin time (s)
GPT     : 0.55    : PT sensitivity exponent to the effective factor level (-)

// ---------------- drug PK ----------------
VH      : 3.5     : UFH volume of distribution (L)
KELH    : 0.70    : UFH elimination (1/h)
KAE     : 0.50    : enoxaparin SC absorption (1/h)
VE      : 5.0     : enoxaparin volume of distribution (L)
KELE    : 0.154   : enoxaparin elimination, t1/2 4.5 h (1/h)
FE      : 0.90    : enoxaparin bioavailability (-)
RENIMP  : 0.55    : maximal renal-failure reduction in enoxaparin clearance (-)
VTM     : 3.5     : thrombomodulin alfa volume of distribution (L)
KELTM   : 0.0347  : thrombomodulin alfa elimination, t1/2 20 h (1/h)
VAPC    : 10.0    : drotrecogin alfa volume of distribution (L)
KELAPC  : 3.73    : drotrecogin alfa elimination (1/h)
VTX     : 12.0    : tranexamic acid volume of distribution (L)
KELTX   : 0.58    : tranexamic acid elimination (1/h)
VATR    : 100     : ATRA volume of distribution (L)
KAATR   : 1.5     : ATRA absorption (1/h)
KELATR  : 0.50    : ATRA elimination (1/h)
VARG    : 12.0    : argatroban volume of distribution (L)
KELARG  : 1.33    : argatroban elimination (1/h)

// ---------------- infusion rates (zero unless a scenario sets them) ----
RATE_UFH  : 0 : UFH infusion (IU/h)
RATE_ENOX : 0 : enoxaparin input into the SC depot (anti-Xa IU/h)
RATE_RTM  : 0 : thrombomodulin alfa infusion (mg/h)
RATE_APC  : 0 : drotrecogin alfa infusion (ug/h)
RATE_TXA  : 0 : tranexamic acid infusion (mg/h)
RATE_ARG  : 0 : argatroban infusion (mg/h)
RATE_ATRA : 0 : ATRA input into the gut depot (mg/h)
FIBINF    : 0 : continuous fibrinogen replacement (mg/dL/h)
PLTINF    : 0 : continuous platelet replacement (10^9/L/h)
FFPINF    : 0 : continuous FFP replacement (fraction of one full dose per h)

$INIT @annotated
PATH  :     0 : pathogen / DAMP load (AU)
BLAST :     0 : APL promyelocyte burden (10^9/L)
IL6   :   5.0 : interleukin-6 (pg/mL)
TNF   :  10.0 : tumour necrosis factor alpha (pg/mL)
NET   :  0.30 : NET-derived cell-free DNA (ug/mL)
HIST  :  0.10 : circulating histones H3/H4 (ug/mL)
TF    :  0.02 : tissue factor activity (nM-eq)
TM    :   1.0 : endothelial thrombomodulin (fraction of normal)
EPCR  :   1.0 : endothelial protein C receptor (fraction of normal)
GLX   :   1.0 : endothelial glycocalyx integrity (fraction of normal)
FII   :   1.0 : prothrombin (fraction of normal)
FV    :   1.0 : factor V (fraction of normal)
FVII  :   1.0 : factor VII (fraction of normal)
FVIII :   1.0 : factor VIII (fraction of normal)
FX    :   1.0 : factor X (fraction of normal)
AT    :   1.0 : antithrombin (fraction of normal)
PC    :   1.0 : protein C zymogen (fraction of normal)
PS    :   1.0 : free protein S (fraction of normal)
TFPI  :   1.0 : tissue factor pathway inhibitor (fraction of normal)
FIB   :   300 : fibrinogen (mg/dL)
THR   :     0 : thrombin (nM)
FXA   :     0 : factor Xa (nM)
APC   :     0 : endogenously activated protein C (nM)
PLT   :   250 : platelet count (10^9/L)
PACT  :     0 : activated platelet fraction (-)
SFM   :     0 : soluble fibrin monomer complex (nM)
FDEP  :     0 : DEPOSITED microvascular fibrin — the STOCK (AU)
TPA   :   5.0 : tissue plasminogen activator (ng/mL)
PAI1  :  20.0 : plasminogen activator inhibitor 1 (ng/mL)
PLG   :   1.0 : plasminogen (fraction of normal)
PLN   :     0 : plasmin (nM)
A2AP  :   1.0 : alpha-2-antiplasmin (fraction of normal)
TAFIA :     0 : activated TAFI (fraction)
DD    :  0.25 : D-dimer (mg/L FEU)
OKID  :     0 : renal injury index (0-1)
OLIV  :     0 : hepatic injury index (0-1)
OLUN  :     0 : pulmonary injury index (0-1)
OCNS  :     0 : cerebral injury index (0-1)
BLEEDC:     0 : cumulative transfusion-requiring bleeding burden (AU)
CUMH  :     0 : cumulative mortality hazard (-)
HEPC  :     0 : UFH central amount (IU)
ENXD  :     0 : enoxaparin SC depot (anti-Xa IU)
ENXC  :     0 : enoxaparin central amount (anti-Xa IU)
RTM   :     0 : thrombomodulin alfa amount (mg)
APCX  :     0 : drotrecogin alfa amount (ug)
TXAC  :     0 : tranexamic acid amount (mg)
ATRAD :     0 : ATRA gut depot (mg)
ATRAC :     0 : ATRA central amount (mg)
ARGA  :     0 : argatroban amount (mg)

// NOTE: initial conditions are declared HERE, not as <CMT>_0 assignments in
// $MAIN.  $MAIN runs before every simulation and would silently overwrite the
// init() values that set up the APL presenting state.  The numbers above are
// the healthy steady state and must stay consistent with IL6B / TNFB / NETB /
// HISTB / TFB / FIB0 / PLT0 / TPAB / PAIB / DDB in $PARAM.  Every pathological
// driver in $ODE acts on the INCREMENT above these values, so this state is an
// stationary state of the system.  Verified numerically in an independent
// Python/scipy re-implementation: over 28 simulated healthy days no state
// drifts by more than 1.2% (worst: platelets 250 -> 247, fibrinogen 300 -> 298).
// The residual drift is not numerical -- it is the constitutive tissue-factor
// term TFB, which sustains a basal thrombin of 1.0 pM, a basal FXa of 0.014 nM
// and a fibrin stock of 1.3e-3 AU.  That basal tone is deliberate: setting
// TFB = 0 would make the healthy state exactly stationary but would also
// remove the physiological low-grade thrombin generation that the protein C
// pathway is there to restrain.

$ODE
#define pos(x) ((x) > 0 ? (x) : 0.0)

double PATHv  = pos(PATH);
double BLASTv = pos(BLAST);

// ---------------- inflammation (increment-driven) ----------------
double ATRAE = EMAXATRA * ATRAC / VATR / (EC50ATRA + ATRAC / VATR);
dxdt_PATH  = -KCLRP * SRCCTL * PATHv * (1.0 - 0.5 * OLIV);
dxdt_BLAST =  KGB * BLASTv * (1.0 - BLASTv / BMAX) - KDIFF * ATRAE * BLASTv;

dxdt_IL6 = KDIL6 * (IL6B + KSIL6 * SIL6 * PATHv + KIL6BL * BLASTv - IL6);
dxdt_TNF = KDTNF * (TNFB + KSTNF * SIL6 * PATHv + KTNFBL * BLASTv - TNF);

double dTNF = pos(TNF - TNFB);
double dIL6 = pos(IL6 - IL6B);
double fTNF = dTNF / (KTNF + dTNF);
double APF  = dIL6 / (KIL6 + dIL6);          // <- CLOCK 1 lives in this term

dxdt_NET  = KSNET * fTNF - KDNET * (NET - NETB);
double dNET = pos(NET - NETB);
dxdt_HIST = KHNET * dNET - KDHIST * (HIST - HISTB);
double dHIST = pos(HIST - HISTB);
double fHIST = dHIST / (KHIST + dHIST);

// ---------------- drug concentrations ----------------
double HNEUT = 1.0 / (1.0 + dHIST / KHEPN);      // PF4/histone heparin neutralisation
double CHEP  = HEPC / (VH * 1000.0) * HNEUT;     // anti-Xa IU/mL
double CENX  = ENXC / (VE * 1000.0) * HNEUT;
double HEPX  = CHEP + CENX;
double CRTM  = RTM / VTM;                        // ug/mL
double CAPCX = APCX / VAPC * 0.017857;           // nM (56 kDa)
double CTXA  = TXAC / VTX;                       // ug/mL
double CARG  = ARGA / VARG;                      // ug/mL

double AMPX = 1.0 + EMAXHX * HEPX / (EC50H + HEPX);
double AMPI = 1.0 + EMAXHI * CHEP / (EC50H + CHEP);
double DTI  = 1.0 + EMAXDTI * CARG / (KDTI50 + CARG);

// ---------------- tissue factor & endothelium ----------------
dxdt_TF = KDTF * (TFB - TF) + KSTF * fTNF + KTFBL * BLASTv + TFEXO * PATHv;
double SHED = fTNF + fHIST * (1.0 - EHMGB * CRTM / (EC50TM + CRTM));
dxdt_TM   = KRTM * (1.0 - TM)   - KSHTM * SHED * TM;
dxdt_EPCR = KREP * (1.0 - EPCR) - KSHEP * SHED * EPCR;
dxdt_GLX  = KRGL * (1.0 - GLX)  - KSHGL * SHED * GLX;

// ---------------- HEPARIN IS A CATALYST ----------------------------
// Rate = KATII * (AT              <- spontaneous, LINEAR in AT
//               + (AMP-1)*ATcat)  <- heparin-catalysed, SATURATING in AT
// so AT depletion lowers the ceiling sub-linearly and AT supplementation
// above normal buys very little on the heparin arm.
double ATv   = pos(AT);
double ATcat = ATv / (KMAT + ATv) * (KMAT + 1.0);
double ATXA  = ATv + (AMPX - 1.0) * ATcat;
double ATII  = ATv + (AMPI - 1.0) * ATcat;

// ---------------- thrombin generation ----------------
double THRv  = pos(THR);
double FXAv  = pos(FXA);
double PSURF = PACT * (pos(PLT) / PLT0) + PSMV * TF / (KTFS + TF);
double APCF  = 1.0 / (1.0 + (APC + CAPCX) / KIAPC);   // FVa/FVIIIa inactivation

double vXa = KXA * TF * pos(FVII) / (1.0 + RTFPI * pos(TFPI))
           + KTEN * pos(FVIII) * APCF * THRv / (KTHRF + THRv) * PSURF;
dxdt_FXA = vXa - (KATXA * ATXA + KDXA) * FXAv;

double vIIa     = KPRO * FXAv * pos(FV) * APCF * PSURF * pos(FII);
double KILL_AT  = KATII * ATII;
double KILL_TOT = (KILL_AT + KDTHR + KTMB * TM) * DTI;
dxdt_THR = vIIa - KILL_TOT * THRv;
double ATSHARE = KILL_AT / (KILL_TOT > 1e-9 ? KILL_TOT : 1e-9);

// ---------------- the ONE consumption clock ----------------
double KCONS = KC0 * THRv / (KTHRC + THRv);
double LIVF  = 1.0 - 0.85 * OLIV;             // hepatic synthetic feedback
double LEAK  = KLEAK * (1.0 - GLX);

#define ZYM(X, TH, W, APMUL, EXTRA) \
  ( (0.693147/(TH)) * (APMUL) * LIVF - (0.693147/(TH)) * pos(X) - ((KCONS)*(W) + (EXTRA)) * pos(X) )

double APCT = APC + CAPCX;
dxdt_FII   = ZYM(FII,   THF2,  WF2,  1.0, 0.0);
dxdt_FV    = ZYM(FV,    THF5,  WF5,  1.0, 0.0) - KAPCV * APCT * pos(FV);
dxdt_FVII  = ZYM(FVII,  THF7,  WF7,  1.0, 0.0);
dxdt_FVIII = ZYM(FVIII, THF8,  WF8,  1.0 + APOS8 * APF, 0.0) - KAPC8 * APCT * pos(FVIII);
dxdt_FX    = ZYM(FX,    THF10, WF10, 1.0, 0.0);
dxdt_PC    = ZYM(PC,    THPC,  WPC,  1.0/(1.0 + ANEGPC * APF), LEAK + KELPC * fHIST);
dxdt_PS    = ZYM(PS,    THPS,  WPS,  1.0/(1.0 + ANEGPS * APF), LEAK);
dxdt_TFPI  = ZYM(TFPI,  THTFPI, WTFPI, 1.0, KELPC * fHIST);

// Antithrombin is consumed STOICHIOMETRICALLY with the inactivation flux
// it actually catalyses — which is why heparin makes AT fall further.
double kdAT = 0.693147 / THAT;
dxdt_AT = kdAT / (1.0 + ANEGAT * APF) * LIVF - kdAT * ATv
        - (LEAK + KELAT * fHIST) * ATv
        - KATC * vIIa * ATSHARE * (1.0 + KATXA * ATXA / (KATXA * ATXA + KDXA));

// ---------------- CLOCK 1 : fibrinogen ----------------
double kdF  = 0.693147 / THFIB;
double PLNv = pos(PLN);
double A2v  = pos(A2AP);
double TXAI = 1.0 / (1.0 + CTXA / IC50TXA);
double FIBv = pos(FIB);
// systemic fibrinogenolysis is gated on alpha-2-antiplasmin EXHAUSTION
double FGNLYS = KFGN * PLNv * (1.0 - A2v) * (1.0 - A2v) * TXAI;
dxdt_FIB = kdF * FIB0 * (1.0 + APOSF * APF) * LIVF - kdF * FIBv
         - KCONS * WFIB * FIBv - FGNLYS * FIBv
         + FIBINF + 0.20 * FFPINF * FIB0;

// ---------------- protein C : a product of FOUR terms ----------------
double RTMC   = EMTM * CRTM / (EC50TM + CRTM);
double TMEFF  = TM + RTMC;            if (TMEFF  > 2.0) TMEFF  = 2.0;
double TMTAFI = TM + FTAFIS * RTMC;   if (TMTAFI > 2.0) TMTAFI = 2.0;
double vAPC = KPCA * THRv * TMEFF * EPCR * pos(PC);
dxdt_APC = vAPC - KDAPC * APC;

// ---------------- platelets ----------------
double PLTv  = pos(PLT);
double PPROD = KDPLT * PLT0 * MSUP * pow(PLT0 / (PLTv > 15.0 ? PLTv : 15.0), ETPO);
double PCONS = KPC0 * THRv / (KTHRP + THRv) + KPFD * FDEP;
dxdt_PLT  = PPROD - KDPLT * PLTv - PCONS * PLTv + PLTINF;
dxdt_PACT = KACT * (WTHRP * THRv / (KTHRP + THRv) + AHISTP * fHIST) * (1.0 - PACT)
          - KDEACT * PACT;

// ---------------- CLOCK 2 : fibrinolysis ----------------
dxdt_TPA  = KDTPA * (TPAB - TPA) + KSTPA * (fTNF + AENDO * (1.0 - GLX))
          + KTPABL * BLASTv;
dxdt_PAI1 = KDPAI * (PAIB - PAI1) + KSPAI * SPAI * APF - KAPCPAI * APCT * pos(PAI1);

double TPAnM = pos(TPA)  * 0.0147;            // 68 kDa
double PAInM = pos(PAI1) * 0.0222 * PAIACT;   // 45 kDa, active fraction only
double ANX   = 1.0 + KANX * BLASTv / (KBL50 + BLASTv);
double FREETPA = TPAnM * (FLOC + (1.0 - FLOC) / (1.0 + PAInM / KPAI)) * ANX;

double FSURF = FDEP / (KFS + FDEP);
double vPLN  = (KPG * FREETPA * pos(PLG) * FSURF
             +  KUPA * BLASTv / (KBL50 + BLASTv) * pos(PLG)) * TXAI;
dxdt_PLN  = vPLN - (KA2AP * A2v + KDPLN) * PLNv;
dxdt_A2AP = KDA2 * (1.0 - A2v) - CA2 * PLNv * A2v;
dxdt_PLG  = (0.693147 / 60.0) * (1.0 - pos(PLG)) - KPLGC * vPLN;

dxdt_TAFIA = KTAFI * THRv * TMTAFI - KDTAFI * TAFIA;
double PLNEFF = PLNv / (1.0 + ATAFI * TAFIA);

// ---------------- fibrin: the RATE and the STOCK ----------------
dxdt_SFM = KFPA * THRv * FIBv * 0.0294 - KSFMD * SFM - KLYSS * PLNEFF * SFM;
double THN  = pow(THRv, NHD);
double vDEP = KDEP * THN / (pow(KTHD, NHD) + THN) * (FIBv / FIB0);
double LYSFLUX = KLYS * PLNEFF * FDEP;
double SFMLYS  = KLYSS * PLNEFF * SFM;
dxdt_FDEP = vDEP - LYSFLUX - KRES * FDEP;
dxdt_DD   = KDDG * (LYSFLUX + KDDS * SFMLYS) - KDDD * (DD - DDB);

// ---------------- organ injury ----------------
double INFL = PATHv / (KHP + PATHv);
dxdt_OKID = (KIKID * FDEP + KIIK * INFL + KIKIDH * dHIST) * (1 - OKID) - KRKID * OKID;
dxdt_OLIV = (KILIV * FDEP + KIIL * INFL) * (1 - OLIV) - KRLIV * OLIV;
dxdt_OLUN = (KILUN * FDEP + KIIU * INFL) * (1 - OLUN) - KRLUN * OLUN;
dxdt_OCNS = (KICNS * FDEP + KIIC * INFL) * (1 - OCNS) - KRCNS * OCNS;

// ---------------- bleeding: a PRODUCT of independent deficits ----------
double FIBf = (FIBv > 1.0 ? FIBv : 1.0) / 100.0;
double feff = 5.0 / (1.0/(pos(FII)  > 1e-3 ? pos(FII)  : 1e-3)
                   + 1.0/(pos(FV)   > 1e-3 ? pos(FV)   : 1e-3)
                   + 1.0/(pos(FVII) > 1e-3 ? pos(FVII) : 1e-3)
                   + 1.0/(pos(FX)   > 1e-3 ? pos(FX)   : 1e-3)
                   + 1.0/(FIBf      > 1e-3 ? FIBf      : 1e-3));
double PTs = PT0 * pow(feff / (5.0 / (4.0 + 1.0/3.0)), -GPT);

double hplt = 1.0 / (1.0 + pow(PLTv / KBPLT, GBP));
double hfib = 1.0 / (1.0 + pow(FIBv / KBFIB, GBF));
double dpt  = PTs - PT0; if (dpt < 0) dpt = 0;
double hpt  = pow(dpt / KBPT, GBT) / (1.0 + pow(dpt / KBPT, GBT));
double hlys = PLNv / (KBLYS + PLNv);
double ACG  = ATv * ((AMPI - 1.0) + (AMPX - 1.0)) + 4.0 * (DTI - 1.0);
double hhep = FHEP * ACG / (KBHEP + ACG);
double HCOMP = (1-hplt) * (1-hfib) * (1-hpt) * (1-hlys) * (1-hhep);
double BIDXv = 1.0 - HCOMP;
double HFAIL = -log(HCOMP > 1e-6 ? HCOMP : 1e-6);
dxdt_BLEEDC = KBLD * BIDXv - KBLR * BLEEDC;

// ---------------- mortality hazard ----------------
// Note the FIRST term: the part of septic mortality that no anticoagulant
// can reach.  It is the structural reason every anticoagulant trial in
// sepsis has produced a small absolute effect at best.
double ORG = OKID + OLIV + OLUN + OCNS;
double hf  = HFAIL - HFAIL0; if (hf < 0) hf = 0;
dxdt_CUMH = HZ0 + HZI * INFL * INFL + HZS * pow(ORG/4.0, 2.0)
          + HZB * hf * hf + HZF * FDEP / (KHF + FDEP);

// ---------------- drug PK ----------------
dxdt_HEPC  = RATE_UFH  - KELH * HEPC;
dxdt_ENXD  = RATE_ENOX - KAE * ENXD;
dxdt_ENXC  = FE * KAE * ENXD - KELE * (1.0 - RENIMP * OKID) * ENXC;
dxdt_RTM   = RATE_RTM  - KELTM * RTM;
dxdt_APCX  = RATE_APC  - KELAPC * APCX;
dxdt_TXAC  = RATE_TXA  - KELTX * TXAC;
dxdt_ATRAD = RATE_ATRA - KAATR * ATRAD;
dxdt_ATRAC = KAATR * ATRAD - KELATR * ATRAC;
dxdt_ARGA  = RATE_ARG  - KELARG * ARGA;

$TABLE
// ---- derived laboratory values ----
double FIBo  = FIB > 1.0 ? FIB : 1.0;
double FIBfo = FIBo / 100.0;
double feffo = 5.0 / (1.0/(FII  > 1e-3 ? FII  : 1e-3) + 1.0/(FV   > 1e-3 ? FV   : 1e-3)
                    + 1.0/(FVII > 1e-3 ? FVII : 1e-3) + 1.0/(FX   > 1e-3 ? FX   : 1e-3)
                    + 1.0/(FIBfo> 1e-3 ? FIBfo: 1e-3));
capture PT  = PT0 * pow(feffo / (5.0 / (4.0 + 1.0/3.0)), -GPT);
capture INR = PT / PT0;
capture PLTc = PLT;
capture FIBc = FIB;
capture DDc  = DD;
capture ATpc = 100.0 * AT;
capture PCpc = 100.0 * PC;
capture A2pc = 100.0 * A2AP;
capture TMpc = 100.0 * TM;

// ---- SOFA (cardiovascular component not modelled) ----
double coag = PLT < 20 ? 4 : (PLT < 50 ? 3 : (PLT < 100 ? 2 : (PLT < 150 ? 1 : 0)));
capture SOFA = 4.0 * (OLUN + OKID + OLIV + OCNS) + coag;
capture ORGF = OKID + OLIV + OLUN + OCNS;

// ---- diagnostic scores are OUTPUTS of the state vector, not inputs ----
double sPLT = PLT < 50 ? 2 : (PLT < 100 ? 1 : 0);
double sDD  = DD  > 5  ? 3 : (DD  > 1   ? 2 : 0);
double sPT  = (PT - 12.0) > 6 ? 2 : ((PT - 12.0) > 3 ? 1 : 0);
double sFIB = FIB < 100 ? 1 : 0;
capture ISTH = sPLT + sDD + sPT + sFIB;             // >= 5 : overt DIC
double cPLT = PLT < 100 ? 2 : (PLT < 150 ? 1 : 0);
double cINR = INR > 1.4 ? 2 : (INR > 1.2 ? 1 : 0);
double cSOF = SOFA >= 2 ? 2 : (SOFA >= 1 ? 1 : 0);
capture SIC  = cPLT + cINR + cSOF;                  // >= 4 : SIC
capture OVERT = ISTH >= 5 ? 1 : 0;

// ---- bleeding index ----
double hpltT = 1.0 / (1.0 + pow(PLT / KBPLT, GBP));
double hfibT = 1.0 / (1.0 + pow(FIBo / KBFIB, GBF));
double dptT  = PT - PT0; if (dptT < 0) dptT = 0;
double hptT  = pow(dptT / KBPT, GBT) / (1.0 + pow(dptT / KBPT, GBT));
double hlysT = PLN / (KBLYS + PLN);
capture BIDX = 1.0 - (1-hpltT) * (1-hfibT) * (1-hptT) * (1-hlysT);
capture MORT = 100.0 * (1.0 - exp(-CUMH));

// ---- the two clocks, reported explicitly ----
double dIL6T = IL6 - IL6B; if (dIL6T < 0) dIL6T = 0;
capture APFo    = dIL6T / (KIL6 + dIL6T);                 // CLOCK 1
capture FIBSYN  = 1.0 + APOSF * APFo;                     // acute-phase multiplier
double TPAnMT   = TPA  * 0.0147;
double PAInMT   = PAI1 * 0.0222 * PAIACT;
capture PAITPA  = (TPA > 0) ? (PAI1 * 0.0222) / (TPA * 0.0147) : 0;   // CLOCK 2, molar
capture FREETPA = TPAnMT * (FLOC + (1.0 - FLOC) / (1.0 + PAInMT / KPAI))
                * (1.0 + KANX * BLAST / (KBL50 + BLAST));
capture APCPROD = TM * EPCR * PC;        // the three fractions that multiply
capture APCTHR  = (THR > 1e-9) ? APC / THR : 0;   // APC generated per unit thrombin
'

mod <- mcode_cache("dic_qsp", dic_code, atol = 1e-9, rtol = 1e-7, maxsteps = 1e6)

WT <- 70   # reference adult body weight (kg)

## =====================================================================
##  AETIOLOGY PRESETS
##  The two clocks are the ONLY things that differ between the phenotypes.
##  Everything else -- every rate constant, every half-life, every
##  stoichiometric weight -- is shared.
## =====================================================================
AETIOLOGY <- list(

  healthy = list(
    par  = list(),
    init = list()),

  # Septic DIC.  CLOCK 1 forward (IL-6 storm -> fibrinogen synthesis x4),
  # CLOCK 2 shut down (PAI-1 induced ~28-fold over t-PA).
  sepsis = list(
    par  = list(SIL6 = 1.00, SPAI = 1.00, MSUP = 1.00),
    init = list(PATH = 1.0)),

  # Identical infection in a host whose monocytes barely induce tissue
  # factor.  Used to isolate the mortality DIC itself contributes, with
  # infection severity held constant.
  sepsis_noDIC = list(
    par  = list(SIL6 = 1.00, SPAI = 1.00, KSTF = 0.06),
    init = list(PATH = 1.0)),

  # APL at diagnosis.  CLOCK 1 nearly off (no IL-6 storm, so nothing
  # replenishes fibrinogen) and CLOCK 2 reversed (PAI-1 not induced,
  # annexin A2 on the blast surface amplifies t-PA).
  apl = list(
    par  = list(SIL6 = 0.10, SPAI = 0.15, MSUP = 0.20),
    init = list(BLAST = 30, PLT = 35, FIB = 110,
                A2AP = 0.40, PLG = 0.60, DD = 22)),

  # Obstetric catastrophe / major trauma: a large one-off tissue-factor
  # load, intermediate acute-phase response, brisk early fibrinolysis.
  obstetric = list(
    par  = list(SIL6 = 0.45, SPAI = 0.40, TFEXO = 2.2, KCLRP = 0.12),
    init = list(PATH = 1.0))
)

## ---------------------------------------------------------------------
##  Dosing helpers
## ---------------------------------------------------------------------

## Combine event objects safely (c() on ev objects does NOT concatenate).
ev_c <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(NULL)
  as.ev(dplyr::bind_rows(lapply(parts, as.data.frame)))
}

bolus <- function(cmt, amt, times) ev(time = times, amt = amt, cmt = cmt)

## An infusion is a piecewise-constant change to a RATE_* parameter, so it
## can be started and stopped at arbitrary times.
infusion <- function(param, rate, start, stop) {
  data.frame(time = c(start, stop), name = param, value = c(rate, 0),
             stringsAsFactors = FALSE)
}
inf_c <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(NULL)
  dplyr::bind_rows(parts)
}

## Antithrombin concentrate: 1 IU/kg raises AT activity by about 1.4%, so
## the dose enters the AT compartment directly as a fraction of normal.
ATC_LOAD  <- 0.014 * (6000 / WT)    # 6000 IU load  -> +1.20  (AT ~ 220%)
ATC_MAINT <- 0.014 * (3000 / WT)    # 3000 IU q12h  -> +0.60 each

## ---------------------------------------------------------------------
##  Simulation driver
## ---------------------------------------------------------------------
run_dic <- function(aet = "sepsis", events = NULL, infusions = NULL,
                    par = list(), end = 672, delta = 0.5) {

  A <- AETIOLOGY[[aet]]
  if (is.null(A)) stop("unknown aetiology: ", aet)

  m <- mod %>% param(modifyList(A$par, par))
  if (length(A$init)) m <- m %>% init(A$init)

  if (is.null(infusions) || !nrow(infusions)) {
    return(as.data.frame(mrgsim(m, events = events, end = end, delta = delta)))
  }

  ## Integrate segment by segment so infusion rates can switch mid-run.
  cmts  <- as.character(mrgsolve::cmt(mod))
  brk   <- sort(unique(c(0, infusions$time[infusions$time < end], end)))
  edf   <- if (is.null(events)) NULL else as.data.frame(events)
  cur   <- m
  res   <- NULL
  carry <- NULL

  for (i in seq_len(length(brk) - 1L)) {
    t0 <- brk[i]; t1 <- brk[i + 1L]

    act <- infusions[infusions$time <= t0, , drop = FALSE]
    if (nrow(act)) {
      last <- act[!duplicated(act$name, fromLast = TRUE), , drop = FALSE]
      cur  <- param(cur, setNames(as.list(last$value), last$name))
    }
    if (!is.null(carry)) cur <- init(cur, carry)

    e <- NULL
    if (!is.null(edf)) {
      sel <- edf[edf$time >= t0 & edf$time < t1, , drop = FALSE]
      if (nrow(sel)) { sel$time <- sel$time - t0; e <- as.ev(sel) }
    }

    seg <- mrgsim(cur, events = e, end = t1 - t0, delta = delta) %>% as.data.frame()
    seg$time <- seg$time + t0
    carry <- as.list(seg[nrow(seg), cmts])
    res   <- rbind(res, if (i == 1L) seg else seg[-1L, , drop = FALSE])
  }
  res
}

## =====================================================================
##  SIXTEEN THERAPEUTIC SCENARIOS
## =====================================================================
SCEN <- list(

  # ---- reference states ------------------------------------------------
  A_healthy = list(aet = "healthy",
    desc = "Healthy control -- must be an exact fixed point of the system"),

  B_sepsis_noDIC = list(aet = "sepsis_noDIC",
    desc = "Same infection, DIC-resistant host -- isolates DIC-attributable death"),

  C_septic_DIC = list(aet = "sepsis",
    desc = "Septic DIC, supportive care only -- the reference arm"),

  # ---- heparins: an effect that is a PRODUCT with a saturable AT term ---
  D_UFH = list(aet = "sepsis",
    inf = infusion("RATE_UFH", 18 * WT, 0, 336),
    desc = "UFH 18 U/kg/h IV for 14 d"),

  E_enoxaparin = list(aet = "sepsis",
    ev = bolus("ENXD", 1 * WT * 100, seq(0, 335, 12)),
    desc = "Enoxaparin 1 mg/kg SC q12h -- renally cleared, accumulates as AKI develops"),

  M_argatroban = list(aet = "sepsis",
    inf = infusion("RATE_ARG", 8.4, 0, 336),
    desc = "Argatroban 2 ug/kg/min -- a DIRECT thrombin inhibitor, with no AT term"),

  # ---- protein therapeutics: each enters a DIFFERENT factor -------------
  F_AT_alone = list(aet = "sepsis",
    ev = ev_c(bolus("AT", ATC_LOAD, 0), bolus("AT", ATC_MAINT, seq(12, 84, 12))),
    desc = "Antithrombin concentrate WITHOUT heparin (KyberSept no-heparin subgroup)"),

  G_AT_plus_UFH = list(aet = "sepsis",
    ev  = ev_c(bolus("AT", ATC_LOAD, 0), bolus("AT", ATC_MAINT, seq(12, 84, 12))),
    inf = infusion("RATE_UFH", 18 * WT, 0, 336),
    desc = "Antithrombin concentrate WITH concomitant heparin (KyberSept as conducted)"),

  H_thrombomodulin = list(aet = "sepsis",
    ev = bolus("RTM", 0.06 * WT, seq(0, 120, 24)),
    desc = "Thrombomodulin alfa (ART-123) 0.06 mg/kg/d x 6 d (SCARLET)"),

  I_rhAPC_early = list(aet = "sepsis",
    inf = infusion("RATE_APC", 24 * WT, 6, 102),
    desc = "Drotrecogin alfa 24 ug/kg/h x 96 h, started at 6 h (PROWESS-like)"),

  J_rhAPC_late = list(aet = "sepsis",
    inf = infusion("RATE_APC", 24 * WT, 30, 126),
    desc = "Drotrecogin alfa started at 30 h -- does TIMING explain PROWESS-SHOCK?"),

  # ---- antifibrinolytics: ONE term whose clinical sign flips ------------
  L_TXA_sepsis = list(aet = "sepsis",
    ev  = bolus("TXAC", 1000, 0),
    inf = infusion("RATE_TXA", 125, 0, 96),
    desc = "Tranexamic acid in the fibrinolytic-SHUTDOWN phenotype -- predicted HARM"),

  # ---- APL: the other setting of the same two clocks --------------------
  N_APL_supportive = list(aet = "apl",
    desc = "APL DIC, supportive care only -- the haemorrhagic phenotype"),

  O_APL_ATRA = list(aet = "apl",
    ev = bolus("ATRAD", 40, seq(0, 660, 12)),
    desc = "ATRA 45 mg/m2/d -- removes the cause of BOTH clock settings"),

  Q_APL_ATRA_TXA = list(aet = "apl",
    ev  = ev_c(bolus("ATRAD", 40, seq(0, 660, 12)),
               bolus("FIB", 100, seq(0, 72, 12)),
               bolus("TXAC", 1000, 0)),
    inf = infusion("RATE_TXA", 125, 0, 96),
    desc = "ATRA + fibrinogen + TXA -- reproduces the ATRA/antifibrinolytic thrombosis risk"),

  R_APL_ATRA_platelets = list(aet = "apl",
    ev = ev_c(bolus("ATRAD", 40, seq(0, 660, 12)), bolus("PLT", 40, seq(0, 144, 24))),
    desc = "ATRA + platelets -- correct the WORST term of the product first")
)

run_scenario <- function(key, end = 672) {
  S   <- SCEN[[key]]
  out <- run_dic(aet = S$aet, events = S$ev, infusions = S$inf,
                 par = if (is.null(S$par)) list() else S$par, end = end)
  out$scenario <- key
  out
}

run_all <- function(end = 672) dplyr::bind_rows(lapply(names(SCEN), run_scenario, end = end))

summarise_dic <- function(df, hours = c(24, 48, 96, 672)) {
  df %>%
    dplyr::filter(time %in% hours) %>%
    dplyr::transmute(scenario, time,
      PLT = round(PLTc), FIB = round(FIBc), INR = round(INR, 2),
      DD = round(DDc, 1), AT = round(ATpc), PC = round(PCpc),
      THR = round(THR, 2), PLN = round(PLN, 2), A2AP = round(A2pc),
      FDEP = round(FDEP, 2), SOFA = round(SOFA, 1), ISTH, SIC,
      BIDX = round(BIDX, 2), MORT = round(MORT, 1))
}

## Start-time sweep used to test whether "given too late" explains
## PROWESS-SHOCK.  (In this model it does not -- see CALIBRATION NOTES.)
apc_timing_sweep <- function(starts = c(0, 6, 12, 24, 30, 48, 72, 96)) {
  do.call(rbind, lapply(starts, function(t0) {
    r <- run_dic("sepsis", infusions = infusion("RATE_APC", 24 * WT, t0, t0 + 96))
    data.frame(start_h = t0, MORT28 = round(tail(r$MORT, 1), 1))
  }))
}

## ---------------------------------------------------------------------
## Example session:
##   res <- run_all()
##   summarise_dic(res) %>% dplyr::filter(time == 48)  %>% print(n = 40)
##   summarise_dic(res) %>% dplyr::filter(time == 672) %>%
##     dplyr::select(scenario, MORT)
##   apc_timing_sweep()
## ---------------------------------------------------------------------
## =====================================================================
##  CALIBRATION NOTES — what this parameter set reproduces
##  (values below are the ACHIEVED simulated results, not targets)
## ---------------------------------------------------------------------
##  Septic DIC, supportive care, 48 h:
##      platelets 78 x10^9/L · fibrinogen 245 mg/dL (NORMAL — clock 1)
##      INR 1.39 · D-dimer 10.3 mg/L · AT 60% · protein C 33%
##      thrombomodulin 32% · ISTH 5 (overt) · SIC 5 · SOFA 8.3
##  At 24 h the same run scores ISTH 2 (NOT overt) while SIC is already 4.
##      -> the model reproduces, from the equations alone, the observation
##         that motivated the SIC score: the ISTH score gives fibrinogen a
##         point only below 100 mg/dL, and in sepsis fibrinogen never gets
##         there because IL-6 triples its synthesis.
##
##  CLOCK 1, same consumption rate, only the acute-phase sign differs:
##      sepsis  (APOSF 3, APF 0.91) -> fibrinogen steady state 281 mg/dL
##      APL     (APOSF 3, APF 0.38) -> fibrinogen steady state  87 mg/dL
##
##  CLOCK 2, simulated at 48 h:
##      healthy     tPA   5.0  PAI-1  20 ng/mL   molar  6.0:1  free tPA 0.0336 nM
##      septic DIC  tPA  20.2  PAI-1 564 ng/mL   molar 42.1:1  free tPA 0.0222 nM
##      APL         tPA  20.9  PAI-1  53 ng/mL   molar  3.9:1  free tPA 0.2882 nM
##      -> in sepsis tPA rose 4-fold and free tPA FELL to 0.66x of healthy.
##      -> APL free tPA is 13x the septic value.
##
##  The four-term APC product at 48 h in septic DIC:
##      TM 0.318 x EPCR 0.559 x PC 0.333 = 0.059 of normal
##      APC generated per unit thrombin falls 17-fold vs healthy.
##
##  Heparin's ceiling (KATII*(AT + (AMP-1)*ATcat), 1/h, UFH 0.5 IU/mL):
##      AT 100% -> 17.7    AT 60% -> 15.0    AT 40% -> 12.8    AT 220% -> 22.3
##      at AT 40% the ceiling is 72% of the AT-100% ceiling;
##      raising AT to 220% buys +26% on the heparin arm but 3.7x on the
##      UNCATALYSED arm — which is why AT concentrate helps heparin-naive
##      patients more (simulated: -7.3 mortality points without heparin,
##      -5.3 with, and a higher bleeding burden when combined).
##
##  28-day mortality, simulated:
##      B DIC-resistant  28.4 %      C septic DIC       43.3 %
##      D UFH            28.3 %      E enoxaparin       36.3 %
##      F AT alone       36.0 %      G AT + UFH         23.0 %
##      H rTM            43.8 %      I rhAPC early      37.4 %
##      J rhAPC late     35.8 %      L TXA              51.1 %
##      M argatroban     32.2 %      N APL supportive   32.8 %
##      O APL + ATRA      7.5 %      Q APL+ATRA+TXA     21.8 %
##      R APL+ATRA+plts   6.7 %
##
##  Agreement and disagreement with the trials (both stated):
##    AGREES  rhAPC early -5.9 points vs PROWESS -6.1
##    AGREES  AT concentrate without heparin -7.3 vs KyberSept subgroup -5.8
##    AGREES  thrombomodulin alfa essentially neutral, as in SCARLET
##    AGREES  ATRA cuts APL early death by ~three quarters
##    AGREES  antifibrinolytics on top of ATRA are harmful in APL
##    DISAGREES  the model makes UFH (-15 points) and AT+UFH (-20 points)
##               far more effective than any trial or meta-analysis has
##               shown.  This is the model's most exposed prediction and
##               should be treated as a falsifiable hypothesis, not a
##               result.  The likeliest missing mechanisms are heparin-
##               induced bleeding into already-injured organs and
##               practical under-dosing from monitoring difficulty in DIC.
##    REFUTED A PRIOR HYPOTHESIS  the "a rate-acting drug given after the
##               stock exists cannot reverse it" argument does NOT explain
##               PROWESS-SHOCK in this model.  A start-time sweep gives
##               37.8 / 37.4 / 36.9 / 36.1 / 35.8 / 35.9 / 37.3 / 39.1 %
##               for starts at 0 / 6 / 12 / 24 / 30 / 48 / 72 / 96 h:
##               benefit PEAKS at 24-48 h, tracking peak fibrin deposition,
##               and only collapses beyond ~72 h.  Timing matters, but not
##               in the direction the "too late" story assumes.
## =====================================================================
