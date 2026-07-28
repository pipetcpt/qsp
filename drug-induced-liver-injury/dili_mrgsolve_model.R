# =============================================================================
# Drug-Induced Liver Injury (DILI) — QSP model, mrgsolve specification
# =============================================================================
#
# THESIS ENCODED IN THIS MODEL
# ----------------------------
#   DILI is a RATE problem, not a DOSE problem.
#
#   Every toxic species below is governed by a balance between a production
#   flux and a disposal flux, and every disposal flux saturates:
#
#       reactive metabolite   vs   glutathione RESYNTHESIS (cysteine-limited)
#       mitochondrial damage  vs   mitochondrial REPAIR
#       bile acids            vs   BSEP EXPORT (itself ATP-dependent)
#       hepatocyte loss       vs   REGENERATION
#
#   And one loop is POSITIVE: p-JNK docks on mitochondrial Sab, inhibits
#   respiration, makes more ROS, which activates more JNK. With a saturating
#   gain this loop is BISTABLE. Consequences that the model therefore
#   *computes* rather than contains as parameters:
#
#     * a sharp dose threshold (there is no "threshold" parameter anywhere)
#     * an antidote window that MOVES with dose and with host state
#     * the hepatocellular / cholestatic R-ratio (no "pattern" parameter)
#     * Hy's Law as a conjunction of a RATE arm (ALT) and a RESERVE arm (TBIL)
#
# STRUCTURE:  33 ODE compartments
#   PK (4)          AGUT AC AP ALIV
#   cofactors (2)   PAPS UDPGA
#   toxicodyn (9)   RM GSH CYS ADD ROS NRF2 MITO ATP JNK
#   bile (2)        BAH BAP
#   tissue (3)      HEP HEPS NECR
#   innate (4)      DAMP KC TNF IL10
#   adaptive (2)    TCELL TREG
#   biomarkers (6)  ALT AST TBIL ALP FV MIR
#   antidote (1)    NACC
#
# UNITS: time = hours; drug amounts = umol; hepatic concentrations = uM;
#        GSH = mM; cell masses = fraction of normal liver; ALT/AST/ALP = U/L;
#        TBIL = mg/dL; bile acids = uM.
#
# CALIBRATION ANCHORS (see dili_references.md for the numbered sources):
#   - APAP clearance ~19-21 L/h with glucuronide:sulfate:CYP ~ 55:30:6  [19]
#   - sulfation saturates first (PAPS-limited) at supratherapeutic dose [20]
#   - hepatic GSH 5-10 mM, cysteine is the rate-limiting substrate      [25,26]
#   - ATP present -> apoptosis; ATP absent -> oncotic necrosis          [31]
#   - JNK-Sab loop is required for injury (Sab-null mice are protected) [36,37]
#   - NAC within ~8-10 h nearly abolishes hepatotoxicity, benefit falls
#     away by 16-24 h                                                   [79]
#   - BSEP IC50 < ~25-50 uM marks cholestatic DILI liability            [48,49]
#   - ALT serum t1/2 ~47 h, AST ~17 h, miR-122 much shorter             [71,77]
#   - Hy's Law: ALT >= 3xULN AND TBIL >= 2xULN -> ~10% death/transplant [10]
#
# VERIFICATION: dili_reference_check.py re-implements this exact system in
# Python/scipy and regenerates every number quoted in README.md. If you edit
# the equations here, edit that file too and re-run it.
#
# Author: QSP-routine 2026-07-28 | requires mrgsolve >= 1.0
# =============================================================================

library(mrgsolve)

dili_code <- '
$PROB
# Drug-Induced Liver Injury (DILI) QSP model — 33 ODEs
# Reactive-metabolite flux vs thiol resynthesis; mitochondrial bioenergetics;
# bistable JNK-Sab amplification; BSEP-mediated cholestasis; DAMP-driven innate
# immunity; HLA-restricted adaptive immunity; regeneration reserve.
# Time unit: hours.

$PARAM @annotated
// ------------------------------------------------ physiology / anatomy ----
WT      :   70.0 : Body weight (kg)
VLIV    :    1.5 : Liver volume (L)
QH      :   80.0 : Hepatic blood flow (L/h) [~1.35 L/min]
VC      :   63.0 : Central volume of distribution (L) [APAP ~0.9 L/kg]
VP      :   25.0 : Peripheral volume (L)
QD      :   25.0 : Intercompartmental clearance (L/h)
CLR     :    1.4 : Renal clearance of parent drug (L/h)
KP      :    1.0 : Liver-to-plasma partition coefficient (OATP uptake raises it)
FULIV   :    0.8 : Unbound fraction in liver
FA      :   0.88 : Oral bioavailability x fraction absorbed
KA      :    1.2 : First-order absorption rate (1/h)

// ---------------------------------- phase II conjugation (cofactor-limited)
VMAX_SULT: 1800.0 : Vmax sulfation (umol/h)
KM_SULT  :  300.0 : Km sulfation (uM) — saturates first
VMAX_UGT :33000.0 : Vmax glucuronidation (umol/h)
KM_UGT   : 3000.0 : Km glucuronidation (uM)
APAPS    :  400.0 : PAPS cofactor pool (umol)
KSP      :    3.0 : PAPS regeneration rate (1/h) -> max sulfation 1200 umol/h
AUD      : 1500.0 : UDPGA cofactor pool (umol)
KSU      :    4.0 : UDPGA regeneration rate (1/h) -> max glucuronidation 6000

// --------------------------------------------------------- bioactivation --
VMAX_CYP : 6000.0 : Vmax CYP bioactivation (umol/h)
KM_CYP   : 6000.0 : Km CYP bioactivation (uM) — near-linear, so it does NOT
                  : saturate while phase II does; the shunt to reactive
                  : metabolite therefore rises with dose ON ITS OWN
FCYP     :    1.0 : CYP2E1 induction multiplier (2.2 = chronic ethanol)
FGSH     :    1.0 : Multiplier on GSH pool & synthesis (0.6 = fasting/ETOH)

// ------------------------------------------------- reactive metabolite ----
KGST     :   77.0 : GSH conjugation rate constant (per mM per h)
KBIND    :    3.0 : Covalent protein-binding rate constant (1/h)
KNQO     :    1.5 : Reductive/alternative detoxification (1/h)
FMITO    :    0.3 : Fraction of adducts formed on mitochondrial proteins

// -------------------------------------------------------- glutathione -----
GSH0     :    6.5 : Baseline hepatic GSH (mM)
VMAX_GSH :   0.65 : Max GSH synthesis rate constant (mM/h)
KCYS     :    1.0 : Km of GCL for cysteine (normalised units)
KDEGG    :   0.05 : Basal net GSH loss (1/h)
KOX      :   0.05 : GSH oxidation by excess ROS (1/h)
KCYSIN   :    0.5 : Cysteine pool replenishment rate (1/h)
CYSBASE  :   1.65 : Cysteine supply set-point (1.05 = malnourished)
ALPHA_CYS:    1.0 : Cysteine consumed per unit GSH synthesised

// -------------------------------------------------------------- adducts ---
KADD_REP :   0.06 : Adduct removal (autophagy/proteolysis) (1/h)

// --------------------------------------------------- ROS and Nrf2 ---------
ROS_BASE :    5.0 : Basal ROS production (au/h)
KROS_EL  :    5.0 : ROS elimination rate constant (1/h)
KROS_MITO:    2.5 : ROS from mitochondrial dysfunction (au/h)
KROS_SAB :   12.0 : ROS from JNK-Sab docking — THE POSITIVE FEEDBACK GAIN
KROS_KC  :    3.0 : ROS from activated Kupffer cells / neutrophils
KROS_BA  :    2.5 : ROS from excess intracellular bile acids
KROS_NAC :  0.004 : Direct ROS scavenging by NAC (per uM per h)
KNRF     :   0.15 : Nrf2 response rate (1/h)
ENRF     :    1.8 : Max fold-induction of Nrf2 above baseline
KNRF_H   :    0.8 : ROS excess giving half-maximal Nrf2 induction

// ---------------------------------------------- mitochondria and ATP ------
KMD_ADD  :  0.016 : Mito damage per unit adduct (1/h)
KMD_ROS  :   0.05 : Mito damage per unit excess ROS (1/h)
KMD_JNK  :   0.18 : Mito damage per unit active JNK (1/h)
KMD_BA   :   0.05 : Mito damage per unit excess bile acid (1/h)
KMR      :   0.10 : Mitochondrial repair/biogenesis rate (1/h)
KATP     :    2.0 : ATP equilibration rate (1/h)
FUNC     :    0.0 : Uncoupler / FAO-block burden (0-1); e.g. valproate

// ------------------------------------------- JNK-Sab positive feedback ----
KJ_ON    :   0.50 : Max JNK activation rate (1/h)
KJ_OFF   :   0.12 : JNK dephosphorylation rate (1/h)
KJ_HILL  :    2.0 : Hill coefficient of ROS -> JNK activation
KJ_K     :    3.0 : Excess ROS giving half-maximal JNK activation
KJ_TNF   :   0.10 : TNF-driven JNK activation (1/h)
KTNF_J   :    0.4 : TNF giving half-maximal JNK activation

// ------------------------------------------------------- bile acids -------
BAH0     :   30.0 : Baseline intrahepatocyte bile acids (uM)
BAP0     :    3.0 : Baseline serum total bile acids (uM)
BACRIT   :  100.0 : Intracellular bile-acid cytotoxicity threshold (uM)
VMAX_BSEP: 1064.0 : Vmax of BSEP canalicular export (umol/h)
KM_BSEP  :   60.0 : Km of BSEP (uM)
KI_BSEP  :  1e9   : Drug Ki for BSEP inhibition (uM); 1e9 = no liability
KMRP     :   0.60 : Basolateral MRP3/4 efflux rate constant (1/h)
EFXR_MRP :    2.0 : FXR-driven induction of MRP3/4
KFXR     :   80.0 : Bile acid giving half-maximal FXR activation (uM)
FMAX_FXR :   0.60 : Max fractional suppression of CYP7A1 by FXR (incomplete!)
KSYN_BA  :  420.0 : Max de novo bile acid synthesis (umol/h)
KUPT     :    3.0 : NTCP re-uptake rate constant (L/h)
KEL_BA   :   0.15 : Systemic bile acid elimination (1/h)
FCONV_BA : 0.0134 : Hepatic efflux -> plasma conversion factor

// -------------------------------------------------- hepatocyte life cycle -
KS_ADD   :  0.012 : Stress transition per unit adduct (1/h)
KS_BA    :  0.010 : Stress transition per unit excess bile acid (1/h)
KS_TNF   :  0.010 : Stress transition per unit TNF (1/h)
KS_KILL  :  0.030 : Stress transition per unit effector T cell (1/h)
KREC     :   0.10 : Recovery of stressed hepatocytes (1/h)
KNEC     :   0.10 : Max necrosis rate of stressed hepatocytes (1/h)
KJ_MPT   :   0.40 : JNK giving half-maximal MPT pore opening
KATP_MPT :   0.45 : ATP at which the MPT ATP-gate is half open
KBA_DEATH:  0.015 : Direct bile-acid-mediated killing
KTC_DEATH:  0.025 : Direct T-cell-mediated killing
KREG     :  0.045 : Regeneration rate constant (1/h)
AGEF     :    1.0 : Age factor on regeneration (0.6 = elderly)

// ------------------------------------------------------ innate immunity ---
KDAMP    :   30.0 : DAMP release per unit necrosis
KDAMP_EL :   0.30 : DAMP clearance (1/h)
KKC      :   0.50 : Kupffer activation rate (1/h)
KKC_OFF  :   0.10 : Kupffer deactivation rate (1/h)
KIL10    :   0.05 : IL-10 potency on Kupffer deactivation
KTNF     :   0.30 : TNF production by Kupffer cells (1/h)
KTNF_EL  :   0.60 : TNF elimination (1/h)
KIL10P   :   0.06 : IL-10 production by Kupffer cells (1/h)
KIL10_EL :   0.10 : IL-10 elimination (1/h)
DANGER0  :  0.004 : Baseline danger signal (low-grade hepatocyte turnover)

// ----------------------------------------------------- adaptive immunity --
HLA      :    0.0 : Risk-allele carriage (0 = none, 1 = e.g. HLA-B*57:01)
KT_ON    :  0.006 : Drug-specific T-cell priming rate
KT_OFF   :  0.010 : Effector T-cell contraction rate (1/h)
KD_DAMP  :   0.05 : DAMP giving half-maximal co-stimulation
KTREG_S  :   0.25 : Treg potency in suppressing priming
KTREG    :   0.02 : Treg pool turnover rate (1/h)
ICI      :    0.0 : Immune-checkpoint blockade (0 = intact, 1 = tolerance off)
KTREG_IL10:  0.05 : IL-10-driven Treg induction
STER     :    0.0 : Corticosteroid given (0/1)
STER_EFF :   0.85 : Max fractional suppression of T-cell effector function
STER_T0  :  1e9   : Time (h) at which corticosteroid starts
STER_TAU :    6.0 : Onset time constant of the steroid effect (h). A smooth
                  : onset is both more physiological than a step and far
                  : kinder to the ODE solver.

// ------------------------------------------------------------ biomarkers --
ALT0     :   25.0 : Baseline serum ALT (U/L)
ULN_ALT  :   40.0 : Upper limit of normal, ALT (U/L)
KALT_REL :12000.0 : ALT released per unit fractional necrosis
KALT_EL  :0.01475 : ALT elimination (1/h) [t1/2 = 47 h]
AST0     :   25.0 : Baseline serum AST (U/L)
KAST_REL :16800.0 : AST released per unit fractional necrosis
KAST_EL  :0.04077 : AST elimination (1/h) [t1/2 = 17 h]
ALP0     :   70.0 : Baseline serum ALP (U/L)
ULN_ALP  :  120.0 : Upper limit of normal, ALP (U/L)
KALP_REL :    9.0 : ALP released per unit cholangiocyte injury
KALP_EL  :0.00413 : ALP elimination (1/h) [t1/2 = 7 days]
KCHOL    :   0.35 : Cholangiocyte injury per unit excess bile acid
TBIL0    :    0.6 : Baseline total bilirubin (mg/dL)
ULN_TBIL :    1.2 : Upper limit of normal, total bilirubin (mg/dL)
KBIL_CLR :   0.15 : Bilirubin clearance rate constant (1/h) at full liver mass
KBIL_BA  :  300.0 : Bile acid causing half-maximal MRP2 impairment (uM)
KBIL_TNF :    0.8 : TNF causing half-maximal inflammatory cholestasis
KBIL_ALT :  0.003 : Non-hepatic (renal/alternative) bilirubin elimination (1/h)
                  : Without this term TBIL diverges as hepatocyte mass -> 0.
KFV_EL   : 0.1155 : Clotting factor turnover (1/h) [factor VII t1/2 ~6 h]
MIR0     :    1.0 : Baseline miR-122 (fold)
KMIR_REL : 9000.0 : miR-122 released per unit fractional necrosis
KMIR_EL  :   0.23 : miR-122 elimination (1/h) [t1/2 ~3 h]

// ------------------------------------------------------------------- NAC --
VNAC     :   33.0 : NAC volume of distribution (L)
CLNAC    :   11.0 : NAC clearance (L/h)
EMAX_NAC :    1.6 : Max fold-increase in cysteine supply from NAC
EC50_NAC :  150.0 : NAC concentration for half-maximal cysteine supply (uM)

$CMT @annotated
AGUT  : Drug in gut lumen (umol)
AC    : Drug in central compartment (umol)
AP    : Drug in peripheral compartment (umol)
ALIV  : Drug in liver (umol)
PAPS  : PAPS cofactor availability (fraction of pool)
UDPGA : UDPGA cofactor availability (fraction of pool)
RM    : Reactive metabolite in liver (uM)
GSH   : Hepatic reduced glutathione (mM)
CYS   : Cysteine availability (normalised, 1 = rest)
ADD   : Mitochondrial protein adducts (arbitrary units)
ROS   : Reactive oxygen species (fold of baseline)
NRF2  : Nrf2/ARE antioxidant capacity (fold of baseline)
MITO  : Mitochondrial functional capacity (fraction)
ATP   : Hepatocyte ATP (fraction of normal)
JNK   : Active phospho-JNK (fraction of maximum)
BAH   : Intrahepatocyte bile acids (uM)
BAP   : Serum total bile acids (uM)
HEP   : Viable hepatocyte mass (fraction of normal liver)
HEPS  : Stressed but recoverable hepatocyte mass (fraction)
NECR  : Cumulative necrotic hepatocyte mass (fraction)
DAMP  : Damage-associated molecular patterns (au)
KC    : Activated Kupffer cells (fraction of maximum)
TNF   : TNF-alpha (au)
IL10  : IL-10 (au)
TCELL : Drug-specific effector T cells (fraction of maximum)
TREG  : Regulatory T-cell / tolerance tone (fraction, 1 = intact)
ALT   : Serum ALT (U/L)
AST   : Serum AST (U/L)
TBIL  : Serum total bilirubin (mg/dL)
ALP   : Serum alkaline phosphatase (U/L)
FV    : Hepatic clotting-factor activity (fraction, 1 = normal)
MIR   : Serum miR-122 (fold of baseline)
NACC  : N-acetylcysteine in central compartment (umol)

$GLOBAL
#define CcP   (AC/VC)
#define CpP   (AP/VP)
#define CLIVP (ALIV/VLIV)
#define CuP   (FULIV*CLIVP)
#define CNAC  (NACC/VNAC)

// Hill function with a hard floor at zero — used for every threshold-like
// term so that no term can go negative when a state overshoots slightly.
double hillf(double x, double k, double n) {
  if (x <= 0.0) return 0.0;
  double xn = pow(x, n);
  return xn / (pow(k, n) + xn);
}

$MAIN
PAPS_0  = 1.0;
UDPGA_0 = 1.0;
GSH_0   = GSH0 * FGSH;
CYS_0   = 1.0;
ROS_0   = 1.0;
NRF2_0  = 1.0;
MITO_0  = 1.0;
ATP_0   = 1.0;
BAH_0   = BAH0;
BAP_0   = BAP0;
HEP_0   = 1.0;
TREG_0  = 1.0;
ALT_0   = ALT0;
AST_0   = AST0;
TBIL_0  = TBIL0;
ALP_0   = ALP0;
FV_0    = 1.0;
MIR_0   = MIR0;

$ODE
// ---- guards: keep states physical without changing any rate law ----------
double gsh  = (GSH  > 1e-9) ? GSH  : 1e-9;
double atp  = (ATP  > 1e-6) ? ATP  : 1e-6;
double mito = (MITO > 0.0 ) ? MITO : 0.0;
double hep  = (HEP  > 0.0 ) ? HEP  : 0.0;
double heps = (HEPS > 0.0 ) ? HEPS : 0.0;

// =========================================================================
// 1. PHARMACOKINETICS
// =========================================================================
double v_sult = VMAX_SULT * CuP / (KM_SULT + CuP) * PAPS;
double v_ugt  = VMAX_UGT  * CuP / (KM_UGT  + CuP) * UDPGA;
double v_cyp  = VMAX_CYP  * CuP / (KM_CYP  + CuP) * FCYP;

dxdt_AGUT = -KA * AGUT;
dxdt_AC   =  KA*AGUT*FA - QH*CcP + QH*CLIVP/KP - CLR*CcP - QD*CcP + QD*CpP;
dxdt_AP   =  QD*CcP - QD*CpP;
dxdt_ALIV =  QH*CcP - QH*CLIVP/KP - (v_sult + v_ugt + v_cyp);

// Cofactor depletion. This is what makes the *fraction* of drug shunted to
// bioactivation rise with dose — without it there is no dose nonlinearity.
dxdt_PAPS  = KSP * (1.0 - PAPS)  - v_sult / APAPS;
dxdt_UDPGA = KSU * (1.0 - UDPGA) - v_ugt  / AUD;

// =========================================================================
// 2. REACTIVE METABOLITE — the competition that decides everything
//    RM partitions between GSH conjugation and protein binding purely by
//    their relative rates. At GSH = 6.5 mM, ~0.6% binds protein; at
//    GSH = 0.1 mM, ~25% does. No switch, no threshold parameter.
// =========================================================================
double prod_RM = v_cyp / VLIV;              // uM/h
double v_conj  = KGST * gsh * RM;           // uM/h
double v_bind  = KBIND * RM;                // uM/h
dxdt_RM = prod_RM - v_conj - v_bind - KNQO * RM;

// =========================================================================
// 3. GLUTATHIONE AND CYSTEINE — the disposal rate, and its ceiling
// =========================================================================
double enac      = EMAX_NAC * CNAC / (EC50_NAC + CNAC);
double v_gsh_syn = VMAX_GSH * FGSH * CYS / (KCYS + CYS) * NRF2;

dxdt_GSH = v_gsh_syn - KDEGG*gsh - v_conj/1000.0
           - KOX * ((ROS > 1.0) ? (ROS - 1.0) : 0.0) * gsh;
dxdt_CYS = KCYSIN * (CYSBASE*(1.0 + enac) - CYS) - ALPHA_CYS * v_gsh_syn;

// =========================================================================
// 4. ADDUCTS, ROS, Nrf2
// =========================================================================
double autoph = 1.0 + 1.5 * ((atp < 1.0) ? (1.0 - atp) : 0.0);
dxdt_ADD = v_bind * FMITO - KADD_REP * autoph * ADD;

double ba_ex = (BAH > BACRIT) ? (BAH - BACRIT)/BACRIT : 0.0;
double ros_prod = ROS_BASE + KROS_MITO*(1.0 - mito) + KROS_SAB*JNK
                  + KROS_KC*KC + KROS_BA*ba_ex;
double redox_cap = (0.30 + 0.70 * gsh/(GSH0*FGSH)) * NRF2;
double ros_el = KROS_EL * ROS * redox_cap + KROS_NAC * CNAC * ROS;
dxdt_ROS = ros_prod - ros_el;

double nrf_drive = 1.0 + ENRF * hillf(ROS - 1.0, KNRF_H, 1.0);
dxdt_NRF2 = KNRF * (nrf_drive - NRF2);

// =========================================================================
// 5. MITOCHONDRIA AND ATP
// =========================================================================
double dmg = KMD_ADD*ADD + KMD_ROS*((ROS>1.0)?(ROS-1.0):0.0)
             + KMD_JNK*JNK + KMD_BA*ba_ex;
dxdt_MITO = KMR * (1.0 - mito) * (0.2 + 0.8*atp) - dmg * mito;
dxdt_ATP  = KATP * (mito * (1.0 - FUNC) - ATP);

// =========================================================================
// 6. THE JNK-Sab LOOP
//    JNK is activated by ROS (Hill, n = 2) and feeds KROS_SAB*JNK straight
//    back into ros_prod above. That single line is the bistability.
// =========================================================================
double jact = KJ_ON * hillf(ROS - 1.0, KJ_K, KJ_HILL) + KJ_TNF*TNF/(KTNF_J + TNF);
dxdt_JNK = jact * (1.0 - JNK) - KJ_OFF * JNK;

// =========================================================================
// 7. BILE ACIDS
//    BSEP flux is scaled by ATP (it is an ABC transporter) and by surviving
//    mass, so hepatocellular injury produces secondary cholestasis without
//    any extra term.
// =========================================================================
double fxr      = hillf(BAH, KFXR, 2.0);
double v_bsep   = VMAX_BSEP * BAH/(KM_BSEP + BAH) / (1.0 + CuP/KI_BSEP) * atp * hep;
double v_mrp    = KMRP * BAH * (1.0 + EFXR_MRP*fxr);
double v_syn_ba = KSYN_BA * (1.0 - FMAX_FXR*fxr) * hep;
double v_upt    = KUPT * BAP * hep;
dxdt_BAH = v_syn_ba + v_upt - v_bsep - v_mrp;
dxdt_BAP = FCONV_BA * v_mrp * VLIV - KEL_BA * BAP;

// =========================================================================
// 8. ADAPTIVE IMMUNITY (idiosyncratic arm)
// =========================================================================
// NOTE: inside $ODE the integration clock is SOLVERTIME, not TIME (TIME is the
// timestamp of the current data record and would make this step fire late).
double ster_on = (STER_T0 < 1e8)
                 ? 1.0/(1.0 + exp(-(SOLVERTIME - STER_T0)/STER_TAU)) : 0.0;
double ster_s = 1.0 - STER * ster_on * STER_EFF;
dxdt_TREG = KTREG * ((1.0 - ICI) - TREG)
            + KTREG_IL10 * IL10 * ((TREG < 1.0) ? (1.0 - TREG) : 0.0);
double t_act = KT_ON * HLA * ADD
               * ((DAMP + DANGER0)/(KD_DAMP + DAMP + DANGER0))
               / (1.0 + TREG/KTREG_S) * ster_s;
dxdt_TCELL = t_act * (1.0 - TCELL) - KT_OFF * TCELL;
double tc_eff = TCELL * ster_s;

// =========================================================================
// 9. HEPATOCYTE LIFE CYCLE
//    MPT opening needs BOTH sustained JNK and ATP depletion. With ATP intact
//    the ATP gate holds the pore shut (and the cell dies by apoptosis
//    instead, which the model treats as a slow, non-DAMP-releasing route).
// =========================================================================
double stress  = KS_ADD*ADD + KS_BA*ba_ex + KS_TNF*TNF + KS_KILL*tc_eff;
double recover = KREC * heps * (0.2 + 0.8*gsh/(GSH0*FGSH)) / (1.0 + JNK/0.2);
double mpt     = hillf(JNK, KJ_MPT, 3.0) / (1.0 + pow(atp/KATP_MPT, 4.0));
double death   = KNEC * heps * (mpt + KBA_DEATH*ba_ex + KTC_DEATH*tc_eff);
double prim    = 0.3 + 0.7 * TNF/(0.3 + TNF);
double space   = 1.0 - hep - heps;
double regen   = KREG * AGEF * hep * ((space > 0.0) ? space : 0.0) * prim;

dxdt_HEP  = -stress*hep + recover + regen;
dxdt_HEPS =  stress*hep - recover - death;
dxdt_NECR =  death;

// =========================================================================
// 10. INNATE IMMUNITY
// =========================================================================
dxdt_DAMP = KDAMP*death - KDAMP_EL*DAMP;
dxdt_KC   = KKC*DAMP*(1.0 - KC) - KKC_OFF*KC*(1.0 + IL10/KIL10);
dxdt_TNF  = KTNF*KC - KTNF_EL*TNF*(1.0 + IL10);
dxdt_IL10 = KIL10P*KC - KIL10_EL*IL10;

// =========================================================================
// 11. BIOMARKERS
//     ALT is a low-pass filter on the RATE of death.
//     TBIL is set by the RESERVE that survives (clearance scales with mass).
//     They are structurally different quantities — that is Hy s Law.
// =========================================================================
double chol_inj = KCHOL * ba_ex;
dxdt_ALT = KALT_EL*ALT0 + KALT_REL*death - KALT_EL*ALT;
dxdt_AST = KAST_EL*AST0 + KAST_REL*death - KAST_EL*AST;
dxdt_ALP = KALP_EL*ALP0 + KALP_REL*chol_inj - KALP_EL*ALP;

double fchol     = (1.0/(1.0 + BAH/KBIL_BA)) * (1.0/(1.0 + TNF/KBIL_TNF));
double kbil_prod = (KBIL_CLR * (1.0/(1.0 + BAH0/KBIL_BA)) + KBIL_ALT) * TBIL0;
dxdt_TBIL = kbil_prod - KBIL_CLR * hep * fchol * TBIL - KBIL_ALT * TBIL;

dxdt_FV  = KFV_EL*hep - KFV_EL*FV;
dxdt_MIR = KMIR_EL*MIR0 + KMIR_REL*death - KMIR_EL*MIR;

// =========================================================================
// 12. ANTIDOTE PK
// =========================================================================
dxdt_NACC = -CLNAC/VNAC * NACC;

$TABLE
double CP      = AC/VC;                       // plasma drug (uM)
double CP_UGML = CP * 151.16 / 1000.0;        // ug/mL if MW = APAP
double CLIVER  = ALIV/VLIV;
double CU      = FULIV * CLIVER;
double GSHPCT  = 100.0 * GSH/(GSH0*FGSH);
double SURV    = HEP + HEPS;                  // surviving hepatocyte mass
double LOST    = 1.0 - SURV;
double INR     = 1.0 + 0.333*(1.0/((FV > 0.001) ? FV : 0.001) - 1.0);
if (INR > 20.0) INR = 20.0;                   // >=20 is not a survivable state
double RRATIO  = (ALT/ULN_ALT) / (ALP/ULN_ALP);
double HYLAW   = ((ALT >= 3.0*ULN_ALT) && (TBIL >= 2.0*ULN_TBIL)) ? 1.0 : 0.0;
double REDOXG  = (0.30 + 0.70*GSH/(GSH0*FGSH)) * NRF2;   // loop redox capacity
double FBIOACT = (VMAX_CYP*CU/(KM_CYP+CU)*FCYP)
                 / (VMAX_CYP*CU/(KM_CYP+CU)*FCYP
                    + VMAX_SULT*CU/(KM_SULT+CU)*PAPS
                    + VMAX_UGT*CU/(KM_UGT+CU)*UDPGA + 1e-12);

$CAPTURE @annotated
CP      : Plasma drug concentration (uM)
CP_UGML : Plasma drug concentration (ug/mL, MW 151.16)
CU      : Unbound intrahepatic drug (uM)
GSHPCT  : Hepatic GSH (% of baseline)
SURV    : Surviving hepatocyte mass (fraction)
LOST    : Lost hepatocyte mass (fraction)
INR     : International normalised ratio (capped at 20)
RRATIO  : R ratio, (ALT/ULN)/(ALP/ULN); >=5 hepatocellular, <=2 cholestatic
HYLAW   : Hy s Law satisfied (1/0)
REDOXG  : Redox capacity g governing JNK-loop bistability
FBIOACT : Fraction of hepatic metabolism going to the reactive metabolite
'

mod <- mcode("dili", dili_code)

# =============================================================================
# Helpers
# =============================================================================
MW_APAP <- 151.16
MW_NAC  <- 163.20

mgkg_to_umol <- function(mgkg, wt = 70, mw = MW_APAP) mgkg * wt * 1000 / mw

# Events are built as a plain data frame (ID/time/amt/rate/cmt/evid) rather
# than by rbind-ing ev() objects, so this works across mrgsolve versions.
CMT_AGUT <- 1L    # position of AGUT in $CMT
CMT_NACC <- 33L   # position of NACC in $CMT

#' Oral dose event(s) into the gut compartment.
oral <- function(mgkg, times = 0, wt = 70, mw = MW_APAP) {
  data.frame(ID = 1, time = times, amt = mgkg_to_umol(mgkg, wt, mw),
             rate = 0, cmt = CMT_AGUT, evid = 1)
}

#' 21-hour Prescott IV N-acetylcysteine regimen: 150 mg/kg over 1 h,
#' then 50 mg/kg over 4 h, then 100 mg/kg over 16 h.  [ref 78]
#' `maint_h` prolongs the maintenance infusion AT THE STANDARD RATE
#' (100 mg/kg per 16 h), i.e. it gives more NAC rather than spreading the
#' same dose thinner. 16 = the standard 21-hour course.
nac_prescott <- function(start, wt = 70, maint_h = 16) {
  if (is.null(start) || length(start) == 0 || is.na(start)) return(NULL)
  a <- function(mgkg) mgkg * wt * 1000 / MW_NAC
  data.frame(ID = 1,
             time = c(start, start + 1, start + 5),
             amt  = c(a(150), a(50),   a(100) * maint_h/16),
             rate = c(a(150)/1, a(50)/4, a(100)/16),
             cmt  = CMT_NACC, evid = 1)
}

#' Run one scenario. Extra named arguments are passed straight to param().
sim_dili <- function(dose_mgkg = 0, dose_times = 0, nac_start = NULL,
                     nac_maint_h = 16, tend = 336, delta = 0.25, wt = 70,
                     ...) {
  e <- if (dose_mgkg > 0)
    oral(dose_mgkg / length(dose_times), dose_times, wt) else NULL
  n <- nac_prescott(nac_start, wt, nac_maint_h)
  e <- if (is.null(e)) n else if (is.null(n)) e else rbind(e, n)
  m <- mod
  pp <- list(...)
  if (length(pp)) m <- param(m, pp)
  if (is.null(e)) mrgsim(m, end = tend, delta = delta)
  else mrgsim(data_set(m, e[order(e$time), ]), end = tend, delta = delta)
}

# =============================================================================
# Drug archetypes — the ONLY thing that differs between injury phenotypes.
# There is no parameter named "pattern", "cholestatic", or "severity".
# =============================================================================
SLOW_CL <- list(VMAX_UGT = 2500, VMAX_SULT = 200, CLR = 0.5)

# Drug B: low-clearance, OATP-concentrated BSEP inhibitor, weak bioactivator
#         (troglitazone / bosentan archetype)
DRUG_B <- c(SLOW_CL, list(KP = 25, FULIV = 0.05, KI_BSEP = 0.5, VMAX_CYP = 150))

# Drug C: strong bioactivator with very slow adduct turnover, no BSEP liability
#         (flucloxacillin / amoxicillin-clavulanate archetype)
DRUG_C <- c(SLOW_CL, list(KP = 6, FULIV = 0.5, VMAX_CYP = 400,
                          KADD_REP = 0.001))

# Vulnerable host: chronic ethanol (CYP2E1 induction) plus fasting
# (low cysteine supply and a reduced GSH pool)
VULNERABLE <- list(FCYP = 2.2, CYSBASE = 1.05, FGSH = 0.6)

# =============================================================================
# THIRTEEN SCENARIOS
# =============================================================================
scenarios <- function() {
  q6h_7d  <- seq(0, 162, by = 6)
  q12_28d <- seq(0, 660, by = 12)
  q24_56d <- seq(0, 1320, by = 24)

  list(
    S1  = list(label = "APAP therapeutic 1 g q6h x 7 d",
               args = list(dose_mgkg = 400, dose_times = q6h_7d, tend = 504)),
    S2  = list(label = "APAP 150 mg/kg acute, untreated",
               args = list(dose_mgkg = 150, tend = 336)),
    S3  = list(label = "APAP 250 mg/kg acute, untreated",
               args = list(dose_mgkg = 250, tend = 336)),
    S4  = list(label = "APAP 350 mg/kg acute, untreated",
               args = list(dose_mgkg = 350, tend = 336)),
    S5  = list(label = "APAP 350 mg/kg + NAC at 8 h",
               args = list(dose_mgkg = 350, nac_start = 8, tend = 336)),
    S6  = list(label = "APAP 350 mg/kg + NAC at 16 h",
               args = list(dose_mgkg = 350, nac_start = 16, tend = 336)),
    S7  = list(label = "APAP 350 mg/kg + NAC at 24 h",
               args = list(dose_mgkg = 350, nac_start = 24, tend = 336)),
    S8  = list(label = "APAP 150 mg/kg, ethanol + fasting host",
               args = c(list(dose_mgkg = 150, tend = 336), VULNERABLE)),
    S9  = list(label = "APAP 150 mg/kg, same host + NAC at 8 h",
               args = c(list(dose_mgkg = 150, nac_start = 8, tend = 336),
                        VULNERABLE)),
    S10 = list(label = "Drug B (BSEP inhibitor) 8 mg/kg q12h x 28 d",
               args = c(list(dose_mgkg = 8*28, dose_times = q12_28d,
                             tend = 1400), DRUG_B)),
    S11 = list(label = "Drug C (HLA risk allele) 10 mg/kg qd x 56 d",
               args = c(list(dose_mgkg = 10*56, dose_times = q24_56d,
                             tend = 2400, HLA = 1), DRUG_C)),
    S12 = list(label = "Drug C + checkpoint inhibitor (tolerance removed)",
               args = c(list(dose_mgkg = 10*56, dose_times = q24_56d,
                             tend = 2400, HLA = 1, ICI = 1), DRUG_C)),
    S13 = list(label = "Drug C + ICI + corticosteroid from day 42",
               args = c(list(dose_mgkg = 10*56, dose_times = q24_56d,
                             tend = 2400, HLA = 1, ICI = 1, STER = 1,
                             STER_T0 = 1008), DRUG_C))
  )
}

summarise_run <- function(s) {
  d <- as.data.frame(s)
  data.frame(
    peak_ALT  = max(d$ALT),
    peak_AST  = max(d$AST),
    peak_ALP  = max(d$ALP),
    peak_TBIL = max(d$TBIL),
    peak_INR  = max(d$INR),
    R_ratio   = max(d$ALT)/40 / (max(d$ALP)/120),
    hy_law    = as.integer(any(d$HYLAW > 0)),
    GSH_nadir = min(d$GSHPCT),
    JNK_peak  = max(d$JNK),
    lost_mass = max(d$LOST),
    peak_BA   = max(d$BAP),
    peak_miR  = max(d$MIR)
  )
}

run_all <- function() {
  sc <- scenarios()
  do.call(rbind, lapply(names(sc), function(k) {
    r <- summarise_run(do.call(sim_dili, sc[[k]]$args))
    cbind(id = k, label = sc[[k]]$label, r)
  }))
}

# =============================================================================
# ANALYSES — the four results the model is built to compute
# =============================================================================

#' [A] Dose-response: locate the OFF -> ON transition. There is no threshold
#'     parameter in the model; the knee comes from GSH exhaustion feeding the
#'     bistable JNK-Sab loop.
analysis_dose_threshold <- function(doses = seq(100, 400, by = 10)) {
  do.call(rbind, lapply(doses, function(d) {
    r <- summarise_run(sim_dili(dose_mgkg = d, tend = 336))
    cbind(dose_mgkg = d, r)
  }))
}

#' [B] Rate vs dose: hold the total dose fixed and vary only the delivery
#'     rate. If DILI were a dose problem these rows would be identical.
analysis_rate_vs_dose <- function(dose = 350,
                                  spread_h = c(0, 3, 6, 12, 18, 24, 36, 48)) {
  do.call(rbind, lapply(spread_h, function(h) {
    tt <- if (h == 0) 0 else seq(0, h, length.out = max(2, h))
    r  <- summarise_run(sim_dili(dose_mgkg = dose, dose_times = tt, tend = 500))
    cbind(spread_h = h, r)
  }))
}

#' [C] NAC window: sweep the start time in three hosts. The "cliff" is not a
#'     parameter; it is the moment the trajectory crosses the separatrix.
analysis_nac_window <- function(times = c(2,4,6,8,10,12,14,16,20,24,32,48)) {
  hosts <- list(
    "350 mg/kg, normal host"        = list(dose_mgkg = 350),
    "500 mg/kg, normal host"        = list(dose_mgkg = 500),
    "150 mg/kg, ethanol + fasting"  = c(list(dose_mgkg = 150), VULNERABLE))
  do.call(rbind, lapply(names(hosts), function(h) {
    do.call(rbind, lapply(c(times, NA), function(tn) {
      a <- c(hosts[[h]], list(nac_start = if (is.na(tn)) NULL else tn,
                              tend = 336))
      cbind(host = h, nac_start = tn, summarise_run(do.call(sim_dili, a)))
    }))
  }))
}

#' [D] Bistability of the isolated JNK-Sab loop as a function of the redox
#'     capacity g = (0.30 + 0.70*GSH/GSH0) * NRF2.  Returns the fixed points.
analysis_bistability <- function(gs = seq(0.6, 1.6, by = 0.05)) {
  p <- as.list(param(mod))
  fp <- function(g) {
    J <- seq(0, 1, length.out = 40001)
    ros <- (p$ROS_BASE + p$KROS_SAB*J) / (p$KROS_EL*g)
    x <- pmax(ros - 1, 0)
    act <- p$KJ_ON * x^p$KJ_HILL / (p$KJ_K^p$KJ_HILL + x^p$KJ_HILL)
    f <- act*(1 - J) - p$KJ_OFF*J
    J[which(diff(sign(f)) != 0)]
  }
  do.call(rbind, lapply(gs, function(g) {
    r <- fp(g)
    data.frame(g = g, n_fixed_points = length(r),
               separatrix = if (length(r) >= 2) r[1] else NA,
               on_state   = if (length(r) >= 1) r[length(r)] else NA,
               regime = if (length(r) >= 2) "bistable"
                        else if (length(r) == 1 && r[1] > 0.3) "monostable-ON"
                        else "monostable-OFF")
  }))
}

#' [E] Hy s Law anatomy: at what lost-mass does each arm cross?
analysis_hy_law <- function(doses = seq(100, 600, by = 10)) {
  d <- do.call(rbind, lapply(doses, function(x) {
    r <- summarise_run(sim_dili(dose_mgkg = x, tend = 336))
    cbind(dose_mgkg = x, r)
  }))
  alt_arm <- d[d$peak_ALT >= 120, ][1, ]
  bil_arm <- d[d$peak_TBIL >= 2.4, ][1, ]
  list(table = d,
       alt_arm_first_at = alt_arm,
       bil_arm_first_at = bil_arm,
       reserve_to_rate_ratio = bil_arm$lost_mass / alt_arm$lost_mass)
}

# =============================================================================
# Example session
# =============================================================================
if (interactive()) {
  print(run_all())
  print(analysis_dose_threshold())
  print(analysis_rate_vs_dose())
  print(analysis_nac_window())
  print(analysis_bistability())
  str(analysis_hy_law()$reserve_to_rate_ratio)

  # a single trajectory
  s <- sim_dili(dose_mgkg = 350, nac_start = 8, tend = 336)
  plot(s, ALT + TBIL + GSHPCT + JNK + ATP + LOST ~ time)
}
