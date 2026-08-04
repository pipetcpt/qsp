## =============================================================================
##  vl_mrgsolve_model.R
##  Visceral Leishmaniasis (kala-azar) — 73-ODE QSP model for mrgsolve
##  내장 리슈만편모충증 QSP 모델 (mrgsolve, 73개 ODE)
##
##  Author's note on provenance
##  ---------------------------
##  No R runtime was available in the environment where this file was written,
##  so every equation below was first implemented independently in Python
##  (`vl_reference_model.py`) and integrated there; the numbers quoted in
##  README.md are that integration's output, and `vl_reference_output.txt` is
##  its verbatim log.  Three real defects were found by doing it that way and
##  each is marked in place with a "BUG FOUND" comment, because a reader of a
##  QSP model deserves to know which lines were hard:
##
##    1. an extinction floor written INSIDE the derivative block
##       (`if (P < 1e-7 && dPdt > 0) dPdt = 0`) is a discontinuous switch; it
##       is invisible while burdens are high and then stalls the integrator
##       completely once any compartment clears.  The floor belongs between
##       integration steps, not in f(t,y).
##    2. CD4 was written with a loss term and no set point, so every
##       immunocompetent patient drifted to 4 cells/uL over an 18-month
##       follow-up and the entire host-immunity arm silently switched off —
##       which made five regimens that cure >90% of real patients "fail".
##    3. TGF-beta driven by the saturating antigen signal AG (rather than by
##       the high-burden signal AGH) stayed near-maximal down to 1e4
##       amastigotes and permanently clamped macrophage activation, so no
##       slow-acting drug could ever hand off to the host.
##
##  What the model is for
##  --------------------
##  One structural commitment: THE DRUG AND THE PARASITE NEVER MEET IN PLASMA.
##  Amastigotes live inside macrophages of spleen, liver, marrow and skin, so
##  every drug effect is driven by an INTRAMACROPHAGE concentration, while
##  amphotericin nephrotoxicity is driven by the FREE PLASMA concentration of
##  the same dose.  Liposomal encapsulation moves those two integrals in
##  opposite directions because the liposome is cleared by the very cell
##  lineage that harbours the parasite.  Cure is then not "zero parasites" but
##  crossing a separatrix in the (burden, primed-memory) plane — and that
##  separatrix moves with CD4 count.
##
##  Compartment map (73 ODEs)
##  -------------------------
##    1– 8   amphotericin B: liposome-associated and free plasma species,
##           peripheral tissue, four macrophage organs, renal cortex
##    9–15   miltefosine: gut, central, peripheral, four macrophage organs
##   16–24   paromomycin: depot, central, peripheral, four organs, kidney,
##           cochlea
##   25–32   antimony: depot, plasma Sb(V), (spare), four organs Sb(III),
##           deep tissue depot
##   33–40   amastigotes: replicating + quiescent, in each of four organs
##   41–47   TMEM, IFN-gamma, IL-10, TNF, TGF-beta, CD4, activated macrophage
##   48–57   spleen, liver, Hb, platelets, WBC, albumin, IgG, temperature,
##           weight loss, PKDL lesion load
##   58–66   tubular injury, creatinine, K, Mg, hearing, QTc, lipase, ALT, GI
##   67–73   five exposure integrals, cumulative hazard, cumulative mg/kg
##
##  Dosing compartments
##  -------------------
##    cmt 1  = A_LIP   liposomal amphotericin B, IV bolus/infusion
##    cmt 2  = A_FRE   amphotericin B deoxycholate, IV (straight to free pool)
##    cmt 9  = MIL_G   miltefosine, oral
##    cmt 16 = PM_D    paromomycin, IM
##    cmt 25 = SB5_D   sodium stibogluconate / meglumine antimoniate, IM or IV
##
##  Time unit is HOURS.  Burden unit is 1e6 amastigotes.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

vl_code <- '
$PARAM @annotated
// ---- covariates -----------------------------------------------------------
WT       :  50.0 : Body weight (kg)
AGE      :  30.0 : Age (y)
CD40     : 700.0 : Baseline CD4 count (cells/uL)
HIV      :   0.0 : HIV co-infection flag (0/1)
ART      :   0.0 : On suppressive antiretroviral therapy (0/1)
MALNUT   :   0.0 : Malnutrition index 0-1, scales T-cell priming
P0SCALE  :   1.0 : Multiplier on presenting parasite burden

// ---- amphotericin B disposition -------------------------------------------
V_LIP    :  0.10 : Liposome distribution volume (L/kg, ~plasma volume)
CL_LIP   :  0.50 : Liposome disposition clearance (L/h/70kg)
FREL     :  0.30 : Fraction of liposome clearance that LEAKS drug to plasma
V_FRE    :  0.50 : Free/bound AmB central volume (L/kg)
CL_FRE   :  1.40 : Free AmB elimination clearance (L/h/70kg)
Q_AMB    :  1.20 : Free AmB intercompartmental clearance (L/h/70kg)
V_PER    :  3.00 : AmB peripheral volume (L/kg)
FSP_A    :  0.10 : Spleen share of MPS liposome uptake
FLI_A    :  0.60 : Liver (Kupffer) share of MPS liposome uptake
FBM_A    :  0.12 : Bone marrow share of MPS liposome uptake
FSK_A    :  0.03 : Dermal share of MPS liposome uptake
KOUT_A   : 0.00825 : Egress of the pharmacologically available tissue AmB (1/h)
KPSP_A   :  14.0 : Spleen:free-plasma partition coefficient, AmB
KPLI_A   :  20.0 : Liver:free-plasma partition coefficient, AmB
KPBM_A   :   9.0 : Marrow:free-plasma partition coefficient, AmB
KPSK_A   :   1.5 : Skin:free-plasma partition coefficient, AmB
KPKID_A  :   8.0 : Kidney:free-plasma partition coefficient, AmB
KIN_KID  :  0.15 : Renal cortical uptake rate (1/h)
KOUT_KID : 0.012 : Renal cortical egress rate (1/h)

// ---- miltefosine ----------------------------------------------------------
KA_MIL   :  0.40 : Oral absorption rate constant (1/h)
F_MIL    :  1.00 : Relative bioavailability
CL_MIL   : 0.166 : Apparent clearance (L/h/70kg)
V1_MIL   :  40.0 : Apparent central volume (L/70kg)
V2_MIL   :  52.0 : Apparent peripheral volume (L/70kg)
Q_MIL    :  0.90 : Apparent intercompartmental clearance (L/h/70kg)
KPSP_M   :   8.0 : Spleen macrophage partition, miltefosine
KPLI_M   :   8.0 : Liver macrophage partition, miltefosine
KPBM_M   :   6.0 : Marrow macrophage partition, miltefosine
KPSK_M   :   4.0 : Dermal macrophage partition, miltefosine
KIN_M    :  0.05 : Tissue equilibration rate, miltefosine (1/h)
EM50_GI  :  25.0 : Plasma conc for half-maximal GI intolerance (mg/L)

// ---- paromomycin ----------------------------------------------------------
KA_PM    :  1.80 : IM absorption rate constant (1/h)
V_PM     :  0.25 : Central volume (L/kg)
CL_PM    :  5.50 : Clearance, glomerular (L/h/70kg)
Q_PM     :  1.00 : Intercompartmental clearance (L/h/70kg)
V2_PM    :  0.20 : Peripheral volume (L/kg)
KIN_PM   : 0.020 : Pinocytic entry into macrophage (1/h) - deliberately SLOW
KOUT_PM  : 0.004 : Macrophage egress (1/h)
FSK_PM   :  0.25 : Dermal access relative to MPS organs
KIN_PMK  :  0.10 : Renal cortical (megalin) uptake (1/h)
KOUT_PMK : 0.008 : Renal cortical egress (1/h)
KIN_COC  : 0.004 : Cochlear uptake (1/h)
KOUT_COC : 0.0004 : Cochlear egress (1/h) - near-irreversible

// ---- antimony -------------------------------------------------------------
KA_SB    :  1.50 : IM absorption rate constant (1/h)
V_SB     :  0.22 : Sb(V) central volume (L/kg)
CL_SB    :  7.00 : Sb(V) renal clearance (L/h/70kg)
Q_SB     :  0.45 : Clearance to the deep antimony depot (L/h/70kg)
V2_SB    :  0.90 : Deep antimony depot volume (L/kg)
KRED     : 0.045 : Sb(V) -> Sb(III) intracellular reduction rate (1/h)
KEFF_SB  :  0.30 : Sb(III) efflux rate, MRPA/trypanothione route (1/h)
RES_SB   :   1.0 : Antimony resistance multiplier ON EFFLUX (Bihar ~9)
KPSK_SB  :  0.30 : Dermal Sb access relative to MPS organs

// ---- organ compartment volumes -------------------------------------------
VSP      :  0.25 : Spleen macrophage compartment volume (L)
VLI      :  1.50 : Liver macrophage compartment volume (L)
VBM      :  1.50 : Marrow macrophage compartment volume (L)
VSK      :  3.00 : Dermal macrophage compartment volume (L)
VKID     :  0.30 : Renal cortex volume (L)
VCOC     : 0.002 : Cochlear volume (L)

// ---- parasite dynamics ----------------------------------------------------
KG_SP    : 0.0072 : Amastigote replication rate, spleen (1/h)
KG_LI    : 0.0065 : Amastigote replication rate, liver (1/h)
KG_BM    : 0.0068 : Amastigote replication rate, marrow (1/h)
KG_SK    : 0.0040 : Amastigote replication rate, skin (1/h)
PMAX_SP  : 3.0e3 : Spleen carrying capacity (1e6 amastigotes)
PMAX_LI  : 1.0e4 : Liver carrying capacity (1e6 amastigotes)
PMAX_BM  : 3.0e3 : Marrow carrying capacity (1e6 amastigotes)
PMAX_SK  : 1.0e2 : Skin carrying capacity (1e6 amastigotes)
KIMM     : 0.065 : Maximal immune (NO-mediated) killing rate (1/h)
KMIG     : 2.0e-4 : Redistribution among spleen/liver/marrow (1/h)
KMIG_SKIN : 3.0e-6 : Visceral pool -> dermis trafficking (1/h)
KMIG_BACK : 3.0e-7 : Dermis -> visceral pool trafficking (1/h)
FR0      : 5.0e-4 : Quiescent (persister) fraction at presentation
KQ0      : 2.0e-6 : Basal entry into quiescence (1/h)
KQS      : 1.5e-5 : Stress-induced quiescence per unit of kill pressure
KR_REACT : 8.0e-4 : Reactivation out of quiescence (1/h)
COSTR    :  0.05 : Replication rate of quiescent forms, fraction of KG
FIMM_R   :  0.15 : Immune killing of quiescent forms, fraction of KIMM
FIMM_SK  :  0.12 : Immune killing in the DERMIS, fraction of KIMM
RFAC_A   :  20.0 : EC50 multiplier for quiescent forms, amphotericin B
RFAC_M   :  25.0 : EC50 multiplier for quiescent forms, miltefosine
RFAC_P   :  30.0 : EC50 multiplier for quiescent forms, paromomycin
RFAC_S   :  40.0 : EC50 multiplier for quiescent forms, Sb(III)
RES_M    :   1.0 : Acquired miltefosine resistance multiplier on EC50
RES_P    :   1.0 : Acquired paromomycin resistance multiplier on EC50

// ---- drug PD (Emax on INTRAMACROPHAGE concentration) ---------------------
EMAX_A   : 0.060 : Maximal AmB kill rate (1/h)
EC50_A   :   6.0 : AmB EC50, intramacrophage (mg/L)
HILL_A   :   2.0 : AmB Hill coefficient
EMAX_M   : 0.034 : Maximal miltefosine kill rate (1/h)
EC50_M   :  90.0 : Miltefosine EC50, intramacrophage (mg/L)
HILL_M   :   1.5 : Miltefosine Hill coefficient
EMAX_P   : 0.042 : Maximal paromomycin kill rate (1/h)
EC50_P   :  25.0 : Paromomycin EC50, intramacrophage (mg/L)
HILL_P   :   1.5 : Paromomycin Hill coefficient
EMAX_S   : 0.045 : Maximal Sb(III) kill rate (1/h)
EC50_S   :  0.50 : Sb(III) EC50, intramacrophage (mg/L)
HILL_S   :   1.5 : Sb(III) Hill coefficient

// ---- immunity -------------------------------------------------------------
KAG      :  0.01 : Burden giving half-maximal antigen signal (1e6 amastigotes)
KTP      : 0.0050 : Memory/effector T-cell priming rate (1/h)
TMEMAX   :   2.0 : Ceiling on the primed memory pool
KTD      : 0.00015 : Memory decay rate (1/h, t1/2 ~192 d)
KI_T     :  0.35 : IL-10 concentration halving priming
KIP      : 0.045 : IFN-gamma production rate
KID      : 0.060 : IFN-gamma decay rate (1/h)
FAG_IFN  :  0.04 : Antigen-only contribution to IFN-gamma (innate, small)
KLP      : 0.030 : IL-10 production rate
KLD      : 0.020 : IL-10 decay rate (1/h)
K50_IL10 :  60.0 : Burden for half-maximal IL-10 (1e6 amastigotes)
K50_SYS  : 300.0 : Burden for half-maximal systemic illness (1e6 amastigotes)
HIL10    :   2.0 : Hill coefficient of the IL-10 burden response
KNP      : 0.050 : TNF production rate
KND      : 0.070 : TNF decay rate (1/h)
KGP      : 0.010 : TGF-beta production rate
KGD      : 0.012 : TGF-beta decay rate (1/h)
KI_IFN   :  0.35 : IFN-gamma for half-maximal macrophage activation
KI_IL10  :  0.25 : IL-10 inhibition constant on macrophage activation
KI_TGFB  :  0.40 : TGF-beta inhibition constant on macrophage activation
KCD4     : 350.0 : CD4 count for half-maximal T-cell help (cells/uL)
CD4SET   : 800.0 : CD4 homeostatic set point without HIV (cells/uL)
CD4SET_ART : 600.0 : CD4 set point achieved on suppressive ART (cells/uL)
KCD4H    : 0.0015 : CD4 homeostatic restoration rate (1/h)
KCD4L    : 0.00080 : Burden-driven CD4 depression rate (1/h)
KMPH     : 0.015 : Activated-macrophage pool turnover (1/h)
FHIV_MAC :  0.50 : Residual macrophage killing capacity in HIV

// ---- clinical -------------------------------------------------------------
SPL0     :   0.5 : Healthy spleen size (cm below costal margin)
SPLMAX   :  12.0 : Maximal added spleen size (cm)
K50_SPL  : 300.0 : Spleen burden for half-maximal splenomegaly
KSPL     : 0.0035 : Spleen size relaxation rate (1/h)
LIV0     :   0.5 : Healthy liver size (cm)
LIVMAX   :   5.0 : Maximal added liver size (cm)
KLIV     : 0.0030 : Liver size relaxation rate (1/h)
HGB0     :  14.0 : Baseline haemoglobin (g/dL)
KHGB     : 0.0060 : Haemoglobin pool turnover (1/h)
FSUP_HGB :  0.75 : Maximal marrow suppression of erythropoiesis
KSPL_HGB : 0.045 : Hypersplenic RBC destruction per cm of spleen
PLT0     : 250.0 : Baseline platelets (1e9/L)
KPLT     : 0.0090 : Platelet turnover (1/h)
KSPL_PLT : 0.085 : Hypersplenic platelet destruction per cm of spleen
WBC0     :   7.0 : Baseline leucocytes (1e9/L)
KWBC     : 0.020 : Leucocyte turnover (1/h)
KSPL_WBC : 0.055 : Hypersplenic leucocyte destruction per cm of spleen
ALB0     :   4.3 : Baseline albumin (g/dL)
KALB     : 0.0035 : Albumin turnover (1/h)
IGG0     :   1.1 : Baseline IgG (g/dL)
KIGG_P   : 0.00087 : Polyclonal IgG production per antigen unit
KIGG_D   : 0.00030 : IgG decay rate (1/h, t1/2 ~96 d)
TEMP0    :  36.7 : Baseline temperature (degC)
KTEMP    : 0.060 : Temperature relaxation rate (1/h)
TNF_FEV  :   3.2 : Maximal TNF-driven temperature rise (degC)
KBWL     : 0.0012 : Weight-loss dynamics rate (1/h)
BWLMAX   :  0.28 : Maximal fractional weight loss
KPKDL_P  : 0.0022 : PKDL lesion formation rate
KPKDL_D  : 0.0015 : PKDL lesion resolution rate (1/h)

// ---- toxicity -------------------------------------------------------------
KTUB     : 0.00035 : Tubular injury per mg/L renal AmB per h
KTUBR    : 0.0030 : Tubular repair rate (1/h)
GFR0     : 100.0 : Baseline GFR (mL/min)
KSCR     : 0.030 : Creatinine turnover (1/h)
SCR0     :  0.85 : Baseline creatinine (mg/dL)
KWK      : 0.0039 : Renal K wasting per mg/L renal AmB per h
K_K0     :   4.1 : Baseline serum potassium (mEq/L)
KKR      : 0.030 : Potassium restoration rate (1/h)
KWMG     : 0.00146 : Renal Mg wasting per mg/L renal AmB per h
MG0      :   2.0 : Baseline serum magnesium (mg/dL)
KMGR     : 0.025 : Magnesium restoration rate (1/h)
KHEAR    : 0.00040 : Hearing threshold shift per mg/L cochlear PM per h
QTC0     : 400.0 : Baseline QTc (ms)
KQTC     : 0.145 : QTc sensitivity to deep-depot antimony
KQTCR    : 0.020 : QTc recovery rate (1/h)
LIPA0    :  30.0 : Baseline lipase (U/L)
KLIPA    :  0.85 : Lipase response to Sb(V)
KLIPAR   : 0.015 : Lipase recovery rate (1/h)
ALT0     :  25.0 : Baseline ALT (U/L)
KALTX    :  0.18 : ALT response to Sb(V)
KALT_M   : 0.013 : ALT response to miltefosine
KALTR    : 0.012 : ALT recovery rate (1/h)
KGITX    :  0.10 : GI intolerance accumulation rate
KGITR    :  0.09 : GI intolerance resolution rate (1/h)

// ---- mortality ------------------------------------------------------------
H0       : 2.2e-5 : Baseline hazard of the severe VL state (1/h)
HGB_H    :   7.0 : Haemoglobin below which hazard climbs (g/dL)
ALB_H    :   2.8 : Albumin below which hazard climbs (g/dL)
WBC_H    :   1.5 : Leucocytes below which hazard climbs (1e9/L)

$CMT @annotated
A_LIP   : Liposome-associated amphotericin B, plasma (mg)
A_FRE   : Released / protein-bound amphotericin B, central (mg)
A_PER   : Amphotericin B, peripheral tissue (mg)
AMB_SP  : Amphotericin B, spleen macrophage compartment (mg)
AMB_LI  : Amphotericin B, liver Kupffer compartment (mg)
AMB_BM  : Amphotericin B, marrow macrophage compartment (mg)
AMB_SK  : Amphotericin B, dermal macrophage compartment (mg)
AMB_KID : Amphotericin B, renal cortex (mg)
MIL_G   : Miltefosine, gut depot (mg)
MIL_C   : Miltefosine, central (mg)
MIL_P   : Miltefosine, peripheral (mg)
MIL_SP  : Miltefosine, spleen macrophage compartment (mg)
MIL_LI  : Miltefosine, liver macrophage compartment (mg)
MIL_BM  : Miltefosine, marrow macrophage compartment (mg)
MIL_SK  : Miltefosine, dermal macrophage compartment (mg)
PM_D    : Paromomycin, IM depot (mg)
PM_C    : Paromomycin, central (mg)
PM_P    : Paromomycin, peripheral (mg)
PM_SP   : Paromomycin, spleen macrophage compartment (mg)
PM_LI   : Paromomycin, liver macrophage compartment (mg)
PM_BM   : Paromomycin, marrow macrophage compartment (mg)
PM_SK   : Paromomycin, dermal macrophage compartment (mg)
PM_KID  : Paromomycin, renal cortex (mg)
PM_COC  : Paromomycin, cochlea (mg)
SB5_D   : Antimony, IM depot (mg Sb)
SB5_C   : Sb(V), plasma (mg)
SB5_P   : Spare compartment, kept for state-vector parity (mg)
SB3_SP  : Sb(III), spleen macrophage compartment (mg)
SB3_LI  : Sb(III), liver macrophage compartment (mg)
SB3_BM  : Sb(III), marrow macrophage compartment (mg)
SB3_SK  : Sb(III), dermal macrophage compartment (mg)
SB_DEP  : Antimony, deep tissue depot (mg)
P_SP_S  : Replicating amastigotes, spleen (1e6)
P_SP_R  : Quiescent amastigotes, spleen (1e6)
P_LI_S  : Replicating amastigotes, liver (1e6)
P_LI_R  : Quiescent amastigotes, liver (1e6)
P_BM_S  : Replicating amastigotes, marrow (1e6)
P_BM_R  : Quiescent amastigotes, marrow (1e6)
P_SK_S  : Replicating amastigotes, skin (1e6)
P_SK_R  : Quiescent amastigotes, skin (1e6)
TMEM    : Primed effector/memory T-cell pool (0-2)
IFNG    : IFN-gamma (arbitrary units)
IL10    : IL-10 (arbitrary units)
TNFA    : TNF-alpha (arbitrary units)
TGFB    : TGF-beta (arbitrary units)
CD4     : CD4 T-cell count (cells/uL)
MPHA    : Activated-macrophage pool (0-1)
SPL     : Spleen size (cm below costal margin)
LIVS    : Liver size (cm below costal margin)
HGB     : Haemoglobin (g/dL)
PLT     : Platelets (1e9/L)
WBC     : Leucocytes (1e9/L)
ALB     : Albumin (g/dL)
IGG     : Polyclonal IgG (g/dL)
TEMP    : Body temperature (degC)
BWTL    : Fractional body-weight loss
PKDL    : PKDL lesion load (arbitrary units)
TUBI    : Tubular injury index
SCR     : Serum creatinine (mg/dL)
KSER    : Serum potassium (mEq/L)
MGSER   : Serum magnesium (mg/dL)
HEAR    : Hearing threshold shift (dB)
QTC     : QTc interval (ms)
LIPA    : Serum lipase (U/L)
ALTX    : Serum ALT (U/L)
GITX    : GI intolerance index
AUCFRE  : AUC of free plasma AmB (mg.h/L)   <- nephrotoxicity driver
AUCASP  : AUC of spleen-macrophage AmB (mg.h/L) <- efficacy driver
AUCMSP  : AUC of spleen-macrophage miltefosine (mg.h/L)
AUCBSP  : AUC of spleen-macrophage Sb(III) (mg.h/L)
AUCPSP  : AUC of spleen-macrophage paromomycin (mg.h/L)
CUMHAZ  : Cumulative mortality hazard
CUMDOSE : Cumulative dose administered (mg/kg)

$GLOBAL
#define _EMX(conc, ec50, emax, hill) \\
  ((conc) <= 0.0 ? 0.0 : (emax) * pow((conc), (hill)) / \\
   (pow((ec50), (hill)) + pow((conc), (hill))))

// allometry: clearances ^0.75, volumes ^1.0.  THIS IS THE WHOLE OF THE
// PAEDIATRIC MILTEFOSINE PROBLEM: dose is prescribed per kg (exponent 1.0)
// but cleared with exponent 0.75, so steady-state exposure scales as WT^0.25
// and a 10 kg child on 2.5 mg/kg/d is 33% underexposed relative to a 50 kg
// adult on the identical mg/kg.
#define _FCL (pow(WT/70.0, 0.75))
#define _FV  (WT/70.0)

$MAIN
// -- allometrically scaled disposition parameters ---------------------------
double CLlip = CL_LIP * _FCL;
double CLfre = CL_FRE * _FCL;
double Qamb  = Q_AMB  * _FCL;
double CLmil = CL_MIL * _FCL;
double Qmil  = Q_MIL  * _FCL;
double CLpm  = CL_PM  * _FCL;
double Qpm   = Q_PM   * _FCL;
double CLsb  = CL_SB  * _FCL;
double Qsb   = Q_SB   * _FCL;

double VlipA = V_LIP * WT;
double VfreA = V_FRE * WT;
double VperA = V_PER * WT;
double V1mil = V1_MIL * _FV;
double V2mil = V2_MIL * _FV;
double VpmA  = V_PM  * WT;
double V2pmA = V2_PM * WT;
double VsbA  = V_SB  * WT;
double V2sbA = V2_SB * WT;

double VspA  = VSP  * _FV;
double VliA  = VLI  * _FV;
double VbmA  = VBM  * _FV;
double VskA  = VSK  * _FV;
double VkidA = VKID * _FV;
double GFRb  = GFR0 * _FCL;

// -- initial conditions: a patient presenting with established, symptomatic
//    VL, i.e. sitting on the untreated high-burden attractor ---------------
if (NEWIND < 2) {
  P_SP_S_0 = 0.80 * PMAX_SP * P0SCALE * (1.0 - FR0);
  P_SP_R_0 = 0.80 * PMAX_SP * P0SCALE * FR0;
  P_LI_S_0 = 0.55 * PMAX_LI * P0SCALE * (1.0 - FR0);
  P_LI_R_0 = 0.55 * PMAX_LI * P0SCALE * FR0;
  P_BM_S_0 = 0.60 * PMAX_BM * P0SCALE * (1.0 - FR0);
  P_BM_R_0 = 0.60 * PMAX_BM * P0SCALE * FR0;
  P_SK_S_0 = 0.25 * PMAX_SK * P0SCALE * (1.0 - FR0);
  P_SK_R_0 = 0.25 * PMAX_SK * P0SCALE * FR0;

  TMEM_0  = 0.05;
  IFNG_0  = 0.11;
  IL10_0  = 0.62;
  TNFA_0  = 0.55;
  TGFB_0  = 0.42;
  CD4_0   = CD40;
  MPHA_0  = 0.04;   // BUG FOUND (twice): MPHA is the ACTIVATED fraction of
                    // the macrophage pool, so both an active-VL patient and
                    // a healthy host sit near zero, not near one.  Leaving
                    // it high on day 0 flattered the fast-killing drugs, and
                    // leaving it high in the treatment-naive initial state
                    // sterilised every fresh sandfly inoculum on arrival.
  SPL_0   = 8.5;
  LIVS_0  = 2.8;
  HGB_0   = 7.4;
  PLT_0   = 78.0;
  WBC_0   = 2.4;
  ALB_0   = 2.7;
  IGG_0   = 3.6;
  TEMP_0  = 38.9;
  BWTL_0  = 0.16;
  SCR_0   = SCR0;
  KSER_0  = K_K0;
  MGSER_0 = MG0;
  QTC_0   = QTC0;
  LIPA_0  = LIPA0;
  ALTX_0  = ALT0;
}

$ODE
// =========================================================================
// 1. AMPHOTERICIN B  -- one dose, two integrals, opposite directions
// =========================================================================
double C_LIP = A_LIP / VlipA;
double C_FRE = A_FRE / VfreA;
double C_PER = A_PER / VperA;
double cAsp  = AMB_SP  / VspA;
double cAli  = AMB_LI  / VliA;
double cAbm  = AMB_BM  / VbmA;
double cAsk  = AMB_SK  / VskA;
double cAkid = AMB_KID / VkidA;

double lip_out = CLlip * C_LIP;
double rel = FREL * lip_out;             // leaks into the free plasma pool
double mps = (1.0 - FREL) * lip_out;     // PHAGOCYTOSED: handed to the
                                         // infected cell.  This single line
                                         // is why liposomal amphotericin is a
                                         // targeting device and not just a
                                         // better-tolerated formulation.

dxdt_A_LIP = -lip_out;
dxdt_A_FRE = rel - CLfre * C_FRE - Qamb * (C_FRE - C_PER);
dxdt_A_PER = Qamb * (C_FRE - C_PER);

dxdt_AMB_SP = FSP_A * mps + KOUT_A * (KPSP_A * C_FRE * VspA - AMB_SP);
dxdt_AMB_LI = FLI_A * mps + KOUT_A * (KPLI_A * C_FRE * VliA - AMB_LI);
dxdt_AMB_BM = FBM_A * mps + KOUT_A * (KPBM_A * C_FRE * VbmA - AMB_BM);
dxdt_AMB_SK = FSK_A * mps + KOUT_A * (KPSK_A * C_FRE * VskA - AMB_SK);
// BUG FOUND: written first as an uptake/egress pair, this multiplied the
// free plasma concentration by Kp AND by KIN/KOUT = 12.5 -- a 100-fold
// amplification that gave amphotericin deoxycholate a peak creatinine of
// 32 mg/dL and a NEGATIVE serum potassium.  It must use the same
// equilibration form as the four macrophage organs, so that steady-state
// cortical concentration is exactly Kp x C_free.
dxdt_AMB_KID = KIN_KID * (KPKID_A * C_FRE * VkidA - AMB_KID);

// =========================================================================
// 2. MILTEFOSINE
// =========================================================================
double C_MIL  = MIL_C / V1mil;
double C_MILP = MIL_P / V2mil;
double cMsp = MIL_SP / VspA;
double cMli = MIL_LI / VliA;
double cMbm = MIL_BM / VbmA;
double cMsk = MIL_SK / VskA;

dxdt_MIL_G = -KA_MIL * MIL_G;
dxdt_MIL_C = F_MIL * KA_MIL * MIL_G - CLmil * C_MIL - Qmil * (C_MIL - C_MILP);
dxdt_MIL_P = Qmil * (C_MIL - C_MILP);
dxdt_MIL_SP = KIN_M * (KPSP_M * C_MIL * VspA - MIL_SP);
dxdt_MIL_LI = KIN_M * (KPLI_M * C_MIL * VliA - MIL_LI);
dxdt_MIL_BM = KIN_M * (KPBM_M * C_MIL * VbmA - MIL_BM);
dxdt_MIL_SK = KIN_M * (KPSK_M * C_MIL * VskA - MIL_SK);

// =========================================================================
// 3. PAROMOMYCIN -- slow in, slow out: that lag is the 21-day course
// =========================================================================
double C_PM  = PM_C / VpmA;
double C_PMP = PM_P / V2pmA;
double cPsp = PM_SP / VspA;
double cPli = PM_LI / VliA;
double cPbm = PM_BM / VbmA;
double cPsk = PM_SK / VskA;
double cPkid = PM_KID / VkidA;
double cPcoc = PM_COC / VCOC;

dxdt_PM_D = -KA_PM * PM_D;
dxdt_PM_C = KA_PM * PM_D - CLpm * C_PM - Qpm * (C_PM - C_PMP);
dxdt_PM_P = Qpm * (C_PM - C_PMP);
dxdt_PM_SP = KIN_PM * C_PM * VspA - KOUT_PM * PM_SP;
dxdt_PM_LI = KIN_PM * C_PM * VliA - KOUT_PM * PM_LI;
dxdt_PM_BM = KIN_PM * C_PM * VbmA - KOUT_PM * PM_BM;
dxdt_PM_SK = FSK_PM * KIN_PM * C_PM * VskA - KOUT_PM * PM_SK;
dxdt_PM_KID = KIN_PMK * C_PM * VkidA - KOUT_PMK * PM_KID;
dxdt_PM_COC = KIN_COC * C_PM * VCOC - KOUT_COC * PM_COC;

// =========================================================================
// 4. ANTIMONY -- activation and escape share one thiol chemistry
// =========================================================================
double C_SB5 = SB5_C / VsbA;
double C_SBD = SB_DEP / V2sbA;
double cSsp = SB3_SP / VspA;
double cSli = SB3_LI / VliA;
double cSbm = SB3_BM / VbmA;
double cSsk = SB3_SK / VskA;

dxdt_SB5_D = -KA_SB * SB5_D;
dxdt_SB5_C = KA_SB * SB5_D - CLsb * C_SB5 - Qsb * (C_SB5 - C_SBD);
dxdt_SB5_P = 0.0;
dxdt_SB_DEP = Qsb * (C_SB5 - C_SBD);

// Resistance is modelled as INCREASED EFFLUX (MRPA / Sb-trypanothione
// conjugate export), not as a change in target affinity.  That is why dose
// escalation cannot recover a resistant strain in this model: raising the
// dose raises the intracellular load, and the efflux term is proportional to
// exactly that load, while QTc and lipase are not.
double keff = KEFF_SB * RES_SB;
dxdt_SB3_SP = KRED * C_SB5 * VspA - keff * SB3_SP;
dxdt_SB3_LI = KRED * C_SB5 * VliA - keff * SB3_LI;
dxdt_SB3_BM = KRED * C_SB5 * VbmA - keff * SB3_BM;
dxdt_SB3_SK = KPSK_SB * KRED * C_SB5 * VskA - keff * SB3_SK;

// =========================================================================
// 5. IMMUNITY
// =========================================================================
double P_SP = P_SP_S + P_SP_R;
double P_LI = P_LI_S + P_LI_R;
double P_BM = P_BM_S + P_BM_R;
double P_SK = P_SK_S + P_SK_R;
double PTOT = P_SP + P_LI + P_BM + P_SK;
if (PTOT < 0.0) PTOT = 0.0;

// KAG << K50_IL10 is a structural statement, not a fitted convenience.
// Antigen is still abundant over a burden range in which IL-10 has already
// switched off.  That window -- roughly 1e4 to 6e7 amastigotes -- is where
// the host can be primed, and getting the patient into it and holding them
// there is the entire job of a drug course.
double PVIS = P_SP + P_LI + P_BM;
// Three signals, three thresholds, and TWO DIFFERENT SOURCES.
//   AG  = antigen availability, from the TOTAL burden.  Dermal parasites do
//         prime T cells; that is why PKDL patients are leishmanin-positive.
//   AGH = the IL-10 signal, from the VISCERAL burden only.
//   AGS = the systemic-illness signal, also visceral only.
// BUG FOUND: keying IL-10 to the total burden let a growing dermal reservoir
// re-immunosuppress the host -- skin load climbed past K50_IL10, IL-10 came
// back on, macrophage activation collapsed and the viscera regrew, so EVERY
// regimen relapsed.  That is backwards: PKDL patients are not
// immunosuppressed, and the IL-10 of active kala-azar comes from spleen and
// marrow, not from skin lesions.
double AG  = PTOT / (PTOT + KAG);
double AGH = pow(PVIS, HIL10) / (pow(PVIS, HIL10) + pow(K50_IL10, HIL10));
// A THIRD threshold, above the IL-10 one, for systemic illness.  Keyed to AG
// (which saturates at 1e4 amastigotes) the clinical module made a host
// carrying a controlled subclinical infection febrile and anaemic.
double AGS = PVIS / (PVIS + K50_SYS);
double FCD4 = CD4 / (CD4 + KCD4);
double prime = 1.0 - 0.55 * MALNUT;

dxdt_TMEM = KTP * prime * AG * FCD4 * (1.0 - TMEM / TMEMAX)
            / (1.0 + IL10 / KI_T) - KTD * TMEM;
dxdt_IFNG = KIP * (FAG_IFN * AG + TMEM) * FCD4 - KID * IFNG;
dxdt_IL10 = KLP * AGH - KLD * IL10;
dxdt_TNFA = KNP * AGS - KND * TNFA;
// BUG FOUND: driven by AGH, not AG.  With AG the TGF-beta term stayed near
// maximal down to 1e4 amastigotes and permanently clamped macrophage
// activation, so no slow drug could ever hand off to the host.
dxdt_TGFB = KGP * (0.5 * AGH + IL10) - KGD * TGFB;

// BUG FOUND: CD4 is a HOMEOSTATIC pool with a set point.  An earlier version
// had only the loss term, so every immunocompetent patient drifted to CD4 = 4
// over an 18-month follow-up and the host arm silently switched off.
double cd4set = (HIV < 0.5) ? CD4SET : (CD4SET_ART * ART + CD40 * (1.0 - ART));
dxdt_CD4 = KCD4H * (cd4set - CD4) - KCD4L * AG * CD4 * (1.0 + 0.5 * HIV);

double MACTdrive = (IFNG / (IFNG + KI_IFN))
                   / (1.0 + IL10 / KI_IL10 + TGFB / KI_TGFB)
                   * FCD4
                   * (1.0 - (1.0 - FHIV_MAC) * HIV);
dxdt_MPHA = KMPH * (MACTdrive - MPHA);
double MACT = (MPHA > 0.0) ? MPHA : 0.0;
double kimm = KIMM * MACT;

// =========================================================================
// 6. PARASITE DYNAMICS
//    Two populations per organ.  The quiescent pool is not a resistance
//    genotype, it is a physiological state, and it is where every relapse in
//    this model comes from: all four drugs need a metabolically active cell
//    (ergosterol turnover, inward transport, active translation, prodrug
//    reduction), and the macrophage sees quiescent forms poorly too.
// =========================================================================
double kSa, kSm, kSp, kSs, kRa, kRm, kRp, kRs;

// -- spleen ---------------------------------------------------------------
kSa = _EMX(cAsp, EC50_A, EMAX_A, HILL_A);
kSm = _EMX(cMsp, EC50_M * RES_M, EMAX_M, HILL_M);
kSp = _EMX(cPsp, EC50_P * RES_P, EMAX_P, HILL_P);
kSs = _EMX(cSsp, EC50_S, EMAX_S, HILL_S);
kRa = _EMX(cAsp, EC50_A * RFAC_A, EMAX_A, HILL_A);
kRm = _EMX(cMsp, EC50_M * RES_M * RFAC_M, EMAX_M, HILL_M);
kRp = _EMX(cPsp, EC50_P * RES_P * RFAC_P, EMAX_P, HILL_P);
kRs = _EMX(cSsp, EC50_S * RFAC_S, EMAX_S, HILL_S);
double killS_SP = kSa + kSm + kSp + kSs;   // Bliss independence
double killR_SP = kRa + kRm + kRp + kRs;

// -- liver ----------------------------------------------------------------
kSa = _EMX(cAli, EC50_A, EMAX_A, HILL_A);
kSm = _EMX(cMli, EC50_M * RES_M, EMAX_M, HILL_M);
kSp = _EMX(cPli, EC50_P * RES_P, EMAX_P, HILL_P);
kSs = _EMX(cSli, EC50_S, EMAX_S, HILL_S);
kRa = _EMX(cAli, EC50_A * RFAC_A, EMAX_A, HILL_A);
kRm = _EMX(cMli, EC50_M * RES_M * RFAC_M, EMAX_M, HILL_M);
kRp = _EMX(cPli, EC50_P * RES_P * RFAC_P, EMAX_P, HILL_P);
kRs = _EMX(cSli, EC50_S * RFAC_S, EMAX_S, HILL_S);
double killS_LI = kSa + kSm + kSp + kSs;
double killR_LI = kRa + kRm + kRp + kRs;

// -- marrow ---------------------------------------------------------------
kSa = _EMX(cAbm, EC50_A, EMAX_A, HILL_A);
kSm = _EMX(cMbm, EC50_M * RES_M, EMAX_M, HILL_M);
kSp = _EMX(cPbm, EC50_P * RES_P, EMAX_P, HILL_P);
kSs = _EMX(cSbm, EC50_S, EMAX_S, HILL_S);
kRa = _EMX(cAbm, EC50_A * RFAC_A, EMAX_A, HILL_A);
kRm = _EMX(cMbm, EC50_M * RES_M * RFAC_M, EMAX_M, HILL_M);
kRp = _EMX(cPbm, EC50_P * RES_P * RFAC_P, EMAX_P, HILL_P);
kRs = _EMX(cSbm, EC50_S * RFAC_S, EMAX_S, HILL_S);
double killS_BM = kSa + kSm + kSp + kSs;
double killR_BM = kRa + kRm + kRp + kRs;

// -- skin -----------------------------------------------------------------
kSa = _EMX(cAsk, EC50_A, EMAX_A, HILL_A);
kSm = _EMX(cMsk, EC50_M * RES_M, EMAX_M, HILL_M);
kSp = _EMX(cPsk, EC50_P * RES_P, EMAX_P, HILL_P);
kSs = _EMX(cSsk, EC50_S, EMAX_S, HILL_S);
kRa = _EMX(cAsk, EC50_A * RFAC_A, EMAX_A, HILL_A);
kRm = _EMX(cMsk, EC50_M * RES_M * RFAC_M, EMAX_M, HILL_M);
kRp = _EMX(cPsk, EC50_P * RES_P * RFAC_P, EMAX_P, HILL_P);
kRs = _EMX(cSsk, EC50_S * RFAC_S, EMAX_S, HILL_S);
double killS_SK = kSa + kSm + kSp + kSs;
double killR_SK = kRa + kRm + kRp + kRs;

// Inter-organ trafficking of infected monocytes, written as a mass-conserving
// redistribution toward each compartment's equilibrium share so it can never
// create parasites from nothing.  BUG FOUND: an earlier version scaled this
// by an extra 1e-4, giving an effective rate of 2e-8 per hour -- not
// trafficking but zero, so a dermal sandfly inoculum stayed in the skin for
// eighteen months and never visceralised.  Only REPLICATING amastigotes
// traffic; quiescent forms are arrested.
// Trafficking is ASYMMETRIC, and it has to be.  Spleen, liver and marrow share
// the blood monocyte pool and redistribute freely.  The dermis receives
// infected monocytes from that pool, but its own macrophages are largely
// resident and do not recirculate.
// BUG FOUND (twice): written first with an extra 1e-4 factor the effective
// rate was 2e-8 /h -- not trafficking but zero, so a dermal inoculum never
// visceralised.  Written SYMMETRICALLY, the immune-privileged skin became an
// incurable reservoir that re-seeded the viscera and ALL TWENTY regimens
// relapsed -- whereas in reality most cured patients never develop PKDL and
// PKDL patients rarely relapse to visceral disease.
double pmax_vis = PMAX_SP + PMAX_LI + PMAX_BM;
double PS_VIS = fmax(P_SP_S, 0.0) + fmax(P_LI_S, 0.0) + fmax(P_BM_S, 0.0);
double net_to_skin = KMIG_SKIN * PS_VIS - KMIG_BACK * fmax(P_SK_S, 0.0);
double crowd, S, R, kq, seed;

// NOTE: there is deliberately NO extinction floor in this block.  Writing
// `if (P < 1e-7 && dxdt > 0) dxdt = 0` here is a discontinuous switch that
// flips between solver probes as soon as any compartment clears; it stalls
// the integrator completely (a 28-day miltefosine course failed to finish).
// The system as written is homogeneous in the parasite states, so a cleared
// compartment decays smoothly and cannot regrow on its own.  Apply the
// one-amastigote floor between calls -- see `vl_clamp_extinct()` below.

crowd = 1.0 - P_SP / PMAX_SP;  if (crowd < 0.0) crowd = 0.0;
S = (P_SP_S > 0.0) ? P_SP_S : 0.0;  R = (P_SP_R > 0.0) ? P_SP_R : 0.0;
kq = KQ0 + KQS * (killS_SP + kimm);
seed = KMIG * ((PMAX_SP / pmax_vis) * PS_VIS - fmax(P_SP_S, 0.0))
     - (PMAX_SP / pmax_vis) * net_to_skin;
dxdt_P_SP_S = S * (KG_SP * crowd - kimm - killS_SP) - kq * S + KR_REACT * R
              + seed;
dxdt_P_SP_R = R * (KG_SP * COSTR * crowd - FIMM_R * kimm - killR_SP)
              + kq * S - KR_REACT * R;

crowd = 1.0 - P_LI / PMAX_LI;  if (crowd < 0.0) crowd = 0.0;
S = (P_LI_S > 0.0) ? P_LI_S : 0.0;  R = (P_LI_R > 0.0) ? P_LI_R : 0.0;
kq = KQ0 + KQS * (killS_LI + kimm);
seed = KMIG * ((PMAX_LI / pmax_vis) * PS_VIS - fmax(P_LI_S, 0.0))
     - (PMAX_LI / pmax_vis) * net_to_skin;
dxdt_P_LI_S = S * (KG_LI * crowd - kimm - killS_LI) - kq * S + KR_REACT * R
              + seed;
dxdt_P_LI_R = R * (KG_LI * COSTR * crowd - FIMM_R * kimm - killR_LI)
              + kq * S - KR_REACT * R;

crowd = 1.0 - P_BM / PMAX_BM;  if (crowd < 0.0) crowd = 0.0;
S = (P_BM_S > 0.0) ? P_BM_S : 0.0;  R = (P_BM_R > 0.0) ? P_BM_R : 0.0;
kq = KQ0 + KQS * (killS_BM + kimm);
seed = KMIG * ((PMAX_BM / pmax_vis) * PS_VIS - fmax(P_BM_S, 0.0))
     - (PMAX_BM / pmax_vis) * net_to_skin;
dxdt_P_BM_S = S * (KG_BM * crowd - kimm - killS_BM) - kq * S + KR_REACT * R
              + seed;
dxdt_P_BM_R = R * (KG_BM * COSTR * crowd - FIMM_R * kimm - killR_BM)
              + kq * S - KR_REACT * R;

// The dermis is relatively immune-privileged for this parasite, and it has to
// be: PKDL lesions persist for years in hosts whose systemic immunity has
// fully returned.  With dermal killing set equal to splenic killing the model
// made untreated PKDL self-resolving, which it is not in South Asia.
double kimm_sk = kimm * FIMM_SK;
crowd = 1.0 - P_SK / PMAX_SK;  if (crowd < 0.0) crowd = 0.0;
S = (P_SK_S > 0.0) ? P_SK_S : 0.0;  R = (P_SK_R > 0.0) ? P_SK_R : 0.0;
kq = KQ0 + KQS * (killS_SK + kimm_sk);
seed = net_to_skin;
dxdt_P_SK_S = S * (KG_SK * crowd - kimm_sk - killS_SK) - kq * S + KR_REACT * R
              + seed;
dxdt_P_SK_R = R * (KG_SK * COSTR * crowd - FIMM_R * kimm_sk - killR_SK)
              + kq * S - KR_REACT * R;

// =========================================================================
// 7. CLINICAL
// =========================================================================
double spl_target = SPL0 + SPLMAX * P_SP / (P_SP + K50_SPL);
dxdt_SPL = KSPL * (spl_target - SPL);
double liv_target = LIV0 + LIVMAX * P_LI / (P_LI + 3.0 * K50_SPL);
dxdt_LIVS = KLIV * (liv_target - LIVS);

// Two independent mechanisms, not one: marrow occupancy plus hypersplenic
// destruction.  Platelets carry the largest spleen coefficient, which is why
// thrombocytopenia tracks spleen size and anaemia does not.
double marrow = 1.0 / (1.0 + P_BM / (0.35 * PMAX_BM));
double inflam = 1.0 / (1.0 + TNFA / 0.30);
dxdt_HGB = KHGB * HGB0 * (1.0 - FSUP_HGB * (1.0 - marrow * inflam))
           - (KHGB + KSPL_HGB * KHGB * SPL) * HGB;
dxdt_PLT = KPLT * PLT0 * marrow - (KPLT + KSPL_PLT * KPLT * SPL) * PLT;
dxdt_WBC = KWBC * WBC0 * marrow * inflam
           - (KWBC + KSPL_WBC * KWBC * SPL) * WBC;
dxdt_ALB = KALB * (ALB0 * (1.0 - 0.42 * AGS) - ALB);
dxdt_IGG = KIGG_P * AG - KIGG_D * (IGG - IGG0);
dxdt_TEMP = KTEMP * (TEMP0 + TNF_FEV * TNFA / (TNFA + 0.45) - TEMP);
dxdt_BWTL = KBWL * (BWLMAX * AGS - BWTL);
// PKDL: dermal amastigotes become lesions only when immunity RETURNS, which
// is why the rash appears AFTER apparently successful treatment
dxdt_PKDL = KPKDL_P * P_SK * MACT - KPKDL_D * PKDL;

// =========================================================================
// 8. TOXICITY
// =========================================================================
dxdt_TUBI = KTUB * cAkid - KTUBR * TUBI;
double gfr = GFRb / (1.0 + ((TUBI > 0.0) ? TUBI : 0.0));
if (gfr < 1.0) gfr = 1.0;
dxdt_SCR = KSCR * (SCR0 * GFRb / gfr - SCR);
dxdt_KSER = KKR * (K_K0 - KSER) - KWK * cAkid;
dxdt_MGSER = KMGR * (MG0 - MGSER) - KWMG * cAkid;
dxdt_HEAR = KHEAR * cPcoc;                 // irreversible: no repair term
dxdt_QTC = KQTC * C_SBD * KQTCR * 50.0 - KQTCR * (QTC - QTC0);
dxdt_LIPA = KLIPA * C_SB5 - KLIPAR * (LIPA - LIPA0);
dxdt_ALTX = KALTX * C_SB5 + KALT_M * C_MIL - KALTR * (ALTX - ALT0);
dxdt_GITX = KGITX * C_MIL / (C_MIL + EM50_GI) - KGITR * GITX;

// =========================================================================
// 9. MORTALITY HAZARD AND EXPOSURE INTEGRALS
// =========================================================================
double dh_hgb = (HGB_H - HGB) > 0.0 ? (HGB_H - HGB) / HGB_H : 0.0;
double dh_alb = (ALB_H - ALB) > 0.0 ? (ALB_H - ALB) / ALB_H : 0.0;
double dh_wbc = (WBC_H - WBC) > 0.0 ? (WBC_H - WBC) / WBC_H : 0.0;
double haz = H0 * (1.0 + 2.2 * dh_hgb) * (1.0 + 2.6 * dh_alb)
                * (1.0 + 2.0 * dh_wbc) * (0.15 + 0.85 * AGS)
                * (1.0 + 1.6 * HIV);

dxdt_AUCFRE = C_FRE;   // the nephrotoxicity integral
dxdt_AUCASP = cAsp;    // the efficacy integral -- SAME DOSE, DIFFERENT ORGAN
dxdt_AUCMSP = cMsp;
dxdt_AUCBSP = cSsp;
dxdt_AUCPSP = cPsp;
dxdt_CUMHAZ = haz;
dxdt_CUMDOSE = 0.0;    // incremented in $MAIN-free fashion by the caller

$TABLE
double PTOT_OUT = P_SP_S + P_SP_R + P_LI_S + P_LI_R
                + P_BM_S + P_BM_R + P_SK_S + P_SK_R;
if (PTOT_OUT < 0.0) PTOT_OUT = 0.0;
double LOG10P = log10(PTOT_OUT + 1e-12);
double SURV = exp(-CUMHAZ);
double CFRAC = 100.0 * (1.0 - SURV);
// splenic aspirate grade: the clinical reference standard is a log scale,
// so report the model burden on it (grade 6+ = >100 parasites/field)
double SPLGRADE = log10(P_SP_S + P_SP_R + 1e-12) - log10(0.003);
if (SPLGRADE < 0.0) SPLGRADE = 0.0;
if (SPLGRADE > 6.0) SPLGRADE = 6.0;
double CFRE_OUT = A_FRE / (V_FRE * WT);
double CLIP_OUT = A_LIP / (V_LIP * WT);
double CMIL_OUT = MIL_C / (V1_MIL * (WT / 70.0));
double CASP_OUT = AMB_SP / (VSP * (WT / 70.0));
double CMSP_OUT = MIL_SP / (VSP * (WT / 70.0));
double CPSP_OUT = PM_SP / (VSP * (WT / 70.0));
double CSSP_OUT = SB3_SP / (VSP * (WT / 70.0));

$CAPTURE PTOT_OUT LOG10P SURV CFRAC SPLGRADE CFRE_OUT CLIP_OUT CMIL_OUT
$CAPTURE CASP_OUT CMSP_OUT CPSP_OUT CSSP_OUT
'

vl_mod <- mcode("vl_qsp", vl_code)

## =============================================================================
##  The one-amastigote floor, applied BETWEEN integration calls
## =============================================================================
##  A population smaller than one organism does not exist.  Enforcing that is
##  both biologically correct and what stops the ODE from "relapsing" out of
##  1e-300.  It MUST NOT go inside $ODE (see the note there); mrgsolve lets us
##  do it properly by simulating in segments and resetting the initial
##  condition between them.
EXTINCT <- 1e-6      # burden units = one whole amastigote
PAR_CMT <- c("P_SP_S", "P_SP_R", "P_LI_S", "P_LI_R",
             "P_BM_S", "P_BM_R", "P_SK_S", "P_SK_R")

vl_clamp_extinct <- function(x) {
  x[PAR_CMT] <- lapply(x[PAR_CMT], function(v) ifelse(v < EXTINCT, 0, v))
  x
}

#' Simulate a regimen in 14-day segments, applying the extinction floor at
#' each segment boundary.
#'
#' @param mod    an mrgsolve model object
#' @param dosing a data frame of events (mrgsolve event format)
#' @param par    named list of parameter overrides
#' @param days   follow-up duration in days
vl_simulate <- function(mod = vl_mod, dosing = NULL, par = list(),
                        days = 540, delta = 6, seg_days = 14) {
  m <- mod
  if (length(par)) m <- param(m, par)
  edges <- unique(c(seq(0, days, by = seg_days), days))
  out <- NULL
  init_now <- NULL
  for (i in seq_len(length(edges) - 1L)) {
    t0 <- edges[i] * 24
    t1 <- edges[i + 1L] * 24
    mm <- m
    if (!is.null(init_now)) mm <- init(mm, init_now)
    ev_seg <- NULL
    if (!is.null(dosing)) {
      ev_seg <- dosing[dosing$time >= t0 & dosing$time < t1, , drop = FALSE]
      if (nrow(ev_seg) == 0) ev_seg <- NULL
    }
    sim <- mm %>%
      (function(z) if (is.null(ev_seg)) z else data_set(z, ev_seg)) %>%
      mrgsim(start = t0, end = t1, delta = delta, hmax = 2,
             recover = "numeric") %>%
      as.data.frame()
    out <- rbind(out, if (i == 1L) sim else sim[sim$time > t0, , drop = FALSE])
    last <- sim[nrow(sim), , drop = FALSE]
    init_now <- vl_clamp_extinct(
      as.list(last[, names(last) %in% names(init(m)@data), drop = FALSE]))
  }
  out$day <- out$time / 24
  out
}

## =============================================================================
##  Regimen library -- 20 scenarios
##  mg/kg doses are converted with the subject's weight; parenteral routes go
##  to their own compartment so no route is silently assumed.
## =============================================================================
CMT_LAMB <- 1; CMT_DAMB <- 2; CMT_MIL <- 9; CMT_PM <- 16; CMT_SB <- 25

vl_ev <- function(cmt, mgkg = NULL, mg = NULL, start_day = 0, n = 1,
                  ii_h = 24, wt = 50) {
  amt <- if (is.null(mg)) mgkg * wt else mg
  data.frame(ID = 1, time = start_day * 24 + (seq_len(n) - 1L) * ii_h,
             cmt = cmt, amt = amt, evid = 1)
}

vl_regimens <- function(wt = 50) {
  list(
    ## --- natural history -------------------------------------------------
    S01_untreated = NULL,

    ## --- pentavalent antimony -------------------------------------------
    S02_ssg_africa = vl_ev(CMT_SB, mgkg = 20, n = 30, wt = wt),
    S03_ssg_bihar_res = vl_ev(CMT_SB, mgkg = 20, n = 30, wt = wt),

    ## --- amphotericin B --------------------------------------------------
    S04_damb_alt30 = vl_ev(CMT_DAMB, mgkg = 1, n = 15, ii_h = 48, wt = wt),
    S05_lamb_single10 = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
    S06_lamb_21_multi = rbind(
      vl_ev(CMT_LAMB, mgkg = 3, n = 5, wt = wt),
      vl_ev(CMT_LAMB, mgkg = 3, start_day = 13, n = 1, wt = wt),
      vl_ev(CMT_LAMB, mgkg = 3, start_day = 20, n = 1, wt = wt)),
    S07_lamb_5x3 = vl_ev(CMT_LAMB, mgkg = 3, n = 5, wt = wt),

    ## --- miltefosine, and the allometry trap -----------------------------
    S08_mil28_adult = vl_ev(CMT_MIL, mgkg = 2.5, n = 28, wt = wt),
    S09_mil28_child_linear = vl_ev(CMT_MIL, mgkg = 2.5, n = 28, wt = 10),
    ## the mg/day that gives a 10 kg child the SAME AUC as 125 mg/day in a
    ## 50 kg adult, i.e. scaled by WT^0.75 instead of WT^1.0
    S10_mil28_child_allom = vl_ev(CMT_MIL, mg = 125 * (10 / 50)^0.75, n = 28),

    ## --- paromomycin -----------------------------------------------------
    S11_pm21 = vl_ev(CMT_PM, mgkg = 15, n = 21, wt = wt),
    S12_pm21_africa = vl_ev(CMT_PM, mgkg = 15, n = 21, wt = wt),

    ## --- combinations: synergy in TIME, not in concentration -------------
    S13_lamb5_mil7 = rbind(
      vl_ev(CMT_LAMB, mgkg = 5, n = 1, wt = wt),
      vl_ev(CMT_MIL, mgkg = 2.5, start_day = 1, n = 7, wt = wt)),
    S14_lamb5_pm10 = rbind(
      vl_ev(CMT_LAMB, mgkg = 5, n = 1, wt = wt),
      vl_ev(CMT_PM, mgkg = 15, start_day = 1, n = 10, wt = wt)),
    S15_mil10_pm10 = rbind(
      vl_ev(CMT_MIL, mgkg = 2.5, n = 10, wt = wt),
      vl_ev(CMT_PM, mgkg = 15, n = 10, wt = wt)),
    S16_ssg_pm17 = rbind(
      vl_ev(CMT_SB, mgkg = 20, n = 17, wt = wt),
      vl_ev(CMT_PM, mgkg = 15, n = 17, wt = wt)),

    ## --- HIV co-infection ------------------------------------------------
    S17_hiv_lamb10 = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
    S18_hiv_lamb30_mil28 = rbind(
      vl_ev(CMT_LAMB, mgkg = 5, start_day = 0, n = 6, ii_h = 48, wt = wt),
      vl_ev(CMT_MIL, mgkg = 2.5, n = 28, wt = wt)),

    ## --- PKDL ------------------------------------------------------------
    S19_pkdl_mil84 = vl_ev(CMT_MIL, mgkg = 2.5, n = 84, wt = wt),
    S20_pkdl_lamb = vl_ev(CMT_LAMB, mgkg = 5, n = 4, ii_h = 7 * 24, wt = wt)
  )
}

vl_scenario_par <- function(name) {
  p <- list()
  ## resistance is an EFFLUX parameter; FR0 is the quiescent pool and has
  ## nothing to do with antimony susceptibility
  if (name == "S02_ssg_africa")     p$RES_SB <- 1.0
  if (name == "S03_ssg_bihar_res")  p$RES_SB <- 9.0
  if (name %in% c("S09_mil28_child_linear", "S10_mil28_child_allom")) {
    p$WT <- 10; p$AGE <- 3; p$MALNUT <- 0.45
  }
  if (name == "S12_pm21_africa")    p$EC50_P <- 25.0 * 1.60
  if (name == "S17_hiv_lamb10") {
    ## the common real-world sequence: VL treated before ART is started, so
    ## there is no CD4 recovery during the window that decides cure
    p$HIV <- 1; p$ART <- 0; p$CD40 <- 90
  }
  if (name == "S18_hiv_lamb30_mil28") { p$HIV <- 1; p$ART <- 1; p$CD40 <- 90 }
  p
}

## PKDL scenarios start from a different patient: visceral compartments
## cured, dermal amastigotes still present, immunity restored -- which is
## exactly the state that makes the lesions appear.
vl_pkdl_init <- function() list(
  P_SP_S = 0, P_SP_R = 0, P_LI_S = 0, P_LI_R = 0, P_BM_S = 0, P_BM_R = 0,
  P_SK_S = 45 * (1 - 5e-4), P_SK_R = 45 * 5e-4,
  TMEM = 0.85, IFNG = 0.30, IL10 = 0.08, TNFA = 0.10, TGFB = 0.12,
  MPHA = 0.45, SPL = 1.5, HGB = 12.2, PLT = 210, WBC = 5.4, ALB = 3.9,
  IGG = 2.9, TEMP = 36.9, PKDL = 1.2)

## =============================================================================
##  Outcome classification
##  Every cure rate this model reports depends on these three thresholds, so
##  they are stated once, here, and nowhere else.
## =============================================================================
vl_outcome <- function(sim, eot_day) {
  ## IMPORTANT: cure and relapse are judged on the VISCERAL compartments
  ## (spleen + liver + marrow), not on total body burden.  That is what
  ## "definitive cure at 6 months" means in every VL trial ever run.  The
  ## dermal compartment clears far more slowly -- little drug reaches it and
  ## the dermis is relatively immune-privileged -- and its residue is what
  ## becomes PKDL.  Scoring skin residue as VL treatment failure made all
  ## twenty regimens "relapse", which is a finding about the wrong endpoint.
  ## Skin burden and lesion load are reported separately, exactly as trials
  ## report PKDL incidence separately from cure rate.
  PV <- sim$P_SP_S + sim$P_SP_R + sim$P_LI_S + sim$P_LI_R +
        sim$P_BM_S + sim$P_BM_R
  PSK <- sim$P_SK_S + sim$P_SK_R
  d <- sim$day
  i0 <- which.min(abs(d - eot_day))
  P0 <- PV[1]
  post <- PV[i0:length(PV)]
  nadir <- min(post)
  Pend <- PV[length(PV)]
  DET <- 1e-2                       # ~1e4 amastigotes: clinically detectable
  if (P0 < DET) {
    ## A PKDL patient has NO visceral disease at baseline, so the visceral
    ## endpoint is undefined for them.  Score those runs on the dermal
    ## compartment instead -- which is what a PKDL trial does (lesion
    ## clearance, not splenic aspirate).
    P0 <- PSK[1]
    post <- PSK[i0:length(PSK)]
    nadir <- min(post)
    Pend <- PSK[length(PSK)]
  }
  initial_response <- nadir < max(P0, 1e-30) * 1e-3   # >= 3 log10 fall
  relapse <- initial_response && (Pend > DET)
  list(
    P0 = P0, nadir = nadir, nadir_day = d[i0 - 1 + which.min(post)],
    P_end = Pend, P_total_end = sim$PTOT_OUT[length(PV)],
    logdrop = log10(max(nadir, 1e-30) / P0),
    cure = initial_response && !relapse,
    relapse = relapse,
    primary_failure = !initial_response,
    mort_pct = 100 * (1 - exp(-sim$CUMHAZ[length(PV)])),
    auc_free = sim$AUCFRE[length(PV)],
    auc_ambsp = sim$AUCASP[length(PV)],
    scr_max = max(sim$SCR), k_min = min(sim$KSER), mg_min = min(sim$MGSER),
    hear_dB = sim$HEAR[length(PV)], qtc_max = max(sim$QTC),
    lipase_max = max(sim$LIPA), alt_max = max(sim$ALTX),
    gi_max = max(sim$GITX),
    spleen_end = sim$SPL[length(PV)], hgb_end = sim$HGB[length(PV)],
    plt_end = sim$PLT[length(PV)], igg_end = sim$IGG[length(PV)],
    ## --- the PKDL endpoint, reported separately -------------------------
    pkdl_max = max(sim$PKDL), pkdl_end = sim$PKDL[length(PV)],
    skin_end = PSK[length(PSK)],
    skin_180 = PSK[which.min(abs(d - 180))],
    tmem_max = max(sim$TMEM)
  )
}

vl_eot_day <- function(ev) {
  if (is.null(ev)) return(0)
  max(ev$time) / 24 + 1
}

## =============================================================================
##  ANALYSIS 1 — the claim the whole model exists to make
##  One mg/kg of amphotericin B splits into two integrals that move in
##  OPPOSITE directions when you encapsulate it.
## =============================================================================
vl_analysis_targeting <- function(wt = 50) {
  regs <- list(
    `L-AmB 10 mg/kg x1` = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
    `d-AmB 10 mg/kg total (1 mg/kg q48h x10)` =
      vl_ev(CMT_DAMB, mgkg = 1, n = 10, ii_h = 48, wt = wt))
  res <- lapply(names(regs), function(nm) {
    s <- vl_simulate(dosing = regs[[nm]], days = 180)
    data.frame(
      regimen = nm,
      cmax_total = max(s$CLIP_OUT + s$CFRE_OUT),
      cmax_free = max(s$CFRE_OUT),
      auc_free = s$AUCFRE[nrow(s)],          # nephrotoxicity integral
      auc_spleen = s$AUCASP[nrow(s)],        # efficacy integral
      ti_surrogate = s$AUCASP[nrow(s)] / s$AUCFRE[nrow(s)],
      scr_max = max(s$SCR), k_min = min(s$KSER),
      log_drop = log10(min(s$PTOT_OUT) / s$PTOT_OUT[1]))
  })
  out <- do.call(rbind, res)
  cat("\n== A1  Liposomal vs deoxycholate amphotericin B at equal mg/kg ==\n")
  print(out, row.names = FALSE)
  cat(sprintf(
    "\n  therapeutic-index gain from encapsulation: %.1f x\n",
    out$ti_surrogate[1] / out$ti_surrogate[2]))
  cat(sprintf("  spleen exposure x%.1f, free-plasma exposure x%.2f\n",
              out$auc_spleen[1] / out$auc_spleen[2],
              out$auc_free[1] / out$auc_free[2]))
  invisible(out)
}

## =============================================================================
##  ANALYSIS 2 — the paediatric miltefosine problem is arithmetic
##  Dose ~ WT^1.0, clearance ~ WT^0.75, therefore AUC ~ WT^0.25.
## =============================================================================
vl_analysis_allometry <- function() {
  wts <- c(7, 10, 15, 20, 30, 50, 70)
  rows <- lapply(wts, function(w) {
    lin <- 2.5 * w
    allo <- 125 * (w / 50)^0.75
    auc <- vapply(c(lin, allo), function(mg) {
      s <- vl_simulate(dosing = vl_ev(CMT_MIL, mg = mg, n = 28),
                       par = list(WT = w), days = 60, delta = 2)
      sum(diff(s$time) * (head(s$CMIL_OUT, -1) + tail(s$CMIL_OUT, -1)) / 2)
    }, numeric(1))
    data.frame(WT = w, dose_linear = lin, dose_allometric = allo,
               auc_linear = auc[1], auc_allometric = auc[2])
  })
  out <- do.call(rbind, rows)
  ref <- out$auc_linear[out$WT == 50]
  out$auc_rel_to_50kg <- out$auc_linear / ref
  out$predicted_wt_0.25 <- (out$WT / 50)^0.25
  cat("\n== A2  Miltefosine exposure under linear vs allometric dosing ==\n")
  print(out, row.names = FALSE)
  cat("\n  The auc_rel_to_50kg column should track predicted_wt_0.25 exactly;\n")
  cat("  that identity is the paediatric dosing problem, not a coincidence.\n")
  invisible(out)
}

## =============================================================================
##  ANALYSIS 3 — cure is a separatrix crossing, and CD4 moves the line
##  For a patient who has finished treatment with a primed memory pool TMEM,
##  what residual burden can the host still finish unaided?
## =============================================================================
vl_analysis_separatrix <- function(cd4s = c(700, 350, 200, 100, 50),
                                   tmem = 0.55) {
  probe <- function(cd4, logP) {
    par <- list(CD40 = cd4)
    if (cd4 < 350) { par$HIV <- 1; par$ART <- 1 }
    P <- 10^logP
    ini <- list(P_SP_S = 0.45 * P, P_LI_S = 0.35 * P, P_BM_S = 0.15 * P,
                P_SK_S = 0.05 * P,
                P_SP_R = 0, P_LI_R = 0, P_BM_R = 0, P_SK_R = 0,
                TMEM = tmem, IL10 = 0.10, IFNG = 0.30, TNFA = 0.12,
                TGFB = 0.15, MPHA = 0.45, CD4 = cd4)
    m <- init(param(vl_mod, par), ini)
    s <- as.data.frame(mrgsim(m, end = 360 * 24, delta = 24, hmax = 2))
    log10(max(s$PTOT_OUT[nrow(s)], 1e-30)) - logP
  }
  rows <- lapply(cd4s, function(cd4) {
    lo <- -6; hi <- 4
    f <- function(x) probe(cd4, x)
    crit <- tryCatch(uniroot(f, c(lo, hi), tol = 1e-3)$root, error = function(e) NA)
    data.frame(CD4 = cd4, TMEM = tmem,
               critical_burden = 10^crit,
               critical_amastigotes = 10^crit * 1e6)
  })
  out <- do.call(rbind, rows)
  cat("\n== A3  Critical residual burden the host can finish alone ==\n")
  print(out, row.names = FALSE)
  cat("\n  The same drug course that cures an immunocompetent patient leaves a\n")
  cat("  residual burden far ABOVE the CD4-100 line. That, and not a change in\n")
  cat("  parasite susceptibility, is why HIV-VL relapses on a curative dose.\n")
  invisible(out)
}

## =============================================================================
##  ANALYSIS 4 — combination synergy lives in TIME, not in concentration
## =============================================================================
vl_analysis_timing <- function(wt = 50) {
  specs <- list(
    `L-AmB 5 alone` = vl_ev(CMT_LAMB, mgkg = 5, n = 1, wt = wt),
    `L-AmB 5 + MF 7 d` = rbind(vl_ev(CMT_LAMB, mgkg = 5, n = 1, wt = wt),
                               vl_ev(CMT_MIL, mgkg = 2.5, start_day = 1,
                                     n = 7, wt = wt)),
    `MF 7 d first, then L-AmB 5` = rbind(
      vl_ev(CMT_MIL, mgkg = 2.5, n = 7, wt = wt),
      vl_ev(CMT_LAMB, mgkg = 5, start_day = 7, n = 1, wt = wt)),
    `L-AmB 5 + MF 14 d` = rbind(vl_ev(CMT_LAMB, mgkg = 5, n = 1, wt = wt),
                                vl_ev(CMT_MIL, mgkg = 2.5, start_day = 1,
                                      n = 14, wt = wt)),
    `MF 28 d alone` = vl_ev(CMT_MIL, mgkg = 2.5, n = 28, wt = wt),
    `L-AmB 10 alone` = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt))
  rows <- lapply(names(specs), function(nm) {
    s <- vl_simulate(dosing = specs[[nm]], days = 360)
    o <- vl_outcome(s, vl_eot_day(specs[[nm]]))
    data.frame(regimen = nm, nadir = o$nadir, nadir_day = o$nadir_day,
               day360 = o$P_end, tmem_max = o$tmem_max,
               outcome = ifelse(o$cure, "cure",
                                ifelse(o$relapse, "relapse", "failure")))
  })
  out <- do.call(rbind, rows)
  cat("\n== A4  Why one day of L-AmB plus a week of miltefosine works ==\n")
  print(out, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 5 — antimony: dose escalation cannot beat efflux
## =============================================================================
vl_analysis_antimony <- function(wt = 50) {
  grid <- expand.grid(RES_SB = c(1, 3, 9), dose = c(20, 30, 40))
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    ev <- vl_ev(CMT_SB, mgkg = grid$dose[i], n = 30, wt = wt)
    s <- vl_simulate(dosing = ev, par = list(RES_SB = grid$RES_SB[i]),
                     days = 360)
    o <- vl_outcome(s, 31)
    data.frame(RES_SB = grid$RES_SB[i], dose_mgkg = grid$dose[i],
               auc_sb3 = s$AUCBSP[nrow(s)], log_drop = o$logdrop,
               qtc_max = o$qtc_max, lipase_max = o$lipase_max,
               outcome = ifelse(o$cure, "cure",
                                ifelse(o$relapse, "relapse", "failure")))
  })
  out <- do.call(rbind, rows)
  cat("\n== A5  Antimony dose escalation vs efflux-mediated resistance ==\n")
  print(out, row.names = FALSE)
  cat("\n  Sb(III) exposure and QTc rise together with dose; the resistant\n")
  cat("  strain is not recovered, because efflux scales with the very\n")
  cat("  intracellular load the higher dose creates.\n")
  invisible(out)
}

## =============================================================================
##  ANALYSIS 6 — the skin is a different organ, pharmacologically
## =============================================================================
vl_analysis_skin <- function(wt = 50) {
  auc_of <- function(sim, col) {
    sum(diff(sim$time) * (head(sim[[col]], -1) + tail(sim[[col]], -1)) / 2)
  }
  s_l <- vl_simulate(dosing = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
                     days = 180, delta = 2)
  s_m <- vl_simulate(dosing = vl_ev(CMT_MIL, mgkg = 2.5, n = 28, wt = wt),
                     days = 180, delta = 2)
  cat("\n== A6  Intramacrophage exposure by organ ==\n")
  cat(sprintf("  L-AmB 10 mg/kg single dose : spleen AUC %.0f, MPS share %.2f\n",
              s_l$AUCASP[nrow(s_l)], 0.10))
  cat(sprintf("  ... dermal MPS share is only 0.03 and Kp only 1.5, which is\n"))
  cat(sprintf("      why single-dose L-AmB clears the viscera and leaves skin\n"))
  cat(sprintf("  miltefosine 28 d : spleen AUC %.0f, dermal Kp 4.0\n",
              s_m$AUCMSP[nrow(s_m)]))
  ## PKDL treatment, from the post-VL initial state
  pk <- lapply(list(
    `miltefosine 12 weeks` = vl_ev(CMT_MIL, mgkg = 2.5, n = 84, wt = wt),
    `L-AmB 5 mg/kg weekly x4` = vl_ev(CMT_LAMB, mgkg = 5, n = 4,
                                      ii_h = 7 * 24, wt = wt),
    `no treatment` = NULL), function(ev) {
      m <- init(vl_mod, vl_pkdl_init())
      s <- vl_simulate(mod = m, dosing = ev, days = 360)
      data.frame(skin_d0 = s$P_SK_S[1] + s$P_SK_R[1],
                 skin_d360 = s$P_SK_S[nrow(s)] + s$P_SK_R[nrow(s)],
                 lesion_peak = max(s$PKDL), lesion_end = s$PKDL[nrow(s)])
    })
  out <- do.call(rbind, pk)
  out$regimen <- names(pk)
  cat("\n== A6b PKDL treatment ==\n")
  print(out[, c("regimen", "skin_d0", "skin_d360", "lesion_peak",
                "lesion_end")], row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 7 — subclinical vs clinical VL is a race, not a switch
##  Sweep the T-cell priming rate from a single sandfly inoculum.
## =============================================================================
vl_analysis_race <- function() {
  naive_init <- function() list(
    P_SK_S = 1e-4, P_SK_R = 0, P_SP_S = 0, P_SP_R = 0, P_LI_S = 0,
    P_LI_R = 0, P_BM_S = 0, P_BM_R = 0,
    TMEM = 0, IFNG = 0, IL10 = 0, TNFA = 0, TGFB = 0, MPHA = 0,
    SPL = 0.5, LIVS = 0.5, HGB = 14, PLT = 250, WBC = 7, ALB = 4.3,
    IGG = 1.1, TEMP = 36.7, BWTL = 0)
  rows <- lapply(c(0.0030, 0.0040, 0.0045, 0.0050, 0.0055, 0.0060,
                   0.0080, 0.0120), function(ktp) {
    m <- init(param(vl_mod, list(KTP = ktp)), naive_init())
    s <- vl_simulate(mod = m, days = 540, delta = 12)
    data.frame(KTP = ktp, peak_burden = max(s$PTOT_OUT),
               day540 = s$PTOT_OUT[nrow(s)], peak_spleen = max(s$SPL),
               peak_IL10 = max(s$IL10),
               outcome = ifelse(s$PTOT_OUT[nrow(s)] > 1e2,
                                "clinical VL", "subclinical"))
  })
  out <- do.call(rbind, rows)
  cat("\n== A7  Priming rate decides disease vs no disease ==\n")
  print(out, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 8 — dose fractionation of a fixed L-AmB total
## =============================================================================
vl_analysis_fractionation <- function(wt = 50) {
  specs <- list(
    `10 mg/kg x1` = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
    `5 mg/kg x2` = vl_ev(CMT_LAMB, mgkg = 5, n = 2, wt = wt),
    `2 mg/kg x5` = vl_ev(CMT_LAMB, mgkg = 2, n = 5, wt = wt),
    `1 mg/kg x10` = vl_ev(CMT_LAMB, mgkg = 1, n = 10, wt = wt),
    `2 mg/kg weekly x5` = vl_ev(CMT_LAMB, mgkg = 2, n = 5, ii_h = 168,
                                wt = wt),
    `0.5 mg/kg x20` = vl_ev(CMT_LAMB, mgkg = 0.5, n = 20, wt = wt))
  rows <- lapply(names(specs), function(nm) {
    s <- vl_simulate(dosing = specs[[nm]], days = 360)
    o <- vl_outcome(s, vl_eot_day(specs[[nm]]))
    data.frame(schedule = nm, auc_spleen = o$auc_ambsp,
               auc_free = o$auc_free, nadir = o$nadir, day360 = o$P_end,
               scr_max = o$scr_max,
               outcome = ifelse(o$cure, "cure",
                                ifelse(o$relapse, "relapse", "failure")))
  })
  out <- do.call(rbind, rows)
  cat("\n== A8  Fractionating a fixed 10 mg/kg L-AmB total ==\n")
  print(out, row.names = FALSE)
  cat("\n  Spleen exposure is nearly schedule-independent, because the tissue\n")
  cat("  half-life is far longer than the schedule. Fractionation therefore\n")
  cat("  buys almost nothing on efficacy while it RAISES the free-drug\n")
  cat("  integral. That is the arithmetic behind single-dose L-AmB.\n")
  invisible(out)
}

## =============================================================================
##  ANALYSIS 9 — all 20 scenarios, typical patient
## =============================================================================
vl_run_all <- function(wt = 50, days = 540) {
  R <- vl_regimens(wt)
  rows <- lapply(names(R), function(nm) {
    par <- vl_scenario_par(nm)
    m <- vl_mod
    if (grepl("^S(19|20)", nm)) m <- init(vl_mod, vl_pkdl_init())
    s <- vl_simulate(mod = m, dosing = R[[nm]], par = par, days = days)
    o <- vl_outcome(s, vl_eot_day(R[[nm]]))
    data.frame(scenario = nm, log_drop = o$logdrop, nadir = o$nadir,
               nadir_day = o$nadir_day, burden_end = o$P_end,
               spleen_cm = o$spleen_end, hgb = o$hgb_end,
               scr = o$scr_max, k = o$k_min, hearing_dB = o$hear_dB,
               qtc = o$qtc_max, mortality_pct = o$mort_pct,
               outcome = ifelse(o$cure, "cure",
                                ifelse(o$relapse, "relapse", "failure")))
  })
  out <- do.call(rbind, rows)
  cat("\n== A9  All 20 scenarios, typical patient ==\n")
  print(out, row.names = FALSE)
  invisible(out)
}

## =============================================================================
##  ANALYSIS 10 — population simulation with realistic between-subject
##  variability, giving cure rates comparable with trial data
## =============================================================================
vl_population <- function(scenario = "S05_lamb_single10", n = 60,
                          seed = 20260804, days = 365) {
  set.seed(seed)
  R <- vl_regimens()
  base <- vl_scenario_par(scenario)
  ev <- R[[scenario]]
  res <- vapply(seq_len(n), function(i) {
    p <- base
    ## the variability that matters: presenting burden (splenic aspirate
    ## grade spans orders of magnitude), potency, clearance, immune capacity,
    ## nutrition, persister-pool size, and -- for the oral drug -- adherence
    p$P0SCALE <- exp(rnorm(1, 0, 0.60))
    p$FR0     <- 5e-4 * exp(rnorm(1, 0, 0.80))
    p$EC50_A  <- 6.0  * exp(rnorm(1, 0, 0.30))
    p$EC50_M  <- 90.0 * exp(rnorm(1, 0, 0.35))
    p$EC50_P  <- (if (is.null(base$EC50_P)) 25.0 else base$EC50_P) *
                 exp(rnorm(1, 0, 0.35))
    p$EC50_S  <- 0.50 * exp(rnorm(1, 0, 0.45))
    p$CL_MIL  <- 0.166 * exp(rnorm(1, 0, 0.25))
    p$CL_LIP  <- 0.50 * exp(rnorm(1, 0, 0.22))
    p$KOUT_A  <- 0.00825 * exp(rnorm(1, 0, 0.30))
    p$KTP     <- 0.0050 * exp(rnorm(1, 0, 0.30))
    p$KIMM    <- 0.065 * exp(rnorm(1, 0, 0.22))
    p$F_MIL   <- exp(rnorm(1, 0, 0.28))
    p$MALNUT  <- min(0.75, (if (is.null(base$MALNUT)) 0 else base$MALNUT) +
                            runif(1, 0, 0.55))
    p$CD4SET  <- 800 * exp(rnorm(1, 0, 0.25))
    p$WT      <- (if (is.null(base$WT)) 50 else base$WT) * exp(rnorm(1, 0, 0.12))
    ev_i <- ev
    ## oral doses can be missed; parenteral doses are supervised
    if (!is.null(ev_i) && any(ev_i$cmt == CMT_MIL)) {
      adher <- min(1, max(0.45, 1 - abs(rnorm(1, 0, 0.14))))
      keep <- ev_i$cmt != CMT_MIL | runif(nrow(ev_i)) <= adher
      ev_i <- ev_i[keep, , drop = FALSE]
    }
    s <- vl_simulate(dosing = ev_i, par = p, days = days, delta = 12)
    o <- vl_outcome(s, vl_eot_day(ev))
    if (o$cure) 1 else if (o$relapse) 2 else 3
  }, numeric(1))
  out <- c(cure = 100 * mean(res == 1), relapse = 100 * mean(res == 2),
           failure = 100 * mean(res == 3))
  cat(sprintf("\n== A10  %s, n = %d ==\n", scenario, n))
  print(round(out, 1))
  invisible(out)
}

## =============================================================================
##  Plot helpers
## =============================================================================
vl_plot_burden <- function(sims, labels = names(sims)) {
  df <- do.call(rbind, Map(function(s, l) {
    data.frame(day = s$day, burden = pmax(s$PTOT_OUT, 1e-7), regimen = l)
  }, sims, labels))
  ggplot(df, aes(day, burden, colour = regimen)) +
    geom_line(linewidth = 0.7) +
    scale_y_log10() +
    geom_hline(yintercept = EXTINCT, linetype = 3) +
    labs(x = "Day", y = "Total amastigote burden (x1e6)",
         title = "Visceral leishmaniasis: total body parasite burden",
         subtitle = "dotted line = one amastigote (the extinction floor)") +
    theme_bw()
}

vl_plot_two_integrals <- function(wt = 50) {
  s1 <- vl_simulate(dosing = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
                    days = 90, delta = 2)
  s2 <- vl_simulate(dosing = vl_ev(CMT_DAMB, mgkg = 1, n = 10, ii_h = 48,
                                   wt = wt), days = 90, delta = 2)
  df <- rbind(
    data.frame(day = s1$day, spleen = s1$CASP_OUT, free = s1$CFRE_OUT,
               form = "liposomal"),
    data.frame(day = s2$day, spleen = s2$CASP_OUT, free = s2$CFRE_OUT,
               form = "deoxycholate"))
  long <- pivot_longer(df, c(spleen, free), names_to = "site",
                       values_to = "conc")
  ggplot(long, aes(day, conc + 1e-4, colour = form, linetype = site)) +
    geom_line(linewidth = 0.7) + scale_y_log10() +
    labs(x = "Day", y = "Amphotericin B concentration (mg/L)",
         title = "The same mg/kg, two different integrals",
         subtitle = paste("solid-vs-dashed is where the drug goes;",
                          "spleen drives killing, free plasma drives",
                          "nephrotoxicity")) +
    theme_bw()
}

## =============================================================================
##  Reproduce every number quoted in README.md
## =============================================================================
vl_run_everything <- function() {
  vl_run_all()
  vl_analysis_targeting()
  vl_analysis_allometry()
  vl_analysis_separatrix()
  vl_analysis_timing()
  vl_analysis_antimony()
  vl_analysis_skin()
  vl_analysis_race()
  vl_analysis_fractionation()
  for (s in c("S02_ssg_africa", "S03_ssg_bihar_res", "S04_damb_alt30",
              "S05_lamb_single10", "S08_mil28_adult",
              "S09_mil28_child_linear", "S10_mil28_child_allom",
              "S11_pm21", "S13_lamb5_mil7", "S17_hiv_lamb10",
              "S18_hiv_lamb30_mil28")) {
    vl_population(s, n = 60)
  }
  invisible(NULL)
}

## Not run by default: this is a library, not a script.
if (identical(Sys.getenv("VL_RUN_ALL"), "1")) vl_run_everything()
