# ============================================================================
# VEXAS Syndrome QSP Model — mrgsolve specification
#
# Disease: VEXAS (Vacuoles, E1 enzyme, X-linked, Autoinflammatory, Somatic)
# Mechanism: somatic UBA1 (codon Met41) mutation in HSC → loss of cytosolic
#            UBA1b → defective cytoplasmic ubiquitination → ER stress, UPR,
#            mitochondrial dysfunction, NF-κB / NLRP3 / type-I IFN
#            hyperactivation → autoinflammation + myelodysplastic cytopenias.
#
# 23 ODE compartments cover:
#   • HSC clone dynamics + VAF (UBA1 mutant fraction)
#   • Misfolded protein burden / ER stress / mitochondrial ROS
#   • Core cytokines: IL-6, IL-1β, TNF-α, IFN-α, CXCL8, CCL2
#   • Acute-phase: CRP, ferritin
#   • Hematology: Hb (macrocytic anemia), Platelets, Neutrophils
#   • Clinical: fever index, skin/chondritis activity, VTE risk
#   • Steroid HPA suppression
#   • Drug PK: prednisone, tocilizumab (TMDD), anakinra, canakinumab,
#             ruxolitinib (1-cpt PK), azacitidine
#
# Therapy scenarios in `scenarios()`:
#   1. Untreated natural history
#   2. Prednisone 1 mg/kg/d → taper
#   3. Tocilizumab 162 mg SC q1w + low-dose prednisone
#   4. Anakinra 100 mg SC q24h
#   5. Ruxolitinib 10 mg PO BID + 10 mg/d prednisone
#   6. Azacitidine 75 mg/m² SC d1–7 q28 + supportive
#   7. Allogeneic HSCT (instantaneous clone replacement)
#
# Author: QSP-routine 2026-06-30  |  Compatible with mrgsolve ≥ 1.5
# ============================================================================

library(mrgsolve)

vexas_code <- '
$PROB
# VEXAS syndrome QSP model (23 ODEs)
# Author: QSP-routine 2026-06-30
# Time units: hours.  Concentration units indicated per compartment.

$PARAM @annotated
// ---- Disease genetics / clone ----
VAF0       :  0.55 : Baseline UBA1 mutant variant allele fraction (monocyte)
k_clone    : 0.0008: Clonal expansion rate (1/h) [age-driven]
k_HSCT     : 0.0   : HSCT-induced clone eradication rate (1/h)
k_Aza      : 0.0   : Azacitidine clone-modulating rate (1/h)

// ---- Proteostasis ----
kp_misf    : 0.05  : Misfolded protein production proportional to VAF (au/h)
kd_misf    : 0.02  : Baseline misfolded clearance (1/h)
hill_misf  : 2     : Hill coefficient (UPR/inflammation)

// ---- Mitochondrial / ER stress ----
k_ROS_in   : 0.03  : ROS generation by misfolded burden (au/h)
k_ROS_out  : 0.30  : ROS scavenging (1/h)
k_ERst_in  : 0.04  : ER-stress accumulation (au/h)
k_ERst_out : 0.10  : ER-stress resolution (1/h)

// ---- Cytokine synthesis (NF-κB / inflammasome / IFN) ----
k_IL1b_syn : 0.25  : IL-1β synthesis rate (pg/mL/h)
k_IL1b_deg : 0.40  : IL-1β degradation (1/h, t1/2 ≈ 1.7h)
k_IL6_syn  : 0.50  : IL-6 synthesis rate (pg/mL/h)
k_IL6_deg  : 0.40  : IL-6 degradation (1/h)
k_TNF_syn  : 0.30  : TNF-α synthesis (pg/mL/h)
k_TNF_deg  : 0.50  : TNF-α degradation (1/h, t1/2 ≈ 80 min)
k_IFN_syn  : 0.15  : IFN-α synthesis (au/h)
k_IFN_deg  : 0.35  : IFN-α degradation (1/h)
k_CXCL8_syn: 0.20  : CXCL8/IL-8 synthesis (pg/mL/h)
k_CXCL8_deg: 0.45  : CXCL8 degradation (1/h)
k_CCL2_syn : 0.10  : CCL2/MCP-1 synthesis (pg/mL/h)
k_CCL2_deg : 0.30  : CCL2 degradation (1/h)

// ---- Acute-phase reactants ----
k_CRP_syn  : 0.08  : CRP synthesis driven by IL-6 (mg/L/h)
EC50_IL6_CRP: 35   : IL-6 EC50 for CRP production (pg/mL)
k_CRP_deg  : 0.038 : CRP first-order elimination (1/h, t1/2 ≈ 18h)
k_FER_syn  : 0.05  : Ferritin synthesis from IL-6/IL-1β (ng/mL/h)
k_FER_deg  : 0.012 : Ferritin elimination (1/h, t1/2 ≈ 60h)
EC50_IL6_FER: 50   : IL-6 EC50 for ferritin (pg/mL)

// ---- Hematology ----
Hb0        : 13.0  : Baseline normal Hb (g/dL)
k_Hb_in    : 0.018 : Hb production rate (g/dL / h)
k_Hb_loss  : 0.0012: Hb intrinsic loss/turnover (1/h ~ 25-d t1/2)
Imax_Hb    : 0.65  : Max % suppression of erythropoiesis by inflammation
EC50_Hb    : 80    : Combined inflam index for half-max Hb suppression
PLT0       : 250   : Baseline platelets (×10⁹/L)
k_PLT_in   : 0.30  : Platelet production rate (×10⁹/L/h)
k_PLT_loss : 0.0058: Platelet elim (1/h, t1/2 ≈ 5d)
Imax_PLT   : 0.55  : Max % suppression of platelet production
EC50_PLT   : 70    : Inflam EC50 for PLT suppression
ANC0       : 4.5   : Baseline ANC (×10⁹/L)
k_ANC_in   : 0.025 : ANC production
k_ANC_loss : 0.10  : ANC turnover (1/h)
Imax_ANC   : 0.40  : Max % suppression
EC50_ANC   : 70    : Inflam EC50 for ANC

// ---- Clinical activity ----
k_Fev_syn  : 0.05  : Fever index buildup rate
k_Fev_dec  : 0.10  : Fever resolution (1/h)
k_Skin_syn : 0.03  : Skin/chondritis activity buildup
k_Skin_dec : 0.05  : Skin activity resolution
k_VTE_syn  : 0.0008: VTE risk buildup (1/h) driven by IL-6+ferritin
k_VTE_dec  : 0.005 : VTE risk decay (1/h)
EC50_VTE   : 60    : Inflam EC50 for VTE risk

// ---- HPA axis (steroid suppression) ----
HPA0       : 1.0   : Baseline endogenous cortisol output (au)
k_HPA_supp : 0.04  : Suppression of HPA by exogenous GC (1/h per mg/mL drug)
k_HPA_rec  : 0.005 : HPA recovery rate (1/h)

// ---- Drug PK: PREDNISONE (oral, 1-cpt) ----
ka_PRED    : 0.9   : Absorption (1/h)
V_PRED     : 60    : Volume (L)
CL_PRED    : 16    : Clearance (L/h, ~3h t1/2)
F_PRED     : 0.85  : Bioavailability

// ---- Drug PK: TOCILIZUMAB (TMDD-lite, SC) ----
ka_TOC     : 0.012 : SC absorption rate (1/h ~ 4d t1/2 abs)
V_TOC      : 5.0   : Vc (L)
CL_TOC     : 0.013 : Linear CL (L/h)
F_TOC      : 0.80  : SC bioavailability
Kd_TOC     : 0.5   : Tocilizumab-IL6R apparent Kd (mg/L)
IL6Rmax    : 5.0   : Apparent IL-6R density (mg/L)

// ---- Drug PK: ANAKINRA (SC, IL-1Ra recombinant) ----
ka_ANA     : 0.60  : Absorption (1/h)
V_ANA      : 12    : Volume (L)
CL_ANA     : 2.8   : CL (L/h, t1/2 ≈ 3-4h)
F_ANA      : 0.95  : SC F
Kd_ANA     : 1.0   : Anakinra-IL1R apparent Kd (mg/L)

// ---- Drug PK: CANAKINUMAB (SC mAb, anti-IL-1β) ----
ka_CAN     : 0.005 : Absorption (1/h)
V_CAN      : 6.0   : V (L)
CL_CAN     : 0.005 : CL (L/h, t1/2 ≈ 26d)
F_CAN      : 0.66  : F
Kd_CAN     : 0.3   : Canakinumab Kd (mg/L)

// ---- Drug PK: RUXOLITINIB (oral, 1-cpt) ----
ka_RUX     : 1.6   : Absorption (1/h)
V_RUX      : 75    : V (L)
CL_RUX     : 22    : CL (L/h, t1/2 ≈ 3h)
F_RUX      : 0.95  : F
Kd_RUX     : 100   : Apparent JAK Kd (μg/L)

// ---- Drug PK: AZACITIDINE (SC, 1-cpt) ----
ka_AZA     : 1.8   : Absorption (1/h)
V_AZA      : 76    : V (L)
CL_AZA     : 250   : CL (L/h, t1/2 ≈ 0.7h)
F_AZA      : 0.89  : SC F

// ---- Body / dosing helpers ----
WT         : 75    : Body weight (kg)
BSA        : 1.85  : Body surface area (m^2)

$CMT @annotated
// Disease/biology compartments
VAF        : UBA1 mutant variant allele fraction (unitless 0-1)
MISF       : Misfolded protein burden (au)
ROS        : ROS / mitochondrial stress (au)
ERST       : ER stress / UPR (au)
IL1B       : IL-1β (pg/mL)
IL6        : IL-6 (pg/mL)
TNFa       : TNF-α (pg/mL)
IFNa       : IFN-α (au)
CXCL8      : CXCL8 (pg/mL)
CCL2       : CCL2 (pg/mL)
CRP        : CRP (mg/L)
FER        : Ferritin (ng/mL)
HB         : Hb (g/dL)
PLT        : Platelet count (×10⁹/L)
ANC        : Absolute neutrophil count (×10⁹/L)
FEV        : Fever activity index (0-100)
SKIN       : Skin/chondritis activity index (0-100)
VTE        : VTE risk index (0-100)
HPA        : HPA axis output (fraction of baseline)
// Drug PK depots & central compartments
AGUT_PRED  : Prednisone gut (mg)
CC_PRED    : Prednisone central (mg)
ADEP_TOC   : Tocilizumab SC depot (mg)
CC_TOC     : Tocilizumab central (mg)
ADEP_ANA   : Anakinra SC depot (mg)
CC_ANA     : Anakinra central (mg)
ADEP_CAN   : Canakinumab SC depot (mg)
CC_CAN     : Canakinumab central (mg)
AGUT_RUX   : Ruxolitinib gut (mg)
CC_RUX     : Ruxolitinib central (mg)
ADEP_AZA   : Azacitidine SC depot (mg)
CC_AZA     : Azacitidine central (mg)

$MAIN
VAF_0   = VAF0;
MISF_0  = 30.0;
ROS_0   = 5.0;
ERST_0  = 6.0;
IL1B_0  = 8.0;
IL6_0   = 120.0;
TNFa_0  = 25.0;
IFNa_0  = 8.0;
CXCL8_0 = 60.0;
CCL2_0  = 280.0;
CRP_0   = 110.0;
FER_0   = 1800.0;
HB_0    = 9.5;
PLT_0   = 130.0;
ANC_0   = 5.5;
FEV_0   = 40.0;
SKIN_0  = 35.0;
VTE_0   = 25.0;
HPA_0   = 1.0;

$ODE
// ---- Drug concentrations ----
double C_PRED  = CC_PRED / V_PRED;        // mg/L prednisone
double C_TOC   = CC_TOC  / V_TOC;
double C_ANA   = CC_ANA  / V_ANA;
double C_CAN   = CC_CAN  / V_CAN;
double C_RUX   = CC_RUX  / V_RUX * 1000;  // μg/L
double C_AZA   = CC_AZA  / V_AZA;

// ---- Drug PD: anti-cytokine occupancy ----
double occ_TOC = C_TOC / (Kd_TOC + C_TOC);
double occ_ANA = C_ANA / (Kd_ANA + C_ANA);
double occ_CAN = C_CAN / (Kd_CAN + C_CAN);
double occ_RUX = C_RUX / (Kd_RUX + C_RUX);

// ---- Effective free cytokine signaling ----
double IL6_eff   = IL6 * (1.0 - 0.95*occ_TOC);
double IL1_eff   = IL1B * (1.0 - 0.90*occ_ANA - 0.85*occ_CAN);
double IFN_eff   = IFNa * (1.0 - 0.70*occ_RUX);
double NFkB_act  = (1 + C_PRED/(0.05 + C_PRED) > 0 ? (1.0 - 0.6*C_PRED/(0.05 + C_PRED)) : 1.0);

// Combined inflammation index (for hematology suppression)
double INFLAM    = IL6_eff + 2.0*IL1_eff + 0.5*TNFa + 0.3*IFN_eff;

// ---- Disease biology ODEs ----
dxdt_VAF  = (k_clone*VAF*(1-VAF)) - k_HSCT*VAF - k_Aza*VAF;

dxdt_MISF = kp_misf*100*VAF - kd_misf*MISF;
dxdt_ROS  = k_ROS_in*MISF - k_ROS_out*ROS;
dxdt_ERST = k_ERst_in*MISF - k_ERst_out*ERST;

// Cytokines: NF-κB driven (steroids ↓), inflammasome (IL-1β), IFN
dxdt_IL1B  = k_IL1b_syn*(ROS/(ROS + 5))*100*NFkB_act - k_IL1b_deg*IL1B;
dxdt_IL6   = k_IL6_syn *(ERST/(ERST + 6))*100*NFkB_act - k_IL6_deg*IL6;
dxdt_TNFa  = k_TNF_syn *(MISF/(MISF + 40))*100*NFkB_act - k_TNF_deg*TNFa;
dxdt_IFNa  = k_IFN_syn *(ROS/(ROS + 4))*100*(1.0 - 0.7*occ_RUX) - k_IFN_deg*IFNa;
dxdt_CXCL8 = k_CXCL8_syn*(IL1_eff/(IL1_eff + 5))*100 - k_CXCL8_deg*CXCL8;
dxdt_CCL2  = k_CCL2_syn *(IL6_eff/(IL6_eff + 40))*100 - k_CCL2_deg*CCL2;

// Acute-phase
dxdt_CRP   = k_CRP_syn*(IL6_eff/(EC50_IL6_CRP + IL6_eff))*100 - k_CRP_deg*CRP;
dxdt_FER   = k_FER_syn*((IL6_eff+IL1_eff)/(EC50_IL6_FER + IL6_eff + IL1_eff))*100 - k_FER_deg*FER;

// Hematology — suppression by inflammation
double sup_Hb  = Imax_Hb  * INFLAM / (EC50_Hb + INFLAM);
double sup_PLT = Imax_PLT * INFLAM / (EC50_PLT + INFLAM);
double sup_ANC = Imax_ANC * INFLAM / (EC50_ANC + INFLAM);

dxdt_HB  = k_Hb_in *(1 - sup_Hb )*Hb0  - k_Hb_loss *HB;
dxdt_PLT = k_PLT_in*(1 - sup_PLT)*PLT0 - k_PLT_loss*PLT;
dxdt_ANC = k_ANC_in*(1 - sup_ANC)*ANC0 - k_ANC_loss*ANC;

// Clinical activity scores
dxdt_FEV  = k_Fev_syn *(IL1_eff/(IL1_eff + 10) + IL6_eff/(IL6_eff + 50))*50 - k_Fev_dec *FEV;
dxdt_SKIN = k_Skin_syn*(CXCL8/(CXCL8 + 30) + IL1_eff/(IL1_eff + 10))*50 - k_Skin_dec*SKIN;
dxdt_VTE  = k_VTE_syn *((IL6_eff+FER/100)/(EC50_VTE + IL6_eff + FER/100))*100 - k_VTE_dec*VTE;

// HPA axis
dxdt_HPA  = -k_HPA_supp*C_PRED*HPA + k_HPA_rec*(1 - HPA);

// ---- Drug PK ODEs ----
dxdt_AGUT_PRED = -ka_PRED*AGUT_PRED;
dxdt_CC_PRED   =  ka_PRED*AGUT_PRED - CL_PRED/V_PRED*CC_PRED;

dxdt_ADEP_TOC  = -ka_TOC*ADEP_TOC;
dxdt_CC_TOC    =  ka_TOC*ADEP_TOC - CL_TOC/V_TOC*CC_TOC;

dxdt_ADEP_ANA  = -ka_ANA*ADEP_ANA;
dxdt_CC_ANA    =  ka_ANA*ADEP_ANA - CL_ANA/V_ANA*CC_ANA;

dxdt_ADEP_CAN  = -ka_CAN*ADEP_CAN;
dxdt_CC_CAN    =  ka_CAN*ADEP_CAN - CL_CAN/V_CAN*CC_CAN;

dxdt_AGUT_RUX  = -ka_RUX*AGUT_RUX;
dxdt_CC_RUX    =  ka_RUX*AGUT_RUX - CL_RUX/V_RUX*CC_RUX;

dxdt_ADEP_AZA  = -ka_AZA*ADEP_AZA;
dxdt_CC_AZA    =  ka_AZA*ADEP_AZA - CL_AZA/V_AZA*CC_AZA;

$CAPTURE @annotated
C_PRED   : Prednisone plasma conc (mg/L)
C_TOC    : Tocilizumab plasma conc (mg/L)
C_ANA    : Anakinra plasma conc (mg/L)
C_CAN    : Canakinumab plasma conc (mg/L)
C_RUX    : Ruxolitinib plasma conc (μg/L)
C_AZA    : Azacitidine plasma conc (mg/L)
INFLAM   : Composite inflammation index
IL6_eff  : Free/biologically-active IL-6 (pg/mL)
IL1_eff  : Free IL-1β (pg/mL)
IFN_eff  : Free IFN-α (au)
occ_TOC  : IL-6R occupancy
occ_ANA  : IL-1R occupancy
occ_RUX  : JAK occupancy
'

# ---- Compile ----------------------------------------------------------------
mod_vexas <- mcode("vexas_qsp", vexas_code)

# ---- Helpers ---------------------------------------------------------------
prednisone_taper <- function(start_mg = 60, taper_weeks = 12, end_mg = 5) {
  weeks <- 0:taper_weeks
  doses <- pmax(seq(start_mg, end_mg, length.out = length(weeks)), end_mg)
  data.frame(time = weeks*168, cmt = "AGUT_PRED", amt = doses, evid = 1, ii = 24, addl = 6)
}

scenarios <- function() {
  list(
    "1_untreated"   = NULL,
    "2_pred_taper"  = prednisone_taper(60, 12, 5),
    "3_toci_lowGC"  = rbind(
      data.frame(time = seq(0, 84*24, by = 7*24), cmt = "ADEP_TOC", amt = 162, evid = 1, ii = 0, addl = 0),
      data.frame(time = seq(0, 84*24, by = 24),  cmt = "AGUT_PRED", amt = 10, evid = 1, ii = 0, addl = 0)
    ),
    "4_anakinra"    = data.frame(time = seq(0, 84*24, by = 24), cmt = "ADEP_ANA", amt = 100, evid = 1, ii = 0, addl = 0),
    "5_ruxo_pred"   = rbind(
      data.frame(time = seq(0, 84*24, by = 12), cmt = "AGUT_RUX", amt = 10, evid = 1, ii = 0, addl = 0),
      data.frame(time = seq(0, 84*24, by = 24), cmt = "AGUT_PRED", amt = 10, evid = 1, ii = 0, addl = 0)
    ),
    "6_azacitidine" = data.frame(time = c(outer(seq(0, 6)*24, seq(0, 168, by = 28*24), `+`)),
                                  cmt = "ADEP_AZA", amt = 75*1.85, evid = 1, ii = 0, addl = 0),
    "7_HSCT"        = data.frame(time = 24, cmt = "VAF", amt = 0, evid = 2, # set k_HSCT high transiently
                                  ii = 0, addl = 0)  # use $TABLE override in user code
  )
}

# ---- Example run ------------------------------------------------------------
if (interactive()) {
  library(dplyr); library(ggplot2); library(tidyr)
  evts <- scenarios()[["3_toci_lowGC"]]
  out  <- mod_vexas %>%
    ev(evts) %>%
    mrgsim(end = 84*24, delta = 6) %>% as.data.frame()
  out %>% select(time, IL6, IL6_eff, CRP, FER, HB, PLT) %>%
    pivot_longer(-time) %>%
    ggplot(aes(time/24, value)) + geom_line() + facet_wrap(~name, scales = "free_y") +
    labs(x = "Days", title = "VEXAS — Tocilizumab + low-dose prednisone")
}

# ---- Parameter sources / calibration notes ---------------------------------
# • UBA1 VAF ranges: Beck et al. NEJM 2020 (median ~75% monocyte);
#   Georgin-Lavialle et al. AJH 2024.
# • IL-6 elevation: median 80–250 pg/mL in active disease (Patel 2022).
# • CRP: typically 80–200 mg/L during flare.
# • Ferritin: 1000–10,000 ng/mL; HLH-like crises >10,000.
# • Hb baseline 8–10 g/dL macrocytic (MCV >100); PLT 80–150×10⁹/L.
# • Tocilizumab PK: V≈5L, CL≈10–15 mL/h, t1/2 11–13 d (SC, EMA/FDA labels).
# • Anakinra PK: t1/2 4–6h, F≈95% SC (Kineret label).
# • Canakinumab PK: t1/2 26 d, F≈66% SC (Ilaris label).
# • Ruxolitinib PK: t1/2 3h, CL 22 L/h (Jakafi label).
# • Azacitidine PK: t1/2 ~0.7h (Vidaza label, oral CC-486 differs).
# • Prednisone PK: F≈0.85, CL 6–8 L/h/70kg, t1/2 3h (active prednisolone).
# • 5-yr survival ~50% untreated, ~60–70% on biologic+HSCT (Hines 2023).
# • Allogeneic HSCT modeled as instantaneous VAF→0 (Hadjadj 2024 cohort).
# ============================================================================
