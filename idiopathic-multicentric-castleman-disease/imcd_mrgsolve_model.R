# =============================================================================
# Idiopathic Multicentric Castleman Disease (iMCD) — mrgsolve QSP Model
# IL-6 cytokine storm with siltuximab, tocilizumab, sirolimus, rituximab,
# anakinra, ruxolitinib, CHOP, prednisone
#
# 19 ODE compartments, 7 treatment scenarios
# Calibration: CONCERT trial (van Rhee 2014 Lancet Oncol) for siltuximab
#              ACTEMRA studies (Nishimoto 2005 Blood) for tocilizumab
#              TAFRO subtype (Iwaki 2016, Fajgenbaum 2019 JCI) for sirolimus
#
# Author: QSP CCR Library · CDCN 2017 diagnostic criteria
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

code <- '
$PROB
# iMCD QSP Model (Idiopathic Multicentric Castleman Disease)
- IL-6 cytokine storm + lymphadenopathy + acute phase
- Drugs : siltuximab, tocilizumab, sirolimus, rituximab, anakinra,
          ruxolitinib, CHOP-doxo, CHOP-cyclo, prednisone
- Endpoints : CRP, IL-6 (total/free), Hb, IgG, LN tumor burden,
              VEGF, platelet, CDCN response, TAFRO severity, OS hazard

$PARAM @annotated
// ---------------- siltuximab (anti-IL-6 mAb, FDA 2014) ------------------
CL_SILT   :  0.232  : siltuximab CL (L/d) van Rhee 2014 PopPK
V1_SILT   :  4.47   : Vc (L)
V2_SILT   :  2.49   : Vp (L)
Q_SILT    :  0.527  : Q (L/d)
Kd_SILT   :  1e-3   : Kd IL-6 binding (nM ≈ 1 pM) Kurzrock 2013
MW_SILT   :  145000 : Da
DOSE_SILT :  11     : mg/kg q3w label dose
// ---------------- tocilizumab (anti-IL-6R mAb) -------------------------
CL_TOCZ   :  0.20   : (L/d) Frey 2010
V1_TOCZ   :  4.1    :
V2_TOCZ   :  2.5    :
Q_TOCZ    :  0.50   :
Kd_TOCZ   :  2.5    : nM (IL-6R)
DOSE_TOCZ :  8      : mg/kg q2w
// ---------------- sirolimus (mTORC1) ----------------------------------
CL_SIRO   :  8.2    : L/h  -> 197 L/d MacDonald 2000
V1_SIRO   :  12     : L
V2_SIRO   :  100    : L
Q_SIRO    :  10     : L/h
F_SIRO    :  0.14   : oral bioavail
Ka_SIRO   :  2.0    : 1/h
DOSE_SIRO :  2      : mg PO QD (target trough 6-14 ng/mL)
EC50_SIRO :  8      : ng/mL
// ---------------- rituximab (anti-CD20) -------------------------------
CL_RTX    :  0.32   : L/d
V1_RTX    :  3.3    :
V2_RTX    :  4.5    :
Q_RTX     :  0.42   :
DOSE_RTX  :  375    : mg/m2 weekly x4 induction
// ---------------- anakinra (IL-1Rα) -----------------------------------
CL_ANA    :  9.0    : L/h
V_ANA     :  20     : L
F_ANA     :  0.95   : SC
Ka_ANA    :  0.40   : 1/h
DOSE_ANA  :  100    : mg SC QD
// ---------------- ruxolitinib (JAK1/2) --------------------------------
CL_RUX    :  18     : L/h
V_RUX     :  72     : L
F_RUX     :  0.95   :
Ka_RUX    :  3.0    : 1/h
DOSE_RUX  :  20     : mg PO BID
IC50_RUX  :  300    : nM JAK1/2
// ---------------- CHOP doxorubicin ------------------------------------
CL_DOXO   :  45     : L/h
V_DOXO    :  500    : L
DOSE_DOXO :  50     : mg/m2 q21d
// ---------------- CHOP cyclophosphamide -------------------------------
CL_CYC    :  6      : L/h
V_CYC     :  35     : L
DOSE_CYC  :  750    : mg/m2 q21d
// ---------------- Prednisone ------------------------------------------
CL_PRED   :  5.6    : L/h
V_PRED    :  35     : L
Ka_PRED   :  2.5    :
F_PRED    :  0.85   :
EC50_PRED :  40     : ng/mL
DOSE_PRED :  60     : mg PO QD (≈1 mg/kg)
//
// ---------------- IL-6 axis & turnover --------------------------------
kIL6_base : 0.20    : baseline IL-6 production (pg/mL/d) — healthy
kIL6_iMCD : 2.5     : iMCD over-production multiplier
kIL6_deg  : 8       : IL-6 elimination 1/d (half-life ~2 h, but in serum apparent)
sIL6R_b   : 50      : ng/mL baseline soluble IL-6R
//
// ---------------- Lymph node / plasmablast ----------------------------
kgrow_LN  : 0.012   : 1/d growth rate of LN size (cm composite)
LN_max    : 25      : cm composite
kshrink_LN: 0.005   : 1/d natural turnover
//
kprol_PB  : 0.20    : plasmablast proliferation 1/d (under IL-6)
kdeath_PB : 0.10    : 1/d basal apoptosis
PB_base   : 0.05    : baseline plasmablast fraction (% of LN)
//
kprol_Bmem: 0.08    : memory B turnover
kdeath_Bmem: 0.04   :
RTX_kill  : 0.030   : 1/(d * (mg/L)) CD20 lysis rate constant
//
// ---------------- Acute phase ----------------------------------------
kCRP_in   : 12      : mg/L/d max CRP synthesis under IL-6
kCRP_out  : 0.55    : 1/d (CRP half-life ~ 19 h)
EC50_CRP  : 5       : pg/mL IL-6
//
kHb_in    : 0.16    : g/dL/d homeostasis
kHb_out   : 0.012   : 1/d (RBC lifespan ~120 d)
Hb_max    : 14.5    : g/dL
Hepcidin_EC50 : 30  : pg/mL IL-6 driving hepcidin
//
kIgG_in   : 0.6     : g/dL/d
kIgG_out  : 0.035   : 1/d (IgG t1/2 ~21 d)
IgG_max   : 5.0     : g/dL polyclonal max under disease
//
// ---------------- VEGF / TAFRO --------------------------------------
kVEGF_in  : 50      : pg/mL/d
kVEGF_out : 1.5     : 1/d
EC50_VEGF : 20      : pg/mL IL-6
//
kAnasarca_in : 0.20 :  units/d
kAnasarca_out: 0.10 :  1/d
EC50_VEGFan  : 200  : pg/mL VEGF
//
// ---------------- Platelet ------------------------------------------
kPlt_in   : 25      : x10^9/L per d
kPlt_out  : 0.10    : 1/d
Plt_max   : 400     :
Plt_TAFRO_kill : 0.15 : platelet consumption rate driven by IL-6/VEGF
//
// ---------------- mTORC1 / TAFRO disease driver ---------------------
kmTOR_in  : 0.10    : arbitrary signal/d
kmTOR_out : 0.05    : 1/d
//
// ---------------- Hazard / safety ----------------------------------
h0_OS     : 5e-5    : 1/d baseline 5-yr OS hazard ~ 0.09
beta_CRP  : 0.0008  : CRP hazard coefficient

$CMT @annotated
SILT_C   : siltuximab central (mg)
SILT_P   : siltuximab peripheral (mg)
TOCZ_C   : tocilizumab central (mg)
TOCZ_P   : tocilizumab peripheral (mg)
SIRO_GUT : sirolimus gut (mg)
SIRO_C   : sirolimus central (mg)
SIRO_P   : sirolimus peripheral (mg)
RTX_C    : rituximab central (mg)
RTX_P    : rituximab peripheral (mg)
ANA_SC   : anakinra SC (mg)
ANA_C    : anakinra central (mg)
RUX_GUT  : ruxolitinib gut (mg)
RUX_C    : ruxolitinib central (mg)
DOXO_C   : doxorubicin central (mg)
CYC_C    : cyclophosphamide central (mg)
PRED_GUT : prednisone gut (mg)
PRED_C   : prednisone central (mg)
IL6_T    : IL-6 serum (pg/mL)  — TOTAL (free + bound)
IL6_F    : IL-6 free  (pg/mL)
LN       : lymph node composite size (cm)
PB       : plasmablast burden  (fraction LN)
Bmem     : memory B-cell pool (% baseline)
CRP      : C-reactive protein (mg/L)
Hb       : Hemoglobin (g/dL)
IgG      : Polyclonal IgG (g/dL)
VEGF     : VEGF-A (pg/mL)
Anasarca : Anasarca / 3rd-space fluid score (0–10)
Plt      : Platelet x10^9/L
mTOR     : mTORC1 activity (0–1)
HAZ      : Cumulative OS hazard (death by AE/disease)

$MAIN
// initial conditions
IL6_T_0   = 60;            // pg/mL active disease (CONCERT median)
IL6_F_0   = 30;
LN_0      = 8;             // cm composite (significant lymphadenopathy)
PB_0      = 0.12;          // 12% of LN are plasmablasts
Bmem_0    = 100;
CRP_0     = 120;           // mg/L active iMCD
Hb_0      = 9.5;           // g/dL anemia of inflammation
IgG_0     = 4.0;           // hypergammaglobulinemia
VEGF_0    = 800;           // pg/mL
Anasarca_0= 2.5;
Plt_0     = 300;
mTOR_0    = 0.55;
HAZ_0     = 0;

$ODE
// ---------- PK ODEs ----------
double CL_SILT_d = CL_SILT;   // L/d
double k10_SILT  = CL_SILT_d / V1_SILT;
double k12_SILT  = Q_SILT    / V1_SILT;
double k21_SILT  = Q_SILT    / V2_SILT;
double Csilt    = SILT_C / V1_SILT;     // mg/L
dxdt_SILT_C = -(k10_SILT + k12_SILT) * SILT_C + k21_SILT * SILT_P;
dxdt_SILT_P = k12_SILT * SILT_C - k21_SILT * SILT_P;

double k10_TOCZ = CL_TOCZ / V1_TOCZ;
double k12_TOCZ = Q_TOCZ  / V1_TOCZ;
double k21_TOCZ = Q_TOCZ  / V2_TOCZ;
double Ctocz   = TOCZ_C / V1_TOCZ;
dxdt_TOCZ_C = -(k10_TOCZ + k12_TOCZ) * TOCZ_C + k21_TOCZ * TOCZ_P;
dxdt_TOCZ_P = k12_TOCZ * TOCZ_C - k21_TOCZ * TOCZ_P;

double k10_SIRO = (CL_SIRO * 24) / V1_SIRO;     // /d
double k12_SIRO = (Q_SIRO  * 24) / V1_SIRO;
double k21_SIRO = (Q_SIRO  * 24) / V2_SIRO;
double Csiro    = SIRO_C / V1_SIRO * 1000;       // ng/mL
dxdt_SIRO_GUT = -Ka_SIRO * SIRO_GUT;
dxdt_SIRO_C   =  Ka_SIRO * SIRO_GUT * F_SIRO - (k10_SIRO + k12_SIRO) * SIRO_C + k21_SIRO * SIRO_P;
dxdt_SIRO_P   =  k12_SIRO * SIRO_C - k21_SIRO * SIRO_P;

double k10_RTX = CL_RTX / V1_RTX;
double k12_RTX = Q_RTX  / V1_RTX;
double k21_RTX = Q_RTX  / V2_RTX;
double Crtx   = RTX_C / V1_RTX;
dxdt_RTX_C = -(k10_RTX + k12_RTX) * RTX_C + k21_RTX * RTX_P;
dxdt_RTX_P = k12_RTX * RTX_C - k21_RTX * RTX_P;

double k10_ANA = (CL_ANA * 24) / V_ANA;
double Cana    = ANA_C / V_ANA * 1000;            // ng/mL
dxdt_ANA_SC = -Ka_ANA * 24 * ANA_SC;
dxdt_ANA_C  =  Ka_ANA * 24 * ANA_SC * F_ANA - k10_ANA * ANA_C;

double k10_RUX = (CL_RUX * 24) / V_RUX;
double Crux    = RUX_C / V_RUX * 1000;
dxdt_RUX_GUT = -Ka_SIRO * 24 * RUX_GUT;
dxdt_RUX_C   =  Ka_SIRO * 24 * RUX_GUT * F_RUX - k10_RUX * RUX_C;

double k10_DOXO = (CL_DOXO * 24) / V_DOXO;
double Cdoxo = DOXO_C / V_DOXO;
dxdt_DOXO_C = -k10_DOXO * DOXO_C;

double k10_CYC = (CL_CYC * 24) / V_CYC;
double Ccyc   = CYC_C / V_CYC;
dxdt_CYC_C = -k10_CYC * CYC_C;

double k10_PRED = (CL_PRED * 24) / V_PRED;
double Cpred = PRED_C / V_PRED * 1000;
dxdt_PRED_GUT = -Ka_PRED * 24 * PRED_GUT;
dxdt_PRED_C   =  Ka_PRED * 24 * PRED_GUT * F_PRED - k10_PRED * PRED_C;

// ---------- IL-6 axis ----------
// Disease drive (IL-6 production from LN plasmablasts + monocytes + mTOR)
double DiseaseDrive = kIL6_iMCD * (LN / LN_0) * (1 + 0.6 * (mTOR - 0.3));
double IL6_synth    = kIL6_base + DiseaseDrive * 1.0;       // pg/mL/d

// Free IL-6 = total / (1 + siltuximab effect)
// siltuximab binds IL-6 directly with Kd ~ 1 pM (negligible)
double Csilt_molar = Csilt / 145000 * 1e9;                  // nM
double IL6_F_calc  = IL6_T * 1.0 / (1.0 + Csilt_molar / Kd_SILT);

// Tocilizumab blocks IL-6R signal (post-receptor block)
double Ctocz_molar = Ctocz / 145000 * 1e9;
double TOCZ_block  = 1.0 / (1.0 + Ctocz_molar / Kd_TOCZ);

// Ruxolitinib blocks JAK1/2
double JAK_block   = 1.0 / (1.0 + Crux / IC50_RUX);

// Prednisone blocks NF-kB → IL-6 production
double Pred_block  = 1.0 / (1.0 + Cpred / EC50_PRED);

dxdt_IL6_T = IL6_synth * Pred_block - kIL6_deg * IL6_T;
dxdt_IL6_F = IL6_synth * Pred_block * 1.0 / (1.0 + Csilt_molar / Kd_SILT) - kIL6_deg * IL6_F;

double IL6_signal = IL6_F * TOCZ_block * JAK_block;          // effective signal at gp130 + post-receptor

// ---------- mTORC1 ----------
double Sirol_block = 1.0 / (1.0 + Csiro / EC50_SIRO);
dxdt_mTOR = kmTOR_in * (1 + 0.4 * (IL6_signal / 50)) * Sirol_block - kmTOR_out * mTOR;

// ---------- Plasmablast & B-cell ----------
double IL6_prol = IL6_signal / (IL6_signal + 10);            // saturating
double Bortz_kill = 0;       // placeholder if bortezomib added
double cyto_kill = 0.05 * (Cdoxo / (Cdoxo + 0.1)) + 0.04 * (Ccyc / (Ccyc + 5));
dxdt_PB = kprol_PB * (IL6_prol + 0.6 * mTOR) * (1 - PB / 0.8) - (kdeath_PB + cyto_kill) * PB;

dxdt_Bmem = kprol_Bmem * (1 - Bmem/120) - (kdeath_Bmem + RTX_kill * Crtx) * Bmem;

// ---------- Lymph node size ----------
double LN_drive = 0.7 * (PB / 0.10) + 0.4 * mTOR + 0.2 * (IL6_signal / 30);
dxdt_LN = kgrow_LN * LN_drive * (1 - LN/LN_max) - (kshrink_LN + 0.5 * cyto_kill) * LN;

// ---------- Acute phase ----------
dxdt_CRP = kCRP_in * IL6_signal / (IL6_signal + EC50_CRP) - kCRP_out * CRP;

// hepcidin / anemia
double Hepc = IL6_signal / (IL6_signal + Hepcidin_EC50);
dxdt_Hb  = kHb_in * (1 - Hepc) * (1 - Hb / Hb_max) - kHb_out * Hb;

// IgG polyclonal
dxdt_IgG = kIgG_in * (PB / 0.10) * (1 - IgG / IgG_max) - (kIgG_out + 0.1 * RTX_kill * Crtx) * IgG;

// VEGF / Anasarca / TAFRO axis (sirolimus reduces)
dxdt_VEGF = kVEGF_in * IL6_signal / (IL6_signal + EC50_VEGF) * (0.7 + 0.6 * mTOR) - kVEGF_out * VEGF;
dxdt_Anasarca = kAnasarca_in * VEGF / (VEGF + EC50_VEGFan) - kAnasarca_out * Anasarca;

// Platelet — TAFRO consumption
double TAFRO_consum = Plt_TAFRO_kill * (VEGF / 800) * (mTOR / 0.5);
dxdt_Plt = kPlt_in * (1 - Plt/Plt_max) - (kPlt_out + TAFRO_consum) * Plt;

// ---------- Hazard ----------
dxdt_HAZ = h0_OS * (1 + beta_CRP * CRP);

$TABLE
double IL6_serum = IL6_T;
double IL6_free  = IL6_F;
double LN_size   = LN;
double CRP_lab   = CRP;
double Hb_lab    = Hb;
double IgG_lab   = IgG;
double VEGF_lab  = VEGF;
double Plt_lab   = Plt;
double Anasarca_lab = Anasarca;
double mTOR_act  = mTOR;
double Survival  = exp(-HAZ);

// CDCN response composite (lower is better)
double CDCN_resp = 0.25*(IL6_serum/60) + 0.25*(CRP/120) + 0.20*(LN/8) +
                   0.10*(1 - Hb/12) + 0.10*(IgG/4) + 0.10*(Anasarca/2.5);

$CAPTURE Csilt Ctocz Csiro Crtx Crux Cdoxo Ccyc Cpred Cana
         IL6_serum IL6_free LN_size CRP_lab Hb_lab IgG_lab VEGF_lab Plt_lab Anasarca_lab mTOR_act
         Survival CDCN_resp
'

mod <- mcode("imcd", code)

# =============================================================================
# Treatment scenarios
# =============================================================================
sim_one <- function(label, ev) {
  mrgsim(mod, events = ev, end = 365, delta = 1) %>%
    as_tibble() %>% mutate(scenario = label)
}

# ---- Scenario 1: Untreated natural history -----
ev_none <- ev(amt = 0, cmt = "SILT_C", time = 0)
s1 <- sim_one("S1: Untreated", ev_none)

# ---- Scenario 2: Siltuximab 11 mg/kg IV q3w (CONCERT trial label) -----
# 70 kg -> 770 mg per dose
ev_silt <- ev(amt = 770, cmt = "SILT_C", time = seq(0, 360, 21))
s2 <- sim_one("S2: Siltuximab 11 mg/kg q3w", ev_silt)

# ---- Scenario 3: Tocilizumab 8 mg/kg IV q2w (Nishimoto regimen) -----
ev_tocz <- ev(amt = 560, cmt = "TOCZ_C", time = seq(0, 360, 14))
s3 <- sim_one("S3: Tocilizumab 8 mg/kg q2w", ev_tocz)

# ---- Scenario 4: Sirolimus PO QD (TAFRO Fajgenbaum 2019 JCI) -----
ev_siro <- ev(amt = 2, cmt = "SIRO_GUT", ii = 1, addl = 364)
s4 <- sim_one("S4: Sirolimus 2 mg QD (TAFRO)", ev_siro)

# ---- Scenario 5: Rituximab 375 mg/m2 weekly x4 + prednisone -----
ev_rtx_pred <- bind_rows(
  ev(amt = 750, cmt = "RTX_C", time = c(0, 7, 14, 21)),     # 1.5 m2 ~750 mg
  ev(amt = 60,  cmt = "PRED_GUT", ii = 1, addl = 27, time = 0)
)
s5 <- sim_one("S5: Rituximab + Prednisone", ev_rtx_pred)

# ---- Scenario 6: CHOP-like + Siltuximab (induction) -----
ev_chop_silt <- bind_rows(
  ev(amt = 75 * 1.5,  cmt = "DOXO_C", time = seq(0, 6*21, 21)),
  ev(amt = 750 * 1.5, cmt = "CYC_C",  time = seq(0, 6*21, 21)),
  ev(amt = 100,       cmt = "PRED_GUT", ii = 1, addl = 4, time = seq(0, 6*21, 21)),
  ev(amt = 770,       cmt = "SILT_C", time = seq(0, 360, 21))
)
s6 <- sim_one("S6: CHOP + Siltuximab", ev_chop_silt)

# ---- Scenario 7: Triple therapy (Sirolimus + Siltuximab + Anakinra) for refractory TAFRO -----
ev_triple <- bind_rows(
  ev(amt = 770, cmt = "SILT_C", time = seq(0, 360, 21)),
  ev(amt = 2,   cmt = "SIRO_GUT", ii = 1, addl = 364),
  ev(amt = 100, cmt = "ANA_SC", ii = 1, addl = 364)
)
s7 <- sim_one("S7: Siltuximab + Sirolimus + Anakinra", ev_triple)

# Combine and quick plot
all_sim <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

# Example plot
if (interactive()) {
  ggplot(all_sim, aes(time, CRP_lab, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(title = "iMCD QSP: CRP trajectory by treatment scenario",
         x = "Day", y = "CRP (mg/L)") +
    theme_bw()
}

# =============================================================================
# Calibration notes (representative literature anchors)
# =============================================================================
# 1) CONCERT trial (van Rhee 2014 Lancet Oncol): siltuximab 11 mg/kg q3w
#    durable tumor + symptomatic response 34% vs 0% placebo (P=0.0012).
#    CRP normalized to < 10 mg/L by week 6; IL-6 increased due to clearance
#    block (total IL-6 rises because free IL-6 is sequestered).
#
# 2) Nishimoto 2005 Blood: tocilizumab 8 mg/kg q2w —
#    CRP < 10 mg/L by week 2; IgG, fibrinogen, ESR normalize by week 4;
#    lymph node shrinkage by 50% at 6 months in ~60% of patients.
#
# 3) Fajgenbaum 2019 JCI Insight: sirolimus rescues 3/3 refractory TAFRO patients,
#    target trough 6-14 ng/mL — rationale for mTORC1 / Tfh hyperactivation.
#
# 4) van Rhee 2018 Blood Adv consensus treatment guidelines — 1st line siltuximab
#    (or tocilizumab where siltuximab unavailable). Sirolimus for refractory.
#
# 5) Dispenzieri & Fajgenbaum 2020 Blood: 5-year OS ≈ 65% (TAFRO worse, IPL better).
# =============================================================================
