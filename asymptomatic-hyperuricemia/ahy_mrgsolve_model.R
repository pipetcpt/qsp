## ============================================================
## Asymptomatic Hyperuricemia (AHY) — QSP mrgsolve Model
## 무증상 고요산혈증 정량적 시스템 약리학 ODE 모델
##
## Compartments (19):
##   1.  UA_plasma      — 혈청 요산 (Serum Uric Acid, mg/dL)
##   2.  UA_tissue      — 조직 요산 풀 (Tissue UA pool, mg)
##   3.  XO_free        — 자유 XO 활성 (Free XO activity, nmol/min/mg)
##   4.  Oxypurinol     — 옥시퓨리놀 농도 (allopurinol active metabolite, mg/L)
##   5.  Febuxostat_C   — 페북소스타트 혈중 농도 (mg/L)
##   6.  Uricosuric_C   — 요산배설촉진제 혈중 농도 (mg/L)
##   7.  URAT1_free     — 활성 URAT1 분획 (Unbound URAT1 fraction)
##   8.  UrinaryUA      — 요중 요산 (Urinary UA pool, mg/day)
##   9.  MSU_depot      — MSU 결정 침착량 (Crystal depot, mg)
##   10. Endothelial_fn — 내피 기능 (Endothelial function, 0-1)
##   11. NO_level       — 산화질소 수준 (NO, rel. units)
##   12. BP             — 혈압 (Mean arterial pressure, mmHg)
##   13. GFR            — 사구체 여과율 (eGFR, mL/min/1.73m²)
##   14. IL1beta        — IL-1β (pg/mL)
##   15. CRP            — 고감도 CRP (hs-CRP, mg/L)
##   16. InsulinResist  — 인슐린 저항성 지수 (HOMA-IR, rel.)
##   17. CV_risk_score  — 누적 심혈관 위험 점수 (rel. units)
##   18. Tophus_vol     — 토파이 부피 (Tophus volume, mm³)
##   19. ABCG2_frac     — 기능성 ABCG2 분획 (Functional ABCG2)
##
## Parameter Calibration Reference:
##   - Dalbeth et al., Nat Rev Dis Primers 2019 (gout/UA physiology)
##   - Perez-Ruiz et al., Ann Rheum Dis 2011 (XO inhibition PK/PD)
##   - Becker et al., NEJM 2005 (febuxostat Phase III)
##   - FitzGerald et al., Arthritis Rheum 2011 (URAT1/ABCG2)
##   - Dalbeth et al., Lancet 2021 (treat-to-target)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- mrgsolve model code block ----
ahy_model_code <- '
$PROB
Asymptomatic Hyperuricemia (AHY) QSP Model
무증상 고요산혈증 ODE 시뮬레이션
19 compartments: UA pool, XO, drug PK, CV/renal effects

$PARAM @annotated
// --- Purine Production Parameters ---
k_prod_base  : 700   : mg/day, 기저 요산 생성 속도 (endogenous + dietary)
k_dietary    : 200   : mg/day, 식이 퓨린에서의 요산 생성 기여도
fructose_load: 0     : 과당 부하 지수 (0=없음, 1=고과당 식이)
alcohol_use  : 0     : 알코올 사용 지수 (0=없음, 1=음주)
k_fructose   : 0.15  : 과당→ 요산 생성 계수
k_alcohol    : 0.10  : 알코올→ 요산 생성 계수

// --- XO Activity Parameters ---
XO_Vmax      : 1.0   : XO 최대 반응속도 (nmol/min/mg, normalized)
XO_Km        : 50    : XO Km (μM, for hypoxanthine)
k_XO_deg     : 0.02  : XO 분해속도 (/day)
k_XO_syn     : 0.02  : XO 합성속도 (/day, homeostatic)

// --- Renal Excretion Parameters ---
k_renal_base : 400   : mg/day, 기저 신장 배설 속도
k_URAT1      : 0.70  : URAT1 재흡수 분획 (S1+S3)
k_OAT_secr   : 0.50  : OAT 분비 분획 (S2)
k_ABCG2_ren  : 0.15  : 신장 ABCG2 분비 분획
GFR_0        : 90    : mL/min/1.73m², 기저 GFR
// --- Intestinal Handling ---
k_gut_ABCG2  : 0.33  : 장관 ABCG2 분비 분획 (총 배설의 33%)

// --- Drug PK Parameters: Allopurinol / Oxypurinol ---
F_allo       : 0.90  : allopurinol 생체이용률
ka_allo      : 3.0   : /h, 흡수속도상수
CL_allo      : 18    : L/h, 알로퓨리놀 청소율
Vd_allo      : 1.6   : L/kg, 알로퓨리놀 분포용적
ke_oxy       : 0.033 : /h, 옥시퓨리놀 제거속도 (t½=21h)
k_form_oxy   : 0.60  : 알로퓨리놀→옥시퓨리놀 전환율
IC50_oxy     : 8.0   : mg/L, XO 50% 억제 옥시퓨리놀 농도
n_oxy        : 1.5   : XO 억제 힐계수

// --- Drug PK Parameters: Febuxostat ---
F_feb        : 0.84  : febuxostat 생체이용률
ka_feb       : 2.5   : /h, 흡수속도상수
CL_feb       : 5.2   : L/h
Vd_feb       : 1.2   : L/kg, 분포용적
ke_feb       : 0.12  : /h, 제거속도 (t½≈5.8h)
IC50_feb     : 0.001 : mg/L, XO 50% 억제 (non-purine, IC50=0.6nM)
n_feb        : 1.0   : 힐계수

// --- Drug PK Parameters: Uricosuric (generic) ---
F_uric       : 0.95  : 요산배설촉진제 생체이용률
ka_uric      : 2.0   : /h
CL_uric      : 15    : L/h
Vd_uric      : 0.8   : L/kg
ke_uric      : 0.14  : /h (t½≈5h, lesinurad)
IC50_URAT1   : 0.05  : mg/L, URAT1 50% 억제
k_uricosuric_eff : 0.35 : 요산배설 촉진 최대 효과 (35% 배설 증가)

// --- MSU Crystal Dynamics ---
k_nucl       : 0.001 : /day, 결정 핵화 속도 (>6.8mg/dL)
k_growth     : 0.005 : /(mg/dL*day), 결정 성장 속도
k_dissolve   : 0.003 : /day, 결정 용해 속도 (when UA<6.0)
SUA_sat      : 6.8   : mg/dL, MSU 포화 농도

// --- Endothelial / NO / BP ---
k_endo_dmg   : 0.005 : /day/(mg/dL), UA 직접 내피 손상 속도
k_endo_rep   : 0.02  : /day, 내피 회복 속도
NO_0         : 1.0   : 기저 NO 수준 (normalized)
k_NO_UA      : 0.008 : /(mg/dL), UA-NO 역관계 계수
k_BP_NO      : 5.0   : mmHg per unit NO decrease
BP_0         : 85    : mmHg, 기저 MAP

// --- GFR decline ---
k_GFR_UA     : 0.0008: /(mg/dL/day), UA 유발 GFR 감소 속도
k_GFR_BP     : 0.0005: /(mmHg/day), 고혈압 유발 GFR 감소
k_GFR_rep    : 0.001 : /day, 치료시 GFR 부분 회복

// --- Inflammation ---
k_IL1_MSU    : 0.02  : pg/mL/mg·day, MSU → IL-1β 유도
k_IL1_deg    : 0.30  : /day, IL-1β 소실
k_CRP_IL1    : 0.50  : mg/L/pg·mL/day, IL-1β → CRP
k_CRP_deg    : 0.15  : /day, CRP 소실

// --- Insulin Resistance ---
k_IR_UA      : 0.015 : /(mg/dL), UA → HOMA-IR 계수
k_IR_base    : 1.0   : 기저 HOMA-IR

// --- CV Risk ---
k_CV_UA      : 0.002 : SUA-CV 위험 계수
k_CV_BP      : 0.003 : BP-CV 위험 계수
k_CV_CRP     : 0.001 : CRP-CV 위험 계수
k_CV_decay   : 0.0001: /day, CV 위험 배경 진행

// --- Body weight / ABCG2 genetics ---
BW           : 70    : kg, 체중
ABCG2_Q141K  : 0     : ABCG2 Q141K 변이 (0=없음, 1=있음)
k_ABCG2_inh  : 0.50  : Q141K에 의한 ABCG2 기능 감소 (50%)

$INIT
UA_plasma    = 7.5   // mg/dL, 기저 고요산혈증
UA_tissue    = 500   // mg, 조직 분포
XO_free      = 1.0   // normalized
Oxypurinol   = 0     // mg/L
Febuxostat_C = 0     // mg/L
Uricosuric_C = 0     // mg/L
URAT1_free   = 1.0   // fraction
UrinaryUA    = 400   // mg/day (running pool)
MSU_depot    = 0     // mg
Endothelial_fn = 1.0 // normalized (1=normal)
NO_level     = 1.0   // normalized
BP           = 85    // MAP mmHg
GFR          = 90    // mL/min/1.73m²
IL1beta      = 5     // pg/mL
CRP          = 2.0   // mg/L
InsulinResist = 1.5  // HOMA-IR
CV_risk_score = 0    // cumulative
Tophus_vol   = 0     // mm³
ABCG2_frac   = 1.0   // functional fraction

$PARAM
// Dosing parameters (overridden in simulation)
DOSE_allo    = 0   // mg/day allopurinol
DOSE_feb     = 0   // mg/day febuxostat
DOSE_uric    = 0   // mg/day uricosuric

$ODE
// ---- 1. ABCG2 function (genetics) ----
double ABCG2_fn = ABCG2_frac * (1.0 - ABCG2_Q141K * k_ABCG2_inh);

// ---- 2. Purine Production (UA generation) ----
double UA_prod = (k_prod_base + k_dietary
                 + k_fructose * fructose_load * k_prod_base
                 + k_alcohol  * alcohol_use  * k_prod_base)
                * XO_free / 24.0; // mg/h

// ---- 3. XO inhibition by drugs ----
double inh_oxy = pow(Oxypurinol, n_oxy) / (pow(IC50_oxy, n_oxy) + pow(Oxypurinol, n_oxy));
double inh_feb = pow(Febuxostat_C, n_feb) / (pow(IC50_feb, n_feb) + pow(Febuxostat_C, n_feb));
double XO_inhibition = 1.0 - std::max(inh_oxy, inh_feb);
dxdt_XO_free = k_XO_syn - k_XO_deg * XO_free;
// XO effective activity modified by drugs
double XO_active = XO_free * XO_inhibition;

// ---- 4. Effective UA production ----
double UA_prod_eff = UA_prod * XO_active; // mg/h

// ---- 5. Renal excretion ----
double GFR_ratio = GFR / GFR_0;
// URAT1 inhibition by uricosuric
double inh_URAT1 = (Uricosuric_C > 0) ?
    Uricosuric_C / (IC50_URAT1 + Uricosuric_C) * k_uricosuric_eff : 0;
double URAT1_eff = k_URAT1 * (1.0 - inh_URAT1);
double ABCG2_ren_eff = k_ABCG2_ren * ABCG2_fn;
double FE_UA = 1.0 - URAT1_eff + k_OAT_secr + ABCG2_ren_eff; // net excretion fraction
FE_UA = std::min(FE_UA, 0.25); // cap at 25%
double k_renal_eff = k_renal_base * GFR_ratio * FE_UA / 12.0; // /h
double k_gut_eff   = k_gut_ABCG2 * ABCG2_fn * k_renal_base / 24.0; // /h

// ---- 6. UA plasma dynamics ----
// Distribution constant (plasma ↔ tissue)
double k_dist  = 0.5;  // /h plasma→tissue
double k_redist = 0.2; // /h tissue→plasma
dxdt_UA_plasma = UA_prod_eff / (0.6 * BW * 10.0)  // production → plasma (mg/dL/h)
                - k_renal_eff * UA_plasma           // renal excretion
                - k_gut_eff * UA_plasma             // gut secretion
                - k_dist * UA_plasma
                + k_redist * (UA_tissue / (0.4 * BW * 10.0)); // tissue redistribution

dxdt_UA_tissue = k_dist * UA_plasma * (0.6 * BW * 10.0)
                - k_redist * UA_tissue;

// ---- 7. Urinary UA ----
dxdt_UrinaryUA = k_renal_eff * UA_plasma * 24.0 - 0.01 * UrinaryUA; // mg/day dynamics

// ---- 8. Drug PK ----
// Allopurinol (simplified one-compartment; dose given via event)
// Oxypurinol: formed from allopurinol (modeled as input)
double k_oxy_form = (DOSE_allo > 0) ? DOSE_allo * F_allo * k_form_oxy / (Vd_allo * BW * 24.0) : 0;
dxdt_Oxypurinol   = k_oxy_form - ke_oxy * Oxypurinol;

// Febuxostat
double k_feb_in = (DOSE_feb > 0) ? DOSE_feb * F_feb / (Vd_feb * BW * 24.0) : 0;
dxdt_Febuxostat_C = k_feb_in - ke_feb * Febuxostat_C;

// Uricosuric
double k_uric_in = (DOSE_uric > 0) ? DOSE_uric * F_uric / (Vd_uric * BW * 24.0) : 0;
dxdt_Uricosuric_C = k_uric_in - ke_uric * Uricosuric_C;

// ---- 9. URAT1 fraction ----
dxdt_URAT1_free = 0.05 * (1.0 - URAT1_free) - 0.05 * inh_URAT1;

// ---- 10. MSU Crystal dynamics ----
double SUA_excess = std::max(UA_plasma - SUA_sat, 0.0); // mg/dL above saturation
double nucl_rate  = k_nucl  * SUA_excess * (MSU_depot < 1.0 ? 1.0 : 0.1);
double growth_rate = k_growth * SUA_excess * MSU_depot;
double dissolve_rate = k_dissolve * std::max(SUA_sat - UA_plasma, 0.0) * MSU_depot;
dxdt_MSU_depot = nucl_rate + growth_rate - dissolve_rate;

// Tophus: sustained crystal accumulation
dxdt_Tophus_vol = std::max(MSU_depot - 10.0, 0.0) * 0.01 - 0.001 * Tophus_vol;

// ---- 11. Endothelial function & NO ----
dxdt_Endothelial_fn = k_endo_rep * (1.0 - Endothelial_fn)
                     - k_endo_dmg * UA_plasma * Endothelial_fn;
dxdt_NO_level = 0.05 * (1.0 - NO_level)
               - k_NO_UA * (UA_plasma - 5.0) * std::max(UA_plasma - 5.0, 0.0) / 10.0;

// ---- 12. Blood pressure ----
double BP_delta_NO = k_BP_NO * (1.0 - NO_level);
dxdt_BP = 0.01 * (BP_0 + BP_delta_NO - BP);

// ---- 13. GFR decline ----
double GFR_dmg = k_GFR_UA * std::max(UA_plasma - 6.0, 0.0)
               + k_GFR_BP * std::max(BP - 90.0, 0.0);
double GFR_rep = (UA_plasma < 6.0) ? k_GFR_rep * (GFR_0 - GFR) : 0;
dxdt_GFR = GFR_rep - GFR_dmg * GFR;

// ---- 14. Inflammation ----
dxdt_IL1beta = k_IL1_MSU * MSU_depot - k_IL1_deg * IL1beta;
dxdt_CRP     = k_CRP_IL1 * IL1beta   - k_CRP_deg * CRP;

// ---- 15. Insulin resistance ----
dxdt_InsulinResist = 0.001 * (k_IR_base + k_IR_UA * UA_plasma - InsulinResist);

// ---- 16. CV risk score ----
dxdt_CV_risk_score = k_CV_UA  * std::max(UA_plasma - 6.0, 0.0)
                   + k_CV_BP  * std::max(BP - 90.0, 0.0)
                   + k_CV_CRP * std::max(CRP - 3.0, 0.0)
                   + k_CV_decay;

// ---- 17. ABCG2 fraction ----
dxdt_ABCG2_frac = 0; // genetic, constant

$TABLE
double SUA_mg_dL = UA_plasma;
double GoutFlareRisk = (UA_plasma > 9.0) ? 3.0 :
                       (UA_plasma > 8.0) ? 2.0 :
                       (UA_plasma > 7.0) ? 1.5 : 1.0;
double UricosuriaUA_g_day = UrinaryUA / 1000.0;
double eGFR_CKD_stage = (GFR >= 90) ? 1 : (GFR >= 60) ? 2 :
                         (GFR >= 30) ? 3 : (GFR >= 15) ? 4 : 5;
double XO_pct_inhibition = (1.0 - XO_active) * 100.0;
double MAP = BP;
double hsCRP = CRP;
double HOMA_IR = InsulinResist;
double Tophus_mm3 = Tophus_vol;
double CV_risk = CV_risk_score;

$CAPTURE SUA_mg_dL GFR MAP hsCRP HOMA_IR
         GoutFlareRisk XO_pct_inhibition Oxypurinol Febuxostat_C
         Uricosuric_C MSU_depot Tophus_mm3 IL1beta CRP CV_risk
         UrinaryUA UricosuriaUA_g_day
'

## Compile model
mod <- mcode("ahy_model", ahy_model_code, quiet = TRUE)

## ============================================================
## TREATMENT SCENARIOS (5 major scenarios)
## ============================================================

## Simulation parameters
sim_days <- 730   # 2 years
dt       <- 1/24  # hourly

run_sim <- function(model, params = list(), label = "Untreated") {
    mod_upd <- param(model, params)
    out <- mrgsim(mod_upd,
                  end   = sim_days,
                  delta = 1,          # daily output
                  digits = 4)
    df <- as.data.frame(out)
    df$scenario <- label
    df
}

## ---- Scenario 1: Untreated AHY (SUA=7.5mg/dL) ----
sc1 <- run_sim(mod, list(), label = "1. 미치료 AHY (SUA=7.5)")

## ---- Scenario 2: Allopurinol 300mg/day ----
sc2 <- run_sim(mod,
               list(DOSE_allo = 300),
               label = "2. 알로퓨리놀 300mg/day")

## ---- Scenario 3: Febuxostat 80mg/day ----
sc3 <- run_sim(mod,
               list(DOSE_feb = 80),
               label = "3. 페북소스타트 80mg/day")

## ---- Scenario 4: Lesinurad/Uricosuric 200mg + Allopurinol 300mg ----
sc4 <- run_sim(mod,
               list(DOSE_allo = 300, DOSE_uric = 200),
               label = "4. 알로퓨리놀+요산배설촉진제 병용")

## ---- Scenario 5: High Fructose Diet + AHY (no treatment) ----
sc5 <- run_sim(mod,
               list(fructose_load = 1.0, alcohol_use = 0.5),
               label = "5. 고과당식이+알코올 (생활습관 위험)")

## ---- Scenario 6: ABCG2 Q141K variant + Febuxostat 120mg ----
sc6 <- run_sim(mod,
               list(ABCG2_Q141K = 1, DOSE_feb = 120),
               label = "6. ABCG2 Q141K 변이+페북소스타트120")

## ---- Scenario 7: Treat-to-Target (<6mg/dL) Allopurinol + Lifestyle ----
sc7 <- run_sim(mod,
               list(DOSE_allo = 600,
                    fructose_load = 0,
                    alcohol_use   = 0),
               label = "7. 적극적 치료목표 SUA<6mg/dL")

all_sc <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6, sc7)

## ============================================================
## VISUALIZATION
## ============================================================

scenario_colors <- c(
  "1. 미치료 AHY (SUA=7.5)"                  = "#E53935",
  "2. 알로퓨리놀 300mg/day"                   = "#1E88E5",
  "3. 페북소스타트 80mg/day"                  = "#43A047",
  "4. 알로퓨리놀+요산배설촉진제 병용"         = "#8E24AA",
  "5. 고과당식이+알코올 (생활습관 위험)"      = "#F4511E",
  "6. ABCG2 Q141K 변이+페북소스타트120"       = "#00ACC1",
  "7. 적극적 치료목표 SUA<6mg/dL"             = "#039BE5"
)

p1 <- ggplot(all_sc, aes(time, SUA_mg_dL, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 6.8, linetype = "dashed", color = "red", size = 0.5) +
  geom_hline(yintercept = 6.0, linetype = "dotted", color = "blue", size = 0.5) +
  annotate("text", x = 600, y = 7.0, label = "MSU 포화점 (6.8)", size = 3, color = "red") +
  annotate("text", x = 600, y = 6.15, label = "치료 목표 (<6.0)", size = 3, color = "blue") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "혈청 요산 (SUA) 시뮬레이션",
       x = "시간 (일)", y = "SUA (mg/dL)", color = "시나리오") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

p2 <- ggplot(all_sc, aes(time, GFR, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "사구체 여과율 (eGFR)", x = "시간 (일)", y = "eGFR (mL/min/1.73m²)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p3 <- ggplot(all_sc, aes(time, MSU_depot, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "MSU 결정 침착량", x = "시간 (일)", y = "MSU 결정 (mg)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p4 <- ggplot(all_sc, aes(time, CV_risk, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "누적 심혈관 위험 점수", x = "시간 (일)", y = "CV Risk Score (rel.)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p5 <- ggplot(all_sc, aes(time, hsCRP, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 3.0, linetype = "dashed", color = "red") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "고감도 CRP (hs-CRP)", x = "시간 (일)", y = "hs-CRP (mg/L)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p6 <- ggplot(all_sc, aes(time, XO_pct_inhibition, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "XO 억제율 (%)", x = "시간 (일)", y = "XO 억제율 (%)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

combined_plot <- (p1 / (p2 + p3) / (p4 + p5 + p6)) +
  plot_annotation(
    title = "무증상 고요산혈증 (AHY) QSP 시뮬레이션 — 7개 치료 시나리오",
    subtitle = "2년간 시뮬레이션 | mrgsolve ODE 모델",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(combined_plot)

## ============================================================
## SENSITIVITY ANALYSIS: SUA target vs. GFR at 2 years
## ============================================================
cat("\n=== 2년 시뮬레이션 결과 요약 ===\n")
summary_table <- all_sc %>%
  filter(time == sim_days) %>%
  select(scenario, SUA_mg_dL, GFR, MSU_depot, CV_risk,
         hsCRP, HOMA_IR, Tophus_mm3, GoutFlareRisk) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
print(summary_table, width = 120)

## ============================================================
## DOSE-RESPONSE: Allopurinol dose vs. SUA at steady state
## ============================================================
cat("\n=== 알로퓨리놀 용량-반응 곡선 (정상상태 SUA) ===\n")
allo_doses <- c(0, 100, 200, 300, 400, 600, 900)
ss_SUA <- sapply(allo_doses, function(d) {
  sc <- run_sim(mod, list(DOSE_allo = d), label = as.character(d))
  tail(sc$SUA_mg_dL, 1)
})
dose_resp <- data.frame(
  Dose_mg   = allo_doses,
  SUA_mg_dL = round(ss_SUA, 2),
  pct_reduction = round((1 - ss_SUA / ss_SUA[1]) * 100, 1)
)
print(dose_resp)

## ============================================================
## PARAMETER REFERENCE TABLE
## ============================================================
param_ref <- data.frame(
  Parameter = c("기저 요산 생성", "신장 기저 배설", "XO IC50 (옥시퓨리놀)",
                "XO IC50 (페북소스타트)", "URAT1 재흡수 분획",
                "장관 ABCG2 분비 분획", "MSU 포화 농도", "GFR 기저치",
                "알로퓨리놀 생체이용률", "페북소스타트 생체이용률"),
  Value     = c("700 mg/day", "400 mg/day", "8.0 mg/L",
                "0.001 mg/L", "70%", "33%",
                "6.8 mg/dL", "90 mL/min/1.73m²",
                "90%", "84%"),
  Reference = c("Dalbeth 2019 NRP", "Reginato 2012 NRP",
                "Perez-Ruiz 2011 ARD", "Becker 2005 NEJM",
                "FitzGerald 2011 AR", "Mandal 2021 JCI",
                "Khanna 2012 ACR", "KDIGO 2012",
                "Hande 1984 AJM", "Schumacher 2008 NEJM")
)
cat("\n=== 주요 모델 파라미터 참고문헌 ===\n")
print(param_ref, row.names = FALSE)
