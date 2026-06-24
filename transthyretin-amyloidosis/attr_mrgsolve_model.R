##############################################################################
# Transthyretin Amyloidosis (ATTR) — Quantitative Systems Pharmacology Model
# mrgsolve ODE Model
#
# 질환: 트랜스티레틴(TTR) 아밀로이드증
#   · ATTRwt  (야생형 — 주로 심근병증, ≥60세)
#   · ATTRv   (변이형 — V30M 다발신경병증, V122I 심근병증 등)
#
# 주요 치료제 PK/PD 포함:
#   1. 타파미디스 (Tafamidis)   — TTR 사량체 동력학적 안정화
#   2. 파티시란  (Patisiran)    — siRNA, TTR mRNA↓80% (APOLLO NEJM 2018)
#   3. 뷔트리시란(Vutrisiran)   — GalNAc-siRNA, TTR↓87% (HELIOS-A NEJM 2022)
#   4. 이노테르센(Inotersen)    — ASO, TTR↓75% (NEURO-TTR NEJM 2018)
#
# 참조 임상시험:
#   · ATTR-ACT   (Maurer et al., NEJM 2018)  — 타파미디스 심장 ATTR
#   · APOLLO     (Adams et al., NEJM 2018)   — 파티시란 FAP
#   · HELIOS-A   (Solomon et al., NEJM 2022) — 뷔트리시란 FAP
#   · NEURO-TTR  (Benson et al., NEJM 2018)  — 이노테르센 FAP
#
# ODE 구획 (20개):
#   TTR PK: TTR_C, TTR_P (2구획)
#   응집 경로: TTR_MF, TTR_OL (2구획)
#   아밀로이드 침착: AMY_H, AMY_N, AMY_GI (3구획)
#   심장 병리: LV_THICK, BNP_C, TROP_C, CARD_FUNC (4구획)
#   신경 병리: NIS_TOT, AUTO_NP (2구획)
#   약물 PK: TAF_GUT, TAF_C, PAT_EFF, VUT_EFF, INO_EFF (5구획)
#   복합 결과: SixMWT, KCCQ_IDX (2구획)
#   총 20구획
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

##############################################################################
# mrgsolve 모델 정의
##############################################################################
code_attr <- '
$PROB
// Transthyretin Amyloidosis (ATTR) — QSP Model
// 20-compartment ODE system
// Units: time=days, concentration=ug/mL (or normalized AU)

$PARAM
// ============================================================
// TTR 생산·소실 파라미터
// ============================================================
kprod_TTR  = 104.0   // ug/mL/day: 간 TTR 생산율 (Css=300→kel×300=104)
kel_TTR    = 0.347   // /day: TTR 일차 소실 (t1/2=2d)
k12_TTR    = 0.1     // /day: 중심→말초 분포
k21_TTR    = 0.05    // /day: 말초→중심 분포

// ============================================================
// TTR 응집·아밀로이드 형성 파라미터
// ============================================================
// 정상(mut_f=1.0): 야생형, V30M·V122I = 더 높은 불안정성
mut_factor = 1.0     // 1.0=WT, 1.5=Val122I, 2.0=Val30M (더 불안정)
kdiss_base = 5e-4    // /day: 기저 사량체 해리 속도 (야생형)
kconf      = 0.05    // /day: 단량체 → 잘못폴딩 전환율
kolig      = 0.10    // /day: 잘못폴딩 → 올리고머
kfib       = 0.02    // /day: 올리고머 → 피브릴
frac_heart = 0.45    // 피브릴 심장 침착 분율
frac_nerve = 0.35    // 피브릴 신경 침착 분율
frac_GI    = 0.20    // 피브릴 GI 침착 분율
kclear_amyl = 0.003  // /day: 아밀로이드 내인성 제거율 (매우 느림)

// ============================================================
// 심장 아밀로이드 → 임상 파라미터
// ============================================================
LV_base    = 10.0    // mm: 기저 LV 벽두께 (정상 8-12mm)
LV_kgrowth = 0.8     // mm/AU: 아밀로이드/LV비후 스케일
BNP_base   = 80.0    // pg/mL: 기저 NT-proBNP
BNP_scale  = 15.0    // pg/mL/AU: 아밀로이드→BNP 증가
TROP_base  = 5.0     // ng/L: 기저 트로포닌T
TROP_scale = 0.8     // ng/L/AU: 아밀로이드→트로포닌 증가
CARD_base  = 1.0     // AU: 기저 심장 기능 (1=정상)
CARD_decline = 0.15  // /AU: 아밀로이드→기능 저하율

// ============================================================
// 신경 아밀로이드 → 임상 파라미터
// ============================================================
NIS_base   = 10.0    // 기저 NIS 점수 (정상 0, 진단시 평균 ~80)
NIS_scale  = 5.0     // NIS/AU: 신경 아밀로이드→NIS 증가

// ============================================================
// 타파미디스 PK 파라미터 (1구획 경구)
// ============================================================
ka_TAF   = 0.347    // /h → /day: 흡수속도 (Tmax~4h)
CL_TAF   = 7.1      // L/day: 겉보기 청소율 (t1/2=59h = 2.46d → CL=Vd*ke)
Vd_TAF   = 18.0     // L: 분포용적
F_TAF    = 1.0      // 경구 생체이용률 (~100%)

// 타파미디스 PD (Emax 모델 — 사량체 해리 억제)
EC50_TAF = 2.0      // ug/mL: 해리 억제 EC50
Emax_TAF = 0.80     // 최대 억제율 (80%)

// ============================================================
// siRNA/ASO 효과 파라미터 (효과 구획 모델)
// ============================================================
// 각 약물: 간 내 효과 구획에서 서서히 TTR 생산 억제
// 파티시란 Q3W IV: 최대 억제 80%, ke_eff≈0.05/day
kout_PAT = 0.035    // /day: 파티시란 효과 소실 속도 (지속 효과)
Emax_PAT = 0.80     // 최대 억제율
EC50_PAT = 0.5      // 정규화된 효과 구획 농도

// 뷔트리시란 Q3M SC: 최대 억제 87%
kout_VUT = 0.025    // /day: 더 긴 지속성 (Q3M 투여 정당화)
Emax_VUT = 0.87
EC50_VUT = 0.5

// 이노테르센 QW SC: 최대 억제 75%
kout_INO = 0.045    // /day: 약간 빠른 소실
Emax_INO = 0.75
EC50_INO = 0.5

// ============================================================
// 6분보행·삶의 질 파라미터
// ============================================================
SixMWT_base = 420.0 // m: 기저 6분보행거리 (ATTRwt 진단시 평균)
SixMWT_k    = 0.5   // m/AU: 아밀로이드 증가에 따른 거리 감소율
KCCQ_base   = 70.0  // 기저 KCCQ-OS (0-100)
KCCQ_k      = 2.5   // KCCQ/AU: 아밀로이드→삶의 질 저하

$FIXED
Vd_TAF_fixed = 18.0

$CMT
// === TTR PK (2구획) ===
TTR_C         // 혈장 TTR 농도 (ug/mL)
TTR_P         // 말초 TTR (ug/mL×L 단위)

// === 응집 경로 (2구획) ===
TTR_MF        // 잘못폴딩된 단량체 (AU)
TTR_OL        // 가용성 올리고머 (AU) - 세포독성!

// === 아밀로이드 침착 (3구획) ===
AMY_H         // 심장 아밀로이드 (g, 상대 단위)
AMY_N         // 신경 아밀로이드 (AU)
AMY_GI        // 위장관 아밀로이드 (AU)

// === 심장 병리 출력 (4구획) ===
LV_THICK      // LV 벽두께 (mm)
BNP_C         // NT-proBNP (pg/mL)
TROP_C        // 고감도 트로포닌T (ng/L)
CARD_FUNC     // 심장 기능 지수 (1=정상, 0=중증 기능저하)

// === 신경 병리 (2구획) ===
NIS_TOT       // NIS 총점
AUTO_NP       // 자율신경 병증 점수 (0-1)

// === 약물 PK (5구획) ===
TAF_GUT       // 타파미디스 위장관 (mg)
TAF_C         // 타파미디스 혈장 (ug/mL)
PAT_EFF       // 파티시란 효과 구획 (0-1 정규화)
VUT_EFF       // 뷔트리시란 효과 구획 (0-1 정규화)
INO_EFF       // 이노테르센 효과 구획 (0-1 정규화)

// === 복합 결과 (2구획) ===
SixMWT        // 6분보행거리 (m)
KCCQ_IDX      // KCCQ-OS 지수 (0-100)

$INIT
TTR_C   = 300.0   // 정상 혈장 TTR (μg/mL = 30 mg/dL)
TTR_P   = 150.0   // 말초 분포 초기값
TTR_MF  = 0.001   // 매우 낮은 기저 잘못폴딩
TTR_OL  = 0.0001  // 기저 올리고머 (극미량)
AMY_H   = 0.0     // 아밀로이드 축적 없음 (진단 전 또는 정상)
AMY_N   = 0.0
AMY_GI  = 0.0
LV_THICK  = 10.0  // 정상 LV 벽두께 (mm)
BNP_C     = 80.0  // 정상 NT-proBNP (pg/mL)
TROP_C    = 5.0   // 정상 트로포닌
CARD_FUNC = 1.0   // 정상 심장 기능
NIS_TOT   = 10.0  // 경미한 기저 NIS
AUTO_NP   = 0.05  // 거의 정상 자율신경
TAF_GUT   = 0.0
TAF_C     = 0.0
PAT_EFF   = 0.0
VUT_EFF   = 0.0
INO_EFF   = 0.0
SixMWT    = 420.0 // 정상 6분보행거리
KCCQ_IDX  = 70.0  // 기저 KCCQ-OS

$MAIN
// ============================================================
// 타파미디스 약동학 파라미터 계산
// ============================================================
double ke_TAF = CL_TAF / Vd_TAF;  // /day: 소실 속도
double ka_TAF_day = ka_TAF * 24.0; // /h → /day 변환

// ============================================================
// 타파미디스 PD: Emax 모델 (사량체 해리 억제)
// ============================================================
double Inh_TAF = (Emax_TAF * TAF_C) / (EC50_TAF + TAF_C + 1e-10);

// ============================================================
// 유효 사량체 해리 속도 (돌연변이 인자 × 약물 억제)
// ============================================================
double kdiss_eff = kdiss_base * mut_factor * (1.0 - Inh_TAF);

// ============================================================
// siRNA/ASO TTR 생산 억제율 계산 (가산 억제)
// ============================================================
double Inh_PAT = (Emax_PAT * PAT_EFF) / (EC50_PAT + PAT_EFF + 1e-10);
double Inh_VUT = (Emax_VUT * VUT_EFF) / (EC50_VUT + VUT_EFF + 1e-10);
double Inh_INO = (Emax_INO * INO_EFF) / (EC50_INO + INO_EFF + 1e-10);
// 결합 억제 (Bliss 독립성 모델)
double Inh_total = 1.0 - (1.0 - Inh_PAT)*(1.0 - Inh_VUT)*(1.0 - Inh_INO);
double kprod_eff = kprod_TTR * (1.0 - Inh_total);

// ============================================================
// 응집 경로에서 분리된 TTR 공급량 추정
// ============================================================
// 사량체 해리 → 단량체 → 잘못폴딩 경로의 구동력
double TTR_avail = (TTR_C > 0) ? TTR_C : 0.0;
double rate_dissoc = kdiss_eff * TTR_avail;

$ODE
// ============================================================
// TTR PK (2구획 모델)
// ============================================================
dxdt_TTR_C = kprod_eff - kel_TTR * TTR_C
             - k12_TTR * TTR_C + k21_TTR * TTR_P
             - rate_dissoc;  // 해리로 인한 사량체 손실

dxdt_TTR_P = k12_TTR * TTR_C - k21_TTR * TTR_P;

// ============================================================
// 응집 연쇄반응 (TTR 미스폴딩 → 올리고머)
// ============================================================
// 잘못폴딩 단량체 생성: 사량체 해리 후 구조 변형
double input_MF = kconf * rate_dissoc;
dxdt_TTR_MF = input_MF - kolig * TTR_MF;

// 올리고머 형성: 잘못폴딩에서 핵형성
double input_OL = kolig * TTR_MF;
dxdt_TTR_OL = input_OL - kfib * TTR_OL;

// ============================================================
// 아밀로이드 장기 침착 (피브릴 → 침착)
// ============================================================
double rate_fib = kfib * TTR_OL;  // 피브릴 생성 속도
dxdt_AMY_H = frac_heart * rate_fib - kclear_amyl * AMY_H;
dxdt_AMY_N = frac_nerve * rate_fib - kclear_amyl * AMY_N;
dxdt_AMY_GI = frac_GI  * rate_fib - kclear_amyl * AMY_GI;

// ============================================================
// 심장 병리 동역학 (아밀로이드 부담 → 임상 지표)
// ============================================================
// LV 벽두께: 아밀로이드 침착에 비례하여 증가
double LV_target = LV_base + LV_kgrowth * AMY_H;
dxdt_LV_THICK = 0.05 * (LV_target - LV_THICK);  // 일차 반응 속도

// NT-proBNP: 충만압 상승 반영
double BNP_target = BNP_base + BNP_scale * AMY_H
                    + 50.0 * TTR_OL;  // 올리고머 직접독성 기여
dxdt_BNP_C = 0.1 * (BNP_target - BNP_C);

// 고감도 트로포닌T: 심근 손상 반영
double TROP_target = TROP_base + TROP_scale * AMY_H
                     + 0.5 * TTR_OL;
dxdt_TROP_C = 0.1 * (TROP_target - TROP_C);

// 심장 기능 지수: 아밀로이드 증가 시 감소 (LVEF 보존→저하)
double CARD_target = CARD_base - CARD_decline * AMY_H;
if(CARD_target < 0.1) CARD_target = 0.1;  // 최소값 설정
dxdt_CARD_FUNC = 0.02 * (CARD_target - CARD_FUNC);

// ============================================================
// 신경 병리 동역학
// ============================================================
double NIS_target = NIS_base + NIS_scale * AMY_N;
dxdt_NIS_TOT = 0.05 * (NIS_target - NIS_TOT);

double AUTO_target = 0.05 + 0.3 * AMY_N;
if(AUTO_target > 1.0) AUTO_target = 1.0;
dxdt_AUTO_NP = 0.05 * (AUTO_target - AUTO_NP);

// ============================================================
// 타파미디스 PK (1구획 경구 모델)
// ============================================================
// 투여: 61mg QD oral (dose=61 mg → TAF_GUT에 투여)
dxdt_TAF_GUT = -ka_TAF_day * TAF_GUT;
dxdt_TAF_C   = (ka_TAF_day * TAF_GUT * F_TAF) / Vd_TAF
               - ke_TAF * TAF_C;

// ============================================================
// siRNA/ASO 효과 구획 (간 축적/효과 지속 모델)
// ============================================================
// kin_PAT 등: 투여 이벤트에서 직접 ADD로 증가
// kout_PAT: 서서히 효과 소실
dxdt_PAT_EFF = -kout_PAT * PAT_EFF;   // INPUT은 $ODE event로 처리
dxdt_VUT_EFF = -kout_VUT * VUT_EFF;
dxdt_INO_EFF = -kout_INO * INO_EFF;

// ============================================================
// 복합 결과 지표
// ============================================================
// 6분보행거리: 심장 기능 및 심장 아밀로이드 부담 반영
double SixMWT_target = SixMWT_base * CARD_FUNC
                       - SixMWT_k * AMY_H;
if(SixMWT_target < 50) SixMWT_target = 50;  // 최소값
dxdt_SixMWT = 0.05 * (SixMWT_target - SixMWT);

// KCCQ-OS: 심장 기능 및 증상에 따라 변화
double KCCQ_target = KCCQ_base * CARD_FUNC * 1.1
                     - KCCQ_k * AMY_H;
if(KCCQ_target < 0) KCCQ_target = 0;
if(KCCQ_target > 100) KCCQ_target = 100;
dxdt_KCCQ_IDX = 0.05 * (KCCQ_target - KCCQ_IDX);

$TABLE
// ============================================================
// 유도 변수 (표에 출력할 임상 지표들)
// ============================================================
capture TTR_plasma     = TTR_C;          // 혈장 TTR (ug/mL)
capture TTR_mg_dL      = TTR_C / 10.0;  // 혈장 TTR (mg/dL)
capture Agg_misfold    = TTR_MF;         // 잘못폴딩 단량체 (AU)
capture Oligomers      = TTR_OL;         // 올리고머 (AU)
capture CardiacAmyloid = AMY_H;          // 심장 아밀로이드 (AU)
capture NerveAmyloid   = AMY_N;          // 신경 아밀로이드 (AU)
capture LV_wall_mm     = LV_THICK;       // LV 벽두께 (mm)
capture NT_proBNP      = BNP_C;          // NT-proBNP (pg/mL)
capture TroponinT_ngL  = TROP_C;         // hsTnT (ng/L)
capture CardFunc_idx   = CARD_FUNC;      // 심장 기능 (1=정상)
capture NIS_score      = NIS_TOT;        // NIS 총점
capture Walk_6min_m    = SixMWT;         // 6MWT 거리 (m)
capture KCCQ_score     = KCCQ_IDX;       // KCCQ-OS (0-100)
capture Tafamidis_ug   = TAF_C;          // 타파미디스 혈장 (ug/mL)
capture Inh_TAF_frac   = (Emax_TAF * TAF_C) / (EC50_TAF + TAF_C + 1e-10);
capture Inh_siRNA_frac = 1.0 - (1.0 - (Emax_PAT*PAT_EFF)/(EC50_PAT+PAT_EFF+1e-10))
                              *(1.0 - (Emax_VUT*VUT_EFF)/(EC50_VUT+VUT_EFF+1e-10))
                              *(1.0 - (Emax_INO*INO_EFF)/(EC50_INO+INO_EFF+1e-10));

// NAC 병기 분류 (임상 편의용 출력)
// 1: BNP<3000 & TropT<50 ng/L
// 2: 하나 초과
// 3: 둘 다 초과
capture NAC_stage = (NT_proBNP >= 3000 && TroponinT_ngL >= 50) ? 3.0 :
                    (NT_proBNP >= 3000 || TroponinT_ngL >= 50) ? 2.0 : 1.0;

$CAPTURE TTR_plasma TTR_mg_dL Agg_misfold Oligomers CardiacAmyloid NerveAmyloid
         LV_wall_mm NT_proBNP TroponinT_ngL CardFunc_idx NIS_score
         Walk_6min_m KCCQ_score Tafamidis_ug Inh_TAF_frac Inh_siRNA_frac NAC_stage
'

##############################################################################
# 모델 컴파일
##############################################################################
mod_attr <- mcode("ATTR_QSP", code_attr)
mod_attr <- mod_attr %>% param(ka_TAF = 0.347 * 24)

##############################################################################
# 헬퍼 함수: siRNA 투여 일정 생성
##############################################################################
gen_siRNA_events <- function(n_doses, interval_days, eff_bolus = 1.0,
                             cmt_name = "PAT_EFF") {
  data.frame(
    time  = seq(0, by = interval_days, length.out = n_doses),
    amt   = eff_bolus,
    cmt   = cmt_name,
    evid  = 1,
    rate  = 0
  )
}

##############################################################################
# 시뮬레이션 1: 자연 경과 (ATTRwt 심근병증, 무치료)
# 참조: Ruberg et al., Circ Heart Fail 2019
##############################################################################
sim_natural <- function(mod, end_days = 1825) {
  mod %>%
    param(mut_factor = 1.0) %>%
    mrgsim(end = end_days, delta = 7) %>%
    as.data.frame() %>%
    mutate(scenario = "1. 자연 경과 (무치료)")
}

##############################################################################
# 시뮬레이션 2: 타파미디스 61mg QD (ATTR-ACT 설계)
# 참조: Maurer et al., NEJM 2018 — 30개월 추적
##############################################################################
sim_tafamidis <- function(mod, end_days = 1825) {
  # 타파미디스 61mg QD 투여 일정
  taf_dosing <- ev(
    time = 0, amt = 61, cmt = "TAF_GUT", evid = 1, rate = 0,
    ii = 1, addl = end_days - 1
  )
  mod %>%
    param(mut_factor = 1.0) %>%
    mrgsim(events = taf_dosing, end = end_days, delta = 7) %>%
    as.data.frame() %>%
    mutate(scenario = "2. 타파미디스 61mg QD (ATTR-ACT)")
}

##############################################################################
# 시뮬레이션 3: 파티시란 0.3mg/kg Q3W IV (APOLLO 설계, V30M FAP)
# 참조: Adams et al., NEJM 2018 — 18개월 추적
##############################################################################
sim_patisiran <- function(mod, end_days = 540, body_wt = 70) {
  n_doses <- floor(end_days / 21) + 1
  pat_events <- data.frame(
    time  = seq(0, by = 21, length.out = n_doses),
    amt   = 1.0,          # 효과 구획 정규화 입력
    cmt   = "PAT_EFF",
    evid  = 1,
    rate  = 0,
    ID    = 1
  )
  mod %>%
    param(mut_factor = 2.0) %>%  # V30M: 더 불안정한 사량체
    mrgsim_df(data = pat_events, end = end_days, delta = 7,
              add = list(ID = 1)) %>%
    mutate(scenario = "3. 파티시란 Q3W IV (APOLLO)")
}

##############################################################################
# 시뮬레이션 4: 뷔트리시란 25mg Q3M SC (HELIOS-A 설계)
# 참조: Solomon et al., NEJM 2022 — 9개월 추적 (모델: 24개월)
##############################################################################
sim_vutrisiran <- function(mod, end_days = 720) {
  n_doses <- floor(end_days / 91) + 1
  vut_events <- data.frame(
    time  = seq(0, by = 91, length.out = n_doses),
    amt   = 1.0,
    cmt   = "VUT_EFF",
    evid  = 1,
    rate  = 0,
    ID    = 1
  )
  mod %>%
    param(mut_factor = 2.0) %>%
    mrgsim_df(data = vut_events, end = end_days, delta = 7,
              add = list(ID = 1)) %>%
    mutate(scenario = "4. 뷔트리시란 Q3M SC (HELIOS-A)")
}

##############################################################################
# 시뮬레이션 5: 이노테르센 284mg QW SC (NEURO-TTR 설계)
# 참조: Benson et al., NEJM 2018 — 15개월 추적
##############################################################################
sim_inotersen <- function(mod, end_days = 540) {
  n_doses <- floor(end_days / 7) + 1
  ino_events <- data.frame(
    time  = seq(0, by = 7, length.out = n_doses),
    amt   = 1.0,
    cmt   = "INO_EFF",
    evid  = 1,
    rate  = 0,
    ID    = 1
  )
  mod %>%
    param(mut_factor = 2.0) %>%
    mrgsim_df(data = ino_events, end = end_days, delta = 7,
              add = list(ID = 1)) %>%
    mutate(scenario = "5. 이노테르센 QW SC (NEURO-TTR)")
}

##############################################################################
# 시뮬레이션 6: 타파미디스 + 적극적 HF 관리 병용
# (이뇨제·ACEi/ARB 등 증상 치료 병행 모델)
##############################################################################
sim_combo <- function(mod, end_days = 1825) {
  taf_dosing <- ev(
    time = 0, amt = 61, cmt = "TAF_GUT", evid = 1, rate = 0,
    ii = 1, addl = end_days - 1
  )
  # 적극적 HF 관리: BNP_scale 파라미터 감소 (이뇨제 효과 간접 반영)
  mod %>%
    param(mut_factor = 1.0, BNP_scale = 10.0) %>%
    mrgsim(events = taf_dosing, end = end_days, delta = 7) %>%
    as.data.frame() %>%
    mutate(scenario = "6. 타파미디스 + 적극적 HF 관리 병용")
}

##############################################################################
# 전체 시나리오 실행 및 결합
##############################################################################
run_all_scenarios <- function(mod) {
  message(">> 시나리오 1: 자연 경과 시뮬레이션 중...")
  s1 <- sim_natural(mod)

  message(">> 시나리오 2: 타파미디스 61mg QD 시뮬레이션 중...")
  s2 <- sim_tafamidis(mod)

  message(">> 시나리오 3: 파티시란 Q3W IV 시뮬레이션 중...")
  s3 <- tryCatch(sim_patisiran(mod), error = function(e) {
    message("파티시란 시뮬 오류: ", e$message)
    NULL
  })

  message(">> 시나리오 4: 뷔트리시란 Q3M SC 시뮬레이션 중...")
  s4 <- tryCatch(sim_vutrisiran(mod), error = function(e) {
    message("뷔트리시란 시뮬 오류: ", e$message)
    NULL
  })

  message(">> 시나리오 5: 이노테르센 QW SC 시뮬레이션 중...")
  s5 <- tryCatch(sim_inotersen(mod), error = function(e) {
    message("이노테르센 시뮬 오류: ", e$message)
    NULL
  })

  message(">> 시나리오 6: 타파미디스 + HF 관리 병용...")
  s6 <- sim_combo(mod)

  # 성공한 시뮬레이션만 결합
  result_list <- list(s1, s2, s3, s4, s5, s6)
  result_list <- result_list[!sapply(result_list, is.null)]

  bind_rows(result_list) %>%
    mutate(time_years = time / 365.25)
}

##############################################################################
# 시각화 함수 (주요 결과 그림)
##############################################################################
plot_attr_results <- function(results) {
  # 색상 팔레트
  pal <- c(
    "1. 자연 경과 (무치료)"           = "#E74C3C",
    "2. 타파미디스 61mg QD (ATTR-ACT)" = "#27AE60",
    "3. 파티시란 Q3W IV (APOLLO)"      = "#2980B9",
    "4. 뷔트리시란 Q3M SC (HELIOS-A)"  = "#8E44AD",
    "5. 이노테르센 QW SC (NEURO-TTR)"  = "#D35400",
    "6. 타파미디스 + 적극적 HF 관리 병용" = "#16A085"
  )

  # 1) 혈장 TTR 농도
  p1 <- ggplot(results, aes(time_years, TTR_mg_dL, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal) +
    labs(title = "혈장 TTR 농도", x = "시간 (년)", y = "TTR (mg/dL)",
         color = "치료 시나리오") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  # 2) 심장 아밀로이드 부담
  p2 <- ggplot(results, aes(time_years, CardiacAmyloid, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal) +
    labs(title = "심장 아밀로이드 부담", x = "시간 (년)",
         y = "아밀로이드 (AU)", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  # 3) NT-proBNP
  p3 <- ggplot(results, aes(time_years, NT_proBNP, color = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 3000, linetype = "dashed", color = "gray40") +
    annotate("text", x = 0.2, y = 3200, label = "NAC 3기 경계 (3000 pg/mL)",
             size = 3, hjust = 0) +
    scale_color_manual(values = pal) +
    labs(title = "NT-proBNP", x = "시간 (년)", y = "NT-proBNP (pg/mL)",
         color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  # 4) 6분보행거리
  p4 <- ggplot(results, aes(time_years, Walk_6min_m, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal) +
    labs(title = "6분보행거리 (6MWT)", x = "시간 (년)", y = "거리 (m)",
         color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  # 5) LV 벽두께
  p5 <- ggplot(results, aes(time_years, LV_wall_mm, color = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 12, linetype = "dotted", color = "blue") +
    annotate("text", x = 0.1, y = 12.3, label = "정상 상한 (12mm)",
             size = 3, color = "blue", hjust = 0) +
    scale_color_manual(values = pal) +
    labs(title = "LV 벽두께", x = "시간 (년)", y = "두께 (mm)",
         color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  # 6) KCCQ 점수 (삶의 질)
  p6 <- ggplot(results, aes(time_years, KCCQ_score, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = "KCCQ-OS (삶의 질)", x = "시간 (년)", y = "점수 (0-100)",
         color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "none")

  (p1 + p2) / (p3 + p4) / (p5 + p6) +
    plot_annotation(
      title = "ATTR 아밀로이드증 QSP 모델 — 치료 시나리오 비교",
      subtitle = paste0("시나리오 1: 자연 경과 (ATTRwt 무치료) | 시나리오 2-6: 치료 효과 비교\n",
                        "참조: ATTR-ACT·APOLLO·HELIOS-A·NEURO-TTR 임상 데이터"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 9))
    )
}

##############################################################################
# 파라미터 민감도 분석 (타파미디스 Emax·EC50 효과)
##############################################################################
run_sensitivity_taf <- function(mod) {
  ec50_vals <- c(0.5, 1.0, 2.0, 4.0, 8.0)
  emax_vals <- c(0.60, 0.70, 0.80, 0.90)

  results_list <- list()
  for (ec50 in ec50_vals) {
    taf_ev <- ev(time = 0, amt = 61, cmt = "TAF_GUT", evid = 1,
                 rate = 0, ii = 1, addl = 1824)
    out <- mod %>%
      param(EC50_TAF = ec50, mut_factor = 1.0) %>%
      mrgsim(events = taf_ev, end = 1825, delta = 7) %>%
      as.data.frame() %>%
      mutate(EC50 = ec50)
    results_list[[length(results_list) + 1]] <- out
  }
  bind_rows(results_list)
}

##############################################################################
# 가상 환자 집단 시뮬레이션 (이분산성 분석)
##############################################################################
run_virtual_population <- function(mod, n_patients = 100, end_days = 1825) {
  set.seed(42)
  # 환자 특성 샘플링 (로그정규 분포 가정)
  pop_params <- data.frame(
    ID         = 1:n_patients,
    kprod_TTR  = rlnorm(n_patients, log(104),  0.15),
    kel_TTR    = rlnorm(n_patients, log(0.347), 0.12),
    mut_factor = rlnorm(n_patients, log(1.0),   0.10),
    kdiss_base = rlnorm(n_patients, log(5e-4),  0.30),
    kfib       = rlnorm(n_patients, log(0.02),  0.25)
  )

  # 타파미디스 투여 데이터셋
  taf_data <- do.call(rbind, lapply(1:n_patients, function(i) {
    data.frame(
      ID = i, time = 0, amt = 61, cmt = "TAF_GUT",
      evid = 1, rate = 0, ii = 1, addl = end_days - 1
    )
  }))

  # 시뮬레이션
  out <- mod %>%
    data_set(taf_data) %>%
    idata_set(pop_params) %>%
    mrgsim(end = end_days, delta = 30) %>%
    as.data.frame()

  out
}

##############################################################################
# 실행 예시 (직접 실행 시)
##############################################################################
if (interactive()) {
  # 모델 로드 확인
  cat("=== ATTR QSP 모델 파라미터 ===\n")
  print(param(mod_attr))
  cat("\n=== ATTR QSP 모델 초기 상태 ===\n")
  print(init(mod_attr))

  # 모든 치료 시나리오 실행
  cat("\n=== 치료 시나리오 시뮬레이션 시작 ===\n")
  results <- run_all_scenarios(mod_attr)

  # 5년 시점 요약
  cat("\n=== 5년 시점 주요 결과 요약 ===\n")
  summary_5y <- results %>%
    filter(abs(time - 1825) < 14) %>%  # 5년 시점 (±2주)
    group_by(scenario) %>%
    slice(1) %>%
    ungroup() %>%
    select(scenario, TTR_mg_dL, CardiacAmyloid, NT_proBNP,
           Walk_6min_m, KCCQ_score, NAC_stage) %>%
    as.data.frame()
  print(summary_5y)

  # 그래프 출력
  p_final <- plot_attr_results(results)
  print(p_final)

  # 민감도 분석
  cat("\n=== 타파미디스 EC50 민감도 분석 ===\n")
  sens_results <- run_sensitivity_taf(mod_attr)
  p_sens <- ggplot(
    sens_results %>% filter(time %% 30 < 8),
    aes(time / 365.25, CardiacAmyloid, color = factor(EC50))
  ) +
    geom_line() +
    labs(title = "타파미디스 EC50 민감도", x = "시간 (년)",
         y = "심장 아밀로이드 (AU)", color = "EC50 (μg/mL)") +
    theme_bw()
  print(p_sens)

  # 가상 환자 집단 시뮬레이션
  cat("\n=== 가상 환자 집단 시뮬레이션 (N=100) ===\n")
  pop_results <- run_virtual_population(mod_attr, n_patients = 50)
  p_pop <- ggplot(pop_results, aes(time / 365.25, CardiacAmyloid, group = ID)) +
    geom_line(alpha = 0.3, color = "steelblue") +
    stat_summary(aes(group = 1), fun = median, geom = "line",
                 color = "red", linewidth = 1.2) +
    labs(title = "가상 환자 집단 — 타파미디스 심장 아밀로이드 반응",
         subtitle = "파란 선: 개별 환자 (N=50) | 빨간 선: 중앙값",
         x = "시간 (년)", y = "심장 아밀로이드 (AU)") +
    theme_bw()
  print(p_pop)
}
