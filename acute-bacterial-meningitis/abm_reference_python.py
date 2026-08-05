#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
급성 세균성 수막염 (Acute Bacterial Meningitis, pneumococcal) — QSP 참조 구현
===========================================================================

이 파일은 abm_mrgsolve_model.R 의 63개 ODE 를 **의존성 없는 순수 Python RK4** 로
다시 구현한 것이다.  이 환경에는 R 런타임이 없으므로, mrgsolve 모델의 방정식이
실제로 물리적·약리학적으로 성립하는지를 여기서 수치적으로 검증하고 그 결과를
abm_reference_output.txt 로 남긴다.  두 파일의 파라미터와 방정식은 1:1 대응한다.

-------------------------------------------------------------------------------
모델이 계산으로 주장하는 것
-------------------------------------------------------------------------------
(1) 손상 flux 는 균 수가 아니라 **죽이는 속도 × 균 수 × 한 마리당 염증 수확량**
    이라는 곱이다.  이 곱은 첫 항생제 투여 직후에 최대가 된다 (N 은 그 순간
    가장 크고 k_kill 은 0 → Emax 로 점프한다).  그래서 스테로이드는 "언젠가"
    주는 약이 아니라 그 봉우리보다 **먼저** 켜져 있어야 하는 방패다.

(2) 친수성 항생제의 CSF 침투는 염증(장벽 투과성 Pb)의 증가함수인데, 덱사메타손은
    염증을 끈다.  즉 **덱사메타손은 항생제가 들어갈 문을 함께 닫는다.**
    세프트리악손은 C/MIC 여유가 100배 이상이라 순효과가 이득이지만, 반코마이신
    (MIC 1 mg/L · 침투 2-13 %)은 여유가 없어 같은 조작이 살균속도를 반으로
    깎는다 — 부호가 약에 따라 갈린다.

(3) β-락탐은 용해성 살균이라 한 마리당 수확량 Y 가 크고 리팜핀은 비용해성이라
    작다.  **같은 log-kill 에서 손상 적분이 다르다.**

-------------------------------------------------------------------------------
이 수치검증이 잡아낸 실제 결함 (v1 → v2)
-------------------------------------------------------------------------------
[F1]  포도당 소비항에 기질 의존성이 없어, CSF 포도당이 t=2 h 에 0 으로 붕괴한
      뒤에도 호중구·세균이 계속 "없는 포도당"을 먹었다.  그 결과 락테이트가
      243 mmol/L (실측 6-12) 까지 올라갔다.  → 모든 소비항에 Michaelis-Menten
      (Km_use = 5 mg/dL) 를 곱했다.  락테이트는 이제 포도당 유입속도에
      상한이 묶인다 (락테이트는 포도당에서만 나오므로 이것이 물리적 진실이다).
[F2]  호중구 유입상수 k_influx = 1.4e4 가 CSF 백혈구를 121,000/µL 로 만들었다
      (폐렴균 수막염 실측 1,000-5,000/µL).  → 500 으로 교정, 최고 ~4,000/µL.
[F3]  CSF 단백 투과계수가 커서 1,517 mg/dL 가 나왔다 (실측 100-500).
      → PS_prot 0.30 → 0.09.  Pb=10 에서 ~300, Pb=19 에서 ~520 mg/dL.
[F4]  CSF 알부민 계수는 정상 상태에서도 Q_alb 25 를 만들었다 (정상 <8).
      → PS_alb 0.55 → 0.105.  정상 Q_alb 5, 수막염 48-87 (실측 30-100+).
[F5]  MAP 감쇠상수 k_map = 0.32 이 평균혈압을 30 mmHg (하한 clamp) 로 떨어뜨려
      CPP 가 1.3 mmHg 가 되었다 — 모든 시나리오가 뇌사로 끝났다.
      → 0.05.  패혈증성 저혈압이 MAP 65-70 으로 수렴 (실측과 정합).
[F6]  MMP-9 생성상수가 CSF MMP-9 를 6,207 ng/mL 로 만들었다 (실측 100-1,000).
      → kMMP 900 → 60, 최고 ~410 ng/mL.
[F7]  장벽 투과성 Pb 가 상한 20 에 붙어버려(19.4) 염증의 세기를 더 이상
      표현하지 못했고, 그 결과 (2)의 반코마이신 이야기가 사라졌다.
      → k_pb 0.55 → 0.020, 최고 Pb ~9.  덱사메타손의 직접 장벽복구 효과
      (k_pb_off × (1+1.5·Idex)) 를 추가하니 침투율이 0.29 → 0.16 (-46 %),
      동물실험 보고와 같은 크기가 되었다.
[F8]  세균 포도당 소비계수 q_bact = 6 mg/dL/h per 1e7 CFU/mL 은 폐렴균
      바이오매스 수율(건조중량 0.3 pg/세포, 균질발효 수율 ~25 g/mol)에서
      계산하면 1.8 이어야 한다.  → 2.0.  덕분에 "호중구가 세균보다 포도당을
      10배 더 먹는다"는 결론이 가정이 아니라 계산 결과로 나온다.
[F9]  세포벽 신호의 K_cw = 0.6 CWU/mL 이 너무 작아 항생제 투여 전부터
      Mg(대식세포 활성)가 포화되어 있었고, 그래서 **투여 직후 사이토카인
      폭발이 관찰되지 않았다**(TNF 531 → 756, 1.4배).  이 모델의 요점이
      수치적으로 지워져 있었다.  → K_cw 10.0 으로 올려 기저 CW 3.3 이
      민감구간에 놓이게 했다.  이제 폭발이 3-5배로 보인다.
[F10] 발열항 k_temp = 1.10 이 체온을 40.3 °C 로 수렴시켰다.  → 0.85, 39.6 °C.
[F11] rhs() 안에서 d[0](세프트리악손 중앙구획)과 Isch(허혈지수)가 각각 두 번
      계산되어 앞의 식이 조용히 덮어써지고 있었다 — 사실상 죽은 코드였다.
      → 제거.
[F12] 항생제 T>MIC 추적자는 세프트리악손에서 언제나 100 % 라 아무것도
      구분하지 못했다.  → 살균 역치인 **T > 4×MIC** (Emax 모델의 EC50)
      로 바꿨다.  그러자 내성균(MIC 4 mg/L)에서 세프트리악손이 16 mg/L 를
      한 번도 넘지 못한다는 사실 — 반코마이신을 더하는 임상적 이유 — 이
      출력으로 드러난다.
[F13] 적분 안정성: ICP 방정식의 가장 빠른 고유값은 K_el·ICP/R_out 로
      회복기(R_out 0.167)에 최대가 된다.  ICP 60 · R_out 0.167 에서
      λ = 33/h 이므로 RK4 안정조건 dt < 2.785/λ = 0.084 h.  후기 스텝을
      0.04 h 로 잡고 dt 절반(0.02 h) 검산에서 일치를 확인했다 (아래 [8]).

-------------------------------------------------------------------------------
2차 수치검증에서 추가로 잡아낸 결함 (v2 → v3)
-------------------------------------------------------------------------------
[F14] **혈중 세균에 수용력이 없었다.**  dN_b/dt = (μ_b − clearance)·N_b 가
      순증식이면 N_b 가 336 h 동안 e^16.8 = 2×10^10 CFU/mL 까지 자란다.
      그 값이 k_seed 로 CSF 에 역주입되어 CSF 균수를 N_max 를 넘겨 밀어붙이고,
      세균 포도당 소비항이 폭주해 **CSF 락테이트 292 mmol/L** (실측 6-12) 가
      출력되었다.  → N_b 에 logistic 수용력 (N_b,max = 10^8 CFU/mL) 추가.
      이후 락테이트는 포도당 유입 상한이 자동으로 묶어 ≤9 mmol/L.
[F15] **삼투압 단위가 10배 틀렸다.**  mg/L → mOsm/kg 환산에 만니톨 MW 182,
      글리세롤 MW 92 대신 18.2 / 9.2 를 썼다.  0.5 g/kg 만니톨이 혈장
      삼투압을 113 mOsm/kg 올린다는 (생리적으로 불가능한) 값이 나왔다.
      → 182 / 92 로 교정하니 11.3 mOsm/kg (문헌 10-15) 가 되었고, 그만큼
      효과가 사라지므로 k_osm 을 0.045 → 1.2 mL/h per mOsm 으로 재보정했다
      (0.5 g/kg 이 뇌수분 ~13 mL/h 를 빼내는, 보고된 크기).
[F16] **위험함수가 영구손상을 시간에 곱하고 있었다.**  Haz 에 h_cort·(1−N_cort)
      항이 들어 있어, 회복한 환자도 남은 피질손상 때문에 14일 내내 위험을
      적립했다.  그 결과 24 시나리오 중 22개가 사망확률 99 % 이상으로 나와
      **모든 치료 비교가 무의미해졌다.**  → 위험함수는 급성 생리이탈
      (ICP·CPP·SOFA·균혈증)만 적분하고, 구조적 손상은 종료시점에 1회
      평가하는 항으로 분리했다.
[F17] **자동조절 공식이 자동조절을 해로운 것으로 만들었다.**  CBF =
      CBF0·(AutoR·plateau + (1−AutoR)·CPP/CPP0) 에서 plateau 가 CPP<45 에
      0 이 되므로, 자동조절이 남아 있는(AutoR 큰) 환자가 압력수동
      환자보다 CBF 가 더 낮게 계산되었다 (CPP 38 에서 Isch 0.54 vs 0.07).
      → CBF = CBF0·[f_passive + AutoR·max(0, f_auto − f_passive)] 로 바꿔
      자동조절이 결코 손해가 되지 않도록 했다.
[F18] **손상 상수들이 두 자릿수 컸다.**  k_apo_dg 0.026 /h 는 치상핵을 30 h
      시간상수로 없애 모든 시나리오의 인지 z 를 −3 으로 만들었고, k_hc 0.024
      는 전원을 90 dB 이상 난청으로 만들었다 (실측: 어떤 난청 20-30 %,
      중증 ~10 %).  → k_apo_dg 0.0064, k_hc 0.0040.  청력역치 변환도
      단순 거듭제곱에서 **와우 예비능을 반영한 문턱 곡선**으로 교체했다:
      dB = 120·clamp((L−0.15)/0.55, 0, 1)^0.9 (L = 유모세포 소실분율).
[F19] **멸균시간이 부착 아집단에 의해 결정되고 있었다.**  임상에서 재검
      요추천자 배양이 보는 것은 CSF 유리균이다.  부착·격리 아집단(kill×0.35)
      을 함께 세니 감수성 균주에서 42.6 h 가 나왔다 (실측 4-24 h).
      → 멸균시간은 N_c 기준, 부착 아집단 제거시간은 별도 엔드포인트로 분리
      (완전 치료기간이 필요한 이유가 이 간극으로 표현된다).
[F20] 세프트리악손 CSF 침투 파라미터 PS_cef = 4.0 은 중증 염증(Pb 9)에서
      CSF 농도를 20 mg/L 이상으로 올렸다 (보고 1-12 mg/L).  → 3.0.
[F21] 자연 자가용해율 k_autolysis 0.020 /h 는 항생제 투여 전 세포벽 부하를
      이미 높게 만들어 **투여 직후 폭발이 1.4배로만 보였다.**  → 0.008 /h
      (시간당 0.8 % 자가용해).  이제 CW 1.3 → 17 (13배) 로 보인다.
[F22] 반코마이신 Emax 0.90 /h 는 실험적 폐렴균 수막염의 보고 살균속도
      (−0.3~−0.4 log10 CFU/mL/h, 증식 상쇄 후)와 맞지 않았다.  → 1.20.
[F23] 숙주방어 지수 host_def 가 CSF 탐식(원래 작은 항)만 조절해서
      **면역저하 시나리오가 정상 숙주와 완전히 동일한 출력**을 냈다.
      → 억제 실패를 반영해 순증식도 함께 조절 (μ_scale = 1 + 0.5·(1−host_def)).
[F24] 지속주입 시나리오가 간헐투여와 동일했다 — 감수성 균주에서는 C/MIC
      여유가 200배라 T>4×MIC 가 어느 쪽이든 100 % 이기 때문이다.  이것은
      결함이 아니라 결과이므로, 여유가 없는 **내성균에서의 지속주입**
      시나리오(S25)를 추가해 구분이 생기는 자리를 보이도록 했다.
[F33] **t < 0 에 예정된 투여가 조용히 사라졌다 — 그리고 그것이 이 모델의
      중심 주장을 검증 불가능하게 만들고 있었다.**  "덱사메타손을 첫 항생제보다
      20분 먼저"는 t = −0.33 h 의 주입으로 쓰기 자연스럽지만, 적분은 t = 0 에서
      시작하므로 그 주입 구간은 한 번도 평가되지 않는다.  그래서 −20분 팔의
      실제 첫 투여는 두 번째 용량이 들어가는 **+5.67 h** 였다.  즉 "먼저 준
      팔"이 "+4 h 팔"보다 늦게 받고 있었고, 타이밍이 무관해 보이던 것은
      생리가 아니라 일정표 버그였다.  같은 이유로 "리팜핀 2 h 선행" 시나리오는
      리팜핀을 **전혀** 받지 않았다 — 그 증거는 출력에 있었다: S16 의 최고
      세포벽·용소·TNF 가 S02 와 소수점 두 자리까지 동일했다 (17.4 / 5.58 / 737).
      → 전체 일정을 가장 이른 사건이 t=0 이 되도록 평행이동한다 (물리적으로도
      이것이 맞는 해석이다: 발현 시점에 스테로이드를 먼저 주고 20분 뒤 항생제).
[F25] **피질 신경손실이 ROS 항에 의해 만성적으로 진행했다.**  k_ros_cort 0.010 /h
      × hill(ROS) 0.73 이 100 h 넘게 꺼지지 않아, 최적 치료 환자도 피질의 56 %
      를 잃고 사망확률이 87 % 로 나왔다.  → 0.0008 /h.  실제 폐렴균 수막염의
      피질 손상은 반점상이고 대부분 관류 실패에 따른 것이다.
[F26] 치상핵·유모세포 상수도 여전히 커서 최적 치료에서 인지 z −1.6, 청력
      48 dB 이 나왔다 (실측: 생존자 인지장애 ~30 %, 어떤 난청 20-30 %).
      → k_apo_dg 0.0064 → 0.0020, k_hc 0.0040 → 0.0022.
[F27] 위험함수 계수를 눈으로 맞추는 것을 그만두고, **성분별 적분을 상태변수로
      추가**했다 (I_icp, I_cpp, I_sofa, I_isch, I_bact — 상태 57개 → 63개).
      이제 Haz = h0·T + Σ h_i·I_i 로 사후 분해되므로, 계수를 임상시험
      목표치(폐렴균 사망 34 %/14 %)에 맞춰 산술로 풀 수 있다 (아래 [10]).

단위: 시간 h · 약물 mg/L · 세균 CFU/mL · 사이토카인 pg/mL · PMN cells/µL ·
      포도당 mg/dL · 락테이트 mmol/L · 압력 mmHg · 용적 mL.
"""

import math
import sys
from bisect import bisect_right

# ===========================================================================
# 상태변수 (63 ODE)
# ===========================================================================
STATES = [
    # --- 세프트리악손 2-구획 + CSF (0-2)
    "Cef_c", "Cef_p", "Cef_csf",
    # --- 반코마이신 2-구획 + CSF (3-5)
    "Van_c", "Van_p", "Van_csf",
    # --- 리팜핀 1-구획 + CSF (6-7)
    "Rif_c", "Rif_csf",
    # --- 덱사메타손 + CSF + 전사효과구획 (8-10)
    "Dex_c", "Dex_csf", "Dex_TR",
    # --- 삼투요법: 만니톨, 글리세롤(장·혈장), 뇌내 축적 삼투질 (11-14)
    "Mann_c", "Gly_gut", "Gly_c", "Osm_br",
    # --- 세균과 그 산물 (15-19)
    "Nc", "Nadh", "Nb", "CW", "PLY",
    # --- 염증 (20-29)
    "Mg", "TNF", "IL1", "IL6", "IL10", "CXCL8", "Comp", "PMN", "MMP9", "ROS",
    # --- 장벽·CSF 역학 (30-39)
    "Pb", "Alb_csf", "Prot_csf", "Glc_csf", "Lac_csf", "Vcsf_net", "R_out",
    "ICP", "Vbr", "AutoR",
    # --- 전신 (40-43)
    "MAP", "Temp", "SOFA", "Vol",
    # --- 손상·엔드포인트 (44-49)
    "Ncort", "Ndg", "HC", "Oss", "Sz", "Haz",
    # --- 추적자 (50-56)
    "AUC_cef", "AUC_van", "AUC_lysis", "AUC_TNF", "T_cef", "T_van", "AUC_ICP",
    # --- 위험함수 성분별 적분 (57-62) — 보정을 추측이 아니라 산술로 하기 위해
    "I_icp", "I_cpp", "I_sofa", "I_isch", "I_bact", "I_ncsf",
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)

# ===========================================================================
# 파라미터
# ===========================================================================
P = dict(
    # ---------------- 세프트리악손 PK (2 g IV q12h, 0.5 h 주입) --------------
    # Cmax(총) ~250 mg/L · t1/2 6-9 h · 단백결합 85-95 % 가 포화성
    V1_cef=8.0, V2_cef=8.0, CL_cef=1.0, Q_cef=1.0,
    Bmax_cef=333.0, Kd_cef=25.0,               # 알부민 Langmuir (mg/L)
    PS_cef=3.0, Eff_cef=0.0,                   # CSF 확산 clearance (mL/h) @ Pb=1  [F20]
    MIC_cef=0.03,
    Emax_cef=1.40, EC50r_cef=4.0, h_cef=2.0,   # 1.40/h = 0.61 log10 CFU/mL/h
    Y_cef=1.00,                                # 용해성 살균 → 세포벽 수확량 최대
    # ---------------- 반코마이신 PK (15 mg/kg q6h, 1 h 주입) ---------------
    V1_van=20.0, V2_van=30.0, CL_van=4.0, Q_van=8.0, fu_van=0.50,
    PS_van=1.4, Eff_van=10.0,                  # 능동유출 있음
    MIC_van=1.0,
    Emax_van=1.20, EC50r_van=4.0, h_van=1.5,   # [F22]
    Y_van=0.90,
    # ---------------- 리팜핀 PK (600 mg IV q12h) ---------------------------
    V_rif=50.0, CL_rif=12.0, fu_rif=0.20,
    PS_rif=30.0, a_rif=0.10, Eff_rif=0.0,      # 친유성 → Pb 민감도 낮음
    MIC_rif=0.03,
    Emax_rif=0.70, EC50r_rif=2.0, h_rif=1.0,
    Y_rif=0.15,                                # 비용해성 → 수확량 작다
    # ---------------- 덱사메타손 PK/PD (0.15 mg/kg q6h) --------------------
    V_dex=70.0, CL_dex=16.0, fu_dex=0.32, PS_dex=32.0,
    ktr_on=0.46, ktr_off=0.040,                # 비대칭 (켜짐 t1/2 1.5 h / 꺼짐 17 h)
    Imax_dex=0.80, IC50_TR=0.010,
    # ---------------- 삼투요법 ---------------------------------------------
    V_mann=17.0, CL_mann=6.0,
    ka_gly=1.2, V_gly=42.0, CL_gly=8.0,
    MW_mann=182.0, MW_gly=92.0,                # [F15] 10배 단위오류 교정
    k_osm=1.20, k_osm_leak=0.010, k_osm_br_out=0.030,
    # ---------------- CSF 물리 -------------------------------------------
    Vcsf=150.0, Qf0=21.0, P_ss=6.0,            # mL · mL/h(0.35 mL/min) · mmHg
    R_out0=0.167, R_out_max=1.20,              # mmHg/(mL/h)
    k_ro_on=0.09, k_ro_off=0.020,
    K_el=0.092,                                # 탄성계수 1/mL  (PVI 25 mL)
    # ---------------- 세균 -------------------------------------------------
    mu_max=0.85, K_glc=10.0, Nmax=1.0e9,
    k_adh=0.030, k_des=0.010, prot_adh=0.35, mu_adh=0.30,
    k_shed=2.0e-4, k_seed=0.010,
    mu_b=0.45, k_clr_b=0.40, kill_b_boost=1.0, Nb_max=1.0e8,   # [F14]
    kphag_max=0.15, K_pmn_ph=2500.0, K_comp=1.0,
    k_autolysis=0.008, Y_auto=1.0, Y_phag=0.30,                # [F21]
    # ---------------- 세포벽·용소 -----------------------------------------
    kCW_cl=0.060,                              # t1/2 ~12 h (세포벽은 오래 남는다)
    yPLY=0.50, kPLY=0.35, ply_secr=0.006,
    # ---------------- 자연면역 --------------------------------------------
    kmg_on=1.10, kmg_off=0.070,
    K_cw=10.0, K_nl=1.20, K_ply=1.50,          # [F9]
    kTNF=1400.0, kel_TNF=0.70, K_tnf=250.0,
    kIL1=320.0, kel_IL1=0.35, K_il1=120.0,
    kIL6=9.0e3, kel_IL6=0.25, K_il6=4000.0,
    kIL10=260.0, kel_IL10=0.20, K_il10=600.0,
    kC8=5.0e3, kel_C8=0.30, K_c8=1500.0,
    kcomp=1.20, kel_comp=0.25,
    k_influx=500.0, k_egress=0.055, k_apop=0.045,   # [F2]
    kMMP=60.0, kel_MMP=0.087, K_mmp=350.0,          # [F6]
    kROS=1.30, kel_ROS=0.90, K_ros=0.60,
    # ---------------- 장벽 -------------------------------------------------
    Pb_max=20.0, k_pb=0.020, k_pb_off=0.030, dex_pb=1.5,   # [F7]
    PS_alb=0.105, Alb_ser=42000.0,                          # [F4]
    PS_prot=0.090, Prot_ser=70000.0,                        # [F3]
    # ---------------- 포도당·락테이트 ------------------------------------
    Tmax_glc=74.0, Km_glut=90.0, Glc_pl=100.0,
    q_pmn=0.0050, q_bact=2.0, q_brain=1.0,                  # [F8]
    Km_use=5.0,                                             # [F1]
    inflam_glut=0.30,
    k_lac=0.500, Lac_base=1.6,
    # ---------------- 뇌부종·관류 ----------------------------------------
    k_vas=1.05, Vbr_max=60.0, k_cyt=2.2, k_vbr_res=0.020,
    k_ar_loss=0.16, k_ar_rec=0.012,
    CBF0=50.0, CPP0=75.0, CBF_crit=0.55,
    # ---------------- 전신 -------------------------------------------------
    MAP0=88.0, k_map=0.050, k_map_rec=0.10, k_cush=0.30,    # [F5]
    k_temp=0.85, k_temp_off=0.35, Temp0=37.0,               # [F10]
    k_sofa=0.55, k_sofa_rec=0.10,
    # ---------------- 손상 -------------------------------------------------
    k_isch=0.055, k_ros_cort=0.0008,                           # [F25]
    k_apo_dg=0.0020, k_ply_dg=0.5, k_ros_dg=0.4,               # [F18][F26]
    k_hc=0.0022, w_hc_ply=0.5, w_hc_ros=0.7, w_hc_pmn=0.4,     # [F18][F26]
    K_ply2=1.0, k_oss=0.0020,
    # [F32] 청력 문턱곡선.  예비능 0.15/구간 0.55 로는 모든 환자가 20-35 dB 의
    #       중간 난청에 뭉쳐, 실측 분포(어떤 난청 20-30 % · 중증 ~10 %)의
    #       꼬리가 재현되지 않았다.  예비능을 늘리고 구간을 좁혀 곡선을 세웠다.
    hc_reserve=0.22, hc_span=0.40, hc_exp=0.90,
    k_sz=0.10, k_sz_off=0.12,
    # ---------------- 사망 위험함수 --------------------------------------
    # [F16] 급성 생리이탈만 적분한다.  구조적 손상은 종료시점 1회 평가.
    # [F30] 계수는 [10] 절의 성분별 적분과 임상시험 목표치(폐렴균 사망 34 %/14 %)
    # 로부터 산술적으로 결정했다 — 눈으로 맞춘 값이 아니다.
    #   1차 시도: 모든 계수 s배 → 코호트 무-DEX 1.256·s = 0.416 → s = 0.331
    #                              코호트 DEX   0.538·s = 0.151 → s = 0.281
    #   모델의 위험비 2.33 이 목표 2.76 보다 약간 작아 두 s 가 갈리므로,
    #   타협값 s = 0.30 을 모든 계수에 적용했다 (예측 31 % / 15 %).
    h0=1.2e-5, h_icp=6.0e-4, h_cpp=6.0e-4, h_sofa=6.0e-5, h_bact=4.5e-5,
    h_cort_final=0.48, h_isch=6.0e-4,
    # [F31] CSF 내 지속감염 자체의 위험.  이 항이 없으면 무치료 자연경과가
    #       14일 사망 70 % 로 나온다 (실제로는 거의 100 %).  치료군에서는
    #       N_c 가 24 h 안에 사라지므로 기여가 무시할 만하다.
    h_ncsf=1.5e-4,
)

# ===========================================================================
# 보조 함수
# ===========================================================================
def free_saturable(Ctot, Bmax, Kd):
    """포화성 알부민 결합에서 유리농도: Cf + Bmax·Cf/(Kd+Cf) = Ctot 의 양근."""
    if Ctot <= 0.0:
        return 0.0
    b = Kd + Bmax - Ctot
    return (-b + math.sqrt(b * b + 4.0 * Kd * Ctot)) / 2.0


def hill(x, K, n=1.0):
    if x <= 0.0:
        return 0.0
    xn = x if n == 1.0 else x ** n
    Kn = K if n == 1.0 else K ** n
    return xn / (Kn + xn)


def emax_kill(C, MIC, Emax, EC50r, h):
    """C/MIC 비에 대한 Emax 살균속도 (/h).  '10× MBC' 경험칙과 정합."""
    if C <= 0.0 or MIC <= 0.0:
        return 0.0
    r = C / MIC
    rh = r ** h
    return Emax * rh / (EC50r ** h + rh)


def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


# ===========================================================================
# 초기조건: 발현(입원) 시점.  증상 발생 12-36 h 후이므로 이미 염증이 진행되어 있다
# ===========================================================================
def init_state(N0=1.0e7, Nb0=1.0e3):
    y = [0.0] * NS
    y[IX["Nc"]] = N0
    y[IX["Nadh"]] = N0 * 0.05
    y[IX["Nb"]] = Nb0
    y[IX["Mg"]] = 0.55
    y[IX["TNF"]] = 420.0
    y[IX["IL1"]] = 180.0
    y[IX["IL6"]] = 2.5e4
    y[IX["IL10"]] = 260.0
    y[IX["CXCL8"]] = 4.0e3
    y[IX["Comp"]] = 0.9
    y[IX["PMN"]] = 1200.0
    y[IX["MMP9"]] = 160.0
    y[IX["ROS"]] = 0.55
    y[IX["Pb"]] = 6.0
    y[IX["Alb_csf"]] = 1400.0        # Q_alb ~33
    y[IX["Prot_csf"]] = 180.0        # mg/dL
    y[IX["Glc_csf"]] = 34.0
    y[IX["Lac_csf"]] = 5.2
    y[IX["R_out"]] = 0.48
    y[IX["ICP"]] = 18.0
    y[IX["Vbr"]] = 4.0
    y[IX["AutoR"]] = 0.75
    y[IX["MAP"]] = 88.0
    y[IX["Temp"]] = 39.1
    y[IX["SOFA"]] = 2.0
    y[IX["Vol"]] = 1.0
    y[IX["Ncort"]] = 1.0
    y[IX["Ndg"]] = 1.0
    y[IX["HC"]] = 1.0
    y[IX["Oss"]] = 0.0
    return y


# ===========================================================================
# 투여 계획 (구간 정속 IV 주입 + 경구 볼루스)
# ===========================================================================
class Regimen:
    def __init__(self):
        self.by_ix = {}
        self.boluses = []
        self._sorted = None

    def _add(self, t0, t1, rate, ix):
        self.by_ix.setdefault(ix, []).append((t0, t1, rate))
        self._sorted = None

    def iv(self, state, dose_mg, t_first, interval, n_dose, dur=0.5):
        for k in range(n_dose):
            t0 = t_first + k * interval
            self._add(t0, t0 + dur, dose_mg / dur, IX[state])

    def infuse(self, state, rate, t0, t1):
        self._add(t0, t1, rate, IX[state])

    def po(self, state, dose_mg, t_first, interval, n_dose):
        for k in range(n_dose):
            self.boluses.append((t_first + k * interval, dose_mg, IX[state]))

    def shift_to_zero(self):
        """[F33] t<0 에 예정된 투여를 살린다.

        "덱사메타손을 첫 항생제보다 20분 먼저"는 t = −0.33 h 의 주입으로
        쓰기 쉽지만, 적분이 t = 0 에서 시작하므로 그 주입은 **조용히
        사라진다**.  그래서 −20분 팔은 실제로는 두 번째 용량이 들어가는
        +5.67 h 가 첫 투여가 되어 있었고, 타이밍 비교 전체가 무의미했다.
        (같은 이유로 리팜핀 2 h 선행 시나리오는 리팜핀을 전혀 받지 않았다 —
        그래서 S16 의 세포벽·용소·TNF 최고치가 S02 와 소수점까지 같았다.)

        물리적으로 맞는 해석은 "발현(t=0)에 스테로이드를 먼저 주고 20분 뒤
        항생제"이므로, 전체 일정을 가장 이른 사건이 t=0 이 되도록 평행이동한다.
        """
        times = [t0 for lst in self.by_ix.values() for (t0, _, _) in lst]
        times += [tb for (tb, _, _) in self.boluses]
        if not times:
            return 0.0
        tmin = min(times)
        if tmin >= 0.0:
            return 0.0
        sh = -tmin
        self.by_ix = {ix: [(t0 + sh, t1 + sh, r) for (t0, t1, r) in lst]
                      for ix, lst in self.by_ix.items()}
        self.boluses = [(tb + sh, a, ix) for (tb, a, ix) in self.boluses]
        self._sorted = None
        return sh

    def _prep(self):
        self._sorted = {}
        for ix, lst in self.by_ix.items():
            lst = sorted(lst)
            self._sorted[ix] = (lst, [e[0] for e in lst])

    def rate(self, t, idx):
        if self._sorted is None:
            self._prep()
        ent = self._sorted.get(idx)
        if not ent:
            return 0.0
        lst, starts = ent
        j = bisect_right(starts, t) - 1
        r = 0.0
        while j >= 0:
            t0, t1, rt = lst[j]
            if t < t1:
                r += rt
            if t0 < t - 48.0:          # 48 h 이전 주입은 끝났음이 보장됨
                break
            j -= 1
        return r

    def due(self, t, dt):
        return [(a, ix) for (tb, a, ix) in self.boluses if t <= tb < t + dt]


# ===========================================================================
# 우변 (RHS)
# ===========================================================================
def rhs(t, y, p, reg, cfg):
    d = [0.0] * NS

    Cef_c, Cef_p, Cef_csf = y[0], y[1], y[2]
    Van_c, Van_p, Van_csf = y[3], y[4], y[5]
    Rif_c, Rif_csf = y[6], y[7]
    Dex_c, Dex_csf, Dex_TR = y[8], y[9], y[10]
    Mann_c, Gly_gut, Gly_c, Osm_br = y[11], y[12], y[13], y[14]
    Nc, Nadh, Nb, CW, PLY = y[15], y[16], y[17], y[18], y[19]
    Mg, TNF, IL1, IL6, IL10, CXCL8, Comp, PMN, MMP9, ROS = y[20:30]
    Pb, Alb, Prot, Glc, Lac, _Vnet, R_out, ICP, Vbr, AutoR = y[30:40]
    MAP, Temp, SOFA, Vol = y[40:44]
    Ncort, Ndg, HC, Oss, Sz, _Haz = y[44:50]

    # ---------- 혈장 유리농도 ----------
    Cef_pl = Cef_c / p["V1_cef"]
    Cef_f = free_saturable(Cef_pl, p["Bmax_cef"], p["Kd_cef"])
    Van_pl = Van_c / p["V1_van"]
    Van_f = p["fu_van"] * Van_pl
    Rif_pl = Rif_c / p["V_rif"]
    Rif_f = p["fu_rif"] * Rif_pl
    Dex_pl = Dex_c / p["V_dex"]
    Dex_f = p["fu_dex"] * Dex_pl

    Pb_hyd = Pb                                  # 친수성 약물
    Pb_lip = 1.0 + p["a_rif"] * (Pb - 1.0)       # 친유성 약물
    Qbulk = p["Qf0"]                             # CSF bulk flow = sink (mL/h)

    # ---------- 약물 PK ----------
    d[0] = reg.rate(t, 0) - p["CL_cef"] * Cef_pl - p["Q_cef"] * (Cef_pl - Cef_p / p["V2_cef"])
    d[1] = p["Q_cef"] * (Cef_pl - Cef_p / p["V2_cef"])
    d[2] = (p["PS_cef"] * Pb_hyd * (Cef_f - Cef_csf) - (Qbulk + p["Eff_cef"]) * Cef_csf) / p["Vcsf"]

    d[3] = reg.rate(t, 3) - p["CL_van"] * Van_pl - p["Q_van"] * (Van_pl - Van_p / p["V2_van"])
    d[4] = p["Q_van"] * (Van_pl - Van_p / p["V2_van"])
    d[5] = (p["PS_van"] * Pb_hyd * (Van_f - Van_csf) - (Qbulk + p["Eff_van"]) * Van_csf) / p["Vcsf"]

    d[6] = reg.rate(t, 6) - p["CL_rif"] * Rif_pl
    d[7] = (p["PS_rif"] * Pb_lip * (Rif_f - Rif_csf) - (Qbulk + p["Eff_rif"]) * Rif_csf) / p["Vcsf"]

    d[8] = reg.rate(t, 8) - p["CL_dex"] * Dex_pl
    d[9] = (p["PS_dex"] * Pb_lip * (Dex_f - Dex_csf) - Qbulk * Dex_csf) / p["Vcsf"]
    ktr = p["ktr_on"] if Dex_csf > Dex_TR else p["ktr_off"]
    d[10] = ktr * (Dex_csf - Dex_TR)
    Idex = p["Imax_dex"] * hill(Dex_TR, p["IC50_TR"])

    # ---------- 삼투요법 ----------
    d[11] = reg.rate(t, 11) - p["CL_mann"] * (Mann_c / p["V_mann"])
    gly_abs = p["ka_gly"] * Gly_gut
    d[12] = -gly_abs
    d[13] = gly_abs - p["CL_gly"] * (Gly_c / p["V_gly"])
    # [F15] mg/L ÷ MW(g/mol) = mmol/L = mOsm/kg 등가 (10배 단위오류 교정)
    Osm_pl = (Mann_c / p["V_mann"]) / p["MW_mann"] + (Gly_c / p["V_gly"]) / p["MW_gly"]
    Osm_grad = max(0.0, Osm_pl - Osm_br)
    d[14] = p["k_osm_leak"] * (1.0 + 0.35 * (Pb - 1.0)) * Osm_pl - p["k_osm_br_out"] * Osm_br

    # ---------- 세균 ----------
    kill_cef = emax_kill(Cef_csf, cfg["MIC_cef"], p["Emax_cef"], p["EC50r_cef"], p["h_cef"])
    kill_van = emax_kill(Van_csf, cfg["MIC_van"], p["Emax_van"], p["EC50r_van"], p["h_van"])
    kill_rif = emax_kill(Rif_csf, cfg["MIC_rif"], p["Emax_rif"], p["EC50r_rif"], p["h_rif"])
    kill_tot = kill_cef + kill_van + kill_rif

    f_glc = hill(Glc, p["Km_use"])                     # [F1] 기질 의존성
    mu = p["mu_max"] * hill(Glc, p["K_glc"]) * (1.0 - Nc / p["Nmax"]) * cfg["mu_scale"]
    k_phag = (p["kphag_max"] * hill(PMN, p["K_pmn_ph"])
              * (0.2 + 0.8 * hill(Comp, p["K_comp"])) * cfg["host_def"])

    d[15] = ((mu - kill_tot - k_phag) * Nc + p["k_des"] * Nadh - p["k_adh"] * Nc
             + p["k_seed"] * Nb - p["k_shed"] * Nc)
    d[16] = (p["k_adh"] * Nc - p["k_des"] * Nadh
             + (p["mu_adh"] * hill(Glc, p["K_glc"]) - kill_tot * p["prot_adh"]) * Nadh)
    # [F14] 혈중 세균에도 수용력이 있어야 한다 (없으면 N_b 가 발산해 CSF 를 역오염)
    d[17] = ((p["mu_b"] * (1.0 - Nb / p["Nb_max"]) - p["k_clr_b"] * cfg["host_def"]
              - p["kill_b_boost"] * kill_tot) * Nb + p["k_shed"] * Nc)

    # ---------- 용해 flux = 죽이는 속도 × 균 수 × 수확량 ----------
    Y_adh = p["Y_cef"] if (kill_cef + kill_van) > kill_rif else p["Y_rif"]
    lysis_flux = ((kill_cef * p["Y_cef"] + kill_van * p["Y_van"] + kill_rif * p["Y_rif"]) * Nc
                  + kill_tot * p["prot_adh"] * Y_adh * Nadh
                  + p["k_autolysis"] * p["Y_auto"] * Nc
                  + k_phag * p["Y_phag"] * Nc)
    d[18] = lysis_flux / 1.0e6 - p["kCW_cl"] * CW
    d[19] = p["yPLY"] * (lysis_flux + p["ply_secr"] * Nc) / 1.0e6 - p["kPLY"] * PLY

    # ---------- 관류 (ROS·손상에 필요) ----------
    # [F17] 자동조절은 결코 손해가 될 수 없다: 압력수동 흐름을 하한으로 둔다
    CPP = MAP - ICP
    f_passive = clamp(CPP / p["CPP0"], 0.0, 1.6)
    f_auto = clamp((CPP - 25.0) / 25.0, 0.0, 1.0)
    CBF = p["CBF0"] * (f_passive + AutoR * max(0.0, f_auto - f_passive))
    Isch = clamp(1.0 - CBF / (p["CBF_crit"] * p["CBF0"]), 0.0, 1.0)

    # ---------- 자연면역 ----------
    S_pamp = (hill(CW, p["K_cw"]) + 0.6 * hill(Nc / 1.0e6, p["K_nl"])
              + 0.5 * hill(PLY, p["K_ply"]))
    d[20] = p["kmg_on"] * S_pamp * (1.0 - Mg) - p["kmg_off"] * Mg
    f10 = 1.0 / (1.0 + IL10 / p["K_il10"])
    d[21] = p["kTNF"] * Mg * (1.0 - Idex) * f10 - p["kel_TNF"] * TNF
    d[22] = (p["kIL1"] * Mg * (0.6 + 0.4 * hill(PLY, p["K_ply"])) * (1.0 - Idex) * f10
             - p["kel_IL1"] * IL1)
    d[23] = (p["kIL6"] * (0.5 * Mg + hill(TNF, p["K_tnf"]) + hill(IL1, p["K_il1"]))
             * (1.0 - 0.6 * Idex) - p["kel_IL6"] * IL6)
    d[24] = (p["kIL10"] * (Mg + hill(TNF, p["K_tnf"])) * (1.0 + 0.5 * Idex)
             - p["kel_IL10"] * IL10)
    d[25] = (p["kC8"] * (hill(TNF, p["K_tnf"]) + hill(IL1, p["K_il1"]) + 0.3 * Mg)
             * (1.0 - 0.7 * Idex) - p["kel_C8"] * CXCL8)
    d[26] = p["kcomp"] * hill(Nc / 1.0e6, 2.0) - p["kel_comp"] * Comp
    adh_mol = hill(TNF, p["K_tnf"]) * (1.0 - 0.6 * Idex)
    d[27] = (p["k_influx"] * hill(CXCL8, p["K_c8"]) * adh_mol
             * (1.0 + 0.5 * (Pb - 1.0) / 19.0) - (p["k_egress"] + p["k_apop"]) * PMN)
    d[28] = p["kMMP"] * hill(PMN, 2500.0) * (1.0 - 0.5 * Idex) - p["kel_MMP"] * MMP9
    d[29] = p["kROS"] * (hill(PMN, 2500.0) + 0.5 * Mg + 0.3 * Isch) - p["kel_ROS"] * ROS

    # ---------- 장벽 ----------
    d[30] = (p["k_pb"] * (hill(MMP9, p["K_mmp"]) + 0.5 * hill(TNF, p["K_tnf"])
                          + 0.3 * hill(PMN, 2500.0)) * (p["Pb_max"] - Pb)
             - p["k_pb_off"] * (1.0 + p["dex_pb"] * Idex) * (Pb - 1.0))
    d[31] = (p["PS_alb"] * Pb * (p["Alb_ser"] - Alb) - Qbulk * Alb) / p["Vcsf"]
    d[32] = (p["PS_prot"] * Pb * (p["Prot_ser"] / 10.0 - Prot) - Qbulk * Prot) / p["Vcsf"]

    # ---------- 포도당·락테이트 ----------
    glut = 1.0 - p["inflam_glut"] * hill(Pb - 1.0, 8.0)
    influx = p["Tmax_glc"] * glut * (p["Glc_pl"] / (p["Km_glut"] + p["Glc_pl"])
                                     - Glc / (p["Km_glut"] + Glc))
    use_pmn = p["q_pmn"] * PMN * f_glc
    use_bact = p["q_bact"] * (Nc / 1.0e7) * f_glc
    use_brain = p["q_brain"] * (1.0 + 0.3 * (Temp - 37.0)) * f_glc
    d[33] = influx - use_pmn - use_bact - use_brain - (Qbulk / p["Vcsf"]) * Glc
    d[34] = (2.0 * 0.0556 * (use_pmn + use_bact + use_brain * (0.3 + 0.7 * Isch))
             - p["k_lac"] * (Lac - p["Lac_base"]))

    # ---------- CSF 역학 ----------
    Qf = p["Qf0"] * (1.0 - 0.3 * clamp((ICP - 20.0) / 40.0, 0.0, 1.0))
    d[36] = (p["k_ro_on"] * (hill(Prot, 200.0) + hill(PMN, 2000.0))
             * (p["R_out_max"] - R_out) - p["k_ro_off"] * (R_out - p["R_out0"]))
    absorb = max(0.0, (ICP - p["P_ss"]) / max(R_out, 1e-3))
    drain = cfg["csf_drain"](t) if cfg["csf_drain"] else 0.0
    osm_shrink = p["k_osm"] * Osm_grad / (1.0 + 0.3 * (Pb - 1.0))
    dVbr = (p["k_vas"] * hill(Pb - 1.0, 6.0) * (1.0 - Vbr / p["Vbr_max"])
            + p["k_cyt"] * Isch - osm_shrink - p["k_vbr_res"] * Vbr)
    d[38] = dVbr
    d[35] = Qf - absorb - drain                      # 추적자: CSF 순 축적속도
    d[37] = p["K_el"] * ICP * (Qf + dVbr - absorb - drain)
    d[39] = -p["k_ar_loss"] * hill(TNF, p["K_tnf"]) * AutoR + p["k_ar_rec"] * (1.0 - AutoR)

    # ---------- 전신 ----------
    d[40] = (-p["k_map"] * (hill(TNF, p["K_tnf"]) + 0.4 * hill(IL6, p["K_il6"])) * MAP
             + p["k_map_rec"] * (p["MAP0"] * Vol - MAP)
             + p["k_cush"] * max(0.0, ICP - 30.0) + cfg["vasopressor"])
    d[41] = (p["k_temp"] * (hill(IL1, p["K_il1"]) + 0.5 * hill(IL6, p["K_il6"]))
             * (1.0 - cfg["antipyretic"]) - p["k_temp_off"] * (Temp - p["Temp0"]))
    d[42] = (p["k_sofa"] * (hill(math.log10(1.0 + max(Nb, 0.0)), 2.0)
                            + 0.5 * hill(TNF, p["K_tnf"])) - p["k_sofa_rec"] * SOFA)
    d[43] = 0.0

    # ---------- 손상 ----------
    d[44] = -(p["k_isch"] * Isch * Isch + p["k_ros_cort"] * hill(ROS, p["K_ros"])) * Ncort
    d[45] = -(p["k_apo_dg"] * (hill(CW, p["K_cw"]) + p["k_ply_dg"] * hill(PLY, p["K_ply2"])
                               + p["k_ros_dg"] * hill(ROS, p["K_ros"]))) * Ndg
    d[46] = -(p["k_hc"] * (p["w_hc_ply"] * hill(PLY, p["K_ply2"])
                           + p["w_hc_ros"] * hill(ROS, p["K_ros"])
                           + p["w_hc_pmn"] * hill(PMN, 2500.0))) * HC
    d[47] = p["k_oss"] * (1.0 - HC) * hill(PMN, 2000.0) * (1.0 - Oss)
    d[48] = (p["k_sz"] * ((1.0 - Ncort) + hill(ROS, p["K_ros"])
                          + 0.4 * max(0.0, 1.0 - Glc / 25.0))
             * (1.0 - cfg["anticonvulsant"]) - p["k_sz_off"] * Sz)
    # [F16] 급성 생리이탈만 적분 — 영구손상(N_cort)은 종료시점에 1회 평가한다
    # [F28] ICP·CPP 항은 2차식이다.  선형이면 "CPP 45 를 100 h" 와 "CPP 5 를 10 h"
    #       가 같은 위험이 되어, 경증의 지속적 이탈이 위험함수를 지배해버린다.
    z_icp = max(0.0, ICP - 25.0) / 10.0
    z_cpp = max(0.0, 50.0 - CPP) / 10.0
    d[49] = (p["h0"] + p["h_icp"] * z_icp * z_icp
             + p["h_cpp"] * z_cpp * z_cpp
             + p["h_sofa"] * SOFA + p["h_isch"] * Isch * Isch
             + p["h_bact"] * math.log10(1.0 + max(Nb, 0.0))
             + p["h_ncsf"] * math.log10(1.0 + max(Nc, 0.0)))

    # ---------- 추적자 ----------
    d[50] = Cef_csf
    d[51] = Van_csf
    d[52] = lysis_flux / 1.0e6
    d[53] = TNF
    d[54] = 1.0 if Cef_csf > 4.0 * cfg["MIC_cef"] else 0.0      # [F12] T>4×MIC
    d[55] = 1.0 if Van_csf > 4.0 * cfg["MIC_van"] else 0.0
    d[56] = ICP
    # 위험함수 성분별 적분 (Haz = h0·T + Σ h_i · I_i 로 사후 분해가 가능해진다)
    d[57] = z_icp * z_icp
    d[58] = z_cpp * z_cpp
    d[59] = SOFA
    d[60] = Isch * Isch
    d[61] = math.log10(1.0 + max(Nb, 0.0))
    d[62] = math.log10(1.0 + max(Nc, 0.0))
    return d


NONNEG = [IX[s] for s in (
    "Cef_c", "Cef_p", "Cef_csf", "Van_c", "Van_p", "Van_csf", "Rif_c", "Rif_csf",
    "Dex_c", "Dex_csf", "Dex_TR", "Mann_c", "Gly_gut", "Gly_c", "Osm_br", "CW",
    "PLY", "TNF", "IL1", "IL6", "IL10", "CXCL8", "Comp", "PMN", "MMP9", "ROS",
    "Alb_csf", "Prot_csf", "Glc_csf", "Lac_csf", "Vbr", "SOFA", "Sz", "Haz")]


# ===========================================================================
# 적분기 (RK4 · 단계 스케줄)
# ===========================================================================
def simulate(cfg, tmax=336.0, record=None, dt_scale=1.0):
    p = P
    reg = cfg["reg"]
    y = init_state(cfg["N0"], cfg["Nb0"])
    t = 0.0
    out = []
    rec = sorted(record) if record else []
    ri = 0
    sterile_t = None      # [F19] N_c 기준 (재검 요추천자 배양이 보는 것)
    adh_t = None          #       부착 아집단 제거 시점 (완전 치료기간의 근거)
    peak = dict(ICP=y[IX["ICP"]], TNF=y[IX["TNF"]], PMN=y[IX["PMN"]], Nc=y[IX["Nc"]],
                PLY=0.0, Pb=y[IX["Pb"]], Prot=y[IX["Prot_csf"]], Lac=y[IX["Lac_csf"]],
                Temp=y[IX["Temp"]], Qalb=0.0, lysis=0.0, Vbr=y[IX["Vbr"]])
    trough = dict(CPP=y[IX["MAP"]] - y[IX["ICP"]], Glc=y[IX["Glc_csf"]],
                  MAP=y[IX["MAP"]], CBF=50.0)

    while t < tmax - 1e-9:
        dt = (0.005 if t < 4.0 else (0.02 if t < 48.0 else 0.04)) * dt_scale
        if rec and ri < len(rec) and t + dt > rec[ri]:
            dt = max(rec[ri] - t, 1e-6)

        for (amt, ix) in reg.due(t, dt):
            y[ix] += amt

        k1 = rhs(t, y, p, reg, cfg)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NS)]
        k2 = rhs(t + 0.5 * dt, y2, p, reg, cfg)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NS)]
        k3 = rhs(t + 0.5 * dt, y3, p, reg, cfg)
        y4 = [y[i] + dt * k3[i] for i in range(NS)]
        k4 = rhs(t + dt, y4, p, reg, cfg)
        for i in range(NS):
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

        for i in NONNEG:
            if y[i] < 0.0:
                y[i] = 0.0
        y[IX["Mg"]] = clamp(y[IX["Mg"]], 0.0, 1.0)
        y[IX["Pb"]] = clamp(y[IX["Pb"]], 1.0, P["Pb_max"])
        y[IX["AutoR"]] = clamp(y[IX["AutoR"]], 0.0, 1.0)
        y[IX["ICP"]] = clamp(y[IX["ICP"]], 2.0, 120.0)
        y[IX["MAP"]] = clamp(y[IX["MAP"]], 30.0, 160.0)
        for s in ("Ncort", "Ndg", "HC", "Oss"):
            y[IX[s]] = clamp(y[IX[s]], 0.0, 1.0)
        for s in ("Nc", "Nadh", "Nb"):
            if y[IX[s]] < 1.0e-2:          # 검출한계 이하 → 소멸
                y[IX[s]] = 0.0
        t += dt

        if sterile_t is None and y[IX["Nc"]] < 10.0:
            sterile_t = t
        if adh_t is None and y[IX["Nc"]] + y[IX["Nadh"]] < 10.0:
            adh_t = t
        cpp = y[IX["MAP"]] - y[IX["ICP"]]
        peak["ICP"] = max(peak["ICP"], y[IX["ICP"]])
        peak["TNF"] = max(peak["TNF"], y[IX["TNF"]])
        peak["PMN"] = max(peak["PMN"], y[IX["PMN"]])
        peak["Nc"] = max(peak["Nc"], y[IX["Nc"]])
        peak["PLY"] = max(peak["PLY"], y[IX["PLY"]])
        peak["Pb"] = max(peak["Pb"], y[IX["Pb"]])
        peak["Prot"] = max(peak["Prot"], y[IX["Prot_csf"]])
        peak["Lac"] = max(peak["Lac"], y[IX["Lac_csf"]])
        peak["Temp"] = max(peak["Temp"], y[IX["Temp"]])
        peak["Vbr"] = max(peak["Vbr"], y[IX["Vbr"]])
        peak["Qalb"] = max(peak["Qalb"], 1000.0 * y[IX["Alb_csf"]] / P["Alb_ser"])
        trough["CPP"] = min(trough["CPP"], cpp)
        trough["Glc"] = min(trough["Glc"], y[IX["Glc_csf"]])
        trough["MAP"] = min(trough["MAP"], y[IX["MAP"]])

        if rec and ri < len(rec) and abs(t - rec[ri]) < 1e-6:
            out.append((t, list(y)))
            ri += 1

    return dict(y=y, series=out, peak=peak, trough=trough, sterile_t=sterile_t,
                endpoints=endpoints(y, peak, trough, sterile_t, adh_t))


def hearing_dB(HC, Oss):
    """[F18] 와우 예비능(reserve)을 반영한 문턱 곡선.
    유모세포 소실 15 % 까지는 역치이동이 거의 없고, 그 뒤 급격히 오른다."""
    L = 1.0 - HC
    x = clamp((L - P["hc_reserve"]) / P["hc_span"], 0.0, 1.0)
    return min(120.0, 120.0 * (x ** P["hc_exp"]) + 15.0 * Oss)


def endpoints(y, peak, trough, sterile_t, adh_t=None):
    HC, Oss = y[IX["HC"]], y[IX["Oss"]]
    dB = hearing_dB(HC, Oss)
    # [F16] 급성 위험 적분 + 구조적 손상의 1회 평가
    haz_tot = y[IX["Haz"]] + P["h_cort_final"] * (1.0 - y[IX["Ncort"]])
    pdeath = 1.0 - math.exp(-haz_tot)
    cogz = -(2.5 * (1.0 - y[IX["Ndg"]]) + 1.5 * (1.0 - y[IX["Ncort"]]))
    focal = 1.0 - math.exp(-3.0 * (1.0 - y[IX["Ncort"]]))
    seq = min(1.0, (1.0 if dB > 60.0 else 0.0) * 0.5
              + max(0.0, -cogz - 1.0) * 0.35 + focal * 0.4)
    return dict(
        death=pdeath, hear_dB=dB, hear_any=dB > 25.0, hear_severe=dB > 60.0,
        cog_z=cogz, focal=focal, sterile_t=sterile_t, adh_t=adh_t,
        haz_acute=y[IX["Haz"]], haz_struct=P["h_cort_final"] * (1.0 - y[IX["Ncort"]]),
        peak_ICP=peak["ICP"], min_CPP=trough["CPP"], min_MAP=trough["MAP"],
        peak_Vbr=peak["Vbr"], peak_Qalb=peak["Qalb"], peak_Temp=peak["Temp"],
        I_icp=y[IX["I_icp"]], I_cpp=y[IX["I_cpp"]], I_sofa=y[IX["I_sofa"]],
        I_isch=y[IX["I_isch"]], I_bact=y[IX["I_bact"]], I_ncsf=y[IX["I_ncsf"]],
        Ncort_f=y[IX["Ncort"]],
        AUC_lysis=y[IX["AUC_lysis"]], AUC_TNF=y[IX["AUC_TNF"]],
        AUC_cef=y[IX["AUC_cef"]], AUC_van=y[IX["AUC_van"]],
        T_cef=y[IX["T_cef"]], T_van=y[IX["T_van"]],
        peak_TNF=peak["TNF"], peak_PMN=peak["PMN"], peak_Pb=peak["Pb"],
        peak_PLY=peak["PLY"], min_Glc=trough["Glc"], peak_Lac=peak["Lac"],
        peak_Prot=peak["Prot"], Ndg=y[IX["Ndg"]], Ncort=y[IX["Ncort"]],
        unfav=min(1.0, pdeath + (1.0 - pdeath) * seq),
    )


# ===========================================================================
# 시나리오 구성기
# ===========================================================================
def make_cfg(name, cef=True, cef_delay=0.0, cef_ci=False, van=False, rif=False,
             rif_lead=0.0, dex=False, dex_time=0.0, dex_days=4,
             mannitol=False, glycerol=False, MIC_cef=0.03, MIC_van=1.0,
             MIC_rif=0.03, N0=1.0e7, Nb0=1.0e3, host_def=1.0, mu_scale=1.0,
             antipyretic=0.0, anticonvulsant=0.0, vasopressor=0.0,
             drain=False, wt=70.0, tmax=336.0):
    reg = Regimen()
    if cef:
        if cef_ci:
            reg.iv("Cef_c", 2000.0, cef_delay, 12.0, 1, dur=0.5)      # 부하
            for k in range(14):                                        # 4 g/일 지속주입
                reg.infuse("Cef_c", 4000.0 / 24.0,
                           cef_delay + 0.5 + 24.0 * k, cef_delay + 0.5 + 24.0 * (k + 1))
        else:
            reg.iv("Cef_c", 2000.0, cef_delay, 12.0, 28, dur=0.5)
    if van:
        reg.iv("Van_c", 15.0 * wt, cef_delay, 6.0, 56, dur=1.0)
    if rif:
        reg.iv("Rif_c", 600.0, cef_delay - rif_lead, 12.0, 28, dur=0.5)
    if dex:
        reg.iv("Dex_c", 0.15 * wt, dex_time, 6.0, 4 * dex_days, dur=0.25)
    if mannitol:
        reg.iv("Mann_c", 0.5 * wt * 1000.0, 1.0, 6.0, 12, dur=0.5)     # 0.5 g/kg q6h ×3 d
    if glycerol:
        reg.po("Gly_gut", 0.25 * wt * 1000.0, 0.0, 6.0, 16)            # 6 g/kg/일 ×4 d
    shift = reg.shift_to_zero()          # [F33] t<0 투여를 살린다
    drain_fn = (lambda t: 12.0 if 2.0 <= t <= 96.0 else 0.0) if drain else None
    # [F23] 숙주방어 저하는 탐식뿐 아니라 억제 실패(순증식 상승)로도 나타난다
    mu_scale = mu_scale * (1.0 + 0.5 * (1.0 - host_def))
    return dict(name=name, reg=reg, MIC_cef=MIC_cef, MIC_van=MIC_van, MIC_rif=MIC_rif,
                N0=N0, Nb0=Nb0, host_def=host_def, mu_scale=mu_scale,
                antipyretic=antipyretic, anticonvulsant=anticonvulsant,
                vasopressor=vasopressor, csf_drain=drain_fn, tmax=tmax,
                shift=shift)


SCENARIOS = [
    ("S01 무치료 자연경과", dict(cef=False)),
    ("S02 CEF 단독 (t=0)", dict()),
    ("S03 CEF + DEX 20분 전", dict(dex=True, dex_time=-0.33)),
    ("S04 CEF + DEX 동시 (0 h)", dict(dex=True, dex_time=0.0)),
    ("S05 CEF + DEX +2 h", dict(dex=True, dex_time=2.0)),
    ("S06 CEF + DEX +4 h", dict(dex=True, dex_time=4.0)),
    ("S07 CEF + DEX +12 h", dict(dex=True, dex_time=12.0)),
    ("S08 CEF + DEX 1일만", dict(dex=True, dex_time=0.0, dex_days=1)),
    ("S09 항생제 3 h 지연", dict(cef_delay=3.0, dex=True, dex_time=3.0)),
    ("S10 항생제 6 h 지연", dict(cef_delay=6.0, dex=True, dex_time=6.0)),
    ("S11 항생제 12 h 지연", dict(cef_delay=12.0, dex=True, dex_time=12.0)),
    ("S12 내성균(MIC 4) CEF 단독", dict(MIC_cef=4.0)),
    ("S13 내성균 + VAN, DEX 없음", dict(MIC_cef=4.0, van=True)),
    ("S14 내성균 + VAN + DEX", dict(MIC_cef=4.0, van=True, dex=True, dex_time=0.0)),
    ("S15 내성균 + VAN + RIF + DEX", dict(MIC_cef=4.0, van=True, rif=True, dex=True, dex_time=0.0)),
    ("S16 RIF 2 h 선행 → CEF", dict(rif=True, rif_lead=2.0)),
    ("S17 CEF + DEX + 만니톨", dict(dex=True, dex_time=0.0, mannitol=True)),
    ("S18 CEF + DEX + 글리세롤", dict(dex=True, dex_time=0.0, glycerol=True)),
    ("S19 고균량 1e8 + DEX", dict(N0=1.0e8, dex=True, dex_time=0.0)),
    ("S20 조기 저균량 1e5 + DEX", dict(N0=1.0e5, dex=True, dex_time=0.0)),
    ("S21 고령/면역저하 + DEX", dict(host_def=0.4, dex=True, dex_time=0.0)),
    ("S22 CEF 지속주입 + DEX", dict(cef_ci=True, dex=True, dex_time=0.0)),
    ("S23 해열제만 (항생제 없음)", dict(cef=False, antipyretic=0.7)),
    ("S24 표준치료 CEF+VAN+DEX", dict(van=True, dex=True, dex_time=-0.33)),
    ("S25 내성균 CEF지속주입+VAN+DEX", dict(MIC_cef=4.0, cef_ci=True, van=True,
                                       dex=True, dex_time=-0.33)),   # [F24]
    ("S26 표준치료 + EVD 배액", dict(van=True, dex=True, dex_time=-0.33, drain=True)),
]

COHORT = [
    (1.0e6, 0.0, 1.00), (3.0e6, 2.0, 1.00), (1.0e7, 1.0, 0.90),
    (3.0e7, 4.0, 0.85), (1.0e8, 6.0, 0.70), (1.0e7, 3.0, 0.55),
    (5.0e6, 1.0, 1.00), (2.0e7, 8.0, 0.60), (8.0e7, 2.0, 0.80),
    (3.0e7, 24.0, 0.40),          # [F32] 늦게 온 환자 — 중증 난청의 꼬리를 만든다
]


def fmt(x, n=2):
    if x is None:
        return "n/a"
    return f"{x:.{n}f}"


def log10s(x):
    return math.log10(x) if x > 1.0 else 0.0


# ===========================================================================
# 메인
# ===========================================================================
def main():
    L = []
    W = L.append
    W("=" * 112)
    W("급성 세균성 수막염 (폐렴균) QSP 모델 — Python RK4 참조 구현 출력")
    W("63 ODE · 26 시나리오 · 가상코호트 10명 × 2 · abm_mrgsolve_model.R 와 동일 방정식")
    W("=" * 112)
    W("")

    # ---------------- [0] 구조적 정합성 ----------------
    W("[0] 손으로 검산되는 구조적 정합성 (모델을 돌리기 전에 맞아야 하는 것들)")
    W("-" * 112)
    W(f"  CSF 항정 ICP = P_ss + Qf·R_out")
    for r, tag in [(P["R_out0"], "정상"), (0.48, "발현시"), (1.20, "중증")]:
        W(f"     R_out={r:5.3f} ({tag:4s}) → ICP = {P['P_ss'] + P['Qf0']*r:5.1f} mmHg"
          f"   [문헌: 정상 7-15 · 수막염 20-40]")
    cf250, cf30 = (free_saturable(250.0, P["Bmax_cef"], P["Kd_cef"]),
                   free_saturable(30.0, P["Bmax_cef"], P["Kd_cef"]))
    W(f"  세프트리악손 포화성 결합: C_tot 30 → fu {cf30/30:.3f} ; 250 mg/L → fu {cf250/250:.3f}"
      f"   [보고: 0.05→0.20 상승]")
    for nm, PS, Eff, fu in (("CEF", P["PS_cef"], P["Eff_cef"], cf250 / 250),
                            ("VAN", P["PS_van"], P["Eff_van"], P["fu_van"]),
                            ("RIF", P["PS_rif"], P["Eff_rif"], P["fu_rif"]),
                            ("DEX", P["PS_dex"], 0.0, P["fu_dex"])):
        lip = (nm == "RIF" or nm == "DEX")
        def ratio(pb):
            pbe = 1.0 + P["a_rif"] * (pb - 1.0) if lip else pb
            return PS * pbe / (PS * pbe + P["Qf0"] + Eff)
        W(f"  {nm} CSF/유리혈장 항정비: Pb=1 {ratio(1):.3f} · Pb=4 {ratio(4):.3f} · Pb=9 {ratio(9):.3f}"
          f"  → 총농도 침투 {ratio(1)*fu:.3f}/{ratio(9)*fu:.3f}")
    W("     [보고 총농도 침투율: CEF 0.015-0.10 · VAN 0.01-0.15 · RIF 0.10-0.20 · DEX 0.15-0.30]")
    g = 60.0
    W(f"  CSF 포도당 60 mg/dL 균형: 반출유입 {P['Tmax_glc']*(P['Glc_pl']/(P['Km_glut']+P['Glc_pl'])-g/(P['Km_glut']+g)):.1f}"
      f" = bulk {g*P['Qf0']/P['Vcsf']:.1f} + 뇌 {P['q_brain']:.1f} mg/dL/h")
    W(f"  포도당 소비 경쟁 (PMN 4000/µL vs 세균 1e7/mL): {P['q_pmn']*4000:.1f} vs {P['q_bact']:.1f} mg/dL/h"
      f"  → 호중구가 {P['q_pmn']*4000/P['q_bact']:.0f}배 더 먹는다")
    W(f"  RK4 안정조건: ICP 방정식 λ_max = K_el·ICP/R_out = {P['K_el']*60/P['R_out0']:.1f}/h (ICP 60·회복기)"
      f" → dt < 2.785/λ = {2.785/(P['K_el']*60/P['R_out0']):.3f} h ; 사용 dt 0.04 h")
    W("")

    # ---------------- 시나리오 실행 ----------------
    rec_times = [0.5, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 48, 72, 120, 168, 240, 336]
    R = {}
    for label, kw in SCENARIOS:
        cfg = make_cfg(label, **kw)
        R[label] = simulate(cfg, tmax=cfg["tmax"], record=rec_times)
        print(f"  ... {label}", file=sys.stderr)

    # ---------------- [1] 자연경과 ----------------
    W("[1] S01 무치료 자연경과 — 진단 지표가 교과서 범위로 가는가")
    W("-" * 112)
    W("   t[h] log10Nc  PMN/µL 단백mg/dL 포도당 락테이트  Q_alb   TNF  IL6/1e3   Pb   ICP   CPP   MAP  체온")
    for (t, y) in R["S01 무치료 자연경과"]["series"]:
        if t > 120:
            continue
        W(f"  {t:5.0f} {log10s(y[IX['Nc']]):7.2f} {y[IX['PMN']]:7.0f} {y[IX['Prot_csf']]:8.0f}"
          f" {y[IX['Glc_csf']]:6.1f} {y[IX['Lac_csf']]:8.1f} {1000*y[IX['Alb_csf']]/P['Alb_ser']:6.1f}"
          f" {y[IX['TNF']]:5.0f} {y[IX['IL6']]/1000:8.1f} {y[IX['Pb']]:4.1f}"
          f" {y[IX['ICP']]:5.1f} {y[IX['MAP']]-y[IX['ICP']]:5.1f} {y[IX['MAP']]:5.1f} {y[IX['Temp']]:5.1f}")
    W("  [문헌 폐렴균 수막염: CSF 백혈구 1,000-5,000/µL · 단백 100-500 mg/dL · 포도당 <40(혈당비<0.4)")
    W("   · 락테이트 >3.5 mmol/L · Q_alb 30-100 · CSF TNF 100-1,000 pg/mL · ICP 종종 >20 mmHg]")
    W("")
    W("[1b] S24 표준치료 — 회복하는가 (위험함수가 만성적으로 적립되지 않는지 확인)")
    W("-" * 112)
    W("   t[h] log10Nc  PMN/µL 단백mg/dL 포도당 락테이트   TNF   Mg   Pb R_out   ICP   CPP   MAP  SOFA")
    for (t, y) in R["S24 표준치료 CEF+VAN+DEX"]["series"]:
        W(f"  {t:5.0f} {log10s(y[IX['Nc']]):7.2f} {y[IX['PMN']]:7.0f} {y[IX['Prot_csf']]:8.0f}"
          f" {y[IX['Glc_csf']]:6.1f} {y[IX['Lac_csf']]:8.1f} {y[IX['TNF']]:5.0f}"
          f" {y[IX['Mg']]:4.2f} {y[IX['Pb']]:4.1f} {y[IX['R_out']]:5.2f}"
          f" {y[IX['ICP']]:5.1f} {y[IX['MAP']]-y[IX['ICP']]:5.1f} {y[IX['MAP']]:5.1f}"
          f" {y[IX['SOFA']]:5.2f}")
    W("")

    # ---------------- [2] 용해 폭발 ----------------
    W("[2] 손상 flux 의 봉우리는 '첫 항생제 투여 직후'다  (S02 CEF 단독, DEX 없음)")
    W("-" * 112)
    W("   t[h] CEF_CSF C/MIC kkill/h log10Nc 용해flux  CW   PLY   Mg    TNF   PMN/µL")
    for (t, y) in R["S02 CEF 단독 (t=0)"]["series"]:
        if t > 24:
            continue
        kk = emax_kill(y[IX["Cef_csf"]], 0.03, P["Emax_cef"], P["EC50r_cef"], P["h_cef"])
        lf = (kk * P["Y_cef"] * y[IX["Nc"]] + P["k_autolysis"] * y[IX["Nc"]]) / 1e6
        W(f"  {t:5.1f} {y[IX['Cef_csf']]:7.2f} {y[IX['Cef_csf']]/0.03:6.0f} {kk:7.2f}"
          f" {log10s(y[IX['Nc']]):7.2f} {lf:8.2f} {y[IX['CW']]:5.1f} {y[IX['PLY']]:5.2f}"
          f" {y[IX['Mg']]:5.2f} {y[IX['TNF']]:6.0f} {y[IX['PMN']]:7.0f}")
    e2 = R["S02 CEF 단독 (t=0)"]["endpoints"]
    W(f"  최고 용해flux 는 투여 후 첫 시간 안에 나타나고, TNF 는 {e2['peak_TNF']:.0f} pg/mL 로 올라간다.")
    W("  [문헌: Mustafa 1989 등 — 항생제 투여 후 CSF TNF/IL-1 이 오히려 상승한다]")
    W("")

    # ---------------- [3] 덱사메타손 타이밍 ----------------
    W("[3] 덱사메타손 타이밍 — 방패는 봉우리보다 먼저 켜져 있어야 한다")
    W("-" * 112)
    W("  시나리오                    AUC_lysis AUC_TNF/1e3 peakTNF peakPMN peakPb 청력dB 인지z 사망% 불량%")
    for lab in ["S02 CEF 단독 (t=0)", "S03 CEF + DEX 20분 전", "S04 CEF + DEX 동시 (0 h)",
                "S05 CEF + DEX +2 h", "S06 CEF + DEX +4 h", "S07 CEF + DEX +12 h",
                "S08 CEF + DEX 1일만"]:
        e = R[lab]["endpoints"]
        W(f"  {lab:26s} {e['AUC_lysis']:9.1f} {e['AUC_TNF']/1000:11.1f} {e['peak_TNF']:7.0f}"
          f" {e['peak_PMN']:7.0f} {e['peak_Pb']:6.1f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    b, a = R["S02 CEF 단독 (t=0)"]["endpoints"], R["S03 CEF + DEX 20분 전"]["endpoints"]
    c4 = R["S06 CEF + DEX +4 h"]["endpoints"]
    W(f"  → 스테로이드 자체의 이득은 크다: 사망 {b['death']*100:.1f} % → {a['death']*100:.1f} %"
      f" (AUC_TNF −{(1-a['AUC_TNF']/b['AUC_TNF'])*100:.0f} %).")
    W(f"     그런데 **투여 시점의 효과는 이 조건에서 거의 없다**:"
      f" −20분 {a['death']*100:.1f} % vs +4 h {c4['death']*100:.1f} %.")
    W(f"     반대로 **투여 기간**은 분명하다: 4일 {a['death']*100:.1f} % vs 1일"
      f" {R['S08 CEF + DEX 1일만']['endpoints']['death']*100:.1f} %.")
    W("")
    W("  같은 스윕을 중증 환자에게 (균량 1e8 · 항생제 6 h 지연) — 폭발이 클 때는 어떤가")
    W("  DEX 시점(h)  AUC_lysis AUC_TNF/1e3 peakTNF peakPMN peakICP 청력dB 인지z 사망% 불량%")
    for dt in [None, -0.33, 0.0, 2.0, 4.0, 8.0, 12.0, 24.0]:
        if dt is None:
            cfg = make_cfg("sev-nodex", cef_delay=6.0, N0=1.0e8, van=True)
            tag = "  없음"
        else:
            cfg = make_cfg(f"sev{dt}", cef_delay=6.0, N0=1.0e8, van=True,
                           dex=True, dex_time=6.0 + dt)
            tag = f"{dt:+6.2f}"
        e = simulate(cfg, tmax=336.0)["endpoints"]
        W(f"  {tag:>11s} {e['AUC_lysis']:10.1f} {e['AUC_TNF']/1000:11.1f}"
          f" {e['peak_TNF']:7.0f} {e['peak_PMN']:7.0f} {e['peak_ICP']:7.1f}"
          f" {e['hear_dB']:6.1f} {e['cog_z']:6.2f} {e['death']*100:5.1f}"
          f" {e['unfav']*100:5.1f}")
        print(f"  ... severe sweep {tag}", file=sys.stderr)
    W("")
    W("  해석 — 이 모델이 예측하는 것과 예측하지 못하는 것:")
    W("  · 스테로이드를 주는 것 자체의 이득은 임상시험 크기로 재현된다 (아래 [9]).")
    W("  · 그러나 '첫 항생제 투여 **전**이어야 한다'는 부분은 이 모델에서 약하게")
    W("    나온다.  이유는 구조적이다: 사이토카인 폭발은 2-6 h 폭인데, 청력·인지·")
    W("    피질 손상 적분은 수십~수백 시간 규모다.  즉 봉우리가 14일 적분에서")
    W("    차지하는 비중이 작다.  모델은 대신 **투여 기간**(1일 vs 4일)과")
    W("    **항생제 지연**에 훨씬 민감하다.")
    W("  · 이것은 편리한 결론이 아니라 반증 가능한 예측이다.  빠진 기전의 가장 유력한")
    W("    후보는 **손상항의 문턱**이다: 이 모델의 모든 손상항은 구동자에 1차")
    W("    (dHC/dt ∝ driver)이므로, '3 h 동안 2배 높은 봉우리'와 '100 h 동안 1.06배")
    W("    높은 수준'이 같은 손상을 낸다.  실제 유모세포·뉴런 사멸에 문턱이 있다면")
    W("    봉우리가 불균형하게 비싸지고, 그때 비로소 '항생제 전'이 산술적으로")
    W("    중요해진다.  문턱을 넣어 권고를 재현하는 것은 쉽지만, 그것은 결론에")
    W("    맞춰 구조를 고르는 일이므로 여기서는 하지 않고 예측으로 남긴다.")
    W("")

    # ---------------- [4] 항생제 지연 ----------------
    W("[4] 항생제 지연 — 늦게 죽이면 더 많은 균을 죽여야 하므로 폭발이 커진다")
    W("-" * 112)
    W("  시나리오                    AUC_lysis peakICP minCPP 부종mL peakLac 청력dB 사망% 불량%")
    for lab in ["S04 CEF + DEX 동시 (0 h)", "S09 항생제 3 h 지연",
                "S10 항생제 6 h 지연", "S11 항생제 12 h 지연"]:
        e = R[lab]["endpoints"]
        W(f"  {lab:26s} {e['AUC_lysis']:9.1f} {e['peak_ICP']:7.1f} {e['min_CPP']:6.1f}"
          f" {e['peak_Vbr']:6.1f} {e['peak_Lac']:7.1f} {e['hear_dB']:6.1f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    W("  [문헌: door-to-antibiotic 지연 1시간마다 예후 악화 — Proulx 2005 (>6 h OR 8.4),")
    W("   Køster-Rasmussen 2008.  모델은 이것을 '더 큰 용해 폭발'로 기전적으로 재현한다]")
    W("")

    # ---------------- [5] 반코마이신 문 닫기 ----------------
    W("[5] 부호가 갈리는 곳 — 덱사메타손은 반코마이신이 들어갈 문을 함께 닫는다")
    W("-" * 112)
    W("  시나리오                        AUC_van_CSF C_van(48h) T>4×MIC T_cef>4×MIC peakPb"
      " 멸균t 부착제거t 사망%")
    for lab in ["S13 내성균 + VAN, DEX 없음", "S14 내성균 + VAN + DEX",
                "S15 내성균 + VAN + RIF + DEX", "S25 내성균 CEF지속주입+VAN+DEX",
                "S12 내성균(MIC 4) CEF 단독", "S24 표준치료 CEF+VAN+DEX"]:
        e = R[lab]["endpoints"]
        c48 = next((y[IX["Van_csf"]] for (tt, y) in R[lab]["series"] if abs(tt - 48) < 1e-6), 0.0)
        W(f"  {lab:30s} {e['AUC_van']:11.1f} {c48:10.2f} {e['T_van']:7.1f} {e['T_cef']:11.1f}"
          f" {e['peak_Pb']:6.1f} {fmt(e['sterile_t'],1):>6s} {fmt(e['adh_t'],1):>9s}"
          f" {e['death']*100:5.1f}")
    v13, v14 = R["S13 내성균 + VAN, DEX 없음"]["endpoints"], R["S14 내성균 + VAN + DEX"]["endpoints"]
    W(f"  → DEX 는 장벽을 조여(최고 Pb {v13['peak_Pb']:.1f} → {v14['peak_Pb']:.1f}) 반코마이신 CSF AUC 를"
      f" {(1-v14['AUC_van']/max(v13['AUC_van'],1e-9))*100:.0f} % 깎고,")
    W(f"     그 결과 CSF 멸균이 {fmt(v13['sterile_t'],1)} h → {fmt(v14['sterile_t'],1)} h 로 늦어진다."
      f"  같은 조작이 세프트리악손에서는 아무 해가 없다 (C/MIC 여유 200배).")
    W(f"     내성균에서 세프트리악손 단독은 4×MIC(16 mg/L)를 "
      f"{R['S12 내성균(MIC 4) CEF 단독']['endpoints']['T_cef']:.1f} h 만 넘고 멸균에 실패한다")
    W("     — 반코마이신을 더하는 임상 근거가 출력으로 나온다.")
    W("  [문헌: Paris 1994 (토끼, DEX 가 VAN CSF 침투 감소) vs Ricard 2007 (사람, 고용량에서 유지)")
    W("   → 모델은 '여유가 없다'는 쪽을 정량화한다]")
    W("")

    # ---------------- [6] 용해성 vs 비용해성 ----------------
    W("[6] 같은 log-kill, 다른 손상 적분 — 비용해성 선행투여(리팜핀)")
    W("-" * 112)
    W("  시나리오                    멸균t[h] AUC_lysis peakCW peakPLY peakTNF 청력dB 인지z 사망%")
    for lab in ["S02 CEF 단독 (t=0)", "S16 RIF 2 h 선행 → CEF", "S03 CEF + DEX 20분 전"]:
        e, pk = R[lab]["endpoints"], R[lab]["peak"]
        W(f"  {lab:26s} {fmt(e['sterile_t'],1):>8s} {e['AUC_lysis']:9.1f} "
          f"{max(y[IX['CW']] for _, y in R[lab]['series']):6.1f} {e['peak_PLY']:7.2f}"
          f" {e['peak_TNF']:7.0f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f} {e['death']*100:5.1f}")
    W("  [문헌: Nau/Böttcher — 비용해성 항생제 선행투여가 세포벽 방출과 CSF 염증을 줄인다]")
    W("")

    # ---------------- [7] 전체 시나리오 ----------------
    W("[7] 26 시나리오 종합")
    W("-" * 112)
    W("  시나리오                        멸균t 부착제거 peakICP minCPP 최저포도당 최고Lac Q_alb"
      " 청력dB 인지z 급성위험 구조위험 사망% 불량%")
    for label, _ in SCENARIOS:
        e = R[label]["endpoints"]
        W(f"  {label:30s} {fmt(e['sterile_t'],1):>6s} {fmt(e['adh_t'],1):>7s}"
          f" {e['peak_ICP']:7.1f} {e['min_CPP']:6.1f} {e['min_Glc']:10.1f} {e['peak_Lac']:7.1f}"
          f" {e['peak_Qalb']:5.0f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f}"
          f" {e['haz_acute']:8.3f} {e['haz_struct']:8.3f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    W("")

    # ---------------- [8] dt 절반 검산 ----------------
    W("[8] 적분 수렴 검산 (dt 절반) — 같은 시나리오, 스텝만 절반")
    W("-" * 112)
    W("  시나리오                       엔드포인트      dt=1.0     dt=0.5    상대차%")
    for lab, kw in [("S02 CEF 단독 (t=0)", dict()),
                    ("S24 표준치료 CEF+VAN+DEX", dict(van=True, dex=True, dex_time=-0.33))]:
        h = simulate(make_cfg(lab, **kw), tmax=336.0, dt_scale=0.5)["endpoints"]
        f = R[lab]["endpoints"]
        for key, nm in [("peak_ICP", "peak ICP"), ("hear_dB", "청력 dB"),
                        ("death", "사망확률"), ("AUC_lysis", "AUC_lysis")]:
            rel = abs(h[key] - f[key]) / max(abs(f[key]), 1e-12) * 100
            W(f"  {lab:26s} {nm:12s} {f[key]:10.4f} {h[key]:10.4f} {rel:9.3f}")
    W("")

    # ---------------- [9] 가상 코호트 ----------------
    W("[9] 가상 코호트 10명 × (DEX 유/무) — European Dexamethasone Study 대조")
    W("-" * 112)
    W("  환자 균량log10 지연h 숙주방어  사망%(DEX-) 사망%(DEX+) 청력dB(-) 청력dB(+) 불량%(-) 불량%(+)")
    dm, dp, hm, hp, um, up = [], [], [], [], [], []
    for i, (n0, dl, hd) in enumerate(COHORT, 1):
        a = simulate(make_cfg("-", cef_delay=dl, N0=n0, host_def=hd, van=True), 336.0)["endpoints"]
        b = simulate(make_cfg("+", cef_delay=dl, N0=n0, host_def=hd, van=True,
                              dex=True, dex_time=dl - 0.33), 336.0)["endpoints"]
        dm.append(a["death"]); dp.append(b["death"])
        hm.append(a["hear_dB"]); hp.append(b["hear_dB"])
        um.append(a["unfav"]); up.append(b["unfav"])
        W(f"  {i:4d} {math.log10(n0):8.1f} {dl:6.1f} {hd:8.2f} {a['death']*100:11.1f}"
          f" {b['death']*100:11.1f} {a['hear_dB']:9.1f} {b['hear_dB']:9.1f}"
          f" {a['unfav']*100:8.1f} {b['unfav']*100:8.1f}")
        print(f"  ... cohort {i}", file=sys.stderr)
    n = len(COHORT)
    W(f"  평균                        {sum(dm)/n*100:11.1f} {sum(dp)/n*100:11.1f}"
      f" {sum(hm)/n:9.1f} {sum(hp)/n:9.1f} {sum(um)/n*100:8.1f} {sum(up)/n*100:8.1f}")
    W("  [대조 de Gans & van de Beek NEJM 2002 폐렴균 아군: 사망 34 % → 14 %,")
    W("   불량결과 52 % → 26 %.  소아 Hib 중증 난청은 ~15 % → ~5 %]")
    W("")

    # ---------------- [10] 위험함수 분해 ----------------
    W("[10] 사망 위험함수의 성분별 분해 — 어느 생리이탈이 실제로 사람을 죽이는가")
    W("-" * 112)
    W("  Haz = h0·T + h_icp·I_icp + h_cpp·I_cpp + h_sofa·I_sofa + h_isch·I_isch")
    W("        + h_bact·I_bact + h_ncsf·I_ncsf  (급성)  +  h_cort·(1−N_cort)  (구조)")
    W("")
    W("  시나리오                        ICP항   CPP항  SOFA항  허혈항 균혈증항 CSF감염   기저   구조항"
      "   합계  사망%")
    for label in ["S01 무치료 자연경과", "S02 CEF 단독 (t=0)", "S03 CEF + DEX 20분 전",
                  "S11 항생제 12 h 지연", "S14 내성균 + VAN + DEX",
                  "S20 조기 저균량 1e5 + DEX", "S24 표준치료 CEF+VAN+DEX",
                  "S26 표준치료 + EVD 배액"]:
        e = R[label]["endpoints"]
        c = dict(icp=P["h_icp"] * e["I_icp"], cpp=P["h_cpp"] * e["I_cpp"],
                 sofa=P["h_sofa"] * e["I_sofa"], isch=P["h_isch"] * e["I_isch"],
                 bact=P["h_bact"] * e["I_bact"], ncsf=P["h_ncsf"] * e["I_ncsf"],
                 base=P["h0"] * 336.0,
                 struct=P["h_cort_final"] * (1.0 - e["Ncort_f"]))
        tot = sum(c.values())
        W(f"  {label:30s} {c['icp']:6.3f} {c['cpp']:7.3f} {c['sofa']:7.3f} {c['isch']:7.3f}"
          f" {c['bact']:8.3f} {c['ncsf']:7.3f} {c['base']:6.3f} {c['struct']:7.3f} {tot:6.3f}"
          f" {e['death']*100:6.1f}")
    W("  → 잘 치료된 환자에서 남는 위험은 대부분 전신 장기부전(SOFA)이고,")
    W("     치료가 늦거나 실패한 환자에서는 관류(CPP·허혈)와 구조손상 항이 지배한다.")
    W("")

    W("=" * 112)
    W("끝.  이 출력은 abm_reference_python.py 실행 결과이며 abm_mrgsolve_model.R 의")
    W("방정식·파라미터와 1:1 대응한다.  파일 머리말의 [F1]-[F32] 는 이 수치검증이")
    W("실제로 잡아낸 결함과 그 교정 내역이다.")
    W("=" * 112)

    txt = "\n".join(L)
    print(txt)
    with open("abm_reference_output.txt", "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
