## =============================================================================
##  abm_mrgsolve_model.R
##  Acute Bacterial Meningitis (Streptococcus pneumoniae)
##  급성 세균성 수막염 — 63-ODE QSP 모델 (mrgsolve)
##
##  -----------------------------------------------------------------------------
##  Provenance / 이 파일이 만들어진 방식
##  -----------------------------------------------------------------------------
##  이 저장소 환경에는 R 런타임이 없었다.  그래서 아래 63개 방정식은 먼저
##  의존성 없는 순수 Python RK4 로 독립 구현되어 (`abm_reference_python.py`)
##  실제로 적분되었고, README.md 에 인용된 모든 숫자는 그 적분의 출력이며
##  `abm_reference_output.txt` 는 그 실행 로그 전문이다.
##
##  그렇게 했더니 **32개의 실제 결함**이 드러났다.  각각을 고친 자리에는
##  "BUG FOUND" 주석을 남겼다 — QSP 모델을 읽는 사람은 어느 줄이 어려웠는지
##  알 권리가 있다.  전체 목록은 abm_reference_python.py 머리말의 [F1]-[F32] 에
##  있고, 그중 방정식 구조를 바꾼 것들만 여기 옮긴다:
##
##    F1   포도당 소비항에 기질 의존성이 없어 CSF 포도당이 0 으로 붕괴한 뒤에도
##         호중구와 세균이 "없는 포도당"을 계속 먹었다 → 락테이트 243 mmol/L.
##         락테이트는 포도당에서만 나온다.  소비항에 Michaelis-Menten 필수.
##    F14  혈중 세균에 수용력이 없어 N_b 가 336 h 동안 2e10 CFU/mL 로 발산하고
##         그것이 CSF 를 역오염시켰다.  logistic 수용력 추가.
##    F15  삼투압 환산에 만니톨 MW 182 대신 18.2 를 써서 0.5 g/kg 이 혈장
##         삼투압을 113 mOsm/kg 올렸다 (실제 11.3).  10배 단위오류.
##    F16  위험함수가 영구손상 (1 − N_cort) 을 14일 내내 적분해서, 회복한
##         환자도 사망확률 99 % 가 되었다.  구조손상은 종료시점 1회 평가.
##    F17  자동조절 공식이 CPP<45 에서 자동조절이 남은 환자를 압력수동 환자보다
##         더 허혈하게 만들었다.  자동조절은 결코 손해일 수 없다.
##    F28  ICP·CPP 위험항이 선형이면 "CPP 45 를 100 h" 와 "CPP 5 를 10 h" 가
##         같은 위험이 된다.  경증의 지속적 이탈이 위험함수를 지배해버렸다.
##         2차식으로 교체.
##    F33  t < 0 에 예정된 투여가 조용히 사라졌다 - 그리고 그것이 이 모델의
##         중심 주장을 검증 불가능하게 만들고 있었다.  "덱사메타손을 첫 항생제보다
##         20분 먼저"를 time = -0.33 으로 쓰면, 적분이 t = 0 에서 시작하므로 그
##         용량은 평가되지 않는다.  그래서 "먼저 준 팔"의 실제 첫 투여는 두 번째
##         용량이 들어가는 +5.67 h 였고, "타이밍이 무관하다"는 관찰은 생리가 아니라
##         일정표 버그였다.  증거는 출력에 있었다: 리팜핀 2 h 선행 시나리오의 최고
##         세포벽/용소/TNF 가 리팜핀 없는 시나리오와 소수점 두 자리까지 같았다.
##         -> zero_shift() 로 전체 일정을 가장 이른 사건이 t=0 이 되도록 옮긴다.
##    F29  세균 소멸 하한(floor)을 도함수 블록 안에 쓰면 안 된다 — 그것은
##         불연속 스위치이고 적분기를 멈춘다.  Python 참조구현은 적분 스텝
##         **사이**에서 하한을 적용하고, mrgsolve 판에서는 아예 상태를 건드리지
##         않고 $TABLE 의 보고 임계값(N_c < 10 CFU/mL = 배양 검출한계)으로만
##         멸균을 판정한다.  LSODA 는 하한 없이도 안정적으로 0 에 수렴한다.
##
##  -----------------------------------------------------------------------------
##  What the model is for / 이 모델의 구조적 약속
##  -----------------------------------------------------------------------------
##  단 하나의 구조적 약속: **뇌를 다치게 하는 양은 균 수가 아니라 곱이다.**
##
##      injury flux  =  k_kill(t) · N(t) · Y_lysis · (1 − E_dex(t))
##
##  세 인자 모두가 시간의 함수이고, k_kill 은 첫 항생제 투여에서 0 → E_max 로
##  점프하는데 N 은 정확히 그때 최대다.  따라서 이 곱은 **치료의 첫 몇 시간에
##  최대가 된다.**  이 하나의 사실이 치료 지도 전체를 정렬한다:
##
##    · 스테로이드는 그 봉우리보다 먼저 켜져 있어야 한다 (그래서 "항생제 전 또는
##      함께").  전사효과 구획의 비대칭(켜짐 t½ 1.5 h / 꺼짐 17 h)이 그 창을 만든다.
##    · 항생제를 늦추면 N 이 지수적으로 커져 있으므로 봉우리가 커진다.
##      "빨리 주라"와 "천천히 죽여라"는 충돌하지 않는다 — 최적해는
##      **방패를 올린 상태로 최대한 빨리** 다.
##    · 같은 log-kill 이라도 Y 가 다르면 손상 적분이 다르다 (β-락탐 vs 리팜핀).
##
##  그리고 두 개의 문은 사실 같은 문이다 — 혈액-CSF 장벽 Pb.  염증이 열고
##  덱사메타손이 닫는다.  열리면 친수성 항생제가 들어오고(cluster 13) 동시에
##  알부민·호중구·물이 들어온다(cluster 6-9).  그래서 덱사메타손의 **부호는
##  약에 따라 갈린다**: 감수성 균주에 대한 세프트리악손은 C/MIC 여유가 200배라
##  순이득이지만, 세팔로스포린 내성균에 대한 반코마이신은 여유가 없어 같은
##  조작이 CSF 멸균을 2배 늦춘다.
##
##  -----------------------------------------------------------------------------
##  Compartment map (63 ODEs)
##  -----------------------------------------------------------------------------
##     1– 3  ceftriaxone: central, peripheral, CSF
##     4– 6  vancomycin: central, peripheral, CSF
##     7– 8  rifampin: central, CSF
##     9–11  dexamethasone: plasma, CSF, transcriptional effect compartment
##    12–15  osmotherapy: mannitol plasma, glycerol gut + plasma, brain osmoles
##    16–20  bacteria: CSF free, adherent/sequestered, blood, cell wall, pneumolysin
##    21–30  innate: microglia, TNF, IL-1b, IL-6, IL-10, CXCL8, complement,
##           neutrophils, MMP-9, ROS
##    31–40  barrier and CSF: Pb, CSF albumin, CSF protein, glucose, lactate,
##           net CSF accumulation, outflow resistance, ICP, brain water,
##           autoregulation
##    41–44  systemic: MAP, temperature, SOFA, volume status
##    45–50  injury: cortical neurones, dentate gyrus, hair cells, ossificans,
##           seizure burden, cumulative acute hazard
##    51–57  exposure trackers: AUC_cef, AUC_van, AUC_lysis, AUC_TNF,
##           T>4xMIC cef, T>4xMIC van, AUC_ICP
##    58–63  hazard component integrals (I_icp, I_cpp, I_sofa, I_isch, I_bact,
##           I_ncsf) — added so the hazard coefficients can be SOLVED
##           arithmetically against the trial targets (pneumococcal mortality
##           34% without and 14% with dexamethasone) instead of eyeballed.
##           The solution was a single scale factor of 0.30 on every weight.
##
##  Dosing compartments
##  -------------------
##    cmt 1  = CEF_C   ceftriaxone,   IV infusion (2 g over 0.5 h, q12h)
##    cmt 4  = VAN_C   vancomycin,    IV infusion (15 mg/kg over 1 h, q6h)
##    cmt 7  = RIF_C   rifampin,      IV infusion (600 mg over 0.5 h, q12h)
##    cmt 9  = DEX_C   dexamethasone, IV (0.15 mg/kg over 0.25 h, q6h x 4 d)
##    cmt 12 = MANN_C  mannitol,      IV (0.5 g/kg over 0.5 h, q6h)
##    cmt 13 = GLY_GUT glycerol,      oral (1.5 g/kg q6h)
##
##  Time unit is HOURS throughout.  t = 0 is PRESENTATION (hospital arrival),
##  which is already 12–36 h into the illness — hence the non-naive initial
##  conditions in $MAIN.
## =============================================================================

library(mrgsolve)
library(dplyr)

abm_code <- '
$PARAM @annotated
// ---------------------------- ceftriaxone PK -------------------------------
V1_cef   :  8.0  : Central volume (L)
V2_cef   :  8.0  : Peripheral volume (L)
CL_cef   :  1.0  : Clearance (L/h)
Q_cef    :  1.0  : Intercompartmental clearance (L/h)
Bmax_cef : 333.0 : Albumin binding capacity (mg/L)
Kd_cef   : 25.0  : Albumin binding Kd (mg/L)
PS_cef   :  3.0  : CSF diffusion clearance at Pb=1 (mL/h)
Eff_cef  :  0.0  : CSF active efflux (mL/h)
MIC_cef  :  0.03 : MIC (mg/L) - susceptible pneumococcus
Emax_cef :  1.40 : Max kill rate (1/h) = 0.61 log10 CFU/mL/h
EC50r_cef:  4.0  : C/MIC at half-maximal kill
h_cef    :  2.0  : Hill coefficient
Y_cef    :  1.00 : Lysis yield (reference, lytic killing)

// ---------------------------- vancomycin PK --------------------------------
V1_van   : 20.0  : Central volume (L)
V2_van   : 30.0  : Peripheral volume (L)
CL_van   :  4.0  : Clearance (L/h)
Q_van    :  8.0  : Intercompartmental clearance (L/h)
fu_van   :  0.50 : Free fraction
PS_van   :  1.4  : CSF diffusion clearance at Pb=1 (mL/h)
Eff_van  : 10.0  : CSF active efflux (mL/h)
MIC_van  :  1.0  : MIC (mg/L)
Emax_van :  1.20 : Max kill rate (1/h)
EC50r_van:  4.0  : C/MIC at half-maximal kill
h_van    :  1.5  : Hill coefficient
Y_van    :  0.90 : Lysis yield

// ---------------------------- rifampin PK ----------------------------------
V_rif    : 50.0  : Volume (L)
CL_rif   : 12.0  : Clearance (L/h)
fu_rif   :  0.20 : Free fraction
PS_rif   : 30.0  : CSF diffusion clearance (mL/h) - lipophilic
a_rif    :  0.10 : Sensitivity of rifampin CSF entry to Pb (low)
Eff_rif  :  0.0  : CSF active efflux (mL/h)
MIC_rif  :  0.03 : MIC (mg/L)
Emax_rif :  0.70 : Max kill rate (1/h)
EC50r_rif:  2.0  : C/MIC at half-maximal kill
h_rif    :  1.0  : Hill coefficient
Y_rif    :  0.15 : Lysis yield - NON-lytic killing

// ---------------------------- dexamethasone --------------------------------
V_dex    : 70.0  : Volume (L)
CL_dex   : 16.0  : Clearance (L/h)
fu_dex   :  0.32 : Free fraction
PS_dex   : 32.0  : CSF diffusion clearance (mL/h)
ktr_on   :  0.46 : Transcriptional effect onset rate (1/h) - t1/2 1.5 h
ktr_off  :  0.040: Transcriptional effect offset rate (1/h) - t1/2 17 h
Imax_dex :  0.80 : Max transcriptional inhibition
IC50_TR  :  0.010: Effect-compartment IC50 (mg/L)
dex_pb   :  1.5  : Direct barrier-tightening effect of dexamethasone

// ---------------------------- osmotherapy ----------------------------------
V_mann   : 17.0  : Mannitol volume (L)
CL_mann  :  6.0  : Mannitol clearance (L/h)
ka_gly   :  1.2  : Glycerol absorption rate (1/h)
V_gly    : 42.0  : Glycerol volume (L)
CL_gly   :  8.0  : Glycerol clearance (L/h)
MW_mann  : 182.0 : Mannitol molecular weight (g/mol)
MW_gly   : 92.0  : Glycerol molecular weight (g/mol)
k_osm    :  1.20 : Brain water removal per mOsm gradient (mL/h per mOsm)
k_osm_leak: 0.010: Osmole leak into brain (1/h)
k_osm_br_out: 0.030 : Brain osmole washout (1/h)

// ---------------------------- CSF physics ----------------------------------
Vcsf     : 150.0 : CSF volume (mL)
Qf0      : 21.0  : CSF formation rate (mL/h) = 0.35 mL/min
P_ss     :  6.0  : Sagittal sinus pressure (mmHg)
R_out0   :  0.167: Baseline CSF outflow resistance (mmHg per mL/h)
R_out_max:  1.20 : Maximal outflow resistance (mmHg per mL/h)
k_ro_on  :  0.09 : Outflow resistance rise rate
k_ro_off :  0.020: Outflow resistance recovery rate (1/h)
K_el     :  0.092: Intracranial elastance coefficient (1/mL) - PVI 25 mL

// ---------------------------- bacteria -------------------------------------
mu_max   :  0.85 : Max growth rate (1/h) - doubling 49 min
K_glc    : 10.0  : Glucose Km for growth (mg/dL)
Nmax     :  1e9  : CSF carrying capacity (CFU/mL)
k_adh    :  0.030: Adhesion rate to ependyma/plexus (1/h)
k_des    :  0.010: Desorption rate (1/h)
prot_adh :  0.35 : Fraction of kill rate reaching adherent pool
mu_adh   :  0.30 : Growth rate of adherent pool (1/h)
k_shed   :  2e-4 : CSF to blood shedding (1/h)
k_seed   :  0.010: Blood to CSF seeding (per CFU/mL per h)
mu_b     :  0.45 : Blood growth rate (1/h)
k_clr_b  :  0.40 : Splenic/hepatic clearance (1/h)
Nb_max   :  1e8  : Bloodstream carrying capacity (CFU/mL)
kill_b_boost: 1.0 : Antibiotic kill applied to bloodstream
kphag_max:  0.15 : Max CSF phagocytic kill (1/h) - opsonin-poor space
K_pmn_ph : 2500.0: PMN for half-maximal phagocytosis (cells/uL)
K_comp   :  1.0  : Complement for half-maximal opsonisation
k_autolysis: 0.008 : Spontaneous LytA autolysis rate (1/h)
Y_auto   :  1.00 : Lysis yield of autolysis
Y_phag   :  0.30 : Lysis yield of phagocytic killing

// ---------------------------- bacterial products ---------------------------
kCW_cl   :  0.060: Cell wall clearance (1/h) - t1/2 12 h
yPLY     :  0.50 : Pneumolysin yield per lysis unit
kPLY     :  0.35 : Pneumolysin clearance (1/h)
ply_secr :  0.006: Pneumolysin release by live organisms

// ---------------------------- innate immunity ------------------------------
kmg_on   :  1.10 : Microglial activation rate
kmg_off  :  0.070: Microglial deactivation rate (1/h)
K_cw     : 10.0  : Cell wall for half-maximal PAMP signal (CWU/mL)
K_nl     :  1.20 : Live organism signal constant (1e6 CFU/mL)
K_ply    :  1.50 : Pneumolysin signal constant
kTNF     : 1400.0: TNF production rate (pg/mL/h)
kel_TNF  :  0.70 : TNF elimination (1/h)
K_tnf    : 250.0 : TNF EC50 (pg/mL)
kIL1     : 320.0 : IL-1b production rate
kel_IL1  :  0.35 : IL-1b elimination (1/h)
K_il1    : 120.0 : IL-1b EC50 (pg/mL)
kIL6     : 9000.0: IL-6 production rate
kel_IL6  :  0.25 : IL-6 elimination (1/h)
K_il6    : 4000.0: IL-6 EC50 (pg/mL)
kIL10    : 260.0 : IL-10 production rate
kel_IL10 :  0.20 : IL-10 elimination (1/h)
K_il10   : 600.0 : IL-10 inhibition constant (pg/mL)
kC8      : 5000.0: CXCL8 production rate
kel_C8   :  0.30 : CXCL8 elimination (1/h)
K_c8     : 1500.0: CXCL8 EC50 (pg/mL)
kcomp    :  1.20 : Complement activation rate
kel_comp :  0.25 : Complement decay (1/h)
k_influx : 500.0 : PMN influx rate constant (cells/uL/h)
k_egress :  0.055: PMN egress (1/h)
k_apop   :  0.045: PMN apoptosis (1/h)
kMMP     : 60.0  : MMP-9 production (ng/mL/h)
kel_MMP  :  0.087: MMP-9 elimination (1/h) - t1/2 8 h
K_mmp    : 350.0 : MMP-9 EC50 (ng/mL)
kROS     :  1.30 : ROS production rate
kel_ROS  :  0.90 : ROS clearance (1/h)
K_ros    :  0.60 : ROS EC50

// ---------------------------- barrier --------------------------------------
Pb_max   : 20.0  : Maximum permeability multiplier
k_pb     :  0.020: Barrier opening rate
k_pb_off :  0.030: Barrier repair rate (1/h)
PS_alb   :  0.105: CSF albumin transfer clearance (mL/h at Pb=1)
Alb_ser  : 42000.0 : Serum albumin (mg/L)
PS_prot  :  0.090: CSF total protein transfer clearance (mL/h at Pb=1)
Prot_ser : 70000.0 : Serum total protein (mg/L)

// ---------------------------- glucose / lactate ----------------------------
Tmax_glc : 74.0  : GLUT1 max transport (mg/dL/h)
Km_glut  : 90.0  : GLUT1 Km (mg/dL)
Glc_pl   : 100.0 : Plasma glucose (mg/dL)
q_pmn    :  0.0050 : Glucose use per PMN (mg/dL/h per cell/uL)
q_bact   :  2.0  : Glucose use per 1e7 CFU/mL (mg/dL/h)
q_brain  :  1.0  : Brain glucose draw from CSF (mg/dL/h)
Km_use   :  5.0  : Michaelis constant for glucose utilisation (mg/dL)
inflam_glut: 0.30: Fractional reduction of GLUT1 by inflammation
k_lac    :  0.500: CSF lactate clearance (1/h)
Lac_base :  1.6  : Baseline CSF lactate (mmol/L)

// ---------------------------- oedema / perfusion ---------------------------
k_vas    :  1.05 : Vasogenic oedema rate (mL/h)
Vbr_max  : 60.0  : Maximum brain water excess (mL)
k_cyt    :  2.20 : Cytotoxic oedema rate (mL/h per unit ischaemia)
k_vbr_res:  0.020: Oedema resolution (1/h)
k_ar_loss:  0.16 : Autoregulation loss rate
k_ar_rec :  0.012: Autoregulation recovery rate (1/h)
CBF0     : 50.0  : Baseline CBF (mL/100g/min)
CPP0     : 75.0  : Reference CPP (mmHg)
CBF_crit :  0.55 : Ischaemic CBF threshold as fraction of CBF0

// ---------------------------- systemic -------------------------------------
MAP0     : 88.0  : Baseline MAP (mmHg)
k_map    :  0.050: Sepsis-driven MAP fall
k_map_rec:  0.10 : MAP recovery rate (1/h)
k_cush   :  0.30 : Cushing response gain
k_temp   :  0.85 : Fever production
k_temp_off: 0.35 : Temperature recovery (1/h)
Temp0    : 37.0  : Baseline temperature (C)
k_sofa   :  0.55 : SOFA accrual
k_sofa_rec: 0.10 : SOFA recovery (1/h)

// ---------------------------- injury ---------------------------------------
k_isch   :  0.055 : Ischaemic cortical loss rate (1/h)
k_ros_cort: 0.0008: ROS-driven cortical loss rate (1/h)
k_apo_dg :  0.0020: Dentate apoptosis rate (1/h)
k_ply_dg :  0.5   : Pneumolysin weight on dentate apoptosis
k_ros_dg :  0.4   : ROS weight on dentate apoptosis
k_hc     :  0.0022: Hair cell loss rate (1/h)
w_hc_ply :  0.5   : Pneumolysin weight on hair cell loss
w_hc_ros :  0.7   : ROS weight on hair cell loss
w_hc_pmn :  0.4   : PMN weight on hair cell loss
K_ply2   :  1.0   : Pneumolysin EC50 for injury
k_oss    :  0.0020: Labyrinthitis ossificans rate (1/h)
hc_reserve: 0.22  : Cochlear reserve before threshold shift
hc_span  :  0.40  : Hair cell loss span mapping to 0-120 dB
hc_exp   :  0.90  : Audiometric curve exponent
k_sz     :  0.10  : Seizure burden accrual
k_sz_off :  0.12  : Seizure burden decay (1/h)

// ---------------------------- mortality hazard -----------------------------
h0       :  1.2e-5 : Baseline hazard (1/h)
h_icp    :  6.0e-4 : ICP hazard weight (quadratic)
h_cpp    :  6.0e-4 : CPP hazard weight (quadratic)
h_sofa   :  6.0e-5 : SOFA hazard weight
h_isch   :  6.0e-4 : Ischaemia hazard weight (quadratic)
h_bact   :  4.5e-5 : Bacteraemia hazard weight
h_ncsf   :  1.5e-4 : Persistent CSF infection hazard weight
h_cort_final: 0.48 : Structural (cortical loss) hazard, evaluated once at end

// ---------------------------- scenario switches ----------------------------
HOST_DEF :  1.0  : Host defence index 0.4-1.0
MU_SCALE :  1.0  : Growth scaling (raised when containment fails)
ANTIPYR  :  0.0  : Antipyretic effect 0-1
ANTICONV :  0.0  : Anticonvulsant effect 0-1
VASOPRESS:  0.0  : Vasopressor support (mmHg/h)
CSF_DRAIN:  0.0  : External ventricular drainage (mL/h)
N0       :  1e7  : Initial CSF bacterial density (CFU/mL)
NB0      :  1e3  : Initial bacteraemia (CFU/mL)

$CMT @annotated
CEF_C    : Ceftriaxone central (mg)
CEF_P    : Ceftriaxone peripheral (mg)
CEF_CSF  : Ceftriaxone CSF (mg/L)
VAN_C    : Vancomycin central (mg)
VAN_P    : Vancomycin peripheral (mg)
VAN_CSF  : Vancomycin CSF (mg/L)
RIF_C    : Rifampin central (mg)
RIF_CSF  : Rifampin CSF (mg/L)
DEX_C    : Dexamethasone plasma (mg)
DEX_CSF  : Dexamethasone CSF (mg/L)
DEX_TR   : Dexamethasone transcriptional effect (mg/L equivalent)
MANN_C   : Mannitol plasma (mg)
GLY_GUT  : Glycerol gut (mg)
GLY_C    : Glycerol plasma (mg)
OSM_BR   : Brain accumulated osmoles (mOsm/kg equivalent)
NC       : CSF free bacteria (CFU/mL)
NADH     : Adherent/sequestered bacteria (CFU/mL)
NB       : Bloodstream bacteria (CFU/mL)
CW       : CSF cell wall pool (CWU/mL)
PLY      : CSF pneumolysin (ng/mL)
MG       : Microglial activation (0-1)
TNF      : CSF TNF-alpha (pg/mL)
IL1      : CSF IL-1beta (pg/mL)
IL6      : CSF IL-6 (pg/mL)
IL10     : CSF IL-10 (pg/mL)
CXCL8    : CSF CXCL8/IL-8 (pg/mL)
COMP     : CSF complement activation (arb)
PMN      : CSF neutrophils (cells/uL)
MMP9     : CSF MMP-9 (ng/mL)
ROS      : CSF oxidative burden (arb)
PB       : Barrier permeability multiplier
ALB_CSF  : CSF albumin (mg/L)
PROT_CSF : CSF total protein (mg/dL)
GLC_CSF  : CSF glucose (mg/dL)
LAC_CSF  : CSF lactate (mmol/L)
VCSF_NET : Cumulative net CSF accumulation (mL) - tracker
R_OUT    : CSF outflow resistance (mmHg per mL/h)
ICP      : Intracranial pressure (mmHg)
VBR      : Brain water excess (mL)
AUTOR    : Autoregulation integrity (0-1)
MAP      : Mean arterial pressure (mmHg)
TEMP     : Core temperature (C)
SOFA     : Organ dysfunction score
VOL      : Volume status multiplier
NCORT    : Cortical neurone viability (0-1)
NDG      : Dentate gyrus neurone viability (0-1)
HC       : Cochlear hair cell viability (0-1)
OSS      : Labyrinthitis ossificans (0-1)
SZ       : Seizure burden
HAZ      : Cumulative acute mortality hazard
AUC_CEF  : AUC of CSF ceftriaxone (mg*h/L)
AUC_VAN  : AUC of CSF vancomycin (mg*h/L)
AUC_LYS  : AUC of lysis flux (CWU/mL)
AUC_TNF  : AUC of CSF TNF (pg*h/mL)
T_CEF    : Time ceftriaxone CSF > 4xMIC (h)
T_VAN    : Time vancomycin CSF > 4xMIC (h)
AUC_ICP  : AUC of ICP (mmHg*h)
I_ICP    : Integral of squared ICP excess
I_CPP    : Integral of squared CPP deficit
I_SOFA   : Integral of SOFA
I_ISCH   : Integral of squared ischaemia
I_BACT   : Integral of log10 bacteraemia
I_NCSF   : Integral of log10 CSF bacterial density

$GLOBAL
#define SQ(x) ((x)*(x))

// Saturable albumin binding: solve Cf + Bmax*Cf/(Kd+Cf) = Ctot for the
// positive root.  Ceftriaxone binding is saturable, so free fraction RISES
// with dose (0.076 at 30 mg/L to 0.167 at 250 mg/L) - this matters because
// only free drug crosses into CSF.
double free_sat(double Ctot, double Bmax, double Kd) {
  if (Ctot <= 0.0) return 0.0;
  double b = Kd + Bmax - Ctot;
  return (-b + sqrt(b*b + 4.0*Kd*Ctot)) / 2.0;
}

double hillf(double x, double K) {
  if (x <= 0.0) return 0.0;
  return x / (K + x);
}

double hilln(double x, double K, double n) {
  if (x <= 0.0) return 0.0;
  double xn = pow(x, n);
  return xn / (pow(K, n) + xn);
}

// Emax kill on the C/MIC ratio.  EC50r = 4 reproduces the clinical
// "concentration must exceed ~10x MBC in CSF" rule of thumb.
double emax_kill(double C, double MIC, double Emax, double EC50r, double h) {
  if (C <= 0.0 || MIC <= 0.0) return 0.0;
  double r = C / MIC;
  double rh = pow(r, h);
  return Emax * rh / (pow(EC50r, h) + rh);
}

double clampd(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

$MAIN
// t = 0 is PRESENTATION, not infection.  Symptoms have been running 12-36 h,
// so the patient already has an established purulent exudate.  These initial
// conditions are set to the reported CSF panel of pneumococcal meningitis.
NC_0       = N0;
NADH_0     = 0.05 * N0;
NB_0       = NB0;
MG_0       = 0.55;
TNF_0      = 420.0;
IL1_0      = 180.0;
IL6_0      = 2.5e4;
IL10_0     = 260.0;
CXCL8_0    = 4.0e3;
COMP_0     = 0.9;
PMN_0      = 1200.0;
MMP9_0     = 160.0;
ROS_0      = 0.55;
PB_0       = 6.0;
ALB_CSF_0  = 1400.0;     // Q_alb about 33
PROT_CSF_0 = 180.0;
GLC_CSF_0  = 34.0;
LAC_CSF_0  = 5.2;
R_OUT_0    = 0.48;
ICP_0      = 18.0;
VBR_0      = 4.0;
AUTOR_0    = 0.75;
MAP_0      = 88.0;
TEMP_0     = 39.1;
SOFA_0     = 2.0;
VOL_0      = 1.0;
NCORT_0    = 1.0;
NDG_0      = 1.0;
HC_0       = 1.0;
OSS_0      = 0.0;

$ODE
// ===========================================================================
// Plasma free concentrations
// ===========================================================================
double Cef_pl = CEF_C / V1_cef;
double Cef_f  = free_sat(Cef_pl, Bmax_cef, Kd_cef);
double Van_pl = VAN_C / V1_van;
double Van_f  = fu_van * Van_pl;
double Rif_pl = RIF_C / V_rif;
double Rif_f  = fu_rif * Rif_pl;
double Dex_pl = DEX_C / V_dex;
double Dex_f  = fu_dex * Dex_pl;

// Hydrophilic drugs see the full barrier state; lipophilic drugs barely do.
// This asymmetry is the whole reason dexamethasone can shut vancomycin out
// while leaving rifampin (and itself) untouched.
double Pb_hyd = PB;
double Pb_lip = 1.0 + a_rif * (PB - 1.0);
double Qbulk  = Qf0;                       // CSF bulk flow acts as the sink

// ===========================================================================
// Drug PK
// ===========================================================================
dxdt_CEF_C   = -CL_cef*Cef_pl - Q_cef*(Cef_pl - CEF_P/V2_cef);
dxdt_CEF_P   =  Q_cef*(Cef_pl - CEF_P/V2_cef);
dxdt_CEF_CSF = (PS_cef*Pb_hyd*(Cef_f - CEF_CSF) - (Qbulk + Eff_cef)*CEF_CSF)/Vcsf;

dxdt_VAN_C   = -CL_van*Van_pl - Q_van*(Van_pl - VAN_P/V2_van);
dxdt_VAN_P   =  Q_van*(Van_pl - VAN_P/V2_van);
dxdt_VAN_CSF = (PS_van*Pb_hyd*(Van_f - VAN_CSF) - (Qbulk + Eff_van)*VAN_CSF)/Vcsf;

dxdt_RIF_C   = -CL_rif*Rif_pl;
dxdt_RIF_CSF = (PS_rif*Pb_lip*(Rif_f - RIF_CSF) - (Qbulk + Eff_rif)*RIF_CSF)/Vcsf;

dxdt_DEX_C   = -CL_dex*Dex_pl;
dxdt_DEX_CSF = (PS_dex*Pb_lip*(Dex_f - DEX_CSF) - Qbulk*DEX_CSF)/Vcsf;

// Asymmetric transcriptional effect compartment.  Glucocorticoid repression
// switches ON fast (mRNA destabilisation, IkB induction: t1/2 ~1.5 h) and OFF
// slowly (t1/2 ~17 h).  That asymmetry IS the timing window: a dose given
// before the lysis burst is already working when the burst arrives, and a dose
// given after it cannot retroactively cover it.
double ktr = (DEX_CSF > DEX_TR) ? ktr_on : ktr_off;
dxdt_DEX_TR = ktr * (DEX_CSF - DEX_TR);
double Idex = Imax_dex * hillf(DEX_TR, IC50_TR);

// ===========================================================================
// Osmotherapy
// ===========================================================================
dxdt_MANN_C  = -CL_mann*(MANN_C/V_mann);
double gly_abs = ka_gly * GLY_GUT;
dxdt_GLY_GUT = -gly_abs;
dxdt_GLY_C   =  gly_abs - CL_gly*(GLY_C/V_gly);

// BUG FOUND (F15): this conversion originally divided by 18.2 and 9.2 instead
// of the molecular weights 182 and 92, so 0.5 g/kg of mannitol raised plasma
// osmolality by 113 mOsm/kg - physiologically impossible.  mg/L / (g/mol) =
// mmol/L = mOsm/kg equivalent.  Correct value is 11.3 mOsm/kg.
double Osm_pl   = (MANN_C/V_mann)/MW_mann + (GLY_C/V_gly)/MW_gly;
double Osm_grad = (Osm_pl > OSM_BR) ? (Osm_pl - OSM_BR) : 0.0;
dxdt_OSM_BR = k_osm_leak*(1.0 + 0.35*(PB - 1.0))*Osm_pl - k_osm_br_out*OSM_BR;

// ===========================================================================
// Bacteria and the kill term
// ===========================================================================
double kill_cef = emax_kill(CEF_CSF, MIC_cef, Emax_cef, EC50r_cef, h_cef);
double kill_van = emax_kill(VAN_CSF, MIC_van, Emax_van, EC50r_van, h_van);
double kill_rif = emax_kill(RIF_CSF, MIC_rif, Emax_rif, EC50r_rif, h_rif);
double kill_tot = kill_cef + kill_van + kill_rif;

// BUG FOUND (F1): every glucose-consuming term needs this substrate factor.
// Without it, consumption continued after CSF glucose hit zero and the model
// manufactured lactate out of nothing (243 mmol/L against a measured 6-12).
double f_glc = hillf(GLC_CSF, Km_use);

double mu = mu_max * hillf(GLC_CSF, K_glc) * (1.0 - NC/Nmax) * MU_SCALE;
double k_phag = kphag_max * hillf(PMN, K_pmn_ph)
                * (0.2 + 0.8*hillf(COMP, K_comp)) * HOST_DEF;

dxdt_NC   = (mu - kill_tot - k_phag)*NC + k_des*NADH - k_adh*NC
            + k_seed*NB - k_shed*NC;
dxdt_NADH = k_adh*NC - k_des*NADH
            + (mu_adh*hillf(GLC_CSF, K_glc) - kill_tot*prot_adh)*NADH;

// BUG FOUND (F14): without the logistic capacity term, bloodstream bacteria
// grew without bound (2e10 CFU/mL by day 14) and the k_seed term re-infected
// the CSF from an impossible reservoir.
dxdt_NB   = (mu_b*(1.0 - NB/Nb_max) - k_clr_b*HOST_DEF - kill_b_boost*kill_tot)*NB
            + k_shed*NC;

// ---------------------------------------------------------------------------
// THE PRODUCT.  This single expression is the model's central claim: what
// injures the brain is the RATE of bacterial destruction times the yield per
// destroyed organism, not the bacterial count.  k_kill jumps from 0 to E_max
// at the first dose while NC is still maximal, so this flux PEAKS IN THE
// FIRST HOURS OF TREATMENT.  Y_rif = 0.15 x Y_cef is why a non-lytic agent
// buys the same log-kill at one seventh of the inflammatory cost.
// ---------------------------------------------------------------------------
double Y_adh = ((kill_cef + kill_van) > kill_rif) ? Y_cef : Y_rif;
double lysis_flux = (kill_cef*Y_cef + kill_van*Y_van + kill_rif*Y_rif)*NC
                  + kill_tot*prot_adh*Y_adh*NADH
                  + k_autolysis*Y_auto*NC
                  + k_phag*Y_phag*NC;

dxdt_CW  = lysis_flux/1.0e6 - kCW_cl*CW;
dxdt_PLY = yPLY*(lysis_flux + ply_secr*NC)/1.0e6 - kPLY*PLY;

// ===========================================================================
// Perfusion (needed by ROS and injury terms)
// ===========================================================================
double CPP = MAP - ICP;
// BUG FOUND (F17): the original form was
//   CBF = CBF0*(AUTOR*plateau + (1-AUTOR)*CPP/CPP0)
// with plateau = 0 below CPP 45.  That made autoregulation HARMFUL: a patient
// with intact autoregulation computed a lower CBF than a pressure-passive one
// (ischaemia index 0.54 vs 0.07 at CPP 38).  Autoregulation can only ever add
// flow, so the pressure-passive curve must be the floor.
double f_passive = clampd(CPP/CPP0, 0.0, 1.6);
double f_auto    = clampd((CPP - 25.0)/25.0, 0.0, 1.0);
double CBF  = CBF0 * (f_passive + AUTOR*((f_auto > f_passive) ? (f_auto - f_passive) : 0.0));
double Isch = clampd(1.0 - CBF/(CBF_crit*CBF0), 0.0, 1.0);

// ===========================================================================
// Innate immunity
// ===========================================================================
double S_pamp = hillf(CW, K_cw) + 0.6*hillf(NC/1.0e6, K_nl) + 0.5*hillf(PLY, K_ply);
dxdt_MG = kmg_on*S_pamp*(1.0 - MG) - kmg_off*MG;

double f10 = 1.0/(1.0 + IL10/K_il10);
dxdt_TNF   = kTNF*MG*(1.0 - Idex)*f10 - kel_TNF*TNF;
dxdt_IL1   = kIL1*MG*(0.6 + 0.4*hillf(PLY, K_ply))*(1.0 - Idex)*f10 - kel_IL1*IL1;
dxdt_IL6   = kIL6*(0.5*MG + hillf(TNF, K_tnf) + hillf(IL1, K_il1))*(1.0 - 0.6*Idex)
             - kel_IL6*IL6;
dxdt_IL10  = kIL10*(MG + hillf(TNF, K_tnf))*(1.0 + 0.5*Idex) - kel_IL10*IL10;
dxdt_CXCL8 = kC8*(hillf(TNF, K_tnf) + hillf(IL1, K_il1) + 0.3*MG)*(1.0 - 0.7*Idex)
             - kel_C8*CXCL8;
dxdt_COMP  = kcomp*hillf(NC/1.0e6, 2.0) - kel_comp*COMP;

double adh_mol = hillf(TNF, K_tnf)*(1.0 - 0.6*Idex);
dxdt_PMN  = k_influx*hillf(CXCL8, K_c8)*adh_mol*(1.0 + 0.5*(PB - 1.0)/19.0)
            - (k_egress + k_apop)*PMN;
dxdt_MMP9 = kMMP*hillf(PMN, 2500.0)*(1.0 - 0.5*Idex) - kel_MMP*MMP9;
dxdt_ROS  = kROS*(hillf(PMN, 2500.0) + 0.5*MG + 0.3*Isch) - kel_ROS*ROS;

// ===========================================================================
// Barrier - the door that works both ways
// ===========================================================================
// Dexamethasone tightens the barrier both indirectly (less MMP-9, less TNF)
// and directly (tight-junction protein expression), the dex_pb term.  Both
// arms reduce hydrophilic antibiotic entry.
dxdt_PB = k_pb*(hillf(MMP9, K_mmp) + 0.5*hillf(TNF, K_tnf) + 0.3*hillf(PMN, 2500.0))
          *(Pb_max - PB)
          - k_pb_off*(1.0 + dex_pb*Idex)*(PB - 1.0);
dxdt_ALB_CSF  = (PS_alb*PB*(Alb_ser - ALB_CSF) - Qbulk*ALB_CSF)/Vcsf;
dxdt_PROT_CSF = (PS_prot*PB*(Prot_ser/10.0 - PROT_CSF) - Qbulk*PROT_CSF)/Vcsf;

// ===========================================================================
// CSF metabolism.  The GLUT1 carrier is BIDIRECTIONAL, which is why CSF
// glucose almost never reaches zero: as it falls, net influx rises.
// ===========================================================================
double glut = 1.0 - inflam_glut*hillf(PB - 1.0, 8.0);
double influx = Tmax_glc*glut*(Glc_pl/(Km_glut + Glc_pl) - GLC_CSF/(Km_glut + GLC_CSF));
double use_pmn   = q_pmn*PMN*f_glc;
double use_bact  = q_bact*(NC/1.0e7)*f_glc;
double use_brain = q_brain*(1.0 + 0.3*(TEMP - 37.0))*f_glc;
dxdt_GLC_CSF = influx - use_pmn - use_bact - use_brain - (Qbulk/Vcsf)*GLC_CSF;
// 2 lactate per glucose consumed.  Because production is now tied to actual
// glucose consumption, CSF lactate is bounded by the glucose influx rate.
dxdt_LAC_CSF = 2.0*0.0556*(use_pmn + use_bact + use_brain*(0.3 + 0.7*Isch))
               - k_lac*(LAC_CSF - Lac_base);

// ===========================================================================
// CSF hydrodynamics.  Exponential compliance (dP/dV = K_el * P) from the
// Marmarou pressure-volume index; at steady state this reduces to the Davson
// relation ICP = P_ss + Qf * R_out, which is how R_out was calibrated.
// ===========================================================================
double Qf = Qf0*(1.0 - 0.3*clampd((ICP - 20.0)/40.0, 0.0, 1.0));
dxdt_R_OUT = k_ro_on*(hillf(PROT_CSF, 200.0) + hillf(PMN, 2000.0))*(R_out_max - R_OUT)
             - k_ro_off*(R_OUT - R_out0);
double absorb = (ICP > P_ss) ? (ICP - P_ss)/R_OUT : 0.0;

// Osmotherapy efficacy is inversely proportional to barrier leak: the very
// disruption that causes the oedema also lets the osmotic agent equilibrate.
double osm_shrink = k_osm*Osm_grad/(1.0 + 0.3*(PB - 1.0));
double dVbr = k_vas*hillf(PB - 1.0, 6.0)*(1.0 - VBR/Vbr_max)
              + k_cyt*Isch - osm_shrink - k_vbr_res*VBR;
dxdt_VBR      = dVbr;
dxdt_VCSF_NET = Qf - absorb - CSF_DRAIN;
dxdt_ICP      = K_el*ICP*(Qf + dVbr - absorb - CSF_DRAIN);
dxdt_AUTOR    = -k_ar_loss*hillf(TNF, K_tnf)*AUTOR + k_ar_rec*(1.0 - AUTOR);

// ===========================================================================
// Systemic
// ===========================================================================
dxdt_MAP  = -k_map*(hillf(TNF, K_tnf) + 0.4*hillf(IL6, K_il6))*MAP
            + k_map_rec*(MAP0*VOL - MAP)
            + k_cush*((ICP > 30.0) ? (ICP - 30.0) : 0.0) + VASOPRESS;
dxdt_TEMP = k_temp*(hillf(IL1, K_il1) + 0.5*hillf(IL6, K_il6))*(1.0 - ANTIPYR)
            - k_temp_off*(TEMP - Temp0);
dxdt_SOFA = k_sofa*(hillf(log10(1.0 + NB), 2.0) + 0.5*hillf(TNF, K_tnf))
            - k_sofa_rec*SOFA;
dxdt_VOL  = 0.0;

// ===========================================================================
// Injury - two anatomically separate arms with different drivers
// ===========================================================================
dxdt_NCORT = -(k_isch*SQ(Isch) + k_ros_cort*hillf(ROS, K_ros))*NCORT;
dxdt_NDG   = -(k_apo_dg*(hillf(CW, K_cw) + k_ply_dg*hillf(PLY, K_ply2)
                         + k_ros_dg*hillf(ROS, K_ros)))*NDG;
dxdt_HC    = -(k_hc*(w_hc_ply*hillf(PLY, K_ply2) + w_hc_ros*hillf(ROS, K_ros)
                     + w_hc_pmn*hillf(PMN, 2500.0)))*HC;
dxdt_OSS   = k_oss*(1.0 - HC)*hillf(PMN, 2000.0)*(1.0 - OSS);
dxdt_SZ    = k_sz*((1.0 - NCORT) + hillf(ROS, K_ros)
                   + 0.4*((GLC_CSF < 25.0) ? (1.0 - GLC_CSF/25.0) : 0.0))
             *(1.0 - ANTICONV) - k_sz_off*SZ;

// ===========================================================================
// Mortality hazard.
// BUG FOUND (F16): the h_cort*(1-NCORT) term used to live here, inside the
// integral.  A permanent injury multiplied by elapsed time gives an unbounded
// hazard, so recovered patients accrued risk for 14 days and 22 of 26
// scenarios returned >99% mortality - which erased every treatment comparison.
// Structural damage is now scored ONCE, in $TABLE.
// BUG FOUND (F28): the ICP and CPP terms are quadratic.  Linear terms made
// "CPP 45 for 100 h" equivalent to "CPP 5 for 10 h", so mild sustained
// deviation dominated the hazard and no intervention could move it.
// ===========================================================================
double z_icp = ((ICP > 25.0) ? (ICP - 25.0) : 0.0)/10.0;
double z_cpp = ((CPP < 50.0) ? (50.0 - CPP) : 0.0)/10.0;
dxdt_HAZ = h0 + h_icp*SQ(z_icp) + h_cpp*SQ(z_cpp) + h_sofa*SOFA
           + h_isch*SQ(Isch) + h_bact*log10(1.0 + NB)
           + h_ncsf*log10(1.0 + NC);

// ===========================================================================
// Trackers
// ===========================================================================
dxdt_AUC_CEF = CEF_CSF;
dxdt_AUC_VAN = VAN_CSF;
dxdt_AUC_LYS = lysis_flux/1.0e6;
dxdt_AUC_TNF = TNF;
dxdt_T_CEF   = (CEF_CSF > 4.0*MIC_cef) ? 1.0 : 0.0;
dxdt_T_VAN   = (VAN_CSF > 4.0*MIC_van) ? 1.0 : 0.0;
dxdt_AUC_ICP = ICP;
dxdt_I_ICP   = SQ(z_icp);
dxdt_I_CPP   = SQ(z_cpp);
dxdt_I_SOFA  = SOFA;
dxdt_I_ISCH  = SQ(Isch);
dxdt_I_BACT  = log10(1.0 + NB);
dxdt_I_NCSF  = log10(1.0 + NC);

$TABLE
// ---------------------------------------------------------------------------
// Derived, reportable quantities
// ---------------------------------------------------------------------------
double CPPo = MAP - ICP;
double fp   = clampd(CPPo/CPP0, 0.0, 1.6);
double fa   = clampd((CPPo - 25.0)/25.0, 0.0, 1.0);
double CBFo = CBF0*(fp + AUTOR*((fa > fp) ? (fa - fp) : 0.0));
double ISCHo = clampd(1.0 - CBFo/(CBF_crit*CBF0), 0.0, 1.0);

double QALB  = 1000.0*ALB_CSF/Alb_ser;                 // albumin quotient
double GLCR  = GLC_CSF/Glc_pl;                         // CSF/plasma glucose ratio
double LOGNC = log10(1.0 + NC);
double LOGNB = log10(1.0 + NB);

// Culture-detectable sterility uses a REPORTING threshold, not a state floor.
// BUG FOUND (F29): an extinction floor written into the derivative block is a
// discontinuous switch that stalls the integrator once any compartment clears.
double STERILE = (NC < 10.0) ? 1.0 : 0.0;
double ADH_CLEAR = ((NC + NADH) < 10.0) ? 1.0 : 0.0;

// Audiometric mapping with cochlear reserve: the first ~15% of hair cell loss
// costs almost no threshold, then the curve turns sharply.
double L_hc = 1.0 - HC;
double xdb  = clampd((L_hc - hc_reserve)/hc_span, 0.0, 1.0);
double HEAR_DB = 120.0*pow(xdb, hc_exp) + 15.0*OSS;
if (HEAR_DB > 120.0) HEAR_DB = 120.0;

// Acute hazard integral plus the one-shot structural term.
double HAZ_TOT = HAZ + h_cort_final*(1.0 - NCORT);
double P_DEATH = 1.0 - exp(-HAZ_TOT);
double COG_Z   = -(2.5*(1.0 - NDG) + 1.5*(1.0 - NCORT));
double FOCAL   = 1.0 - exp(-3.0*(1.0 - NCORT));

// Penetration ratios, reported so the delivery step stays visible
double Cef_plt = CEF_C/V1_cef;
double PEN_CEF = (Cef_plt > 1e-9) ? CEF_CSF/Cef_plt : 0.0;
double Van_plt = VAN_C/V1_van;
double PEN_VAN = (Van_plt > 1e-9) ? VAN_CSF/Van_plt : 0.0;

// Kill rates and the product, for the diagnostic plots
double KCEF = emax_kill(CEF_CSF, MIC_cef, Emax_cef, EC50r_cef, h_cef);
double KVAN = emax_kill(VAN_CSF, MIC_van, Emax_van, EC50r_van, h_van);
double KRIF = emax_kill(RIF_CSF, MIC_rif, Emax_rif, EC50r_rif, h_rif);
double KTOT = KCEF + KVAN + KRIF;
double LYSFLUX = ((KCEF*Y_cef + KVAN*Y_van + KRIF*Y_rif)*NC
                  + k_autolysis*Y_auto*NC)/1.0e6;
double IDEX = Imax_dex*hillf(DEX_TR, IC50_TR);

$CAPTURE
CPPo CBFo ISCHo QALB GLCR LOGNC LOGNB STERILE ADH_CLEAR HEAR_DB HAZ_TOT
P_DEATH COG_Z FOCAL PEN_CEF PEN_VAN KCEF KVAN KRIF KTOT LYSFLUX IDEX
'

mod <- mcode_cache("abm", abm_code)

## =============================================================================
##  Dosing helpers
## =============================================================================
WT <- 70

ev_cef <- function(t0 = 0, n = 28, dose = 2000, dur = 0.5)
  ev(time = t0, amt = dose, cmt = 1, rate = dose / dur, ii = 12, addl = n - 1)

## Continuous infusion arm: 2 g load, then 4 g/day.  For a susceptible strain
## this changes nothing (T>4xMIC is already 100%); it only matters when the
## C/MIC margin is thin, which is the point of scenario S25.
ev_cef_ci <- function(t0 = 0, days = 14)
  ev(time = t0, amt = 2000, cmt = 1, rate = 4000) +
  ev(time = t0 + 0.5, amt = 4000 * days, cmt = 1, rate = 4000 / 24)

ev_van <- function(t0 = 0, n = 56, mgkg = 15, dur = 1)
  ev(time = t0, amt = mgkg * WT, cmt = 4, rate = mgkg * WT / dur,
     ii = 6, addl = n - 1)

ev_rif <- function(t0 = 0, n = 28, dose = 600, dur = 0.5)
  ev(time = t0, amt = dose, cmt = 7, rate = dose / dur, ii = 12, addl = n - 1)

## Dexamethasone.  t0 is the whole story: negative means before the first
## antibiotic dose, positive means after it.
ev_dex <- function(t0 = -0.33, days = 4, mgkg = 0.15, dur = 0.25)
  ev(time = t0, amt = mgkg * WT, cmt = 9, rate = mgkg * WT / dur,
     ii = 6, addl = 4 * days - 1)

ev_mann <- function(t0 = 1, n = 12, gkg = 0.5, dur = 0.5)
  ev(time = t0, amt = gkg * WT * 1000, cmt = 12,
     rate = gkg * WT * 1000 / dur, ii = 6, addl = n - 1)

ev_gly <- function(t0 = 0, n = 16, gkg = 1.5)
  ev(time = t0, amt = gkg * WT * 1000, cmt = 13, ii = 6, addl = n - 1)

## =============================================================================
##  26 therapeutic scenarios
##  (identical to the Python reference; abm_reference_output.txt holds the
##   numbers these produce)
## =============================================================================
## BUG FOUND (F33): an event at a negative time is never evaluated, because
## integration starts at t = 0.  Every "give X before Y" scenario silently lost
## its first dose.  Shift the whole schedule so the earliest event sits at t = 0
## - which is also the physically correct reading: at presentation you give the
## steroid, and the antibiotic follows 20 minutes later.
zero_shift <- function(e) {
  if (is.null(e)) return(NULL)
  d <- as.data.frame(e)
  tmin <- min(d$time)
  if (tmin < 0) d$time <- d$time - tmin
  as.ev(d)
}

run_abm <- function(dosing = NULL, end = 336, delta = 0.25, ...) {
  p <- list(...)
  m <- mod
  if (length(p)) m <- param(m, p)
  dosing <- zero_shift(dosing)
  if (is.null(dosing)) {
    mrgsim(m, end = end, delta = delta, hmax = 0.05)
  } else {
    mrgsim(m, events = dosing, end = end, delta = delta, hmax = 0.05)
  }
}

scenarios <- list(
  S01 = list(lab = "무치료 자연경과",             ev = NULL,                                     p = list()),
  S02 = list(lab = "CEF 단독 (t=0)",              ev = ev_cef(),                                 p = list()),
  S03 = list(lab = "CEF + DEX 20분 전",            ev = ev_cef() + ev_dex(-0.33),                 p = list()),
  S04 = list(lab = "CEF + DEX 동시 (0 h)",         ev = ev_cef() + ev_dex(0),                     p = list()),
  S05 = list(lab = "CEF + DEX +2 h",               ev = ev_cef() + ev_dex(2),                     p = list()),
  S06 = list(lab = "CEF + DEX +4 h",               ev = ev_cef() + ev_dex(4),                     p = list()),
  S07 = list(lab = "CEF + DEX +12 h",              ev = ev_cef() + ev_dex(12),                    p = list()),
  S08 = list(lab = "CEF + DEX 1일만",              ev = ev_cef() + ev_dex(0, days = 1),           p = list()),
  S09 = list(lab = "항생제 3 h 지연",               ev = ev_cef(3) + ev_dex(3),                    p = list()),
  S10 = list(lab = "항생제 6 h 지연",               ev = ev_cef(6) + ev_dex(6),                    p = list()),
  S11 = list(lab = "항생제 12 h 지연",              ev = ev_cef(12) + ev_dex(12),                  p = list()),
  S12 = list(lab = "내성균(MIC 4) CEF 단독",        ev = ev_cef(),                                 p = list(MIC_cef = 4)),
  S13 = list(lab = "내성균 + VAN, DEX 없음",        ev = ev_cef() + ev_van(),                      p = list(MIC_cef = 4)),
  S14 = list(lab = "내성균 + VAN + DEX",           ev = ev_cef() + ev_van() + ev_dex(0),          p = list(MIC_cef = 4)),
  S15 = list(lab = "내성균 + VAN + RIF + DEX",     ev = ev_cef() + ev_van() + ev_rif() + ev_dex(0), p = list(MIC_cef = 4)),
  S16 = list(lab = "RIF 2 h 선행 → CEF",           ev = ev_cef() + ev_rif(-2),                    p = list()),
  S17 = list(lab = "CEF + DEX + 만니톨",           ev = ev_cef() + ev_dex(0) + ev_mann(),         p = list()),
  S18 = list(lab = "CEF + DEX + 글리세롤",          ev = ev_cef() + ev_dex(0) + ev_gly(),          p = list()),
  S19 = list(lab = "고균량 1e8 + DEX",             ev = ev_cef() + ev_dex(0),                     p = list(N0 = 1e8)),
  S20 = list(lab = "조기 저균량 1e5 + DEX",         ev = ev_cef() + ev_dex(0),                     p = list(N0 = 1e5)),
  S21 = list(lab = "고령/면역저하 + DEX",           ev = ev_cef() + ev_dex(0),                     p = list(HOST_DEF = 0.4, MU_SCALE = 1.3)),
  S22 = list(lab = "CEF 지속주입 + DEX",           ev = ev_cef_ci() + ev_dex(0),                  p = list()),
  S23 = list(lab = "해열제만 (항생제 없음)",         ev = NULL,                                     p = list(ANTIPYR = 0.7)),
  S24 = list(lab = "표준치료 CEF+VAN+DEX",         ev = ev_cef() + ev_van() + ev_dex(-0.33),      p = list()),
  S25 = list(lab = "내성균 CEF지속주입+VAN+DEX",    ev = ev_cef_ci() + ev_van() + ev_dex(-0.33),   p = list(MIC_cef = 4)),
  S26 = list(lab = "표준치료 + EVD 배액",           ev = ev_cef() + ev_van() + ev_dex(-0.33),      p = list(CSF_DRAIN = 12))
)

simulate_scenario <- function(s, end = 336) {
  args <- c(list(dosing = s$ev, end = end), s$p)
  do.call(run_abm, args)
}

## -----------------------------------------------------------------------------
##  Endpoint extraction
## -----------------------------------------------------------------------------
endpoints_abm <- function(sim) {
  d <- as.data.frame(sim)
  last <- d[nrow(d), ]
  st <- d$time[d$STERILE > 0.5]
  ad <- d$time[d$ADH_CLEAR > 0.5]
  data.frame(
    sterile_h   = if (length(st)) min(st) else NA_real_,
    adh_clear_h = if (length(ad)) min(ad) else NA_real_,
    peak_ICP    = max(d$ICP),
    min_CPP     = min(d$CPPo),
    min_MAP     = min(d$MAP),
    peak_TNF    = max(d$TNF),
    peak_PMN    = max(d$PMN),
    peak_Pb     = max(d$PB),
    peak_Qalb   = max(d$QALB),
    min_Glc     = min(d$GLC_CSF),
    peak_Lac    = max(d$LAC_CSF),
    peak_Prot   = max(d$PROT_CSF),
    AUC_lysis   = last$AUC_LYS,
    AUC_TNF     = last$AUC_TNF,
    AUC_van     = last$AUC_VAN,
    T_cef_4MIC  = last$T_CEF,
    T_van_4MIC  = last$T_VAN,
    hear_dB     = last$HEAR_DB,
    cog_z       = last$COG_Z,
    haz_acute   = last$HAZ,
    p_death     = last$P_DEATH
  )
}

run_all_scenarios <- function(end = 336) {
  out <- lapply(names(scenarios), function(k) {
    e <- endpoints_abm(simulate_scenario(scenarios[[k]], end = end))
    cbind(id = k, label = scenarios[[k]]$lab, e)
  })
  do.call(rbind, out)
}

## =============================================================================
##  Virtual cohort: 10 patients x (with / without dexamethasone)
##  Compared against de Gans & van de Beek NEJM 2002 pneumococcal subgroup
##  (death 34% -> 14%, unfavourable outcome 52% -> 26%).
## =============================================================================
cohort <- data.frame(
  id       = 1:10,
  N0       = c(1e6, 3e6, 1e7, 3e7, 1e8, 1e7, 5e6, 2e7, 8e7, 3e7),
  delay_h  = c(0, 2, 1, 4, 6, 3, 1, 8, 2, 24),
  host_def = c(1.00, 1.00, 0.90, 0.85, 0.70, 0.55, 1.00, 0.60, 0.80, 0.40)
)

run_cohort <- function(end = 336) {
  res <- lapply(seq_len(nrow(cohort)), function(i) {
    r <- cohort[i, ]
    mu <- 1 + 0.5 * (1 - r$host_def)
    base <- ev_cef(r$delay_h) + ev_van(r$delay_h)
    a <- endpoints_abm(run_abm(base, end = end, N0 = r$N0,
                               HOST_DEF = r$host_def, MU_SCALE = mu))
    b <- endpoints_abm(run_abm(base + ev_dex(r$delay_h - 0.33), end = end,
                               N0 = r$N0, HOST_DEF = r$host_def, MU_SCALE = mu))
    data.frame(id = r$id, logN0 = log10(r$N0), delay_h = r$delay_h,
               host_def = r$host_def,
               death_no_dex = a$p_death, death_dex = b$p_death,
               dB_no_dex = a$hear_dB, dB_dex = b$hear_dB)
  })
  do.call(rbind, res)
}

## =============================================================================
##  Focused analyses that make the model's three claims visible
## =============================================================================

## (1) The product peaks at the first dose.
##     Plot LYSFLUX and TNF over the first 24 h of S02.
claim1_burst <- function() {
  d <- as.data.frame(run_abm(ev_cef(), end = 24, delta = 0.05))
  d[, c("time", "CEF_CSF", "KCEF", "LOGNC", "LYSFLUX", "CW", "PLY", "MG", "TNF", "PMN")]
}

## (2) Dexamethasone timing sweep.  The benefit is a step function of whether
##     the shield was up when the burst happened, not a smooth dose-response.
claim2_timing <- function(times = c(-0.33, 0, 1, 2, 4, 8, 12, 24)) {
  do.call(rbind, lapply(times, function(tt) {
    e <- endpoints_abm(run_abm(ev_cef() + ev_dex(tt)))
    data.frame(dex_time_h = tt, AUC_TNF = e$AUC_TNF, peak_PMN = e$peak_PMN,
               peak_Pb = e$peak_Pb, hear_dB = e$hear_dB, p_death = e$p_death)
  }))
}

## (3) The sign flip.  Same steroid, opposite consequence, depending on how
##     much C/MIC margin the antibiotic has.
claim3_signflip <- function() {
  rbind(
    cbind(arm = "susceptible, CEF, no DEX",
          endpoints_abm(run_abm(ev_cef()))),
    cbind(arm = "susceptible, CEF + DEX",
          endpoints_abm(run_abm(ev_cef() + ev_dex(-0.33)))),
    cbind(arm = "resistant, CEF+VAN, no DEX",
          endpoints_abm(run_abm(ev_cef() + ev_van(), MIC_cef = 4))),
    cbind(arm = "resistant, CEF+VAN + DEX",
          endpoints_abm(run_abm(ev_cef() + ev_van() + ev_dex(-0.33), MIC_cef = 4)))
  )
}

## (4) Antibiotic delay sweep: the injury integral grows because N grows.
claim4_delay <- function(delays = c(0, 1, 3, 6, 12, 24)) {
  do.call(rbind, lapply(delays, function(dl) {
    e <- endpoints_abm(run_abm(ev_cef(dl) + ev_van(dl) + ev_dex(dl - 0.33)))
    data.frame(delay_h = dl, AUC_lysis = e$AUC_lysis, peak_Lac = e$peak_Lac,
               sterile_h = e$sterile_h, hear_dB = e$hear_dB, p_death = e$p_death)
  }))
}

## (5) Hazard decomposition: which physiological derangement actually kills.
hazard_decomposition <- function(ids = c("S01", "S02", "S03", "S11", "S14",
                                         "S20", "S24", "S26")) {
  p <- as.list(param(mod))
  do.call(rbind, lapply(ids, function(k) {
    d <- as.data.frame(simulate_scenario(scenarios[[k]]))
    L <- d[nrow(d), ]
    data.frame(
      id = k, label = scenarios[[k]]$lab,
      icp    = p$h_icp  * L$I_ICP,
      cpp    = p$h_cpp  * L$I_CPP,
      sofa   = p$h_sofa * L$I_SOFA,
      isch   = p$h_isch * L$I_ISCH,
      bact   = p$h_bact * L$I_BACT,
      ncsf   = p$h_ncsf * L$I_NCSF,
      base   = p$h0 * max(d$time),
      struct = p$h_cort_final * (1 - L$NCORT),
      total  = L$HAZ_TOT, p_death = L$P_DEATH
    )
  }))
}

## =============================================================================
##  Structural checks that must pass before any of the above means anything.
##  These are hand-verifiable identities, not fitted results.
## =============================================================================
structural_checks <- function() {
  p <- as.list(param(mod))
  cat("CSF steady-state ICP = P_ss + Qf * R_out\n")
  for (r in c(p$R_out0, 0.48, p$R_out_max))
    cat(sprintf("   R_out = %5.3f -> ICP = %5.1f mmHg\n",
                r, p$P_ss + p$Qf0 * r))
  cat("   [literature: normal 7-15, meningitis 20-40]\n\n")

  fs <- function(Ct) { b <- p$Kd_cef + p$Bmax_cef - Ct
                       (-b + sqrt(b^2 + 4 * p$Kd_cef * Ct)) / 2 }
  cat(sprintf("Ceftriaxone saturable binding: fu(30) = %.3f, fu(250) = %.3f\n",
              fs(30) / 30, fs(250) / 250))
  cat("   [reported: free fraction rises 0.05 -> 0.20 with dose]\n\n")

  pen <- function(PS, Eff, pb, lip = FALSE) {
    pbe <- if (lip) 1 + p$a_rif * (pb - 1) else pb
    PS * pbe / (PS * pbe + p$Qf0 + Eff)
  }
  cat("CSF/free-plasma steady-state ratio\n")
  cat(sprintf("   CEF Pb=1 %.3f  Pb=9 %.3f  (x fu -> total %.3f / %.3f)\n",
              pen(p$PS_cef, p$Eff_cef, 1), pen(p$PS_cef, p$Eff_cef, 9),
              pen(p$PS_cef, p$Eff_cef, 1) * fs(250) / 250,
              pen(p$PS_cef, p$Eff_cef, 9) * fs(250) / 250))
  cat(sprintf("   VAN Pb=1 %.3f  Pb=9 %.3f  (x fu -> total %.3f / %.3f)\n",
              pen(p$PS_van, p$Eff_van, 1), pen(p$PS_van, p$Eff_van, 9),
              pen(p$PS_van, p$Eff_van, 1) * p$fu_van,
              pen(p$PS_van, p$Eff_van, 9) * p$fu_van))
  cat(sprintf("   RIF Pb=1 %.3f  Pb=9 %.3f  (lipophilic, Pb-insensitive)\n",
              pen(p$PS_rif, p$Eff_rif, 1, TRUE), pen(p$PS_rif, p$Eff_rif, 9, TRUE)))
  cat("   [reported total penetration: CEF 1.5-10%, VAN 1-15%, RIF 10-20%]\n\n")

  g <- 60
  cat(sprintf("CSF glucose balance at 60 mg/dL: net influx %.1f = bulk %.1f + brain %.1f\n",
              p$Tmax_glc * (p$Glc_pl / (p$Km_glut + p$Glc_pl) - g / (p$Km_glut + g)),
              g * p$Qf0 / p$Vcsf, p$q_brain))
  cat(sprintf("Glucose competition: PMN at 4000/uL uses %.1f, bacteria at 1e7 use %.1f mg/dL/h\n",
              p$q_pmn * 4000, p$q_bact))
  cat("   -> neutrophils outconsume the bacteria about tenfold, which is why\n")
  cat("      hypoglycorrhachia persists after the CSF is culture-negative\n\n")

  cat(sprintf("RK4/LSODA stability: ICP eigenvalue K_el*ICP/R_out = %.1f /h at ICP 60,\n",
              p$K_el * 60 / p$R_out0))
  cat(sprintf("   recovery-phase R_out -> explicit step must stay below %.3f h\n",
              2.785 / (p$K_el * 60 / p$R_out0)))
}

## =============================================================================
##  Usage
## =============================================================================
if (FALSE) {
  structural_checks()

  ## the three claims
  head(claim1_burst(), 40)
  claim2_timing()
  claim3_signflip()
  claim4_delay()
  hazard_decomposition()

  ## everything
  all_sc <- run_all_scenarios()
  print(all_sc, digits = 3)

  co <- run_cohort()
  cat(sprintf("cohort mean mortality: no DEX %.1f%%, DEX %.1f%%\n",
              100 * mean(co$death_no_dex), 100 * mean(co$death_dex)))
  cat(sprintf("cohort mean hearing threshold shift: no DEX %.1f dB, DEX %.1f dB\n",
              mean(co$dB_no_dex), mean(co$dB_dex)))

  ## a look at the burst itself
  s <- run_abm(ev_cef() + ev_dex(-0.33), end = 48, delta = 0.05)
  plot(s, LYSFLUX + CW + PLY + TNF + PMN + IDEX ~ time)

  ## and at the door
  s2 <- run_abm(ev_cef() + ev_van() + ev_dex(-0.33), end = 168)
  plot(s2, PB + PEN_VAN + VAN_CSF + LOGNC ~ time)
}
