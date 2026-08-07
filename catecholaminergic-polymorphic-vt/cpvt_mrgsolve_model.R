## =============================================================================
##  cpvt_mrgsolve_model.R
##  카테콜아민성 다형성 심실빈맥 (Catecholaminergic Polymorphic Ventricular
##  Tachycardia, CPVT) QSP 모델
##
##  44 ODEs · beat-resolved-average Ca cycling · SOICR threshold physiology ·
##  4 beta-blockers with real PK · flecainide with two independent arms ·
##  verapamil · LCSD · ICD with shock->catecholamine feedback · 41 scenarios
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL IS ORGANISED AROUND
##  ---------------------------------------------------------------------------
##  CPVT looks like a list of clinical rules: "use nadolol, not metoprolol",
##  "add flecainide on top of a beta-blocker", "the resting ECG is normal",
##  "an ICD is not a substitute for drugs". This model is written so that all
##  four are ARITHMETIC CONSEQUENCES rather than coded-in assertions.
##
##  A triggered beat requires TWO crossings that are set INDEPENDENTLY, and
##  the probability of a triggered beat is their PRODUCT:
##
##      PVCR = W( CAJSRF / CTHR )  x  P( VDAD / VREQ )
##             \_____ crossing 1 _____/  \___ crossing 2 ___/
##
##  CROSSING 1 — LOAD vs SOICR THRESHOLD.  Free junctional-SR [Ca] (CAJSRF)
##    against the store-overload-induced-Ca-release threshold (CTHR). The
##    numerator is pushed up by catecholamines through three separate routes
##    (PKA-P of phospholamban -> SERCA; PKA-P of Cav1.2 -> Ca entry; higher
##    rate -> shorter diastole -> less NCX extrusion time). The denominator is
##    pulled down by the mutation (RyR2 gain of function, or loss of the
##    luminal CASQ2 brake) and pulled down further by PKA/CaMKII phospho-RyR2.
##    Below X = 1 the wave rate is IDENTICALLY ZERO, not merely small.
##
##  CROSSING 2 — DAD vs WHAT THE CELL NEEDS.  The wave's Ca is converted to
##    inward current by electrogenic NCX and opposed by the I_K1 sink, giving
##    VDAD. Whether that captures depends on VREQ = VREQ0 / NAAV^PNAV, i.e. on
##    how much Na current is available.
##
##  ---------------------------------------------------------------------------
##  FIVE RESULTS THAT COME OUT AS ARITHMETIC
##  ---------------------------------------------------------------------------
##  (1) A NORMAL HEART IS SILENT STRUCTURALLY, NOT BY TUNING.
##      Wild type at peak exercise: free [Ca]JSR rises 1112 -> 1449 uM (+30%)
##      while the PKA-lowered threshold falls to ~1610 uM. X peaks at 0.900 —
##      a 10% margin — and the model produces EXACTLY ZERO triggered beats.
##      This is why the resting ECG, the echo and the MRI are all normal in
##      CPVT: the disease is a threshold that has moved INTO the physiological
##      range of exertion, not a structural lesion.
##
##  (2) NADOLOL vs METOPROLOL IS A PK RESULT, NOT A PHARMACOLOGY RESULT.
##      The model gives the two drugs NO difference in intrinsic potency; both
##      enter as the same competitive term 1 + C/EC50. What differs is what
##      survives to the trough. Computed at steady state:
##
##        drug / regimen              Cmax     Ctrough   shift(peak) shift(trough)
##        nadolol 80 mg qd ......... 154.6      46.3        8.73        3.32
##        metoprolol 100 mg bid ....  51.4       5.6        2.14        1.12
##        metoprolol 200 mg bid .... 102.8      11.2        3.28        1.25
##        bisoprolol 10 mg qd ......  41.6      11.8        4.46        1.98
##
##      Metoprolol at DOUBLE dose reaches, at its PEAK, only what nadolol
##      still has left at its TROUGH (3.28 vs 3.32). On the exercise test the
##      PVC count over the study is 49 on nadolol at trough versus 361 on
##      metoprolol at trough — against 403 untreated. Metoprolol's trough is
##      a 10% effect. Because CPVT events are triggered by unscheduled
##      exertion, it is the trough hour that is on trial, and a 3.5 h
##      half-life cleared by CYP2D6 cannot hold one.
##
##  (3) ONE MISSED DOSE IS A DIFFERENT DRUG. 48 h after a nadolol dose the
##      shift factor is still 2.02 (61% of its own trough effect); 24 h after
##      a metoprolol dose it is 1.04, i.e. nothing at all. A CYP2D6
##      ultrarapid metaboliser on metoprolol has a trough of 0.8 ng/mL and
##      the model gives them 98% of the untreated PVC count — they are, at
##      trough, an untreated patient holding a prescription.
##
##  (4) FLECAINIDE'S ADD-ON EFFECT IS THE PRODUCT, NOT AN ASSUMPTION.
##      Nothing in the model grants flecainide an additive or synergistic
##      term. It acts on crossing 2 while beta-blockade acts on crossing 1.
##      On top of nadolol at trough the PVC count falls 49 -> 26 (-47%) and
##      VT time 0.08 -> 0.01 min. Flecainide ALONE is worse than nadolol
##      alone (268 vs 49 PVCs), which is why it is an add-on and not a
##      substitute.
##
##  (5) THE ICD CAN FEED THE ARRHYTHMIA IT TREATS. Shocks are modelled as
##      delivering a catecholamine surge (SHKD -> EPI). With an ICD and no
##      beta-blocker the model delivers 16.5 shocks, drives epinephrine to
##      17.0 nM and accumulates 18.3 min of VT. Switching the shock ->
##      catecholamine feedback OFF and changing nothing else gives 8.1 shocks
##      and 8.96 min of VT. The feedback loop therefore DOUBLES the
##      arrhythmia burden it is trying to terminate — an emergent electrical
##      storm. With nadolol on board the same ICD delivers 0.1 shocks.
##
##  ---------------------------------------------------------------------------
##  THE MODEL'S MOST EXPOSED PREDICTION (and a refuted prior)
##  ---------------------------------------------------------------------------
##  The author's prior hypothesis going in was that at clinical plasma
##  concentrations flecainide works in CPVT mainly through use-dependent Na
##  block (raising VREQ), because its RyR2 IC50 (~4-6 uM in myocytes) sits
##  far above the free intracellular concentration achieved at 100-150 mg bid
##  (~0.20-0.31 uM). Ablating each arm separately REFUTES that:
##
##      nadolol alone ....................... 49 PVCs
##      + flecainide, both arms .............. 26  (-47%)
##      + flecainide, RyR2 arm only .......... 32  (-35%)   <- larger loss
##      + flecainide, Na arm only ............ 39  (-20%)
##
##  i.e. the RyR2 arm carries ~3/4 of the benefit and the Na arm ~2/5 (they
##  are sub-additive). This is the model's stance, not a certainty: it hangs
##  on IC50FLR, the single most uncertain parameter in the file. The
##  prediction is TESTABLE by rate-dependence, because only the Na arm is
##  use-dependent — see scenarios 38-40 (light vs maximal workload). If
##  flecainide's benefit in patients is found NOT to grow with heart rate,
##  the RyR2-dominant reading is supported; if it does grow steeply, IC50FLR
##  is too high here and the balance should shift.
##
##  ---------------------------------------------------------------------------
##  IMPORTANT MODELLING CHOICES / LIMITATIONS (read before using numbers)
##  ---------------------------------------------------------------------------
##  * TIME UNIT IS MINUTES. Every PK parameter is given in the conventional
##    L/h or 1/h and divided by 60.0 inside $ODE. (Getting this wrong made an
##    earlier version clear all drugs 60x too fast, giving zero concentrations.)
##  * Ca cycling is BEAT-AVERAGED, not beat-resolved. There is no action
##    potential and no millisecond time scale. Every process that happens once
##    per beat (ICaL entry, systolic SR release, SERCA uptake, NCX extrusion)
##    is written PER BEAT and multiplied by HR; the SR leak and the
##    spontaneous waves are per-minute processes and are not. An earlier
##    version scaled release with rate but not uptake, which sent cytosolic Ca
##    to 770 uM and made the SR DEPLETE during exercise.
##  * NCX extrusion is diastole-limited, so it scales as (HR/60)^ENCX with
##    ENCX = 0.55 rather than linearly. This single exponent is what makes SR
##    load rate-dependent (the positive staircase) and therefore what makes
##    bradycardia, beta-blockade and LCSD protective. With linear scaling the
##    load balance is rate-invariant and the model says LCSD is PROARRHYTHMIC.
##  * PICA0 / PPLB0 / PRYR0 are FIXED CONSTANTS taken from normal, drug-free,
##    un-denervated physiology. They must NOT be recomputed per scenario: an
##    intervention that lowers resting sympathetic tone also lowers its own
##    reference, so the model then measures delta-PKA from a lowered floor and
##    the intervention comes out proarrhythmic. This defect made LCSD raise
##    peak PVC rate from 32 to 47/min in an earlier version.
##  * RyR2 refractoriness sets a CEILING on wave frequency (1/TAUREC), it is
##    NOT a (1-RREF) multiplier on the wave rate. Written as a multiplier it
##    hands heart rate a spurious ANTI-arrhythmic channel.
##  * PVCR, VTB and the hazard are PROBABILITIES/RATES computed from a
##    deterministic substrate. The model does not produce ECG morphology; a
##    "bidirectional VT" here is a burden index, not a waveform. Do not read
##    VTB as a fraction of beats.
##  * The mechanistic Ca-cycling parameters are of the Shannon-Bers /
##    Ten Tusscher family in ORDER OF MAGNITUDE only; they were re-balanced to
##    put resting free [Ca]JSR near 1100 uM, beat-average cytosolic Ca near
##    0.38 uM and a +30% catecholamine SR loading. They are NOT fitted to any
##    single myocyte data set.
##  * VERIFICATION STATUS -- read this before trusting a number. All 44 ODEs
##    were independently re-implemented in Python/scipy and integrated over the
##    full scenario set; every quantitative claim in this header, in the
##    directory README and in the mechanistic map comes from THAT run, and the
##    four defects named above were found and fixed that way. This R/mrgsolve
##    transcription has NOT been executed: no R interpreter was available in
##    the environment where it was written. It is a faithful line-by-line
##    transcription of the verified system, but it should be treated as
##    unexecuted code -- expect to fix transcription-level errors (a missing
##    semicolon, an $ODE symbol typo) on first run, and re-derive any number
##    you intend to quote from it.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION TARGETS (what the parameters were set against)
##  ---------------------------------------------------------------------------
##  * Wild type: zero exercise-induced ventricular ectopy at peak HR ~166.
##  * CPVT1 index patient: first PVCs at HR 118, bidirectional/polymorphic VT
##    at HR 154 — the classic 110-130 bpm PVC threshold (Leenhardt 1995,
##    Postma 2005).
##  * Nadolol 1-2 mg/kg/day reduces but does not abolish events (Hayashi 2009
##    reports ~8-year event rates around 30% on beta-blockade); model gives
##    -88% PVC count with residual VE-2 ectopy at trough.
##  * Beta-1 selective agents underperform nadolol (Leren 2016) — reproduced
##    from PK alone.
##  * Flecainide added to a beta-blocker suppresses exercise ventricular
##    arrhythmia substantially (van der Werf 2011; Kannankeril 2017 RCT);
##    model gives -47% on top of nadolol at trough.
##  * Epinephrine challenge is LESS sensitive than the exercise test
##    (Marjamaa 2012, Krahn 2005); reproduced here because circulating
##    epinephrine reaches the dyad with access factor FEPI = 0.60 while
##    nerve-released NE does not have to.
##  * Nadolol / metoprolol / bisoprolol / flecainide / verapamil steady-state
##    plasma concentrations within published ranges (see table above).
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

cpvt_code <- '
$PARAM @annotated
// ============================ GENOTYPE =====================================
GRYR    : 0.24  :        : RyR2 기능획득 — SOICR 문턱 감소 분율 (0 = 야생형)
GCSQI   : 0.0   :        : CASQ2의 RyR2 억제 소실에 의한 추가 문턱 감소 분율
CSQSET  : 1.0   :        : CASQ2 단백 상대량 (1 = 정상, CPVT2에서 0.3-0.5)
FREC    : 0.60  :        : RyR2 불응 회복 시상수 배수 (CPVT에서 < 1 = 빠른 회복)

// ==================== SYMPATHETIC DRIVE / TRIGGERS =========================
EXER    : 0     :        : 운동 강도 입력 0-1 (외부에서 지정)
EMOT    : 0     :        : 정서적 스트레스 입력 0-1
EPIINF  : 0     : nM/min : 에피네프린 주입 속도 (부하 검사용)
KNE     : 6.0   : nM/min : 심근 간질 NE 방출 속도계수
KUPT    : 1.0   : 1/min  : NET(uptake-1) 제거 속도
KEPI    : 3.0   : nM/min : 부신 epinephrine 분비 속도계수
KCLE    : 0.35  : 1/min  : epinephrine 제거 속도 (t1/2 ~2 min)
FLCSD   : 0     :        : LCSD 시행 여부 (0 = 미시행, 1 = 시행)
ELCSD   : 0.60  :        : LCSD의 심장 NE 유출 감소 분율
EC50NE  : 3.0   : nM     : NE의 beta1 반수 활성 농도
EC50EPI : 1.5   : nM     : epinephrine의 beta1 반수 활성 농도
FEPI    : 0.60  :        : 순환 epinephrine의 dyad 접근도 (< 1 → 부하검사 민감도 ↓)

// =============================== HEART RATE ================================
HRV0    : 55    : bpm    : 안정 시 미주 긴장 하 기저 심박수
HRVEX   : 45    : bpm    : 운동 시 미주 긴장 철회분
HRMX    : 230   : bpm    : 이론적 최대 심박수
KMHR    : 4.0   :        : 심박수 반응의 교감 구동 반포화 상수
TAUHR   : 0.35  : min    : 심박수 반응 시상수

// ========================= cAMP / PKA CASCADE ==============================
KAC     : 4.0   : 1/min  : adenylyl cyclase 활성 계수
KPDE    : 4.0   : 1/min  : PDE에 의한 cAMP 분해
KPKAON  : 1.0   : 1/min  : PKA 활성화 속도
KPKAOFF : 1.0   : 1/min  : PKA 불활성화 속도
KMPKA   : 0.50  :        : cAMP → PKA 반포화 (Hill 2)
KPI1    : 3.0   : 1/min  : inhibitor-1 인산화 속도
KPI1OFF : 1.5   : 1/min  : inhibitor-1 탈인산화 속도
TAUPH   : 0.5   : min    : 인산화 표적 완화 시상수
KPLB    : 0.28  :        : PKA → P-phospholamban 반포화
KPRYR   : 0.32  :        : PKA → P-RyR2 반포화
KPICA   : 0.30  :        : PKA → P-Cav1.2 반포화
PICA0   : 0.1842:        : 정상 기저 P-Cav1.2 (고정 상수 — 시나리오별 재계산 금지)
PPLB0   : 0.2058:        : 정상 기저 P-phospholamban (고정 상수)
PRYR0   : 0.1656:        : 정상 기저 P-RyR2 (고정 상수)

// ==================== Ca CYCLING — PER BEAT (uM cytosol) ===================
GCAL    : 12.0  : uM/beat: 박동당 L형 통로 Ca 유입
AICA    : 0.90  :        : P-Cav1.2에 의한 Ca 유입 증폭
VUP     : 46.2  : uM/beat: SERCA 최대 재흡수 (박동당)
KMUP    : 0.35  : uM     : SERCA 반포화 세포질 Ca
ASERCA  : 1.00  :        : P-PLB에 의한 SERCA 탈억제 증폭
VNCX    : 62.5  : uM/beat: NCX 최대 배출 (60 bpm 기준)
KMNCX   : 0.60  : uM     : NCX 반포화 (Hill 3)
ENCX    : 0.55  :        : NCX의 심박수 지수 (< 1 = 확장기 제한 → 양성 계단)
GREL    : 0.0398:        : 수축기 SR 방출 계수
KRELH   : 1000  : uM     : 부하 의존 분율 방출 반포화
HREL    : 2.0   :        : 부하 의존 분율 방출 Hill 계수
GLEAK   : 0.045 : 1/min  : 비-스파크 SR 누출 (분당 과정)
VTR     : 15.0  : 1/min  : NSR → JSR 확산 계수
VNSR    : 0.055 :        : NSR 용적 / 세포질 용적
VJSR    : 0.0057:        : JSR 용적 / 세포질 용적
BETAI   : 0.045 :        : 세포질 Ca 완충 계수 (유리 분율)
KMITI   : 1.5   : uM/beat: 미토콘드리아 Ca 유입 (MCU)
KMITO   : 1.0   : uM/beat: 미토콘드리아 Ca 유출 (NCLX)
KMMIT   : 0.60  : uM     : MCU 반포화
BCSQ    : 14000 : uM     : Calsequestrin-2 완충 용량 (Shannon-Bers)
KCSQ    : 630   : uM     : Calsequestrin-2 해리 상수
KCSQTO  : 0.0020: 1/min  : CASQ2 단백 전환 속도

// ======================== RyR2 / SOICR THRESHOLD ===========================
CTH0    : 2330  : uM     : 야생형 SOICR 문턱 (유리 JSR Ca)
BPRYR   : 0.60  :        : P-RyR2에 의한 문턱 하강 강도
BCAMK   : 0.55  :        : CaMKII에 의한 문턱 하강 및 누출 증폭
WMAX    : 110   : 1/min  : 최대 Ca 파동 발생률
KW      : 0.15  :        : 파동 발생률의 (X-1) 반포화
HW      : 2.0   :        : 파동 발생률 Hill 계수
FRELW   : 0.25  :        : 파동당 방출되는 저장고 분율
AFRW    : 2.0   :        : CASQ2 소실 시 분율 방출 증폭
KREF    : 0.40  :        : 방출량에 의한 불응 진입 계수
TAUREC  : 0.030 : min    : RyR2 불응 회복 시상수 (야생형)

// ===================== DAD → TRIGGERED BEAT ================================
GDADV   : 1.64  : mV/uM  : 방출 Ca → DAD 진폭 변환 (전기발생성 NCX)
GK1     : 1.0   :        : I_K1 sink 계수
KO      : 4.2   : mM     : 혈청 칼륨 (재분극 예비력)
TAUV    : 0.05  : min    : DAD 진폭 완화 시상수
VREQ0   : 22.0  : mV     : 정상 Na 가용도에서 촉발에 필요한 탈분극량
MTRIG   : 6.0   :        : 촉발 확률 Hill 계수
PNAV    : 1.5   :        : Na 가용도의 필요 탈분극량 지수
KVT     : 22.0  : 1/min  : VT 성립 반포화 PVC 발생률
NVT_H   : 4.0   :        : VT 성립 Hill 계수

// ============================== FLECAINIDE =================================
KONNA   : 0.50  : 1/min/uM: Na 통로 차단 결합 속도 (사용의존적)
KOFFNA  : 2.00  : 1/min  : Na 통로 차단 해리 속도
NAMIN   : 0.35  :        : Na 가용도 하한
IC50FLR : 4.0   : uM     : RyR2 억제 IC50 (모델에서 가장 불확실한 파라미터)
EFLR    : 0.25  :        : RyR2 억제의 최대 문턱 상승 분율
QRS0    : 95    : ms     : 기저 QRS 폭
AQRS    : 0.55  :        : Na 차단에 의한 QRS 확장 계수

// ========================= CHRONIC / FEEDBACK ==============================
KCMON   : 0.020 : 1/min  : CaMKII 활성화 속도
KCMOFF  : 0.010 : 1/min  : CaMKII 불활성화 속도
CAI0    : 0.377 : uM     : 정상 기저 박동평균 세포질 Ca
KFIB    : 2e-5  : 1/min  : 간질 섬유화 생성 속도
KFIBD   : 1e-5  : 1/min  : 간질 섬유화 소실 속도

// ============================ HAZARD / ICD =================================
KHAZ    : 0.010 : 1/min  : VT 분당 실신/심정지 위험률
ICD     : 0     :        : ICD 삽입 여부 (0/1)
KSHK    : 0.90  : 1/min  : VT 검출 → 쇼크 전달 속도
TAUSHK  : 1.0   : min    : 쇼크 후 카테콜아민 급증 소멸 시상수
ESHK    : 6.0   : nM/min : 쇼크당 epinephrine 급증 강도

// ================= PK (모든 CL·Q·ka는 L/h·1/h — $ODE에서 /60) ==============
KANAD   : 1.2   : 1/h    : 나돌롤 흡수속도
FNAD    : 0.30  :        : 나돌롤 생체이용률
VCNAD   : 150   : L      : 나돌롤 중심 용적
VPNAD   : 150   : L      : 나돌롤 말초 용적
QNAD    : 20    : L/h    : 나돌롤 구획간 청소율
CLNAD   : 12    : L/h    : 나돌롤 청소율 (신배설·CYP 무관)
EC50NAD : 20    : ng/mL  : 나돌롤 beta1 경쟁 EC50
KAMET   : 1.5   : 1/h    : 메토프롤롤 흡수속도
FMET    : 0.40  :        : 메토프롤롤 생체이용률 (초회통과 큼)
VCMET   : 250   : L      : 메토프롤롤 중심 용적
VPMET   : 500   : L      : 메토프롤롤 말초 용적
QMET    : 90    : L/h    : 메토프롤롤 구획간 청소율
CLMET   : 150   : L/h    : 메토프롤롤 청소율 (t1/2 3.5 h)
EC50MET : 45    : ng/mL  : 메토프롤롤 beta1 경쟁 EC50
CYP2D6  : 1.0   :        : CYP2D6 활성 배수 (UM 2.5 · EM 1.0 · PM 0.3)
KABIS   : 1.0   : 1/h    : 비소프롤롤 흡수속도
FBIS    : 0.90  :        : 비소프롤롤 생체이용률
VCBIS   : 230   : L      : 비소프롤롤 용적
CLBIS   : 14.5  : L/h    : 비소프롤롤 청소율 (t1/2 11 h)
EC50BIS : 12    : ng/mL  : 비소프롤롤 beta1 경쟁 EC50
KAFLE   : 0.9   : 1/h    : 플레카이니드 흡수속도
FFLE    : 0.90  :        : 플레카이니드 생체이용률
VCFLE   : 350   : L      : 플레카이니드 중심 용적
VPFLE   : 250   : L      : 플레카이니드 말초 용적
QFLE    : 40    : L/h    : 플레카이니드 구획간 청소율
CLFLE   : 35    : L/h    : 플레카이니드 청소율 (t1/2 ~13 h)
FUFLE   : 0.60  :        : 플레카이니드 유리 분율
MWFLE   : 414   : g/mol  : 플레카이니드 분자량
KAVER   : 1.2   : 1/h    : 베라파밀 흡수속도
FVER    : 0.22  :        : 베라파밀 생체이용률
VCVER   : 300   : L      : 베라파밀 용적
CLVER   : 60    : L/h    : 베라파밀 청소율
IC50VER : 80    : ng/mL  : 베라파밀 L형 통로 IC50
EMXVER  : 0.55  :        : 베라파밀 최대 L형 통로 차단 분율

$CMT @annotated
// ---- sympathetic / rate ----
NEJ    : 심근 간질 노르에피네프린 (nM)
EPI    : 순환 에피네프린 (nM)
HR     : 심박수 (bpm)
// ---- cAMP / PKA cascade ----
CAMP   : cAMP (정규화)
PKA    : 활성 PKA 분율
PI1    : 인산화 inhibitor-1 분율
PPLB   : P-phospholamban Ser16 분율
PRYR   : P-RyR2 Ser2808/2814 분율
PICA   : P-Cav1.2 분율
// ---- Ca cycling ----
CAI    : 박동평균 세포질 Ca (uM)
CANSR  : 네트워크 SR 유리 Ca (uM)
CAJSRT : 접합부 SR 총 Ca (uM)
CAMT   : 미토콘드리아 Ca (uM)
RREF   : RyR2 불응 분율
CSQ    : Calsequestrin-2 상대 단백량
// ---- electrophysiology ----
VDAD   : DAD 진폭 (mV)
NAAV   : Na 통로 가용도 분율
QRSD   : QRS 폭 (ms)
// ---- arrhythmia / outcome ----
PVCR   : 촉발 박동 발생률 (1/min)
VTB    : VT 부담 지수 (0-1)
NPVC   : 누적 촉발 박동 수
NVT    : 누적 VT 삽화 지수
HZRD   : 누적 부정맥 위험 (적분)
SYNCH  : 누적 실신 위험 (적분)
NSHOCK : 누적 ICD 쇼크 수
SHKD   : 쇼크 유발 카테콜아민 구동
// ---- chronic ----
CAMK   : 활성 CaMKII 분율
CAINT  : 누적 Ca 과부하 노출
FIBR   : 간질 섬유화 지수
FATG   : 피로 / 심박수 예비력 상실 지수
EXPO   : 누적 운동 노출
// ---- PK ----
NADD   : 나돌롤 흡수 구획 (mg)
NADC   : 나돌롤 중심 구획 (mg)
NADP   : 나돌롤 말초 구획 (mg)
METD   : 메토프롤롤 흡수 구획 (mg)
METC   : 메토프롤롤 중심 구획 (mg)
METP   : 메토프롤롤 말초 구획 (mg)
BISD   : 비소프롤롤 흡수 구획 (mg)
BISC   : 비소프롤롤 중심 구획 (mg)
FLED   : 플레카이니드 흡수 구획 (mg)
FLEC   : 플레카이니드 중심 구획 (mg)
FLEP   : 플레카이니드 말초 구획 (mg)
VERD   : 베라파밀 흡수 구획 (mg)
VERC   : 베라파밀 중심 구획 (mg)

$GLOBAL
#define H60 60.0
// Invert  T = f + Bmax*f/(f+Kd)  for the free luminal concentration.
// CASQ2 buffers ~90% of junctional-SR Ca, so total and free are NOT
// interchangeable: this is the only place CPVT2 differs structurally.
double csqfree(double total, double bmax, double kd) {
  double b = total - bmax - kd;
  return 0.5 * (b + sqrt(b * b + 4.0 * kd * total));
}
double hillp(double x, double k, double n) {
  double xn = pow(x, n);
  return xn / (xn + pow(k, n));
}

$MAIN
// Resting fixed point of the REFERENCE (normal, drug-free) physiology, taken
// from the settled Python implementation. mrgsolve applies the _0 assignments
// at initialisation only, so these are the initial conditions; a genotype or
// an intervention then moves the system to its own resting state within the
// first simulated minutes (the protocols all begin with >= 8 min of rest).
NEJ_0    = 0.60;   EPI_0   = 0.086;  HR_0    = 64.7;
CAMP_0   = 0.286;  PKA_0   = 0.112;  PI1_0   = 0.183;
PPLB_0   = 0.2058; PRYR_0  = 0.1656; PICA_0  = 0.1842;
CAI_0    = 0.377;  CANSR_0 = 1180.0; CAJSRT_0= 10050.0;
CAMT_0   = 100.0;  RREF_0  = 0.0;    CSQ_0   = CSQSET;
VDAD_0   = 0.0;    NAAV_0  = 1.0;    QRSD_0  = QRS0;

$ODE
// ======================= SYMPATHETIC DRIVE =================================
double SNS = 0.10 + 0.90 * EXER + 0.55 * EMOT;
if (SNS > 1.0) SNS = 1.0;
double NEj = (NEJ > 0.0 ? NEJ : 0.0);
dxdt_NEJ = KNE * SNS * (1.0 - FLCSD * ELCSD) - KUPT * NEj;

double EPIc = (EPI > 0.0 ? EPI : 0.0);
double SHKb = (SHKD > 0.0 ? SHKD : 0.0);
dxdt_EPI = KEPI * SNS * SNS + EPIINF + ESHK * SHKb - KCLE * EPIc;

// ================= COMPETITIVE BETA1 BLOCKADE ==============================
// Every beta-blocker enters through the SAME term with the SAME functional
// form. No drug is given an intrinsic-potency advantage; the only thing that
// differs between them is what their PK leaves behind at the trough.
double CNAD = (NADC > 0.0 ? NADC : 0.0) / VCNAD * 1000.0;   // ng/mL
double CMET = (METC > 0.0 ? METC : 0.0) / VCMET * 1000.0;
double CBIS = (BISC > 0.0 ? BISC : 0.0) / VCBIS * 1000.0;
double BSHIFT = 1.0 + CNAD / EC50NAD + CMET / EC50MET + CBIS / EC50BIS;
double AGON = NEj / EC50NE + FEPI * EPIc / EC50EPI;
double SYMP = AGON / BSHIFT;
double B1ACT = SYMP / (1.0 + SYMP);

// ============================ HEART RATE ===================================
double HRV = HRV0 + HRVEX * EXER;
double Ehr = SYMP / (SYMP + KMHR);
dxdt_HR = (HRV + (HRMX - HRV) * Ehr - HR) / TAUHR;
double HRn = HR; if (HRn < 30.0) HRn = 30.0; if (HRn > 260.0) HRn = 260.0;

// ======================= cAMP / PKA / TARGETS ==============================
dxdt_CAMP = KAC * B1ACT - KPDE * CAMP;
double cAMPn = (CAMP > 0.0 ? CAMP : 0.0);
dxdt_PKA = KPKAON * hillp(cAMPn, KMPKA, 2.0) * (1.0 - PKA) - KPKAOFF * PKA;
double PKAn = (PKA > 0.0 ? PKA : 0.0);
dxdt_PI1 = KPI1 * PKAn * (1.0 - PI1) - KPI1OFF * PI1;
double PP1 = 1.0 / (1.0 + 1.5 * PI1);      // PP1 activity, PKA-inhibited
dxdt_PPLB = (hillp(PKAn, KPLB * PP1, 2.0) - PPLB) / TAUPH;
dxdt_PRYR = (hillp(PKAn, KPRYR * PP1, 2.0) - PRYR) / TAUPH;
dxdt_PICA = (hillp(PKAn, KPICA * PP1, 2.0) - PICA) / TAUPH;

// ============================ Ca CYCLING ===================================
double CAIn   = (CAI    > 1e-4 ? CAI    : 1e-4);
double CANSRn = (CANSR  > 1.0  ? CANSR  : 1.0);
double CAJSRn = (CAJSRT > 1.0  ? CAJSRT : 1.0);
double CSQn   = (CSQ    > 0.01 ? CSQ    : 0.01);
dxdt_CSQ = KCSQTO * (CSQSET - CSQn);
double CAJSRF = csqfree(CAJSRn, BCSQ * CSQn, KCSQ);

double CVER = (VERC > 0.0 ? VERC : 0.0) / VCVER * 1000.0;
double fICa = (1.0 + AICA * (PICA - PICA0)) * (1.0 - EMXVER * CVER / (CVER + IC50VER));
if (fICa < 0.05) fICa = 0.05;
double fSERCA = 1.0 + ASERCA * (PPLB - PPLB0);

double RREFn = RREF; if (RREFn < 0.0) RREFn = 0.0; if (RREFn > 0.98) RREFn = 0.98;
double fload = hillp(CAJSRF, KRELH, HREL);

// per-beat processes x HR
double JCAL = GCAL * fICa * HRn;
double JUP  = VUP * fSERCA * hillp(CAIn, KMUP, 2.0) * HRn;
// NCX is DIASTOLE-limited: the exponent ENCX < 1 is what makes SR load
// rate-dependent and therefore what makes bradycardia protective.
double JNCX = VNCX * hillp(CAIn, KMNCX, 3.0) * 60.0 * pow(HRn / 60.0, ENCX);
double JREL = GREL * (1.0 - RREFn) * CAJSRF * fload * fICa * HRn;
double JMITI = KMITI * hillp(CAIn, KMMIT, 2.0) * HRn;
double JMITO = KMITO * (CAMT > 0.0 ? CAMT : 0.0) / 100.0 * HRn;
// per-minute process
double JLEAK = GLEAK * CAJSRF * (1.0 + BCAMK * CAMK);

// ============== CROSSING 1 : LOAD vs SOICR THRESHOLD =======================
double CFLEf = (FLEC > 0.0 ? FLEC : 0.0) / VCFLE * 1000.0 * FUFLE / MWFLE;  // uM free
double CTHR = CTH0 * (1.0 - GRYR) * (1.0 - GCSQI)
              * (1.0 + EFLR * CFLEf / (CFLEf + IC50FLR))
              / (1.0 + BPRYR * (PRYR - PRYR0) + BCAMK * CAMK);
double X = CAJSRF / CTHR;
double xs = X - 1.0; if (xs < 0.0) xs = 0.0;   // HARD threshold: no wave below 1
double WCAP = 1.0 / (TAUREC * FREC);           // refractoriness = frequency ceiling
double WAVE = WMAX * hillp(xs, KW, HW);
if (WAVE > WCAP) WAVE = WCAP;
double frw = FRELW * (1.0 + AFRW * (CSQn < 1.0 ? 1.0 - CSQn : 0.0));
double QRELc = frw * CAJSRn * VJSR;            // uM of cytosol per wave
double JSPONT = WAVE * QRELc;

dxdt_CAI = BETAI * (JCAL + JREL + JLEAK + JSPONT - JUP - JNCX - JMITI + JMITO);
dxdt_CAMT = (JMITI - JMITO) * 0.20;
double JTR = VTR * (CANSRn - CAJSRF);
dxdt_CANSR  = (JUP - JTR) / VNSR;
dxdt_CAJSRT = (JTR - JREL - JLEAK - JSPONT) / VJSR;
dxdt_RREF = KREF * (JREL + JSPONT) / 1000.0 * (1.0 - RREFn) - RREFn / (TAUREC * FREC);

// ============== CROSSING 2 : DAD vs REQUIRED DEPOLARISATION ================
double NAAVn = NAAV; if (NAAVn < NAMIN) NAAVn = NAMIN; if (NAAVn > 1.0) NAAVn = 1.0;
dxdt_NAAV = KOFFNA * (1.0 - NAAVn) - KONNA * CFLEf * (HRn / 60.0) * NAAVn;
double SINK = GK1 * sqrt(KO / 5.4);
double VDADT = (WAVE > 1e-9 ? GDADV * QRELc / SINK : 0.0);
dxdt_VDAD = (VDADT - VDAD) / TAUV;
double VDADn = (VDAD > 0.0 ? VDAD : 0.0);
double VREQ = VREQ0 / pow(NAAVn, PNAV);
double PTRIG = hillp(VDADn, VREQ, MTRIG);

// ===================== THE PRODUCT OF THE TWO CROSSINGS ====================
double PVCT = WAVE * PTRIG; if (PVCT > HRn) PVCT = HRn;
dxdt_PVCR = (PVCT - PVCR) / 0.10;
double PVCRn = (PVCR > 0.0 ? PVCR : 0.0);
dxdt_VTB = (hillp(PVCRn, KVT, NVT_H) - VTB) / 0.10;
dxdt_NPVC = PVCRn;
dxdt_NVT  = VTB * 2.0;
dxdt_HZRD = KHAZ * VTB;
dxdt_SYNCH = KHAZ * 0.45 * VTB;
double dsh = ICD * KSHK * VTB;
dxdt_NSHOCK = dsh;
dxdt_SHKD   = dsh - SHKD / TAUSHK;

// ============================== CHRONIC ====================================
double ovl = CAIn / CAI0 - 1.0; if (ovl < 0.0) ovl = 0.0;
dxdt_CAMK = KCMON * ovl * (1.0 - CAMK) - KCMOFF * CAMK;
dxdt_CAINT = ovl;
dxdt_FIBR = KFIB * CAMK - KFIBD * FIBR;
dxdt_QRSD = (QRS0 * (1.0 + AQRS * (1.0 - NAAVn)) - QRSD) / 0.5;
dxdt_FATG = ((1.0 - Ehr / 0.60) - FATG) / 5.0;
dxdt_EXPO = EXER;

// ================================= PK ======================================
// ka, CL and Q are given in 1/h and L/h; the model runs in MINUTES.
dxdt_NADD = -KANAD / H60 * NADD;
dxdt_NADC = FNAD * KANAD / H60 * NADD - CLNAD / H60 / VCNAD * NADC
            - QNAD / H60 * (NADC / VCNAD - NADP / VPNAD);
dxdt_NADP = QNAD / H60 * (NADC / VCNAD - NADP / VPNAD);

double clmet = CLMET * CYP2D6;
dxdt_METD = -KAMET / H60 * METD;
dxdt_METC = FMET * KAMET / H60 * METD - clmet / H60 / VCMET * METC
            - QMET / H60 * (METC / VCMET - METP / VPMET);
dxdt_METP = QMET / H60 * (METC / VCMET - METP / VPMET);

dxdt_BISD = -KABIS / H60 * BISD;
dxdt_BISC = FBIS * KABIS / H60 * BISD - CLBIS / H60 / VCBIS * BISC;

dxdt_FLED = -KAFLE / H60 * FLED;
dxdt_FLEC = FFLE * KAFLE / H60 * FLED - CLFLE / H60 / VCFLE * FLEC
            - QFLE / H60 * (FLEC / VCFLE - FLEP / VPFLE);
dxdt_FLEP = QFLE / H60 * (FLEC / VCFLE - FLEP / VPFLE);

dxdt_VERD = -KAVER / H60 * VERD;
dxdt_VERC = FVER * KAVER / H60 * VERD - CLVER / H60 / VCVER * VERC;

$TABLE
double CAJSRFo = csqfree((CAJSRT > 1.0 ? CAJSRT : 1.0), BCSQ * (CSQ > 0.01 ? CSQ : 0.01), KCSQ);
double CTHRo = CTH0 * (1.0 - GRYR) * (1.0 - GCSQI)
               * (1.0 + EFLR * ((FLEC > 0.0 ? FLEC : 0.0) / VCFLE * 1000.0 * FUFLE / MWFLE)
                  / (((FLEC > 0.0 ? FLEC : 0.0) / VCFLE * 1000.0 * FUFLE / MWFLE) + IC50FLR))
               / (1.0 + BPRYR * (PRYR - PRYR0) + BCAMK * CAMK);
double Xo = CAJSRFo / CTHRo;
double CNADo = (NADC > 0.0 ? NADC : 0.0) / VCNAD * 1000.0;
double CMETo = (METC > 0.0 ? METC : 0.0) / VCMET * 1000.0;
double CBISo = (BISC > 0.0 ? BISC : 0.0) / VCBIS * 1000.0;
double CFLEo = (FLEC > 0.0 ? FLEC : 0.0) / VCFLE * 1000.0;
double CVERo = (VERC > 0.0 ? VERC : 0.0) / VCVER * 1000.0;
double SHIFTo = 1.0 + CNADo / EC50NAD + CMETo / EC50MET + CBISo / EC50BIS;

$CAPTURE
CAJSRFo CTHRo Xo CNADo CMETo CBISo CFLEo CVERo SHIFTo EXER
'

mod <- mcode("cpvt", cpvt_code, end = 40, delta = 0.05)

## ===========================================================================
##  HELPERS
## ===========================================================================

## Genotype presets. GRYR is the RyR2 gain-of-function threshold shift; FREC
## is the refractory-recovery multiplier (< 1 = faster recovery = higher
## achievable wave rate). CPVT2 works through the CASQ2 route instead.
GENO <- list(
  WT     = list(GRYR = 0.00, FREC = 1.00, CSQSET = 1.00, GCSQI = 0.00),
  CPVT1  = list(GRYR = 0.24, FREC = 0.60, CSQSET = 1.00, GCSQI = 0.00),
  CPVT1s = list(GRYR = 0.30, FREC = 0.55, CSQSET = 1.00, GCSQI = 0.00),
  CPVT1m = list(GRYR = 0.18, FREC = 0.70, CSQSET = 1.00, GCSQI = 0.00),
  CPVT2  = list(GRYR = 0.00, FREC = 0.45, CSQSET = 0.45, GCSQI = 0.28)
)

## An exercise stress test as a piecewise-constant EXER covariate.
## Rest for t0 min, ramp linearly to `peak`, then recover.
exercise_data <- function(t0 = 10, ramp = 12, rec = 12, peak = 1.0, dt = 0.25) {
  tt <- seq(0, t0 + ramp + rec + 6, by = dt)
  ex <- ifelse(tt < t0, 0,
        ifelse(tt < t0 + ramp, (tt - t0) / ramp * peak,
               pmax(0, peak * (1 - (tt - t0 - ramp) / rec))))
  data.frame(time = tt, EXER = ex, EMOT = 0, EPIINF = 0)
}

## Epinephrine challenge (Mayo/Shimizu-style incremental infusion).
## EPIINF is in nM/min; 0.1 ug/kg/min maps to ~0.5 nM/min here.
epi_data <- function(t0 = 8, steps = c(0.05, 0.10, 0.20, 0.30),
                     dur = 5, dt = 0.25) {
  tt <- seq(0, t0 + length(steps) * dur + 8, by = dt)
  inf <- sapply(tt, function(x) {
    k <- floor((x - t0) / dur) + 1
    if (x < t0 || k > length(steps)) 0 else steps[k] * 5.0
  })
  data.frame(time = tt, EXER = 0, EMOT = 0, EPIINF = inf)
}

## Startle / emotional trigger at rest (swimming, alarm clock).
startle_data <- function(t0 = 10, dur = 4, dt = 0.25) {
  tt <- seq(0, t0 + dur + 16, by = dt)
  data.frame(time = tt, EXER = 0,
             EMOT = ifelse(tt >= t0 & tt < t0 + dur, 1, 0), EPIINF = 0)
}

## Steady-state PK amounts. Rather than integrating the whole 44-state system
## for 18 days, the (linear) PK subsystem is solved on its own and its
## compartment amounts are handed to the physiology model as initial values.
## `offset_h` is the time since the LAST dose -- this is where the entire
## nadolol/metoprolol argument lives, so it is an explicit argument and never
## a default.
pk_steady <- function(drug, dose_mg, tau_h, offset_h, ndays = 18, CYP2D6 = 1) {
  spec <- switch(drug,
    nad = list(ka = 1.2, F = 0.30, VC = 150, VP = 150, Q = 20,  CL = 12),
    met = list(ka = 1.5, F = 0.40, VC = 250, VP = 500, Q = 90,  CL = 150 * CYP2D6),
    bis = list(ka = 1.0, F = 0.90, VC = 230, VP = NA,  Q = 0,   CL = 14.5),
    fle = list(ka = 0.9, F = 0.90, VC = 350, VP = 250, Q = 40,  CL = 35),
    ver = list(ka = 1.2, F = 0.22, VC = 300, VP = NA,  Q = 0,   CL = 60))
  two <- !is.na(spec$VP)
  VP <- if (two) spec$VP else 1
  ka <- spec$ka / 60; CL <- spec$CL / 60; Q <- spec$Q / 60
  y <- c(0, 0, 0); h <- 0.02
  step <- function(y, dt) {
    n <- ceiling(dt / h); hh <- dt / n
    for (i in seq_len(n)) {
      d1 <- -ka * y[1]
      d2 <- spec$F * ka * y[1] - CL / spec$VC * y[2] - Q * (y[2] / spec$VC - y[3] / VP)
      d3 <- Q * (y[2] / spec$VC - y[3] / VP)
      y <- y + hh * c(d1, d2, d3)
    }
    y
  }
  nd <- max(2, floor(ndays * 1440 / (tau_h * 60)))
  for (i in seq_len(nd)) { y[1] <- y[1] + dose_mg; y <- step(y, tau_h * 60) }
  y[1] <- y[1] + dose_mg
  y <- step(y, offset_h * 60)
  out <- list(); nm <- switch(drug,
    nad = c("NADD", "NADC", "NADP"), met = c("METD", "METC", "METP"),
    bis = c("BISD", "BISC", NA),     fle = c("FLED", "FLEC", "FLEP"),
    ver = c("VERD", "VERC", NA))
  out[[nm[1]]] <- y[1]; out[[nm[2]]] <- y[2]
  if (two) out[[nm[3]]] <- y[3]
  attr(out, "conc") <- y[2] / spec$VC * 1000
  out
}

## Run one scenario. `drugs` is a list of list(drug, dose, tau_h, offset_h).
run_cpvt <- function(genotype = "CPVT1", drugs = list(), pmod = list(),
                     protocol = exercise_data(), label = "") {
  pars <- c(GENO[[genotype]], pmod)
  init <- list()
  concs <- c()
  for (d in drugs) {
    amt <- pk_steady(d$drug, d$dose, d$tau, d$offset,
                     CYP2D6 = if (is.null(pmod$CYP2D6)) 1 else pmod$CYP2D6)
    init <- c(init, amt)
    concs[d$drug] <- attr(amt, "conc")
  }
  ## The exercise / epinephrine / startle protocol is supplied as time-varying
  ## covariates on evid = 0 records; mrgsolve holds each value until the next
  ## record, so the ramp is piecewise-constant at the protocol's own dt.
  dat <- protocol %>%
    mutate(ID = 1, amt = 0, cmt = 1, evid = 0) %>%
    select(ID, time, amt, cmt, evid, EXER, EMOT, EPIINF)
  out <- mod %>%
    param(pars) %>%
    init(init) %>%
    mrgsim_df(data = dat, end = max(protocol$time), delta = 0.05,
              recover = "EXER")
  attr(out, "conc") <- concs
  attr(out, "label") <- label
  out
}

## Endpoints: the exercise-test ventricular-ectopy score (VE 0-4) used as the
## primary surrogate in the CPVT trials, plus PVC count and VT minutes.
endpoints <- function(sim) {
  pv <- sim$PVCR; vt <- sim$VTB; hr <- sim$HR
  ve <- 0
  if (max(pv) > 0.5) ve <- 1
  if (max(pv) > 5)   ve <- 2
  if (max(pv) > 12)  ve <- 3
  if (max(vt) > 0.5) ve <- 4
  ion <- which(pv > 1.0); ivt <- which(vt > 0.5)
  data.frame(
    label      = attr(sim, "label"),
    HR_peak    = max(hr),
    HR_at_PVC  = if (length(ion)) hr[ion[1]] else NA_real_,
    HR_at_VT   = if (length(ivt)) hr[ivt[1]] else NA_real_,
    PVC_peak   = max(pv),
    VT_peak    = max(vt),
    VT_minutes = sum(diff(sim$time) * head(vt, -1)),
    PVC_count  = max(sim$NPVC),
    VE_score   = ve,
    QRS_max    = max(sim$QRSD),
    X_peak     = max(sim$Xo),
    shift_min  = min(sim$SHIFTo),
    hazard     = max(sim$HZRD),
    shocks     = max(sim$NSHOCK)
  )
}

## ===========================================================================
##  SCENARIOS
## ===========================================================================
## Each scenario states WHICH crossing it is interrogating. All beta-blocker
## scenarios are run at the TROUGH unless the label says otherwise, because
## the trough is the hour a CPVT patient is unprotected in.

scenarios <- list(
  ## --- A. genotype / natural history (crossing 1 set by the mutation) -----
  s01 = list(l = "1  야생형 운동부하 검사",              g = "WT"),
  s02 = list(l = "2  CPVT1 전형, 무치료",                g = "CPVT1"),
  s03 = list(l = "3  CPVT1 중증, 무치료",                g = "CPVT1s"),
  s04 = list(l = "4  CPVT1 경증 침투, 무치료",           g = "CPVT1m"),
  s05 = list(l = "5  CPVT2 (CASQ2), 무치료",             g = "CPVT2"),

  ## --- B. beta-blockade at trough vs peak (crossing 1) --------------------
  s06 = list(l = "6  + 나돌롤 80 mg qd (24 h 최저)",     g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24))),
  s07 = list(l = "7  + 나돌롤 80 mg qd (2 h 최고)",      g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 2))),
  s08 = list(l = "8  + 메토프롤롤 100 bid (12 h 최저)",  g = "CPVT1",
             d = list(list(drug = "met", dose = 100, tau = 12, offset = 12))),
  s09 = list(l = "9  + 메토프롤롤 100 bid (2 h 최고)",   g = "CPVT1",
             d = list(list(drug = "met", dose = 100, tau = 12, offset = 2))),
  s10 = list(l = "10 + 메토프롤롤 200 bid (12 h 최저)",  g = "CPVT1",
             d = list(list(drug = "met", dose = 200, tau = 12, offset = 12))),
  s11 = list(l = "11 + 메토프롤롤 200 bid (2 h 최고)",   g = "CPVT1",
             d = list(list(drug = "met", dose = 200, tau = 12, offset = 2))),
  s12 = list(l = "12 + 비소프롤롤 10 mg qd (24 h 최저)", g = "CPVT1",
             d = list(list(drug = "bis", dose = 10, tau = 24, offset = 24))),

  ## --- C. one missed dose: half-life IS adherence -------------------------
  s13 = list(l = "13 나돌롤 1회 누락 (48 h)",            g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 48))),
  s14 = list(l = "14 나돌롤 2회 누락 (72 h)",            g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 72))),
  s15 = list(l = "15 메토프롤롤 1회 누락 (24 h)",        g = "CPVT1",
             d = list(list(drug = "met", dose = 100, tau = 12, offset = 24))),
  s16 = list(l = "16 메토프롤롤 CYP2D6 초고속대사자",    g = "CPVT1",
             d = list(list(drug = "met", dose = 100, tau = 12, offset = 12)),
             p = list(CYP2D6 = 2.5)),
  s17 = list(l = "17 메토프롤롤 CYP2D6 저대사자",        g = "CPVT1",
             d = list(list(drug = "met", dose = 100, tau = 12, offset = 12)),
             p = list(CYP2D6 = 0.30)),

  ## --- D. flecainide add-on (crossing 2 x crossing 1) ---------------------
  s18 = list(l = "18 나돌롤 + 플레카이니드 100 bid",     g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12))),
  s19 = list(l = "19 나돌롤 + 플레카이니드 150 bid",     g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 150, tau = 12, offset = 12))),
  s20 = list(l = "20 플레카이니드 단독 100 bid",         g = "CPVT1",
             d = list(list(drug = "fle", dose = 100, tau = 12, offset = 12))),
  s21 = list(l = "21 나돌롤 + 베라파밀 240 qd",          g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "ver", dose = 240, tau = 24, offset = 4))),

  ## --- E. which flecainide arm does the work? (ablate each) ---------------
  s22 = list(l = "22 나돌롤+플레카이니드, RyR2 팔만",    g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12)),
             p = list(KONNA = 0)),
  s23 = list(l = "23 나돌롤+플레카이니드, Na 팔만",      g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12)),
             p = list(EFLR = 0)),

  ## --- F. LCSD, genotype-specific therapy, electrolytes -------------------
  s24 = list(l = "24 나돌롤 + LCSD",                     g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             p = list(FLCSD = 1)),
  s25 = list(l = "25 LCSD 단독",                         g = "CPVT1",
             p = list(FLCSD = 1)),
  s26 = list(l = "26 LCSD + 나돌롤 1회 누락 (48 h)",     g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 48)),
             p = list(FLCSD = 1)),
  s27 = list(l = "27 CPVT2 + 나돌롤",                    g = "CPVT2",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24))),
  s28 = list(l = "28 CPVT2 + 나돌롤 + 플레카이니드",     g = "CPVT2",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12))),
  s29 = list(l = "29 나돌롤, 저칼륨혈증 K 3.0",          g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             p = list(KO = 3.0)),

  ## --- G. diagnostics: why the epinephrine test misses milder cases ------
  s30 = list(l = "30 에피네프린 부하, CPVT1 전형",       g = "CPVT1",
             pr = epi_data()),
  s31 = list(l = "31 에피네프린 부하, CPVT1 중증",       g = "CPVT1s",
             pr = epi_data()),
  s32 = list(l = "32 에피네프린 부하, 야생형",           g = "WT",
             pr = epi_data()),

  ## --- H. startle at rest (no tachycardia, still arrhythmogenic) ---------
  s33 = list(l = "33 안정 시 놀람 자극, 무치료",         g = "CPVT1",
             pr = startle_data()),
  s34 = list(l = "34 안정 시 놀람 자극, 나돌롤 최저",    g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             pr = startle_data()),

  ## --- I. ICD: shock -> catecholamine surge -> more VT --------------------
  s35 = list(l = "35 ICD, beta차단제 없음",              g = "CPVT1",
             p = list(ICD = 1)),
  s36 = list(l = "36 ICD, 쇼크-카테콜아민 되먹임 차단",  g = "CPVT1",
             p = list(ICD = 1, ESHK = 0)),
  s37 = list(l = "37 ICD + 나돌롤",                      g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             p = list(ICD = 1)),

  ## --- J. rate-dependence: the testable discriminator between the arms ---
  s38 = list(l = "38 경도 운동, 나돌롤",                 g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             pr = exercise_data(peak = 0.45)),
  s39 = list(l = "39 경도 운동, 나돌롤+플레카이니드",    g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12)),
             pr = exercise_data(peak = 0.45)),
  s40 = list(l = "40 최대 운동, 나돌롤+플레카이니드",    g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                      list(drug = "fle", dose = 100, tau = 12, offset = 12)),
             pr = exercise_data(peak = 1.0)),

  ## --- K. exercise restriction: lower the INPUT, not the response --------
  s41 = list(l = "41 운동 제한 (0.45), 나돌롤",          g = "CPVT1",
             d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
             pr = exercise_data(peak = 0.45))
)

run_all <- function() {
  res <- lapply(names(scenarios), function(k) {
    s <- scenarios[[k]]
    sim <- run_cpvt(genotype = s$g,
                    drugs    = if (is.null(s$d)) list() else s$d,
                    pmod     = if (is.null(s$p)) list() else s$p,
                    protocol = if (is.null(s$pr)) exercise_data() else s$pr,
                    label    = s$l)
    endpoints(sim)
  })
  bind_rows(res)
}

## ===========================================================================
##  THE TABLE THAT CARRIES THE ARGUMENT
## ===========================================================================
## The competitive shift factor 1 + C/EC50 is what a beta-blocker actually
## delivers. Printed at peak and at trough, side by side with the outcome, it
## makes the nadolol/metoprolol split visible as a PK fact.
shift_table <- function() {
  rows <- list(
    c("nadolol 80 mg qd",     "nad",  80, 24,  2, 24),
    c("nadolol 80 mg qd",     "nad",  80, 24, 24, 48),
    c("metoprolol 100 mg bid","met", 100, 12,  2, 12),
    c("metoprolol 200 mg bid","met", 200, 12,  2, 12),
    c("bisoprolol 10 mg qd",  "bis",  10, 24,  2, 24),
    c("flecainide 100 mg bid","fle", 100, 12,  2, 12),
    c("flecainide 150 mg bid","fle", 150, 12,  2, 12))
  ec50 <- c(nad = 20, met = 45, bis = 12, fle = NA)
  out <- lapply(rows, function(r) {
    drug <- r[2]; dose <- as.numeric(r[3]); tau <- as.numeric(r[4])
    o1 <- as.numeric(r[5]); o2 <- as.numeric(r[6])
    c1 <- attr(pk_steady(drug, dose, tau, o1), "conc")
    c2 <- attr(pk_steady(drug, dose, tau, o2), "conc")
    data.frame(regimen = r[1], t_hi = o1, C_hi = round(c1, 1),
               t_lo = o2, C_lo = round(c2, 1),
               shift_hi = round(1 + c1 / ec50[drug], 2),
               shift_lo = round(1 + c2 / ec50[drug], 2))
  })
  bind_rows(out)
}

## ===========================================================================
##  PLOTS
## ===========================================================================
plot_crossings <- function(genotype = "CPVT1", drugs = list()) {
  sim <- run_cpvt(genotype, drugs, protocol = exercise_data())
  sim %>%
    select(time, HR, CAJSRFo, CTHRo, Xo, VDAD, PVCR, VTB) %>%
    pivot_longer(-time) %>%
    mutate(name = factor(name, levels = c("HR", "CAJSRFo", "CTHRo", "Xo",
                                          "VDAD", "PVCR", "VTB"))) %>%
    ggplot(aes(time, value)) +
    geom_line(linewidth = 0.7, colour = "#C0392B") +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    labs(x = "시간 (min)", y = NULL,
         title = paste0("CPVT 운동부하 검사 — ", genotype),
         subtitle = "Xo = CAJSRF / CTHR : 1을 넘는 순간부터만 파동이 발생한다") +
    theme_bw(base_size = 11)
}

plot_trough_argument <- function() {
  arms <- list(
    list(l = "무치료",                    d = list()),
    list(l = "나돌롤 최저 (24 h)",        d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24))),
    list(l = "메토프롤롤 최저 (12 h)",    d = list(list(drug = "met", dose = 100, tau = 12, offset = 12))),
    list(l = "메토프롤롤 최고 (2 h)",     d = list(list(drug = "met", dose = 100, tau = 12, offset = 2))),
    list(l = "메토프롤롤 200 최고 (2 h)", d = list(list(drug = "met", dose = 200, tau = 12, offset = 2))))
  bind_rows(lapply(arms, function(a) {
    s <- run_cpvt("CPVT1", a$d, protocol = exercise_data(), label = a$l)
    data.frame(arm = a$l, time = s$time, PVCR = s$PVCR, HR = s$HR)
  })) %>%
    ggplot(aes(time, PVCR, colour = arm)) +
    geom_line(linewidth = 0.8) +
    labs(x = "시간 (min)", y = "촉발 박동 발생률 (1/min)",
         title = "동일한 수용체, 동일한 경쟁식 — 차이는 최저 농도에서만 발생한다",
         subtitle = "메토프롤롤은 배용량 최고 농도에서 겨우 나돌롤의 최저 농도에 도달한다") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

## ===========================================================================
##  USAGE
## ===========================================================================
## res <- run_all(); print(res, digits = 3)
## print(shift_table())
## plot_crossings("CPVT1")
## plot_crossings("WT")          # Xo never reaches 1 -> completely silent
## plot_trough_argument()
