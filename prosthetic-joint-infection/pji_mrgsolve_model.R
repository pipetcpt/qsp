## =============================================================================
##  pji_mrgsolve_model.R
##  인공관절 감염 (Prosthetic Joint Infection, PJI) QSP 모델
##  Implant-associated Staphylococcus aureus osteomyelitis
##
##  53 ODEs · 6 antibiotics with explicit bone penetration · 5 bacterial
##  phenotypic states · 2 resistant lineages · local elution from an
##  antibiotic-loaded cement spacer · innate immunity with frustrated
##  phagocytosis · RANKL/OPG osteolysis · 2018 ICM biomarkers · drug toxicity
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL IS ORGANISED AROUND
##  ---------------------------------------------------------------------------
##  PJI doctrine looks like a list of empirical rules ("you must debride",
##  "rifampicin is special", "never give rifampicin alone", "DAIR only works
##  early"). This model is written so that all four fall out of arithmetic
##  rather than being coded in as assertions.
##
##  PILLAR 1 — THE FOREIGN BODY MOVES THE INFECTIVE DOSE, NOT THE ORGANISM.
##    Two mechanisms, both explicit: (a) FBFACT scales down phagocyte killing
##    next to metal, and (b) KATT moves planktonic cells onto a surface where
##    phagocytosis is "frustrated" (FRUST = 0.03). Neither is a threshold
##    parameter. Simulated: an implant-free niche clears up to ~10^7 CFU;
##    with an implant, ~10^1 CFU establishes. That ~10^5-10^6 fold shift is
##    the model's version of Elek & Conen's suture experiment (1957) and of
##    Zimmerli's tissue-cage model (1982).
##
##  PILLAR 2 — EVERY ANTIBIOTIC IS ITS FREE BONE CONCENTRATION DIVIDED BY MBEC.
##    The model never uses the planktonic MIC to decide what happens in the
##    biofilm. For each drug it carries (i) AUC_bone/AUC_plasma, (ii) the free
##    fraction, (iii) the planktonic MIC and (iv) the biofilm MBEC. The ratio
##    C_bone,free / MBEC at standard adult doses then computes to:
##
##       vancomycin 1 g q12h ...... 2.083 / 512 = 0.0041
##       cefazolin 2 g q8h (MSSA).. 2.500 / 256 = 0.0098
##       daptomycin 8 mg/kg qd .... 0.350 / 32  = 0.0109
##       linezolid 600 mg q12h .... 2.218 / 128 = 0.0173
##       levofloxacin 750 mg qd ... 1.203 / 64  = 0.0188
##       RIFAMPICIN 450 mg q12h ... 0.281 / 1.0 = 0.281 (week 1)
##                                  0.141 / 1.0 = 0.141 (autoinduced steady state)
##
##    (all six computed by ratio_table(); rifampicin sits 7-68x above every
##     other option and is still BELOW 1 -- which is why it needs surgery.)
##
##    Rifampicin is not a better antibiotic. It is the only one whose target
##    (RNA polymerase) is still load-bearing in an adherent, barely dividing
##    cell, so its MBEC never climbs 2-3 logs above its MIC the way a cell-wall
##    agent's does. The entire clinical special status of rifampicin in
##    implant infection is that one column of the table.
##
##  PILLAR 3 — SURGERY IS A MUTANT-SUPPLY OPERATION, NOT A CLEANING OPERATION.
##    P(a pre-existing rpoB mutant) = 1 - exp(-mu * N), mu ~ 1e-8 per division.
##    A mature biofilm carries N ~ 1e9-1e10, so that probability is ~1.000 and
##    rifampicin is dead on arrival. A 4-log debridement takes N to ~1e5-1e6
##    and the probability to ~1e-3-1e-2. Debridement does not "help" the
##    antibiotic; it is the step that makes Pillar 2 usable at all. The model
##    computes this probability continuously (PRPOB) and also carries explicit
##    rpoB (RP/RB) and gyrA (QP/QB) lineages generated from the growth flux.
##
##  PILLAR 4 — BIOFILM MATURITY IS A CLOCK, AND DAIR RACES IT.
##    PHIB (the MBEC/MIC tolerance multiplier) is not constant: it is scaled by
##    matrix maturity EPS, which rises with a ~28-day time constant. A biofilm
##    debrided at day 7 is still ~40% "planktonic-like"; at day 90 it is fully
##    tolerant. The 3-week symptom-duration rule for DAIR is therefore an
##    emergent property of an EPS differential equation, not an if-statement.
##
##  ---------------------------------------------------------------------------
##  IMPORTANT MODELLING CHOICES / LIMITATIONS (read before using numbers)
##  ---------------------------------------------------------------------------
##  * Bone compartments are FREE-concentration effect compartments driven by
##    total plasma concentration through FB_x = (AUC_bone/AUC_plasma) x fu.
##    They do not feed mass back to plasma (periprosthetic mass is negligible).
##  * The model is deterministic. Bacterial extinction is enforced by a smooth
##    sub-single-cell decay term (KEXT below 1 CFU) so the ODE cannot carry
##    "0.001 of a bacterium" and regrow it, which is the classic deterministic
##    artefact in antimicrobial models. Stochastic outcomes (ID50, resistance
##    emergence, cure) are reported as PROBABILITIES computed from the
##    deterministic burden, not as trajectories.
##  * Cure probability is mapped as PCURE = exp(-PSEED * N_end): a surviving
##    burden of ~1 CFU next to a retained implant is roughly a coin flip.
##  * Cytokine and biomarker units are local/periprosthetic, not serum, except
##    CRP / ESR / PLT / SCR / ALT which are systemic.
##  * Calibration is to published trial-level endpoints (Zimmerli 1998 JAMA;
##    OVIVA 2019; DATIPO 2021; Lora-Tamayo 2013) and to standard PK/bone
##    penetration reviews (Landersdorfer 2009; Thabit 2019). It has NOT been
##    fitted to individual patient data. Educational / research use only.
##
##  ---------------------------------------------------------------------------
##  VERIFICATION (all 53 ODEs were re-implemented independently in
##  Python/scipy and run before this file was finalised). That pass exposed
##  and fixed six defects, each of which had silently broken one of the
##  pillars:
##
##  1. IMMUNE CAPACITY WRITTEN PER COMPARTMENT. A Michaelis-Menten sink on
##     each bacterial state made a rare rpoB clone face the ENTIRE phagocyte
##     capacity as first-order clearance at ~19 /h, sterilising exactly the
##     subpopulation Pillar 3 is about. Rewritten as one shared capacity
##     allocated in proportion to accessibility (NPHAG / SHR).
##  2. EXTINCTION FLOOR APPLIED PER STATE. Forcing every compartment below
##     1 CFU to zero annihilated the rpoB lineage on every dip. Extinction is
##     now a property of the CLONE. Before the fix an untreated 1e10 lesion
##     carried 0.01 mutants; after it, 63 -- against the analytic mu*N = 100.
##  3. DEBRIDEMENT RESET BIOFILM MATURITY. EPS was stripped in proportion to
##     the surgical log-kill, so a 2.5-log DAIR dropped EPS 0.62 -> 0.026 and
##     handed the residual biofilm near-planktonic susceptibility. Result:
##     levofloxacin MONOTHERAPY cured. A retained implant keeps its adherent
##     matrix (EPSCUT retention factor 0.55 -> 0.05).
##  4. NO CARRYING CAPACITY ON THE FREE POOL. Dispersed organisms grew to
##     1.8e10 and outnumbered the biofilm, making the infection spuriously
##     vancomycin-curable. NPCAP (fluid/abscess) and NICCAP (bounded by host
##     cell number) added.
##  5. UNBOUNDED LOOSENING AND OSTEOBLAST COLLAPSE. LOOSEN ran to 266 on a
##     0-100 scale and KAPOP drove osteoblasts to 8% of baseline.
##  6. PHAGOCYTE SATURATION CONSTANT INCONSISTENT WITH THE ID50. KMIMM = 1e4
##     against a target ID50 of ~1e7 made clearance first-order over three
##     logs with a large rate constant, so the immune system held the biofilm
##     at 1e4.8 CFU and no PJI ever established.
##
##  KNOWN CALIBRATION GAP (stated, not tuned away): the model's DAIR-timing
##  crossover falls between 3 and 6 months of infection duration, whereas the
##  clinical rule is ~3 weeks of symptoms. Debridement at day 21 sterilises on
##  day 60 of an 84-day course; at day 90 it sterilises on day 84 of 84 (zero
##  margin); at day 180 it fails outright. The direction and the mechanism
##  (EPS maturity ~10 d, sequestrum shielding ~3 mo) are right; the absolute
##  clock is too slow because the model has no soft-tissue abscess or
##  mechanical-loosening channel, which is how late DAIR actually fails.
##
##  Author: QSP Disease Model Library (Claude Code Routine)
## =============================================================================

library(mrgsolve)
library(dplyr)

pji_code <- '
$PROB
# Prosthetic Joint Infection (PJI) QSP model
# Implant-associated Staphylococcus aureus osteomyelitis
# 53 ODEs | biofilm PK/PD | mutant supply | surgery | osteolysis

$GLOBAL
// Emax kill helper: EMAX * C / (EC50eff + C)
#define EKILL(EMX, CC, EC) ((EMX)*(CC)/((EC)+(CC)))
// Biofilm tolerance multiplier scaled by matrix maturity (Pillar 4).
// EPSM in [0,1); a bare, freshly attached cell is nearly planktonic.
#define PHIMAT(PHI, EPSM) (1.0 + ((PHI) - 1.0)*(EPSM))

$PARAM @annotated
// ---------------------------------------------------------------- SETUP ----
IMPL    :  1    : 임플란트 존재 (1 = prosthesis in situ, 0 = implant-free bone)
MSSA    :  0    : 메티실린 감수성 (1 = MSSA, cefazolin active; 0 = MRSA)
IMSUP   :  1    : 숙주 면역 계수 (1 = 정상, 0.5 = 당뇨/스테로이드/RA)
INOC0   :  100  : CFU    : 접종 균량 (유주균으로 t=0에 투입)

// ------------------------------------------------- VANCOMYCIN PK (IV) ------
VCVAN   : 40    : L      : VAN 중심 분포용적
VPVAN   : 40    : L      : VAN 말초 분포용적
QVAN    :  8    : L/h    : VAN 구획간 청소율
CLVAN   :  4    : L/h    : VAN 전신 청소율 (CrCl ~100 mL/min)
FUVAN   :  0.50 :        : VAN 유리 분율
PENVAN  :  0.20 :        : VAN AUC_bone / AUC_plasma
KEQVAN  :  0.30 : 1/h    : VAN 골 구획 평형 속도

// ------------------------------------------------ RIFAMPICIN PK (PO) -------
KARIF   :  1.20 : 1/h    : RIF 흡수속도
FRIF    :  0.90 :        : RIF 생체이용률
VCRIF   : 55    : L      : RIF 분포용적
CLRIF0  : 12    : L/h    : RIF 기저 청소율 (자가유도 전)
FURIF   :  0.20 :        : RIF 유리 분율
PENRIF  :  0.50 :        : RIF AUC_bone / AUC_plasma
KEQRIF  :  0.50 : 1/h    : RIF 골 구획 평형 속도
KENZ    :  0.0072: 1/h   : CYP3A4 유도 전환속도 (t1/2 ~ 4일)
EIND    :  1.00 :        : RIF 자가유도 최대 배수 (청소율 x2)

// ---------------------------------------------- LEVOFLOXACIN PK (PO) -------
KALVX   :  1.50 : 1/h    : LVX 흡수속도
FLVX    :  0.99 :        : LVX 생체이용률
VCLVX   : 90    : L      : LVX 분포용적
CLLVX   :  9    : L/h    : LVX 청소율
FULVX   :  0.70 :        : LVX 유리 분율
PENLVX  :  0.50 :        : LVX AUC_bone / AUC_plasma
KEQLVX  :  0.50 : 1/h    : LVX 골 구획 평형 속도

// ------------------------------------------------ DAPTOMYCIN PK (IV) -------
VCDAP   : 10    : L      : DAP 분포용적
CLDAP   :  0.80 : L/h    : DAP 청소율
FUDAP   :  0.08 :        : DAP 유리 분율 (92% 알부민 결합)
PENDAP  :  0.15 :        : DAP AUC_bone / AUC_plasma
KEQDAP  :  0.25 : 1/h    : DAP 골 구획 평형 속도

// ------------------------------------------------- LINEZOLID PK (PO) -------
KALZD   :  1.80 : 1/h    : LZD 흡수속도
FLZD    :  1.00 :        : LZD 생체이용률
VCLZD   : 45    : L      : LZD 분포용적
CLLZD0  :  7    : L/h    : LZD 기저 청소율
FULZD   :  0.69 :        : LZD 유리 분율
PENLZD  :  0.45 :        : LZD AUC_bone / AUC_plasma
KEQLZD  :  0.50 : 1/h    : LZD 골 구획 평형 속도
RIFLZD  :  0.50 :        : RIF에 의한 LZD 청소율 유도 강도 (AUC -32%)

// -------------------------------------------------- CEFAZOLIN PK (IV) ------
VCCFZ   : 15    : L      : CFZ 분포용적
CLCFZ   :  4    : L/h    : CFZ 청소율
FUCFZ   :  0.20 :        : CFZ 유리 분율
PENCFZ  :  0.20 :        : CFZ AUC_bone / AUC_plasma
KEQCFZ  :  0.50 : 1/h    : CFZ 골 구획 평형 속도

// --------------------------------- LOCAL ELUTION (ALBC / SPACER) -----------
VJOINT  :  0.05 : L      : 유효 관절강·계면 용적
KFAST   :  0.05 : 1/h    : 표면층 폭발 방출 속도
KSLOW   :  0.0015: 1/h   : 기질 확산 방출 속도 (t1/2 ~ 19일)
KCLJ    :  0.35 : 1/h    : 관절액 약물 소실 속도
FLOCP   :  0.50 :        : 국소 약물의 유주균 접근도
FLOCB   :  0.15 :        : 국소 약물의 잔존 바이오필름 접근도

// -------------------------------------------- MIC / MBEC (S. aureus) -------
MICVAN  :  1.0  : mg/L   : VAN 유주균 MIC
MICRIF  :  0.012: mg/L   : RIF 유주균 MIC
MICLVX  :  0.25 : mg/L   : LVX 유주균 MIC
MICDAP  :  0.50 : mg/L   : DAP 유주균 MIC
MICLZD  :  2.0  : mg/L   : LZD 유주균 MIC
MICCFZ  :  1.0  : mg/L   : CFZ 유주균 MIC (MSSA)
MBCVAN  : 512   : mg/L   : VAN 바이오필름 MBEC
MBCRIF  :  1.0  : mg/L   : RIF 바이오필름 MBEC
MBCLVX  : 64    : mg/L   : LVX 바이오필름 MBEC
MBCDAP  : 32    : mg/L   : DAP 바이오필름 MBEC
MBCLZD  : 128   : mg/L   : LZD 바이오필름 MBEC
MBCCFZ  : 256   : mg/L   : CFZ 바이오필름 MBEC

// ------------------------------- MAXIMUM KILL RATES (planktonic) -----------
EMXVAN  :  0.45 : 1/h    : VAN 최대 살균속도 (시간의존·완만)
EMXRIF  :  0.80 : 1/h    : RIF 최대 살균속도
EMXLVX  :  1.20 : 1/h    : LVX 최대 살균속도 (농도의존)
EMXDAP  :  2.00 : 1/h    : DAP 최대 살균속도 (급속 농도의존)
EMXLZD  :  0.12 : 1/h    : LZD 최대 살균속도 (정균제)
EMXCFZ  :  0.80 : 1/h    : CFZ 최대 살균속도
EMXLZDI :  0.85 :        : LZD 최대 증식억제 분율 (정균 작용의 본체)

// ---- STATE TOLERANCE MULTIPLIERS (applied ON TOP of the biofilm PHI) ------
PHPVAN  : 10.0  :        : VAN 지속균 추가 내약 배수
PHPRIF  : 25.0  :        : RIF 지속균 추가 내약 배수
PHPLVX  : 25.0  :        : LVX 지속균 추가 내약 배수
PHPDAP  :  8.0  :        : DAP 지속균 추가 내약 배수
PHPLZD  : 20.0  :        : LZD 지속균 추가 내약 배수
PHPCFZ  : 40.0  :        : CFZ 지속균 추가 내약 배수 (증식 의존 최대)
PHSVAN  :  3.0  :        : VAN SCV 추가 내약 배수
PHSRIF  :  2.0  :        : RIF SCV 추가 내약 배수
PHSLVX  :  4.0  :        : LVX SCV 추가 내약 배수
PHSDAP  :  3.0  :        : DAP SCV 추가 내약 배수
PHSLZD  :  1.5  :        : LZD SCV 추가 내약 배수
PHSCFZ  :  6.0  :        : CFZ SCV 추가 내약 배수
PHIVAN  : 100.0 :        : VAN 세포내 추가 내약 배수
PHIRIF  : 200.0 :        : RIF 세포내 추가 내약 배수 (세포내는 비증식·SCV형)
PHILVX  :  5.0  :        : LVX 세포내 추가 내약 배수
PHIDAP  : 50.0  :        : DAP 세포내 추가 내약 배수 (파고솜에서 불활성)
PHILZD  :  3.0  :        : LZD 세포내 추가 내약 배수
PHICFZ  : 60.0  :        : CFZ 세포내 추가 내약 배수
CARVAN  :  0.20 :        : VAN 세포내 축적비
CARRIF  :  6.00 :        : RIF 세포내 축적비
CARLVX  :  5.00 :        : LVX 세포내 축적비
CARDAP  :  0.50 :        : DAP 세포내 축적비
CARLZD  :  1.00 :        : LZD 세포내 축적비
CARCFZ  :  0.10 :        : CFZ 세포내 축적비

// ------------------------------------------ BACTERIAL POPULATION ----------
MUP     :  0.35 : 1/h    : 유주균 최대 증식속도
MUB     :  0.080: 1/h    : 바이오필름 증식속도
MUS     :  0.050: 1/h    : SCV 증식속도
MUI     :  0.005: 1/h    : 세포내 균 증식속도
NMAX    :  1e10 : CFU    : 병소 전체 최대 수용력
NBCAP   :  3e9  : CFU    : 임플란트 표면 최대 정착 수용력
NPCAP   :  1e9  : CFU    : 비부착(관절액·농양) 균 최대 수용력
NICCAP  :  2e8  : CFU    : 세포내 저장고 최대 수용력 (숙주세포 수 제한)
KDNAT   :  0.005: 1/h    : 자연 사멸속도
KATT    :  0.40 : 1/h    : 부착속도 계수 (유주균 -> 바이오필름)
KDISP   :  0.020: 1/h    : 분산속도 계수 (agr 의존)
KPER0   :  0.0008:1/h    : 기저 지속균 형성속도
APER    :  6.0  :        : 성숙·스트레스에 의한 지속균 형성 증폭
KWAKE   :  0.0040:1/h    : 지속균 재활성(wake-up) 속도
KSCV0   :  0.00050:1/h   : SCV 형성속도
ASCV    :  8.0  :        : 국소 항생제 노출에 의한 SCV 유도 증폭
KREVS   :  0.0040:1/h    : SCV 정상형 복귀속도
KINT0   :  0.0060:1/h    : 골모세포 내재화 속도
KRELIC  :  0.0100:1/h    : 세포내 균 방출속도 (숙주세포 사멸 포함)
KEXT    :  2.0  : 1/h    : 1 CFU 미만 소멸 강제 속도 (클론 단위 적용)
NEXT    :  1.0  : CFU    : 소멸 문턱 (1 세포)

// ------------------------------------------------- BIOFILM MATRIX ---------
KSEPS   :  0.0035: 1/h   : EPS 생성속도 (성숙 시상수 ~10일)
KDEPS   :  0.0002: 1/h   : EPS 소실속도 (유지된 임플란트의 기질은 오래 남음)
KEPSN   :  3e7  : CFU    : EPS 생성 반포화 표면균량
KMAT    :  0.35 :        : EPS 성숙 반포화 상수
KSAIP   :  0.30 : 1/h    : AIP 생성속도
KDAIP   :  0.20 : 1/h    : AIP 소실속도
KAIPN   :  1e8  : CFU    : AIP 생성 반포화 균량
KAIP    :  0.40 :        : agr 활성 반포화 AIP

// --------------------------------------------- RESISTANCE (mutant supply) --
MURPOB  :  1e-8 :        : rpoB 변이 빈도 (분열당)
MUGYRA  :  3e-9 :        : gyrA 변이 빈도 (분열당)
FITRPOB :  0.15 :        : rpoB 변이체 적합도 비용
FITGYRA :  0.10 :        : gyrA 변이체 적합도 비용

// --------------------------------------------------- INNATE IMMUNITY ------
KIMM    :  5e6  : CFU/h  : 최대 식균 제거 용량 (완전 동원 시)
KMIMM   :  1e6  : CFU    : 식균 제거 포화 상수 (용량 공유)
FBFACT  :  0.12 :        : 이물 관련 국소 면역결핍 계수
FRUST   :  0.030:        : 좌절 식균작용 잔여 효율 (바이오필름)
FICIMM  :  0.020:        : 세포내 균에 대한 면역 효율
KPMN    :  500  : cells/uL : 면역 동원 반포화 PMN 밀도
PMN0    : 50    : cells/uL : 정상 관절액 호중구 밀도
MONO0   : 150   : cells/uL : 정상 관절액 단핵구 밀도
KINPMN  :  7.5  : cells/uL/h : PMN 유입 기저속도
KOUTPMN :  0.15 : 1/h    : PMN 소실속도
EMXCHEM : 1200  :        : 주화성 최대 동원 배수
KCHEM   :  0.35 :        : 주화성 반포화 신호
KTOXP   :  0.35 : 1/h    : 백혈구독소에 의한 PMN 사멸 최대속도
KINMON  :  1.5  : cells/uL/h : 단핵구 유입 기저속도
KOUTMON :  0.010: 1/h    : 단핵구 소실속도
EMXMON  : 18    :        : 단핵구 동원 최대 배수
KBACS   :  1e7  : CFU    : 세균 PAMP 신호 반포화 균량
KSMDSC  :  0.030: 1/h    : MDSC 동원 속도 (바이오필름 특이)
KDMDSC  :  0.020: 1/h    : MDSC 소실속도
KEPSM   :  0.30 :        : MDSC 동원 반포화 EPS
EMDSC   :  0.70 :        : MDSC 최대 면역억제 분율
KMD     :  0.50 :        : MDSC 억제 반포화
KSMAC   :  0.10 : 1/h    : 대식세포 활성화 속도
KDMAC   :  0.050: 1/h    : 대식세포 비활성화 속도

// ------------------------------------------------------- CYTOKINES --------
IL1B0   :  2.0  : pg/mL  : IL-1beta 기저
TNFA0   :  5.0  : pg/mL  : TNF-alpha 기저
IL60    :  3.0  : pg/mL  : IL-6 기저
IL100   :  4.0  : pg/mL  : IL-10 기저
KDIL1   :  0.70 : 1/h    : IL-1beta 소실
KDTNF   :  0.50 : 1/h    : TNF-alpha 소실
KDIL6   :  0.25 : 1/h    : IL-6 소실
KDIL10  :  0.30 : 1/h    : IL-10 소실
EIL1    : 60    :        : IL-1beta 최대 유도 배수
ETNF    : 40    :        : TNF-alpha 최대 유도 배수
EIL6    : 120   :        : IL-6 최대 유도 배수
EIL10   : 30    :        : IL-10 최대 유도 배수

// -------------------------------------------------- BIOMARKERS ------------
CRP0    :  3.0  : mg/L   : CRP 기저
KDCRP   :  0.0365:1/h    : CRP 소실 (t1/2 19 h)
ECRP    : 40    :        : CRP 최대 유도 배수
KCRP    : 60    : pg/mL  : CRP 유도 반포화 IL-6
ESR0    : 10    : mm/h   : ESR 기저
KDESR   :  0.00413:1/h   : ESR 소실 (t1/2 7일)
EESR    :  8.0  :        : ESR 최대 유도 배수
ADEF0   :  0.30 : S/CO   : alpha-defensin 기저
KDADEF  :  0.10 : 1/h    : alpha-defensin 소실
EADEF   : 12    :        : alpha-defensin 최대 유도 배수
KADEF   : 4000  : cells/uL : alpha-defensin 유도 반포화 PMN

// ----------------------------------------------------- BONE / OSTEOLYSIS --
KSRANK  :  0.10 : 1/h    : RANKL 생성 기저
KDRANK  :  0.10 : 1/h    : RANKL 소실
ERANK   :  9.0  :        : 사이토카인에 의한 RANKL 최대 유도
KSOPG   :  0.10 : 1/h    : OPG 생성 기저
KDOPG   :  0.10 : 1/h    : OPG 소실
KINOBL  :  0.40 : 1/h    : 골모세포 유입
KOUTOBL :  0.0040:1/h    : 골모세포 소실
KAPOP   :  0.0040:1/h    : 골모세포 아포토시스 최대속도
KINOCL  :  1.20 : 1/h    : 파골세포 유입 최대
KRR     :  1.00 :        : RANKL/OPG 반포화
KOUTOCL :  0.0060:1/h    : 파골세포 소실
KFORM   :  0.020: %/h    : 골 형성 속도 계수
KRES    :  0.020: %/h    : 골 흡수 속도 계수
KHOM    :  0.00050:1/h   : 골량 항상성 복귀
KLOOSE  :  0.0030:1/h    : 계면 이완 진행 계수
KREPAIR :  0.00050:1/h   : 계면 회복 계수
KSSEQ   :  0.0020:1/h    : 부골(sequestrum) 형성 계수
KCSEQ   :  0.00030:1/h   : 부골 자연 소실
KSEQA   :  0.30 :        : 부골에 의한 약물 접근 차단 반포화
SEQMAX  :  3.0  :        : 부골 최대 축적량

// ---------------------------------------------------------- TOXICITY ------
SCR0    :  0.90 : mg/dL  : 혈청 크레아티닌 기저
KELCR   :  0.030: 1/h    : 크레아티닌 소실
ENEPH   :  0.55 :        : VAN 최대 신독성 분율
AUCTHR  : 600   : mg*h/L : VAN AUC24 신독성 문턱
KNEPH   : 250   : mg*h/L : VAN 신독성 반포화
PLT0    : 250   : 1e9/L  : 혈소판 기저
KOUTPLT :  0.00417:1/h   : 혈소판 소실 (수명 10일)
EPLT    :  0.60 :        : LZD 최대 혈소판 생성 억제
KPLT    :  9.0  : mg/L   : LZD 혈소판독성 반포화 농도
ALT0    : 25    : U/L    : ALT 기저
KDALT   :  0.0144: 1/h   : ALT 소실 (t1/2 48 h)
EALT    :  4.0  :        : RIF 최대 ALT 유도
KALT    :  6.0  : mg/L   : RIF 간독성 반포화 농도

// ------------------------------------------------ SURGERY / STRATEGY ------
TSURG1  : 1e9   : h      : 1차 수술 시각 (기본: 시행 안 함)
LOGK1   :  0.0  : log10  : 1차 수술 균 제거량 (DAIR 2.5 / 1단계 4.5 / 2단계 6)
EXCH1   :  0    :        : 1차 수술 임플란트 교환 여부 (1 = 교환, EPS 초기화)
TSURG2  : 1e9   : h      : 2차 수술 시각 (2단계 재치환)
LOGK2   :  0.0  : log10  : 2차 수술 균 제거량
EXCH2   :  0    :        : 2차 수술 임플란트 교환 여부
TSDUR   :  1.0  : h      : 수술 제거 적용 구간 길이
SPCF0   :  0.0  : mg     : 스페이서 빠른 방출 풀 초기량
SPCS0   :  0.0  : mg     : 스페이서 느린 방출 풀 초기량

// -------------------------------------------------------- ENDPOINTS -------
PSEED   :  0.35 :        : 생존 1 CFU 당 재발 확률 (임플란트 유지 시)
PSEEDX  :  0.10 :        : 생존 1 CFU 당 재발 확률 (임플란트 교환 시)

$CMT @annotated
VANC   : VAN 중심구획 (mg)
VANP   : VAN 말초구획 (mg)
VANB   : VAN 골 유리농도 (mg/L)
RIFG   : RIF 위장관 흡수구획 (mg)
RIFC   : RIF 중심구획 (mg)
RIFB   : RIF 골 유리농도 (mg/L)
ENZR   : CYP3A4 유도 상태 (배수)
LVXG   : LVX 위장관 흡수구획 (mg)
LVXC   : LVX 중심구획 (mg)
LVXB   : LVX 골 유리농도 (mg/L)
DAPC   : DAP 중심구획 (mg)
DAPB   : DAP 골 유리농도 (mg/L)
LZDG   : LZD 위장관 흡수구획 (mg)
LZDC   : LZD 중심구획 (mg)
LZDB   : LZD 골 유리농도 (mg/L)
CFZC   : CFZ 중심구획 (mg)
CFZB   : CFZ 골 유리농도 (mg/L)
SPCF   : 스페이서 빠른 방출 풀 (mg)
SPCS   : 스페이서 느린 방출 풀 (mg)
LOCJ   : 관절액 국소 항생제 농도 (mg/L)
NP     : 유주균 (CFU)
NB     : 바이오필름 정착균 (CFU)
NPER   : 지속균 persisters (CFU)
NSCV   : 소집락변이체 SCV (CFU)
NIC    : 세포내 균 (CFU)
RP     : rpoB 변이 유주균 (CFU)
RB     : rpoB 변이 바이오필름균 (CFU)
QP     : gyrA 변이 유주균 (CFU)
QB     : gyrA 변이 바이오필름균 (CFU)
EPS    : 바이오필름 기질 성숙도 (0-1)
AIP    : agr 자가유도 펩티드 (정규화)
PMN    : 관절액 호중구 (cells/uL)
MONO   : 관절액 단핵구 (cells/uL)
MDSC   : 골수유래 억제세포 (정규화)
MAC    : 활성 대식세포 (정규화)
IL1B   : IL-1beta (pg/mL)
TNFA   : TNF-alpha (pg/mL)
IL6    : IL-6 (pg/mL)
IL10   : IL-10 (pg/mL)
CRP    : 혈청 CRP (mg/L)
ESR    : 혈청 ESR (mm/h)
ADEF   : 관절액 alpha-defensin (S/CO)
RANKL  : RANKL (정규화)
OPG    : OPG (정규화)
OBL    : 골모세포 (정규화 100)
OCL    : 파골세포 (정규화 100)
BONEV  : 임플란트 주위 골량 (% 기저)
LOOSEN : 계면 이완 지수 (0-100)
SEQ    : 부골 / 무혈관 격막 (정규화)
SCR    : 혈청 크레아티닌 (mg/dL)
PLT    : 혈소판 (1e9/L)
ALT    : 혈청 ALT (U/L)
AUCF   : VAN AUC24 이동추정 (mg*h/L)

$MAIN
VANB_0   = 0.0;   RIFB_0 = 0.0;  LVXB_0 = 0.0;
DAPB_0   = 0.0;   LZDB_0 = 0.0;  CFZB_0 = 0.0;
ENZR_0   = 1.0;
SPCF_0   = SPCF0;
SPCS_0   = SPCS0;
LOCJ_0   = 0.0;

// Inoculum. NOTE: this is driven by the PARAMETER INOC0, not by init().
// $MAIN runs on every record, so writing a literal here would silently
// override anything passed to mrgsolve::init() -- change INOC0 instead.
NP_0     = INOC0;
NB_0     = 0.0;
NPER_0   = 0.0;
NSCV_0   = 0.0;
NIC_0    = 0.0;
RP_0     = 0.0;  RB_0 = 0.0;  QP_0 = 0.0;  QB_0 = 0.0;

EPS_0    = 0.0;
AIP_0    = 0.0;
PMN_0    = PMN0;
MONO_0   = MONO0;
MDSC_0   = 0.02;
MAC_0    = 0.10;
IL1B_0   = IL1B0;
TNFA_0   = TNFA0;
IL6_0    = IL60;
IL10_0   = IL100;
CRP_0    = CRP0;
ESR_0    = ESR0;
ADEF_0   = ADEF0;
RANKL_0  = 1.0;
OPG_0    = 1.0;
OBL_0    = 100.0;
OCL_0    = 100.0;
BONEV_0  = 100.0;
LOOSEN_0 = 0.0;
SEQ_0    = 0.0;
SCR_0    = SCR0;
PLT_0    = PLT0;
ALT_0    = ALT0;
AUCF_0   = 0.0;

$ODE
// Guard: the stiff 10-log bacterial collapse can leave states a hair below
// zero, and pow(negative, 3) then poisons the extinction term.
if(NP < 0)   NP   = 0.0;   if(NB < 0)   NB   = 0.0;
if(NPER < 0) NPER = 0.0;   if(NSCV < 0) NSCV = 0.0;
if(NIC < 0)  NIC  = 0.0;   if(RP < 0)   RP   = 0.0;
if(RB < 0)   RB   = 0.0;   if(QP < 0)   QP   = 0.0;
if(QB < 0)   QB   = 0.0;   if(EPS < 0)  EPS  = 0.0;
if(AIP < 0)  AIP  = 0.0;   if(SEQ < 0)  SEQ  = 0.0;

// =====================================================================
// 0. Surgical burden-reduction pulses (Pillar 3)
//    A pulse of height ln(10)*LOGK/TSDUR applied for TSDUR hours multiplies
//    every bacterial state by exactly 10^-LOGK. No event plumbing needed and
//    the arithmetic is exact.
// =====================================================================
double SR1 = 0.0, SR2 = 0.0, WIPE = 0.0;
if(SOLVERTIME >= TSURG1 && SOLVERTIME < (TSURG1 + TSDUR)) {
  SR1 = 2.302585093*LOGK1/TSDUR;
  if(EXCH1 > 0.5) WIPE = 1.0;
}
if(SOLVERTIME >= TSURG2 && SOLVERTIME < (TSURG2 + TSDUR)) {
  SR2 = 2.302585093*LOGK2/TSDUR;
  if(EXCH2 > 0.5) WIPE = 1.0;
}
double SURG = SR1 + SR2;
// Debridement also strips matrix; an implant exchange strips it completely.
// A RETAINED implant keeps its adherent matrix: debridement removes pus and
// necrotic soft tissue, not the biofilm bonded to metal. Only an exchange
// wipes the matrix. This asymmetry is what makes DAIR a race against EPS.
double EPSCUT = SURG*(WIPE > 0.5 ? 1.0 : 0.05);
double SEQCUT = SURG*(WIPE > 0.5 ? 1.0 : 0.10);

// =====================================================================
// 1. Pharmacokinetics
// =====================================================================
double CVAN = VANC/VCVAN;
double CVANP= VANP/VPVAN;
double CRIF = RIFC/VCRIF;
double CLVX = LVXC/VCLVX;
double CDAP = DAPC/VCDAP;
double CLZD = LZDC/VCLZD;
double CCFZ = CFZC/VCCFZ;

// Rifampicin autoinduction and its induction of linezolid clearance
double CLRIF = CLRIF0*ENZR;
double CLLZD = CLLZD0*(1.0 + RIFLZD*(ENZR - 1.0));

dxdt_VANC = -CLVAN*CVAN - QVAN*CVAN + QVAN*CVANP;
dxdt_VANP =  QVAN*CVAN - QVAN*CVANP;
dxdt_RIFG = -KARIF*RIFG;
dxdt_RIFC =  KARIF*FRIF*RIFG - CLRIF*CRIF;
dxdt_ENZR =  KENZ*((1.0 + EIND*CRIF/(3.0 + CRIF)) - ENZR);
dxdt_LVXG = -KALVX*LVXG;
dxdt_LVXC =  KALVX*FLVX*LVXG - CLLVX*CLVX;
dxdt_DAPC = -CLDAP*CDAP;
dxdt_LZDG = -KALZD*LZDG;
dxdt_LZDC =  KALZD*FLZD*LZDG - CLLZD*CLZD;
dxdt_CFZC = -CLCFZ*CCFZ;

// Free BONE concentrations. FB = (AUC_bone/AUC_plasma) x fu  (Pillar 2)
double FBVAN = PENVAN*FUVAN;
double FBRIF = PENRIF*FURIF;
double FBLVX = PENLVX*FULVX;
double FBDAP = PENDAP*FUDAP;
double FBLZD = PENLZD*FULZD;
double FBCFZ = PENCFZ*FUCFZ;

dxdt_VANB = KEQVAN*(FBVAN*CVAN - VANB);
dxdt_RIFB = KEQRIF*(FBRIF*CRIF - RIFB);
dxdt_LVXB = KEQLVX*(FBLVX*CLVX - LVXB);
dxdt_DAPB = KEQDAP*(FBDAP*CDAP - DAPB);
dxdt_LZDB = KEQLZD*(FBLZD*CLZD - LZDB);
dxdt_CFZB = KEQCFZ*(FBCFZ*CCFZ - CFZB);

// Local elution from antibiotic-loaded cement (two-pool, burst + matrix)
dxdt_SPCF = -KFAST*SPCF;
dxdt_SPCS = -KSLOW*SPCS;
dxdt_LOCJ = (KFAST*SPCF + KSLOW*SPCS)/VJOINT - KCLJ*LOCJ;

// Rolling 24 h vancomycin AUC estimate (for nephrotoxicity)
dxdt_AUCF = (24.0*CVAN - AUCF)/24.0;

// =====================================================================
// 2. Biofilm state variables and drug access
// =====================================================================
double EPSM  = EPS/(EPS + KMAT);                    // maturity 0-1
double AGRA  = AIP/(AIP + KAIP);                    // agr activity 0-1
double ACCESS= 1.0/(1.0 + SEQ/KSEQA);               // sequestrum shielding

// Effective drug concentrations seen by each niche
double VANEP = (VANB + FLOCP*LOCJ);
double VANEB = (VANB + FLOCB*LOCJ)*ACCESS;
double RIFEB = RIFB*ACCESS;
double LVXEB = LVXB*ACCESS;
double DAPEB = DAPB*ACCESS;
double LZDEB = LZDB*ACCESS;
double CFZEB = CFZB*ACCESS;

// Maturity-scaled biofilm tolerance multipliers (Pillar 4)
double PBVAN = PHIMAT(MBCVAN/MICVAN, EPSM);
double PBRIF = PHIMAT(MBCRIF/MICRIF, EPSM);
double PBLVX = PHIMAT(MBCLVX/MICLVX, EPSM);
double PBDAP = PHIMAT(MBCDAP/MICDAP, EPSM);
double PBLZD = PHIMAT(MBCLZD/MICLZD, EPSM);
double PBCFZ = PHIMAT(MBCCFZ/MICCFZ, EPSM);

double MS = (MSSA > 0.5) ? 1.0 : 0.0;

// =====================================================================
// 3. Antibiotic kill by bacterial state
//    Every term is EMAX * C / (EC50_state + C) with
//    EC50_state = MIC x (biofilm PHI) x (state PHI).
// =====================================================================
// --- planktonic (rifampicin-susceptible, FQ-susceptible) ---
double KLP = EKILL(EMXVAN, VANEP, MICVAN)
           + EKILL(EMXRIF, RIFB , MICRIF)
           + EKILL(EMXLVX, LVXB , MICLVX)
           + EKILL(EMXDAP, DAPB , MICDAP)
           + EKILL(EMXLZD, LZDB , MICLZD)
           + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double INHP = EMXLZDI*LZDB/(MICLZD + LZDB);

// --- biofilm-embedded ---
double KLB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
           + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF)
           + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX)
           + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
           + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
           + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);
double INHB = EMXLZDI*LZDEB/(MICLZD*PBLZD + LZDEB);

// --- persisters ---
double KLPER = EKILL(EMXVAN, VANEB, MICVAN*PBVAN*PHPVAN)
             + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF*PHPRIF)
             + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX*PHPLVX)
             + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP*PHPDAP)
             + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD*PHPLZD)
             + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ*PHPCFZ);

// --- small colony variants ---
double KLSCV = EKILL(EMXVAN, VANEB, MICVAN*PBVAN*PHSVAN)
             + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF*PHSRIF)
             + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX*PHSLVX)
             + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP*PHSDAP)
             + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD*PHSLZD)
             + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ*PHSCFZ);

// --- intracellular: cellular accumulation ratio x free bone concentration,
//     shielded by the same avascular sequestrum (the canalicular reservoir
//     sits inside dead cortical bone) ---
double VANIC = CARVAN*VANB*ACCESS;  double RIFIC = CARRIF*RIFB*ACCESS;
double LVXIC = CARLVX*LVXB*ACCESS;  double DAPIC = CARDAP*DAPB*ACCESS;
double LZDIC = CARLZD*LZDB*ACCESS;  double CFZIC = CARCFZ*CFZB*ACCESS;
double KLIC = EKILL(EMXVAN, VANIC, MICVAN*PBVAN*PHIVAN)
            + EKILL(EMXRIF, RIFIC, MICRIF*PBRIF*PHIRIF)
            + EKILL(EMXLVX, LVXIC, MICLVX*PBLVX*PHILVX)
            + EKILL(EMXDAP, DAPIC, MICDAP*PBDAP*PHIDAP)
            + EKILL(EMXLZD, LZDIC, MICLZD*PBLZD*PHILZD)
            + MS*EKILL(EMXCFZ, CFZIC, MICCFZ*PBCFZ*PHICFZ);

// --- rpoB mutants: identical EXCEPT rifampicin does nothing ---
double KLRP = EKILL(EMXVAN, VANEP, MICVAN)
            + EKILL(EMXLVX, LVXB , MICLVX)
            + EKILL(EMXDAP, DAPB , MICDAP)
            + EKILL(EMXLZD, LZDB , MICLZD)
            + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double KLRB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
            + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX)
            + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
            + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
            + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);

// --- gyrA mutants: identical EXCEPT levofloxacin does nothing ---
double KLQP = EKILL(EMXVAN, VANEP, MICVAN)
            + EKILL(EMXRIF, RIFB , MICRIF)
            + EKILL(EMXDAP, DAPB , MICDAP)
            + EKILL(EMXLZD, LZDB , MICLZD)
            + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double KLQB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
            + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF)
            + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
            + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
            + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);

// =====================================================================
// 4. Innate immunity (Pillar 1)
// =====================================================================
double NTOT  = NP + NB + NPER + NSCV + NIC + RP + RB + QP + QB;
double NSURF = NB + NPER + RB + QB;                  // surface-attached pool
double CROWD = 1.0 - NTOT/NMAX;  if(CROWD < 0.0) CROWD = 0.0;

double PMNF  = PMN/(PMN + KPMN);
double MDSUP = 1.0 - EMDSC*MDSC/(KMD + MDSC);
double FB    = (IMPL > 0.5) ? FBFACT : 1.0;
double IMMCAP= KIMM*PMNF*FB*IMSUP*MDSUP;

// Phagocyte capacity is a SHARED resource, allocated across the accessible
// burden in proportion to how reachable each state is. Writing one
// Michaelis-Menten sink per compartment (the obvious but wrong way) makes a
// rare mutant clone face the entire immune system as if it were alone, which
// silently sterilises exactly the subpopulation the model exists to track.
double NPHAG = (NP + RP + QP)
             + FRUST*(NB + RB + QB + NSCV)
             + FICIMM*NIC;
double IMMTOT= IMMCAP*NPHAG/(KMIMM + NPHAG);
double SHR   = IMMTOT/(NPHAG + 1e-30);
double IMMP  = SHR*NP;
double IMMB  = SHR*FRUST*NB;                          // frustrated phagocytosis
double IMMS  = SHR*FRUST*NSCV;
double IMMI  = SHR*FICIMM*NIC;
double IMMRP = SHR*RP;
double IMMRB = SHR*FRUST*RB;
double IMMQP = SHR*QP;
double IMMQB = SHR*FRUST*QB;

// =====================================================================
// 5. Bacterial state transitions
// =====================================================================
double SFREE = (IMPL > 0.5) ? (1.0 - NSURF/NBCAP) : 0.0;
if(SFREE < 0.0) SFREE = 0.0;
double ATT   = KATT*SFREE;
double DISP  = KDISP*AGRA;
double STRESS= KLB/(KLB + 0.05);
double KPER  = KPER0*(1.0 + APER*(EPSM + STRESS));
double KSCV  = KSCV0*(1.0 + ASCV*LOCJ/(LOCJ + 20.0));
double KINT  = KINT0*(OBL/100.0);

// Smooth sub-single-cell extinction. Applied PER CLONE, not per phenotypic
// state: a cell can move between planktonic / sessile / persister / SCV /
// intracellular, so "fewer than one organism left" is a property of the
// lineage. Applying it per compartment annihilates the rare rpoB clone the
// instant it dips below 1 CFU in any one state and destroys Pillar 3.
double NSUSC = NP + NB + NPER + NSCV + NIC;
double NRPOB = RP + RB;
double NGYRA = QP + QB;
double EXS = KEXT/(1.0 + pow(NSUSC/NEXT, 3.0));
double EXR = KEXT/(1.0 + pow(NRPOB/NEXT, 3.0));
double EXQ = KEXT/(1.0 + pow(NGYRA/NEXT, 3.0));

// Carrying capacities: the free (joint fluid / abscess) pool and the
// intracellular reservoir are each bounded by their own niche, not only by
// the lesion total. Without this the dispersed planktonic pool overruns the
// biofilm and the infection becomes spuriously vancomycin-curable.
double NFREE = NP + RP + QP;
double CROWDP= CROWD*(1.0 - NFREE/NPCAP);  if(CROWDP < 0.0) CROWDP = 0.0;
double CROWDI= CROWD*(1.0 - NIC/NICCAP);   if(CROWDI < 0.0) CROWDI = 0.0;

// Growth fluxes (needed for mutant supply: mutants arise at DIVISION)
double GFP = MUP*(1.0 - INHP)*CROWDP*NP;
double GFB = MUB*(1.0 - INHB)*CROWD*NB;
double GFS = MUS*CROWD*NSCV;
double GFI = MUI*CROWDI*NIC;
double MUTR = MURPOB*(GFP + GFB);      // rpoB mutant supply  (Pillar 3)
double MUTQ = MUGYRA*(GFP + GFB);      // gyrA mutant supply

dxdt_NP   = GFP - MUTR - MUTQ
          - (KDNAT + KLP + EXS + SURG)*NP
          - IMMP - ATT*NP + DISP*NB - KINT*NP + KRELIC*NIC;

dxdt_NB   = GFB
          + ATT*NP - DISP*NB
          - (KDNAT + KLB + EXS + SURG)*NB
          - IMMB - KPER*NB + KWAKE*NPER - KSCV*NB + KREVS*NSCV;

dxdt_NPER = KPER*NB - KWAKE*NPER - (KLPER + EXS + SURG)*NPER;

dxdt_NSCV = GFS + KSCV*NB - KREVS*NSCV
          - (KDNAT + KLSCV + EXS + SURG)*NSCV - IMMS;

dxdt_NIC  = GFI + KINT*NP - KRELIC*NIC
          - (KLIC + EXS + 0.5*SURG)*NIC - IMMI;

dxdt_RP   = MUTR + MUP*(1.0 - FITRPOB)*(1.0 - INHP)*CROWDP*RP
          - (KDNAT + KLRP + EXR + SURG)*RP - IMMRP - ATT*RP + DISP*RB;

dxdt_RB   = MUB*(1.0 - FITRPOB)*(1.0 - INHB)*CROWD*RB
          + ATT*RP - DISP*RB
          - (KDNAT + KLRB + EXR + SURG)*RB - IMMRB;

dxdt_QP   = MUTQ + MUP*(1.0 - FITGYRA)*(1.0 - INHP)*CROWDP*QP
          - (KDNAT + KLQP + EXQ + SURG)*QP - IMMQP - ATT*QP + DISP*QB;

dxdt_QB   = MUB*(1.0 - FITGYRA)*(1.0 - INHB)*CROWD*QB
          + ATT*QP - DISP*QB
          - (KDNAT + KLQB + EXQ + SURG)*QB - IMMQB;

// =====================================================================
// 6. Biofilm matrix and quorum sensing
// =====================================================================
dxdt_EPS = KSEPS*(NSURF/(NSURF + KEPSN))*(1.0 - EPS) - KDEPS*EPS - EPSCUT*EPS;
dxdt_AIP = KSAIP*(NP + NB)/((NP + NB) + KAIPN) - KDAIP*AIP;

// =====================================================================
// 7. Innate immune cells, MDSC and cytokines
// =====================================================================
double BACSIG = NTOT/(NTOT + KBACS);
double ACTSIG = BACSIG*MDSUP;
double TOXL   = AGRA*BACSIG;

dxdt_PMN  = KINPMN*(1.0 + EMXCHEM*ACTSIG/(KCHEM + ACTSIG))
          - KOUTPMN*PMN - KTOXP*TOXL*PMN;
dxdt_MONO = KINMON*(1.0 + EMXMON*BACSIG/(KCHEM + BACSIG))
          - KOUTMON*MONO;
dxdt_MDSC = KSMDSC*EPS/(EPS + KEPSM) - KDMDSC*MDSC;
dxdt_MAC  = KSMAC*BACSIG - KDMAC*MAC;

dxdt_IL1B = KDIL1*IL1B0*(1.0 + EIL1*ACTSIG)  - KDIL1*IL1B;
dxdt_TNFA = KDTNF*TNFA0*(1.0 + ETNF*ACTSIG)  - KDTNF*TNFA;
dxdt_IL6  = KDIL6*IL60 *(1.0 + EIL6*ACTSIG)  - KDIL6*IL6;
dxdt_IL10 = KDIL10*IL100*(1.0 + EIL10*MDSC/(KMD + MDSC)) - KDIL10*IL10;

// =====================================================================
// 8. Systemic biomarkers (2018 ICM criteria read-outs)
// =====================================================================
dxdt_CRP  = KDCRP*CRP0*(1.0 + ECRP*(IL6 - IL60)/(KCRP + (IL6 - IL60))) - KDCRP*CRP;
dxdt_ESR  = KDESR*ESR0*(1.0 + EESR*(IL6 - IL60)/(KCRP + (IL6 - IL60))) - KDESR*ESR;
dxdt_ADEF = KDADEF*ADEF0*(1.0 + EADEF*(PMN - PMN0)/(KADEF + (PMN - PMN0))) - KDADEF*ADEF;

// =====================================================================
// 9. Bone: RANKL/OPG, osteoclasts, osteolysis, loosening, sequestrum
// =====================================================================
double CYTB = (TNFA - TNFA0)/(50.0 + (TNFA - TNFA0))
            + (IL1B - IL1B0)/(30.0 + (IL1B - IL1B0))
            + (IL6  - IL60 )/(80.0 + (IL6  - IL60 ));
if(CYTB < 0.0) CYTB = 0.0;
double APOS = 0.6*AGRA*BACSIG + NIC/(NIC + 1e6);

dxdt_RANKL = KSRANK*(1.0 + ERANK*CYTB/3.0) - KDRANK*RANKL;
dxdt_OPG   = KSOPG*(OBL/100.0)             - KDOPG*OPG;
dxdt_OBL   = KINOBL - KOUTOBL*OBL - KAPOP*APOS*OBL;
double RR  = RANKL/(OPG + 1e-6);
dxdt_OCL   = KINOCL*RR/(KRR + RR) - KOUTOCL*OCL;
dxdt_BONEV = KFORM*(OBL/100.0) - KRES*(OCL/100.0) + KHOM*(100.0 - BONEV);
double BLOSS = 100.0 - BONEV;  if(BLOSS < 0.0) BLOSS = 0.0;
dxdt_LOOSEN= KLOOSE*BLOSS*(1.0 - LOOSEN/100.0) - KREPAIR*LOOSEN;
dxdt_SEQ   = KSSEQ*(KRES*(OCL/100.0))*BACSIG*10.0*(1.0 - SEQ/SEQMAX)
           - KCSEQ*SEQ - SEQCUT*SEQ;

// =====================================================================
// 10. Drug toxicity
// =====================================================================
double EXAUC = AUCF - AUCTHR;  if(EXAUC < 0.0) EXAUC = 0.0;
double NEPH  = ENEPH*EXAUC/(KNEPH + EXAUC);
dxdt_SCR  = KELCR*SCR0 - KELCR*(1.0 - NEPH)*SCR;
dxdt_PLT  = KOUTPLT*PLT0*(1.0 - EPLT*CLZD/(KPLT + CLZD)) - KOUTPLT*PLT;
dxdt_ALT  = KDALT*ALT0*(1.0 + EALT*CRIF/(KALT + CRIF)) - KDALT*ALT;

$TABLE
double NTOTAL = NP + NB + NPER + NSCV + NIC + RP + RB + QP + QB;
if(NTOTAL < 0.0) NTOTAL = 0.0;
double NSUSC  = NP + NB + NPER + NSCV + NIC;
double NRES   = RP + RB + QP + QB;
capture LOGNTOT = log10(NTOTAL + 1e-12);
capture LOGNP   = log10(NP  + 1e-12);
capture LOGNB   = log10(NB  + 1e-12);
capture LOGNPER = log10(NPER+ 1e-12);
capture LOGNSCV = log10(NSCV+ 1e-12);
capture LOGNIC  = log10(NIC + 1e-12);
capture LOGNRES = log10(NRES+ 1e-12);
capture RESFRAC = NRES/(NTOTAL + 1e-12);

// --- Pillar 2 read-out: achievable free bone concentration over MBEC -------
capture CPVAN = VANC/VCVAN;
capture CPRIF = RIFC/VCRIF;
capture CPLVX = LVXC/VCLVX;
capture CPDAP = DAPC/VCDAP;
capture CPLZD = LZDC/VCLZD;
capture CPCFZ = CFZC/VCCFZ;
capture RATVAN = (VANB + FLOCB*LOCJ)/MBCVAN;
capture RATRIF = RIFB/MBCRIF;
capture RATLVX = LVXB/MBCLVX;
capture RATDAP = DAPB/MBCDAP;
capture RATLZD = LZDB/MBCLZD;
capture RATCFZ = CFZB/MBCCFZ;
capture RATMAX = RATVAN;
if(RATRIF > RATMAX) RATMAX = RATRIF;
if(RATLVX > RATMAX) RATMAX = RATLVX;
if(RATDAP > RATMAX) RATMAX = RATDAP;
if(RATLZD > RATMAX) RATMAX = RATLZD;
if(MSSA > 0.5 && RATCFZ > RATMAX) RATMAX = RATCFZ;

// --- Pillar 3 read-out: mutant supply --------------------------------------
capture PRPOB = 1.0 - exp(-MURPOB*NTOTAL);   // P(pre-existing rpoB mutant)
capture PGYRA = 1.0 - exp(-MUGYRA*NTOTAL);
capture EMUTR = MURPOB*NTOTAL;               // expected rpoB mutants present

// --- Pillar 4 read-out: biofilm maturity ------------------------------------
capture EPSMAT = EPS/(EPS + KMAT);

// --- clinical / diagnostic read-outs ---------------------------------------
capture SYNWBC = PMN + MONO;
capture PMNPCT = 100.0*PMN/(PMN + MONO + 1e-9);
// 2018 ICM minor-criteria score (CRP 2, D-dimer 2, synovial WBC 3,
// PMN% 2, alpha-defensin 3, LE 3): >=6 infected, 2-5 inconclusive
double SC = 0.0;
if(CRP > 10.0)     SC += 2.0;
if(SYNWBC > 3000.0)SC += 3.0;
if(PMNPCT > 80.0)  SC += 2.0;
if(ADEF > 1.0)     SC += 3.0;
capture ICMSC = SC;

// Cure probability: PSEED is the chance that ONE surviving CFU regrows.
double PS = (EXCH1 > 0.5 || EXCH2 > 0.5) ? PSEEDX : PSEED;
capture PCURE = exp(-PS*NTOTAL);
capture PFAIL = 1.0 - exp(-PS*NTOTAL);

capture CRCL  = 100.0*(SCR0/SCR);           // crude creatinine-clearance proxy
capture AKIFL = (SCR > 1.5*SCR0) ? 1.0 : 0.0;
capture TCPFL = (PLT < 100.0)    ? 1.0 : 0.0;
capture HEPFL = (ALT > 3.0*ALT0) ? 1.0 : 0.0;
'

pji_mod <- mrgsolve::mcode_cache("pji", pji_code)
pji_mod <- mrgsolve::update(pji_mod, end = 4320, delta = 4,
                            atol = 1e-8, rtol = 1e-6, maxsteps = 500000)

## =============================================================================
##  DOSING BUILDERS
##  Standard adult regimens; amt in mg, time in hours.
## =============================================================================
DAY <- 24

ev_van  <- function(start = 0, weeks = 6, amt = 1000, ii = 12)
  ev(time = start, cmt = "VANC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_rif  <- function(start = 0, weeks = 12, amt = 450, ii = 12)
  ev(time = start, cmt = "RIFG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_lvx  <- function(start = 0, weeks = 12, amt = 750, ii = 24)
  ev(time = start, cmt = "LVXG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_dap  <- function(start = 0, weeks = 6, amt = 560, ii = 24)
  ev(time = start, cmt = "DAPC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_lzd  <- function(start = 0, weeks = 6, amt = 600, ii = 12)
  ev(time = start, cmt = "LZDG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_cfz  <- function(start = 0, weeks = 6, amt = 2000, ii = 8)
  ev(time = start, cmt = "CFZC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)

## =============================================================================
##  SCENARIOS
##  All start from a 100 CFU intra-operative contamination of a hip/knee
##  prosthesis at t = 0 and are simulated for 180 days (4320 h).
##  Diagnosis is assumed at day 21 (acute post-operative PJI) unless stated.
## =============================================================================
DX  <- 21*DAY     # diagnosis / surgery time for the acute presentations
DXC <- 90*DAY     # delayed (chronic) presentation

scen <- list(

  ## 1. Natural history, no treatment ---------------------------------------
  "01_무치료 자연경과" = list(
    par = list(),
    ev  = ev(time = 0, cmt = "NP", amt = 0)),

  ## 2. Antibiotics only, implant retained, NO surgery -----------------------
  ##    ("medical management" — the arm that shows Pillar 2 alone is not enough)
  "02_수술없이 반코마이신 단독" = list(
    par = list(),
    ev  = ev_van(DX, 6)),

  ## 3. DAIR + vancomycin monotherapy 6 weeks --------------------------------
  "03_DAIR + 반코마이신 6주" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_van(DX, 6)),

  ## 4. DAIR + vancomycin + rifampicin 12 weeks ------------------------------
  "04_DAIR + 반코+리팜피신 12주" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_van(DX, 6) + ev_rif(DX, 12)),

  ## 5. DAIR + levofloxacin + rifampicin (Zimmerli 1998 JAMA analogue) -------
  "05_DAIR + 레보플록사신+리팜피신 12주" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 6. DAIR + levofloxacin alone (the Zimmerli control arm) -----------------
  "06_DAIR + 레보플록사신 단독 12주" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12)),

  ## 7. DAIR + RIFAMPICIN MONOTHERAPY (Pillar 3 demonstration) ---------------
  "07_DAIR + 리팜피신 단독 12주" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_rif(DX, 12)),

  ## 8. NO debridement + rifampicin combination -----------------------------
  ##    (mature 1e10 biofilm: the rpoB mutant is already there)
  "08_변연절제 없이 리팜피신 병용" = list(
    par = list(),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 9. One-stage exchange + levofloxacin/rifampicin -------------------------
  "09_1단계 교환 + LVX/RIF 12주" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 10. Two-stage exchange with an antibiotic-loaded spacer -----------------
  "10_2단계 교환 + 스페이서" = list(
    par = list(TSURG1 = DX,        LOGK1 = 6.0, EXCH1 = 1,
               TSURG2 = DX + 8*7*DAY, LOGK2 = 1.5, EXCH2 = 1,
               SPCF0  = 1200, SPCS0 = 700),
    ev  = ev_van(DX, 6) + ev_rif(DX, 12)),

  ## 11. OVIVA analogue: fully oral vs IV backbone ---------------------------
  "11_OVIVA 경구(LVX+RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_lvx(DX, 6) + ev_rif(DX, 6)),
  "12_OVIVA 정주(VAN)+RIF" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_van(DX, 6) + ev_rif(DX, 6)),

  ## 13-14. DATIPO analogue: 6 vs 12 weeks -----------------------------------
  "13_DATIPO 6주 (DAIR+LVX/RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 6) + ev_rif(DX, 6)),
  "14_DATIPO 12주 (DAIR+LVX/RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 15-16. Late DAIR on a mature biofilm (Pillar 4) -------------------------
  ##    Read these against scenario 05 (same regimen, DAIR on day 21):
  ##    day  21 -> sterile on day 60 of an 84-day course (24 days of margin)
  ##    day  90 -> sterile on day 84 of 84 (zero margin)
  ##    day 180 -> never sterile; relapses after antibiotics stop
  "15_지연 진단(90일) + DAIR + LVX/RIF" = list(
    par = list(TSURG1 = DXC, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DXC, 12) + ev_rif(DXC, 12)),
  "16_만성(180일) + DAIR + LVX/RIF" = list(
    par = list(TSURG1 = 180*DAY, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(180*DAY, 12) + ev_rif(180*DAY, 12)),

  ## 16. Daptomycin + rifampicin --------------------------------------------
  "17_DAIR + 답토마이신+리팜피신" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_dap(DX, 6) + ev_rif(DX, 12)),

  ## 17. Linezolid + rifampicin (with the CYP/P-gp interaction) --------------
  "18_DAIR + 리네졸리드+리팜피신" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lzd(DX, 6) + ev_rif(DX, 12)),

  ## 18. MSSA: cefazolin + rifampicin ---------------------------------------
  "19_MSSA DAIR + 세파졸린+리팜피신" = list(
    par = list(MSSA = 1, TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_cfz(DX, 6) + ev_rif(DX, 12)),

  ## 19. Immunocompromised host (diabetes / anti-TNF) ------------------------
  "20_면역저하 숙주 DAIR + LVX/RIF" = list(
    par = list(IMSUP = 0.5, TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 20. Chronic oral suppression (inoperable patient) -----------------------
  "21_만성 경구 억제요법" = list(
    par = list(),
    ev  = ev_lvx(DX, 26))
)

run_scen <- function(mod, s, end = 4320) {
  p <- s$par
  m <- mod
  if(length(p)) m <- mrgsolve::param(m, p)
  mrgsolve::mrgsim(m, events = s$ev, end = end, delta = 4) %>% as.data.frame()
}

run_all <- function(mod = pji_mod, end = 4320) {
  out <- lapply(names(scen), function(nm) {
    d <- run_scen(mod, scen[[nm]], end = end)
    d$scenario <- nm
    d
  })
  dplyr::bind_rows(out)
}

## =============================================================================
##  PILLAR 1 CHECK — the foreign body moves the infective dose
##  Run a dose-escalation with and without an implant and find the inoculum at
##  which the model still has > 1e6 CFU at day 60.
## =============================================================================
id50_scan <- function(mod = pji_mod,
                      inoc = 10^seq(0, 9, by = 0.5)) {
  f <- function(impl) {
    sapply(inoc, function(n0) {
      m <- mrgsolve::param(mod, IMPL = impl, INOC0 = n0)
      o <- mrgsolve::mrgsim(m, end = 60*24, delta = 24) %>% as.data.frame()
      max(tail(o$LOGNTOT, 1), -12)
    })
  }
  data.frame(inoculum = inoc,
             log10_burden_d60_with_implant    = f(1),
             log10_burden_d60_without_implant = f(0))
}

## =============================================================================
##  PILLAR 2 CHECK — the C_bone,free / MBEC table at steady state
## =============================================================================
ratio_table <- function(mod = pji_mod) {
  p <- as.list(mrgsolve::param(mod))
  auc24 <- c(
    VAN = 2*1000/p$CLVAN,
    RIF = 2*450*p$FRIF/p$CLRIF0/2,       # /2 approximates the autoinduced state
    LVX = 750*p$FLVX/p$CLLVX,
    DAP = 560/p$CLDAP,
    LZD = 2*600*p$FLZD/p$CLLZD0,
    CFZ = 3*2000/p$CLCFZ)
  cavg  <- auc24/24
  fbone <- c(VAN = p$PENVAN*p$FUVAN, RIF = p$PENRIF*p$FURIF,
             LVX = p$PENLVX*p$FULVX, DAP = p$PENDAP*p$FUDAP,
             LZD = p$PENLZD*p$FULZD, CFZ = p$PENCFZ*p$FUCFZ)
  mbec  <- c(VAN = p$MBCVAN, RIF = p$MBCRIF, LVX = p$MBCLVX,
             DAP = p$MBCDAP, LZD = p$MBCLZD, CFZ = p$MBCCFZ)
  mic   <- c(VAN = p$MICVAN, RIF = p$MICRIF, LVX = p$MICLVX,
             DAP = p$MICDAP, LZD = p$MICLZD, CFZ = p$MICCFZ)
  data.frame(drug = names(auc24),
             AUC24_plasma = round(auc24, 1),
             Cavg_plasma  = round(cavg, 2),
             f_bone_free  = round(fbone, 3),
             Cbone_free   = round(cavg*fbone, 3),
             MIC          = mic,
             MBEC         = mbec,
             ratio_MIC    = round(cavg*fbone/mic, 1),
             ratio_MBEC   = round(cavg*fbone/mbec, 4),
             row.names = NULL)
}

## =============================================================================
##  PILLAR 3 CHECK — mutant supply before and after debridement
## =============================================================================
mutant_supply <- function(N = 10^c(4:10), mu = 1e-8) {
  data.frame(burden_CFU        = N,
             expected_mutants  = signif(mu*N, 3),
             P_preexisting     = signif(1 - exp(-mu*N), 3))
}

## =============================================================================
##  SUMMARY TABLE
## =============================================================================
summarise_scen <- function(df) {
  df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      burden_end_log10 = round(dplyr::last(LOGNTOT), 2),
      resistant_log10  = round(dplyr::last(LOGNRES), 2),
      P_cure           = round(dplyr::last(PCURE), 3),
      CRP_end          = round(dplyr::last(CRP), 1),
      CRP_peak         = round(max(CRP), 1),
      synWBC_peak      = round(max(SYNWBC), 0),
      bone_volume_end  = round(dplyr::last(BONEV), 1),
      loosening_end    = round(dplyr::last(LOOSEN), 1),
      max_bone_MBEC_ratio = round(max(RATMAX), 3),
      # Mutant supply at the LOWEST burden the regimen ever reaches: this is
      # the post-debridement number that decides whether rifampicin is usable.
      P_rpoB_at_nadir  = round(1 - exp(-1e-8*min(10^LOGNTOT)), 4),
      AKI              = max(AKIFL),
      thrombocytopenia = max(TCPFL),
      hepatotoxicity   = max(HEPFL),
      .groups = "drop")
}

## =============================================================================
##  EXAMPLE USAGE
## -----------------------------------------------------------------------------
##  out <- run_all()
##  summarise_scen(out)
##  ratio_table()          # Pillar 2 — the table that explains rifampicin
##  mutant_supply()        # Pillar 3 — why debridement precedes rifampicin
##  id50_scan()            # Pillar 1 — the foreign-body ID50 shift
## =============================================================================
