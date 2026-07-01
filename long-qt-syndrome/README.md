# 선천성 QT 연장 증후군 (Congenital Long QT Syndrome, LQTS) — QSP Model

> 심실 재분극 지연(QTc 연장)을 유발하는 심장 이온통로병증으로, KCNQ1(IKs,
> LQT1)·KCNH2/hERG(IKr, LQT2)·SCN5A(late-INa, LQT3) 등 원인유전자에 따라
> 서로 다른 유발인자(운동/청각자극/수면)와 베타차단제 반응성을 보이며,
> 조기후탈분극(EAD)·재분극 이질성(dispersion)을 매개로 Torsades de Pointes
> (TdP) -> 심실세동(VF) -> 급사(SCD)로 진행할 수 있다. 후천성(약물유발)
> LQTS와 CiPA(Comprehensive in vitro Proarrhythmia Assay) 다중이온채널
> 차단 패러다임도 함께 모델링하였다.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`lqts_qsp_model.dot`](lqts_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`lqts_qsp_model.svg`](lqts_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`lqts_qsp_model.png`](lqts_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`lqts_mrgsolve_model.R`](lqts_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`lqts_shiny_app.R`](lqts_shiny_app.R) |
| 📚 References             | [`lqts_references.md`](lqts_references.md) |

---

## 1. 배경 및 역학 (Background & Epidemiology)

선천성 LQTS의 유병률은 약 1/2,000명으로 추정되며(Schwartz 2009
*Circulation*), 상염색체 우성 유전(Romano-Ward 증후군)이 대다수이고,
KCNQ1/KCNE1 이중대립유전자 변이로 인한 상염색체 열성형(Jervell-Lange-
Nielsen 증후군)은 선천성 감각신경성 난청을 동반한다. 원인유전자에 따라
LQT1(KCNQ1, ~35-40%), LQT2(KCNH2/hERG, ~25-30%), LQT3(SCN5A, ~5-10%)로
분류되며, 나머지는 ANK2(LQT4), CACNA1C(LQT8, Timothy 증후군),
CALM1-3(LQT14-16), TRDN, AKAP9 등 희귀 아형이다. 국제 LQTS 레지스트리
자료에 따르면 치료받지 않은 고위험군의 심장 사건(실신·심정지·급사) 발생률은
유전형·QTc 길이·성별에 따라 크게 다르며(Goldenberg 2008 *Circulation*;
Sauer 2007 *JACC*), 특히 QTc>500ms, 여성 LQT2, 실신 과거력이 대표적
고위험 인자이다. 약물유발(후천성) LQTS는 hERG 채널을 차단하는 다양한
비심장용 약물(항생제, 항정신병제, 항히스타민제 등)에 의해 발생하며,
CiPA 이니셔티브는 hERG 단일채널 차단이 아닌 다중이온채널(hERG/ICaL/
late-INa) 균형을 반영하는 qNet 유사 지표로 부정맥 위험을 평가한다.

## 2. 발병기전 (Pathophysiology)

* **유전형별 채널 결함**: LQT1은 KCNQ1(+ KCNE1 베타소단위) 기능소실로
  IKs(느린 지연정류 K+ 전류) 감소, 특히 교감신경 항진 시 IKs를 증강시키는
  "안전판(safety-valve)" 기전이 작동하지 않아 운동/수영 중 QTc가 역설적으로
  더 연장된다. LQT2는 KCNH2(hERG)의 이동장애 또는 관통결함으로 IKr(빠른
  지연정류 K+ 전류) 감소, 청각자극·감정·산후기에 취약하다. LQT3는 SCN5A의
  불완전한 불활성화로 지속성/후기 INa가 증가해 활동전위 고평부(phase 2)를
  연장시키며, 안정/수면/서맥 시 위험이 커진다.
* **활동전위 및 재분극 예비능(repolarization reserve)**: IKs·IKr 감소 또는
  late-INa 증가는 모두 재분극 예비능을 감소시켜 활동전위 지속시간(APD)을
  연장시킨다. 심외막·M세포·심내막 간 APD 기울기가 확대되면 재분극
  이질성(TDR)이 커지고, 연장된 phase 2/3 동안 ICaL의 창전류(window current)
  재활성화로 조기후탈분극(EAD)이 발생한다. EAD가 성공적으로 전파되면
  기능적 재진입 회로(phase-2 reentry)가 형성되어 Torsades de Pointes로
  이어질 수 있다.
* **자율신경계 유발인자**: 교감신경 활성은 정상 심장에서 ICaL과 IKs를 함께
  증가시켜 QT를 상대적으로 안정시키지만, LQT1에서는 변이 KCNQ1 채널이 PKA
  인산화에 반응하지 못해 이 균형이 깨진다. LQT2는 청각/정서적 자극에,
  LQT3는 서맥·수면 등 저교감신경 상태에 취약하다.
* **후천성/약물유발 LQTS (CiPA)**: hERG 채널을 차단하는 약물은 선천성
  LQT2와 유사한 표현형을 유발하며, 저칼륨혈증은 세포외 K+ 의존적 IKr
  전도도 감소를 통해 이 효과를 증폭시킨다. CiPA의 qNet 개념은 hERG 차단
  단독이 아니라 ICaL·late-INa 동시차단 정도를 종합해 실제 TdP 위험을
  더 정확히 예측하고자 한다.

## 3. 모델 구조 (Model Structure)

### 3.1 기계론적 지도 — 14개 클러스터, 129개 노드

1. 유전학 — 이온채널 유전형 (KCNQ1/KCNH2/SCN5A/희귀유전자)
2. 이온채널 기능이상 (분자 수준: trafficking, gating, inactivation defect)
3. 심실 활동전위 및 이온전류 (Phase 0-4, INa/Ito/ICaL/IKs/IKr/IK1/INaK/INaCa)
4. 경벽 이질성 및 재분극 분산 (심외막/M세포/심내막/Purkinje, TDR, Tpeak-Tend)
5. EAD/부정맥 기질 및 촉발 (ICaL 재활성화, Ca2+ 과부하, CaMKII, 재진입)
6. 자율신경계 유발인자 (교감신경 항진, 운동/청각/수면 유발인자별 특이성)
7. TdP -> VF -> SCD 케스케이드
8. CiPA/hERG 약물차단 약리학 (후천성 LQTS, qNet 지표)
9. 베타차단제 PK/PD (propranolol, nadolol)
10. Mexiletine PK/PD (late-INa 차단, LQT3 표적치료)
11. 칼륨/전해질 조절 (경구 K+, spironolactone, 저칼륨혈증 위험)
12. 시술/기기 중재 (LCSD, ICD)
13. 임상 진단 (QTc 보정공식, T파 형태, 운동부하검사, Schwartz 점수)
14. 임상 엔드포인트 및 위험도 층화 (QTc, TdP 확률, SCD 위험, 레지스트리)

### 3.2 mrgsolve 모델 — 23개 ODE 구획

* **약물 PK (14개 구획)** — Propranolol(2-구획: gut/central/peripheral),
  Nadolol(2-구획, 신장 청소), Mexiletine(1-구획, CYP2D6[주]+CYP1A2[부]
  대사), 경구 KCl(2-구획), Spironolactone(2-구획), 가상 QT연장 약물
  X(2-구획, CiPA 스타일 hERG 차단제).
* **질환 PD (9개 구획)** — 혈청 K+, 유전형별 채널전도도 지수(GKs/GKr/
  GNa-late), 교감신경 구동 지수(SYMP_DRIVE), 축약형 재분극예비능/QTc
  대리모델(QTC), EAD 확률/기질 지수(EAD_SUBSTRATE), 누적 TdP 위험
  해저드(TDP_HAZARD), 누적 기대 TdP 사건수(TDP_EVENTS, 시간-사건 생존
  분석의 Poisson 강도 적분).
* QTc는 단순화된 lumped turnover 모델로, 유전형별 IKs/IKr 결핍과
  late-INa 과잉을 가중합하여 목표 QTc를 계산하고 1차 turnover로
  수렴시킨다 (전체 40개 상태의 O'Hara-Rudy 활동전위 모델이 아닌, 기전에
  기반한 축소모델).

### 3.3 10개 치료/시뮬레이션 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 무치료 LQT1 + 운동 유발 | Priori 2004 NEJM, Ackerman 1999 Mayo Clin Proc |
| 2 | LQT1 + propranolol 2 mg/kg/day | Vincent 2009 Circulation |
| 3 | LQT1 + nadolol 1 mg/kg/day | Vincent 2009 Circulation (nadolol 비교군) |
| 4 | 무치료 LQT2 + 청각/정서 유발, 여성 | Schwartz 2001 Circulation, Priori 2003 Circulation |
| 5 | LQT2 + 베타차단제 + K+/spironolactone | Schwartz 2001 Circulation (유전형별 효과) |
| 6 | 무치료 LQT3 + 수면/서맥 유발 | Schwartz 1995 Circulation, Zareba 1998 NEJM |
| 7 | LQT3 + mexiletine | Moss 2000/2008 계열, Ruan 2007 Circulation |
| 8 | LQT3 + mexiletine + propranolol 병용 | Schwartz 1995 Circulation (유전형-특이 반응) |
| 9 | 임의 유전형(LQT2) + QT연장 약물 X 추가 (후천성-선천성 중첩) | CiPA (Colatsky 2016, Vicente 2019) |
| 10 | 고위험 LQT2 심정지 후: LCSD + 베타차단제 + ICD 구조 | Schwartz 2004 Circulation, Zareba 2003 |

## 4. Shiny 대시보드 (8개 탭)

1. **환자/유전형 프로파일** — 유전형, 유발인자, 성별, 실신력 설정.
2. **약물 PK** — Propranolol/Nadolol/Mexiletine/QT연장약물 X 혈중농도.
3. **이온채널/PD** — GKs/GKr/GNa-late 전도도 지수, 교감신경 구동, qNet 유사지표.
4. **QTc & ECG 대리지표** — 시간에 따른 QTc 궤적, 위험 역치선.
5. **TdP/부정맥 위험** — EAD 기질 지수, 누적 위험 해저드, 누적 TdP 사건확률.
6. **임상 엔드포인트** — 실신/TdP/SCD 대리 엔드포인트 요약표.
7. **시나리오 비교** — 8개 시나리오 동시 비교(QTc, TdP 확률).
8. **바이오마커/위험도 층화** — Schwartz 위험점수 계산기, 위험 바이오마커 표.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg lqts_qsp_model.dot -o lqts_qsp_model.svg
dot -Tpng -Gdpi=150 lqts_qsp_model.dot -o lqts_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","shinydashboard","DT"))
library(mrgsolve)
mod <- mread("lqts_mrgsolve_model.R")
# 시나리오 헬퍼는 파일 하단 주석의 LQTS_simulate_scenarios() 참고
out <- mod %>% param(GENOTYPE=1, TRIGGER=1) %>%
  ev(amt=40, cmt="GUT_PROP", ii=8/24, addl=270) %>%
  mrgsim(end=2160, delta=6)
plot(out, "QTc_ms,TdP_event_probability,GKs_idx")

# 3) Shiny 대시보드 실행
shiny::runApp("lqts_shiny_app.R")
```

## 6. 참고문헌 요약

`lqts_references.md`에 64편의 문헌을 10개 섹션(유전학, 이온채널
전기생리학, 역학/레지스트리, 유전형-표현형/유발인자, 베타차단제,
mexiletine/LQT3, ICD/LCSD, CiPA/약물유발 LQTS, 진단/위험도층화, 종설)으로
분류하여 수록하였다. 핵심 문헌: Priori 2004 JAMA(유전형별 베타차단제
반응), Vincent 2009 Circulation(propranolol vs nadolol), Moss/Ruan
mexiletine 계열(LQT3 표적치료), Schwartz 2001 Circulation(유전형-유발인자
특이성), Colatsky 2016·Vicente 2019(CiPA qNet 개념).

## 7. 한계점 (Limitations)

* 연구/교육/가설생성용이며 임상적 의사결정에 사용할 수 없다.
* QTc/재분극예비능 모델은 완전한 다중구획 심근세포 활동전위 모델(예:
  O'Hara-Rudy, ten Tusscher)이 아닌 축소된 lumped turnover 모델이다.
  실제 EAD 발생은 확률적(stochastic)이고 세포간 결합(cell-to-cell
  coupling), 조직 수준 파동전파(wave propagation)에 의존하나, 본 모델은
  이를 EAD_SUBSTRATE라는 단일 결정론적 지수로 단순화하였다.
  TdP로의 전환은 Poisson 강도 적분(누적 해저드)으로 근사한 것으로, 실제
  임상시험/레지스트리에서 관찰되는 사건율과 정성적 방향성(유전형별,
  치료별 상대적 위험 감소)은 반영하지만 정량적 절대값은 검증되지
  않았다.
* Mexiletine의 late-INa 차단 EC50/Emax, QT연장 약물 X의 hERG-IC50 등은
  대표적 문헌값에서 유추한 근사 파라미터이며, 개별 환자 데이터로
  적합(fitting)되지 않았다.
* 개체간 변동성(IIV)은 포함하지 않은 typical-value 모델이다.
* 본 컨테이너에는 R/mrgsolve 실행환경이 설치되어 있지 않아 mrgsolve 코드는
  문헌 기반 설계 및 자체 코드 검토(구획/차원/경계값 점검) 단계까지
  완료되었으며, 실제 컴파일·적분 실행으로 수치를 검증하지는 못했다. `.dot`
  파일은 Graphviz `dot`으로 실제 렌더링하여 SVG/PNG 생성을 확인했다.
