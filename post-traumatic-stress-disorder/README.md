# Post-Traumatic Stress Disorder (외상후 스트레스 장애, PTSD) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking trauma exposure
> and risk/resilience factors (FKBP5 x childhood-adversity interaction,
> peritraumatic dissociation, social support) to noradrenergic/glucocorticoid
> fear-memory consolidation, HPA axis dysregulation (enhanced glucocorticoid
> negative feedback → paradoxically low/normal cortisol), locus
> coeruleus-noradrenergic hyperarousal, and amygdala-hippocampus-vmPFC
> fear-extinction circuit failure (amygdala hyperreactivity, vmPFC
> hypoactivation, impaired extinction recall) driving the four DSM-5
> symptom clusters (intrusion, avoidance, negative cognition/mood,
> hyperarousal) and a composite CAPS-5-like endpoint — coupled to SSRI
> (sertraline/paroxetine), prazosin (alpha-1 antagonist, nightmares),
> ketamine/esketamine (rapid NMDA-antagonist extinction facilitation), and
> MDMA-assisted / trauma-focused psychotherapy (session-dose driven
> extinction-learning boost) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ptsd_qsp_model.dot`](ptsd_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`ptsd_qsp_model.svg`](ptsd_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ptsd_qsp_model.png`](ptsd_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ptsd_mrgsolve_model.R`](ptsd_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ptsd_shiny_app.R`](ptsd_shiny_app.R) |
| 📚 References             | [`ptsd_references.md`](ptsd_references.md) |

---

## 1. 질환 개요 (Disease in one paragraph)

외상후 스트레스 장애(PTSD)는 생명을 위협하거나 심각한 신체적 위해를
동반하는 사건(전투, 성폭력, 교통사고, 자연재해 등)에 노출된 후 발생하는
만성 신경정신질환으로, 편도체를 중심으로 한 노르아드레날린 및
당질코르티코이드 매개 기억 강화 기전이 지나치게 강력한 공포기억을
고착시키는 것에서 시작한다. 핵심 병태생리는 **HPA 축의 역설적 조절이상**
— 증강된 당질코르티코이드 수용체(GR) 음성 되먹임으로 인해 급성기 스트레스에도
불구하고 기저 코르티솔이 오히려 낮거나 정상으로 유지되는 현상(Yehuda
가설) — 과 **청반(locus coeruleus)-노르아드레날린 과각성**의 병존이다.
신경영상 연구는 **편도체 과반응성**과 **복내측전전두엽(vmPFC)
저활성화**가 결합되어 하향식 공포 억제(소거)가 실패하는 회로 기능장애를
일관되게 보여주며, 반복적 스트레스에 의한 **해마 부피 감소**는 맥락 기반
공포 구별 능력을 저하시켜 안전한 상황에서도 위협을 느끼는 과잉일반화를
유발한다. 이러한 회로 이상은 DSM-5의 4개 증상 클러스터 — 침습적 재경험,
회피, 인지/기분의 부정적 변화, 과각성/반응성 — 으로 발현되며 수면장애
(악몽, REM 분절화)와 우울증·물질사용장애·자살 위험 등 동반질환을 동반한다.
1차 약물치료는 FDA 승인 SSRI(**Sertraline, Paroxetine**)이며, 악몽에는
**Prazosin**이 부가요법으로 사용된다. 최근에는 **케타민/에스케타민**의
급속 소거학습 촉진 효과와, 3상 임상에서 유효성이 확인된 **MDMA 보조
심리치료**가 표준 노출기반 심리치료(지속노출치료, 인지처리치료, EMDR)를
보완하는 새로운 축으로 부상하고 있다.

## 2. 기계론적 지도 클러스터 (15개 클러스터, 127개 노드)

1. 위험요인/유전-환경 상호작용 (트라우마 유형, FKBP5 x 아동기 트라우마,
   5-HTTLPR, 주변외상성 해리, 사회적 지지, 회복탄력성)
2. 급성 스트레스 반응 및 공포 습득 (BLA 위협탐지, 노르아드레날린/코르티솔
   매개 기억 강화, 급성 스트레스 장애, 과잉일반화)
3. HPA 축 조절이상 (CRH-ACTH-코르티솔, 증강된 GR 되먹임, 역설적 저코르티솔,
   FKBP5 탈메틸화)
4. 청반-노르아드레날린 과각성 (LC 발화, NE 방출, 알파1/베타 수용체, 심박수/
   HRV, 경악반사, REM 교감 급증)
5. 편도체-해마-vmPFC 회로 기능장애 (BLA 과반응성, vmPFC 저활성, 해마
   위축/신경생성 억제, 뇌섬엽 과활성)
6. 공포기억 습득/소거 신경가소성 (고착, 재고착, 소거학습, BDNF-TrkB/mTOR,
   맥락의존 재발)
7. 신경염증 및 산화스트레스 (IL-6, TNF-α, CRP, 미세아교세포, 키뉴레닌 경로)
8. DSM-5 증상 클러스터 B-E (침습, 회피, 부정적 인지/기분, 과각성, 해리
   하위유형)
9. 수면장애 (REM 분절화, 악몽, 불면증, 서파수면 감소)
10. 동반질환 및 기능적 결과 (우울증, 물질사용장애, 자살위험, 만성통증,
    심혈관위험, 기능장애)
11. SSRI(Sertraline/Paroxetine) PK/PD
12. Prazosin(알파1 차단제) PK/PD
13. Ketamine/Esketamine PK/PD
14. MDMA 보조 심리치료 및 외상중심 심리치료(PE/CPT/EMDR) PK/PD
15. 임상 평가/바이오마커/엔드포인트 (CAPS-5, PCL-5, 관해기준, HRV,
    경악습관화, fMRI)

## 3. mrgsolve 모델 (24 ODE 구획)

* **약물 PK (10개 구획)** — Sertraline 위장관데포/혈장(2), Paroxetine
  위장관데포/혈장(2), Prazosin 위장관데포/혈장(2), Ketamine 혈장/효과구획(2,
  생체상 지연), MDMA 위장관데포/혈장(2).
* **질환/PD (13개 구획)** — 코르티솔(CORTISOL), 청반/NE 긴장도(NE_TONE),
  편도체 반응성 지수(AMYG_REACT), vmPFC 억제긴장도(VMPFC_TONE), 공포기억
  강도(FEAR_MEM), 소거기억 강도(EXT_MEM), 누적 심리치료 용량(THERAPY_CUM),
  DSM-5 4개 클러스터 중증도(INTRUSION/AVOIDANCE/NEGCOG/HYPERAROUSE), 수면장애
  지수(SLEEP_DIST), 복합 CAPS-5 점수(CAPS5).
* **시간 추적 (1개 구획)** — FX_WEEKS.
* 핵심 설계: **SSRI**는 편도체 반응성 목표치(amyg_target)를 직접 낮추고,
  **Prazosin**은 청반 목표치(ne_target)와 수면장애 목표치(sleep_target)만
  낮춰 야간 증상에 선택적으로 작용하도록 분리했다. **Ketamine**은 효과구획
  (생체상, mTOR/시냅스형성 지연 반영)을 통해 소거기억 형성 속도(ext_drive)를
  일시적으로 가속하고, **MDMA**는 급성 처리세션 동안 편도체 목표치를 직접
  둔화시키면서 동시에 소거기억 형성도 가속하는 이중 작용으로 구현했다.
  **누적 심리치료 용량(THERAPY_CUM)**은 비약물성 사건(session dose event)
  으로 증가하며 소거기억 형성속도(ext_boost)를 선형적으로 가속해, 약물과
  심리치료의 상호보완적 병용 효과를 기전적으로 표현한다.

### 10개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 자연경과 - 중등도 트라우마(SEVERITY=1.0), 미치료 | Kessler 1995 Arch Gen Psychiatry |
| 2 | 자연경과 - 중증/반복 트라우마 + FKBP5 + 해리형(SEVERITY=1.5) | Binder 2008 JAMA; Klengel 2013 Nat Neurosci |
| 3 | Sertraline 100-200mg/day | Brady 2000 JAMA |
| 4 | Paroxetine 20-50mg/day | Marshall 2001 Am J Psychiatry |
| 5 | Prazosin 1-15mg qhs 부가요법 | Raskind 2013 Am J Psychiatry; Raskind 2018 NEJM |
| 6 | 외상중심 심리치료 주1회 x12(PE/CPT/EMDR) | Foa 2005 JAMA; Resick 2002 JCCP |
| 7 | Ketamine 0.5mg/kg IV x6(2주) | Feder 2014 JAMA Psychiatry; Feder 2021 Am J Psychiatry |
| 8 | MDMA 보조 심리치료(3세션+12회 준비/통합) | Mitchell 2021 Nat Med; Mitchell 2023 Nat Med |
| 9 | 고회복탄력성 + 경증 트라우마(자연관해) | Southwick 2014 Eur J Psychotraumatol |
| 10 | 병용: SSRI + 주간 심리치료 + Prazosin | Krystal 2017 Biol Psychiatry (병용 근거) |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — 트라우마 중증도, FKBP5 위험, 해리 하위유형,
   회복탄력성, 시뮬레이션 기간.
2. **PK** — 약물별(Sertraline/Paroxetine/Prazosin/Ketamine/MDMA) 혈장 농도.
3. **PD 주요지표(공포회로/HPA)** — 편도체/vmPFC/공포-소거기억, 코르티솔/NE
   긴장도.
4. **임상 엔드포인트** — CAPS-5 궤적(관해 임계선 포함), 증상 클러스터별
   추이, 관해도달시간.
5. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
6. **바이오마커** — 코르티솔, NE 긴장도, 편도체 반응성, vmPFC 긴장도, 소거
   기억, 수면장애 지수.
7. **수면/야간증상** — 수면장애 지수 궤적(Prazosin 효과 확인).
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg ptsd_qsp_model.dot -o ptsd_qsp_model.svg
dot -Tpng -Gdpi=150 ptsd_qsp_model.dot -o ptsd_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("ptsd_mrgsolve_model.R") %>% param(TRAUMA_SEVERITY = 1.0)
e_sert <- ev(amt = 100, cmt = "SERT_GUT", time = 0, ii = 24, addl = 180)
e_pe   <- ev(amt = 1, cmt = "THERAPY_CUM", time = 168, ii = 168, addl = 11)
out <- mod %>% ev(e_sert) %>% ev(e_pe) %>% mrgsim(end = 24*7*52, delta = 24)  # 52주 추적
plot(out, c("AMYG_REACT", "VMPFC_TONE", "EXT_MEM", "CAPS5"))

# 3) Shiny 대시보드 실행
shiny::runApp("ptsd_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| 생애 유병률 및 만성화 경과 | 전국 역학조사(NCS) | Kessler 1995 Arch Gen Psychiatry |
| Sertraline RCT 반응률 | 무작위 대조 임상시험 | Brady 2000 JAMA |
| Paroxetine 고정용량 RCT | 무작위 대조 임상시험 | Marshall 2001 Am J Psychiatry |
| Prazosin 악몽/수면 개선(1차) vs 무효(대규모 확증, 이질적 반응) | RCT 2건 | Raskind 2013 Am J Psychiatry; Raskind 2018 NEJM |
| Ketamine 단회/반복 투여 급속 증상 감소 | RCT 2건 | Feder 2014 JAMA Psychiatry; Feder 2021 Am J Psychiatry |
| MDMA 보조 심리치료 3상 유효성(2건 확증) | 무작위 대조 3상 | Mitchell 2021 Nat Med; Mitchell 2023 Nat Med |
| 지속노출/인지처리치료 효과크기 | RCT 및 메타분석 | Foa 2005 JAMA; Watts 2013 J Clin Psychiatry |
| FKBP5 x 아동기 트라우마 위험 상호작용 | 유전-환경 연관 연구 | Binder 2008 JAMA; Klengel 2013 Nat Neurosci |
| 편도체 과반응성/vmPFC 저활성 회로 | 신경영상 메타분석 | Shin & Liberzon 2010 Neuropsychopharmacology |
| 해마 부피 감소 | 신경영상 메타분석 | Woon 2010 Prog Neuropsychopharmacol Biol Psychiatry |

## 7. 모델 검증 상태

이 컨테이너에는 초기 R/mrgsolve 실행환경이 설치되어 있지 않아(`Rscript`
부재), mrgsolve 모델은 **문헌 기반 파라미터 설계 및 코드 자체검토(24개
컴파트먼트 전부 `$CMT`/`$ODE`/`$INIT`와 1:1 대응 확인, SSRI/Prazosin/
Ketamine/MDMA 각각의 작용점을 편도체·청반·효과구획·급성세션효과로 의도적으로
분리한 로직 검토)** 단계까지 완료되었으며 실제 컴파일·적분 실행으로 수치를
검증하지는 못했다. `.dot` 파일은 `apt-get install graphviz`로 설치한
Graphviz `dot`으로 실제 렌더링해 SVG/PNG를 생성·확인했다(127 노드, 15
클러스터). 참고문헌은 WebSearch로 저자/연도/저널 정보를 PubMed와 대조하여
PMID를 개별 검증했다. mrgsolve/R 환경이 있는 곳에서 위 "실행 방법"대로
실행해 수치 적분 결과를 확인할 것을 권장한다.
