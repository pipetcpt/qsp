# Trigeminal Neuralgia (삼차신경통, TN) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking neurovascular
> compression (NVC) of the trigeminal root entry zone (REZ) to focal
> demyelination, voltage-gated sodium channel (Nav1.3/1.6/1.7)
> upregulation, ephaptic crosstalk, ectopic/afterdischarge generation, and
> central trigeminal-nucleus sensitization that together produce
> paroxysmal lancinating facial pain — coupled to anticonvulsant PK/PD
> (carbamazepine with autoinduction, oxcarbazepine prodrug/MHD, baclofen,
> gabapentin, pregabalin) and interventional therapy (microvascular
> decompression, percutaneous radiofrequency rhizotomy) with time-dependent
> recurrence.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`tn_qsp_model.dot`](tn_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`tn_qsp_model.svg`](tn_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`tn_qsp_model.png`](tn_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`tn_mrgsolve_model.R`](tn_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`tn_shiny_app.R`](tn_shiny_app.R) |
| 📚 References             | [`tn_references.md`](tn_references.md) |

---

## 1.질환 개요 (Disease in one paragraph)

삼차신경통은 삼차신경 근진입부(REZ)에서 발생하는 박동성 동맥(주로 상소뇌동맥)
또는 정맥의 만성 압박이 국소 탈수초화를 유발하는 것이 핵심 병인으로, 노출된
나출축삭(naked axon)에서 Nav1.7/1.3/1.6 나트륨채널의 이상 발현·재분포가
발생해 이소성 자발방전과 후발방전(afterdischarge)을 생성한다(Devor의
"ignition hypothesis"). 탈수초 부위의 인접 축삭 간 전기적 신호누출
(ephaptic crosstalk)로 경촉각(Aβ) 신호가 통각(Aδ/C) 경로로 오배선되어
얼굴을 살짝 만지거나 씹기·양치 등의 무해자극이 전기충격 같은 통증 발작
(paroxysm)을 유발한다. 반복적 구심입력은 삼차신경핵 미부(Vc)에서 NMDA
매개 wind-up과 중추감작을 일으켜 통증의 강도·지속성·이질통을 증폭시킨다.
1차 약물치료는 나트륨채널 차단제인 카바마제핀·옥스카르바제핀이며, 반응
불충분 시 바클로펜(GABA-B)·가바펜틴/프레가발린(α2δ) 병용요법을 추가한다.
약물난치성 환자는 미세혈관감압술(MVD, 병인 제거)이나 감마나이프·경피
고주파 신경근절제술·풍선압박술 등 파괴적/비파괴적 시술로 전환하며, 시술
후에도 시간에 따라 재발할 수 있다.

## 2. 기계론적 지도 클러스터 (19개 클러스터, 108개 노드)

1. 병인/위험인자 (고령·고혈압·혈관확장·유전소인·MS/종양/AVM 동반)
2. 신경혈관 압박 (SCA/AICA/정맥 압박, REZ, 박동성 압박)
3. 국소 탈수초화 (슈반세포 손상, 나출축삭, MRI 소견)
4. Nav 채널병증 (Nav1.7/1.3/1.6/1.8, 이소성방전, 후발방전, 임계값저하)
5. 말초 신호누출/과민화 (ephaptic crosstalk, Aβ→Aδ/C 오배선, 국소 신경염증)
6. 중추 경로/전달 (삼차신경절→Vc→삼차시상로→VPM→체감각피질, 하행조절)
7. 중추감작 (NMDA, wind-up, 아교세포 활성화, 이질통)
8. 이차성 TN 병인 (MS 플라크, 뇌교소각종양, 양측성/감각소실 red flag)
9. Carbamazepine PK (자가유도, 활성대사체 epoxide, HLA-B*1502)
10. Oxcarbazepine PK (전구약물→MHD, 저나트륨혈증)
11. 병용/2차 약물 PK (바클로펜, 라모트리진, 가바펜틴, 프레가발린, BoNT-A)
12. 약물 PD (Nav 차단, 시냅스 글루타메이트 감소, α2δ/GABA-B, SNARE 절단)
13. 시술/수술 (MVD, 감마나이프, RF rhizotomy, 풍선압박, 글리세롤)
14. 이상반응/안전성 (졸림, 간독성, 저나트륨혈증, SJS/TEN, 무과립구증, 시술후 무감각통)
15. 임상 엔드포인트 (발작빈도, NRS, BNI, 완전관해율, 재발까지 시간, QoL)
16. 자연경과 (관해-재발성 경과, 진행성 악화, 비전형 TN)
17. 삼차신경 분지/임상양상 (V1/V2/V3, ICHD-3, 고해상도 MRI, 역학)
18. 정신건강/기능 동반영향 (우울·불안·수면장애·기능장애·돌봄부담)
19. 2차 약물요법 (토피라메이트, 페니토인, 레베티라세탐, 수마트립탄)

## 3. mrgsolve 모델 (17 ODE 구획)

* **약물 PK (5종, 12개 구획)** — Carbamazepine(위장관/중심/epoxide/효소유도상태 4구획,
  1구획 turnover 모델로 자가유도 반영), Oxcarbazepine(위장관/MHD 2구획),
  Baclofen(위장관/중심 2구획), Gabapentin(포화성 흡수 위장관/중심 2구획),
  Pregabalin(선형 위장관/중심 2구획).
* **질환/PD (5개 구획)** — Nav 채널 상향발현 지수(NAV_UPREG, 잔존 NVC 구동),
  이소성방전(ECTOPIC), 중추감작지수(CENTSENS), 발작빈도(PAROX, 회/일),
  통증강도(PAIN, NRS 0-10).
* **안전성/시술 (5개 구획)** — 혈장 Na+(저나트륨혈증), 졸림/운동실조 점수,
  MVD_STATE·RF_STATE(시술 후 관해→재발 진행 상태변수).
* Nav 차단은 CBZ(모약+epoxide 결합노출)와 OXC(MHD)의 Emax 모델을 결합
  (`1-(1-block_cbz)(1-block_oxc)`)하며, 중추감작 억제는 GABA-B(바클로펜)와
  α2δ(가바펜틴+프레가발린) 기전의 가산적 억제로 모델링.

### 7개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 무치료 자연경과 | Maarbjerg 2014 Headache 자연경과 코호트 |
| 2 | 카바마제핀 단독요법 (200mg BID→400mg TID 적정) | Zakrzewska 1989 JNNP 교차연구, Wiffen 2014 Cochrane |
| 3 | 옥스카르바제핀 단독요법 (300mg BID→1200mg/day) | Zakrzewska 1989 JNNP, Besi 2015 |
| 4 | 카바마제핀+바클로펜 병용 (난치성) | Fromm 1984 Ann Neurol RCT |
| 5 | MVD 시술 (14일차, 술후 CBZ 감량) | Barker 1996 NEJM, Sindou 2006 |
| 6 | CBZ 불내성→가바펜틴+프레가발린 전환 | Al-Quliti 2015 review |
| 7 | 경피 고주파 신경근절제술(RF, 30일차)·재발 추적 | Tronnier 2001 Neurosurgery |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — NVC 중증도, 침범 분지, 이차성(MS) 여부 조절.
2. **PK** — 약물별 혈중농도(CBZ/epoxide/OXC-MHD/바클로펜/가바펜틴/프레가발린).
3. **경로 PD** — Nav 채널차단 분율, 중추감작지수.
4. **임상 엔드포인트** — 발작빈도, 통증 NRS, BNI 근사 등급.
5. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
6. **바이오마커** — Nav 채널 발현/이소성흥분성 지수(영상 소견 연계 개념).
7. **안전성** — 혈장 Na+(저나트륨혈증 경고선), 졸림/운동실조 점수.
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg tn_qsp_model.dot -o tn_qsp_model.svg
dot -Tpng -Gdpi=150 tn_qsp_model.dot -o tn_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("tn_mrgsolve_model.R")
source("tn_mrgsolve_model.R")  # run_scenarios() 헬퍼 로드용(선택)
results <- run_scenarios(mod)
plot(results$cbz_mono %>% mrgsolve::filter_sims(time <= 720), c("PAROX","PAIN","CBZ_conc"))

# 3) Shiny 대시보드 실행
shiny::runApp("tn_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| CBZ 자가유도 (t1/2 36h→12-17h) | 반복투여 후 청소율 증가 | Bertilsson & Tomson 1986 Clin Pharmacokinet |
| OXC 저나트륨혈증 발생률 ~2.7% | CBZ 대비 상대위험 | Dong 2005 Neurology |
| 바클로펜 병용 효과 | CBZ 단독 대비 발작감소 | Fromm 1984 Ann Neurol |
| MVD 장기 무통생존 | ~70% 10년 무통, 연간 재발 ~1-4% | Barker 1996 NEJM; Sindou 2006 |
| RF rhizotomy/감마나이프 재발 | MVD 대비 빠른 재발 경향 | Tronnier 2001; Kondziolka 1996 |
| HLA-B*1502 SJS/TEN 위험 | 아시아인 CBZ 개시 전 검사 권고 | Chung 2004 Nature; Ferrell 2008 |

## 7. 모델 검증 상태

이 컨테이너에는 R/mrgsolve 실행환경이 설치되어 있지 않아 (`Rscript` 부재),
mrgsolve 모델은 **문헌 기반 파라미터 설계 및 코드 자체검토(차원/한계값
검사)** 단계까지 완료되었으며 실제 컴파일·적분 실행으로 수치를 검증하지는
못했다. `.dot` 파일은 Graphviz `dot`으로 렌더링해 실제로 SVG/PNG를
생성·확인했다. mrgsolve/R 환경이 있는 곳에서 위 "실행 방법"대로 실행해
`PAROX`·`PAIN`·`NA_PLASMA` 등의 궤적이 기대한 방향(치료 시 발작빈도·통증
감소, OXC에서 저나트륨혈증 저하 등)으로 움직이는지 확인 권장.

## 8. 한계 (Caveats)

* 연구/교육/가설생성용이며 임상적 의사결정에 사용할 수 없음.
* Nav 채널 차단·중추감작 억제의 EC50/Emax는 문헌에서 직접 추정되지 않은
  대리(surrogate) 파라미터로, 상대적 방향성/순위를 보이기 위한 것.
* MVD/RF 재발률은 문헌의 생존곡선을 단일 지수함수로 단순화한 근사치.
* 개체간 변동성(IIV)은 포함하지 않은 typical-value 모델.
