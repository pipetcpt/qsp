# 간신 증후군 (Hepatorenal Syndrome, HRS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도 / 신장

[![HRS QSP Model](hrs_qsp_model.png)](hrs_qsp_model.svg)

## 개요 (Overview)
간신 증후군(Hepatorenal Syndrome, HRS)은 진행된 간경변 및 문맥압 항진증에서 발생하는 기능적, 잠재적으로 가역적인 급성 신손상 형태입니다. 내장(splanchnic) 혈관의 강력한 확장이 유효 순환 혈액량(effective arterial blood volume, EABV)의 감소를 유발하고, 이에 대한 보상으로 RAAS·교감신경계·바소프레신(AVP)이 극대로 활성화되어 신장 혈관을 강하게 수축시키는 것이 핵심 병태생리입니다. 2015년 ICA(International Club of Ascites) 기준과 2019년 개정판에서는 HRS-AKI(급성)와 HRS-NAKI(non-AKI, 아만성/만성)로 재분류되었으며, 30일 사망률이 40–90%에 달하는 치명적 합병증입니다. 2022년 FDA 승인된 테를리프레신을 포함한 혈관수축제와 알부민의 병용요법이 표준치료이며, 근본치료는 간이식(± 신-간 동시이식, SLK)입니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 내장 혈관 확장 | eNOS/iNOS·NO·CO·PGI₂·글루카곤 ↑ | SVR ↓, MAP ↓, splanchnic pooling |
| EABV 감소 | 유효 혈량 저하 → 압수용체 unloading | 심장 하이퍼다이내믹, RAP ↓ |
| 신경호르몬 활성화 | RAAS, 교감신경, AVP ↑↑ | 신 혈관 수축·나트륨/수분 저류 |
| 신장 혈관 수축 | AngII, NE, ET-1, TXA₂, 아데노신, 엔도카나비노이드 | RBF ↓, GFR ↓, sCr ↑ |
| 신장 프로스타글란딘 소실 | PGE₂/PGI₂ ↓ (NSAID 취약) | Afferent 수축 악화 |
| 심경변성 심근병증 | β-adrenergic 하향조절, 확장기 기능부전 | 스트레스 시 CO 부족, HRS 유발 |
| 전신 염증 | 장 세균 전위 · LPS · TNF-α · IL-6 · CRP | SIRS/ACLF, 다장기 기능부전 |
| 유발인자 | SBP, GI 출혈, LVP without albumin, NSAID | 급성 EABV 감소, HRS-AKI 발생 |

## 주요 약물 표적 (Drug Targets)
- **테를리프레신 (Terlipressin)**: 리실-바소프레신 전구체, V1a 아고니스트 → 내장 혈관 수축, MAP·RBF·GFR ↑ (CONFIRM/REVERSE/OT-0401 RCT). 볼루스(1 mg q6h) 또는 지속주입(2 mg/24h).
- **노르에피네프린 (Norepinephrine)**: α1 아드레날린 수용체 → 전신 혈관 수축, ICU 세팅에서 대체 옵션. Terli-nonresponder 또는 부족 상황.
- **미도드린 + 옥트레오타이드 (Midodrine + Octreotide)**: 미도드린(α1, PO) + 옥트레오타이드(글루카곤·VIP 억제, SC) + 알부민. 외래·비ICU 세팅에서 사용.
- **알부민 25% (Albumin)**: 삼투 확장 + 항염증 + 내피 안정화. 1 g/kg (D1) → 20–40 g/day (ATTIRE, ANSWER, SBP 예방 데이터).
- **항생제**: SBP 치료 (세프트리악손), SBP 예방 (노플록사신), 리팍시민 (장세균 전위 억제).
- **TIPS · RRT · 간이식**: 절차적 개입. Bridge-to-transplant 개념.

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [hrs_qsp_model.dot](hrs_qsp_model.dot) | Graphviz 기계론적 지도 소스 (~148 노드 / 15 클러스터) |
| [hrs_qsp_model.svg](hrs_qsp_model.svg) | SVG 벡터 이미지 |
| [hrs_qsp_model.png](hrs_qsp_model.png) | PNG 이미지 (150 dpi) |
| [hrs_mrgsolve_model.R](hrs_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획 / 10 시나리오) |
| [hrs_shiny_app.R](hrs_shiny_app.R) | Shiny 대시보드 (8 탭) |
| [hrs_references.md](hrs_references.md) | 참고문헌 (82편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**:
  - 약물 PK (10): 테를리프레신 2-cpt + 리실-바소프레신, 노르에피네프린 1-cpt, 미도드린→데스글리미도드린, 옥트레오타이드 SC/central, 알부민
  - 신경호르몬 (4): PRA, aldosterone, endogenous NE, AVP
  - 혈역학 (4): MAP, SVR, RBF, cardiac output
  - 신장/임상 (6): GFR, sCr, urine Na, urine output, serum Na, systemic inflammation
  - 결과 (3): 30일 위험도 적분, terlipressin 허혈 AUC, bilirubin
- **주요 치료 시나리오 (10)**:
  1. 자연경과 (natural history)
  2. 테를리프레신 볼루스 (1 mg q6h) + 알부민 40 g/일
  3. 테를리프레신 지속주입 (2 mg/24h) + 알부민
  4. 노르에피네프린 CIV + 알부민
  5. 미도드린 + 옥트레오타이드 + 알부민
  6. 알부민 단독 (ATTIRE-like 비교군)
  7. SBP 유발 후 terli 반응
  8. LVP-without-albumin PICD
  9. NSAID 유발 악화
  10. TIPS 혈역학 대체 시나리오
- **보정/근거**: CONFIRM (Wong 2021, NEJM), REVERSE (Boyer 2016), OT-0401 (Sanyal 2008), ATTIRE (China 2021), ANSWER (Caraceni 2018), Cavallin 2015/2016 (CI vs bolus).

## Shiny 대시보드 (Dashboard, 8 tabs)
1. **Overview** — MAP·RBF·GFR·sCr 트렌드 + 모델 스키마
2. **Drug PK** — 테를리/리실-VP, NE, 미도드린, 옥트레오타이드, 알부민
3. **Neurohormonal** — RAAS (renin, aldosterone), SNS, AVP
4. **Renal & urinary** — urine Na, urine output, serum Na
5. **Clinical endpoints** — HRS 반응, 30/90일 생존, MELD 궤적
6. **Scenario comparison** — 7종 요법 비교 표/그래프
7. **Safety** — Terli ischemic AUC, MAP overshoot 위험
8. **References** — 주요 참고문헌 링크

## 실행 방법 (Usage)
```r
library(mrgsolve); library(shiny)
source("hrs_mrgsolve_model.R")   # builds hrs_mod
shiny::runApp("hrs_shiny_app.R")
```

## 참고문헌 (References)
[hrs_references.md](hrs_references.md) — ICA-2015/2019, CONFIRM, REVERSE, OT-0401, ATTIRE, ANSWER, PREDICT, CANONIC 등 82편.

## 라이선스 & 면책 (License · Disclaimer)
- 상위 저장소 라이선스([../LICENSE](../LICENSE))를 따릅니다.
- 이 모델은 교육·연구용이며, 임상 의사결정을 대체하지 않습니다. 실제 진료는 최신 가이드라인과 전문가 판단을 우선합니다.
