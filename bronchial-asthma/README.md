# 기관지 천식 (Bronchial Asthma, BA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![BA QSP Model](ba_qsp_model.png)](ba_qsp_model.svg)

## 개요 (Overview)
기관지 천식은 전 세계 약 3억 명이 이환된 가장 흔한 만성 기도 질환으로, 알레르기 항원·바이러스·오염물질 노출에 의해 유발되는 가역적 기도 폐쇄, 기도 과민성, 만성 기도 염증을 특징으로 합니다. Th2/호산구 중심의 type 2 염증(IL-4, IL-5, IL-13, TSLP)이 알레르기 천식의 핵심이며, 흡입 코르티코스테로이드(ICS)+LABA가 1차 치료 표준입니다. 중증 호산구 천식에는 IL-5(메폴리주맙·벤라리주맙), IL-4Rα(두필루맙), TSLP(테제펠루맙) 표적 생물학제제가 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TSLP-Th2 활성화 | 상피세포 → TSLP → DC → Th2 분화 | IgE 생성, 알레르기 염증 시작 |
| IL-5 매개 호산구 생성/활성화 | IL-5 → 골수 호산구 분화·혈중 이동 | 혈액·조직 호산구증가증 |
| IL-13 매개 기도 리모델링 | IL-13 → 배세포 과형성, 점액 과분비, AHR | FEV1 저하, 만성 기도 폐쇄 |
| IgE-비만세포 축 | FcεRI → 비만세포 탈과립 → 히스타민·류코트리엔 | 조기 알레르기 반응, 기관지경련 |
| 기도 평활근 수축 | β2-AR 신호(cAMP), 류코트리엔 수용체 | 가역적 기관지 수축 |
| 기도 리모델링 | TGF-β, 기저막 비후, 평활근 비대 | 비가역적 기도 폐쇄 |

## 주요 약물 표적 (Drug Targets)
- **ICS (흡입 코르티코스테로이드)**: GR 매개 기도 염증 전반 억제 — 모든 단계 천식 기본 치료
- **LABA (지속성 β2 작용제)**: β2-AR → cAMP → 기도 평활근 이완 — ICS와 병합
- **메폴리주맙 (항-IL-5)**: IL-5 중화 → 호산구 감소 — 중증 호산구 천식 (MENSA/SIRIUS)
- **벤라리주맙 (항-IL-5Rα)**: IL-5Rα 차단 → 호산구 ADCC 고갈 — q4w → q8w 투여
- **두필루맙 (항-IL-4Rα)**: IL-4/IL-13 공동 차단 — type 2 고가 천식 (QUEST/VENTURE)
- **테제펠루맙 (항-TSLP)**: 상류 TSLP 차단 → type 2 전체 cascade 억제 (NAVIGATOR)
- **오말리주맙 (항-IgE)**: 유리 IgE 중화, FcεRI 하향 조절 — 중증 알레르기 천식

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [ba_qsp_model.dot](ba_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 215 노드 / 16 클러스터) |
| [ba_qsp_model.svg](ba_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ba_qsp_model.png](ba_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ba_mrgsolve_model.R](ba_mrgsolve_model.R) | mrgsolve ODE 모델 (약 28 구획 / 5개 치료 시나리오) |
| [ba_shiny_app.R](ba_shiny_app.R) | Shiny 대시보드 |
| [ba_references.md](ba_references.md) | 참고문헌 (약 41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 생물학제제 PK(메폴리주맙·벤라리주맙·두필루맙·테제펠루맙 각 2구획 SC PK, 오말리주맙 TMDD 포함), ICS/LABA PK(폐·전신), 면역 PD(TSLP·IL-5·IL-13·혈액 호산구·조직 호산구·ASM 긴장도·점액·FEV1)
- **주요 치료 시나리오**: ① ICS/LABA 단독, ② +메폴리주맙, ③ +벤라리주맙, ④ +두필루맙, ⑤ +테제펠루맙
- **보정/근거**: MENSA(메폴리주맙), SIROCCO(벤라리주맙), QUEST(두필루맙), NAVIGATOR(테제펠루맙) 임상시험 FEV1 및 호산구 감소 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(호산구 수치·IgE·천식 중증도·아형 선택), 약동학(각 생물학제제 혈장 농도), 기도 염증 PD(TSLP·IL-5·IL-13·호산구), 임상 엔드포인트(FEV1·악화율·증상 점수), 치료 시나리오 비교(5개 생물학제제 오버레이), 바이오마커 예측(바이오마커별 약물 반응) 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("ba_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ba_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ba_qsp_model.dot -o ba_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [ba_references.md](ba_references.md) 참조 (약 41편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
