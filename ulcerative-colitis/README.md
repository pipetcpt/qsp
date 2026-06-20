# 궤양성 대장염 (Ulcerative Colitis, UC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![UC QSP Model](uc_qsp_model.png)](uc_qsp_model.svg)

## 개요 (Overview)

궤양성 대장염(UC)은 대장 점막에 국한된 만성 재발성 염증성 장질환으로, 직장에서 시작하여 근위부 대장으로 연속적으로 확장된다. 서양권 유병률은 인구 10만 명당 약 100~300명이며 아시아·한국에서도 꾸준히 증가 추세이다. 핵심 발병기전은 장벽 기능 장애 → 장내 미생물 항원 노출 → Th2/Th17 주도의 점막 염증 → TNF-α·IL-17·IL-13 과분비로 요약된다. 5-ASA, 생물학적 제제(항TNF·항인테그린·항IL-12/23), JAK 억제제, S1P 조절제 등 다양한 기전의 치료제가 단계적으로 사용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 장벽 기능 장애 | Claudin·Occludin 발현 저하, mucin 층 감소 | 장내 세균 항원 점막 침투 |
| Th2 주도 염증 | TSLP·IL-33 → ILC2 → IL-13·IL-5 과분비 | 점막 상피 손상, 호산구 침윤 |
| TNF-α 경로 | 대식세포·Th1 → TNF-α → NF-κB → 점막 궤양 | 혈변·점액 변·복통, Mayo 점수 상승 |
| IL-17 경로 | Th17 → IL-17A → 중성구 동원, CXCL8 유도 | 장 점막 호중구 침윤, 대장염 악화 |
| 조절 T세포 부전 | Treg→IL-10 경로 약화, 면역 관용 파탄 | 만성 재발성 경과 |
| 장 귀소 경로 | α4β7-MAdCAM-1 축, CCR9-CCL25 → 림프구 대장 귀소 | 표적 치료 가능 경로(베돌리주맙) |

## 주요 약물 표적 (Drug Targets)

- **항TNF-α — 인플릭시맙(Infliximab)**: TNF-α 중화, 점막 치유 유도 (ACT-1/2 trials)
- **항α4β7 인테그린 — 베돌리주맙(Vedolizumab)**: 장 선택적 림프구 귀소 차단 (GEMINI 1/2)
- **JAK 억제제 — 토파시티닙(Tofacitinib)**: JAK1/3 억제, 사이토카인 신호 차단 (OCTAVE Induction/Sustain)
- **항IL-12/23 — 우스테키누맙(Ustekinumab)**: IL-12·IL-23 p40 차단, Th1/Th17 억제 (UNIFI)
- **S1P 수용체 조절제 — 오자니모드(Ozanimod)**: 림프구 이차 림프기관 억류, 장 귀소 차단 (TRUE NORTH)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [uc_qsp_model.dot](uc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 422 노드 / 12 클러스터) |
| [uc_qsp_model.svg](uc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [uc_qsp_model.png](uc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [uc_mrgsolve_model.R](uc_mrgsolve_model.R) | mrgsolve ODE 모델 (약 35 구획 / 6 치료 시나리오) |
| [uc_shiny_app.R](uc_shiny_app.R) | Shiny 대시보드 |
| [uc_references.md](uc_references.md) | 참고문헌 (약 47편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(인플릭시맙 2구획+TMDD, 베돌리주맙 1구획+TMDD, 토파시티닙 경구 2구획, 오자니모드+활성대사체 3구획, 우스테키누맙 SC 2구획) + 사이토카인(TNFα·IL-17·IL-13·IL-10) + 면역세포(Th2·Th17·Treg·중성구) + 질환 지표(Mayo 점수·점막 치유 지수·CRP·대변 칼프로텍틴)으로 총 약 35개 구획
- **주요 치료 시나리오**: (1) 위약, (2) 인플릭시맙 5 mg/kg IV 유도·유지, (3) 베돌리주맙 SC, (4) 토파시티닙 10 mg BID 유도 → 5 mg BID 유지, (5) 우스테키누맙 IV 유도 → SC 유지, (6) 오자니모드 0.92 mg/day
- **보정/근거**: ACT-1/2(인플릭시맙 점막 치유율), GEMINI(베돌리주맙 임상 관해율), OCTAVE(토파시티닙 Mayo 점수 반응), UNIFI(우스테키누맙), TRUE NORTH(오자니모드) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·Mayo 점수 기저치·CRP·대변 칼프로텍틴 설정) · PK 시각화(각 생물학적 제제·소분자 혈중 농도-시간 곡선) · 점막 PD 지표(Mayo 점수·점막 치유 지수 시계열) · 임상 엔드포인트(임상 관해·내시경 반응률 비교) · 치료 시나리오 비교(6개 요법 장기 효능) · 바이오마커 패널(사이토카인·면역세포 프로파일) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("uc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("uc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg uc_qsp_model.dot -o uc_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [uc_references.md](uc_references.md) 참조 (약 47편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
