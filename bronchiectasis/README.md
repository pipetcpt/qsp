# 기관지 확장증 (Bronchiectasis, BEX) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![BEX QSP Model](bex_qsp_model.png)](bex_qsp_model.svg)

## 개요 (Overview)
기관지 확장증은 기도 감염과 염증의 악순환(Cole의 vicious cycle)에 의해 기관지벽이 비가역적으로 확장·파괴되는 만성 기도 질환입니다. 전 세계 유병률은 인구 10만 명당 약 50~500명이며, 가장 흔한 원인은 소아기 폐렴 후 손상, 원발성 섬모 운동이상증, 면역결핍, 낭포성 섬유증입니다. 호중구 엘라스타제에 의한 기도 구조 손상과 녹농균(P. aeruginosa) 등의 만성 감염·바이오필름 형성이 진행의 핵심입니다. 기도 청결 요법, 흡입 항생제, 아지스로마이신 유지 치료가 표준 관리법입니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 점막 섬모 청소 기능 손상 | 점액 과분비, MCC 감소, CFTR 이상 | 세균 정착, 반복 감염 |
| 세균 과증식 및 바이오필름 | P. aeruginosa alginate, 형광성 바이오필름 | 항생제 내성, 만성 감염 |
| 호중구 주도 염증 | IL-8 → 호중구 모집 → NE, MMP 분비 | 기도벽 구조 파괴 |
| 호중구 엘라스타제(NE) | 엘라스틴·콜라겐 분해, SLPI 무력화 | 비가역적 기관지 확장 |
| 악순환 강화 | 세균 → 더 많은 호중구 → 더 많은 NE | 진행성 폐기능 저하 |
| 급성 악화 | 세균 부하 급증, 사이토카인 폭풍 | 입원, FEV1 가속 손실 |

## 주요 약물 표적 (Drug Targets)
- **아지스로마이신 (장기 저용량)**: 항균 + 항염 (IL-8↓, 중성구 이동↓) — 악화 빈도 감소 (EMBRACE 시험)
- **흡입 토브라마이신 (TIP)**: 폐 내 고농도 항균 → P. aeruginosa 억제 — RESPIRE 시험 기반
- **시프로플록사신 (경구)**: 급성 악화 시 전신 항균 — 14일 요법
- **도르나제 알파 (흡입 DNase)**: 점성 DNA 분해 → 점액 점도 감소, MCC 개선 — 주로 낭포성 섬유증
- **흡입 고장성 식염수/만니톨**: 기도 수분화 → MCC 향상 — 비낭포성 섬유증 BEX

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [bex_qsp_model.dot](bex_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 172 노드 / 15 클러스터) |
| [bex_qsp_model.svg](bex_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [bex_qsp_model.png](bex_qsp_model.png) | PNG 이미지 (150 dpi) |
| [bex_mrgsolve_model.R](bex_mrgsolve_model.R) | mrgsolve ODE 모델 (약 17 구획 / 6개 치료 시나리오) |
| [bex_shiny_app.R](bex_shiny_app.R) | Shiny 대시보드 |
| [bex_references.md](bex_references.md) | 참고문헌 (약 50편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK(아지스로마이신 3구획, 흡입 토브라마이신 2구획, 경구 시프로플록사신 2구획, 도르나제 알파 폐 구획), 질환 상태(세균 부하·바이오필름·호중구·IL-8·NE·MCC·기도 손상·악화 상태)
- **주요 치료 시나리오**: ① 비치료(자연경과), ② 아지스로마이신 250 mg 주3회 유지(EMBRACE), ③ 흡입 토브라마이신 300 mg BID 사이클(RESPIRE), ④ AZM+흡입 토브라마이신 병합, ⑤ AZM+TIP+도르나제 알파 3제, ⑥ 급성 악화 시 시프로플록사신 500 mg BID×14일
- **보정/근거**: EMBRACE(아지스로마이신), RESPIRE-1/2(시프로플록사신 흡입), Chalmers 악화 예측 모델 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(기저질환·세균 배양·FEV1 기저치·악화 빈도), 약동학(각 항생제 혈장·폐 농도), 감염 PD 탭(세균 부하·바이오필름·NE), 기도 염증·손상 탭(IL-8·호중구·MCC·기도 손상 점수), 치료 시나리오 비교(6개), 악화 예측 바이오마커 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("bex_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("bex_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg bex_qsp_model.dot -o bex_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [bex_references.md](bex_references.md) 참조 (약 50편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
