# 요로결석 (만성 재발성) (Urolithiasis, URI) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![URI QSP Model](uri_qsp_model.png)](uri_qsp_model.svg)

## 개요 (Overview)

만성 재발성 요로결석은 소변 내 칼슘·수산염·요산·구연산의 과포화 → 결정 핵형성 → Randall 플라크 기저 부착 → 결석 성장의 경로로 발생하며, 평생 재발률이 50%에 달하는 흔한 비뇨기 질환이다. 전 세계 유병률은 약 10~15%이며 선진국에서 증가 추세이다. 원발성 고수산뇨증(PH1)과 같은 유전성 원인부터 고칼슘뇨증·고요산뇨증 등 대사 이상까지 원인이 다양하다. 수분 섭취 증가와 저염·저단백 식이가 기본이며, 약물로는 티아지드계 이뇨제·구연산칼륨·알로푸리놀·탐수로신 등이 결석 유형에 따라 선택된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 소변 과포화 | 칼슘·수산염·인산·요산 과잉 배설 → 이온 활동도 곱 상승 | CaOx·CaP·요산 결정 핵형성 |
| Randall 플라크 | 신유두 기저막 인회석 침착 → 수산칼슘 결석 부착 거점 | 재발성 결석 성장 위험 증가 |
| 구연산 억제 감소 | 저구연산뇨증 → 결정 억제인자 부족 | 결석 성장 촉진, 재발 위험 증가 |
| PH1 고수산뇨증 | AGXT 돌연변이 → 간 수산 과잉 생산 → 심한 신장 수산 침착 | 조기 신부전, 신장·전신 수산증 |
| 신장 염증·손상 | 결정 침착 → 산화 스트레스 → 세뇨관 상피 손상 → 염증 | eGFR 감소, 만성 신부전 진행 |
| 요 pH 조절 | pH 5.5 이하 → 요산 결정화; pH 7 이상 → CaP 침전 위험 | 결석 유형에 따른 pH 목표 상이 |

## 주요 약물 표적 (Drug Targets)

- **티아지드계 이뇨제 — HCTZ 25 mg/day**: 원위 세뇨관 칼슘 재흡수 촉진 → 요칼슘 감소 (고칼슘뇨증형)
- **구연산칼륨 — 60 mEq/day**: 요 구연산 증가 + 요 알칼리화 → CaOx·요산 결석 억제
- **알로푸리놀 — 300 mg/day**: XO 억제 → 요산 생산 감소 → 요산 결석·고요산뇨성 CaOx 결석 예방
- **루마시란(Lumasiran)** — siRNA 기반 AGXT 보완 치료: PH1에서 간 HAO1 억제 → 수산 생산 감소 (ILLUMINATE-A/B)
- **탐수로신 — α1차단제**: 요관 평활근 이완 → 결석 자연 배출 촉진

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [uri_qsp_model.dot](uri_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 431 노드 / 10 클러스터) |
| [uri_qsp_model.svg](uri_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [uri_qsp_model.png](uri_qsp_model.png) | PNG 이미지 (150 dpi) |
| [uri_mrgsolve_model.R](uri_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 6 치료 시나리오) |
| [uri_shiny_app.R](uri_shiny_app.R) | Shiny 대시보드 |
| [uri_references.md](uri_references.md) | 참고문헌 (약 38편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(HCTZ 3구획, 알로푸리놀+옥시푸리놀 3구획, 구연산칼륨 2구획, 탐수로신 2구획) + 생리 구획(혈장 Ca·수산·요산·구연산 4구획) + 결석 크기·신장 염증·eGFR 3구획으로 총 약 19개 구획
- **주요 치료 시나리오**: (1) 무치료(CaOx 결석 형성자), (2) HCTZ 25 mg/day, (3) 구연산칼륨 60 mEq/day, (4) 알로푸리놀 300 mg/day, (5) PH1 + 루마시란, (6) 생활 습관 개선 + 병용 요법
- **보정/근거**: ILLUMINATE-A/B(루마시란 PH1), AUA 결석 가이드라인 24시간 요 분석 데이터, 티아지드 메타분석(요칼슘 감소 효과) 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·수분 섭취량·결석 유형·기저 eGFR·유전성 원인 여부 설정) · PK 시각화(각 약물 혈중 농도-시간 곡선) · 소변 화학 PD 지표(요칼슘·요수산·요구연산·요 pH 시계열) · 임상 엔드포인트(결석 크기 궤적·eGFR 변화·재발률) · 치료 시나리오 비교(6개 요법 장기 결석 성장 억제) · 바이오마커 패널(혈장 Ca·요산·구연산·신장 염증 지수) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("uri_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("uri_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg uri_qsp_model.dot -o uri_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [uri_references.md](uri_references.md) 참조 (약 38편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
