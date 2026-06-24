# 에반스 증후군 (Evans Syndrome, ES) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈액

[![ES QSP Model](es_qsp_model.png)](es_qsp_model.svg)

## 개요 (Overview)

에반스 증후군(Evans Syndrome)은 자가면역 용혈성 빈혈(AIHA)과 면역성 혈소판감소증(ITP)이 동시 또는 순차적으로 발생하는 희귀 다계열 자가면역 혈구감소 질환입니다. 유병률은 인구 10만 명당 약 0.02명으로 매우 드물며, 어린이와 성인 모두에서 발생합니다. 자가반응 B세포의 항적혈구(항IgG 온난형) 및 항혈소판 자가항체 생성, 비장에서의 항체 매개 혈구 파괴, Treg 기능 결손이 핵심 병태생리입니다. 치료는 스테로이드·리툭시맙·IVIG·혈소판 생성 촉진제(TPO-RA)를 단계적으로 사용합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가반응 B세포 | 자가항원에 반응하는 B세포 클론 증식 | 항RBC·항혈소판 자가항체 생성 |
| Treg 결손 | FoxP3+ Treg 수·기능 감소 → 말초 내성 이상 | 자가면역 반응 지속 |
| Fc수용체 매개 혈구 파괴 | IgG 피복 적혈구·혈소판 → 비장 대식세포 탐식 | 용혈·혈소판감소 |
| 보체 매개 용혈 | 항IgM(한랭형) → 보체 → 혈관내 용혈 | 급성 빈혈, 혈색소뇨 |
| 혈소판 생성 억제 | 자가항체·T세포 → 거핵구 성숙 억제 | 혈소판감소 가중 |
| PI3K-δ 경로 | (일부 ALPS 연관) → 림프구 자멸사 저항 | 자가면역 림프구 축적 |
| 이차 면역결핍 | 면역억제 치료 합병증 | 감염 위험 증가 |

## 주요 약물 표적 (Drug Targets)

- **코르티코스테로이드** (프레드니솔론): 자가항체 생성 억제, 비장 탐식 억제, 1차 치료
- **리툭시맙(Rituximab)**: 항CD20 → B세포 고갈 → 자가항체 감소, 2차 치료
- **IVIG**: Fc수용체 봉쇄 → 항체 매개 혈구 파괴 차단, 급성기 조절
- **TPO 수용체 작용제** (엘트롬보팍, 로미플로스팀): 거핵구 자극 → 혈소판 생성 촉진
- **시롤리무스(mTOR 억제제)**: ALPS 연관 에반스 증후군에서 자가반응 림프구 억제
- **MMF(마이코페놀레이트 모페틸)**: 림프구 증식 억제, 스테로이드 감량 유지

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [es_qsp_model.dot](es_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 163 노드 / 10 클러스터) |
| [es_qsp_model.svg](es_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [es_qsp_model.png](es_qsp_model.png) | PNG 이미지 (150 dpi) |
| [es_mrgsolve_model.R](es_mrgsolve_model.R) | mrgsolve ODE 모델 (약 26 구획 / 5 치료 시나리오) |
| [es_shiny_app.R](es_shiny_app.R) | Shiny 대시보드 |
| [es_references.md](es_references.md) | 참고문헌 (약 56편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(프레드니솔론·MPA·엘트롬보팍·시롤리무스·IVIG 각 depot/중심) + PD 구획(자가반응 B세포, 형질세포, 자가항체(항RBC·항혈소판), Treg/Teff 비율, 적혈구 수, 혈색소, 혈소판 수, 거핵구, 비장 탐식 지수, 보체 활성화)
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 코르티코스테로이드 단독(1차), ③ 리툭시맙 + 프레드니솔론(2차), ④ IVIG(급성기 구제), ⑤ 시롤리무스 + 스테로이드(불응성 ALPS 연관)
- **보정/근거**: Michel 등 에반스 증후군 코호트, Aladjidi 소아 AIHA 데이터, 리툭시맙 AIHA/ITP 메타분석 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 혈색소·혈소판·Coombs 검사·Treg 수준 설정), ② PK 탭(스테로이드·리툭시맙·IVIG·TPO-RA 혈중 농도), ③ 면역 PD 탭(B세포·자가항체·Treg/Teff 비율·보체 추이), ④ 임상 엔드포인트(혈색소·혈소판·Hb 관해 달성률), ⑤ 시나리오 비교(5개 치료 전략 동시 비교), ⑥ 바이오마커(형질세포·거핵구·비장 탐식 지수 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("es_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("es_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg es_qsp_model.dot -o es_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [es_references.md](es_references.md) 참조 (약 56편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
