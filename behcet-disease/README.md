# 베체트병 (Behçet's Disease, BD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![BD QSP Model](bd_qsp_model.png)](bd_qsp_model.svg)

## 개요 (Overview)
베체트병(BD)은 반복성 구강 궤양, 생식기 궤양, 포도막염, 피부 병변을 특징으로 하는 전신 가변혈관 염증 질환으로, 실크로드를 따라 터키·이란·일본 등에서 유병률이 높습니다(10만 명당 1~300명). 발병기전의 핵심은 호중구 과활성화와 IL-1β/TNF-α 매개 내피세포 손상이며, Th1/Th17 편향과 Treg 기능 부전도 중요합니다. 콜히친이 점막 궤양 예방의 기본 치료이며, 안구 및 혈관 침범 시 TNF 억제제, IL-1 차단제, PDE4 억제제(아프레밀라스트)가 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 호중구 과활성화 | NLRP3 인플라마솜, NET 형성, IL-8 유도 | 피부 과민성, 무균 농포진 |
| IL-1β 매개 염증 | NLRP3 → caspase-1 → IL-1β 성숙 | 구강/생식기 궤양, 발열 |
| TNF-α/IL-6 경로 | Th1 세포, 단핵구 활성화 | 포도막염, 혈관 내피 손상 |
| Th17/IL-17A 경로 | IL-23 → RORγt → IL-17A | 점막 궤양, 호중구 모집 |
| 내피세포 활성화 | E-셀렉틴, VCAM-1 발현 증가 | 혈관 혈전·동맥류 위험 |
| 미세생물군-면역 축 | HLA-B51, 장내 미생물 이상 | 재발 트리거, 지역 차이 |

## 주요 약물 표적 (Drug Targets)
- **콜히친**: 튜불린 중합 억제 → 호중구 주화성 차단 — 구강 궤양 예방 1차 치료
- **프레드니솔론**: GR 매개 광범위 항염 — 급성 눈·혈관 침범 유도 치료
- **아달리무맙/인플릭시맙 (항-TNF)**: TNF-α 중화 → Th1 염증 억제 — 포도막염·혈관 베체트
- **아프레밀라스트 (PDE4 억제제)**: cAMP 증가 → IL-17/TNF 감소 — 구강 궤양, RELIEF 시험 승인
- **카나키누맙 (항-IL-1β)**: NLRP3 하류 차단 — 난치성 BD, 호중구 활성화 억제
- **인터페론-α**: NK세포·Treg 조절 — 포도막염, 일부 지역 사용

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [bd_qsp_model.dot](bd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 158 노드 / 10 클러스터) |
| [bd_qsp_model.svg](bd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [bd_qsp_model.png](bd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [bd_mrgsolve_model.R](bd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 8개 치료 시나리오) |
| [bd_shiny_app.R](bd_shiny_app.R) | Shiny 대시보드 |
| [bd_references.md](bd_references.md) | 참고문헌 (약 42편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK(콜히친 3구획, 프레드니솔론 2구획, 아달리무맙 2구획, 아프레밀라스트, 카나키누맙 2구획), 면역 세포(호중구·Th1·Th17·Treg), 사이토카인(TNF-α·IL-1β·IL-6·IL-17A), 장기 침범 지수(내피세포 활성화, 구강 궤양, 안구 염증, BDCAF 복합지수)
- **주요 치료 시나리오**: ① 비치료, ② 콜히친 단독(1 mg/일), ③ 프레드니솔론 단독, ④ 아달리무맙(항-TNF), ⑤ 아프레밀라스트(PDE4i), ⑥ 카나키누맙(항-IL-1β), ⑦ 콜히친+프레드니솔론 병합, ⑧ 아달리무맙+아프레밀라스트(난치성)
- **보정/근거**: RELIEF 시험(아프레밀라스트), Vitale·Hatemi 메타분석, EULAR 베체트병 권고안 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(침범 장기·HLA-B51·재발 빈도), 약동학(각 약물 혈장 농도-시간), 면역 PD 지표(호중구·사이토카인 추이), 임상 활성 지표(BDCAF, 구강 궤양 빈도, 포도막염 재발), 치료 시나리오 비교(8개 오버레이), 장기 침범 바이오마커 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("bd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("bd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg bd_qsp_model.dot -o bd_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [bd_references.md](bd_references.md) 참조 (약 42편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
