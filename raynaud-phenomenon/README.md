# 레이노 현상 (Raynaud's Phenomenon) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관·결합조직

[![Raynaud QSP Model](raynaud_qsp_model.png)](raynaud_qsp_model.svg)

## 개요 (Overview)

레이노 현상(Raynaud's Phenomenon)은 한랭 자극이나 감정적 스트레스에 의해 유발되는 소동맥 및 세동맥의 간헐적 혈관경련으로, 손가락 및 발가락에 삼색 변화(창백→청색→홍조)를 일으킵니다. 일차성(특발성)과 이차성(전신경화증·루푸스·MCTD 동반)으로 분류되며, 전 세계 유병률은 3-5%입니다. 핵심 기전은 α2C-아드레날린 수용체 과민성, 내피세포 기능 장애(eNOS↓·ET-1↑·PGI2↓), RhoA/ROCK 경로를 통한 칼슘 감수성 증가, 신경펩타이드 불균형(CGRP↓·NPY↑)이며, 이차성에서는 혈관 리모델링과 디지털 궤양이 추가됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| α2C-AR 활성화 | 한랭→α2C-AR 상향조절→NE→혈관수축 | 혈관경련 발작 |
| RhoA/ROCK 경로 | RhoGEF→RhoA-GTP→ROCK→MYPT1 인산화→MLC 인산화(Ca²⁺ 비의존성) | 칼슘 감수성 증가, 지속적 수축 |
| 내피세포 NO 결핍 | eNOS↓→NO↓→sGC/cGMP↓→MLCP↓ | 혈관 이완 손상 |
| 엔도텔린-1 과잉 | ETA-R→Gq/Gα12→IP3/RhoA→VSMC 수축 | 혈관수축 + 리모델링 |
| PGI2 결핍 | COX-2↓→PGI2↓→IP-R/cAMP/PKA↓ | 혈소판 응집↑, 혈관이완↓ |
| 신경펩타이드 불균형 | CGRP↓(한랭 차단), NPY↑ co-release | 신경성 혈관이완 소실 |
| 산화 스트레스 | NOX2/XO→ROS→NF-κB→ET-1↑·eNOS 탈결합 | 악순환 산화 손상 |
| 혈소판 활성화 | TXA2↑·5-HT·ADP → GP IIb/IIIa | 미세혈전, 청색증 지속 |

## 주요 약물 표적 (Drug Targets)

- **칼슘 채널 차단제 (니페디핀, 암로디핀)**: L형 VGCC 차단 → [Ca²⁺]i↓; 1차 치료
- **PDE5 억제제 (실데나필, 타달라필)**: cGMP 분해 억제 → 혈관이완; 2차 치료 또는 이차성
- **엔도텔린 수용체 길항제 (보센탄)**: ETA/ETB 차단 → ET-1 효과 차단; 이차성·디지털 궤양 예방
- **프로스타사이클린 유사체 (일로프로스트)**: IP 수용체 작용제 → cAMP↑; 이차성 중증
- **α1-차단제 (프라조신)**: α1-AR 차단 → 교감신경성 혈관수축↓
- **기타 (플루옥세틴, 로사르탄)**: 간접 혈관활성 효과

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [raynaud_qsp_model.dot](raynaud_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드 / 13 클러스터) |
| [raynaud_qsp_model.svg](raynaud_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [raynaud_qsp_model.png](raynaud_qsp_model.png) | PNG 이미지 (150 dpi) |
| [raynaud_mrgsolve_model.R](raynaud_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획 / 9개 시나리오) |
| [raynaud_shiny_app.R](raynaud_shiny_app.R) | Shiny 대시보드 (7탭) |
| [raynaud_references.md](raynaud_references.md) | 참고문헌 (62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 니페디핀/실데나필/보센탄/일로프로스트/프라조신의 PK 10구획 + NE·RhoA·Cai·cGMP·cAMP·ET-1·ROS·DBF의 PD 8구획 (총 18 ODE)
- **주요 치료 시나리오**: 무치료 일차성, 니페디핀 30mg QD, 실데나필 50mg BID, 보센탄 125mg BID(이차성), 일로프로스트 IV 5일(이차성), 프라조신 1mg BID, 니페디핀+실데나필 병용, 무치료 이차성(SSc), 한랭 자극 시험
- **보정/근거**: Thompson & Pope *Rheumatology* 2005(CCB 메타분석), Fries 등 *Circulation* 2005(실데나필 RCT), Matucci-Cerinic 등 *Ann Rheum Dis* 2011(보센탄 RAPIDS-2), Belch 등 *Ann Rheum Dis* 1995(일로프로스트 RCT) 기반 파라미터 설정

## Shiny 대시보드 (Dashboard)

7개 탭 구성: (1) **환자 프로파일** — 아형·α2-AR 감수성·유발 인자·중증도 설정; (2) **약물 PK** — 5개 약물 혈중 농도 시간 경과; (3) **혈관활성 매개인자** — ET-1·NE·ROS·RhoA 동태; (4) **혈관운동 반응** — 디지털 혈류·한랭 자극 시험·발작 빈도; (5) **임상 엔드포인트** — RCS·혈관경련 빈도·VAS·디지털 궤양 위험; (6) **시나리오 비교** — 9개 치료 간 비교 분석; (7) **바이오마커** — 신호 매개인자 상관성·모세혈관경 지수.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("raynaud_mrgsolve_model.R")
out <- mrgsim(mod, end = 2016, delta = 2)
plot(out)
# Shiny 대시보드:
shiny::runApp("raynaud_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg raynaud_qsp_model.dot -o raynaud_qsp_model.svg
dot -Tpng -Gdpi=150 raynaud_qsp_model.dot -o raynaud_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [raynaud_references.md](raynaud_references.md) 참조 (62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
