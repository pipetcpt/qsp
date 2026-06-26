# IgA 혈관염 (HSP) (IgA Vasculitis, IgAV) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![IgAV QSP Model](igav_qsp_model.png)](igav_qsp_model.svg)

## 개요 (Overview)
IgA 혈관염(IgAV, 구 헤노흐-쇤라인 자반증/HSP)은 IgA1 면역복합체가 소혈관 벽에 침착되어 발생하는 전신성 혈관염입니다. 소아에서 가장 흔한 전신 혈관염으로, 연간 발생률은 소아 10만 명당 약 10~22명이며 성인에서는 드물지만 더 중증으로 경과합니다. 촉지성 자반(100%), 관절통(75%), 복통(65%), 신염(40~50%)의 4대 임상 증상을 특징으로 합니다. Gd-IgA1 기반 면역복합체 형성, 보체(대체경로) 활성화, 내피세포 손상, 호중구·마크로파지 유입이 핵심 병인이며, 치료는 스테로이드 및 중증 신침범 시 면역억제제(MMF, 사이클로포스파마이드, 리툭시맙)를 사용합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Gd-IgA1 면역복합체 형성 | 갈락토스 결핍 IgA1 + 자가항체(IgG) → IC 형성 | 혈관벽 침착 |
| 보체 대체경로 활성화 | C3, C5a 생성, MAC 형성 | 내피세포 손상·혈관 투과성 증가 |
| 호중구 침윤 | C5a·IL-8에 의한 호중구 유입, NETs 형성 | 백혈구파괴성 혈관염 |
| 피부 소혈관염 | IgA·C3 피부 침착, 진피 소혈관 손상 | 촉지성 자반 |
| 신장 침범 | 메산지움 IgA 침착, ANCA 유사 반응 | 혈뇨·단백뇨·신부전 |
| 장관 침범 | 장 소혈관 염증·출혈 | 복통·장중첩·혈변 |
| 관절 침범 | 활막 소혈관 염증 | 관절통·관절염 |

## 주요 약물 표적 (Drug Targets)
- **코르티코스테로이드** (프레드니솔론): 광범위 면역억제, NF-κB 억제, 장·관절 증상 완화
- **마이코페놀레이트 모페틸(MMF)**: B세포 증식 억제 → 자가항체·IgA 생성 감소, 중증 신염
- **리툭시맙** (항CD20): 중증 또는 재발 신염에 B세포 고갈
- **사이클로포스파마이드**: 급속진행성 사구체신염 동반 중증 신침범
- **다프손**: 피부 혈관염에 보조 사용
- **ACEi/ARB**: 단백뇨·혈압 관리

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [igav_qsp_model.dot](igav_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 200 노드 / 10 클러스터) |
| [igav_qsp_model.svg](igav_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [igav_qsp_model.png](igav_qsp_model.png) | PNG 이미지 (150 dpi) |
| [igav_mrgsolve_model.R](igav_mrgsolve_model.R) | mrgsolve ODE 모델 (약 28 구획 / 5개 치료 시나리오) |
| [igav_shiny_app.R](igav_shiny_app.R) | Shiny 대시보드 |
| [igav_references.md](igav_references.md) | 참고문헌 (약 70편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(프레드니솔론 3구획, MMF 3구획, 리툭시맙 3구획 + CD20 결합) + 질환 PD 구획(Gd-IgA1, 자가항체, 면역복합체, 보체 C3/C5, 메산지움, 신장 염증, 단백뇨, eGFR, 피부 자반, 소화기 염증, BAFF, IL-6, TNF, 섬유화)
- **주요 치료 시나리오**: ① 자연경과(무치료), ② 스테로이드 단독(저~중용량), ③ MMF + 스테로이드, ④ 리툭시맙 단독(중증 신염), ⑤ 스테로이드 + MMF + ACEi 병용
- **보정/근거**: SHARE 이니셔티브 권고안, KDIGO 2021 소아 혈관염 지침, 단일기관 코호트 연구 데이터를 파라미터 보정 기준으로 참조

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(나이·체중·신장 침범 여부), 약물 PK 농도 추이(스테로이드·MMF·리툭시맙), 주요 PD 바이오마커(IgA·보체·CD20 B세포), 임상 엔드포인트(자반 회복·단백뇨·eGFR), 치료 시나리오 비교, 장기 신기능 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("igav_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("igav_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg igav_qsp_model.dot -o igav_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [igav_references.md](igav_references.md) 참조 (약 70편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
