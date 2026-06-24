# 미세변화 신증후군 (Minimal Change Disease, MCD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![MCD QSP Model](mcd_qsp_model.png)](mcd_qsp_model.svg)

## 개요 (Overview)
미세변화 신증후군(MCD)은 광학현미경 상 정상 소견이나 전자현미경에서 족세포 족돌기 소실이 특징적인 사구체 질환으로, 소아 신증후군의 약 70~80%, 성인의 약 15~25%를 차지합니다. 특히 소아에서 연간 발생률은 10만 명당 약 2~7명입니다. 핵심 병인은 T세포 면역조절 이상으로 분비되는 순환 투과 인자(permeability factor)가 족세포의 nephrin·podocin 등 슬릿막 단백을 손상시켜 대량 단백뇨를 유발하는 것으로 이해되며, 최근 anti-nephrin 항체도 일부 환자에서 발견되어 B세포 역할도 주목받고 있습니다. 스테로이드(프레드니솔론)에 대부분 반응하나 약 50%에서 재발하며, 스테로이드 의존성·빈발 재발 시 칼시뉴린 억제제(타크로리무스·사이클로스포린), 리툭시맙을 사용합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| T세포 면역조절 이상 | CD4+ Th2 우세·Treg 감소, IL-13·IL-4 과잉 생성 | 순환 투과 인자 방출 |
| 순환 투과 인자 | 특성 미확정(hemopexin·CLCF-1·CD80?), 족세포 손상 유도 | 슬릿막 단백 손상·족돌기 소실 |
| anti-Nephrin 항체 | IgG 항nephrin 자가항체(일부 환자) → 직접 슬릿막 손상 | 급성 대량 단백뇨 |
| 족세포 슬릿막 손상 | nephrin·podocin·synaptopodin 재배열·소실 | 선택적 단백뇨(알부민 > IgG) |
| 저알부민혈증 | 대량 알부민 소실 → 혈청 알부민 <2.5 g/dL | 부종·혈전·저칼슘혈증 |
| 이차 고지혈증 | LDL 생성 증가·HDL 소실로 총 콜레스테롤 상승 | 심혈관 위험·지방요 |
| 세관 손상 | 과부하 단백뇨에 의한 세관 독성 | 드물게 ATN 동반 |

## 주요 약물 표적 (Drug Targets)
- **프레드니솔론** (고용량 스테로이드): T세포·B세포 억제, NF-κB 차단 → 투과 인자 생성 억제 (1차, 소아 2 mg/kg/일)
- **타크로리무스** (칼시뉴린 억제제): IL-2 억제 → T세포 활성화 차단, 족세포 직접 보호 효과; 스테로이드 의존성/빈발 재발
- **사이클로스포린** (칼시뉴린 억제제): 타크로리무스와 유사 기전, 2차 선택
- **리툭시맙** (항CD20): B세포·T세포 조절, 빈발 재발·스테로이드 의존성 MCD에 점증적 사용
- **마이코페놀레이트 모페틸(MMF)**: 스테로이드 절약 목적, 성인 스테로이드 의존성
- **사이클로포스파마이드**: 장기 관해 유도, 빈발 재발 소아에 역사적으로 사용

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [mcd_qsp_model.dot](mcd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 216 노드 / 17 클러스터) |
| [mcd_qsp_model.svg](mcd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mcd_qsp_model.png](mcd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [mcd_mrgsolve_model.R](mcd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [mcd_shiny_app.R](mcd_shiny_app.R) | Shiny 대시보드 |
| [mcd_references.md](mcd_references.md) | 참고문헌 (약 46편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(프레드니솔론 3구획, 사이클로스포린 2구획, 타크로리무스 2구획, 리툭시맙 2구획) + 질환 PD 구획(CD4 효과 T세포, Treg, B세포, 투과 인자, anti-nephrin 항체, 족세포 손상, 슬릿막 손상, 단백뇨, 혈청 알부민, 부종, 콜레스테롤, eGFR)
- **주요 치료 시나리오**: ① 자연경과(무치료), ② 프레드니솔론 고용량 유도, ③ 프레드니솔론 + 사이클로스포린(빈발 재발), ④ 프레드니솔론 + 타크로리무스, ⑤ 리툭시맙 + 저용량 스테로이드, ⑥ MMF + 저용량 스테로이드(스테로이드 절약)
- **보정/근거**: KDIGO 2021 신증후군 가이드라인, MENTOR MCD 코호트, ICON 소아 신증후군 연구 데이터 기반 파라미터 설정

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(나이·기저 단백뇨·혈청 알부민·재발 이력), 약물 PK 농도 추이(스테로이드·칼시뉴린 억제제·리툭시맙), 주요 PD 바이오마커(족세포 손상·CD20 B세포·T세포), 임상 엔드포인트(단백뇨 관해·재발 횟수·eGFR 변화), 치료 시나리오 비교, 재발 위험 및 스테로이드 누적 용량 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("mcd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("mcd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg mcd_qsp_model.dot -o mcd_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [mcd_references.md](mcd_references.md) 참조 (약 46편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
