# 현미경적 다발혈관염 (MPA) (Microscopic Polyangiitis, MPA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![MPA QSP Model](mpa_qsp_model.png)](mpa_qsp_model.svg)

## 개요 (Overview)
현미경적 다발혈관염(MPA)은 MPO-ANCA(항미엘로페록시다제 항중성구 세포질 항체) 양성이 특징적인 소혈관 괴사성 혈관염입니다. 연간 발생률은 100만 명당 약 3~7명으로 드문 질환이지만, 진단 지연 시 급속진행성 사구체신염(RPGN)·폐모세혈관염(폐출혈)으로 이어지는 생명 위협적 상태를 초래합니다. MPO-ANCA가 과활성화된 호중구를 직접 자극해 내피세포에 달라붙어 산화 버스트·세포독성 효소를 분비함으로써 괴사성 혈관염을 유발합니다. 보체 대체경로(C5a-C5aR1 축)도 호중구 활성화를 증폭시킵니다. 치료는 유도(리툭시맙 또는 사이클로포스파마이드 + 스테로이드)와 유지(리툭시맙 또는 아자티오프린) 단계로 나뉘며, 아바코판(C5aR1 억제제)이 스테로이드 절약 전략으로 2021년 승인되었습니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| MPO-ANCA 생성 | B세포 이상 활성화, CD4 Th17 증가, Treg 감소 | 혈청 MPO-ANCA 양성 |
| 호중구 과활성화 | MPO-ANCA가 primed 호중구의 FcγRIIa·PR3 자극 | 산화 버스트·NET 형성·내피 손상 |
| 보체 C5a-C5aR1 축 | 보체 대체경로 활성 → C5a 생성 → 호중구 C5aR1 자극 | 호중구 모집·활성화 피드백 증폭 |
| 신장 괴사성 사구체신염 | 사구체 모세혈관 내피 손상·반월체 형성 | 급속진행성 신부전 |
| 폐 모세혈관염 | 폐 모세혈관 내피 손상, 적혈구 유출 | 폐포 출혈·DAH |
| 전신 염증 | CRP·ESR 상승, 발열·체중감소 | BVAS 점수 반영 |
| 조직 섬유화 | 반월체의 섬유화 전환, TGF-β 매개 | 만성 신기능 저하 |

## 주요 약물 표적 (Drug Targets)
- **리툭시맙** (항CD20): B세포 고갈 → MPO-ANCA 생성 억제, RAVE/RITUXVAS 시험에서 유도·유지 효과 입증
- **사이클로포스파마이드** (알킬화제): B·T세포 고갈, 중증 폐·신 침범 유도 요법
- **프레드니솔론** (고용량 스테로이드): 급성 염증 억제, 유도 기간 병용
- **아바코판 (Avacopan, C5aR1 억제제)**: 보체 C5a-C5aR1 축 차단 → 호중구 활성화 억제, 스테로이드 절약; ADVOCATE 시험에서 완전 관해 비열등성 입증
- **아자티오프린 / 메토트렉세이트**: 완화 후 유지 면역억제
- **벨리무맙 + 리툭시맙**: BAFF 차단 병용 전략 연구 중

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [mpa_qsp_model.dot](mpa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 156 노드 / 11 클러스터) |
| [mpa_qsp_model.svg](mpa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mpa_qsp_model.png](mpa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [mpa_mrgsolve_model.R](mpa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 5개 치료 시나리오) |
| [mpa_shiny_app.R](mpa_shiny_app.R) | Shiny 대시보드 |
| [mpa_references.md](mpa_references.md) | 참고문헌 (약 62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(사이클로포스파마이드 3구획 + 4-OH-CY 대사산물, 리툭시맙 2구획, 프레드니솔론 2구획, 아바코판 2구획) + 질환 PD 구획(B세포, 형질세포, MPO-ANCA, C5a, 활성화 호중구, 내피 손상, 신장 염증/섬유화, GFR, 폐 염증, CRP)
- **주요 치료 시나리오**: ① 자연경과, ② 사이클로포스파마이드 + 스테로이드(유도), ③ 리툭시맙 + 스테로이드(유도), ④ 아바코판 + 리툭시맙(스테로이드 절약), ⑤ 유지: 리툭시맙 6개월 간격
- **보정/근거**: RAVE(리툭시맙 vs CY), ADVOCATE(아바코판), MEPEX(혈장교환) 시험 데이터 기반 파라미터 설정

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(MPO-ANCA 역가·기저 eGFR·BVAS·폐 침범 유무), 약물 PK 농도 추이, 주요 PD 바이오마커(B세포·MPO-ANCA·C5a·CRP), 임상 엔드포인트(BVAS·eGFR·단백뇨 변화), 치료 시나리오 비교, 관해 유지 및 재발 위험 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("mpa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("mpa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg mpa_qsp_model.dot -o mpa_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [mpa_references.md](mpa_references.md) 참조 (약 62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
