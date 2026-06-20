# 막성 신병증 (Membranous Nephropathy, MN) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![MN QSP Model](mn_qsp_model.png)](mn_qsp_model.svg)

## 개요 (Overview)
막성 신병증(MN)은 성인 원발성 신증후군의 가장 흔한 원인 중 하나로, 전 세계 사구체 신염의 약 20~37%를 차지합니다. 원발성 MN의 약 70~80%에서 족세포 표면 항원인 M형 포스포리파제 A2 수용체(PLA2R1)에 대한 자가항체(항PLA2R 항체)가 발견됩니다. 항PLA2R 항체-항원 면역복합체가 사구체 기저막 상피하에 침착되어 보체(막공격복합체, MAC) 형성 → 족세포 손상 → 사구체 기저막 비후 → 대량 단백뇨(신증후군)로 이어지는 기전이 핵심입니다. 리툭시맙(항CD20) B세포 고갈이 현재 원발성 MN의 표준 치료로 자리잡았으며, 칼시뉴린 억제제(타크로리무스·사이클로스포린), 사이클로포스파마이드(폰티첼리 요법)도 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항PLA2R 항체 생성 | B세포·형질세포의 항PLA2R1 IgG4 생성 (HLA-DQA1 연관) | 혈청 항PLA2R 항체 양성 |
| 상피하 면역복합체 침착 | 항PLA2R IgG가 족세포 표면 PLA2R1에 결합 → 상피하 IC 형성 | GBM 비후·스파이크 형성 |
| 보체 활성화 | 렉틴경로/고전경로 → C3·C5b-9(MAC) 생성 | 족세포 직접 손상 |
| 족세포 손상 | MAC에 의한 산화 스트레스·액틴 재배열·슬릿막 단백 소실 | 단백뇨 증가 |
| 사구체 기저막 변화 | 족세포 손상 → GBM 재형성 이상 | GBM 비후·상피하 침착 |
| 알부민 소실 | 대량 단백뇨 → 혈청 알부민 감소 | 부종·혈전 위험 |
| 레닌-안지오텐신 활성화 | 저알부민혈증 → RAAS 활성 → 수분·나트륨 저류 | 부종 악화·신기능 저하 |

## 주요 약물 표적 (Drug Targets)
- **리툭시맙** (항CD20): B세포 고갈 → 항PLA2R 항체 생성 억제, GEMRITUX·MENTOR 시험에서 효과 입증
- **타크로리무스** (칼시뉴린 억제제): T세포 의존적 B세포 활성 억제, 족세포 직접 보호
- **사이클로포스파마이드** (알킬화제): 폰티첼리 요법(사이클로포스파마이드 + 스테로이드 교대), 중증 MN
- **ACEi/ARB**: 단백뇨 감소·사구체 고혈압 완화, 모든 환자 1차 보존 치료
- **오파투무맙 / 오비누투주맙** (차세대 항CD20): 리툭시맙 불응 시 대안
- **벨리무맙** (항BAFF): 형질세포 생존 신호 억제 연구 중

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [mn_qsp_model.dot](mn_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 212 노드 / 13 클러스터) |
| [mn_qsp_model.svg](mn_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mn_qsp_model.png](mn_qsp_model.png) | PNG 이미지 (150 dpi) |
| [mn_mrgsolve_model.R](mn_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [mn_shiny_app.R](mn_shiny_app.R) | Shiny 대시보드 |
| [mn_references.md](mn_references.md) | 참고문헌 (약 53편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(리툭시맙 2구획 + CD20 결합, 타크로리무스 2구획, 사이클로포스파마이드 2구획 + 대사산물) + 질환 PD 구획(CD20 B세포·형질세포, 항PLA2R1 항체, IgG 침착, 보체 MAC, 족세포 손상, GBM 비후, 단백뇨, 혈청 알부민, eGFR, 안지오텐신 II, 알도스테론)
- **주요 치료 시나리오**: ① 자연경과(자발적 관해 가능성 포함), ② 리툭시맙 단독(375 mg/m² × 1~2회), ③ 타크로리무스 + ACEi, ④ 폰티첼리 요법(사이클로포스파마이드 + 스테로이드), ⑤ 리툭시맙 반복 투여(B세포 재증식 시), ⑥ 타크로리무스 + 리툭시맙 병용
- **보정/근거**: MENTOR 시험(리툭시맙 vs 사이클로스포린), GEMRITUX 시험, RI-CYCLO 시험 데이터 기반 파라미터 설정

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(기저 항PLA2R 항체 역가·단백뇨·eGFR·혈청 알부민), 약물 PK 농도 추이(리툭시맙·타크로리무스), 주요 PD 바이오마커(CD20 B세포·항PLA2R 항체·보체), 임상 엔드포인트(단백뇨 완전/부분 관해·eGFR 변화), 치료 시나리오 비교, 면역억제 강도별 반응 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("mn_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("mn_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg mn_qsp_model.dot -o mn_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [mn_references.md](mn_references.md) 참조 (약 53편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
