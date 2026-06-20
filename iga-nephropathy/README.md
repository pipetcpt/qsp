# IgA 신병증 (IgA Nephropathy, IgAN) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![IgAN QSP Model](igan_qsp_model.png)](igan_qsp_model.svg)

## 개요 (Overview)
IgA 신병증(IgAN)은 전 세계에서 가장 흔한 원발성 사구체 신염으로, 성인 원발성 사구체 신염의 약 30~40%를 차지합니다. 핵심 발병기전은 소위 "4단계 타격(four-hit hypothesis)"으로 설명됩니다: 비정상적으로 갈락토스화가 부족한 IgA1(Gd-IgA1) 과잉 생성, Gd-IgA1에 대한 자가항체(IgG/IgA) 형성, 면역복합체 형성, 메산지움 침착 및 보체·메산지움 세포 활성화에 의한 신손상입니다. 치료는 스파르센탄(SPARTAN 시험), 타르게탄(타겟팅 BAFF/APRIL), 이프타코판(보체 C3 억제), 부데소니드 장용 제제(네페콘), 그리고 전통적인 스테로이드·ACEi/ARB까지 다양합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Gd-IgA1 과잉 생성 | 점막 B세포의 갈락토실트랜스퍼라제(GALNT2/C1GALT1) 이상 | 혈중 Gd-IgA1 상승 |
| 자가항체 형성 | 항Gd-IgA1 IgG 항체, BAFF·APRIL 증가 | 면역복합체(IC) 형성 |
| 메산지움 침착 | IC의 메산지움 수용체(CD71 등) 결합 | 메산지움 세포 활성화·증식 |
| 보체 활성화 | 대체경로/렉틴경로 활성, C3·C5b-9 생성 | 족세포 손상·GFR 저하 |
| 족세포 손상 | 보체·사이토카인(IL-6, TNF-α)에 의한 슬릿막 손상 | 단백뇨 증가 |
| 세관간질 섬유화 | TGF-β 매개 근섬유모세포 활성화 | eGFR 점진적 저하 |
| 혈압·레닌-안지오텐신계 | 단백뇨·사구체 고혈압 | 신기능 악화 가속 |

## 주요 약물 표적 (Drug Targets)
- **RAAS 억제제** (ACEi/ARB): 사구체 고혈압 완화, 단백뇨 감소
- **부데소니드 장용 제제 (네페콘)**: 점막 IgA 생성 억제 → Gd-IgA1 감소
- **스파르센탄** (이중 안지오텐신/엔도텔린 길항제): 단백뇨 감소
- **이프타코판 (Iptacopan)**: 보체 인자 B(CFB) 억제 → 대체경로 차단
- **시베팔리맙 / 폴로카나맙** (APRIL/BAFF 억제): 자가항체 생성 억제
- **리툭시맙** (항CD20): B세포 고갈, 중증 단백뇨 환자 선별 사용

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [igan_qsp_model.dot](igan_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 169 노드 / 12 클러스터) |
| [igan_qsp_model.svg](igan_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [igan_qsp_model.png](igan_qsp_model.png) | PNG 이미지 (150 dpi) |
| [igan_mrgsolve_model.R](igan_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 5개 치료 시나리오) |
| [igan_shiny_app.R](igan_shiny_app.R) | Shiny 대시보드 |
| [igan_references.md](igan_references.md) | 참고문헌 (약 46편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(부데소니드 2구획, 스파르센탄 2구획, 이프타코판 2구획, 시베팔리맙 3구획) + 질환 PD 구획(Gd-IgA1, 자가항체, 면역복합체, 보체, 메산지움, 족세포, 세관간질섬유화, UPCR, eGFR, 혈압)
- **주요 치료 시나리오**: ① 무치료 자연경과, ② 부데소니드(네페콘) 단독, ③ 스파르센탄 단독, ④ 이프타코판 + ACEi, ⑤ 시베팔리맙 + 스테로이드 병용
- **보정/근거**: PROTECT 시험(스파르센탄), APPLAUSE-IgAN 시험(이프타코판), NefIgArd 시험(부데소니드 장용 제제) 데이터를 파라미터 보정 기준으로 참조

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(기저 eGFR·UPCR·혈압), 약물 PK 농도 추이, 주요 PD 바이오마커(Gd-IgA1·자가항체·보체), 임상 엔드포인트(단백뇨·eGFR 변화), 치료 시나리오 비교, 장기 신기능 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("igan_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("igan_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg igan_qsp_model.dot -o igan_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [igan_references.md](igan_references.md) 참조 (약 46편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
