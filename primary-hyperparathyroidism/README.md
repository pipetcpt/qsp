# 원발성 부갑상선 기능 항진증 (PHPT) (Primary Hyperparathyroidism, PHPT) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![PHPT QSP Model](phpt_qsp_model.png)](phpt_qsp_model.svg)

## 개요 (Overview)
원발성 부갑상선 기능 항진증(PHPT)은 부갑상선 선종(85%)·증식증·드물게 암에 의한 PTH 자율 과다 분비로 만성 고칼슘혈증이 유발되는 내분비 질환이다. 유병률은 인구 1,000명당 1~2명이며, 50~60대 여성에서 가장 흔하다. 칼슘 감지 수용체(CaSR)의 PTH-칼슘 음성 피드백이 결손되어 PTH가 지속 분비되고, PTH는 신장 칼슘 재흡수·1,25-비타민D 활성화·파골세포 자극(RANK-RANKL 경로)을 통해 골소실(특히 피질골)과 신석회화·신결석을 유발한다. 근치는 부갑상선절제술이며, 수술 불가 또는 무증상 PHPT에는 CaSR 양성조절제(시나칼셋)와 데노수맙이 약물 치료로 사용된다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| PTH-CaSR 피드백 이탈 경로 | CaSR 민감도 저하·선종 자율 분비 | PTH 과다·고칼슘혈증 |
| PTH-골 경로 | RANKL 증가 → 파골세포 활성화 → 피질골 소실 | BMD 감소·골절 위험 |
| PTH-신장 경로 | 원위세뇨관 칼슘 재흡수 증가·FGF23 조절 | 고칼슘뇨·신결석·신기능 저하 |
| 비타민D 활성화 경로 | PTH → CYP27B1 활성 → 1,25-(OH)₂D 증가 | 장 칼슘 흡수 증가 |
| 뼈 재형성 불균형 경로 | 조골세포(OB)·파골세포(OC) 분리 활성 | 섬유성 골염(심한 경우) |
| FGF23-Klotho 축 | PTH에 의한 FGF23 간접 조절 | 인산 조절 이상 |

## 주요 약물 표적 (Drug Targets)
- **시나칼셋 (calcimimetic)**: CaSR 양성조절제 → PTH 분비 억제·혈청 칼슘 감소 (FDA 승인, 수술 불가 환자)
- **데노수맙**: 항RANKL 단클론항체 → 파골세포 활성 억제 → BMD 보전
- **비스포스포네이트 (알렌드로네이트 등)**: 파골세포 아포토시스 유도 → 피질골 보호
- **부갑상선절제술**: PTH 과분비 선종 제거 — 유일한 근치적 치료
- **수분 공급 및 루프이뇨제**: 고칼슘혈증 급성기 대증 (입원 설정)

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [phpt_qsp_model.dot](phpt_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 128 노드 / 12 클러스터) |
| [phpt_qsp_model.svg](phpt_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [phpt_qsp_model.png](phpt_qsp_model.png) | PNG 이미지 (150 dpi) |
| [phpt_mrgsolve_model.R](phpt_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 8 치료 시나리오) |
| [phpt_shiny_app.R](phpt_shiny_app.R) | Shiny 대시보드 |
| [phpt_references.md](phpt_references.md) | 참고문헌 (약 62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 시나칼셋 1구획·데노수맙 2구획+RANKL 결합·알렌드로네이트(골 결합) PK 구획, PTH·혈청Ca·혈청PO4·비타민D(25·1,25)·조골세포·파골세포·RANKL·요중Ca·GFR·요추BMD·대퇴경부BMD PD 구획
- **주요 치료 시나리오**: ① 정상(건강), ② 무치료 경증 PHPT, ③ 무치료 중증 PHPT, ④ 시나칼셋 60 mg/일, ⑤ 데노수맙 60 mg q6mo, ⑥ 부갑상선절제술(90일째), ⑦ 시나칼셋+데노수맙, ⑧ CKD 동반 PHPT + 시나칼셋 90 mg
- **보정/근거**: 국제 PHPT 가이드라인(5th International Workshop, 2022)과 시나칼셋 임상시험(SHOPPE, PRIMARY 등) 데이터를 참고하여 혈청 Ca·PTH 정상화 시간 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(증상 유무, 골다공증 등급, 신결석 이력) · 시나칼셋/데노수맙 PK 프로파일 · PTH-칼슘-비타민D PD 축 · BMD 및 골절 위험 임상 엔드포인트 · 치료 시나리오 비교(혈청 Ca·BMD 변화) · 바이오마커(PTH·24h 요칼슘·GFR) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("phpt_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("phpt_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg phpt_qsp_model.dot -o phpt_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [phpt_references.md](phpt_references.md) 참조 (약 62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
