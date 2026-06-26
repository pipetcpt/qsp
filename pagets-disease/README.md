# 파젯병 (골) (Paget's Disease of Bone, PBD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![PBD QSP Model](pbd_qsp_model.png)](pbd_qsp_model.svg)

## 개요 (Overview)

골 파젯병(PBD)은 특정 부위의 골 개조(bone remodeling)가 비정상적으로 항진되는 만성 대사 골질환으로, 50세 이상 북유럽계에서 약 2–3%의 유병률을 보인다. SQSTM1(p62) 유전자 변이와 파라믹소바이러스 감염이 핵심 위험 인자로, 과활성화된 파골세포가 비정상 골흡수를 주도하고 이에 반응한 조골세포의 과도한 골형성이 혼재한다. 결과적으로 무질서하게 재구성된 직물뼈(woven bone)가 형성되어 골 변형(두개골 비대·장골 만곡)·통증·청력 소실·신경 압박·고심박출 심부전 등을 유발한다. 졸레드론산(단회 5mg IV)이 현재 최강의 1차 치료제이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| SQSTM1/NF-κB 과활성 | SQSTM1 변이 → p62 응집 → TRAF6/RANK 신호 과잉 → NF-κB↑ | 파골세포 과활성, 비정상 골흡수 |
| RANKL 과잉 발현 | 골세포·조골세포 RANKL↑ → 파골세포 분화·활성↑ | 골흡수 표지자(CTX, NTX) 상승 |
| 파골세포 수 증가 | 거핵세포 유사 과형성 파골세포 → 골 과흡수 | bsALP·NTX 수배 상승 |
| 반응적 조골세포 활성 | 흡수 후 과도한 골형성 → 직물뼈 생성 | bsALP↑, 비정상 골 구조 |
| 혈관 과형성 | 파젯 병변 내 혈관화↑ → 고심박출 증가 | 소수에서 고박출성 심부전 |
| 신경·혈관 압박 | 비대해진 뼈 → 청신경공·척추공 협착 | 청력 소실, 신경근 압박 |
| 통증 기전 | 골막 신장·미세골절·신경 압박 | 국소 통증, 관절 이차 변성 |

## 주요 약물 표적 (Drug Targets)

- **졸레드론산 (Zoledronic acid)**: 3세대 N-비스포스포네이트; 파골세포 메발로네이트 경로(FPP 합성효소) 억제 → 파골세포 아폽토시스 (5mg IV 단회, 수년간 관해 유지; HORIZON-PDB)
- **알렌드로네이트 (Alendronate)**: 경구 비스포스포네이트; 골흡수 60–70% 억제 (40mg/일 6개월)
- **칼시토닌 (Calcitonin SC)**: 파골세포 억제; 역사적 표준 치료, 졸레드론산 등장 후 2차 대안
- **데노수맙 (Denosumab)**: 항RANKL → 파골세포 분화 차단; PBD에서의 역할 연구 중
- **지지 치료**: 적절한 칼슘·비타민D 보충; 수술(변형 교정·청력 보조기)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [pbd_qsp_model.dot](pbd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 157 노드 / 12 클러스터) |
| [pbd_qsp_model.svg](pbd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pbd_qsp_model.png](pbd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pbd_mrgsolve_model.R](pbd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 7개 치료 시나리오) |
| [pbd_shiny_app.R](pbd_shiny_app.R) | Shiny 대시보드 |
| [pbd_references.md](pbd_references.md) | 참고문헌 (약 53편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 졸레드론산 IV 2구획+골 결합, 알렌드로네이트 경구 흡수·중심·말초, 칼시토닌 SC 2구획, 데노수맙 SC+중심+말초; RANKL·OPG, 파골세포 전구체·파골세포·조골세포 전구체·조골세포, BMD, bsALP·NTX·CTX(골전환 표지자), 통증 점수
- **주요 치료 시나리오**: ① 무치료 ② 졸레드론산 5mg IV 단회 ③ 알렌드로네이트 40mg/일 ×6개월 ④ 칼시토닌 100IU SC 매일 ×6개월 ⑤ 데노수맙 60mg SC Q6M ⑥ 졸레드론산+지지 치료 ⑦ 알렌드로네이트→졸레드론산 순차요법
- **보정/근거**: HORIZON-PDB(졸레드론산 단회 5mg vs 리세드로네이트), PRISM(파젯병 경과 추적) 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (이환 부위 수·기저 ALP·통증 점수·SQSTM1 변이 유무 설정), ② **PK** (혈장 및 골 내 약물 농도-시간 곡선), ③ **PD 주요지표** (파골세포·조골세포 활성, RANKL/OPG 추이), ④ **임상 엔드포인트** (ALP/bsALP 정상화, NTX/CTX 감소, 통증 개선), ⑤ **시나리오 비교** (7개 치료 전략 직접 비교), ⑥ **바이오마커** (BMD·골전환 마커 장기 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("pbd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pbd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pbd_qsp_model.dot -o pbd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [pbd_references.md](pbd_references.md) 참조 (약 53편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
