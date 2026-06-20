# 면역혈소판감소자반증 (ITP) (Immune Thrombocytopenic Purpura, ITP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈액

[![ITP QSP Model](itp_qsp_model.png)](itp_qsp_model.svg)

## 개요 (Overview)
면역혈소판감소자반증(ITP)은 항혈소판 자가항체에 의한 혈소판 파괴 촉진 및 거핵구-혈소판 생성 억제로 발생하는 자가면역 출혈 질환입니다. 유병률은 성인 약 10만 명당 9.5명이며, 발생은 소아와 60대 이상 여성에서 이중 정점을 보입니다. 항혈소판 항체(주로 항GP IIb/IIIa, Ib/IX)가 혈소판을 opsonize하여 비장·간 대식세포에 의한 파괴를 촉진하며, T세포 면역조절 이상(Th17 증가/Treg 감소), 세포독성 T세포에 의한 직접 혈소판 파괴도 기전에 포함됩니다. 치료는 1차(스테로이드·IVIG), 2차(TPO 수용체 작용제·리툭시맙), 3차(비장 절제·FcRn 억제제·Syk 억제제)로 단계적으로 접근합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항혈소판 항체 생성 | B세포 활성화(Th2/Th17 지원), 항GP IIb/IIIa IgG 생성 | 혈소판 opsonization |
| 비장 파괴 | FcγR+ 대식세포에 의한 opsonized 혈소판 탐식 | 혈소판 수 감소 |
| 거핵구 억제 | 항혈소판 항체가 거핵구 성숙·혈소판 방출 억제 | 생성 저하 |
| Treg/Th17 불균형 | Treg 감소, Th17 증가 → B세포 조절 이상 | 자가항체 지속 생성 |
| TPO 피드백 파괴 | 혈소판 감소 시 TPO 증가 불충분(비장에서 제거) | 생성 보상 실패 |
| 세포독성 T세포 | 혈소판 직접 파괴 | 항체 음성 ITP 기전 |

## 주요 약물 표적 (Drug Targets)
- **코르티코스테로이드** (프레드니솔론): FcγR 발현 억제, 항체 생성 억제, 비장 파괴 감소 (1차)
- **IVIG / 항D 면역글로불린**: FcγR 포화·차단, 혈소판 파괴 일시 억제 (응급)
- **TPO 수용체 작용제** (로미플로스팀, 엘트롬보파그): 거핵구 증식·혈소판 생성 촉진 (2차)
- **리툭시맙** (항CD20): B세포 고갈, 자가항체 생성 억제 (2차)
- **포스타마티닙 (Syk 억제제)**: 대식세포 FcγR 신호 억제, 비장 파괴 차단 (3차)
- **에파가르티지모드 / 니포칼리맙 (FcRn 억제제)**: IgG 반감기 단축 → 항체 소모 (3차)

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [itp_qsp_model.dot](itp_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 110 노드 / 10 클러스터) |
| [itp_qsp_model.svg](itp_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [itp_qsp_model.png](itp_qsp_model.png) | PNG 이미지 (150 dpi) |
| [itp_mrgsolve_model.R](itp_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 6개 치료 시나리오) |
| [itp_shiny_app.R](itp_shiny_app.R) | Shiny 대시보드 |
| [itp_references.md](itp_references.md) | 참고문헌 (약 51편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(프레드니솔론 1구획, IVIG 1구획, 로미플로스팀 2구획, 리툭시맙 2구획, 포스타마티닙 2구획, 에파가르티지모드 1구획, 엘트롬보파그 1구획) + 질환 PD 구획(혈소판·비장 혈소판, TPO, 거핵구 전구세포·거핵구, 자가항체, B세포, Treg, Th17, 대식세포)
- **주요 치료 시나리오**: ① 자연경과, ② 스테로이드 단독, ③ IVIG + 스테로이드, ④ 로미플로스팀 단독, ⑤ 엘트롬보파그 단독, ⑥ 리툭시맙 단독, ⑦ 포스타마티닙, ⑧ 에파가르티지모드
- **보정/근거**: EXTEND 시험(엘트롬보파그), 로미플로스팀 3상 시험, ASH 2019 ITP 가이드라인 데이터 기반 파라미터 설정

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(기저 혈소판 수·출혈 중증도·이전 치료력), 약물 PK 농도 추이, 혈소판·거핵구·자가항체 PD 바이오마커, 임상 엔드포인트(혈소판 수 ≥30 × 10⁹/L 달성·출혈 위험), 치료 시나리오 비교, 장기 관해 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("itp_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("itp_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg itp_qsp_model.dot -o itp_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [itp_references.md](itp_references.md) 참조 (약 51편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
