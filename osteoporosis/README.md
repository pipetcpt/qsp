# 골다공증 (Osteoporosis, OP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![OP QSP Model](op_qsp_model.png)](op_qsp_model.svg)

## 개요 (Overview)

골다공증은 골량 감소(T-점수 ≤ −2.5)와 골 미세구조 이상으로 골절 위험이 높아지는 만성 대사 골질환으로, 50세 이상 여성의 약 30%, 남성의 약 12%에 영향을 미친다. 에스트로겐 결핍(폐경)은 OPG/RANKL 비율을 낮춰 파골세포 활성을 급격히 증가시키고, 골형성·골흡수의 불균형이 지속되어 BMD가 저하된다. 글루코코르티코이드 유발 골다공증(GIOP)과 부갑상선호르몬(PTH) 과잉 분비도 중요한 이차 원인이다. 비스포스포네이트(골흡수 억제), 데노수맙(항RANKL), 테리파라타이드·로모소주맙(골형성 촉진)이 현재 주요 치료제이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| RANKL-RANK-OPG 축 | 에스트로겐↓ → RANKL↑/OPG↓ → 파골세포 분화·활성 증가 | 골흡수 촉진, BMD 감소 |
| 에스트로겐 결핍 | 에스트로겐 수용체 신호↓ → 조골세포 아폽토시스↑ | 골형성 감소, 폐경 후 급속 골 소실 |
| PTH 신호 | 간헐적 PTH↑ → cAMP → Wnt 신호 강화 → 조골세포 증식 | 테리파라타이드에 의한 골형성 촉진 |
| Wnt/스클레로스틴 | 스클레로스틴(스클레로스틴-Wnt 길항) → LRP5/6 차단 → 조골세포 기능 억제 | 골형성 감소; 로모소주맙 표적 |
| 글루코코르티코이드 과잉 | GC → 조골세포 아폽토시스, RANKL↑ → 조골세포·파골세포 불균형 | GIOP (급성 골 소실) |
| 칼슘-비타민D 결핍 | PTH↑ → 신장 Ca 재흡수, 장 Ca 흡수↓ → 이차성 부갑상선 기능 항진 | 골연화증, 골절 위험 가중 |
| 골 전환 마커 | CTX(파골세포), P1NP(조골세포), bsALP | 치료 반응 모니터링 |

## 주요 약물 표적 (Drug Targets)

- **알렌드로네이트/졸레드로네이트 (비스포스포네이트)**: 파골세포 메발로네이트 경로 억제 → 파골세포 아폽토시스 → 골흡수 50–70% 감소
- **데노수맙 (Denosumab)**: 항RANKL 단클론항체 → 파골세포 분화 차단 (6개월 피하주사)
- **테리파라타이드 (Teriparatide, PTH1–34)**: 간헐적 PTH 작용 → Wnt 신호 → 조골세포 증식·활성; BMD 10% 이상 증가
- **로모소주맙 (Romosozumab)**: 항스클레로스틴 → Wnt 신호 탈억제 → 골형성↑/골흡수↓ 이중 효과 (ARCH, FRAME 임상시험)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [op_qsp_model.dot](op_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 189 노드 / 10 클러스터) |
| [op_qsp_model.svg](op_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [op_qsp_model.png](op_qsp_model.png) | PNG 이미지 (150 dpi) |
| [op_mrgsolve_model.R](op_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [op_shiny_app.R](op_shiny_app.R) | Shiny 대시보드 |
| [op_references.md](op_references.md) | 참고문헌 (약 54편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 알렌드로네이트 골 결합 1구획, 졸레드로네이트 골 결합 1구획, 데노수맙 SC+중심 2구획, 테리파라타이드 중심 1구획, 로모소주맙 SC+중심+말초 3구획; 에스트로겐·PTH·칼슘·RANKL·OPG·스클레로스틴·조골세포 전구체·조골세포·파골세포 전구체·파골세포·BMD·CTX·P1NP·bsALP·10년 골절 위험
- **주요 치료 시나리오**: S1 폐경 후 무치료, S2 알렌드로네이트 70mg/주, S3 졸레드로네이트 5mg/년, S4 데노수맙 60mg/6개월, S5 테리파라타이드 20μg/일, S6 로모소주맙→데노수맙 순차요법
- **보정/근거**: FIT(알렌드로네이트), HORIZON(졸레드로네이트), FREEDOM(데노수맙), ARCH(로모소주맙) 임상시험 BMD 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (폐경 여부·기저 BMD·T-점수·GIOP 여부 설정), ② **PK** (혈장/골 약물 농도-시간 곡선), ③ **PD 주요지표** (RANKL/OPG 비율, 조골세포·파골세포 활성 추이), ④ **임상 엔드포인트** (BMD g/cm², T-점수, 10년 골절 위험 변화), ⑤ **시나리오 비교** (6개 치료 전략 직접 비교), ⑥ **바이오마커** (CTX, P1NP, bsALP 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("op_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("op_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg op_qsp_model.dot -o op_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [op_references.md](op_references.md) 참조 (약 54편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
