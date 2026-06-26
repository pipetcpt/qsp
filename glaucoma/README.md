# 녹내장 (Glaucoma, GLAUC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 안과

[![Glaucoma QSP Model](glauc_qsp_model.png)](glauc_qsp_model.svg)

## 개요 (Overview)

녹내장은 전 세계 비가역적 실명의 주요 원인으로, 2040년까지 약 1억 1,180만 명이 이환될 것으로 예측됩니다. 핵심 병태생리는 안압(IOP) 상승 또는 혈류 조절 이상에 의한 시신경유두(ONH) 구조 손상 → 망막신경절세포(RGC) 진행성 소실 → 시야 결손입니다. Goldman 방정식(`IOP = (F − Fu) / C + Pv`)이 모델의 수력학적 핵심으로, 섬유주(trabecular)와 포도막-공막(uveoscleral) 유출로, 그리고 섬모체에서의 방수 생성이 모두 약물 표적입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 방수 생성 | 섬모체 상피 CA-II·β2-AR·cAMP-PKA | F 증가 → IOP 상승 |
| 섬유주 유출 | ROCK1/2 → 액틴 세포골격·ECM → 유출 저항 | C 감소 → IOP 상승 |
| 포도막-공막 유출 | FP 수용체 → ECM 리모델링(MMP) | Fu 감소 → IOP 상승 |
| ONH 구조 손상 | 사상판(LC) 변형·지지 조직 압박 | 축삭 교통 차단 |
| RGC 세포사멸 | 미토콘드리아 스트레스·산화 스트레스·Caspase-3 | RGC 비가역적 소실 |
| RNFL 변성 | RGC 축삭 소실 → OCT RNFL 두께 감소 | 시야 손상 전구 지표 |
| 시야 손상 | RNFL 임계 소실(~20%) 이후 MD 악화 | 시각 기능 저하 |
| 안관류압(OPP) | OPP = MAP × 2/3 − IOP | 정상안압 녹내장 관련 |

## 주요 약물 표적 (Drug Targets)

- **라타노프로스트 (FP 효현제 PG 유사체)**: FP 수용체 → PKC·PI3K → MMP 분비 → ECM 리모델링 → 포도막-공막 유출 ↑ (~25–35%); IOP −25~30%; QD 점안
- **티몰롤 (β-차단제)**: β2-AR 길항 → cAMP↓ → PKA↓ → Na/K-ATPase 활성↓ → 방수 생성 ↓ (~20–25%); IOP −20~25%; BID 점안
- **도르졸라미드 (CA-II 억제제)**: CA-II 억제 → HCO₃⁻↓ → Na+ 공동수송↓ → 방수 생성 ↓ (~15–20%); IOP −15~20%; TID 점안
- **브리모니딘 (α2 효현제)**: α2-AR → cAMP↓ → 방수 생성 ↓ + 포도막-공막 유출 ↑; IOP −20%; BID 점안 + 신경보호 효과 (BDNF↑)
- **네타수딜 (ROCK 억제제)**: ROCK1/2 억제 → 액틴 세포골격 이완 → 섬유주 유출 ↑ (~50%); IOP −20%; QD 점안; ROCKET-1/2 (Serle 2018)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [glauc_qsp_model.dot](glauc_qsp_model.dot) | Graphviz 기계론적 지도 (147+ 노드, 12 클러스터) |
| [glauc_qsp_model.svg](glauc_qsp_model.svg) | SVG 벡터 이미지 |
| [glauc_qsp_model.png](glauc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [glauc_mrgsolve_model.R](glauc_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획, 7 시나리오) |
| [glauc_shiny_app.R](glauc_shiny_app.R) | Shiny 대시보드 (8 탭) |
| [glauc_references.md](glauc_references.md) | 참고문헌 (62편+, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조 (26 ODEs)**:
  - PK (10 구획): LAT_DEPOT·TIM_DEPOT·DZL_DEPOT·BRI_DEPOT·NET_DEPOT (점안 흡수), LAT_AC·TIM_AC·DZL_AC·BRI_AC·NET_AC (전방 챔버)
  - 신호전달 (5 구획): FP_RO (FP 수용체 점유), BETA_RO (β-AR 점유), ALPHA_RO (α2-AR 점유), CA_INH (CA-II 억제), ROCK_INH (ROCK 억제)
  - 방수역학 (5 구획): ECM_CM (ECM 리모델링), cAMP_CIL (섬모체 cAMP), AQ_PROD (방수 생성율), C_TRAB (섬유주 유출 용량), FU_UVEA (포도막-공막 유출)
  - 질환 (6 구획): IOP_cmt (안압 mmHg), ONH_STRESS (ONH 응력), RGC_PCT (RGC 생존율 %), RNFL_UM (RNFL 두께 µm), VF_MD (시야 MD dB), OPP_cmt (안관류압)
- **Goldman 방정식**: `IOP = (AQ_PROD − FU_UVEA) / C_TRAB + EVP`
- **RGC 소실 ODE**: `dRGC/dt = −kRGC × (IOP_excess² + 0.3 × OPP_deficit) × RGC_PCT`
- **RNFL 추적**: RGC 소실에 90일 지연(tau) 적용
- **주요 치료 시나리오 (7가지)**:
  1. 무치료 자연 경과
  2. 라타노프로스트 단독 QD
  3. 티몰롤 단독 BID
  4. 라타노프로스트 + 티몰롤 (Xalacom® 병합)
  5. 3제 병용 (라타노 + 티몰롤 + 도르졸라미드)
  6. 네타수딜 단독 QD (ROCKET-1/2)
  7. 네타수딜 + 라타노프로스트 (Rocklatan®)
- **보정 근거**: OHTS (Gordon 2002), EMGT (Heijl 2002), LiGHT (Gazzard 2019), CIGTS (Musch 1999), AGIS (AGIS investigators 2000), ROCKET-1/2 (Serle 2018)

## Shiny 대시보드 (Dashboard)

8개 탭:
1. **환자 프로파일** — 기저 IOP·위험인자·녹내장 유형 설정; OHTS 5년 위험도 계산기
2. **Drug PK** — 전방 챔버 농도·수용체 점유율·PK 파라미터 표
3. **방수역학 & IOP** — AQ 생성/유출, IOP 경시 변화, Goldman 산포도
4. **RGC & 신경보호** — RGC 생존율·ONH 응력·OPP·5년 RGC 막대 차트
5. **임상 엔드포인트** — RNFL 두께·시야 MD·VFI 추이·엔드포인트 요약표
6. **시나리오 비교** — 7가지 치료 전략 동시 비교 (IOP·RNFL·MD·5년 요약표)
7. **바이오마커 패널** — cAMP·AQ 생성·ECM 상태·레이더차트·IOP 피크/트로프
8. **민감도 분석** — kRGC·dFu_ECM·dC_ROCK·EVP·MAP 슬라이더로 토네이도 차트

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("glauc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365 * 24 * 5, delta = 24)
plot(out)
# Shiny 대시보드:
shiny::runApp("glauc_shiny_app.R")
```

## 참고 임상시험 (Key Clinical Trials)

| 시험 | n | 중재 | 결과 |
|------|---|------|------|
| OHTS (Gordon 2002) | 1,636 | 안압 하강 치료 vs 관찰 | 5년 녹내장 발생 4.4% vs 9.5% |
| EMGT (Heijl 2002) | 255 | 레이저+약물 vs 무치료 | 시야 악화 45% vs 62% (HR 0.50) |
| LiGHT (Gazzard 2019) | 718 | SLT 우선 vs 약물 우선 | 3년 IOP 목표달성 74% vs 71% |
| CIGTS (Musch 1999) | 607 | 수술 vs 약물 | 5년 MD −2.24 vs −1.83 dB |
| AGIS (2000) | 591 | 수술 순서 비교 | IOP <18 mmHg 유지 시 시야 안정 |
| ROCKET-1/2 (Serle 2018) | 1,613 | 네타수딜 0.02% vs 티몰롤 | IOP −3.3 vs −3.0 mmHg |

## 주요 파라미터 (Key Parameters)

| 파라미터 | 값 | 단위 | 출처 |
|----------|-----|------|------|
| 정상 방수 생성율 (F₀) | 2.5 | µL/min | Brubaker 1982 |
| 섬유주 유출 용량 (C₀) | 0.25 | µL/min/mmHg | Grant 1958 |
| 포도막-공막 유출 (Fu₀) | 0.5 | µL/min | Bill 1975 |
| 공막 정맥압 (EVP) | 9.0 | mmHg | Brubaker 1982 |
| 기저 IOP | 16.0 | mmHg | OHTS baseline |
| RGC 소실 속도 (kRGC) | 0.0005 | /IOP²·day | Quigley 2006 모델 보정 |
| RNFL 시간 상수 (τ) | 90 | 일 | Miki 2014 OCT |
| VF-RNFL 임계 | 20 | % RNFL 소실 | Hood & Kardon 2007 |
