# 건선성 관절염 (Psoriatic Arthritis, PsA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![PsA QSP Model](psa_qsp_model.png)](psa_qsp_model.svg)

## 개요 (Overview)

건선성 관절염(PsA)은 건선 피부 병변을 동반하거나 선행하는 만성 염증성 관절염으로, 건선 환자의 약 25~30%에서 발생하며 전 세계 유병률은 0.1~0.25%입니다. TNF-α·IL-23/IL-17 축이 부착부염(enthesitis), 활막 관절염, 손발가락 염(dactylitis), 척추 염증(척추염) 및 손발톱 건선을 동시에 매개하는 것이 핵심 발병기전입니다. IL-23p19·IL-17A·TNF·JAK 경로가 주요 치료 표적이며, 생물학제제(항TNF·IL-17 억제제·IL-12/23·IL-23p19 억제제)와 JAK 억제제(토파시티닙·우파다시티닙·필고티닙)가 표준 치료를 구성합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| IL-23/Th17 축 | DC/대식세포 IL-23 → Th17 분화 → IL-17A·IL-17F·IL-22 분비 | 피부 건선 플라크, 관절 염증 |
| TNF-α 경로 | 시노비오사이트·대식세포 TNF → NF-κB → 염증성 사이토카인 cascade | 활막 비후, 골미란, CRP 상승 |
| 부착부염 경로 | 역학적 스트레스 + IL-17/TNF → 부착부 섬유아세포 활성화 | 부착부 통증, 골형성·골흡수 불균형 |
| JAK-STAT 신호 | IL-6/IL-2/IFN-γ → JAK1/3-STAT3/STAT1 활성화 | 활막 증식, 면역세포 분화 |
| RANKL/OPG 경로 | TNF·IL-17 → 파골세포 활성화 → 골미란 | 방사선 골 손상 진행 |
| IL-12/Th1 축 | IL-12 → Th1 → IFN-γ → 대식세포 M1 분극 | 만성 육아종 염증, 전신 증상 |
| 각질세포 과활성화 | IL-17·IL-22 → 각질세포 증식·CXCL 분비 | 피부 건선 플라크 형성 |

## 주요 약물 표적 (Drug Targets)

- **항TNF 제제 (아달리무맙, 에타너셉트, 세르톨리주맙)**: TNF-α 직접 중화; 관절·피부·방사선 진행 억제
- **IL-17A 억제제 (세쿠키누맙, 익세키주맙)**: IL-17A 차단; 피부·부착부염·축성 병변에 우수
- **IL-23p19 억제제 (구셀쿠맙, 리산키주맙, 틸드라키주맙)**: IL-23 선택 억제; 장기 관해 유도
- **IL-12/23p40 억제제 (우스테키누맙)**: IL-12·IL-23 이중 차단; 피부·관절 병변 동시 조절
- **JAK 억제제 (토파시티닙, 우파다시티닙)**: 인트라셀룰러 JAK-STAT 신호 억제; 광범위 사이토카인 차단
- **DMARD (메토트렉세이트, 설파살라진, 레플루노미드)**: 말초 관절 위주; 피부 병변에도 부분 효과

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [psa_qsp_model.dot](psa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 15 클러스터) |
| [psa_qsp_model.svg](psa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [psa_qsp_model.png](psa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [psa_mrgsolve_model.R](psa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 25 구획 / 약 21개 시나리오) |
| [psa_shiny_app.R](psa_shiny_app.R) | Shiny 대시보드 |
| [psa_references.md](psa_references.md) | 참고문헌 (약 41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 항TNF/항IL-17/항IL-23/JAK억제제의 2구획 PK + TNF·IL-17A·IL-23 생성-제거 PD 모듈, 활막 염증 지수(SII), PASI 스코어 구획, 관절 손상 진행(mTSS) 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, 메토트렉세이트, 아달리무맙, 세쿠키누맙, 구셀쿠맙, 우스테키누맙, 토파시티닙, 아달리무맙+메토트렉세이트 병용 등
- **보정/근거**: FUTURE 임상시험(세쿠키누맙), DISCOVER(구셀쿠맙), PSUMMIT(우스테키맙), OPAL BROADEN(토파시티닙), ADEPT(아달리무맙) 데이터 기반 파라미터 보정

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 관절 서브타입·피부 침범·동반 질환 설정; (2) **PK 프로파일** — 생물학제제/소분자 혈중 농도 경시 변화; (3) **PD 주요지표** — IL-17A·TNF·IL-23 억제율; (4) **임상 엔드포인트** — ACR20/50/70, PASI 75/90/100, mTSS 변화; (5) **시나리오 비교** — 단일 및 병용 치료 비교; (6) **바이오마커** — CRP, ESR, 부착부 초음파 스코어.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("psa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("psa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg psa_qsp_model.dot -o psa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [psa_references.md](psa_references.md) 참조 (약 41편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
