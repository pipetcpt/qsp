# 이상지질혈증 (Dyslipidemia, DYSLIP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![DYSLIP QSP Model](dyslip_qsp_model.png)](dyslip_qsp_model.svg)

## 개요 (Overview)

이상지질혈증(Dyslipidemia)은 LDL 콜레스테롤 상승, HDL 감소, 중성지방 상승 등 지질 항상성의 이상으로 정의되며, 죽상동맥경화증과 심혈관질환의 핵심 위험인자입니다. 전 세계 성인의 약 40%가 이상지질혈증을 보유하며, 고LDL혈증은 심근경색 위험을 2~3배 높입니다. 핵심 병태생리는 간의 HMG-CoA 환원효소 활성, LDL 수용체(LDLR) 발현 조절, PCSK9에 의한 LDLR 분해, 역콜레스테롤 수송(RCT) 경로의 항상성 이상입니다. 스타틴·에제티미브·PCSK9 억제제·인클리시란이 주요 치료제입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| HMG-CoA 환원효소 경로 | 간 콜레스테롤 합성 → VLDL 분비 | LDL-C 상승 |
| LDL 수용체 경로 | LDLR → LDL 간 흡수 → PCSK9에 의한 분해 | 혈중 LDL-C 농도 결정 |
| PCSK9 경로 | PCSK9 → LDLR 리소좀 분해 촉진 | LDL-C 상승, 심혈관 위험 증가 |
| 역콜레스테롤 수송 | ApoA1/HDL → 말초 콜레스테롤 → 간 배출 | HDL 보호 기전 |
| 장관 콜레스테롤 흡수 | NPC1L1 → 식이·담즙 콜레스테롤 흡수 | LDL-C 기여 |
| 죽상경화 진행 | 산화LDL → 대식세포 거품세포 → 플라크 | 심근경색·뇌졸중 |
| 중성지방 경로 | VLDL 과분비·LPL 감소 → TG 상승 | 췌장염, 잔류 위험 |

## 주요 약물 표적 (Drug Targets)

- **스타틴** (아토르바스타틴, 로수바스타틴): HMG-CoA 환원효소 억제 → 간 콜레스테롤 합성 감소, LDLR 발현 증가
- **에제티미브**: NPC1L1 억제 → 장관 콜레스테롤 흡수 억제
- **PCSK9 억제제** (에볼로쿠맙, 알리로쿠맙): PCSK9 결합 차단 → LDLR 보호 → LDL-C 50~60% 추가 감소
- **인클리시란**: PCSK9 mRNA siRNA → 간 PCSK9 합성 억제 (2회/년 피하주사)
- **피브레이트**: PPARα 활성화 → TG 감소, HDL 증가

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [dyslip_qsp_model.dot](dyslip_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 214 노드 / 13 클러스터) |
| [dyslip_qsp_model.svg](dyslip_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [dyslip_qsp_model.png](dyslip_qsp_model.png) | PNG 이미지 (150 dpi) |
| [dyslip_mrgsolve_model.R](dyslip_mrgsolve_model.R) | mrgsolve ODE 모델 (약 27 구획 / 6 치료 시나리오) |
| [dyslip_shiny_app.R](dyslip_shiny_app.R) | Shiny 대시보드 |
| [dyslip_references.md](dyslip_references.md) | 참고문헌 (약 62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(아토르바스타틴·에볼로쿠맙·에제티미브·인클리시란 각 depot/중심/말초) + PD 구획(HMG-CoA, 메발론산, 간 콜레스테롤, 지단백 동역학, LDL-C, HDL-C, TG, PCSK9, LDLR, 죽상경화 지표)
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 스타틴 단독(아토르바스타틴 40 mg), ③ 스타틴 + 에제티미브, ④ 스타틴 + PCSK9 억제제(에볼로쿠맙 140 mg Q2W), ⑤ 3제 병합(스타틴 + 에제티미브 + 에볼로쿠맙), ⑥ 인클리시란 + 스타틴(ORION 시험 기반)
- **보정/근거**: FOURIER(에볼로쿠맙), ODYSSEY OUTCOMES(알리로쿠맙), ORION-10/11(인클리시란), 4S/WOSCOPS 스타틴 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 LDL-C·HDL-C·TG·심혈관 위험도 설정), ② PK 탭(4개 약물 혈중 농도 및 PCSK9 수준), ③ 지질 PD 탭(LDL-C·HDL-C·TG·non-HDL-C 추이), ④ 임상 엔드포인트(10년 ASCVD 위험 변화·LDL 목표 달성률), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(PCSK9·LDLR·죽상경화 플라크 부피 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("dyslip_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("dyslip_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg dyslip_qsp_model.dot -o dyslip_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [dyslip_references.md](dyslip_references.md) 참조 (약 62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
