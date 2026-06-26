# 제2형 당뇨병 (Type 2 Diabetes, T2DM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![T2DM QSP Model](t2dm_qsp_model.png)](t2dm_qsp_model.svg)

## 개요 (Overview)

제2형 당뇨병(T2DM)은 말초 인슐린 저항성과 췌장 베타세포 기능 부전이 복합적으로 작용하여 발생하는 가장 흔한 대사 질환으로, 전 세계 약 5억 3천만 명이 이환되어 있다(IDF 2021). 초기에는 보상성 인슐린 과분비로 혈당이 유지되나, 베타세포 소진이 진행되면 결국 현성 고혈당으로 이어진다. 과잉 지방산(FFA)과 이소성 지방 축적이 간 및 말초 인슐린 저항성의 핵심 기전이며, GLP-1 분비 감소·DPP-4 과활성도 혈당 조절 장애를 악화시킨다. 메트포르민·GLP-1 수용체 작용제·SGLT-2 억제제·인슐린 등 다양한 약물 클래스가 상호보완적 기전으로 사용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 말초 인슐린 저항성 | IRS-1 세린 인산화, GLUT4 전위 감소, FFA·세라마이드 독성 | 골격근 포도당 흡수 저하 |
| 간 인슐린 저항성 | 간 지방 축적 → 인슐린 신호 억제 → EGP 증가 | 공복 고혈당, 간 글루코스 과잉 생산 |
| 베타세포 기능 부전 | 당·지질 독성, 산화 스트레스, 아밀로이드 침착 → 베타세포 소실 | 인슐린 분비 점진적 감소 |
| GLP-1·인크레틴 저하 | 장 L세포 GLP-1 분비 저하, DPP-4 분해 증가 | 식후 인슐린 반응 저하, 글루카곤 억제 실패 |
| 신장 포도당 재흡수 증가 | SGLT2 과발현, 요당 역치 상승 | 소변 포도당 손실 감소 → 혈당 상승 기여 |
| 지방·신장 손상 합병증 | 만성 고혈당 → AGE 축적, 사구체 여과율 저하, 알부민뇨 | eGFR 감소, UACR 증가, 심혈관 위험 |

## 주요 약물 표적 (Drug Targets)

- **메트포르민** — 간 EGP 억제(AMPK 활성화), 체중 중립: 1차 치료제
- **GLP-1 수용체 작용제** — 세마글루타이드(Semaglutide): 인슐린 분비 촉진·글루카곤 억제·체중 감소·심혈관 보호 (SUSTAIN, PIONEER)
- **SGLT-2 억제제** — 엠파글리플로진(Empagliflozin): 신장 포도당 재흡수 차단, 심부전·신장 보호 (EMPA-REG OUTCOME)
- **DPP-4 억제제** — 시타글립틴: GLP-1 분해 억제, 혈당 강하, 체중 중립
- **설포닐우레아** — 글리메피리드: 베타세포 인슐린 분비 촉진 (저혈당 위험)
- **인슐린** — 데글루덱: 기저 혈당 조절, 진행성 베타세포 부전 시 필수

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [t2dm_qsp_model.dot](t2dm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 331 노드 / 11 클러스터) |
| [t2dm_qsp_model.svg](t2dm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [t2dm_qsp_model.png](t2dm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [t2dm_mrgsolve_model.R](t2dm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 27 구획 / 7 치료 시나리오) |
| [t2dm_shiny_app.R](t2dm_shiny_app.R) | Shiny 대시보드 |
| [t2dm_references.md](t2dm_references.md) | 참고문헌 (약 40편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(메트포르민 3구획, 엠파글리플로진 1구획, 세마글루타이드 SC 2구획, DPP-4 억제제 1구획, 설포닐우레아 1구획, 인슐린 SC 2구획, 피오글리타존 1구획) + 포도당-인슐린 동역학(Gp·Gt·Ip·X-action) + 내분비(글루카곤·GLP-1·베타세포 질량) + 인슐린 저항성(IR_H·IR_P·FFA) + 체중 + 임상 엔드포인트(HbA1c·eGFR·UACR)로 총 약 27개 구획
- **주요 치료 시나리오**: (1) 무치료, (2) 메트포르민 단독, (3) 메트포르민 + 엠파글리플로진, (4) 메트포르민 + 세마글루타이드, (5) DPP-4 억제제, (6) 인슐린, (7) 3제 병용(Met + Empa + Sema)
- **보정/근거**: EMPA-REG OUTCOME(엠파글리플로진 심신 보호), SUSTAIN-6(세마글루타이드 MACE 감소), UKPDS(메트포르민 장기 이점), UKPDS 베타세포 기능 감소 모델 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·HbA1c 기저치·eGFR·인슐린 저항성 수준 설정) · PK 시각화(각 약물 혈중 농도-시간 곡선) · 혈당 PD 지표(공복 혈당·HbA1c 시계열) · 임상 엔드포인트(체중 변화·eGFR·UACR) · 치료 시나리오 비교(7개 요법 HbA1c 궤적) · 바이오마커 패널(GLP-1·글루카곤·베타세포 질량) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("t2dm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("t2dm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg t2dm_qsp_model.dot -o t2dm_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [t2dm_references.md](t2dm_references.md) 참조 (약 40편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
