# 그레이브스병 (Graves' Disease, GD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![GD QSP Model](gd_qsp_model.png)](gd_qsp_model.svg)

## 개요 (Overview)

그레이브스병은 TSH 수용체(TSHR)를 자극하는 자가항체(TRAb)에 의해 갑상선이 지속적으로 과자극되어 갑상선기능항진증이 발생하는 자가면역 내분비질환입니다. 갑상선기능항진증의 가장 흔한 원인으로, 유병률은 여성 10명 중 1명, 남성은 그보다 낮으며 30~40대에 호발합니다. B세포 유래 TRAb가 TSHR에 결합하여 cAMP 신호를 지속적으로 활성화시키면 T4·T3 과생성, 갑상선 종대, 빈맥 등의 전형적 증상이 나타납니다. 특징적인 안구 병변(그레이브스 안병증, GO)은 안와 섬유아세포의 IGF-1R·TSHR 공동 자극을 통해 발생합니다. 항갑상선제(메티마졸·PTU), 방사성요오드(RAI), 수술이 주요 치료 선택지입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TRAb 형성 | HLA-DR3/DQ2 제한 항원 제시 → B세포 자극 → 자극성 TRAb IgG 분비 | TSHR 지속 자극 |
| cAMP/PKA 신호 과활성 | TSHR-Gs-adenylyl cyclase → cAMP → TPO·NIS 상향 → T4·T3 과합성 | 갑상선기능항진 |
| HPT 축 억제 | fT3·fT4 과잉 → 시상하부·뇌하수체 음성 피드백 | TSH 억제(<0.01 mIU/L) |
| 갑상선 종대 | cAMP → VEGF·IGF-1 → 갑상선 세포 증식 | 갑상선종, 혈관 증가 |
| 심혈관 영향 | fT3 과잉 → 심박수 증가·심근 수축력 증가·QTc 단축 | 빈맥, 심방세동 위험 |
| 그레이브스 안병증 | TRAb + IGF-1R 활성 → 안와 섬유아세포·지방세포 GAG 축적 | 안구 돌출, 복시 |
| 골대사 영향 | fT3 과잉 → 파골세포 활성 증가 → 골흡수 | 골밀도 감소 |

## 주요 약물 표적 (Drug Targets)

- **메티마졸 (Methimazole, MMI)**: TPO(갑상선과산화효소) 억제 → T4·T3 합성 차단; 1차 항갑상선제
- **프로필티오우라실 (PTU)**: TPO 억제 + 말초 T4→T3 전환 억제(D1 억제); 임신 1삼분기·갑상선 폭풍에 선호
- **방사성요오드 (RAI, ¹³¹I)**: NIS를 통한 甲상선 내 집적 → β선에 의한 갑상선 세포 파괴
- **프로프라놀롤 (Propranolol)**: 비선택적 β차단제 → 빈맥·진전 즉시 조절 + 일부 T4→T3 전환 억제
- **블록-앤-리플레이스 (Block-and-Replace)**: 고용량 MMI + LT4 보충; 갑상선 기능 안정화
- **리툭시맙 (Rituximab)**: B세포 고갈; 중증 불응성 그레이브스 안병증에 적용

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gd_qsp_model.dot](gd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 9 클러스터) |
| [gd_qsp_model.svg](gd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gd_qsp_model.png](gd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gd_mrgsolve_model.R](gd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 6개 치료 시나리오) |
| [gd_shiny_app.R](gd_shiny_app.R) | Shiny 대시보드 |
| [gd_references.md](gd_references.md) | 참고문헌 (약 43편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(MMI, PTU, RAI 혈청·갑상선 내 방사선량, 프로프라놀롤) + 질환 PD(TPO 억제, 갑상선 질량, T4·T3·fT4·fT3·rT3·TSH 구획, B세포, TRAb, 골흡수, 심박수, 그레이브스 안병증 활성도) 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② 메티마졸 30 mg/day, ③ 방사성요오드 15 mCi, ④ MMI + 프로프라놀롤 병합, ⑤ 블록-앤-리플레이스(MMI + LT4), ⑥ 리툭시맙 + MMI(불응성 GO)
- **보정/근거**: Ross et al. (ATA 2016 가이드라인), Cooper et al. (NEJM 2005), Nakamura et al. (JCEM 2007) 데이터를 기반으로 TSH·fT4 정상화 시간 경과를 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 TRAb, fT4, 갑상선 크기, GO 중증도) 설정 탭, 약물 PK(MMI·RAI 동역학), 갑상선 호르몬 PD(TSH·fT4·fT3 궤적), 심혈관·골 임상 엔드포인트(심박수, 골흡수 지표), 그레이브스 안병증 활성도, 6개 치료 시나리오 비교, 바이오마커(TRAb·TSH·fT4) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gd_qsp_model.dot -o gd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gd_references.md](gd_references.md) 참조 (약 43편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
