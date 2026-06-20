# 만성 갑상선 기능 저하증 (Chronic Hypothyroidism, HYPO) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![HYPO QSP Model](hypo_qsp_model.png)](hypo_qsp_model.svg)

## 개요 (Overview)

만성 갑상선 기능 저하증은 성인 여성의 약 5%, 남성의 약 1~2%에서 발생하는 흔한 내분비 질환으로, 하시모토 갑상선염이 가장 흔한 원인이다. 시상하부-뇌하수체-갑상선(HPT) 축의 피드백 조절이 와해되어 TSH 상승·FT4 저하가 특징적으로 나타나며, 대사율 저하, 심박수 감소, LDL 상승, 골 무기질 밀도 이상 등 전신 효과가 복합적으로 발생한다. 레보티록신(LT4) 단독 또는 LT4+리오티로닌(LT3) 병용 보충요법이 표준 치료이며, 용량 적정이 핵심 관건이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| HPT 축 피드백 파괴 | T4/T3 저하 → TRH↑ → TSH↑ | 갑상선 자극 과잉, 갑상선 종대 |
| 말초 T4→T3 전환 | 5'-탈요오드효소(DIO1/DIO2) 활성 저하 | 조직 내 활성형 T3 부족 |
| 대사 저하 | 기초대사율(BMR) 감소, 산소 소비량 저하 | 체중 증가, 피로, 저체온 |
| 심혈관 효과 | 심박수·심박출량 감소, 심근 수축력 저하 | 서맥, 심낭 삼출 |
| 지질 대사 이상 | LDL 수용체 감소 → LDL 상승 | 심혈관 위험 증가 |
| 골 대사 이상 | 골 재형성 속도 저하 → BMD 변화 | 장기 갑상선 기능 항진증 과다 치료시 BMD 감소 |
| 신경·근육 효과 | Na-K-ATPase 활성 저하, 신경전도 속도 저하 | 인지 장애, 근육 경련 |

## 주요 약물 표적 (Drug Targets)

- **레보티록신(LT4)**: T4 보충 → 말초 T4→T3 전환으로 효과 발휘, 긴 반감기(7일)로 1일 1회 복용
- **리오티로닌(LT3)**: 직접 활성형 T3 보충, 조직 반응 빠르지만 반감기 짧음(1일) → 병용 시 분할 투여
- **탈요오드효소 조절**: 셀레늄 등 보조영양소로 DIO1/DIO2 활성 지지 (보조요법)
- **TSH 수용체**: LT4 적정 용량 피드백의 핵심 표적 — 과잉 억제 시 골 소실·심방세동 위험

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [hypo_qsp_model.dot](hypo_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 421 노드 / 13 클러스터) |
| [hypo_qsp_model.svg](hypo_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [hypo_qsp_model.png](hypo_qsp_model.png) | PNG 이미지 (150 dpi) |
| [hypo_mrgsolve_model.R](hypo_mrgsolve_model.R) | mrgsolve ODE 모델 (약 16 구획 / 7개 치료 시나리오) |
| [hypo_shiny_app.R](hypo_shiny_app.R) | Shiny 대시보드 |
| [hypo_references.md](hypo_references.md) | 참고문헌 (약 38편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: TRH, TSH, 총T4(TT4), 총T3(TT3), rT3 내인성 호르몬 동태 + LT4/LT3 gut·중심(2구획) PK + 심박수(HR), LDL, 기초대사율(BMR), 증상 점수, 골밀도(BMD) PD 구획
- **주요 치료 시나리오**: (1) 정상 甲, (2) 무치료 갑상선 기능 저하, (3) LT4 100 μg QD, (4) 부분적 갑상선 잔여 + LT4 50 μg, (5) LT4 175 μg (과잉 치료), (6) LT4 100 μg + LT3 10 μg 병용, (7) LT4 125 μg 고용량
- **보정/근거**: Jonklaas J et al. (Thyroid 2014) 권고안 및 Hoermann et al. TSH-FT4 비선형 관계 데이터 기반

## Shiny 대시보드 (Dashboard)

7개 탭으로 구성: (1) 환자 프로파일 — 갑상선 잔여 기능·체중·연령 설정; (2) PK 탭 — LT4/LT3 혈청 농도·반감기 시각화; (3) HPT 축 — TRH·TSH·FT4·FT3 시계열; (4) PD 주요지표 — HR, BMR, LDL; (5) 임상 엔드포인트 — 증상 점수, BMD; (6) 시나리오 비교 — 7가지 용량 전략 비교; (7) 바이오마커 — TSH 목표 범위 달성, T3/T4 비율

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("hypo_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("hypo_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg hypo_qsp_model.dot -o hypo_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [hypo_references.md](hypo_references.md) 참조 (약 38편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
