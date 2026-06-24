# 상염색체 우성 다낭신 (Autosomal Dominant PKD, ADPKD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![ADPKD QSP Model](adpkd_qsp_model.png)](adpkd_qsp_model.svg)

## 개요 (Overview)

상염색체 우성 다낭신(ADPKD)은 PKD1(폴리시스틴-1) 또는 PKD2(폴리시스틴-2) 유전자 변이에 의해 발생하는 가장 흔한 단일유전자 신장질환으로, 출생 1,000명당 약 1명의 빈도를 보인다. 폴리시스틴 기능 소실은 세포내 Ca²⁺ 저하와 cAMP 상승을 유발하며, 이는 낭종 상피세포의 증식과 체액 분비를 촉진하여 낭종이 수십~수천 개로 증식·확장된다. 총 신장 용적(TKV)은 질환 진행의 핵심 대리지표로, 연간 약 5–6% 증가하며 TKV > 750 mL(신장 길이 > 16.5 cm)가 급속 진행자 기준이다. 무치료 시 ESRD로의 진행 중앙값은 PKD1 변이에서 58세, PKD2 변이에서 79세이다. 톨밥탄(V2 수용체 길항제)이 cAMP 경로를 억제하여 TKV 증가와 eGFR 저하를 유의하게 늦추는 최초의 승인 약물이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| cAMP 과잉 | Ca²⁺↓ → AC 활성화 → cAMP↑ | 세포 증식, CFTR 매개 체액 분비 |
| mTOR 경로 | PC1 소실 → mTOR 과활성 → 세포 성장 | 낭종 상피 증식 |
| RAAS 활성화 | 낭종에 의한 허혈 → 레닌 증가 | 고혈압, 사구체 고혈압 |
| AVP/V2 수용체 | 집합관 V2R → cAMP → AQP2 발현 및 낭종 체액 분비 | TKV 증가 가속 |
| 섬유화 | TGF-β·PDGF → 낭종 주변 섬유화 | 기능 신원 소실 |
| 낭종내 압력 | 낭종 확장 → 인근 혈관·신원 압박 | eGFR 저하 |
| 안지오텐신 II | 산화스트레스, 세포자멸 촉진 | 신기능 악화 |

## 주요 약물 표적 (Drug Targets)

- **톨밥탄 (Tolvaptan, V2R 길항제)**: TEMPO 3:4 및 REPRISE 임상시험에서 TKV 증가율 약 49% 감소, eGFR 저하 지연 입증; 간독성 모니터링 필요
- **에베롤리무스 (Everolimus, mTOR 억제제)**: TKV 억제 효과 있으나 신기능 보호 미입증; 임상 사용 제한적
- **소마토스타틴 유사체 (Octreotide LAR)**: 간낭종 및 신낭종 성장 억제; 간다낭신에 더 효과적
- **ACE 억제제/ARB**: 혈압 조절 및 신보호, 현재 표준 치료
- **연구 중 표적**: MEK 억제제, CFTR 억제제, 항섬유화제

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [adpkd_qsp_model.dot](adpkd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 184 노드 / 10 클러스터) |
| [adpkd_qsp_model.svg](adpkd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [adpkd_qsp_model.png](adpkd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [adpkd_mrgsolve_model.R](adpkd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 다수 치료 시나리오) |
| [adpkd_shiny_app.R](adpkd_shiny_app.R) | Shiny 대시보드 |
| [adpkd_references.md](adpkd_references.md) | 참고문헌 (약 53편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 톨밥탄(경구 3구획 PK) + 에베롤리무스(2구획) + 옥트레오타이드 LAR 데포 + ACEi 구획 + AVP 수준·집합관 cAMP·mTOR 활성·안지오텐신 II·혈압 PD 구획 + TKV·eGFR·뇨삼투압·기능 신원 분율 결과 구획
- **주요 치료 시나리오**: ① 무치료, ② 톨밥탄 단독(저위험), ③ 톨밥탄 단독(고위험/급속 진행), ④ 에베롤리무스 + 톨밥탄 병용, ⑤ ACEi + 톨밥탄, ⑥ 옥트레오타이드 병용
- **보정/근거**: TEMPO 3:4 임상시험(TKV, eGFR 데이터), REPRISE 연구, Torres et al. NEJM 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(PKD1/PKD2·초기 TKV·eGFR·연령) 탭, 톨밥탄 PK 및 V2R 점유율 탭, TKV 성장 곡선 탭, eGFR 저하 예측 탭, 치료 시나리오 비교 탭, 뇨삼투압·바이오마커 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("adpkd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("adpkd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg adpkd_qsp_model.dot -o adpkd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [adpkd_references.md](adpkd_references.md) 참조 (약 53편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
