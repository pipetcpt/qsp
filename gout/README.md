# 통풍 (Gout, GOUT) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![GOUT QSP Model](gout_qsp_model.png)](gout_qsp_model.svg)

## 개요 (Overview)

통풍은 혈청 요산(uric acid)의 과포화로 인해 관절 내 MSU(monosodium urate) 결정이 침착되고, NLRP3 인플라마솜-IL-1β 축을 통한 급성 염증성 관절염 발작이 반복되는 대사질환입니다. 전 세계 성인 유병률은 약 1~4%이며, 서구에서는 가장 흔한 염증성 관절염입니다. 요산 생성 과잉(퓨린 대사 과잉)과 신장 배설 감소(URAT1, ABCG2 변이) 모두 과요산혈증의 원인이 되며, 만성화 시 통풍결절(tophus) 형성과 신장 손상이 동반됩니다. 요산저하제(알로푸리놀·페북소스타트)와 급성 발작 억제제(콜히친·NSAIDs·IL-1 차단제)가 치료의 양축을 이룹니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 퓨린 대사·XO 경로 | 잔틴산화효소(XO) → 요산 생성 과잉; URAT1/ABCG2 수송체 변이 → 배설 감소 | 혈청 요산 상승(>6.8 mg/dL) |
| MSU 결정 형성 | 과포화 → 핵 생성; 온도·pH·결합단백 조건에 따른 결정화 | 관절 내 결정 침착 |
| NLRP3 인플라마솜 활성 | MSU → TLR2/4 신호 + NLRP3 조립 → Caspase-1 → IL-1β 성숙·분비 | 급성 통풍 발작 |
| 호중구 동원 | IL-1β·TNF-α → E-셀렉틴·ICAM-1 상향 → 호중구 침윤 | 관절 부종·통증·발열 |
| 통풍결절 형성 | 만성 MSU 축적 → 대식세포·섬유모세포 반응 → 결절 | 관절 손상, 요로결석 |
| 신기능 손상 | 고요산혈증 → 세뇨관 결정 침착 → 사이토카인 지속 | eGFR 감소 |

## 주요 약물 표적 (Drug Targets)

- **알로푸리놀 (Allopurinol)**: XO 경쟁적 억제 → 요산 생성 감소; 옥시퓨리놀로 전환되어 지속 억제
- **페북소스타트 (Febuxostat)**: XO 비경쟁적 억제제; 신기능 저하 환자에서 알로푸리놀 대안
- **레시누라드/벤즈브로마론 (Lesinurad/Benzbromarone)**: URAT1 억제 → 요산 신배설 촉진
- **콜히친 (Colchicine)**: 튜불린 중합 억제 → 호중구 이동 차단; 급성 발작 예방·치료
- **NSAIDs (인도메타신 등)**: COX-1/2 억제 → 프로스타글란딘 감소; 급성 발작 1차 치료
- **아나킨라/카나키누맙 (Anakinra/Canakinumab)**: IL-1 수용체 길항제 / 항-IL-1β 단클론항체; 난치성 발작

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gout_qsp_model.dot](gout_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 120+ 노드 / 13 클러스터) |
| [gout_qsp_model.svg](gout_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gout_qsp_model.png](gout_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gout_mrgsolve_model.R](gout_mrgsolve_model.R) | mrgsolve ODE 모델 (약 37 구획 / 8개 치료 시나리오) |
| [gout_shiny_app.R](gout_shiny_app.R) | Shiny 대시보드 |
| [gout_references.md](gout_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 요산 PK(장·중심·말초·활막 구획), MSU 결정·결절, 염증 매개체(IL-1β·TNF-α·호중구·통증·관절 손상·eGFR), 약물 PK(알로푸리놀·옥시퓨리놀, 페북소스타트, 레시누라드, 콜히친, 인도메타신, 아나킨라, 카나키누맙-IL-1β 결합체) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② 알로푸리놀 300 mg/day, ③ 페북소스타트 80 mg/day, ④ 알로푸리놀+레시누라드 병합, ⑤ 콜히친 급성 발작 프로토콜, ⑥ 인도메타신 급성 발작, ⑦ 카나키누맙 150 mg SC(난치성 발작), ⑧ 페북소스타트+콜히친 장기 예방
- **보정/근거**: CONFIRM/CONFIRMS (페북소스타트), ACR 2012 통풍 지침, Tardif et al. COLCOT 등 임상시험 데이터를 기반으로 요산 목표치(<6 mg/dL) 달성 곡선을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(혈청 요산 기저치, eGFR, 동반 질환) 설정 탭, 요산 PK 동역학, 결정·결절 형성 PD, 염증 지표(IL-1β·CRP·호중구), 급성 발작 및 임상 엔드포인트(발작 빈도, 관절 손상), 8개 치료 시나리오 비교, 바이오마커(혈청 요산, eGFR 추이) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gout_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gout_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gout_qsp_model.dot -o gout_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gout_references.md](gout_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
