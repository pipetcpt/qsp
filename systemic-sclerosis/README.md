# 전신경화증 (Systemic Sclerosis, SSc) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![SSc QSP Model](ssc_qsp_model.png)](ssc_qsp_model.svg)

## 개요 (Overview)

전신경화증(공피증, Scleroderma)은 혈관병증(vasculopathy), 자가면역 활성화, 그리고 진행성 섬유화가 복합적으로 나타나는 희귀 결합조직 질환으로 전 세계 유병률은 인구 100만 명당 약 150~300명으로 추정된다. 피부·폐·심장·신장·소화관 등 다기관에 콜라겐이 과잉 침착되어 비가역적 기관 손상을 초래한다. 핵심 발병기전은 내피세포 손상 → TGF-β/IL-6 과분비 → 섬유아세포 활성화 → 근섬유아세포 전환 → ECM 축적의 연쇄로 요약된다. 현재 치료 전략은 면역억제(마이코페놀레이트·토실리주맙), 항섬유화(닌테다닙), 혈관확장(보센탄·일로프로스트) 세 축으로 구성된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 혈관 내피 손상 | 항내피세포 항체, Endothelin-1 상승, NO·PGI2 감소 | 레이노 현상, 디지털 궤양, PAH |
| TGF-β 과분비 | 섬유아세포 TGF-βRII 과발현, Smad2/3 인산화 | 콜라겐 I·III 합성 증가, ECM 축적 |
| IL-6/Th17 축 | IL-6 자가분비 루프, Th17:Treg 불균형 | 전신 염증, 자가항체 생산 |
| 근섬유아세포 활성화 | FAct → α-SMA+ Myo 전환, PDGFR·FGFR 신호 | 피부 및 폐 섬유화(mRSS 상승, FVC 저하) |
| 폐혈관 리모델링 | ET-1↑, PVR 증가, RV 비대 | 폐동맥 고혈압(PAH), 6MWT 단축 |
| B세포 자가면역 | 항Scl-70, 항ACA 항체, B세포 활성화 | 면역복합체 침착, 기관 손상 증폭 |

## 주요 약물 표적 (Drug Targets)

- **FGFR/PDGFR/VEGFR 억제제** — 닌테다닙(Nintedanib): 섬유아세포 증식·콜라겐 합성 억제 (SENSCIS trial)
- **IL-6R 차단제** — 토실리주맙(Tocilizumab): IL-6 신호 차단, FVC 보존 (faSScinate, focuSSced trial)
- **T세포/B세포 억제제** — 마이코페놀레이트 모페틸(MMF/MPA): 면역억제, 피부·폐 안정화
- **엔도텔린 수용체 길항제(ERA)** — 보센탄(Bosentan): ET-1 차단, 디지털 궤양 예방 (RAPIDS-1/2)
- **프로스타사이클린 유사체** — 이로프로스트(Iloprost, 흡입): PGI2 유사 작용, PVR 감소 (STEP trial)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ssc_qsp_model.dot](ssc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 367 노드 / 15 클러스터) |
| [ssc_qsp_model.svg](ssc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ssc_qsp_model.png](ssc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ssc_mrgsolve_model.R](ssc_mrgsolve_model.R) | mrgsolve ODE 모델 (약 30 구획 / 5 치료 시나리오) |
| [ssc_shiny_app.R](ssc_shiny_app.R) | Shiny 대시보드 |
| [ssc_references.md](ssc_references.md) | 참고문헌 (약 49편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(닌테다닙 2구획 경구, 토실리주맙 2구획 IV, MMF/MPA 1구획, 보센탄 1구획, 이로프로스트 흡입 1구획) + 질환 PD 구획(TGF-β, IL-6, Th17, Bnv, 섬유아세포활성화, 근섬유아세포, Col1/Col3, ECM, mRSS, FVC, DLCO, ET-1, NO, PGI2, 내피세포 완전성, PVR, 6MWD)으로 총 약 30개 구획
- **주요 치료 시나리오**: (1) 닌테다닙 150 mg BID, (2) 토실리주맙 8 mg/kg IV q4w, (3) MMF 3,000 mg/day, (4) 보센탄 125 mg BID, (5) 닌테다닙 + 토실리주맙 + MMF 병용
- **보정/근거**: SENSCIS 임상시험(닌테다닙 FVC 보존), faSScinate/focuSSced(토실리주맙 mRSS), RAPIDS-1/2(보센탄 디지털 궤양), STEP trial(이로프로스트) 데이터를 근거로 파라미터 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중, 질환 아형, 기저 mRSS·FVC·PVR 설정) · PK 시각화(각 약물 혈중 농도-시간 곡선) · 섬유화 PD 지표(mRSS·FVC·DLCO 시계열) · 혈관 PD 지표(ET-1·PVR·6MWD) · 치료 시나리오 비교(단일요법 vs 병용요법) · 바이오마커 패널(TGF-β·IL-6·Th17) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ssc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ssc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ssc_qsp_model.dot -o ssc_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ssc_references.md](ssc_references.md) 참조 (약 49편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
