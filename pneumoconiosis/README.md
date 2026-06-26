# 진폐증 (Pneumoconiosis, PNM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![PNM QSP Model](pnm_qsp_model.png)](pnm_qsp_model.svg)

## 개요 (Overview)
진폐증은 실리카(규소), 석탄 분진, 석면 등 무기 분진의 장기 흡입으로 유발되는 직업성 폐질환으로, 전 세계적으로 수백만 명의 광산·건설·제조업 종사자가 노출되어 있으며 진행성 폐섬유화가 특징이다. 핵심 발병기전은 폐포 대식세포의 분진 탐식 → NLRP3 인플라마솜 활성화 → IL-1β·TNF-α 방출 → TGF-β 매개 근섬유아세포 증식·콜라겐 침착으로 이어지는 불가역적 섬유화 과정이다. 현재 노출 차단이 가장 중요하며, 피르페니돈·닌테다닙 등 항섬유화제와 NAC 항산화 치료가 보조적으로 사용된다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 분진 침착·과부하 경로 | 폐포 내 분진 축적, 섬모 청소 과부하 | 폐포 대식세포 활성화 |
| NLRP3 인플라마솜 경로 | 식세포 불안정화 → 카스파제-1·IL-1β·IL-18 방출 | 급성·만성 폐 염증 |
| 산화 스트레스 경로 | ROS 생성, GSH 고갈, 미토콘드리아 손상 | 세포 사멸·섬유화 촉진 |
| 대식세포 극화 경로 | M1→M2 전환, 대식세포 피로사이토시스 | TGF-β·IL-10 분비 증가 |
| TGF-β 섬유화 경로 | Smad2/3 신호, 근섬유아세포 활성화 | 콜라겐 침착·폐 경직 |
| 폐 기능 저하 경로 | FVC·DLCO 감소, 폐혈관 저항 증가 | 호흡 부전·폐동맥 고혈압 |

## 주요 약물 표적 (Drug Targets)
- **피르페니돈**: TGF-β 신호 억제·항섬유화·항산화 — IPF 적응증에서 효과 확인, 진폐증 적용 연구 중
- **닌테다닙**: PDGFR·FGFR·VEGFR 억제 — 섬유화 진행 억제
- **N-아세틸시스테인 (NAC)**: 항산화(GSH 전구체) — ROS 매개 손상 경감
- **테트란드린**: 칼슘 채널 차단·항염 — 진폐증 대상 임상 연구 진행 중
- **실리카 노출 차단**: 분진 방호구, 환기 시스템 — 근본적 예방

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pnm_qsp_model.dot](pnm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 144 노드 / 11 클러스터) |
| [pnm_qsp_model.svg](pnm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pnm_qsp_model.png](pnm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pnm_mrgsolve_model.R](pnm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 5 치료 시나리오) |
| [pnm_shiny_app.R](pnm_shiny_app.R) | Shiny 대시보드 |
| [pnm_references.md](pnm_references.md) | 참고문헌 (약 51편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 폐포 분진 침착·청소 구획, 대식세포(휴지·활성), NLRP3·IL-1β·TNF-α·TGF-β·IL-10·ROS·GSH·호중구 염증 구획, 근섬유아세포·콜라겐·FVC·폐혈관저항 섬유화 PD 구획, NAC/피르페니돈/닌테다닙/테트란드린 PK 구획
- **주요 치료 시나리오**: ① 노출 지속(무치료), ② NAC 600 mg TID, ③ 피르페니돈 801 mg TID, ④ 닌테다닙 150 mg BID, ⑤ NAC+피르페니돈 병용
- **보정/근거**: CAPACITY, ASCEND(IPF 피르페니돈) 및 INPULSIS(닌테다닙) 임상시험의 FVC 감소 억제 데이터 참고 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(분진 종류, 노출 기간, 폐기능 등급) · 분진 청소·염증 PK/PD · 섬유화 진행(FVC·DLCO) PD · 임상 엔드포인트(폐기능 감소율, 폐혈관 저항) · 치료 시나리오 비교(섬유화 억제 효과) · 산화 스트레스·사이토카인 바이오마커 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pnm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pnm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pnm_qsp_model.dot -o pnm_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pnm_references.md](pnm_references.md) 참조 (약 51편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
