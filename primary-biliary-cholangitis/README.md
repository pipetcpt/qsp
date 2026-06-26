# 원발성 담즙성 담관염 (PBC) (Primary Biliary Cholangitis, PBC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![PBC QSP Model](pbc_qsp_model.png)](pbc_qsp_model.svg)

## 개요 (Overview)
원발성 담즙성 담관염(PBC)은 간내 소담관(small intrahepatic bile ducts)에 대한 자가면역 공격으로 담즙정체·간섬유화가 진행되는 만성 담즙성 간질환이다. 중년 여성에서 압도적으로 호발하며(여:남 = 10:1), 유병률은 10만 명당 약 40명이다. 항미토콘드리아 항체(AMA, 특히 anti-PDC-E2)가 95% 이상에서 양성이며, 담관 상피세포(cholangiocytes)에 대한 CD4/CD8 T세포·NK세포 매개 손상이 핵심 기전이다. UDCA가 표준 1차 치료이며, UDCA 불충분 반응자에 FXR 작용제(오베티콜산) 및 PPAR 작용제(엘라피브라노르·셀라델파르·베자피브레이트)가 병용된다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가면역 담관 손상 경로 | AMA → PDC-E2 인식, CD4/CD8·NK세포 침윤 | 소담관 파괴·담즙정체 |
| 담즙산 독성 경로 | 소수성 담즙산 축적 → 세포 사멸·산화 스트레스 | 간세포·담관 상피 손상 |
| FXR-FGF19 경로 | FXR 활성화 → FGF19 방출 → 담즙산 합성 억제 | 담즙산 부하 감소 |
| PPAR 항염·항섬유화 경로 | PPARα/δ 활성 → 담즙산 독성 감소·염증 억제 | ALP·GGT 감소, 섬유화 억제 |
| 간 섬유화 경로 | 간성상세포(HSC) TGF-β 자극 → 콜라겐 침착 | 간경변·문맥 고혈압 |
| 담즙산 장간 순환 | 장내 담즙산 재흡수·엔테로헤파틱 순환 | 담즙산 풀(pool) 조절 |

## 주요 약물 표적 (Drug Targets)
- **UDCA (우르소데옥시콜산)**: 친수성 담즙산 대체·세포 보호·면역조절 — 1차 표준치료
- **오베티콜산 (OCA, FXR 작용제)**: FXR 활성화 → FGF19 증가 → 담즙산 합성 억제 (POISE 시험)
- **엘라피브라노르 (PPARα/δ)**: 담즙산 독성·염증 감소 (ELATIVE Phase 3)
- **셀라델파르 (PPARδ)**: 담즙산 재흡수 조절·항염 (RESPONSE Phase 3)
- **베자피브레이트 (PPARα)**: 담즙산 합성 억제·ALP 정상화 (BEZURSO 시험)

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pbc_qsp_model.dot](pbc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 181 노드 / 10 클러스터) |
| [pbc_qsp_model.svg](pbc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pbc_qsp_model.png](pbc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pbc_mrgsolve_model.R](pbc_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 7 치료 시나리오) |
| [pbc_shiny_app.R](pbc_shiny_app.R) | Shiny 대시보드 |
| [pbc_references.md](pbc_references.md) | 참고문헌 (약 57편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: UDCA·OCA·엘라피브라노르·셀라델파르·베자피브레이트 PK 구획(위장관·중심), AMA·Th1세포·담관 상피 손상·담즙산 독성·FGF19·ALP·빌리루빈·GGT·섬유화·IgM PD 구획
- **주요 치료 시나리오**: ① 무치료, ② UDCA 단독(표준), ③ UDCA+OCA(POISE 요법), ④ UDCA+엘라피브라노르(ELATIVE), ⑤ UDCA+셀라델파르(RESPONSE), ⑥ UDCA+베자피브레이트(BEZURSO), ⑦ 삼중요법(UDCA+ELF+SEL)
- **보정/근거**: POISE(OCA), ELATIVE(엘라피브라노르), RESPONSE(셀라델파르), BEZURSO(베자피브레이트) Phase 3 임상시험의 ALP 반응률 및 GLOBE 점수 개선을 참고하여 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(UDCA 반응성, PBC-40 증상 점수, 간섬유화 단계) · 약물 PK 프로파일 · 담즙산·FXR/PPAR PD · 간 생화학 임상 엔드포인트(ALP·빌리루빈·GGT) · 치료 시나리오 비교(ALP 정상화율·GLOBE 점수) · 자가면역·섬유화 바이오마커(AMA·IgM·Fibroscan) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pbc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pbc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pbc_qsp_model.dot -o pbc_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pbc_references.md](pbc_references.md) 참조 (약 57편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
