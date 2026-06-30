# 미세 대장염 (Microscopic Colitis, MC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 만성 위장관·면역

[![MC QSP Model](mc_qsp_model.png)](mc_qsp_model.svg)

## 개요 (Overview)
미세 대장염(Microscopic Colitis, MC)은 만성 또는 간헐적 수양성 비혈변 설사를 특징으로 하는 만성 염증성 장 질환으로, 대장 내시경 소견은 정상이거나 거의 정상이지만 조직검사에서 특징적 이상을 보입니다. 두 가지 주된 아형으로 **림프구성 대장염(LC; CD8+ 상피내림프구 ≥ 20/100 EC)**과 **콜라겐성 대장염(CC; 상피하 콜라겐대 ≥ 10 μm)**가 있습니다. 60대 이상 여성에서 흔하며, NSAID/PPI/SSRI 등 약물 노출이 강력한 유발 인자입니다. 셀리악병·자가면역 갑상선염 등과 자주 동반되고, 부데소니드 MMX가 표준 치료입니다.

## 핵심 병태생리 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 약물 트리거 | NSAID(COX·점막 손상), PPI(담즙·미생물군 변동), SSRI(점막 5-HT) | 발병 위험 ↑ 2–7배 |
| Th1/CD8+ IEL 면역 | IFN-γ, IL-15, T-bet 매개 상피내 림프구 증식 | 상피 세포 사멸 ↑ |
| 사이토카인 폭풍 | TNF-α, IL-6, IL-17A, IL-23 | 장벽 단백질(claudin/ZO-1) 감소 |
| 콜라겐대 침착(CC) | 페리크립트 근섬유모세포의 TGF-β/MMP-9 불균형, type I/III/VI 콜라겐 | 상피하 콜라겐대 비후 |
| 담즙산 흡수 장애 | 회장 ASBT 기능 저하, FGF19 ↓, CYP7A1 ↑ | 대장 BA 부하 ↑ → 분비성 설사 |
| 점막 수송 변화 | NHE3 ↓, NKCC1 ↑, CFTR/cAMP 활성화 | 흡수 → 분비 전환 |
| 미생물군 변화 | F. prausnitzii/Akkermansia ↓, Proteobacteria ↑ | SCFA ↓ · LPS ↑ |
| 점막 장벽 약화 | TJ 단백질·MUC2·Paneth defensin 변화 | 항원 통과 → 면역 활성화 |

## 주요 약물 표적 (Drug Targets)
- **부데소니드 MMX (Budesonide MMX)** — 국소 GR 활성, NF-κB 전사 억제로 유도/유지 표준
- **콜레스티라민 / 콜레세벨람** — 담즙산 격리, BAM 동반 환자에 효과
- **메살라민 (5-ASA)** — PPAR-γ 활성, 보조 치료
- **비스무트 서브살리실레이트** — 항분비·항감염
- **로페라마이드** — μ-opioid 매개 점막 분비/연동 감소 (대증)
- **아자티오프린 / 메토트렉세이트** — 부데소니드 의존/불응 환자
- **인플릭시맙 / 아달리무맙** — 항-TNF, 중증 불응성
- **베돌리주맙** — 항-α4β7 통합소, 장 선택적

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [mc_qsp_model.dot](mc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 130 노드 / 13 클러스터) |
| [mc_qsp_model.svg](mc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mc_qsp_model.png](mc_qsp_model.png) | PNG 썸네일 (150 dpi) |
| [mc_mrgsolve_model.R](mc_mrgsolve_model.R) | mrgsolve ODE 모델 (33 구획 / 10 치료 시나리오) |
| [mc_shiny_app.R](mc_shiny_app.R) | Shiny 대시보드 (7 탭) |
| [mc_references.md](mc_references.md) | 참고문헌 81편 (PubMed) |

## mrgsolve 모델 (ODE Model)
- **PK 구획**: 부데소니드 3구획, 콜레스티라민 1구획(루멘), 메살라민 2구획, 비스무트·로페라마이드, 아자티오프린→6-TGN, 인플릭시맙·아달리무맙·베돌리주맙 mAb PK
- **PD/질환 구획**: CD8+ IEL, IFN-γ, TNF-α, IL-6, IL-17, TGF-β, 상피 장벽, 콜라겐대 두께(μm), 담즙산 부하, 순 수분 흐름, 배변 횟수, Hjortswang 점수, HRQoL, BMD, HPA
- **치료 시나리오 (10개)**:
  1. 자연 경과
  2. 부데소니드 MMX 9 mg 8주 유도 + 6 mg 6개월 유지
  3. 콜레스티라민 4 g TID 6개월
  4. 메살라민 2.4 g/일 8주
  5. 비스무트 524 mg QID 8주
  6. 로페라마이드 4 mg BID
  7. 아자티오프린 2 mg/kg
  8. 인플릭시맙 5 mg/kg q8w
  9. 베돌리주맙 300 mg q8w
  10. 부데소니드 + 콜레스티라민 병용
- **보정 근거**: BUC-60/BUC-63(부데소니드), Ung 2000(콜레스티라민), VICTORIA-MC·Riviere 2019(베돌리주맙), Esteve 2011(항-TNF)

## Shiny 대시보드 (Dashboard)
환자 프로파일 → 약물 PK → 점막 면역(사이토카인/IEL) → 조직학(콜라겐대/장벽) → 임상(배변·Hjortswang·HRQoL·수분) → 시나리오 비교 → 장기(BMD·HPA) 의 7개 탭.

## 실행 방법 (Usage)
```r
library(mrgsolve)
source("mc_mrgsolve_model.R")
out <- mc_model %>% mrgsim(events = mc_scenarios[["02_budesonide_taper"]],
                            end = 365*24, delta = 24)
plot(out, STOOL_CLIP ~ time)
shiny::runApp("mc_shiny_app.R")
```
```bash
dot -Tsvg mc_qsp_model.dot -o mc_qsp_model.svg
dot -Tpng -Gdpi=150 mc_qsp_model.dot -o mc_qsp_model.png
```

## 임상적 의의 (Clinical Relevance)
- MC는 만성 설사 환자의 약 10–20%를 차지하며 65세 이상 여성에서 진단율이 가장 높습니다.
- 부데소니드는 단기 임상 관해율 80–90%, 그러나 중단 후 60–80% 재발하여 저용량 유지 요법의 정량적 설계가 임상적으로 중요합니다.
- 본 QSP 모델은 약물·담즙산·점막 장벽·면역·콜라겐대를 통합하여 ① 유도-유지 전략, ② 담즙산 격리제 병용, ③ 베돌리주맙·항-TNF 의 효율-위험 시뮬레이션을 지원합니다.

---
*Generated: 2026-06-30 · Claude Code Routine · QSP Disease Model Library*
