# 가려운 결절 (Prurigo Nodularis, PN) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 피부과 / 신경면역

[![PN QSP Model](pn_qsp_model.png)](pn_qsp_model.svg)

## 개요 (Overview)

가려운 결절(Prurigo Nodularis, PN)은 만성 신경면역성 피부질환으로, 극심한 가려움증(평균 NRS 7–8/10)과 전신에 분포하는 딱딱한 과각화성 결절(수십~수백 개)이 특징이다. 미국 내 유병률은 약 720만 명으로 추정되며, 어두운 피부 톤을 가진 환자에서 현저히 높은 빈도를 보인다. 핵심 발병기전은 **IL-31/IL-31RA 축을 통한 신경 가려움 신호**, **Th2/ILC2 기반 IL-4·IL-13·IL-31 과잉 생성**, **피부 장벽 기능 저하**, **진피 신경 과증식(NGF 매개)**, **중추 감작(central sensitization)**, 그리고 **긁기-결절 형성 악순환**이다. 2022년 dupilumab(IL-4Rα 차단)이 최초로 FDA 승인을 받았고, 2023년 nemolizumab(IL-31Rα 차단)도 승인되어 새로운 표준 치료제로 자리 잡았다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Th2/ILC2 활성화 | IL-4, IL-13, IL-31, TSLP, IL-33 | 전신 Th2 편향, 피부 염증 |
| IL-31/IL-31RA 신호 | IL-31RA, OSMR, JAK1/2, STAT3/5 | 가려움 NRS↑, 수면 장애 |
| 피부 장벽 손상 | Filaggrin↓, Claudin-1↓, TEWL↑ | 알레르겐 침투, 염증 촉진 |
| IgE·비만세포 축 | FcεRI crosslink, 히스타민, PGD₂ | 가려움 급성 악화 |
| 진피 신경 과증식 | NGF (비만세포 유래), SP, CGRP | 신경 밀도↑, 이질통 |
| 중추 감작 | μ/κ 오피오이드 불균형, TRPV1, TRPA1 | 만성 가려움 지속 |
| 섬유화·결절 | TGF-β, IL-13 → 섬유모세포 | 과각화성 결절 형성 |

## 주요 약물 표적 (Drug Targets)

- **Dupilumab (IL-4Rα 차단)**: IL-4·IL-13 동시 차단; PRIME 1/2 Phase 3 — Week 24 IGA 0/1 달성률 44.1% vs 16.0% (위약)
- **Nemolizumab (IL-31Rα 차단)**: IL-31 직접 차단; ARCADIA 1/2 Phase 3 — Week 16 Peak Pruritus NRS ≥4점 감소 56.3% vs 20.9%
- **Tralokinumab (IL-13 중화)**: TRALooPN Phase 2b — Week 16 IGA 0/1 달성률 37.1% vs 11.8%
- **Cyclosporine (칼시뉴린 억제)**: Th1/Th2 전반 억제; 소규모 RCT에서 가려움 60% 감소
- **Nalbuphine ER (κ-작용/μ-길항 오피오이드)**: PRIME PNT01-02 — 가려움 NRS 유의한 감소
- **Topical corticosteroids (TCS)**: 국소 항염증; 보조 치료

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [pn_qsp_model.dot](pn_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드 / 10 클러스터) |
| [pn_qsp_model.svg](pn_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pn_qsp_model.png](pn_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pn_mrgsolve_model.R](pn_mrgsolve_model.R) | mrgsolve ODE 모델 (27 구획 / 7개 치료 시나리오) |
| [pn_shiny_app.R](pn_shiny_app.R) | Shiny 대시보드 (7개 탭) |
| [pn_references.md](pn_references.md) | 참고문헌 (45편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조 (27개 ODE)**:  
  - PK: Dupilumab 2구획 SC + TMDD (IL-4Rα), Nemolizumab 2구획 SC + TMDD (IL-31Rα), Tralokinumab 2구획 SC, Cyclosporine 1구획 경구, Nalbuphine ER 1구획 경구, TCS 피부 데포 2구획
  - 면역 PD: Th2 세포, IL-4, IL-13, IL-31, IgE, 비만세포, 호산구
  - 피부/신경 PD: TEWL, 진피 신경 밀도, 결절 수, 중추 감작, DNRS
- **주요 치료 시나리오 (7개)**:
  1. 위약 (무치료)
  2. Dupilumab 300mg SC Q2W (PRIME trial)
  3. Nemolizumab 60mg SC Q4W (ARCADIA trial)
  4. Tralokinumab 300mg SC Q2W (TRALooPN trial)
  5. Cyclosporine 5mg/kg/day PO
  6. Nalbuphine ER 54mg PO BID
  7. Dupilumab + TCS 병용
- **보정 근거**: PRIME/PRIME2(dupilumab), ARCADIA 1/2(nemolizumab), TRALooPN(tralokinumab), Siepmann 2013(cyclosporine), Zeidler 2018(nalbuphine ER)

## Shiny 대시보드 (Dashboard)

7개 탭으로 구성:
1. **환자 프로파일**: 체중·나이·성별·기저 IL-31·IgE 설정, 치료 시나리오 선택
2. **약동학 (PK)**: 혈장 약물 농도-시간 곡선, 수용체 점유율 (IL-4Rα, IL-31Rα)
3. **면역 PD**: IL-4/IL-13/IL-31 변화, Th2 세포·IgE, 비만세포·호산구
4. **피부/신경**: TEWL, 진피 신경 밀도, 결절 수, 중추 감작 상태
5. **임상 엔드포인트**: Itch NRS 시간 경과, IGA 점수, DNRS
6. **시나리오 비교**: 전체 7개 시나리오 Itch NRS 비교, Waterfall, IGA 막대 그래프
7. **바이오마커**: 혈청 IL-31, IgE, 호산구, TEWL-NRS 상관 산점도

## 실행 방법 (Usage)

```r
library(mrgsolve)
# ODE 모델 실행:
source("pn_mrgsolve_model.R")

# Shiny 대시보드:
shiny::runApp("pn_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링 (sfdp: 사이클 그래프 지원)
sfdp -Tsvg -Goverlap=prism pn_qsp_model.dot -o pn_qsp_model.svg
sfdp -Tpng -Gdpi=150 -Goverlap=prism pn_qsp_model.dot -o pn_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [pn_references.md](pn_references.md) 참조 (45편, PubMed 링크 포함).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
