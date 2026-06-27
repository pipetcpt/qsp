# 방광암 (Bladder Cancer / Urothelial Carcinoma, BLCA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 비뇨기종양

[![BLCA QSP Model](blca_qsp.png)](blca_qsp.svg)

## 개요 (Overview)

방광암(요로상피암, Urothelial Carcinoma)은 방광 요로상피세포에서 발생하는 악성 종양으로, 남성에서 여성보다 약 3-4배 더 흔하며 전 세계 매년 약 57만 명의 신환이 발생한다. 비침습성(NMIBC: Ta/T1/CIS)과 근층침습성(MIBC: T2-T4) 및 전이성(mUBC)으로 분류된다. 주요 분자 병태는 FGFR3 돌연변이(저등급 NMIBC ~60%), TP53·RB1·CDKN2A 변이(고등급 MIBC), PI3K/AKT/mTOR 경로 활성화이며, PD-L1 발현을 통한 면역 회피가 핵심 기전이다. 치료는 비침습성 단계의 BCG 방광내 주입 면역요법부터, 전이성 단계의 GC 화학요법(시스플라틴+젬시타빈), PD-1/PD-L1 면역관문억제제(펨브롤리주맙, 아테졸리주맙), FGFR 억제제(에르다피티닙), 항체-약물 복합체(엔포르투맙 베도틴, 사시투주맙 고비테칸)까지 다양하다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 의의 |
|------|----------------|-----------|
| FGFR3 → RAS/MAPK | FGFR3 돌연변이/융합 (~60% NMIBC) | 에르다피티닙 표적; ORR 40.4% (BLC2001) |
| PI3K/AKT/mTOR | PIK3CA, ERBB2 증폭 | 세포 증식·생존 |
| TP53/RB1/CDKN2A | 세포 주기 탈조절 | 고등급 MIBC, 불량 예후 |
| PD-L1/PD-1 면역 회피 | IFN-γ 유도 PD-L1 상향 | 펨브롤리주맙/아테졸리주맙 표적 |
| TGF-β/Treg/MDSC 면역억제 | 종양 미세환경 면역 회피 | BCG 효능 제한 |
| Nectin-4/TROP2 표면 항원 | 요로상피암에 높은 발현 | EV/SG ADC 치료 표적 |
| BCG → TLR2/4 → CD8+ CTL | TH1 반응, 직접 종양 사멸 | NMIBC 표준치료 |

## 주요 약물 표적 (Drug Targets)

- **BCG (방광내 주입)**: TLR2/4 자극 → DC 성숙 → Th1/CD8 면역활성; NMIBC CRR 55–70% (SWOG S8507)
- **시스플라틴+젬시타빈 (GC)**: DNA 가교 + 뉴클레오시드 유사체; mUBC ORR 49%, mOS 13.8개월 (von der Maase 2000)
- **펨브롤리주맙 (Pembrolizumab)**: 항PD-1; 2L mUBC ORR 21.1%, OS HR 0.73 (KEYNOTE-045)
- **아테졸리주맙 (Atezolizumab)**: 항PD-L1; 1L 시스플라틴 부적합 ORR 15% (IMvigor210/211)
- **에르다피티닙 (Erdafitinib)**: FGFR1-4 억제; FGFR 변이 mUBC ORR 40.4% (BLC2001)
- **엔포르투맙 베도틴 (Enfortumab Vedotin)**: 항Nectin-4 ADC→MMAE; 3L mUBC ORR 40.6%, OS HR 0.70 (EV-301)
- **사시투주맙 고비테칸 (Sacituzumab Govitecan)**: 항TROP2 ADC→SN-38; mUBC 구제치료 (TROPHY-U-01)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [blca_qsp.dot](blca_qsp.dot) | Graphviz 기계론적 지도 소스 (112 노드 / 11 클러스터) |
| [blca_qsp.svg](blca_qsp.svg) | SVG 벡터 이미지 (확대 가능) |
| [blca_qsp.png](blca_qsp.png) | PNG 이미지 (150 dpi) |
| [blca_mrgsolve_model.R](blca_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획 / 7개 치료 시나리오) |
| [blca_shiny_app.R](blca_shiny_app.R) | Shiny 대시보드 (7개 탭) |
| [blca_references.md](blca_references.md) | 참고문헌 (41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조 (20개)**:
  - PK: BCG_depot/BCG_eff (방광내), Cis_plasm (시스플라틴), Gem_plasm (젬시타빈), Pembro_c/p (펨브롤리주맙 2구획), Atezo_c/p (아테졸리주맙 2구획), Erda_dep/Erda_plasm (에르다피티닙 경구 2구획), EnFV_c (엔포르투맙 베도틴)
  - 면역: CD8_eff (세포독성 T세포), Treg, MDSC_cmt
  - 종양: TumBurd (종양 세포 부담), FGFR3act (FGFR3 신호), PDL1_lvl (PD-L1 발현), IFNg_cmt (IFN-γ)
  - 바이오마커: NMP22_cmt (소변 NMP22), SLD_cmt (종양 크기, mm)
- **주요 치료 시나리오 (7개)**:
  ① 무치료(Placebo) ② BCG 81mg 방광내 주입 주1×6→q3w 유지 ③ 시스플라틴 70mg/m²+젬시타빈 1000mg/m² IV q3w ④ 펨브롤리주맙 200mg IV q3w ⑤ 아테졸리주맙 1200mg IV q3w ⑥ 에르다피티닙 8mg 경구 QD (FGFR 변이 환자) ⑦ 엔포르투맙 베도틴 1.25mg/kg IV d1,8,15 q28d
- **보정 근거**: SWOG S8507(BCG), von der Maase 2000(GC), KEYNOTE-045(펨브롤리주맙), IMvigor210/211(아테졸리주맙), BLC2001(에르다피티닙), EV-301(엔포르투맙 베도틴)

## Shiny 대시보드 (Dashboard)

7개 탭: ① **환자 프로파일** (성별·나이·체중·병기·FGFR3/PD-L1 상태), ② **약동학 PK** (혈중 농도, 수용체 점유율), ③ **면역 바이오마커** (CD8/Treg/MDSC/IFN-γ/PD-L1), ④ **종양 역학** (종양 세포 수, SLD, FGFR3 신호), ⑤ **임상 엔드포인트** (워터폴 플롯, RECIST, 요약표), ⑥ **시나리오 비교** (7개 치료 직접 비교), ⑦ **바이오마커 분석** (NMP22, PD-L1 vs IFN-γ 위상공간)

## 실행 방법 (Usage)

```r
library(mrgsolve)
source("blca_mrgsolve_model.R")
# Shiny 대시보드:
shiny::runApp("blca_shiny_app.R")
```

```bash
# 기계론적 지도 렌더링
dot -Tsvg blca_qsp.dot -o blca_qsp.svg
dot -Tpng -Gdpi=150 blca_qsp.dot -o blca_qsp.png
```

## 보정 요약 (Calibration Summary)

| 치료 | 임상시험 | 주요 엔드포인트 | 보정 파라미터 |
|------|---------|----------------|--------------|
| BCG | SWOG S8507 | CRR(CIS) 55–70% | Emax_BCG_kill=0.65, Emax_BCG_CD8=3.0 |
| GC | von der Maase 2000 | ORR 49%, mOS 13.8mo | Emax_Cis=0.55, Emax_Gem=0.45 |
| Pembrolizumab | KEYNOTE-045 | ORR 21.1%, OS HR=0.73 | IC50_Pembro=0.40, Emax_IO_kill=0.30 |
| Atezolizumab | IMvigor210/211 | ORR 15%, OS 11.1mo | IC50_Atezo=0.60 |
| Erdafitinib | BLC2001 | ORR 40.4% (FGFR-altered) | EC50_Erda=0.12, Emax_Erda_kill=0.40 |
| Enfortumab Vedotin | EV-301 | ORR 40.6%, OS HR=0.70 | EC50_EV=0.08, Emax_EV_kill=0.42 |

## 참고문헌 (References)

자세한 인용은 [blca_references.md](blca_references.md) 참조 (41편).

---

*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
