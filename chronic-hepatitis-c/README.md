# 만성 C형 간염 (Chronic Hepatitis C, CHC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 감염성/간담도

[![HCV QSP Model](HCV_qsp_model.png)](HCV_qsp_model.svg)

## 개요 (Overview)

C형 간염 바이러스(HCV)는 전 세계 약 5,800만 명이 만성 감염 상태에 있는 혈액 매개 플라비바이러스로, 치료받지 않으면 20–30년에 걸쳐 간섬유화 → 간경변 → 간세포암(HCC)으로 진행한다. 직접작용 항바이러스제(DAA)의 등장으로 현재는 8–12주 투약으로 95% 이상의 지속바이러스반응(SVR12)을 달성할 수 있다. 본 QSP 모델은 Perelson/Neumann의 표적세포 제한 바이러스 동역학 틀에 현대 DAA PK/PD, 숙주 면역 모듈, 간섬유화 동역학을 통합하여 치료 전략별 바이러스학적 반응과 장기 간 예후를 시뮬레이션한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 바이러스 복제 | NS5B RNA-의존 RNA 중합효소 → 네거티브-strand 복제 | 고바이러스혈증(log10 6–7 IU/mL) |
| NS5A 복제 복합체 | NS5A 단백질 → 복제소체(membranous web) 형성 | RNA 복제 허브 |
| NS3/4A 프로테아제 | 폴리단백 처리 → 구조·비구조 단백 성숙 | 면역 회피(MAVS·TRIF 절단) |
| 표적세포 감염 | HCV E2–CD81·SR-BI–CLDN1–OCLN 수용체 복합체 | 간세포 감염 확산 |
| 선천면역 회피 | NS3/4A → MAVS/TRIF 절단 → IFN-β 억제 | 지속 감염 확립 |
| T세포 소진 | 만성 항원 노출 → PD-1/Tim-3 상향조절 | CD8+ CTL 기능부전 |
| 간섬유화 | 만성 염증 → 간성상세포(HSC) 활성화 → TGF-β/콜라겐 축적 | F0→F4 진행, 간경변 |
| HCC 발생 | 산화 스트레스·TP53/β-catenin 변이 | 연 1–5% (간경변에서) |

## 주요 약물 표적 (Drug Targets)

- **소포스부비르(SOF, NS5B 억제제)**: 뉴클레오시드 유사체 → 간세포 내 트리포스페이트(SOF-TP)로 활성화 → RNA 복제 종결
- **레디파스비르(LED) / 벨파타스비르(VEL) / 피브렌타스비르(PIB) (NS5A 억제제)**: 복제소체 조립 및 바이러스 분비 억제; 피코몰 수준 효능
- **글레카프레비르(GLE) (NS3/4A 프로테아제 억제제)**: 폴리단백 처리 차단 → 바이러스 성숙 억제
- **리바비린(RBV)**: 면역 조절·불충실 복제 유도, IFN 상승 효과; GT3/간경변 병용에 사용
- **페그인터페론-α(PEG-IFN)**: JAK-STAT → ISG 유도 → 범항바이러스; 현재는 DAA 불내성 특수 상황에만 사용

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [HCV_qsp_model.dot](HCV_qsp_model.dot) | Graphviz 기계론적 지도 소스 (130+ 노드 / 12 클러스터) |
| [HCV_qsp_model.svg](HCV_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [HCV_qsp_model.png](HCV_qsp_model.png) | PNG 이미지 (150 dpi) |
| [HCV_mrgsolve_model.R](HCV_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획 / 7 치료 시나리오) |
| [HCV_shiny_app.R](HCV_shiny_app.R) | Shiny 대시보드 (6 탭) |
| [HCV_references.md](HCV_references.md) | 참고문헌 (65편, 17 섹션) |

## 기계론적 지도 구성 (Mechanistic Map Clusters)

| 클러스터 | 내용 |
|----------|------|
| Virion Entry & Uncoating | CD81·SR-BI·CLDN1·OCLN 수용체, IRES 번역, 내포작용 |
| HCV Replication Complex | NS5B RNA 중합효소, 복제소체, NS5A 도메인 |
| Assembly & Secretion | 코어·E1/E2 조립, VLDL 경로, NS2/NS3 처리 |
| Innate Immune Evasion | MDA5/RIG-I-MAVS 경로, NS3/4A 절단, STAT1 억제 |
| Adaptive Immunity | CD4 Th1, CD8 CTL, B세포, PD-1/Tim-3 소진 |
| Hepatic Pathology | HSC 활성화, TGF-β, Metavir F-score, HCC 위험 |
| DAA PK Compartments | SOF 간내 트리포스페이트, NS5A 억제제 혈장, GLE/PIB PK |
| Drug PD Effects | εp (복제 억제), εi (감염성 억제), 병합 효능 |
| Perelson Viral Kinetics | T(표적세포), I(감염세포), V(바이러스 RNA) ODE |
| Clinical Outcomes | SVR12, 재발, 내성 관련 치환(RAS), 간외 합병증 |
| Viral Dynamics Parameters | β, δ, p, c 파라미터 구성 |
| DAA Regimens | SOF/LED, SOF/VEL, GLE/PIB, PEG-IFN/RBV 등 |

## mrgsolve ODE 모델 (ODE Model)

### 구획 구조 (20 Compartments)

| 범주 | 구획 | 설명 |
|------|------|------|
| DAA PK | `SOF_Tp` | 소포스부비르 간내 활성 트리포스페이트 |
| DAA PK | `LED_p` | 레디파스비르 혈장 농도 |
| DAA PK | `VEL_p` | 벨파타스비르 혈장 농도 |
| DAA PK | `NS5A_i` | NS5A 억제제 복합 효과 |
| DAA PK | `GLE_p` | 글레카프레비르 혈장 농도 |
| DAA PK | `PIB_p` | 피브렌타스비르 혈장 농도 |
| DAA PK | `RBV_p` | 리바비린 혈장 농도 |
| DAA PK | `RBV_RBC` | 리바비린 적혈구 내 축적 |
| IFN | `PEGIFN_p` | 페그인터페론-α 혈장 농도 |
| 바이러스 동역학 | `T_cell` | 표적 간세포 수 |
| 바이러스 동역학 | `I_cell` | 감염 간세포 수 |
| 바이러스 동역학 | `V_rna` | 혈장 HCV RNA (IU/mL) |
| 바이러스 동역학 | `V_def` | 결함 바이러스 입자 |
| 면역 | `CTL` | CD8+ CTL (소진 모델링) |
| 면역 | `NK_cell` | NK 세포 활성 |
| 면역 | `Treg_HCV` | HCV 특이 조절 T세포 |
| 간 병리 | `ALT` | ALT (간세포 손상 대리) |
| 간 병리 | `Fibro_met` | Metavir 섬유화 점수 (F0–F4) |
| 간 병리 | `HSC_act` | 간성상세포 활성화 지수 |
| 간 병리 | `HCC_idx` | HCC 누적 위험 지수 |

### 주요 치료 시나리오

| 시나리오 | 요법 | 기간 | 근거 임상시험 |
|----------|------|------|---------------|
| 1 | SOF/LED 병합 | 12주 | ION-1/2/3 |
| 2 | SOF/VEL 병합 | 12주 | ASTRAL-1/2/3 |
| 3 | GLE/PIB 병합 | 8주 | ENDURANCE, EXPEDITION |
| 4 | PEG-IFN + RBV | 48주 | ADVANCE, ILLUMINATE |
| 5 | SOF/VEL (간경변) | 24주 | ASTRAL-4 |
| 6 | SOF/LED + RBV (GT3/RAS) | 24주 | LONESTAR |
| 7 | 무치료 (자연 경과) | — | 자연사 코호트 |

### 주요 파라미터 (Perelson Viral Kinetics)

| 파라미터 | 기호 | 값 | 출처 |
|----------|------|-----|------|
| 바이러스 소실율 | c | 22 /일 | Neumann et al., Science 1998 |
| 감염 간세포 사망율 | δ | 0.08–0.15 /일 | Perelson et al., J Theor Biol 2013 |
| 바이러스 생성율 | p | 100 virions/cell/일 | Dahari et al., Hepatology 2007 |
| 감염 전달 계수 | β | 1.5×10⁻⁷ /virion/일 | Rong et al., Biophys J 2010 |
| SOF 생산 억제율 | εp_SOF | 0.999 | INSPIRE 모델 |
| NS5A 억제율 | εp_NS5A | 0.9999 | 피코몰 IC50 기반 |

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성:

1. **환자 프로파일**: 유전형(GT1–6)·간섬유화 등급(F0–F4)·기저 바이러스 부하·IL-28B 유전형 설정
2. **약동학(PK)**: DAA 혈장 농도 및 간내 SOF-TP 농도 시계열; εp·εi 동역학
3. **PD 주요지표**: 혈장 HCV RNA log10 IU/mL, ALT, 감염 간세포 수
4. **임상 엔드포인트**: 바이러스학적 반응 타임라인(RVR/EVR/SVR12), 섬유화 진행, HCC 위험
5. **시나리오 비교**: 7가지 치료 전략의 1년 결과 비교, DataTable 요약
6. **면역 지형**: CTL 동역학(소진 포함), NK/Treg 비율, CTL:Treg 비

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("HCV_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("HCV_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg HCV_qsp_model.dot -o HCV_qsp_model.svg
dot -Tpng -Gdpi=150 HCV_qsp_model.dot -o HCV_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [HCV_references.md](HCV_references.md) 참조 (65편, 17 섹션):

- 역학 · 바이러스학 · 바이러스 동역학 모델링
- 선천/적응 면역 · CTL 소진 · NK세포
- DAA 각 계열 (NS5B / NS5A / NS3 프로테아제 / RBV / PEG-IFN)
- 내성 관련 치환(RAS) · SVR 결과
- 섬유화 동역학 · 간외 합병증 · WHO 박멸 목표 · IL-28B 유전형

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
