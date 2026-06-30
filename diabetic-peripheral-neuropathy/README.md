# Diabetic Peripheral Neuropathy (DPN) — QSP Model

> 당뇨병성 말초신경병증의 4중 대사 손상(폴리올·AGE·PKC·헥소사민) → 미세혈관/내막 저관류 → 신경섬유 손상 → 통증·감각저하 →
> 족부궤양으로 이어지는 인과 사슬을 통합한 정량적 시스템 약리학 모델.

| 항목 | 값 |
|------|-----|
| 디렉토리 | `diabetic-peripheral-neuropathy/` |
| 약어 | `dpn` |
| 분류 | 만성질환 · 신경/내분비 |
| 핵심 ODE 구획 | 26 (PK 10 + 질환/통증 16) |
| 약물 시나리오 | 8 (무치료·프레가발린·둘록세틴·병용·α-리포산·에팔레스타트·캡사이신 8%·집중혈당+병용) |
| Shiny 탭 | 8 |
| 기계론적 지도 노드 | 130+ (10 클러스터) |
| 참고문헌 | 60편 (18개 섹션) |

## 1. 산출물 (Deliverables)
| 파일 | 설명 |
|------|------|
| [`dpn_qsp_model.dot`](dpn_qsp_model.dot) / [`.svg`](dpn_qsp_model.svg) / [`.png`](dpn_qsp_model.png) | Graphviz 기계론적 지도 |
| [`dpn_mrgsolve_model.R`](dpn_mrgsolve_model.R) | mrgsolve QSP ODE 모델 |
| [`dpn_shiny_app.R`](dpn_shiny_app.R) | Shiny 대시보드 (8 탭) |
| [`dpn_references.md`](dpn_references.md) | PubMed 인용 60편 |

## 2. 병태생리 요약 (Pathophysiology)
1. **고혈당 → 세포내 포도당 축적**: 인슐린 비의존성 신경/슈반/내피세포에서 GLUT1로 유입 → 4가지 손상 경로 동시 활성화.
2. **폴리올 경로**: aldose reductase가 glucose→sorbitol→fructose 환원, NADPH 소모로 GSH 결핍 및 삼투압 스트레스.
3. **AGE/RAGE**: 비효소적 단백질 글리케이션 → RAGE → NF-κB → TNF-α·IL-6·VCAM-1.
4. **PKC-β**: DAG↑ → eNOS↓·ET-1↑·VEGF 이상.
5. **헥소사민**: GFAT를 통해 UDP-GlcNAc 생성 → O-GlcNAc 변형 → GAPDH·PAI-1 조절.
6. **미토콘드리아 과산화물·PARP**: NAD⁺ 고갈로 4 경로가 추가 증폭(Brownlee unified hypothesis).
7. **미세혈관 손상**: vasa nervorum 협착·기저막 비후·내막저혈류 → 내막 저산소증.
8. **신경섬유 손상**: 원위부 die-back, IENFD 감소, NCV 저하, NGF↓, Schwann 손상.
9. **통증 기전**: Nav1.7/1.8 과흥분, TRPV1/TRPA1 감작, 척수후각 NMDA·미세아교세포 활성 → 중추 감작.
10. **임상 결과**: 통증 NRS·MNSI·TCNS·BPI 악화, 균형/수면 손상, 족부궤양→절단→사망률.

## 3. 약물 PK/PD (Drug PK/PD)
| 약물 | 표적 | 모델링 깊이 |
|------|------|------------|
| Pregabalin | Cav α2δ-1 → 글루타메이트 방출↓ | 1차흡수·신장 CL eGFR 보정, EC50_pain=4 mg/L, Emax=0.45 |
| Duloxetine | SNRI → 하행 5-HT/NE | CYP2D6 표현형 조정, EC50=0.08 mg/L, Emax=0.35 |
| α-Lipoic acid | ROS 청소·GSH 재생 | 경구 F=0.30, IV F=1.0, EC50=0.5 mg/L |
| Epalrestat | Aldose reductase | IC50=5 mg/L, 폴리올 억제 |
| Capsaicin 8% | TRPV1 desensitisation | 효과구획, τ=90일 |
| Lidocaine 5% | Nav1.7/1.8 차단 | 국소 효과구획 |
| Ruboxistaurin*, Aminoguanidine* | 실패 약물 (참고 라인) | 음성 통제 |

(*) Ph3 음성, 음성 대조 목적으로 포함.

## 4. 보정 앵커 (Calibration anchors)
- **DCCT/EDIC**: 강력한 혈당 조절 → 신경병증 5년 발병 60% 감소.
- **SENZA-PDN HF10 SCS**: 3개월 통증 50% 이상 감소 76% vs CMM 5%.
- **ALADIN/NATHAN-1**: ALA 600 mg IV 3주 → TSS Δ −2.7; 4년 NATHAN-1은 일부 서브그룹에서 NIS-LL 개선.
- **Pregabalin 300-600 mg/d**: NRS Δ −1.3 vs 위약.
- **Duloxetine 60-120 mg/d**: NRS Δ −1.4.
- **Capsaicin 8% (STEP)**: NRS Δ −1.0, 12주 지속.
- **Epalrestat (ADCT)**: 3년 NCV 안정화.

## 5. 실행 (Usage)
```bash
# 기계론적 지도 렌더링
dot -Tsvg dpn_qsp_model.dot -o dpn_qsp_model.svg
dot -Tpng -Gdpi=150 dpn_qsp_model.dot -o dpn_qsp_model.png
```
```r
# mrgsolve 모델 (R)
install.packages(c("mrgsolve","dplyr","ggplot2"))
library(mrgsolve)
mod <- mread("dpn_mrgsolve_model.R")
ev  <- ev(amt=150, ii=12, addl=730, cmt="GUT_PG")   # Pregabalin 300 mg BID
out <- mrgsim(mod, events=ev, end=365, delta=1)
plot(out, NRS+IENFD+NCV~time)

# Shiny 대시보드
install.packages(c("shiny","shinydashboard","DT","tidyr","scales"))
shiny::runApp("dpn_shiny_app.R")
```

## 6. 한계 (Limitations)
- 평균 환자 표현형 기준 — 인구 시뮬레이션(가상 환자군)은 추후 확장 필요.
- 자율신경병증(CAN), 통증 없는 대섬유 변형, 운동·치료적 신경자극은 정성적으로만 표현.
- 캡사이신·리도카인 패치는 단순 효과구획으로 추상화 — 약동학적 국소 노출은 표현하지 않음.
- 모든 파라미터는 문헌·임상시험 평균치 기반의 출발점이며, 실데이터 적합 시 재추정 필요.

## 7. 변경 이력
- **2026-06-30** v1.0 — 초기 모델 (Claude Code Routine).
