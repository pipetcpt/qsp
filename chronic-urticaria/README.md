# 만성 자발성 두드러기 (Chronic Spontaneous Urticaria, CSU) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 알레르기·면역

[![CSU QSP Model](csu_qsp_model.png)](csu_qsp_model.svg)

## 개요 (Overview)

만성 자발성 두드러기(CSU)는 전 세계 유병률 약 1~2%의 만성 피부 알레르기 질환으로, 명확한 외부 유발 요인 없이 6주 이상 팽진(wheal)·홍반(flare)·소양감이 지속됩니다. 핵심 발병기전은 피부 비만세포(mast cell)의 FcεRI 결합 IgE 또는 자가항체(anti-FcεRI, anti-IgE)에 의한 자발적 활성화로, 히스타민·사이토카인(IL-31, IL-33, TSLP)을 방출하여 두드러기 삼극(Wheal·Flare·Itch)을 유도합니다. 혈청 총 IgE 상승 및 thyroid peroxidase(TPO) 자가항체 양성이 특정 표현형을 정의합니다. H1-항히스타민제(1차), 오말리주맙(2차), 두필루맙(3차/최신 승인)이 표준 치료 사다리를 구성합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| IgE/FcεRI 교차결합 | 자가항원-IgE → FcεRI 다량체 → 비만세포 탈과립 | 히스타민·프로스타글란딘·류코트리엔 방출 → 팽진·홍반 |
| 자가항체 경로 | Anti-FcεRIα/Anti-IgE IgG → 보체 활성화 → Mas-related GPR X2 | 비면역학적 MC 활성화 |
| Type-2 사이토카인 | IL-4·IL-13·IL-31·IL-33·TSLP → Th2/ILC2 증폭 | 만성 염증, IgE 클래스 전환 |
| JAK-STAT 경로 | IL-4/IL-13 → JAK1/TYK2 → STAT6 → FcεRI 발현 증가 | 비만세포 감수성 상승 |
| 신경면역 교차반응 | IL-31 → OSMRβ/IL-31RA → C 신경섬유 | 소양감(itch) |
| BTK 신호전달 | IgE-FcεRI → BTK → PLCγ → IP3/DAG → Ca²⁺ 유입 | 탈과립 핵심 단계 — BTKi 표적 |
| RANKL/IL-17(2차) | 호산구·Th17 → 연조직 리모델링 | 만성 피부 장벽 손상 |

## 주요 약물 표적 (Drug Targets)

- **H1-항히스타민제 (세티리진, 빌라스틴, 펙소페나딘)**: H1R 역효현제; 비만세포 탈과립 차단 → 소양감·팽진 억제
- **오말리주맙 (항-IgE)**: 유리 IgE와 결합 → FcεRI 발현 및 무장 비만세포 감소; ASTERIA·GLACIAL 시험 근거
- **리겔리주맙 (항-IgE, 고친화도)**: 오말리주맙 대비 ~50배 높은 IgE 친화도; 임상 3상 진행
- **두필루맙 (항-IL-4Rα)**: IL-4/IL-13 공통 수용체 차단; LIBERTY-CSU CUPID A·B 시험 (2023 NEJM) 승인 근거
- **리겔리주맙·베르세필리맙 (2세대 항-IgE)**: 고친화도 IgE 차단으로 더 빠른 IgE 억제
- **BTK 억제제 (리미브루티닙, 에보브루티닙)**: BTK 표적 → FcεRI 하위 신호 차단; 경구 투여 편의성
- **사이클로스포린**: 비만세포/T세포 calcineurin 억제; 항히스타민 불응 중증 CSU

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [csu_qsp_model.dot](csu_qsp_model.dot) | Graphviz 기계론적 지도 소스 (110+ 노드 / 9 클러스터) |
| [csu_qsp_model.svg](csu_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [csu_qsp_model.png](csu_qsp_model.png) | PNG 이미지 (150 dpi) |
| [csu_mrgsolve_model.R](csu_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획 / 7 치료 시나리오) |
| [csu_shiny_app.R](csu_shiny_app.R) | Shiny 대시보드 (8탭) |
| [csu_references.md](csu_references.md) | 참고문헌 (45편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조** (18개):
  - *PK 구획 (10개)*: 항히스타민제 GI·혈장, 오말리주맙 depot·중심·말초, 두필루맙 depot·중심·말초, BTKi GI·혈장
  - *PD 구획 (8개)*: 유리 IgE, IgE-오말리주맙 복합체, 비만세포 priming·활성화 지수, 피부·혈장 히스타민, 피부 IL-31, 피부 IL-33
- **주요 치료 시나리오 (7가지)**:
  1. 무치료 자연 경과
  2. 세티리진 10 mg QD
  3. 고용량 항히스타민제 40 mg/day
  4. 오말리주맙 300 mg q4wk
  5. 오말리주맙 300 mg + 항히스타민제 병용
  6. 두필루맙 300 mg q2wk (600 mg 부하 후)
  7. BTK 억제제 25 mg QD
- **임상 보정 데이터**: ASTERIA I·II (오말리주맙), GLACIAL (오말리주맙 고용량), LIBERTY-CSU CUPID A·B (두필루맙), Lowe 2014 (오말리주맙 PopPK)
- **UAS7 대리 지표**: 비만세포 활성화·피부 히스타민·IL-31 동태에서 산출

## Shiny 대시보드 (Dashboard)

8개 탭으로 구성:
1. **환자 프로파일** — 기저 IgE, 중증도, 자가항체 상태, 치료 선택
2. **약물 PK 프로파일** — 항히스타민제·오말리주맙·두필루맙·BTKi 혈중 농도
3. **IgE & 비만세포** — 유리 IgE 억제, 비만세포 priming·활성화, 히스타민 동태
4. **사이토카인 네트워크** — IL-31(소양감)·IL-33(알라민) 동태
5. **임상 엔드포인트** — UAS7 시간 경과, 두드러기 조절 율(WCU ≤ 6)
6. **시나리오 비교** — 7가지 치료 전략 UAS7 비교 및 요약 표
7. **바이오마커** — 혈청 IgE, IgE 억제율, 비만세포 활성화 대리지표, 호산구 지수
8. **References** — 모델 근거 및 주요 임상시험 참고문헌

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("csu_mrgsolve_model")
out <- mrgsim(mod, end = 24 * 168)  # 24주 시뮬레이션
plot(out)
# Shiny 대시보드:
shiny::runApp("csu_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
sfdp -x -Goverlap=scale -Tsvg csu_qsp_model.dot -o csu_qsp_model.svg
sfdp -x -Goverlap=scale -Tpng -Gdpi=150 csu_qsp_model.dot -o csu_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [csu_references.md](csu_references.md) 참조 (45편, EAACI 가이드라인·ASTERIA·GLACIAL·LIBERTY-CSU CUPID·BTKi 포함).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
