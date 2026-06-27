# Waldenström's Macroglobulinemia (WM) — QSP Model

[![WM](wm_qsp_model.png)](wm_qsp_model.svg)

**분류:** 혈액종양학 (Hematologic Oncology) · **디렉토리:** [`waldenstrom-macroglobulinemia/`](.)

---

## 병태생리 (Pathophysiology)

발덴스트롬 거대글로불린혈증(WM)은 골수(BM)와 림프절에서 **림프형질세포성 림프종(LPL)**으로 분류되는 희귀 B세포 종양입니다. 진단 기준은 ①골수에서 림프형질세포 ≥10% 침윤, ②단클론 IgM 파라단백 검출이며, 전체 비호지킨 림프종의 약 1-2%를 차지합니다.

### 핵심 발병기전

| 경로 | 세부 메커니즘 |
|------|------------|
| **MYD88 L265P 돌연변이** (~95%) | MyD88 자발 이량체화 → IRAK4/IRAK1 → TRAF6 → IKKβ → NF-κB 구성적 활성화; MYD88→HCK→BTK 직접 축 |
| **CXCR4 WHIM 돌연변이** (~35%) | 수용체 내재화 감소 → 지속적 CXCL12/SDF-1 신호 → PI3K/AKT → BM 체류 증가; 이브루티닙 저항성과 연관 |
| **BCR/BTK 신호전달** | 항원→BCR→LYN/SYK→BTK→PLCγ2→PKCβ/DAG/IP3→CBM 복합체→IKK→NF-κB |
| **PI3K/AKT/mTOR** | BCR 및 CXCR4 → PI3K-δ/γ → PIP3 → AKT(pT308) → mTORC1 → 단백 합성·증식 |
| **BM 미세환경** | 기질세포 CXCL12/BAFF/APRIL → 종양세포 생존; 비만세포 IL-6 파라크린; Treg NK 억제 |
| **IgM 과다 생성** | LPC → 형질세포 분화(BLIMP1/IRF4) → 오량체 IgM 분비 → 과점도 증후군·한랭글로불린·항MAG 신경병증 |

---

## 모델 사양 (Model Specifications)

| 구성 요소 | 내용 |
|----------|------|
| **기계론적 지도** | 12 서브그래프 클러스터, 120+ 노드 (Drug PK·BTK/BCR·MYD88/NF-κB·PI3K/AKT/mTOR·CXCR4/BM·LPC 분화·TME·세포사멸·IgM 합병증·약물 PD·임상 결과·유전체) |
| **mrgsolve ODE** | 20 구획 (Ibrutinib/Zanubrutinib/Rituximab/Venetoclax PK + BTK 점유율·NF-κB·LPC·PC·IgM·Hgb·점도·BM 침윤·BCL-2·CD20·Apoptosis·Proteasome·NK) |
| **치료 시나리오** | 7개 (Watch & Wait · 이브루티닙 · 이브루티닙+리툭시맙 · 자누브루티닙 · R-벤다무스틴 · BDR · 베네토클락스) |
| **Shiny 앱** | 7탭 (환자 프로파일/IPSSWM · Drug PK · PD 핵심 · 임상 엔드포인트 · 시나리오 비교 · 바이오마커 · About) |
| **참고문헌** | 62개 (진단·역학·분자생물학·임상시험·PK/PD·합병증·신약) |

---

## 임상시험 보정 데이터 (Clinical Trial Calibration)

| 시험 | 요법 | ORR | VGPR/CR | PFS |
|------|------|-----|---------|-----|
| iNNOVATOR (Treon 2015) | Ibrutinib 420 mg/d | 91.5% | 30.4% | 69% @ 2년 |
| INNOVATE (Dimopoulos 2018) | Ibrutinib + Rituximab | 92% | 43% | 82% @ 30개월 |
| ASPEN (Tam 2020) | Zanubrutinib 160 mg BID | 93.7% | 28.4% | 84% @ 18개월 |
| Rummel 2013 | R-Bendamustine | 96% | 44% | mPFS 69 mo |
| Dimopoulos 2013 | BDR | 83% | 22% | mPFS 43 mo |
| Castillo 2018 | Venetoclax | 84% | 36% | 미도달 |

---

## 주요 약물 기전 (Drug Mechanisms)

| 약물 | 표적 | 기전 |
|------|------|------|
| **Ibrutinib** | BTK C481 (covalent) | BCR/MYD88→BTK 경로 차단 → NF-κB ↓, 증식 억제 |
| **Zanubrutinib** | BTK C481 (covalent) | 이브루티닙보다 선택적; EGFR/ITK 부작용 최소화 |
| **Rituximab** | CD20 | ADCC·CDC·직접 세포사멸; IgM flare 주의 |
| **Bortezomib** | 26S Proteasome β5 | IκBα 안정화 → NF-κB ↓; UPR → ER stress 세포사멸 |
| **Venetoclax** | BCL-2 BH3 (Ki <1 nM) | BAX/BAK 방출 → MOMP → Caspase cascade |
| **Bendamustine** | DNA 이중가닥 | 알킬화 + 퓨린 아날로그 효과 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| [`wm_qsp_model.dot`](wm_qsp_model.dot) | Graphviz 기계론적 지도 소스 |
| [`wm_qsp_model.svg`](wm_qsp_model.svg) | 벡터 지도 (인터랙티브 뷰) |
| [`wm_qsp_model.png`](wm_qsp_model.png) | 래스터 지도 (150 dpi) |
| [`wm_mrgsolve_model.R`](wm_mrgsolve_model.R) | mrgsolve ODE 모델 + 7개 시나리오 시뮬레이션 |
| [`wm_shiny_app.R`](wm_shiny_app.R) | 7탭 인터랙티브 Shiny 대시보드 |
| [`wm_references.md`](wm_references.md) | 62개 참고문헌 (섹션별 분류) |

---

## 실행 방법 (Usage)

```r
# mrgsolve 모델
library(mrgsolve)
source("wm_mrgsolve_model.R")   # 모델 컴파일 + 7개 시나리오 자동 실행

# Shiny 대시보드
shiny::runApp("waldenstrom-macroglobulinemia/wm_shiny_app.R")
```

```bash
# Graphviz 렌더링
dot -Tsvg wm_qsp_model.dot -o wm_qsp_model.svg
dot -Tpng -Gdpi=150 wm_qsp_model.dot -o wm_qsp_model.png
```

---

*생성일: 2026-06-27 · Claude Code Routine (CCR)*
