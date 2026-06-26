# Hemophilia A (혈우병 A) — QSP Model

**디렉토리:** `hemophilia-a/` | **약어:** HA | **모델 버전:** 1.0.0 | **작성일:** 2026-06-23

---

## 질환 개요 (Disease Overview)

혈우병 A(Hemophilia A)는 응고인자 VIII(FVIII) 결핍에 의한 X염색체 연관 열성 출혈 질환이다. 중증(FVIII <1 IU/dL), 중등도(1–5 IU/dL), 경증(>5 IU/dL)으로 분류되며, 전 세계 남성 약 10,000명 중 1명 빈도로 발생한다. 심각한 관절 내 출혈(혈관절증)이 반복되면 비가역적인 혈우병성 관절병증(Hemophilic Arthropathy)으로 진행하여 삶의 질이 극도로 저하된다.

핵심 병태생리:
- **내인성 경로 결함**: FVIII가 FIXa의 보조인자로 tenase 복합체를 구성; 결핍 시 FXa 생성 감소 → 트롬빈 생성 불충분
- **불안정 피브린 클롯**: 트롬빈 급증(burst) 부재로 섬유소 중합체 강도 저하
- **억제항체(Inhibitor)**: 중증 HA의 약 30%에서 항FVIII IgG 항체 발생 → 대체 치료 필요

---

## 핵심 치료 패러다임 (Treatment Paradigms)

| 전략 | 약물 | 기전 | 주요 임상시험 |
|------|------|------|-------------|
| **FVIII 대체** | 표준 반감기(SHL) FVIII | 직접 FVIII 보충 | Manco-Johnson 2007 NEJM |
| **연장 반감기** | Fc 융합·PEG화 FVIII | FcRn 재활용으로 t½ ↑ | A-LONG (Mahlangu 2014) |
| **비대체 예방** | **에미시주맙** (Hemlibra) | FIXa–FX 이중특이항체, FVIII 모방 | HAVEN 1/3/4 (Oldenburg 2017) |
| **항트롬빈 억제** | **피투시란** (fitusiran) | siRNA로 AT mRNA 녹다운 → 트롬빈↑ | ATLAS-INH (Young 2023) |
| **TFPI 억제** | **마스타시맙** (marstacimab) | 항TFPI 항체 → 외인성 경로 증폭 | BASIS 연구 |
| **유전자치료** | valoctocogene roxaparvovec | AAV5 벡터 FVIII 발현 | HOPE-B (Ozelo 2022) |

---

## 기계론적 지도 (Mechanistic Map)

[![Hemophilia A QSP Map](ha_qsp_model.png)](ha_qsp_model.svg)

> 클릭하면 고해상도 SVG 파일을 볼 수 있습니다.

**구성:** 167 노드 · 10 서브그래프 클러스터

| 클러스터 | 핵심 내용 | 치료 표적 |
|---------|---------|---------|
| **혈관 손상·혈소판** | TF 노출·vWF·GPIb·GPIIbIIIa·TXA2·ADP | 혈소판 활성화 초기 플러그 |
| **외인성 경로** | TF·FVIIa·TF/FVIIa 복합체·TFPI 억제 | 마스타시맙(TFPI 중화) |
| **FVIII 생물학·내인성 경로** | F8 유전자·VWF 보호·FVIIIa·Xase 복합체 | FVIII 대체요법·EHL |
| **공통 경로·트롬빈 생성** | Prothrombinase·Thrombin burst·FXIIIa·피브린 | TGA/ETP 임상지표 |
| **자연 항응고 기전** | AT·Protein C/S·EPCR·Thrombomodulin·tPA | 피투시란(AT 감소) |
| **억제항체 면역학** | TH세포·B세포·항FVIII IgG·BU 역가 | 면역관용 유도(ITI)·에미시주맙 |
| **약물 PK/PD** | SHL/EHL FVIII·에미시주맙·피투시란·데스모프레신 | 각 약물 노출-반응 |
| **출혈 표현형** | 혈관절증·근육내·두개내·GI 출혈·ABR | 예방요법 목표 ABR <3 |
| **혈우병성 관절병증** | 활막 철 침착·ROS·연골파괴·Pettersson 점수 | 관절 보호 예방 |
| **임상 엔드포인트** | FVIII trough ≥1%·ETP·ABR·EQ-5D·관절점수 | Zero-bleed 표현형 |

---

## mrgsolve ODE 모델 사양 (Model Specifications)

**파일:** [`ha_mrgsolve_model.R`](ha_mrgsolve_model.R)

### 구획 (Compartments) — 16개

| # | 구획 | 단위 | 설명 |
|---|------|------|------|
| 1 | `FVIII_C` | IU/dL | FVIII 중심 구획 |
| 2 | `FVIII_P` | IU/dL | FVIII 말초 구획 |
| 3 | `EMIC_SC` | mg | 에미시주맙 SC 데포 |
| 4 | `EMIC_C` | mg/L | 에미시주맙 중심 구획 |
| 5 | `EMIC_P` | mg/L | 에미시주맙 말초 구획 |
| 6 | `FITU_SC` | mg | 피투시란 SC 데포 |
| 7 | `FITU_C` | mg/L | 피투시란 중심 구획 |
| 8 | `AT_mRNA` | rel. | 항트롬빈 mRNA (기저치=1) |
| 9 | `AT_prot` | rel. | 항트롬빈 단백질 (기저치=1) |
| 10 | `Inhibitor` | BU/mL | FVIII 억제항체 역가 |
| 11 | `Thrombin_ETP` | norm. | 트롬빈 생성 포텐셜 |
| 12 | `CumBleeds` | count | 누적 출혈 건수 |
| 13 | `JointScore` | 0-100 | Pettersson 관절 점수 |
| 14 | `QoL` | 0-1 | 삶의 질 (EQ-5D) |
| 15 | `Synovitis` | 0-1 | 활막 염증 지수 |
| 16 | `FVIII_eff` | IU/dL | 유효 FVIII 활성 (FVIII + 에미시주맙) |

### 치료 시나리오 (Treatment Scenarios) — 7개

| 시나리오 | 요법 | 투여 방법 | 임상 보정 |
|---------|------|---------|---------|
| 1 | 예방요법 없음 | 온디맨드 | ABR ~30 (무치료 중증 HA) |
| 2 | SHL-FVIII 예방 | 25 IU/kg 3×/주 IV | Manco-Johnson 2007; ABR ~3–4 |
| 3 | EHL-FVIII 예방 | 50 IU/kg Q3-4일 IV | A-LONG 2014; ABR ~2–3 |
| 4 | 에미시주맙 Q1W | 1.5 mg/kg SC Q1W (3 mg/kg 부하 ×4) | HAVEN 3 2018; ABR 1.5 |
| 5 | 에미시주맙 Q4W | 6 mg/kg SC Q4W (부하 후) | HAVEN 4 2019; ABR 2.4 |
| 6 | 피투시란 Q1M | 80 mg SC Q1M | ATLAS-INH 2023; ABR ~0 |
| 7 | FVIII + 에미시주맙 | 병용 | 수술·고위험 시기 |

---

## Shiny 대시보드 (Interactive Dashboard)

**파일:** [`ha_shiny_app.R`](ha_shiny_app.R) | **탭: 6개**

| 탭 | 주요 기능 |
|----|---------|
| **1. Patient Profile** | 체중·FVIII 중증도·억제항체·치료 선택; Value Box 요약 |
| **2. FVIII PK** | FVIII 농도-시간 곡선 (선형/로그); 1%·15% trough 참조선; PK 요약 통계 |
| **3. PD Core Metrics** | ETP 트롬빈 생성 지표; 피투시란에 의한 AT 단백질 감소 추적 |
| **4. Bleed Risk & ABR** | 순간 ABR 시간경과; 관절 점수 + QoL 복합 플롯; Value Box |
| **5. Scenario Comparison** | 6개 시나리오 동시 비교; ABR + 관절 점수 장기 추이 |
| **6. Biomarkers** | FVIII 활성 vs. ETP 산포도; 억제항체 역가 동태; 임상 결과 요약표 |

```r
# 실행 방법
install.packages(c("shiny", "bslib", "plotly", "dplyr", "tidyr", "ggplot2"))
shiny::runApp("hemophilia-a/ha_shiny_app.R")
```

---

## 참고문헌 (References)

**파일:** [`ha_references.md`](ha_references.md)

55개 PubMed 인용 — 주요 섹션:
- WFH 가이드라인·역학 (5)
- FVIII 생물학·분자병리 (5)
- 트롬빈 생성 (5)
- FVIII 약동학 (5)
- SHL/EHL FVIII 임상시험 (4)
- 에미시주맙 HAVEN 1/3/4 (6)
- 피투시란 ATLAS-INH (4)
- 마스타시맙·기타 (3)
- 억제항체 발생 (4)
- 혈우병성 관절병증 (4)
- QoL·환자보고결과 (3)
- 유전자치료 (2)
- QSP/PK-PD 모델링 (5)

---

## 주요 임상 파라미터 요약

| 파라미터 | 중증 HA (미치료) | SHL-FVIII 예방 | 에미시주맙 Q1W |
|---------|----------------|--------------|-------------|
| ABR | ~30/년 | ~3–4/년 | ~1.5/년 |
| FVIII trough | <1 IU/dL | 1–5 IU/dL | ~15 IU/dL 등가 |
| ETP (정상 대비) | ~15% | ~40–60% | ~70–80% |
| QoL (EQ-5D) | 0.55–0.65 | 0.75–0.85 | 0.85–0.92 |
| 10년 관절 점수 | >50 | 15–25 | 8–15 |

---

## 산출물 요약

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`ha_qsp_model.dot/.svg/.png`](ha_qsp_model.svg) | **167 노드, 10 클러스터** (혈관손상·외인성경로·FVIII/내인성경로·공통경로·항응고기전·억제항체면역·약물PK/PD·출혈표현형·관절병증·임상엔드포인트) |
| ⚙️ mrgsolve ODE | [`ha_mrgsolve_model.R`](ha_mrgsolve_model.R) | **16구획 ODE** (FVIII 2구획·에미시주맙 3구획·피투시란 2구획·AT mRNA/단백질·억제항체·ETP·CumBleeds·관절점수·QoL·Synovitis·FVIII_eff), **7 치료 시나리오** |
| 📊 Shiny 앱 | [`ha_shiny_app.R`](ha_shiny_app.R) | **6탭** (환자 프로파일·FVIII PK·PD 핵심지표·출혈위험·시나리오 비교·바이오마커), bslib darkly, plotly, 내장 ODE 시뮬레이터 |
| 📚 참고문헌 | [`ha_references.md`](ha_references.md) | **55개 PubMed 인용** (HAVEN 1/3/4·ATLAS-INH·A-LONG·HOPE-B·Manco-Johnson 2007 등) |
