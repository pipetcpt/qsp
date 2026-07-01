# Hemophilia B (혈우병 B, Christmas Disease) — QSP Model

**디렉토리:** `hemophilia-b/` | **약어:** HB | **모델 버전:** 1.0.0 | **작성일:** 2026-07-01

---

## 질환 개요 (Disease Overview)

혈우병 B(Hemophilia B)는 응고인자 IX(FIX) 결핍에 의한 X염색체 연관 열성 출혈 질환으로, 혈우병 A보다 약 4-5배 드물게 발생한다(남성 신생아 약 25,000-30,000명 중 1명). 중증(FIX <1 IU/dL), 중등도(1-5 IU/dL), 경증(>5 IU/dL)으로 분류되며, 반복적인 관절 내 출혈(혈관절증)이 비가역적인 혈우병성 관절병증으로 진행할 수 있다.

핵심 병태생리:
- **F9 유전자 결함**: Xq27.1에 위치한 *F9* 유전자의 missense/nonsense/large deletion 변이로 간세포에서 FIX 합성·분비 장애 발생; Hemophilia B Leiden 변이는 사춘기 이후 자연 호전이 특징적
- **비타민 K 의존성 감마-카르복실화**: GGCX/VKORC1 경로로 FIX Gla 도메인이 카르복실화되어야 인지질 표면 결합 및 응고 활성이 가능; 와파린은 이 경로를 차단
- **내인성 텐아제 결핍**: FIXa가 FVIIIa 보조인자와 Ca2+/인지질 표면에서 텐아제 복합체를 구성하지 못해 FX 활성화 및 트롬빈 생성 저하
- **억제항체(Inhibitor)**: 중증 HB의 약 1-3%(HA보다 낮음)에서 항FIX IgG 항체 발생하나, null-mutation 유전형에서 위험이 높고 **신증후군·아나필락시스**라는 HB 특유의 면역관용유도(ITI) 합병증 위험이 있음
- **FIX-Padua (R338L) 변이**: 자연 발생 초고활성 변이(비교 활성 ~8배)로 최신 AAV 유전자치료 벡터의 트랜스진으로 활용됨

---

## 핵심 치료 패러다임 (Treatment Paradigms)

| 전략 | 약물 | 기전 | 주요 임상시험 |
|------|------|------|-------------|
| **표준 반감기(SHL) FIX 대체** | 노나콕 알파(BeneFIX)·Rixubis | 직접 FIX 보충, IV | 등록연구 population PK |
| **연장 반감기(EHL) — Fc융합** | 에프트레노네콕 알파(Alprolix) | FcRn 재활용 → t½ ~82h | B-LONG (Powell 2013 NEJM) |
| **연장 반감기(EHL) — 알부민융합** | 알부트레페노네콕 알파(Idelvion) | 알부민 융합 → t½ ~102-104h | PROLONG-9FP (Santagostino 2016) |
| **연장 반감기(EHL) — 당페길화** | 노나콕 베타 페골(Refixia) | GlycoPEGylation → t½ ~93h | pathfinder2 (Collins 2014) |
| **TFPI 억제(비대체)** | **콘시주맙**(Concizumab) | 항TFPI 항체 → 외인성 경로 증폭 | explorer7/8 (Shapiro/Chowdary) |
| **TFPI 억제(비대체)** | **마스타시맙**(Marstacimab) | 항TFPI 항체, 주 1회 SC | BASIS (Pipe 2023 NEJM) |
| **항트롬빈 억제(비대체)** | **피투시란**(Fitusiran) | siRNA로 AT mRNA 녹다운 → 트롬빈↑ | ATLAS-A/B, ATLAS-INH (Young 2023) |
| **유전자치료** | 에트라나코진 데자파르보벡(Hemgenix) | AAV5 벡터, FIX-Padua 발현 | HOPE-B (Pipe 2023 NEJM) |
| **유전자치료** | 피다나코진 엘라파르보벡(Beqvez) | AAV-Rh74var 벡터, FIX-Padua 발현 | BENEGENE-2 (Cuker 2024) |

---

## 기계론적 지도 (Mechanistic Map)

[![Hemophilia B QSP Map](hb_qsp_model.png)](hb_qsp_model.svg)

> 클릭하면 고해상도 SVG 파일을 볼 수 있습니다.

**구성:** 125 노드 · 13 서브그래프 클러스터

| 클러스터 | 핵심 내용 | 치료 표적 |
|---------|---------|---------|
| **F9 유전자·간 합성·비타민K 회로** | F9 변이 유형·GGCX·VKORC1·Gla 카르복실화·FIX-Padua | 유전자치료 트랜스진 설계 |
| **내인성 경로·텐아제 복합체** | FXII/FXI/FIXa·FVIIIa 보조인자·텐아제·FX 활성화 | FIX 대체요법 표적 |
| **공통 경로·트롬빈 생성·피브린** | Prothrombinase·트롬빈 burst·FXIIIa·TAFI·플라스민 | ETP 임상지표 |
| **출혈 표현형·혈우병성 관절병증** | 자연 혈관절증·활막 철 침착·MMP·target joint | 예방요법 목표 ABR <3 |
| **FIX 대체요법 PK(SHL/EHL)** | BeneFIX/Alprolix/Idelvion/Refixia·FcRn 재활용 | trough 유지 전략 |
| **비대체 재균형요법** | 콘시주맙·마스타시맙(TFPI)·피투시란(AT) | 인자-비의존적 지혈 |
| **AAV 유전자치료** | 간세포 형질도입·에피솜 발현·캡시드 면역반응 | 1회 투여 근치적 접근 |
| **억제항체 면역학** | 항FIX IgG·Bethesda 역가·ITI·신증후군·아나필락시스 | HB 특유 ITI 위험관리 |
| **약력학·검사 바이오마커** | one-stage/chromogenic assay·TGA·시약 편차 | EHL 검사 오차 보정 |
| **임상 엔드포인트·PRO** | ABR/AJBR·HJHS·Haem-A-QoL·수술 지혈 | Zero-bleed 표현형 |
| **생애주기·특수집단 보정인자** | 신생아 생리적 저FIX·소아 Vd·간/신기능 | 용량 개별화 |
| **장기 동반질환·과거력 위험** | 만성 통증·정형외과 수술·과거 수혈력·사망률 | 삶의 질·기대수명 |
| **보조적 지혈 조치** | 트라넥삼산·데스모프레신(무효)·물리치료 | 국소/보조 관리 |

---

## mrgsolve ODE 모델 사양 (Model Specifications)

**파일:** [`hb_mrgsolve_model.R`](hb_mrgsolve_model.R)

### 구획 (Compartments) — 24개

| # | 구획 | 단위 | 설명 |
|---|------|------|------|
| 1 | `FIX_C` | IU/dL | FIX 중심 구획 (SHL) |
| 2 | `FIX_P` | IU/dL | FIX 말초 구획 (SHL) |
| 3 | `FIXe_C` | IU/dL | FIX 중심 구획 (EHL: Fc/알부민/당페길화 선택) |
| 4 | `FIXe_P` | IU/dL | FIX 말초 구획 (EHL) |
| 5 | `AAV_Vector` | rel. | 순환 AAV 벡터 게놈 |
| 6 | `Transduced_Hep` | 0-1 | 형질도입된 간세포 분획 |
| 7 | `Capsid_Antigen` | rel. | 일시적 캡시드 항원 제시 풀(면역반응만 구동, ~5주 반감기) |
| 8 | `Transgene_Expr` | IU/dL | 내인성 FIX-Padua 트랜스진 발현 |
| 9 | `Capsid_Immune` | 0-1 | 캡시드 특이 면역 활성화 |
| 10 | `ALT_level` | fold | ALT 배수 상승(간독성 지표) |
| 11 | `CONC_SC` | mg | 콘시주맙 SC 데포 |
| 12 | `CONC_C` | ng/mL | 콘시주맙 중심 구획 |
| 13 | `MARS_SC` | mg | 마스타시맙 SC 데포 |
| 14 | `MARS_C` | ng/mL | 마스타시맙 중심 구획 |
| 15 | `FITU_SC` | mg | 피투시란 SC 데포 |
| 16 | `FITU_C` | mg/L | 피투시란 중심 구획 |
| 17 | `AT_mRNA` | rel. | 항트롬빈 mRNA (기저치=1) |
| 18 | `AT_prot` | rel. | 항트롬빈 단백질 (기저치=1) |
| 19 | `Inhibitor` | BU/mL | 항FIX 억제항체 역가 |
| 20 | `Thrombin_ETP` | norm. | 트롬빈 생성 포텐셜 |
| 21 | `CumBleeds` | count | 누적 출혈 건수 |
| 22 | `JointScore` | 0-100 | 혈우병성 관절병증 점수 |
| 23 | `Synovitis` | 0-1 | 활막 염증 지수 |
| 24 | `QoL` | 0-1 | 삶의 질 (Haem-A-QoL 유래 효용) |

### 치료 시나리오 (Treatment Scenarios) — 10개

| 시나리오 | 요법 | 투여 방법 | 임상 보정 |
|---------|------|---------|---------|
| 1 | 예방요법 없음 | 온디맨드 | ABR ~28 (무치료 중증 HB) |
| 2 | SHL-rFIX 예방 | 40 IU/kg 2×/주 IV | 등록연구 population PK |
| 3 | EHL-rFIX-Fc 예방 | 50 IU/kg Q7-10일 IV | B-LONG 2013; t½ ~82h |
| 4 | EHL-rFIX-알부민 예방 | 75 IU/kg Q14일 IV | PROLONG-9FP 2016; t½ ~102-104h |
| 5 | 당페길화-rFIX 예방 | 40 IU/kg Q7일 IV | pathfinder2 2014; t½ ~93h |
| 6 | 콘시주맙 매일 SC | 부하 210mg + 유지 15mg/일 SC | explorer7/8 |
| 7 | 마스타시맙 주간 SC | 부하 300mg + 유지 150mg/주 SC | BASIS 2023 |
| 8 | 피투시란 월간 SC | 50 mg/월 SC | ATLAS-A/B 2023 |
| 9 | AAV 유전자치료 1회 | 벡터 IV 단회 투여 + 스테로이드 반응성 조절 | HOPE-B 2023 |
| 10 | 억제항체 양성 환자 | ITI 프로토콜 + 우회인자 온디맨드 관리 | Null-mutation 고위험 시나리오 |

---

## Shiny 대시보드 (Interactive Dashboard)

**파일:** [`hb_shiny_app.R`](hb_shiny_app.R) | **탭: 8개**

| 탭 | 주요 기능 |
|----|---------|
| **1. Patient Profile** | 체중·FIX 중증도·유전형(null mutation)·억제항체·치료 선택; Value Box 요약 |
| **2. FIX PK** | SHL/EHL FIX 농도-시간 곡선 (선형/로그); 1%·15% trough 참조선 |
| **3. Non-Factor Rebalancing** | 콘시주맙/마스타시맙 농도; 피투시란에 의한 AT 단백질 감소 추적 |
| **4. Gene Therapy** | AAV 트랜스진 발현 경과; 캡시드 면역반응 및 ALT 간독성 동태 |
| **5. PD Core Metrics** | ETP 트롬빈 생성 지표; 총 유효 FIX-당량(대체+트랜스진+TFPI) |
| **6. Bleed Risk & Arthropathy** | 순간 ABR 시간경과; 관절 점수 + QoL 복합 플롯; Value Box |
| **7. Scenario Comparison** | 9개 치료 옵션 동시 비교; ABR + 관절 점수 장기 추이 |
| **8. Inhibitor & Biomarkers** | FIX 활성 vs. ETP 산포도; 억제항체 역가 동태; 임상 결과 요약표 |

```r
# 실행 방법
install.packages(c("shiny", "bslib", "plotly", "dplyr", "tidyr", "ggplot2"))
shiny::runApp("hemophilia-b/hb_shiny_app.R")
```

---

## 참고문헌 (References)

**파일:** [`hb_references.md`](hb_references.md)

71개의 검증된 PubMed 인용 — 주요 섹션:
- 질환 개요·역학 (WFH 가이드라인)
- F9 유전자·분자유전학·FIX 생물학 (FIX-Padua 포함)
- 응고 캐스케이드·트롬빈 생성
- FIX 약동학 — SHL 제제
- FIX 약동학 — EHL 제제 (B-LONG·PROLONG-9FP·pathfinder2)
- 비대체 재균형요법 (콘시주맙·마스타시맙·피투시란)
- AAV 유전자치료 (HOPE-B·BENEGENE-2)
- 억제항체·면역 합병증 (신증후군·아나필락시스)
- 혈우병성 관절병증
- 임상 엔드포인트·QoL·가이드라인

---

## 주요 임상 파라미터 요약

| 파라미터 | 중증 HB (미치료) | SHL-rFIX 예방 | EHL-rFIX-Fc 예방 | AAV 유전자치료 |
|---------|----------------|--------------|-----------------|--------------|
| ABR | ~28/년 | ~13-16/년 | ~7-11/년 | ~4-6/년 |
| FIX trough/발현 | <1 IU/dL | 1-5 IU/dL | 5-10 IU/dL | ~15-30 IU/dL(지속) |
| ETP (정상 대비) | ~15% | ~30-40% | ~40-55% | ~60-75% |
| QoL (효용) | 0.55-0.65 | 0.70-0.80 | 0.75-0.85 | 0.80-0.90 |

---

## 산출물 요약

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`hb_qsp_model.dot/.svg/.png`](hb_qsp_model.svg) | **125 노드, 13 클러스터** (F9유전자·내인성경로·공통경로·출혈표현형·FIX대체PK·비대체재균형·AAV유전자치료·억제항체면역·바이오마커·임상엔드포인트·생애주기·장기동반질환·보조지혈) |
| ⚙️ mrgsolve ODE | [`hb_mrgsolve_model.R`](hb_mrgsolve_model.R) | **24구획 ODE** (FIX SHL/EHL 4구획·AAV유전자치료 5구획·콘시주맙/마스타시맙 4구획·피투시란/AT 4구획·억제항체·ETP·CumBleeds·관절점수·QoL·Synovitis), **10 치료 시나리오** |
| 📊 Shiny 앱 | [`hb_shiny_app.R`](hb_shiny_app.R) | **8탭** (환자프로파일·FIX PK·비대체재균형·유전자치료·PD핵심지표·출혈위험/관절병증·시나리오비교·억제항체/바이오마커), bslib darkly, plotly, 내장 ODE 시뮬레이터 |
| 📚 참고문헌 | [`hb_references.md`](hb_references.md) | **71개 PubMed 인용** (HOPE-B·B-LONG·PROLONG-9FP·pathfinder2·BASIS·ATLAS-A/B 등) |
