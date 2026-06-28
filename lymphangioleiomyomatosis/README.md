# Lymphangioleiomyomatosis (LAM) — QSP Model

[![LAM QSP Model](lam_qsp_model.png)](lam_qsp_model.svg)

**분류**: 희귀 폐질환 / 종양 억제 유전자 질환 (Rare Lung Disease / Tumor Suppressor Gene Disorder)  
**약어**: LAM  
**디렉토리**: [`lymphangioleiomyomatosis/`](.)

---

## 질환 개요 (Disease Overview)

**림프관평활근종증(LAM)**은 주로 가임기 여성에서 발생하는 희귀 진행성 낭성 폐질환입니다.  
TSC1 또는 TSC2 유전자의 기능 소실 변이로 mTORC1이 과활성화되어 LAM 세포(평활근 유사 신생물 세포)가 폐 실질을 침윤하고 파괴하며, 양측 폐에 얇은 벽의 낭종을 형성합니다.

| 항목 | 내용 |
|------|------|
| 유병률 | 여성 10만 명당 1-9명 (전세계 약 3,000-6,000명) |
| 주요 발병 연령 | 가임기 여성 (20-40대) |
| 유전 배경 | TSC2 변이(散發性 80%), TSC1/2 배선 변이(TSC-LAM) |
| 핵심 경로 | TSC1/2 → Rheb-GTP → mTORC1 과활성 |
| 진단 바이오마커 | 혈청 VEGF-D >800 pg/mL (민감도 ~73%, 특이도 ~100%) |
| 주요 동반질환 | 신장 혈관근지방종(AML, 50-60%), 기흉, 유미흉 |
| 1차 치료 | 시롤리무스 2mg/day (ERS 2022, ATS 2017) |
| 예후 | 연간 FEV1 감소 ~120 mL/yr (미치료); 치료 시 안정화 |

---

## QSP 모델 핵심 기전 (Key Mechanisms)

### 1. TSC1/TSC2 이중 타격 모델
- **TSC-LAM**: TSC1 또는 TSC2 배선 변이 + 체세포 LOH → TSC 복합체 완전 소실
- **散發性 LAM**: TSC2 체세포 변이(de novo) → tuberin(TSC2 단백) 기능 소실
- TSC 복합체(hamartin-tuberin)는 Rheb의 GTPase 활성화 단백(GAP) — 소실 시 Rheb-GTP ↑↑

### 2. mTOR 신호 과활성 경로
```
TSC2 기능 소실
    ↓
Rheb-GTP↑ (정상의 ~2.5배)
    ↓
mTORC1 과활성 (정상의 ~4배)
    ↓
S6K1-pT389↑ / 4E-BP1-P↑
    ↓
LAM 세포 증식 · 생존 · 침윤 ↑
```

### 3. LAM 세포 생물학
- **기원**: 자궁/골반 LAM 세포 → 혈행/림프 전파 → 폐 정착
- **표현형**: HMB-45+, α-SMA+, desmin+, ERα+, PR+ (평활근 + 멜라닌세포 특성)
- **침윤**: MMP-2, MMP-9, MMP-13 분비 → ECM 분해 → 기도 주변 낭종 형성
- **VEGF-D**: 혈청 VEGF-D 분비 → VEGFR-3 → 림프관신생 → 유미흉/후복막 림프관종

### 4. 에스트로겐 조절
- ERα 발현 → E2가 PI3K/Akt 비게놈 경로 활성 → TSC2 Akt-인산화 억제 → mTORC1 추가 상향
- 폐경 후 진행 속도 감소, 임신 시 악화 경향 → 에스트로겐 의존성 확인

### 5. 약물 PK/PD (mTOR 억제제)

| 약물 | 기전 | 표적 혈중 농도 | 주요 임상 근거 |
|------|------|----------------|----------------|
| **시롤리무스** 2mg/day | FKBP12 결합 → mTORC1 알로스테릭 억제 | 5-15 ng/mL (trough) | MILES (NEJM 2011) |
| **에버롤리무스** 10mg/day | FKBP12 결합 → mTORC1 억제 | 5-10 ng/mL | EXIST-2 (Lancet 2016) |

---

## ODE 모델 구조 (18 compartments)

| 번호 | 구획 | 설명 |
|------|------|------|
| 1-3 | SIRO_GUT/C/P | 시롤리무스 PK (2-구획 구강) |
| 4-6 | EVER_GUT/C/P | 에버롤리무스 PK |
| 7 | RHEB_GTP | Rheb-GTP 분율 |
| 8 | MTORC1 | mTORC1 활성 (정규화) |
| 9 | S6K1_P | S6K1-pT389 인산화 |
| 10 | EBPP1 | 4E-BP1 인산화 |
| 11 | LAM_CELLS | LAM 세포 부담 |
| 12 | VEGFD | 혈청 VEGF-D (pg/mL) |
| 13 | MMP_ACT | MMP 활성 (정규화) |
| 14 | ESTROGEN | 에스트로겐 수준 |
| 15 | CYST_VOL | 폐 낭종 부피 (%) |
| 16 | FEV1_PCT | FEV1 (%예측치) |
| 17 | DLCO_PCT | DLCO (%예측치) |
| 18 | AML_VOL | 신장 AML 부피 (mL) |

---

## 치료 시나리오 (5개)

1. **무치료** — 자연 경과 (~120 mL/yr FEV1 감소)
2. **시롤리무스 2mg/day** — MILES 임상시험 보정
3. **에버롤리무스 10mg/day** — EXIST-2 보정 (AML 50% 감소)
4. **시롤리무스 12개월 후 중단** — MILES 중단 후 FEV1 재감소 모사
5. **에버롤리무스 + GnRH 작용제** — mTOR 억제 + 에스트로겐 차단 병용

---

## 임상시험 보정 (Clinical Calibration)

| 임상시험 | 핵심 결과 | 모델 반영 |
|---------|---------|---------|
| **MILES** (McCormack, NEJM 2011) | 시롤리무스: FEV1 +153mL vs 대조군, VEGF-D ~30%↓ | FEV1 안정화, VEGF-D 감소 PD |
| **Bissler et al.** (NEJM 2008) | 시롤리무스: AML -47%, 중단 후 재증가 | AML ODE 구획 |
| **Johnson et al.** (NEJM 2010) | 자연 경과: FEV1 ~117 mL/yr 감소 | kFEV1_decline 파라미터 |
| **EXIST-2** (Kingswood, Lancet 2016) | 에버롤리무스: AML >50% 감소 | kAML_shrink 파라미터 |
| **Young et al.** (Ann Int Med 2011) | VEGF-D >800 pg/mL: 민감도 73%, 특이도 100% | VEGFD_LAM=1500 초기값 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| [`lam_qsp_model.dot`](lam_qsp_model.dot) | Graphviz 기계론적 지도 소스 |
| [`lam_qsp_model.svg`](lam_qsp_model.svg) | 벡터 기계론적 지도 |
| [`lam_qsp_model.png`](lam_qsp_model.png) | 래스터 기계론적 지도 (150 dpi) |
| [`lam_mrgsolve_model.R`](lam_mrgsolve_model.R) | 18구획 ODE PK/PD 모델 |
| [`lam_shiny_app.R`](lam_shiny_app.R) | 6탭 인터랙티브 Shiny 대시보드 |
| [`lam_references.md`](lam_references.md) | PubMed 참고문헌 50편 |

---

## Shiny 앱 탭 구조 (6 tabs)

| 탭 | 주요 내용 |
|----|---------|
| 1. Patient Profile | 환자 프로파일, 중증도 게이지, 질환 개요 |
| 2. Drug PK | 농도-시간 곡선, 정상상태 PK, PK 파라미터 표 |
| 3. mTOR Pathway (PD) | mTORC1 억제율, S6K1-pT389, 4E-BP1, Rheb-GTP |
| 4. Clinical Endpoints | FEV1, DLCO, 낭종 부피, 추정 6분 보행 거리 |
| 5. Scenario Comparison | 5개 치료 시나리오 비교, 12/24개월 결과 표 |
| 6. Biomarker Dashboard | VEGF-D, S6K1, MMP, AML 부피 패널 |

---

## 기계론적 지도 클러스터 (14 subgraphs, 120+ nodes)

1. Genetic Basis (TSC1/TSC2)
2. LAM Cell Origin & Phenotype
3. Upstream mTOR Signaling Inputs (PI3K/Akt, AMPK, HIF-1α)
4. TSC/Rheb/mTOR Core Axis
5. mTORC1 Downstream Substrates (S6K1, 4E-BP1, autophagy, TFEB)
6. LAM Cell Biology (proliferation, migration, invasion, MMP)
7. Lung Pathology (cysts, ECM degradation, airflow obstruction)
8. Lymphatic & Extra-Pulmonary Involvement (VEGF-D, AML, chylothorax)
9. Hormonal Regulation (E2, ERα, GnRH)
10. Drug PK (Sirolimus & Everolimus)
11. Drug PD (mTOR inhibition, PD biomarkers)
12. Adverse Effects (stomatitis, pneumonitis, hyperlipidemia)
13. Disease Biomarkers (VEGF-D, S6K1, CT score, spirometry)
14. Clinical Endpoints & Management (FEV1, transplant, guidelines)

---

## 참고문헌 요약

- McCormack FX et al. **MILES Trial**. *N Engl J Med* 2011;364:1595-1606.
- Kingswood JC et al. **EXIST-2**. *Lancet* 2016;387:1629-1638.
- Young LR et al. **VEGF-D Biomarker**. *Ann Intern Med* 2011;154:743-751.
- Johnson SR et al. **Natural History**. *N Engl J Med* 2010;363:950-959.
- Gupta N et al. **ATS/JRS Guidelines**. *Am J Respir Crit Care Med* 2017;196:1337-1348.

전체 참고문헌 50편: [`lam_references.md`](lam_references.md)
