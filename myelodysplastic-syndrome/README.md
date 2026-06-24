# Myelodysplastic Syndrome (MDS) QSP Model
## 골수이형성 증후군 정량적 시스템 약리학 모델

---

## Disease Overview / 질환 개요

**English:**
Myelodysplastic Syndrome (MDS) is a heterogeneous group of clonal hematopoietic stem cell disorders characterized by ineffective hematopoiesis, peripheral blood cytopenias, and risk of transformation to acute myeloid leukemia (AML). MDS affects approximately 3-5 per 100,000 persons annually, with incidence rising sharply with age (median age 70-75 years). The 2022 WHO classification distinguishes key subtypes: MDS with low blasts (MDS-LB), MDS with increased blasts 1/2 (MDS-IB1, MDS-IB2), MDS with ring sideroblasts (MDS-RS), MDS with del(5q), and others.

**한국어:**
골수이형성 증후군(MDS)은 클론성 조혈모세포 장애로, 무효 조혈, 말초혈액 혈구감소증, 급성 골수성 백혈병(AML)으로의 전환 위험을 특징으로 하는 이질적인 질환군입니다. MDS는 연간 10만 명당 약 3-5명에게 발생하며, 고령(중앙 연령 70-75세)에서 발생률이 급격히 증가합니다. 2022 WHO 분류는 MDS-LB, MDS-IB1/IB2, MDS-RS, MDS-del(5q) 등 주요 아형을 구분합니다.

---

## Key Pathophysiology / 핵심 병태생리

**Molecular Mechanisms:**
- **Splicing factor mutations**: SF3B1 (ring sideroblasts), SRSF2, U2AF1, ZRSR2
- **Epigenetic dysregulation**: TET2, DNMT3A, IDH1/2, ASXL1, EZH2
- **Transcription factors**: RUNX1, TP53, ETV6
- **Chromosomal abnormalities**: del(5q), monosomy 7, complex karyotype

**Ineffective Hematopoiesis Mechanisms:**
- Excess apoptosis in BM progenitors (TNF-α, IFN-γ, Fas/TRAIL)
- GDF11/GDF15 overproduction → Smad2/3 activation → erythroid differentiation block
- TGF-β superfamily dysregulation → ineffective erythropoiesis
- Hepcidin dysregulation → iron overload (from transfusions)

**Clonal Evolution:**
- CHIP → clonal cytopenia of undetermined significance (CCUS) → MDS → AML
- Clonal selection under hypomethylating agent pressure

---

## Model Files / 모델 파일

| File | Description | Size |
|------|-------------|------|
| `mds_qsp_model.dot` | Graphviz mechanistic map source | ~10 KB |
| `mds_qsp_model.svg` | Rendered vector image | ~200 KB |
| `mds_qsp_model.png` | Rendered raster image (150 dpi) | ~500 KB |
| `mds_mrgsolve_model.R` | mrgsolve ODE model (18 compartments) | ~15 KB |
| `mds_shiny_app.R` | Shiny interactive dashboard (6 tabs) | ~25 KB |
| `mds_references.md` | 35+ curated PubMed references | ~8 KB |

---

## Model Specifications / 모델 사양

| Component | Specification |
|-----------|---------------|
| Mechanistic Map Nodes | 130+ |
| Mechanistic Map Clusters | 10 |
| ODE Compartments | 18 |
| Drug PK Compartments | 6 (AZA, DEC, LEN, LUSP, EPO, VEN) |
| Disease PD Variables | 12 |
| Treatment Scenarios | 7 |
| Shiny App Tabs | 6 |
| References | 35+ |

---

## Treatment Scenarios / 치료 시나리오

| # | Treatment | Target Population | Key Trial | Primary Endpoint |
|---|-----------|------------------|-----------|-----------------|
| 1 | BSC (Best Supportive Care) | All MDS | - | Disease progression baseline |
| 2 | AZA 75 mg/m² SC d1-7 q28d | Int/High/VHigh risk | AZA-001 (Fenaux 2009) | OS benefit vs CCR |
| 3 | DEC 20 mg/m² IV d1-5 q28d | Int/High risk | DACO-020 | CR + PR rate |
| 4 | Oral-DEC/Cedazuridine 35/100 mg d1-5 q28d | Int/High risk | ASTX727 (Garcia-Manero 2020) | AUC bioequivalence to IV DEC |
| 5 | Lenalidomide 10 mg QD d1-21 q28d | del(5q) MDS | MDS-003 (List 2006) | TI rate 67% |
| 6 | Luspatercept 1.0 mg/kg SC q21d | MDS-RS, TD | COMMANDS (Platzbecker 2023) | TI 59% vs 31% (ESA) |
| 7 | VEN 400 mg QD + AZA 75 mg/m² | High-risk/AML transition | VIALE-A (DiNardo 2020) | CR rate 24% |

---

## Mechanistic Map Preview / 기계론적 지도 미리보기

[![MDS QSP Mechanistic Map](mds_qsp_model.png)](mds_qsp_model.svg)

---

## IPSS-R Risk Stratification / IPSS-R 위험도 분류

| Risk Group | Score | Median OS | AML 25% |
|------------|-------|-----------|---------|
| Very Low | ≤1.5 | 8.8 years | NR |
| Low | >1.5-3 | 5.3 years | 10.8 years |
| Intermediate | >3-4.5 | 3.0 years | 3.2 years |
| High | >4.5-6 | 1.6 years | 1.4 years |
| Very High | >6 | 0.8 years | 0.7 years |

---

## Usage Instructions / 사용 방법

### Running the mrgsolve Model / mrgsolve 모델 실행

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))

# Run model
source("mds_mrgsolve_model.R")
```

### Running the Shiny App / Shiny 앱 실행

```r
# Install dependencies
install.packages(c("shiny", "mrgsolve", "dplyr", "ggplot2", "plotly", "bslib", "DT"))

# Launch app
shiny::runApp("mds_shiny_app.R")
```

### Rendering the Mechanistic Map / 기계론적 지도 렌더링

```bash
# SVG (vector)
dot -Tsvg mds_qsp_model.dot -o mds_qsp_model.svg

# PNG (150 dpi)
dot -Tpng -Gdpi=150 mds_qsp_model.dot -o mds_qsp_model.png
```

---

## Key Clinical Findings Modeled / 모델링된 주요 임상 소견

1. **AZA vs CCR**: AZA significantly extends OS in high-risk MDS (24.5 vs 15.0 months, AZA-001 trial)
2. **Luspatercept superiority**: COMMANDS trial showed TI rate 59% vs 31% for ESA-naive MDS-RS patients
3. **Lenalidomide del5q**: 67% TI rate and 45% cytogenetic remission (MDS-003)
4. **Iron overload**: Each RBC transfusion adds ~200mg iron; ferritin >1000 ng/mL in transfusion-dependent patients
5. **AML transformation**: Risk increases exponentially with blast % and IPSS-R score
6. **VEN+AZA**: 24% CR rate in AML; data in high-risk MDS emerging

---

## 참고문헌 / References

주요 참고문헌은 `mds_references.md` 파일을 참조하세요.
See `mds_references.md` for complete bibliography (35+ PubMed citations).

---

*Generated by Claude Code Routine (CCR) | Date: 2026-06-23*
