# 갈색세포종/부신경절종 (Pheochromocytoma/Paraganglioma, PPGL) — QSP 모델

> **디렉토리:** `pheochromocytoma/` | **약어:** PPGL | **날짜:** 2026-06-25

[![PPGL QSP 기계론적 지도](ppgl_qsp_model.png)](ppgl_qsp_model.svg)

---

## 질환 개요 (Disease Overview)

갈색세포종(Pheochromocytoma)은 **부신 수질의 크롬친화세포(chromaffin cell)**에서 발생하는 카테콜아민 분비 신경내분비종양으로, 부신 외 교감신경 절에서 발생하면 **부신경절종(Paraganglioma)**으로 분류됩니다. 두 종양을 통칭하여 **PPGL**이라 합니다.

| 특성 | 내용 |
|------|------|
| **유병률** | 인구 10만 명당 2–8명; 우연발견 부신종양의 5% |
| **유전성** | 전체의 40%는 생식세포 변이 (SDH, VHL, RET, NF1, MAX 등) |
| **악성 비율** | 전체의 10–17%; SDHB 변이 시 최대 40–80% |
| **주요 증상** | 두통–발한–심계항진 삼주증(WHO 기준); 저혈압·고혈압·위기 발작 |
| **진단 바이오마커** | 혈장 유리 메타네프린(NMN/MN), 24h 소변 카테콜아민, 크로모그라닌-A |
| **수술 치유율** | 국소성 PPGL: R0 절제 후 생화학적 완치 80–95% |

### 병태생리 핵심

1. **카테콜아민 과분비**: TH(tyrosine hydroxylase) 의존적 생합성 경로를 통해 NE/EPI가 과량 분비 → α₁-수용체를 통한 혈관수축 → **고혈압 위기** 유발
2. **유전자 군집화**: SDHB/C/D 변이(가성저산소 군집) vs. RET/NF1/TMEM127 변이(키나제 신호 군집) — 임상 표현형과 악성도가 다름
3. **심혈관 독성**: 만성 카테콜아민 과다 → 카테콜아민 심근병증, 타코츠보 유사 심근병증, 부정맥

---

## 기계론적 지도 (Mechanistic Map)

| 파일 | 형식 |
|------|------|
| `ppgl_qsp_model.dot` | Graphviz 소스 |
| `ppgl_qsp_model.svg` | 벡터 (고해상도) |
| `ppgl_qsp_model.png` | 래스터 (150 dpi) |

### 11개 클러스터 (130+ 노드)

| 클러스터 | 주요 구성 요소 |
|---------|------------|
| ① 유전·분자 드라이버 | SDHB/C/D/A·VHL·RET·NF1·MAX·EPAS1/HIF-2α·가성저산소·키나제 클러스터·악성도 |
| ② 종양 생물학 | 크롬친화세포·종양 증식·세포자멸·HIF-1α·VEGF-A·신생혈관·CgA·NSE·전이 |
| ③ 카테콜아민 생합성 | 타이로신→TH→DOPA→AADC→도파민→DBH→NE→PNMT→EPI; COMT·MAO 대사 |
| ④ 저장·분비 | 크롬친화 과립·VMAT2·Ca²⁺ 유입·엑소사이토시스·SNARE·NET 재흡수 |
| ⑤ 아드레날린 수용체 신호 | α₁/α₂/β₁/β₂/β₃-AR·Gq/Gi/Gs·PLC→IP3/DAG→PKC·cAMP→PKA·MAPK |
| ⑥ 심혈관 효과 | SBP·DBP·MAP·HR·CO·SVR·LVEF·카테콜아민 심근병증·고혈압 위기·부정맥 |
| ⑦ 대사 효과 | 간 당분해·글루카곤↑·인슐린 억제·고혈당·FFA·BAT 열생성·BMR↑·체중감소 |
| ⑧ α-차단제 PK | 페녹시벤자민(비가역적) 2구획·독사조신(선택적 α₁) 1구획 PK |
| ⑨ 전신 치료 PK | 수니티닙 2구획·¹³¹I-MIBG·¹⁷⁷Lu-DOTATATE·CVD 화학요법 |
| ⑩ 약력학·바이오마커 | α₁/β 수용체 점유율·TH 억제·VEGFR 억제·혈장 NMN/MN·CgA·RECIST 반응 |
| ⑪ 수술·주술기 관리 | 수술 전 알파·베타 차단·메티로신·IV 수액·복강경/개복 수술·수술 중 위기·수술 후 저혈압 |

---

## mrgsolve ODE 모델

**파일:** `ppgl_mrgsolve_model.R`

### ODE 구획 (20개)

| # | 구획 | 생물학적 의미 |
|---|------|------------|
| 1–3 | `PHE_gut`, `PHE_C`, `PHE_P` | 페녹시벤자민 GI·중심·말초 구획 |
| 4 | `DOX_C` | 독사조신 중심 구획 |
| 5 | `MET_C` | 메티로신 중심 구획 |
| 6 | `BB_C` | 베타차단제(프로프라놀롤) 중심 |
| 7–8 | `SUNIT_C`, `SUNIT_P` | 수니티닙 중심·말초 (악성 PPGL) |
| 9 | `TH_act` | TH 효소 활성 (메티로신에 의해 억제) |
| 10 | `NE_store` | 크롬친화 과립 NE 저장 풀 |
| 11 | `NE_plasma` | 혈장 NE |
| 12 | `EPI_plasma` | 혈장 EPI |
| 13 | `TUMvol` | 종양 부피 (mL) |
| 14 | `VEGF_tum` | 혈장 VEGF (pg/mL) |
| 15 | `SBP` | 수축기 혈압 (mmHg) |
| 16 | `DBP` | 이완기 혈압 (mmHg) |
| 17 | `HR` | 심박수 (bpm) |
| 18 | `GLU` | 혈장 포도당 (mmol/L) |
| 19 | `FFA` | 유리지방산 (mmol/L) |
| 20 | `CgA_plasma` | 혈장 크로모그라닌-A (ng/mL) |

### 치료 시나리오 (6개)

| 시나리오 | 요법 | 임상 근거 |
|---------|------|---------|
| S0 | 무치료 (자연 경과) | — |
| S1 | **페녹시벤자민** 60 mg/d × 14일 → 수술 | Kinney 2002 J Cardiothorac Vasc Anesth |
| S2 | **독사조신** 16 mg/d × 14일 → 수술 | Shao 2016 World J Surg (meta-analysis) |
| S3 | **PHE + 메티로신 2 g/d + 프로프라놀롤** 3제 병용 | Steinsapir 1997 Arch Intern Med |
| S4 | **수니티닙** 37.5 mg/d (전이성 악성 PPGL) | Niemeijer 2014 J Clin Endocrinol Metab |
| S5 | **메티로신 단독** (수술 불가, 증상 조절) | Engelman 1968 NEJM |

### 임상 보정 (Calibration)

| 임상시험/연구 | 보정 표적 |
|-------------|---------|
| Lentschener 2009 Hypertension | 페녹시벤자민 vs 독사조신 수술 전후 혈압 조절 동등성 |
| Steinsapir 1997 | 메티로신 TH 억제 40–80%, 카테콜아민 합성 감소 |
| Niemeijer 2014 | 수니티닙 malignant PPGL: ORR 25%, 종양 안정화 |
| Averbuch 1988 | CVD 화학요법 ORR 37%, 임상 반응률 79% |

---

## Shiny 대시보드

**파일:** `ppgl_shiny_app.R`

| 탭 | 주요 내용 |
|----|---------|
| 1. 환자 프로파일 | PPGL 유형·종양 부피·기저 NE/EPI·SBP/HR·CgA·유전자 변이 |
| 2. 약물 PK | PHE/DOX/MET/BB/수니티닙 혈중농도 시뮬레이션 |
| 3. 카테콜아민 | 혈장 NE·EPI·NMN·MN 동태; TH 활성; CgA 종양 지표 |
| 4. 심혈관 | SBP·DBP·HR·α-차단 비율·대사 효과(혈당·FFA) |
| 5. 종양 (악성) | 종양 부피·VEGF·RECIST 반응 평가; 치료 옵션 비교표 |
| 6. 시나리오 비교 | 6개 시나리오 병렬 비교; 주요 지표 요약 테이블 |

---

## 실행 방법 (Usage)

```bash
# 기계론적 지도 렌더링
dot -Tsvg ppgl_qsp_model.dot -o ppgl_qsp_model.svg
dot -Tpng -Gdpi=150 ppgl_qsp_model.dot -o ppgl_qsp_model.png
```

```r
# mrgsolve 모델 실행
install.packages(c("mrgsolve","dplyr","ggplot2"))
source("ppgl_mrgsolve_model.R")

# Shiny 대시보드 실행
install.packages(c("shiny","shinydashboard","tidyr"))
shiny::runApp("ppgl_shiny_app.R")
```

---

## 참고문헌 (References)

**파일:** `ppgl_references.md` — 45편 PubMed 인용 (12개 섹션)

| 섹션 | 편수 |
|------|------|
| 랜드마크 리뷰·가이드라인 | 5 |
| 역학·유전학 | 4 |
| 카테콜아민 생합성·대사 | 5 |
| 수술 전 관리·알파차단 | 6 |
| 메티로신 PK/PD | 3 |
| 악성 PPGL 전신 치료 | 6 |
| 심혈관·혈역학 효과 | 4 |
| 생화학적 진단·바이오마커 | 3 |
| 영상·위치 확인 | 2 |
| 분자 병태생리 | 3 |
| 수니티닙 PK 모델링 | 2 |
| QSP/PK-PD 모델링 | 2 |
