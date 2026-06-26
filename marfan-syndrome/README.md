# Marfan Syndrome (마르팡 증후군) — QSP Model

**디렉토리:** `marfan-syndrome/`  
**약어 (Abbreviation):** MFS  
**추가일:** 2026-06-25  
**카테고리:** 결합조직질환 / 심혈관 유전질환

---

## 개요 (Overview)

마르팡 증후군(Marfan Syndrome, MFS)은 **FBN1 유전자(Chr 15q21.1)**의 생식세포 돌연변이로 인해 발생하는 상염색체 우성 결합조직질환이다.  
전 세계 유병률은 1/5,000~1/10,000이며, 남녀 동등하게 발생한다.  

결함 있는 피브릴린-1(fibrillin-1) 단백질은 세포외기질(ECM) 내 마이크로피브릴 네트워크를 파괴하여:

- **TGF-β 서열화 손상** → 자유 TGF-β1/2 증가 → SMAD2/3 + ERK/MAPK 과활성
- **ECM 분해** (MMP-2/9 상향 조절) → 대동맥 중막 변성
- **대동맥 근부 확장** → 대동맥판 역류, 박리 위험
- **골격 특징:** 큰 키, 거미손가락증, 척추측만증, 오목/새가슴
- **안과 특징:** 수정체 탈구(60-70%), 근시, 망막박리

---

## 기계론적 지도 (Mechanistic Map)

[![MFS QSP Mechanistic Map](mfs_qsp_model.png)](mfs_qsp_model.svg)

> 클릭하면 전체 해상도 SVG를 볼 수 있습니다.

### 14개 클러스터 구성

| 클러스터 | 내용 |
|----------|------|
| ① Genetic & Molecular Foundation | FBN1 돌연변이, 피브릴린-1 단백질, LTBP, 마이크로피브릴 |
| ② TGF-β Signalling | SMAD2/3, SMAD4, SMAD7, 정규/비정규 경로 |
| ③ MAPK/ERK Pathway | RAS-RAF-MEK-ERK, p38, JNK, PI3K/AKT |
| ④ ECM Remodelling | MMP-2/9/13, TIMP-1/2, 엘라스틴 분해, 콜라겐 분절 |
| ⑤ Vascular SMC Biology | VSMC 표현형 전환, 세포자멸사, NOX4, NF-κB |
| ⑥ Aortic Pathology | 대동맥 근부·발살바동·STJ·상행 대동맥, 라플라스 응력 |
| ⑦ Cardiac Manifestations | 대동맥판 역류, 승모판 탈출, LV 확장, 박리 |
| ⑧ Haemodynamic Parameters | HR, SBP, dP/dt_max, PWV, 맥파 속도 |
| ⑨ Skeletal System | 큰 키, 거미손가락증, 오목/새가슴, 척추측만증 |
| ⑩ Ocular System | 수정체 탈구, 근시, 망막박리, 녹내장 |
| ⑪ Other Systemic Features | 경막 확장, 기흉, 탈장, 수면무호흡 |
| ⑫ Drug PK | 아테노롤 2-구획, 로사르탄/EXP-3174 PK |
| ⑬ Drug PD | β1-차단, AT1R 차단, TGF-β 억제 기전 |
| ⑭ Clinical Endpoints & Biomarkers | Z-점수, 연간 성장률, AR 등급, 혈장 TGF-β |

**총 노드 수:** 130+ nodes | **클러스터:** 14개

---

## mrgsolve ODE 모델 (Compartmental Model)

**파일:** `mfs_mrgsolve_model.R`

### 20개 구획 (20 Compartments)

| 구획 | 설명 |
|------|------|
| DEPOT_ATN | 아테노롤 장관 흡수 구획 |
| C1_ATN | 아테노롤 중심 구획 |
| C2_ATN | 아테노롤 말초 구획 |
| DEPOT_LOS | 로사르탄 장관 흡수 구획 |
| C1_LOS | 로사르탄 중심 구획 |
| C_EXP3174 | EXP-3174 활성 대사물 |
| TGFb | 혈장 자유 TGF-β1 [ng/mL] |
| pSMAD | 인산화 SMAD2/3 (배수) |
| pERK | 인산화 ERK1/2 (배수) |
| MMP | 순환 MMP 활성 [U/mL] |
| Ao_Diam | 대동맥 근부 직경 [mm] |
| AR_Grade | 대동맥판 역류 등급 (0-4) |
| HR | 심박수 [bpm] |
| SBP | 수축기 혈압 [mmHg] |
| dPdt | dP/dt_max [mmHg/s] |
| NT_proBNP | NT-proBNP [pg/mL] |
| LVEDD | 좌심실 이완기말 직경 [mm] |
| TGFb_plasma_obs | 관측 혈장 TGF-β1 |
| Systemic_score | 겐트 전신 점수 (0-20) |

### 6가지 치료 시나리오

| 시나리오 | 치료 | 임상 근거 |
|----------|------|----------|
| 1. 무치료 | — | 자연경과 (Salim 1994; Rossig 2019) |
| 2. Atenolol 50mg QD | β1-차단 | PHN RCT — Lacro et al. NEJM 2014 |
| 3. Atenolol 100mg QD | β1-차단 (고용량) | Shores et al. NEJM 1994 |
| 4. Losartan 50mg QD | ARB | PHN RCT — Lacro et al. NEJM 2014 |
| 5. Losartan 100mg QD | ARB (고용량) | COMPARE — Radonic et al. EHJ 2010 |
| 6. Atenolol + Losartan | 병용 요법 | AIMS — Forteza et al. JACC 2016 |

### 주요 임상시험 보정 데이터

| 임상시험 | 약물 | 기간 | 주요 결과 |
|----------|------|------|----------|
| Shores 1994 NEJM | 프로프라노롤 | 10년 | 아테노롤군 vs 무치료: 대동맥 성장률 50% ↓ |
| Lacro 2014 NEJM (PHN) | 아테노롤 vs 로사르탄 | 3년 | 두 군 유사 (aortic root Z-score ↓ ~0.12) |
| Radonic 2010 EHJ (COMPARE) | 로사르탄 100mg | 2년 | 아테노롤 대비 대동맥 성장률 유의한 차이 없음 |
| Forteza 2016 JACC (AIMS) | 이르베사르탄 vs 아테노롤 | 3년 | 대동맥 성장률 동등 (0.48 vs 0.52 mm/yr) |
| Brooke 2008 NEJM | 로사르탄 (소아) | ~2.3년 | 로사르탄군 성장률 현저 ↓ (hist. control) |
| Habashi 2006 Science | 로사르탄 (Fbn1+/- 마우스) | — | 대동맥 확장 및 TGF-β 신호 완전 억제 |

---

## Shiny 대시보드 (Interactive Dashboard)

**파일:** `mfs_shiny_app.R`

### 7개 탭 구성

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | 겐트 기준표, 병태생리 개요, 기계론적 지도 미리보기 |
| ② Drug PK | 아테노롤/로사르탄/EXP-3174 농도-시간 곡선, PK 요약 |
| ③ TGF-β / Mol. PD | TGF-β1, p-SMAD2/3, p-ERK1/2, MMP 동태 |
| ④ 심혈관 엔드포인트 | 대동맥 근부 직경, Z-점수, AR 등급, HR/dP/dt, LVEDD, SBP |
| ⑤ 시나리오 비교 | 6개 치료 시나리오 병렬 비교, 5년 요약 테이블 |
| ⑥ 바이오마커 | NT-proBNP, TGF-β, 겐트 점수, 연간 성장률, 임계값 참조표 |
| ⑦ 수술 결정 지원 | 수술 역치 시각화, ESC/AHA 가이드라인, 치료별 역치 도달 시간 |

---

## 참고문헌 (References)

**파일:** `mfs_references.md` — 총 50건

| 섹션 | 인용 수 |
|------|---------|
| 1. 유전학 & 분자 병인 | 8 |
| 2. TGF-β 신호 | 6 |
| 3. 대동맥 병리 & 자연경과 | 6 |
| 4. 임상시험 — β-차단제 | 4 |
| 5. 임상시험 — ARB | 6 |
| 6. 약동학 | 5 |
| 7. MMP / ECM | 4 |
| 8. 안과 | 2 |
| 9. 골격 & 전신 | 3 |
| 10. 수술 | 3 |
| 11. QSP / 모델링 | 3 |

---

## 실행 방법 (Usage)

```bash
# 1) 기계론적 지도 렌더링 (Graphviz 필요)
dot -Tsvg mfs_qsp_model.dot -o mfs_qsp_model.svg
dot -Tpng -Gdpi=150 mfs_qsp_model.dot -o mfs_qsp_model.png
```

```r
# 2) mrgsolve ODE 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("mfs_mrgsolve_model.R")
# → 6 시나리오 × 5년 시뮬레이션 실행, 요약 테이블 출력

# 3) Shiny 대시보드 실행
install.packages(c("shiny", "shinydashboard", "DT"))
shiny::runApp("mfs_shiny_app.R")
```

---

## 파일 목록 (File Structure)

```
marfan-syndrome/
├── README.md                  ← 이 파일
├── mfs_qsp_model.dot          ← Graphviz 기계론적 지도 소스
├── mfs_qsp_model.svg          ← SVG 렌더링
├── mfs_qsp_model.png          ← PNG 렌더링 (150 dpi)
├── mfs_mrgsolve_model.R       ← mrgsolve ODE QSP 모델
├── mfs_shiny_app.R            ← Shiny 인터랙티브 대시보드
└── mfs_references.md          ← 참고문헌 (50건)
```
