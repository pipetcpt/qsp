# Lambert-Eaton Myasthenic Syndrome (LEMS) QSP Model

> **분류**: 자가면역 신경근육질환 (Autoimmune Neuromuscular Disease)  
> **약어**: LEMS  
> **모델 버전**: v1.0 | 2026-06-27

---

## 질환 개요 (Disease Overview)

**람베르트-이튼 근무력증후군 (Lambert-Eaton Myasthenic Syndrome)**은 전시냅스 신경근육 접합부(NMJ)의 P/Q형 전압개폐칼슘채널(VGCC)을 표적으로 하는 자가항체에 의해 유발되는 희귀 자가면역 질환입니다.

### 핵심 병태생리

```
SCLC 종양 발현 VGCC
        ↓  (분자 모방, molecular mimicry)
  항VGCC IgG 자가항체 생성
        ↓  (P/Q형 VGCC 차단 및 내재화)
  전시냅스 Ca²⁺ 유입 감소
        ↓
  ACh 소포 융합 및 방출 감소
        ↓
  종판전위(EPP) 진폭 감소 → 안전계수(Safety Factor) < 1
        ↓
  근육 활동전위 불발 → 근력 저하 / CMAP 감소
```

### 주요 임상 특징
| 특징 | 설명 |
|------|------|
| **근위부 하지 근력 저하** | 상지보다 하지 우세 |
| **무반사/반사저하** | 건반사 감소 또는 소실 |
| **자율신경 기능 이상** | 구강건조, 변비, 발기부전, 기립저혈압 |
| **촉진 현상 (Facilitation)** | 단기 운동 후 근력 및 CMAP 일시적 개선 |
| **항VGCC 항체** | P/Q형 VGCC에 대한 자가항체 양성 (>85%) |

### 역학 (Epidemiology)
- 유병률: 약 3.4/100만 명
- 약 50-60%: 소세포폐암(SCLC) 동반 (부종양 LEMS)
- 나머지: 비종양성 자가면역 LEMS
- HLA-B8, DR3 연관 (비종양성 LEMS)

---

## 기계론적 지도 (Mechanistic Map)

[![LEMS QSP 기계론적 지도](lems_qsp_model.png)](lems_qsp_model.svg)

> 클릭하면 벡터 형식(SVG)의 고해상도 맵이 열립니다.

### 포함된 서브그래프 (12개 클러스터, 100+ 노드)

| 클러스터 | 내용 |
|----------|------|
| SCLC 종양 | 종양 성장, VGCC 항원 발현, 부종양 기전 |
| 종양학 치료 | 항암화학요법, 면역항암제, 방사선치료 |
| 면역계 활성화 | DC, CD4+ Th2, Tfh, 배중심 반응, B세포, 형질세포 |
| VGCC 자가항체 | 항체 동태, FcRn 재순환, 보체 활성화 |
| 전시냅스 말단 | SNARE 복합체, ACh 소포 풀, Ca²⁺ 동태, 촉진 |
| NMJ/후시냅스 | nAChR, EPP, 안전계수, 근육 AP, 수축 |
| 자율신경계 | 교감/부교감 VGCC 차단 → 자율신경 기능이상 |
| 아미파리딘 PK | 2구획 PK (흡수, 중앙, 말초) |
| 면역억제제 PK | 프레드니솔론, 아자티오프린, MMF PK |
| 아미파리딘 PD | K⁺ 채널 차단 → AP 연장 → Ca²⁺ 유입 증가 |
| 면역억제 PD | GR 활성화, NF-κB 억제, B세포 억제 |
| 임상 엔드포인트 | CMAP, QMG 점수, MRC, 촉진비, 안전계수 |

---

## mrgsolve ODE 모델 (`lems_mrgsolve_model.R`)

### 모델 구획 (15개 ODE 구획)

```r
CMT: A_gut, A_central, A_periph          # 아미파리딘 PK (2-구획)
     P_gut, P_central, P_periph          # 프레드니솔론 PK (2-구획)
     Ab_VGCC                             # 항VGCC 항체 동태
     VGCC_free                           # VGCC 기능적 분율
     RRP                                 # ACh 즉시방출 소포 풀
     EPP_amp                             # 종판전위 진폭
     CMAP                                # CMAP 진폭
     QMG                                 # QMG 점수
     Bcell                               # B세포 동태
     Tumor                               # SCLC 종양 크기
     Facil                               # 촉진 상태 변수
```

### 주요 방정식 요약

| 방정식 | 설명 |
|--------|------|
| `VGCC_free' = k_recov×(1-f) - k_block×Ab×f` | 항체 의존성 VGCC 차단 |
| `Ab_VGCC' = kin×Bcell×(1-Imax_pred) - kout×Ab` | 항체 생성/소거 |
| `ACh_release ∝ Ca_pre^Hill_Ca` | Ca²⁺ 의존성 소포 방출 (Hill 방정식) |
| `CMAP = f(EPP/EPP_thresh)` | 안전계수 기반 CMAP |
| `K_block = Emax×C^n/(EC50^n + C^n)` | 아미파리딘 PD (Emax 모델) |

### 치료 시나리오 (6개)

1. **자연경과 (무치료)**: 항VGCC 항체 고수준 → CMAP 지속 감소
2. **아미파리딘 단독요법** (15 mg TID): 빠른 CMAP 개선 (1-2일 내)
3. **프레드니솔론 단독** (40 mg/일): 느린 항체 감소, 장기적 CMAP 개선
4. **병합요법** (아미파리딘 + 프레드니솔론): 조기+장기 효과 결합
5. **혈장교환 + 아미파리딘**: 빠른 항체 제거 후 아미파리딘으로 유지
6. **부종양 LEMS — 항암 + 아미파리딘**: 종양 축소 → 항원 감소 → 항체 감소

---

## Shiny 앱 (`lems_shiny_app.R`)

### 탭 구성 (7개 탭)

| 탭 | 주요 기능 |
|----|-----------|
| **환자 프로파일** | 항VGCC 항체 초기값, 부종양 여부, PE 세션 수 설정 |
| **약물 PK** | 아미파리딘/프레드니솔론 PK 프로파일, K⁺ 채널 차단 Emax 곡선 |
| **VGCC & 항체** | 항체 동태, VGCC 기능 상태, 항체-CMAP 산점도 |
| **NMJ / CMAP** | CMAP 시계열, 안전계수, EPP 진폭 |
| **임상 엔드포인트** | QMG 점수, 촉진 현상, 임상 반응 요약표 |
| **시나리오 비교** | 6개 치료 시나리오 동시 비교 시각화 |
| **바이오마커** | 항체-CMAP 관계, 용량-반응, B세포 동태, 바이오마커 요약표 |

### 실행 방법

```r
library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)

shiny::runApp("lems_shiny_app.R")
```

---

## 파라미터 보정 (Parameter Calibration)

| 파라미터 | 값 | 참고문헌 |
|----------|-----|----------|
| 아미파리딘 Ka | 0.693 h⁻¹ (Tmax ~1h) | Haroldsen et al. (2015) |
| 아미파리딘 CL | 18.0 L/h | Oh et al. (2009) |
| 아미파리딘 EC50 (K⁺ 차단) | 120 ng/mL | Sanders et al. (2018) |
| IgG 반감기 | ~14일 (kout=0.05 h⁻¹) | FcRn 모델 표준값 |
| VGCC 차단 속도 | 0.001 L/pmol/h | Nagel et al. (1988) |
| EPP 안전계수 역치 | 1.0 (EPP/threshold) | Wood & Slater (2001) |
| CMAP (정상) | 5.0 mV | 표준 전기진단 범위 |
| QMG 최대 | 39점 | LEMS 임상시험 척도 |

---

## 파일 목록 (File Index)

| 파일 | 설명 |
|------|------|
| `lems_qsp_model.dot` | Graphviz 기계론적 지도 소스 (12 클러스터, 100+ 노드) |
| `lems_qsp_model.svg` | 벡터 형식 기계론적 지도 |
| `lems_qsp_model.png` | 래스터 형식 기계론적 지도 (150 dpi) |
| `lems_mrgsolve_model.R` | mrgsolve ODE PK/PD 모델 (15구획, 6 시나리오) |
| `lems_shiny_app.R` | Shiny 인터랙티브 대시보드 (7탭) |
| `lems_references.md` | 참고문헌 (56개 PubMed 링크) |
| `README.md` | 이 파일 |

---

## 약어 목록 (Abbreviations)

| 약어 | 설명 |
|------|------|
| LEMS | Lambert-Eaton Myasthenic Syndrome |
| VGCC | Voltage-Gated Calcium Channel |
| SCLC | Small Cell Lung Cancer |
| NMJ | Neuromuscular Junction |
| CMAP | Compound Muscle Action Potential |
| EPP | Endplate Potential |
| DAP | Diaminopyridine (3,4-DAP = amifampridine) |
| QMG | Quantitative Myasthenia Gravis score |
| RRP | Readily Releasable Pool (ACh vesicles) |
| PK/PD | Pharmacokinetics/Pharmacodynamics |
| PE | Plasma Exchange |
| IVIG | Intravenous Immunoglobulin |
| GR | Glucocorticoid Receptor |
| FcRn | Neonatal Fc Receptor |
| RNS | Repetitive Nerve Stimulation |

---

*모델 생성: Claude Code Routine (CCR) | 날짜: 2026-06-27*
