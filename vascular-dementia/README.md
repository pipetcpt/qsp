# 혈관성 치매 (Vascular Dementia) QSP 모델

> **분류**: 신경계 / 뇌혈관 질환  
> **약어**: VaD  
> **빌드 날짜**: 2026-06-27  
> **모델 버전**: 1.0.0

---

## 질환 개요 (Disease Overview)

**혈관성 치매(Vascular Dementia, VaD)**는 뇌혈관 질환(cerebrovascular disease)에 의한 뇌 손상으로 발생하는 인지 장애로, **알츠하이머병 다음으로 흔한 치매의 원인**(전체 치매의 약 20–30%)입니다. 소혈관 질환(small vessel disease)에 의한 백질변성(WMH)과 열공 경색(lacunar infarct), 대혈관 죽상경화증, 심장색전증에 의한 피질 경색 등 다양한 병리 아형을 포함합니다.

### 핵심 병태생리 경로 (Key Pathophysiological Pathways)

```
고혈압/당뇨/이상지질혈증 → 소혈관 질환(SVD)
     ↓
 백질변성(WMH) + 열공경색 + 피질 미세경색
     ↓
혈뇨장벽(BBB) 손상 → DAMPs 유입
     ↓
미세아교세포 M1 활성화 → NF-κB → IL-1β/TNF-α/IL-6
     ↓
산화 스트레스(ROS) → 뉴런 손상
     ↓
콜린성 기능 저하 + 시냅스 소실
     ↓
인지 저하 (MMSE↓, 실행 기능↓, 처리 속도↓)
```

---

## 모델 구성 (Model Structure)

### 1. 기계론적 지도 (Mechanistic Map)

| 항목 | 세부 |
|------|------|
| 노드 수 | **114개** |
| 서브그래프 클러스터 | **9개** |
| 클러스터 목록 | Drug PK/PD · 혈관 위험인자 · 뇌혈관 병리 · 혈뇨장벽 · 신경염증 · 산화 스트레스 · 신경전달물질 · 뇌 구조 · 임상 결과 |

[![기계론적 지도 미리보기](vad_qsp_model.png)](vad_qsp_model.svg)

### 2. mrgsolve ODE 모델

| 항목 | 세부 |
|------|------|
| ODE 구획 수 | **18개** |
| 약물 PK 구획 | 7개 (AHT depot/central, APT, Statin, AChEI brain, Memantine brain, Cilostazol) |
| 생리·PD 구획 | 11개 (BP, LDL, CBF, WMH, Infarct, Microglia, Cytokine, ROS, ACh, Synapse, MMSE) |
| 치료 시나리오 | **6개** |
| 임상시험 보정 | PROGRESS, SCOPE, SPRINT-MIND, SIGNAL2, Black et al. 2003 (donepezil) |

**6가지 치료 시나리오:**

| # | 시나리오 | 약물 |
|---|---------|------|
| 1 | 무치료 | — |
| 2 | 강압제 단독 | 강압제 (ARB/ACEi) |
| 3 | 혈관 위험인자 병합 | 강압제 + 항혈소판제 + 스타틴 |
| 4 | 증상 치료 | AChEI (도네페질) + 메만틴 |
| 5 | 포괄적 치료 | 혈관 + 증상 + 실로스타졸 |
| 6 | 최적+ (고용량 스타틴) | 시나리오 5 + 스타틴 80 mg/day |

### 3. Shiny 대시보드

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | MMSE, WMH, BP, CBF, LDL, ACh, 시냅스 밀도 기저값 설정 및 위험인자 요약 |
| 2. Drug PK/PD | 혈중 농도-시간 프로파일 (48h), 정상상태 PD 효과 요약 |
| 3. 혈관·뇌관류 | SBP, CBF, WMH 진행, 미세경색 누적 시뮬레이션 |
| 4. 신경생물학적 기전 | 미세아교세포/사이토카인, ROS, ACh 톤, 시냅스 밀도 |
| 5. 임상 엔드포인트 | MMSE 궤적, 인지 도메인, VaD 단계 분류 |
| 6. 시나리오 비교 | 6가지 치료 전략 MMSE/WMH/CBF 2년 비교 |
| 7. 바이오마커 | MRI/CSF/혈액 바이오마커 패널 (정상범위 vs VaD 범위) |

### 4. 참고문헌

| 항목 | 세부 |
|------|------|
| 총 참고문헌 수 | **70개** |
| 분류 | 진단/역학 · 소혈관 질환 · 혈뇨장벽 · 신경염증 · 산화 스트레스 · 콜린성 · 임상시험(BP/스타틴/항혈소판/실로스타졸/도네페질/메만틴) · QSP 방법론 · 바이오마커 · 예방/위험인자 |

---

## 모델 파라미터 주요 근거 (Key Parameter Calibration)

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| Emax_AHT_BP | 25 mmHg | PROGRESS trial (perindopril-based Rx) |
| WMH 진행 속도 | ~0.65 mL/yr (untreated) | PROGRESS MRI substudy, Dufouil 2005 |
| AChEI 효과 | +1.0 MMSE @ 24주 | Black et al. 2003 (donepezil in VaD) |
| 실로스타졸 WMH 감소 | −11.5% | SIGNAL2 trial |
| E_ST_LDL (스타틴) | 50% LDL 감소 | Multiple statin RCTs |
| Emax_AChEI | 75% AChE 억제 | Pharmacodynamic modeling |
| CBF 기저값 (VaD) | ~45-55 mL/100g/min | ASL-MRI, xenon-CT 데이터 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| `vad_qsp_model.dot` | Graphviz 기계론적 지도 소스 |
| `vad_qsp_model.svg` | SVG 벡터 이미지 |
| `vad_qsp_model.png` | PNG 래스터 이미지 (150 dpi) |
| `vad_mrgsolve_model.R` | mrgsolve ODE 모델 (18 구획, 6 시나리오) |
| `vad_shiny_app.R` | Shiny 인터랙티브 대시보드 (7 탭) |
| `vad_references.md` | 참고문헌 (70개) |
| `README.md` | 이 파일 |

---

## 실행 방법 (How to Run)

### 기계론적 지도 렌더링 (Graphviz 필요)
```bash
dot -Tsvg vad_qsp_model.dot -o vad_qsp_model.svg
dot -Tpng -Gdpi=150 vad_qsp_model.dot -o vad_qsp_model.png
```

### mrgsolve 모델 실행 (R 필요)
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("vad_mrgsolve_model.R")
```

### Shiny 앱 실행
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "DT", "tidyr"))
shiny::runApp("vad_shiny_app.R")
```

---

## 핵심 약물 작용기전 요약 (Drug Mechanism Summary)

| 약물 | 작용기전 | 주요 효과 |
|------|---------|----------|
| **강압제** (ARB/ACEi) | RAAS 억제 → SVD 예방, BBB 보호 | WMH 진행 억제, 뇌관류 개선 |
| **스타틴** | HMG-CoA 억제 + 다면효과 (eNOS↑, NF-κB↓, BDNF↑) | LDL↓, 신경염증↓, 혈관 기능 개선 |
| **항혈소판제** | TXA2 합성 억제 (COX-1) | 혈소판 응집↓, 소혈관 혈전↓ |
| **AChEI (도네페질)** | AChE 억제 → 시냅스 ACh↑ | 콜린성 톤 회복, MMSE 안정화 |
| **메만틴** | 저친화성 NMDA-R 차단 | 흥분독성↓, 시냅스 소실 억제 |
| **실로스타졸** | PDE3 억제 → cAMP↑ → 혈관 이완 | CBF↑, WMH 진행↓, 혈소판 응집↓ |

---

## 한계 및 주의 사항 (Limitations)

- 본 모델은 **교육 및 연구 목적**의 반정량적 QSP 모델입니다.
- 임상 환자 데이터에 대한 개별 보정 없이 직접 임상 적용 불가.
- 혼합형 VaD(알츠하이머 공존)의 특성은 별도 모델 확장이 필요합니다.
- 파라미터는 공개 문헌의 대표값(mean)으로, 개인 변이를 충분히 반영하지 않을 수 있습니다.

---

## 관련 모델 (Related Models in This Repository)

- [알츠하이머병](../alzheimers-disease/) — Aβ/Tau 캐스케이드 중심 모델
- [본태성 고혈압](../essential-hypertension/) — RAAS/SNS 중심 BP 조절 모델
- [뇌졸중](../ischemic-stroke/) — 급성 허혈 및 재관류 손상 모델
- [당뇨병성 신병증](../diabetic-nephropathy/) — 대사-혈관 기전 모델
- [다발성 경화증](../multiple-sclerosis/) — CNS 신경염증 모델
