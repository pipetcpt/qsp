# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Model Library

| # | Date | Category | Disease | Mechanistic Map | mrgsolve | Shiny App | References |
|---|------|----------|---------|-----------------|----------|-----------|------------|
| 1 | 2026-06-19 | 만성질환 | **복부 대동맥류 (Abdominal Aortic Aneurysm, AAA)** | [![AAA](abdominal-aortic-aneurysm/aaa_qsp_model.png)](abdominal-aortic-aneurysm/aaa_qsp_model.svg) | [aaa_mrgsolve_model.R](abdominal-aortic-aneurysm/aaa_mrgsolve_model.R) | [aaa_shiny_app.R](abdominal-aortic-aneurysm/aaa_shiny_app.R) | [aaa_references.md](abdominal-aortic-aneurysm/aaa_references.md) |

---

## Abdominal Aortic Aneurysm (AAA) — Model Details

### 질환 개요 (Disease Overview)
복부 대동맥류(AAA)는 복부 대동맥 직경이 3 cm 이상으로 확장된 상태로, 파열 시 사망률이 80% 이상에 달하는 치명적인 혈관 질환입니다. 주요 위험인자는 흡연, 고령(>65세), 남성, 가족력, 고혈압입니다.

### 병태생리 핵심 (Key Pathophysiology)
1. **MMP 매개 ECM 분해**: MMP-2(젤라티나아제 A), MMP-9(젤라티나아제 B), MMP-12(대식세포 엘라스타아제)가 탄력소(elastin)와 콜라겐을 분해하여 대동맥 벽 구조를 약화시킴
2. **대식세포 극성화**: M1 대식세포가 TNF-α, IL-1β, IL-6 등 염증성 사이토카인과 MMP를 분비하여 염증 반응을 증폭
3. **혈관 평활근세포(VSMC) 소실**: MMP-9, TNF-α, 활성산소종(ROS)에 의한 VSMC 세포사멸(anoikis)이 대동맥 벽 완전성을 손상
4. **산화 스트레스**: NADPH 산화효소가 생성하는 초과산화물이 ROS 풀을 형성하여 MMP 활성화, NF-κB 신호 증폭
5. **내강내 혈전(ILT)**: 난류 혈류로 혈소판 활성화 → 피브린 형성 → ILT 성장 → 국소 저산소증 및 추가 MMP 방출

### 약물 PK/PD (Drug PK/PD Mechanisms)
| 약물 | 기전 | 주요 파라미터 | 임상 근거 |
|------|------|--------------|-----------|
| **독시사이클린 (Doxycycline)** | MMP-2, MMP-9, MMP-12 직접 억제 (zinc chelation) | Imax=0.85, IC50=0.20 mg/L (조직), t½=18h | PHAST-1/2 (Lindeman 2009, Meijer 2013) |
| **스타틴 (Simvastatin)** | HMG-CoA 억제 + 다면발현 효과 (NF-κB↓, ROS↓, MMP-9↓) | MMP-9 Imax=0.40, IC50=0.05 mg/L | Brady 2004, Meijer 2012 |
| **프로프라놀롤 (Propranolol)** | β1 차단 → 심박수↓, 혈압↓ → 대동맥 벽 스트레스 감소 | BP Imax=0.20, IC50=0.02 mg/L | UK Small Aneurysm Trial 2002 |

### 모델 구조 (Model Structure)
- **ODE 구획**: 총 20개 (PK 9개 + PD 11개)
  - 독시사이클린: GUT → CENTRAL → PERIPHERAL + AORTIC_TISSUE (4구획)
  - 스타틴: GUT → LIVER → CENTRAL (3구획)
  - 프로프라놀롤: GUT → CENTRAL → PERIPHERAL (3구획)
  - PD: MAC · TNF · ROS · MMP9 · MMP2 · ELASTIN · COLLAGEN · VSMC · ILT · DIAMETER (10구획)
- **치료 시나리오**: 6개 (무치료, 독시사이클린, 스타틴, 프로프라놀롤, 복합 2제, 3제 병용요법)

### 기계론적 지도 특징 (Mechanistic Map Features)
- **노드 수**: 130+ (13개 서브그래프 클러스터)
- **클러스터**: 위험인자 · 혈역학적 스트레스 · 내피세포 기능이상 · 면역세포 침윤 · 사이토카인 네트워크 · MMP/TIMP 축 · 산화 스트레스 · VSMC 생물학 · ECM 재형성 · 내강내 혈전 · 임상 결과 · 약물 PK · 약물 PD
- **주요 표적 강조**: Active MMP-2, Active MMP-9, VSMC Apoptosis, Rupture Risk

### Shiny 앱 탭 구성 (Shiny App Tabs)
1. **Patient Profile** — 환자 특성, 위험 점수, AAA 병기 분류표
2. **Drug PK** — 독시사이클린/스타틴/프로프라놀롤 혈장 및 조직 농도 시뮬레이션
3. **MMP Biomarkers** — MMP-9/2, TNF-α, ROS, 대식세포 활성도 동역학
4. **Aortic Wall Dynamics** — 대동맥 직경, 탄력소/콜라겐 함량, VSMC 밀도, ILT 용적
5. **Scenario Comparison** — 6개 치료 시나리오 직접 비교 (직경, MMP-9, 파열 위험)
6. **Rupture Risk Assessment** — 현재 환자 상태 기반 파열 위험 추정 + 5년 장기 예측

### 핵심 임상시험 (Key Clinical Trials)
- **PHAST-1** (Lindeman 2009, *Circulation*): 독시사이클린이 AAA 조직의 호중구 및 세포독성 T세포를 선택적으로 감소시킴
- **PHAST-2** (Meijer 2013, *Ann Intern Med*): 독시사이클린 100mg/day 18개월 투여 → MMP-9 유의미하게 감소 (그러나 직경 성장 억제는 통계적 유의성 미달)
- **UK Small Aneurysm Trial** (2002, *NEJM*): 직경 4.0–5.5cm AAA에서 즉각 수술 vs 감시 비교

