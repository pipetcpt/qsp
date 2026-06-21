# 낭성 섬유증 (Cystic Fibrosis, CF) — QSP Model

## 개요

낭성 섬유증은 **CFTR(Cystic Fibrosis Transmembrane conductance Regulator)** 유전자 돌연변이로 인한 상염색체 열성 유전 질환으로, 전 세계 약 10만 명이 이환된 희귀질환입니다. CFTR 단백질의 기능 이상은 외분비선에서 두꺼운 점액 분비를 초래하며, 폐·소화관·생식기관 등 여러 장기를 침범합니다. 최근 CFTR 모듈레이터(Trikafta = ELX/TEZ/IVA)의 등장으로 치료 패러다임이 근본적으로 바뀌었습니다.

---

## 핵심 기전 (Pathophysiology)

| 단계 | 핵심 기전 | 임상 결과 |
|------|----------|---------|
| CFTR 유전자 돌연변이 | 6가지 클래스 (I–VI) 기능 이상 | CFTR 단백질 부재/기능부전 |
| ΔF508 처리 장애 (Class II) | ER 품질 관리 → ERAD (~99% 분해) | 막에 도달하는 CFTR <1% |
| ASL 탈수 | ENaC 과활성 + CFTR 분비 ↓ → 삼투압 이상 | 기도 표면 액체 붕괴 (<4 μm) |
| 점액 이동 장애 | 점액 섬모 이동속도 ↓, 점도↑ | 세균 (Pseudomonas) 집락화 |
| 만성 감염 | Pa 바이오필름 형성, 항생제 저항 | 중성구 침윤, 기도 파괴 |
| 전신 염증 | IL-8/TNF-α/NE 과다 분비 | FEV1 감소, 폐 섬유화 |

---

## 모델 구성 (Model Components)

### 1. 기계론적 지도 (`cf_qsp_model.dot`)

| 클러스터 | 주요 내용 | 노드 수 |
|---------|---------|--------|
| ① CFTR Gene & mRNA | 유전자, mRNA, 6가지 돌연변이 클래스 | 12 |
| ② CFTR Protein Processing | ER QC, 샤페론, ERAD, 교정제 결합 | 14 |
| ③ CFTR Channel (Membrane) | 골지 이동, 막 위치, 채널 개폐 | 15 |
| ④ ENaC & Ion Transport | ENaC, NKCC1, Na/K-ATPase, TMEM16A | 11 |
| ⑤ ASL & Mucociliary Clearance | PCL, 점액층, 섬모 박동, MCC | 13 |
| ⑥ Airway Inflammation | NF-κB, IL-8, TNF-α, 중성구, NETs, ROS | 25 |
| ⑦ Infection & Biofilm | P. aeruginosa, S. aureus, 바이오필름, 항생제 | 15 |
| ⑧ Lung Function | FEV1, LCI, 기관지 확장증, 이식 기준 | 13 |
| ⑨ CFTR Modulator PK/PD | Ivacaftor/ELX/TEZ 약동학 구획 | 16 |
| ⑩ Other Therapies | DNase, 고장식염수, PERT, 유전자치료 | 10 |
| ⑪ Systemic Effects | CFRD, 췌장 기능, 간, 골밀도 | 11 |
| ⑫ Clinical Endpoints | 땀 염소, NPD, ppFEV1, CFQ-R | 10 |
| **총계** | | **165 노드, 247 엣지** |

### 2. mrgsolve ODE 모델 (`cf_mrgsolve_model.R`)

**구획 (25 ODEs):**
- 약물 PK: Ivacaftor (3구획) + Elexacaftor (2구획) + Tezacaftor (2구획)
- CFTR 생물학: Band B 처리, 막 위치 CFTR
- 기도 표면 액체 (ASL)
- 감염: Pa 부유균 + 바이오필름
- 염증: IL-8, 중성구, 누적 손상 점수
- 폐 기능: ppFEV1, 누적 악화 횟수
- 전신: BMI, 췌장 기능

**치료 시나리오 (7개):**

| 시나리오 | 약물 | 돌연변이 | 목표 임상 지표 |
|---------|-----|---------|-------------|
| 1 | 무치료 | ΔF508/ΔF508 | 기준선 |
| 2 | Ivacaftor (Kalydeco) | G551D | ppFEV1 +10.6 pp |
| 3 | LUM/IVA (Orkambi) | ΔF508/ΔF508 | ppFEV1 +2.6 pp |
| 4 | TEZ/IVA (Symdeko) | ΔF508/ΔF508 | ppFEV1 +3.4 pp |
| 5 | ETI (Trikafta) | ΔF508/ΔF508 | ppFEV1 +**14.3 pp** |
| 6 | ETI + Tobramycin | ΔF508/ΔF508 | 감염 병용 제어 |
| 7 | Early ETI (age 6) | ΔF508/ΔF508 | 조기 개입 효과 |

### 3. Shiny 앱 (`cf_shiny_app.R`)

6개 탭 구성:
1. **Patient Profile** — 돌연변이 클래스 선택, 임상시험 벤치마크
2. **CFTR Modulator PK** — Ivacaftor/ELX/TEZ 혈장 농도 프로파일
3. **CFTR Function** — 교정율(correction), 증강율(potentiation), 땀 염소
4. **Lung Function** — ppFEV1 궤적, 누적 악화 횟수
5. **Scenario Comparison** — 5가지 치료 비교, 52주 엔드포인트 테이블
6. **ASL & Infection** — 기도 표면 액체 높이, 세균 부담, 염증 마커

---

## 약물 PK/PD 파라미터 요약

| 약물 | 경로 | F (%) | Vc (L) | t½ (h) | 주요 메커니즘 |
|-----|------|--------|--------|--------|------------|
| **Ivacaftor** (VX-770) | p.o. 150mg q12h | 67 | 97 | 12 | CFTR 전위상태 안정화 (↑Po) |
| **Elexacaftor** (VX-445) | p.o. 200mg q24h | 80 | 193 | 27 | NBD1 결합 → CFTR 교정 |
| **Tezacaftor** (VX-661) | p.o. 100mg q24h | 70 | 271 | 14 | MSD1 결합 → 이동 교정 |
| **Lumacaftor** (VX-809) | p.o. 200mg q12h | 65 | ~300 | 26 | MSD1 결합 (단독 효과 제한) |
| **Tobramycin** (TOBI) | 흡입 300mg bid | — | — | — | Pa 살균 (aminoglycoside) |

---

## 임상시험 보정 기준 (Calibration Targets)

| 임상시험 | 치료 | 주요 지표 | 결과 |
|---------|-----|---------|-----|
| STRIVE (NEJM 2011) | Ivacaftor | ppFEV1 | +10.6 pp (G551D) |
| TRAFFIC/TRANSPORT (NEJM 2015) | LUM/IVA | ppFEV1 | +2.6–3.0 pp (ΔF508) |
| EVOLENT (NEJM 2017) | TEZ/IVA | ppFEV1 | +3.4 pp (ΔF508) |
| VX-445-102 (NEJM 2019) | **ETI (Trikafta)** | ppFEV1 | **+14.3 pp** |
| VX-445-103 (AURORA) | **ETI** | Sweat Cl | **-41.8 mmol/L** |

---

## 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [cf_qsp_model.dot](cf_qsp_model.dot) | Graphviz 기계론적 지도 (165 노드, 12 클러스터) |
| [cf_qsp_model.svg](cf_qsp_model.svg) | SVG 벡터 이미지 |
| [cf_qsp_model.png](cf_qsp_model.png) | PNG 150 dpi |
| [cf_mrgsolve_model.R](cf_mrgsolve_model.R) | mrgsolve ODE 모델 (25 구획, 7 시나리오) |
| [cf_shiny_app.R](cf_shiny_app.R) | Shiny 대시보드 (6탭) |
| [cf_references.md](cf_references.md) | 참고문헌 (53편, PubMed 링크) |

---

## 미리보기

[![CF QSP Model](cf_qsp_model.png)](cf_qsp_model.svg)
