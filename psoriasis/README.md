# 건선 (Psoriasis) — QSP Model

> **분류:** 만성 자가면역 피부질환 / Chronic Autoimmune Skin Disease  
> **핵심 경로:** IL-23/IL-17 axis · TNF-α · Keratinocyte hyperproliferation  
> **추가일:** 2026-06-20

---

## 개요 (Overview)

건선(Psoriasis vulgaris)은 전 세계 인구의 약 2-3%에서 발생하는 만성 면역 매개 피부 질환으로, 
각질화 세포의 과다 증식과 면역세포 침윤으로 특징적인 인설성 홍반 판(plaque)을 형성합니다.

**핵심 발병 기전:**
- 수지상세포(mDC) → **IL-23** 분비 → **Th17** 분화
- Th17 → **IL-17A** (주요 효과기 사이토카인) 분비
- IL-17A → 각질형성세포 **NF-κB** 활성화 → 과다증식 (acanthosis)
- TNF-α → 혈관 내피세포 활성화 → 면역세포 피부 유입 증폭

---

## 기계론적 지도 (Mechanistic Map)

[![Psoriasis QSP Mechanistic Map](pso_qsp_model.png)](pso_qsp_model.svg)

*클릭하면 SVG 벡터 이미지로 확대 가능*

### 주요 클러스터 구성
| 클러스터 | 구성 요소 |
|----------|-----------|
| 유전·환경 유발인자 | HLA-Cw6, IL23R, CARD14, 스트렙토콕커스, 물리적 외상, 약물, 비만 |
| 선천 면역 활성화 | pDC, mDC, 랑게르한스세포, 대식세포, 호중구, LL-37–DNA 복합체, TLR7/9 |
| 적응 면역(T세포) | Th17(★), Th1, Th22, Treg, ILC3, 조직상주기억T세포 |
| 핵심 사이토카인 | IL-17A(★), IL-22, TNF-α(★), IFN-γ, IL-36, TSLP |
| 각질형성세포 병리 | NF-κB/STAT3 활성화, 과증식, 이상분화, 장벽기능 손상 |
| 혈관·기질 재형성 | VEGF-A 유도 혈관신생, ICAM-1/VCAM-1, 백혈구 유출 |
| 임상 지표 | PASI, BSA, PGA, DLQI, 건선관절염 |
| 생물학적 제제 PK | 아달리무맙·세쿠키누맙·리산키주맙·우스테키누맙 |
| 소분자약물 PK | 아프레밀라스트·토파시티닙·메토트렉세이트 |

---

## mrgsolve ODE 모델 (Pharmacokinetic/Pharmacodynamic Model)

**파일:** [`pso_mrgsolve_model.R`](pso_mrgsolve_model.R)

### 모델 구조 (25개 ODE 구획)

| 구획 범주 | 구획 수 | 주요 구획 |
|-----------|---------|-----------|
| 면역·세포 | 8 | DC, IL-23, Th17, IL-17A, TNF-α, IFN-γ, KC지수, PASI |
| 아달리무맙 PK | 4 | SC저장소, 중심, 말초, TNF-결합 복합체 |
| 세쿠키누맙 PK | 4 | SC저장소, 중심, 말초, IL-17A-결합 |
| 리산키주맙 PK | 3 | SC저장소, 중심, IL-23-결합 |
| 우스테키누맙 PK | 3 | SC저장소, 중심, p40-결합 |
| 아프레밀라스트 PK | 3 | GI, 중심, 말초 |
| 토파시티닙 PK | 2 | GI, 중심 |
| 메토트렉세이트 PK | 3 | GI, 중심, 폴리글루타메이트 |

### 시뮬레이션 시나리오 (7개)

| 시나리오 | 약물 | 용법 | 기전 |
|----------|------|------|------|
| 1 | 무치료 | — | 자연경과 |
| 2 | 아달리무맙 | 40mg SC q2w | TNF-α 중화 |
| 3 | 세쿠키누맙 | 300mg SC wk0-4, q4w | IL-17A 중화 |
| 4 | 리산키주맙 | 150mg SC wk0,4, q12w | IL-23 p19 억제 |
| 5 | 아프레밀라스트 | 30mg PO BID | PDE4 억제 → cAMP↑ |
| 6 | 토파시티닙 | 10mg PO BID | JAK1/3 억제 → STAT3/1↓ |
| 7 | 메토트렉세이트 | 20mg PO 주1회 | DHFR 억제 → 폴레이트↓ |

### 주요 PK 파라미터

| 약물 | 생체이용률(F) | CL (L/h) | V₁ (L) | t½ (일) | Kd / IC50 |
|------|-------------|---------|-------|---------|-----------|
| 아달리무맙 | 0.64 | 0.247 | 7.0 | ~14d | Kd=0.1 nM (TNF-α) |
| 세쿠키누맙 | 0.73 | 0.191 | 7.1 | ~27d | Kd=0.08 nM (IL-17A) |
| 리산키주맙 | 0.89 | 0.078 | 11.2 | ~28d | Kd=0.06 nM (IL-23p19) |
| 우스테키누맙 | 0.57 | 0.252 | 15.1 | ~21d | Kd=0.9 nM (p40) |
| 아프레밀라스트 | 0.73 | 9.5 | 86.6 | 9h | IC50=74 nM (PDE4) |
| 토파시티닙 | 0.74 | 22.8 | 87.0 | 3h | IC50=1-5 nM (JAK1/3) |
| 메토트렉세이트 | 0.70 | 4.8 | 24.0 | 3-10h | IC50=1 nM (DHFR) |

### 임상 보정 데이터 (PASI 반응률 Wk16 기준)

| 약물 | PASI75 (실제) | PASI90 (실제) | 임상시험 |
|------|-------------|-------------|---------|
| 아달리무맙 40mg q2w | 71% | 45% | CHAMPION |
| 세쿠키누맙 300mg q4w | 77-80% | 59-67% | FIXTURE |
| 리산키주맙 150mg q12w | 88-91% | 72-75% | UltIMMa |
| 우스테키누맙 45mg q12w | 67-71% | 42% | PHOENIX-1 |
| 아프레밀라스트 30mg BID | 33-40% | 18% | ESTEEM-1 |
| 토파시티닙 10mg BID | 39-46% | 22% | OPT Pivotal |
| 메토트렉세이트 20mg qw | 26-36% | — | Heydendael 2003 |

---

## Shiny 인터랙티브 앱 (Interactive Dashboard)

**파일:** [`pso_shiny_app.R`](pso_shiny_app.R)

### 6개 탭 구성

| 탭 | 내용 |
|----|------|
| **환자 프로파일** | 기저 PASI, IL-17A, TNF-α, Th17 세포 수 설정; 질환 표현형 선택 |
| **PK 프로파일** | 생물학적제제 및 소분자약물 혈중 농도 시간 곡선 |
| **사이토카인 PD** | IL-17A, TNF-α, Th17, IL-23 동역학 |
| **PASI 엔드포인트** | PASI 점수 변화, PASI75/90/100 반응률 |
| **시나리오 비교** | 7가지 치료법 동시 비교 |
| **바이오마커 대시보드** | 반응률 ValueBox, 히트맵, 바이오마커 시계열 |

---

## 참고문헌 (References)

**파일:** [`pso_references.md`](pso_references.md)

- 총 42개 문헌 (PubMed 링크 포함)
- 섹션 분류: 발병기전 · IL-23/Th17 축 · IL-17A/각질형성세포 · TNF-α · 임상시험 · PK/PD 모델 · 바이오마커 · TYK2 억제제(신규)

---

## 파일 목록 (File List)

| 파일 | 설명 |
|------|------|
| [`pso_qsp_model.dot`](pso_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 11 클러스터) |
| [`pso_qsp_model.svg`](pso_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [`pso_qsp_model.png`](pso_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [`pso_mrgsolve_model.R`](pso_mrgsolve_model.R) | mrgsolve ODE 모델 (25 구획, 7 시나리오) |
| [`pso_shiny_app.R`](pso_shiny_app.R) | Shiny 대시보드 (6탭) |
| [`pso_references.md`](pso_references.md) | 참고문헌 42편 (PubMed 링크) |

---

## 주요 임상 맥락 (Clinical Context)

- **유병률:** 전 세계 2-3%, 한국 약 1.5%
- **중등도-중증 기준:** PASI ≥ 10 또는 BSA ≥ 10% 또는 DLQI ≥ 10
- **치료 목표:** PASI90 (피부 90% 개선) 또는 IGA 0/1 (완전/거의 완전 소실)
- **건선관절염 동반:** 약 30%에서 발생; IL-17A, TNF-α 억제제 병용 효과적
- **동반질환:** 심혈관질환(CVD), 대사증후군, 염증성장질환(IBD), 우울증 위험 증가

---

## 신규 치료 표적 (Emerging Targets)

| 표적 | 약물 | 개발 단계 |
|------|------|-----------|
| TYK2 (IL-23/IFN-γ 경로) | Deucravacitinib | FDA 승인 (2022) |
| IL-17A + IL-17F | Bimekizumab | FDA 승인 (2023) |
| IL-13 | Tralokinumab | 탐색 임상 |
| OX40L | Amlitelimab | Phase 3 진행 중 |
| IL-4Rα | Dupilumab | 건선 적응증 탐색 |
