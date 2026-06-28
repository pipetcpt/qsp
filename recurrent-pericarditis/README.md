# Recurrent Pericarditis (RP) — QSP Model

[![RP Mechanistic Map](rp_qsp_model.png)](rp_qsp_model.svg)

**분류 (Classification):** 심혈관 / 자가염증 (Cardiovascular · Autoinflammatory)  
**디렉토리 (Directory):** `recurrent-pericarditis/`

---

## 개요 (Overview)

반복 심낭염(Recurrent Pericarditis, RP)은 급성 심낭염 초기 에피소드 이후 증상이 없는 기간(≥4-6주)을 거쳐 재발하는 염증성 심낭 질환으로, 치료하지 않으면 **30-50%**의 환자에서 재발합니다.  
핵심 병리기전은 **NLRP3 인플라마솜 → IL-1β/IL-18 축**이며, 콜히친·NSAIDs·코르티코스테로이드·아나킨라(IL-1Ra)·리로나셉트(IL-1 trap)가 주요 치료제입니다.

Recurrent Pericarditis affects ~30-50% of patients after a first episode without prophylaxis. The central driver is the **NLRP3 inflammasome → IL-1β/IL-18 axis**, making it an excellent target for IL-1 blockers. This QSP model integrates all major pharmacological pathways to simulate disease dynamics and treatment response.

---

## 파일 목록 (File List)

| 파일 | 설명 |
|------|------|
| `rp_qsp_model.dot` | Graphviz 기계론적 지도 (14 클러스터, 162+ 노드) |
| `rp_qsp_model.svg` | 벡터 그래픽 (인터랙티브 브라우저 뷰) |
| `rp_qsp_model.png` | PNG 섬네일 (150 dpi) |
| `rp_mrgsolve_model.R` | mrgsolve ODE PK/PD 모델 (24구획, 7 시나리오) |
| `rp_shiny_app.R` | Shiny 대시보드 (7탭) |
| `rp_references.md` | 참고문헌 57편 |

---

## 병태생리 (Pathophysiology)

| 핵심 병리 경로 | 세부 메커니즘 |
|--------------|-------------|
| **NLRP3 인플라마솜** | Signal 1 (TLR→NF-κB→NLRP3/pro-IL-1β 전사) + Signal 2 (ATP, 결정체, K⁺ efflux, ROS) → ASC 스펙 → Casp-1 → IL-1β/IL-18 절단·분비 |
| **IL-1β 증폭 루프** | IL-1β → IL-1R1 → IRAK4 → NF-κB → TNF/IL-6/CXCL8 → 추가 NLRP3 활성화 (양성 피드백) |
| **중성구 침윤** | CXCL8/LTB4 → 심낭 중성구 유입 → NET 형성 → DAMP 방출 → 인플라마솜 재활성화 |
| **섬유화 경로** | M2 대식세포 → TGF-β1 → 심낭 근섬유모세포 → Collagen I/III → 심낭 비후 → 수축성 심낭염 |
| **적응면역** | MEMORY T cells → 재활성화 (재발 유발); B cell → 항심장 항체 → 면역복합체 → 보체 활성화 |
| **에이코사노이드** | AA → COX-2 → PGE2 (발열, 통증, 혈관투과성 증가) |

---

## 모델 사양 (Model Specifications)

### 기계론적 지도 (Mechanistic Map)
- **14 서브그래프 클러스터**: 병인/Etiology · DAMP-PAMP 인식 · NLRP3 인플라마솜 · NF-κB · 선천면역세포 · 적응면역 · 에이코사노이드 · 심낭 병리 · 임상 바이오마커 · 콜히친 PK/PD · NSAID PK/PD · 코르티코스테로이드 PK/PD · 생물학적 제제 PK/PD · 임상 결과
- **162+ 노드**: 모든 주요 분자·세포·임상 구성 요소 포함

### mrgsolve ODE 모델 (24 구획)

| 구획 그룹 | 구획 수 | 내용 |
|-----------|---------|------|
| Colchicine PK | 3 | 소화관·중심·말초 |
| Ibuprofen PK | 2 | 소화관·중심 |
| Prednisone PK | 2 | 소화관·중심 |
| Anakinra PK | 2 | SC depot·중심 |
| Rilonacept PK | 2 | SC depot·중심 |
| 면역·염증 | 7 | NLRP3_ACT·IL1B·IL18·TNF·IL6·NEUTRO·M1_MACRO |
| 심낭 병리 | 4 | INFLAM·EFFUSION·FIBRIN·FIBROSIS |
| 임상 바이오마커 | 2 | CRP·PAIN |

### 치료 시나리오 (7가지)

| # | 시나리오 | 근거 임상시험 |
|---|---------|------------|
| 1 | 무치료 (자연경과) | — |
| 2 | Ibuprofen 600 mg TID × 4주 | ESC 2015 |
| 3 | Colchicine 0.5 mg BID × 3개월 | COPE, ICAP |
| 4 | Colchicine + Ibuprofen (병합) | COPE, ICAP |
| 5 | Prednisone 0.5 mg/kg/d → 점진 감량 | ESC 2015 |
| 6 | Anakinra 100 mg/d SC × 6개월 | AIRTRIP 2016 |
| 7 | Rilonacept 320→160 mg SC qw | RHAPSODY 2021 |

### 임상시험 보정 데이터 (Calibration)

| 임상시험 | 약물 | 대조군 재발률 | 치료군 재발률 | RRR |
|---------|------|------------|------------|-----|
| COPE (2005) | Colchicine + ASA | 45% | 24% | 47% |
| ICAP (2013) | Colchicine 0.5 mg BID | 32.3% | 16.7% | 48% |
| CORP (2011) | Colchicine (2nd 에피소드) | 45.5% | 19.2% | 58% |
| AIRTRIP (2016) | Anakinra 100 mg/d | 90.9% | 18.2% | 80% |
| RHAPSODY (2021) | Rilonacept 320→160 mg qw | 74.4% | 8.8% | HR 0.04 |

---

## Shiny 앱 탭 구성 (7탭)

| 탭 | 내용 |
|----|------|
| ① 환자 & 시나리오 | 환자 프로파일, 치료 선택, ESC 진단기준, 개요 플롯 |
| ② Drug PK | 농도-시간 곡선, PD 억제율, PK 파라미터 표 |
| ③ 인플라마솜/사이토카인 | NLRP3·IL-1β·IL-18·TNF·IL-6·면역세포 dynamics |
| ④ 심낭 병리 | 염증·삼출액·피브린·섬유화·위험 계층화 표 |
| ⑤ 임상 엔드포인트 | 통증 VAS·CRP·요약 표·치료 목표 |
| ⑥ 시나리오 비교 | 7개 시나리오 염증/CRP/IL-1β/삼출액 비교 + 임상시험 벤치마크 |
| ⑦ 바이오마커 | 시점별 스냅샷, 궤적, Emax 곡선, 위험 레이더 차트 |

---

## 핵심 약물 작용 기전 요약

```
Colchicine:   β-tubulin 결합 → 미세소관 분해
              → NLRP3 ASC 조립 억제 (IC50 ≈ 0.5 ng/mL)
              → 중성구 이동 억제, L-selectin shed 억제

Ibuprofen:    COX-1/2 억제 → PGH2 감소 → PGE2↓
              → 발열·통증·혈관투과성 감소

Prednisone:   GRα → NF-κB 전사억제 → IL-1β/TNF/IL-6↓
              → 단, 급격한 테이퍼 시 반동 재발 위험

Anakinra:     IL-1R1 경쟁적 차단 (IC50 ≈ 0.1 nM)
              → IL-1β 신호 하부 차단

Rilonacept:   IL-1α/β 이중 포획 (KD < 1 pM)
              → 순환 IL-1β 중화 → NLRP3 양성 루프 차단
```

---

## 참고문헌 (References)

총 57편 — [rp_references.md](rp_references.md) 참조

주요 참고:
- Imazio M, et al. NEJM 2016 (AIRTRIP): https://pubmed.ncbi.nlm.nih.gov/27668557/
- Klein AL, et al. NEJM 2021 (RHAPSODY): https://pubmed.ncbi.nlm.nih.gov/33405895/
- Imazio M, et al. NEJM 2013 (ICAP): https://pubmed.ncbi.nlm.nih.gov/24131175/
- Adler Y, et al. Eur Heart J 2015 (ESC Guidelines): https://pubmed.ncbi.nlm.nih.gov/26320112/

---

## 면책 조항 (Disclaimer)

본 모델은 교육 및 연구 목적의 QSP 모델이며, 임상 의사결정에 직접 사용해서는 안 됩니다.
