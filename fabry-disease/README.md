# 파브리병 QSP 모델 (Fabry Disease Quantitative Systems Pharmacology Model)

> **디렉토리:** `fabry-disease/` | **약어:** FBR | **날짜:** 2026-06-24  
> **분류:** 희귀·리소소말 저장 질환 (X-연관 리소소말 저장 질환)

[![FBR QSP 기계론적 지도](fbr_qsp_model.png)](fbr_qsp_model.svg)

---

## 질환 개요

**파브리병(Fabry disease)** 은 *GLA* 유전자(Xq22.1) 변이로 인한 X-연관 리소소말 저장 질환으로, α-갈락토시다제 A(α-Gal A) 효소의 결핍 또는 기능 저하로 인해 당스핑고지질인 **글로보트리아오실세라미드(Gb3/GL-3)** 와 그 탈아실화 유도체인 **lyso-Gb3** 가 신장, 심장, 뇌, 피부, 말초신경계 등 전신에 점진적으로 축적됩니다.

| 특성 | 내용 |
|------|------|
| 유병률 | 1:40,000–1:117,000 (고전형) / 1:3,000–1:10,000 (후기 발현형 포함) |
| 유전 방식 | X-연관 (반성 유전) — 반접합 남성 중증, 이형접합 여성 가변 |
| 주요 표현형 | 고전형(잔여 효소 <1%), 후기 발현형 심장/신장형(잔여 1–30%) |
| 진단 지연 | 평균 10–20년 (증상 발현 → 진단) |
| 사망 원인 | 신부전, 심장 사건(급사·심부전), 뇌졸중 |

---

## 핵심 기전 (14개 클러스터)

| 클러스터 | 핵심 내용 |
|---------|----------|
| 1. 유전적 기반 | GLA 유전자 변이(Xq22.1): 미스센스 ~60%, 넌센스/프레임시프트, 스플라이싱; 적합 변이 (~40%, 미갈라스타트 대상) |
| 2. α-Gal A 효소 생물학 | ER 합성 → 골지체 M6P 인산화 → M6P 수용체 → 리소솜 전달 → pH 4.5–5.0에서 Gb3 가수분해 |
| 3. 당스핑고지질 대사 | 세라미드 → GCS → GlcCer → LacCer → A4GalT → Gb3 생합성; Gb3 탈아실화 → lyso-Gb3(독성 신호분자) |
| 4. 신장 병리 | 족세포 Gb3 축적 → 족돌기 소실 → 단백뇨(UPCR↑) → TGF-β 섬유화 → FSGS → eGFR 저하 → ESRD |
| 5. 심장 병리 | 심근세포 Gb3 → 좌심실 비후(LVMi↑) → 이완기 기능 장애 → 심근 섬유화(LGE) → 부정맥 → 급사 |
| 6. 신경계 병리 | DRG Gb3 → 소섬유 신경병증 → 신경병성 통증(BPI-SF) + CNS 혈관내피 → 백질 병변 → TIA/뇌졸중 |
| 7. 기타 장기 | 혈관각화종, 무한증, 각막 소용돌이, 위장관 운동장애, 감각신경성 난청 |
| 8. 염증 폭포 | Lyso-Gb3 → TLR4/NF-κB → IL-6/TNF-α → NLRP3 인플라마좀 → eNOS↓ → 내피세포 활성화 |
| 9. ERT PK/PD | 아갈시다제 베타(1 mg/kg Q2W, t½ ~45min), 알파(0.2 mg/kg Q2W), 페구니알시다제 알파(1 mg/kg Q4W, t½ ~80h) → M6P → 리소솜 Gb3 분해(Emax ~80%) |
| 10. 샤페론 (미갈라스타트) | 150 mg QOD 경구; α-Gal A 잘못 접힌 단백질 안정화; 적합 변이 전용; ATTRACT 임상 (vs ERT 비열등) |
| 11. 기질감소요법 (SRT) | 루세라스탓 1000 mg TID (GCS 억제, Gb3 전구체↓); MODIFY 임상 (BPI-SF −1.5점); 벵글루스탓 (CNS 투과 가능) |
| 12. 바이오마커 | 혈장 lyso-Gb3 (μg/L, 가장 민감), 소변 Gb3 (nmol/mg Cr), α-Gal A 활성 (nmol/h/mg), DBS 신생아 선별 |
| 13. 임상 엔드포인트 | eGFR 기울기, UPCR, LVMi, BPI-SF 통증, EQ-5D QoL, MSSI 중증도 |
| 14. 자연 경과 | 고전형 남성(소아기 발현) vs 후기 발현형 vs 여성 보인자; Fabry Registry 데이터; 진단 지연 평균 10–20년 |

---

## mrgsolve ODE 모델 (22 구획)

| 모듈 | 구획 | 핵심 동역학 |
|------|------|------------|
| 아갈시다제 베타 PK | A_AGAB_C, A_AGAB_P, A_AGAB_LYS | 2구획 + 리소솜 전달; CL=0.42 L/h; t½ ~45min |
| 아갈시다제 알파 PK | A_AGAA_C, A_AGAA_P, A_AGAA_LYS | 2구획; CL=0.55 L/h |
| 미갈라스타트 PK | A_MIG_GUT, A_MIG_C | 1구획 경구; ka=0.82/h, F=75%, t½~3.5h |
| 루세라스탓 PK | A_LUC_GUT, A_LUC_C | 1구획 경구; IC50_GCS=0.18 μg/mL, Emax=42% |
| α-Gal A 효소 | E_GalA | ERT(Emax 70 nmol/h/mg) + 미갈라스타트(Emax 6) + 기저 잔여 활성 |
| 당스핑고지질 | GB3_PLM, GB3_KID, GB3_HRT, LGB3_PLM | 합성 – 효소 의존 분해; SRT 상류 감소 |
| 염증 | INFLAM | lyso-Gb3 구동 k_in – 제거 k_out ODE |
| 신장 기능 | eGFR, UPCR | Gb3_KID 누적 의존성 감소; ERT 보호 |
| 심장 기능 | LVMi | Gb3_HRT 구동 비후; ERT → 역전 |
| 신경병성 통증 | PAIN | BPI-SF; 효소 활성 의존 감소 |

---

## 6가지 치료 시나리오 임상 근거

| 시나리오 | 치료법 | 임상시험 | 주요 결과 |
|---------|--------|---------|---------|
| S1: 자연경과 | 없음 | Mehta 2009 Eur J Clin Invest | eGFR −3–12/yr, LVMi 매년 증가, lyso-Gb3 30–80 μg/L |
| S2: 아갈시다제 베타 | 1 mg/kg IV Q2W | FABRY-001 (Eng 2001 NEJM); Banikazemi 2007 AIM | 신장·심장·뇌 복합 사건 61% 감소 |
| S3: 아갈시다제 알파 | 0.2 mg/kg IV Q2W | Schiffmann 2001 Ann Intern Med | 신경병성 통증 개선, 신기능 안정 |
| S4: 미갈라스타트 | 150 mg PO QOD | ATTRACT (Germain 2016 NEJM); Hughes 2017 Lancet | ERT와 비열등; eGFR 기울기 −0.3 ml/min/yr |
| S5: 페구니알시다제 알파 | 1 mg/kg IV Q4W | BRIGHT (Schiffmann 2021 JAMA Intern Med) | eGFR 안정, lyso-Gb3 −50%, Q4W 편의성 |
| S6: ERT + 루세라스탓 | 아갈시다제 베타 + 1000 mg TID | MODIFY (Lenders 2022 Lancet DE) | Gb3 추가 감소, BPI-SF −1.5점 통증 개선 |

---

## QSP 모델 구성 파일

| 산출물 | 파일 | 사양 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`fbr_qsp_model.dot`](fbr_qsp_model.dot) · [`fbr_qsp_model.svg`](fbr_qsp_model.svg) · [`fbr_qsp_model.png`](fbr_qsp_model.png) | **138 노드, 14 클러스터** |
| ⚙️ mrgsolve ODE | [`fbr_mrgsolve_model.R`](fbr_mrgsolve_model.R) | **22 구획 ODE**, **6 치료 시나리오** |
| 📊 Shiny 앱 | [`fbr_shiny_app.R`](fbr_shiny_app.R) | **8 탭** (환자 프로파일·PK/효소·Gb3 동역학·신장·심장·시나리오 비교·바이오마커·가상 집단) |
| 📚 참고문헌 | [`fbr_references.md`](fbr_references.md) | **60개 PubMed 인용** (14개 섹션) |

---

## 실행 방법

```bash
# 1. 기계론적 지도 렌더링 (Graphviz)
dot -Tsvg fbr_qsp_model.dot -o fbr_qsp_model.svg
dot -Tpng -Gdpi=150 fbr_qsp_model.dot -o fbr_qsp_model.png
```

```r
# 2. mrgsolve ODE 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))
source("fbr_mrgsolve_model.R")

# 3. Shiny 대시보드 실행
shiny::runApp("fbr_shiny_app.R")
```

---

## 주요 파라미터 요약

| 약물 | 용법 | t½ | 작용기전 | 임상 효능 |
|------|------|-----|---------|----------|
| 아갈시다제 베타 (Fabrazyme) | 1 mg/kg IV Q2W | ~45분 | M6P→리소솜 Gb3 분해 (Emax ~80%) | 복합 사건 61% 감소 |
| 아갈시다제 알파 (Replagal) | 0.2 mg/kg IV Q2W | ~45–110분 | M6P→리소솜 Gb3 분해 (Emax ~70%) | 신경병성 통증·신기능 개선 |
| 페구니알시다제 알파 (Elfabrio) | 1 mg/kg IV Q4W | ~80시간 | PEGylation 연장 t½ (Emax ~85%) | eGFR 안정, Q4W 편의성 |
| 미갈라스타트 (Galafold) | 150 mg PO QOD | ~3.5시간 | α-Gal A 샤페론 (EC50 ~0.25 μg/mL) | 적합 변이에서 ERT 비열등 |
| 루세라스탓 (조합) | 1000 mg PO TID | ~8시간 | GCS 억제 IC50 ~0.18 μg/mL (Emax 42%) | Gb3↓, BPI-SF −1.5점 |

---

*Claude Code Routine (CCR) — 자동 생성 QSP 모델 | 2026-06-24*
