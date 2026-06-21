# 만성 신우신염 (Chronic Pyelonephritis) — QSP Model

> **분류**: 만성질환 / 신장-요로계  
> **날짜**: 2026-06-20  
> **약어**: CPN

---

## 개요 (Overview)

만성 신우신염은 **반복적인 상부 요로 세균 감염** 또는 **방광요관역류(VUR)**로 인해 신장 실질 및 신배(renal calyx)에 영구적인 반흔(scar)이 형성되는 만성 신장 질환입니다. 최종적으로 만성 신부전(CKD)으로 진행할 수 있습니다.

| 항목 | 내용 |
|------|------|
| 주요 원인균 | *E. coli* (UPEC, 75%), *Klebsiella*, *Proteus*, *Pseudomonas* |
| 발병 기전 | VUR → 세균 역류 → TLR4/NFκB 활성화 → 염증 → TGF-β1 → 섬유화 → 신장 반흔 |
| 주요 위험인자 | VUR, 당뇨병, 면역억제, 임신, 요로폐색, 신경인성방광 |
| 임상 양상 | 고열, 옆구리 통증, 농뇨/혈뇨, 만성 고혈압, 단백뇨, GFR 저하 |
| 재발률 | 5년 재발 ~50%(위험인자 있는 경우) |

---

## 주요 병태생리 경로 (Key Pathophysiological Pathways)

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|----------|
| UPEC 부착 & 침입 | Type I Fimbriae(FimH), P-Fimbriae(PapG), α-Hemolysin | 신우 집락화 → 신실질 침입 |
| VUR 매개 역류 | Grade I–V VUR, 복합 유두 | 신내 역류 → 반복 감염 |
| TLR4/NFκB 활성화 | LPS → TLR4 → MyD88 → IKK → NFκB | IL-1β, IL-6, IL-8, TNF-α 과다 생성 |
| NLRP3 인플라마좀 | HlyA → K⁺ efflux → NLRP3 → Caspase-1 → IL-1β, Pyroptosis | 세뇨관 상피 괴사 |
| 보체 활성화 | LPS → C3 → C5a, MAC | 세포 용해, 호중구 동원 |
| EMT & TGF-β1 섬유화 | TGF-β1 → EMT → 근섬유모세포 → Collagen I/III | 신장 간질 섬유화, DMSA 흡수 감소 |
| 레닌-안지오텐신 | Scar → RAAS 활성화 → Ang II → TGF-β1 | 이차성 고혈압, GFR 저하 |

---

## 항생제 PK/PD 파라미터

| 항생제 | 용법 | 생체이용률 | t₁/₂ | PK/PD 지표 | 임상 목표 |
|-------|------|----------|------|------------|----------|
| Ciprofloxacin | 500 mg PO BID | 70% | ~4 h | fAUC/MIC | > 125 |
| TMP-SMX | 160/800 mg PO BID | 95% / 85% | ~10 h | T>MIC | > 40% |
| Nitrofurantoin (예방) | 100 mg PO QD | 75% | 0.3–1 h | 요중 Cmax/MIC | > 4× |
| Fosfomycin | 3 g PO 단회 | ~40%(요중) | 4–8 h | Cmax/MIC (요중) | > 8× |
| Gentamicin (중증) | 5 mg/kg IV QD | — (IV) | ~2 h | Cmax/MIC | > 10 |

---

## 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`cpn_qsp_model.dot`](cpn_qsp_model.dot) | Graphviz 기계론적 지도 소스 (10 클러스터, 140+ 노드) |
| [`cpn_qsp_model.svg`](cpn_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [`cpn_qsp_model.png`](cpn_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [`cpn_mrgsolve_model.R`](cpn_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획, 7 시나리오) |
| [`cpn_shiny_app.R`](cpn_shiny_app.R) | Shiny 대시보드 (7탭: 환자/PK/세균/신기능/시나리오/바이오마커/참고문헌) |
| [`cpn_references.md`](cpn_references.md) | 참고문헌 38편 (PubMed 링크 포함) |

---

## mrgsolve ODE 구획 구조 (18 compartments)

```
PK  : Cipro_gut → Cipro_C ↔ Cipro_P → [renal conc.]
      TMP_gut   → TMP_C
      NIT_gut   → NIT_urine

Disease :
  Bacteria   (log₁₀ CFU/g)       — 세균 부하 (성장/사멸/면역)
  Biofilm    (0–1)                — 생물막 분율
  Neutrophil (normalised)        — 호중구 (급성 염증)
  Macrophage (normalised)        — 대식세포 (만성 염증)
  IL6        (normalised)        — 사이토카인 (급성기)
  TGFb1      (normalised)        — 섬유화 구동인자
  Collagen   (normalised)        — 간질 콜라겐
  RenalScar  (0–1)               — 신장 반흔 (비가역)
  GFR        (mL/min/1.73m²)     — 신사구체 여과율
```

---

## 치료 시나리오 (7가지)

| # | 시나리오 | 결과 요약 |
|---|---------|----------|
| S1 | 미치료 (No Antibiotics) | 세균 부하 지속 → 빠른 GFR 저하, 반흔 형성 |
| S2 | Ciprofloxacin 500 mg BID × 14일 | 7–10일 내 세균 청소, GFR 보존 |
| S3 | TMP-SMX 160/800 mg BID × 14일 | 감수성 균주에서 S2와 유사한 효과 |
| S4 | Ciprofloxacin 500 mg BID × 7일 | 14일 대비 효과 약간 감소, 내성 선택압 낮음 |
| S5 | Nitrofurantoin 100 mg QD × 6개월 예방 | 재발률 60% 감소, GFR 보존 효과 |
| S6 | Cipro 14일 → Nitrofurantoin 6개월 | 가장 우수한 GFR 보존 및 반흔 억제 |
| S7 | TMP-SMX + 내성균 (MIC × 4) | 치료 실패 시뮬레이션 → 세균 지속, 빠른 섬유화 |

---

## 기계론적 지도 미리보기

[![CPN QSP Mechanistic Map](cpn_qsp_model.png)](cpn_qsp_model.svg)

*클릭하면 확대 가능한 SVG로 이동합니다.*

---

## Shiny 앱 탭 구성 (7탭)

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | GFR, VUR 등급, 동반질환, 항생제 선택 |
| ② PK 모니터링 | Cipro/TMP/NIT 혈중 농도, fAUC/MIC |
| ③ 세균 동태 | log₁₀ CFU, 생물막, Kill Rate, 면역세포 |
| ④ 신기능 (GFR) | GFR 추이, 크레아티닌, 반흔, CKD 단계 전환 |
| ⑤ 시나리오 비교 | 5가지 항생제 전략 동시 비교 |
| ⑥ 바이오마커 | IL-6, TGF-β1, 콜라겐, 요로패혈증 위험 |
| ⑦ 참고문헌 | 섹션별 분류 PubMed 링크 |
