# 무증상 고요산혈증 (Asymptomatic Hyperuricemia, AHY)
## QSP Disease Model — 정량적 시스템 약리학 모델

---

## 개요 (Overview)

**무증상 고요산혈증(Asymptomatic Hyperuricemia, AHY)**은 혈청 요산 농도가 남성 ≥7.0mg/dL, 여성 ≥6.0mg/dL를 초과하지만, 통풍 발작·토파이·요산 신석증 등 임상 증상이 없는 상태를 말합니다.

전 세계 성인의 15-25%, 아시아인에서는 더 높은 유병률을 보이며, 장기적으로 통풍, 만성 신부전(CKD), 고혈압, 심혈관 질환의 위험인자로 작용합니다.

---

## 기계론적 지도 (Mechanistic Map)

[![AHY QSP Model](ahy_qsp_model.png)](ahy_qsp_model.svg)

> 클릭하면 전체 SVG 벡터 이미지가 열립니다 (100+ 노드, 15 서브그래프 클러스터).

---

## 병태생리 요약 (Pathophysiology Summary)

### 요산 대사 핵심 경로

| 경로 | 주요 효소/수송체 | 임상 의미 |
|------|-----------------|----------|
| 퓨린 생합성 (De Novo) | PRPP 합성효소, 복합 효소군 | 과생성 → SUA↑ |
| 퓨린 이화 | 잔틴 산화효소(XO/XDH) | 알로퓨리놀·페북소스타트 표적 |
| 신장 재흡수 | URAT1(SLC22A12), GLUT9(SLC2A9) | 요산배설 감소 (배설 감소형) |
| 신장 분비 | OAT1, OAT3, ABCG2, NPT1/4 | 분비 감소 → SUA↑ |
| 장관 분비 | ABCG2(장관) | Q141K 변이 → 분비↓ |
| MSU 결정 | 포화점 6.8mg/dL | 핵화·성장·용해 동적 균형 |

### 고요산혈증의 장기 영향

```
SUA↑ (≥7.0mg/dL)
├─ 직접 독성 → 내피세포 손상 → NO↓ → 혈압↑ → 심혈관 위험↑
├─ XO 부산물 (H₂O₂, O₂⁻) → 산화스트레스 → 죽상경화↑
├─ MSU 결정 형성 (>6.8mg/dL) → NLRP3 → IL-1β → 통풍 발작
├─ 신세뇨관 손상 → 간질 섬유화 → eGFR 감소 → CKD 악화
└─ 인슐린 저항성↑ → 대사증후군 → NAFLD, T2DM
```

---

## 주요 병인 유형 (Subtypes)

| 유형 | 빈도 | 원인 | 특징 |
|------|------|------|------|
| 배설 감소형 | ~85% | URAT1↑, ABCG2↓, GFR↓ | 요중 UA <600mg/day |
| 과생성형 | ~10% | XO 과활성, HGPRT 결핍 | 요중 UA >800mg/day |
| 혼합형 | ~5% | 복합 원인 | 식이+유전 복합 |

---

## 유전적 위험 인자 (Key Genetic Variants)

| 유전자 | 변이 | 효과 | 빈도(아시아) |
|--------|------|------|------------|
| ABCG2 | Q141K (rs2231142) | 장관+신장 분비 50%↓ | 34% |
| SLC22A12 | rs11602903 | URAT1 재흡수 기능 변화 | 5-20% |
| SLC2A9 | rs16890979 | GLUT9 재흡수↑ | 10-25% |
| XDH | 다수 | XO 활성↑ | 드묾 |
| HLA-B*58:01 | — | 알로퓨리놀 독성 위험 | 6-8% (한국/중국) |

---

## 약물 PK/PD 파라미터 (Drug PK/PD Parameters)

### XO 억제제

| 파라미터 | 알로퓨리놀 300mg | 페북소스타트 80mg |
|---------|----------------|-----------------|
| 생체이용률 (F) | 90% | 84% |
| 최고혈중농도 | 2-3 μg/mL | 3.4 μg/mL |
| 반감기 (t½) | 1-2h (옥시퓨리놀: 18-30h) | 5-8h |
| XO IC50 | 옥시퓨리놀: 8.0mg/L | 0.001mg/L (non-purine) |
| SUA 감소 효과 | 30-40% | 40-53% |
| 신장 배설 | 주 경로 | 이중 경로(간/신) |

### 요산배설촉진제

| 파라미터 | 레시누라드 200mg | 도티누라드 4mg |
|---------|----------------|--------------|
| 생체이용률 | ≈100% | ≈100% |
| t½ | 5h | 14-17h |
| 표적 수송체 | URAT1+OAT4 | URAT1 선택적 |
| SUA 감소 | 병용 시 15-30% 추가 | 단독 15-20% |

### 페글로티카제 (생물학적 제제)

| 파라미터 | 값 |
|---------|---|
| 투여 경로 | IV 8mg q2주 |
| t½ | 6-14일 (PEG화) |
| 기전 | UA → 알란토인 (용해도 10배↑) |
| SUA 감소 | >80% (24시간 내) |
| 실패율 | 40-50% (항PEG 항체) |

---

## mrgsolve ODE 모델 구조 (Model Structure)

### 19개 구획 (Compartments)

| # | 구획 | 설명 |
|---|------|------|
| 1 | UA_plasma | 혈청 요산 (mg/dL) |
| 2 | UA_tissue | 조직 요산 풀 (mg) |
| 3 | XO_free | 잔틴 산화효소 활성 (normalized) |
| 4 | Oxypurinol | 옥시퓨리놀 혈중 농도 (mg/L) |
| 5 | Febuxostat_C | 페북소스타트 혈중 농도 (mg/L) |
| 6 | Uricosuric_C | 요산배설촉진제 혈중 농도 (mg/L) |
| 7 | URAT1_free | 활성 URAT1 분획 |
| 8 | UrinaryUA | 요중 요산 (mg/day) |
| 9 | MSU_depot | MSU 결정 침착량 (mg) |
| 10 | Endothelial_fn | 내피 기능 (0-1) |
| 11 | NO_level | 산화질소 수준 (rel.) |
| 12 | BP | 평균 동맥압 (mmHg) |
| 13 | GFR | 사구체 여과율 (mL/min/1.73m²) |
| 14 | IL1beta | IL-1β (pg/mL) |
| 15 | CRP | hs-CRP (mg/L) |
| 16 | InsulinResist | HOMA-IR |
| 17 | CV_risk_score | 누적 심혈관 위험 점수 |
| 18 | Tophus_vol | 토파이 부피 (mm³) |
| 19 | ABCG2_frac | 기능성 ABCG2 분획 |

---

## 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 설명 | 예상 SUA 변화 (2년) |
|---------|------|-------------------|
| 1. 미치료 AHY | SUA=7.5 기저 | →8.2mg/dL (점진적 상승) |
| 2. 알로퓨리놀 300mg | XO 억제 30-40% | →5.2mg/dL (목표 달성) |
| 3. 페북소스타트 80mg | XO 억제 40-53% | →4.5mg/dL (목표 달성) |
| 4. 알로+요산배설촉진제 병용 | 이중 기전 | →4.0mg/dL (적극적 치료) |
| 5. 고과당식이+알코올 | 생활습관 위험 | →9.5mg/dL (급격 상승) |
| 6. ABCG2 Q141K+페북소120 | 유전변이+최대용량 | →4.8mg/dL |
| 7. 적극적 목표 SUA<6 | 알로600mg+생활습관 교정 | →5.0mg/dL |

---

## Shiny 대시보드 탭 구조 (Shiny App Tabs)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 개요, 초기 위험 평가, SUA 시계열 |
| 2. 약물 PK | 혈중 농도-시간 곡선, XO 억제율 |
| 3. 요산 역학 | SUA, 요중 UA, MSU 결정 침착, 통풍 발작 위험 |
| 4. 심혈관·신장 영향 | eGFR, 혈압, 심혈관 위험, 내피기능/NO |
| 5. 시나리오 비교 | 7개 치료 시나리오 동시 비교, 2년 결과표 |
| 6. 바이오마커 | IL-1β, hs-CRP, HOMA-IR, SUA-CRP 산점도 |

---

## 파일 목록 (Model Files)

| 파일 | 설명 |
|------|------|
| [ahy_qsp_model.dot](ahy_qsp_model.dot) | Graphviz 기계론적 지도 소스 (15 클러스터, 100+ 노드) |
| [ahy_qsp_model.svg](ahy_qsp_model.svg) | SVG 벡터 이미지 (고해상도 확대 가능) |
| [ahy_qsp_model.png](ahy_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [ahy_mrgsolve_model.R](ahy_mrgsolve_model.R) | mrgsolve ODE 모델 (19구획, 7시나리오, 용량-반응 분석) |
| [ahy_shiny_app.R](ahy_shiny_app.R) | Shiny 대시보드 (6탭 인터랙티브 시뮬레이터) |
| [ahy_references.md](ahy_references.md) | 참고문헌 45편 (PubMed 링크, 섹션별 분류) |

---

## 임상 치료 가이드라인 (Treatment Guidelines)

| 학회 | 치료 권고 (AHY) |
|------|----------------|
| ACR 2020 | SUA ≥9mg/dL + 동반질환 → 치료 고려 (약한 권고) |
| EULAR 2016 | 무증상 AHY에서 ULT는 일반적으로 권장하지 않음 |
| 중국 류마티스학회 2023 | SUA ≥8mg/dL + 심혈관/신장 위험 → 치료 고려 |
| 일본 통풍학회 | SUA ≥8mg/dL 지속 → 치료 권고 |

> **핵심 논란**: AHY에서 요산 강하 치료가 심혈관·신장 예후를 개선하는지에 대한 무작위 대조시험 증거가 아직 불충분함.

---

*생성일: 2026-06-20 | Claude Code Routine (CCR) | 질환 분류: 만성질환/대사*
