# 자가면역 다발성 내분비병증 (APS / APECED)
## Autoimmune Polyendocrinopathy Syndrome Type 1 — QSP Model

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![APS QSP Model](aps_qsp_model.png)](aps_qsp_model.svg)

---

## 개요 (Overview)

자가면역 다발성 내분비병증-칸디다증-외배엽이영양증(APECED, Autoimmune Polyendocrinopathy-Candidiasis-Ectodermal Dystrophy)은 **AIRE** (autoimmune regulator) 유전자의 돌연변이로 인해 중추 면역관용이 파괴되어 여러 내분비 기관이 자가면역 공격을 받는 희귀 단일 유전자 질환입니다. APS Type 1으로도 불리며, 3대 주요 구성 요소는 만성 피부칸디다증(CMC), 부갑상선기능저하증, 애디슨병입니다.

- **유병률**: 핀란드 1/25,000, 사르데냐 1/14,400, 이란 유대인 1/9,000
- **발병 연령**: 주로 소아기 (CMC 첫 징후: 평균 5세)
- **유전**: 상염색체 열성 (AIRE 유전자, 21q22.3), >300 종의 돌연변이 확인
- **진단 기준**: 3개 주 증상 중 2개 이상 또는 가족 내 APECED 병력

---

## 발병기전 요약 (Pathophysiology Summary)

| 단계 | 기전 | 결과 |
|------|------|------|
| AIRE 기능 소실 | 흉선 수질 상피세포(mTEC)에서 조직 특이 항원(TSA) 발현 불능 | 자가반응 T세포 음성선택 실패 |
| 중추 관용 붕괴 | 자가반응 CD4+/CD8+ T세포 말초로 방출 | 다장기 자가면역 공격 |
| 자가항체 생성 | 자가반응 T세포 도움 → B세포 → IgG 자가항체 | 항21-OH(부신), 항NALP5(부갑상선), 항GAD65(췌장), 항TPO(갑상선) |
| 항IFN-α 항체 | APS1 특이 병리 (항IFN-α/IFN-ω 중화항체) | CMC 위험 증가, 바이러스 감염 취약 |
| 장기 파괴 | CTL 직접 손상 + 보체 매개 파괴 | 부신피질, 부갑상선, 베타세포, 갑상선 기능 상실 |

---

## 모델 구조 (Model Architecture)

### 기계론적 지도 (Mechanistic Map)

- **클러스터**: 11개 서브그래프
- **노드 수**: 140+ (면역계, 4개 표적 장기, 약물 PK/PD, 임상 끝점)
- **주요 경로**: AIRE→흉선 음성선택→자가반응 T세포→자가항체→장기 파괴→임상 표현형

| 클러스터 | 내용 |
|---------|------|
| 1. AIRE & 흉선 | AIRE 유전자/단백질, mTEC, TSA 발현, 음성선택, Treg 생성 |
| 2. 말초 관용 | Treg, CTLA-4, PD-1, B 조절세포, 에너지 |
| 3. 자가항체 | B세포 활성화, GC 반응, 형질세포, 10종 주요 자가항체 |
| 4. 사이토카인 | IFN-γ, IL-12/23, IL-17, TNF-α, JAK-STAT, NF-κB |
| 5. 부신 병리 | 부신피질 파괴, 코르티솔, ACTH, HPA축, 알도스테론 |
| 6. 부갑상선 병리 | PTG 파괴, PTH, Ca²⁺, 비타민D 활성화, 신경근육 흥분성 |
| 7. 췌장 병리 | 베타세포 파괴, 인슐린, 혈당, HbA1c |
| 8. 갑상선 병리 | 갑상선 파괴, TSH, FT4, 대사 효과 |
| 9. 호르몬 대체 PK | 히드로코르티손, 플루드로코르티손, 칼시트리올, 레보티록신 |
| 10. 면역억제제 PK | 사이클로스포린 A, 아바타셉트, 리툭시맙, 토파시티닙, IVIG |
| 11. 임상 끝점 | 장기 기능 점수, 자가항체 역가, QoL, APS 복합 점수 |

---

### mrgsolve ODE 모델 (20 구획)

| # | 구획 | 단위 | 정상값 |
|---|------|------|--------|
| 1 | AIRE_func | 0-1 scale | 1.0 (정상) |
| 2 | AutoT_pool | cells/µL | 2.0 |
| 3 | Treg_pool | cells/µL | 15.0 |
| 4 | AutoAb_adren | U/mL (항21-OH) | 1.0 |
| 5 | AutoAb_PTG | U/mL (항NALP5) | 1.0 |
| 6 | AutoAb_beta | U/mL (항GAD65) | 1.0 |
| 7 | AutoAb_thy | U/mL (항TPO) | 1.0 |
| 8 | Adrenal_fn | % (0-100) | 100 |
| 9 | Cortisol_c | µg/dL | 12.0 |
| 10 | PTG_fn | % (0-100) | 100 |
| 11 | PTH_plasma | pg/mL | 40.0 |
| 12 | Ca_serum | mg/dL | 9.4 |
| 13 | Beta_mass | % (0-100) | 100 |
| 14 | Insulin_p | pmol/L | 60.0 |
| 15 | Glucose_p | mg/dL | 90.0 |
| 16 | Thyroid_fn | % (0-100) | 100 |
| 17 | TSH_plasma | mIU/L | 2.0 |
| 18 | FT4_plasma | ng/dL | 1.2 |
| 19-23 | Drug_CsA/Aba/RTX/JAKi/HC | ng/mL 또는 µg/dL | 0 |

---

### 치료 시나리오 (7개)

| # | 시나리오 | AIRE 중증도 | 치료 |
|---|---------|------------|------|
| 1 | 자연경과 (중증) | 90% 기능 소실 | 없음 |
| 2 | HRT 단독 | 90% 기능 소실 | HC 20mg/day |
| 3 | HRT + CsA | 90% 기능 소실 | HC + 사이클로스포린 A 3.5mg/kg/day |
| 4 | HRT + 아바타셉트 | 90% 기능 소실 | HC + Abatacept 125mg SC/주 |
| 5 | HRT + 리툭시맙 | 90% 기능 소실 | HC + RTX 375mg/m² q6mo |
| 6 | HRT + 토파시티닙 | 90% 기능 소실 | HC + JAKi 10mg/day |
| 7 | 조기 개입 (경증) | 30% 기능 소실 | HC 15mg/day + 조기 시작 |

---

## 주요 약물 PK/PD 파라미터

| 약물 | 기전 | 주요 PK | 임상 효과 |
|------|------|---------|---------|
| 히드로코르티손 20mg/day | GR 작용제 | F=0.96, t½=1.5h, CL=1.1L/min | AM 코르티솔 정상화 (8-20 µg/dL) |
| 플루드로코르티손 100µg/day | MR 작용제 | F=0.99, t½=3.5h | Na+/K+ 균형 정상화 |
| 칼시트리올 0.5µg/day | VDR 작용제 | t½=5-8h | Ca²⁺ 정상화 (8.5-10.5 mg/dL) |
| 레보티록신 75-100µg/day | 갑상선호르몬 대체 | F=0.75-0.80, t½=9일 | TSH 정상화 (0.4-4.0 mIU/L) |
| 사이클로스포린 A 3.5mg/kg/day | 칼시뉴린 억제 | F=0.30(가변), t½=8-24h | AutoT 억제 ~50% at Cp=150ng/mL |
| 아바타셉트 125mg SC/주 | CD28/B7 차단 | F=0.79, t½=13일 | Treg/AutoT 비율 개선 |
| 리툭시맙 375mg/m² q6mo | 항CD20 (B세포 고갈) | IV, t½=21일 | 자가항체 ~85% 감소 |
| 토파시티닙 10mg/day | JAK1/3 억제 | F=0.74, t½=3h | IFN-γ/IL-17 신호 70% 차단 |

---

## 임상 끝점 (Clinical Endpoints)

| 끝점 | 정상 범위 | APS1 목표 |
|------|----------|----------|
| AM 코르티솔 (8시) | 8-20 µg/dL | ≥8 µg/dL on HRT |
| 혈청 Ca²⁺ (교정) | 8.5-10.5 mg/dL | 8.0-9.0 mg/dL (안전 범위) |
| 혈청 PTH | 15-65 pg/mL | 측정 불가 (기능 소실) |
| 공복 혈당 | 70-100 mg/dL | <130 mg/dL (T1DM) |
| HbA1c | <5.7% | <7.0% (T1DM 목표) |
| TSH | 0.4-4.0 mIU/L | 0.5-2.5 mIU/L (LT4 치료 중) |
| Free T4 | 0.8-1.8 ng/dL | 1.0-1.5 ng/dL |
| Anti-21-OH Ab | <1 U/mL | 모니터링 (>1 → 부신 위험) |
| 장기 기능 점수 | 100% | >50% 유지 목표 |

---

## 연간 스크리닝 패널 (APS1 Lifetime Monitoring)

APS1 환자는 평생 새로운 구성 요소가 추가될 수 있어 매년 자가항체 스크리닝이 필수입니다.

| 항체 | 대상 장기 | 임계값 | 조치 |
|------|---------|--------|------|
| 항-21-hydroxylase IgG | 부신 | >1 U/mL | AM 코르티솔, ACTH 자극시험 |
| 항-NALP5 IgG | 부갑상선 | >10 U/mL | 혈청 Ca, PTH 측정 |
| 항-GAD65 IgG | 췌장 β세포 | >5 U/mL | 공복혈당, OGTT |
| 항-TPO/TG IgG | 갑상선 | >34 IU/mL | TSH 측정 |
| 항-IFN-α 중화항체 | 전신 바이러스 방어 | 양성 | APS1 확진 마커 |
| 항-위벽세포/내인성인자 Ab | 위장 | 양성 | 비타민 B12 수준 |
| 항-17α-OH IgG | 생식샘 | 양성 | LH/FSH/에스트라디올 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| [aps_qsp_model.dot](aps_qsp_model.dot) | Graphviz 기계론적 지도 소스 (140+ 노드, 11 클러스터) |
| [aps_qsp_model.svg](aps_qsp_model.svg) | SVG 벡터 이미지 (확대 가능, 인터랙티브) |
| [aps_qsp_model.png](aps_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [aps_mrgsolve_model.R](aps_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 7 시나리오) |
| [aps_shiny_app.R](aps_shiny_app.R) | Shiny 대시보드 (7탭: 환자/면역/PK/내분비/끝점/시나리오/바이오마커) |
| [aps_references.md](aps_references.md) | 참고문헌 60편 (PubMed 링크 포함) |

---

## 주요 참고문헌

- Husebye ES et al. *N Engl J Med.* 2018;378:1132–1141
- Perheentupa J. *J Clin Endocrinol Metab.* 2006;91:2843–2850
- Anderson MS et al. *Science.* 2002;298:1395–1401
- Alimohammadi M et al. *N Engl J Med.* 2008;358:1018–1028
- Landegren N et al. *Sci Rep.* 2016;6:20104
