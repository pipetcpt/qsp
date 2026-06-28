# 급성 신손상 (Acute Kidney Injury, AKI) QSP Model

[![AKI Map](aki_qsp_model.png)](aki_qsp_model.svg)

**분류:** 신장학 / 중환자의학 (Nephrology / Critical Care)  
**디렉토리:** [`acute-kidney-injury/`](.)

---

## 병태생리 개요 (Pathophysiology)

급성 신손상(AKI)은 **48시간 이내 혈청 크레아티닌 ≥ 0.3 mg/dL 상승, 또는 7일 이내 기저치의 1.5배 이상 상승, 또는 6시간 이상 소변량 < 0.5 mL/kg/h** 로 정의되는 임상 증후군입니다(KDIGO 2012). 전 세계 입원 환자의 10–15%, ICU 환자의 50% 이상에서 발생하며 원내 사망률을 최대 5–10배 높입니다.

### 세 가지 주요 아형 (Three Subtypes)

| 아형 | 주요 기전 | 대표 원인 |
|------|----------|----------|
| **허혈-재관류 손상 (IRI)** | ATP 고갈 → 미토콘드리아 기능부전 → 아폽토시스/괴사/페롭토시스 | 대수술, 심인성 쇼크, 저혈량증 |
| **신독성 (NTX)** | 직접 세뇨관 세포 독성, 산화스트레스, GSH 고갈 | 시스플라틴, 아미노글리코시드, 조영제 |
| **패혈증 연관 (SA-AKI)** | LPS/DAMPs → TLR4 → NF-κB → 미세혈관 장애 + 사이토카인 폭풍 | 패혈증, 패혈성 쇼크 |

### 핵심 병태생리 경로

| 경로 | 메커니즘 |
|------|---------|
| **미세혈관 기능부전** | ET-1↑, Ang II↑, NO↓, TXA2↑ → 구심성 세동맥 수축 → GFR↓ → 수질 저산소증 |
| **세뇨관 세포 손상** | ATP 고갈 → mPTP 개방 → Cyt-c 방출 → Casp-9/3 → 아폽토시스; GPx4↓ → 지질과산화 → 페롭토시스 |
| **산화스트레스** | ROS↑ (미토콘드리아·NOX) → GSH 고갈 → Nrf2/HO-1 대항 부족 → 세포사멸 |
| **염증 폭풍** | DAMPs → TLR4 → NF-κB → IL-6, TNF-α, IL-1β↑ → 중성구 침윤 → NET형성 |
| **세뇨관 폐색 및 역류** | 솔가장자리 탈락 → 관내 캐스트 → 폐색 → filtrate backleak |
| **세뇨관 수복 실패** | G2/M 정지 → TGF-β1↑ → 근섬유모세포 활성화 → 섬유화 → AKI-to-CKD |

---

## KDIGO AKI 병기 (Staging)

| Stage | 혈청 크레아티닌 | 소변량 |
|-------|--------------|-------|
| **1** | ×1.5–1.9 기저치 또는 +0.3 mg/dL | < 0.5 mL/kg/h × 6–12h |
| **2** | ×2.0–2.9 | < 0.5 mL/kg/h × ≥12h |
| **3** | ×3.0 이상 또는 투석 | < 0.3 mL/kg/h × ≥24h 또는 무뇨 × ≥12h |

---

## 모델 산출물 (Deliverables)

| 파일 | 설명 |
|------|------|
| `aki_qsp_model.dot` | Graphviz 기계론적 지도 (12 클러스터, 100+ 노드) |
| `aki_qsp_model.svg` | 벡터 그래픽 (인터랙티브 뷰) |
| `aki_qsp_model.png` | 래스터 이미지 (150 dpi) |
| `aki_mrgsolve_model.R` | mrgsolve ODE 모델 (20 구획, 7 시나리오) |
| `aki_shiny_app.R` | Shiny 인터랙티브 대시보드 (6 탭) |
| `aki_references.md` | 참고문헌 55개 (섹션별 분류) |

---

## 모델 사양 (Model Specifications)

### 기계론적 지도
- **12 서브그래프 클러스터**: ①유발 인자·위험 요인, ②신장 미세혈관 기능부전, ③세뇨관 세포 손상 (아폽토시스/괴사/페롭토시스), ④산화스트레스 (ROS·GSH·Nrf2), ⑤선천 면역·염증 (NF-κB·NLRP3·IL-6·TNF-α), ⑥사구체 여과 및 소변 역학, ⑦AKI 바이오마커, ⑧Furosemide PK (OAT1/3·NKCC2), ⑨Vasopressors·NAC·CRRT PK/PD, ⑩세뇨관 수복·재생, ⑪부적응적 수복 → AKI-to-CKD, ⑫임상 결과
- **100+ 노드** 포함

### mrgsolve ODE 모델 (20 구획)
| 모듈 | 구획 |
|------|------|
| **Drug PK** | Furosemide central/peripheral/gut · Norepinephrine central · NAC central |
| **AKI 병태생리** | ATP · ROS · GSH · 세뇨관 세포 활력(TCV) · GFR · IL-6 · TNF-α |
| **바이오마커** | NGAL · KIM-1 · 혈청 크레아티닌 · Cystatin C |
| **수복·섬유화** | Repair Capacity · TGF-β1 · Myofibroblast · Fibrosis Index |

### 치료 시나리오 (7개)
1. **IRI — 무치료**: 자연 경과 대조군
2. **IRI + Furosemide 40mg IV q12h**: 이뇨제 전략 (NKCC2 억제)
3. **IRI + Norepinephrine**: 혈압 지지 → RBF 자동조절 회복
4. **IRI + NAC 예방투여**: GSH 보충 → ROS 제거 → 산화손상↓
5. **SA-AKI — 무치료**: 패혈증 자연경과
6. **SA-AKI + NE + Furosemide + CRRT**: 복합 중재 전략
7. **신독성 AKI (시스플라틴)**: 장기 AKI-to-CKD 전환 모니터링

### 임상시험 보정 파라미터
| 참고 임상시험 | 보정 데이터 |
|-------------|-----------|
| Mehta et al. Lancet 2015 | KDIGO staging criteria (sCr×1.5/2.0/3.0) |
| Gaudry et al. NEJM 2016 | Early CRRT vs delayed — 90d mortality ~48% |
| Zarbock et al. JAMA 2016 | Early RRT initiation benefit in AKI-3 |
| Mishra et al. JASN 2003 | NGAL rise within 2h of IRI |
| Han et al. Kidney Int 2002 | KIM-1 shedding kinetics |
| Meersch et al. ICM 2017 | TIMP-2·IGFBP7 prediction cutoff |
| Felker et al. NEJM 2011 | Furosemide dose-response in decongestion |

### Shiny 대시보드 (6 탭)
| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | AKI 아형 선택, 위험인자, KDIGO 병기 기준표, 개요 그래프 |
| ② Drug PK | Furosemide 혈중농도, NKCC2 억제율, NE/NAC 농도, 이뇨반응 |
| ③ 신장 바이오마커 | NGAL, KIM-1, Cystatin C vs Cr, 세뇨관 세포 활력, ROS/GSH |
| ④ 임상 엔드포인트 | eGFR 궤적, AKI Stage 시계열, 소변량, IL-6/TNF-α |
| ⑤ 시나리오 비교 | 전체 7개 시나리오 GFR·sCr·NGAL 비교 + 임상요약표 |
| ⑥ AKI-to-CKD | TGF-β1, 근섬유모세포, 섬유화 지수, 장기(30일) GFR, CKD 위험도 |

---

## 실행 방법 (Usage)

```r
# 1) mrgsolve 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("aki_mrgsolve_model.R")

# 2) Shiny 앱 실행
install.packages(c("shiny", "shinydashboard", "mrgsolve", "DT"))
shiny::runApp("aki_shiny_app.R")
```

```bash
# 3) 기계론적 지도 재렌더링
sfdp -Tsvg -Goverlap=prism aki_qsp_model.dot -o aki_qsp_model.svg
sfdp -Tpng -Goverlap=prism -Gdpi=150 aki_qsp_model.dot -o aki_qsp_model.png
```

---

## 핵심 약물 PK/PD 요약

| 약물 | 모델 특성 |
|------|---------|
| **Furosemide** | 2구획 PK · OAT1/3 활성분비 → 세뇨관내 농도 결정 · NKCC2 IC50=0.5 mg/L · AKI 시 OAT 발현↓ → 이뇨 저항 |
| **Norepinephrine** | 1구획 PK (t½ ~2.5분) · α1-수용체 → MAP↑ → 구심성 세동맥 자동조절 → RBF 회복 |
| **N-Acetylcysteine** | 1구획 PK · Cysteine 전구체 → GSH 보충 · thiol-NO 교환 → eNOS↑ |
| **CRRT** | 3 L/h (~50 mL/min) 크레아티닌 제거 · IL-6 흡착(0.5 L/h) · 체액 균형 조절 |

---

## 참고문헌 요약 (References Summary)

55개 PubMed 문헌, 16개 섹션으로 분류:
- 역학·부담(4), KDIGO 기준(2), IRI 병태생리(5), 신독성(3), 패혈증-AKI(3),
  산화스트레스·미토콘드리아(4), 염증(3), 바이오마커(5), Furosemide(3),
  Vasopressors(2), NAC(2), CRRT(3), AKI-to-CKD(5), QSP 모델링(4),
  임상시험·가이드라인(4), 바이오마커-유도 예방(3)

전체 목록: [`aki_references.md`](aki_references.md)

---

## 면책 조항

본 모델은 교육 및 연구 목적의 정성적·반정량적 QSP 모델이며, 임상 의사결정에 직접 사용해서는 안 됩니다.
