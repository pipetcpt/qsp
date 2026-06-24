# 자궁내막증 (Endometriosis) QSP Model

## Overview / 개요

### English

Endometriosis is a chronic inflammatory estrogen-dependent disease affecting approximately **10% of reproductive-age women** (approximately 190 million women worldwide). It is characterised by the presence of endometrial-like tissue outside the uterine cavity, most commonly on the ovaries, fallopian tubes, peritoneum, and rectovaginal septum. The disease causes significant morbidity through chronic pelvic pain, dysmenorrhea, dyspareunia, and infertility. Diagnosis is typically delayed 7–10 years from symptom onset.

The pathophysiology involves a complex interplay of retrograde menstruation, impaired immune surveillance, local estrogen biosynthesis via aromatase upregulation, progesterone resistance, neuroangiogenesis, and a self-sustaining PGE2–aromatase positive feedback loop. The QSP model in this directory captures these mechanisms quantitatively using an ODE-based framework spanning the HPO axis, lesion dynamics, inflammation, and pain signalling.

### 한국어

자궁내막증은 만성 염증성 에스트로겐 의존성 질환으로, 가임기 여성의 약 **10%** (전 세계 약 1억 9천만 명)에서 발생합니다. 자궁 내막과 유사한 조직이 자궁강 밖(난소, 복막, 직장-질 사이 등)에 존재하며, 만성 골반통·월경통·성교통·불임을 야기합니다. 증상 발현부터 진단까지 평균 7~10년이 소요됩니다.

발병기전은 역행성 월경, 면역 감시 기능 장애, 방향화효소(aromatase) 과발현에 의한 국소 에스트로겐 생합성, 프로게스테론 저항성, 신경혈관신생, PGE2–aromatase 양성 피드백 루프 등이 복합적으로 작용합니다.

---

## Mechanistic Map / 기계론적 지도

| 경로 (Pathway) | 핵심 분자/세포 (Key Nodes) | 임상 결과 (Clinical Outcome) |
|---|---|---|
| HPO 축 (HPO Axis) | GnRH → LH/FSH → E2 → 음성 피드백 | 배란 주기 조절, E2 공급 |
| 역행성 월경 (Retrograde Menstruation) | 내막 세포 → 복강 → 착상 | 병소 형성 |
| 방향화효소 과발현 (Aromatase Upregulation) | SF-1 ↑, PGE2 → CYP19A1 | 국소 E2 과잉생성 |
| 프로게스테론 저항성 (Progesterone Resistance) | PR-B ↓, HSD17B2 ↓ | 내막 증식 억제 실패 |
| PGE2 피드백 (PGE2 Loop) | COX-2 → PGE2 → EP2/EP4 → aromatase | 염증 지속, 증식 촉진 |
| NK 세포 기능 이상 (NK Dysfunction) | NK 세포 감소/기능 저하 | 이소 세포 clearance 실패 |
| 신경침투/NGF (Neural Invasion) | NGF → TrkA → 통증 신호 | 이상통증, 중추 감작 |
| 산화 스트레스 (Oxidative Stress) | ROS → NF-κB → IL-6/IL-8 | 염증 증폭, 세포 생존 촉진 |
| 혈관신생 (Angiogenesis) | VEGF → 신생혈관 | 병소 영양 공급 |
| 거대세포 침윤 (Macrophage Infiltration) | M2 대식세포 → TGF-β | 면역 억제, 섬유화 |

---

## Drug PK/PD Parameters / 약물 PK/PD 파라미터

| 약물 (Drug) | 기전 (Mechanism) | 주요 PK 파라미터 | 임상 효과 (Clinical Effect) |
|---|---|---|---|
| **Leuprolide depot 3.75 mg/월** | GnRH 작용제 → 탈감작 → 뇌하수체 억제 | t½ = 3 h (활성형); depot 방출 ~4주; Tmax ~4 h | E2 < 20 pg/mL (거세 수준); flare-up 1–2주 |
| **Elagolix 150 mg/일** | GnRH 길항제 → 즉각적 FSH/LH 억제 (부분 억제) | t½ = 4–6 h; Tmax = 1 h; 선형 PK | E2 부분 억제 (~50%); 월경 유지 가능 |
| **Elagolix 200 mg BID** | GnRH 길항제 → 완전 억제 | t½ = 4–6 h; BID 투여 | E2 완전 억제 (< 12 pg/mL); 무월경 유발 |
| **Dienogest 2 mg/일** | 선택적 프로게스틴 → 내막 위축, E2 억제 (부분) | t½ = 9–10 h; Tmax = 1.5 h; F ≈ 91% | E2 ~60% 감소; 병소 위축; 통증 감소 |
| **Letrozole 2.5 mg/일 + Add-back** | 비스테로이드 방향화효소 저해제 | t½ = 48 h; Tmax = 1 h; 간 대사 | 말초 E2 > 90% 감소; 골 손실 예방 위해 add-back 필요 |
| **Combined OCP (EE/DRSP)** | 시상하부-뇌하수체 억제; 내막 위축 | EE: t½ = 6–12 h; DRSP: t½ = 30 h | E2 부분 감소; 통증 완화; 내막 얇아짐 |
| **NSAIDs (Naproxen 500 mg BID)** | COX-1/2 억제 → PGE2 감소 | t½ = 14 h; Tmax = 2–4 h | 통증 완화; 병소 증식에는 영향 제한적 |

---

## Pathophysiology Table / 발병기전 요약

| 경로 | 핵심 메커니즘 | 임상 이상 |
|---|---|---|
| HPO 축 조절 이상 | GnRH pulse 변화 → FSH/LH 과잉 → 난소 E2 과잉 생성 | 과도한 내막 증식, 불규칙 월경 |
| 역행성 월경 | 월경혈 내 생존 내막 세포 → 복강 내 착상 및 증식 | 복막/난소/직장 병소 형성 |
| 방향화효소 과발현 | SF-1 이소성 발현 → CYP19A1 유도 → 국소 E2 생합성 | 시스테믹 치료에 불완전 반응 |
| 프로게스테론 저항성 | PR-B 발현 저하, MAPK 과활성화 → progestin 효과 감소 | 프로게스틴 치료 반응 저하 |
| PGE2–Aromatase 양성 피드백 | PGE2 → EP2/EP4 → cAMP → SF-1 → aromatase → E2 → COX-2 → PGE2 ↑ | 자기 증폭 염증/증식 루프 |
| NK 세포 기능 이상 | NK 세포 수 감소, cytotoxicity 저하 → 이소 세포 생존 | 면역 제거 실패 → 병소 지속 |
| 신경침투/NGF 과발현 | 이소 내막 → NGF → TrkA 발현 감각신경 침투 → 중추 감작 | 만성 통증, 이상통증, 통각과민 |
| 산화 스트레스 | 복막액 내 철 과잉, ROS → NF-κB 활성화 | IL-6, IL-8, TNF-α 증가 → 병소 성장 |

---

## Treatment Scenarios / 치료 시나리오

| 시나리오 | 약물 | 목표 | 예측 효과 |
|---|---|---|---|
| 1 | 치료 없음 (No Treatment) | 기저 질환 경과 관찰 | 병소 서서히 증가; 통증 지속 |
| 2 | Leuprolide depot (GnRH 작용제) | 완전 E2 억제, 병소 위축 | E2 < 20 pg/mL; 병소 40–60% 감소; BMD 손실 위험 |
| 3 | Elagolix 150 mg/일 | 부분 E2 억제, 통증 조절 | E2 ~50% 감소; 월경통 개선; BMD 영향 최소 |
| 4 | Elagolix 200 mg BID | 완전 E2 억제, 강력한 통증 조절 | E2 < 12 pg/mL; 무월경; BMD 손실 주의 (2년 미만 사용) |
| 5 | Dienogest 2 mg/일 | 장기 유지, 내막 위축 | E2 부분 억제; 병소 감소; BMD 상대적 보존 |
| 6 | Letrozole 2.5 mg/일 + Add-back | 말초 aromatase 차단 + 골 보호 | E2 > 90% 감소; 재발성/수술 후 잔류 병소 치료 |
| 7 | 복합 경구피임약 (Combined OCP) | 호르몬 주기 억제, 통증 완화 | E2 중등도 감소; 장기 사용 가능; 임신 원하지 않는 경우 |

---

## File List / 파일 목록

| 파일 | 설명 |
|---|---|
| `endo_shiny_app.R` | Shiny 인터랙티브 대시보드 (6탭) |
| `endo_references.md` | 참고문헌 목록 (44개, 섹션별 분류) |
| `README.md` | 이 파일 — 모델 개요 및 파라미터 요약 |

> **참고**: DOT/SVG/PNG 기계론적 지도 파일 및 mrgsolve ODE 모델 파일은 추후 추가 예정입니다.

---

## Key References / 주요 참고문헌

1. **Giudice LC, Kao LC.** "Endometriosis." *Lancet* 2004;364(9447):1789–1799. [PMID: 15488215](https://pubmed.ncbi.nlm.nih.gov/15488215/)

2. **Zondervan KT, Becker CM, Missmer SA.** "Endometriosis." *N Engl J Med* 2020;382(13):1244–1256. [PMID: 31777792](https://pubmed.ncbi.nlm.nih.gov/31777792/)

3. **Taylor HS et al.** "Treatment of Endometriosis-Associated Pain with Elagolix, an Oral GnRH Antagonist." *N Engl J Med* 2017;377(1):28–40. [PMID: 28985706](https://pubmed.ncbi.nlm.nih.gov/28985706/)

4. **Strowitzki T et al.** "Dienogest is as effective as leuprolide acetate in treating the painful symptoms of endometriosis: a 24-week, randomized, multicentre, open-label trial." *Hum Reprod* 2010;25(3):633–641. [PMID: 20188814](https://pubmed.ncbi.nlm.nih.gov/20188814/)

5. **Bulun SE et al.** "Molecular basis for treating endometriosis with aromatase inhibitors." *Hum Reprod Update* 2000;6(5):413–418. [PMID: 12065431](https://pubmed.ncbi.nlm.nih.gov/12065431/)
