# 심인성 쇼크 (Cardiogenic Shock, CS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관·중환자

[![Cardiogenic Shock QSP Model](cs_qsp_model.png)](cs_qsp_model.svg)

## 개요 (Overview)

심인성 쇼크(Cardiogenic Shock, CS)는 심장 펌프 기능의 1차적 부전으로 인해 조직 관류가 부족해지는 임상 증후군입니다. SHOCK 등록 연구에 따르면 약 **80%가 급성 심근경색(특히 STEMI)** 으로부터 발생하며, 나머지는 전격성 심근염·Takotsubo·말기 DCM 악화·치명적 부정맥·기계적 합병증(유두근 파열·VSR·자유벽 파열)에서 비롯됩니다. 30일 사망률은 SCAI 분류상 **E 단계에서 ~67%, B 단계에서 ~5%**, 전체적으로 약 40-50%에 이릅니다.

핵심 병태생리는 **저심박출(CI<2.2) → 평균동맥압↓ → 조직 저관류 → 젖산↑ + 신경호르몬 폭주(NE 10-20×, RAAS) + DAMP-매개 SIRS 염증(NO 과잉으로 vasoplegia)** 이라는 다단계 악순환입니다. 최근 SCAI 2022 업데이트는 A(at risk)→E(extremis) 5단계로 분류해 치료 강도와 MCS 적용을 권고합니다. 치료의 중심은 **신속한 재관류(SHOCK trial 6년 RRR ~30%) · 노르에피네프린(1차 vasopressor, SOAP-II) · 강심제(dobutamine/milrinone DOREMI 비열등) · 기계적 순환보조(IABP-SHOCK II 음성, DanGer Shock 2024: Impella CP가 STEMI-CS에서 사망률 감소 입증)** 입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| LV 펌프 실패 | 심근 괴사·기절(stunning) → dP/dt↓ → SV<35mL → CI<2.2 | SCAI C-E 진입, 혈역학적 부전 |
| LVEDP/PCWP 상승 | 후방 부전 → 폐 충혈 → ARDS | 저산소혈증, 호흡부전 |
| 신경호르몬 폭주 | NE 10-20×, RAAS·AngII↑, β1 GRK2 하향조절 | 심근 산소 요구↑, cAMP 반응↓ |
| 미세순환 부전 | De Backer DBI 이질성·기능적 단락 | StO2↓, 젖산↑, DO2/VO2 불균형 |
| SIRS / iNOS-NO | DAMP·TNF-α·IL-6 → iNOS → NO 과잉 → SVR↓ (혼합형 쇼크) | 후기 vasoplegia, 30-40% |
| 다발성 장기부전(MOF) | AKI(KDIGO 2-3)·쇼크 간·DIC·허혈성 장염 | SCAI E 67% mortality |
| RV 부전 (하벽 STEMI 25%) | TAPSE↓·septal bowing → LV 전부하↓ | LV under-filling, 양심실 부전 |
| 보상 실패의 악순환 | 산증 + 베타1 desensitization → 강심제 효과↓ | 약물 내성, 진행성 사망 |

## 주요 약물 표적 (Drug Targets)

- **노르에피네프린(NE)**: α1 강력 + 약한 β1 → MAP·관류 회복 1차 선택(2018 AHA/ACC). EC50_α1 ≈ 0.04 µg/min
- **에피네프린**: α·β1·β2 균형. 심정지·아나필락시스/난치성 쇼크에 사용. Levy 2018 cohort에서 사망률 증가 신호
- **도파민**: 저용량 D1(신·내장 혈관확장)·고용량 α1. **SOAP-II에서 NE 대비 부정맥·사망률 증가 → 1차 선택에서 후퇴**
- **도부타민(Dobutamine)**: β1>β2 → cAMP·PKA·SERCA2a → 수축력↑. 2.5-20 µg/kg/min. β1 desensitization 제한
- **밀리논(Milrinone)**: PDE3 억제 → cAMP↑ → 강심·혈관 확장. **DOREMI(NEJM 2021): 도부타민과 비열등**
- **레보시멘단(Levosimendan)**: 트로포닌 C Ca²⁺ 감수성↑(에너지 비용 없음) + K-ATP 개방 → 혈관 확장. 활성 대사체 OR-1896 (반감기 ~80 h)
- **바소프레신(AVP)**: V1aR 작용. 카테콜아민 저항성 vasoplegia 또는 폐고혈압 동반 RV 부전에서 보조
- **오메캄티브 메카르빌**: 심근 미오신 활성제(GALACTIC-HF). 현재 만성 HFrEF 대상이며 CS 적응증 미확립

### 기계적 순환보조 (MCS)
- **IABP**: 이완기 풍선 팽창 → 관상 관류 ↑·후부하 ↓. IABP-SHOCK II에서 사망률 차이 없음(주의)
- **Impella CP/5.5**: 미세축류 펌프, LV 직접 unload 2.5-5.5 L/min. **DanGer Shock 2024(NEJM): 6개월 사망 26→18% (HR 0.74)**
- **VA-ECMO**: 4-7 L/min 전신 지원. **ECLS-SHOCK(2023): 30일 사망률 차이 없음**, ECMO-CS·EURO-SHOCK 역시 음성
- **장기 LVAD / 심장이식**: SCAI D-E 또는 회복 불가시 다리(bridge) 또는 종착(destination)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cs_qsp_model.dot](cs_qsp_model.dot) | Graphviz 기계론적 지도 소스 (120+ 노드 · 12 클러스터) |
| [cs_qsp_model.svg](cs_qsp_model.svg) | SVG 벡터 이미지 |
| [cs_qsp_model.png](cs_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cs_mrgsolve_model.R](cs_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획 · 7 치료 시나리오 · 약물 PK/PD + MCS) |
| [cs_shiny_app.R](cs_shiny_app.R) | Shiny 대시보드 (6 탭) |
| [cs_references.md](cs_references.md) | 참고문헌 (63편 PubMed) |

## mrgsolve 모델 (ODE Model)

- **구획 구조 (총 26)**: 약물 PK 8개(NE·EPI·Dobu·DA·Mil·Levo + OR-1896 + AVP), 혈역학·생리 12개(EF·CO·MAP·SVR·PCWP·HR·LAC·SvO2·CR·ALT·UO·STUN), 매개체 6개(NO 과잉·TNF-α·DAMP·RAAS·SNS·CUMHAZ)
- **수용체 표적**: α1·β1·β2·D1·V1a·PDE3·TroponinC·KATP. 각 수용체 점유율 Emax/EC50 모형
- **MCS**: IABP(+0.5 L/min 보조), Impella CP(+3.5 L/min, LV unload 0.25), VA-ECMO(+4.0 L/min, LV distension 위험)
- **시나리오 (7)**: (1) 무치료 기저선, (2) NE+Dobu(표준), (3) NE+Milrinone(DOREMI), (4) Levosimendan, (5) NE+Dobu+IABP, (6) NE+Impella CP, (7) NE+VA-ECMO
- **보정/근거**: SOAP-II, IABP-SHOCK II, ECLS-SHOCK, DanGer Shock, DOREMI, SURVIVE, SHOCK trial을 정성·반정량 보정
- **SCAI 단계 출력**: LAC·MAP·CO 기반 5단계(A→E) 산출. 단순화 Cox-like 사망률 해저드(NO 과잉·CO·MAP·LAC 가중)

## Shiny 대시보드 (Dashboard)

6개 탭: ① 환자 프로파일·SCAI 단계 · ② 약물 PK(8개 약물 농도) · ③ 혈역학 PD(MAP·CO·SVR·PCWP·HR·SvO2) · ④ 임상 엔드포인트(젖산·요량·Cr·ALT) + 생존 곡선 · ⑤ 7개 시나리오 비교 · ⑥ 염증·NO 바이오마커.

## 실행 방법 (Usage)

```r
library(mrgsolve)
source("cs_mrgsolve_model.R")
mod <- build_cs()
sim <- run_cs_scenarios(end = 72)
ggplot(sim, aes(time, MAP, colour = Scenario)) + geom_line()
# Shiny:
shiny::runApp("cs_shiny_app.R")
```

```bash
# 기계론적 지도 렌더링
dot -Tsvg cs_qsp_model.dot -o cs_qsp_model.svg
dot -Tpng -Gdpi=150 cs_qsp_model.dot -o cs_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [cs_references.md](cs_references.md) 참조 (총 63편 PubMed).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
