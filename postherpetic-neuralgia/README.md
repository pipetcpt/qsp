# Postherpetic Neuralgia (PHN) — QSP Model

대상포진 후 신경통(Postherpetic Neuralgia, PHN)에 대한 정량적 시스템 약리학(QSP)
모델 패키지입니다. 수두-대상포진 바이러스(VZV)의 재활성화 → 신경절 손상 → 말초·중추
민감화 → 다양한 통증 표현형으로 이어지는 병태생리와, 항바이러스/예방백신 (Shingrix)
부터 가바펜티노이드, TCA·SNRI, 5% 리도카인 패치, 8% 캡사이신 패치, NMDA 차단제,
오피오이드, 그리고 NaV1.7/NGF/P2X3 등 미래 표적까지 약리학을 한 곳에 통합합니다.

## 파일 구성 / Files

| 파일 | 설명 |
|---|---|
| `phn_qsp_model.dot` | Graphviz 기계론적 지도 — 13개 클러스터, 170+ 노드 (VZV 잠복-재활성, 신경 손상, 말초/중추 민감화, 신경염증, 임상 표현형, 백신·항바이러스·다중약물 약리, PK·안전성·환자 공변량) |
| `phn_qsp_model.svg` / `.png` | DOT 렌더링 결과 (svg + 150 dpi png) |
| `phn_mrgsolve_model.R` | mrgsolve QSP 모델 — 24개 ODE 컴파트먼트 (PK: GBP·PGB·AMI/NOR·DLX·LIDO·CAP·VAL·TRA·RZV; PD: VZV_LOAD, GANG_INJ, IENF, NAV_ACT, CSEN, MICROG, KCC2, NMDA_TONE, NGF, CMI, PAIN, ALLO, SLEEP, MOOD, AE_SED). 10개 시나리오 사전 정의 |
| `phn_shiny_app.R` | Shiny 대시보드 — 8개 탭 (Patient · PK · PD physiology · Clinical endpoints · Scenarios · Vaccine · Safety · Biomarkers/QST) |
| `phn_references.md` | PubMed 참고문헌 80개 (자연사·기전·백신·항바이러스·약물별 RCT·QSP/PK 모델·QST 표현형) |

## 사용법 / How to use

```r
# 1) DOT 렌더링 (이미 SVG/PNG 포함):
# dot -Tsvg phn_qsp_model.dot -o phn_qsp_model.svg
# dot -Tpng -Gdpi=150 phn_qsp_model.dot -o phn_qsp_model.png

# 2) mrgsolve 모델 로드 & 시나리오 시뮬레이션
source("phn_mrgsolve_model.R")
out <- mrgsim_e(mod, e_combo, end = 180*24, delta = 24) |> as.data.frame()

# 3) Shiny 대시보드
shiny::runApp("phn_shiny_app.R")
```

## 시나리오 (10) / Scenario library

1. **Placebo** — 자연 경과
2. **Valaciclovir 1 g q8h × 7 d** — 급성 zoster 항바이러스
3. **RZV (Shingrix)** — 0, 2-6개월 2회 접종 예방
4. **Gabapentin titration → 3600 mg/d**
5. **Pregabalin 75 → 150 bid (300 mg/d)**
6. **Amitriptyline 25 → 75 mg HS**
7. **Duloxetine 30 → 60 mg qd**
8. **Lidocaine 5% patch** 일 1회 12 h 적용
9. **Capsaicin 8% patch** 단회 60-min 적용, 90일 간격 재적용
10. **Combo** (AV + RZV + PGB + lidocaine + AMI) — 통합 관리 표준화

## 모델 보정 앵커 / Calibration anchors

- RZV 효능: HZ 97% (50-69세), 91% (≥70세), PHN 88-91% (ZOE-50/70)
- Valaciclovir 1 g tid → 급성 통증 회복 시간 ~30% 단축 (Beutner 1995, Tyring 1995)
- Pregabalin 300-600 mg → 50% 통증감소 책임자율 ~35-50% (Dworkin 2003)
- Gabapentin 1800-3600 mg → NRS ~30% 감소 (Rice 2001)
- Amitriptyline 25-100 mg HS → NNT ~2.7 (Watson 1982)
- Capsaicin 8%: 단회 적용 후 12 주간 NRS ~30% 감소 (Backonja 2008 / STRIDE)
- Lidocaine 5% patch: allodynia 30-50% 감소, 매우 낮은 systemic Css (Galer 2002)

## 모델 범위 한계 / Caveats

- 통증 ODE는 phenomenological — in-silico 시나리오 비교 및 trial design용
- PK는 단일/2-구획 단순화, popPK 변동성 미포함
- VZV reactivation drive는 단일 first-order 변수로 추상화
- 캡사이신 IENF 영향은 약리학적 defunctionalization 평균 효과로 가정
