# 지역사회획득 폐렴 (Community-Acquired Pneumonia, CAP) — QSP 모델

<p align="center">
  <a href="cap_qsp_model.svg"><img src="cap_qsp_model.png" width="720" alt="CAP QSP mechanistic map"></a>
</p>

폐렴구균을 원형으로 하는 지역사회획득 폐렴의 정량적 시스템 약리학(QSP) 모델입니다.
세균 부하, **살균 방식에 따라 달라지는 PAMP 방출**, 숙주 염증 증폭, 폐포-모세혈관 장벽
파괴와 폐내 션트, 전신 장기부전, 그리고 세프트리악손·아목시실린·아지트로마이신·
레보플록사신·하이드로코르티손의 PK/PD를 하나의 미분방정식계로 묶었습니다.

| 파일 | 내용 |
|---|---|
| [`cap_qsp_model.dot`](cap_qsp_model.dot) | 기계론적 지도 소스 (14개 클러스터, 150+ 노드, 270+ 엣지) |
| [`cap_qsp_model.svg`](cap_qsp_model.svg) / [`cap_qsp_model.png`](cap_qsp_model.png) | 렌더링된 지도 (Graphviz `dot`, PNG 150 dpi) |
| [`cap_mrgsolve_model.R`](cap_mrgsolve_model.R) | mrgsolve 모델 — **38개 ODE 구획**, 9개 사전정의 시나리오 |
| [`cap_shiny_app.R`](cap_shiny_app.R) | Shiny 대시보드 — 8개 탭 |
| [`cap_references.md`](cap_references.md) | PubMed 검증 문헌 71편 (PMID 전수 대조) |

---

## 이 모델이 하려는 말 (The three structural commitments)

교과서적 CAP 모델은 "항생제가 균을 죽이면 환자가 좋아진다"로 끝납니다. 그 그림 하나로는
문헌이 요구하는 사실 여섯 가지를 동시에 만족시킬 수 없습니다. 그래서 세 가지 구조적
선택을 했습니다.

### ① fT>MIC와 fAUC/MIC는 입력이 아니라 **적분 결과**다

모델 어디에도 "β-락탐은 시간의존적, 퀴놀론은 농도의존적"이라고 적혀 있지 않습니다.
두 약 모두 **같은** 포화형 Emax 함수를 씁니다 — 유리 ELF 농도에 대해, EC50은 그 균주의
MIC에 비례. 유일한 차이는 Hill 계수입니다(세프트리악손 1.5, 레보플록사신 2.6).
`TAM`(fT>MIC)과 `AUCF`(fAUC/MIC)는 그 곡선을 시간 적분하는 **누적 구획**이며,
따라서 두 지표는 파라미터가 아니라 관측량입니다.

### ② 살균은 PAMP의 **원천**이지 소멸원이 아니다

β-락탐은 PBP를 억제하고 LytA 자가용해효소를 작동시켜 균을 터뜨립니다. 그 결과 세포벽
(LTA·펩티도글리칸)과 뉴몰라이신이 폐포로 한꺼번에 쏟아집니다. 비용균성 살균
(마크로라이드, 그리고 정도는 덜하지만 퀴놀론)은 훨씬 적게 방출합니다.

```
dPAMP/dt = YLIVE·(살아있는 균의 증식) + YHOST·(탐식 살균)
         + (YLYS_BL·k_βlactam + YLYS_FQ·k_quinolone + YLYS_MAC·k_macrolide)·BE
         − KPAMP·PAMP
```

`YLYS_BL = 1.00`, `YLYS_FQ = 0.35`, `YLYS_MAC = 0.12`.
이 항 하나가 **첫 투여 후 6-24시간의 염증 급등과 일시적 임상 악화**를 만들고,
동시에 **마크로라이드 역설**을 설명합니다 — 중증 폐렴구균 CAP에서 균주가 마크로라이드
**내성**이어도 병용 시 생존 이득이 관찰되는 이유는, 그 이득이 살균 경로가 아니라
독소/사이토카인 경로(`MAC_PLY`, `MAC_ANTI`, `MAC_MIG`)로 흐르기 때문입니다.
`MIC_AZI = 64`로 두면 `kazi ≈ 0`이 되지만 세 개의 면역조절 항은 그대로 남습니다.

이 구조는 반증 가능합니다: `YLYS_BL = 0`으로 두면 초기 급등도, 마크로라이드 이득의
상당 부분도 함께 사라집니다.

### ③ 코르티코스테로이드는 두 팔을 가지며, **순효과의 부호가 조건부**다

하나의 효과구획 `HCE`에서 두 방향이 갈라집니다.

| 팔 | 항 | 방향 |
|---|---|---|
| NF-κB transrepression | `GC_ANTI` (`IMAX_GC` 0.72) | 사이토카인 증폭기·장벽 손상·션트·혈관마비 ↓ (**이득**) |
| 옵소닌 탐식 억제 | `GC_PHAG` (`IMAX_GCP` 0.45) | 숙주 살균능 ↓ (**손해**) |

따라서 순효과는 "항생제가 이미 균을 없애서 숙주 살균이 **불필요해졌는가**"에 달려 있습니다.
감수성 균주 + 유효 β-락탐이면 스테로이드는 보호적이고(CAPE-COD 방향), 내성 균주이거나
바이러스 폐렴이라 항생제가 겨눌 표적이 없으면 **같은 용량이 해롭습니다**.

시나리오 6과 7은 스테로이드 용량·숙주·모든 파라미터가 동일하고 오직 MIC만 다릅니다.
그런데도 `Mortality_pct`의 부호가 뒤집힙니다. `IMAX_GCP = 0`으로 두면 이 주장 자체가
사라지므로, 역시 반증 가능합니다.

---

## 기계론적 지도 (14개 클러스터)

1. 숙주 감수성·침입 경로 — 연령·동반질환·흡인·점액섬모 청소·분비형 IgA·백신
2. 병원체·병독인자 — 협막, 뉴몰라이신, LytA, PspA/PspC, NanA, 비정형균, 바이러스
3. 폐포 상피·장벽 — AT1/AT2, 밀착연접, 계면활성제, ENaC/Na-K-ATPase, 유리질막
4. 선천면역 인지 — TLR2/4/9, NOD2, NLRP3, cGAS-STING, MyD88-NF-κB, 가스더민-D
5. 세포 효과기 — 폐포대식세포, 호중구, NET, MPO/엘라스타제, 단핵구 M1/M2, 효페로사이토시스
6. 사이토카인 네트워크 — TNF-α, IL-1β, IL-6, IL-8, CXCL1/2, CCL2, IL-17A, IL-10, TGF-β
7. 급성기 반응·바이오마커 — CRP, PCT, MR-proADM, 림프구감소
8. 가스교환·호흡역학 — V/Q 불균형, 션트, HPV, 사강, 순응도, PaO₂/FiO₂
9. 전신 파급·장기부전 — 균혈증, 내피 활성화, 혈관마비, 젖산, AKI, DIC, 심근손상, SOFA
10. 항생제 PK — 흡수, 단백결합, ELF 이행, 폐조직 축적, 신/담즙 청소
11. 항생제 PD·내성 — PBP2x/2b, ermB/mefA, GyrA/ParC, 유출펌프, 잠복균, **용균 방식**
12. 보조요법 — 스테로이드(GR·전사억제·면역억제), 산소/HFNC/NIV, 수액/승압제, 백신
13. 해소·복구 — 효페로사이토시스, SPM, AT2→AT1, 기질 재형성, 기질화 폐렴
14. 임상 결과 — 발열 소실, 임상안정(Halm), 치료실패, ICU 이송, 재원기간, 28일 사망, 재발

---

## mrgsolve 모델 구조 (38 ODE)

| 구획군 | 개수 | 상태변수 |
|---|---|---|
| 약동학 | 11 | `CEFC` `CEFP` `AMXD` `AMXC` `AZID` `AZIC` `AZIT` `LVXD` `LVXC` `HCC` `HCE` |
| 병원체 | 5 | `BE` 폐포 세포외균 · `BP` 잠복/세포내균 · `BB` 균혈증 · `PAMP` · `VIR` |
| 면역 | 9 | `AM` `NB` `NL` `TNF` `IL1` `IL6` `IL8` `IL10` `SPM` |
| 바이오마커 | 3 | `CRP` `PCT` `LYM` |
| 폐 | 4 | `PERM` `EDEM` `SURF` `FIBP` |
| 전신 | 3 | `MAPD` `LACT` `AKID` |
| 누적기 | 3 | `TAM` (fT>MIC) · `AUCF` (fAUC) · `CHZ` (누적 사망 위험) |

시간 단위는 **시간(hour)** 입니다 — 투여간격 내부의 노출 변화와 투여 직후 PAMP 급등을
분해하는 것이 모델의 요점이기 때문입니다. 기본 시뮬레이션 구간은 336 h(14일)입니다.

산출 지표: `PaO2_FiO2`, `SpO2_pct`, `Shunt_fraction`, `Temperature`, `Resp_rate`,
`Heart_rate`, `MAP_mmHg`, `Lactate`, `SOFA_score`, `CURB65`, `Clinically_stable`
(Halm 기준 그대로 구현), `Mortality_pct = 100·(1 − e^(−CHZ))`.

### 9개 시나리오

| # | 시나리오 | 보여주는 것 |
|---|---|---|
| 1 | 무치료 자연경과 | 숙주 단독 청소 → 7-9일 "위기", 높은 누적 위험 |
| 2 | 아목시실린 1 g q8h (경증 외래) | 경구 β-락탐으로도 충분한 경우 |
| 3 | 세프트리악손 2 g q24h (병동) | 기준 단독요법, fT>MIC ~100% |
| 4 | 세프트리악손 + 아지트로마이신 | 균 곡선은 거의 동일한데 PAMP·IL-8·션트가 낮아짐 |
| 5 | 레보플록사신 750 mg qd | 같은 세균학적 종점, 다른 PK/PD 경로 (fAUC/MIC ~70) |
| 6 | 중증 CAP + 하이드로코르티손 (CAPE-COD) | 스테로이드 **보호** 방향 |
| 7 | 동일 스테로이드 + 페니실린 비감수성 균주 | 스테로이드 **유해** 방향 — 부호 반전 |
| 8 | 세프트리악손 5일 단기요법 | 안정화 후 중단의 안전성, 재발은 `BP`가 결정 |
| 9 | 항생제 18시간 지연 (중증) | 지연 중 누적된 위험은 이후 회복되지 않음 |

### 돌려볼 만한 민감도 손잡이

| 파라미터 | 조작 | 결과 |
|---|---|---|
| `YLYS_BL` | 1.0 → 0.2 | 투여 직후 염증 급등 소멸, 마크로라이드 이득 대부분 소실 |
| `IMAX_GCP` | 0.45 → 0 | 스테로이드 유해 팔 제거 → 시나리오 7이 더 이상 해롭지 않음 |
| `MIC_AZI` | 0.12 → 64 | ermB 내성 — 살균은 죽고 면역조절은 살아남음 |
| `KPERS` | 0.010 → 0 | 단기요법 후 재발 자체가 사라짐 |
| `CRCL` | 90 → 25 | 레보플록사신 축적, 세프트리악손은 담즙 경로로 빠져나가 거의 불변 |
| `VIRAL` | 0 → 1 | PCT는 낮은데 CRP는 상승 — 두 마커가 갈라짐 |

---

## 사용법

```r
# 1) 모델 시뮬레이션
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread_cache("cap_mrgsolve_model.R")

# 세프트리악손 2 g IV q24h × 7일 + 아지트로마이신 500 mg × 5일, 6시간째 첫 투여
ev_abx <- ev(amt = 2000, cmt = "CEFC", time = 6, ii = 24, addl = 6) +
          ev(amt = 500,  cmt = "AZID", time = 6, ii = 24, addl = 4)

out <- mod |> mrgsim(events = ev_abx, end = 336, delta = 0.25) |> as_tibble()
ggplot(out, aes(time/24, Bacteria_log10)) + geom_line()

# 2) 9개 시나리오 일괄 실행 — 모델 파일 하단 주석의 CAP_simulate_scenarios() 참조

# 3) Shiny 대시보드
shiny::runApp("cap_shiny_app.R")
```

지도 재렌더링:

```bash
dot -Tsvg cap_qsp_model.dot -o cap_qsp_model.svg
dot -Tpng -Gdpi=150 cap_qsp_model.dot -o cap_qsp_model.png
```

---

## 검증 앵커

- 무치료 자연경과가 항생제 이전 시대 폐렴구균 폐렴의 7-9일 "위기"를 재현하도록
  `KPHAG_AM`/`KPHAG_N`/`KG`를 조정
- 세프트리악손 2 g q24h: Cmax ~150 mg/L, trough ~10 mg/L, fu ~0.10, ELF/혈장 유리분율 ~0.4
  → MIC 0.25 mg/L에 대해 fT>MIC ~100%
- 레보플록사신 750 mg qd: AUC24 ~100 mg·h/L, fu 0.70, ELF/혈장 ~1.2 → fAUC/MIC ~70
- 아지트로마이신: 혈청 Cmax ~0.4 mg/L인데 ELF 1-3 mg/L (심부 폐 구획 분리 이유)
- Halm 1998: 유효 치료 시 임상안정까지 중앙값 ~3일
- CAPE-COD (Dequin 2023): 중증 CAP에서 28일 사망 6.2% vs 11.9%

---

## ⚠️ 면책 조항

본 모델은 **교육 및 연구 목적의 정성적·반정량적 QSP 모델**입니다. 공개 문헌에 근거해
구성되었으나 환자 수준 데이터에 적합·검증되지 않았으며, **임상 의사결정·처방·규제 제출에
사용해서는 안 됩니다.**
