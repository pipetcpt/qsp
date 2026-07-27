# Clostridioides difficile Infection (CDI) — QSP Model

**클로스트리디오이데스 디피실 감염 — 정량적 시스템 약리학 모델**

<a href="cdi_qsp_model.svg"><img src="cdi_qsp_model.png" width="100%" alt="CDI QSP mechanistic map"></a>

*(그림 클릭 → 확대 가능한 SVG / click the map for the zoomable SVG)*

---

## 1. 이 모델의 출발점 (The premise)

CDI는 감염질환 중에서 **병원체를 죽이는 것이 가장 쉬운 부분**인 드문 경우입니다.
경구 반코마이신은 며칠 안에 대변에서 영양형 *C. difficile*을 사실상 멸균하고 약
80%의 환자를 치료합니다 — 그리고 그중 **네 명 중 한 명이 재발**합니다. 병원체를
제거한 그 약이, 약이 닿지 않는 포자 저장고가 기다리는 동안 정착 저항성
(colonization resistance)을 담당하는 미생물 길드를 바닥에 눌러 놓았기 때문입니다.

지난 15년의 모든 치료적 진보 — 피닥소마이신, 베즐로톡수맙, FMT, SER-109,
RBX2660 — 는 살균 속도가 아니라 **바로 이 간극**을 겨냥한 것입니다.

CDI is unusual among infections in that killing the organism is the easy part.
Oral vancomycin sterilises the stool within days and cures ~80% of episodes —
and then one in four patients relapses, because the same drug that cleared the
organism also held down the microbial guild that provides colonization
resistance while a drug-proof spore reservoir waited. Every therapeutic advance
of the last 15 years targets that gap rather than the kill rate.

**그래서 이 모델은 병원체가 아니라 재발 루프(recurrence loop)를 중심으로
설계되었습니다.**

```
항생제 노출
  → 절대 혐기성 길드 붕괴 (특히 bai⁺ 7α-dehydroxylating Clostridia)
  → 2차 담즙산(DCA/LCA) 소실 + 결합형 1차 담즙산(taurocholate) 상승
     + 시알산·Stickland 아미노산 영양 니치 개방
  → CspC 매개 포자 발아 + 무제한 영양형 증식   (= 정착 저항성 상실)
  → PaLoc 탈억제 (CodY/CcpA off, SigD/Spo0A on)
  → TcdA/TcdB가 RhoA/Rac1/Cdc42를 글루코실화
  → 밀착연접 붕괴 · 대장상피 사멸 · FZD 차단에 의한 줄기세포 재생 실패
     · pyrin 인플라마좀 IL-1β · 호중구 침윤 · 위막 형성
  → 설사 / 백혈구증가 / 저알부민혈증 / 신손상
  → 그리고 항생제가 끝난 뒤에도 포자 저장고는 살아남고 미생물군집은 아직
     회복되지 않았으므로 → 재발
```

---

## 2. 산출물 (Deliverables)

| 파일 | 내용 |
|------|------|
| [`cdi_qsp_model.dot`](cdi_qsp_model.dot) | Graphviz 기계론적 지도 소스 — **20개 모듈 · 242 노드 · 378 엣지** |
| [`cdi_qsp_model.svg`](cdi_qsp_model.svg) / [`.png`](cdi_qsp_model.png) | 렌더링된 지도 (PNG는 150 dpi) |
| [`cdi_mrgsolve_model.R`](cdi_mrgsolve_model.R) | **61 ODE 구획 · 240 파라미터 · 18개 치료 시나리오** mrgsolve 모델 + 검증 리포트 |
| [`cdi_shiny_app.R`](cdi_shiny_app.R) | **9개 탭** 인터랙티브 대시보드 (regimen builder 포함) |
| [`cdi_references.md`](cdi_references.md) | **129편** — 모든 PMID를 PubMed E-utilities로 조회·제목 대조 검증 |

재현:

```bash
dot -Tsvg cdi_qsp_model.dot -o cdi_qsp_model.svg
dot -Tpng -Gdpi=150 cdi_qsp_model.dot -o cdi_qsp_model.png
Rscript cdi_mrgsolve_model.R        # 18개 시나리오 실행 + 검증 리포트 출력
Rscript -e 'shiny::runApp("cdi_shiny_app.R")'
```

---

## 3. 모델 구조 (Model structure)

세 개의 층이 동시에 돌아갑니다.

### (1) 생태 층 (Ecology) — 6개 길드

| 상태 | 길드 | 건강 시 상대풍부도 | 모델 내 역할 |
|---|---|---|---|
| `MB_SBA` | bai⁺ 7α-dehydroxylating Clostridia (*C. scindens* 등) | 0.10 | **2차 담즙산의 유일한 공급원 — 재발의 결정 인자** |
| `MB_BUT` | butyrogenic Lachnospiraceae / Ruminococcaceae | 0.30 | 부티르산 생산, 장벽 지지, 독소 발현 억제 |
| `MB_BAC` | Bacteroidetes | 0.42 | BSH 탈결합 + **시알산을 유리시키는 공급원** |
| `MB_BIF` | Bifidobacterium / Actinobacteria | 0.08 | BSH, 아미노산 경쟁 |
| `MB_ENT` | Enterobacteriaceae | 0.004 | 니치가 비면 폭발적 증식 → LPS |
| `MB_ENC` | Enterococcus | 0.002 | VRE 선택, *C. difficile* 교차 급이 |

이 층에서 가장 중요한 구조적 선택: **영양 경쟁 길드 집합에서 Bacteroidetes를
제외**했습니다. Bacteroidetes는 뮤신 시알리다아제로 시알산을 *유리시키는* 쪽이고
소비하는 쪽이 아니기 때문입니다 (Ng 2013). 이것이 **Bacteroides를 살려두는
반코마이신이 그럼에도 니치를 활짝 열어 두는 이유**이며, 모델에서 반코마이신과
피닥소마이신 팔이 갈라지는 기계론적 지점입니다.

길드 회복 동역학도 in vitro 배가시간이 아니라 **관찰된 항생제 후 회복 시간
(수 주~수 개월)** 에 맞췄고, 여기에 0차 재정착(re-seeding) 플럭스를 더했습니다.
그래야 회복 시간이 nadir 깊이에 로그 선형으로 좌우되지 않고, 대신 **약제별 nadir
깊이가 kill/재정착 균형으로 결정**됩니다. 재정착의 브레이크는 총 항균 활성이
아니라 **길드별 kill rate** 입니다 — 협범위 대 광범위 논쟁이 한 줄의 수식으로
표현되는 지점입니다.

### (2) 담즙산·영양 니치 (The germination switch)

$MAIN에서 **건강한 정상상태를 파라미터로부터 대수적으로 풀어냅니다** —
결합형 풀 → BSH 탈결합 → bai 7α-dehydroxylation 3단계 캐스케이드를,
cholate와 chenodeoxycholate의 Michaelis–Menten 균형에 대한 **닫힌 형태의 근**까지
포함해서. 따라서 `K7A`나 BSH 파라미터를 바꾸면 baseline이 조용히 흘러가는 대신
같이 움직입니다.

- `BA_TCA` (taurocholate 등 결합형 1차) → CspC 수용체 작용제 → **발아**
- `BA_CDCA` (유리 chenodeoxycholate) → 경쟁적 **발아 억제제**
- `BA_DCA` / `BA_LCA` (2차) → **영양형 증식·포자형성·독소 발현 억제**
- `NUT_SIA` (유리 시알산 + succinate), `NUT_AA` (proline/glycine/leucine),
  `SCFA_BUT` (부티르산)

### (3) 병원체 · 숙주 (Pathogen and host)

- 4개 *C. difficile* 구획: 관내 포자 · 관내 영양형 · 점막 부착 · **점막/biofilm
  포자 저장고**(= 재발의 씨앗)
- PaLoc 조절: CodY/CcpA 영양 억제 때문에 **독소가 균 정점보다 늦게 올라갑니다**
  (집단이 스스로의 니치를 소모하면서 유도됨) — in vivo 관찰과 일치
- TcdB → CSPG4/FZD/Nectin-3 → 자가절단 → Rho GTPase 글루코실화 → 액틴 붕괴
  → 밀착연접 실패 + **FZD 차단으로 Lgr5⁺ 줄기세포 재생 자체가 표적이 됨**
- pyrin/NLRP3 → IL-1β → IL-8 → 호중구 → 위막; IL-22는 보호적 limb
- 임상 판독: `STOOL`(무형변/일), `WBC`, `CRE`, `ALB` + IDSA/SHEA 중증도 플래그

### 구획 요약 (61 ODEs)

| 그룹 | 개수 | 상태 |
|---|---|---|
| 미생물 길드 | 6 | `MB_SBA` `MB_BUT` `MB_BAC` `MB_BIF` `MB_ENT` `MB_ENC` |
| 담즙산 | 5 | `BA_TCA` `BA_CA` `BA_CDCA` `BA_DCA` `BA_LCA` |
| 영양소 · SCFA | 3 | `NUT_SIA` `NUT_AA` `SCFA_BUT` |
| *C. difficile* | 4 | `CD_SPORE_L` `CD_VEG` `CD_MUC` `CD_SPORE_B` |
| 독소 | 6 | `TCDA` `TCDB` `TOX_CPLX` `TCDA_MUC` `TCDB_MUC` `CDT` |
| 상피 · 장벽 | 5 | `EPI` `EPI_SC` `EPI_TJ` `EPI_MUCUS` `EPI_PERM` |
| 면역 | 7 | `IM_IL8` `IM_IL1B` `IM_TNF` `IM_NEUT` `IM_IL22` `IM_PSM` `AB_IGG` |
| 임상 | 5 | `EPI_H2O` `STOOL` `WBC` `ALB` `CRE` |
| 약물 PK | 20 | 반코마이신 3 · 피닥소마이신+OP-1118 4 · 메트로니다졸 3 · 베즐로톡수맙 3 · 리팍시민 2 · 리디닐라졸 2 · 유발 항생제 2 · 생균치료제 1 |

---

## 4. 치료 시나리오 (18)

| # | 시나리오 |
|---|---|
| S01 | 건강 대조군 (self-calibrated baseline 검증) |
| S02 | 항생제만 — 병원체 없는 허용 상태 |
| S03 | 미치료 CDI (자연 경과) |
| S04 | 무증상 보균 (기존 항독소 IgG 높음) |
| S05 | 메트로니다졸 500 mg q8h × 10일 |
| S06 | 반코마이신 125 mg q6h × 10일 (표준 1차) |
| S07 | 피닥소마이신 200 mg q12h × 10일 |
| S08 | 피닥소마이신 extended-pulsed (EXTEND) |
| S09 | 반코마이신 + 베즐로톡수맙 (MODIFY 전략) |
| S10 | 피닥소마이신 + 베즐로톡수맙 |
| S11 | 반코마이신 → FMT |
| S12 | 반코마이신 → SER-109 (경구 Firmicutes 포자 ×3일) |
| S13 | 반코마이신 → RBX2660 (직장 현탁액) |
| S14 | 리디닐라졸 200 mg q12h × 10일 |
| S15 | 반코마이신 + 리팍시민 chaser |
| S16 | 반코마이신 taper/pulse (IDSA 다회 재발 요법) |
| S17 | 리보타입 027 · 면역저하 · 5일 지연 치료 (fulminant) |
| S18 | 유발 항생제를 끊지 않은 최악의 경우 |

---

## 5. 검증 결과 (Validation)

### Baseline은 가정이 아니라 풀린 고정점입니다

무약물 90일 시뮬레이션에서 **최대 표류 0.000%**:

| 상태 | Day 0 | Day 90 | 표류 |
|---|---|---|---|
| `MB_SBA` | 0.100 | 0.100 | 5×10⁻⁵ % |
| `BA_TCA` | 12.0 µM | 12.0 µM | 0 % |
| `BA_DCA` | 449.8 µM | 449.8 µM | 1×10⁻⁵ % |
| `BA_LCA` | 300.0 µM | 300.0 µM | 9×10⁻⁶ % |
| `NUT_SIA` | 0.050 mM | 0.050 mM | −4×10⁻⁵ % |
| `SCFA_BUT` | 15.0 mM | 15.0 mM | 4×10⁻⁵ % |
| `STOOL` / `WBC` / `ALB` / `CRE` | 0.8 / 7.0 / 4.20 / 0.90 | 동일 | 0 % |

### 약물 노출은 문헌값을 재현합니다

| 약제 | 모델 | 보고값 |
|---|---|---|
| 반코마이신 (대변) | 929 µg/g | 500–3000 µg/g |
| 피닥소마이신 + OP-1118 (대변) | 630 + 288 = 918 µg/g | ~1000–1400 µg/g (Sears 2012) |
| 메트로니다졸 (대변) | 10.4 µg/g | 9.3 µg/g 수양성 변 (Bolton 1986) |
| 메트로니다졸 (혈장) | 14.3 mg/L | trough ~10, Cmax ~25 mg/L |
| 리디닐라졸 (대변) | 863 µg/g | ~1000 µg/g |
| 베즐로톡수맙 (혈장) | Cmax ~230 mg/L, t½ 19일 | Yee 2019 popPK |

메트로니다졸은 **염증 점막을 통한 분비로만 관내에 도달**하도록 모델링했으므로,
환자가 좋아지면 대변 농도가 함께 떨어집니다 (Bolton 1986: 수양성 변 9.3 →
정상 변 1.2 µg/g). 용량이 아니라 점막 염증이 노출을 결정합니다.

### 급성기 및 반응

| 시나리오 | 정점 무형변/일 | 정점 WBC | 최대 Cr | 최저 Alb | TTROD | 모델 내 재발 |
|---|---|---|---|---|---|---|
| S03 미치료 | 7.9 | 14.4 | 1.52 | 3.81 | 12일 | — (자연 소실) |
| S05 메트로니다졸 | 4.9 | 12.9 | 1.45 | 3.85 | **4일** | **22일** |
| S06 반코마이신 | 4.4 | 15.3 | 1.55 | 3.61 | **4일** | **19일** |
| S07 피닥소마이신 | 4.0 | 11.0 | 1.26 | 4.12 | **3일** | 없음 |
| S11 반코 → FMT | 4.4 | 11.3 | 1.29 | 3.96 | **4일** | 중단 후 1–2일 flicker |
| S14 리디닐라졸 | 4.2 | 11.1 | 1.27 | 4.12 | **4일** | 없음 |
| S17 027 fulminant | 9.5 | 16.4 | 1.58 | 3.63 | 10일 | 33일 |
| S18 항생제 지속 | 4.7 | 16.2 | 1.58 | 3.03 | 4일 | 21일 |

치료군 TTROD 3–4일 (임상시험 중앙값 2–4일). 미치료군은 12일에 정점을 지나
군집이 회복되면서 자연 소실 — 관찰된 자연 경과와 일치.

### 치료 종료 시점의 생태 상태가 다음을 결정합니다

| 시나리오 | `MB_SBA` (정상 대비) | `BA_DCA` (µM) | 저장고 (log CFU/g) | RRI |
|---|---|---|---|---|
| S06 반코마이신 | **2.8 %** | **11** | 6.99 | 0.560 |
| S05 메트로니다졸 | 9.1 % | 57 | 6.56 | 0.397 |
| S07 피닥소마이신 | **8.3 %** | **68** | 5.28 | 0.250 |
| S08 피닥소 EXTEND | 12.5 % | 95 | 3.83 | 0.017 |
| S14 리디닐라졸 | 19.7 % | 156 | 5.34 | 0.162 |
| S11 반코 → FMT | **39.4 %** | **180** | 6.50 | 0.159 |
| S12 반코 → SER-109 | **106.8 %** | **424** | 6.49 | 0.000 |

### 재발률 매핑은 6개 3상 앵커에 대해 R² ≈ 0.73

결정론적 궤적은 재발하거나 하지 않을 뿐이지만, 실제 코호트는 갈라집니다.
그래서 기계론적 **재발 위험 지수(RRI)** — 살아남은 포자 저장고 × 회복되지 않은
복원 길드 ÷ 이용 가능한 중화 항체 — 를 관찰된 8주 재발률에 로지스틱으로
명시적으로(그리고 감사 가능하게) 매핑합니다:

```
logit(p) = −2.24 + 1.99 × RRI        (n = 6 앵커, R² = 0.725)
```

| 시나리오 | RRI | 관찰 8주 재발 | 예측 | 출처 |
|---|---|---|---|---|
| 메트로니다졸 | 0.397 | 23.0 % | 19.1 % | Johnson 2014 CID |
| 반코마이신 | 0.560 | 25.3 % | 24.6 % | Louie 2011 / Cornely 2012 |
| 피닥소마이신 | 0.250 | 15.4 % | 15.0 % | Louie 2011 / Cornely 2012 |
| 반코 + 베즐로톡수맙 | 0.386 | 16.5 % | 18.7 % | Wilcox 2017 MODIFY I/II |
| 반코 → FMT | 0.159 | 9.0 % | 12.8 % | van Nood 2013 / Kelly 2016 |
| 반코 → SER-109 | 0.000 | 12.0 % | 9.7 % | Feuerstadt 2022 ECOSPOR III |

순서가 옳고, 예측이 관찰 범위(9–25 %) 안에 들어옵니다. RBX2660(PUNCH CD3)은
대조군 재발률 42.5 %가 Louie 2011과 동일한 모집단이 아니므로 **의도적으로 앵커에서
제외**했습니다.

---

## 6. 모델이 말하는 것 (What the model says)

두 개의 치료적 양(quantity)이 구조적으로 분리되어 있습니다:
**`KILLCD`**(병원체가 얼마나 빨리 죽는가)와 **`MB_SBA`에 대한 부수적 kill**
(복원 길드가 얼마나 심하게 맞는가). 첫 번째만 잘하고 두 번째에서 실패하는 약제는
**에피소드를 치료하고 재발을 사들입니다.**

- **피닥소마이신 · 리디닐라졸** — 둘 다 잘함 → 모델 궤적에서 재발하지 않음
- **반코마이신 · 메트로니다졸** — 치료하고 재발 (19–22일)
- **FMT · SER-109 · RBX2660** — 길드를 다시 심어 재발을 중단시킴
- **베즐로톡수맙** — 어느 쪽도 하지 않지만 취약 구간 동안 독소를 중화
- **유발 항생제를 끊지 않는 것** — 가장 나쁨 (중증 기준 충족 32일)

`S02`(항생제만, 병원체 없음)와 `S04`(같은 노출, 항독소 IgG 3.2)가 왜 이 라이브러리
안에 있는지도 여기서 설명됩니다: 허용 상태와 병원체와 숙주 항체는 서로 다른
축이고, 무증상 보균은 그 세 번째 축의 결과입니다 (Kyne 2000).

---

## 7. 한계 (Limitations)

- **준정량적(semi-quantitative)이며 결정론적.** 파라미터는 문헌에 근거한
  자릿수 수준 추정치이고, 적합된 집단 모델이 아닙니다. 팔 수준 재발률은 단일
  실행에서 재발을 세는 대신 명시적 RRI→확률 매핑에서 나옵니다.
- **FMT / SER-109 앵커는 재발성 CDI 모집단**에서 나왔으므로 1차 에피소드 팔보다
  기저 위험이 높습니다 → 이 두 앵커는 보수적입니다.
- **독소 농도 스케일**은 자릿수 수준입니다. 발표된 분석법이 단위와 (세포독성
  역가 대 면역측정 질량) 측정 대상이 다릅니다. 보정된 특징은 절대값이 아니라
  *상대적* 동역학(CodY/CcpA 영양 억제로 독소가 균 정점에 뒤처지는 것)입니다.
- **반코마이신 taper/pulse**(S16)는 모델이 잘 다루지 못합니다. 실제 taper의
  이점 중 일부는 투여 간격 사이에 포자가 발아하고 다음 pulse에 죽는 것인데,
  구획 모델은 이 시간적 미세구조를 온전히 담지 못합니다.
- **리보타입 027 파라미터**(`RT027` `RTTOXF` `RTSPOR`)는 BI/NAP1/027 계통에서
  보고된 집합적 표현형이며 특정 측정 균주가 아닙니다.

---

## 8. 면책 (Disclaimer)

교육 및 연구 목적의 QSP 모델입니다. 공개 문헌과 임상시험 데이터로 구성했으나
독립적으로 검증·인증되지 않았습니다. **임상 의사결정, 처방, 규제 제출에 직접
사용해서는 안 됩니다.**

Educational / research QSP model. Not independently validated. Not for clinical
decision-making, prescribing, or regulatory submission.

전체 인용 목록: [`cdi_references.md`](cdi_references.md) (129편, PubMed 검증).
