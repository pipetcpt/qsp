# Alport Syndrome (알포트 증후군, AS) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking COL4A3/A4/A5
> genotype-dependent type-IV collagen network failure (defective α-chain
> folding/ER retention, absent α3α4α5(IV) GBM network, mechanically inferior
> compensatory α1α2(IV) network) to glomerular basement membrane (GBM)
> structural progression (thinning → lamellation/basket-weave splitting →
> segmental/global glomerulosclerosis), podocyte foot-process effacement and
> apoptosis, proteinuria/hematuria, compensatory hyperfiltration and RAAS/
> endothelin-1-driven glomerular hypertension, and a TGF-β1/CTGF/miR-21
> fibrotic cascade culminating in progressive eGFR decline and ESRD — coupled
> to cochlear (progressive sensorineural hearing loss) and ocular (anterior
> lenticonus, dot-and-fleck retinopathy) basement-membrane phenotypes, and to
> RAAS blockade (ramipril/losartan), sparsentan (dual ETA/AT1 antagonist),
> bardoxolone methyl (Nrf2 activator), lademirsen (anti-miR-21 antisense
> oligonucleotide), dapagliflozin (SGLT2 inhibitor), and finerenone
> (nonsteroidal mineralocorticoid receptor antagonist) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`alp_qsp_model.dot`](alp_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`alp_qsp_model.svg`](alp_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`alp_qsp_model.png`](alp_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`alp_mrgsolve_model.R`](alp_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`alp_shiny_app.R`](alp_shiny_app.R) |
| 📚 References             | [`alp_references.md`](alp_references.md) |

---

## 1. 질환 개요 (Disease in one paragraph)

알포트 증후군은 IV형 콜라겐 α3, α4, α5 사슬을 코딩하는 **COL4A5**(Xq22.3,
X연관, 전체의 ~80%), **COL4A3/COL4A4**(2q36.3, 상염색체) 유전자의 변이로
발생하는 유전성 사구체기저막(GBM) 질환이다. 결함이 있는 α사슬은 소포체에서
삼중나선 형성에 실패해 ERAD로 분해되며, 그 결과 성숙한 GBM에서 정상적인
α3α4α5(IV) network가 형성되지 못하고 태아형 α1α1α2(IV) network가 잔류한다.
이 대체 네트워크는 기계적 강도가 낮아 GBM이 점진적으로 얇아지고(thin GBM) →
층판화("basket-weave" 분리) → 미세 파열(microtear)로 진행하며, 이는 족세포
발돌기 소실과 사구체경화로 이어진다. 사구체 손실은 대상성 과여과와 사구체
내압 상승을 유발하고, 이는 RAAS(안지오텐신 II)와 엔도텔린-1(ET-1) 경로를
활성화시켜 메산지움세포 활성화 → TGF-β1/CTGF/miR-21 매개 섬유화 캐스케이드를
증폭시키는 악순환을 형성한다. 유전형에 따라 중증도가 크게 다르며, **X연관
남성**(반접합, 가장 중증, 평균 20-30대 말기신부전)과 **상염색체 열성**(양쪽
대립유전자 변이) 환자는 가장 빠르게 진행하고, **X연관 여성**(모자이크,
X-비활성화에 따라 다양)과 **상염색체 우성 이형접합**(thin basement membrane
nephropathy 스펙트럼)은 상대적으로 경미하다. GBM과 동일한 IV형 콜라겐이
와우(기저막/혈관조)와 수정체 전낭에도 존재하여, 진행성 고음역 감각신경성
난청과 전방수정체원추(anterior lenticonus)·점상-반점상 망막병증 등 신장외
표현형을 동반하는 것이 특징이다. 표준 치료는 **RAAS 차단제**(ACE
억제제/ARB)로 조기 투여 시 말기신부전 도달을 10년 이상 지연시키며, 최근에는
이중 ETA/AT1 길항제 **sparsentan**, Nrf2 활성화제 **bardoxolone methyl**,
항-miR-21 항센스올리고 **lademirsen**, SGLT2 억제제, 비스테로이드성
미네랄로코르티코이드 수용체 길항제(**finerenone**) 등 표적치료가 임상에서
평가되고 있다.

## 2. 기계론적 지도 클러스터 (12개 클러스터, 118개 노드)

1. 유전적 병인 (COL4A5 X연관, COL4A3/COL4A4 상염색체 — XLAS 남성/여성,
   ARAS, ADAS/thin-BM 스펙트럼, digenic modifier, truncating vs missense
   변이군, 가족 선별검사·유전상담)
2. 분자병태생리: IV형 콜라겐 네트워크 (α3α4α5(IV) 삼중나선 형성 실패 → ERAD
   → 태아형 α1α2(IV) 잔류 → NC1 sulfilimine crosslink 소실 → GBM 기계적
   강도 저하, 렌즈낭/와우 기저막 조립 결함, 이식 후 항-GBM 알로면역 위험)
3. GBM 구조 진행 (Thin GBM → basket-weave 층판화 → 불규칙 비후/박화 반복 →
   진행성 splitting → 음전하장벽 소실 → 미세파열 → 분절/전체 사구체경화 →
   네프론 손실)
4. 족세포 손상 및 단백뇨 (기계적 스트레스 전달 → 발돌기 소실 → nephrin/
   podocin 슬릿막 붕괴 → 액틴 재배열 → 세포자멸/탈락 → 알부민 누출 →
   미세알부민뇨 → 현성 단백뇨 → 신증후군범위 단백뇨, 지속성 혈뇨)
5. 사구체 혈류역학 및 RAAS (네프론손실 대상성 과여과 → 단일신원GFR 증가 →
   사구체모세관고혈압 → 유입/유출세동맥 긴장도 변화 → RAAS 활성화 →
   안지오텐신II/AT1 신호 → 알도스테론 섬유화 증폭, ET-1/ETA 매개 혈관수축·
   메산지움 증식, 사구체내압 악순환 피드백, 이차성 고혈압)
6. 섬유화 캐스케이드 (메산지움세포 활성화 → 매트릭스 확장, TGF-β1/CTGF
   유도, 세뇨관상피 EMT/EndMT, 간질 섬유모세포 증식·콜라겐 침착, 세뇨관
   위축, 모세혈관 희박화, IFTA 복합지표, miR-21 증폭 루프)
7. 염증 증폭 (대식세포 침윤·M1/M2 분극, 보체 대체경로 활성화, 세뇨관 NLRP3
   인플라마솜, MCP-1/CCL2 화학주성, T세포 침윤, 활성산소종 증폭)
8. 이신장외 병변 — 와우 (기저막/혈관조 콜라겐IV 결함 → 코르티기관 유모세포
   취약성 → 진행성 고음역 감각신경성 난청 → 어음영역 침범, 전정 침범(희귀))
9. 이신장외 병변 — 안구 (전방수정체낭 콜라겐IV 결함 → 전방수정체원추 →
   낭 미세파열("oil-droplet" 반사) → 점상-반점상 망막병증, 후다형성각막
   이상증(희귀), 재발성각막상피미란, 굴절변화/백내장 위험)
10. 약물 약동학 (Ramipril/Ramiprilat, Losartan/E-3174, Sparsentan,
    Bardoxolone methyl, Lademirsen SC/혈장/신장조직 3구획, Dapagliflozin,
    Finerenone)
11. 약물 작용기전 (ACE억제·ARB → AngII/AT1 억제 → 병용 RAAS 차단 → 사구체
    내압 저하, Sparsentan 이중 ETA/AT1 차단, Bardoxolone Nrf2/Keap1 →
    항산화·항섬유화 및 급성 eGFR 상승 안전성 신호, Lademirsen anti-miR-21
    → CTGF 억제, SGLT2 억제 → 관구사구체피드백 회복, Finerenone MR차단)
12. 임상 엔드포인트 (eGFR 궤적/감소율, UACR, 40% eGFR감소 도달시간, ESRD,
    투석/이식, 이식 후 항-GBM 알로면역신염 위험, 청력역치 진행, 안구병기
    복합지표, 신장/전체 생존)

## 3. mrgsolve 모델 (24 ODE 구획)

* **약물 PK (15개 구획)** — Ramipril 위장관데포/Ramiprilat 혈장(2), Losartan
  데포/E-3174 활성대사체 혈장(2), Sparsentan 데포/혈장(2), Bardoxolone
  methyl 데포/혈장(2), Lademirsen 피하데포/혈장/신장조직(3), Dapagliflozin
  데포/혈장(2), Finerenone 데포/혈장(2).
* **질환/PD (7개 구획)** — GBM 구조건전성지표(GBM_INTEG), 생존족세포분율
  (PODO_FRAC), 잔존 기능 네프론분율(NEPHRON_FRAC), 섬유화/IFTA지표
  (FIBROSIS), miR-21 상대활성도(MIR21), UACR, eGFR.
* **신장외 임상 엔드포인트 (2개 구획, 지연 이완모델)** — 청력역치변화
  (HEARING_LOSS), 안구중증도지표(OCULAR_SCORE).
* GBM 손상은 유전형별 SEVERITY 배율과 사구체내압 악순환 피드백(vicious
  feedback)에 의해 가속되며, 이는 족세포 소실률에 직접 반영된다. 네프론
  손실이 진행할수록 단일신원GFR과 사구체내압이 상승해 RAAS/ET-1 축을
  자극하고, 이는 다시 TGF-β1/CTGF 섬유화 구동력을 증폭시키는 양성 피드백
  루프로 구현했다. RAAS 차단제(ACEi/ARB)와 sparsentan의 AT1 차단 성분은
  경쟁적 억제로 통합(at1_block_total)되며, sparsentan의 ETA 차단과 직접
  항섬유화 성분, bardoxolone의 Nrf2 항산화/항섬유화 성분, lademirsen의
  anti-miR-21 성분, finerenone의 MR차단 항섬유화 성분이 각각 독립적인
  Emax 항으로 fibrotic_drive를 저감시킨다. SGLT2 억제제는 관구사구체
  피드백 회복을 통해 단일신원GFR 상승을 직접 저감한다. Bardoxolone은
  급성 크레아티닌 비의존적 eGFR 상승(CARDINAL 임상의 안전성 신호)을 별도
  항으로 반영했다.

### 10개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 자연경과 - XLAS 남성 (SEVERITY=1.0) | Jais 2000 J Am Soc Nephrol |
| 2 | 자연경과 - ARAS (SEVERITY=1.15) | Storey 2013 J Am Soc Nephrol |
| 3 | 자연경과 - ADAS/이형접합 thin-BM (SEVERITY=0.4) | Kamiyoshi 2016 Clin J Am Soc Nephrol |
| 4 | Ramipril, 단백뇨 발생 후 시작 | Gross 2012 Kidney Int |
| 5 | Ramipril 조기(무증상기) 시작 | Gross 2020 Kidney Int (EARLY PRO-TECT Alport) |
| 6 | Losartan (ARB) 단독 | Temme 2012 Kidney Int (RAAS 억제 외삽) |
| 7 | Sparsentan (dual ETA/AT1) 단독 | Komers 2023 Lancet (DUPLEX FSGS 기전 외삽) |
| 8 | Bardoxolone methyl | Chertow 2021 Am J Nephrol (CARDINAL) |
| 9 | Lademirsen (RG-012, anti-miR-21) | Gomez 2015 J Clin Invest (HERA 감쇄반영) |
| 10 | 병용: RAAS 최대차단 + Dapagliflozin + Sparsentan | Heerspink 2020 NEJM 외삽 |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — 유전형(중증도) 선택, 연령, 시뮬레이션 기간.
2. **PK** — 약물별(ACEi/ARB/Sparsentan/Bardoxolone/Lademirsen/SGLT2i/
   Finerenone) 혈장/조직 농도.
3. **PD 주요지표 (신장)** — GBM 구조건전성·족세포분율, 섬유화/혈류역학.
4. **임상 엔드포인트** — eGFR 궤적(ESRD 임계선 포함), UACR, ESRD 도달시간.
5. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
6. **바이오마커** — UACR, GBM 구조지표, miR-21 활성도, IFTA 지표.
7. **이신장외 병변 (와우·안구)** — 청력역치 진행, 안구중증도지표.
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg alp_qsp_model.dot -o alp_qsp_model.svg
dot -Tpng -Gdpi=150 alp_qsp_model.dot -o alp_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("alp_mrgsolve_model.R") %>% param(SEVERITY = 1.0, IS_MALE_XL = 1)
e_ram <- ev(amt = 5, cmt = "RAM_GUT", time = 0, ii = 24, addl = 365*15-1)
out <- mod %>% ev(e_ram) %>% mrgsim(end = 24*365*15, delta = 24)  # 15년 추적
plot(out, c("GBM_INTEG", "PODO_FRAC", "EGFR", "UACR"))

# 3) Shiny 대시보드 실행
shiny::runApp("alp_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| XLAS 남성 자연경과 ESRD 도달연령 (~20대 후반) | 195가족 코호트 | Jais 2000 J Am Soc Nephrol |
| Ramipril 조기투여 ESRD 지연 (~13년) | 후향적 코호트 | Gross 2012 Kidney Int |
| Ramipril 무증상기 투여 알부민뇨/GBM 이득 | 무작위 대조 3상 | Gross 2020 Kidney Int (EARLY PRO-TECT Alport) |
| Sparsentan 단백뇨 감소 (FSGS 기전 외삽) | 무작위 대조 3상 | Komers 2023 Lancet (DUPLEX) |
| Bardoxolone methyl 급성 eGFR 상승 및 만성기 신호 | Alport 특이적 3상 | Chertow 2021 Am J Nephrol (CARDINAL) |
| Lademirsen anti-miR-21 항섬유화 (전임상 근거) | 마우스 모델 | Gomez 2015 J Clin Invest |
| Dapagliflozin 단백뇨/과여과 감소 (비당뇨 CKD 외삽) | 무작위 대조 3상 | Heerspink 2020 NEJM (DAPA-CKD) |
| 진행성 감각신경성 난청 자연경과 | 측두골 병리 | Merchant 2004 Laryngoscope |

## 7. 모델 검증 상태

이 컨테이너에는 초기 R/mrgsolve 실행환경이 설치되어 있지 않아(`Rscript`
부재), mrgsolve 모델은 **문헌 기반 파라미터 설계 및 코드 자체검토(구획/
파라미터 일치성 — 25개 컴파트먼트 전부 `$ODE`/`$INIT`와 1:1 대응 확인,
괄호·중괄호 균형 검사)** 단계까지 완료되었으며 실제 컴파일·적분 실행으로
수치를 검증하지는 못했다. `.dot` 파일은 `apt-get install graphviz`로 설치한
Graphviz `dot`으로 실제 렌더링해 SVG/PNG를 생성·확인했다(118 노드, 12
클러스터). 40편의 PubMed 문헌은 실제 PubMed 페이지의 저자/연도/저널 정보와
대조해 PMID를 개별 검증했다. mrgsolve/R 환경이 있는 곳에서 위 "실행 방법"
대로 실행해 수치 적분 결과를 확인할 것을 권장한다.
