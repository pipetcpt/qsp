# Retinitis Pigmentosa (색소성 망막염, RP) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking RHO/RPGR/USH2A/
> PDE6/RPE65 genotype-dependent primary rod photoreceptor apoptosis
> (rhodopsin misfolding/ER stress, PDE6-cGMP-Ca2+ excitotoxicity, RPGR-
> ciliopathy transport failure) to secondary cone death (RdCVF loss,
> oxidative stress) and microglial/gliotic amplification — producing
> progressive night blindness, visual field constriction, and central
> vision loss — coupled to voretigene neparvovec (AAV2-RPE65 subretinal
> gene therapy), investigational RPGR gene augmentation (AAV8/AAV5-RPGR),
> MCO-010 optogenetic gene therapy (AAV2 multi-characteristic opsin,
> intravitreal), CNTF encapsulated-cell neuroprotection, N-acetylcysteine
> (antioxidant), and vitamin A palmitate/DHA supplementation PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`rp_qsp_model.dot`](rp_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`rp_qsp_model.svg`](rp_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`rp_qsp_model.png`](rp_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`rp_mrgsolve_model.R`](rp_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`rp_shiny_app.R`](rp_shiny_app.R) |
| 📚 References             | [`rp_references.md`](rp_references.md) |

---

## 1. 질환 개요 (Disease in one paragraph)

색소성 망막염은 90여 개 이상의 유전자 변이로 발생하는 가장 흔한 유전성 망막
이영양증으로, 간상세포(rod)의 원발성 사멸이 원추세포(cone)의 이차적 사멸로
이어져 실명에 이르는 진행성 질환이다. 상염색체 우성 **RHO**(로돕신 오접힘 →
소포체 스트레스/UPR), X연관 **RPGR**(섬모 수송 결함), 상염색체 열성
**USH2A**(Usher 증후군 동반 가능), **PDE6A/B**(cGMP 가수분해 실패 → CNG 채널
과활성 → Ca2+ 독성), **RPE65**(시각회로 이성질화효소 결핍, Leber 선천성
흑암시/조기발현중증망막이영양증) 등이 주요 원인이다. 간상세포가 소실되면
망막 산소소비가 줄어 원추세포가 상대적 고산소 환경에 노출되고, 간상세포유래
원추세포생존인자(RdCVF/NXNL1)가 고갈되어 산화스트레스 매개 이차 원추세포
사멸이 촉발된다. 미세아교세포 활성화와 뮐러세포 반응성 신경교증이 이 과정을
증폭시키며, 임상적으로는 야맹증(최초 증상) → 진행성 주변시야 소실(고리암점 →
터널시야) → 말기 중심시력 소실 순으로 진행한다. RPE65 이대립유전자 변이
환자에는 **voretigene neparvovec**(AAV2-RPE65 망막하 유전자치료, 2017년 3상
Russell Lancet 연구 기반 FDA 승인)가 시각회로를 회복시키며, X연관 RP에는
AAV8/AAV5-RPGR 유전자치료가 임상 개발 중이다(Cehajic-Kapetanovic 2020 Nat
Med, XIRIUS 3상). 광수용체 생존과 무관하게 작동하는 **MCO-010 광유전학
치료**(다중특성 옵신, 잔존 양극/신경절세포 발현)는 말기·유전형 무관 환자에
적용 가능하며, CNTF 서방출 임플란트(신경보호, paradox한 ERG 억제 동반),
N-아세틸시스테인(항산화), 비타민A palmitate+DHA(경도 진행지연, Berson 1993/
2004)가 보조/신경보호 요법으로 사용된다. 말기 환자는 Argus II 망막보철
(전기자극 우회경로, 2020년 상업적 단종)이나 저시력 보조기기의 대상이 된다.

## 2. 기계론적 지도 클러스터 (12개 클러스터, 136개 노드)

1. 유전적 병인 (RHO 우성 · RPGR/RP2 X연관 · USH2A/PDE6A/B/CRB1/NR2E3/RP1/
   EYS 열성 · PRPF 스플라이싱인자 우성 · BBS 섬모병증 · Usher I/II/III)
2. 분자병태생리: 광수용체 단백질 이상 (로돕신 오접힘/UPR, PDE6 기능소실 →
   cGMP 축적 → CNG채널 과활성 → Ca2+ 독성 → calpain/caspase, RPGR-ORF15/
   BBSome 섬모수송결함, PRPF 스플라이싱 이상, PARP1/AIF 세포사)
3. 간상세포-원추세포 상호작용 및 이차 사멸 (RdCVF/NXNL1 소실, 망막 산소소비
   감소 → 고산소 노출 → 산화스트레스 → 원추세포 사멸, 포도당/젖산 셔틀 붕괴,
   mTOR/자가포식 이상, 중심와 원추세포 상대적 보존)
4. 망막색소상피/시각회로 (RPE65 이성질화효소, LRAT, 11-cis-레티날 재생,
   아포옵신 축적/구성적 신호, RPE 식세포작용, 색소이주/골편양색소침착, 혈관
   협착, 시신경위축)
5. 신경염증/망막 리모델링 (뮐러세포 반응성 신경교증, 미세아교세포 활성화,
   보체 활성화, 양극/수평세포 수상돌기 위축, 이상 신경돌기 발아, 신경절세포
   상대적 보존)
6. 임상 표현형 및 병기 (야맹증 → 고리암점 → 터널시야 → 중심시력소실, Usher
   난청, Bardet-Biedl 비만/다지증/신장질환, 후낭하백내장, 유전형별 진행속도)
7. 임상 엔드포인트 및 바이오마커 (ERG 간상/원추세포 진폭, Goldmann 시야,
   BCVA, OCT EZ폭, 안저자가형광 고리, FST, MLMT, CST, NEI-VFQ-25)
8. 유전자치료 PK/PD (Voretigene neparvovec 망막하 AAV2-RPE65, Botaretigene
   sparoparvovec/Cotoretigene toliparvovec AAV8/AAV5-RPGR, CRISPR 편집)
9. 광유전학 치료 PK/PD (MCO-010 유리체강내 AAV2 다중특성옵신, 잔존
   양극/신경절세포 발현 우회경로)
10. 신경보호/항산화 치료 PK/PD (CNTF 캡슐화세포 임플란트, N-아세틸시스테인,
    비타민A palmitate+DHA, 재조합 RdCVF 투자적, 발프로산 역사적)
11. 합병증 관리 (낭포황반부종-도르졸라마이드/아세타졸아마이드, 후낭하백내장
    수술, 저시력 재활)
12. 망막 보철/보조기기 (Argus II 망막상 전극배열, 저시력 보조기기)

## 3. mrgsolve 모델 (23 ODE 구획)

* **약물/벡터 PK (10개 구획)** — Voretigene neparvovec(망막하 벡터게놈/RPE65
  발현 2구획), 투자적 RPGR 유전자치료(벡터/발현 2구획), MCO-010(벡터/옵신
  발현 2구획), CNTF 임플란트(조직농도 1구획), N-아세틸시스테인(경구 데포/
  혈장 2구획), 비타민A/DHA 저장고(1구획).
* **질환/PD (6개 구획)** — 간상세포 생존분율(ROD_FRAC), 원추세포 생존분율
  (CONE_FRAC), 산화스트레스지수(ROS), 미세아교세포 활성지수(MICROGLIA),
  망막신경절세포 생존분율(RGC_FRAC).
* **임상 엔드포인트 (7개 구획, 지연 이완모델)** — ERG 간상/원추세포 진폭,
  시야면적, BCVA, FST, MLMT, 낭포황반부종 중심망막두께(CME_CST).
* RPE65 유전형에서는 유효 시각회로 플럭스가 낮을수록 아포옵신 독성이 증가해
  간상세포 사멸률이 가중되며, voretigene neparvovec 발현이 이를 Emax
  모델로 회복시킨다. XLRP 유전형에서는 RPGR 유전자치료 발현이 섬모수송결함
  기여분을 경쟁적으로 교정한다. MCO-010은 광수용체 생존과 독립적으로 작동
  하나 신경절세포 생존분율에 비례하는 우회 신호를 제공한다. CNTF는 간상/
  원추세포 사멸을 완화하지만 동시에 가역적인 ERG 진폭 억제(Birch 2013 역설
  현상)를 유발하도록 모델링했다.

### 10개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 자연경과 - RHO-adRP (SEVERITY=1.0) | Berson 1985 Am J Ophthalmol |
| 2 | 자연경과 - RPGR-XLRP (SEVERITY=1.6) | Hartong 2006 Lancet |
| 3 | 자연경과 - RPE65-LCA/EOSRD (SEVERITY=2.5) | Cideciyan 2013 PNAS |
| 4 | Voretigene neparvovec (RPE65-이대립유전자) | Russell 2017 Lancet, Maguire 2019 Ophthalmology |
| 5 | 투자적 RPGR 유전자치료 (XLRP) | Cehajic-Kapetanovic 2020 Nat Med, XIRIUS 2024 |
| 6 | MCO-010 광유전학 치료 (말기, 유전형 무관) | Busskamp 2010 Science, Sahel 2021 Nat Med |
| 7 | CNTF 캡슐화세포 임플란트 | Sieving 2006 PNAS, Birch 2013 Am J Ophthalmol |
| 8 | N-아세틸시스테인 경구 | Campochiaro 2020 J Clin Invest |
| 9 | 비타민A palmitate + DHA 경구 | Berson 1993/2004 Arch Ophthalmol |
| 10 | Voretigene neparvovec + NAC 병용 | 병용 신경보호 가설 (외삽) |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — 유전형(중증도) 선택, 연령, 시뮬레이션 기간.
2. **PK** — 유전자치료 발현(RPE65/RPGR/MCO 옵신), CNTF 조직농도, NAC 혈장농도,
   비타민A 저장고.
3. **PD 주요지표** — 간상/원추세포 생존분율, 산화스트레스·미세아교세포 활성.
4. **임상 엔드포인트** — ERG, 시야면적/BCVA, FST/MLMT.
5. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
6. **바이오마커** — 신경절세포 생존분율, 낭포황반부종 중심망막두께(CAI 병용
   옵션 포함).
7. **유전자치료·광유전학** — 시각회로 회복 지표(RPE65), 광유전학 우회신호
   (MCO-010).
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg rp_qsp_model.dot -o rp_qsp_model.svg
dot -Tpng -Gdpi=150 rp_qsp_model.dot -o rp_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("rp_mrgsolve_model.R") %>% param(SEVERITY = 2.5, IS_RPE65 = 1)
e_gt <- ev(amt = 100, cmt = "GT65_VG", time = 0)  # voretigene neparvovec 1회 망막하 투여
out <- mod %>% ev(e_gt) %>% mrgsim(end = 24*365*10, delta = 24)  # 10년 추적
plot(out, c("ROD_FRAC", "CONE_FRAC", "ERG_ROD", "MLMT"))

# 3) Shiny 대시보드 실행
shiny::runApp("rp_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| 자연경과 ERG 진폭 감소율 (~연 15-20%) | 3년 추적 자연경과 | Berson 1985 Am J Ophthalmol |
| Voretigene neparvovec MLMT/FST 개선 | 1년 3상, 4년 지속성 | Russell 2017 Lancet, Maguire 2019 Ophthalmology |
| RPGR 유전자치료 용량의존 망막 위축 | 1/2상 고용량군 | Cehajic-Kapetanovic 2020 Nat Med |
| 광유전학 부분 시기능 회복 | 단일 환자 증례 | Sahel 2021 Nat Med |
| CNTF 역설적 ERG 억제 | 무작위 3상 | Birch 2013 Am J Ophthalmol |
| 비타민A 경구보충 ERG 감소 완화(~20%) | DBA 무작위시험 | Berson 1993 Arch Ophthalmol |
| N-아세틸시스테인 원추세포 기능 개선 | 1상 용량증량 | Campochiaro 2020 J Clin Invest |
| CAI(도르졸라마이드/아세타졸아마이드) 황반부종 반응률 | 개방표지/후향 연구 | Fishman 1989, Grover 1997 |

## 7. 모델 검증 상태

이 컨테이너에는 R/mrgsolve 실행환경이 설치되어 있지 않아(`Rscript` 부재),
mrgsolve 모델은 **문헌 기반 파라미터 설계, 코드 자체검토(구획/파라미터
일치성, 괄호·차원 검사)** 단계까지 완료되었으며 실제 컴파일·적분 실행으로
수치를 검증하지는 못했다. `.dot` 파일은 Graphviz `dot`(로컬 설치)으로
렌더링해 실제로 SVG/PNG를 생성·확인했다(136 노드, 12 클러스터). 45편의
PubMed 문헌은 모두 실제 PubMed 페이지 대조로 PMID를 개별 검증했다(MCO-010
RESTORE 임상시험은 이 문서 작성 시점 동료검토 논문이 없어 학술대회 초록
수준으로만 언급하고 별도 문헌으로 등재하지 않음). mrgsolve/R 환경이 있는
곳에서 위 "실행 방법"대로 실행해 수치 적분 결과를 확인할 것을 권장한다.
