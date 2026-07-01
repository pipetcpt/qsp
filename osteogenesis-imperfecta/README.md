# Osteogenesis Imperfecta (골형성부전증, OI) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking COL1A1/COL1A2
> (and recessive CRTAP/P3H1/PPIB/SERPINH1/FKBP10/WNT1) collagen-biosynthesis
> defects to osteoblast/osteocyte dysfunction, sclerostin-mediated Wnt
> suppression, RANKL/OPG-driven high bone turnover, and excess TGF-beta
> bioavailability — producing low bone mineral density and recurrent
> low-trauma fracture — coupled to bisphosphonate (pamidronate, zoledronic
> acid), denosumab (anti-RANKL), teriparatide (PTH1R anabolic, adult type I
> only), setrusumab (anti-sclerostin, investigational), and fresolimumab
> (anti-TGF-beta, investigational) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`oi_qsp_model.dot`](oi_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`oi_qsp_model.svg`](oi_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`oi_qsp_model.png`](oi_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`oi_mrgsolve_model.R`](oi_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`oi_shiny_app.R`](oi_shiny_app.R) |
| 📚 References             | [`oi_references.md`](oi_references.md) |

---

## 1. 질환 개요 (Disease in one paragraph)

골형성부전증은 제1형 콜라겐을 암호화하는 **COL1A1/COL1A2** 유전자의 상염색체
우성 돌연변이(대부분 글리신 치환)가 원인인 가장 흔한 유전성 골취약질환으로,
드물게는 콜라겐 3-수산화 복합체(CRTAP·P3H1·PPIB)나 샤페론(SERPINH1/HSP47,
FKBP10), WNT1 등의 상염색체 열성 변이로도 발생한다. 정량적(quantitative,
null allele) 결함은 비교적 경증인 type I을, 구조적(qualitative) 결함은
주산기 치명적인 type II부터 진행성 변형을 동반하는 type III, 중등도의
type IV까지 다양한 중증도를 유발한다(Sillence 분류). 잘못 접힌 콜라겐은
소포체 스트레스와 미접힘단백질반응(UPR)을 유발해 조골세포 사멸을 증가시키고,
붕괴된 골세포(osteocyte) 네트워크는 **스클레로스틴(SOST)** 과다분비로 Wnt/
베타카테닌 신호를 억제해 골형성을 더 저하시킨다. 동시에 **RANKL/OPG 비율
상승**으로 파골세포 분화가 촉진되고, 무질서한 기질에서 유리되는 **과잉
활성 TGF-베타**가 조골세포 성숙을 저해하며 파골세포형성을 상승적으로
자극한다. 그 결과 골형성-골흡수 커플링이 파괴되어 저골량·고회전율 골대사와
반복적인 저외상성 골절, 진행성 장골 만곡, 척추압박골절, 치아형성부전,
공막의 청회색 변화, 청력저하, 관절과이완 등 골외 증상이 동반된다. 표준
치료는 정맥용 비스포스포네이트(파미드로네이트 주기요법, 졸레드론산)로
파골세포 자멸을 유도해 골흡수를 억제하며, 중증 열성형에는 데노수맙이
시도되고, 성인 type I 한정으로 테리파라타이드가 사용되나 성장판이 열려있는
소아에서는 골육종 위험으로 금기이다. 항-스클레로스틴 항체 세트루수맙과
항-TGF-베타 항체 프레솔리무맙은 병인 표적치료로 임상 개발 중이다.

## 2. 기계론적 지도 클러스터 (15개 클러스터, 106개 노드)

1. 유전적 병인 (COL1A1/COL1A2 우성 · CRTAP/P3H1/PPIB 3-수산화복합체 ·
   SERPINH1/FKBP10 샤페론 · IFITM5/WNT1/SP7 열성형)
2. 콜라겐 생합성/분비 경로 (삼중나선 형성, P4H/LH 변형, HSP47, ER 스트레스/
   UPR, 프로콜라겐 절단, 원섬유 조립, LOX 가교)
3. 골기질/조골세포 기능이상 (기질 침착 결함, 과무기질화, 골강도 저하,
   조골세포 자멸, 골세포 네트워크 붕괴)
4. 골재형성 커플링: WNT/스클레로스틴 · RANKL-OPG (골세포 SOST 분비,
   LRP5/6 억제, RANKL/OPG 비율 상승, 파골세포 활성화)
5. TGF-베타 신호 과활성 (잠재형 저장고 유리, SMAD2/3, 조골세포 성숙 저해,
   파골세포형성 상승, 근육약화)
6. 골격 표현형 및 Sillence 분류 (Type I-VIII, 반복골절, 진행성 만곡,
   척추압박골절/척추측만, 저신장, 두개저함몰)
7. 골외 증상 (치아형성부전, 청회색 공막, 청력저하, 심장판막/대동맥근 확장,
   관절과이완, 쉬운 멍, 제한성 폐질환)
8. 임상 엔드포인트 (연간 골절율, DXA BMD Z-score, 골대사표지자, 성장 Z-score,
   기능/이동성, 만성통증/QoL, 척추압박골절 수)
9. 비스포스포네이트 PK/PD (파미드로네이트/졸레드론산/경구제, 골결합,
   FPPS 억제, 파골세포 자멸, 급성기반응, 저칼슘혈증, ONJ/비전형골절)
10. 데노수맙 PK/PD (항-RANKL 항체, 표적매개약물동태, 중단 시 반동성 골흡수)
11. 테리파라타이드/아나볼릭 (PTH1R 간헐자극, 골형성 증가, 성장판 개방 금기)
12. 항-스클레로스틴 항체 (세트루수맙, ASTEROID/ORBIT 임상시험)
13. 수술/재활 중재 (골수강내 로드삽입, 척추유합, 물리치료, 비타민D/칼슘)
14. 약물-질환 통합 피드백 (골재형성단위 활성빈도, 커플링 균형, 골무기질
    함량 궤적)
15. (5번 클러스터 내) 항-TGF-베타 항체 프레솔리무맙 투자적 중재

## 3. mrgsolve 모델 (22 ODE 구획)

* **약물 PK (5종, 11개 구획)** — Pamidronate(중심/골결합 2구획),
  Zoledronic acid(중심/골결합 2구획), Denosumab(SC 데포/중심 2구획,
  비선형 Michaelis-Menten 표적매개소실 포함), Teriparatide(SC 데포/중심
  2구획), Setrusumab(중심/말초 2구획 선형 mAb), Fresolimumab(중심 1구획
  선형 mAb).
* **질환/PD (11개 구획)** — 스클레로스틴(SOST), RANKL, OPG, 활성
  TGF-베타(TGFB), 조골세포활성지수(OB), 파골세포활성지수(OC), 골무기질
  함량지수(BMC), 골형성표지자 P1NP, 골흡수표지자 CTX, 누적골절수(FX_CUM),
  신장 Z-score(HEIGHT_Z).
* 세트루수맙은 스클레로스틴에 대해 Emax 경쟁 중화 모델로 유효 SOST 신호를
  낮추고, 프레솔리무맙은 동일한 방식으로 유효 TGF-베타를 낮춘다. 데노수맙은
  유리 RANKL을 Emax 모델로 중화하여 파골세포분화 신호(RANKL/OPG 비율)를
  낮춘다. 비스포스포네이트의 골결합량은 파골세포 자멸 확률(Emax, 약물별
  개별 IC50 — 졸레드론산이 파미드로네이트보다 저용량에서 강력)로 변환된다.
  골절 위험은 기저 대비 BMC 결손에 비례하는 순간 위험도(hazard)로 적분되어
  누적 골절수를 산출한다.

### 7개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 자연경과 (미치료, type III/IV) | Marini 2017 Nat Rev Dis Primers, Rauch 2004 Lancet |
| 2 | IV 파미드로네이트 주기요법 (1mg/kg/day x3일 q3개월) | Glorieux 1998 NEJM |
| 3 | IV 졸레드론산 (0.05mg/kg 6개월마다) | Vuorimies 2017, Barros 2012 |
| 4 | 데노수맙 SC (중증 열성형, off-label) | Hoyer-Kuhn 2014 JBMR |
| 5 | 테리파라타이드 SC (성인 type I 한정, 20ug/day) | Orwoll 2014 JBMR |
| 6 | 세트루수맙 (항-스클레로스틴, 투자적) | Glorieux/ASTEROID·ORBIT 임상시험 |
| 7 | 프레솔리무맙 (항-TGF-베타, 투자적) | Grafe 2014 Nat Med 전임상 근거 |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — Sillence 분류 중증도, 연령, 체중, 성장판 개방 여부.
2. **PK** — 약물별 혈중농도(파미드로네이트/졸레드론산/데노수맙/테리파라
   타이드/세트루수맙/프레솔리무맙).
3. **PD 주요지표** — 유효 스클레로스틴·유리 RANKL, 조골/파골세포 활성지수.
4. **임상 엔드포인트** — BMC 지수(Z-score 변화), 누적 골절수.
5. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
6. **바이오마커** — 골형성표지자(P1NP)·골흡수표지자(CTX) 추이.
7. **성장/골격** — 신장 Z-score 궤적, 테리파라타이드 성장판 금기 경고.
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg oi_qsp_model.dot -o oi_qsp_model.svg
dot -Tpng -Gdpi=150 oi_qsp_model.dot -o oi_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("oi_mrgsolve_model.R")
e_pam <- ev(amt = 30, cmt = "PAM_CENT", ii = 24, addl = 2, rate = 30/2, time = 0) %>%
  ev_repeat(ii = 24*30*3, addl = 12)  # 파미드로네이트 3개월 주기 x4년
out <- mod %>% ev(e_pam) %>% mrgsim(end = 24*365*4, delta = 24)
plot(out, c("BMC", "P1NP", "CTX", "FX_CUM"))

# 3) Shiny 대시보드 실행
shiny::runApp("oi_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| 주기적 IV 파미드로네이트 BMD 상승 | 4년간 요추 BMD Z-score 개선 | Glorieux 1998 NEJM |
| 비스포스포네이트 종합효과 | 골절율/BMD 메타분석 | Dwan 2016 Cochrane (CD005088) |
| 경구 리세드로네이트 12개월 무효과 | 골절 예방효과 미입증 (PLUTO) | Bishop 2013 Lancet |
| 졸레드론산 vs 파미드로네이트 비열등성 | 연 2회 vs 3개월 주기 | Vuorimies 2017 Horm Res Paediatr |
| 데노수맙 중증 열성형 사례군 | 골절 감소, 중단 시 반동성 흡수 | Hoyer-Kuhn 2014 JBMR |
| 테리파라타이드 성인 type I | BMD 증가(성장판 폐쇄 후 한정) | Orwoll 2014 JBMR |
| 항-TGF-베타 중화 전임상 근거 | oim/oim 마우스 골량 정상화 | Grafe 2014 Nat Med |

## 7. 모델 검증 상태

이 컨테이너에는 R/mrgsolve 실행환경이 설치되어 있지 않아(`Rscript` 부재),
mrgsolve 모델은 **문헌 기반 파라미터 설계 및 코드 자체검토(차원/한계값
검사)** 단계까지 완료되었으며 실제 컴파일·적분 실행으로 수치를 검증하지는
못했다. `.dot` 파일은 Graphviz `dot`으로 렌더링해 실제로 SVG/PNG를
생성·확인했다(106 노드, 15 클러스터). mrgsolve/R 환경이 있는 곳에서 위
"실행 방법"대로 실행해 수치 적분 결과를 확인할 것을 권장한다.
