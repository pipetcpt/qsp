# C3 Glomerulopathy (C3 신병증, C3G: DDD & C3GN) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking genetic (CFH/
> CFI/CFB/C3/MCP variants, CFHR5 duplication) and acquired (C3 Nephritic
> Factor [C3NeF], anti-Factor H autoantibody) drivers of uncontrolled
> alternative-pathway (AP) C3-convertase (C3bBb) amplification to glomerular
> C3b deposition (mesangial/subendothelial in C3 glomerulonephritis [C3GN],
> intramembranous ribbon-like dense deposits in Dense Deposit Disease [DDD]),
> mesangial proliferation, podocyte injury, proteinuria, TGF-β-driven
> interstitial fibrosis, nephron loss, eGFR decline, ESKD and kidney
> transplant recurrence — coupled to upstream alternative-pathway inhibitor
> (iptacopan [Factor B], pegcetacoplan [C3/C3b], danicopan [Factor D]) and
> terminal-pathway inhibitor (eculizumab, ravulizumab [C5]) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`c3g_qsp_model.dot`](c3g_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`c3g_qsp_model.svg`](c3g_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`c3g_qsp_model.png`](c3g_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`c3g_mrgsolve_model.R`](c3g_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`c3g_shiny_app.R`](c3g_shiny_app.R) |
| 📚 References             | [`c3g_references.md`](c3g_references.md) |

---

## 1. 질환 개요 (Disease in one paragraph)

C3 신병증(C3 Glomerulopathy, C3G)은 보체 **대안경로(Alternative Pathway,
AP)**의 유전적 또는 후천적 조절 이상으로 인해 C3 전환효소(C3bBb)가
과도하게 안정화되어 사구체에 C3 분해산물이 우세하게 침착되는 희귀 신장
질환군으로, 조직학적으로 막내 리본형 고전자밀도 침착을 특징으로 하는
**치밀침착병(Dense Deposit Disease, DDD)**과 사구체계막/내피하 침착을
특징으로 하는 **C3 사구체신염(C3GN)**의 두 아형으로 나뉜다. 핵심 병인은
**Factor H(CFH)/Factor I(CFI)의 기능소실 변이**로 인한 조절 저하,
**Factor B(CFB)/C3의 기능획득 변이**로 인한 전환효소 안정화, 또는 가장
흔하게는 **C3 신염인자(C3 Nephritic Factor, C3NeF)** — C3bBb 전환효소에
결합해 Factor H 매개 분해를 차단하는 자가항체 — 에 의한 후천적 조절이상이다.
이렇게 안정화된 전환효소는 유체상과 사구체 표면에서 C3를 지속적으로
분해하여 C3b를 침착시키고, 하류로 C5 전환효소 형성 및 막공격복합체
(C5b-9)를 통한 말단경로 활성화로 이어져 메산지움 증식, 족세포 손상,
단백뇨를 유발한다. 병리적으로 축적된 사구체 손상은 TGF-β 매개 세뇨관간질
섬유화와 네프론 소실을 가속화하여 10년 내 약 절반의 환자가 말기신부전
(ESKD)에 도달하며, 신장이식 후에도 최대 50-90%에서 조직학적으로
재발한다. 특징적인 신장외 병변으로 안구의 **드루젠(drusen)** 및 **후천성
부분 지방이상증(acquired partial lipodystrophy)**이 동반될 수 있다. 표준
치료(RAAS 차단제, 면역억제제)는 효과가 제한적이었으나, 최근 대안경로
상류를 직접 차단하는 경구 **Factor B 억제제 iptacopan**과 피하주사
**C3/C3b 억제제 pegcetacoplan**이 3상 임상에서 단백뇨를 유의하게
감소시켰으며, 말단경로만 차단하는 **C5 억제제(eculizumab, ravulizumab)**는
상류 C3 침착이 지속되어 반응이 이질적이라는 점이 기전적으로 중요한
차별점이다.

## 2. 기계론적 지도 클러스터 (16개 클러스터, 120개 노드)

1. 유전적/후천적 소인 (CFH/CFI/CFB/C3/MCP 변이, CFHR5 중복, 단클론감마병증,
   감염 유발 인자, 유전자 패널검사)
2. 자가항체 매개 기전 (C3NeF, 항-Factor H/B/C3b 자가항체, 형질모세포 클론,
   지속성 vs 일시성 역가)
3. 대안경로 증폭 루프 (C3 tick-over → Factor B/D → C3bBb 전환효소 →
   Properdin 안정화 → C3b 침착/C3a)
4. 보체 조절인자 네트워크 (Factor H/I, MCP/CD46, DAF/CD55, CR1, 순 조절
   예비능)
5. 말단경로/막공격복합체 (C5 전환효소 → C5a → C5b-9 MAC → 용해성/비용해성
   신호, sC5b-9)
6. 사구체 침착 조직병리 (메산지움 C3 침착, MPGN 양상, C3GN vs DDD EM 분류,
   반월체, GBM 이중윤곽)
7. 사구체 세포 손상 (메산지움 증식, 내피 활성화, 족세포 C5aR1/C3aR 신호,
   발돌기 소실, 백혈구 침윤)
8. 세뇨관간질 손상 및 섬유화 진행 (RAAS 활성화, TGF-β1/CTGF, 근섬유모세포,
   IF/TA, 사구체경화, 네프론 소실 피드백)
9. 신장외 임상 발현 (드루젠/망막병증, 후천성 부분 지방이상증, 감염 취약성)
10. 신장이식 및 재발 (지속적 전신 AP 조절이상, 동종이식 내 재발, 이식 실패)
11. Iptacopan PK/PD (경구 Factor B 억제제)
12. Pegcetacoplan PK/PD (C3/C3b 억제제)
13. C5 억제제 PK/PD (Eculizumab/Ravulizumab)
14. Danicopan PK/PD (경구 Factor D 억제제, 병용)
15. 바이오마커 (혈청 C3/C4, Factor B/H/I, sC5b-9, C3NeF 역가, AP 기능검사)
16. 임상 엔드포인트 (완전/부분 관해, 단백뇨, eGFR 궤적, ESKD, 이식 생존,
    조직 활성도)

## 3. mrgsolve 모델 (21 ODE 구획)

* **약물 PK (10개 구획)** — Iptacopan 위장관데포/혈장(2), Pegcetacoplan
  피하데포/혈장(2), Eculizumab 중심/말초(2, 2-구획 단클론항체), Ravulizumab
  중심/말초(2, FcRn 재순환), Danicopan 위장관데포/혈장(2).
* **질환/PD (10개 구획)** — 대안경로 전환효소 활성도지수(AP_ACTIVITY),
  혈청 C3(C3_LEVEL), 가용성 말단복합체(SC5B9), 사구체 침착부담
  (GLOM_DEPOSIT), 메산지움 활성도지수(MESANGIAL), 생존 족세포분율
  (PODO_FRAC), 섬유화지표(FIBROSIS), 잔존 네프론분율(NEPHRON_FRAC),
  단백뇨(UPCR), eGFR.
* **시간 추적 (1개 구획)** — FX_YEARS.
* 핵심 설계: **상류(Factor B/D, C3) 억제제**(iptacopan/pegcetacoplan/
  danicopan)는 AP_ACTIVITY 자체를 낮춰 사구체 침착(GLOM_DEPOSIT)과
  말단복합체(SC5B9) 양쪽을 모두 감소시키지만, **말단경로(C5) 억제제**
  (eculizumab/ravulizumab)는 c5_block 항을 통해 SC5B9만 차단하고 상류 C3b
  침착 경로는 그대로 유지되도록 구현했다 — 이는 실제 임상에서 C5 억제제의
  C3G 치료 반응이 이질적이고 제한적인 이유를 기전적으로 반영한다.
  Danicopan은 C5 억제제 배경 위에서 잔여 AP 활성(돌파 활성화)을 부분적으로
  추가 억제하는 병용요법으로 설계했다.

### 10개 시나리오

| # | 시나리오 | 보정 근거 |
|---|---|---|
| 1 | 자연경과 - C3GN 표현형 (SEVERITY=1.0) | Servais 2012 Kidney Int |
| 2 | 자연경과 - DDD 표현형, 고티터 C3NeF (SEVERITY=1.35) | Smith 2019 Nat Rev Nephrol |
| 3 | Iptacopan 200mg BID | Bomback 2025 Lancet (APPEAR-C3G) |
| 4 | Pegcetacoplan 1080mg SC 주2회 | Dixon 2023 Kidney Int Rep |
| 5 | Eculizumab 900mg IV q2w (허가외) | Bomback 2012 CJASN |
| 6 | Ravulizumab 체중기반 IV q8w | Lee 2019 Blood (PK 가교) |
| 7 | Danicopan 150mg TID + Eculizumab 병용 | Risitano 2021 Br J Haematol (PNH EVH 외삽) |
| 8 | CFH 기능소실 변이 + Iptacopan | Gale 2010 Lancet (CFHR5 신병증) |
| 9 | 신장이식 후 재발 + Pegcetacoplan 예방 | Zand 2014 J Am Soc Nephrol |
| 10 | Iptacopan 2년 후 중단 → 재발 | APPEAR-C3G 연장기 임상 설계 외삽 |

## 4. Shiny 대시보드 (8탭)

1. **환자 프로파일** — 표현형(DDD/C3GN 중증도), C3NeF 역가, CFH 변이,
   연령, 시뮬레이션 기간.
2. **PK** — 약물별(Iptacopan/Pegcetacoplan/Eculizumab/Ravulizumab/Danicopan)
   혈장 농도.
3. **PD 주요지표 (보체계)** — AP 활성도지수, 혈청 C3 및 sC5b-9.
4. **임상 엔드포인트** — eGFR 궤적(ESKD 임계선 포함), UPCR, ESKD 도달시간.
5. **사구체 조직학** — 침착부담·메산지움 활성도, 족세포분율·섬유화.
6. **시나리오 비교** — 다중 시나리오 중첩 비교 및 요약표.
7. **바이오마커** — 혈청 C3, sC5b-9, AP 활성도, 침착부담, 메산지움 지표.
8. **참고문헌** — 전체 문헌 목록.

## 5. 실행 방법

```bash
# 1) 기계론적 지도 렌더링
dot -Tsvg c3g_qsp_model.dot -o c3g_qsp_model.svg
dot -Tpng -Gdpi=150 c3g_qsp_model.dot -o c3g_qsp_model.png
```

```r
# 2) R/mrgsolve 시뮬레이션
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("c3g_mrgsolve_model.R") %>% param(SEVERITY = 1.0, C3NEF_TITER = 1.0)
e_ipta <- ev(amt = 200, cmt = "IPTA_GUT", time = 0, ii = 12, addl = 2*365*2-1)
out <- mod %>% ev(e_ipta) %>% mrgsim(end = 24*365*10, delta = 24)  # 10년 추적
plot(out, c("AP_ACTIVITY", "C3_LEVEL", "UPCR", "EGFR"))

# 3) Shiny 대시보드 실행
shiny::runApp("c3g_shiny_app.R")
```

## 6. 주요 임상 보정 근거

| 엔드포인트 | 비교대상 | 근거 |
|---|---|---|
| 자연경과 10년 ESKD 도달률 (~50%) | 등록 코호트 | Servais 2012 Kidney Int |
| Iptacopan 6개월 단백뇨 감소 (35.1% vs 위약) | 무작위 대조 3상 | Bomback 2025 Lancet (APPEAR-C3G) |
| Pegcetacoplan 48주 단백뇨 ≥50% 감소 | 단일군 2상 | Dixon 2023 Kidney Int Rep |
| Eculizumab 반응 이질성 (sC5b-9 반응군에서만 이득) | 개방표지 개념증명 | Bomback 2012 CJASN |
| Ravulizumab C5 포화 및 PK 가교 | PNH 프로그램 | Lee 2019 Blood |
| Danicopan 부가요법 돌파 활성 억제 (PNH EVH 외삽) | 무작위 대조 3상(ALPHA) | Risitano 2021 Br J Haematol |
| 신장이식 후 재발률 (최대 50-90%) | 후향적 코호트 | Zand 2014 J Am Soc Nephrol |
| CFHR5 신병증 상염색체 우성 특이 유전자 중복 | 가족 코호트 | Gale 2010 Lancet |

## 7. 모델 검증 상태

이 컨테이너에는 초기 R/mrgsolve 실행환경이 설치되어 있지 않아(`Rscript`
부재), mrgsolve 모델은 **문헌 기반 파라미터 설계 및 코드 자체검토(21개
컴파트먼트 전부 `$CMT`/`$ODE`/`$INIT`와 1:1 대응 확인, 상류 AP 억제제와
말단(C5) 억제제의 작용점을 의도적으로 분리한 로직 검토)** 단계까지
완료되었으며 실제 컴파일·적분 실행으로 수치를 검증하지는 못했다. `.dot`
파일은 `apt-get install graphviz`로 설치한 Graphviz `dot`으로 실제
렌더링해 SVG/PNG를 생성·확인했다(120 노드, 16 클러스터). 32편의 PubMed
문헌은 WebSearch로 저자/연도/저널 정보를 실제 PubMed·저널 페이지와
대조하여 PMID를 개별 검증했다(학회 초록 등 PMID 미색인 4건은 저널
웹페이지 링크로 대체 표기). mrgsolve/R 환경이 있는 곳에서 위 "실행 방법"
대로 실행해 수치 적분 결과를 확인할 것을 권장한다.
