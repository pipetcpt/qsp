# 부신피질암 (ACC) QSP 모델 — 참고문헌
# Adrenocortical Carcinoma QSP Model — References

모든 PMID 는 작성 시점(2026-08)에 PubMed 에서 직접 조회하여 확인했습니다.
링크는 `https://pubmed.ncbi.nlm.nih.gov/<PMID>/` 형식입니다.

> **모델 파라미터의 근거 표기 규칙.** 아래 각 항목의 *"모델에서의 역할"* 은 그
> 문헌이 `acc_mrgsolve_model.R` 의 어떤 파라미터·구조를 지지하는지를 명시합니다.
> 문헌에서 직접 측정되지 않고 **추정(ASSUMED)** 한 결합은 그 사실을 그대로
> 적었습니다. 특히 **에토포시드 청소율 유도 계수**는 미토테인–에토포시드
> 상호작용 연구에서 측정된 값이 아니라, ① CYP3A4 유도 크기(midazolam probe)와
> ② 에토포시드의 CYP3A4 대사 분율(문헌상 대략 0.3–0.4)에서 유도한 값입니다.

---

## 1. 핵심 임상시험 (Pivotal Clinical Trials)

1. **Fassnacht M, Terzolo M, Allolio B, Baudin E, Haak H, Berruti A, et al.; FIRM-ACT Study Group.**
   Combination chemotherapy in advanced adrenocortical carcinoma.
   *N Engl J Med.* 2012;366(23):2189-97. PMID: **22551107**
   <https://pubmed.ncbi.nlm.nih.gov/22551107/>
   *모델에서의 역할:* EDP-M vs Sz-M 두 시나리오(07, 10)의 보정 기준.
   ORR 23.2% vs 9.2%, PFS 5.0 vs 2.1개월, OS 14.8 vs 12.0개월(유의차 없음).
   모델의 `SLPETO`/`SLPDOX`/`SLPCIS`/`SLPSZ` 는 이 순서를 재현하도록 설정.

2. **Terzolo M, Angeli A, Fassnacht M, Daffara F, Tauchmanova L, Conton PA, et al.**
   Adjuvant mitotane treatment for adrenocortical carcinoma.
   *N Engl J Med.* 2007;356(23):2372-80. PMID: **17554118**
   <https://pubmed.ncbi.nlm.nih.gov/17554118/>
   *모델에서의 역할:* 보조요법 세팅의 근거. RFS 엔드포인트.

3. **Terzolo M, Fassnacht M, Perotti P, Libe R, Kastelan D, Lacroix A, et al.**
   Adjuvant mitotane versus surveillance in low-grade, localised adrenocortical
   carcinoma (ADIUVO): an international, randomised, phase 3 trial.
   *Lancet Diabetes Endocrinol.* 2023;11(10):720-730. PMID: **37619579**
   <https://pubmed.ncbi.nlm.nih.gov/37619579/>
   *모델에서의 역할:* 저위험군에서 보조 미토테인의 이익이 없다는 결과 —
   모델에서 `KI67`/`ENSATst` 에 따라 보조요법 이익이 달라지는 구조의 근거.

4. **Fassnacht M, Berruti A, Baudin E, Demeure MJ, Gilbert J, Haak H, et al.**
   Linsitinib (OSI-906) versus placebo for patients with locally advanced or
   metastatic adrenocortical carcinoma: a double-blind, randomised, phase 3 study.
   *Lancet Oncol.* 2015;16(4):426-35. PMID: **25795408**
   <https://pubmed.ncbi.nlm.nih.gov/25795408/>
   *모델에서의 역할:* 시나리오 17. IGF2 가 ~90%에서 과발현되는데도 IGF1R
   억제가 실패한 사실 — 모델은 이를 ① `APROLIF` 가 작다(ACC 성장이 IGF에
   의해 제한되지 않음)와 ② 인슐린-IR-A 구제 경로의 조합으로 설명.

5. **Jones RL, Kim ES, Nava-Parada P, Alam S, Johnson FM, Stephens AW, et al.**
   Phase I study of intermittent oral dosing of the insulin-like growth factor-1
   and insulin receptors inhibitor OSI-906 in patients with advanced solid tumors.
   *Clin Cancer Res.* 2015;21(4):693-700. PMID: **25208878**
   <https://pubmed.ncbi.nlm.nih.gov/25208878/>
   *모델에서의 역할:* 린시티닙 PK(`CLLIN`, `V1LIN`, `FLIN`)와 용량제한
   고혈당의 근거.

6. **Habra MA, Stephen B, Campbell M, Hess K, Tapia C, Xu M, et al.**
   Phase II clinical trial of pembrolizumab efficacy and safety in advanced
   adrenocortical carcinoma.
   *J Immunother Cancer.* 2019;7(1):253. PMID: **31533818**
   <https://pubmed.ncbi.nlm.nih.gov/31533818/>
   *모델에서의 역할:* 시나리오 14–16. ORR ~14–23% 수준을 재현하도록
   `SLPIMM`/`EMAXPD1` 설정.

7. **Raj N, Zheng Y, Kelly V, Katz SS, Chou J, Do RKG, et al.**
   PD-1 Blockade in Advanced Adrenocortical Carcinoma.
   *J Clin Oncol.* 2020;38(1):71-80. PMID: **31644329**
   <https://pubmed.ncbi.nlm.nih.gov/31644329/>
   *모델에서의 역할:* 위와 동일. 스테로이드 과다 환자에서 반응이 낮다는
   임상 관찰이 결과 F 의 출발점.

8. **Berruti A, Terzolo M, Pia A, Angeli A, Dogliotti L.**
   Mitotane associated with etoposide, doxorubicin, and cisplatin in the
   treatment of advanced adrenocortical carcinoma.
   *Cancer.* 1998;83(10):2194-200. PMID: **9827725**
   <https://pubmed.ncbi.nlm.nih.gov/9827725/>
   *모델에서의 역할:* EDP-M 용법(에토포시드 100 mg/m² d2-4, 독소루비신
   40 mg/m² d1, 시스플라틴 40 mg/m² d3-4, q28d)의 근거.

9. **Khan TS, Imam H, Juhlin C, Skogseid B, Grondal S, Tibblin S, et al.**
   Streptozocin and o,p'DDD in the treatment of adrenocortical cancer patients:
   long-term survival in its adjuvant use.
   *Ann Oncol.* 2000;11(10):1281-7. PMID: **11106117**
   <https://pubmed.ncbi.nlm.nih.gov/11106117/>
   *모델에서의 역할:* Sz-M 용법(1 g d1-5 유도 후 2 g q3wk)의 근거.

10. **Megerle F, Kroiss M, Hahner S, Fassnacht M, et al.**
    Mitotane Monotherapy in Patients With Advanced Adrenocortical Carcinoma.
    *J Clin Endocrinol Metab.* 2018;103(4):1686-1695. PMID: **29452402**
    <https://pubmed.ncbi.nlm.nih.gov/29452402/>
    *모델에서의 역할:* 시나리오 05. 저종양부하에서 단독요법의 이익 —
    모델의 `SLPMIT`(미토테인 자체 세포독성이 약함)의 근거.

---

## 2. 미토테인 약동학 · 치료약물모니터링 (Mitotane PK & TDM) — 모델의 심장부

11. **Arshad U, Taubert M, Kurlbaum M, Frechen S, Herterich S, Megerle F, et al.**
    Enzyme autoinduction by mitotane supported by population pharmacokinetic
    modelling in a large cohort of adrenocortical carcinoma patients.
    *Eur J Endocrinol.* 2018;179(5):287-297. PMID: **30087117**
    <https://pubmed.ncbi.nlm.nih.gov/30087117/>
    *모델에서의 역할:* **이 모델의 중심 근거.** ① 겉보기 분포용적 ~6,086 L
    (개체간 변동 **CV 81.5%**) → `V1MIT + VLEANTIS + KPFAT*FATKG` = 5,966 L
    로 재현. ② **BMI 가 미토테인 처분에 영향**을 준다 → 지방량을 저장고
    용적의 결정 변수로 둔 근거. ③ 청소율의 **선형 효소 자기유도** →
    `FMIND`. ④ 고용량 시작 요법과 **16일째 첫 TDM** 권고 → 시나리오 01–04.

12. **Terzolo M, Baudin AE, Ardito A, Kroiss M, Leboulleux S, Daffara F, et al.**
    Mitotane levels predict the outcome of patients with adrenocortical
    carcinoma treated adjuvantly following radical resection.
    *Eur J Endocrinol.* 2013;169(3):263-70. PMID: **23704714**
    <https://pubmed.ncbi.nlm.nih.gov/23704714/>
    *모델에서의 역할:* 치료 목표 농도 **≥14 mg/L** 와 노출–반응 관계.
    모델의 `INWIN`/`TIW` 출력과 결과 A·B 의 기준선.

13. **Kerkhofs TM, Baudin E, Terzolo M, Allolio B, Chadarevian R, Mueller HH, et al.**
    Comparison of two mitotane starting dose regimens in patients with advanced
    adrenocortical carcinoma.
    *J Clin Endocrinol Metab.* 2013;98(12):4759-67. PMID: **24057287**
    <https://pubmed.ncbi.nlm.nih.gov/24057287/>
    *모델에서의 역할:* **결과 B 의 핵심.** 고용량 시작이 미토테인 농도나
    이상반응률에서 유의한 차이를 보이지 않았다는 RCT 결과. 모델은 이것을
    "용량은 랜덤화했지만 체성분은 층화하지 않았다 + 보고된 엔드포인트(목표
    도달 시간)가 공변량이 가장 크게 해치는 환자를 검열해 버린다"로 설명하며,
    popPK(11번)의 고용량 권고와 모순되지 않음을 보입니다.

14. **Kerkhofs TM, Derijks LJ, Ettaieb MH, Eekhoff EM, Neef C, Gelderblom H, et al.**
    Short-term variation in plasma mitotane levels confirms the importance of
    trough level monitoring.
    *Eur J Endocrinol.* 2014;171(6):677-83. PMID: **25201518**
    <https://pubmed.ncbi.nlm.nih.gov/25201518/>
    *모델에서의 역할:* 일중 변동이 제한적이라는 관찰 → 중심구획을 혈액 +
    급속평형 제지방 조직으로 크게(`V1MIT` 400 L) 잡아 일중 진폭을 억제한
    구조적 결정의 근거.

15. **Puglisi S, Calabrese A, Basile V, Ceccato F, Scaroni C, Altieri B, et al.**
    Mitotane Concentrations Influence Outcome in Patients with Advanced
    Adrenocortical Carcinoma.
    *Cancers (Basel).* 2020;12(3):740. PMID: **32245135**
    <https://pubmed.ncbi.nlm.nih.gov/32245135/>
    *모델에서의 역할:* 진행성 질환에서도 농도–결과 관계가 성립.

16. **Puglisi S, Perotti P, Cosentini D, Fiorentini C, Basile V, et al.**
    Decreased Plasma Mitotane Levels Are Associated with Improved
    Recurrence-Free Survival in Adrenocortical Carcinoma Patients Treated
    Adjuvantly — target attainment analysis.
    *J Clin Med.* 2019;8(11):1850. PMID: **31684071**
    <https://pubmed.ncbi.nlm.nih.gov/31684071/>
    *모델에서의 역할:* 보조요법에서 목표 농도 도달 여부와 재발 위험의 관계.

17. **Corso CR, Acco A, Bach C, Bonatto SJR, de Figueiredo BC, de Souza LM.**
    Pharmacological profile and effects of mitotane in adrenocortical carcinoma.
    *Br J Clin Pharmacol.* 2021;87(7):2698-2710. PMID: **33382119**
    <https://pubmed.ncbi.nlm.nih.gov/33382119/>
    *모델에서의 역할:* 종말 반감기 범위(대략 18–160일), 극단적 친유성,
    o,p'-DDA / o,p'-DDE 대사체. 모델의 유도 전(86일) / 완전 유도 후(42일)
    반감기가 이 범위 안에 들어오는지 확인하는 기준.

18. **Hescot S, Paci A, Seck A, Slama A, Viengchareun S, Trabado S, et al.**
    Lipoprotein-Free Mitotane Exerts High Cytotoxic Activity in Adrenocortical
    Carcinoma.
    *J Clin Endocrinol Metab.* 2015;100(8):2890-2898. PMID: **26120791**
    <https://pubmed.ncbi.nlm.nih.gov/26120791/>
    *모델에서의 역할:* 지단백에 결합하지 않은 분획이 세포독성을 낸다는 점 —
    "혈장 총 농도가 아니라 지질상과의 분배가 작용을 결정한다"는 저장고
    관점을 세포 수준에서 지지.

---

## 3. 미토테인 작용기전 — 부신파괴 (Mechanism of Adrenolysis)

19. **Sbiera S, Leich E, Liebisch G, Sbiera I, Schirbel A, Wiemer L, et al.**
    Mitotane Inhibits Sterol-O-Acyl Transferase 1 Triggering Lipid-Mediated
    Endoplasmic Reticulum Stress and Apoptosis in Adrenocortical Carcinoma Cells.
    *Endocrinology.* 2015;156(11):3895-908. PMID: **26305886**
    <https://pubmed.ncbi.nlm.nih.gov/26305886/>
    *모델에서의 역할:* 부신파괴 경로의 분자 기전. 모델 cluster 6 의
    `SOAT1INH → OXYSTER → ERSTRESS → CASP → ADRLYSIS` 사슬과
    `KLYS`/`EC50LYS`/`HLYS` 의 근거.

20. **Weigand I, Altieri B, Lacombe AMF, Basile V, Kircher S, Landwehr LS, et al.**
    Expression of SOAT1 in Adrenocortical Carcinoma and Response to Mitotane
    Monotherapy.
    *J Clin Endocrinol Metab.* 2020;105(8):dgaa293. PMID: **32449514**
    <https://pubmed.ncbi.nlm.nih.gov/32449514/>
    *모델에서의 역할:* 표적(SOAT1) 발현량이 반응과 관련 → 모델의
    `SOAT1exp` 노드와 반응 예측 변동성의 근거.

21. **Smith DC, Kroiss M, Kebebew E, Habra MA, Chugh R, Schneider BJ, et al.**
    A phase 1 study of nevanimibe HCl, a novel SOAT1 inhibitor, in
    adrenocortical carcinoma.
    *Invest New Drugs.* 2020;38(5):1421-1429. PMID: **31984451**
    <https://pubmed.ncbi.nlm.nih.gov/31984451/>
    *모델에서의 역할:* 선택적 SOAT1 억제제가 부신파괴 축만 자극하는 대조
    분기(지도의 `NEVANIM` 노드) — "미토테인의 이익 축을 유도 없이 얻으면
    어떻게 되는가"라는 반사실 실험의 임상적 근거.

22. **Weigand I, Ronchi CL, Rizk-Rabin M, Dalmazi GD, Wild V, Bathon K, et al.**
    Active steroid hormone synthesis renders adrenocortical cells highly
    susceptible to type II ferroptosis induction.
    *Cell Death Dis.* 2020;11(3):192. PMID: **32184394**
    <https://pubmed.ncbi.nlm.nih.gov/32184394/>
    *모델에서의 역할:* 활발히 스테로이드를 합성하는 세포가 더 취약하다는 점 —
    모델에서 `CORT_SYN → FERROPT` 로 연결(분비형 종양이 미토테인에 더
    민감할 수 있는 기전).

23. **Ronchi CL, Sbiera S, Volante M, Steinhauer S, Scott-Wild V, Altieri B, et al.**
    CYP2W1 is highly expressed in adrenal glands and is positively associated
    with the response to mitotane in adrenocortical carcinoma.
    *PLoS One.* 2014;9(8):e105855. PMID: **25144458**
    <https://pubmed.ncbi.nlm.nih.gov/25144458/>
    *모델에서의 역할:* 반응 예측 약물유전체 표지자(지도의 `CYP2W1` 노드).

24. **van Koetsveld PM, Creemers SG, Dogan F, Franssen GJH, de Herder WW, Hofland LJ, et al.**
    The Efficacy of Mitotane in Human Primary Adrenocortical Carcinoma Cultures.
    *J Clin Endocrinol Metab.* 2020;105(2):407-417. PMID: **31586196**
    <https://pubmed.ncbi.nlm.nih.gov/31586196/>
    *모델에서의 역할:* in vitro 감수성의 개체간 차이 — 가상 집단 레이어에서
    반응 변동을 두는 근거.

25. **Krüger AF, et al.**
    Adrenocortical Mitochondria-Associated Membranes and the Lipidoproteomic
    Response to Mitotane.
    *J Endocr Soc.* 2025;10(1):bvaf198. PMID: **41472891**
    <https://pubmed.ncbi.nlm.nih.gov/41472891/>
    *모델에서의 역할:* 미토콘드리아-ER 접촉면과 지질 재배치 — 지도의
    `MITOCHOND`/`ERSTRESS` 결합의 최신 근거.

---

## 4. CYP3A4 유도와 약물상호작용 (Induction & DDI) — 해악 축 #1

26. **Kroiss M, Quinkler M, Lutz WK, Allolio B, Fassnacht M.**
    Drug interactions with mitotane by induction of CYP3A4 metabolism in the
    clinical management of adrenocortical carcinoma.
    *Clin Endocrinol (Oxf).* 2011;75(5):585-91. PMID: **21883349**
    <https://pubmed.ncbi.nlm.nih.gov/21883349/>
    *모델에서의 역할:* `EMAXIND`/`EC50IND` 및 유도의 임상적 광범위성
    (지도 cluster 8 의 피해약물 목록).

27. **van Erp NP, Guchelaar HJ, Ploeger BA, Romijn JA, den Hartigh J, Gelderblom H.**
    Mitotane has a strong and a durable inducing effect on CYP3A4 activity.
    *Eur J Endocrinol.* 2011;164(4):621-6. PMID: **21220434**
    <https://pubmed.ncbi.nlm.nih.gov/21220434/>
    *모델에서의 역할:* midazolam probe 로 측정한 **유도의 크기와 지속성**.
    모델에서 ENZ 가 ~4배까지 오르고 투약 중단 후에도 수개월 유지되는
    (결과 H) 거동의 직접 근거.

28. **Theile D, Haefeli WE, Weiss J.**
    Effects of adrenolytic mitotane on drug elimination pathways assessed
    in vitro.
    *Endocrine.* 2015;49(3):842-53. PMID: **25542188**
    <https://pubmed.ncbi.nlm.nih.gov/25542188/>
    *모델에서의 역할:* CYP3A4 외 CYP2B6·UGT·SULT·P-gp 동반 유도 →
    `ENZ2B6`, `PGP`/`FPGPIND`, `DDI_LEVO` 노드.

29. **Kroiss M, Quinkler M, Johanssen S, van Erp NP, Lankheet N, Pöllinger A, et al.**
    Sunitinib in refractory adrenocortical carcinoma: a phase II, single-arm,
    open-label trial (SIRAC).
    *J Clin Endocrinol Metab.* 2012;97(10):3495-503. PMID: **22851488**
    <https://pubmed.ncbi.nlm.nih.gov/22851488/>
    *모델에서의 역할:* 미토테인 유도가 병용 TKI 노출을 크게 낮춘다는
    임상 증거(지도의 `DDI_SUNI`). 결과 D 의 일반화 사례.

---

## 5. 스테로이드 보충 · 결합단백 · 검사 함정 (Replacement, Binding Proteins, Assay)

30. **Chortis V, Taylor AE, Schneider P, Tomlinson JW, Hughes BA, O'Neil DM, et al.**
    Mitotane therapy in adrenocortical cancer induces CYP3A4 and inhibits
    5alpha-reductase, explaining the need for personalized glucocorticoid and
    androgen replacement.
    *J Clin Endocrinol Metab.* 2013;98(1):161-71. PMID: **23162091**
    <https://pubmed.ncbi.nlm.nih.gov/23162091/>
    *모델에서의 역할:* **결과 E 의 직접 근거.** 코르티솔 청소율 유도로
    하이드로코르티손 필요량이 대략 두 배가 된다(`FMHC`), 그리고
    5α-환원효소 억제로 남성 안드로겐 보충이 필요해진다(`INH5AR`).

31. **Nader N, Raverot G, Emptoz-Bonneton A, Déchaud H, Bonnay M, Baudin E, Pugeat M.**
    Mitotane has an estrogenic effect on sex hormone-binding globulin and
    corticosteroid-binding globulin in humans.
    *J Clin Endocrinol Metab.* 2006;91(6):2165-70. PMID: **16551731**
    <https://pubmed.ncbi.nlm.nih.gov/16551731/>
    *모델에서의 역할:* **결과 E 의 나머지 절반.** CBG·SHBG 가 대략 두 배로
    상승(`EMAXCBG`, `EMAXSHBG`). 모델은 이것이 유리 코르티솔을 크게 바꾸지는
    않지만 **총/유리 비를 16.9 → 37.2 로 올려** 총 코르티솔 판독을 오도한다는
    것을 정량적으로 보입니다 — 즉 "환자의 오차"(유도)와 "검사의 오차"(CBG)는
    다른 종류이며 같은 방향을 가리킵니다.

32. **Daffara F, De Francia S, Reimondo G, Zaggia B, Aroasio E, Porpiglia F, et al.**
    Prospective evaluation of mitotane toxicity in adrenocortical cancer
    patients treated adjuvantly.
    *Endocr Relat Cancer.* 2008;15(4):1043-53. PMID: **18824557**
    <https://pubmed.ncbi.nlm.nih.gov/18824557/>
    *모델에서의 역할:* 독성 프로파일 — GGT/ALT 상승(`AALT`), 유리 T4 저하
    (`AFT4`), 소화기 불내성, 부신 기능저하의 사실상 보편성(`KLYS` 결과로
    부신 피질이 거의 완전히 소실되는 것과 대조).

---

## 6. 코르티솔과 면역회피 (Glucocorticoid-Mediated Immune Escape) — 결과 F

33. **Landwehr LS, Altieri B, Schreiner J, Sbiera I, Weigand I, Kroiss M, et al.**
    Interplay between glucocorticoids and tumor-infiltrating lymphocytes on the
    prognosis of adrenocortical carcinoma.
    *J Immunother Cancer.* 2020;8(1):e000469. PMID: **32474412**
    <https://pubmed.ncbi.nlm.nih.gov/32474412/>
    *모델에서의 역할:* **결과 F 의 직접 근거.** 글루코코르티코이드 과다가
    종양침윤림프구를 고갈시키고 예후를 악화. 모델의
    `GR → GCIMSUP → TEFF(억제)` 결합과 `GCIC50`/`HGC`.

34. **Baechle JJ, Hanna DN, Sekhar KR, Rathmell JC, Rathmell WK, Baregamian N.**
    Integrative computational immunogenomic profiling of cortisol-secreting
    adrenocortical carcinoma.
    *J Cell Mol Med.* 2021;25(21):10061-10072. PMID: **34664400**
    <https://pubmed.ncbi.nlm.nih.gov/34664400/>
    *모델에서의 역할:* 분비형 ACC 의 면역 배제 표현형.

35. **Wu K, Liu X, Wang Y, Liu Y, Zhang Y, Wang Z, et al.**
    Discovery of a glucocorticoid receptor (GR) activity signature correlates
    with immune cell infiltration in adrenocortical carcinoma.
    *J Immunother Cancer.* 2023;11(10):e007528. PMID: **37793855**
    <https://pubmed.ncbi.nlm.nih.gov/37793855/>
    *모델에서의 역할:* GR 활성 자체가 침윤과 상관 → 모델이 총 코르티솔이
    아니라 **GR 점유율(`GRO`)** 을 면역 억제의 구동 변수로 쓴 근거.

36. **Landwehr LS, Schreiner J, Appenzeller S, Kircher S, Herterich S, Sbiera S, et al.**
    Expression and Prognostic Relevance of PD-1, PD-L1, and CTLA-4 Immune
    Checkpoints in Adrenocortical Carcinoma.
    *J Clin Endocrinol Metab.* 2024;109(9):2325-2334. PMID: **38415841**
    <https://pubmed.ncbi.nlm.nih.gov/38415841/>
    *모델에서의 역할:* 관문 분자 발현과 예후 — `PD1`/`CTLA4` 노드.

37. **Maier T, et al.**
    Wnt/beta-catenin pathway activation is associated with glucocorticoid
    secretion in adrenocortical carcinoma, but not directly with immune cell
    infiltration.
    *Front Endocrinol (Lausanne).* 2025;16:1502117. PMID: **40130164**
    <https://pubmed.ncbi.nlm.nih.gov/40130164/>
    *모델에서의 역할:* Wnt 활성과 분비 표현형의 연관(지도의
    `WNT → STEROIDOG_PHEN`), 그리고 면역 침윤과는 **직접** 연결되지 않는다는
    점 — 모델이 Wnt→면역을 직접 연결하지 않고 코르티솔을 매개로만 연결한
    구조적 선택의 근거.

---

## 7. 유전체 · 예후 (Genomics & Prognosis)

38. **Assié G, Letouzé E, Fassnacht M, Jouinot A, Luscap W, Barreau O, et al.**
    Integrated genomic characterization of adrenocortical carcinoma.
    *Nat Genet.* 2014;46(6):607-12. PMID: **24747642**
    <https://pubmed.ncbi.nlm.nih.gov/24747642/>
    *모델에서의 역할:* 드라이버 빈도(지도 cluster 1)와 분자 아형.

39. **Zheng S, Cherniack AD, Dewal N, Moffitt RA, Danilova L, Murray BA, et al.;
    Cancer Genome Atlas Research Network.**
    Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma.
    *Cancer Cell.* 2016;29(5):723-736. PMID: **27165744**
    <https://pubmed.ncbi.nlm.nih.gov/27165744/>
    *모델에서의 역할:* TP53 ~21%, ZNRF3 ~19%, CTNNB1 ~16%, IGF2 과발현
    ~90% 등 지도의 유전체 노드 수치.

40. **Beuschlein F, Weigel J, Saeger W, Kroiss M, Wild V, Daffara F, et al.**
    Major prognostic role of Ki67 in localized adrenocortical carcinoma after
    complete resection.
    *J Clin Endocrinol Metab.* 2015;100(3):841-9. PMID: **25559399**
    <https://pubmed.ncbi.nlm.nih.gov/25559399/>
    *모델에서의 역할:* `KI67 → KGROW` 결합과 보조요법 대상 선정.

41. **Arlt W, Biehl M, Taylor AE, Hahner S, Libé R, Hughes BA, et al.**
    Urine steroid metabolomics as a biomarker tool for detecting malignancy in
    adrenal tumors.
    *J Clin Endocrinol Metab.* 2011;96(12):3775-3784. PMID: **21917861**
    <https://pubmed.ncbi.nlm.nih.gov/21917861/>
    *모델에서의 역할:* 지도의 `URSTER` 노드. 모델이 DHEAS·17-OHP·
    11-deoxycortisol 을 별도 출력으로 둔 이유(효소 차단 시 기질이 쌓이는
    패턴 자체가 판독 정보).

---

## 8. 종설 · 진료지침 (Reviews & Guidelines)

42. **Fassnacht M, Dekkers OM, Else T, Baudin E, Berruti A, de Krijger RR, et al.**
    European Society of Endocrinology Clinical Practice Guidelines on the
    management of adrenocortical carcinoma in adults, in collaboration with the
    European Network for the Study of Adrenal Tumors.
    *Eur J Endocrinol.* 2018;179(4):G1-G46. PMID: **30299884**
    <https://pubmed.ncbi.nlm.nih.gov/30299884/>
    *모델에서의 역할:* 치료 알고리즘, TDM 목표, 스테로이드 보충 권고 —
    시나리오 설계 전반의 임상 골격.

43. **Fassnacht M, Assié G, Baudin E, Eisenhofer G, de la Fouchardiere C, Haak HR, et al.**
    Adrenocortical carcinomas and malignant phaeochromocytomas: ESMO-EURACAN
    Clinical Practice Guidelines for diagnosis, treatment and follow-up.
    *Ann Oncol.* 2020;31(11):1476-1490. PMID: **32861807**
    <https://pubmed.ncbi.nlm.nih.gov/32861807/>
    *모델에서의 역할:* 진행성 질환 치료 순서(EDP-M 1차)의 근거.

44. **Else T, Kim AC, Sabolch A, Raymond VM, Kandathil A, Caoili EM, et al.**
    Adrenocortical carcinoma.
    *Endocr Rev.* 2014;35(2):282-326. PMID: **24423978**
    <https://pubmed.ncbi.nlm.nih.gov/24423978/>
    *모델에서의 역할:* 병태생리 종설 — 지도 전반의 배경.

45. **Fassnacht M, Kroiss M, Allolio B.**
    Update in adrenocortical carcinoma.
    *J Clin Endocrinol Metab.* 2013;98(12):4551-64. PMID: **24081734**
    <https://pubmed.ncbi.nlm.nih.gov/24081734/>
    *모델에서의 역할:* 임상 관리 종설.

46. **Allolio B, Fassnacht M.**
    Clinical review: Adrenocortical carcinoma: clinical update.
    *J Clin Endocrinol Metab.* 2006;91(6):2027-37. PMID: **16551738**
    <https://pubmed.ncbi.nlm.nih.gov/16551738/>
    *모델에서의 역할:* 호르몬 과다 분비 표현형의 빈도(약 60%가 기능성) —
    `SECRETOR` 스위치의 근거.

---

## 9. 모델에서 명시적으로 **추정(ASSUMED)** 한 결합

투명성을 위해, 문헌에서 직접 측정되지 않은 결합을 모두 여기에 모았습니다.

| 모델 요소 | 값 | 근거의 성격 |
|---|---|---|
| `FM3A4ETO` (에토포시드의 CYP3A4 대사 분율) | 0.35 | **ASSUMED.** 미토테인–에토포시드 상호작용 연구는 확인하지 못했습니다. 26·27번의 CYP3A4 유도 크기와 에토포시드가 부분적으로 CYP3A4 대사·부분적으로 신배설된다는 일반 약리 지식에서 유도. 결과: ENZ≈4 에서 청소율 약 2배 → 노출 약 1/2. 이 계수는 결과 D 의 크기를 직접 정하므로 **가장 먼저 실측으로 대체되어야 하는 파라미터**입니다. |
| `KPFAT` (지방 1 kg 당 분배용적 203 L/kg) | 203 | **부분 추정.** 11번의 겉보기 Vss ~6,086 L 와 BMI 공변량 효과를 재현하도록 역산. 조직 분배계수 실측치가 아닙니다. |
| `AFAT` (코르티솔 과다 → 지방 증가) | 0.35 | **ASSUMED.** 쿠싱증후군의 체성분 변화 크기에서 대략 설정. 결과 C 의 크기를 정합니다(모델은 이 되먹임이 치료 중에는 거의 작동하지 않고 진단 전 기간에 집중된다고 예측). |
| `IC50RB` < `IC50RA` (대사 IR-B 가 종양 IR-A 보다 린시티닙에 더 민감) | 0.35 vs 1.50 | **ASSUMED, 기전적으로 동기부여됨.** 5번의 용량제한 고혈당이 항종양 효과보다 먼저 나타난다는 관찰을 "대사 수용체가 먼저 차단된다"로 해석. 결과 G 의 인슐린 구제 크기를 정합니다. |
| `SLPIMM`, `EMAXPD1` | 0.010, 2.2 | **보정값.** 6·7번의 ORR ~14–23% 를 재현하도록 설정. |
| `CNSTHR` = 24 mg/L (Hill EC50) | 24 | **구조적 선택.** 임상 경보 수준 20 mg/L 를 Hill 곡선의 EC50 이 아니라 **상승 구간**에 놓기 위해 EC50 을 위로 옮겼습니다. 20 mg/L 에서 최대 손상의 약 17%. |
| 종양 성장 `KGROW`, `TVMAX` | 0.0047, 20,000 mL | **보정값.** 미치료 시 부피 배가시간 약 50–60일이 되도록 설정. ACC 의 실제 성장률은 매우 이질적입니다. |

---

## 10. 모델 구조와 문헌의 대응 요약

| 모델 결과 | 문헌 근거 번호 | 성격 |
|---|---|---|
| Vss ~6,000 L, 대부분 지방 | 11, 17, 18 | 직접 지지 |
| 치료창 14–20 mg/L | 12, 15, 16 | 직접 지지 |
| A. 지방이 τ 만 바꾸고 Css 는 안 바꾼다 | 11 (BMI 공변량) | 문헌 + 모델 연역 |
| B. 용량 RCT vs 층화되지 않은 체성분 | 11 (CV 81.5%), 13 (RCT 무의차) | **모델의 새로운 설명** |
| C. 진단 전 쿠싱이 저장고를 미리 키운다 | 46 (기능성 60%), 32 | **모델의 새로운 예측** |
| D. EDP-M 은 자기 파트너를 대사시킨다 | 26, 27, 28, 29 + ASSUMED 계수 | 기전은 확립, 크기는 추정 |
| E. 총 코르티솔이 거짓말한다 | 30 (유도), 31 (CBG) | 직접 지지, 정량화는 모델 |
| F. 종양이 스스로 스테로이드 전처치를 한다 | 33, 34, 35, 37 | 직접 지지 |
| G. 표적은 있는데 약은 실패했다 | 4, 5 | 문헌 + 모델의 기전적 설명 |
| H. washout 이 없다 | 27, 17 | 직접 지지 |

---

## 11. 도구 (Tools)

- **mrgsolve** — Baron KT. *mrgsolve: Simulate from ODE-Based Models.*
  <https://mrgsolve.org/>
- **Graphviz** — <https://graphviz.org/>
- **Shiny** — <https://shiny.posit.co/>
- Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO.
  Model of chemotherapy-induced myelosuppression with parameter consistency
  across drugs. *J Clin Oncol.* 2002;20(24):4713-21. PMID: **12488418**
  <https://pubmed.ncbi.nlm.nih.gov/12488418/>
  *모델에서의 역할:* 호중구 4-transit 구획 구조(`PROLN`–`TR3N`–`CIRCN`),
  `MTTN`/`GAMN`.

---

> **면책.** 본 모델과 참고문헌 정리는 교육·연구 목적입니다. 임상 의사결정,
> 처방, 규제 제출에 사용하지 마십시오. 파라미터는 설명을 위한 근사치이며
> 실제 환자 데이터에 대한 적합·검증이 별도로 필요합니다.
