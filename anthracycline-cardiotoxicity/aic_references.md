# 안트라사이클린 심장독성 (Anthracycline-Induced Cardiotoxicity) QSP 모델 — 참고문헌 (References)

안트라사이클린 심장독성의 기계론적 지도(`aic_qsp_model.dot`), mrgsolve ODE 모델
(`aic_mrgsolve_model.R`), Shiny 대시보드(`aic_shiny_app.R`), 독립 수치검증
스크립트(`aic_reference_check.py`) 구축에 사용한 문헌 목록입니다. 총 62편.

> **PMID 검증에 관한 註**: 아래 62편 중 **61편은 PubMed E-utilities
> (`esearch`/`esummary`) 로 조회하여 제목·저자·저널·연도를 실제 레코드와 대조
> 확인한 PMID 직접 링크**입니다. 추측으로 기입한 식별자는 없습니다. 나머지
> 1편(51번 Slamon 2001)은 자동 조회에서 단일 레코드로 확정되지 않아 PMID를
> 추측하는 대신 **검색 질의 URL**로 두었습니다 — 잘못된 PMID는 전혀 다른 논문을
> 인용하게 되기 때문입니다.
>
> **중요**: `README.md` 와 모델 파일에 "모델:" 로 표시된 수치는 모두 이 모델의
> **계산 결과**이며 문헌 값이 아닙니다. 문헌은 파라미터 보정(calibration)의
> 근거 및 목표값(anchor)으로만 사용되었고, 그 대응 관계는
> `aic_mrgsolve_model.R` 상단의 CALIBRATION ANCHORS 절에 정리되어 있습니다.

---

## 1. 역학 · 누적용량-반응 · 자연경과 (Epidemiology, dose-response, natural history)

모델의 A1 분석(누적용량-발생률 곡선)과 후기 진행 성분(KDNH, FIB)의 보정 근거.

1. Von Hoff DD, et al. **Risk factors for doxorubicin-induced congestive heart
   failure.** *Ann Intern Med.* 1979. — 누적용량-심부전 관계의 원전; 모델
   A1의 1차 anchor. <https://pubmed.ncbi.nlm.nih.gov/496103/>
2. Swain SM, Whaley FS, Ewer MS. **Congestive heart failure in patients treated
   with doxorubicin: a retrospective analysis of three trials.** *Cancer.* 2003.
   — 630 mg/m²에서 심부전 약 5%, 용량-반응의 현대적 재분석.
   <https://pubmed.ncbi.nlm.nih.gov/12767102/>
3. Cardinale D, et al. **Early detection of anthracycline cardiotoxicity and
   improvement with heart failure therapy.** *Circulation.* 2015. — n=2625,
   CTRCD 9%, 98%가 1년 내 발생, 중앙 3.5개월. 모델의 시간척도 anchor.
   <https://pubmed.ncbi.nlm.nih.gov/25948538/>
4. Cardinale D, et al. **Anthracycline-induced cardiomyopathy: clinical
   relevance and response to pharmacologic therapy.** *J Am Coll Cardiol.*
   2010. — 치료 개시 시점별 LVEF 회복률(2개월 내 64% → 6개월 이후 ~0%);
   모델 A5(가역성 창)의 anchor. <https://pubmed.ncbi.nlm.nih.gov/20117401/>
5. Lipshultz SE, et al. **Chronic progressive cardiac dysfunction years after
   doxorubicin therapy for childhood acute lymphoblastic leukemia.**
   *J Clin Oncol.* 2005. — 소아 생존자의 지속적 진행; 비가역 풀 + 성장
   불균형 개념의 근거. <https://pubmed.ncbi.nlm.nih.gov/15837978/>
6. Mulrooney DA, et al. **Major cardiac events for adult survivors of childhood
   cancer diagnosed between 1970 and 1999.** *BMJ.* 2020. — 장기 코호트의
   후기 심장사건. <https://pubmed.ncbi.nlm.nih.gov/31941657/>
7. Cardinale D, Iacopo F, Cipolla CM. **Cardiotoxicity of anthracyclines.**
   *Front Cardiovasc Med.* 2020. — 임상 종합 리뷰.
   <https://pubmed.ncbi.nlm.nih.gov/32258060/>
8. Camilli M, et al. **Anthracyclines, diastolic dysfunction and the road to
   heart failure in cancer survivors.** *Prog Cardiovasc Dis.* 2024. — 이완기
   기능장애·섬유화 경로. <https://pubmed.ncbi.nlm.nih.gov/39025347/>

## 2. Top2b 의존 DNA 손상 축 (Topoisomerase IIβ arm) — 첨두 구동 축의 근거

모델의 `TOXN → DSB → P53` 축과 덱스라족산 Top2b 분해 기전의 근거.

9. Zhang S, et al. **Identification of the molecular basis of doxorubicin-
   induced cardiotoxicity.** *Nat Med.* 2012. — 심근 특이적 Top2b 결손 마우스가
   심장독성에 저항; 본 모델에서 Top2b를 독립 상태변수로 둔 직접적 근거.
   <https://pubmed.ncbi.nlm.nih.gov/23104132/>
10. Lyu YL, et al. **Topoisomerase IIβ mediated DNA double-strand breaks:
    implications in doxorubicin cardiotoxicity and prevention by dexrazoxane.**
    *Cancer Res.* 2007. <https://pubmed.ncbi.nlm.nih.gov/17875725/>
11. Vejpongsa P, Yeh ETH. **Topoisomerase 2β: a promising molecular target for
    primary prevention of anthracycline-induced cardiotoxicity.**
    *Clin Pharmacol Ther.* 2014. <https://pubmed.ncbi.nlm.nih.gov/24091715/>
12. Deng S, et al. **Dexrazoxane may prevent doxorubicin-induced DNA damage via
    depleting both topoisomerase II isoforms.** *BMC Cancer.* 2014. — 모델에서
    덱스라족산 효과를 철 킬레이션이 아니라 Top2b 단백 분해로 구현한 근거.
    <https://pubmed.ncbi.nlm.nih.gov/25406834/>
13. Martin E, et al. **Evaluation of the topoisomerase II-inactive
    bisdioxopiperazine ICRF-161 as a protectant against doxorubicin-induced
    cardiomyopathy.** *Toxicology.* 2009. — 철 킬레이션만 하고 Top2를
    억제하지 않는 유사체는 보호효과가 약하다 → 기전 분리의 실험적 근거.
    <https://pubmed.ncbi.nlm.nih.gov/19010377/>
14. Sawyer DB, et al. **Mechanisms of anthracycline cardiac injury: can we
    identify strategies for cardioprotection?** *Prog Cardiovasc Dis.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20728697/>

## 3. 철 · ROS · 미토콘드리아 · 페롭토시스 (Iron-redox arm) — AUC 구동 축의 근거

모델의 `TOXR → ROS → MITOD → 사멸` 축, LIP 되먹임, 준임계 루프 이득의 근거.

15. Ichikawa Y, et al. **Cardiotoxicity of doxorubicin is mediated through
    mitochondrial iron accumulation.** *J Clin Invest.* 2014. — 미토콘드리아
    철 축적이 핵심; 모델의 LIP 풀과 FEAMP 항.
    <https://pubmed.ncbi.nlm.nih.gov/24382354/>
16. Wallace KB, Sardão VA, Oliveira PJ. **Mitochondrial determinants of
    doxorubicin-induced cardiomyopathy.** *Circ Res.* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/32213135/>
17. Fang X, et al. **Ferroptosis as a target for protection against
    cardiomyopathy.** *Proc Natl Acad Sci USA.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30692261/>
18. Tadokoro T, et al. **Mitochondria-dependent ferroptosis plays a pivotal
    role in doxorubicin cardiotoxicity.** *JCI Insight.* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/36946465/>
19. Minotti G, et al. **Anthracyclines: molecular advances and pharmacologic
    developments in antitumor activity and cardiotoxicity.**
    *Pharmacol Rev.* 2004. — 퀴논 산화환원 순환·대사체·철 화학의 표준 리뷰.
    <https://pubmed.ncbi.nlm.nih.gov/15169927/>

## 4. 대사 · 독소루비시놀 · 약물유전체 (Metabolism, doxorubicinol, pharmacogenomics)

모델에서 대사체를 별도 심근 풀(CHM, t½ 28 d)로 둔 근거와 FM 변이(A7 분석).

20. Olson RD, et al. **Doxorubicin cardiotoxicity may be caused by its
    metabolite, doxorubicinol.** *Proc Natl Acad Sci USA.* 1988. — 대사체
    가설의 원전; 모델의 WM 가중과 CHM 축적.
    <https://pubmed.ncbi.nlm.nih.gov/2897122/>
21. Forrest GL, et al. **Human carbonyl reductase overexpression in the heart
    advances the development of doxorubicin-induced cardiotoxicity in
    transgenic mice.** *Cancer Res.* 2000. — CBR 과발현이 독성을 앞당김
    → FM 파라미터의 인과적 근거. <https://pubmed.ncbi.nlm.nih.gov/11016643/>
22. Blanco JG, et al. **Anthracycline-related cardiomyopathy after childhood
    cancer: role of polymorphisms in carbonyl reductase genes.**
    *J Clin Oncol.* 2012. <https://pubmed.ncbi.nlm.nih.gov/22124095/>
23. Aminkeng F, et al. **A coding variant in RARG confers susceptibility to
    anthracycline-induced cardiotoxicity in childhood cancer.**
    *Nat Genet.* 2015. <https://pubmed.ncbi.nlm.nih.gov/26237429/>
24. Visscher H, et al. **Pharmacogenomic prediction of anthracycline-induced
    cardiotoxicity in children.** *J Clin Oncol.* 2012. — SLC28A3 등 수송체
    변이. <https://pubmed.ncbi.nlm.nih.gov/21900104/>

## 5. 약동학 (Pharmacokinetics)

모델 PK 블록(3-구획 + 대사체 + 리포조말 담체)의 파라미터 근거.

25. Speth PAJ, et al. **Clinical pharmacokinetics of doxorubicin.**
    *Clin Pharmacokinet.* 1988. — CL·Vss·종말 t½ anchor.
    <https://pubmed.ncbi.nlm.nih.gov/3042244/>
26. Callies S, et al. **A population pharmacokinetic model for doxorubicin and
    doxorubicinol in patients.** *Cancer Chemother Pharmacol.* 2003. — 모약물-
    대사체 동시 모델; 본 모델 FM·CLM·VM의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/12647011/>
27. Gabizon A, Shmeeda H, Barenholz Y. **Pharmacokinetics of pegylated
    liposomal doxorubicin: review of animal and human studies.**
    *Clin Pharmacokinet.* 2003. — PLD의 긴 순환 반감기와 낮은 유리약물 첨두;
    모델의 KEL_LIP·KREL_LIP. <https://pubmed.ncbi.nlm.nih.gov/12739982/>

## 6. 제형 · 투여 스케줄 (Formulation and schedule) — 첨두 가설의 임상 검증

모델 A3(노출 지표 비교)의 핵심 근거. 누적용량이 동일한데 결과가 다르다는
관찰은, 손상축이 AUC가 아닌 첨두를 읽는다는 주장의 임상적 대응물이다.

28. Legha SS, et al. **Reduction of doxorubicin cardiotoxicity by prolonged
    continuous intravenous infusion.** *Ann Intern Med.* 1982. — 지속주입의
    원전. <https://pubmed.ncbi.nlm.nih.gov/7059060/>
29. van Dalen EC, et al. **Different dosage schedules for reducing
    cardiotoxicity in people with cancer receiving anthracycline
    chemotherapy.** *Cochrane Database Syst Rev.* 2016. — 지속주입 심부전
    RR≈0.27; 모델 A3/A8의 anchor.
    <https://pubmed.ncbi.nlm.nih.gov/26938118/>
30. O'Brien MER, et al. **Reduced cardiotoxicity and comparable efficacy in a
    phase III trial of pegylated liposomal doxorubicin versus conventional
    doxorubicin for first-line treatment of metastatic breast cancer.**
    *Ann Oncol.* 2004. <https://pubmed.ncbi.nlm.nih.gov/14998846/>
31. Batist G, et al. **Reduced cardiotoxicity and preserved antitumor efficacy
    of liposome-encapsulated doxorubicin and cyclophosphamide compared with
    conventional doxorubicin in metastatic breast cancer.** *J Clin Oncol.*
    2001. — 항종양 효과 유지 + 심장독성 감소 = 치료지수 분리의 임상 증거.
    <https://pubmed.ncbi.nlm.nih.gov/11230490/>

## 7. 덱스라족산 (Dexrazoxane)

모델 A4(용량 절감 등가), 결과 5·6의 근거.

32. de Baat EC, et al. **Dexrazoxane for preventing or reducing cardiotoxicity
    in adults and children with cancer receiving anthracyclines.**
    *Cochrane Database Syst Rev.* 2022. — 심부전 RR≈0.29; 모델 anchor.
    <https://pubmed.ncbi.nlm.nih.gov/36162822/>
33. Swain SM, et al. **Delayed administration of dexrazoxane provides
    cardioprotection for patients with advanced breast cancer treated with
    doxorubicin-containing therapy.** *J Clin Oncol.* 1997.
    <https://pubmed.ncbi.nlm.nih.gov/9193324/>
34. Lipshultz SE, et al. **Assessment of dexrazoxane as a cardioprotectant in
    doxorubicin-treated children with high-risk acute lymphoblastic
    leukaemia: long-term follow-up of a prospective, randomised, multicentre
    trial.** *Lancet Oncol.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20850381/>

## 8. 바이오마커 · 변형률 영상 (Biomarkers and strain imaging)

모델 A2(선행시간)와 TNI/BNP/GLS 관측식의 근거.

35. Cardinale D, et al. **Prognostic value of troponin I in cardiac risk
    stratification of cancer patients undergoing high-dose chemotherapy.**
    *Circulation.* 2004. — 트로포닌의 예측력, 특히 높은 음성예측도.
    <https://pubmed.ncbi.nlm.nih.gov/15148277/>
36. Cardinale D, et al. **Prevention of high-dose chemotherapy-induced
    cardiotoxicity in high-risk patients by enalapril.** *Circulation.* 2006.
    — 트로포닌 유도 조기개입 전략. <https://pubmed.ncbi.nlm.nih.gov/17101852/>
37. Lipshultz SE, et al. **Predictive value of cardiac troponin T in pediatric
    patients at risk for myocardial injury.** *Circulation.* 1997.
    <https://pubmed.ncbi.nlm.nih.gov/9355905/>
38. Ky B, et al. **Early increases in multiple biomarkers predict subsequent
    cardiotoxicity in patients with breast cancer treated with doxorubicin,
    taxanes, and trastuzumab.** *J Am Coll Cardiol.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24291281/>
39. Oikonomou EK, et al. **Assessment of prognostic value of left ventricular
    global longitudinal strain for early prediction of chemotherapy-induced
    cardiotoxicity: a systematic review and meta-analysis.**
    *JAMA Cardiol.* 2019. — GLS 상대 감소의 예측력; 모델 GLS 가중치 근거.
    <https://pubmed.ncbi.nlm.nih.gov/31433450/>
40. Negishi K, et al. **Independent and incremental value of deformation
    indices for prediction of trastuzumab-induced cardiotoxicity.**
    *J Am Soc Echocardiogr.* 2013. — GLS 상대 15% 역치의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/23562088/>
41. Cardinale D, et al. **Trastuzumab-induced cardiotoxicity: clinical and
    prognostic implications of troponin I evaluation.** *J Clin Oncol.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20679614/>

## 9. 심장보호 약물 임상시험 (Cardioprotection trials)

모델 A9(상류 vs 하류)와 결과 7·8의 근거. 모델은 스타틴을 ROS 생성 차단
(노출 중에만 유효), ACEi/BB/ARNI를 가역 결손·섬유화 축(노출 후에도 유효)에
배치하며, 그 구분은 아래 시험들의 상반된 결과를 설명하기 위한 것이다.

42. Neilan TG, et al. **Atorvastatin for anthracycline-associated cardiac
    dysfunction: the STOP-CA randomized clinical trial.** *JAMA.* 2023.
    — ≥10% LVEF 감소 9% vs 22%; 모델 A9 anchor.
    <https://pubmed.ncbi.nlm.nih.gov/37552303/>
43. Acar Z, et al. **Efficiency of atorvastatin in the protection of
    anthracycline-induced cardiomyopathy.** *J Am Coll Cardiol.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21851890/>
44. Bosch X, et al. **Enalapril and carvedilol for preventing chemotherapy-
    induced left ventricular systolic dysfunction (OVERCOME).**
    *J Am Coll Cardiol.* 2013. <https://pubmed.ncbi.nlm.nih.gov/23583763/>
45. Heck SL, et al. **Prevention of cardiac dysfunction during adjuvant breast
    cancer therapy (PRADA): extended follow-up of a 2×2 factorial,
    randomized, placebo-controlled, double-blind clinical trial of
    candesartan and metoprolol.** *Circulation.* 2021. — 효과가 작고 지속성이
    제한적이라는 점이 모델의 보수적 하류 가중치 근거.
    <https://pubmed.ncbi.nlm.nih.gov/33993702/>
46. Avila MS, et al. **Carvedilol for prevention of chemotherapy-related
    cardiotoxicity: the CECCY trial.** *J Am Coll Cardiol.* 2018. — LVEF
    일차종점은 음성이나 트로포닌은 감소 → 모델에서 BB를 부분적 ROS 소거 +
    신경호르몬 축으로 나눈 근거. <https://pubmed.ncbi.nlm.nih.gov/29540327/>
47. Cardinale D, et al. **Anthracycline-induced cardiotoxicity: a multicenter
    randomised trial comparing two strategies for guiding prevention with
    enalapril (ICOS-ONE).** *Eur J Cancer.* 2018. — 예방적 투여와 트로포닌
    유도 투여가 동등. <https://pubmed.ncbi.nlm.nih.gov/29567630/>
48. Livi L, et al. **Cardioprotective strategy for patients with nonmetastatic
    breast cancer who are receiving an anthracycline-based chemotherapy: the
    SAFE randomized clinical trial.** *JAMA Oncol.* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34436523/>
49. Gregorietti V, et al. **Use of sacubitril/valsartan in patients with
    cardiotoxicity and heart failure due to chemotherapy.**
    *Cardiooncology.* 2020. — ARNI 구제요법; 모델의 최대 회복 레버.
    <https://pubmed.ncbi.nlm.nih.gov/33292750/>
50. Quagliariello V, et al. **The SGLT-2 inhibitor empagliflozin improves
    myocardial strain, reduces cardiac fibrosis and pro-inflammatory
    cytokines in non-diabetic mice treated with doxorubicin.**
    *Cardiovasc Diabetol.* 2021. <https://pubmed.ncbi.nlm.nih.gov/34301253/>

## 10. 트라스투주맙 · ErbB2 (HER2 blockade interaction)

모델 A6(상호작용)과 결과 9의 근거. ErbB2 차단을 "복구능 결손"으로 구현하고,
전통적 Type I/II 이분법을 두 결손의 비율로 재해석하였다.

51. Slamon DJ, et al. **Use of chemotherapy plus a monoclonal antibody against
    HER2 for metastatic breast cancer that overexpresses HER2.**
    *N Engl J Med.* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Slamon+chemotherapy+plus+monoclonal+antibody+against+HER2+metastatic+breast+cancer+overexpresses+HER2+2001>
52. Romond EH, et al. **Trastuzumab plus adjuvant chemotherapy for operable
    HER2-positive breast cancer.** *N Engl J Med.* 2005. — NSABP B-31/N9831;
    동시 투여 심부전율. <https://pubmed.ncbi.nlm.nih.gov/16236738/>
53. Perez EA, et al. **Trastuzumab plus adjuvant chemotherapy for human
    epidermal growth factor receptor 2-positive breast cancer: planned joint
    analysis of overall survival from NSABP B-31 and NCCTG N9831.**
    *J Clin Oncol.* 2014. <https://pubmed.ncbi.nlm.nih.gov/25332249/>
54. Crone SA, et al. **ErbB2 is essential in the prevention of dilated
    cardiomyopathy.** *Nat Med.* 2002. — ErbB2 생존 신호의 필수성; 모델의
    WTRD·WTR_REP·WTR_REG. <https://pubmed.ncbi.nlm.nih.gov/11984589/>
55. Ewer MS, Lippman SM. **Type II chemotherapy-related cardiac dysfunction:
    time to recognize a new entity.** *J Clin Oncol.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15860848/>
56. Telli ML, et al. **Trastuzumab-related cardiotoxicity: calling into
    question the concept of reversibility.** *J Clin Oncol.* 2007. — Type II의
    "완전 가역성" 반박; 모델이 가역/비가역을 이분법이 아닌 비율로 둔 근거.
    <https://pubmed.ncbi.nlm.nih.gov/17687157/>

## 11. 심근세포 재생 (Cardiomyocyte renewal) — 비가역성의 정량적 근거

모델에서 KREG = 2×10⁻⁵/day (≈0.7 %/yr)로 둔 직접적 근거.

57. Bergmann O, et al. **Evidence for cardiomyocyte renewal in humans.**
    *Science.* 2009. — ¹⁴C 연대측정; 성체 심근세포 연간 갱신율 ~1%.
    <https://pubmed.ncbi.nlm.nih.gov/19342590/>
58. Bergmann O, et al. **Dynamics of cell generation and turnover in the human
    heart.** *Cell.* 2015. <https://pubmed.ncbi.nlm.nih.gov/26073943/>

## 12. 진료지침 · 정의 · 감시 전략 (Guidelines, definitions, surveillance)

CTRCD 정의(LVEF 10점 하락 & <50%), GLS 상대 15% 역치, 감시 주기의 근거.

59. Lyon AR, et al. **2022 ESC Guidelines on cardio-oncology developed in
    collaboration with the European Hematology Association (EHA), the
    European Society for Therapeutic Radiology and Oncology (ESTRO) and the
    International Cardio-Oncology Society (IC-OS).** *Eur Heart J.* 2022.
    <https://pubmed.ncbi.nlm.nih.gov/36017575/>
60. Plana JC, et al. **Expert consensus for multimodality imaging evaluation of
    adult patients during and after cancer therapy: a report from the American
    Society of Echocardiography and the European Association of Cardiovascular
    Imaging.** 2014. — CTRCD 정의의 표준.
    <https://pubmed.ncbi.nlm.nih.gov/25239940/>
61. Curigliano G, et al. **Management of cardiac disease in cancer patients
    throughout oncological treatment: ESMO consensus recommendations.**
    *Ann Oncol.* 2020. <https://pubmed.ncbi.nlm.nih.gov/31959335/>

## 13. 전임상 · in vitro 모델 (Preclinical and in vitro systems)

개인차(KDROS·TOXN50의 개체간 변이) 구조의 근거.

62. Burridge PW, et al. **Human induced pluripotent stem cell-derived
    cardiomyocytes recapitulate the predilection of breast cancer patients to
    doxorubicin-induced cardiotoxicity.** *Nat Med.* 2016. — 환자별 감수성이
    세포 수준에서 재현됨 → 모델의 개체간 변이가 약동학이 아니라 손상
    감수성에 있다는 가정의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/27089514/>

---

## 문헌과 모델 구조의 대응 관계 (How the literature maps onto the model)

| 모델 구조 | 대응 문헌 |
|---|---|
| Top2b를 독립 상태변수로 둠 (첨두 구동, DSB → p53 기억) | 9, 10, 11, 14 |
| 덱스라족산 = Top2b 프로테아좀 분해 (철 킬레이션이 아님) | 12, 13, 32 |
| 느린 잔류 풀 + 대사체 축적 (AUC 구동 산화환원 축) | 15, 16, 19, 20, 21 |
| 대사체를 별도 심근 풀로 (t½ 28 d) · FM 변이 | 20, 21, 22, 26 |
| 두 축이 서로 다른 노출 지표를 읽음 (스케줄·제형 효과) | 28, 29, 30, 31 |
| 근세포 소실의 비가역성 (KREG ≈ 0.7 %/yr) | 57, 58 |
| 가역적 기능 결손 + 가역성 창 (섬유화 시계) | 4, 36, 49 |
| 보상성 비대에 의한 LVEF 가면 → GLS·트로포닌 선행 | 35, 37, 38, 39, 40 |
| 상류(스타틴) vs 하류(ACEi/BB/ARNI) 비대칭 | 42, 43, 44, 45, 46, 47, 48, 49 |
| ErbB2 차단 = 복구능 결손 (초가산적 상호작용) | 51, 52, 53, 54 |
| Type I/II 이분법의 재해석 (두 결손의 비율) | 55, 56 |
| CTRCD 정의 및 감시 | 59, 60, 61 |
| 개체간 변이를 손상 감수성에 배치 | 23, 24, 62 |

## 모델이 문헌보다 낙관적인 부분 (Where this model overstates benefit)

정직한 기록을 위해 남긴다. 두 항목은 상한(upper bound)으로 취급해야 한다.

- **페길화 리포조말 독소루비신**: 모델의 CTRCD 상대위험 ≈ 0.08 이나 Cochrane
  메타분석의 임상 심부전 RR은 0.20 이다(문헌 29·30·31). 모델에서 유리약물
  첨두가 핵 손상축의 지배적 구동변수이므로 캡슐화의 이득이 과대평가된다.
- **보호제 병용(덱스라족산 + 스타틴 + ACEi/BB)**: 모델은 480 mg/m²에서 CTRCD를
  거의 0까지 낮추지만, 이 조합을 검증한 무작위 시험은 존재하지 않는다.
