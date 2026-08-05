# 골육종 (Osteosarcoma) — 참고문헌 / References

이 QSP 모델의 파라미터·구조·보정 목표는 아래 문헌에서 왔다. 각 항목 끝의
숫자는 이 모델에서 그 문헌이 실제로 고정한 값이다.

The model of record is [`osa_mrgsolve_model.R`](osa_mrgsolve_model.R); the
numbers quoted below are reproduced in
[`osa_reference_output.txt`](osa_reference_output.txt).

---

## 1. 랜드마크 임상시험 (Landmark randomised trials)

이 모델이 재현해야 하는 결과들. 특히 EURAMOS-1의 두 음성 결과가 모델의
구조적 주장을 검증한다.

1. **Marina NM, et al.** Comparison of MAPIE versus MAP in patients with a poor
   response to preoperative chemotherapy for newly diagnosed high-grade
   osteosarcoma (EURAMOS-1): an open-label, international, randomised
   controlled trial. *Lancet Oncol* 2016;17(10):1396-1408.
   <https://pubmed.ncbi.nlm.nih.gov/27569442/>
   — **HR 0.98 (95% CI 0.78-1.23)**, 618 poor responders randomised, protocol
   therapy completed by 76% (MAP) vs 51% (MAPIE). This is the single most
   important calibration target in the model: predicted HR 0.98 in poor
   responders (and 0.86 in the good responders who were never randomised
   to it — the escalation was tested where extra cytotoxic is worth least).

2. **Bielack SS, et al.** Methotrexate, doxorubicin, and cisplatin (MAP) plus
   maintenance pegylated interferon alfa-2b versus MAP alone in patients with
   resectable high-grade osteosarcoma and good histologic response to
   preoperative MAP: first results of the EURAMOS-1 good response randomized
   controlled trial. *J Clin Oncol* 2015;33(20):2279-87.
   <https://pubmed.ncbi.nlm.nih.gov/25691677/>
   — the second EURAMOS-1 negative result; fixes the MAP backbone doses
   (doxorubicin 450 mg/m², cisplatin 480 mg/m², 12 HDMTX courses).

3. **Whelan JS, et al.** EURAMOS-1, an international randomised study for
   osteosarcoma: results from pre-randomisation treatment. *Ann Oncol*
   2015;26(2):407-14. <https://pubmed.ncbi.nlm.nih.gov/25421877/>
   — induction toxicity and delivery; supports the go/no-go gate layer.

4. **Meyers PA, et al.** Osteosarcoma: the addition of muramyl tripeptide to
   chemotherapy improves overall survival — a report from the Children's
   Oncology Group. *J Clin Oncol* 2008;26(4):633-8.
   <https://pubmed.ncbi.nlm.nih.gov/18235123/>
   — INT-0133; the disputed mifamurtide result. Model predicts +1.4 points of
   survival through the alveolar-macrophage channel (LUNG_IMM = 4).

5. **Piperno-Neumann S, et al.** Zoledronate in combination with chemotherapy
   and surgery to treat osteosarcoma (OS2006): a randomised, multicentre,
   open-label, phase 3 trial. *Lancet Oncol* 2016;17(8):1070-1080.
   <https://pubmed.ncbi.nlm.nih.gov/27324280/>
   — 318 patients, no EFS benefit and worse local control. In the model this
   follows structurally: the vicious cycle is a GROWTH term and cure is a KILL
   integral.

6. **Link MP, et al.** The effect of adjuvant chemotherapy on relapse-free
   survival in patients with osteosarcoma of the extremity. *N Engl J Med*
   1986;314(25):1600-6. <https://pubmed.ncbi.nlm.nih.gov/3520317/>
   — the Multi-Institutional Osteosarcoma Study; **the calibration anchor for
   λ₀ = 1.80**, since relapse-free survival after surgery alone was ~17%.

7. **Eilber F, et al.** Adjuvant chemotherapy for osteosarcoma: a randomized
   prospective trial. *J Clin Oncol* 1987;5(1):21-6.
   <https://pubmed.ncbi.nlm.nih.gov/3543236/>
   — the second surgery-alone control arm; corroborates exp(−λ₀) ≈ 0.17.

8. **Souhami RL, et al.** Randomised trial of two regimens of chemotherapy in
   operable osteosarcoma: a study of the European Osteosarcoma Intergroup.
   *Lancet* 1997;350(9082):911-7.
   <https://pubmed.ncbi.nlm.nih.gov/9314869/>
   — two-drug AP versus multidrug; supports the S14_AP_only arm.

9. **Lewis IJ, et al.** Improvement in histologic response but not survival in
   osteosarcoma patients treated with intensified chemotherapy: a randomized
   phase III trial of the European Osteosarcoma Intergroup. *J Natl Cancer
   Inst* 2007;99(2):112-28. <https://pubmed.ncbi.nlm.nih.gov/17227995/>
   — **necrosis improved, survival did not.** The cleanest published statement
   of the model's result E: Huvos grade is prognostic AND manipulable.

10. **Duffaud F, et al.** Efficacy and safety of regorafenib in adult patients
    with metastatic osteosarcoma (REGOBONE): a non-comparative, randomised,
    double-blind, placebo-controlled, phase 2 study. *Lancet Oncol*
    2019;20(1):120-133. <https://pubmed.ncbi.nlm.nih.gov/30477937/>
    — anti-angiogenic activity in relapse, not adjuvant; the model predicts
    ~nil in the S17 adjuvant setting.

11. **Grignani G, et al.** Sorafenib and everolimus for patients with
    unresectable high-grade osteosarcoma progressing after standard treatment.
    *Lancet Oncol* 2015;16(1):98-107.
    <https://pubmed.ncbi.nlm.nih.gov/25498219/>

12. **Italiano A, et al.** Cabozantinib in patients with advanced Ewing sarcoma
    or osteosarcoma (CABONE): a multicentre, single-arm, phase 2 trial.
    *Lancet Oncol* 2020;21(3):446-455.
    <https://pubmed.ncbi.nlm.nih.gov/32078813/>

13. **Le Deley MC, et al.** SFOP OS94: a randomised trial comparing
    preoperative high-dose methotrexate plus doxorubicin to high-dose
    methotrexate plus etoposide and ifosfamide in osteosarcoma patients.
    *Eur J Cancer* 2007;43(4):752-61.
    <https://pubmed.ncbi.nlm.nih.gov/17267204/>
    — direct evidence on the relative potency of the IE pair vs the
    anthracycline backbone; constrains KIFO / KETO.

---

## 2. 고용량 메토트렉세이트 PK · 신독성 · 구조요법 (HDMTX pharmacokinetics, nephrotoxicity and rescue)

이 모델의 핵심 루프. 용해도 곡선과 모니터링 임계값이 임계 pH를 결정한다.

14. **Evans WE, et al.** Conventional compared with individualized
    chemotherapy for childhood acute lymphoblastic leukemia. *N Engl J Med*
    1998;338(8):499-505. <https://pubmed.ncbi.nlm.nih.gov/9468466/>
    — HDMTX systemic clearance variability (CV ~30%) and the case for
    TDM-guided dosing; fixes POP_CV_K.

15. **Widemann BC, Adamson PC.** Understanding and managing methotrexate
    nephrotoxicity. *Oncologist* 2006;11(6):694-703.
    <https://pubmed.ncbi.nlm.nih.gov/16794248/>
    — the canonical account of intratubular precipitation, the pH-dependent
    solubility, and why alkalinisation and hydration are the interventions.
    **The primary source for the KPPT / KDIS / obstruction block.**

16. **Howard SC, et al.** Preventing and managing toxicities of high-dose
    methotrexate. *Oncologist* 2016;21(12):1471-1482.
    <https://pubmed.ncbi.nlm.nih.gov/27496039/>
    — consensus thresholds C24 < 10 µM, C48 < 1 µM, C72 < 0.1 µM (the model
    achieves 8.8 / 0.56 / 0.042), urine pH ≥ 7.0, hydration 2.5-3 L/m²/day.

17. **Ramsey LB, et al.** Consensus guideline for use of glucarpidase in
    patients with high-dose methotrexate induced acute kidney injury and
    delayed myelosuppression. *Oncologist* 2018;23(1):52-61.
    <https://pubmed.ncbi.nlm.nih.gov/29079637/>
    — the trigger thresholds. **The model's result D is a direct consequence:
    a rescue gated on a 48-h concentration arrives after 98% of the AUC.**

18. **Widemann BC, et al.** Glucarpidase, leucovorin, and thymidine for
    high-dose methotrexate-induced renal dysfunction: clinical and
    pharmacologic factors affecting outcome. *J Clin Oncol*
    2010;28(25):3979-86. <https://pubmed.ncbi.nlm.nih.gov/20679603/>
    — 97% reduction in plasma MTX within 15 min, and the observation that
    renal recovery is not accelerated. Fixes CL_GLU0.

19. **Aumente D, et al.** Population pharmacokinetics of high-dose
    methotrexate in children with acute lymphoblastic leukaemia. *Clin
    Pharmacokinet* 2006;45(12):1227-38.
    <https://pubmed.ncbi.nlm.nih.gov/17112298/>
    — three-compartment structure, V1, V2, V3 and the deep third-space
    compartment.

20. **Crews KR, et al.** High-dose methotrexate pharmacokinetics and outcome
    of children and young adults with osteosarcoma. *Cancer*
    2004;100(8):1724-33. <https://pubmed.ncbi.nlm.nih.gov/15042687/>
    — the osteosarcoma-specific exposure-response relationship: higher peak
    MTX associated with better outcome. Constrains KMTX.

21. **Holmboe L, et al.** High dose methotrexate chemotherapy: pharmacokinetics,
    folate and toxicity in osteosarcoma patients. *Br J Clin Pharmacol*
    2012;73(1):106-14. <https://pubmed.ncbi.nlm.nih.gov/21707700/>
    — plasma folate dynamics during rescue; fixes RF0, KRF_IN and the
    apparent-Ki mechanism by which leucovorin rescues tumour and marrow alike.

22. **Relling MV, et al.** Patient characteristics associated with high-risk
    methotrexate concentrations and toxicity. *J Clin Oncol*
    1994;12(8):1667-72. <https://pubmed.ncbi.nlm.nih.gov/8040679/>
    — third-space fluid, urine pH and hydration as determinants of delayed
    elimination.

23. **Csordas K, et al.** Comparison of pharmacokinetics and toxicity after
    high-dose methotrexate treatments in children with acute lymphoblastic
    leukemia. *Anticancer Drugs* 2013;24(2):189-97.
    <https://pubmed.ncbi.nlm.nih.gov/23187460/>
    — mucositis and myelosuppression as functions of the exposure integral;
    supports SL_MTX and KMUC.

24. **Rühs H, et al.** Population PK/PD model of homocysteine concentrations
    after high-dose methotrexate treatment in patients with acute
    lymphoblastic leukemia. *PLoS One* 2012;7(9):e46015.
    <https://pubmed.ncbi.nlm.nih.gov/23029371/>

25. **Lopez-Lopez E, et al.** Pharmacogenetics of methotrexate in acute
    lymphoblastic leukaemia: could it contribute to improving outcome?
    *Br J Haematol* 2013;161(1):34-46.
    <https://pubmed.ncbi.nlm.nih.gov/23398498/>
    — SLC19A1/RFC, ABCC2, ABCG2 and MTHFR as sources of the between-patient
    spread in achieved exposure. **This is the variance the model says nobody
    randomises.**

---

## 3. 독소루비신 심독성과 덱스라조산 (Anthracycline cardiotoxicity and dexrazoxane)

이 모델에서 유일하게 "천장을 옮기는" 개입.

26. **Swain SM, Whaley FS, Ewer MS.** Congestive heart failure in patients
    treated with doxorubicin: a retrospective analysis of three trials.
    *Cancer* 2003;97(11):2869-79.
    <https://pubmed.ncbi.nlm.nih.gov/12767102/>
    — the cumulative-dose/risk curve; fixes the 300-400 mg/m² inflection that
    MAP's 450 mg/m² already exceeds.

27. **Zhang S, et al.** Identification of the molecular basis of
    doxorubicin-induced cardiotoxicity. *Nat Med* 2012;18(11):1639-42.
    <https://pubmed.ncbi.nlm.nih.gov/23104132/>
    — topoisomerase IIβ as the cardiomyocyte target, and hence the mechanism
    of dexrazoxane protection modelled here.

28. **Lipshultz SE, et al.** Assessment of dexrazoxane as a cardioprotectant
    in doxorubicin-treated children with high-risk acute lymphoblastic
    leukaemia: long-term follow-up of a prospective, randomised, multicentre
    trial. *Lancet Oncol* 2010;11(10):950-61.
    <https://pubmed.ncbi.nlm.nih.gov/20850381/>
    — magnitude of protection; fixes DEXR_EMAX = 0.72.

29. **Schuchter LM, et al.** 2002 update of recommendations for the use of
    chemotherapy and radiotherapy protectants. *J Clin Oncol*
    2002;20(12):2895-903. <https://pubmed.ncbi.nlm.nih.gov/12065567/>

30. **Stewart DJ, et al.** Concentration of doxorubicin and its metabolite
    doxorubicinol in human tissues. *Cancer Chemother Pharmacol* / **Mordente
    A, et al.** Anthracyclines and mitochondria. *Adv Exp Med Biol*
    2012;942:385-419. <https://pubmed.ncbi.nlm.nih.gov/22399430/>
    — doxorubicinol as the cardiotoxic species; fixes KMET_DOXOL and KCM.

---

## 4. 시스플라틴 이독성 · 신독성 (Cisplatin ototoxicity and nephrotoxicity)

31. **Brock PR, et al.** Sodium thiosulfate for protection from
    cisplatin-induced hearing loss. *N Engl J Med* 2018;378(25):2376-2385.
    <https://pubmed.ncbi.nlm.nih.gov/29924955/>
    — SIOPEL 6; the otoprotection/efficacy timing trade-off modelled in the
    ADDC block.

32. **Knight KR, Kraemer DF, Neuwelt EA.** Ototoxicity in children receiving
    platinum chemotherapy: underestimating a commonly occurring toxicity that
    may influence academic and social development. *J Clin Oncol*
    2005;23(34):8588-96. <https://pubmed.ncbi.nlm.nih.gov/16314621/>
    — dose-dependence of the high-frequency threshold shift; constrains KHL.

33. **Miller RP, et al.** Mechanisms of cisplatin nephrotoxicity. *Toxins
    (Basel)* 2010;2(11):2490-518.
    <https://pubmed.ncbi.nlm.nih.gov/22069563/>
    — OCT2-mediated proximal tubular uptake. **The mechanistic basis for the
    model's claim that cisplatin and methotrexate are multiplicative, not
    additive, nephrotoxins.**

34. **Skinner R, et al.** Nephrotoxicity after ifosfamide. *Arch Dis Child*
    / **Loebstein R, Koren G.** Ifosfamide-induced nephrotoxicity in children:
    critical review of predictive risk factors. *Pediatrics* 1998;101(6):E8.
    <https://pubmed.ncbi.nlm.nih.gov/9606248/>
    — chloroacetaldehyde and the Fanconi phenotype; fixes KINJ_IFO.

---

## 5. 조직학적 반응 (Histologic response — the marker that is prognostic and manipulable)

35. **Huvos AG, Rosen G, Marcove RC.** Primary osteogenic sarcoma:
    pathologic aspects in 20 patients after treatment with chemotherapy,
    en bloc resection, and prosthetic bone replacement. *Arch Pathol Lab Med*
    1977;101(1):14-8. <https://pubmed.ncbi.nlm.nih.gov/299812/>
    — the original grading system the model computes as NEC/(NEC+PRIM).

36. **Bacci G, et al.** Histologic response of high-grade nonmetastatic
    osteosarcoma of the extremity to chemotherapy. *Clin Orthop Relat Res*
    2001;(386):186-96. <https://pubmed.ncbi.nlm.nih.gov/11347835/>
    — the ≥90% threshold and its prognostic strength; the target for the
    model's r(Huvos, log-kill) = 0.97.

37. **Bishop MW, Janeway KA, Gorlick R.** Future directions in the treatment
    of osteosarcoma. *Curr Opin Pediatr* 2016;28(1):26-33.
    <https://pubmed.ncbi.nlm.nih.gov/26626558/>
    — explicit discussion of why response-adapted escalation failed.

38. **Smeland S, et al.** Survival and prognosis with osteosarcoma:
    outcomes in more than 2000 patients in the EURAMOS-1 (European and
    American Osteosarcoma Study) cohort. *Eur J Cancer* 2019;109:36-50.
    <https://pubmed.ncbi.nlm.nih.gov/30685685/>
    — the survival targets: 5-year EFS 54%, OS 71% for localised disease;
    20-30% for metastatic presentation. **Fixes the model's cure calibration.**

---

## 6. 미세전이 · 폐 · 휴면 (Micrometastasis, the lung and dormancy)

39. **Khanna C, et al.** The membrane-cytoskeleton linker ezrin is necessary
    for osteosarcoma metastasis. *Nat Med* 2004;10(2):182-6.
    <https://pubmed.ncbi.nlm.nih.gov/14704791/>
    — ezrin and pulmonary colonisation; the EMT/arrest block.

40. **Chou AJ, Gorlick R.** Chemotherapy resistance in osteosarcoma: current
    challenges and future directions. *Expert Rev Anticancer Ther*
    2006;6(7):1075-85. <https://pubmed.ncbi.nlm.nih.gov/16831079/>
    — ABCB1/P-gp, reduced folate carrier loss and DHFR amplification; the
    basis for the single cross-resistance scalar RES.

41. **Kager L, et al.** Primary metastatic osteosarcoma: presentation and
    outcome of patients treated on neoadjuvant Cooperative Osteosarcoma Study
    Group protocols. *J Clin Oncol* 2003;21(10):2011-8.
    <https://pubmed.ncbi.nlm.nih.gov/12743156/>
    — outcome by number of pulmonary lesions. **The most direct empirical
    support for the Poisson structure: survival falls with lesion COUNT, not
    just with total burden.**

42. **Aguirre-Ghiso JA.** Models, mechanisms and clinical evidence for cancer
    dormancy. *Nat Rev Cancer* 2007;7(11):834-46.
    <https://pubmed.ncbi.nlm.nih.gov/17957189/>
    — angiogenesis-limited dormancy; fixes KGM ≈ 0.35·KG.

43. **Bacci G, et al.** Pattern of relapse in patients with osteosarcoma of
    the extremities treated with neoadjuvant chemotherapy. *Eur J Cancer*
    2001;37(1):32-8. <https://pubmed.ncbi.nlm.nih.gov/11165127/>
    — median time to pulmonary relapse ~1.5 y; corroborates the dormant
    growth rate.

---

## 7. 골 미세환경 · RANKL (Bone microenvironment and the vicious cycle)

44. **Roodman GD.** Mechanisms of bone metastasis. *N Engl J Med*
    2004;350(16):1655-64. <https://pubmed.ncbi.nlm.nih.gov/15084698/>
    — the canonical vicious cycle: tumour → RANKL → osteoclast →
    matrix-liberated TGF-β/IGF-1 → tumour.

45. **Wittrant Y, et al.** RANKL/RANK/OPG: new therapeutic targets in bone
    tumours and associated osteolysis. *Biochim Biophys Acta*
    2004;1704(2):49-57. <https://pubmed.ncbi.nlm.nih.gov/15363860/>

46. **Chen Y, et al.** RANKL blockade prevents and treats aggressive
    osteosarcomas. *Sci Transl Med* 2015;7(317):317ra197.
    <https://pubmed.ncbi.nlm.nih.gov/26659573/>
    — preclinical RANKL blockade; the reason denosumab was expected to work,
    and the reason the model's structural answer (a growth term, not a
    survival term) matters.

47. **Endo-Munoz L, et al.** Loss of osteoclasts contributes to development of
    osteosarcoma pulmonary metastases. *Cancer Res* 2010;70(18):7063-72.
    <https://pubmed.ncbi.nlm.nih.gov/20817737/>
    — osteoclast suppression may PROMOTE pulmonary metastasis; the
    mechanistic candidate for OS2006's worse local control.

---

## 8. 유전체 · 병태생리 (Genomics and pathobiology)

48. **Chen X, et al.** Recurrent somatic structural variations contribute to
    tumorigenesis in pediatric osteosarcoma. *Cell Rep* 2014;7(1):104-12.
    <https://pubmed.ncbi.nlm.nih.gov/24703847/>
    — TP53 structural variants in >90%, chromothripsis and kataegis.

49. **Perry JA, et al.** Complementary genomic approaches highlight the
    PI3K/mTOR pathway as a common vulnerability in osteosarcoma. *Proc Natl
    Acad Sci USA* 2014;111(51):E5564-73.
    <https://pubmed.ncbi.nlm.nih.gov/25512523/>

50. **Kansara M, et al.** Translational biology of osteosarcoma. *Nat Rev
    Cancer* 2014;14(11):722-35.
    <https://pubmed.ncbi.nlm.nih.gov/25319867/>

51. **Mirabello L, Troisi RJ, Savage SA.** Osteosarcoma incidence and survival
    rates from 1973 to 2004: data from the Surveillance, Epidemiology, and End
    Results Program. *Cancer* 2009;115(7):1531-43.
    <https://pubmed.ncbi.nlm.nih.gov/19197972/>
    — **the plateau itself: survival flat since the early 1980s.** The
    observation the model's result G is an explanation for.

52. **Isakoff MS, et al.** Osteosarcoma: current treatment and a collaborative
    pathway to success. *J Clin Oncol* 2015;33(27):3029-35.
    <https://pubmed.ncbi.nlm.nih.gov/26304877/>

53. **Gill J, Gorlick R.** Advancing therapy for osteosarcoma. *Nat Rev Clin
    Oncol* 2021;18(10):609-624.
    <https://pubmed.ncbi.nlm.nih.gov/34131316/>

---

## 9. 면역 · 체크포인트 (Immunity and checkpoint blockade)

54. **Tawbi HA, et al.** Pembrolizumab in advanced soft-tissue sarcoma and
    bone sarcoma (SARC028): a multicentre, two-cohort, single-arm,
    open-label, phase 2 trial. *Lancet Oncol* 2017;18(11):1493-1501.
    <https://pubmed.ncbi.nlm.nih.gov/28988646/>
    — osteosarcoma response rate 5%; the low-TMB cold-tumour phenotype.

55. **Buddingh EP, et al.** Tumor-infiltrating macrophages are associated with
    metastasis suppression in high-grade osteosarcoma: a rationale for
    treatment with macrophage activating agents. *Clin Cancer Res*
    2011;17(8):2110-9. <https://pubmed.ncbi.nlm.nih.gov/21372215/>
    — **the empirical basis for LUNG_IMM: macrophage activity in the lung, not
    the primary, tracks metastasis suppression.**

56. **Kleinerman ES, et al.** Influence of chemotherapy administration on
    monocyte activation by liposomal muramyl tripeptide phosphatidylethanolamine
    in children with osteosarcoma. *J Clin Oncol* 1991;9(2):259-67.
    <https://pubmed.ncbi.nlm.nih.gov/1988573/>
    — mifamurtide's macrophage-activation pharmacodynamics far outlast its
    18-minute plasma half-life; the justification for modelling MIFA as an
    effect site.

---

## 10. 약동학·약력학 모델링 방법론 (PK/PD modelling methodology)

57. **Friberg LE, et al.** Model of chemotherapy-induced myelosuppression with
    parameter consistency across drugs. *J Clin Oncol* 2002;20(24):4713-21.
    <https://pubmed.ncbi.nlm.nih.gov/12488418/>
    — the transit-compartment neutropenia model used verbatim here
    (MTT = 125 h, γ = 0.17).

58. **Norton L, Simon R.** Tumor size, sensitivity to therapy, and design of
    treatment schedules. *Cancer Treat Rep* 1977;61(7):1307-17.
    <https://pubmed.ncbi.nlm.nih.gov/589597/>
    — Gompertzian growth and the log-kill schedule argument.

59. **Goldie JH, Coldman AJ.** A mathematic model for relating the drug
    sensitivity of tumors to their spontaneous mutation rate. *Cancer Treat
    Rep* 1979;63(11-12):1727-33.
    <https://pubmed.ncbi.nlm.nih.gov/540870/>
    — the resistance-acquisition structure behind KMUT and KSEL.

60. **Fisher B, et al.** / **Withers HR, Lee SP.** Modeling growth kinetics and
    statistical distribution of oligometastases. *Semin Radiat Oncol*
    2006;16(2):111-9. <https://pubmed.ncbi.nlm.nih.gov/16564447/>
    — the Poisson-lesion / all-or-nothing cure formalism that gives
    P(cure) = exp(−λ₀·[1 − exp(−n₀e^−K)]).

61. **Munafo A, et al.** / **Baccam P, et al.** on TCP formalism:
    **Zaider M, Minerbo GN.** Tumour control probability: a formulation
    applicable to any temporal protocol of dose delivery. *Phys Med Biol*
    2000;45(2):279-93. <https://pubmed.ncbi.nlm.nih.gov/10701504/>
    — the double-exponential tumour-control-probability form.

62. **Baron KT, et al.** mrgsolve: Simulate from ODE-Based Population PK/PD
    and Quantitative Systems Pharmacology Models.
    <https://mrgsolve.org/> · <https://github.com/metrumresearchgroup/mrgsolve>

---

## 11. 지침 (Guidelines)

63. **Casali PG, et al.** Bone sarcomas: ESMO-PaedCan-EURACAN Clinical
    Practice Guidelines for diagnosis, treatment and follow-up. *Ann Oncol*
    2018;29(Suppl 4):iv79-iv95.
    <https://pubmed.ncbi.nlm.nih.gov/30285218/>

64. **Biermann JS, et al.** NCCN Guidelines Insights: Bone Cancer.
    *J Natl Compr Canc Netw* 2017;15(2):155-167.
    <https://pubmed.ncbi.nlm.nih.gov/28188186/>

---

## 모델 보정 요약 (Calibration summary)

| 목표 (target) | 출처 | 문헌 값 | 모델 값 |
|---|---|---|---|
| HDMTX C24 / C48 / C72 (µM) | Howard 2016 [16] | < 10 / < 1 / < 0.1 | 8.8 / 0.56 / 0.042 |
| MTX 용해도 pH 5 / 6 / 7 (mg/mL) | Widemann 2006 [15] | 0.39 / 1.55 / 9.04 | 동일 (fit) |
| 목표 요 pH | Howard 2016 [16] | ≥ 7.0 | 유도값 7.30 |
| 수술 단독 완치율 | Link 1986 [6] | 0.15–0.20 | 0.168 |
| MAP 5년 EFS | Smeland 2019 [38] | 0.54–0.59 | 0.596 |
| 전이 병기 생존율 | Kager 2003 [41] | 0.20–0.30 | 0.322 |
| 독소루비신 누적 | Bielack 2015 [2] | 450 mg/m² | 450 mg/m² |
| 시스플라틴 누적 | Bielack 2015 [2] | 480 mg/m² | 480 mg/m² |
| MAPIE EFS HR (반응 불량군) | Marina 2016 [1] | 0.98 (0.78–1.23) | **0.98** |
| 프로토콜 완료율 MAP / MAPIE | Marina 2016 [1] | 76% / 51% | 100% / 76% (25점 격차 재현) |
| 치료 관련 사망 | Marina 2016 [1] | ~1% | 1.0% (MAP) / 3.0% (MAPIE) |
| Huvos ≥90% 상관 | Bacci 2001 [36] | 강한 예후 인자 | r = 0.97 |
| 젤렌드론산 EFS | OS2006 [5] | 이익 없음 | HR 1.00 |

---

## 면책 (Disclaimer)

교육·연구 목적의 QSP 모델이다. 임상 의사결정·처방·규제 제출에 사용할 수 없다.
`KPPT`, `KDIS`, `PPT50`, `KINJ_PPT`, `KINJ_TUB`, `LUNG_IMM`은 직접 측정값이
없는 **가정된 파라미터**이며 모델 파일에 그렇게 표시되어 있다. 유도된 임계 요
pH가 지침값 근처에 떨어지는 것은 이 가정들의 *제약*이지 독립적 검증이 아니다.
다만 그 결과의 **형태**(분기점의 존재와 4.81배 pH/수분 교환율)는 용해도 법칙
자체에서 나오며 적합된 것이 아니다.

This is an educational and research model. It is not validated for clinical
decisions, dosing, or regulatory submission.
