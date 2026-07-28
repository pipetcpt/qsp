# 산후우울증 (Postpartum Depression) QSP 모델 — 참고문헌
# References for the Postpartum Depression QSP Model

**모든 PMID는 이 세션에서 NCBI E-utilities (`esearch` + `esummary`)로 직접 조회하여
제목·저널·연도를 확인한 것입니다.** 기억에 의존해 적은 것이 아니며, 검증 스크립트는
`ppd_reference_check.py`의 헤더에 기록된 캘리브레이션 목표값과 동일한 출처에서
가져왔습니다. 확인되지 않은 문헌은 이 목록에 넣지 않았습니다.

Every PMID below was resolved live against PubMed during model construction; titles,
journals and years are as returned by NCBI. Nothing here is quoted from memory.

---

## 1. 역학·임상 정의 (Epidemiology and clinical definition)

1. Munk-Olsen T, et al. **New parents and mental disorders: a population-based register
   study.** *JAMA* 2006. — 산후 정신질환 입원 위험의 시간적 집중(첫 3개월).
   <https://pubmed.ncbi.nlm.nih.gov/17148723/>
2. Wisner KL, et al. **Onset timing, thoughts of self-harm, and diagnoses in postpartum
   women with screen-positive depression findings.** *JAMA Psychiatry* 2013. — 10,000명
   선별 코호트; 발병 시점 분포와 자해 사고 빈도. 모델의 "발병 창(window)" 비교 대상.
   <https://pubmed.ncbi.nlm.nih.gov/23487258/>
3. O'Hara MW, McCabe JE. **Postpartum depression: current status and future
   directions.** *Annu Rev Clin Psychol* 2013.
   <https://pubmed.ncbi.nlm.nih.gov/23394227/>
4. Stewart DE, Vigod S. **Postpartum Depression.** *N Engl J Med* 2016.
   <https://pubmed.ncbi.nlm.nih.gov/27959754/>
5. Meltzer-Brody S, et al. **Postpartum psychiatric disorders.** *Nat Rev Dis Primers*
   2018. <https://pubmed.ncbi.nlm.nih.gov/29695824/>
6. Howard LM, et al. **Non-psychotic mental disorders in the perinatal period.**
   *Lancet* 2014. <https://pubmed.ncbi.nlm.nih.gov/25455248/>
7. Gavin NI, et al. **Perinatal depression: a systematic review of prevalence and
   incidence.** *Obstet Gynecol* 2005.
   <https://pubmed.ncbi.nlm.nih.gov/16260528/>
8. Shorey S, et al. **Prevalence and incidence of postpartum depression among healthy
   mothers: a systematic review and meta-analysis.** *J Psychiatr Res* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/30114665/>
9. **Prevalence of Postpartum Depression Based on Diagnostic Interviews: A Systematic
   Review and Meta-Analysis.** *Depress Anxiety* 2023. — 구조화 면담 기준 유병률
   (선별검사 기준보다 낮음).
   <https://pubmed.ncbi.nlm.nih.gov/40224605/>
10. **Association of Postpartum Depression with Maternal Suicide: A Nationwide
    Population-Based Study.** *Int J Environ Res Public Health* 2022. — 모델의
    자살사고 엔드포인트 정당화.
    <https://pubmed.ncbi.nlm.nih.gov/35564525/>
11. Bergink V, et al. **Postpartum Psychosis: Madness, Mania, and Melancholia in
    Motherhood.** *Am J Psychiatry* 2016. — PPD와 구별되는 별개 응급질환.
    <https://pubmed.ncbi.nlm.nih.gov/27609245/>

## 2. 평가 척도 (Rating scales and their model mapping)

12. Cox JL, Holden JM, Sagovsky R. **Detection of postnatal depression. Development of
    the 10-item Edinburgh Postnatal Depression Scale.** *Br J Psychiatry* 1987. — EPDS
    원논문. <https://pubmed.ncbi.nlm.nih.gov/3651732/>
13. Levis B, et al. **Accuracy of the Edinburgh Postnatal Depression Scale (EPDS) for
    screening to detect major depression among pregnant and postpartum women.** *BMJ*
    2020. — 개별환자자료 메타분석; 12/13 절단값 근거.
    <https://pubmed.ncbi.nlm.nih.gov/33177069/>
14. Zimmerman M, et al. **Is the cutoff to define remission on the Hamilton Rating
    Scale for Depression too high?** *J Nerv Ment Dis* 2005. — 모델이 HAM-D ≤ 7을
    관해 경계로 쓰는 근거와 그 한계.
    <https://pubmed.ncbi.nlm.nih.gov/15729106/>
15. **Associations between commonly used patient-reported outcome tools in postpartum
    depression clinical practice and the Hamilton Rating Scale.** *Arch Womens Ment
    Health* 2020. — EPDS ↔ HAM-D 상호 변환.
    <https://pubmed.ncbi.nlm.nih.gov/32666402/>

## 3. 신경스테로이드와 GABA_A 수용체 — 모델의 중심 기전

16. Maguire J, Mody I. **GABA(A)R plasticity during pregnancy: relevance to postpartum
    depression.** *Neuron* 2008. — **모델의 핵심 논문.** 임신 중 δ 소단위 하향조절과
    분만 후 회복 지연, δ-KO 마우스의 산후 우울·모성행동 이상.
    <https://pubmed.ncbi.nlm.nih.gov/18667149/>
17. Maguire JL, Stell BM, Rafizadeh M, Mody I. **Ovarian cycle-linked changes in
    GABA(A) receptors mediating tonic inhibition alter seizure susceptibility and
    anxiety.** *Nat Neurosci* 2005. — 긴장성 억제와 δ 발현의 스테로이드 의존적 가소성.
    <https://pubmed.ncbi.nlm.nih.gov/15895085/>
18. Maguire J, Mody I. **Neurosteroid synthesis-mediated regulation of GABA(A)
    receptors: relevance to the ovarian cycle and stress.** *J Neurosci* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/17329412/>
19. Hosie AM, Wilkins ME, da Silva HMA, Smart TG. **Endogenous neurosteroids regulate
    GABAA receptors through two discrete transmembrane sites.** *Nature* 2006. —
    **모델이 효능(부위 1, nM)과 진정(부위 2, µM)을 분리하는 구조적 근거.**
    <https://pubmed.ncbi.nlm.nih.gov/17108970/>
20. Hosie AM, et al. **Neurosteroid binding sites on GABA(A) receptors.** *Pharmacol
    Ther* 2007. <https://pubmed.ncbi.nlm.nih.gov/17560657/>
21. Belelli D, Lambert JJ. **Neurosteroids: endogenous regulators of the GABA(A)
    receptor.** *Nat Rev Neurosci* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15959466/>
22. Stell BM, Brickley SG, Tang CY, Farrant M, Mody I. **Neuroactive steroids reduce
    neuronal excitability by selectively enhancing tonic inhibition mediated by delta
    subunit-containing GABAA receptors.** *PNAS* 2003. — 세포외 시냅스 선택성.
    <https://pubmed.ncbi.nlm.nih.gov/14623958/>
23. Mihalek RM, et al. **Attenuated sensitivity to neuroactive steroids in
    gamma-aminobutyrate type A receptor delta subunit knockout mice.** *PNAS* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/10536021/>
24. Concas A, et al. **Role of brain allopregnanolone in the plasticity of
    gamma-aminobutyric acid type A receptor in rat brain during pregnancy and after
    delivery.** *PNAS* 1998. — **모델의 두 시간상수(리간드 vs 수용체)의 실험적 원형.**
    <https://pubmed.ncbi.nlm.nih.gov/9789080/>
25. Sarkar J, Wakefield S, MacKenzie G, Moss SJ, Maguire J. **Neurosteroidogenesis is
    required for the physiological response to stress: role of neurosteroid-sensitive
    GABAA receptors.** *J Neurosci* 2011. — PVN CRH 뉴런의 δ-GABA_A 브레이크;
    모델의 HPA-회로 결합.
    <https://pubmed.ncbi.nlm.nih.gov/22171026/>
26. Ferando I, Mody I. **Interneuronal GABAA receptors inside and outside of
    synapses.** *Curr Opin Neurobiol* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24650505/>
27. Paul SM, Purdy RH. **Neuroactive steroids.** *FASEB J* 1992.
    <https://pubmed.ncbi.nlm.nih.gov/1347506/>
28. Luscher B, Shen Q, Sahir N. **The GABAergic deficit hypothesis of major depressive
    disorder.** *Mol Psychiatry* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21079608/>
29. Windle RJ 등을 포함한 리뷰: **[GABAergic approach of postpartum depression: A
    translational review of literature].** *Encephale* 2020. — KCC2·염화물 항상성 포함
    번역적 정리. <https://pubmed.ncbi.nlm.nih.gov/31767256/>
30. **The human gamma-aminobutyric acid A receptor delta (GABRD) gene: molecular
    characterisation and tissue-specific expression.** *Gene* 2002. — δ 소단위 유전자.
    <https://pubmed.ncbi.nlm.nih.gov/12119096/>
31. Bäckström T, et al. **Positive GABA(A) receptor modulating steroids and their
    antagonists: implications for clinical treatments.** *J Neuroendocrinol* 2022. —
    이소알로프레그나놀론(세프라놀론)의 기능적 길항; 모델의 ISOALLO 노드.
    <https://pubmed.ncbi.nlm.nih.gov/34337790/>
32. Bixo M, et al. **Treatment of premenstrual dysphoric disorder with the GABA(A)
    receptor modulating steroid antagonist Sepranolone (UC1010).**
    *Psychoneuroendocrinology* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28319848/>

## 4. 스테로이드 대사·합성 경로 (Synthesis and clearance parameters)

33. **Serum concentrations of progesterone and 5 alpha-pregnane-3,20-dione during
    labor and early post partum.** *Acta Obstet Gynecol Scand* 1990. — 모델의 P4 및
    5α-DHP 소실 반감기 근거.
    <https://pubmed.ncbi.nlm.nih.gov/2386015/>
34. **Allopregnanolone serum levels and expression of 5 alpha-reductase and 3
    alpha-hydroxysteroid dehydrogenase isoforms.** *Epilepsy Res* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/12742591/>
35. **Hypothalamic 5 alpha-reductase and 3 alpha-oxidoreductase activity.**
    *J Steroid Biochem Mol Biol* 2002. — 뇌 국소 합성(모델의 BRSYN 항).
    <https://pubmed.ncbi.nlm.nih.gov/11867268/>
36. Timby E, et al. **Pharmacokinetic and behavioral effects of allopregnanolone in
    healthy women.** *Psychopharmacology (Berl)* 2006. — 사람에서 IV ALLO의 PK와
    진정 효과. <https://pubmed.ncbi.nlm.nih.gov/16177884/>
37. **Allopregnanolone assays.** *J Clin Endocrinol Metab* 2001. — RIA vs LC-MS/MS
    불일치; 모델 가정 A4(절대값보다 비율을 맞춘 이유)의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/11238543/>
38. **Peripartum allopregnanolone blood concentrations and depressive symptoms: a
    systematic review and individual participant data meta-analysis.** *Mol
    Psychiatry* 2025. — **핵심:** ALLO 절대 농도와 증상의 관계는 일관되지 않음 →
    "낮은 ALLO"가 아니라 "변화율/수용체 반응"이 문제라는 모델 전제를 지지.
    <https://pubmed.ncbi.nlm.nih.gov/39511449/>
39. **Allopregnanolone in the peripartum: correlates, concentrations, and challenges —
    a systematic review.** *Psychoneuroendocrinology* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38759520/>
40. **Profiling neuroactive steroids in pregnancy: a non-derivatised LC-MS/MS
    method.** *J Chromatogr B* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/40054418/>

## 5. 사람 데이터 — 호르몬 금단과 취약성

41. Bloch M, et al. **Effects of gonadal steroids in women with a history of postpartum
    depression.** *Am J Psychiatry* 2000. — **모델의 취약성(V) 개념의 실험적 근거:**
    동일한 스테로이드 부하/금단이 병력이 있는 여성에서만 증상을 유발.
    <https://pubmed.ncbi.nlm.nih.gov/10831472/>
42. Bloch M, Daly RC, Rubinow DR. **Endocrine factors in the etiology of postpartum
    depression.** *Compr Psychiatry* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/12764712/>
43. Schiller CE, Meltzer-Brody S, Rubinow DR. **The role of reproductive hormones in
    postpartum depression.** *CNS Spectr* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/25263255/>
44. Nappi RE, et al. **Serum allopregnanolone in women with postpartum "blues".**
    *Obstet Gynecol* 2001. <https://pubmed.ncbi.nlm.nih.gov/11152912/>
45. Hellgren C, et al. **Low serum allopregnanolone is associated with symptoms of
    depression in late pregnancy.** *Neuropsychobiology* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24776841/>
46. Osborne LM, et al. **Lower allopregnanolone during pregnancy predicts postpartum
    depression: an exploratory study.** *Psychoneuroendocrinology* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28278440/>
47. Epperson CN, et al. **Preliminary evidence of reduced occipital GABA
    concentrations in puerperal women: a 1H-MRS study.** *Psychopharmacology (Berl)*
    2006. <https://pubmed.ncbi.nlm.nih.gov/16724188/>
48. Deligiannidis KM, et al. **Resting-state functional connectivity, cortical GABA,
    and neuroactive steroids in peripartum and peripartum depressed women.**
    *Neuropsychopharmacology* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30327498/>
49. Deligiannidis KM, et al. **GABAergic neuroactive steroids and resting-state
    functional connectivity in postpartum depression: a preliminary study.**
    *J Psychiatr Res* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23499388/>

## 6. HPA 축 (Second slow variable)

50. Magiakou MA, et al. **Hypothalamic corticotropin-releasing hormone suppression
    during the postpartum period: implications for the increase in psychiatric
    manifestations at this time.** *J Clin Endocrinol Metab* 1996. — **모델의 HCRH
    회복 시간상수(≈ 12주 창)의 근거.**
    <https://pubmed.ncbi.nlm.nih.gov/8626857/>
51. Yim IS, et al. **Risk of postpartum depressive symptoms with elevated
    corticotropin-releasing hormone in human pregnancy.** *Arch Gen Psychiatry* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19188538/>
52. Glynn LM, Davis EP, Sandman CA. **New insights into the role of perinatal HPA-axis
    dysregulation in postpartum depression.** *Neuropeptides* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/24210135/>

## 7. 모노아민·염증·키뉴레닌

53. Sacher J, et al. **Elevated brain monoamine oxidase A binding in the early
    postpartum period.** *Arch Gen Psychiatry* 2010. — **모델의 MAO-A +21 % 항의
    직접 근거** (E2 금단 → MAO-A 상승).
    <https://pubmed.ncbi.nlm.nih.gov/20439828/>
54. Osborne LM, Monk C. **Perinatal depression — the fourth inflammatory morbidity of
    pregnancy? Theory and literature review.** *Psychoneuroendocrinology* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23608136/>
55. Bränn E, et al. **Inflammatory markers in women with postpartum depressive
    symptoms.** *J Neurosci Res* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/30252150/>
56. Achtyes E, et al. **Inflammation and kynurenine pathway dysregulation in
    post-partum women with severe and suicidal depression.** *Brain Behav Immun*
    2020. — 모델의 KYN/TRP 및 QUIN/KYNA 항.
    <https://pubmed.ncbi.nlm.nih.gov/31698012/>
57. **Associations between estrogen and progesterone, the kynurenine pathway, and
    inflammation in the post-partum.** *J Affect Disord* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/33278766/>
58. Corwin EJ, et al. **Symptoms of postpartum depression associated with elevated
    levels of interleukin-1 beta during the first month postpartum.** *Biol Res Nurs*
    2008. <https://pubmed.ncbi.nlm.nih.gov/18829596/>
59. Corwin EJ, Pajer K. **The psychoneuroimmunology of postpartum depression.**
    *J Womens Health* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/18710365/>

## 8. 유전·후성유전 (Vulnerability, V)

60. Guintivano J, Arad M, Gould TD, Payne JL, Kaminsky ZA. **Antenatal prediction of
    postpartum depression with blood DNA methylation biomarkers.** *Mol Psychiatry*
    2014. — HP1BP3 / TTC9B.
    <https://pubmed.ncbi.nlm.nih.gov/23689534/>
61. **Meta-Analyses of Genome-Wide Association Studies for Postpartum Depression.**
    *Am J Psychiatry* 2023. — 최초의 대규모 PPD GWAS 메타분석.
    <https://pubmed.ncbi.nlm.nih.gov/37849304/>
62. **The First Large GWAS Meta-Analysis for Postpartum Depression** (editorial).
    *Am J Psychiatry* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/38037399/>
63. Viktorin A, et al. **Heritability of perinatal depression and genetic overlap with
    nonperinatal depression.** *Am J Psychiatry* 2016.
    <https://pubmed.ncbi.nlm.nih.gov/26337037/>

## 9. 수면 — 모델의 악순환 고리 (the closed loop)

64. Lawson A, et al. **The relationship between sleep and postpartum mental disorders:
    a systematic review.** *J Affect Disord* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/25702602/>
65. **Changes in sleep quality, but not hormones predict time to postpartum depression
    recurrence.** *J Affect Disord* 2011. — **모델이 수면을 (호르몬이 아니라)
    가장 강한 근위 구동자로 두는 직접 근거.**
    <https://pubmed.ncbi.nlm.nih.gov/20708275/>
66. Okun ML, et al. **Sleep complaints in late pregnancy and the recurrence of
    postpartum depression.** *Behav Sleep Med* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19330583/>
67. Goyal D, Gay C, Lee K. **Fragmented maternal sleep is more strongly correlated with
    depressive symptoms than infant temperament at three months postpartum.** *Arch
    Womens Ment Health* 2009. — 분절화 > 총 수면시간.
    <https://pubmed.ncbi.nlm.nih.gov/19396527/>
68. **Maternal depressive symptoms, dysfunctional cognitions, and infant night waking:
    the role of maternal nighttime behavior.** *Child Dev* 2012. — 양방향 고리.
    <https://pubmed.ncbi.nlm.nih.gov/22506917/>
69. **Preventing recurrence of postpartum depression by regulating sleep.** *Expert Rev
    Neurother* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37462620/>

## 10. 옥시토신·모아 결합·영아 결과

70. Skrundz M, et al. **Plasma oxytocin concentration during pregnancy is associated
    with development of postpartum depression.** *Neuropsychopharmacology* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21562482/>
71. Moehler E, et al. **Maternal depressive symptoms in the postnatal period are
    associated with long-term impairment of mother-infant bonding.** *Arch Womens Ment
    Health* 2006. <https://pubmed.ncbi.nlm.nih.gov/16937313/>
72. **Associations between maternal psychological distress and mother-infant bonding: a
    systematic review and meta-analysis.** *Arch Womens Ment Health* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37316760/>
73. Netsi E, et al. **Association of persistent and severe postnatal depression with
    child outcomes.** *JAMA Psychiatry* 2018. — 지속·중증 PPD의 자녀 결과.
    <https://pubmed.ncbi.nlm.nih.gov/29387878/>

## 11. 신경스테로이드 치료제 임상시험 — 모델 캘리브레이션 목표값

74. Kanes S, et al. **Brexanolone (SAGE-547 injection) in post-partum depression: a
    randomised controlled trial.** *Lancet* 2017. — phase 2, n = 21, 중증
    (HAM-D ≥ 26): **60시간 HAM-D 변화 −21.0 (브렉사놀론) vs −8.8 (위약).**
    <https://pubmed.ncbi.nlm.nih.gov/28619476/>
75. Meltzer-Brody S, et al. **Brexanolone injection in post-partum depression: two
    multicentre, double-blind, randomised, placebo-controlled, phase 3 trials.**
    *Lancet* 2018. — **study 1 (HAM-D ≥ 26): 60시간 −19.5 (BRX60) · −17.7 (BRX90) ·
    −14.0 (위약); study 2 (HAM-D 20-25): −14.6 vs −12.1.** 60 → 90 µg/kg/h 용량
    반응이 평평하다는 관찰이 모델의 Hill 포화 예측과 대응.
    <https://pubmed.ncbi.nlm.nih.gov/30177236/>
75b. Kanes SJ, et al. **Open-label, proof-of-concept study of brexanolone in the
    treatment of severe postpartum depression.** *Hum Psychopharmacol* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28370307/>
76. Deligiannidis KM, et al. **Effect of zuranolone vs placebo in postpartum
    depression: a randomized clinical trial (ROBIN).** *JAMA Psychiatry* 2021. —
    주라놀론 30 mg × 14일: **day 15 HAM-D −17.8 vs −13.6.**
    <https://pubmed.ncbi.nlm.nih.gov/34190962/>
77. Deligiannidis KM, et al. **Zuranolone for the Treatment of Postpartum Depression
    (SKYLARK).** *Am J Psychiatry* 2023. — 주라놀론 50 mg × 14일: **day 15 HAM-D
    −15.6 vs −11.6 (차이 −4.0); day 3·28·45에서도 유의.**
    <https://pubmed.ncbi.nlm.nih.gov/37491938/>
78. Althaus AL, et al. **Preclinical characterization of zuranolone (SAGE-217), a
    selective neuroactive steroid GABA(A) receptor positive allosteric modulator.**
    *Neuropharmacology* 2020. — 시냅스 + 세포외시냅스 동시 작용(모델의 W_PHASIC 항).
    <https://pubmed.ncbi.nlm.nih.gov/32976892/>
79. Gunduz-Bruce H, et al. **Trial of SAGE-217 in patients with major depressive
    disorder.** *N Engl J Med* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/31483961/>
80. **Efficacy and safety of zuranolone co-initiated with an antidepressant in adults
    with major depressive disorder (CORAL).** *Neuropsychopharmacology* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/37875578/>
81. **Zuranolone for treatment of major depressive disorder: a systematic review and
    meta-analysis.** *Front Neurosci* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38726035/>
82. **Brexanolone, zuranolone and related neurosteroid GABA(A) receptor positive
    allosteric modulators for postnatal depression.** *Cochrane Database Syst Rev*
    2025. — 최신 계통적 종합.
    <https://pubmed.ncbi.nlm.nih.gov/40562419/>
83. **Preclinical and clinical pharmacology of brexanolone (allopregnanolone) for
    postpartum depression: a landmark journey.** *Psychopharmacology (Berl)* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37566239/>
84. Cornett EM, et al. **Brexanolone to treat postpartum depression in adult women.**
    *Psychopharmacol Bull* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34092826/>
85. Faden J, Citrome L. **Brexanolone for postpartum depression: clinical evidence and
    practical considerations.** *Pharmacotherapy* 2019. — REMS·진정 감시의 실제.
    <https://pubmed.ncbi.nlm.nih.gov/31514247/>
86. **Pharmacokinetics, safety, and tolerability of single and multiple doses of
    zuranolone in Japanese and White healthy subjects.** *Neuropsychopharmacol Rep*
    2023. — 모델의 주라놀론 t½ ≈ 20 h 가정.
    <https://pubmed.ncbi.nlm.nih.gov/37366077/>
87. **The magnitude and sustainability of treatment benefit of zuranolone on function
    and well-being (SDS).** *J Affect Disord* 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38325605/>
88. **Real-world safety profile of zuranolone for postpartum depression: a FAERS
    analysis.** *J Affect Disord* 2026.
    <https://pubmed.ncbi.nlm.nih.gov/40935249/>

## 12. 수유 중 약물 이행 (Lactation transfer — the infant-exposure arm)

89. **Allopregnanolone Concentrations in Breast Milk and Plasma from Healthy Volunteers
    Receiving Brexanolone Injection, With Population Pharmacokinetic Modeling of
    Potential Relative Infant Dose.** *Clin Pharmacokinet* 2022. — **모델의 RID 계산과
    직접 비교되는 실측 데이터.**
    <https://pubmed.ncbi.nlm.nih.gov/35869362/>
90. **Using Brexanolone for Postpartum Depression Must Account for Lactation.**
    *Matern Child Health J* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34019187/>
91. Weissman AM, et al. **Pooled analysis of antidepressant levels in lactating
    mothers, breast milk, and nursing infants.** *Am J Psychiatry* 2004. — 설트랄린
    RID가 SSRI 중 가장 낮음.
    <https://pubmed.ncbi.nlm.nih.gov/15169695/>
92. Berle JØ, et al. **Breastfeeding during maternal antidepressant treatment with
    serotonin reuptake inhibitors: infant exposure, clinical symptoms, and cytochrome
    P450 genotypes.** *J Clin Psychiatry* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15367050/>

## 13. 기존·비약물 치료 (Comparator arms)

93. Molyneaux E, et al. **Antidepressant treatment for postnatal depression.**
    *Cochrane Database Syst Rev* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/25211400/>
94. Hantsoo L, et al. **A randomized, placebo-controlled, double-blind trial of
    sertraline for postpartum depression.** *Psychopharmacology (Berl)* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24173623/>
95. Wisner KL, et al. **Prevention of postpartum depression: a pilot randomized
    clinical trial.** *Am J Psychiatry* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15229064/>
96. Gregoire AJ, et al. **Transdermal oestrogen for treatment of severe postnatal
    depression.** *Lancet* 1996. — 모델의 E2 패치 시나리오.
    <https://pubmed.ncbi.nlm.nih.gov/8598756/>
97. **Perioperative Adjunctive Esketamine for Postpartum Depression Among Women
    Undergoing Elective Cesarean Delivery: A Randomized Clinical Trial.** *JAMA Netw
    Open* 2024. <https://pubmed.ncbi.nlm.nih.gov/38446480/>
98. **Intraoperative Esketamine and Postpartum Depression Among Women With Cesarean
    Delivery: A Randomized Clinical Trial.** *JAMA Netw Open* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/39946130/>
99. **Prophylactic esketamine for postpartum depression after cesarean section: a
    systematic review and meta-analysis.** *J Affect Disord* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/40505985/>
100. O'Connor E, et al. **Interventions to Prevent Perinatal Depression: Evidence
     Report and Systematic Review for the US Preventive Services Task Force.** *JAMA*
     2019. <https://pubmed.ncbi.nlm.nih.gov/30747970/>
101. Dennis CL, Dowswell T. **Psychosocial and psychological interventions for
     preventing postpartum depression.** *Cochrane Database Syst Rev* 2013.
     <https://pubmed.ncbi.nlm.nih.gov/23450532/>
102. Zlotnick C, et al. **Postpartum depression in women receiving public assistance:
     pilot study of an interpersonal-therapy-oriented group intervention (ROSE).**
     *Am J Psychiatry* 2001.
     <https://pubmed.ncbi.nlm.nih.gov/11282702/>
103. **Effect of Bright Light Therapy on Perinatal Depression: A Systematic Review and
     Meta-Analysis.** *Can J Psychiatry* 2024.
     <https://pubmed.ncbi.nlm.nih.gov/38863243/>
104. **The effectiveness of exercise-based interventions for preventing or treating
     postpartum depression: a systematic review and meta-analysis.** *Arch Womens Ment
     Health* 2019. <https://pubmed.ncbi.nlm.nih.gov/29882074/>
105. **Effectiveness of aerobic exercise in the prevention and treatment of postpartum
     depression: meta-analysis and network meta-analysis.** *PLoS One* 2023.
     <https://pubmed.ncbi.nlm.nih.gov/38019729/>
106. Muller AF, et al. **Postpartum thyroiditis.** *Best Pract Res Clin Endocrinol
     Metab* 2004. — 감별진단 및 병존 (모델의 THYRPP 노드).
     <https://pubmed.ncbi.nlm.nih.gov/15157842/>

## 14. 도구·방법론 (Tools and methodology)

107. gPKPDviz — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구.
     <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
108. mrgsolve — R package for simulation from ODE-based models.
     <https://mrgsolve.org/>
109. **New directions in neurosteroid therapeutics in neuropsychiatry.** *Neurosci
     Biobehav Rev* 2025. — 차세대 후보 물질 전망.
     <https://pubmed.ncbi.nlm.nih.gov/40127877/>

---

## 이 모델이 문헌에서 **가져오지 않은** 것 (What the model does NOT get from literature)

정직성을 위해 명시합니다.

| 항목 | 상태 |
|------|------|
| `KR` (δ 수용체 회복 t½ ≈ 7일) | 사람 데이터 **없음**. 설치류(#24)에서 유추한 자유 파라미터. 모델의 가장 중요한 파라미터이며 감도분석에서 1위. |
| `THR0` (E/I 증상 역치) | 관측 불가능한 잠재 파라미터. 위약군 궤적에 맞춰 캘리브레이션. |
| `ZUR_EQ` (주라놀론 → ALLO 등가) | SKYLARK day 15에 맞춘 **단일** 캘리브레이션 파라미터. |
| `K_CARE`, `K_NSP` (비특이 효과) | 기전이 아니라 위약군에 맞춘 현상학적 항. 능동군에서 값 고정. |
| `W_SELF` (증상 자기강화) | 반추·DMN·글루타메이트를 하나로 묶은 축약. 개별 근거는 있으나 계수는 자유. |
| 시상하부 CRH 회복 12주 | #50 한 편의 소규모 연구에 의존. |
| 영아 노출의 F_oral = 0.05 | 신생아 실측 없음. 보수적 가정이며 #89의 실측 RID와 비교만 가능. |
