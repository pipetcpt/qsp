# 정상압 수두증 (iNPH) QSP 모델 — 참고문헌

> **검증 방법 (verification).** 아래 1–18절의 모든 PMID는 NCBI PubMed
> E-utilities API(`esearch` + `esummary`)로 조회하여 **저자·연도·학술지·제목이
> 실제 레코드와 일치함을 확인**한 것입니다. 기억에 의존해 작성한 인용은
> 없습니다. 19절은 반대로 **모델의 어떤 부분이 문헌으로 뒷받침되지 않는지**를
> 명시합니다 — 민감도 분석에서 답을 실제로 소유하는 것으로 밝혀진 파라미터
> (`kW_out`, `Hcol_cm`, `kAQ_rec`)가 바로 그 목록에 있습니다.
>
> 형식: `저자 연도 학술지 — 제목 [PMID]`

---

## 1. 질환 개념·진단기준·가이드라인

1. Adams RD 1965 *N Engl J Med* — Symptomatic occult hydrocephalus with "normal" cerebrospinal-fluid pressure: a treatable syndrome. [PMID 14303656](https://pubmed.ncbi.nlm.nih.gov/14303656/) — 원 기술. 모델의 출발점인 "압력은 정상인데 질환은 있다"는 역설 그 자체.
2. Relkin N 2005 *Neurosurgery* — Diagnosing idiopathic normal-pressure hydrocephalus. [PMID 16160425](https://pubmed.ncbi.nlm.nih.gov/16160425/)
3. Marmarou A 2005 *Neurosurgery* — The value of supplemental prognostic tests for the preoperative assessment of idiopathic normal-pressure hydrocephalus. [PMID 16160426](https://pubmed.ncbi.nlm.nih.gov/16160426/)
4. Marmarou A 2005 *Acta Neurochir Suppl* — Guidelines for management of idiopathic normal pressure hydrocephalus: progress to date. [PMID 16463856](https://pubmed.ncbi.nlm.nih.gov/16463856/)
5. Klinge P 2005 *Neurosurgery* — Outcome of shunting in idiopathic normal-pressure hydrocephalus and the value of outcome assessment in shunted patients. [PMID 16160428](https://pubmed.ncbi.nlm.nih.gov/16160428/)
6. Ishikawa M 2008 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus. [PMID 18408356](https://pubmed.ncbi.nlm.nih.gov/18408356/)
7. Mori E 2012 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus: second edition. [PMID 23183074](https://pubmed.ncbi.nlm.nih.gov/23183074/)
8. Nakajima M 2021 *Neurol Med Chir (Tokyo)* — Guidelines for management of idiopathic normal pressure hydrocephalus (third edition). [PMID 33455998](https://pubmed.ncbi.nlm.nih.gov/33455998/) — 진단 기준·검사 알고리즘의 1차 출처.
9. Tsakanikas D 2007 *Semin Neurol* — Normal pressure hydrocephalus. [PMID 17226742](https://pubmed.ncbi.nlm.nih.gov/17226742/)
10. Rovira À 2026 *Can Assoc Radiol J* — Idiopathic normal pressure hydrocephalus: a comprehensive review. [PMID 42239999](https://pubmed.ncbi.nlm.nih.gov/42239999/)
11. Schmidt E 2025 *J Neurosurg Sci* — Treatment of iNPH: novel insights. [PMID 40045806](https://pubmed.ncbi.nlm.nih.gov/40045806/)
12. Taddei G 2026 *Eur J Health Econ* — Economic analysis and healthcare implications of underdiagnosed idiopathic normal pressure hydrocephalus. [PMID 42201617](https://pubmed.ncbi.nlm.nih.gov/42201617/)

## 2. CSF 역학의 정량 이론 — Davson 식, Marmarou 압력-용적 곡선, 유출저항

13. Marmarou A 1975 *J Neurosurg* — Compartmental analysis of compliance and outflow resistance of the cerebrospinal fluid system. [PMID 1181384](https://pubmed.ncbi.nlm.nih.gov/1181384/) — **모델의 수력학 블록 전체가 이 논문의 형식을 따른다**: 지수 압력-용적 곡선 `C = 1/(E1·(P−P0))`와 유출저항 `Rout`.
14. Yamada S 2023 *Neurol Med Chir (Tokyo)* — Cerebrospinal fluid production and absorption and ventricular enlargement mechanisms in hydrocephalus. [PMID 36858632](https://pubmed.ncbi.nlm.nih.gov/36858632/) — `If0 = 0.35 mL/min`, 흡수 경로 분포.
15. Hansson W 2026 *Neurosurgery* — Intracranial pressure dynamics and cerebrospinal fluid outflow resistance in the elderly. [PMID 41854330](https://pubmed.ncbi.nlm.nih.gov/41854330/) — 고령 정상군의 `Rout` 분포. `Rout_norm = 9`의 근거.
16. Qvarlander S 2014 *Med Biol Eng Comput* — CSF dynamic analysis of a predictive pulsatility-based infusion test for normal pressure hydrocephalus. [PMID 24151060](https://pubmed.ncbi.nlm.nih.gov/24151060/) — 주입검사로 `Rout`와 박동성을 동시에 얻는 방법.
17. Lalou AD 2020 *Fluids Barriers CNS* — Cerebrospinal fluid dynamics in non-acute post-traumatic ventriculomegaly. [PMID 32228689](https://pubmed.ncbi.nlm.nih.gov/32228689/) — `Rout`·RAP·보상예비능의 임상 해석.
18. Zeiler FA 2019 *Acta Neurochir (Wien)* — Compensatory-reserve-weighted intracranial pressure versus intracranial pressure for outcome association. [PMID 31053909](https://pubmed.ncbi.nlm.nih.gov/31053909/) — 평균압보다 보상예비능이 더 정보량이 많다는 논지.
19. Lindstrøm EK 2018 *Neuroimage Clin* — Cerebrospinal fluid volumetric net flow rate and direction in idiopathic normal pressure hydrocephalus. [PMID 30238917](https://pubmed.ncbi.nlm.nih.gov/30238917/)
20. Rosenberg GA 1982 *Brain Res* — The effect of increased CSF pressure on interstitial fluid flow during ventriculocisternal perfusion. [PMID 7055690](https://pubmed.ncbi.nlm.nih.gov/7055690/) — 압력 상승 → 간질액 유동, `Wpv` 구획의 실험적 근거.
21. James AE Jr 1975 *J Neurol Sci* — Pathophysiology of chronic communicating hydrocephalus in dogs. [PMID 1172782](https://pubmed.ncbi.nlm.nih.gov/1172782/) — 경상피(transependymal) 유출과 뇌실 확장의 고전적 실험.

## 3. 박동성 가설 — 평균압이 아니라 파형이 병인이라는 주장

22. Eide PK 2006 *Acta Neurochir (Wien)* — Intracranial pulse pressure amplitude levels determined during preoperative assessment of subjects with possible idiopathic normal pressure hydrocephalus. [PMID 17039303](https://pubmed.ncbi.nlm.nih.gov/17039303/) — **`E1_init`와 `Vp_init`를 이 문헌의 AMP 분포(정상 ~2, iNPH >4 mmHg)에 맞춰 보정했다.**
23. Eide PK 2010 *Cerebrospinal Fluid Res* — Cerebrospinal fluid pulse pressure amplitude during lumbar infusion in idiopathic normal pressure hydrocephalus. [PMID 20205911](https://pubmed.ncbi.nlm.nih.gov/20205911/)
24. Eide PK 2008 *Med Eng Phys* — Comparison of simultaneous continuous intracranial pressure signals from sensors placed within brain parenchyma and epidural space. [PMID 17336574](https://pubmed.ncbi.nlm.nih.gov/17336574/) — AMP 측정 자체의 부위 의존성. 이 모델이 AMP를 "누운 자세"에서만 보고하는 이유.
25. Eide PK 2018 *J Neurosurg* — The pathophysiology of chronic noncommunicating hydrocephalus: lessons from continuous intracranial pressure monitoring. [PMID 28799879](https://pubmed.ncbi.nlm.nih.gov/28799879/)
26. Park EH 2012 *J Neurosurg* — Impaired pulsation absorber mechanism in idiopathic normal pressure hydrocephalus. [PMID 23061391](https://pubmed.ncbi.nlm.nih.gov/23061391/) — 두개내 박동 흡수 실패. `kappa_tm`이 표현하려는 기전.
27. Bateman GA 2002 *Neuroradiology* — Pulse-wave encephalopathy: a comparative study of the hydrodynamics of leukoaraiosis and normal-pressure hydrocephalus. [PMID 12221445](https://pubmed.ncbi.nlm.nih.gov/12221445/)
28. Bateman GA 2003 *Neuroradiology* — The reversibility of reduced cortical vein compliance in normal-pressure hydrocephalus following shunt insertion. [PMID 12592485](https://pubmed.ncbi.nlm.nih.gov/12592485/) — **션트 후 순응도가 수개월에 걸쳐 회복된다는 관찰. `kE1_rec`의 근거.**
29. Bateman GA 2004 *Med Hypotheses* — Pulse wave encephalopathy: a spectrum hypothesis incorporating Alzheimer's disease, vascular dementia and normal pressure hydrocephalus. [PMID 14962623](https://pubmed.ncbi.nlm.nih.gov/14962623/)
30. Bateman GA 2016 *Fluids Barriers CNS* — A comparison between the pathophysiology of multiple sclerosis and normal pressure hydrocephalus. [PMID 27658732](https://pubmed.ncbi.nlm.nih.gov/27658732/)
31. Murambi RT 2025 *Fluids Barriers CNS* — Deformation of brain in normal pressure hydrocephalus is more readily associated with slow vasomotion than with the cardiac pulse. [PMID 40533775](https://pubmed.ncbi.nlm.nih.gov/40533775/) — **모델과 부분적으로 상충하는 문헌**: 이 모델은 심박동 성분에 `kappa_tm`을 걸었다. 19절 참조.

## 4. 글림파틱 배출·AQP4

32. Hasan-Olive MM 2019 *Glia* — Loss of perivascular aquaporin-4 in idiopathic normal pressure hydrocephalus. [PMID 30306658](https://pubmed.ncbi.nlm.nih.gov/30306658/) — **`AQ_init = 0.45` (정상 0.85)의 직접 근거.**
33. Eide PK 2020 *Brain Commun* — MRI biomarkers of cerebrospinal fluid tracer dynamics in idiopathic normal pressure hydrocephalus. [PMID 33381757](https://pubmed.ncbi.nlm.nih.gov/33381757/)
34. Jacobsen HH 2020 *Invest Ophthalmol Vis Sci* — In vivo evidence for impaired glymphatic function in the visual pathway of patients with normal pressure hydrocephalus. [PMID 33201186](https://pubmed.ncbi.nlm.nih.gov/33201186/)
35. Eide PK 2019 *Gerontol Geriatr Med* — In vivo imaging of molecular clearance from human entorhinal cortex. [PMID 31819895](https://pubmed.ncbi.nlm.nih.gov/31819895/)
36. Eide PK 2022 *J Cereb Blood Flow Metab* — Altered glymphatic enhancement of cerebrospinal fluid tracer in individuals with chronic poor sleep quality. [PMID 35350917](https://pubmed.ncbi.nlm.nih.gov/35350917/) — **수면 → 글림파틱 결합(`sleep_base`, 멜라토닌 시나리오 S19)의 인체 근거.**
37. Yang Y 2025 *Eur J Neurol* — Alterations of glymphatic system before and after shunt surgery in patients with idiopathic normal pressure hydrocephalus. [PMID 40365713](https://pubmed.ncbi.nlm.nih.gov/40365713/) — 션트 후 글림파틱 회복. 9절 예측의 검증 대상.
38. Mossige I 2026 *Radiology* — Comparing glymphatic function measures: DTI-ALPS and CSF tracer dynamics. [PMID 41631990](https://pubmed.ncbi.nlm.nih.gov/41631990/)
39. Broggi M 2025 *Neurol Sci* — Implications of the glymphatic system in the diagnostic and surgical workup of normal pressure hydrocephalus. [PMID 40524080](https://pubmed.ncbi.nlm.nih.gov/40524080/)
40. He W 2025 *Neurosurg Rev* — Association between perivascular spaces, DTI-derived indices, and choroid plexus with ventriculomegaly. [PMID 41081971](https://pubmed.ncbi.nlm.nih.gov/41081971/)
41. Zahran A 2026 *CNS Neurosci Ther* — Glymphatic system dysfunction in central nervous system diseases. [PMID 41792880](https://pubmed.ncbi.nlm.nih.gov/41792880/)
42. Herz J 2018 *Methods Mol Biol* — Morphological and functional analysis of CNS-associated lymphatics. [PMID 30242757](https://pubmed.ncbi.nlm.nih.gov/30242757/) — 수막 림프관 유출 경로.

## 5. 영상 형태학 — DESH, Evans index, 뇌각(callosal angle)

43. Yamada S 2026 *Fluids Barriers CNS* — CT-based automatic segmentation of key CSF regions for detecting DESH. [PMID 42337617](https://pubmed.ncbi.nlm.nih.gov/42337617/)
44. Skalický P 2021 *J Clin Neurosci* — Role of DESH, callosal angle and cingulate sulcus sign in prediction of gait responsiveness after shunting. [PMID 33334664](https://pubmed.ncbi.nlm.nih.gov/33334664/) — `DESH` 배수의 임상적 의미.
45. Rohatgi S 2024 *J Comput Assist Tomogr* — Correlating Evans index, callosal angle, and lateral ventricle volume with gait response outcomes. [PMID 38657140](https://pubmed.ncbi.nlm.nih.gov/38657140/)
46. Rohatgi S 2024 *Cureus* — Predicting gait speed improvement in idiopathic normal pressure hydrocephalus patients: the role of Evans index. [PMID 38975434](https://pubmed.ncbi.nlm.nih.gov/38975434/)
47. Selcuk M 2026 *Surg Radiol Anat* — Normative Evans index values in the adult population: an age- and sex-stratified MRI study. [PMID 42262516](https://pubmed.ncbi.nlm.nih.gov/42262516/) — `evans_index()` 정상 기준선.
48. Neikter J 2020 *AJNR Am J Neuroradiol* — Ventricular volume is more strongly associated with clinical improvement than the Evans index after shunting. [PMID 32527841](https://pubmed.ncbi.nlm.nih.gov/32527841/) — **모델이 Evans index를 예후 지표로 쓰지 않는 이유.**
49. Virhammar J 2019 *J Neurosurg* — Increase in callosal angle and decrease in ventricular volume after shunt surgery. [PMID 29393749](https://pubmed.ncbi.nlm.nih.gov/29393749/) — 션트 후 형태 변화의 크기. `f_plastic`(가소성 확장 비율)의 보정 표적.
50. Holmgren RT 2026 *Fluids Barriers CNS* — Ventricular volumetry in relation to clinical response and overdrainage after shunt surgery. [PMID 42381066](https://pubmed.ncbi.nlm.nih.gov/42381066/)
51. Yeh PY 2026 *J Neuroradiol* — Predicting ventricular volume reduction and overdrainage in idiopathic normal pressure hydrocephalus. [PMID 41213359](https://pubmed.ncbi.nlm.nih.gov/41213359/)
52. Sahuquillo J 2026 *Biomedicines* — Redefining idiopathic normal pressure hydrocephalus using AI-driven brain volumetry. [PMID 41898322](https://pubmed.ncbi.nlm.nih.gov/41898322/)
53. Barough SS 2025 *medRxiv* — Disproportionately elevated sulcal index (DESI). [PMID 41404281](https://pubmed.ncbi.nlm.nih.gov/41404281/) — *프리프린트(peer-review 전)*.

## 6. 뇌혈류·뇌실주위 백질

54. Virhammar J 2017 *AJNR Am J Neuroradiol* — Arterial spin-labeling perfusion MR imaging demonstrates regional CBF decrease in idiopathic normal pressure hydrocephalus. [PMID 28860216](https://pubmed.ncbi.nlm.nih.gov/28860216/) — **`CBF0`와 iNPH의 뇌실주위 관류 저하폭 보정 근거.**
55. Ziegelitz D 2014 *J Magn Reson Imaging* — Cerebral perfusion measured by DSC MRI is reduced in patients with idiopathic normal pressure hydrocephalus. [PMID 24006249](https://pubmed.ncbi.nlm.nih.gov/24006249/)
56. Ziegelitz D 2016 *J Cereb Blood Flow Metab* — Pre- and postoperative cerebral blood flow changes in patients with idiopathic normal pressure hydrocephalus. [PMID 26661191](https://pubmed.ncbi.nlm.nih.gov/26661191/) — 션트 후 관류 회복.
57. Tanaka A 1997 *Neurosurgery* — Cerebral blood flow and autoregulation in normal pressure hydrocephalus. [PMID 9179888](https://pubmed.ncbi.nlm.nih.gov/9179888/) — `Autoreg` 상태변수의 근거.
58. Tuniz F 2017 *Fluids Barriers CNS* — The role of perfusion and diffusion MRI in the assessment of patients affected by probable idiopathic normal pressure hydrocephalus. [PMID 28899431](https://pubmed.ncbi.nlm.nih.gov/28899431/)
59. Ades-Aron B 2018 *AJNR Am J Neuroradiol* — Diffusional kurtosis along the corticospinal tract in adult normal pressure hydrocephalus. [PMID 30385473](https://pubmed.ncbi.nlm.nih.gov/30385473/) — **하지 피질척수로가 가장 뇌실 인접 경로라는 해부학적 근거(`CST_leg`).**
60. Kanno S 2017 *Fluids Barriers CNS* — A change in brain white matter after shunt surgery in idiopathic normal pressure hydrocephalus: a tract-based spatial statistics study. [PMID 28132644](https://pubmed.ncbi.nlm.nih.gov/28132644/) — 회복 가능한 성분의 존재.
61. Di Curzio DL 2013 *Exp Neurol* — Reduced subventricular zone proliferation and white matter damage in juvenile ferrets with kaolin-induced hydrocephalus. [PMID 23769908](https://pubmed.ncbi.nlm.nih.gov/23769908/)
62. Zadka Y 2023 *J Appl Physiol* — Mechanisms of reduced cerebral blood flow in cerebral edema and elevated intracranial pressure. [PMID 36603049](https://pubmed.ncbi.nlm.nih.gov/36603049/) — `kCBF_W`(간질 부종 → 관류 저하)의 기전적 근거.

## 7. CSF 바이오마커 — 그리고 희석 인공물

63. Jingami N 2015 *J Alzheimers Dis* — Idiopathic normal pressure hydrocephalus has a different cerebrospinal fluid biomarker profile from Alzheimer's disease. [PMID 25428256](https://pubmed.ncbi.nlm.nih.gov/25428256/) — **9절의 핵심 관찰: iNPH에서 Aβ42와 p-tau가 *모두* 낮다.**
64. Jingami N 2019 *J Alzheimers Dis* — Two-point dynamic observation of Alzheimer's disease cerebrospinal fluid biomarkers in idiopathic normal pressure hydrocephalus. [PMID 31561378](https://pubmed.ncbi.nlm.nih.gov/31561378/)
65. Miyajima M 2013 *PLoS One* — Leucine-rich α2-glycoprotein is a novel biomarker of neurodegenerative disease in human cerebrospinal fluid. [PMID 24058569](https://pubmed.ncbi.nlm.nih.gov/24058569/) — `A_lrg` 구획.
66. Grønning R 2023 *Fluids Barriers CNS* — Association between ventricular CSF biomarkers and outcome after shunt surgery in idiopathic normal pressure hydrocephalus. [PMID 37880775](https://pubmed.ncbi.nlm.nih.gov/37880775/) — **측정 부위(뇌실 vs 요추)에 따라 값이 달라진다 — 모델이 단일 CSF 구획을 쓰는 한계.**
67. Migliorati K 2021 *Neurol Res* — P-tau as prognostic marker in long-term follow-up for patients with shunted iNPH. [PMID 33059546](https://pubmed.ncbi.nlm.nih.gov/33059546/)
68. Kanemoto H 2023 *Int Psychogeriatr* — Cerebrospinal fluid amyloid beta and response of cognition to a tap test in idiopathic normal pressure hydrocephalus. [PMID 34399871](https://pubmed.ncbi.nlm.nih.gov/34399871/)
69. Cihlo M 2025 *Neurosurg Rev* — Value of biomarkers in the prediction of shunt responsivity in patients with normal pressure hydrocephalus. [PMID 40457021](https://pubmed.ncbi.nlm.nih.gov/40457021/)
70. Kwiecień A 2026 *Molecules* — Decoding the CSF proteomic signature of idiopathic normal pressure hydrocephalus: a systematic review. [PMID 42451686](https://pubmed.ncbi.nlm.nih.gov/42451686/)
71. Lolansen SD 2021 *Dis Markers* — Inflammatory markers in cerebrospinal fluid from patients with hydrocephalus: a systematic literature review. [PMID 33613789](https://pubmed.ncbi.nlm.nih.gov/33613789/)

## 8. 동반 알츠하이머 병리 — 반응자/비반응자의 실체

72. Ye BS 2025 *Alzheimers Dement* — Degenerative pathologies on cortical biopsy, dopaminergic depletion, and shunt efficacy in iNPH. [PMID 41366839](https://pubmed.ncbi.nlm.nih.gov/41366839/) — **`APOE`/`Ab_plq` 항이 모델의 유일한 비반응자 기전인 근거.**
73. Luikku AJ 2024 *J Neuropathol Exp Neurol* — Deep learning assisted quantitative analysis of Aβ and microglia in patients with idiopathic normal pressure hydrocephalus. [PMID 39101555](https://pubmed.ncbi.nlm.nih.gov/39101555/)
74. Greenberg ABW 2024 *Cereb Cortex* — Utility of cortical tissue analysis in normal pressure hydrocephalus. [PMID 38275188](https://pubmed.ncbi.nlm.nih.gov/38275188/)
75. Mattoli MV 2020 *Int J Mol Sci* — Usefulness of brain PET with different tracers in the evaluation of patients with idiopathic normal pressure hydrocephalus. [PMID 32906629](https://pubmed.ncbi.nlm.nih.gov/32906629/)

## 9. 진단적 섭동 — tap test, 지속 요추 배액, 주입검사

76. Rydja J 2021 *Fluids Barriers CNS* — Evaluating the cerebrospinal fluid tap test with the Hellström iNPH scale. [PMID 33827613](https://pubmed.ncbi.nlm.nih.gov/33827613/)
77. Liu C 2021 *Clin Neurol Neurosurg* — A pilot study of multiple time points and multidomain assessment in cerebrospinal fluid tap test. [PMID 34749022](https://pubmed.ncbi.nlm.nih.gov/34749022/) — **평가 시점이 결과를 바꾼다 — 6절의 시간상수 불일치 논지와 직접 대응.**
78. Mládek A 2022 *Neurosurgery* — Prediction of shunt responsiveness in suspected patients with normal pressure hydrocephalus using the lumbar infusion test. [PMID 35080523](https://pubmed.ncbi.nlm.nih.gov/35080523/)
79. Akar K 2026 *Fluids Barriers CNS* — Sensitivity of physiotherapy-based clinical tests in detecting change in gait and balance performance. [PMID 41703573](https://pubmed.ncbi.nlm.nih.gov/41703573/) — **검사 민감도가 관찰자의 판정 역치에 달려 있다는 근거(6절).**
80. Wolfsegger T 2017 *J Neurol Sci* — Cognitive impairment predicts worse short-term response to spinal tap test in normal pressure hydrocephalus. [PMID 28716246](https://pubmed.ncbi.nlm.nih.gov/28716246/)
81. Kudelić N 2023 *Front Neurol* — Predictive value of spinal CSF volume in the preoperative assessment of patients with idiopathic normal-pressure hydrocephalus. [PMID 37869132](https://pubmed.ncbi.nlm.nih.gov/37869132/)
82. Cai H 2025 *J Med Internet Res* — Predictive value of digital neuropsychological and gait assessments on shunt outcome. [PMID 41289581](https://pubmed.ncbi.nlm.nih.gov/41289581/)
83. Guarracino I 2025 *Brain Sci* — Real-time neuropsychological testing during infusion in hydrocephalus. [PMID 39851404](https://pubmed.ncbi.nlm.nih.gov/39851404/)

## 10. 션트 무작위 시험 — 밸브 압력과 수술 자체

84. Boon AJ 1998 *J Neurosurg* — Dutch Normal-Pressure Hydrocephalus Study: randomized comparison of low- and medium-pressure shunts. [PMID 9488303](https://pubmed.ncbi.nlm.nih.gov/9488303/) — **5절 적정 지도가 재현하려는 시험. 낮은 개방압이 더 효과적이지만 과배액 합병증이 늘어난다는 양면성.**
85. Lemcke J 2013 *J Neurol Neurosurg Psychiatry* — Safety and efficacy of gravitational shunt valves in patients with idiopathic normal pressure hydrocephalus (SVASONA). [PMID 23457222](https://pubmed.ncbi.nlm.nih.gov/23457222/) — **중력식 보조기 무작위 시험. 3절·5절의 핵심 근거.**
86. Lemcke J 2012 *Acta Neurochir Suppl* — On the method of a randomised comparison of programmable valves with and without gravitational units. [PMID 22327702](https://pubmed.ncbi.nlm.nih.gov/22327702/)
87. Lemcke J 2010 *Acta Neurochir Suppl* — Is it possible to minimize overdrainage complications with gravitational units in patients with idiopathic normal pressure hydrocephalus? [PMID 19812931](https://pubmed.ncbi.nlm.nih.gov/19812931/)
88. Luciano M 2023 *Neurosurgery* — Placebo-controlled effectiveness of idiopathic normal pressure hydrocephalus shunting: a randomized pilot. [PMID 36700738](https://pubmed.ncbi.nlm.nih.gov/36700738/)
89. Luciano MG 2025 *N Engl J Med* — A randomized trial of shunting for idiopathic normal-pressure hydrocephalus. [PMID 40960253](https://pubmed.ncbi.nlm.nih.gov/40960253/) — **위약(밸브 폐쇄) 대조 시험. 모델이 션트 효과를 과대평가한다는 12절 자기비판의 기준점.**
90. Saper CB 2026 *Ann Neurol* — Time to reconsider the value of shunting procedures for "idiopathic normal pressure hydrocephalus". [PMID 41741941](https://pubmed.ncbi.nlm.nih.gov/41741941/) — **모델의 결론과 정면으로 긴장 관계에 있는 논평. 반드시 함께 읽어야 한다.**
91. Klinge P 2012 *Acta Neurol Scand* — One-year outcome in the European multicentre study on iNPH. [PMID 22571428](https://pubmed.ncbi.nlm.nih.gov/22571428/) — 실제 반응률 60–80%.
92. Salih A 2024 *EClinicalMedicine* — The effectiveness of various CSF diversion surgeries in idiopathic normal pressure hydrocephalus: a systematic review. [PMID 39539993](https://pubmed.ncbi.nlm.nih.gov/39539993/)
93. Brenner LBO 2026 *Neurosurgery* — The Lumboperitoneal Shunt Study: a systematic review and single-arm meta-analysis of 2696 patients. [PMID 42294933](https://pubmed.ncbi.nlm.nih.gov/42294933/)
94. Burrows EJ 2026 *Neurosurg Pract* — Clinical subgroups and treatment outcomes in idiopathic normal pressure hydrocephalus. [PMID 42282951](https://pubmed.ncbi.nlm.nih.gov/42282951/)

## 11. 션트 수력학·정수압·과배액 — 이 모델의 3절

95. de Jong DA 2000 *Acta Neurochir (Wien)* — Hydrostatic and hydrodynamic considerations in shunted normal pressure hydrocephalus. [PMID 10819253](https://pubmed.ncbi.nlm.nih.gov/10819253/) — **3절의 정수압 컬럼 계산이 따르는 문헌. `Hcol_cm`의 근거.**
96. Kajimoto Y 2000 *J Neurosurg* — Posture-related changes in the pressure environment of the ventriculoperitoneal shunt system. [PMID 11014539](https://pubmed.ncbi.nlm.nih.gov/11014539/) — **직립 시 구동압이 수십 cmH₂O 증가한다는 직접 측정. `f_up` 가중 방식의 정당화.**
97. Medow JE 2012 *J Neurosurg Pediatr* — Posture-independent piston valve: a novel valve mechanism that actuates based on intracranial pressure alone. [PMID 22208323](https://pubmed.ncbi.nlm.nih.gov/22208323/)
98. Corin AS 2026 *Neurochirurgie* — Effectiveness of anti-siphon devices in CSF shunts for preventing overdrainage in normal pressure hydrocephalus. [PMID 42486335](https://pubmed.ncbi.nlm.nih.gov/42486335/) — `asd_eff`.
99. Hafizka Y 2026 *World Neurosurg* — Effect of anti-siphon devices on postoperative outcomes in idiopathic normal pressure hydrocephalus. [PMID 41941959](https://pubmed.ncbi.nlm.nih.gov/41941959/)
100. Meier U 2013 *Neurosurgery* — Predictors of subsequent overdrainage and clinical outcomes after ventriculoperitoneal shunting for idiopathic normal pressure hydrocephalus. [PMID 24257332](https://pubmed.ncbi.nlm.nih.gov/24257332/) — 위축(`atrophy`)이 과배액 위험인자.
101. Chung K 2026 *JAAPA* — Ventriculoperitoneal shunt management in patients with normal pressure hydrocephalus and cerebral atrophy. [PMID 41662144](https://pubmed.ncbi.nlm.nih.gov/41662144/)
102. Kawahara T 2022 *Surg Neurol Int* — Dural sac shrinkage signs on spinal MRI indicate overdrainage after lumboperitoneal shunt. [PMID 35855156](https://pubmed.ncbi.nlm.nih.gov/35855156/)
103. Younes B 2026 *Childs Nerv Syst* — Incidence and management of overdrainage in hydrocephalus patients treated with adjustable valves. [PMID 41758240](https://pubmed.ncbi.nlm.nih.gov/41758240/)
104. Chelmis F 2026 *Cureus* — Diagnostic challenges and surgical outcomes of Miyazaki syndrome. [PMID 42460193](https://pubmed.ncbi.nlm.nih.gov/42460193/) — 과배액의 극단 표현형.
105. Iimori T 2026 *Clin Neurol Neurosurg* — Safety of continuing antithrombotic therapy during lumboperitoneal shunting for idiopathic normal pressure hydrocephalus. [PMID 42068899](https://pubmed.ncbi.nlm.nih.gov/42068899/) — `antithrombotic` 위험 계수.

## 12. 밸브 선택·적정·재수술

106. Reis RC 2026 *J Neurol Surg A Cent Eur Neurosurg* — Optimizing initial shunt pressure in idiopathic normal pressure hydrocephalus. [PMID 41571242](https://pubmed.ncbi.nlm.nih.gov/41571242/) — **5절 적정 지도의 임상적 대응물.**
107. Colonna S 2026 *World Neurosurg* — Valve selection and long-term outcomes in idiopathic normal pressure hydrocephalus. [PMID 41720264](https://pubmed.ncbi.nlm.nih.gov/41720264/)
108. Ahmed M 2023 *Clin Neurol Neurosurg* — Fixed versus adjustable differential pressure valves in idiopathic normal pressure hydrocephalus. [PMID 37209623](https://pubmed.ncbi.nlm.nih.gov/37209623/)
109. Katiyar V 2021 *Neurol India* — Comparison of programmable and non-programmable shunts for normal pressure hydrocephalus: a meta-analysis. [PMID 35102997](https://pubmed.ncbi.nlm.nih.gov/35102997/)
110. Oksa S 2026 *Acta Neurochir (Wien)* — Reduced risk of shunt revision with adjustable valves: a population-based cohort study over three decades. [PMID 41634436](https://pubmed.ncbi.nlm.nih.gov/41634436/) — 재수술률(`k_occl` 시나리오 S18).

## 13. 제3뇌실 천공술 (ETV)

111. Scalia G 2025 *Brain Sci* — The effects of endoscopic third ventriculostomy versus ventriculoperitoneal shunt on neuropsychological outcomes. [PMID 40426679](https://pubmed.ncbi.nlm.nih.gov/40426679/)
112. Harbaugh TD 2025 *Cureus* — Ventriculoperitoneal shunting versus endoscopic third ventriculostomy for the surgical management of idiopathic normal pressure hydrocephalus. [PMID 40034632](https://pubmed.ncbi.nlm.nih.gov/40034632/) — **S20(ETV 무효 예측)의 임상적 대조.**
113. Mohamed B 2024 *Cureus* — A single-centre experience of the management and surgical outcomes of late-onset idiopathic aqueductal stenosis. [PMID 38868257](https://pubmed.ncbi.nlm.nih.gov/38868257/) — 폐쇄성 수두증에서는 ETV가 작동한다는 대비.

## 14. 약물치료 — acetazolamide의 PK/PD

114. Kunka RL 1979 *J Pharm Sci* — Nonlinear model for acetazolamide. [PMID 423125](https://pubmed.ncbi.nlm.nih.gov/423125/) — **포화 적혈구 결합(`az_Bmax`, `az_kon`, `az_koff`) 구조의 직접 근거.**
115. Maren TH 1979 *J Appl Physiol* — Effect of varying CO₂ equilibria on rates of HCO₃⁻ formation in cerebrospinal fluid. [PMID 118142](https://pubmed.ncbi.nlm.nih.gov/118142/) — 탄산탈수효소–CSF 형성 결합.
116. Vogh BP 1975 *Am J Physiol* — Sodium, chloride, and bicarbonate movement from plasma to cerebrospinal fluid. [PMID 803792](https://pubmed.ncbi.nlm.nih.gov/803792/) — `az_Emax = 0.50`의 동물 근거.
117. Alperin N 2014 *Neurology* — Low-dose acetazolamide reverses periventricular white matter hyperintensities in iNPH. [PMID 24634454](https://pubmed.ncbi.nlm.nih.gov/24634454/)
118. Gilbert GJ 2014 *Neurology* — Low-dose acetazolamide reverses periventricular white matter hyperintensities in iNPH (correspondence). [PMID 25367060](https://pubmed.ncbi.nlm.nih.gov/25367060/)
119. Virhammar J 2026 *Lancet Neurol* — Safety, tolerability, and efficacy of acetazolamide in idiopathic normal pressure hydrocephalus (DRAIN). [PMID 42127932](https://pubmed.ncbi.nlm.nih.gov/42127932/) — **8절 결론(만성 효과는 작고 부작용은 용량 의존적)의 임상 검증 대상. 이 무작위 시험 결과가 모델의 S12/S13 해석을 판정한다.**

## 15. CSF 분비 생리 — 분자 수준

120. Jensen DB 2025 *Adv Sci* — The Na⁺,K⁺,2Cl⁻ cotransporter, not aquaporin 1, sustains cerebrospinal fluid secretion. [PMID 39692709](https://pubmed.ncbi.nlm.nih.gov/39692709/) — **`NKCC1_cp` 노드와 루프계 이뇨제(`bu_Emax`) 표적의 근거. 지도에서 AQP1보다 NKCC1을 상위에 둔 이유.**
121. Deffner F 2022 *Cell Mol Life Sci* — Aquaporin-4 expression in the human choroid plexus. [PMID 35072772](https://pubmed.ncbi.nlm.nih.gov/35072772/)

## 16. 임상 증상·평가 척도

122. Hellström P 2012 *Acta Neurol Scand* — A new scale for assessment of severity and outcome in iNPH. [PMID 22587624](https://pubmed.ncbi.nlm.nih.gov/22587624/)
123. Passaretti M 2023 *Mov Disord Clin Pract* — Gait analysis in idiopathic normal pressure hydrocephalus: a meta-analysis. [PMID 38026510](https://pubmed.ncbi.nlm.nih.gov/38026510/) — **`G_max`와 iNPH 보행속도 기준선(≈0.6 m/s)의 보정 근거.**
124. Mills R 2026 *Fluids Barriers CNS* — Three-dimensional kinematic gait signatures of idiopathic normal pressure hydrocephalus. [PMID 42174608](https://pubmed.ncbi.nlm.nih.gov/42174608/)
125. Ekblom M 2025 *Fluids Barriers CNS* — Blinded gait assessment in idiopathic normal pressure hydrocephalus: reliability and correlation with clinical scales. [PMID 40926257](https://pubmed.ncbi.nlm.nih.gov/40926257/)
126. Lander J 2026 *Fluids Barriers CNS* — Slow unsteady gait in a population-based cohort: links to ventriculomegaly and iNPH-related imaging markers. [PMID 42286688](https://pubmed.ncbi.nlm.nih.gov/42286688/)
127. Sakakibara R 2008 *Neurourol Urodyn* — Mechanism of bladder dysfunction in idiopathic normal pressure hydrocephalus. [PMID 18092331](https://pubmed.ncbi.nlm.nih.gov/18092331/) — **`Urin` 구획의 전두–교뇌 배뇨 회로 근거.**
128. Krzastek SC 2017 *Neurourol Urodyn* — Characterization of lower urinary tract symptoms in patients with idiopathic normal pressure hydrocephalus. [PMID 27490149](https://pubmed.ncbi.nlm.nih.gov/27490149/)
129. Rovčanin B 2025 *Diagnostics* — Time course of symptoms in normal-pressure hydrocephalus: a systematic review. [PMID 40722526](https://pubmed.ncbi.nlm.nih.gov/40722526/) — 7절(조기 vs 지연 수술)의 임상적 배경.
130. Kimura T 2021 *World Neurosurg* — Preoperative predictive factors of short-term outcome in idiopathic normal pressure hydrocephalus. [PMID 33895373](https://pubmed.ncbi.nlm.nih.gov/33895373/)
131. Uchigami H 2022 *Clin Neurol Neurosurg* — Preoperative factors associated with shunt responsiveness in patients with idiopathic normal-pressure hydrocephalus. [PMID 36049404](https://pubmed.ncbi.nlm.nih.gov/36049404/) — **증상 지속기간이 예후를 예측한다는 관찰. `k_perm`(비가역 전환률)이 표현하려는 현상.**
132. Mansour MA 2025 *J Neurosurg Case Lessons* — Beyond the triad: akinetic mutism in idiopathic normal pressure hydrocephalus. [PMID 40720908](https://pubmed.ncbi.nlm.nih.gov/40720908/)

## 17. 병인·유전·상피

133. Sato H 2016 *PLoS One* — A segmental copy number loss of the SFMBT1 gene is a genetic risk for shunt-responsive idiopathic normal pressure hydrocephalus. [PMID 27861535](https://pubmed.ncbi.nlm.nih.gov/27861535/)
134. Tipton PW 2023 *Neurol Genet* — CWH43 variants are associated with disease risk and clinical phenotypic measures in patients with normal pressure hydrocephalus. [PMID 37476022](https://pubmed.ncbi.nlm.nih.gov/37476022/)
135. Yang HW 2023 *Proc Natl Acad Sci USA* — A role for mutations in AK9 and other genes affecting ependymal cells in idiopathic normal pressure hydrocephalus. [PMID 38100419](https://pubmed.ncbi.nlm.nih.gov/38100419/) — `Epend_loss` 노드.
136. Piccinin CC 2025 *Mov Disord* — Genetic risk factors in normal pressure hydrocephalus: what we know and what is next. [PMID 40266017](https://pubmed.ncbi.nlm.nih.gov/40266017/)
137. Yang D 2021 *Proc Natl Acad Sci USA* — Increased plasmin-mediated proteolysis of L1CAM in a mouse model of idiopathic normal pressure hydrocephalus. [PMID 34380733](https://pubmed.ncbi.nlm.nih.gov/34380733/)
138. Botfield H 2013 *Brain* — Decorin prevents the development of juvenile communicating hydrocephalus. [PMID 23983032](https://pubmed.ncbi.nlm.nih.gov/23983032/) — **TGF-β1 → 지주막 섬유화 → `Rout` 상승 경로의 인과적 실험 근거.**
139. Yan H 2016 *Brain Res* — Decorin alleviated chronic hydrocephalus via inhibiting TGF-β1/Smad/CTGF pathway after subarachnoid hemorrhage. [PMID 26556770](https://pubmed.ncbi.nlm.nih.gov/26556770/)
140. Tan Q 2017 *Brain Res* — Cannabinoid receptor 2 activation restricts fibrosis and alleviates hydrocephalus after intraventricular hemorrhage. [PMID 27769788](https://pubmed.ncbi.nlm.nih.gov/27769788/)
141. Feng Z 2020 *Neurosci Lett* — uPA alleviates kaolin-induced hydrocephalus. [PMID 32497735](https://pubmed.ncbi.nlm.nih.gov/32497735/)
142. Ishikawa M 2008 *Brain Nerve* — [Idiopathic normal pressure hydrocephalus — overviews and pathogenesis]. [PMID 18402067](https://pubmed.ncbi.nlm.nih.gov/18402067/) — *일본어*.

## 18. 수면호흡장애·감별진단·수학적 모델링

143. Riedel CS 2022 *Sleep* — Sleep-disordered breathing is frequently associated with idiopathic normal pressure hydrocephalus. [PMID 34739077](https://pubmed.ncbi.nlm.nih.gov/34739077/) — `OSA_node` → `Sleep_SWS` → 글림파틱.
144. Román GC 2019 *Curr Neurol Neurosci Rep* — Sleep-disordered breathing and idiopathic normal-pressure hydrocephalus: recent pathophysiological advances. [PMID 31144048](https://pubmed.ncbi.nlm.nih.gov/31144048/)
145. Riedel CS 2024 *Fluids Barriers CNS* — Transient intracranial pressure elevations (B waves) associated with sleep. [PMID 39702226](https://pubmed.ncbi.nlm.nih.gov/39702226/)
146. Yun S 2025 *Neuroradiology* — Data-driven differentiation of idiopathic normal-pressure hydrocephalus and progressive supranuclear palsy. [PMID 41204957](https://pubmed.ncbi.nlm.nih.gov/41204957/)
147. Deng Z 2024 *Clin Neurol Neurosurg* — Evaluation of imaging indicators in differentiating idiopathic normal pressure hydrocephalus from Alzheimer's disease. [PMID 38823198](https://pubmed.ncbi.nlm.nih.gov/38823198/)
148. Netterwala A 2026 *BMJ Case Rep* — Recognising CADASIL in adults with NPH-like syndrome. [PMID 42236108](https://pubmed.ncbi.nlm.nih.gov/42236108/)
149. Suppa A 2025 *Mov Disord* — Neurophysiology of atypical parkinsonian syndromes: a study group position paper. [PMID 40356334](https://pubmed.ncbi.nlm.nih.gov/40356334/)
150. Dreyer LW 2024 *Fluids Barriers CNS* — Modeling CSF circulation and the glymphatic system during infusion using subject-specific intracranial pressures. [PMID 39407250](https://pubmed.ncbi.nlm.nih.gov/39407250/) — **가장 가까운 선행 모델링 연구. 이 모델과의 차이는 약물 PK/PD와 션트 하드웨어를 함께 갖는다는 점.**
151. Chu KH 2024 *Brain Spine* — Mathematical modelling of cerebral haemodynamics and their effects on ICP. [PMID 38510619](https://pubmed.ncbi.nlm.nih.gov/38510619/)
152. Doron O 2021 *Fluids Barriers CNS* — Interactions of brain, blood, and CSF: a novel mathematical model of cerebral edema. [PMID 34530863](https://pubmed.ncbi.nlm.nih.gov/34530863/)
153. Ursino M 2010 *Ann Biomed Eng* — A model of cerebrovascular reactivity including the circle of Willis and cortical anastomoses. [PMID 20094916](https://pubmed.ncbi.nlm.nih.gov/20094916/) — 뇌혈관 순응도·자동조절 모델링의 고전적 계보.
154. Lee KJ 2015 *Biomed Eng Online* — Non-invasive detection of intracranial hypertension using a simplified intracranial hemo- and hydro-dynamic model. [PMID 26024843](https://pubmed.ncbi.nlm.nih.gov/26024843/)
155. Gadda G 2015 *Am J Physiol Heart Circ Physiol* — A new hemodynamic model for the study of cerebral venous outflow. [PMID 25398980](https://pubmed.ncbi.nlm.nih.gov/25398980/)
156. Lampe R 2014 *Comput Math Methods Med* — Mathematical modelling of cerebral blood circulation and cerebral autoregulation. [PMID 25126111](https://pubmed.ncbi.nlm.nih.gov/25126111/)

---

## 19. 문헌으로 뒷받침되지 **않는** 부분 (what is NOT sourced)

QSP 모델에서 위험한 것은 인용이 없는 파라미터가 아니라 **인용이 없다는
사실이 감춰진 파라미터**입니다. 아래는 위 문헌 어디에서도 직접 측정값을
얻을 수 없었던 항목이며, `inph_model_report.txt` 12절과 동일한 목록입니다.

| 파라미터 | 값 | 상태 |
|---|---|---|
| `kW_out` | 1.300 /day | 뇌실주위 간질 수분 배출률. **민감도 분석에서 1위(탄력도 +2.27) — 답이 실제로 매달려 있는 파라미터인데 자릿수 수준의 추정이다.** 20·21번 문헌이 경상피 유동의 존재를 지지하지만 속도는 지지하지 않는다. |
| `kAQ_rec` | 0.008 /day | 감압 후 AQP4 재극성화 속도. **민감도 3위.** 32번 문헌이 극성 소실을 보이지만 회복 속도의 인체 측정값은 없다. |
| `kappa_tm` | 0.060 | **인체 뇌실에서 측정된 적이 없다.** ICP 박동 중 mantle을 횡단하는 분율. **모델 전체에서 가장 영향력이 클 것으로 예상했으나 민감도 분석에서 14개 중 10위(−0.049)로 나왔다 — 예상과 반박을 둘 다 보고서에 남겼다.** 26·27번 문헌이 기전의 존재를 지지하지만 크기는 지지하지 않는다. 31번 문헌은 심박동보다 느린 vasomotion이 변형을 더 잘 설명한다고 보고하여 이 항의 **형태 자체와 긴장 관계**에 있다. |
| `DESH` | 1.55 | `kappa_tm`에 곱하는 형태학 배수. 44번 문헌이 DESH가 예후와 연관됨을 보이지만 배수값은 순수 추정. |
| `k_perm` | 1.0e-4 /day | 회복 가능 손상 → 영구 손상 전환률. **7절(조기 vs 지연 수술) 결과를 이 파라미터 하나가 결정한다.** 131·129번 문헌이 방향을 지지하나 속도는 아니다. |
| `kflux_ab`, `kflux_tau`, `kdeg_csf` | — | 뇌→CSF 바이오마커 유출 상수. 9절의 예측 **방향**이 이 비율에 달려 있고, 본문에서 그렇게 명시했다. |
| `Hcol_cm` | 45 cmH₂O | 환자 신장에 비례하는 뇌실–복강 컬럼. 95·96번 문헌이 크기의 자릿수를 지지하지만 개별값은 환자마다 다르다. 3절의 지배항이며 **민감도 분석에서 2위 — 약물 다섯 개가 들어 있는 모델에서 환자의 키가 그 어떤 약보다 2년 후 보행에 더 중요하게 나왔다.** |
| `k_haz` | 1.5e-5 | 수액낭 부피 → 증상성 사건 위험률 변환. **따라서 5절의 경막하 발생률(%)은 서열척도이며 보정된 발생률이 아니다.** |
| `etv_dRout` | 0.06 | 소통성 수두증에서 ETV가 `Rout`를 거의 낮추지 못한다는 주장은 **구조적 논증**(저항이 stoma 하류에 있다)이며 측정값이 아니다. 112번 문헌의 임상 결과와 방향은 일치한다. |
| `az_tau_esc`, `az_f_esc` | 18 d, 0.78 | acetazolamide 만성 효과 소실. 기전(신장 산-염기 보상, 융모 수송체 상향조절)은 널리 인정되지만 iNPH에서 이 속도로 측정된 바 없다. 119번(DRAIN 시험)이 판정 근거가 된다. |

### 모델이 실제 문헌과 어긋나는 지점 (misses)

1. **반응률.** 이 모델은 `Rout`가 상승한 환자를 거의 모두 반응자로 만든다. 91번 문헌의 실측 반응률은 60–80%이고, 89·90번은 그보다도 회의적이다. → **모델은 선별되지 않은 집단에서 션트의 가치를 과대평가한다.**
2. **뇌실 크기.** 모델은 성공적 션트 후 Evans index를 0.01–0.02 감소시킨다. 48·49번 문헌은 Evans index가 임상 개선과 잘 대응하지 않으며 변화도 작다고 보고한다 — 모델은 여전히 크기를 압력에 과도하게 결합시키고 있다.
3. **Acetazolamide.** 8절은 급성 효과 39–43%, 만성 12%를 준다. 117·118번은 백질 병변 호전을 보고했으나 표본이 작고, 119번(무작위 시험)이 최종 판정자다. S12/S13은 **상한**으로 읽어야 한다.
4. **박동성의 주파수 대역.** 31번 문헌은 뇌 변형이 심박동보다 느린 vasomotion과 더 잘 연관된다고 보고한다. 이 모델은 `AMP`(심박동 성분)에 `kappa_tm`을 걸었다. 두 해석 중 어느 쪽이 옳은지 모델은 판정하지 않으며, 판정할 데이터도 갖고 있지 않다.
5. **단일 CSF 구획.** 66번 문헌은 뇌실 CSF와 요추 CSF의 바이오마커가 다르다고 보고한다. 이 모델은 하나의 CSF 구획만 갖고 있어 tap test(요추)와 션트(뇌실)를 같은 농도로 읽는다.

---

## 20. 도구·방법론

- mrgsolve (R): <https://mrgsolve.org/>
- Graphviz: <https://graphviz.org/>
- Shiny (R): <https://shiny.posit.co/>
- NCBI E-utilities (본 문헌 목록의 검증에 사용): <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

> **면책.** 본 모델과 문헌 목록은 교육·연구 목적입니다. 임상 의사결정,
> 처방, 규제 제출에 사용할 수 없습니다. 특히 5절의 밸브 적정 지도는
> 보정되지 않은 위험 계수를 포함하므로 **실제 밸브 설정의 근거가 될 수
> 없습니다.**
