# 방사선 유발 폐손상 (RILI) QSP 모델 — 참고문헌

# Radiation-Induced Lung Injury — References


> **PMID 검증 방법 (how these PMIDs were obtained).** 이 목록의 어떤 항목도
> 기억에서 작성되지 않았습니다. 각 절의 주제별 질의를 NCBI E-utilities
> (`esearch` + `esummary`)에 직접 보내고, PubMed 가 실제로 반환한 레코드의
> 저자·연도·저널·제목·PMID 를 그대로 옮겼습니다. 따라서 모든 PMID 링크는
> 해석됩니다. 질의는 관련도(relevance) 정렬을 사용했으므로 일부 항목은
> 주변 주제이며, 그런 경우에도 반환된 레코드를 임의로 바꾸지 않았습니다.
>
> 총 **236편**. 모델 파라미터의 근거는 `rili_reference_model.py` 의
> 해당 파라미터 주석에, 보정 목표는 `rili_calibration.py` 에 있습니다.


---

## 1. 선량-체적 관계와 방사선 폐렴 위험 (Dose-Volume Relationships & RP Risk)

MLD·V20·V5 가 폐렴 위험과 어떻게 연결되는지에 대한 임상 근거. 모델의 NTCP 보정 목표(QUANTEC: MLD 13/20/24 Gy → 10/20/30%)가 여기서 온다.

1. Rodrigues G et al. (2004). Prediction of radiation pneumonitis by dose - volume histogram parameters in lung cancer--a systematic review. *Radiother Oncol* 71:127-38. [PMID 15110445](https://pubmed.ncbi.nlm.nih.gov/15110445/)
2. Ma Y et al. (2025). Optimal dose-volume histogram thresholds for radiation pneumonitis prevention in lung cancer patients receiving immunotherapy. *Radiat Oncol* 20:60. [PMID 40264199](https://pubmed.ncbi.nlm.nih.gov/40264199/)
3. Zhang Y et al. (2025). Predicting Severe Radiation Pneumonitis in Patients With Locally-Advanced Non-Small-Cell Lung Cancer After Thoracic Radiotherapy: Development and Validation of a Nomogram Based on the Clinical, Hematological, and Dose-Volume Histogram Parameters. *Clin Lung Cancer* 26:393-406. [PMID 40087057](https://pubmed.ncbi.nlm.nih.gov/40087057/)
4. Kirakli EK et al. (2023). Ipsilateral lung dose as a correlative measure for radiation pneumonitis in patients treated with definitive concurrent radiochemotherapy. *J Cancer Res Ther* 19:1153-1159. [PMID 37787278](https://pubmed.ncbi.nlm.nih.gov/37787278/)
5. Milano MT et al. (2024). Comparison of Risks of Late Effects From Radiation Therapy in Children Versus Adults: Insights From the QUANTEC, HyTEC, and PENTEC Efforts. *Int J Radiat Oncol Biol Phys* 119:387-400. [PMID 38069917](https://pubmed.ncbi.nlm.nih.gov/38069917/)
6. Kirli Bolukbas M et al. (2020). Effect of lung volume on helical radiotherapy in esophageal cancer: are there predictive factors to achieve acceptable lung doses?. *Strahlenther Onkol* 196:805-812. [PMID 31980833](https://pubmed.ncbi.nlm.nih.gov/31980833/)
7. He R et al. (2023). Model development of dose and volume predictors for esophagitis induced during chemoradiotherapy for lung cancer as a step towards radiobiological treatment planning. *BMC Pulm Med* 23:379. [PMID 37814254](https://pubmed.ncbi.nlm.nih.gov/37814254/)
8. Kirli Bolukbas M et al. (2023). Cardiac protective techniques in left breast radiotherapy: rapid selection criteria for routine clinical decision making. *Eur J Med Res* 28:504. [PMID 37941070](https://pubmed.ncbi.nlm.nih.gov/37941070/)
9. Frédéric-Moreau T et al. (2026). Cumulative lung dose-volume predictors of radiation pneumonitis in thoracic reirradiation: A systematic scoping review. *Crit Rev Oncol Hematol* 222:105287. [PMID 41864297](https://pubmed.ncbi.nlm.nih.gov/41864297/)
10. Tucker SL et al. (2019). Validation of Effective Dose as a Better Predictor of Radiation Pneumonitis Risk Than Mean Lung Dose: Secondary Analysis of a Randomized Trial. *Int J Radiat Oncol Biol Phys* 103:403-410. [PMID 30291994](https://pubmed.ncbi.nlm.nih.gov/30291994/)
11. Flakus MJ et al. (2023). Metrics of dose to highly ventilated lung are predictive of radiation-induced pneumonitis in lung cancer patients. *Radiother Oncol* 182:109553. [PMID 36813178](https://pubmed.ncbi.nlm.nih.gov/36813178/)

---

## 2. NTCP 형식론과 체적효과 지수 (NTCP Formalism & the Volume-Effect Exponent)

LKB 모델·gEUD·기능적 하부단위(FSU) 병렬 구조. 모델은 n 을 입력하지 않고 300개 무작위 DVH에서 측정한다(§2 결과: n = 1.43 vs 문헌 0.99).

12. Wang Z et al. (2020). Lyman-Kutcher-Burman normal tissue complication probability modeling for radiation-induced esophagitis in non-small cell lung cancer patients receiving proton radiotherapy. *Radiother Oncol* 146:200-204. [PMID 32220701](https://pubmed.ncbi.nlm.nih.gov/32220701/)
13. Semenenko VA et al. (2008). Lyman-Kutcher-Burman NTCP model parameters for radiation pneumonitis and xerostomia based on combined analysis of published clinical data. *Phys Med Biol* 53:737-55. [PMID 18199912](https://pubmed.ncbi.nlm.nih.gov/18199912/)
14. Kavousi N et al. (2020). Ipsilateral lung normal tissue complication probability parameters for different dose calculation algorithms in radiotherapy of breast cancer. *J Cancer Res Ther* 16:1323-1330. [PMID 33342791](https://pubmed.ncbi.nlm.nih.gov/33342791/)
15. Ming X et al. (2024). Intensity-modulated proton and carbon-ion radiotherapy using a fixed-beam system for locally advanced lung cancer: dosimetric comparison with x-ray radiotherapy and normal tissue complication probability (NTCP) evaluation. *Phys Med Biol* 69. [PMID 38064747](https://pubmed.ncbi.nlm.nih.gov/38064747/)
16. Anetai Y et al. (2023). Effective optimization strategy for large optimization volume object, remaining volume at risk (RVR):α-value selection and usage from generalized equivalent uniform dose (gEUD) curve deviation perspective. *Phys Med Biol* 68. [PMID 36745933](https://pubmed.ncbi.nlm.nih.gov/36745933/)
17. Wu Q et al. (2005). Dose sculpting with generalized equivalent uniform dose. *Med Phys* 32:1387-96. [PMID 15984690](https://pubmed.ncbi.nlm.nih.gov/15984690/)
18. Wu Q et al. (2003). Intensity-modulated radiotherapy optimization with gEUD-guided dose-volume objectives. *Phys Med Biol* 48:279-91. [PMID 12608607](https://pubmed.ncbi.nlm.nih.gov/12608607/)
19. Haghbin A et al. (2022). Prediction of chronic kidney disease in abdominal cancers radiation therapy using the functional assays of normal tissue complication probability models. *J Cancer Res Ther* 18:718-724. [PMID 35900545](https://pubmed.ncbi.nlm.nih.gov/35900545/)
20. Jackson A et al. (1993). Probability of radiation-induced complications for normal tissues with parallel architecture subject to non-uniform irradiation. *Med Phys* 20:613-25. [PMID 8350812](https://pubmed.ncbi.nlm.nih.gov/8350812/)
21. Niemierko A et al. (1993). Modeling of normal tissue response to radiation: the critical volume model. *Int J Radiat Oncol Biol Phys* 25:135-45. [PMID 8416870](https://pubmed.ncbi.nlm.nih.gov/8416870/)
22. Gagliardi G et al. (2000). Radiation pneumonitis after breast cancer irradiation: analysis of the complication probability using the relative seriality model. *Int J Radiat Oncol Biol Phys* 46:373-81. [PMID 10661344](https://pubmed.ncbi.nlm.nih.gov/10661344/)
23. Tajiki S et al. (2023). A systematic review of the normal tissue complication probability models and parameters: Head and neck cancers treated with conformal radiotherapy. *Head Neck* 45:3146-3156. [PMID 37767820](https://pubmed.ncbi.nlm.nih.gov/37767820/)
24. Dell'Oro M et al. (2022). Normal tissue complication probability modeling to guide individual treatment planning in pediatric cranial proton and photon radiotherapy. *Med Phys* 49:742-755. [PMID 34796509](https://pubmed.ncbi.nlm.nih.gov/34796509/)
25. Cheraghi S et al. (2017). Normal tissue complication probability modeling of radiation-induced sensorineural hearing loss after head-and-neck radiation therapy. *Int J Radiat Biol* 93:1327-1333. [PMID 28967273](https://pubmed.ncbi.nlm.nih.gov/28967273/)

---

## 3. 선형-이차 모델과 분할 효과 (Linear-Quadratic Model & Fractionation)

폐 α/β ≈ 3 Gy, NSCLC α/β ≈ 10 Gy. 두 값의 차이가 저분할 조사에서 치료비를 어느 방향으로 움직이는지를 결정한다.

26. Bentzen SM et al. (2000). Quantitative clinical radiobiology of early and late lung reactions. *Int J Radiat Biol* 76:453-62. [PMID 10815624](https://pubmed.ncbi.nlm.nih.gov/10815624/)
27. Gay HA et al. (2009). Isodose-based methodology for minimizing the morbidity and mortality of thoracic hypofractionated radiotherapy. *Radiother Oncol* 91:369-78. [PMID 19058867](https://pubmed.ncbi.nlm.nih.gov/19058867/)
28. Van Dyk J et al. (1989). Radiation-induced lung damage: dose-time-fractionation considerations. *Radiother Oncol* 14:55-69. [PMID 2928557](https://pubmed.ncbi.nlm.nih.gov/2928557/)
29. Fowler JF et al. (2004). A challenge to traditional radiation oncology. *Int J Radiat Oncol Biol Phys* 60:1241-56. [PMID 15519797](https://pubmed.ncbi.nlm.nih.gov/15519797/)
30. Ohtakara K et al. (2023). 5-Fraction Re-radiosurgery for Progression Following 8-Fraction Radiosurgery of Brain Metastases From Lung Adenocarcinoma: Importance of Gross Tumor Coverage With Biologically Effective Dose ≥80 Gy and Internal Dose Increase. *Cureus* 15:e42299. [PMID 37609081](https://pubmed.ncbi.nlm.nih.gov/37609081/)
31. Denekamp J (1986). Cell kinetics and radiation biology. *Int J Radiat Biol Relat Stud Phys Chem Med* 49:357-80. [PMID 3510997](https://pubmed.ncbi.nlm.nih.gov/3510997/)
32. Marks LB et al. (2010). Radiation dose-volume effects in the lung. *Int J Radiat Oncol Biol Phys* 76:S70-6. [PMID 20171521](https://pubmed.ncbi.nlm.nih.gov/20171521/)
33. Klement RJ et al. (2020). Estimation of the α/β ratio of non-small cell lung cancer treated with stereotactic body radiotherapy. *Radiother Oncol* 142:210-216. [PMID 31431371](https://pubmed.ncbi.nlm.nih.gov/31431371/)
34. Lee P et al. (2021). Local Control After Stereotactic Body Radiation Therapy for Stage I Non-Small Cell Lung Cancer. *Int J Radiat Oncol Biol Phys* 110:160-171. [PMID 30954520](https://pubmed.ncbi.nlm.nih.gov/30954520/)
35. Chi A et al. (2013). What would be the most appropriate α/β ratio in the setting of stereotactic body radiation therapy for early stage non-small cell lung cancer. *Biomed Res Int* 2013:391021. [PMID 24350266](https://pubmed.ncbi.nlm.nih.gov/24350266/)
36. Zhou S (2026). A modified lethal-potentially lethal model for oxygen-mediated FLASH sparing in stem cell niches. *Med Phys* 53:e70469. [PMID 42108222](https://pubmed.ncbi.nlm.nih.gov/42108222/)
37. Wein LM et al. (2000). Dynamic optimization of a linear-quadratic model with incomplete repair and volume-dependent sensitivity and repopulation. *Int J Radiat Oncol Biol Phys* 47:1073-83. [PMID 10863081](https://pubmed.ncbi.nlm.nih.gov/10863081/)
38. Millar WT et al. (2007). Effects of very low dose-rate (90)Sr/(90)Y exposure on the acute moist desquamation response of pig skin. *Radiother Oncol* 83:187-95. [PMID 17467835](https://pubmed.ncbi.nlm.nih.gov/17467835/)

---

## 4. 폐포 상피·II형 폐포세포·계면활성제 (Alveolar Epithelium, AT2 & Surfactant)

치사손상 세포가 분열을 시도할 때까지 계속 기능한다는 사실이 이 모델의 잠재기(4–12주)를 만든다. AT2 회전율·KL-6·SP-D 근거.

39. Liu X et al. (2021). β-Catenin/Lin28/let-7 regulatory network determines type II alveolar epithelial stem cell differentiation phenotypes following thoracic irradiation. *J Radiat Res* 62:119-132. [PMID 33302295](https://pubmed.ncbi.nlm.nih.gov/33302295/)
40. Coggle JE et al. (1986). Radiation effects in the lung. *Environ Health Perspect* 70:261-91. [PMID 3549278](https://pubmed.ncbi.nlm.nih.gov/3549278/)
41. Beach TA et al. (2023). Epithelial Responses in Radiation-Induced Lung Injury (RILI) Allow Chronic Inflammation and Fibrogenesis. *Radiat Res* 199:439-451. [PMID 37237442](https://pubmed.ncbi.nlm.nih.gov/37237442/)
42. Finkelstein JN (1990). Physiologic and toxicologic responses of alveolar type II cells. *Toxicology* 60:41-52. [PMID 2180132](https://pubmed.ncbi.nlm.nih.gov/2180132/)
43. Kohno N (1999). Serum marker KL-6/MUC1 for the diagnosis and management of interstitial pneumonitis. *J Med Invest* 46:151-8. [PMID 10687309](https://pubmed.ncbi.nlm.nih.gov/10687309/)
44. Goto K et al. (2001). Serum levels of KL-6 are useful biomarkers for severe radiation pneumonitis. *Lung Cancer* 34:141-8. [PMID 11557124](https://pubmed.ncbi.nlm.nih.gov/11557124/)
45. Hara R et al. (2004). Serum levels of KL-6 for predicting the occurrence of radiation pneumonitis after stereotactic radiotherapy for lung tumors. *Chest* 125:340-4. [PMID 14718465](https://pubmed.ncbi.nlm.nih.gov/14718465/)
46. Mahmoud Hewala TI et al. (2024). Dosimetry and Biochemical Comparison of Early Radiation-Induced Lung Toxicity in Breast Cancer Patients Treated with 3D-CRT and IMRT: the Role of Serum Interleukin-6 and Pulmonary Surfactant Protein-D. *Asian Pac J Cancer Prev* 25:1707-1713. [PMID 38809643](https://pubmed.ncbi.nlm.nih.gov/38809643/)
47. Malaviya R et al. (2015). Radiation-induced lung injury and inflammation in mice: role of inducible nitric oxide synthase and surfactant protein D. *Toxicol Sci* 144:27-38. [PMID 25552309](https://pubmed.ncbi.nlm.nih.gov/25552309/)
48. Xu L et al. (2019). Genetic variants of SP-D confer susceptibility to radiation pneumonitis in lung cancer patients undergoing thoracic radiation therapy. *Cancer Med* 8:2599-2611. [PMID 30897289](https://pubmed.ncbi.nlm.nih.gov/30897289/)
49. Banerjee ER et al. (2012). Characterization of lung stem cell niches in a mouse model of bleomycin-induced fibrosis. *Stem Cell Res Ther* 3:21. [PMID 22643035](https://pubmed.ncbi.nlm.nih.gov/22643035/)
50. Mukhametshina RT et al. (2013). Quantitative proteome analysis of alveolar type-II cells reveals a connection of integrin receptor subunits beta 2/6 and WNT signaling. *J Proteome Res* 12:5598-608. [PMID 24175614](https://pubmed.ncbi.nlm.nih.gov/24175614/)
51. Yang J et al. (2025). Fibrotic extracellular matrix microenvironment induces alveolar type II epithelial cell senescence via integrin-β1/FAK/YAP signaling pathway. *Int Immunopharmacol* 162:115169. [PMID 40639050](https://pubmed.ncbi.nlm.nih.gov/40639050/)
52. Han L et al. (2026). TMEM131-Mediated Soluble TRAIL Triggered Type II Alveolar Epithelial Cell Senescence in Radiation-Induced Lung Injury. *Adv Sci (Weinh)* 13:e09973. [PMID 41293955](https://pubmed.ncbi.nlm.nih.gov/41293955/)
53. Citrin DE et al. (2013). Role of type II pneumocyte senescence in radiation-induced lung fibrosis. *J Natl Cancer Inst* 105:1474-84. [PMID 24052614](https://pubmed.ncbi.nlm.nih.gov/24052614/)

---

## 5. 미세혈관 내피 손상과 모세혈관 소실 (Microvascular Endothelium & Capillary Rarefaction)

ASMase–ceramide 매개 내피 세포사, 내피 전구세포 소진, 영구적 수복 천장 EC_max = 1/(1+KECIRR·BED) 의 근거.

54. Boittin FX et al. (2024). The Protein Kinase A Inhibitor KT5720 Prevents Endothelial Dysfunctions Induced by High-Dose Irradiation. *Int J Mol Sci* 25. [PMID 38396945](https://pubmed.ncbi.nlm.nih.gov/38396945/)
55. Wynn TA (2008). Cellular and molecular mechanisms of fibrosis. *J Pathol* 214:199-210. [PMID 18161745](https://pubmed.ncbi.nlm.nih.gov/18161745/)
56. Choi KJ et al. (2020). Endothelial-to-mesenchymal transition in anticancer therapy and normal tissue damage. *Exp Mol Med* 52:781-792. [PMID 32467609](https://pubmed.ncbi.nlm.nih.gov/32467609/)
57. Nam JK et al. (2025). Combined HIF-1α blockade and CHIR99021 treatment reverses pulmonary fibrosis via modulation endothelial-to-mesenchymal transition. *iScience* 28:114028. [PMID 41399517](https://pubmed.ncbi.nlm.nih.gov/41399517/)
58. Abe JI et al. (2024). Radiation-Induced Macrovessel/Microvessel Disease. *Arterioscler Thromb Vasc Biol* 44:2407-2415. [PMID 39445428](https://pubmed.ncbi.nlm.nih.gov/39445428/)
59. Peterson LM et al. (1992). Vascular response to radiation injury in the rat lung. *Radiat Res* 129:139-48. [PMID 1734443](https://pubmed.ncbi.nlm.nih.gov/1734443/)
60. Jagtap J et al. (2021). A rapid dynamic in vivo near-infrared fluorescence imaging assay to track lung vascular permeability after acute radiation injury. *Am J Physiol Lung Cell Mol Physiol* 320:L436-L450. [PMID 33404364](https://pubmed.ncbi.nlm.nih.gov/33404364/)
61. Chatterjee S et al. (2019). LGM2605 Reduces Space Radiation-Induced NLRP3 Inflammasome Activation and Damage in In Vitro Lung Vascular Networks. *Int J Mol Sci* 20. [PMID 30621290](https://pubmed.ncbi.nlm.nih.gov/30621290/)
62. Leonetti D et al. (2020). Secretion of Acid Sphingomyelinase and Ceramide by Endothelial Cells Contributes to Radiation-Induced Intestinal Toxicity. *Cancer Res* 80:2651-2662. [PMID 32291318](https://pubmed.ncbi.nlm.nih.gov/32291318/)
63. Kolesnick R et al. (2003). Radiation and ceramide-induced apoptosis. *Oncogene* 22:5897-906. [PMID 12947396](https://pubmed.ncbi.nlm.nih.gov/12947396/)
64. Niaudet C et al. (2017). Plasma membrane reorganization links acid sphingomyelinase/ceramide to p38 MAPK pathways in endothelial cells apoptosis. *Cell Signal* 33:10-21. [PMID 28179144](https://pubmed.ncbi.nlm.nih.gov/28179144/)
65. Adzraku SY et al. (2024). Endothelial Robo4 suppresses endothelial-to-mesenchymal transition induced by irradiation and improves hematopoietic reconstitution. *Cell Death Dis* 15:159. [PMID 38383474](https://pubmed.ncbi.nlm.nih.gov/38383474/)
66. Turchan WT et al. (2016). Irradiated human endothelial progenitor cells induce bystander killing in human non-small cell lung and pancreatic cancer cells. *Int J Radiat Biol* 92:427-33. [PMID 27258472](https://pubmed.ncbi.nlm.nih.gov/27258472/)
67. Li Y et al. (2024). Silencing endomucin in bone marrow sinusoids improves hematopoietic stem and progenitor cell homing during transplantation. *Stem Cells* 42:889-901. [PMID 38995653](https://pubmed.ncbi.nlm.nih.gov/38995653/)

---

## 6. DAMP 방출과 선천면역 감지 (DAMP Release & Innate Sensing)

HMGB1 · cGAS–STING · NLRP3. 모델에서 DAMP 원천은 선량이 아니라 사멸 플럭스 DTH = KMIT × DOOM 이다.

68. Ge X et al. (2024). STING facilitates the development of radiation-induced lung injury via regulating the PERK/eIF2α pathway. *Transl Lung Cancer Res* 13:3010-3025. [PMID 39670000](https://pubmed.ncbi.nlm.nih.gov/39670000/)
69. Yang C et al. (2023). Role of the cGAS-STING pathway in radiotherapy for non-small cell lung cancer. *Radiat Oncol* 18:145. [PMID 37667279](https://pubmed.ncbi.nlm.nih.gov/37667279/)
70. Zhang Y et al. (2023). STING-Dependent Sensing of Self-DNA Driving Pyroptosis Contributes to Radiation-Induced Lung Injury. *Int J Radiat Oncol Biol Phys* 117:928-941. [PMID 37230431](https://pubmed.ncbi.nlm.nih.gov/37230431/)
71. Zhang M et al. (2025). NLRP3 inflammasome mediates pyroptosis of alveolar macrophages to induce radiation lung injury. *J Hazard Mater* 484:136740. [PMID 39642726](https://pubmed.ncbi.nlm.nih.gov/39642726/)
72. Jiang YC et al. (2025). Raspberry ketone alleviates radiation-induced lung injury through the STAT2-P2X7r/NLRP3 signaling pathway. *Phytomedicine* 145:156984. [PMID 40544736](https://pubmed.ncbi.nlm.nih.gov/40544736/)
73. Wu X et al. (2019). Melatonin Alleviates Radiation-Induced Lung Injury via Regulation of miR-30e/NLRP3 Axis. *Oxid Med Cell Longev* 2019:4087298. [PMID 30755784](https://pubmed.ncbi.nlm.nih.gov/30755784/)

---

## 7. 사이토카인 네트워크와 빠른 고리 (Cytokine Network & the Fast Loop)

TNF-α · IL-6 · IL-1β · TGF-β1 의 혈중 농도와 폐렴 위험. 빠른 고리의 이득 예산(0.072 < KCE 0.150)이 폐렴이 자연 소실되는 이유다.

74. Chen Y et al. (2001). Circulating IL-6 as a predictor of radiation pneumonitis. *Int J Radiat Oncol Biol Phys* 49:641-8. [PMID 11172944](https://pubmed.ncbi.nlm.nih.gov/11172944/)
75. Rübe CE et al. (2008). Cytokine plasma levels: reliable predictors for radiation pneumonitis?. *PLoS One* 3:e2898. [PMID 18682839](https://pubmed.ncbi.nlm.nih.gov/18682839/)
76. Stenmark MH et al. (2012). Combining physical and biologic parameters to predict radiation-induced lung toxicity in patients with non-small-cell lung cancer treated with definitive radiation therapy. *Int J Radiat Oncol Biol Phys* 84:e217-22. [PMID 22935395](https://pubmed.ncbi.nlm.nih.gov/22935395/)
77. Anscher MS et al. (1998). Plasma transforming growth factor beta1 as a predictor of radiation pneumonitis. *Int J Radiat Oncol Biol Phys* 41:1029-35. [PMID 9719112](https://pubmed.ncbi.nlm.nih.gov/9719112/)
78. Zhou X et al. (2022). 3,3'-Diindolylmethane attenuates inflammation and fibrosis in radiation-induced lung injury by regulating NF-κB/TGF-β/Smad signaling pathways. *Exp Lung Res* 48:103-113. [PMID 35594367](https://pubmed.ncbi.nlm.nih.gov/35594367/)
79. Chen G et al. (2023). Cepharanthine Ameliorates Pulmonary Fibrosis by Inhibiting the NF-κB/NLRP3 Pathway, Fibroblast-to-Myofibroblast Transition and Inflammation. *Molecules* 28. [PMID 36677811](https://pubmed.ncbi.nlm.nih.gov/36677811/)
80. Gandhi KA et al. (2019). Oral administration of 3,3'-diselenodipropionic acid prevents thoracic radiation induced pneumonitis in mice by suppressing NF-kB/IL-17/G-CSF/neutrophil axis. *Free Radic Biol Med* 145:8-19. [PMID 31521664](https://pubmed.ncbi.nlm.nih.gov/31521664/)

---

## 8. TGF-β1 활성화와 근섬유아세포 (TGF-β1 Activation & the Myofibroblast)

잠재형 TGF-β1 은 integrin αvβ6 가 전달하는 힘으로 풀린다 — 강성이 문턱을 넘어야 활성화된다는 점이 이 모델의 이중안정성의 근거.

81. Munger JS et al. (1999). The integrin alpha v beta 6 binds and activates latent TGF beta 1: a mechanism for regulating pulmonary inflammation and fibrosis. *Cell* 96:319-28. [PMID 10025398](https://pubmed.ncbi.nlm.nih.gov/10025398/)
82. Puthawala K et al. (2008). Inhibition of integrin alpha(v)beta6, an activator of latent transforming growth factor-beta, prevents radiation-induced lung fibrosis. *Am J Respir Crit Care Med* 177:82-90. [PMID 17916808](https://pubmed.ncbi.nlm.nih.gov/17916808/)
83. Sheppard D (2015). Epithelial-mesenchymal interactions in fibrosis and repair. Transforming growth factor-β activation by epithelial cells and fibroblasts. *Ann Am Thorac Soc* 12 Suppl 1:S21-3. [PMID 25830829](https://pubmed.ncbi.nlm.nih.gov/25830829/)
84. Ding H et al. (2023). ROS-responsive microneedles loaded with integrin avβ6-blocking antibodies for the treatment of pulmonary fibrosis. *J Control Release* 360:365-375. [PMID 37331606](https://pubmed.ncbi.nlm.nih.gov/37331606/)
85. Santos A et al. (2018). Matrix Stiffness: the Conductor of Organ Fibrosis. *Curr Rheumatol Rep* 20:2. [PMID 29349703](https://pubmed.ncbi.nlm.nih.gov/29349703/)
86. Hinz B (2009). Tissue stiffness, latent TGF-beta1 activation, and mechanical signal transduction: implications for the pathogenesis and treatment of fibrosis. *Curr Rheumatol Rep* 11:120-6. [PMID 19296884](https://pubmed.ncbi.nlm.nih.gov/19296884/)
87. Vissers G et al. (2024). The role of fibrosis in endometriosis: a systematic review. *Hum Reprod Update* 30:706-750. [PMID 39067455](https://pubmed.ncbi.nlm.nih.gov/39067455/)
88. Freeberg MAT et al. (2025). Piezo2 Is a Key Mechanoreceptor in Lung Fibrosis that Drives Myofibroblast Differentiation. *Am J Pathol* 195:626-638. [PMID 39855300](https://pubmed.ncbi.nlm.nih.gov/39855300/)
89. Yue B et al. (2024). SPP1 induces idiopathic pulmonary fibrosis and NSCLC progression via the PI3K/Akt/mTOR pathway. *Respir Res* 25:362. [PMID 39369217](https://pubmed.ncbi.nlm.nih.gov/39369217/)
90. Tu J et al. (2024). Nintedanib Mitigates Radiation-Induced Pulmonary Fibrosis by Suppressing Epithelial Cell Inflammatory Response and Inhibiting Fibroblast-to-Myofibroblast Transition. *Int J Biol Sci* 20:3353-3371. [PMID 38993568](https://pubmed.ncbi.nlm.nih.gov/38993568/)
91. Chen F et al. (2025). The pivotal role of TGF-β/Smad pathway in fibrosis pathogenesis and treatment. *Front Oncol* 15:1649179. [PMID 40969268](https://pubmed.ncbi.nlm.nih.gov/40969268/)
92. Kim H et al. (2020). LXA(4)-FPR2 signaling regulates radiation-induced pulmonary fibrosis via crosstalk with TGF-β/Smad signaling. *Cell Death Dis* 11:653. [PMID 32811815](https://pubmed.ncbi.nlm.nih.gov/32811815/)
93. Ma HY et al. (2023). LOXL4, but not LOXL2, is the critical determinant of pathological collagen cross-linking and fibrosis in the lung. *Sci Adv* 9:eadf0133. [PMID 37235663](https://pubmed.ncbi.nlm.nih.gov/37235663/)
94. Chen W et al. (2020). Lysyl Oxidase (LOX) Family Members: Rationale and Their Potential as Therapeutic Targets for Liver Fibrosis. *Hepatology* 72:729-741. [PMID 32176358](https://pubmed.ncbi.nlm.nih.gov/32176358/)
95. Peng L et al. (2025). NUDT21 regulates lysyl oxidase-like 2(LOXL2) to influence ECM protein cross-linking in silica-induced pulmonary fibrosis. *Ecotoxicol Environ Saf* 290:117572. [PMID 39700768](https://pubmed.ncbi.nlm.nih.gov/39700768/)

---

## 9. 이중안정성·역학적 기억·수학적 모델 (Bistability, Mechanical Memory & Models)

느린 고리의 세 고정점(정상 1.000 · 분리선 1.586 · 섬유화 2.9)과 LOX 가교결합이 만드는 비가역성의 근거.

96. Mu M et al. (2026). Multi-omics reveals a novel Cxcr4(+) subpopulation of alveolar macrophages and therapeutic effect of AMD3100 in mice with advanced silicosis. *Clin Transl Med* 16:e70705. [PMID 42204838](https://pubmed.ncbi.nlm.nih.gov/42204838/)
97. Xiong W et al. (2003). A positive-feedback-based bistable 'memory module' that governs a cell fate decision. *Nature* 426:460-5. [PMID 14647386](https://pubmed.ncbi.nlm.nih.gov/14647386/)
98. Tian XJ et al. (2019). Modeling ncRNA-Mediated Circuits in Cell Fate Decision. *Methods Mol Biol* 1912:411-426. [PMID 30635903](https://pubmed.ncbi.nlm.nih.gov/30635903/)
99. Hsu C et al. (2016). Protein Dimerization Generates Bistability in Positive Feedback Loops. *Cell Rep* 16:1204-1210. [PMID 27425609](https://pubmed.ncbi.nlm.nih.gov/27425609/)
100. Hari K et al. (2022). Emergent properties of coupled bistable switches. *J Biosci* 47. [PMID 36550692](https://pubmed.ncbi.nlm.nih.gov/36550692/)
101. Larreta-Garde V et al. (2002). Modeling extracellular matrix degradation balance with proteinase/transglutaminase cycle. *J Theor Biol* 217:105-24. [PMID 12183135](https://pubmed.ncbi.nlm.nih.gov/12183135/)

---

## 10. 코르티코스테로이드 치료 (Corticosteroid Treatment)

빠른 고리에만 도달하는 약. 모델은 증상 지수를 낮추지만 섬유화 전환 체적을 전혀 바꾸지 못한다고 예측한다.

102. Arroyo-Hernández M et al. (2021). Radiation-induced lung injury: current evidence. *BMC Pulm Med* 21:9. [PMID 33407290](https://pubmed.ncbi.nlm.nih.gov/33407290/)
103. Doshita K et al. (2023). Incidence and Treatment Outcome of Radiation Pneumonitis in Patients With Limited-stage Small Cell Lung Cancer Treated With Concurrent Accelerated Hyperfractionated Radiation Therapy and Chemotherapy. *Adv Radiat Oncol* 8:101129. [PMID 36845617](https://pubmed.ncbi.nlm.nih.gov/36845617/)
104. Forschner A et al. (2014). Radiation recall dermatitis and radiation pneumonitis during treatment with vemurafenib. *Melanoma Res* 24:512-6. [PMID 24743051](https://pubmed.ncbi.nlm.nih.gov/24743051/)
105. Muraoka T et al. (2002). Corticosteroid refractory radiation pneumonitis that remarkably responded to cyclosporin A. *Intern Med* 41:730-3. [PMID 12322802](https://pubmed.ncbi.nlm.nih.gov/12322802/)
106. Zafar B et al. (2026). Pulmonary Complications of Cancer Therapy: Clinical Presentations, Imaging Patterns, and Management Strategies. *Medicina (Kaunas)* 62. [PMID 41901659](https://pubmed.ncbi.nlm.nih.gov/41901659/)
107. Naidoo J et al. (2020). Chronic immune checkpoint inhibitor pneumonitis. *J Immunother Cancer* 8. [PMID 32554618](https://pubmed.ncbi.nlm.nih.gov/32554618/)
108. Zhao J et al. (2026). Interim Positron Emission Tomography-Guided Dose-Adapted Residual Site Radiation Therapy Improves Survival in Diffuse Large B-Cell Lymphoma Patients With Partial Metabolic Response After R-CHOP: A Retrospective Cohort Analysis. *Adv Radiat Oncol* 11:102078. [PMID 42440462](https://pubmed.ncbi.nlm.nih.gov/42440462/)
109. Tsung I et al. (2021). A Pilot Study of Checkpoint Inhibitors in Solid Organ Transplant Recipients with Metastatic Cutaneous Squamous Cell Carcinoma. *Oncologist* 26:133-138. [PMID 32969143](https://pubmed.ncbi.nlm.nih.gov/32969143/)

---

## 11. 아미포스틴과 방사선 보호 (Amifostine & Radioprotection)

WR-1065 혈장 반감기 ≈ 8분. 보호 = PFAMI·exp(−ln2·Δt/8min) 이라는 시간적 동시성 요건이 상반된 시험 결과의 기전적 설명이 된다.

110. Zhang XJ et al. (2012). Prediction of radiation pneumonitis in lung cancer patients: a systematic review. *J Cancer Res Clin Oncol* 138:2103-16. [PMID 22842662](https://pubmed.ncbi.nlm.nih.gov/22842662/)
111. Antonadou D et al. (2001). Randomized phase III trial of radiation treatment +/- amifostine in patients with advanced-stage lung cancer. *Int J Radiat Oncol Biol Phys* 51:915-22. [PMID 11704311](https://pubmed.ncbi.nlm.nih.gov/11704311/)
112. Wang S et al. (2012). [Effect of amifostine on locally advanced non-small cell lung cancer patients treated with radiotherapy: a meta-analysis of randomized controlled trials]. *Zhongguo Fei Ai Za Zhi* 15:539-44. [PMID 22989457](https://pubmed.ncbi.nlm.nih.gov/22989457/)
113. Komaki R et al. (2002). Randomized phase III study of chemoradiation with or without amifostine for patients with favorable performance status inoperable stage II-III non-small cell lung cancer: preliminary results. *Semin Radiat Oncol* 12:46-9. [PMID 11917284](https://pubmed.ncbi.nlm.nih.gov/11917284/)
114. Lawrence YR et al. (2013). The addition of amifostine to carboplatin and paclitaxel based chemoradiation in locally advanced non-small cell lung cancer: long-term follow-up of Radiation Therapy Oncology Group (RTOG) randomized trial 9801. *Lung Cancer* 80:298-305. [PMID 23477890](https://pubmed.ncbi.nlm.nih.gov/23477890/)
115. Werner-Wasik M et al. (2004). Randomized phase II study of amifostine mucosal protection by either subcutaneous injection or rapid IV bolus for patients with inoperable stage II-IIIA/B or stage IV non-small cell lung cancer with oligometastases receiving concurrent radiochemotherapy with carboplatin and paclitaxel followed by optional consolidative chemotherapy: a follow-up study after RTOG 98-01. *Semin Oncol* 31:47-51. [PMID 15726523](https://pubmed.ncbi.nlm.nih.gov/15726523/)
116. Cassatt DR et al. (2003). Effects of dose and schedule on the efficacy of ethyol: preclinical studies. *Semin Oncol* 30:31-9. [PMID 14727238](https://pubmed.ncbi.nlm.nih.gov/14727238/)
117. van der Vijgh WJ et al. (1996). Amifostine (Ethyol): pharmacokinetic and pharmacodynamic effects in vivo. *Eur J Cancer* 32A Suppl 4:S26-30. [PMID 8976819](https://pubmed.ncbi.nlm.nih.gov/8976819/)
118. Cassatt DR et al. (2002). Preclinical studies on the radioprotective efficacy and pharmacokinetics of subcutaneously administered amifostine. *Semin Oncol* 29:2-8. [PMID 12577236](https://pubmed.ncbi.nlm.nih.gov/12577236/)

---

## 12. SOD 모방체·아바소파셈 (SOD Mimetics & Avasopasem)

초과산화물 불균화. 아미포스틴과 같은 시간 창 논리를 따르지만 반감기가 27분으로 더 길다.

119. Sonis ST (2021). Superoxide Dismutase as an Intervention for Radiation Therapy-Associated Toxicities: Review and Profile of Avasopasem Manganese as a Treatment Option for Radiation-Induced Mucositis. *Drug Des Devel Ther* 15:1021-1029. [PMID 33716500](https://pubmed.ncbi.nlm.nih.gov/33716500/)
120. Sonis ST et al. (2023). Avasopasem for the treatment of radiotherapy-induced severe oral mucositis. *Expert Opin Investig Drugs* 32:463-470. [PMID 37365149](https://pubmed.ncbi.nlm.nih.gov/37365149/)
121. Mapuskar KA et al. (2023). Avasopasem manganese (GC4419) protects against cisplatin-induced chronic kidney disease: An exploratory analysis of renal metrics from a randomized phase 2b clinical trial in head and neck cancer patients. *Redox Biol* 60:102599. [PMID 36640725](https://pubmed.ncbi.nlm.nih.gov/36640725/)
122. Cui Z et al. (2025). Radioprotection redefined: drug discovery at the intersection of tardigrade biology and translational pharmacology. *Front Pharmacol* 16:1713914. [PMID 41341030](https://pubmed.ncbi.nlm.nih.gov/41341030/)
123. Ashcraft KA et al. (2015). Novel Manganese-Porphyrin Superoxide Dismutase-Mimetic Widens the Therapeutic Margin in a Preclinical Head and Neck Cancer Model. *Int J Radiat Oncol Biol Phys* 93:892-900. [PMID 26530759](https://pubmed.ncbi.nlm.nih.gov/26530759/)
124. Thompson JS et al. (2010). The manganese superoxide dismutase mimetic, M40403, protects adult mice from lethal total body irradiation. *Free Radic Res* 44:529-40. [PMID 20298121](https://pubmed.ncbi.nlm.nih.gov/20298121/)
125. Moulder JE et al. (2007). Future strategies for mitigation and treatment of chronic radiation-induced normal tissue injury. *Semin Radiat Oncol* 17:141-8. [PMID 17395044](https://pubmed.ncbi.nlm.nih.gov/17395044/)
126. Anderson CM et al. (2022). Two-Year Tumor Outcomes of a Phase 2B, Randomized, Double-Blind Trial of Avasopasem Manganese (GC4419) Versus Placebo to Reduce Severe Oral Mucositis Owing to Concurrent Radiation Therapy and Cisplatin for Head and Neck Cancer. *Int J Radiat Oncol Biol Phys* 114:416-421. [PMID 35724774](https://pubmed.ncbi.nlm.nih.gov/35724774/)
127. Anderson C et al. (2025). Avasopasem manganese treatment for severe oral mucositis from chemoradiotherapy for locally advanced head and neck cancer: phase 3 randomized controlled trial (ROMAN). *EClinicalMedicine* 89:103539. [PMID 41127563](https://pubmed.ncbi.nlm.nih.gov/41127563/)
128. Anderson CM et al. (2019). Phase IIb, Randomized, Double-Blind Trial of GC4419 Versus Placebo to Reduce Severe Oral Mucositis Due to Concurrent Radiotherapy and Cisplatin For Head and Neck Cancer. *J Clin Oncol* 37:3256-3265. [PMID 31618127](https://pubmed.ncbi.nlm.nih.gov/31618127/)

---

## 13. 항섬유화제·ACE 억제제 (Antifibrotics & ACE Inhibitors)

피르페니돈·닌테다닙·ACE 억제제·펜톡시필린. 모두 느린 고리에만 도달하므로 폐렴 엔드포인트로는 검출될 수 없다.

129. Ying H et al. (2021). Pirfenidone modulates macrophage polarization and ameliorates radiation-induced lung fibrosis by inhibiting the TGF-β1/Smad3 pathway. *J Cell Mol Med* 25:8662-8675. [PMID 34327818](https://pubmed.ncbi.nlm.nih.gov/34327818/)
130. Hou Z et al. (2025). Pirfenidone for grade 2 and grade 3 radiation-induced lung injury: a multicentre, open-label, randomised, phase 2 trial. *Lancet Oncol* 26:1552-1562. [PMID 41207313](https://pubmed.ncbi.nlm.nih.gov/41207313/)
131. Chen C et al. (2022). Pirfenidone for the prevention of radiation-induced lung injury in patients with locally advanced oesophageal squamous cell carcinoma: a protocol for a randomised controlled trial. *BMJ Open* 12:e060619. [PMID 36302570](https://pubmed.ncbi.nlm.nih.gov/36302570/)
132. Park HR et al. (2025). A human lung organoid platform for studying radiation-induced pulmonary fibrosis and antifibrotic drug screening. *Sci Rep* 16:1905. [PMID 41466095](https://pubmed.ncbi.nlm.nih.gov/41466095/)
133. Dabholkar S et al. (2022). Nintedanib-A case of treating concurrent idiopathic pulmonary fibrosis and non-small cell lung cancer. *Respirol Case Rep* 10:e0902. [PMID 35059200](https://pubmed.ncbi.nlm.nih.gov/35059200/)
134. De Ruysscher D et al. (2017). Nintedanib reduces radiation-induced microscopic lung fibrosis but this cannot be monitored by CT imaging: A preclinical study with a high precision image-guided irradiator. *Radiother Oncol* 124:482-487. [PMID 28774597](https://pubmed.ncbi.nlm.nih.gov/28774597/)
135. Moore ZR et al. (2024). Biomarkers associated with pulmonary exacerbations in a randomized trial of nintedanib for radiation pneumonitis. *Radiother Oncol* 196:110320. [PMID 38740091](https://pubmed.ncbi.nlm.nih.gov/38740091/)
136. Harder EM et al. (2015). Angiotensin-converting enzyme inhibitors decrease the risk of radiation pneumonitis after stereotactic body radiation therapy. *Pract Radiat Oncol* 5:e643-9. [PMID 26412341](https://pubmed.ncbi.nlm.nih.gov/26412341/)
137. Sun F et al. (2018). Angiotensin-converting Enzyme Inhibitors Decrease the Incidence of Radiation-induced Pneumonitis Among Lung Cancer Patients: A Systematic Review and Meta-analysis. *J Cancer* 9:2123-2131. [PMID 29937931](https://pubmed.ncbi.nlm.nih.gov/29937931/)
138. Wang LW et al. (2000). Can angiotensin-converting enzyme inhibitors protect against symptomatic radiation pneumonitis?. *Radiat Res* 153:405-10. [PMID 10761000](https://pubmed.ncbi.nlm.nih.gov/10761000/)
139. Medhora M et al. (2012). Radiation damage to the lung: mitigation by angiotensin-converting enzyme (ACE) inhibitors. *Respirology* 17:66-71. [PMID 22023053](https://pubmed.ncbi.nlm.nih.gov/22023053/)
140. Kaidar-Person O et al. (2018). Pentoxifylline and vitamin E for treatment or prevention of radiation-induced fibrosis in patients with breast cancer. *Breast J* 24:816-819. [PMID 29687536](https://pubmed.ncbi.nlm.nih.gov/29687536/)
141. Chiao TB et al. (2005). Role of pentoxifylline and vitamin E in attenuation of radiation-induced fibrosis. *Ann Pharmacother* 39:516-22. [PMID 15701781](https://pubmed.ncbi.nlm.nih.gov/15701781/)
142. Harpsø M et al. (2025). Pentoxifylline and vitamin E for treating radiation-induced fibrosis in breast and head and neck cancer patients. *Strahlenther Onkol* 201:1044-1048. [PMID 40542130](https://pubmed.ncbi.nlm.nih.gov/40542130/)

---

## 14. 면역관문 억제제와 방사선 폐렴 (Checkpoint Inhibitors & Pneumonitis)

PACIFIC 공고요법. 모델에서 더발루맙은 항을 더하지 않고 빠른 고리의 이득 GIMM 을 올리므로, 초과 위험이 MLD 에 따라 커진다.

143. Antonia SJ et al. (2017). Durvalumab after Chemoradiotherapy in Stage III Non-Small-Cell Lung Cancer. *N Engl J Med* 377:1919-1929. [PMID 28885881](https://pubmed.ncbi.nlm.nih.gov/28885881/)
144. Garassino MC et al. (2022). Durvalumab After Sequential Chemoradiotherapy in Stage III, Unresectable NSCLC: The Phase 2 PACIFIC-6 Trial. *J Thorac Oncol* 17:1415-1427. [PMID 35961520](https://pubmed.ncbi.nlm.nih.gov/35961520/)
145. Bradley JD et al. (2025). Simultaneous Durvalumab and Platinum-Based Chemoradiotherapy in Unresectable Stage III Non-Small Cell Lung Cancer: The Phase III PACIFIC-2 Study. *J Clin Oncol* 43:3610-3621. [PMID 41082707](https://pubmed.ncbi.nlm.nih.gov/41082707/)
146. Girard N et al. (2023). Treatment Characteristics and Real-World Progression-Free Survival in Patients With Unresectable Stage III NSCLC Who Received Durvalumab After Chemoradiotherapy: Findings From the PACIFIC-R Study. *J Thorac Oncol* 18:181-193. [PMID 36307040](https://pubmed.ncbi.nlm.nih.gov/36307040/)
147. Zhou P et al. (2022). Risk Factors for Immune Checkpoint Inhibitor-Related Pneumonitis in Cancer Patients: A Systemic Review and Meta-Analysis. *Respiration* 101:1035-1050. [PMID 36108598](https://pubmed.ncbi.nlm.nih.gov/36108598/)
148. Chen X et al. (2021). Radiation Versus Immune Checkpoint Inhibitor Associated Pneumonitis: Distinct Radiologic Morphologies. *Oncologist* 26:e1822-e1832. [PMID 34251728](https://pubmed.ncbi.nlm.nih.gov/34251728/)
149. Chao Y et al. (2022). Risk factors for immune checkpoint inhibitor-related pneumonitis in non-small cell lung cancer. *Transl Lung Cancer Res* 11:295-306. [PMID 35280322](https://pubmed.ncbi.nlm.nih.gov/35280322/)
150. Tan P et al. (2023). Risk Factors for Refractory Immune Checkpoint Inhibitor-related Pneumonitis in Patients With Lung Cancer. *J Immunother* 46:64-73. [PMID 36637978](https://pubmed.ncbi.nlm.nih.gov/36637978/)
151. Ran X et al. (2025). PARP inhibitor radiosensitization enhances anti-PD-L1 immunotherapy through stabilizing chemokine mRNA in small cell lung cancer. *Nat Commun* 16:2166. [PMID 40038278](https://pubmed.ncbi.nlm.nih.gov/40038278/)
152. Kordbacheh T et al. (2018). Radiotherapy and anti-PD-1/PD-L1 combinations in lung cancer: building better translational research platforms. *Ann Oncol* 29:301-310. [PMID 29309540](https://pubmed.ncbi.nlm.nih.gov/29309540/)
153. Shu Z et al. (2024). PD-L1 deglycosylation promotes its nuclear translocation and accelerates DNA double-strand-break repair in cancer. *Nat Commun* 15:6830. [PMID 39122729](https://pubmed.ncbi.nlm.nih.gov/39122729/)

---

## 15. 정위체부방사선치료(SBRT) 독성 (SBRT Toxicity)

말초 SBRT 의 낮은 MLD 와 급격한 선량 감쇠. 모델이 SBRT 폐렴을 과소예측하는 지점이기도 하다(§한계).

154. Hanania AN et al. (2019). Radiation-Induced Lung Injury: Assessment and Management. *Chest* 156:150-162. [PMID 30998908](https://pubmed.ncbi.nlm.nih.gov/30998908/)
155. Yan M et al. (2023). Stereotactic body radiotherapy for Ultra-Central lung Tumors: A systematic review and Meta-Analysis and International Stereotactic Radiosurgery Society practice guidelines. *Lung Cancer* 182:107281. [PMID 37393758](https://pubmed.ncbi.nlm.nih.gov/37393758/)
156. Chang JY et al. (2021). Stereotactic ablative radiotherapy for operable stage I non-small-cell lung cancer (revised STARS): long-term results of a single-arm, prospective trial with prespecified comparison to surgery. *Lancet Oncol* 22:1448-1457. [PMID 34529930](https://pubmed.ncbi.nlm.nih.gov/34529930/)
157. Lindberg K et al. (2021). The HILUS-Trial-a Prospective Nordic Multicenter Phase 2 Study of Ultracentral Lung Tumors Treated With Stereotactic Body Radiotherapy. *J Thorac Oncol* 16:1200-1210. [PMID 33823286](https://pubmed.ncbi.nlm.nih.gov/33823286/)
158. Blais E et al. (2017). [Lung dose constraints for normo-fractionated radiotherapy and for stereotactic body radiation therapy]. *Cancer Radiother* 21:584-596. [PMID 28886981](https://pubmed.ncbi.nlm.nih.gov/28886981/)
159. Hoffmann L et al. (2022). Thorough design and pre-trial quality assurance (QA) decrease dosimetric impact of delineation and dose planning variability in the STRICTLUNG and STARLUNG trials for stereotactic body radiotherapy (SBRT) of central and ultra-central lung tumours. *Radiother Oncol* 171:53-61. [PMID 35421513](https://pubmed.ncbi.nlm.nih.gov/35421513/)
160. Levy A et al. (2024). Stereotactic Body Radiotherapy for Centrally Located Inoperable Early-Stage NSCLC: EORTC 22113-08113 LungTech Phase II Trial Results. *J Thorac Oncol* 19:1297-1309. [PMID 38788924](https://pubmed.ncbi.nlm.nih.gov/38788924/)

---

## 16. 주요 임상시험과 조사 기법 비교 (Key Trials & Delivery Techniques)

RTOG 0617 · 양성자 vs IMRT · 기능적 폐 회피 · 저분할 요법.

161. Movsas B et al. (2016). Quality of Life Analysis of a Radiation Dose-Escalation Study of Patients With Non-Small-Cell Lung Cancer: A Secondary Analysis of the Radiation Therapy Oncology Group 0617 Randomized Clinical Trial. *JAMA Oncol* 2:359-67. [PMID 26606200](https://pubmed.ncbi.nlm.nih.gov/26606200/)
162. Ma L et al. (2019). A current review of dose-escalated radiotherapy in locally advanced non-small cell lung cancer. *Radiol Oncol* 53:6-14. [PMID 30840594](https://pubmed.ncbi.nlm.nih.gov/30840594/)
163. McNew LK et al. (2017). The relationship between cardiac radiation dose and mediastinal lymph node involvement in stage III non-small cell lung cancer patients. *Adv Radiat Oncol* 2:192-196. [PMID 28740931](https://pubmed.ncbi.nlm.nih.gov/28740931/)
164. Hudson A et al. (2018). Is heterogeneity in stage 3 non-small cell lung cancer obscuring the potential benefits of dose-escalated concurrent chemo-radiotherapy in clinical trials?. *Lung Cancer* 118:139-147. [PMID 29571993](https://pubmed.ncbi.nlm.nih.gov/29571993/)
165. Zou Z et al. (2020). Scanning Beam Proton Therapy versus Photon IMRT for Stage III Lung Cancer: Comparison of Dosimetry, Toxicity, and Outcomes. *Adv Radiat Oncol* 5:434-443. [PMID 32529138](https://pubmed.ncbi.nlm.nih.gov/32529138/)
166. Yu NY et al. (2022). Cardiopulmonary Toxicity Following Intensity-Modulated Proton Therapy (IMPT) Versus Intensity-Modulated Radiation Therapy (IMRT) for Stage III Non-Small Cell Lung Cancer. *Clin Lung Cancer* 23:e526-e535. [PMID 36104272](https://pubmed.ncbi.nlm.nih.gov/36104272/)
167. Vogelius IR et al. (2011). Estimated radiation pneumonitis risk after photon versus proton therapy alone or combined with chemotherapy for lung cancer. *Acta Oncol* 50:772-6. [PMID 21767173](https://pubmed.ncbi.nlm.nih.gov/21767173/)
168. Ireland RH et al. (2016). Functional Image-guided Radiotherapy Planning for Normal Lung Avoidance. *Clin Oncol (R Coll Radiol)* 28:695-707. [PMID 27637724](https://pubmed.ncbi.nlm.nih.gov/27637724/)
169. Eslick EM et al. (2019). SPECT V/Q in Lung Cancer Radiotherapy Planning. *Semin Nucl Med* 49:31-36. [PMID 30545514](https://pubmed.ncbi.nlm.nih.gov/30545514/)
170. Bedford JL et al. (2023). Functional lung avoidance in radiotherapy using optimisation of biologically effective dose with non-coplanar beam orientations. *Phys Imaging Radiat Oncol* 28:100518. [PMID 38077270](https://pubmed.ncbi.nlm.nih.gov/38077270/)
171. Swaminath A et al. (2024). Stereotactic vs Hypofractionated Radiotherapy for Inoperable Stage I Non-Small Cell Lung Cancer: The LUSTRE Phase 3 Randomized Clinical Trial. *JAMA Oncol* 10:1571-1575. [PMID 39298144](https://pubmed.ncbi.nlm.nih.gov/39298144/)
172. Iyengar P et al. (2021). Accelerated Hypofractionated Image-Guided vs Conventional Radiotherapy for Patients With Stage II/III Non-Small Cell Lung Cancer and Poor Performance Status: A Randomized Clinical Trial. *JAMA Oncol* 7:1497-1505. [PMID 34383006](https://pubmed.ncbi.nlm.nih.gov/34383006/)
173. Iyengar P et al. (2018). Consolidative Radiotherapy for Limited Metastatic Non-Small-Cell Lung Cancer: A Phase 2 Randomized Clinical Trial. *JAMA Oncol* 4:e173501. [PMID 28973074](https://pubmed.ncbi.nlm.nih.gov/28973074/)

---

## 17. 환자 위험인자와 예측모델 (Patient Risk Factors & Prediction Models)

기존 간질성 폐질환·COPD·유전다형성. 모델에서는 기저 예비능 RF = √(DLREF/DLCO₀) 로 들어가 같은 손상을 다른 등급으로 만든다.

174. Yamaguchi T et al. (2021). Pre-existing interstitial lung disease is associated with onset of nivolumab-induced pneumonitis in patients with solid tumors: a retrospective analysis. *BMC Cancer* 21:924. [PMID 34399710](https://pubmed.ncbi.nlm.nih.gov/34399710/)
175. Lu H et al. (2026). Integrated machine learning risk model for predicting radiation pneumonitis in lung cancer patients with interstitial lung disease. *Ann Med* 58:2634444. [PMID 41755502](https://pubmed.ncbi.nlm.nih.gov/41755502/)
176. Walls GM et al. (2023). Clinicoradiological outcomes after radical radiotherapy for lung cancer in patients with interstitial lung disease. *BJR Open* 5:20220049. [PMID 37389005](https://pubmed.ncbi.nlm.nih.gov/37389005/)
177. Lu H et al. (2026). Time-of-day radiotherapy alters radiation pneumonitis risk but not survival in lung cancer patients with interstitial lung disease. *Ther Adv Med Oncol* 18:17588359261446805. [PMID 42147351](https://pubmed.ncbi.nlm.nih.gov/42147351/)
178. Ishijima M et al. (2015). Patients with severe emphysema have a low risk of radiation pneumonitis following stereotactic body radiotherapy. *Br J Radiol* 88:20140596. [PMID 25490255](https://pubmed.ncbi.nlm.nih.gov/25490255/)
179. Zhou Z et al. (2017). Pulmonary emphysema is a risk factor for radiation pneumonitis in NSCLC patients with squamous cell carcinoma after thoracic radiation therapy. *Sci Rep* 7:2748. [PMID 28584268](https://pubmed.ncbi.nlm.nih.gov/28584268/)
180. Uchida Y et al. (2017). Exclusion of emphysematous lung from dose-volume estimates of risk improves prediction of radiation pneumonitis. *Radiat Oncol* 12:160. [PMID 28969651](https://pubmed.ncbi.nlm.nih.gov/28969651/)
181. Wang Y et al. (2015). Effect of transforming growth factor-β1 869C/T polymorphism and radiation pneumonitis. *Int J Clin Exp Pathol* 8:2835-9. [PMID 26045792](https://pubmed.ncbi.nlm.nih.gov/26045792/)
182. Alsbeih G et al. (2010). Association between normal tissue complications after radiotherapy and polymorphic variations in TGFB1 and XRCC1 genes. *Radiat Res* 173:505-11. [PMID 20334523](https://pubmed.ncbi.nlm.nih.gov/20334523/)
183. Fernet M et al. (2004). Genetic biomarkers of therapeutic radiation sensitivity. *DNA Repair (Amst)* 3:1237-43. [PMID 15279812](https://pubmed.ncbi.nlm.nih.gov/15279812/)
184. Yiu WS et al. (2024). DNA Repair Genetics and the Risk of Radiation Pneumonitis in Patients With Lung Cancer: A Systematic Review and Meta-analysis. *Clin Oncol (R Coll Radiol)* 36:e182-e196. [PMID 38653664](https://pubmed.ncbi.nlm.nih.gov/38653664/)
185. Andreassen CN et al. (2006). Risk of radiation-induced subcutaneous fibrosis in relation to single nucleotide polymorphisms in TGFB1, SOD2, XRCC1, XRCC3, APEX and ATM--a study based on DNA from formalin fixed paraffin embedded tissue samples. *Int J Radiat Biol* 82:577-86. [PMID 16966185](https://pubmed.ncbi.nlm.nih.gov/16966185/)
186. Kawahara D et al. (2025). Prediction of radiation pneumonitis after CRT in patients with advanced NSCLC using multi-region radiomics and attention-based ensemble learning. *Med Phys* 52:e70140. [PMID 41261069](https://pubmed.ncbi.nlm.nih.gov/41261069/)
187. Yu H et al. (2019). Machine Learning to Build and Validate a Model for Radiation Pneumonitis Prediction in Patients with Non-Small Cell Lung Cancer. *Clin Cancer Res* 25:4343-4350. [PMID 30992302](https://pubmed.ncbi.nlm.nih.gov/30992302/)
188. Ai Y et al. (2025). Integrating deep learning and multi-omics features in radiation pneumonitis prediction for lung cancer patients using PET/CT. *BMC Med Imaging* 25:426. [PMID 41146084](https://pubmed.ncbi.nlm.nih.gov/41146084/)

---

## 18. 폐기능 변화와 영상 정량화 (Pulmonary Function & Imaging Quantification)

흉부 방사선치료 후 DLCO 감소, CT 섬유화 정량화.

189. Niezink AGH et al. (2017). Pulmonary Function Changes After Radiotherapy for Lung or Esophageal Cancer: A Systematic Review Focusing on Dose-Volume Parameters. *Oncologist* 22:1257-1264. [PMID 28550029](https://pubmed.ncbi.nlm.nih.gov/28550029/)
190. Karlsen J et al. (2024). Pulmonary Function and Lung Fibrosis up to 12 Years After Breast Cancer Radiotherapy. *Int J Radiat Oncol Biol Phys* 118:1066-1077. [PMID 38099884](https://pubmed.ncbi.nlm.nih.gov/38099884/)
191. Schröder C et al. (2017). Changes in pulmonary function and influencing factors after high-dose intrathoracic radio(chemo)therapy. *Strahlenther Onkol* 193:125-131. [PMID 27783103](https://pubmed.ncbi.nlm.nih.gov/27783103/)
192. Henderson M et al. (2008). Baseline pulmonary function as a predictor for survival and decline in pulmonary function over time in patients undergoing stereotactic body radiotherapy for the treatment of stage I non-small-cell lung cancer. *Int J Radiat Oncol Biol Phys* 72:404-9. [PMID 18394819](https://pubmed.ncbi.nlm.nih.gov/18394819/)
193. Mahmutovic Persson I et al. (2020). Imaging Biomarkers and Pathobiological Profiling in a Rat Model of Drug-Induced Interstitial Lung Disease Induced by Bleomycin. *Front Physiol* 11:584. [PMID 32636756](https://pubmed.ncbi.nlm.nih.gov/32636756/)
194. Zaher A et al. (2024). Exploratory Analysis of Image-Guided Ionizing Radiation Delivery to Induce Long-Term Iron Accumulation and Ferritin Expression in a Lung Injury Model: Preliminary Results. *Bioengineering (Basel)* 11. [PMID 38391668](https://pubmed.ncbi.nlm.nih.gov/38391668/)
195. Chen L et al. (2019). The impact of right ventricular function on prognosis in patients with stage III non-small cell lung cancer after concurrent chemoradiotherapy. *Int J Cardiovasc Imaging* 35:1009-1017. [PMID 30941563](https://pubmed.ncbi.nlm.nih.gov/30941563/)
196. Torre-Bouscoulet L et al. (2018). Longitudinal Evaluation of Lung Function in Patients With Advanced Non-Small Cell Lung Cancer Treated With Concurrent Chemoradiation Therapy. *Int J Radiat Oncol Biol Phys* 101:910-918. [PMID 29976503](https://pubmed.ncbi.nlm.nih.gov/29976503/)
197. Wada S et al. (2026). Four-Dimensional Computed Tomography of Respiratory Function Changes Post-Radiotherapy for Lung Cancer. *Cureus* 18:e103971. [PMID 41727802](https://pubmed.ncbi.nlm.nih.gov/41727802/)

---

## 19. QSP·정량 모델링 방법론 (QSP & Quantitative Modelling Methodology)

mrgsolve 및 모델 기반 신약개발(MIDD) 방법론.

198. Klionsky DJ et al. (2021). Guidelines for the use and interpretation of assays for monitoring autophagy (4th edition)(1). *Autophagy* 17:1-382. [PMID 33634751](https://pubmed.ncbi.nlm.nih.gov/33634751/)
199. Kosinsky Y et al. (2018). Radiation and PD-(L)1 treatment combinations: immune response and dose optimization via a predictive systems model. *J Immunother Cancer* 6:17. [PMID 29486799](https://pubmed.ncbi.nlm.nih.gov/29486799/)
200. Zhang LW et al. (2025). IL-23 Receptor Agonism by Mulberroside C Activates the RASGRP1/RAS/ERK Pathway Contributing to Leukopenia Treatments. *Phytother Res* 39:3578-3600. [PMID 40619172](https://pubmed.ncbi.nlm.nih.gov/40619172/)
201. Elmokadem A et al. (2019). Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol* 8:883-893. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
202. Jiang J et al. (2025). Naringenin inhibits ferroptosis to reduce radiation-induced lung injury: insights from network Pharmacology and molecular docking. *Pharm Biol* 63:1-10. [PMID 39969099](https://pubmed.ncbi.nlm.nih.gov/39969099/)
203. Zhi Y et al. (2025). Novel NLRP3 inhibitors mitigate acute radiation-induced lung injury by suppressing pyroptosis in alveolar epithelial cells. *Toxicol Appl Pharmacol* 502:117458. [PMID 40618790](https://pubmed.ncbi.nlm.nih.gov/40618790/)
204. Li J et al. (2026). Mai-Men-Dong decoction alleviates radiation-induced lung injury by regulating ferroptosis through modulation of the STAT6/p53/SLC7A11 axis. *J Ethnopharmacol* 360:121072. [PMID 41448369](https://pubmed.ncbi.nlm.nih.gov/41448369/)
205. GBD 2021 Stroke Risk Factor Collaborators (2024). Global, regional, and national burden of stroke and its risk factors, 1990-2021: a systematic analysis for the Global Burden of Disease Study 2021. *Lancet Neurol* 23:973-1003. [PMID 39304265](https://pubmed.ncbi.nlm.nih.gov/39304265/)
206. Bai JPF et al. (2024). Creating a Roadmap to Quantitative Systems Pharmacology-Informed Rare Disease Drug Development: A Workshop Report. *Clin Pharmacol Ther* 115:201-205. [PMID 37984065](https://pubmed.ncbi.nlm.nih.gov/37984065/)
207. Zhu AZX et al. (2022). Applications of Quantitative System Pharmacology Modeling to Model-Informed Drug Development. *Methods Mol Biol* 2486:71-86. [PMID 35437719](https://pubmed.ncbi.nlm.nih.gov/35437719/)
208. Lesko LJ (2021). Perspective on model-informed drug development. *CPT Pharmacometrics Syst Pharmacol* 10:1127-1129. [PMID 34404115](https://pubmed.ncbi.nlm.nih.gov/34404115/)
209. Cogno N et al. (2022). An Agent-Based Model of Radiation-Induced Lung Fibrosis. *Int J Mol Sci* 23. [PMID 36430398](https://pubmed.ncbi.nlm.nih.gov/36430398/)
210. Cogno N et al. (2024). Mechanistic model of radiotherapy-induced lung fibrosis using coupled 3D agent-based and Monte Carlo simulations. *Commun Med (Lond)* 4:16. [PMID 38336802](https://pubmed.ncbi.nlm.nih.gov/38336802/)
211. Wang YC et al. (2024). Hyaluronic acid-based injectable formulation developed to mitigate metastasis and radiation-induced skin fibrosis in breast cancer treatment. *Carbohydr Polym* 336:122136. [PMID 38670762](https://pubmed.ncbi.nlm.nih.gov/38670762/)

---

## 20. 그 밖의 기전·현상 (Other Mechanisms & Phenomena)

노화·NOX4·대식세포 분극·림프구감소·재조사·FLASH·기질화 폐렴·야외(out-of-field) 손상 등.

212. Dörr W et al. (2001). Consequential late effects in normal tissues. *Radiother Oncol* 61:223-31. [PMID 11730991](https://pubmed.ncbi.nlm.nih.gov/11730991/)
213. Dörr W (2015). Radiobiology of tissue reactions. *Ann ICRP* 44:58-68. [PMID 25816259](https://pubmed.ncbi.nlm.nih.gov/25816259/)
214. Palmer JD et al. (2021). Late effects of radiation therapy in pediatric patients and survivorship. *Pediatr Blood Cancer* 68 Suppl 2:e28349. [PMID 33818893](https://pubmed.ncbi.nlm.nih.gov/33818893/)
215. Cheng B et al. (2022). Anti-PD-L1/TGF-βR fusion protein (SHR-1701) overcomes disrupted lymphocyte recovery-induced resistance to PD-1/PD-L1 inhibitors in lung cancer. *Cancer Commun (Lond)* 42:17-36. [PMID 34981670](https://pubmed.ncbi.nlm.nih.gov/34981670/)
216. Prades-Sagarra È et al. (2025). L19-IL2 reverts radiation-induced lymphopenia in a mouse model of lung cancer. *Radiother Oncol* 208:110908. [PMID 40288691](https://pubmed.ncbi.nlm.nih.gov/40288691/)
217. van Rossum PSN et al. (2023). Severe radiation-induced lymphopenia during concurrent chemoradiotherapy for stage III non-small cell lung cancer: external validation of two prediction models. *Front Oncol* 13:1278723. [PMID 38023221](https://pubmed.ncbi.nlm.nih.gov/38023221/)
218. He Y et al. (2019). Cellular senescence and radiation-induced pulmonary fibrosis. *Transl Res* 209:14-21. [PMID 30981698](https://pubmed.ncbi.nlm.nih.gov/30981698/)
219. Su L et al. (2021). Potential role of senescent macrophages in radiation-induced pulmonary fibrosis. *Cell Death Dis* 12:527. [PMID 34023858](https://pubmed.ncbi.nlm.nih.gov/34023858/)
220. Zhou S et al. (2022). Alveolar type 2 epithelial cell senescence and radiation-induced pulmonary fibrosis. *Front Cell Dev Biol* 10:999600. [PMID 36407111](https://pubmed.ncbi.nlm.nih.gov/36407111/)
221. Zhang Y et al. (2022). SIRT1 prevents cigarette smoking-induced lung fibroblasts activation by regulating mitochondrial oxidative stress and lipid metabolism. *J Transl Med* 20:222. [PMID 35568871](https://pubmed.ncbi.nlm.nih.gov/35568871/)
222. Zheng M et al. (2024). Traditional Chinese medicine inspired dual-drugs loaded inhalable nano-therapeutics alleviated idiopathic pulmonary fibrosis by targeting early inflammation and late fibrosis. *J Nanobiotechnology* 22:14. [PMID 38166847](https://pubmed.ncbi.nlm.nih.gov/38166847/)
223. Fang L et al. (2021). Osthole Attenuates Bleomycin-Induced Pulmonary Fibrosis by Modulating NADPH Oxidase 4-Derived Oxidative Stress in Mice. *Oxid Med Cell Longev* 2021:3309944. [PMID 34527170](https://pubmed.ncbi.nlm.nih.gov/34527170/)
224. Ni J et al. (2023). STING signaling activation modulates macrophage polarization via CCL2 in radiation-induced lung injury. *J Transl Med* 21:590. [PMID 37667317](https://pubmed.ncbi.nlm.nih.gov/37667317/)
225. Favaudon V et al. (2014). Ultrahigh dose-rate FLASH irradiation increases the differential response between normal and tumor tissue in mice. *Sci Transl Med* 6:245ra93. [PMID 25031268](https://pubmed.ncbi.nlm.nih.gov/25031268/)
226. Lee SE et al. (2026). Localized normal tissue-sparing effects of proton FLASH radiotherapy in a preclinical lung irradiation model. *Br J Radiol* 99:459-467. [PMID 41564308](https://pubmed.ncbi.nlm.nih.gov/41564308/)
227. Dubail M et al. (2025). Sparing effects of FLASH irradiation in patient-derived lung tissue. *Radiother Oncol* 212:111126. [PMID 40921334](https://pubmed.ncbi.nlm.nih.gov/40921334/)
228. Hong JH et al. (2019). High-Dose Thoracic Re-irradiation of Lung Cancer Using Highly Conformal Radiotherapy Is Effective with Acceptable Toxicity. *Cancer Res Treat* 51:1156-1166. [PMID 30514067](https://pubmed.ncbi.nlm.nih.gov/30514067/)
229. Rulach R et al. (2021). An International Expert Survey on the Indications and Practice of Radical Thoracic Reirradiation for Non-Small Cell Lung Cancer. *Adv Radiat Oncol* 6:100653. [PMID 33851065](https://pubmed.ncbi.nlm.nih.gov/33851065/)
230. Rulach R et al. (2026). Cumulative oesophageal dose and risk of high-grade toxicity in thoracic re-irradiation: a dose/toxicity analysis. *Clin Transl Radiat Oncol* 57:101108. [PMID 41583535](https://pubmed.ncbi.nlm.nih.gov/41583535/)
231. Kuriakose J et al. (2023). Osimertinib-induced radiation recall pneumonitis. *BMJ Support Palliat Care*. [PMID 37116943](https://pubmed.ncbi.nlm.nih.gov/37116943/)
232. Kalisz KR et al. (2019). Immune Checkpoint Inhibitor Therapy-related Pneumonitis: Patterns and Management. *Radiographics* 39:1923-1937. [PMID 31584861](https://pubmed.ncbi.nlm.nih.gov/31584861/)
233. Grassi F et al. (2023). Radiation Recall Pneumonitis: The Open Challenge in Differential Diagnosis of Pneumonia Induced by Oncological Treatments. *J Clin Med* 12. [PMID 36835977](https://pubmed.ncbi.nlm.nih.gov/36835977/)
234. Epler GR et al. (2020). Post-Breast Cancer Radiotherapy Bronchiolitis Obliterans Organizing Pneumonia. *Respir Care* 65:686-692. [PMID 31892515](https://pubmed.ncbi.nlm.nih.gov/31892515/)
235. Ailloud A et al. (2024). [Bronchiolitis obliterans organizing pneumonia after radiotherapy: A systematic review and case report]. *Cancer Radiother* 28:707-718. [PMID 39581827](https://pubmed.ncbi.nlm.nih.gov/39581827/)
236. Caroprese M et al. (2024). Bronchiolitis Obliterans Organizing Pneumonia After Breast Radiation Therapy. *Pract Radiat Oncol* 14:e443-e448. [PMID 39032596](https://pubmed.ncbi.nlm.nih.gov/39032596/)

---

## 인용 대상이 없는 구조적 선택 (Structural choices with no citation)

아래 항목은 문헌에서 직접 가져온 것이 아니라 이 모델의 **구조적 결정**이며,
근거가 아니라 가정으로 읽어야 합니다.

- **6개 DVH 구간으로의 이산화.** 실제 DVH는 연속이며, 6개 구간은 체적효과를
  물을 수 있는 최소 해상도로 선택했습니다. 구간 수를 늘리면 유도된 n 값이
  달라질 수 있습니다.
- **BED 전달률의 평활화.** 분할 조사를 하루 단위 볼루스가 아니라
  `RB_b = BED_b / T_course` 의 연속 전달률로 표현했습니다. 총 BED와 분할
  의존성은 정확히 보존되지만, 하루 이내의 동태(아분치사 손상 회복)는
  표현되지 않습니다. 하류의 모든 시간상수가 ≥ 일 단위이므로 채택했습니다.
- **PNI50 · PNISL.** 기전으로 식별되지 않는 두 개의 척도 상수이며 QUANTEC
  용량-반응에 적합시켰습니다(`rili_calibration.py`).
- **CTCAE 등급 문턱.** 연속 지수를 등급으로 바꾸는 절단점은 관측이 아니라
  모델의 산출 분포에 맞춘 규약입니다.
- **경구 약물의 연속 투여율 근사.** 프레드니솔론·피르페니돈·닌테다닙·
  리시노프릴은 mg/day 연속 입력으로 표현했습니다. 더발루맙만 실제 q28d
  볼루스입니다.
- **아미포스틴의 대수적 처리.** 반감기 8분을 ODE로 적분하는 대신
  보호계수를 투여-조사 간격의 함수로 대수적으로 계산했습니다.

