# 관상동맥 미세혈관 기능장애 (CMD / ANOCA–INOCA) — 참고문헌

> Coronary Microvascular Dysfunction · Angina/Ischaemia with No Obstructive Coronary Arteries
> QSP 모델(`cmd_qsp_model.dot`, `cmd_mrgsolve_model.R`, `cmd_reference_model.py`)의
> 모든 파라미터·구조·보정 목표의 근거.

**인용 원칙.** 이 목록의 모든 PMID는 이 세션에서 PubMed E-utilities로 직접
조회하여 제목·저널·연도를 확인한 것입니다. 모델의 정량적 주장(§2, §3, §11)에
쓰인 문헌은 초록 전문을 조회해 수치를 확인했으며, 해당 항목에는 **모델에서
사용한 수치**를 함께 적었습니다. 확인하지 못한 문헌은 넣지 않았습니다.

---

## 1. 정의 · 명명법 · 진료지침 (Definitions, nomenclature, guidelines)

1. Ong P, et al. **International standardization of diagnostic criteria for microvascular angina.** Int J Cardiol. 2018. [PMID 29031990](https://pubmed.ncbi.nlm.nih.gov/29031990/) — COVADIS 미세혈관 협심증 진단기준. 모델의 endotype 분류 축(CFR·미세혈관저항·ACh 반응)이 이 기준을 따릅니다.
2. Beltrame JF, et al. **International standardization of diagnostic criteria for vasospastic angina.** Eur Heart J. 2017. [PMID 26245334](https://pubmed.ncbi.nlm.nih.gov/26245334/) — COVADIS 경련성 협심증 기준.
3. Vrints C, et al. **2024 ESC Guidelines for the management of chronic coronary syndromes.** Eur Heart J. 2024. [PMID 39210710](https://pubmed.ncbi.nlm.nih.gov/39210710/)
4. Vrints C, et al. **[2024 ESC 만성 관상동맥증후군 지침, 이탈리아어판].** G Ital Cardiol. 2024. [PMID 39611224](https://pubmed.ncbi.nlm.nih.gov/39611224/)
5. Gulati M, et al. **2021 AHA/ACC/ASE/CHEST/SAEM/SCCT/SCMR Guideline for the Evaluation and Diagnosis of Chest Pain.** Circulation. 2021. [PMID 34709879](https://pubmed.ncbi.nlm.nih.gov/34709879/)
6. Writing Committee Members. **2021 AHA/ACC Chest Pain Guideline (JACC).** J Am Coll Cardiol. 2021. [PMID 34756653](https://pubmed.ncbi.nlm.nih.gov/34756653/)
7. Gulati M, et al. **2021 Chest Pain Guideline: Executive Summary.** Circulation. 2021. [PMID 34709928](https://pubmed.ncbi.nlm.nih.gov/34709928/)
8. Boden WE, et al. **Myocardial Ischemic Syndromes: A New Nomenclature to Harmonize Evolving International Clinical Practice Guidelines.** Circulation. 2024. [PMID 39210827](https://pubmed.ncbi.nlm.nih.gov/39210827/)
9. Kaur G, et al. **Chest Pain in Women: Considerations From the 2021 AHA/ACC Chest Pain Guideline.** Curr Probl Cardiol. 2023. [PMID 36921653](https://pubmed.ncbi.nlm.nih.gov/36921653/)
10. Pepine CJ, et al. **ANOCA/INOCA/MINOCA: Open artery ischemia.** Am Heart J Plus. 2023. [PMID 37064505](https://pubmed.ncbi.nlm.nih.gov/37064505/)
11. Ashokprabhu ND, et al. **INOCA/ANOCA: Mechanisms and novel treatments.** Am Heart J Plus. 2023. [PMID 37377840](https://pubmed.ncbi.nlm.nih.gov/37377840/)
12. Rinaldi R, et al. **Management of angina pectoris.** Trends Cardiovasc Med. 2025. [PMID 40086653](https://pubmed.ncbi.nlm.nih.gov/40086653/)
13. Burgess S, et al. **Challenges in diagnosing coronary microvascular dysfunction and coronary vasospasm.** Cardiovasc Revasc Med. 2025. [PMID 40312200](https://pubmed.ncbi.nlm.nih.gov/40312200/)
14. Scarica V, et al. **Coronary microvascular dysfunction: pathophysiology, diagnosis, and therapeutic strategies across cardiovascular diseases.** EXCLI J. 2025. [PMID 40376434](https://pubmed.ncbi.nlm.nih.gov/40376434/)
15. Ya'Qoub L, et al. **Non-obstructive Plaque and Treatment of INOCA: More to Be Learned.** Curr Atheroscler Rep. 2022. [PMID 35781776](https://pubmed.ncbi.nlm.nih.gov/35781776/)
16. Vancheri F, et al. **Coronary Microvascular Dysfunction.** J Clin Med. 2020. [PMID 32899944](https://pubmed.ncbi.nlm.nih.gov/32899944/)
17. Smilowitz NR, et al. **Coronary Microvascular Disease in Contemporary Clinical Practice.** Circ Cardiovasc Interv. 2023. [PMID 37259860](https://pubmed.ncbi.nlm.nih.gov/37259860/)
18. Parwani P, et al. **Contemporary Diagnosis and Management of Patients with MINOCA.** Curr Cardiol Rep. 2023. [PMID 37067753](https://pubmed.ncbi.nlm.nih.gov/37067753/)
19. Boden WE, et al. **Evolving Management Paradigm for Stable Ischemic Heart Disease Patients.** J Am Coll Cardiol. 2023. [PMID 36725179](https://pubmed.ncbi.nlm.nih.gov/36725179/)

---

## 2. 이 모델의 척추 — Endotype과 침습적 생리 (The model's backbone)

이 절의 첫 문헌이 모델 전체의 구조적 출발점입니다. **CFR은 비(ratio)이며,
같은 낮은 값이 분모(안정 유량 과다)에서 오기도 하고 분자(최대충혈 유량 한계)에서
오기도 합니다.**

20. **Rahman H, et al. Coronary Microvascular Dysfunction Is Associated With Myocardial Ischemia and Abnormal Coronary Perfusion During Exercise. Circulation. 2019;140:1805–1816.** [PMID 31707835](https://pubmed.ncbi.nlm.nih.gov/31707835/)
    **모델에서 사용한 수치(보정 목표 T1–T5):** ANOCA 85명(78% 여성, 57±10세). CFR<2.5를 MVD로 정의 → 45명(53%). MVD 중 62%는 **기능적**(최대충혈 미세혈관저항 <2.5 mmHg/cm/s), 38%는 **구조적**(≥2.5). 안정 미세혈관저항 기능적 4.2±1.0, 구조적 6.9±1.7, 대조 7.3±2.2 mmHg/(cm/s). MVD군 82%에서 유발성 허혈(최대충혈 심내막하/심표면 관류비 <1.0) 대 대조군 22%. 전체 심근 관류예비능 2.01±0.41 대 2.68±0.49. 관류 효율은 대조군에서 운동 시 59±11→65±14%로 상승하나 MVD에서는 61±12→44±10%로 하락. 운동 시 수축기압 구조적 188±25 대 기능적 161±27 대 대조 156±30 mmHg(모델 파라미터 `KSBP_STRUCT/FUNC/CTRL`의 근거). **두 endotype의 부하 시 심근 관류와 운동 관류효율은 서로 비슷했다** — 모델이 재현해야 할, 그리고 §4의 "reserve exhaustion workload"로 설명하는 긴장(tension).
21. Sinha A, et al. **Rethinking False Positive Exercise Electrocardiographic Stress Tests by Assessing Coronary Microvascular Function.** J Am Coll Cardiol. 2024. [PMID 38199706](https://pubmed.ncbi.nlm.nih.gov/38199706/) — ANOCA 102명. 운동 ECG에서 허혈이 나타나면 CMD에 대해 **특이도 100%**, ACh 유량예비능이 운동 허혈의 최강 예측인자. 모델의 "운동 부하검사는 위양성이 아니라 미세혈관 기질을 보고 있다"는 전제.
22. Rahman H, et al. **Optimal Use of Vasodilators for Diagnosis of Microvascular Angina in the Cardiac Catheterization Laboratory.** Circ Cardiovasc Interv. 2020. [PMID 32519879](https://pubmed.ncbi.nlm.nih.gov/32519879/)
23. Ford TJ, et al. **Assessment of Vascular Dysfunction in Patients Without Obstructive Coronary Artery Disease: Why, How, and When.** JACC Cardiovasc Interv. 2020. [PMID 32819476](https://pubmed.ncbi.nlm.nih.gov/32819476/)
24. Suda A, et al. **Coronary Functional Abnormalities in Patients With Angina and Nonobstructive Coronary Artery Disease.** J Am Coll Cardiol. 2019. [PMID 31699275](https://pubmed.ncbi.nlm.nih.gov/31699275/)
25. Ang DTY, et al. **Phenotype-based management of coronary microvascular dysfunction.** J Nucl Cardiol. 2022. [PMID 35672569](https://pubmed.ncbi.nlm.nih.gov/35672569/)
26. Lanza GA, Crea F. **Primary coronary microvascular dysfunction: clinical presentation, pathophysiology, and management.** Circulation. 2010. [PMID 20516386](https://pubmed.ncbi.nlm.nih.gov/20516386/)
27. Crea F, et al. **Coronary microvascular dysfunction: an update.** Eur Heart J. 2014. [PMID 24366916](https://pubmed.ncbi.nlm.nih.gov/24366916/)
28. Crea F, et al. **Pathophysiology of Coronary Microvascular Dysfunction.** Circ J. 2022. [PMID 34759123](https://pubmed.ncbi.nlm.nih.gov/34759123/)
29. Del Buono MG, et al. **Coronary Microvascular Dysfunction Across the Spectrum of Cardiovascular Diseases: JACC State-of-the-Art Review.** J Am Coll Cardiol. 2021. [PMID 34556322](https://pubmed.ncbi.nlm.nih.gov/34556322/)
30. Bairey Merz CN, et al. **Treatment of coronary microvascular dysfunction.** Cardiovasc Res. 2020. [PMID 32087007](https://pubmed.ncbi.nlm.nih.gov/32087007/)
31. Nogami K, et al. **Chest pain patterns and coronary microvascular function in non-obstructive coronary artery disease.** EuroIntervention. 2025. [PMID 40887985](https://pubmed.ncbi.nlm.nih.gov/40887985/)
32. Beck S, et al. **Invasive Diagnosis of Coronary Functional Disorders Causing Angina Pectoris.** Eur Cardiol. 2021. [PMID 34276812](https://pubmed.ncbi.nlm.nih.gov/34276812/)
33. Chalikias G, et al. **Slow Coronary Flow: Pathophysiology, Clinical Implications, and Therapeutic Management.** Angiology. 2021. [PMID 33779300](https://pubmed.ncbi.nlm.nih.gov/33779300/)
34. Lee SH, et al. **Clinical Relevance of Ischemia with Nonobstructive Coronary Arteries According to Coronary Microvascular Dysfunction.** J Am Heart Assoc. 2022. [PMID 35475358](https://pubmed.ncbi.nlm.nih.gov/35475358/)

---

## 3. 지표 — CFR · IMR · MRR · 절대혈류 (Indices)

모델은 **비(ratio)와 저항(resistance)을 구별**합니다. 치료 중에도 해석이
유지되는 지표는 최대충혈 저항이며, 이 절의 문헌들이 그 지표들의 정의와 한계를
제공합니다.

35. Fearon WF, et al. **Novel index for invasively assessing the coronary microcirculation.** Circulation. 2003. [PMID 12821539](https://pubmed.ncbi.nlm.nih.gov/12821539/) — IMR의 원 정의(Pd × Tmn). 모델의 `IMR = MR_hyp × 8.35` 스케일링은 정상 territory가 18 U를 읽도록 맞춘 것.
36. Martínez GJ, et al. **The index of microcirculatory resistance in the physiologic assessment of the coronary microcirculation.** Coron Artery Dis. 2015. [PMID 26247265](https://pubmed.ncbi.nlm.nih.gov/26247265/)
37. Yong AS, et al. **Calculation of the index of microcirculatory resistance without coronary wedge pressure measurement in the presence of epicardial stenosis.** JACC Cardiovasc Interv. 2013. [PMID 23347861](https://pubmed.ncbi.nlm.nih.gov/23347861/)
38. De Bruyne B, et al. **Coronary thermodilution to assess flow reserve: experimental validation.** Circulation. 2001. [PMID 11673336](https://pubmed.ncbi.nlm.nih.gov/11673336/)
39. Collet C, et al. **A Systematic Approach to the Evaluation of the Coronary Microcirculation Using Bolus Thermodilution: CATH CMD.** J Soc Cardiovasc Angiogr Interv. 2024. [PMID 39131992](https://pubmed.ncbi.nlm.nih.gov/39131992/)
40. Gutiérrez-Barrios A, et al. **Continuous Thermodilution Method to Assess Coronary Flow Reserve.** Am J Cardiol. 2021. [PMID 33220317](https://pubmed.ncbi.nlm.nih.gov/33220317/)
41. Candreva A, et al. **Automation of intracoronary continuous thermodilution for absolute coronary flow and microvascular resistance measurements.** Catheter Cardiovasc Interv. 2022. [PMID 35723684](https://pubmed.ncbi.nlm.nih.gov/35723684/)
42. Belmonte M, et al. **Measuring Absolute Coronary Flow and Microvascular Resistance by Thermodilution: JACC Review Topic of the Week.** J Am Coll Cardiol. 2024. [PMID 38325996](https://pubmed.ncbi.nlm.nih.gov/38325996/)
43. Pijls NHJ, et al. **Absolute coronary blood flow measurement and the principle of microvascular resistance reserve.** Cardiovasc Interv Ther. 2026. [PMID 41432885](https://pubmed.ncbi.nlm.nih.gov/41432885/) — MRR의 원리. 모델은 MRR이 **약물의 혈압 강하를 보정하지 못한다**고 예측합니다(안정·최대충혈 양쪽에서 같은 압력이 내려가므로 보정항이 상쇄되어 MRR = CFR/FFR로 항등). MRR이 값을 발휘하는 대상은 심표면 협착이며, 그 두 문제는 자주 혼동됩니다.
44. Sakai K, et al. **Impact of vessel volume on thermodilution measurements in patients with coronary microvascular dysfunction.** Catheter Cardiovasc Interv. 2024. [PMID 38566527](https://pubmed.ncbi.nlm.nih.gov/38566527/)
45. Gallinoro E, et al. **Microvascular Dysfunction in Patients With Type II Diabetes Mellitus: Invasive Assessment of Absolute Coronary Blood Flow and Microvascular Resistance Reserve.** Front Cardiovasc Med. 2021. [PMID 34738020](https://pubmed.ncbi.nlm.nih.gov/34738020/)
46. Paolisso P, et al. **Absolute coronary flow and microvascular resistance reserve in patients with severe aortic stenosis.** Heart. 2022. [PMID 35977812](https://pubmed.ncbi.nlm.nih.gov/35977812/)
47. Johnson NP, et al. **Invasive FFR and Noninvasive CFR in the Evaluation of Ischemia: What Is the Future?** J Am Coll Cardiol. 2016. [PMID 27282899](https://pubmed.ncbi.nlm.nih.gov/27282899/)
48. Simova I. **Coronary Flow Velocity Reserve Assessment with Transthoracic Doppler Echocardiography.** Eur Cardiol. 2015. [PMID 30310417](https://pubmed.ncbi.nlm.nih.gov/30310417/)
49. Olsen RH, et al. **Coronary flow velocity reserve by echocardiography: feasibility, reproducibility and agreement with PET.** Cardiovasc Ultrasound. 2016. [PMID 27267255](https://pubmed.ncbi.nlm.nih.gov/27267255/)
50. Loftspring E, et al. **Angiography-Derived Versus Coronary Guidewire-Derived Index of Microcirculatory Resistance in Patients With INOCA.** J Soc Cardiovasc Angiogr Interv. 2025. [PMID 41324045](https://pubmed.ncbi.nlm.nih.gov/41324045/)
51. Gao B, et al. **Quantitative Flow Ratio-Derived Index of Microcirculatory Resistance as a Novel Tool to Identify Microcirculatory Function in Patients with INOCA.** Cardiology. 2024. [PMID 37839404](https://pubmed.ncbi.nlm.nih.gov/37839404/)
52. Wang S, et al. **Myocardial Blood Flow Quantification Using Stress Cardiac Magnetic Resonance Improves Detection of Coronary Artery Disease.** JACC Cardiovasc Imaging. 2024. [PMID 39297850](https://pubmed.ncbi.nlm.nih.gov/39297850/)
53. Rasmussen LD, et al. **Myocardial Blood Flow by Magnetic Resonance in Patients With Suspected Coronary Stenosis: Comparison to PET and Invasive Physiology.** Circ Cardiovasc Imaging. 2024. [PMID 38889213](https://pubmed.ncbi.nlm.nih.gov/38889213/)
54. Rasmussen LD, et al. **Impact of Absolute Myocardial Blood Flow Quantification on the Diagnostic Performance of PET-Based Perfusion Scans Using 82Rubidium.** Circ Cardiovasc Imaging. 2024. [PMID 38227687](https://pubmed.ncbi.nlm.nih.gov/38227687/)
55. Soman P, et al. **Absolute Myocardial Blood Flow Quantification With PET: Should Diagnostic Cutoffs Be Tracer Specific?** Circ Cardiovasc Imaging. 2025. [PMID 40438936](https://pubmed.ncbi.nlm.nih.gov/40438936/)
56. Klein R, et al. **Selection of PET Camera and Implications on the Reliability and Accuracy of Absolute Myocardial Blood Flow Quantification.** Curr Cardiol Rep. 2020. [PMID 32770426](https://pubmed.ncbi.nlm.nih.gov/32770426/)
57. Hagemann CE, et al. **Quantitative myocardial blood flow with Rubidium-82 PET: a clinical perspective.** Am J Nucl Med Mol Imaging. 2015. [PMID 26550537](https://pubmed.ncbi.nlm.nih.gov/26550537/)
58. Bhave NM, et al. **Considerations when measuring myocardial perfusion reserve by cardiovascular magnetic resonance using regadenoson.** J Cardiovasc Magn Reson. 2012. [PMID 23272658](https://pubmed.ncbi.nlm.nih.gov/23272658/)
59. Wöhrle J, et al. **Myocardial perfusion reserve in cardiovascular magnetic resonance: Correlation to coronary microvascular dysfunction.** J Cardiovasc Magn Reson. 2006. [PMID 17060099](https://pubmed.ncbi.nlm.nih.gov/17060099/)
60. Taqueti VR, et al. **Clinical significance of noninvasive coronary flow reserve assessment in patients with ischemic heart disease.** Curr Opin Cardiol. 2016. [PMID 27652814](https://pubmed.ncbi.nlm.nih.gov/27652814/)
61. Taqueti VR, et al. **The role of positron emission tomography in the evaluation of myocardial ischemia in women.** J Nucl Cardiol. 2016. [PMID 27488383](https://pubmed.ncbi.nlm.nih.gov/27488383/)
62. Gunasekaran V, et al. **Assessment of coronary microvascular dysfunction in INOCA using 13N-ammonia PET: Lack of correlation with angiographic flow grades.** J Nucl Cardiol. 2026. [PMID 41407151](https://pubmed.ncbi.nlm.nih.gov/41407151/)
63. Singh H, et al. **Potential Role of 13N-NH3 Cardiac PET in Monitoring Treatment Response in Patients with Microvascular Angina.** Indian J Nucl Med. 2025. [PMID 40735747](https://pubmed.ncbi.nlm.nih.gov/40735747/)
64. Zaman MU, et al. **Cardiac Positron Emission Tomography Myocardial Perfusion Imaging: Seeing Beyond Perfusion.** World J Nucl Med. 2026. [PMID 42395164](https://pubmed.ncbi.nlm.nih.gov/42395164/)

---

## 4. 예후 (Prognosis)

65. **Kelshiker MA, et al. Coronary flow reserve and cardiovascular outcomes: a systematic review and meta-analysis. Eur Heart J. 2022.** [PMID 34849697](https://pubmed.ncbi.nlm.nih.gov/34849697/)
    **모델에서 사용한 수치:** 비정상 CFR에서 MACE HR 3.42(95% CI 2.92–3.99). **CFR 0.1 단위 감소당 사망 HR 1.16(1.04–1.29), MACE HR 1.08.** 비폐색성 관상동맥 환자에서 비정상 CFR의 사망 HR 5.44(3.78–7.83). 모델의 위험함수 `dCHMORT/dt = H0_MORT·exp(ln1.16·(2.5−CFR)/0.1)`가 이 기울기를 그대로 씁니다.
66. Gdowski MA, et al. **Association of Isolated Coronary Microvascular Dysfunction With Mortality and Major Adverse Cardiac Events: A Systematic Review and Meta-Analysis.** J Am Heart Assoc. 2020. [PMID 32345133](https://pubmed.ncbi.nlm.nih.gov/32345133/)
67. Luo X, et al. **Impact of Isolated Coronary Microvascular Disease Diagnosed Using Various Measurement Modalities on Prognosis: An Updated Systematic Review and Meta-Analysis.** Cardiology. 2024. [PMID 37708863](https://pubmed.ncbi.nlm.nih.gov/37708863/)
68. Gallinoro E, et al. **Prognostic Value of Microvascular Resistance Reserve in Coronary Artery Disease: A Systematic Review and Meta-Analysis.** JACC Cardiovasc Interv. 2026. [PMID 41881651](https://pubmed.ncbi.nlm.nih.gov/41881651/)
69. Seitz A, et al. **Prognostic implications of coronary artery stenosis and coronary spasm in patients with stable angina: 5-year follow-up of the ACOVA study.** Coron Artery Dis. 2020. [PMID 32168049](https://pubmed.ncbi.nlm.nih.gov/32168049/)
70. Weber BN, et al. **Impaired Coronary Vasodilator Reserve and Adverse Prognosis in Patients With Systemic Inflammatory Disorders.** JACC Cardiovasc Imaging. 2021. [PMID 33744132](https://pubmed.ncbi.nlm.nih.gov/33744132/)
71. Shah NR, et al. **Prognostic Value of Coronary Flow Reserve in Patients with Dialysis-Dependent ESRD.** J Am Soc Nephrol. 2016. [PMID 26459635](https://pubmed.ncbi.nlm.nih.gov/26459635/)
72. Huck DM, et al. **Comparative effectiveness of PET and SPECT myocardial perfusion imaging for predicting risk in patients with cardiometabolic disease.** J Nucl Cardiol. 2024. [PMID 38996910](https://pubmed.ncbi.nlm.nih.gov/38996910/)
73. Gould KL, et al. **Subendocardial and Transmural Myocardial Ischemia: Clinical Characteristics, Prevalence, and Outcomes With and Without Revascularization.** JACC Cardiovasc Imaging. 2023. [PMID 36599572](https://pubmed.ncbi.nlm.nih.gov/36599572/)
74. Eftekhari A, et al. **Changes in microvascular resistance following percutaneous coronary intervention — From the ILIAS global registry.** Int J Cardiol. 2023. [PMID 37633364](https://pubmed.ncbi.nlm.nih.gov/37633364/)
75. Jeyaprakash P, et al. **Index of Microcirculatory Resistance to predict microvascular obstruction in STEMI: systematic review and meta-analysis.** Catheter Cardiovasc Interv. 2024. [PMID 38179600](https://pubmed.ncbi.nlm.nih.gov/38179600/)
76. Zhang Y, et al. **Prognostic Value of Coronary Angiography-Derived Index of Microcirculatory Resistance in NSTEMI Patients.** JACC Cardiovasc Interv. 2024. [PMID 39115479](https://pubmed.ncbi.nlm.nih.gov/39115479/)
77. Zheng Y, et al. **Prognostic Value of Angiography-Derived IMR in Patients With Intermediate Coronary Stenosis.** JACC Cardiovasc Interv. 2025. [PMID 39880572](https://pubmed.ncbi.nlm.nih.gov/39880572/)
78. Chen D, et al. **Combined risk estimates of diabetes and angiography-derived IMR in NSTEMI.** Cardiovasc Diabetol. 2024. [PMID 39152477](https://pubmed.ncbi.nlm.nih.gov/39152477/)
79. Zhang Y, et al. **Prognostic Value of Angiography-Derived Index of Microvascular Resistance in Hypertrophic Cardiomyopathy.** MedComm. 2025. [PMID 40717902](https://pubmed.ncbi.nlm.nih.gov/40717902/)

---

## 5. 자동조절과 층별 관류 — 모델의 물리 (Autoregulation & transmural perfusion)

**"공급은 이완기에만 산다"** 는 모델의 두 번째 축이며, 다음 문헌들이 DPTI/SPTI
형식과 심박수의 이중 진입(demand와 관류시간)을 제공합니다.

80. **Duncker DJ, Bache RJ. Regulation of coronary blood flow during exercise. Physiol Rev. 2008;88:1009–86.** [PMID 18626066](https://pubmed.ncbi.nlm.nih.gov/18626066/) — 모델의 대사성 자동조절 구조, 안정 시 산소추출률 ~70%, 최대추출 상한, K_ATP·아데노신 매개 확장의 근거.
81. Duncker DJ, et al. **Regulation of coronary resistance vessel tone in response to exercise.** J Mol Cell Cardiol. 2012. [PMID 22037538](https://pubmed.ncbi.nlm.nih.gov/22037538/)
82. Duncker DJ, et al. **Role of K+ATP channels in coronary vasodilation during exercise.** Circulation. 1993. [PMID 8353886](https://pubmed.ncbi.nlm.nih.gov/8353886/)
83. Duncker DJ, et al. **Role of K+ATP channels and adenosine in the regulation of coronary blood flow during exercise with normal and restricted coronary blood flow.** J Clin Invest. 1996. [PMID 8613554](https://pubmed.ncbi.nlm.nih.gov/8613554/)
84. Hoffman JI, Buckberg GD. **Transmural myocardial perfusion.** Prog Cardiovasc Dis. 1987. [PMID 2953043](https://pubmed.ncbi.nlm.nih.gov/2953043/) — 심내막하가 이완기에만 관류되고 혈관외 압박을 시리즈 저항으로 받는다는 구조(모델 `W_ENDO`, `RCOMP_K`, `PHI_SYS`).
85. **Buckberg GD, et al. Ischemia in aortic stenosis: hemodynamic prediction. Am J Cardiol. 1975.** [PMID 1130286](https://pubmed.ncbi.nlm.nih.gov/1130286/) — DPTI/SPTI(=SEVR) 형식의 원전.
86. Brazier JR, et al. **Effects of tachycardia on the adequacy of subendocardial oxygen delivery in experimental aortic stenosis.** Am Heart J. 1975. [PMID 1155327](https://pubmed.ncbi.nlm.nih.gov/1155327/)
87. Brazier JR, et al. **Papillary muscle ischemia with patent coronary arteries.** Surgery. 1975. [PMID 1166409](https://pubmed.ncbi.nlm.nih.gov/1166409/)
88. Canty JM Jr, et al. **Effect of tachycardia on regional function and transmural myocardial perfusion during graded coronary pressure reduction in conscious dogs.** Circulation. 1990. [PMID 2225378](https://pubmed.ncbi.nlm.nih.gov/2225378/)
89. Buck JD, et al. **Changes in ischemic blood flow distribution and dynamic severity of a coronary stenosis induced by beta blockade in the canine heart.** Circulation. 1981. [PMID 6115724](https://pubmed.ncbi.nlm.nih.gov/6115724/)
90. Buck JD, et al. **Effects of sotalol and vagal stimulation on ischemic myocardial blood flow distribution in the canine heart.** J Pharmacol Exp Ther. 1981. [PMID 7463353](https://pubmed.ncbi.nlm.nih.gov/7463353/)
91. Chemla D, et al. **Subendocardial viability ratio estimated by arterial tonometry: a critical evaluation in elderly hypertensive patients with increased aortic stiffness.** Clin Exp Pharmacol Physiol. 2008. [PMID 18346166](https://pubmed.ncbi.nlm.nih.gov/18346166/)
92. Reitan JA, et al. **A computer evaluation of the ratio of the diastolic pressure-time index to the time-tension index from three arterial sites in dogs.** J Clin Monit. 1986. [PMID 3711953](https://pubmed.ncbi.nlm.nih.gov/3711953/)
93. Kissling G, et al. **Mechanical determinants of myocardial oxygen consumption with special reference to external work and efficiency.** Cardiovasc Res. 1992. [PMID 1451165](https://pubmed.ncbi.nlm.nih.gov/1451165/) — 모델 MVO2 식의 장력-시간 항.
94. Balady GJ, et al. **Comparison of determinants of myocardial oxygen consumption during arm and leg exercise in normal persons.** Am J Cardiol. 1986. [PMID 3717042](https://pubmed.ncbi.nlm.nih.gov/3717042/)
95. Richalet JP, et al. **Myocardial oxygen extraction and oxygen-hemoglobin equilibrium curve during moderate exercise.** Eur J Appl Physiol. 1981. [PMID 7197622](https://pubmed.ncbi.nlm.nih.gov/7197622/) — 모델의 `E_REST = 0.70`, `E_MAX = 0.80`.

---

## 6. 내피 생물학 — NO · BH4 · ROS · ADMA (Endothelial biology)

96. Yuyun MF, et al. **Endothelial dysfunction, endothelial nitric oxide bioavailability, tetrahydrobiopterin, and 5-methyltetrahydrofolate in cardiovascular disease.** Microvasc Res. 2018. [PMID 29596860](https://pubmed.ncbi.nlm.nih.gov/29596860/)
97. Bendall JK, et al. **Tetrahydrobiopterin in cardiovascular health and disease.** Antioxid Redox Signal. 2014. [PMID 24294830](https://pubmed.ncbi.nlm.nih.gov/24294830/) — 모델의 BH4↔BH2 산화와 eNOS 탈짝지음(`KBH4_OX`).
98. Cherng TW, et al. **Mechanisms of diesel-induced endothelial nitric oxide synthase dysfunction in coronary arterioles.** Environ Health Perspect. 2011. [PMID 20870565](https://pubmed.ncbi.nlm.nih.gov/20870565/)
99. Landim MB, et al. **Asymmetric dimethylarginine (ADMA) and endothelial dysfunction: implications for atherogenesis.** Clinics. 2009. [PMID 19488614](https://pubmed.ncbi.nlm.nih.gov/19488614/)
100. Böger RH. **Association of asymmetric dimethylarginine and endothelial dysfunction.** Clin Chem Lab Med. 2003. [PMID 14656027](https://pubmed.ncbi.nlm.nih.gov/14656027/)
101. Bełtowski J, et al. **Asymmetric dimethylarginine (ADMA) as a target for pharmacotherapy.** Pharmacol Rep. 2006. [PMID 16702618](https://pubmed.ncbi.nlm.nih.gov/16702618/)
102. Chan NN, et al. **ADMA: a potential link between endothelial dysfunction and cardiovascular diseases in insulin resistance syndrome?** Diabetologia. 2002. [PMID 12488950](https://pubmed.ncbi.nlm.nih.gov/12488950/)
103. Thengchaisri N, et al. **H2O2 Mediates VEGF- and Flow-Induced Dilations of Coronary Arterioles in Early Type 1 Diabetes: Role of Vascular Arginase and PI3K-Linked eNOS Uncoupling.** Int J Mol Sci. 2022. [PMID 36613929](https://pubmed.ncbi.nlm.nih.gov/36613929/)
104. Mahmoud AM, et al. **Nox2 contributes to hyperinsulinemia-induced redox imbalance and impaired vascular function.** Redox Biol. 2017. [PMID 28600985](https://pubmed.ncbi.nlm.nih.gov/28600985/)
105. Younis W, et al. **Soluble guanylyl cyclase, the NO receptor, drives vasorelaxation via endothelial S-nitrosation.** Proc Natl Acad Sci USA. 2025. [PMID 41037641](https://pubmed.ncbi.nlm.nih.gov/41037641/)
106. Friebe A, et al. **NO-GC in cells 'off the beaten track'.** Nitric Oxide. 2018. [PMID 29626542](https://pubmed.ncbi.nlm.nih.gov/29626542/)
107. Xiao S, et al. **Soluble Guanylate Cyclase Stimulators and Activators: Where are We and Where to Go?** Mini Rev Med Chem. 2019. [PMID 31362687](https://pubmed.ncbi.nlm.nih.gov/31362687/)
108. Torfgård KE, Ahlner J. **Mechanisms of action of nitrates.** Cardiovasc Drugs Ther. 1994. [PMID 7873467](https://pubmed.ncbi.nlm.nih.gov/7873467/)
109. Al-Badri A, et al. **Peripheral Microvascular Function Reflects Coronary Vascular Function.** Arterioscler Thromb Vasc Biol. 2019. [PMID 31018659](https://pubmed.ncbi.nlm.nih.gov/31018659/)
110. McChord J, et al. **Coronary Endothelial Dysfunction: Diagnostic Necessity or Futile Effort in Patients With Non-Obstructive Angina?** Catheter Cardiovasc Interv. 2025. [PMID 40745893](https://pubmed.ncbi.nlm.nih.gov/40745893/)

---

## 7. 엔도텔린 축 (Endothelin)

111. **Ford TJ, et al. Genetic dysregulation of endothelin-1 is implicated in coronary microvascular dysfunction. Eur Heart J. 2020.** [PMID 31972008](https://pubmed.ncbi.nlm.nih.gov/31972008/) — 협심증 391명에서 rs9349379-G가 혈중 ET-1 상승과 CMD에 연관. 모델 파라미터 `GENO`, `A_GENE_E = 0.30`의 근거이자 PRIZE의 정밀의료 가설.
112. Feng J, et al. **Endothelin-1-induced contractile responses of human coronary arterioles via endothelin-A receptors and PKC-alpha signaling pathways.** Surgery. 2010. [PMID 20079914](https://pubmed.ncbi.nlm.nih.gov/20079914/) — ETA 매개 세동맥 수축(모델 `F_ETA = 0.75`, `F_ETB2 = 0.25`).
113. Dashwood MR, et al. **Regional variations in endothelin-1 and its receptor subtypes in human coronary vasculature.** Endothelium. 1998. [PMID 9832333](https://pubmed.ncbi.nlm.nih.gov/9832333/)
114. DeFily DV, et al. **Endothelin antagonists block alpha1-adrenergic constriction of coronary arterioles.** Am J Physiol. 1999. [PMID 10070088](https://pubmed.ncbi.nlm.nih.gov/10070088/) — α1과 ET 축의 결합(모델 `K_A1_TN`).
115. Lamping KG, et al. **Effects of 17 beta-estradiol on coronary microvascular responses to endothelin-1.** Am J Physiol. 1996. [PMID 8853349](https://pubmed.ncbi.nlm.nih.gov/8853349/)
116. Sauvageau S, et al. **Evaluation of endothelin-1-induced pulmonary vasoconstriction following myocardial infarction.** Exp Biol Med. 2006. [PMID 16741009](https://pubmed.ncbi.nlm.nih.gov/16741009/)

---

## 8. Rho-kinase와 경련 (Rho-kinase & spasm)

117. Shimokawa H. **2014 Williams Harvey Lecture: importance of coronary vasomotion abnormalities—from bench to bedside.** Eur Heart J. 2014. [PMID 25354517](https://pubmed.ncbi.nlm.nih.gov/25354517/)
118. Shimokawa H. **Cellular and molecular mechanisms of coronary artery spasm.** Jpn Circ J. 2000. [PMID 10651199](https://pubmed.ncbi.nlm.nih.gov/10651199/) — 모델의 ROCK→MLC 탈인산화 억제→Ca 비의존 감작.
119. Yoo SY, Kim JY. **Recent insights into the mechanisms of vasospastic angina.** Korean Circ J. 2009. [PMID 20049135](https://pubmed.ncbi.nlm.nih.gov/20049135/)
120. Oi K, et al. **Remnant lipoproteins from patients with sudden cardiac death enhance coronary vasospastic activity through upregulation of Rho-kinase.** Arterioscler Thromb Vasc Biol. 2004. [PMID 15044207](https://pubmed.ncbi.nlm.nih.gov/15044207/) — 모델에서 LDL/잔여지단백이 `ROCK_D`를 올리는 경로.
121. **Masumoto A, et al. Suppression of coronary artery spasm by the Rho-kinase inhibitor fasudil in patients with vasospastic angina. Circulation. 2002.** [PMID 11927519](https://pubmed.ncbi.nlm.nih.gov/11927519/) — 모델 파수딜 팔의 근거(`EM_FAS_RK = 0.72`).
122. Otsuka T, et al. **Administration of the Rho-kinase inhibitor fasudil following nitroglycerin additionally dilates the site of coronary spasm.** Coron Artery Dis. 2008. [PMID 18300747](https://pubmed.ncbi.nlm.nih.gov/18300747/) — 나이트레이트로 풀리지 않는 경련 성분이 남는다는 관찰(모델에서 나이트레이트가 심표면에만 작용하는 이유).
123. Mohri M, et al. **Angina pectoris caused by coronary microvascular spasm.** Lancet. 1998. [PMID 9643687](https://pubmed.ncbi.nlm.nih.gov/9643687/)
124. Sun H, et al. **Coronary microvascular spasm causes myocardial ischemia in patients with vasospastic angina.** J Am Coll Cardiol. 2002. [PMID 11869851](https://pubmed.ncbi.nlm.nih.gov/11869851/) — 미세혈관 경련이 협착 없이 허혈을 만든다 → 모델이 경련을 "긴장도"가 아니라 **혈관 폐쇄 + 구동압 손실**로 쓰는 이유(bug B20).
125. Ong P, et al. **High prevalence of a pathological response to acetylcholine testing in patients with stable angina pectoris and unobstructed coronary arteries. The ACOVA Study.** J Am Coll Cardiol. 2012. [PMID 22322081](https://pubmed.ncbi.nlm.nih.gov/22322081/)
126. Seitz A, et al. **Characterization and implications of intracoronary hemodynamic assessment during coronary spasm provocation testing.** Clin Res Cardiol. 2023. [PMID 37195455](https://pubmed.ncbi.nlm.nih.gov/37195455/)
127. Feenstra RGT, et al. **Haemodynamic characterisation of different endotypes in coronary artery vasospasm in reaction to acetylcholine.** Int J Cardiol Heart Vasc. 2022. [PMID 36017267](https://pubmed.ncbi.nlm.nih.gov/36017267/)
128. Feenstra RGT, et al. **Post-spastic flow recovery time to document vasospasm induced ischemia during acetylcholine provocation testing.** Int J Cardiol Heart Vasc. 2023. [PMID 37275626](https://pubmed.ncbi.nlm.nih.gov/37275626/)
129. Feenstra RGT, et al. **Do ECG changes induced during intracoronary vasospasm provocation testing reflect those during spontaneous angina episodes in vasospastic angina?** Eur Heart J Case Rep. 2024. [PMID 39161720](https://pubmed.ncbi.nlm.nih.gov/39161720/)
130. Huang J, et al. **Invasive Evaluation for Coronary Vasospasm.** US Cardiol. 2023. [PMID 39493950](https://pubmed.ncbi.nlm.nih.gov/39493950/)
131. Aswathappa S, et al. **A Comprehensive Literature Review Discussing Diagnostic Challenges of Prinzmetal or Vasospastic Angina.** Cureus. 2025. [PMID 40486459](https://pubmed.ncbi.nlm.nih.gov/40486459/)
132. Mehta HH, et al. **The Spontaneous Coronary Slow-Flow Phenomenon: Reversal by Intracoronary Nicardipine.** J Invasive Cardiol. 2019. [PMID 30555052](https://pubmed.ncbi.nlm.nih.gov/30555052/)
133. Sykes R, et al. **Myocardial Bridging Independently Associates With Coronary Artery Spasm.** JACC Cardiovasc Interv. 2026. [PMID 42508854](https://pubmed.ncbi.nlm.nih.gov/42508854/)
134. Toya T, et al. **Coronary Endothelial Dysfunction and Vasomotor Dysregulation in Myocardial Bridging.** J Cardiovasc Dev Dis. 2025. [PMID 39997488](https://pubmed.ncbi.nlm.nih.gov/39997488/)

---

## 9. 아데노신 · K_ATP · 통각 (Adenosine, K_ATP, nociception)

이 절은 모델의 가장 비자명한 예측을 지탱합니다. **기능적 endotype은 최소저항이
정상이므로 어떤 부하에서도 심내막하 공급-요구 결손이 사실상 0이고, 따라서 그
협심증은 허혈이 아니라 구심성(A1-아데노신) 신호와 중추 감작이 운반한다.**

135. **Elliott PM, et al. Effect of oral aminophylline in patients with angina and normal coronary arteriograms (cardiac syndrome X). Heart. 1997;77:523.** [PMID 9227295](https://pubmed.ncbi.nlm.nih.gov/9227295/) — 모델이 재현하는 관찰: 아데노신 수용체 차단이 CFR을 거의 움직이지 않으면서 운동시간을 늘린다. 모델은 이 효과가 **기능적 endotype에 국한**되고 구조적 endotype에서는 A2A 차단으로 확장 예비력을 잃어 해로울 수 있다고 예측합니다.
136. Zhou X, et al. **A1 adenosine receptor negatively modulates coronary reactive hyperemia via counteracting A2A-mediated H2O2 production and KATP opening.** Am J Physiol Heart Circ Physiol. 2013. [PMID 24043252](https://pubmed.ncbi.nlm.nih.gov/24043252/) — A1과 A2A의 반대 작용(모델 `EM_AMI_A1` 대 `EM_AMI_A2`).
137. Peart JN, Headrick JP. **Adenosinergic cardioprotection: multiple receptors, multiple pathways.** Pharmacol Ther. 2007. [PMID 17408751](https://pubmed.ncbi.nlm.nih.gov/17408751/)
138. Riou LM, et al. **Influence of propranolol, enalaprilat, verapamil, and caffeine on adenosine A2A-receptor-mediated coronary vasodilation.** J Am Coll Cardiol. 2002. [PMID 12427424](https://pubmed.ncbi.nlm.nih.gov/12427424/) — 카페인/메틸잔틴이 아데노신 확장을 방해한다 → 진단 검사와 치료의 충돌.
139. Niiya K, et al. **Glibenclamide reduces the coronary vasoactivity of adenosine receptor agonists.** J Pharmacol Exp Ther. 1994. [PMID 7965706](https://pubmed.ncbi.nlm.nih.gov/7965706/)
140. Lanza GA, et al. **Effect of spinal cord stimulation on spontaneous and stress-induced angina and 'ischemia-like' ST-segment depression in patients with cardiac syndrome X.** Eur Heart J. 2005. [PMID 15642701](https://pubmed.ncbi.nlm.nih.gov/15642701/)
141. Lanza GA, et al. **Spinal cord stimulation in patients with refractory anginal pain and normal coronary arteries.** Ital Heart J. 2001. [PMID 11214698](https://pubmed.ncbi.nlm.nih.gov/11214698/)
142. Eliasson T, et al. **Spinal cord stimulation in angina pectoris with normal coronary arteriograms.** Coron Artery Dis. 1993. [PMID 8287216](https://pubmed.ncbi.nlm.nih.gov/8287216/)
143. Lanza GA, et al. **Management of microvascular angina pectoris.** Am J Cardiovasc Drugs. 2014. [PMID 24174173](https://pubmed.ncbi.nlm.nih.gov/24174173/)
144. de Silva R, et al. **Refractory angina: mechanisms and stratified treatment in obstructive and non-obstructive chronic myocardial ischaemic syndromes.** Eur Heart J. 2025. [PMID 40590516](https://pubmed.ncbi.nlm.nih.gov/40590516/)
145. Tyrer P, et al. **Cognitive behaviour therapy for non-cardiac pain in the chest (COPIC): a multicentre randomized controlled trial with economic evaluation.** BMC Psychol. 2015. [PMID 26596540](https://pubmed.ncbi.nlm.nih.gov/26596540/) — 모델 `CBT` 스위치.
146. Eriksson-Liebon M, et al. **Long-term effects and predictors of change of internet-delivered CBT on cardiac anxiety in patients with non-cardiac chest pain: RCT.** BMC Psychiatry. 2024. [PMID 38504157](https://pubmed.ncbi.nlm.nih.gov/38504157/)
147. Thesen T, et al. **Patients with depression symptoms are more likely to experience improvements of internet-based CBT: secondary analysis in non-cardiac chest pain.** BMC Psychiatry. 2023. [PMID 37838653](https://pubmed.ncbi.nlm.nih.gov/37838653/)
148. Achem SR. **Recent developments in chest pain of undetermined origin.** Curr Gastroenterol Rep. 2000. [PMID 10957931](https://pubmed.ncbi.nlm.nih.gov/10957931/)
149. Shrestha S, Pasricha PJ. **Update on noncardiac chest pain.** Dig Dis. 2000. [PMID 11279332](https://pubmed.ncbi.nlm.nih.gov/11279332/)

---

## 10. 구조 리모델링 · 동반질환 (Structural remodelling & comorbidity)

150. Camici PG, et al. **Coronary microvascular dysfunction in hypertrophy and heart failure.** Cardiovasc Res. 2020. [PMID 31999329](https://pubmed.ncbi.nlm.nih.gov/31999329/) — 모델의 LVH·모세혈관 밀도 불균형(`CAPD`, `LVH`).
151. Paulus WJ, Tschöpe C. **A novel paradigm for heart failure with preserved ejection fraction: comorbidities drive myocardial dysfunction and remodeling through coronary microvascular endothelial inflammation.** J Am Coll Cardiol. 2013. [PMID 23684677](https://pubmed.ncbi.nlm.nih.gov/23684677/)
152. **Shah SJ, et al. Prevalence and correlates of coronary microvascular dysfunction in heart failure with preserved ejection fraction: PROMIS-HFpEF. Eur Heart J. 2018.** [PMID 30165580](https://pubmed.ncbi.nlm.nih.gov/30165580/) — CMD와 충만압의 결합. 모델의 운동 시 LVEDP 상승(`K_LVDP_W`, `K_ICF_LW`)이 심내막하 관류를 깎는 경로의 임상적 근거.
153. Sinha A, et al. **Coronary microvascular dysfunction and heart failure with preserved ejection fraction: what are the mechanistic links?** Curr Opin Cardiol. 2023. [PMID 37668191](https://pubmed.ncbi.nlm.nih.gov/37668191/)
154. Erhardsson M, et al. **Regional differences and coronary microvascular dysfunction in heart failure with preserved ejection fraction.** ESC Heart Fail. 2023. [PMID 37920127](https://pubmed.ncbi.nlm.nih.gov/37920127/)
155. Chandramouli C, et al. **Sex differences in proteomic correlates of coronary microvascular dysfunction among patients with HFpEF.** Eur J Heart Fail. 2022. [PMID 35060248](https://pubmed.ncbi.nlm.nih.gov/35060248/)
156. Venkateshvaran A, et al. **Association of epicardial adipose tissue with proteomics, coronary flow reserve, cardiac structure and function, and quality of life in HFpEF: PROMIS-HFpEF.** Eur J Heart Fail. 2022. [PMID 36196462](https://pubmed.ncbi.nlm.nih.gov/36196462/)
157. Mahmoud I, et al. **Epicardial adipose tissue differentiates in patients with and without coronary microvascular dysfunction.** Int J Obes. 2021. [PMID 34172829](https://pubmed.ncbi.nlm.nih.gov/34172829/)
158. Patel NH, et al. **Epicardial adipose tissue attenuation on CT in women with coronary microvascular dysfunction.** Atherosclerosis. 2024. [PMID 38944545](https://pubmed.ncbi.nlm.nih.gov/38944545/)
159. Agabiti-Rosei E, Rizzoni D. **[Structural and functional changes of the microcirculation in hypertension: influence of pharmacological therapy].** Drugs. 2003. [PMID 12708883](https://pubmed.ncbi.nlm.nih.gov/12708883/) — media/lumen 내향 리모델링과 ACE 억제제에 의한 부분 역전(모델 `EM_RAM_ML`, `TAU_ML`).
160. Feihl F, et al. **The macrocirculation and microcirculation of hypertension.** Curr Hypertens Rep. 2009. [PMID 19442327](https://pubmed.ncbi.nlm.nih.gov/19442327/)
161. Agabiti-Rosei E, Rizzoni D. **From macro- to microcirculation: benefits in hypertension and diabetes.** J Hypertens Suppl. 2008. [PMID 19363848](https://pubmed.ncbi.nlm.nih.gov/19363848/)
162. Sezer M, et al. **Bimodal Pattern of Coronary Microvascular Involvement in Diabetes Mellitus.** J Am Heart Assoc. 2016. [PMID 27930353](https://pubmed.ncbi.nlm.nih.gov/27930353/)
163. Niewiara Ł, et al. **Impaired coronary flow reserve in patients with poor type 2 diabetes control.** Cardiol J. 2024. [PMID 36342032](https://pubmed.ncbi.nlm.nih.gov/36342032/)
164. Huang R, et al. **Relationship between glycosylated hemoglobin A1c and coronary flow reserve in patients with type 2 diabetes.** Expert Rev Cardiovasc Ther. 2015. [PMID 25695762](https://pubmed.ncbi.nlm.nih.gov/25695762/)
165. Y-Hassan S, et al. **Coronary microvascular dysfunction in Takotsubo syndrome: cause or consequence.** Am J Cardiovasc Dis. 2021. [PMID 34084653](https://pubmed.ncbi.nlm.nih.gov/34084653/)
166. Castaldi G, et al. **Angiography-derived index of microvascular resistance in takotsubo syndrome.** Int J Cardiovasc Imaging. 2023. [PMID 36336756](https://pubmed.ncbi.nlm.nih.gov/36336756/)
167. Chitturi KR, et al. **Coronary microvascular dysfunction and cancer therapy-related cardiovascular toxicity.** Cardiovasc Revasc Med. 2024. [PMID 38789343](https://pubmed.ncbi.nlm.nih.gov/38789343/)
168. Türkoğlu C, et al. **The Relationship Between H2FPEF Score and Coronary Slow Flow Phenomenon.** Turk Kardiyol Dern Ars. 2022. [PMID 35695359](https://pubmed.ncbi.nlm.nih.gov/35695359/)
169. Amirzadegan A, et al. **Coronary slow flow phenomenon and microalbuminuria.** Turk Kardiyol Dern Ars. 2019. [PMID 31802772](https://pubmed.ncbi.nlm.nih.gov/31802772/)

---

## 11. 무작위 임상시험과 약물치료 (Randomised trials & therapeutics)

모델이 재현해야 하는 다섯 개의 임상시험 정박점(anchor)입니다. 세 개는 **음성**
결과이고, 모델이 하는 일의 절반은 그 음성 결과들이 왜 음성이었는지를
정량적으로 분해하는 것입니다.

### 11.1 계층화 치료 — CorMicA

170. **Ford TJ, et al. Stratified Medical Therapy Using Invasive Coronary Function Testing in Angina: The CorMicA Trial. J Am Coll Cardiol. 2018;72:2841–2855.** [PMID 30266608](https://pubmed.ncbi.nlm.nih.gov/30266608/)
     **모델에서 사용한 수치:** 391명 등록, 관상동맥조영에서 폐색성 병변 206명(53.7%), 비폐색성 151명(39%)을 1:1 무작위(중재 76 / 눈가림 대조 75). 중재는 guidewire 기반 CFR·IMR·FFR + ACh 혈관반응성 검사에 연동된 계층화 치료. **6개월 SAQ 요약점수 평균 +11.7 U(95% CI 5.0–18.4, p=0.001)**, EQ-5D +0.10(0.01–0.18), VAS +14.5(7.8–21.3). 6개월 MACE 차이 없음(2.6% 대 2.6%). 모델의 §X가 이 수치를 목표로 삼고 +5.4 U를 산출합니다(방향 일치, 크기 과소).
171. Ford TJ, et al. **Rationale and design of the BHF CorMicA stratified medicine clinical trial.** Am Heart J. 2018. [PMID 29803987](https://pubmed.ncbi.nlm.nih.gov/29803987/)
172. Ford TJ, et al. **How to Diagnose and Manage Angina Without Obstructive Coronary Artery Disease: Lessons from CorMicA.** Interv Cardiol. 2019. [PMID 31178933](https://pubmed.ncbi.nlm.nih.gov/31178933/)
173. Heggie R, et al. **Stratified medicine using invasive coronary function testing in angina: A cost-effectiveness analysis of the BHF CorMicA trial.** Int J Cardiol. 2021. [PMID 33992700](https://pubmed.ncbi.nlm.nih.gov/33992700/)

### 11.2 라놀라진 — RWISE

174. **Bairey Merz CN, et al. A randomized, placebo-controlled trial of late Na current inhibition (ranolazine) in coronary microvascular dysfunction (CMD): impact on angina and myocardial perfusion reserve. Eur Heart J. 2016;37:1504–13.** [PMID 26614823](https://pubmed.ncbi.nlm.nih.gov/26614823/)
     **모델에서 사용한 수치:** 128명(96% 여성), 라놀라진 500–1000 mg 1일 2회 2주, 이중맹검 교차. **전체적으로 SAQ 차이 없음.** 부하 시 최고 심박수 −3.55 bpm(p<0.001) → 모델 `EM_RAN_HR = 0.052`. SAQ-7 변화와 MPRI 변화의 상관 0.25(p=0.005). **CFR<2.5 하위군에서만** MPRI(p=0.014)·협심증 빈도(p=0.027)·SAQ-7(p=0.041) 개선. 모델의 §VII이 재현하는 희석: CMD 계층 +2.3 U 대 전체 코호트 +1.2 U(둘 다 MCID 미만).
175. Hampilos KE, et al. **Myocardial biomarkers in coronary microvascular dysfunction: Response to ranolazine.** Am Heart J Plus. 2025. [PMID 40093309](https://pubmed.ncbi.nlm.nih.gov/40093309/)
176. Zhu H, et al. **Effects of the Antianginal Drugs Ranolazine, Nicorandil, and Ivabradine on Coronary Microvascular Function in Patients With Nonobstructive Coronary Artery Disease: A Meta-analysis of RCTs.** Clin Ther. 2019. [PMID 31548105](https://pubmed.ncbi.nlm.nih.gov/31548105/)
177. Patel S, et al. **Contemporary Antianginal Therapy.** Am J Cardiovasc Drugs. 2026. [PMID 40999181](https://pubmed.ncbi.nlm.nih.gov/40999181/)

### 11.3 지보텐탄 — PRIZE

178. **Morrow A, et al. Zibotentan in Microvascular Angina: A Randomized, Placebo-Controlled, Crossover Trial. Circulation. 2024.** [PMID 39217504](https://pubmed.ncbi.nlm.nih.gov/39217504/)
     **모델에서 사용한 수치:** 미세혈관 협심증 118명(63.5±9.2세, 60.2% 여성, 21.2% 당뇨), rs9349379-G 대립유전자 빈도 50%로 강화, 지보텐탄 10 mg/일 12주, 순차 교차. 완전자료 103명에서 **트레드밀(Bruce) 지속시간 차이 −4.26초(95% CI −19.60 ~ +11.06, p=0.5871)**, 2차 지표 모두 개선 없음. **지보텐탄은 혈압을 낮추고 혈중 ET-1을 올렸다.** 이상반응 60.2% 대 위약 14.4%(p<0.001), 체액 보유 우세. 모델의 §VI가 −4.53초를 산출하고(관측 신뢰구간 내부) 그 전부를 체액 보유(−6.99초)에 귀속시킵니다. 혈압 강하는 오히려 +2.41초(요구 감소가 구동압 손실을 앞선다).
179. Morrow AJ, et al. **Rationale and design of the MRC's Precision Medicine with Zibotentan in Microvascular Angina (PRIZE) trial.** Am Heart J. 2020. [PMID 32942043](https://pubmed.ncbi.nlm.nih.gov/32942043/)
180. Morrow A, et al. **Exercise treadmill testing for efficacy evaluation in randomized, controlled trials.** Am Heart J. 2026. [PMID 41687797](https://pubmed.ncbi.nlm.nih.gov/41687797/) — 모델의 §VI 결론(트레드밀 지속시간은 이 모집단에서 둔한 기구)과 직접 관련.
181. Pasupathy S, et al. **Anti-Anginal Efficacy of Zibotentan in the Coronary Slow-Flow Phenomenon.** J Clin Med. 2024. [PMID 38592159](https://pubmed.ncbi.nlm.nih.gov/38592159/)

### 11.4 집중적 약물치료 — WARRIOR

182. **Pepine CJ, et al. Women's IschemiA TRial to Reduce Events In Non-ObstRuctive CAD (WARRIOR): a randomised controlled trial. Open Heart. 2026;13:e004115.** [PMID 41932694](https://pubmed.ncbi.nlm.nih.gov/41932694/)
     **모델에서 사용한 수치:** ANOCA/INOCA 의심 여성 2476명, 71개 기관. 집중치료(고강도 스타틴 + ACEi/ARB + 아스피린) 대 통상치료. 2.5년에 421건(집중 221 / 통상 200), **1차 종점 HR 1.13(95% CI 0.94–1.37, p=0.20)**, 2차 종점도 차이 없음. **협심증 입원이 MACE의 주 구성요소.** 오염 보정 민감도 분석 HR 0.74(0.352–1.558, p=0.43). 등록이 계획보다 적어 고령(평균 64세)·혈압/LDL 양호·스타틴·ACEi/ARB 사용률 높은 집단. 모델의 §IX가 참 HR 0.90을 내고 배경 오염 80%에서 관측 HR 0.98이 됨을 보입니다.
183. Lakshmanan S, et al. **Comparison of risk profiles of participants in the WARRIOR trial, using CCTA vs invasive coronary angiography.** Prog Cardiovasc Dis. 2024. [PMID 38547955](https://pubmed.ncbi.nlm.nih.gov/38547955/)

### 11.5 심박수 조절 · 칼슘차단제 · 기타 약물

184. Camici PG, et al. **Ivabradine in chronic stable angina: Effects by and beyond heart rate reduction.** Int J Cardiol. 2016. [PMID 27104917](https://pubmed.ncbi.nlm.nih.gov/27104917/) — 모델의 §III(심박수 감소 이익의 43%가 이완기 관류시간)와 직접 대응.
185. Heusch G. **Ivabradine: Cardioprotection By and Beyond Heart Rate Reduction.** Drugs. 2016. [PMID 27041289](https://pubmed.ncbi.nlm.nih.gov/27041289/)
186. Giavarini A, et al. **The Role of Ivabradine in the Management of Angina Pectoris.** Cardiovasc Drugs Ther. 2016. [PMID 27475447](https://pubmed.ncbi.nlm.nih.gov/27475447/)
187. Bucchi A, et al. **Heart rate reduction via selective 'funny' channel blockers.** Curr Opin Pharmacol. 2007. [PMID 17267284](https://pubmed.ncbi.nlm.nih.gov/17267284/)
188. Borer JS, et al. **Characterization of the heart rate-lowering action of ivabradine, a selective I(f) current inhibitor.** Am J Ther. 2008. [PMID 18806523](https://pubmed.ncbi.nlm.nih.gov/18806523/) — 모델 이바브라딘 PK/PD.
189. Chaudhary R, et al. **Ivabradine: Heart Failure and Beyond.** J Cardiovasc Pharmacol Ther. 2016. [PMID 26721645](https://pubmed.ncbi.nlm.nih.gov/26721645/)
190. Doesch AO, et al. **Heart rate reduction after heart transplantation with beta-blocker versus the selective If channel antagonist ivabradine.** Transplantation. 2007. [PMID 17989604](https://pubmed.ncbi.nlm.nih.gov/17989604/)
191. Rognoni A, et al. **Ivabradine: cardiovascular effects.** Recent Pat Cardiovasc Drug Discov. 2009. [PMID 19149708](https://pubmed.ncbi.nlm.nih.gov/19149708/)
192. Chen JW, et al. **Effects of short-term treatment of nicorandil on exercise-induced myocardial ischemia and abnormal cardiac autonomic activity in microvascular angina.** Am J Cardiol. 1997. [PMID 9205016](https://pubmed.ncbi.nlm.nih.gov/9205016/)
193. Hirohata A, et al. **Nicorandil prevents microvascular dysfunction resulting from PCI in patients with stable angina pectoris: a randomised study.** EuroIntervention. 2014. [PMID 24457276](https://pubmed.ncbi.nlm.nih.gov/24457276/)
194. Zhang Y, et al. **The effectiveness and safety of nicorandil in the treatment of patients with microvascular angina: protocol for systematic review and meta-analysis.** Medicine. 2021. [PMID 33466132](https://pubmed.ncbi.nlm.nih.gov/33466132/)
195. Jia Q, et al. **The effect of nicorandil in patients with cardiac syndrome X: a meta-analysis of RCTs.** Medicine. 2020. [PMID 32925783](https://pubmed.ncbi.nlm.nih.gov/32925783/)
196. Pavão RB, et al. **Aspirin plus verapamil relieves angina and perfusion abnormalities in patients with coronary microvascular dysfunction and Chagas disease.** Rev Soc Bras Med Trop. 2021. [PMID 34787258](https://pubmed.ncbi.nlm.nih.gov/34787258/)
197. Denardo SJ, et al. **Effect of phosphodiesterase type 5 inhibition on microvascular coronary dysfunction in women: a WISE ancillary study.** Clin Cardiol. 2011. [PMID 21780138](https://pubmed.ncbi.nlm.nih.gov/21780138/) — 모델 실데나필 팔.
198. Guarini G, et al. **Trimetazidine and Other Metabolic Modifiers.** Eur Cardiol. 2018. [PMID 30697354](https://pubmed.ncbi.nlm.nih.gov/30697354/) — 모델 `TMZ` 스위치(같은 ATP를 더 적은 산소로).
199. Nalbantgil S, et al. **The Effect of Trimetazidine in the Treatment of Microvascular Angina.** Int J Angiol. 1999. [PMID 9826407](https://pubmed.ncbi.nlm.nih.gov/9826407/)
200. Lanza GA. **[Therapy of microvascular angina].** Cardiologia. 1993. [PMID 7912650](https://pubmed.ncbi.nlm.nih.gov/7912650/)
201. Ferrara L, et al. **[Syndrome X and microvascular angina].** Minerva Cardioangiol. 1998. [PMID 9882962](https://pubmed.ncbi.nlm.nih.gov/9882962/)
202. Rusticali G, et al. **[The noninvasive identification of patients with angina and normal coronary arteries].** G Ital Cardiol. 1995. [PMID 8529853](https://pubmed.ncbi.nlm.nih.gov/8529853/)
203. Hinoi T, et al. **Acute effect of atorvastatin on coronary circulation measured by transthoracic Doppler echocardiography in patients without coronary artery disease by angiography.** Am J Cardiol. 2005. [PMID 15979441](https://pubmed.ncbi.nlm.nih.gov/15979441/)
204. Bountioukos M, et al. **Effect of atorvastatin on myocardial contractile reserve assessed by tissue Doppler imaging in moderately hypercholesterolemic patients without heart disease.** Am J Cardiol. 2003. [PMID 12943890](https://pubmed.ncbi.nlm.nih.gov/12943890/)
205. Mortensen MB, et al. **Influence of intensive lipid-lowering on CT-derived fractional flow reserve in patients with stable chest pain: FLOWPROMOTE.** Clin Cardiol. 2022. [PMID 36056636](https://pubmed.ncbi.nlm.nih.gov/36056636/)

---

## 12. 비약물·기구 치료 (Non-pharmacological & device therapy)

206. Kissel CK, Nikoletou D. **Cardiac Rehabilitation and Exercise Prescription in Symptomatic Patients with Non-Obstructive Coronary Artery Disease—a Systematic Review.** Curr Treat Options Cardiovasc Med. 2018. [PMID 30121850](https://pubmed.ncbi.nlm.nih.gov/30121850/) — 모델 `REHAB` 스위치.
207. Hausvater A, et al. **Cardiac Rehabilitation for Patients With INOCA and MINOCA: A Review.** J Cardiopulm Rehabil Prev. 2025. [PMID 40476778](https://pubmed.ncbi.nlm.nih.gov/40476778/)
208. Carvalho EE, et al. **Improved endothelial function and reversal of myocardial perfusion defects after aerobic physical training in a patient with microvascular myocardial ischemia.** Am J Phys Med Rehabil. 2011. [PMID 20531160](https://pubmed.ncbi.nlm.nih.gov/20531160/)
209. Akyuz A. **Exercise and Coronary Heart Disease.** Adv Exp Med Biol. 2020. [PMID 32342457](https://pubmed.ncbi.nlm.nih.gov/32342457/)
210. Hong J, et al. **Exercise training mitigates ER stress and UCP2 deficiency-associated coronary vascular dysfunction in atherosclerosis.** Sci Rep. 2021. [PMID 34326395](https://pubmed.ncbi.nlm.nih.gov/34326395/)
211. Tryon D, et al. **Coronary Sinus Reducer Improves Angina, Quality of Life, and Coronary Flow Reserve in Microvascular Dysfunction.** JACC Cardiovasc Interv. 2024. [PMID 39520443](https://pubmed.ncbi.nlm.nih.gov/39520443/) — 관상정맥압을 올려 심내막하 관류를 재분배한다는 발상. 모델의 `PV`와 심내막하 구동압 항이 이 기전을 담고 있습니다.
212. Tebaldi M, et al. **Coronary Sinus Narrowing Improves Coronary Microcirculation Function in Patients With Refractory Angina: INROAD.** Circ Cardiovasc Interv. 2024. [PMID 38227697](https://pubmed.ncbi.nlm.nih.gov/38227697/)
213. Konigstein M, et al. **Coronary Sinus Narrowing for the Treatment of Patients With Angina and Evidence of Microvascular Dysfunction.** Can J Cardiol. 2026. [PMID 42309353](https://pubmed.ncbi.nlm.nih.gov/42309353/)
214. Tomaniak M, et al. **Coronary Sinus Reduction for Refractory Angina Caused by Microvascular Dysfunction—A Systematic Review.** J Clin Med. 2025. [PMID 41517541](https://pubmed.ncbi.nlm.nih.gov/41517541/)
215. Ashokprabhu ND, et al. **Enhanced External Counterpulsation for the Treatment of Angina With Nonobstructive Coronary Artery Disease.** Am J Cardiol. 2024. [PMID 37890564](https://pubmed.ncbi.nlm.nih.gov/37890564/)
216. Kronhaus KD, Lawson WE. **Enhanced external counterpulsation is an effective treatment for Syndrome X.** Int J Cardiol. 2009. [PMID 18590931](https://pubmed.ncbi.nlm.nih.gov/18590931/)
217. Bondesson SM, et al. **Reduced peripheral vascular reactivity in refractory angina pectoris: Effect of enhanced external counterpulsation.** J Geriatr Cardiol. 2011. [PMID 22783308](https://pubmed.ncbi.nlm.nih.gov/22783308/)

---

## 13. 성별 · 폐경 · 집단 (Sex, menopause, population)

218. Waheed N, et al. **Sex differences in non-obstructive coronary artery disease.** Cardiovasc Res. 2020. [PMID 31958135](https://pubmed.ncbi.nlm.nih.gov/31958135/)
219. Jansen TPJ, et al. **Sex Differences in Coronary Function Test Results in Patients With Angina and Nonobstructive Disease.** Front Cardiovasc Med. 2021. [PMID 34722680](https://pubmed.ncbi.nlm.nih.gov/34722680/)
220. Steinberg RR, et al. **Coronary microvascular disease in women: epidemiology, mechanisms, evaluation, and treatment.** Can J Physiol Pharmacol. 2024. [PMID 38728748](https://pubmed.ncbi.nlm.nih.gov/38728748/)
221. Mathew D, et al. **Coronary microvascular dysfunction in menopausal women.** Heart. 2026. [PMID 42331611](https://pubmed.ncbi.nlm.nih.gov/42331611/) — 모델 `RF_MENO`(에스트로겐 소실 → eNOS).
222. Shufelt CL, et al. **Sex-Specific Physiology and Cardiovascular Disease.** Adv Exp Med Biol. 2018. [PMID 30051400](https://pubmed.ncbi.nlm.nih.gov/30051400/)
223. SenthilKumar G, et al. **17β-Estradiol promotes sex-specific dysfunction in isolated human arterioles.** Am J Physiol Heart Circ Physiol. 2023. [PMID 36607795](https://pubmed.ncbi.nlm.nih.gov/36607795/)
224. Lam CSP, et al. **Sex differences in heart failure.** Eur Heart J. 2019. [PMID 31800034](https://pubmed.ncbi.nlm.nih.gov/31800034/)
225. Carlini NA, et al. **Vascular function in women with heart failure with preserved ejection fraction: a mismatch beyond diastole.** J Appl Physiol. 2025. [PMID 40839391](https://pubmed.ncbi.nlm.nih.gov/40839391/)
226. Ryk-Adamska M, et al. **Ophthalmological Microvascular Changes in ANOCA/INOCA Disease and Ophthalmological Methods to Detect Them—A Systematic Review.** J Clin Med. 2026. [PMID 41753032](https://pubmed.ncbi.nlm.nih.gov/41753032/)

---

## 14. 삶의 질 척도와 임상시험 종점 (Patient-reported endpoints)

227. Spertus JA, et al. **Minimally Important Kansas City Cardiomyopathy Questionnaire Changes Across the Spectrum of Heart Failure Severity.** JACC Heart Fail. 2025. [PMID 40908080](https://pubmed.ncbi.nlm.nih.gov/40908080/) — 모델의 `MCID_SAQ = 10 U` 설정 근거로 참고했습니다. **주의: 이 논문은 KCCQ에 대한 것이며 SAQ가 아닙니다.** SAQ 요약점수의 MCID는 이 세션에서 직접 확인하지 못했으므로 모델은 CorMicA가 달성한 +11.7 U를 임상적 중요성의 실용적 기준으로 삼고 10 U를 임계값으로 씁니다. 이 값은 모델의 결론(하위군 효과가 임계값 미만) 방향에는 영향을 주지만, +2.3 U라는 산출값 자체는 임계값 선택과 무관합니다.
228. Cao H, et al. **Use of comparative effectiveness research for similar Chinese patent medicine for angina pectoris: a new approach based on patient-important outcomes.** Trials. 2014. [PMID 24641790](https://pubmed.ncbi.nlm.nih.gov/24641790/)
229. Picano E, et al. **The clinical use of stress echocardiography in chronic coronary syndromes and beyond coronary artery disease: EACVI clinical consensus statement.** Eur Heart J Cardiovasc Imaging. 2024. [PMID 37798126](https://pubmed.ncbi.nlm.nih.gov/37798126/)

---

## 15. 모델의 반증 가능한 예측 (Falsifiable predictions, and where to look)

모델이 문헌에서 확인하지 못한 채 내놓은 예측들입니다. 각각을 검증하거나
반증할 수 있는 측정을 함께 적었습니다.

| # | 예측 | 검증 방법 | 근거 문헌 |
|---|------|-----------|-----------|
| P1 | 기능적 endotype의 **안정 시 심근 산소추출률은 정상의 약 56%**로 낮고, 관상정맥 산소포화도가 그만큼 높다 | 관상정맥동 채혈 1회 | §5(80, 95)로부터의 순수 연역. 안정 유량 상승 자체는 20에서 측정됨 |
| P2 | 최대충혈 미세혈관저항 ≥2.5를 "구조적"으로 읽는 것은 **리모델링과 아데노신 내성 수축긴장을 구별하지 못한다**. ROCK 억제 후 재측정하면 일부가 기준 아래로 내려간다 | 급성 Rho-kinase 또는 ETA 차단 후 최대충혈 측정 반복 | 121, 122, 124(파수딜이 남은 수축 성분을 푼다) |
| P3 | 이바브라딘/β차단의 이익은 **가장 낮은 CFR이 아니라 가장 높은 안정 심박수**를 가진 환자에서 가장 크다 | 계층화 무작위배정 또는 사후 상호작용 검정 | §III, 184–186 |
| P4 | 아미노필린은 **기능적 endotype에만** 유효하고 구조적 endotype에서는 A2A 차단으로 해로울 수 있다 | endotype 계층화 교차시험 | 135, 136, 138 |
| P5 | 지보텐탄의 음성 결과는 표적 실패가 아니라 **체액 보유**다. 이뇨제 병용 또는 저용량에서 운동시간 신호가 되살아난다 | PRIZE 설계에 이뇨제 병용군 추가 | 178 |
| P6 | 치료 반응 종점으로 CFR 대신 **최대충혈 미세혈관저항**을 미리 정해야 한다. 혈압을 낮추는 약은 CFR을 낮추지만 최대충혈 저항은 낮추지 않는다 | 기존 시험 자료의 재분석(두 지표를 함께 보고한 시험) | 43, 178 |

---

## 16. 방법론 · 도구 (Methods & tooling)

- mrgsolve: <https://mrgsolve.org/>
- Graphviz: <https://graphviz.org/>
- PubMed E-utilities: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- 이 디렉토리의 `cmd_reference_model.py`는 R 런타임 없이 모든 방정식을 먼저
  실행하기 위한 의존성 없는 Python 참조 구현이며, 수치 작업에서 드러난 **25개의
  실제 결함(B1–B25)** 을 파일 상단 BUG LOG에 남기고 각 결함이 있던 줄에 주석을
  달아두었습니다. `cmd_reference_output.txt`가 모든 산출 수치의 원본입니다.

---

*이 문서의 모든 PMID는 PubMed E-utilities로 직접 조회해 제목·저널·연도를
확인했습니다. 정량적 정박점(20, 65, 170, 174, 178, 182)은 초록 전문을 조회해
수치를 확인했습니다. 교육·연구 목적의 모델이며 임상 의사결정에 사용할 수
없습니다.*
