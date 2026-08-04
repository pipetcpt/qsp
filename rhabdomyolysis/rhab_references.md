# 횡문근분해증 유발 급성 신손상 — 참고문헌
# Rhabdomyolysis-Induced Acute Kidney Injury — References

본 QSP 모델의 구조와 파라미터 근거가 되는 문헌 80편입니다. 모든 PubMed 식별자(PMID)는 NCBI E-utilities로 조회하여 제목·학술지·연도를 확인한 것이며, 추정하거나 임의로 생성한 식별자는 없습니다.

All 80 PubMed identifiers below were resolved through the NCBI E-utilities API and their titles, journals and years verified against the returned records; none were reconstructed from memory. Where a citation supports a specific parameter or structural choice in `rhab_mrgsolve_model.R`, the connection is stated in the section preamble rather than claimed for individual papers.

---

## A. 정의 · 역학 · 임상 스펙트럼 · Definitions, epidemiology and clinical spectrum

*How rhabdomyolysis is defined, how often it causes AKI, and the historical origin of crush syndrome.*

1. **COVID-19-Associated Rhabdomyolysis in Hospitalized Patients: Retrospective Cohort Study.** Saudi Med J · 2026 · Esquerdo-Serrano P. [PMID 42238209](https://pubmed.ncbi.nlm.nih.gov/42238209/)
2. **Rhabdomyolysis and acute kidney injury.** N Engl J Med · 2009 · Bosch X. [PMID 19571284](https://pubmed.ncbi.nlm.nih.gov/19571284/)
3. **A risk prediction score for kidney failure or mortality in rhabdomyolysis.** JAMA Intern Med · 2013 · McMahon GM. [PMID 24000014](https://pubmed.ncbi.nlm.nih.gov/24000014/)
4. **Beyond muscle destruction: a systematic review of rhabdomyolysis for clinical practice.** Crit Care · 2016 · Chavez LO. [PMID 27301374](https://pubmed.ncbi.nlm.nih.gov/27301374/)
5. **Crush injuries with impairment of renal function. 1941.** J Am Soc Nephrol · 1998 · Bywaters EG. [PMID 9527411](https://pubmed.ncbi.nlm.nih.gov/9527411/)
6. **Rhabdomyolysis.** 2026 · Rout P. [PMID 28846335](https://pubmed.ncbi.nlm.nih.gov/28846335/)

## B. 미오글로빈 신독성 기전 · Mechanisms of myoglobin nephrotoxicity

*Heme iron, ferrihaemate chemistry, oxidative injury and the megalin/cubilin route into the proximal tubule — the model's ARM A.*

7. **L-Carnitine ameliorates glycerol-induced myoglobinuric acute renal failure in rats.** Ren Fail · 2009 · Ustundag S. [PMID 19212909](https://pubmed.ncbi.nlm.nih.gov/19212909/)
8. **Myoglobin inhibits proliferation of cultured human proximal tubular (HK-2) cells.** Kidney Int · 1996 · Iwata M. [PMID 8872953](https://pubmed.ncbi.nlm.nih.gov/8872953/)
9. **Hemoglobin and myoglobin associated oxidative stress: from molecular mechanisms to disease States.** Curr Med Chem · 2005 · Reeder BJ. [PMID 16305469](https://pubmed.ncbi.nlm.nih.gov/16305469/)
10. **The effect of hyperbaric oxygen therapy on rhabdomyolysis-induced myoglobinuric acute renal failure in rats.** Ren Fail · 2016 · Cebi G. [PMID 27765004](https://pubmed.ncbi.nlm.nih.gov/27765004/)
11. **Renal uptake of myoglobin is mediated by the endocytic receptors megalin and cubilin.** Am J Physiol Renal Physiol · 2003 · Gburek J. [PMID 12724130](https://pubmed.ncbi.nlm.nih.gov/12724130/)
12. **Regulation and immunohistochemical analysis of stress protein heme oxygenase-1 in rat kidney with myoglobinuric acute renal failure.** Biochem Biophys Res Commun · 1997 · Ishizuka S. [PMID 9367889](https://pubmed.ncbi.nlm.nih.gov/9367889/)
13. **Hemoglobin- and myoglobin-induced acute renal failure in rats: role of iron in nephrotoxicity.** Am J Physiol · 1988 · Paller MS. [PMID 3414810](https://pubmed.ncbi.nlm.nih.gov/3414810/)
14. **Redox reactions of myoglobin.** Antioxid Redox Signal · 2013 · Richards MP. [PMID 22900975](https://pubmed.ncbi.nlm.nih.gov/22900975/)
15. **Iron, heme oxygenase, and glutathione: effects on myohemoglobinuric proximal tubular injury.** Kidney Int · 1995 · Zager RA. [PMID 8544424](https://pubmed.ncbi.nlm.nih.gov/8544424/)

## C. 세관내 캐스트 형성 · Intratubular cast formation

*Tamm-Horsfall/uromodulin co-aggregation and the pH- and concentration-dependence behind the model's squared cast term — ARM B.*

16. **Mechanisms of intranephronal proteinaceous cast formation by low molecular weight proteins.** J Clin Invest · 1990 · Sanders PW. [PMID 2298921](https://pubmed.ncbi.nlm.nih.gov/2298921/)
17. **Myoglobin cast nephropathy following multiple bee stings.** Indian J Pathol Microbiol · 2023 · Madireddy N. [PMID 36656236](https://pubmed.ncbi.nlm.nih.gov/36656236/)
18. **Calcium oxalate monohydrate aggregation induced by aggregation of desialylated Tamm-Horsfall protein.** Urol Res · 2011 · Viswanathan P. [PMID 21229239](https://pubmed.ncbi.nlm.nih.gov/21229239/)

## D. 신장 혈관수축 · 허혈 · Renal vasoconstriction and ischaemia

*Nitric-oxide scavenging, endothelin and thromboxane, medullary hypoxia — ARM C.*

19. **Unmasking the Janus face of myoglobin in health and disease.** J Exp Biol · 2010 · Hendgen-Cotta UB. [PMID 20675542](https://pubmed.ncbi.nlm.nih.gov/20675542/)
20. **Renal cortical blood flow in glycerol-induced acute renal failure in the rat.** Circ Res · 1976 · Kurtz TW. [PMID 1244225](https://pubmed.ncbi.nlm.nih.gov/1244225/)

## E. 수액 요법 · 압사 증후군 · 재난 신장학 · Fluid resuscitation, crush syndrome and disaster nephrology

*The evidence that fluid TIMING dominates fluid DOSE, including the field-start protocol the model reproduces as scenario 6.*

21. **Early management of shock and prophylaxis of acute renal failure in traumatic rhabdomyolysis.** N Engl J Med · 1990 · Better OS. [PMID 2407958](https://pubmed.ncbi.nlm.nih.gov/2407958/)
22. **A Retrospctive Study of 377 Patients Admitted as an Emergency with Crush Syndrome After the Türkiye-Syria Earthquakes.** Med Sci Monit · 2024 · Görmeli Kurt N. [PMID 39127884](https://pubmed.ncbi.nlm.nih.gov/39127884/)
23. **Pediatric kidney care experience after the 2023 Türkiye earthquake.** Nephrol Dial Transplant · 2024 · Bakkaloğlu SA. [PMID 38327222](https://pubmed.ncbi.nlm.nih.gov/38327222/)
24. **Management of crush victims in mass disasters: highlights from recently published recommendations.** Clin J Am Soc Nephrol · 2013 · Sever MS. [PMID 23024157](https://pubmed.ncbi.nlm.nih.gov/23024157/)
25. **Rhabdomyolysis.** Acta Clin Belg · 2007 · Sever MS. [PMID 18284003](https://pubmed.ncbi.nlm.nih.gov/18284003/)
26. **The clinical features and outcome of crush patients with acute kidney injury after the Wenchuan earthquake: differences between elderly and younger adults.** Injury · 2012 · Zhang L. [PMID 21144512](https://pubmed.ncbi.nlm.nih.gov/21144512/)
27. **Balanced Crystalloid Solutions.** Am J Respir Crit Care Med · 2019 · Semler MW. [PMID 30407838](https://pubmed.ncbi.nlm.nih.gov/30407838/)
28. **Acute exertional rhabdomyolysis.** Am Fam Physician · 1995 · Line RL. [PMID 7625324](https://pubmed.ncbi.nlm.nih.gov/7625324/)
29. **Disaster nephrology: crush injury and beyond.** Kidney Int · 2014 · Gibney RT. [PMID 24107850](https://pubmed.ncbi.nlm.nih.gov/24107850/)

## F. 중탄산염 · 만니톨 · Bicarbonate and mannitol

*The negative trials that the model's product structure reconciles with the positive mechanism, plus osmotic nephrosis.*

30. **Preventing renal failure in patients with rhabdomyolysis: do bicarbonate and mannitol make a difference?.** J Trauma · 2004 · Brown CV. [PMID 15211124](https://pubmed.ncbi.nlm.nih.gov/15211124/)
31. **Prophylaxis of acute renal failure in patients with rhabdomyolysis.** Ren Fail · 1997 · Homsi E. [PMID 9101605](https://pubmed.ncbi.nlm.nih.gov/9101605/)
32. **[Rhabdomyolysis].** Acta Med Port · 2005 · Rosa NG. [PMID 16584660](https://pubmed.ncbi.nlm.nih.gov/16584660/)
33. **Mannitol-induced acute renal failure.** Neth J Med · 1997 · van Hengel P. [PMID 9038039](https://pubmed.ncbi.nlm.nih.gov/9038039/)
34. **Osmotic nephrosis with mannitol: review article.** Ren Fail · 2014 · Nomani AZ. [PMID 24941319](https://pubmed.ncbi.nlm.nih.gov/24941319/)
35. **The Role of Urine Alkalinization in Preventing Rhabdomyolysis-Induced Acute Kidney Injury and Need for Dialysis: A Systematic Review and Meta-Analysis.** Bull Emerg Trauma · 2025. [PMID 41268470](https://pubmed.ncbi.nlm.nih.gov/41268470/)

## G. 정질액 선택 · 염소 부하 · Crystalloid choice and chloride load

*Why saline buys the flow term and spends the pH term.*

36. **Balanced Crystalloids versus Saline in Critically Ill Adults.** N Engl J Med · 2018 · Semler MW. [PMID 29485925](https://pubmed.ncbi.nlm.nih.gov/29485925/)
37. **[Hyperchloremic acidosis druing plasma volume replacement].** Ann Fr Anesth Reanim · 2002 · Blanloeil Y. [PMID 11963385](https://pubmed.ncbi.nlm.nih.gov/11963385/)
38. **A Crossover Trial of Hospital-Wide Lactated Ringer's Solution versus Normal Saline.** N Engl J Med · 2025 · McIntyre L. [PMID 40503714](https://pubmed.ncbi.nlm.nih.gov/40503714/)
39. **A prospective, randomized, comparison study on effect of perioperative use of chloride liberal intravenous fluids versus chloride restricted intravenous fluids on postoperative acute kidney injury in patients undergoing off-pump coronary artery bypass grafting surgeries.** Ann Card Anaesth · 2018 · Bhaskaran K. [PMID 30333337](https://pubmed.ncbi.nlm.nih.gov/30333337/)

## H. 표지자 동역학 — 두 개의 시계 · Marker kinetics — the two clocks

*Creatine kinase and myoglobin half-lives, renal clearance of myoglobin, and why the dipstick is often negative.*

40. **Prognostic value, kinetics and effect of CVVHDF on serum of the myoglobin and creatine kinase in critically ill patients with rhabdomyolysis.** Acta Anaesthesiol Scand · 2005 · Mikkelsen TS. [PMID 15954972](https://pubmed.ncbi.nlm.nih.gov/15954972/)
41. **Rapid renal clearance of immunoreactive canine plasma myoglobin.** Circulation · 1982 · Klocke FJ. [PMID 7074811](https://pubmed.ncbi.nlm.nih.gov/7074811/)
42. **Investigating skeletal muscle biomarkers for the early detection of Australian myotoxic snake envenoming: an animal model pilot study.** Clin Toxicol (Phila) · 2024 · Johnston CI. [PMID 38804832](https://pubmed.ncbi.nlm.nih.gov/38804832/)
43. **Mild rhabdomyolysis in a child with fever and "hematuria".** Pediatr Nephrol · 2003 · Tasic V. [PMID 12736809](https://pubmed.ncbi.nlm.nih.gov/12736809/)
44. **Exertional rhabdomyolysis: a case series of 30 hospitalized patients.** Mil Med · 2015 · Oh RC. [PMID 25643388](https://pubmed.ncbi.nlm.nih.gov/25643388/)
45. **Rhabdomyolysis.** Dis Mon · 2020 · Cabral BMI. [PMID 32532456](https://pubmed.ncbi.nlm.nih.gov/32532456/)

## I. 전해질 · 산염기 · Electrolytes and acid-base

*Potassium, the biphasic calcium-phosphate course, urate and the rebound hypercalcaemia of the recovery phase.*

46. **[Crush syndrome].** Minerva Chir · 2007 · Scapellato S. [PMID 17641588](https://pubmed.ncbi.nlm.nih.gov/17641588/)
47. **Incidence of Electrolyte Imbalances Following Traumatic Rhabdomyolysis: A Systematic Review and Meta-Analysis.** Cureus · 2024 · Safari S. [PMID 38817473](https://pubmed.ncbi.nlm.nih.gov/38817473/)
48. **Rhabdomyolysis-Induced Resistant Hypercalcemia During the Recovery Phase of Acute Kidney Injury.** Cureus · 2026 · Das A. [PMID 42306390](https://pubmed.ncbi.nlm.nih.gov/42306390/)
49. **Hyperuricemic nephropathies.** Nephron · 1999 · Steele TH. [PMID 9873214](https://pubmed.ncbi.nlm.nih.gov/9873214/)
50. **Potassium Disorders: Hypokalemia and Hyperkalemia.** Am Fam Physician · 2015 · Viera AJ. [PMID 26371733](https://pubmed.ncbi.nlm.nih.gov/26371733/)
51. **Acute hyperkalemia in the emergency department: a summary from a Kidney Disease: Improving Global Outcomes conference.** Eur J Emerg Med · 2020 · Lindner G. [PMID 32852924](https://pubmed.ncbi.nlm.nih.gov/32852924/)

## J. 구획 증후군 · 근막절개 · Compartment syndrome and fasciotomy

*The delta-pressure threshold and the limb-versus-kidney trade-off of reperfusion.*

52. **Do one-time intracompartmental pressure measurements have a high false-positive rate in diagnosing compartment syndrome?.** J Trauma Acute Care Surg · 2014 · Whitney A. [PMID 24458053](https://pubmed.ncbi.nlm.nih.gov/24458053/)
53. **Long-term physical outcome of patients who suffered crush syndrome after the 1995 Hanshin-Awaji earthquake: prognostic indicators in retrospect.** J Trauma · 2002 · Matsuoka T. [PMID 11791049](https://pubmed.ncbi.nlm.nih.gov/11791049/)
54. **Acute compartment syndrome.** Injury · 2017 · Schmidt AH. [PMID 28449851](https://pubmed.ncbi.nlm.nih.gov/28449851/)

## K. 체외 순환 요법 · Extracorporeal therapy

*High cut-off membranes, myoglobin clearance and the arithmetic ceiling on removing a molecule the kidney clears faster.*

55. **High cut-off renal replacement therapy for removal of myoglobin in severe rhabdomyolysis and acute kidney injury: a case series.** Nephron Clin Pract · 2012 · Heyne N. [PMID 23327834](https://pubmed.ncbi.nlm.nih.gov/23327834/)
56. **Myoglobin clearance by continuous venous-venous haemofiltration in rhabdomyolysis with acute kidney injury: a case series.** Injury · 2012 · Zhang L. [PMID 20843513](https://pubmed.ncbi.nlm.nih.gov/20843513/)
57. **Factors predicting kidney replacement therapy in pediatric earthquake victims with crush syndrome in the first week following rescue.** Eur J Pediatr · 2023 · Atmis B. [PMID 37804325](https://pubmed.ncbi.nlm.nih.gov/37804325/)
58. **Rhabdomyolysis and acute renal failure in severely burned patients.** Burns · 2011 · Stollwerck PL. [PMID 20965664](https://pubmed.ncbi.nlm.nih.gov/20965664/)

## L. 원인별 병인 — 약물 · 운동 · 유전 · Cause-specific aetiology — drugs, exertion, genetics

*Statin myopathy, exertional rhabdomyolysis, malignant hyperthermia and dantrolene.*

59. **Statin-Associated Myopathy: Emphasis on Mechanisms and Targeted Therapy.** Int J Mol Sci · 2021 · Vinci P. [PMID 34769118](https://pubmed.ncbi.nlm.nih.gov/34769118/)
60. **Cerivastatin and gemfibrozil-associated rhabdomyolysis.** Ann Pharmacother · 2001 · Bruno-Joyce J. [PMID 11573847](https://pubmed.ncbi.nlm.nih.gov/11573847/)
61. **Twelve cases of exertional rhabdomyolysis in college football players from the same institution over a 23-year span: a descriptive study.** Phys Sportsmed · 2018 · Thompson TL. [PMID 29855209](https://pubmed.ncbi.nlm.nih.gov/29855209/)
62. **Dantrolene requires Mg(2+) to arrest malignant hyperthermia.** Proc Natl Acad Sci U S A · 2017 · Choi RH. [PMID 28373535](https://pubmed.ncbi.nlm.nih.gov/28373535/)
63. **Pharmacokinetics of Dantrolene in the Plasma Exchange Treatment of Malignant Hyperthermia in a 14-Year-Old Chinese Boy: A Case Report and Literature Review.** Front Med (Lausanne) · 2022 · Li X. [PMID 36035384](https://pubmed.ncbi.nlm.nih.gov/36035384/)
64. **Genome-Wide Analysis of Exertional Rhabdomyolysis in Sickle Cell Trait Positive African Americans.** Genes (Basel) · 2024 · Ren M. [PMID 38674343](https://pubmed.ncbi.nlm.nih.gov/38674343/)
65. **Drug induced rhabdomyolysis.** Curr Opin Pharmacol · 2012 · Hohenegger M. [PMID 22560920](https://pubmed.ncbi.nlm.nih.gov/22560920/)
66. **Delayed Presentation of Acute Gluteal Compartment Syndrome.** Am J Case Rep · 2016 · Tasch JJ. [PMID 27432320](https://pubmed.ncbi.nlm.nih.gov/27432320/)

## M. 예측 점수 · 결과 · 회복 · Prediction scores, outcomes and recovery

*The McMahon score, KDIGO staging, and the unusually good long-term renal recovery of rhabdomyolysis-AKI.*

67. **Rhabdomyolysis and acute kidney injury: creatine kinase as a prognostic marker and validation of the McMahon Score in a 10-year cohort: A retrospective observational evaluation.** Eur J Anaesthesiol · 2016 · Simpson JP. [PMID 27259093](https://pubmed.ncbi.nlm.nih.gov/27259093/)
68. **The Role of Scoring Systems and Urine Dipstick in Prediction of Rhabdomyolysis-induced Acute Kidney Injury: a Systematic Review.** Iran J Kidney Dis · 2016 · Safari S. [PMID 27225716](https://pubmed.ncbi.nlm.nih.gov/27225716/)
69. **Specific macrophage subtypes influence the progression of rhabdomyolysis-induced kidney injury.** J Am Soc Nephrol · 2015 · Belliere J. [PMID 25270069](https://pubmed.ncbi.nlm.nih.gov/25270069/)
70. **Diagnosis, evaluation, and management of acute kidney injury: a KDIGO summary (Part 1).** Crit Care · 2013 · Kellum JA. [PMID 23394211](https://pubmed.ncbi.nlm.nih.gov/23394211/)
71. **Development and validation of an interpretable multi-task model to predict outcomes in patients with rhabdomyolysis: a multicenter retrospective cohort study.** EClinicalMedicine · 2025 · Liu C. [PMID 40896465](https://pubmed.ncbi.nlm.nih.gov/40896465/)
72. **Early start of hemoadsorption is associated with improved kidney recovery in ICU patients with ischemic/reperfusion cause of rhabdomyolysis.** PLoS One · 2026 · Premuzic V. [PMID 42378236](https://pubmed.ncbi.nlm.nih.gov/42378236/)

## N. 정량 모델링 방법론 · Quantitative modelling methodology

*mrgsolve, systems models of tubular injury, potassium and volume homeostasis, and the crystallisation chemistry used here.*

73. **Multiscale Mathematical Model of Drug-Induced Proximal Tubule Injury: Linking Urinary Biomarkers to Epithelial Cell Injury and Renal Dysfunction.** Toxicol Sci · 2018 · Gebremichael Y. [PMID 29126144](https://pubmed.ncbi.nlm.nih.gov/29126144/)
74. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** CPT Pharmacometrics Syst Pharmacol · 2019 · Elmokadem A. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
75. **Dialysis and Acid-Base Balance: A Comparative Physiological Analysis of Boston and Stewart Models.** J Clin Med · 2025 · Kroustalakis N. [PMID 41303241](https://pubmed.ncbi.nlm.nih.gov/41303241/)
76. **A mathematical model of potassium homeostasis: Effect of feedforward and feedback controls.** PLoS Comput Biol · 2022 · Stadt MM. [PMID 36538563](https://pubmed.ncbi.nlm.nih.gov/36538563/)
77. **Role of arginine vasopressin in the regulation of extracellular fluid volume.** Med Sci Sports Exerc · 1996 · Norsk P. [PMID 8897402](https://pubmed.ncbi.nlm.nih.gov/8897402/)
78. **Urine alkalinization for dissolution of uric acid stones and treatment of other urological diseases with a treatment combining potassium magnesium citrate and theobromine.** Arch Ital Urol Androl · 2025 · Rodriguez-Hesles CA. [PMID 40162813](https://pubmed.ncbi.nlm.nih.gov/40162813/)
79. **[Tumor lysis syndrome: risk factors and treatment].** Wien Klin Wochenschr · 2005 · Hörl WH. [PMID 15986584](https://pubmed.ncbi.nlm.nih.gov/15986584/)
80. **Modeling the transfer of low molecular weight plasma proteins during hemodialysis and online hemodiafiltration.** Artif Organs · 2021 · Glancey G. [PMID 33001450](https://pubmed.ncbi.nlm.nih.gov/33001450/)

---

## 인용 원칙 (How these references are used)

이 모델은 위 문헌들의 **중심 경향치**를 이용해 파라미터를 설정했으며, 특정 환자 데이터셋에 적합(fitting)한 것이 아닙니다. 모델이 재현하도록 보정된 임상 정박점(anchor)은 다음과 같습니다.

| 정박점 (anchor) | 근거 섹션 | 모델 출력 |
|---|---|---|
| 중증 압사 CK 최고치 수만 U/L, 24–72시간에 정점 | A, H | 63,900 U/L @ 38–39 h |
| 미오글로빈 반감기 약 2–3시간, 조기 정점 | H | 정점 11.8 h, t½ 약 2.4 h |
| CK 반감기 약 36–48시간 | H | 36 h (설정값), 정점 지연 25.7 h |
| 현장에서 수액을 시작하면 신부전이 거의 발생하지 않음 | E | KDIGO 1, eGFR d90 113.8 |
| 24시간 지연 시 투석이 필요한 AKI | E | KDIGO 3, 무뇨 14 h, McMahon 11 |
| 충분한 수액 위에 추가한 중탄산염의 이득은 검출되지 않음 | F | 절대 이득 43 → 18 단위 (100 → 1000 mL/h) |
| 만니톨 100 g/일은 중립, 누적 고용량은 신독성 | F | 크레아티닌 2.17 vs 2.15 / 5.18 |
| 미치료 압사에서 K⁺ 7–9 mmol/L | I | 8.07 mmol/L |
| 회복기 반동성 고칼슘혈증 | I | 이온화 Ca 0.88 → 1.41 mmol/L |
| 구획압 ΔP < 30 mmHg가 근막절개 기준 | J | ΔP 게이트 중심값 30 mmHg |
| 고차단막 CVVH의 미오글로빈 제거 효과는 제한적 | K | 총 이득의 약 1/10만이 Mb 제거 |
| 횡문근분해증 AKI의 신기능 회복은 상대적으로 양호 | M | eGFR d90 84.8–113.8 |

정박점을 맞추지 못했거나 의도적으로 과장한 시나리오(만니톨 과다투여, 근막절개 미시행 등)는 `rhab_mrgsolve_model.R` 하단의 CALIBRATION NOTES와 `README.md`의 한계 항목에 명시했습니다.
