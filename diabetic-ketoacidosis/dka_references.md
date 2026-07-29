# 참고문헌 — 당뇨병성 케톤산증 / 고혈당 고삼투압 상태 QSP 모델
# References — Diabetic Ketoacidosis & Hyperglycaemic Hyperosmolar State QSP Model

각 항목은 모델의 어느 부분을 근거하는지 함께 표시했습니다.
Each entry notes which part of the model it supports. 102 references.

---

## 1. 정의·역학·진료지침 (Definitions, epidemiology, guidelines)

1. Kitabchi AE, Umpierrez GE, Miles JM, Fisher JN. **Hyperglycemic crises in adult patients with diabetes.** Diabetes Care 2009;32:1335-43. — the source of the deficit table (water 100 mL/kg, Na 7-10, Cl 3-5, K 3-5 mmol/kg) that the model's 24 h lead-in is calibrated against. https://pubmed.ncbi.nlm.nih.gov/19564476/
2. Umpierrez GE, Davis GM, ElSayed NA, et al. **Hyperglycaemic crises in adults with diabetes: a consensus report.** Diabetologia 2024;67:1455-79. — current diagnostic thresholds, including the shift to β-hydroxybutyrate ≥3.0 mmol/L and the recognition of euglycaemic presentations. https://pubmed.ncbi.nlm.nih.gov/38907161/
3. Dhatariya KK, et al (Joint British Diabetes Societies). **The management of diabetic ketoacidosis in adults.** 2021 (JBDS 02). — the ketone-driven protocol, target BHB fall >0.5 mmol/L/h, fixed-rate 0.1 U/kg/h. Underlies the model's resolution criteria.
4. Dhatariya KK, Glaser NS, Codner E, Umpierrez GE. **Diabetic ketoacidosis.** Nat Rev Dis Primers 2020;6:40. https://pubmed.ncbi.nlm.nih.gov/32409703/
5. Wolfsdorf JI, Glaser N, Agus M, et al. **ISPAD Clinical Practice Consensus Guidelines 2018: Diabetic ketoacidosis and the hyperglycemic hyperosmolar state.** Pediatr Diabetes 2018;19(Suppl 27):155-77. — paediatric fluid rates used in the fluid-rate arms. https://pubmed.ncbi.nlm.nih.gov/29900641/
6. Glaser N, Fritsch M, Priyambada L, et al. **ISPAD Clinical Practice Consensus Guidelines 2022: Diabetic ketoacidosis and hyperglycemic hyperosmolar state.** Pediatr Diabetes 2022;23:835-56. https://pubmed.ncbi.nlm.nih.gov/36250645/
7. Benoit SR, Zhang Y, Geiss LS, Gregg EW, Albright A. **Trends in diabetic ketoacidosis hospitalizations and in-hospital mortality — United States, 2000-2014.** MMWR 2018;67:362-5. https://pubmed.ncbi.nlm.nih.gov/29596401/
8. Desai D, Mehta D, Mathias P, Menon G, Schubart UK. **Health care utilization and burden of diabetic ketoacidosis in the U.S. over the past decade.** Diabetes Care 2018;41:1631-8. https://pubmed.ncbi.nlm.nih.gov/29773640/
9. Pasquel FJ, Umpierrez GE. **Hyperosmolar hyperglycemic state: a historic review of the clinical presentation, diagnosis, and treatment.** Diabetes Care 2014;37:3124-31. — the HHS phenotype the model reproduces with two parameter changes. https://pubmed.ncbi.nlm.nih.gov/25342831/
10. Fayfman M, Pasquel FJ, Umpierrez GE. **Management of hyperglycemic crises: diabetic ketoacidosis and hyperglycemic hyperosmolar state.** Med Clin North Am 2017;101:587-606. https://pubmed.ncbi.nlm.nih.gov/28372715/

---

## 2. 인슐린 용량-반응: 모델의 핵심 비대칭 (Insulin dose–response — the four IC50s the model turns on)

11. Zierler KL, Rabinowitz D. **Effect of very small concentrations of insulin on forearm metabolism.** J Clin Invest 1964;43:950-62. — the original demonstration that antilipolysis occurs at insulin concentrations far below those needed for glucose uptake. The model's IC50_LIP = 15 vs EC50_UP = 60 µU/mL. https://pubmed.ncbi.nlm.nih.gov/14169520/
12. Rizza RA, Mandarino LJ, Gerich JE. **Dose-response characteristics for effects of insulin on production and utilisation of glucose in man.** Am J Physiol 1981;240:E630-9. — half-maximal suppression of hepatic glucose output near 30 µU/mL, half-maximal stimulation of utilisation near 60. https://pubmed.ncbi.nlm.nih.gov/7018429/
13. Nurjhan N, Campbell PJ, Kennedy FP, Miles JM, Gerich JE. **Insulin dose-response characteristics for suppression of glycerol release and conversion to glucose in humans.** Diabetes 1986;35:1326-31. https://pubmed.ncbi.nlm.nih.gov/3530852/
14. Jensen MD, Caruso M, Heiling V, Miles JM. **Insulin regulation of lipolysis in nondiabetic and IDDM subjects.** Diabetes 1989;38:1595-601. https://pubmed.ncbi.nlm.nih.gov/2573554/
15. Groop LC, Bonadonna RC, DelPrato S, et al. **Glucose and free fatty acid metabolism in non-insulin-dependent diabetes mellitus: evidence for multiple sites of insulin resistance.** J Clin Invest 1989;84:205-13. https://pubmed.ncbi.nlm.nih.gov/2661589/
16. Kitabchi AE, Ayyagari V, Guerra SM. **The efficacy of low-dose versus conventional therapy of insulin for treatment of diabetic ketoacidosis.** Ann Intern Med 1976;84:633-8. — the trial the model *derives* rather than assumes: an 8-fold reduction in insulin dose did not slow recovery. https://pubmed.ncbi.nlm.nih.gov/818050/
17. Alberti KG, Hockaday TD, Turner RC. **Small doses of intramuscular insulin in the treatment of diabetic "coma".** Lancet 1973;2:515-22. https://pubmed.ncbi.nlm.nih.gov/4125959/
18. Padilla AJ, Loeb JN. **"Low-dose" versus "high-dose" insulin regimens in the management of uncontrolled diabetes.** Am J Med 1977;63:843-8. https://pubmed.ncbi.nlm.nih.gov/22195/
19. Ferrannini E, Barrett EJ, Bevilacqua S, DeFronzo RA. **Effect of fatty acids on glucose production and utilization in man.** J Clin Invest 1983;72:1737-47. — the Randle-type lipid resistance the model uses for the falling insulin requirement during treatment. https://pubmed.ncbi.nlm.nih.gov/6358256/
20. Boden G, Chen X, Ruiz J, White JV, Rossetti L. **Mechanisms of fatty acid-induced inhibition of glucose uptake.** J Clin Invest 1994;93:2438-46. https://pubmed.ncbi.nlm.nih.gov/8200979/

---

## 3. 문맥 인슐린 우선노출 — HHS가 케톤산증이 아닌 이유 (Portal insulin privilege)

21. Field JB. **Extraction of insulin by liver.** Annu Rev Med 1973;24:309-14. — hepatic first-pass extraction of ~50%, the basis of the model's PORTF parameter. https://pubmed.ncbi.nlm.nih.gov/4267050/
22. Eaton RP, Allen RC, Schade DS. **Hepatic removal of insulin in normal man: dose response to endogenous insulin secretion.** J Clin Endocrinol Metab 1983;56:1294-300. https://pubmed.ncbi.nlm.nih.gov/6341391/
23. Meier JJ, Veldhuis JD, Butler PC. **Pulsatile insulin secretion dictates systemic insulin delivery by regulating hepatic insulin extraction in humans.** Diabetes 2005;54:1649-56. https://pubmed.ncbi.nlm.nih.gov/15919785/
24. Chap Z, Ishida T, Chou J, et al. **First-pass hepatic extraction and metabolic effects of insulin and insulin analogues.** Am J Physiol 1987;252:E209-17. https://pubmed.ncbi.nlm.nih.gov/3548518/
25. Edgerton DS, Lautz M, Scott M, et al. **Insulin's direct effects on the liver dominate the control of hepatic glucose production.** J Clin Invest 2006;116:521-7. https://pubmed.ncbi.nlm.nih.gov/16453026/

---

## 4. 케톤 생성: CPT-1 관문과 malonyl-CoA (Ketogenesis — the gate)

26. McGarry JD, Foster DW. **Regulation of hepatic fatty acid oxidation and ketone body production.** Annu Rev Biochem 1980;49:395-420. — the substrate-times-gate architecture the model implements literally. https://pubmed.ncbi.nlm.nih.gov/6157353/
27. McGarry JD, Mannaerts GP, Foster DW. **A possible role for malonyl-CoA in the regulation of hepatic fatty acid oxidation and ketogenesis.** J Clin Invest 1977;60:265-70. https://pubmed.ncbi.nlm.nih.gov/874089/
28. McGarry JD, Woeltje KF, Kuwajima M, Foster DW. **Regulation of ketogenesis and the renaissance of carnitine palmitoyltransferase.** Diabetes Metab Rev 1989;5:271-84. https://pubmed.ncbi.nlm.nih.gov/2656156/
29. Foster DW, McGarry JD. **The metabolic derangements and treatment of diabetic ketoacidosis.** N Engl J Med 1983;309:159-69. https://pubmed.ncbi.nlm.nih.gov/6408658/
30. Hegardt FG. **Mitochondrial 3-hydroxy-3-methylglutaryl-CoA synthase: a control enzyme in ketogenesis.** Biochem J 1999;338:569-82. https://pubmed.ncbi.nlm.nih.gov/10051425/
31. Puchalska P, Crawford PA. **Multi-dimensional roles of ketone bodies in fuel metabolism, signaling, and therapeutics.** Cell Metab 2017;25:262-84. — the definitive review of ketone body flux, including the SCOT/OXCT1 asymmetry that makes the liver unable to oxidise its own product. https://pubmed.ncbi.nlm.nih.gov/28178565/
32. Cotter DG, Schugar RC, Crawford PA. **Ketone body metabolism and cardiovascular disease.** Am J Physiol Heart Circ Physiol 2013;304:H1060-76. https://pubmed.ncbi.nlm.nih.gov/23396451/
33. Miles JM, Haymond MW, Nissen SL, Gerich JE. **Effects of free fatty acid availability, glucagon excess, and insulin deficiency on ketone body production in postabsorptive man.** J Clin Invest 1983;71:1554-61. — the experiment that separates the substrate arm from the gate. https://pubmed.ncbi.nlm.nih.gov/6863539/
34. Keller U, Lustenberger M, Müller-Brand J, Gerber PP, Stauffacher W. **Human ketone body production and utilization studied using tracer techniques: regulation by free fatty acids, insulin, catecholamines, and thyroid hormones.** Diabetes Metab Rev 1989;5:285-98. https://pubmed.ncbi.nlm.nih.gov/2656157/
35. Balasse EO, Féry F. **Ketone body production and disposal: effects of fasting, diabetes, and exercise.** Diabetes Metab Rev 1989;5:247-70. — the saturable disposal that makes the ketone arm bounded, and the linear cerebral arm. https://pubmed.ncbi.nlm.nih.gov/2656155/
36. Owen OE, Morgan AP, Kemp HG, Sullivan JM, Herrera MG, Cahill GF. **Brain metabolism during fasting.** J Clin Invest 1967;46:1589-95. — concentration-driven (non-saturable) cerebral ketone uptake, the KLIN_KOX term. https://pubmed.ncbi.nlm.nih.gov/6061736/
37. Hall SE, Wastney ME, Bolton TM, Braaten JT, Berman M. **Ketone body kinetics in humans: the effects of insulin-dependent diabetes, obesity, and starvation.** J Lipid Res 1984;25:1184-94. https://pubmed.ncbi.nlm.nih.gov/6438334/
38. Fery F, Balasse EO. **Ketone body turnover during and after exercise in overnight-fasted and starved humans.** Am J Physiol 1983;245:E318-25. https://pubmed.ncbi.nlm.nih.gov/6412636/

---

## 5. BHB/AcAc 비와 산화환원 상태 (Hepatic redox and the nitroprusside paradox)

39. Williamson DH, Lund P, Krebs HA. **The redox state of free nicotinamide-adenine dinucleotide in the cytoplasm and mitochondria of rat liver.** Biochem J 1967;103:514-27. — the equilibrium that sets BHB/AcAc from NADH/NAD⁺. https://pubmed.ncbi.nlm.nih.gov/4291787/
40. Sherwin RS, Hendler RG, Felig P. **Effect of ketone infusions on amino acid and nitrogen metabolism in man.** J Clin Invest 1975;55:1382-90. https://pubmed.ncbi.nlm.nih.gov/1133182/
41. Umpierrez GE, Watts NB, Phillips LS. **Clinical utility of β-hydroxybutyrate determined by reflectance meter in the management of diabetic ketoacidosis.** Diabetes Care 1995;18:137-8. https://pubmed.ncbi.nlm.nih.gov/7538325/
42. Sheikh-Ali M, Karon BS, Basu A, et al. **Can serum β-hydroxybutyrate be used to diagnose diabetic ketoacidosis?** Diabetes Care 2008;31:643-7. https://pubmed.ncbi.nlm.nih.gov/18184896/
43. Wallace TM, Matthews DR. **Recent advances in the monitoring and management of diabetic ketoacidosis.** QJM 2004;97:773-80. — why the urine nitroprusside reaction can worsen while the patient improves. https://pubmed.ncbi.nlm.nih.gov/15569808/
44. Csako G. **False-positive results for ketone with the drug mesna and other free-sulfhydryl compounds.** Clin Chem 1987;33:289-92. https://pubmed.ncbi.nlm.nih.gov/3802509/
45. Reichard GA Jr, Skutches CL, Hoeldtke RD, Owen OE. **Acetone metabolism in humans during diabetic ketoacidosis.** Diabetes 1986;35:668-74. — acetone's long half-life, hence days of positive strips and ketotic breath after resolution. https://pubmed.ncbi.nlm.nih.gov/3086240/
46. Owen OE, Trapp VE, Skutches CL, et al. **Acetone metabolism during diabetic ketoacidosis.** Diabetes 1982;31:242-8. https://pubmed.ncbi.nlm.nih.gov/6802673/

---

## 6. 산-염기: Stewart / 강이온차 접근 (Acid–base: the strong-ion framework)

47. Stewart PA. **Modern quantitative acid-base chemistry.** Can J Physiol Pharmacol 1983;61:1444-61. — the framework in which bicarbonate is a dependent variable, as implemented here. https://pubmed.ncbi.nlm.nih.gov/6423247/
48. Figge J, Rossing TH, Fencl V. **The role of serum proteins in acid-base equilibria.** J Lab Clin Med 1991;117:453-67. — the weak-acid coefficients used verbatim in the pH solve. https://pubmed.ncbi.nlm.nih.gov/2045713/
49. Figge J, Mydosh T, Fencl V. **Serum proteins and acid-base equilibria: a follow-up.** J Lab Clin Med 1992;120:713-9. https://pubmed.ncbi.nlm.nih.gov/1431500/
50. Kellum JA. **Determinants of blood pH in health and disease.** Crit Care 2000;4:6-14. https://pubmed.ncbi.nlm.nih.gov/11094491/
51. Story DA, Morimatsu H, Bellomo R. **Strong ions, weak acids and base excess: a simplified Fencl-Stewart approach to clinical acid-base disorders.** Br J Anaesth 2004;92:54-60. https://pubmed.ncbi.nlm.nih.gov/14742325/
52. Morgan TJ. **The Stewart approach — one clinician's perspective.** Clin Biochem Rev 2009;30:41-54. https://pubmed.ncbi.nlm.nih.gov/19565024/
53. Albert MS, Dell RB, Winters RW. **Quantitative displacement of acid-base equilibrium in metabolic acidosis.** Ann Intern Med 1967;66:312-22. — the source of the respiratory-compensation relation, and of the observation that Winter's line overestimates compensation at very low bicarbonate. https://pubmed.ncbi.nlm.nih.gov/6016545/
54. Fulop M. **Ventilatory response in patients with acute severe metabolic acidosis.** Clin Sci 1976;50:367-73. — measured PCO₂ in DKA, which is why the model uses PCO₂ = 13 + 1.15·HCO₃ rather than Winter's formula. https://pubmed.ncbi.nlm.nih.gov/819308/
55. Adrogué HJ, Wilson H, Boyd AE, Suki WN, Eknoyan G. **Plasma acid-base patterns in diabetic ketoacidosis.** N Engl J Med 1982;307:1603-10. — the pattern the model reproduces, including the frequency of a superimposed hyperchloraemic component. https://pubmed.ncbi.nlm.nih.gov/6815530/

---

## 7. 치료 후 고염소성 산증 (Post-treatment hyperchloraemic acidosis)

56. Adrogué HJ, Eknoyan G, Suki WK. **Diabetic ketoacidosis: role of the kidney in the acid-base homeostasis re-evaluated.** Kidney Int 1984;25:591-8. — urinary loss of ketoanion with sodium as loss of *potential bicarbonate*: the mechanism the model separates from the chloride load. https://pubmed.ncbi.nlm.nih.gov/6431025/
57. Oh MS, Carroll HJ, Goldstein DA, Fein IA. **Hyperchloremic acidosis during the recovery phase of diabetic ketosis.** Ann Intern Med 1978;89:925-7. https://pubmed.ncbi.nlm.nih.gov/102229/
58. Oh MS, Banerji MA, Carroll HJ. **The mechanism of hyperchloremic acidosis during the recovery phase of diabetic ketoacidosis.** Diabetes 1981;30:310-3. https://pubmed.ncbi.nlm.nih.gov/6782024/
59. Halperin ML, Kamel KS. **D-lactic acidosis and the ketoacidoses: the importance of the urine.** Kidney Int 1996;49:1-8. https://pubmed.ncbi.nlm.nih.gov/8770942/
60. Taylor D, Durward A, Tibby SM, et al. **The influence of hyperchloraemia on acid base interpretation in diabetic ketoacidosis.** Intensive Care Med 2006;32:295-301. https://pubmed.ncbi.nlm.nih.gov/16450092/
61. Chua HR, Venkatesh B, Stachowski E, et al. **Plasma-Lyte 148 vs 0.9% saline for fluid resuscitation in diabetic ketoacidosis.** J Crit Care 2012;27:138-45. https://pubmed.ncbi.nlm.nih.gov/22440386/
62. Van Zyl DG, Rheeder P, Delport E. **Fluid management in diabetic-acidosis — Ringer's lactate versus normal saline: a randomized controlled trial.** QJM 2012;105:337-43. https://pubmed.ncbi.nlm.nih.gov/22109683/
63. Self WH, Evans CS, Jenkins CA, et al. **Clinical effects of balanced crystalloids vs saline in adults with diabetic ketoacidosis: a subgroup analysis of cluster randomized clinical trials (SMART/SALT-ED).** JAMA Netw Open 2020;3:e2024596. https://pubmed.ncbi.nlm.nih.gov/33196806/
64. Ramanan M, Attokaran A, Murray L, et al. **Sodium chloride or Plasmalyte-148 evaluation in severe diabetic ketoacidosis (SCOPE-DKA): a cluster, crossover, randomized, controlled trial.** Intensive Care Med 2021;47:1248-57. https://pubmed.ncbi.nlm.nih.gov/34609547/
65. Semler MW, Self WH, Wanderer JP, et al. **Balanced crystalloids versus saline in critically ill adults (SMART).** N Engl J Med 2018;378:829-39. https://pubmed.ncbi.nlm.nih.gov/29485925/
66. Yunos NM, Bellomo R, Hegarty C, et al. **Association between a chloride-liberal vs chloride-restrictive intravenous fluid administration strategy and kidney injury in critically ill adults.** JAMA 2012;308:1566-72. https://pubmed.ncbi.nlm.nih.gov/23073953/

---

## 8. 신장: 삼투성 다뇨, 포도당 문턱, 암모니아 생성 (Renal handling)

67. DeFronzo RA, Hompesch M, Kasichayanula S, et al. **Characterization of renal glucose reabsorption in response to dapagliflozin in healthy subjects and subjects with type 2 diabetes.** Diabetes Care 2013;36:3169-76. — Tm glucose, threshold and splay, and the ~50-60% reduction with SGLT2 inhibition. https://pubmed.ncbi.nlm.nih.gov/23735727/
68. Rave K, Nosek L, Posner J, Heise T, Roggen K, van Hoogdalem EJ. **Renal glucose excretion as a function of blood glucose concentration in subjects with type 2 diabetes.** Diabetologia 2006;49:1274-82. https://pubmed.ncbi.nlm.nih.gov/16596360/
69. Vallon V, Thomson SC. **Targeting renal glucose reabsorption to treat hyperglycaemia: the pleiotropic effects of SGLT2 inhibition.** Diabetologia 2017;60:215-25. https://pubmed.ncbi.nlm.nih.gov/27878313/
70. Thomson SC, Blantz RC. **Glomerulotubular balance, tubuloglomerular feedback, and salt homeostasis.** J Am Soc Nephrol 2008;19:2272-5. — the reason Tm falls less than filtration does in a prerenal state (FTM_GFR in the model). https://pubmed.ncbi.nlm.nih.gov/19005009/
71. Gennari FJ, Kassirer JP. **Osmotic diuresis.** N Engl J Med 1974;291:714-20. — the fractional-excretion multipliers the model uses for the sodium and potassium losses of DKA. https://pubmed.ncbi.nlm.nih.gov/4606663/
72. Halperin ML, Cheema-Dhadli S, Lin SH, Kamel KS. **Properties permitting the renal cortex to be the oxygen sensor for the release of erythropoietin: clinical implications.** Clin J Am Soc Nephrol 2006;1:1049-53 (for the ammoniagenesis-PEPCK coupling discussed therein). https://pubmed.ncbi.nlm.nih.gov/17699325/
73. Halperin ML, Jungas RL. **Metabolic production and renal disposal of hydrogen ions.** Kidney Int 1983;24:709-13. — the ledger the model keeps: ketoanion with NH4⁺ preserves base, ketoanion with Na⁺ loses it. https://pubmed.ncbi.nlm.nih.gov/6323704/
74. Nissim I, Yudkoff M, Segal S. **Metabolic fate of glutamine in rat renal tubules: a study with 15N.** Metabolism 1987;36:1057-63. https://pubmed.ncbi.nlm.nih.gov/3670093/
75. Weiner ID, Verlander JW. **Renal ammonia metabolism and transport.** Compr Physiol 2013;3:201-20. — the days-long time course of ammoniagenic adaptation (the model's TNH4 = 12 h and NH4MAX ceiling). https://pubmed.ncbi.nlm.nih.gov/23720285/
76. Orth SR, Ritz E. **Adaptation of renal ammoniagenesis to metabolic acidosis.** Miner Electrolyte Metab 1996;22:118-25. https://pubmed.ncbi.nlm.nih.gov/8676808/

---

## 9. 수분·나트륨·삼투압 (Water, sodium and osmolality)

77. Katz MA. **Hyperglycemia-induced hyponatremia — calculation of expected serum sodium depression.** N Engl J Med 1973;289:843-4. https://pubmed.ncbi.nlm.nih.gov/4763428/
78. Hillier TA, Abbott RD, Barrett EJ. **Hyponatremia: evaluating the correction factor for hyperglycemia.** Am J Med 1999;106:399-403. — the 0.024 mmol/L per mg/dL correction the model uses for corrected sodium. https://pubmed.ncbi.nlm.nih.gov/10225241/
79. Adrogué HJ, Madias NE. **Hypernatremia.** N Engl J Med 2000;342:1493-9. https://pubmed.ncbi.nlm.nih.gov/10816188/
80. Kreisberg RA. **Diabetic ketoacidosis: new concepts and trends in pathogenesis and treatment.** Ann Intern Med 1978;88:681-95. https://pubmed.ncbi.nlm.nih.gov/417653/
81. Waldhäusl W, Kleinberger G, Korn A, Dudczak R, Bratusch-Marrain P, Nowotny P. **Severe hyperglycemia: effects of rehydration on endocrine derangements and blood glucose concentration.** Diabetes 1979;28:577-84. — fluid alone lowers glucose substantially while leaving ketone bodies almost unchanged: the two-lever result the model reproduces. https://pubmed.ncbi.nlm.nih.gov/109354/
82. Owen OE, Licht JH, Sapir DG. **Renal function and effects of partial rehydration during diabetic ketoacidosis.** Diabetes 1981;30:510-8. https://pubmed.ncbi.nlm.nih.gov/6785122/
83. Luzi L, Barrett EJ, Groop LC, Ferrannini E, DeFronzo RA. **Metabolic effects of low-dose insulin therapy on glucose metabolism in diabetic ketoacidosis.** Diabetes 1988;37:1470-7. — measured glucose production and disposal during treatment, against which the model's flux decomposition is calibrated. https://pubmed.ncbi.nlm.nih.gov/3141232/

---

## 10. 칼륨·인·마그네슘 (Potassium, phosphate, magnesium)

84. Adrogué HJ, Lederer ED, Suki WN, Eknoyan G. **Determinants of plasma potassium levels in diabetic ketoacidosis.** Medicine 1986;65:163-72. — the finding the model's potassium set-point equation encodes: insulin lack and hypertonicity, not the acid itself, drive the hyperkalaemia. https://pubmed.ncbi.nlm.nih.gov/3084904/
85. Adrogué HJ, Madias NE. **Changes in plasma potassium concentration during acute acid-base disturbances.** Am J Med 1981;71:456-67. — organic acidoses shift potassium far less than mineral acidoses do. https://pubmed.ncbi.nlm.nih.gov/7025622/
86. Aronson PS, Giebisch G. **Effects of pH on potassium: new explanations for old observations.** J Am Soc Nephrol 2011;22:1981-9. https://pubmed.ncbi.nlm.nih.gov/22034506/
87. Murthy K, Harrington JT, Siegel RD. **Profound hypokalemia in diabetic ketoacidosis: a therapeutic challenge.** Endocr Pract 2005;11:331-4. https://pubmed.ncbi.nlm.nih.gov/16191494/
88. Fisher JN, Kitabchi AE. **A randomized study of phosphate therapy in the treatment of diabetic ketoacidosis.** J Clin Endocrinol Metab 1983;57:177-80. https://pubmed.ncbi.nlm.nih.gov/6408100/
89. Wilson HK, Keuer SP, Lea AS, Boyd AE, Eknoyan G. **Phosphate therapy in diabetic ketoacidosis.** Arch Intern Med 1982;142:517-20. https://pubmed.ncbi.nlm.nih.gov/6802093/
90. Ditzel J, Lundgaard Hansen H. **Disturbance of inorganic phosphate metabolism in diabetes mellitus: clinical implications.** Diabetes Metab Syndr Obes 2010;3:319-24. https://pubmed.ncbi.nlm.nih.gov/21437101/

---

## 11. 중탄산염 치료 (Bicarbonate therapy)

91. Morris LR, Murphy MB, Kitabchi AE. **Bicarbonate therapy in severe diabetic ketoacidosis.** Ann Intern Med 1986;105:836-40. — the randomised trial finding no benefit, which the model reproduces mechanistically. https://pubmed.ncbi.nlm.nih.gov/3096181/
92. Chua HR, Schneider A, Bellomo R. **Bicarbonate in diabetic ketoacidosis — a systematic review.** Ann Intensive Care 2011;1:23. https://pubmed.ncbi.nlm.nih.gov/21906367/
93. Duhon B, Attridge RL, Franco-Martinez AC, Maxwell PR, Hughes DW. **Intravenous sodium bicarbonate therapy in severely acidotic diabetic ketoacidosis.** Ann Pharmacother 2013;47:970-5. https://pubmed.ncbi.nlm.nih.gov/23821610/
94. Ohman JL Jr, Marliss EB, Aoki TT, Munichoodappa CS, Khanna VV, Kozak GP. **The cerebrospinal fluid in diabetic ketoacidosis.** N Engl J Med 1971;284:283-90. — the paradoxical CSF acidification after bicarbonate. https://pubmed.ncbi.nlm.nih.gov/4992707/

---

## 12. 뇌부종 (Cerebral oedema)

95. Glaser N, Barnett P, McCaslin I, et al. **Risk factors for cerebral edema in children with diabetic ketoacidosis.** N Engl J Med 2001;344:264-9. — the risk factors the model reproduces (low PCO₂, high urea at presentation) and the ones it does not need (rate of fluid administration). https://pubmed.ncbi.nlm.nih.gov/11172152/
96. Kuppermann N, Ghetti S, Schunk JE, et al (PECARN DKA FLUID Study Group). **Clinical trial of fluid infusion rates for pediatric diabetic ketoacidosis.** N Engl J Med 2018;378:2275-87. — the negative randomised trial the model reproduces: a 6-fold range of fluid rate barely moves brain volume. https://pubmed.ncbi.nlm.nih.gov/29897851/
97. Glaser NS, Wootton-Gorges SL, Marcin JP, et al. **Mechanism of cerebral edema in children with diabetic ketoacidosis.** J Pediatr 2004;145:164-71. — the shift away from a purely osmotic explanation towards hypoperfusion and reperfusion injury. https://pubmed.ncbi.nlm.nih.gov/15289761/
98. Glaser NS, Marcin JP, Wootton-Gorges SL, et al. **Correlation of clinical and biochemical findings with DKA-related cerebral edema in children using magnetic resonance diffusion-weighted imaging.** J Pediatr 2008;153:541-6. https://pubmed.ncbi.nlm.nih.gov/18589447/
99. Gullans SR, Verbalis JG. **Control of brain volume during hyperosmolar and hypoosmolar conditions.** Annu Rev Med 1993;44:289-301. — the organic-osmolyte store whose slow washout the model represents as OSMB. https://pubmed.ncbi.nlm.nih.gov/8476253/
100. Lien YH, Shapiro JI, Chan L. **Effects of hypernatremia on organic brain osmoles.** J Clin Invest 1990;85:1427-35. — accumulation and washout time constants for taurine, myo-inositol and glycerophosphocholine. https://pubmed.ncbi.nlm.nih.gov/2332498/
101. Verbalis JG, Gullans SR. **Hyponatremia causes large sustained reductions in brain content of multiple organic osmolytes in rats.** Brain Res 1991;567:274-82. https://pubmed.ncbi.nlm.nih.gov/1817731/
102. Silver SM, Clark EC, Schroeder BM, Sterns RH. **Pathogenesis of cerebral edema after treatment of diabetic ketoacidosis.** Kidney Int 1997;51:1237-44. https://pubmed.ncbi.nlm.nih.gov/9083293/

---

## 13. 정상혈당 케톤산증과 SGLT2 억제제 (Euglycaemic ketoacidosis)

103. Peters AL, Buschur EO, Buse JB, Cohan P, Diner JC, Hirsch IB. **Euglycemic diabetic ketoacidosis: a potential complication of treatment with sodium-glucose cotransporter 2 inhibition.** Diabetes Care 2015;38:1687-93. https://pubmed.ncbi.nlm.nih.gov/26078479/
104. Ogawa W, Sakaguchi K. **Euglycemic diabetic ketoacidosis induced by SGLT2 inhibitors: possible mechanism and contributing factors.** J Diabetes Investig 2016;7:135-8. https://pubmed.ncbi.nlm.nih.gov/27042263/
105. Taylor SI, Blau JE, Rother KI. **SGLT2 inhibitors may predispose to ketoacidosis.** J Clin Endocrinol Metab 2015;100:2849-52. https://pubmed.ncbi.nlm.nih.gov/26086329/
106. Ferrannini E, Baldi S, Frascerra S, et al. **Shift to fatty substrate utilization in response to sodium-glucose cotransporter 2 inhibition in subjects without diabetes and patients with type 2 diabetes.** Diabetes 2016;65:1190-5. https://pubmed.ncbi.nlm.nih.gov/26861783/
107. Douros A, Lix LM, Fralick M, et al. **Sodium-glucose cotransporter-2 inhibitors and the risk for diabetic ketoacidosis: a multicenter cohort study.** Ann Intern Med 2020;173:417-25. https://pubmed.ncbi.nlm.nih.gov/32716707/
108. Munro JF, Campbell IW, McCuish AC, Duncan LJ. **Euglycaemic diabetic ketoacidosis.** BMJ 1973;2:578-80. — the syndrome described before SGLT2 inhibitors existed. https://pubmed.ncbi.nlm.nih.gov/4197425/

---

## 14. 알코올성 케톤산증과 임신 (Alcoholic ketoacidosis, pregnancy)

109. Wrenn KD, Slovis CM, Minion GE, Rutkowski R. **The syndrome of alcoholic ketoacidosis.** Am J Med 1991;91:119-28. — the high BHB:AcAc ratio at a normal glucose that the model's redox term produces. https://pubmed.ncbi.nlm.nih.gov/1867237/
110. McGuire LC, Cruickshank AM, Munro PT. **Alcoholic ketoacidosis.** Emerg Med J 2006;23:417-20. https://pubmed.ncbi.nlm.nih.gov/16714496/
111. Sibai BM, Viteri OA. **Diabetic ketoacidosis in pregnancy.** Obstet Gynecol 2014;123:167-78. — accelerated starvation and a lower buffer reserve, the two changes used in the pregnancy phenotype. https://pubmed.ncbi.nlm.nih.gov/24463676/
112. Felig P, Lynch V. **Starvation in human pregnancy: hypoalaninemia and hypoglycemia accompanied by hyperketonemia.** Science 1970;170:990-2. https://pubmed.ncbi.nlm.nih.gov/5529067/

---

## 15. 피하 인슐린 프로토콜과 전환 (Subcutaneous protocols and the transition)

113. Umpierrez GE, Latif K, Stoever J, et al. **Efficacy of subcutaneous insulin lispro versus continuous intravenous regular insulin for the treatment of patients with diabetic ketoacidosis.** Am J Med 2004;117:291-6. https://pubmed.ncbi.nlm.nih.gov/15336577/
114. Umpierrez GE, Cuervo R, Karabell A, Latif K, Freire AX, Kitabchi AE. **Treatment of diabetic ketoacidosis with subcutaneous insulin aspart.** Diabetes Care 2004;27:1873-8. https://pubmed.ncbi.nlm.nih.gov/15277410/
115. Andrade-Castellanos CA, Colunga-Lozano LE, Delgado-Figueroa N, Gonzalez-Padilla DA. **Subcutaneous rapid-acting insulin analogues for diabetic ketoacidosis.** Cochrane Database Syst Rev 2016;1:CD011281. https://pubmed.ncbi.nlm.nih.gov/26798030/
116. Homko C, Deluzio A, Jimenez C, Kolaczynski JW, Boden G. **Comparison of insulin aspart and lispro: pharmacokinetic and metabolic effects.** Diabetes Care 2003;26:2027-31. — the absorption rate constants used for the s.c. depot. https://pubmed.ncbi.nlm.nih.gov/12832307/
117. Mudaliar SR, Lindberg FA, Joyce M, et al. **Insulin aspart (B28 asp-insulin): a fast-acting analog of human insulin.** Diabetes Care 1999;22:1501-6. https://pubmed.ncbi.nlm.nih.gov/10480516/
118. Hsia E, Seggelke S, Gibbs J, et al. **Subcutaneous administration of glargine to diabetic patients receiving insulin infusion prevents rebound hyperglycemia.** J Clin Endocrinol Metab 2012;97:3132-7. — the overlap the model shows is not a formality. https://pubmed.ncbi.nlm.nih.gov/22685233/

---

## 16. 대사반응 조절 호르몬과 인슐린 저항성 (Counter-regulatory physiology)

119. Gerich JE, Lorenzi M, Bier DM, et al. **Prevention of human diabetic ketoacidosis by somatostatin: evidence for an essential role of glucagon.** N Engl J Med 1975;292:985-9. https://pubmed.ncbi.nlm.nih.gov/804137/
120. Miles JM, Rizza RA, Haymond MW, Gerich JE. **Effects of acute insulin deficiency on glucose and ketone body turnover in man.** Diabetes 1980;29:926-30. https://pubmed.ncbi.nlm.nih.gov/6773829/
121. Barrett EJ, DeFronzo RA, Bevilacqua S, Ferrannini E. **Insulin resistance in diabetic ketoacidosis.** Diabetes 1982;31:923-8. — the resistance that resolves as the acidosis and the NEFA concentration fall. https://pubmed.ncbi.nlm.nih.gov/6813166/
122. Hirsch IB, Emmett M. **Diabetic ketoacidosis and hyperosmolar hyperglycemic state in adults: epidemiology and pathogenesis.** UpToDate (continuously updated) — used for cross-checking conventional teaching against model output.
123. Kitabchi AE, Umpierrez GE, Fisher JN, Murphy MB, Stentz FB. **Thirty years of personal experience in hyperglycemic crises: DKA and HHS.** J Clin Endocrinol Metab 2008;93:1541-52. https://pubmed.ncbi.nlm.nih.gov/18270259/
124. Stentz FB, Umpierrez GE, Cuervo R, Kitabchi AE. **Proinflammatory cytokines, markers of cardiovascular risks, oxidative stress, and lipid peroxidation in patients with hyperglycemic crises.** Diabetes 2004;53:2079-86. https://pubmed.ncbi.nlm.nih.gov/15277389/

---

## 17. 측정 인공물 (Measurement artefacts the model reproduces)

125. Cotten SW, Willis MS. **Interference of acetoacetate and other compounds in creatinine assays.** (Jaffe reaction interference; see also Owen 1982 above.) Clin Chim Acta reviews of picrate-method interference. https://pubmed.ncbi.nlm.nih.gov/25445815/
126. Molitch ME, Rodman E, Hirsch CA, Dubinsky E. **Spurious serum creatinine elevations in ketoacidosis.** Ann Intern Med 1980;93:280-1. — the artefact captured by CREA_JAFFE in the model. https://pubmed.ncbi.nlm.nih.gov/7406382/
127. Kemperman FA, Weber JA, Gorgels J, van Zanten AP, Krediet RT, Arisz L. **The influence of ketoacids on plasma creatinine assays in diabetic ketoacidosis.** J Intern Med 2000;248:511-7. https://pubmed.ncbi.nlm.nih.gov/11155144/

---

## 18. 모델링 방법론 (Modelling methodology)

128. Baron KT, Gastonguay MR. **mrgsolve: Simulate from ODE-based population PK/PD and systems pharmacology models.** R package. https://mrgsolve.org
129. Bergman RN, Ider YZ, Bowden CR, Cobelli C. **Quantitative estimation of insulin sensitivity.** Am J Physiol 1979;236:E667-77. — the minimal-model tradition this model deliberately departs from by adding the renal and acid-base arms. https://pubmed.ncbi.nlm.nih.gov/443421/
130. Dalla Man C, Rizza RA, Cobelli C. **Meal simulation model of the glucose-insulin system.** IEEE Trans Biomed Eng 2007;54:1740-9. https://pubmed.ncbi.nlm.nih.gov/17926672/
131. Hovorka R, Canonico V, Chassin LJ, et al. **Nonlinear model predictive control of glucose concentration in subjects with type 1 diabetes.** Physiol Meas 2004;25:905-20. — the effect-compartment structure adopted here for peripheral insulin action. https://pubmed.ncbi.nlm.nih.gov/15382830/
132. Guyton AC, Coleman TG, Granger HJ. **Circulation: overall regulation.** Annu Rev Physiol 1972;34:13-46. — the whole-body fluid and electrolyte accounting tradition the volume block follows. https://pubmed.ncbi.nlm.nih.gov/4334846/

---

## 근거의 한계 (Where the evidence is thin)

- **뇌부종의 정량적 기전.** The osmotic and ischaemic contributions to cerebral oedema are not separately measurable in patients; the model's split between them (KOSMB/COMPL versus KINJ/KVASO) reproduces the *reported risk factors* and the negative fluid-rate trial, but the partition itself is an inference, not a measurement.
- **케톤 처리 용량의 개인차.** Ketone oxidation Vmax in acutely ill humans is known only from small tracer studies; the saturable/non-saturable split (VMAX_KOX vs KLIN_KOX) is constrained by steady-state concentrations rather than measured directly.
- **문맥 인슐린 비율.** PORTF is set from first-pass extraction studies in the fasting state; whether the same ratio holds during a hyperglycaemic crisis with hepatic congestion is unknown, and this parameter carries the whole DKA-versus-HHS distinction.
- **중탄산염의 조건수.** As documented in the model output, the computed bicarbonate is a small residual of large strong-ion terms and is correspondingly sensitive to renal chloride handling. Its predictions should be treated as less reliable than those for the anion gap and β-hydroxybutyrate.
