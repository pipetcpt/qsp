# 살리실산(아스피린) 중독 QSP 모델 — 참고문헌
# Salicylate (Aspirin) Poisoning QSP Model — References

모든 PMID는 PubMed E-utilities로 조회하여 확인했습니다.
Every PMID below was resolved against PubMed rather than written from memory.
Entries marked **[CAL]** supplied a number that was *fitted*; entries marked
**[VAL]** supplied a number the model was *tested against* but not fitted to.

---

## 1. 총설 · 역학 (Reviews and epidemiology)

1. Palmer BF, Clegg DJ. **Salicylate Toxicity.** *N Engl J Med* 2020;382:2544-2555. — the modern reference review; the source of the management structure used in the scenarios. [PMID 32579814](https://pubmed.ncbi.nlm.nih.gov/32579814/)
2. Temple AR. **Pathophysiology of aspirin overdosage toxicity, with implications for management.** *Pediatrics* 1978;62(5 Pt 2):873-876. — still the clearest short statement of why the acid-base disturbance and the tissue distribution are the same problem. [PMID 364398](https://pubmed.ncbi.nlm.nih.gov/364398/)
3. Temple AR. **Acute and chronic effects of aspirin toxicity and their treatment.** *Arch Intern Med* 1981;141:364-369. [PMID 7469627](https://pubmed.ncbi.nlm.nih.gov/7469627/)
4. Hill JB. **Salicylate intoxication.** *N Engl J Med* 1973;288:1110-1113. [PMID 4572648](https://pubmed.ncbi.nlm.nih.gov/4572648/)
5. Krause DS, Wolf BA, Shaw LM. **Acute aspirin overdose: mechanisms of toxicity.** *Ther Drug Monit* 1992;14:441-451. [PMID 1485363](https://pubmed.ncbi.nlm.nih.gov/1485363/)
6. Meredith TJ, Vale JA. **Non-narcotic analgesics. Problems of overdosage.** *Drugs* 1986;32(Suppl 4):177-205. [PMID 3552583](https://pubmed.ncbi.nlm.nih.gov/3552583/)
7. O'Malley GF. **Emergency department management of the salicylate-poisoned patient.** *Emerg Med Clin North Am* 2007;25:333-346. [PMID 17482023](https://pubmed.ncbi.nlm.nih.gov/17482023/)
8. Runde TJ, Nappe TM. **Salicylates Toxicity.** *StatPearls*, updated 2024. [PMID 29763054](https://pubmed.ncbi.nlm.nih.gov/29763054/)
9. Reingardiene D, Lazauskas R. **Acute salicylate poisoning.** *Medicina (Kaunas)* 2006;42:79-88. [PMID 16467617](https://pubmed.ncbi.nlm.nih.gov/16467617/)
10. Gummin DD, et al. **2023 Annual Report of the National Poison Data System (NPDS): 41st Annual Report.** *Clin Toxicol (Phila)* 2024;62:793-1027. — exposure and fatality denominators. [PMID 39688840](https://pubmed.ncbi.nlm.nih.gov/39688840/)
11. Shively RM, Hoffman RS, Manini AF. **Acute salicylate poisoning: risk factors for severe outcome.** *Clin Toxicol (Phila)* 2017;55:175-180. **[VAL]** — the model's severity index is checked against this risk structure. [PMID 28064509](https://pubmed.ncbi.nlm.nih.gov/28064509/)
12. McGuigan MA. **A two-year review of salicylate deaths in Ontario.** *Arch Intern Med* 1987;147:510-512. — the chronic/elderly death phenotype that scenario S12 reproduces. [PMID 3827428](https://pubmed.ncbi.nlm.nih.gov/3827428/)

---

## 2. 약동학 — 포화 대사 (Pharmacokinetics: capacity-limited metabolism)

13. Levy G. **Pharmacokinetics of salicylate elimination in man.** *J Pharm Sci* 1965;54:959-967. **[CAL]** — the Michaelis-Menten parameterisation of salicylurate and phenolic glucuronide formation used in `$PARAM` (VMAXSU/KMSU, VMAXPG/KMPG). [PMID 5862532](https://pubmed.ncbi.nlm.nih.gov/5862532/)
14. Levy G. **Clinical pharmacokinetics of salicylates: a re-assessment.** *Br J Clin Pharmacol* 1980;10(Suppl 2):285S-290S. **[CAL]** [PMID 7437270](https://pubmed.ncbi.nlm.nih.gov/7437270/)
15. Levy G, Tsuchiya T. **Salicylate accumulation kinetics in man.** *N Engl J Med* 1972;287:430-432. **[CAL]** — the non-linear rise of steady-state concentration with dose; anchors the 3.9 g/day → ~200 mg/L check. [PMID 5044917](https://pubmed.ncbi.nlm.nih.gov/5044917/)
16. Bedford C, Cummings AJ, Martin BK. **A kinetic study of the elimination of salicylate in man.** *Br J Pharmacol Chemother* 1965;24:418-431. [PMID 14320855](https://pubmed.ncbi.nlm.nih.gov/14320855/)
17. Needs CJ, Brooks PM. **Clinical pharmacokinetics of the salicylates.** *Clin Pharmacokinet* 1985;10:164-177. **[CAL]** — half-life, bioavailability and metabolite fractions. [PMID 3888490](https://pubmed.ncbi.nlm.nih.gov/3888490/)
18. Furst DE, Tozer TN, Melmon KL. **Salicylate clearance, the resultant of protein binding and metabolism.** *Clin Pharmacol Ther* 1979;26:380-389. — the paper that states the model's central PK point: clearance is a *product* of binding and enzyme terms, so neither can be interpreted alone. [PMID 466931](https://pubmed.ncbi.nlm.nih.gov/466931/)
19. Roberts DM, Buckley NA. **Pharmacokinetic considerations in clinical toxicology: clinical applications.** *Clin Pharmacokinet* 2007;46:897-939. [PMID 17922558](https://pubmed.ncbi.nlm.nih.gov/17922558/)

---

## 3. 단백결합 포화 — 곱셈 인자 1 (Protein binding: MULTIPLIER 1)

20. Ekstrand R, Alván G, Orme ML, et al. **Concentration dependent plasma protein binding of salicylate in rheumatoid patients.** *Clin Pharmacokinet* 1979;4:137-143. **[CAL]** — the concentration-dependence of the free fraction that the two-class binding model (B1/KD1, B2/KD2) was fitted to. [PMID 378501](https://pubmed.ncbi.nlm.nih.gov/378501/)
21. Furst DE, Tozer TN, Melmon KL. *(as above, ref 18)* — free fraction rising from ~8% to >40% across the toxic range.
22. Bannwarth B, Netter P, Pourel J, et al. **Clinical pharmacokinetics of nonsteroidal anti-inflammatory drugs in the cerebrospinal fluid.** *Biomed Pharmacother* 1989;43:121-126. — CSF concentrations track the *unbound* plasma concentration, not the total. [PMID 2660917](https://pubmed.ncbi.nlm.nih.gov/2660917/)

---

## 4. pH 분배와 조직 분포 — 곱셈 인자 2 (pH partition: MULTIPLIER 2)

23. Hill JB. **Experimental salicylate poisoning: observations on the effects of altering blood pH on tissue and plasma salicylate concentrations.** *Pediatrics* 1971;47:658-665. **[CAL]** — the key experiment behind the whole model: lowering blood pH *raises* tissue and brain salicylate while *lowering* the plasma concentration. Scenario S07 vs S08 is this experiment transposed to the ventilator. [PMID 5089754](https://pubmed.ncbi.nlm.nih.gov/5089754/)
24. Macpherson CR, Milne MD, Evans BM. **The excretion of salicylate.** *Br J Pharmacol Chemother* 1955;10:484-489. **[CAL]** — non-ionic back-diffusion of a weak acid across the tubular epithelium; the origin of the `AREAB` term. [PMID 13276608](https://pubmed.ncbi.nlm.nih.gov/13276608/)
25. Morgan AG, Polak A. **The excretion of salicylate in salicylate poisoning.** *Clin Sci* 1971;41:475-484. **[VAL]** — the measured dependence of renal salicylate clearance on urine pH. [PMID 5123233](https://pubmed.ncbi.nlm.nih.gov/5123233/)

### 4a. 세포내 pH 완충 — 왜 CO2와 HCO3가 대칭이 아닌가

26. Portman MA, Lassen NA, Cooper TG, et al. **Intra- and extracellular pH of the brain in vivo studied by 31P-NMR during hyper- and hypocapnia.** *J Appl Physiol* 1991;71:2168-2172. **[CAL]** — the observation that fixes `BETAB`: brain intracellular pH is defended against acute PaCO2 change while plasma pH is not. Without this the entire ventilator effect disappears. [PMID 1778908](https://pubmed.ncbi.nlm.nih.gov/1778908/)
27. Cadoux-Hudson TA, Blackledge MJ, Rajagopalan B, et al. **Response of the human brain to a hypercapnic acid load in vivo.** *Clin Sci (Lond)* 1990;79:1-3. **[CAL]** [PMID 2167784](https://pubmed.ncbi.nlm.nih.gov/2167784/)
28. Mellergård P, Ouyang YB, Siesjö BK. **Intracellular pH regulation in cultured rat astrocytes in CO2/HCO3-containing media.** *Exp Brain Res* 1993;95:371-380. — the transport machinery that makes the brain buffer base follow plasma bicarbonate slowly (`TAUBB`). [PMID 8224063](https://pubmed.ncbi.nlm.nih.gov/8224063/)

---

## 5. 미토콘드리아 탈공역 (Mitochondrial uncoupling)

29. Haas R, Parker WD Jr, Stumpf D, Eguren LA. **Salicylate-induced loose coupling: protonmotive force measurements.** *Biochem Pharmacol* 1985;34:900-902. **[CAL]** — the protonophore mechanism and its concentration dependence (`EC50U`). [PMID 3977963](https://pubmed.ncbi.nlm.nih.gov/3977963/)
30. Kaplan EH, Kennedy J, Davis J. **Effects of salicylate and other benzoates on oxidative enzymes of the tricarboxylic acid cycle in rat tissue homogenates.** *Arch Biochem Biophys* 1954;51:47-61. [PMID 13181459](https://pubmed.ncbi.nlm.nih.gov/13181459/)
31. Trost LC, Lemasters JJ. **Role of the mitochondrial permeability transition in salicylate toxicity to cultured rat hepatocytes: implications for the pathogenesis of Reye's syndrome.** *Toxicol Appl Pharmacol* 1997;147:431-441. [PMID 9439738](https://pubmed.ncbi.nlm.nih.gov/9439738/)
32. Trost LC, Lemasters JJ. **The mitochondrial permeability transition: a new pathophysiological mechanism for Reye's syndrome and toxic liver injury.** *J Pharmacol Exp Ther* 1996;278:1000-1005. [PMID 8819478](https://pubmed.ncbi.nlm.nih.gov/8819478/)
33. Schrör K. **Aspirin and Reye syndrome: a review of the evidence.** *Paediatr Drugs* 2007;9:195-204. [PMID 17523700](https://pubmed.ncbi.nlm.nih.gov/17523700/)
34. Glasgow JF. **Reye's syndrome: the case for a causal link with aspirin.** *Drug Saf* 2006;29:1111-1121. [PMID 17147458](https://pubmed.ncbi.nlm.nih.gov/17147458/)

---

## 6. 산-염기 · 호흡 (Acid-base and ventilation)

35. Gabow PA, Anderson RJ, Potts DE, Schrier RW. **Acid-base disturbances in the salicylate-intoxicated adult.** *Arch Intern Med* 1978;138:1481-1484. **[VAL]** — the mixed respiratory-alkalosis-plus-metabolic-acidosis pattern, and the observation that most adults present *alkalaemic*. Scenario S02 is checked against this. [PMID 708168](https://pubmed.ncbi.nlm.nih.gov/708168/)
36. Tenney SM, Miller RM. **The respiratory and circulatory actions of salicylate.** *Am J Med* 1955;19:498-508. **[CAL]** — direct medullary stimulation independent of acid-base status; the source of `EMAXR`/`EC50R`. [PMID 13258583](https://pubmed.ncbi.nlm.nih.gov/13258583/)
37. Anderson RJ, Potts DE, Gabow PA, Rumack BH, Schrier RW. **Unrecognized adult salicylate intoxication.** *Ann Intern Med* 1976;85:745-748. — the classic series in which the diagnosis was missed for a mean of several days. [PMID 999110](https://pubmed.ncbi.nlm.nih.gov/999110/)

### 6a. 인공호흡의 위험 — 이 모델의 중심 실험

38. Stolbach AI, Hoffman RS, Nelson LS. **Mechanical ventilation was associated with acidemia in a case series of salicylate-poisoned patients.** *Acad Emerg Med* 2008;15:866-869. **[VAL]** — the observational counterpart of scenario S07: intubated salicylate-poisoned patients became acidaemic. The model reproduces the direction and the approximate magnitude without being fitted to it. [PMID 18821862](https://pubmed.ncbi.nlm.nih.gov/18821862/)
39. Greenberg MI, Hendrickson RG, Hofman M. **Deleterious effects of endotracheal intubation in salicylate poisoning.** *Ann Emerg Med* 2003;41:583-584. **[VAL]** [PMID 12705252](https://pubmed.ncbi.nlm.nih.gov/12705252/)
40. McCabe DJ, Lu JJ. **The association of hemodialysis and survival in intubated salicylate-poisoned patients.** *Am J Emerg Med* 2017;35:899-903. **[VAL]** — the interaction the model predicts: once the patient is intubated, only extracorporeal removal is fast enough. [PMID 28438446](https://pubmed.ncbi.nlm.nih.gov/28438446/)

---

## 7. 신장 처리 · 요 알칼리화 (Renal handling and urinary alkalinisation)

41. Proudfoot AT, Krenzelok EP, Vale JA. **Position Paper on urine alkalinization.** *J Toxicol Clin Toxicol* 2004;42:1-26. **[CAL]** — the regimen implemented as `ALKALI_1`/`ALKALI_2`, the urine pH target of 7.5-8.5, the arterial pH ceiling, and the instruction that potassium must be replaced. [PMID 15083932](https://pubmed.ncbi.nlm.nih.gov/15083932/)
42. Prescott LF, Balali-Mood M, Critchley JA, Johnstone AF, Proudfoot AT. **Diuresis or urinary alkalinisation for salicylate poisoning?** *Br Med J (Clin Res Ed)* 1982;285:1383-1386. **[VAL]** — the comparison the model must get right: alkalinisation, not volume, is what increases clearance. Reproduced as S04 vs S05. *(Search: "Diuresis or urinary alkalinisation for salicylate poisoning" — BMJ 1982;285:1383)*
43. Berg KJ. **Acute acetylsalicylic acid poisoning: treatment with forced alkaline diuresis and diuretics.** *Eur J Clin Pharmacol* 1977;12:111-116. [PMID 21797](https://pubmed.ncbi.nlm.nih.gov/21797/)
44. Ullmann E. **Factors modifying renal tubular bicarbonate reabsorption in the dog.** *J Physiol* 1968;198:1-14. **[CAL]** — the dependence of the bicarbonate threshold on PCO2, which in this model is why the patient's own hyperventilation makes the urine easier to alkalinise (and why intubation makes it harder). [PMID 5636989](https://pubmed.ncbi.nlm.nih.gov/5636989/)
45. Cogan MG, Liu FY. **Metabolic alkalosis in the rat. Evidence that reduced glomerular filtration rather than enhanced tubular bicarbonate reabsorption is responsible for maintaining the alkalotic state.** *J Clin Invest* 1983;71:1141-1160. **[CAL]** — the volume/potassium dependence of bicarbonate handling behind `THK` and `THVOL`. [PMID 6853706](https://pubmed.ncbi.nlm.nih.gov/6853706/)
46. Higgins RM, Connolly JO, Hendry BM. **Alkalinization and hemodialysis in severe salicylate poisoning: comparison of elimination techniques in the same patient.** *Clin Nephrol* 1998;50:178-183. **[VAL]** — a within-patient comparison of the two clearances; the model's 6 L/h vs 2.7 L/h ratio is checked against it. [PMID 9776422](https://pubmed.ncbi.nlm.nih.gov/9776422/)

---

## 8. 체외 제거 (Extracorporeal removal)

47. Juurlink DN, Gosselin S, Kielstein JT, et al. (EXTRIP Workgroup). **Extracorporeal Treatment for Salicylate Poisoning: Systematic Review and Recommendations From the EXTRIP Workgroup.** *Ann Emerg Med* 2015;66:165-181. **[CAL]** — the indication thresholds encoded in the Shiny app's dialysis tab. [PMID 25986310](https://pubmed.ncbi.nlm.nih.gov/25986310/)
48. Jacobsen D, Wiik-Larsen E, Bredesen JE. **Haemodialysis or haemoperfusion in severe salicylate poisoning?** *Hum Toxicol* 1988;7:161-163. **[CAL]** — the extracorporeal clearance value (`CLHD`). [PMID 3378803](https://pubmed.ncbi.nlm.nih.gov/3378803/)
49. Wrathall G, Sinclair R, Moore A, Pogson D. **Three case reports of the use of haemodiafiltration in the treatment of salicylate overdose.** *Hum Exp Toxicol* 2001;20:491-495. [PMID 11776412](https://pubmed.ncbi.nlm.nih.gov/11776412/)
50. Fertel BS, Nelson LS, Goldfarb DS. **The underutilization of hemodialysis in patients with salicylate poisoning.** *Kidney Int* 2009;75:1349-1353. **[VAL]** — the clinical cost of the delay that scenario S10 quantifies. [PMID 18716600](https://pubmed.ncbi.nlm.nih.gov/18716600/)

---

## 9. 위장관 제염 · 지연 흡수 (Decontamination and delayed absorption)

51. Wortzman DJ, Grunfeld A. **Delayed absorption following enteric-coated aspirin overdose.** *Ann Emerg Med* 1987;16:434-436. **[VAL]** — the late peak reproduced in scenario S11. [PMID 3826813](https://pubmed.ncbi.nlm.nih.gov/3826813/)
52. Pierce RP, Gazewood J, Blake RL Jr. **Salicylate poisoning from enteric-coated aspirin. Delayed absorption may complicate management.** *Postgrad Med* 1991;89:61-64. [PMID 2008403](https://pubmed.ncbi.nlm.nih.gov/2008403/)
53. Bogacz K, Caldron P. **Enteric-coated aspirin bezoar: elevation of serum salicylate level by barium study.** *Am J Med* 1987;83:783-786. — the concretion compartment (`ACONC`). [PMID 2823600](https://pubmed.ncbi.nlm.nih.gov/2823600/)
54. Rivera W, Kleinschmidt KC, Velez LI, et al. **Delayed salicylate toxicity at 35 hours without early manifestations following a single salicylate ingestion.** *Ann Pharmacother* 2004;38:1186-1188. **[VAL]** [PMID 15173556](https://pubmed.ncbi.nlm.nih.gov/15173556/)
55. Moss MJ, Ruha AM, Gerkin RD, et al. **Salicylate toxicity after undetectable serum salicylate concentration: a retrospective cohort study.** *Clin Toxicol (Phila)* 2019;57:832-837. **[VAL]** — the single-level trap the model reproduces as the "3 h level understates the peak" check. [PMID 30306804](https://pubmed.ncbi.nlm.nih.gov/30306804/)
56. Hillman RJ, Prescott LF. **Treatment of salicylate poisoning with repeated oral charcoal.** *Br Med J (Clin Res Ed)* 1985;291:1472. [PMID 3933714](https://pubmed.ncbi.nlm.nih.gov/3933714/)
57. Barone JA, Raia JJ, Huang YC. **Evaluation of the effects of multiple-dose activated charcoal on the absorption of orally administered salicylate in a simulated toxic ingestion model.** *Ann Emerg Med* 1988;17:34-37. **[CAL]** [PMID 3337412](https://pubmed.ncbi.nlm.nih.gov/3337412/)
58. Kirshenbaum LA, Mathews SC, Sitar DS, Tenenbein M. **Does multiple-dose charcoal therapy enhance salicylate excretion?** *Arch Intern Med* 1990;150:1281-1283. — the honest negative result: MDAC works on absorption, not on elimination. The model reproduces this asymmetry. [PMID 2191636](https://pubmed.ncbi.nlm.nih.gov/2191636/)
59. Vertrees JE, McWilliams BC, Kelly HW. **Repeated oral administration of activated charcoal for treating aspirin overdose in young children.** *Pediatrics* 1990;85:594-598. [PMID 2314974](https://pubmed.ncbi.nlm.nih.gov/2314974/)

---

## 10. 중추신경 · 저혈당뇌증 (CNS injury and neuroglycopenia)

60. Thurston JH, Pollock PG, Warren SK, Jones EM. **Reduced brain glucose with normal plasma glucose in salicylate poisoning.** *J Clin Invest* 1970;49:2139-2145. **[CAL]** — the direct measurement behind the brain-glucose ODE and the clinical instruction to give dextrose whatever the serum glucose says. [PMID 4319971](https://pubmed.ncbi.nlm.nih.gov/4319971/)
61. Kuzak N, Brubacher JR, Kennedy JR. **Reversal of salicylate-induced euglycemic delirium with dextrose.** *Clin Toxicol (Phila)* 2007;45:526-529. **[VAL]** — scenario S14. [PMID 17503260](https://pubmed.ncbi.nlm.nih.gov/17503260/)
62. Cotton EK, Fahlberg VI. **Hypoglycemia with salicylate poisoning.** *Am J Dis Child* 1964;108:171-173. [PMID 14159936](https://pubmed.ncbi.nlm.nih.gov/14159936/)
63. Rauschka H, Aboul-Enein F, Bauer J, et al. **Acute cerebral white matter damage in lethal salicylate intoxication.** *Neurotoxicology* 2007;28:33-37. — post-mortem confirmation that the brain is the target organ. [PMID 16930716](https://pubmed.ncbi.nlm.nih.gov/16930716/)
64. Rosenfeld RG, Liebhaber MI. **Acute encephalopathy in siblings. Reye syndrome vs salicylate intoxication.** *Am J Dis Child* 1976;130:295-297. [PMID 1258838](https://pubmed.ncbi.nlm.nih.gov/1258838/)

---

## 11. 이독성 (Ototoxicity — the earliest clinical signal)

65. Cazals Y. **Auditory sensori-neural alterations induced by salicylate.** *Prog Neurobiol* 2000;62:583-631. **[CAL]** — the concentration-effect relationship behind the tinnitus readout. [PMID 10880852](https://pubmed.ncbi.nlm.nih.gov/10880852/)
66. Jung TT, Rhee CK, Lee CS, Park YS, Choi DC. **Ototoxicity of salicylate, nonsteroidal antiinflammatory drugs, and quinine.** *Otolaryngol Clin North Am* 1993;26:791-810. [PMID 8233489](https://pubmed.ncbi.nlm.nih.gov/8233489/)

---

## 12. 폐 · 기타 장기 (Pulmonary oedema and other organ effects)

67. Heffner JE, Sahn SA. **Salicylate-induced pulmonary edema. Clinical features and prognosis.** *Ann Intern Med* 1981;95:405-409. **[CAL]** — the age and chronicity dependence encoded in `AGEF`. [PMID 7283290](https://pubmed.ncbi.nlm.nih.gov/7283290/)
68. Fisher CJ Jr, Albertson TE, Foulke GE. **Salicylate-induced pulmonary edema: clinical characteristics in children.** *Am J Emerg Med* 1985;3:33-37. [PMID 3970751](https://pubmed.ncbi.nlm.nih.gov/3970751/)
69. Walters JS, Woodring JH, Stelling CB, Rosenbaum HD. **Salicylate-induced pulmonary edema.** *Radiology* 1983;146:289-293. [PMID 6849076](https://pubmed.ncbi.nlm.nih.gov/6849076/)
70. Goldsweig HG, Kapusta M, Schwartz J. **Bleeding, salicylates, and prolonged prothrombin time: three case reports and a review of the literature.** *J Rheumatol* 1976;3:37-42. [PMID 1271387](https://pubmed.ncbi.nlm.nih.gov/1271387/)

---

## 13. 만성 살리실산증 (Chronic salicylism — the same molecule, a different disease)

71. Bailey RB, Jones SR. **Chronic salicylate intoxication. A common cause of morbidity in the elderly.** *J Am Geriatr Soc* 1989;37:556-561. **[VAL]** — the phenotype scenario S12 is built to produce: severe toxicity at plasma levels a nomogram calls moderate. [PMID 2715563](https://pubmed.ncbi.nlm.nih.gov/2715563/)
72. Gaudreault P, Temple AR, Lovejoy FH Jr. **The relative severity of acute versus chronic salicylate poisoning in children: a clinical comparison.** *Pediatrics* 1982;70:566-569. **[VAL]** — the quantitative statement that the same level means more toxicity in chronic exposure. [PMID 7122154](https://pubmed.ncbi.nlm.nih.gov/7122154/)
73. Halani S, Chu W. **Salicylate toxicity from chronic bismuth subsalicylate use.** *BMJ Case Rep* 2020;13:e235500. [PMID 33257373](https://pubmed.ncbi.nlm.nih.gov/33257373/)
74. Wolowich WR, Hadley CM, Kelley MT, Walson PD, Casavant MJ. **Plasma salicylate from methyl salicylate cream compared to oil of wintergreen.** *J Toxicol Clin Toxicol* 2003;41:355-358. [PMID 12870876](https://pubmed.ncbi.nlm.nih.gov/12870876/)
75. Davis JE. **Are one or two dangerous? Methyl salicylate exposure in toddlers.** *J Emerg Med* 2007;32:63-69. [PMID 17239735](https://pubmed.ncbi.nlm.nih.gov/17239735/)

---

## 14. 노모그램과 측정의 실패 (The nomogram and the limits of measurement)

76. Done AK. **Salicylate intoxication. Significance of measurements of salicylate in blood in cases of acute ingestion.** *Pediatrics* 1960;26:800-807. — the original nomogram, evaluated as a model *output* in section 8 of the verification report rather than assumed. [PMID 13723722](https://pubmed.ncbi.nlm.nih.gov/13723722/)
77. Dugandzic RM, Tierney MG, Dickinson GE, Dolan MC, McKnight DR. **Evaluation of the validity of the Done nomogram in the management of acute salicylate intoxication.** *Ann Emerg Med* 1989;18:1186-1190. **[VAL]** — the empirical refutation the model reproduces mechanistically. [PMID 2817562](https://pubmed.ncbi.nlm.nih.gov/2817562/)
78. Done AK. **Treatment of salicylate poisoning.** *Mod Treat* 1971;8:528-551. [PMID 4941269](https://pubmed.ncbi.nlm.nih.gov/4941269/)

---

## 15. 의원성 위해 (Iatrogenic hazards)

79. Anderson CJ, Kaufman PL, Sturm RJ. **Toxicity of combined therapy with carbonic anhydrase inhibitors and aspirin.** *Am J Ophthalmol* 1978;86:516-519. — acetazolamide alkalinises the urine while acidifying the blood: in this model that is the single worst possible combination of the two multipliers. [PMID 707596](https://pubmed.ncbi.nlm.nih.gov/707596/)
80. Hazouard E, Ferrandiere M, Rateau H, et al. **[Salicylism and glaucoma: reciprocal augmentation of the toxicity of acetazolamide and acetylsalicylic acid].** *J Fr Ophtalmol* 1999;22:73-75. [PMID 10221197](https://pubmed.ncbi.nlm.nih.gov/10221197/)

---

## 16. 진료 지침 (Guidelines and management pathways)

81. Dargan PI, Wallace CI, Jones AL. **An evidence based flowchart to guide the management of acute salicylate (aspirin) overdose.** *Emerg Med J* 2002;19:206-209. [PMID 11971828](https://pubmed.ncbi.nlm.nih.gov/11971828/)
82. Chyka PA, Erdman AR, Christianson G, et al. **Salicylate poisoning: an evidence-based consensus guideline for out-of-hospital management.** *Clin Toxicol (Phila)* 2007;45:95-131. [PMID 17364628](https://pubmed.ncbi.nlm.nih.gov/17364628/)
83. Thisted B, Krantz T, Strøom J, Sørensen MB. **Acute salicylate self-poisoning in 177 consecutive patients treated in ICU.** *Acta Anaesthesiol Scand* 1987;31:312-316. **[VAL]** — the ICU denominator: mortality, ventilation rates, and the acid-base distribution at presentation. [PMID 3591255](https://pubmed.ncbi.nlm.nih.gov/3591255/)

---

## 17. 인용은 되었으나 PubMed 색인에서 확인하지 못한 항목

정직성을 위해 별도로 표시합니다. 아래 두 편은 본문에서 인용했으나 위의
E-utilities 조회로 PMID를 확정하지 못했으므로 서지사항만 기재합니다.
(Listed separately for honesty: cited in the text but their PMIDs could not be
confirmed by the E-utilities lookups above, so only the bibliographic details
are given.)

- Prescott LF, Balali-Mood M, Critchley JA, Johnstone AF, Proudfoot AT. **Diuresis or urinary alkalinisation for salicylate poisoning?** *Br Med J (Clin Res Ed)* 1982;285:1383-1386.
- Chapman BJ, Proudfoot AT. **Adverse effects of forced alkaline diuresis in the treatment of acute salicylate poisoning.** *Q J Med* 1989;72:699-707.

---

## 18. 모델링 방법론 (Modelling methodology)

- mrgsolve: <https://mrgsolve.org/> — the ODE simulation engine used by `sal_mrgsolve_model.R`.
- Shiny: <https://shiny.posit.co/> — the dashboard framework used by `sal_shiny_app.R`.
- Graphviz: <https://graphviz.org/> — used to render `sal_qsp_model.dot`.
- QSP in R: <https://vantage-research.net/qsp-in-r/>
- gPKPDviz (mrgsolve-based PK/PD Shiny tool): <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · <https://github.com/Genentech/gPKPDviz/>

---

## 단위 환산 (Unit conversions used throughout)

| 표기 | 환산 |
|------|------|
| mg/dL → mg/L | × 10 |
| mg/L → mmol/L | ÷ 138.12 (salicylate) |
| mg/dL → mmol/L | × 0.0724 |
| 치료 항염증 범위 | 150–300 mg/L (15–30 mg/dL) |
| 이명 역치 | 대략 200–300 mg/L |
| EXTRIP 투석 고려 (신기능 정상) | > 900 mg/L (90 mg/dL, 6.6 mmol/L) |
| 아스피린 → 살리실산 몰 환산 | × 138.12 / 180.16 = 0.767 |
