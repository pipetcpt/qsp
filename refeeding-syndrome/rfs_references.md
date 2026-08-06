# Refeeding Syndrome — References

References supporting the QSP model in this directory. They are grouped by the
part of the model each one constrains, and the annotation says **what the model
took from it** — because the calibration claim in `rfs_mrgsolve_model.R` is that
only three numbers were fitted to refeeding syndrome itself and everything else
came from normal physiology or non-refeeding pharmacology. These are the sources
for "everything else".

---

## 1. Definition, incidence and the guideline disagreement

1. Friedli N, Stanga Z, Sobotka L, et al. **Revisiting the refeeding syndrome: results of a systematic review.** *Nutrition* 2017;35:151–160. — <https://pubmed.ncbi.nlm.nih.gov/28241991/>
2. da Silva JSV, Seres DS, Sabino K, et al. **ASPEN Consensus Recommendations for Refeeding Syndrome.** *Nutr Clin Pract* 2020;35(2):178–195. — the 2020 consensus definition, risk criteria, and the explicit position that hypophosphataemia should be treated rather than used as a reason to withhold feeding. <https://pubmed.ncbi.nlm.nih.gov/32115791/>
3. National Institute for Health and Care Excellence. **Nutrition support for adults (CG32).** 2006, updated 2017. — the 10 kcal/kg/d starting rate, the 5 kcal/kg/d extreme-risk rate, and the history-based risk criteria the model reproduces in `rfs_sweep_history()`. <https://www.nice.org.uk/guidance/cg32>
4. Mehanna HM, Moledina J, Travis J. **Refeeding syndrome: what it is, and how to prevent and treat it.** *BMJ* 2008;336(7659):1495–1498. — <https://pubmed.ncbi.nlm.nih.gov/18583681/>
5. Friedli N, Stanga Z, Culkin A, et al. **Management and prevention of refeeding syndrome in medical inpatients: an evidence-based and consensus-supported algorithm.** *Nutrition* 2018;47:13–20. — <https://pubmed.ncbi.nlm.nih.gov/29429529/>
6. Doig GS, Simpson F, Heighes PT, et al. **Restricted versus continued standard caloric intake during the management of refeeding syndrome in critically ill adults: a randomised, parallel-group, multicentre, single-blind controlled trial.** *Lancet Respir Med* 2015;3(12):943–952. — the principal randomised evidence on caloric restriction; the comparison `rfs_sweep_gir()` and scenarios 01/02/07 are aimed at it. <https://pubmed.ncbi.nlm.nih.gov/26597128/>
7. Olthof LE, Koekkoek WACK, van Setten C, et al. **Impact of caloric intake in critically ill patients with, and without, refeeding syndrome: a retrospective study.** *Clin Nutr* 2018;37(5):1609–1617. — <https://pubmed.ncbi.nlm.nih.gov/28866138/>
8. Marvin VA, Brown D, Portlock J, Livingstone C. **Factors contributing to the development of hypophosphataemia when refeeding using parenteral nutrition.** *Pharm World Sci* 2008;30(4):329–335. — <https://pubmed.ncbi.nlm.nih.gov/18204995/>
9. Rio A, Whelan K, Goff L, et al. **Occurrence of refeeding syndrome in adults started on artificial nutrition support: prospective cohort study.** *BMJ Open* 2013;3(1):e002173. — incidence in an unselected artificial-nutrition cohort. <https://pubmed.ncbi.nlm.nih.gov/23315514/>
10. Crook MA, Hally V, Panteli JV. **The importance of the refeeding syndrome.** *Nutrition* 2001;17(7-8):632–637. — <https://pubmed.ncbi.nlm.nih.gov/11448586/>
11. Boateng AA, Sriram K, Meguid MM, Crook M. **Refeeding syndrome: treatment considerations based on collective analysis of literature case reports.** *Nutrition* 2010;26(2):156–167. — <https://pubmed.ncbi.nlm.nih.gov/19913393/>
12. Schnitker MA, Mattman PE, Bliss TL. **A clinical study of malnutrition in Japanese prisoners of war.** *Ann Intern Med* 1951;35(1):69–96. — the original description; the historical mortality the hazard scale is anchored against. <https://pubmed.ncbi.nlm.nih.gov/14847450/>

---

## 2. Body phosphorus: how much, and where — the model's central ratio

13. Bansal VK. **Serum Inorganic Phosphorus.** In: Walker HK, Hall WD, Hurst JW, eds. *Clinical Methods*, 3rd ed. Boston: Butterworths; 1990, ch. 198. — total body phosphorus ~700 g, the 85 % bone / 14 % intracellular / <1 % extracellular partition that gives the 0.06 % figure the whole model turns on. <https://www.ncbi.nlm.nih.gov/books/NBK311/>
14. Berndt T, Kumar R. **Phosphatonins and the regulation of phosphate homeostasis.** *Annu Rev Physiol* 2007;69:341–359. — <https://pubmed.ncbi.nlm.nih.gov/17002593/>
15. Institute of Medicine. **Dietary Reference Intakes for Calcium, Phosphorus, Magnesium, Vitamin D, and Fluoride.** Washington DC: National Academies Press; 1997. — dietary phosphorus intake, fractional absorption 0.55–0.70, and the balance the model's `FABSP` reproduces. <https://www.ncbi.nlm.nih.gov/books/NBK109825/>
16. Marks J, Debnam ES, Unwin RJ. **Phosphate homeostasis and the renal-gastrointestinal axis.** *Am J Physiol Renal Physiol* 2010;299(2):F285–F296. — <https://pubmed.ncbi.nlm.nih.gov/20534871/>
17. Sabbagh Y, O'Brien SP, Song W, et al. **Intestinal NaPi-IIb plays a major role in phosphate absorption and homeostasis.** *J Am Soc Nephrol* 2009;20(11):2348–2358. — <https://pubmed.ncbi.nlm.nih.gov/19729436/>
18. Forster IC, Hernando N, Biber J, Murer H. **Phosphate transporters of the SLC20 and SLC34 families.** *Mol Aspects Med* 2013;34(2-3):386–395. — the PiT-1/PiT-2 and NaPi-II kinetics behind `KMPUP`. <https://pubmed.ncbi.nlm.nih.gov/23506879/>
19. Payne RB. **Renal tubular reabsorption of phosphate (TmP/GFR): indications and interpretation.** *Ann Clin Biochem* 1998;35(Pt 2):201–206. — the threshold formulation `GFR × max(0, Pser − TmP/GFR)` and the 0.80–1.35 mmol/L normal range the model derives 0.96 for. <https://pubmed.ncbi.nlm.nih.gov/9547891/>
20. Walton RJ, Bijvoet OL. **Nomogram for derivation of renal threshold phosphate concentration.** *Lancet* 1975;2(7929):309–310. — <https://pubmed.ncbi.nlm.nih.gov/50513/>
21. Blaine J, Chonchol M, Levi M. **Renal control of calcium, phosphate, and magnesium homeostasis.** *Clin J Am Soc Nephrol* 2015;10(7):1257–1272. — <https://pubmed.ncbi.nlm.nih.gov/25287933/>

---

## 3. Hypophosphataemia: consequences, thresholds, and the energetic mechanism

22. Amanzadeh J, Reilly RF Jr. **Hypophosphatemia: an evidence-based approach to its clinical consequences and management.** *Nat Clin Pract Nephrol* 2006;2(3):136–148. — the <0.32 mmol/L severity threshold used for the cardiac and respiratory terms. <https://pubmed.ncbi.nlm.nih.gov/16932412/>
23. Knochel JP. **The pathophysiology and clinical characteristics of severe hypophosphatemia.** *Arch Intern Med* 1977;137(2):203–220. — <https://pubmed.ncbi.nlm.nih.gov/836118/>
24. O'Connor LR, Wheeler WS, Bethune JE. **Effect of hypophosphatemia on myocardial performance in man.** *N Engl J Med* 1977;297(17):901–903. — the myocardial contractility term. <https://pubmed.ncbi.nlm.nih.gov/909464/>
25. Aubier M, Murciano D, Lecocguic Y, et al. **Effect of hypophosphatemia on diaphragmatic contractility in patients with acute respiratory failure.** *N Engl J Med* 1985;313(7):420–424. — the diaphragm arm and the reason `RESULT 10` attributes most of the PaCO2 gap to the pump rather than the CO2 load. <https://pubmed.ncbi.nlm.nih.gov/4022081/>
26. Lichtman MA, Miller DR, Cohen J, Waterhouse C. **Reduced red cell glycolysis, 2,3-diphosphoglycerate and adenosine triphosphate concentration, and increased hemoglobin-oxygen affinity caused by hypophosphatemia.** *Ann Intern Med* 1971;74(4):562–568. — the 2,3-DPG compartment and the left-shifted oxygen curve. <https://pubmed.ncbi.nlm.nih.gov/5551159/>
27. Travis SF, Sugerman HJ, Ruberg RL, et al. **Alterations of red-cell glycolytic intermediates and oxygen transport as a consequence of hypophosphatemia in patients receiving intravenous hyperalimentation.** *N Engl J Med* 1971;285(14):763–768. — the original parenteral-nutrition observation. <https://pubmed.ncbi.nlm.nih.gov/5570340/>
28. Davis SV, Olichwier KK, Chakko SC. **Reversible depression of myocardial performance in hypophosphatemia.** *Am J Med Sci* 1988;295(3):183–187. — <https://pubmed.ncbi.nlm.nih.gov/3348135/>
29. Subramanian R, Khardori R. **Severe hypophosphatemia: pathophysiologic implications, clinical presentations, and treatment.** *Medicine (Baltimore)* 2000;79(1):1–8. — <https://pubmed.ncbi.nlm.nih.gov/10670405/>
30. Brautbar N, Baczynski R, Carpenter C, et al. **Impaired energy metabolism in rat myocardium during phosphate depletion.** *Am J Physiol* 1982;242(4):F699–F704. — the ATP/phosphocreatine dependence on inorganic phosphate. <https://pubmed.ncbi.nlm.nih.gov/7072654/>
31. Cirillo M, Ciacci C, De Santo NG. **Age, renal tubular phosphate reabsorption, and serum phosphate levels in adults.** *N Engl J Med* 2008;359(8):864–866. — <https://pubmed.ncbi.nlm.nih.gov/18716307/>
32. Boyd JW. **The relationships between blood haemoglobin concentration, packed cell volume and plasma protein concentration in dehydration.** *Br Vet J* 1981;137(2):166–172. — background for the haemoconcentration correction on measured concentrations.

---

## 4. Glucose, insulin, the incretin effect, and the transcellular shift

33. Bergman RN, Ider YZ, Bowden CR, Cobelli C. **Quantitative estimation of insulin sensitivity.** *Am J Physiol* 1979;236(6):E667–E677. — the minimal-model structure (remote insulin compartment X) used for the glucose–insulin submodel. <https://pubmed.ncbi.nlm.nih.gov/443421/>
34. Nauck MA, Homberger E, Siegel EG, et al. **Incretin effects of increasing glucose loads in man calculated from venous insulin and C-peptide responses.** *J Clin Endocrinol Metab* 1986;63(2):492–498. — the ~35 % incretin share of fed insulin secretion, and the oral-vs-intravenous asymmetry the model uses to separate enteral feed from IV dextrose. <https://pubmed.ncbi.nlm.nih.gov/3522621/>
35. Nauck MA, Meier JJ. **Incretin hormones: their role in health and disease.** *Diabetes Obes Metab* 2018;20(Suppl 1):5–21. — <https://pubmed.ncbi.nlm.nih.gov/29364588/>
36. DeFronzo RA, Felig P, Ferrannini E, Wahren J. **Effect of graded doses of insulin on splanchnic and peripheral potassium metabolism in man.** *Am J Physiol* 1980;238(5):E421–E427. — the insulin–potassium shift and its saturation, which set `EMAXK`. <https://pubmed.ncbi.nlm.nih.gov/6990783/>
37. Sterns RH, Grieff M, Bernstein PL. **Treatment of hyperkalemia: something old, something new.** *Kidney Int* 2016;89(3):546–554. — the 0.6–1.0 mmol/L fall in serum potassium from insulin-dextrose, used as a non-refeeding calibration point. <https://pubmed.ncbi.nlm.nih.gov/26880451/>
38. Newsholme EA, Start C. **Regulation in Metabolism.** London: Wiley; 1973. — glycolytic intermediate pools and the phosphate they sequester.
39. Klein S, Sakurai Y, Romijn JA, Carroll RM. **Progressive alterations in lipid and glucose metabolism during short-term fasting in young adult men.** *Am J Physiol* 1993;265(5 Pt 1):E801–E806. — <https://pubmed.ncbi.nlm.nih.gov/8238506/>
40. Soeters MR, Soeters PB, Schooneman MG, et al. **Adaptive reciprocity of lipid and glucose metabolism in human short-term starvation.** *Am J Physiol Endocrinol Metab* 2012;303(12):E1397–E1407. — <https://pubmed.ncbi.nlm.nih.gov/23074240/>
41. Newman JC, Verdin E. **Ketone bodies as signaling metabolites.** *Trends Endocrinol Metab* 2014;25(1):42–52. — <https://pubmed.ncbi.nlm.nih.gov/24140022/>
42. Boden G. **Effects of free fatty acids on gluconeogenesis and glycogenolysis.** *Life Sci* 2003;72(9):977–988. — <https://pubmed.ncbi.nlm.nih.gov/12495777/>

---

## 5. Starvation physiology, body composition and the metabolic adaptation

43. Keys A, Brozek J, Henschel A, Mickelsen O, Taylor HL. **The Biology of Human Starvation.** Minneapolis: University of Minnesota Press; 1950. — the Minnesota Starvation Experiment: the reference dataset for the starvation phase, including the body-composition and cardiac changes.
44. Cahill GF Jr. **Fuel metabolism in starvation.** *Annu Rev Nutr* 2006;26:1–22. — <https://pubmed.ncbi.nlm.nih.gov/16848698/>
45. Forbes GB. **Lean body mass-body fat interrelationships in humans.** *Nutr Rev* 1987;45(8):225–231. — the fat/lean partition rule used in the starvation phase. <https://pubmed.ncbi.nlm.nih.gov/3306482/>
46. Hall KD. **What is the required energy deficit per unit weight loss?** *Int J Obes (Lond)* 2008;32(3):573–576. — energy densities of adipose and lean tissue (9440 and 1800 kcal/kg). <https://pubmed.ncbi.nlm.nih.gov/17848938/>
47. Cunningham JJ. **A reanalysis of the factors influencing basal metabolic rate in normal adults.** *Am J Clin Nutr* 1980;33(11):2372–2374. — the REE equation. <https://pubmed.ncbi.nlm.nih.gov/7435418/>
48. Rosenbaum M, Leibel RL. **Adaptive thermogenesis in humans.** *Int J Obes (Lond)* 2010;34(Suppl 1):S47–S55. — the adaptive-thermogenesis floor and its slow recovery, reused as the marker of the insulin-resistant adapted state. <https://pubmed.ncbi.nlm.nih.gov/20935667/>
49. Pierson RN Jr, Wang J, Colt EW, Neumann P. **Body composition measurements in normal man: the potassium, sodium, sulfate and tritium spaces.** *J Chronic Dis* 1982;35(6):419–428. — the ~70 mmol total body potassium per kg fat-free mass the model uses (72). <https://pubmed.ncbi.nlm.nih.gov/7042769/>
50. Wang Z, St-Onge MP, Lecumberri B, et al. **Body cell mass: model development and validation at the cellular level of body composition.** *Am J Physiol Endocrinol Metab* 2004;286(1):E123–E128. — <https://pubmed.ncbi.nlm.nih.gov/14532167/>
51. Björntorp P. **Effects of energy reduction on metabolic rates and body composition.** *Int J Obes* 1981;5(Suppl 1):119–124.
52. Newton JL, Travis SPL. **Nutrition and the gut: refeeding.** *Clin Med (Lond)* 2006;6(2):141–146. — <https://pubmed.ncbi.nlm.nih.gov/16688969/>

---

## 6. Thiamine: stores, kinetics, transport, and Wernicke encephalopathy

53. Ariaey-Nejad MR, Balaghi M, Baker EM, Sauberlich HE. **Thiamin metabolism in man.** *Am J Clin Nutr* 1970;23(6):764–778. — the whole-body store (~25–30 mg) and the biological half-life (9–18 d) from which the model *predicts* the 1.1–1.4 mg/d requirement rather than assuming it. <https://pubmed.ncbi.nlm.nih.gov/5445068/>
54. Institute of Medicine. **Dietary Reference Intakes for Thiamin, Riboflavin, Niacin, Vitamin B6, Folate, Vitamin B12, Pantothenic Acid, Biotin, and Choline.** Washington DC: National Academies Press; 1998. — the published RDA the emergent check is scored against. <https://www.ncbi.nlm.nih.gov/books/NBK114310/>
55. Smithline HA, Donnino M, Greenblatt DJ. **Pharmacokinetics of high-dose oral thiamine hydrochloride in healthy subjects.** *BMC Clin Pharmacol* 2012;12:4. — the saturation of oral absorption, and the data behind the model's refuted prediction that oral cannot replete. <https://pubmed.ncbi.nlm.nih.gov/22305197/>
56. Thomson AD, Baker H, Leevy CM. **Patterns of 35S-thiamine hydrochloride absorption in the malnourished alcoholic patient.** *J Lab Clin Med* 1970;76(1):34–45. — ethanol impairment of intestinal thiamine uptake. <https://pubmed.ncbi.nlm.nih.gov/5449328/>
57. Said HM. **Intestinal absorption of water-soluble vitamins in health and disease.** *Biochem J* 2011;437(3):357–372. — ThTR-1/ThTR-2 kinetics and the passive route at high concentrations, which is the mechanistic basis of the "two clocks" result. <https://pubmed.ncbi.nlm.nih.gov/21749321/>
58. Rindi G, Laforenza U. **Thiamine intestinal transport and related issues: recent aspects.** *Proc Soc Exp Biol Med* 2000;224(4):246–255. — <https://pubmed.ncbi.nlm.nih.gov/10964259/>
59. Donnino MW, Vega J, Miller J, Walsh M. **Myths and misconceptions of Wernicke's encephalopathy: what every emergency physician should know.** *Ann Emerg Med* 2007;50(6):715–721. — including the "thiamine before glucose" question the timing sweep addresses directly. <https://pubmed.ncbi.nlm.nih.gov/17681641/>
60. Schabelman E, Kuo D. **Glucose before thiamine for Wernicke encephalopathy: a literature review.** *J Emerg Med* 2012;42(4):488–494. — the evidence base for the rule the model only partly supports. <https://pubmed.ncbi.nlm.nih.gov/22104258/>
61. Galvin R, Bråthen G, Ivashynka A, et al. **EFNS guidelines for diagnosis, therapy and prevention of Wernicke encephalopathy.** *Eur J Neurol* 2010;17(12):1408–1418. — the 200–500 mg IV dosing used in the thiamine arms. <https://pubmed.ncbi.nlm.nih.gov/20642790/>
62. Sechi G, Serra A. **Wernicke's encephalopathy: new clinical settings and recent advances in diagnosis and management.** *Lancet Neurol* 2007;6(5):442–455. — <https://pubmed.ncbi.nlm.nih.gov/17434099/>
63. Klein M, Weksler N, Gurman GM. **Fatal metabolic acidosis caused by thiamine deficiency.** *J Emerg Med* 2004;26(3):301–303. — type B lactic acidosis from a blocked pyruvate gate. <https://pubmed.ncbi.nlm.nih.gov/15028328/>
64. Attaluri P, Castillo A, Edriss H, Nugent K. **Thiamine deficiency: an important consideration in critically ill patients.** *Am J Med Sci* 2018;356(4):382–390. — <https://pubmed.ncbi.nlm.nih.gov/30360807/>
65. Patel MS, Korotchkina LG. **Regulation of the pyruvate dehydrogenase complex.** *Biochem Soc Trans* 2006;34(Pt 2):217–222. — TPP dependence of the PDH gate. <https://pubmed.ncbi.nlm.nih.gov/16545080/>
66. Frank LL. **Thiamin in clinical practice.** *JPEN J Parenter Enteral Nutr* 2015;39(5):503–520. — <https://pubmed.ncbi.nlm.nih.gov/25564426/>

---

## 7. Potassium and magnesium: renal handling and the coupling between them

67. Youn JH, McDonough AA. **Recent advances in understanding integrative control of potassium homeostasis.** *Annu Rev Physiol* 2009;71:381–401. — <https://pubmed.ncbi.nlm.nih.gov/18759644/>
68. Palmer BF. **Regulation of potassium homeostasis.** *Clin J Am Soc Nephrol* 2015;10(6):1050–1060. — the aldosterone/ROMK control the model uses, including potassium as the dominant aldosterone secretagogue. <https://pubmed.ncbi.nlm.nih.gov/24721891/>
69. Squires RD, Huth EJ. **Experimental potassium depletion in normal human subjects. I. Relation of ionic intakes to the renal conservation of potassium.** *J Clin Invest* 1959;38(7):1134–1148. — the incomplete renal adaptation to potassium depletion (a floor of roughly 10–20 mmol/d) that the model's `KOBL` and ROMK term reproduce; the absence of this was one of the defects the Python verification exposed. <https://pubmed.ncbi.nlm.nih.gov/13664789/>
70. Huang CL, Kuo E. **Mechanism of hypokalemia in magnesium deficiency.** *J Am Soc Nephrol* 2007;18(10):2649–2652. — the intracellular magnesium block on ROMK, i.e. the mechanism behind refractory hypokalaemia. <https://pubmed.ncbi.nlm.nih.gov/17804670/>
71. Whang R, Whang DD, Ryan MP. **Refractory potassium repletion. A consequence of magnesium deficiency.** *Arch Intern Med* 1992;152(1):40–45. — the clinical phenomenon the model reproduces only modestly, and says so. <https://pubmed.ncbi.nlm.nih.gov/1728927/>
72. Agus ZS. **Hypomagnesemia.** *J Am Soc Nephrol* 1999;10(7):1616–1622. — <https://pubmed.ncbi.nlm.nih.gov/10405219/>
73. de Baaij JHF, Hoenderop JGJ, Bindels RJM. **Magnesium in man: implications for health and disease.** *Physiol Rev* 2015;95(1):1–46. — total body magnesium, the exchangeable fraction, and the renal threshold that makes a bolus wasteful. <https://pubmed.ncbi.nlm.nih.gov/25540137/>
74. Elin RJ. **Assessment of magnesium status for diagnosis and therapy.** *Magnes Res* 2010;23(4):S194–S198. — <https://pubmed.ncbi.nlm.nih.gov/20736141/>
75. Rude RK, Oldham SB, Singer FR. **Functional hypoparathyroidism and parathyroid hormone end-organ resistance in human magnesium deficiency.** *Clin Endocrinol (Oxf)* 1976;5(3):209–224. — magnesium dependence of PTH secretion, which in this model is what disables the skeletal phosphate supply line. <https://pubmed.ncbi.nlm.nih.gov/947478/>
76. Vormann J. **Magnesium: nutrition and metabolism.** *Mol Aspects Med* 2003;24(1-3):27–37. — <https://pubmed.ncbi.nlm.nih.gov/12537987/>

---

## 8. Calcium, PTH, vitamin D and FGF23

77. Brown EM. **Role of the calcium-sensing receptor in extracellular calcium homeostasis.** *Best Pract Res Clin Endocrinol Metab* 2013;27(3):333–343. — the sigmoidal PTH–calcium relationship. <https://pubmed.ncbi.nlm.nih.gov/23856263/>
78. Shimada T, Hasegawa H, Yamazaki Y, et al. **FGF-23 is a potent regulator of vitamin D metabolism and phosphate homeostasis.** *J Bone Miner Res* 2004;19(3):429–435. — <https://pubmed.ncbi.nlm.nih.gov/15040831/>
79. Kuro-o M. **The FGF23 and Klotho system beyond mineral metabolism.** *Clin Exp Nephrol* 2017;21(Suppl 1):64–69. — <https://pubmed.ncbi.nlm.nih.gov/28062938/>
80. Block GA, Hulbert-Shearon TE, Levin NW, Port FK. **Association of serum phosphorus and calcium × phosphate product with mortality risk in chronic hemodialysis patients.** *Am J Kidney Dis* 1998;31(4):607–617. — the calcium–phosphate product threshold used in the precipitation term. <https://pubmed.ncbi.nlm.nih.gov/9531176/>
81. Bushinsky DA, Monk RD. **Electrolyte quintet: calcium.** *Lancet* 1998;352(9124):306–311. — <https://pubmed.ncbi.nlm.nih.gov/9690425/>

---

## 9. Repletion pharmacology: doses, rates and their hazards

82. Taylor BE, Huey WY, Buchman TG, et al. **Treatment of hypophosphatemia using a protocol based on patient weight and serum phosphorus level in a surgical intensive care unit.** *J Am Coll Surg* 2004;198(2):198–204. — the weight-based repletion protocol the model's `P_DOSE` band is taken from. <https://pubmed.ncbi.nlm.nih.gov/14759775/>
83. Charron T, Bernard F, Skrobik Y, et al. **Intravenous phosphate in the intensive care unit: more aggressive repletion regimens for moderate and severe hypophosphatemia.** *Intensive Care Med* 2003;29(8):1273–1278. — <https://pubmed.ncbi.nlm.nih.gov/12845429/>
84. Brown KA, Dickerson RN, Morgan LM, et al. **A new graduated dosing regimen for phosphorus replacement in patients receiving nutrition support.** *JPEN J Parenter Enteral Nutr* 2006;30(3):209–214. — <https://pubmed.ncbi.nlm.nih.gov/16639067/>
85. Bugg NC, Jones JA. **Hypophosphataemia: pathophysiology, effects and management on the intensive care unit.** *Anaesthesia* 1998;53(9):895–902. — <https://pubmed.ncbi.nlm.nih.gov/9849285/>
86. Kraft MD, Btaiche IF, Sacks GS, Kudsk KA. **Treatment of electrolyte disorders in adult patients in the intensive care unit.** *Am J Health Syst Pharm* 2005;62(16):1663–1682. — potassium, magnesium and phosphate repletion rates. <https://pubmed.ncbi.nlm.nih.gov/16085929/>
87. Hébert P, Mehta N, Wang J, et al. **Functional magnesium deficiency in critically ill patients identified using a magnesium-loading test.** *Crit Care Med* 1997;25(5):749–755. — the renal-threshold behaviour behind the bolus-versus-infusion result. <https://pubmed.ncbi.nlm.nih.gov/9187592/>
88. Ayuk J, Gittoes NJL. **How should hypomagnesaemia be investigated and treated?** *Clin Endocrinol (Oxf)* 2011;75(6):743–746. — <https://pubmed.ncbi.nlm.nih.gov/21569071/>

---

## 10. Route, formula composition and the respiratory cost of carbohydrate

89. Askanazi J, Rosenbaum SH, Hyman AI, et al. **Respiratory changes induced by the large glucose loads of total parenteral nutrition.** *JAMA* 1980;243(14):1444–1447. — the CO2 load of glucose-based feeding. <https://pubmed.ncbi.nlm.nih.gov/6767014/>
90. Talpers SS, Romberger DJ, Bunce SB, Pingleton SK. **Nutritionally associated increased carbon dioxide production. Excess total calories vs high proportion of carbohydrate calories.** *Chest* 1992;102(2):551–555. — the finding that total calories matter more than carbohydrate fraction for CO2, which the model's RESULT 10 revisits and partly reinterprets. <https://pubmed.ncbi.nlm.nih.gov/1643945/>
91. Wolfe RR, O'Donnell TF Jr, Stone MD, et al. **Investigation of factors determining the optimal glucose infusion rate in total parenteral nutrition.** *Metabolism* 1980;29(9):892–900. — the ~4–5 mg/kg/min ceiling on glucose oxidation used in the gas-exchange block. <https://pubmed.ncbi.nlm.nih.gov/6774136/>
92. Singer P, Blaser AR, Berger MM, et al. **ESPEN guideline on clinical nutrition in the intensive care unit.** *Clin Nutr* 2019;38(1):48–79. — <https://pubmed.ncbi.nlm.nih.gov/30348463/>
93. Elia M, Cummings JH. **Physiological aspects of energy metabolism and gastrointestinal effects of carbohydrates.** *Eur J Clin Nutr* 2007;61(Suppl 1):S40–S74. — <https://pubmed.ncbi.nlm.nih.gov/17992186/>

---

## 11. Population-specific evidence: anorexia nervosa, alcohol, and critical illness

94. Garber AK, Cheng J, Accurso EC, et al. **Short-term outcomes of the study of refeeding to optimize inpatient gains for patients with anorexia nervosa: a randomized clinical trial.** *JAMA Pediatr* 2021;175(1):19–27. — the higher-calorie refeeding trial in anorexia nervosa. <https://pubmed.ncbi.nlm.nih.gov/33074282/>
95. O'Connor G, Nicholls D. **Refeeding hypophosphatemia in adolescents with anorexia nervosa: a systematic review.** *Nutr Clin Pract* 2013;28(3):358–364. — <https://pubmed.ncbi.nlm.nih.gov/23515708/>
96. Whitelaw M, Gilbertson H, Lam PY, Sawyer SM. **Does aggressive refeeding in hospitalized adolescents with anorexia nervosa result in increased hypophosphatemia?** *J Adolesc Health* 2010;46(6):577–582. — <https://pubmed.ncbi.nlm.nih.gov/20472215/>
97. Society for Adolescent Health and Medicine. **Refeeding hypophosphatemia in hospitalized adolescents with anorexia nervosa: a position statement.** *J Adolesc Health* 2022;71(4):517–520. — <https://pubmed.ncbi.nlm.nih.gov/35970149/>
98. Sacks GS, Walker J, Dickerson RN, et al. **Observations of hypophosphatemia and its management in nutrition support.** *Nutr Clin Pract* 1994;9(3):105–108. — <https://pubmed.ncbi.nlm.nih.gov/8065226/>
99. Palmese S, Scarano V, Bruno G. **Hypophosphatemia and metabolic acidosis: a case of severe refeeding syndrome.** *Nutr Ther Metab* 2010;28(2):89–92.
100. Elnenaei MO, Alaghband-Zadeh J, Sherwood R, et al. **Leptin and insulin growth factor 1: diagnostic markers of the refeeding syndrome and mortality.** *Br J Nutr* 2011;106(6):906–912. — <https://pubmed.ncbi.nlm.nih.gov/21736847/>
101. Pourhassan M, Cuvelier I, Gehrke I, et al. **Risk factors of refeeding syndrome in malnourished older hospitalized patients.** *Clin Nutr* 2018;37(4):1354–1359. — <https://pubmed.ncbi.nlm.nih.gov/28647294/>

---

## 12. QSP, mrgsolve and modelling method

102. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn* 2015;42:S84–S85. — <https://mrgsolve.org>
103. Peterson MC, Riggs MM. **A physiologically based mathematical model of integrated calcium homeostasis and bone remodeling.** *Bone* 2010;46(1):49–63. — the reference QSP treatment of the calcium–phosphate–PTH–vitamin D system, and the structural model for this file's bone block. <https://pubmed.ncbi.nlm.nih.gov/19664732/>
104. Riggs MM, Peterson MC, Gastonguay MR. **Multiscale physiology-based modeling of mineral bone disorder in patients with impaired kidney function.** *J Clin Pharmacol* 2012;52(1 Suppl):45S–53S. — <https://pubmed.ncbi.nlm.nih.gov/22232755/>
105. Nijhout HF, Best JA, Reed MC. **Systems biology of robustness and homeostatic mechanisms.** *Wiley Interdiscip Rev Syst Biol Med* 2019;11(3):e1440. — on why a model whose healthy state is an exact steady state behaves differently from one that is merely initialised near it. <https://pubmed.ncbi.nlm.nih.gov/30576081/>
106. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2016;5(5):235–249. — <https://pubmed.ncbi.nlm.nih.gov/27299936/>

---

### How these were used

The calibration policy stated in `rfs_mrgsolve_model.R` is that **three numbers**
were fitted to refeeding syndrome itself: the size of the flux-driven rise in the
cellular organic-phosphate set-point, the rate at which cells approach it, and one
global mortality-hazard scale. Sections 2, 5, 6, 7 and 8 above supply the normal
physiology; sections 4, 9 and 10 supply non-refeeding pharmacology; section 1 and
section 11 supply the clinical outcomes the model's predictions are **compared
against** rather than fitted to.

The strongest internal check comes from section 6: references 53 and 54 give the
whole-body thiamine store and its half-life independently of any requirement
figure, and the model then computes a dietary requirement of **1.42 mg/d** against
a published RDA of **1.1–1.4 mg/d** — a number it was never shown.
