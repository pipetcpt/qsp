# Anorexia Nervosa (AN) — References

Curated bibliography supporting the mechanistic map (`an_qsp_model.dot`), the
mrgsolve QSP model (`an_mrgsolve_model.R`) and the Shiny dashboard
(`an_shiny_app.R`) in this directory.

Every numbered entry below was resolved against PubMed and links to a
**verified PMID**; the two entries without a PMID are books/monographs and link
to a PubMed search instead. Organized by the mechanistic cluster each source
supports.

---

## 1. Epidemiology, Natural History & Mortality

1. Arcelus J, Mitchell AJ, Wales J, et al. Mortality rates in patients with anorexia nervosa and other eating disorders. A meta-analysis of 36 studies. *Arch Gen Psychiatry*. 2011;68:724-31. [PMID 21727255](https://pubmed.ncbi.nlm.nih.gov/21727255/) — SMR 5.86; the anchor for the model's mortality framing.
2. Zipfel S, Giel KE, Bulik CM, et al. Anorexia nervosa: aetiology, assessment, and treatment. *Lancet Psychiatry*. 2015;2:1099-111. [PMID 26514083](https://pubmed.ncbi.nlm.nih.gov/26514083/)
3. Treasure J, Zipfel S, Micali N, et al. Anorexia nervosa. *Nat Rev Dis Primers*. 2015;1:15074. [PMID 27189821](https://pubmed.ncbi.nlm.nih.gov/27189821/)
4. Treasure J, Duarte TA, Schmidt U. Eating disorders. *Lancet*. 2020;395:899-911. [PMID 32171414](https://pubmed.ncbi.nlm.nih.gov/32171414/)
5. van Eeden AE, van Hoeken D, Hoek HW. Incidence, prevalence and mortality of anorexia nervosa and bulimia nervosa. *Curr Opin Psychiatry*. 2021;34:515-524. [PMID 34419970](https://pubmed.ncbi.nlm.nih.gov/34419970/)

## 2. Genetic & Metabo-Psychiatric Architecture (map cluster 1)

6. Watson HJ, Yilmaz Z, Thornton LM, et al. Genome-wide association study identifies eight risk loci and implicates metabo-psychiatric origins for anorexia nervosa. *Nat Genet*. 2019;51:1207-1214. [PMID 31308545](https://pubmed.ncbi.nlm.nih.gov/31308545/) — the basis for treating AN as partly metabolic, not purely psychiatric.
7. Duncan L, Yilmaz Z, Gaspar H, et al. Significant locus and metabolic genetic correlations revealed in genome-wide association study of anorexia nervosa. *Am J Psychiatry*. 2017;174:850-858. [PMID 28494655](https://pubmed.ncbi.nlm.nih.gov/28494655/)
8. Hübel C, Gaspar HA, Coleman JRI, et al. Genomics of body fat percentage may contribute to sex bias in anorexia nervosa. *Am J Med Genet B Neuropsychiatr Genet*. 2019;180:428-438. [PMID 30593698](https://pubmed.ncbi.nlm.nih.gov/30593698/)
9. Bulik CM, Blake L, Austin J. Genetics of eating disorders: what the clinician needs to know. *Psychiatr Clin North Am*. 2019;42:59-73. [PMID 30704640](https://pubmed.ncbi.nlm.nih.gov/30704640/)

## 3. Neurocircuitry, Habit & Neurotransmitters (map clusters 2-3)

10. Kaye WH, Fudge JL, Paulus M. New insights into symptoms and neurocircuit function of anorexia nervosa. *Nat Rev Neurosci*. 2009;10:573-84. [PMID 19603056](https://pubmed.ncbi.nlm.nih.gov/19603056/)
11. Foerde K, Steinglass JE, Shohamy D, et al. Neural mechanisms supporting maladaptive food choices in anorexia nervosa. *Nat Neurosci*. 2015;18:1571-3. [PMID 26457555](https://pubmed.ncbi.nlm.nih.gov/26457555/) — dorsal-striatal engagement during food choice, the empirical basis for the habit node.
12. Steinglass JE, Walsh BT. Neurobiological model of the persistence of anorexia nervosa. *J Eat Disord*. 2016;4:19. [PMID 27195123](https://pubmed.ncbi.nlm.nih.gov/27195123/)
13. Uniacke B, Walsh BT, Foerde K, et al. The role of habits in anorexia nervosa: where we are and where to go from here? *Curr Psychiatry Rep*. 2018;20:61. [PMID 30039342](https://pubmed.ncbi.nlm.nih.gov/30039342/)
14. Frank GKW, Shott ME, DeGuzman MC. The neurobiology of eating disorders. *Child Adolesc Psychiatr Clin N Am*. 2019;28:629-640. [PMID 31443880](https://pubmed.ncbi.nlm.nih.gov/31443880/)
15. Bailer UF, Frank GK, Henry SE, et al. Altered brain serotonin 5-HT1A receptor binding after recovery from anorexia nervosa measured by positron emission tomography. *Arch Gen Psychiatry*. 2005;62:1032-41. [PMID 16143735](https://pubmed.ncbi.nlm.nih.gov/16143735/)

## 4. Starvation Physiology, Energy Balance & Body-Composition Modeling (map cluster 6; model FM/FFM/ADAPT)

16. Keys A, Brožek J, Henschel A, Mickelsen O, Taylor HL. *The Biology of Human Starvation*. University of Minnesota Press; 1950. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Minnesota+semi-starvation+experiment+Keys+biology+of+human+starvation) — semi-starvation in healthy volunteers reproduced food preoccupation, ritualization and cognitive rigidity; the source of the model's starvation-perpetuation term.
17. Hall KD, Sacks G, Chandramohan D, et al. Quantification of the effect of energy imbalance on bodyweight. *Lancet*. 2011;378:826-37. [PMID 21872751](https://pubmed.ncbi.nlm.nih.gov/21872751/) — energy-partition framework behind the FM/FFM ODEs.
18. Hall KD. Body fat and fat-free mass inter-relationships: Forbes's theory revisited. *Br J Nutr*. 2007;97:1059-63. [PMID 17367567](https://pubmed.ncbi.nlm.nih.gov/17367567/) — the dFFM/dFM = C/FM rule with C ≈ 10.4 kg implemented in `[ODE]` section 5.
19. Dulloo AG, Jacquet J, Montani JP, et al. How dieting makes the lean fatter: from a perspective of body composition autoregulation. *Obes Rev*. 2015;16 Suppl 1:25-35. [PMID 25614201](https://pubmed.ncbi.nlm.nih.gov/25614201/) — adaptive thermogenesis (`ADAPT`, `ADAPT_MAX`).
20. Krahn DD, Rock C, Dechert RE, et al. Changes in resting energy expenditure and body composition in anorexia nervosa patients during refeeding. *J Am Diet Assoc*. 1993;93:434-8. [PMID 8454812](https://pubmed.ncbi.nlm.nih.gov/8454812/) — refeeding hypermetabolism (`K_HYPERMET`).
21. Marzola E, Nasser JA, Hashim SA, et al. Nutritional rehabilitation in anorexia nervosa: review of the literature and implications for treatment. *BMC Psychiatry*. 2013;13:290. [PMID 24200367](https://pubmed.ncbi.nlm.nih.gov/24200367/)
22. Schebendach J, Mayer LE, Devlin MJ, et al. Dietary energy density and diet variety as risk factors for relapse in anorexia nervosa: a replication. *Int J Eat Disord*. 2012;45:79-84. [PMID 21448937](https://pubmed.ncbi.nlm.nih.gov/21448937/)

## 5. Neuroendocrine Adaptation — Leptin, Ghrelin, Thyroid, HPA, GH/IGF-1 (map clusters 4, 7, 8, 10, 11)

23. Misra M, Klibanski A. Endocrine consequences of anorexia nervosa. *Lancet Diabetes Endocrinol*. 2014;2:581-92. [PMID 24731664](https://pubmed.ncbi.nlm.nih.gov/24731664/) — the umbrella reference for the model's endocrine block.
24. Schorr M, Miller KK. The endocrine manifestations of anorexia nervosa: mechanisms and management. *Nat Rev Endocrinol*. 2017;13:174-186. [PMID 27811940](https://pubmed.ncbi.nlm.nih.gov/27811940/)
25. Misra M, Miller KK, Kuo K, et al. Secretory dynamics of leptin in adolescent girls with anorexia nervosa and healthy adolescents. *Am J Physiol Endocrinol Metab*. 2005;289:E373-81. [PMID 15811876](https://pubmed.ncbi.nlm.nih.gov/15811876/)
26. Germain N, Galusca B, Grouselle D, et al. Ghrelin/obestatin ratio in two populations with low bodyweight: constitutional thinness and anorexia nervosa. *Psychoneuroendocrinology*. 2009;34:413-9. [PMID 18995969](https://pubmed.ncbi.nlm.nih.gov/18995969/)
27. Croxson MS, Ibbertson HK. Low serum triiodothyronine (T3) and hypothyroidism in anorexia nervosa. *J Clin Endocrinol Metab*. 1977;44:167-74. [PMID 401822](https://pubmed.ncbi.nlm.nih.gov/401822/) — the low-T3 (euthyroid-sick) state reproduced by `T3_SUPP`.
28. Misra M, Miller KK, Almazan C, et al. Alterations in cortisol secretory dynamics in adolescent girls with anorexia nervosa and effects on bone metabolism. *J Clin Endocrinol Metab*. 2004;89:4972-80. [PMID 15472193](https://pubmed.ncbi.nlm.nih.gov/15472193/)
29. Lawson EA, Donoho D, Miller KK, et al. Hypercortisolemia is associated with severity of bone loss and depression in hypothalamic amenorrhea and anorexia nervosa. *J Clin Endocrinol Metab*. 2009;94:4710-6. [PMID 19837921](https://pubmed.ncbi.nlm.nih.gov/19837921/) — the cortisol → bone arm (`K_P1NP_CORT`, `K_CTX_CORT`).
30. Misra M, Miller KK, Bjornson J, et al. Alterations in growth hormone secretory dynamics in adolescent girls with anorexia nervosa and effects on bone metabolism. *J Clin Endocrinol Metab*. 2003;88:5615-23. [PMID 14671143](https://pubmed.ncbi.nlm.nih.gov/14671143/) — GH resistance with low IGF-1 (`IGF_SUPP`).
31. Modan-Moses D, Yaroslavsky A, Kochavi B, et al. Linear growth and final height characteristics in adolescent females with anorexia nervosa. *PLoS One*. 2012;7:e45504. [PMID 23029058](https://pubmed.ncbi.nlm.nih.gov/23029058/)

## 6. Reproductive Axis & Functional Hypothalamic Amenorrhea (map cluster 9)

32. Köpp W, Blum WF, von Prittwitz S, et al. Low leptin levels predict amenorrhea in underweight and eating disordered females. *Mol Psychiatry*. 1997;2:335-40. [PMID 9246675](https://pubmed.ncbi.nlm.nih.gov/9246675/)
33. Audi L, Mantzoros CS, Vidal-Puig A, et al. Leptin in relation to resumption of menses in women with anorexia nervosa. *Mol Psychiatry*. 1998;3:544-7. [PMID 9857982](https://pubmed.ncbi.nlm.nih.gov/9857982/) — source of the ~1.85 ng/mL leptin threshold used as `LEP50_GNRH`.

## 7. Bone — Low Formation / High Resorption & Bone-Directed Therapy (map clusters 12, 19)

34. Misra M, Klibanski A. Anorexia nervosa and bone. *J Endocrinol*. 2014;221:R163-76. [PMID 24898127](https://pubmed.ncbi.nlm.nih.gov/24898127/)
35. Fazeli PK, Klibanski A. Bone metabolism in anorexia nervosa. *Curr Osteoporos Rep*. 2014;12:82-9. [PMID 24419863](https://pubmed.ncbi.nlm.nih.gov/24419863/) — the uncoupled low-P1NP/high-CTX phenotype the BMD ODE encodes.
36. Misra M, Katzman D, Miller KK, et al. Physiologic estrogen replacement increases bone density in adolescent girls with anorexia nervosa. *J Bone Miner Res*. 2011;26:2430-8. [PMID 21698665](https://pubmed.ncbi.nlm.nih.gov/21698665/) — transdermal, not oral, estradiol; the `E2_PATCH` arm.
37. Grinspoon S, Thomas L, Miller K, et al. Effects of recombinant human IGF-I and oral contraceptive administration on bone density in anorexia nervosa. *J Clin Endocrinol Metab*. 2002;87:2883-91. [PMID 12050268](https://pubmed.ncbi.nlm.nih.gov/12050268/) — the IGF-1 → osteoblast arm.
38. Fazeli PK, Wang IS, Miller KK, et al. Teriparatide increases bone formation and bone mineral density in adult women with anorexia nervosa. *J Clin Endocrinol Metab*. 2014;99:1322-9. [PMID 24456286](https://pubmed.ncbi.nlm.nih.gov/24456286/) — the `TPTD` arm (`EMAX_TPTD`).
39. Miller KK, Meenaghan E, Lawson EA, et al. Effects of risedronate and low-dose transdermal testosterone on bone mineral density in women with anorexia nervosa: a randomized, placebo-controlled study. *J Clin Endocrinol Metab*. 2011;96:2081-8. [PMID 21525157](https://pubmed.ncbi.nlm.nih.gov/21525157/) — the `BIS` arm (`EMAX_BIS`).
40. Bredella MA, Fazeli PK, Daley SM, et al. Marrow fat composition in anorexia nervosa. *Bone*. 2014;66:199-204. [PMID 24953711](https://pubmed.ncbi.nlm.nih.gov/24953711/) — paradoxical marrow adiposity behind `K_P1NP_LEP`.

## 8. Cardiovascular & Autonomic Complications (map cluster 13)

41. Mehler PS, Brown C. Anorexia nervosa — medical complications. *J Eat Disord*. 2015;3:11. [PMID 25834735](https://pubmed.ncbi.nlm.nih.gov/25834735/)
42. Mehler PS, Krantz MJ, Sachs KV. Treatments of medical complications of anorexia nervosa and bulimia nervosa. *J Eat Disord*. 2015;3:15. [PMID 25874112](https://pubmed.ncbi.nlm.nih.gov/25874112/)
43. Sachs KV, Harnke B, Mehler PS, et al. Cardiovascular complications of anorexia nervosa: a systematic review. *Int J Eat Disord*. 2016;49:238-48. [PMID 26710932](https://pubmed.ncbi.nlm.nih.gov/26710932/) — bradycardia/QTc behavior of the `HR` compartment.
44. Casiero D, Frishman WH. Cardiovascular complications of eating disorders. *Cardiol Rev*. 2006;14:227-31. [PMID 16924163](https://pubmed.ncbi.nlm.nih.gov/16924163/)
45. Giovinazzo S, Sukkar SG, Rosa GM, et al. Anorexia nervosa and heart disease: a systematic review. *Eat Weight Disord*. 2019;24:199-207. [PMID 30173377](https://pubmed.ncbi.nlm.nih.gov/30173377/)

## 9. Gastrointestinal Physiology, Gut Peptides & Microbiome (map cluster 5)

46. Norris ML, Harrison ME, Isserlin L, et al. Gastrointestinal complications associated with anorexia nervosa: a systematic review. *Int J Eat Disord*. 2016;49:216-37. [PMID 26407541](https://pubmed.ncbi.nlm.nih.gov/26407541/)
47. Benini L, Todesco T, Dalle Grave R, et al. Gastric emptying in patients with restricting and binge/purging subtypes of anorexia nervosa. *Am J Gastroenterol*. 2004;99:1448-54. [PMID 15307858](https://pubmed.ncbi.nlm.nih.gov/15307858/) — delayed emptying behind `GI_SAT_MAX`.
48. Kleiman SC, Watson HJ, Bulik-Sullivan EC, et al. The intestinal microbiota in acute anorexia nervosa and during renourishment. *Psychosom Med*. 2015;77:969-81. [PMID 26428446](https://pubmed.ncbi.nlm.nih.gov/26428446/)
49. Fetissov SO. Role of the gut microbiota in host appetite control: bacterial growth to animal feeding behaviour. *Nat Rev Endocrinol*. 2017;13:11-25. [PMID 27616451](https://pubmed.ncbi.nlm.nih.gov/27616451/)

## 10. Refeeding Syndrome & Nutritional Rehabilitation Strategy (map cluster 15; model PHOS/POT/MG/THIA)

50. Crook MA, Hally V, Panteli JV. The importance of the refeeding syndrome. *Nutrition*. 2001;17:632-7. [PMID 11448586](https://pubmed.ncbi.nlm.nih.gov/11448586/)
51. Whitelaw M, Gilbertson H, Lam PY, et al. Does aggressive refeeding in hospitalized adolescents with anorexia nervosa result in increased hypophosphatemia? *J Adolesc Health*. 2010;46:577-82. [PMID 20472215](https://pubmed.ncbi.nlm.nih.gov/20472215/)
52. Golden NH, Keane-Miller C, Sainani KL, et al. Higher caloric intake in hospitalized adolescents with anorexia nervosa is associated with reduced length of stay and no increased rate of refeeding syndrome. *J Adolesc Health*. 2013;53:573-8. [PMID 23830088](https://pubmed.ncbi.nlm.nih.gov/23830088/)
53. O'Connor G, Nicholls D, Hudson L, et al. Refeeding low weight hospitalized adolescents with anorexia nervosa: a multicenter randomized controlled trial. *Nutr Clin Pract*. 2016;31:681-9. [PMID 26869609](https://pubmed.ncbi.nlm.nih.gov/26869609/)
54. Friedli N, Stanga Z, Sobotka L, et al. Revisiting the refeeding syndrome: results of a systematic review. *Nutrition*. 2017;35:151-160. [PMID 28087222](https://pubmed.ncbi.nlm.nih.gov/28087222/)
55. da Silva JSV, Seres DS, Sabino K, et al. ASPEN consensus recommendations for refeeding syndrome. *Nutr Clin Pract*. 2020;35:178-195. [PMID 32115791](https://pubmed.ncbi.nlm.nih.gov/32115791/) — definition of the phosphate/potassium/magnesium thresholds used in the composite risk index.
56. Garber AK, Cheng J, Accurso EC, et al. Short-term outcomes of the Study of Refeeding to Optimize Inpatient Gains for patients with anorexia nervosa (StRONG): a multicenter randomized clinical trial. *JAMA Pediatr*. 2021;175:19-27. [PMID 33074282](https://pubmed.ncbi.nlm.nih.gov/33074282/) — higher-calorie refeeding restored weight faster with no excess electrolyte events; scenarios 2 vs 3.

## 11. Pharmacotherapy (map clusters 17-18)

57. Attia E, Haiman C, Walsh BT, et al. Does fluoxetine augment the inpatient treatment of anorexia nervosa? *Am J Psychiatry*. 1998;155:548-51. [PMID 9546003](https://pubmed.ncbi.nlm.nih.gov/9546003/) — no benefit in the underweight state.
58. Walsh BT, Kaplan AS, Attia E, et al. Fluoxetine after weight restoration in anorexia nervosa: a randomized controlled trial. *JAMA*. 2006;295:2605-12. [PMID 16772623](https://pubmed.ncbi.nlm.nih.gov/16772623/) — together with #57, the evidence encoded as the model's nutritional gate on SSRI efficacy.
59. Bissada H, Tasca GA, Barber AM, et al. Olanzapine in the treatment of low body weight and obsessive thinking in women with anorexia nervosa: a randomized, double-blind, placebo-controlled trial. *Am J Psychiatry*. 2008;165:1281-8. [PMID 18558642](https://pubmed.ncbi.nlm.nih.gov/18558642/)
60. Attia E, Steinglass JE, Walsh BT, et al. Olanzapine versus placebo in adult outpatients with anorexia nervosa: a randomized clinical trial. *Am J Psychiatry*. 2019;176:449-456. [PMID 30654643](https://pubmed.ncbi.nlm.nih.gov/30654643/) — modest weight benefit with no psychopathology benefit; the reason `EMAX_OLZ_APP` >> `EMAX_OLZ_COG`.
61. Andries A, Frystyk J, Flyvbjerg A, et al. Dronabinol in severe, enduring anorexia nervosa: a randomized controlled trial. *Int J Eat Disord*. 2014;47:18-23. [PMID 24105610](https://pubmed.ncbi.nlm.nih.gov/24105610/)
62. Blanchet C, Guillaume S, Bat-Pitault F, et al. Medication in AN: a multidisciplinary overview of meta-analyses and systematic reviews. *J Clin Med*. 2019;8:278. [PMID 30823566](https://pubmed.ncbi.nlm.nih.gov/30823566/)

## 12. Psychological Treatment & Clinical Guidelines (map cluster 20)

63. Lock J, Le Grange D, Agras WS, et al. Randomized clinical trial comparing family-based treatment with adolescent-focused individual therapy for adolescents with anorexia nervosa. *Arch Gen Psychiatry*. 2010;67:1025-32. [PMID 20921118](https://pubmed.ncbi.nlm.nih.gov/20921118/) — the FBT advantage encoded as `FBT_BONUS`.
64. Fairburn CG, Cooper Z, Doll HA, et al. Enhanced cognitive behaviour therapy for adults with anorexia nervosa: a UK-Italy study. *Behav Res Ther*. 2013;51:R2-8. [PMID 23084515](https://pubmed.ncbi.nlm.nih.gov/23084515/)
65. Schmidt U, Magill N, Renwick B, et al. The Maudsley Outpatient Study of Treatments for Anorexia Nervosa and Related Conditions (MOSAIC): MANTRA versus specialist supportive clinical management. *J Consult Clin Psychol*. 2015;83:796-807. [PMID 25984803](https://pubmed.ncbi.nlm.nih.gov/25984803/)
66. Le Grange D, Lock J, Agras WS, et al. Randomized clinical trial of family-based treatment and cognitive-behavioral therapy for adolescent bulimia nervosa. *J Am Acad Child Adolesc Psychiatry*. 2015;54:886-94. [PMID 26506579](https://pubmed.ncbi.nlm.nih.gov/26506579/)
67. Hay P, Chinn D, Forbes D, et al. Royal Australian and New Zealand College of Psychiatrists clinical practice guidelines for the treatment of eating disorders. *Aust N Z J Psychiatry*. 2014;48:977-1008. [PMID 25351912](https://pubmed.ncbi.nlm.nih.gov/25351912/)
68. Hilbert A, Hoek HW, Schmidt R. Evidence-based clinical guidelines for eating disorders: international comparison. *Curr Opin Psychiatry*. 2017;30:423-437. [PMID 28777107](https://pubmed.ncbi.nlm.nih.gov/28777107/)
69. Society for Adolescent Health and Medicine. Medical management of restrictive eating disorders in adolescents and young adults. *J Adolesc Health*. 2022;71:648-654. [PMID 36058805](https://pubmed.ncbi.nlm.nih.gov/36058805/) — medical-instability criteria behind the `Medically_unstable` flag.

## 13. Systems / QSP Methodology

70. Nijhout HF, Best JA, Reed MC. Systems biology of robustness and homeostatic mechanisms. *Wiley Interdiscip Rev Syst Biol Med*. 2019;11:e1440. [PMID 30371009](https://pubmed.ncbi.nlm.nih.gov/30371009/)
71. Baron KT, Gastonguay MR. *mrgsolve: Simulate from ODE-Based Population PK/PD and Systems Pharmacology Models*. R package. [https://mrgsolve.org](https://mrgsolve.org) — the simulation engine used here.

---

## How the evidence maps onto the model

| Model element | Primary sources |
|---|---|
| FM/FFM energy partition (`FORBES_C`, `RHO_*`) | 17, 18 |
| Adaptive thermogenesis (`ADAPT_MAX`) & refeeding hypermetabolism (`K_HYPERMET`) | 19, 20 |
| Starvation-perpetuated cognition (`K_STARVE_DRV`) | 16, 12, 13 |
| Hypoleptinemia (`K_LEP_FM`, `LEP_ACUTE`) | 25, 23 |
| Low-T3 syndrome (`T3_SUPP`) → bradycardia (`HR_BRADY_MAX`) | 27, 43 |
| Hypercortisolemia (`K_CORT_DEF`) → lean & bone catabolism | 28, 29 |
| GH resistance / low IGF-1 (`IGF_SUPP`) | 30, 37 |
| Leptin-gated GnRH pulse generator (`LEP50_GNRH` = 1.85 ng/mL) | 32, 33 |
| Low-formation / high-resorption bone (`P1NP`, `CTX`, `BMD`) | 34, 35, 40 |
| Transdermal E2 vs oral contraceptive (`E2_PATCH_GAIN`, `OCP_IGF_SUPP`) | 36, 37 |
| Teriparatide (`EMAX_TPTD`) / risedronate (`EMAX_BIS`) | 38, 39 |
| Refeeding insulin surge, PO4/K/Mg shift, thiamine | 50, 51, 54, 55 |
| Higher-calorie vs conservative refeeding (scenarios 2-5) | 52, 53, 56 |
| SSRI nutritional gate (`GATE_MBMI50`, `GATE_HILL`) | 57, 58 |
| Olanzapine appetite >> cognition (`EMAX_OLZ_APP` vs `EMAX_OLZ_COG`) | 59, 60 |
| FBT advantage in adolescents (`FBT_BONUS`) | 63 |
| Medical-instability criteria (`Medically_unstable`) | 69, 41, 42 |

---

## Disclaimer

This bibliography supports an **educational and research** QSP model. The model
is not validated for clinical decision-making, and nothing in it constitutes a
refeeding protocol. Refeeding a severely malnourished patient is an inpatient
medical procedure requiring specialist supervision and electrolyte monitoring.
