# Lipodystrophy syndromes — annotated reference list

**Every PubMed identifier below was resolved against the NCBI E-utilities API,
and the title, first author, journal and year printed here are what PubMed
returns for that identifier** — not what this file was drafted with. That matters:
of eight identifiers written from memory during drafting, four turned out to
point at unrelated papers (a herpesvirus dihydrofolate reductase, a linkage
analysis program, and two others). They were dropped rather than corrected by
guesswork, and no entry appears here whose metadata was not fetched.

The annotation after each entry says what the reference is used FOR in
`lipo_mrgsolve_model.R` and `lipo_qsp_model.dot`: a parameter, a calibration
target, or a structural claim. Two annotations record DISAGREEMENT between the
model and the literature; those are the entries to read first.

## 1. Classification, epidemiology and diagnosis

1. Brown RJ. *The Diagnosis and Management of Lipodystrophy Syndromes: A Multi-Society Practice Guideline*. J Clin Endocrinol Metab 2016. [PMID 27710244](https://pubmed.ncbi.nlm.nih.gov/27710244/) — consensus guideline
2. Feingold KR. *Lipodystrophy Syndromes: Presentation and Treatment*. 2000. [PMID 29989768](https://pubmed.ncbi.nlm.nih.gov/29989768/) — clinical overview
3. Araújo-Vilar D. *Lipodystrophies: A diagnostic and therapeutic challenge*. Endocrinol Diabetes Nutr (Engl Ed) 2026. [PMID 42120113](https://pubmed.ncbi.nlm.nih.gov/42120113/) — current classification, diagnosis and management review
4. Akinci B. *Natural History of Congenital Generalized Lipodystrophy: A Nationwide Study From Turkey*. J Clin Endocrinol Metab 2016. [PMID 27144933](https://pubmed.ncbi.nlm.nih.gov/27144933/) — natural history
5. Santos MG. *[Generalized congenital lipodystrophy: correlation with leptin and other biochemical parameters]*. Acta Cir Bras 2005. [PMID 17768808](https://pubmed.ncbi.nlm.nih.gov/17768808/) — leptin in CGL
6. Cook K. *Effects of Metreleptin on Patient Outcomes and Quality of Life in Generalized and Partial Lipodystrophy*. J Endocr Soc 2021. [PMID 33817539](https://pubmed.ncbi.nlm.nih.gov/33817539/) — quality-of-life outcomes on metreleptin

## 2. Genetics of generalised lipodystrophy

7. Agarwal AK. *AGPAT2 is mutated in congenital generalized lipodystrophy linked to chromosome 9q34*. Nat Genet 2002. [PMID 11967537](https://pubmed.ncbi.nlm.nih.gov/11967537/) — CGL1
8. Magré J. *Identification of the gene altered in Berardinelli-Seip congenital lipodystrophy on chromosome 11q13*. Nat Genet 2001. [PMID 11479539](https://pubmed.ncbi.nlm.nih.gov/11479539/) — CGL2 / seipin
9. Agarwal AK. *Phenotypic and genetic heterogeneity in congenital generalized lipodystrophy*. J Clin Endocrinol Metab 2003. [PMID 14557463](https://pubmed.ncbi.nlm.nih.gov/14557463/) — CGL1 vs CGL2
10. Savage DB. *Human metabolic syndrome resulting from dominant-negative mutations in the nuclear receptor peroxisome proliferator-activated receptor-gamma*. Diabetes 2003. [PMID 12663460](https://pubmed.ncbi.nlm.nih.gov/12663460/) — PPARG
11. Mancioppi V. *A new mutation in the CAVIN1/PTRF gene in two siblings with congenital generalized lipodystrophy type 4: case reports and review of the literature*. Front Endocrinol (Lausanne) 2023. [PMID 37501786](https://pubmed.ncbi.nlm.nih.gov/37501786/) — CGL4 / CAVIN1-PTRF, with myopathy and arrhythmia
12. Wang H. *Seipin is required for converting nascent to mature lipid droplets*. Elife 2016. [PMID 27564575](https://pubmed.ncbi.nlm.nih.gov/27564575/) — seipin function
13. Combot Y. *Seipin localizes at endoplasmic-reticulum-mitochondria contact sites to control mitochondrial calcium import and metabolism in adipocytes*. Cell Rep 2022. [PMID 35021082](https://pubmed.ncbi.nlm.nih.gov/35021082/) — seipin biology
14. Softic S. *Lipodystrophy Due to Adipose Tissue-Specific Insulin Receptor Knockout Results in Progressive NAFLD*. Diabetes 2016. [PMID 27207510](https://pubmed.ncbi.nlm.nih.gov/27207510/) — mouse fat-specific model
15. Rajan R. *A series of genetically confirmed congenital lipodystrophy and diabetes in adult southern Indian patients*. Sci Rep 2024. [PMID 39550450](https://pubmed.ncbi.nlm.nih.gov/39550450/) — genetically confirmed congenital lipodystrophy series

## 3. Genetics of partial lipodystrophy

16. Cao H. *Nuclear lamin A/C R482Q mutation in canadian kindreds with Dunnigan-type familial partial lipodystrophy*. Hum Mol Genet 2000. [PMID 10587585](https://pubmed.ncbi.nlm.nih.gov/10587585/) — FPLD2
17. O'Rahilly S. *Human genetics illuminates the paths to metabolic disease*. Nature 2009. [PMID 19924209](https://pubmed.ncbi.nlm.nih.gov/19924209/) — context
18. Gandotra S. *Perilipin deficiency and autosomal dominant partial lipodystrophy*. N Engl J Med 2011. [PMID 21345103](https://pubmed.ncbi.nlm.nih.gov/21345103/) — FPLD4 / PLIN1
19. Visser ME. *Characterisation of non-obese diabetic patients with marked insulin resistance identifies a novel familial partial lipodystrophy-associated PPARγ mutation (Y151C)*. Diabetologia 2011. [PMID 21479595](https://pubmed.ncbi.nlm.nih.gov/21479595/) — severe insulin resistance without obesity -- the phenotype the capacity argument predicts
20. Fernández-Pombo A. *A cohort analysis of familial partial lipodystrophy from two Mediterranean countries*. Diabetes Obes Metab 2024. [PMID 39171574](https://pubmed.ncbi.nlm.nih.gov/39171574/) — FPLD2 genotype-phenotype
21. Jéru I. *Diagnostic Challenge in PLIN1-Associated Familial Partial Lipodystrophy*. J Clin Endocrinol Metab 2019. [PMID 31504636](https://pubmed.ncbi.nlm.nih.gov/31504636/) — PLIN1-associated FPLD4, diagnostic difficulty
22. George S. *A family with severe insulin resistance and diabetes due to a mutation in AKT2*. Science 2004. [PMID 15166380](https://pubmed.ncbi.nlm.nih.gov/15166380/) — AKT2
23. Farhan SM. *A novel LIPE nonsense mutation found using exome sequencing in siblings with late-onset familial partial lipodystrophy*. Can J Cardiol 2014. [PMID 25475467](https://pubmed.ncbi.nlm.nih.gov/25475467/) — LIPE
24. Shackleton S. *LMNA, encoding lamin A/C, is mutated in partial lipodystrophy*. Nat Genet 2000. [PMID 10655060](https://pubmed.ncbi.nlm.nih.gov/10655060/) — FPLD2: LMNA is the partial-lipodystrophy gene

## 4. Adipose tissue expandability, lipid droplets and the buffer hypothesis

25. Carobbio S. *Adipose Tissue Function and Expandability as Determinants of Lipotoxicity and the Metabolic Syndrome*. Adv Exp Med Biol 2017. [PMID 28585199](https://pubmed.ncbi.nlm.nih.gov/28585199/) — the expandability hypothesis
26. Frayn KN. *Adipose tissue as a buffer for daily lipid flux*. Diabetologia 2002. [PMID 12242452](https://pubmed.ncbi.nlm.nih.gov/12242452/) — the phrase this model is built on: adipose tissue as a buffer for daily lipid flux
27. Chait A. *Adipose Tissue Distribution, Inflammation and Its Metabolic Consequences, Including Diabetes and Cardiovascular Disease*. Front Cardiovasc Med 2020. [PMID 32158768](https://pubmed.ncbi.nlm.nih.gov/32158768/) — depot biology
28. Lefebvre AM. *Depot-specific differences in adipose tissue gene expression in lean and obese subjects*. Diabetes 1998. [PMID 9421381](https://pubmed.ncbi.nlm.nih.gov/9421381/) — depot heterogeneity
29. Nielsen S. *Splanchnic lipolysis in human obesity*. J Clin Invest 2004. [PMID 15173884](https://pubmed.ncbi.nlm.nih.gov/15173884/) — portal delivery
30. Spalding KL. *Dynamics of fat cell turnover in humans*. Nature 2008. [PMID 18454136](https://pubmed.ncbi.nlm.nih.gov/18454136/) — adipocyte turnover
31. Lee JE. *Transcriptional and epigenetic regulation of PPARγ expression during adipogenesis*. Cell Biosci 2014. [PMID 24904744](https://pubmed.ncbi.nlm.nih.gov/24904744/) — PPARg adipogenesis
32. Bergman RN. *Free fatty acids and pathogenesis of type 2 diabetes mellitus*. Trends Endocrinol Metab 2000. [PMID 11042464](https://pubmed.ncbi.nlm.nih.gov/11042464/) — free fatty acids and the portal hypothesis
33. Virtue S. *Adipose tissue expandability, lipotoxicity and the Metabolic Syndrome--an allostatic perspective*. Biochim Biophys Acta 2010. [PMID 20056169](https://pubmed.ncbi.nlm.nih.gov/20056169/) — adipose expandability and lipotoxicity -- the buffer thesis

## 5. Leptin biology and the hypothalamic transducer

34. Zhang Y. *Positional cloning of the mouse obese gene and its human homologue*. Nature 1994. [PMID 7984236](https://pubmed.ncbi.nlm.nih.gov/7984236/) — leptin discovery
35. Tartaglia LA. *Identification and expression cloning of a leptin receptor, OB-R*. Cell 1995. [PMID 8548812](https://pubmed.ncbi.nlm.nih.gov/8548812/) — LepR
36. Considine RV. *Serum immunoreactive-leptin concentrations in normal-weight and obese humans*. N Engl J Med 1996. [PMID 8532024](https://pubmed.ncbi.nlm.nih.gov/8532024/) — leptin vs adiposity
37. Friedman JM. *Leptin and the regulation of body weight in mammals*. Nature 1998. [PMID 9796811](https://pubmed.ncbi.nlm.nih.gov/9796811/) — physiology review
38. Cowley MA. *Leptin activates anorexigenic POMC neurons through a neural network in the arcuate nucleus*. Nature 2001. [PMID 11373681](https://pubmed.ncbi.nlm.nih.gov/11373681/) — POMC/AgRP
39. Howard JK. *Enhanced leptin sensitivity and attenuation of diet-induced obesity in mice with haploinsufficiency of Socs3*. Nat Med 2004. [PMID 15220914](https://pubmed.ncbi.nlm.nih.gov/15220914/) — SOCS3 brake
40. Zabolotny JM. *PTP1B regulates leptin signal transduction in vivo*. Dev Cell 2002. [PMID 11970898](https://pubmed.ncbi.nlm.nih.gov/11970898/) — PTP1B brake
41. Gruzdeva O. *Leptin resistance: underlying mechanisms and diagnosis*. Diabetes Metab Syndr Obes 2019. [PMID 30774404](https://pubmed.ncbi.nlm.nih.gov/30774404/) — why leptin fails in obesity
42. Heymsfield SB. *Recombinant leptin for weight loss in obese and lean adults: a randomized, controlled, dose-escalation trial*. JAMA 1999. [PMID 10546697](https://pubmed.ncbi.nlm.nih.gov/10546697/) — leptin fails in common obesity
43. Farooqi IS. *Beneficial effects of leptin on obesity, T cell hyporesponsiveness, and neuroendocrine/metabolic dysfunction of human congenital leptin deficiency*. J Clin Invest 2002. [PMID 12393845](https://pubmed.ncbi.nlm.nih.gov/12393845/) — leptin deficiency treated
44. Davis JF. *Leptin regulates energy balance and motivation through action at distinct neural circuits*. Biol Psychiatry 2011. [PMID 21035790](https://pubmed.ncbi.nlm.nih.gov/21035790/) — reward
45. Cheung CC. *Leptin is a metabolic gate for the onset of puberty in the female rat*. Endocrinology 1997. [PMID 9003028](https://pubmed.ncbi.nlm.nih.gov/9003028/) — reproductive permissiveness
46. Lawler K. *Leptin-Mediated Changes in the Human Metabolome*. J Clin Endocrinol Metab 2020. [PMID 32392278](https://pubmed.ncbi.nlm.nih.gov/32392278/) — metabolomic signature of leptin replacement
47. Ahima RS. *Role of leptin in the neuroendocrine response to fasting*. Nature 1996. [PMID 8717038](https://pubmed.ncbi.nlm.nih.gov/8717038/) — leptin as a starvation signal

## 6. Metreleptin: pharmacology and clinical trials

48. Vahidi Rad M. *Leptin Therapy Improves Metabolic Dysfunction in Immune Checkpoint Inhibitor-induced Lipodystrophy*. JCEM Case Rep 2026. [PMID 41700130](https://pubmed.ncbi.nlm.nih.gov/41700130/) — acquired lipodystrophy after immune-checkpoint inhibition, treated with leptin
49. Petersen KF. *Leptin reverses insulin resistance and hepatic steatosis in patients with severe lipodystrophy*. J Clin Invest 2002. [PMID 12021250](https://pubmed.ncbi.nlm.nih.gov/12021250/) — KEY calibration target: hepatic triglyceride falls ~86% on leptin
50. Javor ED. *Long-term efficacy of leptin replacement in patients with generalized lipodystrophy*. Diabetes 2005. [PMID 15983199](https://pubmed.ncbi.nlm.nih.gov/15983199/) — Javor 2005
51. Brush M. *Effects of Metreleptin in Patients With Generalized Lipodystrophy Before vs After the Onset of Severe Metabolic Disease*. J Clin Endocrinol Metab 2025. [PMID 38757950](https://pubmed.ncbi.nlm.nih.gov/38757950/) — metreleptin outcomes, treated vs untreated periods
52. Diker-Cohen T. *Partial and generalized lipodystrophy: comparison of baseline characteristics and response to metreleptin*. J Clin Endocrinol Metab 2015. [PMID 25734254](https://pubmed.ncbi.nlm.nih.gov/25734254/) — KEY calibration target: the generalised/partial dissociation in response
53. Brown RJ. *Metreleptin-mediated improvements in insulin sensitivity are independent of food intake in humans with lipodystrophy*. J Clin Invest 2018. [PMID 29723161](https://pubmed.ncbi.nlm.nih.gov/29723161/) — KEY and the reference this model's inference 2 DISAGREES with: human data suggest a larger intake-independent component than the model reproduces
54. Ajluni N. *Efficacy and Safety of Metreleptin in Patients with Partial Lipodystrophy: Lessons from an Expanded Access Program*. J Diabetes Metab 2016. [PMID 27642538](https://pubmed.ncbi.nlm.nih.gov/27642538/) — PL expanded access
55. Akinci G. *Metreleptin Treatment in Patients with Non-HIV Associated Lipodystrophy*. Recent Pat Endocr Metab Immune Drug Discov 2015. [PMID 26556498](https://pubmed.ncbi.nlm.nih.gov/26556498/)
56. Akinci B. *The complicated clinical course in a case of atypical lipodystrophy after development of neutralizing antibody to metreleptin: treatment with setmelanotide*. Endocrinol Diabetes Metab Case Rep 2020. [PMID 32213649](https://pubmed.ncbi.nlm.nih.gov/32213649/)
57. Baykal AP. *Leptin decreases de novo lipogenesis in patients with lipodystrophy*. JCI Insight 2020. [PMID 32573497](https://pubmed.ncbi.nlm.nih.gov/32573497/) — KEY: leptin lowers de novo lipogenesis in patients -- the KLDNL term
58. Oral EA. *Effect of leptin replacement on pituitary hormone regulation in patients with severe lipodystrophy*. J Clin Endocrinol Metab 2002. [PMID 12107209](https://pubmed.ncbi.nlm.nih.gov/12107209/) — endocrine effects
59. Magkos F. *Leptin replacement improves postprandial glycemia and insulin sensitivity in human immunodeficiency virus-infected lipoatrophic men treated with pioglitazone: a pilot study*. Metabolism 2011. [PMID 21081243](https://pubmed.ncbi.nlm.nih.gov/21081243/)
60. Muniyappa R. *Metreleptin therapy lowers plasma angiopoietin-like protein 3 in patients with generalized lipodystrophy*. J Clin Lipidol 2017. [PMID 28502512](https://pubmed.ncbi.nlm.nih.gov/28502512/)
61. Chan JL. *Immunogenicity associated with metreleptin treatment in patients with obesity or lipodystrophy*. Clin Endocrinol (Oxf) 2016. [PMID 26589105](https://pubmed.ncbi.nlm.nih.gov/26589105/) — immunogenicity: binding and neutralising anti-metreleptin antibodies (scenario 11)
62. Perakakis N. *Evidence from clinical studies of leptin: current and future clinical applications in humans*. Metabolism 2024. [PMID 39490439](https://pubmed.ncbi.nlm.nih.gov/39490439/) — current and future clinical applications of leptin
63. Oral EA. *Leptin-replacement therapy for lipodystrophy*. N Engl J Med 2002. [PMID 11856796](https://pubmed.ncbi.nlm.nih.gov/11856796/) — KEY: the pivotal observation that leptin replacement transforms generalised lipodystrophy
64. Chan JL. *Clinical effects of long-term metreleptin treatment in patients with lipodystrophy*. Endocr Pract 2011. [PMID 22068254](https://pubmed.ncbi.nlm.nih.gov/22068254/) — long-term metreleptin cohort, HbA1c and triglyceride targets

## 7. Ectopic lipid, DAG-PKC and selective hepatic insulin resistance

65. Samuel VT. *Inhibition of protein kinase Cepsilon prevents hepatic insulin resistance in nonalcoholic fatty liver disease*. J Clin Invest 2007. [PMID 17318260](https://pubmed.ncbi.nlm.nih.gov/17318260/)
66. Brown MS. *Selective versus total insulin resistance: a pathogenic paradox*. Cell Metab 2008. [PMID 18249166](https://pubmed.ncbi.nlm.nih.gov/18249166/) — Brown & Goldstein
67. Krssak M. *Intramyocellular lipid concentrations are correlated with insulin sensitivity in humans: a 1H NMR spectroscopy study*. Diabetologia 1999. [PMID 10027589](https://pubmed.ncbi.nlm.nih.gov/10027589/)
68. Donnelly KL. *Sources of fatty acids stored in liver and secreted via lipoproteins in patients with nonalcoholic fatty liver disease*. J Clin Invest 2005. [PMID 15864352](https://pubmed.ncbi.nlm.nih.gov/15864352/) — Donnelly 2005, DNL contribution
69. Adiels M. *Overproduction of large VLDL particles is driven by increased liver fat content in man*. Diabetologia 2006. [PMID 16463046](https://pubmed.ncbi.nlm.nih.gov/16463046/) — Adiels, substrate-driven VLDL
70. Bikman BT. *Ceramides as modulators of cellular and whole-body metabolism*. J Clin Invest 2011. [PMID 22045572](https://pubmed.ncbi.nlm.nih.gov/22045572/)
71. Petersen KF. *Reversal of muscle insulin resistance by weight reduction in young, lean, insulin-resistant offspring of parents with type 2 diabetes*. Proc Natl Acad Sci U S A 2012. [PMID 22547801](https://pubmed.ncbi.nlm.nih.gov/22547801/) — reversibility of IMCL-driven IR
72. Samuel VT. *Mechanisms for insulin resistance: common threads and missing links*. Cell 2012. [PMID 22385956](https://pubmed.ncbi.nlm.nih.gov/22385956/) — the DAG-PKC account of lipid-induced insulin resistance

## 8. Liver disease in lipodystrophy

73. Safar Zadeh E. *The liver diseases of lipodystrophy: the long-term effect of leptin treatment*. J Hepatol 2013. [PMID 23439261](https://pubmed.ncbi.nlm.nih.gov/23439261/) — Safar Zadeh 2013
74. Singh S. *Fibrosis progression in nonalcoholic fatty liver vs nonalcoholic steatohepatitis: a systematic review and meta-analysis of paired-biopsy studies*. Clin Gastroenterol Hepatol 2015. [PMID 24768810](https://pubmed.ncbi.nlm.nih.gov/24768810/) — fibrosis rates
75. Friedman SL. *Hepatic stellate cells: protean, multifunctional, and enigmatic cells of the liver*. Physiol Rev 2008. [PMID 18195085](https://pubmed.ncbi.nlm.nih.gov/18195085/) — fibrogenesis
76. Ajluni N. *Spectrum of disease associated with partial lipodystrophy: lessons from a trial cohort*. Clin Endocrinol (Oxf) 2017. [PMID 28199729](https://pubmed.ncbi.nlm.nih.gov/28199729/) — spectrum of partial lipodystrophy in a translational cohort
77. Brown RJ. *Effects of Metreleptin in Pediatric Patients With Lipodystrophy*. J Clin Endocrinol Metab 2017. [PMID 28324110](https://pubmed.ncbi.nlm.nih.gov/28324110/) — paediatric metreleptin
78. Qu Y. *Diagnostic accuracy of hepatic proton density fat fraction measured by magnetic resonance imaging for the evaluation of liver steatosis with histology as reference standard: a meta-analysis*. Eur Radiol 2019. [PMID 30877459](https://pubmed.ncbi.nlm.nih.gov/30877459/) — MRI-PDFF

## 9. Triglyceride metabolism, APOC3, ANGPTL3 and pancreatitis

79. Young SG. *GPIHBP1 and Lipoprotein Lipase, Partners in Plasma Triglyceride Metabolism*. Cell Metab 2019. [PMID 31269429](https://pubmed.ncbi.nlm.nih.gov/31269429/)
80. Jong MC. *Apolipoprotein C-III deficiency accelerates triglyceride hydrolysis by lipoprotein lipase in wild-type and apoE knockout mice*. J Lipid Res 2001. [PMID 11590213](https://pubmed.ncbi.nlm.nih.gov/11590213/) — apoC-III restrains LPL-mediated hydrolysis (F_APOC3 -> F_LPL)
81. Musunuru K. *Exome sequencing, ANGPTL3 mutations, and familial combined hypolipidemia*. N Engl J Med 2010. [PMID 20942659](https://pubmed.ncbi.nlm.nih.gov/20942659/)
82. Volpi L. *Evinacumab Treatment of Homozygous Familial Hypercholesterolemia Patients in a Real-World Setting: An Overview After the First Year of Therapy*. AACE Endocrinol Diabetes 2026. [PMID 42491512](https://pubmed.ncbi.nlm.nih.gov/42491512/) — evinacumab in routine practice
83. Rosenson RS. *Longer-Term Efficacy and Safety of Evinacumab in Patients With Refractory Hypercholesterolemia*. JAMA Cardiol 2023. [PMID 37703006](https://pubmed.ncbi.nlm.nih.gov/37703006/)
84. Kota SK. *Hypertriglyceridemia-induced recurrent acute pancreatitis: A case-based review*. Indian J Endocrinol Metab 2012. [PMID 22276267](https://pubmed.ncbi.nlm.nih.gov/22276267/)
85. Hansen SEJ. *Genetic Variants Associated With Increased Plasma Levels of Triglycerides, via Effects on the Lipoprotein Lipase Pathway, Increase Risk of Acute Pancreatitis*. Clin Gastroenterol Hepatol 2021. [PMID 32801009](https://pubmed.ncbi.nlm.nih.gov/32801009/) — triglyceride-raising variants and pancreatitis risk
86. Mozaffarian D. *Omega-3 fatty acids and cardiovascular disease: effects on risk factors, molecular pathways, and clinical events*. J Am Coll Cardiol 2011. [PMID 22051327](https://pubmed.ncbi.nlm.nih.gov/22051327/)
87. Liu C. *Lipoprotein lipase transporter GPIHBP1 and triglyceride-rich lipoprotein metabolism*. Clin Chim Acta 2018. [PMID 30218660](https://pubmed.ncbi.nlm.nih.gov/30218660/) — GPIHBP1 and triglyceride-rich lipoprotein metabolism
88. Oral EA. *Assessment of efficacy and safety of volanesorsen for treatment of metabolic complications in patients with familial partial lipodystrophy: Results of the BROADEN study: Volanesorsen in FPLD; The BROADEN Study*. J Clin Lipidol 2022. [PMID 36402670](https://pubmed.ncbi.nlm.nih.gov/36402670/) — KEY calibration target: volanesorsen in FPLD -- triglyceride collapses, hepatic fat barely moves (inference 3)
89. Olesen SS. *Hypertriglyceridemia is often under recognized as an aetiologic risk factor for acute pancreatitis: A population-based cohort study*. Pancreatology 2021. [PMID 33608229](https://pubmed.ncbi.nlm.nih.gov/33608229/) — TG and pancreatitis risk
90. Lefebvre AM. *Regulation of lipoprotein metabolism by thiazolidinediones occurs through a distinct but complementary mechanism relative to fibrates*. Arterioscler Thromb Vasc Biol 1997. [PMID 9327774](https://pubmed.ncbi.nlm.nih.gov/9327774/) — thiazolidinedione effects on lipoprotein metabolism

## 10. Thiazolidinediones and capacity-directed therapy

91. Hwang YC. *Effects of rosiglitazone on body fat distribution and insulin sensitivity in Korean type 2 diabetes mellitus patients*. Metabolism 2008. [PMID 18328348](https://pubmed.ncbi.nlm.nih.gov/18328348/)
92. Simha V. *Prolonged thiazolidinedione therapy does not reverse fat loss in patients with familial partial lipodystrophy, Dunnigan variety*. Diabetes Obes Metab 2008. [PMID 19040647](https://pubmed.ncbi.nlm.nih.gov/19040647/) — KEY for inference 4: prolonged thiazolidinedione therapy does NOT reverse fat loss where no expandable depot remains
93. Sharma AM. *Review: Peroxisome proliferator-activated receptor gamma and adipose tissue--understanding obesity-related changes in regulation of lipid and glucose metabolism*. J Clin Endocrinol Metab 2007. [PMID 17148564](https://pubmed.ncbi.nlm.nih.gov/17148564/)
94. Yki-Järvinen H. *Thiazolidinediones and the liver in humans*. Curr Opin Lipidol 2009. [PMID 19779336](https://pubmed.ncbi.nlm.nih.gov/19779336/) — TZD review
95. Palavicini JP. *The Insulin-Sensitizer Pioglitazone Remodels Adipose Tissue Phospholipids in Humans*. Front Physiol 2021. [PMID 34925073](https://pubmed.ncbi.nlm.nih.gov/34925073/) — pioglitazone remodels adipose tissue -- the capacity term
96. Luedtke A. *Thiazolidinedione response in familial lipodystrophy patients with LMNA mutations: a case series*. Horm Metab Res 2012. [PMID 22274718](https://pubmed.ncbi.nlm.nih.gov/22274718/) — TZD in FPLD
97. Slama L. *Effect of pioglitazone on HIV-1-related lipodystrophy: a randomized double-blind placebo-controlled trial (ANRS 113)*. Antivir Ther 2008. [PMID 18389900](https://pubmed.ncbi.nlm.nih.gov/18389900/) — TZD in acquired PL
98. Christensen ML. *Single- and multiple-dose pharmacokinetics of pioglitazone in adolescents with type 2 diabetes*. J Clin Pharmacol 2005. [PMID 16172178](https://pubmed.ncbi.nlm.nih.gov/16172178/) — pioglitazone pharmacokinetics (KAPIO, VPIO, CLPIO)

## 11. Glycaemia, incretins and other adjunctive therapy

99. Lamothe S. *Safety and effectiveness in an uncontrolled setting of glucagon-like-peptide-1 receptor agonists in patients with familial partial lipodystrophy: Real-life experience from a national reference network*. Diabetes Obes Metab 2025. [PMID 39829337](https://pubmed.ncbi.nlm.nih.gov/39829337/)
100. Li Y. *Liraglutide use in pediatric type 2 familial partial lipodystrophy caused by LMNA mutation: a case report*. BMC Pediatr 2025. [PMID 40619352](https://pubmed.ncbi.nlm.nih.gov/40619352/) — incretin in lipodystrophy
101. Zamora JM. *Case Series of U-500 Insulin Use in Adults With Type 2 Diabetes and Severe Insulin Resistance*. Can J Diabetes 2021. [PMID 32847768](https://pubmed.ncbi.nlm.nih.gov/32847768/) — U-500 insulin
102. Simha V. *Inherited lipodystrophies and hypertriglyceridemia*. Curr Opin Lipidol 2009. [PMID 19494770](https://pubmed.ncbi.nlm.nih.gov/19494770/) — inherited lipodystrophy and hypertriglyceridaemia; the dietary-fat argument (scenario 18)

## 12. Renal, cardiac, reproductive and other complications

103. Hoff FW. *Genotype-phenotype heterogeneity among patients with lipodystrophy harboring rare POLD1 variants*. J Clin Endocrinol Metab 2026. [PMID 41742372](https://pubmed.ncbi.nlm.nih.gov/41742372/)
104. Liberato CBR. *Early Left Ventricular Systolic Dysfunction Detected by Two-Dimensional Speckle-Tracking Echocardiography in Young Patients with Congenital Generalized Lipodystrophy*. Diabetes Metab Syndr Obes 2020. [PMID 32021357](https://pubmed.ncbi.nlm.nih.gov/32021357/)
105. Hoff FW. *A Novel Subtype of Acquired Generalized Lipodystrophy Associated With Subcutaneous Panniculitis-Like T-cell Lymphoma*. JCEM Case Rep 2024. [PMID 38681964](https://pubmed.ncbi.nlm.nih.gov/38681964/)
106. Brănişteanu DD. *Barraquer-Simons syndrome. Report of a case and review of the literature*. Rev Med Chir Soc Med Nat Iasi 2000. [PMID 12089983](https://pubmed.ncbi.nlm.nih.gov/12089983/)
107. Power DA. *Familial incidence of C3 nephritic factor, partial lipodystrophy and membranoproliferative glomerulonephritis*. Q J Med 1990. [PMID 2385743](https://pubmed.ncbi.nlm.nih.gov/2385743/)
108. Mosbah H. *Health-related Quality of Life, Social, and Psychological Well-Being of 109 Adult Patients With Genetic Lipodystrophy*. J Clin Endocrinol Metab 2025. [PMID 39657019](https://pubmed.ncbi.nlm.nih.gov/39657019/)
109. Morguetti MJ. *Podocytopathies associated with familial partial lipodystrophy due to LMNA variants: report of two cases*. Arch Endocrinol Metab 2024. [PMID 38739524](https://pubmed.ncbi.nlm.nih.gov/38739524/) — renal phenotype
110. Nguyen ML. *Leptin Attenuates Cardiac Hypertrophy in Patients With Generalized Lipodystrophy*. J Clin Endocrinol Metab 2021. [PMID 34223895](https://pubmed.ncbi.nlm.nih.gov/34223895/) — cardiac phenotype
111. Corvillo F. *An overview of lipodystrophy and the role of the complement system*. Mol Immunol 2019. [PMID 31177059](https://pubmed.ncbi.nlm.nih.gov/31177059/) — APL
112. Brown RJ. *Lymphoma in acquired generalized lipodystrophy*. Leuk Lymphoma 2016. [PMID 25864863](https://pubmed.ncbi.nlm.nih.gov/25864863/) — lymphoma association

## 13. Modelling methods and quantitative pharmacology

113. Musante CJ. *Quantitative Systems Pharmacology: A Case for Disease Models*. Clin Pharmacol Ther 2017. [PMID 27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)
114. Hall KD. *Quantification of the effect of energy imbalance on bodyweight*. Lancet 2011. [PMID 21872751](https://pubmed.ncbi.nlm.nih.gov/21872751/) — Hall energy-balance model
115. Sayed A. *Translating the HbA1c assay into estimated average glucose values in children and adolescents with type 1 diabetes mellitus*. Acta Biomed 2018. [PMID 30049928](https://pubmed.ncbi.nlm.nih.gov/30049928/) — secondary HbA1c-glucose calibration
116. Visentin R. *Hepatic insulin sensitivity in healthy and prediabetic subjects: from a dual- to a single-tracer oral minimal model*. Am J Physiol Endocrinol Metab 2015. [PMID 25991649](https://pubmed.ncbi.nlm.nih.gov/25991649/) — dual-tracer hepatic insulin sensitivity, used for the selective-resistance term
117. Upton RN. *Basic concepts in population modeling, simulation, and model-based drug development: part 3-introduction to pharmacodynamic modeling methods*. CPT Pharmacometrics Syst Pharmacol 2014. [PMID 24384783](https://pubmed.ncbi.nlm.nih.gov/24384783/)
118. Bai JP. *Quantitative Systems Pharmacology for Rare Disease Drug Development*. J Pharm Sci 2023. [PMID 37422281](https://pubmed.ncbi.nlm.nih.gov/37422281/) — QSP
119. Hall KD. *Obesity Energetics: Body Weight Regulation and the Effects of Diet Composition*. Gastroenterology 2017. [PMID 28193517](https://pubmed.ncbi.nlm.nih.gov/28193517/) — energy-balance model
120. Hijmans BS. *A systems biology approach reveals the physiological origin of hepatic steatosis induced by liver X receptor activation*. FASEB J 2015. [PMID 25477282](https://pubmed.ncbi.nlm.nih.gov/25477282/) — hepatic lipid model
121. Nathan DM. *Translating the A1C assay into estimated average glucose values*. Diabetes Care 2008. [PMID 18540046](https://pubmed.ncbi.nlm.nih.gov/18540046/) — the HbA1c-mean-glucose map used in dxdt_A1C

---

**121 references, all with PubMed identifiers verified against E-utilities.**

Retrieval note: identifiers were resolved by title search, then re-fetched by
identifier so that the printed metadata comes from PubMed rather than from the
query. Papers dated 2026 are recent additions to the index; the sections they
appear in also cite the original work where one exists.
