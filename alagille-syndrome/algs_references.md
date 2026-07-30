# Alagille Syndrome (ALGS) — QSP model references

154 references, every PMID resolved through the NCBI E-utilities API on
2026-07-30 rather than written from memory. Twelve further queries returned
no confident match and were dropped rather than guessed at.

## The five publications this model is calibrated to

Everything numeric in `algs_mrgsolve_model.R` traces back to one of these.
They are listed first because the model has no other quantitative input.

| # | What it supplies | Numbers used |
|---|---|---|
| 1 | **ICONIC** — maralixibat 380 µg/kg/day, n=31, randomised drug-withdrawal phase 2b in ALGS | withdrawal-period serum bile acid difference −117 µmol/L; week-48 sBA −96 µmol/L and ItchRO(Obs) −1.6 points; improvements maintained to week 204 in the 15 who continued |
| 2 | **ASSERT** — odevixibat 120 µg/kg/day, n=52 (35 drug / 17 placebo), phase 3 in ALGS | sBA 237→149 on drug vs 246→271 on placebo, LS-mean difference −113 µmol/L (95% CI −179 to −47), p=0.0012; PRUCISION scratching 2.8→1.1 vs 3.0→2.2, difference −0.9 (95% CI −1.4 to −0.3), p=0.0024 |
| 3 | **GALA** — 1 433 children, 67 centres, 29 countries, natural history | native liver survival 66.8% / 54.4% / 40.3% at 5 / 10 / 18 yr; adverse liver event 51.5% by 10 yr and 66.0% by 18 yr; relative to total bilirubin <5.0 mg/dL at 6–12 months, transplant hazard 4.8× at 5–<10 and 15.6× at ≥10; post-transplant patient survival 91.0% at 10 yr |
| 4 | **Maralixibat vs GALA event-free survival** — 84 treated to 6 yr vs 469 aligned GALA controls | HR 0.305 (95% CI 0.189–0.491), p<0.0001 |
| 5 | **MRGPRX4 / autotaxin pruritus mechanism** | the receptor identity behind the bile-acid arm of the itch model, and the reason TGR5 is carried as a null branch |

Two things follow from that list and are worth stating plainly. First, the
model has exactly **one** placebo-controlled bile-acid number and **one**
placebo-controlled itch number to calibrate against, which is why
`FAILURE 1` (non-identifiability of duct capacity, blockade and synthesis
reserve from a single ratio) is structural rather than fixable. Second, the
survival comparison in row 4 is not a randomised comparison, and the model's
Axis 4 is entirely about how much of it the drug's own bile-acid effect can
carry.

## How to read the rest

Sections follow the mechanistic map (`algs_qsp_model.dot`) cluster order.
A reference appearing here means the model uses the mechanism, not that the
model reproduces that paper's numbers — only the five above are calibration
targets.



## 1. Genetics — JAG1, NOTCH2 and the Notch pathway

1. Li L, et al. *Alagille syndrome is caused by mutations in human Jagged1, which encodes a ligand for Notch1.* Nat Genet 1997. [PMID 9207788](https://pubmed.ncbi.nlm.nih.gov/9207788/)
2. Oda T, et al. *Mutations in the human Jagged1 gene are responsible for Alagille syndrome.* Nat Genet 1997. [PMID 9207787](https://pubmed.ncbi.nlm.nih.gov/9207787/)
3. McDaniell R, et al. *NOTCH2 mutations cause Alagille syndrome, a heterogeneous disorder of the notch signaling pathway.* Am J Hum Genet 2006. [PMID 16773578](https://pubmed.ncbi.nlm.nih.gov/16773578/)
4. Gilbert MA, et al. *Alagille syndrome mutation update: Comprehensive overview of JAG1 and NOTCH2 mutation frequencies and insight into missense variant classification.* Hum Mutat 2019. [PMID 31343788](https://pubmed.ncbi.nlm.nih.gov/31343788/)
5. Spinner NB, et al. *Jagged1 mutations in alagille syndrome.* Hum Mutat 2001. [PMID 11139239](https://pubmed.ncbi.nlm.nih.gov/11139239/)
6. Warthen DM, et al. *Jagged1 (JAG1) mutations in Alagille syndrome: increasing the mutation detection rate.* Hum Mutat 2006. [PMID 16575836](https://pubmed.ncbi.nlm.nih.gov/16575836/)
7. Kopan R, et al. *The canonical Notch signaling pathway: unfolding the activation mechanism.* Cell 2009. [PMID 19379690](https://pubmed.ncbi.nlm.nih.gov/19379690/)
8. Falo-Sanjuan J, et al. *Decoding the Notch signal.* Dev Growth Differ 2020. [PMID 31886523](https://pubmed.ncbi.nlm.nih.gov/31886523/)
9. Hayeck TJ, et al. *Likelihood-based calibration improves the clinical utility of JAG1 functional data for variant classification.* Am J Hum Genet 2026. [PMID 42442366](https://pubmed.ncbi.nlm.nih.gov/42442366/)
10. Leonard LD, et al. *Clinical utility gene card for: Alagille Syndrome (ALGS).* Eur J Hum Genet 2014. [PMID 23881058](https://pubmed.ncbi.nlm.nih.gov/23881058/)
11. Adams JM, et al. *The Roles of Notch Signaling in Liver Development and Disease.* Biomolecules 2019. [PMID 31615106](https://pubmed.ncbi.nlm.nih.gov/31615106/)
12. Kamath BM, et al. *Renal involvement and the role of Notch signalling in Alagille syndrome.* Nat Rev Nephrol 2013. [PMID 23752887](https://pubmed.ncbi.nlm.nih.gov/23752887/)
13. Turnpenny PD, et al. *Alagille syndrome: pathogenesis, diagnosis and management.* Eur J Hum Genet 2012. [PMID 21934706](https://pubmed.ncbi.nlm.nih.gov/21934706/)


## 2. Ductal plate malformation and postnatal biliary repair

14. Libbrecht L, et al. *Peripheral bile duct paucity and cholestasis in the liver of a patient with Alagille syndrome: further evidence supporting a lack of postnatal bile duct branching and elongation.* Am J Surg Pathol 2005. [PMID 15897750](https://pubmed.ncbi.nlm.nih.gov/15897750/)
15. Fabris L, et al. *Analysis of liver repair mechanisms in Alagille syndrome and biliary atresia reveals a role for notch signaling.* Am J Pathol 2007. [PMID 17600123](https://pubmed.ncbi.nlm.nih.gov/17600123/)
16. Lemaigre FP, et al. *Mechanisms of liver development: concepts for understanding liver disorders and design of novel therapies.* Gastroenterology 2009. [PMID 19328801](https://pubmed.ncbi.nlm.nih.gov/19328801/)
17. Antoniou A, et al. *Intrahepatic bile ducts develop according to a new mode of tubulogenesis regulated by the transcription factor SOX9.* Gastroenterology 2009. [PMID 19403103](https://pubmed.ncbi.nlm.nih.gov/19403103/)
18. Sparks EE, et al. *Notch signaling regulates formation of the three-dimensional architecture of intrahepatic bile ducts in mice.* Hepatology 2010. [PMID 20069650](https://pubmed.ncbi.nlm.nih.gov/20069650/)
19. Zong Y, et al. *Notch signaling controls liver development by regulating biliary differentiation.* Development 2009. [PMID 19369401](https://pubmed.ncbi.nlm.nih.gov/19369401/)
20. Falix FA, et al. *Hepatic Notch2 deficiency leads to bile duct agenesis perinatally and secondary bile duct formation after weaning.* Dev Biol 2014. [PMID 25446530](https://pubmed.ncbi.nlm.nih.gov/25446530/)
21. Andersson ER, et al. *Mouse Model of Alagille Syndrome and Mechanisms of Jagged1 Missense Mutations.* Gastroenterology 2018. [PMID 29162437](https://pubmed.ncbi.nlm.nih.gov/29162437/)
22. Raynaud P, et al. *Biliary differentiation and bile duct morphogenesis in development and disease.* Int J Biochem Cell Biol 2011. [PMID 19735739](https://pubmed.ncbi.nlm.nih.gov/19735739/)
23. Shneider BL, et al. *Plasma proteome correlations with liver stiffness in pediatric cholestasis implicate epithelial to mesenchymal transition.* Hepatol Commun 2025. [PMID 41021277](https://pubmed.ncbi.nlm.nih.gov/41021277/)
24. Russo P, et al. *Key Histopathologic Features of Liver Biopsies That Distinguish Biliary Atresia From Other Causes of Infantile Cholestasis and Their Correlation With Outcome: A Multicenter Study.* Am J Surg Pathol 2016. [PMID 27776008](https://pubmed.ncbi.nlm.nih.gov/27776008/)


## 3. Natural history — the GALA cohort and long-term outcome

25. Kamath BM, et al. *Outcomes of Childhood Cholestasis in Alagille Syndrome: Results of a Multicenter Observational Study.* Hepatol Commun 2020. [PMID 33313463](https://pubmed.ncbi.nlm.nih.gov/33313463/)
26. Emerick KM, et al. *Features of Alagille syndrome in 92 patients: frequency and relation to prognosis.* Hepatology 1999. [PMID 10051485](https://pubmed.ncbi.nlm.nih.gov/10051485/)
27. Lykavieris P, et al. *Outcome of liver disease in children with Alagille syndrome: a study of 163 patients.* Gut 2001. [PMID 11511567](https://pubmed.ncbi.nlm.nih.gov/11511567/)
28. Hoffenberg EJ, et al. *Outcome of syndromic paucity of interlobular bile ducts (Alagille syndrome) with onset of cholestasis in infancy.* J Pediatr 1995. [PMID 7636645](https://pubmed.ncbi.nlm.nih.gov/7636645/)
29. Quiros-Tejeira RE, et al. *Variable morbidity in alagille syndrome: a review of 43 cases.* J Pediatr Gastroenterol Nutr 1999. [PMID 10512403](https://pubmed.ncbi.nlm.nih.gov/10512403/)
30. Kamath BM, et al. *Systematic Review: The Epidemiology, Natural History, and Burden of Alagille Syndrome.* J Pediatr Gastroenterol Nutr 2018. [PMID 29543694](https://pubmed.ncbi.nlm.nih.gov/29543694/)
31. Mitchell E, et al. *Alagille Syndrome.* Clin Liver Dis 2018. [PMID 30266153](https://pubmed.ncbi.nlm.nih.gov/30266153/)
32. Perez CFM, et al. *Elevated Serum Bile Acids Predict Poor Liver Outcomes in Children With Alagille Syndrome: Results From the GALA Study Group.* Liver Int 2025. [PMID 41250932](https://pubmed.ncbi.nlm.nih.gov/41250932/)


## 4. Diagnosis, bile duct paucity and biochemical phenotype

33. Nischal KK, et al. *Ocular ultrasound in Alagille syndrome: a new sign.* Ophthalmology 1997. [PMID 9022108](https://pubmed.ncbi.nlm.nih.gov/9022108/)
34. Alagille D, et al. *Hepatic ductular hypoplasia associated with characteristic facies, vertebral malformations, retarded physical, mental, and sexual development, and cardiac murmur.* J Pediatr 1975. [PMID 803282](https://pubmed.ncbi.nlm.nih.gov/803282/)
35. Subramaniam P, et al. *Diagnosis of Alagille syndrome-25 years of experience at King's College Hospital.* J Pediatr Gastroenterol Nutr 2011. [PMID 21119543](https://pubmed.ncbi.nlm.nih.gov/21119543/)
36. Neimark E, et al. *Novel surgical and pharmacological approaches to chronic cholestasis in children: partial external biliary diversion for intractable pruritus and xanthomas in Alagille Syndrome.* J Pediatr Gastroenterol Nutr 2003. [PMID 12593401](https://pubmed.ncbi.nlm.nih.gov/12593401/)
37. Minhaj A, et al. *Abstract: Normal Serum IgG4 in Biopsy Proven IgG4 Related Hypophysitis: A Case Report.* S D Med 2026. [PMID 42526010](https://pubmed.ncbi.nlm.nih.gov/42526010/)
38. Liebe R, et al. *Use of genetic analysis in adult cholestatic liver disease: lessons from progressive paediatric syndromes and cohort studies.* Gut 2026. [PMID 42309808](https://pubmed.ncbi.nlm.nih.gov/42309808/)
39. Xuan W, et al. *Alagille syndrome caused by p.L2014Vfs*10 in NOTCH2: a case report and review of the literature.* J Med Case Rep 2026. [PMID 41928326](https://pubmed.ncbi.nlm.nih.gov/41928326/)


## 5. Bile acid physiology and the enterohepatic circulation

40. Hofmann AF, et al. *Key discoveries in bile acid chemistry and biology and their clinical applications: history of the last eight decades.* J Lipid Res 2014. [PMID 24838141](https://pubmed.ncbi.nlm.nih.gov/24838141/)
41. Trauner M, et al. *Benefits and challenges to therapeutic targeting of bile acid circulation in cholestatic liver disease.* Hepatology 2025. [PMID 40601300](https://pubmed.ncbi.nlm.nih.gov/40601300/)
42. Lan T, et al. *Mouse organic solute transporter alpha deficiency alters FGF15 expression and bile acid metabolism.* J Hepatol 2012. [PMID 22542490](https://pubmed.ncbi.nlm.nih.gov/22542490/)
43. Li T, et al. *Bile Acid Signaling in Metabolic and Inflammatory Diseases and Drug Development.* Pharmacol Rev 2024. [PMID 38977324](https://pubmed.ncbi.nlm.nih.gov/38977324/)
44. Craddock AL, et al. *Expression and transport properties of the human ileal and renal sodium-dependent bile acid transporter.* Am J Physiol 1998. [PMID 9458785](https://pubmed.ncbi.nlm.nih.gov/9458785/)
45. Oelkers P, et al. *Primary bile acid malabsorption caused by mutations in the ileal sodium-dependent bile acid transporter gene (SLC10A2).* J Clin Invest 1997. [PMID 9109432](https://pubmed.ncbi.nlm.nih.gov/9109432/)
46. Hofmann AF, et al. *Bile acids: chemistry, pathochemistry, biology, pathobiology, and therapeutics.* Cell Mol Life Sci 2008. [PMID 18488143](https://pubmed.ncbi.nlm.nih.gov/18488143/)
47. Dawson PA, et al. *Targeted deletion of the ileal bile acid transporter eliminates enterohepatic cycling of bile acids in mice.* J Biol Chem 2003. [PMID 12819193](https://pubmed.ncbi.nlm.nih.gov/12819193/)
48. Ridlon JM, et al. *Consequences of bile salt biotransformations by intestinal bacteria.* Gut Microbes 2016. [PMID 26939849](https://pubmed.ncbi.nlm.nih.gov/26939849/)
49. HOFMANN AF, et al. *THE INTRALUMINAL PHASE OF FAT DIGESTION IN MAN: THE LIPID CONTENT OF THE MICELLAR AND OIL PHASES OF INTESTINAL CONTENT OBTAINED DURING FAT DIGESTION AND ABSORPTION.* J Clin Invest 1964. [PMID 14162533](https://pubmed.ncbi.nlm.nih.gov/14162533/)
50. Carey MC, et al. *The characteristics of mixed micellar solutions with particular reference to bile.* Am J Med 1970. [PMID 4924587](https://pubmed.ncbi.nlm.nih.gov/4924587/)
51. Angelin B, et al. *Hepatic uptake of bile acids in man. Fasting and postprandial concentrations of individual bile acids in portal venous and systemic blood serum.* J Clin Invest 1982. [PMID 7119112](https://pubmed.ncbi.nlm.nih.gov/7119112/)


## 6. FXR, SHP and the FGF19 feedback arm

52. Holt JA, et al. *Definition of a novel growth factor-dependent signal cascade for the suppression of bile acid biosynthesis.* Genes Dev 2003. [PMID 12815072](https://pubmed.ncbi.nlm.nih.gov/12815072/)
53. Goodwin B, et al. *A regulatory cascade of the nuclear receptors FXR, SHP-1, and LRH-1 represses bile acid biosynthesis.* Mol Cell 2000. [PMID 11030332](https://pubmed.ncbi.nlm.nih.gov/11030332/)
54. Lu TT, et al. *Molecular basis for feedback regulation of bile acid synthesis by nuclear receptors.* Mol Cell 2000. [PMID 11030331](https://pubmed.ncbi.nlm.nih.gov/11030331/)
55. Sauter G, et al. *Serum concentrations of 7alpha-hydroxy-4-cholesten-3-one reflect bile acid synthesis in humans.* Hepatology 1996. [PMID 8707250](https://pubmed.ncbi.nlm.nih.gov/8707250/)
56. Wong MH, et al. *Identification of a mutation in the ileal sodium-dependent bile acid transporter gene that abolishes transport activity.* J Biol Chem 1995. [PMID 7592981](https://pubmed.ncbi.nlm.nih.gov/7592981/)
57. Kim I, et al. *Differential regulation of bile acid homeostasis by the farnesoid X receptor in liver and intestine.* J Lipid Res 2007. [PMID 17720959](https://pubmed.ncbi.nlm.nih.gov/17720959/)
58. Zhang JH, et al. *Potent stimulation of fibroblast growth factor 19 expression in the human ileum by bile acids.* Am J Physiol Gastrointest Liver Physiol 2013. [PMID 23518683](https://pubmed.ncbi.nlm.nih.gov/23518683/)
59. Al-Dury S, et al. *Ileal Bile Acid Transporter Inhibition for the Treatment of Chronic Constipation, Cholestatic Pruritus, and NASH.* Front Pharmacol 2018. [PMID 30186169](https://pubmed.ncbi.nlm.nih.gov/30186169/)


## 7. Hepatocyte bile acid transporters and adaptive responses

60. Trauner M, et al. *Bile salt transporters: molecular characterization, function, and regulation.* Physiol Rev 2003. [PMID 12663868](https://pubmed.ncbi.nlm.nih.gov/12663868/)
61. Stieger B, et al. *Genetic variations of bile salt transporters as predisposing factors for drug-induced cholestasis, intrahepatic cholestasis of pregnancy and therapeutic response of viral hepatitis.* Expert Opin Drug Metab Toxicol 2011. [PMID 21320040](https://pubmed.ncbi.nlm.nih.gov/21320040/)
62. Krishna M, et al. *Ileal Bile Acid Transporter Inhibitors in Cholestasis: Potential for More Than Just Paediatrics?.* Liver Int 2026. [PMID 41645895](https://pubmed.ncbi.nlm.nih.gov/41645895/)
63. Pollheimer MJ, et al. *Lysyl oxidase-like protein 2 (LOXL2) modulates barrier function in cholangiocytes in cholestasis.* J Hepatol 2018. [PMID 29709678](https://pubmed.ncbi.nlm.nih.gov/29709678/)
64. Denk GU, et al. *Multidrug resistance-associated protein 4 is up-regulated in liver but down-regulated in kidney in obstructive cholestasis in the rat.* J Hepatol 2004. [PMID 15030973](https://pubmed.ncbi.nlm.nih.gov/15030973/)
65. Boyer JL, et al. *Upregulation of a basolateral FXR-dependent bile acid efflux transporter OSTalpha-OSTbeta in cholestasis in humans and rodents.* Am J Physiol Gastrointest Liver Physiol 2006. [PMID 16423920](https://pubmed.ncbi.nlm.nih.gov/16423920/)
66. Keitel V, et al. *Expression and localization of hepatobiliary transport proteins in progressive familial intrahepatic cholestasis.* Hepatology 2005. [PMID 15841457](https://pubmed.ncbi.nlm.nih.gov/15841457/)
67. Halilbasic E, et al. *Bile acid transporters and regulatory nuclear receptors in the liver and beyond.* J Hepatol 2013. [PMID 22885388](https://pubmed.ncbi.nlm.nih.gov/22885388/)


## 8. Maralixibat — trials and long-term outcome

68. Gonzales E, et al. *Efficacy and safety of maralixibat treatment in patients with Alagille syndrome and cholestatic pruritus (ICONIC): a randomised phase 2 study.* Lancet 2021. [PMID 34755627](https://pubmed.ncbi.nlm.nih.gov/34755627/)
69. Shneider BL, et al. *Placebo-Controlled Randomized Trial of an Intestinal Bile Salt Transport Inhibitor for Pruritus in Alagille Syndrome.* Hepatol Commun 2018. [PMID 30288474](https://pubmed.ncbi.nlm.nih.gov/30288474/)
70. Hansen BE, et al. *Event-free survival of maralixibat-treated patients with Alagille syndrome compared to a real-world cohort from GALA.* Hepatology 2024. [PMID 38146932](https://pubmed.ncbi.nlm.nih.gov/38146932/)
71. Kamath BM, et al. *Potential of ileal bile acid transporter inhibition as a therapeutic target in Alagille syndrome and progressive familial intrahepatic cholestasis.* Liver Int 2020. [PMID 32492754](https://pubmed.ncbi.nlm.nih.gov/32492754/)
72. Himes R, et al. *Real-world experience of maralixibat in Alagille syndrome: Novel findings outside of clinical trials.* J Pediatr Gastroenterol Nutr 2024. [PMID 38334237](https://pubmed.ncbi.nlm.nih.gov/38334237/)
73. Loomes KM, et al. *Maralixibat for the treatment of PFIC: Long-term, IBAT inhibition in an open-label, Phase 2 study.* Hepatol Commun 2022. [PMID 35507739](https://pubmed.ncbi.nlm.nih.gov/35507739/)
74. *.*  2026. [PMID 42330149](https://pubmed.ncbi.nlm.nih.gov/42330149/)
75. Shirley M, et al. *Maralixibat: First Approval.* Drugs 2022. [PMID 34813049](https://pubmed.ncbi.nlm.nih.gov/34813049/)


## 9. Odevixibat and other IBAT inhibitors

76. Ovchinsky N, et al. *Efficacy and safety of odevixibat in patients with Alagille syndrome (ASSERT): a phase 3, double-blind, randomised, placebo-controlled trial.* Lancet Gastroenterol Hepatol 2024. [PMID 38670135](https://pubmed.ncbi.nlm.nih.gov/38670135/)
77. Kaltenbach E, et al. *Letter to The Editor, Regarding "Indirect Comparison of Maralixibat and Odevixibat for the Treatment of Progressive Familial Intrahepatic Cholestasis" Recently Published by Lacey and Colleagues.* Clin Ther 2026. [PMID 41298181](https://pubmed.ncbi.nlm.nih.gov/41298181/)
78. Thompson RJ, et al. *Interim results from an ongoing, open-label, single-arm trial of odevixibat in progressive familial intrahepatic cholestasis.* JHEP Rep 2023. [PMID 37456676](https://pubmed.ncbi.nlm.nih.gov/37456676/)
79. Baumann U, et al. *Effects of odevixibat on pruritus and bile acids in children with cholestatic liver disease: Phase 2 study.* Clin Res Hepatol Gastroenterol 2021. [PMID 34182185](https://pubmed.ncbi.nlm.nih.gov/34182185/)
80. Gelhorn H, et al. *Psychometric validation of the Worst Itch Numerical Rating Scale (WI-NRS) and other patient-reported outcome measures for assessing severity and impact of pruritus in patients with primary biliary cholangitis.* Orphanet J Rare Dis 2025. [PMID 40745549](https://pubmed.ncbi.nlm.nih.gov/40745549/)
81. Sohal A, et al. *New therapies for primary biliary cholangitis.* Hepatol Int 2026. [PMID 42101777](https://pubmed.ncbi.nlm.nih.gov/42101777/)
82. Karpen SJ, et al. *Use of a Comprehensive 66-Gene Cholestasis Sequencing Panel in 2171 Cholestatic Infants, Children, and Young Adults.* J Pediatr Gastroenterol Nutr 2021. [PMID 33720099](https://pubmed.ncbi.nlm.nih.gov/33720099/)
83. Deeks ED, et al. *Correction to: Odevixibat: First Approval.* Drugs 2021. [PMID 34554439](https://pubmed.ncbi.nlm.nih.gov/34554439/)


## 10. Cholestatic pruritus — mechanism

84. Yang J, et al. *Development of a clinically viable MRGPRX4 inverse agonist for cholestatic itch treatment.* Nat Chem Biol 2026. [PMID 41957282](https://pubmed.ncbi.nlm.nih.gov/41957282/)
85. Kremer AE, et al. *Lysophosphatidic acid is a potential mediator of cholestatic pruritus.* Gastroenterology 2010. [PMID 20546739](https://pubmed.ncbi.nlm.nih.gov/20546739/)
86. Zhang P, et al. *The role of autotaxin in pruritus and preterm labor in early-onset intrahepatic cholestasis of pregnancy.* BMC Pregnancy Childbirth 2026. [PMID 41620689](https://pubmed.ncbi.nlm.nih.gov/41620689/)
87. Beuers U, et al. *Pruritus in cholestasis: facts and fiction.* Hepatology 2014. [PMID 24807046](https://pubmed.ncbi.nlm.nih.gov/24807046/)
88. Bergasa NV, et al. *The pruritus of cholestasis: From bile acids to opiate agonists: Relevant after all these years.* Med Hypotheses 2018. [PMID 29317077](https://pubmed.ncbi.nlm.nih.gov/29317077/)
89. Alemi F, et al. *The TGR5 receptor mediates bile acid-induced itch and analgesia.* J Clin Invest 2013. [PMID 23524965](https://pubmed.ncbi.nlm.nih.gov/23524965/)
90. Thornton JR, et al. *Opioid peptides and primary biliary cirrhosis.* BMJ 1988. [PMID 3147046](https://pubmed.ncbi.nlm.nih.gov/3147046/)
91. Invernizzi P, et al. *Primary biliary cholangitis (PBC): evolving approaches and expert perspectives.* Expert Rev Gastroenterol Hepatol 2026. [PMID 41399111](https://pubmed.ncbi.nlm.nih.gov/41399111/)
92. Kremer AE, et al. *Receptors, cells and circuits involved in pruritus of systemic disorders.* Biochim Biophys Acta 2014. [PMID 24568861](https://pubmed.ncbi.nlm.nih.gov/24568861/)
93. Romeo M, et al. *Peroxisome Proliferator-Activated Receptor (PPAR) Agonists in Chronic Liver Diseases: Translating Mechanistic Insights into Clinical Practice and Future Perspectives.* Cells 2026. [PMID 42505401](https://pubmed.ncbi.nlm.nih.gov/42505401/)
94. Lieu T, et al. *The bile acid receptor TGR5 activates the TRPA1 channel to induce itch in mice.* Gastroenterology 2014. [PMID 25194674](https://pubmed.ncbi.nlm.nih.gov/25194674/)


## 11. Antipruritic therapy other than IBAT inhibitors

95. Ghent CN, et al. *Treatment of pruritus in primary biliary cirrhosis with rifampin. Results of a double-blind, crossover, randomized trial.* Gastroenterology 1988. [PMID 3275568](https://pubmed.ncbi.nlm.nih.gov/3275568/)
96. Wolfhagen FH, et al. *Oral naltrexone treatment for cholestatic pruritus: a double-blind, placebo-controlled study.* Gastroenterology 1997. [PMID 9322521](https://pubmed.ncbi.nlm.nih.gov/9322521/)
97. Mayo MJ, et al. *Sertraline as a first-line treatment for cholestatic pruritus.* Hepatology 2007. [PMID 17326161](https://pubmed.ncbi.nlm.nih.gov/17326161/)
98. Bachs L, et al. *Effects of long-term rifampicin administration in primary biliary cirrhosis.* Gastroenterology 1992. [PMID 1587427](https://pubmed.ncbi.nlm.nih.gov/1587427/)
99. Datta DV, et al. *Cholestyramine for long term relief of the pruritus complicating intrahepatic cholestasis.* Gastroenterology 1966. [PMID 5905351](https://pubmed.ncbi.nlm.nih.gov/5905351/)
100. Wu WN, et al. *A Rare Coexistence of Biliary Atresia and Alagille Syndrome in a Neonate: Clinical Implications of Dual Etiology in Neonatal Cholestasis.* Diagnostics (Basel) 2026. [PMID 42351413](https://pubmed.ncbi.nlm.nih.gov/42351413/)


## 12. Surgical biliary diversion and liver transplantation

101. Whitington PF, et al. *Partial external diversion of bile for the treatment of intractable pruritus associated with intrahepatic cholestasis.* Gastroenterology 1988. [PMID 3371608](https://pubmed.ncbi.nlm.nih.gov/3371608/)
102. Modi BP, et al. *Ileal exclusion for refractory symptomatic cholestasis in Alagille syndrome.* J Pediatr Surg 2007. [PMID 17502187](https://pubmed.ncbi.nlm.nih.gov/17502187/)
103. Wang KS, et al. *Analysis of surgical interruption of the enterohepatic circulation as a treatment for pediatric cholestasis.* Hepatology 2017. [PMID 28027587](https://pubmed.ncbi.nlm.nih.gov/28027587/)
104. Kamath BM, et al. *Outcomes of liver transplantation for patients with Alagille syndrome: the studies of pediatric liver transplantation experience.* Liver Transpl 2012. [PMID 22454296](https://pubmed.ncbi.nlm.nih.gov/22454296/)
105. Carlos EC, et al. *Wilms Tumor After Orthotopic Liver Transplant in a Patient With Alagille Syndrome.* Urology 2018. [PMID 29879405](https://pubmed.ncbi.nlm.nih.gov/29879405/)
106. Kamath BM, et al. *Alagille syndrome and liver transplantation.* J Pediatr Gastroenterol Nutr 2010. [PMID 19949348](https://pubmed.ncbi.nlm.nih.gov/19949348/)
107. Liu WW, et al. *(Clinical characteristics and follow-up in children with Alagille syndrome).* Zhonghua Gan Zang Bing Za Zhi 2026. [PMID 42036231](https://pubmed.ncbi.nlm.nih.gov/42036231/)
108. Jarasvaraparn C, et al. *Exploring odevixibat's efficacy in alagille syndrome: insights from recent clinical trials and IBAT inhibitor experiences.* Expert Opin Pharmacother 2024. [PMID 39155775](https://pubmed.ncbi.nlm.nih.gov/39155775/)


## 13. Cholestatic liver injury, fibrosis and portal hypertension

109. Hasan MN, et al. *Gly-βMCA modulates bile acid metabolism to reduce hepatobiliary injury in Mdr2 KO mice.* Am J Physiol Gastrointest Liver Physiol 2025. [PMID 40418643](https://pubmed.ncbi.nlm.nih.gov/40418643/)
110. Woolbright BL, et al. *Therapeutic targets for cholestatic liver injury.* Expert Opin Ther Targets 2016. [PMID 26479335](https://pubmed.ncbi.nlm.nih.gov/26479335/)
111. Tsuchida T, et al. *Mechanisms of hepatic stellate cell activation.* Nat Rev Gastroenterol Hepatol 2017. [PMID 28487545](https://pubmed.ncbi.nlm.nih.gov/28487545/)
112. Pinzani M, et al. *Liver fibrosis: from the bench to clinical targets.* Dig Liver Dis 2004. [PMID 15115333](https://pubmed.ncbi.nlm.nih.gov/15115333/)
113. Penz-Österreicher M, et al. *Fibrosis in autoimmune and cholestatic liver disease.* Best Pract Res Clin Gastroenterol 2011. [PMID 21497742](https://pubmed.ncbi.nlm.nih.gov/21497742/)
114. Sato K, et al. *Ductular Reaction in Liver Diseases: Pathological Mechanisms and Translational Significances.* Hepatology 2019. [PMID 30070383](https://pubmed.ncbi.nlm.nih.gov/30070383/)
115. Berzigotti A, et al. *Effect of meal ingestion on liver stiffness in patients with cirrhosis and portal hypertension.* PLoS One 2013. [PMID 23520531](https://pubmed.ncbi.nlm.nih.gov/23520531/)
116. Fabris L, et al. *Emerging concepts in biliary repair and fibrosis.* Am J Physiol Gastrointest Liver Physiol 2017. [PMID 28526690](https://pubmed.ncbi.nlm.nih.gov/28526690/)


## 14. Malabsorption, fat-soluble vitamins and growth

117. Kamath BM, et al. *Fat Soluble Vitamin Assessment and Supplementation in Cholestasis.* Clin Liver Dis 2022. [PMID 35868689](https://pubmed.ncbi.nlm.nih.gov/35868689/)
118. Sokol RJ, et al. *Frequency and clinical progression of the vitamin E deficiency neurologic disorder in children with prolonged neonatal cholestasis.* Am J Dis Child 1985. [PMID 4061425](https://pubmed.ncbi.nlm.nih.gov/4061425/)
119. Feranchak AP, et al. *Comparison of indices of vitamin A status in children with chronic liver disease.* Hepatology 2005. [PMID 16175620](https://pubmed.ncbi.nlm.nih.gov/16175620/)
120. Shneider BL, et al. *Efficacy of fat-soluble vitamin supplementation in infants with biliary atresia.* Pediatrics 2012. [PMID 22891232](https://pubmed.ncbi.nlm.nih.gov/22891232/)
121. Huysentruyt K, et al. *Condition-Specific Growth Charts for Children With Alagille Syndrome.* JAMA Netw Open 2025. [PMID 41284294](https://pubmed.ncbi.nlm.nih.gov/41284294/)
122. Quiros-Tejeira RE, et al. *Does liver transplantation affect growth pattern in Alagille syndrome?.* Liver Transpl 2000. [PMID 10980057](https://pubmed.ncbi.nlm.nih.gov/10980057/)
123. Bucuvalas JC, et al. *Growth hormone insensitivity associated with elevated circulating growth hormone-binding protein in children with Alagille syndrome and short stature.* J Clin Endocrinol Metab 1993. [PMID 8501153](https://pubmed.ncbi.nlm.nih.gov/8501153/)
124. Rovner AJ, et al. *Rethinking growth failure in Alagille syndrome: the role of dietary intake and steatorrhea.* J Pediatr Gastroenterol Nutr 2002. [PMID 12394373](https://pubmed.ncbi.nlm.nih.gov/12394373/)
125. Ramaccioni V, et al. *Nutritional aspects of chronic liver disease and liver transplantation in children.* J Pediatr Gastroenterol Nutr 2000. [PMID 10776944](https://pubmed.ncbi.nlm.nih.gov/10776944/)


## 15. Cholestatic dyslipidaemia, lipoprotein X and xanthoma

126. Seidel D, et al. *A lipoprotein characterizing obstructive jaundice. II. Isolation and partial characterization of the protein moieties of low density lipoproteins.* J Clin Invest 1970. [PMID 5480863](https://pubmed.ncbi.nlm.nih.gov/5480863/)
127. Nemes K, et al. *Cholesterol metabolism in cholestatic liver disease and liver transplantation: From molecular mechanisms to clinical implications.* World J Hepatol 2016. [PMID 27574546](https://pubmed.ncbi.nlm.nih.gov/27574546/)
128. Fellin R, et al. *Lipoprotein-X fifty years after its original discovery.* Nutr Metab Cardiovasc Dis 2019. [PMID 30503707](https://pubmed.ncbi.nlm.nih.gov/30503707/)
129. Davit-Spraul A, et al. *Abnormal lipoprotein pattern in patients with Alagille syndrome depends on Icterus severity.* Gastroenterology 1996. [PMID 8831598](https://pubmed.ncbi.nlm.nih.gov/8831598/)
130. Garcia MA, et al. *Alagille syndrome: cutaneous manifestations in 38 children.* Pediatr Dermatol 2005. [PMID 15660889](https://pubmed.ncbi.nlm.nih.gov/15660889/)


## 16. Cardiac and vascular disease in Alagille syndrome

131. Kulikauskas MR, et al. *Three high throughput compatible cell-based assays for identifying small molecule JAG1 upregulators for Alagille syndrome.* SLAS Discov 2026. [PMID 42498207](https://pubmed.ncbi.nlm.nih.gov/42498207/)
132. Kamath BM, et al. *Vascular anomalies in Alagille syndrome: a significant cause of morbidity and mortality.* Circulation 2004. [PMID 14993126](https://pubmed.ncbi.nlm.nih.gov/14993126/)
133. Kamath BM, et al. *Renal anomalies in Alagille syndrome: a disease-defining feature.* Am J Med Genet A 2012. [PMID 22105858](https://pubmed.ncbi.nlm.nih.gov/22105858/)
134. Emerick KM, et al. *Intracranial vascular abnormalities in patients with Alagille syndrome.* J Pediatr Gastroenterol Nutr 2005. [PMID 15990638](https://pubmed.ncbi.nlm.nih.gov/15990638/)
135. Woolfenden AR, et al. *Moyamoya syndrome in children with Alagille syndrome: additional evidence of a vasculopathy.* Pediatrics 1999. [PMID 9925853](https://pubmed.ncbi.nlm.nih.gov/9925853/)
136. Mohiaddin H, et al. *Does severe pulmonary hypertension affect long-term survival after TAVI for severe aortic stenosis?.* Br J Cardiol 2026. [PMID 42524360](https://pubmed.ncbi.nlm.nih.gov/42524360/)
137. Ragheb DK, et al. *Durability of Aortic Homografts in Pulmonary Atresia and Major Aortopulmonary Collateral Arteries.* World J Pediatr Congenit Heart Surg 2024. [PMID 39166263](https://pubmed.ncbi.nlm.nih.gov/39166263/)
138. Salem JE, et al. *Hypertension and aortorenal disease in Alagille syndrome.* J Hypertens 2012. [PMID 22525199](https://pubmed.ncbi.nlm.nih.gov/22525199/)
139. High FA, et al. *The multifaceted role of Notch in cardiac development and disease.* Nat Rev Genet 2008. [PMID 18071321](https://pubmed.ncbi.nlm.nih.gov/18071321/)


## 17. Renal, skeletal and ocular involvement

140. Sanderson E, et al. *Vertebral anomalies in children with Alagille syndrome: an analysis of 50 consecutive patients.* Pediatr Radiol 2002. [PMID 11819079](https://pubmed.ncbi.nlm.nih.gov/11819079/)
141. Bales CB, et al. *Pathologic lower extremity fractures in children with Alagille syndrome.* J Pediatr Gastroenterol Nutr 2010. [PMID 20453673](https://pubmed.ncbi.nlm.nih.gov/20453673/)
142. Hingorani M, et al. *Ocular abnormalities in Alagille syndrome.* Ophthalmology 1999. [PMID 9951486](https://pubmed.ncbi.nlm.nih.gov/9951486/)
143. Varol Fİ, et al. *Alagille syndrome case series: five new variants and two large deletions.* Eur J Pediatr 2026. [PMID 42377562](https://pubmed.ncbi.nlm.nih.gov/42377562/)
144. Kohut TJ, et al. *Alagille Syndrome: A Focused Review on Clinical Features, Genetics, and Treatment.* Semin Liver Dis 2021. [PMID 34215014](https://pubmed.ncbi.nlm.nih.gov/34215014/)
145. Ayoub MD, et al. *Alagille Syndrome: Current Understanding of Pathogenesis, and Challenges in Diagnosis and Management.* Clin Liver Dis 2022. [PMID 35868679](https://pubmed.ncbi.nlm.nih.gov/35868679/)


## 18. Modelling methods — QSP, bile acid and hazard models

146. Baier V, et al. *A Physiology-Based Model of Human Bile Acid Metabolism for Predicting Bile Acid Tissue Levels After Drug Administration in Healthy Subjects and BRIC Type 2 Patients.* Front Physiol 2019. [PMID 31611804](https://pubmed.ncbi.nlm.nih.gov/31611804/)
147. Zhang Y, et al. *Effect of bile duct ligation on bile acid composition in mouse serum and liver.* Liver Int 2012. [PMID 22098667](https://pubmed.ncbi.nlm.nih.gov/22098667/)
148. Brathovde M, et al. *A lean additive frailty model: With an application to clustering of melanoma in Norwegian families.* Stat Med 2023. [PMID 37527835](https://pubmed.ncbi.nlm.nih.gov/37527835/)
149. Bjarnason H, et al. *Fisher information for two gamma frailty bivariate Weibull models.* Lifetime Data Anal 2000. [PMID 10763561](https://pubmed.ncbi.nlm.nih.gov/10763561/)
150. Balan TA, et al. *A tutorial on frailty models.* Stat Methods Med Res 2020. [PMID 32466712](https://pubmed.ncbi.nlm.nih.gov/32466712/)
151. Buchsbaum M, et al. *Neural events and psychophysical law.* Science 1971. [PMID 5550509](https://pubmed.ncbi.nlm.nih.gov/5550509/)
152. Carstens E, et al. *Brain Processing of Itch and Scratching.*  2014. [PMID 24830004](https://pubmed.ncbi.nlm.nih.gov/24830004/)
153. Musante CJ, et al. *Quantitative Systems Pharmacology: A Case for Disease Models.* Clin Pharmacol Ther 2017. [PMID 27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)
154. Goryanin I, et al. *Validation of AI-enabled surrogate models in quantitative systems pharmacology: a practical, context-of-use-driven review.* Drug Discov Today 2026. [PMID 42409163](https://pubmed.ncbi.nlm.nih.gov/42409163/)


---

## Sources deliberately NOT cited as evidence

- **Regulatory label text and press releases.** Approval wording for
  maralixibat and odevixibat is indication language, not data, and the
  trials behind it are cited directly above.
- **Any number for the ASSERT population's mean baseline total bilirubin.**
  The model needs it (it anchors the duct-capacity→bilirubin map that
  produces the comparison with the GALA cut-points) and could not find it
  in the primary publication's abstract. It is therefore an assumption,
  swept in scenario 6b, and flagged as `FAILURE 2`. The convergence with
  GALA's 5.0 and 10.0 mg/dL thresholds should be read with that in mind —
  only the *ratio* of the two model boundaries is anchor-independent.
- **Serial-biopsy duct:portal-tract ratio series.** The postnatal ductular
  repair rate would be identifiable from one and is not identifiable
  without one; no such series appears to exist, so that rate is descriptive
  (`FAILURE 4`).

## Disclaimer

Research and teaching model. Not validated for clinical use, not a decision
support tool, and no parameter in it should be used to manage a patient.

