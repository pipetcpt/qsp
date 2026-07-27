# Autoimmune Pulmonary Alveolar Proteinosis (aPAP) — References

158 numbered entries supporting the mechanistic map (`pap_qsp_model.dot`), the
mrgsolve model (`pap_mrgsolve_model.R`) and the Shiny dashboard
(`pap_shiny_app.R`).

**Every PMID below was resolved through the NCBI E-utilities API and the
returned title, journal and year checked against the citation as written.**
This is not a formality. A first draft of this file cited the PAGE trial as
PMID 31483961. That identifier resolves to a paper on a different subject
entirely (Gunduz-Bruce, *N Engl J Med* 2019); the PAGE trial is 31483963. A
second recalled PMID, for the rituximab series, turned out to belong to an
unrelated respiratory paper and was removed rather than guessed at again. Two
errors in two attempts is the reason nothing here is cited unverified.

Links are `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`.

---

## 1. The autoantibody: discovery, and why it is a buffer rather than an inhibitor

The quantitative core of the model is here. Alveolar GM-CSF sits at tens of
pg/mL while neutralising antibody sits at µg/mL, so free ligand is set by a
binding equilibrium in the excess-antibody regime. That single fact produces
both the titre threshold and the absence of a titre–severity relationship
above it.

1. Kitamura T, et al. Serological diagnosis of idiopathic pulmonary alveolar proteinosis. *Am J Respir Crit Care Med* 2000. [PMID 10934102](https://pubmed.ncbi.nlm.nih.gov/10934102/)
2. Uchida K, et al. High-affinity autoantibodies specifically eliminate granulocyte-macrophage colony-stimulating factor activity in the lungs of patients with idiopathic pulmonary alveolar proteinosis. *Blood* 2004. [PMID 14512323](https://pubmed.ncbi.nlm.nih.gov/14512323/) — the paper the buffer module is built on: the antibody removes GM-CSF *activity* in the lung.
3. Uchida K, et al. Granulocyte/macrophage-colony-stimulating factor autoantibodies and myeloid cell immune functions in healthy subjects. *Blood* 2009. [PMID 19282464](https://pubmed.ncbi.nlm.nih.gov/19282464/) — low-level GMAb in healthy people; the basis for the sub-threshold region of the titration curve.
4. Uchida K, et al. Response: Granulocyte/macrophage colony-stimulating factor autoantibodies and myeloid cell immune functions in healthy persons. *Blood* 2010. [PMID 20075172](https://pubmed.ncbi.nlm.nih.gov/20075172/)
5. Sakagami T, et al. Patient-derived granulocyte/macrophage colony-stimulating factor autoantibodies reproduce pulmonary alveolar proteinosis in nonhuman primates. *Am J Respir Crit Care Med* 2010. [PMID 20224064](https://pubmed.ncbi.nlm.nih.gov/20224064/) — the 5 µg/mL critical threshold, identical in lung and blood, and the anchor for diagnostic D04.
6. Kusakabe Y, et al. A standardized blood test for the routine clinical diagnosis of impaired GM-CSF signaling using flow cytometry. *J Immunol Methods* 2014. [PMID 25068538](https://pubmed.ncbi.nlm.nih.gov/25068538/)
7. Bradhurst P, et al. A bioluminescent inhibition immunoassay for detecting GM-CSF inhibitory activity in serum. *Clin Exp Immunol* 2026. [PMID 42132313](https://pubmed.ncbi.nlm.nih.gov/42132313/) — directly relevant to D05: an assay of *neutralisation* rather than binding is what the model predicts should correlate with severity.
8. Salvator H, et al. Neutralizing GM-CSF autoantibodies in pulmonary alveolar proteinosis, cryptococcal meningitis and severe nocardiosis. *Respir Res* 2022. [PMID 36221098](https://pubmed.ncbi.nlm.nih.gov/36221098/)
9. Hirose M, et al. B cell-activating factors in autoimmune pulmonary alveolar proteinosis. *Orphanet J Rare Dis* 2021. [PMID 33653382](https://pubmed.ncbi.nlm.nih.gov/33653382/)
10. Akasaka K, et al. Cytokine profiles associated with disease severity and prognosis of autoimmune pulmonary alveolar proteinosis. *Respir Investig* 2024. [PMID 38705133](https://pubmed.ncbi.nlm.nih.gov/38705133/)
11. Fujishiro E, et al. Comorbidity of autoimmune diseases in patients with autoimmune pulmonary alveolar proteinosis. *Respir Investig* 2026. [PMID 42250564](https://pubmed.ncbi.nlm.nih.gov/42250564/)

## 2. GM-CSF, its receptor, and the alveolar macrophage differentiation program

The model treats the disease as loss of a *signal*: receptor occupancy drives
PU.1 → PPARγ → lysosomal catabolic machinery, with a floor that is
GM-CSF-independent.

12. Dranoff G, et al. Involvement of granulocyte-macrophage colony-stimulating factor in pulmonary homeostasis. *Science* 1994. [PMID 8171324](https://pubmed.ncbi.nlm.nih.gov/8171324/)
13. Stanley E, et al. Granulocyte/macrophage colony-stimulating factor-deficient mice show no major perturbation of hematopoiesis but develop a characteristic pulmonary pathology. *Proc Natl Acad Sci USA* 1994. [PMID 8202532](https://pubmed.ncbi.nlm.nih.gov/8202532/) — the two 1994 knockouts that made PAP a signalling disease rather than a secretion disease.
14. Dranoff G, Mulligan RC. Activities of granulocyte-macrophage colony-stimulating factor revealed by gene transfer and gene knockout studies. *Stem Cells* 1994. [PMID 7696961](https://pubmed.ncbi.nlm.nih.gov/7696961/)
15. Shibata Y, et al. GM-CSF regulates alveolar macrophage differentiation and innate immunity in the lung through PU.1. *Immunity* 2001. [PMID 11672538](https://pubmed.ncbi.nlm.nih.gov/11672538/) — PU.1 as the master regulator; the transducer in `$ODE`.
16. Bonfield TL, et al. Peroxisome proliferator-activated receptor-gamma is deficient in alveolar macrophages from patients with alveolar proteinosis. *Am J Respir Cell Mol Biol* 2003. [PMID 12805087](https://pubmed.ncbi.nlm.nih.gov/12805087/)
17. Bonfield TL, et al. Peroxisome proliferator-activated receptor-gamma regulates the expression of alveolar macrophage macrophage colony-stimulating factor. *J Immunol* 2008. [PMID 18566389](https://pubmed.ncbi.nlm.nih.gov/18566389/)
18. Baker AD, et al. Targeted PPARγ deficiency in alveolar macrophages disrupts surfactant catabolism. *J Lipid Res* 2010. [PMID 20064973](https://pubmed.ncbi.nlm.nih.gov/20064973/) — PPARγ loss alone impairs catabolism: the model's link from signalling to digestive capacity.
19. Malur A, et al. Restoration of PPARγ reverses lipid accumulation in alveolar macrophages of GM-CSF knockout mice. *Am J Physiol Lung Cell Mol Physiol* 2011. [PMID 21036914](https://pubmed.ncbi.nlm.nih.gov/21036914/)
20. Gschwend J, et al. Alveolar macrophages rely on GM-CSF from alveolar epithelial type 2 cells before and after birth. *J Exp Med* 2021. [PMID 34431978](https://pubmed.ncbi.nlm.nih.gov/34431978/) — AEC2 as the source of the alveolar GM-CSF the model produces at a constant rate.
21. Trapnell BC, Whitsett JA. Pulmonary alveolar proteinosis, a primary immunodeficiency of impaired GM-CSF stimulation of macrophages. *Curr Opin Immunol* 2009. [PMID 19796925](https://pubmed.ncbi.nlm.nih.gov/19796925/)
22. Westermann F, et al. Adenophages are an atypical macrophage population in exocrine glands sustained by ILC2-derived GM-CSF. *Nat Immunol* 2026. [PMID 41461985](https://pubmed.ncbi.nlm.nih.gov/41461985/)
23. Boussetta T, et al. The peptidyl-prolyl isomerase Pin1 controls GM-CSF-induced priming of NADPH oxidase in human neutrophils. *Int Immunopharmacol* 2024. [PMID 38851160](https://pubmed.ncbi.nlm.nih.gov/38851160/) — GM-CSF priming of neutrophils, the host-defence arm of the map.

## 3. The cholesterol arm

Why the macrophage stalls: it is the cholesterol, not the phospholipid, that
it cannot dispose of. This is the module the statin acts on, and the reason a
systemically delivered agent has leverage an aerosol does not.

24. Thomassen MJ, et al. ABCG1 is deficient in alveolar macrophages of GM-CSF knockout mice and patients with pulmonary alveolar proteinosis. *J Lipid Res* 2007. [PMID 17848583](https://pubmed.ncbi.nlm.nih.gov/17848583/)
25. McCarthy C, et al. Statin as a novel pharmacotherapy of pulmonary alveolar proteinosis. *Nat Commun* 2018. [PMID 30087322](https://pubmed.ncbi.nlm.nih.gov/30087322/) — the mechanism and the first clinical observations behind scenario S26.
26. Huang J, et al. Causal role of lipid metabolism in pulmonary alveolar proteinosis: an observational and mendelian randomisation study. *Thorax* 2024. [PMID 38124156](https://pubmed.ncbi.nlm.nih.gov/38124156/)
27. Shi S, et al. Assessment of Statin Treatment for Pulmonary Alveolar Proteinosis without Hypercholesterolemia: A 12-Month Prospective, Longitudinal, and Observational Study. *Biomed Res Int* 2022. [PMID 36330458](https://pubmed.ncbi.nlm.nih.gov/36330458/)
28. Shi S, et al. Long-term follow-up and successful treatment of pulmonary alveolar proteinosis without hypercholesterolemia with statin therapy: a case report. *J Int Med Res* 2021. [PMID 33926277](https://pubmed.ncbi.nlm.nih.gov/33926277/)
29. Takano T, et al. A case of autoimmune pulmonary alveolar proteinosis with severe respiratory failure treated with segmental lung lavage and oral statin therapy. *Respir Med Case Rep* 2022. [PMID 35707405](https://pubmed.ncbi.nlm.nih.gov/35707405/)
30. Venkatraman R, et al. Engineered biodegradable polymeric nanoparticles injectable system of atorvastatin for improved therapeutic effect. *Sci Rep* 2025. [PMID 41361564](https://pubmed.ncbi.nlm.nih.gov/41361564/) — atorvastatin disposition, used for the statin PK block.

## 4. Surfactant metabolism: the mass balance, and why recycling is not clearance

The single most consequential structural decision in the model. Type-II
re-uptake and re-secretion is a large flux with zero net effect at steady
state; counting it as clearance predicts a pool that settles at 2–3× normal
instead of the 30–100× that patients reach.

31. Wright JR. Clearance and recycling of pulmonary surfactant. *Am J Physiol* 1990. [PMID 2200279](https://pubmed.ncbi.nlm.nih.gov/2200279/) — the review that separates the recycling loop from true catabolism.
32. Bruni R, et al. Postnatal transformations of alveolar surfactant in the rabbit: changes in pool size, pool morphology and isoforms. *Biochim Biophys Acta* 1988. [PMID 3337839](https://pubmed.ncbi.nlm.nih.gov/3337839/)
33. Vedovelli L, et al. Simultaneous measurement of phosphatidylglycerol and disaturated-phosphatidylcholine palmitate kinetics from alveolar surfactant. *J Mass Spectrom* 2011. [PMID 22012664](https://pubmed.ncbi.nlm.nih.gov/22012664/) — stable-isotope surfactant kinetics; the basis for the turnover rate.
34. Cogo PE, et al. Surfactant disaturated-phosphatidylcholine kinetics in acute respiratory distress syndrome by stable isotopes and a two compartment model. *Respir Res* 2007. [PMID 17313681](https://pubmed.ncbi.nlm.nih.gov/17313681/)
35. Savov J, et al. Incorporation of biotinylated SP-A into rat lung surfactant layer, type II cells, and clara cells. *Am J Physiol Lung Cell Mol Physiol* 2000. [PMID 10893210](https://pubmed.ncbi.nlm.nih.gov/10893210/)
36. Pinto RA, et al. Association of surfactant protein C with isolated alveolar type II cells. *Biochim Biophys Acta* 1995. [PMID 7893733](https://pubmed.ncbi.nlm.nih.gov/7893733/)
37. Uy CC, et al. Granulocyte-macrophage colony-stimulating factor increases surfactant phospholipid in premature rabbits. *Pediatr Res* 1999. [PMID 10541327](https://pubmed.ncbi.nlm.nih.gov/10541327/)
38. Herbein JF, Wright JR. Enhanced clearance of surfactant protein D during LPS-induced acute inflammation in rat lung. *Am J Physiol Lung Cell Mol Physiol* 2001. [PMID 11404270](https://pubmed.ncbi.nlm.nih.gov/11404270/)
39. Phelps DS, et al. Comparison of the Toponomes of Alveolar Macrophages From Wild Type and Surfactant Protein A Knockout Mice. *Front Immunol* 2022. [PMID 35572576](https://pubmed.ncbi.nlm.nih.gov/35572576/)
40. Aono Y, et al. Surfactant protein-D regulates effector cell function and fibrotic lung remodeling in response to bleomycin injury. *Am J Respir Crit Care Med* 2012. [PMID 22198976](https://pubmed.ncbi.nlm.nih.gov/22198976/)
41. Laycock A, et al. Knock-out mouse models and single particle ICP-MS reveal that SP-D and SP-A deficiency reduces agglomeration of inhaled particles. *Nanotoxicology* 2025. [PMID 39868723](https://pubmed.ncbi.nlm.nih.gov/39868723/)

## 5. Congenital surfactant dysfunction (the differential the model does not treat as PAP)

42. Hamvas A. Inherited surfactant protein-B deficiency and surfactant protein-C associated disease: clinical features and evaluation. *Semin Perinatol* 2006. [PMID 17142157](https://pubmed.ncbi.nlm.nih.gov/17142157/)
43. Peca D, et al. Altered surfactant homeostasis and recurrent respiratory failure secondary to TTF-1 nuclear targeting defect. *Respir Res* 2011. [PMID 21867529](https://pubmed.ncbi.nlm.nih.gov/21867529/)
44. Beverstock AM, et al. Progressive respiratory failure in a term neonate with ABCA3 surfactant deficiency. *J Neonatal Perinatal Med* 2026. [PMID 40852890](https://pubmed.ncbi.nlm.nih.gov/40852890/)
45. Kang MH, et al. A lung tropic AAV vector improves survival in a mouse model of surfactant B deficiency. *Nat Commun* 2020. [PMID 32764559](https://pubmed.ncbi.nlm.nih.gov/32764559/)
46. Barnett RC, et al. Electroporation-mediated gene delivery of surfactant protein B (SP-B) restores expression and improves survival in mouse model of SP-B deficiency. *Exp Biol Med* 2017. [PMID 28581337](https://pubmed.ncbi.nlm.nih.gov/28581337/)

## 6. Hereditary PAP: the receptor mutations that make GM-CSF a structural null

Diagnostic D13 requires inhaled GM-CSF to do *nothing* when the receptor is
absent. These are the papers that make that a real requirement.

47. Suzuki T, et al. Hereditary pulmonary alveolar proteinosis: pathogenesis, presentation, diagnosis, and therapy. *Am J Respir Crit Care Med* 2010. [PMID 20622029](https://pubmed.ncbi.nlm.nih.gov/20622029/)
48. Tanaka T, et al. Adult-onset hereditary pulmonary alveolar proteinosis caused by a single-base deletion in CSF2RB. *J Med Genet* 2011. [PMID 21075760](https://pubmed.ncbi.nlm.nih.gov/21075760/)
49. Shima K, et al. A murine model of hereditary pulmonary alveolar proteinosis caused by homozygous Csf2ra gene disruption. *Am J Physiol Lung Cell Mol Physiol* 2022. [PMID 35043685](https://pubmed.ncbi.nlm.nih.gov/35043685/)
50. Chen Q, et al. A Hereditary Pulmonary Alveolar Proteinosis Caused by a Novel Hemizygous Variation of the CSF2RA Gene. *Mol Genet Genomic Med* 2025. [PMID 41311221](https://pubmed.ncbi.nlm.nih.gov/41311221/)
51. Sivasubramanian D, et al. Hereditary pulmonary alveolar proteinosis in a 5-year-old child: Diagnostic insights and therapeutic approach. *Radiol Case Rep* 2025. [PMID 40486150](https://pubmed.ncbi.nlm.nih.gov/40486150/)
52. Griese M, et al. Autoimmune pulmonary alveolar proteinosis in children. *ERJ Open Res* 2022. [PMID 35350279](https://pubmed.ncbi.nlm.nih.gov/35350279/)

## 7. Cell and gene therapy: replacing capacity instead of restoring signal

Scenario S30 — the only intervention in the model that helps hereditary PAP.

53. Suzuki T, Arumugam P, et al. Long-Term Safety and Efficacy of Gene-Pulmonary Macrophage Transplantation Therapy of PAP in Csf2ra−/− mice. *Mol Ther* 2019. [PMID 31326401](https://pubmed.ncbi.nlm.nih.gov/31326401/)
54. Arumugam P, et al. A toxicology study of Csf2ra complementation and pulmonary macrophage transplantation therapy of hereditary PAP in mice. *Mol Ther Methods Clin Dev* 2024. [PMID 38596536](https://pubmed.ncbi.nlm.nih.gov/38596536/)
55. Hetzel M, et al. Function and Safety of Lentivirus-Mediated Gene Transfer for CSF2RA-Deficiency. *Hum Gene Ther Methods* 2017. [PMID 28854814](https://pubmed.ncbi.nlm.nih.gov/28854814/)
56. Li F, et al. Gene therapy of Csf2ra deficiency in mouse fetal monocyte precursors restores alveolar macrophage development and function. *JCI Insight* 2022. [PMID 35393945](https://pubmed.ncbi.nlm.nih.gov/35393945/)
57. Mishra-Sopori V, et al. Restitutio ad integrum: Rescuing the Alveolar Macrophage Function with HSCT in Pulmonary Alveolar Proteinosis Due to CSF2Rα Deficiency. *J Clin Immunol* 2024. [PMID 39621143](https://pubmed.ncbi.nlm.nih.gov/39621143/)

## 8. Secondary PAP: a cell-number disease with intact signalling

Scenario S29 is built with the autoantibody switched off and monocyte supply
reduced, which is why inhaled GM-CSF has almost nothing left to restore.

58. Lim J, et al. Secondary Pulmonary Alveolar Proteinosis. *Semin Respir Crit Care Med* 2025. [PMID 41213627](https://pubmed.ncbi.nlm.nih.gov/41213627/)
59. Shimaya M, et al. Autoimmune Pulmonary Alveolar Proteinosis Complicated by Myelodysplastic Syndrome. *Intern Med* 2024. [PMID 37839886](https://pubmed.ncbi.nlm.nih.gov/37839886/)
60. Raffáč Š, et al. GATA2 deficiency in an adult with alveolar proteinosis, infections, lymphadenopathy with granulomatosis. *Front Immunol* 2025. [PMID 41322420](https://pubmed.ncbi.nlm.nih.gov/41322420/)
61. Gao WJ, et al. A case of pulmonary alveolar proteinosis secondary to GATA2 deficiency combined with splenic M. kansasii infection. *Zhonghua Jie He He Hu Xi Za Zhi* 2025. [PMID 41362143](https://pubmed.ncbi.nlm.nih.gov/41362143/)
62. Marcu AD, et al. Insights into Pediatric GATA2-Related MDS. *Biomedicines* 2025. [PMID 40299403](https://pubmed.ncbi.nlm.nih.gov/40299403/)
63. Dournes G, et al. CT features of genetic mutation-related pulmonary alveolar proteinosis (CCR2 and GATA2 deficiency). *Diagn Interv Imaging* 2025. [PMID 40340131](https://pubmed.ncbi.nlm.nih.gov/40340131/)
64. Haraguchi M, et al. Disseminated nontuberculous mycobacteriosis and fungemia after second delivery in a patient with MonoMAC syndrome. *BMC Infect Dis* 2021. [PMID 34051752](https://pubmed.ncbi.nlm.nih.gov/34051752/)
65. Katakura T, et al. Secondary Pulmonary Alveolar Proteinosis Complicated by Hemophagocytic Syndrome in a Patient with Adult-onset Still's disease. *Intern Med* 2026. [PMID 40603093](https://pubmed.ncbi.nlm.nih.gov/40603093/)
66. Torén K, et al. Occupational exposures and risk of pulmonary alveolar proteinosis (PAP). *Scand J Work Environ Health* 2026. [PMID 41582849](https://pubmed.ncbi.nlm.nih.gov/41582849/)
67. Liu N, et al. Mechanism of nano-indium-tin oxide inducing pulmonary alveolar proteinosis in Sprague-Dawley rats. *Zhonghua Lao Dong Wei Sheng Zhi Ye Bing Za Zhi* 2020. [PMID 33287472](https://pubmed.ncbi.nlm.nih.gov/33287472/)
68. Hagmeyer L, Randerath W. Smoking-related interstitial lung disease. *Dtsch Arztebl Int* 2015. [PMID 25797422](https://pubmed.ncbi.nlm.nih.gov/25797422/)
69. Hou SX, et al. Autoimmune pulmonary alveolar proteinosis induced by brigatinib. *Zhonghua Jie He He Hu Xi Za Zhi* 2026. [PMID 42236457](https://pubmed.ncbi.nlm.nih.gov/42236457/)

## 9. Natural history, cohorts and the severity score

The model's presenting patient is an *event* — the crossing of the enrolment
criterion — rather than a fitted equilibrium, and patients whose catabolic
floor holds them below it never present. Inoue's 31.8% asymptomatic fraction
is what that construction is answerable to (diagnostic D16).

70. Seymour JF, Presneill JJ. Pulmonary alveolar proteinosis: progress in the first 44 years. *Am J Respir Crit Care Med* 2002. [PMID 12119235](https://pubmed.ncbi.nlm.nih.gov/12119235/)
71. Inoue Y, et al. Characteristics of a large cohort of patients with autoimmune pulmonary alveolar proteinosis in Japan. *Am J Respir Crit Care Med* 2008. [PMID 18202348](https://pubmed.ncbi.nlm.nih.gov/18202348/) — 223 patients; incidence 0.49/10⁶; 31.8% asymptomatic; the disease severity score; and the finding that DSS does **not** correlate with GMAb titre.
72. Trapnell BC, et al. Pulmonary alveolar proteinosis. *Nat Rev Dis Primers* 2019. [PMID 30846703](https://pubmed.ncbi.nlm.nih.gov/30846703/)
73. Lettieri S, et al. Pathogenesis-driven treatment of primary pulmonary alveolar proteinosis. *Eur Respir Rev* 2024. [PMID 39142709](https://pubmed.ncbi.nlm.nih.gov/39142709/)
74. Rønnov-Jessen I, et al. Pulmonary alveolar proteinosis in Denmark: a retrospective cohort study. *Eur Clin Respir J* 2025. [PMID 41333616](https://pubmed.ncbi.nlm.nih.gov/41333616/)
75. Papiris SA, et al. Pulmonary Alveolar Proteinosis in Greece-Türkiye-Cyprus: Answers in a Real-Life Comparison. *Respirology* 2026. [PMID 42457194](https://pubmed.ncbi.nlm.nih.gov/42457194/)
76. Zhong YF, et al. Clinical characterization of 39 patients with autoimmune pulmonary alveolar proteinosis. *J Int Med Res* 2026. [PMID 41795808](https://pubmed.ncbi.nlm.nih.gov/41795808/)
77. Oh JH, et al. A Multicenter Registry for Rare Interstitial Lung Diseases in Korea: Baseline Characteristics and Clinical Outcomes. *Tuberc Respir Dis (Seoul)* 2026. [PMID 42338237](https://pubmed.ncbi.nlm.nih.gov/42338237/)
78. Xiong J, et al. Outcomes and outcome measures in studies of pulmonary alveolar proteinosis: a scoping review. *Eur Respir Rev* 2026. [PMID 41638878](https://pubmed.ncbi.nlm.nih.gov/41638878/)
79. Bai J, et al. A New Scale to Assess the Severity and Prognosis of Pulmonary Alveolar Proteinosis. *Can Respir J* 2016. [PMID 27635117](https://pubmed.ncbi.nlm.nih.gov/27635117/)
80. Sahoo S, et al. Spontaneous Resolution in Autoimmune Pulmonary Alveolar Proteinosis: A Case Series. *Chest* 2025. [PMID 40506130](https://pubmed.ncbi.nlm.nih.gov/40506130/) — the phenomenon reproduced in scenario S06.
81. Carrington JM, Hershberger DM. Pulmonary Alveolar Proteinosis. *StatPearls* 2026. [PMID 29493933](https://pubmed.ncbi.nlm.nih.gov/29493933/)
82. Huaringa AJ, Malek AO. Pulmonary alveolar proteinosis: a case report and world literature review. *Respirol Case Rep* 2016. [PMID 28031836](https://pubmed.ncbi.nlm.nih.gov/28031836/)

## 10. Biomarkers and imaging

Serum markers in the model are driven by alveolar-to-blood leak proportional
to the burden, which is why they track severity while the antibody titre does
not.

83. Campo I, et al. An exploratory study investigating biomarkers associated with autoimmune pulmonary alveolar proteinosis. *Sci Rep* 2022. [PMID 35610268](https://pubmed.ncbi.nlm.nih.gov/35610268/)
84. Bai JW, et al. CYFRA21-1 is a more sensitive biomarker to assess the severity of pulmonary alveolar proteinosis. *BMC Pulm Med* 2022. [PMID 34980056](https://pubmed.ncbi.nlm.nih.gov/34980056/)
85. Kimura T, et al. Extracellular DNA in bronchoalveolar lavage fluid as a candidate biomarker of disease severity in autoimmune pulmonary alveolar proteinosis. *Respir Investig* 2026. [PMID 41619656](https://pubmed.ncbi.nlm.nih.gov/41619656/)
86. Huang H, et al. Geospatial and temporal trends of interstitial lung disease subtypes and KL-6 biomarker in China. *Respir Res* 2026. [PMID 41699651](https://pubmed.ncbi.nlm.nih.gov/41699651/)
87. Tokura S, et al. A Semiquantitative Computed Tomographic Grading System for Evaluating Therapeutic Response in Pulmonary Alveolar Proteinosis. *Ann Am Thorac Soc* 2017. [PMID 28489417](https://pubmed.ncbi.nlm.nih.gov/28489417/) — CT grading; the model outputs mean lung density in HU because PAGE used it as an endpoint.
88. Chang Q, et al. Discriminating nodular pulmonary alveolar proteinosis from lung adenocarcinoma with radiomic model. *Ther Adv Respir Dis* 2026. [PMID 41905780](https://pubmed.ncbi.nlm.nih.gov/41905780/)

## 11. Whole lung lavage

The model treats lavage as an instantaneous mechanical removal that does not
touch the accumulation rate — so time to recurrence is mass removed divided by
net imbalance, and is a *prediction* rather than a parameter (D12).

89. Beccaria M, et al. Long-term durable benefit after whole lung lavage in pulmonary alveolar proteinosis. *Eur Respir J* 2004. [PMID 15083749](https://pubmed.ncbi.nlm.nih.gov/15083749/) — >70% recurrence-free at 7 years; DLCO still 75 ± 19 %pred and A-aDO2 27 ± 11 mmHg at 5 years. The residual-deficit ceiling in the model answers to this.
90. Campo I, et al. Whole lung lavage therapy for pulmonary alveolar proteinosis: a global survey of current practices and procedures. *Orphanet J Rare Dis* 2016. [PMID 27577926](https://pubmed.ncbi.nlm.nih.gov/27577926/)
91. Ataya A, et al. Whole Lung Lavage in Pulmonary Alveolar Proteinosis. *Chest* 2026. [PMID 41072904](https://pubmed.ncbi.nlm.nih.gov/41072904/)
92. Luo W, et al. Gravity-driven modified whole-lung lavage improves procedural efficiency in pulmonary alveolar proteinosis. *J Thorac Dis* 2026. [PMID 42444971](https://pubmed.ncbi.nlm.nih.gov/42444971/)
93. Ju R, et al. High-frequency jet ventilation in managing airway during whole-lung lavage under general anesthesia. *Respir Med Case Rep* 2025. [PMID 40290789](https://pubmed.ncbi.nlm.nih.gov/40290789/)
94. Pu J, et al. Timing of Whole Lung Lavage in Autoimmune Pulmonary Alveolar Proteinosis with Concurrent Opportunistic Infection. *Infect Drug Resist* 2026. [PMID 41852848](https://pubmed.ncbi.nlm.nih.gov/41852848/)
95. Le H, et al. Successful Whole-Lung Lavage in Anti-GM-CSF-Negative Pulmonary Alveolar Proteinosis. *Am J Case Rep* 2026. [PMID 42130041](https://pubmed.ncbi.nlm.nih.gov/42130041/)
96. Huang L, et al. Nursing care for patients with pulmonary alveolar proteinosis undergoing whole-lung lavage therapy. *Medicine (Baltimore)* 2026. [PMID 42499099](https://pubmed.ncbi.nlm.nih.gov/42499099/)
97. Iwase A, et al. A case of recurrent pulmonary alveolar proteinosis treated by pulmonary lavage. *Nihon Kyobu Shikkan Gakkai Zasshi* 1992. [PMID 1630060](https://pubmed.ncbi.nlm.nih.gov/1630060/)
98. Tsuchiyama T, et al. Pulmonary alveolar proteinosis with coincidental thoracic injury: successful bronchoalveolar lavage. *Nihon Kyobu Shikkan Gakkai Zasshi* 1995. [PMID 7609340](https://pubmed.ncbi.nlm.nih.gov/7609340/)
99. Schoch OD, et al. BAL findings in a patient with pulmonary alveolar proteinosis successfully treated with GM-CSF. *Thorax* 2002. [PMID 11867836](https://pubmed.ncbi.nlm.nih.gov/11867836/)
100. Herron M, et al. Optimising bronchoalveolar lavage: lessons from alpha-1 antitrypsin deficiency. *Thorax* 2024. [PMID 39586664](https://pubmed.ncbi.nlm.nih.gov/39586664/)
101. Bäckström E, et al. Possible Extraction of Drugs from Lung Tissue During Broncho-alveolar Lavage. *J Pharm Sci* 2022. [PMID 34890629](https://pubmed.ncbi.nlm.nih.gov/34890629/) — BAL dilution, which is why ELF and BAL concentrations differ by roughly two orders of magnitude in the model.

## 12. GM-CSF therapy: the trials the model is validated against

102. Kavuru MS, et al. Exogenous granulocyte-macrophage colony-stimulating factor administration for pulmonary alveolar proteinosis. *Am J Respir Crit Care Med* 2000. [PMID 10764303](https://pubmed.ncbi.nlm.nih.gov/10764303/)
103. Seymour JF, et al. Therapeutic efficacy of granulocyte-macrophage colony-stimulating factor in patients with idiopathic acquired alveolar proteinosis. *Am J Respir Crit Care Med* 2001. [PMID 11179134](https://pubmed.ncbi.nlm.nih.gov/11179134/) — subcutaneous GM-CSF, partial responses at escalating dose; scenarios S12–S13.
104. Tazawa R, et al. Inhaled granulocyte/macrophage-colony stimulating factor as therapy for pulmonary alveolar proteinosis. *Am J Respir Crit Care Med* 2010. [PMID 20167854](https://pubmed.ncbi.nlm.nih.gov/20167854/) — phase 2: A-aDO2 −12.3 mmHg, 62% response, and 29 of 35 stable for a year *after stopping* (the observation that argues against the flat-recovery reading in D08).
105. Tazawa R, et al. Inhaled GM-CSF for Pulmonary Alveolar Proteinosis. *N Engl J Med* 2019. [PMID 31483963](https://pubmed.ncbi.nlm.nih.gov/31483963/) — PAGE: A-aDO2 −4.50 vs +0.17 mmHg; CT density −36.08 HU. The placebo arm that did *not* move, in the trial that excluded improvers during a 12-week run-in (diagnostic D10).
106. Trapnell BC, et al. Inhaled Molgramostim Therapy in Autoimmune Pulmonary Alveolar Proteinosis. *N Engl J Med* 2020. [PMID 32897035](https://pubmed.ncbi.nlm.nih.gov/32897035/) — IMPALA: A-aDO2 −12.8 vs −6.6 mmHg at 24 weeks, n=138; continuous better than intermittent.
107. Trapnell BC, et al. Phase 3 Trial of Inhaled Molgramostim in Autoimmune Pulmonary Alveolar Proteinosis. *N Engl J Med* 2025. [PMID 40834301](https://pubmed.ncbi.nlm.nih.gov/40834301/) — IMPALA-2: DLCO +9.8 vs +3.8 %pred at 24 weeks, +11.6 vs +4.7 at 48 weeks, SGRQ-T −11.5 vs −4.9. The primary validation target (D07) and the source of the 48-week failure (D08).
108. Trapnell BC, et al. Inhaled molgramostim therapy for the treatment of autoimmune pulmonary alveolar proteinosis: a plain language summary. *Hosp Pract* 2024. [PMID 39165153](https://pubmed.ncbi.nlm.nih.gov/39165153/)
109. Tazawa R, et al. Inhaled granulocyte-macrophage colony-stimulating factor for autoimmune pulmonary alveolar proteinosis: from patients to practice. *Respir Investig* 2026. [PMID 42462405](https://pubmed.ncbi.nlm.nih.gov/42462405/)
110. Dang M, et al. Efficacy and Safety of Inhaled GM-CSF in Autoimmune Pulmonary Alveolar Proteinosis: A Systematic Review and Meta-analysis. *Lung* 2026. [PMID 42045620](https://pubmed.ncbi.nlm.nih.gov/42045620/)
111. Naeem U, et al. Efficacy of recombinant human GM-CSF compared to placebo for autoimmune pulmonary alveolar proteinosis: a meta-analysis. *BMC Pulm Med* 2026. [PMID 41742124](https://pubmed.ncbi.nlm.nih.gov/41742124/)
112. Higgins M, et al. Targeting autoimmune pulmonary alveolar proteinosis with GM-CSF: insights from clinical trials and emerging strategies. *Expert Opin Biol Ther* 2026. [PMID 42466621](https://pubmed.ncbi.nlm.nih.gov/42466621/)
113. Matsumoto T, et al. A Marked Improvement in Mild but Progressive Autoimmune Pulmonary Alveolar Proteinosis Treated with Inhaled GM-CSF. *Intern Med* 2026. [PMID 42366042](https://pubmed.ncbi.nlm.nih.gov/42366042/)
114. McCarthy J, et al. Recombinant GM-CSF drug evaluation review. *Immunotherapy* 2025. [PMID 41085041](https://pubmed.ncbi.nlm.nih.gov/41085041/)
115. Takazoe M, et al. Sargramostim in patients with Crohn's disease: results of a phase 1-2 study. *J Gastroenterol* 2009. [PMID 19352588](https://pubmed.ncbi.nlm.nih.gov/19352588/) — systemic sargramostim dosing and tolerability.
116. Kelsen JR, et al. Phase I trial of sargramostim in pediatric Crohn's disease. *Inflamm Bowel Dis* 2010. [PMID 20052780](https://pubmed.ncbi.nlm.nih.gov/20052780/)

## 13. Guidelines and current management

117. McCarthy C, et al. European Respiratory Society guidelines for the diagnosis and management of pulmonary alveolar proteinosis. *Eur Respir J* 2024. [PMID 39147411](https://pubmed.ncbi.nlm.nih.gov/39147411/)
118. Alfaro T, et al. Summary for clinicians: ERS guidelines on pulmonary alveolar proteinosis. *Breathe (Sheff)* 2025. [PMID 40365091](https://pubmed.ncbi.nlm.nih.gov/40365091/)
119. Huang JF, et al. Chinese translation of European Respiratory Society guidelines for pulmonary alveolar proteinosis. *Zhonghua Jie He He Hu Xi Za Zhi* 2025. [PMID 41073306](https://pubmed.ncbi.nlm.nih.gov/41073306/)
120. Jouneau S, et al. Pharmacotherapy for Autoimmune Pulmonary Alveolar Proteinosis. *Drugs* 2025. [PMID 40866780](https://pubmed.ncbi.nlm.nih.gov/40866780/)
121. Bendstrup E, et al. Recent advances in the diagnosis and management of pulmonary alveolar proteinosis. *Expert Rev Respir Med* 2025. [PMID 40772394](https://pubmed.ncbi.nlm.nih.gov/40772394/)

## 14. Antibody-directed therapy (the model's untested predictions)

Diagnostic D19 flags all of this as hypothesis. Rituximab cannot reach the
CD20-negative long-lived plasma cell that supplies most of the antibody;
plasmapheresis removes only the intravascular pool; FcRn blockade lowers total
IgG enough to cross the critical titre from a median baseline.

122. The rituximab evidence in aPAP is a small number of open-label series rather than a controlled trial, and is summarised in refs. 117, 120 and 121 (ERS guidelines 2024; Jouneau 2025; Bendstrup 2025). A specific citation was cut from this list rather than guessed: the PMID first written here from memory turned out to belong to an unrelated paper, which is the second time that happened while compiling this file.
123. Springer JM, et al. Dose-dependent Pharmacological Response to Rituximab in ANCA-associated vasculitis. *J Rheumatol* 2021. [PMID 34334366](https://pubmed.ncbi.nlm.nih.gov/34334366/) — rituximab exposure–B-cell response, used for the depletion module.
124. Liang P, et al. Targeting FcRn for immunomodulation: a promising therapy in autoimmune inflammatory rheumatic diseases. *Inflamm Res* 2026. [PMID 42257855](https://pubmed.ncbi.nlm.nih.gov/42257855/)
125. Khateb M, et al. Nipocalimab and other FcRn blockers in neuromuscular disorders. *Pharmacol Ther* 2026. [PMID 42398839](https://pubmed.ncbi.nlm.nih.gov/42398839/)
126. Jiang L, et al. Efgartigimod for generalized myasthenia gravis: a comprehensive review of clinical evidence. *Front Neurol* 2026. [PMID 42358938](https://pubmed.ncbi.nlm.nih.gov/42358938/) — the magnitude of IgG reduction achievable, which is what scenarios S24–S25 turn on.
127. Sirina J, et al. The uptake, intracellular trafficking and recycling of FcRn-blocking therapeutics in human endothelial cells. *Cell Immunol* 2026. [PMID 42385333](https://pubmed.ncbi.nlm.nih.gov/42385333/)
128. Khumalo A, et al. Systemic Inflammation as a Modulator of FcRn-dependent IgG Pharmacokinetics. *Curr HIV/AIDS Rep* 2026. [PMID 42430031](https://pubmed.ncbi.nlm.nih.gov/42430031/)
129. Komatsu Y, et al. Anti-LGI1 encephalitis during FcRn inhibition with efgartigimod for myasthenia gravis. *Front Immunol* 2026. [PMID 42233017](https://pubmed.ncbi.nlm.nih.gov/42233017/)

## 15. Inhaled delivery, and the reach problem

The model's central pharmacological claim is that aerosol follows ventilation
while the burden sits in units that do not ventilate — so *reach*, not dose,
limits inhaled GM-CSF (diagnostic D17).

130. Hertel S, et al. Prediction of protein degradation during vibrating mesh nebulization. *Eur J Pharm Biopharm* 2014. [PMID 24709473](https://pubmed.ncbi.nlm.nih.gov/24709473/)
131. Mahri S, et al. Nebulization of PEGylated recombinant human deoxyribonuclease I using vibrating membrane nebulizers. *Eur J Pharm Sci* 2023. [PMID 37423579](https://pubmed.ncbi.nlm.nih.gov/37423579/)
132. Li C, et al. Gelatin Stabilizes Nebulized Proteins in Pulmonary Drug Delivery. *ACS Biomater Sci Eng* 2022. [PMID 35608934](https://pubmed.ncbi.nlm.nih.gov/35608934/)
133. Boger E, Fridén M. Assessment of Epithelial Lining Fluid Partitioning of Systemically Administered Monoclonal Antibodies. *J Pharm Sci* 2023. [PMID 36632919](https://pubmed.ncbi.nlm.nih.gov/36632919/) — ELF/plasma partitioning of IgG, the source of the transudation ratio κ.
134. Out TA, et al. Permeability or local production of immunoglobulins and other inflammatory proteins in asthma. *Eur Respir J Suppl* 1991. [PMID 1953911](https://pubmed.ncbi.nlm.nih.gov/1953911/)
135. Florio G, et al. Impact of Positive End-Expiratory Pressure and FiO₂ on Lung Mechanics and Intrapulmonary Shunt. *J Intensive Care Med* 2024. [PMID 37926984](https://pubmed.ncbi.nlm.nih.gov/37926984/)

## 16. Gas exchange: the physics the endpoints are computed from

A-aDO2 and DLCO are not scores in this model. They come from a shunt equation,
the alveolar gas equation, the Severinghaus dissociation curve (and its
closed-form inverse) and the Fick relation.

136. Severinghaus JW. Simple, accurate equations for human blood O₂ dissociation computations. *J Appl Physiol* 1979. [PMID 35496](https://pubmed.ncbi.nlm.nih.gov/35496/) — the curve and inverse implemented in `$GLOBAL`.
137. Roughton FJ, Forster RE. Relative importance of diffusion and chemical reaction rates in determining rate of exchange of gases in the human lung. *J Appl Physiol* 1957. [PMID 13475180](https://pubmed.ncbi.nlm.nih.gov/13475180/) — the membrane/capillary partition of DLCO.
138. D'Souza AW, et al. A comparison of pulmonary capillary blood volume and membrane diffusing capacity assessed via two pulmonary function techniques. *ERJ Open Res* 2026. [PMID 41878275](https://pubmed.ncbi.nlm.nih.gov/41878275/)
139. Collins SÉ, et al. Relationship of Pulmonary Vascular Structure and Function With Exercise Capacity in Health and COPD. *Chest* 2025. [PMID 39368737](https://pubmed.ncbi.nlm.nih.gov/39368737/)
140. Wells AU, et al. Functional impairment in lone cryptogenic fibrosing alveolitis and fibrosing alveolitis associated with systemic sclerosis. *Am J Respir Crit Care Med* 1997. [PMID 9154872](https://pubmed.ncbi.nlm.nih.gov/9154872/)
141. Satrell E, et al. Development of lung diffusion to adulthood following extremely preterm birth. *Eur Respir J* 2022. [PMID 34625479](https://pubmed.ncbi.nlm.nih.gov/34625479/)

## 17. Host defence, infection and the complications of a lipid-rich alveolus

142. Rosen LB, et al. Nocardia-induced granulocyte macrophage colony-stimulating factor is neutralized by autoantibodies in disseminated/extrapulmonary nocardiosis. *Clin Infect Dis* 2015. [PMID 25472947](https://pubmed.ncbi.nlm.nih.gov/25472947/)
143. Lo YF, et al. The Pathogenic Role of Anti-Granulocyte-Macrophage Colony-Stimulating Factor Autoantibodies in Nocardiosis. *J Clin Immunol* 2024. [PMID 39133333](https://pubmed.ncbi.nlm.nih.gov/39133333/)
144. Menon V, et al. Concurrent Nocardia, Cryptococcus and Mycobacterium Infections Unmask Anti-GM-CSF Antibody-Associated Immunodeficiency. *Case Rep Infect Dis* 2026. [PMID 42358565](https://pubmed.ncbi.nlm.nih.gov/42358565/)
145. Coirier V, et al. A case report of Covid-19 in an autoimmune pulmonary alveolar proteinosis. *Respir Med Case Rep* 2023. [PMID 36874265](https://pubmed.ncbi.nlm.nih.gov/36874265/)
146. El Mawla Z, et al. A Rare Case of Pulmonary Alveolar Proteinosis Superimposed by Severe COVID-19 Pneumonia. *Respirol Case Rep* 2026. [PMID 41551542](https://pubmed.ncbi.nlm.nih.gov/41551542/)
147. Costiniuk CT, et al. Potential role of alveolar macrophages in HIV persistence and lung disease. *Virologie* 2024. [PMID 39248668](https://pubmed.ncbi.nlm.nih.gov/39248668/)

## 18. Irreversible remodelling and transplantation

The ceiling on recovery. Beccaria's 5-year post-lavage DLCO of 75 ± 19 %pred is
the observation the fibrosis/structural term answers to.

148. Guirriec Y, et al. Pulmonary fibrosis in patients with autoimmune pulmonary alveolar proteinosis: a retrospective nationwide study. *ERJ Open Res* 2024. [PMID 39624377](https://pubmed.ncbi.nlm.nih.gov/39624377/)
149. Sugino K, et al. Pleuroparenchymal Fibroelastosis Complicated by Pulmonary Alveolar Proteinosis After Peripheral Blood Stem Cell Transplantation. *Respirol Case Rep* 2026. [PMID 41815743](https://pubmed.ncbi.nlm.nih.gov/41815743/)
150. Schwarz S, et al. Lung Transplantation for Pulmonary Alveolar Proteinosis — a Retrospective International Multi-Center Analysis. *Eur Respir J* 2026. [PMID 42425729](https://pubmed.ncbi.nlm.nih.gov/42425729/)
151. Takaki M, et al. Recurrence of pulmonary alveolar proteinosis after bilateral lung transplantation. *Respir Med Case Rep* 2016. [PMID 27595063](https://pubmed.ncbi.nlm.nih.gov/27595063/) — recurrence in the graft: the autoantibody is systemic, so replacing the lung does not remove the cause.
152. Lopez O, et al. Diagnostic and Therapeutic Dilemmas in Recurrent Pulmonary Alveolar Proteinosis After Bilateral Lung Transplantation. *Respirol Case Rep* 2026. [PMID 42293034](https://pubmed.ncbi.nlm.nih.gov/42293034/)
153. Sharma S, et al. Atrial septal defect and pulmonary alveolar proteinosis in an adult. *Indian J Thorac Cardiovasc Surg* 2026. [PMID 41835822](https://pubmed.ncbi.nlm.nih.gov/41835822/)
154. Zheng Q, et al. Rare case of autoimmune pulmonary alveolar proteinosis with ANCA positivity. *BMC Pulm Med* 2026. [PMID 41735921](https://pubmed.ncbi.nlm.nih.gov/41735921/)
155. Ali MB, et al. A Case of Autoimmune Pulmonary Alveolar Proteinosis. *Clin Case Rep* 2026. [PMID 41716453](https://pubmed.ncbi.nlm.nih.gov/41716453/)

## 19. QSP methodology

156. Rao R, et al. A quantitative systems pharmacology model of the pathophysiology and treatment of COVID-19. *NPJ Syst Biol Appl* 2023. [PMID 37059734](https://pubmed.ncbi.nlm.nih.gov/37059734/)
157. Braniff N, et al. An integrated quantitative systems pharmacology virtual population approach for calibration. *CPT Pharmacometrics Syst Pharmacol* 2025. [PMID 39508122](https://pubmed.ncbi.nlm.nih.gov/39508122/) — virtual-population construction, the approach used for the trial replications.
158. Lu T, et al. gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve. *CPT Pharmacometrics Syst Pharmacol* 2024. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)

---

## Note on two citations used for numbers this model depends on

**The 5 µg/mL critical threshold** (ref. 5, Sakagami 2010) is quoted here as a
*serum* concentration that is "similar in lungs and blood". The model does not
assume that: it carries an explicit ELF/serum transudation ratio and a
neutralising fraction, and diagnostic D04 shows the half-signal titre landing
at 5.6 µg/mL. Those two parameters are the calibration levers. The *shape* of
the curve — steep, then flat — is not a lever; it is what a stoichiometric
buffer does.

**The alveolar surfactant pool and its turnover** (refs. 31–34) are the
weakest quantitative link in the model. Human alveolar phospholipid pool size
and de novo synthesis rate are known to roughly an order of magnitude, mostly
from animal and neonatal stable-isotope work. The model fixes the healthy pool
at 300 mg with 200 mg/d of net turnover and *derives* every Vmax from that, so
the flux split is exact by construction — but the absolute scale carries that
uncertainty, and any statement in milligrams should be read as an order of
magnitude. What the model is actually sensitive to is the *ratio* of the
macrophage sink to production, not the absolute size of either.
