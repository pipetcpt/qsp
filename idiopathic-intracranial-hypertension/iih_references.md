# Idiopathic Intracranial Hypertension (IIH) — QSP 모델 참고문헌

> **Verification note.** Every PMID below was resolved and its title,
> journal and year read back from the PubMed E-utilities API at build
> time by `scripts`-style title-field queries, not recalled from memory.
> Entries whose returned title did not overlap the requested title were
> discarded rather than guessed at. 112 references.

The sections follow the structure of the model: the first five build
the fixed-point equation, the next four supply the drives that enter it,
and the last five supply the endpoints the equation is read out through.


## 1. 정의 · 진단기준 · 역학 (Definition, diagnostic criteria, epidemiology)

1. Whiteley W, Al-Shahi R, Warlow CP et al. **CSF opening pressure: reference interval and the effect of body mass index.** *Neurology* 2006. [PMID 17101909](https://pubmed.ncbi.nlm.nih.gov/17101909/)
2. Raoof N, Sharrack B, Pepper IM et al. **The incidence and prevalence of idiopathic intracranial hypertension in Sheffield, UK.** *European journal of neurology* 2011. [PMID 21426442](https://pubmed.ncbi.nlm.nih.gov/21426442/)
3. Friedman DI, Liu GT, Digre KB. **Revised diagnostic criteria for the pseudotumor cerebri syndrome in adults and children.** *Neurology* 2013. [PMID 23966248](https://pubmed.ncbi.nlm.nih.gov/23966248/)
4. McGeeney BE, Friedman DI. **Pseudotumor cerebri pathophysiology.** *Headache* 2014. [PMID 24433163](https://pubmed.ncbi.nlm.nih.gov/24433163/)
5. Markey KA, Mollan SP, Jensen RH et al. **Understanding idiopathic intracranial hypertension: mechanisms, management, and future directions.** *The Lancet. Neurology* 2016. [PMID 26700907](https://pubmed.ncbi.nlm.nih.gov/26700907/)
6. Sundholm A, Burkill S, Sveinsson O et al. **Population-based incidence and clinical characteristics of idiopathic intracranial hypertension.** *Acta neurologica Scandinavica* 2017. [PMID 28244170](https://pubmed.ncbi.nlm.nih.gov/28244170/)
7. Mollan SP, Davies B, Silver NC et al. **Idiopathic intracranial hypertension: consensus guidelines on management.** *Journal of neurology, neurosurgery, and psychiatry* 2018. [PMID 29903905](https://pubmed.ncbi.nlm.nih.gov/29903905/)
8. Mollan SP, Hornby C, Mitchell J et al. **Evaluation and management of adult idiopathic intracranial hypertension.** *Practical neurology* 2018. [PMID 30154235](https://pubmed.ncbi.nlm.nih.gov/30154235/)
9. Adderley NJ, Subramanian A, Nirantharakumar K et al. **Association Between Idiopathic Intracranial Hypertension and Risk of Cardiovascular Diseases in Women in the United Kingdom.** *JAMA neurology* 2019. [PMID 31282950](https://pubmed.ncbi.nlm.nih.gov/31282950/)
10. Bsteh G, Pemp B, Marik W et al. **Diagnosis, treatment and monitoring of idiopathic intracranial hypertension: Consensus recommendations of the Austrian IIH network (AN4IH).** *Cephalalgia : an international journal of headache* 2025. [PMID 41223050](https://pubmed.ncbi.nlm.nih.gov/41223050/)


## 2. 뇌척수액 동역학과 Davson 관계식 (CSF hydrodynamics and the Davson relation)

11. Cutler RW, Page L, Galicich J et al. **Formation and absorption of cerebrospinal fluid in man.** *Brain : a journal of neurology* 1968. [PMID 5304069](https://pubmed.ncbi.nlm.nih.gov/5304069/)
12. Davson H, Hollingsworth G, Segal MB. **The mechanism of drainage of the cerebrospinal fluid.** *Brain : a journal of neurology* 1970. [PMID 5490270](https://pubmed.ncbi.nlm.nih.gov/5490270/)
13. Marmarou A, Shulman K, LaMorgese J. **Compartmental analysis of compliance and outflow resistance of the cerebrospinal fluid system.** *Journal of neurosurgery* 1975. [PMID 1181384](https://pubmed.ncbi.nlm.nih.gov/1181384/)
14. Marmarou A, Shulman K, Rosende RM. **A nonlinear analysis of the cerebrospinal fluid system and intracranial pressure dynamics.** *Journal of neurosurgery* 1978. [PMID 632857](https://pubmed.ncbi.nlm.nih.gov/632857/)
15. Huang TY, Chung HW, Chen MY et al. **Supratentorial cerebrospinal fluid production rate in healthy adults: quantification with two-dimensional cine phase-contrast MR imaging with high temporal and spatial resolution.** *Radiology* 2004. [PMID 15516623](https://pubmed.ncbi.nlm.nih.gov/15516623/)
16. Takahashi H, Tanaka H, Fujita N et al. **Variation in supratentorial cerebrospinal fluid production rate in one day: measurement by nontriggered phase-contrast magnetic resonance imaging.** *Japanese journal of radiology* 2011. [PMID 21359936](https://pubmed.ncbi.nlm.nih.gov/21359936/)
17. Bateman GA, Siddique SH. **Cerebrospinal fluid absorption block at the vertex in chronic hydrocephalus: obstructed arachnoid granulations or elevated venous pressure?.** *Fluids and barriers of the CNS* 2014. [PMID 24955236](https://pubmed.ncbi.nlm.nih.gov/24955236/)
18. Leinonen V, Vanninen R, Rauramaa T. **Cerebrospinal fluid circulation and hydrocephalus.** *Handbook of clinical neurology* 2017. [PMID 28987185](https://pubmed.ncbi.nlm.nih.gov/28987185/)
19. Ahn JH, Cho H, Kim JH et al. **Meningeal lymphatic vessels at the skull base drain cerebrospinal fluid.** *Nature* 2019. [PMID 31341278](https://pubmed.ncbi.nlm.nih.gov/31341278/)
20. Du T, Raghunandan A, Mestre H et al. **Restoration of cervical lymphatic vessel function in aging rescues cerebrospinal fluid drainage.** *Nature aging* 2024. [PMID 39147980](https://pubmed.ncbi.nlm.nih.gov/39147980/)


## 3. 맥락총 분비 기전과 약물 표적 (Choroid plexus secretion and its drug targets)

21. McCarthy KD, Reed DJ. **The effect of acetazolamide and furosemide on cerebrospinal fluid production and choroid plexus carbonic anhydrase activity.** *The Journal of pharmacology and experimental therapeutics* 1974. [PMID 4207244](https://pubmed.ncbi.nlm.nih.gov/4207244/)
22. Vogh BP, Langham MR Jr. **The effect of furosemide and bumetanide on cerebrospinal fluid formation.** *Brain research* 1981. [PMID 6791768](https://pubmed.ncbi.nlm.nih.gov/6791768/)
23. Kumpulainen T, Korhonen LK. **Immunohistochemical localization of carbonic anhydrase isoenzyme C in the central and peripheral nervous system of the mouse.** *The journal of histochemistry and cytochemistry : official journal of the Histochemistry Society* 1982. [PMID 6801110](https://pubmed.ncbi.nlm.nih.gov/6801110/)
24. Nógrádi A, Kelly C, Carter ND. **Localization of acetazolamide-resistant carbonic anhydrase III in human and rat choroid plexus by immunocytochemistry and in situ hybridisation.** *Neuroscience letters* 1993. [PMID 8506074](https://pubmed.ncbi.nlm.nih.gov/8506074/)
25. Fisone G, Snyder GL, Aperia A et al. **Na+,K(+)-ATPase phosphorylation in the choroid plexus: synergistic regulation by serotonin/protein kinase C and isoproterenol/cAMP-PK/PP-1 pathways.** *Molecular medicine (Cambridge, Mass.)* 1998. [PMID 9606178](https://pubmed.ncbi.nlm.nih.gov/9606178/)
26. Oshio K, Watanabe H, Song Y et al. **Reduced cerebrospinal fluid production and intracranial pressure in mice lacking choroid plexus water channel Aquaporin-1.** *FASEB journal : official publication of the Federation of American Societies for Experimental Biology* 2005. [PMID 15533949](https://pubmed.ncbi.nlm.nih.gov/15533949/)
27. Damkier HH, Brown PD, Praetorius J. **Cerebrospinal fluid secretion by the choroid plexus.** *Physiological reviews* 2013. [PMID 24137023](https://pubmed.ncbi.nlm.nih.gov/24137023/)
28. Lehtinen MK, Bjornsson CS, Dymecki SM et al. **The choroid plexus and cerebrospinal fluid: emerging roles in development, disease, and therapy.** *The Journal of neuroscience : the official journal of the Society for Neuroscience* 2013. [PMID 24198345](https://pubmed.ncbi.nlm.nih.gov/24198345/)
29. Uldall M, Botfield H, Jansen-Olesen I et al. **Acetazolamide lowers intracranial pressure and modulates the cerebrospinal fluid secretion pathway in healthy rats.** *Neuroscience letters* 2017. [PMID 28219789](https://pubmed.ncbi.nlm.nih.gov/28219789/)
30. Sadegh C, Xu H, Sutin J et al. **Choroid plexus-targeted NKCC1 overexpression to treat post-hemorrhagic hydrocephalus.** *Neuron* 2023. [PMID 36893755](https://pubmed.ncbi.nlm.nih.gov/36893755/)
31. Wang Q, Liu F, Li Y et al. **Choroid plexus CCL2‒CCR2 signaling orchestrates macrophage recruitment and cerebrospinal fluid hypersecretion in hydrocephalus.** *Acta pharmaceutica Sinica. B* 2024. [PMID 39525574](https://pubmed.ncbi.nlm.nih.gov/39525574/)
32. Jensen MN, Israelsen IME, Wardman JH et al. **Glucagon-like peptide-1 receptor modulates cerebrospinal fluid secretion and intracranial pressure in rats.** *Fluids and barriers of the CNS* 2025. [PMID 40275284](https://pubmed.ncbi.nlm.nih.gov/40275284/)


## 4. 정맥동 협착 — 증폭기 (Venous sinus stenosis: the Starling resistor that closes the loop)

33. Kotani J, Sugioka S, Momota Y et al. **Effect of sevoflurane on intracranial pressure, sagittal sinus pressure, and the intracranial volume-pressure relation in cats.** *Journal of neurosurgical anesthesiology* 1992. [PMID 15815463](https://pubmed.ncbi.nlm.nih.gov/15815463/)
34. McGonigal A, Bone I, Teasdale E. **Resolution of transverse sinus stenosis in idiopathic intracranial hypertension after L-P shunt.** *Neurology* 2004. [PMID 14872049](https://pubmed.ncbi.nlm.nih.gov/14872049/)
35. Higgins JN, Gillard JH, Owler BK et al. **MR venography in idiopathic intracranial hypertension: unappreciated and misunderstood.** *Journal of neurology, neurosurgery, and psychiatry* 2004. [PMID 15026510](https://pubmed.ncbi.nlm.nih.gov/15026510/)
36. De Simone R, Ranieri A, Cardillo G et al. **High prevalence of bilateral transverse sinus stenosis-associated IIHWOP in unresponsive chronic headache sufferers: pathogenetic implications in primary headache progression.** *Cephalalgia : an international journal of headache* 2011. [PMID 21493643](https://pubmed.ncbi.nlm.nih.gov/21493643/)
37. DE Simone R, Ranieri A, Bonavita V. **Starling resistors, autoregulation of cerebral perfusion and the pathogenesis of idiopathic intracranial hypertension.** *Panminerva medica* 2017. [PMID 27598891](https://pubmed.ncbi.nlm.nih.gov/27598891/)
38. Giridharan N, Patel SK, Ojugbeli A et al. **Understanding the complex pathophysiology of idiopathic intracranial hypertension and the evolving role of venous sinus stenting: a comprehensive review of the literature.** *Neurosurgical focus* 2018. [PMID 29961379](https://pubmed.ncbi.nlm.nih.gov/29961379/)
39. Buell T, Ding D, Raper D et al. **Resolution of venous pressure gradient in a patient with idiopathic intracranial hypertension after ventriculoperitoneal shunt placement: A proof of secondary cerebral sinovenous stenosis.** *Surgical neurology international* 2021. [PMID 33500829](https://pubmed.ncbi.nlm.nih.gov/33500829/)
40. Fargen KM, Midtlien JP, Margraf CR et al. **Idiopathic intracranial hypertension pathogenesis: The jugular hypothesis.** *Interventional neuroradiology : journal of peritherapeutic neuroradiology, surgical procedures and related neurosciences* 2024. [PMID 39113487](https://pubmed.ncbi.nlm.nih.gov/39113487/)
41. Abou-Mrad T, Alaraj A. **Venous sinus stenosis intracranial hypertension, rethinking idiopathic intracranial hypertension in the setting of venous sinus stenosis: A call for new nomenclature and diagnostic precision.** *Interventional neuroradiology : journal of peritherapeutic neuroradiology, surgical procedures and related neurosciences* 2024. [PMID 39311024](https://pubmed.ncbi.nlm.nih.gov/39311024/)
42. Cagnazzo F, Villain M, van Dokkum LE et al. **Concordance between venous sinus pressure and intracranial pressure in patients investigated for idiopathic intracranial hypertension.** *The journal of headache and pain* 2024. [PMID 39289632](https://pubmed.ncbi.nlm.nih.gov/39289632/)
43. White T, Shah K, Ryu B et al. **Diagnostic accuracy of venous manometry to predict elevated intracranial pressure.** *Frontiers in neurology* 2026. [PMID 41743050](https://pubmed.ncbi.nlm.nih.gov/41743050/)


## 5. 두개내 순응도 · ICP 계측 (Intracranial compliance and the measurement of ICP)

44. Poca MA, Sahuquillo J, Topczewski T et al. **Posture-induced changes in intracranial pressure: a comparative study in patients with and without a cerebrospinal fluid block at the craniovertebral junction.** *Neurosurgery* 2006. [PMID 16639324](https://pubmed.ncbi.nlm.nih.gov/16639324/)
45. Velazquez Sanchez VF, Al Dayri G, Tschan CA. **Long-term telemetric intracranial pressure monitoring for diagnosis and therapy optimisation of idiopathic intracranial hypertension.** *BMC neurology* 2021. [PMID 34493231](https://pubmed.ncbi.nlm.nih.gov/34493231/)
46. Mitchell JL, Buckham R, Lyons H et al. **Evaluation of diurnal and postural intracranial pressure employing telemetric monitoring in idiopathic intracranial hypertension.** *Fluids and barriers of the CNS* 2022. [PMID 36320018](https://pubmed.ncbi.nlm.nih.gov/36320018/)
47. Maroufi SF, Katkade O, Puppalla P et al. **Lumbar Puncture Opening Pressure and Continuous Intracranial Pressure Monitoring: Concordance and Clinical Implications in Idiopathic Intracranial Hypertension.** *Neurosurgery* 2026. [PMID 42084378](https://pubmed.ncbi.nlm.nih.gov/42084378/)


## 6. 비만 · 복강내압 · 중심정맥압 (Obesity, intra-abdominal pressure, central venous pressure)

48. Sugerman HJ, DeMaria EJ, Felton WL 3rd et al. **Increased intra-abdominal pressure and cardiac filling pressures in obesity-associated pseudotumor cerebri.** *Neurology* 1997. [PMID 9270586](https://pubmed.ncbi.nlm.nih.gov/9270586/)
49. Ko MW, Chang SC, Ridha MA et al. **Weight gain and recurrence in idiopathic intracranial hypertension: a case-control study.** *Neurology* 2011. [PMID 21536635](https://pubmed.ncbi.nlm.nih.gov/21536635/)
50. Schwartz R, Kliper E, Stern N et al. **The obesity pattern of idiopathic intracranial hypertension in men.** *Graefe's archive for clinical and experimental ophthalmology = Albrecht von Graefes Archiv fur klinische und experimentelle Ophthalmologie* 2013. [PMID 23955783](https://pubmed.ncbi.nlm.nih.gov/23955783/)
51. Karam M, Alsaif A, Alroumi D et al. **Obstructive Sleep Apnea in Idiopathic Intracranial Hypertension: Systematic Review and Meta-Analysis.** *Neuro-ophthalmology (Aeolus Press)* 2025. [PMID 40919091](https://pubmed.ncbi.nlm.nih.gov/40919091/)


## 7. 내분비 축: 11β-HSD1 · 안드로겐 (The adipose–endocrine drive on the choroid plexus)

52. Sinclair AJ, Onyimba CU, Khosla P et al. **Corticosteroids, 11beta-hydroxysteroid dehydrogenase isozymes and the rabbit choroid plexus.** *Journal of neuroendocrinology* 2007. [PMID 17620103](https://pubmed.ncbi.nlm.nih.gov/17620103/)
53. Quintela T, Alves CH, Gonçalves I et al. **5Alpha-dihydrotestosterone up-regulates transthyretin levels in mice and rat choroid plexus via an androgen receptor independent pathway.** *Brain research* 2008. [PMID 18634756](https://pubmed.ncbi.nlm.nih.gov/18634756/)
54. Alves CH, Gonçalves I, Socorro S et al. **Androgen receptor is expressed in murine choroid plexus and downregulated by 5alpha-dihydrotestosterone in male and female mice.** *Journal of molecular neuroscience : MN* 2009. [PMID 19015999](https://pubmed.ncbi.nlm.nih.gov/19015999/)
55. Sinclair AJ, Walker EA, Burdon MA et al. **Cerebrospinal fluid corticosteroid levels and cortisol metabolism in patients with idiopathic intracranial hypertension: a link between 11beta-HSD1 and intracranial pressure regulation?.** *The Journal of clinical endocrinology and metabolism* 2010. [PMID 20826586](https://pubmed.ncbi.nlm.nih.gov/20826586/)
56. O'Reilly MW, Westgate CS, Hornby C et al. **A unique androgen excess signature in idiopathic intracranial hypertension is linked to cerebrospinal fluid dynamics.** *JCI insight* 2019. [PMID 30753168](https://pubmed.ncbi.nlm.nih.gov/30753168/)
57. Hardy RS, Botfield H, Markey K et al. **11βHSD1 Inhibition with AZD4017 Improves Lipid Profiles and Lean Muscle Mass in Idiopathic Intracranial Hypertension.** *The Journal of clinical endocrinology and metabolism* 2021. [PMID 33098644](https://pubmed.ncbi.nlm.nih.gov/33098644/)
58. Savaşcı D, Toydemir HE, Ekin M et al. **Hyperandrogenism and polycystic ovary syndrome phenotypes in idiopathic intracranial hypertension.** *Irish journal of medical science* 2026. [PMID 42458175](https://pubmed.ncbi.nlm.nih.gov/42458175/)


## 8. 임상시험과 약물치료 (Randomised trials and drug therapy)

59. Bietti G, Virno M, Pecori-Giraldi J. **Acetazolamide, metabolic acidosis, and intraocular pressure.** *American journal of ophthalmology* 1975. [PMID 1163584](https://pubmed.ncbi.nlm.nih.gov/1163584/)
60. Vorstrup S, Henriksen L, Paulson OB. **Effect of acetazolamide on cerebral blood flow and cerebral metabolic rate for oxygen.** *The Journal of clinical investigation* 1984. [PMID 6501565](https://pubmed.ncbi.nlm.nih.gov/6501565/)
61. Yano I, Takayama A, Takano M et al. **Pharmacokinetics and pharmacodynamics of acetazolamide in patients with transient intraocular pressure elevation.** *European journal of clinical pharmacology* 1998. [PMID 9591933](https://pubmed.ncbi.nlm.nih.gov/9591933/)
62. NORDIC Idiopathic Intracranial Hypertension Study Group Writing Committee, Wall M, McDermott MP et al. **Effect of acetazolamide on visual function in patients with idiopathic intracranial hypertension and mild visual loss: the idiopathic intracranial hypertension treatment trial.** *JAMA* 2014. [PMID 24756514](https://pubmed.ncbi.nlm.nih.gov/24756514/)
63. Wall M, Kupersmith MJ, Kieburtz KD et al. **The idiopathic intracranial hypertension treatment trial: clinical profile at baseline.** *JAMA neurology* 2014. [PMID 24756302](https://pubmed.ncbi.nlm.nih.gov/24756302/)
64. ten Hove MW, Friedman DI, Patel AD et al. **Safety and Tolerability of Acetazolamide in the Idiopathic Intracranial Hypertension Treatment Trial.** *Journal of neuro-ophthalmology : the official journal of the North American Neuro-Ophthalmology Society* 2016. [PMID 26587993](https://pubmed.ncbi.nlm.nih.gov/26587993/)
65. Smith SV, Friedman DI. **The Idiopathic Intracranial Hypertension Treatment Trial: A Review of the Outcomes.** *Headache* 2017. [PMID 28758206](https://pubmed.ncbi.nlm.nih.gov/28758206/)
66. Scotton WJ, Botfield HF, Westgate CS et al. **Topiramate is more effective than acetazolamide at lowering intracranial pressure.** *Cephalalgia : an international journal of headache* 2019. [PMID 29898611](https://pubmed.ncbi.nlm.nih.gov/29898611/)
67. Gulati S, Aref AA. **Oral acetazolamide for intraocular pressure lowering: balancing efficacy and safety in ophthalmic practice.** *Expert review of clinical pharmacology* 2021. [PMID 34003717](https://pubmed.ncbi.nlm.nih.gov/34003717/)
68. Grech O, Mitchell JL, Lyons HS et al. **Effect of glucagon like peptide-1 receptor agonist exenatide, used as an intracranial pressure lowering agent, on cognition in Idiopathic Intracranial Hypertension.** *Eye (London, England)* 2024. [PMID 38212401](https://pubmed.ncbi.nlm.nih.gov/38212401/)
69. Almaqhawi A, Alokley A, Alamri R et al. **Effectiveness of Topiramate Versus Acetazolamide in the Management of Idiopathic Intracranial Hypertension: ASystematic Review and Meta-Analysis.** *Medicina (Kaunas, Lithuania)* 2025. [PMID 40142261](https://pubmed.ncbi.nlm.nih.gov/40142261/)
70. Mitchell JL, Lyons HS, Walker JK et al. **A randomized sequential cross-over trial evaluating five purportedly ICP-lowering drugs in idiopathic intracranial hypertension.** *Headache* 2025. [PMID 39853738](https://pubmed.ncbi.nlm.nih.gov/39853738/)
71. Ahmed W, Gandhi OH, Yu N et al. **Efficacy of glucagon-like peptide-1 receptor agonists in idiopathic intracranial hypertension: A systematic review and meta-analysis.** *Journal of the neurological sciences* 2026. [PMID 41468715](https://pubmed.ncbi.nlm.nih.gov/41468715/)


## 9. 체중감량 · 대사수술 (Weight loss and bariatric surgery)

72. Sinclair AJ, Burdon MA, Nightingale PG et al. **Low energy diet and intracranial pressure in women with idiopathic intracranial hypertension: prospective cohort study.** *BMJ (Clinical research ed.)* 2010. [PMID 20610512](https://pubmed.ncbi.nlm.nih.gov/20610512/)
73. Sun WYL, Switzer NJ, Dang JT et al. **Idiopathic intracranial hypertension and bariatric surgery: a systematic review.** *Canadian journal of surgery. Journal canadien de chirurgie* 2020. [PMID 32195557](https://pubmed.ncbi.nlm.nih.gov/32195557/)
74. Mollan SP, Mitchell JL, Ottridge RS et al. **Effectiveness of Bariatric Surgery vs Community Weight Management Intervention for the Treatment of Idiopathic Intracranial Hypertension: A Randomized Clinical Trial.** *JAMA neurology* 2021. [PMID 33900360](https://pubmed.ncbi.nlm.nih.gov/33900360/)
75. Rubino D, Abrahamsson N, Davies M et al. **Effect of Continued Weekly Subcutaneous Semaglutide vs Placebo on Weight Loss Maintenance in Adults With Overweight or Obesity: The STEP 4 Randomized Clinical Trial.** *JAMA* 2021. [PMID 33755728](https://pubmed.ncbi.nlm.nih.gov/33755728/)


## 10. 유두부종 · OCT · 시기능 (Papilloedema, OCT and visual function)

76. Radius RL, Anderson DR. **Fast axonal transport in early experimental disc edema.** *Investigative ophthalmology & visual science* 1980. [PMID 6153175](https://pubmed.ncbi.nlm.nih.gov/6153175/)
77. Fard MA, Fakhree S, Abdi P et al. **Quantification of peripapillary total retinal volume in pseudopapilledema and mild papilledema using spectral-domain optical coherence tomography.** *American journal of ophthalmology* 2014. [PMID 24727146](https://pubmed.ncbi.nlm.nih.gov/24727146/)
78. Optical Coherence Tomography Substudy Committee, NORDIC Idiopathic Intracranial Hypertension Study Group. **Papilledema Outcomes from the Optical Coherence Tomography Substudy of the Idiopathic Intracranial Hypertension Treatment Trial.** *Ophthalmology* 2015. [PMID 26198807](https://pubmed.ncbi.nlm.nih.gov/26198807/)
79. Zhang Z, Wu S, Jonas JB et al. **Dynein, kinesin and morphological changes in optic nerve axons in a rat model with cerebrospinal fluid pressure reduction: the Beijing Intracranial and Intraocular Pressure (iCOP) study.** *Acta ophthalmologica* 2016. [PMID 26178710](https://pubmed.ncbi.nlm.nih.gov/26178710/)
80. Frisén L. **Swelling of the Optic Nerve Head: A Backstage View of a Staging Scheme.** *Journal of neuro-ophthalmology : the official journal of the North American Neuro-Ophthalmology Society* 2017. [PMID 28187078](https://pubmed.ncbi.nlm.nih.gov/28187078/)
81. Price DA, Harris A, Siesky B et al. **The Influence of Translaminar Pressure Gradient and Intracranial Pressure in Glaucoma: A Review.** *Journal of glaucoma* 2020. [PMID 31809396](https://pubmed.ncbi.nlm.nih.gov/31809396/)
82. Hao J, Pircher A, Miller NR et al. **Cerebrospinal fluid and optic nerve sheath compartment syndrome: A common pathophysiological mechanism in five different cases?.** *Clinical & experimental ophthalmology* 2020. [PMID 31648390](https://pubmed.ncbi.nlm.nih.gov/31648390/)
83. Sibony PA, Kupersmith MJ, Kardon RH. **Optical Coherence Tomography Neuro-Toolbox for the Diagnosis and Management of Papilledema, Optic Disc Edema, and Pseudopapilledema.** *Journal of neuro-ophthalmology : the official journal of the North American Neuro-Ophthalmology Society* 2021. [PMID 32909979](https://pubmed.ncbi.nlm.nih.gov/32909979/)
84. Behbehani R, Ali A, Al-Moosa A. **Course and Predictors of Visual Outcome of Idiopathic Intracranial Hypertension.** *Neuro-ophthalmology (Aeolus Press)* 2022. [PMID 35273409](https://pubmed.ncbi.nlm.nih.gov/35273409/)
85. Costello F, Hamann S. **Advantages and Pitfalls of the Use of Optical Coherence Tomography for Papilledema.** *Current neurology and neuroscience reports* 2024. [PMID 38261144](https://pubmed.ncbi.nlm.nih.gov/38261144/)


## 11. 두통과 중추 감작 (Headache, central sensitisation, medication overuse)

86. deSouza RM, Toma A, Watkins L. **Medication overuse headache - An under-diagnosed problem in shunted idiopathic intracranial hypertension patients.** *British journal of neurosurgery* 2015. [PMID 25136917](https://pubmed.ncbi.nlm.nih.gov/25136917/)
87. Wattiez AS, Sowers LP, Russo AF. **Calcitonin gene-related peptide (CGRP): role in migraine pathophysiology and therapeutic targeting.** *Expert opinion on therapeutic targets* 2020. [PMID 32003253](https://pubmed.ncbi.nlm.nih.gov/32003253/)
88. Mollan SP, Grech O, Sinclair AJ. **Headache attributed to idiopathic intracranial hypertension and persistent post-idiopathic intracranial hypertension headache: A narrative review.** *Headache* 2021. [PMID 34106464](https://pubmed.ncbi.nlm.nih.gov/34106464/)
89. Yiangou A, Mitchell JL, Fisher C et al. **Erenumab for headaches in idiopathic intracranial hypertension: A prospective open-label evaluation.** *Headache* 2021. [PMID 33316102](https://pubmed.ncbi.nlm.nih.gov/33316102/)
90. Bsteh G, Krajnc N, Zaic S et al. **Acute headache treatment in idiopathic intracranial hypertension: treating to the phenotype?.** *The journal of headache and pain* 2025. [PMID 41299250](https://pubmed.ncbi.nlm.nih.gov/41299250/)
91. Mandloi S, Nisar A, Shing SR et al. **The Impact of Venous Stenting on Symptoms and Quality of Life in Patients with Idiopathic Intracranial Hypertension and Spontaneous Cerebrospinal Fluid Leak.** *Journal of neurological surgery. Part B, Skull base* 2025. [PMID 41140429](https://pubmed.ncbi.nlm.nih.gov/41140429/)


## 12. 시술: 스텐트 · 션트 · ONSF (Procedures: stenting, shunting, sheath fenestration)

92. Manfield JH, Yu KK, Efthimiou E et al. **Bariatric Surgery or Non-surgical Weight Loss for Idiopathic Intracranial Hypertension? A Systematic Review and Comparison of Meta-analyses.** *Obesity surgery* 2017. [PMID 27981458](https://pubmed.ncbi.nlm.nih.gov/27981458/)
93. Asif H, Craven CL, Siddiqui AH et al. **Idiopathic intracranial hypertension: 120-day clinical, radiological, and manometric outcomes after stent insertion into the dural venous sinus.** *Journal of neurosurgery* 2018. [PMID 28984521](https://pubmed.ncbi.nlm.nih.gov/28984521/)
94. Nicholson P, Brinjikji W, Radovanovic I et al. **Venous sinus stenting for idiopathic intracranial hypertension: a systematic review and meta-analysis.** *Journal of neurointerventional surgery* 2019. [PMID 30166333](https://pubmed.ncbi.nlm.nih.gov/30166333/)
95. Yiangou A, Mitchell J, Markey KA et al. **Therapeutic lumbar puncture for headache in idiopathic intracranial hypertension: Minimal gain, is it worth the pain?.** *Cephalalgia : an international journal of headache* 2019. [PMID 29911422](https://pubmed.ncbi.nlm.nih.gov/29911422/)
96. Greener DL, Akarca D, Durnford AJ et al. **Idiopathic Intracranial Hypertension: Shunt Failure and the Role of Obesity.** *World neurosurgery* 2020. [PMID 31954904](https://pubmed.ncbi.nlm.nih.gov/31954904/)
97. Fiani B, Kondilis A, Doan T et al. **Venous sinus stenting for intractable pulsatile tinnitus: A review of indications and outcomes.** *Surgical neurology international* 2021. [PMID 33767885](https://pubmed.ncbi.nlm.nih.gov/33767885/)
98. Kaur N, Patro SK, Gupta AK et al. **Idiopathic Intracranial Hypertension and Endoscopic Optic Nerve Sheath Fenestration.** *Indian journal of otolaryngology and head and neck surgery : official publication of the Association of Otolaryngologists of India* 2022. [PMID 36452828](https://pubmed.ncbi.nlm.nih.gov/36452828/)
99. Subramanian PS, Miller NR. **Optic nerve sheath fenestration: Does it still have a role in treating patients with elevated intracranial pressure?.** *Clinical & experimental ophthalmology* 2023. [PMID 37314300](https://pubmed.ncbi.nlm.nih.gov/37314300/)
100. El-Hajj VG, Roy J, Musmar B et al. **Venous sinus stenting versus ventriculoperitoneal shunting for idiopathic intracranial hypertension: propensity score weighted, cost consequence analysis.** *Journal of neurointerventional surgery* 2026. [PMID 41781209](https://pubmed.ncbi.nlm.nih.gov/41781209/)
101. Mehta A, Goldman D, Philbrick B et al. **Time is vision: a systematic review of urgent venous sinus stenting for fulminant idiopathic intracranial hypertension.** *Frontiers in neurology* 2026. [PMID 42147834](https://pubmed.ncbi.nlm.nih.gov/42147834/)


## 13. 전격형 · 특수 집단 · 이차성 (Fulminant disease, special populations, secondary causes)

102. Lee AG. **Pseudotumor cerebri after treatment with tetracycline and isotretinoin for acne.** *Cutis* 1995. [PMID 7634848](https://pubmed.ncbi.nlm.nih.gov/7634848/)
103. Dinkin M, Oliveira C. **Men Are from Mars, Idiopathic Intracranial Hypertension Is from Venous: The Role of Venous Sinus Stenosis and Stenting in Idiopathic Intracranial Hypertension.** *Seminars in neurology* 2019. [PMID 31847040](https://pubmed.ncbi.nlm.nih.gov/31847040/)
104. Gaier ED, Heidary G. **Pediatric Idiopathic Intracranial Hypertension.** *Seminars in neurology* 2019. [PMID 31847041](https://pubmed.ncbi.nlm.nih.gov/31847041/)
105. Bouffard MA. **Fulminant Idiopathic Intracranial Hypertension.** *Current neurology and neuroscience reports* 2020. [PMID 32219578](https://pubmed.ncbi.nlm.nih.gov/32219578/)
106. Scott C, Kaliaperumal C. **Idiopathic intracranial hypertension and pregnancy: A comprehensive review of management.** *Clinical neurology and neurosurgery* 2022. [PMID 35461091](https://pubmed.ncbi.nlm.nih.gov/35461091/)
107. Angelette AL, Rando LL, Wadhwa RD et al. **Tetracycline-, Doxycycline-, Minocycline-Induced Pseudotumor Cerebri and Esophageal Perforation.** *Advances in therapy* 2023. [PMID 36763302](https://pubmed.ncbi.nlm.nih.gov/36763302/)


## 14. QSP 모델링 방법론 (QSP and modelling methodology)

108. Giulioni M, Ursino M, Alvisi C. **Correlations among intracranial pulsatility, intracranial hemodynamics, and transcranial Doppler wave form: literature review and hypothesis for future studies.** *Neurosurgery* 1988. [PMID 3288898](https://pubmed.ncbi.nlm.nih.gov/3288898/)
109. Lakin WD, Stevens SA, Tranmer BI et al. **A whole-body mathematical model for intracranial pressure dynamics.** *Journal of mathematical biology* 2003. [PMID 12673511](https://pubmed.ncbi.nlm.nih.gov/12673511/)
110. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT: pharmacometrics & systems pharmacology* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
111. Bai JP, Wang J, Zhang Y et al. **Quantitative Systems Pharmacology for Rare Disease Drug Development.** *Journal of pharmaceutical sciences* 2023. [PMID 37422281](https://pubmed.ncbi.nlm.nih.gov/37422281/)
112. Lesko LJ, van der Graaf PH. **Reflections on Model-Informed Drug Development.** *Clinical pharmacology and therapeutics* 2024. [PMID 39012331](https://pubmed.ncbi.nlm.nih.gov/39012331/)


---

## 모델 파라미터가 어느 문헌에서 왔는가 (Parameter provenance)

| 파라미터 | 값 | 근거 |
|---|---|---|
| `FFORM0` | 0.35 mL/min (504 mL/day) | 고전적 인체 CSF 생성률; phase-contrast MRI 측정치와 정합 (Cutler 1968; Radiology 2004) |
| `E1` | 0.20 /mL | Marmarou 지수 압력-용적 관계 (PVI); ICP 10 mmHg에서 순응도 0.5 mL/mmHg |
| `ROUTB`, `KROUTP` | 0, 1.69 mmHg·min/mL per mmHg | R_out을 독립 병변이 아니라 정맥압 의존 저항으로 둠 → 마른 사람 8.8, BMI 38에서 17.1 (주입검사 보고치 범위) |
| `RTS0` | 0.00333 → 2 mmHg 기울기 | 정상 상시상정맥동–경정맥 압력차 < 4 mmHg |
| `PCRIT`, `PSTIFF` | 10.95, 4.14 mmHg | 탐색이 아니라 **연립해**: BMI 38에서 ICP 25.9 mmHg, 루프이득 0.53, 기울기 8.4 mmHg 세 조건에서 유도 |
| `KIAP` | 0.35 mmHg per BMI unit | 복강내압 → 흉강내압 → 중심정맥압 전달 (Sugerman 1997) |
| `EMAXCA` | 0.42 | 아세타졸아미드의 CSF 생성 억제 상한 (30–50 % 보고 범위의 중앙) |
| `EMAXGLP` | 0.30 | GLP-1R 작용제의 맥락총 분비 억제 (Jensen 2025 rat; Mitchell/Grech 2024 human) |
| `EMAXGC`, `KANDR` | 0.10, 0.012/BMI unit | 지방–내분비 구동을 11β-HSD1 의존분(차단 가능)과 안드로겐 의존분(차단 불가)으로 분리 → AZD4017 음성 결과가 **입력이 아니라 예측**이 됨 |
| `EMAXIOP` | 2.5 mmHg | 탄산탈수효소 억제의 안압 강하; 녹내장 고용량치보다 보수적으로 |
| `TLPDTHR` | 3 mmHg | 축삭형질 수송 정체가 시작되는 층판간 압력차 (Radius & Anderson 1980) |
| `KAXL`, `KISCH` | 0.0012, 0.006 /day | 만성 정체성 축삭 소실과 TLPD > 22 mmHg에서의 허혈성 소실을 분리 → 전격형이 별도 파라미터 없이 생성됨 |
| `FMOHDES` | 0.60 | 약물과용이 중추 감작의 소실을 차단 → ICP 정상화 후에도 두통이 남는 구조 |

## 보정 앵커 (Calibration anchors actually used)

| 시험 | 발표값 | 모델값 |
|---|---|---|
| Sinclair 2010 BMJ, 저열량식 −15.7 %, 3개월 | ICP −8.4 cmH2O | −6.0 cmH2O (90일) / −7.1 (180일) |
| Mitchell·Grech 2024 (exenatide), 12주 | ICP −5.6 cmH2O | −4.9 cmH2O (동시대조 대비 −6.0) |
| Hardy·Markey 2021 (AZD4017), 12주 | 유의한 ICP 변화 없음 | −0.3 cmH2O |
| IIHTT 2014 (Wall), 6개월 | 유두부종 등급 차 −0.70 | −0.50 |
| IIHTT 2014 (Wall), 6개월 | 시야 MD 차 +0.71 dB | +0.89 dB |
| IIHTT 2014 (Wall), 6개월 | ICP 차 −4.4 cmH2O | −6.1 cmH2O (**약 1.4배 과대예측 — 보정하지 않고 기록**) |
| IIH:WT 2021 (Mollan), 12개월 | 대사수술 vs 식이 ICP 차 −6.0 cmH2O | 대사수술 −10.0 cmH2O (제시시점 대비) |

