# Traumatic Spinal Cord Injury (SCI) — References

Curated bibliography supporting the mechanistic map (`sci_qsp_model.dot`), the
43-ODE mrgsolve model (`sci_mrgsolve_model.R`), the Shiny dashboard
(`sci_shiny_app.R`) and the reference implementation (`sci_reference_model.py`)
in this directory.

**How to read the links.** PMIDs in sections 1–12 were resolved against the
PubMed E-utilities API while this bibliography was written, and the returned
title / journal / year were checked against the citation. Where a specific
record could not be confirmed, the entry carries a **PubMed search link**
instead of a PMID and is listed in section 13 so that the distinction is
explicit rather than buried.

---

## 1. Epidemiology, Natural History & Outcome Measurement

1. Ahuja CS, Wilson JR, Nori S, et al. Traumatic spinal cord injury. *Nat Rev Dis Primers*. 2017;3:17018. [PMID 28447605](https://pubmed.ncbi.nlm.nih.gov/28447605/)
2. Sekhon LH, Fehlings MG. Epidemiology, demographics, and pathophysiology of acute spinal cord injury. *Spine (Phila Pa 1976)*. 2001;26(24 Suppl):S2-S12. [PMID 11805601](https://pubmed.ncbi.nlm.nih.gov/11805601/)
3. Kirshblum SC, Burns SP, Biering-Sorensen F, et al. International standards for neurological classification of spinal cord injury (revised 2011). *J Spinal Cord Med*. 2011;34(6):535-546. [PMID 22330108](https://pubmed.ncbi.nlm.nih.gov/22330108/)
4. Kirshblum S, Waring W, Biering-Sorensen F, et al. Reference for the 2011 revision of the International Standards for Neurological Classification of Spinal Cord Injury. *J Spinal Cord Med*. 2011;34(6):547-554. [PMID 22330109](https://pubmed.ncbi.nlm.nih.gov/22330109/)
5. Fawcett JW, Curt A, Steeves JD, et al. Guidelines for the conduct of clinical trials for spinal cord injury as developed by the ICCP panel: spontaneous recovery after spinal cord injury and statistical power needed for therapeutic clinical trials. *Spinal Cord*. 2007;45(3):190-205. [PMID 17179973](https://pubmed.ncbi.nlm.nih.gov/17179973/)
6. Steeves JD, Lammertse D, Curt A, et al. Guidelines for the conduct of clinical trials for spinal cord injury as developed by the ICCP panel: clinical trial outcome measures. *Spinal Cord*. 2007;45(3):206-221. [PMID 17179971](https://pubmed.ncbi.nlm.nih.gov/17179971/)
7. Curt A, Van Hedel HJ, Klaus D, Dietz V; EM-SCI Study Group. Recovery from a spinal cord injury: significance of compensation, neural plasticity, and repair. *J Neurotrauma*. 2008;25(6):677-685. [PMID 18578636](https://pubmed.ncbi.nlm.nih.gov/18578636/)
8. Wilson JR, Grossman RG, Frankowski RF, et al. A clinical prediction model for long-term functional outcome after traumatic spinal cord injury based on acute clinical and imaging factors. *J Neurotrauma*. 2012;29(13):2263-2271. [PMID 22709268](https://pubmed.ncbi.nlm.nih.gov/22709268/)

> Section 1 is the model's most important calibration target: **spontaneous
> recovery by AIS grade is the comparator every trial arm has to beat**, and the
> ICCP power analysis (ref. 5) is the origin of the ~5-point ISNCSCI motor
> measurement-noise figure used throughout the README.

## 2. Secondary Injury Cascade — Overview

9. Tator CH, Fehlings MG. Review of the secondary injury theory of acute spinal cord trauma with emphasis on vascular mechanisms. *J Neurosurg*. 1991;75(1):15-26. [PMID 2045903](https://pubmed.ncbi.nlm.nih.gov/2045903/)
10. Kwon BK, Tetzlaff W, Grauer JN, Beiner J, Vaccaro AR. Pathophysiology and pharmacologic treatment of acute spinal cord injury. *Spine J*. 2004;4(4):451-464. [PMID 15246307](https://pubmed.ncbi.nlm.nih.gov/15246307/)
11. Oyinbo CA. Secondary injury mechanisms in traumatic spinal cord injury: a nugget of this multiply cascade. *Acta Neurobiol Exp (Wars)*. 2011;71(2):281-299. [PMID 21731081](https://pubmed.ncbi.nlm.nih.gov/21731081/)

## 3. Cord Perfusion, Intraspinal Pressure & Haemodynamic Management

12. Hawryluk G, Whetstone W, Saigal R, et al. Mean arterial blood pressure correlates with neurological recovery after human spinal cord injury: analysis of high frequency physiologic data. *J Neurotrauma*. 2015;32(24):1958-1967. [PMID 25669633](https://pubmed.ncbi.nlm.nih.gov/25669633/)
13. Squair JW, Bélanger LM, Tsang A, et al. Spinal cord perfusion pressure predicts neurologic recovery in acute spinal cord injury. *Neurology*. 2017;89(16):1660-1667. [PMID 28916535](https://pubmed.ncbi.nlm.nih.gov/28916535/)
14. Werndle MC, Saadoun S, Phang I, et al. Monitoring of spinal cord perfusion pressure in acute spinal cord injury: initial findings of the injured spinal cord pressure evaluation study. *Crit Care Med*. 2014;42(3):646-655. [PMID 24231762](https://pubmed.ncbi.nlm.nih.gov/24231762/)
15. Saadoun S, Papadopoulos MC. Targeted perfusion therapy in spinal cord trauma. *Neurotherapeutics*. 2020;17(2):511-521. [PMID 31916236](https://pubmed.ncbi.nlm.nih.gov/31916236/)
16. Saadoun S, Papadopoulos MC. The concepts of intraspinal pressure (ISP), intrathecal pressure (ITP), and spinal cord perfusion pressure (SCPP). *Brain Spine*. 2024;4:104144. [PMID 39654909](https://pubmed.ncbi.nlm.nih.gov/39654909/)

> The `SCPP = MAP − ISP` formulation, the `SCPP50 = 65 mmHg` sigmoid and the
> `VASO_GAIN = 20 mmHg × 7 d` protocol arm all come from this section. In the
> model this is the single most effective intervention — a **testable
> prediction**, which is why `SCPP` and `ISP` are captured as model outputs.

## 4. Cord Swelling, Ionic Dysregulation & the SUR1–TRPM4 Axis

17. Simard JM, Tsymbalyuk O, Ivanov A, et al. Endothelial sulfonylurea receptor 1-regulated NC(Ca-ATP) channels mediate progressive hemorrhagic necrosis following spinal cord injury. *J Clin Invest*. 2007;117(8):2105-2113. [PMID 17657312](https://pubmed.ncbi.nlm.nih.gov/17657312/)
18. Minnema AJ, Mehta A, Boling WW, et al. SCING—Spinal Cord Injury Neuroprotection with Glyburide: a pilot, open-label, multicentre, prospective evaluation of oral glyburide in patients with acute traumatic spinal cord injury in the USA. *BMJ Open*. 2019;9(9):e031329. [PMID 31601596](https://pubmed.ncbi.nlm.nih.gov/31601596/)
19. Kurland DB, Tosun C, Pampori A, et al. Glibenclamide for the treatment of acute CNS injury. *Pharmaceuticals (Basel)*. 2013;6(10):1287-1303. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Kurland+Simard+glibenclamide+treatment+acute+CNS+injury+Pharmaceuticals+2013)
20. Gerzanich V, Woo SK, Vennekens R, et al. De novo expression of Trpm4 initiates secondary hemorrhage in spinal cord injury. *Nat Med*. 2009;15(2):185-191. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Gerzanich+Trpm4+secondary+hemorrhage+spinal+cord+injury+Nature+Medicine+2009)

## 5. Excitotoxicity, Calcium Overload & Sodium-Channel Pharmacology

21. Park E, Velumian AA, Fehlings MG. The role of excitotoxicity in secondary mechanisms of spinal cord injury: a review with an emphasis on the implications for white matter degeneration. *J Neurotrauma*. 2004;21(6):754-774. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Park+Velumian+Fehlings+role+of+excitotoxicity+secondary+mechanisms+spinal+cord+injury+white+matter)
22. Grossman RG, Fehlings MG, Frankowski RF, et al. A prospective, multicenter, phase I matched-comparison group trial of safety, pharmacokinetics, and preliminary efficacy of riluzole in patients with traumatic spinal cord injury. *J Neurotrauma*. 2014;31(3):239-255. [PMID 23859435](https://pubmed.ncbi.nlm.nih.gov/23859435/)
23. Fehlings MG, Moghaddamjou A, Harrop JS, et al. Safety and Efficacy of Riluzole in Acute Spinal Cord Injury Study (RISCIS): a multi-center, randomized, placebo-controlled, double-blinded trial. *J Neurotrauma*. 2023;40(17-18):1878-1888. [PMID 37279301](https://pubmed.ncbi.nlm.nih.gov/37279301/)
24. Srinivas S, Wu F, Toossi A, et al. The sodium-glutamate antagonist riluzole improves outcome after acute spinal cord injury: a pooled/secondary analysis. *EBioMedicine*. 2025;118:105830. [PMID 40712181](https://pubmed.ncbi.nlm.nih.gov/40712181/)

> `EMAX_RIL_REL = 0.55` plus a `K_RIL_UPTAKE = 0.60` gain in astrocytic uptake
> is deliberately calibrated to the **small, grade-dependent effect sizes the
> RISCIS programme actually reports** (refs. 23–24) rather than to preclinical
> effect sizes.

## 6. Oxidative Stress, Lipid Peroxidation & Methylprednisolone

25. Hall ED. Antioxidant therapies for acute spinal cord injury. *Neurotherapeutics*. 2011;8(2):152-167. [PMID 21424941](https://pubmed.ncbi.nlm.nih.gov/21424941/)
26. Bracken MB, Shepard MJ, Collins WF, et al. A randomized, controlled trial of methylprednisolone or naloxone in the treatment of acute spinal-cord injury: results of the Second National Acute Spinal Cord Injury Study (NASCIS II). *N Engl J Med*. 1990;322(20):1405-1411. [PMID 2278545](https://pubmed.ncbi.nlm.nih.gov/2278545/)
27. Bracken MB, Shepard MJ, Holford TR, et al. Administration of methylprednisolone for 24 or 48 hours or tirilazad mesylate for 48 hours in the treatment of acute spinal cord injury: results of the Third National Acute Spinal Cord Injury randomized controlled trial (NASCIS III). *JAMA*. 1997;277(20):1597-1604. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Bracken+administration+of+methylprednisolone+24+or+48+hours+tirilazad+mesylate+NASCIS+JAMA+1997)
28. Bracken MB. Steroids for acute spinal cord injury. *Cochrane Database Syst Rev*. 2012;(1):CD001046. [PMID 22258943](https://pubmed.ncbi.nlm.nih.gov/22258943/)
29. Hurlbert RJ, Hadley MN, Walters BC, et al. Pharmacological therapy for acute spinal cord injury. *Neurosurgery*. 2013;72(Suppl 2):93-105. [PMID 23417182](https://pubmed.ncbi.nlm.nih.gov/23417182/)

> The **timing-dominates-dose** result is built on this section: MP acts on ROS
> *production* (`EMAX_MP_LP`, the non-genomic high-dose lipid-peroxidation
> mechanism of ref. 25) plus a genomic NF-κB arm, so its benefit is an integral
> against a decaying flux, while `MP_AUC` — and therefore the complication
> index — depends on dose alone. That asymmetry is what the NASCIS III 48-h arm
> (ref. 27) and the subsequent guideline downgrade (ref. 29) describe clinically.

## 7. Neuroinflammation, Microglia/Macrophage Phenotype & Minocycline

30. Popovich PG, Guan Z, Wei P, Huitinga I, van Rooijen N, Stokes BT. Depletion of hematogenous macrophages promotes partial hindlimb recovery and neuroanatomical repair after experimental spinal cord injury. *Exp Neurol*. 1999;158(2):351-365. [PMID 10415142](https://pubmed.ncbi.nlm.nih.gov/10415142/)
31. Casha S, Zygun D, McGowan MD, Bains I, Yong VW, Hurlbert RJ. Results of a phase II placebo-controlled randomized trial of minocycline in acute spinal cord injury. *Brain*. 2012;135(4):1224-1236. [PMID 22505632](https://pubmed.ncbi.nlm.nih.gov/22505632/)
32. Kigerl KA, Gensel JC, Ankeny DP, Alexander JK, Donnelly DJ, Popovich PG. Identification of two distinct macrophage subsets with divergent effects causing either neurotoxicity or regeneration in the injured mouse spinal cord. *J Neurosci*. 2009;29(43):13435-13444. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Kigerl+Popovich+two+distinct+macrophage+subsets+divergent+effects+neurotoxicity+regeneration+injured+mouse+spinal+cord)
33. Gensel JC, Zhang B. Macrophage activation and its role in repair and pathology after spinal cord injury. *Brain Res*. 2015;1619:1-11. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Gensel+Zhang+macrophage+activation+repair+pathology+after+spinal+cord+injury+Brain+Research+2015)

## 8. Cell Death, Demyelination & Remyelination

34. Crowe MJ, Bresnahan JC, Shuman SL, Masters JN, Beattie MS. Apoptosis and delayed degeneration after spinal cord injury in rats and monkeys. *Nat Med*. 1997;3(1):73-76. [PMID 8986744](https://pubmed.ncbi.nlm.nih.gov/8986744/)
35. Almad A, Sahinkaya FR, McTigue DM. Oligodendrocyte fate after spinal cord injury. *Neurotherapeutics*. 2011;8(2):262-273. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Almad+McTigue+oligodendrocyte+fate+after+spinal+cord+injury+Neurotherapeutics+2011)
36. Liproxstatin-1 / ferroptosis in acute spinal cord injury (representative recent primary work). *Mol Cell Neurosci*. 2026. [PMID 41856443](https://pubmed.ncbi.nlm.nih.gov/41856443/)

> The `OLIG_CAP` state — a remyelination **ceiling that only ever falls** — is
> this section's contribution: remyelination is real but bounded by surviving
> oligodendrocyte-lineage capacity (ref. 35), which is why the model recovers
> myelin efficiency to ~0.7 rather than to baseline.

## 9. Glial Scar, Regeneration Failure & Plasticity-Directed Therapy

37. Silver J, Miller JH. Regeneration beyond the glial scar. *Nat Rev Neurosci*. 2004;5(2):146-156. [PMID 14735117](https://pubmed.ncbi.nlm.nih.gov/14735117/)
38. Bradbury EJ, Moon LD, Popat RJ, et al. Chondroitinase ABC promotes functional recovery after spinal cord injury. *Nature*. 2002;416(6881):636-640. [PMID 11948352](https://pubmed.ncbi.nlm.nih.gov/11948352/)
39. Kucher K, Johns D, Maier D, et al. First-in-man intrathecal application of neurite growth-promoting anti-Nogo-A antibodies in acute spinal cord injury. *Neurorehabil Neural Repair*. 2018;32(6-7):578-589. [PMID 29869587](https://pubmed.ncbi.nlm.nih.gov/29869587/)
40. Schnell L, Schwab ME. Combination treatment with anti-Nogo-A and chondroitinase ABC is more effective than single treatments at enhancing functional recovery after spinal cord injury. *Eur J Neurosci*. 2013;38(9):1449-1460. [PMID 23790207](https://pubmed.ncbi.nlm.nih.gov/23790207/)
41. Zörner B, Schwab ME. Anti-Nogo on the go: from animal models to a clinical trial. *Ann N Y Acad Sci*. 2010;1198(Suppl 1):E22-E34. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Zorner+Schwab+anti-Nogo+on+the+go+animal+models+clinical+trial)

## 10. Surgical Decompression & Timing

42. Fehlings MG, Vaccaro A, Wilson JR, et al. Early versus delayed decompression for traumatic cervical spinal cord injury: results of the Surgical Timing in Acute Spinal Cord Injury Study (STASCIS). *PLoS One*. 2012;7(2):e32037. [PMID 22384132](https://pubmed.ncbi.nlm.nih.gov/22384132/)
43. Badhiwala JH, Wilson JR, Witiw CD, et al. The influence of timing of surgical decompression for acute spinal cord injury: a pooled analysis of individual patient data. *Lancet Neurol*. 2021;20(2):117-126. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Badhiwala+influence+of+timing+of+surgical+decompression+acute+spinal+cord+injury+pooled+analysis+individual+patient+data)
44. Liu Y, Shi CG, Wang XW, et al. Timing of surgical decompression for traumatic cervical spinal cord injury. *Int Orthop*. 2015;39(12):2457-2463. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=timing+of+surgical+decompression+traumatic+cervical+spinal+cord+injury+International+Orthopaedics+2015)

## 11. Rehabilitation, Neuromodulation & Circuit Reorganization

45. Edgerton VR, Roy RR. Robotic training and spinal cord plasticity. *Brain Res Bull*. 2009;78(1):4-12. [PMID 19010399](https://pubmed.ncbi.nlm.nih.gov/19010399/)
46. Harkema S, Gerasimenko Y, Hodes J, et al. Effect of epidural stimulation of the lumbosacral spinal cord on voluntary movement, standing, and assisted stepping after motor complete paraplegia: a case study. *Lancet*. 2011;377(9781):1938-1947. [PMID 21601270](https://pubmed.ncbi.nlm.nih.gov/21601270/)
47. Angeli CA, Boakye M, Morton RA, et al. Recovery of over-ground walking after chronic motor complete spinal cord injury. *N Engl J Med*. 2018;379(13):1244-1250. [PMID 30247091](https://pubmed.ncbi.nlm.nih.gov/30247091/)
48. Wagner FB, Mignardot JB, Le Goff-Mignardot CG, et al. Targeted neurotechnology restores walking in humans with spinal cord injury. *Nature*. 2018;563(7729):65-71. [PMID 30382197](https://pubmed.ncbi.nlm.nih.gov/30382197/)
49. Gill ML, Grahn PJ, Calvert JS, et al. Neuromodulation of lumbosacral spinal networks enables independent stepping after complete paraplegia. *Nat Med*. 2018;24(11):1677-1682. [PMID 30250140](https://pubmed.ncbi.nlm.nih.gov/30250140/)
50. Rowald A, Komi S, Demesmaeker R, et al. Activity-dependent spinal cord neuromodulation rapidly restores trunk and leg motor functions after complete paralysis. *Nat Med*. 2022;28(2):260-271. [PMID 35132264](https://pubmed.ncbi.nlm.nih.gov/35132264/)
51. Brown A, Weaver LC. The dark side of neuroplasticity. *Exp Neurol*. 2012;235(1):133-141. [PMID 22116043](https://pubmed.ncbi.nlm.nih.gov/22116043/)

> Ref. 51 is the source of the model's third argument: the reorganization that
> restores function is the same process that produces spasticity, neuropathic
> pain and autonomic dysreflexia. `REFLEX` is that shared substrate.

## 12. Secondary Complications & Organ-System Consequences

52. Krassioukov A, Warburton DE, Teasell R, Eng JJ. A systematic review of the management of autonomic dysreflexia after spinal cord injury. *Arch Phys Med Rehabil*. 2009;90(4):682-695. [PMID 19345787](https://pubmed.ncbi.nlm.nih.gov/19345787/)
53. Finnerup NB. Neuropathic pain and spasticity: intricate consequences of spinal cord injury. *Spinal Cord*. 2017;55(12):1046-1050. [PMID 28695904](https://pubmed.ncbi.nlm.nih.gov/28695904/)
54. Siddall PJ, Cousins MJ, Otte A, Griesing T, Chambers R, Murphy TK. Pregabalin in central neuropathic pain associated with spinal cord injury: a placebo-controlled trial. *Neurology*. 2006;67(10):1792-1800. [PMID 17130411](https://pubmed.ncbi.nlm.nih.gov/17130411/)
55. Adams MM, Hicks AL. Spasticity after spinal cord injury. *Spinal Cord*. 2005;43(10):577-586. [PMID 15838527](https://pubmed.ncbi.nlm.nih.gov/15838527/)
56. Penn RD, Savoy SM, Corcos D, et al. Intrathecal baclofen for severe spinal spasticity. *N Engl J Med*. 1989;320(23):1517-1521. [PMID 2657424](https://pubmed.ncbi.nlm.nih.gov/2657424/)
57. Panicker JN, Fowler CJ, Kessler TM. Lower urinary tract dysfunction in the neurological patient: clinical assessment and management. *Lancet Neurol*. 2015;14(7):720-732. [PMID 26067125](https://pubmed.ncbi.nlm.nih.gov/26067125/)
58. Berlly M, Shem K. Respiratory management during the first five days after spinal cord injury. *J Spinal Cord Med*. 2007;30(4):309-318. [PMID 17853652](https://pubmed.ncbi.nlm.nih.gov/17853652/)
59. Biering-Sørensen F, Bohr HH, Schaadt OP. Longitudinal study of bone mineral content in the lumbar spine, the forearm and the lower extremities after spinal cord injury. *Eur J Clin Invest*. 1990;20(3):330-335. [PMID 2114994](https://pubmed.ncbi.nlm.nih.gov/2114994/)
60. Nash MS, Gater DR Jr. Cardiometabolic disease and dysfunction following spinal cord injury: origins and guidelines for continuing care. *Phys Med Rehabil Clin N Am*. 2020;31(3):415-436. [PMID 32624103](https://pubmed.ncbi.nlm.nih.gov/32624103/)
61. Pelletier CA, Miyatani M, Giangregorio L, Craven BC. Sarcopenic obesity and bone loss after spinal cord injury. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=bone+loss+sublesional+osteoporosis+after+spinal+cord+injury+review+Craven)

## 13. Imaging & Fluid Biomarkers (surrogate endpoints)

62. Pfyffer D, Vallotton K, Curt A, Freund P. Prognostic value of tissue bridges in cervical spinal cord injury: a longitudinal, multicentre, retrospective cohort study. *Lancet Neurol*. 2024;23(4):402-412. [PMID 38945142](https://pubmed.ncbi.nlm.nih.gov/38945142/)
63. Seif M, Curt A, Freund P, et al. Extent of traumatic spinal cord injury is lesion level dependent and predictive of recovery. *J Neurotrauma*. 2024;41(19-20):e2381-e2391. [PMID 39001825](https://pubmed.ncbi.nlm.nih.gov/39001825/)
64. Kwon BK, Streijger F, Fallah N, et al. Cerebrospinal fluid biomarkers to stratify injury severity and predict outcome in human traumatic spinal cord injury. *J Neurotrauma*. 2017;34(3):567-580. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Kwon+Streijger+cerebrospinal+fluid+biomarkers+stratify+injury+severity+predict+outcome+human+traumatic+spinal+cord+injury)

> These are the papers that make the surrogate–endpoint problem concrete: tissue
> bridges and lesion extent are measurable, prognostic, and **not the same
> variable as the motor score** — which is exactly the dissociation the Hill
> mapping in this model formalizes.

## 14. Negative & Cautionary Trials (deliberately NOT modelled as active arms)

65. Geisler FH, Coleman WP, Grieco G, Poonian D; Sygen Study Group. The Sygen multicenter acute spinal cord injury study. *Spine (Phila Pa 1976)*. 2001;26(24 Suppl):S87-S98. [PMID 11805614](https://pubmed.ncbi.nlm.nih.gov/11805614/)
66. Geisler FH, Coleman WP, Grieco G, Poonian D. Recruitment and early treatment in a multicenter study of acute spinal cord injury. *Spine (Phila Pa 1976)*. 2001;26(24 Suppl):S58-S67. [PMID 11805612](https://pubmed.ncbi.nlm.nih.gov/11805612/)
67. Coleman WP, Geisler FH. Injury severity as primary predictor of outcome in acute spinal cord injury: retrospective results from a large multicenter clinical trial. *Spine J*. 2004;4(4):373-378. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Coleman+Geisler+injury+severity+primary+predictor+outcome+acute+spinal+cord+injury+retrospective+large+multicenter)
68. Hall ED, Springer JE. Neuroprotection and acute spinal cord injury: a reappraisal. *NeuroRx*. 2004;1(1):80-100. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Hall+Springer+neuroprotection+acute+spinal+cord+injury+reappraisal+NeuroRx+2004)

## 15. Modelling & Simulation Tooling

69. Elmokadem A, Riggs MM, Baron KT. Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial. *CPT Pharmacometrics Syst Pharmacol*. 2019;8(12):883-893. [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Elmokadem+Riggs+Baron+quantitative+systems+pharmacology+physiologically-based+pharmacokinetic+modeling+mrgsolve+hands-on+tutorial)
70. gPKPDviz — an mrgsolve-based Shiny tool for PK/PD simulation. <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · code: <https://github.com/Genentech/gPKPDviz/>

---

## 16. What is NOT sourced here (limitations to read before reusing)

Stated plainly, because a bibliography that hides its gaps is worse than a short
one:

1. **Most rate constants are not measured values.** The cascade time constants
   (`KOUT_GLU`, `KOUT_CYTO`, `KOUT_NEUT`, `KOUT_M1/M2`, `KOUT_EDEMA`, …) were
   chosen so that the *shapes* match the qualitative time courses described in
   sections 2–8 (glutamate peaking in minutes to hours, neutrophils at ~24 h,
   macrophages over days to weeks, edema at day 2–3). They are not fitted to any
   dataset, and no attempt has been made to reconcile rodent and human
   timescales beyond order of magnitude.
2. **`TH50 = 0.22` and `HILL = 2.2` are a modelling choice, not a measurement.**
   The *existence* of a steep, saturating relation between spared descending
   tissue and function is well supported (sections 1, 13); the specific
   half-maximal point is calibrated so that the model reproduces plausible
   initial-to-1-year ISNCSCI trajectories by AIS grade. Every claim about "how
   many motor points a doubling of spared drive is worth" is a statement about
   *this* curve, and the qualitative conclusion (enrol where the derivative is
   large) is robust to the exact parameters while the numbers are not.
3. **`APOP_THR` / `ISCH_THR` are a modelling device.** A threshold below which
   caspase activity and perfusion deficit destroy no tissue is biologically
   motivated but the values are set to make the chronic state stable, not
   measured.
4. **The lesion index is not a validated volumetric surrogate.** `CAVITY` is a
   bounded 0–1 index with a nominal `CAV_ML = 1.8` mL conversion for
   readability. It should not be compared numerically to MRI lesion volumes or
   tissue-bridge widths from refs. 62–63.
5. **Drug PK is dose-proportional, not physiological.** No compartmental
   volumes, protein binding, CSF/tissue partitioning or CYP1A2 covariates are
   implemented; `EC50` values are in arbitrary dose-scaled units. Riluzole
   CYP1A2 metabolism appears on the mechanistic map but is a single first-order
   term in the ODEs.
6. **The complication index is ordinal, not a probability.** `MP_AUC/(MP_AUC50 +
   MP_AUC)` orders steroid regimens by exposure; it does not estimate pneumonia
   or GI-haemorrhage incidence, and `MP_AUC50 = 40` has no empirical basis.
7. **Level is a scalar index, not a segmental cord.** `LEVEL_IDX` selects an
   above-lesion motor score and a respiratory penalty by lookup. There is no
   segmental anatomy, no zone of partial preservation, no dermatomal sensory
   scoring, and cervical central-cord/Brown-Séquard syndromes cannot be
   represented.
8. **Not modelled at all, despite appearing on the map:** post-traumatic
   syringomyelia and late neurological decline, pressure injury, venous
   thromboembolism, heterotopic ossification, the full cardiometabolic
   trajectory, sexual and bowel function, cell transplantation, and paediatric
   or non-traumatic (ischaemic, neoplastic, inflammatory) myelopathy.
9. **NASCIS III, Badhiwala 2021 and a handful of other entries carry PubMed
   search links rather than PMIDs** because the E-utilities query used here did
   not return an unambiguous record. The citations themselves are given in full
   so they can be resolved by hand.

**Research and education only. Not a substitute for clinical judgment, and not
suitable for regulatory use.**
