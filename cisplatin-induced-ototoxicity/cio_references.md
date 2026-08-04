# 시스플라틴 유발 이독성 (Cisplatin-Induced Ototoxicity, CIO) — 참고문헌

> 이 목록의 모든 항목은 NCBI E-utilities 로 PubMed 를 직접 조회하여 반환된 레코드이며,
> PMID · 제목 · 저널 · 연도는 조회 결과를 그대로 옮긴 것입니다 (조회일 2026-08-04).
> 기억에 의존해 PMID 를 적은 항목은 없습니다.

각 절의 머리말은 그 문헌군이 모델의 **어느 파라미터 또는 어느 주장**을 떠받치는지를 밝힙니다.
모델 자체는 [`cio_mrgsolve_model.R`](cio_mrgsolve_model.R), 검증 실행은
[`cio_reference_output.txt`](cio_reference_output.txt) 를 참조하십시오.

---

## A. 백금 축적과 와우 약동학 (Cochlear platinum kinetics)

> 모델의 첫 번째 구조적 선택(백금은 농도가 아니라 저장고이며, 혈장→혈관조→외림프 캐스케이드가 6시간 창을 만든다)의 근거.

1. Breglio AM et al. *Cisplatin is retained in the cochlea indefinitely following chemotherapy*. Nat Commun 2017. [PMID 29162831](https://pubmed.ncbi.nlm.nih.gov/29162831/)
2. Saito T et al. *The effect of sodium thiosulfate on ototoxicity and pharmacokinetics after cisplatin treatment in guinea pigs*. Eur Arch Otorhinolaryngol 1997. [PMID 9248736](https://pubmed.ncbi.nlm.nih.gov/9248736/)
3. Hellberg V et al. *Cochlear pharmacokinetics of cisplatin: an in vivo study in the guinea pig*. Laryngoscope 2013. [PMID 23754209](https://pubmed.ncbi.nlm.nih.gov/23754209/)
4. Miettinen S et al. *Blood flow-independent accumulation of cisplatin in the guinea pig cochlea*. Acta Otolaryngol 1997. [PMID 9039482](https://pubmed.ncbi.nlm.nih.gov/9039482/)
5. Laurell G et al. *Distribution of cisplatin in perilymph and cerebrospinal fluid after intravenous administration in the guinea pig*. Cancer Chemother Pharmacol 1995. [PMID 7720182](https://pubmed.ncbi.nlm.nih.gov/7720182/)
6. Lyu AR et al. *CORM‑2 reduces cisplatin accumulation in the mouse inner ear and protects against cisplatin-induced ototoxicity*. J Adv Res 2024. [PMID 38030129](https://pubmed.ncbi.nlm.nih.gov/38030129/)
7. Thomas JP et al. *High accumulation of platinum-DNA adducts in strial marginal cells of the cochlea is an early event in cisplatin but not carboplatin ototoxicity*. Mol Pharmacol 2006. [PMID 16569706](https://pubmed.ncbi.nlm.nih.gov/16569706/)
8. Ruggiero A et al. *Progressive Sensorineural Hearing Loss Following Cisplatin Chemotherapy: Mechanisms Underlying Cochlear Retention and Long-Term Ototoxicity*. Pharmaceuticals (Basel) 2026. [PMID 42198452](https://pubmed.ncbi.nlm.nih.gov/42198452/)
9. Köppen C et al. *Quantitative imaging of platinum based on laser ablation-inductively coupled plasma-mass spectrometry to investigate toxic side effects of cisplatin*. Metallomics 2015. [PMID 26477751](https://pubmed.ncbi.nlm.nih.gov/26477751/)
10. Ke Y et al. *The Breakdown of Blood-Labyrinth Barrier Makes it Easier for Drugs to Enter the Inner Ear*. Laryngoscope 2024. [PMID 37987231](https://pubmed.ncbi.nlm.nih.gov/37987231/)
11. Sekulic M et al. *Human blood-labyrinth barrier model to study the effects of cytokines and inflammation*. Front Mol Neurosci 2023. [PMID 37808472](https://pubmed.ncbi.nlm.nih.gov/37808472/)

## B. 수송체와 세포 내 유입 (Transporters and cellular entry)

> 와우 유입 경로와 KINBL·KUPT 파라미터의 근거.

12. Schoeberl A et al. *The copper transporter CTR1 and cisplatin accumulation at the single-cell level by LA-ICP-TOFMS*. Front Mol Biosci 2022. [PMID 36518851](https://pubmed.ncbi.nlm.nih.gov/36518851/)
13. Lin X et al. *The copper transporter CTR1 regulates cisplatin uptake in Saccharomyces cerevisiae*. Mol Pharmacol 2002. [PMID 12391279](https://pubmed.ncbi.nlm.nih.gov/12391279/)
14. Wang X et al. *Copper transporter Ctr1 contributes to enhancement of the sensitivity of cisplatin in esophageal squamous cell carcinoma*. Transl Oncol 2023. [PMID 36689863](https://pubmed.ncbi.nlm.nih.gov/36689863/)
15. Pabla N et al. *The copper transporter Ctr1 contributes to cisplatin uptake by renal tubular cells during cisplatin nephrotoxicity*. Am J Physiol Renal Physiol 2009. [PMID 19144690](https://pubmed.ncbi.nlm.nih.gov/19144690/)
16. Lanvers-Kaminsky C et al. *Human OCT2 variant c.808G>T confers protection effect against cisplatin-induced ototoxicity*. Pharmacogenomics 2015. [PMID 25823781](https://pubmed.ncbi.nlm.nih.gov/25823781/)
17. Maruyama A et al. *Susceptibility of mouse cochlear hair cells to cisplatin ototoxicity largely depends on sensory mechanoelectrical transduction channels both Ex Vivo and In Vivo*. Hear Res 2024. [PMID 38718672](https://pubmed.ncbi.nlm.nih.gov/38718672/)
18. Hellberg V et al. *Immunohistochemical localization of OCT2 in the cochlea of various species*. Laryngoscope 2015. [PMID 25892279](https://pubmed.ncbi.nlm.nih.gov/25892279/)
19. Ciarimboli G et al. *Organic cation transporter 2 mediates cisplatin-induced oto- and nephrotoxicity and is a target for protective interventions*. Am J Pathol 2010. [PMID 20110413](https://pubmed.ncbi.nlm.nih.gov/20110413/)
20. Lin LY et al. *Extracellular Ca(2+) and Mg(2+) modulate aminoglycoside blockade of mechanotransducer channel-mediated Ca(2+) entry in zebrafish hair cells: an in vivo study with the SIET*. Am J Physiol Cell Physiol 2013. [PMID 24005042](https://pubmed.ncbi.nlm.nih.gov/24005042/)
21. Alharazneh A et al. *Functional hair cell mechanotransducer channels are required for aminoglycoside ototoxicity*. PLoS One 2011. [PMID 21818312](https://pubmed.ncbi.nlm.nih.gov/21818312/)
22. Gale JE et al. *FM1-43 dye behaves as a permeant blocker of the hair-cell mechanotransducer channel*. J Neurosci 2001. [PMID 11549711](https://pubmed.ncbi.nlm.nih.gov/11549711/)

## C. 산화 스트레스와 항산화 예비능 (Oxidative stress and antioxidant reserve)

> 글루타티온 예비능 기울기(BGSH)와 포화 소모 항(KCON·KMG)의 근거.

23. Mukherjea D et al. *Transtympanic administration of short interfering (si)RNA for the NOX3 isoform of NADPH oxidase protects against cisplatin-induced hearing loss in the rat*. Antioxid Redox Signal 2010. [PMID 20214492](https://pubmed.ncbi.nlm.nih.gov/20214492/)
24. Mukherjea D et al. *NOX3 NADPH oxidase couples transient receptor potential vanilloid 1 to signal transducer and activator of transcription 1-mediated inflammation and hearing loss*. Antioxid Redox Signal 2011. [PMID 20712533](https://pubmed.ncbi.nlm.nih.gov/20712533/)
25. Rybak LP et al. *siRNA-mediated knock-down of NOX3: therapy for hearing loss?*. Cell Mol Life Sci 2012. [PMID 22562580](https://pubmed.ncbi.nlm.nih.gov/22562580/)
26. Sha SH et al. *Differential vulnerability of basal and apical hair cells is based on intrinsic susceptibility to free radicals*. Hear Res 2001. [PMID 11335071](https://pubmed.ncbi.nlm.nih.gov/11335071/)
27. Wang X et al. *Cisplatin-induced ototoxicity: From signaling network to therapeutic targets*. Biomed Pharmacother 2023. [PMID 36455457](https://pubmed.ncbi.nlm.nih.gov/36455457/)
28. Liu Z et al. *Inhibition of Gpx4-mediated ferroptosis alleviates cisplatin-induced hearing loss in C57BL/6 mice*. Mol Ther 2024. [PMID 38414247](https://pubmed.ncbi.nlm.nih.gov/38414247/)
29. Tan WJT et al. *Role of mitochondrial dysfunction and oxidative stress in sensorineural hearing loss*. Hear Res 2023. [PMID 37167889](https://pubmed.ncbi.nlm.nih.gov/37167889/)
30. Rybak LP et al. *Ototoxicity*. Kidney Int 2007. [PMID 17653135](https://pubmed.ncbi.nlm.nih.gov/17653135/)

## D. 세포사 경로 (Hair cell death pathways)

> 사멸 위험함수(KAP·OXC·HILL)와 염증 증폭항(GINFL), 이차 신경절 퇴행(KDEAFF)의 근거.

31. Ramkumar V et al. *Transient Receptor Potential Channels and Auditory Functions*. Antioxid Redox Signal 2022. [PMID 34465184](https://pubmed.ncbi.nlm.nih.gov/34465184/)
32. Schmitt NC et al. *Cisplatin-induced hair cell death requires STAT1 and is attenuated by epigallocatechin gallate*. J Neurosci 2009. [PMID 19321781](https://pubmed.ncbi.nlm.nih.gov/19321781/)
33. Levano S et al. *Loss of STAT1 protects hair cells from ototoxicity through modulation of STAT3, c-Jun, Akt, and autophagy factors*. Cell Death Dis 2015. [PMID 26673664](https://pubmed.ncbi.nlm.nih.gov/26673664/)
34. Kaur T et al. *Adenosine A1 Receptor Protects Against Cisplatin Ototoxicity by Suppressing the NOX3/STAT1 Inflammatory Pathway in the Cochlea*. J Neurosci 2016. [PMID 27053204](https://pubmed.ncbi.nlm.nih.gov/27053204/)
35. Ding D et al. *Cell death after co-administration of cisplatin and ethacrynic acid*. Hear Res 2007. [PMID 16978814](https://pubmed.ncbi.nlm.nih.gov/16978814/)
36. Ding D et al. *Review: ototoxic characteristics of platinum antitumor drugs*. Anat Rec (Hoboken) 2012. [PMID 23044998](https://pubmed.ncbi.nlm.nih.gov/23044998/)
37. Al Aameri RFH et al. *Role of RGS17 in cisplatin-induced cochlear inflammation and ototoxicity via caspase-3 activation*. Front Immunol 2025. [PMID 40061942](https://pubmed.ncbi.nlm.nih.gov/40061942/)
38. Kim J et al. *Alpha-lipoic acid protects against cisplatin-induced ototoxicity via the regulation of MAPKs and proinflammatory cytokines*. Biochem Biophys Res Commun 2014. [PMID 24796665](https://pubmed.ncbi.nlm.nih.gov/24796665/)
39. Cheng AG et al. *Mechanisms of hair cell death and protection*. Curr Opin Otolaryngol Head Neck Surg 2005. [PMID 16282762](https://pubmed.ncbi.nlm.nih.gov/16282762/)
40. Febles NK et al. *A combinatorial approach to protect sensory tissue against cisplatin-induced ototoxicity*. Hear Res 2022. [PMID 35051751](https://pubmed.ncbi.nlm.nih.gov/35051751/)
41. Jiang W et al. *BRCA1 deficiency promotes DNA damage in cochlear hair cells with activation of ATM-p53 pathway independent of CHK2*. Hear Res 2025. [PMID 40684720](https://pubmed.ncbi.nlm.nih.gov/40684720/)
42. Bedeir MM et al. *Multiplex immunohistochemistry reveals cochlear macrophage heterogeneity and local auditory nerve inflammation in cisplatin-induced hearing loss*. Front Neurol 2022. [PMID 36341090](https://pubmed.ncbi.nlm.nih.gov/36341090/)
43. Wood MB et al. *The Contribution of Immune Infiltrates to Ototoxicity and Cochlear Hair Cell Loss*. Front Cell Neurosci 2017. [PMID 28446866](https://pubmed.ncbi.nlm.nih.gov/28446866/)
44. Ito T et al. *Tissue-Resident Macrophages in the Stria Vascularis*. Front Neurol 2022. [PMID 35185769](https://pubmed.ncbi.nlm.nih.gov/35185769/)
45. Moser T et al. *Auditory neuropathy--neural and synaptic mechanisms*. Nat Rev Neurol 2016. [PMID 26891769](https://pubmed.ncbi.nlm.nih.gov/26891769/)

## E. 음조지형과 기저부 취약성 (Tonotopy and basal vulnerability)

> 8개 음조지형 밴드의 Greenwood 위치와 기저부 취약성 기울기의 근거.

46. Sridhar D et al. *A frequency-position function for the human cochlear spiral ganglion*. Audiol Neurootol 2006. [PMID 17063006](https://pubmed.ncbi.nlm.nih.gov/17063006/)
47. Greenwood DD et al. *A cochlear frequency-position function for several species--29 years later*. J Acoust Soc Am 1990. [PMID 2373794](https://pubmed.ncbi.nlm.nih.gov/2373794/)
48. Greenwood DD et al. *Comparing octaves, frequency ranges, and cochlear-map curvature across species*. Hear Res 1996. [PMID 8789821](https://pubmed.ncbi.nlm.nih.gov/8789821/)
49. Dalian D et al. *OTOTOXIC EFFECTS OF CARBOPLATIN IN ORGANOTYPIC CULTURES IN CHINCHILLAS AND RATS*. J Otol 2012. [PMID 25593588](https://pubmed.ncbi.nlm.nih.gov/25593588/)
50. Hibino H et al. *How is the highly positive endocochlear potential formed? The specific architecture of the stria vascularis and the roles of the ion-transport apparatus*. Pflugers Arch 2010. [PMID 20012478](https://pubmed.ncbi.nlm.nih.gov/20012478/)

## F. 티오황산나트륨 임상시험 (Sodium thiosulfate trials)

> S03·S04 및 Part 3 지연 스윕의 임상 대조군. SIOPEL-6와 ACCL0431의 상반된 결과가 모델의 핵심 검증 대상이다.

51. Freyer DR et al. *Effects of sodium thiosulfate versus observation on development of cisplatin-induced hearing loss in children with cancer (ACCL0431): a multicentre, randomised, controlled, open-label, phase 3 trial*. Lancet Oncol 2017. [PMID 27914822](https://pubmed.ncbi.nlm.nih.gov/27914822/)
52. Brock PR et al. *Sodium Thiosulfate for Protection from Cisplatin-Induced Hearing Loss*. N Engl J Med 2018. [PMID 29924955](https://pubmed.ncbi.nlm.nih.gov/29924955/)
53. Freyer DR et al. *Prevention of cisplatin-induced ototoxicity in children and adolescents with cancer: a clinical practice guideline*. Lancet Child Adolesc Health 2020. [PMID 31866182](https://pubmed.ncbi.nlm.nih.gov/31866182/)
54. Kallenberger EM et al. *Preventing Hearing Loss in Children Receiving Cisplatin: A Systematic Review and Meta-Analysis*. Laryngoscope 2025. [PMID 40165641](https://pubmed.ncbi.nlm.nih.gov/40165641/)
55. Freyer DR et al. *Interventions for cisplatin-induced hearing loss in children and adolescents with cancer*. Lancet Child Adolesc Health 2019. [PMID 31160205](https://pubmed.ncbi.nlm.nih.gov/31160205/)
56. Orgel E et al. *Reevaluation of sodium thiosulfate otoprotection using the consensus International Society of Paediatric Oncology Ototoxicity Scale: A report from the Children's Oncology Group study ACCL0431*. Pediatr Blood Cancer 2023. [PMID 37416942](https://pubmed.ncbi.nlm.nih.gov/37416942/)
57. Ohlsen TJD et al. *Otoprotective Effects of Sodium Thiosulfate by Demographic and Clinical Characteristics: A Report From Children's Oncology Group Study ACCL0431*. Pediatr Blood Cancer 2025. [PMID 39654065](https://pubmed.ncbi.nlm.nih.gov/39654065/)
58. Orgel E et al. *Assessment of provider perspectives on otoprotection research for children and adolescents: A Children's Oncology Group Cancer Control and Supportive Care Committee survey*. Pediatr Blood Cancer 2020. [PMID 32886425](https://pubmed.ncbi.nlm.nih.gov/32886425/)
59. Meijer AJM et al. *Use of Sodium Thiosulfate as an Otoprotectant in Patients With Cancer Treated With Platinum Compounds: A Review of the Literature*. J Clin Oncol 2024. [PMID 38648563](https://pubmed.ncbi.nlm.nih.gov/38648563/)
60. Farese S et al. *Sodium thiosulfate pharmacokinetics in hemodialysis patients and healthy volunteers*. Clin J Am Soc Nephrol 2011. [PMID 21566113](https://pubmed.ncbi.nlm.nih.gov/21566113/)
61. Goel R et al. *Effect of sodium thiosulfate on the pharmacokinetics and toxicity of cisplatin*. J Natl Cancer Inst 1989. [PMID 2552131](https://pubmed.ncbi.nlm.nih.gov/2552131/)
62. Sooriyaarachchi M et al. *Chemical basis for the detoxification of cisplatin-derived hydrolysis products by sodium thiosulfate*. J Inorg Biochem 2016. [PMID 27324827](https://pubmed.ncbi.nlm.nih.gov/27324827/)
63. Verschraagen M et al. *The chemical reactivity of BNP7787 and its metabolite mesna with the cytostatic agent cisplatin: comparison with the nucleophiles thiosulfate, DDTC, glutathione and its disulfide GSSG*. Cancer Chemother Pharmacol 2003. [PMID 12715205](https://pubmed.ncbi.nlm.nih.gov/12715205/)
64. Treskes M et al. *The chemical reactivity of the modulating agent WR2721 (ethiofos) and its main metabolites with the antitumor agents cisplatin and carboplatin*. Biochem Pharmacol 1991. [PMID 1659819](https://pubmed.ncbi.nlm.nih.gov/1659819/)

## G. 국소 투여 및 기타 보호제 (Local delivery and other otoprotectants)

> S05(고실내 투여)와 Part 8 경로 논증의 근거.

65. Dubashynskaya NV et al. *Advanced delivery systems for sodium thiosulfate in cisplatin otoprotection*. Eur J Pharm Biopharm 2026. [PMID 42437593](https://pubmed.ncbi.nlm.nih.gov/42437593/)
66. Videhult Pierre P et al. *Impact of pH on Intratympanic Sodium Thiosulfate-hyaluronan Gel in Preventing Cisplatin-induced Ototoxicity*. Otol Neurotol 2026. [PMID 41493913](https://pubmed.ncbi.nlm.nih.gov/41493913/)
67. Viglietta V et al. *Phase 1 study to evaluate safety, tolerability and pharmacokinetics of a novel intra-tympanic administered thiosulfate to prevent cisplatin-induced hearing loss in cancer patients*. Invest New Drugs 2020. [PMID 32157599](https://pubmed.ncbi.nlm.nih.gov/32157599/)
68. Tawalbeh M et al. *Intratympanic N-acetylcysteine in the prevention of cisplatin-induced ototoxicity: a systematic review and meta-analysis of randomized controlled trials*. BMC Pharmacol Toxicol 2025. [PMID 39905500](https://pubmed.ncbi.nlm.nih.gov/39905500/)
69. Rosas-Gutiérrez GDC et al. *[Efficacy of intratympanic infiltration of N-acetyl cysteine in cisplatin ototoxicity]*. Rev Med Inst Mex Seguro Soc 2023. [PMID 38016189](https://pubmed.ncbi.nlm.nih.gov/38016189/)
70. Gausterer JC et al. *Intratympanic application of poloxamer 407 hydrogels results in sustained N-acetylcysteine delivery to the inner ear*. Eur J Pharm Biopharm 2020. [PMID 32173603](https://pubmed.ncbi.nlm.nih.gov/32173603/)
71. Salt AN et al. *Principles of local drug delivery to the inner ear*. Audiol Neurootol 2009. [PMID 19923805](https://pubmed.ncbi.nlm.nih.gov/19923805/)
72. Borkholder DA et al. *Round window membrane intracochlear drug delivery enhanced by induced advection*. J Control Release 2014. [PMID 24291333](https://pubmed.ncbi.nlm.nih.gov/24291333/)
73. Hammer DR et al. *Novel dual-lumen microneedle delivers adeno-associated viral vectors in the guinea pig inner ear via the round window membrane*. Biomed Microdevices 2025. [PMID 40493265](https://pubmed.ncbi.nlm.nih.gov/40493265/)
74. Wang X et al. *A prestin-targeting peptide-guided drug delivery system rearranging concentration gradient in the inner ear: An improved strategy against hearing loss*. Eur J Pharm Sci 2023. [PMID 37295658](https://pubmed.ncbi.nlm.nih.gov/37295658/)
75. Fouladi M et al. *Amifostine protects against cisplatin-induced ototoxicity in children with average-risk medulloblastoma*. J Clin Oncol 2008. [PMID 18669462](https://pubmed.ncbi.nlm.nih.gov/18669462/)
76. Fisher MJ et al. *Amifostine for children with medulloblastoma treated with cisplatin-based chemotherapy*. Pediatr Blood Cancer 2004. [PMID 15390300](https://pubmed.ncbi.nlm.nih.gov/15390300/)
77. Araujo AGFS et al. *Cisplatin and ototoxicity in childhood: the perspective of supporting otoprotective agentes*. Braz J Biol 2024. [PMID 39140499](https://pubmed.ncbi.nlm.nih.gov/39140499/)
78. Ekborn A et al. *D-Methionine and cisplatin ototoxicity in the guinea pig: D-methionine influences cisplatin pharmacokinetics*. Hear Res 2002. [PMID 12031515](https://pubmed.ncbi.nlm.nih.gov/12031515/)
79. Campbell KC et al. *Oral D-methionine protects against cisplatin-induced hearing loss in humans: phase 2 randomized clinical trial in India*. Int J Audiol 2022. [PMID 34622731](https://pubmed.ncbi.nlm.nih.gov/34622731/)
80. Campbell KC et al. *D-methionine provides excellent protection from cisplatin ototoxicity in the rat*. Hear Res 1996. [PMID 8951454](https://pubmed.ncbi.nlm.nih.gov/8951454/)

## H. 임상 역학과 위험인자 (Clinical epidemiology and risk factors)

> 위험인자 파라미터(AGE·AMGLY·NOISE·FURO·GFR0)와 Part 9 분해의 근거.

81. Moke DJ et al. *Prevalence and risk factors for cisplatin-induced hearing loss in children, adolescents, and young adults: a multi-institutional North American cohort study*. Lancet Child Adolesc Health 2021. [PMID 33581749](https://pubmed.ncbi.nlm.nih.gov/33581749/)
82. Castelán-Martínez OD et al. *Hearing loss in Mexican children treated with cisplatin*. Int J Pediatr Otorhinolaryngol 2014. [PMID 25037447](https://pubmed.ncbi.nlm.nih.gov/25037447/)
83. Siemens A et al. *Role of Cisplatin Dose Intensity and TPMT Variation in the Development of Hearing Loss in Children*. Ther Drug Monit 2023. [PMID 36917731](https://pubmed.ncbi.nlm.nih.gov/36917731/)
84. Abu-Arja MH et al. *The cochlear dose and the age at radiotherapy predict severe hearing loss after passive scattering proton therapy and cisplatin in children with medulloblastoma*. Neuro Oncol 2024. [PMID 38916058](https://pubmed.ncbi.nlm.nih.gov/38916058/)
85. Brock PR et al. *Cisplatin ototoxicity in children: a practical grading system*. Med Pediatr Oncol 1991. [PMID 2056973](https://pubmed.ncbi.nlm.nih.gov/2056973/)
86. Yancey A et al. *Risk factors for cisplatin-associated ototoxicity in pediatric oncology patients*. Pediatr Blood Cancer 2012. [PMID 22431292](https://pubmed.ncbi.nlm.nih.gov/22431292/)
87. Turan C et al. *Cisplatin ototoxicity in children: risk factors and its relationship with polymorphisms of DNA repair genes ERCC1, ERCC2, and XRCC1*. Cancer Chemother Pharmacol 2019. [PMID 31586226](https://pubmed.ncbi.nlm.nih.gov/31586226/)
88. Fung C et al. *Toxicities Associated with Cisplatin-Based Chemotherapy and Radiotherapy in Long-Term Testicular Cancer Survivors*. Adv Urol 2018. [PMID 29670654](https://pubmed.ncbi.nlm.nih.gov/29670654/)
89. Sanchez VA et al. *Impact of cisplatin dose, renal function, and other factors on audiometrically-assessed ototoxicity in more than 1400 adult-onset cancer survivors from The Platinum Study: a multicentre cohort study*. EClinicalMedicine 2026. [PMID 41969333](https://pubmed.ncbi.nlm.nih.gov/41969333/)
90. Sprauten M et al. *Impact of long-term serum platinum concentrations on neuro- and ototoxicity in Cisplatin-treated survivors of testicular cancer*. J Clin Oncol 2012. [PMID 22184390](https://pubmed.ncbi.nlm.nih.gov/22184390/)
91. Sanchez VA et al. *Comprehensive Audiologic Analyses After Cisplatin-Based Chemotherapy*. JAMA Oncol 2024. [PMID 38842797](https://pubmed.ncbi.nlm.nih.gov/38842797/)
92. Skalleberg J et al. *The Relationship Between Cisplatin-related and Age-related Hearing Loss During an Extended Follow-up*. Laryngoscope 2020. [PMID 32065408](https://pubmed.ncbi.nlm.nih.gov/32065408/)
93. Fetoni AR et al. *Long-term auditory follow-up in the management of pediatric platinum-induced ototoxicity*. Eur Arch Otorhinolaryngol 2022. [PMID 35024956](https://pubmed.ncbi.nlm.nih.gov/35024956/)
94. Zsiros J et al. *Dose-dense cisplatin-based chemotherapy and surgery for children with high-risk hepatoblastoma (SIOPEL-4): a prospective, single-arm, feasibility study*. Lancet Oncol 2013. [PMID 23831416](https://pubmed.ncbi.nlm.nih.gov/23831416/)
95. Kohn S et al. *Ototoxicity resulting from combined administration of cisplatin and gentamicin*. Laryngoscope 1997. [PMID 9121325](https://pubmed.ncbi.nlm.nih.gov/9121325/)
96. Steyger PS et al. *Potentiation of Chemical Ototoxicity by Noise*. Semin Hear 2009. [PMID 20523755](https://pubmed.ncbi.nlm.nih.gov/20523755/)
97. Rybak LP et al. *Local Drug Delivery for Prevention of Hearing Loss*. Front Cell Neurosci 2019. [PMID 31338024](https://pubmed.ncbi.nlm.nih.gov/31338024/)
98. Chen Y et al. *Early Physiological and Cellular Indicators of Cisplatin-Induced Ototoxicity*. J Assoc Res Otolaryngol 2021. [PMID 33415542](https://pubmed.ncbi.nlm.nih.gov/33415542/)
99. Rybak LP et al. *Pathophysiology of furosemide ototoxicity*. J Otolaryngol 1982. [PMID 7042998](https://pubmed.ncbi.nlm.nih.gov/7042998/)
100. Rybak LP et al. *Ototoxicity of loop diuretics*. Otolaryngol Clin North Am 1993. [PMID 8233492](https://pubmed.ncbi.nlm.nih.gov/8233492/)
101. Rybak LP et al. *Comparative acute ototoxicity of loop diuretic compounds*. Eur Arch Otorhinolaryngol 1991. [PMID 1930985](https://pubmed.ncbi.nlm.nih.gov/1930985/)
102. Maideen NMP et al. *A Comprehensive Review of the Pharmacologic Perspective on Loop Diuretic Drug Interactions with Therapeutically Used Drugs*. Curr Drug Metab 2022. [PMID 35366769](https://pubmed.ncbi.nlm.nih.gov/35366769/)
103. Gehin W et al. *Sensorineural hearing loss after pediatric cranial radiotherapy: a multicenter analysis of dosimetric and clinical risk factors*. Support Care Cancer 2026. [PMID 42104091](https://pubmed.ncbi.nlm.nih.gov/42104091/)
104. Cohen-Cutler S et al. *Hearing Loss Risk in Pediatric Patients Treated with Cranial Irradiation and Cisplatin-Based Chemotherapy*. Int J Radiat Oncol Biol Phys 2021. [PMID 33677052](https://pubmed.ncbi.nlm.nih.gov/33677052/)

## I. 등급 체계와 청력 모니터링 (Grading scales and audiological monitoring)

> Brock·SIOP Boston·ASHA·CTCAE 등급을 모델의 출력으로 계산하기 위한 정의.

105. Landier W et al. *Ototoxicity in children with high-risk neuroblastoma: prevalence, risk factors, and concordance of grading scales--a report from the Children's Oncology Group*. J Clin Oncol 2014. [PMID 24419114](https://pubmed.ncbi.nlm.nih.gov/24419114/)
106. Clemens E et al. *A comparison of the Muenster, SIOP Boston, Brock, Chang and CTCAEv4.03 ototoxicity grading scales applied to 3,799 audiograms of childhood cancer patients treated with platinum-based chemotherapy*. PLoS One 2019. [PMID 30763334](https://pubmed.ncbi.nlm.nih.gov/30763334/)
107. Brock PR et al. *Platinum-induced ototoxicity in children: a consensus review on mechanisms, predisposition, and protection, including a new International Society of Pediatric Oncology Boston ototoxicity scale*. J Clin Oncol 2012. [PMID 22547603](https://pubmed.ncbi.nlm.nih.gov/22547603/)
108. Chawla A et al. *Incidence and Severity of Carboplatin-Associated Hearing Loss in Children With Cancer Assessed by the SIOP Boston 2012 Ototoxicity Criteria*. Pediatr Blood Cancer 2026. [PMID 42504125](https://pubmed.ncbi.nlm.nih.gov/42504125/)
109. Gertson K et al. *Prevalence of Ototoxicity Following Hematopoietic Stem Cell Transplantation in Pediatric Patients*. Biol Blood Marrow Transplant 2020. [PMID 31494228](https://pubmed.ncbi.nlm.nih.gov/31494228/)
110. Stevenson LJ et al. *Extended High-Frequency Audiometry for Ototoxicity Monitoring: A Longitudinal Evaluation of Drug-Resistant Tuberculosis Treatment*. Am J Audiol 2023. [PMID 36490390](https://pubmed.ncbi.nlm.nih.gov/36490390/)
111. Lough M et al. *Extended high-frequency audiometry in research and clinical practice*. J Acoust Soc Am 2022. [PMID 35364938](https://pubmed.ncbi.nlm.nih.gov/35364938/)
112. Fausti SA et al. *High-frequency monitoring for early detection of cisplatin ototoxicity*. Arch Otolaryngol Head Neck Surg 1993. [PMID 8499098](https://pubmed.ncbi.nlm.nih.gov/8499098/)
113. Bhagat SP et al. *Monitoring carboplatin ototoxicity with distortion-product otoacoustic emissions in children with retinoblastoma*. Int J Pediatr Otorhinolaryngol 2010. [PMID 20667604](https://pubmed.ncbi.nlm.nih.gov/20667604/)
114. Reavis KM et al. *Distortion-product otoacoustic emission test performance for ototoxicity monitoring*. Ear Hear 2011. [PMID 20625302](https://pubmed.ncbi.nlm.nih.gov/20625302/)
115. Konrad-Martin D et al. *Long-Term Variability of Distortion-Product Otoacoustic Emissions in Infants and Children and Its Relation to Pediatric Ototoxicity Monitoring*. Ear Hear 2020. [PMID 29280917](https://pubmed.ncbi.nlm.nih.gov/29280917/)
116. McMillan GP et al. *Accuracy of distortion-product otoacoustic emissions-based ototoxicity monitoring using various primary frequency step-sizes*. Int J Audiol 2012. [PMID 22676700](https://pubmed.ncbi.nlm.nih.gov/22676700/)

## J. 신독성과 전신 약동학 (Nephrotoxicity and systemic PK)

> 신독성 하위모델(KTUB·KLOSS·KPERM)과 GFR→청력 되먹임(Part 7)의 근거.

117. Fillastre JP et al. *Cisplatin nephrotoxicity*. Toxicol Lett 1989. [PMID 2650023](https://pubmed.ncbi.nlm.nih.gov/2650023/)
118. Reece PA et al. *Creatinine clearance as a predictor of ultrafilterable platinum disposition in cancer patients treated with cisplatin: relationship between peak ultrafilterable platinum plasma levels and nephrotoxicity*. J Clin Oncol 1987. [PMID 3806171](https://pubmed.ncbi.nlm.nih.gov/3806171/)
119. Bonetti A et al. *Cisplatin pharmacokinetics in elderly patients*. Ther Drug Monit 1994. [PMID 7846745](https://pubmed.ncbi.nlm.nih.gov/7846745/)
120. Dumas M et al. *Evaluation of the effect of furosemide on ultrafilterable platinum kinetics in patients treated with cis-diamminedichloroplatinum*. Cancer Chemother Pharmacol 1989. [PMID 2909288](https://pubmed.ncbi.nlm.nih.gov/2909288/)
121. Goel R et al. *Comparison of the pharmacokinetics of ultrafilterable cisplatin species detectable by derivatization with diethyldithiocarbamate or atomic absorption spectroscopy*. Eur J Cancer 1990. [PMID 2156545](https://pubmed.ncbi.nlm.nih.gov/2156545/)
122. Crona DJ et al. *A Systematic Review of Strategies to Prevent Cisplatin-Induced Nephrotoxicity*. Oncologist 2017. [PMID 28438887](https://pubmed.ncbi.nlm.nih.gov/28438887/)
123. Li J et al. *A systematic review for prevention of cisplatin-induced nephrotoxicity using different hydration protocols and meta-analysis for magnesium hydrate supplementation*. Clin Exp Nephrol 2024. [PMID 37530867](https://pubmed.ncbi.nlm.nih.gov/37530867/)
124. Suppadungsuk S et al. *Preloading magnesium attenuates cisplatin-associated nephrotoxicity: pilot randomized controlled trial (PRAGMATIC study)*. ESMO Open 2022. [PMID 34953401](https://pubmed.ncbi.nlm.nih.gov/34953401/)
125. Matsui M et al. *Magnesium supplementation therapy to prevent cisplatin-induced acute nephrotoxicity in pediatric cancer: a randomized phase-2 trial*. Int J Clin Oncol 2024. [PMID 38564107](https://pubmed.ncbi.nlm.nih.gov/38564107/)
126. McKeage MJ et al. *Comparative adverse effect profiles of platinum drugs*. Drug Saf 1995. [PMID 8573296](https://pubmed.ncbi.nlm.nih.gov/8573296/)
127. Schweitzer VG et al. *Cisplatin-induced ototoxicity: the effect of pigmentation and inhibitory agents*. Laryngoscope 1993. [PMID 8464301](https://pubmed.ncbi.nlm.nih.gov/8464301/)
128. Schell MJ et al. *Hearing loss in children and young adults receiving cisplatin with or without prior cranial irradiation*. J Clin Oncol 1989. [PMID 2715805](https://pubmed.ncbi.nlm.nih.gov/2715805/)
129. Osanto S et al. *Long-term effects of chemotherapy in patients with testicular cancer*. J Clin Oncol 1992. [PMID 1372350](https://pubmed.ncbi.nlm.nih.gov/1372350/)

## K. 약물유전체학 (Pharmacogenomics)

> 모델이 개인차로 남겨둔 부분 — 현재 결정론적 파라미터가 대표하는 집단 변이의 원천.

130. Thiesen S et al. *TPMT, COMT and ACYP2 genetic variants in paediatric cancer patients with cisplatin-induced ototoxicity*. Pharmacogenet Genomics 2017. [PMID 28445188](https://pubmed.ncbi.nlm.nih.gov/28445188/)
131. Hagleitner MM et al. *Influence of genetic variants in TPMT and COMT associated with cisplatin induced hearing loss in patients with cancer: two new cohorts and a meta-analysis reveal significant heterogeneity between cohorts*. PLoS One 2014. [PMID 25551397](https://pubmed.ncbi.nlm.nih.gov/25551397/)
132. Carleton BC et al. *Role of TPMT and COMT genetic variation in cisplatin-induced ototoxicity*. Clin Pharmacol Ther 2014. [PMID 24193170](https://pubmed.ncbi.nlm.nih.gov/24193170/)
133. Vos HI et al. *Replication of a genetic variant in ACYP2 associated with cisplatin-induced hearing loss in patients with osteosarcoma*. Pharmacogenet Genomics 2016. [PMID 26928270](https://pubmed.ncbi.nlm.nih.gov/26928270/)
134. Xu H et al. *Common variants in ACYP2 influence susceptibility to cisplatin-induced hearing loss*. Nat Genet 2015. [PMID 25665007](https://pubmed.ncbi.nlm.nih.gov/25665007/)
135. Tserga E et al. *The genetic vulnerability to cisplatin ototoxicity: a systematic review*. Sci Rep 2019. [PMID 30837596](https://pubmed.ncbi.nlm.nih.gov/30837596/)
136. Scott EN et al. *Systematic Critical Review of Genetic Factors Associated with Cisplatin-induced Ototoxicity: Canadian Pharmacogenomics Network for Drug Safety 2022 Update*. Ther Drug Monit 2023. [PMID 37726872](https://pubmed.ncbi.nlm.nih.gov/37726872/)
137. Oldenburg J et al. *Genetic variants associated with cisplatin-induced ototoxicity*. Pharmacogenomics 2008. [PMID 18855538](https://pubmed.ncbi.nlm.nih.gov/18855538/)

## L. 종양 효능과 백금-DNA 부가체 (Tumour efficacy and Pt-DNA adducts)

> 종양 효능 하위모델(KADF·KREPAIR·TUMPERF)과 Part 3의 효능 손실 계산의 근거.

138. Reed E et al. *Platinum-DNA adduct, nucleotide excision repair and platinum based anti-cancer chemotherapy*. Cancer Treat Rev 1998. [PMID 9861196](https://pubmed.ncbi.nlm.nih.gov/9861196/)
139. Olaussen KA et al. *PARP1 impact on DNA repair of platinum adducts: preclinical and clinical read-outs*. Lung Cancer 2013. [PMID 23410825](https://pubmed.ncbi.nlm.nih.gov/23410825/)
140. Galluzzi L et al. *Molecular mechanisms of cisplatin resistance*. Oncogene 2012. [PMID 21892204](https://pubmed.ncbi.nlm.nih.gov/21892204/)
141. Tan WJT et al. *Molecular Characteristics of Cisplatin-Induced Ototoxicity and Therapeutic Interventions*. Int J Mol Sci 2023. [PMID 38003734](https://pubmed.ncbi.nlm.nih.gov/38003734/)
142. Elmorsy EA et al. *Advances in understanding cisplatin-induced toxicity: Molecular mechanisms and protective strategies*. Eur J Pharm Sci 2024. [PMID 39423903](https://pubmed.ncbi.nlm.nih.gov/39423903/)
143. Trimmer EE et al. *Cisplatin*. Essays Biochem 1999. [PMID 10730196](https://pubmed.ncbi.nlm.nih.gov/10730196/)
144. Au JL et al. *Determinants of drug delivery and transport to solid tumors*. J Control Release 2001. [PMID 11489481](https://pubmed.ncbi.nlm.nih.gov/11489481/)
145. Salavati H et al. *Drug transport modeling in solid tumors: A computational exploration of spatial heterogeneity of biophysical properties*. Comput Biol Med 2023. [PMID 37392620](https://pubmed.ncbi.nlm.nih.gov/37392620/)
146. Kushner BH et al. *Ototoxicity from high-dose use of platinum compounds in patients with neuroblastoma*. Cancer 2006. [PMID 16779793](https://pubmed.ncbi.nlm.nih.gov/16779793/)
147. Jillella AP et al. *Ototoxicity after high-dose chemotherapy with cyclophosphamide, thiotepa and carboplatin followed by stem cell transplantation in patients with breast cancer*. Med Oncol 2000. [PMID 11114707](https://pubmed.ncbi.nlm.nih.gov/11114707/)

## M. 결과, 삶의 질, 임상 지침 (Outcomes, quality of life, guidelines)

> 임상 종말점(PTA·언어발달·삶의 질)과 등급-치료변경 연결의 근거.

148. Rajput K et al. *Ototoxicity-induced hearing loss and quality of life in survivors of paediatric cancer*. Int J Pediatr Otorhinolaryngol 2020. [PMID 33152988](https://pubmed.ncbi.nlm.nih.gov/33152988/)
149. Miaskowski C et al. *Associations among hearing loss, multiple co-occurring symptoms, and quality of life outcomes in cancer survivors*. J Cancer Surviv 2023. [PMID 36454519](https://pubmed.ncbi.nlm.nih.gov/36454519/)
150. Weiss A et al. *Hearing loss and quality of life in survivors of paediatric CNS tumours and other cancers*. Qual Life Res 2019. [PMID 30306534](https://pubmed.ncbi.nlm.nih.gov/30306534/)
151. Neuwelt EA et al. *Critical need for international consensus on ototoxicity assessment criteria*. J Clin Oncol 2010. [PMID 20194840](https://pubmed.ncbi.nlm.nih.gov/20194840/)

## N. QSP 방법론 (QSP methodology)

> 모델링 방법론.

152. Hardiansyah D et al. *Quantitative Systems Pharmacology Model of Chimeric Antigen Receptor T-Cell Therapy*. Clin Transl Sci 2019. [PMID 30990958](https://pubmed.ncbi.nlm.nih.gov/30990958/)
153. Bai JPF et al. *Creating a Roadmap to Quantitative Systems Pharmacology-Informed Rare Disease Drug Development: A Workshop Report*. Clin Pharmacol Ther 2024. [PMID 37984065](https://pubmed.ncbi.nlm.nih.gov/37984065/)
154. Zhu AZX et al. *Applications of Quantitative System Pharmacology Modeling to Model-Informed Drug Development*. Methods Mol Biol 2022. [PMID 35437719](https://pubmed.ncbi.nlm.nih.gov/35437719/)
155. Lesko LJ et al. *Perspective on model-informed drug development*. CPT Pharmacometrics Syst Pharmacol 2021. [PMID 34404115](https://pubmed.ncbi.nlm.nih.gov/34404115/)
156. Elmokadem A et al. *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial*. CPT Pharmacometrics Syst Pharmacol 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
157. Lu T et al. *gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve*. CPT Pharmacometrics Syst Pharmacol 2024. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
158. Hooijmaijers R et al. *Building an adaptive dose simulation framework to aid dose and schedule selection*. CPT Pharmacometrics Syst Pharmacol 2023. [PMID 37574587](https://pubmed.ncbi.nlm.nih.gov/37574587/)

---

## 문헌 총계

총 **158편**, 14개 주제 절. 모든 링크는 PubMed 원문 페이지로 연결됩니다.

## 모델이 문헌과 어긋나는 지점 (Disagreements, stated rather than tuned away)

1. **푸로세미드 병용** — 모델은 성인 600 mg/m² 배경에서 PTA +17.7 dB 를 예측한다.
   이는 고리이뇨제-백금 병용의 위험을 보고한 문헌(위 H절)의 방향과는 일치하지만 크기는
   사람 데이터로 검증된 적이 없다. 모델에서 가장 노출된 예측이다.
2. **분할 투여** — 같은 누적용량을 5일에 나누어도 모델은 PTA 0.50 dB 밖에 못 벌어준다.
   외림프가 혈장 첨두를 완전히 평활화하기 때문이며, 분할이 이독성을 의미 있게 줄인다는
   임상적 통념과 충돌한다. 모델은 그 통념이 참이라면 사멸함수의 Hill 계수가 1.7 보다
   훨씬 커야 함을 함의한다(Part 4 의 Hill 스윕).
3. **주간 저용량 요법(S20)** — 누적 360 mg/m² 인데도 모델의 GFR 이 60.7 mL/min 까지
   떨어진다(600 mg/m² 3주 간격은 79.9). 세뇨관 회복 반감기(KREPT)가 7일 간격보다 길기
   때문이며, 이 예측은 별도 검증이 필요하다.
4. **ACCL0431 해석** — 관류가 나쁜 전이성 종양에서 6시간 지연 티오황산염이 부가체
   형성의 3.7%를 잠식한다는 계산은 모델 수준의 가설이지 관측이 아니다. 다만 4시간으로
   당기면 14.5% 로 급증한다는 비선형성은 검증 가능한 형태의 예측이다.
