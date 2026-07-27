# Clostridioides difficile Infection (CDI) — References

Literature underpinning `cdi_qsp_model.dot` (mechanistic map) and
`cdi_mrgsolve_model.R` (61-compartment ODE model).

**129 references, every PMID resolved and title-verified against the PubMed
E-utilities API** (`esearch` + `esummary`) rather than transcribed from memory.
Sections follow the module structure of the mechanistic map.

---

## How these sources map onto the model

| Model element | Anchoring evidence |
|---|---|
| `BA_TCA → CD_GERM` (taurocholate germinant, CspC receptor) | Sorg & Sonenshein 2008; Francis 2013; Kochan 2017 |
| `BA_CDCA ⊣ CD_GERM` (competitive inhibition) | Sorg & Sonenshein 2010 |
| `MB_SBA → BA_DCA/BA_LCA` (bai 7α-dehydroxylation) | Ridlon 2006; Buffie 2015; Studer 2016; Kang 2019 |
| `BA_DCA/BA_LCA ⊣ CD_VEG` (outgrowth inhibition) | Thanissery 2017; Winston 2016; Łukawska 2022 |
| Antibiotic loss of secondary bile acids | Theriot 2014; Weingarden 2014 |
| `NUT_SIA → CD_VEG` (microbiota-liberated sialic acid) | Ng 2013; Almagro-Moreno 2009 |
| `NUT_AA` Stickland fermentation | Bouillaut 2013; Battaglioli 2018 |
| `SCFA_BUT` barrier support and virulence damping | Fachi 2019, 2020; Pensinger 2023 |
| `CD_SPORE_B` reservoir, `Spo0A` sporulation | Underwood 2009; Deakin 2012; Fimlaid 2013; Semenyuk 2014 |
| PaLoc control (`TcdR`, `TcdC`, `CodY`, `CcpA`, `SigD`) | Mani 2001; Dineen 2007; Dupuy 2008; Carter 2011; Antunes 2012; El Meouche 2013 |
| TcdB receptors (CSPG4, FZD, Nectin-3) | Yuan 2015; LaFrance 2015; Tao 2016 |
| Rho-GTPase glucosylation → barrier failure | Just 1995; Hecht 1988; Nusrat 2001 |
| `TCDB_MUC ⊣ EPI_SC` (FZD blockade of renewal) | Tao 2016; Mileto 2020 |
| Pyrin/NLRP3 → `IM_IL1B` | Ng 2010; Xu 2014; Hasegawa 2012 |
| `IM_IL22` protective limb | Hasegawa 2014; Frisbee 2019 |
| `AB_IGG` determines carriage vs recurrence | Kyne 2000, 2001 |
| Vancomycin vs fidaxomicin faecal PK, OP-1118 | Sears 2012; Babakhani 2011 |
| Fidaxomicin blocks sporulation and toxin output | Babakhani 2012; Louie 2012 |
| Fidaxomicin RNA-polymerase switch-region MoA | Boyaci 2018; Lin 2018 |
| Metronidazole faecal exposure falls as mucosa heals | Bolton 1986 |
| Bezlotoxumab PK and recurrence effect | Yee 2019; Wilcox 2017; Gerding 2018 |
| FMT / SER-109 / RBX2660 / VE303 engraftment | van Nood 2013; Kelly 2016; Feuerstadt 2022; Khanna 2022; Louie 2023 |
| Recurrence-rate calibration anchors | Louie 2011; Cornely 2012; Johnson 2014; Guery 2018 |
| Severity thresholds, fulminant disease, colectomy | McDonald 2018; Johnson 2021; Dallal 2002; Neal 2011 |

---

### 1. Epidemiology, burden, risk factors and natural history

1. Bartlett JG et al. *Antibiotic-associated pseudomembranous colitis due to toxin-producing clostridia.* N Engl J Med. 1978. [PMID 625309](https://pubmed.ncbi.nlm.nih.gov/625309/)
2. Bartlett JG et al. *Clinical practice. Antibiotic-associated diarrhea.* N Engl J Med. 2002. [PMID 11821511](https://pubmed.ncbi.nlm.nih.gov/11821511/)
3. Dallal RM et al. *Fulminant Clostridium difficile: an underappreciated and increasing cause of death and complications.* Ann Surg. 2002. [PMID 11882758](https://pubmed.ncbi.nlm.nih.gov/11882758/)
4. Loo VG et al. *A predominantly clonal multi-institutional outbreak of Clostridium difficile-associated diarrhea with high morbidity and mortality.* N Engl J Med. 2005. [PMID 16322602](https://pubmed.ncbi.nlm.nih.gov/16322602/)
5. McDonald LC et al. *An epidemic, toxin gene-variant strain of Clostridium difficile.* N Engl J Med. 2005. [PMID 16322603](https://pubmed.ncbi.nlm.nih.gov/16322603/)
6. Rupnik M et al. *Clostridium difficile infection: new developments in epidemiology and pathogenesis.* Nat Rev Microbiol. 2009. [PMID 19528959](https://pubmed.ncbi.nlm.nih.gov/19528959/)
7. Tleyjeh IM et al. *Association between proton pump inhibitor therapy and clostridium difficile infection: a contemporary systematic review and meta-analysis.* PLoS One. 2012. [PMID 23236397](https://pubmed.ncbi.nlm.nih.gov/23236397/)
8. Brown KA et al. *Meta-analysis of antibiotics and the risk of community-associated Clostridium difficile infection.* Antimicrob Agents Chemother. 2013. [PMID 23478961](https://pubmed.ncbi.nlm.nih.gov/23478961/)
9. He M et al. *Emergence and global spread of epidemic healthcare-associated Clostridium difficile.* Nat Genet. 2013. [PMID 23222960](https://pubmed.ncbi.nlm.nih.gov/23222960/)
10. Abou Chakra CN et al. *Risk factors for recurrence, complications and mortality in Clostridium difficile infection: a systematic review.* PLoS One. 2014. [PMID 24897375](https://pubmed.ncbi.nlm.nih.gov/24897375/)
11. Slimings C et al. *Antibiotics and hospital-acquired Clostridium difficile infection: update of systematic review and meta-analysis.* J Antimicrob Chemother. 2014. [PMID 24324224](https://pubmed.ncbi.nlm.nih.gov/24324224/)
12. Lessa FC et al. *Burden of Clostridium difficile infection in the United States.* N Engl J Med. 2015. [PMID 26061850](https://pubmed.ncbi.nlm.nih.gov/26061850/)
13. Rao K et al. *Clostridium difficile ribotype 027: relationship to age, detectability of toxins A or B in stool with rapid testing, severe infection, and mortality.* Clin Infect Dis. 2015. [PMID 25828993](https://pubmed.ncbi.nlm.nih.gov/25828993/)
14. Crobach MJT et al. *Understanding Clostridium difficile Colonization.* Clin Microbiol Rev. 2018. [PMID 29540433](https://pubmed.ncbi.nlm.nih.gov/29540433/)
15. Dieterle MG et al. *Systemic Inflammatory Mediators Are Effective Biomarkers for Predicting Adverse Outcomes in Clostridioides difficile Infection.* mBio. 2020. [PMID 32371595](https://pubmed.ncbi.nlm.nih.gov/32371595/)
16. Guh AY et al. *Trends in U.S. Burden of Clostridioides difficile Infection and Outcomes.* N Engl J Med. 2020. [PMID 32242357](https://pubmed.ncbi.nlm.nih.gov/32242357/)

### 2. Guidelines and consensus documents

17. McDonald LC et al. *Clinical Practice Guidelines for Clostridium difficile Infection in Adults and Children: 2017 Update by the Infectious Diseases Society of America (IDSA) and Society for Healthcare Epidemiology of America (SHEA).* Clin Infect Dis. 2018. [PMID 29562266](https://pubmed.ncbi.nlm.nih.gov/29562266/)
18. Sartelli M et al. *2019 update of the WSES guidelines for management of Clostridioides (Clostridium) difficile infection in surgical patients.* World J Emerg Surg. 2019. [PMID 30858872](https://pubmed.ncbi.nlm.nih.gov/30858872/)
19. Kelly CR et al. *ACG Clinical Guidelines: Prevention, Diagnosis, and Treatment of Clostridioides difficile Infections.* Am J Gastroenterol. 2021. [PMID 34003176](https://pubmed.ncbi.nlm.nih.gov/34003176/)
20. van Prehn J et al. *European Society of Clinical Microbiology and Infectious Diseases: 2021 update on the treatment guidance document for Clostridioides difficile infection in adults.* Clin Microbiol Infect. 2021. [PMID 34678515](https://pubmed.ncbi.nlm.nih.gov/34678515/)

### 3. Microbiota, colonization resistance and community recovery

21. Chang JY et al. *Decreased diversity of the fecal Microbiome in recurrent Clostridium difficile-associated diarrhea.* J Infect Dis. 2008. [PMID 18199029](https://pubmed.ncbi.nlm.nih.gov/18199029/)
22. Theriot CM et al. *Antibiotic-induced shifts in the mouse gut microbiome and metabolome increase susceptibility to Clostridium difficile infection.* Nat Commun. 2014. [PMID 24445449](https://pubmed.ncbi.nlm.nih.gov/24445449/)
23. Smith AB et al. *Enterococci enhance Clostridioides difficile pathogenesis.* Nature. 2022. [PMID 36385534](https://pubmed.ncbi.nlm.nih.gov/36385534/)

### 4. Bile acid metabolism and spore germination

24. Ridlon JM et al. *Bile salt biotransformations by human intestinal bacteria.* J Lipid Res. 2006. [PMID 16299351](https://pubmed.ncbi.nlm.nih.gov/16299351/)
25. Sorg JA et al. *Bile salts and glycine as cogerminants for Clostridium difficile spores.* J Bacteriol. 2008. [PMID 18245298](https://pubmed.ncbi.nlm.nih.gov/18245298/)
26. Sorg JA et al. *Inhibiting the initiation of Clostridium difficile spore germination using analogs of chenodeoxycholic acid, a bile acid.* J Bacteriol. 2010. [PMID 20675492](https://pubmed.ncbi.nlm.nih.gov/20675492/)
27. Francis MB et al. *Bile acid recognition by the Clostridium difficile germinant receptor, CspC, is important for establishing infection.* PLoS Pathog. 2013. [PMID 23675301](https://pubmed.ncbi.nlm.nih.gov/23675301/)
28. Weingarden AR et al. *Microbiota transplantation restores normal fecal bile acid composition in recurrent Clostridium difficile infection.* Am J Physiol Gastrointest Liver Physiol. 2014. [PMID 24284963](https://pubmed.ncbi.nlm.nih.gov/24284963/)
29. Buffie CG et al. *Precision microbiome reconstitution restores bile acid mediated resistance to Clostridium difficile.* Nature. 2015. [PMID 25337874](https://pubmed.ncbi.nlm.nih.gov/25337874/)
30. Studer N et al. *Functional Intestinal Bile Acid 7α-Dehydroxylation by Clostridium scindens Associated with Protection from Clostridium difficile Infection in a Gnotobiotic Mouse Model.* Front Cell Infect Microbiol. 2016. [PMID 28066726](https://pubmed.ncbi.nlm.nih.gov/28066726/)
31. Winston JA et al. *Impact of microbial derived secondary bile acids on colonization resistance against Clostridium difficile in the gastrointestinal tract.* Anaerobe. 2016. [PMID 27163871](https://pubmed.ncbi.nlm.nih.gov/27163871/)
32. Thanissery R et al. *Inhibition of spore germination, growth, and toxin activity of clinically relevant C. difficile strains by gut microbiota derived secondary bile acids.* Anaerobe. 2017. [PMID 28279860](https://pubmed.ncbi.nlm.nih.gov/28279860/)
33. Kang JD et al. *Bile Acid 7α-Dehydroxylating Gut Bacteria Secrete Antibiotics that Inhibit Clostridium difficile: Role of Secondary Bile Acids.* Cell Chem Biol. 2019. [PMID 30482679](https://pubmed.ncbi.nlm.nih.gov/30482679/)
34. Aguirre AM et al. *Bile acid-independent protection against Clostridioides difficile infection.* PLoS Pathog. 2021. [PMID 34665847](https://pubmed.ncbi.nlm.nih.gov/34665847/)
35. Łukawska A et al. *Impact of Primary and Secondary Bile Acids on Clostridioides difficile Infection.* Pol J Microbiol. 2022. [PMID 35635171](https://pubmed.ncbi.nlm.nih.gov/35635171/)

### 5. Nutrient niche, Stickland metabolism and short-chain fatty acids

36. Almagro-Moreno S et al. *Sialic acid catabolism confers a competitive advantage to pathogenic vibrio cholerae in the mouse intestine.* Infect Immun. 2009. [PMID 19564383](https://pubmed.ncbi.nlm.nih.gov/19564383/)
37. Bouillaut L et al. *Proline-dependent regulation of Clostridium difficile Stickland metabolism.* J Bacteriol. 2013. [PMID 23222730](https://pubmed.ncbi.nlm.nih.gov/23222730/)
38. Ng KM et al. *Microbiota-liberated host sugars facilitate post-antibiotic expansion of enteric pathogens.* Nature. 2013. [PMID 23995682](https://pubmed.ncbi.nlm.nih.gov/23995682/)
39. Battaglioli EJ et al. *Clostridioides difficile uses amino acids associated with gut microbial dysbiosis in a subset of patients with diarrhea.* Sci Transl Med. 2018. [PMID 30355801](https://pubmed.ncbi.nlm.nih.gov/30355801/)
40. Collins J et al. *Dietary trehalose enhances virulence of epidemic Clostridium difficile.* Nature. 2018. [PMID 29310122](https://pubmed.ncbi.nlm.nih.gov/29310122/)
41. Fachi JL et al. *Butyrate Protects Mice from Clostridium difficile-Induced Colitis through an HIF-1-Dependent Mechanism.* Cell Rep. 2019. [PMID 30995474](https://pubmed.ncbi.nlm.nih.gov/30995474/)
42. Fachi JL et al. *Acetate coordinates neutrophil and ILC3 responses against C. difficile through FFAR2.* J Exp Med. 2020. [PMID 31876919](https://pubmed.ncbi.nlm.nih.gov/31876919/)
43. Pensinger DA et al. *Butyrate Differentiates Permissiveness to Clostridioides difficile Infection and Influences Growth of Diverse C. difficile Isolates.* Infect Immun. 2023. [PMID 36692308](https://pubmed.ncbi.nlm.nih.gov/36692308/)

### 6. Life cycle: sporulation, biofilm and persistence

44. Underwood S et al. *Characterization of the sporulation initiation pathway of Clostridium difficile and its role in toxin production.* J Bacteriol. 2009. [PMID 19783633](https://pubmed.ncbi.nlm.nih.gov/19783633/)
45. Deakin LJ et al. *The Clostridium difficile spo0A gene is a persistence and transmission factor.* Infect Immun. 2012. [PMID 22615253](https://pubmed.ncbi.nlm.nih.gov/22615253/)
46. Fimlaid KA et al. *Global analysis of the sporulation pathway of Clostridium difficile.* PLoS Genet. 2013. [PMID 23950727](https://pubmed.ncbi.nlm.nih.gov/23950727/)
47. Edwards AN et al. *Conserved oligopeptide permeases modulate sporulation initiation in Clostridium difficile.* Infect Immun. 2014. [PMID 25069979](https://pubmed.ncbi.nlm.nih.gov/25069979/)
48. Semenyuk EG et al. *Spore formation and toxin production in Clostridium difficile biofilms.* PLoS One. 2014. [PMID 24498186](https://pubmed.ncbi.nlm.nih.gov/24498186/)

### 7. Adhesion, surface architecture and quorum sensing

49. Schwan C et al. *Clostridium difficile toxin CDT induces formation of microtubule-based protrusions and increases adherence of bacteria.* PLoS Pathog. 2009. [PMID 19834554](https://pubmed.ncbi.nlm.nih.gov/19834554/)
50. Ryan A et al. *A role for TLR4 in Clostridium difficile infection and the recognition of surface layer proteins.* PLoS Pathog. 2011. [PMID 21738466](https://pubmed.ncbi.nlm.nih.gov/21738466/)
51. Spigaglia P et al. *Surface layer protein A variant of Clostridium difficile PCR-ribotype 027.* Emerg Infect Dis. 2011. [PMID 21291621](https://pubmed.ncbi.nlm.nih.gov/21291621/)
52. Aubry A et al. *Modulation of toxin production by the flagellar regulon in Clostridium difficile.* Infect Immun. 2012. [PMID 22851750](https://pubmed.ncbi.nlm.nih.gov/22851750/)
53. Yoshino Y et al. *Clostridium difficile flagellin stimulates toll-like receptor 5, and toxin B promotes flagellin-induced chemokine production via TLR5.* Life Sci. 2013. [PMID 23261530](https://pubmed.ncbi.nlm.nih.gov/23261530/)
54. Darkoh C et al. *Toxin synthesis by Clostridium difficile is regulated through quorum signaling.* mBio. 2015. [PMID 25714717](https://pubmed.ncbi.nlm.nih.gov/25714717/)

### 8. PaLoc regulation of toxin expression

55. Mani N et al. *Regulation of toxin synthesis in Clostridium difficile by an alternative RNA polymerase sigma factor.* Proc Natl Acad Sci U S A. 2001. [PMID 11320220](https://pubmed.ncbi.nlm.nih.gov/11320220/)
56. Dineen SS et al. *Repression of Clostridium difficile toxin gene expression by CodY.* Mol Microbiol. 2007. [PMID 17725558](https://pubmed.ncbi.nlm.nih.gov/17725558/)
57. Dupuy B et al. *Clostridium difficile toxin synthesis is negatively regulated by TcdC.* J Med Microbiol. 2008. [PMID 18480323](https://pubmed.ncbi.nlm.nih.gov/18480323/)
58. Carter GP et al. *The anti-sigma factor TcdC modulates hypervirulence in an epidemic BI/NAP1/027 clinical isolate of Clostridium difficile.* PLoS Pathog. 2011. [PMID 22022270](https://pubmed.ncbi.nlm.nih.gov/22022270/)
59. Antunes A et al. *Global transcriptional control by glucose and carbon regulator CcpA in Clostridium difficile.* Nucleic Acids Res. 2012. [PMID 22989714](https://pubmed.ncbi.nlm.nih.gov/22989714/)
60. El Meouche I et al. *Characterization of the SigD regulon of C. difficile and its positive control of toxin production through the regulation of tcdR.* PLoS One. 2013. [PMID 24358307](https://pubmed.ncbi.nlm.nih.gov/24358307/)
61. Fletcher JR et al. *Shifts in the Gut Metabolome and Clostridium difficile Transcriptome throughout Colonization and Infection in a Mouse Model.* mSphere. 2018. [PMID 29600278](https://pubmed.ncbi.nlm.nih.gov/29600278/)

### 9. Toxin structure, receptors and intoxication mechanism

62. Hecht G et al. *Clostridium difficile toxin A perturbs cytoskeletal structure and tight junction permeability of cultured human intestinal epithelial monolayers.* J Clin Invest. 1988. [PMID 3141478](https://pubmed.ncbi.nlm.nih.gov/3141478/)
63. Just I et al. *Glucosylation of Rho proteins by Clostridium difficile toxin B.* Nature. 1995. [PMID 7777059](https://pubmed.ncbi.nlm.nih.gov/7777059/)
64. Nusrat A et al. *Clostridium difficile toxins disrupt epithelial barrier function by altering membrane microdomain localization of tight junction proteins.* Infect Immun. 2001. [PMID 11179295](https://pubmed.ncbi.nlm.nih.gov/11179295/)
65. Reineke J et al. *Autocatalytic cleavage of Clostridium difficile toxin B.* Nature. 2007. [PMID 17334356](https://pubmed.ncbi.nlm.nih.gov/17334356/)
66. Pruitt RN et al. *Structural organization of the functional domains of Clostridium difficile toxins A and B.* Proc Natl Acad Sci U S A. 2010. [PMID 20624955](https://pubmed.ncbi.nlm.nih.gov/20624955/)
67. Chumbler NM et al. *Clostridium difficile Toxin B causes epithelial cell necrosis through an autoprocessing-independent mechanism.* PLoS Pathog. 2012. [PMID 23236283](https://pubmed.ncbi.nlm.nih.gov/23236283/)
68. Farrow MA et al. *Clostridium difficile toxin B-induced necrosis is mediated by the host epithelial cell NADPH oxidase complex.* Proc Natl Acad Sci U S A. 2013. [PMID 24167244](https://pubmed.ncbi.nlm.nih.gov/24167244/)
69. Gerding DN et al. *Clostridium difficile binary toxin CDT: mechanism, epidemiology, and potential clinical importance.* Gut Microbes. 2014. [PMID 24253566](https://pubmed.ncbi.nlm.nih.gov/24253566/)
70. Xu H et al. *Innate immune sensing of bacterial modifications of Rho GTPases by the Pyrin inflammasome.* Nature. 2014. [PMID 24919149](https://pubmed.ncbi.nlm.nih.gov/24919149/)
71. LaFrance ME et al. *Identification of an epithelial cell receptor responsible for Clostridium difficile TcdB-induced cytotoxicity.* Proc Natl Acad Sci U S A. 2015. [PMID 26038560](https://pubmed.ncbi.nlm.nih.gov/26038560/)
72. Cowardin CA et al. *The binary toxin CDT enhances Clostridium difficile virulence by suppressing protective colonic eosinophilia.* Nat Microbiol. 2016. [PMID 27573114](https://pubmed.ncbi.nlm.nih.gov/27573114/)
73. Tao L et al. *Frizzled proteins are colonic epithelial receptors for C. difficile toxin B.* Nature. 2016. [PMID 27680706](https://pubmed.ncbi.nlm.nih.gov/27680706/)
74. Aktories K et al. *Clostridium difficile Toxin Biology.* Annu Rev Microbiol. 2017. [PMID 28657883](https://pubmed.ncbi.nlm.nih.gov/28657883/)
75. Chandrasekaran R et al. *The role of toxins in Clostridium difficile infection.* FEMS Microbiol Rev. 2017. [PMID 29048477](https://pubmed.ncbi.nlm.nih.gov/29048477/)
76. Mileto SJ et al. *Clostridioides difficile infection damages colonic stem cells via TcdB, impairing epithelial repair and recovery from disease.* Proc Natl Acad Sci U S A. 2020. [PMID 32198200](https://pubmed.ncbi.nlm.nih.gov/32198200/)

### 10. Host immune and inflammatory response

77. Kelly CP et al. *Neutrophil recruitment in Clostridium difficile toxin A enteritis in the rabbit.* J Clin Invest. 1994. [PMID 7907603](https://pubmed.ncbi.nlm.nih.gov/7907603/)
78. Pothoulakis C et al. *CP-96,345, a substance P antagonist, inhibits rat intestinal responses to Clostridium difficile toxin A but not cholera toxin.* Proc Natl Acad Sci U S A. 1994. [PMID 7508124](https://pubmed.ncbi.nlm.nih.gov/7508124/)
79. Sun X et al. *Essential role of the glucosyltransferase activity in Clostridium difficile toxin-induced secretion of TNF-alpha by macrophages.* Microb Pathog. 2009. [PMID 19324080](https://pubmed.ncbi.nlm.nih.gov/19324080/)
80. Ng J et al. *Clostridium difficile toxin-induced inflammation and intestinal injury are mediated by the inflammasome.* Gastroenterology. 2010. [PMID 20398664](https://pubmed.ncbi.nlm.nih.gov/20398664/)
81. Hasegawa M et al. *Protective role of commensals against Clostridium difficile infection via an IL-1β-mediated positive-feedback loop.* J Immunol. 2012. [PMID 22888139](https://pubmed.ncbi.nlm.nih.gov/22888139/)
82. Hasegawa M et al. *Interleukin-22 regulates the complement system to promote resistance against pathobionts after pathogen-induced intestinal damage.* Immunity. 2014. [PMID 25367575](https://pubmed.ncbi.nlm.nih.gov/25367575/)
83. Frisbee AL et al. *IL-33 drives group 2 innate lymphoid cell-mediated protection during Clostridium difficile infection.* Nat Commun. 2019. [PMID 31221971](https://pubmed.ncbi.nlm.nih.gov/31221971/)

### 11. Antibacterial therapy: randomised trials and comparative effectiveness

84. Louie TJ et al. *Tolevamer, a novel nonantibiotic polymer, compared with vancomycin in the treatment of mild to moderately severe Clostridium difficile-associated diarrhea.* Clin Infect Dis. 2006. [PMID 16838228](https://pubmed.ncbi.nlm.nih.gov/16838228/)
85. Musher DM et al. *Nitazoxanide for the treatment of Clostridium difficile colitis.* Clin Infect Dis. 2006. [PMID 16838229](https://pubmed.ncbi.nlm.nih.gov/16838229/)
86. Johnson S et al. *Interruption of recurrent Clostridium difficile-associated diarrhea episodes by serial therapy with vancomycin and rifaximin.* Clin Infect Dis. 2007. [PMID 17304459](https://pubmed.ncbi.nlm.nih.gov/17304459/)
87. Zar FA et al. *A comparison of vancomycin and metronidazole for the treatment of Clostridium difficile-associated diarrhea, stratified by disease severity.* Clin Infect Dis. 2007. [PMID 17599306](https://pubmed.ncbi.nlm.nih.gov/17599306/)
88. Garey KW et al. *A randomized, double-blind, placebo-controlled pilot study to assess the ability of rifaximin to prevent recurrent diarrhoea in patients with Clostridium difficile infection.* J Antimicrob Chemother. 2011. [PMID 21948965](https://pubmed.ncbi.nlm.nih.gov/21948965/)
89. Louie TJ et al. *Fidaxomicin versus vancomycin for Clostridium difficile infection.* N Engl J Med. 2011. [PMID 21288078](https://pubmed.ncbi.nlm.nih.gov/21288078/)
90. Cornely OA et al. *Fidaxomicin versus vancomycin for infection with Clostridium difficile in Europe, Canada, and the USA: a double-blind, non-inferiority, randomised controlled trial.* Lancet Infect Dis. 2012. [PMID 22321770](https://pubmed.ncbi.nlm.nih.gov/22321770/)
91. Johnson S et al. *Vancomycin, metronidazole, or tolevamer for Clostridium difficile infection: results from two multinational, randomized, controlled trials.* Clin Infect Dis. 2014. [PMID 24799326](https://pubmed.ncbi.nlm.nih.gov/24799326/)
92. Freeman J et al. *Pan-European longitudinal surveillance of antibiotic resistance among prevalent Clostridium difficile ribotypes.* Clin Microbiol Infect. 2015. [PMID 25701178](https://pubmed.ncbi.nlm.nih.gov/25701178/)
93. Stevens VW et al. *Comparative Effectiveness of Vancomycin and Metronidazole for the Prevention of Recurrence and Death in Patients With Clostridium difficile Infection.* JAMA Intern Med. 2017. [PMID 28166328](https://pubmed.ncbi.nlm.nih.gov/28166328/)
94. Vickers RJ et al. *Efficacy and safety of ridinilazole compared with vancomycin for the treatment of Clostridium difficile infection: a phase 2, randomised, double-blind, active-controlled, non-inferiority study.* Lancet Infect Dis. 2017. [PMID 28461207](https://pubmed.ncbi.nlm.nih.gov/28461207/)
95. Guery B et al. *Extended-pulsed fidaxomicin versus vancomycin for Clostridium difficile infection in patients 60 years and older (EXTEND): a randomised, controlled, open-label, phase 3b/4 trial.* Lancet Infect Dis. 2018. [PMID 29273269](https://pubmed.ncbi.nlm.nih.gov/29273269/)
96. Carlson TJ et al. *Ridinilazole for the treatment of Clostridioides difficile infection.* Expert Opin Investig Drugs. 2019. [PMID 30767587](https://pubmed.ncbi.nlm.nih.gov/30767587/)
97. Boekhoud IM et al. *Plasmid-mediated metronidazole resistance in Clostridioides difficile.* Nat Commun. 2020. [PMID 32001686](https://pubmed.ncbi.nlm.nih.gov/32001686/)

### 12. Antibacterial pharmacokinetics and pharmacodynamics

98. Bolton RP et al. *Faecal metronidazole concentrations during oral and intravenous therapy for antibiotic associated colitis due to Clostridium difficile.* Gut. 1986. [PMID 3781329](https://pubmed.ncbi.nlm.nih.gov/3781329/)
99. Babakhani F et al. *Postantibiotic effect of fidaxomicin and its major metabolite, OP-1118, against Clostridium difficile.* Antimicrob Agents Chemother. 2011. [PMID 21709084](https://pubmed.ncbi.nlm.nih.gov/21709084/)
100. Babakhani F et al. *Fidaxomicin inhibits spore production in Clostridium difficile.* Clin Infect Dis. 2012. [PMID 22752866](https://pubmed.ncbi.nlm.nih.gov/22752866/)
101. Louie TJ et al. *Fidaxomicin preserves the intestinal microbiome during and after treatment of Clostridium difficile infection (CDI) and reduces both toxin reexpression and recurrence of CDI.* Clin Infect Dis. 2012. [PMID 22752862](https://pubmed.ncbi.nlm.nih.gov/22752862/)
102. Sears P et al. *Fidaxomicin attains high fecal concentrations with minimal plasma concentrations following oral administration in patients with Clostridium difficile infection.* Clin Infect Dis. 2012. [PMID 22752859](https://pubmed.ncbi.nlm.nih.gov/22752859/)
103. Boyaci H et al. *Fidaxomicin jams Mycobacterium tuberculosis RNA polymerase motions needed for initiation via RbpA contacts.* Elife. 2018. [PMID 29480804](https://pubmed.ncbi.nlm.nih.gov/29480804/)
104. Lin W et al. *Structural Basis of Transcription Inhibition by Fidaxomicin (Lipiarmycin A3).* Mol Cell. 2018. [PMID 29606590](https://pubmed.ncbi.nlm.nih.gov/29606590/)
105. de Gunzburg J et al. *Protection of the Human Gut Microbiome From Antibiotics.* J Infect Dis. 2018. [PMID 29186529](https://pubmed.ncbi.nlm.nih.gov/29186529/)

### 13. Faecal microbiota transplantation and live biotherapeutics

106. McFarland LV et al. *A randomized placebo-controlled trial of Saccharomyces boulardii in combination with standard antibiotics for Clostridium difficile disease.* JAMA. 1994. [PMID 8201735](https://pubmed.ncbi.nlm.nih.gov/8201735/)
107. van Nood E et al. *Duodenal infusion of donor feces for recurrent Clostridium difficile.* N Engl J Med. 2013. [PMID 23323867](https://pubmed.ncbi.nlm.nih.gov/23323867/)
108. Seekatz AM et al. *Recovery of the gut microbiome following fecal microbiota transplantation.* mBio. 2014. [PMID 24939885](https://pubmed.ncbi.nlm.nih.gov/24939885/)
109. Cammarota G et al. *Randomised clinical trial: faecal microbiota transplantation by colonoscopy vs. vancomycin for the treatment of recurrent Clostridium difficile infection.* Aliment Pharmacol Ther. 2015. [PMID 25728808](https://pubmed.ncbi.nlm.nih.gov/25728808/)
110. Kelly CR et al. *Effect of Fecal Microbiota Transplantation on Recurrence in Multiply Recurrent Clostridium difficile Infection: A Randomized Trial.* Ann Intern Med. 2016. [PMID 27547925](https://pubmed.ncbi.nlm.nih.gov/27547925/)
111. Goldenberg JZ et al. *Probiotics for the prevention of Clostridium difficile-associated diarrhea in adults and children.* Cochrane Database Syst Rev. 2017. [PMID 29257353](https://pubmed.ncbi.nlm.nih.gov/29257353/)
112. Kao D et al. *Effect of Oral Capsule- vs Colonoscopy-Delivered Fecal Microbiota Transplantation on Recurrent Clostridium difficile Infection: A Randomized Clinical Trial.* JAMA. 2017. [PMID 29183074](https://pubmed.ncbi.nlm.nih.gov/29183074/)
113. Dubberke ER et al. *Results From a Randomized, Placebo-Controlled Clinical Trial of a RBX2660-A Microbiota-Based Drug for the Prevention of Recurrent Clostridium difficile Infection.* Clin Infect Dis. 2018. [PMID 29617739](https://pubmed.ncbi.nlm.nih.gov/29617739/)
114. DeFilipp Z et al. *Drug-Resistant E. coli Bacteremia Transmitted by Fecal Microbiota Transplant.* N Engl J Med. 2019. [PMID 31665575](https://pubmed.ncbi.nlm.nih.gov/31665575/)
115. Kelly CR et al. *Fecal Microbiota Transplantation Is Highly Effective in Real-World Practice: Initial Results From the FMT National Registry.* Gastroenterology. 2021. [PMID 33011173](https://pubmed.ncbi.nlm.nih.gov/33011173/)
116. Feuerstadt P et al. *SER-109, an Oral Microbiome Therapy for Recurrent Clostridioides difficile Infection.* N Engl J Med. 2022. [PMID 35045228](https://pubmed.ncbi.nlm.nih.gov/35045228/)
117. Khanna S et al. *Efficacy and Safety of RBX2660 in PUNCH CD3, a Phase III, Randomized, Double-Blind, Placebo-Controlled Trial with a Bayesian Primary Analysis for the Prevention of Recurrent Clostridioides difficile Infection.* Drugs. 2022. [PMID 36287379](https://pubmed.ncbi.nlm.nih.gov/36287379/)
118. Louie T et al. *VE303, a Defined Bacterial Consortium, for Prevention of Recurrent Clostridioides difficile Infection: A Randomized Clinical Trial.* JAMA. 2023. [PMID 37060545](https://pubmed.ncbi.nlm.nih.gov/37060545/)
119. Bhat A et al. *Safety and efficacy of fecal microbiota transplantation versus antibiotics for treating clostridioides difficile infection: systematic review and meta-analysis.* Eur J Clin Microbiol Infect Dis. 2026. [PMID 41081988](https://pubmed.ncbi.nlm.nih.gov/41081988/)

### 14. Toxin-directed immunotherapy, antibody and vaccines

120. Kyne L et al. *Asymptomatic carriage of Clostridium difficile and serum levels of IgG antibody against toxin A.* N Engl J Med. 2000. [PMID 10666429](https://pubmed.ncbi.nlm.nih.gov/10666429/)
121. Kotloff KL et al. *Safety and immunogenicity of increasing doses of a Clostridium difficile toxoid vaccine administered to healthy adults.* Infect Immun. 2001. [PMID 11159994](https://pubmed.ncbi.nlm.nih.gov/11159994/)
122. Kyne L et al. *Association between antibody response to toxin A and protection against recurrent Clostridium difficile diarrhoea.* Lancet. 2001. [PMID 11213096](https://pubmed.ncbi.nlm.nih.gov/11213096/)
123. Lowy I et al. *Treatment with monoclonal antibodies against Clostridium difficile toxins.* N Engl J Med. 2010. [PMID 20089970](https://pubmed.ncbi.nlm.nih.gov/20089970/)
124. Wilcox MH et al. *Bezlotoxumab for Prevention of Recurrent Clostridium difficile Infection.* N Engl J Med. 2017. [PMID 28121498](https://pubmed.ncbi.nlm.nih.gov/28121498/)
125. Gerding DN et al. *Bezlotoxumab for Prevention of Recurrent Clostridium difficile Infection in Patients at Increased Risk for Recurrence.* Clin Infect Dis. 2018. [PMID 29538686](https://pubmed.ncbi.nlm.nih.gov/29538686/)
126. Yee KL et al. *Population Pharmacokinetics and Pharmacodynamics of Bezlotoxumab in Adults with Primary and Recurrent Clostridium difficile Infection.* Antimicrob Agents Chemother. 2019. [PMID 30455246](https://pubmed.ncbi.nlm.nih.gov/30455246/)
127. de Bruyn G et al. *Safety, immunogenicity, and efficacy of a Clostridioides difficile toxoid vaccine candidate: a phase 3 multicentre, observer-blind, randomised, controlled trial.* Lancet Infect Dis. 2021. [PMID 32946836](https://pubmed.ncbi.nlm.nih.gov/32946836/)

### 15. Surgery and salvage in fulminant disease

128. Neal MD et al. *Diverting loop ileostomy and colonic lavage: an alternative to total abdominal colectomy for the treatment of severe, complicated Clostridium difficile associated disease.* Ann Surg. 2011. [PMID 21865943](https://pubmed.ncbi.nlm.nih.gov/21865943/)

### 16. Additional primary sources and reviews

129. Smits WK et al. *Clostridium difficile infection.* Nat Rev Dis Primers. 2016. [PMID 27158839](https://pubmed.ncbi.nlm.nih.gov/27158839/)

---

## Verification method

Each entry was resolved programmatically:

```
esearch(db=pubmed, term=<title>[ti])  →  esummary(db=pubmed, id=<candidates>)
```

and accepted only when the returned record's title matched the queried title
above a symmetric token-overlap threshold. Author, year, journal and title in
this file are the values PubMed returned, not values written from memory. No
PMID appears here that was not resolved this way.

## Notes and caveats

- **Population differences between the recurrence anchors.** The vancomycin,
  fidaxomicin, metronidazole and bezlotoxumab arms come from first-episode or
  mixed populations; the FMT and SER-109 trials enrolled recurrent-CDI
  populations with a higher baseline recurrence hazard. The model's
  RRI→probability calibration uses six anchors and notes this explicitly in
  `cdi_mrgsolve_model.R`; RBX2660 (PUNCH CD3) is deliberately excluded as an
  anchor because its 42.5% control-arm recurrence is not the same population as
  Louie 2011.
- **Metronidazole.** Its documented deficit against vancomycin in the pooled
  phase 3 analysis is chiefly in *cure* (72.7% vs 81.1%), not recurrence
  (~23% vs ~21%). The model reproduces the mechanism proposed for this — colonic
  drug exposure that depends on mucosal inflammation and therefore disappears as
  the patient improves — rather than assuming a fixed potency deficit.
- **Toxin concentration scale.** Faecal TcdA/TcdB are reported in this model in
  ng/mL of faecal supernatant. Published assays differ in units and in whether
  they report cytotoxin titre or immunoassay mass, so the absolute scale is
  order-of-magnitude rather than exact; the *relative* dynamics (toxin lagging
  the bacterial peak because of CodY/CcpA nutrient repression) are the
  calibrated feature.
- **Ribotype 027 parameters** (`RT027`, `RTTOXF`, `RTSPOR`) represent the
  aggregate phenotype reported for the BI/NAP1/027 lineage — truncated `tcdC`,
  higher toxin output, higher sporulation, fluoroquinolone resistance — not a
  single measured strain.

## Disclaimer

Educational and research QSP model. Semi-quantitative and not independently
validated. Not for clinical decision-making, prescribing, or regulatory
submission.
