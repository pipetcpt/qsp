# Head and Neck Squamous Cell Carcinoma — Reference List

Literature underpinning `hnscc_qsp_model.dot`, `hnscc_mrgsolve_model.R` and
`hnscc_shiny_app.R`. Every parameter value or structural choice in the model
that is not a free fitted quantity traces to one of the entries below.

References marked **[ANCHOR]** are the six clinical results used to fit the
five free tumour-control parameters. References marked **[HELD-OUT]** were
deliberately *not* used in fitting and are the tests the model had to pass
afterwards. That split is the whole basis for claiming the model has any
predictive content, so it is recorded here as well as in `README.md`.

---

## 1 · Epidemiology, aetiology and field cancerisation

1. Sung H, et al. Global Cancer Statistics 2020: GLOBOCAN estimates of incidence and mortality. *CA Cancer J Clin* 2021. <https://pubmed.ncbi.nlm.nih.gov/33538338/>
2. Chaturvedi AK, et al. Human papillomavirus and rising oropharyngeal cancer incidence in the United States. *J Clin Oncol* 2011. <https://pubmed.ncbi.nlm.nih.gov/21969503/>
3. Gillison ML, et al. Distinct risk factor profiles for human papillomavirus type 16-positive and -negative head and neck cancers. *J Natl Cancer Inst* 2008. <https://pubmed.ncbi.nlm.nih.gov/18334711/>
4. Hashibe M, et al. Interaction between tobacco and alcohol use and the risk of head and neck cancer (INHANCE consortium). *Cancer Epidemiol Biomarkers Prev* 2009. <https://pubmed.ncbi.nlm.nih.gov/19190158/>
5. Seitz HK, Stickel F. Molecular mechanisms of alcohol-mediated carcinogenesis. *Nat Rev Cancer* 2007. <https://pubmed.ncbi.nlm.nih.gov/17646865/>
6. Warnakulasuriya S, et al. Areca nut and oral cancer. *Addict Biol* 2002 / IARC Monograph 85. <https://pubmed.ncbi.nlm.nih.gov/12006218/>
7. Slaughter DP, Southwick HW, Smejkal W. Field cancerization in oral stratified squamous epithelium. *Cancer* 1953. <https://pubmed.ncbi.nlm.nih.gov/13094644/>
8. Braakhuis BJM, et al. A genetic explanation of Slaughter's concept of field cancerization. *Cancer Res* 2003. <https://pubmed.ncbi.nlm.nih.gov/12727832/>
9. Kutler DI, et al. High incidence of head and neck squamous cell carcinoma in patients with Fanconi anemia. *Arch Otolaryngol Head Neck Surg* 2003. <https://pubmed.ncbi.nlm.nih.gov/12874067/>
10. Leemans CR, Snijders PJF, Brakenhoff RH. The molecular landscape of head and neck cancer. *Nat Rev Cancer* 2018. <https://pubmed.ncbi.nlm.nih.gov/29497144/>

## 2 · HPV16 oncogenic program

11. Scheffner M, et al. The HPV-16 E6 and E6-AP complex functions as a ubiquitin-protein ligase in the ubiquitination of p53. *Cell* 1993. <https://pubmed.ncbi.nlm.nih.gov/8402913/>
12. Dyson N, et al. The human papilloma virus-16 E7 oncoprotein is able to bind to the retinoblastoma gene product. *Science* 1989. <https://pubmed.ncbi.nlm.nih.gov/2537532/>
13. Klaes R, et al. Overexpression of p16INK4A as a specific marker for dysplastic and neoplastic epithelial lesions of the cervix uteri. *Int J Cancer* 2001. <https://pubmed.ncbi.nlm.nih.gov/11340556/>
14. Lewis JS Jr, et al. p16 immunohistochemistry as a standalone test for risk stratification in oropharyngeal squamous cell carcinoma. *Head Neck Pathol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20386991/>
15. Park JW, et al. Deregulation of DNA damage response signalling by human papillomavirus E7. *J Virol* 2014 / Wallace NA, Galloway DA. *Annu Rev Virol* 2015. <https://pubmed.ncbi.nlm.nih.gov/26958917/>
16. Rieckmann T, et al. HNSCC cell lines positive for HPV and p16 possess higher cellular radiosensitivity due to an impaired DSB repair capacity. *Radiother Oncol* 2013. <https://pubmed.ncbi.nlm.nih.gov/23830196/>
17. Kimple RJ, et al. Enhanced radiation sensitivity in HPV-positive head and neck cancer. *Cancer Res* 2013. <https://pubmed.ncbi.nlm.nih.gov/23744740/>
18. Weaver AN, et al. DNA double strand break repair defect and sensitivity to poly ADP-ribose polymerase inhibition in HPV-positive HNSCC. *Oncotarget* 2015. <https://pubmed.ncbi.nlm.nih.gov/26356814/>
19. Klingenberg B, et al. E6/E7 and the tumour immune microenvironment in HPV-driven oropharyngeal cancer. *Br J Cancer* 2010 / Wansom D, et al. *Arch Otolaryngol* 2012. <https://pubmed.ncbi.nlm.nih.gov/22986739/>
20. Henderson S, et al. APOBEC-mediated cytosine deamination links PIK3CA helical domain mutations to human papillomavirus-driven tumour development. *Cell Rep* 2014. <https://pubmed.ncbi.nlm.nih.gov/25131207/>

## 3 · Genomics of HPV-negative disease

21. The Cancer Genome Atlas Network. Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature* 2015. <https://pubmed.ncbi.nlm.nih.gov/25631445/>
22. Stransky N, et al. The mutational landscape of head and neck squamous cell carcinoma. *Science* 2011. <https://pubmed.ncbi.nlm.nih.gov/21798893/>
23. Agrawal N, et al. Exome sequencing of head and neck squamous cell carcinoma reveals inactivating mutations in NOTCH1. *Science* 2011. <https://pubmed.ncbi.nlm.nih.gov/21798897/>
24. Poeta ML, et al. TP53 mutations and survival in squamous-cell carcinoma of the head and neck. *N Engl J Med* 2007. <https://pubmed.ncbi.nlm.nih.gov/18160686/>
25. Grandis JR, Tweardy DJ. Elevated levels of transforming growth factor alpha and epidermal growth factor receptor messenger RNA are early markers of carcinogenesis in head and neck cancer. *Cancer Res* 1993. <https://pubmed.ncbi.nlm.nih.gov/8339283/>
26. Chung CH, et al. Increased epidermal growth factor receptor gene copy number is associated with poor prognosis in head and neck squamous cell carcinomas. *J Clin Oncol* 2006. <https://pubmed.ncbi.nlm.nih.gov/16769987/>
27. Lui VWY, et al. Frequent mutation of the PI3K pathway in head and neck cancer defines predictive biomarkers. *Cancer Discov* 2013. <https://pubmed.ncbi.nlm.nih.gov/23619167/>
28. Martin D, et al. The head and neck cancer cell oncogenome: a platform for the development of precision molecular therapies. *Oncotarget* 2014. <https://pubmed.ncbi.nlm.nih.gov/25277197/>
29. Namani A, et al. NRF2-regulated metabolic gene signature in HNSCC / Kim YR, et al. Oncogenic NRF2 mutations in squamous cell carcinomas. *J Pathol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20521232/>
30. Morris LGT, et al. Recurrent somatic mutation of FAT1 in multiple human cancers leads to aberrant Wnt activation. *Nat Genet* 2013. <https://pubmed.ncbi.nlm.nih.gov/23354438/>

## 4 · EGFR signalling and its nuclear/repair role

31. Grandis JR, et al. Levels of TGF-alpha and EGFR protein in head and neck squamous cell carcinoma and patient survival. *J Natl Cancer Inst* 1998. <https://pubmed.ncbi.nlm.nih.gov/9625170/>
32. Dittmann K, Mayer C, Rodemann HP. Radiation-induced epidermal growth factor receptor nuclear import is linked to activation of DNA-dependent protein kinase. *J Biol Chem* 2005. <https://pubmed.ncbi.nlm.nih.gov/15879598/>
33. Dittmann K, et al. Radiation-induced caveolin-1 associated EGFR internalization is linked with nuclear EGFR transport and activation of DNA-PK. *Mol Cancer* 2008. <https://pubmed.ncbi.nlm.nih.gov/18644176/>
34. Huang SM, Harari PM. Modulation of radiation response after epidermal growth factor receptor blockade in squamous cell carcinomas. *Clin Cancer Res* 2000. <https://pubmed.ncbi.nlm.nih.gov/10741744/>
35. Milas L, et al. In vivo enhancement of tumor radioresponse by C225 antiepidermal growth factor receptor antibody. *Clin Cancer Res* 2000. <https://pubmed.ncbi.nlm.nih.gov/10741709/>
36. Bonner JA, et al. Enhanced apoptosis with combination C225/radiation treatment serves as the impetus for clinical investigation. *Semin Radiat Oncol* 2002. <https://pubmed.ncbi.nlm.nih.gov/12174342/>
37. Sen M, et al. Targeting STAT3 abrogates EGFR inhibitor resistance in cancer. *Clin Cancer Res* 2012. <https://pubmed.ncbi.nlm.nih.gov/22553343/>
38. Madoz-Gúrpide J, et al. Proteomics-based validation of genomic data / Knowles LM, et al. HGF and c-Met participate in paracrine tumorigenic pathways in head and neck squamous cell cancer. *Clin Cancer Res* 2009. <https://pubmed.ncbi.nlm.nih.gov/19351759/>

## 5 · Radiobiology — the linear-quadratic model and the five R's

39. Withers HR. The four R's of radiotherapy. *Adv Radiat Biol* 1975 (see also Steel GG, et al. *Int J Radiat Oncol Biol Phys* 1989 for the fifth). <https://pubmed.ncbi.nlm.nih.gov/2646263/>
40. Fowler JF. The linear-quadratic formula and progress in fractionated radiotherapy. *Br J Radiol* 1989. <https://pubmed.ncbi.nlm.nih.gov/2670032/>
41. Brenner DJ. The linear-quadratic model is an appropriate methodology for determining isoeffective doses at large doses per fraction. *Semin Radiat Oncol* 2008. <https://pubmed.ncbi.nlm.nih.gov/18725110/>
42. Bentzen SM, et al. Clinical radiobiology of squamous cell carcinoma of the oropharynx. *Int J Radiat Oncol Biol Phys* 1991. <https://pubmed.ncbi.nlm.nih.gov/2032890/>
43. Withers HR, Taylor JM, Maciejewski B. The hazard of accelerated tumor clonogen repopulation during radiotherapy. *Acta Oncol* 1988. <https://pubmed.ncbi.nlm.nih.gov/3390344/>
44. Bentzen SM, Thames HD. Clinical evidence for tumor clonogen regeneration: interpretations of the data. *Radiother Oncol* 1991. <https://pubmed.ncbi.nlm.nih.gov/1808721/>
45. Fowler JF, Lindstrom MJ. Loss of local control with prolongation in radiotherapy. *Int J Radiat Oncol Biol Phys* 1992. <https://pubmed.ncbi.nlm.nih.gov/1526880/>
46. **[HELD-OUT]** Hansen O, et al. Importance of overall treatment time for the outcome of radiotherapy of advanced head and neck carcinoma. *Radiother Oncol* 1997. <https://pubmed.ncbi.nlm.nih.gov/9106920/>
47. **[HELD-OUT]** Bese NS, Hendry J, Jeremic B. Effects of prolongation of overall treatment time due to unplanned interruptions during radiotherapy. *Int J Radiat Oncol Biol Phys* 2007. <https://pubmed.ncbi.nlm.nih.gov/17721981/>
48. Begg AC, et al. Cell proliferation as predictor of response to radiotherapy: Tpot measurement in head and neck cancer. *Radiother Oncol* 1992 / *Int J Radiat Oncol Biol Phys* 1990. <https://pubmed.ncbi.nlm.nih.gov/1631160/>
49. Baumann M, Krause M, Hill R. Exploring the role of cancer stem cells in radioresistance. *Nat Rev Cancer* 2008. <https://pubmed.ncbi.nlm.nih.gov/18784658/>
50. Prince ME, et al. Identification of a subpopulation of cells with cancer stem cell properties in head and neck squamous cell carcinoma. *Proc Natl Acad Sci USA* 2007. <https://pubmed.ncbi.nlm.nih.gov/17210912/>
51. **[HELD-OUT]** Horiot JC, et al. Hyperfractionation versus conventional fractionation in oropharyngeal carcinoma: EORTC 22791. *Radiother Oncol* 1992. <https://pubmed.ncbi.nlm.nih.gov/1281137/>
52. **[HELD-OUT]** Lacas B, et al. Meta-analysis of radiotherapy in carcinomas of head and neck (MARCH): updated individual-patient-data analysis. *Lancet Oncol* 2017. <https://pubmed.ncbi.nlm.nih.gov/28757375/>
53. Overgaard J, et al. Five compared with six fractions per week of conventional radiotherapy (DAHANCA 6 & 7). *Lancet* 2003. <https://pubmed.ncbi.nlm.nih.gov/14522594/>

## 6 · Hypoxia, reoxygenation and hypoxic radiosensitisers

54. Gray LH, et al. The concentration of oxygen dissolved in tissues at the time of irradiation as a factor in radiotherapy. *Br J Radiol* 1953. <https://pubmed.ncbi.nlm.nih.gov/13106296/>
55. Nordsmark M, et al. Prognostic value of tumor oxygenation in 397 head and neck tumors after primary radiation therapy: an international multi-center study. *Radiother Oncol* 2005. <https://pubmed.ncbi.nlm.nih.gov/16098619/>
56. Brizel DM, et al. Tumor oxygenation predicts for the likelihood of distant metastases in human soft tissue sarcoma / head and neck cancer. *Cancer Res* 1996. <https://pubmed.ncbi.nlm.nih.gov/8674032/>
57. **[ANCHOR]** Overgaard J, et al. A randomized double-blind phase III study of nimorazole as a hypoxic radiosensitizer of primary radiotherapy in supraglottic larynx and pharynx carcinoma (DAHANCA 5-85). *Radiother Oncol* 1998. <https://pubmed.ncbi.nlm.nih.gov/9635689/>
58. Overgaard J. Hypoxic modification of radiotherapy in squamous cell carcinoma of the head and neck — a systematic review and meta-analysis. *Radiother Oncol* 2011. <https://pubmed.ncbi.nlm.nih.gov/21684026/>
59. Toustrup K, et al. Development of a hypoxia gene expression classifier with predictive impact for hypoxic modification of radiotherapy in head and neck cancer. *Cancer Res* 2011. <https://pubmed.ncbi.nlm.nih.gov/21659506/>
60. Rischin D, et al. Tirapazamine, cisplatin, and radiation versus cisplatin and radiation for advanced squamous cell carcinoma of the head and neck (TROG 02.02, HeadSTART). *J Clin Oncol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20516441/>
61. Lee N, et al. Prospective trial incorporating pre-/mid-treatment [18F]-misonidazole PET for head-and-neck cancer. *Int J Radiat Oncol Biol Phys* 2009. <https://pubmed.ncbi.nlm.nih.gov/18929449/>
62. Semenza GL. Targeting HIF-1 for cancer therapy. *Nat Rev Cancer* 2003. <https://pubmed.ncbi.nlm.nih.gov/13130303/>
63. Hoff CM, et al. Importance of hemoglobin concentration and its modification for the outcome of head and neck cancer patients treated with radiotherapy. *Acta Oncol* 2012. <https://pubmed.ncbi.nlm.nih.gov/22150165/>
64. Browman GP, et al. Influence of cigarette smoking on the efficacy of radiation therapy in head and neck cancer. *N Engl J Med* 1993. <https://pubmed.ncbi.nlm.nih.gov/8426619/>

## 7 · Chemoradiotherapy — the platinum backbone

65. **[ANCHOR]** Pignon JP, et al. Meta-analysis of chemotherapy in head and neck cancer (MACH-NC): an update on 93 randomised trials and 17,346 patients. *Radiother Oncol* 2009. <https://pubmed.ncbi.nlm.nih.gov/19446902/>
66. Lacas B, et al. MACH-NC update: individual patient data network meta-analysis. *Lancet Oncol* 2021. <https://pubmed.ncbi.nlm.nih.gov/34049003/>
67. Adelstein DJ, et al. An intergroup phase III comparison of standard radiation therapy and two schedules of concurrent chemoradiotherapy in unresectable head and neck cancer. *J Clin Oncol* 2003. <https://pubmed.ncbi.nlm.nih.gov/12506176/>
68. Forastiere AA, et al. Concurrent chemotherapy and radiotherapy for organ preservation in advanced laryngeal cancer (RTOG 91-11). *N Engl J Med* 2003 / long-term *J Clin Oncol* 2013. <https://pubmed.ncbi.nlm.nih.gov/14645636/>
69. Cooper JS, et al. Postoperative concurrent radiotherapy and chemotherapy for high-risk squamous-cell carcinoma of the head and neck (RTOG 9501). *N Engl J Med* 2004. <https://pubmed.ncbi.nlm.nih.gov/15128893/>
70. Bernier J, et al. Postoperative irradiation with or without concomitant chemotherapy for locally advanced head and neck cancer (EORTC 22931). *N Engl J Med* 2004. <https://pubmed.ncbi.nlm.nih.gov/15128894/>
71. Bernier J, et al. Defining risk levels in locally advanced head and neck cancers: a comparative analysis of concurrent chemoradiation in the RTOG 9501 and EORTC 22931 trials. *Head Neck* 2005. <https://pubmed.ncbi.nlm.nih.gov/16161069/>
72. Ang KK, et al. Human papillomavirus and survival of patients with oropharyngeal cancer (RTOG 0129). *N Engl J Med* 2010. <https://pubmed.ncbi.nlm.nih.gov/20530316/>
73. Strojan P, et al. Cumulative cisplatin dose in concurrent chemoradiotherapy for head and neck cancer: a systematic review. *Head Neck* 2016. <https://pubmed.ncbi.nlm.nih.gov/26894744/>
74. **[HELD-OUT]** Noronha V, et al. Once-a-week versus once-every-3-weeks cisplatin chemoradiation for locally advanced head and neck cancer: a phase III randomized noninferiority trial. *J Clin Oncol* 2018. <https://pubmed.ncbi.nlm.nih.gov/29331199/>
75. **[HELD-OUT]** Kiyota N, et al. Weekly cisplatin plus radiation for postoperative head and neck cancer (JCOG1008): a multicenter, noninferiority, phase II/III trial. *J Clin Oncol* 2022. <https://pubmed.ncbi.nlm.nih.gov/35357907/>
76. Wilkins AC, et al. Cisplatin plus radiotherapy: mechanisms of interaction / Boeckman HJ, et al. Cisplatin sensitizes cancer cells to ionizing radiation via inhibition of nonhomologous end joining. *Mol Cancer Res* 2005. <https://pubmed.ncbi.nlm.nih.gov/15956241/>
77. Sears CR, Turchi JJ. Complex cisplatin-double strand break (DSB) lesions directly impair cellular non-homologous end-joining (NHEJ) independent of downstream damage response signaling. *J Biol Chem* 2012. <https://pubmed.ncbi.nlm.nih.gov/22992740/>
78. Olaussen KA, et al. DNA repair by ERCC1 in non-small-cell lung cancer and cisplatin-based adjuvant chemotherapy. *N Engl J Med* 2006 (ERCC1 principle applied to HNSCC in Handra-Luca A, et al. *Clin Cancer Res* 2007). <https://pubmed.ncbi.nlm.nih.gov/17671124/>

## 8 · Cetuximab and EGFR-directed de-escalation

79. **[ANCHOR]** Bonner JA, et al. Radiotherapy plus cetuximab for squamous-cell carcinoma of the head and neck. *N Engl J Med* 2006. <https://pubmed.ncbi.nlm.nih.gov/16467544/>
80. Bonner JA, et al. Radiotherapy plus cetuximab for locoregionally advanced head and neck cancer: 5-year survival data from a phase 3 randomised trial. *Lancet Oncol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20004617/>
81. **[ANCHOR]** Gillison ML, et al. Radiotherapy plus cetuximab or cisplatin in human papillomavirus-positive oropharyngeal cancer (NRG Oncology RTOG 1016): a randomised, multicentre, non-inferiority trial. *Lancet* 2019. <https://pubmed.ncbi.nlm.nih.gov/30449625/>
82. **[ANCHOR]** Mehanna H, et al. Radiotherapy plus cisplatin or cetuximab in low-risk human papillomavirus-positive oropharyngeal cancer (De-ESCALaTE HPV): an open-label randomised controlled phase 3 trial. *Lancet* 2019. <https://pubmed.ncbi.nlm.nih.gov/30449623/>
83. **[HELD-OUT]** Ang KK, et al. Randomized phase III trial of concurrent accelerated radiation plus cisplatin with or without cetuximab for stage III to IV head and neck carcinoma (RTOG 0522). *J Clin Oncol* 2014. <https://pubmed.ncbi.nlm.nih.gov/25154822/>
84. Vermorken JB, et al. Platinum-based chemotherapy plus cetuximab in head and neck cancer (EXTREME). *N Engl J Med* 2008. <https://pubmed.ncbi.nlm.nih.gov/18784101/>
85. Bibeau F, et al. Impact of FcγRIIa-FcγRIIIa polymorphisms and KRAS mutations on the clinical outcome of patients treated with cetuximab. *J Clin Oncol* 2009. <https://pubmed.ncbi.nlm.nih.gov/19204207/>
86. Taylor RJ, et al. FcγRIIIa polymorphisms and cetuximab-induced cytotoxicity in squamous cell carcinoma of the head and neck. *Cancer Immunol Immunother* 2009. <https://pubmed.ncbi.nlm.nih.gov/18987856/>
87. Dirix P, Vanstraelen B, Nuyts S. Cetuximab-induced skin rash and survival / Bonner JA, et al. Association of acneiform rash severity with outcome. *Ann Oncol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20130225/>
88. Schrag D, et al. Cetuximab therapy and symptomatic hypomagnesemia. *J Natl Cancer Inst* 2005 (mechanism: Groenestege WMT, et al. *J Clin Invest* 2007). <https://pubmed.ncbi.nlm.nih.gov/17671655/>
89. Chung CH, et al. Cetuximab-induced anaphylaxis and IgE specific for galactose-alpha-1,3-galactose. *N Engl J Med* 2008. <https://pubmed.ncbi.nlm.nih.gov/18334654/>
90. Dirix P, et al. Population pharmacokinetics of cetuximab / Fracasso PM, et al. A phase 1 escalating single-dose and weekly fixed-dose study of cetuximab. *Clin Cancer Res* 2007. <https://pubmed.ncbi.nlm.nih.gov/17332291/>

## 9 · Immune microenvironment and checkpoint blockade

91. Ferris RL, et al. Nivolumab for recurrent squamous-cell carcinoma of the head and neck (CheckMate 141). *N Engl J Med* 2016. <https://pubmed.ncbi.nlm.nih.gov/27718784/>
92. Burtness B, et al. Pembrolizumab alone or with chemotherapy versus cetuximab with chemotherapy for recurrent or metastatic squamous cell carcinoma of the head and neck (KEYNOTE-048). *Lancet* 2019. <https://pubmed.ncbi.nlm.nih.gov/31679945/>
93. Cohen EEW, et al. Pembrolizumab versus methotrexate, docetaxel, or cetuximab for recurrent head and neck squamous cell carcinoma (KEYNOTE-040). *Lancet* 2019. <https://pubmed.ncbi.nlm.nih.gov/30509740/>
94. Kulangara K, et al. Clinical utility of the combined positive score for PD-L1 expression and the approval of pembrolizumab for treatment of gastric cancer / PD-L1 CPS scoring in HNSCC. *Arch Pathol Lab Med* 2019. <https://pubmed.ncbi.nlm.nih.gov/30407860/>
95. Lee NY, et al. Avelumab plus standard-of-care chemoradiotherapy versus chemoradiotherapy alone in patients with locally advanced squamous cell carcinoma of the head and neck (JAVELIN Head and Neck 100). *Lancet Oncol* 2021. <https://pubmed.ncbi.nlm.nih.gov/33002432/>
96. Machiels JP, et al. Pembrolizumab plus concurrent chemoradiotherapy versus placebo plus chemoradiotherapy for locally advanced HNSCC (KEYNOTE-412). *Lancet Oncol* 2024. <https://pubmed.ncbi.nlm.nih.gov/38848741/>
97. Uppaluri R, et al. Neoadjuvant and adjuvant pembrolizumab in resectable locally advanced head and neck squamous cell carcinoma. *Clin Cancer Res* 2020. <https://pubmed.ncbi.nlm.nih.gov/32398345/>
98. Mandal R, et al. The head and neck cancer immune landscape and its immunotherapeutic implications. *JCI Insight* 2016. <https://pubmed.ncbi.nlm.nih.gov/27942583/>
99. Cillo AR, et al. Immune landscape of viral- and carcinogen-driven head and neck cancer. *Immunity* 2020. <https://pubmed.ncbi.nlm.nih.gov/31924475/>
100. Ferris RL, et al. Tumor antigen-targeted, monoclonal antibody-based immunotherapy: clinical response, cellular immunity and immunoescape. *J Clin Oncol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20940192/>
101. Leibowitz MS, et al. Deficiency of activated STAT1 in head and neck cancer cells mediates TAP1-dependent escape from cytotoxic T lymphocytes. *Cancer Immunol Immunother* 2011. <https://pubmed.ncbi.nlm.nih.gov/21519829/>
102. Ogino T, et al. HLA class I antigen down-regulation in primary laryngeal squamous cell carcinoma lesions as a poor prognostic marker. *Cancer Res* 2006. <https://pubmed.ncbi.nlm.nih.gov/16982766/>
103. Golden EB, Apetoh L. Radiotherapy and immunogenic cell death. *Semin Radiat Oncol* 2015. <https://pubmed.ncbi.nlm.nih.gov/25732004/>
104. Deng L, et al. STING-dependent cytosolic DNA sensing promotes radiation-induced type I interferon-dependent antitumor immunity. *Immunity* 2014. <https://pubmed.ncbi.nlm.nih.gov/25517616/>
105. Tang C, et al. Ipilimumab with stereotactic ablative radiation therapy / Twyman-Saint Victor C, et al. Radiation and dual checkpoint blockade activate non-redundant immune mechanisms. *Nature* 2015. <https://pubmed.ncbi.nlm.nih.gov/25754329/>
106. Ahn MJ, et al. Nivolumab plus ipilimumab in recurrent/metastatic HNSCC (CheckMate 651). *J Clin Oncol* 2023. <https://pubmed.ncbi.nlm.nih.gov/36867686/>
107. Postow MA, Sidlow R, Hellmann MD. Immune-related adverse events associated with immune checkpoint blockade. *N Engl J Med* 2018. <https://pubmed.ncbi.nlm.nih.gov/29320654/>
108. Wang PF, et al. Immune-related adverse events associated with anti-PD-1/PD-L1 treatment for malignancies: a meta-analysis. *Front Pharmacol* 2017. <https://pubmed.ncbi.nlm.nih.gov/29075197/>

## 10 · Induction chemotherapy, taxanes and 5-fluorouracil

109. Vermorken JB, et al. Cisplatin, fluorouracil, and docetaxel in unresectable head and neck cancer (TAX 323). *N Engl J Med* 2007. <https://pubmed.ncbi.nlm.nih.gov/17960013/>
110. Posner MR, et al. Cisplatin and fluorouracil alone or with docetaxel in head and neck cancer (TAX 324). *N Engl J Med* 2007. <https://pubmed.ncbi.nlm.nih.gov/17960012/>
111. Haddad R, et al. Induction chemotherapy followed by concurrent chemoradiotherapy versus concurrent chemoradiotherapy alone (PARADIGM). *Lancet Oncol* 2013. <https://pubmed.ncbi.nlm.nih.gov/23578174/>
112. Cohen EEW, et al. Phase III randomized trial of induction chemotherapy in patients with N2 or N3 locally advanced head and neck cancer (DeCIDE). *J Clin Oncol* 2014. <https://pubmed.ncbi.nlm.nih.gov/24958820/>
113. Longley DB, Harkin DP, Johnston PG. 5-Fluorouracil: mechanisms of action and clinical strategies. *Nat Rev Cancer* 2003. <https://pubmed.ncbi.nlm.nih.gov/12724731/>
114. Amstutz U, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) guideline for dihydropyrimidine dehydrogenase genotype and fluoropyrimidine dosing. *Clin Pharmacol Ther* 2018. <https://pubmed.ncbi.nlm.nih.gov/29152729/>
115. Bruno R, et al. Population pharmacokinetics/pharmacodynamics of docetaxel in phase II studies in patients with cancer. *J Clin Oncol* 1998. <https://pubmed.ncbi.nlm.nih.gov/9508171/>

## 11 · Platinum pharmacokinetics and organ toxicity

116. Urien S, Lokiec F. Population pharmacokinetics of total and unbound plasma cisplatin in adult patients. *Br J Clin Pharmacol* 2004. <https://pubmed.ncbi.nlm.nih.gov/15606443/>
117. Vermorken JB, et al. Pharmacokinetics of free and total platinum species after short-term infusion of cisplatin. *Cancer Treat Rep* 1984 / Erdlenbruch B, et al. *Br J Clin Pharmacol* 2001. <https://pubmed.ncbi.nlm.nih.gov/11736867/>
118. Calvert AH, et al. Carboplatin dosage: prospective evaluation of a simple formula based on renal function. *J Clin Oncol* 1989. <https://pubmed.ncbi.nlm.nih.gov/2681557/>
119. Ciarimboli G, et al. Organic cation transporter 2 mediates cisplatin-induced oto- and nephrotoxicity. *Am J Pathol* 2010. <https://pubmed.ncbi.nlm.nih.gov/20431035/>
120. Miller RP, et al. Mechanisms of cisplatin nephrotoxicity. *Toxins* 2010. <https://pubmed.ncbi.nlm.nih.gov/22069601/>
121. Rybak LP, et al. Mechanisms of cisplatin-induced ototoxicity and prevention. *Hear Res* 2007. <https://pubmed.ncbi.nlm.nih.gov/17706109/>
122. Breglio AM, et al. Cisplatin is retained in the cochlea indefinitely following chemotherapy. *Nat Commun* 2017. <https://pubmed.ncbi.nlm.nih.gov/29162831/>
123. Ross CJD, et al. Genetic variants in TPMT and COMT are associated with hearing loss in children receiving cisplatin chemotherapy. *Nat Genet* 2009 / Xu H, et al. ACYP2. *Nat Genet* 2015. <https://pubmed.ncbi.nlm.nih.gov/26301491/>
124. Friberg LE, et al. Model of chemotherapy-induced myelosuppression with parameter consistency across drugs. *J Clin Oncol* 2002. <https://pubmed.ncbi.nlm.nih.gov/12488408/>
125. Frame D. Best practice management of CINV in oncology patients / Hesketh PJ. Chemotherapy-induced nausea and vomiting. *N Engl J Med* 2008. <https://pubmed.ncbi.nlm.nih.gov/18525043/>

## 12 · Normal-tissue injury: mucositis, xerostomia, dysphagia, late effects

126. Sonis ST. The pathobiology of mucositis. *Nat Rev Cancer* 2004. <https://pubmed.ncbi.nlm.nih.gov/15057287/>
127. Sonis ST. Pathobiology of oral mucositis: novel insights and opportunities. *J Support Oncol* 2007. <https://pubmed.ncbi.nlm.nih.gov/17958279/>
128. Trotti A, et al. Mucositis incidence, severity and associated outcomes in patients with head and neck cancer receiving radiotherapy with or without chemotherapy: a systematic literature review. *Radiother Oncol* 2003. <https://pubmed.ncbi.nlm.nih.gov/12873674/>
129. Elting LS, et al. Risk, outcomes, and costs of radiation-induced oral mucositis among patients with head-and-neck malignancies. *Int J Radiat Oncol Biol Phys* 2007. <https://pubmed.ncbi.nlm.nih.gov/17544601/>
130. Dörr W, Hendry JH. Consequential late effects in normal tissues. *Radiother Oncol* 2001. <https://pubmed.ncbi.nlm.nih.gov/11395258/>
131. Eisbruch A, et al. Dose, volume, and function relationships in parotid salivary glands following conformal and intensity-modulated irradiation of head and neck cancer. *Int J Radiat Oncol Biol Phys* 1999. <https://pubmed.ncbi.nlm.nih.gov/10487552/>
132. Deasy JO, et al. Radiotherapy dose-volume effects on salivary gland function (QUANTEC). *Int J Radiat Oncol Biol Phys* 2010. <https://pubmed.ncbi.nlm.nih.gov/20171519/>
133. **[HELD-OUT]** Nutting CM, et al. Parotid-sparing intensity modulated versus conventional radiotherapy in head and neck cancer (PARSPORT): a phase 3 multicentre randomised controlled trial. *Lancet Oncol* 2011. <https://pubmed.ncbi.nlm.nih.gov/21236730/>
134. Konings AWT, Coppes RP, Vissink A. On the mechanism of salivary gland radiosensitivity. *Int J Radiat Oncol Biol Phys* 2005. <https://pubmed.ncbi.nlm.nih.gov/16226396/>
135. Eisbruch A, et al. Dysphagia and aspiration after chemoradiotherapy for head-and-neck cancer: which anatomic structures are affected and can they be spared by IMRT? *Int J Radiat Oncol Biol Phys* 2004. <https://pubmed.ncbi.nlm.nih.gov/15519788/>
136. Hutcheson KA, et al. Dynamic Imaging Grade of Swallowing Toxicity (DIGEST): scale development and validation. *Cancer* 2017. <https://pubmed.ncbi.nlm.nih.gov/27808424/>
137. Straub JM, et al. Radiation-induced fibrosis: mechanisms and implications for therapy. *J Cancer Res Clin Oncol* 2015. <https://pubmed.ncbi.nlm.nih.gov/25573644/>
138. Nabil S, Samman N. Incidence and prevention of osteoradionecrosis after dental extraction in irradiated patients: a systematic review. *Int J Oral Maxillofac Surg* 2011. <https://pubmed.ncbi.nlm.nih.gov/21396799/>
139. Boomsma MJ, Bijl HP, Langendijk JA. Radiation-induced hypothyroidism in head and neck cancer patients: a systematic review. *Radiother Oncol* 2011. <https://pubmed.ncbi.nlm.nih.gov/21255848/>
140. Dorresteijn LDA, et al. Increased risk of ischemic stroke after radiotherapy on the neck in patients younger than 60 years. *J Clin Oncol* 2002. <https://pubmed.ncbi.nlm.nih.gov/11844823/>
141. Bjordal K, et al. A 12 country field study of the EORTC QLQ-C30 and the head and neck cancer specific module (EORTC QLQ-H&N35). *Eur J Cancer* 2000. <https://pubmed.ncbi.nlm.nih.gov/10974628/>

## 13 · Host state, nutrition and cachexia

142. Fearon K, et al. Definition and classification of cancer cachexia: an international consensus. *Lancet Oncol* 2011. <https://pubmed.ncbi.nlm.nih.gov/21296615/>
143. Baracos VE, et al. Cancer-associated cachexia. *Nat Rev Dis Primers* 2018. <https://pubmed.ncbi.nlm.nih.gov/29345251/>
144. Johnen H, et al. Tumor-induced anorexia and weight loss are mediated by the TGF-beta superfamily cytokine MIC-1 (GDF-15). *Nat Med* 2007. <https://pubmed.ncbi.nlm.nih.gov/17982465/>
145. Langius JAE, et al. Critical weight loss is a major prognostic indicator for disease-specific survival in patients with head and neck cancer receiving radiotherapy. *Br J Cancer* 2013. <https://pubmed.ncbi.nlm.nih.gov/23695023/>
146. Wendrich AW, et al. Low skeletal muscle mass is a predictive factor for chemotherapy dose-limiting toxicity in patients with locally advanced head and neck cancer. *Oral Oncol* 2017. <https://pubmed.ncbi.nlm.nih.gov/28688688/>
147. Piccirillo JF, et al. Development of a new head and neck cancer-specific comorbidity index (ACE-27). *Arch Otolaryngol Head Neck Surg* 2002. <https://pubmed.ncbi.nlm.nih.gov/12365042/>

## 14 · Staging, biomarkers and response assessment

148. Lydiatt WM, et al. Head and Neck cancers — major changes in the American Joint Committee on Cancer eighth edition cancer staging manual. *CA Cancer J Clin* 2017. <https://pubmed.ncbi.nlm.nih.gov/28128848/>
149. O'Sullivan B, et al. Development and validation of a staging system for HPV-related oropharyngeal cancer by the ICON-S group. *Lancet Oncol* 2016. <https://pubmed.ncbi.nlm.nih.gov/26670617/>
150. Eisenhauer EA, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). *Eur J Cancer* 2009. <https://pubmed.ncbi.nlm.nih.gov/19097774/>
151. Chera BS, et al. Plasma circulating tumor HPV DNA for the surveillance of cancer recurrence in HPV-associated oropharyngeal cancer. *J Clin Oncol* 2020. <https://pubmed.ncbi.nlm.nih.gov/32348703/>
152. Chera BS, et al. Rapid clearance profile of plasma circulating tumor HPV type 16 DNA during chemoradiotherapy correlates with disease control. *Clin Cancer Res* 2019. <https://pubmed.ncbi.nlm.nih.gov/31088830/>
153. Ferris RL, et al. Extranodal extension and its prognostic significance / Amin MB. AJCC Cancer Staging Manual 8th ed. *Head Neck* 2019. <https://pubmed.ncbi.nlm.nih.gov/31179584/>

## 15 · QSP / modelling methodology

154. Fowler JF. Rapid repopulation in radiotherapy: a debate on mechanism. *Radiother Oncol* 1991 / Fowler JF. Biological factors influencing optimum fractionation. *Acta Oncol* 2001. <https://pubmed.ncbi.nlm.nih.gov/11441937/>
155. Webb S, Nahum AE. A model for calculating tumour control probability in radiotherapy including the effects of inhomogeneous distributions of dose and clonogenic cell density. *Phys Med Biol* 1993. <https://pubmed.ncbi.nlm.nih.gov/8346284/>
156. Zaider M, Minerbo GN. Tumour control probability: a formulation applicable to any temporal protocol of dose delivery. *Phys Med Biol* 2000. <https://pubmed.ncbi.nlm.nih.gov/10701509/>
157. Bentzen SM, Tucker SL. Quantifying the position and steepness of radiation dose-response curves. *Int J Radiat Biol* 1997. <https://pubmed.ncbi.nlm.nih.gov/9246191/>
158. Curtis SB. Lethal and potentially lethal lesions induced by radiation — a unified repair model. *Radiat Res* 1986. <https://pubmed.ncbi.nlm.nih.gov/3737748/>
159. Brenner DJ, et al. A convenient extension of the linear-quadratic model to include redistribution and reoxygenation. *Int J Radiat Oncol Biol Phys* 1995. <https://pubmed.ncbi.nlm.nih.gov/7635772/>
160. Baker RE, Peña JM, Jayamohan J, Jérusalem A. Mechanistic models versus machine learning, a fight worth fighting for the biological community? *Biol Lett* 2018. <https://pubmed.ncbi.nlm.nih.gov/29769297/>
161. Elmokadem A, et al. Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial. *CPT Pharmacometrics Syst Pharmacol* 2019. <https://pubmed.ncbi.nlm.nih.gov/31654579/>
162. Baron KT, et al. mrgsolve: Simulate from ODE-Based Models. R package. <https://mrgsolve.org/>
163. Chelliah V, et al. Quantitative systems pharmacology approaches for immuno-oncology: adding virtual patients to the development paradigm. *Clin Pharmacol Ther* 2021. <https://pubmed.ncbi.nlm.nih.gov/32433790/>
164. Jafarnejad M, et al. A computational model of neoadjuvant PD-1 inhibition in non-small cell lung cancer. *AAPS J* 2019. <https://pubmed.ncbi.nlm.nih.gov/30984968/>
165. Wang H, et al. In silico simulation of a clinical trial with anti-CTLA-4 and anti-PD-L1 immunotherapies in metastatic breast cancer using a systems pharmacology model. *R Soc Open Sci* 2019. <https://pubmed.ncbi.nlm.nih.gov/31417737/>

---

## How the anchors and held-out tests are used

| Role | Trial / source | Quantity | Used for |
|------|----------------|----------|----------|
| **[ANCHOR]** | Bonner 2006 / 2010 (79, 80) | 3-yr LRC 34 % (RT) and 47 % (cetuximab-RT), HPV-unselected | fits `ALPHA0`, `FCSC0` |
| **[ANCHOR]** | MACH-NC (65, 66), RTOG 0129 (72) | LRC ≈ 58 % for cisplatin-RT, HPV-negative | fits `PHIPT` |
| **[ANCHOR]** | RTOG 1016 (81), De-ESCALaTE (82) | 5-yr LRC 90 % (cisplatin) vs 83 % (cetuximab), HPV-positive | fits `KHRA` |
| **[ANCHOR]** | DAHANCA 5-85 (57) | 5-yr LRC 33 % → 49 % with nimorazole | fits `KOXTR` |
| **[HELD-OUT]** | RTOG 0522 (83) | cetuximab added to cisplatin-RT: no PFS/OS benefit | **model over-predicts** — see `README.md` |
| **[HELD-OUT]** | Hansen 1997 (46), Bese 2007 (47), Withers 1988 (43) | ≈ 1 % LRC lost per day of prolongation; 0.6-0.9 Gy/day | reproduced, and predicted to be position-dependent |
| **[HELD-OUT]** | EORTC 22791 (51), MARCH (52) | hyperfractionation improves LRC | reproduced |
| **[HELD-OUT]** | Noronha 2018 (74), JCOG1008 (75) | weekly vs 3-weekly cisplatin roughly comparable | direction reproduced, magnitude wrong |
| **[HELD-OUT]** | PARSPORT (133), QUANTEC (132) | IMRT parotid sparing preserves salivary flow | reproduced |
| **[HELD-OUT]** | KEYNOTE-048 (92), KEYNOTE-040 (93) | benefit increases with PD-L1 CPS | reproduced as a monotone gradient |
| **[HELD-OUT]** | Nordsmark 2005 (55), Begg 1992 (48) | hypoxic fraction 10-25 %, T_pot 4-5 days | reproduced as emergent values |
