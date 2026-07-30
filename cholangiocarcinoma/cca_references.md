# Cholangiocarcinoma QSP Model — References

Every PMID below was resolved against PubMed E-utilities while this model was
being built; the title shown is the title PubMed returned. References are
grouped by the part of the model they constrain, and the model parameter or
structural choice each one supports is named explicitly where there is one.

Trial-level anchors used to calibrate or to test the model are marked
**[ANCHOR]** (used to fit) or **[TEST]** (predicted, not fitted).

---

## 1. Pivotal systemic-therapy trials — the anchors

1. Valle J, Wasan H, Palmer DH, et al. **Cisplatin plus gemcitabine versus
   gemcitabine for biliary tract cancer.** *N Engl J Med* 2010;362:1273-81.
   PMID [20375404](https://pubmed.ncbi.nlm.nih.gov/20375404/) — **[ANCHOR]**
   ABC-02. Median OS 11.7 mo, median PFS 8.0 mo, ORR 26.1%. The reference
   point for `K_GEM`, `K_CIS` and the two hazard slopes.
2. Bridgewater J, Lopes A, Palmer D, et al. **Quality of life, long-term
   survivors and long-term outcome from the ABC-02 study.** *Br J Cancer*
   2016;114:965-71. PMID [27115567](https://pubmed.ncbi.nlm.nih.gov/27115567/)
   — the long-tail behaviour of the control arm.
3. Oh DY, Ruth He A, Qin S, et al. **Durvalumab plus gemcitabine and cisplatin
   in advanced biliary tract cancer.** *NEJM Evid* 2022;1(8).
   PMID [38319896](https://pubmed.ncbi.nlm.nih.gov/38319896/) — **[TEST]**
   TOPAZ-1. Median OS 12.8 vs 11.5 mo; 24-month OS 24.9% vs 10.4%. The
   observation that forced the mixture parameter `PI_IMMUNE`: a mean-effect
   model cannot move a tail without moving a median.
4. Oh DY, et al. **Durvalumab plus chemotherapy for advanced biliary tract
   cancer: a post hoc analysis of the TOPAZ-1 randomized clinical trial.**
   *JAMA Oncol* 2026. PMID [42424063](https://pubmed.ncbi.nlm.nih.gov/42424063/)
5. Kelley RK, Ueno M, Yoo C, et al. **Pembrolizumab in combination with
   gemcitabine and cisplatin compared with gemcitabine and cisplatin alone for
   patients with advanced biliary tract cancer (KEYNOTE-966).** *Lancet*
   2023;401:1853-65. PMID [37075781](https://pubmed.ncbi.nlm.nih.gov/37075781/)
   — the same small-median / separated-tail signature in a second trial.
6. Lamarca A, Palmer DH, Wasan HS, et al. **Second-line FOLFOX chemotherapy
   versus active symptom control for advanced biliary tract cancer (ABC-06).**
   *Lancet Oncol* 2021;22:690-701.
   PMID [33798493](https://pubmed.ncbi.nlm.nih.gov/33798493/) — the size of
   the second-line effect, which the model represents as `K_2L`.
7. Yoo C, Kim KP, Jeong JH, et al. **Liposomal irinotecan plus fluorouracil and
   leucovorin versus fluorouracil and leucovorin for metastatic biliary tract
   cancer after progression on gemcitabine plus cisplatin (NIFTY).**
   *Lancet Oncol* 2021;22:1560-72.
   PMID [34656226](https://pubmed.ncbi.nlm.nih.gov/34656226/)
8. Yoo C, et al. **Liposomal irinotecan for previously treated patients with
   biliary tract cancer: a pooled analysis of NIFTY and NALIRICC.**
   *J Hepatol* 2025. PMID [40147791](https://pubmed.ncbi.nlm.nih.gov/40147791/)
9. Primrose JN, Fox RP, Palmer DH, et al. **Capecitabine compared with
   observation in resected biliary tract cancer (BILCAP).** *Lancet Oncol*
   2019;20:663-73. PMID [30922733](https://pubmed.ncbi.nlm.nih.gov/30922733/)
   — **[TEST]** the adjuvant scenario. RFS 24.4 vs 17.5 mo.
10. Stein A, Arnold D, Bridgewater J, et al. **Adjuvant chemotherapy with
    gemcitabine and cisplatin compared to observation after curative intent
    resection of cholangiocarcinoma (ACTICCA-1).** *BMC Cancer* 2015;15:564.
    PMID [26228433](https://pubmed.ncbi.nlm.nih.gov/26228433/)
11. Eckel F, Schmid RM. **Chemotherapy in advanced biliary tract carcinoma: a
    pooled analysis of clinical trials.** *Br J Cancer* 2007;96:896-902.
    PMID [17325704](https://pubmed.ncbi.nlm.nih.gov/17325704/)
12. Glimelius B, Hoffman K, Sjödén PO, et al. **Chemotherapy improves survival
    and quality of life in advanced pancreatic and biliary cancer.**
    *Ann Oncol* 1996;7:593-600.
    PMID [8879373](https://pubmed.ncbi.nlm.nih.gov/8879373/) — the
    supportive-care comparator the S1 scenario is checked against.
13. Rimini M, et al. **Factors associated with reaching maintenance therapy in
    patients with advanced biliary tract cancer treated with durvalumab:
    real-world data.** *Int J Cancer* 2025.
    PMID [40387725](https://pubmed.ncbi.nlm.nih.gov/40387725/) — real-world
    attrition, which is what the dose gate is trying to reproduce.

## 2. Targeted therapy — FGFR2, IDH1, HER2, BRAF, NTRK, MSI

14. Abou-Alfa GK, Sahai V, Hollebecque A, et al. **Pemigatinib for previously
    treated, locally advanced or metastatic cholangiocarcinoma (FIGHT-202).**
    *Lancet Oncol* 2020;21:671-84.
    PMID [32203698](https://pubmed.ncbi.nlm.nih.gov/32203698/) — **[TEST]**
    ORR 35.5%, median PFS 6.9 mo, median OS 21.1 mo in the fusion cohort.
15. Goyal L, Meric-Bernstam F, Hollebecque A, et al. **Futibatinib for
    FGFR2-rearranged intrahepatic cholangiocarcinoma (FOENIX-CCA2).**
    *N Engl J Med* 2023;388:228-39.
    PMID [36652354](https://pubmed.ncbi.nlm.nih.gov/36652354/) — **[TEST]**
    the covalent inhibitor. Sets `COVAL = 1` and `RHO = RHO_FUT`.
16. Javle M, Roychowdhury S, Kelley RK, et al. **Infigratinib (BGJ398) in
    previously treated patients with advanced or metastatic
    cholangiocarcinoma with FGFR2 fusions or rearrangements.**
    *Lancet Gastroenterol Hepatol* 2021;6:803-15.
    PMID [34358484](https://pubmed.ncbi.nlm.nih.gov/34358484/)
17. Goyal L, Saha SK, Liu LY, et al. **Polyclonal secondary FGFR2 mutations
    drive acquired resistance to FGFR inhibition in patients with
    FGFR2 fusion-positive cholangiocarcinoma.** *Cancer Discov*
    2017;7:252-63. PMID [28034880](https://pubmed.ncbi.nlm.nih.gov/28034880/)
    — the observation behind structural commitment 2. Resistance mutations
    are POLYCLONAL, which is what a pre-existing-clone model predicts and an
    induced-resistance model does not.
18. Goyal L, Shi L, Liu LY, et al. **TAS-120 overcomes resistance to
    ATP-competitive FGFR inhibitors in patients with FGFR2 fusion-positive
    intrahepatic cholangiocarcinoma.** *Cancer Discov* 2019;9:1064-79.
    PMID [31109923](https://pubmed.ncbi.nlm.nih.gov/31109923/) — the
    cross-resistance asymmetry encoded as `RHO_PEM = 0.05` vs `RHO_FUT = 0.45`.
19. Arai Y, Totoki Y, Hosoda F, et al. **Fibroblast growth factor receptor 2
    tyrosine kinase fusions define a unique molecular subtype of
    cholangiocarcinoma.** *Hepatology* 2014;59:1427-34.
    PMID [24122810](https://pubmed.ncbi.nlm.nih.gov/24122810/)
20. Wu YM, Su F, Kalyana-Sundaram S, et al. **Identification of targetable FGFR
    gene fusions in diverse cancers.** *Cancer Discov* 2013;3:636-47.
    PMID [23558953](https://pubmed.ncbi.nlm.nih.gov/23558953/)
21. Abou-Alfa GK, Macarulla T, Javle MM, et al. **Ivosidenib in IDH1-mutant,
    chemotherapy-refractory cholangiocarcinoma (ClarIDHy).** *Lancet Oncol*
    2020;21:796-807. PMID [32416072](https://pubmed.ncbi.nlm.nih.gov/32416072/)
    — **[TEST]** ORR ~2% with a PFS hazard ratio of 0.37. This is the
    observation that forces a drug node with a growth-suppression term and a
    kill term of zero.
22. Zhu AX, Macarulla T, Javle MM, et al. **Final overall survival efficacy
    results of ivosidenib for patients with advanced cholangiocarcinoma with
    IDH1 mutation (ClarIDHy).** *JAMA Oncol* 2021;7:1669-77.
    PMID [34554208](https://pubmed.ncbi.nlm.nih.gov/34554208/)
23. Saha SK, Parachoniak CA, Ghanta KS, et al. **Mutant IDH inhibits HNF-4α to
    block hepatocyte differentiation and promote biliary cancer.** *Nature*
    2014;513:110-4. PMID [25043045](https://pubmed.ncbi.nlm.nih.gov/25043045/)
    — the differentiation block that makes the IDH1 axis cytostatic.
24. Wu MJ, Shi L, Dubrot J, et al. **Mutant IDH inhibits IFNγ-TET2 signaling to
    promote immunoevasion and tumor maintenance in cholangiocarcinoma.**
    *Cancer Discov* 2022;12:812-35.
    PMID [34848557](https://pubmed.ncbi.nlm.nih.gov/34848557/)
25. Harding JJ, Fan J, Oh DY, et al. **Zanidatamab for HER2-amplified,
    unresectable, locally advanced or metastatic biliary tract cancer
    (HERIZON-BTC-01).** *Lancet Oncol* 2023;24:772-82.
    PMID [37276871](https://pubmed.ncbi.nlm.nih.gov/37276871/)
26. Meric-Bernstam F, Hanna DL, El-Khoueiry AB, et al. **Zanidatamab, a novel
    bispecific antibody, for HER2-expressing or HER2-amplified cancers.**
    *Lancet Oncol* 2022;23:1558-70.
    PMID [36400106](https://pubmed.ncbi.nlm.nih.gov/36400106/)
27. Ohba A, Morizane C, Kawamoto Y, et al. **Trastuzumab deruxtecan in
    HER2-expressing biliary tract cancer (HERB; NCCH1805).** *J Clin Oncol*
    2024;42:3207-17. PMID [39102634](https://pubmed.ncbi.nlm.nih.gov/39102634/)
28. Subbiah V, Lassen U, Élez E, et al. **Dabrafenib plus trametinib in
    patients with BRAF^V600E-mutated biliary tract cancer (ROAR).**
    *Lancet Oncol* 2020;21:1234-43.
    PMID [32818466](https://pubmed.ncbi.nlm.nih.gov/32818466/)
29. Drilon A, Laetsch TW, Kummar S, et al. **Efficacy of larotrectinib in TRK
    fusion-positive cancers in adults and children.** *N Engl J Med*
    2018;378:731-9. PMID [29466156](https://pubmed.ncbi.nlm.nih.gov/29466156/)
30. Marabelle A, et al. **Pembrolizumab in microsatellite-instability-high and
    mismatch-repair-deficient advanced solid tumors: updated KEYNOTE-158
    results.** *Nat Cancer* 2025.
    PMID [39979665](https://pubmed.ncbi.nlm.nih.gov/39979665/)
31. Piha-Paul SA, Oh DY, Ueno M, et al. **Efficacy and safety of pembrolizumab
    for advanced biliary cancer: KEYNOTE-158 and KEYNOTE-028.**
    *Int J Cancer* 2020;147:2190-8.
    PMID [32359091](https://pubmed.ncbi.nlm.nih.gov/32359091/) — single-agent
    checkpoint blockade ORR ~6-13%, i.e. the size of `PI_IMMUNE`.
32. Ueno M, Ikeda M, Morizane C, et al. **Nivolumab alone or in combination
    with cisplatin plus gemcitabine in Japanese patients with unresectable or
    recurrent biliary tract cancer.** *Lancet Gastroenterol Hepatol*
    2019;4:611-21. PMID [31109808](https://pubmed.ncbi.nlm.nih.gov/31109808/)
33. Kommalapati A, Tella SH, Borad M, et al. **FGFR inhibitors in oncology:
    insight on the management of toxicities in clinical practice.**
    *Cancers (Basel)* 2021;13:2968.
    PMID [34199304](https://pubmed.ncbi.nlm.nih.gov/34199304/) — the
    hyperphosphataemia class effect, modelled as `EMAX_PHOS`.
34. Gattineni J, Baum M. **Genetic disorders of phosphate regulation.**
    *Pediatr Nephrol* 2012;27:1477-87.
    PMID [22350303](https://pubmed.ncbi.nlm.nih.gov/22350303/) — FGF23-Klotho
    control of NaPi-2a/2c, the on-target route to the phosphate rise.

## 3. Genomics, subtype and cell of origin

35. Nakamura H, Arai Y, Totoki Y, et al. **Genomic spectra of biliary tract
    cancer.** *Nat Genet* 2015;47:1003-10.
    PMID [26258846](https://pubmed.ncbi.nlm.nih.gov/26258846/)
36. Jusakul A, Cutcutache I, Yong CH, et al. **Whole-genome and epigenomic
    landscapes of etiologically distinct subtypes of cholangiocarcinoma.**
    *Cancer Discov* 2017;7:1116-35.
    PMID [28667006](https://pubmed.ncbi.nlm.nih.gov/28667006/)
37. Farshidfar F, Zheng S, Gingras MC, et al. **Integrative genomic analysis of
    cholangiocarcinoma identifies distinct IDH-mutant molecular profiles.**
    *Cell Rep* 2017;18:2780-94.
    PMID [28658632](https://pubmed.ncbi.nlm.nih.gov/28658632/)
38. Lowery MA, Ptashkin R, Jordan E, et al. **Comprehensive molecular profiling
    of intrahepatic and extrahepatic cholangiocarcinomas.** *Clin Cancer Res*
    2018;24:4154-61. PMID [29848569](https://pubmed.ncbi.nlm.nih.gov/29848569/)
    — the anatomic confinement of FGFR2 fusions and IDH1 mutations to
    intrahepatic disease, which the model encodes as the `FHILAR` covariate.
39. Boerner T, Drill E, Pak LM, et al. **Genetic determinants of outcome in
    intrahepatic cholangiocarcinoma.** *Hepatology* 2021;74:1429-44.
    PMID [33765338](https://pubmed.ncbi.nlm.nih.gov/33765338/)
40. Banales JM, Marin JJG, Lamarca A, et al. **Cholangiocarcinoma 2020: the
    next horizon in mechanisms and management.** *Nat Rev Gastroenterol
    Hepatol* 2020;17:557-88.
    PMID [32606456](https://pubmed.ncbi.nlm.nih.gov/32606456/)
41. Rimassa L, et al. **Mapping the landscape of biliary tract cancer in
    Europe: challenges and controversies.** *Lancet Reg Health Eur* 2025.
    PMID [40093398](https://pubmed.ncbi.nlm.nih.gov/40093398/)

## 4. Chronic biliary injury, aetiology, bile-acid signalling

42. Boonstra K, Weersma RK, van Erpecum KJ, et al. **Population-based
    epidemiology, malignancy risk, and outcome of primary sclerosing
    cholangitis.** *Hepatology* 2013;58:2045-55.
    PMID [23775876](https://pubmed.ncbi.nlm.nih.gov/23775876/)
43. Isomoto H, Mott JL, Kobayashi S, et al. **Sustained IL-6/STAT-3 signaling
    in cholangiocarcinoma cells due to SOCS-3 epigenetic silencing.**
    *Gastroenterology* 2007;132:384-96.
    PMID [17241887](https://pubmed.ncbi.nlm.nih.gov/17241887/) — the IL-6 node
    that drives both the stromal and the cachexia arms of the model.
44. Jaiswal M, LaRusso NF, Burgart LJ, Gores GJ. **Inflammatory cytokines
    induce DNA damage and inhibit DNA repair in cholangiocarcinoma cells by a
    nitric oxide-dependent mechanism.** *Cancer Res* 2000;60:184-90.
    PMID [10646872](https://pubmed.ncbi.nlm.nih.gov/10646872/) — the
    inflammation-to-mutation-rate edge that justifies a non-trivial `MU_RES`.
45. Wu T. **Cyclooxygenase-2 and prostaglandin signaling in
    cholangiocarcinoma.** *Biochim Biophys Acta* 2005;1755:135-50.
    PMID [15921858](https://pubmed.ncbi.nlm.nih.gov/15921858/)
46. Liu R, Li X, Hylemon PB, Zhou H. **Conjugated bile acids promote invasive
    growth via sphingosine 1-phosphate receptor 2.** *Am J Pathol*
    2018;188:2042-58. PMID [29963993](https://pubmed.ncbi.nlm.nih.gov/29963993/)
47. Yu H, Zhao T, Liu S, et al. **MRGPRX4 is a bile acid receptor for human
    cholestatic itch.** *eLife* 2019;8:e48431.
    PMID [31500698](https://pubmed.ncbi.nlm.nih.gov/31500698/) — the pruritus
    endpoint on the map.

## 5. Desmoplastic stroma and the drug-penetration barrier

48. Sirica AE. **The role of cancer-associated myofibroblasts in intrahepatic
    cholangiocarcinoma.** *Nat Rev Gastroenterol Hepatol* 2011;9:44-54.
    PMID [22143274](https://pubmed.ncbi.nlm.nih.gov/22143274/)
49. Fabris L, Perugorria MJ, Mertens J, et al. **The tumour microenvironment
    and immune milieu of cholangiocarcinoma.** *Liver Int* 2019;39(Suppl 1):
    63-78. PMID [30907492](https://pubmed.ncbi.nlm.nih.gov/30907492/)
50. Cadamuro M, Brivio S, Mertens J, et al. **Platelet-derived growth factor-D
    enables liver myofibroblasts to promote tumor lymphangiogenesis in
    cholangiocarcinoma.** *J Hepatol* 2019;70:700-9.
    PMID [30553841](https://pubmed.ncbi.nlm.nih.gov/30553841/)
51. Manzanares MÁ, Campbell DJW, Maldonado GT, Sirica AE. **Overexpression of
    periostin and distinct mesothelin forms predict malignant progression in a
    rat cholangiocarcinoma model.** *Hepatol Commun* 2018;2:155-72.
    PMID [29404524](https://pubmed.ncbi.nlm.nih.gov/29404524/)
52. Heldin CH, Rubin K, Pietras K, Östman A. **High interstitial fluid
    pressure — an obstacle in cancer therapy.** *Nat Rev Cancer* 2004;4:806-13.
    PMID [15510161](https://pubmed.ncbi.nlm.nih.gov/15510161/) — the physical
    argument for `F_PEN`, the penetration factor that multiplies every
    cytotoxic exposure at the tumour.
53. Provenzano PP, Cuevas C, Chang AE, et al. **Enzymatic targeting of the
    stroma ablates physical barriers to treatment of pancreatic ductal
    adenocarcinoma.** *Cancer Cell* 2012;21:418-29.
    PMID [22439937](https://pubmed.ncbi.nlm.nih.gov/22439937/)
54. Olive KP, Jacobetz MA, Davidson CJ, et al. **Inhibition of Hedgehog
    signaling enhances delivery of chemotherapy in a mouse model of pancreatic
    cancer.** *Science* 2009;324:1457-61.
    PMID [19460966](https://pubmed.ncbi.nlm.nih.gov/19460966/) — the
    experiment the F4 falsifier is a simulation of.

## 6. Biliary drainage, stents and cholangitis — the gate

55. Sangchan A, Kongkasame W, Pugkhem A, et al. **Efficacy of metal and plastic
    stents in unresectable complex hilar cholangiocarcinoma: a randomized
    controlled trial.** *Gastrointest Endosc* 2012;76:93-9.
    PMID [22595446](https://pubmed.ncbi.nlm.nih.gov/22595446/) — **[TEST]**
    the S5 scenario. Sets `KOCC_SEMS` and `KOCC_PLAS`.
56. Chen X, et al. **Efficacy and safety of preoperative biliary drainage in
    patients with hilar cholangiocarcinoma: a systematic review and
    meta-analysis.** *Int J Surg* 2025.
    PMID [40072352](https://pubmed.ncbi.nlm.nih.gov/40072352/)
57. Rousian M, et al. **Primary percutaneous stenting above the ampulla versus
    endoscopic drainage for unresectable malignant hilar biliary obstruction.**
    *BMC Cancer* 2025. PMID [40346549](https://pubmed.ncbi.nlm.nih.gov/40346549/)
58. Yokoe M, Hata J, Takada T, et al. **Tokyo Guidelines 2018: diagnostic
    criteria and severity grading of acute cholecystitis.** *J Hepatobiliary
    Pancreat Sci* 2018;25:41-54.
    PMID [29032636](https://pubmed.ncbi.nlm.nih.gov/29032636/) — the companion
    Tokyo Guidelines severity framework for biliary infection; the model's
    `CHOLI` state is a continuous analogue of that grading.
59. Ortner ME, Caca K, Berr F, et al. **Successful photodynamic therapy for
    nonresectable cholangiocarcinoma: a randomized prospective study.**
    *Gastroenterology* 2003;125:1355-63.
    PMID [14598251](https://pubmed.ncbi.nlm.nih.gov/14598251/) — a
    patency-directed intervention with a survival benefit and no systemic
    antitumour action, which is the clinical form of the model's central claim.
60. Johnson PJ, Berhane S, Kagebayashi C, et al. **Assessment of liver function
    in patients with hepatocellular carcinoma: a new evidence-based approach —
    the ALBI grade.** *J Clin Oncol* 2015;33:550-8.
    PMID [25512453](https://pubmed.ncbi.nlm.nih.gov/25512453/) — the exact
    formula used for `ALBI` in `$TABLE`, and the driver of `h_biliary`.

## 7. Surgery and locoregional therapy

61. Nagino M, Ebata T, Yokoyama Y, et al. **Evolution of surgical treatment for
    perihilar cholangiocarcinoma: a single-center 34-year review of 574
    consecutive resections.** *Ann Surg* 2013;258:129-40.
    PMID [23059502](https://pubmed.ncbi.nlm.nih.gov/23059502/)
62. Ribero D, Pinna AD, Guglielmi A, et al. **Surgical approach for long-term
    survival of patients with intrahepatic cholangiocarcinoma: a
    multi-institutional analysis of 434 patients.** *Arch Surg* 2012;147:1107-13.
    PMID [22910846](https://pubmed.ncbi.nlm.nih.gov/22910846/)
63. Groot Koerkamp B, Wiggers JK, Allen PJ, et al. **Survival after resection of
    perihilar cholangiocarcinoma — development and external validation of a
    prognostic nomogram.** *Ann Oncol* 2016;27:753.
    PMID [26920702](https://pubmed.ncbi.nlm.nih.gov/26920702/)
64. Weber SM, Ribero D, O'Reilly EM, et al. **Intrahepatic cholangiocarcinoma:
    expert consensus statement.** *HPB (Oxford)* 2015;17:669-80.
    PMID [26172134](https://pubmed.ncbi.nlm.nih.gov/26172134/)
65. Darwish Murad S, Kim WR, Harnois DM, et al. **Efficacy of neoadjuvant
    chemoradiation, followed by liver transplantation, for perihilar
    cholangiocarcinoma at 12 US centers.** *Gastroenterology* 2012;143:88-98.
    PMID [22504095](https://pubmed.ncbi.nlm.nih.gov/22504095/)
66. Franssen S, Soares KC, Jolissaint JS, et al. **Comparison of hepatic
    arterial infusion pump chemotherapy vs resection for patients with
    multifocal intrahepatic cholangiocarcinoma.** *JAMA Surg* 2022;157:590-6.
    PMID [35544131](https://pubmed.ncbi.nlm.nih.gov/35544131/)
67. Edeline J, Touchefeu Y, Guiu B, et al. **Radioembolization plus chemotherapy
    for first-line treatment of locally advanced intrahepatic
    cholangiocarcinoma (MISPHEC).** *JAMA Oncol* 2020;6:51-9.
    PMID [31670746](https://pubmed.ncbi.nlm.nih.gov/31670746/)
68. Grassberger C, Hong TS, Hato T, et al. **Differential association between
    circulating lymphocyte populations and outcome after radiation therapy in
    subtypes of liver cancer.** *Int J Radiat Oncol Biol Phys* 2018;101:1222-5.
    PMID [29859792](https://pubmed.ncbi.nlm.nih.gov/29859792/)

## 8. Drug pharmacology — the PK/PD submodels

69. Abbruzzese JL, Grunewald R, Weeks EA, et al. **A phase I clinical, plasma,
    and cellular pharmacology study of gemcitabine.** *J Clin Oncol*
    1991;9:491-8. PMID [1999720](https://pubmed.ncbi.nlm.nih.gov/1999720/) —
    `CL_GEM`, `V1_GEM`, and the saturable dFdCTP accumulation that is
    represented by `VMAX_TP`/`KM_TP`.
70. Plunkett W, Huang P, Xu YZ, et al. **Gemcitabine: metabolism, mechanisms of
    action, and self-potentiation.** *Semin Oncol* 1995;22(4 Suppl 11):3-10.
    PMID [7481842](https://pubmed.ncbi.nlm.nih.gov/7481842/) — the RRM1
    self-potentiation loop drawn on the map.
71. Bergman AM, Pinedo HM, Peters GJ. **Determinants of resistance to
    2',2'-difluorodeoxycytidine (gemcitabine).** *Drug Resist Updat*
    2002;5:19-33. PMID [12127861](https://pubmed.ncbi.nlm.nih.gov/12127861/)
72. Vos LJ, Yusuf D, Lui A, et al. **Predictive and prognostic properties of
    human equilibrative nucleoside transporter 1 expression in
    gemcitabine-treated pancreatobiliary cancer.** *JCO Precis Oncol* 2019.
    PMID [35100740](https://pubmed.ncbi.nlm.nih.gov/35100740/) — the hENT1
    node.
73. Kelland L. **The resurgence of platinum-based cancer chemotherapy.**
    *Nat Rev Cancer* 2007;7:573-84.
    PMID [17625587](https://pubmed.ncbi.nlm.nih.gov/17625587/)
74. Siddik ZH. **Cisplatin: mode of cytotoxic action and molecular basis of
    resistance.** *Oncogene* 2003;22:7265-79.
    PMID [14576837](https://pubmed.ncbi.nlm.nih.gov/14576837/) — adduct
    formation and ERCC1-dependent removal, the `KFORM_PT`/`KREP_PT` pair.
75. Friberg LE, Henningsson A, Maas H, et al. **Model of chemotherapy-induced
    myelosuppression with parameter consistency across drugs.** *J Clin Oncol*
    2002;20:4713-21. PMID [12488418](https://pubmed.ncbi.nlm.nih.gov/12488418/)
    — the exact structure of the ANC submodel: proliferating pool, three
    transit compartments, `MTT` and the feedback exponent `GAMMA`.
76. Terranova N, Germani M, Del Bene F, Magni P. **A predictive
    pharmacokinetic-pharmacodynamic model of tumor growth kinetics in xenograft
    mice after administration of anticancer agents.** *Cancer Chemother
    Pharmacol* 2013;72:471-82.
    PMID [23812004](https://pubmed.ncbi.nlm.nih.gov/23812004/) — the
    Simeoni-type exponential-to-linear growth switch and the damaged-cell
    transit chain used for `TD1`-`TD3`.

## 9. Resistance as clonal selection, and drug-tolerant persistence

77. Goldie JH, Coldman AJ. **A mathematic model for relating the drug
    sensitivity of tumors to their spontaneous mutation rate.** *Cancer Treat
    Rep* 1979;63:1727-33.
    PMID [526911](https://pubmed.ncbi.nlm.nih.gov/526911/) — the seeding rule
    `T_R(0) = T(0) · μ · ln(N)` used verbatim. There is no fitted
    "time to resistance" anywhere in the model because of this reference.
78. Sharma SV, Lee DY, Li B, et al. **A chromatin-mediated reversible
    drug-tolerant state in cancer cell subpopulations.** *Cell*
    2010;141:69-80. PMID [20371346](https://pubmed.ncbi.nlm.nih.gov/20371346/)
    — the persister compartment `T_P`: reversible, slow-cycling, NOT mutated.
79. Vasan N, Baselga J, Hyman DM. **A view on drug resistance in cancer.**
    *Nature* 2019;575:299-309.
    PMID [31723286](https://pubmed.ncbi.nlm.nih.gov/31723286/)
80. Bozic I, Reiter JG, Allen B, et al. **Evolutionary dynamics of cancer in
    response to targeted combination therapy.** *eLife* 2013;2:e00747.
    PMID [23805382](https://pubmed.ncbi.nlm.nih.gov/23805382/) — the
    arithmetic that shows a clone seeded at 10⁻⁵ cannot dominate inside a
    six-month window, which is why this model reports that negative result
    instead of tuning around it.
81. Diaz LA Jr, Williams RT, Wu J, et al. **The molecular evolution of acquired
    resistance to targeted EGFR blockade in colorectal cancers.** *Nature*
    2012;486:537-40. PMID [22722843](https://pubmed.ncbi.nlm.nih.gov/22722843/)
82. Bettegowda C, Sausen M, Leary RJ, et al. **Detection of circulating tumor
    DNA in early- and late-stage human malignancies.** *Sci Transl Med*
    2014;6:224ra24. PMID [24553385](https://pubmed.ncbi.nlm.nih.gov/24553385/)
    — the shedding-versus-burden relation behind the `CTDNA` compartment.
83. Cowzer D, et al. **Clinical utility and prognostic implications of
    circulating cell-free DNA in biliary tract cancer.** *JCO Precis Oncol*
    2025. PMID [40956997](https://pubmed.ncbi.nlm.nih.gov/40956997/)

## 10. Guidelines and reviews used for regimen definitions

84. Vogel A, Bridgewater J, Edeline J, et al. **Biliary tract cancer: ESMO
    Clinical Practice Guideline for diagnosis, treatment and follow-up.**
    *Ann Oncol* 2023;34:127-40.
    PMID [36372281](https://pubmed.ncbi.nlm.nih.gov/36372281/) — the source of
    the bilirubin and performance-status thresholds encoded in rules 1 and 4.
85. Bridgewater J, Galle PR, Khan SA, et al. **Guidelines for the diagnosis and
    management of intrahepatic cholangiocarcinoma.** *J Hepatol*
    2014;60:1268-89. PMID [24681130](https://pubmed.ncbi.nlm.nih.gov/24681130/)
86. Kam AE, Masood A, Shroff RT. **Current and emerging therapies for advanced
    biliary tract cancers.** *Lancet Gastroenterol Hepatol* 2021;6:956-69.
    PMID [34626563](https://pubmed.ncbi.nlm.nih.gov/34626563/)
87. Lamarca A, Barriuso J, McNamara MG, Valle JW. **Molecular targeted
    therapies: ready for "prime time" in biliary tract cancer.** *J Hepatol*
    2020;73:170-85. PMID [32171892](https://pubmed.ncbi.nlm.nih.gov/32171892/)

---

## How the references map onto the three structural commitments

| Commitment | Load-bearing references |
|---|---|
| **I — the delivered dose is an output** | 1, 6, 9, 13, 55-60, 84. Rules 1-4 are the ESMO/label dose-modification thresholds; the stent-patency constants come from the randomised metal-versus-plastic data; ALBI is Johnson's published formula. Reference 59 is the clinical form of the claim: an intervention with no antitumour action that changes survival. |
| **II — resistance is selected, not induced** | 17, 18, 77-83. Polyclonality at progression (17) is the observation; Goldie-Coldman (77) supplies the seeding rule; Bozic (80) supplies the arithmetic that limits what that rule can explain; Sharma (78) supplies the mechanism the model uses instead for median PFS. |
| **III — two competing hazards on different clocks** | 3, 5, 31, 32, 60. The tail-without-median signature in TOPAZ-1 (3) and KEYNOTE-966 (5), plus the ~10% single-agent checkpoint response rate (31, 32) that sizes the mixture fraction, plus ALBI (60) for the biliary clock. |
