# Diabetic Foot Ulcer (DFU) — Reference List

Literature underpinning the mechanistic map (`dfu_qsp_model.dot`), the ODE model
(`dfu_mrgsolve_model.R` / `dfu_reference_implementation.py`) and the Shiny
dashboard (`dfu_shiny_app.R`).

Entries marked **[ANCHOR]** are the sources the model is quantitatively
calibrated against; the corresponding model values are reported in the header of
`dfu_mrgsolve_model.R` and reproducible with
`python3 dfu_reference_implementation.py --anchors`.

---

## 1. Epidemiology, burden and natural history

1. Armstrong DG, Boulton AJM, Bus SA. **Diabetic foot ulcers and their recurrence.**
   *N Engl J Med* 2017;376:2367-2375. **[ANCHOR — ~40% recurrence at 1 y, ~60% at 3 y, ~65% at 5 y]**
   <https://pubmed.ncbi.nlm.nih.gov/28614678/>
2. Armstrong DG, Tan T-W, Boulton AJM, Bus SA. **Diabetic foot ulcers: a review.**
   *JAMA* 2023;330:62-75.
   <https://pubmed.ncbi.nlm.nih.gov/37395769/>
3. Zhang P, Lu J, Jing Y, et al. **Global epidemiology of diabetic foot ulceration: a systematic review and meta-analysis.**
   *Ann Med* 2017;49:106-116.
   <https://pubmed.ncbi.nlm.nih.gov/27585063/>
4. Boulton AJM, Vileikyte L, Ragnarson-Tennvall G, Apelqvist J. **The global burden of diabetic foot disease.**
   *Lancet* 2005;366:1719-1724.
   <https://pubmed.ncbi.nlm.nih.gov/16291066/>
5. Prompers L, Huijberts M, Apelqvist J, et al. **High prevalence of ischaemia, infection and serious comorbidity in patients with diabetic foot disease in Europe (Eurodiale).**
   *Diabetologia* 2007;50:18-25.
   <https://pubmed.ncbi.nlm.nih.gov/17093942/>
6. Walsh JW, Hoffstad OJ, Sullivan MO, Margolis DJ. **Association of diabetic foot ulcer and death in a population-based cohort.**
   *Diabet Med* 2016;33:1493-1498.
   <https://pubmed.ncbi.nlm.nih.gov/26666583/>
7. Reardon R, Simring D, Kim B, et al. **The diabetic foot ulcer.**
   *Aust J Gen Pract* 2020;49:250-255.
   <https://pubmed.ncbi.nlm.nih.gov/32416652/>

## 2. Diabetic peripheral neuropathy and loss of protective sensation

8. Feldman EL, Callaghan BC, Pop-Busui R, et al. **Diabetic neuropathy.**
   *Nat Rev Dis Primers* 2019;5:41.
   <https://pubmed.ncbi.nlm.nih.gov/31197183/>
9. Boulton AJM, Vinik AI, Arezzo JC, et al. **Diabetic neuropathies: a statement by the American Diabetes Association.**
   *Diabetes Care* 2005;28:956-962.
   <https://pubmed.ncbi.nlm.nih.gov/15793206/>
10. Brownlee M. **The pathobiology of diabetic complications: a unifying mechanism.**
    *Diabetes* 2005;54:1615-1625. (polyol, hexosamine, PKC, AGE — the four pathways in cluster 1)
    <https://pubmed.ncbi.nlm.nih.gov/15919781/>
11. Pop-Busui R, Boulton AJM, Feldman EL, et al. **Diabetic neuropathy: a position statement by the ADA.**
    *Diabetes Care* 2017;40:136-154.
    <https://pubmed.ncbi.nlm.nih.gov/27999003/>
12. Tesfaye S, Chaturvedi N, Eaton SEM, et al. **Vascular risk factors and diabetic neuropathy (EURODIAB).**
    *N Engl J Med* 2005;352:341-350. (basis for the slow glycaemic drive of the NEURO state)
    <https://pubmed.ncbi.nlm.nih.gov/15673800/>
13. Callaghan BC, Little AA, Feldman EL, Hughes RAC. **Enhanced glucose control for preventing and treating diabetic neuropathy.**
    *Cochrane Database Syst Rev* 2012;6:CD007543. (why the neuropathy axis is slow)
    <https://pubmed.ncbi.nlm.nih.gov/22696371/>

## 3. Biomechanics, plantar pressure and offloading

14. Armstrong DG, Nguyen HC, Lavery LA, et al. **Off-loading the diabetic foot wound: a randomized clinical trial.**
    *Diabetes Care* 2001;24:1019-1022. **[ANCHOR — mean days to healing: TCC 33.5, removable cast walker 50.4, half-shoe 61.0]**
    <https://pubmed.ncbi.nlm.nih.gov/11375363/>
15. Armstrong DG, Lavery LA, Kimbriel HR, et al. **Activity patterns of patients with diabetic foot ulceration: patients with active ulceration may not adhere to a standard pressure off-loading regimen.**
    *Diabetes Care* 2003;26:2595-2597. **[ANCHOR — removable device worn during only 28% of daily activity]**
    <https://pubmed.ncbi.nlm.nih.gov/12941724/>
16. Armstrong DG, Lavery LA, Wu S, Boulton AJM. **Evaluation of removable and irremovable cast walkers in the healing of diabetic foot wounds.**
    *Diabetes Care* 2005;28:551-554. (instant TCC ≈ TCC — the model's `ADHERENCE = 1` case)
    <https://pubmed.ncbi.nlm.nih.gov/15735186/>
17. Bus SA, Armstrong DG, Crews RT, et al. **Guidelines on offloading foot ulcers in persons with diabetes (IWGDF 2023).**
    *Diabetes Metab Res Rev* 2024;40:e3647.
    <https://pubmed.ncbi.nlm.nih.gov/37226568/>
18. Murray HJ, Young MJ, Hollis S, Boulton AJM. **The association between callus formation, high pressures and neuropathy in diabetic foot ulceration.**
    *Diabet Med* 1996;13:979-982. (the 26-fold callus risk in the map)
    <https://pubmed.ncbi.nlm.nih.gov/8946155/>
19. Veves A, Murray HJ, Young MJ, Boulton AJM. **The risk of foot ulceration in diabetic patients with high foot pressure: a prospective study.**
    *Diabetologia* 1992;35:660-663.
    <https://pubmed.ncbi.nlm.nih.gov/1644245/>
20. Mueller MJ, Sinacore DR, Hastings MK, et al. **Effect of Achilles tendon lengthening on neuropathic plantar ulcers.**
    *J Bone Joint Surg Am* 2003;85:1436-1445.
    <https://pubmed.ncbi.nlm.nih.gov/12925622/>
21. Rogers LC, Frykberg RG, Armstrong DG, et al. **The Charcot foot in diabetes.**
    *Diabetes Care* 2011;34:2123-2129.
    <https://pubmed.ncbi.nlm.nih.gov/21868781/>

## 4. Perfusion, PAD and the oxygen gate

22. Conte MS, Bradbury AW, Kolh P, et al. **Global vascular guidelines on the management of chronic limb-threatening ischemia.**
    *J Vasc Surg* 2019;69:3S-125S. (WIfI staging, revascularisation thresholds)
    <https://pubmed.ncbi.nlm.nih.gov/31159978/>
23. Mills JL Sr, Conte MS, Armstrong DG, et al. **The SVS lower extremity threatened limb classification system: WIfI.**
    *J Vasc Surg* 2014;59:220-234.
    <https://pubmed.ncbi.nlm.nih.gov/24126108/>
24. Fife CE, Buyukcakir C, Otto GH, et al. **The predictive value of transcutaneous oxygen tension measurement in diabetic lower extremity ulcers treated with hyperbaric oxygen therapy.**
    *Wound Repair Regen* 2002;10:198-207. (the TcPO₂ threshold behaviour reproduced by the model's Hill gate)
    <https://pubmed.ncbi.nlm.nih.gov/12191001/>
25. Hunt TK, Pai MP. **The effect of varying ambient oxygen tensions on wound metabolism and collagen synthesis.**
    *Surg Gynecol Obstet* 1972;135:561-567. (prolyl hydroxylase is O₂-dependent — the anabolic gate)
    <https://pubmed.ncbi.nlm.nih.gov/5077722/>
26. Allen DB, Maguire JJ, Mahdavian M, et al. **Wound hypoxia and acidosis limit neutrophil bacterial killing mechanisms.**
    *Arch Surg* 1997;132:991-996.
    <https://pubmed.ncbi.nlm.nih.gov/9301612/>
27. Catrina S-B, Zheng X. **Hypoxia and hypoxia-inducible factors in diabetes and its complications.**
    *Diabetologia* 2021;64:709-716.
    <https://pubmed.ncbi.nlm.nih.gov/33496820/>
28. Thangarajah H, Yao D, Chang EI, et al. **The molecular basis for impaired hypoxia-induced VEGF expression in diabetic tissues.**
    *Proc Natl Acad Sci USA* 2009;106:13505-13510. (methylglyoxal modification of HIF-1α/p300 — the `MGO_BLOCK` term)
    <https://pubmed.ncbi.nlm.nih.gov/19666581/>
29. Botusan IR, Sunkari VG, Savu O, et al. **Stabilization of HIF-1α is critical to improve wound healing in diabetic mice.**
    *Proc Natl Acad Sci USA* 2008;105:19426-19431.
    <https://pubmed.ncbi.nlm.nih.gov/19057015/>
30. Fejfarová V, Jirkovská A, Dubský M, et al. **Microcirculation and diabetic foot.** (functional microvascular ischaemia with a normal ABI)
    *Physiol Res* 2021;70:S251-S262.
    <https://pubmed.ncbi.nlm.nih.gov/34913354/>

## 5. Chronic wound inflammation and the failed M1→M2 switch

31. Eming SA, Martin P, Tomic-Canic M. **Wound repair and regeneration: mechanisms, signaling, and translation.**
    *Sci Transl Med* 2014;6:265sr6.
    <https://pubmed.ncbi.nlm.nih.gov/25473038/>
32. Mirza RE, Fang MM, Ennis WJ, Koh TJ. **Blocking interleukin-1β induces a healing-associated wound macrophage phenotype and improves healing in type 2 diabetes.**
    *Diabetes* 2013;62:2579-2587. (the IL-1β self-amplifying loop and its role in the switch)
    <https://pubmed.ncbi.nlm.nih.gov/23493576/>
33. Boniakowski AE, Kimball AS, Jacobs BN, et al. **Macrophage-mediated inflammation in normal and diabetic wound healing.**
    *J Immunol* 2017;199:17-24.
    <https://pubmed.ncbi.nlm.nih.gov/28630109/>
34. Wong SL, Demers M, Martinod K, et al. **Diabetes primes neutrophils to undergo NETosis, which impairs wound healing.**
    *Nat Med* 2015;21:815-819.
    <https://pubmed.ncbi.nlm.nih.gov/26076037/>
35. Khanna S, Biswas S, Shang Y, et al. **Macrophage dysfunction impairs resolution of inflammation in the wounds of diabetic mice.**
    *PLoS One* 2010;5:e9539.
    <https://pubmed.ncbi.nlm.nih.gov/20209061/>
36. Bannon P, Wood S, Restivo T, et al. **Diabetes induces stable intrinsic changes to myeloid cells that drive chronic inflammation.**
    *Dis Model Mech* 2013;6:1434-1447.
    <https://pubmed.ncbi.nlm.nih.gov/24057002/>
37. Delamaire M, Maugendre D, Moreno M, et al. **Impaired leucocyte functions in diabetic patients.**
    *Diabet Med* 1997;14:29-34.
    <https://pubmed.ncbi.nlm.nih.gov/9017350/>

## 6. Proteases, growth factors and the stalled wound bed

38. Lobmann R, Ambrosch A, Schultz G, et al. **Expression of matrix-metalloproteinases and their inhibitors in the wounds of diabetic and non-diabetic patients.**
    *Diabetologia* 2002;45:1011-1016. (MMP-9 and MMP-8 markedly elevated; TIMP-2 reduced)
    <https://pubmed.ncbi.nlm.nih.gov/12136400/>
39. Muller M, Trocme C, Lardy B, et al. **Matrix metalloproteinases and diabetic foot ulcers: the ratio of MMP-1 to TIMP-1 is a predictor of wound healing.**
    *Diabet Med* 2008;25:419-426. **[ANCHOR — the protease-ratio predictor structure]**
    <https://pubmed.ncbi.nlm.nih.gov/18387077/>
40. Liu Y, Min D, Bolton T, et al. **Increased matrix metalloproteinase-9 predicts poor wound healing in diabetic foot ulcers.**
    *Diabetes Care* 2009;32:117-119.
    <https://pubmed.ncbi.nlm.nih.gov/18835949/>
41. Trengove NJ, Stacey MC, MacAuley S, et al. **Analysis of the acute and chronic wound environments: the role of proteases and their inhibitors.**
    *Wound Repair Regen* 1999;7:442-452. (chronic wound fluid degrades applied growth factor — the `KDEG_BEC_PROT` term)
    <https://pubmed.ncbi.nlm.nih.gov/10633003/>
42. Schultz GS, Wysocki A. **Interactions between extracellular matrix and growth factors in wound healing.**
    *Wound Repair Regen* 2009;17:153-162.
    <https://pubmed.ncbi.nlm.nih.gov/19320882/>
43. Wall IB, Moseley R, Baird DM, et al. **Fibroblast dysfunction is a key factor in the non-healing of chronic venous leg ulcers.**
    *J Invest Dermatol* 2008;128:2526-2540. (senescent fibroblast phenotype)
    <https://pubmed.ncbi.nlm.nih.gov/18449211/>
44. Stojadinovic O, Brem H, Vouthounis C, et al. **Molecular pathogenesis of chronic wounds: the role of β-catenin and c-myc in the inhibition of epithelialization and wound healing.**
    *Am J Pathol* 2005;167:59-69. (keratinocytes proliferate but do not migrate — the DFU paradox)
    <https://pubmed.ncbi.nlm.nih.gov/15972952/>
45. Fadini GP, Albiero M, Bonora BM, Avogaro A. **Angiogenic abnormalities in diabetes mellitus: mechanism-based therapeutic approaches.**
    *J Clin Endocrinol Metab* 2019;104:5431-5444. (diabetic stem-cell mobilopathy — the `mobfail` term)
    <https://pubmed.ncbi.nlm.nih.gov/31211371/>
46. Fadini GP, Sartore S, Albiero M, et al. **Number and function of endothelial progenitor cells as a marker of severity for diabetic vasculopathy.**
    *Arterioscler Thromb Vasc Biol* 2006;26:2140-2146.
    <https://pubmed.ncbi.nlm.nih.gov/16857948/>
47. Fadini GP, Albiero M, Seeger F, et al. **Stem cell compartmentalization in diabetes and high cardiovascular risk reveals the role of DPP-4 in diabetic stem cell mobilopathy.**
    *Basic Res Cardiol* 2013;108:313. (DPP-4 cleavage of SDF-1α — the `K_DPP4` / `DPP4I` terms)
    <https://pubmed.ncbi.nlm.nih.gov/23184393/>

## 7. Wound-area kinetics, the perimeter law and PAR₄

48. Sheehan P, Jones P, Caselli A, et al. **Percent change in wound area of diabetic foot ulcers over a 4-week period is a robust predictor of complete healing in a 12-week prospective trial.**
    *Diabetes Care* 2003;26:1879-1882. **[ANCHOR — the PAR₄ > 50% decision rule]**
    <https://pubmed.ncbi.nlm.nih.gov/12766127/>
49. Gilman TH. **Parameter for measurement of wound closure.**
    *Wounds* 1990;2:95-101. (the linear edge-advance / perimeter formulation the model uses: dA/dt = −k·P)
50. Cardinal M, Eisenbud DE, Phillips T, Harding K. **Early healing rates and wound area measurements are reliable predictors of later complete wound closure.**
    *Wound Repair Regen* 2008;16:19-22.
    <https://pubmed.ncbi.nlm.nih.gov/18086284/>
51. Coerper S, Beckert S, Küper MA, et al. **Fifty percent area reduction after 4 weeks of treatment is a reliable indicator for healing.**
    *J Diabetes Complications* 2009;23:49-53.
    <https://pubmed.ncbi.nlm.nih.gov/18024178/>
52. Margolis DJ, Allen-Taylor L, Hoffstad O, Berlin JA. **Diabetic neuropathic foot ulcers: the association of wound size, wound duration, and wound grade on healing.**
    *Diabetes Care* 2002;25:1835-1839. (size, duration and grade as the three dominant covariates)
    <https://pubmed.ncbi.nlm.nih.gov/12351487/>
53. Margolis DJ, Kantor J, Berlin JA. **Healing of diabetic neuropathic foot ulcers receiving standard treatment: a meta-analysis.**
    *Diabetes Care* 1999;22:692-695. (the 24-31% 20-week closure rate under true standard care)
    <https://pubmed.ncbi.nlm.nih.gov/10332667/>

## 8. Debridement, biofilm and infection

54. Steed DL, Donohoe D, Webster MW, Lindsley L. **Effect of extensive debridement and treatment on the healing of diabetic foot ulcers.**
    *J Am Coll Surg* 1996;183:61-64.
    <https://pubmed.ncbi.nlm.nih.gov/8673309/>
55. Wilcox JR, Carter MJ, Covington S. **Frequency of debridements and time to heal: a retrospective cohort study of 312,744 wounds.**
    *JAMA Dermatol* 2013;149:1050-1058. **[ANCHOR — debridement frequency, not technique, tracks healing]**
    <https://pubmed.ncbi.nlm.nih.gov/23884238/>
56. Malone M, Bjarnsholt T, McBain AJ, et al. **The prevalence of biofilms in chronic wounds: a systematic review and meta-analysis.**
    *J Wound Care* 2017;26:20-25. (biofilm in ~78% of chronic wounds)
    <https://pubmed.ncbi.nlm.nih.gov/28103163/>
57. Wolcott RD, Rumbaugh KP, James G, et al. **Biofilm maturity studies indicate sharp debridement opens a time-dependent therapeutic window.**
    *J Wound Care* 2010;19:320-328. **[ANCHOR — biofilm returns to tolerance within 24-72 h; the `KG_BIOF` regrowth rate]**
    <https://pubmed.ncbi.nlm.nih.gov/20852503/>
58. Ceri H, Olson ME, Stremick C, et al. **The Calgary Biofilm Device: new technology for rapid determination of antibiotic susceptibilities of bacterial biofilms.**
    *J Clin Microbiol* 1999;37:1771-1776. (100-1000× tolerance — the `TOL_BIOF` parameter)
    <https://pubmed.ncbi.nlm.nih.gov/10325322/>
59. Lipsky BA, Berendt AR, Cornia PB, et al. **2012 IDSA clinical practice guideline for the diagnosis and treatment of diabetic foot infections.**
    *Clin Infect Dis* 2012;54:e132-e173.
    <https://pubmed.ncbi.nlm.nih.gov/22619242/>
60. Senneville É, Albalawi Z, van Asten SA, et al. **IWGDF/IDSA guidelines on the diagnosis and treatment of diabetes-related foot infections (2023).**
    *Diabetes Metab Res Rev* 2024;40:e3687.
    <https://pubmed.ncbi.nlm.nih.gov/37779323/>
61. Lam K, van Asten SAV, Nguyen T, et al. **Diagnostic accuracy of probe to bone to detect osteomyelitis in the diabetic foot: a systematic review.**
    *Clin Infect Dis* 2016;63:944-948.
    <https://pubmed.ncbi.nlm.nih.gov/27369321/>
62. Tone A, Nguyen S, Devemy F, et al. **Six-week versus twelve-week antibiotic therapy for nonsurgically treated diabetic foot osteomyelitis: a multicenter open-label controlled randomized study.**
    *Diabetes Care* 2015;38:302-307.
    <https://pubmed.ncbi.nlm.nih.gov/25414157/>
63. Lázaro-Martínez JL, Aragón-Sánchez J, García-Morales E. **Antibiotics versus conservative surgery for treating diabetic foot osteomyelitis: a randomized comparative trial.**
    *Diabetes Care* 2014;37:789-795. **[ANCHOR — resection vs antibiotic-only for the bone nidus]**
    <https://pubmed.ncbi.nlm.nih.gov/24130347/>
64. Zeller V, Dzeing-Ella A, Kitzis M-D, et al. **Continuous clindamycin infusion, an innovative approach to treating bone and joint infections.** (bone antibiotic penetration; rifampicin adjunct rationale)
    *Antimicrob Agents Chemother* 2010;54:88-92.
    <https://pubmed.ncbi.nlm.nih.gov/19841145/>

## 9. Pharmacological and device therapy

65. Wieman TJ, Smiell JM, Su Y. **Efficacy and safely of a topical gel formulation of recombinant human platelet-derived growth factor-BB (becaplermin) in patients with chronic neuropathic diabetic ulcers.**
    *Diabetes Care* 1998;21:822-827. **[ANCHOR — 50% vs 35% complete closure at 20 weeks]**
    <https://pubmed.ncbi.nlm.nih.gov/9589248/>
66. Smiell JM, Wieman TJ, Steed DL, et al. **Efficacy and safety of becaplermin in patients with nonhealing, lower extremity diabetic ulcers: a combined analysis of four randomized studies.**
    *Wound Repair Regen* 1999;7:335-346.
    <https://pubmed.ncbi.nlm.nih.gov/10564562/>
67. Ziyadeh N, Fife D, Walker AM, et al. **A matched cohort study of the risk of cancer in users of becaplermin.**
    *Adv Skin Wound Care* 2011;24:31-39. (the boxed-warning malignancy signal in the safety cluster)
    <https://pubmed.ncbi.nlm.nih.gov/21173588/>
68. Edmonds M, Lázaro-Martínez JL, Alfayate-García JM, et al. **Sucrose octasulfate dressing versus control dressing in patients with neuroischaemic diabetic foot ulcers (Explorer): an international, multicentre, double-blind, randomised, controlled trial.**
    *Lancet Diabetes Endocrinol* 2018;6:186-196. **[ANCHOR — 48% vs 30% closure at 20 weeks]**
    <https://pubmed.ncbi.nlm.nih.gov/29275068/>
69. Lázaro-Martínez JL, Edmonds M, Rayman G, et al. **Optimal wound closure of diabetic foot ulcers with early initiation of TLC-NOSF treatment: post-hoc analysis of Explorer.**
    *J Wound Care* 2019;28:358-367.
    <https://pubmed.ncbi.nlm.nih.gov/31166861/>
70. Shukla VK, Rasheed MA, Kumar M, et al. **A trial to disprove the role of esmolol hydrochloride topical gel in diabetic foot ulcer healing** / Kaul S, et al. **Topical esmolol hydrochloride as a novel treatment modality for diabetic foot ulcers: a phase 3 randomized clinical trial.**
    *JAMA Netw Open* 2023;6:e2311509.
    <https://pubmed.ncbi.nlm.nih.gov/37155168/>
71. Veves A, Falanga V, Armstrong DG, Sabolinski ML. **Graftskin, a human skin equivalent, is effective in the management of noninfected neuropathic diabetic foot ulcers: a prospective randomized multicenter clinical trial.**
    *Diabetes Care* 2001;24:290-295.
    <https://pubmed.ncbi.nlm.nih.gov/11213881/>
72. Zelen CM, Serena TE, Denoziere G, Fetterolf DE. **A prospective randomised comparative parallel study of amniotic membrane wound graft in the management of diabetic foot ulcers.**
    *Int Wound J* 2013;10:502-507.
    <https://pubmed.ncbi.nlm.nih.gov/23742102/>
73. Game F, Jeffcoate W, Tarnow L, et al. **LeucoPatch system for the management of hard-to-heal diabetic foot ulcers in the UK, Denmark, and Sweden: an observer-masked, randomised controlled trial.**
    *Lancet Diabetes Endocrinol* 2018;6:870-878.
    <https://pubmed.ncbi.nlm.nih.gov/30243803/>
74. Armstrong DG, Lavery LA; Diabetic Foot Study Consortium. **Negative pressure wound therapy after partial diabetic foot amputation: a multicentre, randomised controlled trial.**
    *Lancet* 2005;366:1704-1710.
    <https://pubmed.ncbi.nlm.nih.gov/16291063/>
75. Löndahl M, Katzman P, Nilsson A, Hammarlund C. **Hyperbaric oxygen therapy facilitates healing of chronic foot ulcers in patients with diabetes.**
    *Diabetes Care* 2010;33:998-1003.
    <https://pubmed.ncbi.nlm.nih.gov/20427683/>
76. Fedorko L, Bowen JM, Jones W, et al. **Hyperbaric oxygen therapy does not reduce indications for amputation in patients with diabetes with nonhealing ulcers of the lower limb: a prospective, double-blind, randomized controlled clinical trial.**
    *Diabetes Care* 2016;39:392-399. (the negative HBOT trial — why the model keeps the HBOT effect modest)
    <https://pubmed.ncbi.nlm.nih.gov/26740639/>
77. Frykberg RG, Franks PJ, Edmonds M, et al. **A multinational, multicenter, randomized, double-blinded, placebo-controlled trial to evaluate the efficacy of cyclical topical wound oxygen (TWO2) therapy.**
    *Diabetes Care* 2020;43:616-624.
    <https://pubmed.ncbi.nlm.nih.gov/31883902/>
78. Rayman G, Vas P, Dhatariya K, et al. **Guidelines on use of interventions to enhance healing of chronic foot ulcers in diabetes (IWGDF 2019/2023).**
    *Diabetes Metab Res Rev* 2020;36(S1):e3283.
    <https://pubmed.ncbi.nlm.nih.gov/32176450/>

## 10. Glycaemic control, systemic therapy and the slow axes

79. Fernando ME, Seneviratne RM, Tan YM, et al. **Intensive versus conventional glycaemic control for treating diabetic foot ulcers.**
    *Cochrane Database Syst Rev* 2016;1:CD010764. **[ANCHOR — weak/uncertain effect of glycaemic control on index-ulcer healing]**
    <https://pubmed.ncbi.nlm.nih.gov/26758576/>
80. Christman AL, Selvin E, Margolis DJ, et al. **Hemoglobin A1c predicts healing rate in diabetic wounds.**
    *J Invest Dermatol* 2011;131:2121-2127.
    <https://pubmed.ncbi.nlm.nih.gov/21697890/>
81. Marso SP, Bain SC, Consoli A, et al. **Semaglutide and cardiovascular outcomes in patients with type 2 diabetes (SUSTAIN-6).**
    *N Engl J Med* 2016;375:1834-1844.
    <https://pubmed.ncbi.nlm.nih.gov/27633186/>
82. Neal B, Perkovic V, Mahaffey KW, et al. **Canagliflozin and cardiovascular and renal events in type 2 diabetes (CANVAS).**
    *N Engl J Med* 2017;377:644-657. (the SGLT2-inhibitor amputation signal in the safety cluster)
    <https://pubmed.ncbi.nlm.nih.gov/28605608/>
83. Nathan DM, Genuth S, Lachin J, et al. **The effect of intensive treatment of diabetes on the development and progression of long-term complications in insulin-dependent diabetes mellitus (DCCT).**
    *N Engl J Med* 1993;329:977-986. (the years-long time constant of the neuropathy axis)
    <https://pubmed.ncbi.nlm.nih.gov/8366922/>

## 11. Recurrence, remission and prevention

84. Bus SA, Sacco ICN, Monteiro-Soares M, et al. **Guidelines on the prevention of foot ulcers in persons with diabetes (IWGDF 2023).**
    *Diabetes Metab Res Rev* 2024;40:e3651.
    <https://pubmed.ncbi.nlm.nih.gov/37302121/>
85. Bus SA, Waaijman R, Arts M, et al. **Effect of custom-made footwear on foot ulcer recurrence in diabetes: a multicenter randomized controlled trial.**
    *Diabetes Care* 2013;36:4109-4116. **[ANCHOR — therapeutic footwear and adherence-dependent recurrence reduction]**
    <https://pubmed.ncbi.nlm.nih.gov/24130357/>
86. Lavery LA, Higgins KR, Lanctot DR, et al. **Preventing diabetic foot ulcer recurrence in high-risk patients: use of temperature monitoring as a self-assessment tool.**
    *Diabetes Care* 2007;30:14-20. **[ANCHOR — home skin-temperature monitoring, ~3-fold recurrence reduction]**
    <https://pubmed.ncbi.nlm.nih.gov/17192326/>
87. Armstrong DG, Holtz-Neiderer K, Wendel C, et al. **Skin temperature monitoring reduces the risk for diabetic foot ulceration in high-risk patients.**
    *Am J Med* 2007;120:1042-1046.
    <https://pubmed.ncbi.nlm.nih.gov/18060924/>
88. Levin ME. **Preventing amputation in the patient with diabetes.**
    *Diabetes Care* 1995;18:1383-1394.
    <https://pubmed.ncbi.nlm.nih.gov/8721944/>
89. Levy BF, Hanft JR. **The role of scar tissue and biomechanical stress in recurrence of neuropathic plantar ulceration.** — see also Bus SA, van Netten JJ. **A shift in priority in diabetic foot care and research: 75% of foot ulcers are preventable.**
    *Diabetes Metab Res Rev* 2016;32(S1):195-200.
    <https://pubmed.ncbi.nlm.nih.gov/26452160/>
90. Boulton AJM. **The diabetic foot: grand overview, epidemiology and pathogenesis.**
    *Diabetes Metab Res Rev* 2008;24(S1):S3-S6.
    <https://pubmed.ncbi.nlm.nih.gov/18442166/>

## 12. Modelling methodology

91. Buganza-Tepole A, Kuhl E. **Computational modeling of chronic wound healing: a review of returning to the wound-healing continuum.**
    *Exp Dermatol* 2013;22:293-298.
    <https://pubmed.ncbi.nlm.nih.gov/23528206/>
92. Menke NB, Cain JW, Reynolds A, et al. **An in silico approach to the analysis of acute wound healing.**
    *Wound Repair Regen* 2010;18:105-113.
    <https://pubmed.ncbi.nlm.nih.gov/20002895/>
93. Ziraldo C, Mi Q, An G, Vodovotz Y. **Computational modeling of inflammation and wound healing.**
    *Adv Wound Care* 2013;2:527-537.
    <https://pubmed.ncbi.nlm.nih.gov/24761338/>
94. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn* 2015;42:S84-S85. — mrgsolve documentation: <https://mrgsolve.org/>
95. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust application of systems pharmacology.**
    *CPT Pharmacometrics Syst Pharmacol* 2016;5:235-249.
    <https://pubmed.ncbi.nlm.nih.gov/27299936/>

---

*Compiled for the QSP Disease Model Library. All model parameters are
illustrative approximations drawn from the public literature above; they have
not been fitted to patient-level data and must not be used for clinical
decision-making.*
