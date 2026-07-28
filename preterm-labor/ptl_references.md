# Preterm Labor / Spontaneous Preterm Birth — Reference Library

Curated evidence base for the QSP model in this directory
(`ptl_qsp_model.dot` · `ptl_mrgsolve_model.R` · `ptl_shiny_app.R`).
References are grouped by the model block they support, so that any parameter
or structural choice can be traced back to a source.

**100 references.** Links are PubMed or, for guidelines and tools, the issuing
body's canonical page.

---

## 1. Epidemiology, definitions and the "syndrome" framing

The single most important conceptual reference for this model is Romero's
argument that preterm birth is not a disease but a **syndrome** — several
distinct pathological processes converging on one common terminal pathway.
The entire architecture of the mechanistic map (four ignition arms, three
effector limbs) is that argument drawn as a graph.

1. Romero R, Dey SK, Fisher SJ. **Preterm labor: one syndrome, many causes.**
   *Science* 2014;345(6198):760-5. <https://pubmed.ncbi.nlm.nih.gov/25124429/>
2. Goldenberg RL, Culhane JF, Iams JD, Romero R. **Epidemiology and causes of
   preterm birth.** *Lancet* 2008;371(9606):75-84.
   <https://pubmed.ncbi.nlm.nih.gov/18177778/>
3. Blencowe H, Cousens S, Oestergaard MZ, et al. **National, regional, and
   worldwide estimates of preterm birth rates in the year 2010.**
   *Lancet* 2012;379(9832):2162-72.
   <https://pubmed.ncbi.nlm.nih.gov/22682464/>
4. Chawanpaiboon S, Vogel JP, Moller AB, et al. **Global, regional, and
   national estimates of levels of preterm birth in 2014.**
   *Lancet Glob Health* 2019;7(1):e37-e46.
   <https://pubmed.ncbi.nlm.nih.gov/30389451/>
5. Ohuma EO, Moller AB, Bradley E, et al. **National, regional, and global
   estimates of preterm birth in 2020.** *Lancet* 2023;402(10409):1261-71.
   <https://pubmed.ncbi.nlm.nih.gov/37805217/>
6. Vogel JP, Chawanpaiboon S, Moller AB, et al. **The global epidemiology of
   preterm birth.** *Best Pract Res Clin Obstet Gynaecol* 2018;52:3-12.
   <https://pubmed.ncbi.nlm.nih.gov/29779863/>
7. Manuck TA, Rice MM, Bailit JL, et al. **Preterm neonatal morbidity and
   mortality by gestational age: a contemporary cohort.**
   *Am J Obstet Gynecol* 2016;215(1):103.e1-103.e14.
   <https://pubmed.ncbi.nlm.nih.gov/26772790/>
8. Stoll BJ, Hansen NI, Bell EF, et al. **Trends in care practices, morbidity,
   and mortality of extremely preterm neonates, 1993-2012.**
   *JAMA* 2015;314(10):1039-51.
   <https://pubmed.ncbi.nlm.nih.gov/26348753/>

---

## 2. Arm 1 — the placental CRH clock

Sets `CRH_TDOUBLE`, `CRH_MAX`, `FB_CORT`, `CLOCK_ADV` and the fetal
adrenal / estriol cascade (`DHEAS`, `CORTF`, `E3`).

9. McLean M, Bisits A, Davies J, et al. **A placental clock controlling the
   length of human pregnancy.** *Nat Med* 1995;1(5):460-3.
   <https://pubmed.ncbi.nlm.nih.gov/7585095/>
10. Challis JR, Matthews SG, Gibb W, Lye SJ. **Endocrine and paracrine
    regulation of birth at term and preterm.** *Endocr Rev* 2000;21(5):514-50.
    <https://pubmed.ncbi.nlm.nih.gov/11041447/>
11. Smith R. **Parturition.** *N Engl J Med* 2007;356(3):271-83.
    <https://pubmed.ncbi.nlm.nih.gov/17229954/>
12. Wadhwa PD, Garite TJ, Porto M, et al. **Placental corticotropin-releasing
    hormone, spontaneous preterm birth, and fetal growth restriction: a
    prospective investigation.** *Am J Obstet Gynecol* 2004;191(4):1063-9.
    <https://pubmed.ncbi.nlm.nih.gov/15507922/>
13. Grammatopoulos DK. **Placental corticotrophin-releasing hormone and its
    receptors in human pregnancy and labour.**
    *J Neuroendocrinol* 2008;20(4):432-8.
    <https://pubmed.ncbi.nlm.nih.gov/18266948/>
14. Mesiano S, Jaffe RB. **Developmental and functional biology of the primate
    fetal adrenal cortex.** *Endocr Rev* 1997;18(3):378-403.
    <https://pubmed.ncbi.nlm.nih.gov/9183569/>

---

## 3. Arm 2 — functional progesterone withdrawal (PR-A / PR-B)

The reason `PRW` is a state variable rather than a function of plasma
progesterone. In humans there is **no systemic progesterone fall** before
labour; withdrawal is a receptor-isoform event inside the myometrium.

15. Mesiano S, Chan EC, Fitter JT, et al. **Progesterone withdrawal and
    estrogen activation in human parturition are coordinated by progesterone
    receptor A expression in the myometrium.**
    *J Clin Endocrinol Metab* 2002;87(6):2924-30.
    <https://pubmed.ncbi.nlm.nih.gov/12050275/>
16. Merlino AA, Welsh TN, Tan H, et al. **Nuclear progesterone receptors in the
    human pregnancy myometrium: evidence that parturition involves functional
    progesterone withdrawal mediated by increased expression of progesterone
    receptor-A.** *J Clin Endocrinol Metab* 2007;92(5):1927-33.
    <https://pubmed.ncbi.nlm.nih.gov/17341556/>
17. Renthal NE, Chen CC, Williams KC, et al. **miR-200 family and targets,
    ZEB1 and ZEB2, modulate uterine quiescence and contractility during
    pregnancy and labor.** *Proc Natl Acad Sci USA* 2010;107(48):20828-33.
    <https://pubmed.ncbi.nlm.nih.gov/21079000/>
18. Zakar T, Hertelendy F. **Progesterone withdrawal: key to parturition.**
    *Am J Obstet Gynecol* 2007;196(4):289-96.
    <https://pubmed.ncbi.nlm.nih.gov/17403397/>
19. Nadeem L, Shynlova O, Matysiak-Zablocki E, et al. **Molecular evidence of
    functional progesterone withdrawal in human myometrium.**
    *Nat Commun* 2016;7:11565.
    <https://pubmed.ncbi.nlm.nih.gov/27220952/>
20. Hardy DB, Janowski BA, Corey DR, Mendelson CR. **Progesterone receptor
    plays a major antiinflammatory role in human myometrial cells by
    antagonism of NF-kB activation of cyclooxygenase 2 expression.**
    *Mol Endocrinol* 2006;20(11):2724-33.
    <https://pubmed.ncbi.nlm.nih.gov/16772530/>
21. Condon JC, Jeyasuria P, Faust JM, et al. **A decline in the levels of
    progesterone receptor coactivators in the pregnant uterus at term may
    antagonize progesterone receptor function and contribute to the initiation
    of parturition.** *Proc Natl Acad Sci USA* 2003;100(16):9518-23.
    <https://pubmed.ncbi.nlm.nih.gov/12886011/>

---

## 4. Arm 3 — intra-amniotic infection, sterile inflammation and FIRS

Sets `BACT`, `NFKB`, `IL6`, `FIRS`, `IL6_THRESH` and `PR_ON_NFKB`.
Note reference 24: **sterile** intra-amniotic inflammation is roughly twice as
common as culture-proven MIAC, which is why `STERILE_DRIVE` is a separate
parameter and is non-zero by default.

22. Romero R, Miranda J, Chaiworapongsa T, et al. **Prevalence and clinical
    significance of sterile intra-amniotic inflammation in patients with
    preterm labor and intact membranes.**
    *Am J Reprod Immunol* 2014;72(5):458-74.
    <https://pubmed.ncbi.nlm.nih.gov/25078709/>
23. Romero R, Espinoza J, Kusanovic JP, et al. **The preterm parturition
    syndrome.** *BJOG* 2006;113(Suppl 3):17-42.
    <https://pubmed.ncbi.nlm.nih.gov/17206962/>
24. Romero R, Miranda J, Chaemsaithong P, et al. **Sterile and microbial-
    associated intra-amniotic inflammation in preterm prelabor rupture of
    membranes.** *J Matern Fetal Neonatal Med* 2015;28(12):1394-409.
    <https://pubmed.ncbi.nlm.nih.gov/25190175/>
25. Gomez R, Romero R, Ghezzi F, et al. **The fetal inflammatory response
    syndrome.** *Am J Obstet Gynecol* 1998;179(1):194-202.
    <https://pubmed.ncbi.nlm.nih.gov/9704787/>
26. Yoon BH, Romero R, Park JS, et al. **Fetal exposure to an intra-amniotic
    inflammation and the development of cerebral palsy at the age of three
    years.** *Am J Obstet Gynecol* 2000;182(3):675-81.
    <https://pubmed.ncbi.nlm.nih.gov/10739529/>
27. Combs CA, Gravett M, Garite TJ, et al. **Amniotic fluid infection,
    inflammation, and colonization in preterm labor with intact membranes.**
    *Am J Obstet Gynecol* 2014;210(2):125.e1-125.e15.
    <https://pubmed.ncbi.nlm.nih.gov/24274987/>
28. DiGiulio DB, Callahan BJ, McMurdie PJ, et al. **Temporal and spatial
    variation of the human microbiota during pregnancy.**
    *Proc Natl Acad Sci USA* 2015;112(35):11060-5.
    <https://pubmed.ncbi.nlm.nih.gov/26283357/>
29. Fettweis JM, Serrano MG, Brooks JP, et al. **The vaginal microbiome and
    preterm birth.** *Nat Med* 2019;25(6):1012-21.
    <https://pubmed.ncbi.nlm.nih.gov/31142849/>
30. Sweeney EL, Dando SJ, Kallapur SG, Knox CL. **The human Ureaplasma species
    as causative agents of chorioamnionitis.**
    *Clin Microbiol Rev* 2016;30(1):349-79.
    <https://pubmed.ncbi.nlm.nih.gov/27974410/>
31. Gomez-Lopez N, Romero R, Xu Y, et al. **Are amniotic fluid neutrophils in
    women with intraamniotic infection and/or inflammation of fetal or
    maternal origin?** *Am J Obstet Gynecol* 2017;217(6):693.e1-693.e16.
    <https://pubmed.ncbi.nlm.nih.gov/28964823/>
32. Gomez-Lopez N, Motomura K, Miller D, et al. **Inflammasomes: their role in
    normal and complicated pregnancies.**
    *J Immunol* 2019;203(11):2757-69.
    <https://pubmed.ncbi.nlm.nih.gov/31740550/>

---

## 5. Arm 4 — decidual activation, senescence, thrombin and overdistension

Sets `THROMBIN_DRV`, `STRETCH`, the SASP/DAMP term and the PAR-1 uterotonic
contribution to `CAI`.

33. Menon R, Behnia F, Polettini J, et al. **Placental membrane aging and HMGB1
    signaling associated with human parturition.**
    *Aging (Albany NY)* 2016;8(2):216-30.
    <https://pubmed.ncbi.nlm.nih.gov/26851389/>
34. Behnia F, Taylor BD, Woodson M, et al. **Chorioamniotic membrane
    senescence: a signal for parturition?**
    *Am J Obstet Gynecol* 2015;213(3):359.e1-359.e16.
    <https://pubmed.ncbi.nlm.nih.gov/26025293/>
35. Elovitz MA, Baron J, Phillippe M. **The role of thrombin in preterm
    parturition.** *Am J Obstet Gynecol* 2001;185(5):1059-63.
    <https://pubmed.ncbi.nlm.nih.gov/11717633/>
36. Rosen T, Kuczynski E, O'Neill LM, et al. **Plasma levels of
    thrombin-antithrombin complexes predict preterm premature rupture of the
    fetal membranes.** *J Matern Fetal Med* 2001;10(5):297-300.
    <https://pubmed.ncbi.nlm.nih.gov/11730490/>
37. Lyall F, Lye S, Teoh T, et al. **Expression of Gsalpha, connexin-43,
    connexin-26, and EP1, 3, and 4 receptors in myometrium of prelabor
    singleton versus multiple gestations and the effects of mechanical
    stretch and steroids on Gsalpha.**
    *J Soc Gynecol Investig* 2002;9(5):299-307.
    <https://pubmed.ncbi.nlm.nih.gov/12383914/>
38. Shynlova O, Tsui P, Jaffer S, Lye SJ. **Integration of endocrine and
    mechanical signals in the regulation of myometrial functions during
    pregnancy and labour.**
    *Eur J Obstet Gynecol Reprod Biol* 2009;144(Suppl 1):S2-10.
    <https://pubmed.ncbi.nlm.nih.gov/19299064/>

---

## 6. Prostaglandins, the 15-PGDH barrier and the CAP programme

Sets `COX2`, `PG`, `PGDH_FRAC`, `OXTR`, `CX43` and their maxima.

39. Olson DM. **The role of prostaglandins in the initiation of parturition.**
    *Best Pract Res Clin Obstet Gynaecol* 2003;17(5):717-30.
    <https://pubmed.ncbi.nlm.nih.gov/12972010/>
40. Challis JR, Sloboda DM, Alfaidy N, et al. **Prostaglandins and mechanisms
    of preterm birth.** *Reproduction* 2002;124(1):1-17.
    <https://pubmed.ncbi.nlm.nih.gov/12090913/>
41. Fuchs AR, Fuchs F, Husslein P, Soloff MS. **Oxytocin receptors in the human
    uterus during pregnancy and parturition.**
    *Am J Obstet Gynecol* 1984;150(6):734-41.
    <https://pubmed.ncbi.nlm.nih.gov/6093521/>
42. Kimura T, Takemura M, Nomura S, et al. **Expression of oxytocin receptor in
    human pregnant myometrium.** *Endocrinology* 1996;137(2):780-5.
    <https://pubmed.ncbi.nlm.nih.gov/8593830/>
43. Chow L, Lye SJ. **Expression of the gap junction protein connexin-43 is
    increased in the human myometrium toward term and with the onset of
    labor.** *Am J Obstet Gynecol* 1994;170(3):788-95.
    <https://pubmed.ncbi.nlm.nih.gov/8141201/>
44. Wray S, Prendergast C. **The myometrium: from excitation to contractions
    and labour.** *Adv Exp Med Biol* 2019;1124:233-63.
    <https://pubmed.ncbi.nlm.nih.gov/31183830/>
45. Aguilar HN, Mitchell BF. **Physiological pathways and molecular mechanisms
    regulating uterine contractility.**
    *Hum Reprod Update* 2010;16(6):725-44.
    <https://pubmed.ncbi.nlm.nih.gov/20551073/>

---

## 7. Cervical remodelling and cervical length

Sets `COLL`, `HA`, `MMP`, `CLEN`, `CL_REF`, `CL_STEEP` and the `f_cerv` hazard
limb. Reference 46 is the source of the model's cervical-length distribution
and of the continuous, non-linear CL-to-risk relation.

46. Iams JD, Goldenberg RL, Meis PJ, et al. **The length of the cervix and the
    risk of spontaneous premature delivery.** *N Engl J Med* 1996;334(9):567-72.
    <https://pubmed.ncbi.nlm.nih.gov/8569824/>
47. Word RA, Li XH, Hnat M, Carrick K. **Dynamics of cervical remodeling during
    pregnancy and parturition: mechanisms and current concepts.**
    *Semin Reprod Med* 2007;25(1):69-79.
    <https://pubmed.ncbi.nlm.nih.gov/17205425/>
48. Timmons B, Akins M, Mahendroo M. **Cervical remodeling during pregnancy and
    parturition.** *Trends Endocrinol Metab* 2010;21(6):353-61.
    <https://pubmed.ncbi.nlm.nih.gov/20172738/>
49. Vink J, Feltovich H. **Cervical etiology of spontaneous preterm birth.**
    *Semin Fetal Neonatal Med* 2016;21(2):106-12.
    <https://pubmed.ncbi.nlm.nih.gov/26776146/>
50. Myers KM, Feltovich H, Mazza E, et al. **The mechanical role of the cervix
    in pregnancy.** *J Biomech* 2015;48(9):1511-23.
    <https://pubmed.ncbi.nlm.nih.gov/25841293/>

---

## 8. Fetal membranes and PPROM

Sets `MEMB`, `KDEG_MEMB`, the ZAM weak-zone concept and the `f_memb` limb.

51. Menon R, Richardson LS. **Preterm prelabor rupture of the membranes: a
    disease of the fetal membranes.**
    *Semin Perinatol* 2017;41(7):409-19.
    <https://pubmed.ncbi.nlm.nih.gov/28912049/>
52. Parry S, Strauss JF 3rd. **Premature rupture of the fetal membranes.**
    *N Engl J Med* 1998;338(10):663-70.
    <https://pubmed.ncbi.nlm.nih.gov/9486996/>
53. Kumar D, Moore RM, Mercer BM, et al. **The physiology of fetal membrane
    weakening and rupture: insights gained from the determination of physical
    properties revisited.** *Placenta* 2016;42:59-73.
    <https://pubmed.ncbi.nlm.nih.gov/27238715/>

---

## 9. Biomarkers and prediction

Sets `FFN`, `FFN_POS` and the QUiPP-style composite risk readout.

54. Lockwood CJ, Senyei AE, Dische MR, et al. **Fetal fibronectin in cervical
    and vaginal secretions as a predictor of preterm delivery.**
    *N Engl J Med* 1991;325(10):669-74.
    <https://pubmed.ncbi.nlm.nih.gov/1870640/>
55. Berghella V, Saccone G. **Fetal fibronectin testing for prevention of
    preterm birth in singleton pregnancies with threatened preterm labor.**
    *Cochrane Database Syst Rev* 2019;7:CD006843.
    <https://pubmed.ncbi.nlm.nih.gov/31356681/>
56. Watson HA, Carter J, Seed PT, et al. **The QUiPP App: a safe alternative to
    a treat-all strategy for threatened preterm labor.**
    *Ultrasound Obstet Gynecol* 2017;50(3):342-6.
    <https://pubmed.ncbi.nlm.nih.gov/28067433/>
57. Ngo TTM, Moufarrej MN, Rasmussen MH, et al. **Noninvasive blood tests for
    fetal development predict gestational age and preterm delivery.**
    *Science* 2018;360(6393):1133-6.
    <https://pubmed.ncbi.nlm.nih.gov/29880692/>

---

## 10. Progestogen prevention

Sets `EMAX_P4_PRW`, `EC50_P4_LOC`, `FUP_RATIO`, `OHPC_LOCAL`. Note the
deliberate tension between references 58 and 60 — the model encodes it rather
than resolving it.

58. Meis PJ, Klebanoff M, Thom E, et al. **Prevention of recurrent preterm
    delivery by 17 alpha-hydroxyprogesterone caproate.**
    *N Engl J Med* 2003;348(24):2379-85.
    <https://pubmed.ncbi.nlm.nih.gov/12802023/>
59. Hassan SS, Romero R, Vidyadhari D, et al. **Vaginal progesterone reduces
    the rate of preterm birth in women with a sonographic short cervix: a
    multicenter, randomized, double-blind, placebo-controlled trial.**
    *Ultrasound Obstet Gynecol* 2011;38(1):18-31.
    <https://pubmed.ncbi.nlm.nih.gov/21472815/>
60. Blackwell SC, Gyamfi-Bannerman C, Biggio JR Jr, et al. **17-OHPC to prevent
    recurrent preterm birth in singleton gestations (PROLONG study).**
    *Am J Perinatol* 2020;37(2):127-36.
    <https://pubmed.ncbi.nlm.nih.gov/31600795/>
61. EPPPIC Group. **Evaluating Progestogens for Preventing Preterm birth
    International Collaborative (EPPPIC): meta-analysis of individual
    participant data from randomised controlled trials.**
    *Lancet* 2021;397(10280):1183-94.
    <https://pubmed.ncbi.nlm.nih.gov/33773630/>
62. Norman JE, Marlow N, Messow CM, et al. **Vaginal progesterone prophylaxis
    for preterm birth (the OPPTIMUM study): a multicentre, randomised,
    double-blind trial.** *Lancet* 2016;387(10033):2106-16.
    <https://pubmed.ncbi.nlm.nih.gov/26921136/>
63. Cicinelli E, de Ziegler D, Bulletti C, et al. **Direct transport of
    progesterone from vagina to uterus.**
    *Obstet Gynecol* 2000;95(3):403-6.
    <https://pubmed.ncbi.nlm.nih.gov/10711552/>
64. Caritis SN, Sharma S, Venkataramanan R, et al. **Pharmacokinetics of
    17-hydroxyprogesterone caproate in multifetal gestation.**
    *Am J Obstet Gynecol* 2011;205(1):40.e1-8.
    <https://pubmed.ncbi.nlm.nih.gov/21620357/>

---

## 11. Tocolysis — pharmacology and the efficacy ceiling

The core evidence behind the model's central structural claim: tocolytics
delay birth by 48 hours to 7 days and do **not** improve neonatal outcome.

65. Flenady V, Wojcieszek AM, Papatsonis DN, et al. **Calcium channel blockers
    for inhibiting preterm labour and birth.**
    *Cochrane Database Syst Rev* 2014;(6):CD002255.
    <https://pubmed.ncbi.nlm.nih.gov/24901312/>
66. Flenady V, Reinebrant HE, Liley HG, et al. **Oxytocin receptor antagonists
    for inhibiting preterm labour.**
    *Cochrane Database Syst Rev* 2014;(6):CD004452.
    <https://pubmed.ncbi.nlm.nih.gov/24903678/>
67. Worldwide Atosiban versus Beta-agonists Study Group. **Effectiveness and
    safety of the oxytocin antagonist atosiban versus beta-adrenergic agonists
    in the treatment of preterm labour.** *BJOG* 2001;108(2):133-42.
    <https://pubmed.ncbi.nlm.nih.gov/11236112/>
68. Goodwin TM, Millar L, North L, et al. **The pharmacokinetics of the
    oxytocin antagonist atosiban in pregnant women with preterm uterine
    contractions.** *Am J Obstet Gynecol* 1995;173(3 Pt 1):913-7.
    <https://pubmed.ncbi.nlm.nih.gov/7573267/>
69. Haas DM, Caldwell DM, Kirkpatrick P, et al. **Tocolytic therapy for preterm
    delivery: systematic review and network meta-analysis.**
    *BMJ* 2012;345:e6226.
    <https://pubmed.ncbi.nlm.nih.gov/23048010/>
70. Moise KJ Jr, Huhta JC, Sharif DS, et al. **Indomethacin in the treatment of
    premature labor: effects on the fetal ductus arteriosus.**
    *N Engl J Med* 1988;319(6):327-31.
    <https://pubmed.ncbi.nlm.nih.gov/3292537/>
71. Reinebrant HE, Pileggi-Castro C, Romero CL, et al. **Cyclo-oxygenase (COX)
    inhibitors for treating preterm labour.**
    *Cochrane Database Syst Rev* 2015;(6):CD001992.
    <https://pubmed.ncbi.nlm.nih.gov/26057712/>
72. Neilson JP, West HM, Dowswell T. **Betamimetics for inhibiting preterm
    labour.** *Cochrane Database Syst Rev* 2014;(2):CD004352.
    <https://pubmed.ncbi.nlm.nih.gov/24578236/>
73. Crowther CA, Brown J, McKinlay CJ, Middleton P. **Magnesium sulphate for
    preventing preterm birth in threatened preterm labour.**
    *Cochrane Database Syst Rev* 2014;(8):CD001060.
    <https://pubmed.ncbi.nlm.nih.gov/25126773/>
74. Vogel JP, Nardin JM, Dowswell T, et al. **Combination of tocolytic agents
    for inhibiting preterm labour.**
    *Cochrane Database Syst Rev* 2014;(7):CD006169.
    <https://pubmed.ncbi.nlm.nih.gov/25010869/>

---

## 12. Antenatal corticosteroids — the intervention that changes outcome

Sets `EMAX_SURF`, `EC50_BETF`, `KSURF_ACS`, `ACS_DECAY` and the RDS logistic.

75. Liggins GC, Howie RN. **A controlled trial of antepartum glucocorticoid
    treatment for prevention of the respiratory distress syndrome in premature
    infants.** *Pediatrics* 1972;50(4):515-25.
    <https://pubmed.ncbi.nlm.nih.gov/4561295/>
76. Roberts D, Brown J, Medley N, Dalziel SR. **Antenatal corticosteroids for
    accelerating fetal lung maturation for women at risk of preterm birth.**
    *Cochrane Database Syst Rev* 2020;12:CD004454.
    <https://pubmed.ncbi.nlm.nih.gov/33368142/>
77. Gyamfi-Bannerman C, Thom EA, Blackwell SC, et al. **Antenatal betamethasone
    for women at risk for late preterm delivery (ALPS).**
    *N Engl J Med* 2016;374(14):1311-20.
    <https://pubmed.ncbi.nlm.nih.gov/26842679/>
78. WHO ACTION Trials Collaborators. **Antenatal dexamethasone for early
    preterm birth in low-resource countries.**
    *N Engl J Med* 2020;383(26):2514-25.
    <https://pubmed.ncbi.nlm.nih.gov/33095526/>
79. Ballard PL, Ballard RA. **Scientific basis and therapeutic regimens for use
    of antenatal glucocorticoids.**
    *Am J Obstet Gynecol* 1995;173(1):254-62.
    <https://pubmed.ncbi.nlm.nih.gov/7631700/>
80. Crowther CA, McKinlay CJ, Middleton P, Harding JE. **Repeat doses of
    prenatal corticosteroids for women at risk of preterm birth.**
    *Cochrane Database Syst Rev* 2015;(7):CD003935.
    <https://pubmed.ncbi.nlm.nih.gov/26142898/>

---

## 13. Magnesium sulfate for fetal neuroprotection

A separate indication from tocolysis, with a separate exposure-response.
Sets `EMAX_MGN`, `EC50_MGN`, `KMG_BRAIN`.

81. Rouse DJ, Hirtz DG, Thom E, et al. **A randomized, controlled trial of
    magnesium sulfate for the prevention of cerebral palsy (BEAM).**
    *N Engl J Med* 2008;359(9):895-905.
    <https://pubmed.ncbi.nlm.nih.gov/18753646/>
82. Doyle LW, Crowther CA, Middleton P, et al. **Magnesium sulphate for women
    at risk of preterm birth for neuroprotection of the fetus.**
    *Cochrane Database Syst Rev* 2009;(1):CD004661.
    <https://pubmed.ncbi.nlm.nih.gov/19160238/>
83. Crowther CA, Hiller JE, Doyle LW, Haslam RR. **Effect of magnesium sulfate
    given for neuroprotection before preterm birth (ACTOMgSO4).**
    *JAMA* 2003;290(20):2669-76.
    <https://pubmed.ncbi.nlm.nih.gov/14645308/>
84. Marret S, Marpeau L, Zupan-Simunek V, et al. **Magnesium sulphate given
    before very-preterm birth to protect infant brain: the randomised
    controlled PREMAG trial.** *BJOG* 2007;114(3):310-8.
    <https://pubmed.ncbi.nlm.nih.gov/17169012/>
85. Lu JF, Nightingale CH. **Magnesium sulfate in eclampsia and pre-eclampsia:
    pharmacokinetic principles.**
    *Clin Pharmacokinet* 2000;38(4):305-14.
    <https://pubmed.ncbi.nlm.nih.gov/10803454/>

---

## 14. Antibiotics, PPROM management and mechanical adjuncts

86. Kenyon SL, Taylor DJ, Tarnow-Mordi W. **Broad-spectrum antibiotics for
    preterm, prelabour rupture of fetal membranes: the ORACLE I randomised
    trial.** *Lancet* 2001;357(9261):979-88.
    <https://pubmed.ncbi.nlm.nih.gov/11293640/>
87. Kenyon S, Boulvain M, Neilson JP. **Antibiotics for preterm rupture of
    membranes.** *Cochrane Database Syst Rev* 2013;(12):CD001058.
    <https://pubmed.ncbi.nlm.nih.gov/24297389/>
88. Berghella V, Rafael TJ, Szychowski JM, et al. **Cerclage for short cervix on
    ultrasonography in women with singleton gestations and previous preterm
    birth: a meta-analysis.** *Obstet Gynecol* 2011;117(3):663-71.
    <https://pubmed.ncbi.nlm.nih.gov/21446209/>
89. Goya M, Pratcorona L, Merced C, et al. **Cervical pessary in pregnant women
    with a short cervix (PECEP): an open-label randomised controlled trial.**
    *Lancet* 2012;379(9828):1800-6.
    <https://pubmed.ncbi.nlm.nih.gov/22475493/>
90. Nicolaides KH, Syngelaki A, Poon LC, et al. **A randomized trial of a
    cervical pessary to prevent preterm singleton birth.**
    *N Engl J Med* 2016;374(11):1044-52.
    <https://pubmed.ncbi.nlm.nih.gov/26981934/>
91. Conde-Agudelo A, Romero R, Nicolaides KH. **Cervical pessary to prevent
    preterm birth in asymptomatic high-risk women: a systematic review and
    meta-analysis.** *Am J Obstet Gynecol* 2020;223(1):42-65.e2.
    <https://pubmed.ncbi.nlm.nih.gov/31932170/>

---

## 15. Guidelines and clinical practice statements

92. American College of Obstetricians and Gynecologists. **Practice Bulletin
    No. 234: Prediction and Prevention of Spontaneous Preterm Birth.**
    *Obstet Gynecol* 2021;138(2):e65-e90.
    <https://pubmed.ncbi.nlm.nih.gov/34293771/>
93. American College of Obstetricians and Gynecologists. **Practice Bulletin
    No. 217: Prelabor Rupture of Membranes.**
    *Obstet Gynecol* 2020;135(3):e80-e97.
    <https://pubmed.ncbi.nlm.nih.gov/32080050/>
94. National Institute for Health and Care Excellence. **Preterm labour and
    birth (NG25).** 2015, updated 2022.
    <https://www.nice.org.uk/guidance/ng25>
95. WHO. **WHO recommendations on interventions to improve preterm birth
    outcomes.** Geneva: World Health Organization, 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26447264/>

---

## 16. QSP methodology and modelling of parturition

96. Bhattacharya S, Zhang Q, Andersen ME. **A deterministic map of Waddington's
    epigenetic landscape for cell fate decision** — used here only as a
    methodological reference for bistable-switch framing of the
    quiescence-to-activation transition.
    *BMC Syst Biol* 2011;5:85.
    <https://pubmed.ncbi.nlm.nih.gov/21619617/>
97. Maner WL, Garfield RE. **Identification of human term and preterm labor
    using artificial neural networks on uterine electromyography data.**
    *Ann Biomed Eng* 2007;35(3):465-73.
    <https://pubmed.ncbi.nlm.nih.gov/17111212/>
98. Vink J, Myers K. **Cervical alterations in pregnancy.**
    *Best Pract Res Clin Obstet Gynaecol* 2018;52:88-102.
    <https://pubmed.ncbi.nlm.nih.gov/30314740/>
99. Baron KT, Gastonguay MR. **mrgsolve: Simulate from ODE-based population
    PK/PD and systems pharmacology models.** R package.
    <https://mrgsolve.org/>
100. Bhattaram VA, Bonapace C, Chilukuri DM, et al. **Impact of pharmacometric
     reviews on new drug approval and labeling decisions** — background on
     model-informed drug development, the regulatory frame this library sits
     in. *Clin Pharmacol Ther* 2007;81(2):213-21.
     <https://pubmed.ncbi.nlm.nih.gov/17235330/>

---

## How the evidence maps onto the model

| Model element | Primary references |
|---|---|
| Four-arm / one-terminal-pathway architecture | 1, 2, 23 |
| `CRH`, `CRH_TDOUBLE`, `FB_CORT` | 9, 10, 11, 12, 13 |
| `PRW` as a state (not plasma P4) | 15, 16, 17, 18, 19 |
| `PR_ON_NFKB` trans-repression | 20, 21 |
| `BACT`, `NFKB`, `IL6`, `IL6_THRESH` | 22, 23, 24, 27, 30 |
| `FIRS` and its outcome coupling | 25, 26, 31, 32 |
| `THROMBIN_DRV`, senescence/SASP | 33, 34, 35, 36 |
| `STRETCH` mechanotransduction | 37, 38 |
| `COX2`, `PG`, `PGDH_FRAC` | 39, 40 |
| `OXTR_MAX`, `CX43_MAX` | 41, 42, 43 |
| E-C coupling (`CAI`, `CONTR`) | 44, 45 |
| `CLEN`, `CL_REF`, `f_cerv` | 46, 47, 48, 49, 50 |
| `MEMB`, `f_memb`, PPROM | 51, 52, 53 |
| `FFN`, `FFN_POS`, composite risk | 54, 55, 56, 57 |
| Vaginal P4 PD and `FUP_RATIO` | 59, 61, 62, 63 |
| 17-OHPC and the MEIS/PROLONG tension | 58, 60, 64 |
| Atosiban PK and `KI_ATO` | 67, 68 |
| Tocolytic efficacy ceiling (48 h, no outcome shift) | 65, 66, 69, 71, 72, 73, 74 |
| Indomethacin `DUCT_GA50` | 70, 71 |
| `EMAX_SURF`, ACS window, RDS logistic | 75, 76, 78, 79, 80 |
| Late-preterm steroids and hypoglycaemia | 77 |
| `EMAX_MGN`, `EC50_MGN`, CP logistic | 81, 82, 83, 84 |
| Magnesium PK and CrCl scaling | 85 |
| Erythromycin latency, co-amoxiclav NEC | 86, 87 |
| Cerclage / pessary offsets | 88, 89, 90, 91 |
| Neonatal outcome logistics vs GA | 7, 8 |

---

*Compiled for the QSP Disease Model Library. Research and education only —
not clinical guidance. Where a reference supports a parameter value, the value
is an interpretation of the published data for simulation purposes and has not
been independently validated.*
