# Primary Postpartum Haemorrhage — References

Literature underpinning the QSP model in this directory (`pph_qsp_model.dot`,
`pph_mrgsolve_model.R`, `pph_shiny_app.R`). Grouped by the part of the model each
source constrains, so that any parameter can be traced back to its evidence.

Links are PubMed (`https://pubmed.ncbi.nlm.nih.gov/<PMID>/`). Every PMID below
was resolved against the PubMed E-utilities API while writing this file, so the
author / title / journal / year shown is the record's own metadata rather than a
recollection. Where a source is a guideline without a PMID the canonical URL is
given.

---

## 1. Epidemiology and the clinical shape of the disease

1. Say L, Chou D, Gemmill A, et al. **Global causes of maternal death: a WHO
   systematic analysis.** Lancet Glob Health 2014;2:e323-33.
   Haemorrhage is the single largest direct cause of maternal death worldwide —
   the reason this model exists.
   <https://pubmed.ncbi.nlm.nih.gov/25103301/>
2. Cresswell JA, Alexander M, Chong MYC, et al. **Global and regional causes of
   maternal deaths 2009-20: a WHO systematic analysis.** Lancet Glob Health 2025.
   <https://pubmed.ncbi.nlm.nih.gov/40064189/>
3. Bateman BT, Berman MF, Riley LE, Leffert LR. **The epidemiology of postpartum
   hemorrhage in a large, nationwide sample of deliveries.** Anesth Analg
   2010;110:1368-73. Atony dominates the aetiological mix — the justification for
   building the model around a single contractile-capacity parameter.
   <https://pubmed.ncbi.nlm.nih.gov/20237047/>
4. Kramer MS, Berg C, Abenhaim H, et al. **Incidence, risk factors, and temporal
   trends in severe postpartum hemorrhage.** Am J Obstet Gynecol 2013;209:449.e1-7.
   <https://pubmed.ncbi.nlm.nih.gov/23871950/>
5. Callaghan WM, Kuklina EV, Berg CJ. **Trends in postpartum hemorrhage: United
   States, 1994-2006.** Am J Obstet Gynecol 2010;202:353.e1-6.
   <https://pubmed.ncbi.nlm.nih.gov/20350642/>
6. Mehrabadi A, Hutcheon JA, Lee L, et al. **Temporal trends in postpartum
   hemorrhage and severe postpartum hemorrhage in Canada from 2003 to 2010.**
   J Obstet Gynaecol Can 2014;36:21-33.
   <https://pubmed.ncbi.nlm.nih.gov/24444284/>
7. Bateman BT, Tsen LC, Liu J, Butwick AJ, Huybrechts KF. **Patterns of
   second-line uterotonic use in a large sample of hospitalizations for childbirth
   in the United States: 2007-2011.** Anesth Analg 2014;119:1344-9.
   Real-world escalation practice, against which the model's ladder is checked.
   <https://pubmed.ncbi.nlm.nih.gov/25166464/>

## 2. Uterine blood flow and the volume the model bleeds from

8. Palmer SK, Zamudio S, Coffin C, et al. **Quantitative estimation of human
   uterine artery blood flow and pelvic blood flow redistribution in pregnancy.**
   Obstet Gynecol 1992;80:1000-6. The primary constraint on `UBF0` (750 mL/min).
   <https://pubmed.ncbi.nlm.nih.gov/1448242/>
9. Hytten F. **Blood volume changes in normal pregnancy.** Clin Haematol
   1985;14:601-12. The ~2100 mL expansion that the model converts into 2.8 minutes.
   <https://pubmed.ncbi.nlm.nih.gov/4075604/>
10. Pritchard JA. **Changes in the blood volume during pregnancy and delivery.**
    Anesthesiology 1965;26:393-9. Term blood volume ~100 mL/kg and the volume
    actually lost at delivery.
    <https://pubmed.ncbi.nlm.nih.gov/14313451/>
11. Pritchard JA, Wiggins KM, Dickey JC. **Blood volume changes in pregnancy and
    the puerperium.** Am J Obstet Gynecol 1960;80:956-63.
    <https://pubmed.ncbi.nlm.nih.gov/13738089/>

## 3. How much blood is actually lost (and how badly it is measured)

12. Bose P, Regan F, Paterson-Brown S. **Improving the accuracy of estimated blood
    loss at obstetric haemorrhage using clinical reconstructions.** BJOG
    2006;113:919-24. Visual estimation understates loss by 30-50%; the reason the
    model reports a continuous loss integral rather than a category.
    <https://pubmed.ncbi.nlm.nih.gov/16907938/>
13. Sloan NL, Durocher J, Aldrich T, Blum J, Winikoff B. **What measured blood loss
    tells us about postpartum bleeding: a systematic review.** BJOG 2010;117:788-800.
    <https://pubmed.ncbi.nlm.nih.gov/20406227/>
14. Diaz V, Abalos E, Carroli G. **Methods for blood loss estimation after vaginal
    birth.** Cochrane Database Syst Rev 2018;9:CD010980.
    <https://pubmed.ncbi.nlm.nih.gov/30211952/>

## 4. Haemodynamic read-outs: shock index and the compensation cliff

15. Le Bas A, Chandraharan E, Addei A, Arulkumaran S. **Use of the "obstetric shock
    index" as an adjunct in identifying significant blood loss in patients with
    massive postpartum hemorrhage.** Int J Gynaecol Obstet 2014;124:253-5.
    Calibrates the model's `ShockIx` thresholds.
    <https://pubmed.ncbi.nlm.nih.gov/24373705/>
16. Nathan HL, Seed PT, Hezelgrave NL, et al. **Shock index thresholds to predict
    adverse outcomes in maternal hemorrhage and sepsis.** Acta Obstet Gynecol Scand
    2019;98:1178-86.
    <https://pubmed.ncbi.nlm.nih.gov/31001814/>
17. El Ayadi AM, Nathan HL, Seed PT, et al. **Vital sign prediction of adverse
    maternal outcomes in women with hypovolemic shock: the role of shock index.**
    PLoS One 2016;11:e0148729.
    <https://pubmed.ncbi.nlm.nih.gov/26901161/>
18. Siegel JH, Fabian M, Smith JA, et al. **Oxygen debt criteria quantify the
    effectiveness of early partial resuscitation after hypovolemic hemorrhagic
    shock.** J Trauma 2003;54:862-80. The basis of the model's mortality
    criterion (accumulated oxygen debt, 120 mL O₂/kg) rather than a volume rule.
    <https://pubmed.ncbi.nlm.nih.gov/12777899/>
19. Bilkovski RN, Rivers EP, Horst HM. **Targeted resuscitation strategies after
    injury.** Curr Opin Crit Care 2004;10:529-38.
    <https://pubmed.ncbi.nlm.nih.gov/15616397/>

## 5. The oxytocin receptor: desensitisation and its consequences

20. Phaneuf S, Rodríguez Liñares B, TambyRaja RL, MacKenzie IZ, López Bernal A.
    **Loss of myometrial oxytocin receptors during oxytocin-induced and
    oxytocin-augmented labour.** J Reprod Fertil 2000;120:91-7. The primary
    constraint on `KDES`/`KREC` and on the OTR ≈ 0.39 target after augmentation.
    <https://pubmed.ncbi.nlm.nih.gov/11006150/>
21. Balki M, Ronayne M, Davies S, et al. **Minimum oxytocin dose requirement after
    cesarean delivery for labor arrest.** Obstet Gynecol 2006;107:45-50. The
    ED90 for oxytocin is several-fold higher after labour augmentation — the
    clinical signature of the EC50 right-shift the model implements.
    <https://pubmed.ncbi.nlm.nih.gov/16394038/>
22. Magalhaes JK, Carvalho JC, Parkes RK, Kingdom J, Li Y, Balki M. **Oxytocin
    pretreatment decreases oxytocin-induced myometrial contractions in pregnant
    rats in a concentration-dependent but not time-dependent manner.** Reprod Sci
    2009;16:501-8.
    <https://pubmed.ncbi.nlm.nih.gov/19164477/>
23. Grotegut CA, Paglia MJ, Johnson LNC, Thames B, James AH. **Oxytocin exposure
    during labor among women with postpartum hemorrhage secondary to uterine
    atony.** Am J Obstet Gynecol 2011;204:56.e1-6. The epidemiological counterpart:
    cumulative oxytocin exposure predicts atonic PPH.
    <https://pubmed.ncbi.nlm.nih.gov/21047614/>
24. Rydén G, Sjöholm I. **Half-life of oxytocin in blood of pregnant and
    non-pregnant women.** Acta Endocrinol (Copenh) 1969;61:425-31. Source of the
    ~4 min half-life used for `CL_OXY`.
    <https://pubmed.ncbi.nlm.nih.gov/5820054/>
25. Rydén G, Sjöholm I. **The metabolism of oxytocin in pregnant and non-pregnant
    women.** Acta Obstet Gynecol Scand Suppl 1971;9:37-8.
    <https://pubmed.ncbi.nlm.nih.gov/5287107/>

## 6. Uterotonic pharmacology and comparative trials

26. Widmer M, Piaggio G, Nguyen TMH, et al. (CHAMPION Trial Group). **Heat-stable
    carbetocin versus oxytocin to prevent hemorrhage after vaginal birth.** N Engl
    J Med 2018;379:743-52. Non-inferiority anchor for the carbetocin arm.
    <https://pubmed.ncbi.nlm.nih.gov/29949473/>
27. Rath W. **Prevention of postpartum haemorrhage with the oxytocin analogue
    carbetocin.** Eur J Obstet Gynecol Reprod Biol 2009;147:15-20. Source of the
    ~40 min half-life used for `KE_CBT`.
    <https://pubmed.ncbi.nlm.nih.gov/19616358/>
28. Vernekar SS, Bhandari S, Datta S, et al. **Effect of heat-stable carbetocin vs
    oxytocin for preventing postpartum haemorrhage on post-delivery haemoglobin.**
    2021.
    <https://pubmed.ncbi.nlm.nih.gov/34763599/>
29. Blum J, Winikoff B, Raghavan S, et al. **Treatment of post-partum haemorrhage
    with sublingual misoprostol versus oxytocin in women receiving prophylactic
    oxytocin: a double-blind, randomised, non-inferiority trial.** Lancet
    2010;375:217-23. Constrains misoprostol's weaker, slower `EMAX_MSO`/`KA_MSO`.
    <https://pubmed.ncbi.nlm.nih.gov/20060162/>
30. Winikoff B, Dabash R, Durocher J, et al. **Treatment of post-partum haemorrhage
    with sublingual misoprostol versus oxytocin in women not exposed to oxytocin
    during labour.** Lancet 2010;375:210-6. The companion trial in agonist-naive
    women — the closest published analogue of the model's `OTR0` contrast.
    <https://pubmed.ncbi.nlm.nih.gov/20060161/>
31. Gülmezoglu AM, Villar J, Ngoc NTN, et al. **WHO multicentre randomised trial of
    misoprostol in the management of the third stage of labour.** Lancet
    2001;358:689-95.
    <https://pubmed.ncbi.nlm.nih.gov/11551574/>
32. Lumbiganon P, Hofmeyr J, Gülmezoglu AM, Pinol A, Villar J. **Misoprostol
    dose-related shivering and pyrexia in the third stage of labour.** Br J Obstet
    Gynaecol 1999;106:304-8. The toxicity that caps misoprostol dosing.
    <https://pubmed.ncbi.nlm.nih.gov/10426235/>
33. Cole NM, Carvalho JCA, Downey K, et al. **Second-line uterotonics for uterine
    atony: a randomized controlled trial.** Obstet Gynecol 2024;144:566-74.
    Directly relevant to the model's central pharmacological claim — that after
    oxytocin failure the choice of receptor class, not the dose, is what matters.
    <https://pubmed.ncbi.nlm.nih.gov/39326051/>
34. Mousa HA, Blum J, Abou El Senoun G, Shakur H, Alfirevic Z. **Treatment for
    primary postpartum haemorrhage.** Cochrane Database Syst Rev 2014;2:CD003249.
    <https://pubmed.ncbi.nlm.nih.gov/24523225/>
35. Begley CM, Gyte GML, Devane D, McGuire W, Weeks A, Biesty LM. **Active versus
    expectant management for women in the third stage of labour.** Cochrane
    Database Syst Rev 2019;2:CD007412. The AMTSL effect size the model's
    prophylaxis block is compared against.
    <https://pubmed.ncbi.nlm.nih.gov/30754073/>
36. Yaju Y, Kataoka Y, Eto H, Horiuchi S, Mori R. **Prophylactic interventions after
    delivery of placenta for reducing bleeding during the postnatal period.**
    Cochrane Database Syst Rev 2013;11:CD009328.
    <https://pubmed.ncbi.nlm.nih.gov/24277681/>

## 7. Fibrinogen, dilution, and which factor fails first

37. Charbit B, Mandelbrot L, Samain E, et al. **The decrease of fibrinogen is an
    early predictor of the severity of postpartum hemorrhage.** J Thromb Haemost
    2007;5:266-73. The 2 g/L threshold that the model's `KFIB` Hill function is
    built to reproduce.
    <https://pubmed.ncbi.nlm.nih.gov/17087729/>
38. Cortet M, Deneux-Tharaux C, Dupont C, et al. **Association between fibrinogen
    level and severity of postpartum haemorrhage: secondary analysis of a
    prospective trial.** Br J Anaesth 2012;108:984-9.
    <https://pubmed.ncbi.nlm.nih.gov/22490316/>
39. Hiippala ST, Myllylä GJ, Vahtera EM. **Hemostatic factors and replacement of
    major blood loss with plasma-poor red cell concentrates.** Anesth Analg
    1995;81:360-5. The dilutional ordering the model derives from the ratio of
    critical to starting concentration — fibrinogen first.
    <https://pubmed.ncbi.nlm.nih.gov/7542432/>
40. de Lloyd L, Bovington R, Kaye A, et al. **Standard haemostatic tests following
    major obstetric haemorrhage.** Int J Obstet Anesth 2011;20:135-41.
    <https://pubmed.ncbi.nlm.nih.gov/21439811/>
41. Huissoud C, Carrabin N, Audibert F, et al. **Bedside assessment of fibrinogen
    level in postpartum haemorrhage by thrombelastometry.** BJOG 2009;116:1097-102.
    The FIBTEM read-out the model's `Fib` capture maps onto.
    <https://pubmed.ncbi.nlm.nih.gov/19459866/>
42. Gillissen A, van den Akker T, Caram-Deelder C, et al. **Coagulation parameters
    during the course of severe postpartum hemorrhage: a nationwide retrospective
    cohort study.** Blood Adv 2018;2:2433-42.
    <https://pubmed.ncbi.nlm.nih.gov/30266818/>
43. Collins PW, Cannings-John R, Bruynseels D, et al. **Viscoelastometric-guided
    early fibrinogen concentrate replacement during postpartum haemorrhage: OBS2, a
    double-blind randomized controlled trial.** Br J Anaesth 2017;119:411-21.
    The negative trial the model reproduces (empiric 4 g changes loss by 16 mL when
    fibrinogen is already above the cliff).
    <https://pubmed.ncbi.nlm.nih.gov/28969312/>
44. Wikkelsø AJ, Edwards HM, Afshari A, et al. **Pre-emptive treatment with
    fibrinogen concentrate for postpartum haemorrhage: randomized controlled trial
    (FIB-PPH).** Br J Anaesth 2015;114:623-33. The same null result independently.
    <https://pubmed.ncbi.nlm.nih.gov/25586727/>
45. Collins P, Abdul-Kadir R, Thachil J. **Management of postpartum haemorrhage:
    from research into practice, a narrative review of the literature and the
    Cardiff experience.** Int J Obstet Anesth 2019;37:106-17.
    <https://pubmed.ncbi.nlm.nih.gov/30322667/>

## 8. Fibrinolysis and tranexamic acid

46. WOMAN Trial Collaborators. **Effect of early tranexamic acid administration on
    mortality, hysterectomy, and other morbidities in women with post-partum
    haemorrhage (WOMAN): an international, randomised, double-blind,
    placebo-controlled trial.** Lancet 2017;389:2105-16. Death due to bleeding 1.5%
    vs 1.9% (RR 0.81), benefit confined to treatment within 3 h.
    <https://pubmed.ncbi.nlm.nih.gov/28456509/>
47. Gayet-Ageron A, Prieto-Merino D, Ker K, et al. **Effect of treatment delay on
    the effectiveness and safety of antifibrinolytics in acute severe haemorrhage:
    a meta-analysis of individual patient-level data from 40 138 bleeding
    patients.** Lancet 2018;391:125-32. The ~10%-per-15-minutes decay in benefit
    that the model re-derives from plasmin kinetics rather than assuming.
    <https://pubmed.ncbi.nlm.nih.gov/29126600/>
48. Shakur H, Beaumont D, Pavord S, Gayet-Ageron A, Ker K, Mousa HA.
    **Antifibrinolytic drugs for treating primary postpartum haemorrhage.**
    Cochrane Database Syst Rev 2018;2:CD012964.
    <https://pubmed.ncbi.nlm.nih.gov/29462500/>
49. Ducloy-Bouthors AS, Jude B, Duhamel A, et al. **High-dose tranexamic acid
    reduces blood loss in postpartum haemorrhage.** Crit Care 2011;15:R117.
    <https://pubmed.ncbi.nlm.nih.gov/21496253/>
50. Grassin-Delyle S, Semeraro M, Foissac F, et al. **Pharmacokinetics of tranexamic
    acid after intravenous, intramuscular, and oral routes: a prospective,
    randomised, crossover trial.** Br J Anaesth 2022;128:465-72. Constrains
    `V1_TX`, `CL` and the ~10 mg/L antifibrinolytic threshold.
    <https://pubmed.ncbi.nlm.nih.gov/34998508/>
51. Shakur-Still H, Roberts I, Grassin-Delyle S, et al. **Alternative routes for
    tranexamic acid treatment in obstetric bleeding (WOMAN-PharmacoTXA trial).**
    2023.
    <https://pubmed.ncbi.nlm.nih.gov/37019443/>
52. Dunn A, et al. **Evaluating tranexamic acid dosing strategies for postpartum
    hemorrhage: a population pharmacokinetic approach in pregnant individuals.**
    2024.
    <https://pubmed.ncbi.nlm.nih.gov/40384366/>
53. Iwamoto M. **Plasminogen-plasmin system IX. Specific binding of tranexamic acid
    to plasmin.** Thromb Diath Haemorrh 1975;33:573-85. The lysine-binding-site
    mechanism the model encodes as two separate effects (activation inhibition and
    clot protection).
    <https://pubmed.ncbi.nlm.nih.gov/125463/>
54. Suenson E, Thorsen S. **Initial plasmin-degradation of fibrin as the basis of a
    positive feed-back mechanism in fibrinolysis.** Eur J Biochem 1984;140:513-22.
    <https://pubmed.ncbi.nlm.nih.gov/6233145/>
55. Ker K, Roberts I, Shakur H, Coats TJ. **Antifibrinolytic drugs for acute
    traumatic injury.** Cochrane Database Syst Rev 2015;5:CD004896.
    <https://pubmed.ncbi.nlm.nih.gov/25956410/>

## 9. Transfusion, and what the products themselves cost

56. Pacheco LD, Saade GR, Costantine MM, Clark SL, Hankins GDV. **An update on the
    use of massive transfusion protocols in obstetrics.** Am J Obstet Gynecol
    2016;214:340-4. The loss-triggered (not haemoglobin-triggered) MTP the model
    implements.
    <https://pubmed.ncbi.nlm.nih.gov/26348379/>
57. Butwick AJ, Goodnough LT. **Transfusion and coagulation management in major
    obstetric hemorrhage.** Curr Opin Anaesthesiol 2015;28:275-84.
    <https://pubmed.ncbi.nlm.nih.gov/25812005/>
58. Butwick A. **How do I manage severe postpartum hemorrhage?** Transfusion
    2020;60:897-9.
    <https://pubmed.ncbi.nlm.nih.gov/32319687/>
59. Pasquier P, Gayat E, Rackelboom T, et al. **An observational study of the fresh
    frozen plasma : red blood cell ratio in postpartum hemorrhage.** Anesth Analg
    2013;116:155-61.
    <https://pubmed.ncbi.nlm.nih.gov/23223094/>
60. Williams CR, et al. **Transfusion of blood and blood products for the management
    of postpartum haemorrhage.** Cochrane Database Syst Rev 2025.
    <https://pubmed.ncbi.nlm.nih.gov/39911088/>
61. Giancarelli A, Birrer KL, Alban RF, Hobbs BP, Liu-DeRyke X. **Hypocalcemia in
    trauma patients receiving massive transfusion.** J Surg Res 2016;202:182-7.
    The prevalence and depth of ionised hypocalcaemia the citrate/calcium block is
    calibrated to.
    <https://pubmed.ncbi.nlm.nih.gov/27083965/>
62. Lier H, Krep H, Schroeder S, Stuber F. **Preconditions of hemostasis in trauma:
    a review. The influence of acidosis, hypocalcemia, anemia, and hypothermia on
    functional hemostasis.** J Trauma 2008;65:951-60. The source of the model's
    multiplicative acidosis / hypocalcaemia / hypothermia modifiers.
    <https://pubmed.ncbi.nlm.nih.gov/18849817/>
63. Meng ZH, Wolberg AS, Monroe DM 3rd, Hoffman M. **The effect of temperature and
    pH on the activity of factor VIIa: implications for the efficacy of high-dose
    factor VIIa in hypothermic and acidotic patients.** J Trauma 2003;55:886-91.
    Constrains `M_PHC` and `M_TMPC`.
    <https://pubmed.ncbi.nlm.nih.gov/14608161/>
64. Wolberg AS, Meng ZH, Monroe DM 3rd, Hoffman M. **A systematic evaluation of the
    effect of temperature on coagulation enzyme activity and platelet function.**
    J Trauma 2004;56:1221-8. The ~10%-per-degree rule behind the fluid-warmer
    result.
    <https://pubmed.ncbi.nlm.nih.gov/15211129/>
65. Rossaint R, Afshari A, Bouillon B, et al. **The European guideline on management
    of major bleeding and coagulopathy following trauma: sixth edition.** Crit Care
    2023;27:80. Source of the goal-directed targets used in the resuscitation arms.
    <https://pubmed.ncbi.nlm.nih.gov/36859355/>
66. Rajesh A, et al. **Aggressive calcium chloride dosing reduces early mortality in
    trauma patients receiving whole blood resuscitation.** J Trauma Acute Care Surg
    2026.
    <https://pubmed.ncbi.nlm.nih.gov/41995161/>
67. Shandaliy Y, et al. **Impact of a calcium replacement protocol during massive
    transfusion in trauma patients at a level 2 trauma center.** Am J Health Syst
    Pharm 2024.
    <https://pubmed.ncbi.nlm.nih.gov/38578328/>
68. Carson JL, et al. **Transfusion thresholds and other strategies for guiding red
    blood cell transfusion.** Cochrane Database Syst Rev 2025.
    <https://pubmed.ncbi.nlm.nih.gov/41114449/>

## 10. Fluid strategy: the arm the model finds most costly

69. Gillissen A, Henriquez DDCA, van den Akker T, et al. **Association between fluid
    management and dilutional coagulopathy in severe postpartum haemorrhage: a
    nationwide retrospective cohort study.** BMC Pregnancy Childbirth 2018;18:398.
    The observational counterpart of the model's crystalloid-first arm.
    <https://pubmed.ncbi.nlm.nih.gov/30305108/>
70. Henriquez DDCA, Bloemenkamp KWM, Loeff RM, et al. **Fluid resuscitation during
    persistent postpartum haemorrhage and maternal outcome: a nationwide cohort
    study.** Eur J Obstet Gynecol Reprod Biol 2019;235:49-56.
    <https://pubmed.ncbi.nlm.nih.gov/30784827/>

## 11. Mechanical and surgical control

71. Suarez S, Conde-Agudelo A, Borovac-Pinheiro A, et al. **Uterine balloon tamponade
    for the treatment of postpartum hemorrhage: a systematic review and
    meta-analysis.** Am J Obstet Gynecol 2020;222:293.e1-52. The ~85% success rate
    the balloon arm is checked against.
    <https://pubmed.ncbi.nlm.nih.gov/31917139/>
72. Doumouchtsis SK, Papageorghiou AT, Arulkumaran S. **Systematic review of
    conservative management of postpartum hemorrhage: what to do when medical
    treatment fails.** Obstet Gynecol Surv 2007;62:540-7. The escalation sequence
    the model's mechanical block encodes.
    <https://pubmed.ncbi.nlm.nih.gov/17634155/>
73. Kellie FJ, Wandabwa JN, Mousa HA, Weeks AD. **Mechanical and surgical
    interventions for treating primary postpartum haemorrhage.** Cochrane Database
    Syst Rev 2020;7:CD013663.
    <https://pubmed.ncbi.nlm.nih.gov/32609374/>
74. El-Sokkary M, Wahba K, El-Shahawy Y. **Uterine salvage management for atonic
    postpartum hemorrhage using "modified Lynch suture".** BMC Pregnancy Childbirth
    2016;16:251.
    <https://pubmed.ncbi.nlm.nih.gov/27567670/>
75. Aoki M, et al. **Effects of descending aortic occlusion for massive obstetric
    hemorrhage: nationwide analysis of maternal death in Japan.** Int J Gynaecol
    Obstet 2024. Evidence for the aortic-occlusion bridge that the model finds
    worth more than a later definitive procedure.
    <https://pubmed.ncbi.nlm.nih.gov/42101038/>
76. van de Voort JC, et al. **Resuscitative endovascular balloon occlusion of the
    aorta (REBOA) for non-trauma patients in an urban hospital.** 2024.
    <https://pubmed.ncbi.nlm.nih.gov/39351589/>
77. Zheng TQ, et al. **Innovative approach for rapid control of postpartum
    hemorrhage with abdominal pressure.** Matern Fetal Med 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42051647/>
78. Koyama E, Naruse K, Shigetomi H, et al. **Combination of B-Lynch brace suture
    and uterine artery embolization for atonic bleeding after cesarean section.**
    J Obstet Gynaecol Res 2012;38:345-8.
    <https://pubmed.ncbi.nlm.nih.gov/22136878/>

## 12. Guidelines and quality-improvement programmes

79. Committee on Practice Bulletins-Obstetrics. **Practice Bulletin No. 183:
    Postpartum Hemorrhage.** Obstet Gynecol 2017;130:e168-86.
    <https://pubmed.ncbi.nlm.nih.gov/28937571/>
80. **Prevention and Management of Postpartum Haemorrhage: Green-top Guideline
    No. 52.** BJOG 2017;124:e106-49. (RCOG)
    <https://pubmed.ncbi.nlm.nih.gov/27981719/>
81. Escobar MF, Nassar AH, Theron G, et al. **FIGO recommendations on the management
    of postpartum hemorrhage 2022.** Int J Gynaecol Obstet 2022;157(Suppl 1):3-50.
    <https://pubmed.ncbi.nlm.nih.gov/35297039/>
82. Collis RE, Collins PW. **Haemostatic management of obstetric haemorrhage.**
    Anaesthesia 2015;70(Suppl 1):78-86.
    <https://pubmed.ncbi.nlm.nih.gov/25440400/>
83. Bell SF, Watkins A, John M, et al. **Designing and implementing an all Wales
    postpartum haemorrhage quality improvement project: OBS Cymru.** BMJ Open Qual
    2020;9:e000854. The clearest demonstration that the model's central claim —
    time-to-intervention dominates choice-of-intervention — holds in practice.
    <https://pubmed.ncbi.nlm.nih.gov/32273281/>
84. Dale M, et al. **What is the economic cost of providing an all Wales postpartum
    haemorrhage quality improvement initiative (OBS Cymru)?** 2022.
    <https://pubmed.ncbi.nlm.nih.gov/36066836/>
85. **WHO recommendations for the prevention and treatment of postpartum
    haemorrhage.** World Health Organization.
    <https://www.who.int/publications/i/item/9789241548502>

## 13. Consequences for the survivor

86. Karpati PCJ, Rossignol M, Pirot M, et al. **High incidence of myocardial
    ischemia during postpartum hemorrhage.** Anesthesiology 2004;100:30-6.
    <https://pubmed.ncbi.nlm.nih.gov/14695721/>
87. Mehrabadi A, Liu S, Bartholomew S, et al. **Investigation of a rise in obstetric
    acute renal failure in the United States, 1999-2011.** Obstet Gynecol
    2016;127:899-906.
    <https://pubmed.ncbi.nlm.nih.gov/27054929/>
88. Schury MP, Adigun R. **Sheehan Syndrome.** StatPearls. The pituitary
    consequence of the model's "minutes below MAP 50" integral.
    <https://pubmed.ncbi.nlm.nih.gov/29083621/>
89. Jeon C, et al. **A case report of acute Sheehan syndrome with a review of 29
    existing reports from 1990 to 2024.** 2025.
    <https://pubmed.ncbi.nlm.nih.gov/41089439/>
90. Froeliger A, Deneux-Tharaux C, Loussert L, et al. **Trajectories of
    childbirth-related posttraumatic stress symptoms after a vaginal delivery.**
    Am J Obstet Gynecol 2026. A reminder that blood loss is not the only endpoint.
    <https://pubmed.ncbi.nlm.nih.gov/42361947/>

---

### Notes on what the literature does NOT settle

Three of the model's structural commitments are, as far as this reference list
goes, **assumptions rather than measurements**, and the sensitivity notes at the
end of `pph_mrgsolve_model.R` say how to attack each one:

- the exponent linking myometrial tone to placental-bed patency (`NP = 4`, from
  Poiseuille) — no direct human measurement of flow-versus-tone at the placental
  bed appears to exist;
- the ceiling on what a clot can seal without mechanical apposition
  (`CSEAL0 = 0.35`) — inferred from the clinical existence of uterotonic-
  refractory atony, not measured;
- the time constant of oxytocin-receptor recovery in vivo (`KREC`, t½ 2.6 h) —
  extrapolated from binding-site and in-vitro contractility data (refs 20-22),
  which do not measure recovery directly in the postpartum uterus.

Everything the model claims about therapy follows from those three numbers plus
the shear gate, so they are where a falsification effort should start.
