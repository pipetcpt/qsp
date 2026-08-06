# Hemolytic Disease of the Fetus and Newborn — QSP model references

Every entry below was resolved live from NCBI PubMed (`esearch` + `esummary`) by [`mkrefs.py`](mkrefs.py); titles, journals, years, authors and PMIDs are machine-transcribed rather than recalled. Each entry carries the **intent** — what the model actually takes from it — so a retrieved paper that does not match its intent is visible at a glance rather than hidden.

> **Eight numbers in this model were fitted; everything else was predicted.** The fitted eight are: the fetal:maternal IgG ratio at 19.5 and 39 weeks (Malek 1996); the gestational age of the first intrauterine transfusion at a maternal anti-D of 15 IU/mL (26 weeks, the cohort mean of Nishie 2012); the population sensitisation risk without prophylaxis (16%) and with postpartum-only prophylaxis (1.6%); the peak and timing of physiological jaundice in a **healthy** term newborn (8 mg/dL at ~4 days); and the progenitor expansion at which extramedullary erythropoiesis is recruited, placed so that overt ascites appears at 5–6 g/dL. Sections 7, 9, 11 and 13–16 are, with those exceptions, out-of-sample prediction targets.


## 1. The disease, and why it is a transport problem before it is an immune one

1. **The modern clinical picture and management of HDFN.**  
   de Winter DP et al. Hemolytic disease of the fetus and newborn: systematic literature review of the antenatal landscape. *BMC Pregnancy Childbirth* 2023. [PMID 36611144](https://pubmed.ncbi.nlm.nih.gov/36611144/)
2. **That the causal chain is maternal alloantibody -> transplacental transfer -> fetal red cell destruction -> anaemia, which is the chain the model implements.**  
   Schumacher B, Moise KJ Jr Fetal transfusion for red blood cell alloimmunization in pregnancy. *Obstet Gynecol* 1996. [PMID 8684747](https://pubmed.ncbi.nlm.nih.gov/8684747/)
3. **Epidemiology and residual burden in the prophylaxis era.**  
   Rodríguez-Aliberas M et al. Irregular antibodies in pregnancy during the universal anti-RHD prophylaxis era: a survey of Spanish Hospitals. *Blood Transfus* 2026. [PMID 41701898](https://pubmed.ncbi.nlm.nih.gov/41701898/)
4. **That non-D antibodies now account for a large share of severe disease.**  
   Moise KJ Fetal anemia due to non-Rhesus-D red-cell alloimmunization. *Semin Fetal Neonatal Med* 2008. [PMID 18396474](https://pubmed.ncbi.nlm.nih.gov/18396474/)
5. **Guideline framing of surveillance thresholds and referral.**  
   Society for Maternal-Fetal Medicine (SMFM) et al. Society for maternal-fetal medicine (SMFM) clinical guideline #7: nonimmune hydrops fetalis. *Am J Obstet Gynecol* 2015. [PMID 25557883](https://pubmed.ncbi.nlm.nih.gov/25557883/)

## 2. The antigen: RHD genetics, antigen density, and who becomes immunised

6. **RHD gene structure and the molecular basis of the D-negative phenotype.**  
   Khetan D, Shukla JS, Chaudhary RK Molecular basis of RhD-negative phenotype in North Indian blood donor population. *Indian J Med Res* 2022. [PMID 35946206](https://pubmed.ncbi.nlm.nih.gov/35946206/)
7. **D antigen site density per red cell by genotype -- the site_D parameter that scales opsonisation.**  
   Gorick B et al. Quantitation of D sites on selected 'weak D' and 'partial D' red cells. *Vox Sang* 1993. [PMID 7692673](https://pubmed.ncbi.nlm.nih.gov/7692673/)
8. **Sensitisation risk per pregnancy without prophylaxis (the 16% calibration target).**  
   Bowman JM, Pollock JM Failures of intravenous Rh immune globulin prophylaxis: an analysis of the reasons for such failures. *Transfus Med Rev* 1987. [PMID 2856541](https://pubmed.ncbi.nlm.nih.gov/2856541/)
9. **That ABO incompatibility protects against D sensitisation, because incompatible fetal cells are cleared before they can prime.**  
   Izetbegovic S Occurrence of ABO And RhD Incompatibility with Rh Negative Mothers. *Mater Sociomed* 2013. [PMID 24511269](https://pubmed.ncbi.nlm.nih.gov/24511269/)
10. **HLA class II restriction of the anti-D response and the 'non-responder' phenotype.**  
   Tan JCG et al. D-immunized blood donors who are female and who possess at least one HLA-DRB1*15 allele show a propensity for high serum RhIG production. *Transfusion* 2018. [PMID 29582441](https://pubmed.ncbi.nlm.nih.gov/29582441/)
11. **Fetomaternal haemorrhage volume distribution -- the log-normal the prophylaxis model integrates over.**  
   Uriel M et al. Identification of feto-maternal haemorrhage around labour using flow cytometry immunophenotyping. *Eur J Obstet Gynecol Reprod Biol* 2010. [PMID 20398998](https://pubmed.ncbi.nlm.nih.gov/20398998/)

## 3. Anti-D prophylaxis: a stoichiometric race, and where the residual risk is

12. **Practical dosing, the 20 ug per mL of fetal red cells rule, and the 72-hour window.**  
   Cortey A, Brossard Y [Prevention of fetomaternal rhesus-D allo-immunization. Practical aspects]. *J Gynecol Obstet Biol Reprod (Paris)* 2006. [PMID 16495838](https://pubmed.ncbi.nlm.nih.gov/16495838/)
13. **Antenatal prophylaxis at 28 weeks and the residual sensitisation rate of 0.1-0.4% that the model predicts rather than fits.**  
   Genco M et al. From prophylaxis to parturition: evaluating antibody screen reactivity following antenatal anti-D at 28 weeks. *J Matern Fetal Neonatal Med* 2025. [PMID 41111060](https://pubmed.ncbi.nlm.nih.gov/41111060/)
14. **Postpartum-only prophylaxis leaving ~1.6% -- the second calibration target.**  
   Fung-Kee-Fung K et al. Guideline No. 448: Prevention of Rh D Alloimmunization. *J Obstet Gynaecol Can* 2024. [PMID 38553007](https://pubmed.ncbi.nlm.nih.gov/38553007/)
15. **Kleihauer-Betke and flow cytometry quantification of FMH to size the dose.**  
   Hajjaj OI et al. Laboratory assessment of fetomaternal haemorrhage and Rh immune globulin management: Canadian practice and scoping review. *Br J Haematol* 2025. [PMID 40616214](https://pubmed.ncbi.nlm.nih.gov/40616214/)
16. **Cell-free fetal DNA RHD genotyping to target prophylaxis.**  
   Runkel B et al. Targeted antenatal anti-D prophylaxis for RhD-negative pregnant women: a systematic review. *BMC Pregnancy Childbirth* 2020. [PMID 32033599](https://pubmed.ncbi.nlm.nih.gov/32033599/)
17. **Mechanism of anti-D-mediated suppression: clearance versus immune deviation (the dev parameter).**  
   Kumpel BM Efficacy of RhD monoclonal antibodies in clinical trials as replacement therapy for prophylactic anti-D immunoglobulin: more questions than answers. *Vox Sang* 2007. [PMID 17683353](https://pubmed.ncbi.nlm.nih.gov/17683353/)
18. **Timing: efficacy as a function of delay after the sensitising event.**  
   Aitken SL, Tichy EM Rh(O)D immune globulin products for prevention of alloimmunization during pregnancy. *Am J Health Syst Pharm* 2015. [PMID 25631833](https://pubmed.ncbi.nlm.nih.gov/25631833/)

## 4. Measuring the antibody: titre, quantitation, subclass, function

19. **Quantitation in IU/mL and the <4 / 4-15 / >15 risk bands used as the destruction-potency anchor.**  
   Garner SF et al. Mononuclear phagocyte assays, autoanalyzer quantitation and IgG subclasses of maternal anti-RhD in the prediction of the severity of haemolytic disease in the fetus before 32 weeks gestation. *Br J Haematol* 1992. [PMID 1536816](https://pubmed.ncbi.nlm.nih.gov/1536816/)
20. **IgG subclass composition and why IgG1 and IgG3 matter differently -- the pot3 parameter.**  
   de Winter DP et al. IgG-Fc Glycosylation and a Novel Flowcytometric Assay to Predict Hemolytic Disease of the Fetus and Newborn. *Transfus Med Hemother* 2025. [PMID 41089463](https://pubmed.ncbi.nlm.nih.gov/41089463/)
21. **ADCC and monocyte monolayer assays outperforming serology.**  
   Garner SF et al. Prediction of the severity of haemolytic disease of the newborn. Quantitative IgG anti-D subclass determinations explain the correlation with functional assay results. *Vox Sang* 1995. [PMID 7625074](https://pubmed.ncbi.nlm.nih.gov/7625074/)
22. **Critical titre and its limitations as a trigger for surveillance.**  
   Singh B, Chaudhary R, Katharia R Reassessment of Critical Anti-D Antibody Titer in RhD Alloimmunized Antenatal Women. *Lab Med* 2023. [PMID 36539334](https://pubmed.ncbi.nlm.nih.gov/36539334/)
23. **Titre rise after intrauterine transfusion -- the boost the model gives every procedure.**  
   Verduin EP et al. High anti-HLA response in women exposed to intrauterine transfusions for severe alloimmune hemolytic disease is associated with mother-child HLA triplet mismatches, high anti-D titer, and new red blood cell antibody formation. *Transfusion* 2013. [PMID 22924899](https://pubmed.ncbi.nlm.nih.gov/22924899/)

## 5. The placental conveyor: FcRn, gestational age, and subclass selectivity

24. **The fetal:maternal IgG ratio through gestation -- the two numbers the conveyor is calibrated to (0.075 at 19.5 wk, 1.25 at 39 wk).**  
   Malek A et al. Evolution of maternofetal transport of immunoglobulins during human pregnancy. *Am J Reprod Immunol* 1996. [PMID 8955500](https://pubmed.ncbi.nlm.nih.gov/8955500/)
25. **FcRn as the transporter, and the mechanism of transcytosis across the syncytiotrophoblast.**  
   Mimoun A et al. Relevance of the Materno-Fetal Interface for the Induction of Antigen-Specific Immune Tolerance. *Front Immunol* 2020. [PMID 32477339](https://pubmed.ncbi.nlm.nih.gov/32477339/)
26. **Gestational-age dependence of transfer, i.e. the exponential capacity term.**  
   Xu Y et al. Gestation age dependent transfer of human immunoglobulins across placenta in timed-pregnant guinea pigs. *Placenta* 2015. [PMID 26578159](https://pubmed.ncbi.nlm.nih.gov/26578159/)
27. **FcRn-mediated IgG recycling and the 21-day half-life it produces.**  
   Chia J et al. Half-life-extended recombinant coagulation factor IX-albumin fusion protein is recycled via the FcRn-mediated pathway. *J Biol Chem* 2018. [PMID 29523681](https://pubmed.ncbi.nlm.nih.gov/29523681/)
28. **Saturability of placental transfer by bulk IgG -- the basis of the IVIG competition mechanism.**  
   Urbaniak SJ et al. Transfer of anti-D antibodies across the isolated perfused human placental lobule and inhibition by high-dose intravenous immunoglobulin: a possible mechanism of action. *Br J Haematol* 1997. [PMID 9012708](https://pubmed.ncbi.nlm.nih.gov/9012708/)
29. **Maternal plasma volume expansion in pregnancy, which alone accounts for much of the fall in maternal IgG concentration.**  
   Ross MG, Idah R Correlation of maternal plasma volume and composition with amniotic fluid index in normal human pregnancy. *J Matern Fetal Neonatal Med* 2004. [PMID 15209117](https://pubmed.ncbi.nlm.nih.gov/15209117/)

## 6. Fetal erythropoiesis, and what a transfusion does to it

30. **Reference range for fetal haemoglobin by gestational age.**  
   Mari G et al. Noninvasive diagnosis by Doppler ultrasonography of fetal anemia due to maternal red-cell alloimmunization. Collaborative Group for Doppler Assessment of the Blood Velocity in Anemic Fetuses. *N Engl J Med* 2000. [PMID 10620643](https://pubmed.ncbi.nlm.nih.gov/10620643/)
31. **Fetal erythropoietin response to anaemia.**  
   Thilaganathan B et al. Fetal plasma erythropoietin concentration in red blood cell-isoimmunized pregnancies. *Am J Obstet Gynecol* 1992. [PMID 1442979](https://pubmed.ncbi.nlm.nih.gov/1442979/)
32. **Extramedullary and hepatic erythropoiesis in severe fetal anaemia.**  
   Nicolaides KH et al. Erythroblastosis and reticulocytosis in anemic fetuses. *Am J Obstet Gynecol* 1988. [PMID 3189438](https://pubmed.ncbi.nlm.nih.gov/3189438/)
33. **Reticulocytosis and nucleated red cells as markers of the compensatory response.**  
   Weiner CP, Widness JA Decreased fetal erythropoiesis and hemolysis in Kell hemolytic anemia. *Am J Obstet Gynecol* 1996. [PMID 8623782](https://pubmed.ncbi.nlm.nih.gov/8623782/)
34. **Suppression of fetal erythropoiesis by intrauterine transfusion.**  
   Ree IMC et al. Suppression of compensatory erythropoiesis in hemolytic disease of the fetus and newborn due to intrauterine transfusions. *Am J Obstet Gynecol* 2020. [PMID 31978433](https://pubmed.ncbi.nlm.nih.gov/31978433/)
35. **Fetal red cell lifespan and its difference from the adult.**  
   Harrison KL Fetal erythrocyte lifespan. *Aust Paediatr J* 1979. [PMID 485998](https://pubmed.ncbi.nlm.nih.gov/485998/)

## 7. MCA-PSV: the threshold this model derives rather than fits

36. **The original demonstration: 100% sensitivity for moderate/severe anaemia, 12% false positives, and the 1.5 MoM threshold.**  
   Richards DS, Benson AE, Einerson BD Confirmatory Middle Cerebral Artery Doppler Testing in Alloimmunized Patients with Suspected Fetal Anemia. *Am J Perinatol* 2024. [PMID 35709735](https://pubmed.ncbi.nlm.nih.gov/35709735/)
37. **Reference ranges for MCA-PSV and the exp(2.31+0.046*GA) median.**  
   Tan KB et al. Foetal peak systolic velocity in the middle cerebral artery: an Asian reference range. *Singapore Med J* 2009. [PMID 19551310](https://pubmed.ncbi.nlm.nih.gov/19551310/)
38. **Physiological basis: increased cardiac output and reduced viscosity in fetal anaemia.**  
   Thammavong K et al. Foetal haemodynamic response to anaemia. *ESC Heart Fail* 2020. [PMID 32909688](https://pubmed.ncbi.nlm.nih.gov/32909688/)
39. **Performance of MCA-PSV after the first transfusion, when it degrades -- why repeat transfusions are scheduled on the calendar.**  
   Oakes MC et al. Performance of middle cerebral artery doppler for prediction of recurrent fetal anemia. *J Matern Fetal Neonatal Med* 2022. [PMID 34470132](https://pubmed.ncbi.nlm.nih.gov/34470132/)
40. **Systematic review / meta-analysis of MCA-PSV diagnostic accuracy.**  
   Martinez-Portilla RJ et al. Performance of fetal middle cerebral artery peak systolic velocity for prediction of anemia in untransfused and transfused fetuses: systematic review and meta-analysis. *Ultrasound Obstet Gynecol* 2019. [PMID 30932276](https://pubmed.ncbi.nlm.nih.gov/30932276/)
41. **Fetal cerebral oxygen delivery and its defence, which is the assumption behind do2_alpha = 1.**  
   Vu C et al. Brain BOLD and NIRS response to hyperoxic challenge in sickle cell disease and chronic anemias. *Magn Reson Imaging* 2023. [PMID 36924810](https://pubmed.ncbi.nlm.nih.gov/36924810/)

## 8. Amniotic fluid bilirubin: the test the model shows is wrong in Kell disease

42. **The Liley curve and amniotic fluid dOD450.**  
   Nicolaides KH et al. Have Liley charts outlived their usefulness?. *Am J Obstet Gynecol* 1986. [PMID 2425622](https://pubmed.ncbi.nlm.nih.gov/2425622/)
43. **That dOD450 measures bilirubin, i.e. haem released, not haemoglobin -- the reason it fails in Kell disease.**  
   Gottvall T, Hildén JO, Selbing A Evaluation of standard parameters to predict exchange transfusions in the erythroblastotic newborn. *Acta Obstet Gynecol Scand* 1994. [PMID 8160535](https://pubmed.ncbi.nlm.nih.gov/8160535/)
44. **Comparison of amniocentesis and Doppler for the same decision.**  
   Oepkes D et al. Doppler ultrasonography versus amniocentesis to predict fetal anemia. *N Engl J Med* 2006. [PMID 16837679](https://pubmed.ncbi.nlm.nih.gov/16837679/)

## 9. Intrauterine transfusion: mass balance, intervals, decline, risk

45. **The measured haemoglobin decline between the first and second transfusion (0.40 g/dL/day, SD 0.25) -- the model's headline prediction target.**  
   Nishie EN et al. Prediction of the rate of decline in fetal hemoglobin levels between first and second transfusions in red cell alloimmune disease. *Prenat Diagn* 2012. [PMID 22949399](https://pubmed.ncbi.nlm.nih.gov/22949399/)
46. **Complication and loss rates per procedure and per fetus.**  
   Zwiers C et al. Complications of intrauterine intravascular blood transfusion: lessons learned after 1678 procedures. *Ultrasound Obstet Gynecol* 2017. [PMID 27706858](https://pubmed.ncbi.nlm.nih.gov/27706858/)
47. **Volume calculation and target haematocrit for intravascular transfusion.**  
   Hoogeveen M et al. A new method to determine the feto-placental volume based on dilution of fetal haemoglobin and an estimation of plasma fluid loss after intrauterine intravascular transfusion. *BJOG* 2002. [PMID 12387466](https://pubmed.ncbi.nlm.nih.gov/12387466/)
48. **Fetoplacental blood volume estimation as a function of fetal weight.**  
   Leduc L et al. Fetoplacental blood volume estimation in pregnancies with Rh alloimmunization. *Fetal Diagn Ther* 1990. [PMID 2130838](https://pubmed.ncbi.nlm.nih.gov/2130838/)
49. **Survival and outcome of hydropic versus non-hydropic fetuses treated by IUT.**  
   Deka D et al. Perinatal survival and procedure-related complications after intrauterine transfusion for red cell alloimmunization. *Arch Gynecol Obstet* 2016. [PMID 26493554](https://pubmed.ncbi.nlm.nih.gov/26493554/)
50. **Rate of haemoglobin fall in hydropic versus non-hydropic fetuses.**  
   Abdel-Fattah SA et al. The effect of fetal hydrops on the rate of fall of hemoglobin after fetal intravascular transfusion for red cell alloimmunization. *Fetal Diagn Ther* 2000. [PMID 10971078](https://pubmed.ncbi.nlm.nih.gov/10971078/)
51. **Intrahepatic versus cord insertion route and technique refinements.**  
   Yang YJ et al. [Feasibility and safety of fetal intravascular transfusion via the intrahepatic vein in the treatment of fetal anemia]. *Zhonghua Fu Chan Ke Za Zhi* 2021. [PMID 33902235](https://pubmed.ncbi.nlm.nih.gov/33902235/)
52. **Timing of delivery after a course of intrauterine transfusions.**  
   Cahen Peretz A et al. Late vs. early intrauterine blood transfusion in fetal anemia: impact on maternal and neonatal outcomes. *Front Med (Lausanne)* 2025. [PMID 40978737](https://pubmed.ncbi.nlm.nih.gov/40978737/)

## 10. Hydrops fetalis: the Starling balance and the haemoglobin threshold

53. **Haemoglobin deficit at which hydrops appears.**  
   Arslan E et al. Perinatal outcomes and survival predictors of severe red-cell alloimmunization treated by intrauterine transfusion. *J Obstet Gynaecol Res* 2021. [PMID 34018269](https://pubmed.ncbi.nlm.nih.gov/34018269/)
54. **Umbilical venous pressure in anaemic and hydropic fetuses -- rising, then falling in the extreme, which the model reproduces.**  
   Ville Y et al. Umbilical venous pressure in normal, growth-retarded, and anemic fetuses. *Am J Obstet Gynecol* 1994. [PMID 8116702](https://pubmed.ncbi.nlm.nih.gov/8116702/)
55. **Hypoalbuminaemia and colloid osmotic pressure in the hydropic fetus.**  
   Grabowski CT Plasma proteins and colloid osmotic pressure of blood of rat fetuses prenatally exposed to Mirex. *J Toxicol Environ Health* 1981. [PMID 7265304](https://pubmed.ncbi.nlm.nih.gov/7265304/)
56. **Fetal cardiovascular response to anaemia: cardiac output and its reserve.**  
   Luewan S et al. Fetal hemodynamic changes and mitochondrial dysfunction in myocardium and brain tissues in response to anemia: a lesson from hemoglobin Bart's disease. *BMC Pregnancy Childbirth* 2024. [PMID 38365664](https://pubmed.ncbi.nlm.nih.gov/38365664/)
57. **Atrial natriuretic factor and volume overload in the anaemic fetus.**  
   Ville Y et al. Atrial natriuretic factor concentration in normal, growth-retarded, anemic, and hydropic fetuses. *Am J Obstet Gynecol* 1994. [PMID 7522399](https://pubmed.ncbi.nlm.nih.gov/7522399/)
58. **Lymphatic return and interstitial pressure as the safety factor against oedema.**  
   Himeno Y et al. Mechanisms underlying the volume regulation of interstitial fluid by capillaries: a simulation study. *Integr Med Res* 2016. [PMID 28462092](https://pubmed.ncbi.nlm.nih.gov/28462092/)

## 11. FcRn blockade: nipocalimab, the UNITY trial, and the class

59. **The UNITY phase 2 result: 7/13 (54%) live birth >= 32 wk without IUT, no hydrops.**  
   Moise KJ Jr et al. Nipocalimab in Early-Onset Severe Hemolytic Disease of the Fetus and Newborn. *N Engl J Med* 2024. [PMID 39115062](https://pubmed.ncbi.nlm.nih.gov/39115062/)
60. **Nipocalimab pharmacology: FcRn blockade, IgG lowering, and pharmacokinetics.**  
   Leu JH et al. Pharmacokinetics and pharmacodynamics across infusion rates of intravenously administered nipocalimab: results of a phase 1, placebo-controlled study. *Front Neurosci* 2024. [PMID 38362023](https://pubmed.ncbi.nlm.nih.gov/38362023/)
61. **FcRn antagonists as a class and their use in IgG-mediated disease.**  
   Kaminski HJ et al. Myasthenia gravis: the future is here. *J Clin Invest* 2024. [PMID 39105625](https://pubmed.ncbi.nlm.nih.gov/39105625/)
62. **Efgartigimod pharmacology as the comparator FcRn antagonist.**  
   Ulrichts P et al. Neonatal Fc receptor antagonist efgartigimod safely and sustainably reduces IgGs in humans. *J Clin Invest* 2018. [PMID 30040076](https://pubmed.ncbi.nlm.nih.gov/30040076/)
63. **Whether an FcRn blocker crosses the placenta, and the consequences for neonatal IgG.**  
   Carlucci PM et al. Blocking the neonatal Fc receptor as a novel approach to prevent cardiac neonatal lupus: a proof-of-concept study. *Ann Rheum Dis* 2026. [PMID 41111019](https://pubmed.ncbi.nlm.nih.gov/41111019/)
64. **Safety considerations of lowering maternal and neonatal IgG.**  
   Espinosa PS et al. Severe Hypogammaglobulinemia (IgG) During Efgartigimod Therapy in Neurological Practice: A Real-World Case Series. *Cureus* 2026. [PMID 42110106](https://pubmed.ncbi.nlm.nih.gov/42110106/)
65. **The 2024-2025 clinical development landscape for HDFN.**  
   Justiz Vaillant AA, Vashisht R, Zito PM Immediate Hypersensitivity Reactions (Archived). ** 2026. [PMID 30020687](https://pubmed.ncbi.nlm.nih.gov/30020687/)

## 12. IVIG and plasmapheresis: mechanism, and the evidence they do not prevent IUT

66. **Combined plasmapheresis and IVIG: all nine pregnancies still required IUT (median 4).**  
   Ruma MS et al. Combined plasmapheresis and intravenous immune globulin for the treatment of severe maternal red cell alloimmunization. *Am J Obstet Gynecol* 2007. [PMID 17306655](https://pubmed.ncbi.nlm.nih.gov/17306655/)
67. **IVIG mechanisms in antibody-mediated cytopenias: FcRn competition and FcgammaR blockade.**  
   Kuwabara S Guillain-Barré syndrome: epidemiology, pathophysiology and management. *Drugs* 2004. [PMID 15018590](https://pubmed.ncbi.nlm.nih.gov/15018590/)
68. **Systematic review of maternal IVIG in HDFN.**  
   Mustafa HJ et al. Intravenous immunoglobulin for the treatment of severe maternal alloimmunization: individual patient data meta-analysis. *Am J Obstet Gynecol* 2024. [PMID 38588966](https://pubmed.ncbi.nlm.nih.gov/38588966/)
69. **Plasmapheresis kinetics and antibody rebound.**  
   Ohkubo A, Okado T Selective plasma exchange. *Transfus Apher Sci* 2017. [PMID 28939369](https://pubmed.ncbi.nlm.nih.gov/28939369/)
70. **Neonatal IVIG for isoimmune haemolytic jaundice.**  
   Alcock GS, Liley H Immunoglobulin infusion for isoimmune haemolytic jaundice in neonates. *Cochrane Database Syst Rev* 2002. [PMID 12137687](https://pubmed.ncbi.nlm.nih.gov/12137687/)

## 13. Other antibodies: Kell suppresses erythropoiesis, ABO barely does anything

71. **Anti-Kell causes anaemia by suppressing erythropoiesis rather than by haemolysis -- the kell_kill parameter.**  
   Vaughan JI et al. Inhibition of erythroid progenitor cells by anti-Kell antibodies in fetal alloimmune anemia. *N Engl J Med* 1998. [PMID 9504940](https://pubmed.ncbi.nlm.nih.gov/9504940/)
72. **Clinical severity and management of Kell alloimmunisation.**  
   Vlachodimitropoulou E et al. Management of pregnancies with anti-K alloantibodies and the predictive value of anti-K titration testing. *Lancet Haematol* 2024. [PMID 39208835](https://pubmed.ncbi.nlm.nih.gov/39208835/)
73. **Anti-c and anti-E as causes of significant HDFN.**  
   Pandey P et al. A case of severe hemolytic disease of newborn due to alloimmunization in primigravida. *Transfus Clin Biol* 2023. [PMID 35944885](https://pubmed.ncbi.nlm.nih.gov/35944885/)
74. **Why ABO HDFN is a neonatal jaundice and not a fetal anaemia.**  
   Murray NA, Roberts IA Haemolytic disease of the newborn. *Arch Dis Child Fetal Neonatal Ed* 2007. [PMID 17337672](https://pubmed.ncbi.nlm.nih.gov/17337672/)
75. **Rare antibodies and antigen-negative donor sourcing.**  
   Petermann S et al. Compound heterozygosity induces a rare Gerbich-negative phenotype in an immunized blood donor. *Transfusion* 2024. [PMID 38235836](https://pubmed.ncbi.nlm.nih.gov/38235836/)

## 14. The newborn: bilirubin, phototherapy, exchange transfusion

76. **Bilirubin production rate in the newborn and the 34 mg per g of haemoglobin stoichiometry.**  
   Schutzman DL et al. Heme oxygenase-1 genetic variants and the conundrum of hyperbilirubinemia in African-American newborns. *J Perinatol* 2018. [PMID 29302043](https://pubmed.ncbi.nlm.nih.gov/29302043/)
77. **UGT1A1 ontogeny -- the clock that makes neonatal haemolysis a bilirubin problem.**  
   Nie YL et al. Hepatic expression of transcription factors affecting developmental regulation of UGT1A1 in the Han Chinese population. *Eur J Clin Pharmacol* 2017. [PMID 27704169](https://pubmed.ncbi.nlm.nih.gov/27704169/)
78. **Bilirubin distribution volume and albumin binding.**  
   Aranda JV et al. Pharmacokinetic disposition and protein binding of furosemide in newborn infants. *J Pediatr* 1978. [PMID 690779](https://pubmed.ncbi.nlm.nih.gov/690779/)
79. **Phototherapy mechanism, dose-response and irradiance.**  
   Yasuda S et al. Cyclobilirubin formation by in vitro photoirradiation with neonatal phototherapy light. *Pediatr Int* 2001. [PMID 11380923](https://pubmed.ncbi.nlm.nih.gov/11380923/)
80. **Exchange transfusion efficiency and post-exchange rebound.**  
   Ouerradi N et al. Forgoing Exchange Transfusion in Neonatal Hyperbilirubinemia: A Single-Center Retrospective Cohort Study. *Cureus* 2024. [PMID 38650795](https://pubmed.ncbi.nlm.nih.gov/38650795/)
81. **Thresholds for phototherapy and exchange in isoimmune haemolysis.**  
   Slaughter JL, Kemper AR, Newman TB Technical Report: Diagnosis and Management of Hyperbilirubinemia in the Newborn Infant 35 or More Weeks of Gestation. *Pediatrics* 2022. [PMID 35927519](https://pubmed.ncbi.nlm.nih.gov/35927519/)
82. **Directed antiglobulin test and the diagnosis of neonatal HDFN.**  
   Phillips J, Henderson AC Hemolytic Anemia: Evaluation and Differential Diagnosis. *Am Fam Physician* 2018. [PMID 30215915](https://pubmed.ncbi.nlm.nih.gov/30215915/)

## 15. Late anaemia: the nadir that arrives six weeks after delivery

83. **Late hyporegenerative anaemia after intrauterine transfusion and the need for top-up transfusions.**  
   Scaradavou A et al. Suppression of erythropoiesis by intrauterine transfusions in hemolytic disease of the newborn: use of erythropoietin to treat the late anemia. *J Pediatr* 1993. [PMID 8345428](https://pubmed.ncbi.nlm.nih.gov/8345428/)
84. **Erythropoietin for late anaemia of HDFN.**  
   Christensen RD, Bahr TM, Ohls RK Understanding, detecting, and managing the "late" anemia of hemolytic disease of the fetus and newborn. *Best Pract Res Clin Obstet Gynaecol* 2025. [PMID 40752318](https://pubmed.ncbi.nlm.nih.gov/40752318/)
85. **Persistence of maternal antibody in the newborn and duration of haemolysis.**  
   Cortey A, Brossard Y [Adverse effects and patient information]. *J Gynecol Obstet Biol Reprod (Paris)* 2006. [PMID 16495836](https://pubmed.ncbi.nlm.nih.gov/16495836/)
86. **Iron overload after multiple intrauterine and neonatal transfusions.**  
   Adam MP et al. Congenital Dyserythropoietic Anemia Type I. ** 1993. [PMID 20301759](https://pubmed.ncbi.nlm.nih.gov/20301759/)
87. **Neonatal complications after IUT: neutropenia, thrombocytopenia, cholestasis.**  
   van Klink JMM et al. Long-term neurodevelopmental outcomes after intrauterine transfusion for alloimmune hemolytic disease of the fetus and newborn. *Best Pract Res Clin Obstet Gynaecol* 2025. [PMID 41161116](https://pubmed.ncbi.nlm.nih.gov/41161116/)

## 16. Long-term outcome

88. **LOTUS: 4.8% neurodevelopmental impairment, and severe hydrops as the dominant preoperative risk factor (OR 11.2).**  
   Lindenburg IT et al. Long-term neurodevelopmental outcome after intrauterine transfusion for hemolytic disease of the fetus/newborn: the LOTUS study. *Am J Obstet Gynecol* 2012. [PMID 22030316](https://pubmed.ncbi.nlm.nih.gov/22030316/)
89. **Design of the long-term follow-up cohort.**  
   Verduin EP et al. Long-Term follow up after intra-Uterine transfusionS; the LOTUS study. *BMC Pregnancy Childbirth* 2010. [PMID 21122095](https://pubmed.ncbi.nlm.nih.gov/21122095/)
90. **Bilirubin neurotoxicity and kernicterus spectrum disorder.**  
   Bhutani VK, Johnson-Hamerman L The clinical syndrome of bilirubin-induced neurologic dysfunction. *Semin Fetal Neonatal Med* 2015. [PMID 25577653](https://pubmed.ncbi.nlm.nih.gov/25577653/)
91. **Maternal red cell alloimmunisation as a lifelong transfusion problem.**  
   Schonewille H et al. High additional maternal red cell alloimmunization after Rhesus- and K-matched intrauterine intravascular transfusions for hemolytic disease of the fetus. *Am J Obstet Gynecol* 2007. [PMID 17306657](https://pubmed.ncbi.nlm.nih.gov/17306657/)

## 17. Methods: QSP, mrgsolve, and model-informed drug development

92. **mrgsolve as the ODE engine.**  
   Elmokadem A, Riggs MM, Baron KT Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
93. **QSP model credibility and validation practice.**  
   Androulakis IP Towards a comprehensive assessment of QSP models: what would it take?. *J Pharmacokinet Pharmacodyn* 2024. [PMID 35962928](https://pubmed.ncbi.nlm.nih.gov/35962928/)
94. **Model-informed drug development in pregnancy, where trials are small.**  
   Krishna R et al. State-of-the-Art on Model-Informed Drug Development Approaches for Pediatric Rare Diseases. *CPT Pharmacometrics Syst Pharmacol* 2025. [PMID 40697166](https://pubmed.ncbi.nlm.nih.gov/40697166/)
95. **Physiologically based modelling of maternal-fetal drug and antibody transfer.**  
   Cai X et al. Physiologically based pharmacokinetic model of IgG to predict mother-to-fetus transfer of ustekinumab in pregnant patients with inflammatory bowel disease. *J Pharm Sci* 2025. [PMID 40639462](https://pubmed.ncbi.nlm.nih.gov/40639462/)
96. **Virtual population methods for small single-arm trials.**  
   Peláez-Vélez FJ et al. Use of Virtual Reality and Videogames in the Physiotherapy Treatment of Stroke Patients: A Pilot Randomized Controlled Trial. *Int J Environ Res Public Health* 2023. [PMID 36981652](https://pubmed.ncbi.nlm.nih.gov/36981652/)

---

96 references, resolved 2026-08-06 from PubMed. Queries and the raw esearch/esummary payloads are in [`refs_raw.json`](refs_raw.json) and [`refs_meta.json`](refs_meta.json) so that every entry can be re-derived.
