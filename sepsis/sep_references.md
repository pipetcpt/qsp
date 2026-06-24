# Sepsis QSP Model — References

Comprehensive reference list for the Sepsis Quantitative Systems Pharmacology (QSP) model,
organized by topic. Covers landmark clinical trials, guidelines, pathophysiology, cytokine biology,
coagulation, hemodynamics, antibiotic PK/PD, mathematical/QSP models, and biomarkers.

---

## 1. Landmark Clinical Trials

1. **Rivers E, Nguyen B, Havstad S, Ressler J, Muzzin A, Knoblich B, Peterson E, Tomlanovich M; Early Goal-Directed Therapy Collaborative Group**. Early Goal-Directed Therapy in the Treatment of Severe Sepsis and Septic Shock. *New England Journal of Medicine*. 2001;345(19):1368–1377. PMID: 11794169
   > *Relevance: Seminal EGDT trial establishing hemodynamic resuscitation targets (CVP, MAP, ScvO₂); parameterizes the hemodynamic sub-model including fluid, vasopressor, and transfusion decision thresholds.*

2. **ProCESS Investigators; Yealy DM, Kellum JA, Huang DT, Barnato AE, Weissfeld LA, Pike F, Terndrup T, Wang HE, Hou PC, LoVecchio F, Filbin MR, Shapiro NI, Angus DC**. A Randomized Trial of Protocol-Based Care for Early Septic Shock. *New England Journal of Medicine*. 2014;370(18):1683–1693. PMID: 24635773
   > *Relevance: Showed protocolized resuscitation did not improve survival over usual care; informs model uncertainty around ScvO₂-guided therapy effectiveness.*

3. **ARISE Investigators; ANZICS Clinical Trials Group; Peake SL, Delaney A, Bailey M, Bellomo R, Cameron PA, Cooper DJ, Higgins AM, Holdgate A, Howe BD, Webb SA, Williams P**. Goal-Directed Resuscitation for Patients with Early Septic Shock. *New England Journal of Medicine*. 2014;371(16):1496–1506. PMID: 25272316
   > *Relevance: Confirms EGDT non-superiority across an independent cohort; provides mortality rates and hemodynamic response data for model validation.*

4. **Mouncey PR, Osborn TM, Power GS, Harrison DA, Sadique MZ, Grieve RD, Jahan R, Harvey SE, Bell D, Bion JF, Coats TJ, Singer M, Young JD, Rowan KM; ProMISe Trial Investigators**. Trial of Early, Goal-Directed Resuscitation for Septic Shock. *New England Journal of Medicine*. 2015;372(14):1301–1311. PMID: 25776532
   > *Relevance: European counterpart to ProCESS/ARISE; mortality and organ-failure endpoint data anchor the model's clinical outcome calibration.*

5. **Russell JA, Walley KR, Singer J, Gordon AC, Hébert PC, Cooper DJ, Holmes CL, Mehta S, Granton JT, Storms MM, Cook DJ, Presneill JJ, Ayers D; VASST Investigators**. Vasopressin versus Norepinephrine Infusion in Patients with Septic Shock. *New England Journal of Medicine*. 2008;358(9):877–887. PMID: 18305265
   > *Relevance: VASST trial; provides dose–response parameters for vasopressin and norepinephrine in the vasopressor PK/PD sub-model.*

6. **Venkatesh B, Finfer S, Cohen J, Rajbhandari D, Arabi Y, Bellomo R, Billot L, Correa M, Glass P, Harward M, Joyce C, Li Q, McArthur C, Perner A, Rhodes A, Thompson K, Webb S, Myburgh J; ADRENAL Trial Investigators and the Australian–New Zealand Intensive Care Society Clinical Trials Group**. Adjunctive Glucocorticoid Therapy in Patients with Septic Shock. *New England Journal of Medicine*. 2018;378(9):797–808. PMID: 29347874
   > *Relevance: ADRENAL trial; hydrocortisone effect on shock reversal and time-to-vasopressor cessation informs the corticosteroid PD module.*

7. **Annane D, Renault A, Brun-Buisson C, Megarbane B, Quenot JP, Siami S, Cariou A, Forceville X, Schwebel C, Martin C, Timsit JF, Misset B, Ali Benali M, Colin G, Souweine B, Asehnoune K, Mercier E, Chimot L, Charpentier C, François B, Boulain T, Brault C, Mayaux J, Nseir S, Dehoux M, Helms J, Rigaud JP, Azoulay E, Santré C, Morin-Longuet P, Moreau D, Putaud C, Avargues P, Leblanc G, Morel J, Sapin E, Guérin C, Algotsson L, Cariou A, Jaber S, Galoisy-Guibal L, Faller JP, Czosnyka M, Mira JP, Dhainaut JF; APROCCHSS Trial Investigators**. Hydrocortisone plus Fludrocortisone for Adults with Septic Shock. *New England Journal of Medicine*. 2018;378(9):809–818. PMID: 29373074
   > *Relevance: APROCCHSS trial; combined corticosteroid/mineralocorticoid therapy; informs adrenal-axis feedback node in the neuroendocrine sub-model.*

8. **Bernard GR, Vincent JL, Laterre PF, LaRosa SP, Dhainaut JF, Lopez-Rodriguez A, Steingrub JS, Garber GE, Helterbrand JD, Ely EW, Fisher CJ Jr; Recombinant Human Protein C Worldwide Evaluation in Severe Sepsis (PROWESS) Study Group**. Efficacy and Safety of Recombinant Human Activated Protein C for Severe Sepsis. *New England Journal of Medicine*. 2001;344(10):699–709. PMID: 11236773
   > *Relevance: PROWESS trial (drotrecogin alfa); provides anticoagulant-pathway effect estimates for the protein C node in the coagulation sub-model.*

9. **SAFE Study Investigators; Finfer S, Bellomo R, Boyce N, French J, Myburgh J, Norton R**. A Comparison of Albumin and Saline for Fluid Resuscitation in the Intensive Care Unit. *New England Journal of Medicine*. 2004;350(22):2247–2256. PMID: 15163774
   > *Relevance: SAFE trial; colloid vs. crystalloid osmotic effects on intravascular volume and microvascular leak parameters in the fluid-balance module.*

10. **Semler MW, Self WH, Wanderer JP, Ehrenfeld JM, Wang L, Byrne DW, Stollings JL, Kumar AB, Hughes CG, Hernandez A, Guillamondegui OD, May AK, Weavind L, Casey JD, Siew ED, Shaw AD, Bernard GR, Rice TW; SMART Investigators and the Pragmatic Critical Care Research Group**. Balanced Crystalloids versus Saline in Critically Ill Adults. *New England Journal of Medicine*. 2018;378(9):829–839. PMID: 29485925
   > *Relevance: SMART trial; hyperchloremic acidosis from normal saline informs the chloride-bicarbonate acid–base and renal sub-model nodes.*

11. **Holst LB, Haase N, Wetterslev J, Wernerman J, Guttormsen AB, Karlsson S, Johansson PI, Åneman A, Vang ML, Winding R, Nebrich L, Nibro HL, Rasmussen BS, Lauridsen JR, Nielsen JS, Oldner A, Pettilä V, Cronhjort MB, Andersen LH, Pedersen UG, Reiter N, Wiis J, White JO, Russell L, Thornberg KJ, Hjortrup PB, Müller RG, Møller MH, Steensen M, Tjäder I, Kilsand K, Odeberg-Wernerman S, Sjøbø B, Bundgaard H, Thyø MA, Lodahl D, Mærkedahl R, Albeck C, Illum D, Kruse M, Winkel P, Perner A; TRISS Trial Group; Scandinavian Critical Care Trials Group**. Lower versus Higher Hemoglobin Threshold for Transfusion in Septic Shock. *New England Journal of Medicine*. 2014;371(15):1381–1391. PMID: 25270275
    > *Relevance: TRISS trial; establishes transfusion trigger thresholds (7 vs 9 g/dL) for the oxygen-delivery component of the hemodynamic sub-model.*

---

## 2. Surviving Sepsis Campaign Guidelines

12. **Rhodes A, Evans LE, Alhazzani W, Levy MM, Antonelli M, Ferrer R, Kumar A, Sevransky JE, Sprung CL, Nunnally ME, Rochwerg B, Rubenfeld GD, Angus DC, Annane D, Beale RJ, Bellinghan GJ, Bernard GR, Chiche JD, Coopersmith C, De Backer DP, French CJ, Fujishima S, Gerlach H, Hidalgo JL, Hollenberg SM, Jones AE, Karnad DR, Kleinpell RM, Koh Y, Lisboa TC, Machado FR, Marini JJ, Marshall JC, Mazuski JE, McIntyre LA, McLean AS, Mehta S, Moreno RP, Myburgh J, Navalesi P, Nishida O, Osborn TM, Perner A, Plunkett CM, Ranieri M, Schorr CA, Seckel MA, Seymour CW, Shieh L, Shukri KA, Simpson SQ, Singer M, Thompson BT, Townsend SR, Van der Poll T, Vincent JL, Wiersinga WJ, Zimmerman JL, Dellinger RP**. Surviving Sepsis Campaign: International Guidelines for Management of Sepsis and Septic Shock 2016. *Intensive Care Medicine*. 2017;43(3):304–377. PMID: 28101605
    > *Relevance: SSC 2016 guidelines; primary source for resuscitation bundle components, antibiotic timing, vasopressor targets, and corticosteroid indications modeled in treatment scenarios.*

13. **Evans L, Rhodes A, Alhazzani W, Antonelli M, Coopersmith CM, French C, Machado FR, Mcintyre L, Ostermann M, Prescott HC, Schorr C, Simpson S, Wiersinga WJ, Alshamsi F, Angus DC, Arabi Y, Azevedo L, Beale R, Beilman G, Belley-Cote E, Burry L, Cecconi M, Centofanti J, Coz Yataco A, De Waele J, Dellinger RP, Doi K, Du B, Estenssoro E, Ferrer R, Gomersall C, Hodgson C, Hylander Møller M, Iwashyna T, Jacob S, Kleinpell R, Klompas M, Koh Y, Kumar A, Kwizera A, Lobo S, Masur H, McGloughlin S, Mehta S, Mehta Y, Mer M, Nunnally M, Oczkowski S, Osborn T, Papathanassoglou E, Perner A, Puskarich M, Roberts J, Schweickert W, Seckel M, Sevransky J, Sprung CL, Welte T, Zimmerman J, Levy M**. Surviving Sepsis Campaign: International Guidelines for Management of Sepsis and Septic Shock 2021. *Critical Care Medicine*. 2021;49(11):e1063–e1143. PMID: 34605781
    > *Relevance: SSC 2021 update; revised fluid resuscitation recommendations, vasopressor dosing, and antibiotic stewardship guidance reflected in the model's treatment decision nodes.*

14. **Levy MM, Evans LE, Rhodes A**. The Surviving Sepsis Campaign Bundle: 2018 Update. *Intensive Care Medicine*. 2018;44(6):925–928. PMID: 29675566
    > *Relevance: Hour-1 Bundle specification (lactate, blood cultures, broad-spectrum antibiotics, fluids, vasopressors); defines the time-to-treatment parameters in the early resuscitation module.*

---

## 3. Pathophysiology & Mechanisms

15. **Hotchkiss RS, Karl IE**. The Pathophysiology and Treatment of Sepsis. *New England Journal of Medicine*. 2003;348(2):138–150. PMID: 12519925
    > *Relevance: Foundational review of pro-inflammatory, anti-inflammatory, and apoptotic pathways; provides mechanistic architecture for the innate immunity and immune-suppression nodes.*

16. **Singer M, Deutschman CS, Seymour CW, Shankar-Hari M, Annane D, Bauer M, Bellomo R, Bernard GR, Chiche JD, Coopersmith CM, Hotchkiss RS, Levy MM, Marshall JC, Martin GS, Opal SM, Rubenfeld GD, van der Poll T, Vincent JL, Angus DC**. The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). *JAMA*. 2016;315(8):801–810. PMID: 26903338
    > *Relevance: Sepsis-3 organ-dysfunction definition; SOFA-based clinical endpoint definition used for model output validation and parameter fitting.*

17. **Vincent JL, Moreno R, Takala J, Willatts S, De Mendonça A, Bruining H, Reinhart CK, Suter PM, Thijs LG**. The SOFA (Sepsis-related Organ Failure Assessment) Score to Describe Organ Dysfunction/Failure. On Behalf of the Working Group on Sepsis-Related Problems of the European Society of Intensive Care Medicine. *Intensive Care Medicine*. 1996;22(7):707–710. PMID: 8844239
    > *Relevance: SOFA score components (respiratory, coagulation, liver, cardiovascular, CNS, renal) directly map onto the model's six-organ-system outputs.*

18. **Angus DC, van der Poll T**. Severe Sepsis and Septic Shock. *New England Journal of Medicine*. 2013;369(9):840–851. PMID: 23984731
    > *Relevance: Comprehensive review of immune dysregulation, endothelial damage, and mitochondrial dysfunction; supports parameterization of the inflammatory cascade and microcirculation modules.*

19. **van der Poll T, van de Veerdonk FL, Scicluna BP, Netea MG**. The Immunopathology of Sepsis and Potential Therapeutic Targets. *Nature Reviews Immunology*. 2017;17(7):407–420. PMID: 28436424
    > *Relevance: Detailed innate/adaptive immune crosstalk; T-cell exhaustion and immunosuppression dynamics inform the late-phase immune-suppression state variables.*

20. **Seymour CW, Liu VX, Iwashyna TJ, Brunkhorst FM, Rea TD, Scherag A, Rubenfeld G, Kahn JM, Shankar-Hari M, Singer M, Deutschman CS, Escobar GJ, Angus DC**. Assessment of Clinical Criteria for Sepsis: For the Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). *JAMA*. 2016;315(8):762–774. PMID: 26903335
    > *Relevance: Clinical validation of qSOFA as a screening score; provides population-level transition probabilities for the model's patient-state definitions.*

21. **Gotts JE, Matthay MA**. Sepsis: Pathophysiology and Clinical Management. *BMJ*. 2016;353:i1585. PMID: 27217054
    > *Relevance: Overview of alveolar-capillary barrier disruption and ARDS in sepsis; underpins the pulmonary compartment equations and oxygenation endpoints.*

---

## 4. Cytokine Biology in Sepsis

22. **Beutler B, Cerami A**. Cachectin: More Than a Tumor Necrosis Factor. *New England Journal of Medicine*. 1987;316(7):379–385. PMID: 3027565
    > *Relevance: Original description of TNF-α (cachectin) as a primary mediator of septic shock; provides the basis for the TNF-α production and receptor-signaling equations.*

23. **Dinarello CA**. Infection, Fever, and Exogenous and Endogenous Pyrogens: Some Concepts Have Changed. *Journal of Endotoxin Research*. 2004;10(4):201–222. PMID: 15373964
    > *Relevance: IL-1β role in fever, coagulation activation, and endothelial dysfunction; parameterizes IL-1 production rates and downstream effector pathways in the model.*

24. **Medzhitov R**. Origin and Physiological Roles of Inflammation. *Nature*. 2008;454(7203):428–435. PMID: 18650913
    > *Relevance: Conceptual framework for pattern recognition receptors (TLR/NLR) and initiation of innate inflammatory response; underpins the PAMP-sensing nodes at the model's entry point.*

25. **Liu Q, Zhou YH, Yang ZQ**. The Cytokine Storm of COVID-19: A Flashpoint of Sepsis Research. *Lancet*. 2020;395(10229):1033–1034. PMID: 32192578
    > *Relevance: Cytokine storm dynamics and IL-6/IL-10 feedback loops; informs the anti-inflammatory counter-regulatory branch of the cytokine sub-model.*

26. **Schulte W, Bernhagen J, Bucala R**. Cytokines in Sepsis: Potent Immunoregulators and Potential Therapeutic Targets — An Updated View. *Mediators of Inflammation*. 2013;2013:165974. PMID: 23853427
    > *Relevance: Cytokine biomarker kinetics (TNF-α, IL-6, IL-8, IL-10, MIF); provides concentration–time profiles used to fit model output to clinical observations.*

---

## 5. Coagulation & DIC

27. **Levi M, Ten Cate H**. Disseminated Intravascular Coagulation. *New England Journal of Medicine*. 1999;341(8):586–592. PMID: 10451465
    > *Relevance: DIC pathophysiology: simultaneous coagulation activation and fibrinolysis inhibition; foundation for the coagulation cascade and anticoagulant pathway sub-model.*

28. **Gando S, Levi M, Toh CH**. Disseminated Intravascular Coagulation. *Nature Reviews Disease Primers*. 2016;2:16037. PMID: 27250996
    > *Relevance: Immunothrombosis and endothelial crosstalk; connects inflammatory cytokine state variables to thrombin generation rate equations.*

29. **Iba T, Levy JH, Warkentin TE, Thachil J, van der Poll T, Levi M; Scientific and Standardization Committee on DIC, and the Scientific and Standardization Committee on Perioperative and Critical Care of the International Society on Thrombosis and Haemostasis**. Diagnosis and Management of Sepsis-Induced Coagulopathy and Disseminated Intravascular Coagulation. *Journal of Thrombosis and Haemostasis*. 2019;17(11):1989–1994. PMID: 31410983
    > *Relevance: SIC (sepsis-induced coagulopathy) diagnostic criteria and management; DIC scoring integrated into organ-failure composite endpoint calculations.*

---

## 6. Hemodynamics & Vasopressors

30. **Hollenberg SM, Ahrens TS, Annane D, Astiz ME, Chalfin DB, Dasta JF, Heard SO, Martin C, Napolitano LM, Reilly JM, Richard C, Svensson LG, Totaro R, Vincent JL, Zimmerman JL**. Practice Parameters for Hemodynamic Support of Sepsis in Adult Patients: 2004 Update. *Critical Care Medicine*. 2004;32(9):1928–1948. PMID: 15343024
    > *Relevance: Comprehensive hemodynamic monitoring and vasopressor practice parameters; provides the target MAP, CO/CI, and SVRI thresholds for the cardiovascular sub-model.*

31. **Levy B, Dusang B, Annane D, Bauer C, Cholley B, Chouihed T, Eon B, Francoz C, Garçon P, Houessou G, Jolly A, Kaci R, Levy P, Miette A, Mira JP, Pastré J, Pavot A, Piquet J, Plaisance P, Ricard JD, Richard JC, Rolin N, Teboul JL, Ternacle J, Titeca-Beauport D, Valette X, Vieillard-Baron A, Zafrani L, Meziani F**. Comparison of Norepinephrine–Dobutamine to Epinephrine for Hemodynamics, Lactate Metabolism, and Organ Function Variables in Cardiogenic Shock. *Critical Care Medicine*. 2011;39(3):450–455. PMID: 21057313
    > *Relevance: Norepinephrine dose–response and catecholamine effects on lactate metabolism; parameterizes the vasopressor-dose–MAP relationship and lactate clearance equations.*

32. **Russell JA, Walley KR, Singer J, Gordon AC, Hébert PC, Cooper DJ, Holmes CL, Mehta S, Granton JT, Storms MM, Cook DJ, Presneill JJ, Ayers D; VASST Investigators**. Vasopressin versus Norepinephrine Infusion in Patients with Septic Shock. *New England Journal of Medicine*. 2008;358(9):877–887. PMID: 18305265
    > *Relevance: VASST pharmacodynamic data for vasopressin (0.03 U/min fixed dose); V1 receptor-mediated vasoconstriction kinetics in the vasopressor module.*

33. **De Backer D, Biston P, Devriendt J, Madl C, Chochrad D, Aldecoa C, Brasseur A, Defrance P, Gottignies P, Vincent JL; SOAP II Investigators**. Comparison of Dopamine and Norepinephrine in the Treatment of Shock. *New England Journal of Medicine*. 2010;362(9):779–789. PMID: 20200382
    > *Relevance: Comparative vasopressor PD data (dopamine vs. norepinephrine); arrhythmia risk and mortality data inform the adverse-effect branch of the vasopressor sub-model.*

---

## 7. Antibiotic PK/PD in Sepsis

34. **Roberts JA, Paul SK, Akova M, Bassetti M, De Waele JJ, Dimopoulos G, Kaukonen KM, Koulenti D, Martin C, Montravers P, Rello J, Rhodes A, Starr T, Wallis SC, Lipman J; DALI Study**. DALI: Defining Antibiotic Levels in Intensive Care Unit Patients: Are Current ß-Lactam Antibiotic Doses Sufficient for Critically Ill Patients? *Clinical Infectious Diseases*. 2014;58(8):1072–1083. PMID: 24457344
    > *Relevance: Real-world PK data from ICU patients showing subtherapeutic beta-lactam exposures; motivates the augmented-clearance correction in the antibiotic PK sub-model.*

35. **Kumar A, Roberts D, Wood KE, Light B, Parrillo JE, Sharma S, Suppes R, Feinstein D, Zanotti S, Taiberg L, Gurka D, Kumar A, Cheang M**. Duration of Hypotension Before Initiation of Effective Antimicrobial Therapy Is the Critical Determinant of Survival in Human Septic Shock. *Critical Care Medicine*. 2006;34(6):1589–1596. PMID: 16625125
    > *Relevance: Every hour delay in antibiotics associated with ~7% mortality increase; provides the time-to-antibiotic effectiveness function for the bacterial kill-rate module.*

36. **Udy AA, Roberts JA, Boots RJ, Paterson DL, Lipman J**. Augmented Renal Clearance: Implications for Antibacterial Dosing in the Critically Ill. *Clinical Pharmacokinetics*. 2010;49(1):1–16. PMID: 20000886
    > *Relevance: Augmented renal clearance (ARC) increases drug elimination in hyperdynamic sepsis; provides CrCl-adjustment equations for the antibiotic clearance model.*

37. **Craig WA**. Pharmacokinetic/Pharmacodynamic Parameters: Rationale for Antibacterial Dosing of Mice and Men. *Clinical Infectious Diseases*. 1998;26(1):1–10. PMID: 9455502
    > *Relevance: Classic exposition of time-dependent (beta-lactams: %T>MIC), concentration-dependent (aminoglycosides: Cmax/MIC), and AUC-dependent (fluoroquinolones: AUC/MIC) PK/PD indices; defines killing-rate equations for all three antibiotic classes in the model.*

38. **Rybak MJ, Lomaestro BM, Rotschafer JC, Moellering RC Jr, Craig WA, Billeter M, Dalovisio JR, Levine DP**. Therapeutic Monitoring of Vancomycin in Adult Patients: A Consensus Review of the American Society of Health-System Pharmacists, the Infectious Diseases Society of America, and the Society of Infectious Diseases Pharmacists. *American Journal of Health-System Pharmacy*. 2009;66(1):82–98. PMID: 19106348
    > *Relevance: Vancomycin AUC/MIC target (400–600 mg·h/L) and nephrotoxicity thresholds; parameterizes the vancomycin PK/PD and renal safety constraints in the antibiotic dosing module.*

---

## 8. QSP / Mathematical Models of Sepsis

39. **Chow CC, Clermont G, Kumar R, Lagoa C, Tawadrous Z, Gallo D, Betten B, Bartels J, Constantine G, Fink MP, Billiar TR, Vodovotz Y**. The Acute Inflammatory Response in Diverse Shock States. *Shock*. 2005;24(1):74–84. PMID: 15988324
    > *Relevance: Foundational ODE model of acute inflammation (pathogen load, early/late mediators, damage); the structural skeleton for the present QSP model's inflammatory state variables.*

40. **An G, Nieman G, Vodovotz Y**. Computational and Systems Biology in Trauma and Sepsis: Current State and Future Perspectives. *International Journal of Burns and Trauma*. 2012;2(1):1–10. PMID: 22928162
    > *Relevance: Agent-based modeling of macrophage–neutrophil–epithelial interactions; provides cellular-level mechanistic insights that inform the cell-population sub-model.*

41. **Day J, Rubin J, Vodovotz Y, Chow CC, Reynolds A, Clermont G**. A Reduced Mathematical Model of the Acute Inflammatory Response: II. Capturing Scenarios of Repeated Endotoxin Administration. *Journal of Theoretical Biology*. 2006;242(1):237–256. PMID: 16701699
    > *Relevance: Endotoxin tolerance and repeated-insult modeling; informs LPS dose–response nonlinearity and macrophage desensitization dynamics in the PAMP-sensing module.*

42. **Clermont G, Bartels J, Kumar R, Constantine G, Vodovotz Y, Chow C**. In Silico Design of Clinical Trials: A Method Coming of Age. *Critical Care Medicine*. 2004;32(10):2061–2070. PMID: 15483413
    > *Relevance: Demonstrates in silico trial design for anti-inflammatory therapies in sepsis; the computational framework concept directly underpins the model's virtual-patient cohort simulation capability.*

43. **Vodovotz Y, Clermont G, Chow C, An G**. Mathematical Models of the Acute Inflammatory Response. *Current Opinion in Critical Care*. 2004;10(5):383–390. PMID: 15385762
    > *Relevance: Review of ODE-based acute inflammation models; provides benchmark validation criteria and explains the mathematical rationale for the three-variable (pathogen, early mediator, late mediator) core.*

44. **Torres A, Ferrer M, Badia JR**. Treatment Guidelines and Outcomes of Hospital-Acquired and Ventilator-Associated Pneumonia. *Clinical Infectious Diseases*. 2010;51(Suppl 1):S48–S53. PMID: 20597671
    > *Relevance: Mechanistic model structure for hospital-acquired infection progression; provides transition rates between infection severity states used in the bacterial-clearance sub-model.*

45. **Cockrell C, An G**. Sepsis Reconsidered: Identifying Novel Metrics for Behavioral Landscape Characterization with a High-Performance Computing Implementation of an Agent-Based Model. *Journal of Theoretical Biology*. 2017;430:157–168. PMID: 28716441
    > *Relevance: High-performance agent-based sepsis simulation; provides emergent behavior benchmarks (e.g., bifurcation between survival/death attractors) against which the ODE model's phase-plane behavior is validated.*

---

## 9. Biomarkers

46. **Becker KL, Snider R, Nylen ES**. Procalcitonin Assay in Systemic Inflammation, Infection, and Sepsis: Clinical Utility and Limitations. *Critical Care Medicine*. 2008;36(3):941–952. PMID: 18431284
    > *Relevance: Procalcitonin (PCT) kinetics during infection and clearance; PCT production and half-life equations calibrate the biomarker output module for antibiotic stewardship scenarios.*

47. **Brunkhorst FM, Engel C, Bloos F, Meier-Hellmann A, Ragaller M, Weiler N, Moerer O, Gruendling M, Oppert M, Grond S, Olthoff D, Jaschinski U, John S, Rossaint R, Welte T, Schaefer M, Kern P, Kuhnt E, Kiehntopf M, Hartog C, Natanson C, Loeffler M, Reinhart K; German Competence Network Sepsis (SepNet)**. Intensive Insulin Therapy and Pentastarch Resuscitation in Severe Sepsis. *New England Journal of Medicine*. 2008;358(2):125–139. PMID: 18184958
    > *Relevance: Lactate and PCT as severity biomarkers alongside intensive insulin therapy; provides paired lactate–mortality data for calibrating the lactate compartment and glycemic control sub-model.*

48. **Kellum JA, Pike F, Yealy DM, Huang DT, Shapiro NI, Angus DC; ProCESS Investigators**. Relationship Between Alternative Resuscitation Strategies, Host Response, and Injury Biomarkers, and Outcome in Septic Shock: Analysis of the ProCESS Trial. *Critical Care Medicine*. 2017;45(1):438–445. PMID: 27643972
    > *Relevance: Lactate clearance trajectories in ProCESS trial; provides lactate clearance rate constants (k_lac) for the metabolic sub-model and its link to mortality prediction.*

49. **Wacker C, Prkno A, Brunkhorst FM, Schlattmann P**. Procalcitonin as a Diagnostic Marker for Sepsis: A Systematic Review and Meta-Analysis. *Lancet Infectious Diseases*. 2013;13(5):426–435. PMID: 23375419
    > *Relevance: Pooled PCT sensitivity (77%) and specificity (79%) for sepsis diagnosis; provides Bayesian prior probabilities for the PCT-based diagnostic node in the model's virtual-patient initialization module.*

50. **Shapiro NI, Howell MD, Talmor D, Nathanson LA, Lisbon A, Wolfe RE, Weiss JW**. Serum Lactate as a Predictor of Mortality in Emergency Department Patients with Infection. *Annals of Emergency Medicine*. 2005;45(5):524–528. PMID: 15855951
    > *Relevance: Lactate stratification (< 2, 2–4, > 4 mmol/L) and 60-day mortality; provides the lactate-to-outcome probability mapping for the model's risk-stratification output.*

---

## Additional Key References

51. **Annane D, Sébille V, Charpentier C, Bollaert PE, François B, Korach JM, Capellier G, Cohen Y, Azoulay E, Troché G, Chaumet-Riffaud P, Bellissant E**. Effect of Treatment with Low Doses of Hydrocortisone and Fludrocortisone on Mortality in Patients with Septic Shock. *JAMA*. 2002;288(7):862–871. PMID: 12186604
    > *Relevance: Early RCT of low-dose corticosteroids; baseline cortisol response curves and survival data used to parameterize adrenal insufficiency sub-model prior to ADRENAL/APROCCHSS.*

52. **Sprung CL, Annane D, Keh D, Moreno R, Singer M, Freivogel K, Weiss YG, Benbenishty J, Kalenka A, Forst H, Laterre PF, Reinhart K, Cuthbertson BH, Payen D, Briegel J; CORTICUS Study Group**. Hydrocortisone Therapy for Patients with Septic Shock. *New England Journal of Medicine*. 2008;358(2):111–124. PMID: 18184957
    > *Relevance: CORTICUS trial negative result for hydrocortisone; used together with ADRENAL/APROCCHSS to define subgroup-specific steroid responder probability in the model.*

53. **Reinhart K, Meisner M**. Biomarkers in the Critically Ill Patient: Procalcitonin. *Critical Care Clinics*. 2011;27(2):253–263. PMID: 21440202
    > *Relevance: PCT pharmacokinetics (peak at 24 h, t½ ~24 h), correlation with bacteremia and fungemia; informs the PCT differential equation parameterization and antibiotic de-escalation trigger logic.*

54. **Vincent JL, De Backer D**. Circulatory Shock. *New England Journal of Medicine*. 2013;369(18):1726–1734. PMID: 24171518
    > *Relevance: Classification of distributive, cardiogenic, hypovolemic, and obstructive shock; provides the clinical phenotype branching logic in the model's hemodynamic state-transition structure.*

55. **Marik PE, Khangoora V, Rivera R, Hooper MH, Catravas J**. Hydrocortisone, Vitamin C, and Thiamine for the Treatment of Severe Sepsis and Septic Shock: A Retrospective Before-After Study. *Chest*. 2017;151(6):1229–1238. PMID: 27940189
    > *Relevance: HAT (hydrocortisone-ascorbate-thiamine) protocol pilot data; ascorbate antioxidant effects on vascular permeability and catecholamine sensitivity inform the adjunctive-therapy nodes in the treatment scenario module.*

---

*File generated: 2026-06-24*
*Model version: Sepsis QSP v1.0*
*Total references: 55*
