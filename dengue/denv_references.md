# Dengue / Severe Dengue — References

References for the QSP model in this directory. Grouped by the part of the
model each one constrains. Where a reference fixed a parameter rather than
merely motivating a mechanism, the parameter is named.

---

## 1. Burden, classification and natural history

1. Bhatt S, et al. *The global distribution and burden of dengue.* Nature 2013;496:504–7. — <https://pubmed.ncbi.nlm.nih.gov/23563266/>
2. Messina JP, et al. *The current and future global distribution and population at risk of dengue.* Nat Microbiol 2019;4:1508–15. — <https://pubmed.ncbi.nlm.nih.gov/31182801/>
3. World Health Organization. *Dengue: Guidelines for Diagnosis, Treatment, Prevention and Control.* Geneva, 2009. — <https://pubmed.ncbi.nlm.nih.gov/23762963/>
4. Simmons CP, Farrar JJ, Nguyen VV, Wills B. *Dengue.* N Engl J Med 2012;366:1423–32. — <https://pubmed.ncbi.nlm.nih.gov/22494122/>
5. Wilder-Smith A, Ooi EE, Horstick O, Wills B. *Dengue.* Lancet 2019;393:350–63. — <https://pubmed.ncbi.nlm.nih.gov/30696575/>
6. Halstead SB. *Dengue.* Lancet 2007;370:1644–52. — <https://pubmed.ncbi.nlm.nih.gov/17993365/>
7. Guzman MG, et al. *Dengue infection.* Nat Rev Dis Primers 2016;2:16055. — <https://pubmed.ncbi.nlm.nih.gov/27534439/>
8. Barniol J, et al. *Usefulness and applicability of the revised dengue case classification.* BMC Infect Dis 2011;11:106. — <https://pubmed.ncbi.nlm.nih.gov/21510901/>
9. Alexander N, et al. *Multicentre prospective study on dengue classification in four South-east Asian and three Latin American countries.* Trop Med Int Health 2011;16:936–48. — <https://pubmed.ncbi.nlm.nih.gov/21624014/>
10. Yacoub S, Wills B. *Predicting outcome from dengue.* BMC Med 2014;12:147. — <https://pubmed.ncbi.nlm.nih.gov/25259615/>

---

## 2. Antibody-dependent enhancement — the centre of the model

The entry factor `E(A) = (1−N)·[1 + (Φ−1)·A/(A+K)]` and its parameters
(`FCROSS = 0.12`, `NT50 = 30`, `HNEUT = 2.5`, `KOPS = 4`, `Φ = 14`) come from
this section.

11. Halstead SB, O'Rourke EJ. *Dengue viruses and mononuclear phagocytes. I. Infection enhancement by non-neutralizing antibody.* J Exp Med 1977;146:201–17. — <https://pubmed.ncbi.nlm.nih.gov/406347/>
12. **Salje H, et al.** *Reconstruction of antibody dynamics and infection histories to evaluate dengue risk.* Nature 2018;557:719–23. — <https://pubmed.ncbi.nlm.nih.gov/29795354/> — **the model's primary ADE calibration target**: risk raised at 1:21–1:80, protection above 1:1280.
13. **Katzelnick LC, et al.** *Antibody-dependent enhancement of severe dengue disease in humans.* Science 2017;358:929–32. — <https://pubmed.ncbi.nlm.nih.gov/29097492/> — the peak-risk titre window measured prospectively.
14. Kliks SC, Nimmanitya S, Nisalak A, Burke DS. *Evidence that maternal dengue antibodies are important in the development of dengue hemorrhagic fever in infants.* Am J Trop Med Hyg 1988;38:411–19. — <https://pubmed.ncbi.nlm.nih.gov/3354774/> — sets the infant validation target.
15. Chau TNB, et al. *Dengue in Vietnamese infants — results of infection-enhancement assays correlate with age-related disease severity.* J Infect Dis 2008;198:516–24. — <https://pubmed.ncbi.nlm.nih.gov/18598194/>
16. Chau TNB, et al. *Clinical and virological features of dengue in Vietnamese infants.* PLoS Negl Trop Dis 2010;4:e657. — <https://pubmed.ncbi.nlm.nih.gov/20405057/>
17. Simmons CP, et al. *Maternal antibody and viral factors in the pathogenesis of dengue virus in infants.* J Infect Dis 2007;196:416–24. — <https://pubmed.ncbi.nlm.nih.gov/17597456/> — maternal IgG decay half-life (`dABH`, t½ 43 d).
18. Clapham HE, et al. *Modelling virus and antibody dynamics during dengue virus infection suggests a role for antibody in virus clearance.* PLoS Comput Biol 2016;12:e1004951. — <https://pubmed.ncbi.nlm.nih.gov/27213556/>
19. Boonnak K, et al. *Cell type specificity and host genetic polymorphisms influence antibody-dependent enhancement of dengue virus infection.* J Virol 2011;85:1671–83. — <https://pubmed.ncbi.nlm.nih.gov/21123382/> — magnitude of Φ.
20. Chareonsirisuthigul T, Kalayanarooj S, Ubol S. *Dengue virus (DENV) antibody-dependent enhancement of infection upregulates the anti-inflammatory, reduces the pro-inflammatory... response in THP-1 cells.* J Gen Virol 2007;88:365–75. — <https://pubmed.ncbi.nlm.nih.gov/17251552/> — intrinsic ADE.
21. Ubol S, Halstead SB. *How innate immune mechanisms contribute to antibody-enhanced viral infections.* Clin Vaccine Immunol 2010;17:1829–35. — <https://pubmed.ncbi.nlm.nih.gov/20876821/> — basis for `KINT` and `PINT`.
22. Chan KR, et al. *Ligation of Fc gamma receptor IIB inhibits antibody-dependent enhancement of dengue virus infection.* Proc Natl Acad Sci USA 2011;108:12479–84. — <https://pubmed.ncbi.nlm.nih.gov/21746897/>
23. Chan KR, et al. *Leukocyte immunoglobulin-like receptor B1 is critical for antibody-dependent dengue.* Proc Natl Acad Sci USA 2014;111:2722–7. — <https://pubmed.ncbi.nlm.nih.gov/24550301/>
24. Dejnirattisai W, et al. *Cross-reacting antibodies enhance dengue virus infection in humans.* Science 2010;328:745–8. — <https://pubmed.ncbi.nlm.nih.gov/20448183/> — poorly neutralising anti-prM antibodies.
25. de Alwis R, et al. *Identification of human neutralizing antibodies that bind to complex epitopes on dengue virions.* Proc Natl Acad Sci USA 2012;109:7439–44. — <https://pubmed.ncbi.nlm.nih.gov/22499787/> — basis for the homotypic avidity term `AVIDN = 3`.
26. Dejnirattisai W, et al. *A new class of highly potent, broadly neutralizing antibodies isolated from viremic patients infected with dengue virus.* Nat Immunol 2015;16:170–7. — <https://pubmed.ncbi.nlm.nih.gov/25501631/> — envelope dimer epitope (EDE).
27. Pierson TC, et al. *The stoichiometry of antibody-mediated neutralization and enhancement of West Nile virus infection.* Cell Host Microbe 2007;1:135–45. — <https://pubmed.ncbi.nlm.nih.gov/18005691/> — the multi-hit basis for `HNEUT = 2.5`.
28. Guzman MG, Alvarez M, Halstead SB. *Secondary infection as a risk factor for dengue hemorrhagic fever/dengue shock syndrome: an historical perspective.* Arch Virol 2013;158:1445–59. — <https://pubmed.ncbi.nlm.nih.gov/23471635/>
29. Endy TP, et al. *Determinants of inapparent and symptomatic dengue infection in a prospective study of primary school children.* PLoS Negl Trop Dis 2011;5:e975. — <https://pubmed.ncbi.nlm.nih.gov/21390156/>
30. Waggoner JJ, et al. *Homotypic dengue virus reinfections in Nicaraguan children.* J Infect Dis 2016;214:986–93. — <https://pubmed.ncbi.nlm.nih.gov/27274182/>

---

## 3. Within-host virology, NS1 and viral kinetics

31. Vaughn DW, et al. *Dengue viremia titer, antibody response pattern, and virus serotype correlate with disease severity.* J Infect Dis 2000;181:2–9. — <https://pubmed.ncbi.nlm.nih.gov/10608744/> — peak-viraemia targets.
32. Libraty DH, et al. *High circulating levels of the dengue virus nonstructural protein NS1 early in dengue illness correlate with the development of dengue hemorrhagic fever.* J Infect Dis 2002;186:1165–8. — <https://pubmed.ncbi.nlm.nih.gov/12355369/> — NS1 scale (`KNS1`).
33. Duyen HTL, et al. *Kinetics of plasma viremia and soluble nonstructural protein 1 concentrations in dengue.* J Infect Dis 2011;203:1292–1300. — <https://pubmed.ncbi.nlm.nih.gov/21335562/> — NS1 and viraemia time courses.
34. Clapham HE, et al. *Within-host viral dynamics of dengue serotype 1 infection.* J R Soc Interface 2014;11:20140094. — <https://pubmed.ncbi.nlm.nih.gov/24829280/> — growth and clearance rates.
35. Ben-Shachar R, Koelle K. *Minimal within-host dengue models highlight the specific roles of the immune response in primary and secondary dengue infections.* J R Soc Interface 2015;12:20140886. — <https://pubmed.ncbi.nlm.nih.gov/25519990/>
36. Nguyen NM, et al. *Host and viral features of human dengue cases shape the population of infected and infectious Aedes aegypti mosquitoes.* Proc Natl Acad Sci USA 2013;110:9072–7. — <https://pubmed.ncbi.nlm.nih.gov/23674683/>
37. Muller DA, Young PR. *The flavivirus NS1 protein: molecular and structural biology, immunology, role in pathogenesis and application as a diagnostic biomarker.* Antiviral Res 2013;98:192–208. — <https://pubmed.ncbi.nlm.nih.gov/23523765/>
38. Alcon S, et al. *Enzyme-linked immunosorbent assay specific to dengue virus type 1 nonstructural protein NS1 reveals circulation of the antigen in the blood during the acute phase of disease.* J Clin Microbiol 2002;40:376–81. — <https://pubmed.ncbi.nlm.nih.gov/11825945/>
39. Rico-Hesse R. *Microevolution and virulence of dengue viruses.* Adv Virus Res 2003;59:315–41. — <https://pubmed.ncbi.nlm.nih.gov/14696333/>
40. OhAinle M, et al. *Dynamics of dengue disease severity determined by the interplay between viral genetics and serotype-specific immunity.* Sci Transl Med 2011;3:114ra128. — <https://pubmed.ncbi.nlm.nih.gov/22190239/>

---

## 4. Innate immunity, interferon and viral antagonism

41. Ashour J, et al. *NS5 of dengue virus mediates STAT2 binding and degradation.* J Virol 2009;83:5408–18. — <https://pubmed.ncbi.nlm.nih.gov/19279106/>
42. Muñoz-Jordán JL, et al. *Inhibition of interferon signaling by dengue virus.* Proc Natl Acad Sci USA 2003;100:14333–8. — <https://pubmed.ncbi.nlm.nih.gov/14612562/>
43. Green AM, Beatty PR, Hadjilaou A, Harris E. *Innate immunity to dengue virus infection and subversion of antiviral responses.* J Mol Biol 2014;426:1148–60. — <https://pubmed.ncbi.nlm.nih.gov/24316047/>
44. Kurane I, et al. *Activation of T lymphocytes in dengue virus infections.* J Clin Invest 1991;88:1473–80. — <https://pubmed.ncbi.nlm.nih.gov/1939640/>
45. Nasirudeen AMA, et al. *RIG-I, MDA5 and TLR3 synergistically play an important role in restriction of dengue virus infection.* PLoS Negl Trop Dis 2011;5:e926. — <https://pubmed.ncbi.nlm.nih.gov/21245912/>
46. Avirutnan P, et al. *Vascular leakage in severe dengue virus infections: a potential role for the nonstructural viral protein NS1 and complement.* J Infect Dis 2006;193:1078–88. — <https://pubmed.ncbi.nlm.nih.gov/16544248/>
47. Avirutnan P, et al. *Antagonism of the complement component C4 by flavivirus nonstructural protein NS1.* J Exp Med 2010;207:793–806. — <https://pubmed.ncbi.nlm.nih.gov/20308361/>
48. Nascimento EJM, et al. *Alternative complement pathway deregulation is correlated with dengue severity.* PLoS One 2009;4:e6782. — <https://pubmed.ncbi.nlm.nih.gov/19707565/>

---

## 5. Adaptive immunity, original antigenic sin and immune complexes

49. Mongkolsapaya J, et al. *Original antigenic sin and apoptosis in the pathogenesis of dengue hemorrhagic fever.* Nat Med 2003;9:921–7. — <https://pubmed.ncbi.nlm.nih.gov/12808447/> — basis for `XAVID = 0.45` and `XTNF = 3.1`.
50. Mongkolsapaya J, et al. *T cell responses in dengue hemorrhagic fever: are cross-reactive T cells suboptimal?* J Immunol 2006;176:3821–9. — <https://pubmed.ncbi.nlm.nih.gov/16517753/>
51. Weiskopf D, Sette A. *T-cell immunity to infection with dengue virus in humans.* Front Immunol 2014;5:93. — <https://pubmed.ncbi.nlm.nih.gov/24639680/>
52. Wrammert J, et al. *Rapid and massive virus-specific plasmablast responses during acute dengue virus infection in humans.* J Virol 2012;86:2911–8. — <https://pubmed.ncbi.nlm.nih.gov/22238318/> — the plasmablast/ASC cascade and its lag (`kMAT`).
53. Zompi S, Harris E. *Original antigenic sin in dengue revisited.* Proc Natl Acad Sci USA 2013;110:8761–2. — <https://pubmed.ncbi.nlm.nih.gov/23686584/>
54. Lin CF, et al. *Antibodies from dengue patient sera cross-react with endothelial cells and induce damage.* J Med Virol 2003;69:82–90. — <https://pubmed.ncbi.nlm.nih.gov/12436482/>
55. Chuang YC, et al. *Dengue virus nonstructural protein 1-induced antibodies cross-react with human plasminogen and enhance fibrinolysis.* J Immunol 2016;196:1218–26. — <https://pubmed.ncbi.nlm.nih.gov/26712948/>
56. Sun DS, et al. *Antiplatelet autoantibodies elicited by dengue virus non-structural protein 1 cause thrombocytopenia and mortality in mice.* J Thromb Haemost 2007;5:2291–9. — <https://pubmed.ncbi.nlm.nih.gov/17958746/> — the `APLT` compartment.

---

## 6. Endothelium, the glycocalyx and vascular leak

The glycocalyx equations, the `sigma` map and the damage weights come from
this section.

57. **Puerta-Guardo H, Glasner DR, Harris E.** *Dengue virus NS1 disrupts the endothelial glycocalyx, leading to hyperpermeability.* PLoS Pathog 2016;12:e1005738. — <https://pubmed.ncbi.nlm.nih.gov/27416066/> — `wNS1`.
58. Glasner DR, et al. *Dengue virus NS1 cytokine-independent vascular leak is dependent on endothelial glycocalyx components.* PLoS Pathog 2017;13:e1006673. — <https://pubmed.ncbi.nlm.nih.gov/29121099/>
59. Puerta-Guardo H, et al. *Flavivirus NS1 triggers tissue-specific vascular endothelial dysfunction reflecting disease tropism.* Cell Rep 2019;26:1598–1613. — <https://pubmed.ncbi.nlm.nih.gov/30726741/>
60. **Suwarto S, et al.** *Association of endothelial glycocalyx and tight and adherens junctions with severity of plasma leakage in dengue infection.* J Infect Dis 2017;215:992–9. — <https://pubmed.ncbi.nlm.nih.gov/28453844/>
61. Tang THC, et al. *Diagnosis of severe dengue: challenges, needs and opportunities.* J Infect Public Health 2020;13:193–8. — <https://pubmed.ncbi.nlm.nih.gov/31405790/>
62. Lin CF, et al. *Endothelial cell apoptosis induced by antibodies against dengue virus nonstructural protein 1 via production of nitric oxide.* J Immunol 2002;169:657–64. — <https://pubmed.ncbi.nlm.nih.gov/12097367/>
63. Michels M, et al. *Imaging of the sublingual microcirculation in dengue.* Am J Trop Med Hyg 2013;89:139–43. — <https://pubmed.ncbi.nlm.nih.gov/23716404/>
64. Reitsma S, et al. *The endothelial glycocalyx: composition, functions, and visualization.* Pflugers Arch 2007;454:345–59. — <https://pubmed.ncbi.nlm.nih.gov/17256154/>
65. Woodcock TE, Woodcock TM. *Revised Starling equation and the glycocalyx model of transvascular fluid exchange.* Br J Anaesth 2012;108:384–94. — <https://pubmed.ncbi.nlm.nih.gov/22290457/>
66. Levick JR, Michel CC. *Microvascular fluid exchange and the revised Starling principle.* Cardiovasc Res 2010;87:198–210. — <https://pubmed.ncbi.nlm.nih.gov/20200043/>
67. Bhatt P, et al. *Current understanding of the pathogenesis of dengue virus infection.* Curr Microbiol 2021;78:17–32. — <https://pubmed.ncbi.nlm.nih.gov/33231723/>

---

## 7. Mediators: TNF, VEGF, IL-10 and mast cells

68. Bethell DB, et al. *Pathophysiologic and prognostic role of cytokines in dengue hemorrhagic fever.* J Infect Dis 1998;177:778–82. — <https://pubmed.ncbi.nlm.nih.gov/9498463/> — TNF scale.
69. Srikiatkhachorn A, et al. *Virus-induced decline in soluble vascular endothelial growth receptor 2 is associated with plasma leakage in dengue hemorrhagic fever.* J Virol 2007;81:1592–1600. — <https://pubmed.ncbi.nlm.nih.gov/17151115/> — `wVEG`.
70. Malavige GN, et al. *Serum IL-10 as a marker of severe dengue infection.* BMC Infect Dis 2013;13:341. — <https://pubmed.ncbi.nlm.nih.gov/23883139/> — `kIL10`.
71. **Tissera H, et al.** *Chymase level is a predictive biomarker of dengue hemorrhagic fever in pediatric and adult patients.* J Infect Dis 2017;216:1112–21. — <https://pubmed.ncbi.nlm.nih.gov/28968807/> — `wCHY`, `KCHY`.
72. St John AL, et al. *Immune surveillance by mast cells during dengue infection promotes natural killer cell recruitment and viral clearance.* Proc Natl Acad Sci USA 2011;108:9190–5. — <https://pubmed.ncbi.nlm.nih.gov/21576486/>
73. Rathore APS, et al. *Serum chymase levels correlate with severe dengue warning signs and clinical fluid accumulation in hospitalized pediatric patients.* Sci Rep 2020;10:11856. — <https://pubmed.ncbi.nlm.nih.gov/32678240/>
74. Srikiatkhachorn A, Kelley JF. *Endothelial cells in dengue hemorrhagic fever.* Antiviral Res 2014;109:160–70. — <https://pubmed.ncbi.nlm.nih.gov/25025934/>
75. Rathakrishnan A, et al. *Cytokine expression profile of dengue patients at different phases of illness.* PLoS One 2012;7:e52215. — <https://pubmed.ncbi.nlm.nih.gov/23284941/>

---

## 8. Fluid physiology, effusions and the serosal bed

The second Starling bed — the one with a drain two orders of magnitude
smaller than systemic lymph — comes from this section.

76. **Miserocchi G.** *Physiology and pathophysiology of pleural fluid turnover.* Eur Respir J 1997;10:219–25. — <https://pubmed.ncbi.nlm.nih.gov/9032518/> — `DSERMAX` (~0.65 mL/kg/h).
77. Wiig H, Swartz MA. *Interstitial fluid and lymph formation and transport.* Physiol Rev 2012;92:1005–60. — <https://pubmed.ncbi.nlm.nih.gov/22811424/> — lymph recruitment (`kLYMPH`, `JLMAX`).
78. Guyton AC, Hall JE. *Textbook of Medical Physiology*, chapters on capillary exchange and the lymphatic system. — baseline Pc 17.3, Pi −3.0, Πp 27.4, Πi 8.0 mmHg; `Kf` and `Kfs` follow.
79. Landis EM, Pappenheimer JR. *Exchange of substances through the capillary walls.* Handbook of Physiology, 1963. — the colloid-osmotic-pressure polynomial used on TOTAL plasma protein.
80. Srikiatkhachorn A, et al. *Natural history of plasma leakage in dengue hemorrhagic fever: a serial ultrasonographic study.* Pediatr Infect Dis J 2007;26:283–90. — <https://pubmed.ncbi.nlm.nih.gov/17414388/> — effusion volumes and timing.
81. Michels M, et al. *Long-term ultrasonographic follow-up of dengue.* Am J Trop Med Hyg 2013;88:876–81. — <https://pubmed.ncbi.nlm.nih.gov/23509120/>
82. Wills BA, et al. *Coagulation abnormalities in dengue hemorrhagic fever: serial investigations in 167 Vietnamese children.* Clin Infect Dis 2002;35:277–85. — <https://pubmed.ncbi.nlm.nih.gov/12115093/>
83. Gamble J, et al. *Age-related changes in microvascular permeability: a significant factor in the susceptibility of children to shock?* Clin Sci 2000;98:211–6. — <https://pubmed.ncbi.nlm.nih.gov/10657278/>

---

## 9. Fluid therapy trials

84. **Wills BA, et al.** *Comparison of three fluid solutions for resuscitation in dengue shock syndrome.* N Engl J Med 2005;353:877–89. — <https://pubmed.ncbi.nlm.nih.gov/16135832/> — the crystalloid/colloid comparison the model reproduces.
85. Dung NM, et al. *Fluid replacement in dengue shock syndrome: a randomized, double-blind comparison of four intravenous-fluid regimens.* Clin Infect Dis 1999;29:787–94. — <https://pubmed.ncbi.nlm.nih.gov/10589889/>
86. Ngo NT, et al. *Acute management of dengue shock syndrome: a randomized double-blind comparison of 4 intravenous fluid regimens in the first hour.* Clin Infect Dis 2001;32:204–13. — <https://pubmed.ncbi.nlm.nih.gov/11170909/>
87. Nhan NT, et al. *Acute management of dengue shock syndrome: a randomized double-blind comparison of 4 intravenous fluid regimens.* Clin Infect Dis 2001;32:204–13. — <https://pubmed.ncbi.nlm.nih.gov/11170909/>
88. Rosenberger KD, et al. *Early diagnostic indicators of dengue versus other febrile illnesses in Asia and Latin America (IDAMS study).* Lancet Glob Health 2023;11:e361–72. — <https://pubmed.ncbi.nlm.nih.gov/36796983/>
89. Yacoub S, et al. *Cardio-haemodynamic assessment and venous lactate in severe dengue.* PLoS Negl Trop Dis 2017;11:e0005740. — <https://pubmed.ncbi.nlm.nih.gov/28723901/> — lactate and pulse-pressure targets.
90. Kalayanarooj S. *Clinical manifestations and management of dengue/DHF/DSS.* Trop Med Health 2011;39(4 Suppl):83–7. — <https://pubmed.ncbi.nlm.nih.gov/22500140/>

---

## 10. Haematology, thrombocytopenia and bleeding

91. **Lye DC, et al.** *Prophylactic platelet transfusion plus supportive care versus supportive care alone in adults with dengue and thrombocytopenia (AAPT): a multicentre, open-label, randomised, superiority trial.* Lancet 2017;389:1611–8. — <https://pubmed.ncbi.nlm.nih.gov/28283286/> — the transfusion-futility result the model reproduces structurally.
92. Lye DC, et al. *Lack of efficacy of prophylactic platelet transfusion for severe thrombocytopenia in adults with acute uncomplicated dengue infection.* Clin Infect Dis 2009;48:1262–5. — <https://pubmed.ncbi.nlm.nih.gov/19292665/>
93. de Azeredo EL, Monteiro RQ, de-Oliveira Pinto LM. *Thrombocytopenia in dengue: interrelationship between virus and the imbalance between coagulation and fibrinolysis and inflammatory mediators.* Mediators Inflamm 2015;2015:313842. — <https://pubmed.ncbi.nlm.nih.gov/25999666/>
94. Hottz ED, et al. *Platelets in dengue infection.* Drug Discov Today Dis Mech 2011;8:e33–8. — <https://pubmed.ncbi.nlm.nih.gov/22368665/>
95. Noisakran S, et al. *Role of CD61+ cells in thrombocytopenia of dengue patients.* Int J Hematol 2012;96:600–10. — <https://pubmed.ncbi.nlm.nih.gov/23076877/>
96. Wills B, et al. *Hemostatic changes in Vietnamese children with mild dengue correlate with the severity of vascular leakage rather than bleeding.* Am J Trop Med Hyg 2009;81:638–44. — <https://pubmed.ncbi.nlm.nih.gov/19815879/> — the basis for treating bleeding as a product in which the vessel wall dominates.

---

## 11. Liver, organ involvement and severity

97. Trung DT, et al. *Liver involvement associated with dengue infection in adults in Vietnam.* Am J Trop Med Hyg 2010;83:774–80. — <https://pubmed.ncbi.nlm.nih.gov/20889864/> — AST/ALT targets.
98. Samanta J, Sharma V. *Dengue and its effects on liver.* World J Clin Cases 2015;3:125–31. — <https://pubmed.ncbi.nlm.nih.gov/25685758/>
99. Lee LK, et al. *Clinical relevance and discriminatory value of elevated liver aminotransferase levels for dengue severity.* PLoS Negl Trop Dis 2012;6:e1676. — <https://pubmed.ncbi.nlm.nih.gov/22679523/>
100. Póvoa TF, et al. *The pathology of severe dengue in multiple organs of human fatal cases.* PLoS One 2014;9:e83386. — <https://pubmed.ncbi.nlm.nih.gov/24736395/>
101. Yacoub S, et al. *Cardiac function in Vietnamese patients with different dengue severity grades.* Crit Care Med 2012;40:477–83. — <https://pubmed.ncbi.nlm.nih.gov/22020238/> — myocardial depression (`kMYOC`).
102. Oliveira Neto AF, et al. *Kidney injury in dengue.* Rev Inst Med Trop Sao Paulo 2019;61:e11. — <https://pubmed.ncbi.nlm.nih.gov/30916234/>

---

## 12. Antivirals, adjunctive therapy and negative trials

103. **Nguyen NM, et al.** *A randomized, double-blind placebo controlled trial of balapiravir, a polymerase inhibitor, in adult dengue patients.* J Infect Dis 2013;207:1442–50. — <https://pubmed.ncbi.nlm.nih.gov/23453731/> — enrolment at ≤72 h; the timing result the model explains.
104. Low JG, et al. *Efficacy and safety of celgosivir in patients with dengue fever (CELADEN): a phase 1b, randomised, double-blind, placebo-controlled, proof-of-concept trial.* Lancet Infect Dis 2014;14:706–15. — <https://pubmed.ncbi.nlm.nih.gov/24877997/>
105. Tricou V, et al. *A randomized controlled trial of chloroquine for the treatment of dengue in Vietnamese adults.* PLoS Negl Trop Dis 2010;4:e785. — <https://pubmed.ncbi.nlm.nih.gov/20706626/>
106. **Whitehorn J, et al.** *Lovastatin for the treatment of adult patients with dengue: a randomized, double-blind, placebo-controlled trial.* Clin Infect Dis 2016;62:468–76. — <https://pubmed.ncbi.nlm.nih.gov/26565005/>
107. **Tam DTH, et al.** *Effects of short-course oral corticosteroid therapy in early dengue infection in Vietnamese patients: a randomized, placebo-controlled trial.* Clin Infect Dis 2012;55:1216–24. — <https://pubmed.ncbi.nlm.nih.gov/22865871/>
108. Panpanich R, Sornchai P, Kanjanaratanakorn K. *Corticosteroids for treating dengue shock syndrome.* Cochrane Database Syst Rev 2006;(3):CD003488. — <https://pubmed.ncbi.nlm.nih.gov/16856011/>
109. Kaptein SJF, et al. *A pan-serotype dengue virus inhibitor targeting the NS3–NS4B interaction.* Nature 2021;598:504–9. — <https://pubmed.ncbi.nlm.nih.gov/34616043/> — JNJ-1802/mosnodenvir, the antiviral modelled here.
110. Goethals O, et al. *Blocking NS3–NS4B interaction inhibits dengue virus in non-human primates.* Nature 2023;615:678–86. — <https://pubmed.ncbi.nlm.nih.gov/36922586/>
111. Whitehorn J, et al. *Prophylactic platelets in dengue: survey responses highlight lack of an evidence base.* PLoS Negl Trop Dis 2012;6:e1716. — <https://pubmed.ncbi.nlm.nih.gov/22745844/>
112. Zhang H, et al. *Effect of paracetamol on liver injury in dengue: a review.* J Clin Pharm Ther 2020;45:1157–64. — <https://pubmed.ncbi.nlm.nih.gov/32671868/>

---

## 13. Vaccines and the policy consequence of the bell curve

113. **Sridhar S, et al.** *Effect of dengue serostatus on dengue vaccine safety and efficacy.* N Engl J Med 2018;379:327–40. — <https://pubmed.ncbi.nlm.nih.gov/29897841/> — CYD-TDV harms the seronegative.
114. Hadinegoro SR, et al. *Efficacy and long-term safety of a dengue vaccine in regions of endemic disease.* N Engl J Med 2015;373:1195–1206. — <https://pubmed.ncbi.nlm.nih.gov/26214039/>
115. Halstead SB, Russell PK. *Protective and immunological behavior of chimeric yellow fever dengue vaccine.* Vaccine 2016;34:1643–7. — <https://pubmed.ncbi.nlm.nih.gov/26873054/>
116. Biswal S, et al. *Efficacy of a tetravalent dengue vaccine in healthy children and adolescents (TIDES/TAK-003).* N Engl J Med 2019;381:2009–19. — <https://pubmed.ncbi.nlm.nih.gov/31693803/>
117. Rivera L, et al. *Three-year efficacy and safety of Takeda's dengue vaccine candidate (TAK-003).* Clin Infect Dis 2022;75:107–17. — <https://pubmed.ncbi.nlm.nih.gov/34606595/>
118. Kallás EG, et al. *Live, attenuated, tetravalent Butantan-Dengue vaccine in children and adults.* N Engl J Med 2024;390:397–408. — <https://pubmed.ncbi.nlm.nih.gov/38294972/>
119. World Health Organization. *Dengue vaccine: WHO position paper, September 2018.* Wkly Epidemiol Rec 2018;93:457–76. — <https://pubmed.ncbi.nlm.nih.gov/30264474/>
120. Wilder-Smith A, et al. *Deliberations of the Strategic Advisory Group of Experts on Immunization on the use of CYD-TDV dengue vaccine.* Lancet Infect Dis 2019;19:e31–8. — <https://pubmed.ncbi.nlm.nih.gov/30195995/>
121. Flasche S, et al. *The long-term safety, public health impact, and cost-effectiveness of routine vaccination with a recombinant, live-attenuated dengue vaccine (Dengvaxia): a model comparison study.* PLoS Med 2016;13:e1002181. — <https://pubmed.ncbi.nlm.nih.gov/27898668/>

---

## 14. Modelling method, QSP and tooling

122. Baron KT, et al. *mrgsolve: Simulate from ODE-based population PK/PD and systems pharmacology models.* — <https://mrgsolve.org/>
123. Virtanen P, et al. *SciPy 1.0: fundamental algorithms for scientific computing in Python.* Nat Methods 2020;17:261–72. — <https://pubmed.ncbi.nlm.nih.gov/32015543/>
124. Hindmarsh AC, Petzold LR. *LSODA, ordinary differential equation solver for stiff or non-stiff system.* — the integrator used for the reference implementation.
125. Sorger PK, et al. *Quantitative and systems pharmacology in the post-genomic era.* NIH White Paper, 2011. — <https://www.nigms.nih.gov/training/documents/systemspharmawpsorger2011.pdf>
126. Gadkar K, et al. *A six-stage workflow for robust application of systems pharmacology.* CPT Pharmacometrics Syst Pharmacol 2016;5:235–49. — <https://pubmed.ncbi.nlm.nih.gov/27299936/>

---

## Notes on what is *not* in the model

* **Spatial heterogeneity of the leak.** The model has one systemic bed and
  one serosal bed. Real leakage is regional and the pleural and peritoneal
  cavities differ in compliance; here they are lumped into one compartment
  with one compliance (`CSER`).
* **The incubation period.** `t = 0` is fever onset. Every scenario starts
  from the same viraemia so the antibody state, not a different starting
  point, is the only thing that differs. The epidemiological titre sweep is
  the one analysis run from inoculation instead, precisely because there a
  protected host must be allowed not to become infected.
* **Serotype identity.** DENV-1 to DENV-4 differ in replicative fitness and in
  the cross-reactivity structure of the antibody response. The model carries a
  single scalar titre and a single `FCROSS`; it cannot say which serotype
  sequence is worst.
* **Age.** Children leak more than adults at the same mediator load
  (ref. 83). The model is parameterised for a 70 kg adult.
* **Coagulation in detail.** Fibrinogen is a single state; there is no
  explicit thrombin, no protein C axis and no fibrinolysis. Bleeding is a
  product of four requirements rather than a coagulation cascade.
