# IgE-mediated food allergy and anaphylaxis — references

Literature underpinning `fana_qsp_model.dot`, `fana_mrgsolve_model.R` and
`fana_shiny_app.R`. Sections follow the clusters of the mechanistic map.
Every entry carries a PubMed link. Where a reference supplies a number that
is actually used as a parameter or as a calibration target, that is stated.

---

## 1. Epidemiology, natural history and the burden

1. Gupta RS, Warren CM, Smith BM, et al. **Prevalence and severity of food allergies among US adults.** *JAMA Netw Open.* 2019;2(1):e185630. — https://pubmed.ncbi.nlm.nih.gov/30646188/
2. Sicherer SH, Sampson HA. **Food allergy: a review and update on epidemiology, pathogenesis, diagnosis, prevention, and management.** *J Allergy Clin Immunol.* 2018;141(1):41-58. — https://pubmed.ncbi.nlm.nih.gov/29157945/
3. Savage J, Sicherer S, Wood R. **The natural history of food allergy.** *J Allergy Clin Immunol Pract.* 2016;4(2):196-203. — https://pubmed.ncbi.nlm.nih.gov/26968958/
4. Skolnick HS, Conover-Walker MK, Koerner CB, et al. **The natural history of peanut allergy.** *J Allergy Clin Immunol.* 2001;107(2):367-374. — https://pubmed.ncbi.nlm.nih.gov/11174206/ — *source of the ~20% resolution rate used for `NATRES = 0` in peanut.*
5. Wood RA. **The natural history of food allergy.** *Pediatrics.* 2003;111(6 Pt 3):1631-1637. — https://pubmed.ncbi.nlm.nih.gov/12777603/ — *milk/egg resolution, the basis of `NATRES ≈ 0.9`.*
6. Turner PJ, Baumert JL, Beyer K, et al. **Can we identify patients at risk of life-threatening allergic reactions to food?** *Allergy.* 2016;71(9):1241-1255. — https://pubmed.ncbi.nlm.nih.gov/27138061/
7. Turner PJ, Jerschow E, Umasunthar T, et al. **Fatal anaphylaxis: mortality rate and risk factors.** *J Allergy Clin Immunol Pract.* 2017;5(5):1169-1178. — https://pubmed.ncbi.nlm.nih.gov/28888247/ — *the ~0.03-0.3 per million person-years figure in cluster 19.*
8. Pumphrey RS. **Lessons for management of anaphylaxis from a study of fatal reactions.** *Clin Exp Allergy.* 2000;30(8):1144-1150. — https://pubmed.ncbi.nlm.nih.gov/10931122/ — *median time to arrest by trigger (food ~30 min); the posture deaths in cluster 8.*

## 2. Sensitisation, the dual-allergen-exposure hypothesis and prevention

9. Lack G. **Epidemiologic risks for food allergy.** *J Allergy Clin Immunol.* 2008;121(6):1331-1336. — https://pubmed.ncbi.nlm.nih.gov/18539191/
10. Du Toit G, Roberts G, Sayre PH, et al. **Randomized trial of peanut consumption in infants at risk for peanut allergy (LEAP).** *N Engl J Med.* 2015;372(9):803-813. — https://pubmed.ncbi.nlm.nih.gov/25705822/ — *13.7% → 1.9% at 5 years; the tolerance node of cluster 1.*
11. Du Toit G, Sayre PH, Roberts G, et al. **Effect of avoidance on peanut allergy after early peanut consumption (LEAP-On).** *N Engl J Med.* 2016;374(15):1435-1443. — https://pubmed.ncbi.nlm.nih.gov/26942922/
12. Brough HA, Liu AH, Sicherer S, et al. **Atopic dermatitis increases the effect of exposure to peanut antigen in dust on peanut sensitization and likely food allergy.** *J Allergy Clin Immunol.* 2015;135(1):164-170. — https://pubmed.ncbi.nlm.nih.gov/25457149/
13. Palmer CN, Irvine AD, Terron-Kwiatkowski A, et al. **Common loss-of-function variants of the epidermal barrier protein filaggrin are a major predisposing factor for atopic dermatitis.** *Nat Genet.* 2006;38(4):441-446. — https://pubmed.ncbi.nlm.nih.gov/16550169/
14. Tordesillas L, Berin MC, Sampson HA. **Immunology of food allergy.** *Immunity.* 2017;47(1):32-50. — https://pubmed.ncbi.nlm.nih.gov/28723552/
15. Stefka AT, Feehley T, Tripathi P, et al. **Commensal bacteria protect against food allergen sensitization.** *Proc Natl Acad Sci USA.* 2014;111(36):13145-13150. — https://pubmed.ncbi.nlm.nih.gov/25157157/

## 3. The allergen: molecular identity, digestion resistance, valency

16. Koppelman SJ, Hefle SL, Taylor SL, de Jong GA. **Digestion of peanut allergens Ara h 1, Ara h 2, Ara h 3 and Ara h 6.** *Mol Nutr Food Res.* 2010;54(12):1711-1721. — https://pubmed.ncbi.nlm.nih.gov/20603830/ — *the pepsin resistance of the 2S albumins that motivates `KDEG_GUT` and `FARAH2`.*
17. Blanc F, Adel-Patient K, Drumare MF, et al. **Capacity of purified peanut allergens to induce degranulation in a functional in vitro assay: Ara h 2 and Ara h 6 are the most efficient elicitors.** *Clin Exp Allergy.* 2009;39(8):1277-1285. — https://pubmed.ncbi.nlm.nih.gov/19538351/
18. Klemans RJ, van Os-Medendorp H, Blankestijn M, et al. **Diagnostic accuracy of specific IgE to components in diagnosing peanut allergy: a systematic review.** *Clin Exp Allergy.* 2015;45(4):720-730. — https://pubmed.ncbi.nlm.nih.gov/25226880/
19. Suárez-Fariñas M, Suprun M, Chang HL, et al. **Predicting development of sustained unresponsiveness to milk oral immunotherapy using epitope-specific antibody binding profiles.** *J Allergy Clin Immunol.* 2019;143(3):1038-1046. — https://pubmed.ncbi.nlm.nih.gov/30528770/
20. Nowak-Węgrzyn A, Bloom KA, Sicherer SH, et al. **Tolerance to extensively heated milk in children with cow's milk allergy.** *J Allergy Clin Immunol.* 2008;122(2):342-347. — https://pubmed.ncbi.nlm.nih.gov/18620748/ — *the baked-tolerance node of cluster 18.*
21. Untersmayr E, Jensen-Jarolim E. **The role of protein digestibility and antacids on food allergy outcomes.** *J Allergy Clin Immunol.* 2008;121(6):1301-1308. — https://pubmed.ncbi.nlm.nih.gov/18539189/ — *the PPI cofactor `FPPI`.*
22. Grimshaw KE, King RM, Nordlee JA, et al. **Presentation of allergen in different food preparations affects the nature of the allergic reaction — a case series.** *Clin Exp Allergy.* 2003;33(11):1581-1585. — https://pubmed.ncbi.nlm.nih.gov/14616872/ — *the fat-matrix effect (`FATMATRIX`).*

## 4. FcεRI, the surface, and the receptor as a state variable

23. Kinet JP. **The high-affinity IgE receptor (FcεRI): from physiology to pathology.** *Annu Rev Immunol.* 1999;17:931-972. — https://pubmed.ncbi.nlm.nih.gov/10358778/ — *KD ~1e-10 M, the basis of `KD_FCERI`.*
24. MacGlashan DW Jr, Bochner BS, Adelman DC, et al. **Down-regulation of FcεRI expression on human basophils during in vivo treatment of atopic patients with anti-IgE antibody.** *J Immunol.* 1997;158(3):1438-1445. — https://pubmed.ncbi.nlm.nih.gov/9013989/ — *the slow arm (S2); the ~10-30x basophil fall that the model deliberately does NOT apply to the tissue mast cell.*
25. Beck LA, Marcotte GV, MacGlashan D, et al. **Omalizumab-induced reductions in mast cell FcεRI expression and function.** *J Allergy Clin Immunol.* 2004;114(3):527-530. — https://pubmed.ncbi.nlm.nih.gov/15356552/ — *skin mast cell FcεRI falls less and later than basophil; the source of `RHO_FLOOR = 0.18`.*
26. MacGlashan D Jr. **IgE receptor and signal transduction in mast cells and basophils.** *Curr Opin Immunol.* 2008;20(6):717-723. — https://pubmed.ncbi.nlm.nih.gov/18822373/
27. Borkowski TA, Jouvin MH, Lin SY, Kinet JP. **Minimal requirements for IgE-mediated regulation of surface FcεRI.** *J Immunol.* 2001;167(3):1290-1296. — https://pubmed.ncbi.nlm.nih.gov/11466345/ — *free IgE stabilises the receptor: the positive feedback in cluster 4.*
28. Kepley CL, Youssef L, Andrews RP, et al. **Syk deficiency in nonreleaser basophils.** *J Allergy Clin Immunol.* 1999;104(2 Pt 1):279-284. — https://pubmed.ncbi.nlm.nih.gov/10452746/ — *the `RELEASE` parameter.*

## 5. Cross-linking: the quadratic law, aggregation requirements, the hook effect

29. Dembo M, Goldstein B. **Theory of equilibrium binding of symmetric bivalent haptens to cell surface antibody: application to histamine release from basophils.** *J Immunol.* 1978;121(1):345-353. — https://pubmed.ncbi.nlm.nih.gov/78964/ — *the foundational bell-shaped cross-linking theory; the origin of the square law used in cluster 5.*
30. Goldstein B, Perelson AS. **Equilibrium theory for the clustering of bivalent cell surface receptors by trivalent ligands.** *Biophys J.* 1984;45(6):1109-1123. — https://pubmed.ncbi.nlm.nih.gov/6743742/
31. Faeder JR, Hlavacek WS, Reischl I, et al. **Investigation of early events in FcεRI-mediated signaling using a detailed mathematical model.** *J Immunol.* 2003;170(7):3769-3781. — https://pubmed.ncbi.nlm.nih.gov/12646643/ — *rule-based FcεRI aggregation model; the quantitative ancestor of the engine in cluster 5.*
32. Posner RG, Geng D, Haymore S, et al. **Trivalent antigens for degranulation of mast cells.** *Org Lett.* 2007;9(18):3551-3554. — https://pubmed.ncbi.nlm.nih.gov/17691792/
33. Christensen LH, Holm J, Lund G, et al. **Several distinct properties of the IgE repertoire determine effector cell degranulation in response to allergen challenge.** *J Allergy Clin Immunol.* 2008;122(2):298-304. — https://pubmed.ncbi.nlm.nih.gov/18572230/ — *epitope diversity and the specific-IgE FRACTION, not the titre, drive degranulation. The single best experimental support for S1.*
34. Handlogten MW, Kiziltepe T, Serezani AP, et al. **Inhibition of weak-affinity epitope-IgE interactions prevents mast cell degranulation.** *Nat Chem Biol.* 2013;9(12):789-795. — https://pubmed.ncbi.nlm.nih.gov/24141197/

## 6. Signalling and its brakes

35. Gilfillan AM, Tkaczyk C. **Integrated signalling pathways for mast-cell activation.** *Nat Rev Immunol.* 2006;6(3):218-230. — https://pubmed.ncbi.nlm.nih.gov/16470226/
36. Huber M, Helgason CD, Damen JE, et al. **The src homology 2-containing inositol phosphatase (SHIP) is the gatekeeper of mast cell degranulation.** *Proc Natl Acad Sci USA.* 1998;95(19):11330-11335. — https://pubmed.ncbi.nlm.nih.gov/9736736/
37. Daëron M, Malbec O, Latour S, et al. **Regulation of high-affinity IgE receptor-mediated mast cell activation by murine low-affinity IgG receptors.** *J Clin Invest.* 1995;95(2):577-585. — https://pubmed.ncbi.nlm.nih.gov/7532187/ — *FcγRIIb co-aggregation; the second IgG4 mechanism in cluster 14.*
38. Dispenza MC, Krier-Burris RA, Chhiba KD, et al. **Bruton's tyrosine kinase inhibition effectively protects against human IgE-mediated anaphylaxis.** *J Clin Invest.* 2020;130(9):4759-4770. — https://pubmed.ncbi.nlm.nih.gov/32484802/ — *the `BTKI` node: two days of dosing abolishes reactivity without touching IgE.*
39. MacGlashan DW Jr. **Basophil activation testing.** *J Allergy Clin Immunol.* 2013;132(4):777-787. — https://pubmed.ncbi.nlm.nih.gov/24001573/

## 7. Mediators

40. Vadas P, Gold M, Perelman B, et al. **Platelet-activating factor, PAF acetylhydrolase, and severe anaphylaxis.** *N Engl J Med.* 2008;358(1):28-35. — https://pubmed.ncbi.nlm.nih.gov/18172172/ — *the `PAFAH` parameter and the dominant weight `WP_LEAK`.*
41. Vadas P, Perelman B, Liss G. **Platelet-activating factor, histamine, and tryptase levels in human anaphylaxis.** *J Allergy Clin Immunol.* 2013;131(1):144-149. — https://pubmed.ncbi.nlm.nih.gov/23040367/
42. Schwartz LB, Metcalfe DD, Miller JS, et al. **Tryptase levels as an indicator of mast-cell activation in systemic anaphylaxis and mastocytosis.** *N Engl J Med.* 1987;316(26):1622-1626. — https://pubmed.ncbi.nlm.nih.gov/3295549/
43. Valent P, Akin C, Arock M, et al. **Definitions, criteria and global classification of mast cell disorders with special reference to mast cell activation syndromes: a consensus proposal.** *Int Arch Allergy Immunol.* 2012;157(3):215-225. — https://pubmed.ncbi.nlm.nih.gov/22041891/ — *the (1.2 × baseline + 2) tryptase rise criterion used in diagnostic D12.*
44. Lieberman P, Garvey LH. **Mast cells and anaphylaxis.** *Curr Allergy Asthma Rep.* 2016;16(3):20. — https://pubmed.ncbi.nlm.nih.gov/26857018/
45. Lyons JJ, Yu X, Hughes JD, et al. **Elevated basal serum tryptase identifies a multisystem disorder associated with increased TPSAB1 copy number.** *Nat Genet.* 2016;48(12):1564-1569. — https://pubmed.ncbi.nlm.nih.gov/27749843/ — *hereditary alpha-tryptasaemia; the `MCBURDEN` severity modifier.*

## 8. Cardiovascular physiology of anaphylaxis

46. Fisher MM. **Clinical observations on the pathophysiology and treatment of anaphylactic cardiovascular collapse.** *Anaesth Intensive Care.* 1986;14(1):17-21. — https://pubmed.ncbi.nlm.nih.gov/2869715/ — *up to 35% of plasma volume shifts in ~10 minutes: the calibration target for `KLEAK`, `FRAC50` and the whole of S4.*
47. Brown SG. **Cardiovascular aspects of anaphylaxis: implications for treatment and diagnosis.** *Curr Opin Allergy Clin Immunol.* 2005;5(4):359-364. — https://pubmed.ncbi.nlm.nih.gov/15985820/
48. Brown SG. **Clinical features and severity grading of anaphylaxis.** *J Allergy Clin Immunol.* 2004;114(2):371-376. — https://pubmed.ncbi.nlm.nih.gov/15316518/ — *the severity grading whose ordinal structure cluster 10 criticises.*
49. Kounis NG. **Kounis syndrome: an update on epidemiology, pathogenesis, diagnosis and therapeutic management.** *Clin Chem Lab Med.* 2016;54(10):1545-1559. — https://pubmed.ncbi.nlm.nih.gov/26966931/
50. Pumphrey RS. **Fatal posture in anaphylactic shock.** *J Allergy Clin Immunol.* 2003;112(2):451-452. — https://pubmed.ncbi.nlm.nih.gov/12897756/ — *the "empty ventricle" deaths; the `SUPINE` parameter.*

## 9. Cofactors and the variability of the threshold

51. Niggemann B, Beyer K. **Factors augmenting allergic reactions.** *Allergy.* 2014;69(12):1582-1587. — https://pubmed.ncbi.nlm.nih.gov/25306896/
52. Wölbing F, Fischer J, Köberle M, et al. **About the role and underlying mechanisms of cofactors in anaphylaxis.** *Allergy.* 2013;68(9):1085-1092. — https://pubmed.ncbi.nlm.nih.gov/23909934/ — *the mechanistic basis for treating cofactors as DELIVERY multipliers (S5).*
53. Christensen MJ, Eller E, Mortz CG, et al. **Exercise lowers threshold and increases severity, but wheat-dependent, exercise-induced anaphylaxis can be elicited at rest.** *J Allergy Clin Immunol Pract.* 2018;6(2):514-520. — https://pubmed.ncbi.nlm.nih.gov/28734860/ — *`FEX`.*
54. Dua S, Ruiz-Garcia M, Bond S, et al. **Effect of sleep deprivation and exercise on reaction threshold in adults with peanut allergy: a randomized controlled study.** *J Allergy Clin Immunol.* 2019;144(6):1584-1594. — https://pubmed.ncbi.nlm.nih.gov/31319102/ — *the definitive threshold-shift experiment; the calibration target for the cofactor block.*
55. Hourihane JO, Kilburn SA, Nordlee JA, et al. **An evaluation of the sensitivity of subjects with peanut allergy to very low doses of peanut protein: a randomized, double-blind, placebo-controlled food challenge study.** *J Allergy Clin Immunol.* 1997;100(5):596-600. — https://pubmed.ncbi.nlm.nih.gov/9389287/
56. Taylor SL, Baumert JL, Kruizinga AG, et al. **Establishment of reference doses for residues of allergenic foods: report of the VITAL Expert Panel.** *Food Chem Toxicol.* 2014;63:9-17. — https://pubmed.ncbi.nlm.nih.gov/24184597/ — *the ED01/ED05 population thresholds in cluster 12.*
57. Remington BC, Westerhout J, Meima MY, et al. **Updated population minimal eliciting dose distributions for use in risk assessment of 14 priority food allergens.** *Food Chem Toxicol.* 2020;139:111259. — https://pubmed.ncbi.nlm.nih.gov/32179163/ — *the log-normal population ED distribution used by `vpop()`.*
58. Blumchen K, Beder A, Beschorner J, et al. **Modified oral food challenge used with sensitization biomarkers provides more real-life clinical thresholds for peanut allergy.** *J Allergy Clin Immunol.* 2014;134(2):390-398. — https://pubmed.ncbi.nlm.nih.gov/24831437/

## 10. Diagnostics

59. Sampson HA, Gerth van Wijk R, Bindslev-Jensen C, et al. **Standardizing double-blind, placebo-controlled oral food challenges: American Academy of Allergy, Asthma & Immunology–European Academy of Allergy and Clinical Immunology PRACTALL consensus report.** *J Allergy Clin Immunol.* 2012;130(6):1260-1274. — https://pubmed.ncbi.nlm.nih.gov/23195525/ — *the escalation schedule implemented in `ev_challenge_practall()`.*
60. Santos AF, Douiri A, Bécares N, et al. **Basophil activation test discriminates between allergy and tolerance in peanut-sensitized children.** *J Allergy Clin Immunol.* 2014;134(3):645-652. — https://pubmed.ncbi.nlm.nih.gov/25065721/ — *BAT measures the engine directly; the best single correlate of threshold.*
61. Santos AF, Du Toit G, Douiri A, et al. **Distinct parameters of the basophil activation test reflect the severity and threshold of allergic reactions to peanut.** *J Allergy Clin Immunol.* 2015;135(1):179-186. — https://pubmed.ncbi.nlm.nih.gov/25567045/
62. Bahri R, Custovic A, Korosec P, et al. **Mast cell activation test in the diagnosis of allergic disease and anaphylaxis.** *J Allergy Clin Immunol.* 2018;142(2):485-496. — https://pubmed.ncbi.nlm.nih.gov/29518421/
63. Sindher SB, Long AJ, Purington N, et al. **Analysis of a large standardized food challenge data set to determine predictors of positive outcome across multiple allergens.** *Front Immunol.* 2018;9:2689. — https://pubmed.ncbi.nlm.nih.gov/30524436/ — *the sIgE/total-IgE ratio outperforming sIgE alone: the clinical observation that S1 explains.*
64. Gupta RS, Lau CH, Hamilton RG, et al. **Predicting outcomes of oral food challenges by using the allergen-specific IgE-total IgE ratio.** *J Allergy Clin Immunol Pract.* 2014;2(3):300-305. — https://pubmed.ncbi.nlm.nih.gov/24811021/ — *the direct empirical support for the ratio f; the target of diagnostic D2.*

## 11. Anti-IgE pharmacology

65. Lowe PJ, Tannenbaum S, Gautier A, Jimenez P. **Relationship between omalizumab pharmacokinetics, IgE pharmacodynamics and symptoms in patients with severe persistent allergic (IgE-mediated) asthma.** *Br J Clin Pharmacol.* 2009;68(1):61-76. — https://pubmed.ncbi.nlm.nih.gov/19660004/ — *the source of `CL_OMA`, `VC_OMA`, `KA_OMA` and the binding structure of the TMDD block.*
66. Hayashi N, Tsukamoto Y, Sallas WM, Lowe PJ. **A mechanism-based binding model for the population pharmacokinetics and pharmacodynamics of omalizumab.** *Br J Clin Pharmacol.* 2007;63(5):548-561. — https://pubmed.ncbi.nlm.nih.gov/17096680/ — *the 1:1 quasi-equilibrium approximation used in `bind11()` and the complex-clearance ratio behind `KEL_CPX`.*
67. Wood RA, Togias A, Sicherer SH, et al. **Omalizumab for the treatment of multiple food allergies (OUtMATCH).** *N Engl J Med.* 2024;390(10):889-899. — https://pubmed.ncbi.nlm.nih.gov/38407394/ — *67% vs 7% tolerating ≥600 mg peanut protein at week 16-20; the primary validation target of diagnostics D4 and D15.*
68. Sampson HA, Leung DY, Burks AW, et al. **A phase II, randomized, double-blind, parallel-group, placebo-controlled oral food challenge trial of Xolair (omalizumab) in peanut allergy.** *J Allergy Clin Immunol.* 2011;127(5):1309-1310. — https://pubmed.ncbi.nlm.nih.gov/21397314/
69. Savage JH, Courneya JP, Sterba PM, et al. **Kinetics of mast cell, basophil, and oral food challenge responses in omalizumab-treated adults with peanut allergy.** *J Allergy Clin Immunol.* 2012;130(5):1123-1129. — https://pubmed.ncbi.nlm.nih.gov/23021878/ — *the time course that motivates separating the fast and slow arms (S2).*
70. Gasser P, Tarchevskaya SS, Guntern P, et al. **The mechanistic and functional profile of the therapeutic anti-IgE antibody ligelizumab differs from omalizumab.** *Nat Commun.* 2020;11(1):165. — https://pubmed.ncbi.nlm.nih.gov/31913280/ — *higher affinity, non-proportional clinical gain; consistent with saturation of the occupancy term.*
71. Corren J, Casale TB, Lanier B, et al. **Safety and tolerability of omalizumab.** *Clin Exp Allergy.* 2009;39(6):788-797. — https://pubmed.ncbi.nlm.nih.gov/19302249/
72. Fiocchi A, Artesani MC, Riccardi C, et al. **Impact of omalizumab on food allergy in patients treated for asthma: a real-life study.** *J Allergy Clin Immunol Pract.* 2019;7(6):1901-1909. — https://pubmed.ncbi.nlm.nih.gov/30797776/

## 12. Allergen immunotherapy

73. PALISADE Group of Clinical Investigators; Vickery BP, Vereda A, Casale TB, et al. **AR101 oral immunotherapy for peanut allergy.** *N Engl J Med.* 2018;379(21):1991-2001. — https://pubmed.ncbi.nlm.nih.gov/30449234/ — *67.2% vs 4.0% tolerating a single 600 mg dose; the OIT calibration target.*
74. Jones SM, Kim EH, Nadeau KC, et al. **Efficacy and safety of oral immunotherapy in children aged 1-3 years with peanut allergy (IMPACT): a randomised placebo-controlled study.** *Lancet.* 2022;399(10322):359-371. — https://pubmed.ncbi.nlm.nih.gov/35065784/ — *remission 20% vs 7%, strongly age-dependent: the source of `AGE50_TOL` and diagnostic D14.*
75. Chinthrajah RS, Purington N, Andorf S, et al. **Sustained outcomes in oral immunotherapy for peanut allergy (POISED study).** *Lancet.* 2019;394(10207):1437-1449. — https://pubmed.ncbi.nlm.nih.gov/31522849/ — *loss of desensitisation on discontinuation; the two decay rates in D13.*
76. Vickery BP, Berglund JP, Burk CM, et al. **Early oral immunotherapy in peanut-allergic preschool children is safe and highly effective.** *J Allergy Clin Immunol.* 2017;139(1):173-181. — https://pubmed.ncbi.nlm.nih.gov/27522159/
77. Burks AW, Jones SM, Wood RA, et al. **Oral immunotherapy for treatment of egg allergy in children.** *N Engl J Med.* 2012;367(3):233-243. — https://pubmed.ncbi.nlm.nih.gov/22808958/
78. Santos AF, James LK, Bahnson HT, et al. **IgG4 inhibits peanut-induced basophil and mast cell activation in peanut-tolerant children sensitized to peanut major allergens.** *J Allergy Clin Immunol.* 2015;135(5):1249-1256. — https://pubmed.ncbi.nlm.nih.gov/25670011/ — *the direct demonstration of IgG4 interception; the mechanism behind `KI_G4` (S3).*
79. Shamji MH, Valenta R, Jardetzky T, et al. **The role of allergen-specific IgE, IgG and IgA in allergic disease.** *Allergy.* 2021;76(12):3627-3641. — https://pubmed.ncbi.nlm.nih.gov/33999439/
80. Kim EH, Bird JA, Kulis M, et al. **Sublingual immunotherapy for peanut allergy: clinical and immunologic evidence of desensitization.** *J Allergy Clin Immunol.* 2011;127(3):640-646. — https://pubmed.ncbi.nlm.nih.gov/21281959/
81. Fleischer DM, Greenhawt M, Sussman G, et al. **Effect of epicutaneous immunotherapy vs placebo on reaction to peanut protein ingestion among children with peanut allergy (PEPITES).** *JAMA.* 2019;321(10):946-955. — https://pubmed.ncbi.nlm.nih.gov/30794314/
82. Greenhawt M, Sindher SB, Wang J, et al. **Phase 3 trial of epicutaneous immunotherapy in toddlers with peanut allergy (EPITOPE).** *N Engl J Med.* 2023;388(19):1755-1766. — https://pubmed.ncbi.nlm.nih.gov/37163622/
83. Chu DK, Wood RA, French S, et al. **Oral immunotherapy for peanut allergy (PACE): a systematic review and meta-analysis of efficacy and safety.** *Lancet.* 2019;393(10187):2222-2232. — https://pubmed.ncbi.nlm.nih.gov/31030987/ — *the finding that OIT INCREASES anaphylaxis and adrenaline use; the `OITAE` node.*
84. Lucendo AJ, Arias Á, Tenias JM. **Relation between eosinophilic esophagitis and oral immunotherapy for food allergy: a systematic review with meta-analysis.** *Ann Allergy Asthma Immunol.* 2014;113(6):624-629. — https://pubmed.ncbi.nlm.nih.gov/25216976/
85. Akdis CA, Akdis M. **Mechanisms of allergen-specific immunotherapy and immune tolerance to allergens.** *World Allergy Organ J.* 2015;8(1):17. — https://pubmed.ncbi.nlm.nih.gov/26023323/
86. Tang ML, Ponsonby AL, Orsini F, et al. **Administration of a probiotic with peanut oral immunotherapy: a randomized trial.** *J Allergy Clin Immunol.* 2015;135(3):737-744. — https://pubmed.ncbi.nlm.nih.gov/25592987/

## 13. Type 2 blockade and other pharmacology

87. Rial MJ, Barroso B, Sastre J. **Dupilumab for treatment of food allergy.** *J Allergy Clin Immunol Pract.* 2019;7(2):673-674. — https://pubmed.ncbi.nlm.nih.gov/30075339/
88. Spekhorst LS, de Graaf M, van der Rijst LP, et al. **Dupilumab has a profound effect on specific IgE levels of several food allergens in atopic dermatitis patients.** *Allergy.* 2023;78(3):875-878. — https://pubmed.ncbi.nlm.nih.gov/36239235/ — *sIgE and total IgE fall TOGETHER: the observation that makes the model's dupilumab prediction (D6) testable.*
89. Wood RA, Chinthrajah RS, Eggel A, et al. **The rationale for development of ligelizumab in food allergy.** *World Allergy Organ J.* 2022;15(3):100638. — https://pubmed.ncbi.nlm.nih.gov/35341022/
90. Dispenza MC, Bochner BS, MacGlashan DW Jr. **Targeting the FcεRI pathway as a potential strategy to prevent food-induced anaphylaxis.** *Front Immunol.* 2020;11:614402. — https://pubmed.ncbi.nlm.nih.gov/33391286/
91. Maurer M, Berger W, Giménez-Arnau A, et al. **Remibrutinib, a novel BTK inhibitor, in chronic spontaneous urticaria.** *J Allergy Clin Immunol.* 2022;150(6):1498-1506. — https://pubmed.ncbi.nlm.nih.gov/36096203/
92. Terhorst-Molawi D, Hawro T, Grekowitz E, et al. **Anti-KIT antibody barzolvolimab reduces skin mast cells and disease activity in chronic inducible urticaria.** *Allergy.* 2023;78(5):1269-1279. — https://pubmed.ncbi.nlm.nih.gov/36385701/

## 14. Acute management: adrenaline and the rest

93. Simons FE, Gu X, Simons KJ. **Epinephrine absorption in adults: intramuscular versus subcutaneous injection.** *J Allergy Clin Immunol.* 2001;108(5):871-873. — https://pubmed.ncbi.nlm.nih.gov/11692118/ — *the IM-thigh PK reproduced by `KA_EPI_IM`, `KEL_EPI`, `V_EPI`.*
94. Simons FE, Roberts JR, Gu X, Simons KJ. **Epinephrine absorption in children with a history of anaphylaxis.** *J Allergy Clin Immunol.* 1998;101(1 Pt 1):33-37. — https://pubmed.ncbi.nlm.nih.gov/9449498/
95. Dworaczyk D, Hunt A, Casale TB, et al. **Pharmacokinetics and pharmacodynamics of intranasal epinephrine (neffy) compared with intramuscular epinephrine.** *J Allergy Clin Immunol Pract.* 2023;11(9):2760-2768. — https://pubmed.ncbi.nlm.nih.gov/37301433/ — *the `EPIIN` route.*
96. Cardona V, Ansotegui IJ, Ebisawa M, et al. **World Allergy Organization anaphylaxis guidance 2020.** *World Allergy Organ J.* 2020;13(10):100472. — https://pubmed.ncbi.nlm.nih.gov/33204386/
97. Muraro A, Worm M, Alviani C, et al. **EAACI guidelines: anaphylaxis (2021 update).** *Allergy.* 2022;77(2):357-377. — https://pubmed.ncbi.nlm.nih.gov/34343358/
98. Shaker MS, Wallace DV, Golden DBK, et al. **Anaphylaxis — a 2020 practice parameter update, systematic review, and Grading of Recommendations, Assessment, Development and Evaluation (GRADE) analysis.** *J Allergy Clin Immunol.* 2020;145(4):1082-1123. — https://pubmed.ncbi.nlm.nih.gov/32001253/ — *the weak evidence base for corticosteroids and antihistamines; the dashed edges in clusters 11 and 16.*
99. Liyanage CK, Galappatthy P, Seneviratne SL. **Corticosteroids in management of anaphylaxis: a systematic review of evidence.** *Eur Ann Allergy Clin Immunol.* 2017;49(5):196-207. — https://pubmed.ncbi.nlm.nih.gov/28884986/
100. Lee S, Bellolio MF, Hess EP, et al. **Time of onset and predictors of biphasic anaphylactic reactions: a systematic review and meta-analysis.** *J Allergy Clin Immunol Pract.* 2015;3(3):408-416. — https://pubmed.ncbi.nlm.nih.gov/25680923/ — *the ~4-5% biphasic rate and its timing; the late-phase block.*
101. Fromer L. **Prevention of anaphylaxis: the role of the epinephrine auto-injector.** *Am J Med.* 2016;129(12):1244-1250. — https://pubmed.ncbi.nlm.nih.gov/27555094/
102. Kawano T, Scheuermeyer FX, Stenstrom R, et al. **Epinephrine use in older patients with anaphylaxis: clinical outcomes and cardiovascular complications.** *Resuscitation.* 2017;112:53-58. — https://pubmed.ncbi.nlm.nih.gov/28069483/
103. Anagnostou K, Turner PJ. **Myths, facts and controversies in the diagnosis and management of anaphylaxis.** *Arch Dis Child.* 2019;104(1):83-90. — https://pubmed.ncbi.nlm.nih.gov/29909382/

## 15. Quality of life, risk perception and what the threshold is worth

104. Cummings AJ, Knibb RC, King RM, Lucas JS. **The psychosocial impact of food allergy and food hypersensitivity in children, adolescents and their families: a review.** *Allergy.* 2010;65(8):933-945. — https://pubmed.ncbi.nlm.nih.gov/20180792/
105. DunnGalvin A, Blumchen K, Timmermans F, et al. **APPEAL-1: a pan-European survey of patient/caregiver perceptions of peanut allergy management.** *Allergy.* 2020;75(11):2899-2908. — https://pubmed.ncbi.nlm.nih.gov/32441772/
106. Greenhawt M, Shaker M, Winders T, et al. **Development and acceptability of a shared decision-making tool for commercial peanut allergy therapies.** *Ann Allergy Asthma Immunol.* 2020;125(1):90-96. — https://pubmed.ncbi.nlm.nih.gov/32247741/

## 16. QSP and modelling methodology

107. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn.* 2015;42:S84-S85. (mrgsolve) — https://mrgsolve.org
108. Mager DE, Jusko WJ. **General pharmacokinetic model for drugs exhibiting target-mediated drug disposition.** *J Pharmacokinet Pharmacodyn.* 2001;28(6):507-532. — https://pubmed.ncbi.nlm.nih.gov/11999290/ — *the TMDD framework behind the omalizumab block.*
109. Gibiansky L, Gibiansky E, Kakkar T, Ma P. **Approximations of the target-mediated drug disposition model and identifiability of model parameters.** *J Pharmacokinet Pharmacodyn.* 2008;35(5):573-591. — https://pubmed.ncbi.nlm.nih.gov/18989757/ — *justification for the quasi-equilibrium closed form used in `bind11()`.*
110. Chen X, Hickling TP, Vicini P. **A mechanistic, multiscale mathematical model of immunogenicity for therapeutic proteins.** *CPT Pharmacometrics Syst Pharmacol.* 2014;3:e133. — https://pubmed.ncbi.nlm.nih.gov/25184733/
111. Kim EH, Burks AW. **Food allergy immunotherapy: oral immunotherapy and epicutaneous immunotherapy.** *Allergy.* 2020;75(6):1337-1346. — https://pubmed.ncbi.nlm.nih.gov/31840823/
112. Turner PJ, Patel N, Ballmer-Weber BK, et al. **Peanut can be used as a reference allergen for hazard characterization in food allergen risk management: a rapid evidence assessment and meta-analysis.** *J Allergy Clin Immunol Pract.* 2022;10(1):59-70. — https://pubmed.ncbi.nlm.nih.gov/34506967/ — *the pooled ED distribution against which `vpop()` is checked.*

---

## What this model does not cover, and where to read instead

| Excluded | Why | Read instead |
|---|---|---|
| FPIES | non-IgE, T-cell/innate, delayed emesis; no IgE anywhere in the mechanism | Nowak-Węgrzyn A, et al. *J Allergy Clin Immunol.* 2017;139(4):1111-1126 — https://pubmed.ncbi.nlm.nih.gov/28167094/ |
| Eosinophilic oesophagitis | food-triggered but not IgE-effector-mediated | Dellon ES, Hirano I. *Gastroenterology.* 2018;154(2):319-332 — https://pubmed.ncbi.nlm.nih.gov/28774845/ |
| Alpha-gal syndrome | IgE-mediated but delayed 3-6 h; the kinetics are chylomicron trafficking | Platts-Mills TAE, et al. *J Allergy Clin Immunol.* 2020;145(4):1061-1071 — https://pubmed.ncbi.nlm.nih.gov/32057766/ |
| MRGPRX2 / pseudoallergy | IgE-independent mast cell activation; a parallel pathway, drawn but not modelled | McNeil BD, et al. *Nature.* 2015;519(7542):237-241 — https://pubmed.ncbi.nlm.nih.gov/25517090/ |
| Clonal mast cell disease | enters here only as `MCBURDEN`; as a primary driver it needs its own model | Valent P, et al. *Blood.* 2017;129(11):1420-1427 — https://pubmed.ncbi.nlm.nih.gov/28031180/ |
