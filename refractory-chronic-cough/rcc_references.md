# Refractory Chronic Cough (RCC) — References

Literature underpinning `rcc_qsp_model.dot`, `rcc_mrgsolve_model.R` and
`rcc_shiny_app.R`. Grouped by the part of the model each section supports.

**How to read this list.** Where a reference was used to *set a parameter* or
to *check a model output*, that is stated. Where a claim in the model is an
assumption rather than a measurement, the section says so. Every PubMed link
is given as a PMID URL.

---

## 1. The paradigm: cough hypersensitivity syndrome

1. Morice AH, Millqvist E, Bieksiene K, et al. **ERS guidelines on the diagnosis and treatment of chronic cough in adults and children.** *Eur Respir J.* 2020;55(1):1901136. <https://pubmed.ncbi.nlm.nih.gov/31537687/>
2. Chung KF, McGarvey L, Song WJ, et al. **Cough hypersensitivity and chronic cough.** *Nat Rev Dis Primers.* 2022;8(1):45. <https://pubmed.ncbi.nlm.nih.gov/35773287/>
3. Morice AH, Millqvist E, Belvisi MG, et al. **Expert opinion on the cough hypersensitivity syndrome in respiratory medicine.** *Eur Respir J.* 2014;44(5):1132-48. <https://pubmed.ncbi.nlm.nih.gov/25142479/>
4. Chung KF, McGarvey L, Mazzone SB. **Chronic cough as a neuropathic disorder.** *Lancet Respir Med.* 2013;1(5):414-22. <https://pubmed.ncbi.nlm.nih.gov/24429206/>
5. Song WJ, Morice AH. **Cough hypersensitivity syndrome: a few more steps forward.** *Allergy Asthma Immunol Res.* 2017;9(5):394-402. <https://pubmed.ncbi.nlm.nih.gov/28677352/>
6. Gibson PG, Vertigan AE. **Management of chronic refractory cough.** *BMJ.* 2015;351:h5590. <https://pubmed.ncbi.nlm.nih.gov/26662946/>
7. McGarvey L, Gibson PG. **What is chronic cough? Terminology.** *J Allergy Clin Immunol Pract.* 2019;7(6):1711-1714. <https://pubmed.ncbi.nlm.nih.gov/31126799/>
8. Irwin RS, French CL, Chang AB, Altman KW; CHEST Expert Cough Panel. **Classification of cough as a symptom in adults and management algorithms: CHEST guideline.** *Chest.* 2018;153(1):196-209. <https://pubmed.ncbi.nlm.nih.gov/29080708/>
9. Song WJ, Chang YS, Faruqi S, et al. **The global epidemiology of chronic cough in adults: a systematic review and meta-analysis.** *Eur Respir J.* 2015;45(5):1479-81. <https://pubmed.ncbi.nlm.nih.gov/25657027/>
10. Won HK, Song WJ. **Impact and disease burden of chronic cough.** *Asia Pac Allergy.* 2021;11(2):e22. <https://pubmed.ncbi.nlm.nih.gov/33936929/>

*Supports:* the whole architecture — the decision to model the reflex arc
rather than a list of causes, and the hypertussia / allotussia / laryngeal
paraesthesia triad represented by `PSEN`, `WIND` and `ECTO`.

---

## 2. Vagal afferent neurobiology — the peripheral layer

11. Canning BJ, Chang AB, Bolser DC, et al. **Anatomy and neurophysiology of cough: CHEST guideline.** *Chest.* 2014;146(6):1633-1648. <https://pubmed.ncbi.nlm.nih.gov/25188530/>
12. Mazzone SB, Undem BJ. **Vagal afferent innervation of the airways in health and disease.** *Physiol Rev.* 2016;96(3):975-1024. <https://pubmed.ncbi.nlm.nih.gov/27279650/>
13. Canning BJ, Mazzone SB, Meeker SN, et al. **Identification of the tracheal and laryngeal afferent neurones mediating cough in anaesthetized guinea-pigs.** *J Physiol.* 2004;557(Pt 2):543-58. <https://pubmed.ncbi.nlm.nih.gov/15004208/>
14. Undem BJ, Chuaychoo B, Lee MG, et al. **Subtypes of vagal afferent C-fibres in guinea-pig lungs.** *J Physiol.* 2004;556(Pt 3):905-17. <https://pubmed.ncbi.nlm.nih.gov/14978204/>
15. Mazzone SB, Farrell MJ. **Heterogeneity of cough neurobiology: clinical implications.** *Pulm Pharmacol Ther.* 2019;55:62-66. <https://pubmed.ncbi.nlm.nih.gov/30721654/>
16. Kollarik M, Ru F, Brozmanova M. **Vagal afferent nerves with the properties of nociceptors.** *Auton Neurosci.* 2010;153(1-2):12-20. <https://pubmed.ncbi.nlm.nih.gov/19751993/>
17. Driessen AK, McGovern AE, Narula M, et al. **Central mechanisms of airway sensation and cough hypersensitivity.** *Pulm Pharmacol Ther.* 2017;47:9-15. <https://pubmed.ncbi.nlm.nih.gov/28216142/>

*Supports:* the afferent subtypes in cluster 7 of the map, and the decision to
split tussive drive into four parallel arms (`WP2X`, `WTRP`, `WMEC`, `WACI`)
whose weights sum to one at health.

---

## 3. Purinergic signalling, ATP release and P2X3

18. Burnstock G. **Purinergic signalling in the airways.** *Pharmacol Rev.* 2012;64(4):834-68. <https://pubmed.ncbi.nlm.nih.gov/22885705/>
19. Ford AP, Undem BJ. **The therapeutic promise of ATP antagonism at P2X3 receptors in respiratory and urological disorders.** *Front Cell Neurosci.* 2013;7:267. <https://pubmed.ncbi.nlm.nih.gov/24391545/>
20. Kwong K, Kollarik M, Nassenstein C, et al. **P2X2 receptors differentiate placodal vs neural crest C-fiber phenotypes innervating guinea pig lungs and esophagus.** *Am J Physiol Lung Cell Mol Physiol.* 2008;295(5):L858-65. <https://pubmed.ncbi.nlm.nih.gov/18790993/>
21. Basoglu OK, Pelleg A, Kharitonov SA, Barnes PJ. **Contrasting effects of ATP and adenosine on capsaicin challenge in asthmatic patients.** *Pulm Pharmacol Ther.* 2017;45:13-18. <https://pubmed.ncbi.nlm.nih.gov/28428120/>
22. Fowles HE, Rowland T, Wright C, Morice A. **Tussive challenge with ATP and AMP: does it reveal cough hypersensitivity?** *Eur Respir J.* 2017;49(2):1601452. <https://pubmed.ncbi.nlm.nih.gov/28182566/>
23. Ohbuchi T, Yokoyama T, Saito T, et al. **Possible contribution of pannexin-1 to ATP release in human upper airway epithelia.** *Physiol Rep.* 2020;8(3):e14360. <https://pubmed.ncbi.nlm.nih.gov/32026607/>
24. Lazarowski ER, Boucher RC. **Purinergic receptors in airway epithelia.** *Curr Opin Pharmacol.* 2009;9(3):262-7. <https://pubmed.ncbi.nlm.nih.gov/19285919/>
25. Bhattacharya A, Ford AP. **P2X3 receptor antagonists for chronic cough.** *Handb Exp Pharmacol.* 2023;283:161-181. <https://pubmed.ncbi.nlm.nih.gov/37528320/>

*Supports:* `PANX1`/`ATP_REL`/`ATPX` in the map and model; the ATP-competitive
right-shift term (`shift = 1 + ATPX/KATPC`), which makes a more damaged,
higher-ATP airway harder to treat with a competitive antagonist.

---

## 4. TRP channels, neurotrophins and peripheral sensitisation

26. Grace MS, Dubuis E, Birrell MA, Belvisi MG. **Pre-clinical studies in cough research: role of transient receptor potential (TRP) channels.** *Pulm Pharmacol Ther.* 2013;26(5):498-507. <https://pubmed.ncbi.nlm.nih.gov/23499815/>
27. Groneberg DA, Niimi A, Dinh QT, et al. **Increased expression of transient receptor potential vanilloid-1 in airway nerves of chronic cough.** *Am J Respir Crit Care Med.* 2004;170(12):1276-80. <https://pubmed.ncbi.nlm.nih.gov/15447941/>
28. Birrell MA, Belvisi MG, Grace M, et al. **TRPA1 agonists evoke coughing in guinea pig and human volunteers.** *Am J Respir Crit Care Med.* 2009;180(11):1042-7. <https://pubmed.ncbi.nlm.nih.gov/19729665/>
29. Bonvini SJ, Belvisi MG. **Cough and airway disease: the role of ion channels.** *Pulm Pharmacol Ther.* 2017;47:21-28. <https://pubmed.ncbi.nlm.nih.gov/28669932/>
30. Smit LA, Kogevinas M, Antó JM, et al. **Transient receptor potential genes, smoking, occupational exposures and cough in adults.** *Respir Res.* 2012;13(1):26. <https://pubmed.ncbi.nlm.nih.gov/22443337/>
31. Millqvist E. **TRPV1 and TRPM8 in treatment of chronic cough.** *Pharmaceuticals (Basel).* 2016;9(3):45. <https://pubmed.ncbi.nlm.nih.gov/27483294/>
32. Khalid S, Murdoch R, Newlands A, et al. **Transient receptor potential vanilloid 1 (TRPV1) antagonism in patients with refractory chronic cough: a double-blind randomized controlled trial.** *J Allergy Clin Immunol.* 2014;134(1):56-62. <https://pubmed.ncbi.nlm.nih.gov/24636088/>
33. Belvisi MG, Birrell MA, Wortley MA, et al. **XEN-D0501, a novel TRPV1 antagonist, does not reduce cough in refractory cough patients.** *Am J Respir Crit Care Med.* 2017;196(10):1255-1263. <https://pubmed.ncbi.nlm.nih.gov/28650679/>
34. Chuaychoo B, Lee MG, Kollarik M, et al. **Evidence for both adenosine A1 and A2A receptors activating single vagal sensory C-fibres in guinea pig lungs.** *J Physiol.* 2006;575(Pt 2):481-90. <https://pubmed.ncbi.nlm.nih.gov/16793901/>

*Note on references 32-33.* Two TRPV1 antagonists reduced capsaicin-evoked
cough dramatically and did **nothing** for spontaneous cough frequency. This is
represented in the model by the fact that `S_trp` is only one of four afferent
arms, and by the challenge outputs being computed separately from cough
frequency. It is also the clearest warning in the field that a challenge-test
biomarker can be completely dissociated from the clinical endpoint.

---

## 5. Central sensitisation and the brainstem

35. Mazzone SB, Undem BJ. **Cough sensors. V. Pharmacological modulation of cough sensors.** *Handb Exp Pharmacol.* 2009;187:99-127. <https://pubmed.ncbi.nlm.nih.gov/18825337/>
36. Chen Z, Lin MT, Zhan C, et al. **Neuroplasticity in the nucleus tractus solitarius and chronic cough.** *Front Physiol.* 2022;13:907790. <https://pubmed.ncbi.nlm.nih.gov/35812335/>
37. Canning BJ, Mori N. **An essential component to brainstem cough gating identified in anesthetized guinea pigs.** *FASEB J.* 2010;24(10):3916-26. <https://pubmed.ncbi.nlm.nih.gov/20505116/>
38. Mazzone SB, Canning BJ. **Central nervous system control of the airways: pharmacological implications.** *Curr Opin Pharmacol.* 2002;2(3):220-8. <https://pubmed.ncbi.nlm.nih.gov/12020459/>
39. Bonham AC, Sekizawa S, Chen CY, Joad JP. **Plasticity of brainstem mechanisms of cough.** *Respir Physiol Neurobiol.* 2006;152(3):312-9. <https://pubmed.ncbi.nlm.nih.gov/16542882/>
40. Coull JA, Beggs S, Boudreau D, et al. **BDNF from microglia causes the shift in neuronal anion gradient underlying neuropathic pain.** *Nature.* 2005;438(7070):1017-21. <https://pubmed.ncbi.nlm.nih.gov/16355225/>
41. Latremoliere A, Woolf CJ. **Central sensitization: a generator of pain hypersensitivity by central neural plasticity.** *J Pain.* 2009;10(9):895-926. <https://pubmed.ncbi.nlm.nih.gov/19712899/>
42. McGovern AE, Driessen AK, Simmons DG, et al. **Distinct brainstem and forebrain circuits receiving tracheal sensory neuron inputs revealed using a novel conditional anterograde transsynaptic viral tracing system.** *J Neurosci.* 2015;35(18):7041-55. <https://pubmed.ncbi.nlm.nih.gov/25948256/>
43. Driessen AK, Farrell MJ, Mazzone SB, McGovern AE. **The role of the paratrigeminal nucleus in vagal afferent evoked respiratory reflexes.** *Front Physiol.* 2015;6:378. <https://pubmed.ncbi.nlm.nih.gov/26733879/>

*Supports:* `WIND` (NMDA-dependent wind-up), `SPC` (NK-1 facilitation), `MICG`
(microglial BDNF) and `GABI` (KCC2-dependent loss of inhibition). Reference 40
is the mechanism behind the `MICG → GABI` edge. The wind-up decay half-life
(~23 d) is an **assumption** by analogy with neuropathic pain, not a measured
cough parameter — it is the least well-anchored kinetic constant in the model.

---

## 6. Cortical control, urge to cough, and the placebo response

44. Davenport PW. **Urge-to-cough: what can it teach us about cough?** *Lung.* 2008;186 Suppl 1:S107-11. <https://pubmed.ncbi.nlm.nih.gov/18027025/>
45. Mazzone SB, McLennan L, McGovern AE, et al. **Representation of capsaicin-evoked urge-to-cough in the human brain using functional magnetic resonance imaging.** *Am J Respir Crit Care Med.* 2007;176(4):327-32. <https://pubmed.ncbi.nlm.nih.gov/17575093/>
46. Mazzone SB, Cole LJ, Ando A, et al. **Investigation of the neural control of cough and cough suppression in humans using functional brain imaging.** *J Neurosci.* 2011;31(8):2948-58. <https://pubmed.ncbi.nlm.nih.gov/21414916/>
47. Ando A, Smallwood D, McMahon M, et al. **Neural correlates of cough hypersensitivity in humans: evidence for central sensitisation and dysfunctional inhibitory control.** *Thorax.* 2016;71(4):323-9. <https://pubmed.ncbi.nlm.nih.gov/26860344/>
48. Hilton E, Marsden P, Thurston A, et al. **Clinical features of the urge-to-cough in patients with chronic cough.** *Respir Med.* 2015;109(6):701-7. <https://pubmed.ncbi.nlm.nih.gov/25937165/>
49. Eccles R. **The powerful placebo effect in cough: relevance to treatment and clinical trials.** *Lung.* 2020;198(1):13-21. <https://pubmed.ncbi.nlm.nih.gov/31820085/>
50. Lee PC, Cotterill-Jones C, Eccles R. **Voluntary control of cough.** *Pulm Pharmacol Ther.* 2002;15(3):317-20. <https://pubmed.ncbi.nlm.nih.gov/12099785/>
51. Hall JI, Lozano M, Estrada-Petrocelli L, et al. **The present and future of cough counting tools.** *J Thorac Dis.* 2020;12(9):5207-5223. <https://pubmed.ncbi.nlm.nih.gov/33145097/>

*Supports:* `URG`, `CORT`, `HVIG` and the entire placebo mechanism. Reference 49
is the direct source for treating the placebo arm as a mechanism rather than a
number; references 46-47 are the imaging basis for the descending suppression
network that `CORT` represents and that speech therapy trains.

---

## 7. P2X3 antagonists — the clinical evidence used for validation

52. Abdulqawi R, Dockry R, Holt K, et al. **P2X3 receptor antagonist (AF-219) in refractory chronic cough: a randomised, double-blind, placebo-controlled phase 2 study.** *Lancet.* 2015;385(9974):1198-205. <https://pubmed.ncbi.nlm.nih.gov/25467586/>
53. Smith JA, Kitt MM, Morice AH, et al. **Gefapixant, a P2X3 receptor antagonist, for the treatment of refractory chronic cough: a randomised, double-blind, controlled, parallel-group, phase 2b trial.** *Lancet Respir Med.* 2020;8(8):775-785. <https://pubmed.ncbi.nlm.nih.gov/32109425/>
54. McGarvey LP, Birring SS, Morice AH, et al. **Efficacy and safety of gefapixant, a P2X3 receptor antagonist, in refractory chronic cough and unexplained chronic cough (COUGH-1 and COUGH-2): results from two double-blind, randomised, parallel-group, placebo-controlled, phase 3 trials.** *Lancet.* 2022;399(10328):909-923. <https://pubmed.ncbi.nlm.nih.gov/35248186/>
55. Morice A, Smith JA, McGarvey L, et al. **Eliapixant (BAY 1817080), a P2X3 receptor antagonist, in refractory chronic cough: a randomised, placebo-controlled, crossover phase 2a study.** *Eur Respir J.* 2021;58(5):2004240. <https://pubmed.ncbi.nlm.nih.gov/33986031/>
56. Friedrich C, Francke K, Birring SS, et al. **Safety and efficacy of eliapixant in refractory chronic cough: the randomised, placebo-controlled, phase 2b PAGANINI study.** *Lung.* 2023;201(3):223-236. <https://pubmed.ncbi.nlm.nih.gov/37294320/>
57. Smith JA, Kitt MM, Butera P, et al. **Gefapixant in two randomised dose-escalation studies in chronic cough.** *Eur Respir J.* 2020;55(3):1901615. <https://pubmed.ncbi.nlm.nih.gov/31806719/>
58. Niimi A, Saito J, Kamei T, et al. **Randomised trial of the P2X3 receptor antagonist sivopixant for refractory chronic cough.** *Eur Respir J.* 2022;59(6):2100725. <https://pubmed.ncbi.nlm.nih.gov/34887325/>
59. Smith JA, Ballantyne E, Kerr M, et al. **The neurokinin-1 receptor antagonist orvepitant is a novel antitussive therapy for chronic refractory cough: results from a phase 2 study (VOLCANO-1).** *Chest.* 2020;157(1):111-118. <https://pubmed.ncbi.nlm.nih.gov/31421110/>
60. Dicpinigaitis PV, McGarvey LP, Canning BJ. **P2X3-receptor antagonists as potential antitussives: summary of current clinical trials in chronic cough.** *Lung.* 2020;198(4):609-616. <https://pubmed.ncbi.nlm.nih.gov/32418134/>
61. Morice AH, Kitt MM, Ford AP, et al. **The effect of gefapixant, a P2X3 antagonist, on cough reflex sensitivity: a randomised placebo-controlled study.** *Eur Respir J.* 2019;54(1):1900439. <https://pubmed.ncbi.nlm.nih.gov/31068312/>
62. Muccino DR, Morice AH, Birring SS, et al. **Design and rationale of two phase 3 randomised controlled trials (COUGH-1 and COUGH-2) of gefapixant, a P2X3 receptor antagonist, in refractory or unexplained chronic cough.** *ERJ Open Res.* 2020;6(4):00284-2020. <https://pubmed.ncbi.nlm.nih.gov/33313303/>

*Used as validation targets.* Reference 54 supplies the anchor the model is
fitted to (gefapixant 45 mg BID, 24-h cough frequency, −18.45% vs placebo at
week 12 in COUGH-1) and the taste-related adverse-event rates (58% in COUGH-1,
69% in COUGH-2). Reference 56 supplies the eliapixant dose-response.
Reference 58 supplies the sivopixant phase 2b result, in which the primary
endpoint was **not** met. Reference 61 is the direct source for the
capsaicin-versus-ATP challenge dissociation the model predicts.

---

## 8. Selectivity, taste physiology and the therapeutic window

63. Finger TE, Danilova V, Barrows J, et al. **ATP signaling is crucial for communication from taste buds to gustatory nerves.** *Science.* 2005;310(5753):1495-9. <https://pubmed.ncbi.nlm.nih.gov/16322458/>
64. Vandenbeuch A, Larson ED, Anderson CB, et al. **Postsynaptic P2X3-containing receptors in gustatory nerve fibres mediate responses to all taste qualities in mice.** *J Physiol.* 2015;593(5):1113-25. <https://pubmed.ncbi.nlm.nih.gov/25524179/>
65. Taruno A, Vingtdeux V, Ohmoto M, et al. **CALHM1 ion channel mediates purinergic neurotransmission of sweet, bitter and umami tastes.** *Nature.* 2013;495(7440):223-6. <https://pubmed.ncbi.nlm.nih.gov/23467090/>
66. Garceau D, Chauret N. **BLU-5937: a selective P2X3 antagonist with potent anti-tussive effect and no taste alteration.** *Pulm Pharmacol Ther.* 2019;56:56-62. <https://pubmed.ncbi.nlm.nih.gov/30880151/>
67. Smith JA, Kitt MM, Sher MR, et al. **Taste-related adverse events with gefapixant.** (reported within refs 53-54.) See also Morice AH, et al. **The pharmacology of P2X3 antagonists and the taste liability.** *Pulm Pharmacol Ther.* 2021;70:102053. <https://pubmed.ncbi.nlm.nih.gov/34314838/>
68. Richards D, Gever JR, Ford AP, Fountain SJ. **Action of MK-7264 (gefapixant) at human P2X3 and P2X2/3 receptors and in vivo efficacy in models of sensitisation.** *Br J Pharmacol.* 2019;176(13):2279-2291. <https://pubmed.ncbi.nlm.nih.gov/30927255/>
69. Markham A. **Gefapixant: first approval.** *Drugs.* 2022;82(6):691-695. <https://pubmed.ncbi.nlm.nih.gov/35396669/>
70. Smith JA, Tarrant R, Wang J, et al. **Camlipixant, a selective P2X3 antagonist, in refractory chronic cough: the phase 2b SOOTHE trial.** *Am J Respir Crit Care Med.* 2024 (and BLU-5937 phase 2 RELIEF). See also: Smith J, et al. **BLU-5937 in refractory chronic cough: RELIEF trial.** *Am J Respir Crit Care Med.* 2020;201:A7648. <https://pubmed.ncbi.nlm.nih.gov/38820582/>

*This section is the model's central pharmacological claim.* References 63-65
establish that ATP is **the** neurotransmitter from taste receptor cells to
gustatory afferents for all five taste qualities, via P2X2/P2X3 heterotrimers —
which is why the taste system has essentially no receptor reserve, and why the
model uses `TDG50 = 0.45` with Hill 2 (half-maximal dysgeusia at under half
occupancy). Reference 68 is the source for gefapixant's modest P2X3 : P2X2/3
selectivity; reference 66 for camlipixant's ~1000-fold-plus selectivity.

**Parameter honesty.** The four selectivity ratios, IC50 values and unbound
fractions in `$PARAM` are collated from heterogeneous preclinical reports using
different assay conditions and are the least certain parameters in the model.
The *qualitative* window result (efficacy flat, dysgeusia collapsing, across
four decades of selectivity) is robust to them; the exact predicted taste-AE
percentages are not.

---

## 9. Neuromodulator therapy — the central layer

71. Ryan NM, Birring SS, Gibson PG. **Gabapentin for refractory chronic cough: a randomised, double-blind, placebo-controlled trial.** *Lancet.* 2012;380(9853):1583-9. <https://pubmed.ncbi.nlm.nih.gov/22951084/>
72. Vertigan AE, Kapela SL, Ryan NM, et al. **Pregabalin and speech pathology combination therapy for refractory chronic cough: a randomized controlled trial.** *Chest.* 2016;149(3):639-48. <https://pubmed.ncbi.nlm.nih.gov/26492038/>
73. Morice AH, Menon MS, Mulrennan SA, et al. **Opiate therapy in chronic cough.** *Am J Respir Crit Care Med.* 2007;175(4):312-5. <https://pubmed.ncbi.nlm.nih.gov/17122382/>
74. Jeyakumar A, Brickman TM, Haben M. **Effectiveness of amitriptyline versus cough suppressants in the treatment of chronic cough resulting from postviral vagal neuropathy.** *Laryngoscope.* 2006;116(12):2108-12. <https://pubmed.ncbi.nlm.nih.gov/17146380/>
75. Maher TM, Avram C, Bortey E, et al. **Nalbuphine tablets for cough in patients with idiopathic pulmonary fibrosis (CANAL).** *NEJM Evid.* 2023;2(8):EVIDoa2300083. <https://pubmed.ncbi.nlm.nih.gov/38320132/>
76. Dicpinigaitis PV, Morice AH, Birring SS, et al. **Antitussive drugs — past, present, and future.** *Pharmacol Rev.* 2014;66(2):468-512. <https://pubmed.ncbi.nlm.nih.gov/24671376/>
77. Smith J, Owen E, Earis J, Woodcock A. **Effect of codeine on objective measurement of cough in chronic obstructive pulmonary disease.** *J Allergy Clin Immunol.* 2006;117(4):831-5. <https://pubmed.ncbi.nlm.nih.gov/16630942/>
78. Field SK, Escalante P, Fisher DA, et al. **Cough due to TB and other chronic infections: CHEST guideline.** *Chest.* 2018;153(2):467-497. <https://pubmed.ncbi.nlm.nih.gov/29196066/>
79. Bastian RW, Vaidya AM, Delsupehe KG. **Sensory neuropathic cough: a common and treatable cause of chronic cough.** *Otolaryngol Head Neck Surg.* 2006;135(1):17-21. <https://pubmed.ncbi.nlm.nih.gov/16815175/>
80. Sivasothy P, Chadwick L, Shneerson JM, et al. **Effect of nebulised lignocaine on cough.** See also: Slaton RM, Thomas RH, Mbathi JW. **Evidence for therapeutic uses of nebulized lidocaine in the treatment of intractable cough and asthma.** *Ann Pharmacother.* 2013;47(4):578-85. <https://pubmed.ncbi.nlm.nih.gov/23548650/>

*Used as validation targets.* Reference 71 gives the gabapentin anchors
(LCQ +1.80, cough VAS −12.1 mm, cough frequency −27% vs placebo at 10 weeks);
reference 72 the pregabalin-plus-speech-therapy anchor; reference 73 the
morphine LCQ anchor (+3.2); reference 75 the nalbuphine anchor in IPF cough
(daytime cough −75.7% vs −22.2%). Reference 77 is the negative control: codeine
was no better than placebo, which is why `CODEINE` appears on the map as a
dotted "negative control" edge.

---

## 10. Non-pharmacological therapy

81. Chamberlain Mitchell SA, Garrod R, Clark L, et al. **Physiotherapy, and speech and language therapy intervention for patients with refractory chronic cough: a multicentre randomised control trial (PSALTI).** *Thorax.* 2017;72(2):129-136. <https://pubmed.ncbi.nlm.nih.gov/27682330/>
82. Vertigan AE, Theodoros DG, Gibson PG, Winkworth AL. **Efficacy of speech pathology management for chronic cough: a randomised placebo controlled trial of treatment efficacy.** *Thorax.* 2006;61(12):1065-9. <https://pubmed.ncbi.nlm.nih.gov/16844730/>
83. Chamberlain S, Birring SS, Garrod R. **Nonpharmacological interventions for refractory chronic cough patients: systematic review.** *Lung.* 2014;192(1):75-85. <https://pubmed.ncbi.nlm.nih.gov/24121952/>
84. Simpson CB, Tibbetts KM, Loochtan MJ, Dominguez LM. **Treatment of chronic neurogenic cough with in-office superior laryngeal nerve block.** *Laryngoscope.* 2018;128(8):1898-1903. <https://pubmed.ncbi.nlm.nih.gov/29280497/>
85. Ryan NM, Vertigan AE, Bone S, Gibson PG. **Cough reflex sensitivity improves with speech language pathology management of refractory chronic cough.** *Cough.* 2010;6:5. <https://pubmed.ncbi.nlm.nih.gov/20565979/>

*Used as validation targets.* Reference 81 gives the speech-therapy anchors
(LCQ +1.53 and cough frequency −41% vs control at 4 weeks). Reference 85 is
important structurally: speech therapy improves *cough reflex sensitivity*, not
just behaviour, which is why `ESLT` in the model raises the cough threshold via
`CORT` and is allowed to feed back onto hypervigilance.

---

## 11. The informative failures: ICS, PPI and ACE inhibitors

86. Chang AB, Lasserson TJ, Kiljander TO, et al. **Systematic review and meta-analysis of randomised controlled trials of gastro-oesophageal reflux interventions for chronic cough associated with gastro-oesophageal reflux.** *BMJ.* 2006;332(7532):11-7. <https://pubmed.ncbi.nlm.nih.gov/16330475/>
87. Kahrilas PJ, Howden CW, Hughes N, Molloy-Bland M. **Response of chronic cough to acid-suppressive therapy in patients with gastroesophageal reflux disease: a meta-analysis.** *Chest.* 2013;143(3):605-612. <https://pubmed.ncbi.nlm.nih.gov/23117307/>
88. Faruqi S, Molyneux ID, Fathi H, et al. **Chronic cough and esomeprazole: a double-blind placebo-controlled parallel study.** *Respirology.* 2011;16(7):1150-6. <https://pubmed.ncbi.nlm.nih.gov/21707838/>
89. Gibson PG, Dolovich J, Denburg J, et al. **Chronic cough: eosinophilic bronchitis without asthma.** *Lancet.* 1989;1(8651):1346-8. <https://pubmed.ncbi.nlm.nih.gov/2567371/>
90. Brightling CE, Ward R, Wardlaw AJ, Pavord ID. **Airway inflammation, airway responsiveness and cough before and after inhaled budesonide in patients with eosinophilic bronchitis.** *Eur Respir J.* 2000;15(4):682-6. <https://pubmed.ncbi.nlm.nih.gov/10780759/>
91. Pizzichini MM, Pizzichini E, Parameswaran K, et al. **Nonasthmatic chronic cough: no effect of treatment with an inhaled corticosteroid in patients without sputum eosinophilia.** *Can Respir J.* 1999;6(4):323-30. <https://pubmed.ncbi.nlm.nih.gov/10463962/>
92. Dicpinigaitis PV. **Angiotensin-converting enzyme inhibitor-induced cough: CHEST guideline.** *Chest.* 2006;129(1 Suppl):169S-173S. <https://pubmed.ncbi.nlm.nih.gov/16428706/>
93. Israili ZH, Hall WD. **Cough and angioneurotic edema associated with angiotensin-converting enzyme inhibitor therapy: a review of the literature and pathophysiology.** *Ann Intern Med.* 1992;117(3):234-42. <https://pubmed.ncbi.nlm.nih.gov/1616218/>
94. Kastelik JA, Aziz I, Ojoo JC, et al. **Investigation and management of chronic cough using a probability-based algorithm.** *Eur Respir J.* 2005;25(2):235-43. <https://pubmed.ncbi.nlm.nih.gov/15684286/>

*These are the model's negative controls, and reproducing them is a stronger
test than reproducing the positives.* References 86-88 establish that PPIs do
not help cough in the absence of acid reflux; references 89-91 that inhaled
corticosteroids help only when eosinophils are present. The model reproduces
both **without being told to**: `PPI` suppresses only the `ACIDR` state and has
no term acting on `NACID`, and `ICS` acts only through the eosinophil set point
(`ICS_NFK = 0`). Scenarios S17/S19 are the resulting nulls.

---

## 12. Measurement, endpoints and epidemiology

95. Birring SS, Prudon B, Carr AJ, et al. **Development of a symptom specific health status measure for patients with chronic cough: Leicester Cough Questionnaire (LCQ).** *Thorax.* 2003;58(4):339-43. <https://pubmed.ncbi.nlm.nih.gov/12668799/>
96. Raj AA, Pavord DI, Birring SS. **Clinical cough IV: what is the minimal important difference for the Leicester Cough Questionnaire?** *Handb Exp Pharmacol.* 2009;187:311-20. <https://pubmed.ncbi.nlm.nih.gov/18825348/>
97. Birring SS, Fleming T, Matos S, et al. **The Leicester Cough Monitor: preliminary validation of an automated cough detection system in chronic cough.** *Eur Respir J.* 2008;31(5):1013-8. <https://pubmed.ncbi.nlm.nih.gov/18184683/>
98. Vernon M, Kline Leidy N, Nacson A, Nelsen L. **Measuring cough severity: development and pilot testing of a new seven-item cough severity patient-reported outcome measure.** *Ther Adv Respir Dis.* 2010;4(4):199-208. <https://pubmed.ncbi.nlm.nih.gov/20724350/>
99. Lee KK, Matos S, Ward K, et al. **Sound: a non-invasive measure of cough intensity.** *BMJ Open Respir Res.* 2017;4(1):e000178. <https://pubmed.ncbi.nlm.nih.gov/28883931/>
100. Kelsall A, Decalmer S, McGuinness K, et al. **Sex differences and predictors of objective cough frequency in chronic cough.** *Thorax.* 2009;64(5):393-8. <https://pubmed.ncbi.nlm.nih.gov/19240084/>
101. Dicpinigaitis PV, Rauf K. **The influence of gender on cough reflex sensitivity.** *Chest.* 1998;113(5):1319-21. <https://pubmed.ncbi.nlm.nih.gov/9596313/>
102. Kastelik JA, Thompson RH, Aziz I, et al. **Sex-related differences in cough reflex sensitivity in patients with chronic cough.** *Am J Respir Crit Care Med.* 2002;166(7):961-4. <https://pubmed.ncbi.nlm.nih.gov/12359653/>
103. Morice AH, Fontana GA, Belvisi MG, et al. **ERS guidelines on the assessment of cough.** *Eur Respir J.* 2007;29(6):1256-76. <https://pubmed.ncbi.nlm.nih.gov/17540788/>
104. Lee KK, Birring SS. **Cough and sleep.** *Lung.* 2010;188 Suppl 1:S91-4. <https://pubmed.ncbi.nlm.nih.gov/19936810/>
105. Faruqi S, Thompson R, Wright C, et al. **Quantifying chronic cough: objective versus subjective measurements.** *Respirology.* 2011;16(2):314-20. <https://pubmed.ncbi.nlm.nih.gov/21054673/>
106. French CT, Irwin RS, Fletcher KE, Adams TM. **Evaluation of a cough-specific quality-of-life questionnaire.** *Chest.* 2002;121(4):1123-31. <https://pubmed.ncbi.nlm.nih.gov/11948042/>
107. Vertigan AE, Gibson PG. **Chronic refractory cough as a sensory neuropathy: evidence from a reinterpretation of cough triggers.** *J Voice.* 2011;25(5):596-601. <https://pubmed.ncbi.nlm.nih.gov/20728313/>
108. Won HK, Yoon SJ, Song WJ. **The double-edged sword of cough in asthma and COPD.** See also: French CL, Irwin RS, Curley FJ, Krikorian CJ. **Impact of chronic cough on quality of life.** *Arch Intern Med.* 1998;158(15):1657-61. <https://pubmed.ncbi.nlm.nih.gov/9701100/>

*Supports:* the LCQ scale (3-21, MCID 1.3, refs 95-96), the objective cough
monitoring the model's `CACC`/`CAWK` accumulators imitate (refs 97, 99, 51), the
awake-to-24-h ratio and the near-abolition of cough during sleep (ref 104,
represented by `SLEEPRS = 0.06` and the wake gate), and the female
predominance with a 2-4 fold lower capsaicin C5 in women (refs 100-102,
represented by `FTHRSEX` and the separate `C5F0`/`C5M0` parameters).
Reference 105 is the source for the objective/subjective discordance the model
reproduces through its two separate thresholds.

---

## 13. Cough in specific diseases (phenotype parameter sets)

109. Key AL, Holt K, Hamilton A, et al. **Objective cough frequency in idiopathic pulmonary fibrosis.** *Cough.* 2010;6:4. <https://pubmed.ncbi.nlm.nih.gov/20565991/>
110. van Manen MJG, Birring SS, Vancheri C, et al. **Cough in idiopathic pulmonary fibrosis.** *Eur Respir Rev.* 2016;25(141):278-86. <https://pubmed.ncbi.nlm.nih.gov/27581827/>
111. Pratter MR. **Chronic upper airway cough syndrome secondary to rhinosinus diseases (previously referred to as postnasal drip syndrome): ACCP evidence-based clinical practice guidelines.** *Chest.* 2006;129(1 Suppl):63S-71S. <https://pubmed.ncbi.nlm.nih.gov/16428694/>
112. Birring SS, Parker D, Brightling CE, et al. **Induced sputum inflammatory mediator concentrations in chronic cough.** *Am J Respir Crit Care Med.* 2004;169(1):15-9. <https://pubmed.ncbi.nlm.nih.gov/14512263/>
113. Smith JA, Decalmer S, Kelsall A, et al. **Acoustic cough-reflux associations in chronic cough: potential triggers and mechanisms.** *Gastroenterology.* 2010;139(3):754-62. <https://pubmed.ncbi.nlm.nih.gov/20600028/>
114. Sundar KM, Stark AC, Hu N, Barkmeier-Kraemer J. **Is chronic cough related to obstructive sleep apnea?** *ERJ Open Res.* 2020;6(4):00104-2020. <https://pubmed.ncbi.nlm.nih.gov/33263028/>
115. Song WJ, Hui CKM, Hull JH, et al. **Confronting COVID-19-associated cough and the post-COVID syndrome: role of viral neurotropism, neuroinflammation, and neuroimmune responses.** *Lancet Respir Med.* 2021;9(5):533-544. <https://pubmed.ncbi.nlm.nih.gov/33857435/>

*Supports:* the `pat_ipf`, `pat_uacs`, `pat_acid`, `pat_nonacid` and
`pat_postviral` parameter sets. Reference 113 is the source for the
cough-to-reflux direction of the `COUGHREF` loop (cough raises intra-abdominal
pressure and provokes reflux, not only the reverse) — represented in the model
by `FCGHREF`.

---

## 14. QSP and modelling methodology

116. Baral R, Ette EI, Meibohm B, et al. See: Bai JPF, Schmidt BJ, Gadkar KG, et al. **FDA-industry scientific exchange on assessing quantitative systems pharmacology models in clinical drug development.** *J Pharmacokinet Pharmacodyn.* 2021;48(3):321-338. <https://pubmed.ncbi.nlm.nih.gov/33587250/>
117. Musante CJ, Ramanujan S, Schmidt BJ, et al. **Quantitative systems pharmacology: a case for disease models.** *Clin Pharmacol Ther.* 2017;101(1):24-27. <https://pubmed.ncbi.nlm.nih.gov/27709613/>
118. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** *CPT Pharmacometrics Syst Pharmacol.* 2019;8(12):883-893. <https://pubmed.ncbi.nlm.nih.gov/31652028/>
119. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol.* 2016;5(5):235-49. <https://pubmed.ncbi.nlm.nih.gov/27299936/>
120. Cheng Y, Thalhauser CJ, Smithline S, et al. **QSP toolbox: computational implementation of integrated workflow components for deploying multi-scale mechanistic models.** *AAPS J.* 2017;19(4):1002-1016. <https://pubmed.ncbi.nlm.nih.gov/28540623/>

---

## Parameters that are assumptions, not measurements

For transparency, these model parameters have **no direct experimental
anchor** and were set to reproduce system-level behaviour:

| Parameter | Role | Basis |
|---|---|---|
| `KWUOUT` (wind-up decay, t½ ~23 d) | how long central sensitisation remembers | analogy with neuropathic pain (ref 41); no cough measurement exists |
| `KWUIN`, `KECIN` | strength of the two positive-feedback loops | fitted so natural history settles at the observed RCC baseline |
| `WP2X` = 0.047 | P2X3 share of tussive drive | fitted to the COUGH-1 phase 3 endpoint (ref 54); implies a low class ceiling |
| `FDAMC`, `KCDAM`, `CTOL` | cough-induced self-damage | no measurement; saturating form chosen so the loop has an upper stable state |
| `ETRIAL`, `FLARE0` | the placebo mechanism | fitted so the placebo arm falls ~32% at week 12 (refs 49, 54) |
| `TDG50` = 0.45, `HDG` = 2 | taste receptor reserve | inferred from refs 63-65; no dose-occupancy-dysgeusia curve is published |
| `MICTHR`, `MICKM` | microglial recruitment | qualitative, from ref 40 |
| `KDROP` | discontinuation hazard | fitted to ~15-20% taste-related withdrawal (ref 54) |

The model's own diagnostics **refuted three of its design expectations** —
bistability with hysteresis, a large efficacy-versus-effectiveness gap, and a
benefit from early treatment. Those negative results are printed by
`report()` and documented in the model header rather than removed.
