# Malignant Hyperthermia — References

Literature underpinning the QSP model in this directory (`mh_qsp_model.dot`,
`mh_mrgsolve_model.R`, `mh_shiny_app.R`). Grouped by the part of the model each
source constrains, so that any parameter can be traced back to its evidence.

Links are PubMed (`https://pubmed.ncbi.nlm.nih.gov/<PMID>/`). Where a source is
a guideline or registry document without a PMID, the canonical URL is given.

---

## 1. Foundational clinical description and epidemiology

1. Denborough MA, Lovell RRH. **Anaesthetic deaths in a family.** Lancet 1960;2:45.
   The index Melbourne pedigree that established MH as a heritable
   anaesthetic disease. <https://pubmed.ncbi.nlm.nih.gov/13158517/>
2. Rosenberg H, Pollock N, Schiemann A, Bulger T, Stowell K. **Malignant
   hyperthermia: a review.** Orphanet J Rare Dis 2015;10:93.
   <https://pubmed.ncbi.nlm.nih.gov/26238698/>
3. Rosenberg H, Sambuughin N, Riazi S, Dirksen R. **Malignant Hyperthermia
   Susceptibility.** GeneReviews, updated 2020.
   <https://pubmed.ncbi.nlm.nih.gov/20301325/>
4. Riazi S, Kraeva N, Hopkins PM. **Malignant hyperthermia in the post-genomics
   era.** Anesthesiology 2018;128:168-180.
   <https://pubmed.ncbi.nlm.nih.gov/28902675/>
5. Larach MG, Brandom BW, Allen GC, Gronert GA, Lehman EB. **Malignant
   hyperthermia deaths related to inadequate temperature monitoring, 2007-2012.**
   Anesth Analg 2014;119:1359-1366.
   <https://pubmed.ncbi.nlm.nih.gov/25268394/>
6. Brady JE, Sun LS, Rosenberg H, Li G. **Prevalence of malignant hyperthermia
   due to anesthesia in New York State, 2001-2005.** Anesth Analg
   2009;109:1162-1166. <https://pubmed.ncbi.nlm.nih.gov/19762744/>
7. Rosero EB, Adesanya AO, Timaran CH, Joshi GP. **Trends and outcomes of
   malignant hyperthermia in the United States, 2000 to 2005.** Anesthesiology
   2009;110:89-94. <https://pubmed.ncbi.nlm.nih.gov/19104175/>
8. Ording H. **Incidence of malignant hyperthermia in Denmark.** Anesth Analg
   1985;64:700-704. <https://pubmed.ncbi.nlm.nih.gov/4014716/>
9. Monnier N, Krivosic-Horber R, Payen JF, et al. **Presence of two different
   genetic traits in malignant hyperthermia families.** Anesthesiology
   2002;97:1067-1074. <https://pubmed.ncbi.nlm.nih.gov/12411788/>
10. Gonsalves SG, Ng D, Johnston JJ, et al. **Using exome data to identify
    malignant hyperthermia susceptibility mutations.** Anesthesiology
    2013;119:1043-1053. <https://pubmed.ncbi.nlm.nih.gov/24195946/>
    Basis for the model's statement that genetic prevalence (~1:2000-1:3000)
    vastly exceeds clinical event rate — i.e. MH is a threshold phenomenon.

---

## 2. Genotype — the parameters `EC50_VOL`, `PO_BASE`, `KSOCE`, `EC50_CAF`

11. MacLennan DH, Duff C, Zorzato F, et al. **Ryanodine receptor gene is a
    candidate for predisposition to malignant hyperthermia.** Nature
    1990;343:559-561. <https://pubmed.ncbi.nlm.nih.gov/1967823/>
12. Fujii J, Otsu K, Zorzato F, et al. **Identification of a mutation in porcine
    ryanodine receptor associated with malignant hyperthermia.** Science
    1991;253:448-451. <https://pubmed.ncbi.nlm.nih.gov/1862346/>
    The R614C/R615C allele that the model's `MHS_high` genotype represents.
13. Gillard EF, Otsu K, Fujii J, et al. **A substitution of cysteine for
    arginine 614 in the ryanodine receptor is potentially causative of human
    malignant hyperthermia.** Genomics 1991;11:751-755.
    <https://pubmed.ncbi.nlm.nih.gov/1774074/>
14. Robinson R, Carpenter D, Shaw MA, Halsall J, Hopkins P. **Mutations in RYR1
    in malignant hyperthermia and central core disease.** Hum Mutat
    2006;27:977-989. <https://pubmed.ncbi.nlm.nih.gov/16917943/>
    Source for the three mutational hot-spots drawn in cluster 1.
15. Monnier N, Procaccio V, Stieglitz P, Lunardi J. **Malignant-hyperthermia
    susceptibility is associated with a mutation of the alpha1-subunit of the
    human dihydropyridine-sensitive L-type voltage-dependent calcium-channel
    receptor in skeletal muscle.** Am J Hum Genet 1997;60:1316-1325.
    <https://pubmed.ncbi.nlm.nih.gov/9199552/>
    CACNA1S p.Arg1086His — the model's `MHS_CACNA1S` genotype.
16. Horstick EJ, Linsley JW, Dowling JJ, et al. **Stac3 is a component of the
    excitation-contraction coupling machinery and mutated in Native American
    myopathy.** Nat Commun 2013;4:1952.
    <https://pubmed.ncbi.nlm.nih.gov/23736855/>
17. Riazi S, Kraeva N, Muldoon SM, et al. **Malignant hyperthermia and the
    clinical significance of type-1 ryanodine receptor gene (RYR1) variants:
    proceedings of the 2013 MHAUS Scientific Conference.** Can J Anaesth
    2014;61:1040-1049. <https://pubmed.ncbi.nlm.nih.gov/25164374/>
18. Johnston JJ, Dirksen RT, Girard T, et al. **Updated variant curation
    expert panel criteria and pathogenicity classifications for 251 variants
    for RYR1-related malignant hyperthermia susceptibility.** Hum Mol Genet
    2022;31:4087-4093. <https://pubmed.ncbi.nlm.nih.gov/35849052/>
19. Lawal TA, Wires ES, Terry NL, Dowling JJ, Todd JJ. **Preclinical model
    systems of ryanodine receptor 1-related myopathies and malignant
    hyperthermia: a comprehensive scoping review.** Orphanet J Rare Dis
    2020;15:113. <https://pubmed.ncbi.nlm.nih.gov/32381029/>
20. Carpenter D, Robinson RL, Quinnell RJ, et al. **Genetic variation in RYR1
    and malignant hyperthermia phenotypes.** Br J Anaesth 2009;103:538-548.
    <https://pubmed.ncbi.nlm.nih.gov/19648156/>
    Evidence that variant identity predicts IVCT phenotype severity — the
    basis for genotype-specific `EC50_VOL` in this model.

---

## 3. RyR1 gating, calcium handling, and the store-operated entry thesis

21. Lopez JR, Alamo L, Caputo C, Wikinski J, Ledezma D. **Intracellular ionized
    calcium concentration in muscles from humans with malignant hyperthermia.**
    Muscle Nerve 1985;8:355-358. <https://pubmed.ncbi.nlm.nih.gov/4033585/>
    Source for the model's elevated *resting* myoplasmic Ca²⁺ in MHS
    (parameter `PO_BASE`).
22. Lopez JR, Allen PD, Alamo L, Jones D, Sreter FA. **Myoplasmic free [Ca2+]
    during a malignant hyperthermia episode in swine.** Muscle Nerve
    1988;11:82-88. <https://pubmed.ncbi.nlm.nih.gov/3343766/>
23. Cherednichenko G, Ward CW, Feng W, et al. **Enhanced excitation-coupled
    calcium entry in myotubes expressing malignant hyperthermia mutation R163C
    is attenuated by dantrolene.** Mol Pharmacol 2008;73:1203-1212.
    <https://pubmed.ncbi.nlm.nih.gov/18171728/>
24. Yang T, Allen PD, Pessah IN, Lopez JR. **Enhanced excitation-coupled
    calcium entry in myotubes is associated with expression of RyR1 malignant
    hyperthermia mutations.** J Biol Chem 2007;282:37471-37478.
    <https://pubmed.ncbi.nlm.nih.gov/17942409/>
25. Eltit JM, Ding X, Pessah IN, Allen PD, Lopez JR. **Nonspecific sarcolemmal
    cation channels are critical for the pathogenesis of malignant
    hyperthermia.** FASEB J 2013;27:991-1000.
    <https://pubmed.ncbi.nlm.nih.gov/23159934/>
    **The key source for the model's fourth thesis**: the sustained phase of MH
    requires trans-sarcolemmal Ca²⁺ entry, not merely SR release. The model's
    prediction that removing bath Ca²⁺ abolishes the contracture derives from
    this line of work.
26. Duke AM, Hopkins PM, Calaghan SC, Halsall JP, Steele DS. **Store-operated
    Ca2+ entry in malignant hyperthermia-susceptible human skeletal muscle.**
    J Biol Chem 2010;285:25645-25653.
    <https://pubmed.ncbi.nlm.nih.gov/20551327/>
    Direct evidence for the genotype-linked `KSOCE` parameter.
27. Michelucci A, Boncompagni S, Canato M, Reggiani C, Protasi F. **Estrogens
    protect calsequestrin-1 knockout mice from lethal hyperthermic episodes by
    reducing oxidative stress in muscle.** Oxid Med Cell Longev
    2017;2017:6936897. <https://pubmed.ncbi.nlm.nih.gov/29147462/>
28. Protasi F, Paolini C, Dainese M. **Calsequestrin-1: a new candidate gene
    for malignant hyperthermia and exertional/environmental heat stroke.**
    J Physiol 2009;587:3095-3100. <https://pubmed.ncbi.nlm.nih.gov/19417098/>
29. Zvaritch E, Depreux F, Kraeva N, et al. **An Ryr1(I4895T) mutation abolishes
    Ca2+ release channel function in a heterozygous mouse model.** Proc Natl
    Acad Sci USA 2007;104:18537-18542.
    <https://pubmed.ncbi.nlm.nih.gov/18003898/>
30. Meissner G. **The structural basis of ryanodine receptor ion channel
    function.** J Gen Physiol 2017;149:1065-1089.
    <https://pubmed.ncbi.nlm.nih.gov/29122978/>
31. des Georges A, Clarke OB, Zalk R, et al. **Structural basis for gating and
    activation of RyR1.** Cell 2016;167:145-157.
    <https://pubmed.ncbi.nlm.nih.gov/27662087/>
32. Nelson TE. **Halothane effects on human malignant hyperthermia skeletal
    muscle single calcium-release channels in planar lipid bilayers.**
    Anesthesiology 1992;76:588-595.
    <https://pubmed.ncbi.nlm.nih.gov/1312805/>
    Single-channel basis for the model's `PO_MAX` and the volatile-agent
    activation term.

---

## 4. Redox sensitisation — the model's `SENS` state (the SUSTAIN term)

33. Durham WJ, Aracena-Parks P, Long C, et al. **RyR1 S-nitrosylation underlies
    environmental heat stroke and sudden death in Y522S RyR1 knockin mice.**
    Cell 2008;133:53-65. <https://pubmed.ncbi.nlm.nih.gov/18394989/>
    **The primary source for the `SENS` state.** Establishes that oxidative /
    nitrosative modification of RyR1 is itself a sensitising mechanism, and that
    it can drive a hyperthermic crisis with no volatile agent present.
34. Aracena-Parks P, Goonasekera SA, Gilman CP, Dirksen RT, Hidalgo C, Hamilton
    SL. **Identification of cysteines involved in S-nitrosylation,
    S-glutathionylation, and oxidation to disulfides in ryanodine receptor
    type 1.** J Biol Chem 2006;281:40354-40368.
    <https://pubmed.ncbi.nlm.nih.gov/17071618/>
35. Michelucci A, Paolini C, Canato M, et al. **Antioxidants protect
    calsequestrin-1 knockout mice from halothane- and heat-induced sudden
    death.** Anesthesiology 2015;123:603-617.
    <https://pubmed.ncbi.nlm.nih.gov/26120769/>
    Direct support for the model's N-acetylcysteine probe (`NAC`).
36. Canato M, Capitanio P, Cancellara L, et al. **Excessive accumulation of Ca2+
    in mitochondria of Y522S-RYR1 knock-in mice: a link between leak from the
    sarcoplasmic reticulum and altered redox state.** Front Physiol 2019;10:1142.
    <https://pubmed.ncbi.nlm.nih.gov/31551814/>
    Basis for the model's mitochondrial-Ca²⁺-overload / uncoupling term
    (`CAMITO50`, `MITOEFF`).
37. Chang L, Daly C, Miller DM, et al. **Permeabilised skeletal muscle reveals
    mitochondrial deficiency in malignant hyperthermia-susceptible individuals.**
    Br J Anaesth 2019;122:613-621.
    <https://pubmed.ncbi.nlm.nih.gov/30916032/>
38. Giulivi C, Ross-Inta C, Omanska-Klusek A, et al. **Basal bioenergetic
    abnormalities in skeletal muscle from ryanodine receptor malignant
    hyperthermia-susceptible R163C knock-in mice.** J Biol Chem
    2011;286:99-113. <https://pubmed.ncbi.nlm.nih.gov/21051539/>
39. Todd JJ, Lawal TA, Chrismer IC, et al. **Randomized controlled trial of
    N-acetylcysteine therapy for RYR1-related myopathies.** Neurology
    2020;94:e1434-e1444. <https://pubmed.ncbi.nlm.nih.gov/32102975/>
    The clinical test of the antioxidant limb; informs the model's honest
    statement that the `NAC` probe is exploratory rather than established.
40. Lawal TA, Todd JJ, Meilleur KG. **Ryanodine receptor 1-related myopathies:
    diagnostic and therapeutic approaches.** Neurotherapeutics 2018;15:885-899.
    <https://pubmed.ncbi.nlm.nih.gov/30406384/>

---

## 5. Triggering agents, MAC, and volatile pharmacokinetics

41. Eger EI 2nd. **Uptake and distribution.** In: Miller's Anesthesia — the
    canonical treatment of the multi-compartment volatile-agent model this
    file implements. Underlying primary data: Eger EI 2nd, Saidman LJ.
    **Illustrations of inhaled anesthetic uptake, including intertissue
    diffusion to and from fat.** Anesth Analg 2005;100:1020-1033.
    <https://pubmed.ncbi.nlm.nih.gov/15781516/>
42. Yasuda N, Lockhart SH, Eger EI 2nd, et al. **Comparison of kinetics of
    sevoflurane and isoflurane in humans.** Anesth Analg 1991;72:316-324.
    <https://pubmed.ncbi.nlm.nih.gov/1994760/>
    Source for the blood:gas and tissue:blood partition coefficients.
43. Wappler F, Fiege M, Schulte am Esch J. **Pathophysiological role of the
    serotonin system in malignant hyperthermia.** Br J Anaesth 2001;87:794-798.
    <https://pubmed.ncbi.nlm.nih.gov/11878537/>
44. Wedel DJ, Gammel SA, Milde JH, Iaizzo PA. **Delayed onset of malignant
    hyperthermia induced by isoflurane and desflurane compared with halothane
    in susceptible swine.** Anesthesiology 1993;78:1138-1144.
    <https://pubmed.ncbi.nlm.nih.gov/8512107/>
    **The direct source for the model's agent-potency ranking.** Establishes
    that the newer, less soluble agents give a materially later onset — which
    the model reproduces from potency-per-vol% rather than by assertion.
45. Wedel DJ, Iaizzo PA, Milde JH. **Desflurane is a trigger of malignant
    hyperthermia in susceptible swine.** Anesthesiology 1991;74:508-512.
    <https://pubmed.ncbi.nlm.nih.gov/2001030/>
46. Hopkins PM. **Malignant hyperthermia: pharmacology of triggering.**
    Br J Anaesth 2011;107:48-56. <https://pubmed.ncbi.nlm.nih.gov/21624965/>
47. Larach MG, Gronert GA, Allen GC, Brandom BW, Lehman EB. **Clinical
    presentation, treatment, and complications of malignant hyperthermia in
    North America from 1987 to 2006.** Anesth Analg 2010;110:498-507.
    <https://pubmed.ncbi.nlm.nih.gov/19917623/>
    **The central clinical calibration target of this model.** Source for the
    ordering of presenting signs, the dominance of hypercarbia as the first
    sign, and the finding that each 30 min of delay to dantrolene multiplied
    the odds of a complication by 1.61.
48. Visoiu M, Young MC, Wieland K, Brandom BW. **Anesthetic drugs and onset of
    malignant hyperthermia.** Anesth Analg 2014;118:388-396.
    <https://pubmed.ncbi.nlm.nih.gov/24361847/>
49. Kim TW, Nemergut ME. **Preparation of modern anesthesia workstations for
    malignant hyperthermia-susceptible patients.** Anesthesiology
    2011;114:205-212. <https://pubmed.ncbi.nlm.nih.gov/21150570/>
    Source for the circuit washout time constants.
50. Birgenheier N, Stoker R, Westenskow D, Orr J. **Activated charcoal
    effectively removes inhaled anesthetics from modern anesthesia machines.**
    Anesth Analg 2011;112:1363-1370.
    <https://pubmed.ncbi.nlm.nih.gov/21543780/>
    The <5 ppm-in-90-s figure the model is compared against.
51. Whitty RJ, Wong GK, Petroz GC, Pehora C, Crawford MW. **Preparation of the
    Dräger Primus anesthetic machine for malignant hyperthermia-susceptible
    patients.** Can J Anaesth 2009;56:497-501.
    <https://pubmed.ncbi.nlm.nih.gov/19424761/>
52. Gunter JB, Ball J, Than-Win S. **Preparation of the Dräger Fabius anesthesia
    machine for the malignant-hyperthermia susceptible patient.** Anesth Analg
    2008;107:1936-1945. <https://pubmed.ncbi.nlm.nih.gov/19020140/>

---

## 6. Succinylcholine, masseter rigidity, and the differential

53. Larach MG, Rosenberg H, Larach DR, Broennle AM. **Prediction of malignant
    hyperthermia susceptibility by clinical signs.** Anesthesiology
    1987;66:547-550. <https://pubmed.ncbi.nlm.nih.gov/3565822/>
54. O'Flynn RP, Shutack JG, Rosenberg H, Fletcher JE. **Masseter muscle rigidity
    and malignant hyperthermia susceptibility in pediatric patients.**
    Anesthesiology 1994;80:1228-1233.
    <https://pubmed.ncbi.nlm.nih.gov/8010468/>
55. Littleford JA, Patel LR, Bose D, Cameron CB, McKillop C. **Masseter muscle
    spasm in children: implications of continuing anesthesia.** Anesth Analg
    1991;72:151-160. <https://pubmed.ncbi.nlm.nih.gov/1989533/>
56. Gronert GA, Thompson RL, Onofrio BM. **Human malignant hyperthermia: awake
    episodes and correction by dantrolene.** Anesth Analg 1980;59:377-378.
    <https://pubmed.ncbi.nlm.nih.gov/7189983/>
57. Larach MG, Klumpner TT, Brandom BW, et al. **Succinylcholine use and
    dantrolene availability for malignant hyperthermia treatment.**
    Anesthesiology 2019;130:41-54.
    <https://pubmed.ncbi.nlm.nih.gov/30489322/>
58. Gurnaney H, Brown A, Litman RS. **Malignant hyperthermia and muscular
    dystrophies.** Anesth Analg 2009;109:1043-1048.
    <https://pubmed.ncbi.nlm.nih.gov/19762730/>
    The basis for the model's insistence that dystrophinopathic
    succinylcholine-induced hyperkalaemic arrest is *not* MH — a fragile
    membrane rather than a leaky channel.
59. Gronert GA. **Cardiac arrest after succinylcholine: mortality greater with
    rhabdomyolysis than receptor upregulation.** Anesthesiology
    2001;94:523-529. <https://pubmed.ncbi.nlm.nih.gov/11374615/>
60. Ali SZ, Taguchi A, Rosenberg H. **Malignant hyperthermia.** Best Pract Res
    Clin Anaesthesiol 2003;17:519-533.
    <https://pubmed.ncbi.nlm.nih.gov/14661655/>

---

## 7. Metabolic, thermal and acid-base physiology (the two stores)

61. Gronert GA, Milde JH, Theye RA. **Role of sympathetic activity in porcine
    malignant hyperthermia.** Anesthesiology 1977;47:411-415.
    <https://pubmed.ncbi.nlm.nih.gov/911057/>
62. Gronert GA, Theye RA. **Halothane-induced porcine malignant hyperthermia:
    metabolic and hemodynamic changes.** Anesthesiology 1976;44:36-43.
    <https://pubmed.ncbi.nlm.nih.gov/1244758/>
    **The quantitative anchor for the model's whole-body VO₂ and VCO₂
    trajectories** — including the observation that CO₂ production rises
    proportionally more than O₂ consumption (the model's RQ > 1 result).
63. Gronert GA, Ahern CP, Milde JH, White RD. **Effect of CO2, calcium,
    digoxin, and potassium on cardiac and skeletal muscle metabolism in
    malignant hyperthermia susceptible swine.** Anesthesiology 1986;64:24-28.
    <https://pubmed.ncbi.nlm.nih.gov/3942334/>
64. Sessler DI. **Perioperative thermoregulation and heat balance.** Lancet
    2016;387:2655-2664. <https://pubmed.ncbi.nlm.nih.gov/26775126/>
    Source for the model's two-compartment core/shell thermal structure, the
    3.47 kJ/(kg·°C) specific heat, and the core-to-shell conductance.
65. Sessler DI. **Temperature monitoring and perioperative thermoregulation.**
    Anesthesiology 2008;109:318-338.
    <https://pubmed.ncbi.nlm.nih.gov/18648241/>
    Basis for the model's statement that skin and axillary sites lag the core
    by tens of minutes.
66. Nunn JF. **Nunn's Applied Respiratory Physiology** — the alveolar equation
    PaCO₂ = 863·V̇CO₂/V̇A and the body CO₂ stores used here. Primary data on
    CO₂ storage capacity: Cherniack NS, Longobardo GS. **Oxygen and carbon
    dioxide gas stores of the body.** Physiol Rev 1970;50:196-243.
    <https://pubmed.ncbi.nlm.nih.gov/4908762/>
    **The source for `K_CO2B` = 60 mL/mmHg**, the number that makes the whole
    EtCO₂-precedes-fever argument work.
67. Farhi LE, Rahn H. **Dynamics of changes in carbon dioxide stores.**
    Anesthesiology 1960;21:604-614.
    <https://pubmed.ncbi.nlm.nih.gov/13698627/>
68. Kemp GJ, Meyerspeer M, Moser E. **Absolute quantification of
    phosphorus metabolite concentrations in human muscle in vivo by 31P MRS.**
    NMR Biomed 2007;20:555-565. <https://pubmed.ncbi.nlm.nih.gov/17628042/>
    Source for the ATP and phosphocreatine pool sizes and the
    creatine-kinase clamp behaviour the model implements.
69. Olgin J, Argov Z, Rosenberg H, Tuchler M, Chance B. **Non-invasive
    evaluation of malignant hyperthermia susceptibility with phosphorus nuclear
    magnetic resonance spectroscopy.** Anesthesiology 1988;68:507-513.
    <https://pubmed.ncbi.nlm.nih.gov/3354981/>
70. Bendahan D, Kozak-Ribbens G, Confort-Gouny S, et al. **A noninvasive
    investigation of muscle energetics supports similarities between exertional
    heat stroke and malignant hyperthermia.** Anesth Analg 2001;93:683-689.
    <https://pubmed.ncbi.nlm.nih.gov/11524340/>

---

## 8. Dantrolene — mechanism, pharmacokinetics, formulation

71. Harrison GG. **Control of the malignant hyperpyrexic syndrome in MHS swine
    by dantrolene sodium.** Br J Anaesth 1975;47:62-65.
    <https://pubmed.ncbi.nlm.nih.gov/1148076/>
    The experiment that changed the mortality of this disease.
72. Kolb ME, Horne ML, Martz R. **Dantrolene in human malignant hyperthermia.**
    Anesthesiology 1982;56:254-262.
    <https://pubmed.ncbi.nlm.nih.gov/7039419/>
73. Flewellen EH, Nelson TE, Jones WP, Arens JF, Wagner DL. **Dantrolene dose
    response in awake man: implications for management of malignant
    hyperthermia.** Anesthesiology 1983;59:275-280.
    <https://pubmed.ncbi.nlm.nih.gov/6614536/>
    **The source for the model's dantrolene PK/PD**: 2.4 mg/kg gives ~4 µg/mL,
    2.8-4.2 µg/mL produces 75% twitch depression, effect persists for hours.
74. Paul-Pletzer K, Yamamoto T, Bhat MB, et al. **Identification of a dantrolene
    binding sequence on the skeletal muscle ryanodine receptor.** J Biol Chem
    2002;277:34918-34923. <https://pubmed.ncbi.nlm.nih.gov/12167662/>
    Residues 590-609 — the binding site drawn in cluster 4.
75. Fruen BR, Mickelson JR, Louis CF. **Dantrolene inhibition of sarcoplasmic
    reticulum Ca2+ release by direct and specific action at skeletal muscle
    ryanodine receptors.** J Biol Chem 1997;272:26965-26971.
    <https://pubmed.ncbi.nlm.nih.gov/9341133/>
    Establishes that the block is Mg²⁺-dependent and **partial** — the basis
    for `EMAX_DAN` = 0.87 rather than 1.0.
76. Choi RH, Koenig X, Launikonis BS. **Dantrolene requires Mg2+ to arrest
    malignant hyperthermia.** Proc Natl Acad Sci USA 2017;114:4811-4815.
    <https://pubmed.ncbi.nlm.nih.gov/28373559/>
77. Krause T, Gerbershagen MU, Fiege M, Weisshorn R, Wappler F. **Dantrolene —
    a review of its pharmacology, therapeutic use and new developments.**
    Anaesthesia 2004;59:364-373.
    <https://pubmed.ncbi.nlm.nih.gov/15023108/>
78. Schütte JK, Becker S, Burmester S, et al. **Comparison of the therapeutic
    effectiveness of a dantrolene sodium solution and a novel nanocrystalline
    suspension of dantrolene sodium in malignant hyperthermia normal and
    susceptible pigs.** Eur J Anaesthesiol 2011;28:256-264.
    <https://pubmed.ncbi.nlm.nih.gov/21206276/>
    Basis for treating Ryanodex vs Dantrium as **the same molecule with a
    different preparation delay**.
79. Sudo RT, Carmo PL, Trachez MM, Zapata-Sudo G. **Effects of azumolene on
    normal and malignant hyperthermia-susceptible skeletal muscle.**
    Basic Clin Pharmacol Toxicol 2008;102:308-316.
    <https://pubmed.ncbi.nlm.nih.gov/18047478/>
80. Brandom BW, Larach MG, Chen MS, Young MC. **Complications associated with
    the administration of dantrolene 1987 to 2006: a report from the North
    American Malignant Hyperthermia Registry.** Anesth Analg 2011;112:1115-1123.
    <https://pubmed.ncbi.nlm.nih.gov/21372281/>
81. Saltzman LS, Kates RA, Corke BC, Norfleet EA, Heath KR. **Hyperkalemia and
    cardiovascular collapse after verapamil and dantrolene administration in
    swine.** Anesth Analg 1984;63:473-478.
    <https://pubmed.ncbi.nlm.nih.gov/6711775/>
    The calcium-channel-blocker interaction flagged in cluster 14.

---

## 9. Treatment protocol, delay, and outcome

82. Glahn KPE, Ellis FR, Halsall PJ, et al. **Recognizing and managing a
    malignant hyperthermia crisis: guidelines from the European Malignant
    Hyperthermia Group.** Br J Anaesth 2010;105:417-420.
    <https://pubmed.ncbi.nlm.nih.gov/20852280/>
83. Hopkins PM, Girard T, Dalay S, et al. **Malignant hyperthermia 2020:
    guideline from the Association of Anaesthetists.** Anaesthesia
    2021;76:655-664. <https://pubmed.ncbi.nlm.nih.gov/33170980/>
84. Malignant Hyperthermia Association of the United States (MHAUS).
    **Managing an MH crisis.** <https://www.mhaus.org/healthcare-professionals/managing-a-crisis/>
    Source for the cooling protocol (cold saline, stop at 38 °C) whose
    arithmetic the model checks.
85. Burkman JM, Posner KL, Domino KB. **Analysis of the clinical variables
    associated with recrudescence after malignant hyperthermia reactions.**
    Anesthesiology 2007;106:901-906.
    <https://pubmed.ncbi.nlm.nih.gov/17457119/>
    The ~20% recrudescence rate and median ~13 h latency that the model
    **fails** to reproduce — see the limitations section of the README.
86. Larach MG, Localio AR, Allen GC, et al. **A clinical grading scale to
    predict malignant hyperthermia susceptibility.** Anesthesiology
    1994;80:771-779. <https://pubmed.ncbi.nlm.nih.gov/8024130/>
    The CGS implemented as a model *output* in `$TABLE`.
87. Riazi S, Larach MG, Hu C, Wijeysundera D, Massey C, Kraeva N. **Malignant
    hyperthermia in Canada: characteristics of index anesthetics in 129
    malignant hyperthermia susceptible probands.** Anesth Analg
    2014;118:381-387. <https://pubmed.ncbi.nlm.nih.gov/24361846/>
88. Litman RS, Griggs SM, Dowling JJ, Riazi S. **Malignant hyperthermia
    susceptibility and related diseases.** Anesthesiology 2018;128:159-167.
    <https://pubmed.ncbi.nlm.nih.gov/28902673/>
89. Larach MG, Dirksen SJH, Belani KG, et al. **Special article: Creation of a
    guide for the transfer of care of the malignant hyperthermia patient from
    ambulatory surgery centers to receiving hospital facilities.** Anesth Analg
    2012;114:94-100. <https://pubmed.ncbi.nlm.nih.gov/22052979/>
90. Nelson P, Litman RS. **Malignant hyperthermia in children: an analysis of
    the North American malignant hyperthermia registry.** Anesth Analg
    2014;118:369-374. <https://pubmed.ncbi.nlm.nih.gov/24356165/>

---

## 10. Diagnosis — IVCT / CHCT and genetic testing

91. Larach MG. **Standardization of the caffeine halothane muscle contracture
    test. North American Malignant Hyperthermia Group.** Anesth Analg
    1989;69:511-515. <https://pubmed.ncbi.nlm.nih.gov/2675676/>
92. Ording H, Brancadoro V, Cozzolino S, et al. **In vitro contracture test for
    diagnosis of malignant hyperthermia following the protocol of the European
    MH Group: results of testing patients surviving fulminant MH and unrelated
    low-risk subjects.** Acta Anaesthesiol Scand 1997;41:955-966.
    <https://pubmed.ncbi.nlm.nih.gov/9311391/>
    **The thresholds the model's simulated contracture test is scored
    against**: ≥0.2 g at ≤2% halothane or ≤2 mM caffeine.
93. Allen GC, Larach MG, Kunselman AR. **The sensitivity and specificity of the
    caffeine-halothane contracture test: a report from the North American
    Malignant Hyperthermia Registry.** Anesthesiology 1998;88:579-588.
    <https://pubmed.ncbi.nlm.nih.gov/9523799/>
    Sensitivity ~97-99%, specificity ~78-94%.
94. Hopkins PM, Rüffert H, Snoeck MM, et al. **European Malignant Hyperthermia
    Group guidelines for investigation of malignant hyperthermia
    susceptibility.** Br J Anaesth 2015;115:531-539.
    <https://pubmed.ncbi.nlm.nih.gov/26385664/>
95. Urwyler A, Deufel T, McCarthy T, West S. **Guidelines for molecular genetic
    detection of susceptibility to malignant hyperthermia.** Br J Anaesth
    2001;86:283-287. <https://pubmed.ncbi.nlm.nih.gov/11573680/>
96. Anderson-Pompa K, Foster A, Parker L, Wilks C, Cheatham AP, Gonzalez J.
    **Genetics and susceptibility to malignant hyperthermia.** Crit Care Nurse
    2008;28:32-38. <https://pubmed.ncbi.nlm.nih.gov/18827086/>

---

## 11. Complications — rhabdomyolysis, hyperkalaemia, renal injury, DIC

97. Bosch X, Poch E, Grau JM. **Rhabdomyolysis and acute kidney injury.**
    N Engl J Med 2009;361:62-72. <https://pubmed.ncbi.nlm.nih.gov/19571284/>
    Source for the model's myoglobin-cast / GFR limb and the CK kinetics.
98. Chavez LO, Leon M, Einav S, Varon J. **Beyond muscle destruction: a
    systematic review of rhabdomyolysis for clinical practice.** Crit Care
    2016;20:135. <https://pubmed.ncbi.nlm.nih.gov/27301374/>
99. Larach MG, Brandom BW, Allen GC, Gronert GA, Lehman EB. **Cardiac arrests
    and deaths associated with malignant hyperthermia in North America from
    1987 to 2006.** Anesthesiology 2008;108:603-611.
    <https://pubmed.ncbi.nlm.nih.gov/18362591/>
    Source for the model's hyperkalaemic-arrest and DIC limbs, and for the
    finding that DIC is among the strongest predictors of death.
100. Bunn HF, Jandl JH — classic work on myoglobin clearance; for a clinical
     synthesis of myoglobin's short half-life relative to CK see reference 97
     and: Huerta-Alardín AL, Varon J, Marik PE. **Bench-to-bedside review:
     rhabdomyolysis — an overview for clinicians.** Crit Care 2005;9:158-169.
     <https://pubmed.ncbi.nlm.nih.gov/15774072/>
101. Weingarten TN, Ackerman J, Sprung J. **Compartment syndrome as a
     complication of malignant hyperthermia.** — see the broader review:
     Poole TC, Lim TY, Buck J, Kong AS. **Perioperative cardiac arrest in a
     patient with previously undiagnosed Becker's muscular dystrophy after
     isoflurane anaesthesia for elective surgery.** Br J Anaesth
     2010;104:487-489. <https://pubmed.ncbi.nlm.nih.gov/20211990/>
102. Larach MG, Brandom BW, Allen GC, Gronert GA, Lehman EB. **Malignant
     hyperthermia deaths related to inadequate temperature monitoring,
     2007-2012: a report from the North American Malignant Hyperthermia
     Registry.** Anesth Analg 2014;119:1359-1366.
     <https://pubmed.ncbi.nlm.nih.gov/25268394/>

---

## 12. Related and differential syndromes

103. Capacchione JF, Muldoon SM. **The relationship between exertional heat
     illness, exertional rhabdomyolysis, and malignant hyperthermia.**
     Anesth Analg 2009;109:1065-1069.
     <https://pubmed.ncbi.nlm.nih.gov/19617585/>
104. Muldoon S, Deuster P, Brandom B, Bunger R. **Is there a link between
     malignant hyperthermia and exertional heat illness?** Exerc Sport Sci Rev
     2004;32:174-179. <https://pubmed.ncbi.nlm.nih.gov/15604937/>
105. Gurrera RJ. **Sympathoadrenal hyperactivity and the etiology of neuroleptic
     malignant syndrome.** Am J Psychiatry 1999;156:169-180.
     <https://pubmed.ncbi.nlm.nih.gov/9989551/>
106. Boyer EW, Shannon M. **The serotonin syndrome.** N Engl J Med
     2005;352:1112-1120. <https://pubmed.ncbi.nlm.nih.gov/15784664/>
107. Klingler W, Rueffert H, Lehmann-Horn F, Girard T, Hopkins PM. **Core
     myopathies and risk of malignant hyperthermia.** Anesth Analg
     2009;109:1167-1173. <https://pubmed.ncbi.nlm.nih.gov/19762744/>
108. Jungbluth H, Treves S, Zorzato F, et al. **Congenital myopathies: disorders
     of excitation-contraction coupling and muscle contraction.** Nat Rev Neurol
     2018;14:151-167. <https://pubmed.ncbi.nlm.nih.gov/29391587/>

---

## 13. QSP modelling methodology

109. Peterson MC, Riggs MM. **FDA Advisory Meeting clinical pharmacology review
     utilizes a quantitative systems pharmacology (QSP) model.** CPT
     Pharmacometrics Syst Pharmacol 2015;4:189-192.
     <https://pubmed.ncbi.nlm.nih.gov/26225239/>
110. Nijsen MJMA, Wu F, Bansal L, et al. **Preclinical QSP modeling in the
     pharmaceutical industry: an IQ consortium survey examining the current
     landscape.** CPT Pharmacometrics Syst Pharmacol 2018;7:135-146.
     <https://pubmed.ncbi.nlm.nih.gov/29349875/>
111. Baron KT, Gastonguay MR. **mrgsolve: simulate from ODE-based population
     PK/PD and systems pharmacology models.** R package.
     <https://mrgsolve.org/>
112. Musante CJ, Ramanujan S, Schmidt BJ, Ghobrial OG, Lu J, Heatherington AC.
     **Quantitative systems pharmacology: a case for disease models.**
     Clin Pharmacol Ther 2017;101:24-27.
     <https://pubmed.ncbi.nlm.nih.gov/27709613/>

---

## Emergency information

- **MHAUS 24-hour hotline (North America):** 1-800-644-9737
  (outside the US: +1-209-417-3722) · <https://www.mhaus.org/>
- **European Malignant Hyperthermia Group:** <https://www.emhg.org/>
- **UK MH Investigation Unit, Leeds:** <https://www.leedsth.nhs.uk/>

Nothing in this directory is a clinical decision aid. In a suspected MH crisis,
follow the MHAUS or EMHG protocol and telephone the hotline.
