# Familial Chylomicronemia Syndrome (FCS) — Reference Library

Curated evidence base for the QSP model in this directory
(`fcs_qsp_model.dot` · `fcs_mrgsolve_model.R` · `fcs_shiny_app.R`).
References are grouped by the model block they support, so that every
structural choice and every calibration target can be traced to a source.

**48 references.** Links are PubMed **title queries**
(`pubmed.ncbi.nlm.nih.gov/?term=...`) rather than hard-coded PMIDs, so that
each link resolves to the indexed record even if an identifier is mistyped.
Guidelines and regulatory documents link to the issuing body's page.

---

## 1. Definition, epidemiology and the FCS-versus-MCS distinction

The first structural decision in this model is that FCS and multifactorial
chylomicronemia (MCS) are *different diseases with the same laboratory value*.
FCS is a monogenic loss of the lipolysis machinery; MCS is polygenic risk plus
a secondary factor acting on a machinery that still works. The model expresses
that difference as a single parameter, `F_GENO`, and everything else follows.

1. Brahm AJ, Hegele RA. **Chylomicronaemia — current diagnosis and future
   therapies.** *Nat Rev Endocrinol* 2015;11(6):352-62.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Chylomicronaemia+current+diagnosis+and+future+therapies>
2. Hegele RA, Berberich AJ, Ban MR, et al. **Clinical and biochemical features
   of different molecular etiologies of familial chylomicronemia.**
   *J Clin Lipidol* 2018;12(4):920-927.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+and+biochemical+features+of+different+molecular+etiologies+of+familial+chylomicronemia>
3. Moulin P, Dufour R, Averna M, et al. **Identification and diagnosis of
   patients with familial chylomicronaemia syndrome (FCS): expert panel
   recommendations and proposal of an "FCS score".**
   *Atherosclerosis* 2018;275:265-272.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Identification+and+diagnosis+of+patients+with+familial+chylomicronaemia+syndrome+FCS+score>
4. Paquette M, Bernard S, Hegele RA, Baass A. **Chylomicronemia: differences
   between familial chylomicronemia syndrome and multifactorial
   chylomicronemia.** *Atherosclerosis* 2019;283:137-142.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Chylomicronemia+differences+between+familial+chylomicronemia+syndrome+and+multifactorial+chylomicronemia>
5. Chait A, Eckel RH. **The chylomicronemia syndrome is most often multifactorial:
   a narrative review of causes and treatment.** *Ann Intern Med*
   2019;170(9):626-634.
   <https://pubmed.ncbi.nlm.nih.gov/?term=The+chylomicronemia+syndrome+is+most+often+multifactorial+narrative+review>
6. Stroes ESG, Moulin P, Parhofer KG, et al. **Diagnostic algorithm for familial
   chylomicronemia syndrome.** *Atheroscler Suppl* 2017;23:1-7.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Diagnostic+algorithm+for+familial+chylomicronemia+syndrome>
7. Falko JM. **Familial chylomicronemia syndrome: a clinical guide for
   endocrinologists.** *Endocr Pract* 2018;24(8):756-763.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Familial+chylomicronemia+syndrome+a+clinical+guide+for+endocrinologists>

## 2. Genetic architecture — five genes, one enzyme complex

Supports cluster 1 of the map and the `geno` presets in the R model. The single
most important quantitative statement in this section is that phenotype tracks
**residual post-heparin lipase activity**, and collapses only below roughly
5-10% of normal — which is why `FCS_genotype_gradient()` exists.

8. Rahalkar AR, Giffen F, Har B, et al. **Novel LPL mutations associated with
   lipoprotein lipase deficiency: two case reports and a literature review.**
   *Can J Physiol Pharmacol* 2009;87(3):151-60.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Novel+LPL+mutations+associated+with+lipoprotein+lipase+deficiency+two+case+reports+and+a+literature+review>
9. Beigneux AP, Davies BSJ, Gin P, et al. **Glycosylphosphatidylinositol-anchored
   high-density lipoprotein-binding protein 1 plays a critical role in the
   lipolytic processing of chylomicrons.** *Cell Metab* 2007;5(4):279-91.
   <https://pubmed.ncbi.nlm.nih.gov/?term=GPIHBP1+plays+a+critical+role+in+the+lipolytic+processing+of+chylomicrons>
10. Beigneux AP, Miyashita K, Ploug M, et al. **Autoantibodies against GPIHBP1
    as a cause of hypertriglyceridemia.** *N Engl J Med* 2017;376(17):1647-1658.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Autoantibodies+against+GPIHBP1+as+a+cause+of+hypertriglyceridemia>
11. Péterfy M, Ben-Zeev O, Mao HZ, et al. **Mutations in LMF1 cause combined
    lipase deficiency and severe hypertriglyceridemia.** *Nat Genet*
    2007;39(12):1483-7.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Mutations+in+LMF1+cause+combined+lipase+deficiency+and+severe+hypertriglyceridemia>
12. Priore Oliva C, Pisciotta L, Li Volti G, et al. **Inherited apolipoprotein
    A-V deficiency in severe hypertriglyceridemia.** *Arterioscler Thromb Vasc
    Biol* 2005;25(2):411-7.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Inherited+apolipoprotein+A-V+deficiency+in+severe+hypertriglyceridemia>
13. Dron JS, Hegele RA. **Genetics of hypertriglyceridemia.** *Front Endocrinol*
    2020;11:455.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Genetics+of+hypertriglyceridemia+Dron+Hegele+Frontiers+in+Endocrinology>
14. Dron JS, Wang J, Cao H, et al. **Severe hypertriglyceridemia is primarily
    polygenic.** *J Clin Lipidol* 2019;13(1):80-88.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Severe+hypertriglyceridemia+is+primarily+polygenic>

## 3. The LPL machinery — limb 1 of the clearance sum

Supports cluster 5. The map draws LPL as an enzyme that is *made* in one cell,
*folded* by LMF1, *moved* by GPIHBP1 and *switched on* by apoC-II, because
losing any one of those steps produces the same phenotype — the structural
justification for a single multiplicative genotype factor in the ODEs.

15. Young SG, Fong LG, Beigneux AP, et al. **GPIHBP1 and lipoprotein lipase,
    partners in plasma triglyceride metabolism.** *Cell Metab* 2019;30(1):51-65.
    <https://pubmed.ncbi.nlm.nih.gov/?term=GPIHBP1+and+lipoprotein+lipase+partners+in+plasma+triglyceride+metabolism>
16. Birrane G, Beigneux AP, Dwyer B, et al. **Structure of the lipoprotein
    lipase-GPIHBP1 complex that mediates plasma triglyceride hydrolysis.**
    *Proc Natl Acad Sci U S A* 2019;116(5):1723-1732.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Structure+of+the+lipoprotein+lipase-GPIHBP1+complex+that+mediates+plasma+triglyceride+hydrolysis>
17. Kristensen KK, Leth-Espensen KZ, Mertens HDT, et al. **Unfolding of
    monomeric lipoprotein lipase by ANGPTL4: insight into the regulation of
    plasma triglyceride metabolism.** *Proc Natl Acad Sci U S A*
    2020;117(8):4337-4346.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Unfolding+of+monomeric+lipoprotein+lipase+by+ANGPTL4>
18. Wolska A, Dunbar RL, Freeman LA, et al. **Apolipoprotein C-II: new findings
    related to genetics, biochemistry, and role in triglyceride metabolism.**
    *Atherosclerosis* 2017;267:49-60.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Apolipoprotein+C-II+new+findings+related+to+genetics+biochemistry+and+role+in+triglyceride+metabolism>

## 4. ANGPTL3 / 4 / 8 — the physiological brake on limb 1

Supports cluster 16 and the model's central falsifiable prediction: an ANGPTL3
drug is a **limb-1** drug and therefore cannot work when limb 1 is zero.

19. Chi X, Britt EC, Shows HW, et al. **ANGPTL8 promotes the ability of ANGPTL3
    to inhibit lipoprotein lipase.** *Mol Metab* 2017;6(10):1137-1149.
    <https://pubmed.ncbi.nlm.nih.gov/?term=ANGPTL8+promotes+the+ability+of+ANGPTL3+to+inhibit+lipoprotein+lipase>
20. Musunuru K, Pirruccello JP, Do R, et al. **Exome sequencing, ANGPTL3
    mutations, and familial combined hypolipidemia.** *N Engl J Med*
    2010;363(23):2220-7.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Exome+sequencing+ANGPTL3+mutations+and+familial+combined+hypolipidemia>
21. Ahmad Z, Banerjee P, Hamon S, et al. **Inhibition of angiopoietin-like
    protein 3 with a monoclonal antibody reduces triglycerides in
    hypertriglyceridemia.** *Circulation* 2019;140(6):470-486.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Inhibition+of+angiopoietin-like+protein+3+with+a+monoclonal+antibody+reduces+triglycerides+in+hypertriglyceridemia>
22. Rosenson RS, Gaudet D, Ballantyne CM, et al. **Evinacumab in severe
    hypertriglyceridemia with or without lipoprotein lipase pathway
    mutations: a phase 2 randomized trial.** *Nat Med* 2023;29(3):729-737.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Evinacumab+in+severe+hypertriglyceridemia+with+or+without+lipoprotein+lipase+pathway+mutations>

> Reference 22 is the closest thing to a direct test of the model's second
> claim: in patients **with** biallelic LPL-pathway mutations the ANGPTL3
> antibody produced little triglyceride lowering, while patients **without**
> them responded. That is limb arithmetic, observed.

## 5. apoC-III — the node that sits on both limbs

Supports clusters 6, 7 and the `GC3` amplifier in the ODEs. The decisive
mechanistic point for FCS is that apoC-III inhibits **hepatic remnant uptake**
independently of lipolysis, which is why knocking it down works in a patient
with no lipase at all.

23. Gordts PLSM, Nock R, Son NH, et al. **ApoC-III inhibits clearance of
    triglyceride-rich lipoproteins through LDL family receptors.**
    *J Clin Invest* 2016;126(8):2855-66.
    <https://pubmed.ncbi.nlm.nih.gov/?term=ApoC-III+inhibits+clearance+of+triglyceride-rich+lipoproteins+through+LDL+family+receptors>
24. Ramms B, Patel S, Nora C, et al. **ApoC-III ASO promotes tissue LPL-independent
    clearance of triglyceride-rich lipoproteins in humans and mice.**
    *J Lipid Res* 2019;60(8):1379-1395.
    <https://pubmed.ncbi.nlm.nih.gov/?term=ApoC-III+ASO+promotes+tissue+LPL-independent+clearance+of+triglyceride-rich+lipoproteins>
25. Pollin TI, Damcott CM, Shen H, et al. **A null mutation in human APOC3
    confers a favorable plasma lipid profile and apparent cardioprotection.**
    *Science* 2008;322(5908):1702-5.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+null+mutation+in+human+APOC3+confers+a+favorable+plasma+lipid+profile+and+apparent+cardioprotection>
26. TG and HDL Working Group of the Exome Sequencing Project. **Loss-of-function
    mutations in APOC3, triglycerides, and coronary disease.** *N Engl J Med*
    2014;371(1):22-31.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Loss-of-function+mutations+in+APOC3+triglycerides+and+coronary+disease>
27. Norata GD, Tsimikas S, Pirillo A, Catapano AL. **Apolipoprotein C-III:
    from pathophysiology to pharmacology.** *Trends Pharmacol Sci*
    2015;36(10):675-687.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Apolipoprotein+C-III+from+pathophysiology+to+pharmacology>
28. Taskinen MR, Packard CJ, Borén J. **Emerging evidence that ApoC-III inhibitors
    provide novel options to reduce the residual CVD.** *Curr Atheroscler Rep*
    2019;21(8):27.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emerging+evidence+that+ApoC-III+inhibitors+provide+novel+options+to+reduce+the+residual+CVD>

## 6. Chylomicron assembly, postprandial lipemia and saturable clearance

Supports clusters 3, 4, 6 and the Michaelis-Menten structure of limb 2.

29. Hussain MM. **Intestinal lipid absorption and lipoprotein formation.**
    *Curr Opin Lipidol* 2014;25(3):200-6.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Intestinal+lipid+absorption+and+lipoprotein+formation+Hussain>
30. Brunzell JD, Hazzard WR, Porte D Jr, Bierman EL. **Evidence for a common,
    saturable, triglyceride removal mechanism for chylomicrons and very low
    density lipoproteins in man.** *J Clin Invest* 1973;52(7):1578-85.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Evidence+for+a+common+saturable+triglyceride+removal+mechanism+for+chylomicrons+and+very+low+density+lipoproteins+in+man>
31. Borén J, Taskinen MR, Björnson E, Packard CJ. **Metabolism of triglyceride-rich
    lipoproteins in health and dyslipidaemia.** *Nat Rev Cardiol*
    2022;19(9):577-592.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Metabolism+of+triglyceride-rich+lipoproteins+in+health+and+dyslipidaemia>

> Reference 30 is the empirical origin of this model's most important
> structural choice. Brunzell's 1973 demonstration that triglyceride removal is
> a **common and saturable** capacity shared by chylomicrons and VLDL is what
> makes the steady-state relationship between dietary fat and plasma TG a
> hyperbola rather than a line, and therefore what makes the FCS diet
> prescription a cliff edge rather than a dose-response.

## 7. Hypertriglyceridemic acute pancreatitis — mechanism and threshold

Supports cluster 10 and the convex hazard function. The model's Hill exponent
of 1.7 and the 880 mg/dL (10 mmol/L) threshold both come from this section.

32. Yang AL, McNabb-Baltar J. **Hypertriglyceridemia and acute pancreatitis.**
    *Pancreatology* 2020;20(5):795-800.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hypertriglyceridemia+and+acute+pancreatitis+Yang+McNabb-Baltar>
33. Pedersen SB, Langsted A, Nordestgaard BG. **Nonfasting mild-to-moderate
    hypertriglyceridemia and risk of acute pancreatitis.** *JAMA Intern Med*
    2016;176(12):1834-1842.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nonfasting+mild-to-moderate+hypertriglyceridemia+and+risk+of+acute+pancreatitis>
34. Saharia P, Margolis S, Zuidema GD, Cameron JL. **Acute pancreatitis with
    hyperlipemia: studies with an isolated perfused canine pancreas.**
    *Surgery* 1977;82(1):60-7.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Acute+pancreatitis+with+hyperlipemia+studies+with+an+isolated+perfused+canine+pancreas>
35. Criddle DN. **The role of fat and alcohol in acute pancreatitis: a dangerous
    liaison.** *Pancreatology* 2015;15(4 Suppl):S6-S12.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+role+of+fat+and+alcohol+in+acute+pancreatitis+a+dangerous+liaison>
36. Huang W, Booth DM, Cane MC, et al. **Fatty acid ethyl ester synthase
    inhibition ameliorates ethanol-induced Ca2+-dependent mitochondrial
    dysfunction and acute pancreatitis.** *Gut* 2014;63(8):1313-24.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fatty+acid+ethyl+ester+synthase+inhibition+ameliorates+ethanol-induced+Ca2+-dependent+mitochondrial+dysfunction+and+acute+pancreatitis>
37. Berglund L, Brunzell JD, Goldberg AC, et al. **Evaluation and treatment of
    hypertriglyceridemia: an Endocrine Society clinical practice guideline.**
    *J Clin Endocrinol Metab* 2012;97(9):2969-89.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Evaluation+and+treatment+of+hypertriglyceridemia+an+Endocrine+Society+clinical+practice+guideline>

## 8. Diet, MCT and non-pharmacological management

Supports cluster 3 and the `FMCT` parameter.

38. Williams L, Rhodes KS, Karmally W, et al. **Familial chylomicronemia
    syndrome: bringing to life dietary recommendations throughout the life
    span.** *J Clin Lipidol* 2018;12(4):908-919.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Familial+chylomicronemia+syndrome+bringing+to+life+dietary+recommendations+throughout+the+life+span>
39. Davidson M, Stevenson M, Hsieh A, et al. **The burden of familial
    chylomicronemia syndrome: results from the global IN-FOCUS study.**
    *J Clin Lipidol* 2018;12(4):898-907.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+burden+of+familial+chylomicronemia+syndrome+results+from+the+global+IN-FOCUS+study>
40. Gelrud A, Digenio A, Alexander V, Williams KR. **Treatment burden and
    disease impact in familial chylomicronemia syndrome: patient survey.**
    *J Clin Lipidol* 2017;11(3):814.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+burden+and+disease+impact+in+familial+chylomicronemia+syndrome+patient+survey>

## 9. Conventional lipid drugs and why they fail in true FCS

Supports cluster 12 and the `FCS_limb_decomposition()` output.

41. Gaudet D, Méthot J, Déry S, et al. **Efficacy and long-term safety of
    alipogene tiparvovec (AAV1-LPL S447X) gene therapy for lipoprotein lipase
    deficiency: an open-label trial.** *Gene Ther* 2013;20(4):361-9.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Efficacy+and+long-term+safety+of+alipogene+tiparvovec+AAV1-LPL+S447X+gene+therapy+for+lipoprotein+lipase+deficiency>
42. Gaudet D, Brisson D, Tremblay K, et al. **Targeting APOC3 in the familial
    chylomicronemia syndrome.** *N Engl J Med* 2014;371(23):2200-6.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Targeting+APOC3+in+the+familial+chylomicronemia+syndrome>

## 10. apoC-III antisense — volanesorsen

Supports cluster 13 and calibration targets D.

43. Witztum JL, Gaudet D, Freedman SD, et al. **Volanesorsen and triglyceride
    levels in familial chylomicronemia syndrome (APPROACH).** *N Engl J Med*
    2019;381(6):531-542.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Volanesorsen+and+triglyceride+levels+in+familial+chylomicronemia+syndrome>
44. Gouni-Berthold I, Alexander VJ, Yang Q, et al. **Efficacy and safety of
    volanesorsen in patients with multifactorial chylomicronaemia (COMPASS).**
    *Lancet Diabetes Endocrinol* 2021;9(5):264-275.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Efficacy+and+safety+of+volanesorsen+in+patients+with+multifactorial+chylomicronaemia+COMPASS>
45. European Medicines Agency. **Waylivra (volanesorsen): EPAR — product
    information.** Conditional marketing authorisation, 2019.
    <https://www.ema.europa.eu/en/medicines/human/EPAR/waylivra>

## 11. GalNAc-conjugated apoC-III therapeutics — olezarsen and plozasiran

Supports clusters 14, 15, 17 and calibration targets E and F. This is where
the "GalNAc dividend" of the model comes from: identical target, ~20-30-fold
lower systemic exposure, and the platelet signal disappears.

46. Stroes ESG, Alexander VJ, Karwatowska-Prokopczuk E, et al. **Olezarsen,
    acute pancreatitis, and familial chylomicronemia syndrome (Balance).**
    *N Engl J Med* 2024;390(19):1781-1792.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Olezarsen+acute+pancreatitis+and+familial+chylomicronemia+syndrome>
47. Watts GF, Rosenson RS, Hegele RA, et al. **Plozasiran for managing
    persistent chylomicronemia and pancreatitis risk (PALISADE).**
    *N Engl J Med* 2024/2025.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Plozasiran+for+managing+persistent+chylomicronemia+and+pancreatitis+risk>
48. Gaudet D, Pall D, Watts GF, et al. **Plozasiran (ARO-APOC3) for severe
    hypertriglyceridemia: the SHASTA-2 randomized clinical trial.**
    *JAMA Cardiol* 2024;9(7):620-630.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Plozasiran+ARO-APOC3+for+severe+hypertriglyceridemia+the+SHASTA-2+randomized+clinical+trial>

---

## How each reference enters the model

| Model element | Parameter(s) | Anchoring references |
|---|---|---|
| Two-limb clearance sum | `CL_LPL_MAX`, `VMAX_IND`, `KM_IND` | 30, 31, 24 |
| Genotype factor and its cliff | `F_GENO`, `A_C2`, `A_A5` | 2, 8, 13, 15 |
| GPIHBP1 / LMF1 / apoA-V as one machine | map cluster 5 | 9, 10, 11, 12, 16, 18 |
| ANGPTL3 brake, evinacumab prediction | `IC50_ANG`, `IMAX_EVI` | 17, 19, 20, 21, 22 |
| apoC-III raises limb-2 Vmax | `IMAX_C3`, `GC3` | 23, 24, 25, 26, 27 |
| Chylomicron assembly and lymph delay | `KA_FAT`, `K_ENT`, `K_LYM` | 29, 31 |
| Convex AP hazard and 880 mg/dL threshold | `LAM_MAX`, `TG50_AP`, `HILL_AP` | 32, 33, 37 |
| Unbound-FFA acinar injury mechanism | `K_PFFA`, `PFFA_THR`, `K_INJ` | 34, 35, 36 |
| Diet as the input term; MCT bypass | `F_ABS`, `FMCT` | 38 |
| Disease burden, PROs, brain fog | `KF_IN`, `KX_IN` | 39, 40 |
| Fibrate/omega-3 futility in LPL-null FCS | `EFIB_LPL`, `IFIB_C3`, `I_OM3` | 5, 7, 42 |
| Volanesorsen PK/PD and thrombocytopenia | `EC50_VOL`, `IPLT_S`, `KDES_PLT` | 42, 43, 44, 45 |
| Olezarsen GalNAc PK/PD | `EC50_OLE`, `KL_OLE` | 46 |
| Plozasiran siRNA and RISC persistence | `EC50_SI`, `KRISC_IN/OUT` | 47, 48 |
| Gene therapy as the only limb-1 repair | map cluster 16 | 41 |

---

## A note on what the model does *not* claim

The mechanistic map draws atherosclerotic risk in FCS as **not clearly
increased**, on the grounds that chylomicrons are too large to enter the
arterial intima — in deliberate contrast to remnant-rich multifactorial
chylomicronemia, where remnant cholesterol is atherogenic. That contrast is
still debated, and the model contains no cardiovascular endpoint. Anyone
extending it to cardiovascular outcomes should treat references 25, 26 and 31
as the starting point and should not assume that the FCS and MCS populations
share a risk function, because in this model they do not even share a disease.

---

*Research and education only. Not for clinical decision-making.*
