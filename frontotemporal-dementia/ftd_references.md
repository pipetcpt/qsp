# Frontotemporal Dementia (FTD / FTLD) — QSP Model References

Supporting literature for `ftd_qsp_model.dot`, `ftd_mrgsolve_model.R`,
`ftd_reference_model.py` and `ftd_shiny_app.R`.

## How to read this list — PMID verification status

Two classes of entry, kept deliberately separate:

- **✅ VERIFIED** — the PMID was resolved against the NCBI E-utilities API during
  the session that produced this model, and the returned title / journal /
  year are reproduced below as returned. These you can trust as identifiers.
- **🔍 SEARCH** — no PMID is asserted. A PubMed search URL is given instead.
  These are topics where a specific paper is cited for context but the exact
  identifier was not machine-verified here, so quoting a PMID would be a guess
  dressed up as a fact.

Values used to calibrate the model are marked **[ANCHOR]** and the specific
number taken from the paper is stated, so any disagreement is auditable.

---

## 1. Genetics — the three major genes

1. ✅ **Baker M, et al. Mutations in progranulin cause tau-negative
   frontotemporal dementia linked to chromosome 17.** *Nature* 2006 Aug 24.
   PMID **16862116** — <https://pubmed.ncbi.nlm.nih.gov/16862116/>
2. ✅ **Cruts M, et al. Null mutations in progranulin cause ubiquitin-positive
   frontotemporal dementia linked to chromosome 17q21.** *Nature* 2006 Aug 24.
   PMID **16862115** — <https://pubmed.ncbi.nlm.nih.gov/16862115/>
3. ✅ **DeJesus-Hernandez M, et al. Expanded GGGGCC hexanucleotide repeat in
   noncoding region of C9ORF72 causes chromosome 9p-linked FTD and ALS.**
   *Neuron* 2011 Oct 20. PMID **21944778** —
   <https://pubmed.ncbi.nlm.nih.gov/21944778/>
4. ✅ **Renton AE, et al. A hexanucleotide repeat expansion in C9ORF72 is the
   cause of chromosome 9p21-linked ALS-FTD.** *Neuron* 2011 Oct 20.
   PMID **21944779** — <https://pubmed.ncbi.nlm.nih.gov/21944779/>
5. ✅ **Hutton M, et al. Association of missense and 5'-splice-site mutations in
   tau with the inherited dementia FTDP-17.** *Nature* 1998 Jun 18.
   PMID **9641683** — <https://pubmed.ncbi.nlm.nih.gov/9641683/>
   Basis for `MUT_TAU_BOOST` and the exon-10 / 4R:3R node in the map.
6. ✅ **Van Deerlin VM, et al. Common variants at 7p21 are associated with
   frontotemporal lobar degeneration with TDP-43 inclusions.** *Nat Genet*
   2010 Mar. PMID **20154673** — <https://pubmed.ncbi.nlm.nih.gov/20154673/>
   The TMEM106B locus; basis for `TMEM106B_PROT` / `TMEM_PEN`.
7. ✅ **Pottier C, et al. Potential genetic modifiers of disease risk and age at
   onset in patients with frontotemporal lobar degeneration and C9orf72
   repeat expansions.** *Lancet Neurol* 2018 Jun. PMID **29724592** —
   <https://pubmed.ncbi.nlm.nih.gov/29724592/>
8. ✅ **Seelaar H, et al. Clinical, genetic and pathological heterogeneity of
   frontotemporal dementia: a review.** *J Neurol Neurosurg Psychiatry*
   2011 May. PMID **20971753** — <https://pubmed.ncbi.nlm.nih.gov/20971753/>
9. 🔍 TBK1 / VCP / CHMP2B / SQSTM1 / OPTN / CHCHD10 in FTD-ALS —
   <https://pubmed.ncbi.nlm.nih.gov/?term=TBK1+VCP+CHMP2B+frontotemporal+dementia+ALS+genetics>

## 2. Onset age, duration and survival — **[ANCHOR]**

10. ✅ **[ANCHOR] Moore KM, et al. Age at symptom onset and death and disease
    duration in genetic frontotemporal dementia: an international
    retrospective cohort study.** *Lancet Neurol* 2020 Feb.
    PMID **31810826** — <https://pubmed.ncbi.nlm.nih.gov/31810826/>
    **The single most important calibration target in this model.** From 3403
    individuals in 1492 families (1433 C9orf72, 1179 GRN, 791 MAPT), mean age
    at symptom onset / at death / disease duration:
    - **MAPT 49.5 y (SD 10.0) / 58.5 y / 9.3 y**
    - **C9orf72 58.2 y (SD 9.8) / 65.3 y / 6.4 y**
    - **GRN 61.3 y (SD 8.8) / 68.8 y / 7.1 y**

    The model predicts 50.5 / 54.3 / 55.8 y — MAPT within 1 y, C9orf72 ~4 y
    early, GRN ~5.5 y early. The GRN miss is reported as a *structural* failure
    in section A of `ftd_model_report.txt`: GRN onset is measured to be LATER
    than sporadic bvFTD (~58 y), and no model in which haploinsufficiency only
    ever adds damage can reproduce that ordering.
11. ✅ **Rohrer JD, et al. Presymptomatic cognitive and neuroanatomical changes
    in genetic frontotemporal dementia in the Genetic Frontotemporal dementia
    Initiative (GENFI) study.** *Lancet Neurol* 2015 Mar. PMID **25662776** —
    <https://pubmed.ncbi.nlm.nih.gov/25662776/>

## 3. Progranulin biology and the sortilin axis — the drug-target block

12. ✅ **[ANCHOR] Hu F, et al. Sortilin-mediated endocytosis determines levels of
    the frontotemporal dementia protein, progranulin.** *Neuron* 2010 Nov 18.
    PMID **21092856** — <https://pubmed.ncbi.nlm.nih.gov/21092856/>
    The foundational result the whole model turns on: sortilin is the
    high-affinity receptor that *sets* extracellular PGRN concentration by
    rapid endocytic clearance. Basis for `KCL_SORT_PL` / `KCL_SORT_CSF`.
13. ✅ **Zhou X, et al. Prosaposin facilitates sortilin-independent lysosomal
    trafficking of progranulin.** *J Cell Biol* 2015 Sep 14. PMID **26370502**
    — <https://pubmed.ncbi.nlm.nih.gov/26370502/>
    **This paper is why the model can say anything at all.** It establishes a
    *sortilin-independent* route into the lysosome, which is exactly the
    `1 - P_LYS_SORT` term. If lysosomal delivery were 100% sortilin-dependent
    (p = 1), anti-sortilin therapy would necessarily *reduce* the lysosomal
    pool. The existence of the prosaposin route is what makes benefit possible
    — and the unmeasured *split* between the two routes is what decides the
    sign of the effect.
14. ✅ **Paushter DH, et al. The lysosomal function of progranulin, a guardian
    against neurodegeneration.** *Acta Neuropathol* 2018 Jul. PMID **29744576**
    — <https://pubmed.ncbi.nlm.nih.gov/29744576/>
    Basis for `W_ACT_LYS` (the lysosomal-action hypothesis) and the
    PGRN → cathepsin D / GCase / BMP arm of the map.
15. ✅ **Valdez C, et al. Progranulin-mediated deficiency of cathepsin D results
    in FTD and NCL-like phenotypes in neurons derived from FTD patients.**
    *Hum Mol Genet* 2017 Dec 15. PMID **29036611** —
    <https://pubmed.ncbi.nlm.nih.gov/29036611/>
16. ✅ **Ward ME, et al. Individuals with progranulin haploinsufficiency exhibit
    features of neuronal ceroid lipofuscinosis.** *Sci Transl Med* 2017 Apr 12.
    PMID **28404863** — <https://pubmed.ncbi.nlm.nih.gov/28404863/>
    Human evidence that *heterozygous* PGRN loss already produces a lysosomal
    storage phenotype — the dose-dependence argument for `LIPO` / `K_LIPO_DAM`.
17. ✅ **[ANCHOR] Swift IJ, et al. A systematic review of progranulin
    concentrations in biofluids in over 7,000 people.** *Alzheimers Res Ther*
    2024 Mar 28. PMID **38539243** — <https://pubmed.ncbi.nlm.nih.gov/38539243/>
    Pooled data from 7071 individuals across 75 publications.
    **Plasma PGRN cut-off between pathogenic GRN mutation carriers and
    non-carriers 74.8 ng/mL (Adipogen assay); CSF cut-off 3.43 ng/mL.**
    The model's carrier plasma value (~70 ng/mL against a 200 ng/mL control) sits
    at that cut-off. Also reports **no significant association between plasma
    PGRN and TMEM106B rs1990622 genotype** — which the model reproduces without
    being fitted to it, because TMEM106B is wired to the lysosomal/age
    vulnerability axis and not to PGRN levels.
18. ✅ **Lui H, et al. Progranulin deficiency promotes circuit-specific synaptic
    pruning by microglia via complement activation.** *Cell* 2016 May 5.
    PMID **27114033** — <https://pubmed.ncbi.nlm.nih.gov/27114033/>
    The direct basis for `W_MG_PGRN` → `C1Q` → `C3` → `KPRUNE_C3`, i.e. the
    entire complement-mediated synapse-elimination arm, and for the model's
    prediction that anti-sortilin therapy lowers CSF C1q/C3.
19. 🔍 Anti-sortilin antibody preclinical proof of concept (raising PGRN by
    blocking sortilin) —
    <https://pubmed.ncbi.nlm.nih.gov/?term=anti-sortilin+antibody+progranulin+elevation+preclinical>

## 4. Progranulin-directed therapeutics

20. ✅ **[ANCHOR] Kurnellas M, et al. Latozinemab, a novel progranulin-elevating
    therapy for frontotemporal dementia.** *J Transl Med* 2023 Jun 15.
    PMID **37322482** — <https://pubmed.ncbi.nlm.nih.gov/37322482/>
21. ✅ **[ANCHOR] Ward M, et al. Phase 1 study of latozinemab in
    progranulin-associated frontotemporal dementia.** *Alzheimers Dement (N Y)* 2024 Jan-Mar.
    PMID **38356474** — <https://pubmed.ncbi.nlm.nih.gov/38356474/>
    Together the basis for the modelled 60 mg/kg IV q4w regimen and the target
    plasma/CSF PGRN elevations (plasma restored to the control range, CSF
    roughly doubled), which the model reproduces at 2.70× and 2.08×.
22. ✅ **[ANCHOR] Sevigny J, et al. Progranulin AAV gene therapy for
    frontotemporal dementia: translational studies and phase 1/2 trial
    interim results.** *Nat Med*
    2024 May. PMID **38745011** — <https://pubmed.ncbi.nlm.nih.gov/38745011/>
    Basis for the `AAVVG` → `TRANSD` gene-therapy arm and for the model's
    central differential prediction: gene therapy raises production and so
    escapes the receptor-blockade ceiling that limits the antibody.
23. ✅ **Peripheral expression of brain-penetrant progranulin rescues
    pathologies in mouse models of frontotemporal lobar degeneration.**
    *Sci Transl Med* 2024 Jun 5. PMID **38838131** —
    <https://pubmed.ncbi.nlm.nih.gov/38838131/>
24. ✅ **An anti-sortilin affibody-peptide fusion inhibits sortilin-mediated
    progranulin degradation.** *Front Immunol* 2024. PMID **39185427** —
    <https://pubmed.ncbi.nlm.nih.gov/39185427/>
25. 🔍 INFRONT-2 / INFRONT-3 (NCT04374136) phase 2/3 anti-sortilin results in
    GRN-FTD, including the primary CDR plus NACC-FTLD SB endpoint —
    <https://pubmed.ncbi.nlm.nih.gov/?term=INFRONT-3+latozinemab+GRN+frontotemporal+dementia+phase+3>
    and <https://clinicaltrials.gov/study/NCT04374136>
    The model's section D predicts a 0.65% 96-week slowing on this endpoint,
    i.e. a miss coexisting with unambiguous target engagement.

## 5. TDP-43 proteinopathy and cryptic splicing

26. ✅ **Neumann M, et al. Ubiquitinated TDP-43 in frontotemporal lobar
    degeneration and amyotrophic lateral sclerosis.** *Science* 2006 Oct 6.
    PMID **17023659** — <https://pubmed.ncbi.nlm.nih.gov/17023659/>
27. ✅ **Mackenzie IR, et al. A harmonized classification system for FTLD-TDP
    pathology.** *Acta Neuropathol* 2011 Jul. PMID **21644037** —
    <https://pubmed.ncbi.nlm.nih.gov/21644037/>
    Basis for the TDP type A (GRN) / B (C9orf72) / C (svPPA) nodes.
28. ✅ **[ANCHOR] Melamed Z, et al. Premature polyadenylation-mediated loss of
    stathmin-2 is a hallmark of TDP-43-dependent neurodegeneration.**
    *Nat Neurosci* 2019 Feb. PMID **30643298** —
    <https://pubmed.ncbi.nlm.nih.gov/30643298/>
29. ✅ **[ANCHOR] Klim JR, et al. ALS-implicated protein TDP-43 sustains levels
    of STMN2, a mediator of motor neuron growth and repair.** *Nat Neurosci*
    2019 Feb. PMID **30643292** — <https://pubmed.ncbi.nlm.nih.gov/30643292/>
    28 and 29 together are the basis for the `STMN2T` compartment and for
    `W_TOX_STMN2` / the axonal contribution to NfL release.
30. ✅ **[ANCHOR] Brown AL, et al. TDP-43 loss and ALS-risk SNPs drive
    mis-splicing and depletion of UNC13A.** *Nature* 2022 Mar.
    PMID **35197628** — <https://pubmed.ncbi.nlm.nih.gov/35197628/>
    Basis for the `UNC13AT` compartment and `KTOX_UNC13A`, i.e. cryptic
    splicing acting directly on synaptic function rather than via cell death.
31. 🔍 HDGFL2 cryptic peptide as a presymptomatic CSF marker of TDP-43
    dysfunction —
    <https://pubmed.ncbi.nlm.nih.gov/?term=HDGFL2+cryptic+exon+CSF+TDP-43+presymptomatic>

## 6. C9orf72 repeat toxicity

32. ✅ **Mori K, et al. The C9orf72 GGGGCC repeat is translated into aggregating
    dipeptide-repeat proteins in FTLD/ALS.** *Science* 2013 Mar 15.
    PMID **23393093** — <https://pubmed.ncbi.nlm.nih.gov/23393093/>
33. ✅ **Ash PE, et al. Unconventional translation of C9ORF72 GGGGCC expansion
    generates insoluble polypeptides specific to c9FTD/ALS.** *Neuron*
    2013 Feb 20. PMID **23415312** — <https://pubmed.ncbi.nlm.nih.gov/23415312/>
34. ✅ **Zu T, et al. RAN proteins and RNA foci from antisense transcripts in
    C9ORF72 ALS and frontotemporal dementia.** *PNAS* 2013 Dec 17.
    PMID **24248382** — <https://pubmed.ncbi.nlm.nih.gov/24248382/>
    32–34 are the basis for `RAN_translation` → `POLYGP` / `POLYGR` and for
    treating poly-GP as the CSF-measurable target-engagement marker while
    poly-GR carries the toxicity (`W_TOX_GR`). Section J of the model report
    flags that if the causal DPR is poly-GA instead, the ASO prediction changes
    and the marker trials use would be tracking the wrong arm.
35. 🔍 Arginine-rich dipeptide repeats (poly-GR / poly-PR), nucleolar stress and
    ribosomal toxicity —
    <https://pubmed.ncbi.nlm.nih.gov/?term=poly-GR+poly-PR+nucleolar+stress+C9orf72+toxicity>
36. 🔍 C9orf72 protein haploinsufficiency, the C9-SMCR8-WDR41 complex and
    autophagy/lysosome regulation — basis for `C9_LYSO_PEN` —
    <https://pubmed.ncbi.nlm.nih.gov/?term=C9orf72+SMCR8+WDR41+autophagy+lysosome+haploinsufficiency>
37. 🔍 WVE-004 / FOCUS-C9 — C9orf72 variant-selective ASO, CSF poly-GP lowering
    without clinical benefit; programme discontinued —
    <https://pubmed.ncbi.nlm.nih.gov/?term=WVE-004+C9orf72+antisense+oligonucleotide+poly-GP>
    and <https://clinicaltrials.gov/study/NCT04931862>

## 7. Tau pathology and tau-directed therapy

38. ✅ **Shulman M, et al. Exploratory analyses of clinical outcomes from the
    BIIB080 phase 1b study in mild Alzheimer's disease.** *Nat Aging* 2026 Feb.
    PMID **41673497** — <https://pubmed.ncbi.nlm.nih.gov/41673497/>
    The MAPT-lowering ASO class modelled in `ASOTCSF` / `ASOTTIS`.
39. 🔍 MAPT exon-10 splicing and the 4R:3R tau ratio in FTLD-tau (PSP, CBD,
    Pick's disease) —
    <https://pubmed.ncbi.nlm.nih.gov/?term=MAPT+exon+10+splicing+4R+3R+tau+ratio+FTLD>
40. 🔍 Tau oligomers as the synaptotoxic species, and prion-like trans-synaptic
    tau propagation — basis for `KTOX_TAU_SYN` and `KSEED_TAU` —
    <https://pubmed.ncbi.nlm.nih.gov/?term=tau+oligomers+synaptotoxicity+prion-like+propagation>
41. 🔍 Flortaucipir and the poor performance of first-generation tau-PET ligands
    for 4R tauopathies (the measurement gap noted in the map) —
    <https://pubmed.ncbi.nlm.nih.gov/?term=flortaucipir+4R+tauopathy+PSP+CBD+tau+PET+binding>

## 8. Selective vulnerability and network degeneration

42. ✅ **[ANCHOR] Seeley WW, et al. Early frontotemporal dementia targets neurons
    unique to apes and humans.** *Ann Neurol* 2006 Dec. PMID **17187353** —
    <https://pubmed.ncbi.nlm.nih.gov/17187353/>
    Von Economo neurons as the earliest cell type lost; basis for the
    `NEURSN` compartment being the salience-network pool and for `A_CDR_SN`
    carrying the largest CDR weight.
43. ✅ **Zhou J, et al. Divergent network connectivity changes in behavioural
    variant frontotemporal dementia and Alzheimer's disease.** *Brain*
    2010 May. PMID **20410145** — <https://pubmed.ncbi.nlm.nih.gov/20410145/>
    Salience-network breakdown with relative default-mode sparing — the
    mechanistic contrast with AD, and the `DMN_relative` node in the map.

## 9. Diagnostic criteria and clinical phenotypes

44. ✅ **Rascovsky K, et al. Sensitivity of revised diagnostic criteria for the
    behavioural variant of frontotemporal dementia.** *Brain* 2011 Sep.
    PMID **21810890** — <https://pubmed.ncbi.nlm.nih.gov/21810890/>
45. ✅ **Gorno-Tempini ML, et al. Classification of primary progressive aphasia
    and its variants.** *Neurology* 2011 Mar 15. PMID **21325651** —
    <https://pubmed.ncbi.nlm.nih.gov/21325651/>
    Basis for the svPPA / nfvPPA phenotype nodes and the `TEMP_SHIFT`
    regional-weighting parameter.

## 10. Biomarkers — **[ANCHOR]**

46. ✅ **[ANCHOR] Meeter LH, et al. Neurofilament light chain: a biomarker for
    genetic frontotemporal dementia.** *Ann Clin Transl Neurol* 2016 Aug.
    PMID **27606344** — <https://pubmed.ncbi.nlm.nih.gov/27606344/>
    Basis for `KPROD_NFL`, `KTR_NFL` and the `BH_NFL` survival term. The model
    is calibrated to plasma NfL ~10 pg/mL in controls versus ~50–80 pg/mL
    symptomatic, and CSF ~700 versus ~3000–5000 pg/mL.
47. 🔍 Serum/plasma NfL as a predictor of survival and progression rate in FTD —
    <https://pubmed.ncbi.nlm.nih.gov/?term=serum+neurofilament+light+frontotemporal+dementia+survival+prognosis>
48. 🔍 Low CSF pTau181/tTau ratio as an FTLD-TDP signature distinguishing it
    from Alzheimer's disease —
    <https://pubmed.ncbi.nlm.nih.gov/?term=CSF+p-tau+total+tau+ratio+FTLD-TDP+Alzheimer+differential>
49. 🔍 Plasma GFAP in frontotemporal dementia —
    <https://pubmed.ncbi.nlm.nih.gov/?term=plasma+GFAP+frontotemporal+dementia+biomarker>
50. 🔍 BMP (bis(monoacylglycero)phosphate) and glucosylsphingosine as lysosomal
    pharmacodynamic markers in GRN carriers — the readout the model argues
    should be measured alongside plasma PGRN —
    <https://pubmed.ncbi.nlm.nih.gov/?term=bis(monoacylglycero)phosphate+BMP+progranulin+GRN+lysosomal+biomarker>
51. ✅ **Analytical and clinical validation of a blood progranulin ELISA in
    frontotemporal dementias.** *Clin Chem Lab Med* 2023 Nov 27.
    PMID **37476993** — <https://pubmed.ncbi.nlm.nih.gov/37476993/>
    Relevant caveat: PGRN concentrations are assay-dependent, so the absolute
    ng/mL values in this model should be read as within-assay ratios.

## 11. Symptomatic pharmacology — including the negative trials

52. ✅ **[ANCHOR] Lebert F, et al. Frontotemporal dementia: a randomised,
    controlled trial with trazodone.** *Dement Geriatr Cogn Disord* 2004.
    PMID **15178953** — <https://pubmed.ncbi.nlm.nih.gov/15178953/>
    The one clearly positive behavioural RCT in FTD. Basis for `EMAX_TRZ`;
    the model reproduces ΔNPI = −5.8 against a reported ~6–10 point
    improvement.
53. ✅ **[ANCHOR] Mendez MF, et al. Preliminary findings: behavioral worsening on
    donepezil in patients with frontotemporal dementia.** *Am J Geriatr
    Psychiatry* 2007 Jan. PMID **17194818** —
    <https://pubmed.ncbi.nlm.nih.gov/17194818/>
    Basis for the `K_AGIT_DNP` harm term. The model reproduces ΔNPI = +3.2
    with ΔCDR of exactly zero.
54. ✅ **[ANCHOR] Boxer AL, et al. Memantine in patients with frontotemporal
    lobar degeneration: a multicentre, randomised, double-blind,
    placebo-controlled trial.** *Lancet Neurol* 2013 Feb. PMID **23290598** —
    <https://pubmed.ncbi.nlm.nih.gov/23290598/>
    The negative phase 3 that anchors the `Memantine_null` node in the map.
55. 🔍 Preserved cholinergic system in FTD (normal/near-normal ChAT and nucleus
    basalis) — the mechanistic reason AChEIs lack a rationale, encoded as
    `ACH_INTEGRITY` = 0.95 —
    <https://pubmed.ncbi.nlm.nih.gov/?term=cholinergic+system+preserved+frontotemporal+dementia+ChAT+nucleus+basalis>
56. 🔍 Serotonergic deficits in FTD: frontal 5-HT2A receptor loss and SSRI
    trials for behavioural symptoms —
    <https://pubmed.ncbi.nlm.nih.gov/?term=serotonergic+deficit+5-HT2A+frontotemporal+dementia+SSRI+behaviour>
57. 🔍 Intranasal oxytocin for social cognition and apathy in bvFTD —
    <https://pubmed.ncbi.nlm.nih.gov/?term=intranasal+oxytocin+behavioural+variant+frontotemporal+dementia+trial>

## 12. Trial design and outcome measures

58. ✅ **[ANCHOR] Boxer AL, et al. New directions in clinical trials for
    frontotemporal lobar degeneration: Methods and outcome measures.**
    *Alzheimers Dement* 2020 Jan. PMID **31668596** —
    <https://pubmed.ncbi.nlm.nih.gov/31668596/>
    Source for the CDR plus NACC-FTLD sum-of-boxes (0–24) primary endpoint and
    its progression rate; basis for `CDR_MAX` and `KCDR`.
59. ✅ **Silverman HE, et al. The contribution of behavioral features to
    caregiver burden in FTLD spectrum disorders.** *Alzheimers Dement* 2022 Sep. PMID **34854532** —
    <https://pubmed.ncbi.nlm.nih.gov/34854532/>
60. 🔍 ALLFTD / ARTFL-LEFFTDS longitudinal natural-history cohort —
    <https://pubmed.ncbi.nlm.nih.gov/?term=ALLFTD+ARTFL+LEFFTDS+longitudinal+frontotemporal+lobar+degeneration+cohort>
61. 🔍 FTLD-CDR / CDR plus NACC-FTLD psychometrics and annual progression rates —
    <https://pubmed.ncbi.nlm.nih.gov/?term=CDR+plus+NACC+FTLD+sum+of+boxes+progression+rate>
62. 🔍 Frontotemporal Dementia Rating Scale (FRS), Rasch-derived logit scoring —
    <https://pubmed.ncbi.nlm.nih.gov/?term=frontotemporal+dementia+rating+scale+FRS+Rasch+logit>

## 13. Methods — QSP / mrgsolve

63. 🔍 mrgsolve: simulation from ODE-based PK/PD and QSP models in R —
    <https://mrgsolve.org/> and
    <https://pubmed.ncbi.nlm.nih.gov/?term=mrgsolve+quantitative+systems+pharmacology+R>
64. 🔍 Target-mediated drug disposition and the quasi-equilibrium approximation
    (the `THETA_PL` / `THETA_CNS` formulation used here) —
    <https://pubmed.ncbi.nlm.nih.gov/?term=target-mediated+drug+disposition+quasi-equilibrium+approximation+Mager+Gibiansky>
65. 🔍 Antibody CSF/plasma exposure ratios in the CNS (~0.1–0.3%) — the basis
    for `KIN_CSF` / `KOUT_CSF` giving the model's computed 0.21% —
    <https://pubmed.ncbi.nlm.nih.gov/?term=monoclonal+antibody+CSF+plasma+ratio+0.1%25+brain+exposure>
66. 🔍 Intrathecal ASO CNS distribution and the rostro-caudal exposure gradient
    (`ASO9_ROSTRAL`) —
    <https://pubmed.ncbi.nlm.nih.gov/?term=intrathecal+antisense+oligonucleotide+CNS+distribution+rostrocaudal+gradient>

---

## Anchor summary — what was fitted to what

| Model quantity | Value used | Source | Model result |
|---|---|---|---|
| MAPT onset age | 49.5 y | #10 | 50.5 y (+1.0) |
| C9orf72 onset age | 58.2 y | #10 | 54.3 y (−3.9) |
| GRN onset age | 61.3 y | #10 | 55.8 y (−5.5, reported as structural failure) |
| Plasma PGRN carrier cut-off | 74.8 ng/mL | #17 | ~70 ng/mL vs 200 control |
| CSF PGRN cut-off | 3.43 ng/mL | #17 | 3.0 ng/mL control (assay-dependent) |
| Plasma PGRN rise on anti-sortilin | to control range | #20, #21 | 2.70× → 95% of control |
| CSF PGRN rise on anti-sortilin | ~2× | #20, #21 | 2.08× |
| CSF C1q / C3 fall on treatment | reduced | #18, #21 | −0.039 au each |
| Antibody CSF/plasma ratio | 0.1–0.3% | #65 | 0.21% |
| Plasma NfL, symptomatic | 50–80 pg/mL | #46 | 57.5 (GRN) |
| CSF NfL, symptomatic | 3000–5000 pg/mL | #46 | 3290 (GRN) |
| CDR+NACC-FTLD SB slope | 1.5–2.5 /y | #58, #61 | 1.58–1.76 /y |
| Brain atrophy rate | 2–3 %/y | #11 | 2.16–2.25 %/y |
| CSF poly-GP lowering, C9 ASO | up to ~50% | #37 | 50.1% |
| MAPT mRNA / tau lowering | 30–50% | #38 | 45.6% |
| Trazodone ΔNPI | −6 to −10 | #52 | −5.8 |
| Donepezil behavioural effect | worse | #53 | ΔNPI +3.2, ΔCDR 0.000 |
| TMEM106B vs plasma PGRN | no association | #17 | no association (not fitted) |

**Not anchored to anything, and decisive:** `P_LYS_SORT` — the sortilin share of
lysosomal progranulin delivery. No human measurement exists. It sets the *sign*
of the predicted anti-sortilin effect, and the plasma-PGRN biomarker is
mathematically incapable of measuring it. See section E of
`ftd_model_report.txt`.

---

*Compiled for the QSP Disease Model Library. Research and education only — not
for clinical or regulatory use. PMIDs marked ✅ were resolved against the NCBI
E-utilities API during model construction; entries marked 🔍 deliberately assert
no PMID.*
