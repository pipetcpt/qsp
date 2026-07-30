# High-Risk Neuroblastoma — QSP Model References

모든 PMID는 PubMed E-utilities API로 조회하여 제목·저자·저널·연도를 확인했습니다.
Every PMID below was resolved through the PubMed E-utilities API and the
author / journal / year / title were read back before being cited here.

References are grouped by the part of the model they constrain. Where a model
parameter is an **assumption** rather than a measured value, that is stated
explicitly — the reference then supports the *mechanism*, not the number.

---

## 1. Disease framework, risk stratification and outcomes

These set the clinical targets the model is calibrated against and the endpoint
definitions used in `nb_mrgsolve_model.R`.

1. Maris JM. **Recent advances in neuroblastoma.** *N Engl J Med* 2010.
   [PMID 20558371](https://pubmed.ncbi.nlm.nih.gov/20558371/)
   — the framing review; source of the ~50% long-term survival figure for
   high-risk disease that the model's overall-survival context assumes.
2. Cohn SL, Pearson ADJ, London WB, et al. **The International Neuroblastoma
   Risk Group (INRG) classification system: an INRG Task Force report.**
   *J Clin Oncol* 2009. [PMID 19047291](https://pubmed.ncbi.nlm.nih.gov/19047291/)
   — defines the high-risk group the model simulates and the response
   categories used for the `INRG` readout.
3. Ambros PF, Ambros IM, Brodeur GM, et al. **International consensus for
   neuroblastoma molecular diagnostics: report from the INRG Biology
   Committee.** *Br J Cancer* 2009.
   [PMID 19401703](https://pubmed.ncbi.nlm.nih.gov/19401703/)
4. Cohn SL, et al. **International neuroblastoma risk group consortium: a model
   of networking for rare cancers.** *J Natl Cancer Inst* 2026.
   [PMID 40854111](https://pubmed.ncbi.nlm.nih.gov/40854111/)
5. London WB, et al. **Adaptive Clinical Neuroblastoma Risk Groups — tailoring
   treatment in low- and middle-income countries: an INRG project.**
   *JCO Glob Oncol* 2025. [PMID 41435213](https://pubmed.ncbi.nlm.nih.gov/41435213/)
6. Perwein T, Lackner H, Sovinz P, et al. **Survival and late effects in
   children with stage 4 neuroblastoma.** *Pediatr Blood Cancer* 2011.
   [PMID 21319289](https://pubmed.ncbi.nlm.nih.gov/21319289/)
   — the late-effects burden (hypothyroidism, growth failure, gonadal failure)
   the model reports as `THY` and `ENDOLATE`.
7. Laverdière C, Cheung NK, Kushner BH, et al. **Long-term complications in
   survivors of advanced stage neuroblastoma.** *Pediatr Blood Cancer* 2005.
   [PMID 15714447](https://pubmed.ncbi.nlm.nih.gov/15714447/)

---

## 2. Genomic drivers — the parameters `MYCN`, `KPROL_MYCN`, `IC50ALK`, `SCA`

8. Mossé YP, Laudenslager M, Longo L, et al. **Identification of ALK as a major
   familial neuroblastoma predisposition gene.** *Nature* 2008.
   [PMID 18724359](https://pubmed.ncbi.nlm.nih.gov/18724359/)
   — the ALK F1174L / R1275Q distinction that drives the `IC50ALK` scenario.
9. Rajbhandari P, Lopez G, Capdevila C, et al. **Cross-cohort analysis
   identifies a TEAD4-MYCN positive feedback loop as the core regulatory
   element of high-risk neuroblastoma.** *Cancer Discov* 2018.
   [PMID 29510988](https://pubmed.ncbi.nlm.nih.gov/29510988/)
10. Wei JS, Kuznetsov IB, Zhang S, et al. **Clinically relevant cytotoxic
    immune cell signatures and clonal expansion of T-cell receptors in
    high-risk MYCN-not-amplified human neuroblastoma.** *Clin Cancer Res* 2018.
    [PMID 29784674](https://pubmed.ncbi.nlm.nih.gov/29784674/)
11. Ambros IM, Tonini GP, Pötschger U, et al. **Age dependency of the prognostic
    impact of tumor genomics in localized resectable MYCN-non-amplified
    neuroblastomas (SIOPEN LNESG / COG validation).** *J Clin Oncol* 2020.
    [PMID 32903140](https://pubmed.ncbi.nlm.nih.gov/32903140/)
    — segmental chromosomal aberration burden, the biological content of the
    map's `SCA` node.
12. Meeser A, Bartenhagen C, Werr L, et al. **Reliable assessment of telomere
    maintenance mechanisms in neuroblastoma.** *Cell Biosci* 2022.
    [PMID 36153564](https://pubmed.ncbi.nlm.nih.gov/36153564/)
13. Duan K, Dickson BC, Marrano P, et al. **Adult-onset neuroblastoma: report of
    seven cases with molecular genetic characterization.** *Genes Chromosomes
    Cancer* 2020. [PMID 31749253](https://pubmed.ncbi.nlm.nih.gov/31749253/)
    — the ATRX / ALT-positive, indolent-but-incurable phenotype in older
    patients; the model represents this only as a slower `KPROL`.
14. Liu Y, et al. **N6-methyladenosine-mediated overexpression of lncRNA
    ADAMTS9-AS2 triggers neuroblastoma differentiation via regulating
    LIN28B/let-7/MYCN signaling.** *JCI Insight* 2023.
    [PMID 37991019](https://pubmed.ncbi.nlm.nih.gov/37991019/)

---

## 3. Cell identity and the retinoid-responsive differentiation axis

These support the `TP → TD` transition and the claim that differentiated cells
are post-mitotic and therefore invisible to S/M-phase cytotoxics.

15. van Groningen T, Koster J, Valentijn LJ, et al. **Neuroblastoma is composed
    of two super-enhancer-associated differentiation states.** *Nat Genet* 2017.
    [PMID 28650485](https://pubmed.ncbi.nlm.nih.gov/28650485/)
    — the ADRN / MES axis; the model's `FMES` and `GD2_MES` parameters, i.e. the
    antigen-low escape route.
16. Boeva V, Louis-Brennetot C, Peltier A, et al. **Heterogeneity of
    neuroblastoma cell identity defined by transcriptional circuitries.**
    *Nat Genet* 2017. [PMID 28740262](https://pubmed.ncbi.nlm.nih.gov/28740262/)
17. Decaesteker B, Denecker G, Van Neste C, et al. **TBX2 is a neuroblastoma
    core regulatory circuitry component enhancing MYCN/FOXM1 reactivation of
    DREAM targets.** *Nat Commun* 2018.
    [PMID 30451831](https://pubmed.ncbi.nlm.nih.gov/30451831/)

---

## 4. GD2 antigen — density, biosynthesis, and modulation

18. van den Bijgaart RJE, Kroesen M, Wassink M, et al. **Combined sialic acid
    and histone deacetylase (HDAC) inhibitor treatment up-regulates the
    neuroblastoma antigen GD2.** *J Biol Chem* 2019.
    [PMID 30670592](https://pubmed.ncbi.nlm.nih.gov/30670592/)
    — supports the map's `VORINO ⊣ EZH2 → GD2S` edge and the possibility of
    rescuing antigen-low disease.
19. Kroesen M, Bull C, Gielen PR, et al. **Anti-GD2 mAb and vorinostat
    synergize in the treatment of neuroblastoma.** *Oncoimmunology* 2016.
    [PMID 27471639](https://pubmed.ncbi.nlm.nih.gov/27471639/)

> **Parameter note — `GD2DENS = 8×10⁶ molecules/cell`.** A direct, quantitative
> per-cell GD2 copy-number measurement in primary neuroblastoma is not something
> we could resolve to a specific PubMed record, and the model therefore treats
> the value as an **assumption within the commonly quoted 5–10×10⁶ range**. The
> conclusions that depend on it are *robust to it*: what matters for the model is
> only that antigen per gram (13.3 nmol/g at 8×10⁶/cell) massively exceeds what
> permeability can deliver per course. That inequality holds across the whole
> 10⁶–10⁷ range, and it is the inequality — not the number — that makes bound
> IgG per cell rather than fractional occupancy the correct PD driver. The
> `GD2DENS_NRV = 0.5×10⁶` value for nerve is a weaker assumption still and is
> flagged as such in `README.md`.

---

## 5. Anti-GD2 antibody — pharmacokinetics and clinical efficacy

20. Yu AL, Gilman AL, Ozkaynak MF, et al. **Anti-GD2 antibody with GM-CSF,
    interleukin-2, and isotretinoin for neuroblastoma.** *N Engl J Med* 2010.
    [PMID 20879881](https://pubmed.ncbi.nlm.nih.gov/20879881/)
    — **the primary calibration target (ANBL0032):** 2-year EFS 66% with
    immunotherapy + isotretinoin versus 46% with isotretinoin alone.
21. Marachelian A, Desai A, Balis F, et al. **Comparative pharmacokinetics,
    safety, and tolerability of two sources of ch14.18 in pediatric patients
    with high-risk neuroblastoma following myeloablative therapy.**
    *Cancer Chemother Pharmacol* 2016.
    [PMID 26791869](https://pubmed.ncbi.nlm.nih.gov/26791869/)
    — the PK target for `V1AB`, `CLAB`: Cmax ≈ 11.5 µg/mL and terminal
    t½ ≈ 10 d at 17.5 mg/m²/d × 4 d.
22. Siebert N, Eger C, Seidel D, et al. **Pharmacokinetics and pharmacodynamics
    of ch14.18/CHO in relapsed/refractory high-risk neuroblastoma patients
    treated by long-term infusion in combination with IL-2.** *MAbs* 2016.
    [PMID 26785755](https://pubmed.ncbi.nlm.nih.gov/26785755/)
    — the long-term-infusion schedule; the model's `ab_hours` parameter spans
    this to the naxitamab short-infusion extreme.
23. Siebert N, Jensen C, Troschke-Meurer S, et al. **Impact of HACA on
    immunomodulation and treatment toxicity following ch14.18/CHO long-term
    infusion with interleukin-2: results from a SIOPEN phase 2 trial.**
    *Cancers* 2018. [PMID 30336605](https://pubmed.ncbi.nlm.nih.gov/30336605/)
    — anti-drug antibody and accelerated clearance (map node `HAHA`; not
    implemented as an ODE, listed as a model limitation).
24. Ladenstein R, Pötschger U, Valteau-Couanet D, et al. **Interleukin 2 with
    anti-GD2 antibody ch14.18/CHO (dinutuximab beta) in patients with high-risk
    neuroblastoma (HR-NBL1/SIOPEN): a multicentre, randomised, phase 3 trial.**
    *Lancet Oncol* 2018. [PMID 30442501](https://pubmed.ncbi.nlm.nih.gov/30442501/)
    — **the second calibration target:** adding IL-2 produced no EFS benefit and
    more toxicity. Reproduced by the model as a −1.4% change in ADCC exposure.
25. Mora J, Chan GC, Morgenstern DA, et al. **The anti-GD2 monoclonal antibody
    naxitamab plus GM-CSF for relapsed or refractory high-risk neuroblastoma:
    a phase 2 clinical trial.** *Nat Commun* 2025.
    [PMID 39952926](https://pubmed.ncbi.nlm.nih.gov/39952926/)
26. Lode HN, et al. **Dinutuximab beta versus naxitamab in relapsed/refractory
    neuroblastoma with disease in bone or bone marrow: systematic review and
    matching-adjusted indirect comparison.** *Cancers* 2025.
    [PMID 40940820](https://pubmed.ncbi.nlm.nih.gov/40940820/)
27. Kushner BH, Ostrovnaya I, Cheung IY, et al. **Prolonged progression-free
    survival after consolidating second or later remissions of neuroblastoma
    with anti-GD2 immunotherapy and isotretinoin: a prospective phase II
    study.** *Oncoimmunology* 2015.
    [PMID 26140243](https://pubmed.ncbi.nlm.nih.gov/26140243/)
28. Cheung IY, Cheung NV, Modak S, et al. **Bone marrow minimal residual disease
    was an early response marker and a consistent independent predictor of
    survival after anti-GD2 immunotherapy.** *J Clin Oncol* 2015.
    [PMID 25559819](https://pubmed.ncbi.nlm.nih.gov/25559819/)
    — direct support for the model's central prediction that the **marrow** is
    the compartment anti-GD2 therapy actually clears.
29. Cheung NK, Cheung IY, Kushner BH, et al. **Key role for myeloid cells:
    phase II results of anti-GD2 antibody 3F8 plus GM-CSF for chemoresistant
    osteomedullary neuroblastoma.** *Int J Cancer* 2014.
    [PMID 24644014](https://pubmed.ncbi.nlm.nih.gov/24644014/)
    — the clinical counterpart of the model's granulocyte ADCC arm (`WG`,
    `EMAX_GM_ANC`): the myeloid compartment, not the NK compartment, is what
    GM-CSF recruits.
30. Oesterheld J, Ferguson W, Kraveka JM, et al. **Eflornithine as
    post-immunotherapy maintenance in high-risk neuroblastoma.**
    *J Clin Oncol* 2024. [PMID 37883734](https://pubmed.ncbi.nlm.nih.gov/37883734/)
    — a maintenance modality the model does **not** implement; noted as a gap.

---

## 6. The two Fc effector arms — ADCC (Hill n = 1) versus CDC (Hill n = 2)

This is the structural core of the model. The Hill-coefficient asymmetry is not
fitted; it follows from the stoichiometry of the two effector complexes.

31. Sopp JM, Peters SJ, Rowley TF, et al. **On-target IgG hexamerisation driven
    by a C-terminal IgM tail-piece fusion variant confers augmented complement
    activation.** *Commun Biol* 2021.
    [PMID 34475514](https://pubmed.ncbi.nlm.nih.gov/34475514/)
    — C1q engagement requires an ordered **multimer** of surface-bound IgG, the
    direct justification for `HCDC = 2` (a *pair* of adjacent Fc domains) and
    hence for CDC being superlinear in surface density while ADCC is linear.
32. van Osch TLJ, Nouta J, Derksen NIL, et al. **Fc galactosylation promotes
    hexamerization of human IgG1, leading to enhanced classical complement
    activation.** *J Immunol* 2021.
    [PMID 34408013](https://pubmed.ncbi.nlm.nih.gov/34408013/)
33. Hezareh M, Hessell AJ, Jensen RC, et al. **Effector function activities of a
    panel of mutants of a broadly neutralizing antibody against HIV-1.**
    *J Virol* 2001. [PMID 11711607](https://pubmed.ncbi.nlm.nih.gov/11711607/)
    — the **K322A** IgG1 Fc mutation: C1q binding and CDC are lost while
    FcγR-mediated ADCC is retained. The single most important reference for the
    model's therapeutic-index conclusion, and the source of `C1QEFF = 0.10`.
34. Navid F, Sondel PM, Barfield R, et al. **Phase I trial of a novel anti-GD2
    monoclonal antibody, hu14.18K322A, designed to decrease toxicity in children
    with refractory or recurrent neuroblastoma.** *J Clin Oncol* 2014.
    [PMID 24711551](https://pubmed.ncbi.nlm.nih.gov/24711551/)
    — the clinical realisation of the same idea: an Fc engineered to reduce
    complement activation, with reduced pain.
35. Bishop MW, Hutson PR, Hank JA, et al. **A phase 1 and pharmacokinetic study
    evaluating daily or weekly schedules of the humanized anti-GD2 antibody
    hu14.18K322A in recurrent/refractory solid tumors.** *MAbs* 2020.
    [PMID 32643524](https://pubmed.ncbi.nlm.nih.gov/32643524/)
36. Mise N, Takami M, Suzuki A, et al. **Antibody-dependent cellular
    cytotoxicity toward neuroblastoma enhanced by activated invariant natural
    killer T cells.** *Cancer Sci* 2016.
    [PMID 26749374](https://pubmed.ncbi.nlm.nih.gov/26749374/)
37. Delgado DC, Hank JA, Kolesar J, et al. **Genotypes of NK cell KIR receptors,
    their ligands, and Fcγ receptors in the response of neuroblastoma patients
    to hu14.18-IL2 immunotherapy.** *Cancer Res* 2010.
    [PMID 20935224](https://pubmed.ncbi.nlm.nih.gov/20935224/)
    — the basis for the `FCG` genotype multiplier.
38. Lode HN, et al. **Fcγ receptor polymorphism in patients with
    relapsed/refractory high-risk neuroblastoma correlates with outcomes in the
    SIOPEN dinutuximab beta long-term infusion trial.** *Clin Cancer Res* 2025.
    [PMID 40627545](https://pubmed.ncbi.nlm.nih.gov/40627545/)
    — the FcγR-genotype/outcome association the model's `FCG` axis predicts.
39. Brandetti E, Veneziani I, Melaiu O, et al. **MYCN is an immunosuppressive
    oncogene dampening the expression of ligands for NK-cell-activating
    receptors in human high-risk neuroblastoma.** *Oncoimmunology* 2017.
    [PMID 28680748](https://pubmed.ncbi.nlm.nih.gov/28680748/)
    — why an *antibody-directed* effector mechanism still works in MYCN-amplified
    disease when native NK recognition does not.
40. Veneziani I, Infante P, Ferretti E, et al. **The BET-bromodomain inhibitor
    JQ1 renders neuroblastoma cells more resistant to NK cell-mediated
    recognition and killing by downregulating ligands for NKG2D and DNAM-1.**
    *Oncotarget* 2019. [PMID 31040907](https://pubmed.ncbi.nlm.nih.gov/31040907/)

---

## 7. Delivery barriers — the 200-fold permeability spread that the model rests on

41. Jain RK. **Haemodynamic and transport barriers to the treatment of solid
    tumours.** *Int J Radiat Biol* 1991.
    [PMID 1678003](https://pubmed.ncbi.nlm.nih.gov/1678003/)
    — elevated interstitial fluid pressure and the loss of convective
    macromolecule influx; the mechanism behind the model's `VIFP` term, which is
    what makes anti-GD2 an MRD therapy rather than a debulking agent.
42. Flessner MF, Choi J, Vanpelt H, et al. **Tissue-level transport mechanisms
    of intraperitoneally-administered monoclonal antibodies.**
    *J Control Release* 1998. [PMID 9741914](https://pubmed.ncbi.nlm.nih.gov/9741914/)
    — antibody tissue-level transport, the order-of-magnitude anchor for
    `PSG_TU`.
43. Bush MS, Allt G. **Blood-nerve barrier: distribution of anionic sites on the
    endothelial plasma membrane and basal lamina of dorsal root ganglia.**
    *J Neurocytol* 1991. [PMID 1960538](https://pubmed.ncbi.nlm.nih.gov/1960538/)
    — the dorsal root ganglion microvasculature differs from that of the nerve
    trunk; the DRG is the model's high-permeability toxicity compartment
    (`PSG_NRV`), and this anatomical asymmetry is why the nociceptor sees more
    antibody per cell than a solid tumour does.
44. Mühlethaler-Mottet A, Liberman J, Ćetković H, et al. **The CXCR4/CXCR7/CXCL12
    axis is involved in a secondary but complex control of neuroblastoma
    metastatic cell homing.** *PLoS One* 2015.
    [PMID 25955316](https://pubmed.ncbi.nlm.nih.gov/25955316/)
    — the marrow niche (`BMNICHE` in the map).
45. Coniglio SJ. **Role of tumor-derived chemokines in osteolytic bone
    metastasis.** *Front Endocrinol* 2018.
    [PMID 29930538](https://pubmed.ncbi.nlm.nih.gov/29930538/)

> **Parameter note — `PSG_TU`, `PSG_BM`, `PSG_NRV`.** Only the tumour value has a
> published order-of-magnitude anchor. The marrow and DRG values are **assumed**
> to be 15× and 200× the tumour value respectively, on the qualitative grounds
> that both are barrier-free vascular beds while a solid tumour is not. The
> model's central conclusions are therefore conditional on that *ordering*, not
> on the numbers; `README.md` states the measurement that would falsify them.

---

## 8. Anti-GD2 pain — the dose-limiting toxicity

46. Bertolizio G, Ingelmo P, Cohen S, et al. **Multimodal analgesic plan for
    children undergoing chimeric 14.18 immunotherapy.**
    *J Pediatr Hematol Oncol* 2021.
    [PMID 31972721](https://pubmed.ncbi.nlm.nih.gov/31972721/)
    — the mandatory opioid/gabapentin co-therapy; the clinical scale of the
    `PAIN` state.
47. Mora J, et al. **Desensitizing the autonomic nervous system to mitigate
    anti-GD2 monoclonal antibody side effects.** *Front Oncol* 2024.
    [PMID 38812778](https://pubmed.ncbi.nlm.nih.gov/38812778/)

> **Mechanism note.** That anti-GD2 pain is *complement*-mediated rather than
> ADCC-mediated is the model's most consequential mechanistic commitment, and it
> rests on the K322A evidence above (33–35) — an Fc change that removes C1q
> binding reduces pain — rather than on a direct nerve-complement measurement,
> which we could not resolve to a specific record. This is listed as the model's
> first-order structural assumption in `README.md`.

---

## 9. Cytotoxic therapy, consolidation and stem-cell rescue

48. Matthay KK, Reynolds CP, Seeger RC, et al. **Long-term results for children
    with high-risk neuroblastoma treated on a randomized trial of myeloablative
    therapy followed by 13-cis-retinoic acid: a Children's Oncology Group
    study.** *J Clin Oncol* 2009.
    [PMID 19171716](https://pubmed.ncbi.nlm.nih.gov/19171716/)
49. Ladenstein R, Pötschger U, Pearson ADJ, et al. **Busulfan and melphalan
    versus carboplatin, etoposide, and melphalan as high-dose chemotherapy for
    high-risk neuroblastoma (HR-NBL1/SIOPEN).** *Lancet Oncol* 2017.
    [PMID 28259608](https://pubmed.ncbi.nlm.nih.gov/28259608/)
    — the basis for the model's `hdct = "BuMel"` deeper-log-kill setting.
50. Park JR, Kreissman SG, London WB, et al. **Effect of tandem autologous stem
    cell transplant vs single transplant on event-free survival in patients with
    high-risk neuroblastoma: a randomized clinical trial.** *JAMA* 2019.
    [PMID 31454045](https://pubmed.ncbi.nlm.nih.gov/31454045/)
    — the `tandem = TRUE` scenario.
51. Garaventa A, Poetschger U, Valteau-Couanet D, et al. **Randomized trial of
    two induction therapy regimens for high-risk neuroblastoma: HR-NBL1.5
    SIOPEN study.** *J Clin Oncol* 2021.
    [PMID 34152804](https://pubmed.ncbi.nlm.nih.gov/34152804/)
52. Park JR, Villablanca JG, London WB, et al. **Outcome of high-risk stage 3
    neuroblastoma with myeloablative therapy and 13-cis-retinoic acid.**
    *Pediatr Blood Cancer* 2009.
    [PMID 18937318](https://pubmed.ncbi.nlm.nih.gov/18937318/)

---

## 10. Myelosuppression — the Friberg structure

53. Friberg LE, Henningsson A, Maas H, et al. **Model of chemotherapy-induced
    myelosuppression with parameter consistency across drugs.**
    *J Clin Oncol* 2002. [PMID 12488418](https://pubmed.ncbi.nlm.nih.gov/12488418/)
    — the proliferative-pool + three-transit + feedback structure used verbatim
    for `PROL → TR1 → TR2 → TR3 → ANC`, including the `(CIRC0/ANC)^GAM` term.
54. Hansson EK, Friberg LE. **Limited inter-occasion variability in relation to
    inter-individual variability in chemotherapy-induced myelosuppression.**
    *Cancer Chemother Pharmacol* 2010.
    [PMID 19680655](https://pubmed.ncbi.nlm.nih.gov/19680655/)
55. Quartino AL, Friberg LE, Karlsson MO. **Characterization of endogenous G-CSF
    and the inverse correlation to chemotherapy-induced neutropenia using
    population modeling.** *Pharm Res* 2014.
    [PMID 24919931](https://pubmed.ncbi.nlm.nih.gov/24919931/)
    — the growth-factor feedback that the model collapses into `EMAX_GM_ANC`.
56. Guo Y, et al. **Optimization of clinical dosing schedule to manage
    neutropenia: learnings from a semi-mechanistic modeling simulation
    approach.** *J Pharmacokinet Pharmacodyn* 2020.
    [PMID 31853740](https://pubmed.ncbi.nlm.nih.gov/31853740/)

---

## 11. Isotretinoin — exposure, autoinduction and formulation

57. Veal GJ, Errington J, Sastry J, et al. **Pharmacokinetics and safety of a
    novel oral liquid formulation of 13-cis retinoic acid in children with
    neuroblastoma: a randomized crossover clinical trial.** *Cancers* 2021.
    [PMID 33919763](https://pubmed.ncbi.nlm.nih.gov/33919763/)
    — the formulation/administration problem the model encodes as `FREL`: how
    the drug is physically given to a small child changes its exposure.
58. Gota V, Chinnaswamy G, Vora T, et al. **Pharmacokinetics and
    pharmacogenetics of 13-cis retinoic acid in Indian high-risk neuroblastoma
    patients.** *Cancer Chemother Pharmacol* 2016.
    [PMID 27541143](https://pubmed.ncbi.nlm.nih.gov/27541143/)
    — the ~2–4 µM peak concentration range that `VRA`/`CLRA` are set to
    reproduce.
59. Nelson CH, Buttrick BR, Isoherranen N. **Therapeutic potential of the
    inhibition of the retinoic acid hydroxylases CYP26A1 and CYP26B1 by
    xenobiotics.** *Curr Top Med Chem* 2013.
    [PMID 23688132](https://pubmed.ncbi.nlm.nih.gov/23688132/)
60. Innes J, et al. **New insights into the role of CYP26 in retinoic acid
    clearance.** *Expert Opin Drug Metab Toxicol* 2026.
    [PMID 42023600](https://pubmed.ncbi.nlm.nih.gov/42023600/)
    — the autoinduction mechanism behind `FIND`, `EMAXI`, `KENZ`: a 14-day
    course progressively destroys its own exposure.

---

## 12. ¹³¹I-MIBG — transporter-limited delivery and dosimetry

61. Vallabhajosula S, Nikolopoulou A. **Radioiodinated
    metaiodobenzylguanidine (MIBG): radiochemistry, biology, and pharmacology.**
    *Semin Nucl Med* 2011. [PMID 21803182](https://pubmed.ncbi.nlm.nih.gov/21803182/)
    — the specific-activity / carrier-added versus no-carrier-added distinction
    that the model turns into a quantitative competition at the transporter
    (`SA0`, `KMNET`).
62. Gaze MN, Chang YC, Flux GD, et al. **Feasibility of dosimetry-based
    high-dose ¹³¹I-meta-iodobenzylguanidine with topotecan as a radiosensitizer
    in children with metastatic neuroblastoma.** *Cancer Biother Radiopharm*
    2005. [PMID 15869455](https://pubmed.ncbi.nlm.nih.gov/15869455/)
    — dosimetry-guided prescribing to a whole-body absorbed-dose target, the
    approach the model's `axis7_dosimetry` scenario quantifies.
63. Cash T, Yu AL, Weiss BD, et al. **Phase I study of ¹³¹I-metaiodobenzyl-
    guanidine with dinutuximab ± vorinostat for patients with relapsed or
    refractory neuroblastoma: a NANT trial.** *J Clin Oncol* 2025.
    [PMID 40549985](https://pubmed.ncbi.nlm.nih.gov/40549985/)
    — the combination the model is structured to simulate (radiopharmaceutical
    plus anti-GD2 antibody plus antigen up-regulation).
64. Campbell K, Shusterman S, Chi YY, et al. **Modulation of radiation
    biomarkers in a randomized phase II study of ¹³¹I-MIBG with or without
    radiation sensitizers for relapsed or refractory neuroblastoma.**
    *Int J Radiat Oncol Biol Phys* 2023.
    [PMID 36526235](https://pubmed.ncbi.nlm.nih.gov/36526235/)
65. Temple W, Mendelsohn L, Kim GE, et al. **Vesicular monoamine transporter
    protein expression correlates with clinical features, tumor biology, and
    MIBG avidity in neuroblastoma: a report from the Children's Oncology
    Group.** *Eur J Nucl Med Mol Imaging* 2016.
    [PMID 26338179](https://pubmed.ncbi.nlm.nih.gov/26338179/)
    — VMAT-dependent vesicular retention, i.e. the model's `KWASH`: uptake sets
    the delivered activity, retention sets the absorbed dose.
66. Batra V, Samanta M, Makvandi M, et al. **Dexmedetomidine does not interfere
    with meta-iodobenzylguanidine (MIBG) uptake at clinically relevant
    concentrations.** *Pediatr Blood Cancer* 2017.
    [PMID 27654664](https://pubmed.ncbi.nlm.nih.gov/27654664/)
    — a worked negative example of the drug-interference question the model's
    `NETX` parameter represents.
67. Das S, et al. **Synthesis and evaluation of ⁹⁹ᵐTc-analogues of
    [¹²³/¹³¹I]mIBG for targeting the norepinephrine transporter.**
    *Nucl Med Biol* 2019. [PMID 30770228](https://pubmed.ncbi.nlm.nih.gov/30770228/)
68. Pandit-Taskar N, Zanzonico PB, Staton KD, et al. **Biodistribution and
    dosimetry of ¹⁸F-meta-fluorobenzylguanidine: a first-in-human PET/CT imaging
    study of patients with neuroendocrine malignancies.** *J Nucl Med* 2018.
    [PMID 28705916](https://pubmed.ncbi.nlm.nih.gov/28705916/)

> **Parameter note — `SWB = 5.54×10⁻⁵ Gy/(MBq·d)`.** This is *calibrated*, not
> predicted: it is set so that the whole-body absorbed dose comes out at
> 0.20 mGy/MBq at the reference biological clearance. The model's dosimetry
> scenario therefore does not "predict" that number — it predicts how the
> required **activity** changes when biological clearance varies, which is a
> different and testable claim.

---

## 13. ALK inhibition — a free-fraction problem

69. Foster JH, Voss SD, Hall DC, et al. **Activity of crizotinib in patients
    with ALK-aberrant relapsed/refractory neuroblastoma: a Children's Oncology
    Group study (ADVL0912).** *Clin Cancer Res* 2021.
    [PMID 33568345](https://pubmed.ncbi.nlm.nih.gov/33568345/)
    — **the calibration target for the ALK axis:** crizotinib had activity in
    only a small minority of ALK-aberrant neuroblastomas.
70. Balis FM, Thompson PA, Mosse YP, et al. **First-dose and steady-state
    pharmacokinetics of orally administered crizotinib in children with solid
    tumors: a report on ADVL0912.** *Cancer Chemother Pharmacol* 2017.
    [PMID 28032129](https://pubmed.ncbi.nlm.nih.gov/28032129/)
    — the paediatric crizotinib exposure the model reproduces; combined with
    ~91% protein binding this is what puts free drug below the cellular IC50.
71. Greengard E, Mosse YP, Liu X, et al. **Safety, tolerability and
    pharmacokinetics of crizotinib in combination with cytotoxic chemotherapy
    for pediatric patients with refractory solid tumors or ALCL (ADVL1212).**
    *Cancer Chemother Pharmacol* 2020.
    [PMID 33095287](https://pubmed.ncbi.nlm.nih.gov/33095287/)
72. Goldsmith KC, Park JR, Kayser K, et al. **Lorlatinib with or without
    chemotherapy in ALK-driven refractory/relapsed neuroblastoma: phase 1 trial
    results.** *Nat Med* 2023.
    [PMID 37012551](https://pubmed.ncbi.nlm.nih.gov/37012551/)
    — the successful second-generation comparator; lower protein binding and a
    lower IC50 together, which is why the model separates the two effects.
73. Berko ER, Witek GM, Matkar S, et al. **Circulating tumor DNA reveals
    mechanisms of lorlatinib resistance in patients with relapsed/refractory
    ALK-driven neuroblastoma.** *Nat Commun* 2023.
    [PMID 37147298](https://pubmed.ncbi.nlm.nih.gov/37147298/)
74. Yang F, et al. **Population pharmacokinetic modeling and simulation of
    TQ-B3101 to inform dosing in pediatric patients with solid tumors.**
    *Front Pharmacol* 2021.
    [PMID 35115931](https://pubmed.ncbi.nlm.nih.gov/35115931/)

---

## 14. Platinum ototoxicity and the sodium-thiosulfate trade-off

75. Orgel E, Villaluna D, Krailo MD, et al. **Reevaluation of sodium thiosulfate
    otoprotection using the consensus International Society of Paediatric
    Oncology Ototoxicity Scale: a report from Children's Oncology Group study
    ACCL0431.** *Pediatr Blood Cancer* 2023.
    [PMID 37416942](https://pubmed.ncbi.nlm.nih.gov/37416942/)
    — the source of both halves of the model's `ESTS` / `FTUMSTS` trade-off:
    cochlear protection, against the concern that a thiol scavenger may also
    protect disseminated tumour.
76. Ohlsen TJD, et al. **Otoprotective effects of sodium thiosulfate by
    demographic and clinical characteristics: a report from Children's Oncology
    Group study ACCL0431.** *Pediatr Blood Cancer* 2025.
    [PMID 39654065](https://pubmed.ncbi.nlm.nih.gov/39654065/)
77. Freyer DR, Brock PR, Chang KW, et al. **Interventions for cisplatin-induced
    hearing loss in children and adolescents with cancer.**
    *Lancet Child Adolesc Health* 2019.
    [PMID 31160205](https://pubmed.ncbi.nlm.nih.gov/31160205/)

---

## 15. Biomarkers and response assessment

78. Grèze V, Chambon F, Merlin E, et al. **RT-qPCR for PHOX2B mRNA is a highly
    specific and sensitive method to assess neuroblastoma minimal residual
    disease.** *Oncol Lett* 2017.
    [PMID 28693243](https://pubmed.ncbi.nlm.nih.gov/28693243/)
    — the assay that reads the model's `TM` compartment directly, and therefore
    the measurement against which the model's central marrow prediction is
    testable.
79. Riaz S, Bashir H, Hassan A, et al. **I-131 mIBG scintigraphy Curie versus
    SIOPEN scoring: prognostic value in stage 4 neuroblastoma.**
    *Mol Imaging Radionucl Ther* 2018.
    [PMID 30317848](https://pubmed.ncbi.nlm.nih.gov/30317848/)
    — the semiquantitative imaging endpoint the model's `VTU`/`CURIE` readout
    stands in for.

---

## Software and methods

- **mrgsolve** — Baron KT. *mrgsolve: Simulate from ODE-Based Models.*
  <https://mrgsolve.org>
- **Graphviz** — Gansner ER, North SC. *An open graph visualization system and
  its applications to software engineering.* Softw Pract Exper 2000.
- **Shiny** — Chang W, et al. *shiny: Web Application Framework for R.*
- PubMed E-utilities — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

---

## 면책 (Disclaimer)

이 참고문헌 목록과 모델은 **연구·교육 목적**입니다. 임상 진료나 투약 결정에
사용해서는 안 됩니다. 파라미터 값 중 일부는 문헌 기반이지만 상당수는
가정(assumption)이며, 각 절의 *Parameter note*에 그 사실과 근거의 한계를
명시했습니다.

This reference list and the accompanying model are for **research and
educational purposes only** and must not be used for clinical care or dosing
decisions. Several parameters are assumptions rather than measured values; each
such case is flagged in a *Parameter note* above and in `README.md`.
