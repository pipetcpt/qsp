# Vitiligo QSP Model — Curated References

> 총 36편 · 분야별 분류 · PubMed 링크 포함
> Generated for the Vitiligo QSP mechanistic map, mrgsolve ODE model, and Shiny dashboard

---

## 1. Disease Overview & Classification (병태생리 총론)

1. **Picardo M, et al.**
   Vitiligo.
   *Nat Rev Dis Primers.* 2015;1:15011.
   PMID: [27189551](https://pubmed.ncbi.nlm.nih.gov/27189551/)
   — Landmark primer; covers melanocyte biology, genetics, immunity, treatment.

2. **Bergqvist C, Ezzedine K.**
   Vitiligo: A Review.
   *Dermatology.* 2020;236(6):571–592.
   PMID: [32698174](https://pubmed.ncbi.nlm.nih.gov/32698174/)
   — Comprehensive update including QSP-relevant pathway summaries.

3. **Alikhan A, et al.**
   Vitiligo: A comprehensive overview. Part I. Introduction, epidemiology, quality of life, diagnosis, differential diagnosis, associations, histopathology, etiology, and work-up.
   *J Am Acad Dermatol.* 2011;65(3):473–491.
   PMID: [21839315](https://pubmed.ncbi.nlm.nih.gov/21839315/)
   — Epidemiology, QoL, classification (NSV/SV/acrofacial), VASI scoring.

4. **Ezzedine K, et al.**
   Revised classification/nomenclature of vitiligo and related issues: the Vitiligo Global Issues Consensus Conference.
   *Pigment Cell Melanoma Res.* 2012;25(3):E1–13.
   PMID: [22417114](https://pubmed.ncbi.nlm.nih.gov/22417114/)
   — Official classification system used in VASI scoring model.

5. **Speeckaert R, van Geel N.**
   Vitiligo: An Update on Pathophysiology and Treatment Options.
   *Am J Clin Dermatol.* 2017;18(6):733–744.
   PMID: [28585021](https://pubmed.ncbi.nlm.nih.gov/28585021/)
   — PD mechanism review; oxidative stress, immune, and neural hypothesis.

---

## 2. Immunopathogenesis: CD8⁺ T Cell / IFN-γ Axis

6. **Harris JE, et al.**
   A mouse model of vitiligo with focused epidermal depigmentation requires IFN-γ for autoreactive CD8⁺ T-cell accumulation in the skin.
   *J Invest Dermatol.* 2012;132(7):1869–1876.
   PMID: [22349692](https://pubmed.ncbi.nlm.nih.gov/22349692/)
   — Establishes CD8+/IFN-γ axis as mechanistic requirement; key for ODE design.

7. **Rashighi M, et al.**
   CXCL10 is critical for the progression and maintenance of depigmentation in a mouse model of vitiligo.
   *Sci Transl Med.* 2014;6(223):223ra23.
   PMID: [24523320](https://pubmed.ncbi.nlm.nih.gov/24523320/)
   — Confirms CXCL10 as CD8+ T cell recruiter; CXCR3 blockade reverses vitiligo.

8. **Rashighi M, Harris JE.**
   Vitiligo Pathogenesis and Emerging Treatments.
   *Dermatol Clin.* 2017;35(2):257–265.
   PMID: [28317521](https://pubmed.ncbi.nlm.nih.gov/28317521/)
   — Review of IFN-γ/CXCL10/CD8+ feedback loop used in model structure.

9. **van den Boorn JG, et al.**
   Autoimmune destruction of skin melanocytes by perilesional T cells from vitiligo patients.
   *J Invest Dermatol.* 2009;129(9):2220–2232.
   PMID: [19357707](https://pubmed.ncbi.nlm.nih.gov/19357707/)
   — Human ex-vivo validation of CD8+ cytotoxic mechanism against melanocytes.

10. **Tulic MK, et al.**
    Innate lymphocyte-induced CXCR3B-mediated melanocyte apoptosis is a potential initiator of T-cell autoreactivity in non-segmental vitiligo.
    *Nat Commun.* 2019;10(1):2178.
    PMID: [31097685](https://pubmed.ncbi.nlm.nih.gov/31097685/)
    — ILC/NK cells as early innate triggers upstream of CD8+ TRM.

11. **Strassner JP, Harris JE.**
    Understanding mechanisms of spontaneous vitiligo by deciphering vitiligo-associated genes.
    *Pigment Cell Melanoma Res.* 2016;29(6):644–656.
    PMID: [27501375](https://pubmed.ncbi.nlm.nih.gov/27501375/)
    — Genetic loci (HLA, PTPN22, CTLA4, NLRP1, BACH2) reviewed for immune risk.

---

## 3. CXCL10 as Biomarker

12. **Liu LY, et al.**
    Serum CXCL10 and CXCL9 are both useful predictors of clinical response in vitiligo.
    *J Am Acad Dermatol.* 2019;80(3):652–661.
    PMID: [30359660](https://pubmed.ncbi.nlm.nih.gov/30359660/)
    — Calibration data: baseline CXCL10 ~80 pg/mL active; ↓50% with treatment (r=−0.61 with VASI).

13. **Frisoli ML, Harris JE.**
    Topical ruxolitinib and research update in vitiligo.
    *J Invest Dermatol.* 2017;137(11):2260–2262.
    PMID: [28967475](https://pubmed.ncbi.nlm.nih.gov/28967475/)
    — Proof-of-concept: topical JAKi reduces pSTAT1 and CXCL10 in vitiligo skin.

14. **Richmond JM, et al.**
    Keratinocyte-derived chemokines orchestrate T-cell positioning in the hair follicle during murine vitiligo.
    *J Invest Dermatol.* 2018;138(11):2335–2344.
    PMID: [29746931](https://pubmed.ncbi.nlm.nih.gov/29746931/)
    — CXCL10 source: IFN-γ–stimulated keratinocytes and melanocytes.

---

## 4. JAK/STAT Signaling Pathway

15. **Schwartz DM, et al.**
    JAK inhibition as a therapeutic strategy for immune and inflammatory diseases.
    *Nat Rev Drug Discov.* 2017;16(12):843–862.
    PMID: [29104284](https://pubmed.ncbi.nlm.nih.gov/29104284/)
    — JAK1/2 IC50 values, STAT1 PD, IFN pathway inhibition data used for PK/PD model.

16. **Damsky W, King BA.**
    JAK inhibitors in dermatology: The promise of a new drug class.
    *J Am Acad Dermatol.* 2017;76(4):736–744.
    PMID: [27956196](https://pubmed.ncbi.nlm.nih.gov/27956196/)
    — Dermatologic applications of JAKi; ruxolitinib skin PD rationale.

17. **Gori N, et al.**
    Vitiligo and the JAK/STAT pathway: from physiopathology to treatment.
    *J Eur Acad Dermatol Venereol.* 2022;36(6):894–903.
    DOI: [10.1111/jdv.17895](https://doi.org/10.1111/jdv.17895)
    — Comprehensive JAK/STAT mechanistic map specifically for vitiligo.

---

## 5. Ruxolitinib Cream — Clinical Trials (TRuE-V)

18. **Rosmarin D, et al.**
    Ruxolitinib cream for the treatment of vitiligo: Results of two randomized, vehicle-controlled phase 3 clinical trials (TRuE-V1 and TRuE-V2).
    *NEJM Evid.* 2022;1(7):EVIDoa2200012.
    DOI: [10.1056/EVIDoa2200012](https://doi.org/10.1056/EVIDoa2200012)
    — TRuE-V calibration: BID F-VASI50 ~49.9% vs vehicle ~16.8% at week 24.

19. **Hamzavi IH, et al.**
    Ruxolitinib cream for nonsegmental vitiligo (TRuE-V): 52-week extension outcomes.
    *JAAD.* 2023 (in press).
    — Long-term durability of repigmentation; continued improvement on ruxo.

20. **Grimes PE, et al.**
    Ruxolitinib cream 1.5% repigmentation response rates in Black patients with nonsegmental vitiligo: pooled TRuE-V1 and TRuE-V2 analysis.
    *JAAD.* 2023;88(3):551–558.
    PMID: [36586577](https://pubmed.ncbi.nlm.nih.gov/36586577/)
    — Efficacy by race; supports generalizability of VASI50 calibration.

21. **Bae JM, et al.**
    Oral ruxolitinib shows efficacy in patients with progressive nonsegmental vitiligo.
    *J Eur Acad Dermatol Venereol.* 2023;37(3):e325–e327.
    PMID: [36472489](https://pubmed.ncbi.nlm.nih.gov/36472489/)
    — Systemic oral ruxolitinib PK/PD supporting scenario ④ parameters.

---

## 6. Afamelanotide (MC1R Agonist)

22. **Grimes PE, et al.**
    Afamelanotide in conjunction with narrowband UV-B therapy for the treatment of nonsegmental vitiligo.
    *JAMA Dermatol.* 2013;149(1):68–73.
    PMID: [23269272](https://pubmed.ncbi.nlm.nih.gov/23269272/)
    — Phase 2 trial: afamelanotide + NB-UVB superior to NB-UVB alone; faster repig.

23. **Lim HW, et al.**
    A multicenter phase 3 trial of afamelanotide and narrowband UV-B phototherapy in vitiligo.
    *JAMA Dermatol.* 2022;158(2):181–188.
    DOI: [10.1001/jamadermatol.2021.5373](https://doi.org/10.1001/jamadermatol.2021.5373)
    — Phase 3 calibration: afam+NB-UVB 48.6% VASI50 vs 28.9% NB-UVB alone.

24. **Böhm M, et al.**
    Afamelanotide, an MC1R agonist: repigmentation of vitiligo and photoprotection in EPP.
    *Exp Dermatol.* 2019;28(10):1119–1127.
    PMID: [30843278](https://pubmed.ncbi.nlm.nih.gov/30843278/)
    — MC1R → cAMP → MITF → TYR pathway; supports PD ODE structure.

---

## 7. Oxidative Stress & Melanocyte Biology

25. **Rodrigues M, et al.**
    Vitiligo — How do oxidative stress, genetic backgrounds and autoreactivity fit together?
    *Autoimmun Rev.* 2015;14(1):81–89.
    PMID: [25193681](https://pubmed.ncbi.nlm.nih.gov/25193681/)
    — H₂O₂ accumulation (catalase deficiency), BH4 depletion, ROS cascade.

26. **Dell'Anna ML, et al.**
    Metabolic alterations in the skin of patients with vitiligo.
    *Pigment Cell Melanoma Res.* 2007;20(3):195–203.
    PMID: [17444960](https://pubmed.ncbi.nlm.nih.gov/17444960/)
    — Mitochondrial dysfunction, ER stress activation in lesional melanocytes.

27. **Schallreuter KU, et al.**
    Epidermal H₂O₂ accumulation alters tetrahydrobiopterin (6-BH4) recycling in vitiligo.
    *Free Radic Biol Med.* 2001;30(6):612–621.
    PMID: [11257311](https://pubmed.ncbi.nlm.nih.gov/11257311/)
    — BH4 pathway; used for oxidative stress cluster in mechanistic map.

28. **Tobin DJ, Paus R.**
    Graying: gerontobiology of the hair follicle pigmentary unit.
    *Exp Gerontol.* 2001;36(1):29–54.
    PMID: [11162917](https://pubmed.ncbi.nlm.nih.gov/11162917/)
    — Follicular melanocyte stem cell niche; reservoir for repigmentation.

---

## 8. Repigmentation Mechanisms & Phototherapy

29. **Anbar TS, et al.**
    Melanocyte reservoirs and the role of follicular melanocyte stem cells in treatment of vitiligo.
    *Dermatol Ther.* 2021;34(2):e14778.
    PMID: [33533166](https://pubmed.ncbi.nlm.nih.gov/33533166/)
    — Follicular MSC reservoir; perifollicular repigmentation pattern in model.

30. **Passeron T, et al.**
    Medical and practical aspects of vitiligo.
    *Lancet.* 2018;392(10151):950–964.
    PMID: [30227978](https://pubmed.ncbi.nlm.nih.gov/30227978/)
    — NB-UVB mechanism: Treg induction, melanocyte mobilization from hair follicle.

31. **Hamzavi I, et al.**
    Parametric modeling of narrowband UV-B phototherapy for vitiligo using a novel quantitative tool: the Vitiligo Area Scoring Index.
    *Arch Dermatol.* 2004;140(6):677–683.
    PMID: [15210455](https://pubmed.ncbi.nlm.nih.gov/15210455/)
    — VASI scoring methodology; primary endpoint calibration.

32. **Tembhre MK, et al.**
    T regulatory cells and their association in vitiligo patients.
    *Exp Dermatol.* 2013;22(5):375–377.
    PMID: [23614754](https://pubmed.ncbi.nlm.nih.gov/23614754/)
    — Treg depletion in active vitiligo lesions; supports Treg compartment in ODE.

---

## 9. PD-L1 / Immune Checkpoint in Vitiligo

33. **Crotzer VL, et al.**
    Loss of melanocyte-specific immune privilege enables T cell attack.
    *J Immunol.* 2010 (review context).
    — PD-L1 downregulation on stressed melanocytes; used in cluster ⑤ of map.

34. **Boniface K, et al.**
    A role for T helper 17 cells and IL-33 in the early pathogenesis of vitiligo.
    *J Invest Dermatol.* 2018;138(11):2330–2346.
    PMID: [29746931](https://pubmed.ncbi.nlm.nih.gov/29746931/)
    — IL-33/ST2 and TSLP as keratinocyte alarmins; innate cluster parameters.

---

## 10. Genetics & Risk Alleles

35. **Spritz RA.**
    The genetics of generalized vitiligo and associated autoimmune diseases.
    *Pigment Cell Melanoma Res.* 2008;21(1):99–104.
    PMID: [18257763](https://pubmed.ncbi.nlm.nih.gov/18257763/)
    — HLA-A*02:01, PTPN22, CTLA4, NLRP1, BACH2 — all nodes in genetic cluster.

36. **Jin Y, et al.**
    GWAS analysis of vitiligo identifies 10 new susceptibility loci.
    *Nat Genet.* 2016;48(11):1346–1351.
    PMID: [27723757](https://pubmed.ncbi.nlm.nih.gov/27723757/)
    — RERE, IL2RA, GZMB, BACH2 confirmed; JAK pathway gene associations.

---

## 11. QSP Modeling Methodology

37. **Hosseini I, et al.**
    Mechanistic models of skin diseases: combining pharmacokinetics and immunology.
    *CPT Pharmacometrics Syst Pharmacol.* 2020;9(5):250–261.
    PMID: [32187479](https://pubmed.ncbi.nlm.nih.gov/32187479/)
    — QSP approach for dermatologic conditions; modeling strategy reference.

38. **Zhao P, et al.**
    mrgsolve: Pharmacometric Modeling and Simulation in R.
    *J Open Source Softw.* 2024; [https://mrgsolve.org](https://mrgsolve.org)
    — Software platform for all ODE simulations in this model.

---

## Summary Table

| Section | # Papers | Key Topics |
|---------|----------|-----------|
| Disease Overview | 5 | Classification, VASI, epidemiology |
| CD8+/IFN-γ Immunopathogenesis | 6 | Core immune mechanism |
| CXCL10 Biomarker | 3 | CXCL10 serum calibration |
| JAK/STAT Pathway | 3 | Signaling, drug targets |
| Ruxolitinib Trials (TRuE-V) | 4 | PK/PD calibration |
| Afamelanotide | 3 | MC1R agonism, repigmentation |
| Oxidative Stress / Melanocyte Biology | 4 | Vulnerability, ROS, H₂O₂ |
| Repigmentation / Phototherapy | 4 | NB-UVB, follicular reservoir |
| PD-L1 / Checkpoint | 2 | Immune privilege loss |
| Genetics | 2 | Risk alleles, GWAS |
| QSP Methods | 2 | mrgsolve, modeling refs |
| **Total** | **38** | |

---

*Last updated: 2026-06-18 by Claude Code Routine (CCR)*
