# GIST QSP model — references

Every PMID below was resolved against PubMed (E-utilities `esearch`/`esummary`)
while this model was built, so the identifier, author, journal and year match the
paper named. Each entry says what the model actually takes from it: a parameter,
a structural commitment, or a calibration/validation anchor.

Notation used in the annotations:
**[PARAM]** a number in `$PARAM` / `P` · **[STRUCT]** a reason the equations are
shaped the way they are · **[ANCHOR]** a clinical quantity the model is fitted to
· **[TEST]** a clinical quantity the model predicts and is checked against ·
**[MISS]** an observation the model does not fully reproduce, reported as such.

---

## 1 · Disease biology, cell of origin and molecular classification

1. Corless CL, Barnett CM, Heinrich MC. **Gastrointestinal stromal tumours:
   origin and molecular oncology.** *Nat Rev Cancer* 2011.
   [PMID 22089421](https://pubmed.ncbi.nlm.nih.gov/22089421/)
   — [STRUCT] mutually exclusive driver genotypes; the frequencies used in the
   `GENO` switch (exon 11 ~65-70%, exon 9 ~8-10%, PDGFRA ~5-8%, wild type ~10-15%).

2. Chi P, Chen Y, Zhang L, et al. **ETV1 is a lineage survival factor that
   cooperates with KIT in gastrointestinal stromal tumours.** *Nature* 2010.
   [PMID 20927104](https://pubmed.ncbi.nlm.nih.gov/20927104/)
   — [STRUCT] the interstitial-cell-of-Cajal lineage programme and the
   ERK→ETV1→KIT feed-forward loop drawn in cluster 1 of the map.

3. Ou WB, Zhu MJ, Demetri GD, Fletcher CDM, Fletcher JA. **Protein kinase
   C-theta regulates KIT expression and proliferation in gastrointestinal
   stromal tumors.** *Oncogene* 2008.
   [PMID 18521081](https://pubmed.ncbi.nlm.nih.gov/18521081/)
   — [STRUCT] PLCγ→PKCθ→NF-κB arm and the transcriptional support of KIT.

4. Boikos SA, Pappo AS, Killian JK, et al. **Molecular subtypes of KIT/PDGFRA
   wild-type gastrointestinal stromal tumors: a report from the NIH GIST
   Clinic.** *JAMA Oncol* 2016.
   [PMID 27011036](https://pubmed.ncbi.nlm.nih.gov/27011036/)
   — [PARAM] the `wt_sdh` phenotype: `kitdep = 0.15` (largely KIT-independent)
   and `kpscale = 0.65` (indolent), i.e. no TKI target and slow growth.

5. Zhao R, Wang Y, Huang Y, et al. **Molecular landscape and clinical
   significance of exon 11 mutations in KIT gene among patients with
   gastrointestinal stromal tumours.** *Front Oncol* 2023.
   [PMID 37901323](https://pubmed.ncbi.nlm.nih.gov/37901323/)
   — [STRUCT] heterogeneity within exon 11 (del557-558 vs point mutations);
   represented here by a single exon 11 EC50 and noted as a simplification.

6. Brahmi M, Alberti L, Tirode F, et al. **KIT exon 10 variant (c.1621 A > C)
   single nucleotide polymorphism as predictor of GIST patient outcome.**
   *BMC Cancer* 2015.
   [PMID 26498480](https://pubmed.ncbi.nlm.nih.gov/26498480/)
   — [STRUCT] germline modifiers of outcome; not in the model, listed as a
   known omission.

---

## 2 · KIT/PDGFRA structure, conformational switch control and drug classes

7. Smith BD, Kaufman MD, Lu WP, et al. **Ripretinib (DCC-2618) is a switch
   control kinase inhibitor of a broad spectrum of oncogenic and drug-resistant
   KIT and PDGFRA variants.** *Cancer Cell* 2019.
   [PMID 31085175](https://pubmed.ncbi.nlm.nih.gov/31085175/)
   — [PARAM] the ripretinib EC50 column: potent against the activation loop
   (`C2_RI = 165`) and comparatively weak against the ATP-binding pocket
   (`C1_RI = 2200`). This asymmetry is the mechanism of the INTRIGUE crossover
   and was taken from here, **not** fitted to INTRIGUE.

8. Serrano C, Mariño-Enríquez A, Tao DL, et al. **Complementary activity of
   tyrosine kinase inhibitors against secondary KIT mutations in
   imatinib-resistant gastrointestinal stromal tumours.** *Br J Cancer* 2019.
   [PMID 30792533](https://pubmed.ncbi.nlm.nih.gov/30792533/)
   — [PARAM] the whole drug × clone potency matrix pattern: sunitinib covers
   exon 13/14, regorafenib and ripretinib cover exon 17/18, and no single agent
   covers both. This is the single most important non-clinical input.

9. Heinrich MC, Maki RG, Corless CL, et al. **Primary and secondary kinase
   genotypes correlate with the biological and clinical activity of sunitinib in
   imatinib-resistant gastrointestinal stromal tumor.** *J Clin Oncol* 2008.
   [PMID 18955458](https://pubmed.ncbi.nlm.nih.gov/18955458/)
   — [PARAM] `C1_SU = 12` vs `C2_SU = 380`: sunitinib's activity is confined to
   ATP-pocket secondary mutants. Also the exon 9 > exon 11 ordering of sunitinib
   benefit (`E9_SU = 15` < `E11_SU = 20`).

10. Antonescu CR, Besmer P, Guo T, et al. **Acquired resistance to imatinib in
    gastrointestinal stromal tumor occurs through secondary gene mutation.**
    *Clin Cancer Res* 2005.
    [PMID 15930355](https://pubmed.ncbi.nlm.nih.gov/15930355/)
    — [STRUCT] resistance is a *mutation in the same gene*, not loss of the
    target; sets the two secondary classes and their relative frequency
    (`MU_AL > MU_ATP`).

---

## 3 · Resistance architecture: polyclonality and ctDNA  [C1]

11. Wardelmann E, Merkelbach-Bruse S, Pauls K, et al. **Polyclonal evolution of
    multiple secondary KIT mutations in gastrointestinal stromal tumors under
    treatment with imatinib mesylate.** *Clin Cancer Res* 2006.
    [PMID 16551858](https://pubmed.ncbi.nlm.nih.gov/16551858/)
    — [STRUCT] the reason the state variable is a clone *vector*: different
    secondary mutations in different nodules of the same patient.

12. Liegl B, Kepten I, Le C, et al. **Heterogeneity of kinase inhibitor
    resistance mechanisms in GIST.** *J Pathol* 2008.
    [PMID 18623623](https://pubmed.ncbi.nlm.nih.gov/18623623/)
    — [STRUCT] ≥2 resistance mutations within a single nodule; no single agent
    can cover one patient's tumour, which is what makes efficacy a set-cover
    problem rather than a potency problem.

13. Heinrich MC, Jones RL, George S, et al. **Ripretinib versus sunitinib in
    gastrointestinal stromal tumor: ctDNA biomarker analysis of the phase 3
    INTRIGUE trial.** *Nat Med* 2024.
    [PMID 38182785](https://pubmed.ncbi.nlm.nih.gov/38182785/)
    — [TEST] the model's headline prediction. Observed: KIT exon 11+13/14 →
    sunitinib 15.0 vs ripretinib 4.0 months; exon 11+17/18 → ripretinib 14.2 vs
    sunitinib 1.5 months; ctDNA detected in 280/362 (77%) with KIT mutations in
    213/362 (59%); the two subgroups were *mutually exclusive*. Nothing in the
    model was fitted to these numbers.

14. Bauer S, Jones RL, Blay JY, et al. **Ripretinib versus sunitinib in patients
    with advanced gastrointestinal stromal tumor after treatment with imatinib
    (INTRIGUE): a randomized, open-label, phase III trial.** *J Clin Oncol* 2022.
    [PMID 35947817](https://pubmed.ncbi.nlm.nih.gov/35947817/)
    — [TEST] the ITT null (8.0 vs 8.3 months; KIT exon 11 ITT 8.3 vs 7.0;
    ORR 23.9% vs 14.6%) that the model explains as the average of two
    opposite-signed genotype effects.

15. Kang YK, Ryu MH, Yoo C, et al. **Resumption of imatinib to control
    metastatic or unresectable gastrointestinal stromal tumours after failure of
    imatinib and sunitinib (RIGHT): a randomised, placebo-controlled, phase 3
    trial.** *Lancet Oncol* 2013.
    [PMID 24140183](https://pubmed.ncbi.nlm.nih.gov/24140183/)
    — [TEST] PFS 1.8 vs 0.9 months. The authors' own conclusion — "the disease
    continues to harbour many clones that are sensitive to kinase inhibitors" —
    is commitment C1 stated clinically. The model reproduces the direction
    (17% of viable cells still carry only the primary mutation at fourth line)
    with a **[MISS]** on effect size.

16. Serrano C, Leal A, Kuang Y, et al. **Phase I study of rapid alternation of
    sunitinib and regorafenib for the treatment of tyrosine kinase inhibitor
    refractory gastrointestinal stromal tumors.** *Clin Cancer Res* 2019.
    [PMID 31471313](https://pubmed.ncbi.nlm.nih.gov/31471313/)
    — [STRUCT] the clinical translation of complementary coverage: alternating
    two agents with orthogonal clone spectra. Scenario S21's upfront-combination
    prediction is the same idea with simultaneous rather than alternating dosing.

---

## 4 · First-line imatinib: efficacy, genotype and dose  [C3]

17. Demetri GD, von Mehren M, Blanke CD, et al. **Efficacy and safety of
    imatinib mesylate in advanced gastrointestinal stromal tumors.**
    *N Engl J Med* 2002.
    [PMID 12181401](https://pubmed.ncbi.nlm.nih.gov/12181401/)
    — [ANCHOR] B2222: the response rate and the fact that complete responses are
    essentially absent, which the model reproduces as a partial response of about
    −45% SLD with a persistent quiescent reservoir.

18. Demetri GD, von Mehren M, Blanke CD, et al. **Identification and treatment
    of chemoresistant inoperable or metastatic GIST: experience with the
    selective tyrosine kinase inhibitor imatinib mesylate.**
    *Eur J Cancer* 2002.
    [PMID 12528773](https://pubmed.ncbi.nlm.nih.gov/12528773/)
    — [ANCHOR] early dose-finding context for the 400-600 mg range.

19. Heinrich MC, Corless CL, Demetri GD, et al. **Kinase mutations and imatinib
    response in patients with metastatic gastrointestinal stromal tumor.**
    *J Clin Oncol* 2003.
    [PMID 14645423](https://pubmed.ncbi.nlm.nih.gov/14645423/)
    — [ANCHOR] response by genotype: exon 11 ≫ exon 9 ≫ wild type. Sets the
    ordering of `E11_IM` < `E9_IM` ≪ `D842_IM`.

20. Debiec-Rychter M, Sciot R, Le Cesne A, et al. **KIT mutations and dose
    selection for imatinib in patients with advanced gastrointestinal stromal
    tumours.** *Eur J Cancer* 2006.
    [PMID 16624552](https://pubmed.ncbi.nlm.nih.gov/16624552/)
    — [ANCHOR] EORTC 62005: exon 9 raises the relative risk of progression by
    171% versus exon 11, and 800 mg reduces the relative risk of progression by
    61% **in exon 9 only**. Combined with an exon 11 median PFS of ~25 months
    these give the derived targets used here: exon 9 at 400 mg ≈ 9 months, exon 9
    at 800 mg ≈ 23 months.

21. Heinrich MC, Owzar K, Corless CL, et al. **Correlation of kinase genotype
    and clinical outcome in the North American Intergroup Phase III Trial of
    imatinib mesylate for treatment of advanced gastrointestinal stromal tumor:
    CALGB 150105 study by Cancer and Leukemia Group B and Southwest Oncology
    Group.** *J Clin Oncol* 2008.
    [PMID 18955451](https://pubmed.ncbi.nlm.nih.gov/18955451/)
    — [ANCHOR] S0033: exon 11 median PFS ~25 months, the model's first-line
    calibration target for `TURNOVER`.

22. Gastrointestinal Stromal Tumor Meta-Analysis Group (MetaGIST).
    **Comparison of two doses of imatinib for the treatment of unresectable or
    metastatic gastrointestinal stromal tumors: a meta-analysis of 1,640
    patients.** *J Clin Oncol* 2010.
    [PMID 20124181](https://pubmed.ncbi.nlm.nih.gov/20124181/)
    — [TEST] "a small PFS advantage of high-dose imatinib, essentially among
    patients with KIT exon 9 mutations, but no OS advantage." The model produces
    this as a threshold-geometry consequence: exon 11 is already on the plateau
    at 400 mg, exon 9 crosses the proliferation threshold between 400 and 800 mg.

23. Gronchi A, Blay JY, Trent JC. **The role of high-dose imatinib in the
    management of patients with gastrointestinal stromal tumor.** *Cancer* 2010.
    [PMID 20166214](https://pubmed.ncbi.nlm.nih.gov/20166214/)
    — [STRUCT] review of when escalation is and is not rational.

24. Judson I. **Imatinib in advanced gastrointestinal stromal tumour: when is
    800 mg the correct dose?** *Curr Opin Oncol* 2008.
    [PMID 18525340](https://pubmed.ncbi.nlm.nih.gov/18525340/)
    — [STRUCT] the clinical framing of C3.

25. Blay JY, Shen L, Kang YK, et al. **Nilotinib versus imatinib as first-line
    therapy for patients with unresectable or metastatic gastrointestinal
    stromal tumours (ENESTg1): a randomised phase 3 trial.**
    *Lancet Oncol* 2015.
    [PMID 25882987](https://pubmed.ncbi.nlm.nih.gov/25882987/)
    — [TEST] a first-line agent with a *narrower* clone spectrum than imatinib
    loses; consistent with coverage rather than potency being the currency.

26. Zhou Y, Zhang X, Wu X, et al. **A prospective multicenter phase II study on
    the efficacy and safety of dasatinib in the treatment of metastatic
    gastrointestinal stromal tumor.** *Cancer Med* 2020.
    [PMID 32677196](https://pubmed.ncbi.nlm.nih.gov/32677196/)
    — [TEST] same argument for a different off-spectrum agent.

---

## 5 · Imatinib pharmacokinetics, exposure–response and protein binding

27. Demetri GD, Wang Y, Wehrle E, et al. **Imatinib plasma levels are correlated
    with clinical benefit in patients with unresectable/metastatic
    gastrointestinal stromal tumors.** *J Clin Oncol* 2009.
    [PMID 19451435](https://pubmed.ncbi.nlm.nih.gov/19451435/)
    — [ANCHOR/MISS] median TTP 11.3 months in the lowest trough quartile
    (Cmin < 1110 ng/mL) versus > 30 months in quartiles 2-4. The model's
    steady-state trough at 400 mg is 1205 ng/mL (parent + CGP74588), which
    matches the measured distribution; but the model makes exon 11 first-line PFS
    almost exposure-*independent* above ~250 ng/mL, so it under-predicts this
    association. Reported as a **[MISS]** in `README.md` §9.

28. Widmer N, Decosterd LA, Leyvraz S, et al. **Relationship of imatinib-free
    plasma levels and target genotype with efficacy and tolerability.**
    *Br J Cancer* 2008.
    [PMID 18475296](https://pubmed.ncbi.nlm.nih.gov/18475296/)
    — [PARAM] the free-fraction/genotype interaction that motivates
    `EC50_eff = EC50 · (AGP/AGP0)^HAGP`, `HAGP = 0.90`.

29. Haouala A, Widmer N, Guidi M, et al. **Prediction of free imatinib
    concentrations based on total plasma concentrations in patients with
    gastrointestinal stromal tumours.** *Br J Clin Pharmacol* 2013.
    [PMID 22891806](https://pubmed.ncbi.nlm.nih.gov/22891806/)
    — [PARAM] the inverse relation between AGP and free drug; sets `WAGP = 0.25`
    (the fractional AGP rise with tumour burden) and the exponent above.

30. Gandia P, Arellano C, Lafont T, et al. **Should therapeutic drug monitoring
    of the unbound fraction of imatinib and its main active metabolite
    N-desmethyl-imatinib be developed?** *Cancer Chemother Pharmacol* 2013.
    [PMID 23183914](https://pubmed.ncbi.nlm.nih.gov/23183914/)
    — [PARAM] CGP74588 is treated as equipotent with the parent and its
    concentration is added (`FM_IM = 0.12`, `CLM_IM = 250 L/day`, giving ~14% of
    parent exposure).

31. Bouchet S, Poulette S, Titier K, et al. **Relationship between imatinib
    trough concentration and outcomes in the treatment of advanced
    gastrointestinal stromal tumours in a real-life setting.**
    *Eur J Cancer* 2016.
    [PMID 26851399](https://pubmed.ncbi.nlm.nih.gov/26851399/)
    — [ANCHOR] the real-world trough distribution the model's 1205 ng/mL sits in.

---

## 6 · Second, third and fourth line

32. Demetri GD, van Oosterom AT, Garrett CR, et al. **Efficacy and safety of
    sunitinib in patients with advanced gastrointestinal stromal tumour after
    failure of imatinib: a randomised controlled trial.** *Lancet* 2006.
    [PMID 17046465](https://pubmed.ncbi.nlm.nih.gov/17046465/)
    — [TEST] time to progression 27.3 weeks versus 6.4 weeks on placebo
    (HR 0.33), and the 4-weeks-on/2-weeks-off schedule that the model's
    `on_schedule("4/2")` reproduces, including regrowth during the off-weeks.

33. Demetri GD, Reichardt P, Kang YK, et al. **Efficacy and safety of
    regorafenib for advanced gastrointestinal stromal tumours after failure of
    imatinib and sunitinib (GRID): an international, multicentre, randomised,
    placebo-controlled, phase 3 trial.** *Lancet* 2013.
    [PMID 23177515](https://pubmed.ncbi.nlm.nih.gov/23177515/)
    — [TEST] PFS 4.8 versus 0.9 months (HR 0.27) in third line, and the
    3-weeks-on/1-week-off schedule.

34. Blay JY, Serrano C, Heinrich MC, et al. **Ripretinib in patients with
    advanced gastrointestinal stromal tumours (INVICTUS): a double-blind,
    randomised, placebo-controlled, phase 3 trial.** *Lancet Oncol* 2020.
    [PMID 32511981](https://pubmed.ncbi.nlm.nih.gov/32511981/)
    — [TEST] fourth-line PFS 6.3 versus 1.0 months. A broad-spectrum agent
    working where narrow ones have failed is the set-cover prediction.

35. Heinrich MC, Jones RL, von Mehren M, et al. **Avapritinib in advanced PDGFRA
    D842V-mutant gastrointestinal stromal tumour (NAVIGATOR): a multicentre,
    open-label, phase 1 trial.** *Lancet Oncol* 2020.
    [PMID 32615108](https://pubmed.ncbi.nlm.nih.gov/32615108/)
    — [TEST] the same node, opposite drug ranks: D842V is the most
    imatinib-resistant genotype and the most avapritinib-sensitive
    (`D842_IM = 20000` vs `D842_AV = 40`).

36. Kang YK, George S, Jones RL, et al. **Avapritinib versus regorafenib in
    locally advanced unresectable or metastatic GI stromal tumor: a randomized,
    open-label phase III study (VOYAGER).** *J Clin Oncol* 2021.
    [PMID 34343033](https://pubmed.ncbi.nlm.nih.gov/34343033/)
    — [TEST] a highly potent activation-loop inhibitor does *not* beat
    regorafenib in an unselected later-line population — coverage of the clones
    actually present, not potency against one of them, is what matters.

37. Kurokawa Y, Honma Y, Sawaki A, et al. **Pimitespib in patients with advanced
    gastrointestinal stromal tumor (CHAPTER-GIST-301): a randomized,
    double-blind, placebo-controlled phase III trial.** *Ann Oncol* 2022.
    [PMID 35688358](https://pubmed.ncbi.nlm.nih.gov/35688358/)
    — [STRUCT] an HSP90 inhibitor acts *below* the clone-specific EC50 matrix and
    is therefore genotype-agnostic; a target the model does not yet contain.

38. Boilève A, Faron M, Fromentin AM, et al. **Outcomes of patients with
    metastatic gastrointestinal stromal tumors (GIST) treated with multi-kinase
    inhibitors beyond the third line.** *ESMO Open* 2020.
    [PMID 33246932](https://pubmed.ncbi.nlm.nih.gov/33246932/)
    — [TEST] the shortening of PFS with each successive line, which the model
    produces as the KIT-independent clone becoming dominant.

---

## 7 · Stopping the drug, rechallenge and the reservoir  [C2]

39. Blay JY, Le Cesne A, Cassier PA, et al. **Discontinuation versus continuation
    of imatinib in patients with advanced gastrointestinal stromal tumours
    (BFR14): long-term results of an open-label, randomised phase 3 trial.**
    *Lancet Oncol* 2024.
    [PMID 39127063](https://pubmed.ncbi.nlm.nih.gov/39127063/)
    — [ANCHOR/TEST] the central test of C2. Interruption versus continuation
    median PFS: 6.1 vs 27.8 months after 1 year, 7.0 vs 67.0 after 3 years,
    12.0 vs not reached after 5 years. `KDQ` and `PHI_D` are fitted to the
    1-year arm; the 3-year and 5-year arms are predictions, with the phenotype
    of each randomised cohort taken from its own continuation arm.

40. Le Cesne A, Ray-Coquard I, Bui BN, et al. **Discontinuation of imatinib in
    patients with advanced gastrointestinal stromal tumours after 3 years of
    treatment: an open-label multicentre randomised phase 3 trial.**
    *Lancet Oncol* 2010.
    [PMID 20864406](https://pubmed.ncbi.nlm.nih.gov/20864406/)
    — [TEST] 2-year PFS 80% on continuation versus 16% on interruption. Years of
    deep response leave a reservoir; this is what the drug-insensitive quiescent
    compartment exists to represent.

41. Patrikidou A, Chabaud S, Ray-Coquard I, et al. **Influence of imatinib
    interruption and rechallenge on the residual disease in patients with
    advanced GIST: results of the BFR14 prospective French Sarcoma Group
    randomised phase III trial.** *Ann Oncol* 2013.
    [PMID 23175622](https://pubmed.ncbi.nlm.nih.gov/23175622/)
    — [TEST] re-starting imatinib after interruption re-establishes control in
    most patients: the regrowing cells are the *same sensitive* cells, not new
    resistant ones. Impossible in a single-clone drifting-potency model.

---

## 8 · Adjuvant therapy, risk stratification and surgery

42. Joensuu H, Eriksson M, Sundby Hall K, et al. **Survival outcomes associated
    with 3 years vs 1 year of adjuvant imatinib for patients with high-risk
    gastrointestinal stromal tumors: an analysis of a randomized clinical trial
    after 10-year follow-up.** *JAMA Oncol* 2020.
    [PMID 32469385](https://pubmed.ncbi.nlm.nih.gov/32469385/)
    — [ANCHOR] SSGXVIII: 5-year RFS 71.4% (36 months) versus 53.0% (12 months);
    10-year RFS 52.5% versus 41.8%. `N_MICRO` is fitted to the 12-month arm; the
    convergence of the two curves with time is the model's statement that
    adjuvant therapy delays rather than cures.

43. DeMatteo RP, Ballman KV, Antonescu CR, et al. **Adjuvant imatinib mesylate
    after resection of localised, primary gastrointestinal stromal tumour: a
    randomised, double-blind, placebo-controlled trial.** *Lancet* 2009.
    [PMID 19303137](https://pubmed.ncbi.nlm.nih.gov/19303137/)
    — [ANCHOR] ACOSOG Z9001: the recurrence hazard resumes when the drug stops.

44. Corless CL, Ballman KV, Antonescu CR, et al. **Pathologic and molecular
    features correlate with long-term outcome after adjuvant therapy of resected
    primary GI stromal tumor: the ACOSOG Z9001 trial.** *J Clin Oncol* 2014.
    [PMID 24638003](https://pubmed.ncbi.nlm.nih.gov/24638003/)
    — [STRUCT] adjuvant benefit is genotype-dependent, exactly as the
    metastatic-setting EC50 matrix predicts.

45. Joensuu H, Vehtari A, Riihimäki J, et al. **Risk of recurrence of
    gastrointestinal stromal tumour after surgery: an analysis of pooled
    population-based cohorts.** *Lancet Oncol* 2012.
    [PMID 22153892](https://pubmed.ncbi.nlm.nih.gov/22153892/)
    — [PARAM] size, mitotic index and site map onto the occult residual burden
    `N_MICRO`; the map's cluster 20.

46. Hølmebakk T, Bjerkehagen B, Hompland I, et al. **Recurrence-free survival
    after resection of gastric gastrointestinal stromal tumors classified
    according to a strict definition of tumor rupture: a population-based
    study.** *Ann Surg Oncol* 2018.
    [PMID 29435684](https://pubmed.ncbi.nlm.nih.gov/29435684/)
    — [PARAM] rupture as a large multiplier on `N_MICRO`.

47. Wang D, Zhang Q, Blanke CD, et al. **Phase II trial of neoadjuvant/adjuvant
    imatinib mesylate for advanced primary and metastatic/recurrent operable
    gastrointestinal stromal tumors: long-term follow-up results of Radiation
    Therapy Oncology Group 0132.** *Ann Surg Oncol* 2012.
    [PMID 22203182](https://pubmed.ncbi.nlm.nih.gov/22203182/)
    — [TEST] neoadjuvant downstaging: deep response without eradication, so
    surgery still has to remove the reservoir.

---

## 9 · Imaging: why PET and CT disagree  [C2]

48. Van den Abbeele AD, Gatsonis C, de Vries DJ, et al. **ACRIN 6665/RTOG 0132
    phase II trial of neoadjuvant imatinib mesylate for operable malignant
    gastrointestinal stromal tumor: monitoring with 18F-FDG PET and correlation
    with genotype and GLUT4 expression.** *J Nucl Med* 2012.
    [PMID 22381410](https://pubmed.ncbi.nlm.nih.gov/22381410/)
    — [PARAM] `TAUPET = 0.8 day`: the metabolic response is essentially
    immediate and is a read-out of signalling (glucose-transporter expression
    downstream of KIT), not of mass.

49. Van den Abbeele AD. **The lessons of GIST — PET and PET/CT: a new paradigm
    for imaging.** *Oncologist* 2008.
    [PMID 18434632](https://pubmed.ncbi.nlm.nih.gov/18434632/)
    — [STRUCT] the explicit separation of signalling, viable mass and imaged
    mass into three state variables with three time constants.

---

## 10 · Toxicity, organ systems and dose intensity

50. Jiang X, Zhang X, Wang Y, et al. **Hematologic toxicities of sunitinib in
    patients with gastrointestinal stromal tumors: a systematic review and
    meta-analysis.** *Int J Colorectal Dis* 2022.
    [PMID 35780257](https://pubmed.ncbi.nlm.nih.gov/35780257/)
    — [PARAM] `SLANCSU = 1.1e-3 mL/ng`, giving an ANC nadir of about 2.4 ×10⁹/L
    at 50 mg on a 4/2 schedule.

51. De Leo S, Trevisan M, Fugazzola L. **Endocrine-related adverse conditions
    induced by tyrosine kinase inhibitors.** *Ann Endocrinol (Paris)* 2023.
    [PMID 36963756](https://pubmed.ncbi.nlm.nih.gov/36963756/)
    — [PARAM] the sunitinib thyroid axis: capillary regression → free T4 falls →
    TSH rises (`EC50THY = 60 ng/mL`, `EMAXTHY = 0.80`).

52. Wang E, Bello CL, Kang YK, et al. **Population pharmacokinetics of sunitinib
    and its active metabolite SU012662 in pediatric patients with
    gastrointestinal stromal tumors or other solid tumors.**
    *Eur J Drug Metab Pharmacokinet* 2021.
    [PMID 33852135](https://pubmed.ncbi.nlm.nih.gov/33852135/)
    — [PARAM] sunitinib and SU12662 disposition (`CL_SU`, `V1_SU`, `FM_SU`).

53. Khosravan R, Motzer RJ, Fumagalli E, Rini BI. **Extrapolation of
    pharmacokinetics and pharmacodynamics of sunitinib in children with
    gastrointestinal stromal tumors.** *Cancer Chemother Pharmacol* 2021.
    [PMID 33507338](https://pubmed.ncbi.nlm.nih.gov/33507338/)
    — [PARAM] sunitinib exposure–toxicity (neutropenia, blood pressure) slopes.

54. Fukudo M, Asai K, Tani C, et al. **Pharmacokinetics of the oral multikinase
    inhibitor regorafenib and its association with real-world treatment
    outcomes.** *Invest New Drugs* 2021.
    [PMID 33830408](https://pubmed.ncbi.nlm.nih.gov/33830408/)
    — [PARAM] `CL_RE = 66 L/day`, giving a steady-state parent concentration of
    roughly 2400 ng/mL at 160 mg.

55. Fujita K, Ando Y, Kubota Y, et al. **Association between albumin-bilirubin
    grade and plasma trough concentrations of regorafenib and its metabolites
    M-2 and M-5.** *Invest New Drugs* 2024.
    [PMID 38517650](https://pubmed.ncbi.nlm.nih.gov/38517650/)
    — [PARAM] `FM_RE = 0.40` and the treatment of M-2/M-5 as equipotent, so the
    active exposure is roughly 3600 ng/mL of regorafenib equivalents.

---

## 11 · Guidelines and practice context

56. Casali PG, Blay JY, Abecassis N, et al. **Gastrointestinal stromal tumours:
    ESMO-EURACAN-GENTURIS Clinical Practice Guidelines for diagnosis, treatment
    and follow-up.** *Ann Oncol* 2022.
    [PMID 34560242](https://pubmed.ncbi.nlm.nih.gov/34560242/)
    — [STRUCT] the standard sequence the model's scenarios follow, and the
    genotype-directed exceptions (800 mg for exon 9, avapritinib for D842V).

57. Hirota S, Tateishi U, Nakamoto Y, et al. **English version of Japanese
    Clinical Practice Guidelines 2022 for gastrointestinal stromal tumor (GIST)
    issued by the Japan Society of Clinical Oncology.**
    *Int J Clin Oncol* 2024.
    [PMID 38609732](https://pubmed.ncbi.nlm.nih.gov/38609732/)
    — [STRUCT] response assessment intervals and the role of PET.

---

## What the model deliberately does not contain

- **HSP90 and other target-independent mechanisms** (pimitespib, ref 37): they
  act below the clone × drug EC50 matrix and would need a different layer.
- **Immune contribution.** KIT inhibition alters intratumoural macrophages and
  T cells; the map shows the edges (cluster 11) but the ODEs carry only a fixed
  immune term inside `KD0`.
- **Spatial structure.** Polyclonality is per-patient here, not per-nodule, so
  the model cannot represent focal progression treated by local therapy — which
  is exactly the setting where cytoreductive surgery of a single progressing
  nodule is used clinically.
- **Germline and pharmacogenetic modifiers** (ref 6), and hepatic/renal
  covariates on clearance beyond the single `CLF` multiplier.
- **Overall survival.** The model stops at progression-free survival per line;
  OS requires post-progression therapy assumptions that would not be identifiable
  from the data used here.
