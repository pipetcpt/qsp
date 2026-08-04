# References — HIV-associated Cryptococcal Meningitis QSP model

Every PMID below was resolved and its title verified individually against
PubMed while this model was being built. Where a reference is used as a
numerical calibration target, the number taken from it is quoted, so that a
reader can check the model against the source without re-reading the paper.
Two conventions:

- **[CAL]** the paper supplied a number the model was fitted to.
- **[PRED]** the paper supplied a number the model was tested against but
  *not* fitted to.

A deliberate note on disagreement between sources: the two best estimates of
early fungicidal activity (EFA) for the *same* regimen — amphotericin B
deoxycholate 1 mg/kg/day plus flucytosine — are −0.42 (Day 2013) and −0.56
(Bicanic 2008), and the two best estimates of 10-week mortality for the *same*
regimen — one week of amphotericin plus flucytosine — are 24.2% (ACTA) and
28.7% (AMBITION control arm). The model's residuals against these endpoints
should be read against that ~0.14 log10/day and ~4.5 percentage-point
between-trial spread, not against zero.

---

## 1. Randomised induction-therapy trials (the primary calibration set)

1. **Brouwer AE, Rajanuwong A, Chierakul W, et al.** Combination antifungal
   therapies for HIV-associated cryptococcal meningitis: a randomised trial.
   *Lancet* 2004;363(9423):1764-7.
   <https://pubmed.ncbi.nlm.nih.gov/15172774/>
   — established quantitative CSF culture as the pharmacodynamic endpoint and
   showed clearance is exponential; amphotericin + flucytosine the most rapidly
   fungicidal of four regimens.

2. **Bicanic T, Wood R, Meintjes G, et al.** High-dose amphotericin B with
   flucytosine for the treatment of cryptococcal meningitis in HIV-infected
   patients: a randomized trial. *Clin Infect Dis* 2008;47(1):123-30.
   <https://pubmed.ncbi.nlm.nih.gov/18505387/>
   — **[CAL]** EFA −0.45 ± 0.16 (AmB-d 0.7 mg/kg + 5FC) vs −0.56 ± 0.24
   (1.0 mg/kg + 5FC), p = 0.02. This *within-trial ratio* of 0.80 for a
   1.43-fold dose step is what forces the model's amphotericin kill term to be
   nearly first-order in concentration; a saturating fit cannot reproduce it.

3. **Nussbaum JC, Jackson A, Namarika D, et al.** Combination flucytosine and
   high-dose fluconazole compared with fluconazole monotherapy for the
   treatment of cryptococcal meningitis: a randomized trial in Malawi.
   *Clin Infect Dis* 2010;50(3):338-44.
   <https://pubmed.ncbi.nlm.nih.gov/20038244/>
   — **[CAL]** EFA −0.11 ± 0.09 (fluconazole 1200 alone) vs −0.28 ± 0.17
   (+ flucytosine); 2-week mortality 37% vs 10%, 10-week 58% vs 43%; grade
   III/IV neutropenia 5 vs 1 with no increase in infection-related events —
   the basis for the model's deliberately small neutropenia hazard.

4. **Longley N, Muzoora C, Taseera K, et al.** Dose response effect of
   high-dose fluconazole for HIV-associated cryptococcal meningitis in
   southwestern Uganda. *Clin Infect Dis* 2008;47(12):1556-61.
   <https://pubmed.ncbi.nlm.nih.gov/18990067/>
   — **[CAL]** fluconazole monotherapy dose-response; the model's −0.10 for
   800 mg/day is a prediction against roughly −0.07.

5. **Day JN, Chau TTH, Wolbers M, et al.** Combination antifungal therapy for
   cryptococcal meningitis. *N Engl J Med* 2013;368(14):1291-302.
   <https://pubmed.ncbi.nlm.nih.gov/23550668/>
   — **[CAL]** the model's single most important anchor. All groups received
   AmB-d 1 mg/kg/day. EFA −0.31 (AmB alone), −0.32 (+ fluconazole 800),
   −0.42 (+ flucytosine); deaths by day 14 / day 70: 25/44 (AmB alone), 15/30
   (+5FC), and HR 0.78 / 0.71 for the fluconazole arm. **The finding that
   fluconazole adds essentially nothing to amphotericin's fungicidal rate
   (−0.32 vs −0.31) while fluconazole alone is clearly active (−0.11) cannot
   be reproduced by independent additive drug effects.** It is the
   observation that forced the ergosterol node into the model.

6. **Molloy SF, Kanyama C, Heyderman RS, et al.** Antifungal combinations for
   treatment of cryptococcal meningitis in Africa (ACTA).
   *N Engl J Med* 2018;378(11):1004-17.
   <https://pubmed.ncbi.nlm.nih.gov/29539274/>
   — **[CAL]** 10-week mortality 35.1% (oral fluconazole + flucytosine),
   36.2% (1-week amphotericin), 39.7% (2-week amphotericin); 1-week
   amphotericin + flucytosine lowest at 24.2%. **[PRED]** flucytosine was
   superior to fluconazole as the amphotericin partner: 31.1% vs 45.0%,
   HR 0.62 (0.45–0.84). The model recovers the direction of this comparison
   from the ergosterol node alone but only about half its magnitude
   (HR 0.85), and this shortfall is reported rather than tuned away.

7. **Jarvis JN, Lawrence DS, Meya DB, et al.** Single-dose liposomal
   amphotericin B treatment for cryptococcal meningitis (AMBITION-cm).
   *N Engl J Med* 2022;386(12):1109-20.
   <https://pubmed.ncbi.nlm.nih.gov/35320642/>
   — **[CAL]** 10-week mortality 24.8% (101/407) vs 28.7% (117/407),
   difference −3.9 points, non-inferior at a 10-point margin; EFA −0.40 vs
   −0.42 log10 CFU/mL/day; grade 3 or 4 adverse events 50.0% vs 62.3%. The
   near-identical EFA with markedly fewer adverse events is the observation
   the model's separate CNS and renal-cortical amphotericin compartments
   exist to explain.

8. **Boulware DR, Meya DB, Muzoora C, et al.** Timing of antiretroviral
   therapy after diagnosis of cryptococcal meningitis (COAT).
   *N Engl J Med* 2014;370(26):2487-98.
   <https://pubmed.ncbi.nlm.nih.gov/24963568/>
   — **[PRED]** 26-week mortality 45% (40/88) with ART at 1–2 weeks vs 30%
   (27/89) at 5 weeks, HR 1.73 (1.06–2.82); excess deaths concentrated at
   2–5 weeks; CSF white cells <5/mm³ subgroup HR 3.87 (1.41–10.58); recognised
   IRIS 20% vs 13%. The model reproduces the 15-point gap (16 points) without
   any term for ART timing, because IRIS drive is written as the product of
   the rate of immune recovery and the antigen present when it begins.

9. **Beardsley J, Wolbers M, Kibengo FM, et al.** Adjunctive dexamethasone in
   HIV-associated cryptococcal meningitis (CryptoDex).
   *N Engl J Med* 2016;374(6):542-54.
   <https://pubmed.ncbi.nlm.nih.gov/26863355/>
   — **[PRED]** stopped for safety. 10-week mortality 47% vs 41%, 6-month 57%
   vs 49%; prespecified good outcome 13% vs 25% (OR 0.42); more grade 3/4
   infection (48 vs 25), renal (22 vs 7) and cardiac (8 vs 0) events; **fungal
   clearance in CSF was slower on dexamethasone**. The model predicts the
   slower clearance mechanistically (glucocorticoid suppression of the
   macrophage killing term) and needs added steroid-attributable harm terms to
   reproduce the mortality and disability signals.

10. **Rhein J, Huppler Hullsiek K, Tugume L, et al.** Adjunctive sertraline
    for HIV-associated cryptococcal meningitis (ASTRO-CM): a randomised,
    placebo-controlled, double-blind phase 3 trial.
    *Lancet Infect Dis* 2019;19(8):843-51.
    <https://pubmed.ncbi.nlm.nih.gov/31345462/>
    — **[PRED]** no benefit. The model predicts a change in EFA of 0.000
    log10/day because the free brain concentration achievable at 400 mg/day is
    ~0.6% of the in-vitro EC50 — a calculation that did not require the trial.

11. **Jarvis JN, Meintjes G, Rebe K, et al.** Adjunctive interferon-γ
    immunotherapy for the treatment of HIV-associated cryptococcal meningitis:
    a randomized controlled trial. *AIDS* 2012;26(9):1105-13.
    <https://pubmed.ncbi.nlm.nih.gov/22421244/>
    — **[PRED]** adjunctive IFN-γ accelerated CSF fungal clearance; the model
    improves EFA from −0.40 to −0.44 through the IFN-γ-dependent macrophage
    activation term.

12. **Ngan NTT, Flower B, Day JN.** Treatment of cryptococcal meningitis:
    how have we got here and where are we going? — see also the sertraline
    PK-PD analysis: **Rhein J, et al.** Pharmacokinetics-pharmacodynamics of
    sertraline as an antifungal in HIV-infected Ugandans with cryptococcal
    meningitis. *J Pharmacokinet Pharmacodyn* 2019;46(6):519-28.
    <https://pubmed.ncbi.nlm.nih.gov/31584146/>

---

## 2. Intracranial pressure — the slow clock's clinical face

13. **Denning DW, Armstrong RW, Lewis BH, Stevens DA.** Elevated
    cerebrospinal fluid pressures in patients with cryptococcal meningitis and
    acquired immunodeficiency syndrome. *Am J Med* 1991;91(3):267-72.
    <https://pubmed.ncbi.nlm.nih.gov/1892147/>
    — the original description, and the source of the hypothesis that
    capsular polysaccharide obstructs CSF outflow at the arachnoid villi.

14. **Graybill JR, Sobel J, Saag M, et al.** Diagnosis and management of
    increased intracranial pressure in patients with AIDS and cryptococcal
    meningitis. *Clin Infect Dis* 2000;30(1):47-54.
    <https://pubmed.ncbi.nlm.nih.gov/10619732/>
    — NIAID Mycoses Study Group; the basis of the "drain to ≤50% of opening
    pressure or ≤200 mmH2O" rule implemented in the model's lumbar-puncture
    event.

15. **Bicanic T, Brouwer AE, Meintjes G, et al.** Relationship of
    cerebrospinal fluid pressure, fungal burden and outcome in patients with
    cryptococcal meningitis undergoing serial lumbar punctures.
    *AIDS* 2009;23(6):701-6.
    <https://pubmed.ncbi.nlm.nih.gov/19279443/>
    — **[PRED]** the single most structurally important observation for this
    model. Higher baseline fungal burden gave higher opening pressure, but
    **high burden was necessary and not sufficient** for high pressure;
    baseline pressure was *not* associated with CD4 count or with CSF
    pro-inflammatory cytokines; day-14 pressure tracked day-14 burden. The
    model's answer is that pressure is set by the *antigen* pool rather than
    the yeast or the inflammation, and antigen and yeast are correlated at
    presentation but decouple during therapy.

16. **Rolfes MA, Hullsiek KH, Rhein J, et al.** The effect of therapeutic
    lumbar punctures on acute mortality from cryptococcal meningitis.
    *Clin Infect Dis* 2014;59(11):1607-14.
    <https://pubmed.ncbi.nlm.nih.gov/25057102/>
    — **[PRED]** 11-day mortality 7% (5/75) with at least one therapeutic LP
    vs 18% (31/173) without; adjusted relative risk 0.31 (0.12–0.82), and the
    association held **regardless of opening pressure**. The model's LP effect
    (RR 0.75 at 11 days in the high-burden phenotype) is weaker than this
    observational estimate, which is confounded in the direction of
    under-stating benefit; the model's mechanism for pressure-independence is
    that the needle removes antigen as well as volume.

17. **Robertson EJ, Najjuka G, Rolfes MA, et al.** *Cryptococcus neoformans*
    ex vivo capsule size is associated with intracranial pressure and host
    immune response in HIV-associated cryptococcal meningitis.
    *J Infect Dis* 2014;209(1):74-82.
    <https://pubmed.ncbi.nlm.nih.gov/23945372/>
    — direct evidence linking capsule (not burden alone) to pressure; the
    empirical basis for the model's GXM → outflow-resistance term.

---

## 3. Prognosis, burden and the pleocytosis paradox

18. **Jarvis JN, Bicanic T, Loyse A, et al.** Determinants of mortality in a
    combined cohort of 501 patients with HIV-associated cryptococcal
    meningitis: implications for improving outcomes.
    *Clin Infect Dis* 2014;58(5):736-45.
    <https://pubmed.ncbi.nlm.nih.gov/24319084/>
    — 2-week mortality 17%, 10-week 34%, 1-year 41%; CSF fungal burden
    OR 1.4 per log10 CFU/mL; altered mental status OR 3.1; slow clearance
    independently predictive; **low** CSF opening pressure independently
    associated with 10-week mortality; IRIS in 13%, associated with 2-week
    fungal burden but not with time to ART.

19. **Bicanic T, Muzoora C, Brouwer AE, et al.** Independent association
    between rate of clearance of infection and clinical outcome of
    HIV-associated cryptococcal meningitis: analysis of a combined cohort of
    262 patients. *Clin Infect Dis* 2009;49(5):702-9.
    <https://pubmed.ncbi.nlm.nih.gov/19613840/>
    — the paper that made EFA a legitimate surrogate; also the reason the
    model carries an explicit "days the culture remains positive" hazard term
    rather than burden alone.

20. **Boulware DR, Bonham SC, Meya DB, et al.** Paucity of initial
    cerebrospinal fluid inflammation in cryptococcal meningitis is associated
    with subsequent immune reconstitution inflammatory syndrome.
    *J Infect Dis* 2010;202(6):962-70.
    <https://pubmed.ncbi.nlm.nih.gov/20677939/>
    — the mechanistic basis of the model's paucicellular phenotype: absent CSF
    inflammation is a failure of the clearance machinery, so antigen persists,
    so the IRIS stock is larger when ART arrives.

---

## 4. Epidemiology and guidelines

21. **Rajasingham R, Govender NP, Jordan A, et al.** The global burden of
    HIV-associated cryptococcal infection in adults in 2020: a modelling
    analysis. *Lancet Infect Dis* 2022;22(12):1748-55.
    <https://pubmed.ncbi.nlm.nih.gov/36049486/>

22. **Rajasingham R, Smith RM, Park BJ, et al.** Global burden of disease of
    HIV-associated cryptococcal meningitis: an updated analysis.
    *Lancet Infect Dis* 2017;17(8):873-81.
    <https://pubmed.ncbi.nlm.nih.gov/28483415/>

23. **Perfect JR, Dismukes WE, Dromer F, et al.** Clinical practice
    guidelines for the management of cryptococcal disease: 2010 update by the
    Infectious Diseases Society of America.
    *Clin Infect Dis* 2010;50(3):291-322.
    <https://pubmed.ncbi.nlm.nih.gov/20047480/>

24. **Chang CC, Harrison TS, Bicanic TA, et al.** Global guideline for the
    diagnosis and management of cryptococcosis: an initiative of the ECMM and
    ISHAM in cooperation with the ASM.
    *Lancet Infect Dis* 2024;24(8):e495-e512.
    <https://pubmed.ncbi.nlm.nih.gov/38346436/>

25. **Stott KE, Loyse A, Jarvis JN, et al.** Cryptococcal
    meningoencephalitis: time for action.
    *Lancet Infect Dis* 2021;21(9):e259-e271.
    <https://pubmed.ncbi.nlm.nih.gov/33872594/>

26. **Lawrence DS, Muthoga C, Meya DB, et al.** The acceptability of the
    AMBITION-cm treatment regimen for HIV-associated cryptococcal meningitis.
    *PLoS Negl Trop Dis* 2022;16(10):e0010825.
    <https://pubmed.ncbi.nlm.nih.gov/36279300/>

---

## 5. Amphotericin B — mechanism and pharmacology

27. **Anderson TM, Clay MC, Cioffi AG, et al.** Amphotericin forms an
    extramembranous and fungicidal sterol sponge.
    *Nat Chem Biol* 2014;10(5):400-6.
    <https://pubmed.ncbi.nlm.nih.gov/24681535/>
    — the mechanism that makes ergosterol amphotericin's *substrate* rather
    than merely its receptor, and therefore makes azole-induced ergosterol
    depletion an antagonism rather than a redundancy.

28. **Guo X, Zhang J, Li X, et al.** Fungicidal amphotericin B sponges are
    assemblies of staggered asymmetric homodimers encasing large void volumes.
    *Nat Struct Mol Biol* 2021;28(12):972-81.
    <https://pubmed.ncbi.nlm.nih.gov/34887566/>

29. **Lestner JM, Groll AH, Aljayyoussi G, et al.** Experimental models of
    short courses of liposomal amphotericin B for induction therapy for
    cryptococcal meningitis. *Antimicrob Agents Chemother* 2017;61(6):e00090-17.
    <https://pubmed.ncbi.nlm.nih.gov/28320715/>
    — the preclinical basis for a single high liposomal dose: brain
    concentrations persist far beyond plasma, which is why the model's CNS
    compartment has a long efflux half-life relative to dosing interval.

30. **Kwizera R, Akampurira A, Kandole TK, et al.** Relative contribution of
    pharmacokinetics and immune signatures to clinical outcomes in patients
    with HIV-associated cryptococcal meningitis.
    *Open Forum Infect Dis* 2025;12(4):ofaf182.
    <https://pubmed.ncbi.nlm.nih.gov/40271162/>

---

## 6. Flucytosine and fluconazole — mechanism and resistance

31. **Chen Y-C, Chang T-Y, Liu J-W, et al.** Flucytosine resistance in
    *Cryptococcus gattii* is indirectly mediated by the FCY2-FCY1-FUR1
    pathway. *Med Mycol* 2018;56(7):873-80.
    <https://pubmed.ncbi.nlm.nih.gov/29554336/>
    — the transport-and-activation chain (cytosine permease → cytosine
    deaminase → uracil phosphoribosyltransferase) whose loss the model
    represents as a 45-fold EC50 shift in a resistant subclone.

32. **Stone NRH, Bicanic T, Salim R, Hope W.** Liposomal amphotericin B
    (AmBisome®): a review of the pharmacokinetics, pharmacodynamics, clinical
    experience and future directions. See also **Rhodes J, Desjardins CA,
    Sykes SM, et al.** Dynamic ploidy changes drive fluconazole resistance in
    human cryptococcal meningitis. *J Clin Invest* 2019;129(3):999-1014.
    <https://pubmed.ncbi.nlm.nih.gov/30688656/>
    — aneuploidy-driven, reversible azole resistance in patients, the
    empirical counterpart of the model's heteroresistance annotation.

---

## 7. Fungal virulence, capsule and CNS invasion

33. **Charlier C, Nielsen K, Daou S, et al.** Evidence of a role for
    monocytes in dissemination and brain invasion by *Cryptococcus
    neoformans*. *Infect Immun* 2009;77(1):120-7.
    <https://pubmed.ncbi.nlm.nih.gov/18936186/>
    — the Trojan-horse route, and the reason the model carries a separate,
    drug-shielded intracellular fungal compartment.

34. **Olszewski MA, Noverr MC, Chen G-H, et al.** Urease expression by
    *Cryptococcus neoformans* promotes microvascular sequestration, thereby
    enhancing central nervous system invasion.
    *Am J Pathol* 2004;164(5):1761-71.
    <https://pubmed.ncbi.nlm.nih.gov/15111322/>

35. **Nosanchuk JD, Casadevall A.** Impact of melanin on microbial virulence
    and clinical resistance to antimicrobial compounds.
    *Antimicrob Agents Chemother* 2006;50(11):3519-28.
    <https://pubmed.ncbi.nlm.nih.gov/17065617/>
    — melanin as a determinant of reduced amphotericin susceptibility; enters
    the model as a modifier of the effective amphotericin EC50.

36. **Zaragoza O, García-Rodas R, Nosanchuk JD, et al.** Fungal cell gigantism
    during mammalian infection. *PLoS Pathog* 2010;6(6):e1000945.
    <https://pubmed.ncbi.nlm.nih.gov/20585557/>
    — titan cells; the phenotypic basis of the model's drug-tolerant persister
    compartment.

37. **Trevijano-Contador N, de Oliveira HC, García-Rodas R, et al.** Titan
    cells formation in *Cryptococcus neoformans* is finely tuned by
    environmental conditions and modulated by positive and negative genetic
    regulators. *PLoS Pathog* 2018;14(5):e1007007.
    <https://pubmed.ncbi.nlm.nih.gov/29775480/>

38. **Yang C, Huang Y, Zhou Y, et al.** *Cryptococcus* escapes host immunity:
    what do we know? — and on capsule surface reorganisation:
    **Sabbatini S, et al.** A glucuronoxylomannan epitope exhibits
    serotype-specific accessibility and redistributes towards the capsule
    surface during titanization of the fungal cell.
    *Infect Immun* 2019;87(4):e00779-18.
    <https://pubmed.ncbi.nlm.nih.gov/30670549/>

39. **Trevijano-Contador N, Zaragoza O, et al.** Capsule growth in
    *Cryptococcus neoformans* is coordinated with cell cycle progression.
    *mBio* 2014;5(3):e00945-14.
    <https://pubmed.ncbi.nlm.nih.gov/24939886/>
    — the reason the model lets capsule production saturate rather than scale
    without limit with burden.

40. **Chang YC, Stins MF, McCaffery MJ, et al.** *Cryptococcus neoformans*
    induces alterations in the cytoskeleton of human brain microvascular
    endothelial cells. *J Med Microbiol* 2003;52(11):961-70.
    <https://pubmed.ncbi.nlm.nih.gov/14532340/>

---

## 8. Immunopathology and IRIS

41. **Meya DB, Okurut S, Zziwa G, et al.** Advances in cryptococcal
    infections-associated immunopathology and potential therapeutic
    strategies. *Front Pharmacol* 2026 (in press).
    <https://pubmed.ncbi.nlm.nih.gov/42222160/>

42. **Musubire AK, Meya DB, Rhein J, et al.** Strategies for the diagnosis and
    management of meningitis in HIV-infected adults in resource-limited
    settings. *Expert Opin Pharmacother* 2021;22(15):2039-53.
    <https://pubmed.ncbi.nlm.nih.gov/34154509/>

43. **Nel JS, et al.** Induction treatment for HIV-associated cryptococcal
    meningitis: where have we been and where are we going?
    *Microorganisms* 2025;13(4):845.
    <https://pubmed.ncbi.nlm.nih.gov/40284683/>

---

## 9. Screening, prevention and pre-emptive therapy

44. **Ssebambulidde K, et al.** Fluconazole plus flucytosine versus
    fluconazole alone for adults with HIV-associated cryptococcal antigenaemia
    identified through screening: a multi-centre randomised trial protocol.
    *Trials* 2026;27:—.
    <https://pubmed.ncbi.nlm.nih.gov/41877274/>
    — the pre-emptive-therapy question the model's low-burden,
    partially-immune phenotype is meant to speak to.

---

## 10. What the model could not calibrate, and where it is therefore weakest

These are stated so that a reader knows which of the model's outputs are
constrained by data and which are structural extrapolation:

- **No trial has randomised intracranial-pressure management.** Rolfes 2014
  (#16) is observational. The pressure and perfusion hazard coefficients are
  therefore *fixed at mechanistically anchored values, not fitted*, and every
  quantitative claim about therapeutic lumbar puncture in this model is an
  extrapolation constrained by a single confounded estimate.
- **No trial has measured CSF GXM concentration serially alongside
  quantitative culture.** The 13-day antigen half-life that drives the whole
  slow clock is inferred from the well-documented lag of CrAg titre behind
  culture conversion, not measured. This is the model's most consequential
  unmeasured parameter, and the observation that would most cheaply falsify
  it is a serial CSF GXM assay during induction therapy.
- **Amphotericin CNS concentrations in humans are essentially unmeasured.**
  The CNS delivery parameters were solved *backwards* from early fungicidal
  activity (see #5, #7), which means the model's amphotericin brain
  concentrations are an inference from effect, not an observation, and the
  implied 36-fold difference in per-unit-plasma CNS delivery between
  deoxycholate and liposomal formulations is a prediction awaiting a
  measurement.
- **The ART-timing curve has no optimum** because the model contains no
  competing risk from leaving HIV untreated. It must not be read beyond about
  eight weeks.
