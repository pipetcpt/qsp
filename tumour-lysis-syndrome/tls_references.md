# Tumour Lysis Syndrome — References

References for `tls_qsp_model.dot`, `tls_mrgsolve_model.R`, `tls_shiny_app.R`
and `tls_reference_check.py`. Grouped by the part of the model each one
supports, so that a reader auditing a specific parameter can find its source.

Where the model is **calibrated rather than derived** — the crystal nucleation
rate constants, the obstruction half-max crystal mass, the interstitial
calcium-phosphate deposition constant, the nephron loss/recovery rates, and all
three hazard functions — that is stated in `README.md` and in the model file.
The references below support the *structure* and the *measured* parameters, not
those calibrated ones.

---

## 1. Definition, epidemiology and consensus guidance

1. Cairo MS, Bishop M. **Tumour lysis syndrome: new therapeutic strategies and
   classification.** Br J Haematol. 2004;127(1):3–11. — the laboratory/clinical
   TLS definition implemented in `cairo_bishop()`.
   <https://pubmed.ncbi.nlm.nih.gov/15384972/>
2. Coiffier B, Altman A, Pui CH, Younes A, Cairo MS. **Guidelines for the
   management of pediatric and adult tumor lysis syndrome: an evidence-based
   review.** J Clin Oncol. 2008;26(16):2767–78. — risk stratification and the
   1.5–5% dialysis figure used in the trial ledger.
   <https://pubmed.ncbi.nlm.nih.gov/18509186/>
3. Cairo MS, Coiffier B, Reiter A, Younes A. **Recommendations for the
   evaluation of risk and prophylaxis of tumour lysis syndrome in adults and
   children with malignant diseases: an expert TLS panel consensus.**
   Br J Haematol. 2010;149(4):578–86.
   <https://pubmed.ncbi.nlm.nih.gov/20331465/>
4. Howard SC, Jones DP, Pui CH. **The tumor lysis syndrome.**
   N Engl J Med. 2011;364(19):1844–54. — the standard synthesis; source for the
   general shape of the phosphate and calcium limbs.
   <https://pubmed.ncbi.nlm.nih.gov/21561350/>
5. Jones GL, Will A, Jackson GH, Webb NJA, Rule S. **Guidelines for the
   management of tumour lysis syndrome in adults and children with haematological
   malignancies on behalf of the British Committee for Standards in
   Haematology.** Br J Haematol. 2015;169(5):661–71.
   <https://pubmed.ncbi.nlm.nih.gov/25876990/>
6. Wilson FP, Berns JS. **Onco-nephrology: tumor lysis syndrome.**
   Clin J Am Soc Nephrol. 2012;7(10):1730–9.
   <https://pubmed.ncbi.nlm.nih.gov/22879434/>
7. Darmon M, Vincent F, Camous L, et al. **Tumour lysis syndrome and acute
   kidney injury in high-risk haematology patients in the rasburicase era.**
   Br J Haematol. 2013;162(4):489–97. — that AKI persists after rasburicase,
   which is the observation the model's residual-limb analysis addresses.
   <https://pubmed.ncbi.nlm.nih.gov/23772757/>
8. Belay Y, Yirdaw K, Enawgaw B. **Tumor lysis syndrome in patients with
   hematological malignancies.** J Oncol. 2017;2017:9684909.
   <https://pubmed.ncbi.nlm.nih.gov/29230334/>
9. Barbar T, Jaffer Sathick I. **Tumor lysis syndrome.**
   Adv Chronic Kidney Dis. 2021;28(5):438–46.
   <https://pubmed.ncbi.nlm.nih.gov/35190110/>
10. Williams SM, Killeen AA. **Tumor lysis syndrome.**
    Arch Pathol Lab Med. 2019;143(3):386–93.
    <https://pubmed.ncbi.nlm.nih.gov/30499695/>

---

## 2. Uric acid solubility, speciation and crystallisation

The model's urate solubility law is `S = S_HU · (1 + 10^(pH − pKa))` with
pKa 5.75 and S_HU 0.655 mmol/L, which gives 130 mg/L at pH 5.0 against a
reported ~150 mg/L, and 2070 mg/L at pH 7.0 against a reported ~2000.

11. Kelton J, Kelley WN, Holmes EW. **A rapid method for the diagnosis of acute
    uric acid nephropathy.** Arch Intern Med. 1978;138(4):612–5.
    <https://pubmed.ncbi.nlm.nih.gov/637642/>
12. Wilcox WR, Khalaf A, Weinberger A, Kippen I, Klinenberg JR. **Solubility of
    uric acid and monosodium urate.** Med Biol Eng. 1972;10(4):522–31.
    <https://pubmed.ncbi.nlm.nih.gov/5074854/>
13. Grases F, Villacampa AI, Costa-Bauzá A, Söhnel O. **Uric acid
    calculi: types, etiology and mechanisms of formation.**
    Clin Chim Acta. 2000;302(1–2):89–104.
    <https://pubmed.ncbi.nlm.nih.gov/11074067/>
14. Ngo TC, Assimos DG. **Uric acid nephrolithiasis: recent progress and future
    directions.** Rev Urol. 2007;9(1):17–27.
    <https://pubmed.ncbi.nlm.nih.gov/17396168/>
15. Conger JD, Falk SA. **Intrarenal dynamics in the pathogenesis and
    prevention of acute urate nephropathy.** J Clin Invest. 1977;59(5):786–93.
    — the experimental basis for the obstruction limb, including the effect of
    urine flow on intratubular urate.
    <https://pubmed.ncbi.nlm.nih.gov/845256/>
16. Conger JD, Falk SA, Guggenheim SJ, Burke TJ. **A micropuncture study of the
    early phase of acute urate nephropathy.** J Clin Invest. 1976;58(3):681–9.
    <https://pubmed.ncbi.nlm.nih.gov/956395/>
17. Shimada M, Johnson RJ, May WS, et al. **A novel role for uric acid in acute
    kidney injury associated with tumour lysis syndrome.** Nephrol Dial
    Transplant. 2009;24(10):2960–4. — the crystal-INDEPENDENT soluble-urate
    injury limb (`URATETOX` in the model).
    <https://pubmed.ncbi.nlm.nih.gov/19581334/>
18. Ejaz AA, Mu W, Kang DH, et al. **Could uric acid have a role in acute renal
    failure?** Clin J Am Soc Nephrol. 2007;2(1):16–21.
    <https://pubmed.ncbi.nlm.nih.gov/17699383/>
19. Sanchez-Lozada LG, Tapia E, Santamaría J, et al. **Mild hyperuricemia
    induces vasoconstriction and maintains glomerular hypertension in normal and
    remnant kidney rats.** Kidney Int. 2005;67(1):237–47.
    <https://pubmed.ncbi.nlm.nih.gov/15610247/>
20. Mulay SR, Anders HJ. **Crystallopathies.**
    N Engl J Med. 2016;374(25):2465–76. — the general framework for
    crystal-driven tissue injury including NLRP3 activation.
    <https://pubmed.ncbi.nlm.nih.gov/27332905/>

---

## 3. Renal urate handling and its saturability

The model lets the fractional excretion of urate rise from 8% to ~30% as plasma
urate rises, because URAT1/GLUT9 reabsorption saturates. That is the kidney's
only defence against a purine load, and its price is a higher *tubular* urate
concentration — which is why the defence itself generates the crystal.

21. Enomoto A, Kimura H, Chairoungdua A, et al. **Molecular identification of a
    renal urate–anion exchanger that regulates blood urate levels.**
    Nature. 2002;417(6887):447–52. — URAT1 / SLC22A12.
    <https://pubmed.ncbi.nlm.nih.gov/12024214/>
22. Vitart V, Rudan I, Hayward C, et al. **SLC2A9 is a newly identified urate
    transporter influencing serum urate concentration, urate excretion and
    gout.** Nat Genet. 2008;40(4):437–42.
    <https://pubmed.ncbi.nlm.nih.gov/18327257/>
23. Woodward OM, Köttgen A, Coresh J, et al. **Identification of a urate
    transporter, ABCG2, with a common functional polymorphism causing gout.**
    Proc Natl Acad Sci USA. 2009;106(25):10338–42.
    <https://pubmed.ncbi.nlm.nih.gov/19506252/>
24. Mandal AK, Mount DB. **The molecular physiology of uric acid homeostasis.**
    Annu Rev Physiol. 2015;77:323–45.
    <https://pubmed.ncbi.nlm.nih.gov/25422986/>
25. Maesaka JK, Fishbane S. **Regulation of renal urate excretion: a critical
    review.** Am J Kidney Dis. 1998;32(6):917–33.
    <https://pubmed.ncbi.nlm.nih.gov/9856507/>
26. Sorensen LB. **Role of the intestinal tract in the elimination of uric
    acid.** Arthritis Rheum. 1965;8(5):694–706. — the extrarenal route, ~30% of
    disposal.
    <https://pubmed.ncbi.nlm.nih.gov/5859543/>
27. Sorensen LB. **Degradation of uric acid in man.**
    Metabolism. 1959;8:687–703. — the miscible urate pool (~1.2 g) that the
    lead-time analysis compares against the released load.
    <https://pubmed.ncbi.nlm.nih.gov/13820342/>

---

## 4. Purine catabolism, xanthine oxidase and the xanthine substitution

28. Hille R. **Molybdenum-containing hydroxylases.**
    Arch Biochem Biophys. 2005;433(1):107–16.
    <https://pubmed.ncbi.nlm.nih.gov/15581570/>
29. Harkness RA. **Hypoxanthine, xanthine and uridine in body fluids,
    indicators of ATP depletion.** J Chromatogr. 1988;429:255–78.
    <https://pubmed.ncbi.nlm.nih.gov/3062147/>
30. Band PR, Silverberg DS, Henderson JF, et al. **Xanthine nephropathy in a
    patient with lymphosarcoma treated with allopurinol.**
    N Engl J Med. 1970;283(7):354–7. — the index case for the substitution the
    flux operator makes.
    <https://pubmed.ncbi.nlm.nih.gov/5434139/>
31. LaRosa C, McMullen L, Bakdash S, et al. **Acute renal failure from
    xanthine nephropathy during management of acute leukemia.**
    Pediatr Nephrol. 2007;22(1):132–5.
    <https://pubmed.ncbi.nlm.nih.gov/17024408/>
32. Ichida K, Amaya Y, Kamatani N, Nishino T, Hosoya T, Sakai O. **Identification
    of two mutations in human xanthine dehydrogenase gene responsible for
    classical type I xanthinuria.** J Clin Invest. 1997;99(10):2391–7. — the
    natural experiment for complete XO loss.
    <https://pubmed.ncbi.nlm.nih.gov/9153281/>
33. Gómez-Arbelaez D, et al. **Xanthinuria and xanthine stones.**
    Urolithiasis. 2019;47(3):215–22.
    <https://pubmed.ncbi.nlm.nih.gov/30542938/>

---

## 5. Allopurinol and febuxostat pharmacology (the FLUX operators)

Model values: allopurinol V 90 L, CL 45 L/h (t½ 1.4 h), 75% converted to
oxypurinol; oxypurinol V 40 L, CL 1.2 L/h **scaling with GFR** (t½ 23 h),
XO IC50 1.1 mg/L; febuxostat V 50 L, CL 5 L/h (t½ 6.9 h), hepatic, IC50
0.030 mg/L.

34. Day RO, Graham GG, Hicks M, McLachlan AJ, Stocker SL, Williams KM.
    **Clinical pharmacokinetics and pharmacodynamics of allopurinol and
    oxypurinol.** Clin Pharmacokinet. 2007;46(8):623–44.
    <https://pubmed.ncbi.nlm.nih.gov/17655371/>
35. Turnheim K, Krivanek P, Oberbauer R. **Pharmacokinetics and
    pharmacodynamics of allopurinol in elderly and young subjects.**
    Br J Clin Pharmacol. 1999;48(4):501–9.
    <https://pubmed.ncbi.nlm.nih.gov/10583019/>
36. Graham S, Day RO, Wong H, et al. **Pharmacodynamics of oxypurinol after
    administration of allopurinol to healthy subjects.**
    Br J Clin Pharmacol. 1996;41(4):299–304.
    <https://pubmed.ncbi.nlm.nih.gov/8730975/>
37. Hande KR, Noone RM, Stone WJ. **Severe allopurinol toxicity: description
    and guidelines for prevention in patients with renal insufficiency.**
    Am J Med. 1984;76(1):47–56. — the clinical counterpart of the oxypurinol
    accumulation the model predicts in evolving AKI.
    <https://pubmed.ncbi.nlm.nih.gov/6691361/>
38. Hung SI, Chung WH, Liou LB, et al. **HLA-B*5801 allele as a genetic marker
    for severe cutaneous adverse reactions caused by allopurinol.**
    Proc Natl Acad Sci USA. 2005;102(11):4134–9.
    <https://pubmed.ncbi.nlm.nih.gov/15743917/>
39. Khosravan R, Grabowski BA, Mayer MD, Wu JT, Joseph-Ridge N, Vernillet L.
    **The effect of mild and moderate renal impairment on the pharmacokinetics,
    pharmacodynamics and safety of febuxostat.**
    J Clin Pharmacol. 2008;48(9):1014–24. — the basis for febuxostat needing no
    renal dose adjustment.
    <https://pubmed.ncbi.nlm.nih.gov/18635756/>
40. Spina M, Nagy Z, Ribera JM, et al. **FLORENCE: a randomized, double-blind,
    phase III pivotal study of febuxostat versus allopurinol for the prophylaxis
    of tumor lysis syndrome in patients with haematologic malignancies at
    intermediate to high TLS risk.** Ann Oncol. 2015;26(10):2155–61.
    <https://pubmed.ncbi.nlm.nih.gov/26216382/>
41. Becker MA, Schumacher HR, Wortmann RL, et al. **Febuxostat compared with
    allopurinol in patients with hyperuricemia and gout.**
    N Engl J Med. 2005;353(23):2450–61.
    <https://pubmed.ncbi.nlm.nih.gov/16339094/>
42. Elion GB. **The purine path to chemotherapy.**
    Science. 1989;244(4900):41–7. — allopurinol's origin and the 6-MP
    interaction.
    <https://pubmed.ncbi.nlm.nih.gov/2649979/>

---

## 6. Rasburicase (the POOL operator) and its oxidant limb

Model values: V1 8 L, CL 0.30 L/h (t½ 18 h), Vmax 1.50 mmol/h per mg/L,
Km 0.030 mmol/L — two orders of magnitude below TLS urate, hence the zero-order
behaviour that makes a dose buy a fixed mmol/h.

43. Goldman SC, Holcenberg JS, Finklestein JZ, et al. **A randomized comparison
    between rasburicase and allopurinol in children with lymphoma or leukemia at
    high risk for tumor lysis.** Blood. 2001;97(10):2998–3003. — the −86% at
    4 h and the 2.6-fold AUC difference used in the trial ledger.
    <https://pubmed.ncbi.nlm.nih.gov/11342423/>
44. Pui CH, Mahmoud HH, Wiley JM, et al. **Recombinant urate oxidase for the
    prophylaxis or treatment of hyperuricemia in patients with leukemia or
    lymphoma.** J Clin Oncol. 2001;19(3):697–704.
    <https://pubmed.ncbi.nlm.nih.gov/11157020/>
45. Cortes J, Moore JO, Maziarz RT, et al. **Control of plasma uric acid in
    adults at risk for tumor lysis syndrome: efficacy and safety of rasburicase
    alone and rasburicase followed by allopurinol compared with allopurinol
    alone.** J Clin Oncol. 2010;28(27):4207–13. — the 87% vs 66% response rates
    in the ledger.
    <https://pubmed.ncbi.nlm.nih.gov/20713865/>
46. Vadhan-Raj S, Fayad LE, Fanale MA, et al. **A randomized trial of a
    single-dose rasburicase versus five-daily doses in patients at risk for
    tumor lysis syndrome.** Ann Oncol. 2012;23(6):1640–5. — the dose/schedule
    question the capacity analysis addresses.
    <https://pubmed.ncbi.nlm.nih.gov/22015451/>
47. Ishizawa K, Ogura M, Hamaguchi M, et al. **Safety and efficacy of
    rasburicase (SR29142) in a Japanese phase II study.**
    Cancer Sci. 2009;100(5):357–62.
    <https://pubmed.ncbi.nlm.nih.gov/19154407/>
48. Browning LA, Kruse JA. **Hemolysis and methemoglobinemia secondary to
    rasburicase administration.** Ann Pharmacother. 2005;39(11):1932–5.
    <https://pubmed.ncbi.nlm.nih.gov/16174785/>
49. Sonbol MB, Yadav H, Vaidya R, Rana V, Witzig TE. **Methemoglobinemia and
    hemolysis in a patient with G6PD deficiency treated with rasburicase.**
    Am J Hematol. 2013;88(2):152–4.
    <https://pubmed.ncbi.nlm.nih.gov/22886583/>
50. Cheah CY, Lew TE, Seymour JF, Burbury K. **Rasburicase causing severe
    oxidative hemolysis and methemoglobinemia in a patient with previously
    unrecognized glucose-6-phosphate dehydrogenase deficiency.**
    Acta Haematol. 2013;130(4):254–9.
    <https://pubmed.ncbi.nlm.nih.gov/23860388/>
51. Lopez-Olivo MA, Pratt G, Palla SL, Salahudeen A. **Rasburicase in tumor
    lysis syndrome of the adult: a systematic review and meta-analysis.**
    Am J Kidney Dis. 2013;62(3):481–92. — that the benefit on hard renal
    endpoints is not established, which is the observation the model's
    residual-limb and capacity analyses are about.
    <https://pubmed.ncbi.nlm.nih.gov/23684124/>
52. Cammalleri L, Malaguarnera M. **Rasburicase represents a new tool for
    hyperuricemia in tumor lysis syndrome and in gout.**
    Int J Med Sci. 2007;4(2):83–93.
    <https://pubmed.ncbi.nlm.nih.gov/17396159/>
53. Hummel M, Reiter S, Adam K, Hehlmann R, Buchheidt D. **Effective treatment
    and prophylaxis of hyperuricemia and impaired renal function in tumor lysis
    syndrome with low doses of rasburicase.**
    Eur J Haematol. 2008;80(4):331–6.
    <https://pubmed.ncbi.nlm.nih.gov/18221387/>

---

## 7. Potassium: release, redistribution and removal

54. Allon M, Copkney C. **Albuterol and insulin for treatment of hyperkalemia in
    hemodialysis patients.** Kidney Int. 1990;38(5):869–72. — the −0.6 to
    −1.0 mmol/L shift magnitude used in the ledger.
    <https://pubmed.ncbi.nlm.nih.gov/2266671/>
55. Sterns RH, Grieff M, Bernstein PL. **Treatment of hyperkalemia: something
    old, something new.** Kidney Int. 2016;89(3):546–54. — the distinction
    between moving and removing potassium, which the model's total-body-K
    readout is built to make explicit.
    <https://pubmed.ncbi.nlm.nih.gov/26880451/>
56. Palmer BF, Clegg DJ. **Physiology and pathophysiology of potassium
    homeostasis.** Adv Physiol Educ. 2016;40(4):480–90.
    <https://pubmed.ncbi.nlm.nih.gov/27756725/>
57. Palmer BF. **Regulation of potassium homeostasis.**
    Clin J Am Soc Nephrol. 2015;10(6):1050–60. — flow-dependent distal
    secretion, ROMK/BK, and the aldosterone gain used in `k_capacity()`.
    <https://pubmed.ncbi.nlm.nih.gov/25287933/>
58. Sorensen MV, Matos JE, Praetorius HA, Leipziger J. **Colonic potassium
    handling.** Pflugers Arch. 2010;459(5):645–56. — the GFR-independent
    colonic route that survives anuria.
    <https://pubmed.ncbi.nlm.nih.gov/20143061/>
59. Packham DK, Rasmussen HS, Lavin PT, et al. **Sodium zirconium cyclosilicate
    in hyperkalemia.** N Engl J Med. 2015;372(3):222–31.
    <https://pubmed.ncbi.nlm.nih.gov/25415807/>
60. Weisberg LS. **Management of severe hyperkalemia.**
    Crit Care Med. 2008;36(12):3246–51.
    <https://pubmed.ncbi.nlm.nih.gov/18936701/>
61. Parham WA, Mehdirad AA, Biermann KM, Fredman CS. **Hyperkalemia revisited.**
    Tex Heart Inst J. 2006;33(1):40–7. — the ECG progression in the cardiac
    cluster of the map.
    <https://pubmed.ncbi.nlm.nih.gov/16572868/>
62. Montford JR, Linas S. **How dangerous is hyperkalemia?**
    J Am Soc Nephrol. 2017;28(11):3155–65.
    <https://pubmed.ncbi.nlm.nih.gov/28778861/>

---

## 8. Phosphate, calcium and calcium-phosphate precipitation

Model values: TmP/GFR ≈ 0.90 mmol/L with FGF23-like and PTH-driven suppression;
basal phosphate excretion 23 mmol/day; systemic deposition threshold set at
Ca × PO₄ = 4.84 mmol²/L², i.e. the classical 60 mg²/dL².

63. Murer H, Hernando N, Forster I, Biber J. **Proximal tubular phosphate
    reabsorption: molecular mechanisms.** Physiol Rev. 2000;80(4):1373–409.
    <https://pubmed.ncbi.nlm.nih.gov/11015617/>
64. Payne RB. **Renal tubular reabsorption of phosphate (TmP/GFR): indications
    and interpretation.** Ann Clin Biochem. 1998;35(2):201–6. — the saturable
    reabsorption maximum implemented as `TMP0`.
    <https://pubmed.ncbi.nlm.nih.gov/9547891/>
65. Shimada T, Hasegawa H, Yamazaki Y, et al. **FGF-23 is a potent regulator of
    vitamin D metabolism and phosphate homeostasis.**
    J Bone Miner Res. 2004;19(3):429–35.
    <https://pubmed.ncbi.nlm.nih.gov/15040831/>
66. Kestenbaum B, Sampson JN, Rudser KD, et al. **Serum phosphate levels and
    mortality risk among people with chronic kidney disease.**
    J Am Soc Nephrol. 2005;16(2):520–8.
    <https://pubmed.ncbi.nlm.nih.gov/15615819/>
67. Boyce BF, Rosenberg E, de Papp AE, Duong LT. **The osteoclast, bone
    remodelling and treatment of metabolic bone disease.**
    Eur J Clin Invest. 2012;42(12):1332–41.
    <https://pubmed.ncbi.nlm.nih.gov/22998735/>
68. Brown EM, MacLeod RJ. **Extracellular calcium sensing and extracellular
    calcium signaling.** Physiol Rev. 2001;81(1):239–97. — the CaSR set point
    implemented in the PTH equation.
    <https://pubmed.ncbi.nlm.nih.gov/11152759/>
69. Baird GS. **Ionized calcium.** Clin Chim Acta. 2011;412(9–10):696–701. —
    why total calcium is the wrong measurement in an alkalotic, hypoalbuminaemic
    patient.
    <https://pubmed.ncbi.nlm.nih.gov/21238450/>
70. Kuo IY, Ehrlich BE. **Signaling in muscle contraction.**
    Cold Spring Harb Perspect Biol. 2015;7(2):a006023.
    <https://pubmed.ncbi.nlm.nih.gov/25646377/>
71. Giachelli CM. **Ectopic calcification: gathering hard facts about soft
    tissue mineralization.** Am J Pathol. 1999;154(3):671–5.
    <https://pubmed.ncbi.nlm.nih.gov/10079244/>
72. O'Neill WC. **The fallacy of the calcium-phosphorus product.**
    Kidney Int. 2007;71(8):792–6. — an important caveat on the threshold the
    model uses; read alongside the calibration disclosure in README.md.
    <https://pubmed.ncbi.nlm.nih.gov/17311068/>
73. Evan AP. **Physiopathology and etiology of stone formation in the kidney and
    the urinary tract.** Pediatr Nephrol. 2010;25(5):831–41. — Randall's plaque
    and the pH dependence of calcium-phosphate deposition.
    <https://pubmed.ncbi.nlm.nih.gov/19924445/>
74. Coe FL, Evan AP, Worcester EM, Lingeman JE. **Three pathways for human
    kidney stone formation.** Urol Res. 2010;38(3):147–60.
    <https://pubmed.ncbi.nlm.nih.gov/20411383/>
75. Boonstra AH, Jackson WPU. **Serum calcium survey for hypercalcaemia.**
    Am J Clin Pathol. 1971;55(5):523–6.
    <https://pubmed.ncbi.nlm.nih.gov/5573450/>

---

## 9. Alkalinisation: the trade the model prices on both sides

76. Conger JD. **Acute uric acid nephropathy.**
    Med Clin North Am. 1990;74(4):859–71. — includes the argument for and
    against urinary alkalinisation.
    <https://pubmed.ncbi.nlm.nih.gov/2195257/>
77. Ten Harkel AD, Kist-Van Holthe JE, Van Weel M, Van der Vorst MM.
    **Alkalinization and the prevention of tumour lysis syndrome.**
    Med Pediatr Oncol. 1998;31(1):27–8.
    <https://pubmed.ncbi.nlm.nih.gov/9611927/>
78. Pais VM Jr, Lowe G, Lallas CD, Preminger GM, Assimos DG. **Xanthine
    urolithiasis.** Urology. 2006;67(5):1084–7. — why alkali cannot rescue the
    xanthine limb (pKa 7.7).
    <https://pubmed.ncbi.nlm.nih.gov/16698365/>
79. Pak CY, Sakhaee K, Fuller C. **Successful management of uric acid
    nephrolithiasis with potassium citrate.** Kidney Int. 1986;30(3):422–8.
    <https://pubmed.ncbi.nlm.nih.gov/3784284/>
80. Rodman JS. **Prophylaxis of uric acid stones with alternate day doses of
    alkaline potassium salts.** J Urol. 1991;145(1):97–9.
    <https://pubmed.ncbi.nlm.nih.gov/1985735/>

---

## 10. Venetoclax, dose-ramp design and flux shaping

The venetoclax dose ramp exists because of TLS, which makes it the cleanest
clinical example of a flux-shaping operator: the same molecule and the same
eventual kill, with the release rate as the manipulated variable.

81. Roberts AW, Davids MS, Pagel JM, et al. **Targeting BCL2 with venetoclax in
    relapsed chronic lymphocytic leukemia.** N Engl J Med. 2016;374(4):311–22. —
    the dose-finding experience that produced the ramp.
    <https://pubmed.ncbi.nlm.nih.gov/26639348/>
82. Davids MS, Hallek M, Wierda W, et al. **Comprehensive safety analysis of
    venetoclax monotherapy for patients with relapsed/refractory chronic
    lymphocytic leukemia.** Clin Cancer Res. 2018;24(18):4371–9.
    <https://pubmed.ncbi.nlm.nih.gov/29895707/>
83. Salem AH, Agarwal SK, Dunbar M, Enschede SL, Humerickhouse RA, Wong SL.
    **Pharmacokinetics of venetoclax, a novel BCL-2 inhibitor, in patients with
    relapsed or refractory chronic lymphocytic leukemia or non-Hodgkin
    lymphoma.** J Clin Pharmacol. 2017;57(4):484–92. — the PK values used.
    <https://pubmed.ncbi.nlm.nih.gov/27558232/>
84. Seymour JF, Kipps TJ, Eichhorst B, et al. **Venetoclax–rituximab in
    relapsed or refractory chronic lymphocytic leukemia.**
    N Engl J Med. 2018;378(12):1107–20.
    <https://pubmed.ncbi.nlm.nih.gov/29562156/>
85. DiNardo CD, Jonas BA, Pullarkat V, et al. **Azacitidine and venetoclax in
    previously untreated acute myeloid leukemia.**
    N Engl J Med. 2020;383(7):617–29.
    <https://pubmed.ncbi.nlm.nih.gov/32786187/>
86. Patel VK, Lamothe B, Ayres ML, et al. **Pharmacodynamics and
    pharmacokinetics of venetoclax in CLL.**
    Br J Haematol. 2018;180(4):547–56.
    <https://pubmed.ncbi.nlm.nih.gov/29271470/>

---

## 11. Cytoreduction kinetics, tumour burden and the steroid prephase

87. Magrath I. **The pathogenesis of Burkitt's lymphoma.**
    Adv Cancer Res. 1990;55:133–270. — doubling-time and burden figures for the
    model's index disease.
    <https://pubmed.ncbi.nlm.nih.gov/2166998/>
88. Patte C, Auperin A, Michon J, et al. **The Société Française d'Oncologie
    Pédiatrique LMB89 protocol: highly effective multiagent chemotherapy tailored
    to the tumour burden and initial response in 561 unselected children with
    B-cell lymphomas and L3 leukemia.** Blood. 2001;97(11):3370–9. — the
    cytoreductive prephase in practice.
    <https://pubmed.ncbi.nlm.nih.gov/11369626/>
89. Dunleavy K, Pittaluga S, Shovlin M, et al. **Low-intensity therapy in
    adults with Burkitt lymphoma.** N Engl J Med. 2013;369(20):1915–25.
    <https://pubmed.ncbi.nlm.nih.gov/24224624/>
90. Truong TH, Beyene J, Hitzler J, et al. **Features at presentation predict
    children with acute lymphoblastic leukemia at low risk for tumor lysis
    syndrome.** Cancer. 2007;110(8):1832–9.
    <https://pubmed.ncbi.nlm.nih.gov/17724692/>
91. Montesinos P, Lorenzo I, Martín G, et al. **Tumor lysis syndrome in
    patients with acute myeloid leukemia: identification of risk factors and
    development of a predictive model.** Haematologica. 2008;93(1):67–74.
    <https://pubmed.ncbi.nlm.nih.gov/18166787/>
92. Mirrakhimov AE, Voore P, Khan M, Ali AM. **Tumor lysis syndrome: A clinical
    review.** World J Crit Care Med. 2015;4(2):130–8.
    <https://pubmed.ncbi.nlm.nih.gov/25938028/>
93. Alberts DS, Chen HS, Soehnlen B, Salmon SE. **Cell content and volume of
    human tumour cells.** Cancer Chemother Pharmacol. 1980;4:157–63. — the cell
    volume used to convert intracellular concentrations to per-10¹² content.
    <https://pubmed.ncbi.nlm.nih.gov/7407700/>

---

## 12. Renal replacement therapy and fluid management

94. Ronco C, Bellomo R, Kellum JA. **Acute kidney injury.**
    Lancet. 2019;394(10212):1949–64.
    <https://pubmed.ncbi.nlm.nih.gov/31777389/>
95. Pichette V, Leblanc M, Bonnardeaux A, Ouimet D, Geadah D, Cardinal J.
    **High dialysate flow rate continuous arteriovenous hemodialysis: a new
    approach for the treatment of acute renal failure and tumor lysis
    syndrome.** Am J Kidney Dis. 1994;23(4):591–6. — phosphate clearance and
    rebound.
    <https://pubmed.ncbi.nlm.nih.gov/8154498/>
96. Tan HK, Bellomo R, M'Pisi DA, Ronco C. **Phosphatemic control during acute
    renal failure: intermittent hemodialysis versus continuous hemodiafiltration.**
    Int J Artif Organs. 2001;24(4):186–91.
    <https://pubmed.ncbi.nlm.nih.gov/11394699/>
97. Ho VW, Wu VC, Chen YM, et al. **Continuous renal replacement therapy for
    tumor lysis syndrome.** Nephrology. 2016;21(9):739–46.
    <https://pubmed.ncbi.nlm.nih.gov/26582072/>
98. Ho PJ, Rosner MH. **Tumor lysis syndrome: prevention and treatment.**
    Adv Chronic Kidney Dis. 2022;29(2):162–7.
    <https://pubmed.ncbi.nlm.nih.gov/35817521/>
99. Rose BD, Post TW. **Clinical physiology of acid–base and electrolyte
    disorders.** 5th ed. McGraw-Hill; 2001. — the free-water handling and
    maximal fractional water excretion used to bound urine flow.
100. Ho KM, Sheridan DJ. **Meta-analysis of frusemide to prevent or treat acute
     renal failure.** BMJ. 2006;333(7565):420–5. — why furosemide is a
     dilution operator only in the volume-replete patient, and harmful
     otherwise (the model reproduces the harm through its volume state).
     <https://pubmed.ncbi.nlm.nih.gov/16861256/>

---

## 13. Cytokine release syndrome and the differential

101. Lee DW, Santomasso BD, Locke FL, et al. **ASTCT consensus grading for
     cytokine release syndrome and neurologic toxicity associated with immune
     effector cells.** Biol Blood Marrow Transplant. 2019;25(4):625–38.
     <https://pubmed.ncbi.nlm.nih.gov/30592986/>
102. Neelapu SS, Tummala S, Kebriaei P, et al. **Chimeric antigen receptor
     T-cell therapy — assessment and management of toxicities.**
     Nat Rev Clin Oncol. 2018;15(1):47–62.
     <https://pubmed.ncbi.nlm.nih.gov/28925994/>
103. Howard SC, Trifilio S, Gregory TK, Baxter N, McBride A. **Tumor lysis
     syndrome in the era of novel and targeted agents in patients with
     hematologic malignancies: a systematic review.**
     Ann Hematol. 2016;95(4):563–73.
     <https://pubmed.ncbi.nlm.nih.gov/26758269/>

---

## 14. QSP methodology

104. Nijsen MJMA, Wu F, Bansal L, et al. **Preclinical QSP modeling in the
     pharmaceutical industry: an IQ consortium survey examining the current
     landscape.** CPT Pharmacometrics Syst Pharmacol. 2018;7(3):135–46.
     <https://pubmed.ncbi.nlm.nih.gov/29344814/>
105. Baron KT, Gastonguay MR. **mrgsolve: simulate from ODE-based population
     PK/PD and systems pharmacology models.** — the simulation engine.
     <https://mrgsolve.org/>
106. Musante CJ, Ramanujan S, Schmidt BJ, Ghobrial OG, Lu J, Heatherington AC.
     **Quantitative systems pharmacology: a case for disease models.**
     Clin Pharmacol Ther. 2017;101(1):24–7.
     <https://pubmed.ncbi.nlm.nih.gov/27709613/>
107. Cheng Y, Straube R, Alnaif AE, Huang L, Leil TA, Schmidt BJ. **Virtual
     populations for quantitative systems pharmacology models.**
     Methods Mol Biol. 2022;2486:129–79.
     <https://pubmed.ncbi.nlm.nih.gov/35437723/>

---

## Notes on provenance and on what this model does not have

- **The two threshold concentrations** (`UA_req`, `UA_crit`) are computed from
  the solubility law (refs 11–14), the saturable urate handling (refs 21–27) and
  the urine-flow prescription. They are not fitted to the Cairo–Bishop
  criterion; that the 8 mg/dL criterion falls inside the computed UA_crit range
  is a result, not an input.
- **The crystal nucleation constants, obstruction half-max, nephron
  loss/recovery rates, interstitial deposition constant and all three hazard
  functions are calibrated**, not measured. They were set so that an
  unprophylaxed high-burden Burkitt patient reaches a peak creatinine ratio of
  about 3 with recovery over two weeks. Absolute event probabilities from the
  hazard functions should not be quoted.
- **The model over-separates the rasburicase-versus-allopurinol urate AUC**
  (0.06× versus a reported 0.39×, ref 43). The likely reason is that the trial
  population was mixed-risk while the model's parameterisation is built around
  the extreme high-burden patient; it is recorded as a mismatch rather than
  tuned away.
- **No spatial structure.** Crystal deposition is a lumped renal mass, so the
  model cannot say anything about medullary versus cortical distribution, and
  the medullary concentrating factor (`CF_MED` = 1.4) is a single lumped number
  standing in for a whole axial gradient.
- **No coagulation, no sepsis, no drug-specific nephrotoxicity** (methotrexate,
  aminoglycosides, contrast), all of which co-occur with TLS in real patients
  and all of which would add to the injury drive.
