# Tobacco Use Disorder (Nicotine Dependence) — References

> Curated, PubMed-indexed reference list underpinning the **Tobacco Use
> Disorder / Nicotine Dependence** QSP model (`tud_qsp_model.dot`,
> `tud_mrgsolve_model.R`, `tud_shiny_app.R`).
>
> Every PMID below was resolved against the live PubMed E-utilities API and the
> returned title/journal/year is reproduced verbatim, so each link is a real
> record rather than a reconstructed identifier. References are grouped by the
> part of the model they support, and the **"→ model"** annotations name the
> specific parameter, equation or calibration target each citation anchors.

---

## 1. Burden of disease, dependence constructs & natural history

1. GBD 2019 Risk Factors Collaborators. **Global burden of 87 risk factors in 204 countries and territories, 1990–2019: a systematic analysis for the Global Burden of Disease Study 2019.** Lancet. 2020. <https://pubmed.ncbi.nlm.nih.gov/33069327/>
2. Benowitz NL. **Nicotine addiction.** N Engl J Med. 2010;362(24):2295–2303. <https://pubmed.ncbi.nlm.nih.gov/20554984/> — → model: overall disease framing; the review that motivates modelling dependence as a receptor-adaptation disorder rather than a habit.
3. Jarvis MJ. **Why people smoke.** BMJ. 2004;328(7434):277–279. <https://pubmed.ncbi.nlm.nih.gov/14751901/>
4. Heatherton TF, Kozlowski LT, Frecker RC, Fagerström KO. **The Fagerström Test for Nicotine Dependence: a revision of the Fagerström Tolerance Questionnaire.** Br J Addict. 1991;86(9):1119–1127. <https://pubmed.ncbi.nlm.nih.gov/1932883/> — → model: `FTND`/`HSI` node; time-to-first-cigarette is the emergent consequence of overnight resensitization in the model.
5. Hughes JR. **Effects of abstinence from tobacco: valid symptoms and time course.** Nicotine Tob Res. 2007;9(3):315–327. <https://pubmed.ncbi.nlm.nih.gov/17365764/> — → model: **primary calibration target** for the withdrawal state `WD` — peak at day 1–3, resolution over 2–4 weeks.
6. Shiffman S, Patten C, Gwaltney C, et al. **Natural history of nicotine withdrawal.** Addiction. 2006;101(12):1822–1832. <https://pubmed.ncbi.nlm.nih.gov/17156182/> — → model: shape of the `WD` trajectory and its within-day variability.
7. Piper ME. **Withdrawal: expanding a key addiction construct.** Nicotine Tob Res. 2015;17(12):1405–1415. <https://pubmed.ncbi.nlm.nih.gov/25744958/>
8. Cox LS, Tiffany ST, Christen AG. **Evaluation of the brief questionnaire of smoking urges (QSU-brief) in laboratory and clinical settings.** Nicotine Tob Res. 2001;3(1):7–16. <https://pubmed.ncbi.nlm.nih.gov/11260806/> — → model: `QSU` compartment scale (1–7) and its two-factor (tonic + cue) structure.
9. Hughes JR, Keely J, Naud S. **Shape of the relapse curve and long-term abstinence among untreated smokers.** Addiction. 2004;99(1):29–38. <https://pubmed.ncbi.nlm.nih.gov/14678060/> — → model: justifies the exponential-survival form of `PABST` and the front-loaded hazard.
10. Chaiton M, Diemert L, Cohen JE, et al. **Estimating the number of quit attempts it takes to quit smoking successfully.** BMJ Open. 2016;6(6):e011045. <https://pubmed.ncbi.nlm.nih.gov/27288378/>

## 2. Nicotine pharmacokinetics & delivery kinetics

11. Benowitz NL, Hukkanen J, Jacob P 3rd. **Nicotine chemistry, metabolism, kinetics and biomarkers.** Handb Exp Pharmacol. 2009;192:29–60. <https://pubmed.ncbi.nlm.nih.gov/19184645/> — → model: `CLNP`/`CLN2A6` (total CL ≈ 72 L/h), `V2N`/`V3N` (Vss ≈ 182 L), nicotine t½ ≈ 1.8 h, cotinine t½ 16 h.
12. Hukkanen J, Jacob P 3rd, Benowitz NL. **Metabolism and disposition kinetics of nicotine.** Pharmacol Rev. 2005;57(1):79–115. <https://pubmed.ncbi.nlm.nih.gov/15734728/> — → model: metabolic fractions `FCOT` = 0.70 (nicotine→cotinine) and `FHC` = 0.45 (cotinine→3HC); minor UGT2B10/FMO3 routes.
13. Benowitz NL, Jacob P 3rd. **Intravenous nicotine replacement suppresses nicotine intake from cigarette smoking.** J Pharmacol Exp Ther. 1990;254(3):1000–1005. <https://pubmed.ncbi.nlm.nih.gov/2203896/> — → model: empirical basis for the titration/substitution logic linking plasma nicotine to `CPD`.
14. Shiffman S, Dresler CM, Rohay JM. **Nicotine delivery systems.** Expert Opin Drug Deliv. 2005;2(3):563–577. <https://pubmed.ncbi.nlm.nih.gov/16296775/> — → model: route-specific absorption constants `KASKIN` (patch, t½ ≈ 12 h) and `KAMOU` (buccal), and the phasic-fraction parameters `FPHPATCH` / `FPHORAL` / `FPHCIG`.
15. Benowitz NL, Hall SM, Stewart S, et al. **Nicotine intake and dose response when smoking reduced-nicotine content cigarettes.** Clin Pharmacol Ther. 2006;80(6):703–714. <https://pubmed.ncbi.nlm.nih.gov/17178270/> — → model: `NICCIG` (mg absorbed per cigarette) and the VLNC experiment.
16. Middleton ET, Morice AH. **Breath carbon monoxide as an indication of smoking habit.** Chest. 2000;117(3):758–763. <https://pubmed.ncbi.nlm.nih.gov/10713003/> — → model: `COHB` kinetics (t½ 4 h) and the exhaled-CO ≈ 4 × COHb% conversion used for CO-verified abstinence.

## 3. CYP2A6 pharmacogenetics & the nicotine metabolite ratio (NMR)

17. Dempsey D, Tutka P, Jacob P 3rd, et al. **Nicotine metabolite ratio as an index of cytochrome P450 2A6 metabolic activity.** Clin Pharmacol Ther. 2004;76(1):64–72. <https://pubmed.ncbi.nlm.nih.gov/15229465/> — → model: **defines the `NMR` readout** (3HC/cotinine) and its dependence on a single CYP2A6 activity scalar `F2A6`, which is exactly how the model implements it.
18. Malaiyandi V, Sellers EM, Tyndale RF. **Implications of CYP2A6 genetic variation for smoking behaviors and nicotine dependence.** Clin Pharmacol Ther. 2005;77(3):145–158. <https://pubmed.ncbi.nlm.nih.gov/15735609/> — → model: `F2A6` range (0.35 slow ↔ 1.7–2.0 fast) mapped to *4/*9/*12/*17 alleles.
19. Benowitz NL, Pomerleau OF, Pomerleau CS, Jacob P 3rd. **Nicotine metabolite ratio as a predictor of cigarette consumption.** Nicotine Tob Res. 2003;5(5):621–624. <https://pubmed.ncbi.nlm.nih.gov/14577978/> — → model: the `NMR_FAST → CPD` edge.
20. Bloom AJ, Baker TB, Chen LS, et al. **The contribution of common UGT2B10 and CYP2A6 alleles to variation in nicotine metabolism among European Americans.** Pharmacogenet Genomics. 2013;23(12):706–716. <https://pubmed.ncbi.nlm.nih.gov/24192532/>
21. Lerman C, Schnoll RA, Hawk LW Jr, et al. **Use of the nicotine metabolite ratio as a genetically informed biomarker of response to nicotine patch or varenicline for smoking cessation: a randomised, double-blind placebo-controlled trial.** Lancet Respir Med. 2015;3(2):131–138. <https://pubmed.ncbi.nlm.nih.gov/25588294/> — → model: **the target of experiment 5b.** The model reproduces this NMR × treatment interaction from one PK parameter, with no fitted interaction term.
22. Schnoll RA, Patterson F, Wileyto EP, et al. **Nicotine metabolic rate predicts successful smoking cessation with transdermal nicotine.** Pharmacol Biochem Behav. 2009;92(1):6–11. <https://pubmed.ncbi.nlm.nih.gov/19000709/>
23. Cosgrove KP, Esterlis I, McKee SA, et al. **Sex differences in availability of β2*-nicotinic acetylcholine receptors in recently abstinent tobacco smokers.** Arch Gen Psychiatry. 2012;69(4):418–427. <https://pubmed.ncbi.nlm.nih.gov/22474108/> — → model: the `COV_SEX → CYP2A6` covariate edge and sex differences in receptor availability.

## 4. α4β2* nicotinic receptor pharmacology, occupancy & upregulation

24. Brody AL, Mandelkern MA, London ED, et al. **Cigarette smoking saturates brain α4β2 nicotinic acetylcholine receptors.** Arch Gen Psychiatry. 2006;63(8):907–915. <https://pubmed.ncbi.nlm.nih.gov/16894067/> — → model: **the single most important calibration.** `KDNIC` = 5.4 nM is set so that 0.87 ng/mL plasma nicotine gives 50% β2* occupancy, reproducing this paper's reported EC50; one cigarette then saturates the receptor.
25. Cosgrove KP, Batis J, Bois F, et al. **β2-Nicotinic acetylcholine receptor availability during acute and prolonged abstinence from tobacco smoking.** Arch Gen Psychiatry. 2009;66(6):666–676. <https://pubmed.ncbi.nlm.nih.gov/19487632/> — → model: `KOUTR` (t½ = 14 d) — the receptor-pool time constant that makes the 3–4 week normalization emerge rather than be imposed.
26. Marks MJ, Burch JB, Collins AC. **Effects of chronic nicotine infusion on tolerance development and nicotinic receptors.** J Pharmacol Exp Ther. 1983;226(3):817–825. <https://pubmed.ncbi.nlm.nih.gov/6887012/> — → model: `EMAXR` — desensitization-driven upregulation of the receptor pool.
27. Fenster CP, Rains MF, Noerager B, Quick MW, Lester RA. **Influence of subunit composition on desensitization of neuronal acetylcholine receptors at low concentrations of nicotine.** J Neurosci. 1997;17(15):5747–5759. <https://pubmed.ncbi.nlm.nih.gov/9221773/> — → model: `KOND` (≈7 min on) and `KOFFD` (≈2 h resensitization).
28. Nashmi R, Xiao C, Deshpande P, et al. **Chronic nicotine cell specifically upregulates functional α4* nicotinic receptors: basis for both tolerance in midbrain and enhanced long-term potentiation in perforant path.** J Neurosci. 2007;27(31):8202–8218. <https://pubmed.ncbi.nlm.nih.gov/17670967/> — → model: cell-type-specific upregulation motivating `PHID` (partial, not complete, signal loss per unit desensitization).
29. Mansvelder HD, McGehee DS. **Long-term potentiation of excitatory inputs to brain reward areas by nicotine.** Neuron. 2000;27(2):349–357. <https://pubmed.ncbi.nlm.nih.gov/10985354/> — → model: α7-mediated glutamatergic drive onto VTA DA neurons.
30. Mineur YS, Picciotto MR. **Nicotine receptors and depression: revisiting and revising the cholinergic hypothesis.** Trends Pharmacol Sci. 2010;31(12):580–586. <https://pubmed.ncbi.nlm.nih.gov/20965579/> — → model: `DEPR_COMORB` and the `PSYHX` covariate on lapse hazard.

## 5. Mesolimbic dopamine reward circuitry

31. Di Chiara G, Imperato A. **Drugs abused by humans preferentially increase synaptic dopamine concentrations in the mesolimbic system of freely moving rats.** Proc Natl Acad Sci U S A. 1988;85(14):5274–5278. <https://pubmed.ncbi.nlm.nih.gov/2899326/> — → model: `EMAXDA` — the magnitude of the acute NAc dopamine response to nicotine.
32. Maskos U, Molles BE, Pons S, et al. **Nicotine reinforcement and cognition restored by targeted expression of nicotinic receptors.** Nature. 2005;436(7047):103–107. <https://pubmed.ncbi.nlm.nih.gov/16001069/> — → model: β2-containing receptors in VTA as the necessary node linking `OCC` → `VTA_DA`.
33. Tapper AR, McKinney SL, Nashmi R, et al. **Nicotine activation of α4* receptors: sufficient for reward, tolerance, and sensitization.** Science. 2004;306(5698):1029–1032. <https://pubmed.ncbi.nlm.nih.gov/15528443/> — → model: justifies α4β2* as the single reward-relevant receptor in the activation term.
34. Pons S, Fattore L, Cossu G, et al. **Crucial role of α4 and α6 nicotinic acetylcholine receptor subunits from ventral tegmental area in systemic nicotine self-administration.** J Neurosci. 2008;28(47):12318–12327. <https://pubmed.ncbi.nlm.nih.gov/19020025/> — → model: `A6B2B3 → NAC_DA` presynaptic release edge.
35. Fowler JS, Logan J, Wang GJ, Volkow ND. **Monoamine oxidase and cigarette smoking.** Neurotoxicology. 2003;24(1):75–82. <https://pubmed.ncbi.nlm.nih.gov/12564384/> — → model: `MAOAB` node — a non-nicotine smoke constituent that prolongs dopamine signalling and is absent from all pharmacotherapies.

## 6. Habenula–IPN aversion circuit, α5 and intake titration

36. Fowler CD, Lu Q, Johnson PM, Marks MJ, Kenny PJ. **Habenular α5 nicotinic receptor subunit signalling controls nicotine intake.** Nature. 2011;471(7340):597–601. <https://pubmed.ncbi.nlm.nih.gov/21278726/> — → model: the `MHb → IPN → AVERS_SP → DOSE_CEIL` chain and the inverted-U self-administration ceiling.
37. Frahm S, Ślimak MA, Ferrarese L, et al. **Aversion to nicotine is regulated by the balanced activity of β4 and α5 nicotinic receptor subunits in the medial habenula.** Neuron. 2011;70(3):522–535. <https://pubmed.ncbi.nlm.nih.gov/21555077/>
38. Thorgeirsson TE, Geller F, Sulem P, et al. **A variant associated with nicotine dependence, lung cancer and peripheral arterial disease.** Nature. 2008;452(7187):638–642. <https://pubmed.ncbi.nlm.nih.gov/18385739/> — → model: `A5RISK` covariate (rs16969968) on the lapse hazard.
39. Olfson E, Saccone NL, Johnson EO, et al. **Rare, low frequency and common coding variants in CHRNA5 and their contribution to nicotine dependence in European and African Americans.** Mol Psychiatry. 2016;21(5):601–607. <https://pubmed.ncbi.nlm.nih.gov/26239294/>

## 7. Allostasis, negative affect & the opponent-process account of withdrawal

40. Koob GF, Le Moal M. **Drug addiction, dysregulation of reward, and allostasis.** Neuropsychopharmacology. 2001;24(2):97–129. <https://pubmed.ncbi.nlm.nih.gov/11120394/> — → model: **the structural basis of the `SETP` compartment.** The hedonic set-point as a slow moving average of dopaminergic tone, with withdrawal as the transient mismatch `DEF = SETP − DA`.
41. George O, Ghozland S, Azar MR, et al. **CRF–CRF1 system activation mediates withdrawal-induced increases in nicotine self-administration in nicotine-dependent rats.** Proc Natl Acad Sci U S A. 2007;104(43):17198–17203. <https://pubmed.ncbi.nlm.nih.gov/17921249/> — → model: `ALLO` compartment (CRF/extended-amygdala recruitment) and its slower time constant than `DEF`.
42. Bruijnzeel AW. **κ-Opioid receptor signaling and brain reward function.** Brain Res Rev. 2009;62(1):127–146. <https://pubmed.ncbi.nlm.nih.gov/19804796/> — → model: `DYN_KOR` contribution to the dysphoria component of `WD`.

## 8. Varenicline & cytisinicline — partial agonist pharmacology and trials

43. Rollema H, Chambers LK, Coe JW, et al. **Pharmacological profile of the α4β2 nicotinic acetylcholine receptor partial agonist varenicline, an effective smoking cessation aid.** Neuropharmacology. 2007;52(3):985–994. <https://pubmed.ncbi.nlm.nih.gov/17157884/> — → model: `EMAXV` = 0.45 (intrinsic activity ≈ 45% of nicotine) and the dual relief-plus-blockade mechanism formalised in experiment 5c.
44. Faessel HM, Gibbs MA, Clark DJ, et al. **Multiple-dose pharmacokinetics of the selective nicotinic receptor partial agonist, varenicline, in healthy smokers.** J Clin Pharmacol. 2006;46(12):1439–1448. <https://pubmed.ncbi.nlm.nih.gov/17101743/> — → model: `CLV` = 12 L/h, `VV` = 415 L, t½ 24 h, and the renal-impairment scalar `RFV`.
45. Lotfipour S, Mandelkern M, Alvarez-Estrada M, Brody AL. **A single administration of low-dose varenicline saturates α4β2* nicotinic acetylcholine receptors in the human brain.** Neuropsychopharmacology. 2012;37(7):1738–1748. <https://pubmed.ncbi.nlm.nih.gov/22395733/> — → model: `KDVAR` = 3.6 nM, calibrated so that 1 mg BID yields ≈ 90% β2* occupancy.
46. Gonzales D, Rennard SI, Nides M, et al. **Varenicline, an α4β2 nicotinic acetylcholine receptor partial agonist, vs sustained-release bupropion and placebo for smoking cessation: a randomized controlled trial.** JAMA. 2006;296(1):47–55. <https://pubmed.ncbi.nlm.nih.gov/16820546/>
47. Tonstad S, Tønnesen P, Hajek P, et al. **Effect of maintenance therapy with varenicline on smoking cessation: a randomized controlled trial.** JAMA. 2006;296(1):64–71. <https://pubmed.ncbi.nlm.nih.gov/16820548/> — → model: rationale for the `TSTOPRX` treatment-duration parameter and post-treatment hazard rebound.
48. West R, Baker CL, Cappelleri JC, Bushmakin AG. **Effect of varenicline and bupropion SR on craving, nicotine withdrawal symptoms, and rewarding effects of smoking during a quit attempt.** Psychopharmacology (Berl). 2008;197(3):371–377. <https://pubmed.ncbi.nlm.nih.gov/18084743/> — → model: the empirical separation of *withdrawal relief* from *reward blockade* that `GREL` and `GBLOCK` encode as distinct terms.
49. Anthenelli RM, Benowitz NL, West R, et al. **Neuropsychiatric safety and efficacy of varenicline, bupropion, and nicotine patch in smokers with and without psychiatric disorders (EAGLES): a double-blind, randomised, placebo-controlled clinical trial.** Lancet. 2016;387(10037):2507–2520. <https://pubmed.ncbi.nlm.nih.gov/27116918/> — → model: **primary efficacy calibration.** CAR weeks 9–12 = varenicline 33.5% / bupropion 22.6% / patch 23.4% / placebo 12.5%; also the `NPS_SAFETY` node.
50. Benowitz NL, Pipe A, West R, et al. **Cardiovascular safety of varenicline, bupropion, and nicotine patch in smokers: a randomized clinical trial.** JAMA Intern Med. 2018;178(5):622–631. <https://pubmed.ncbi.nlm.nih.gov/29630702/> — → model: `CV_SAFETY` node (CATS).
51. Koegelenberg CFN, Noor F, Bateman ED, et al. **Efficacy of varenicline combined with nicotine replacement therapy vs varenicline alone for smoking cessation: a randomized clinical trial.** JAMA. 2014;312(2):155–161. <https://pubmed.ncbi.nlm.nih.gov/25005652/> — → model: scenario 9 (varenicline + patch) target.
52. Ebbert JO, Hughes JR, West RJ, et al. **Effect of varenicline on smoking cessation through smoking reduction: a randomized clinical trial.** JAMA. 2015;313(7):687–694. <https://pubmed.ncbi.nlm.nih.gov/25688780/> — → model: scenario 6 (preloading / flexible quit date).
53. Johnstone S, Hughes JR, et al. **Evaluating mediators of the effect of varenicline preloading on smoking cessation.** Addiction. 2025. <https://pubmed.ncbi.nlm.nih.gov/39915904/>
54. Rigotti NA, Benowitz NL, Prochaska JJ, et al. **Cytisinicline for smoking cessation: a randomized clinical trial (ORCA-2).** JAMA. 2023;330(2):152–160. <https://pubmed.ncbi.nlm.nih.gov/37432430/> — → model: scenario 8 target; `KDCYT` = 60 nM reflects cytisine's limited brain penetration.
55. Walker N, Howe C, Glover M, et al. **Cytisine versus nicotine for smoking cessation.** N Engl J Med. 2014;371(25):2353–2362. <https://pubmed.ncbi.nlm.nih.gov/25517706/>
56. Cahill K, Lindson-Hawley N, Thomas KH, Fanshawe TR, Lancaster T. **Nicotine receptor partial agonists for smoking cessation.** Cochrane Database Syst Rev. 2016;(5):CD006103. <https://pubmed.ncbi.nlm.nih.gov/27158893/>

## 9. Bupropion, NRT & other pharmacotherapy

57. Slemmer JE, Martin BR, Damaj MI. **Bupropion is a nicotinic antagonist.** J Pharmacol Exp Ther. 2000;295(1):321–327. <https://pubmed.ncbi.nlm.nih.gov/10991997/> — → model: `KINACHR` — non-competitive nAChR block by hydroxybupropion, the `BLK` term.
58. Hsyu PH, Singh A, Giargiari TD, et al. **Pharmacokinetics of bupropion and its metabolites in cigarette smokers versus nonsmokers.** J Clin Pharmacol. 1997;37(8):737–743. <https://pubmed.ncbi.nlm.nih.gov/9378846/> — → model: `CLB`, `FOH`, `CLOH` — parent and hydroxybupropion exposures (metabolite ≈ 6–8 × parent).
59. Hurt RD, Sachs DP, Glover ED, et al. **A comparison of sustained-release bupropion and placebo for smoking cessation.** N Engl J Med. 1997;337(17):1195–1202. <https://pubmed.ncbi.nlm.nih.gov/9337378/>
60. Hughes JR, Stead LF, Lancaster T. **Antidepressants for smoking cessation.** Cochrane Database Syst Rev. 2003. <https://pubmed.ncbi.nlm.nih.gov/12804385/>
61. Hartmann-Boyce J, Chepkin SC, Ye W, Bullen C, Lancaster T. **Nicotine replacement therapy versus control for smoking cessation.** Cochrane Database Syst Rev. 2018;(5):CD000146. <https://pubmed.ncbi.nlm.nih.gov/29852054/> — → model: patch-arm efficacy target.
62. Theodoulou A, Chepkin SC, Ye W, et al. **Different doses, durations and modes of delivery of nicotine replacement therapy for smoking cessation.** Cochrane Database Syst Rev. 2023;(6):CD013308. <https://pubmed.ncbi.nlm.nih.gov/37335995/> — → model: combination-NRT scenario (patch + PRN) and the phasic/tonic split that explains its advantage.
63. Donny EC, Denlinger RL, Tidey JW, et al. **Randomized trial of reduced-nicotine standards for cigarettes.** N Engl J Med. 2015;373(14):1340–1349. <https://pubmed.ncbi.nlm.nih.gov/26422724/> — → model: VLNC experiment 5e; the ~0.03 mg/cigarette exposure at which dependence variables finally move.
64. Pisinger C, Godtfredsen NS. **Is there a health benefit of reduced tobacco consumption? A systematic review.** Nicotine Tob Res. 2007;9(6):631–646. <https://pubmed.ncbi.nlm.nih.gov/17558820/>

## 10. Organ-system consequences, safety & drug interactions

65. Fletcher C, Peto R. **The natural history of chronic airflow obstruction.** Br Med J. 1977;1(6077):1645–1648. <https://pubmed.ncbi.nlm.nih.gov/871704/> — → model: `FEVSMOKE` = 60 vs `FEVQUIT` = 30 mL/yr.
66. Anthonisen NR, Connett JE, Kiley JP, et al. **Effects of smoking intervention and the use of an inhaled anticholinergic bronchodilator on the rate of decline of FEV1 (Lung Health Study).** JAMA. 1994;272(19):1497–1505. <https://pubmed.ncbi.nlm.nih.gov/7966841/> — → model: post-cessation FEV1 slope recovery.
67. Aubin HJ, Farley A, Lycett D, Lahmek P, Aveyard P. **Weight gain in smokers after quitting cigarettes: meta-analysis.** BMJ. 2012;345:e4439. <https://pubmed.ncbi.nlm.nih.gov/22782848/> — → model: `WTSS` = 5 kg with `KWT` t½ ≈ 90 d.
68. Mineur YS, Abizaid A, Rao Y, et al. **Nicotine decreases food intake through activation of POMC neurons.** Science. 2011;332(6035):1330–1332. <https://pubmed.ncbi.nlm.nih.gov/21659607/> — → model: the `APPETITE → WEIGHT` mechanism.
69. Benowitz NL, Burbank AD. **Cardiovascular toxicity of nicotine: implications for electronic cigarette use.** Trends Cardiovasc Med. 2016;26(6):515–523. <https://pubmed.ncbi.nlm.nih.gov/27079891/> — → model: `HR_BP` node; separates nicotine effects from combustion-product effects.
70. Zevin S, Benowitz NL. **Drug interactions with tobacco smoking. An update.** Clin Pharmacokinet. 1999;36(6):425–438. <https://pubmed.ncbi.nlm.nih.gov/10427467/> — → model: `CYP1A2_IND → DDI` — induction by PAHs, **not** nicotine; the clozapine/olanzapine/theophylline dose trap on quitting.
71. Critchley JA, Capewell S. **Mortality risk reduction associated with smoking cessation in patients with coronary heart disease: a systematic review.** JAMA. 2003;290(1):86–97. <https://pubmed.ncbi.nlm.nih.gov/12837716/>

## 11. QSP / pharmacometric methodology

72. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883–893. <https://pubmed.ncbi.nlm.nih.gov/31652028/> — → model: the mrgsolve implementation conventions used throughout `tud_mrgsolve_model.R`.

### Tooling

- mrgsolve documentation — <https://mrgsolve.org/>
- Graphviz (`dot`) — <https://graphviz.org/>
- Shiny — <https://shiny.posit.co/>

---

## Reference-to-parameter map (quick index)

| Model quantity | Value used | Anchor |
|---|---|---|
| Nicotine CL / Vss / t½ | 72 L/h · 182 L · 1.8 h | ref 11 |
| Cotinine t½ / 3HC route | 16 h · `FHC` 0.45 | refs 11, 12 |
| NMR (3HC/cotinine) | 0.12 (slow) – 0.59 (fast) | refs 17, 18, 21 |
| `KDNIC` (β2* effective Kd) | 5.4 nM ⇒ EC50 0.87 ng/mL | ref 24 |
| `KDVAR` | 3.6 nM ⇒ ~90% occupancy at 1 mg BID | ref 45 |
| `KOND` / `KOFFD` | 7 min on / 2 h off | ref 27 |
| `KOUTR` (receptor pool) | t½ 14 d ⇒ 3–4 wk normalization | refs 25, 26 |
| `EMAXV` / `EMAXC` | 0.45 / 0.30 intrinsic activity | refs 43, 54 |
| `SETP` (hedonic set-point) | t½ 14 d moving average of DA | ref 40 |
| `WD` peak timing | day 1–3, resolves 2–4 wk | refs 5, 6 |
| CAR weeks 9–12 targets | VAR 33.5 / BUP 22.6 / patch 23.4 / PBO 12.5 % | ref 49 |
| Cytisinicline CAR | ~32% vs ~7% placebo | ref 54 |
| Varenicline + patch | superior to varenicline alone | ref 51 |
| Post-cessation weight | ~5 kg | ref 67 |
| FEV1 slope | 60 → 30 mL/yr | refs 65, 66 |
| Exhaled CO conversion | ppm ≈ 4 × COHb% | ref 16 |

---

*This reference list supports an educational QSP model. It is not a clinical
guideline. See the disclaimer in the repository README.*
