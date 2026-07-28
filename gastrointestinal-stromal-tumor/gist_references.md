# Gastrointestinal Stromal Tumor (GIST) — References

Curated bibliography supporting the mechanistic map, mrgsolve QSP model, and
Shiny dashboard in this directory. Organized by section.

**Every PMID below was verified against the NCBI PubMed E-utilities API at the
time this file was generated** — each link resolves to the paper named beside
it. Where a claim in the model rests on a source whose PMID could not be
confirmed, no link is given rather than a guessed one.

---

## 1. Cell of Origin, Lineage Dependency & Molecular Pathogenesis

1. Hirota S, Isozaki K, Moriyama Y, et al. Gain-of-function mutations of c-kit in human gastrointestinal stromal tumors. *Science*. 1998;279(5350):577-580. [PMID 9438854](https://pubmed.ncbi.nlm.nih.gov/9438854/) — the founding observation that GIST is a KIT-driven disease of interstitial-cell-of-Cajal lineage.
2. Heinrich MC, Corless CL, Duensing A, et al. PDGFRA activating mutations in gastrointestinal stromal tumors. *Science*. 2003;299(5607):708-710. [PMID 12522257](https://pubmed.ncbi.nlm.nih.gov/12522257/) — establishes KIT and PDGFRA as mutually exclusive alternative drivers.
3. Corless CL, Barnett CM, Heinrich MC. Gastrointestinal stromal tumours: origin and molecular oncology. *Nat Rev Cancer*. 2011;11(12):865-878. [PMID 22089421](https://pubmed.ncbi.nlm.nih.gov/22089421/) — the standard synthesis of GIST molecular subtypes; the backbone of cluster 2 of the DOT map.
4. Chi P, Chen Y, Zhang L, et al. ETV1 is a lineage survival factor that cooperates with KIT in gastrointestinal stromal tumours. *Nature*. 2010;467(7317):849-853. [PMID 20927104](https://pubmed.ncbi.nlm.nih.gov/20927104/) — the ERK→ETV1-stabilization loop modeled in the `ETV1` compartment.
5. Kawanowa K, Sakuma Y, Sakurai S, et al. High incidence of microscopic gastrointestinal stromal tumors in the stomach. *Hum Pathol*. 2006;37(12):1527-1535. [PMID 16996566](https://pubmed.ncbi.nlm.nih.gov/16996566/) — microscopic "tumorlets" are common; driver mutation alone is not sufficient for clinical disease.

## 2. Risk Stratification & Natural History

6. Miettinen M, Lasota J. Gastrointestinal stromal tumors: pathology and prognosis at different sites. *Semin Diagn Pathol*. 2006;23(2):70-83. [PMID 17193820](https://pubmed.ncbi.nlm.nih.gov/17193820/) — the AFIP size × mitotic index × site risk scheme (`Risk_strat` node).
7. Casali PG, Blay JY, Abecassis N, et al. Gastrointestinal stromal tumours: ESMO-EURACAN-GENTURIS Clinical Practice Guidelines for diagnosis, treatment and follow-up. *Ann Oncol*. 2022;33(1):20-33. [PMID 34560242](https://pubmed.ncbi.nlm.nih.gov/34560242/) — current line-sequencing and genotyping standard encoded in the management cluster.

## 3. Genotype–Response Relationships (first line)

8. Demetri GD, von Mehren M, Blanke CD, et al. Efficacy and safety of imatinib mesylate in advanced gastrointestinal stromal tumors. *N Engl J Med*. 2002;347(7):472-480. [PMID 12181401](https://pubmed.ncbi.nlm.nih.gov/12181401/) — B2222; the response magnitude and time-to-progression the model's first-line scenario is calibrated against.
9. Heinrich MC, Corless CL, Demetri GD, et al. Kinase mutations and imatinib response in patients with metastatic gastrointestinal stromal tumor. *J Clin Oncol*. 2003;21(23):4342-4349. [PMID 14645423](https://pubmed.ncbi.nlm.nih.gov/14645423/) — exon 11 > exon 9 > wild-type response gradient; the basis of `GF_EX9`.
10. Verweij J, Casali PG, Zalcberg J, et al. Progression-free survival in gastrointestinal stromal tumours with high-dose imatinib: randomised trial. *Lancet*. 2004;364(9440):1127-1134. [PMID 15451219](https://pubmed.ncbi.nlm.nih.gov/15451219/) — EORTC 62005, 400 vs 800 mg.
11. Blanke CD, Rankin C, Demetri GD, et al. Phase III randomized, intergroup trial assessing imatinib mesylate at two dose levels in patients with unresectable or metastatic GIST expressing KIT. *J Clin Oncol*. 2008;26(4):626-632. [PMID 18235122](https://pubmed.ncbi.nlm.nih.gov/18235122/) — S0033.
12. Gastrointestinal Stromal Tumor Meta-Analysis Group (MetaGIST). Comparison of two doses of imatinib for the treatment of unresectable or metastatic gastrointestinal stromal tumors: a meta-analysis of 1,640 patients. *J Clin Oncol*. 2010;28(7):1247-1253. [PMID 20124181](https://pubmed.ncbi.nlm.nih.gov/20124181/) — the exon 9 / 800 mg PFS benefit that scenario 4 targets.

## 4. Secondary Resistance & Clonal Evolution

13. Antonescu CR, Besmer P, Guo T, et al. Acquired resistance to imatinib in gastrointestinal stromal tumor occurs through secondary gene mutation. *Clin Cancer Res*. 2005;11(11):4182-4190. [PMID 15930355](https://pubmed.ncbi.nlm.nih.gov/15930355/) — resistance is target-mutation, not drug-delivery, failure.
14. Wardelmann E, Merkelbach-Bruse S, Pauls K, et al. Polyclonal evolution of multiple secondary KIT mutations in gastrointestinal stromal tumors under treatment with imatinib mesylate. *Clin Cancer Res*. 2006;12(6):1743-1749. [PMID 16551858](https://pubmed.ncbi.nlm.nih.gov/16551858/) — **the single most important paper for this model's architecture**: different lesions in one patient carry different secondary mutations, which is why a one-compartment tumor model cannot work.
15. Heinrich MC, Maki RG, Corless CL, et al. Primary and secondary kinase genotypes correlate with the biological and clinical activity of sunitinib in imatinib-resistant gastrointestinal stromal tumor. *J Clin Oncol*. 2008;26(33):5352-5359. [PMID 18955451](https://pubmed.ncbi.nlm.nih.gov/18955451/) — sunitinib is active against ATP-binding-pocket (exon 13/14) but not activation-loop (exon 17/18) mutants; the asymmetry that generates the model's `IC_SUN_A` vs `IC_SUN_L` split.
16. Serrano C, Mariño-Enríquez A, Tao DL, et al. Complementary activity of tyrosine kinase inhibitors against secondary kit mutations in imatinib-resistant gastrointestinal stromal tumours. *Br J Cancer*. 2019;120(6):612-620. [PMID 30792533](https://pubmed.ncbi.nlm.nih.gov/30792533/) — the explicit drug × secondary-mutation potency matrix this model implements.
17. Guo T, Hajdu M, Agaram NP, et al. Mechanisms of sunitinib resistance in gastrointestinal stromal tumors harboring KIT AY502-503 insertion. *Clin Cancer Res*. 2009;15(22):6862-6870. [PMID 19861442](https://pubmed.ncbi.nlm.nih.gov/19861442/).
18. Agaram NP, Wong GC, Guo T, et al. Novel V600E BRAF mutations in imatinib-naive and imatinib-resistant gastrointestinal stromal tumors. *Genes Chromosomes Cancer*. 2008;47(10):853-859. [PMID 18615679](https://pubmed.ncbi.nlm.nih.gov/18615679/) — KIT-independent bypass drivers.
19. Gupta A, Roy S, Lazar AJF, et al. Autophagy inhibition and antimalarial drugs enhance cell death in gastrointestinal stromal tumor. *Proc Natl Acad Sci U S A*. 2010;107(32):14333-14338. [PMID 20660757](https://pubmed.ncbi.nlm.nih.gov/20660757/) — the autophagy-dependent quiescent survivor population represented by the `TQ` compartment.
20. Ng KP, Hillmer AM, Chuah CTH, et al. A common BIM deletion polymorphism mediates intrinsic resistance and inferior responses to tyrosine kinase inhibitors in cancer. *Nat Med*. 2012;18(4):521-528. [PMID 22426421](https://pubmed.ncbi.nlm.nih.gov/22426421/) — the intron-2 2,903-bp deletion (≈12% of East Asians) behind `BIMFAC_DEL`.

## 5. Second-, Third- and Fourth-Line Pharmacology

21. Demetri GD, van Oosterom AT, Garrett CR, et al. Efficacy and safety of sunitinib in patients with advanced gastrointestinal stromal tumour after failure of imatinib: a randomised controlled trial. *Lancet*. 2006;368(9544):1329-1338. [PMID 17046465](https://pubmed.ncbi.nlm.nih.gov/17046465/).
22. Demetri GD, Reichardt P, Kang YK, et al. Efficacy and safety of regorafenib for advanced gastrointestinal stromal tumours after failure of imatinib and sunitinib (GRID): an international, multicentre, randomised, placebo-controlled, phase 3 trial. *Lancet*. 2013;381(9863):295-302. [PMID 23177515](https://pubmed.ncbi.nlm.nih.gov/23177515/).
23. George S, Wang Q, Heinrich MC, et al. Efficacy and safety of regorafenib in patients with metastatic and/or unresectable GI stromal tumor after failure of imatinib and sunitinib: a multicenter phase II trial. *J Clin Oncol*. 2012;30(19):2401-2407. [PMID 22614970](https://pubmed.ncbi.nlm.nih.gov/22614970/).
24. Smith BD, Kaufman MD, Lu WP, et al. Ripretinib (DCC-2618) is a switch control kinase inhibitor of a broad spectrum of oncogenic and drug-resistant KIT and PDGFRA variants. *Cancer Cell*. 2019;35(5):738-751.e9. [PMID 31085175](https://pubmed.ncbi.nlm.nih.gov/31085175/) — the dual switch-pocket/activation-loop binding mode behind ripretinib's broad `IC_RIP_*` row.
25. Blay JY, Serrano C, Heinrich MC, et al. Ripretinib in patients with advanced gastrointestinal stromal tumours (INVICTUS): a double-blind, randomised, placebo-controlled, phase 3 trial. *Lancet Oncol*. 2020;21(7):923-934. [PMID 32511981](https://pubmed.ncbi.nlm.nih.gov/32511981/) — fourth-line registration trial.
26. Bauer S, Jones RL, Blay JY, et al. Ripretinib versus sunitinib in patients with advanced gastrointestinal stromal tumor after treatment with imatinib (INTRIGUE): a randomized, open-label, phase III trial. *J Clin Oncol*. 2022;40(34):3918-3928. [PMID 35947817](https://pubmed.ncbi.nlm.nih.gov/35947817/) — **negative** for the primary endpoint in the second line.
27. Bauer S, Heinrich MC, George S, et al. Ripretinib versus sunitinib in gastrointestinal stromal tumor: ctDNA biomarker analysis of the phase 3 INTRIGUE trial. *Nat Med*. 2024;30(2):498-506. [PMID 38182785](https://pubmed.ncbi.nlm.nih.gov/38182785/) — the ctDNA-defined subgroups (KIT exon 11+13/14 favouring sunitinib, exon 11+17/18 favouring ripretinib) that scenarios 9a–9d reproduce.

## 6. Genotype-Directed Agents (PDGFRA D842V, SDH-deficient, rare drivers)

28. Heinrich MC, Jones RL, von Mehren M, et al. Avapritinib in advanced PDGFRA D842V-mutant gastrointestinal stromal tumour (NAVIGATOR): a multicentre, open-label, phase 1 trial. *Lancet Oncol*. 2020;21(7):935-946. [PMID 32615108](https://pubmed.ncbi.nlm.nih.gov/32615108/) — the type I inhibitor that converts an untreatable genotype into a highly responsive one.
29. Kang YK, George S, Jones RL, et al. Avapritinib versus regorafenib in locally advanced unresectable or metastatic GI stromal tumor: a randomized, open-label phase III study (VOYAGER). *J Clin Oncol*. 2021;39(28):3128-3139. [PMID 34343033](https://pubmed.ncbi.nlm.nih.gov/34343033/) — avapritinib's potency does **not** generalize beyond D842V, consistent with the model's `IC_AVA_A` penalty on V654A.
30. Janeway KA, Kim SY, Lodish M, et al. Defects in succinate dehydrogenase in gastrointestinal stromal tumors lacking KIT and PDGFRA mutations. *Proc Natl Acad Sci U S A*. 2011;108(1):314-318. [PMID 21173220](https://pubmed.ncbi.nlm.nih.gov/21173220/).
31. Killian JK, Kim SY, Miettinen M, et al. Succinate dehydrogenase mutation underlies global epigenomic divergence in gastrointestinal stromal tumor. *Cancer Discov*. 2013;3(6):648-657. [PMID 23550148](https://pubmed.ncbi.nlm.nih.gov/23550148/) — succinate → PHD inhibition → HIF-1α/VEGF; the pseudohypoxia cluster and the high `WVASC_SDH` vascular dependence.

## 7. Clinical Pharmacology, Exposure–Response & TDM

32. Demetri GD, Wang Y, Wehrle E, et al. Imatinib plasma levels are correlated with clinical benefit in patients with unresectable/metastatic gastrointestinal stromal tumors. *J Clin Oncol*. 2009;27(19):3141-3147. [PMID 19451435](https://pubmed.ncbi.nlm.nih.gov/19451435/) — the ~1,100 ng/mL steady-state trough threshold; scenario 5/6.
33. Widmer N, Decosterd LA, Csajka C, et al. Population pharmacokinetics of imatinib and the role of alpha-acid glycoprotein. *Br J Clin Pharmacol*. 2006;62(1):97-112. [PMID 16842382](https://pubmed.ncbi.nlm.nih.gov/16842382/) — the AAG binding term; high tumour burden raises AAG and lowers free drug at unchanged total concentration.
34. Eechoute K, Fransson MN, Reyners AK, et al. A long-term prospective population pharmacokinetic study on imatinib plasma concentrations in GIST patients. *Clin Cancer Res*. 2012;18(20):5780-5787. [PMID 22850565](https://pubmed.ncbi.nlm.nih.gov/22850565/) — the time-dependent rise in apparent clearance implemented as `FCLDRIFT`.
35. Widmer N, Decosterd LA, Csajka C, et al. Therapeutic drug monitoring of imatinib — new data strengthen the case. *Clin Cancer Res*. 2012;18(19):5164-5166. [PMID 22904104](https://pubmed.ncbi.nlm.nih.gov/22904104/).

## 8. Treatment Interruption, Adjuvant Therapy & Surgery

36. Blay JY, Le Cesne A, Ray-Coquard I, et al. Prospective multicentric randomized phase III study of imatinib in patients with advanced gastrointestinal stromal tumors comparing interruption versus continuation of treatment beyond 1 year (BFR14). *J Clin Oncol*. 2007;25(9):1107-1113. [PMID 17369574](https://pubmed.ncbi.nlm.nih.gov/17369574/) — the rapid regrowth on interruption that the quiescent-persister compartment reproduces.
37. Le Cesne A, Ray-Coquard I, Bui BN, et al. Discontinuation of imatinib in patients with advanced gastrointestinal stromal tumours after 3 years of treatment: an open-label multicentre randomised phase 3 trial. *Lancet Oncol*. 2010;11(10):942-949. [PMID 20864406](https://pubmed.ncbi.nlm.nih.gov/20864406/).
38. DeMatteo RP, Ballman KV, Antonescu CR, et al. Adjuvant imatinib mesylate after resection of localised, primary gastrointestinal stromal tumour: a randomised, double-blind, placebo-controlled trial (ACOSOG Z9001). *Lancet*. 2009;373(9669):1097-1104. [PMID 19303137](https://pubmed.ncbi.nlm.nih.gov/19303137/).
39. Joensuu H, Eriksson M, Sundby Hall K, et al. One vs three years of adjuvant imatinib for operable gastrointestinal stromal tumor: a randomized trial (SSGXVIII/AIO). *JAMA*. 2012;307(12):1265-1272. [PMID 22453568](https://pubmed.ncbi.nlm.nih.gov/22453568/) — the "delay, not cure" adjuvant shape reproduced by scenario 12.
40. Raut CP, Posner M, Desai J, et al. Surgical management of advanced gastrointestinal stromal tumors after treatment with targeted systemic therapy using kinase inhibitors. *J Clin Oncol*. 2006;24(15):2325-2331. [PMID 16710031](https://pubmed.ncbi.nlm.nih.gov/16710031/) — debulking a focally resistant nodule.

## 9. Imaging & Molecular Response Assessment

41. Choi H, Charnsangavej C, Faria SC, et al. Correlation of computed tomography and positron emission tomography in patients with metastatic gastrointestinal stromal tumor treated at a single institution with imatinib mesylate: proposal of new computed tomography response criteria. *J Clin Oncol*. 2007;25(13):1753-1759. [PMID 17470865](https://pubmed.ncbi.nlm.nih.gov/17470865/) — the Choi density criteria; the model's `HU`/`MYX` compartments exist to reproduce why RECIST under-calls GIST response.
42. Van den Abbeele AD. The lessons of GIST — PET and PET/CT: a new paradigm for imaging. *Oncologist*. 2008;13(Suppl 2):8-13. [PMID 18434632](https://pubmed.ncbi.nlm.nih.gov/18434632/) — metabolic response within 24-48 h, the fastest readout in the model.
43. Namløs HM, Boye K, Mishkin SJ, et al. Noninvasive detection of ctDNA reveals intratumor heterogeneity and is associated with tumor burden in gastrointestinal stromal tumor. *Mol Cancer Ther*. 2018;17(11):2473-2480. [PMID 30097488](https://pubmed.ncbi.nlm.nih.gov/30097488/) — clone-resolved ctDNA, the basis of the `CTDS`/`CTDR` split.
44. Bhalla S, Boye K, Namløs HM, et al. Cell-free DNA in blood as a noninvasive insight into the sarcoma genome. *Mol Aspects Med*. 2020;72:100827. [PMID 31703948](https://pubmed.ncbi.nlm.nih.gov/31703948/).

## 10. Tumour Microenvironment, Immunity & Toxicity

45. Balachandran VP, Cavnar MJ, Zeng S, et al. Imatinib potentiates antitumor T cell responses in gastrointestinal stromal tumor through the inhibition of Ido. *Nat Med*. 2011;17(9):1094-1100. [PMID 21873989](https://pubmed.ncbi.nlm.nih.gov/21873989/) — the KIT-inhibition → IDO ↓ → CD8 ↑ arm of the immune cluster.
46. Desai J, Yassa L, Marqusee E, et al. Hypothyroidism after sunitinib treatment for patients with gastrointestinal stromal tumors. *Ann Intern Med*. 2006;145(9):660-664. [PMID 17088579](https://pubmed.ncbi.nlm.nih.gov/17088579/) — the `TSH` compartment.

---

## How the references map onto the model

| Model element | Primary anchors |
|---|---|
| Six-clone tumor structure | refs 13, 14, 15, 16 |
| Drug × clone IC50 matrix | refs 9, 15, 16, 24, 28, 29 |
| Imatinib PK, AAG, clearance drift | refs 32, 33, 34, 35 |
| Exon 9 dose effect (`GF_EX9`) | refs 9, 10, 11, 12 |
| Quiescent persister pool (`TQ`) | refs 19, 36, 37 |
| Vascular support (`ANG`, `WVASC_SDH`) | refs 21, 22, 30, 31 |
| SUV / CT density / RECIST ordering | refs 41, 42 |
| Clone-resolved ctDNA | refs 27, 43, 44 |
| INTRIGUE subgroup reversal | refs 15, 16, 24, 26, 27 |
| Adjuvant "delay, not cure" | refs 38, 39 |
| BIM polymorphism | ref 20 |
| Toxicity terms | refs 21, 22, 25, 46 |
