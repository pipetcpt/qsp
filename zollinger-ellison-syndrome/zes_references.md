# Zollinger–Ellison Syndrome / Gastrinoma — QSP Model References

62 PubMed-indexed references supporting the mechanistic map (`zes_qsp_model.dot`),
the mrgsolve model (`zes_mrgsolve_model.R`) and the Shiny dashboard
(`zes_shiny_app.R`).

Every PMID below was resolved against the NCBI E-utilities `esummary` endpoint
while this model was being built, so first author, journal and year are as
PubMed records them. Where a reference supplies a number the model actually
uses, that number is stated next to it, together with the parameter or anchor
it feeds.

---

## 1. Disease definition, cohorts and natural history

The quantitative backbone of this model is the National Institutes of Health
prospective ZES cohort, which is the only series large enough to give
distributions rather than case reports for fasting gastrin and basal acid
output.

1. **Berna MJ et al.** *Serum gastrin in Zollinger-Ellison syndrome: I.
   Prospective study of fasting serum gastrin in 309 patients from the National
   Institutes of Health and comparison with 2229 cases from the literature.*
   Medicine (Baltimore) 2006.
   [PMID 17108778](https://pubmed.ncbi.nlm.nih.gov/17108778/)
   → **anchor A4** (untreated fasting serum gastrin ≈ 900 pg/mL; the model's
   `KSECT0` is fitted to this). Also the source of the observation that
   fasting gastrin is *not* > 1000 pg/mL in most patients, which is why the
   model treats the acid brake and not the gastrin level as diagnostic.

2. **Berna MJ et al.** *Serum gastrin in Zollinger-Ellison syndrome: II.
   Prospective study of gastrin provocative testing in 293 patients from the
   National Institutes of Health and comparison with 537 cases from the
   literature.* Medicine (Baltimore) 2006.
   [PMID 17108779](https://pubmed.ncbi.nlm.nih.gov/17108779/)
   → the secretin-test Δ > 120 pg/mL criterion used in `secretin_test()` and
   validation item **A7**.

3. **Gibril F et al.** *Multiple endocrine neoplasia type 1 and
   Zollinger-Ellison syndrome: a prospective study of 107 cases and comparison
   with 1009 cases from the literature.* Medicine (Baltimore) 2004.
   [PMID 14747767](https://pubmed.ncbi.nlm.nih.gov/14747767/)
   → MEN1 accounts for 20–25% of ZES; multiplicity and duodenal predominance
   of MEN1 gastrinomas (`ZES_MEN1` phenotype: smaller `TUM0`, lower `GRADEF`).

4. **Gibril F et al.** *Prospective study of the natural history of gastrinoma
   in patients with MEN1: definition of an aggressive and a nonaggressive
   form.* J Clin Endocrinol Metab 2001.
   [PMID 11701693](https://pubmed.ncbi.nlm.nih.gov/11701693/)
   → the bimodal growth behaviour that motivates `GRADEF` as a phenotype input
   rather than a fixed constant.

5. **Gibril F, Jensen RT.** *Zollinger-Ellison syndrome revisited: diagnosis,
   biologic markers, associated inherited disorders, and acid hypersecretion.*
   Curr Gastroenterol Rep 2004.
   [PMID 15527675](https://pubmed.ncbi.nlm.nih.gov/15527675/)

6. **Chatzipanagiotou O et al.** *Gastrinoma and Zollinger-Ellison syndrome: A
   thorough update.* J Neuroendocrinol 2023.
   [PMID 37042078](https://pubmed.ncbi.nlm.nih.gov/37042078/)
   → contemporary diagnostic algorithm; the requirement for gastric pH < 2
   alongside hypergastrinaemia, which the model reproduces as an emergent
   discriminator (see `A7`).

7. **Qiu C et al.** *Clinicopathological characteristics and prognosis analysis
   of gastrinoma based on the SEER database.* PLoS One 2025.
   [PMID 41406116](https://pubmed.ncbi.nlm.nih.gov/41406116/)
   → survival stratification by grade and metastatic status; supports treating
   tumour burden, not acid, as the survival driver once acid is controlled.

8. **Pellicano R et al.** *Zollinger-Ellison syndrome in 2006: concepts from a
   clinical point of view.* Panminerva Med 2006.
   [PMID 16633330](https://pubmed.ncbi.nlm.nih.gov/16633330/)

9. **Wilcox CM, Seay T, Arcury JT, et al.** *Treatment strategies for
   Zollinger-Ellison syndrome.* Expert Opin Pharmacother 2009.
   [PMID 19351273](https://pubmed.ncbi.nlm.nih.gov/19351273/)
   → the "titrate until basal acid output is below 10 mEq/h (below 5 after
   previous gastric surgery)" algorithm implemented in `titrate_ppi()`.

10. **Tomassetti P et al.** *Optimal treatment of Zollinger-Ellison syndrome
    and related conditions in elderly patients.* Drugs Aging 2003.
    [PMID 14651442](https://pubmed.ncbi.nlm.nih.gov/14651442/)

11. **Riff BP et al.** *Weight Gain in Zollinger-Ellison Syndrome After Acid
    Suppression.* Pancreas 2016.
    [PMID 26164604](https://pubmed.ncbi.nlm.nih.gov/26164604/)
    → weight loss in untreated ZES is acid-driven (lipase inactivation and
    secretory diarrhoea), which is how `DIARR` and the steatorrhoea path are
    wired: to acid output, not to tumour burden.

---

## 2. Genetics: MEN1, menin, and the parathyroid–gut axis

12. **Adam MP et al. (GeneReviews).** *Multiple Endocrine Neoplasia Type 1.*
    [PMID 20301710](https://pubmed.ncbi.nlm.nih.gov/20301710/)
    → > 90% penetrance of primary hyperparathyroidism in MEN1, usually the
    first manifestation (`GLANDM` as the MEN1 driver).

13. **Huang J et al.** *The same pocket in menin binds both MLL and JUND but
    has opposite effects on transcription.* Nature 2012.
    [PMID 22327296](https://pubmed.ncbi.nlm.nih.gov/22327296/)

14. **Dreijerink KMA et al.** *Multi-omics analyses of MEN1 missense mutations
    identify disruption of menin-MLL and menin-JunD interactions as critical
    requirements for molecular pathogenicity.* Epigenetics Chromatin 2022.
    [PMID 35941657](https://pubmed.ncbi.nlm.nih.gov/35941657/)

15. **Agarwal SK.** *Multiple endocrine neoplasia type 1.* Front Horm Res 2013.
    [PMID 23652667](https://pubmed.ncbi.nlm.nih.gov/23652667/)

16. **Hackeng WM et al.** *A Parathyroid-Gut Axis: Hypercalcemia and the
    Pathogenesis of Gastrinoma in Multiple Endocrine Neoplasia 1.*
    Mol Cancer Res 2021.
    [PMID 33771883](https://pubmed.ncbi.nlm.nih.gov/33771883/)
    → **the reference behind cluster 19 and validation item A8.** This is the
    evidence that the calcium–gastrin coupling is causal rather than
    coincidental, and therefore that `ECA` (CaSR secretagogue gain) and `ACAW`
    (the calcium arm of `DRIVE`) belong in the model as edges.

17. **Sundaresan S et al.** *Deletion of Men1 and somatostatin induces
    hypergastrinemia and gastric carcinoids.* Gut 2017.
    [PMID 26860771](https://pubmed.ncbi.nlm.nih.gov/26860771/)
    → genetic demonstration that losing the somatostatin brake *is* what
    produces unrestrained gastrin; the model's "missing edge" formulation.

18. **van Beek DJ et al.** *Prognostic factors and survival in MEN1 patients
    with gastrinomas: Results from the DutchMEN study group (DMSG).*
    J Surg Oncol 2019.
    [PMID 31401809](https://pubmed.ncbi.nlm.nih.gov/31401809/)

19. **Hoffmann KM et al.** *Patients with multiple endocrine neoplasia type 1
    with gastrinomas have an increased risk of severe esophageal disease
    including stricture and the premalignant condition of Barrett's
    esophagus.* J Clin Endocrinol Metab 2006.
    [PMID 16249283](https://pubmed.ncbi.nlm.nih.gov/16249283/)
    → the `ESOPH` → stricture / Barrett path, and the reason MEN1-ZES is run as
    its own long-horizon scenario (24).

---

## 3. Gastrin: biosynthesis, circulating species, kinetics, feedback

20. **Dimaline R, Varro A.** *Novel roles of gastrin.* J Physiol 2014.
    [PMID 24665102](https://pubmed.ncbi.nlm.nih.gov/24665102/)

21. **Blair AJ 3rd et al.** *Comparison of acid secretory responsiveness to
    gastrin heptadecapeptide and of gastrin heptadecapeptide pharmacokinetics
    in duodenal ulcer patients and normal subjects.* J Clin Invest 1986.
    [PMID 3745438](https://pubmed.ncbi.nlm.nih.gov/3745438/)
    → G17 disappearance kinetics; `KEL17` = 6.93 /h (t½ ≈ 6 min). The much
    slower G34 clearance (`KEL34` = 0.99 /h, t½ ≈ 42 min) is why G34 dominates
    *fasting immunoreactivity* while G17 dominates *bioactivity* — the reason
    the model carries two species and a potency weight (`GBIO = G17 + G34/6`)
    rather than one "gastrin" compartment.

22. **Bramante G et al.** *Ferric ions inhibit proteolytic processing of
    progastrin.* Biochem Biophys Res Commun 2011.
    [PMID 21195058](https://pubmed.ncbi.nlm.nih.gov/21195058/)

23. **Marshall KM et al.** *The C-terminal flanking peptide of progastrin
    induces gastric cell apoptosis and stimulates colonic cell division in
    vivo.* Peptides 2013.
    [PMID 23742999](https://pubmed.ncbi.nlm.nih.gov/23742999/)

24. **Zavros Y et al.** *Hypergastrinemia in response to gastric inflammation
    suppresses somatostatin.* Am J Physiol Gastrointest Liver Physiol 2002.
    [PMID 11751171](https://pubmed.ncbi.nlm.nih.gov/11751171/)
    → the D-cell somatostatin brake and its acid dependence; `HBRK50`,
    `HBRKN`.

25. **Schmidt PT et al.** *Cholecystokinin inhibits gastrin secretion
    independently of paracrine somatostatin secretion in the pig.*
    Scand J Gastroenterol 2004.
    [PMID 15074389](https://pubmed.ncbi.nlm.nih.gov/15074389/)

26. **Bevilacqua M et al.** *Dissimilar PTH, gastrin, and calcitonin responses
    to oral calcium and peptones in hypocalciuric hypercalcemia, primary
    hyperparathyroidism, and normal subjects.* J Bone Miner Res 2006.
    [PMID 16491288](https://pubmed.ncbi.nlm.nih.gov/16491288/)
    → quantitative basis for `ECA`, the calcium sensitivity of gastrin release.

27. **Tzaneva M.** *Effects of duodenogastric reflux on gastrin cells,
    somatostatin cells and serotonin cells in human antral gastric mucosa.*
    Pathol Res Pract 2004.
    [PMID 15310146](https://pubmed.ncbi.nlm.nih.gov/15310146/)

---

## 4. Trophic action of gastrin: ECL cells and parietal cell mass

The key structural claim — that gastrin's trophic action is *selective* — comes
from reference 28, and it is why `EMAXECL` (3.0) exceeds `EMAXPCM` (2.2) in the
model and why the histamine amplifier, not direct parietal CCK2R signalling,
carries most of the drive (`AH` = 6.0 versus `AG` = 3.5).

28. **Bakke I et al.** *Gastrin has a specific proliferative effect on the rat
    enterochromaffin-like cell, but not on the parietal cell: a study by
    elutriation centrifugation.* Acta Physiol Scand 2000.
    [PMID 10759608](https://pubmed.ncbi.nlm.nih.gov/10759608/)

29. **Modlin IM, Kidd M, Tang LH.** *Pathophysiology of the fundic
    enterochromaffin-like (ECL) cell and gastric carcinoid tumours.*
    Ann R Coll Surg Engl 1996.
    [PMID 8678447](https://pubmed.ncbi.nlm.nih.gov/8678447/)

30. **Freston JW et al.** *Effects of hypochlorhydria and hypergastrinemia on
    structure and function of gastrointestinal cells. A review and analysis.*
    Dig Dis Sci 1995.
    [PMID 7859584](https://pubmed.ncbi.nlm.nih.gov/7859584/)
    → time-scales for `KPCM` (t½ ≈ 3 weeks) and `KECL` (t½ ≈ 2 weeks).

31. **Chen D et al.** *Ultrastructure of enterochromaffin-like cells in rat
    stomach: effects of alpha-fluoromethylhistidine-evoked histamine depletion
    and hypergastrinemia.* Cell Tissue Res 1996.
    [PMID 8593676](https://pubmed.ncbi.nlm.nih.gov/8593676/)

32. **Zheng B et al.** *Gastrin-dependent expansion of Cck2r+ corpus
    progenitors accelerates ulcer healing and inhibits gastric dysplasia.*
    Gut 2026.
    [PMID 40983503](https://pubmed.ncbi.nlm.nih.gov/40983503/)

33. **Smith JP, Nadella S, Osborne N.** *Gastrin and Gastric Cancer.*
    Cell Mol Gastroenterol Hepatol 2017.
    [PMID 28560291](https://pubmed.ncbi.nlm.nih.gov/28560291/)

34. **Massironi S et al.** *The Autoimmune Gastritis Puzzle: Emerging Cellular
    Crosstalk and Molecular Pathways Driving Parietal Cell Loss and ECL Cell
    Hyperplasia.* Cells 2025.
    [PMID 41148791](https://pubmed.ncbi.nlm.nih.gov/41148791/)
    → the mirror-image disorder, and the reason `ATRF` exists as a parameter:
    the same acid-brake edge run in the opposite direction.

---

## 5. Parietal cell physiology and the H+/K+-ATPase lifecycle

Cluster 8 of the map and the four-state pump pool (`PUMPI`/`PUMPA`/`PUMPB`/
`PUMPR`) rest on this block.

35. **Engevik AC, Kaji I, Goldenring JR.** *The Physiology of the Gastric
    Parietal Cell.* Physiol Rev 2020.
    [PMID 31670611](https://pubmed.ncbi.nlm.nih.gov/31670611/)
    → tubulovesicle-to-canalicular trafficking as the physiological switch;
    `KACT0`, `KACT1`, `KINACT`.

36. **Shin JM, Munson K, Sachs G.** *Gastric H+,K+-ATPase.*
    Compr Physiol 2011.
    [PMID 23733700](https://pubmed.ncbi.nlm.nih.gov/23733700/)

37. **Sachs G, Shin JM, Vagin O, et al.** *The gastric H,K ATPase as a drug
    target: past, present, and future.* J Clin Gastroenterol 2007.
    [PMID 17575528](https://pubmed.ncbi.nlm.nih.gov/17575528/)

38. **Shin JM, Sachs G.** *Restoration of acid secretion following treatment
    with proton pump inhibitors.* Gastroenterology 2002.
    [PMID 12404233](https://pubmed.ncbi.nlm.nih.gov/12404233/)
    → **the single most load-bearing reference in the file.** Pump protein
    half-life ≈ 50 h sets `KDEGP` = 0.01386 /h, and the *de novo* synthesis
    rate is what makes the dosing interval, rather than the daily dose, the
    thing that fails first in ZES.

39. **Shin JM, Cho YM, Sachs G.** *Chemistry of covalent inhibition of the
    gastric (H+, K+)-ATPase by proton pump inhibitors.* J Am Chem Soc 2004.
    [PMID 15212527](https://pubmed.ncbi.nlm.nih.gov/15212527/)
    → acid-catalysed conversion of the prodrug to the sulfenamide, i.e. the
    acid-activation trap encoded as
    `dPPICAN/dt = KACTP·C·PUMPA/(PUMPA + KMPA) − KOUTCAN·PPICAN`.

40. **Besancon M et al.** *Sites of reaction of the gastric H,K-ATPase with
    extracytoplasmic thiol reagents.* J Biol Chem 1997.
    [PMID 9278394](https://pubmed.ncbi.nlm.nih.gov/9278394/)
    → Cys813/Cys892 adduct formation; irreversibility, hence `PUMPB` as a
    separate pool that only leaves by degradation (`KDEGP`) or slow glutathione
    reduction (`KGSH`).

41. **Shin JM, Kim N.** *Pharmacology of proton pump inhibitors.*
    Curr Gastroenterol Rep 2008.
    [PMID 19006606](https://pubmed.ncbi.nlm.nih.gov/19006606/)

42. **Niyomugabo AP et al.** *Gastric Proton Pumpopathy Associated with Protein
    Sorting Machinery: A Narrative Review.* Clin Exp Gastroenterol 2026.
    [PMID 41889561](https://pubmed.ncbi.nlm.nih.gov/41889561/)

---

## 6. PPI, P-CAB and H2RA pharmacokinetics / pharmacodynamics

43. **Li S et al.** *Prediction of Omeprazole Pharmacokinetics and its
    Inhibition on Gastric Acid Secretion in Humans Using Physiologically Based
    Pharmacokinetic-Pharmacodynamic Modeling.* Pharm Res 2023.
    [PMID 37226024](https://pubmed.ncbi.nlm.nih.gov/37226024/)
    → the closest published analogue of this file's PPI submodel; source of
    `VPPI`, `CLPPI0` and **anchor A5** (≈ 66% acid inhibition at steady state
    on omeprazole 20 mg once daily).

44. **Na JY et al.** *Influence of CYP2C19 Polymorphisms on the
    Pharmacokinetics of Omeprazole in Elderly Subjects.*
    Clin Pharmacol Drug Dev 2021.
    [PMID 34337876](https://pubmed.ncbi.nlm.nih.gov/34337876/)
    → the `CYPF` grid (UM 1.9 / NM 1.0 / IM 0.5 / PM 0.28) used in scenarios
    10–11 and validation item **A5**.

45. **Kagami T et al.** *Potent acid inhibition by vonoprazan in comparison
    with esomeprazole, with reference to CYP2C19 genotype.*
    Aliment Pharmacol Ther 2016.
    [PMID 26991399](https://pubmed.ncbi.nlm.nih.gov/26991399/)
    → the pH > 4 holding-time contrast the model must reproduce, and the
    evidence that a P-CAB is far less CYP2C19-sensitive (`VON2C19` in the map,
    `CLVON` route in the model).

46. **Zerbetto De Palma G et al.** *Protonation of Key Acidic Residues Reveals
    Binding Features of PCABs to Gastric H,K-ATPase.* J Membr Biol 2026.
    [PMID 42217047](https://pubmed.ncbi.nlm.nih.gov/42217047/)
    → K+-competitive, reversible binding without acid activation, and the
    accessibility of *resting* pumps that `FRESTV` encodes.

47. **Wehbe H, Fass R.** *Potassium-competitive acid blockers: Clinical pearls
    and practical insights.* Cleve Clin J Med 2026.
    [PMID 42225381](https://pubmed.ncbi.nlm.nih.gov/42225381/)

48. **Yoon DY et al.** *Effect of meal timing on pharmacokinetics and
    pharmacodynamics of tegoprazan in healthy male volunteers.*
    Clin Transl Sci 2021.
    [PMID 33382926](https://pubmed.ncbi.nlm.nih.gov/33382926/)
    → the meal-timing dependence a PPI has and a P-CAB largely does not; in
    the model this difference is not a rule but a consequence of `FRESTV` and
    of the `PUMPA/(PUMPA+KMPA)` activation term.

49. **Galmiche JP et al.** *Tenatoprazole, a novel proton pump inhibitor with a
    prolonged plasma half-life: effects on intragastric pH and comparison with
    esomeprazole in healthy volunteers.* Aliment Pharmacol Ther 2004.
    [PMID 15023167](https://pubmed.ncbi.nlm.nih.gov/15023167/)
    → the cleanest demonstration that prolonging plasma exposure, not raising
    the dose, is what fixes the trough.

50. **Weigt J et al.** *Nocturnal gastric acid breakthrough is not associated
    with night-time gastroesophageal reflux in GERD patients.* Dig Dis 2009.
    [PMID 19439964](https://pubmed.ncbi.nlm.nih.gov/19439964/)
    → definition of nocturnal acid breakthrough used by `nab_hours()`.

51. **Gillen D, McColl KE.** *Problems related to acid rebound and
    tachyphylaxis.* Best Pract Res Clin Gastroenterol 2001.
    [PMID 11403541](https://pubmed.ncbi.nlm.nih.gov/11403541/)

52. **Clark JH et al.** *Histamine 2-Receptor Antagonists Tachyphylaxis: A
    Scoping Review.* Laryngoscope 2026.
    [PMID 41459837](https://pubmed.ncbi.nlm.nih.gov/41459837/)
    → `RH2`, `TACHM`, `KTACH`: H2-receptor up-regulation under blockade, and
    hence the escape seen in scenarios 14–15.

53. **Namikawa K et al.** *Rebound Acid Hypersecretion after Withdrawal of
    Long-Term Proton Pump Inhibitor (PPI) Treatment — Are PPIs Addictive?*
    Int J Mol Sci 2024.
    [PMID 38791497](https://pubmed.ncbi.nlm.nih.gov/38791497/)
    → validation item **A9**. In the model rebound is not a parameter: it is
    the enlarged `PUMPI` pool discharging onto the enlarged `PCM`.

---

## 7. Acid control in ZES: the therapeutic trials the model is fitted to

54. **Metz DC et al.** *Three-year oral pantoprazole administration is
    effective for patients with Zollinger-Ellison syndrome and other
    hypersecretory conditions.* Aliment Pharmacol Ther 2006.
    [PMID 16423003](https://pubmed.ncbi.nlm.nih.gov/16423003/)
    → **anchor A6**: the maintenance regimen that holds basal acid output
    below 10 mEq/h, and the observation that most patients need twice-daily
    dosing at a total dose well above the reflux dose.

55. **Metz DC et al.** *Maintenance oral pantoprazole therapy is effective for
    patients with Zollinger-Ellison syndrome and idiopathic hypersecretion.*
    Am J Gastroenterol 2003.
    [PMID 12591045](https://pubmed.ncbi.nlm.nih.gov/12591045/)

56. **Hirschowitz BI, Simmons J, Mohnen J.** *Long-term lansoprazole control of
    gastric acid and pepsin secretion in ZE and non-ZE hypersecretors: a
    prospective 10-year study.* Aliment Pharmacol Ther 2001.
    [PMID 11683694](https://pubmed.ncbi.nlm.nih.gov/11683694/)
    → the ten-year stability of the required dose, which the model reproduces
    as a consequence of `PCM` and `ECL` reaching a plateau rather than of any
    explicit tolerance term.

57. **Hirschowitz BI, Mohnen J, Shaw S.** *Long-term treatment with
    lansoprazole for patients with Zollinger-Ellison syndrome.*
    Aliment Pharmacol Ther 1996.
    [PMID 8853754](https://pubmed.ncbi.nlm.nih.gov/8853754/)

58. **Corleto V et al.** *Efficacy of long-term therapy with low doses of
    omeprazole in the control of gastric acid secretion in Zollinger-Ellison
    syndrome patients.* Aliment Pharmacol Ther 1993.
    [PMID 8485270](https://pubmed.ncbi.nlm.nih.gov/8485270/)
    → the counter-observation that a substantial minority are controlled on
    low doses; the model's `A13` virtual population is the mechanism it
    proposes (variability in the two trophic gains, not in drug handling).

59. **Lew EA et al.** *Intravenous pantoprazole rapidly controls gastric acid
    hypersecretion in patients with Zollinger-Ellison syndrome.*
    Gastroenterology 2000.
    [PMID 10734021](https://pubmed.ncbi.nlm.nih.gov/10734021/)
    → `ev_ppi_iv()`.

60. **Pisegna JR et al.** *Inhibition of pentagastrin-induced gastric acid
    secretion by intravenous pantoprazole: a dose-response study.*
    Am J Gastroenterol 1999.
    [PMID 10520836](https://pubmed.ncbi.nlm.nih.gov/10520836/)
    → the pentagastrin-stimulation protocol reproduced by `mao_of()`.

61. **Saeed ZA et al.** *Parenteral antisecretory drug therapy in patients with
    Zollinger-Ellison syndrome.* Gastroenterology 1989.
    [PMID 2565842](https://pubmed.ncbi.nlm.nih.gov/2565842/)
    → the historical H2-antagonist doses (famotidine up to several hundred
    milligrams a day) that scenario 15 tests.

62. **Hirschowitz BI, Simmons J, Mohnen J.** *Risk factors for esophagitis in
    extreme acid hypersecretors with and without Zollinger-Ellison syndrome.*
    Clin Gastroenterol Hepatol 2004.
    [PMID 15017606](https://pubmed.ncbi.nlm.nih.gov/15017606/)

63. **Strader DB et al.** *Esophageal function and occurrence of Barrett's
    esophagus in Zollinger-Ellison syndrome.* Digestion 1995.
    [PMID 8549876](https://pubmed.ncbi.nlm.nih.gov/8549876/)

---

## 8. Diagnosis, and the false positive the model reproduces

64. **Bhattacharya S et al.** *Validity of Secretin Stimulation Testing on
    Proton Pump Inhibitor Therapy for Diagnosis of Zollinger-Ellison
    Syndrome.* Am J Gastroenterol 2021.
    [PMID 34515664](https://pubmed.ncbi.nlm.nih.gov/34515664/)
    → the clinical statement of the problem that validation item **A7** poses
    mechanistically: a PPI raises gastrin by lifting the D-cell acid brake, so
    a healthy subject on a PPI can look like a gastrinoma.

65. **Proye C et al.** *Intraoperative gastrin measurements during surgical
    management of patients with gastrinomas: experience with 20 cases.*
    World J Surg 1998.
    [PMID 9606276](https://pubmed.ncbi.nlm.nih.gov/9606276/)

---

## 9. Surgery

66. **Norton JA et al.** *Does the use of routine duodenotomy (DUODX) affect
    rate of cure, development of liver metastases, or survival in patients with
    Zollinger-Ellison syndrome?* Ann Surg 2004.
    [PMID 15082965](https://pubmed.ncbi.nlm.nih.gov/15082965/)
    → the cure rates behind scenario 18's `TUMCUT` state edit.

67. **Norton JA, Jensen RT.** *Surgical treatment and prognosis of
    gastrinoma.* Best Pract Res Clin Gastroenterol 2005.
    [PMID 16253901](https://pubmed.ncbi.nlm.nih.gov/16253901/)

68. **Norton JA.** *Surgery and prognosis of duodenal gastrinoma as a duodenal
    neuroendocrine tumor.* Best Pract Res Clin Gastroenterol 2005.
    [PMID 16253894](https://pubmed.ncbi.nlm.nih.gov/16253894/)

69. **Albers MB, Manoharan J, Bartsch DK.** *Contemporary surgical management
    of the Zollinger-Ellison syndrome in multiple endocrine neoplasia type 1.*
    Best Pract Res Clin Endocrinol Metab 2019.
    [PMID 31521501](https://pubmed.ncbi.nlm.nih.gov/31521501/)
    → the sequencing rule (parathyroids first) that scenario 17 tests without
    being told.

70. **Gaujoux S et al.** *Surgical Management of Zollinger-Ellison Syndrome in
    Multiple Endocrine Neoplasia Type 1: an AFCE and GTE Cohort Study.*
    World J Surg 2026.
    [PMID 41862416](https://pubmed.ncbi.nlm.nih.gov/41862416/)

71. **Kong W et al.** *Pancreaticoduodenectomy Is the Best Surgical Procedure
    for Zollinger-Ellison Syndrome Associated with Multiple Endocrine Neoplasia
    Type 1.* Cancers (Basel) 2022.
    [PMID 35454834](https://pubmed.ncbi.nlm.nih.gov/35454834/)

72. **Krampitz GW, Norton JA.** *Current management of the Zollinger-Ellison
    syndrome.* Adv Surg 2013.
    [PMID 24298844](https://pubmed.ncbi.nlm.nih.gov/24298844/)

---

## 10. Antitumour therapy: the four systemic arms

73. **Kunz PL et al.** *Randomized Study of Temozolomide or Temozolomide and
    Capecitabine in Patients With Advanced Pancreatic Neuroendocrine Tumors
    (ECOG-ACRIN E2211).* J Clin Oncol 2023.
    [PMID 36260828](https://pubmed.ncbi.nlm.nih.gov/36260828/)
    → the CAPTEM arm (scenario 30); median progression-free survival 22.7
    months for temozolomide + capecitabine.

74. **Evans MG et al.** *Loss of O6-Methylguanine-DNA Methyltransferase Protein
    Expression by Immunohistochemistry Is Associated With Response to
    Capecitabine and Temozolomide.* World J Surg 2025.
    [PMID 39825572](https://pubmed.ncbi.nlm.nih.gov/39825572/)
    → `MGMTM` and `EMGMT`: the model's alkylator kill is multiplied by MGMT
    status rather than assumed uniform.

75. **Fazio N et al.** *Relationship between metabolic toxicity and efficacy of
    everolimus in patients with neuroendocrine tumors: A pooled analysis from
    the randomized, phase 3 RADIANT-3 and RADIANT-4 trials.* Cancer 2021.
    [PMID 33857327](https://pubmed.ncbi.nlm.nih.gov/33857327/)
    → the everolimus arm (scenario 28) and its `EVEADV` toxicity node.

76. **Chan DL et al.** *Markers of Systemic Inflammation in Neuroendocrine
    Tumors: A Pooled Analysis of the RADIANT-3 and RADIANT-4 Studies.*
    Pancreas 2021.
    [PMID 33560090](https://pubmed.ncbi.nlm.nih.gov/33560090/)

77. **Baudin E et al.** *[177Lu]Lu-DOTA-TATE versus sunitinib in patients with
    metastatic progressive neuroendocrine tumours of the pancreas
    (OCLURANDOM): a randomised, controlled, phase 2 trial.*
    Lancet Oncol 2026.
    [PMID 42225102](https://pubmed.ncbi.nlm.nih.gov/42225102/)
    → **the head-to-head that scenarios 29 and 31 reproduce**, and the reason
    those two arms are run from an identical baseline state.

78. **Strosberg JR et al.** *177Lu-Dotatate plus long-acting octreotide versus
    high-dose long-acting octreotide in patients with midgut neuroendocrine
    tumours (NETTER-1): final overall survival and long-term safety results.*
    Lancet Oncol 2021.
    [PMID 34793718](https://pubmed.ncbi.nlm.nih.gov/34793718/)
    → the 7.4 GBq × 4 every-8-weeks schedule in `ev_prrt()`.

79. **Bodei L et al.** *Dosimetry of [177Lu]Lu-DOTATATE in Patients with
    Advanced Midgut Neuroendocrine Tumors: Results from a Substudy of the Phase
    III NETTER-1 Trial.* J Nucl Med 2025.
    [PMID 39947918](https://pubmed.ncbi.nlm.nih.gov/39947918/)
    → `SDOSET`, `SDOSEK`: absorbed dose per unit of retained activity, and the
    23 Gy renal constraint that scenario 33 breaches when amino-acid protection
    is withheld.

80. **Kobayashi N et al.** *Safety and efficacy of peptide receptor
    radionuclide therapy with 177Lu-DOTA0-Tyr3-octreotate in combination with
    amino acid solution infusion in Japanese patients with somatostatin
    receptor positive tumours.* Ann Nucl Med 2021.
    [PMID 34533700](https://pubmed.ncbi.nlm.nih.gov/34533700/)
    → `AALYS`, `FPROT`.

81. **Park EA et al.** *The Impact of Radiopharmaceutical Therapy on Renal
    Function.* Semin Nucl Med 2022.
    [PMID 35314056](https://pubmed.ncbi.nlm.nih.gov/35314056/)

82. **Pavel M et al.** *Efficacy and safety of high-dose lanreotide autogel in
    patients with progressive pancreatic or midgut neuroendocrine tumours:
    CLARINET FORTE phase 2 study results.* Eur J Cancer 2021.
    [PMID 34597974](https://pubmed.ncbi.nlm.nih.gov/34597974/)

83. **Pusceddu S et al.** *Impact of Diabetes and Metformin Use on
    Enteropancreatic Neuroendocrine Tumors: Post Hoc Analysis of the CLARINET
    Study.* Cancers (Basel) 2021.
    [PMID 35008233](https://pubmed.ncbi.nlm.nih.gov/35008233/)

84. **Mujica-Mota R et al.** *Everolimus, lutetium-177 DOTATATE and sunitinib
    for advanced, unresectable or metastatic neuroendocrine tumours with
    disease progression: a systematic review and cost-effectiveness analysis.*
    Health Technol Assess 2018.
    [PMID 30209002](https://pubmed.ncbi.nlm.nih.gov/30209002/)
    → the cross-trial progression-free-survival table used to sanity-check the
    relative ordering of scenarios 27–31.

85. **Gaztambide S, Vazquez JA.** *Short- and long-term effect of a long-acting
    somatostatin analogue, lanreotide (SR-L) on metastatic gastrinoma.*
    J Endocrinol Invest 1999.
    [PMID 10195383](https://pubmed.ncbi.nlm.nih.gov/10195383/)
    → the antisecretory arm `ESSAS`: a somatostatin analogue lowers tumour
    gastrin output, which is the only licensed way to move factors 1 and 3.

86. **Plöckinger U et al.** *Effect of the somatostatin analog octreotide on
    gastric mucosal function and histology during 3 months of preoperative
    treatment in patients with acromegaly.* Eur J Endocrinol 1998.
    [PMID 9820614](https://pubmed.ncbi.nlm.nih.gov/9820614/)
    → `ESSAE`, `ESSAP`: direct SSTR2 effects on ECL mass and on the parietal
    cell, separate from the antisecretory effect on the tumour.

87. **Faggiano A et al.** *The antiproliferative effect of pasireotide LAR
    alone and in combination with everolimus in patients with medullary thyroid
    cancer.* Endocrine 2018.
    [PMID 29572709](https://pubmed.ncbi.nlm.nih.gov/29572709/)

88. **Sawicka-Gutaj N et al.** *Pasireotide — Mechanism of Action and Clinical
    Applications.* Curr Drug Metab 2018.
    [PMID 29595102](https://pubmed.ncbi.nlm.nih.gov/29595102/)

---

## 11. CCK2-receptor antagonism

89. **Boyce M et al.** *Netazepide, a gastrin/cholecystokinin-2 receptor
    antagonist, can eradicate gastric neuroendocrine tumours in patients with
    autoimmune chronic atrophic gastritis.* Br J Clin Pharmacol 2017.
    [PMID 27704617](https://pubmed.ncbi.nlm.nih.gov/27704617/)
    → `KINET`; scenario 19. This is the only agent in the file that blocks the
    gastrin signal itself rather than its consequences.

90. **Boyce M, Moore AR, Sagatun L, et al.** *Potential clinical indications
    for a CCK2 receptor antagonist.* Curr Opin Pharmacol 2016.
    [PMID 27710813](https://pubmed.ncbi.nlm.nih.gov/27710813/)

91. **Lloyd KA et al.** *Netazepide Inhibits Expression of Pappalysin 2 in Type
    1 Gastric Neuroendocrine Tumors.* Cell Mol Gastroenterol Hepatol 2020.
    [PMID 32004755](https://pubmed.ncbi.nlm.nih.gov/32004755/)

---

## 12. Gastric neuroendocrine tumours as a late consequence

92. **Lamberti G et al.** *Gastric neuroendocrine neoplasms.*
    Nat Rev Dis Primers 2024.
    [PMID 38605021](https://pubmed.ncbi.nlm.nih.gov/38605021/)
    → the type-1 / type-2 / type-3 taxonomy; `MENFLG` gates the type-2 path in
    the model because a type-2 gastric NET requires hypergastrinaemia *plus*
    the MEN1 background.

93. **McCarthy DM.** *Proton Pump Inhibitor Use, Hypergastrinemia, and Gastric
    Carcinoids — What Is the Relationship?* Int J Mol Sci 2020.
    [PMID 31963924](https://pubmed.ncbi.nlm.nih.gov/31963924/)

94. **Rais R et al.** *Enterochromaffin-like Cell Hyperplasia-Associated
    Gastric Neuroendocrine Tumors May Arise in the Setting of Proton Pump
    Inhibitor Use.* Arch Pathol Lab Med 2022.
    [PMID 34283890](https://pubmed.ncbi.nlm.nih.gov/34283890/)

95. **Massironi S et al.** *Gastric carcinoids: between underestimation and
    overtreatment.* World J Gastroenterol 2009.
    [PMID 19437556](https://pubmed.ncbi.nlm.nih.gov/19437556/)

96. **Moeinipour Y et al.** *The Association Between Proton Pump Inhibitor Use
    and Gastroenteropancreatic Neuroendocrine Tumors: A Systematic Review.*
    J Gastrointest Cancer 2025.
    [PMID 41191220](https://pubmed.ncbi.nlm.nih.gov/41191220/)

---

## 13. Long-term consequences of profound acid suppression

97. **Choudhury A et al.** *Vitamin B12 deficiency and use of proton pump
    inhibitors: a systematic review and meta-analysis.*
    Expert Rev Gastroenterol Hepatol 2023.
    [PMID 37060552](https://pubmed.ncbi.nlm.nih.gov/37060552/)
    → `PHB12`, `NB12`, `KB12OUT`.

98. **Mumtaz H et al.** *Association of Vitamin B12 deficiency with long-term
    PPIs use: A cohort study.* Ann Med Surg (Lond) 2022.
    [PMID 36268318](https://pubmed.ncbi.nlm.nih.gov/36268318/)

99. **Wang C et al.** *Safety of proton pump inhibitors: an overview of
    systematic reviews and meta-analyses.* BMJ Evid Based Med 2026.
    [PMID 42303374](https://pubmed.ncbi.nlm.nih.gov/42303374/)
    → the umbrella review that sets the *magnitude* of every harm in cluster
    21, and the reason `FMGFL` puts an acid-independent floor under magnesium
    absorption rather than letting it fall without limit.

100. **Sundar C et al.** *Impact of Proton Pump Inhibitor Therapy on Bone
     Mineral Density: An Updated Systematic Review.* Cureus 2026.
     [PMID 42110060](https://pubmed.ncbi.nlm.nih.gov/42110060/)
     → `KBMD` (about 3% bone loss over five years of near-complete
     achlorhydria).

101. **Li P et al.** *Proton pump inhibitors use and the risk of osteoporosis
     and fractures: A two-sample Mendelian randomization study.*
     Medicine (Baltimore) 2026.
     [PMID 42499118](https://pubmed.ncbi.nlm.nih.gov/42499118/)

102. **Azores-Moreno J et al.** *Acute Drug-Induced Tubulointerstitial
     Nephritis: Current Perspectives on Diagnosis and Treatment.*
     Adv Kidney Dis Health 2025.
     [PMID 40947149](https://pubmed.ncbi.nlm.nih.gov/40947149/)
     → the `AKIINT` → `RENALG` → gastrin-clearance loop, i.e. why renal injury
     raises *measured* gastrin without any change in secretion (scenario 35).

103. **Al-Fakhouri Z et al.** *Short-Term Daily Vonoprazan Treatment Is More
     Commonly Associated With Infectious Gastroenteritis Than Short-Term Daily
     PPI Treatment.* Neurogastroenterol Motil 2026.
     [PMID 42260324](https://pubmed.ncbi.nlm.nih.gov/42260324/)
     → the `ENTINF` node, and a caution attached to the model's own finding
     that a P-CAB gives the deepest acid suppression of the three classes.

104. **Uemura N et al.** *Vonoprazan as a Long-term Maintenance Treatment for
     Erosive Esophagitis: VISION, a 5-Year, Randomized, Open-label Study.*
     Clin Gastroenterol Hepatol 2025.
     [PMID 39209187](https://pubmed.ncbi.nlm.nih.gov/39209187/)

105. **Nafea S et al.** *Association of Long-Term Proton Pump Inhibitor Use
     With Nutrient Deficiencies: A Retrospective Cross-Sectional Study.*
     Cureus 2026.
     [PMID 42064495](https://pubmed.ncbi.nlm.nih.gov/42064495/)

---

## 14. Methodology

106. **Elmokadem A, Riggs MM, Baron KT.** mrgsolve — R package for simulation
     from ODE-based population PK/PD and QSP models.
     <https://mrgsolve.org/> (model built and validated under mrgsolve 2.0.1).

107. **gPKPDviz** — an mrgsolve-based Shiny tool for PK/PD simulation.
     Paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> ·
     Code: <https://github.com/Genentech/gPKPDviz/>

---

## Anchor table (what the calibration is fitted to)

| # | Quantity | Target | Source |
|---|----------|--------|--------|
| A1 | Healthy basal acid output | 3.0 mEq/h | refs 35, 60 |
| A2 | Healthy pentagastrin-stimulated maximal acid output | 23 mEq/h | ref 60 |
| A3 | Untreated sporadic ZES basal acid output | 36 mEq/h | refs 1, 5, 9 |
| A4 | Untreated sporadic ZES fasting serum gastrin | 900 pg/mL | ref 1 |
| A5 | Healthy, omeprazole 20 mg od, day 5: acid inhibition | 66% | ref 43 |
| A6 | ZES, omeprazole 30 mg bd, day 10: fasting acid output | 6 mEq/h | refs 54, 55 |

Held-out (not fitted, used only to test the model): the secretin-test
increment and its PPI false positive (refs 2, 64), the CYP2C19 gradient
(ref 44), the vonoprazan pH-holding-time advantage (ref 45), H2-antagonist
escape (refs 52, 61), withdrawal rebound (ref 53), the effect of
parathyroidectomy on gastrin and acid output (refs 16, 69), and the relative
ordering of the four antitumour arms (refs 73, 77, 78, 84).

---

*Educational and research use only. Not validated for clinical or regulatory
use.*
