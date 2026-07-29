# Peripartum Cardiomyopathy (PPCM) — Annotated References

References supporting the mechanistic map (`ppcm_qsp_model.dot`), the mrgsolve
ODE model (`ppcm_mrgsolve_model.R`) and the Shiny dashboard
(`ppcm_shiny_app.R`). Grouped by the role each paper plays in the model, with
the specific model element it justifies noted in *italics*.

---

## 1. Definition, epidemiology and diagnostic framework

1. Sliwa K, Hilfiker-Kleiner D, Petrie MC, et al. **Current state of knowledge on aetiology, diagnosis, management, and therapy of peripartum cardiomyopathy: a position statement from the Heart Failure Association of the ESC Working Group on peripartum cardiomyopathy.** *Eur J Heart Fail* 2010;12:767–78. — <https://pubmed.ncbi.nlm.nih.gov/20675664/>
   *The operational definition used here: LVEF <45%, last month of pregnancy to ~5 months postpartum, no other identifiable cause.*

2. Bauersachs J, König T, van der Meer P, et al. **Pathophysiology, diagnosis and management of peripartum cardiomyopathy: a position statement from the Heart Failure Association of the ESC Study Group on Peripartum Cardiomyopathy.** *Eur J Heart Fail* 2019;21:827–43. — <https://pubmed.ncbi.nlm.nih.gov/31243866/>
   *Source of the BOARD bundle (Bromocriptine, Oral heart-failure therapy, Anticoagulants, vasoRelaxing agents, Diuretics) and of the mandatory-anticoagulation rule modelled by `PROTHR_BRC` + `EAC_THR`.*

3. Arany Z, Elkayam U. **Peripartum cardiomyopathy.** *Circulation* 2016;133:1397–409. — <https://pubmed.ncbi.nlm.nih.gov/27045128/>
   *Integrative review used to structure clusters 1–10 of the mechanistic map.*

4. Davis MB, Arany Z, McNamara DM, Goland S, Elkayam U. **Peripartum cardiomyopathy: JACC state-of-the-art review.** *J Am Coll Cardiol* 2020;75:207–21. — <https://pubmed.ncbi.nlm.nih.gov/31948651/>
   *Contemporary management algorithm; source for the device/mechanical-support cluster (WCD, ICD deferral, IABP/ECMO, LVAD).*

5. Sliwa K, Petrie MC, van der Meer P, et al. **Clinical presentation, management, and 6-month outcomes in women with peripartum cardiomyopathy: an ESC EORP registry.** *Eur Heart J* 2020;41:3787–97. — <https://pubmed.ncbi.nlm.nih.gov/32840318/>
   *739-patient global registry: ~6-month mortality 6%, LVEF ≥50% in ~46%. Anchors the pessimistic tail represented by the untreated scenario.*

6. Sliwa K, Mebazaa A, Hilfiker-Kleiner D, et al. **Clinical characteristics of patients from the worldwide registry on peripartum cardiomyopathy (PPCM): EURObservational Research Programme.** *Eur J Heart Fail* 2017;19:1131–41. — <https://pubmed.ncbi.nlm.nih.gov/28271625/>
   *Regional variation in presentation severity and access to care (`Low_SES_access` node).*

7. Kolte D, Khera S, Aronow WS, et al. **Temporal trends in incidence and outcomes of peripartum cardiomyopathy in the United States: a nationwide population-based study.** *J Am Heart Assoc* 2014;3:e001056. — <https://pubmed.ncbi.nlm.nih.gov/24901108/>
   *Incidence and in-hospital complication rates (cardiogenic shock, thromboembolism) used to scale the endpoint cluster.*

8. Isogai T, Kamiya CA. **Worldwide incidence of peripartum cardiomyopathy and overall maternal mortality.** *Int Heart J* 2019;60:503–11. — <https://pubmed.ncbi.nlm.nih.gov/30982820/>
   *Geographic heterogeneity of incidence.*

---

## 2. The prolactin / cathepsin D / 16-kDa prolactin axis (model backbone)

9. Hilfiker-Kleiner D, Kaminski K, Podewski E, et al. **A cathepsin D-cleaved 16 kDa form of prolactin mediates postpartum cardiomyopathy.** *Cell* 2007;128:589–600. — <https://pubmed.ncbi.nlm.nih.gov/17289576/>
   *The single most important mechanistic source. Cardiomyocyte STAT3 knockout → reduced MnSOD → oxidative stress → cathepsin D activation → cleavage of 23-kDa prolactin to the anti-angiogenic 16-kDa form → PPCM, prevented by bromocriptine. Justifies the `LOSS_STAT3 → amp → ROS → CTSD → PRL16` chain and the `CTSD_TOL` threshold.*

10. Hilfiker-Kleiner D, Sliwa K. **Pathophysiology and epidemiology of peripartum cardiomyopathy.** *Nat Rev Cardiol* 2014;11:364–70. — <https://pubmed.ncbi.nlm.nih.gov/24686946/>
    *Synthesis of the prolactin-cleavage and anti-angiogenic hypotheses into the two-hit framework encoded by `SYN_2HIT`.*

11. Clapp C, Thebault S, Jeziorski MC, Martínez de la Escalera G. **Peptide hormone regulation of angiogenesis.** *Physiol Rev* 2009;89:1177–215. — <https://pubmed.ncbi.nlm.nih.gov/19789380/>
    *Vasoinhibin (16K-PRL) biology: endothelial apoptosis, NF-κB activation, anti-angiogenic action. Justifies `KCAP_OUT` and the `PRL16 → NFkB_endo` edge.*

12. Corbacho AM, Martínez de la Escalera G, Clapp C. **Roles of prolactin and related members of the prolactin/growth hormone/placental lactogen family in angiogenesis.** *J Endocrinol* 2002;173:219–38. — <https://pubmed.ncbi.nlm.nih.gov/12010630/>
    *Opposing actions of full-length versus cleaved prolactin — the reason `PRL23` is modelled as pro-angiogenic through PRLR/JAK2/STAT5 while `PRL16` is anti-angiogenic.*

13. Hilfiker-Kleiner D, Struman I, Hoch M, Podewski E, Sliwa K. **STAT3 and myocardial remodeling in pregnancy and heart disease.** *Heart Fail Rev* 2012;17:583–8. — <https://pubmed.ncbi.nlm.nih.gov/21979759/>
    *Cardioprotective STAT3 signalling during the peripartum haemodynamic load.*

14. Ricke-Hoch M, Bultmann I, Stapel B, et al. **Opposing roles of Akt and STAT3 in the protection of the maternal heart from peripartum stress.** *Cardiovasc Res* 2014;101:587–96. — <https://pubmed.ncbi.nlm.nih.gov/24448317/>
    *Supports treating STAT3 capacity as the key modifiable susceptibility parameter rather than a fixed lesion.*

---

## 3. Anti-angiogenic second hit: PGC-1α, VEGF and placental sFlt-1

15. Patten IS, Rana S, Shahul S, et al. **Cardiac angiogenic imbalance leads to peripartum cardiomyopathy.** *Nature* 2012;485:333–8. — <https://pubmed.ncbi.nlm.nih.gov/22596155/>
    *Cardiac PGC-1α knockout mice develop PPCM; disease requires BOTH a cardiac VEGF deficit AND placental sFlt-1; rescued by bromocriptine or VEGF. This is the direct source of the `PGC_LOSS → VEGF → CAP` path, the `KTRAP` sFlt-1 ligand-trap term, and the `Two_hit_model` node.*

16. Levine RJ, Maynard SE, Qian C, et al. **Circulating angiogenic factors and the risk of preeclampsia.** *N Engl J Med* 2004;350:672–83. — <https://pubmed.ncbi.nlm.nih.gov/14764923/>
    *Kinetics of sFlt-1 rise in late gestation and its fall after delivery — the basis for producing `SFLT` only while antepartum with `KOUT_SFLT` clearance thereafter.*

17. Maynard SE, Min JY, Merchan J, et al. **Excess placental soluble fms-like tyrosine kinase 1 (sFlt1) may contribute to endothelial dysfunction, hypertension, and proteinuria in preeclampsia.** *J Clin Invest* 2003;111:649–58. — <https://pubmed.ncbi.nlm.nih.gov/12618519/>
    *Mechanism of sFlt-1-mediated VEGF sequestration and endothelial dysfunction.*

18. Bello N, Rendon ISH, Arany Z. **The relationship between pre-eclampsia and peripartum cardiomyopathy: a systematic review and meta-analysis.** *J Am Coll Cardiol* 2013;62:1715–23. — <https://pubmed.ncbi.nlm.nih.gov/23916925/>
    *Preeclampsia prevalence in PPCM ~22%, roughly four-fold the background rate. Sets the `PREECL` risk weighting and the `PREECL_SFLT` multiplier.*

19. Damp J, Givertz MM, Semigran M, et al. **Relaxin-2 and soluble Flt1 levels in peripartum cardiomyopathy: results of the multicenter IPAC study.** *JACC Heart Fail* 2016;4:380–8. — <https://pubmed.ncbi.nlm.nih.gov/26970830/>
    *Measured sFlt-1 and relaxin-2 in PPCM patients; higher sFlt-1 associates with worse presentation. Supports the relaxin and sFlt-1 nodes.*

---

## 4. miR-146a exosome axis

20. Halkein J, Tabruyn SP, Ricke-Hoch M, et al. **MicroRNA-146a is a therapeutic target and biomarker for peripartum cardiomyopathy.** *J Clin Invest* 2013;123:2143–54. — <https://pubmed.ncbi.nlm.nih.gov/23619365/>
    *16K-PRL → endothelial NF-κB → miR-146a-5p → exosomal transfer to cardiomyocytes → knockdown of Erbb4, Nras and Notch1; serum miR-146a is elevated specifically in PPCM (not in other cardiomyopathies) and falls with bromocriptine. Source of the entire `MIR → SURV` submodel and of the `miR146a_serum` biomarker node.*

21. Hoes MF, Arany Z, Bauersachs J, et al. **Pathophysiology and risk factors of peripartum cardiomyopathy.** *Nat Rev Cardiol* 2022;19:555–65. — <https://pubmed.ncbi.nlm.nih.gov/35296804/>
    *Current mechanistic synthesis including the exosomal miR-146a pathway and genetic susceptibility.*

22. Bajou K, Herkenne S, Thijssen VL, et al. **PAI-1 mediates the antiangiogenic and profibrinolytic effects of 16K prolactin.** *Nat Med* 2014;20:741–7. — <https://pubmed.ncbi.nlm.nih.gov/24929952/>
    *Links 16K-PRL to PAI-1 and the prothrombotic milieu — mechanistic support for coupling PPCM severity to thrombus propensity.*

---

## 5. Genetic substrate

23. Ware JS, Li J, Mazaika E, et al. **Shared genetic predisposition in peripartum and dilated cardiomyopathies.** *N Engl J Med* 2016;374:233–41. — <https://pubmed.ncbi.nlm.nih.gov/26735901/>
    *Truncating variants (predominantly TTN) in ~15% of PPCM, the same architecture as idiopathic DCM, and associated with lower LVEF at 1 year. Justifies `GEN_TTN` acting through both `GEN_INJ` (less sarcomere reserve) and `GEN_REC` (less recovery capacity).*

24. van Spaendonck-Zwarts KY, Posafalvi A, van den Berg MP, et al. **Titin gene mutations are common in families with both peripartum cardiomyopathy and dilated cardiomyopathy.** *Eur Heart J* 2014;35:2165–73. — <https://pubmed.ncbi.nlm.nih.gov/24558114/>
    *Familial co-segregation of PPCM and DCM through TTN — the `DCM_shared_arch` and `Family_history` nodes.*

25. Goli R, Li J, Brandimarto J, et al. **Genetic and phenotypic landscape of peripartum cardiomyopathy.** *Circulation* 2021;143:1852–62. — <https://pubmed.ncbi.nlm.nih.gov/33874732/>
    *Largest genotype–phenotype study to date: truncating-variant carriers have lower LVEF and worse outcomes; also documents ancestry differences.*

26. Sliwa K, Bauersachs J, Arany Z, Spaendonck-Zwarts KY, Hilfiker-Kleiner D. **Peripartum cardiomyopathy: from genetics to management.** *Eur Heart J* 2021;42:3094–102. — <https://pubmed.ncbi.nlm.nih.gov/34322694/>
    *Ties the genetic substrate to the prolactin/anti-angiogenic mechanism and to management — the "susceptible genotype plus peripartum insult" paradigm that `GEN_INJ` and `GEN_REC` encode.*

---

## 6. Natural history, recovery and prognosis (calibration targets)

27. McNamara DM, Elkayam U, Alharethi R, et al. **Clinical outcomes for peripartum cardiomyopathy in North America: results of the IPAC study (Investigations of Pregnancy-Associated Cardiomyopathy).** *J Am Coll Cardiol* 2015;66:905–14. — <https://pubmed.ncbi.nlm.nih.gov/26293760/>
    *Primary calibration target for the standard-therapy arm: 100 women, entry LVEF ~0.35 rising to ~0.51 at 6 months, 72% reaching LVEF ≥0.50 by 12 months; entry LVEF <0.30 with LVEDD ≥6.0 cm predicts non-recovery (reproduced by the wall-stress and dilation terms); lower recovery in Black women.*

28. Irizarry OC, Levine LD, Lewey J, et al. **Comparison of clinical characteristics and outcomes of peripartum cardiomyopathy between African American and non-African American women.** *JAMA Cardiol* 2017;2:1256–60. — <https://pubmed.ncbi.nlm.nih.gov/28793138/>
    *Later presentation, lower recovery and more frequent worsening in Black women — the `African_ancestry` edges to `LVEF_recovery` and `Mortality`.*

29. Goland S, Bitar F, Modi K, et al. **Evaluation of the clinical relevance of baseline left ventricular ejection fraction as a predictor of recovery or persistence of severe dysfunction in women in the United States with peripartum cardiomyopathy.** *J Card Fail* 2011;17:426–30. — <https://pubmed.ncbi.nlm.nih.gov/21549301/>
    *Baseline LVEF as the dominant predictor of recovery — the `SEV` parameterisation of initial conditions.*

30. Blauwet LA, Delgado-Montero A, Ryo K, et al. **Right ventricular function in peripartum cardiomyopathy at presentation is associated with subsequent left ventricular recovery and clinical outcomes.** *Circ Heart Fail* 2016;9:e002756. — <https://pubmed.ncbi.nlm.nih.gov/27152625/>
    *RV involvement in ~30% and its prognostic weight — the `RV_involvement` node.*

31. Elkayam U, Tummala PP, Rao K, et al. **Maternal and fetal outcomes of subsequent pregnancies in women with peripartum cardiomyopathy.** *N Engl J Med* 2001;344:1567–71. — <https://pubmed.ncbi.nlm.nih.gov/11372007/>
    *Relapse risk in subsequent pregnancy, markedly higher when LVEF has not normalised — the `Subsequent_pregnancy` hexagon and its edge back to peripartum wall stress.*

32. Codsi E, Rose CH, Blauwet LA. **Subsequent pregnancy outcomes in patients with peripartum cardiomyopathy.** *Obstet Gynecol* 2018;131:322–7. — <https://pubmed.ncbi.nlm.nih.gov/29324606/>
    *Contemporary subsequent-pregnancy outcome data supporting the counselling/contraception decision node.*

33. Sliwa K, Förster O, Libhaber E, et al. **Peripartum cardiomyopathy: inflammatory markers as predictors of outcome in 100 prospectively studied patients.** *Eur Heart J* 2006;27:441–6. — <https://pubmed.ncbi.nlm.nih.gov/16143707/>
    *TNF-α, IL-6, CRP and Fas/APO-1 as outcome predictors — the inflammatory cluster and the `sFas_Fas → Mortality` prognostic edge.*

---

## 7. Bromocriptine: trials and pharmacology

34. Sliwa K, Blauwet L, Tibazarwa K, et al. **Evaluation of bromocriptine in the treatment of acute severe peripartum cardiomyopathy: a proof-of-concept pilot study.** *Circulation* 2010;121:1465–73. — <https://pubmed.ncbi.nlm.nih.gov/20308616/>
    *First randomised evidence: bromocriptine + standard therapy raised LVEF from 0.27 to 0.58 versus 0.27 to 0.36 with standard therapy alone. Upper bracket for the modelled bromocriptine effect size.*

35. Hilfiker-Kleiner D, Haghikia A, Berliner D, et al. **Bromocriptine for the treatment of peripartum cardiomyopathy: a multicentre randomized study.** *Eur Heart J* 2017;38:2671–9. — <https://pubmed.ncbi.nlm.nih.gov/28934837/>
    *The 1-week versus 8-week comparison (2.5 mg BID × 1 week versus 2.5 mg BID × 2 weeks then 2.5 mg daily × 6 weeks) on top of standard therapy: full recovery in 52% versus 68%, very low MACE in both arms. Directly defines the two bromocriptine regimens in `ppcm_brc()`. **Note the model's known overshoot:** the trial's mean-LVEF difference between arms was small (0.49 versus 0.51) while the recovery-rate difference was 16 points; the model reproduces the direction and the mechanism (vulnerable-window coverage plus irreversible-scar lock-in) but exaggerates the mean-LVEF gap (47.6% versus 56.5%).*

36. Koenig T, Bauersachs J, Hilfiker-Kleiner D. **Bromocriptine for the treatment of peripartum cardiomyopathy.** *Card Fail Rev* 2018;4:46–9. — <https://pubmed.ncbi.nlm.nih.gov/29892479/>
    *Practical dosing, duration rationale and the requirement for concomitant anticoagulation.*

37. Tremblay-Gravel M, Marquis-Gravel G, Avram R, et al. **The effect of bromocriptine on left ventricular functional recovery in peripartum cardiomyopathy: insights from the BRO-HF retrospective cohort study.** *ESC Heart Fail* 2019;6:27–36. — <https://pubmed.ncbi.nlm.nih.gov/30565877/>
    *A cohort that did NOT find a significant independent bromocriptine benefit — the reason the model's bromocriptine effect is presented as a mechanistic hypothesis rather than settled quantitative fact.*

38. Trongtorsak A, Kittipibul V, Cheungpasitporn W, et al. **Effects of bromocriptine in peripartum cardiomyopathy: a systematic review and meta-analysis.** *Heart Fail Rev* 2022;27:533–43. — <https://pubmed.ncbi.nlm.nih.gov/34129137/>
    *Pooled estimate of LVEF improvement and recovery rate with bromocriptine.*

39. Del Dotto P, Bonuccelli U. **Clinical pharmacokinetics of cabergoline.** *Clin Pharmacokinet* 2003;42:633–45. — <https://pubmed.ncbi.nlm.nih.gov/12844325/>
    *Cited as the closest well-characterised ergot dopamine-agonist PK profile: extensive first-pass metabolism, low oral bioavailability, CYP3A4 clearance, and prolactin suppression that outlasts plasma exposure. Basis for `KA_BRC`, `KOUT_BRC`, `KTR_BRC` and for modelling effect through a smooth `EMAX/EC50` term rather than plasma concentration alone. (Bromocriptine's own labelled PK — oral bioavailability ~6%, CYP3A4 metabolism, biliary excretion — is taken from its product information rather than a primary paper.)*

40. **Bromocriptine thrombotic risk.** No single primary study is cited here because the evidence is regulatory and case-based: bromocriptine was withdrawn from the indication of routine postpartum lactation suppression after reports of stroke, myocardial infarction and seizures, which is precisely why the ESC position statement mandates concomitant anticoagulation whenever it is used for PPCM. See Bauersachs 2019 (ref 2), Koenig 2018 (ref 36) and Sliwa 2010 (ref 34) for the discussion in the PPCM context.
    *Justifies `PROTHR_BRC` and the modelled requirement for concomitant anticoagulation (`EAC_THR`), and is the reason scenarios 8 and 9 are shipped as a matched pair.*

---

## 8. Heart-failure pharmacotherapy in pregnancy and lactation

41. Regitz-Zagrosek V, Roos-Hesselink JW, Bauersachs J, et al. **2018 ESC Guidelines for the management of cardiovascular diseases during pregnancy.** *Eur Heart J* 2018;39:3165–241. — <https://pubmed.ncbi.nlm.nih.gov/30165544/>
    *Source of the antepartum ACEi/ARB fetotoxicity contraindication that the model enforces as a hard gate (`ANTEPARTUM_RASI_BLOCK`), of hydralazine + nitrate as the antepartum-safe vasodilator substitute, and of the breastfeeding-compatibility filter.*

42. McDonagh TA, Metra M, Adamo M, et al. **2021 ESC Guidelines for the diagnosis and treatment of acute and chronic heart failure.** *Eur Heart J* 2021;42:3599–726. — <https://pubmed.ncbi.nlm.nih.gov/34447992/>
    *The four-pillar HFrEF stack (ACEi/ARNI, beta-blocker, MRA, SGLT2 inhibitor) whose effects are parameterised in the HF pharmacotherapy cluster.*

43. Cooper WO, Hernandez-Diaz S, Arbogast PG, et al. **Major congenital malformations after first-trimester exposure to ACE inhibitors.** *N Engl J Med* 2006;354:2443–51. — <https://pubmed.ncbi.nlm.nih.gov/16760444/>
    *Fetal risk of ACE inhibition, underpinning the antepartum gate.*

44. Bateman BT, Patorno E, Desai RJ, et al. **Angiotensin-converting enzyme inhibitors and the risk of congenital malformations.** *Obstet Gynecol* 2017;129:174–84. — <https://pubmed.ncbi.nlm.nih.gov/27984527/>
    *Further characterisation of ACEi fetal risk.*

45. Halpern DG, Weinberg CR, Pinnelas R, et al. **Use of medication for cardiovascular disease during pregnancy: JACC state-of-the-art review.** *J Am Coll Cardiol* 2019;73:457–76. — <https://pubmed.ncbi.nlm.nih.gov/30717713/>
    *Drug-by-drug pregnancy and lactation compatibility — the `Breastfeeding_compat` decision diamond.*

46. Biteker M, Duran NE, Kaya H, et al. **Effect of levosimendan and predictors of recovery in patients with peripartum cardiomyopathy: a randomized clinical trial.** *Clin Res Cardiol* 2011;100:571–7. — <https://pubmed.ncbi.nlm.nih.gov/21212959/>
    *Neutral trial. Deliberately modelled as having NO durable PD effect — levosimendan appears in the map for completeness only.*

47. Sliwa K, Skudicky D, Candy G, et al. **The addition of pentoxifylline to conventional therapy improves outcome in patients with peripartum cardiomyopathy.** *Eur J Heart Fail* 2002;4:305–9. — <https://pubmed.ncbi.nlm.nih.gov/12034156/>
    *Small, unreplicated positive study. Included in the map, deliberately given no PD effect in the ODE model.*

48. Bozkurt B, Villaneuva FS, Holubkov R, et al. **Intravenous immune globulin in the therapy of peripartum cardiomyopathy.** *J Am Coll Cardiol* 1999;34:177–80. — <https://pubmed.ncbi.nlm.nih.gov/10399005/>
    *IVIG pilot data; same treatment in the model as pentoxifylline (mapped, not efficacious).*

49. Duncker D, Haghikia A, König T, et al. **Risk for ventricular fibrillation in peripartum cardiomyopathy with severely reduced left ventricular function — value of the wearable cardioverter/defibrillator.** *Eur J Heart Fail* 2014;16:1331–6. — <https://pubmed.ncbi.nlm.nih.gov/25371320/>
    *Basis for the wearable-defibrillator bridge and the deferral of permanent ICD pending recovery.*

50. Sieweke JT, Pfeffer TJ, Berliner D, et al. **Cardiogenic shock complicating peripartum cardiomyopathy: importance of early left ventricular unloading and bromocriptine therapy.** *Eur Heart J Acute Cardiovasc Care* 2020;9:173–82. — <https://pubmed.ncbi.nlm.nih.gov/30764627/>
    *Management of the cardiogenic-shock phenotype (IABP/ECMO/LVAD) and early bromocriptine in that setting.*

---

## 9. Peripartum cardiovascular physiology (load driver)

51. Sanghavi M, Rutherford JD. **Cardiovascular physiology of pregnancy.** *Circulation* 2014;130:1003–8. — <https://pubmed.ncbi.nlm.nih.gov/25223771/>
    *Plasma-volume expansion of 40–50%, cardiac-output rise of 30–50%, fall in systemic vascular resistance, and the postpartum autotransfusion preload spike. Quantitative basis for `LOAD_PREG`, `LOAD_PP` and `TAU_PP`.*

52. Melchiorre K, Sharma R, Thilaganathan B. **Cardiac structure and function in normal pregnancy.** *Curr Opin Obstet Gynecol* 2012;24:413–21. — <https://pubmed.ncbi.nlm.nih.gov/23000697/>
    *Physiological eccentric LV hypertrophy and the normal ranges used for `LVEDV_N` and `EF_MAX`.*

53. Grattan DR. **60 years of neuroendocrinology: the hypothalamo-prolactin axis.** *J Endocrinol* 2015;226:T101–22. — <https://pubmed.ncbi.nlm.nih.gov/26101377/>
    *Dopaminergic (D2) tonic inhibition of lactotrophs, gestational prolactin rise to ~200 ng/mL and the suckling-induced surge. Basis for `DRIVE_PREG`, `DRIVE_LACT` and the bromocriptine mechanism of action.*

54. Neville MC, Morton J. **Physiology and endocrine changes underlying human lactogenesis II.** *J Nutr* 2001;131:3005S–8S. — <https://pubmed.ncbi.nlm.nih.gov/11694638/>
    *Supply-and-demand autoregulation of milk production, and the practical irreversibility of established involution. Justifies the autocatalytic `LACT` compartment with the absorbing `LACT_MIN` threshold.*

---

## 10. Biomarkers and imaging endpoints

55. Forster O, Hilfiker-Kleiner D, Ansari AA, et al. **Reversal of IFN-γ, oxLDL and prolactin serum levels correlate with clinical improvement in patients with peripartum cardiomyopathy.** *Eur J Heart Fail* 2008;10:861–8. — <https://pubmed.ncbi.nlm.nih.gov/18768352/>
    *Oxidised-LDL and prolactin fall in parallel with recovery — supports using the oxidative-stress and prolactin states as tracked biomarkers.*

56. Hoevelmann J, Engel ME, Muller E, et al. **A global perspective on the management and outcomes of peripartum cardiomyopathy: a systematic review and meta-analysis.** *Eur J Heart Fail* 2022;24:1719–36. — <https://pubmed.ncbi.nlm.nih.gov/35867781/>
    *Pooled recovery and mortality rates across regions; NT-proBNP as an outcome marker.*

57. Arany Z. **Understanding peripartum cardiomyopathy.** *Annu Rev Med* 2018;69:165–76. — <https://pubmed.ncbi.nlm.nih.gov/29029584/>
    *Concise mechanistic review used to cross-check the map's edge directions.*

58. Haghikia A, Podewski E, Libhaber E, et al. **Phenotyping and outcome on contemporary management in a German cohort of patients with peripartum cardiomyopathy.** *Basic Res Cardiol* 2013;108:366. — <https://pubmed.ncbi.nlm.nih.gov/23812247/>
    *German cohort presenting LVEF (~27%), NT-proBNP levels and LVEDD distribution — used to set presenting-state initial conditions.*

59. Schelbert EB, Elkayam U, Cooper LT, et al. **Myocardial damage detected by late gadolinium enhancement cardiac magnetic resonance is uncommon in peripartum cardiomyopathy.** *J Am Heart Assoc* 2017;6:e005472. — <https://pubmed.ncbi.nlm.nih.gov/28468786/>
    *Important counterweight: extensive replacement fibrosis is NOT typical of PPCM. The model's `SCAR` fraction is therefore kept modest in treated arms (0.13–0.25) and grows large only in the untreated, progressively deteriorating arm.*

60. Ersbøll AS, Bojer AS, Hauge MG, et al. **Long-term cardiac function after peripartum cardiomyopathy and preeclampsia: a Danish nationwide, clinical follow-up study using maximal exercise testing and cardiac magnetic resonance imaging.** *J Am Heart Assoc* 2018;7:e008991. — <https://pubmed.ncbi.nlm.nih.gov/30371239/>
    *Long-term functional follow-up including subclinical residual dysfunction — supports the model's incomplete-recovery plateau rather than full normalisation in every arm.*

---

## 11. QSP methodology

61. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus machine learning, a fight worth fighting for the biological community?** *Biol Lett* 2018;14:20170660. — <https://pubmed.ncbi.nlm.nih.gov/29769297/>
    *Rationale for mechanistic, interpretable models of the kind built here.*

62. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019;8:883–93. — <https://pubmed.ncbi.nlm.nih.gov/31674724/>
    *The mrgsolve implementation idiom used in `ppcm_mrgsolve_model.R`.*

63. Bai JPF, Schmidt BJ, Gadkar KG, et al. **FDA-industry scientific exchange on assessing quantitative systems pharmacology models in clinical drug development.** *J Pharmacokinet Pharmacodyn* 2021;48:453–9. — <https://pubmed.ncbi.nlm.nih.gov/33835309/>
    *Expectations for QSP model credibility assessment — the reason this model ships with an explicit verification table and a stated overshoot against its own anchor.*

---

## Model limitations stated explicitly

- **The bromocriptine duration effect is exaggerated.** The model reproduces the *direction* of the German trial's 1-week versus 8-week result and offers a mechanism for it (coverage of the high-oxidative-stress vulnerable window, with early damage locked in as irreversible scar), but the mean-LVEF gap it produces (47.6% versus 56.5% at 6 months) is larger than the trial's (0.49 versus 0.51). Treat the duration comparison as a hypothesis, not a prediction.
- **Bromocriptine efficacy itself is contested.** Ref 37 (BRO-HF) did not find an independent benefit. The model encodes the mechanistic hypothesis of refs 9, 15, 20, 34 and 35.
- **Drug exposures are in dose-proportional arbitrary units**, not validated plasma PK. Only bromocriptine has a structured (two-compartment plus depot) representation; the remaining agents are single exposure compartments.
- **Fibrosis is probably over-weighted as a recovery brake.** Ref 59 shows that late gadolinium enhancement is uncommon in PPCM; the `FIB`/`SCAR` split is the model's attempt to respect that while still producing path-dependent recovery, but the balance between reversible functional impairment and true scar is not identified from data.
- **No mortality/competing-risk model.** `MACE_hazard_yr` is an illustrative algebraic surrogate, not a fitted survival model, and no patient leaves the simulation.
- **Single-patient deterministic simulation.** There is no inter-individual variability block (`$OMEGA`/`$SIGMA`), so scenario outputs are typical-patient trajectories rather than recovery *rates*; comparing the model's LVEF values to trial recovery percentages requires that caveat.
