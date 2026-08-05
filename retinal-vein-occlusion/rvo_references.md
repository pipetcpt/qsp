# 망막정맥폐쇄(RVO) QSP 모델 — 참고문헌
# Retinal Vein Occlusion QSP Model — Annotated References

모든 PMID는 2026-08-05에 PubMed E-utilities로 조회하여 제목·저자·저널·연도를 확인했습니다.
Every PMID below was resolved and verified against the PubMed E-utilities API on
2026-08-05; the title, first author, journal and year shown are the ones PubMed
returned. Each entry states **what the model actually took from it** — a
parameter, a structural assumption, or a calibration target.

---

## 1. 역학과 자연사 (Epidemiology and natural history)

Sets the phenotype prevalences, the baseline central subfield thickness and BCVA
at presentation, and the untreated trajectories the model must reproduce.

1. McIntosh RL, et al. **Natural history of central retinal vein occlusion: an evidence-based systematic review.** *Ophthalmology* 2010;117:1113-23.
   [PMID 20430446](https://pubmed.ncbi.nlm.nih.gov/20430446/)
   → the untreated CRVO trajectory used to check `NAT_NONISCH` / `NAT_ISCH`.
2. Rogers SL, et al. **Natural history of branch retinal vein occlusion: an evidence-based systematic review.** *Ophthalmology* 2010;117:1094-101.
   [PMID 20430447](https://pubmed.ncbi.nlm.nih.gov/20430447/)
   → the BRVO natural history; the reason `OCC_RES` for BRVO is 0.32 and not 0.55.
3. Cugati S, et al. **Ten-year incidence of retinal vein occlusion in an older population: the Blue Mountains Eye Study.** *Arch Ophthalmol* 2006;124:726-32.
   [PMID 16682596](https://pubmed.ncbi.nlm.nih.gov/16682596/)
4. Klein R, et al. **The 15-year cumulative incidence of retinal vein occlusion: the Beaver Dam Eye Study.** *Arch Ophthalmol* 2008;126:513-8.
   [PMID 18413521](https://pubmed.ncbi.nlm.nih.gov/18413521/)
   → fellow-eye risk over 5-15 years, quoted in the decision layer of the map.
5. Hayreh SS, et al. **Incidence of various types of retinal vein occlusion and their recurrence and demographic characteristics.** *Am J Ophthalmol* 1994;117:429-41.
   [PMID 8154523](https://pubmed.ncbi.nlm.nih.gov/8154523/)
6. Central Vein Occlusion Study Group. **Baseline and early natural history report.** *Arch Ophthalmol* 1993;111:1087-95.
   [PMID 7688950](https://pubmed.ncbi.nlm.nih.gov/7688950/)
   → baseline BCVA distribution; the ischaemic / non-ischaemic split.
7. Central Vein Occlusion Study Group. **Natural history and clinical management of central retinal vein occlusion.** *Arch Ophthalmol* 1997;115:486-91.
   [PMID 9109757](https://pubmed.ncbi.nlm.nih.gov/9109757/)
   → the conversion rate from non-ischaemic to ischaemic CRVO; `K_NP`.
8. Branch Vein Occlusion Study Group. **Argon laser photocoagulation for macular edema in branch vein occlusion.** *Am J Ophthalmol* 1984;98:271-82.
   [PMID 6383055](https://pubmed.ncbi.nlm.nih.gov/6383055/)
   → grid laser as the pre-anti-VEGF comparator; the untreated BRVO arm.
9. Cugati S, et al. **Retinal vein occlusion and vascular mortality: pooled data analysis of 2 population-based cohorts.** *Ophthalmology* 2007;114:520-4.
   [PMID 17141315](https://pubmed.ncbi.nlm.nih.gov/17141315/)
   → why the model treats RVO as a systemic vascular event, not only an eye event.

## 2. 압력항의 근거 — 망막 모세혈관 정수압 (The evidence for the pressure arm)

These are the papers that make the pressure arm a measurement rather than a
metaphor. They are the reason `PA = 40`, `RV0 = 0.25` and `PC ≈ 21.6 mmHg` in a
normal eye, and the reason a venous occlusion drives Pc toward the arterial
inlet pressure.

10. Glucksberg MR, Dunn R. **Direct measurement of retinal microvascular pressures in the live, anesthetized cat.** *Microvasc Res* 1993;45:158-65.
    [PMID 8361399](https://pubmed.ncbi.nlm.nih.gov/8361399/)
    → **the source of the pressure-arm numbers.** Servo-null micropipette
    measurements across the retinal microvasculature; sets the arteriolar inlet
    and capillary pressures the model uses.
11. Attariwala R, Giebs CP, Glucksberg MR. **The influence of elevated intraocular pressure on vascular pressures in the cat retina.** *Invest Ophthalmol Vis Sci* 1994;35:1019-25.
    [PMID 8125712](https://pubmed.ncbi.nlm.nih.gov/8125712/)
    → justifies `P_v = IOP + 2 mmHg` and the IOP → Pv → Pc feedback loop.
12. Attariwala R, et al. **The effect of acute experimental retinal vein occlusion on cat retinal vein pressures.** *Invest Ophthalmol Vis Sci* 1997;38:2742-9.
    [PMID 9418726](https://pubmed.ncbi.nlm.nih.gov/9418726/)
    → **direct measurement of what an experimental vein occlusion does to venous
    pressure.** The basis for the `OCC` multiplier magnitudes.
13. Levick JR, Michel CC. **Microvascular fluid exchange and the revised Starling principle.** *Cardiovasc Res* 2010;87:198-210.
    [PMID 20200043](https://pubmed.ncbi.nlm.nih.gov/20200043/)
    → the Starling form used in the model, including the role of the reflection
    coefficient σ and why interstitial oncotic pressure is not free to rise
    without limit.
14. Michel CC, Woodcock TE, Curry FE. **Understanding and extending the Starling principle.** *Acta Anaesthesiol Scand* 2020;64:1032-7.
    [PMID 32270491](https://pubmed.ncbi.nlm.nih.gov/32270491/)
    → the glycocalyx / sub-glycocalyx interpretation of σ that the model encodes
    as `SIGMA = SIGMAX * TJ^0.55`.
15. Ruffini A, et al. **A mathematical model of interstitial fluid flow and retinal tissue deformation in macular edema.** *Invest Ophthalmol Vis Sci* 2024;65:33.
    [PMID 39254963](https://pubmed.ncbi.nlm.nih.gov/39254963/)
    → an independent poroelastic treatment of retinal fluid; the closest published
    precedent for treating macular oedema as a filtration / clearance balance
    rather than as a permeability index.
16. Hayreh SS. **Prevalent misconceptions about acute retinal vascular occlusive disorders.** *Prog Retin Eye Res* 2005;24:493-519.
    [PMID 15845346](https://pubmed.ncbi.nlm.nih.gov/15845346/)
    → the argument that the haemodynamic state, not the mediator profile, defines
    the RVO phenotype. This paper is the closest prior statement of the model's
    central claim.
17. Williamson TH. **Central retinal vein occlusion: what's the story?** *Br J Ophthalmol* 1997;81:698-704.
    [PMID 9349161](https://pubmed.ncbi.nlm.nih.gov/9349161/)
    → the lamina cribrosa anatomy and why CRVO has so few collateral options
    (`COLL_MAX` 0.45-0.70 for CRVO versus 2.2 for BRVO).

## 3. 혈전과 재개통 · 우회로 시계 (Thrombus, recanalisation, collaterals)

18. Rehak M, Wiedemann P. **Retinal vein thrombosis: pathogenesis and management.** *J Thromb Haemost* 2010;8:1886-94.
    [PMID 20492457](https://pubmed.ncbi.nlm.nih.gov/20492457/)
    → the Virchow decomposition in cluster 3 and the `K_LYS` order of magnitude.
19. Tsuboi K, et al. **Collateral vessels in branch retinal vein occlusion: anatomic and functional analyses by OCT angiography.** *Ophthalmol Retina* 2019;3:767-76.
    [PMID 31167729](https://pubmed.ncbi.nlm.nih.gov/31167729/)
    → **the collateral clock made visible.** OCT-A quantification of collateral
    formation in BRVO; the reason the collateral time constant is months
    (`K_COLL = 0.0035 /d`) and the reason BRVO crosses below Pc* and CRVO does not.
20. Tsuboi K, et al. **Longitudinal vasculature changes in branch retinal vein occlusion with projection-resolved optical coherence tomography angiography.** *Graefes Arch Clin Exp Ophthalmol* 2019;257:1875-83.
    [PMID 31165932](https://pubmed.ncbi.nlm.nih.gov/31165932/)
21. Tsuboi K, et al. **Chorioretinal shunt vessel in eyes with central retinal vein occlusion after radial optic neurotomy.** *Ophthalmology* 2018;125:1332-3.
    [PMID 30143093](https://pubmed.ncbi.nlm.nih.gov/30143093/)
    → the surgical attempt to add conductance directly to the venous side.
22. Opremcak EM, Bruce RA, Lomeo MD, et al. **Radial optic neurotomy for central retinal vein occlusion: a retrospective pilot study of 11 consecutive cases.** *Retina* 2001;21:408-15.
    [PMID 11642369](https://pubmed.ncbi.nlm.nih.gov/11642369/)
    → one of the very few interventions that targets `R_v` rather than `L_p`;
    included in the map precisely because that arrow is so rare.
23. Glacet-Bernard A, et al. **Hemodilution therapy using automated erythrocytapheresis in central retinal vein occlusion.** *Graefes Arch Clin Exp Ophthalmol* 2011;249:505-12.
    [PMID 20953877](https://pubmed.ncbi.nlm.nih.gov/20953877/)
    → viscosity as a lever on `R_v`; the other left-half arrow.

## 4. 혈액-망막 장벽과 투과성 (Barrier biology: the permeability arm)

24. Klaassen I, Van Noorden CJ, Schlingemann RO. **Molecular basis of the inner blood-retinal barrier and its breakdown in diabetic macular edema and other pathological conditions.** *Prog Retin Eye Res* 2013;34:19-48.
    [PMID 23416119](https://pubmed.ncbi.nlm.nih.gov/23416119/)
    → the junction complex in cluster 11 and the paracellular / transcellular split.
25. Daruich A, et al. **Mechanisms of macular edema: beyond the surface.** *Prog Retin Eye Res* 2018;63:20-68.
    [PMID 29126927](https://pubmed.ncbi.nlm.nih.gov/29126927/)
    → the single best synthesis of macular oedema as a balance of influx and
    clearance; the structural template for clusters 11-13.
26. Antonetti DA, et al. **Vascular permeability in experimental diabetes is associated with reduced endothelial occludin content.** *Diabetes* 1998;47:1953-9.
    [PMID 9836530](https://pubmed.ncbi.nlm.nih.gov/9836530/)
    → occludin loss as the measurable correlate of `TJ`.
27. Murakami T, et al. **Protein kinase Cβ phosphorylates occludin regulating tight junction trafficking in vascular endothelial growth factor-induced permeability in vivo.** *Diabetes* 2012;61:1573-83.
    [PMID 22438576](https://pubmed.ncbi.nlm.nih.gov/22438576/)
    → the VEGF → PKCβ → occludin Ser490 → internalisation chain; the mechanistic
    content of `KDIS_TJ * DRIVE * TJ`.
28. Bringmann A, et al. **Müller cells in the healthy and diseased retina.** *Prog Retin Eye Res* 2006;25:397-424.
    [PMID 16839797](https://pubmed.ncbi.nlm.nih.gov/16839797/)
    → AQP4 / Kir4.1 water handling; the saturable clearance term `J_out`.
29. Jo AO, et al. **TRPV4 and AQP4 channels synergistically regulate cell volume and calcium homeostasis in retinal Müller glia.** *J Neurosci* 2015;35:13525-37.
    [PMID 26424896](https://pubmed.ncbi.nlm.nih.gov/26424896/)
    → why clearance saturates rather than scaling linearly with oedema.
30. Wang T, et al. **Anti-VEGF therapy prevents Müller intracellular edema by decreasing VEGF-A in diabetic retinopathy.** *Eye Vis (Lond)* 2021;8:13.
    [PMID 33865457](https://pubmed.ncbi.nlm.nih.gov/33865457/)
31. Suzuma I, et al. **Cyclic stretch and hypertension induce retinal expression of vascular endothelial growth factor and vascular endothelial growth factor receptor-2.** *Diabetes* 2001;50:444-54.
    [PMID 11272159](https://pubmed.ncbi.nlm.nih.gov/11272159/)
    → **the justification for the mechanotransduction term** `W_PRESS * (Pc - 28)/12`
    in `DRIVE`: pressure itself induces VEGF and opens junctions, so the pressure
    arm partially drives the permeability arm.
32. Suzuma I, et al. **Stretch-induced retinal vascular endothelial growth factor expression is mediated by phosphatidylinositol 3-kinase and protein kinase C.** *J Biol Chem* 2002;277:1047-57.
    [PMID 11694503](https://pubmed.ncbi.nlm.nih.gov/11694503/)

## 5. 매개체: VEGF, IL-6, Ang-2 (Mediators)

33. Aiello LP, et al. **Vascular endothelial growth factor in ocular fluid of patients with diabetic retinopathy and other retinal disorders.** *N Engl J Med* 1994;331:1480-7.
    [PMID 7526212](https://pubmed.ncbi.nlm.nih.gov/7526212/)
    → the original human intraocular VEGF measurements; anchors `VTONE`.
34. Noma H, et al. **Aqueous humour levels of cytokines are correlated to vitreous levels and severity of macular oedema in branch retinal vein occlusion.** *Eye (Lond)* 2008;22:42-8.
    [PMID 16826241](https://pubmed.ncbi.nlm.nih.gov/16826241/)
    → **the calibration source for `VTONE`, `IL6` and their ratio to CST.**
    Justifies aqueous sampling as a proxy for the vitreous compartment.
35. Noma H, et al. **Vitreous levels of interleukin-6 and vascular endothelial growth factor in macular edema with central retinal vein occlusion.** *Ophthalmology* 2009;116:87-93.
    [PMID 19118700](https://pubmed.ncbi.nlm.nih.gov/19118700/)
    → normal ~50-100 pg/mL VEGF versus 1500-2500 pg/mL in ischaemic CRVO; the
    13-fold induction `1 + E_HIF_V * HIF`. IL-6 5 → ~85 pg/mL sets `E_HIF_I = 16`.
36. Noma H, et al. **Role of inflammation in previously untreated macular edema with branch retinal vein occlusion.** *BMC Ophthalmol* 2014;14:67.
    [PMID 24884703](https://pubmed.ncbi.nlm.nih.gov/24884703/)
    → **the reason IL-6 is a separate state and not folded into VEGF.** IL-6
    correlates with oedema independently of VEGF, so anti-VEGF leaves it standing.
37. Noma H, et al. **Vascular endothelial growth factor and interleukin-6 are correlated with serous retinal detachment in central retinal vein occlusion.** *Curr Eye Res* 2012;37:62-7.
    [PMID 21978265](https://pubmed.ncbi.nlm.nih.gov/21978265/)
    → the separate `SRF` compartment and its distinct mediator dependence.
38. Saharinen P, Eklund L, Alitalo K. **Therapeutic targeting of the angiopoietin-TIE pathway.** *Nat Rev Drug Discov* 2017;16:635-61.
    [PMID 28529319](https://pubmed.ncbi.nlm.nih.gov/28529319/)
    → Ang-2 as a Tie2 antagonist; the `TIE2_BOOST` term that makes Ang-2 blockade
    raise the junction *repair rate constant* rather than simply lower `DRIVE`.
39. Campochiaro PA, Peters KG. **Targeting Tie2 for treatment of diabetic retinopathy and diabetic macular edema.** *Curr Diab Rep* 2016;16:126.
    [PMID 27778249](https://pubmed.ncbi.nlm.nih.gov/27778249/)
40. Shen J, et al. **Targeting VE-PTP activates TIE2 and stabilizes the ocular vasculature.** *J Clin Invest* 2014;124:4564-76.
    [PMID 25180601](https://pubmed.ncbi.nlm.nih.gov/25180601/)
    → the VE-PTP node in cluster 10.

## 6. 안구 약동학 (Ocular pharmacokinetics)

The `TVIT` half-lives, the penetration factor `FPEN`, and the systemic
concentrations all come from this block.

41. Krohne TU, Eter N, Holz FG, Meyer CH. **Intraocular pharmacokinetics of bevacizumab after a single intravitreal injection in humans.** *Am J Ophthalmol* 2008;146:508-12.
    [PMID 18635152](https://pubmed.ncbi.nlm.nih.gov/18635152/)
    → aqueous half-life 9.82 d for bevacizumab.
42. Krohne TU, Liu Z, Holz FG, Meyer CH. **Intraocular pharmacokinetics of ranibizumab following a single intravitreal injection in humans.** *Am J Ophthalmol* 2012;154:682-6.
    [PMID 22818800](https://pubmed.ncbi.nlm.nih.gov/22818800/)
    → aqueous half-life 7.19 d for ranibizumab; the model's `TVIT` for RBZ.
43. Bakri SJ, et al. **Pharmacokinetics of intravitreal ranibizumab (Lucentis).** *Ophthalmology* 2007;114:2179-82.
    [PMID 18054637](https://pubmed.ncbi.nlm.nih.gov/18054637/)
    → **the source of the retina-to-vitreous partition** that the model condenses
    into `FPEN ≈ 0.08`.
44. Bakri SJ, et al. **Pharmacokinetics of intravitreal bevacizumab (Avastin).** *Ophthalmology* 2007;114:855-9.
    [PMID 17467524](https://pubmed.ncbi.nlm.nih.gov/17467524/)
45. Caruso A, et al. **Ocular half-life of intravitreal biologics in humans and other species: meta-analysis and model-based prediction.** *Mol Pharm* 2020;17:695-709.
    [PMID 31876425](https://pubmed.ncbi.nlm.nih.gov/31876425/)
    → the size-dependence of ocular half-life; the basis for brolucizumab 5.1 d
    and faricimab 7.5 d, which are not directly measured in humans.
46. Hutton-Smith LA, et al. **A mechanistic model of the intravitreal pharmacokinetics of large molecules and the pharmacodynamic suppression of ocular VEGF levels by ranibizumab in patients with neovascular age-related macular degeneration.** *Mol Pharm* 2016;13:2941-50.
    [PMID 26726925](https://pubmed.ncbi.nlm.nih.gov/26726925/)
    → **the direct QSP precedent for this model's PD layer**, and the paper that
    frames the discrepancy the `ALPHA` parameter absorbs.
47. Hutton-Smith LA, et al. **Ocular pharmacokinetics of therapeutic antibodies given by intravitreal injection: estimation of retinal permeabilities using a 3-compartment semi-mechanistic model.** *Mol Pharm* 2017;14:2690-6.
    [PMID 28631484](https://pubmed.ncbi.nlm.nih.gov/28631484/)
    → the retinal permeability estimates behind `FPEN`.
48. Avery RL, et al. **Systemic pharmacokinetics and pharmacodynamics of intravitreal aflibercept, bevacizumab, and ranibizumab.** *Retina* 2017;37:1847-58.
    [PMID 28106709](https://pubmed.ncbi.nlm.nih.gov/28106709/)
    → **the calibration target for the systemic arm.** Plasma aflibercept ≈ 0.03
    µg/mL with near-complete plasma free-VEGF suppression; ranibizumab plasma
    concentrations an order of magnitude lower with minimal suppression. The
    model reproduces both from a single flip-flop structure.
49. Chang-Lin JE, et al. **Pharmacokinetics and pharmacodynamics of a sustained-release dexamethasone intravitreal implant.** *Invest Ophthalmol Vis Sci* 2011;52:80-6.
    [PMID 20702826](https://pubmed.ncbi.nlm.nih.gov/20702826/)
    → `K_REL_DEX`, `KEL_DEX` and the ~900 ng/mL vitreous plateau.
50. Chang-Lin JE, et al. **Pharmacokinetics of a sustained-release dexamethasone intravitreal implant in vitrectomized and nonvitrectomized eyes.** *Invest Ophthalmol Vis Sci* 2011;52:4605-9.
    [PMID 21421864](https://pubmed.ncbi.nlm.nih.gov/21421864/)
    → the vitrectomy node in cluster 18.

## 7. 결합 친화도와 억제 지속시간 (Binding affinity and suppression duration)

51. Papadopoulos N, et al. **Binding and neutralization of vascular endothelial growth factor (VEGF) and related ligands by VEGF Trap, ranibizumab and bevacizumab.** *Angiogenesis* 2012;15:171-85.
    [PMID 22302382](https://pubmed.ncbi.nlm.nih.gov/22302382/)
    → **the `KDV` values used verbatim**: aflibercept 0.49 pM, ranibizumab 46 pM,
    bevacizumab 58 pM for VEGF-A165, and aflibercept's PlGF affinity.
52. Stewart MW, Rosenfeld PJ. **Predicted biological activity of intravitreal VEGF Trap.** *Br J Ophthalmol* 2008;92:667-8.
    [PMID 18356264](https://pubmed.ncbi.nlm.nih.gov/18356264/)
    → **the calculation this model corrects.** The affinity-and-half-life argument
    predicts months of VEGF suppression from one aflibercept injection. The model
    reproduces that arithmetic (96 days) and then shows that the retreatment
    interval is set by the pressure arm rather than by loss of suppression.
53. Stewart MW. **Pharmacokinetics, pharmacodynamics and pre-clinical characteristics of ophthalmic drugs that bind VEGF.** *Expert Rev Clin Pharmacol* 2014;7:167-80.
    [PMID 24483136](https://pubmed.ncbi.nlm.nih.gov/24483136/)
54. Stewart MW. **Pharmacokinetic rationale for dosing every 2 weeks versus 4 weeks with intravitreal ranibizumab, bevacizumab, and aflibercept.** *Retina* 2012;32:434-45.
    [PMID 22374154](https://pubmed.ncbi.nlm.nih.gov/22374154/)
55. Blanot M, et al. **Aflibercept off-target effects in diabetic macular edema: an in silico modeling approach.** *Int J Mol Sci* 2024;25:3621.
    [PMID 38612432](https://pubmed.ncbi.nlm.nih.gov/38612432/)
    → an independent in-silico treatment of aflibercept's multi-ligand binding.

## 8. 임상시험 — 항VEGF (Anti-VEGF trials: the calibration targets)

56. Brown DM, et al. **Ranibizumab for macular edema following central retinal vein occlusion: six-month primary end point results of a phase III study (CRUISE).** *Ophthalmology* 2010;117:1124-33.
    [PMID 20381871](https://pubmed.ncbi.nlm.nih.gov/20381871/)
    → **target: +14.9 letters and −452 µm CST at month 6; sham +0.8 letters.**
57. Campochiaro PA, et al. **Ranibizumab for macular edema following branch retinal vein occlusion: six-month primary end point results of a phase III study (BRAVO).** *Ophthalmology* 2010;117:1102-12.
    [PMID 20398941](https://pubmed.ncbi.nlm.nih.gov/20398941/)
    → **target: +18.3 letters at month 6 in BRVO.**
58. Campochiaro PA, et al. **Sustained benefits from ranibizumab for macular edema following central retinal vein occlusion: twelve-month outcomes of a phase III study.** *Ophthalmology* 2011;118:2041-9.
    [PMID 21715011](https://pubmed.ncbi.nlm.nih.gov/21715011/)
    → **target: the deferred (sham → ranibizumab) arm reaches only +7.3 letters at
    month 12 and never catches up.** This is the delay-response calibration.
59. Heier JS, et al. **Ranibizumab for macular edema due to retinal vein occlusions: long-term follow-up in the HORIZON trial.** *Ophthalmology* 2012;119:802-9.
    [PMID 22301066](https://pubmed.ncbi.nlm.nih.gov/22301066/)
    → the loss of gains under reduced monitoring frequency.
60. Campochiaro PA, et al. **Long-term outcomes in patients with retinal vein occlusion treated with ranibizumab: the RETAIN study.** *Ophthalmology* 2014;121:209-19.
    [PMID 24112944](https://pubmed.ncbi.nlm.nih.gov/24112944/)
    → **target: about half of CRVO eyes still require injections at 4 years.** In
    the model these are the eyes whose Pc never crosses below Pc*.
61. Boyer D, et al. **Vascular endothelial growth factor Trap-Eye for macular edema secondary to central retinal vein occlusion: six-month results of the phase 3 COPERNICUS study.** *Ophthalmology* 2012;119:1024-32.
    [PMID 22440275](https://pubmed.ncbi.nlm.nih.gov/22440275/)
    → **target: +17.3 letters and −457 µm at week 24.**
62. Heier JS, et al. **Intravitreal aflibercept injection for macular edema due to central retinal vein occlusion: two-year results from the COPERNICUS study.** *Ophthalmology* 2014;121:1414-20.
    [PMID 24679444](https://pubmed.ncbi.nlm.nih.gov/24679444/)
    → **target: +13.0 letters at week 100 for immediate treatment versus +1.5 for
    the deferred arm.** The single strongest published statement that delay is
    paid permanently.
63. Holz FG, et al. **VEGF Trap-Eye for macular oedema secondary to central retinal vein occlusion: 6-month results of the phase III GALILEO study.** *Br J Ophthalmol* 2013;97:278-84.
    [PMID 23298885](https://pubmed.ncbi.nlm.nih.gov/23298885/)
    → **target: +18.0 letters at week 24.**
64. Korobelnik JF, et al. **Intravitreal aflibercept injection for macular edema resulting from central retinal vein occlusion: one-year results of the phase 3 GALILEO study.** *Ophthalmology* 2014;121:202-8.
    [PMID 24084497](https://pubmed.ncbi.nlm.nih.gov/24084497/)
65. Campochiaro PA, et al. **Intravitreal aflibercept for macular edema following branch retinal vein occlusion: the 24-week results of the VIBRANT study.** *Ophthalmology* 2015;122:538-44.
    [PMID 25315663](https://pubmed.ncbi.nlm.nih.gov/25315663/)
    → **target: +17.0 letters versus grid laser +12.2 at week 24.**
66. Scott IU, et al. **Effect of bevacizumab vs aflibercept on visual acuity among patients with macular edema due to central retinal vein occlusion: the SCORE2 randomized clinical trial.** *JAMA* 2017;317:2072-87.
    [PMID 28492910](https://pubmed.ncbi.nlm.nih.gov/28492910/)
    → **the most important negative result in this model.** Bevacizumab (KD 58 pM)
    is non-inferior to aflibercept (KD 0.49 pM): +18.6 versus +18.9 letters. A
    118-fold affinity difference produces no difference in outcome. Any model in
    which effect scales with affinity is falsified by this trial; the log-dose
    structure of `t_sup` is what survives it.
67. Scott IU, et al. **Comparison of monthly vs treat-and-extend regimens for individuals with macular edema who respond well to anti-VEGF treatment (SCORE2 report 12).** *JAMA Ophthalmol* 2018;136:337-45.
    [PMID 29476687](https://pubmed.ncbi.nlm.nih.gov/29476687/)
    → the treat-and-extend controller in the model, including the ±4-week step.
68. Tadayoni R, et al. **Efficacy and safety of faricimab for macular edema due to retinal vein occlusion: 24-week results from the BALATON and COMINO trials.** *Ophthalmology* 2024;131:950-60.
    [PMID 38280653](https://pubmed.ncbi.nlm.nih.gov/38280653/)
    → **target: faricimab non-inferior to aflibercept (+16.9 vs +17.5 letters).**
    Six monthly loading doses, which is what the model's `FAR_TAE` scenario uses.
69. Danzig CJ, et al. **Faricimab treat-and-extend dosing for macular edema due to retinal vein occlusion: 72-week results.** *Ophthalmol Retina* 2025.
    [PMID 40107501](https://pubmed.ncbi.nlm.nih.gov/40107501/)
    → the achievable extension intervals; the `imax` cap in the controller.
70. Sahni J, et al. **Simultaneous inhibition of angiopoietin-2 and vascular endothelial growth factor-A with faricimab in diabetic macular edema (BOULEVARD).** *Ophthalmology* 2019;126:1155-70.
    [PMID 30905643](https://pubmed.ncbi.nlm.nih.gov/30905643/)
    → the evidence that the Ang-2 arm contributes durability rather than peak effect.
71. Lanzetta P, et al. **Intravitreal aflibercept 8 mg in neovascular age-related macular degeneration (PULSAR): 48-week results.** *Lancet* 2024;403:1141-52.
    [PMID 38461841](https://pubmed.ncbi.nlm.nih.gov/38461841/)
    → the high-molar-dose arm. The model predicts 8 mg extends suppression from 96
    to 114 days — a 19% gain for a 4-fold dose, because duration is logarithmic.
72. Dugel PU, et al. **HAWK and HARRIER: ninety-six-week outcomes from the phase 3 trials of brolucizumab for neovascular age-related macular degeneration.** *Ophthalmology* 2021;128:89-99.
    [PMID 32574761](https://pubmed.ncbi.nlm.nih.gov/32574761/)
    → brolucizumab's molar dose advantage (230.8 nmol per injection).
73. Monés J, et al. **Risk of inflammation, retinal vasculitis, and retinal occlusion-related events with brolucizumab.** *Ophthalmology* 2021;128:1050-9.
    [PMID 33207259](https://pubmed.ncbi.nlm.nih.gov/33207259/)
    → the intraocular inflammation node in cluster 23.

## 9. 스테로이드 (Corticosteroids)

74. Haller JA, et al. **Randomized, sham-controlled trial of dexamethasone intravitreal implant in patients with macular edema due to retinal vein occlusion (GENEVA).** *Ophthalmology* 2010;117:1134-46.
    [Trial report; see also PMID 24447389 for the pooled safety analysis](https://pubmed.ncbi.nlm.nih.gov/24447389/)
    → **target: CST reduction and the ~16% rate of IOP rise ≥10 mmHg.** The model
    reproduces the CST effect only after the steroid is allowed a *direct*
    junction effect (`E_DEX_TJ`) in addition to mediator suppression — VEGF
    suppression alone is not enough to match the trial.
75. Kuppermann BD, et al. **Onset and duration of visual acuity improvement after dexamethasone intravitreal implant in eyes with macular edema due to retinal vein occlusion.** *Retina* 2014;34:1743-9.
    [PMID 24830824](https://pubmed.ncbi.nlm.nih.gov/24830824/)
    → the time course that constrains `K_REL_DEX` and the ~6-month re-dosing interval.
76. Boyer DS, et al. **Three-year, randomized, sham-controlled trial of dexamethasone intravitreal implant in patients with diabetic macular edema (MEAD).** *Ophthalmology* 2014;121:1904-14.
    [PMID 24907062](https://pubmed.ncbi.nlm.nih.gov/24907062/)
    → **the cataract calibration: ~68% cataract-related events over 3 years.** The
    reason the model runs phakic and pseudophakic arms separately.
77. Ip MS, et al. **A randomized trial comparing the efficacy and safety of intravitreal triamcinolone with observation to treat vision loss associated with macular edema secondary to central retinal vein occlusion (SCORE report 5).** *Arch Ophthalmol* 2009;127:1101-14.
    [PMID 19752419](https://pubmed.ncbi.nlm.nih.gov/19752419/)
78. Scott IU, et al. **A randomized trial comparing the efficacy and safety of intravitreal triamcinolone with standard care to treat vision loss associated with macular edema secondary to branch retinal vein occlusion (SCORE report 6).** *Arch Ophthalmol* 2009;127:1115-28.
    [PMID 19752420](https://pubmed.ncbi.nlm.nih.gov/19752420/)
79. Li X, et al. **Safety and efficacy of dexamethasone intravitreal implant for treatment of macular edema secondary to retinal vein occlusion in Chinese patients.** *Graefes Arch Clin Exp Ophthalmol* 2018;256:59-69.
    [PMID 29119239](https://pubmed.ncbi.nlm.nih.gov/29119239/)
80. Wecker T, et al. **Real-life medium term follow-up data for intravitreal dexamethasone implant in retinal vein occlusion.** *Sci Rep* 2021;11:8385.
    [PMID 33859243](https://pubmed.ncbi.nlm.nih.gov/33859243/)

## 10. 허혈 · 신생혈관 · 신생혈관 녹내장 (Ischaemia, neovascularisation, NVG)

81. Central Vein Occlusion Study Group. **A randomized clinical trial of early panretinal photocoagulation for ischemic central vein occlusion (report N).** *Ophthalmology* 1995;102:1434-44.
    [PMID 9097789](https://pubmed.ncbi.nlm.nih.gov/9097789/)
    → **the PRP node.** Prophylactic PRP does not prevent iris neovascularisation
    better than prompt PRP at the first sign of it; the model reproduces this
    because PRP removes the VEGF *source* but not the ischaemia already present.
82. Hayreh SS, Zimmerman MB. **Ocular neovascularization associated with central and hemicentral retinal vein occlusion.** *Retina* 2012;32:1553-65.
    [PMID 22495331](https://pubmed.ncbi.nlm.nih.gov/22495331/)
    → **target: the incidence and timing of INV/ANV in ischaemic CRVO.**
83. Rong AJ, Swaminathan SS, Vanner EA, Parrish RK. **Predictors of neovascular glaucoma in central retinal vein occlusion.** *Am J Ophthalmol* 2019;204:62-9.
    [PMID 31277883](https://pubmed.ncbi.nlm.nih.gov/31277883/)
    → the NVG risk factors; `KM_NVI_NP` and the `NVI > 0.5 & IOP > 25` criterion.
84. Ou WC, et al. **Longitudinal quantification of retinal nonperfusion in the macula of eyes with retinal vein occlusion.** *Ophthalmic Surg Lasers Imaging Retina* 2018;49:249-56.
    [PMID 29664983](https://pubmed.ncbi.nlm.nih.gov/29664983/)
    → **the growth rate of non-perfusion**, which sets `K_NP`.
85. Zhu ZY, et al. **Effect of anti-VEGF treatment on nonperfusion areas in ischemic retinopathy.** *Int J Ophthalmol* 2021;14:1647-54.
    [PMID 34804852](https://pubmed.ncbi.nlm.nih.gov/34804852/)
    → the partial, leukostasis-mediated protection anti-VEGF gives against
    non-perfusion progression, encoded via `W_LEUK_V`.
86. Kadomoto S, et al. **Nonperfusion area quantification in branch retinal vein occlusion: a widefield optical coherence tomography angiography study.** *Retina* 2021;41:1210-7.
    [PMID 33105300](https://pubmed.ncbi.nlm.nih.gov/33105300/)
    → the disc-area units used for `NP`.
87. Hayreh SS. **Photocoagulation for retinal vein occlusion.** *Prog Retin Eye Res* 2021;85:100964.
    [PMID 33713810](https://pubmed.ncbi.nlm.nih.gov/33713810/)

## 11. 구조 바이오마커와 시력 (Structural biomarkers and visual function)

88. Sun JK, et al. **Disorganization of the retinal inner layers as a predictor of visual acuity in eyes with center-involved diabetic macular edema.** *JAMA Ophthalmol* 2014;132:1309-16.
    [PMID 25058813](https://pubmed.ncbi.nlm.nih.gov/25058813/)
    → the `DRIL` state and why thickness alone under-predicts acuity.
89. Sun JK, et al. **Neural retinal disorganization as a robust marker of visual acuity in current and resolved diabetic macular edema.** *Diabetes* 2015;64:2560-70.
    [PMID 25633419](https://pubmed.ncbi.nlm.nih.gov/25633419/)
    → **the evidence that the structural damage persists after the oedema
    resolves.** This is the empirical content of `L_pr` being irreversible.
90. Song Y, et al. **Prognostic biomarkers in treatment-naïve central retinal vein occlusion with macular edema.** *Eur J Med Res* 2025;30:559.
    [PMID 40624719](https://pubmed.ncbi.nlm.nih.gov/40624719/)
91. Hirano Y, et al. **Multimodal imaging of microvascular abnormalities in retinal vein occlusion.** *J Clin Med* 2021;10:405.
    [PMID 33494354](https://pubmed.ncbi.nlm.nih.gov/33494354/)

## 12. 전신 위험인자와 혈압 (Systemic risk and blood pressure)

92. O'Mahoney PR, Wong DT, Ray JG. **Retinal vein occlusion and traditional risk factors for atherosclerosis.** *Arch Ophthalmol* 2008;126:692-9.
    [PMID 18474782](https://pubmed.ncbi.nlm.nih.gov/18474782/)
    → **hypertension in 63-73% of CRVO.** The prior on `PA` in the virtual
    population, and the reason the blood-pressure arrow matters at a population
    level and not just in individuals.
93. Voigt AM, et al. **Incidence of retinal vein occlusion and its association with mortality: results from the Gutenberg Health Study.** *Ophthalmology* 2025;132:—.
    [PMID 39971137](https://pubmed.ncbi.nlm.nih.gov/39971137/)
94. Yau JW, et al. **Retinal vein occlusion: an approach to diagnosis, systemic risk factors and management.** *Intern Med J* 2008;38:904-10.
    [PMID 19120547](https://pubmed.ncbi.nlm.nih.gov/19120547/)

## 13. 실제 진료 결과 (Real-world outcomes and treatment intensity)

95. Hunt A, et al. **Central retinal vein occlusion 36-month outcomes with anti-VEGF: the Fight Retinal Blindness! registry.** *Ophthalmol Retina* 2023;7:317-25.
    [PMID 36371040](https://pubmed.ncbi.nlm.nih.gov/36371040/)
    → **the target for the `REALWORLD` scenario**: injection counts and letter
    gains over 36 months outside a trial.
96. Spooner KL, et al. **Long-term outcomes of anti-VEGF treatment of retinal vein occlusion.** *Eye (Lond)* 2022;36:1194-201.
    [PMID 34117379](https://pubmed.ncbi.nlm.nih.gov/34117379/)
97. Spooner K, et al. **Prospective study of aflibercept for the treatment of persistent macular oedema secondary to retinal vein occlusion.** *Clin Exp Ophthalmol* 2020;48:53-60.
    [PMID 31498950](https://pubmed.ncbi.nlm.nih.gov/31498950/)
    → **the persistent-oedema population.** In the model these are the eyes whose
    residual thickness is a pressure floor, and the trial's modest benefit from
    switching agent is exactly what the model predicts for them.
98. Wecker T, et al. **Five-year visual acuity outcomes and injection patterns in patients with pro-re-nata treatments for AMD, DME and RVO.** *Br J Ophthalmol* 2017;101:353-9.
    [PMID 27215744](https://pubmed.ncbi.nlm.nih.gov/27215744/)
    → the letters-per-injection relationship the dose-response curve reproduces.
99. Bhandari S, et al. **Changes in 12-month outcomes over time for age-related macular degeneration, diabetic macular oedema and retinal vein occlusion.** *Eye (Lond)* 2023;37:1145-54.
    [PMID 35508721](https://pubmed.ncbi.nlm.nih.gov/35508721/)
100. Ciulla TA, et al. **Longer-term anti-VEGF therapy outcomes in neovascular age-related macular degeneration, diabetic macular edema, and retinal vein occlusion.** *Ophthalmol Retina* 2022;6:796-806.
     [PMID 35381391](https://pubmed.ncbi.nlm.nih.gov/35381391/)
     → the treatment-intensity / outcome relationship across three indications.
101. Kawakami S, et al. **Healing rate of macular edema secondary to branch retinal vein occlusion in two years after initiation of anti-VEGF therapy.** *PLoS One* 2023;18:e0280058.
     [PMID 36595494](https://pubmed.ncbi.nlm.nih.gov/36595494/)
     → the fraction of BRVO eyes that become injection-free; in the model, the
     fraction whose Pc crosses below Pc*.
102. Akdemir SC, et al. **Comparison of anti-VEGF results between non-ischemic and ischemic branch retinal vein occlusion.** *J Fr Ophtalmol* 2022;45:1013-20.
     [PMID 36127168](https://pubmed.ncbi.nlm.nih.gov/36127168/)
     → the ischaemic / non-ischaemic response gap the model derives from Pc.

## 14. 주사 관련 안전성 (Injection-related safety)

103. Chang-Lin JE (above, PMID 20702826) — steroid IOP kinetics.
104. Abbas K, et al. **Transient vision and intraocular pressure changes following anti-vascular endothelial growth factor injection.** *Can J Ophthalmol* 2026.
     [PMID 40987347](https://pubmed.ncbi.nlm.nih.gov/40987347/)
     → the transient IOP spike from injecting 0.05 mL into a 4 mL cavity.
105. Zhang C, et al. **Rates of endophthalmitis in prefilled versus nonprefilled syringes for intravitreal injections.** *Ophthalmol Retina* 2026;10:—.
     [PMID 40816641](https://pubmed.ncbi.nlm.nih.gov/40816641/)
     → the per-injection endophthalmitis rate that makes injection count a risk
     as well as a benefit.
106. Hanna RM, et al. **Intravitreal vascular endothelial growth factor inhibitors, hypertension, proteinuria, and renal injury: a concerning trend.** *Curr Opin Nephrol Hypertens* 2022;31:47-56.
     [PMID 34750330](https://pubmed.ncbi.nlm.nih.gov/34750330/)
     → the systemic VEGF-suppression safety node; connects to the model's `ASYS`
     compartment and the aflibercept-versus-ranibizumab plasma difference.

## 15. 모델링 방법 (Modelling methodology)

107. Park J, et al. **Computational fluid dynamics modeling of intravitreal ranibizumab bolus versus subretinal gene therapy.** *Transl Vis Sci Technol* 2025;14:—.
     [PMID 41042032](https://pubmed.ncbi.nlm.nih.gov/41042032/)
     → shows how much the well-mixed-vitreous assumption gives up; part of the
     justification for treating `FPEN` and `ALPHA` as lumped calibration factors.
108. Franz M, et al. **Injected eye, fellow eye and systemic pharmacokinetic modeling of intravitreally administered bevacizumab.** *Mol Pharm* 2026;23:—.
     [PMID 41925314](https://pubmed.ncbi.nlm.nih.gov/41925314/)
109. Chowdhury JM, et al. **Computational modeling of intravitreal ranibizumab kinetics: predicting macular drug concentration.** *PLoS One* 2026;21:e0—.
     [PMID 42096431](https://pubmed.ncbi.nlm.nih.gov/42096431/)
110. Wohlfart S, et al. **Vitrectomy status alters early intravitreal drug distribution: an ex vivo PET-CT study.** *Invest Ophthalmol Vis Sci* 2026;67:—.
     [PMID 42446480](https://pubmed.ncbi.nlm.nih.gov/42446480/)
111. Jing Y, et al. **Faricimab for retinal vein occlusion: a review of current evidence and future perspectives.** *Front Pharmacol* 2025;16:—.
     [PMID 41164764](https://pubmed.ncbi.nlm.nih.gov/41164764/)
112. Chaudhary V, et al. **Emerging clinical evidence of a dual role for Ang-2 and VEGF-A blockade with faricimab in retinal diseases.** *Graefes Arch Clin Exp Ophthalmol* 2025;263:1345-58.
     [PMID 39708087](https://pubmed.ncbi.nlm.nih.gov/39708087/)

---

## 모델 파라미터의 출처 요약 (Where each key parameter comes from)

| 파라미터 | 값 | 출처 / 근거 |
|---|---|---|
| `PA` post-arteriolar inlet pressure | 40 mmHg | Glucksberg 1993 [10], Attariwala 1994 [11] |
| `RV0` / `RA0` | 0.25 / 1.0 | chosen so that normal `Pc` = 21.6 mmHg [10] |
| `PT`, `PI_C`, `PI_T0`, `SIGMA_MAX` | 12, 25, 4 mmHg, 0.92 | Levick & Michel 2010 [13], Michel 2020 [14] |
| `OCC0` per phenotype | 5 – 32 | scaled to Attariwala 1997 experimental occlusion [12] |
| `OCC_RES` residual stenosis | 0.18 – 0.80 | fitted to the natural-history reviews [1,2] and RETAIN [60] |
| `K_COLL`, `COLL_MAX` | 0.0035 /d; 0.45 – 2.2 | Tsuboi 2019 OCT-A collateral kinetics [19,20] |
| `RV_VEGF`, `E_HIF_V` | 24 pM/d, ×12 | Noma 2009 aqueous / vitreous VEGF [35], Aiello 1994 [33] |
| `RI_IL6`, `E_HIF_I` | 30 pg/mL/d, ×16 | Noma 2008/2014 IL-6 [34,36] |
| `RA_ANG2`, `E_HIF_A` | 10 pM/d, ×7 | Saharinen 2017 [38], Sahni 2019 [70] |
| `TVIT` per agent | 5.1 – 9.8 d | Krohne 2008/2012 [41,42], Caruso 2020 [45] |
| `KDV` per agent | 0.49 – 58 pM | Papadopoulos 2012 [51] |
| `FPEN` retina:vitreous | 0.08 | Bakri 2007 [43], Hutton-Smith 2017 [47] |
| `ALPHA` apparent-IC50 multiplier | 120 | **calibrated**, single global value; see [46,52,66] |
| `K_NP` non-perfusion growth | 0.30 DA/d | Ou 2018 [84], Kadomoto 2021 [86] |
| `K_NVI`, NVG criterion | 0.030 /d; NVI>0.5 & IOP>25 | Hayreh 2012 [82], Rong 2019 [83] |
| dexamethasone implant PK | plateau ≈ 900 ng/mL | Chang-Lin 2011 [49] |
| `K_CAT_DEX` | 0.0024 /d | MEAD 3-year cataract rate [76] |
| `E_ED`, `KM_ED`, `E_PR` | 36 letters, 230 µm, 42 letters | CRUISE [56], COPERNICUS [61], Sun 2014/2015 [88,89] |
| `K_EZ` photoreceptor damage | 0.0018 /d | fitted to the deferred arms of CRUISE [58] and COPERNICUS [62] |

## 이 모델이 사용하지 않은 것 (What this model deliberately does not claim)

- **개별 환자 예측 (individual prediction).** The parameters are population-level
  and the calibration targets are trial means. Nothing here is fitted to any
  individual eye.
- **`ALPHA`는 기전이 아니라 보정 상수입니다.** The 120-fold gap between an SPR
  affinity and an apparent in-vivo IC50 is a lumped stand-in for VEGFR2
  competition, local interstitial flux and diffusional limitation. It is applied
  identically to every agent, so it cannot manufacture a difference between
  agents — but it is not a measured quantity, and any absolute suppression
  duration should be read as an order of magnitude.
- **절대 두께 (absolute thickness).** `KF`, `VMAX_PUMP` and `KM_PUMP` are fitted
  so that the CRVO baseline lands near the trial baselines. The *differences*
  between arms are the model's output; the absolute micrometres are calibration.
- **`SIGMA = SIGMAX·TJ^0.55` 관계.** The functional coupling of the reflection
  coefficient to junction integrity is a modelling choice with the right
  qualitative behaviour, not a measured relationship. It matters, because it is
  what sets `Pc*`.
