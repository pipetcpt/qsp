# Arrhythmogenic Right Ventricular Cardiomyopathy — References

Supporting literature for `arvc_qsp_model.dot`, `arvc_mrgsolve_model.R`,
`arvc_shiny_app.R` and `arvc_reference_impl.py`.

**How this list was built.** Every PMID below was resolved against the NCBI
E-utilities API (`esearch` + `esummary`) while the model was being written, so
each link goes to the record whose first author, year and journal are printed
beside it. Where a paper is load-bearing for a specific model choice, the
model file and parameter it feeds are named explicitly. Nothing in the model is
supported by a citation that is not on this list.

**Reading order if you only read five.** Basso 2009 (§1.1) for the disease,
James 2013 (§4.1) for the exercise dose-response that the whole model is built
around, La Gerche 2011 (§4.9) for the two wall-stress numbers that produce RV
selectivity, Fabritz 2011 (§4.7) for the load-reduction experiment the model
was *not* fitted to, and Marcus 2009 (§6.1) for the therapy hierarchy that
forces the two-generator structure.

---

## 1. Disease definition, pathology and natural history

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 1.1 | Basso C, et al. Arrhythmogenic right ventricular cardiomyopathy. *Lancet* 2009 | [19362677](https://pubmed.ncbi.nlm.nih.gov/19362677/) | Overall disease architecture; subepicardial-to-endocardial wavefront of fibrofatty replacement (map cluster 10, `SUBEPI`) |
| 1.2 | Corrado D, et al. Spectrum of clinicopathologic manifestations of ARVC/D: a multicenter study. *J Am Coll Cardiol* 1997 | [9362410](https://pubmed.ncbi.nlm.nih.gov/9362410/) | Regional distribution of disease ("triangle of dysplasia"); node `ANEURYSM` |
| 1.3 | Groeneweg JA, et al. Clinical presentation, long-term follow-up and outcomes of 1001 ARVD/C patients and family members. *Circ Cardiovasc Genet* 2015 | [25820315](https://pubmed.ncbi.nlm.nih.gov/25820315/) | **Calibration target** for `K_INJ`: median age at definite diagnosis; also the incomplete penetrance seen in relatives |
| 1.4 | Sen-Chowdhry S, et al. Clinical and genetic characterization of families with ARVD/C. *Circulation* 2007 | [17372169](https://pubmed.ncbi.nlm.nih.gov/17372169/) | Patterns of disease expression within families; the concealed phase |
| 1.5 | Sen-Chowdhry S, et al. Left-dominant arrhythmogenic cardiomyopathy: an under-recognized clinical entity. *J Am Coll Cardiol* 2008 | [19095136](https://pubmed.ncbi.nlm.nih.gov/19095136/) | Existence of a left-dominant phenotype → motivates `KAPPA_LV` as a genotype-gated switch |
| 1.6 | Marcus FI, et al. Diagnosis of ARVC/D: proposed modification of the Task Force Criteria. *Eur Heart J* 2010 | [20172912](https://pubmed.ncbi.nlm.nih.gov/20172912/) | The 2010 TFC, implemented literally in `observe()` (categories I–VI, definite/borderline/possible logic) |
| 1.7 | Marcus FI, et al. Same, *Circulation* 2010 (companion publication) | [20172911](https://pubmed.ncbi.nlm.nih.gov/20172911/) | Numeric thresholds: RVEDVi ≥110/100 mL/m², RVEF ≤40%, RVOT ≥32 mm, residual myocytes <60%/60–75%, TAD ≥55 ms, PVC >500/24 h |
| 1.8 | Corrado D, et al. Proposed diagnostic criteria for arrhythmogenic cardiomyopathy: European Task Force consensus report. *Int J Cardiol* 2024 | [37844667](https://pubmed.ncbi.nlm.nih.gov/37844667/) | Left-sided tissue characterisation; node `PADUA`; that left dominance tracks DSP/FLNC rather than PKP2 |
| 1.9 | Corrado D. Scarring/arrhythmogenic cardiomyopathy. *Eur Heart J Suppl* 2023 | [37125320](https://pubmed.ncbi.nlm.nih.gov/37125320/) | Reframing of the disease around scar rather than chamber |
| 1.10 | Cipriani A, et al. Cardiac MRI of arrhythmogenic cardiomyopathy: evolving diagnostic perspectives. *Eur Radiol* 2023 | [35788758](https://pubmed.ncbi.nlm.nih.gov/35788758/) | Imaging outputs `CMR_LGE`, `LGE_RING` |
| 1.11 | Arbelo E, et al. 2023 ESC Guidelines for the management of cardiomyopathies. *Eur Heart J* 2023 | [37622657](https://pubmed.ncbi.nlm.nih.gov/37622657/) | Current management framing; exercise recommendations; ICD indications |
| 1.12 | Kaski JP, et al. The 2023 ESC guidelines for cardiomyopathies: the 10 commandments. *Eur Heart J* 2024 | [38289320](https://pubmed.ncbi.nlm.nih.gov/38289320/) | Concise statement of the guideline positions the model is compared against |
| 1.13 | Towbin JA, et al. 2019 HRS expert consensus statement on evaluation, risk stratification and management of arrhythmogenic cardiomyopathy. *Heart Rhythm* 2019 | [31676023](https://pubmed.ncbi.nlm.nih.gov/31676023/) | Genotype-directed management; the ACM umbrella term |
| 1.14 | Corrado D, et al. Treatment of ARVC/D: an international Task Force consensus statement. *Circulation* 2015 | [26216213](https://pubmed.ncbi.nlm.nih.gov/26216213/) | Therapy ladder: exercise restriction → beta-blocker → AAD → ablation → ICD; the ordering the model has to reproduce |

## 2. Genetics and genotype–phenotype

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 2.1 | Gerull B, et al. Mutations in the desmosomal protein plakophilin-2 are common in ARVC. *Nat Genet* 2004 | [15489853](https://pubmed.ncbi.nlm.nih.gov/15489853/) | PKP2 as the dominant genotype; `PKP2_SET = 0.50` for truncating haploinsufficiency |
| 2.2 | McKoy G, et al. Identification of a deletion in plakoglobin in ARVC with palmoplantar keratoderma (Naxos disease). *Lancet* 2000 | [10902626](https://pubmed.ncbi.nlm.nih.gov/10902626/) | `JUP` node; plakoglobin as a disease gene, not only a signalling intermediate |
| 2.3 | Merner ND, et al. ARVC type 5 is a fully penetrant, lethal arrhythmic disorder caused by a missense mutation in TMEM43. *Am J Hum Genet* 2008 | [18313022](https://pubmed.ncbi.nlm.nih.gov/18313022/) | `TMEM43` node; the existence of a *fully* penetrant variant is the boundary case for the penetrance argument in §4 |
| 2.4 | Ortiz-Genga MF, et al. Truncating FLNC mutations are associated with high-risk dilated and arrhythmogenic cardiomyopathies. *J Am Coll Cardiol* 2016 | [27908349](https://pubmed.ncbi.nlm.nih.gov/27908349/) | `FLNC` → `KAPPA_LV` ON; ring-like LGE phenotype |
| 2.5 | Bhonsale A, et al. Impact of genotype on clinical course in ARVD/C-associated mutation carriers. *Eur Heart J* 2015 | [25616645](https://pubmed.ncbi.nlm.nih.gov/25616645/) | Genotype modifies course but does not determine it; supports genotype entering only through the reserve term |
| 2.6 | Carruth ED, et al. Prevalence and EHR-based phenotype of loss-of-function genetic variants in ARVC-associated genes. 2019 | [31638835](https://pubmed.ncbi.nlm.nih.gov/31638835/) | **Key observation (A):** unselected carriers of ARVC-gene LoF variants are mostly phenotype-negative. Falsifies a genotype-driven clock |
| 2.7 | Brandão M, et al. Desmoplakin cardiomyopathy: comprehensive review of an increasingly recognized entity. *J Clin Med* 2023 | [37048743](https://pubmed.ncbi.nlm.nih.gov/37048743/) | DSP as a distinct, LV-first, episodic-injury disease → `DSP_SET = 0.50`, `KAPPA_LV = 5.0` (which has to exceed the RV:LV load ratio of ~3.4 for the LV to actually go first) |
| 2.8 | Hoffman-Andrews L, et al. Desmoplakin cardiomyopathy: recent updates in natural history and management. *Curr Opin Cardiol* 2025 | [40600431](https://pubmed.ncbi.nlm.nih.gov/40600431/) | Current DSP natural history; hot-phase episodes |
| 2.9 | Delmar M, McKenna WJ. The cardiac desmosome and arrhythmogenic cardiomyopathies: from gene to disease. *Circ Res* 2010 | [20847325](https://pubmed.ncbi.nlm.nih.gov/20847325/) | The "connexome" framing of map cluster 2: one protein complex serving mechanical, electrical and signalling roles |

## 3. Molecular mechanism — junction, conduction, calcium, adipogenesis, inflammation

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 3.1 | Sato PY, et al. Loss of plakophilin-2 expression leads to decreased sodium current and slower conduction velocity in cultured cardiac myocytes. *Circ Res* 2009 | [19661460](https://pubmed.ncbi.nlm.nih.gov/19661460/) | `ETA_NAV`, `ECV`: PKP2 loss → I<sub>Na</sub> ↓ ≈30% → CV ↓; the mechanistic basis of Generator II arriving *before* scar |
| 3.2 | Cerrone M, et al. Missense mutations in plakophilin-2 cause sodium current deficit and associate with a Brugada syndrome phenotype. *Circulation* 2014 | [24352520](https://pubmed.ncbi.nlm.nih.gov/24352520/) | Same axis in vivo; PKP2 scaffolds the Nav1.5 complex at the ID crest (`NAV_CPLX`) |
| 3.3 | van Opbergen CJM, et al. Plakophilin-2 haploinsufficiency causes calcium handling deficits and modulates the cardiac response towards stress. *Int J Mol Sci* 2019 | [31438494](https://pubmed.ncbi.nlm.nih.gov/31438494/) | **Generator I in a structurally normal heart:** RyR2 destabilisation and diastolic Ca leak from PKP2 loss alone → `RYR_LEAK`, `CA_DIA` |
| 3.4 | Asimaki A, et al. A new diagnostic test for arrhythmogenic right ventricular cardiomyopathy. *N Engl J Med* 2009 | [19279339](https://pubmed.ncbi.nlm.nih.gov/19279339/) | Plakoglobin redistribution as a diagnostic signal → state `PG_NUC` |
| 3.5 | Garcia-Gras E, et al. Suppression of canonical Wnt/β-catenin signaling by nuclear plakoglobin recapitulates the phenotype of ARVC. *J Clin Invest* 2006 | [16823493](https://pubmed.ncbi.nlm.nih.gov/16823493/) | The Wnt brake on adipogenesis: `WNT_ACT`, and why the replacement tissue is fat |
| 3.6 | Chen SN, et al. The Hippo pathway is activated and is a causal mechanism for adipogenesis in arrhythmogenic cardiomyopathy. *Circ Res* 2014 | [24276085](https://pubmed.ncbi.nlm.nih.gov/24276085/) | `HIPPO`, `YAP_SEQ`, `PPARG`; and the requirement that **both** switches be thrown for adipogenesis (`ADIPO_DRIVE` is a product) |
| 3.7 | Chelko SP, et al. Central role for GSK3β in the pathogenesis of arrhythmogenic cardiomyopathy. *JCI Insight* 2016 | [27170944](https://pubmed.ncbi.nlm.nih.gov/27170944/) | `GSK3I` intervention node and `E_GSK_MAX` |
| 3.8 | Padrón-Barthe L, et al. Severe cardiac dysfunction and death caused by ARVC type 5 are improved by inhibition of glycogen synthase kinase-3β. *Circulation* 2019 | [31567019](https://pubmed.ncbi.nlm.nih.gov/31567019/) | Independent replication of the GSK-3β target in a different genotype (TMEM43) |
| 3.9 | Lombardi R, et al. Cardiac fibro-adipocyte progenitors express desmosome proteins and preferentially differentiate to adipocytes upon deletion of the desmoplakin gene. *Circ Res* 2016 | [27121621](https://pubmed.ncbi.nlm.nih.gov/27121621/) | `FAP` as the cell of origin for both fibrous and fatty replacement — one progenitor, two fates |
| 3.10 | Kim C, et al. Studying arrhythmogenic right ventricular dysplasia with patient-specific iPSCs. *Nature* 2013 | [23354045](https://pubmed.ncbi.nlm.nih.gov/23354045/) | PPARγ-dependent adipogenesis in human cells; metabolic co-requirement |
| 3.11 | Caspi O, et al. Modeling of arrhythmogenic right ventricular cardiomyopathy with human induced pluripotent stem cells. *Circ Cardiovasc Genet* 2013 | [24200905](https://pubmed.ncbi.nlm.nih.gov/24200905/) | Human cellular phenotype including stress-dependence |
| 3.12 | Chelko SP, et al. Therapeutic modulation of the immune response in arrhythmogenic cardiomyopathy. *Circulation* 2019 | [31533459](https://pubmed.ncbi.nlm.nih.gov/31533459/) | NF-κB as an active driver, not a bystander → `INF_STATE`, `PHI_INF`, and the `IL1BLK`/`GC` arms |
| 3.13 | Ariyaratne GHDN, et al. A paradigm shift: arrhythmogenic cardiomyopathy is an inflammatory disease. *Cells* 2026 | [42193878](https://pubmed.ncbi.nlm.nih.gov/42193878/) | Current statement of the inflammatory-amplifier position; hot phase |
| 3.14 | Zarrouk S, et al. Identification of biomarkers of arrhythmogenic cardiomyopathy by plasma proteomics. *Medicina (Kaunas)* 2025 | [39859087](https://pubmed.ncbi.nlm.nih.gov/39859087/) | Candidate circulating readouts for the inflammatory/fibrotic states |

## 4. The load hypothesis — exercise, wall stress, and the experiments that test it

This is the section the model's first commitment stands on.

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 4.1 | James CA, et al. Exercise increases age-related penetrance and arrhythmic risk in ARVD/C-associated desmosomal mutation carriers. *J Am Coll Cardiol* 2013 | [23871885](https://pubmed.ncbi.nlm.nih.gov/23871885/) | **Key observation (B)** and the model's headline prediction check: hazard ratio ≈3.16 for VT/death, highest vs lowest exercise tertile. The model predicts 3.0 from `NMECH = 4` and the two exercise stress coefficients, with nothing fitted |
| 4.2 | Ruwald AC, et al. Association of competitive and recreational sport participation with cardiac events in patients with ARVC. *Eur Heart J* 2015 | [25896080](https://pubmed.ncbi.nlm.nih.gov/25896080/) | Dose-graded risk and the benefit of reducing participation → the `RESTRICT` node acting on the clock |
| 4.3 | Sawant AC, et al. Exercise has a disproportionate role in the pathogenesis of ARVD/C in patients without desmosomal mutations. *J Am Heart Assoc* 2014 | [25516436](https://pubmed.ncbi.nlm.nih.gov/25516436/) | **Key observation (C):** gene-elusive patients had done *more* exercise. Load alone can reach the phenotype → `GENE_ELUSIVE` with normal reserve |
| 4.4 | Saberniak J, et al. Vigorous physical activity impairs myocardial function in ARVC patients and in mutation-positive family members. *Eur J Heart Fail* 2014 | [25319773](https://pubmed.ncbi.nlm.nih.gov/25319773/) | Exercise dose affects *function*, not only events — the structural limb of the clock |
| 4.5 | Lie ØH, et al. Harmful effects of exercise intensity and exercise duration in patients with arrhythmogenic cardiomyopathy. *JACC Clin Electrophysiol* 2018 | [29929667](https://pubmed.ncbi.nlm.nih.gov/29929667/) | Intensity **and** duration both matter — i.e. the exposure is a dose, which is what `f_ex × (1+K_EX)^NMECH` encodes |
| 4.6 | Kirchhof P, et al. Age- and training-dependent development of ARVC in heterozygous plakoglobin-deficient mice. *Circulation* 2006 | [17030684](https://pubmed.ncbi.nlm.nih.gov/17030684/) | Training-dependence of the phenotype in a controlled genetic model |
| 4.7 | Fabritz L, et al. Load-reducing therapy prevents development of arrhythmogenic right ventricular cardiomyopathy in plakoglobin-deficient mice. *J Am Coll Cardiol* 2011 | [21292134](https://pubmed.ncbi.nlm.nih.gov/21292134/) | **The external validation the model was not fitted to.** Furosemide + nitrate — no desmosomal, ion-channel or anti-fibrotic action — prevented RV enlargement and the arrhythmic phenotype. Reproduced by the `LOADRED` node acting only on the RV volume set point |
| 4.8 | Cruz FM, et al. Exercise triggers ARVC phenotype in mice expressing a disease-causing mutated version of human plakophilin-2. *J Am Coll Cardiol* 2015 | [25857910](https://pubmed.ncbi.nlm.nih.gov/25857910/) | Exercise as the trigger in a PKP2-specific model, i.e. the same effect in the model's index genotype |
| 4.9 | La Gerche A, et al. Disproportionate exercise load and remodeling of the athlete's right ventricle. *Med Sci Sports Exerc* 2011 | [21085033](https://pubmed.ncbi.nlm.nih.gov/21085033/) | **The two numbers behind commitment 2:** RV end-systolic wall stress rises far more than LV during exercise → `K_EX_RV = 1.25`, `K_EX_LV = 0.14` |
| 4.10 | La Gerche A, et al. Exercise-induced right ventricular dysfunction and structural remodelling in endurance athletes. *Eur Heart J* 2012 | [22160404](https://pubmed.ncbi.nlm.nih.gov/22160404/) | Acute post-exercise RV dysfunction and chronic remodelling in athletes without a variant |
| 4.11 | La Gerche A, et al. The response of the pulmonary circulation and right ventricle to exercise. *Pulm Circ* 2014 | [25621154](https://pubmed.ncbi.nlm.nih.gov/25621154/) | Why PVR falls less than SVR → `PVR`, `PAP`, `SVR` nodes |
| 4.12 | La Gerche A, Claessen G. Exercise and the right ventricle: a potential Achilles' heel. *Cardiovasc Res* 2017 | [28957535](https://pubmed.ncbi.nlm.nih.gov/28957535/) | Synthesis of the RV-load argument; the athlete's RV as both physiology and risk (`ATH_RV`) |
| 4.13 | Pelliccia A, et al. 2020 ESC Guidelines on sports cardiology and exercise in patients with cardiovascular disease. *Eur Heart J* 2021 | [32860412](https://pubmed.ncbi.nlm.nih.gov/32860412/) | The exercise prescription the model's `RESTRICT` scenario represents |
| 4.14 | Akdis D, et al. Sex hormones affect outcome in ARVC/D. *Eur Heart J* 2017 | [28329361](https://pubmed.ncbi.nlm.nih.gov/28329361/) | Basis for the single female multiplier `SEX_K_FEMALE = 0.72` on the fatigue rate; male predominance is then an output |

## 5. Arrhythmia mechanism, mapping and ablation (Generator II)

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 5.1 | Tschabrunn CM, et al. Isolated critical epicardial arrhythmogenic substrate abnormalities in patients with ARVC. *Heart Rhythm* 2022 | [34883271](https://pubmed.ncbi.nlm.nih.gov/34883271/) | **Why `ABL_ENDO` < `ABL_EPI`:** the critical substrate can be entirely epicardial, so endocardial mapping under-reads it |
| 5.2 | Berruezo A, et al. Combined endocardial and epicardial catheter ablation in ARVD incorporating scar dechanneling. *Circ Arrhythm Electrophysiol* 2012 | [22205683](https://pubmed.ncbi.nlm.nih.gov/22205683/) | Scar dechannelling; `CHANNELS`; the quantitative gap the model predicts between the two approaches |
| 5.3 | Berruezo A, et al. Safety, long-term outcomes and predictors of recurrence after first-line combined endo-epicardial VT substrate ablation. *Europace* 2017 | [28431051](https://pubmed.ncbi.nlm.nih.gov/28431051/) | Durability of the combined approach |
| 5.4 | Fernández-Armenta J, et al. Sinus rhythm detection of conducting channels and VT isthmus in ARVC. *Heart Rhythm* 2014 | [24561159](https://pubmed.ncbi.nlm.nih.gov/24561159/) | The conducting-channel concept behind `SCAR_HET` peaking at intermediate replacement, not at maximal scar |
| 5.5 | Santangeli P, et al. Long-term outcome with catheter ablation of ventricular tachycardia in patients with ARVC. *Circ Arrhythm Electrophysiol* 2015 | [26546346](https://pubmed.ncbi.nlm.nih.gov/26546346/) | Endocardial-only recurrence rates; calibration anchor for `ABL_ENDO = 0.35` |
| 5.6 | Santangeli P, et al. Outcomes of catheter ablation in ARVC without background ICD therapy. *JACC Clin Electrophysiol* 2019 | [30678787](https://pubmed.ncbi.nlm.nih.gov/30678787/) | Ablation as substrate therapy independent of the device |
| 5.7 | Bai R, et al. Ablation of ventricular arrhythmias in ARVD/C: arrhythmia-free survival after endo-epicardial substrate-based mapping and ablation. *Circ Arrhythm Electrophysiol* 2011 | [21665983](https://pubmed.ncbi.nlm.nih.gov/21665983/) | Independent series supporting the endo-vs-endo/epi gap |
| 5.8 | Alahwany SH, et al. Outcomes of VT catheter ablation in paediatric ARVC. *Circ Arrhythm Electrophysiol* 2025 | [41078126](https://pubmed.ncbi.nlm.nih.gov/41078126/) | Substrate-directed ablation in the young, where the model puts Generator I ahead of Generator II |
| 5.9 | Zeppenfeld K, et al. 2022 ESC Guidelines for the management of patients with ventricular arrhythmias and the prevention of sudden cardiac death. *Eur Heart J* 2022 | [36017572](https://pubmed.ncbi.nlm.nih.gov/36017572/) | Guideline placement of AAD, ablation and ICD in ARVC |
| 5.10 | Könemann H, et al. Spotlight on the 2022 ESC guideline management of ventricular arrhythmias and prevention of sudden cardiac death. *Europace* 2023 | [37102266](https://pubmed.ncbi.nlm.nih.gov/37102266/) | Condensed statement of the same recommendations |
| 5.11 | Maharani E, et al. Diagnostic accuracy of electrocardiographic criteria in differentiating arrhythmogenic RV disease. *J Arrhythm* 2026 | [42524098](https://pubmed.ncbi.nlm.nih.gov/42524098/) | ECG discrimination from idiopathic RVOT VT — the `PHENOCOPY` node |

## 6. Drug therapy — the hierarchy the two-generator structure has to reproduce

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 6.1 | Marcus GM, et al. Efficacy of antiarrhythmic drugs in ARVC: a report from the North American ARVC Registry. *J Am Coll Cardiol* 2009 | [19660690](https://pubmed.ncbi.nlm.nih.gov/19660690/) | **Key observation (E)** and the anchor for `LAM2`: amiodarone effective, sotalol not better than nothing, beta-blockers without a demonstrated VT benefit. The sotalol null is a model *prediction* (its IKr and proarrhythmic arms cancel); the beta-blocker result is the model's **stated miss** |
| 6.2 | Wichter T, et al. Efficacy of antiarrhythmic drugs in patients with arrhythmogenic right ventricular disease. *Circulation* 1992 | [1617780](https://pubmed.ncbi.nlm.nih.gov/1617780/) | The earlier, more favourable sotalol experience — the discrepancy the two-generator model is asked to explain |
| 6.3 | Ermakov S, et al. Use of flecainide in combination antiarrhythmic therapy in patients with ARVC. *Heart Rhythm* 2017 | [27939893](https://pubmed.ncbi.nlm.nih.gov/27939893/) | Flecainide as an add-on reduces VT in a structural cardiomyopathy → the RyR2 arm of `FLEC`, and the reason the model gives the drug two signs |
| 6.4 | Ermakov S, et al. Combination drug therapy for patients with intractable VT associated with right ventricular cardiomyopathy. *Pacing Clin Electrophysiol* 2014 | [24102153](https://pubmed.ncbi.nlm.nih.gov/24102153/) | Combination therapy in refractory disease |
| 6.5 | Ermakov S, Scheinman M. Arrhythmogenic right ventricular cardiomyopathy — antiarrhythmic therapy. *Arrhythm Electrophysiol Rev* 2015 | [26835106](https://pubmed.ncbi.nlm.nih.gov/26835106/) | Review of the drug hierarchy, including amiodarone's long-term toxicity burden in young patients (`AMIO_TOX`) |
| 6.6 | Waldo AL, et al. Effect of d-sotalol on mortality in patients with left ventricular dysfunction after recent and remote myocardial infarction (SWORD). *Lancet* 1996 | [8691967](https://pubmed.ncbi.nlm.nih.gov/8691967/) | **Sets the size of the `DISP_REPOL` term.** A randomised trial in which pure IKr block *increased* mortality in a scarred ventricle. Any model making IKr block net-beneficial in scar is contradicted by a trial, not merely by a registry |
| 6.7 | The Cardiac Arrhythmia Suppression Trial (CAST) Investigators. Preliminary report: effect of encainide and flecainide on mortality after myocardial infarction. *N Engl J Med* 1989 | [2473403](https://pubmed.ncbi.nlm.nih.gov/2473403/) | The class-IC caution the flecainide node has to reproduce. In this model it emerges from the conduction-velocity term feeding Generator II, not from a rule forbidding the drug |

## 7. Risk stratification and devices

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 7.1 | Cadrin-Tourigny J, et al. A new prediction model for ventricular arrhythmias in ARVC. *Eur Heart J* 2019 | [30915475](https://pubmed.ncbi.nlm.nih.gov/30915475/) | The seven predictors the model deliberately **emits as outputs** (age, sex, syncope, NSVT, PVC/24 h, number of leads with TWI, RVEF) so its trajectories can be scored externally rather than by an internal score |
| 7.2 | Carrick RT, et al. Longitudinal prediction of ventricular arrhythmic risk in patients with ARVC. *Circ Arrhythm Electrophysiol* 2022 | [36315818](https://pubmed.ncbi.nlm.nih.gov/36315818/) | Risk is a *trajectory*, not a baseline — matching a model whose hazard is a function of moving states |
| 7.3 | Carrick RT, et al. A novel tool for arrhythmic risk stratification in desmoplakin gene variant carriers. *Eur Heart J* 2024 | [39011630](https://pubmed.ncbi.nlm.nih.gov/39011630/) | Genotype-specific risk structure, consistent with `KAPPA_LV` being genotype-gated |
| 7.4 | Bosman LP, et al. Predicting arrhythmic risk in ARVC: a systematic review and meta-analysis. *Heart Rhythm* 2018 | [29408436](https://pubmed.ncbi.nlm.nih.gov/29408436/) | Pooled event rates used as the sanity range for `H0_VA` |
| 7.5 | Orgeron GM, et al. ICD therapy in ARVD/C: predictors of appropriate therapy, outcomes and complications. *J Am Heart Assoc* 2017 | [28588093](https://pubmed.ncbi.nlm.nih.gov/28588093/) | Appropriate-therapy rates and complication burden → `ICD_EFF`, `H_ICD_CX`, `R_ICD_INAPPROP` |
| 7.6 | Bhonsale A, et al. Incidence and predictors of ICD therapy in ARVD/C patients undergoing ICD implantation for primary prevention. *J Am Coll Cardiol* 2011 | [21939834](https://pubmed.ncbi.nlm.nih.gov/21939834/) | Primary- vs secondary-prevention event rates |
| 7.7 | Migliore F, et al. Subcutaneous ICD in patients with ARVC. *Int J Cardiol* 2019 | [30661851](https://pubmed.ncbi.nlm.nih.gov/30661851/) | `SICD` node and its trade-offs in young patients |
| 7.8 | De Marco C, et al. Left ventricular late gadolinium enhancement for arrhythmic risk prediction in ARVC. *Circ Arrhythm Electrophysiol* 2026 | [41608798](https://pubmed.ncbi.nlm.nih.gov/41608798/) | LV substrate as an independent risk axis → the `SCAR_HET_LV` contribution to Generator II |

## 8. Emerging and mechanism-directed therapy

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 8.1 | Bradford WH, et al. Plakophilin 2 gene therapy prevents and rescues ARVC in a mouse model harbouring patient variants. *Nat Cardiovasc Res* 2023 | [39196150](https://pubmed.ncbi.nlm.nih.gov/39196150/) | The `AAV` node: PKP2 restoration prevents *and* partially rescues. In the model it restores the reserve term and therefore the *rate*, which is why timing dominates dose |
| 8.2 | Lin X, et al. Fibroblast growth factor 21 prevents catecholaminergic arrhythmias in a mouse model of PKP2 arrhythmogenic cardiomyopathy. *Heart Rhythm* 2026 | [41759869](https://pubmed.ncbi.nlm.nih.gov/41759869/) | Independent evidence that the *catecholaminergic* (Generator I) limb is separately druggable |
| 8.3 | Chelko SP, et al. (see 3.12) — immune modulation | [31533459](https://pubmed.ncbi.nlm.nih.gov/31533459/) | The `IL1BLK` / `GC` arms as amplifier therapy, explicitly not clock therapy |
| 8.4 | Conte G, et al. Diagnostic, pharmacological and ablation approaches for idiopathic ventricular fibrillation: 2024 EHRA statement. *Europace* 2025 | [40394989](https://pubmed.ncbi.nlm.nih.gov/40394989/) | Context for VF without overt substrate — the phenotype Generator I predicts in concealed-phase carriers |

## 9. Differential diagnosis / phenocopies

| # | Citation | PubMed | What it supports here |
|---|----------|--------|-----------------------|
| 9.1 | Perazzolo Marra M, et al. Morphofunctional abnormalities of mitral annulus and arrhythmic mitral valve prolapse. *Circ Cardiovasc Imaging* 2016 | [27516479](https://pubmed.ncbi.nlm.nih.gov/27516479/) | A structural-arrhythmic phenocopy with LV scar; the `PHENOCOPY` node |
| 9.2 | Corrado D 2024 European Task Force criteria (see 1.8) | [37844667](https://pubmed.ncbi.nlm.nih.gov/37844667/) | Explicit exclusion requirements (sarcoidosis, myocarditis, sport-related remodelling) |
| 9.3 | Cipriani A 2023 (see 1.10) | [35788758](https://pubmed.ncbi.nlm.nih.gov/35788758/) | Imaging discrimination from myocarditis and sarcoidosis |

---

## What the model claims each citation is *not* used for

Being explicit about this matters more than the length of the list:

- The sotalol result is **not** fitted to the ARVC registry. The size of the
  repolarisation-dispersion term is anchored on SWORD (§6.6), a randomised
  mortality trial, and the ARVC registry's null (§6.1) is then a prediction.
  Likewise the flecainide sign flip is anchored on nothing at all — it falls
  out of the same conduction-velocity term that produces Generator II.
- **Nothing** in §4 was used to fit `NMECH`. The fatigue exponent was fixed at 4
  a priori (the Basquin exponent range for load-bearing biological tissue) and
  the ≈3-fold exercise hazard ratio of James 2013 fell out of it. If `NMECH`
  had been fitted, §4.1 would be a calibration target and not a test.
- `K_EX_RV = 1.25` and `K_EX_LV = 0.14` come from exercise CMR (§4.9) and are
  **not** adjusted to make the RV go first. RV selectivity is a consequence.
- Fabritz 2011 (§4.7) is a *prediction check*, not an input. `E_LOADRED_MAX`
  was set to a plausible diuretic/nitrate effect on RV volume, and the
  prevention of the phenotype is what the equations then do.
- Marcus 2009 (§6.1) supplied one fitted number (`LAM2`, the Generator II
  weight, set from the amiodarone-versus-sotalol gap). The sotalol null itself
  is not fitted — it emerges because the drug's two arms cancel. The
  beta-blocker result in the same paper is reported by the model as a **miss**
  rather than absorbed.
- Cadrin-Tourigny 2019 (§7.1) coefficients are deliberately **not** reproduced
  inside the model. The model emits the calculator's seven inputs so that an
  external, independently validated score can be applied to its trajectories;
  transcribing a published linear predictor would have made the comparison
  circular.

## Disclaimer

Educational and research QSP model. Semi-quantitative, assembled from public
literature, not independently validated, and **not for clinical
decision-making**. The AAV-PKP2, GSK-3β-inhibition, IL-1-blockade and
load-reducing-therapy scenarios are investigational or extrapolated from animal
data and are simulated here only to explore model structure.
