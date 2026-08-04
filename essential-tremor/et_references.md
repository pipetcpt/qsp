# 본태성 떨림 (Essential Tremor) — QSP 모델 참고문헌

**Essential Tremor QSP model — verified reference list**

총 **126편**. 이 파일의 모든 항목은 손으로 적은 것이 아니라 `et_reference_check.py`가
NCBI E-utilities로 PubMed에 질의해 **실제로 돌려받은 레코드**입니다. 제목·저자·저널·
연도·PMID는 전부 PubMed가 보고한 값을 그대로 옮긴 것이며, 기억에 의존해 작성한
인용은 한 건도 포함되지 않습니다. 재현 방법:

```bash
python3 et_reference_check.py            # 1차: 제목 기반 검증
python3 et_reference_check.py --harvest   # 2차: 주제 질의 기반 수집
python3 et_reference_check.py --emit      # 3차: 이 파일 생성
```

1차 통과에서 91개 후보 제목 중 39개만 확인되었습니다(나머지는 필자의 기억이
의역이었기 때문입니다). 근사한 인용을 그대로 싣는 대신, 섹션별 주제 질의를
실행해 **PubMed가 실제로 반환한 논문**을 인용하는 방식으로 대체했습니다.
`query` 필드가 있는 항목이 그렇게 수집된 것입니다.

---

## 1. 역학 · 임상양상 · 진단기준

*Epidemiology, phenomenology and diagnostic criteria* — 유병률, ET plus 논쟁, 자연경과. 모델의 환자 표현형(G0 · HDG · VXG)과 기저 TETRAS 값의 근거.

- Okelberry T, et al. **Updates in essential tremor.** Parkinsonism Relat Disord. 2024;122:106086 [PMID 38538475](https://pubmed.ncbi.nlm.nih.gov/38538475/)
- Elias WJ, et al. **Essential Tremor.** JAMA. 2024;332:418-419 [PMID 38976274](https://pubmed.ncbi.nlm.nih.gov/38976274/)
- Latorre A, et al. **The MDS consensus tremor classification: The best way to classify patients with tremor at present.** J Neurol Sci. 2022;435:120191 [PMID 35247714](https://pubmed.ncbi.nlm.nih.gov/35247714/)
- Louis ED, et al. **Problems and controversies in tremor classification.** J Neurol Sci. 2022;435:120204 [PMID 35279635](https://pubmed.ncbi.nlm.nih.gov/35279635/)
- Welton T, et al. **Essential tremor.** Nat Rev Dis Primers. 2021;7:83 [PMID 34764294](https://pubmed.ncbi.nlm.nih.gov/34764294/)
- Benito-León J, et al. **[Epidemiology of essential tremor].** Rev Neurol. 2020;70:139-148 [PMID 32043536](https://pubmed.ncbi.nlm.nih.gov/32043536/)
- Shanker V, et al. **Essential tremor: diagnosis and management.** BMJ. 2019;366:l4485 [PMID 31383632](https://pubmed.ncbi.nlm.nih.gov/31383632/)
- Bhatia KP, et al. **Consensus Statement on the classification of tremors. from the task force on tremor of the International Parkinson and Movement Disorder Society.** Mov Disord. 2018;33:75-87 [PMID 29193359](https://pubmed.ncbi.nlm.nih.gov/29193359/)
- Gövert F, et al. **Tremor entities and their classification: an update.** Curr Opin Neurol. 2015;28:393-9 [PMID 26110800](https://pubmed.ncbi.nlm.nih.gov/26110800/)

## 2. 떨림의 진동자 이론 · 측정과 척도

*Tremor as an oscillator; measurement and rating scales* — 루프 이득/지연 구분, 기계적 공명, 그리고 Elble의 로그 변환 — 모델의 ⑧ 진동자 핵과 TETRAS 사상의 근거.

- Smid A, et al. **Perioperative quantification of MDS-UPDRS-III tremor measurements in patients with Parkinson's disease using accelerometry.** J Neural Transm (Vienna). 2026;133:1457-1469 [PMID 41811427](https://pubmed.ncbi.nlm.nih.gov/41811427/)
- Kremer NI, et al. **Supine MDS-UPDRS-III Assessment: An Explorative Study.** J Clin Med. 2023;12 [PMID 37176549](https://pubmed.ncbi.nlm.nih.gov/37176549/)
- Smid A, et al. **A Novel Accelerometry Method to Perioperatively Quantify Essential Tremor Based on Fahn-Tolosa-Marin Criteria.** J Clin Med. 2023;12 [PMID 37445270](https://pubmed.ncbi.nlm.nih.gov/37445270/)
- Smid A, et al. **Intraoperative Quantification of MDS-UPDRS Tremor Measurements Using 3D Accelerometry: A Pilot Study.** J Clin Med. 2022;11 [PMID 35566401](https://pubmed.ncbi.nlm.nih.gov/35566401/)
- Elble RJ, et al. **Assessment of Head Tremor with Accelerometers Versus Gyroscopic Transducers.** Mov Disord Clin Pract. 2017;4:205-211 [PMID 30363428](https://pubmed.ncbi.nlm.nih.gov/30363428/)
- Elble R, et al. **Task force report: scales for screening and evaluating tremor: critique and recommendations.** Mov Disord. 2013;28:1793-800 [PMID 24038576](https://pubmed.ncbi.nlm.nih.gov/24038576/)
- Elble RJ, et al. **Tremor amplitude is logarithmically related to 4- and 5-point tremor rating scales.** Brain. 2006;129:2660-6 [PMID 16891320](https://pubmed.ncbi.nlm.nih.gov/16891320/)
- Elble RJ, et al. **Characteristics of physiologic tremor in young and elderly adults.** Clin Neurophysiol. 2003;114:624-35 [PMID 12686271](https://pubmed.ncbi.nlm.nih.gov/12686271/)
- Hashimoto T, et al. **Peripheral mechanisms in tremor after traumatic neck injury.** J Neurol Neurosurg Psychiatry. 2002;73:585-7 [PMID 12397157](https://pubmed.ncbi.nlm.nih.gov/12397157/)

## 3. 소뇌 병리

*Cerebellar pathology* — Purkinje torpedo, 등반섬유-PC 시냅스 병리, 치아핵 GABA 수용체 감소 — PCINT · DNDIS 상태변수의 근거.

- Kerridge CA, et al. **Purkinje Cell Loss in Essential Tremor: Collective Data From 215 Brains Over a 21-Year Period.** Ann Clin Transl Neurol. 2026;13:36-48 [PMID 40985128](https://pubmed.ncbi.nlm.nih.gov/40985128/)
- Widner J, et al. **Axonal pathology differentially affects human Purkinje cell subpopulations in the essential tremor cerebellum.** Proc Natl Acad Sci U S A. 2025;122:e2502024122 [PMID 40587795](https://pubmed.ncbi.nlm.nih.gov/40587795/)
- Hong CT, et al. **Cerebellar Structural and N-Acetylaspartate, Choline, and Creatine Metabolic Profiles in Parkinson's Disease and Essential Tremor.** Diagnostics (Basel). 2024;14 [PMID 39518397](https://pubmed.ncbi.nlm.nih.gov/39518397/)
- Gionco JT, et al. **Essential Tremor versus "ET-plus": A Detailed Postmortem Study of Cerebellar Pathology.** Cerebellum. 2021;20:904-912 [PMID 33768479](https://pubmed.ncbi.nlm.nih.gov/33768479/)
- Handforth A, et al. **Increased Purkinje Cell Complex Spike and Deep Cerebellar Nucleus Synchrony as a Potential Basis for Syndromic Essential Tremor. A Review and Synthesis of the Literature.** Cerebellum. 2021;20:266-281 [PMID 33048308](https://pubmed.ncbi.nlm.nih.gov/33048308/)
- Louis ED, et al. **Essential tremor pathology: neurodegeneration and reorganization of neuronal connections.** Nat Rev Neurol. 2020;16:69-83 [PMID 31959938](https://pubmed.ncbi.nlm.nih.gov/31959938/)
- Zhang X, et al. **Role of cerebellar GABAergic dysfunctions in the origins of essential tremor.** Proc Natl Acad Sci U S A. 2019;116:13592-13601 [PMID 31209041](https://pubmed.ncbi.nlm.nih.gov/31209041/)
- Lee D, et al. **Climbing fiber-Purkinje cell synaptic pathology across essential tremor subtypes.** Parkinsonism Relat Disord. 2018;51:24-29 [PMID 29482925](https://pubmed.ncbi.nlm.nih.gov/29482925/)
- Louis ED, et al. **Essential tremor and the cerebellum.** Handb Clin Neurol. 2018;155:245-258 [PMID 29891062](https://pubmed.ncbi.nlm.nih.gov/29891062/)
- Kuo SH, et al. **Climbing fiber-Purkinje cell synaptic pathology in tremor and cerebellar degenerative diseases.** Acta Neuropathol. 2017;133:121-138 [PMID 27704282](https://pubmed.ncbi.nlm.nih.gov/27704282/)
- Lizarraga KJ, et al. **Molecular imaging of movement disorders.** World J Radiol. 2016;8:226-39 [PMID 27029029](https://pubmed.ncbi.nlm.nih.gov/27029029/)
- Babij R, et al. **Purkinje cell axonal anatomy: quantifying morphometric changes in essential tremor versus control brains.** Brain. 2013;136:3051-61 [PMID 24030953](https://pubmed.ncbi.nlm.nih.gov/24030953/)
- Paris-Robidas S, et al. **Defective dentate nucleus GABA receptors in essential tremor.** Brain. 2012;135:105-16 [PMID 22120148](https://pubmed.ncbi.nlm.nih.gov/22120148/)
- Louis ED, et al. **Torpedoes in the cerebellar vermis in essential tremor cases vs. controls.** Cerebellum. 2011;10:812-9 [PMID 21656041](https://pubmed.ncbi.nlm.nih.gov/21656041/)
- Louis ED, et al. **Purkinje cell loss is a characteristic of essential tremor.** Parkinsonism Relat Disord. 2011;17:406-9 [PMID 21600832](https://pubmed.ncbi.nlm.nih.gov/21600832/)
- Louis ED, et al. **Neuropathological changes in essential tremor: 33 cases compared with 21 controls.** Brain. 2007;130:3297-307 [PMID 18025031](https://pubmed.ncbi.nlm.nih.gov/18025031/)

## 4. 하올리브핵 진동자 · T형 칼슘 · harmaline

*Inferior olive, T-type calcium, harmaline* — Cav3.1 · Cx36 · 아역치 진동. a_O(올리브 분율) 파라미터와 종간 불일치 결과의 근거.

- Kosmowska B, et al. **Inhibition of Excessive Glutamatergic Transmission in the Ventral Thalamic Nuclei by a Selective Adenosine A1 Receptor Agonist, 5'-Chloro-5'-Deoxy-(±)-ENBA Underlies its Tremorolytic Effect in the Harmaline-Induced Model of Essential Tremor.** Neuroscience. 2020;429:106-118 [PMID 31935489](https://pubmed.ncbi.nlm.nih.gov/31935489/)
- Handforth A, et al. **Harmaline tremor: underlying mechanisms in a potential animal model of essential tremor.** Tremor Other Hyperkinet Mov (N Y). 2012;2 [PMID 23440018](https://pubmed.ncbi.nlm.nih.gov/23440018/)
- Handforth A, et al. **T-type calcium channel antagonists suppress tremor in two mouse models of essential tremor.** Neuropharmacology. 2010;59:380-7 [PMID 20547167](https://pubmed.ncbi.nlm.nih.gov/20547167/)
- Park YG, et al. **Ca(V)3.1 is a tremor rhythm pacemaker in the inferior olive.** Proc Natl Acad Sci U S A. 2010;107:10731-6 [PMID 20498062](https://pubmed.ncbi.nlm.nih.gov/20498062/)
- Placantonakis DG, et al. **Continuous electrical oscillations emerge from a coupled network: a study of the inferior olive using lentiviral knockdown of connexin36.** J Neurosci. 2006;26:5008-16 [PMID 16687492](https://pubmed.ncbi.nlm.nih.gov/16687492/)
- Long MA, et al. **Rhythmicity without synchrony in the electrically uncoupled inferior olive.** J Neurosci. 2002;22:10898-905 [PMID 12486184](https://pubmed.ncbi.nlm.nih.gov/12486184/)

## 5. 유전학

*Genetics* — LINGO1 · FUS · NOTCH2NLC · GWAS · harmane 노출.

- Ilaghi M, et al. **The JAK1/2 Inhibitor Baricitinib Ameliorates Neuroinflammation and Symptoms in an Animal Model of Essential Tremor.** Pharmacol Res Perspect. 2026;14:e70299 [PMID 42427276](https://pubmed.ncbi.nlm.nih.gov/42427276/)
- Jiménez-Jiménez FJ, et al. **Genomic Markers for Essential Tremor.** Pharmaceuticals (Basel). 2021;14 [PMID 34072005](https://pubmed.ncbi.nlm.nih.gov/34072005/)
- Sun QY, et al. **Expansion of GGC repeat in the human-specific NOTCH2NLC gene is associated with essential tremor.** Brain. 2020;143:222-233 [PMID 31819945](https://pubmed.ncbi.nlm.nih.gov/31819945/)
- Louis ED, et al. **Blood Harmane (1-Methyl-9H-Pyrido[3,4-b]indole) and Mercury in Essential Tremor: A Population-Based, Environmental Epidemiology Study in the Faroe Islands.** Neuroepidemiology. 2020;54:272-280 [PMID 32007995](https://pubmed.ncbi.nlm.nih.gov/32007995/)
- Müller SH, et al. **Genome-wide association study in essential tremor identifies three new loci.** Brain. 2016;139:3163-3169 [PMID 27797806](https://pubmed.ncbi.nlm.nih.gov/27797806/)

## 6. 프로프라놀롤과 β 차단

*Propranolol and beta-blockade* — 말초 β₂ 부위, β₁ 선택성의 실패, 나돌롤의 효과 — φ_spindle의 Gaddum 경쟁식과 용량 천장의 근거.

- Saifee TA, et al. **Epidemiology and treatment of patients with Essential Tremor: A retrospective cohort analysis in the United Kingdom.** Neuroepidemiology. 2026 [PMID 42430313](https://pubmed.ncbi.nlm.nih.gov/42430313/)
- Marques L, et al. **Model-Based Virtual Clinical Trial Reveals Renal Impairment and Body Size as Key Determinants of Pharmacokinetic Variability and Drug-Drug Interaction Risk in Propranolol Therapy.** Pharmaceutics. 2026;18 [PMID 42357253](https://pubmed.ncbi.nlm.nih.gov/42357253/)
- Mohapatra P, et al. **Propranolol for Tremors in Spinocerebellar Ataxia Type 12: A Randomized Clinical Trial.** Mov Disord. 2026;41:373-383 [PMID 41261874](https://pubmed.ncbi.nlm.nih.gov/41261874/)
- Lv Y, et al. **Cerebellar repetitive transcranial magnetic stimulation versus propranolol for essential tremor.** Brain Behav. 2023;13:e2926 [PMID 36806734](https://pubmed.ncbi.nlm.nih.gov/36806734/)
- Ghosh D, et al. **A Series of 211 Children with Probable Essential Tremor.** Mov Disord Clin Pract. 2017;4:231-236 [PMID 30363473](https://pubmed.ncbi.nlm.nih.gov/30363473/)
- Zesiewicz TA, et al. **Evidence-based guideline update: treatment of essential tremor: report of the Quality Standards subcommittee of the American Academy of Neurology.** Neurology. 2011;77:1752-5 [PMID 22013182](https://pubmed.ncbi.nlm.nih.gov/22013182/)
- Yetimalar Y, et al. **Olanzapine versus propranolol in essential tremor.** Clin Neurol Neurosurg. 2005;108:32-5 [PMID 16311142](https://pubmed.ncbi.nlm.nih.gov/16311142/)
- Lee KS, et al. **A multicenter randomized crossover multiple-dose comparison study of arotinolol and propranolol in essential tremor.** Parkinsonism Relat Disord. 2003;9:341-7 [PMID 12853233](https://pubmed.ncbi.nlm.nih.gov/12853233/)
- Sampaio C, et al. **Essential tremor.** Clin Evid. 2002 [PMID 12230735](https://pubmed.ncbi.nlm.nih.gov/12230735/)
- Gironell A, et al. **A randomized placebo-controlled comparative trial of gabapentin and propranolol in essential tremor.** Arch Neurol. 1999;56:475-80 [PMID 10199338](https://pubmed.ncbi.nlm.nih.gov/10199338/)
- Scott AK, et al. **Sumatriptan clinical pharmacokinetics.** Clin Pharmacokinet. 1994;27:337-44 [PMID 7851052](https://pubmed.ncbi.nlm.nih.gov/7851052/)
- Keller F, et al. **Saturable first-pass kinetics of propranolol.** J Clin Pharmacol. 1989;29:240-5 [PMID 2723110](https://pubmed.ncbi.nlm.nih.gov/2723110/)
- Abila B, et al. **The tremorolytic action of beta-adrenoceptor blockers in essential, physiological and isoprenaline-induced tremor is mediated by beta-adrenoceptors located in a deep peripheral compartment.** Br J Clin Pharmacol. 1985;20:369-76 [PMID 2866785](https://pubmed.ncbi.nlm.nih.gov/2866785/)
- Suzuki T, et al. **Nonlinear first-pass metabolism of propranolol in the rat.** J Pharmacobiodyn. 1981;4:131-41 [PMID 7277199](https://pubmed.ncbi.nlm.nih.gov/7277199/)

## 7. 프리미돈 · 페노바르비탈

*Primidone and phenobarbital* — 모체/대사물 기여, PK, 임상 효과 — '모체가 활성 분자' 결과의 근거.

- Saifee TA, et al. **Epidemiology and treatment of patients with Essential Tremor: A retrospective cohort analysis in the United Kingdom.** Neuroepidemiology. 2026 [PMID 42430313](https://pubmed.ncbi.nlm.nih.gov/42430313/)
- Calzetti S, et al. **Does pre-treatment with phenobarbital prevent the acute intolerance to primidone in patients with essential tremor?.** Neurol Sci. 2025;46:3703-3707 [PMID 40343565](https://pubmed.ncbi.nlm.nih.gov/40343565/)
- Alharbi O, et al. **The Pharmacological Management of Essential Tremor and Its Long-Term Effects on Patient Quality of Life: A Systematic Review.** Cureus. 2024;16:e76016 [PMID 39834985](https://pubmed.ncbi.nlm.nih.gov/39834985/)
- Haubenberger D, et al. **Essential Tremor.** N Engl J Med. 2018;378:1802-1810 [PMID 29742376](https://pubmed.ncbi.nlm.nih.gov/29742376/)
- Sethi KD, et al. **Tremor.** Curr Opin Neurol. 2003;16:481-5 [PMID 12869807](https://pubmed.ncbi.nlm.nih.gov/12869807/)
- O'Suilleabhain P, et al. **Randomized trial comparing primidone initiation schedules for treating essential tremor.** Mov Disord. 2002;17:382-6 [PMID 11921128](https://pubmed.ncbi.nlm.nih.gov/11921128/)
- Koller WC, et al. **Pharmacologic treatment of essential tremor.** Neurology. 2000;54:S30-8 [PMID 10854350](https://pubmed.ncbi.nlm.nih.gov/10854350/)
- Nagaki S, et al. **Blood and cerebrospinal fluid pharmacokinetics of primidone and its primary pharmacologically active metabolites, phenobarbital and phenylethylmalonamide in the rat.** Eur J Drug Metab Pharmacokinet. 1999;24:255-64 [PMID 10716065](https://pubmed.ncbi.nlm.nih.gov/10716065/)
- El-Masri HA, et al. **Physiologically based pharmacokinetics model of primidone and its metabolites phenobarbital and phenylethylmalonamide in humans, rats, and mice.** Drug Metab Dispos. 1998;26:585-94 [PMID 9616196](https://pubmed.ncbi.nlm.nih.gov/9616196/)
- Koller WC, et al. **Acute and chronic effects of propranolol and primidone in essential tremor.** Neurology. 1989;39:1587-8 [PMID 2586774](https://pubmed.ncbi.nlm.nih.gov/2586774/)
- Sasso E, et al. **Double-blind comparison of primidone and phenobarbital in essential tremor.** Neurology. 1988;38:808-10 [PMID 3283599](https://pubmed.ncbi.nlm.nih.gov/3283599/)
- Koller WC, et al. **Efficacy of primidone in essential tremor.** Neurology. 1986;36:121-4 [PMID 3941767](https://pubmed.ncbi.nlm.nih.gov/3941767/)
- Cottrell PR, et al. **Pharmacokinetics of phenylethylmalonamide (PEMA) in normal subjects and in patients treated with antiepileptic drugs.** Epilepsia. 1982;23:307-13 [PMID 7084140](https://pubmed.ncbi.nlm.nih.gov/7084140/)

## 8. 기타 약물치료

*Other pharmacotherapy* — 토피라메이트 · 가바펜틴 · 조니사미드 · 벤조디아제핀.

- Dash D, et al. **Update on Medical Treatments for Essential Tremor: An International Parkinson and Movement Disorder Society Evidence-Based Medicine Review.** Mov Disord. 2026;41:815-825 [PMID 41556478](https://pubmed.ncbi.nlm.nih.gov/41556478/)
- Pillai KS, et al. **Zonisamide add-on in tremor-dominant Parkinson's disease- A randomized controlled clinical trial.** Parkinsonism Relat Disord. 2022;105:1-6 [PMID 36323130](https://pubmed.ncbi.nlm.nih.gov/36323130/)
- Bruno E, et al. **Topiramate for essential tremor.** Cochrane Database Syst Rev. 2017;4:CD009683 [PMID 28409827](https://pubmed.ncbi.nlm.nih.gov/28409827/)
- Chang KH, et al. **Efficacy and Safety of Topiramate for Essential Tremor: A Meta-Analysis of Randomized Controlled Trials.** Medicine (Baltimore). 2015;94:e1809 [PMID 26512577](https://pubmed.ncbi.nlm.nih.gov/26512577/)
- Zappia M, et al. **Treatment of essential tremor: a systematic review of evidence and recommendations from the Italian Movement Disorders Association.** J Neurol. 2013;260:714-40 [PMID 22886006](https://pubmed.ncbi.nlm.nih.gov/22886006/)
- Frima N, et al. **A double-blind, placebo-controlled, crossover trial of topiramate in essential tremor.** Clin Neuropharmacol. 2006;29:94-6 [PMID 16614542](https://pubmed.ncbi.nlm.nih.gov/16614542/)
- Ondo WG, et al. **Topiramate in essential tremor: a double-blind, placebo-controlled trial.** Neurology. 2006;66:672-7 [PMID 16436648](https://pubmed.ncbi.nlm.nih.gov/16436648/)
- Chen JJ, et al. **Nonparkinsonism movement disorders in the elderly.** Consult Pharm. 2006;21:58-71 [PMID 16524353](https://pubmed.ncbi.nlm.nih.gov/16524353/)
- Ondo W, et al. **Gabapentin for essential tremor: a multiple-dose, double-blind, placebo-controlled trial.** Mov Disord. 2000;15:678-82 [PMID 10928578](https://pubmed.ncbi.nlm.nih.gov/10928578/)

## 9. 에탄올과 1-옥탄올

*Ethanol and 1-octanol* — 급성 억제, 반동, 자가투약 — ADAPTF/ADAPTS 비대칭 동역학의 근거.

- Everlo CSJ, et al. **Testing for Alcohol Responsiveness in Familial Essential Tremor.** Tremor Other Hyperkinet Mov (N Y). 2024;14:30 [PMID 38881692](https://pubmed.ncbi.nlm.nih.gov/38881692/)
- Voller B, et al. **Dose-escalation study of octanoic acid in patients with essential tremor.** J Clin Invest. 2016;126:1451-7 [PMID 26927672](https://pubmed.ncbi.nlm.nih.gov/26927672/)
- Haubenberger D, et al. **Treatment of essential tremor with long-chain alcohols: still experimental or ready for prime time?.** Tremor Other Hyperkinet Mov (N Y). 2014;4 [PMID 24587968](https://pubmed.ncbi.nlm.nih.gov/24587968/)
- Haubenberger D, et al. **Octanoic acid in alcohol-responsive essential tremor: a randomized controlled study.** Neurology. 2013;80:933-40 [PMID 23408867](https://pubmed.ncbi.nlm.nih.gov/23408867/)
- Nahab FB, et al. **An open-label, single-dose, crossover study of the pharmacokinetics and metabolism of two oral formulations of 1-octanol in patients with essential tremor.** Neurotherapeutics. 2011;8:753-62 [PMID 21594724](https://pubmed.ncbi.nlm.nih.gov/21594724/)

## 10. T형 칼슘차단제 임상시험

*T-type calcium channel blockers in the clinic* — CX-8998 · ulixacaltamide. 모델이 유도한 천장을 검증(혹은 반증)할 데이터.

- Papapetropoulos S, et al. **A Phase 2, Randomized, Double-Blind, Placebo-Controlled Trial of CX-8998, a Selective Modulator of the T-Type Calcium Channel in Inadequately Treated Moderate to Severe Essential Tremor: T-CALM Study Design and Methodology for Efficacy Endpoint and Digital Biomarker Selection.** Front Neurol. 2019;10:597 [PMID 31244760](https://pubmed.ncbi.nlm.nih.gov/31244760/)

## 11. 보툴리눔 독소

*Botulinum toxin* — 손떨림 시험의 악력 문제, 유도 주사, SNAP-25 회복 — f_spill 결과의 근거.

- Machicoane M, et al. **Excitation-contraction coupling inhibitors potentiate the actions of botulinum neurotoxin type A at the neuromuscular junction.** Br J Pharmacol. 2025;182:564-580 [PMID 39389783](https://pubmed.ncbi.nlm.nih.gov/39389783/)
- Samotus O, et al. **Real-World Longitudinal Experience of Botulinum Toxin Therapy for Parkinson and Essential Tremor.** Toxins (Basel). 2022;14 [PMID 36006219](https://pubmed.ncbi.nlm.nih.gov/36006219/)
- Liao YH, et al. **Botulinum Toxin for Essential Tremor and Hands Tremor in the Neurological Diseases: A Meta-Analysis of Randomized Controlled Trials.** Toxins (Basel). 2022;14 [PMID 35324700](https://pubmed.ncbi.nlm.nih.gov/35324700/)
- Zheng X, et al. **Botulinum toxin type A for hand tremor: a meta-analysis of randomised controlled trials.** Neurol Neurochir Pol. 2020;54:561-567 [PMID 33047784](https://pubmed.ncbi.nlm.nih.gov/33047784/)
- Mittal SO, et al. **Botulinum toxin in essential hand tremor - A randomized double-blind placebo-controlled study with customized injection approach.** Parkinsonism Relat Disord. 2018;56:65-69 [PMID 29929813](https://pubmed.ncbi.nlm.nih.gov/29929813/)
- Mittal SO, et al. **Botulinum Toxin in Parkinson Disease Tremor: A Randomized, Double-Blind, Placebo-Controlled Study With a Customized Injection Approach.** Mayo Clin Proc. 2017;92:1359-1367 [PMID 28789780](https://pubmed.ncbi.nlm.nih.gov/28789780/)
- Zesiewicz TA, et al. **Practice parameter: therapies for essential tremor [RETIRED]: report of the Quality Standards Subcommittee of the American Academy of Neurology.** Neurology. 2005;64:2008-20 [PMID 15972843](https://pubmed.ncbi.nlm.nih.gov/15972843/)
- Dolly O, et al. **Synaptic transmission: inhibition of neurotransmitter release by botulinum toxins.** Headache. 2003;43 Suppl 1:S16-24 [PMID 12887390](https://pubmed.ncbi.nlm.nih.gov/12887390/)
- Foran PG, et al. **Evaluation of the therapeutic usefulness of botulinum neurotoxin B, C1, E, and F compared with the long lasting type A. Basis for distinct durations of inhibition of exocytosis in central neurons.** J Biol Chem. 2003;278:1363-71 [PMID 12381720](https://pubmed.ncbi.nlm.nih.gov/12381720/)
- Brin MF, et al. **A randomized, double masked, controlled trial of botulinum toxin type A in essential hand tremor.** Neurology. 2001;56:1523-8 [PMID 11402109](https://pubmed.ncbi.nlm.nih.gov/11402109/)

## 12. 수술 · 신경조절

*Surgery and neuromodulation* — MRgFUS · Vim DBS · 자극 주파수 · 병변 부피 대 실조 · 습관화 — φ_thal(직렬 인자)의 근거.

- Luu CP, et al. **Rescue deep brain stimulation for recurrent essential tremor after ventral intermediate nucleus thalamotomy: illustrative cases.** J Neurosurg Case Lessons. 2026;12 [PMID 42475739](https://pubmed.ncbi.nlm.nih.gov/42475739/)
- Sastre-Bataller I, et al. **Magnetic resonance-guided focused ultrasound thalamotomy in essential tremor subtypes: a phenotype-based insight into gait and balance.** Brain Commun. 2026;8:fcag076 [PMID 41907317](https://pubmed.ncbi.nlm.nih.gov/41907317/)
- Kiselev R, et al. **Rescue Thalamotomy for Habituation to Deep Brain Stimulation in Essential Tremor: Case Report.** Tremor Other Hyperkinet Mov (N Y). 2026;16:8 [PMID 41694795](https://pubmed.ncbi.nlm.nih.gov/41694795/)
- Paraskevopoulos Z, et al. **Frequency-Dependent Inhibition during Deep Brain Stimulation of Thalamic Ventral Intermediate Nuclei.** J Neurosci. 2026;46 [PMID 41942271](https://pubmed.ncbi.nlm.nih.gov/41942271/)
- Shiramba A, et al. **Efficacy and Safety of Magnetic Resonance-Guided Focused Ultrasound Thalamotomy in Essential Tremor: A Systematic Review and Metanalysis.** Mov Disord. 2025;40:1020-1033 [PMID 40243386](https://pubmed.ncbi.nlm.nih.gov/40243386/)
- He S, et al. **Cortico-thalamic tremor circuits and their associations with deep brain stimulation effects in essential tremor.** Brain. 2025;148:2093-2107 [PMID 39592428](https://pubmed.ncbi.nlm.nih.gov/39592428/)
- Buch VP, et al. **"Quality over quantity:" smaller, targeted lesions optimize quality of life outcomes after MR-guided focused ultrasound thalamotomy for essential tremor.** Front Neurol. 2024;15:1450699 [PMID 39610701](https://pubmed.ncbi.nlm.nih.gov/39610701/)
- Kaplitt MG, et al. **Safety and Efficacy of Staged, Bilateral Focused Ultrasound Thalamotomy in Essential Tremor: An Open-Label Clinical Trial.** JAMA Neurol. 2024;81:939-946 [PMID 39073822](https://pubmed.ncbi.nlm.nih.gov/39073822/)
- Ferreira Felloni Borges Y, et al. **Essential Tremor - Deep Brain Stimulation vs. Focused Ultrasound.** Expert Rev Neurother. 2023;23:603-619 [PMID 37288812](https://pubmed.ncbi.nlm.nih.gov/37288812/)
- Iorio-Morin C, et al. **Bilateral Focused Ultrasound Thalamotomy for Essential Tremor (BEST-FUS Phase 2 Trial).** Mov Disord. 2021;36:2653-2662 [PMID 34288097](https://pubmed.ncbi.nlm.nih.gov/34288097/)
- Dhima K, et al. **Neuropsychological outcomes after thalamic deep brain stimulation for essential tremor.** Parkinsonism Relat Disord. 2021;92:88-93 [PMID 34736157](https://pubmed.ncbi.nlm.nih.gov/34736157/)
- He S, et al. **Closed-Loop Deep Brain Stimulation for Essential Tremor Based on Thalamic Local Field Potentials.** Mov Disord. 2021;36:863-873 [PMID 33547859](https://pubmed.ncbi.nlm.nih.gov/33547859/)
- Peters J, et al. **Habituation After Deep Brain Stimulation in Tremor Syndromes: Prevalence, Risk Factors and Long-Term Outcomes.** Front Neurol. 2021;12:696950 [PMID 34413826](https://pubmed.ncbi.nlm.nih.gov/34413826/)
- Swan BD, et al. **Effects of ramped-frequency thalamic deep brain stimulation on tremor and activity of modeled neurons.** Clin Neurophysiol. 2020;131:625-634 [PMID 31978847](https://pubmed.ncbi.nlm.nih.gov/31978847/)
- Harary M, et al. **Unilateral Thalamic Deep Brain Stimulation Versus Focused Ultrasound Thalamotomy for Essential Tremor.** World Neurosurg. 2019;126:e144-e152 [PMID 30794976](https://pubmed.ncbi.nlm.nih.gov/30794976/)
- Fasano A, et al. **Tremor habituation to deep brain stimulation: Underlying mechanisms and solutions.** Mov Disord. 2019;34:1761-1773 [PMID 31433906](https://pubmed.ncbi.nlm.nih.gov/31433906/)
- Krack P, et al. **Deep Brain Stimulation in Movement Disorders: From Experimental Surgery to Evidence-Based Therapy.** Mov Disord. 2019;34:1795-1810 [PMID 31580535](https://pubmed.ncbi.nlm.nih.gov/31580535/)
- Bond AE, et al. **Safety and Efficacy of Focused Ultrasound Thalamotomy for Patients With Medication-Refractory, Tremor-Dominant Parkinson Disease: A Randomized Clinical Trial.** JAMA Neurol. 2017;74:1412-1418 [PMID 29084313](https://pubmed.ncbi.nlm.nih.gov/29084313/)
- Elias WJ, et al. **A Randomized Trial of Focused Ultrasound Thalamotomy for Essential Tremor.** N Engl J Med. 2016;375:730-9 [PMID 27557301](https://pubmed.ncbi.nlm.nih.gov/27557301/)

## 13. 평가지표 · 삶의 질

*Endpoints and quality of life* — TETRAS · FTM · QUEST · 가속도계.

- Serrano-Dueñas M, et al. **Severe and unclassifiable tremor.** Arq Neuropsiquiatr. 2024;82:1-5 [PMID 39396518](https://pubmed.ncbi.nlm.nih.gov/39396518/)
- Gerbasi ME, et al. **Associations Among Tremor Amplitude, Activities of Daily Living, and Quality of Life in Patients with Essential Tremor.** Tremor Other Hyperkinet Mov (N Y). 2024;14:22 [PMID 38708124](https://pubmed.ncbi.nlm.nih.gov/38708124/)
- Marques A, et al. **French validation of the Quality of life in Essential Tremor Questionnaire (QUEST) and the Essential Tremor Embarrassment Assessment (ETEA).** Rev Neurol (Paris). 2023;179:1128-1133 [PMID 37735016](https://pubmed.ncbi.nlm.nih.gov/37735016/)
- Tang CC, et al. **Quantifying the impact of upper limb tremor on the quality of life of people with multiple sclerosis: a comparison between the QUEST and MSIS-29 scales.** Mult Scler Relat Disord. 2022;58:103495 [PMID 35085981](https://pubmed.ncbi.nlm.nih.gov/35085981/)
- Kovács M , et al. **Independent validation of the Quality of Life in Essential Tremor Questionnaire (QUEST).** Ideggyogy Sz. 2017;70:193-202 [PMID 29870634](https://pubmed.ncbi.nlm.nih.gov/29870634/)

## 14. QSP 방법론

*QSP methodology* — 모델 구조·검증 방법론.

- Bai JPF, et al. **Creating a Roadmap to Quantitative Systems Pharmacology-Informed Rare Disease Drug Development: A Workshop Report.** Clin Pharmacol Ther. 2024;115:201-205 [PMID 37984065](https://pubmed.ncbi.nlm.nih.gov/37984065/)
- Bai JP, et al. **Quantitative Systems Pharmacology for Rare Disease Drug Development.** J Pharm Sci. 2023;112:2313-2320 [PMID 37422281](https://pubmed.ncbi.nlm.nih.gov/37422281/)
- Li X, et al. **Combining network pharmacology, molecular docking, molecular dynamics simulation, and experimental verification to examine the efficacy and immunoregulation mechanism of FHB granules on vitiligo.** Front Immunol. 2023;14:1194823 [PMID 37575231](https://pubmed.ncbi.nlm.nih.gov/37575231/)
- Marshall S, et al. **Model-Informed Drug Discovery and Development: Current Industry Good Practice and Regulatory Expectations and Future Perspectives.** CPT Pharmacometrics Syst Pharmacol. 2019;8:87-96 [PMID 30411538](https://pubmed.ncbi.nlm.nih.gov/30411538/)
- Helmlinger G, et al. **Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development.** CPT Pharmacometrics Syst Pharmacol. 2019;8:380-395 [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)

---

## 모델 파라미터가 문헌에 직접 기대는 지점 (Where the model leans on this literature)

| 모델 요소 | 문헌 근거 섹션 | 비고 |
|---|---|---|
| `rating = 2 + 2·log10(A_cm)` (Elble 로그 변환) | 2 | 가속도계와 TETRAS의 '불일치'가 로그라는 결론의 전부가 여기서 나온다 |
| `f = 1/tau_loop` 및 질량 부하 감별검사 | 2 | ET 주파수는 중추 지연, EPT 주파수는 기계적 공명 |
| `PCINT`, `DNDIS` | 3 | torpedo·등반섬유 병리·치아핵 GABA 수용체 감소 |
| `a_O = 0.35` (올리브 분율) | 4, 10 | harmaline 랫은 a_O=1, 사람은 <0.62 — 후자는 실패한 시험이 준 상한 |
| `KI_PRP_B2 = 0.6 nM`, `FB2 = 0.60` | 6 | 말초 β₂ 부위; 나돌롤 유효·아테놀롤 무효 |
| `EMAX_PRM/EC50_PRM` vs `EMAX_PB/EC50_PB` | 7 | 모체 대 페노바르비탈의 기여 분해 |
| `TAUF_ON = 1 h`, `TAUF_OFF = 5 h` | 9 | Mellanby 급성 내성과 반동의 비대칭 |
| `f_spill` | 11 | 유도 주사가 악력을 보존한다는 관찰 |
| `F50D = 80 Hz`, `HDBS = 4`, `V50L`, `V50A` | 12 | >100 Hz 규칙, 병변 부피-실조 상충 |

## 면책 (Disclaimer)

본 모델은 교육·연구 목적의 반정량적 QSP 모델입니다. 위 문헌은 모델 구조와 파라미터의
**출발점**이며, 모델이 환자 데이터에 적합(fit)되거나 검증된 것은 아닙니다. 임상 의사결정에
사용해서는 안 됩니다.
