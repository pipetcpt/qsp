# 다계통 위축 (Multiple System Atrophy, MSA) — 참고문헌

> **PMID 검증 (PMID verification).** 아래 모든 항목은 NCBI E-utilities
> (`esearch` + `esummary`)로 프로그램적으로 조회하여 **PMID · 제1저자 ·
> 저널 · 연도 · 제목이 실제 PubMed 레코드와 일치하는지 확인**한 것입니다.
> 검증되지 않은 인용은 이 파일에 포함하지 않았습니다. 한 건의 예외는
> §12의 verdiperstat 3상(M-STAR)으로, 이는 PubMed 등재 논문이 아니라
> 임상시험 등록번호로 표기했습니다 — 아래에 그 사실을 명시했습니다.
>
> 링크 형식: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`

---

## 1. 진단기준 및 임상 개요 (Diagnostic Criteria & Clinical Overview)

MSA 모델의 임상 엔드포인트(§22 클러스터)와 표현형 분류(MSA-P / MSA-C)는
2022년 MDS 기준을 따릅니다. 인지장애가 더 이상 배제기준이 아니라는 점,
그리고 "prodromal MSA" 범주가 새로 도입된 점이 모델의 조기 개입 시나리오
(시나리오 19, 항-α-synuclein 항체 조기 vs 후기 투여)의 근거입니다.

1. Wenning GK, et al. **The Movement Disorder Society Criteria for the Diagnosis of Multiple System Atrophy.** *Mov Disord* 2022. [PMID 35445419](https://pubmed.ncbi.nlm.nih.gov/35445419/)
2. Gilman S, et al. **Second consensus statement on the diagnosis of multiple system atrophy.** *Neurology* 2008. [PMID 18725592](https://pubmed.ncbi.nlm.nih.gov/18725592/)
3. Sekiya H, et al. **Validation Study of the MDS Criteria for the Diagnosis of Multiple System Atrophy in the Mayo Clinic Brain Bank.** *Neurology* 2023. [PMID 37816641](https://pubmed.ncbi.nlm.nih.gov/37816641/)
4. Sun Y, et al. **Comparison of the second consensus statement with the movement disorder society criteria for multiple system atrophy.** *Parkinsonism Relat Disord* 2023. [PMID 36529110](https://pubmed.ncbi.nlm.nih.gov/36529110/)
5. Lamotte G, et al. **Movement disorder society criteria for the diagnosis of multiple system atrophy — what's new?** *Clin Auton Res* 2022. [PMID 35633428](https://pubmed.ncbi.nlm.nih.gov/35633428/)
6. Stankovic I, et al. **A Review on the Clinical Diagnosis of Multiple System Atrophy.** *Cerebellum* 2023. [PMID 35986227](https://pubmed.ncbi.nlm.nih.gov/35986227/)
7. Wenning GK, et al. **Multiple system atrophy.** *Handb Clin Neurol* 2013. [PMID 24095129](https://pubmed.ncbi.nlm.nih.gov/24095129/)
8. Peeraully T. **Multiple system atrophy.** *Semin Neurol* 2014. [PMID 24963676](https://pubmed.ncbi.nlm.nih.gov/24963676/)
9. Wenning GK, et al. **Multiple system atrophy: a review of 203 pathologically proven cases.** *Mov Disord* 1997. [PMID 9087971](https://pubmed.ncbi.nlm.nih.gov/9087971/)
10. Lin DJ, et al. **The Diagnosis and Natural History of Multiple System Atrophy, Cerebellar Type.** *Cerebellum* 2016. [PMID 26467153](https://pubmed.ncbi.nlm.nih.gov/26467153/)
11. Liu M, et al. **Multiple system atrophy: an update and emerging directions of biomarkers and clinical trials.** *J Neurol* 2024. [PMID 38483626](https://pubmed.ncbi.nlm.nih.gov/38483626/)

---

## 2. 자연경과 · 진행속도 · 생존 (Natural History, Progression Rate & Survival)

**모델 보정 근거 (calibration targets).** mrgsolve 모델의 `KND`
(신경퇴행 구동 계수)와 생존 위험함수(`H0`, `B_U2`, `B_STR`, `B_RESP`)는
아래 코호트에서 보고된 두 수치에 맞추어 보정했습니다:
**UMSARS-II 연간 5–8점 증가**, **발병 후 중앙 생존 6–10년**.
모델 자체 검정(self-test)은 4→8년 구간 기울기와 생존 중앙값을
각각 3–11점/년, 6–11년 범위로 강제합니다.

12. Wenning GK, et al. **The natural history of multiple system atrophy: a prospective European cohort study.** *Lancet Neurol* 2013. [PMID 23391524](https://pubmed.ncbi.nlm.nih.gov/23391524/)
13. Low PA, et al. **Natural history of multiple system atrophy in the USA: a prospective cohort study.** *Lancet Neurol* 2015. [PMID 26025783](https://pubmed.ncbi.nlm.nih.gov/26025783/)
14. Goldstein DS, et al. **Survival in synucleinopathies: A prospective cohort study.** *Neurology* 2015. [PMID 26432848](https://pubmed.ncbi.nlm.nih.gov/26432848/)
15. Wenning GK, et al. **Development and validation of the Unified Multiple System Atrophy Rating Scale (UMSARS).** *Mov Disord* 2004. [PMID 15452868](https://pubmed.ncbi.nlm.nih.gov/15452868/)
16. Palma JA, et al. **Limitations of the Unified Multiple System Atrophy Rating Scale as outcome measure for clinical trials and a roadmap for improvement.** *Clin Auton Res* 2021. [PMID 33554315](https://pubmed.ncbi.nlm.nih.gov/33554315/)
17. Geser F, et al. **The European Multiple System Atrophy-Study Group (EMSA-SG).** *J Neural Transm (Vienna)* 2005. [PMID 16049636](https://pubmed.ncbi.nlm.nih.gov/16049636/)
18. Xiao Y, et al. **Modified version of unified multiple system atrophy rating scale for remote video-based assessments.** *NPJ Parkinsons Dis* 2023. [PMID 37891215](https://pubmed.ncbi.nlm.nih.gov/37891215/)
19. Kaufmann H, et al. **Multiple System Atrophy Combined Outcome Assessment (MuSyCA): process, format, and validation plan.** *Clin Auton Res* 2026. [PMID 41762390](https://pubmed.ncbi.nlm.nih.gov/41762390/)
20. Feng T, et al. **Natural history and 12-month progression of multiple system atrophy in a Chinese cohort.** *BMC Neurol* 2026. [PMID 42286509](https://pubmed.ncbi.nlm.nih.gov/42286509/)

---

## 3. 올리고덴드로글리아 α-시누클레인 병리 (Oligodendroglial α-Synuclein / GCI)

**모델의 네 번째 구조적 주장의 근거.** MSA는 뉴런이 아니라
**올리고덴드로사이트**에서 α-synuclein이 응집하는 질환입니다. 모델의
`ASYNM → ASYNO → GCI` 캐스케이드는 오직 `P25A` (p25α/TPPP의 수초 →
세포체 재분포)에 의해 촉매되며, 그 `P25A`는 다시 `MYE`(수초 온전성)의
감소에 의해 생성됩니다 — 즉 자기증폭 루프입니다. 뉴런 사멸은
`WTROPH`(잃어버린 올리고덴드로글리아 영양지지) 항을 통해 **하류에서**
일어나므로, 뉴런만 표적하는 치료는 원리적으로 잘못된 구획을 겨냥합니다.

21. Wakabayashi K, et al. **Alpha-synuclein immunoreactivity in glial cytoplasmic inclusions in multiple system atrophy.** *Neurosci Lett* 1998. [PMID 9682846](https://pubmed.ncbi.nlm.nih.gov/9682846/)
22. Ndayisaba A, et al. **Multiple System Atrophy: Pathology, Pathogenesis, and Path Forward.** *Annu Rev Pathol* 2025. [PMID 39405585](https://pubmed.ncbi.nlm.nih.gov/39405585/)
23. Reddy K, et al. **Multiple system atrophy: α-Synuclein strains at the neuron-oligodendrocyte crossroad.** *Mol Neurodegener* 2022. [PMID 36435784](https://pubmed.ncbi.nlm.nih.gov/36435784/)
24. Hoffmann A, et al. **Oligodendroglial α-synucleinopathy-driven neuroinflammation in multiple system atrophy.** *Brain Pathol* 2019. [PMID 30444295](https://pubmed.ncbi.nlm.nih.gov/30444295/)
25. Wiseman JA, et al. **Neuronal α-synuclein toxicity is the key driver of neurodegeneration in multiple system atrophy.** *Brain* 2025. [PMID 39908177](https://pubmed.ncbi.nlm.nih.gov/39908177/)
26. Mavroeidi P, et al. **Endogenous oligodendroglial alpha-synuclein and TPPP/p25α orchestrate alpha-synuclein pathology in experimental multiple system atrophy models.** *Acta Neuropathol* 2019. [PMID 31011860](https://pubmed.ncbi.nlm.nih.gov/31011860/)
27. Lindersson E, et al. **p25alpha Stimulates alpha-synuclein aggregation and is co-localized with aggregated alpha-synuclein in alpha-synucleinopathies.** *J Biol Chem* 2005. [PMID 15590652](https://pubmed.ncbi.nlm.nih.gov/15590652/)
28. Ferreira N, et al. **Multiple system atrophy-associated oligodendroglial protein p25α stimulates formation of novel α-synuclein strain with enhanced neurodegenerative potential.** *Acta Neuropathol* 2021. [PMID 33978813](https://pubmed.ncbi.nlm.nih.gov/33978813/)
29. Mavroeidi P, et al. **Autophagy mediates the clearance of oligodendroglial SNCA/alpha-synuclein and TPPP/p25A in multiple system atrophy models.** *Autophagy* 2022. [PMID 35000546](https://pubmed.ncbi.nlm.nih.gov/35000546/)
30. Kragh CL, et al. **FAS-dependent cell death in α-synuclein transgenic oligodendrocyte models of multiple system atrophy.** *PLoS One* 2013. [PMID 23372841](https://pubmed.ncbi.nlm.nih.gov/23372841/)

---

## 4. 선택적 신경세포 소실 — 왜 절전(preganglionic)이 무너지는가 (Selective Neuronal Loss)

**모델의 첫 번째 구조적 주장의 해부학적 근거.** MSA에서 사라지는 것은
척수 중간외측핵(IML)의 **절전** 교감신경세포와 뇌간 심혈관 뉴런이며,
**절후(postganglionic)** 종말은 상대적으로 보존됩니다. 모델은 이 비대칭을
`VULN_PG = 0.08` (다른 집단의 1/10 이하) 단 하나의 파라미터로 표현하고,
`G_CENT`(중추 이득)와 `POSTG`(절후 온전성)를 곱셈 항으로 분리합니다.
같은 방정식에서 `VULN_PG`만 크게 올리면 순수 자율신경실패(PAF) 표현형이
되고, 그것만으로 atomoxetine의 반응이 사라집니다(시나리오 7/8).

31. Oppenheimer DR. **Lateral horn cells in progressive autonomic failure.** *J Neurol Sci* 1980. [PMID 6247458](https://pubmed.ncbi.nlm.nih.gov/6247458/)
32. Terao S, et al. **Disease-specific patterns of neuronal loss in the spinal ventral horn in amyotrophic lateral sclerosis, multiple system atrophy and Werdnig-Hoffmann disease.** *J Neurol* 1994. [PMID 8195817](https://pubmed.ncbi.nlm.nih.gov/8195817/)
33. Cortelli P, et al. **Autonomic blood pressure control.** *Handb Clin Neurol* 2026. [PMID 41896018](https://pubmed.ncbi.nlm.nih.gov/41896018/)

---

## 5. 프리온 유사 전파 (Prion-like Seeding & Strain Biology)

모델의 `SEED` 구획(세포외 seeding-competent α-synuclein)과
`SEED → ASYN_UPTAKE` 양성 피드백은 아래 전파 실험들에 근거합니다. 이
구획은 항-α-synuclein 항체(`MAB`)가 결합할 수 있는 **유일한** 표적이며,
이 때문에 모델에서 항체의 효과는 아직 남은 뉴런 수에 의해 게이팅됩니다.

34. Prusiner SB, et al. **Evidence for α-synuclein prions causing multiple system atrophy in humans with parkinsonism.** *Proc Natl Acad Sci U S A* 2015. [PMID 26324905](https://pubmed.ncbi.nlm.nih.gov/26324905/)
35. Holec SAM, et al. **Multiple system atrophy prions transmit neurological disease to mice expressing wild-type human α-synuclein.** *Acta Neuropathol* 2022. [PMID 36018376](https://pubmed.ncbi.nlm.nih.gov/36018376/)
36. Holec SAM, et al. **α-synuclein prion strains differentially adapt after passage in mice.** *PLoS Pathog* 2024. [PMID 39642110](https://pubmed.ncbi.nlm.nih.gov/39642110/)
37. Dhillon JS, et al. **Dissecting α-synuclein inclusion pathology diversity in multiple system atrophy: implications for the prion-like transmission hypothesis.** *Lab Invest* 2019. [PMID 30737468](https://pubmed.ncbi.nlm.nih.gov/30737468/)
38. Jellinger KA, et al. **Is Multiple System Atrophy a Prion-like Disorder?** *Int J Mol Sci* 2021. [PMID 34576255](https://pubmed.ncbi.nlm.nih.gov/34576255/)

---

## 6. 유전 · CoQ10 (Genetics & the Coenzyme Q10 Axis)

모델 파라미터 `COQ2F`(CoQ10 생합성능)는 COQ2 변이 보유자를 표현하며,
ubiquinol 시나리오(시나리오 17)에서 야생형과 COQ2 결핍형의 반응 차이가
**코딩된 것이 아니라** `CQ → OXS → NDRIVE` 경로에서 유도되어 나옵니다.

39. Multiple-System Atrophy Research Collaboration. **Mutations in COQ2 in familial and sporadic multiple-system atrophy.** *N Engl J Med* 2013. [PMID 23758206](https://pubmed.ncbi.nlm.nih.gov/23758206/)
40. Ogaki K, et al. **Analysis of COQ2 gene in multiple system atrophy.** *Mol Neurodegener* 2014. [PMID 25373618](https://pubmed.ncbi.nlm.nih.gov/25373618/)
41. Porto KJ, et al. **COQ2 V393A confers high risk susceptibility for multiple system atrophy in East Asian population.** *J Neurol Sci* 2021. [PMID 34455210](https://pubmed.ncbi.nlm.nih.gov/34455210/)
42. Procopio R, et al. **Genetic mutation analysis of the COQ2 gene in Italian patients with multiple system atrophy.** *Gene* 2019. [PMID 31398377](https://pubmed.ncbi.nlm.nih.gov/31398377/)
43. Scholz SW, et al. **SNCA variants are associated with increased risk for multiple system atrophy.** *Ann Neurol* 2009. [PMID 19475667](https://pubmed.ncbi.nlm.nih.gov/19475667/)
44. Federoff M, et al. **Multiple system atrophy: the application of genetics in understanding etiology.** *Clin Auton Res* 2015. [PMID 25687905](https://pubmed.ncbi.nlm.nih.gov/25687905/)
45. Bougea A, et al. **Genetics of Multiple System Atrophy and Progressive Supranuclear Palsy: A Systemized Review of the Literature.** *Int J Mol Sci* 2023. [PMID 36982356](https://pubmed.ncbi.nlm.nih.gov/36982356/)
46. Li XY, et al. **Genetic profiles of multiple system atrophy revealed by exome sequencing, long-read sequencing and spinocerebellar ataxia repeat expansion analysis.** *Eur J Neurol* 2024. [PMID 39152783](https://pubmed.ncbi.nlm.nih.gov/39152783/)

---

## 7. 신경염증 · 미엘로퍼옥시다제 (Neuroinflammation & Myeloperoxidase)

**시나리오 18이 왜 null이어야 하는지의 근거.** MSA 형질전환 마우스에서
MPO 억제는 조기에 시작하면 신경보호를 보였지만(#47), **지연 개시**에서는
미세아교세포를 확실히 억제했음에도 신경보호에 실패했습니다(#48). 모델은
이 구분을 `WINDOW` 논리로 재현합니다 — 즉 MPO 억제는 `MGL → MPO → OXS`
경로만 차단하고, 진단 시점에는 이미 `GCI`와 `1-MYE` 항이 `NDRIVE`를
지배하고 있으므로 UMSARS-II 곡선은 거의 움직이지 않습니다.

47. Stefanova N, et al. **Myeloperoxidase inhibition ameliorates multiple system atrophy-like degeneration in a transgenic mouse model.** *Neurotox Res* 2012. [PMID 22161470](https://pubmed.ncbi.nlm.nih.gov/22161470/)
48. Kaindlstorfer C, et al. **Failure of Neuroprotection Despite Microglial Suppression by Delayed-Start Myeloperoxidase Inhibition in a Model of Advanced Multiple System Atrophy.** *Neurotox Res* 2015. [PMID 26194617](https://pubmed.ncbi.nlm.nih.gov/26194617/)

---

## 8. 기립성 저혈압 · 앙와위 고혈압 · 야간 다뇨 (nOH, Supine Hypertension & the Nocturnal Loop)

**모델의 세 번째 구조적 주장의 근거이며, 이 파일에서 가장 중요한 절.**
Goldstein(#49)은 자율신경실패에서 **앙와위 고혈압과 기립성 저혈압이 서로
연관되어 함께 나타난다**는 것을 보였고, Okamoto(#50)는 야간 혈압
dipping의 소실을, Shibao(#51)는 앙와위 고혈압이 **압력 나트륨배설
(pressure natriuresis)** 을 통해 야간 용적 손실을 일으킨다는 것을
보였습니다. Mathias(#54)는 그 루프를 닫는 실험을 했습니다 — 취침 시
desmopressin이 야간 다뇨와 밤사이 체중감소를 줄이고 **다음 날 아침의
기립성 저혈압을 개선**했습니다.

모델에서 `UNAV`와 `UVOL`은 모두 `exp(K*(MAP - MAPNAT))` 항을 가지므로,
누워 있는 밤 동안 상승한 혈압이 그날 아침에 필요한 용적을 미리 써버립니다.
`KPN`(압력 항)은 25 mmHg 상승에 대해 약 1.5–2배의 야간 나트륨배설이
되도록 **의도적으로 완만하게** 두고, 장기 용적 항상성은 별도의 급한
`KVOLN`(용적 오차) 항이 담당합니다. 이 분리가 없으면 루프가 비현실적으로
폭주하거나(급한 압력항) 총 나트륨이 발산합니다(압력항만 있는 경우).

49. Goldstein DS, et al. **Association between supine hypertension and orthostatic hypotension in autonomic failure.** *Hypertension* 2003. [PMID 12835329](https://pubmed.ncbi.nlm.nih.gov/12835329/)
50. Okamoto LE, et al. **Nocturnal blood pressure dipping in the hypertension of autonomic failure.** *Hypertension* 2009. [PMID 19047577](https://pubmed.ncbi.nlm.nih.gov/19047577/)
51. Shibao C, et al. **Clonidine for the treatment of supine hypertension and pressure natriuresis in autonomic failure.** *Hypertension* 2006. [PMID 16391172](https://pubmed.ncbi.nlm.nih.gov/16391172/)
52. Fanciulli A, et al. **Consensus statement on the definition of neurogenic supine hypertension in cardiovascular autonomic failure by the American Autonomic Society and the European Federation of Autonomic Societies.** *Clin Auton Res* 2018. [PMID 29766366](https://pubmed.ncbi.nlm.nih.gov/29766366/)
53. Park JW, et al. **Advances in the Pathophysiology and Management of Supine Hypertension in Patients with Neurogenic Orthostatic Hypotension.** *Curr Hypertens Rep* 2022. [PMID 35230654](https://pubmed.ncbi.nlm.nih.gov/35230654/)
54. Mathias CJ, et al. **The effect of desmopressin on nocturnal polyuria, overnight weight loss, and morning postural hypotension in patients with autonomic failure.** *Br Med J (Clin Res Ed)* 1986. [PMID 3089519](https://pubmed.ncbi.nlm.nih.gov/3089519/)
55. Umbertini E, et al. **Understanding nocturnal polyuria in cardiovascular autonomic failure: Pathophysiological mechanisms and clinical implications.** *Auton Neurosci* 2026. [PMID 42241932](https://pubmed.ncbi.nlm.nih.gov/42241932/)
56. Norcliffe-Kaufmann L, et al. **Orthostatic heart rate changes in patients with autonomic failure caused by neurodegenerative synucleinopathies.** *Ann Neurol* 2018. [PMID 29405350](https://pubmed.ncbi.nlm.nih.gov/29405350/)
57. Pavy-Le Traon A, et al. **New insights into orthostatic hypotension in multiple system atrophy: a European multicentre cohort study.** *J Neurol Neurosurg Psychiatry* 2016. [PMID 25977316](https://pubmed.ncbi.nlm.nih.gov/25977316/)
58. Jiang Q, et al. **Orthostatic Hypotension in Multiple System Atrophy: Related Factors and Disease Prognosis.** *J Parkinsons Dis* 2023. [PMID 38143372](https://pubmed.ncbi.nlm.nih.gov/38143372/)
59. Idiaquez JF, et al. **Neurogenic Orthostatic Hypotension. Lessons From Synucleinopathies.** *Am J Hypertens* 2021. [PMID 33705537](https://pubmed.ncbi.nlm.nih.gov/33705537/)

---

## 9. 신경호르몬 반응의 해리 (Dissociated Neurohormonal Responses: NE, Renin, AVP)

모델이 재현해야 하는 MSA의 신경호르몬 지문은 세 가지입니다:
**(i) 앙와위 혈장 NE는 정상**(절후 종말이 살아 있으므로),
**(ii) 기립 시 NE 증가폭은 둔화**(중추 구동이 없으므로),
**(iii) 레닌은 기립에도 오르지 않는다**(사구체옆장치의 β1 교감신경 자극이
사라졌으므로). 모델에서 (i)–(iii)은 각각 `POSTG`, `G_CENT`,
`KBRENIN*SNA*POSTG` 항에서 자동으로 따라옵니다. AVP의 압수용체 팔
(`KBAVP*G_CENT*...`)이 `G_CENT`로 게이팅되어 있으므로, 저혈압이 더 이상
AVP를 방출시키지 못하는 것도 같은 구조에서 나옵니다.

60. Biaggioni I, et al. **Hyporeninemic normoaldosteronism in severe autonomic failure.** *J Clin Endocrinol Metab* 1993. [PMID 7680352](https://pubmed.ncbi.nlm.nih.gov/7680352/)
61. Giza RJ, et al. **Clinical and neurohormonal characteristics in African Americans with neurogenic orthostatic hypotension.** *Clin Auton Res* 2021. [PMID 33502643](https://pubmed.ncbi.nlm.nih.gov/33502643/)
62. Mendoza-Velásquez JJ, et al. **Autonomic Dysfunction in α-Synucleinopathies.** *Front Neurol* 2019. [PMID 31031694](https://pubmed.ncbi.nlm.nih.gov/31031694/)

---

## 10. 심장 교감신경 영상 및 피부 생검 — MSA와 PD/PAF의 갈림길 (Cardiac Sympathetic Imaging & Skin Biopsy)

`POSTG`가 MSA에서 보존되고 PD/PAF에서 소실된다는 모델의 핵심 가정은
**임상적으로 직접 측정 가능**합니다: ¹²³I-MIBG 심근 섭취가 MSA에서는
정상이고 PD/PAF에서는 감소합니다. 모델은 이를 `MIBG` 노드로 표현하고,
`POSTG`를 통해 atomoxetine 반응과 **같은 상태변수**에 연결합니다 — 즉
모델에서 MIBG는 예측 바이오마커가 됩니다.

63. King AE, et al. **Meta-analysis of 123I-MIBG cardiac scintigraphy for the diagnosis of Lewy body-related disorders.** *Mov Disord* 2011. [PMID 21480373](https://pubmed.ncbi.nlm.nih.gov/21480373/)
64. Alves Do Rego C, et al. **Prospective study of relevance of (123)I-MIBG myocardial scintigraphy and clonidine GH test to distinguish Parkinson's disease and multiple system atrophy.** *J Neurol* 2018. [PMID 29956027](https://pubmed.ncbi.nlm.nih.gov/29956027/)
65. Catalan M, et al. **(123)I-Metaiodobenzylguanidine Myocardial Scintigraphy in Discriminating Degenerative Parkinsonisms.** *Mov Disord Clin Pract* 2021. [PMID 34295947](https://pubmed.ncbi.nlm.nih.gov/34295947/)
66. Yang T, et al. **(131)I-MIBG myocardial scintigraphy for differentiation of Parkinson's disease from multiple system atrophy or essential tremor.** *J Neurol Sci* 2017. [PMID 28131225](https://pubmed.ncbi.nlm.nih.gov/28131225/)
67. Donadio V, et al. **Skin sympathetic fiber α-synuclein deposits: a potential biomarker for pure autonomic failure.** *Neurology* 2013. [PMID 23390175](https://pubmed.ncbi.nlm.nih.gov/23390175/)

---

## 11. 레보도파 무반응의 후시냅스 기원 (Postsynaptic Origin of Levodopa Failure)

**모델의 두 번째 구조적 주장의 결정적 근거.** Churchyard(#68)는 MSA의
도파 저항성이 **후시냅스 D2 수용체 소실**이라는 것을 직접 보였고,
Sawle(#69)은 전/후시냅스 변화를 영상으로 분리했습니다. 모델은
선조체 출력을 `STRIAT = G_POST × DA/(EC50DA + DA)`로 쓰고, MSA에서는
`G_POST = NMSN^GEXP_MSN`이 감소하도록 합니다. 레보도파는 `DA`만 올릴 수
있으므로, **같은 뇌 노출에서 PD는 거의 정상화되고 MSA-P는 점점 반응이
줄어듭니다** — 초기 부분반응 후 1–2년 내 소실, 그리고 이상운동증(dyskinesia)
부재라는 두 red flag가 모두 이 한 줄에서 나옵니다.

68. Churchyard A, et al. **Dopa resistance in multiple-system atrophy: loss of postsynaptic D2 receptors.** *Ann Neurol* 1993. [PMID 8338346](https://pubmed.ncbi.nlm.nih.gov/8338346/)
69. Sawle GV, et al. **Asymmetrical pre-synaptic and post-synaptic changes in the striatal dopamine projection in dopa naïve parkinsonism.** *Brain* 1993. [PMID 8353712](https://pubmed.ncbi.nlm.nih.gov/8353712/)
70. Booij J, et al. **The clinical benefit of imaging striatal dopamine transporters with [123I]FP-CIT SPET in differentiating patients with presynaptic parkinsonism from those with other forms of parkinsonism.** *Eur J Nucl Med* 2001. [PMID 11315592](https://pubmed.ncbi.nlm.nih.gov/11315592/)

---

## 12. 승압 약물 — 병소 위치별 선택성 (Pressor Pharmacology: Selectivity by Lesion Site)

### 12.1 미도드린 (midodrine) — 모든 병소의 하류에서 작용
모델에서 활성대사체 desglymidodrine은 α1 작용제 항 `AGON`에 **직접**
더해지므로, 중추 병소든 절후 병소든 무관하게 작동합니다. 또한
`A1R`(탈신경 과민성)이 이를 곱하므로, **같은 mg이 병이 진행할수록 더 큰
ΔSBP를 만든다**는 임상 관찰이 모델에서 유도됩니다(자체검정 항목).

71. Low PA, et al. **Efficacy of midodrine vs placebo in neurogenic orthostatic hypotension. A randomized, double-blind multicenter study.** *JAMA* 1997. [PMID 9091692](https://pubmed.ncbi.nlm.nih.gov/9091692/)
72. Jankovic J, et al. **Neurogenic orthostatic hypotension: a double-blind, placebo-controlled study with midodrine.** *Am J Med* 1993. [PMID 7687093](https://pubmed.ncbi.nlm.nih.gov/7687093/)
73. Wright RA, et al. **A double-blind, dose-response study of midodrine in neurogenic orthostatic hypotension.** *Neurology* 1998. [PMID 9674789](https://pubmed.ncbi.nlm.nih.gov/9674789/)
74. Fouad-Tarazi FM, et al. **Alpha sympathomimetic treatment of autonomic insufficiency with orthostatic hypotension.** *Am J Med* 1995. [PMID 7503082](https://pubmed.ncbi.nlm.nih.gov/7503082/)

### 12.2 드록시도파 (droxidopa) — 절후 AADC를 **필요로** 한다
모델: `DROXNE = KAADC × CDRX × POSTG × (1 − CARBI)`. 이 한 줄이
(a) 드록시도파가 종말이 살아 있어야 작동한다는 것과 (b) **파킨슨증
때문에 함께 투여하는 카르비도파가 같은 효소를 점유해 서로 길항한다**는
것을 동시에 만들어 냅니다(시나리오 6).

75. Kaufmann H, et al. **Droxidopa for neurogenic orthostatic hypotension: a randomized, placebo-controlled, phase 3 trial.** *Neurology* 2014. [PMID 24944260](https://pubmed.ncbi.nlm.nih.gov/24944260/)
76. Biaggioni I, et al. **Randomized withdrawal study of patients with symptomatic neurogenic orthostatic hypotension responsive to droxidopa.** *Hypertension* 2015. [PMID 25350981](https://pubmed.ncbi.nlm.nih.gov/25350981/)
77. Elgebaly A, et al. **Meta-analysis of the safety and efficacy of droxidopa for neurogenic orthostatic hypotension.** *Clin Auton Res* 2016. [PMID 26951135](https://pubmed.ncbi.nlm.nih.gov/26951135/)
78. Strassheim V, et al. **Droxidopa for orthostatic hypotension: a systematic review and meta-analysis.** *J Hypertens* 2016. [PMID 27442791](https://pubmed.ncbi.nlm.nih.gov/27442791/)
79. Isaacson S, et al. **Long-term safety of droxidopa in patients with symptomatic neurogenic orthostatic hypotension.** *J Am Soc Hypertens* 2016. [PMID 27614923](https://pubmed.ncbi.nlm.nih.gov/27614923/)
80. Chen JJ, et al. **Standing and Supine Blood Pressure Outcomes Associated With Droxidopa and Midodrine in Patients With Neurogenic Orthostatic Hypotension.** *Ann Pharmacother* 2018. [PMID 29972032](https://pubmed.ncbi.nlm.nih.gov/29972032/)

### 12.3 NET 억제제 (atomoxetine, ampreloxetine) — MSA 선택성의 실험적 증거
NET 차단은 **이미 일어나고 있는** NE 방출을 증폭시키므로, 중추 병소 +
온전한 절후 종말(=MSA)에서 크고 절후 병소(=PAF)에서 작습니다. 모델에서
이 선택성은 `NEREL = SNA × NEVES × POSTG`에 NET 차단이 곱해지는 구조
그 자체이며, 어디에도 "atomoxetine은 MSA에 선택적"이라고 쓰여 있지
않습니다(시나리오 7/8, 자체검정에서 강제).

81. Byun JI, et al. **Efficacy of atomoxetine versus midodrine for neurogenic orthostatic hypotension.** *Ann Clin Transl Neurol* 2020. [PMID 31856425](https://pubmed.ncbi.nlm.nih.gov/31856425/)
82. Mwesigwa N, et al. **Atomoxetine on neurogenic orthostatic hypotension: a randomized, double-blind, placebo-controlled crossover trial.** *Clin Auton Res* 2024. [PMID 39294522](https://pubmed.ncbi.nlm.nih.gov/39294522/)
83. Okamoto LE, et al. **Synergistic effect of norepinephrine transporter blockade and α-2 antagonism on blood pressure in autonomic failure.** *Hypertension* 2012. [PMID 22311903](https://pubmed.ncbi.nlm.nih.gov/22311903/)
84. Lo A, et al. **Pharmacokinetics and pharmacodynamics of ampreloxetine, a novel, selective norepinephrine reuptake inhibitor, in symptomatic neurogenic orthostatic hypotension.** *Clin Auton Res* 2021. [PMID 33782836](https://pubmed.ncbi.nlm.nih.gov/33782836/)
85. Kaufmann H, et al. **Safety and efficacy of ampreloxetine in symptomatic neurogenic orthostatic hypotension: a phase 2 trial.** *Clin Auton Res* 2021. [PMID 34657222](https://pubmed.ncbi.nlm.nih.gov/34657222/)
86. Hoxhaj P, et al. **Ampreloxetine Versus Droxidopa in Neurogenic Orthostatic Hypotension: A Comparative Review.** *Cureus* 2023. [PMID 37303338](https://pubmed.ncbi.nlm.nih.gov/37303338/)

### 12.4 α2 길항제 (yohimbine) — 잔여 중추 구동을 필요로 한다
87. Onrot J, et al. **Oral yohimbine in human autonomic failure.** *Neurology* 1987. [PMID 3808301](https://pubmed.ncbi.nlm.nih.gov/3808301/)
88. Biaggioni I, et al. **Manipulation of norepinephrine metabolism with yohimbine in the treatment of autonomic failure.** *J Clin Pharmacol* 1994. [PMID 8089252](https://pubmed.ncbi.nlm.nih.gov/8089252/)

### 12.5 피리도스티그민 (pyridostigmine) — 반사성 신경전달만 증폭
신경절 AChE 억제는 **압수용체가 구동하는** 절전 신경전달을 증폭하므로,
누워 있을 때(반사 오차 ≈ 0)에는 거의 작용하지 않습니다. 모델은
`PYRAMP`을 `SNA_ss`의 **압수용체 항에만** 곱하여 이를 표현하고, 그 결과
"기립 혈압은 올리면서 앙와위 고혈압은 상대적으로 덜 악화시킨다"는
관찰이 유도됩니다(시나리오 13, 자체검정에서 강제).

89. Singer W, et al. **Pyridostigmine treatment trial in neurogenic orthostatic hypotension.** *Arch Neurol* 2006. [PMID 16476804](https://pubmed.ncbi.nlm.nih.gov/16476804/)
90. Okamoto LE, et al. **Clinical Correlates of Efficacy of Pyridostigmine in the Treatment of Orthostatic Hypotension.** *Hypertension* 2025. [PMID 39727053](https://pubmed.ncbi.nlm.nih.gov/39727053/)
91. Holder AC, et al. **Pyridostigmine for the Management of Neurogenic Orthostatic Hypotension: A Systemic Review.** *J Geriatr Psychiatry Neurol* 2025. [PMID 39043171](https://pubmed.ncbi.nlm.nih.gov/39043171/)

### 12.6 용적 확장 · 비약물 요법 (Volume Expansion & Non-Pharmacological Measures)
92. Veazie S, et al. **Fludrocortisone for orthostatic hypotension.** *Cochrane Database Syst Rev* 2021. [PMID 34000076](https://pubmed.ncbi.nlm.nih.gov/34000076/)
93. van Lieshout JJ, et al. **Fludrocortisone and sleeping in the head-up position limit the postural decrease in cardiac output in autonomic failure.** *Clin Auton Res* 2000. [PMID 10750642](https://pubmed.ncbi.nlm.nih.gov/10750642/)
94. May M, et al. **The osmopressor response to water drinking.** *Am J Physiol Regul Integr Comp Physiol* 2011. [PMID 21048076](https://pubmed.ncbi.nlm.nih.gov/21048076/)
95. Okamoto LE, et al. **Efficacy of Servo-Controlled Splanchnic Venous Compression in the Treatment of Orthostatic Hypotension: A Randomized Comparison With Midodrine.** *Hypertension* 2016. [PMID 27271310](https://pubmed.ncbi.nlm.nih.gov/27271310/)
96. van der Stam AH, et al. **The Impact of Head-Up Tilt Sleeping on Orthostatic Tolerance: A Scoping Review.** *Biology (Basel)* 2023. [PMID 37626994](https://pubmed.ncbi.nlm.nih.gov/37626994/)
97. van der Stam AH, et al. **Tolerability and efficacy of full-body head-up tilt sleeping in Parkinson's disease and multiple system atrophy.** *NPJ Parkinsons Dis* 2026. [PMID 42143029](https://pubmed.ncbi.nlm.nih.gov/42143029/)
98. van der Stam AH, et al. **Study protocol for the Heads-Up trial: a phase II randomized controlled trial investigating head-up tilt sleeping to alleviate orthostatic intolerance.** *BMC Neurol* 2024. [PMID 38166676](https://pubmed.ncbi.nlm.nih.gov/38166676/)

### 12.7 식후 저혈압과 옥트레오티드 (Post-prandial Hypotension & Octreotide)
99. Armstrong E, et al. **The effects of the somatostatin analogue, octreotide, on postural hypotension, before and after food ingestion, in primary autonomic failure.** *Clin Auton Res* 1991. [PMID 1822761](https://pubmed.ncbi.nlm.nih.gov/1822761/)
100. Alam M, et al. **Effects of the peptide release inhibitor, octreotide, on daytime hypotension and on nocturnal hypertension in primary autonomic failure.** *J Hypertens* 1995. [PMID 8903629](https://pubmed.ncbi.nlm.nih.gov/8903629/)
101. Smith GD, et al. **Effect of the somatostatin analogue, octreotide, on exercise-induced hypotension in human subjects with chronic sympathetic failure.** *Clin Sci (Lond)* 1995. [PMID 7493436](https://pubmed.ncbi.nlm.nih.gov/7493436/)
102. Jansen RW, et al. **Postprandial hypotension: epidemiology, pathophysiology, and clinical management.** *Ann Intern Med* 1995. [PMID 7825766](https://pubmed.ncbi.nlm.nih.gov/7825766/)
103. Chaudhuri KR, et al. **Alcohol ingestion lowers supine blood pressure, causes splanchnic vasodilatation and worsens postural hypotension in primary autonomic failure.** *J Neurol* 1994. [PMID 8164016](https://pubmed.ncbi.nlm.nih.gov/8164016/)

### 12.8 치료 지침 · 종합 리뷰 (Guidelines & Comprehensive Reviews)
104. Park JW, et al. **Pharmacologic treatment of orthostatic hypotension.** *Auton Neurosci* 2020. [PMID 32979782](https://pubmed.ncbi.nlm.nih.gov/32979782/)
105. Eschlböck S, et al. **Evidence-based treatment of neurogenic orthostatic hypotension and related symptoms.** *J Neural Transm (Vienna)* 2017. [PMID 29058089](https://pubmed.ncbi.nlm.nih.gov/29058089/)
106. Palma JA, et al. **Management of Orthostatic Hypotension.** *Continuum (Minneap Minn)* 2020. [PMID 31996627](https://pubmed.ncbi.nlm.nih.gov/31996627/)
107. Shibao C, et al. **Pharmacotherapy of autonomic failure.** *Pharmacol Ther* 2012. [PMID 21664375](https://pubmed.ncbi.nlm.nih.gov/21664375/)
108. Chen B, et al. **Non-pharmacological and drug treatment of autonomic dysfunction in multiple system atrophy: current status and future directions.** *J Neurol* 2023. [PMID 37477834](https://pubmed.ncbi.nlm.nih.gov/37477834/)
109. Arbique D, et al. **Management of neurogenic orthostatic hypotension.** *J Am Med Dir Assoc* 2014. [PMID 24388946](https://pubmed.ncbi.nlm.nih.gov/24388946/)
110. Vidal-Petiot E, et al. **Orthostatic hypotension: Review and expert position statement.** *Rev Neurol (Paris)* 2024. [PMID 38123372](https://pubmed.ncbi.nlm.nih.gov/38123372/)
111. Fanciulli A, et al. **Management of Orthostatic Hypotension in Parkinson's Disease.** *J Parkinsons Dis* 2020. [PMID 32716319](https://pubmed.ncbi.nlm.nih.gov/32716319/)

---

## 13. 질병조절 치료 시도 (Disease-Modifying Trials — mostly negative)

**⚠️ verdiperstat 3상 (M-STAR) 주의.** 이 시험(NCT03952806, Biohaven)의
1차 종료점 미달 결과는 **PubMed에 등재된 논문으로 확인하지 못했습니다.**
따라서 이 파일은 그 결과를 학술 인용으로 제시하지 않고 등록번호로만
표기하며, 모델의 시나리오 18(verdiperstat null)의 문헌적 근거는 §7의
전임상 지연개시 실험(#48)과 아래 rifampicin 실패(#112)에 둡니다.
시나리오 18은 "3상 결과의 재현"이 아니라 **"진단 시점에 개시된 하류
염증 표적 단일 억제는 이 모델 구조에서 null이 되어야 한다"** 는 구조적
예측으로 읽어야 합니다.

112. Low PA, et al. **Efficacy and safety of rifampicin for multiple system atrophy: a randomised, double-blind, placebo-controlled trial.** *Lancet Neurol* 2014. [PMID 24507091](https://pubmed.ncbi.nlm.nih.gov/24507091/)
113. Singer W, et al. **Optimizing clinical trial design for multiple system atrophy: lessons from the rifampicin study.** *Clin Auton Res* 2015. [PMID 25763826](https://pubmed.ncbi.nlm.nih.gov/25763826/)
114. Mitsui J, et al. **High-dose ubiquinol supplementation in multiple-system atrophy: a multicentre, randomised, double-blinded, placebo-controlled phase 2 trial.** *EClinicalMedicine* 2023. [PMID 37256098](https://pubmed.ncbi.nlm.nih.gov/37256098/)
115. Mitsui J, et al. **Three-Year Follow-Up of High-Dose Ubiquinol Supplementation in a Case of Familial Multiple System Atrophy with Compound Heterozygous COQ2 Mutations.** *Cerebellum* 2017. [PMID 28150130](https://pubmed.ncbi.nlm.nih.gov/28150130/)
116. Singer W, et al. **Intrathecal administration of autologous mesenchymal stem cells in multiple system atrophy.** *Neurology* 2019. [PMID 31152011](https://pubmed.ncbi.nlm.nih.gov/31152011/)
117. Poewe W, et al. **Therapeutic advances in multiple system atrophy and progressive supranuclear palsy.** *Mov Disord* 2015. [PMID 26227071](https://pubmed.ncbi.nlm.nih.gov/26227071/)
118. Eschlböck S, et al. **Interventional trials in atypical parkinsonism.** *Parkinsonism Relat Disord* 2016. [PMID 26421389](https://pubmed.ncbi.nlm.nih.gov/26421389/)
119. Jeong SH, et al. **Drug repurposing for disease-modifying effects in multiple system atrophy.** *Transl Neurodegener* 2026. [PMID 42010648](https://pubmed.ncbi.nlm.nih.gov/42010648/)

---

## 14. 바이오마커 (Biomarkers: NfL, Seed Amplification, MRI)

모델의 `NFL` 구획은 신경 소실 **속도**에 비례하여 방출되고 1차 소실
되므로, 자연경과에서 중기(약 6년)에 최대가 되고 병이 소진되면서 서서히
내려옵니다. 이는 NfL이 **단면 중증도보다 진행 속도의 지표**라는 관찰
(#120–#123)과 부합하며, 동시에 모델이 말기에 NfL을 과대추정하지 않게
합니다.

120. Chelban V, et al. **Neurofilament light levels predict clinical progression and death in multiple system atrophy.** *Brain* 2022. [PMID 35903017](https://pubmed.ncbi.nlm.nih.gov/35903017/)
121. Zhang L, et al. **Neurofilament Light Chain Predicts Disease Severity and Progression in Multiple System Atrophy.** *Mov Disord* 2022. [PMID 34719813](https://pubmed.ncbi.nlm.nih.gov/34719813/)
122. Singer W, et al. **Neurofilament light chain in spinal fluid and plasma in multiple system atrophy: a prospective, longitudinal biomarker study.** *Clin Auton Res* 2023. [PMID 37603107](https://pubmed.ncbi.nlm.nih.gov/37603107/)
123. Huang J, et al. **Plasma Neurofilament Light Chain as a Biomarker for Motor Progression and Disease Milestones in Multiple System Atrophy: An Update.** *Mov Disord* 2026. [PMID 41912379](https://pubmed.ncbi.nlm.nih.gov/41912379/)
124. Jin B, et al. **Plasma neurofilament light chain as a predictor of multiple system atrophy in idiopathic REM sleep behavior disorder.** *J Neurol* 2025. [PMID 41359208](https://pubmed.ncbi.nlm.nih.gov/41359208/)
125. Fernandes Gomes B, et al. **α-Synuclein seed amplification assay as a diagnostic tool for parkinsonian disorders.** *Parkinsonism Relat Disord* 2023. [PMID 37591709](https://pubmed.ncbi.nlm.nih.gov/37591709/)
126. Rossi M, et al. **Comparison of Two α-Synuclein Seed Amplification Assays for Discrimination of Parkinson Disease and Atypical Parkinsonism.** *Mov Disord* 2025. [PMID 40879244](https://pubmed.ncbi.nlm.nih.gov/40879244/)
127. Grossauer A, et al. **α-Synuclein Seed Amplification Assays in the Diagnosis of Synucleinopathies Using Cerebrospinal Fluid — A Systematic Review and Meta-Analysis.** *Mov Disord Clin Pract* 2023. [PMID 37205253](https://pubmed.ncbi.nlm.nih.gov/37205253/)
128. Sugiyama A, et al. **Revisiting 'hot cross bun' sign: a multicentre MRI study of 97 patients with autopsy-confirmed multiple system atrophy.** *J Neurol Neurosurg Psychiatry* 2026. [PMID 41083252](https://pubmed.ncbi.nlm.nih.gov/41083252/)
129. Portet M, et al. **Hot cross bun sign.** *J Neurol* 2019. [PMID 31254063](https://pubmed.ncbi.nlm.nih.gov/31254063/)

---

## 15. 수면 · 협착음 · 돌연사 (Sleep, Stridor & Sudden Death)

모델의 생존 위험함수는 `B_STR`(협착음)과 `B_RESP`(호흡 뉴런 소실)에
독립적인 가중치를 둡니다. 이는 협착음이 UMSARS 운동점수와 **별개로**
생존을 예측하고, 기관절개가 생존을 연장한다는 관찰(#131–#133)을
반영한 것입니다.

130. Silber MH, et al. **Stridor and death in multiple system atrophy.** *Mov Disord* 2000. [PMID 10928581](https://pubmed.ncbi.nlm.nih.gov/10928581/)
131. Cortelli P, et al. **Stridor in multiple system atrophy: Consensus statement on diagnosis, prognosis, and treatment.** *Neurology* 2019. [PMID 31570638](https://pubmed.ncbi.nlm.nih.gov/31570638/)
132. Giannini G, et al. **Early stridor onset and stridor treatment predict survival in 136 patients with MSA.** *Neurology* 2016. [PMID 27566741](https://pubmed.ncbi.nlm.nih.gov/27566741/)
133. Giannini G, et al. **Tracheostomy is associated with increased survival in Multiple System Atrophy patients with stridor.** *Eur J Neurol* 2022. [PMID 35384153](https://pubmed.ncbi.nlm.nih.gov/35384153/)
134. Laga A, et al. **A strategic approach of the management of sleep-disordered breathing in multiple system atrophy.** *J Clin Sleep Med* 2025. [PMID 39539061](https://pubmed.ncbi.nlm.nih.gov/39539061/)
135. Giannini G, et al. **REM Sleep Behaviour Disorder in Multiple System Atrophy: From Prodromal to Progression of Disease.** *Front Neurol* 2021. [PMID 34194385](https://pubmed.ncbi.nlm.nih.gov/34194385/)
136. Postuma RB, et al. **Evolution of Prodromal Multiple System Atrophy from REM Sleep Behavior Disorder: A Descriptive Study.** *J Parkinsons Dis* 2022. [PMID 35094998](https://pubmed.ncbi.nlm.nih.gov/35094998/)
137. Postuma RB, et al. **Risk and predictors of dementia and parkinsonism in idiopathic REM sleep behaviour disorder: a multicentre study.** *Brain* 2019. [PMID 30789229](https://pubmed.ncbi.nlm.nih.gov/30789229/)

---

## 16. 비뇨생식 자율신경 장애 (Urogenital Autonomic Failure — Onuf's Nucleus)

모델에서 `NONUF`(Onuf 핵)는 **가장 취약한** 집단(`VULN_ONUF = 1.05`)으로
설정되어, 잔뇨(PVR)와 발기부전이 운동증상보다 먼저 나타납니다. 이는
방광 기능장애가 MSA의 **첫 증상**일 수 있다는 전향적 관찰(#139)과,
남성에서 발기부전이 가장 이른 증상이라는 보고(#138)를 반영합니다.

138. Papatsoris AG, et al. **Urinary and erectile dysfunction in multiple system atrophy (MSA).** *Neurourol Urodyn* 2008. [PMID 17563111](https://pubmed.ncbi.nlm.nih.gov/17563111/)
139. Sakakibara R, et al. **Bladder dysfunction as the initial presentation of multiple system atrophy: a prospective cohort study.** *Clin Auton Res* 2019. [PMID 30043182](https://pubmed.ncbi.nlm.nih.gov/30043182/)

---

## 17. 모델링 도구 (Modelling Tools)

- **mrgsolve** — ODE 기반 PK/PD·QSP 시뮬레이션: <https://mrgsolve.org/>
- **Graphviz** — 기계론적 지도 렌더링: <https://graphviz.org/>
- **Shiny** — 인터랙티브 대시보드: <https://shiny.posit.co/>
- gPKPDviz (mrgsolve 기반 Shiny 시뮬레이션 도구) 논문:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>

---

## ⚠️ 면책 조항 (Disclaimer)

본 참고문헌 목록과 이에 연결된 QSP 모델은 **교육 및 연구 목적**으로
작성되었습니다. 파라미터는 위 문헌에 보고된 값의 범위에 맞추어 보정한
**설명용 근사치**이며, 개별 환자 데이터에 대한 적합·검증을 거치지
않았습니다. **실제 임상 의사결정, 처방, 규제 제출에 직접 사용해서는 안
됩니다.**
