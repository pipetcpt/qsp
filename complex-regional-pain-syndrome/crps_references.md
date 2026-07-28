# CRPS QSP 모델 — 참고문헌 (References)

복합부위통증증후군(Complex Regional Pain Syndrome, CRPS)의 기계론적 지도,
mrgsolve ODE 모델, Shiny 대시보드 구축에 사용한 문헌 목록입니다.

> **링크 형식에 관한 註**: 아래 링크는 모두 **PubMed 검색 질의(query) URL**입니다.
> PMID를 잘못 기재하면 전혀 다른 논문을 인용하게 되므로, 저자·연도·제목 키워드로
> 검색되는 링크를 사용했습니다. 클릭하면 해당 논문(또는 그 논문이 최상위인 결과
> 목록)으로 이동합니다. 모델 파라미터가 특정 문헌에 정박(anchor)된 부분은
> `crps_mrgsolve_model.R` 상단의 CALIBRATION ANCHORS 절에 대응 관계가 정리되어
> 있습니다.

---

## 1. 진단 기준 · 중증도 척도 (Diagnostic criteria & severity)

1. Harden RN, et al. **Validation of proposed diagnostic criteria (the "Budapest Criteria") for complex regional pain syndrome.** *Pain.* 2010.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Harden+validation+proposed+diagnostic+criteria+Budapest+complex+regional+pain+syndrome+Pain+2010>
2. Harden RN, et al. **A comparison of the CRPS Severity Score (CSS) and diagnostic criteria over time.** *Pain.* 2017.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Harden+CRPS+Severity+Score+responsiveness+Pain+2017>
3. Harden RN, et al. **Complex Regional Pain Syndrome: Practical Diagnostic and Treatment Guidelines, 5th Edition.** *Pain Med.* 2022.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Harden+complex+regional+pain+syndrome+practical+diagnostic+treatment+guidelines+5th+edition>
4. Goebel A, et al. **The Valencia consensus-based adaptation of the IASP complex regional pain syndrome diagnostic criteria.** *Pain.* 2021.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Goebel+Valencia+consensus+adaptation+IASP+complex+regional+pain+syndrome+criteria>

## 2. 역학 · 자연경과 (Epidemiology & natural history)

5. de Mos M, et al. **The incidence of complex regional pain syndrome: a population-based study.** *Pain.* 2007. (26.2/100,000 인-년)
   <https://pubmed.ncbi.nlm.nih.gov/?term=de+Mos+incidence+complex+regional+pain+syndrome+population-based+study+Pain+2007>
6. Sandroni P, et al. **Complex regional pain syndrome type I: incidence and prevalence in Olmsted county.** *Pain.* 2003. (5.46/100,000)
   <https://pubmed.ncbi.nlm.nih.gov/?term=Sandroni+complex+regional+pain+syndrome+incidence+prevalence+Olmsted+county+Pain+2003>
7. Beerthuizen A, et al. **Demographic and medical parameters in the development of CRPS type 1 after fractures.** *Pain.* 2012. (골절 후 발생률)
   <https://pubmed.ncbi.nlm.nih.gov/?term=Beerthuizen+demographic+medical+parameters+development+complex+regional+pain+syndrome+type+1+fractures>
8. Bean DJ, et al. **The outcome of complex regional pain syndrome type 1: a systematic review.** *J Pain.* 2014.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Bean+outcome+complex+regional+pain+syndrome+type+1+systematic+review+J+Pain+2014>
9. Schwartzman RJ, et al. **The natural history of complex regional pain syndrome.** *Clin J Pain.* 2009. (확산·진행)
   <https://pubmed.ncbi.nlm.nih.gov/?term=Schwartzman+natural+history+complex+regional+pain+syndrome+Clin+J+Pain+2009>
10. Bruehl S. **Complex regional pain syndrome.** *BMJ.* 2015. (임상 개요)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bruehl+complex+regional+pain+syndrome+BMJ+2015+review>

## 3. 신경성 염증 · 사이토카인 (Neurogenic inflammation & cytokines)

11. Birklein F, et al. **Neuropeptides and neurogenic inflammation in CRPS.** *Neurology / Neurosci Lett.* 2001-2011 계열.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Birklein+neuropeptide+neurogenic+inflammation+complex+regional+pain+syndrome>
12. Huygen FJPM, et al. **Evidence for local inflammation in complex regional pain syndrome type 1.** *Mediators Inflamm.* 2002. (수포액 IL-6·TNF-α)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Huygen+evidence+local+inflammation+complex+regional+pain+syndrome+blister+fluid+mediators+inflammation>
13. Parkitny L, et al. **Inflammation in complex regional pain syndrome: a systematic review and meta-analysis.** *Neurology.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Parkitny+inflammation+complex+regional+pain+syndrome+systematic+review+meta-analysis+Neurology+2013>
14. Birklein F, Schmelz M. **Neuropeptides, neurogenic inflammation and complex regional pain syndrome.** *Neurosci Lett.* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Birklein+Schmelz+neuropeptides+neurogenic+inflammation+complex+regional+pain+syndrome+Neurosci+Lett>
15. Birklein F, et al. **Activation of cutaneous immune responses in complex regional pain syndrome (keratinocyte proliferation).** *J Pain.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/?term=activation+cutaneous+immune+responses+complex+regional+pain+syndrome+keratinocyte>

## 4. 소섬유 병변 · 말초 신경 (Small-fibre pathology)

16. Oaklander AL, et al. **Evidence of focal small-fiber axonal degeneration in complex regional pain syndrome-I.** *Pain.* 2006. (IENFD 감소)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Oaklander+focal+small+fiber+axonal+degeneration+complex+regional+pain+syndrome+Pain+2006>
17. Albrecht PJ, et al. **Pathologic alterations of cutaneous innervation and vasculature in CRPS.** *Pain.* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Albrecht+pathologic+alterations+cutaneous+innervation+vasculature+complex+regional+pain+syndrome>

## 5. 자가면역 기전 (Autoimmunity)

18. Kohr D, et al. **Autoimmunity against the β2 adrenergic receptor and muscarinic-2 receptor in complex regional pain syndrome.** *Pain.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kohr+autoimmunity+beta2+adrenergic+receptor+muscarinic+2+complex+regional+pain+syndrome+Pain+2011>
19. Tékus V, et al. **A CRPS-IgG-transfer-trauma model reproducing inflammatory and positive sensory signs.** *Pain.* 2014. (수동 이전)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Tekus+CRPS+IgG+transfer+trauma+model+reproducing+inflammatory+sensory+signs+Pain+2014>
20. Helyes Z, et al. **Transfer of complex regional pain syndrome to mice via human autoantibodies.** *PNAS.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Helyes+transfer+complex+regional+pain+syndrome+mice+human+autoantibodies+PNAS+2019>
21. Goebel A, et al. **Intravenous immunoglobulin treatment of the complex regional pain syndrome: a randomized trial.** *Ann Intern Med.* 2010. (소규모 양성)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Goebel+intravenous+immunoglobulin+complex+regional+pain+syndrome+randomized+trial+Annals+2010>
22. Goebel A, et al. **Low-dose intravenous immunoglobulin treatment for long-standing complex regional pain syndrome (LIPS): a randomized trial.** *Ann Intern Med.* 2017. (**음성** — 모델의 자가항체 경로 가중치 제약)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Goebel+low+dose+intravenous+immunoglobulin+long+standing+complex+regional+pain+syndrome+LIPS+randomized>

## 6. 산화 스트레스 · 자유라디칼 소거제 (Oxidative stress & scavengers)

23. Perez RSGM, et al. **The treatment of complex regional pain syndrome type I with free radical scavengers: a randomized controlled study (DMSO vs N-acetylcysteine).** *J Pain Symptom Manage / Pain.* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Perez+treatment+complex+regional+pain+syndrome+free+radical+scavengers+DMSO+N-acetylcysteine+randomized>
24. Zollinger PE, et al. **Can vitamin C prevent complex regional pain syndrome in patients with wrist fractures?** *J Bone Joint Surg Am.* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zollinger+vitamin+C+prevent+complex+regional+pain+syndrome+wrist+fractures+randomized>
25. Eisenberg E, et al. **Evidence for oxidative stress / lipid peroxidation in CRPS.** *Pain Med / Eur J Pain.* 2008-2012 계열.
    <https://pubmed.ncbi.nlm.nih.gov/?term=oxidative+stress+lipid+peroxidation+complex+regional+pain+syndrome>

## 7. 혈관 기능이상 · 온난/한랭 표현형 (Vascular dysfunction, warm vs cold)

26. Groeneweg JG, et al. **Increased endothelin-1 and diminished nitric oxide levels in blister fluids of patients with intermediate cold type CRPS type 1.** *BMC Musculoskelet Disord.* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Groeneweg+endothelin+1+nitric+oxide+blister+fluid+cold+type+complex+regional+pain+syndrome>
27. Groeneweg G, et al. **Effect of tadalafil on blood flow, pain, and function in chronic cold complex regional pain syndrome: a randomized controlled trial.** *BMC Musculoskelet Disord.* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Groeneweg+tadalafil+blood+flow+pain+function+chronic+cold+complex+regional+pain+syndrome+randomized>
28. Wasner G, et al. **Vascular abnormalities in reflex sympathetic dystrophy (CRPS I): mechanisms and diagnostic value.** *Brain.* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Wasner+vascular+abnormalities+reflex+sympathetic+dystrophy+CRPS+mechanisms+diagnostic+Brain>
29. Eberle T, et al. **Warm and cold complex regional pain syndromes: differences beyond skin temperature?** *Neurology.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Eberle+warm+cold+complex+regional+pain+syndromes+differences+beyond+skin+temperature+Neurology>
30. Koban M, et al. **Tissue hypoxia in complex regional pain syndrome.** *Pain.* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Koban+tissue+hypoxia+complex+regional+pain+syndrome+Pain+2003>

## 8. 교감신경-구심신경 결합 (Sympatho-afferent coupling)

31. Baron R, et al. **Causalgia and reflex sympathetic dystrophy: does the sympathetic nervous system contribute to the generation of pain?** *Muscle Nerve / Lancet Neurol.* 계열.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Baron+sympathetic+nervous+system+contribute+generation+pain+causalgia+reflex+sympathetic+dystrophy>
32. Drummond PD, et al. **α1-adrenoceptor expression in hyperalgesic skin of complex regional pain syndrome.** *Pain / Clin J Pain.* 계열.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Drummond+alpha1+adrenoceptor+expression+hyperalgesic+skin+complex+regional+pain+syndrome>
33. O'Connell NE, et al. **Local anaesthetic sympathetic blockade for complex regional pain syndrome (Cochrane review).** 2016.
    <https://pubmed.ncbi.nlm.nih.gov/?term=O%27Connell+local+anaesthetic+sympathetic+blockade+complex+regional+pain+syndrome+Cochrane>

## 9. 중추 감작 · 교세포 (Central sensitisation & glia)

34. Ji RR, et al. **Glia and pain: is chronic pain a gliopathy?** *Pain.* 2013. (교세포 잠금 구조의 근거)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ji+glia+and+pain+is+chronic+pain+a+gliopathy+Pain+2013>
35. Ji RR, et al. **Astrocytes in chronic pain and itch.** *Nat Rev Neurosci.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ji+astrocytes+chronic+pain+itch+Nature+Reviews+Neuroscience>
36. Coull JAM, et al. **BDNF from microglia causes the shift in neuronal anion gradient underlying neuropathic pain.** *Nature.* 2005. (BDNF-KCC2 탈억제)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Coull+BDNF+microglia+shift+neuronal+anion+gradient+neuropathic+pain+Nature+2005>
37. Latremoliere A, Woolf CJ. **Central sensitization: a generator of pain hypersensitivity by central neural plasticity.** *J Pain.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Latremoliere+Woolf+central+sensitization+generator+pain+hypersensitivity+central+neural+plasticity>
38. Seifert F, et al. **Differential endogenous pain modulation in complex regional pain syndrome.** *Brain.* 2009. (하행 조절 손실)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Seifert+differential+endogenous+pain+modulation+complex+regional+pain+syndrome+Brain+2009>

## 10. NMDA 차단 · 케타민 (NMDA antagonism)

39. Sigtermans MJ, et al. **Ketamine produces effective and long-term pain relief in patients with complex regional pain syndrome type 1.** *Pain.* 2009. (100시간 주입, 11주 후 소실)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sigtermans+ketamine+effective+long+term+pain+relief+complex+regional+pain+syndrome+type+1+Pain+2009>
40. Schwartzman RJ, et al. **Outpatient intravenous ketamine for the treatment of complex regional pain syndrome: a double-blind placebo controlled study.** *Pain.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Schwartzman+outpatient+intravenous+ketamine+complex+regional+pain+syndrome+double+blind+placebo>
41. Dahan A, et al. **Population pharmacokinetic-pharmacodynamic modeling of ketamine-induced pain relief of chronic pain.** *Eur J Pain / PLoS One.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dahan+population+pharmacokinetic+pharmacodynamic+ketamine+induced+pain+relief+chronic+pain>
42. Noppers I, et al. **Ketamine for the treatment of chronic non-cancer pain (review).** *Expert Opin Pharmacother.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Noppers+ketamine+treatment+chronic+non+cancer+pain+expert+opinion>

## 11. 골 재형성 · 비스포스포네이트 (Bone remodelling & bisphosphonates)

43. Varenna M, et al. **Treatment of complex regional pain syndrome type I with neridronate: a randomized, double-blind, placebo-controlled study.** *Rheumatology.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Varenna+neridronate+complex+regional+pain+syndrome+type+I+randomized+double+blind+placebo+Rheumatology+2013>
44. Varenna M, et al. **Predictors of responsiveness to bisphosphonate treatment in patients with complex regional pain syndrome type I / long-term follow-up.** *Rheumatology / Clin Exp Rheumatol.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Varenna+predictors+responsiveness+bisphosphonate+complex+regional+pain+syndrome+type+I>
45. Robinson JN, et al. **Efficacy of pamidronate in complex regional pain syndrome type I.** *Pain Med.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Robinson+efficacy+pamidronate+complex+regional+pain+syndrome+type+I+Pain+Medicine+2004>
46. Manicourt DH, et al. **Role of alendronate in therapy for posttraumatic complex regional pain syndrome type I of the lower extremity.** *Arthritis Rheum.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Manicourt+alendronate+posttraumatic+complex+regional+pain+syndrome+type+I+lower+extremity+Arthritis+Rheum>
47. Chevreau M, et al. **Bisphosphonates for treatment of complex regional pain syndrome type 1: a systematic literature review and meta-analysis.** *Joint Bone Spine.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Chevreau+bisphosphonates+complex+regional+pain+syndrome+type+1+systematic+review+meta-analysis>

## 12. 대뇌 재구성 · 재활 · 심리 (Cortical reorganisation, rehabilitation, psychology)

48. Maihöfner C, et al. **Cortical reorganization during recovery from complex regional pain syndrome.** *Neurology.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Maihofner+cortical+reorganization+during+recovery+complex+regional+pain+syndrome+Neurology+2004>
49. Moseley GL. **Graded motor imagery is effective for long-standing complex regional pain syndrome: a randomised controlled trial.** *Pain.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Moseley+graded+motor+imagery+effective+long+standing+complex+regional+pain+syndrome+randomised+controlled+trial>
50. Moseley GL. **Graded motor imagery for pathologic pain: a randomized controlled trial.** *Neurology.* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Moseley+graded+motor+imagery+pathologic+pain+randomized+controlled+trial+Neurology+2006>
51. McCabe CS, et al. **A controlled pilot study of the utility of mirror visual feedback in the treatment of complex regional pain syndrome (type 1).** *Rheumatology.* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/?term=McCabe+controlled+pilot+study+mirror+visual+feedback+complex+regional+pain+syndrome+type+1+Rheumatology>
52. Lewis JS, et al. **Body perception disturbance: a contribution to pain in complex regional pain syndrome.** *Pain.* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lewis+body+perception+disturbance+contribution+pain+complex+regional+pain+syndrome+Pain+2007>
53. de Jong JR, et al. **Reduction of pain-related fear in complex regional pain syndrome type I: the application of graded exposure in vivo.** *Pain.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/?term=de+Jong+reduction+pain+related+fear+complex+regional+pain+syndrome+graded+exposure+in+vivo+Pain+2005>
54. Barnhoorn KJ, et al. **Pain exposure physical therapy versus conventional treatment in complex regional pain syndrome type 1: a randomised controlled trial.** *BMJ Open.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Barnhoorn+pain+exposure+physical+therapy+conventional+treatment+complex+regional+pain+syndrome+randomised+BMJ+Open>
55. Bean DJ, et al. **Do psychological factors influence recovery from complex regional pain syndrome type 1? A prospective study.** *Pain.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bean+psychological+factors+influence+recovery+complex+regional+pain+syndrome+prospective+study+Pain>

## 13. 신경조절 · 중재 (Neuromodulation & interventions)

56. Kemler MA, et al. **Spinal cord stimulation in patients with chronic reflex sympathetic dystrophy.** *N Engl J Med.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kemler+spinal+cord+stimulation+chronic+reflex+sympathetic+dystrophy+New+England+Journal+2000>
57. Kemler MA, et al. **Effect of spinal cord stimulation for chronic complex regional pain syndrome type I: five-year final follow-up of patients in a randomized controlled trial.** *J Neurosurg.* 2008. (5년 차 군간 차이 소실)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kemler+spinal+cord+stimulation+five+year+final+follow+up+complex+regional+pain+syndrome+randomized>
58. Deer TR, et al. **Dorsal root ganglion stimulation yielded higher treatment success rate for CRPS and causalgia (ACCURATE study).** *Pain.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Deer+dorsal+root+ganglion+stimulation+higher+treatment+success+complex+regional+pain+syndrome+causalgia+ACCURATE>
59. van Rijn MA, et al. **Intrathecal baclofen for dystonia of complex regional pain syndrome.** *Pain.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=van+Rijn+intrathecal+baclofen+dystonia+complex+regional+pain+syndrome+Pain+2009>
60. Duong S, et al. **Treatment of complex regional pain syndrome: an updated systematic review and narrative synthesis.** *Can J Anaesth.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Duong+treatment+complex+regional+pain+syndrome+updated+systematic+review+narrative+synthesis>

## 14. 항염 · 대증 약물치료 (Anti-inflammatory & symptomatic pharmacotherapy)

61. Christensen K, et al. **The reflex dystrophy syndrome: response to treatment with systemic corticosteroids.** *Acta Chir Scand.* 1982.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Christensen+reflex+dystrophy+syndrome+response+treatment+systemic+corticosteroids+Acta+Chirurgica>
62. Kalita J, et al. **Comparison of prednisolone with piroxicam in complex regional pain syndrome following stroke: a randomized controlled trial.** *QJM.* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kalita+prednisolone+piroxicam+complex+regional+pain+syndrome+following+stroke+randomized+QJM>
63. van de Vusse AC, et al. **Randomised controlled trial of gabapentin in complex regional pain syndrome type 1.** *BMC Neurol.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/?term=van+de+Vusse+randomised+controlled+trial+gabapentin+complex+regional+pain+syndrome+type+1+BMC+Neurology>
64. Wertli MM, et al. **Rational pain management in complex regional pain syndrome 1 (CRPS 1) — a network meta-analysis.** *Pain Med.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Wertli+rational+pain+management+complex+regional+pain+syndrome+network+meta-analysis+Pain+Medicine>
65. O'Connell NE, et al. **Interventions for treating pain and disability in adults with complex regional pain syndrome (Cochrane overview).** 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=O%27Connell+interventions+treating+pain+disability+adults+complex+regional+pain+syndrome+Cochrane+overview>

## 15. QSP · 모델링 방법론 (QSP methodology)

66. Baral R, et al. / mrgsolve 문서. **mrgsolve: Simulate from ODE-based population PK/PD and QSP models in R.**
    <https://mrgsolve.org/>
67. Gadkar K, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol.* 2016.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Gadkar+six+stage+workflow+robust+application+systems+pharmacology>
68. Ribba B, et al. **Model-informed drug development / QSP for pain indications (review).**
    <https://pubmed.ncbi.nlm.nih.gov/?term=quantitative+systems+pharmacology+model+chronic+pain+central+sensitization>
69. Angeli D, Ferrell JE, Sontag ED. **Detection of multistability, bifurcations, and hysteresis in a large class of biological positive-feedback systems.** *PNAS.* 2004. (다중안정성 형식론)
    <https://pubmed.ncbi.nlm.nih.gov/?term=Angeli+Ferrell+Sontag+detection+multistability+bifurcations+hysteresis+positive+feedback+PNAS+2004>
70. Ferrell JE. **Bistability, bifurcations, and Waddington's epigenetic landscape.** *Curr Biol.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ferrell+bistability+bifurcations+Waddington+epigenetic+landscape+Current+Biology+2012>

---

### 모델 파라미터 정박 요약 (calibration anchors → parameters)

| 문헌 근거 | 모델 파라미터 / 결과 |
|---|---|
| Sigtermans 2009 (케타민 100시간, 주입 중 NRS 7.2→2.7, ~11주 후 소실) | `A_KET`, `EC50_KET`, `EMAX_KET` → 모델 6.04→3.02 (d34), d44에 무치료와 0.3 이내로 복귀 |
| Kemler 2000 / 2008 (SCS 6개월 양성, 5년 차 차이 소실) | `A_SCS`, `SCS_EFF`, `SCS_HAB_K` → 7개월 +3.27, 60개월 +0.15 |
| Varenna 2013 (네리드로네이트 100 mg IV ×4, 1년까지 효과 지속) | `EMAX_NER`, `BONE50_NER`, `W_BONE_PS` → 1년 NRS −1.83, CTX-I −61%, BMD 보존 |
| Goebel 2017 LIPS (장기 CRPS에서 IVIG 음성) | `W_AAB_PS` 작게 유지 → 후기 IVIG 효과 ≈ 0 |
| Zollinger 2007 (비타민 C 예방), Perez 2003 (DMSO/NAC) | `EMAX_NAC`, `EC50_NAC`, `KOUT_ROS` |
| Ji 2013/2019, Coull 2005 (교세포 priming·BDNF-KCC2) | `W_GLIA_SELF`, `GLIA50`, `HILL_GLIA`, `KOUT_GLIA` (잠금 구조) |
| Moseley 2004/2006, Maihöfner 2004 (사용 의존적 피질 재지도화) | `REHAB_CTX`, `KOUT_CTX`, `W_DIS_CTX` |
| de Jong 2005, Bean 2015 (통증 관련 공포가 예후를 결정) | `KFEAR`, `PAIN50_FEAR`, `HILL_FEAR`, `REHAB_FEAR` (분기 파라미터) |
| Eberle 2009, Wasner 2001, Groeneweg 2006 (온난→한랭 전환) | `W_NP_PERF`, `W_A1_PERF`, `W_SS_PERF`, `W_ROS_PERF` |
| Harden 2010/2017 (Budapest 기준, CSS 0-16) | `$TABLE`의 `CSS` 4-도메인 합산 및 `ACTIVE_CRPS` 플래그 |
