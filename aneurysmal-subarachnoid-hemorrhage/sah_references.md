# 참고문헌 — 동맥류성 지주막하출혈 후 지연성 뇌허혈 (aSAH-DCI)
# References — Aneurysmal Subarachnoid Haemorrhage & Delayed Cerebral Ischaemia

이 목록은 `sah_qsp_model.dot` · `sah_mrgsolve_model.R` · `sah_shiny_app.R` 이
주장하는 구조의 근거 문헌입니다. 각 항목의 주석은 "이 논문이 모델의 어느 항을
정당화하는가"를 적었습니다.

**모든 PMID는 PubMed E-utilities(esearch/esummary)로 조회해 제목·저널·연도를
확인한 것이며, 제목은 PubMed가 반환한 원문을 그대로 옮겼습니다.** 조회에
실패한 후보 문헌은 목록에서 제외했습니다(추측한 PMID를 넣지 않았습니다).

## 이 모델의 세 가지 구조적 주장과 그 근거

| 주장 | 내용 | 1차 근거 |
|------|------|----------|
| CLAIM 1 | 하나의 세동맥 예비능을 네 소비자가 나눠 쓴다 — 대혈관 경련 · 미세혈관 긴장 · 미세혈전 · 관류압/SD | Rowland 2012 · Macdonald 2014 · Etminan 2011 |
| CLAIM 2 | 미세순환은 저항의 각주가 아니라 산소 추출 한계를 깎는다 (CTH → OEFmax↓) | Jespersen & Østergaard 2012 · Østergaard 2013 |
| CLAIM 3 | 시간창은 달력이 아니라 헤모글로빈 처리 동역학이다 (전이 사슬 + 유도성 싱크) | Suzuki 1999 (HO-1) · Hugelshofer 2019 (haptoglobin) |

가장 중요한 단일 대조: **CONSCIOUS-1**(혈관경련 65% 감소)과
**CONSCIOUS-2**(결과 변화 없음)를 같은 약물·같은 용량에서 나란히 읽는 것.
이 모델은 그 간격을 '중복성(redundancy) 먼저, 위해 채널(harm channel) 나중'으로
분해해서 계산합니다.

---

**총 89편** (모두 PMID 검증 완료).

## 1. 개요 · 역학 · 정의 및 척도  Overview, Epidemiology, Definitions & Scales

1. **Subarachnoid haemorrhage.** *Lancet* 2007. [PMID 17258671](https://pubmed.ncbi.nlm.nih.gov/17258671/)  
   고전적 총론. 본 모델의 임상 골격(등급, 시간축, 합병증)의 출발점.

2. **Spontaneous subarachnoid haemorrhage.** *Lancet* 2017. [PMID 27637674](https://pubmed.ncbi.nlm.nih.gov/27637674/)  
   현대적 총론. DCI를 혈관경련과 분리해 서술하는 관점의 대표 문헌.

3. **2023 Guideline for the Management of Patients With Aneurysmal Subarachnoid Hemorrhage: A Guideline From the American Heart Association/American Stroke Association.** *Stroke* 2023. [PMID 37212182](https://pubmed.ncbi.nlm.nih.gov/37212182/)  
   2023 AHA/ASA 지침. 니모디핀 권고, 유도 고혈압의 불확실성, 목표 혈압·수액 관리의 근거.

4. **Comment on the 2023 Guidelines for the Management of Patients With Aneurysmal Subarachnoid Hemorrhage.** *Stroke* 2023. [PMID 37581267](https://pubmed.ncbi.nlm.nih.gov/37581267/)  
   위 지침에 대한 논평. 권고 강도의 해석 논쟁을 보여준다.

5. **Critical care management of patients following aneurysmal subarachnoid hemorrhage: recommendations from the Neurocritical Care Society's Multidisciplinary Consensus Conference.** *Neurocrit Care* 2011. [PMID 21773873](https://pubmed.ncbi.nlm.nih.gov/21773873/)  
   Neurocritical Care Society 합의. 모니터링(TCD·PbtO2·미세투석) 사용 근거.

6. **Definition of delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage as an outcome event in clinical trials and observational studies: proposal of a multidisciplinary research group.** *Stroke* 2010. [PMID 20798370](https://pubmed.ncbi.nlm.nih.gov/20798370/)  
   DCI의 표준 정의(신경학적 악화 + 새 경색). 모델의 DCI 판정 규칙이 이 정의를 조작화한 것.

7. **Report of World Federation of Neurological Surgeons Committee on a Universal Subarachnoid Hemorrhage Grading Scale.** *J Neurosurg* 1988. [PMID 3131498](https://pubmed.ncbi.nlm.nih.gov/3131498/)  
   WFNS 등급 원전. 모델의 EBIB(초기 뇌손상 부담) 초기값이 WFNS에서 유도된다.

8. **Prediction of symptomatic vasospasm after subarachnoid hemorrhage: the modified fisher scale.** *Neurosurgery* 2006. [PMID 16823296](https://pubmed.ncbi.nlm.nih.gov/16823296/)  
   modified Fisher 척도. 모델의 CLOT0·IVH0 초기 조건이 이 등급에서 결정된다.

9. **Reliability of the modified Rankin Scale applied by telephone.** *Neurol Int* 2013. [PMID 23717781](https://pubmed.ncbi.nlm.nih.gov/23717781/)  
   mRS 신뢰도. 결과 로짓의 종말점 정의에 대한 측정오차 배경.


## 2. 혈관경련–DCI 해리: 이 모델이 설명하려는 현상  The Vasospasm–DCI Dissociation

10. **Delayed neurological deterioration after subarachnoid haemorrhage.** *Nat Rev Neurol* 2014. [PMID 24323051](https://pubmed.ncbi.nlm.nih.gov/24323051/)  
   지연성 신경학적 악화를 혈관경련 하나로 설명할 수 없다는 핵심 리뷰. 본 모델의 출발 문제의식.

11. **Delayed cerebral ischaemia after subarachnoid haemorrhage: looking beyond vasospasm.** *Br J Anaesth* 2012. [PMID 22879655](https://pubmed.ncbi.nlm.nih.gov/22879655/)  
   제목 그대로 'DCI: 혈관경련을 넘어서'. 4개 소비자 구조의 서술적 선행 문헌.

12. **Relationship between vasospasm, cerebral perfusion, and delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage.** *Neuroradiology* 2009. [PMID 19623472](https://pubmed.ncbi.nlm.nih.gov/19623472/)  
   혈관경련–관류–DCI의 삼각관계가 느슨하다는 것을 정량적으로 보인 연구.

13. **Relationship between angiographic vasospasm and regional hypoperfusion in aneurysmal subarachnoid hemorrhage.** *Stroke* 2012. [PMID 22492520](https://pubmed.ncbi.nlm.nih.gov/22492520/)  
   혈관경련과 국소 저관류의 관계가 부분적임을 보인 연구. 대혈관 경로의 설명력 한계.

14. **The relationship between delayed infarcts and angiographic vasospasm after aneurysmal subarachnoid hemorrhage.** *Neurosurgery* 2013. [PMID 23313984](https://pubmed.ncbi.nlm.nih.gov/23313984/)  
   지연성 경색과 혈관조영 경련의 불일치를 직접 측정. 모델의 r² 계산과 대응.

15. **Angiographic vasospasm is strongly correlated with cerebral infarction after subarachnoid hemorrhage.** *Stroke* 2011. [PMID 21350201](https://pubmed.ncbi.nlm.nih.gov/21350201/)  
   반대 방향의 근거: 경련과 경색의 강한 상관을 보고. 모델은 양쪽을 모두 설명해야 한다.

16. **Predictors of cerebral infarction in aneurysmal subarachnoid hemorrhage.** *Stroke* 2004. [PMID 15218156](https://pubmed.ncbi.nlm.nih.gov/15218156/)  
   뇌경색의 예측인자. 임상 등급·출혈량이 경련보다 강한 예측인자라는 점.

17. **Effect of pharmaceutical treatment on vasospasm, delayed cerebral ischemia, and clinical outcome in patients with aneurysmal subarachnoid hemorrhage: a systematic review and meta-analysis.** *J Cereb Blood Flow Metab* 2011. [PMID 21285966](https://pubmed.ncbi.nlm.nih.gov/21285966/)  
   결정적 메타분석: 약물치료가 혈관경련은 줄이지만 결과는 바꾸지 못한다는 것을 시험 전체에서 보임. 본 모델이 설명하려는 현상 그 자체.

18. **Effect of endothelin-receptor antagonists on delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage remains unclear.** *Stroke* 2009. [PMID 19875735](https://pubmed.ncbi.nlm.nih.gov/19875735/)  
   엔도텔린 수용체 길항제 메타분석. 경련 감소 vs 결과 무변화의 정량적 요약.


## 3. 클라조센탄과 엔도텔린 축: 교훈적 실패  Clazosentan & the Endothelin Axis

19. **Clazosentan to overcome neurological ischemia and infarction occurring after subarachnoid hemorrhage (CONSCIOUS-1): randomized, double-blind, placebo-controlled phase 2 dose-finding trial.** *Stroke* 2008. [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/)  
   CONSCIOUS-1. 클라조센탄이 중등–중증 혈관경련을 약 65% 감소시킴(모델 보정 목표 RR 0.35).

20. **Clazosentan, an endothelin receptor antagonist, in patients with aneurysmal subarachnoid haemorrhage undergoing surgical clipping: a randomised, double-blind, placebo-controlled phase 3 trial (CONSCIOUS-2).** *Lancet Neurol* 2011. [PMID 21640651](https://pubmed.ncbi.nlm.nih.gov/21640651/)  
   CONSCIOUS-2. 같은 약이 3개월 결과를 개선하지 못함(모델 보정 목표 RR ~1.05). 본 모델의 중심 사례.

21. **Randomized trial of clazosentan in patients with aneurysmal subarachnoid hemorrhage undergoing endovascular coiling.** *Stroke* 2012. [PMID 22403047](https://pubmed.ncbi.nlm.nih.gov/22403047/)  
   CONSCIOUS-3. 코일링 환자에서의 확인. 고용량에서 폐부종·저혈압 등 위해 채널이 뚜렷.

22. **Effects of clazosentan on cerebral vasospasm-related morbidity and all-cause mortality after aneurysmal subarachnoid hemorrhage: two randomized phase 3 trials in Japanese patients.** *J Neurosurg* 2022. [PMID 35364589](https://pubmed.ncbi.nlm.nih.gov/35364589/)  
   일본 제3상. 혈관경련 관련 이환·사망을 감소시켜, 종말점 정의에 따라 결론이 달라짐을 보여준다.

23. **The REACT study: design of a randomized phase 3 trial to assess the efficacy and safety of clazosentan for preventing deterioration due to delayed cerebral ischemia after aneurysmal subarachnoid hemorrhage.** *BMC Neurol* 2022. [PMID 36539711](https://pubmed.ncbi.nlm.nih.gov/36539711/)  
   REACT 시험 설계. 클라조센탄 재평가의 현재 위치.

24. **Endothelin and subarachnoid hemorrhage: an overview.** *Neurosurgery* 1998. [PMID 9766314](https://pubmed.ncbi.nlm.nih.gov/9766314/)  
   엔도텔린 축 총론. ET-1/ETA/ETB 구조와 모델의 ETeff 항.

25. **Plasma endothelin-1 as screening marker for cerebral vasospasm after subarachnoid hemorrhage.** *Neurocrit Care* 2014. [PMID 23921571](https://pubmed.ncbi.nlm.nih.gov/23921571/)  
   혈장 ET-1의 선별검사 성능. 효과기 농도와 임상 사건의 느슨한 연결.


## 4. 니모디핀과 칼슘 길항제  Nimodipine and Calcium Antagonists

26. **Effect of oral nimodipine on cerebral infarction and outcome after subarachnoid haemorrhage: British aneurysm nimodipine trial.** *BMJ* 1989. [PMID 2496789](https://pubmed.ncbi.nlm.nih.gov/2496789/)  
   BRANT. 경구 니모디핀이 뇌경색과 불량 결과를 감소. 유일하게 결과를 바꾼 약.

27. **Cerebral arterial spasm--a controlled trial of nimodipine in patients with subarachnoid hemorrhage.** *N Engl J Med* 1983. [PMID 6338383](https://pubmed.ncbi.nlm.nih.gov/6338383/)  
   최초의 니모디핀 RCT. 혈관조영 경련의 변화는 미미했으나 임상 결과가 개선된 원형 사례.

28. **Calcium antagonists for aneurysmal subarachnoid haemorrhage.** *Cochrane Database Syst Rev* 2007. [PMID 17636626](https://pubmed.ncbi.nlm.nih.gov/17636626/)  
   칼슘 길항제 Cochrane 리뷰. 불량 결과 RR 약 0.67(모델 보정 목표).

29. **Nimodipine-Induced Blood Pressure Changes Can Predict Delayed Cerebral Ischemia.** *Front Neurol* 2019. [PMID 31736865](https://pubmed.ncbi.nlm.nih.gov/31736865/)  
   니모디핀에 의한 혈압 변화 자체가 DCI를 예측. 모델의 CPP 위해 채널(DMAP_NIM)의 직접 근거.

30. **NEWTON: Nimodipine Microparticles to Enhance Recovery While Reducing Toxicity After Subarachnoid Hemorrhage.** *Neurocrit Care* 2015. [PMID 25678453](https://pubmed.ncbi.nlm.nih.gov/25678453/)  
   뇌실내 서방형 니모디핀(EG-1962) 제1/2상. 전신 저혈압 없이 국소 노출을 높이는 전략.

31. **NEWTON-2 Cisternal (Nimodipine Microparticles to Enhance Recovery While Reducing Toxicity After Subarachnoid Hemorrhage): A Phase 2, Multicenter, Randomized, Open-Label Safety Study of Intracisternal EG-1962 in Aneurysmal Subarachnoid Hemorrhage.** *Neurosurgery* 2020. [PMID 32985652](https://pubmed.ncbi.nlm.nih.gov/32985652/)  
   NEWTON-2. 국소 전달이 혈관경련을 줄여도 결과는 개선하지 못한 또 하나의 사례.

32. **Nicardipine prolonged-release implants for preventing cerebral vasospasm after subarachnoid hemorrhage: effect and outcome in the first 100 patients.** *Neurol Med Chir (Tokyo)* 2007. [PMID 17895611](https://pubmed.ncbi.nlm.nih.gov/17895611/)  
   니카르디핀 서방형 삽입물. 국소 대혈관 이완의 강한 효과.

33. **Intrathecal nicardipine for cerebral vasospasm after non-traumatic subarachnoid hemorrhage: a meta-analysis.** *Neurosurg Rev* 2025. [PMID 40295406](https://pubmed.ncbi.nlm.nih.gov/40295406/)  
   척수강내 니카르디핀 메타분석. 모델의 S9(척수강내 니카르디핀) 시나리오 근거.


## 5. 기타 개입 임상시험  Other Interventional Trials

34. **Simvastatin in aneurysmal subarachnoid haemorrhage (STASH): a multicentre randomised phase 3 trial.** *Lancet Neurol* 2014. [PMID 24837690](https://pubmed.ncbi.nlm.nih.gov/24837690/)  
   STASH. 심바스타틴이 결과를 개선하지 못함(모델 보정 목표 RR ~1.0).

35. **Magnesium for aneurysmal subarachnoid haemorrhage (MASH-2): a randomised placebo-controlled trial.** *Lancet* 2012. [PMID 22633825](https://pubmed.ncbi.nlm.nih.gov/22633825/)  
   MASH-2. 마그네슘 무효.

36. **Intravenous magnesium sulphate for aneurysmal subarachnoid hemorrhage (IMASH): a randomized, double-blinded, placebo-controlled, multicenter phase III trial.** *Stroke* 2010. [PMID 20378868](https://pubmed.ncbi.nlm.nih.gov/20378868/)  
   IMASH. 마그네슘 무효(확인).

37. **Effects of cilostazol on cerebral vasospasm after aneurysmal subarachnoid hemorrhage: a multicenter prospective, randomized, open-label blinded end point trial.** *J Neurosurg* 2013. [PMID 23039152](https://pubmed.ncbi.nlm.nih.gov/23039152/)  
   실로스타졸 RCT. 항혈소판+미세혈관 확장이라는 다경로 작용.

38. **Efficacy of Cilostazol in Prevention of Delayed Cerebral Ischemia after Aneurysmal Subarachnoid Hemorrhage: A Meta-Analysis.** *J Stroke Cerebrovasc Dis* 2018. [PMID 30093204](https://pubmed.ncbi.nlm.nih.gov/30093204/)  
   실로스타졸 메타분석. DCI RR 약 0.47(모델 보정 목표).

39. **Effectiveness of Lumbar Cerebrospinal Fluid Drain Among Patients With Aneurysmal Subarachnoid Hemorrhage: A Randomized Clinical Trial.** *JAMA Neurol* 2023. [PMID 37330974](https://pubmed.ncbi.nlm.nih.gov/37330974/)  
   EARLYDRAIN. 조기 요추 배액이 불량 결과를 감소. 모델에서 유일하게 '입력 함수를 절단'하는 개입.

40. **Induced Hypertension for Delayed Cerebral Ischemia After Aneurysmal Subarachnoid Hemorrhage: A Randomized Clinical Trial.** *Stroke* 2018. [PMID 29158449](https://pubmed.ncbi.nlm.nih.gov/29158449/)  
   HIMALAIA. 유도 고혈압의 무작위 시험이 무효. 모델은 자동조절능 상태에 따라 효과가 갈리는 것으로 설명.

41. **Ultra-early tranexamic acid after subarachnoid haemorrhage (ULTRA): a randomised controlled trial.** *Lancet* 2021. [PMID 33357465](https://pubmed.ncbi.nlm.nih.gov/33357465/)  
   ULTRA. 초조기 트라넥사민산이 재출혈은 줄이나 결과는 개선하지 못함.

42. **Milrinone and homeostasis to treat cerebral vasospasm associated with subarachnoid hemorrhage: the Montreal Neurological Hospital protocol.** *Neurocrit Care* 2012. [PMID 22528278](https://pubmed.ncbi.nlm.nih.gov/22528278/)  
   Montreal 밀리논 프로토콜. 미세혈관 확장 + 심박출 증가의 구제 요법.

43. **Effect of prophylactic transluminal balloon angioplasty on cerebral vasospasm and outcome in patients with Fisher grade III subarachnoid hemorrhage: results of a phase II multicenter, randomized, clinical trial.** *Stroke* 2008. [PMID 18420953](https://pubmed.ncbi.nlm.nih.gov/18420953/)  
   예방적 풍선 혈관성형술. 대혈관을 물리적으로 넓혀도 결과는 개선되지 않은 사례.

44. **Endovascular Rescue Treatment for Delayed Cerebral Ischemia After Subarachnoid Hemorrhage Is Safe and Effective.** *Front Neurol* 2019. [PMID 30858818](https://pubmed.ncbi.nlm.nih.gov/30858818/)  
   혈관내 구제 치료. 실제 임상에서 대혈관 경로에 개입하는 마지막 수단.


## 6. 헤모글로빈 · 헵토글로빈 · HO-1: 시계  The Haemoglobin Clock

45. **A review of hemoglobin and the pathogenesis of cerebral vasospasm.** *Stroke* 1991. [PMID 1866764](https://pubmed.ncbi.nlm.nih.gov/1866764/)  
   헤모글로빈과 혈관경련 병인 리뷰. 모델의 OXYHB 중심 구조의 근거.

46. **Haptoglobin administration into the subarachnoid space prevents hemoglobin-induced cerebral vasospasm.** *J Clin Invest* 2019. [PMID 31454333](https://pubmed.ncbi.nlm.nih.gov/31454333/)  
   헵토글로빈을 지주막하 공간에 투여하면 혈관경련이 예방됨. 소거 용량이 율속 단계임을 증명.

47. **Haptoglobin genotype and outcome after aneurysmal subarachnoid haemorrhage.** *J Neurol Neurosurg Psychiatry* 2020. [PMID 31937585](https://pubmed.ncbi.nlm.nih.gov/31937585/)  
   헵토글로빈 유전형과 결과. 모델의 HP22 공변량.

48. **Haptoglobin phenotype predicts cerebral vasospasm and clinical deterioration after aneurysmal subarachnoid hemorrhage.** *J Stroke Cerebrovasc Dis* 2013. [PMID 23498376](https://pubmed.ncbi.nlm.nih.gov/23498376/)  
   헵토글로빈 표현형이 혈관경련과 임상 악화를 예측.

49. **Heme oxygenase-1 gene induction as an intrinsic regulation against delayed cerebral vasospasm in rats.** *J Clin Invest* 1999. [PMID 10393699](https://pubmed.ncbi.nlm.nih.gov/10393699/)  
   HO-1 유도가 지연성 혈관경련에 대한 내인성 방어. 모델의 '유도성 싱크'와 지연(t½ 1.9 d).

50. **Bilirubin oxidation products (BOXes) and their role in cerebral vasospasm after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2006. [PMID 16467784](https://pubmed.ncbi.nlm.nih.gov/16467784/)  
   빌리루빈 산화생성물(BOXes)의 혈관수축 작용. 모델의 BOX 상태.

51. **Iron and early brain injury after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2010. [PMID 20736954](https://pubmed.ncbi.nlm.nih.gov/20736954/)  
   철과 초기 뇌손상. Fenton 화학과 ROS 항의 근거.

52. **Ferroptosis mechanisms in early brain injury after subarachnoid hemorrhage.** *Brain Res* 2026. [PMID 41802694](https://pubmed.ncbi.nlm.nih.gov/41802694/)  
   헴/철 매개 페롭토시스. 산화 손상 경로의 최신 정리.


## 7. NO · 엔도텔린 · Rho-kinase: 효과기  Effector Signalling

53. **Acute decrease in cerebral nitric oxide levels after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2000. [PMID 10724124](https://pubmed.ncbi.nlm.nih.gov/10724124/)  
   SAH 직후 뇌 NO 농도의 급성 감소. 모델의 NOB 항.

54. **Delayed cerebral vasospasm and nitric oxide: review, new hypothesis, and proposed treatment.** *Pharmacol Ther* 2005. [PMID 15626454](https://pubmed.ncbi.nlm.nih.gov/15626454/)  
   NO 결핍 가설과 치료 제안. oxyHb의 NO 소거가 중심.

55. **Nitric oxide in early brain injury after subarachnoid hemorrhage.** *Acta Neurochir Suppl* 2011. [PMID 21116923](https://pubmed.ncbi.nlm.nih.gov/21116923/)  
   초기 뇌손상에서의 NO. EBI와 혈관 기능부전의 연결.

56. **Involvement of Rho-kinase-mediated phosphorylation of myosin light chain in enhancement of cerebral vasospasm.** *Circ Res* 2000. [PMID 10926869](https://pubmed.ncbi.nlm.nih.gov/10926869/)  
   Rho-kinase 매개 MLC 인산화가 혈관경련을 증강. 모델의 RHOK(Ca 민감화) 항 — 칼슘 차단만으로 부족한 이유.

57. **Effect of AT877 on cerebral vasospasm after aneurysmal subarachnoid hemorrhage. Results of a prospective placebo-controlled double-blind trial.** *J Neurosurg* 1992. [PMID 1545249](https://pubmed.ncbi.nlm.nih.gov/1545249/)  
   파수딜(AT877) 무작위 시험. Rho-kinase 억제의 임상 근거.


## 8. 미세순환 · 통과시간 이질성 · 미세혈전  Microcirculation, CTH & Microthrombosis

58. **The role of the microcirculation in delayed cerebral ischemia and chronic degenerative changes after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2013. [PMID 24064495](https://pubmed.ncbi.nlm.nih.gov/24064495/)  
   미세순환이 DCI에서 하는 역할. 본 모델 CLAIM 2의 1차 출처.

59. **The roles of cerebral blood flow, capillary transit time heterogeneity, and oxygen tension in brain oxygenation and metabolism.** *J Cereb Blood Flow Metab* 2012. [PMID 22044867](https://pubmed.ncbi.nlm.nih.gov/22044867/)  
   CBF·통과시간 이질성(CTH)·산소분압의 관계. 모델의 OEFmax = 0.85/(1+K·CTH) 형태가 여기서 나온다.

60. **Are We Barking Up the Wrong Vessels? Cerebral Microcirculation After Subarachnoid Hemorrhage.** *Stroke* 2015. [PMID 26152299](https://pubmed.ncbi.nlm.nih.gov/26152299/)  
   '우리가 잘못된 혈관을 보고 있는가?' 미세순환 관점의 선언적 리뷰.

61. **Capillary flow disturbances after experimental subarachnoid hemorrhage: A contributor to delayed cerebral ischemia?.** *Microcirculation* 2019. [PMID 30431201](https://pubmed.ncbi.nlm.nih.gov/30431201/)  
   실험적 SAH 후 모세혈관 흐름 장애가 DCI에 기여. CTH의 실험적 근거.

62. **Microthrombosis after aneurysmal subarachnoid hemorrhage: an additional explanation for delayed cerebral ischemia.** *J Cereb Blood Flow Metab* 2008. [PMID 18628782](https://pubmed.ncbi.nlm.nih.gov/18628782/)  
   미세혈전이 DCI의 추가 설명이라는 대표 문헌. 모델의 MTHR 상태.

63. **Mechanisms of microthrombosis and microcirculatory constriction after experimental subarachnoid hemorrhage.** *Acta Neurochir Suppl* 2013. [PMID 22890667](https://pubmed.ncbi.nlm.nih.gov/22890667/)  
   미세혈전 형성과 미세순환 수축의 기전.

64. **Intraoperative detection of early microvasospasm in patients with subarachnoid hemorrhage by using orthogonal polarization spectral imaging.** *Neurosurgery* 2003. [PMID 12762876](https://pubmed.ncbi.nlm.nih.gov/12762876/)  
   수술 중 조기 미세혈관 경련의 직접 관찰. 대혈관 경련과 독립적인 경로.

65. **Reduction of neutrophil activity decreases early microvascular injury after subarachnoid haemorrhage.** *J Neuroinflammation* 2011. [PMID 21854561](https://pubmed.ncbi.nlm.nih.gov/21854561/)  
   호중구 활성 감소가 초기 미세혈관 손상을 줄임. 염증–미세순환 연결.

66. **Pericyte constriction after stroke: the jury is still out.** *Nat Med* 2010. [PMID 20823870](https://pubmed.ncbi.nlm.nih.gov/20823870/)  
   뇌졸중 후 페리사이트 수축 논쟁.

67. **Capillary pericytes regulate cerebral blood flow in health and disease.** *Nature* 2014. [PMID 24670647](https://pubmed.ncbi.nlm.nih.gov/24670647/)  
   모세혈관 페리사이트가 CBF를 조절한다는 Nature 연구. 미세혈관 저항의 세포학적 기반.


## 9. 확산성 탈분극  Spreading Depolarisation

68. **Delayed ischaemic neurological deficits after subarachnoid haemorrhage are associated with clusters of spreading depolarizations.** *Brain* 2006. [PMID 17067993](https://pubmed.ncbi.nlm.nih.gov/17067993/)  
   SAH 후 지연성 허혈 결손이 확산성 탈분극 군집과 연관. 모델의 SD 경로의 1차 출처.

69. **The role of spreading depression, spreading depolarization and spreading ischemia in neurological disease.** *Nat Med* 2011. [PMID 21475241](https://pubmed.ncbi.nlm.nih.gov/21475241/)  
   확산성 탈분극/확산성 허혈의 총론. 역 신경혈관 결합(inverse coupling)의 개념적 근거.

70. **Recording, analysis, and interpretation of spreading depolarizations in neurointensive care: Review and recommendations of the COSBID research group.** *J Cereb Blood Flow Metab* 2017. [PMID 27317657](https://pubmed.ncbi.nlm.nih.gov/27317657/)  
   COSBID 권고. SD 기록·해석 표준.

71. **Delayed cerebral ischemia and spreading depolarization in absence of angiographic vasospasm after subarachnoid hemorrhage.** *J Cereb Blood Flow Metab* 2012. [PMID 22146193](https://pubmed.ncbi.nlm.nih.gov/22146193/)  
   혈관조영 경련이 없는 상태에서의 DCI와 SD. 모델이 SD를 독립 소비자로 두는 직접 근거.

72. **Spreading Depolarizations: A Therapeutic Target Against Delayed Cerebral Ischemia After Subarachnoid Hemorrhage.** *J Clin Neurophysiol* 2016. [PMID 27258442](https://pubmed.ncbi.nlm.nih.gov/27258442/)  
   SD를 치료 표적으로 보는 관점. 케타민 시나리오의 근거.

73. **Preliminary evidence that ketamine inhibits spreading depolarizations in acute human brain injury.** *Stroke* 2009. [PMID 19520992](https://pubmed.ncbi.nlm.nih.gov/19520992/)  
   케타민이 인간 뇌에서 SD를 억제한다는 예비 근거. 모델의 S12 시나리오.


## 10. 뇌 자동조절능  Cerebral Autoregulation

74. **Clinical relevance of cerebral autoregulation following subarachnoid haemorrhage.** *Nat Rev Neurol* 2013. [PMID 23419369](https://pubmed.ncbi.nlm.nih.gov/23419369/)  
   SAH 후 자동조절능의 임상적 의미. 모델의 AREG 상태와 압력 수동성.

75. **Impairment of cerebral autoregulation predicts delayed cerebral ischemia after subarachnoid hemorrhage: a prospective observational study.** *Stroke* 2012. [PMID 23150652](https://pubmed.ncbi.nlm.nih.gov/23150652/)  
   자동조절능 손상이 DCI를 예측. AREG를 위험 인자로 두는 근거이자 R5의 층화 결과와 대응.


## 11. 모니터링과 바이오마커  Monitoring & Biomarkers

76. **Transcranial Doppler versus angiography in patients with vasospasm due to a ruptured cerebral aneurysm: A systematic review.** *Stroke* 2001. [PMID 11588316](https://pubmed.ncbi.nlm.nih.gov/11588316/)  
   TCD vs 혈관조영 체계적 문헌고찰. 모델의 TCD 진단 성능 계산과 비교 대상.

77. **How is vasospasm screening using transcranial Doppler associated with delayed cerebral ischemia and outcomes in aneurysmal subarachnoid hemorrhage?.** *Acta Neurochir (Wien)* 2019. [PMID 30637487](https://pubmed.ncbi.nlm.nih.gov/30637487/)  
   TCD 선별검사와 DCI·결과의 연관. 유속이 흐름에 교란된다는 점.

78. **Monitoring of brain tissue oxygenation following severe subarachnoid hemorrhage.** *Neurol Res* 2003. [PMID 12866190](https://pubmed.ncbi.nlm.nih.gov/12866190/)  
   중증 SAH에서의 뇌조직 산소 모니터링. PbtO2가 공유 노드를 직접 측정하는 이유.

79. **Cerebral metabolism and intracranial hypertension in high grade aneurysmal subarachnoid haemorrhage patients.** *Acta Neurochir Suppl* 2005. [PMID 16463827](https://pubmed.ncbi.nlm.nih.gov/16463827/)  
   미세투석 대사 지표와 두개내압. LPR 상승의 근거.


## 12. 초기 뇌손상과 전신 합병증  Early Brain Injury & Systemic Complications

80. **Metamorphosis of subarachnoid hemorrhage research: from delayed vasospasm to early brain injury.** *Mol Neurobiol* 2011. [PMID 21161614](https://pubmed.ncbi.nlm.nih.gov/21161614/)  
   연구 패러다임이 지연성 혈관경련에서 초기 뇌손상으로 이동. 모델에서 EBI가 결과의 최대 결정인자인 이유.

81. **Early Brain Injury After Subarachnoid Hemorrhage: Incidence and Mechanisms.** *Stroke* 2023. [PMID 36866673](https://pubmed.ncbi.nlm.nih.gov/36866673/)  
   초기 뇌손상의 발생률과 기전(2023 종합).

82. **The incidence and pathophysiology of hyponatraemia after subarachnoid haemorrhage.** *Clin Endocrinol (Oxf)* 2006. [PMID 16487432](https://pubmed.ncbi.nlm.nih.gov/16487432/)  
   SAH 후 저나트륨혈증의 빈도와 기전. 모델의 NAS 상태 보정(<135 약 30-50%).

83. **Higher hemoglobin is associated with improved outcome after subarachnoid hemorrhage.** *Crit Care Med* 2007. [PMID 17717494](https://pubmed.ncbi.nlm.nih.gov/17717494/)  
   높은 헤모글로빈이 더 좋은 결과와 연관. CaO2가 같은 노드에 들어간다는 임상 근거.

84. **Red blood cell transfusion in patients with subarachnoid hemorrhage: a multidisciplinary North American survey.** *Crit Care* 2011. [PMID 21244675](https://pubmed.ncbi.nlm.nih.gov/21244675/)  
   SAH에서의 적혈구 수혈 실태. S13 시나리오의 임상적 맥락.

85. **Impact of cardiac complications on outcome after aneurysmal subarachnoid hemorrhage: a meta-analysis.** *Neurology* 2009. [PMID 19221297](https://pubmed.ncbi.nlm.nih.gov/19221297/)  
   심장 합병증이 결과에 미치는 영향(메타분석). 신경성 심근 손상 채널.

86. **Fever after subarachnoid hemorrhage: risk factors and impact on outcome.** *Neurology* 2007. [PMID 17314332](https://pubmed.ncbi.nlm.nih.gov/17314332/)  
   발열의 결과 영향. 모델에서 발열은 CMRO2(수요) 측 소비자.

87. **Risk factors for surgery-related cerebral infarction after aneurysmal subarachnoid hemorrhage.** *Neurosurg Rev* 2025. [PMID 40316796](https://pubmed.ncbi.nlm.nih.gov/40316796/)  
   수술 관련 뇌경색의 위험인자. 경색 부담의 비-DCI 성분.


## 13. QSP 방법론  QSP Methodology

88. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)  
   mrgsolve를 이용한 QSP/PBPK 실습 튜토리얼. 본 모델의 구현 도구.

89. **Quantitative Systems Pharmacology: A Case for Disease Models.** *Clin Pharmacol Ther* 2017. [PMID 27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)  
   질환 모델로서의 QSP 논의. 본 라이브러리의 방법론적 위치.



---

## 보정 목표와 모델 산출값 (Calibration targets vs model)

아래 값은 모두 `sah_reference_check.py`(numpy/RK4 독립 구현)가 계산한 것이며,
전체 출력은 `sah_reference_check_output.txt`에 그대로 커밋되어 있습니다.

| 항목 | 문헌 목표 | 출처 |
|------|-----------|------|
| 중등–중증 혈관조영 혈관경련 (대조군) | 약 66% | CONSCIOUS-1 [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/) |
| DCI 발생률 (표준 치료) | 약 30% | Vergouwen 2010 [PMID 20798370](https://pubmed.ncbi.nlm.nih.gov/20798370/) |
| 불량 결과 mRS 4-6 (90일) | 30–35% | Macdonald 2017 [PMID 27637674](https://pubmed.ncbi.nlm.nih.gov/27637674/) |
| 저나트륨혈증 < 135 mmol/L | 30–50% | Sherlock 2006 [PMID 16487432](https://pubmed.ncbi.nlm.nih.gov/16487432/) |
| 클라조센탄 RR 중등–중증 경련 | 0.35 | CONSCIOUS-1 [PMID 18688013](https://pubmed.ncbi.nlm.nih.gov/18688013/) |
| 클라조센탄 RR 불량 결과 | 약 1.05 | CONSCIOUS-2 [PMID 21640651](https://pubmed.ncbi.nlm.nih.gov/21640651/) |
| 니모디핀 RR 불량 결과 | 0.67 | Cochrane [PMID 17636626](https://pubmed.ncbi.nlm.nih.gov/17636626/) |
| 실로스타졸 RR DCI | 약 0.47 | 메타분석 [PMID 30093204](https://pubmed.ncbi.nlm.nih.gov/30093204/) |
| 심바스타틴 RR 불량 결과 | 약 1.00 | STASH [PMID 24837690](https://pubmed.ncbi.nlm.nih.gov/24837690/) |
| 요추 배액 RR 불량 결과 | 0.76 | EARLYDRAIN [PMID 37330974](https://pubmed.ncbi.nlm.nih.gov/37330974/) |

## 모델이 문헌보다 얕거나 다른 지점 (기록해 둔 불일치)

튜닝으로 지우지 않고 명시적으로 남긴 항목입니다. 전체 목록과 설명은 이 디렉토리의
[README.md](README.md) "모델이 문헌과 어긋나거나 실패한 지점" 절에 있습니다.

1. **니모디핀의 결과 개선폭이 얕다.** 모델 RR(불량 결과) 0.80 대 Cochrane 0.67.
2. **혈관조영 경련이 없는 DCI의 비율이 낮다.** 모델 10.9% 대 문헌 보고 20-30%.
3. **modified Fisher 양극단이 벌어져 있다.** mF1 2.3%(문헌 약 6-12%),
   mF4 50.0%(문헌 약 35-40%). 순서는 단조롭게 맞지만 기울기가 급하다.
4. **저나트륨혈증 < 130 mmol/L 이 낮다.** 모델 4.6% 대 문헌 10-20%.
5. **가설이 재현되지 않은 항목.** "유도 고혈압은 압력 수동 환자에서만 듣는다"는
   예측이 계산에서 반대로 나왔다(자동조절 손상군 DCI RR 0.73 vs 온전군 0.62).
   HIMALAIA 무효를 자동조절 층화로 설명하려는 가설은 이 모델에서 지지되지 않는다.

## 방법론 도구

- mrgsolve (R): <https://mrgsolve.org>
- gPKPDviz — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구
  - 논문: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
  - 코드: <https://github.com/Genentech/gPKPDviz/>
- PubMed E-utilities (본 목록의 PMID 검증에 사용): <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
