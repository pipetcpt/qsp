# 자발성 뇌내출혈 (Spontaneous Intracerebral Haemorrhage, ICH) — 참고문헌

QSP 모델(`ich_qsp_model.dot`, `ich_mrgsolve_model.R`, `ich_shiny_app.R`)의 구조·파라미터·
검증 앵커에 사용된 문헌 목록입니다.

**모든 PMID는 PubMed E-utilities(esearch + esummary)로 제목·저자·저널·연도를 실제 조회하여
대조 확인했습니다.** 기억에 의존해 적은 PMID는 없으며, 조회로 확인되지 않은 문헌은
목록에서 제외했습니다.

각 섹션 끝의 **▶ 모델 연결** 항목은 그 문헌군이 모델의 어느 방정식·파라미터·구조적 선택을
뒷받침하는지 명시합니다.

---

## 1. 총론 · 역학 · 진료지침 (Overview, Epidemiology, Guidelines)

1. Sheth KN. **Spontaneous Intracerebral Hemorrhage.** *N Engl J Med.* 2022;387:1589-1596.
   <https://pubmed.ncbi.nlm.nih.gov/36300975/>
2. Puy L, Parry-Jones AR, Sandset EC, Dowlatshahi D, Ziai W, Cordonnier C.
   **Intracerebral haemorrhage.** *Nat Rev Dis Primers.* 2023;9:14.
   <https://pubmed.ncbi.nlm.nih.gov/36928219/>
3. Greenberg SM, Ziai WC, Cordonnier C, et al. **2022 Guideline for the Management of
   Patients With Spontaneous Intracerebral Hemorrhage: A Guideline From the American
   Heart Association/American Stroke Association.** *Stroke.* 2022;53:e282-e361.
   <https://pubmed.ncbi.nlm.nih.gov/35579034/>
4. Cordonnier C, Demchuk A, Ziai W, Anderson CS. **Intracerebral haemorrhage: current
   approaches to acute management.** *Lancet.* 2018;392:1257-1268.
   <https://pubmed.ncbi.nlm.nih.gov/30319113/>
5. Magid-Bernstein J, Girard R, Polster S, et al. **Cerebral Hemorrhage: Pathophysiology,
   Treatment, and Future Directions.** *Circ Res.* 2022;130:1204-1229.
   <https://pubmed.ncbi.nlm.nih.gov/35420918/>
6. Keep RF, Hua Y, Xi G. **Intracerebral haemorrhage: mechanisms of injury and therapeutic
   targets.** *Lancet Neurol.* 2012;11:720-731.
   <https://pubmed.ncbi.nlm.nih.gov/22698888/>
7. Xi G, Keep RF, Hoff JT. **Mechanisms of brain injury after intracerebral haemorrhage.**
   *Lancet Neurol.* 2006;5:53-63. <https://pubmed.ncbi.nlm.nih.gov/16361023/>
8. Aronowski J, Zhao X. **Molecular pathophysiology of cerebral hemorrhage: secondary
   brain injury.** *Stroke.* 2011;42:1781-1786.
   <https://pubmed.ncbi.nlm.nih.gov/21527759/>
9. van Asch CJ, Luitse MJ, Rinkel GJ, van der Tweel I, Algra A, Klijn CJ. **Incidence,
   case fatality, and functional outcome of intracerebral haemorrhage over time,
   according to age, sex, and ethnic origin: a systematic review and meta-analysis.**
   *Lancet Neurol.* 2010;9:167-176. <https://pubmed.ncbi.nlm.nih.gov/20056489/>
10. GBD 2019 Stroke Collaborators. **Global, regional, and national burden of stroke and
    its risk factors, 1990-2019.** *Lancet Neurol.* 2021;20:795-820.
    <https://pubmed.ncbi.nlm.nih.gov/34487721/>
11. Qureshi AI, Tuhrim S, Broderick JP, Batjer HH, Hondo H, Hanley DF. **Spontaneous
    intracerebral hemorrhage.** *N Engl J Med.* 2001;344:1450-1460.
    <https://pubmed.ncbi.nlm.nih.gov/11346811/>
12. Meretoja A, Strbian D, Putaala J, et al. **SMASH-U: a proposal for etiologic
    classification of intracerebral hemorrhage.** *Stroke.* 2012;43:2592-2597.
    <https://pubmed.ncbi.nlm.nih.gov/22858729/>
13. Hostettler IC, Seiffge DJ, Werring DJ. **Intracerebral hemorrhage: an update on
    diagnosis and treatment.** *Expert Rev Neurother.* 2019;19:679-694.
    <https://pubmed.ncbi.nlm.nih.gov/31188036/>

**▶ 모델 연결** — 전체 구조의 골격, 클러스터 1(원인·소혈관 병변)의 심부/뇌엽 구분,
`LOC` 파라미터, 사망률 로짓의 기저 발생률(`MORTA`), 그리고 90일 mRS를 1차 결과로 삼는
설계. 문헌 6-8은 "일차 손상(질량효과) + 이차 손상(헴·철·염증)" 이분법의 근거이며
모델의 **구조적 선택 ②**(두 개의 분리 가능한 손상 채널)가 여기서 나옵니다.

---

## 2. 혈종 확장 — 모델의 핵심 적분량 (Haematoma Expansion)

14. Brott T, Broderick J, Kothari R, et al. **Early hemorrhage growth in patients with
    intracerebral hemorrhage.** *Stroke.* 1997;28:1-5.
    <https://pubmed.ncbi.nlm.nih.gov/8996478/>
15. Davis SM, Broderick J, Hennerici M, et al. **Hematoma growth is a determinant of
    mortality and poor outcome after intracerebral hemorrhage.** *Neurology.*
    2006;66:1175-1181. <https://pubmed.ncbi.nlm.nih.gov/16636233/>
16. Al-Shahi Salman R, Frantzias J, Lee RJ, et al. **Absolute risk and predictors of the
    growth of acute spontaneous intracerebral haemorrhage: a systematic review and
    meta-analysis of individual patient data.** *Lancet Neurol.* 2018;17:885-894.
    <https://pubmed.ncbi.nlm.nih.gov/30120039/>
17. Morotti A, Boulouis G, Dowlatshahi D, et al. **Intracerebral haemorrhage expansion:
    definitions, predictors, and prevention.** *Lancet Neurol.* 2023;22:159-171.
    <https://pubmed.ncbi.nlm.nih.gov/36309041/>
18. Kothari RU, Brott T, Broderick JP, et al. **The ABCs of measuring intracerebral
    hemorrhage volumes.** *Stroke.* 1996;27:1304-1305.
    <https://pubmed.ncbi.nlm.nih.gov/8711791/>
19. Wada R, Aviv RI, Fox AJ, et al. **CT angiography "spot sign" predicts hematoma
    expansion in acute intracerebral hemorrhage.** *Stroke.* 2007;38:1257-1262.
    <https://pubmed.ncbi.nlm.nih.gov/17322083/>
20. Demchuk AM, Dowlatshahi D, Rodriguez-Luna D, et al. **Prediction of haematoma growth
    and outcome in patients with intracerebral haemorrhage using the CT-angiography spot
    sign (PREDICT): a prospective observational study.** *Lancet Neurol.* 2012;11:307-314.
    <https://pubmed.ncbi.nlm.nih.gov/22405630/>
21. Delgado Almandoz JE, Yoo AJ, Stone MJ, et al. **The spot sign score in primary
    intracerebral hemorrhage identifies patients at highest risk of in-hospital mortality
    and poor outcome among survivors.** *Stroke.* 2010;41:54-60.
    <https://pubmed.ncbi.nlm.nih.gov/19910545/>
22. Li Q, Zhang G, Huang YJ, et al. **Blend Sign on Computed Tomography: Novel and
    Reliable Predictor for Early Hematoma Growth in Patients With Intracerebral
    Hemorrhage.** *Stroke.* 2015;46:2119-2123.
    <https://pubmed.ncbi.nlm.nih.gov/26089330/>
23. Li Q, Zhang G, Xiong X, et al. **Black Hole Sign: Novel Imaging Marker That Predicts
    Hematoma Growth in Patients With Intracerebral Hemorrhage.** *Stroke.* 2016;47:1777-1781.
    <https://pubmed.ncbi.nlm.nih.gov/27174523/>
24. Blacquiere D, Demchuk AM, Al-Hazzaa M, et al. **Intracerebral Hematoma Morphologic
    Appearance on Noncontrast Computed Tomography Predicts Significant Hematoma
    Expansion.** *Stroke.* 2015;46:3111-3116.
    <https://pubmed.ncbi.nlm.nih.gov/26451019/>

**▶ 모델 연결** — 이 문헌군이 **구조적 선택 ①**의 정량적 표적입니다. 문헌 14-16은
무치료 확장률 20-38%와 확장이 대부분 첫 수 시간에 일어난다는 사실을 제공하며, 이것이
`KBLEED`(0.085)와 `KSEAL`(0.10, 출혈 창을 6-12시간으로 늘림)의 보정 근거입니다.
초안에서 `KSEAL = 0.80`이었을 때 출혈이 1.5시간 내에 종료되어 강압제가 도달하기 전에
끝나버렸고, BP 효과가 문헌보다 훨씬 작게 나왔습니다 — 이 문헌군 때문에 모델을
고쳤습니다. 문헌 15는 `dV24 -> NIHSS -> mRS` 경로의 근거이고, 문헌 19-24는 지도의
CTA spot sign / blend sign / black hole sign 노드와 `SpotSign` 상태에 대응합니다.

---

## 3. 혈압 조절 — 부호가 반대인 두 경로 (Blood Pressure Management)

25. Anderson CS, Heeley E, Huang Y, et al. **Rapid blood-pressure lowering in patients
    with acute intracerebral hemorrhage (INTERACT2).** *N Engl J Med.* 2013;368:2355-2365.
    <https://pubmed.ncbi.nlm.nih.gov/23713578/>
26. Qureshi AI, Palesch YY, Barsan WG, et al. **Intensive Blood-Pressure Lowering in
    Patients with Acute Cerebral Hemorrhage (ATACH-2).** *N Engl J Med.* 2016;375:1033-1043.
    <https://pubmed.ncbi.nlm.nih.gov/27276234/>
27. Ma L, Hu X, Song L, et al. **The third Intensive Care Bundle with Blood Pressure
    Reduction in Acute Cerebral Haemorrhage Trial (INTERACT3): an international,
    stepped wedge cluster randomised controlled trial.** *Lancet.* 2023;402:27-40.
    <https://pubmed.ncbi.nlm.nih.gov/37245517/>
28. Li G, Lin Y, Yang J, et al. **Intensive Ambulance-Delivered Blood-Pressure Reduction
    in Hyperacute Stroke (INTERACT4).** *N Engl J Med.* 2024;390:1862-1872.
    <https://pubmed.ncbi.nlm.nih.gov/38752650/>
29. Anderson CS, Huang Y, Wang JG, et al. **Intensive blood pressure reduction in acute
    cerebral haemorrhage trial (INTERACT): a randomised pilot trial.** *Lancet Neurol.*
    2008;7:391-399. <https://pubmed.ncbi.nlm.nih.gov/18396107/>
30. Moullaali TJ, Wang X, Martin RH, et al. **Blood pressure control and clinical outcomes
    in acute intracerebral haemorrhage: a preplanned pooled analysis of individual
    participant data.** *Lancet Neurol.* 2019;18:857-864.
    <https://pubmed.ncbi.nlm.nih.gov/31397290/>
31. Manning L, Hirakawa Y, Arima H, et al. **Blood pressure variability and outcome after
    acute intracerebral haemorrhage: a post-hoc analysis of INTERACT2.** *Lancet Neurol.*
    2014;13:364-373. <https://pubmed.ncbi.nlm.nih.gov/24530176/>
32. Qureshi AI, Huang W, Lobanova I, et al. **Outcomes of Intensive Systolic Blood
    Pressure Reduction in Patients With Intracerebral Hemorrhage and Excessively High
    Initial Systolic Blood Pressure: Post Hoc Analysis of a Randomized Clinical Trial.**
    *JAMA Neurol.* 2020;77:1355-1365. <https://pubmed.ncbi.nlm.nih.gov/32897310/>
33. Butcher K, Jeerakathil T, Emery D, et al. **The Intracerebral Haemorrhage Acutely
    Decreasing Arterial Pressure Trial: ICH ADAPT.** *Int J Stroke.* 2010;5:227-233.
    <https://pubmed.ncbi.nlm.nih.gov/20536619/>
34. Gould B, McCourt R, Asdaghi N, et al. **Autoregulation of cerebral blood flow is
    preserved in primary intracerebral hemorrhage.** *Stroke.* 2013;44:1726-1728.
    <https://pubmed.ncbi.nlm.nih.gov/23619129/>
35. Toyoda K, Koga M, Yamamoto H, et al. **Intensive blood pressure lowering with
    nicardipine and outcomes after intracerebral hemorrhage: An individual participant
    data systematic review.** *Int J Stroke.* 2022;17:494-505.
    <https://pubmed.ncbi.nlm.nih.gov/34542358/>

**▶ 모델 연결** — **구조적 선택 ③**의 근거이자 반증 자료입니다. 문헌 25-26이 제공하는
효과 크기(INTERACT2 평균 성장 4.5 vs 5.3 mL, ATACH-2 확장 18.9% vs 24.4%)가
`KBLEED`·`EMXNIC`·`EC5NIC` 보정의 상한을 정합니다. 초안에서 모델이 40% 감소를 내도록
맞춰져 있었는데, 그것은 **시험이 얻지 못한 결과를 재현하는 것**이므로 목표를 문헌 쪽으로
되돌렸습니다. 문헌 32는 `TCPP`(CPP<60 누적시간)·신장 항과 시나리오 4(과도 강압 = 최소
출혈 + 최악 결과)의 근거이며, 문헌 33-34는 `AUTOR`·`CBFREG`·`RSHIFT`(고혈압성 자동조절
우측 이동)의 근거입니다. 문헌 31은 `SBPV` 노드에 대응합니다. 문헌 27은 시나리오 6
(care bundle)의 곱셈적 구조를, 문헌 28은 지도의 `Prehospital` 노드를 지지합니다.

---

## 4. 지혈 · 항응고 역전 (Haemostatic Therapy & Anticoagulation Reversal)

36. Mayer SA, Brun NC, Begtrup K, et al. **Efficacy and safety of recombinant activated
    factor VII for acute intracerebral hemorrhage (FAST).** *N Engl J Med.*
    2008;358:2127-2137. <https://pubmed.ncbi.nlm.nih.gov/18480205/>
37. Sprigg N, Flaherty K, Appleton JP, et al. **Tranexamic acid for hyperacute primary
    IntraCerebral Haemorrhage (TICH-2): an international randomised, placebo-controlled,
    phase 3 superiority trial.** *Lancet.* 2018;391:2107-2115.
    <https://pubmed.ncbi.nlm.nih.gov/29778325/>
38. Baharoglu MI, Cordonnier C, Al-Shahi Salman R, et al. **Platelet transfusion versus
    standard care after acute stroke due to spontaneous cerebral haemorrhage associated
    with antiplatelet therapy (PATCH): a randomised, open-label, phase 3 trial.**
    *Lancet.* 2016;387:2605-2613. <https://pubmed.ncbi.nlm.nih.gov/27178479/>
39. Steiner T, Poli S, Griebe M, et al. **Fresh frozen plasma versus prothrombin complex
    concentrate in patients with intracranial haemorrhage related to vitamin K
    antagonists (INCH): a randomised trial.** *Lancet Neurol.* 2016;15:566-573.
    <https://pubmed.ncbi.nlm.nih.gov/27302126/>
40. Connolly SJ, Sharma M, Cohen AT, et al. **Andexanet for Factor Xa Inhibitor-Associated
    Acute Intracerebral Hemorrhage (ANNEXA-I).** *N Engl J Med.* 2024;390:1745-1755.
    <https://pubmed.ncbi.nlm.nih.gov/38749032/>
41. Connolly SJ, Crowther M, Eikelboom JW, et al. **Full Study Report of Andexanet Alfa
    for Bleeding Associated with Factor Xa Inhibitors (ANNEXA-4).** *N Engl J Med.*
    2019;380:1326-1335. <https://pubmed.ncbi.nlm.nih.gov/30730782/>
42. Pollack CV Jr, Reilly PA, van Ryn J, et al. **Idarucizumab for Dabigatran Reversal —
    Full Cohort Analysis (RE-VERSE AD).** *N Engl J Med.* 2017;377:431-441.
    <https://pubmed.ncbi.nlm.nih.gov/28693366/>
43. Siegal DM, Curnutte JT, Connolly SJ, et al. **Andexanet Alfa for the Reversal of
    Factor Xa Inhibitor Activity.** *N Engl J Med.* 2015;373:2413-2424.
    <https://pubmed.ncbi.nlm.nih.gov/26559317/>
44. Flaherty ML, Tao H, Haverbusch M, et al. **Warfarin use leads to larger intracerebral
    hematomas.** *Neurology.* 2008;71:1084-1089.
    <https://pubmed.ncbi.nlm.nih.gov/18824672/>
45. Kuramatsu JB, Gerner ST, Schellinger PD, et al. **Anticoagulant reversal, blood
    pressure levels, and anticoagulant resumption in patients with anticoagulation-related
    intracerebral hemorrhage.** *JAMA.* 2015;313:824-836.
    <https://pubmed.ncbi.nlm.nih.gov/25710659/>
46. Eddleston M, de la Torre JC, Oldstone MB, Loskutoff DJ, Edgington TS, Mackman N.
    **Astrocytes are the primary source of tissue factor in the murine central nervous
    system.** *J Clin Invest.* 1993;92:349-358.
    <https://pubmed.ncbi.nlm.nih.gov/8326003/>

**▶ 모델 연결** — 클러스터 3-5 전체와 `CLOT`·`FII`·`PLAS`·`PLTF` 방정식.
문헌 46은 뇌가 조직인자(TF)를 매우 풍부하게 발현한다는 근거로 `TF_brain -> FVIIa_TF`
경로와 `KFORM`이 큰 값(3.0)을 갖는 이유입니다. 문헌 44-45는 `FII0 = 20%`(INR 3.0)
시나리오와 OAC-ICH의 더 길고 큰 확장을 지지하지만, **정량적으로는 약 2배**입니다 —
초안에서 단일 곱셈항(`THRGEN*PLTF`)을 쓰자 3.5배가 나왔고, 그래서 `CLOT` 방정식을
혈소판 팔(`WPLT`)과 응고 팔(`WFIB`)로 분리했습니다. 이 분리는 문헌 38(PATCH: 항혈소판
관련 ICH에서 혈소판 수혈이 **해로움**)과도 정합하며, 모델에서 혈소판 수혈·DDAVP가
VKA 출혈에는 아무 효과가 없는 이유입니다. 문헌 40·43은 `KBINDA`·`KOFFA`·`KELCPA`와
안덱사네트 중단 후 **반동**(모델에서 별도 항 없이 해리+재분포로 창발)의 근거이고,
문헌 36은 rFVIIa가 확장은 줄였으나 결과는 개선하지 못한(그리고 혈전증은 늘린) 사례로
**단일 인자 개입의 한계**를 보여주는 반례입니다.

---

## 5. 두개내압 · 관류 · 자동조절 (Intracranial Mechanics, Perfusion, Autoregulation)

47. Mokri B. **The Monro-Kellie hypothesis: applications in CSF volume depletion.**
    *Neurology.* 2001;56:1746-1748. <https://pubmed.ncbi.nlm.nih.gov/11425944/>
48. Marmarou A, Shulman K, LaMorgese J. **Compartmental analysis of compliance and outflow
    resistance of the cerebrospinal fluid system.** *J Neurosurg.* 1975;43:523-534.
    <https://pubmed.ncbi.nlm.nih.gov/1181384/>
49. Davson H, Hollingsworth G, Segal MB. **The mechanism of drainage of the cerebrospinal
    fluid.** *Brain.* 1970;93:665-678. <https://pubmed.ncbi.nlm.nih.gov/5490270/>
50. Zazulia AR, Diringer MN, Videen TO, et al. **Hypoperfusion without ischemia
    surrounding acute intracerebral hemorrhage.** *J Cereb Blood Flow Metab.*
    2001;21:804-810. <https://pubmed.ncbi.nlm.nih.gov/11435792/>
51. Oeinck M, Neunhoeffer F, Buttler KJ, et al. **Dynamic cerebral autoregulation in acute
    intracerebral hemorrhage.** *Stroke.* 2013;44:2722-2728.
    <https://pubmed.ncbi.nlm.nih.gov/23943213/>
52. Reinhard M, Neunhoeffer F, Gerds TA, et al. **Secondary decline of cerebral
    autoregulation is associated with worse outcome after intracerebral hemorrhage.**
    *Intensive Care Med.* 2010;36:264-271.
    <https://pubmed.ncbi.nlm.nih.gov/19838669/>
53. Czosnyka M, Smielewski P, Kirkpatrick P, Laing RJ, Menon D, Pickard JD. **Continuous
    assessment of the cerebral vasomotor reactivity in head injury.** *Neurosurgery.*
    1997;41:11-17. <https://pubmed.ncbi.nlm.nih.gov/9218290/>

**▶ 모델 연결** — 문헌 47이 `MonroKellie` 항등식과 `VADD`를 **상태가 아니라 용적 수지의
잔차**로 다루는 설계의 근거입니다. 문헌 48-49(Davson 식: 흡수 = (ICP − Pss) × 전도도)가
`KCSFC = 4.2`, `PSS = 5`의 출처이며, 이 두 값이 무손상 상태에서 생산(21 mL/h)과 정확히
평형을 이룹니다 — 초안에서는 이 평형이 없어 **CSF가 음수로 발산**했습니다. `VCSFMIN`
(18 mL의 압출 가능 CSF)은 압력-용적 곡선을 지수함수가 아니라 **이상(biphasic)**으로
만들며, 이는 대상성 소진의 실제 생리입니다. 문헌 50이 `CBFISC = 0.70`의 근거입니다 —
혈종주변 조직은 이미 과소관류(oligaemia) 상태이므로 손상 개시 문턱이 고전적 경색
문턱보다 **높은** 곳에 있어야 하고, `CBFISC = 0.55`였을 때는 관류 채널이 전혀 작동하지
않아 U자 곡선이 창발하지 못했습니다. 문헌 51-53은 `AUTOR`의 동적 저하와 `PRx` 노드,
그리고 자동조절 소실이 결과를 악화시킨다는 `KISCH` 항의 근거입니다.

---

## 6. 혈종주변 부종 (Perihaematomal Oedema)

54. Wagner KR, Xi G, Hua Y, et al. **Lobar intracerebral hemorrhage model in pigs: rapid
    edema development in perihematomal white matter.** *Stroke.* 1996;27:490-497.
    <https://pubmed.ncbi.nlm.nih.gov/8610319/>
55. Xi G, Wagner KR, Keep RF, et al. **Role of blood clot formation on early edema
    development after experimental intracerebral hemorrhage.** *Stroke.* 1998;29:2580-2586.
    <https://pubmed.ncbi.nlm.nih.gov/9836771/>
56. Hua Y, Keep RF, Hoff JT, Xi G. **Thrombin preconditioning attenuates brain edema
    induced by erythrocytes and iron.** *J Cereb Blood Flow Metab.* 2003;23:1448-1454.
    <https://pubmed.ncbi.nlm.nih.gov/14663340/>
57. Urday S, Kimberly WT, Beslow LA, et al. **Targeting secondary injury in intracerebral
    haemorrhage — perihaematomal oedema.** *Nat Rev Neurol.* 2015;11:111-122.
    <https://pubmed.ncbi.nlm.nih.gov/25623787/>
58. Volbers B, Giede-Jeppe A, Gerner ST, et al. **Peak perihemorrhagic edema correlates
    with functional outcome in intracerebral hemorrhage.** *Neurology.* 2018;90:e1005-e1012.
    <https://pubmed.ncbi.nlm.nih.gov/29453243/>
59. Simard JM, Chen M, Tarasov KV, Bhatta S, Ivanova S, Melnitchenko L, Tsymbalyuk N,
    West GA, Gerzanich V. **Newly expressed SUR1-regulated NC(Ca-ATP) channel mediates
    cerebral edema after ischemic stroke.** *Nat Med.* 2006;12:433-440.
    <https://pubmed.ncbi.nlm.nih.gov/16550187/>
60. Jiang B, Li L, Chen Q, et al. **Glibenclamide Attenuates Neuroinflammation and
    Promotes Neurological Recovery After Intracerebral Hemorrhage in Aged Rats.**
    *Front Aging Neurosci.* 2021;13:729652.
    <https://pubmed.ncbi.nlm.nih.gov/34512312/>
61. Chen Y, Chen S, Chang J, Wei J, Feng M, Wang R. **Perihematomal Edema After
    Intracerebral Hemorrhage: An Update on Pathogenesis, Risk Factors, and Therapeutic
    Advances.** *Front Immunol.* 2021;12:740632.
    <https://pubmed.ncbi.nlm.nih.gov/34737745/>
62. Gong Y, Xi G, Hu H, et al. **Complement inhibition attenuates brain edema and
    neurological deficits induced by thrombin.** *Acta Neurochir Suppl.* 2005;95:389-392.
    <https://pubmed.ncbi.nlm.nih.gov/16463887/>

**▶ 모델 연결** — `OEDE`(이온성·혈관성, 3-5일 정점)와 `OEDL`(세포독성·철 구동, 7-14일)의
**두 구획 분리**가 문헌 54-57에서 나옵니다. 문헌 55-56·62는 `THRT -> PAR1 -> BBBP -> OEDE`
경로와 `KTHR`·`KTHRI`의 근거이며, 문헌 59-60은 `SUR1_TRPM4` 노드와 `KSUR1` 항입니다.
문헌 58은 부종 정점이 혈종 용적과 **독립적으로** 결과를 예측한다는 근거로, `PHE_rel`
노드와 `STRAIN`이 `VHEM`만이 아니라 `OEDE + OEDL`을 함께 포함하는 이유입니다.
모델 검증에서 총 PHE 정점은 26.2 mL(4.3일)로 문헌의 3-5일 정점과 일치합니다.

---

## 7. 헴 · 철 · 페롭토시스 — 지연성 화학 채널 (Haem, Iron, Ferroptosis)

63. Wu J, Hua Y, Keep RF, Nakamura T, Hoff JT, Xi G. **Iron and iron-handling proteins in
    the brain after intracerebral hemorrhage.** *Stroke.* 2003;34:2964-2969.
    <https://pubmed.ncbi.nlm.nih.gov/14615611/>
64. Pérez de la Ossa N, Sobrino T, Silva Y, et al. **Iron-related brain damage in patients
    with intracerebral hemorrhage.** *Stroke.* 2010;41:810-813.
    <https://pubmed.ncbi.nlm.nih.gov/20185788/>
65. Hua Y, Keep RF, Hoff JT, Xi G. **Deferoxamine therapy for intracerebral hemorrhage.**
    *Acta Neurochir Suppl.* 2008;105:3-6.
    <https://pubmed.ncbi.nlm.nih.gov/19066072/>
66. Hu S, Hua Y, Keep RF, Feng H, Xi G. **Deferoxamine therapy reduces brain hemin
    accumulation after intracerebral hemorrhage in piglets.** *Exp Neurol.* 2019;318:244-250.
    <https://pubmed.ncbi.nlm.nih.gov/31078524/>
67. Selim M, Foster LD, Moy CS, et al. **Deferoxamine mesylate in patients with
    intracerebral haemorrhage (i-DEF): a multicentre, randomised, placebo-controlled,
    double-blind phase 2 trial.** *Lancet Neurol.* 2019;18:428-438.
    <https://pubmed.ncbi.nlm.nih.gov/30898550/>
68. Selim M, Yeatts S, Goldstein JN, et al. **Safety and tolerability of deferoxamine
    mesylate in patients with acute intracerebral hemorrhage.** *Stroke.* 2011;42:3067-3074.
    <https://pubmed.ncbi.nlm.nih.gov/21868742/>
69. Zeng L, Tan L, Li H, Zhang Q, Li Y, Guo J. **Deferoxamine therapy for intracerebral
    hemorrhage: A systematic review.** *PLoS One.* 2018;13:e0193615.
    <https://pubmed.ncbi.nlm.nih.gov/29566000/>
70. Zille M, Karuppagounder SS, Chen Y, et al. **Neuronal Death After Hemorrhagic Stroke
    In Vitro and In Vivo Shares Features of Ferroptosis and Necroptosis.** *Stroke.*
    2017;48:1033-1043. <https://pubmed.ncbi.nlm.nih.gov/28250197/>
71. Dixon SJ, Lemberg KM, Lamprecht MR, et al. **Ferroptosis: an iron-dependent form of
    nonapoptotic cell death.** *Cell.* 2012;149:1060-1072.
    <https://pubmed.ncbi.nlm.nih.gov/22632970/>
72. Zhao X, Song S, Sun G, et al. **Neuroprotective role of haptoglobin after
    intracerebral hemorrhage.** *J Neurosci.* 2009;29:15819-15827.
    <https://pubmed.ncbi.nlm.nih.gov/20016097/>
73. Chen-Roetling J, Regan RF. **Hemopexin increases the neurotoxicity of hemoglobin when
    haptoglobin is absent.** *J Neurochem.* 2018;145:464-473.
    <https://pubmed.ncbi.nlm.nih.gov/29500821/>
74. Wang Y, Kinzie E, Berger FG, Lim SK, Baumann H. **Haptoglobin, an
    inflammation-inducible plasma protein.** *Redox Rep.* 2001;6:379-385.
    <https://pubmed.ncbi.nlm.nih.gov/11865981/>
75. Zhao X, Sun G, Zhang J, et al. **Transcription factor Nrf2 protects the brain from
    damage produced by intracerebral hemorrhage.** *Stroke.* 2007;38:3280-3286.
    <https://pubmed.ncbi.nlm.nih.gov/17962605/>
76. Summers MR, Jacobs A, Tudway D, Perera P, Ricketts C. **Studies in desferrioxamine and
    ferrioxamine metabolism in normal and iron-loaded subjects.** *Br J Haematol.*
    1979;42:547-555. <https://pubmed.ncbi.nlm.nih.gov/476006/>

**▶ 모델 연결** — **구조적 선택 ②**의 화학 채널 전체:
`RBC_lysis -> HEME -> HO1 -> FEII -> Fenton -> LPO -> ferroptosis -> NEUR`.
문헌 63·64가 `KHEMI`·`YFE`·`FERR`(페리틴 격리)와 혈청 페리틴 노드의 근거이고,
문헌 70-71이 `KFERRO2`(지질과산화 → 뉴런 사멸) 항의 기전적 근거입니다.
문헌 72-74는 `HPG` 파라미터(합토글로빈 유전형별 소거능)와 `KSCAV` 항,
문헌 75는 `Nrf2 -> HO1/FERR/MG2` 경로입니다.
문헌 65-69·76은 데페록사민 PK(`CLDFO`, t½ ≈ 18분)와 `KCHEL`의 근거이며,
**모델에서 DFO는 `VHEM`으로 가는 간선이 하나도 없습니다** — 검증 결과 i-DEF 시나리오는
24시간 혈종 용적이 대조군과 **완전히 동일**하면서 14일 철 노출 AUC만 68.8 → 48.1
(−30%) 감소했고, `KCHEL = 0`으로 두면 그 효과가 전부 사라집니다. 이것이 문헌 67이
보고한 "용적 변화 없는 mRS 이동" 패턴을 모델이 구조적으로 재현하는 방식입니다.

---

## 8. 신경염증 · 혈종 청소 (Neuroinflammation & Haematoma Clearance)

77. Mracsko E, Veltkamp R. **Neuroinflammation after intracerebral hemorrhage.**
    *Front Cell Neurosci.* 2014;8:388. <https://pubmed.ncbi.nlm.nih.gov/25477782/>
78. Ma Q, Chen S, Hu Q, Feng H, Zhang JH, Tang J. **NLRP3 inflammasome contributes to
    inflammation after intracerebral hemorrhage.** *Ann Neurol.* 2014;75:209-219.
    <https://pubmed.ncbi.nlm.nih.gov/24273204/>
79. Zhao X, Sun G, Zhang J, et al. **Hematoma resolution as a target for intracerebral
    hemorrhage treatment: role for peroxisome proliferator-activated receptor gamma in
    microglia/macrophages.** *Ann Neurol.* 2007;61:352-362.
    <https://pubmed.ncbi.nlm.nih.gov/17457822/>
80. Gonzales NR, Shah J, Sangha N, et al. **Design of a prospective, dose-escalation study
    evaluating the Safety of Pioglitazone for Hematoma Resolution in Intracerebral
    Hemorrhage (SHRINC).** *Int J Stroke.* 2013;8:388-396.
    <https://pubmed.ncbi.nlm.nih.gov/22340518/>
81. Chang CF, Massey J, Osherov A, Angenendt da Costa LH, Sansing LH. **Bexarotene
    Enhances Macrophage Erythrophagocytosis and Hematoma Clearance in Experimental
    Intracerebral Hemorrhage.** *Stroke.* 2020;51:612-618.
    <https://pubmed.ncbi.nlm.nih.gov/31826730/>
82. Fu Y, Hao J, Zhang N, et al. **Fingolimod for the treatment of intracerebral
    hemorrhage: a 2-arm proof-of-concept study.** *JAMA Neurol.* 2014;71:1092-1101.
    <https://pubmed.ncbi.nlm.nih.gov/25003359/>

**▶ 모델 연결** — `MG1`(전염증)·`MG2`(수복·탐식)·`NEU`·`IL6`·`IL10` 방정식.
문헌 79·81이 모델에서 혈종 흡수 속도를 상수가 아니라 **`KRESE = KRES0 × MG2`**,
즉 수복형 미세아교세포에 의해 속도 제한되게 만든 근거입니다(따라서 항염증 개입은
청소를 늦출 수도 있습니다). 문헌 78은 `NLRP3 -> SUR1_TRPM4`·`IL6` 경로,
문헌 82는 `Fingolimod` 노드입니다.

**모델링 주의사항 하나** — 검증 과정에서 `IL6`↔`NEU` 되먹임을 선형으로 쓰면
루프 이득이 1을 넘어 사이토카인 상태가 **발산**했습니다(24시간에 IL-6이 18배로 폭주).
두 방향 모두 포화(Michaelis) 함수로 바꾸어야 염증 반응이 문헌대로 **유한하고 스스로
해소되는 펄스**가 됩니다. 이는 미용상의 수정이 아니라 구조적 요건입니다.

---

## 9. 뇌실내출혈 · 수두증 (Intraventricular Haemorrhage & Hydrocephalus)

83. Hanley DF, Lane K, McBee N, et al. **Thrombolytic removal of intraventricular
    haemorrhage in treatment of severe stroke: results of the randomised, multicentre,
    multiregion, placebo-controlled CLEAR III trial.** *Lancet.* 2017;389:603-611.
    <https://pubmed.ncbi.nlm.nih.gov/28081952/>
84. Graeb DA, Robertson WD, Lapointe JS, Nugent RA, Harrison PB. **Computed tomographic
    diagnosis of intraventricular hemorrhage. Etiology and prognosis.** *Radiology.*
    1982;143:91-96. <https://pubmed.ncbi.nlm.nih.gov/6977795/>
85. Hallevi H, Albright KC, Aronowski J, et al. **Intraventricular hemorrhage: Anatomic
    relationships and clinical implications.** *Neurology.* 2008;70:848-852.
    <https://pubmed.ncbi.nlm.nih.gov/18332342/>
86. Roh DJ, Boehme A, Mamoon R, et al. **Intraventricular Hemorrhage Expansion in the
    CLEAR III Trial: A Post Hoc Exploratory Analysis.** *Stroke.* 2022;53:1847-1853.
    <https://pubmed.ncbi.nlm.nih.gov/35086362/>
87. Strahle JM, Garton HJ, Maher CO, Muraszko KM, Keep RF, Xi G. **Role of hemoglobin and
    iron in hydrocephalus after neonatal intraventricular hemorrhage.** *Neurosurgery.*
    2014;75:696-705. <https://pubmed.ncbi.nlm.nih.gov/25121790/>

**▶ 모델 연결** — `VIVH`·`VCSF`·`EVD`·`IVTPA` 방정식. 문헌 83이 `EMXTPA`·`KIVHCL`과
시나리오 12의 근거이며, **CLEAR III의 핵심 비대칭(사망률은 감소하나 mRS 0-3은 거의
불변)** 이 모델에서도 재현됩니다. 문헌 84-85는 `IVH_cast`·Graeb 노드,
문헌 87은 IVH 후 철이 수두증에 기여한다는 `Ependymal`·`ShuntDep` 경로입니다.
CSF 흡수 차단(`1 - VIVH/(KIVHB+VIVH)`)이 생산에는 적용되지 않는 **비대칭**이
IVH → 수두증 기전의 전부입니다.

---

## 10. 수술 전략 (Surgical Strategies)

88. Mendelow AD, Gregson BA, Fernandes HM, et al. **Early surgery versus initial
    conservative treatment in patients with spontaneous supratentorial intracerebral
    haematomas in the International Surgical Trial in Intracerebral Haemorrhage (STICH):
    a randomised trial.** *Lancet.* 2005;365:387-397.
    <https://pubmed.ncbi.nlm.nih.gov/15680453/>
89. Mendelow AD, Gregson BA, Rowan EN, et al. **Early surgery versus initial conservative
    treatment in patients with spontaneous supratentorial lobar intracerebral haematomas
    (STICH II): a randomised trial.** *Lancet.* 2013;382:397-408.
    <https://pubmed.ncbi.nlm.nih.gov/23726393/>
90. Hanley DF, Thompson RE, Rosenblum M, et al. **Efficacy and safety of minimally
    invasive surgery with thrombolysis in intracerebral haemorrhage evacuation
    (MISTIE III): a randomised, controlled, open-label, blinded endpoint phase 3 trial.**
    *Lancet.* 2019;393:1021-1032. <https://pubmed.ncbi.nlm.nih.gov/30739747/>
91. Pradilla G, Ratcliff JJ, Hall AJ, et al. **Trial of Early Minimally Invasive Removal
    of Intracerebral Hemorrhage (ENRICH).** *N Engl J Med.* 2024;390:1277-1289.
    <https://pubmed.ncbi.nlm.nih.gov/38598795/>
92. Beck J, Fung C, Strbian D, et al. **Decompressive craniectomy plus best medical
    treatment versus best medical treatment alone for spontaneous severe deep
    supratentorial intracerebral haemorrhage: a randomised controlled clinical trial**
    (SWITCH). *Lancet.* 2024;403:2395-2404.
    <https://pubmed.ncbi.nlm.nih.gov/38761811/>

**▶ 모델 연결** — 문헌 90이 모델의 가장 중요한 구조적 교훈 중 하나입니다: **효과의
매개변수는 절개 방식이 아니라 `EOT_vol`(치료 종료 시 잔여 용적)** 이며, MISTIE III의
이득은 ≤15 mL를 달성한 부분군에 국한되었습니다. 따라서 모델에서 수술은 `FEVAC`으로
`VHEM`을 줄이는 동시에 `SURGREB`로 `NOPEN`을 늘리며(시술 관련 재출혈),
질량효과 채널에만 작용하고 이미 유리된 철에는 손대지 못합니다. 문헌 88-89(STICH I/II의
중성 결과)는 "용적을 줄였다"만으로는 결과가 개선되지 않는다는 반례로서 `STRTHR`
(문턱 이하 변형에서는 질량 손상이 무시할 만함) 항의 존재 이유입니다. 문헌 92는
`DecompCranx -> VDC`(두개 개방 시 대상성 용적 증가)에 대응합니다.

---

## 11. 이차 예방 · 뇌아밀로이드혈관병 (Secondary Prevention & CAA)

93. SPS3 Study Group. **Blood-pressure targets in patients with recent lacunar stroke:
    the SPS3 randomised trial.** *Lancet.* 2013;382:507-515.
    <https://pubmed.ncbi.nlm.nih.gov/23726159/>
94. PROGRESS Collaborative Group. **Randomised trial of a perindopril-based
    blood-pressure-lowering regimen among 6105 individuals with previous stroke or
    transient ischaemic attack.** *Lancet.* 2001;358:1033-1041.
    <https://pubmed.ncbi.nlm.nih.gov/11589932/>
95. RESTART Collaboration. **Effects of antiplatelet therapy after stroke due to
    intracerebral haemorrhage (RESTART): a randomised, open-label trial.** *Lancet.*
    2019;393:2613-2623. <https://pubmed.ncbi.nlm.nih.gov/31128924/>
96. SoSTART Collaboration. **Effects of oral anticoagulation for atrial fibrillation after
    spontaneous intracranial haemorrhage in the UK: a randomised, open-label, assessor-
    masked, pilot-phase, non-inferiority trial.** *Lancet Neurol.* 2021;20:842-853.
    <https://pubmed.ncbi.nlm.nih.gov/34487722/>
97. Schreuder FHBM, van Nieuwenhuizen KM, Hofmeijer J, et al. **Apixaban versus no
    anticoagulation after anticoagulation-associated intracerebral haemorrhage in patients
    with atrial fibrillation (APACHE-AF): a randomised, open-label, phase 2 trial.**
    *Lancet Neurol.* 2021;20:907-916. <https://pubmed.ncbi.nlm.nih.gov/34687635/>
98. Charidimou A, Boulouis G, Frosch MP, et al. **The Boston criteria version 2.0 for
    cerebral amyloid angiopathy: a multicentre, retrospective, MRI-neuropathology
    diagnostic accuracy study.** *Lancet Neurol.* 2022;21:714-725.
    <https://pubmed.ncbi.nlm.nih.gov/35841910/>
99. Charidimou A, Boulouis G, Gurol ME, et al. **Emerging concepts in sporadic cerebral
    amyloid angiopathy.** *Brain.* 2017;140:1829-1850.
    <https://pubmed.ncbi.nlm.nih.gov/28334869/>
100. Charidimou A, Krishnan A, Werring DJ, Rolf Jäger H. **Cerebral microbleeds: a guide to
     detection and clinical relevance in different disease settings.** *Neuroradiology.*
     2013;55:655-674. <https://pubmed.ncbi.nlm.nih.gov/23708941/>

**▶ 모델 연결** — 지도의 클러스터 17(이차 예방) 전체와 `Recurrence` 노드의
심부(약 2%/년) 대 뇌엽·CAA(약 7%/년) 위험 구분. 문헌 98-100은 `CAA`·`Microbleeds`·
`cSS`·`APOE` 노드와 Boston v2.0 기준, 문헌 93-94는 `BP_longterm`,
문헌 95-97은 `OAC_restart`·`Antiplatelet_restart` 딜레마입니다.
(이 클러스터는 지도에 포함되어 있으나 90일 시뮬레이션 지평 밖이므로 ODE에는
구현되지 않았습니다 — 재발은 모델의 상태변수가 아닙니다.)

---

## 12. 예후 · 바이오마커 (Prognosis & Biomarkers)

101. Hemphill JC 3rd, Bonovich DC, Besmertis L, Manley GT, Johnston SC. **The ICH score: a
     simple, reliable grading scale for intracerebral hemorrhage.** *Stroke.*
     2001;32:891-897. <https://pubmed.ncbi.nlm.nih.gov/11283388/>
102. Sembill JA, Castello JP, Sprügel MI, et al. **Multicenter Validation of the max-ICH
     Score in Intracerebral Hemorrhage.** *Ann Neurol.* 2021;89:474-484.
     <https://pubmed.ncbi.nlm.nih.gov/33222266/>
103. Rost NS, Smith EE, Chang Y, et al. **Prediction of functional outcome in patients with
     primary intracerebral hemorrhage: the FUNC score.** *Stroke.* 2008;39:2304-2309.
     <https://pubmed.ncbi.nlm.nih.gov/18556582/>
104. Becker KJ, Baxter AB, Cohen WA, et al. **Withdrawal of support in intracerebral
     hemorrhage may lead to self-fulfilling prophecies.** *Neurology.* 2001;56:766-772.
     <https://pubmed.ncbi.nlm.nih.gov/11274312/>
105. Foerch C, Niessner M, Back T, et al. **Diagnostic accuracy of plasma glial fibrillary
     acidic protein for differentiating intracerebral hemorrhage and cerebral ischemia in
     patients with symptoms of acute stroke.** *Clin Chem.* 2012;58:237-245.
     <https://pubmed.ncbi.nlm.nih.gov/22125303/>
106. Wang X, Moullaali TJ, Li Q, et al. **Utility-Weighted Modified Rankin Scale Scores for
     the Assessment of Stroke Outcome: Pooled Analysis of 20 000+ Patients.** *Stroke.*
     2020;51:2411-2417. <https://pubmed.ncbi.nlm.nih.gov/32640944/>

**▶ 모델 연결** — `$TABLE`의 `ICHSC`(문헌 101의 원 채점 규칙 그대로: GCS·용적 ≥30 mL·
IVH 유무·후두개와·연령 ≥80), `PMORT` 로짓(`MORTA`·`MORTB`), `UWMRS`(문헌 106의
효용 가중치), `GFAP_ratio` 노드(문헌 105). 문헌 104는 지도의 `DNR_bias` 노드로
명시적으로 포함했습니다 — **모든 ICH 시험의 결과 해석을 교란하는 자기실현적 예언**이며,
모델이 재현하는 "치료 효과"와 실제 임상시험 결과를 비교할 때 반드시 감안해야 합니다.

**보정 시 주의** — `ICHSC`는 설계상 **기저(baseline) 점수**입니다. 초안에서 90일
시점의 상태로 ICH score를 계산해 사망률을 산출했는데, 이는 회복된 환자를 기준으로
채점하는 것이어서 의도가 뒤집힙니다. 현재는 24시간 점수 + 누적 두개내압 노출
(`TICP`)로 조합합니다.

---

## 13. 약동학 · 약력학 · QSP 방법론 (PK/PD & QSP Methodology)

107. Cook E, Clifton GG, Vargas R, et al. **Pharmacokinetics, pharmacodynamics, and
     minimum effective clinical dose of intravenous nicardipine.** *Clin Pharmacol Ther.*
     1990;47:706-718. <https://pubmed.ncbi.nlm.nih.gov/2357865/>
108. Keating GM. **Clevidipine: a review of its use for managing blood pressure in
     perioperative and intensive care settings.** *Drugs.* 2014;74:1947-1960.
     <https://pubmed.ncbi.nlm.nih.gov/25312594/>
109. Li S, Ahmadzia HK, Guo D, et al. **Population pharmacokinetics and pharmacodynamics
     of tranexamic acid in women undergoing caesarean delivery.** *Br J Clin Pharmacol.*
     2021;87:3579-3589. <https://pubmed.ncbi.nlm.nih.gov/33576009/>
110. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and
     Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.**
     *CPT Pharmacometrics Syst Pharmacol.* 2019;8:883-893.
     <https://pubmed.ncbi.nlm.nih.gov/31652028/>
111. EFPIA MID3 Workgroup. **Good Practices in Model-Informed Drug Discovery and
     Development: Practice, Application, and Documentation.** *CPT Pharmacometrics Syst
     Pharmacol.* 2016;5:93-122. <https://pubmed.ncbi.nlm.nih.gov/27069774/>
112. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for Quantitative Systems
     Pharmacology (QSP) Models: Design and Estimation.** *CPT Pharmacometrics Syst
     Pharmacol.* 2017;6:496-498. <https://pubmed.ncbi.nlm.nih.gov/28585415/>

**▶ 모델 연결** — 문헌 107이 `CLNIC`·`VCNIC`·`QNIC`·`VPNIC`·`KE0NIC`의 근거이며,
`EC5NIC`는 5-15 mg/h 임상 용량 범위가 Emax 곡선의 **선형 구간**에 놓이도록 보정했습니다
(초안에서 `EC5NIC = 0.035 mg/L`이면 3 mg/h에서 이미 포화되어 guideline·intensive·
overshoot 세 팔이 약리학적으로 동일해졌습니다). 문헌 108은 클레비디핀의 초단시간
에스터 분해 특성, 문헌 109는 TXA의 2구획 PK와 항섬유소분해 EC50입니다.
문헌 110-112는 mrgsolve 구현 관례와 QSP 모델 문서화·검증 원칙입니다.

---

## 부록: 이 모델이 재현하지 않는 것 (What This Model Deliberately Does Not Reproduce)

QSP 모델의 신뢰성은 재현하는 것뿐 아니라 **재현하지 않기로 한 것**에서도 드러납니다.

| 항목 | 이유 |
|---|---|
| 혈압 강하의 큰 mRS 개선 | INTERACT2·ATACH-2가 얻지 못한 결과입니다. 모델도 작게 나옵니다. |
| rFVIIa·혈소판 수혈의 이득 | FAST·PATCH는 각각 중성·유해였습니다. 모델에서도 그렇습니다. |
| 재발·장기 이차예방 | 90일 지평 밖이므로 지도에만 있고 ODE에는 없습니다. |
| 개별 환자 예측 | 환자 수준 데이터로 적합하지 않았습니다. 검증은 시험 수준 방향·크기까지입니다. |
| CAA vs 고혈압성 병태의 분자적 구분 | `LOC`·`CAA` 노드로 표현했으나 Aβ 동역학은 구현하지 않았습니다. |

---

## 조회 방법 (Verification Method)

```bash
# 이 문서의 모든 PMID는 아래와 동일한 방식으로 제목을 대조했습니다.
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi\
?db=pubmed&id=<PMID>&retmode=json" | python3 -m json.tool
```

기억에 의존한 PMID는 오류율이 높습니다. 실제 조회 과정에서 초기 후보 PMID 중
상당수가 **다른 논문을 가리키는 것으로 확인되어 교체**되었고, 조회로 확정하지 못한
문헌은 이 목록에서 제외했습니다.
