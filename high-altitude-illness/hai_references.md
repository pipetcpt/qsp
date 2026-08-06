# High-Altitude Illness (AMS · HACE · HAPE) — References

고산병 QSP 모델의 근거 문헌입니다. 각 절은 모델의 해당 구성요소에
대응하며, 파라미터 보정에 **직접** 사용된 문헌에는 **calibration** 표시를
붙였습니다.

> **PMID 검증.** 아래 모든 PubMed ID 는 NCBI E-utilities (`esearch` +
> `esummary`) 로 조회하여 **제목·저자·저널·연도를 실제로 확인**한 것입니다.
> 기억에 의존해 적은 PMID 는 하나도 없습니다. 조회 결과가 의도한 논문과
> 다른 것으로 확인된 후보들은 목록에서 제외했습니다.

---

## 1. 고도의 물리학 — 기압과 흡입 산소분압 (The physics of altitude)

모델의 유일한 외생 변수는 `PIO2 = FiO2 × (PB(h) − 47)` 입니다. 나머지는
전부 이것의 결과입니다.

1. West JB. **Prediction of barometric pressures at high altitude with the use
   of model atmospheres.** J Appl Physiol. 1996;81(4):1850.
   <https://pubmed.ncbi.nlm.nih.gov/8904608/> — **calibration** (모델의
   기압식 `PB = exp(6.63268 − 0.1112h − 0.00149h²)`)
2. West JB. **Barometric pressures on Mt. Everest: new data and physiological
   significance.** J Appl Physiol. 1999;86(3):1062.
   <https://pubmed.ncbi.nlm.nih.gov/10066724/> — **calibration** (정상 253 mmHg)
3. West JB, Hackett PH, Maret KH, et al. **Maximal exercise at extreme altitudes
   on Mount Everest.** J Appl Physiol. 1983;55(3):688.
   <https://pubmed.ncbi.nlm.nih.gov/6415008/>
4. Grocott MPW, Martin DS, Levett DZH, et al. **Arterial blood gases and oxygen
   content in climbers on Mount Everest.** N Engl J Med. 2009;360(2):140.
   <https://pubmed.ncbi.nlm.nih.gov/19129527/> — **calibration** (8400 m
   PaO2 24.6 mmHg · PaCO2 13.3 mmHg · SaO2 54 %)
5. West JB. **High-altitude medicine.** Am J Respir Crit Care Med.
   2012;186(12):1229. <https://pubmed.ncbi.nlm.nih.gov/23103737/>
6. Grocott M, Montgomery H, Vercueil A. **High-altitude physiology and
   pathophysiology: implications and relevance for intensive care medicine.**
   Crit Care. 2007;11(1):203. <https://pubmed.ncbi.nlm.nih.gov/17291330/>

## 2. 폐 가스교환 — 폐포기체식 · 확산 · 션트 (Pulmonary gas exchange)

7. Severinghaus JW. **Simple, accurate equations for human blood O2 dissociation
   computations.** J Appl Physiol. 1979;46(3):599.
   <https://pubmed.ncbi.nlm.nih.gov/35496/> — **calibration** (산소해리곡선
   `S = 1/(23400/(P³+150P)+1)`)
8. Piiper J, Scheid P. **Model for capillary-alveolar equilibration with special
   reference to O2 uptake in hypoxia.** Respir Physiol. 1981;46(3):193.
   <https://pubmed.ncbi.nlm.nih.gov/6798659/> — **calibration** (확산제한
   `exp(−D_L/βQ̇)`; 모델의 A–a 경사가 해수면 9 mmHg → 정상 3 mmHg 로
   *좁아지는* 이유)
9. Wagner PD, Gale GE, Moon RE, et al. **Pulmonary gas exchange in humans
   exercising at sea level and simulated altitude.** J Appl Physiol.
   1986;61(1):260. <https://pubmed.ncbi.nlm.nih.gov/3090012/>
10. Sutton JR, Reeves JT, Wagner PD, et al. **Operation Everest II: oxygen
    transport during exercise at extreme simulated altitude.** J Appl Physiol.
    1988;64(4):1309. <https://pubmed.ncbi.nlm.nih.gov/3132445/>
11. Calbet JAL, Boushel R, Rådegran G, et al. **Determinants of maximal oxygen
    uptake in severe acute hypoxia.** Am J Physiol Regul Integr Comp Physiol.
    2003;284(2):R291. <https://pubmed.ncbi.nlm.nih.gov/12388461/>
12. Calbet JAL, Lundby C. **Air to muscle O2 delivery during exercise at
    altitude.** High Alt Med Biol. 2009;10(2):123.
    <https://pubmed.ncbi.nlm.nih.gov/19555296/>
13. Wehrlin JP, Hallén J. **Linear decrease in V̇O2max and performance with
    increasing altitude in endurance athletes.** Eur J Appl Physiol.
    2006;96(4):404. <https://pubmed.ncbi.nlm.nih.gov/16311764/>

## 3. 환기 조절과 환기 순응 (Ventilatory control and acclimatisation)

모델의 방어 1. 화학수용체가 아니라 **산-염기**가 한계를 정한다는 것이 핵심.

14. Dempsey JA, Forster HV. **Mediation of ventilatory adaptations.** Physiol
    Rev. 1982;62(1):262. <https://pubmed.ncbi.nlm.nih.gov/6798585/>
15. Dempsey JA, Smith CA. **Pathophysiology of human ventilatory control.**
    Eur Respir J. 2014;44(2):495.
    <https://pubmed.ncbi.nlm.nih.gov/24925922/>
16. Teppema LJ, Dahan A. **The ventilatory response to hypoxia in mammals:
    mechanisms, measurement, and analysis.** Physiol Rev. 2010;90(2):675.
    <https://pubmed.ncbi.nlm.nih.gov/20393196/> — **calibration** (HVR
    개인차 범위, 말초·중추 상호작용의 곱셈적 형태)
17. Weil JV, Byrne-Quinn E, Sodal IE, et al. **Hypoxic ventilatory drive in
    normal man.** J Clin Invest. 1970;49(6):1061.
    <https://pubmed.ncbi.nlm.nih.gov/5422012/> — **calibration**
    (등탄산 HVR 기울기 → 모델 `GP`)
18. Rahn H, Otis AB. **Man's respiratory response during and after
    acclimatization to high altitude.** Am J Physiol. 1949;157(3):445.
    <https://pubmed.ncbi.nlm.nih.gov/18151752/> — **calibration**
    (고전적 PaO2–PaCO2 고도선)
19. Guluzade NA, et al. **A test of the interaction between central and
    peripheral respiratory chemoreflexes in humans.** J Physiol.
    2023;601(19):4341. <https://pubmed.ncbi.nlm.nih.gov/37566804/>
20. Fan JL, Subudhi AW, Evero O, et al. **AltitudeOmics: enhanced
    cerebrovascular reactivity and ventilatory response to CO2 with
    high-altitude acclimatization and re-exposure.** J Appl Physiol.
    2014;116(7):911. <https://pubmed.ncbi.nlm.nih.gov/24356520/>
21. Forster HV, et al. **Control of breathing during exercise.** Compr Physiol.
    2012;2(1):743. <https://pubmed.ncbi.nlm.nih.gov/23728984/>

## 4. 산-염기 — 순응의 율속 단계 (Acid-base: the rate-limiting step)

22. Krapf R, Beeler I, Hertner D, Hulter HN. **Chronic respiratory alkalosis.
    The effect of sustained hyperventilation on renal regulation of acid-base
    equilibrium.** N Engl J Med. 1991;324(20):1394.
    <https://pubmed.ncbi.nlm.nih.gov/1902283/> — **calibration** (만성
    ΔHCO3⁻/ΔPaCO2 ≈ −0.5, 신장 시간상수 `TAU_REN`)
23. Siggaard-Andersen O. **Oxygen and acid-base parameters of arterial and mixed
    venous blood, relevant versus redundant.** Acta Anaesthesiol Scand Suppl.
    1995;107:21. <https://pubmed.ncbi.nlm.nih.gov/8599280/> — **calibration**
    (Van Slyke 식 → 모델이 급성 완충을 파라미터 없이 재현하는 근거)

## 5. 아세타졸아미드 — 무엇을 하는가 (Acetazolamide: what it actually does)

모델의 결론: 아세타졸아미드는 호흡을 자극하지 않는다. **바닥(무호흡 역치)을
작동점보다 더 빨리 내린다.**

24. Swenson ER. **Carbonic anhydrase inhibitors and ventilation: a complex
    interplay of stimulation and suppression.** Eur Respir J. 1998;12(6):1242.
    <https://pubmed.ncbi.nlm.nih.gov/9877470/>
25. Swenson ER, Teppema LJ. **Prevention of acute mountain sickness by
    acetazolamide: as yet an unfinished story.** J Appl Physiol.
    2007;102(4):1305. <https://pubmed.ncbi.nlm.nih.gov/17194729/>
26. Leaf DE, Goldfarb DS. **Mechanisms of action of acetazolamide in the
    prophylaxis and treatment of acute mountain sickness.** J Appl Physiol.
    2007;102(4):1313. <https://pubmed.ncbi.nlm.nih.gov/17023566/>
27. Swenson ER, Hughes JMB. **Effects of acute and chronic acetazolamide on
    resting ventilation and ventilatory responses in men.** J Appl Physiol.
    1993;74(1):230. <https://pubmed.ncbi.nlm.nih.gov/8444696/> —
    **calibration** (HCO3⁻ 하강폭 → `ACZ_EMAX_REN`)
28. Teppema LJ, Balanos GM, Steinback CD, et al. **Effects of acetazolamide on
    ventilatory, cerebrovascular, and pulmonary vascular responses to hypoxia.**
    Am J Respir Crit Care Med. 2007;175(3):277.
    <https://pubmed.ncbi.nlm.nih.gov/17095745/>
29. Basnyat B, Gertsch JH, Johnson EW, et al. **Efficacy of low-dose
    acetazolamide (125 mg BID) for the prophylaxis of acute mountain sickness.**
    High Alt Med Biol. 2003;4(1):45.
    <https://pubmed.ncbi.nlm.nih.gov/12713711/> — **calibration** (125 mg 용량군)
30. Low EV, Avery AJ, Gupta V, et al. **Identifying the lowest effective dose of
    acetazolamide for the prophylaxis of acute mountain sickness: systematic
    review and meta-analysis.** BMJ. 2012;345:e6779.
    <https://pubmed.ncbi.nlm.nih.gov/23081689/> — **calibration** (RR)
31. Kayser B, Dumont L, Lysakowski C, et al. **Reappraisal of acetazolamide for
    the prevention of acute mountain sickness: a systematic review and
    meta-analysis.** High Alt Med Biol. 2012;13(2):82.
    <https://pubmed.ncbi.nlm.nih.gov/22724610/>
32. van Patot MCT, Leadbetter G, Keyes LE, et al. **Prophylactic low-dose
    acetazolamide reduces the incidence and severity of acute mountain
    sickness.** High Alt Med Biol. 2008;9(4):289.
    <https://pubmed.ncbi.nlm.nih.gov/19115912/>
33. Gertsch JH, Basnyat B, Johnson EW, Onopa J, Holck PS. **Randomised, double
    blind, placebo controlled comparison of ginkgo biloba and acetazolamide for
    prevention of acute mountain sickness.** BMJ. 2004;328(7443):797.
    <https://pubmed.ncbi.nlm.nih.gov/15070635/>
34. Richalet JP, Rivera M, Bouchet P, et al. **Acetazolamide: a treatment for
    chronic mountain sickness.** Am J Respir Crit Care Med. 2005;172(11):1427.
    <https://pubmed.ncbi.nlm.nih.gov/16126936/>

## 6. 수면 중 주기성 호흡 — CO2 예비량 (Periodic breathing and the CO2 reserve)

모델의 `CO2_reserve = PaCO2 − Bc` 가 여기서 나옵니다.

35. Dempsey JA. **Crossing the apnoeic threshold: causes and consequences.**
    Exp Physiol. 2005;90(1):13. <https://pubmed.ncbi.nlm.nih.gov/15572458/>
    — **calibration** (해수면 CO2 예비량 3–6 mmHg)
36. Dempsey JA, Veasey SC, Morgan BJ, O'Donnell CP. **Pathophysiology of sleep
    apnea.** Physiol Rev. 2010;90(1):47.
    <https://pubmed.ncbi.nlm.nih.gov/20086074/>
37. Khoo MCK, Kronauer RE, Strohl KP, Slutsky AS. **Factors inducing periodic
    breathing in humans: a general model.** J Appl Physiol. 1982;53(3):644.
    <https://pubmed.ncbi.nlm.nih.gov/7129986/> — (loop gain 형식론)
38. Bloch KE, Latshang TD, Turk AJ, et al. **Nocturnal periodic breathing during
    acclimatization at very high altitude at Mount Muztagh Ata (7,546 m).**
    Am J Respir Crit Care Med. 2010;182(4):562.
    <https://pubmed.ncbi.nlm.nih.gov/20442435/> — **calibration** (AHI 대 고도;
    ★ 모델이 *맞히지 못하는* 야간 추이의 근거이기도 함 — README 한계 절 참조)
39. Nussbaumer-Ochsner Y, Ursprung J, Siebenmann C, Maggiorini M, Bloch KE.
    **Effect of short-term acclimatization to high altitude on sleep and
    nocturnal breathing.** Sleep. 2012;35(3):419.
    <https://pubmed.ncbi.nlm.nih.gov/22379248/>
40. Fischer R, Lang SM, Leitl M, et al. **Theophylline and acetazolamide reduce
    sleep-disordered breathing at high altitude.** Eur Respir J. 2004;23(1):47.
    <https://pubmed.ncbi.nlm.nih.gov/14738230/> — **calibration**
    (아세타졸아미드의 AHI 감소폭)
41. Burgess A, et al. **Loop gain response to increased cerebral blood flow at
    high altitude.** Sleep Breath. 2024;28(2):873.
    <https://pubmed.ncbi.nlm.nih.gov/38085496/>

## 7. 급성 고산병 — 역학·점수·상승 속도 (AMS: epidemiology and scoring)

42. Roach RC, Hackett PH, Oelz O, et al. **The 2018 Lake Louise Acute Mountain
    Sickness Score.** High Alt Med Biol. 2018;19(1):4.
    <https://pubmed.ncbi.nlm.nih.gov/29583031/> — **calibration**
    (LLS 4개 하위영역, AMS 정의: 두통 ≥1 그리고 총점 ≥3)
43. Hackett PH, Rennie D, Levine HD. **The incidence, importance, and
    prophylaxis of acute mountain sickness.** Lancet. 1976;2(7996):1149.
    <https://pubmed.ncbi.nlm.nih.gov/62991/>
44. Hackett PH, Roach RC. **High-altitude illness.** N Engl J Med.
    2001;345(2):107. <https://pubmed.ncbi.nlm.nih.gov/11450659/>
45. Bärtsch P, Swenson ER. **Acute high-altitude illnesses.** N Engl J Med.
    2013;369(17):1666. <https://pubmed.ncbi.nlm.nih.gov/24152275/>
46. Luks AM, Swenson ER, Bärtsch P. **Acute high-altitude sickness.** Eur Respir
    Rev. 2017;26(143):160096. <https://pubmed.ncbi.nlm.nih.gov/28143879/>
47. Maggiorini M, Bühler B, Walter M, Oelz O. **Prevalence of acute mountain
    sickness in the Swiss Alps.** BMJ. 1990;301(6756):853.
    <https://pubmed.ncbi.nlm.nih.gov/2282425/> — **calibration** (고도별 유병률)
48. Schneider M, Bernasch D, Weymann J, Holle R, Bärtsch P. **Acute mountain
    sickness: influence of susceptibility, preexposure, and ascent rate.**
    Med Sci Sports Exerc. 2002;34(12):1886.
    <https://pubmed.ncbi.nlm.nih.gov/12471292/> — **calibration**
    (상승 속도 스윕의 비교 대상)
49. Bloch KE, Turk AJ, Maggiorini M, et al. **Effect of ascent protocol on acute
    mountain sickness and success at Muztagh Ata, 7546 m.** High Alt Med Biol.
    2009;10(1):25. <https://pubmed.ncbi.nlm.nih.gov/19326598/>
50. Bärtsch P, Bailey DM, Berger MM, Knauth M, Baumgartner RW. **Acute mountain
    sickness: controversies and advances.** High Alt Med Biol. 2004;5(2):110.
    <https://pubmed.ncbi.nlm.nih.gov/15265333/>
51. Imray C, Wright A, Subudhi A, Roach R. **Acute mountain sickness:
    pathophysiology, prevention, and treatment.** Prog Cardiovasc Dis.
    2010;52(6):467. <https://pubmed.ncbi.nlm.nih.gov/20417340/>
52. Canoui-Poitrine F, Veerabudun K, Larmignat P, et al. **Risk prediction score
    for severe high altitude illness: a cohort study.** PLoS One.
    2014;9(7):e100642. <https://pubmed.ncbi.nlm.nih.gov/25068815/>
53. Luks AM, Auerbach PS, Freer L, et al. **Wilderness Medical Society Clinical
    Practice Guidelines for the Prevention and Treatment of Acute Altitude
    Illness: 2019 Update.** Wilderness Environ Med. 2019;30(4S):S3.
    <https://pubmed.ncbi.nlm.nih.gov/31248818/> — **calibration**
    (권장 상승 속도, 예방약 용량)

## 8. 뇌 — 혈류·부종·두개내압 (Cerebral circulation, oedema, ICP)

모델의 방어 3, 그리고 "tight-fit" 가설의 정량적 형태.

54. Wilson MH, Newman S, Imray CH. **The cerebral effects of ascent to high
    altitudes.** Lancet Neurol. 2009;8(2):175.
    <https://pubmed.ncbi.nlm.nih.gov/19161909/>
55. Ross RT. **The random nature of cerebral mountain sickness.** Lancet.
    1985;1(8435):990. <https://pubmed.ncbi.nlm.nih.gov/2859454/> —
    (원조 "tight-fit" 가설 → 모델의 `PVI` 개인차)
56. Marmarou A, Shulman K, Rosende RM. **A nonlinear analysis of the
    cerebrospinal fluid system and intracranial pressure dynamics.**
    J Neurosurg. 1978;48(3):332. <https://pubmed.ncbi.nlm.nih.gov/632857/> —
    **calibration** (압력-용적 지수 PVI, `ICP = ICP₀·10^(ΔV/PVI)`)
57. Lawley JS, Levine BD, Williams MA, et al. **Cerebral spinal fluid dynamics:
    effect of hypoxia and implications for high-altitude illness.** J Appl
    Physiol. 2016;120(2):251. <https://pubmed.ncbi.nlm.nih.gov/26494441/>
58. Hackett PH, Yarnell PR, Hill R, et al. **High-altitude cerebral edema
    evaluated with magnetic resonance imaging: clinical correlation and
    pathophysiology.** JAMA. 1998;280(22):1920.
    <https://pubmed.ncbi.nlm.nih.gov/9851477/> — **calibration**
    (혈관성 부종, 뇌량 팽대부)
59. Kallenberg K, Bailey DM, Christ S, et al. **Magnetic resonance imaging
    evidence of cytotoxic cerebral edema in acute mountain sickness.**
    J Cereb Blood Flow Metab. 2007;27(5):1064.
    <https://pubmed.ncbi.nlm.nih.gov/17024110/>
60. Schoonman GG, Sándor PS, Nirkko AC, et al. **Hypoxia-induced acute mountain
    sickness is associated with intracellular cerebral edema: a 3 T MRI study.**
    J Cereb Blood Flow Metab. 2008;28(1):198.
    <https://pubmed.ncbi.nlm.nih.gov/17519973/>
61. Sagoo RS, Hutchinson CE, Wright A, et al. **Magnetic resonance investigation
    into the mechanisms involved in the development of high-altitude cerebral
    edema.** J Cereb Blood Flow Metab. 2017;37(1):319.
    <https://pubmed.ncbi.nlm.nih.gov/26746867/>
62. Bailey DM, Bärtsch P, Knauth M, Baumgartner RW. **Emerging concepts in acute
    mountain sickness and high-altitude cerebral edema: from the molecular to
    the morphological.** Cell Mol Life Sci. 2009;66(22):3583.
    <https://pubmed.ncbi.nlm.nih.gov/19763397/>
63. Van Osta A, Moraine JJ, Mélot C, et al. **Effects of high altitude exposure
    on cerebral hemodynamics in normal subjects.** Stroke. 2005;36(3):557.
    <https://pubmed.ncbi.nlm.nih.gov/15692117/> — **calibration** (급성 CBF 증가폭)
64. Ainslie PN, Subudhi AW. **Cerebral blood flow at high altitude.** High Alt
    Med Biol. 2014;15(2):133. <https://pubmed.ncbi.nlm.nih.gov/24971767/>
65. Willie CK, Macleod DB, Shaw AD, et al. **Regional brain blood flow in man
    during acute changes in arterial blood gases.** J Physiol. 2012;590(14):3261.
    <https://pubmed.ncbi.nlm.nih.gov/22495584/> — **calibration**
    (CO2 반응성의 **포화** 형태 → 모델이 선형 반응성을 쓰지 않는 이유)
66. Ainslie PN, Ogoh S, Burgess K, et al. **Differential effects of acute
    hypoxia and high altitude on cerebral blood flow velocity and dynamic
    cerebral autoregulation.** J Appl Physiol. 2008;104(2):490.
    <https://pubmed.ncbi.nlm.nih.gov/18048592/>
67. Subudhi AW, Fan JL, Evero O, et al. **AltitudeOmics: effect of ascent and
    acclimatization to 5260 m on regional cerebral oxygen delivery.** Exp
    Physiol. 2014;99(5):772. <https://pubmed.ncbi.nlm.nih.gov/24243839/>
68. Wilson MH, Imray CHE, Hargens AR. **The headache of high altitude and
    microgravity — similarities with clinical syndromes of cerebral venous
    hypertension.** High Alt Med Biol. 2011;12(4):379.
    <https://pubmed.ncbi.nlm.nih.gov/22087727/>
69. Wilson MH, Davagnanam I, Holland G, et al. **Cerebral venous system and
    anatomical predisposition to high-altitude headache.** Ann Neurol.
    2013;73(3):381. <https://pubmed.ncbi.nlm.nih.gov/23444324/>

## 9. 고산 폐부종 — 모세혈관 압력과 응력파괴 (HAPE: capillary pressure and stress failure)

모델의 방어 2가 어떻게 병이 되는가. 핵심 문헌은 62번(Maggiorini 2001):
HAPE 는 염증이 아니라 **압력**에서 시작한다.

70. Maggiorini M, Mélot C, Pierre S, et al. **High-altitude pulmonary edema is
    initially caused by an increase in capillary pressure.** Circulation.
    2001;103(16):2078. <https://pubmed.ncbi.nlm.nih.gov/11319198/> —
    **calibration** (HAPE 모세혈관압 19–26 mmHg 대 대조군 ~13 mmHg →
    `PCAP_CRIT`)
71. West JB, Mathieu-Costello O. **Stress failure of pulmonary capillaries: role
    in lung and heart disease.** Lancet. 1992;340(8822):762.
    <https://pubmed.ncbi.nlm.nih.gov/1356184/> — **calibration**
    (응력파괴의 초선형 압력 의존성 → `J ∝ (Pcap − Pcrit)^1.5`)
72. Elliott AR, Fu Z, Tsukimoto K, Prediletto R, Mathieu-Costello O, West JB.
    **Short-term reversibility of ultrastructural changes in pulmonary
    capillaries caused by stress failure.** J Appl Physiol. 1992;73(3):1150.
    <https://pubmed.ncbi.nlm.nih.gov/1400030/>
73. Bärtsch P, Mairbäurl H, Maggiorini M, Swenson ER. **Physiological aspects of
    high-altitude pulmonary edema.** J Appl Physiol. 2005;98(3):1101.
    <https://pubmed.ncbi.nlm.nih.gov/15703168/>
74. Swenson ER, Bärtsch P. **High-altitude pulmonary edema.** Compr Physiol.
    2012;2(4):2753. <https://pubmed.ncbi.nlm.nih.gov/23720264/>
75. Swenson ER, et al. **Early hours in the development of high-altitude
    pulmonary edema: time course and mechanisms.** J Appl Physiol.
    2020;128(6):1539. <https://pubmed.ncbi.nlm.nih.gov/32213112/> —
    (염증은 결과이지 원인이 아니라는 근거)
76. Dehnert C, Berger MM, Mairbäurl H, Bärtsch P. **High altitude pulmonary
    edema: a pressure-induced leak.** Respir Physiol Neurobiol.
    2007;158(2-3):266. <https://pubmed.ncbi.nlm.nih.gov/17602898/>
77. Hopkins SR, Garg J, Bolar DS, Balouch J, Levin DL. **Pulmonary blood flow
    heterogeneity during hypoxia and high-altitude pulmonary edema.** Am J
    Respir Crit Care Med. 2005;171(1):83.
    <https://pubmed.ncbi.nlm.nih.gov/15486339/> — **calibration**
    (★ 모델의 핵심 파라미터 `a` = HPV 불균일성의 직접 근거)
78. Hultgren HN. **High-altitude pulmonary edema: current concepts.** Annu Rev
    Med. 1996;47:267. <https://pubmed.ncbi.nlm.nih.gov/8712781/>
79. Maggiorini M. **Prevention and treatment of high-altitude pulmonary edema.**
    Prog Cardiovasc Dis. 2010;52(6):500.
    <https://pubmed.ncbi.nlm.nih.gov/20417343/>
80. Grünig E, Mereles D, Hildebrandt W, et al. **Stress Doppler
    echocardiography for identification of susceptibility to high altitude
    pulmonary edema.** J Am Coll Cardiol. 2000;35(4):980.
    <https://pubmed.ncbi.nlm.nih.gov/10732898/> — **calibration**
    (HAPE 감수성자의 과장된 HPV → `LAMBDA_MAX_S`)
81. Sartori C, Allemann Y, Trueb L, et al. **Augmented vasoreactivity in adult
    life associated with perinatal vascular insult.** Lancet.
    1999;353(9171):2205. <https://pubmed.ncbi.nlm.nih.gov/10392986/>
82. Allemann Y, Hutter D, Lipp E, et al. **Patent foramen ovale and
    high-altitude pulmonary edema.** JAMA. 2006;296(24):2954.
    <https://pubmed.ncbi.nlm.nih.gov/17190896/>
83. Duplain H, Sartori C, Lepori M, et al. **Exhaled nitric oxide in
    high-altitude pulmonary edema: role in the regulation of pulmonary vascular
    tone.** Am J Respir Crit Care Med. 2000;162(1):221.
    <https://pubmed.ncbi.nlm.nih.gov/10903245/>
84. Berger MM, Dehnert C, Bailey DM, et al. **Transpulmonary plasma ET-1 and
    nitrite differences in high altitude pulmonary hypertension.** High Alt Med
    Biol. 2009;10(1):17. <https://pubmed.ncbi.nlm.nih.gov/19326597/>

## 10. 저산소성 폐혈관 수축 (Hypoxic pulmonary vasoconstriction)

85. Sylvester JT, Shimoda LA, Aaronson PI, Ward JPT. **Hypoxic pulmonary
    vasoconstriction.** Physiol Rev. 2012;92(1):367.
    <https://pubmed.ncbi.nlm.nih.gov/22298659/> — **calibration**
    (HPV 의 PAO2 의존성 곡선 → `P50_HPV`, `N_HPV`)
86. Weir EK, López-Barneo J, Buckler KJ, Archer SL. **Acute oxygen-sensing
    mechanisms.** N Engl J Med. 2005;353(19):2042.
    <https://pubmed.ncbi.nlm.nih.gov/16282179/>
87. Dorrington KL, Clar C, Young JD, et al. **Time course of the human pulmonary
    vascular response to 8 hours of isocapnic hypoxia.** Am J Physiol.
    1997;273(3):H1126. <https://pubmed.ncbi.nlm.nih.gov/9321798/> —
    **calibration** (느린 HPV 성분 → `TAU_HPV_S`, `HPVS_MAX`)
88. Naeije R, Dedobbeleer C. **Pulmonary hypertension and the right ventricle in
    hypoxia.** Exp Physiol. 2013;98(8):1247.
    <https://pubmed.ncbi.nlm.nih.gov/23625956/>
89. Faoro V, et al. **Pulmonary circulation and gas exchange at exercise in
    Sherpas at high altitude.** J Appl Physiol. 2014;116(7):919.
    <https://pubmed.ncbi.nlm.nih.gov/23869067/>

## 11. 폐포액 제거 — 배수구가 닫히는 두 번째 양성 되먹임 (Alveolar fluid clearance)

90. Mairbäurl H. **Role of alveolar epithelial sodium transport in high altitude
    pulmonary edema (HAPE).** Respir Physiol Neurobiol. 2006;151(2-3):178.
    <https://pubmed.ncbi.nlm.nih.gov/16337225/>
91. Vivona ML, Matthay M, Chabaud MB, Friedlander G, Clerici C. **Hypoxia
    reduces alveolar epithelial sodium and fluid transport in rats: reversal by
    beta-adrenergic agonist treatment.** Am J Respir Cell Mol Biol.
    2001;25(5):554. <https://pubmed.ncbi.nlm.nih.gov/11713096/> —
    **calibration** (★ 모델의 두 번째 양성 되먹임 `AFC_HYP50`)
92. Sartori C, Allemann Y, Duplain H, et al. **Salmeterol for the prevention of
    high-altitude pulmonary edema.** N Engl J Med. 2002;346(21):1631.
    <https://pubmed.ncbi.nlm.nih.gov/12023995/> — **calibration**
    (β2 작용제의 제거율 상승 → `SAL_EMAX_AFC`)

## 12. HAPE 예방·치료 임상시험 (HAPE prophylaxis and treatment trials)

93. Bärtsch P, Maggiorini M, Ritter M, Noti C, Vock P, Oelz O. **Prevention of
    high-altitude pulmonary edema by nifedipine.** N Engl J Med.
    1991;325(18):1284. <https://pubmed.ncbi.nlm.nih.gov/1922223/> —
    **calibration** (HAPE 발생률 63 % → 10 %)
94. Oelz O, Maggiorini M, Ritter M, et al. **Prevention and treatment of high
    altitude pulmonary edema by a calcium channel blocker.** Int J Sports Med.
    1992;13 Suppl 1:S65. <https://pubmed.ncbi.nlm.nih.gov/1483797/>
95. Maggiorini M, Brunner-La Rocca HP, Peth S, et al. **Both tadalafil and
    dexamethasone may reduce the incidence of high-altitude pulmonary edema: a
    randomized trial.** Ann Intern Med. 2006;145(7):497.
    <https://pubmed.ncbi.nlm.nih.gov/17015867/> — **calibration**
    (위약 74 % · 타다라필 14 % · 덱사메타손 29 % → 모델 시나리오 09–11 의
    비교 대상)
96. Scherrer U, Vollenweider L, Delabays A, et al. **Inhaled nitric oxide for
    high-altitude pulmonary edema.** N Engl J Med. 1996;334(10):624.
    <https://pubmed.ncbi.nlm.nih.gov/8592525/>
97. Nieto Estrada VH, Molano Franco D, Medina RD, et al. **Interventions for
    preventing high altitude illness: Part 1. Commonly-used classes of drugs.**
    Cochrane Database Syst Rev. 2017;6:CD009761.
    <https://pubmed.ncbi.nlm.nih.gov/28653390/>

## 13. 덱사메타손 · 소염진통제 · 물리적 처치 (Dexamethasone, NSAIDs, physical measures)

98. Ferrazzini G, Maggiorini M, Kriemler S, Bärtsch P, Oelz O. **Successful
    treatment of acute mountain sickness with dexamethasone.** Br Med J.
    1987;294(6584):1380. <https://pubmed.ncbi.nlm.nih.gov/3109663/>
99. Levine BD, Yoshimura K, Kobayashi T, et al. **Dexamethasone in the treatment
    of acute mountain sickness.** N Engl J Med. 1989;321(25):1707.
    <https://pubmed.ncbi.nlm.nih.gov/2687688/> — **calibration**
    (★ 증상은 좋아지고 SaO2 는 움직이지 않는다는 비대칭의 원 자료)
100. Burns P, et al. **Altitude sickness prevention with ibuprofen relative to
     acetazolamide.** Am J Med. 2019;132(2):247.
     <https://pubmed.ncbi.nlm.nih.gov/30419226/>
101. Wang J, et al. **Comparative effects of pharmacological interventions for
     the prevention of acute mountain sickness: a systematic review and network
     meta-analysis.** Travel Med Infect Dis. 2025;66:102850.
     <https://pubmed.ncbi.nlm.nih.gov/40383249/>
102. Bärtsch P, Merki B, Hofstetter D, Maggiorini M, Kayser B, Oelz O.
     **Treatment of acute mountain sickness by simulated descent: a randomised
     controlled trial.** BMJ. 1993;306(6885):1098.
     <https://pubmed.ncbi.nlm.nih.gov/8495155/> — **calibration**
     (가모우백 = 모의 하산; 모델의 "하산 등가 미터" 계산의 임상적 대응물)
103. Freeman K, Shalit M, Stroh G. **Use of the Gamow Bag by EMT-basic park
     rangers for treatment of high-altitude pulmonary edema and high-altitude
     cerebral edema.** Wilderness Environ Med. 2004;15(3):198.
     <https://pubmed.ncbi.nlm.nih.gov/15473460/>
104. Luks AM, Swenson ER. **Travel to high altitude with pre-existing lung
     disease.** Eur Respir J. 2007;29(4):770.
     <https://pubmed.ncbi.nlm.nih.gov/17400877/>

## 14. 적혈구생성·혈장량·점도 (Erythropoiesis, plasma volume, viscosity)

모델의 결론: 1주차의 Hct 상승은 적혈구가 아니라 **물**이다. 그리고 최적
헤마토크리트는 고도와 무관하게 43 % 다.

105. Eckardt KU, Boutellier U, Kurtz A, et al. **Rate of erythropoietin
     formation in humans in response to acute hypobaric hypoxia.** J Appl
     Physiol. 1989;66(4):1785. <https://pubmed.ncbi.nlm.nih.gov/2732171/> —
     **calibration** (EPO 24–48 h 최고치, ×3–5)
106. Ge RL, Witkowski S, Zhang Y, et al. **Determinants of erythropoietin
     release in response to short-term hypobaric hypoxia.** J Appl Physiol.
     2002;92(6):2361. <https://pubmed.ncbi.nlm.nih.gov/12015348/>
107. Siebenmann C, Robach P, Lundby C. **Regulation of blood volume in
     lowlanders exposed to high altitude.** J Appl Physiol. 2017;123(4):957.
     <https://pubmed.ncbi.nlm.nih.gov/28572493/> — **calibration**
     (혈장량 수축 폭 → `PV_CONTRACT`)
108. Rasmussen P, Siebenmann C, Díaz V, Lundby C. **Red cell volume expansion at
     altitude: a meta-analysis and Monte Carlo simulation.** Med Sci Sports
     Exerc. 2013;45(9):1767. <https://pubmed.ncbi.nlm.nih.gov/23502972/> —
     **calibration** (Hb mass 증가 속도 ≈ +1 %/주)
109. Winslow RM, Monge CC, Statham NJ, et al. **Variability of oxygen affinity
     of blood: human subjects native to high altitude.** J Appl Physiol.
     1981;51(6):1411. <https://pubmed.ncbi.nlm.nih.gov/7319874/>
110. Mairbäurl H. **Red blood cells in sports: effects of exercise and training
     on oxygen supply by red blood cells.** Front Physiol. 2013;4:332.
     <https://pubmed.ncbi.nlm.nih.gov/24273518/>

## 15. 만성 고산병과 집단 적응 (Chronic mountain sickness and population adaptation)

111. Villafuerte FC, Corante N. **Chronic mountain sickness: clinical aspects,
     etiology, management, and treatment.** High Alt Med Biol. 2016;17(2):61.
     <https://pubmed.ncbi.nlm.nih.gov/27218284/> — **calibration**
     (Qinghai 점수, Hb 진단역치)
112. León-Velarde F, Maggiorini M, Reeves JT, et al. **Consensus statement on
     chronic and subacute high altitude diseases.** High Alt Med Biol.
     2005;6(2):147. <https://pubmed.ncbi.nlm.nih.gov/16060849/>
113. Beall CM. **Two routes to functional adaptation: Tibetan and Andean
     high-altitude natives.** Proc Natl Acad Sci USA. 2007;104 Suppl 1:8655.
     <https://pubmed.ncbi.nlm.nih.gov/17494744/> — (★ 모델의 최적
     헤마토크리트 결과가 말하는 것의 인류학적 대응물)
114. Beall CM. **Natural selection on EPAS1 (HIF2α) associated with low
     hemoglobin concentration in Tibetan highlanders.** Proc Natl Acad Sci USA.
     2010;107(25):11459. <https://pubmed.ncbi.nlm.nih.gov/20534544/>
115. Huerta-Sánchez E, Jin X, Asan, et al. **Altitude adaptation in Tibetans
     caused by introgression of Denisovan-like DNA.** Nature. 2014;512(7513):194.
     <https://pubmed.ncbi.nlm.nih.gov/25043035/>
116. Simonson TS, Yang Y, Huff CD, et al. **Genetic evidence for high-altitude
     adaptation in Tibet.** Science. 2010;329(5987):72.
     <https://pubmed.ncbi.nlm.nih.gov/20466884/>
117. Beall CM. **Human adaptability studies at high altitude: research designs
     and major concepts during fifty years of discovery.** Am J Hum Biol.
     2013;25(2):141. <https://pubmed.ncbi.nlm.nih.gov/23349118/>
118. Levett DZH, Radford EJ, Menassa DA, et al. **Acclimatization of skeletal
     muscle mitochondria to high-altitude hypoxia during an ascent of Everest.**
     FASEB J. 2012;26(4):1431. <https://pubmed.ncbi.nlm.nih.gov/22186874/>

## 16. 체액·신경호르몬 — AMS 의 항이뇨 표현형 (Fluid and neurohumoral responses)

119. Bärtsch P, Maggiorini M, Schobersberger W, et al. **Enhanced
     exercise-induced rise of aldosterone and vasopressin preceding mountain
     sickness.** J Appl Physiol. 1991;71(1):136.
     <https://pubmed.ncbi.nlm.nih.gov/1917735/> — **calibration**
     (AMS = 항이뇨 표현형 → 모델의 `ADH`–`FLUID` 경로)
120. Loeppky JA, Icenogle MV, Maes D, et al. **Early fluid retention and severe
     acute mountain sickness.** J Appl Physiol. 2005;98(2):591.
     <https://pubmed.ncbi.nlm.nih.gov/15501929/>
121. Swenson ER, Duncan TB, Goldberg SV, Ramirez G, Ahmad S, Schoene RB.
     **Diuretic effect of acute hypoxia in humans: relationship to hypoxic
     ventilatory responsiveness and renal hormones.** J Appl Physiol.
     1995;78(2):377. <https://pubmed.ncbi.nlm.nih.gov/7759405/> —
     ("잘 순응하는 사람은 소변을 본다" 의 원자료)

---

## 모델 구성요소 ↔ 문헌 대응표 (Component-to-reference map)

| 모델 구성요소 | 방정식/파라미터 | 근거 문헌 |
|---|---|---|
| 기압 | `PB = exp(6.63268 − 0.1112h − 0.00149h²)` | 1, 2 |
| 폐포기체식 | `PAO2 = PIO2 − PaCO2·[FiO2+(1−FiO2)/R]` | 5, 18 |
| 산소해리곡선 | Severinghaus + Bohr + 2,3-DPG | 7, 109 |
| 확산제한 | Piiper–Scheid `exp(−D_L/βQ̇)` | 8, 9, 10 |
| 중추 무호흡 역치 `Bc` | `[HCO3⁻]_CSF/(0.03·10^(pH*−6.1))` | 35, 36 |
| 말초 화학반사 이득 `GP` | 곱셈적 O₂×CO₂ 상호작용 | 16, 17, 19 |
| 환기 순응 `VAH` | 경동맥소체 가소성, τ ≈ 60 h | 14, 16, 20 |
| 신장 산-염기 `TAU_REN` | τ ≈ 34 h, ΔHCO3⁻/ΔPaCO2 ≈ −0.5 | 22, 23 |
| 아세타졸아미드 PD | 신장 CA + 맥락총 CA | 24–28 |
| CO₂ 예비량 | `PaCO2 − Bc` | 35, 38, 40 |
| HPV `λ(PAO2)` | Hill, P50 45 mmHg | 85, 86, 87 |
| HPV 불균일성 `a` | 두 구획 분할, 상한 = 1/(1−a) | 77, 80 |
| 모세혈관 응력파괴 | `J ∝ (Pcap − 19.5)^1.5` | 70, 71, 72 |
| 폐포액 제거 `AFC` | 저산소가 ENaC/Na-K-ATPase 억제 | 90, 91, 92 |
| 뇌 CO₂ 반응성 | **포화형** 시그모이드 | 65, 66 |
| Monro–Kellie / PVI | `ICP = ICP₀·10^(ΔV/PVI)` | 55, 56, 57 |
| 혈관성/세포독성 부종 | VEGF·BBB / Na-K-ATPase 실패 | 58, 59, 60, 61 |
| Lake Louise 점수 | 두통·위장·피로·어지럼 (0–3 각) | 42 |
| 혈장량 수축 | −20 %, τ ≈ 40 h | 107 |
| 적혈구생성 | EPO 24–48 h 최고, Hb mass +1 %/주 | 105, 106, 108 |
| 최적 Hct | `Hct* = 1/(k·γ)`, k = 2.31 | 108, 110, 113 |
| 니페디핀 / 타다라필 / 덱사메타손 | HAPE 예방 3-arm | 93, 95 |
| 살메테롤 | 제거율 단독 상승 | 92 |
| 덱사메타손의 비대칭 | 증상↓, PaO2 불변 | 99 |
| 가모우백 | +105 mmHg = 모의 하산 | 102, 103 |
| 만성 고산병 | Qinghai 점수, Hb 역치 | 111, 112 |

---

## 이 모델이 문헌과 어긋나는 지점 (Where the model misses)

정직하게 적어 둡니다. 자세한 수치는 `README.md` 의 「한계」 절에 있습니다.

1. **야간 주기성 호흡의 시간 추이 (문헌 38).** 모델은 CO₂ 예비량 하나로
   불안정성을 설명하므로 **도착 첫날 밤이 가장 나쁘고** 이후 개선된다고
   예측합니다. Bloch 등의 현장 수면다원검사는 4559 m 에서 AHI 가 1일차보다
   3일차에 *더 높아지는* 것을 보고했습니다. 빠진 기전은 순응에 따라
   상승하는 말초 화학반사 이득(loop gain)이며, 모델은 이것을 환기에는
   반영하지만 안정성 지표에는 반영하지 않았습니다.
2. **8400 m PaCO2 (문헌 4).** 모델 14.6 mmHg 대 관측 13.3 mmHg — 극단
   고도에서 환기를 약 10 % 과소예측합니다.
3. **최적 헤마토크리트 (문헌 113).** 모델은 순환이 점도 부담을 전혀
   흡수하지 않는다는 가정(γ = 1)에서 43 % 를 냅니다. 실제 심박출량은
   조절되므로 이 값은 하한으로 읽어야 하며, README 에 γ 민감도 분석을
   함께 제시했습니다.
