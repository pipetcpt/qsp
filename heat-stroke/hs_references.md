# 열사병 (Heat Stroke) QSP 모델 — 참고문헌

**Heat stroke — exertional (EHS) and classic / non-exertional (NEHS)**

이 문서는 `hs_qsp_model.dot`, `hs_mrgsolve_model.R`, `hs_core.py`에 들어간
모든 구조적 가정과 파라미터의 출처를 정리한다. 각 절의 마지막에 **이 모델이
해당 문헌에서 무엇을 가져왔는지**를 명시했다.

> 인용은 근거의 위치를 가리키는 것이지 모델의 타당성을 보증하지 않는다.
> 본 모델은 교육·연구 목적이며 임상 의사결정에 사용해서는 안 된다.

---

## 1. 정의, 역학, 두 가지 열사병 (Definitions & epidemiology)

1. Bouchama A, Knochel JP. **Heat stroke.** N Engl J Med. 2002;346(25):1978-88.
   <https://pubmed.ncbi.nlm.nih.gov/12075060/>
2. Epstein Y, Yanovich R. **Heatstroke.** N Engl J Med. 2019;380(25):2449-2459.
   <https://pubmed.ncbi.nlm.nih.gov/31216400/>
3. Bouchama A, Abuyassin B, Lehe C, et al. **Classic and exertional heatstroke.**
   Nat Rev Dis Primers. 2022;8(1):8. <https://pubmed.ncbi.nlm.nih.gov/35115565/>
4. Leon LR, Bouchama A. **Heat stroke.** Compr Physiol. 2015;5(2):611-47.
   <https://pubmed.ncbi.nlm.nih.gov/25880507/>
5. Sorensen C, Hess J. **Treatment and Prevention of Heat-Related Illness.**
   N Engl J Med. 2022;387(15):1404-1413.
   <https://pubmed.ncbi.nlm.nih.gov/36239647/>
6. Argaud L, Ferry T, Le QH, et al. **Short- and long-term outcomes of heatstroke
   following the 2003 heat wave in Lyon, France.** Arch Intern Med.
   2007;167(20):2177-83. <https://pubmed.ncbi.nlm.nih.gov/17698677/>
7. Wallace RF, Kriebel D, Punnett L, et al. **Prior heat illness hospitalization
   and risk of early death.** Environ Res. 2007;104(2):290-5.
   <https://pubmed.ncbi.nlm.nih.gov/17335799/>
8. DeMartini JK, Casa DJ, Stearns R, et al. **Effectiveness of cold water
   immersion in the treatment of exertional heat stroke at the Falmouth Road
   Race.** Med Sci Sports Exerc. 2015;47(2):240-5.
   <https://pubmed.ncbi.nlm.nih.gov/24983342/>

> **모델 반영:** 두 표현형(EHS/NEHS)을 별개 질환이 아니라 **같은 열평형
> 방정식의 분자(numerator)가 다른 두 경우**로 구현했다. 사망률 차이
> (EHS 즉시 냉각 시 거의 0% — Falmouth 코호트 274/274 생존 — 대 NEHS
> 30–60%)는 내재적 치명도가 아니라 §3의 시계 속도와 §4의 용량으로 유도된다.

---

## 2. 체온조절 생리와 열평형 방정식 (Thermoregulation & heat balance)

9. Gagge AP, Stolwijk JA, Nishi Y. **An effective temperature scale based on a
   simple model of human physiological regulatory response.** ASHRAE Trans.
   1971;77:247-262.
10. Havenith G, Fiala D. **Thermal indices and thermophysiological modeling for
    heat stress.** Compr Physiol. 2015;6(1):255-302.
    <https://pubmed.ncbi.nlm.nih.gov/26756633/>
11. Cramer MN, Jay O. **Biophysical aspects of human thermoregulation during
    heat stress.** Auton Neurosci. 2016;196:3-13.
    <https://pubmed.ncbi.nlm.nih.gov/26971392/>
12. Cramer MN, Gagnon D, Laitano O, Crandall CG. **Human temperature regulation
    under heat stress in health, disease, and injury.** Physiol Rev.
    2022;102(4):1907-1989. <https://pubmed.ncbi.nlm.nih.gov/35679471/>
13. Kenny GP, Jay O. **Thermometry, calorimetry, and mean body temperature
    during heat stress.** Compr Physiol. 2013;3(4):1689-719.
    <https://pubmed.ncbi.nlm.nih.gov/24265242/>
14. Charkoudian N. **Skin blood flow in adult human thermoregulation.** Mayo
    Clin Proc. 2003;78(5):603-12. <https://pubmed.ncbi.nlm.nih.gov/12744548/>
15. Rowell LB. **Human cardiovascular adjustments to exercise and thermal
    stress.** Physiol Rev. 1974;54(1):75-159.
    <https://pubmed.ncbi.nlm.nih.gov/4587247/>
16. Buck AL. **New equations for computing vapor pressure and enhancement
    factor.** J Appl Meteorol. 1981;20:1527-1532.
17. Stull R. **Wet-bulb temperature from relative humidity and air
    temperature.** J Appl Meteorol Climatol. 2011;50:2267-2269.

> **모델 반영:** 3-노드(코어·근육·피부) Gagge 계열 열평형. 건열교환
> `Q_dry = A(T_sk−T_o)/(R_cl+1/f_cl·h)`, 증발능
> `E_max = A(P_sk−P_a)/(R_e,cl+1/f_cl·h_e)`, Lewis 관계 `h_e = 16.5·h_c`.
> 피부혈류 0.3→7.5 L/min (ref 14), 내장혈류 재분배 (ref 15).
> 포화수증기압은 Buck(16), 습구온도는 Stull(17) 근사.

---

## 3. 보상 가능/불가능 경계와 임계 습구온도 (The fixed point)

18. Lind AR. **A physiological criterion for setting thermal environmental
    limits for everyday work.** J Appl Physiol. 1963;18:51-6.
    <https://pubmed.ncbi.nlm.nih.gov/13930723/>
19. Belding HS, Hatch TF. **Index for evaluating heat stress in terms of
    resulting physiological strains.** Heat Pip Air Cond. 1955;27:129-136.
20. Sherwood SC, Huber M. **An adaptability limit to climate change due to heat
    stress.** Proc Natl Acad Sci USA. 2010;107(21):9552-5.
    <https://pubmed.ncbi.nlm.nih.gov/20439769/>
21. Vecellio DJ, Wolf ST, Cottle RM, Kenney WL. **Evaluating the 35°C wet-bulb
    temperature adaptability threshold for young, healthy subjects
    (PSU HEAT Project).** J Appl Physiol. 2022;132(2):340-345.
    <https://pubmed.ncbi.nlm.nih.gov/34913738/>
22. Wolf ST, Kenney WL. **The dry heat is hotter: the critical environmental
    limits of young adults in humid vs dry heat.** Am J Physiol Regul Integr
    Comp Physiol. 2023;324(5):R601-R608.
    <https://pubmed.ncbi.nlm.nih.gov/36939211/>
23. Foster J, Smallcombe JW, Hodder S, et al. **An advanced empirical model for
    quantifying the impact of heat and climate change on human physical work
    capacity.** Int J Biometeorol. 2021;65(7):1215-1229.
    <https://pubmed.ncbi.nlm.nih.gov/33674931/>
24. Vanos J, Guzman-Echavarria G, Baldwin JW, et al. **A physiological approach
    for assessing human survivability and liveability to heat in a changing
    climate.** Nat Commun. 2023;14(1):7653.
    <https://pubmed.ncbi.nlm.nih.gov/37993427/>

> **모델 반영:** 이 모델의 **1번 결과**. 보상가능성은 시뮬레이션 문제가 아니라
> **열평형 방정식에 고정점(fixed point)이 존재하는가**라는 문제로 정의했다.
> 계산 결과 안정시 임계 습구온도 **35.7 °C** — Sherwood-Huber(20)의 35 °C
> 한계가 파라미터로 넣지 않았는데도 방정식에서 나온다 — 이고, 900 W 작업 시
> **약 21 °C**로 떨어진다. 이는 Vecellio(21)·Wolf(22)의 "젊고 건강한 성인의
> 실제 임계 습구온도는 35 °C보다 훨씬 낮다"는 실측을 재현한다.
> **부수적으로 얻은 반증:** 이 경계 위 모든 조건에서 `E_actual = E_max,환경`
> 이고 발한 능력(비순화 668 W / 순화 1233 W)은 한 번도 제한 요인이 되지
> 않는다. 따라서 이 모델에서 **열순화는 경계를 움직이지 못한다** (평균 이동
> −0.35 °C). 열순화의 이득은 경계가 아니라 과도기(transient)와 심혈관
> 예비력에 있다는 검증 가능한 구조적 주장이 된다.

---

## 4. 열용량(thermal dose): CEM43 (Sapareto-Dewey)

25. Sapareto SA, Dewey WC. **Thermal dose determination in cancer therapy.**
    Int J Radiat Oncol Biol Phys. 1984;10(6):787-800.
    <https://pubmed.ncbi.nlm.nih.gov/6547421/>
26. Dewhirst MW, Viglianti BL, Lora-Michiels M, et al. **Basic principles of
    thermal dosimetry and thermal thresholds for tissue damage from
    hyperthermia.** Int J Hyperthermia. 2003;19(3):267-94.
    <https://pubmed.ncbi.nlm.nih.gov/12745972/>
27. van Rhoon GC, Samaras T, Yarmolenko PS, et al. **CEM43°C thermal dose
    thresholds: a potential guide for magnetic resonance radiofrequency
    exposure levels?** Eur Radiol. 2013;23(8):2215-27.
    <https://pubmed.ncbi.nlm.nih.gov/23553588/>
28. Yarmolenko PS, Moon EJ, Landon C, et al. **Thresholds for thermal damage to
    normal tissues: an update.** Int J Hyperthermia. 2011;27(4):320-43.
    <https://pubmed.ncbi.nlm.nih.gov/21591897/>
29. Wright NT. **On a relationship between the Arrhenius parameters from
    thermal damage studies.** J Biomech Eng. 2003;125(2):300-4.
    <https://pubmed.ncbi.nlm.nih.gov/12751294/>

> **모델 반영:** 이 모델의 **2번 결과**. 손상변수를 최고 체온이 아니라
> **누적 열용량 CEM43 = ∫R^(43−T_c)dt** (R = 0.25 미만 43 °C, 0.50 이상)로
> 잡았다. 이 한 줄이 "시간과 온도는 1:1로 교환되지 않는다"를 강제한다 —
> 40 °C에서 0.016/분, 42 °C에서 0.250/분, 43 °C에서 1.000/분.
> 종양온열치료에서 확립된 이 선량 법칙을 열사병에 적용한 것이 이 모델의
> 핵심 이식(transplant)이며, 냉각 지연과 냉각 방식을 **하나의 축** 위에
> 올려놓는 것을 가능하게 한다.

---

## 5. 냉각: 방식과 속도 (Cooling modality and rate)

30. Casa DJ, McDermott BP, Lee EC, et al. **Cold water immersion: the gold
    standard for exertional heatstroke treatment.** Exerc Sport Sci Rev.
    2007;35(3):141-9. <https://pubmed.ncbi.nlm.nih.gov/17620933/>
31. Casa DJ, Armstrong LE, Kenny GP, et al. **Exertional heat stroke: new
    concepts regarding cause and care.** Curr Sports Med Rep. 2012;11(3):115-23.
    <https://pubmed.ncbi.nlm.nih.gov/22580488/>
32. McDermott BP, Casa DJ, Ganio MS, et al. **Acute whole-body cooling for
    exercise-induced hyperthermia: a systematic review.** J Athl Train.
    2009;44(1):84-93. <https://pubmed.ncbi.nlm.nih.gov/19180223/>
33. Zhang Y, Davis JK, Casa DJ, Bishop PA. **Optimizing cold water immersion for
    exercise-induced hyperthermia: a meta-analysis.** Med Sci Sports Exerc.
    2015;47(11):2464-72. <https://pubmed.ncbi.nlm.nih.gov/25848900/>
34. Proulx CI, Ducharme MB, Kenny GP. **Effect of water temperature on cooling
    efficiency during hyperthermia in humans.** J Appl Physiol.
    2003;94(4):1317-23. <https://pubmed.ncbi.nlm.nih.gov/12626467/>
35. Luhring KE, Butts CL, Smith CR, et al. **Cooling effectiveness of a modified
    cold-water immersion method after exercise-induced hyperthermia.**
    J Athl Train. 2016;51(11):946-951.
    <https://pubmed.ncbi.nlm.nih.gov/27834504/>
36. Hosokawa Y, Adams WM, Belval LN, et al. **Tarp-assisted cooling as a method
    of whole-body cooling in hyperthermic individuals.** Ann Emerg Med.
    2017;69(3):347-352. <https://pubmed.ncbi.nlm.nih.gov/27522596/>
37. Butts CL, McDermott BP, Buening BJ, et al. **Physiologic and perceptual
    responses to cold-shower cooling after exercise-induced hyperthermia.**
    J Athl Train. 2016;51(3):252-7.
    <https://pubmed.ncbi.nlm.nih.gov/26967831/>
38. Gaudio FG, Grissom CK. **Cooling methods in heat stroke.** J Emerg Med.
    2016;50(4):607-16. <https://pubmed.ncbi.nlm.nih.gov/26525947/>
39. Filep EM, Murata Y, Endres BD, et al. **Exertional heat stroke, modality
    cooling rate, and survival outcomes: a systematic review.** Medicina
    (Kaunas). 2020;56(11):589. <https://pubmed.ncbi.nlm.nih.gov/33167418/>
40. Weiner JS, Khogali M. **A physiological body-cooling unit for treatment of
    heat stroke.** Lancet. 1980;1(8167):507-9.
    <https://pubmed.ncbi.nlm.nih.gov/6102236/>
41. Hoedemaekers CW, Ezzahti M, Gerritsen A, van der Hoeven JG. **Comparison of
    cooling methods to induce and maintain normo- and hypothermia in ICU
    patients.** Crit Care. 2007;11(4):R91.
    <https://pubmed.ncbi.nlm.nih.gov/17718920/>
42. Kim F, Olsufka M, Longstreth WT Jr, et al. **Pilot randomized clinical trial
    of prehospital induction of mild hypothermia with cold normal saline.**
    Circulation. 2007;115(24):3064-70.
    <https://pubmed.ncbi.nlm.nih.gov/17548731/>
43. Bernard S, Buist M, Monteiro O, Smith K. **Induced hypothermia using large
    volume, ice-cold intravenous fluid in comatose survivors of out-of-hospital
    cardiac arrest.** Resuscitation. 2003;56(1):9-13.
    <https://pubmed.ncbi.nlm.nih.gov/12505732/>
44. Belval LN, Casa DJ, Adams WM, et al. **Consensus statement — prehospital
    care of exertional heat stroke.** Prehosp Emerg Care. 2018;22(3):392-397.
    <https://pubmed.ncbi.nlm.nih.gov/29336710/>

> **모델 반영:** 각 냉각 방식을 **하나의 집중 전도도 UA (W/K)** 로 표현하고,
> 위 문헌의 발표된 코어 냉각속도에 맞춰 `hs_calibrate.py`에서 수치적으로
> 적합시켰다: 빙수침수 2 °C → UA 23.5 (0.22 °C/분), 냉수침수 14 °C → 27.7
> (0.17), 방수포 냉각 → 20.1 (0.14), 냉수 샤워 → 15.9 (0.10), 증발·대류 →
> 16.2 (0.08), 혈관내 카테터 → 4.6 (0.06), 얼음팩 → 1.6 (0.032),
> 수동 → 0.016. 방식을 **하나의 역가 숫자**로 만든 것이 요점이다 —
> "어느 냉각기"와 "언제 시작하는가"를 같은 축에 올린다.
> 4 °C 정질액 2 L는 전도도가 아니라 **일회성 엔탈피 싱크**로 구현했고
> (관류 가중 분배), 계산된 코어 하강 1.4–1.7 °C는 Kim(42)·Bernard(43)의
> 관측치와 일치한다.
> **모델이 과대예측하는 지점:** 2 °C 대 14 °C 물의 냉각속도 비를 1.55로
> 예측하는데 메타분석 실측은 1.16–1.29다. 원인은 전도도가 선형이고 떨림
> 열생산 상쇄와 불완전 침수를 모델링하지 않았기 때문이며,
> `hs_verification_output.txt`에 이 모델의 가장 노출된 열역학적 예측으로
> 표시해 두었다. 시나리오에는 방식별로 개별 적합한 UA를 쓰므로 보고되는
> 냉각속도 자체는 관측치와 일치한다.

---

## 6. 냉각 지연이 결과를 결정한다 (Time to cooling)

45. Heled Y, Rav-Acha M, Shani Y, et al. **The "golden hour" for heatstroke
    treatment.** Mil Med. 2004;169(3):184-6.
    <https://pubmed.ncbi.nlm.nih.gov/15080235/>
46. Vicario SJ, Okabajue R, Haltom T. **Rapid cooling in classic heatstroke:
    effect on mortality rates.** Am J Emerg Med. 1986;4(5):394-8.
    <https://pubmed.ncbi.nlm.nih.gov/3741548/>
47. Demartini JK, Casa DJ, Belval LN, et al. **Environmental conditions and the
    occurrence of exertional heat illnesses and exertional heat stroke at the
    Falmouth Road Race.** J Athl Train. 2014;49(4):478-85.
    <https://pubmed.ncbi.nlm.nih.gov/24960519/>
48. Rublee C, Dresser C, Giudice C, et al. **Evidence-based heatstroke
    management in the emergency department.** West J Emerg Med.
    2021;22(2):186-195. <https://pubmed.ncbi.nlm.nih.gov/33856299/>
49. Pryor RR, Roth RN, Suyama J, Hostler D. **Exertional heat illness:
    emerging concepts and advances in prehospital care.** Prehosp Disaster Med.
    2015;30(3):297-305. <https://pubmed.ncbi.nlm.nih.gov/25868636/>

> **모델 반영:** 이 모델의 **3번 결과 — 교환율(exchange rate)**.
> 42.0 °C에서 쓰러진 환자에서 증발냉각 → 빙수침수로 **방식을 올리는 것**은
> CEM43 기준 **4.6–6.0분의 지연**과 같은 값이며, 이 교환율은 쓰러진 온도
> (41/42/43 °C)에 거의 무관하다. 반면 절대 위험은 지수적으로 커진다.
> "cool first, transport second"(44, 48)의 정량적 내용이 바로 이것이다:
> 방식은 몇 배 수준, 지연은 상한이 없다.

---

## 7. 위장관 장벽, 내독소혈증, 내장 허혈 (Gut barrier & endotoxaemia)

50. Lambert GP. **Stress-induced gastrointestinal barrier dysfunction and its
    inflammatory effects.** J Anim Sci. 2009;87(14 Suppl):E101-8.
    <https://pubmed.ncbi.nlm.nih.gov/18791134/>
51. Lim CL, Mackinnon LT. **The roles of exercise-induced immune system
    disturbances in the pathology of heat stroke.** Sports Med.
    2006;36(1):39-64. <https://pubmed.ncbi.nlm.nih.gov/16445310/>
52. Hall DM, Buettner GR, Oberley LW, et al. **Mechanisms of circulatory and
    intestinal barrier dysfunction during whole body hyperthermia.** Am J
    Physiol Heart Circ Physiol. 2001;280(2):H509-21.
    <https://pubmed.ncbi.nlm.nih.gov/11158946/>
53. Dokladny K, Zuhl MN, Moseley PL. **Intestinal epithelial barrier function
    and tight junction proteins with heat and exercise.** J Appl Physiol.
    2016;120(6):692-701. <https://pubmed.ncbi.nlm.nih.gov/26359485/>
54. Ogden HB, Child RB, Fallowfield JL, et al. **The gastrointestinal exertional
    heat stroke paradigm: pathophysiology, assessment, severity, aetiology and
    nutritional countermeasures.** Nutrients. 2020;12(2):537.
    <https://pubmed.ncbi.nlm.nih.gov/32093001/>
55. Bouchama A, Parhar RS, el-Yazigi A, et al. **Endotoxemia and release of
    tumor necrosis factor and interleukin 1 alpha in acute heatstroke.**
    J Appl Physiol. 1991;70(6):2640-4.
    <https://pubmed.ncbi.nlm.nih.gov/1885459/>
56. Snipe RMJ, Khoo A, Kitic CM, et al. **The impact of exertional-heat stress
    on gastrointestinal integrity, systemic endotoxin and cytokine profile.**
    Eur J Appl Physiol. 2018;118(2):389-400.
    <https://pubmed.ncbi.nlm.nih.gov/29234915/>

> **모델 반영:** 피부혈관확장(열을 버리는 바로 그 기전)이 내장혈류를
> 1.5 → 0.3 L/min으로 빼앗는 **cutaneous steal**을 명시적으로 구현했다.
> 장세포 손상은 허혈(ISCH²)과 직접 열손상(CEM43) 두 항의 합이고, 문 하나가
> 열리는 사건이 **동시에 출구도 닫는다** — Kupffer 세포의 내독소 청소율이
> 간 기능과 내장 관류 양쪽에 비례하기 때문이다.

---

## 8. 사이토카인, HMGB1, 그리고 확정 스위치 (Commitment switch)

57. Bouchama A, Hammami MM, Haq A, et al. **Evidence for endothelial cell
    activation/injury in heatstroke.** Crit Care Med. 1996;24(7):1173-8.
    <https://pubmed.ncbi.nlm.nih.gov/8674331/>
58. Bouchama A, Roberts G, Al Mohanna F, et al. **Inflammatory, hemostatic, and
    clinical changes in a baboon experimental model for heatstroke.**
    J Appl Physiol. 2005;98(2):697-705.
    <https://pubmed.ncbi.nlm.nih.gov/15475605/>
59. Leon LR, Helwig BG. **Heat stroke: role of the systemic inflammatory
    response.** J Appl Physiol. 2010;109(6):1980-8.
    <https://pubmed.ncbi.nlm.nih.gov/20522730/>
60. Tong H, Tang Y, Chen Y, et al. **HMGB1 activity inhibition alleviating liver
    injury in heatstroke.** J Trauma Acute Care Surg. 2013;74(3):801-7.
    <https://pubmed.ncbi.nlm.nih.gov/23425739/>
61. Tong HS, Tang YQ, Chen Y, et al. **Early elevated HMGB1 level predicting the
    outcome in exertional heatstroke.** J Trauma. 2011;71(4):808-14.
    <https://pubmed.ncbi.nlm.nih.gov/21460740/>
62. Andersson U, Tracey KJ. **HMGB1 is a therapeutic target for sterile
    inflammation and infection.** Annu Rev Immunol. 2011;29:139-62.
    <https://pubmed.ncbi.nlm.nih.gov/21219181/>
63. Scaffidi P, Misteli T, Bianchi ME. **Release of chromatin protein HMGB1 by
    necrotic cells triggers inflammation.** Nature. 2002;418(6894):191-5.
    <https://pubmed.ncbi.nlm.nih.gov/12110890/>
64. Dehbi M, Baturcam E, Eldali A, et al. **Hsp-72, a candidate prognostic
    indicator of heatstroke.** Cell Stress Chaperones. 2010;15(5):593-603.
    <https://pubmed.ncbi.nlm.nih.gov/20205026/>
65. Ferat-Osorio E, Sánchez-Anaya A, Gutiérrez-Mendoza M, et al. **Heat shock
    protein 70 down-regulates the production of TLR-mediated inflammatory
    cytokines.** J Inflamm (Lond). 2014;11:19.
    <https://pubmed.ncbi.nlm.nih.gov/25053921/>

> **모델 반영:** 이 모델의 **4번 결과 — regime 3**. HMGB1을 확정변수로 삼아
> **안장-마디 쌍안정(saddle-node bistability)** 으로 구현했다:
> `dH/dt = drive + A·H³/(H³+K³) − c_eff·H`.
> 고정점은 OFF = 0, **불안정 = 10.91 ng/mL**, ON = 50.3 ng/mL.
> 괴사세포가 HMGB1을 방출하고(63) HMGB1이 RAGE/TLR4를 통해 다시 괴사를
> 만드는(60, 62) 양성 되먹임에, 생산은 포화하고 청소는 선형이라는 두 조건이
> 붙으면 반드시 두 개의 안정 상태가 생긴다.
> **핵심 계산:** 확정을 일으키는 **용량은 7.8 ± 0.2 CEM43으로 냉각 방식에
> 거의 무관**한 반면(퍼짐 0.69), 확정을 일으키는 **지연은 17분(얼음팩)에서
> 42분(빙수침수)까지** 벌어진다. 용량은 환자의 성질이고, 감당할 수 있는
> 지연은 냉각기의 성질이다.
> HSP70은 열용량을 나누는 보호항으로 들어간다(64, 65).

---

## 9. 응고장애와 DIC (Coagulopathy)

66. Bouchama A, Bridey F, Hammami MM, et al. **Activation of coagulation and
    fibrinolysis in heatstroke.** Thromb Haemost. 1996;76(6):909-15.
    <https://pubmed.ncbi.nlm.nih.gov/8972008/>
67. Hifumi T, Kondo Y, Shimizu K, Miyake Y. **Heat stroke.** J Intensive Care.
    2018;6:30. <https://pubmed.ncbi.nlm.nih.gov/29850022/>
68. Iba T, Levy JH, Warkentin TE, et al. **Diagnosis and management of
    sepsis-induced coagulopathy and disseminated intravascular coagulation.**
    J Thromb Haemost. 2019;17(11):1989-1994.
    <https://pubmed.ncbi.nlm.nih.gov/31410983/>
69. Taylor FB Jr, Toh CH, Hoots WK, et al. **Towards definition, clinical and
    laboratory criteria, and a scoring system for disseminated intravascular
    coagulation (ISTH).** Thromb Haemost. 2001;86(5):1327-30.
    <https://pubmed.ncbi.nlm.nih.gov/11816725/>
70. Ogura Y, Sakurai R, Kawakami R, et al. **Heat-stroke-induced DIC and
    recombinant thrombomodulin.** Acute Med Surg. 2018;5(4):327-334.
    <https://pubmed.ncbi.nlm.nih.gov/30338078/>
71. Yokobori S, Koido Y, Shishido H, et al. **Feasibility and safety of
    intravascular temperature management for severe heat stroke (ICE-HEAT).**
    Crit Care Med. 2018;46(7):e670-e676.
    <https://pubmed.ncbi.nlm.nih.gov/29652723/>
72. Johansson PI, Stensballe J, Ostrowski SR. **Shock-induced endotheliopathy
    (SHINE) in acute critical illness.** Crit Care. 2017;21(1):25.
    <https://pubmed.ncbi.nlm.nih.gov/28173843/>

> **모델 반영:** 트롬빈 생성 → 피브리노겐·혈소판·protein C 소비 → D-dimer
> 상승의 소비성 응고장애로 구현. **ISTH overt-DIC 점수와 SOFA·GCS는 전부
> 상태변수가 아니라 출력으로만 계산**한다 — 점수를 적분하면 점수가
> 병태생리를 구동하게 되어 순환논증이 된다.
> 피브리노겐은 소비항과 급성기 합성항이 **부호가 반대**인 유일한 인자다.

---

## 10. 약물: 해열제, 단트롤렌, 스테로이드, 트롬보모듈린

73. Bouchama A, Cafege A, Devol EB, et al. **Ineffectiveness of dantrolene
    sodium in the treatment of heatstroke.** Crit Care Med. 1991;19(2):176-80.
    <https://pubmed.ncbi.nlm.nih.gov/1989755/>
74. Krause T, Gerbershagen MU, Fiege M, et al. **Dantrolene — a review of its
    pharmacology, therapeutic use and new developments.** Anaesthesia.
    2004;59(4):364-73. <https://pubmed.ncbi.nlm.nih.gov/15023108/>
75. Rosenberg H, Pollock N, Schiemann A, et al. **Malignant hyperthermia: a
    review.** Orphanet J Rare Dis. 2015;10:93.
    <https://pubmed.ncbi.nlm.nih.gov/26238698/>
76. Prescott LF. **Paracetamol, alcohol and the liver.** Br J Clin Pharmacol.
    2000;49(4):291-301. <https://pubmed.ncbi.nlm.nih.gov/10759684/>
77. Mitchell JR, Jollow DJ, Potter WZ, et al. **Acetaminophen-induced hepatic
    necrosis. IV. Protective role of glutathione.** J Pharmacol Exp Ther.
    1973;187(1):211-7. <https://pubmed.ncbi.nlm.nih.gov/4746329/>
78. Whitcomb DC, Block GD. **Association of acetaminophen hepatotoxicity with
    fasting and ethanol use.** JAMA. 1994;272(23):1845-50.
    <https://pubmed.ncbi.nlm.nih.gov/7990219/>
79. Saito T, Maruyama I, Shimazaki S, et al. **Efficacy and safety of
    recombinant human soluble thrombomodulin (ART-123) in disseminated
    intravascular coagulation.** J Thromb Haemost. 2007;5(1):31-41.
    <https://pubmed.ncbi.nlm.nih.gov/17059423/>
80. Abeyama K, Stern DM, Ito Y, et al. **The N-terminal domain of
    thrombomodulin sequesters high-mobility group-B1 protein: a novel
    antiinflammatory mechanism.** J Clin Invest. 2005;115(5):1267-74.
    <https://pubmed.ncbi.nlm.nih.gov/15841214/>
81. Vaity C, Al-Subaie N, Cecconi M. **Cooling techniques for hyperthermia in
    the intensive care unit.** Crit Care. 2015;19:103.
    <https://pubmed.ncbi.nlm.nih.gov/25886948/>
82. Niven DJ, Stelfox HT, Laupland KB. **Antipyretic therapy in febrile
    critically ill adults: a systematic review and meta-analysis.**
    J Crit Care. 2013;28(3):303-10.
    <https://pubmed.ncbi.nlm.nih.gov/23159136/>
83. Sund-Levander M, Grodzinsky E. **Assessment of body temperature measurement
    options.** Br J Nurs. 2013;22(15):880-8.
    <https://pubmed.ncbi.nlm.nih.gov/24006812/>

> **모델 반영:** 이 모델의 **5번 결과 — 해열제 논증**. 발열은 **설정점**을
> 올리고 몸이 그것을 방어하지만, 열사병은 설정점이 37 °C 그대로이고 효과기가
> 이미 포화되어 있다. 그래서 해열제를 **역치에만** 작용하도록 구현했고,
> dTc/dt의 설정점 민감도는 보상 가능 영역에서 크고 **비보상 영역에서 0으로
> 붕괴한다** — 해열제는 열사병이 아닌 영역에서만 작동한다(82).
> 비용은 따로 청구된다: CYP2E1이 유도되고 글루타치온이 이미 소모되는
> 간에서의 NAPQI(76-78), NSAID의 신장 프로스타글란딘 차단과 혈소판 기능
> 억제.
> 단트롤렌은 RyR1 기전(74, 75) — 즉 **악성 고열의 기전이지 이것의 기전이
> 아니다** — 이므로 근육 열생산의 ≤12%만 줄이고, Bouchama의 RCT(73)에서
> 관찰된 null을 재현한다.
> 트롬보모듈린 알파는 이 모델에서 **확정변수를 건드리는 유일한 약물**이다:
> lectin-유사 도메인이 HMGB1을 분해하므로(80) HMGB1 청소율에 더해지고,
> 안장-마디 계산상 **+0.00192/분을 넘으면 ON 상태 자체가 사라진다**.
> 따라서 그 치료 창은 체온계가 아니라 스위치가 정한다 — 검증 가능한
> 구조적 예측이다(70, 79).

---

## 11. 열순화, 탈수, 위험인자 (Acclimatisation, hydration, risk)

84. Périard JD, Racinais S, Sawka MN. **Adaptations and mechanisms of human heat
    acclimation.** Scand J Med Sci Sports. 2015;25 Suppl 1:20-38.
    <https://pubmed.ncbi.nlm.nih.gov/25943654/>
85. Périard JD, Eijsvogels TMH, Daanen HAM. **Exercise under heat stress:
    thermoregulation, hydration, performance implications, and mitigation
    strategies.** Physiol Rev. 2021;101(4):1873-1979.
    <https://pubmed.ncbi.nlm.nih.gov/33829868/>
86. Sawka MN, Burke LM, Eichner ER, et al. **ACSM position stand: exercise and
    fluid replacement.** Med Sci Sports Exerc. 2007;39(2):377-90.
    <https://pubmed.ncbi.nlm.nih.gov/17277604/>
87. Cheuvront SN, Kenefick RW. **Dehydration: physiology, assessment, and
    performance effects.** Compr Physiol. 2014;4(1):257-85.
    <https://pubmed.ncbi.nlm.nih.gov/24692140/>
88. Kenney WL, Craighead DH, Alexander LM. **Heat waves, aging, and human
    cardiovascular health.** Med Sci Sports Exerc. 2014;46(10):1891-9.
    <https://pubmed.ncbi.nlm.nih.gov/24598696/>
89. Semenza JC, Rubin CH, Falter KH, et al. **Heat-related deaths during the
    July 1995 heat wave in Chicago.** N Engl J Med. 1996;335(2):84-90.
    <https://pubmed.ncbi.nlm.nih.gov/8649494/>
90. Westaway K, Frank O, Husband A, et al. **Medicines can affect
    thermoregulation and accentuate the risk of dehydration and heat-related
    illness during hot weather.** J Clin Pharm Ther. 2015;40(4):363-7.
    <https://pubmed.ncbi.nlm.nih.gov/25917210/>
91. Bongers CCWG, Alsady M, Nijenhuis T, et al. **Impact of acute versus
    prolonged exercise and dehydration on kidney function and injury.**
    Physiol Rep. 2018;6(11):e13734. <https://pubmed.ncbi.nlm.nih.gov/29890037/>

> **모델 반영:** 열순화는 SWMAX ×1.85, 발한 역치 −0.35 °C, 혈장량 +11.5%,
> HSP 기저치 상승으로 한 개의 손잡이(`ACCLIM`)에 묶여 있다(84, 85).
> 탈수는 발한을 체중 결손 1%당 5.5% 줄이고(86, 87) 혈장량·MAP·GFR로
> 이어진다. 항콜린제 부담(90)은 `FSW_DRUG`로, 노화는 `FSW_AGE`·`FVD_AGE`로
> 들어간다 — 시카고 1995 폭염(89)에서 확인된 위험인자 구조 그대로다.
> 모델에서 **항콜린제 하나가 고정점의 존부를 뒤집는다**: 같은 40 °C/55% RH
> 방에서 항콜린제가 없는 노인은 경계 위에서 서서히 기어오르고(72시간에
> 38.1 °C), 있는 노인은 고정점을 잃고 달아난다.

---

## 12. 장기손상, 횡문근융해, 예후 (Organ injury & outcome)

92. Argaud L, Ferry T, Le QH, et al. (ref 6 참조 — 장기 예후)
93. Pease S, Bouadma L, Kermarrec N, et al. **Early organ dysfunction course,
    cooling time and outcome in classic heatstroke.** Intensive Care Med.
    2009;35(8):1454-8. <https://pubmed.ncbi.nlm.nih.gov/19404609/>
94. Varghese GM, John G, Thomas K, et al. **Predictors of multi-organ
    dysfunction in heatstroke.** Emerg Med J. 2005;22(3):185-7.
    <https://pubmed.ncbi.nlm.nih.gov/15735266/>
95. Hifumi T, Kondo Y, Shimazaki J, et al. **Prognostic significance of
    disseminated intravascular coagulation in patients with heat stroke.**
    J Intensive Care. 2018;6:33. <https://pubmed.ncbi.nlm.nih.gov/29942522/>
96. Garcin JM, Bronstein JA, Cremades S, et al. **Acute liver failure is
    frequent during heat stroke.** World J Gastroenterol. 2008;14(1):158-9.
    <https://pubmed.ncbi.nlm.nih.gov/18176982/>
97. Zeller L, Novack V, Barski L, et al. **Exertional heatstroke: clinical
    characteristics, diagnostic and therapeutic considerations.** Eur J Intern
    Med. 2011;22(3):296-9. <https://pubmed.ncbi.nlm.nih.gov/21570652/>
98. Bosch X, Poch E, Grau JM. **Rhabdomyolysis and acute kidney injury.**
    N Engl J Med. 2009;361(1):62-72.
    <https://pubmed.ncbi.nlm.nih.gov/19571284/>
99. Lawton EM, Pearce H, Gabb GM. **Review article: environmental heatstroke and
    long-term clinical neurological outcomes.** Emerg Med Australas.
    2019;31(2):163-173. <https://pubmed.ncbi.nlm.nih.gov/30632286/>
100. Yaqub BA. **Neurologic manifestations of heatstroke at the Mecca
     pilgrimage.** Neurology. 1987;37(6):1004-6.
     <https://pubmed.ncbi.nlm.nih.gov/3587615/>

> **모델 반영:** 간 손상은 CEM43·허혈·NAPQI 세 항의 합이고 AST/ALT가
> 48–72시간에 정점(96)에 이르도록 반감기를 맞췄다. 횡문근융해는 **코어가
> 아니라 근육 노드의 CEM43**로 구동된다 — 운동 중 근육은 직장 탐침이 읽는
> 값보다 1–2 °C 높으므로, 같은 코어 온도에서 EHS가 NEHS보다 CK가 훨씬 높은
> 이유가 자동으로 나온다. 소뇌 Purkinje 세포가 가장 열민감한 뉴런이라는
> 관찰(99, 100)은 CNS 손상항의 낮은 역치로 반영했다.

---

## 13. 진료 지침과 합의문 (Guidelines)

101. Casa DJ, DeMartini JK, Bergeron MF, et al. **National Athletic Trainers'
     Association Position Statement: Exertional Heat Illnesses.** J Athl Train.
     2015;50(9):986-1000. <https://pubmed.ncbi.nlm.nih.gov/26381473/>
102. Roberts WO, Armstrong LE, Sawka MN, et al. **ACSM Expert Consensus
     Statement on Exertional Heat Illness: Recognition, Management, and Return
     to Activity.** Curr Sports Med Rep. 2021;20(9):470-484.
     <https://pubmed.ncbi.nlm.nih.gov/34524191/>
103. Lipman GS, Gaudio FG, Eifling KP, et al. **Wilderness Medical Society
     Clinical Practice Guidelines for the Prevention and Treatment of
     Heat Illness: 2019 Update.** Wilderness Environ Med. 2019;30(4S):S33-S46.
     <https://pubmed.ncbi.nlm.nih.gov/31221601/>
104. Epstein Y, Roberts WO. **The pathophysiology of heat stroke: an integrative
     view of the final common pathway.** Scand J Med Sci Sports.
     2011;21(6):742-8. <https://pubmed.ncbi.nlm.nih.gov/21631615/>
105. Hifumi T, Kondo Y, Shimazaki J, et al. **Japanese Association for Acute
     Medicine Heatstroke and Hypothermia Surveillance Committee: heatstroke
     STUDY.** Acute Med Surg. 2019;6(3):270-278.
     <https://pubmed.ncbi.nlm.nih.gov/31304029/>
106. Périard JD, DeGroot D, Jay O. **Exertional heat stroke in sport and the
     military: epidemiology and mitigation.** Exp Physiol.
     2022;107(10):1111-1121. <https://pubmed.ncbi.nlm.nih.gov/35953078/>

> **모델 반영:** 냉각 중단 역치 38.6 °C(101), "cool first, transport
> second"(101-103), 근무-휴식 주기와 WBGT 기반 관리(106)를 시나리오와
> 파라미터 기본값에 반영했다.

---

## 14. 이 모델의 노출된 예측과 한계 (Exposed predictions & limitations)

모델이 문헌과 **어긋나거나** 검증 불가능한 지점을 명시한다. 이것들이
먼저 틀릴 곳이다.

| # | 예측/한계 | 상태 |
|---|-----------|------|
| 1 | 2 °C 대 14 °C 물의 냉각속도 비 **1.55** (메타분석 1.16–1.29) | **과대예측.** 전도도 선형·떨림 상쇄 미모델링·불완전 침수 미고려. 시나리오는 방식별 개별 적합 UA를 사용하므로 보고 속도 자체는 일치 |
| 2 | 열순화가 보상경계를 **전혀 움직이지 못한다** | 반직관적. 모델 내부에서는 `E_actual = E_max,환경`이므로 필연이지만, WMAX 0.85나 대류 가정이 관대하면 실제로는 건조열에서 발한이 제한요인이 될 수 있다 |
| 3 | 확정 용량 **7.8 CEM43** | 역학(30분 내 냉각 시 생존, 60분 지연 시 아님)에 맞춰 보정한 값이지 독립 측정치가 아니다. HMGB1 절대값 보정은 소규모 코호트(61)에 의존 |
| 4 | rTM이 ON 상태를 **소멸**시킨다 | 안장-마디 계산의 직접 귀결. 임상 근거는 후향적 등록연구(70)뿐이고 sepsis DIC의 SCARLET 시험은 null이었다 |
| 5 | 안정시 코어 평형이 **36.8 °C** (정상 37.0) | 능동적 한랭 방어를 넣지 않았기 때문. NEHS 발현 시점을 약간 늦추는 방향으로 편향 |
| 6 | 개체간 변이 없음 | 단일 결정론적 환자. 유전적 소인(RYR1 변이), 최근 발열, 수면부족, 비만 등 위험인자는 파라미터로만 들어간다 |

---

## 15. 방법론 (QSP methodology)

107. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust
     application of systems pharmacology.** CPT Pharmacometrics Syst Pharmacol.
     2016;5(5):235-49. <https://pubmed.ncbi.nlm.nih.gov/27299936/>
108. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and
     systems pharmacology models in R with mrgsolve.** J Pharmacokinet
     Pharmacodyn. 2015;42:S84-S85.
109. Musante CJ, Ramanujan S, Schmidt BJ, et al. **Quantitative systems
     pharmacology: a case for disease models.** Clin Pharmacol Ther.
     2017;101(1):24-27. <https://pubmed.ncbi.nlm.nih.gov/27709613/>
110. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for quantitative
     systems pharmacology (QSP) models.** CPT Pharmacometrics Syst Pharmacol.
     2017;6(8):496-498. <https://pubmed.ncbi.nlm.nih.gov/28643452/>

---

**총 110개 문헌.** 모든 링크는 PubMed 또는 PMC.
파라미터 출처는 `hs_mrgsolve_model.R`의 `$PARAM` 주석과
`hs_core.py`의 `P0` 딕셔너리에 인라인으로도 표시되어 있다.
