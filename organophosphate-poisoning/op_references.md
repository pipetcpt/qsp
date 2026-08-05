# 급성 유기인계 살충제 중독 — 참고문헌
# Acute Organophosphorus (OP) Insecticide Self-Poisoning — References

이 문서는 `op_qsp_model.dot` (기계론적 지도), `op_mrgsolve_model.R` (51-state ODE
모델), `op_reference_model.py` (검증용 Python 재구현), `op_shiny_app.R`
(대시보드)에 들어간 **모든 파라미터와 구조적 가정의 근거**를 섹션별로 정리한
것입니다. 각 항목 끝의 대괄호는 그 문헌이 모델의 어느 숫자로 들어갔는지를
표시합니다.

> **읽는 순서 제안** — 모델의 핵심 산술(옥심 천장 Ω)은 §2, 임상 결과와의 대조는
> §7·§8, 그리고 "왜 무작위 시험이 음성이었는가"는 §8·§13에 있습니다.

---

## 1. 역학과 질병 부담 (Epidemiology and burden)

유기인계 살충제 자의 음독은 전 세계 자살 사망의 상당 부분을 차지하며, 그
치명률은 약물학이 아니라 **어떤 제형이 시장에 있는가**로 결정된다는 것이 이
분야의 가장 확립된 사실입니다. 이 사실은 모델의 `REG`(규제) 노드와 §13의
가상 시험 결과로 반영되어 있습니다.

1. Mew EJ, Padmanathan P, Konradsen F, et al. The global burden of fatal self-poisoning with pesticides 2006-15: systematic review. *J Affect Disord* 2017;219:93-104. [전 세계 사망 부담 · 모델 서문] <https://pubmed.ncbi.nlm.nih.gov/28535450/>
2. Gunnell D, Eddleston M, Phillips MR, Konradsen F. The global distribution of fatal pesticide self-poisoning: systematic review. *BMC Public Health* 2007;7:357. [연간 사망 추정] <https://pubmed.ncbi.nlm.nih.gov/18154668/>
3. Eddleston M, Buckley NA, Eyer P, Dawson AH. Management of acute organophosphorus pesticide poisoning. *Lancet* 2008;371:597-607. [전반적 임상 프레임, 아트로핀·옥심·환기 우선순위] <https://pubmed.ncbi.nlm.nih.gov/17706760/>
4. Eddleston M, Eyer P, Worek F, et al. Differences between organophosphorus insecticides in human self-poisoning: a prospective cohort study. *Lancet* 2005;366:1452-9. [dimethyl vs diethyl 치명률 차이 · `subclass` 파라미터의 근거] <https://pubmed.ncbi.nlm.nih.gov/16243090/>
5. Eddleston M, Eyer P, Worek F, et al. Predicting outcome using butyrylcholinesterase activity in organophosphorus pesticide self-poisoning. *QJM* 2008;101:467-74. [BChE 예후 가치 · `B_free` 출력] <https://pubmed.ncbi.nlm.nih.gov/18375475/>
6. Dawson AH, Eddleston M, Senarathna L, et al. Acute human lethal toxicity of agricultural pesticides: a prospective cohort study. *PLoS Med* 2010;7:e1000357. [제품별 치명률 순위] <https://pubmed.ncbi.nlm.nih.gov/20968586/>
7. Manuweera G, Eddleston M, Egodage S, Buckley NA. Do targeted bans of insecticides to prevent deaths from self-poisoning result in reduced agricultural output? *Environ Health Perspect* 2008;116:492-5. [규제 개입의 실행 가능성] <https://pubmed.ncbi.nlm.nih.gov/18414632/>
8. Knipe DW, Chang SS, Dawson A, et al. Suicide prevention through means restriction: impact of the 2008-2011 pesticide restrictions on suicide in Sri Lanka. *PLoS One* 2017;12:e0172893. [금지 조치의 국가 단위 효과] <https://pubmed.ncbi.nlm.nih.gov/28264041/>

---

## 2. 에스터라제 반응속도론 — 모델의 핵심 (Esterase kinetics: the core)

모델의 3-상태 스위치(E → EP → EP-aged)에 들어가는 **k_i, k_a, k_s, k_r2 값은 모두
아래 Worek 그룹의 인간 AChE 측정값**에서 왔습니다. 특히 옥심 재활성화가
`k_r2 = k_r_max/(K_D + X)` 형태로 **포화한다**는 사실이 모델의 천장
`Ω/(1+Ω)`을 만들어 냅니다.

9. Worek F, Thiermann H, Szinicz L, Eyer P. Kinetic analysis of interactions between human acetylcholinesterase, structurally different organophosphorus compounds and oximes. *Biochem Pharmacol* 2004;68:2237-48. [**k_i, k_a, k_s, k_r_max, K_D의 1차 출처**] <https://pubmed.ncbi.nlm.nih.gov/15498514/>
10. Worek F, Eyer P, Aurbek N, Szinicz L, Thiermann H. Recent advances in evaluation of oxime efficacy in nerve agent poisoning by in vitro analysis. *Toxicol Appl Pharmacol* 2007;219:226-34. [포화형 재활성화 식의 형식화] <https://pubmed.ncbi.nlm.nih.gov/17169391/>
11. Aurbek N, Thiermann H, Szinicz L, Eyer P, Worek F. Analysis of inhibition, reactivation and aging kinetics of highly toxic organophosphorus compounds with human and pig acetylcholinesterase. *Toxicology* 2006;224:91-9. [노화 반감기 diethyl 33 h / dimethyl 3.7 h] <https://pubmed.ncbi.nlm.nih.gov/16720069/>
12. Worek F, Diepold C, Eyer P. Dimethylphosphoryl-inhibited human cholinesterases: inhibition, reactivation, and aging kinetics. *Arch Toxicol* 1999;73:7-14. [dimethyl OP의 빠른 자발적 재활성화 t½ 0.7 h — 모델이 φ를 diethyl보다 *낮게* 계산하는 이유] <https://pubmed.ncbi.nlm.nih.gov/10207609/>
13. Worek F, Thiermann H, Wille T. Oximes in organophosphate poisoning: 60 years of hope and despair. *Chem Biol Interact* 2016;259:93-8. [옥심 실패의 기전적 정리] <https://pubmed.ncbi.nlm.nih.gov/27125983/>
14. Eyer P. The role of oximes in the management of organophosphorus pesticide poisoning. *Toxicol Rev* 2003;22:165-90. [옥심 약동/약력학 종설 · 재억제 현상] <https://pubmed.ncbi.nlm.nih.gov/15181665/>
15. Mason HJ. The recovery of plasma cholinesterase and erythrocyte acetylcholinesterase activity in workers after over-exposure to dichlorvos. *Occup Med (Lond)* 2000;50:343-7. [적혈구 AChE는 재합성되지 않는다 — `KRBCNEW` 값] <https://pubmed.ncbi.nlm.nih.gov/10975133/>
16. Thiermann H, Szinicz L, Eyer P, et al. Correlation between red blood cell acetylcholinesterase activity and neuromuscular transmission in organophosphate poisoning. *Chem Biol Interact* 2005;157-158:345-7. [RBC AChE와 신경근 전달의 상관 — 그리고 그 한계] <https://pubmed.ncbi.nlm.nih.gov/16289002/>
17. Worek F, Mast U, Kiderlen D, Diepold C, Eyer P. Improved determination of acetylcholinesterase activity in human whole blood. *Clin Chim Acta* 1999;288:73-90. [AChE 측정법 — 바이오마커 출력의 정의] <https://pubmed.ncbi.nlm.nih.gov/10529460/>

---

## 3. 독성동태와 생물학적 활성화 (Toxicokinetics and bioactivation)

모델에서 **AChE를 억제하는 것은 모(母)화합물이 아니라 CYP가 만든 옥손**이며,
그 생성이 포화한다는 점(`VMAXBIO`, `KMBIO`)이 "대량 음독에서 옥손 농도가
용량에 무관해진다"는 결론을 만듭니다.

18. Eyer F, Meischner V, Kiderlen D, et al. Human parathion poisoning: a toxicokinetic analysis. *Toxicol Rev* 2003;22:143-63. [모화합물·옥손 혈중 농도의 실측, 지방 저장고 재분포] <https://pubmed.ncbi.nlm.nih.gov/15181664/>
19. Eyer F, Roberts DM, Buckley NA, et al. Extreme variability in the formation of chlorpyrifos oxon (CPO) in patients poisoned by chlorpyrifos (CPF). *Biochem Pharmacol* 2009;78:531-7. [**옥손 생성의 환자 간 극단적 변이 — 모델에서 `ETA_BIO`, `PON1SC`의 근거이자 Ω의 가장 불확실한 입력**] <https://pubmed.ncbi.nlm.nih.gov/19433067/>
20. Timchalk C, Nolan RJ, Mendrala AL, et al. A physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) model for the organophosphate insecticide chlorpyrifos in rats and humans. *Toxicol Sci* 2002;66:34-53. [클로르피리포스 PBPK 구조 · 지방 구획] <https://pubmed.ncbi.nlm.nih.gov/11861971/>
21. Buratti FM, Volpe MT, Meneguz A, Vittozzi L, Testai E. CYP-specific bioactivation of four organophosphorothioate pesticides by human liver microsomes. *Toxicol Appl Pharmacol* 2003;186:143-54. [CYP2B6/3A4/1A2 desulfuration vs dearylation 분기] <https://pubmed.ncbi.nlm.nih.gov/12620367/>
22. Costa LG, Cole TB, Vitalone A, Furlong CE. Measurement of paraoxonase (PON1) status as a potential biomarker of susceptibility to organophosphate toxicity. *Clin Chim Acta* 2005;352:37-47. [PON1 Q192R — `PON1SC = 0.33` 의 근거] <https://pubmed.ncbi.nlm.nih.gov/15653098/>
23. Furlong CE, Marsillach J, Jarvik GP, Costa LG. Paraoxonases-1, -2 and -3: what are their functions? *Chem Biol Interact* 2016;259:51-62. [PON1 기질 특이성: dimethyl 옥손은 거의 가수분해되지 않는다] <https://pubmed.ncbi.nlm.nih.gov/27238723/>
24. Eddleston M, Worek F, Eyer P, et al. Poisoning with the S-alkyl organophosphorus insecticides profenofos and prothiofos. *QJM* 2009;102:785-92. [비전형 OP의 동태 — 모델 라이브러리 확장 근거] <https://pubmed.ncbi.nlm.nih.gov/19737791/>
25. Eddleston M, Street JM, Self I, et al. A role for solvents in the toxicity of agricultural organophosphorus pesticides. *Toxicology* 2012;294:94-103. [**용매(사이클로헥사논·탄화수소)의 독립적 독성 — 모델의 `SOLVCV`, `SOLVRESP`**] <https://pubmed.ncbi.nlm.nih.gov/22365945/>
26. Roberts DM, Aaron CK. Management of acute organophosphorus pesticide poisoning. *BMJ* 2007;334:629-34. [흡수 창, 활성탄의 시간 의존성] <https://pubmed.ncbi.nlm.nih.gov/17379909/>

---

## 4. 아트로핀 — 폐쇄루프 약물 (Atropine as a closed-loop drug)

모델에서 아트로핀은 **투여 스케줄이 아니라 제어기**로 구현되어 있고, 그 결과
누적 용량이 모델의 *출력*이 됩니다. 아래 문헌이 그 목표치(마른 흉부, HR > 80,
수축기혈압 > 80)와 배가 프로토콜의 근거입니다.

27. Eddleston M, Buckley NA, Checketts H, et al. Speed of initial atropinisation in significant organophosphorus pesticide poisoning — a systematic comparison of recommended regimens. *J Toxicol Clin Toxicol* 2004;42:865-75. [배가(doubling) 프로토콜 vs 관행적 적정 — 모델의 `ATRMODE`] <https://pubmed.ncbi.nlm.nih.gov/15533024/>
28. Eddleston M, Dawson A, Karalliedde L, et al. Early management after self-poisoning with an organophosphorus or carbamate pesticide — a treatment protocol for junior doctors. *Crit Care* 2004;8:R391-7. [아트로핀화 종말점의 조작적 정의] <https://pubmed.ncbi.nlm.nih.gov/15566582/>
29. Abedin MJ, Sayeed AA, Basher A, et al. Open-label randomized clinical trial of atropine bolus injection versus incremental boluses plus infusion for organophosphate poisoning in Bangladesh. *J Med Toxicol* 2012;8:108-17. [주입 vs 볼루스 — `KP_RAPID`/`KP_SLOW`] <https://pubmed.ncbi.nlm.nih.gov/22331327/>
30. Perera PMS, Shahmy S, Gawarammana I, Dawson AH. Comparison of two commonly practiced atropinization regimens in acute organophosphorus and carbamate poisoning. *Hum Exp Toxicol* 2008;27:513-8. [누적 아트로핀 용량의 실측 범위] <https://pubmed.ncbi.nlm.nih.gov/18784205/>
31. Ali-Melkkilä T, Kanto J, Iisalo E. Pharmacokinetics and related pharmacodynamics of anticholinergic drugs. *Acta Anaesthesiol Scand* 1993;37:633-42. [**아트로핀 V_d, CL, 그리고 글리코피롤레이트의 4급 암모늄 구조 → BBB 비투과: 모델의 `BBBF` 0.35 vs 0.02**] <https://pubmed.ncbi.nlm.nih.gov/8249551/>
32. Bardin PG, van Eeden SF. Organophosphate poisoning: grading the severity and comparing treatment between atropine and glycopyrrolate. *Crit Care Med* 1990;18:956-60. [글리코피롤레이트 비교 — 시나리오 S15] <https://pubmed.ncbi.nlm.nih.gov/2394120/>
33. Connors NJ, Harnett ZH, Hoffman RS. Comparison of current recommended regimens of atropinization in organophosphate poisoning. *J Med Toxicol* 2014;10:143-7. [권고안 간 비교] <https://pubmed.ncbi.nlm.nih.gov/24619543/>

---

## 5. 옥심 — 약동학과 임상시험 (Oximes: PK and the trials)

34. Eddleston M, Eyer P, Worek F, et al. Pralidoxime in acute organophosphorus insecticide poisoning — a randomised controlled trial. *PLoS Med* 2009;6:e1000104. [**WHO 용법(30 mg/kg + 8 mg/kg/h)의 무작위 시험, 음성 결과 — 모델 §13의 가상 시험이 재현하려는 대상**] <https://pubmed.ncbi.nlm.nih.gov/19564902/>
35. Buckley NA, Eddleston M, Li Y, Bevan M, Robertson J. Oximes for acute organophosphate pesticide poisoning. *Cochrane Database Syst Rev* 2011;(2):CD005085. [체계적 고찰: 이득 없음] <https://pubmed.ncbi.nlm.nih.gov/21328273/>
36. Pawar KS, Bhoite RR, Pillay CP, Chavan SC, Malshikare DS, Garad SG. Continuous pralidoxime infusion versus repeated bolus injection to treat organophosphorus pesticide poisoning: a randomised controlled trial. *Lancet* 2006;368:2136-41. [지속 주입의 우월성 — 모델의 S02 vs S03 대조] <https://pubmed.ncbi.nlm.nih.gov/17174705/>
37. Medicis JJ, Stork CM, Howland MA, Hoffman RS, Goldfrank LR. Pharmacokinetics following a loading plus a continuous infusion of pralidoxime compared with the traditional short infusion regimen in human volunteers. *J Toxicol Clin Toxicol* 1996;34:289-95. [프랄리독심 PK 파라미터 `OXV1/OXCL`] <https://pubmed.ncbi.nlm.nih.gov/8667465/>
38. Thiermann H, Szinicz L, Eyer F, et al. Modern strategies in therapy of organophosphate poisoning. *Toxicol Lett* 1999;107:233-9. [오비독심 용법 250 mg + 750 mg/24 h] <https://pubmed.ncbi.nlm.nih.gov/10414801/>
39. Thiermann H, Zilker T, Eyer F, Felgenhauer N, Eyer P, Worek F. Monitoring of neuromuscular transmission in organophosphate pesticide-poisoned patients. *Toxicol Lett* 2009;191:297-304. [실제 환자에서 재활성화가 관찰되는 조건] <https://pubmed.ncbi.nlm.nih.gov/19825404/>
40. Thiermann H, Worek F, Kehe K. Limitations and challenges in treatment of acute chemical warfare agent poisoning. *Chem Biol Interact* 2013;206:435-43. [옥심 천장 개념의 정성적 선행 서술] <https://pubmed.ncbi.nlm.nih.gov/23994499/>
41. Kiderlen D, Eyer P, Worek F. Formation and disposition of diethylphosphoryl-obidoxime, a potent anticholinesterase that is hydrolyzed by human paraoxonase (PON1). *Biochem Pharmacol* 2005;69:1853-67. [인산화 옥심의 재억제 — 모델의 `REBOUND` 노드] <https://pubmed.ncbi.nlm.nih.gov/15869746/>
42. Rahimi R, Nikfar S, Abdollahi M. Increased morbidity and mortality in acute human organophosphate-poisoned patients treated by oximes: a meta-analysis of clinical trials. *Hum Exp Toxicol* 2006;25:157-62. [**옥심 사용군에서 사망률이 더 높았다는 메타분석 — 모델이 (의도치 않게) 재현하는 방향과 일치하며, README의 '가장 노출된 예측'으로 표시**] <https://pubmed.ncbi.nlm.nih.gov/16696287/>
43. Peter JV, Moran JL, Graham P. Oxime therapy and outcomes in human organophosphate poisoning: an evaluation using meta-analytic techniques. *Crit Care Med* 2006;34:502-10. [동일 주제의 독립 메타분석] <https://pubmed.ncbi.nlm.nih.gov/16424734/>
44. Syed S, Gurcoo SA, Farooqui AK, Nisa W, Sofi K, Wani TM. Is the World Health Organization-recommended dose of pralidoxime effective in the treatment of organophosphorus poisoning? *Indian J Anaesth* 2015;59:31-6. [WHO 용법의 재검증 시도] <https://pubmed.ncbi.nlm.nih.gov/25684811/>

---

## 6. 니코틴성 증후군과 중간증후군 (Nicotinic limb and intermediate syndrome)

45. Senanayake N, Karalliedde L. Neurotoxic effects of organophosphorus insecticides: an intermediate syndrome. *N Engl J Med* 1987;316:761-3. [중간증후군의 최초 기술 — 모델의 `Rn_des`/`IMS`] <https://pubmed.ncbi.nlm.nih.gov/3029588/>
46. Jayawardane P, Dawson AH, Weerasinghe V, Karalliedde L, Buckley NA, Senanayake N. The spectrum of intermediate syndrome following acute organophosphate poisoning: a prospective cohort study from Sri Lanka. *PLoS Med* 2008;5:e147. [발현 시점 24-96 h, 지속 4-18일 — `KDES`, `KRECDES`] <https://pubmed.ncbi.nlm.nih.gov/18630983/>
47. Jayawardane P, Senanayake N, Buckley NA, Dawson AH. Electrophysiological correlates of the intermediate syndrome in acute organophosphate poisoning. *Clin Toxicol* 2012;50:250-3. [반복자극 감쇠 — 탈분극 차단의 전기생리] <https://pubmed.ncbi.nlm.nih.gov/22385107/>
48. Karalliedde L, Baker D, Marrs TC. Organophosphate-induced intermediate syndrome: aetiology and relationships with myopathy. *Toxicol Rev* 2006;25:1-14. [기전 가설의 정리] <https://pubmed.ncbi.nlm.nih.gov/16856766/>
49. Wood SJ, Slater CR. Safety factor at the neuromuscular junction. *Prog Neurobiol* 2001;64:393-429. [**신경근 접합부 안전계수 — 모델의 `SF_NMJ` = 0.28, 즉 안정 시 환기에는 최대 흡기력의 28%만 필요하다는 구조**] <https://pubmed.ncbi.nlm.nih.gov/11275359/>
50. Yang CC, Deng JF. Intermediate syndrome following organophosphate insecticide poisoning. *J Chin Med Assoc* 2007;70:467-72. [임상 경과와 인공환기 소요] <https://pubmed.ncbi.nlm.nih.gov/18063499/>

---

## 7. 호흡부전과 지지요법 (Respiratory failure and supportive care)

50개가 넘는 임상 관찰이 한 방향을 가리킵니다: **이 병에서 생존을 결정하는
단일 요소는 인공환기의 가용성**이며, 모델의 시나리오 S13(환기기 없음)이 이를
정량화합니다.

51. Eddleston M, Mohamed F, Davies JOJ, et al. Respiratory failure in acute organophosphorus pesticide self-poisoning. *QJM* 2006;99:513-22. [**호흡부전의 원인 분해: 중추 억제 · 분비물 · 신경근 차단 — 모델 §8 구조의 직접 근거**] <https://pubmed.ncbi.nlm.nih.gov/16829539/>
52. Hulse EJ, Davies JOJ, Simpson AJ, Sciuto AM, Eddleston M. Respiratory complications of organophosphorus nerve agent and insecticide poisoning: implications for respiratory and critical care. *Am J Respir Crit Care Med* 2014;190:1342-54. [기도 분비물, 흡인, ARDS] <https://pubmed.ncbi.nlm.nih.gov/25419614/>
53. Gaspari RJ, Paydarfar D. Respiratory recovery following organophosphate poisoning in a rat model is suppressed by isolated hypoxia at the point of apnea. *Toxicology* 2011;290:59-64. [중추성 무호흡의 저산소 악순환] <https://pubmed.ncbi.nlm.nih.gov/21875643/>
54. Carey JL, Dunn C, Gaspari RJ. Central respiratory failure during acute organophosphate poisoning. *Respir Physiol Neurobiol* 2013;189:403-10. [preBötzinger 리듬 발생기의 실패 — 모델의 `PREBOT`/`DMAXNONR`] <https://pubmed.ncbi.nlm.nih.gov/23933009/>
55. Eddleston M, Chowdhury FR. Pharmacological treatment of organophosphorus insecticide poisoning: the old and the (possible) new. *Br J Clin Pharmacol* 2016;81:462-70. [현행 및 후보 치료의 비판적 정리] <https://pubmed.ncbi.nlm.nih.gov/26366467/>
56. Munidasa UADD, Gawarammana IB, Kularatne SAM, Kumarasiri PVR, Goonasekera CDA. Survival pattern in patients with acute organophosphate poisoning receiving intensive care. *J Toxicol Clin Toxicol* 2004;42:343-7. [중환자실 생존 곡선 — 모델 사망률 보정 기준] <https://pubmed.ncbi.nlm.nih.gov/15461240/>
57. Senanayake N, de Silva HJ, Karalliedde L. A scale to assess severity in organophosphorus intoxication: POP scale. *Hum Exp Toxicol* 1993;12:297-9. [**Peradeniya OP Poisoning scale — 모델 `$TABLE`에서 상태벡터로부터 재구성하는 점수**] <https://pubmed.ncbi.nlm.nih.gov/8104007/>
58. Davies JOJ, Eddleston M, Buckley NA. Predicting outcome in acute organophosphorus poisoning with a poison severity score or the Glasgow coma scale. *QJM* 2008;101:371-9. [중증도 점수의 예측력] <https://pubmed.ncbi.nlm.nih.gov/18319295/>

---

## 8. 중추신경계 · 경련 · 벤조디아제핀 (CNS, seizures, benzodiazepines)

59. McDonough JH Jr, Shih TM. Neuropharmacological mechanisms of nerve agent-induced seizure and neuropathology. *Neurosci Biobehav Rev* 1997;21:559-79. [콜린성 → 글루타민성 전환, 30분 창 — 모델의 `SEIZ`→`NMDA`] <https://pubmed.ncbi.nlm.nih.gov/9353792/>
60. Shih TM, Duniho SM, McDonough JH. Control of nerve agent-induced seizures is critical for neuroprotection and survival. *Toxicol Appl Pharmacol* 2003;188:69-80. [디아제팜의 시간 의존적 효과] <https://pubmed.ncbi.nlm.nih.gov/12691725/>
61. Eddleston M, Buckley NA, Checketts H, et al. Diazepam in organophosphorus poisoning. (In: Management reviews above; see also) Marrs TC, Rice P, Vale JA. The role of oximes in the treatment of nerve agent poisoning in civilian casualties. *Toxicol Rev* 2006;25:297-323. [디아제팜 병용의 근거] <https://pubmed.ncbi.nlm.nih.gov/17288500/>
62. Bird SB, Gaspari RJ, Dickson EW. Early death due to severe organophosphate poisoning is a centrally mediated process. *Acad Emerg Med* 2003;10:295-8. [**조기 사망은 중추성이다 — 모델이 아트로핀의 BBB 투과를 결과에 연결하는 이유**] <https://pubmed.ncbi.nlm.nih.gov/12670843/>
63. Chen Y. Organophosphate-induced brain damage: mechanisms, neuropsychiatric and neurological consequences, and potential therapeutic strategies. *Neurotoxicology* 2012;33:391-400. [COPIND] <https://pubmed.ncbi.nlm.nih.gov/22498093/>

---

## 9. 지연성 다발신경병증 (OPIDN) 과 NTE

64. Lotti M, Moretto A. Organophosphate-induced delayed polyneuropathy. *Toxicol Rev* 2005;24:37-49. [**NTE 억제 + 노화 > 70%가 발병 조건 — 모델의 `NTETH`**] <https://pubmed.ncbi.nlm.nih.gov/16042503/>
65. Glynn P. Neuropathy target esterase and phospholipid deacylation. *Biochim Biophys Acta* 2005;1736:87-93. [NTE(PNPLA6)의 생화학] <https://pubmed.ncbi.nlm.nih.gov/16137924/>
66. Jokanović M, Kosanović M, Brkić D, Vukomanović P. Organophosphate induced delayed polyneuropathy in man: an overview. *Clin Neurol Neurosurg* 2011;113:7-10. [원인 화합물 목록 — dimethyl OP는 거의 유발하지 않는다] <https://pubmed.ncbi.nlm.nih.gov/20880635/>

---

## 10. 보조 치료 · 실패한 후보들 (Adjuncts and the failed candidates)

67. Pajoumand A, Shadnia S, Rezaie A, Abdi M, Abdollahi M. Benefits of magnesium sulfate in the management of acute human poisoning by organophosphorus insecticides. *Hum Exp Toxicol* 2004;23:565-9. [황산마그네슘 — 시나리오 S16] <https://pubmed.ncbi.nlm.nih.gov/15688984/>
68. Basher A, Rahman SH, Ghose A, Arif SM, Faiz MA, Dawson AH. Phase II study of magnesium sulfate in acute organophosphate pesticide poisoning. *Clin Toxicol* 2013;51:35-40. [용량 탐색] <https://pubmed.ncbi.nlm.nih.gov/23148565/>
69. Perera PM, Jayamanna SF, Hettiarachchi R, et al. A phase II clinical trial to assess the safety of clonidine in acute organophosphorus pesticide poisoning. *Trials* 2009;10:73. [클로니딘 — 모델의 `CLONID` 노드] <https://pubmed.ncbi.nlm.nih.gov/19678923/>
70. Eddleston M, Juszczak E, Buckley NA, et al. Multiple-dose activated charcoal in acute self-poisoning: a randomised controlled trial. *Lancet* 2008;371:579-87. [**다회 활성탄: 이득 없음 — 모델은 그 이유를 흡수 창의 산술로 제시**] <https://pubmed.ncbi.nlm.nih.gov/18280328/>
71. Li Y, Tse ML, Gawarammana I, Buckley N, Eddleston M. Systematic review of controlled clinical trials of gastric lavage in acute organophosphorus pesticide poisoning. *Clin Toxicol* 2009;47:179-92. [위세척] <https://pubmed.ncbi.nlm.nih.gov/19306191/>
72. Peter JV, Moran JL, Graham PL. Advances in the management of organophosphate poisoning. *Expert Opin Pharmacother* 2007;8:1451-64. [보조요법 종설] <https://pubmed.ncbi.nlm.nih.gov/17661728/>
73. Chowdhary S, Bhattacharyya R, Banerjee D. Acute organophosphorus poisoning. *Clin Chim Acta* 2014;431:66-76. [진단검사 관점] <https://pubmed.ncbi.nlm.nih.gov/24508989/>

---

## 11. 생촉매/생포집제 (Bioscavengers) — 왜 살충제에는 통하지 않는가

모델 §6은 **화학량론적 포집제가 살충제 음독에서는 산술적으로 불가능**함을
보입니다(음독량 mmol vs 투여 가능한 포집제 µmol). 신경작용제 예방에서는 같은
전략이 성립하는데, 그 차이가 정확히 몰수에 있다는 점이 아래 문헌들과 일치합니다.

74. Lenz DE, Yeung D, Smith JR, Sweeney RE, Lumley LA, Cerasoli DM. Stoichiometric and catalytic scavengers as protection against nerve agent toxicity: a mini review. *Toxicology* 2007;233:31-9. [화학량론적 포집의 원리와 몰수 제약] <https://pubmed.ncbi.nlm.nih.gov/17188793/>
75. Nachon F, Brazzolotto X, Trovaslet M, Masson P. Progress in the development of enzyme-based nerve agent bioscavengers. *Chem Biol Interact* 2013;206:536-44. [촉매형 포집제] <https://pubmed.ncbi.nlm.nih.gov/23811386/>
76. Rosenberg YJ, Laube B, Mao L, et al. Pulmonary delivery of an aerosolized recombinant human butyrylcholinesterase pretreatment protects against aerosolized paraoxon in macaques. *Chem Biol Interact* 2013;203:167-71. [rHuBChE 전처치] <https://pubmed.ncbi.nlm.nih.gov/23103586/>
77. Ashani Y, Pistinner S. Estimation of the upper limit of human butyrylcholinesterase dose required for protection against organophosphates toxicity: a mathematically based toxicokinetic model. *Toxicol Sci* 2004;77:358-67. [**필요한 BChE 용량의 수학적 상한 — 모델 §6의 선행 계산**] <https://pubmed.ncbi.nlm.nih.gov/14691206/>

---

## 12. 심혈관 독성 (Cardiovascular toxicity)

78. Ludomirsky A, Klein HO, Sarelli P, et al. Q-T prolongation and polymorphous ("torsade de pointes") ventricular arrhythmias associated with organophosphorus insecticide poisoning. *Am J Cardiol* 1982;49:1654-8. [QTc·염전성 심실빈맥] <https://pubmed.ncbi.nlm.nih.gov/7081053/>
79. Yurumez Y, Yavuz Y, Saglam H, et al. Electrocardiographic findings of acute organophosphate poisoning. *J Emerg Med* 2009;36:39-42. [심전도 소견의 빈도] <https://pubmed.ncbi.nlm.nih.gov/18024061/>
80. Chuang FR, Jang SW, Lin JL, Chern MS, Chen JB, Hsu KT. QTc prolongation indicates a poor prognosis in patients with organophosphate poisoning. *Am J Emerg Med* 1996;14:451-3. [예후 표지자] <https://pubmed.ncbi.nlm.nih.gov/8765108/>
81. Anand S, Singh S, Nahar Saikia U, Bhalla A, Paul Sharma Y, Singh D. Cardiac abnormalities in acute organophosphate poisoning. *Clin Toxicol* 2009;47:230-5. [심근 손상] <https://pubmed.ncbi.nlm.nih.gov/19306191/>

---

## 13. 모델링 방법론 (QSP / PBPK methodology)

82. Baron KT, Gastonguay MR. Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve. *J Pharmacokinet Pharmacodyn* 2015;42:S84-5. [mrgsolve] <https://mrgsolve.org>
83. Timchalk C, Poet TS. Development of a physiologically based pharmacokinetic and pharmacodynamic model to determine dosimetry and cholinesterase inhibition for a binary mixture of chlorpyrifos and diazinon in the rat. *Neurotoxicology* 2008;29:428-43. [OP PBPK/PD 선행 모델] <https://pubmed.ncbi.nlm.nih.gov/18359092/>
84. Gearhart JM, Jepson GW, Clewell HJ 3rd, Andersen ME, Conolly RB. Physiologically based pharmacokinetic and pharmacodynamic model for the inhibition of acetylcholinesterase by diisopropylfluorophosphate. *Toxicol Appl Pharmacol* 1990;106:295-310. [AChE 억제의 PBPK/PD 원형] <https://pubmed.ncbi.nlm.nih.gov/2260096/>
85. Bhattacharya S, Conolly RB, Kaminski NE, Thomas RS, Andersen ME, Zhang Q. A bistable switch underlying B-cell differentiation and its disruption by the environmental contaminant TCDD. *Toxicol Sci* 2010;115:51-65. [독성 시스템의 스위치 해석 방법론 — 본 모델의 3-상태 스위치 해석에 사용한 사고 틀] <https://pubmed.ncbi.nlm.nih.gov/20123757/>
86. Sterner TR, Ruark CD, Covington TR, Yu KO, Gearhart JM. A physiologically based pharmacokinetic model for the oxime TMB-4: simulation of rodent and human data. *Arch Toxicol* 2013;87:661-80. [옥심 PBPK] <https://pubmed.ncbi.nlm.nih.gov/23151815/>

---

## 14. 지침과 종설 (Guidelines and general reviews)

87. World Health Organization. *The WHO Recommended Classification of Pesticides by Hazard and Guidelines to Classification 2019.* Geneva: WHO, 2020. [WHO Class I/II 분류 — 모델의 `REG` 노드] <https://www.who.int/publications/i/item/9789240005662>
88. Eddleston M. Novel clinical toxicology and pharmacology of organophosphorus insecticide self-poisoning. *Annu Rev Pharmacol Toxicol* 2019;59:341-60. [**이 모델이 재현하려는 임상 서사 전체의 가장 좋은 단일 요약**] <https://pubmed.ncbi.nlm.nih.gov/30110584/>
89. Vale A, Lotti M. Organophosphorus and carbamate insecticide poisoning. *Handb Clin Neurol* 2015;131:149-68. [신경학 관점의 종설] <https://pubmed.ncbi.nlm.nih.gov/26563788/>
90. Robb EL, Baker MB. Organophosphate toxicity. *StatPearls* [Internet]. Treasure Island (FL): StatPearls Publishing, 2023. [교육용 요약] <https://www.ncbi.nlm.nih.gov/books/NBK470430/>

---

## 부록 A — 모델 파라미터별 출처 대응표 (Parameter provenance)

| 파라미터 | 값 (기본: 클로르피리포스) | 출처 |
|---|---|---|
| `KI` 옥손의 AChE 억제 상수 | 0.18 nM⁻¹h⁻¹ (= 3.0×10⁶ M⁻¹min⁻¹) | 9, 11 |
| `T12AGE` 노화 반감기 | diethyl 33 h / dimethyl 3.7 h | 9, 11, 12 |
| `T12SPO` 자발적 재활성화 반감기 | diethyl 77 h / dimethyl 0.7 h | 9, 12 |
| `KRMAX`, `KDOX` 옥심 재활성화 포화 | 2-PAM 36 h⁻¹ / 200 µM; 오비독심 48 h⁻¹ / 40 µM | 9, 10, 13 |
| `VMAXBIO`, `KMBIO` CYP 탈황화 | 6.0×10⁴ nmol/h, 20 µM | 19, 20, 21 |
| `CLOXON`, `PON1SC` 옥손 가수분해 | 300 L/h, QQ에서 0.33배 | 22, 23 |
| `VFTH` 지방 저장고 | 700 L (fenthion 1750 L) | 18, 20 |
| `OXV1/OXCL` 프랄리독심 PK | 24 L / 21 L/h | 34, 37 |
| `ATR_V/ATR_CL/BBBF` 아트로핀 | 200 L / 46 L/h / 0.35 | 31 |
| `BBBF` 글리코피롤레이트 | 0.02 | 31, 32 |
| `SF_NMJ` 신경근 안전계수 | 0.28 | 49 |
| `KDES/KRECDES` 중간증후군 | 0.075 / 0.050 h⁻¹ | 45, 46, 47 |
| `NTETH` OPIDN 임계 | 노화 NTE 70% | 64, 66 |
| `HVENT` 인공환기 합병증 위험 | 0.0012 h⁻¹ (≈2.9%/일) | 51, 52, 56 |
| `SOLVCV` 용매 심근 억제 | 사이클로헥사논 1.8 vs 탄화수소 0.4-0.5 | 25 |
| 가상 시험의 음독량 분포 | 로그정규, 중앙값 28 mL | 4, 6, 34 |

## 부록 B — 모델이 재현에 성공/실패한 임상 관찰

| 임상 관찰 | 모델 | 근거 문헌 |
|---|---|---|
| dimethyl OP의 치명률이 diethyl보다 높다 | 재현 (노화 속도 + 옥손 농도) | 4 |
| 옥심은 AChE를 재활성화하지만 사망률을 바꾸지 않는다 | 재현 (§13 가상 시험 RR ≈ 1) | 34, 35, 43 |
| 지속 주입이 볼루스보다 낫다 | 재현 (S02 vs S03) | 36 |
| 빠른 아트로핀 배가가 관행적 적정보다 낫다 | 재현 (S02 vs S14) | 27, 29 |
| 다회 활성탄은 이득이 없다 | 재현 (S17 vs S18: 1시간 이후 효과 소실) | 70 |
| 조기 사망은 중추성이다 | 재현 (S15 글리코피롤레이트에서 악화) | 62 |
| BChE는 예후 표지자다 | 부분 재현 (출력으로 제공, 예후식은 미적합) | 5 |
| 옥심군의 사망률이 오히려 높다는 메타분석 | **의도치 않게 재현** — README의 '가장 노출된 예측' 참조 | 42, 43 |
