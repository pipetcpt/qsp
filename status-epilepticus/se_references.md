# 경련지속상태 QSP 모델 — 참고문헌
# Status Epilepticus QSP Model — References

이 모델의 모든 구조적 주장과 파라미터는 아래 문헌에 대응합니다.
**모든 PMID는 NCBI E-utilities로 조회하여 저자·저널·연도·제목을 확인했습니다.**

각 항목 끝의 `→` 표시는 그 문헌이 모델의 **어느 부분**에 쓰였는지를 가리킵니다.
(`§3` = 기계론적 지도의 클러스터 번호, `PARAM` = mrgsolve 파라미터, `CAL` = 보정 근거)

---

## 1. 정의 · 역학 · 진료지침 (Definition · epidemiology · guidelines)

1. Trinka E, et al. **A definition and classification of status epilepticus — Report of the ILAE Task Force on Classification of Status Epilepticus.** Epilepsia 2015. — [PMID 26336950](https://pubmed.ncbi.nlm.nih.gov/26336950/)
   → t1 = 5분(자가지속 전환), t2 = 30분(손상 시작). 모델의 `N_transit` 노드와 `INJURY` 적분기의 시간 척도가 여기에서 옵니다. §2, §9

2. Glauser T, et al. **Evidence-Based Guideline: Treatment of Convulsive Status Epilepticus in Children and Adults (American Epilepsy Society).** Epilepsy Curr 2016. — [PMID 26900382](https://pubmed.ncbi.nlm.nih.gov/26900382/)
   → 1차/2차/3차 투여 시각과 용량의 기준값. 시나리오 02·21

3. Brophy GM, et al. **Guidelines for the evaluation and management of status epilepticus.** Neurocrit Care 2012. — [PMID 22528274](https://pubmed.ncbi.nlm.nih.gov/22528274/)
   → 3차 마취제 용량 범위(midazolam 0.05-2 mg/kg/h, propofol 2-10 mg/kg/h). §14, PARAM

4. DeLorenzo RJ, et al. **A prospective, population-based epidemiologic study of status epilepticus in Richmond, Virginia.** Neurology 1996. — [PMID 8780085](https://pubmed.ncbi.nlm.nih.gov/8780085/)
   → 발생률 및 병인별 분포. §1

5. Betjemann JP, Lowenstein DH. **Status epilepticus in adults.** Lancet Neurol 2015. — [PMID 25908090](https://pubmed.ncbi.nlm.nih.gov/25908090/)
   → 전반적 임상 프레임. 전 클러스터

6. Neligan A, Shorvon SD. **Frequency and prognosis of convulsive status epilepticus of different causes: a systematic review.** Arch Neurol 2010. — [PMID 20697043](https://pubmed.ncbi.nlm.nih.gov/20697043/)
   → 병인별 사망률 — `EDRIVE`가 결과 분산을 지배한다는 §18 해리 3의 근거

7. Neligan A, et al. **Change in mortality of generalized convulsive status epilepticus in high-income countries over time.** JAMA Neurol 2019. — [PMID 31135807](https://pubmed.ncbi.nlm.nih.gov/31135807/)
   → 치료 발전에도 사망률이 크게 변하지 않았다는 관찰 — 모델의 "선택보다 시계가 지배한다"는 주장과 부합

8. Kantanen AM, et al. **Incidence and mortality of super-refractory status epilepticus in adults.** Epilepsy Behav 2015. — [PMID 26141934](https://pubmed.ncbi.nlm.nih.gov/26141934/)
   → SRSE 빈도 및 사망률. §18 `K_rse`, `K_mort`

9. Kantanen AM, et al. **Long-term outcome of refractory status epilepticus in adults: a retrospective population-based study.** Epilepsy Res 2017. — [PMID 28402834](https://pubmed.ncbi.nlm.nih.gov/28402834/)

10. Kantanen AM, et al. **Predictors of hospital and one-year mortality in intensive care patients with refractory status epilepticus.** Crit Care 2017. — [PMID 28330483](https://pubmed.ncbi.nlm.nih.gov/28330483/)

11. Lowenstein DH, Bleck T, Macdonald RL. **It's time to revise the definition of status epilepticus.** Epilepsia 1999. — [PMID 9924914](https://pubmed.ncbi.nlm.nih.gov/9924914/)
    → 30분 → 5분 정의 전환의 논거가 바로 이 모델의 시계 논리입니다

12. Fountain NB. **Status epilepticus: risk factors and complications.** Epilepsia 2000. — [PMID 10885737](https://pubmed.ncbi.nlm.nih.gov/10885737/) → §10, §11

13. Mayer SA, et al. **Refractory status epilepticus: frequency, risk factors, and impact on outcome.** Arch Neurol 2002. — [PMID 11843690](https://pubmed.ncbi.nlm.nih.gov/11843690/)

14. Novy J, Logroscino G, Rossetti AO. **Refractory status epilepticus: a prospective observational study.** Epilepsia 2010. — [PMID 19817823](https://pubmed.ncbi.nlm.nih.gov/19817823/)

---

## 2. 자가지속(self-sustainment)과 시간 의존성

15. Lothman E. **The biochemical basis and pathophysiology of status epilepticus.** Neurology 1990. — [PMID 2185436](https://pubmed.ncbi.nlm.nih.gov/2185436/)
    → 보상기 → 비보상기 전환(약 30분). §10 → §11

16. Mazarati AM, et al. **Self-sustaining status epilepticus after brief electrical stimulation of the perforant path.** Brain Res 1998. — [PMID 9729413](https://pubmed.ncbi.nlm.nih.gov/9729413/)
    → 유발 자극을 제거해도 발작이 지속된다 = `SEIZ` 상태변수의 쌍안정성. §2

17. Wasterlain CG, et al. **Self-sustaining status epilepticus: a condition maintained by potentiation of glutamate receptors and by plasticity of the GABA-A receptor.** Epilepsia 2000. — [PMID 10999535](https://pubmed.ncbi.nlm.nih.gov/10999535/)
    → **이 모델의 두 시계 개념 자체의 출처.** §3, §5

18. Wasterlain CG, Chen JW. **Mechanistic and pharmacologic aspects of status epilepticus and its treatment with new antiepileptic drugs.** Epilepsia 2008. — [PMID 19087119](https://pubmed.ncbi.nlm.nih.gov/19087119/)

19. Chen JW, Wasterlain CG. **Status epilepticus: pathophysiology and management in adults.** Lancet Neurol 2006. — [PMID 16488380](https://pubmed.ncbi.nlm.nih.gov/16488380/)

20. Lowenstein DH, Alldredge BK. **Status epilepticus.** N Engl J Med 1998. — [PMID 9521986](https://pubmed.ncbi.nlm.nih.gov/9521986/)

21. Shorvon S, Ferlisi M. **The treatment of super-refractory status epilepticus: a critical review of available therapies and a clinical treatment protocol.** Brain 2011. — [PMID 21914716](https://pubmed.ncbi.nlm.nih.gov/21914716/) → §14

22. Ferlisi M, Shorvon S. **The outcome of therapies in refractory and super-refractory convulsive status epilepticus and recommendations for therapy.** Brain 2012. — [PMID 22577217](https://pubmed.ncbi.nlm.nih.gov/22577217/)

23. Bleck TP. **Refractory status epilepticus.** Curr Opin Crit Care 2005. — [PMID 15758590](https://pubmed.ncbi.nlm.nih.gov/15758590/)

---

## 3. **시계 1** — 시냅스 γ2-GABA-A 수용체의 내재화 (모델의 중심축)

24. **Kapur J, Macdonald RL. Rapid seizure-induced reduction of benzodiazepine and Zn²⁺ sensitivity of hippocampal dentate granule cell GABA-A receptors.** J Neurosci 1997. — [PMID 9295398](https://pubmed.ncbi.nlm.nih.gov/9295398/)
    → **diazepam 역가 약 20배 저하 vs phenobarbital 약 3배.** 모델은 이 숫자를 입력하지 않고, 세 개의 풀(§3·§4)에서 `CREQ_BZD` / `CREQ_PB`로 재현합니다 (7.0배 vs 2.6배 @30분 — 방향과 비대칭은 재현, 크기는 과소추정). CAL

25. **Naylor DE, Liu H, Wasterlain CG. Trafficking of GABA(A) receptors, loss of inhibition, and a mechanism for pharmacoresistance in status epilepticus.** J Neurosci 2005. — [PMID 16120773](https://pubmed.ncbi.nlm.nih.gov/16120773/)
    → **SE 1시간 후 표면 γ2 약 47% 감소.** `KENDO`/`KREC`/`KDEGR`은 이 한 점에 맞췄고, 30분 값(0.69)은 결과입니다. CAL, PARAM

26. Goodkin HP, Yeh JL, Kapur J. **Status epilepticus increases the intracellular accumulation of GABA-A receptors.** J Neurosci 2005. — [PMID 15944379](https://pubmed.ncbi.nlm.nih.gov/15944379/)
    → `RENDO` 구획(엔도솜 풀)의 존재 근거

27. **Goodkin HP, et al. Subunit-specific trafficking of GABA(A) receptors during status epilepticus.** J Neurosci 2008. — [PMID 18322097](https://pubmed.ncbi.nlm.nih.gov/18322097/)
    → **β2/3·γ2 함유 수용체는 내재화되지만 δ 함유는 그렇지 않다.** 모델에서 `RSYN`은 내려가고 `REXTRA`는 내려가지 않는 이유. §3 vs §4

28. Terunuma M, et al. **Deficits in phosphorylation of GABA(A) receptors by intimately associated protein kinase C activity underlie compromised synaptic inhibition during status epilepticus.** J Neurosci 2008. — [PMID 18184780](https://pubmed.ncbi.nlm.nih.gov/18184780/)
    → `T_ca` → `T_gephyrin` 인산화 경로

29. Kittler JT, et al. **Phospho-dependent binding of the clathrin AP2 adaptor complex to GABA-A receptors regulates the efficacy of inhibitory synaptic transmission.** PNAS 2005. — [PMID 16192353](https://pubmed.ncbi.nlm.nih.gov/16192353/)
    → `T_ap2` 노드

30. Jacob TC, Moss SJ, Jurd R. **GABA(A) receptor trafficking and its role in the dynamic modulation of neuronal inhibition.** Nat Rev Neurosci 2008. — [PMID 18382465](https://pubmed.ncbi.nlm.nih.gov/18382465/)

31. Feng HJ, Mathews GC, Kao C, Macdonald RL. **Alterations of GABA-A receptor function and allosteric modulation during development of status epilepticus.** J Neurophysiol 2008. — [PMID 18216225](https://pubmed.ncbi.nlm.nih.gov/18216225/)
    → `FBZS`(벤조디아제핀 감수성 분획)라는 **두 번째** 감소 항의 근거

32. Deeb TZ, Maguire J, Moss SJ. **Possible alterations in GABA-A receptor signaling that underlie benzodiazepine-resistant seizures.** Epilepsia 2012. — [PMID 23216581](https://pubmed.ncbi.nlm.nih.gov/23216581/)

33. Goodkin HP, Kapur J. **The impact of diazepam's discovery on the treatment and understanding of status epilepticus.** Epilepsia 2009. — [PMID 19674049](https://pubmed.ncbi.nlm.nih.gov/19674049/)

34. Joshi S, Kapur J. **Mechanisms of status epilepticus: α-amino-3-hydroxy-5-methyl-4-isoxazolepropionic acid receptor hypothesis.** Epilepsia 2018. — [PMID 30159880](https://pubmed.ncbi.nlm.nih.gov/30159880/)
    → `AMPACP`(Ca 투과성 AMPA 증가) 상태변수

---

## 4. **남아 있는 풀** — 시냅스외 δ-GABA-A와 신경스테로이드

35. Farrant M, Nusser Z. **Variations on an inhibitory theme: phasic and tonic activation of GABA(A) receptors.** Nat Rev Neurosci 2005. — [PMID 15738957](https://pubmed.ncbi.nlm.nih.gov/15738957/)
    → phasic(시냅스) vs tonic(시냅스외) 구분 — 모델의 `WSYN`/`WEXT` 분해

36. Rogawski MA, Loya CM, Reddy K, et al. **Neuroactive steroids for the treatment of status epilepticus.** Epilepsia 2013. — [PMID 24001085](https://pubmed.ncbi.nlm.nih.gov/24001085/)
    → δ 수용체는 벤조디아제핀 비감수성 · 신경스테로이드 고효능. `FSYNALLO = 0.40`

37. Reddy DS, Rogawski MA. **Neurosteroid replacement therapy for catamenial epilepsy.** Neurotherapeutics 2009. — [PMID 19332335](https://pubmed.ncbi.nlm.nih.gov/19332335/) → `D_allowith` 노드

38. Zolkowska D, et al. **Intramuscular allopregnanolone and ganaxolone in a mouse model of treatment-resistant status epilepticus.** Epilepsia 2018. — [PMID 29453777](https://pubmed.ncbi.nlm.nih.gov/29453777/)

39. Rosenthal ES, et al. **Brexanolone as adjunctive therapy in super-refractory status epilepticus.** Ann Neurol 2017. — [PMID 28779545](https://pubmed.ncbi.nlm.nih.gov/28779545/)
    → 시나리오 18의 용량(86 µg/kg/h)과 시점

40. Vaitkevicius H, et al. **Intravenous ganaxolone for the treatment of refractory status epilepticus: results from an open-label, dose-finding, phase 2 trial.** Epilepsia 2022. — [PMID 35748707](https://pubmed.ncbi.nlm.nih.gov/35748707/)

41. Singh RK, et al. **Intravenous ganaxolone in pediatric super-refractory status epilepticus.** Epilepsy Behav Rep 2022. — [PMID 36325100](https://pubmed.ncbi.nlm.nih.gov/36325100/)

42. Gasior M, et al. **Intravenous ganaxolone: pharmacokinetics, pharmacodynamics, safety, and tolerability in healthy adults.** Clin Pharmacol Drug Dev 2024. — [PMID 38231434](https://pubmed.ncbi.nlm.nih.gov/38231434/)

---

## 5. **시계 2** — 시냅스 NMDA 수용체의 유입, 그리고 케타민

43. **Naylor DE, et al. Rapid surface accumulation of NMDA receptors increases glutamatergic excitation during status epilepticus.** Neurobiol Dis 2013. — [PMID 23313318](https://pubmed.ncbi.nlm.nih.gov/23313318/)
    → **SE 1시간 후 시냅스 NR1 약 38% 증가.** `KEXO`/`KENDN`/`FEXO0`는 이 한 점에 맞췄습니다. CAL

44. Mazarati AM, Wasterlain CG. **N-methyl-D-aspartate receptor antagonists abolish the maintenance phase of self-sustaining status epilepticus in rat.** Neurosci Lett 1999. — [PMID 10327162](https://pubmed.ncbi.nlm.nih.gov/10327162/)
    → 유지기에는 NMDA 차단이 듣는다 = 케타민의 시간 부호가 양(+)인 이유

45. Rice AC, DeLorenzo RJ. **NMDA receptor activation during status epilepticus is required for the development of epilepsy.** Brain Res 1998. — [PMID 9519269](https://pubmed.ncbi.nlm.nih.gov/9519269/) → §19 `Z_epilep`

46. Fujikawa DG. **Neuroprotective effect of ketamine administered after status epilepticus onset.** Epilepsia 1995. — [PMID 7821277](https://pubmed.ncbi.nlm.nih.gov/7821277/)

47. Gaspard N, et al. **Intravenous ketamine for the treatment of refractory status epilepticus: a retrospective multicenter study.** Epilepsia 2013. — [PMID 23758557](https://pubmed.ncbi.nlm.nih.gov/23758557/)
    → 케타민 용량(1.5-3 mg/kg 볼루스, 1-10 mg/kg/h). PARAM

48. Rosati A, et al. **Efficacy and safety of ketamine in refractory status epilepticus in children.** Neurology 2012. — [PMID 23197747](https://pubmed.ncbi.nlm.nih.gov/23197747/)

49. Rosati A, De Masi S, Guerrini R. **Ketamine for refractory status epilepticus: a systematic review.** CNS Drugs 2018. — [PMID 30232735](https://pubmed.ncbi.nlm.nih.gov/30232735/)

50. Alkhachroum A, et al. **Ketamine to treat super-refractory status epilepticus.** Neurology 2020. — [PMID 32873691](https://pubmed.ncbi.nlm.nih.gov/32873691/)

51. Ilvento L, et al. **Ketamine in refractory convulsive status epilepticus in children avoids endotracheal intubation.** Epilepsy Behav 2015. — [PMID 26189786](https://pubmed.ncbi.nlm.nih.gov/26189786/)
    → 모델에서 케타민만 MAP을 **올리고**(`KETMAP`) 호흡 부담을 늘리지 않는 이유

52. **Niquet J, et al. Simultaneous triple therapy for the treatment of status epilepticus.** Neurobiol Dis 2017. — [PMID 28461248](https://pubmed.ncbi.nlm.nih.gov/28461248/)
    → 두 개의 움직이는 풀을 **동시에** 겨냥하는 다제요법 — 시나리오 15·17의 근거

---

## 6. 염화물 항상성 — GABA가 흥분성이 되는 경로

53. Rivera C, et al. **BDNF-induced TrkB activation down-regulates the K⁺-Cl⁻ cotransporter KCC2 and impairs neuronal Cl⁻ extrusion.** J Cell Biol 2002. — [PMID 12473684](https://pubmed.ncbi.nlm.nih.gov/12473684/) → `C_bdnf` → `KCC2`

54. Pathak HR, et al. **Disrupted dentate granule cell chloride regulation enhances synaptic excitability during development of temporal lobe epilepsy.** J Neurosci 2007. — [PMID 18094240](https://pubmed.ncbi.nlm.nih.gov/18094240/)

55. Kaila K, et al. **Cation-chloride cotransporters in neuronal development, plasticity and disease.** Nat Rev Neurosci 2014. — [PMID 25234263](https://pubmed.ncbi.nlm.nih.gov/25234263/)

56. **Burman RJ, et al. Excitatory GABAergic signalling is associated with benzodiazepine resistance in status epilepticus.** Brain 2019. — [PMID 31553050](https://pubmed.ncbi.nlm.nih.gov/31553050/)
    → 네 번째 실패 모드. 모델에서 **GABA 작용제 자체가 Cl⁻을 싣는다**(`FCLDRG`)는 되먹임의 근거

57. Pressler RM, et al. **Bumetanide for the treatment of seizures in newborn babies with hypoxic ischaemic encephalopathy (NEMO).** Lancet Neurol 2015. — [PMID 25765333](https://pubmed.ncbi.nlm.nih.gov/25765333/)
    → NKCC1 표적화가 임상에서 실패한 근거 — `C_nkcc1` 노드에 명시

---

## 7. 내인성 종결자와 촉진자

58. Young D, Dragunow M. **Status epilepticus may be caused by loss of adenosine anticonvulsant mechanisms.** Neuroscience 1994. — [PMID 8152537](https://pubmed.ncbi.nlm.nih.gov/8152537/) → `ADO` 상태변수

59. Aronica E, et al. **Glial adenosine kinase — a neuropathological marker of the epileptic brain.** Neurochem Int 2013. — [PMID 23385089](https://pubmed.ncbi.nlm.nih.gov/23385089/) → `ADK` 상승

60. Mazarati AM, et al. **Galanin modulation of seizures and seizure modulation of hippocampal galanin in animal models of status epilepticus.** J Neurosci 1998. — [PMID 9822761](https://pubmed.ncbi.nlm.nih.gov/9822761/) → `NETPEP` 고갈

61. Mazarati AM, et al. **Modulation of hippocampal excitability and seizures by galanin.** J Neurosci 2000. — [PMID 10934278](https://pubmed.ncbi.nlm.nih.gov/10934278/)

62. Liu H, et al. **Substance P is expressed in hippocampal principal neurons during status epilepticus and plays a critical role in the maintenance of status epilepticus.** PNAS 1999. — [PMID 10220458](https://pubmed.ncbi.nlm.nih.gov/10220458/)
    → 억제성 펩타이드는 고갈되고 흥분성 펩타이드는 상승 = `NETPEP`의 단조 감소

---

## 8. 신경염증 · 혈액뇌장벽

63. Vezzani A, French J, Bartfai T, Baram TZ. **The role of inflammation in epilepsy.** Nat Rev Neurol 2011. — [PMID 21135885](https://pubmed.ncbi.nlm.nih.gov/21135885/)

64. Maroso M, et al. **Toll-like receptor 4 and high-mobility group box-1 are involved in ictogenesis and can be targeted to reduce seizures.** Nat Med 2010. — [PMID 20348922](https://pubmed.ncbi.nlm.nih.gov/20348922/) → `I_hmgb1`

65. **Balosso S, et al. A novel non-transcriptional pathway mediates the proconvulsive effects of interleukin-1β.** Brain 2008. — [PMID 18952671](https://pubmed.ncbi.nlm.nih.gov/18952671/)
    → IL-1β → Src → NR2B Tyr1472 인산화 → NMDA 전류 증가. 모델의 `GIL` 항

66. Friedman A, Kaufer D, Heinemann U. **Blood-brain barrier breakdown-inducing astrocytic transformation: novel targets for the prevention of epilepsy.** Epilepsy Res 2009. — [PMID 19362806](https://pubmed.ncbi.nlm.nih.gov/19362806/) → 알부민-TGF-β → `Z_epilep`

67. Marchi N, et al. **Seizure-promoting effect of blood-brain barrier disruption.** Epilepsia 2007. — [PMID 17319915](https://pubmed.ncbi.nlm.nih.gov/17319915/)

68. Kenney-Jung DL, et al. **Febrile infection-related epilepsy syndrome treated with anakinra.** Ann Neurol 2016. — [PMID 27770579](https://pubmed.ncbi.nlm.nih.gov/27770579/) → 시나리오 20

69. Lai YC, et al. **Anakinra usage in febrile infection related epilepsy syndrome: an international cohort.** Ann Clin Transl Neurol 2020. — [PMID 33506622](https://pubmed.ncbi.nlm.nih.gov/33506622/)

70. Gaspard N, et al. **New-onset refractory status epilepticus: etiology, clinical features, and outcome.** Neurology 2015. — [PMID 26296517](https://pubmed.ncbi.nlm.nih.gov/26296517/)

71. Sculier C, Gaspard N. **New onset refractory status epilepticus (NORSE).** Seizure 2019. — [PMID 30482654](https://pubmed.ncbi.nlm.nih.gov/30482654/)

---

## 9. 흥분독성 신경 손상 — 두 번째 적분기

72. **Meldrum BS, Brierley JB. Prolonged epileptic seizures in primates: ischemic cell change and its relation to ictal physiological events.** Arch Neurol 1973. — [PMID 4629379](https://pubmed.ncbi.nlm.nih.gov/4629379/)
    → 인공호흡·혈압 유지 하에서도 손상이 생긴다 = `INJURY`가 `SEIZ`와 **별개의** 상태변수여야 하는 이유

73. Meldrum BS, Horton RW. **Physiology of status epilepticus in primates.** Arch Neurol 1973. — [PMID 4629380](https://pubmed.ncbi.nlm.nih.gov/4629380/) → §10 → §11 전환

74. Wasterlain CG, Fujikawa DG, Penix L, Sankar R. **Pathophysiological mechanisms of brain damage from status epilepticus.** Epilepsia 1993. — [PMID 8385002](https://pubmed.ncbi.nlm.nih.gov/8385002/)

75. Fujikawa DG. **The temporal evolution of neuronal damage from pilocarpine-induced status epilepticus.** Brain Res 1996. — [PMID 8828581](https://pubmed.ncbi.nlm.nih.gov/8828581/) → `NEURLOSS` 포화 곡선

76. DeGiorgio CM, et al. **Serum neuron-specific enolase in the major subtypes of status epilepticus.** Neurology 1999. — [PMID 10078721](https://pubmed.ncbi.nlm.nih.gov/10078721/) → `X_nse` 노드

---

## 10. 예후 · 결과 지표

77. Rossetti AO, Logroscino G, Bromfield EB. **Status Epilepticus Severity Score (STESS).** J Neurol 2008. — [PMID 18769858](https://pubmed.ncbi.nlm.nih.gov/18769858/)

78. Leitinger M, et al. **Epidemiology-based mortality score in status epilepticus (EMSE).** Neurocrit Care 2015. — [PMID 25412806](https://pubmed.ncbi.nlm.nih.gov/25412806/)

79. Leitinger M, et al. **Predicting outcome of status epilepticus.** Epilepsy Behav 2015. — [PMID 26071999](https://pubmed.ncbi.nlm.nih.gov/26071999/)

80. Giovannini G, et al. **Mortality, morbidity and refractoriness prediction in status epilepticus: comparison of STESS and EMSE scores.** Seizure 2017. — [PMID 28226274](https://pubmed.ncbi.nlm.nih.gov/28226274/)

81. Hocker S, Tatum WO, LaRoche S, Freeman WD. **Refractory and super-refractory status epilepticus — an update.** Curr Neurol Neurosci Rep 2014. — [PMID 24760477](https://pubmed.ncbi.nlm.nih.gov/24760477/)

82. Pujar SS, et al. **Long-term prognosis after childhood convulsive status epilepticus: a prospective cohort study.** Lancet Child Adolesc Health 2018. — [PMID 30169233](https://pubmed.ncbi.nlm.nih.gov/30169233/) → §19 `Z_epilep`

---

## 11. 1차 치료 임상시험 (first line)

83. **Treiman DM, et al. A comparison of four treatments for generalized convulsive status epilepticus (VA Cooperative Study).** N Engl J Med 1998. — [PMID 9738086](https://pubmed.ncbi.nlm.nih.gov/9738086/)
    → **현성 SE에서 lorazepam 64.9%, subtle SE에서 7.7-24.2%.** 같은 약, 다른 창. 모델의 §18 해리 1 및 `MOTG`/`MOTDRG` 파라미터. CAL

84. **Alldredge BK, et al. A comparison of lorazepam, diazepam, and placebo for the treatment of out-of-hospital status epilepticus (PHTSE).** N Engl J Med 2001. — [PMID 11547716](https://pubmed.ncbi.nlm.nih.gov/11547716/)
    → lorazepam 59.1% / diazepam 42.6% / 위약 21.1%; **호흡 합병증은 위약 22.5% vs lorazepam 10.6%** — 발작이 약보다 호흡을 더 억제한다. `KRESPS > KRESPD`의 근거. CAL

85. **Silbergleit R, et al. Intramuscular versus intravenous therapy for prehospital status epilepticus (RAMPART).** N Engl J Med 2012. — [PMID 22335736](https://pubmed.ncbi.nlm.nih.gov/22335736/)
    → **IM midazolam 73.4% vs IV lorazepam 63.4%.** RAMPART 전용 파라미터는 모델에 하나도 없습니다 — 흡수 지연 < 투여 지연이라는 산술로만 재현합니다 (모델 80.0 vs 78.8, 방향 일치·폭 과소). CAL

86. Welch RD, et al. **Intramuscular midazolam versus intravenous lorazepam for the prehospital treatment of status epilepticus in the pediatric population.** Epilepsia 2015. — [PMID 25597369](https://pubmed.ncbi.nlm.nih.gov/25597369/)

87. Leppik IE, et al. **Double-blind study of lorazepam and diazepam in status epilepticus.** JAMA 1983. — [PMID 6131148](https://pubmed.ncbi.nlm.nih.gov/6131148/) → diazepam의 재분포·조기 재발

88. McMullan J, et al. **Midazolam versus diazepam for the treatment of status epilepticus in children and young adults: a meta-analysis.** Acad Emerg Med 2010. — [PMID 20624136](https://pubmed.ncbi.nlm.nih.gov/20624136/)

89. **Sathe AG, et al. Underdosing of benzodiazepines in patients with status epilepticus enrolled in the Established Status Epilepticus Treatment Trial.** Acad Emerg Med 2019. — [PMID 31161706](https://pubmed.ncbi.nlm.nih.gov/31161706/)
    → 저용량 투여가 예외가 아니라 규칙. 시나리오 05·06의 근거

90. Sathe AG, et al. **Patterns of benzodiazepine underdosing in the Established Status Epilepticus Treatment Trial.** Epilepsia 2021. — [PMID 33567109](https://pubmed.ncbi.nlm.nih.gov/33567109/)

---

## 12. 2차·3차 치료 임상시험

91. **Kapur J, et al. Randomized trial of three anticonvulsant medications for status epilepticus (ESETT).** N Engl J Med 2019. — [PMID 31774955](https://pubmed.ncbi.nlm.nih.gov/31774955/)
    → **levetiracetam 47% / fosphenytoin 45% / valproate 46%.** 모델의 `SLEV`·`SPHT`·`SVPA`(및 대응 EC50)는 이 **동률**에만 맞춘 유일한 효능 파라미터입니다 (모델 46.6 / 50.3 / 47.2). CAL

92. Chamberlain JM, et al. **Efficacy of levetiracetam, fosphenytoin, and valproate for established status epilepticus by age group (ESETT).** Lancet 2020. — [PMID 32203691](https://pubmed.ncbi.nlm.nih.gov/32203691/)

93. Navarro V, et al. **Prehospital treatment with levetiracetam plus clonazepam or placebo plus clonazepam in status epilepticus (SAMUKeppra).** Lancet Neurol 2016. — [PMID 26627366](https://pubmed.ncbi.nlm.nih.gov/26627366/)
    → 병원 전 levetiracetam 추가는 음성 — 모델에서 2차 약제의 이득이 1차 실패자에 국한된다는 예측과 부합

94. Lyttle MD, et al. **Levetiracetam versus phenytoin for second-line treatment of paediatric convulsive status epilepticus (EcLiPSE).** Lancet 2019. — [PMID 31005385](https://pubmed.ncbi.nlm.nih.gov/31005385/)

95. Dalziel SR, et al. **Levetiracetam versus phenytoin for second-line treatment of convulsive status epilepticus in children (ConSEPT).** Lancet 2019. — [PMID 31005386](https://pubmed.ncbi.nlm.nih.gov/31005386/)

96. Claassen J, Hirsch LJ, Emerson RG, Mayer SA. **Treatment of refractory status epilepticus with pentobarbital, propofol, or midazolam: a systematic review.** Epilepsia 2002. — [PMID 11903460](https://pubmed.ncbi.nlm.nih.gov/11903460/) → §14

97. Rossetti AO, et al. **A randomized trial for the treatment of refractory status epilepticus.** Neurocrit Care 2011. — [PMID 20878265](https://pubmed.ncbi.nlm.nih.gov/20878265/)

98. Prabhakar H, et al. **Propofol versus thiopental sodium for the treatment of refractory status epilepticus.** Cochrane Database Syst Rev 2017. — [PMID 28155226](https://pubmed.ncbi.nlm.nih.gov/28155226/)

99. Legriel S, et al. **Hypothermia for neuroprotection in convulsive status epilepticus (HYBERNATUS).** N Engl J Med 2016. — [PMID 28002714](https://pubmed.ncbi.nlm.nih.gov/28002714/)
    → EEG 발작은 줄었으나 기능 예후는 개선 없음 — `COOL` 플래그가 `INJURY`에만 작용하고 `EDRIVE`에는 작용하지 않는 이유

---

## 13. 약동학 (PK 파라미터의 출처)

100. Greenblatt DJ, et al. **Pharmacokinetic comparison of sublingual lorazepam with intravenous, intramuscular, and oral lorazepam.** J Pharm Sci 1982. — [PMID 6121043](https://pubmed.ncbi.nlm.nih.gov/6121043/) → `VLZP`, `CLLZP`

101. Ramsay RE, DeToledo J. **Intravenous administration of fosphenytoin: options for the management of seizures.** Neurology 1996. — [PMID 8649609](https://pubmed.ncbi.nlm.nih.gov/8649609/) → `KCONV`, 주입속도 제한

102. Clements JA, Nimmo WS. **Pharmacokinetics and analgesic effect of ketamine in man.** Br J Anaesth 1981. — [PMID 7459184](https://pubmed.ncbi.nlm.nih.gov/7459184/) → `V1KET`/`V2KET`/`CLKET`

103. Patsalos PN, et al. **Serum protein binding of 25 antiepileptic drugs in a routine clinical setting: a comparison of free non-protein-bound concentrations.** Epilepsia 2017. — [PMID 28542801](https://pubmed.ncbi.nlm.nih.gov/28542801/)
     → `FUPHT0 = 0.10`, valproate의 포화성 결합(`FUVPA0`/`FUVPAM`/`FUVPA50`)

104. Patsalos PN, et al. **Therapeutic drug monitoring of antiepileptic drugs in epilepsy: a 2018 update.** Ther Drug Monit 2018. — [PMID 29957667](https://pubmed.ncbi.nlm.nih.gov/29957667/)

---

## 14. 약제 내성 — 수송체 가설 vs 표적 가설

105. Löscher W, Potschka H. **Drug resistance in brain diseases and the role of drug efflux transporters.** Nat Rev Neurosci 2005. — [PMID 16025095](https://pubmed.ncbi.nlm.nih.gov/16025095/)

106. Remy S, Beck H. **Molecular and cellular mechanisms of pharmacoresistance in epilepsy.** Brain 2006. — [PMID 16317026](https://pubmed.ncbi.nlm.nih.gov/16317026/)
     → 두 가설의 구분. 모델의 `BPRATIO` 출력이 이 둘을 갈라놓습니다

107. Bankstahl JP, Löscher W. **Resistance to antiepileptic drugs and expression of P-glycoprotein in two rat models of status epilepticus.** Epilepsy Res 2008. — [PMID 18760905](https://pubmed.ncbi.nlm.nih.gov/18760905/) → `PGP` 상태변수

108. Bauer B, et al. **Seizure-induced up-regulation of P-glycoprotein at the blood-brain barrier through glutamate and cyclooxygenase-2 signaling.** Mol Pharmacol 2008. — [PMID 18094072](https://pubmed.ncbi.nlm.nih.gov/18094072/) → `I_cox2` → `R_pgp` 간선

109. Zibell G, et al. **Prevention of seizure-induced up-regulation of endothelial P-glycoprotein by COX-2 inhibition.** Neuropharmacology 2009. — [PMID 19371577](https://pubmed.ncbi.nlm.nih.gov/19371577/)

110. Potschka H, Fedrowitz M, Löscher W. **P-glycoprotein-mediated efflux of phenobarbital, lamotrigine, and felbamate at the blood-brain barrier.** Neurosci Lett 2002. — [PMID 12113905](https://pubmed.ncbi.nlm.nih.gov/12113905/)
     → 약물별 P-gp 민감도(`FPGPBZD` ≪ `FPGPPHT`)의 근거

---

## 15. 비경련성 SE와 EEG — 해리 1의 실증

111. Treiman DM, Walton NY, Kendrick C. **A progressive sequence of electroencephalographic changes during generalized convulsive status epilepticus.** Epilepsy Res 1990. — [PMID 2303022](https://pubmed.ncbi.nlm.nih.gov/2303022/)
     → 운동 발현이 사라져도 EEG는 진행한다 — `MOTG` 감쇠의 직접 근거

112. DeLorenzo RJ, et al. **Persistent nonconvulsive status epilepticus after the control of convulsive status epilepticus.** Epilepsia 1998. — [PMID 9701373](https://pubmed.ncbi.nlm.nih.gov/9701373/)
     → 임상적으로 종료된 환자의 48%에서 EEG 발작이 지속

113. Hirsch LJ, et al. **American Clinical Neurophysiology Society's Standardized Critical Care EEG Terminology: 2021 version.** J Clin Neurophysiol 2021. — [PMID 33475321](https://pubmed.ncbi.nlm.nih.gov/33475321/)

114. Claassen J, et al. **Detection of electrographic seizures with continuous EEG monitoring in critically ill patients.** Neurology 2004. — [PMID 15159471](https://pubmed.ncbi.nlm.nih.gov/15159471/)

115. Leitinger M, et al. **Salzburg consensus criteria for non-convulsive status epilepticus — approach to clinical application.** Epilepsy Behav 2015. — [PMID 26092326](https://pubmed.ncbi.nlm.nih.gov/26092326/)

116. Trinka E, Leitinger M. **Which EEG patterns in coma are nonconvulsive status epilepticus?** Epilepsy Behav 2015. — [PMID 26148985](https://pubmed.ncbi.nlm.nih.gov/26148985/)

117. Shneker BF, Fountain NB. **Assessment of acute morbidity and mortality in nonconvulsive status epilepticus.** Neurology 2003. — [PMID 14581666](https://pubmed.ncbi.nlm.nih.gov/14581666/)

---

## 16. 지연 사슬 — 모델에서 가장 강력하고, 유일하게 순수 운영적인 변수

118. **Gaínza-Lein M, et al. Association of time to treatment with short-term outcomes for pediatric patients with refractory convulsive status epilepticus.** JAMA Neurol 2018. — [PMID 29356811](https://pubmed.ncbi.nlm.nih.gov/29356811/)
     → 첫 투여 지연이 사망·발작 지속시간과 독립적으로 연관 — `Q_cliff`의 임상 대응

119. Sánchez Fernández I, et al. **Factors associated with treatment delays in pediatric refractory convulsive status epilepticus.** Neurology 2018. — [PMID 29643084](https://pubmed.ncbi.nlm.nih.gov/29643084/)

120. Kämppi L, Ritvanen J, Mustonen H, Soinila S. **Analysis of the delay components in the treatment of status epilepticus.** Neurocrit Care 2013. — [PMID 23817962](https://pubmed.ncbi.nlm.nih.gov/23817962/)
     → §17 지연 사슬의 실측 구성요소

121. Kämppi L, et al. **The essence of the first 2.5 h in the treatment of generalized convulsive status epilepticus.** Seizure 2018. — [PMID 29306214](https://pubmed.ncbi.nlm.nih.gov/29306214/)

122. Guterman EL, et al. **Prehospital midazolam use and outcomes among patients with out-of-hospital status epilepticus.** Neurology 2020. — [PMID 32943481](https://pubmed.ncbi.nlm.nih.gov/32943481/)
     → 병원 전 벤조디아제핀이 실제로는 소수에게만 투여된다는 실세계 데이터

123. Vasquez A, et al. **Hospital emergency treatment of convulsive status epilepticus: comparison of pathways from ten pediatric research centers.** Pediatr Neurol 2018. — [PMID 30075875](https://pubmed.ncbi.nlm.nih.gov/30075875/)

---

## 17. 전신 생리 — 보상기와 비보상기

124. Simon RP. **Physiologic consequences of status epilepticus.** Epilepsia 1985. — [PMID 3922751](https://pubmed.ncbi.nlm.nih.gov/3922751/)
     → `CATAMP`/`TAUCAT`, 자동조절 실패, 저혈당·고체온의 시간 경과

125. Walton NY. **Systemic effects of generalized convulsive status epilepticus.** Epilepsia 1993. — [PMID 8462491](https://pubmed.ncbi.nlm.nih.gov/8462491/)

126. Manno EM, et al. **Cardiac pathology in status epilepticus.** Ann Neurol 2005. — [PMID 16240367](https://pubmed.ncbi.nlm.nih.gov/16240367/) → `S2_takot` 노드

---

## 이 모델이 반증될 수 있는 네 지점 (Falsifiable predictions)

문헌으로 **아직 검증되지 않은** 모델의 출력들입니다. 이 중 하나라도 틀리면 구조가 틀린 것입니다.

| # | 예측 | 모델이 내놓는 값 | 검증 방법 |
|---|------|-----------------|-----------|
| 1 | 벤조디아제핀 필요 농도의 **수직 점근선** | 기본 환자에서 45-50분 사이 | 투여 시각별 1차 반응률의 급락 지점을 전향적으로 측정 |
| 2 | 2차 약제의 절벽은 **존재하되 약 40분 늦다** | 1차 25분 vs 2차 62분 | ESETT 유형 시험에서 2차 투여 시각을 계층화한 반응률 |
| 3 | 벤조디아제핀 회수 가능 풀과 NMDA 제거 가능 지분의 **교차 시각** | 280분 | 케타민 병용 시점을 무작위 배정한 시험 |
| 4 | 20분 종료 vs 60분 종료의 **누적 흥분독성 부하 차이** | 신경세포 손실 9.3% vs 42.7% | 종료 시각별 NSE/NfL 및 해마 용적 추적 |

---

## 면책 (Disclaimer)

본 모델은 **교육·연구 목적의 정성적·반정량적 QSP 모델**입니다.
공개 문헌에 근거하여 구성되었으나 독립적으로 검증·인증되지 않았으며,
**실제 임상 의사결정, 처방, 또는 규제 제출에 사용해서는 안 됩니다.**
