# 급성 심근경색 (STEMI) QSP 모델 — 참고문헌
## Acute Myocardial Infarction · Annotated Reference List

이 목록은 `ami_qsp_model.dot`의 각 클러스터와 `ami_mrgsolve_model.R`의 각
파라미터가 어디에서 왔는지를 추적할 수 있도록 **기전 블록별로** 정리한
것입니다. 각 항목 뒤의 대괄호는 모델에서 그 문헌이 실제로 무엇을 고정하는지
(어떤 구조 또는 어떤 수치를 정하는지) 표시합니다.

> **읽는 방법.** ★ 표시는 이 모델의 골격을 직접 결정한 문헌입니다. 나머지는
> 파라미터의 크기(order of magnitude)를 잡거나, 모델이 재현해야 할 임상
> 관찰을 제공합니다. 모델이 **재현하지 못하거나 문헌과 어긋나는 지점은
> §12에 따로 모아 두었습니다** — 숨기지 않는 것이 이 저장소의 규칙입니다.

---

## 1 · 괴사 파면과 관통벽 기울기 (The Wavefront) — 모델의 1번 경주

1. ★ Reimer KA, Lowe JE, Rasmussen MM, Jennings RB. **The wavefront phenomenon of ischemic cell death. 1. Myocardial infarct size vs duration of coronary occlusion in dogs.** *Circulation* 1977;56:786-794. <https://pubmed.ncbi.nlm.nih.gov/912839/> — [모델의 층 구조 자체. 심내막→심외막 진행과 40 min/3 h/6 h/24 h 시점의 AAR 대비 괴사 분율이 `KE`·`KNI`·`GC` 보정의 1차 앵커]
2. ★ Reimer KA, Jennings RB. **The "wavefront phenomenon" of myocardial ischemic cell death. II. Transmural progression of necrosis within the framework of ischemic bed size (myocardium at risk) and collateral flow.** *Lab Invest* 1979;40:633-644. <https://pubmed.ncbi.nlm.nih.gov/449273/> — [`GC` 관통벽 부수혈류 기울기와 "생존 테두리"의 근거]
3. Jennings RB, Reimer KA. **The cell biology of acute myocardial ischemia.** *Annu Rev Med* 1991;42:225-246. <https://pubmed.ncbi.nlm.nih.gov/2035969/> — [E(에너지 충전도) 상태변수의 개념적 근거]
4. Schaper W, Görge G, Winkler B, Schaper J. **The collateral circulation of the heart.** *Prog Cardiovasc Dis* 1988;31:57-77. <https://pubmed.ncbi.nlm.nih.gov/3293127/> — [`COLL` 기본값 0.10과 그 종간·개인간 변동폭]
5. Rentrop KP, Cohen M, Blanke H, Phillips RA. **Changes in collateral channel filling immediately after controlled coronary artery occlusion by an angioplasty balloon in human subjects.** *J Am Coll Cardiol* 1985;5:587-592. <https://pubmed.ncbi.nlm.nih.gov/3156171/> — [사람에서 부수혈관의 즉시 동원 = `RENTRY` 노드]
6. Meier P, Gloekler S, Zbinden R, et al. **Beneficial effect of recruitable collaterals: a 10-year follow-up study in patients with stable coronary artery disease undergoing quantitative collateral measurements.** *Circulation* 2007;116:975-983. <https://pubmed.ncbi.nlm.nih.gov/17679611/> — [부수혈류가 예후를 바꾼다는 임상 확인 — 실험 3의 논지]
7. Christian TF, Schwartz RS, Gibbons RJ. **Determinants of infarct size in reperfusion therapy for acute myocardial infarction.** *Circulation* 1992;86:81-90. <https://pubmed.ncbi.nlm.nih.gov/1617791/> — [AAR·부수혈류·재관류 시간의 상대적 기여 — 세 변수 모두 모델에 존재]
8. Hausenloy DJ, Yellon DM. **Myocardial ischemia-reperfusion injury: a neglected therapeutic target.** *J Clin Invest* 2013;123:92-100. <https://pubmed.ncbi.nlm.nih.gov/23281415/> — [허혈 손상과 재관류 손상의 분리 — 모델이 NI와 NR을 따로 적분하는 이유]

## 2 · 동면 심근과 기초 대사 요구 (Hibernation & the BASAL Criterion)

9. ★ Rahimtoola SH. **The hibernating myocardium.** *Am Heart J* 1989;117:211-221. <https://pubmed.ncbi.nlm.nih.gov/2783527/> — [모델의 생존 판정 기준: 기초 요구를 지불할 수 있으면 무운동이더라도 생존]
10. Heusch G. **Hibernating myocardium.** *Physiol Rev* 1998;78:1055-1085. <https://pubmed.ncbi.nlm.nih.gov/9790569/> — [혈류–기능 관계(perfusion-contraction matching) = `CF`(감당 가능한 수축) 식의 직접 근거]
11. Ross J Jr. **Myocardial perfusion-contraction matching. Implications for coronary heart disease and hibernation.** *Circulation* 1991;83:1076-1083. <https://pubmed.ncbi.nlm.nih.gov/1999009/> — [수축이 공급에 맞추어 내려간다는 관측 — `CF = afford × E`]
12. Gibbs CL. **Cardiac energetics.** *Physiol Rev* 1978;58:174-254. <https://pubmed.ncbi.nlm.nih.gov/146862/> — [`DEM_BAS ≈ 0.20`: 정지 심장의 기초 산소소모가 박동 심장의 약 20%]
13. Suga H. **Ventricular energetics.** *Physiol Rev* 1990;70:247-277. <https://pubmed.ncbi.nlm.nih.gov/2181496/> — [PVA–MVO₂ 관계, `DEMAND ∝ rate-pressure product` 근사의 정당화]
14. Braunwald E, Kloner RA. **The stunned myocardium: prolonged, postischemic ventricular dysfunction.** *Circulation* 1982;66:1146-1149. <https://pubmed.ncbi.nlm.nih.gov/6754130/> — [`STUN` 상태변수와 `KST_OFF`(회복 시간상수 ~3.5 d)]
15. Bolli R, Marbán E. **Molecular and cellular mechanisms of myocardial stunning.** *Physiol Rev* 1999;79:609-634. <https://pubmed.ncbi.nlm.nih.gov/10221990/> — [기절의 ROS·Ca²⁺ 이중 원인 → `dSTUN/dt`가 두 항의 합인 이유]

## 3 · 혐기 대사, 산증, 그리고 산물 억제 (Anaerobic Metabolism & Acidosis)

16. Neely JR, Grotyohann LW. **Role of glycolytic products in damage to ischemic myocardium.** *Circ Res* 1984;55:816-824. <https://pubmed.ncbi.nlm.nih.gov/6499136/> — [`HI50`: 해당 산물이 스스로 해당을 멈춘다 — 모델에서 저혈류 층의 공급이 붕괴하는 이유]
17. Rovetto MJ, Lamberton WF, Neely JR. **Mechanisms of glycolytic inhibition in ischemic rat hearts.** *Circ Res* 1975;37:742-751. <https://pubmed.ncbi.nlm.nih.gov/128461/> — [PFK/GAPDH 수준의 억제 기전 — `PRODINH` 노드]
18. Cross HR, Opie LH, Radda GK, Clarke K. **Is a high glycogen content beneficial or detrimental to the ischemic rat heart?** *Circ Res* 1996;78:482-491. <https://pubmed.ncbi.nlm.nih.gov/8593707/> — [`KGLY`와 글리코겐의 양면성: 초기 완충, 후기 산증]
19. Dennis SC, Gevers W, Opie LH. **Protons in ischemia: where do they come from; where do they go to?** *J Mol Cell Cardiol* 1991;23:1077-1086. <https://pubmed.ncbi.nlm.nih.gov/1658348/> — [`KH`·`KHW`: H⁺ 생성과 혈류 의존적 제거]
20. Karmazyn M, Gan XT, Humphreys RA, Yoshida H, Kusumoto K. **The myocardial Na⁺-H⁺ exchange: structure, regulation, and its role in heart disease.** *Circ Res* 1999;85:777-786. <https://pubmed.ncbi.nlm.nih.gov/10532945/> — [`NHE1 → Na⁺ → 역방향 NCX → Ca²⁺`: 재관류 Ca²⁺ 급증(`KCA_RP`)의 기전]

## 4 · 재관류 손상 — 숙신산·ROS·mPTP (Reperfusion Injury) — 모델의 2번 경주

21. ★ Chouchani ET, Pell VR, Gaude E, et al. **Ischaemic accumulation of succinate controls reperfusion injury through mitochondrial ROS.** *Nature* 2014;515:431-435. <https://pubmed.ncbi.nlm.nih.gov/25383517/> — [`SUC` 상태변수와 `KSUC`/`KSUCO`/`KROS`. 이 논문이 없으면 산화제 폭발이 "재관류 사건"이 되는 이유를 쓸 수 없음]
22. ★ Halestrap AP. **A pore way to die: the role of mitochondria in reperfusion injury and cardioprotection.** *Biochem Soc Trans* 2010;38:841-860. <https://pubmed.ncbi.nlm.nih.gov/20658967/> — [mPTP의 Ca²⁺·ROS·pH 삼중 의존성 = `dP/dt` 식의 세 인자]
23. ★ Lemasters JJ, Bond JM, Chacon E, et al. **The pH paradox in ischemia-reperfusion injury to cardiac myocytes.** *EXS* 1996;76:99-114. <https://pubmed.ncbi.nlm.nih.gov/8805791/> — [`PHGATE = 1/(1+H/HP50)`. 모델 전체에서 가장 중요한 한 줄]
24. Bond JM, Herman B, Lemasters JJ. **Protection by acidotic pH against anoxia/reoxygenation injury to rat neonatal cardiac myocytes.** *Biochem Biophys Res Commun* 1991;179:798-803. <https://pubmed.ncbi.nlm.nih.gov/1909130/> — [산성이 보호적이라는 직접 실험 — 실험 14의 `HP50` 절단 대조]
25. Griffiths EJ, Halestrap AP. **Mitochondrial non-specific pores remain closed during cardiac ischaemia, but open upon reperfusion.** *Biochem J* 1995;307:93-98. <https://pubmed.ncbi.nlm.nih.gov/7717999/> — [공극이 허혈 중이 아니라 재관류 시에 열린다는 결정적 관찰]
26. Baines CP, Kaiser RA, Purcell NH, et al. **Loss of cyclophilin D reveals a critical role for mitochondrial permeability transition in cell death.** *Nature* 2005;434:658-662. <https://pubmed.ncbi.nlm.nih.gov/15800627/> — [`CYPD` 노드와 CsA의 분자 표적]
27. Piper HM, García-Dorado D, Ovize M. **A fresh look at reperfusion injury.** *Cardiovasc Res* 1998;38:291-300. <https://pubmed.ncbi.nlm.nih.gov/9709389/> — [과수축 경직(`CONTRACT`)이 독립적 사망 경로라는 논거]
28. Zweier JL, Talukder MAH. **The role of oxidants and free radicals in reperfusion injury.** *Cardiovasc Res* 2006;70:181-190. <https://pubmed.ncbi.nlm.nih.gov/16580655/> — [ROS 공급원의 시간적 분리: 미토콘드리아(분) vs 백혈구(시간) — `KROS` vs `KROS_N`]
29. Hausenloy DJ, Yellon DM. **New directions for protecting the heart against ischaemia-reperfusion injury: targeting the Reperfusion Injury Salvage Kinase (RISK)-pathway.** *Cardiovasc Res* 2004;61:448-460. <https://pubmed.ncbi.nlm.nih.gov/14962476/> — [`RISK` 노드]
30. Lecour S. **Activation of the protective Survivor Activating Factor Enhancement (SAFE) pathway against reperfusion injury.** *J Mol Cell Cardiol* 2009;47:32-40. <https://pubmed.ncbi.nlm.nih.gov/19269290/> — [`SAFE` 노드]
31. Murry CE, Jennings RB, Reimer KA. **Preconditioning with ischemia: a delay of lethal cell injury in ischemic myocardium.** *Circulation* 1986;74:1124-1136. <https://pubmed.ncbi.nlm.nih.gov/3769170/> — [`PRECOND`/`PRECOND_EFF` 공변량]
32. Heusch G, Gersh BJ. **The pathophysiology of acute myocardial infarction and strategies of protection beyond reperfusion: a continual challenge.** *Eur Heart J* 2017;38:774-784. <https://pubmed.ncbi.nlm.nih.gov/27354052/> — [심근보호 임상시험이 반복 실패한 구조적 이유 — 실험 5의 논지와 직접 대응]

## 5 · 미세혈관 폐색 / No-Reflow (Microvascular Obstruction)

33. ★ Kloner RA, Ganote CE, Jennings RB. **The "no-reflow" phenomenon after temporary coronary occlusion in the dog.** *J Clin Invest* 1974;54:1496-1508. <https://pubmed.ncbi.nlm.nih.gov/4140198/> — [`MVO` 상태변수의 원 관찰: 심외막이 열려도 조직은 재관류되지 않는다]
34. ★ de Waha S, Patel MR, Granger CB, et al. **Relationship between microvascular obstruction and adverse events following primary percutaneous coronary intervention: an individual patient data pooled analysis from seven randomized trials.** *Eur Heart J* 2017;38:3502-3510. <https://pubmed.ncbi.nlm.nih.gov/29020248/> — [MVO가 경색 크기와 **독립적으로** 예후를 예측 — 실험 7이 재현해야 하는 사실]
35. Niccoli G, Scalone G, Lerman A, Crea F. **Coronary microvascular obstruction in acute myocardial infarction.** *Eur Heart J* 2016;37:1024-1033. <https://pubmed.ncbi.nlm.nih.gov/26364289/> — [MVO의 네 구성요소(내피 팽윤·백혈구 마개·미세혈전·색전) = 모델의 네 항]
36. Fearon WF, Low AF, Yong AS, et al. **Prognostic value of the index of microcirculatory resistance measured after primary percutaneous coronary intervention.** *Circulation* 2013;127:2436-2441. <https://pubmed.ncbi.nlm.nih.gov/23681066/> — [`IMRN` 노드의 임상적 대응물]
37. Bekkers SC, Yazdani SK, Virmani R, Waltenberger J. **Microvascular obstruction: underlying pathophysiology and clinical diagnosis.** *J Am Coll Cardiol* 2010;55:1649-1660. <https://pubmed.ncbi.nlm.nih.gov/20394867/> — [심근내 출혈(`IMH`)과 철 침착의 후기 역할]
38. Ito H, Maruyama A, Iwakura K, et al. **Clinical implications of the 'no reflow' phenomenon. A predictor of complications and left ventricular remodeling in reperfused anterior wall myocardial infarction.** *Circulation* 1996;93:223-228. <https://pubmed.ncbi.nlm.nih.gov/8548892/> — [MVO → 재형성 연결(실험 7의 6개월 EDV)]

## 6 · 선천면역과 흉터 형성 (Innate Immunity & Scar Formation)

39. ★ Frangogiannis NG. **The inflammatory response in myocardial injury, repair, and remodelling.** *Nat Rev Cardiol* 2014;11:255-265. <https://pubmed.ncbi.nlm.nih.gov/24663091/> — [염증→수복 2상 프로그램: 모델의 M1→M2 전환 구조]
40. ★ Nahrendorf M, Swirski FK, Aikawa E, et al. **The healing myocardium sequentially mobilizes two monocyte subsets with divergent and complementary functions.** *J Exp Med* 2007;204:3037-3047. <https://pubmed.ncbi.nlm.nih.gov/18025128/> — [Ly6Cʰⁱ→Ly6Cˡᵒ 순차 동원 = `KSW`와 `switch` 함수]
40b. Horckmans M, Ring L, Duchene J, et al. **Neutrophils orchestrate post-myocardial infarction healing by polarizing macrophages towards a reparative phenotype.** *Eur Heart J* 2017;38:187-197. <https://pubmed.ncbi.nlm.nih.gov/28158426/> — [호중구 소실(efferocytosis)이 M2 전환의 방아쇠 — `switch = f(1/NEU)`의 근거]
41. Toldo S, Abbate A. **The NLRP3 inflammasome in acute myocardial infarction.** *Nat Rev Cardiol* 2018;15:203-214. <https://pubmed.ncbi.nlm.nih.gov/29143812/> — [`NLRP3 → IL-1β` 축과 콜히친/anakinra의 표적]
42. Frangogiannis NG. **The extracellular matrix in myocardial injury, repair, and remodeling.** *J Clin Invest* 2017;127:1600-1612. <https://pubmed.ncbi.nlm.nih.gov/28459429/> — [MMP/TIMP·콜라겐 성숙 = `MMP`/`COL`/`SCARSTR`]
43. Ducharme A, Frantz S, Aikawa M, et al. **Targeted deletion of matrix metalloproteinase-9 attenuates left ventricular enlargement and collagen accumulation after experimental myocardial infarction.** *J Clin Invest* 2000;106:55-62. <https://pubmed.ncbi.nlm.nih.gov/10880046/> — [MMP-9이 확장을 매개한다는 유전학적 증거 — `KTH`의 MMP 의존성]
44. Cleutjens JP, Kandala JC, Guarda E, Guntaka RV, Weber KT. **Regulation of collagen degradation in the rat myocardium after infarction.** *J Mol Cell Cardiol* 1995;27:1281-1292. <https://pubmed.ncbi.nlm.nih.gov/8531210/> — [흉터 성숙의 시간 척도(약 2주) = `KCO`/`KSC`]
45. Hammerman H, Kloner RA, Hale S, Schoen FJ, Braunwald E. **Dose-dependent effects of short-term methylprednisolone on myocardial infarct extent, scar formation, and ventricular function.** *Circulation* 1983;68:446-452. <https://pubmed.ncbi.nlm.nih.gov/6872182/> — [★ 항염이 과하면 흉터가 늦어져 확장이 악화 — 실험 11의 U자 반응이 재현해야 하는 고전]
46. Ridker PM, Everett BM, Thuren T, et al. **Antiinflammatory therapy with canakinumab for atherosclerotic disease (CANTOS).** *N Engl J Med* 2017;377:1119-1131. <https://pubmed.ncbi.nlm.nih.gov/28845751/> — [`ANTIIL1` PK/PD의 임상 앵커]
47. Tardif JC, Kouz S, Waters DD, et al. **Efficacy and safety of low-dose colchicine after myocardial infarction (COLCOT).** *N Engl J Med* 2019;381:2497-2505. <https://pubmed.ncbi.nlm.nih.gov/31733140/> — [`COLCH` 0.5 mg/일의 임상 근거]
48. Broch K, Anstensrud AK, Woxholt S, et al. **Randomized trial of interleukin-6 receptor inhibition in patients with acute ST-segment elevation myocardial infarction (ASSAIL-MI).** *J Am Coll Cardiol* 2021;77:1845-1855. <https://pubmed.ncbi.nlm.nih.gov/33858621/> — [`ANTIIL6`이 구제 지수를 개선했다는 보고]
49. Abbate A, Trankle CR, Buckley LF, et al. **Interleukin-1 blockade inhibits the acute inflammatory response in patients with ST-segment-elevation myocardial infarction.** *J Am Heart Assoc* 2020;9:e014941. <https://pubmed.ncbi.nlm.nih.gov/32122219/> — [anakinra의 CRP 억제 크기 — `KC1`/`KC2` 검증]

## 7 · 재형성과 Laplace 분기 (Remodelling & the Laplace Bifurcation) — 모델의 분기

50. ★ Pfeffer MA, Braunwald E. **Ventricular remodeling after myocardial infarction. Experimental observations and clinical implications.** *Circulation* 1990;81:1161-1172. <https://pubmed.ncbi.nlm.nih.gov/2138525/> — [경색 크기 → 확장 → 심부전의 연쇄. 모델의 분기 구조 전체]
51. ★ Grossman W, Jones D, McLaurin LP. **Wall stress and patterns of hypertrophy in the human left ventricle.** *J Clin Invest* 1975;56:56-64. <https://pubmed.ncbi.nlm.nih.gov/124746/> — [σ = P·r/(2h), 구심성 비후가 벽응력을 정상화한다 = 모델의 유일한 부귀환]
52. ★ Hutchins GM, Bulkley BH. **Infarct expansion versus extension: two different complications of acute myocardial infarction.** *Am J Cardiol* 1978;41:1127-1132. <https://pubmed.ncbi.nlm.nih.gov/665521/> — [`THIN`/`EXPANSION`을 `NI`/`NR`과 구별해서 따로 둔 이유]
53. Erlebacher JA, Weiss JL, Weisfeldt ML, Bulkley BH. **Early dilation of the infarcted segment in acute transmural myocardial infarction: role of infarct expansion in acute left ventricular enlargement.** *J Am Coll Cardiol* 1984;4:201-208. <https://pubmed.ncbi.nlm.nih.gov/6736460/> — [1일차 확장이 재형성과 다른 기전이라는 관찰 = `EDV`를 `EDVS`와 분리한 근거]
54. White HD, Norris RM, Brown MA, Brandt PW, Whitlock RM, Wild CJ. **Left ventricular end-systolic volume as the major determinant of survival after recovery from myocardial infarction.** *Circulation* 1987;76:44-51. <https://pubmed.ncbi.nlm.nih.gov/3594774/> — [ESV가 EF보다 강한 예후 인자 — 모델이 ESV를 명시적으로 계산하는 이유]
55. Bolognese L, Neskovic AN, Parodi G, et al. **Left ventricular remodeling after primary coronary angioplasty: patterns of left ventricular dilation and long-term prognostic implications.** *Circulation* 2002;106:2351-2357. <https://pubmed.ncbi.nlm.nih.gov/12403666/> — [재관류 시대의 확장 표현형 분포 — 실험 8의 "안정/서행/발산" 3분류의 임상 대응]
56. Gaudron P, Eilles C, Kugler I, Ertl G. **Progressive left ventricular dysfunction and remodeling after myocardial infarction. Potential mechanisms and early predictors.** *Circulation* 1993;87:755-763. <https://pubmed.ncbi.nlm.nih.gov/8443896/> — [발산군과 비발산군이 실제로 나뉜다는 관측 = 분기의 임상적 증거]
57. Sutton MG, Sharpe N. **Left ventricular remodeling after myocardial infarction: pathophysiology and therapy.** *Circulation* 2000;101:2981-2988. <https://pubmed.ncbi.nlm.nih.gov/10869273/> — [종합 리뷰 — 클러스터 11의 구성]
58. Konstam MA, Kramer DG, Patel AR, Maron MS, Udelson JE. **Left ventricular remodeling in heart failure: current concepts in clinical significance and assessment.** *JACC Cardiovasc Imaging* 2011;4:98-108. <https://pubmed.ncbi.nlm.nih.gov/21232712/> — [`EDVI` 종말점의 정의와 임상적 의미]

## 8 · 재관류 전략 — 시간, 용해제, PCI (Reperfusion Strategy)

59. ★ Nallamothu BK, Bates ER. **Percutaneous coronary intervention versus fibrinolytic therapy in acute myocardial infarction: is timing (almost) everything?** *Am J Cardiol* 2003;92:824-826. <https://pubmed.ncbi.nlm.nih.gov/14516884/> — [PCI 관련 지연 60분 역전점 — 실험 6이 계산으로 재현해야 하는 값]
60. ★ Pinto DS, Kirtane AJ, Nallamothu BK, et al. **Hospital delays in reperfusion for ST-elevation myocardial infarction: implications when selecting a reperfusion strategy.** *Circulation* 2006;114:2019-2025. <https://pubmed.ncbi.nlm.nih.gov/17075010/> — [역전점이 환자 특성에 따라 이동한다는 관찰]
61. Keeley EC, Boura JA, Grines CL. **Primary angioplasty versus intravenous thrombolytic therapy for acute myocardial infarction: a quantitative review of 23 randomised trials.** *Lancet* 2003;361:13-20. <https://pubmed.ncbi.nlm.nih.gov/12517460/> — [PCI 우월성의 크기 — 실험 6의 좌측 끝]
62. Boersma E, Maas AC, Deckers JW, Simoons ML. **Early thrombolytic treatment in acute myocardial infarction: reappraisal of the golden hour.** *Lancet* 1996;348:771-775. <https://pubmed.ncbi.nlm.nih.gov/8813982/> — [★ 시간–편익 곡선이 직선이 아니라 쌍곡선이라는 근거 — 실험 2의 핵심]
63. Steg PG, Bonnefoy E, Chabaud S, et al. **Impact of time to treatment on mortality after prehospital fibrinolysis or primary angioplasty (CAPTIM).** *Circulation* 2003;108:2851-2856. <https://pubmed.ncbi.nlm.nih.gov/14657226/> — [2시간 이내에서는 병원전 용해가 PCI에 뒤지지 않음]
64. Armstrong PW, Gershlick AH, Goldstein P, et al. **Fibrinolysis or primary PCI in ST-segment elevation myocardial infarction (STREAM).** *N Engl J Med* 2013;368:1379-1387. <https://pubmed.ncbi.nlm.nih.gov/23473396/> — [약물–침습 병용 전략 = 시나리오 S4]
65. Van de Werf F, Ristić AD, Averkov OV, et al. **STREAM-2: half-dose tenecteplase or primary percutaneous coronary intervention in older patients with ST-elevation myocardial infarction.** *Circulation* 2023;148:753-764. <https://pubmed.ncbi.nlm.nih.gov/37439213/> — [고령 반량 용해 — 실험 13의 임상 대응]
66. Assessment of the Safety and Efficacy of a New Thrombolytic (ASSENT-2) Investigators. **Single-bolus tenecteplase compared with front-loaded alteplase in acute myocardial infarction.** *Lancet* 1999;354:716-722. <https://pubmed.ncbi.nlm.nih.gov/10475182/> — [TNK 용량·출혈 프로파일]
67. Tanswell P, Modi N, Combs D, Danays T. **Pharmacokinetics and pharmacodynamics of tenecteplase in fibrinolytic therapy of acute myocardial infarction.** *Clin Pharmacokinet* 2002;41:1229-1245. <https://pubmed.ncbi.nlm.nih.gov/12452737/> — [`TNK_V1`·`TNK_CL`·`TNK_Q`·`TNK_V2` 전부]
68. Fitzgerald DJ, Catella F, Roy L, FitzGerald GA. **Marked platelet activation in vivo after intravenous streptokinase in patients with acute myocardial infarction.** *Circulation* 1988;77:142-150. <https://pubmed.ncbi.nlm.nih.gov/3121207/> — [★ `KPLT_PLN`: 플라스민이 혈소판을 활성화한다 — 용해 후 재폐색의 기전]
69. Ohman EM, Califf RM, Topol EJ, et al. **Consequences of reocclusion after successful reperfusion therapy in acute myocardial infarction (TAMI).** *Circulation* 1990;82:781-791. <https://pubmed.ncbi.nlm.nih.gov/2118195/> — [재폐색 발생률과 대가]

## 9 · 항혈소판·항응고 약리 (Antiplatelet & Anticoagulant Pharmacology)

70. Wallentin L, Becker RC, Budaj A, et al. **Ticagrelor versus clopidogrel in patients with acute coronary syndromes (PLATO).** *N Engl J Med* 2009;361:1045-1057. <https://pubmed.ncbi.nlm.nih.gov/19717846/> — [`P2Y12I` 선택과 임상 효과]
71. Teng R, Butler K. **Pharmacokinetics, pharmacodynamics, tolerability and safety of single ascending doses of ticagrelor.** *Eur J Clin Pharmacol* 2010;66:487-496. <https://pubmed.ncbi.nlm.nih.gov/20091161/> — [`TIC_KA`·`TIC_F`·`TIC_V`·`TIC_CL`·`TIC_EC50`]
72. Sabatine MS, Cannon CP, Gibson CM, et al. **Addition of clopidogrel to aspirin and fibrinolytic therapy for myocardial infarction with ST-segment elevation (CLARITY-TIMI 28).** *N Engl J Med* 2005;352:1179-1189. <https://pubmed.ncbi.nlm.nih.gov/15758000/> — [★ 용해제에 P2Y₁₂를 더하면 개통이 유지된다 — 실험 6의 재폐색 표]
73. ISIS-2 Collaborative Group. **Randomised trial of intravenous streptokinase, oral aspirin, both, or neither among 17,187 cases of suspected acute myocardial infarction.** *Lancet* 1988;332:349-360. <https://pubmed.ncbi.nlm.nih.gov/2899772/> — [아스피린과 용해제의 상가적 효과 — `ASA_EFF`]
74. Stone GW, Witzenbichler B, Guagliumi G, et al. **Bivalirudin during primary PCI in acute myocardial infarction (HORIZONS-AMI).** *N Engl J Med* 2008;358:2218-2230. <https://pubmed.ncbi.nlm.nih.gov/18499566/> — [`HEPARIN` 노드의 대안과 출혈 상충]

## 10 · 심근보호·요구 감소·GDMT (Cardioprotection, Demand Reduction & GDMT)

75. ★ Piot C, Croisille P, Staat P, et al. **Effect of cyclosporine on reperfusion injury in acute myocardial infarction.** *N Engl J Med* 2008;359:473-481. <https://pubmed.ncbi.nlm.nih.gov/18669426/> — [CsA 예비연구의 양성 결과 — 실험 5의 상한]
76. ★ Cung TT, Morel O, Cayla G, et al. **Cyclosporine before PCI in patients with acute myocardial infarction (CIRCUS).** *N Engl J Med* 2015;373:1021-1031. <https://pubmed.ncbi.nlm.nih.gov/26321103/> — [확증시험은 중립. 모델은 이 결과와 §12에서 명시적으로 대면합니다]
77. Ottani F, Latini R, Staszewsky L, et al. **Cyclosporine A in reperfused myocardial infarction: the multicenter, controlled, open-label CYCLE trial.** *J Am Coll Cardiol* 2016;67:365-374. <https://pubmed.ncbi.nlm.nih.gov/26821623/> — [두 번째 중립 결과]
78. ★ Ibanez B, Macaya C, Sánchez-Brunete V, et al. **Effect of early metoprolol on infarct size in ST-segment-elevation myocardial infarction patients undergoing primary percutaneous coronary intervention (METOCARD-CNIC).** *Circulation* 2013;128:1495-1503. <https://pubmed.ncbi.nlm.nih.gov/24002794/> — [재관류 **전** 정맥 베타차단제가 경색을 줄인다 — 실험 10이 요구 감소만으로 재현해야 하는 결과]
79. Roolvink V, Ibáñez B, Ottervanger JP, et al. **Early intravenous beta-blockers in patients with ST-segment elevation myocardial infarction before primary percutaneous coronary intervention (EARLY-BAMI).** *J Am Coll Cardiol* 2016;67:2705-2715. <https://pubmed.ncbi.nlm.nih.gov/27050189/> — [같은 개입이 더 늦게 주어지면 중립 — 실험 10의 시간 의존성]
80. Pfeffer MA, Braunwald E, Moyé LA, et al. **Effect of captopril on mortality and morbidity in patients with left ventricular dysfunction after myocardial infarction (SAVE).** *N Engl J Med* 1992;327:669-677. <https://pubmed.ncbi.nlm.nih.gov/1386652/> — [ACEi가 확장을 늦춘다 — 실험 9]
81. Pitt B, Remme W, Zannad F, et al. **Eplerenone, a selective aldosterone blocker, in patients with left ventricular dysfunction after myocardial infarction (EPHESUS).** *N Engl J Med* 2003;348:1309-1321. <https://pubmed.ncbi.nlm.nih.gov/12668699/> — [`MRA` 효과 크기]
82. Montalescot G, Pitt B, Lopez de Sa E, et al. **Early eplerenone treatment in patients with acute ST-elevation myocardial infarction without heart failure (REMINDER).** *Eur Heart J* 2014;35:2295-2302. <https://pubmed.ncbi.nlm.nih.gov/24780614/> — [MRA 조기 투여]
83. Butler J, Jones WS, Udell JA, et al. **Empagliflozin after acute myocardial infarction (EMPACT-MI).** *N Engl J Med* 2024;390:1455-1466. <https://pubmed.ncbi.nlm.nih.gov/38587239/> — [`SGLT2` — 1차 종말점 중립이지만 심부전 입원은 감소. §12에서 논의]
84. Pfeffer MA, Claggett B, Lewis EF, et al. **Angiotensin receptor-neprilysin inhibition in acute myocardial infarction (PARADISE-MI).** *N Engl J Med* 2021;385:1845-1855. <https://pubmed.ncbi.nlm.nih.gov/34758252/> — [`ARNI` 플래그]
85. Hofmann R, James SK, Jernberg T, et al. **Oxygen therapy in suspected acute myocardial infarction (DETO2X-AMI).** *N Engl J Med* 2017;377:1240-1249. <https://pubmed.ncbi.nlm.nih.gov/28844200/> — [`O2SUPP`: 저산소가 없으면 산소는 무익 — 모델에서 O₂DEL이 이미 정상인 이유]
86. Møller JE, Engstrøm T, Jensen LO, et al. **Microaxial flow pump or standard care in infarct-related cardiogenic shock (DanGer Shock).** *N Engl J Med* 2024;390:1382-1393. <https://pubmed.ncbi.nlm.nih.gov/38587239/> — [`MCS` 노드]
87. Byrne RA, Rossello X, Coughlan JJ, et al. **2023 ESC Guidelines for the management of acute coronary syndromes.** *Eur Heart J* 2023;44:3720-3826. <https://pubmed.ncbi.nlm.nih.gov/37622654/> — [시나리오 S1–S7의 임상적 타당성 근거]

## 11 · 표지자, 영상, 종말점 (Biomarkers, Imaging & Endpoints)

88. ★ Thygesen K, Alpert JS, Jaffe AS, et al. **Fourth universal definition of myocardial infarction (2018).** *Eur Heart J* 2019;40:237-269. <https://pubmed.ncbi.nlm.nih.gov/30165617/> — [경색 정의와 표지자 판정 기준]
89. Katus HA, Remppis A, Scheffold T, Diederich KW, Kuebler W. **Intracellular compartmentation of cardiac troponin T and its release kinetics in patients with reperfused and nonreperfused myocardial infarction.** *Am J Cardiol* 1991;67:1360-1367. <https://pubmed.ncbi.nlm.nih.gov/1904190/> — [★ 세척 현상: 재관류가 표지자 곡선의 **모양**을 바꾼다 — 실험 12]
90. Vasile VC, Jaffe AS. **High-sensitivity cardiac troponin for the diagnosis of patients with acute coronary syndromes.** *Curr Cardiol Rep* 2017;19:92. <https://pubmed.ncbi.nlm.nih.gov/28840493/> — [hs-cTn 동역학 — `KTND`]
91. Chia S, Senatore F, Raffel OC, Lee H, Wackers FJ, Jang IK. **Utility of cardiac biomarkers in predicting infarct size, left ventricular function, and clinical outcome after primary percutaneous coronary intervention.** *JACC Cardiovasc Interv* 2008;1:415-423. <https://pubmed.ncbi.nlm.nih.gov/19463339/> — [AUC vs 정점의 예측력 차이 — 실험 12의 결론]
92. Stone GW, Selker HP, Thiele H, et al. **Relationship between infarct size and outcomes following primary PCI: patient-level analysis from 10 randomized trials.** *J Am Coll Cardiol* 2016;67:1674-1683. <https://pubmed.ncbi.nlm.nih.gov/27056772/> — [★ 경색 크기 %LV와 사망·심부전의 정량적 연결 — 모델 종말점 스케일의 앵커]
93. Eitel I, de Waha S, Wöhrle J, et al. **Comprehensive prognosis assessment by CMR imaging after ST-segment elevation myocardial infarction.** *J Am Coll Cardiol* 2014;64:1217-1226. <https://pubmed.ncbi.nlm.nih.gov/25236513/> — [IS·MVO·MSI의 독립적 예후 기여]
94. Ibanez B, Aletras AH, Arai AE, et al. **Cardiac MRI endpoints in myocardial infarction experimental and clinical trials: JACC Scientific Expert Panel.** *J Am Coll Cardiol* 2019;74:238-256. <https://pubmed.ncbi.nlm.nih.gov/31236842/> — [구제 지수(MSI) 정의 — 모델의 `SALV` 계산식]
95. Francone M, Bucciarelli-Ducci C, Carbone I, et al. **Impact of primary coronary angioplasty delay on myocardial salvage, infarct size, and microvascular damage in patients with ST-segment elevation myocardial infarction.** *J Am Coll Cardiol* 2009;54:2145-2153. <https://pubmed.ncbi.nlm.nih.gov/19942086/> — [★ 사람에서 지연 시간별 MSI — 실험 2의 보정 표적 그 자체]

## 12 · 전기생리·합병증 (Electrophysiology & Complications)

96. Janse MJ, Wit AL. **Electrophysiological mechanisms of ventricular arrhythmias resulting from myocardial ischemia and infarction.** *Physiol Rev* 1989;69:1049-1169. <https://pubmed.ncbi.nlm.nih.gov/2678165/> — [손상 전류·회귀 — 클러스터 13]
97. Manning AS, Hearse DJ. **Reperfusion-induced arrhythmias: mechanisms and prevention.** *J Mol Cell Cardiol* 1984;16:497-518. <https://pubmed.ncbi.nlm.nih.gov/6086901/> — [`REPARR` 노드]
98. Hochman JS, Sleeper LA, Webb JG, et al. **Early revascularization in acute myocardial infarction complicated by cardiogenic shock (SHOCK).** *N Engl J Med* 1999;341:625-634. <https://pubmed.ncbi.nlm.nih.gov/10460813/> — [`SHOCK` 종말점]
99. Figueras J, Alcalde O, Barrabés JA, et al. **Changes in hospital mortality rates in 425 patients with acute ST-elevation myocardial infarction and cardiac rupture over a 30-year period.** *Circulation* 2008;118:2783-2789. <https://pubmed.ncbi.nlm.nih.gov/19064683/> — [`RUPT` 빈도와 흉터 시기의 관계]

## 13 · QSP 방법론과 도구 (QSP Methodology & Tooling)

100. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus machine learning, a fight worth fighting for the biological community?** *Biol Lett* 2018;14:20170660. <https://pubmed.ncbi.nlm.nih.gov/29769297/> — [기전 모델을 쓰는 이유]
101. Musante CJ, Ramanujan S, Schmidt BJ, Ghobrial OG, Lu J, Heatherington AC. **Quantitative systems pharmacology: a case for disease models.** *Clin Pharmacol Ther* 2017;101:24-27. <https://pubmed.ncbi.nlm.nih.gov/27709613/> — [질환 모델 라이브러리의 근거]
102. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019;8:883-893. <https://pubmed.ncbi.nlm.nih.gov/31674729/> — [mrgsolve 구현 관례]
103. Zhang XY, Trame MN, Lesko LJ, Schmidt S. **Sobol sensitivity analysis: a tool to guide the development and evaluation of systems pharmacology models.** *CPT Pharmacometrics Syst Pharmacol* 2015;4:69-79. <https://pubmed.ncbi.nlm.nih.gov/27548289/> — [실험 14의 구조적 민감도 분석 방법론적 배경]
104. Bradshaw EL, Spilker ME, Zang R, et al. **Applications of quantitative systems pharmacology in model-informed drug discovery: perspective on impact and opportunities.** *CPT Pharmacometrics Syst Pharmacol* 2019;8:777-791. <https://pubmed.ncbi.nlm.nih.gov/31535440/> — [QSP의 의사결정 활용]
105. Ottesen JT, Andreassen KV, Christiansen C, et al. **Age-related changes in the bone-vascular axis: a systems approach.** *J Theor Biol* 2016. <https://pubmed.ncbi.nlm.nih.gov/27208429/> — [다장기 QSP의 구조적 선례]

---

## 14 · 이 모델이 문헌과 어긋나는 지점 (Where This Model Disagrees With the Literature)

QSP 모델의 가치는 맞는 곳보다 **틀리는 곳을 명시하는 데** 있습니다. 아래는
`ami_reference_check.py`를 돌려 확인된 불일치이며, 매끄럽게 다듬지 않고
그대로 남겨 둡니다.

| # | 모델의 예측 | 문헌 | 어떻게 다루었는가 |
|---|---|---|---|
| 1 | mPTP 억제(CsA)가 재관류 직전 투여 시 경색을 줄인다 (실험 5) | Piot 2008 [75]은 양성, **CIRCUS [76]·CYCLE [77]은 중립** | 모델은 CsA가 "효과가 있다"를 주장하지 않습니다. 실험 5는 예측된 효과가 **허혈 시간에 대해 비단조**임을 보여 주며, 넓은 허혈 시간대를 등록한 실용시험이 스스로 효과크기를 희석한다는 **구조적 설명**을 제시합니다. 이것은 중립 결과의 변호가 아니라, 중립 결과와 양성 예비연구가 동시에 참일 수 있는 조건의 진술입니다. |
| 2 | SGLT2 억제제가 6개월 EDV를 줄인다 (실험 9) | **EMPACT-MI [83]의 1차 복합 종말점은 중립** (심부전 입원은 감소) | 모델의 SGLT2 효과는 용적(`GVOL_S`)과 확장 이득(`0.15×sglt`) 두 항뿐이며, 두 항 모두 **재형성 대리지표**에 작용합니다. 사망을 포함한 복합 종말점을 모델은 계산하지 않으므로, 모델과 EMPACT-MI는 서로 다른 것을 측정하고 있습니다. 모델이 임상시험을 "예측했다"고 읽어서는 안 됩니다. |
| 3 | 항IL-1β 강도를 높이면 흉터가 약해져 확장이 악화 (실험 11) | Hammerman 1983 [45]은 스테로이드에서 이를 확인, 그러나 **CANTOS [46]·COLCOT [47]에서는 흉터 약화 신호 없음** | 모델의 U자 반응은 IL-1β를 **거의 완전히** 차단할 때만 나타납니다(`CAN_EMAX > 0.95`). 임상 용량의 canakinumab/colchicine은 그 영역에 도달하지 않으며, 모델도 그 구간에서는 해를 예측하지 않습니다. 즉 U자의 오른쪽 가지는 **임상적으로 접근되지 않은 외삽**이며, 그렇게 표시했습니다. |
| 4 | 완전 폐색 24–48 h에서 AAR의 약 80%가 괴사, 심외막 테두리 생존 | Reimer 1977 [1]의 개 실험은 약 75%. 사람 데이터는 AAR 정의에 따라 60–90%로 넓게 분포 | 개 데이터에 맞추었습니다. 사람의 부수혈관은 개보다 빈약하므로 `COLL`을 낮추면(실험 3의 0.02–0.05행) 생존 테두리가 사라집니다. 종간 차이를 파라미터 한 개로 흡수한 것이며, 이는 **가정이지 결과가 아닙니다.** |
| 5 | 용해–PCI 역전점이 PCI 관련 지연 약 60–120 min | Nallamothu [59] 60 min, Pinto [60] 환자군에 따라 40–180 min | 일치하지만, 모델의 역전점은 `KLYS_PL`·`KREG`·`KPCI` 세 상수에 민감합니다. 이 값들은 개통률(90분 TIMI 3 ≈ 60% vs ≈ 95%)에 맞추어 역산한 것이므로, 역전점은 **독립적 예측이 아니라 개통률 데이터의 재표현**에 가깝습니다. 실험 6의 서술에서 이 점을 명시했습니다. |
| 6 | 사망률·심부전 입원 등 임상 사건 발생률 | Stone 2016 [92], Eitel 2014 [93] | **모델은 사건율을 계산하지 않습니다.** 클러스터 15의 종말점 노드는 지도에서 인과 경로를 보여 주기 위한 것이며, R 모델의 출력은 IS·MVO·EF·EDV·BNP·흉터강도까지입니다. 대리지표에서 사건율로 넘어가는 단계는 이 모델에 **구현되어 있지 않습니다.** |
| 7 | 실험 13의 두개내 출혈 지수 | 실제 ICH 위험 모형 (연령·체중·혈압·과거 뇌졸중) | 실험 13의 `ICH index`는 플라스민 AUC와 연령만으로 만든 **보고용 지표**이며 검증된 위험 예측 모형이 아닙니다. 절대 위험도로 읽어서는 안 됩니다. 코드와 출력에 그렇게 표시했습니다. |
| 8 | **베타차단제 단독**이 6개월 재형성을 **악화**시킨다 (실험 9: 성장 380%, 무치료 205%) | CAPRICORN [80 계열]·BHAT 등에서 경색 후 베타차단제는 사망과 재형성을 **개선** | **모델이 틀리는 지점입니다.** 이 모델에서 베타차단제의 유일한 경로는 심박수·수축력 감소이며, 그것이 심박출량을 낮추어 RAAS 구동을 키웁니다(Ang II 2.65 → 3.13). 급성 혈역학으로는 옳지만 만기 결과와는 반대입니다. 모델에 베타차단제의 **항부정맥·항산화·β1 수용체 상향회복·직접 역재형성** 효과가 없기 때문입니다. ACEi와 병용하면(ACEi 단독 61.7% → ACEi+BB 48.2%) 방향이 옳아지므로, 결함은 RAAS 차단이 없을 때의 단독 요법 궤적에 국한됩니다. 보정하지 않고 그대로 남겨 둡니다. |
| 8b | 정맥 베타차단제의 경색 감소 **크기** (모델 2.0%) | METOCARD-CNIC [78]은 약 20% 감소 (25.6 vs 32.0 g) | **방향과 시점은 맞지만 크기가 작습니다.** 이 모델에서 심박수 조절이 파면에 닿는 경로는 **확장기 연장 → 부수혈류 증가**(`KDIAST`) 하나뿐입니다. 기초 대사 요구는 심박수에 의존하지 않으므로(그렇게 두는 것이 정의상 옳습니다) 요구 감소가 이미 기초 이하인 조직을 구할 수는 없습니다. 또 Ibanez 그룹이 뒤에 보고한 **호중구–혈소판 응집 억제** 기전이 모델에 없습니다. `KDIAST`를 올려 20%에 맞출 수 있었지만, 그러면 시험 결과에 파라미터를 맞추는 것이 되므로 하지 않았습니다. |
| 8c | MVO가 경색 크기와 **독립적으로** 6개월 EDV를 예측하지 못한다 (MVO 3.17 vs 13.33 %LV에서 EDV 차이 2.5 mL) | de Waha 2017 [34]·Eitel 2014 [93]: MVO는 경색 크기 보정 후에도 독립적 예후 인자 | **모델이 재현하지 못하는 지점입니다.** MVO의 *발생*은 경색 크기와 잘 분리됩니다(MVO/경색 비 0.23–0.96). 문제는 *지속*입니다: `KMVO_R` = 0.012/h(반감기 약 2.4일)로 MVO가 재형성의 시간 척도에 닿기 전에 해소되고, 모델에는 **심근내 출혈의 철 침착**이나 **영구적 미세혈관 소실(capillary rarefaction)** 같은 비가역 경로가 없습니다. 그 경로를 넣으면 재현될 것으로 예상되지만, 데이터에 맞추기 위해 `KMVO_R`을 임의로 낮추지는 않았습니다. |
| 8d | **IL-1β 차단이 아무 효과를 내지 못한다** (투여 시점·차단 강도를 모두 훑어도 흉터·박화·EDV·EF 불변) | CANTOS [46]는 심혈관 사건 감소, COLCOT [47]도 감소, Abbate [49]는 CRP 억제 | **명백한 음성 결과이며 구조적 결함입니다.** 모델에서 흉터를 짓는 사슬은 M1 → M2 → TGF-β → 근섬유모세포 → 콜라겐이며 **IL-1β를 경유하지 않습니다.** IL-1β는 MMP 활성화의 20%만 기여하므로 차단해도 역학이 움직이지 않습니다. 임상 효과를 재현하려면 IL-1β가 (i) M1 유지·M2 전환 지연, (ii) 심근세포 음성 수축력, (iii) 경계 영역 세포자멸사에 직접 작용하는 경로가 필요합니다. 실험을 삭제하지 않고 **음성 결과로 보고**합니다. |
| 8e | 정점 트로포닌이 경색 크기를 **잘못된 순서로** 세우지 않는다 (39.7 → 50.2 → 53.5, 실제 경색과 동순) | Katus 1991 [89]·Chia 2008 [91]: 재관류된 경색은 더 작은 경색에서 더 높은 정점 | **시간 축은 재현하고 크기 축은 재현하지 못합니다.** 정점 *시각*은 재관류로 6.1 → 3.1시간 앞당겨집니다(옳음). 정점 *높이*는 실제 경색 크기 순서를 유지합니다(관찰과 다름). 또 AUC/경색이 135 → 68로 2배 변하므로 AUC도 큰 경색을 과소평가합니다. 원인은 `WASH_FLOOR = 0.15`로 추정합니다 — 폐색된 심근에서도 표지자의 15%가 유리되므로 대비가 충분히 벌어지지 않습니다. 맞추기 위해 낮추지 않았습니다. |
| 9 | 재형성 분리선의 **절대 위치** | 임계 경색 크기 약 18–20 %LV [50, 55, 56] | 모델의 **무치료** 분리선은 약 5 %LV, **전체 GDMT** 분리선은 약 18–20 %LV입니다. 인용되는 수치가 지침 치료 코호트에서 나온 것이므로 치료 분리선은 일치합니다. 그러나 무치료 분리선이 그만큼 낮은 것이 ACE 억제제 이전 시대의 실제 자연사인지, 아니면 `KDIL`/`KHYP_E` 비의 보정 오차인지 이 모델만으로는 구별할 수 없습니다. **두 해석을 구별하는 데이터로 보정하지 않았습니다.** |
| 10 | 발산 분기의 용적·BNP 수치 (EDV 600–1100 mL, BNP 70–730) | 실제 심실은 그렇게 커지지 않고, 그 전에 사망하거나 다른 기전이 개입 | 모델에는 **사망도, 다른 포화 기전도 없습니다.** `EDVMAX`(450 mL)는 생물학이 아니라 **수치적 천장**이며, 그 지점에서 적분을 멈춥니다. 분리선을 넘은 궤적은 "이 심실은 안정되지 않는다"는 **판정**을 보고하는 것이며 밀리리터가 아닙니다. 출력의 모든 해당 행에 `DIVERGED`와 도달 일자를 표시했습니다. |

## 15 · 모델이 의도적으로 포함하지 않은 것 (Deliberate Omissions)

- **관상동맥 해부의 공간 구조.** 위험 영역은 `AAR` 하나의 스칼라이며, LAD/RCA/LCx의 차이나 다혈관 질환은 없습니다.
- **우심실 경색.** 하벽 경색의 RV 침범과 그 특유의 전부하 의존적 혈역학은 없습니다.
- **당뇨·만성 신장병 등 동반질환.** `AGE`·`PRECOND`·`COLL` 세 공변량만 있습니다.
- **재경색과 판 진행.** 지도에는 있으나 ODE에는 없습니다.
- **혈소판 수·응고인자의 개체 변동**, 유전형(CYP2C19 등)에 따른 P2Y₁₂ 반응 차이.
- **부정맥의 동역학.** 클러스터 13은 지도에만 존재하고 ODE에는 없습니다 — VF는 사건이지 상태변수가 아니며, 이 모델의 시간 척도(시–월)와 맞지 않습니다.

이 목록은 모델을 쓸 때 **무엇을 물어봐서는 안 되는지**를 정하기 위한 것입니다.
