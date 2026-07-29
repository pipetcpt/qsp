# 만성 고칼륨혈증 QSP 모델 — 참고문헌 (References)

**Chronic hyperkalaemia in CKD and heart failure — the RAAS-inhibitor potassium dilemma**

이 문서는 `hk_qsp_model.dot`, `hk_mrgsolve_model.R`, `hk_reference_model.py`,
`hk_shiny_app.R` 에 사용된 **모든 구조적 가정과 파라미터의 출처**를 정리한 것입니다.
각 항목은 모델의 어느 부분에 쓰였는지 함께 표시합니다.

> 표기: **[F]** = 이 문헌으로 파라미터를 *적합(fit)* 함 · **[V]** = 적합하지 않고
> *검증(validation)* 에만 사용 · **[S]** = 구조(structure)의 근거 ·
> **[C]** = 임상 맥락(context)

---

## 1. 칼륨 항상성의 기본 구조 — 세포내외 분포와 완충 (Potassium homeostasis)

모델의 핵심 관계식
`K_total = Ce·V_ECF + Ci0·LAMrel·(Ce/Ce0)^α·V_ICF` 와 α = 0.25 의 근거.

1. **[S]** Palmer BF. *Regulation of Potassium Homeostasis.* Clin J Am Soc Nephrol. 2015;10(6):1050-60. — https://pubmed.ncbi.nlm.nih.gov/24721891/
2. **[S]** Palmer BF, Clegg DJ. *Physiology and Pathophysiology of Potassium Homeostasis: Core Curriculum 2019.* Am J Kidney Dis. 2019;74(5):682-95. — https://pubmed.ncbi.nlm.nih.gov/31227226/
3. **[S]** Youn JH, McDonough AA. *Recent advances in understanding integrative control of potassium homeostasis.* Annu Rev Physiol. 2009;71:381-401. — https://pubmed.ncbi.nlm.nih.gov/18759636/
4. **[V]** Sterns RH, Cox M, Feig PU, Singer I. *Internal potassium balance and the control of the plasma potassium concentration.* Medicine (Baltimore). 1981;60(5):339-54. — https://pubmed.ncbi.nlm.nih.gov/6268928/ — 혈청 K 3.0 에서 총체내 결핍 200-400 mmol. 모델은 이 값을 **보지 않고** 302 mmol 을 예측 (α 로부터 유도).
5. **[S]** Clausen T. *Quantification of Na+,K+ pumps and their transport rate in skeletal muscle: functional significance.* J Gen Physiol. 2013;142(4):327-45. — https://pubmed.ncbi.nlm.nih.gov/24081980/ — 골격근이 총 체내 K 의 ~75% 를 보유하며 완충의 주체.
6. **[S]** McDonough AA, Youn JH. *Potassium Homeostasis: The Knowns, the Unknowns, and the Health Benefits.* Physiology. 2017;32(2):100-11. — https://pubmed.ncbi.nlm.nih.gov/28202621/
7. **[S]** Gumz ML, Rabinowitz L, Wingo CS. *An Integrated View of Potassium Homeostasis.* N Engl J Med. 2015;373(1):60-72. — https://pubmed.ncbi.nlm.nih.gov/26132942/

## 2. 세포내 이동 조절 — 인슐린 · β2 · pH · 삼투압 (Transcellular shift)

`LAMrel = f_ins · f_β2 · f_pH · f_glu · f_aldo` 의 각 항.

8. **[F]** Nguyen TQ, Maalouf NM, Sakhaee K, Moe OW. *Comparison of insulin action on glucose versus potassium uptake in humans.* Clin J Am Soc Nephrol. 2011;6(7):1533-9. — https://pubmed.ncbi.nlm.nih.gov/21734082/ — 인슐린의 K 흡수는 포도당 흡수보다 낮은 농도에서 포화 → 기저 인슐린이 이미 대부분의 작용을 함 (모델의 `f_ins` 를 전 범위 함수로 둔 근거).
9. **[V]** Allon M, Copkney C. *Albuterol and insulin for treatment of hyperkalemia in hemodialysis patients.* Kidney Int. 1990;38(5):869-72. — https://pubmed.ncbi.nlm.nih.gov/2266671/ — 인슐린 -0.6 ~ -1.0, 살부타몰 -0.6 ~ -1.0 mmol/L.
10. **[V]** Ngugi NN, McLigeyo SO, Kayima JK. *Treatment of hyperkalaemia by altering the transcellular gradient in patients with renal failure.* East Afr Med J. 1997;74(8):503-9. — https://pubmed.ncbi.nlm.nih.gov/9487416/
11. **[S]** Aronson PS, Giebisch G. *Effects of pH on potassium: new explanations for old observations.* J Am Soc Nephrol. 2011;22(11):1981-9. — https://pubmed.ncbi.nlm.nih.gov/21980112/ — 무기산증만이 K 를 유의하게 이동시키며 유기산증(락트산·케톤산)은 그렇지 않다는 점.
12. **[S]** Adrogué HJ, Madias NE. *Changes in plasma potassium concentration during acute acid-base disturbances.* Am J Med. 1981;71(3):456-67. — https://pubmed.ncbi.nlm.nih.gov/7025622/
13. **[S]** Clausen T. *Hormonal and pharmacological modification of plasma potassium homeostasis.* Fundam Clin Pharmacol. 2010;24(5):595-605. — https://pubmed.ncbi.nlm.nih.gov/20618871/
14. **[C]** Kamel KS, Wei C. *Controversial issues in the treatment of hyperkalaemia.* Nephrol Dial Transplant. 2003;18(11):2215-8. — https://pubmed.ncbi.nlm.nih.gov/14551344/

## 3. 신장의 칼륨 처리와 원위 세뇨관 (Renal handling, the ASDN)

`E_ren = FD·filt + S_cap·RASDN·(Ce/(Km+Ce))·q^0.5·(HCO3/24)^n` 의 구조.

15. **[S]** Palmer BF, Clegg DJ. *Hyperkalemia across the Continuum of Kidney Function.* Clin J Am Soc Nephrol. 2018;13(1):155-7. — https://pubmed.ncbi.nlm.nih.gov/29114006/
16. **[F]** Hayes CP Jr, McLeod ME, Robinson RR. *An extravenal [sic] mechanism for the maintenance of potassium balance in severe chronic renal failure.* Trans Assoc Am Physicians. 1967;80:207-16. — https://pubmed.ncbi.nlm.nih.gov/6082243/ — 신부전에서 대장 분비의 상향 조절(모델의 `KC_COL`).
17. **[F]** Schultze RG, Taggart DD, Shapiro H, et al. *On the adaptation in potassium excretion associated with nephron reduction in the dog.* J Clin Invest. 1971;50(5):1061-8. — https://pubmed.ncbi.nlm.nih.gov/5552407/ — 잔존 네프론당 분비능 상향(모델의 `ADAPT_P`, `ADAPT_MX`).
18. **[S]** Welling PA. *Roles and Regulation of Renal K Channels.* Annu Rev Physiol. 2016;78:415-35. — https://pubmed.ncbi.nlm.nih.gov/26654186/ — ROMK 와 BK(Maxi-K)의 역할 분담.
19. **[S]** Wang WH, Giebisch G. *Regulation of potassium (K) handling in the renal collecting duct.* Pflugers Arch. 2009;458(1):157-68. — https://pubmed.ncbi.nlm.nih.gov/18839206/
20. **[S]** Terker AS, Zhang C, McCormick JA, et al. *Potassium modulates electrolyte balance and blood pressure through effects on distal cell voltage and chloride.* Cell Metab. 2015;21(1):39-50. — https://pubmed.ncbi.nlm.nih.gov/25565204/ — Kir4.1/5.1–WNK–SPAK–NCC 스위치(지도의 cluster 5).
21. **[S]** Cuevas CA, Su XT, Wang MX, et al. *Potassium Sensing by Renal Distal Tubules Requires Kir4.1.* J Am Soc Nephrol. 2017;28(6):1814-25. — https://pubmed.ncbi.nlm.nih.gov/28052988/
22. **[S]** Sorensen MV, Grossmann S, Roesinger M, et al. *Rapid dephosphorylation of the renal sodium chloride cotransporter in response to oral potassium intake in mice.* Kidney Int. 2013;83(5):811-24. — https://pubmed.ncbi.nlm.nih.gov/23447069/
23. **[S]** Preston RA, Afshartous D, Rodco R, et al. *Evidence for a gastrointestinal-renal kaliuretic signaling axis in humans.* Kidney Int. 2015;88(6):1383-91. — https://pubmed.ncbi.nlm.nih.gov/26308672/ — 장-신장 되먹임(지도의 `gutsensor`).
24. **[S]** Rabelink TJ, Koomans HA, Hené RJ, Dorhout Mees EJ. *Early and late adjustment to potassium loading in humans.* Kidney Int. 1990;38(5):942-7. — https://pubmed.ncbi.nlm.nih.gov/2266680/
25. **[S]** DuBose TD Jr. *Regulation of Potassium Homeostasis in CKD.* Adv Chronic Kidney Dis. 2017;24(5):305-14. — https://pubmed.ncbi.nlm.nih.gov/29031357/

## 4. 산-염기와 칼륨 (Acid-base)

모델에서 **알칼리 요법은 세포 이동이 아니라 신장 분비를 통해** 작용한다는 결론의 근거.

26. **[F]** de Brito-Ashurst I, Varagunam M, Raftery MJ, Yaqoob MM. *Bicarbonate supplementation slows progression of CKD and improves nutritional status.* J Am Soc Nephrol. 2009;20(9):2075-84. — https://pubmed.ncbi.nlm.nih.gov/19608703/
27. **[F]** Di Iorio BR, Bellasi A, Raphael KL, et al. (UBI Study) *Treatment of metabolic acidosis with sodium bicarbonate delays progression of chronic kidney disease.* J Nephrol. 2019;32(6):989-1001. — https://pubmed.ncbi.nlm.nih.gov/31598912/
28. **[S]** Wesson DE, Buysse JM, Bushinsky DA. *Mechanisms of Metabolic Acidosis-Induced Kidney Injury in CKD.* J Am Soc Nephrol. 2020;31(3):469-82. — https://pubmed.ncbi.nlm.nih.gov/31988269/
29. **[V]** Blumberg A, Weidmann P, Shaw S, Gnädinger M. *Effect of various therapeutic approaches on plasma potassium and major regulating factors in terminal renal failure.* Am J Med. 1988;85(4):507-12. — https://pubmed.ncbi.nlm.nih.gov/3052050/ — **중탄산염은 급성 고칼륨혈증 치료제가 아니다**는 고전적 음성 결과.
30. **[V]** Allon M, Shanklin N. *Effect of bicarbonate administration on plasma potassium in dialysis patients.* Am J Kidney Dis. 1996;28(4):508-14. — https://pubmed.ncbi.nlm.nih.gov/8840939/

## 5. RAAS 축과 무기질코르티코이드 수용체 (RAAS and the MR)

31. **[S]** Arroyo JP, Ronzaud C, Lagnaz D, et al. *Aldosterone paradox: differential regulation of ion transport in distal nephron.* Physiology. 2011;26(2):115-23. — https://pubmed.ncbi.nlm.nih.gov/21487030/ — 알도스테론 역설(Na 보존 vs K 분비)의 구조.
32. **[S]** Himathongkam T, Dluhy RG, Williams GH. *Potassim[sic]-aldosterone-renin interrelationships.* J Clin Endocrinol Metab. 1975;41(1):153-9. — https://pubmed.ncbi.nlm.nih.gov/1167307/ — 혈청 K 의 직접적 알도스테론 자극(모델의 `SK_ALDO`, 알도스테론 탈출 경로).
33. **[C]** DeFronzo RA. *Hyperkalemia and hyporeninemic hypoaldosteronism.* Kidney Int. 1980;17(1):118-34. — https://pubmed.ncbi.nlm.nih.gov/6990088/ — 제4형 신세뇨관성 산증.
34. **[F]** Pitt B, Zannad F, Remme WJ, et al. (RALES) *The effect of spironolactone on morbidity and mortality in patients with severe heart failure.* N Engl J Med. 1999;341(10):709-17. — https://pubmed.ncbi.nlm.nih.gov/10471456/ — 스피로노락톤 25 mg 의 K 상승 +0.30 mmol/L (모델의 `KI_MRA` 적합 앵커).
35. **[V]** Pitt B, Remme W, Zannad F, et al. (EPHESUS) *Eplerenone, a selective aldosterone blocker, in patients with left ventricular dysfunction after myocardial infarction.* N Engl J Med. 2003;348(14):1309-21. — https://pubmed.ncbi.nlm.nih.gov/12668699/
36. **[V]** Bakris GL, Agarwal R, Anker SD, et al. (FIDELIO-DKD) *Effect of Finerenone on Chronic Kidney Disease Outcomes in Type 2 Diabetes.* N Engl J Med. 2020;383(23):2219-29. — https://pubmed.ncbi.nlm.nih.gov/33264825/ — 관찰된 K 상승 +0.23; 모델은 가정된 MR 부하에서 +0.10 으로 **불일치**하며 이를 그대로 보고함.
37. **[V]** Pitt B, Filippatos G, Agarwal R, et al. (FIGARO-DKD) *Cardiovascular Events with Finerenone in Kidney Disease and Type 2 Diabetes.* N Engl J Med. 2021;385(24):2252-63. — https://pubmed.ncbi.nlm.nih.gov/34449181/
38. **[C]** Juurlink DN, Mamdani MM, Lee DS, et al. *Rates of hyperkalemia after publication of the Randomized Aldactone Evaluation Study.* N Engl J Med. 2004;351(6):543-51. — https://pubmed.ncbi.nlm.nih.gov/15295047/ — 시험 밖 인구에서의 고칼륨혈증 급증: 모델이 강조하는 "시험이 배제한 eGFR 구간" 문제의 실제 사례.

## 6. 칼륨 결합제 (Potassium binders)

39. **[V]** Weir MR, Bakris GL, Bushinsky DA, et al. (OPAL-HK) *Patiromer in patients with kidney disease and hyperkalemia receiving RAAS inhibitors.* N Engl J Med. 2015;372(3):211-21. — https://pubmed.ncbi.nlm.nih.gov/25415805/ — 4주 ΔK -1.01 (모델 -0.95, 적합하지 않음).
40. **[V]** Bakris GL, Pitt B, Weir MR, et al. (AMETHYST-DN) *Effect of Patiromer on Serum Potassium Level in Patients With Hyperkalemia and Diabetic Kidney Disease.* JAMA. 2015;314(2):151-61. — https://pubmed.ncbi.nlm.nih.gov/26172895/
41. **[V]** Kosiborod M, Rasmussen HS, Lavin P, et al. (HARMONIZE) *Effect of sodium zirconium cyclosilicate on potassium lowering for 28 days among outpatients with hyperkalemia.* JAMA. 2014;312(21):2223-33. — https://pubmed.ncbi.nlm.nih.gov/25402495/ — 5/10/15 g 유지용량에서 K 4.8/4.5/4.4 (모델 4.89/4.63/4.49).
42. **[V]** Packham DK, Rasmussen HS, Lavin PT, et al. (ZS-003) *Sodium zirconium cyclosilicate in hyperkalemia.* N Engl J Med. 2015;372(3):222-31. — https://pubmed.ncbi.nlm.nih.gov/25415807/
43. **[S]** Stavros F, Yang A, Leon A, et al. *Characterization of structure and function of ZS-9, a K selective ion trap.* PLoS One. 2014;9(12):e114686. — https://pubmed.ncbi.nlm.nih.gov/25531770/ — SZC 가 위/근위 장관에서 작용하는 구조적 근거(모델에서 patiromer 와 작용 부위를 다르게 둔 이유).
44. **[S]** Li L, Harrison SD, Cope MJ, Park C, et al. *Mechanism of Action and Pharmacology of Patiromer, a Nonabsorbed Cross-Linked Polymer That Lowers Serum Potassium Concentration.* J Cardiovasc Pharmacol Ther. 2016;21(5):456-65. — https://pubmed.ncbi.nlm.nih.gov/26856345/ — 원위 결장에서 Ca-K 교환.
45. **[C]** Agarwal R, Rossignol P, Romero A, et al. (AMBER) *Patiromer versus placebo to enable spironolactone use in patients with resistant hypertension and CKD.* Lancet. 2019;394(10208):1540-50. — https://pubmed.ncbi.nlm.nih.gov/31533906/ — "RAASi enablement" 개념의 무작위 근거.
46. **[C]** Butler J, Anker SD, Lund LH, et al. (DIAMOND) *Patiromer for the management of hyperkalaemia in heart failure with reduced ejection fraction.* Eur Heart J. 2022;43(41):4362-73. — https://pubmed.ncbi.nlm.nih.gov/35900838/
47. **[C]** Sterns RH, Rojas M, Bernstein P, Chennupati S. *Ion-exchange resins for the treatment of hyperkalemia: are they safe and effective?* J Am Soc Nephrol. 2010;21(5):733-5. — https://pubmed.ncbi.nlm.nih.gov/20167700/ — SPS(폴리스티렌설폰산)의 근거 부족과 대장 괴사 신호.
48. **[C]** Natale P, Palmer SC, Ruospo M, et al. *Potassium binders for chronic hyperkalaemia in people with CKD.* Cochrane Database Syst Rev. 2020;6:CD013165. — https://pubmed.ncbi.nlm.nih.gov/32588430/

## 7. 급성 고칼륨혈증 치료 (Acute management)

49. **[C]** Lindner G, Burdmann EA, Clase CM, et al. *Acute hyperkalemia in the emergency department: a summary from a Kidney Disease: Improving Global Outcomes conference.* Eur J Emerg Med. 2020;27(5):329-37. — https://pubmed.ncbi.nlm.nih.gov/32852924/
50. **[V]** Harel Z, Kamel KS. *Optimal Dose and Method of Administration of Intravenous Insulin in the Management of Emergency Hyperkalemia: A Systematic Review.* PLoS One. 2016;11(5):e0154963. — https://pubmed.ncbi.nlm.nih.gov/27148740/
51. **[C]** Coca A, Valencia AL, Bustamante J, et al. *Hypoglycemia following intravenous insulin plus glucose for hyperkalemia in patients with impaired renal function.* PLoS One. 2017;12(2):e0172961. — https://pubmed.ncbi.nlm.nih.gov/28245289/ — 모델이 **재현하지 못하는** 15-20% 저혈당 발생률(한계로 명시).
52. **[S]** Parham WA, Mehdirad AA, Biermann KM, Fredman CS. *Hyperkalemia revisited.* Tex Heart Inst J. 2006;33(1):40-7. — https://pubmed.ncbi.nlm.nih.gov/16572868/
53. **[S]** Diercks DB, Shumaik GM, Harrigan RA, et al. *Electrocardiographic manifestations: electrolyte abnormalities.* J Emerg Med. 2004;27(2):153-60. — https://pubmed.ncbi.nlm.nih.gov/15261358/
54. **[V]** Montague BT, Ouellette JR, Buller GK. *Retrospective review of the frequency of ECG changes in hyperkalemia.* Clin J Am Soc Nephrol. 2008;3(2):324-30. — https://pubmed.ncbi.nlm.nih.gov/18235147/ — 첨예 T 파의 민감도가 낮다는 점(모델이 T 파를 진단 지표로 쓰지 않는 이유).
55. **[S]** Weiss JN, Qu Z, Shivkumar K. *Electrophysiology of Hypokalemia and Hyperkalemia.* Circ Arrhythm Electrophysiol. 2017;10(3):e004667. — https://pubmed.ncbi.nlm.nih.gov/28314851/ — Em, Na 채널 가용성, 전도속도의 정량적 관계(모델의 cluster 12).

## 8. 역학과 임상 결과 (Epidemiology and outcomes)

56. **[F]** Kovesdy CP, Matsushita K, Sang Y, et al. (CKD Prognosis Consortium) *Serum potassium and adverse outcomes across the range of kidney function: a CKD Prognosis Consortium meta-analysis.* Eur Heart J. 2018;39(17):1535-42. — https://pubmed.ncbi.nlm.nih.gov/29554312/ — U 자형 위험곡선(모델의 `BK_HI`, `BK_LO` 적합 근거).
57. **[F]** Collins AJ, Pitt B, Reaven N, et al. *Association of Serum Potassium with All-Cause Mortality in Patients with and without Heart Failure, CKD, and/or Diabetes.* Am J Nephrol. 2017;46(3):213-21. — https://pubmed.ncbi.nlm.nih.gov/28866674/
58. **[C]** Einhorn LM, Zhan M, Hsu VD, et al. *The frequency of hyperkalemia and its significance in chronic kidney disease.* Arch Intern Med. 2009;169(12):1156-62. — https://pubmed.ncbi.nlm.nih.gov/19546417/
59. **[C]** Luo J, Brunelli SM, Jensen DE, Yang A. *Association between Serum Potassium and Outcomes in Patients with Reduced Kidney Function.* Clin J Am Soc Nephrol. 2016;11(1):90-100. — https://pubmed.ncbi.nlm.nih.gov/26500246/
60. **[C]** James G, Kim J, Mellström C, Ford KL, et al. *Serum potassium variability as a predictor of clinical outcomes in patients with cardiorenal disease or diabetes.* Clin Kidney J. 2022;15(4):758-70. — https://pubmed.ncbi.nlm.nih.gov/35371436/

## 9. RAASi 감량·중단의 대가 (The cost of down-titration) — 모델의 핵심 결과

61. **[F]** Epstein M, Reaven NL, Funk SE, et al. *Evaluation of the treatment gap between clinical guidelines and the utilization of renin-angiotensin-aldosterone system inhibitors.* Am J Manag Care. 2015;21(11 Suppl):S212-20. — https://pubmed.ncbi.nlm.nih.gov/26619183/ — 고칼륨혈증 후 RAASi 감량/중단과 사망률 증가의 연관.
62. **[V]** Xie X, Liu Y, Perkovic V, et al. *Renin-Angiotensin System Inhibitors and Kidney and Cardiovascular Outcomes in Patients With CKD: A Bayesian Network Meta-analysis.* Am J Kidney Dis. 2016;67(5):728-41. — https://pubmed.ncbi.nlm.nih.gov/26597926/ — 모델의 `LNHR_RD` 및 eGFR 기울기 효과.
63. **[V]** Brenner BM, Cooper ME, de Zeeuw D, et al. (RENAAL) *Effects of losartan on renal and cardiovascular outcomes in patients with type 2 diabetes and nephropathy.* N Engl J Med. 2001;345(12):861-9. — https://pubmed.ncbi.nlm.nih.gov/11565518/
64. **[V]** Fu EL, Evans M, Clase CM, et al. *Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes.* J Am Soc Nephrol. 2021;32(2):424-35. — https://pubmed.ncbi.nlm.nih.gov/33372009/
65. **[C]** Bhandari S, Mehta S, Khwaja A, et al. (STOP-ACEi) *Renin-Angiotensin System Inhibition in Advanced Chronic Kidney Disease.* N Engl J Med. 2022;387(22):2021-32. — https://pubmed.ncbi.nlm.nih.gov/36326117/ — 진행된 CKD 에서 RAASi 중단이 eGFR 을 개선하지 않았다는 무작위 결과 — 모델의 hazard 층이 관찰연구에 근거함을 상기시키는 반례.
66. **[C]** Leon SJ, Whitlock R, Rigatto C, et al. *Hyperkalemia-Related Discontinuation of Renin-Angiotensin-Aldosterone System Inhibitors and Clinical Outcomes in CKD.* Am J Kidney Dis. 2022;80(2):164-73. — https://pubmed.ncbi.nlm.nih.gov/35085685/

## 10. 식이 · 대장 분비 · 병용약물 (Diet, colon, co-medication)

67. **[C]** Cupisti A, Kovesdy CP, D'Alessandro C, Kalantar-Zadeh K. *Dietary Approach to Recurrent or Chronic Hyperkalaemia in Patients with Decreased Kidney Function.* Nutrients. 2018;10(3):261. — https://pubmed.ncbi.nlm.nih.gov/29495340/
68. **[C]** Bernier-Jean A, Wong G, Saglimbene V, et al. *Dietary Potassium Intake and All-Cause Mortality in Adults Treated with Hemodialysis.* Clin J Am Soc Nephrol. 2021;16(12):1851-61. — https://pubmed.ncbi.nlm.nih.gov/34853064/ — 식이 K 제한의 근거가 생각보다 약하다는 점.
69. **[C]** Neal B, Wu Y, Feng X, et al. (SSaSS) *Effect of Salt Substitution on Cardiovascular Events and Death.* N Engl J Med. 2021;385(12):1067-77. — https://pubmed.ncbi.nlm.nih.gov/34459569/ — 저나트륨 대체염 = KCl 부하(지도의 `saltsub`).
70. **[S]** Sandle GI, Gaiger E, Tapster S, Goodship TH. *Enhanced rectal potassium secretion in chronic renal insufficiency: evidence for large intestinal potassium adaptation in man.* Clin Sci (Lond). 1986;71(4):393-401. — https://pubmed.ncbi.nlm.nih.gov/3757437/
71. **[C]** Perazella MA. *Drug-induced hyperkalemia: old culprits and new offenders.* Am J Med. 2000;109(4):307-14. — https://pubmed.ncbi.nlm.nih.gov/10996582/
72. **[C]** Antoniou T, Gomes T, Juurlink DN, et al. *Trimethoprim-sulfamethoxazole-induced hyperkalemia in patients receiving inhibitors of the renin-angiotensin system.* Arch Intern Med. 2010;170(12):1045-9. — https://pubmed.ncbi.nlm.nih.gov/20585070/
73. **[C]** Neuen BL, Oshima M, Agarwal R, et al. *Sodium-Glucose Cotransporter 2 Inhibitors and Risk of Hyperkalemia in People With Type 2 Diabetes: A Meta-Analysis.* Circulation. 2022;145(19):1460-70. — https://pubmed.ncbi.nlm.nih.gov/35394821/ — SGLT2i 의 고칼륨혈증 위험 감소(모델의 `SGLT2I` 스위치).

## 11. 당뇨병성 케톤산증 — 같은 숫자, 반대의 저장고 (The partition patient)

74. **[C]** Kitabchi AE, Umpierrez GE, Miles JM, Fisher JN. *Hyperglycemic crises in adult patients with diabetes.* Diabetes Care. 2009;32(7):1335-43. — https://pubmed.ncbi.nlm.nih.gov/19564476/ — DKA 의 총 체내 K 결핍 3-5 mmol/kg 에도 혈청 K 는 정상~상승.
75. **[S]** Adrogué HJ, Lederer ED, Suki WN, Eknoyan G. *Determinants of plasma potassium levels in diabetic ketoacidosis.* Medicine (Baltimore). 1986;65(3):163-72. — https://pubmed.ncbi.nlm.nih.gov/3084904/
76. **[S]** Nicolis GL, Kahn T, Sanchez A, Gabrilove JL. *Glucose-induced hyperkalemia in diabetic subjects.* Arch Intern Med. 1981;141(1):49-53. — https://pubmed.ncbi.nlm.nih.gov/7447584/ — 고삼투압에 의한 K 이동(모델의 `K_GLU`).

## 12. 임상 진료지침 (Guidelines)

77. **[C]** Kidney Disease: Improving Global Outcomes (KDIGO) CKD Work Group. *KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease.* Kidney Int. 2024;105(4S):S117-S314. — https://pubmed.ncbi.nlm.nih.gov/38490803/
78. **[C]** Clase CM, Carrero JJ, Ellison DH, et al. *Potassium homeostasis and management of dyskalemia in kidney diseases: conclusions from a KDIGO Controversies Conference.* Kidney Int. 2020;97(1):42-61. — https://pubmed.ncbi.nlm.nih.gov/31706619/
79. **[C]** McDonagh TA, Metra M, Adamo M, et al. *2021 ESC Guidelines for the diagnosis and treatment of acute and chronic heart failure.* Eur Heart J. 2021;42(36):3599-3726. — https://pubmed.ncbi.nlm.nih.gov/34447992/
80. **[C]** Heidenreich PA, Bozkurt B, Aguilar D, et al. *2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure.* Circulation. 2022;145(18):e895-e1032. — https://pubmed.ncbi.nlm.nih.gov/35363499/

## 13. QSP 방법론 (Modelling methodology)

81. **[S]** Elmokadem A, Riggs MM, Baron KT. *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve.* CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883-93. — https://pubmed.ncbi.nlm.nih.gov/31652028/
82. **[S]** Guyton AC, Coleman TG, Granger HJ. *Circulation: overall regulation.* Annu Rev Physiol. 1972;34:13-46. — https://pubmed.ncbi.nlm.nih.gov/4334846/ — 체액 구획 부피(V_ECF 0.20·BW, V_ICF 0.36·BW)의 고전적 출처.
83. **[S]** Hallow KM, Gebremichael Y. *A quantitative systems physiology model of renal function and blood pressure regulation: Model description.* CPT Pharmacometrics Syst Pharmacol. 2017;6(6):383-92. — https://pubmed.ncbi.nlm.nih.gov/28548387/ — 신장 QSP 모델링의 구조적 참고.
84. **[S]** Hallow KM, Gebremichael Y. *A quantitative systems physiology model of renal function and blood pressure regulation: Application in salt-sensitive hypertension.* CPT Pharmacometrics Syst Pharmacol. 2017;6(6):393-400. — https://pubmed.ncbi.nlm.nih.gov/28556624/

---

## 부록 A — 파라미터별 출처 요약 (Parameter provenance)

| 파라미터 | 값 | 출처 유형 | 근거 |
|---|---|---|---|
| `ALPHA` | 0.25 | 구조 + 검증 | #4, #5 — 만성 완충 224 mmol/(mmol/L), 급성 66 |
| `S0` | 75.18 mmol/day | **적합** | eGFR 100 에서 K 4.20 (#1, #15) |
| `ADAPT_P` | 0.879 | **적합** | eGFR 20 에서 K 5.00 (#17, #58) |
| `ADAPT_MX` | 6.36 | **적합** | eGFR 12 에서 K 5.30 (#17, #25) |
| `KI_MRA` | 0.1029 mg/L | **적합** | RALES ΔK +0.30 (#34) |
| `N_HCO3` | 0.435 | **적합** | HCO3 18→24 에서 ΔK -0.30 (#26, #27) |
| `KC_COL` | 1.00 | 구조 | 신부전 대장 분비 상향 (#16, #70) |
| `PHIMAX_P` / `D50_PAT` | 0.42 / 12 g | 구조 | 작용 부위·기전 (#44); 검증은 #39 |
| `PHIMAX_S` / `D50_SZC` | 0.45 / 8 g | 구조 | 작용 부위·기전 (#43); 검증은 #41, #42 |
| `EMAX_INS` | 0.120 | 검증 정합 | 인슐린 ΔK -0.6~-1.0 (#9, #50) |
| `EMAX_B2` | 0.070 | 검증 정합 | 살부타몰 ΔK -0.6~-1.0 (#9) |
| `K_PH` | 0.190 | 구조 | 무기산증 ΔK +0.3/0.1 pH (#11, #12) |
| `K_GLU` | 0.0025 | 구조 | 고혈당 유발 K 이동 (#76) |
| `SK_ALDO` | 1.30 | 구조 | K 의 직접적 알도스테론 자극 (#32) |
| `BK_HI`,`BK_LO` | 0.642 / 0.916 | **적합** | U 자형 사망 위험 (#56, #57) |
| `LNHR_RD` | -0.248 | 문헌값 | RAASi HR 0.78 (#62, #63) |
| `SLOPE_RD` | 1.40 mL/min/yr | 문헌값 | (#62, #63) |
| `VTH0`,`H_MID`,`S_CA` | -70, -78, 15 | 구조 | 막전위·Na 채널 가용성 (#55) |

## 부록 B — 모델이 다루지 않는 것 (Explicitly out of scope)

아래는 임상적으로 중요하지만 이 모델에 **의도적으로 포함하지 않은** 것들이며,
각각의 참고문헌을 함께 둡니다. 모델을 이 영역에 사용해서는 안 됩니다.

- **급성 신손상(AKI on CKD)** — 실제 중증 고칼륨혈증의 최다 원인 (#49)
- **횡문근융해증·종양용해증** — ICF 저장고 자체로부터의 K 방출; 본 모델은 ICF 를 수동적 저장고로만 취급 (#52)
- **가성 고칼륨혈증** — 용혈·백혈구증가·혈소판증가 (#52)
- **투석** — 간헐적·초고청소율 제거항 (#78)
- **디곡신 중독** — Na-K-ATPase 직접 억제 (#71)
- **운동 유발 일과성 고칼륨혈증** — 작업근 간질강 축적 (#5)
- **소아·임신** — 파라미터가 성인 70 kg 기준

## 부록 C — 인용 수 (Citation count)

총 **84편**의 1차 문헌 및 지침 (PubMed 링크 포함), 13개 섹션으로 분류.
