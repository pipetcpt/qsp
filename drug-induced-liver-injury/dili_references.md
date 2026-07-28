# 약물유발 간손상 (DILI) — 참고문헌
# Drug-Induced Liver Injury — Annotated Reference List

이 목록은 `dili_qsp_model.dot` (기계론적 지도), `dili_mrgsolve_model.R` (ODE 모델),
`dili_shiny_app.R` (대시보드)에 사용된 구조·파라미터·임상 앵커의 근거입니다.
모델의 각 서브그래프 클러스터 및 각 ODE 블록에 대응하도록 섹션을 나누었습니다.

> **읽는 법.** ★ 표시는 모델 구조 자체(어떤 항이 존재하는가)를 결정한 문헌,
> ▲ 표시는 특정 파라미터 값 또는 임상 캘리브레이션 앵커를 제공한 문헌입니다.

---

## 1. 총론 · 역학 · 인과성 평가 (Overview, Epidemiology, Causality)

1. ★ Andrade RJ, Chalasani N, Björnsson ES, et al. **Drug-induced liver injury.** *Nat Rev Dis Primers.* 2019;5(1):58. — DILI 전체 병태생리의 표준 리뷰. 본 모델의 14개 클러스터 구획은 이 리뷰의 개념 틀을 따랐습니다. <https://pubmed.ncbi.nlm.nih.gov/31439850/>
2. ▲ Chalasani N, Bonkovsky HL, Fontana R, et al. **Features and Outcomes of 899 Patients With Drug-Induced Liver Injury: The DILIN Prospective Study.** *Gastroenterology.* 2015;148(7):1340-52.e7. — 만성화 비율(약 10–20%), 사망/이식률, 패턴 분포의 앵커. <https://pubmed.ncbi.nlm.nih.gov/25754159/>
3. ▲ Björnsson ES, Bergmann OM, Björnsson HK, et al. **Incidence, presentation, and outcomes in patients with drug-induced liver injury in the general population of Iceland.** *Gastroenterology.* 2013;144(7):1419-25. — 인구 기반 발생률 19.1/100,000/년. <https://pubmed.ncbi.nlm.nih.gov/23419359/>
4. ★ Danan G, Benichou C. **Causality assessment of adverse reactions to drugs—I. A novel method based on the conclusions of international consensus meetings: application to drug-induced liver injuries.** *J Clin Epidemiol.* 1993;46(11):1323-30. — RUCAM 원전. <https://pubmed.ncbi.nlm.nih.gov/8229110/>
5. Hayashi PH, Lucena MI, Fontana RJ, et al. **A revised electronic version of RUCAM for the diagnosis of DILI.** *Hepatology.* 2022;76(1):18-31. — RECAM. <https://pubmed.ncbi.nlm.nih.gov/35014066/>
6. ★ Aithal GP, Watkins PB, Andrade RJ, et al. **Case definition and phenotype standardization in drug-induced liver injury.** *Clin Pharmacol Ther.* 2011;89(6):806-15. — R-비(R ratio) 및 간세포형/담즙정체형/혼합형 정의의 표준. 모델의 `R` 출력 정의가 여기서 옵니다. <https://pubmed.ncbi.nlm.nih.gov/21544079/>
7. European Association for the Study of the Liver. **EASL Clinical Practice Guidelines: Drug-induced liver injury.** *J Hepatol.* 2019;70(6):1222-1261. <https://pubmed.ncbi.nlm.nih.gov/30926241/>
8. Fontana RJ, Liou I, Reuben A, et al. **AASLD practice guidance on drug, herbal, and dietary supplement-induced liver injury.** *Hepatology.* 2023;77(3):1036-1065. <https://pubmed.ncbi.nlm.nih.gov/35899384/>
9. ▲ Navarro VJ, Khan I, Björnsson E, et al. **Liver injury from herbal and dietary supplements.** *Hepatology.* 2017;65(1):363-373. — HDS가 미국 DILIN 등록의 20% 이상을 차지. <https://pubmed.ncbi.nlm.nih.gov/27677775/>

---

## 2. Hy's Law · 예후 · 규제 프레임워크 (Hy's Law, Prognosis, Regulatory)

10. ★▲ Temple R. **Hy's law: predicting serious hepatotoxicity.** *Pharmacoepidemiol Drug Saf.* 2006;15(4):241-3. — ALT ≥3×ULN **AND** TBIL ≥2×ULN → 사망/이식 약 10%. 본 모델의 핵심 주장(“Hy's Law는 속도×예비력의 결합 검정”)이 겨냥하는 규칙. <https://pubmed.ncbi.nlm.nih.gov/16552790/>
11. ★ Robles-Diaz M, Lucena MI, Kaplowitz N, et al. **Use of Hy's law and a new composite algorithm to predict acute liver failure in patients with drug-induced liver injury.** *Gastroenterology.* 2014;147(1):109-118.e5. — Hy's Law의 예측 성능과 nR 기반 개선. <https://pubmed.ncbi.nlm.nih.gov/24704526/>
12. ▲ U.S. FDA. **Guidance for Industry: Drug-Induced Liver Injury — Premarketing Clinical Evaluation.** 2009. — 규제상 Hy's Law 사례 정의 및 중단 규칙(stopping rule). <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/drug-induced-liver-injury-premarketing-clinical-evaluation>
13. ▲ O'Grady JG, Alexander GJ, Hayllar KM, Williams R. **Early indicators of prognosis in fulminant hepatic failure.** *Gastroenterology.* 1989;97(2):439-45. — King's College 기준 원전(APAP형: pH<7.30, 또는 INR>6.5 + Cr>3.4 mg/dL + 뇌증 III–IV). <https://pubmed.ncbi.nlm.nih.gov/2490426/>
14. Bernal W, Wendon J. **Acute liver failure.** *N Engl J Med.* 2013;369(26):2525-34. <https://pubmed.ncbi.nlm.nih.gov/24369077/>
15. ▲ Lee WM. **Acetaminophen (APAP) hepatotoxicity—Isn't it time for APAP to go away?** *J Hepatol.* 2017;67(6):1324-1331. — 미국 ALF의 최대 원인(약 46%). <https://pubmed.ncbi.nlm.nih.gov/28734939/>

---

## 3. 아세트아미노펜 생체변환 · 용량 임계 (APAP Biotransformation & Dose Thresholds)

16. ★▲ Mitchell JR, Jollow DJ, Potter WZ, Gillette JR, Brodie BB. **Acetaminophen-induced hepatic necrosis. IV. Protective role of glutathione.** *J Pharmacol Exp Ther.* 1973;187(1):211-7. — GSH 고갈이 괴사의 전제조건임을 확립한 고전. 모델의 “GSH 게이트” 구조(`v_bind`가 GSH에 의해 경쟁적으로 억제됨)의 근거. <https://pubmed.ncbi.nlm.nih.gov/4746329/>
17. ★ Jollow DJ, Mitchell JR, Potter WZ, et al. **Acetaminophen-induced hepatic necrosis. II. Role of covalent binding in vivo.** *J Pharmacol Exp Ther.* 1973;187(1):195-202. <https://pubmed.ncbi.nlm.nih.gov/4746327/>
18. ★▲ Dahlin DC, Miwa GT, Lu AY, Nelson SD. **N-acetyl-p-benzoquinone imine: a cytochrome P-450-mediated oxidation product of acetaminophen.** *Proc Natl Acad Sci USA.* 1984;81(5):1327-31. — NAPQI 동정. <https://pubmed.ncbi.nlm.nih.gov/6424115/>
19. ▲ Prescott LF. **Kinetics and metabolism of paracetamol and phenacetin.** *Br J Clin Pharmacol.* 1980;10(Suppl 2):291S-298S. — 치료용량에서 glucuronide ≈ 50–60%, sulfate ≈ 25–35%, CYP 경로 ≈ 5–10%, 신 배설 ≈ 5%. 모델의 Vmax/Km 배분과 총 청소율(약 19–21 L/h)의 앵커. <https://pubmed.ncbi.nlm.nih.gov/7002186/>
20. ★▲ Slattery JT, Wilson JM, Kalhorn TF, Nelson SD. **Dose-dependent pharmacokinetics of acetaminophen: evidence of glutathione depletion in humans.** *Clin Pharmacol Ther.* 1987;41(4):413-8. — 인체에서 용량 의존적 sulfation 포화. 모델의 **PAPS/UDPGA 보조인자 고갈 ODE**가 여기서 옵니다(이 항이 있어야 “고용량에서 활성화 분율이 스스로 상승”이 계산으로 나옵니다). <https://pubmed.ncbi.nlm.nih.gov/3829578/>
21. ▲ Rumack BH, Matthew H. **Acetaminophen poisoning and toxicity.** *Pediatrics.* 1975;55(6):871-6. — Rumack-Matthew 노모그램 원전(4시간 200 µg/mL 선, 미국은 150 µg/mL '치료선'). <https://pubmed.ncbi.nlm.nih.gov/1134886/>
22. ▲ Rumack BH. **Acetaminophen hepatotoxicity: the first 35 years.** *J Toxicol Clin Toxicol.* 2002;40(1):3-20. — 150 mg/kg / 7.5 g 임계, 250 mg/kg 이상에서 위험 급증. 모델 용량–반응 knee 위치의 임상 앵커. <https://pubmed.ncbi.nlm.nih.gov/11990202/>
23. ▲ Zhao L, Pickering G. **Paracetamol metabolism and related genetic differences.** *Drug Metab Rev.* 2011;43(1):41-52. <https://pubmed.ncbi.nlm.nih.gov/20977384/>
24. Chiew AL, Reith D, Pomerleau A, et al. **Updated guidelines for the management of paracetamol poisoning in Australia and New Zealand.** *Med J Aust.* 2020;212(4):175-183. <https://pubmed.ncbi.nlm.nih.gov/31786822/>

---

## 4. 글루타티온 · 시스테인 · 티올 방어 (Glutathione, Cysteine, Thiol Defence)

25. ★ Lu SC. **Glutathione synthesis.** *Biochim Biophys Acta.* 2013;1830(5):3143-53. — GCL이 율속 효소이고 **시스테인이 율속 기질**임. 모델에서 NAC 효과가 `CYS` 상태변수를 통해서만 들어가는 이유. <https://pubmed.ncbi.nlm.nih.gov/22995213/>
26. ★▲ Lu SC. **Regulation of hepatic glutathione synthesis: current concepts and controversies.** *FASEB J.* 1999;13(10):1169-83. — 간 GSH 농도 5–10 mM, 회전 반감기 2–4시간. <https://pubmed.ncbi.nlm.nih.gov/10385608/>
27. ★ Fernández-Checa JC, Kaplowitz N. **Hepatic mitochondrial glutathione: transport and role in disease and toxicity.** *Toxicol Appl Pharmacol.* 2005;204(3):263-73. — 미토콘드리아 GSH 풀이 별도 구획이며 그 고갈이 치명적. <https://pubmed.ncbi.nlm.nih.gov/15845418/>
28. ▲ Chen Y, Dong H, Thompson DC, et al. **Glutathione defense mechanism in liver injury: insights from animal models.** *Food Chem Toxicol.* 2013;60:38-44. <https://pubmed.ncbi.nlm.nih.gov/23856494/>
29. ★ Yuan L, Kaplowitz N. **Glutathione in liver diseases and hepatotoxicity.** *Mol Aspects Med.* 2009;30(1-2):29-41. <https://pubmed.ncbi.nlm.nih.gov/18786561/>

---

## 5. 미토콘드리아 독성 · 생체에너지학 (Mitochondrial Toxicity & Bioenergetics)

30. ★ Jaeschke H, McGill MR, Ramachandran A. **Oxidant stress, mitochondria, and cell death mechanisms in drug-induced liver injury: lessons learned from acetaminophen hepatotoxicity.** *Drug Metab Rev.* 2012;44(1):88-106. — 미토콘드리아 부가체 → ROS → MPT → 괴사 경로의 표준 서술. 모델 클러스터 4의 뼈대. <https://pubmed.ncbi.nlm.nih.gov/22229890/>
31. ★▲ Kon K, Kim JS, Jaeschke H, Lemasters JJ. **Mitochondrial permeability transition in acetaminophen-induced necrosis and apoptosis of cultured mouse hepatocytes.** *Hepatology.* 2004;40(5):1170-9. — **ATP가 있으면 아폽토시스, 없으면 괴사**. 모델의 MPT 항에 있는 ATP 게이트 `1/(1+(ATP/K)^4)`의 근거. <https://pubmed.ncbi.nlm.nih.gov/15486922/>
32. ★ Ramachandran A, Jaeschke H. **Acetaminophen Hepatotoxicity.** *Semin Liver Dis.* 2019;39(2):221-234. <https://pubmed.ncbi.nlm.nih.gov/30849782/>
33. ★ Pessayre D, Fromenty B, Berson A, et al. **Central role of mitochondria in drug-induced liver injury.** *Drug Metab Rev.* 2012;44(1):34-87. — 탈공역·β-산화 억제·mtDNA 고갈의 세 가지 독립적 미토콘드리아 손상 축. 모델의 `FUNC`(탈공역/FAO 차단 부담) 파라미터. <https://pubmed.ncbi.nlm.nih.gov/21892896/>
34. ▲ Nadanaciva S, Will Y. **New insights in drug-induced mitochondrial toxicity.** *Curr Pharm Des.* 2011;17(20):2100-12. <https://pubmed.ncbi.nlm.nih.gov/21718246/>
35. ★ Masubuchi Y, Suda C, Horie T. **Involvement of mitochondrial permeability transition in acetaminophen-induced liver injury in mice.** *J Hepatol.* 2005;42(1):110-6. <https://pubmed.ncbi.nlm.nih.gov/15629515/>

---

## 6. JNK–Sab 증폭 고리 (JNK–Sab Amplification Loop) — 쌍안정성의 근원

36. ★★ Hanawa N, Shinohara M, Saberi B, et al. **Role of JNK translocation to mitochondria leading to inhibition of mitochondria bioenergetics in acetaminophen-induced liver injury.** *J Biol Chem.* 2008;283(20):13565-77. — **p-JNK가 미토콘드리아로 전위하여 호흡을 억제하고 ROS를 더 만든다** = 본 모델의 정귀환 고리 그 자체. <https://pubmed.ncbi.nlm.nih.gov/18337250/>
37. ★★ Win S, Than TA, Han D, Petrovic LM, Kaplowitz N. **c-Jun N-terminal kinase (JNK)-dependent acute liver injury from acetaminophen or tumor necrosis factor (TNF) requires mitochondrial Sab protein expression in mice.** *J Biol Chem.* 2011;286(40):35071-8. — Sab(SH3BP5)이 고리의 도킹 단백질임을 증명. Sab 결손이면 손상이 사라짐 → 고리의 **필수성**. <https://pubmed.ncbi.nlm.nih.gov/21844199/>
38. ★★ Win S, Than TA, Zhang J, Oo C, Min RWM, Kaplowitz N. **New insights into the role and mechanism of c-Jun-N-terminal kinase signaling in the pathobiology of liver diseases.** *Hepatology.* 2018;67(5):2013-2024. — 일시적 JNK 활성화는 생존, **지속적** 활성화는 사망 신호. 모델이 “쌍안정”으로 표현하는 임상적 실체. <https://pubmed.ncbi.nlm.nih.gov/29194686/>
39. ★ Gunawan BK, Liu ZX, Han D, Hanawa N, Gaarde WA, Kaplowitz N. **c-Jun N-terminal kinase plays a major role in murine acetaminophen hepatotoxicity.** *Gastroenterology.* 2006;131(1):165-78. <https://pubmed.ncbi.nlm.nih.gov/16831600/>
40. ★ Nakagawa H, Maeda S, Hikiba Y, et al. **Deletion of apoptosis signal-regulating kinase 1 attenuates acetaminophen-induced liver injury by inhibiting c-Jun N-terminal kinase activation.** *Gastroenterology.* 2008;135(4):1311-21. — ASK1이 ROS→JNK 전달자. <https://pubmed.ncbi.nlm.nih.gov/18700144/>
41. Saito C, Lemasters JJ, Jaeschke H. **c-Jun N-terminal kinase modulates oxidant stress and peroxynitrite formation independent of inducible nitric oxide synthase in acetaminophen hepatotoxicity.** *Toxicol Appl Pharmacol.* 2010;246(1-2):8-17. <https://pubmed.ncbi.nlm.nih.gov/20423716/>
42. ★ Win S, Min RWM, Zhang J, et al. **Hepatic Mitochondrial SAB Deletion or Knockdown Alleviates Diet-Induced Metabolic Syndrome, Steatohepatitis, and Hepatic Fibrosis.** *Hepatology.* 2021;74(6):3127-3145. — 고리가 APAP 외 다른 간질환에도 공통. <https://pubmed.ncbi.nlm.nih.gov/34272738/>

---

## 7. Nrf2 적응 반응 · 오토파지 (Nrf2 Adaptation & Autophagy)

43. ★ Enomoto A, Itoh K, Nagayoshi E, et al. **High sensitivity of Nrf2 knockout mice to acetaminophen hepatotoxicity associated with decreased expression of ARE-regulated drug metabolizing enzymes and antioxidant genes.** *Toxicol Sci.* 2001;59(1):169-77. — Nrf2 결손 시 감수성 급증 → 모델에서 `NRF2`가 GSH 합성능과 ROS 제거능을 동시에 곱셈으로 올리는 근거. <https://pubmed.ncbi.nlm.nih.gov/11134556/>
44. ★ Goldring CE, Kitteringham NR, Elsby R, et al. **Activation of hepatic Nrf2 in vivo by acetaminophen in CD-1 mice.** *Hepatology.* 2004;39(5):1267-76. <https://pubmed.ncbi.nlm.nih.gov/15122754/>
45. ★ Ni HM, Bockus A, Boggess N, Jaeschke H, Ding WX. **Activation of autophagy protects against acetaminophen-induced hepatotoxicity.** *Hepatology.* 2012;55(1):222-32. — 오토파지가 손상 미토콘드리아·부가체를 제거. 모델의 `KADD_REP × autophagy(ATP)` 항. <https://pubmed.ncbi.nlm.nih.gov/21932416/>
46. ▲ Dinkova-Kostova AT, Kostov RV, Kazantsev AG. **The role of Nrf2 signaling in counteracting neurodegenerative diseases.** *FEBS J.* 2018;285(19):3576-3590. — KEAP1 시스테인 센서 기전 일반론. <https://pubmed.ncbi.nlm.nih.gov/29323772/>
47. ★ Kaplowitz N. **Adaptation vs. injury: what determines the outcome?** In: *Drug-Induced Liver Disease* (3rd ed). — “적응(adaptation)”이 대다수 노출자의 귀결이라는 개념. (개관: Watkins PB. **Idiosyncratic liver injury: challenges and approaches.** *Toxicol Pathol.* 2005;33(1):1-5. <https://pubmed.ncbi.nlm.nih.gov/15805049/>)

---

## 8. 담즙산 · BSEP · 담즙정체성 DILI (Bile Acids, BSEP, Cholestatic DILI)

48. ★★ Morgan RE, Trauner M, van Staden CJ, et al. **Interference with bile salt export pump function is a susceptibility factor for human liver injury in drug development.** *Toxicol Sci.* 2010;118(2):485-500. — BSEP 억제 IC50이 낮을수록 DILI 위험이 높음. 모델의 `KI_BSEP` 경쟁적 억제 항. <https://pubmed.ncbi.nlm.nih.gov/20829430/>
49. ★▲ Dawson S, Stahl S, Paul N, Barber J, Kenna JG. **In vitro inhibition of the bile salt export pump correlates with risk of cholestatic drug-induced liver injury in humans.** *Drug Metab Dispos.* 2012;40(1):130-8. — BSEP IC50 < 25–50 µM 인 약물에서 담즙정체형 DILI 위험 상승. <https://pubmed.ncbi.nlm.nih.gov/21965623/>
50. ★ Stieger B, Fattinger K, Madon J, Kullak-Ublick GA, Meier PJ. **Drug- and estrogen-induced cholestasis through inhibition of the hepatocellular bile salt export pump (Bsep) of rat liver.** *Gastroenterology.* 2000;118(2):422-30. <https://pubmed.ncbi.nlm.nih.gov/10648470/>
51. ★ Woolbright BL, Jaeschke H. **Novel insight into mechanisms of cholestatic liver injury.** *World J Gastroenterol.* 2012;18(36):4985-93. — 담즙산의 세포독성이 단순 세정 작용이 아니라 염증 매개(CXCL1/2 유도)임. 모델에서 `BADETER → CXCL` 간선. <https://pubmed.ncbi.nlm.nih.gov/23049205/>
52. ★ Chiang JYL. **Bile acid metabolism and signaling.** *Compr Physiol.* 2013;3(3):1191-212. — FXR-SHP-CYP7A1 음성 되먹임 및 장 FGF19 축. 모델의 `FMAX_FXR` 부분 억제. <https://pubmed.ncbi.nlm.nih.gov/23897684/>
53. ▲ Slijepcevic D, van de Graaf SFJ. **Bile Acid Uptake Transporters as Targets for Therapy.** *Dig Dis.* 2017;35(3):251-258. <https://pubmed.ncbi.nlm.nih.gov/28249287/>
54. ★ Aleo MD, Luo Y, Swiss R, Bonin PD, Potter DM, Will Y. **Human drug-induced liver injury severity is highly associated with dual inhibition of liver mitochondrial function and bile salt export pump.** *Hepatology.* 2014;60(3):1015-22. — **미토 독성 + BSEP 억제의 이중 타격**이 최악의 조합. 모델에서 BSEP 활성이 `× ATP` 로 스케일되는 이유(간세포 손상이 2차 담즙정체를 만든다). <https://pubmed.ncbi.nlm.nih.gov/24799086/>
55. Padda MS, Sanchez M, Akhtar AJ, Boyer JL. **Drug-induced cholestasis.** *Hepatology.* 2011;53(4):1377-87. <https://pubmed.ncbi.nlm.nih.gov/21480339/>

---

## 9. 선천면역 · 무균성 염증 (Innate Immunity & Sterile Inflammation)

56. ★ Jaeschke H, Ramachandran A. **Mechanisms and pathophysiological significance of sterile inflammation during acetaminophen hepatotoxicity.** *Food Chem Toxicol.* 2020;138:111240. — DAMP → 쿠퍼세포 → 사이토카인 축이 손상을 **확대하기도 하고 수복을 개시하기도 함**(TNF의 양날). <https://pubmed.ncbi.nlm.nih.gov/32145352/>
57. ★ Scaffidi P, Misteli T, Bianchi ME. **Release of chromatin protein HMGB1 by necrotic cells triggers inflammation.** *Nature.* 2002;418(6894):191-5. — HMGB1 = 대표적 DAMP. <https://pubmed.ncbi.nlm.nih.gov/12110890/>
58. ★ Imaeda AB, Watanabe A, Sohail MA, et al. **Acetaminophen-induced hepatotoxicity in mice is dependent on Tlr9 and the Nalp3 inflammasome.** *J Clin Invest.* 2009;119(2):305-14. <https://pubmed.ncbi.nlm.nih.gov/19164858/>
59. ★ Bourdi M, Masubuchi Y, Reilly TP, et al. **Protection against acetaminophen-induced liver injury and lethality by interleukin 10: role of inducible nitric oxide synthase.** *Hepatology.* 2002;35(2):289-98. — IL-10의 보호 효과. 모델의 `IL10 ⊣ KC, TNF` 억제 간선. <https://pubmed.ncbi.nlm.nih.gov/11826401/>
60. ★ Yamada Y, Kirillova I, Peschon JJ, Fausto N. **Initiation of liver growth by tumor necrosis factor: deficient liver regeneration in mice lacking type I tumor necrosis factor receptor.** *Proc Natl Acad Sci USA.* 1997;94(4):1441-6. — TNF/IL-6가 재생 프라이밍. 모델의 `prim` 항. <https://pubmed.ncbi.nlm.nih.gov/9037072/>
61. Antoniades CG, Quaglia A, Taams LS, et al. **Source and characterization of hepatic macrophages in acetaminophen-induced acute liver failure in humans.** *Hepatology.* 2012;56(2):735-46. <https://pubmed.ncbi.nlm.nih.gov/22334567/>
62. Woolbright BL, Jaeschke H. **Role of the inflammasome in acetaminophen-induced liver injury and acute liver failure.** *J Hepatol.* 2017;66(4):836-848. <https://pubmed.ncbi.nlm.nih.gov/27913221/>

---

## 10. 적응면역 · HLA · 특이체질 DILI (Adaptive Immunity, HLA, Idiosyncratic DILI)

63. ★★ Uetrecht J, Naisbitt DJ. **Idiosyncratic adverse drug reactions: current concepts.** *Pharmacol Rev.* 2013;65(2):779-808. — 합텐 가설, p-i 개념, 위험 신호(danger signal) 요구성. 모델의 T세포 활성화 항이 **부가체 × 위험신호 × (1/관용)** 의 곱인 이유. <https://pubmed.ncbi.nlm.nih.gov/23476052/>
64. ★★ Daly AK, Donaldson PT, Bhatnagar P, et al. **HLA-B*5701 genotype is a major determinant of drug-induced liver injury due to flucloxacillin.** *Nat Genet.* 2009;41(7):816-9. — 오즈비 약 80. 모델의 `HLA` 스위치. <https://pubmed.ncbi.nlm.nih.gov/19483685/>
65. ★▲ Lucena MI, Molokhia M, Shen Y, et al. **Susceptibility to amoxicillin-clavulanate-induced liver injury is influenced by multiple HLA class I and II alleles.** *Gastroenterology.* 2011;141(1):338-47. — HLA-DRB1*15:01 등. <https://pubmed.ncbi.nlm.nih.gov/21570397/>
66. ★ Urban TJ, Nicoletti P, Chalasani N, et al. **Minocycline hepatotoxicity: Clinical characterization and identification of HLA-B*35:02 as a risk factor.** *J Hepatol.* 2017;67(1):137-144. <https://pubmed.ncbi.nlm.nih.gov/28323125/>
67. ★★ Metushi IG, Hayes MA, Uetrecht J. **Treatment of PD-1−/− mice with anti-CTLA-4 antibody unmasks liver injury from amodiaquine.** *Hepatology.* 2015;61(4):1332-42. — **면역 관용을 제거하면 잠재적 부가체 부하가 임상 간염으로 전환된다**는 직접 증명. 모델의 ICI 시나리오 구조가 여기서 옵니다. <https://pubmed.ncbi.nlm.nih.gov/25482010/>
68. ★ Chakraborty M, Fullerton AM, Semple K, et al. **Drug-induced allergic hepatitis develops in mice when myeloid-derived suppressor cells are depleted prior to halothane treatment.** *Hepatology.* 2015;62(2):546-57. <https://pubmed.ncbi.nlm.nih.gov/25712247/>
68b. Cho T, Uetrecht J. **How Reactive Metabolites Induce an Immune Response That Sometimes Leads to an Idiosyncratic Drug Reaction.** *Chem Res Toxicol.* 2017;30(1):295-314. <https://pubmed.ncbi.nlm.nih.gov/27775332/>
69. ▲ De Martin E, Michot JM, Papouin B, et al. **Characterization of liver injury induced by cancer immunotherapy using immune checkpoint inhibitors.** *J Hepatol.* 2018;68(6):1181-1190. — ICI 간염의 임상·병리 표현형과 스테로이드 반응. <https://pubmed.ncbi.nlm.nih.gov/29427729/>
70. ▲ Peeraphatdit TB, Wang J, Odenwald MA, Hu S, Hart J, Charlton MR. **Hepatotoxicity From Immune Checkpoint Inhibitors: A Systematic Review and Management Recommendation.** *Hepatology.* 2020;72(1):315-329. — grade 3–4 간독성 발생률 및 스테로이드/MMF 단계적 치료. <https://pubmed.ncbi.nlm.nih.gov/32167613/>

---

## 11. 바이오마커 (Mechanistic Biomarkers)

71. ★★ Antoine DJ, Dear JW, Lewis PS, et al. **Mechanistic biomarkers provide early and sensitive detection of acetaminophen-induced acute liver injury at first presentation to hospital.** *Hepatology.* 2013;58(2):777-87. — miR-122·HMGB1·K18이 ALT보다 먼저 상승. 모델에서 miR-122의 짧은 반감기 때문에 **선행 상승이 계산으로 나옵니다**. <https://pubmed.ncbi.nlm.nih.gov/23390034/>
72. ★▲ Starkey Lewis PJ, Dear J, Platt V, et al. **Circulating microRNAs as potential markers of human drug-induced liver injury.** *Hepatology.* 2011;54(5):1767-76. <https://pubmed.ncbi.nlm.nih.gov/22045675/>
73. ▲ Church RJ, Kullak-Ublick GA, Aubrecht J, et al. **Candidate biomarkers for the diagnosis and prognosis of drug-induced liver injury: An international collaborative effort.** *Hepatology.* 2019;69(2):760-773. — 다기관 검증. <https://pubmed.ncbi.nlm.nih.gov/29357190/>
74. ★▲ McGill MR, Sharpe MR, Williams CD, Taha M, Curry SC, Jaeschke H. **The mechanism underlying acetaminophen-induced hepatotoxicity in humans and mice involves mitochondrial damage and nuclear DNA fragmentation.** *J Clin Invest.* 2012;122(4):1574-83. — 인체에서 GLDH·mtDNA·nDNA 단편의 상승. <https://pubmed.ncbi.nlm.nih.gov/22378043/>
75. ▲ Dear JW, Clarke JI, Francis B, et al. **Risk stratification after paracetamol overdose using mechanistic biomarkers: results from two prospective cohort studies.** *Lancet Gastroenterol Hepatol.* 2018;3(2):104-113. <https://pubmed.ncbi.nlm.nih.gov/29146439/>
76. ▲ James LP, Letzig L, Simpson PM, et al. **Pharmacokinetics of acetaminophen-protein adducts in adults with acetaminophen overdose and acute liver failure.** *Drug Metab Dispos.* 2009;37(8):1779-84. — 혈중 APAP-Cys 부가체(노출 확진 바이오마커)의 동태. <https://pubmed.ncbi.nlm.nih.gov/19439490/>
77. ▲ Ozer J, Ratner M, Shaw M, Bailey W, Schomaker S. **The current state of serum biomarkers of hepatotoxicity.** *Toxicology.* 2008;245(3):194-205. — ALT 혈청 반감기 약 47시간, AST 약 17시간. 모델의 `KALT_EL`·`KAST_EL`. <https://pubmed.ncbi.nlm.nih.gov/18291570/>

---

## 12. N-아세틸시스테인 및 기타 치료 (NAC and Other Therapeutics)

78. ★★ Prescott LF, Illingworth RN, Critchley JA, Stewart MJ, Adam RD, Proudfoot AT. **Intravenous N-acetylcystine: the treatment of choice for paracetamol poisoning.** *Br Med J.* 1979;2(6198):1097-100. — Prescott 요법(150 → 50 → 100 mg/kg) 원전. 모델의 NAC 주입 스케줄. <https://pubmed.ncbi.nlm.nih.gov/519312/>
79. ★★ Smilkstein MJ, Knapp GL, Kulig KW, Rumack BH. **Efficacy of oral N-acetylcysteine in the treatment of acetaminophen overdose. Analysis of the national multicenter study (1976 to 1985).** *N Engl J Med.* 1988;319(24):1557-62. — **10시간 이내 투여 시 간독성 거의 0%, 16–24시간에서 급격히 악화.** 본 모델의 NAC 시간창 시뮬레이션이 재현해야 하는 핵심 임상 앵커. <https://pubmed.ncbi.nlm.nih.gov/3059186/>
80. ▲ Bateman DN, Dear JW, Thanacoody HK, et al. **Reduction of adverse effects from intravenous acetylcysteine treatment for paracetamol poisoning: a randomised controlled trial.** *Lancet.* 2014;383(9918):697-704. — SNAP 12시간 요법. <https://pubmed.ncbi.nlm.nih.gov/24290406/>
81. ★ Lauterburg BH, Corcoran GB, Mitchell JR. **Mechanism of action of N-acetylcysteine in the protection against the hepatotoxicity of acetaminophen in rats in vivo.** *J Clin Invest.* 1983;71(4):980-91. — NAC의 주 기전은 **시스테인 공급을 통한 GSH 재합성**. <https://pubmed.ncbi.nlm.nih.gov/6833497/>
82. ★ Lee WM, Hynan LS, Rossaro L, et al. **Intravenous N-acetylcysteine improves transplant-free survival in early stage non-acetaminophen acute liver failure.** *Gastroenterology.* 2009;137(3):856-64. — 비-APAP ALF에서도 유익(미세순환/산소전달 기전). <https://pubmed.ncbi.nlm.nih.gov/19524577/>
83. ▲ Akakpo JY, Ramachandran A, Kandel SE, et al. **4-Methylpyrazole protects against acetaminophen hepatotoxicity in mice and in primary human hepatocytes.** *Hum Exp Toxicol.* 2018;37(12):1310-1322. — 포메피졸(CYP2E1 + JNK 억제). <https://pubmed.ncbi.nlm.nih.gov/29768939/>
84. ▲ Morrison EE, Oatey K, Gallagher B, et al. **Principal results of a randomised open label exploratory, safety and tolerability study with calmangafodipir in patients treated with a 12 h regimen of N-acetylcysteine for paracetamol overdose (POP trial).** *EBioMedicine.* 2019;46:423-430. <https://pubmed.ncbi.nlm.nih.gov/31351929/>
85. ▲ Larsen FS, Schmidt LE, Bernsmeier C, et al. **High-volume plasma exchange in patients with acute liver failure: An open randomised controlled trial.** *J Hepatol.* 2016;64(1):69-78. — 고용량 혈장교환이 ALF 생존 개선. <https://pubmed.ncbi.nlm.nih.gov/26325537/>
86. ▲ Lheureux PE, Penaloza A, Zahir S, Gris M. **Science review: carnitine in the treatment of valproic acid-induced toxicity.** *Crit Care.* 2005;9(5):431-40. <https://pubmed.ncbi.nlm.nih.gov/16277730/>

---

## 13. 숙주 위험 인자 · 약물 물성 (Host Risk Factors & Physicochemical Rules)

87. ★▲ Chen M, Borlak J, Tong W. **High lipophilicity and high daily dose of oral medications are associated with significant risk for drug-induced liver injury.** *Hepatology.* 2013;58(1):388-96. — “Rule-of-Two”(logP > 3 그리고 일일용량 > 100 mg). 모델 지도의 `LIPOPHIL` 노드. <https://pubmed.ncbi.nlm.nih.gov/23258593/>
88. ★▲ Lammert C, Einarsson S, Saha C, Niklasson A, Bjornsson E, Chalasani N. **Relationship between daily dose of oral medications and idiosyncratic drug-induced liver injury: search for signals.** *Hepatology.* 2008;47(6):2003-9. — 일일용량 ≥50 mg 에서 위험 집중. <https://pubmed.ncbi.nlm.nih.gov/18454504/>
89. ★ Whitcomb DC, Block GD. **Association of acetaminophen hepatotoxicity with fasting and ethanol use.** *JAMA.* 1994;272(23):1845-50. — 금식·음주의 이중 타격. 모델의 취약 숙주 시나리오(`FCYP`↑, `FGSH`↓, `CYSBASE`↓). <https://pubmed.ncbi.nlm.nih.gov/7990219/>
90. ★ Zimmerman HJ, Maddrey WC. **Acetaminophen (paracetamol) hepatotoxicity with regular intake of alcohol: analysis of instances of therapeutic misadventure.** *Hepatology.* 1995;22(3):767-73. <https://pubmed.ncbi.nlm.nih.gov/7657281/>
91. ▲ Michaut A, Moreau C, Robin MA, Fromenty B. **Acetaminophen-induced liver injury in obesity and nonalcoholic fatty liver disease.** *Liver Int.* 2014;34(7):e171-9. — 비만/MASLD의 CYP2E1 상승과 미토 예비력 감소. <https://pubmed.ncbi.nlm.nih.gov/24575897/>
92. ▲ Nicoletti P, Aithal GP, Bjornsson ES, et al. **Association of Liver Injury From Specific Drugs, or Groups of Drugs, With Polymorphisms in HLA and Other Genes in a Genome-Wide Association Study.** *Gastroenterology.* 2017;152(5):1078-1089. — PTPN22 등 비-HLA 인자. <https://pubmed.ncbi.nlm.nih.gov/28043905/>

---

## 14. 재생 · 예비력 · 만성화 (Regeneration, Functional Reserve, Chronicity)

93. ★ Michalopoulos GK, Bhushan B. **Liver regeneration: biological and pathological mechanisms and implications.** *Nat Rev Gastroenterol Hepatol.* 2021;18(1):40-55. — HGF/c-Met 주도 재생, TGF-β1의 종료 신호. 모델의 `regen` 항 구조. <https://pubmed.ncbi.nlm.nih.gov/32764740/>
94. ★ Fausto N, Campbell JS, Riehle KJ. **Liver regeneration.** *Hepatology.* 2006;43(2 Suppl 1):S45-53. <https://pubmed.ncbi.nlm.nih.gov/16447274/>
95. ▲ Fontana RJ, Hayashi PH, Barnhart H, et al. **Persistent liver biochemistry abnormalities are more common in older patients and those with cholestatic drug induced liver injury.** *Am J Gastroenterol.* 2015;110(10):1450-9. — 담즙정체형·고령에서 만성화가 흔함. <https://pubmed.ncbi.nlm.nih.gov/26346867/>
96. ▲ Bonkovsky HL, Kleiner DE, Gu J, et al. **Clinical presentations and outcomes of bile duct loss caused by drugs and herbal and dietary supplements.** *Hepatology.* 2017;65(4):1267-1277. — 소실성 담관 증후군(VBDS). <https://pubmed.ncbi.nlm.nih.gov/27981596/>

---

## 15. QSP 모델링 방법론 · DILIsym (QSP Methodology & Prior Models)

97. ★★ Howell BA, Yang Y, Kumar R, et al. **In vitro to in vivo extrapolation and species response comparisons for drug-induced liver injury (DILI) using DILIsym™: a mechanistic, mathematical model of DILI.** *J Pharmacokinet Pharmacodyn.* 2012;39(5):527-41. — DILI QSP의 표준 플랫폼. 본 모델은 DILIsym의 공개된 구조 개념(GSH·미토·담즙산·생명주기 서브모델)을 참조하되 독립적으로 단순화·재구현한 것입니다. <https://pubmed.ncbi.nlm.nih.gov/22996471/>
98. ★★ Watkins PB. **DILIsym: Quantitative systems toxicology impacting drug development.** *Curr Opin Toxicol.* 2020;23-24:67-73. <https://doi.org/10.1016/j.cotox.2020.06.003>
99. ★ Woodhead JL, Howell BA, Yang Y, et al. **An analysis of N-acetylcysteine treatment for acetaminophen overdose using a systems model of drug-induced liver injury.** *J Pharmacol Exp Ther.* 2012;342(2):529-40. — NAC 시간창을 QSP로 설명한 선행 연구. 본 모델 [4]번 분석의 직접적 선례. <https://pubmed.ncbi.nlm.nih.gov/22645248/>
100. ★ Woodhead JL, Yang K, Brouwer KLR, et al. **Mechanistic modeling reveals the critical knowledge gaps in bile acid-mediated DILI.** *CPT Pharmacometrics Syst Pharmacol.* 2014;3(3):e123. <https://pubmed.ncbi.nlm.nih.gov/24646538/>
101. ★ Longo DM, Yang Y, Watkins PB, Howell BA, Siler SQ. **Elucidating Differences in the Hepatotoxic Potential of Tolcapone and Entacapone With DILIsym, a Mechanistic Model of Drug-Induced Liver Injury.** *CPT Pharmacometrics Syst Pharmacol.* 2016;5(1):31-9. <https://pubmed.ncbi.nlm.nih.gov/26844013/>
102. ★ Shoda LKM, Woodhead JL, Siler SQ, Watkins PB, Howell BA. **Linking physiology to toxicity using DILIsym®, a mechanistic mathematical model of drug-induced liver injury.** *Biopharm Drug Dispos.* 2014;35(1):33-49. <https://pubmed.ncbi.nlm.nih.gov/24214486/>
103. ★ Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn.* 2015;42:S84-85. — mrgsolve. <https://github.com/metrumresearchgroup/mrgsolve>
104. Bhattacharya S, Shoda LKM, Zhang Q, et al. **Modeling drug- and chemical-induced hepatotoxicity with systems biology approaches.** *Front Physiol.* 2012;3:462. <https://pubmed.ncbi.nlm.nih.gov/23248599/>

---

## 16. 쌍안정성 · 임계 현상 (Bistability & Threshold Behaviour) — 모델 구조의 이론적 근거

105. ★ Ferrell JE Jr. **Self-perpetuating states in signal transduction: positive feedback, double-negative feedback and bistability.** *Curr Opin Cell Biol.* 2002;14(2):140-8. — 포화형 이득을 갖는 정귀환 고리가 쌍안정을 만드는 조건. 모델의 JNK–Sab 고리 해석 틀. <https://pubmed.ncbi.nlm.nih.gov/11891111/>
106. ★ Tyson JJ, Chen KC, Novak B. **Sniffers, buzzers, toggles and blinkers: dynamics of regulatory and signaling pathways in the cell.** *Curr Opin Cell Biol.* 2003;15(2):221-31. <https://pubmed.ncbi.nlm.nih.gov/12648679/>
107. ★ Bagci EZ, Vodovotz Y, Billiar TR, Ermentrout GB, Bahar I. **Bistability in apoptosis: roles of bax, bcl-2, and mitochondrial permeability transition pores.** *Biophys J.* 2006;90(5):1546-59. — MPT를 포함한 세포사 결정의 쌍안정 모델 선례. <https://pubmed.ncbi.nlm.nih.gov/16339882/>

---

### 인용 관례 (Citation Conventions)
- PubMed 링크는 `https://pubmed.ncbi.nlm.nih.gov/<PMID>/` 형식입니다.
- ★/▲ 표시의 의미는 문서 상단 참조.
- 총 문헌 수: 108건 (★ 구조 결정 문헌 다수 포함).
- 본 문헌 목록은 교육·연구 목적의 모델 문서화를 위한 것이며, 임상 진료 지침을
  대체하지 않습니다.
