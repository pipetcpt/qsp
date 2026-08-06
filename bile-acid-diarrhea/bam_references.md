# 담즙산 설사 / 담즙산 흡수장애 (Bile Acid Diarrhoea / Malabsorption) — 참고문헌

> **모든 PMID는 NCBI E-utilities로 조회하여 검증한 실제 문헌입니다.**
> 링크 형식: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`
> 각 섹션 머리말은 그 문헌군이 이 QSP 모델의 **어느 방정식·어느 파라미터**를
> 뒷받침하는지를 밝힙니다. 모델의 주장은 문헌에서 나온 것과 모델이 스스로
> 유도한 것을 구분해서 표시했습니다 (⚙ = 모델이 유도한 것, 📖 = 문헌이 준 것).

---

## 목차

| 섹션 | 주제 | 편수 |
|---|---|---:|
| A | 장간순환의 정량 생리 — 풀 크기·순환 횟수·합성속도 | 14 |
| B | ASBT / 회장 수송 — `φ` 축의 분자적 실체 | 12 |
| C | FXR–FGF19–CYP7A1 되먹임 — `κ` 축과 되먹임 이득 | 13 |
| D | 담즙산 설사의 병태생리와 아형 | 12 |
| E | 대장 분비·운동 — 3 mM 역치와 양성 되먹임 | 11 |
| F | 진단검사 — SeHCAT · C4 · FGF19 · 대변 담즙산 | 14 |
| G | 담즙산 결합수지 (class A) | 8 |
| H | FXR 작용제 (class B) | 9 |
| I | ASBT 억제제 (class C) — 거울상 | 8 |
| J | 회장절제·100 cm 규칙·지방변 | 7 |
| K | 미생물총 — 탈포합과 7α-탈수산화 (class E) | 8 |
| L | IBS-D 중첩과 대증치료 (class D) | 8 |
| M | 담즙산 QSP·PBPK 모델링과 도구 | 10 |
| | **합계** | **134** |

---

## A. 장간순환의 정량 생리 — 풀 크기·순환 횟수·합성속도

> 모델의 정상 정상상태 보정값이 여기서 나옵니다: 총 풀 3 g(6,000 µmol),
> 하루 4–6 회 순환, 십이지장 전달 20–40 mmol/day, 간 합성 0.35 g/day.
> 특히 Hofmann 1983은 이 저장소 모델의 직계 조상에 해당하는 최초의
> 장간순환 시뮬레이션 모델입니다.

1. Hofmann AF, et al. *Description and simulation of a physiological pharmacokinetic model for the metabolism and enterohepatic circulation of bile acids in man.* J Clin Invest 1983. [PMID 6682120](https://pubmed.ncbi.nlm.nih.gov/6682120/)
2. Chiang JY. *Bile acid metabolism and signaling.* Compr Physiol 2013. [PMID 23897684](https://pubmed.ncbi.nlm.nih.gov/23897684/)
3. Chiang JYL. *Bile Acid Metabolism in Liver Pathobiology.* Gene Expr 2018. [PMID 29325602](https://pubmed.ncbi.nlm.nih.gov/29325602/)
4. Chiang JY. *Bile acids: regulation of synthesis.* J Lipid Res 2009. [PMID 19346330](https://pubmed.ncbi.nlm.nih.gov/19346330/)
5. Di Ciaula A, et al. *Bile Acid Physiology.* Ann Hepatol 2017. [PMID 29080336](https://pubmed.ncbi.nlm.nih.gov/29080336/)
6. Roberts MS, et al. *Enterohepatic circulation: physiological, pharmacokinetic and clinical implications.* Clin Pharmacokinet 2002. [PMID 12162761](https://pubmed.ncbi.nlm.nih.gov/12162761/)
7. Vantrappen G, et al. *A new method for the measurement of bile acid turnover and pool size by a double label, single intubation technique.* J Lipid Res 1981. [PMID 7017051](https://pubmed.ncbi.nlm.nih.gov/7017051/)
8. Hulzebos CV, et al. *Measurement of parameters of cholic acid kinetics in plasma using a microscale stable isotope dilution technique.* J Lipid Res 2001. [PMID 11714862](https://pubmed.ncbi.nlm.nih.gov/11714862/)
9. Stellaard F, et al. *Determination of deoxycholic acid pool size and input rate using [24-13C]deoxycholic acid.* J Lipid Res 1986. [PMID 3559388](https://pubmed.ncbi.nlm.nih.gov/3559388/)
10. Schwartz CC, et al. *Cholesterol kinetics in subjects with bile fistula. Positive relationship between size of the bile acid precursor pool and bile acid synthetic rate.* J Clin Invest 1993. [PMID 8450070](https://pubmed.ncbi.nlm.nih.gov/8450070/)
11. Balistreri WF, et al. *Validation of use of 11,12-2H-labeled chenodeoxycholic acid in isotope dilution measurements of bile acid kinetics.* Pediatr Res 1975. [PMID 1103070](https://pubmed.ncbi.nlm.nih.gov/1103070/)
12. Javitt NB. *Bile acid synthesis from cholesterol: regulatory and auxiliary pathways.* FASEB J 1994. [PMID 8001744](https://pubmed.ncbi.nlm.nih.gov/8001744/)
13. Schwarz M, et al. *Two 7 alpha-hydroxylase enzymes in bile acid biosynthesis.* Curr Opin Lipidol 1998. [PMID 9559267](https://pubmed.ncbi.nlm.nih.gov/9559267/)
14. Sonne DP, et al. *Postprandial gallbladder emptying in patients with type 2 diabetes.* Eur J Endocrinol 2014. [PMID 24986531](https://pubmed.ncbi.nlm.nih.gov/24986531/)

---

## B. ASBT / 회장 수송 — `φ` 축의 분자적 실체

> 모델의 첫 번째 병변축 `φ`는 **회장 ASBT 흡수능의 생존 분율**입니다.
> Oelkers 1997의 SLC10A2 돌연변이는 `φ → 0`의 인간 실험이고,
> Dawson 2003의 Slc10a2 녹아웃은 그 동물판입니다. Hruz 2006은 회장
> 질환에서 ASBT가 실제로 감소함을 사람에서 보였습니다.

15. Oelkers P, et al. *Primary bile acid malabsorption caused by mutations in the ileal sodium-dependent bile acid transporter gene (SLC10A2).* J Clin Invest 1997. [PMID 9109432](https://pubmed.ncbi.nlm.nih.gov/9109432/)
16. Dawson PA, et al. *Targeted deletion of the ileal bile acid transporter eliminates enterohepatic cycling of bile acids in mice.* J Biol Chem 2003. [PMID 12819193](https://pubmed.ncbi.nlm.nih.gov/12819193/)
17. Dawson PA, et al. *Bile acid transporters.* J Lipid Res 2009. [PMID 19498215](https://pubmed.ncbi.nlm.nih.gov/19498215/)
18. Hruz P, et al. *Adaptive regulation of the ileal apical sodium dependent bile acid transporter (ASBT) in patients with obstructive cholestasis.* Gut 2006. [PMID 16150853](https://pubmed.ncbi.nlm.nih.gov/16150853/)
19. Balesaria S, et al. *Exploring possible mechanisms for primary bile acid malabsorption: evidence for different regulation of ileal bile acid transporter transcripts.* Eur J Gastroenterol Hepatol 2008. [PMID 18403943](https://pubmed.ncbi.nlm.nih.gov/18403943/)
20. Renner O, et al. *Mutation screening of apical sodium-dependent bile acid transporter (SLC10A2).* Hum Genet 2009. [PMID 19184108](https://pubmed.ncbi.nlm.nih.gov/19184108/)
21. Dawson PA, et al. *The heteromeric organic solute transporter alpha-beta, Ostalpha-Ostbeta, is an ileal basolateral bile acid transporter.* J Biol Chem 2005. [PMID 15563450](https://pubmed.ncbi.nlm.nih.gov/15563450/)
22. Rao A, et al. *The organic solute transporter alpha-beta, Ostalpha-Ostbeta, is essential for intestinal bile acid transport and homeostasis.* PNAS 2008. [PMID 18292224](https://pubmed.ncbi.nlm.nih.gov/18292224/)
23. Ballatori N, et al. *OSTalpha-OSTbeta: a major basolateral bile acid and steroid transporter.* Hepatology 2005. [PMID 16317684](https://pubmed.ncbi.nlm.nih.gov/16317684/)
24. Sultan M, et al. *Organic solute transporter-β (SLC51B) deficiency in two brothers with congenital diarrhea.* Hepatology 2018. [PMID 28898457](https://pubmed.ncbi.nlm.nih.gov/28898457/)
25. Wang L, et al. *Mechanism of Asbt (Slc10a2)-related bile acid malabsorption in diarrhea after pelvic radiation.* Int J Radiat Biol 2020. [PMID 31900034](https://pubmed.ncbi.nlm.nih.gov/31900034/)
26. Chothe PP, et al. *Human bile acid transporter ASBT (SLC10A2) forms functional non-covalent homodimers.* Biochim Biophys Acta 2018. [PMID 29198943](https://pubmed.ncbi.nlm.nih.gov/29198943/)

---

## C. FXR–FGF19–CYP7A1 되먹임 — `κ` 축과 되먹임 이득

> 이 모델의 중심 주장 — **센서는 흡수된 flux를 재고, 증상은 흘러넘친
> flux가 만든다** — 은 이 문헌군에서 나옵니다. Inagaki 2005가 회장
> FGF15/19가 CYP7A1을 억제하는 내분비 신호임을 확립했고, Lundåsen 2006이
> 사람에서 FGF19의 뚜렷한 일주기 변동과 담즙산 합성 조절을 보였습니다.
> Lu 2000 / Kerr 2002는 SHP 팔이 **보조적**(제거해도 되먹임이 남음)임을
> 보였고, 이는 모델에서 FGF19 팔에 더 큰 가중치(Rf가 지배)를 준 근거입니다.

27. Inagaki T, et al. *Fibroblast growth factor 15 functions as an enterohepatic signal to regulate bile acid homeostasis.* Cell Metab 2005. [PMID 16213224](https://pubmed.ncbi.nlm.nih.gov/16213224/)
28. Lundåsen T, et al. *Circulating intestinal fibroblast growth factor 19 has a pronounced diurnal variation and modulates hepatic bile acid synthesis in man.* J Intern Med 2006. [PMID 17116003](https://pubmed.ncbi.nlm.nih.gov/17116003/)
29. Lu TT, et al. *Molecular basis for feedback regulation of bile acid synthesis by nuclear receptors.* Mol Cell 2000. [PMID 11030331](https://pubmed.ncbi.nlm.nih.gov/11030331/)
30. Kerr TA, et al. *Loss of nuclear receptor SHP impairs but does not eliminate negative feedback regulation of bile acid synthesis.* Dev Cell 2002. [PMID 12062084](https://pubmed.ncbi.nlm.nih.gov/12062084/)
31. Makishima M, et al. *Identification of a nuclear receptor for bile acids.* Science 1999. [PMID 10334992](https://pubmed.ncbi.nlm.nih.gov/10334992/)
32. Fu T, et al. *FXR Primes the Liver for Intestinal FGF15 Signaling by Transient Induction of β-Klotho.* Mol Endocrinol 2016. [PMID 26505219](https://pubmed.ncbi.nlm.nih.gov/26505219/)
33. Byun S, et al. *Postprandial FGF19-induced phosphorylation by Src is critical for FXR function in bile acid homeostasis.* Nat Commun 2018. [PMID 29968724](https://pubmed.ncbi.nlm.nih.gov/29968724/)
34. Wang C, et al. *Hepatocyte FRS2α is essential for the endocrine fibroblast growth factor to limit the amplitude of bile acid production.* Curr Mol Med 2014. [PMID 25056539](https://pubmed.ncbi.nlm.nih.gov/25056539/)
35. Bouju A, et al. *A primer on the pleiotropic endocrine fibroblast growth factor FGF19/FGF15.* Differentiation 2024. [PMID 39500656](https://pubmed.ncbi.nlm.nih.gov/39500656/)
36. Ferrell JM, Chiang JYL. *Circadian rhythms in liver metabolism and disease.* Acta Pharm Sin B 2015. [PMID 26579436](https://pubmed.ncbi.nlm.nih.gov/26579436/)
37. Al-Khaifi A, et al. *Asynchronous rhythms of circulating conjugated and unconjugated bile acids.* J Intern Med 2018. [PMID 29964306](https://pubmed.ncbi.nlm.nih.gov/29964306/)
38. Haeusler RA, et al. *Impaired generation of 12-hydroxylated bile acids links hepatic insulin signaling with dyslipidemia.* Cell Metab 2012. [PMID 22197325](https://pubmed.ncbi.nlm.nih.gov/22197325/)
39. Patankar JV, et al. *Genetic ablation of Cyp8b1 preserves host metabolic function.* Am J Physiol Endocrinol Metab 2018. [PMID 29066462](https://pubmed.ncbi.nlm.nih.gov/29066462/)

---

## D. 담즙산 설사의 병태생리와 아형

> Type 1 (회장질환/절제, `φ↓`) · Type 2 (특발성/일차성, `κ↓`) ·
> Type 3 (담낭절제·SIBO·소아지방변증 등) · Type 4 (약물유발).
> Walters 2009 계열의 연구가 특발성 BAD에서 **공복 FGF19가 낮다**는
> 것을 보인 것이 `κ` 축을 도입한 직접적 근거입니다 (섹션 F의 Pattni 2012,
> Borup 2015 참조).

40. Camilleri M. *Bile Acid Diarrhea: Prevalence, Pathogenesis, and Therapy* / *The Role of Bile Acids in Chronic Diarrhea.* Am J Gastroenterol 2020. [PMID 32558690](https://pubmed.ncbi.nlm.nih.gov/32558690/)
41. Camilleri M, Vijayvargiya P. *Bile Acid Diarrhea in Adults and Adolescents.* Neurogastroenterol Motil 2022. [PMID 34751982](https://pubmed.ncbi.nlm.nih.gov/34751982/)
42. Walters JRF. *Managing bile acid diarrhea: aspects of contention.* Expert Rev Gastroenterol Hepatol 2024. [PMID 39264409](https://pubmed.ncbi.nlm.nih.gov/39264409/)
43. Barbara G, et al. *Bile acid diarrhea in patients with chronic diarrhea. Current appraisal and recommendations for clinical practice.* Dig Liver Dis 2025. [PMID 39827025](https://pubmed.ncbi.nlm.nih.gov/39827025/)
44. Di Ciaula A, et al. *Advances in the pathophysiology, diagnosis and management of chronic diarrhoea from bile acid malabsorption.* Eur J Intern Med 2024. [PMID 39069430](https://pubmed.ncbi.nlm.nih.gov/39069430/)
45. Barkun AN, et al. *Bile acid malabsorption in chronic diarrhea: pathophysiology and treatment.* Can J Gastroenterol 2013. [PMID 24199211](https://pubmed.ncbi.nlm.nih.gov/24199211/)
46. Gracie DJ, et al. *Prevalence of, and predictors of, bile acid malabsorption in outpatients with chronic diarrhea.* Neurogastroenterol Motil 2012. [PMID 22765392](https://pubmed.ncbi.nlm.nih.gov/22765392/)
47. Ruiz-Campos L, et al. *Systematic review with meta-analysis: the prevalence of bile acid malabsorption and response to colestyramine.* Aliment Pharmacol Ther 2019. [PMID 30585336](https://pubmed.ncbi.nlm.nih.gov/30585336/)
48. Farrugia A, et al. *Rates of Bile Acid Diarrhoea After Cholecystectomy: A Multicentre Audit.* World J Surg 2021. [PMID 33982189](https://pubmed.ncbi.nlm.nih.gov/33982189/)
49. Farrugia A, et al. *Bile acid diarrhoea and metabolic changes after cholecystectomy: a prospective case-control study.* BMC Gastroenterol 2024. [PMID 39174936](https://pubmed.ncbi.nlm.nih.gov/39174936/)
50. Barrera F, et al. *Effect of cholecystectomy on bile acid synthesis and circulating levels of fibroblast growth factor 19.* Ann Hepatol 2016. [PMID 26256900](https://pubmed.ncbi.nlm.nih.gov/26256900/)
51. Bouchoucha M, et al. *Metformin and digestive disorders.* Diabetes Metab 2011. [PMID 21236717](https://pubmed.ncbi.nlm.nih.gov/21236717/)

---

## E. 대장 분비·운동 — 3 mM 역치와 양성 되먹임

> 모델의 `EC50_G = 3 mM CDCA-equivalent`는 **Mekjian 1971**의 사람 대장
> 관류 실험에서 직접 옵니다: 이수산기 담즙산(CDCA·DCA)이 3 mM을 넘으면
> 순 분비로 전환됩니다. 삼수산기(CA)는 훨씬 약합니다 → 모델의 potency
> 가중치 (CA 0.15 / CDCA 1.00 / DCA 1.15 / LCA 0.30).
> Farack 1984는 loperamide가 DCA 유발 분비 자체를 억제함을 보였고,
> 이는 모델에서 loperamide를 "대증"이 아니라 **loop gain 저하제**로
> 분류한 근거입니다.

52. Mekjian HS, Phillips SF, Hofmann AF. *Colonic secretion of water and electrolytes induced by bile acids: perfusion studies in man.* J Clin Invest 1971. [PMID 4938344](https://pubmed.ncbi.nlm.nih.gov/4938344/)
53. Rao AS, et al. *Chenodeoxycholate in females with irritable bowel syndrome-constipation: a pharmacodynamic and pharmacogenetic analysis.* Gastroenterology 2010. [PMID 20691689](https://pubmed.ncbi.nlm.nih.gov/20691689/)
54. Odunsi-Shiyanbade ST, et al. *Effects of chenodeoxycholate and a bile acid sequestrant, colesevelam, on intestinal transit and bowel function.* Clin Gastroenterol Hepatol 2010. [PMID 19879973](https://pubmed.ncbi.nlm.nih.gov/19879973/)
55. Alemi F, et al. *The receptor TGR5 mediates the prokinetic actions of intestinal bile acids and is required for normal defecation.* Gastroenterology 2013. [PMID 23041323](https://pubmed.ncbi.nlm.nih.gov/23041323/)
56. Ward JB, et al. *The bile acid receptor, TGR5, regulates basal and cholinergic-induced secretory responses in rat colon.* Neurogastroenterol Motil 2013. [PMID 23634890](https://pubmed.ncbi.nlm.nih.gov/23634890/)
57. Poole DP, et al. *Expression and function of the bile acid receptor GpBAR1 (TGR5) in the murine enteric nervous system.* Neurogastroenterol Motil 2010. [PMID 20236244](https://pubmed.ncbi.nlm.nih.gov/20236244/)
58. Bunnett NW. *Neuro-humoral signalling by bile acids and the TGR5 receptor in the gastrointestinal tract.* J Physiol 2014. [PMID 24614746](https://pubmed.ncbi.nlm.nih.gov/24614746/)
59. Yu Y, et al. *Deoxycholic acid activates colonic afferent nerves via 5-HT3 receptor-dependent and -independent mechanisms.* Am J Physiol Gastrointest Liver Physiol 2019. [PMID 31216174](https://pubmed.ncbi.nlm.nih.gov/31216174/)
60. Camilleri M, et al. *Bile acid detergency: permeability, inflammation, and effects of sulfation.* Am J Physiol Gastrointest Liver Physiol 2022. [PMID 35258349](https://pubmed.ncbi.nlm.nih.gov/35258349/)
61. Zeng H, et al. *Deoxycholic Acid Modulates Cell-Junction Gene Expression and Increases Intestinal Barrier Dysfunction.* Molecules 2022. [PMID 35163990](https://pubmed.ncbi.nlm.nih.gov/35163990/)
62. Farack UM, Loeschke K. *Inhibition by loperamide of deoxycholic acid induced intestinal secretion.* Naunyn Schmiedebergs Arch Pharmacol 1984. [PMID 6328335](https://pubmed.ncbi.nlm.nih.gov/6328335/)

---

## F. 진단검사 — SeHCAT · C4 · FGF19 · 대변 담즙산

> **모델의 가장 검정 가능한 예측이 이 섹션과 충돌하거나 부합합니다.**
> ⚙ 모델은 SeHCAT이 `φ`만을 (그것도 매우 좁은 구간에서 압축적으로) 읽고,
> C4·대변 담즙산은 `S = f(φ, κ)`를 읽는다고 주장합니다. 그래서 두 검사가
> 불완전하게만 일치하고 (Valentin 2016 메타분석), SeHCAT 정상인데
> C4가 높은 환자가 존재해야 합니다 (Vijayvargiya 2017·2020, Borup 2024).

63. Riemsma R, et al. *SeHCAT [tauroselcholic (selenium-75) acid] for the investigation of bile acid malabsorption: a systematic review and cost-effectiveness analysis.* Health Technol Assess 2013. [PMID 24351663](https://pubmed.ncbi.nlm.nih.gov/24351663/)
64. Wedlake L, et al. *Systematic review: the prevalence of idiopathic bile acid malabsorption as diagnosed by SeHCAT scanning in patients with diarrhoea-predominant IBS.* Aliment Pharmacol Ther 2009. [PMID 19570102](https://pubmed.ncbi.nlm.nih.gov/19570102/)
65. van Tilburg AJ, et al. *The selenium-75-homocholic acid taurine test reevaluated: combined measurement of fecal selenium-75 activity and 3 alpha-hydroxy bile acids.* J Nucl Med 1991. [PMID 2045936](https://pubmed.ncbi.nlm.nih.gov/2045936/)
66. Pattni SS, et al. *Fibroblast Growth Factor 19 and 7α-Hydroxy-4-Cholesten-3-one in the Diagnosis of Patients With Possible Bile Acid Diarrhea.* Clin Transl Gastroenterol 2012. [PMID 23238290](https://pubmed.ncbi.nlm.nih.gov/23238290/)
67. Vijayvargiya P, et al. *Performance characteristics of serum C4 and FGF19 measurements to exclude the diagnosis of bile acid diarrhoea.* Aliment Pharmacol Ther 2017. [PMID 28691284](https://pubmed.ncbi.nlm.nih.gov/28691284/)
68. Vijayvargiya P, et al. *Fecal Bile Acid Testing in Assessing Patients With Chronic Unexplained Diarrhea.* Am J Gastroenterol 2020. [PMID 32618660](https://pubmed.ncbi.nlm.nih.gov/32618660/)
69. Valentin N, et al. *Biomarkers for bile acid diarrhoea in functional bowel disorder with diarrhoea: a systematic review and meta-analysis.* Gut 2016. [PMID 26347530](https://pubmed.ncbi.nlm.nih.gov/26347530/)
70. Borup C, et al. *Diagnosis of bile acid diarrhoea by fasting and postprandial measurements of fibroblast growth factor 19.* Eur J Gastroenterol Hepatol 2015. [PMID 26426834](https://pubmed.ncbi.nlm.nih.gov/26426834/)
71. Borup C, et al. *Prospective comparison of diagnostic tests for bile acid diarrhoea.* Aliment Pharmacol Ther 2024. [PMID 37794830](https://pubmed.ncbi.nlm.nih.gov/37794830/)
72. Battat R, et al. *Serum Concentrations of 7α-hydroxy-4-cholesten-3-one Are Associated With Bile Acid Diarrhea in Patients With Crohn's Disease.* Clin Gastroenterol Hepatol 2019. [PMID 30448597](https://pubmed.ncbi.nlm.nih.gov/30448597/)
73. Lyutakov I, et al. *Methods for diagnosing bile acid malabsorption: a systematic review.* BMC Gastroenterol 2019. [PMID 31726982](https://pubmed.ncbi.nlm.nih.gov/31726982/)
74. Dilmaghani S, et al. *Simplifying Diagnosis of Bile Acid Diarrhea With Clinical and Biochemical Measurements on Blood and Single Stool Sample.* Clin Gastroenterol Hepatol 2025. [PMID 40378992](https://pubmed.ncbi.nlm.nih.gov/40378992/)
75. Peleman C, et al. *Colonic Transit and Bile Acid Synthesis or Excretion in Patients With Irritable Bowel Syndrome-Diarrhea Without Bile Acid Malabsorption.* Clin Gastroenterol Hepatol 2017. [PMID 27856362](https://pubmed.ncbi.nlm.nih.gov/27856362/)
76. Lupianez-Merly C, Camilleri M. *Recent developments in diagnosing bile acid diarrhea.* Expert Rev Gastroenterol Hepatol 2023. [PMID 38086533](https://pubmed.ncbi.nlm.nih.gov/38086533/)

---

## G. 담즙산 결합수지 (class A — 관강 부하를 묶는다)

> ⚙ 모델은 결합수지가 **자기가 치료하려는 구동력을 스스로 올린다**고
> 예측합니다: 관강 담즙산을 묶으면 센서 리간드도 사라져 FGF19가 떨어지고
> CYP7A1/C4가 올라갑니다. 이 "구조적 탈출"은 결합수지의 용량-반응이
> 얕은 이유이고, FXR 작용제와의 병용이 상승적인 이유입니다.
> 📖 결합수지가 담즙산 합성을 실제로 증가시킨다는 것은 지질 문헌에서
> 오래된 사실입니다 (LDL 강하 기전 자체가 그것입니다).

77. Borup C, et al. *Efficacy and safety of colesevelam for the treatment of bile acid diarrhoea: a double-blind, randomised, placebo-controlled, phase 4 clinical trial.* Lancet Gastroenterol Hepatol 2023. [PMID 36758570](https://pubmed.ncbi.nlm.nih.gov/36758570/)
78. Camilleri M, et al. *Effect of colesevelam on faecal bile acids and bowel functions in diarrhoea-predominant irritable bowel syndrome.* Aliment Pharmacol Ther 2015. [PMID 25594801](https://pubmed.ncbi.nlm.nih.gov/25594801/)
79. Beigel F, et al. *Colesevelam for the treatment of bile acid malabsorption-associated diarrhea in patients with Crohn's disease.* J Crohns Colitis 2014. [PMID 24953836](https://pubmed.ncbi.nlm.nih.gov/24953836/)
80. Soares GAR, et al. *Efficacy of Bile Acid Sequestrants in the Treatment of Bile Acid Diarrhea: A Meta-Analysis of Randomized Controlled Trials.* J Clin Pharmacol 2025. [PMID 39428959](https://pubmed.ncbi.nlm.nih.gov/39428959/)
81. Wilcox C, et al. *Systematic review: the management of chronic diarrhoea due to bile acid malabsorption.* Aliment Pharmacol Ther 2014. [PMID 24602022](https://pubmed.ncbi.nlm.nih.gov/24602022/)
82. Kårhus ML, et al. *Safety and efficacy of liraglutide versus colesevelam for the treatment of bile acid diarrhoea: a randomised, double-blind, active-comparator, non-inferiority clinical trial.* Lancet Gastroenterol Hepatol 2022. [PMID 35868334](https://pubmed.ncbi.nlm.nih.gov/35868334/)
83. Ellegaard AM, et al. *Liraglutide and Colesevelam Change Serum and Fecal Bile Acid Levels in a Randomized Trial.* Clin Transl Gastroenterol 2024. [PMID 39602188](https://pubmed.ncbi.nlm.nih.gov/39602188/)
84. Scaldaferri F, et al. *Use and indications of cholestyramine and bile acid sequestrants.* Intern Emerg Med 2013. [PMID 21739227](https://pubmed.ncbi.nlm.nih.gov/21739227/)

---

## H. FXR 작용제 (class B — 사라진 신호를 공급한다)

> ⚙ 모델의 분류: FXR 작용제는 수용체 점유에서 **덧셈**이므로 `φ`와
> 무관하게 작동하며, 유일하게 `S`(=대장 부하) 자체를 낮춥니다.
> ⚙ 예측: FXR 작용제는 회장을 고치지 않고도 SeHCAT 수치를 올립니다
> (풀 회전율 `S/P`가 떨어지므로). 📖 Walters 2015(OBADIAH)와
> Camilleri 2020(tropifexor)이 C4 급감·FGF19 상승·대장통과 지연을
> 사람에서 보였습니다. 📖 대가는 LDL-C 상승입니다 (Siddiqui 2020, Pencek 2016).

85. Walters JR, et al. *The response of patients with bile acid diarrhoea to the farnesoid X receptor agonist obeticholic acid.* Aliment Pharmacol Ther 2015. [PMID 25329562](https://pubmed.ncbi.nlm.nih.gov/25329562/)
86. Camilleri M, et al. *Randomised clinical trial: significant biochemical and colonic transit effects of the farnesoid X receptor agonist tropifexor in patients with primary bile acid diarrhoea.* Aliment Pharmacol Ther 2020. [PMID 32702169](https://pubmed.ncbi.nlm.nih.gov/32702169/)
87. Pellicciari R, et al. *6alpha-ethyl-chenodeoxycholic acid (6-ECDCA), a potent and selective FXR agonist.* J Med Chem 2002. [PMID 12166927](https://pubmed.ncbi.nlm.nih.gov/12166927/)
88. Jiang L, et al. *Structural basis of tropifexor as a potent and selective agonist of farnesoid X receptor.* Biochem Biophys Res Commun 2021. [PMID 33121679](https://pubmed.ncbi.nlm.nih.gov/33121679/)
89. Gege C, et al. *Nonsteroidal FXR Ligands: Current Status and Clinical Applications.* Handb Exp Pharmacol 2019. [PMID 31197565](https://pubmed.ncbi.nlm.nih.gov/31197565/)
90. Fiorucci S, et al. *Obeticholic Acid: An Update of Its Pharmacological Activities in Liver Disorders.* Handb Exp Pharmacol 2019. [PMID 31201552](https://pubmed.ncbi.nlm.nih.gov/31201552/)
91. Siddiqui MS, et al. *Impact of obeticholic acid on the lipoprotein profile in patients with non-alcoholic steatohepatitis.* J Hepatol 2020. [PMID 31634532](https://pubmed.ncbi.nlm.nih.gov/31634532/)
92. Pencek R, et al. *Effects of obeticholic acid on lipoprotein metabolism in healthy volunteers.* Diabetes Obes Metab 2016. [PMID 27109453](https://pubmed.ncbi.nlm.nih.gov/27109453/)
93. Kremoser C. *FXR agonists for NASH: How are they different and what difference do they make?* J Hepatol 2021. [PMID 33985820](https://pubmed.ncbi.nlm.nih.gov/33985820/)

---

## I. ASBT 억제제 (class C — `φ`를 일부러 낮춘다) · 거울상

> ⚙ 모델의 **위조 가능성 시험**: ASBT 억제제는 정상 장에서 완전한 BAM
> 표현형을 만들어야 하며, 그것이 elobixibat이 만성변비 치료제인 이유여야
> 합니다. 📖 Simrén 2011 · Nakajima/Miner 계열 시험이 그대로입니다.
> 📖 Carreño 2025는 IBAT 억제제 관련 설사를 **C4 농도로 예측**하는
> 약동-약력 분석으로, 이 모델의 "C4는 `S`를 읽는다"는 축과 정확히 같은
> 논리를 임상시험 데이터에서 사용합니다.

94. Simrén M, et al. *Randomised clinical trial: the ileal bile acid transporter inhibitor A3309 vs. placebo in patients with chronic idiopathic constipation.* Aliment Pharmacol Ther 2011. [PMID 21545606](https://pubmed.ncbi.nlm.nih.gov/21545606/)
95. Miner PB Jr. *Elobixibat, the first-in-class Ileal Bile Acid Transporter inhibitor, for the treatment of Chronic Idiopathic Constipation.* Expert Opin Pharmacother 2018. [PMID 30129377](https://pubmed.ncbi.nlm.nih.gov/30129377/)
96. Taniguchi S, et al. *Elobixibat, an ileal bile acid transporter inhibitor, induces giant migrating contractions during natural defecation.* Neurogastroenterol Motil 2018. [PMID 30129138](https://pubmed.ncbi.nlm.nih.gov/30129138/)
97. Manabe N, et al. *Elobixibat Improves Stool/Gas Distribution and Fecal Bile Acids in Older Adults With Chronic Constipation.* JGH Open 2025. [PMID 40832007](https://pubmed.ncbi.nlm.nih.gov/40832007/)
98. Carreño F, et al. *Analysis of C4 Concentrations to Predict Impact of Patient-Reported Diarrhea Associated With the Ileal Bile Acid Transporter Inhibitor...* CPT Pharmacometrics Syst Pharmacol 2025. [PMID 39945351](https://pubmed.ncbi.nlm.nih.gov/39945351/)
99. Al-Dury S, Marschall HU. *Ileal Bile Acid Transporter Inhibition for the Treatment of Chronic Constipation, Cholestatic Pruritus, and NASH.* Front Pharmacol 2018. [PMID 30186169](https://pubmed.ncbi.nlm.nih.gov/30186169/)
100. Shirley M. *Maralixibat: First Approval.* Drugs 2022. [PMID 34813049](https://pubmed.ncbi.nlm.nih.gov/34813049/)
101. Peverelle M, et al. *Review Article: Ileal Bile Acid Transport (IBAT) Inhibitors as an Emerging Treatment for Cholestatic Liver Disease.* Aliment Pharmacol Ther 2026. [PMID 41953994](https://pubmed.ncbi.nlm.nih.gov/41953994/)

---

## J. 회장절제 · 100 cm 규칙 · 지방변

> ⚙ 모델은 **100 cm 규칙을 데이터에 맞추지 않고 유도**합니다:
> ASBT 밀도의 원위 편중(감쇠거리 60 cm) + CYP7A1 보상 천장(6배) +
> 임계 미셀 농도(1.5 mM) 세 가지만으로 담즙산설사→지방변 전환점이
> 나옵니다. 📖 Hofmann & Poley 1972 / Poley & Hofmann 1976이 그 임상
> 원전이며, Skouras 2019는 절제 길이와 BAM 중증도의 상관을 보였습니다.

102. Hofmann AF, Poley JR. *Role of bile acid malabsorption in pathogenesis of diarrhea and steatorrhea in patients with ileal resection.* Gastroenterology 1972. [PMID 5029077](https://pubmed.ncbi.nlm.nih.gov/5029077/)
103. Poley JR, Hofmann AF. *Role of fat maldigestion in pathogenesis of steatorrhea in ileal resection.* Gastroenterology 1976. [PMID 6360](https://pubmed.ncbi.nlm.nih.gov/6360/)
104. Skouras T, et al. *Brief report: length of ileal resection correlates with severity of bile acid malabsorption in Crohn's disease.* Int J Colorectal Dis 2019. [PMID 30116880](https://pubmed.ncbi.nlm.nih.gov/30116880/)
105. Akerlund JE, et al. *Hepatic metabolism of cholesterol in Crohn's disease. Effect of partial resection of ileum.* Gastroenterology 1991. [PMID 2001802](https://pubmed.ncbi.nlm.nih.gov/2001802/)
106. Färkkilä MA, et al. *Plasma lathosterol as a screening test for bile acid malabsorption due to ileal resection.* Clin Sci (Lond) 1996. [PMID 8777839](https://pubmed.ncbi.nlm.nih.gov/8777839/)
107. Siener R, et al. *Intestinal Oxalate Absorption, Enteric Hyperoxaluria, and Risk of Urinary Stone Formation.* Nutrients 2024. [PMID 38257157](https://pubmed.ncbi.nlm.nih.gov/38257157/)
108. Larsen HM, et al. *Chronic loose stools following right-sided hemicolectomy for colon cancer and the association with bile acid malabsorption.* Colorectal Dis 2023. [PMID 36347822](https://pubmed.ncbi.nlm.nih.gov/36347822/)

---

## K. 미생물총 — 탈포합과 7α-탈수산화 (class E — 부하가 아니라 potency)

> ⚙ 모델의 두 번째 불일치 설명: 미생물총은 대장에 도달하는 담즙산의
> **양**을 (거의) 바꾸지 않고 **분비 효력**을 바꿉니다. CA(w 0.15)를
> DCA(w 1.15)로 바꾸는 7α-탈수산화가 율속 단계이므로, 같은 대변 담즙산
> 배설량에서도 증상이 3배 차이날 수 있습니다.

109. Ridlon JM, et al. *Bile salt biotransformations by human intestinal bacteria.* J Lipid Res 2006. [PMID 16299351](https://pubmed.ncbi.nlm.nih.gov/16299351/)
110. Wells JE, Hylemon PB. *Identification and characterization of a bile acid 7alpha-dehydroxylation operon in Clostridium sp. strain TO-931.* Appl Environ Microbiol 2000. [PMID 10698778](https://pubmed.ncbi.nlm.nih.gov/10698778/)
111. Ridlon JM, et al. *Identification and characterization of two bile acid coenzyme A transferases from Clostridium scindens.* J Lipid Res 2012. [PMID 22021638](https://pubmed.ncbi.nlm.nih.gov/22021638/)
112. Ridlon JM, et al. *Consequences of bile salt biotransformations by intestinal bacteria.* Gut Microbes 2016. [PMID 26939849](https://pubmed.ncbi.nlm.nih.gov/26939849/)
113. Ridlon JM, et al. *The human gut sterolbiome: bile acid-microbiome endocrine aspects and therapeutics.* Acta Pharm Sin B 2015. [PMID 26579434](https://pubmed.ncbi.nlm.nih.gov/26579434/)
114. Guzior DV, et al. *Bile salt hydrolase acyltransferase activity expands bile acid diversity.* Nature 2024. [PMID 38326608](https://pubmed.ncbi.nlm.nih.gov/38326608/)
115. Sinha SR, et al. *Dysbiosis-Induced Secondary Bile Acid Deficiency Promotes Intestinal Inflammation.* Cell Host Microbe 2020. [PMID 32101703](https://pubmed.ncbi.nlm.nih.gov/32101703/)
116. Camilleri M, et al. *Comparison of biochemical, microbial and mucosal mRNA expression in bile acid diarrhoea and irritable bowel syndrome with diarrhoea.* Gut 2023. [PMID 35580964](https://pubmed.ncbi.nlm.nih.gov/35580964/)

---

## L. IBS-D 중첩과 대증치료 (class D — loop gain을 낮춘다)

> ⚙ 모델에서 loperamide·ondansetron은 단순 대증요법이 아닙니다:
> 통과속도를 늦추면 회장 접촉시간이 늘어 `f`가 올라가고 유출 자체가
> 줄어듭니다. 그래서 모델은 loperamide 하에서 **대변 담즙산 배설량도**
> (수분뿐 아니라) 줄어들 것을 예측합니다.

117. Camilleri M. *Irritable bowel syndrome: treatment based on pathophysiology and biomarkers.* Gut 2023. [PMID 36307180](https://pubmed.ncbi.nlm.nih.gov/36307180/)
118. Garsed K, et al. *A randomised trial of ondansetron for the treatment of irritable bowel syndrome with diarrhoea.* Gut 2014. [PMID 24334242](https://pubmed.ncbi.nlm.nih.gov/24334242/)
119. Gunn D, et al. *Randomised, placebo-controlled trial and meta-analysis show benefit of ondansetron for irritable bowel syndrome with diarrhoea (TRITON).* Aliment Pharmacol Ther 2023. [PMID 36866724](https://pubmed.ncbi.nlm.nih.gov/36866724/)
120. Sun WM, et al. *Effects of loperamide oxide on gastrointestinal transit time and anorectal function in patients with chronic diarrhoea.* Scand J Gastroenterol 1997. [PMID 9018764](https://pubmed.ncbi.nlm.nih.gov/9018764/)
121. Mainguet P, Fiasse R. *Double-blind placebo-controlled study of loperamide (Imodium) in chronic diarrhoea caused by ileocolic disease or resection.* Gut 1977. [PMID 326642](https://pubmed.ncbi.nlm.nih.gov/326642/)
122. Vijayvargiya P, et al. *Safety and Efficacy of Eluxadoline in Patients with Irritable Bowel Syndrome-Diarrhea With or Without Bile Acid Diarrhea.* Dig Dis Sci 2022. [PMID 35122592](https://pubmed.ncbi.nlm.nih.gov/35122592/)
123. Shin A, et al. *Bowel functions, fecal unconjugated primary and secondary bile acids, and colonic transit in patients with irritable bowel syndrome.* Clin Gastroenterol Hepatol 2013. [PMID 23639599](https://pubmed.ncbi.nlm.nih.gov/23639599/)
124. Vijayvargiya P, et al. *Bile and fat excretion are biomarkers of clinically significant diarrhoea and constipation in irritable bowel syndrome.* Aliment Pharmacol Ther 2019. [PMID 30740753](https://pubmed.ncbi.nlm.nih.gov/30740753/)

---

## M. 담즙산 QSP·PBPK 모델링과 도구

> 이 모델의 방법론적 선조들. Hofmann 1983(섹션 A-1)이 원형이고,
> Voronova 2020이 현대적 PBPK 재구성, Guiastrennec 2018이 인체 혈장
> 담즙산의 장간순환 모델, Meessen 2020이 개인별 식후 반응 모델입니다.

125. Voronova V, et al. *A Physiology-Based Model of Bile Acid Distribution and Metabolism Under Healthy and Pathologic Conditions.* Cell Mol Gastroenterol Hepatol 2020. [PMID 32112828](https://pubmed.ncbi.nlm.nih.gov/32112828/)
126. Guiastrennec B, et al. *Model-Based Prediction of Plasma Concentration and Enterohepatic Circulation of Total Bile Acids in Humans.* CPT Pharmacometrics Syst Pharmacol 2018. [PMID 30070437](https://pubmed.ncbi.nlm.nih.gov/30070437/)
127. Meessen ECE, et al. *Model-based data analysis of individual human postprandial plasma bile acid responses.* Physiol Rep 2020. [PMID 32170845](https://pubmed.ncbi.nlm.nih.gov/32170845/)
128. Mayo AK, et al. *Quantitative Systems Toxicology Model Predicts Obeticholic Acid-Associated Liver Injury.* Clin Pharmacol Ther 2026. [PMID 42332345](https://pubmed.ncbi.nlm.nih.gov/42332345/)
129. Generaux G, et al. *Quantitative systems toxicology (QST) reproduces species differences in PF-04895162 liver safety due to combined mitochondrial and bile acid toxicity.* Pharmacol Res Perspect 2019. [PMID 31624633](https://pubmed.ncbi.nlm.nih.gov/31624633/)
130. Nigam SK. *The Systems Biology of Drug Metabolizing Enzymes and Transporters: Relevance to Quantitative Systems Pharmacology.* Clin Pharmacol Ther 2020. [PMID 32119114](https://pubmed.ncbi.nlm.nih.gov/32119114/)
131. Elmokadem A, et al. *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.* CPT Pharmacometrics Syst Pharmacol 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
132. Lu T, et al. *gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve.* CPT Pharmacometrics Syst Pharmacol 2024. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
133. Ghallab A, et al. *Enteronephrohepatic circulation of bile acids and therapeutic potential of systemic bile acid transporter inhibition.* J Hepatol 2025. [PMID 40414504](https://pubmed.ncbi.nlm.nih.gov/40414504/)
134. De Bruijn VMP, et al. *Intestinal in vitro transport assay combined with physiologically based kinetic modeling.* ALTEX 2024. [PMID 37528756](https://pubmed.ncbi.nlm.nih.gov/37528756/)

---

## 모델의 주장과 그 근거의 대응표

| 모델의 주장 | 근거 유형 | 문헌 |
|---|---|---|
| 정상 풀 3 g, 하루 4–6회 순환, 합성 0.35 g/day | 📖 보정 | 1, 7, 8, 10 |
| 회장 보존율 f ≈ 0.97, ASBT는 원위 편중 | 📖 보정 | 15–18, 102 |
| 대장 분비 역치 ≈ 3 mM 이수산기 담즙산 | 📖 파라미터 직접 | **52** (Mekjian 1971) |
| CA는 약하고 DCA가 가장 강한 분비자극제 | 📖 potency 가중치 | 52, 55, 60 |
| FGF19가 CYP7A1 억제의 주 팔, SHP는 보조 | 📖 구조 | 27, 28, 29, 30 |
| 특발성(2형) BAD = 낮은 공복 FGF19 (`κ↓`) | 📖 병변축 | 66, 67, 70 |
| **정상상태에서 대장 부하 = 간 합성** | ⚙ 모델이 유도 | (질량보존) |
| **되먹임을 얼리면 흡수장애만으로는 설사가 없다** | ⚙ 모델이 유도 | 검정 대상 |
| **SeHCAT은 `φ`만, C4/대변담즙산은 `S`를 읽는다** | ⚙ 모델이 유도 | 69, 67, 68과 부합 |
| **SeHCAT 음성 담즙산설사가 존재해야 한다** | ⚙ 예측 | 68, 71과 부합 |
| 결합수지는 C4를 올린다 (구조적 탈출) | ⚙ + 📖 | 77, 78 |
| **FXR 작용제는 회장을 고치지 않고 SeHCAT을 올린다** | ⚙ 예측 (미검증) | 85, 86 |
| ASBT 억제제 = 의원성 BAM | ⚙ + 📖 | 94–98 |
| **100 cm 규칙은 세 상수에서 유도된다** | ⚙ 모델이 유도 | 102, 103, 104와 부합 |
| loperamide는 대변 담즙산도 줄인다 | ⚙ 예측 | 62, 120, 121과 부합 |
| 미생물총은 부하가 아니라 potency를 바꾼다 | ⚙ 예측 | 109–113, 116 |

---

*문헌 134편 · PubMed 링크 134개 · 모든 PMID는 NCBI E-utilities로 조회·검증됨.*
