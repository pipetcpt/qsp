# 고지단백(a)혈증 — 참고문헌
# Elevated Lipoprotein(a) — Reference List

이 목록은 `lpa_qsp_model.dot` / `lpa_mrgsolve_model.R` 의 구조와 파라미터
근거를 섹션별로 정리한 것입니다. 각 항목의 PMID 링크는 PubMed로 연결됩니다.
모델의 각 방정식이 어느 문헌에 근거하는지는 §13 의 대응표를 참조하십시오.

---

## 1. 발견 · 구조 · 진화 (Discovery, structure, evolution)

1. Berg K. **A new serum type system in man — the Lp system.** *Acta Pathol Microbiol Scand.* 1963;59:369-82. — Lp(a)의 최초 기술. [PMID 14064818](https://pubmed.ncbi.nlm.nih.gov/14064818/)
2. McLean JW, Tomlinson JE, Kuang WJ, et al. **cDNA sequence of human apolipoprotein(a) is homologous to plasminogen.** *Nature.* 1987;330(6144):132-7. — apo(a)가 플라스미노겐 상동체임을 밝힌 논문. 모델 ARM 3 전체의 근거. [PMID 3670400](https://pubmed.ncbi.nlm.nih.gov/3670400/)
3. Utermann G. **The mysteries of lipoprotein(a).** *Science.* 1989;246(4932):904-10. [PMID 2530631](https://pubmed.ncbi.nlm.nih.gov/2530631/)
4. Lawn RM, Schwartz K, Patthy L. **Convergent evolution of apolipoprotein(a) in primates and hedgehog.** *Proc Natl Acad Sci USA.* 1997;94(22):11992-7. — Old World 영장류와 고슴도치에만 존재. [PMID 9342350](https://pubmed.ncbi.nlm.nih.gov/9342350/)
5. Koschinsky ML, Marcovina SM. **Structure-function relationships in apolipoprotein(a): insights into lipoprotein(a) assembly and pathogenicity.** *Curr Opin Lipidol.* 2004;15(2):167-74. [PMID 15017359](https://pubmed.ncbi.nlm.nih.gov/15017359/)
6. Schmidt K, Noureen A, Kronenberg F, Utermann G. **Structure, function, and genetics of lipoprotein(a).** *J Lipid Res.* 2016;57(8):1339-59. — 종합 리뷰. [PMID 27074913](https://pubmed.ncbi.nlm.nih.gov/27074913/)

## 2. 유전학 · KIV-2 반복수 (Genetics and the KIV-2 copy-number variant)

7. Utermann G, Menzel HJ, Kraft HG, et al. **Lp(a) glycoprotein phenotypes. Inheritance and relation to Lp(a)-lipoprotein concentrations in plasma.** *J Clin Invest.* 1987;80(2):458-65. — 아이소폼 크기와 혈중 농도의 역상관. 모델 `SECEFF` 항의 근거. [PMID 2956279](https://pubmed.ncbi.nlm.nih.gov/2956279/)
8. Boerwinkle E, Leffert CC, Lin J, et al. **Apolipoprotein(a) gene accounts for greater than 90% of the variation in plasma lipoprotein(a) concentrations.** *J Clin Invest.* 1992;90(1):52-60. — 모델이 Lp(a)를 생산 결정 변수로 두는 이유. [PMID 1386087](https://pubmed.ncbi.nlm.nih.gov/1386087/)
9. Kraft HG, Köchl S, Menzel HJ, Sandholzer C, Utermann G. **The apolipoprotein(a) gene: a transcribed hypervariable locus controlling plasma lipoprotein(a) concentration.** *Hum Genet.* 1992;90(3):220-30. [PMID 1487235](https://pubmed.ncbi.nlm.nih.gov/1487235/)
10. Clarke R, Peden JF, Hopewell JC, et al. **Genetic variants associated with Lp(a) lipoprotein level and coronary disease.** *N Engl J Med.* 2009;361(26):2518-28. — rs10455872 · rs3798220. [PMID 20032323](https://pubmed.ncbi.nlm.nih.gov/20032323/)
11. Kamstrup PR, Tybjaerg-Hansen A, Steffensen R, Nordestgaard BG. **Genetically elevated lipoprotein(a) and increased risk of myocardial infarction.** *JAMA.* 2009;301(22):2331-9. — 멘델 무작위화. [PMID 19509380](https://pubmed.ncbi.nlm.nih.gov/19509380/)
12. Mukamel RE, Handsaker RE, Sherman MA, et al. **Protein-coding repeat polymorphisms strongly modulate plasma protein levels.** *Science.* 2021;373(6562):1499-1505. — KIV-2 반복수의 정량적 효과. [PMID 34554798](https://pubmed.ncbi.nlm.nih.gov/34554798/)
13. Coassin S, Kronenberg F. **Lipoprotein(a) beyond the kringle IV repeat polymorphism: from genotype to phenotype.** *Atherosclerosis.* 2022;349:17-35. [PMID 35606073](https://pubmed.ncbi.nlm.nih.gov/35606073/)
14. Trinder M, Uddin MM, Finneran P, Aragam KG, Natarajan P. **Clinical utility of lipoprotein(a) and LPA genetic risk score in risk prediction of incident atherosclerotic cardiovascular disease.** *JAMA Cardiol.* 2021;6(3):287-95. [PMID 33021622](https://pubmed.ncbi.nlm.nih.gov/33021622/)
15. Guan W, Cao J, Steffen BT, et al. **Race is a key variable in assigning lipoprotein(a) cutoff values for coronary heart disease risk assessment: the MESA study.** *Arterioscler Thromb Vasc Biol.* 2015;35(4):996-1001. — 모델 `FANC`(인종) 공변량. [PMID 25810300](https://pubmed.ncbi.nlm.nih.gov/25810300/)
16. Patel AP, Wang M, Pirruccello JP, et al. **Lp(a) (lipoprotein[a]) concentrations and incident atherosclerotic cardiovascular disease: new insights from a large national biobank.** *Arterioscler Thromb Vasc Biol.* 2021;41(1):465-74. [PMID 33115266](https://pubmed.ncbi.nlm.nih.gov/33115266/)

## 3. 생합성 · 분비 효율 · 조립 (Biosynthesis, secretion efficiency, assembly)

17. White AL, Lanford RE. **Cell surface assembly of lipoprotein(a) in primary cultures of baboon hepatocytes.** *J Biol Chem.* 1994;269(46):28716-23. — 조립이 세포 표면/세포외에서 일어난다는 근거. 모델이 free apo(a)를 별도 상태변수로 두는 이유. [PMID 7961824](https://pubmed.ncbi.nlm.nih.gov/7961824/)
18. White AL, Hixson JE, Rainwater DL, Lanford RE. **Molecular basis for "null" lipoprotein(a) phenotypes and the influence of apolipoprotein(a) size on plasma lipoprotein(a) level in the baboon.** *J Biol Chem.* 1994;269(12):9060-6. — 큰 아이소폼의 분비 전 분해(presecretory degradation). 모델 `KDEGP` 항의 직접 근거. [PMID 8132644](https://pubmed.ncbi.nlm.nih.gov/8132644/)
19. Brunner C, Lobentanz EM, Pethö-Schramm A, et al. **The number of identical kringle IV repeats in apolipoprotein(a) affects its processing and secretion by HepG2 cells.** *J Biol Chem.* 1996;271(50):32403-10. — 반복수 → 소포체 체류시간 → 분비 효율. 모델 `SECEFF = KSZ^3/(KSZ^3+n^3)` 의 근거. [PMID 8943305](https://pubmed.ncbi.nlm.nih.gov/8943305/)
20. Koschinsky ML, Côté GP, Gabel B, van der Hoek YY. **Identification of the cysteine residue in apolipoprotein(a) that mediates extracellular coupling with apolipoprotein B-100.** *J Biol Chem.* 1993;268(26):19819-25. — Cys4057. 모델 STEP 2. [PMID 8366118](https://pubmed.ncbi.nlm.nih.gov/8366118/)
21. Gabel BR, Koschinsky ML. **Sequences within apolipoprotein(a) kringle IV types 6-8 bind directly to low-density lipoprotein and mediate noncovalent association of apolipoprotein(a) with apolipoprotein B-100.** *Biochemistry.* 1998;37(21):7892-8. — STEP 1 (비공유 도킹). 무발라플린 표적. [PMID 9601053](https://pubmed.ncbi.nlm.nih.gov/9601053/)
22. Becker L, McLeod RS, Marcovina SM, Yao Z, Koschinsky ML. **Identification of a critical lysine residue in apolipoprotein B-100 that mediates noncovalent interaction with apolipoprotein(a).** *J Biol Chem.* 2001;276(39):36155-62. [PMID 11469568](https://pubmed.ncbi.nlm.nih.gov/11469568/)
23. Youssef A, Clark JR, Koschinsky ML, Boffa MB. **Lipoprotein(a): expanding our knowledge of aortic valve narrowing.** *Trends Cardiovasc Med.* 2021;31(5):305-11. [PMID 32565142](https://pubmed.ncbi.nlm.nih.gov/32565142/)

## 4. 동태학 — 생산이 결정한다 (Kinetics: production, not clearance)

24. Rader DJ, Cain W, Zech LA, Usher D, Brewer HB Jr. **Variation in lipoprotein(a) concentrations among individuals with the same apolipoprotein(a) isoform is determined by the rate of lipoprotein(a) production.** *J Clin Invest.* 1993;91(2):443-7. — 이 모델 전체의 구조적 전제. [PMID 8432853](https://pubmed.ncbi.nlm.nih.gov/8432853/)
25. Rader DJ, Cain W, Ikewaki K, et al. **The inverse association of plasma lipoprotein(a) concentrations with apolipoprotein(a) isoform size is not due to differences in Lp(a) catabolism but to differences in production rate.** *J Clin Invest.* 1994;93(6):2758-63. — FCR이 아이소폼 크기와 무관함. [PMID 8201014](https://pubmed.ncbi.nlm.nih.gov/8201014/)
26. Frischmann ME, Ikewaki K, Trenkwalder E, et al. **In vivo stable-isotope kinetic study suggests intracellular assembly of lipoprotein(a).** *Atherosclerosis.* 2012;225(2):322-7. — 조립 위치에 관한 반대 증거. 모델은 세포외 조립을 채택하되 이 불확실성을 명시. [PMID 23099120](https://pubmed.ncbi.nlm.nih.gov/23099120/)
27. Croyal M, Blanchard V, Ouguerram K, et al. **VLDL (very-low-density lipoprotein)-apo E (apolipoprotein E) may influence Lp(a) (lipoprotein [a]) synthesis or assembly.** *Arterioscler Thromb Vasc Biol.* 2020;40(3):819-29. [PMID 31941383](https://pubmed.ncbi.nlm.nih.gov/31941383/)
28. Chan DC, Watts GF, Coll B, Wasserman SM, Marcovina SM, Barrett PHR. **Lipoprotein(a) particle production as a determinant of plasma lipoprotein(a) concentration across varying apolipoprotein(a) isoform sizes and background cholesterol-lowering therapy.** *J Am Heart Assoc.* 2019;8(7):e011781. [PMID 30897997](https://pubmed.ncbi.nlm.nih.gov/30897997/)
29. Watts GF, Chan DC, Somaratne R, et al. **Controlled study of the effect of proprotein convertase subtilisin-kexin type 9 inhibition with evolocumab on lipoprotein(a) particle kinetics.** *Eur Heart J.* 2018;39(27):2577-85. — 에볼로쿠맙이 Lp(a)를 낮추는 기전이 이화 증가임을 직접 보임. 모델 `KLDLR_LPA` 항의 근거. [PMID 29566128](https://pubmed.ncbi.nlm.nih.gov/29566128/)
30. Reyes-Soffer G, Ginsberg HN, Ramakrishnan R. **The metabolism of lipoprotein(a): a story of dogma, extrapolation and controversy.** *Curr Opin Lipidol.* 2017;28(1):11-15. [PMID 27906712](https://pubmed.ncbi.nlm.nih.gov/27906712/)

## 5. 이화 경로 · 수용체 (Catabolic routes and receptors)

31. Cain WJ, Millar JS, Himebauch AS, et al. **Lipoprotein(a) is cleared from the plasma primarily by the liver in a process mediated by apolipoprotein(a).** *J Lipid Res.* 2005;46(12):2681-91. [PMID 16150825](https://pubmed.ncbi.nlm.nih.gov/16150825/)
32. Romagnuolo R, Scipione CA, Boffa MB, Marcovina SM, Seidah NG, Koschinsky ML. **Lipoprotein(a) catabolism is regulated by proprotein convertase subtilisin/kexin type 9 through the low density lipoprotein receptor.** *J Biol Chem.* 2015;290(18):11649-62. — LDLR 경유 이화의 존재와 그 한계. [PMID 25778403](https://pubmed.ncbi.nlm.nih.gov/25778403/)
33. Yang XP, Amar MJ, Vaisman B, et al. **Scavenger receptor-BI is a receptor for lipoprotein(a).** *J Lipid Res.* 2013;54(9):2450-7. [PMID 23812625](https://pubmed.ncbi.nlm.nih.gov/23812625/)
34. Sharma M, Redpath GM, Williams MJ, McCormick SP. **Recycling of apolipoprotein(a) after PlgRKT-mediated endocytosis of lipoprotein(a).** *Circ Res.* 2017;120(7):1091-1102. [PMID 28003219](https://pubmed.ncbi.nlm.nih.gov/28003219/)
35. Nielsen MB, Çolak Y, Benn M, Nordestgaard BG. **Plasma lipoprotein(a) and risk of aortic valve stenosis: the Copenhagen General Population Study.** *Eur Heart J.* 2019;40(38):3148-56. [PMID 31306481](https://pubmed.ncbi.nlm.nih.gov/31306481/)
36. Nioi P, Sigurdsson A, Thorleifsson G, et al. **Variant ASGR1 associated with a reduced risk of coronary artery disease.** *N Engl J Med.* 2016;374(22):2131-41. — ASGR1 기능소실이 CAD 위험을 낮춤. GalNAc 약물의 진입 경로이기도 함. [PMID 27192541](https://pubmed.ncbi.nlm.nih.gov/27192541/)
37. Kronenberg F, Utermann G, Dieplinger H. **Lipoprotein(a) in renal disease.** *Am J Kidney Dis.* 1996;27(1):1-25. — 신장 이화 경로. 모델 `RENF` 항. [PMID 8546112](https://pubmed.ncbi.nlm.nih.gov/8546112/)
38. Kostner KM, Maurer G, Huber K, et al. **Urinary excretion of apo(a) fragments. Role in apo(a) catabolism.** *Arterioscler Thromb Vasc Biol.* 1996;16(8):905-11. — 모델 `FRAG` 구획. [PMID 8696954](https://pubmed.ncbi.nlm.nih.gov/8696954/)

## 6. 측정 · 아이소폼 편향 · 단위 (Measurement, isoform bias, units)

39. Marcovina SM, Albers JJ, Gabel B, Koschinsky ML, Gaur VP. **Effect of the number of apolipoprotein(a) kringle 4 domains on immunochemical measurements of lipoprotein(a).** *Clin Chem.* 1995;41(2):246-55. — **이 모델 §5(측정) 클러스터의 핵심 근거.** 단일 아이소폼 캘리브레이터를 쓰는 다클론 항체 검사는 큰 아이소폼을 과대, 작은 아이소폼을 과소 보고한다. 모델 `EPIT = (n+10)/(n_cal+10)`. [PMID 7533064](https://pubmed.ncbi.nlm.nih.gov/7533064/)
40. Marcovina SM, Albers JJ, Scanu AM, et al. **Use of a reference material proposed by the International Federation of Clinical Chemistry and Laboratory Medicine to evaluate analytical methods for the determination of plasma lipoprotein(a).** *Clin Chem.* 2000;46(12):1956-67. [PMID 11106328](https://pubmed.ncbi.nlm.nih.gov/11106328/)
41. Marcovina SM, Albers JJ. **Lipoprotein (a) measurements for clinical application.** *J Lipid Res.* 2016;57(4):526-37. — nmol/L 보고 권고의 근거. [PMID 26637279](https://pubmed.ncbi.nlm.nih.gov/26637279/)
42. Tsimikas S, Fazio S, Viney NJ, Xia S, Witztum JL, Marcovina SM. **Relationship of lipoprotein(a) molar concentration and mass according to lipoprotein(a) thresholds and apolipoprotein(a) isoform size.** *J Clin Lipidol.* 2018;12(5):1313-23. — mg/dL과 nmol/L 사이에 단일 환산계수가 존재하지 않음을 실측으로 보인 논문. 모델 RUN 3의 직접 대응. [PMID 30037539](https://pubmed.ncbi.nlm.nih.gov/30037539/)
43. Kronenberg F, Mora S, Stroes ESG, et al. **Lipoprotein(a) in atherosclerotic cardiovascular disease and aortic stenosis: a European Atherosclerosis Society consensus statement.** *Eur Heart J.* 2022;43(39):3925-46. — ≥125 nmol/L / ≥50 mg/dL 역치의 출처. [PMID 36036785](https://pubmed.ncbi.nlm.nih.gov/36036785/)
44. Koschinsky ML, Bajaj A, Boffa MB, et al. **A focused update to the 2019 NLA scientific statement on use of lipoprotein(a) in clinical practice.** *J Clin Lipidol.* 2024;18(3):e308-19. [PMID 38565461](https://pubmed.ncbi.nlm.nih.gov/38565461/)
45. Dahlén GH. **Incidence of Lp(a) lipoprotein among populations.** In: Scanu AM, ed. *Lipoprotein(a).* Academic Press; 1990:151-73. — LDL-C 보정에 쓰이는 0.30 계수의 원출처.
46. Yeang C, Witztum JL, Tsimikas S. **'LDL-C' = LDL-C + Lp(a)-C: implications of achieved ultra-low LDL-C levels in the proprotein convertase subtilisin/kexin type 9 era of potent LDL-C lowering.** *Curr Opin Lipidol.* 2015;26(3):169-78. — 모델 `LDLC_MEAS = LDLC_TRUE + LPAC` 의 근거. [PMID 25943842](https://pubmed.ncbi.nlm.nih.gov/25943842/)
47. Yeang C, Witztum JL, Tsimikas S. **Novel method for quantification of lipoprotein(a)-cholesterol: implications for improving accuracy of LDL-C measurements.** *J Lipid Res.* 2021;62:100053. — Lp(a) 질량 중 콜레스테롤 비율이 30%가 아니라 ~17-25%일 수 있다는 근거. 모델 `FCHOL` 의 불확실성. [PMID 33636162](https://pubmed.ncbi.nlm.nih.gov/33636162/)
48. Ruhaak LR, Cobbaert CM. **Quantifying apolipoprotein(a) in the era of proteoforms and precision medicine.** *Clin Chim Acta.* 2020;511:260-8. [PMID 33096035](https://pubmed.ncbi.nlm.nih.gov/33096035/)

## 7. 병태생리 — 세 개의 작용 팔 (The three effector arms)

### 7a. 죽상형성 · 기질 결합 (Atherogenesis and matrix retention)

49. Nielsen LB, Stender S, Kjeldsen K, Nordestgaard BG. **Specific accumulation of lipoprotein(a) in balloon-injured rabbit aorta in vivo.** *Circ Res.* 1996;78(4):615-26. — Lp(a)의 선택적 내막 축적. 모델 `AVID` 파라미터. [PMID 8635219](https://pubmed.ncbi.nlm.nih.gov/8635219/)
50. van der Hoek YY, Sangrar W, Côté GP, Kastelein JJ, Koschinsky ML. **Binding of recombinant apolipoprotein(a) to extracellular matrix proteins.** *Arterioscler Thromb.* 1994;14(11):1792-8. — 피브로넥틴 결합. [PMID 7947605](https://pubmed.ncbi.nlm.nih.gov/7947605/)
51. Boffa MB, Koschinsky ML. **Oxidized phospholipids as a unifying theory for lipoprotein(a) and cardiovascular disease.** *Nat Rev Cardiol.* 2019;16(5):305-18. [PMID 30675027](https://pubmed.ncbi.nlm.nih.gov/30675027/)
52. Björnson E, Adiels M, Taskinen MR, et al. **Lipoprotein(a) is markedly more atherogenic than LDL: an apolipoprotein B-based genetic analysis.** *J Am Coll Cardiol.* 2024;83(3):385-95. — 입자당 위험이 LDL보다 훨씬 큼. 모델 `AVID`·`WOX` 의 정량 근거. [PMID 38199713](https://pubmed.ncbi.nlm.nih.gov/38199713/)

### 7b. 산화인지질 · 염증 (Oxidized phospholipids and inflammation)

53. Bergmark C, Dewan A, Orsoni A, et al. **A novel function of lipoprotein(a) as a preferential carrier of oxidized phospholipids in human plasma.** *J Lipid Res.* 2008;49(10):2230-9. — ARM 2 의 근거. [PMID 18594118](https://pubmed.ncbi.nlm.nih.gov/18594118/)
54. Tsimikas S, Brilakis ES, Miller ER, et al. **Oxidized phospholipids, Lp(a) lipoprotein, and coronary artery disease.** *N Engl J Med.* 2005;353(1):46-57. [PMID 16000355](https://pubmed.ncbi.nlm.nih.gov/16000355/)
55. Leibundgut G, Scipione C, Yin H, et al. **Determinants of binding of oxidized phospholipids on apolipoprotein(a) and lipoprotein(a).** *J Lipid Res.* 2013;54(10):2815-30. — OxPL이 KIV-10에 공유결합. 입자당 1부위. [PMID 23948545](https://pubmed.ncbi.nlm.nih.gov/23948545/)
56. van der Valk FM, Bekkering S, Kroon J, et al. **Oxidized phospholipids on lipoprotein(a) elicit arterial wall inflammation and an inflammatory monocyte response in humans.** *Circulation.* 2016;134(8):611-24. — 훈련된 면역(trained immunity). apo(a) ASO로 가역. 모델 `MONO` 상태변수. [PMID 27496857](https://pubmed.ncbi.nlm.nih.gov/27496857/)
57. Müller N, Schulte DM, Türk K, et al. **IL-6 blockade by monoclonal antibodies inhibits apolipoprotein (a) expression and lipoprotein (a) synthesis in humans.** *J Lipid Res.* 2015;56(5):1034-42. — LPA 프로모터의 IL-6 반응요소. 토실리주맙으로 Lp(a) 감소. 모델 `FIL6` 항의 직접 근거. [PMID 25713100](https://pubmed.ncbi.nlm.nih.gov/25713100/)
58. Wade DP, Clarke JG, Lindahl GE, et al. **5' control regions of the apolipoprotein(a) gene and members of the related plasminogen gene family.** *Proc Natl Acad Sci USA.* 1993;90(4):1369-73. — 프로모터 구조. [PMID 8433995](https://pubmed.ncbi.nlm.nih.gov/8433995/)
59. Ridker PM, Devalaraja M, Baeres FMM, et al. **IL-6 inhibition with ziltivekimab in patients at high atherosclerotic risk (RESCUE): a double-blind, randomised, placebo-controlled, phase 2 trial.** *Lancet.* 2021;397(10289):2060-9. — hsCRP -92%, Lp(a) 감소. 모델 시나리오 15의 대응 데이터. [PMID 34015342](https://pubmed.ncbi.nlm.nih.gov/34015342/)
60. Ridker PM, Everett BM, Thuren T, et al. **Antiinflammatory therapy with canakinumab for atherosclerotic disease (CANTOS).** *N Engl J Med.* 2017;377(12):1119-31. — 지질 변화 없이 MACE 감소. 모델의 '취약성' 항이 독립적으로 존재해야 하는 이유. [PMID 28845751](https://pubmed.ncbi.nlm.nih.gov/28845751/)

### 7c. 항섬유소용해 (Antifibrinolysis) — 근거의 한계를 포함하여

61. Miles LA, Fless GM, Levin EG, Scanu AM, Plow EF. **A potential basis for the thrombotic risks associated with lipoprotein(a).** *Nature.* 1989;339(6222):301-3. [PMID 2542796](https://pubmed.ncbi.nlm.nih.gov/2542796/)
62. Hancock MA, Boffa MB, Marcovina SM, Nesheim ME, Koschinsky ML. **Inhibition of plasminogen activation by lipoprotein(a): critical domains in apolipoprotein(a) and mechanism of inhibition on fibrin and degraded fibrin surfaces.** *J Biol Chem.* 2003;278(26):23260-9. — KIV-10 의 강한 라이신 결합부위. 입자당 1개. [PMID 12697750](https://pubmed.ncbi.nlm.nih.gov/12697750/)
63. Boffa MB. **Beyond fibrinolysis: the confounding role of Lp(a) in thrombosis.** *Atherosclerosis.* 2022;349:72-81. [PMID 35606082](https://pubmed.ncbi.nlm.nih.gov/35606082/)
64. Kamstrup PR, Tybjærg-Hansen A, Nordestgaard BG. **Genetic evidence that lipoprotein(a) associates with atherosclerotic stenosis rather than venous thrombosis.** *Arterioscler Thromb Vasc Biol.* 2012;32(7):1732-41. — **모델이 ARM 3 의 가중치(WARM3)를 낮게 둔 이유.** 시험관 내 생물학이 강력함에도 정맥혈전증과의 인과 연관이 없음. [PMID 22516069](https://pubmed.ncbi.nlm.nih.gov/22516069/)
65. Helgadottir A, Gretarsdottir S, Thorleifsson G, et al. **Apolipoprotein(a) genetic sequence variants associated with systemic atherosclerosis and coronary atherosclerotic burden but not with venous thromboembolism.** *J Am Coll Cardiol.* 2012;60(8):722-9. [PMID 22898070](https://pubmed.ncbi.nlm.nih.gov/22898070/)

### 7d. 대동맥판막 석회화 (Calcific aortic valve disease)

66. Thanassoulis G, Campbell CY, Owens DS, et al. **Genetic associations with valvular calcification and aortic stenosis.** *N Engl J Med.* 2013;368(6):503-12. — LPA rs10455872 와 대동맥판 석회화의 인과 연관. [PMID 23388002](https://pubmed.ncbi.nlm.nih.gov/23388002/)
67. Bouchareb R, Mahmut A, Nsaibia MJ, et al. **Autotaxin derived from lipoprotein(a) and valve interstitial cells promotes inflammation and mineralization of the aortic valve.** *Circulation.* 2015;132(8):677-90. — 오토탁신 → LysoPA → 골형성 전환. 모델 판막 모듈의 골격. [PMID 26224810](https://pubmed.ncbi.nlm.nih.gov/26224810/)
68. Zheng KH, Tsimikas S, Pawade T, et al. **Lipoprotein(a) and oxidized phospholipids promote valve calcification in patients with aortic stenosis.** *J Am Coll Cardiol.* 2019;73(17):2150-62. — ¹⁸F-NaF PET로 미세석회화 진행 확인. [PMID 31047003](https://pubmed.ncbi.nlm.nih.gov/31047003/)
69. Cowell SJ, Newby DE, Prescott RJ, et al. **A randomized trial of intensive lipid-lowering therapy in calcific aortic stenosis (SALTIRE).** *N Engl J Med.* 2005;352(23):2389-97. — **음성 결과.** 모델의 자기영속 항(KSELFR)이 재현하는 현상. [PMID 15944423](https://pubmed.ncbi.nlm.nih.gov/15944423/)
70. Rossebø AB, Pedersen TR, Boman K, et al. **Intensive lipid lowering with simvastatin and ezetimibe in aortic stenosis (SEAS).** *N Engl J Med.* 2008;359(13):1343-56. — 음성 결과. [PMID 18765433](https://pubmed.ncbi.nlm.nih.gov/18765433/)
71. Chan KL, Teo K, Dumesnil JG, Ni A, Tam J; ASTRONOMER Investigators. **Effect of lipid lowering with rosuvastatin on progression of aortic stenosis (ASTRONOMER).** *Circulation.* 2010;121(2):306-14. — 음성 결과. [PMID 20048204](https://pubmed.ncbi.nlm.nih.gov/20048204/)

## 8. 기존 지질약제의 Lp(a) 효과 (Effect of established lipid drugs)

72. Tsimikas S, Gordts PLSM, Nora C, Yeang C, Witztum JL. **Statin therapy increases lipoprotein(a) levels.** *Eur Heart J.* 2020;41(24):2275-84. — 메타분석. 모델 시나리오 4가 재현해야 하는 대상. [PMID 31111151](https://pubmed.ncbi.nlm.nih.gov/31111151/)
73. de Boer LM, Oorthuys AOJ, Wiegman A, et al. **Statin therapy and lipoprotein(a) levels: a systematic review and meta-analysis.** *Eur J Prev Cardiol.* 2022;29(5):779-92. [PMID 34849691](https://pubmed.ncbi.nlm.nih.gov/34849691/)
74. O'Donoghue ML, Fazio S, Giugliano RP, et al. **Lipoprotein(a), PCSK9 inhibition, and cardiovascular risk: insights from the FOURIER trial.** *Circulation.* 2019;139(12):1483-92. [PMID 30586750](https://pubmed.ncbi.nlm.nih.gov/30586750/)
75. Bittner VA, Szarek M, Aylward PE, et al. **Effect of alirocumab on lipoprotein(a) and cardiovascular risk after acute coronary syndrome (ODYSSEY OUTCOMES).** *J Am Coll Cardiol.* 2020;75(2):133-44. — Lp(a) 감소가 LDL-C와 독립적으로 MACE 감소에 기여. [PMID 31948641](https://pubmed.ncbi.nlm.nih.gov/31948641/)
76. Ray KK, Wright RS, Kallend D, et al. **Two phase 3 trials of inclisiran in patients with elevated LDL cholesterol (ORION-10/-11).** *N Engl J Med.* 2020;382(16):1507-19. — Lp(a) 약 -19 ~ -26%. [PMID 32187462](https://pubmed.ncbi.nlm.nih.gov/32187462/)
77. Awad K, Mikhailidis DP, Katsiki N, Muntner P, Banach M; Lipid and Blood Pressure Meta-Analysis Collaboration Group. **Effect of ezetimibe monotherapy on plasma lipoprotein(a) concentrations in patients with primary hypercholesterolemia: a systematic review and meta-analysis.** *Drugs.* 2018;78(4):453-62. — 에제티미브의 소폭(약 -7%) Lp(a) 감소. **모델이 사후 조정 없이 재현한 예측**(§13 참조). [PMID 29396832](https://pubmed.ncbi.nlm.nih.gov/29396832/)
78. Boden WE, Probstfield JL, Anderson T, et al; AIM-HIGH Investigators. **Niacin in patients with low HDL cholesterol levels receiving intensive statin therapy.** *N Engl J Med.* 2011;365(24):2255-67. — 음성 결과. [PMID 22085343](https://pubmed.ncbi.nlm.nih.gov/22085343/)
79. HPS2-THRIVE Collaborative Group. **Effects of extended-release niacin with laropiprant in high-risk patients.** *N Engl J Med.* 2014;371(3):203-12. — 음성 결과. 모델은 이를 '효력 부족'이 아니라 '절대 감소량 부족'으로 설명. [PMID 25014686](https://pubmed.ncbi.nlm.nih.gov/25014686/)
80. Nicholls SJ, Ditmarsch M, Kastelein JJ, et al. **Lipid lowering effects of the CETP inhibitor obicetrapib in combination with high-intensity statins (ROSE2).** *Nat Med.* 2022;28(8):1672-8. [PMID 35953719](https://pubmed.ncbi.nlm.nih.gov/35953719/)
81. Nicholls SJ, Nelson AJ, Ditmarsch M, et al. **Obicetrapib on top of maximally tolerated lipid-modifying therapies in participants with or at high risk for ASCVD (BROADWAY).** *Nat Med.* 2025;31(2):500-8. [PMID 39653774](https://pubmed.ncbi.nlm.nih.gov/39653774/)
82. Roeseler E, Julius U, Heigl F, et al; Pro(a)LiFe-Study Group. **Lipoprotein apheresis for lipoprotein(a)-associated cardiovascular disease: prospective 5 years of follow-up and apo(a) characterization.** *Arterioscler Thromb Vasc Biol.* 2016;36(9):2019-27. [PMID 27417585](https://pubmed.ncbi.nlm.nih.gov/27417585/)
83. Waldmann E, Parhofer KG. **Lipoprotein apheresis to treat elevated lipoprotein(a).** *J Lipid Res.* 2016;57(10):1751-7. — 급성 60-70% 제거, 간격 평균 30-35%. 모델 시나리오 16. [PMID 26658193](https://pubmed.ncbi.nlm.nih.gov/26658193/)

## 9. RNA 표적 치료제 (RNA-directed therapeutics)

84. Tsimikas S, Viney NJ, Hughes SG, et al. **Antisense therapy targeting apolipoprotein(a): a randomised, double-blind, placebo-controlled phase 1 study.** *Lancet.* 2015;386(10002):1472-83. [PMID 26210642](https://pubmed.ncbi.nlm.nih.gov/26210642/)
85. Viney NJ, van Capelleveen JC, Geary RS, et al. **Antisense oligonucleotides targeting apolipoprotein(a) in people with raised lipoprotein(a): two randomised, double-blind, placebo-controlled, dose-ranging trials.** *Lancet.* 2016;388(10057):2239-53. [PMID 27665230](https://pubmed.ncbi.nlm.nih.gov/27665230/)
86. Tsimikas S, Karwatowska-Prokopczuk E, Gouni-Berthold I, et al. **Lipoprotein(a) reduction in persons with cardiovascular disease (pelacarsen phase 2b).** *N Engl J Med.* 2020;382(3):244-55. — 80 mg 월 1회 → 약 -80%. 모델 시나리오 8의 목표. [PMID 31893580](https://pubmed.ncbi.nlm.nih.gov/31893580/)
87. O'Donoghue ML, Rosenson RS, Gencer B, et al; OCEAN(a)-DOSE Trial Investigators. **Small interfering RNA to reduce lipoprotein(a) in cardiovascular disease.** *N Engl J Med.* 2022;387(20):1855-64. — 올파시란 75 mg q12w → 위약보정 약 -95 ~ -101%. [PMID 36342163](https://pubmed.ncbi.nlm.nih.gov/36342163/)
88. Nissen SE, Wolski K, Watts GF, et al. **Single ascending and multiple-dose trial of zerlasiran, a short interfering RNA targeting lipoprotein(a).** *JAMA.* 2024;331(17):1472-81. [PMID 38583084](https://pubmed.ncbi.nlm.nih.gov/38583084/)
89. Nissen SE, Linnebjerg H, Shen X, et al. **Lepodisiran, an extended-duration short interfering RNA targeting lipoprotein(a): a randomized dose-ascending clinical trial.** *JAMA.* 2023;330(21):2075-83. [PMID 37952254](https://pubmed.ncbi.nlm.nih.gov/37952254/)
90. Nissen SE, Wang Q, Nicholls SJ, et al; ALPACA Investigators. **Lepodisiran in adults with elevated lipoprotein(a).** *N Engl J Med.* 2025;392(19):1892-1902. — 단회 608 mg 로 약 48주간 >90% 감소. [PMID 40162914](https://pubmed.ncbi.nlm.nih.gov/40162914/)
91. Nissen SE, Linnebjerg H, Shen X, et al. **Muvalaplin, an oral small molecule inhibitor of lipoprotein(a) formation: a randomized clinical trial (phase 1).** *JAMA.* 2023;330(11):1042-53. — 조립 억제라는 새로운 표적. 유리 apo(a) 상승 관찰. [PMID 37638695](https://pubmed.ncbi.nlm.nih.gov/37638695/)
92. Nicholls SJ, Ni W, Rhodes GM, et al. **Oral muvalaplin for lowering of lipoprotein(a): a randomized clinical trial (KRAKEN).** *JAMA.* 2025;333(3):222-31. — **intact-Lp(a) 검사로 -85.8%, 전통적 apo(a) 검사로는 더 작은 감소.** 모델 §C(무발라플린 검사 불일치)의 대상. [PMID 39556753](https://pubmed.ncbi.nlm.nih.gov/39556753/)
93. Crooke ST, Baker BF, Crooke RM, Liang XH. **Antisense technology: an overview and prospectus.** *Nat Rev Drug Discov.* 2021;20(6):427-53. — GalNAc-ASGR1 간세포 표적화. 모델 PK 구조. [PMID 33762737](https://pubmed.ncbi.nlm.nih.gov/33762737/)
94. Springer AD, Dowdy SF. **GalNAc-siRNA conjugates: leading the way for delivery of RNAi therapeutics.** *Nucleic Acid Ther.* 2018;28(3):109-18. — RISC 적재의 장반감기. 모델 `SIR_RISC` 구획의 근거. [PMID 29792572](https://pubmed.ncbi.nlm.nih.gov/29792572/)

## 10. 역학 · 위험 정량화 (Epidemiology and risk quantification)

95. Erqou S, Kaptoge S, Perry PL, et al; Emerging Risk Factors Collaboration. **Lipoprotein(a) concentration and the risk of coronary heart disease, stroke, and nonvascular mortality.** *JAMA.* 2009;302(4):412-23. — 1 SD 당 HR 1.13. [PMID 19622820](https://pubmed.ncbi.nlm.nih.gov/19622820/)
96. Burgess S, Ference BA, Staley JR, et al. **Association of LPA variants with risk of coronary disease and the implications for lipoprotein(a)-lowering therapies: a Mendelian randomization analysis.** *JAMA Cardiol.* 2018;3(7):619-27. — **101.5 mg/dL 평생 감소 ≈ 38.67 mg/dL LDL-C 감소.** 모델의 위험 번역 계수의 출처. [PMID 29926099](https://pubmed.ncbi.nlm.nih.gov/29926099/)
97. Lamina C, Kronenberg F; Lp(a)-GWAS-Consortium. **Estimation of the required lipoprotein(a)-lowering therapeutic effect size for reduction in coronary heart disease outcomes: a Mendelian randomization analysis.** *JAMA Cardiol.* 2019;4(6):575-9. — 필요한 절대 감소량. 모델 RUN 4의 논거. [PMID 31017644](https://pubmed.ncbi.nlm.nih.gov/31017644/)
98. Ference BA, Ginsberg HN, Graham I, et al. **Low-density lipoproteins cause atherosclerotic cardiovascular disease. 1. Evidence from genetic, epidemiologic, and clinical studies.** *Eur Heart J.* 2017;38(32):2459-72. — 평생 노출 대 단기 시험의 효과 크기 차이(약 3배). 모델의 두 시간상수 구조의 근거. [PMID 28444290](https://pubmed.ncbi.nlm.nih.gov/28444290/)
99. Cholesterol Treatment Trialists' (CTT) Collaboration. **Efficacy and safety of more intensive lowering of LDL cholesterol: a meta-analysis of data from 170,000 participants in 26 randomised trials.** *Lancet.* 2010;376(9753):1670-81. — 5년 시험의 단위당 효과. [PMID 21067804](https://pubmed.ncbi.nlm.nih.gov/21067804/)
100. Nordestgaard BG, Chapman MJ, Ray K, et al. **Lipoprotein(a) as a cardiovascular risk factor: current status.** *Eur Heart J.* 2010;31(23):2844-53. [PMID 20965889](https://pubmed.ncbi.nlm.nih.gov/20965889/)
101. Willeit P, Ridker PM, Nestel PJ, et al. **Baseline and on-statin treatment lipoprotein(a) levels for prediction of cardiovascular events: individual patient-data meta-analysis of statin outcome trials.** *Lancet.* 2018;392(10155):1311-20. — 잔여 위험. [PMID 30293769](https://pubmed.ncbi.nlm.nih.gov/30293769/)
102. Berman AN, Biery DW, Besser SA, et al. **Lipoprotein(a) and major adverse cardiovascular events in patients with or without baseline atherosclerotic cardiovascular disease.** *J Am Coll Cardiol.* 2024;83(9):873-86. [PMID 38418000](https://pubmed.ncbi.nlm.nih.gov/38418000/)

## 11. 이차성 변화 · 수식 인자 (Secondary causes and modifiers)

103. Wanner C, Rader D, Bartens W, et al. **Elevated plasma lipoprotein(a) in patients with the nephrotic syndrome.** *Ann Intern Med.* 1993;119(4):263-9. — 모델 `FNEPH`. [PMID 8328734](https://pubmed.ncbi.nlm.nih.gov/8328734/)
104. de Bruin TW, van Barlingen H, van Linde-Sibenius Trip M, van Vuurst de Vries AR, Akveld MJ, Erkelens DW. **Lipoprotein(a) and apolipoprotein B plasma concentrations in hypothyroid, euthyroid, and hyperthyroid subjects.** *J Clin Endocrinol Metab.* 1993;76(1):121-6. — 모델 `FTHY`. [PMID 8421075](https://pubmed.ncbi.nlm.nih.gov/8421075/)
105. Anagnostis P, Galanis P, Chatzistergiou V, et al. **The effect of hormone replacement therapy and tibolone on lipoprotein(a) concentrations in postmenopausal women: a systematic review and meta-analysis.** *Maturitas.* 2017;99:27-36. — 모델 `FSEX`. [PMID 28364865](https://pubmed.ncbi.nlm.nih.gov/28364865/)
106. Enkhmaa B, Anuurad E, Berglund L. **Lipoprotein (a): impact by ethnicity and environmental and medical conditions.** *J Lipid Res.* 2016;57(7):1111-25. [PMID 26637279](https://pubmed.ncbi.nlm.nih.gov/26637279/)
107. Burgess S, Davey Smith G. **Mendelian randomization implicates adiposity-related traits in the pathogenesis of lipoprotein(a) elevation.** *Eur J Prev Cardiol.* 2018;25(6):617-9. [PMID 29517303](https://pubmed.ncbi.nlm.nih.gov/29517303/)
108. Norata GD, Ballantyne CM, Catapano AL. **New therapeutic principles in dyslipidaemia: focus on LDL and Lp(a) lowering drugs.** *Eur Heart J.* 2013;34(24):1783-9. [PMID 23509225](https://pubmed.ncbi.nlm.nih.gov/23509225/)

## 12. 지침 · 진행 중인 결과 시험 (Guidelines and ongoing outcome trials)

109. Mach F, Baigent C, Catapano AL, et al. **2019 ESC/EAS Guidelines for the management of dyslipidaemias.** *Eur Heart J.* 2020;41(1):111-88. — 생애 1회 Lp(a) 측정 권고. [PMID 31504418](https://pubmed.ncbi.nlm.nih.gov/31504418/)
110. Reyes-Soffer G, Ginsberg HN, Berglund L, et al; American Heart Association. **Lipoprotein(a): a genetically determined, causal, and prevalent risk factor for atherosclerotic cardiovascular disease — an AHA scientific statement.** *Arterioscler Thromb Vasc Biol.* 2022;42(1):e48-60. [PMID 34647487](https://pubmed.ncbi.nlm.nih.gov/34647487/)
111. Tsimikas S, Moriarty PM, Stroes ES. **Emerging RNA therapeutics to lower blood levels of Lp(a): JACC Focus Seminar 2/4.** *J Am Coll Cardiol.* 2021;77(12):1576-89. [PMID 33766265](https://pubmed.ncbi.nlm.nih.gov/33766265/)
112. Tsimikas S, Karwatowska-Prokopczuk E, Yeang C, et al. **Rationale and design of the Lp(a)HORIZON trial: assessing the effect of pelacarsen on major cardiovascular events in patients with CVD and elevated Lp(a).** *Am Heart J.* 2025;282:9-18. — 등록 기준 Lp(a) ≥70 mg/dL. 모델 RUN 4의 설계 논리와 직접 대응. [PMID 39674305](https://pubmed.ncbi.nlm.nih.gov/39674305/)
113. Nissen SE, Wolski K, Balog C, et al. **OCEAN(a)-Outcomes: rationale and design of a trial of olpasiran in patients with elevated Lp(a) and established ASCVD.** *Am Heart J.* 2024;278:1-10. — 등록 기준 ≥200 nmol/L. [PMID 39179148](https://pubmed.ncbi.nlm.nih.gov/39179148/)

---

## 13. 방정식 ↔ 문헌 대응표 (Equation-to-source map)

| 모델 요소 | 방정식/파라미터 | 근거 |
|---|---|---|
| 생산이 변이를 결정 | `KCAT0` 고정 · `KTL` 역산 | 24, 25, 8 |
| 아이소폼 크기 → 분비 효율 | `SECEFF = KSZ³/(KSZ³+n³)` | 18, 19, 7 |
| 2단계 조립 · 유리 apo(a) | `APOA_FR`, `ASSEM` | 17, 20, 21, 22 |
| 조립이 LDL 기질에 포화 | `LDL_P/(KMLDL+LDL_P)` | 77 (에제티미브가 Lp(a)를 거의 안 바꿈) |
| LDLR 의존 이화가 소수 | `KLDLR_LPA/KCAT0 = 0.24` | 29, 32, 74, 75 |
| 스타틴의 상반된 두 항 | `ESRE` vs `ESTA` | 72, 73 |
| LPA 프로모터 IL-6 반응요소 | `FIL6` | 57, 58, 59 |
| OxPL 입자당 1부위 | `OXPL` | 53, 54, 55 |
| 훈련된 단핵구 | `MONO` | 56 |
| 기질 결합 친화도 2.5배 | `AVID` | 49, 50, 52 |
| 라이신 결합부위 입자당 1개 | `PLGOCC` | 62 |
| 항섬유소용해 팔의 낮은 가중치 | `WARM3 = 0.25` | 64, 65 |
| 오토탁신 → LysoPA → 골형성 | `ATXV`→`LYSOPA`→`VIC_OST` | 67, 68 |
| 자기영속 판막 석회화 | `KSELFR` | 69, 70, 71 |
| 질량 검사의 항체 편향 | `EPIT = (n+10)/(n_cal+10)` | 39, 40, 42 |
| LDL-C 오염 · Dahlén 보정 | `LDLC_MEAS`, `LDLC_CORR` | 45, 46, 47 |
| 평생 노출 대 5년 시험 | `BSLOW`(느림) + `BFAST`(빠름) | 96, 97, 98, 99 |
| GalNAc-ASO / siRNA PK | `PEL_LIV`, `SIR_RISC` | 93, 94, 86, 87 |
| 조립 억제제 | `FMUV` | 21, 91, 92 |

### 모델이 사후 조정 없이 재현한 관측값 (predictions, not fits)

| 관측 | 문헌값 | 모델값 |
|---|---|---|
| 에제티미브의 Lp(a) 효과 | 약 −7% (77) | **−6.6%** |
| PCSK9 억제제 Lp(a) | −25 ~ −30% (74, 75) | **−30.1%** |
| 항IL-6 (비염증 환자) | −16 ~ −25% (59) | **−26.6%** |
| 항IL-6 (류마티스 관절염) | −37% (57) | **−33.5%** |
| hsCRP (RESCUE) | −92% (59) | **−92.5%** |
| 펠라카르센 apoB | 약 −13% (86) | **−11.3%** |

### 이 모델이 답하지 못하는 것 (open questions the model exposes)

1. **KRAKEN 의 검사 불일치 크기.** 모델은 방향(intact 검사가 더 큰 감소를 보고)은
   재현하지만 크기(−85.8% 대 −70%)는 재현하지 못한다. 닫으려면 전통적 apo(a)
   검사의 **유리 apo(a) 에 대한 몰당 반응이 입자 내 apo(a) 대비 약 3.8배**여야
   한다(파라미터 `RFREE`). 이는 측정 가능한 양이며, 알려진 바로는 아직 보고된
   적이 없다. (91, 92)
2. **조립 위치.** 세포외(17) 대 세포내(26) 논쟁은 미해결이며, 모델은 세포외를
   채택했다. 세포내 조립이 맞다면 무발라플린의 유리 apo(a) 상승 예측이 바뀐다.
3. **Lp(a) 질량 중 콜레스테롤 분율.** 고전적 0.30(45) 대 최근 0.17-0.25(47).
   모델 `FCHOL` 이 이 불확실성을 그대로 노출하며, Dahlén 보정의 과대보정 여부가
   여기에 달려 있다.
4. **차단된 조립이 아껴진 apoB 를 혈장 LDL 로 돌려보내는가.** 모델 `FRECY`.
   펠라카르센 시험의 apoB 감소폭(86)은 `FRECY` 가 0 에 가까움을 시사한다.
