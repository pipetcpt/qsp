# 임신성 당뇨병 (GDM) QSP 모델 — 참고문헌
# Gestational Diabetes Mellitus — Annotated Reference List

> **PMID 검증 방법.** 이 목록의 모든 PMID는 NCBI E-utilities
> (`ecitmatch` 인용 매칭 + `esummary`)로 저널·연도·권·페이지·제1저자를
> 대조하여 확인했습니다. 제목과 서지정보는 PubMed 레코드에서 그대로
> 가져온 것이며, 확인되지 않은 인용은 이 목록에 넣지 않았습니다.
>
> **All PMIDs below were machine-verified** against PubMed via E-utilities
> (`ecitmatch` citation matching, then `esummary` to confirm journal, year,
> volume, pages and first author). Titles are taken verbatim from the PubMed
> record. Citations that failed verification were dropped rather than guessed.

총 **85편** · 12개 섹션. 각 항목 끝의 대괄호는 이 문헌이 모델의 어느 부분을
보정하는지를 표시합니다 (예: `[→ CID 가중치]`).

---

## 1. 기준 임상시험 및 대규모 관찰연구 (Landmark trials & cohorts)

모델의 임상 엔드포인트 보정은 거의 전부 이 섹션에서 나옵니다.

1. HAPO Study Cooperative Research Group. **Hyperglycemia and adverse pregnancy outcomes.** N Engl J Med 2008;358:1991-2002.
   <https://pubmed.ncbi.nlm.nih.gov/18463375/>
   *25,505명 대상 맹검 관찰연구. 공복혈당 7분위에 걸쳐 LGA 5.3% → 26.3%, 제대 C-peptide >90백분위 3.7% → 32.4%. 역치가 아니라 연속적 기울기임을 확립.*
   `[→ P(LGA) 보정의 1차 기준, 시나리오 S9]`

2. Crowther CA, et al. **Effect of treatment of gestational diabetes mellitus on pregnancy outcomes.** N Engl J Med 2005;352:2477-86. (ACHOIS)
   <https://pubmed.ncbi.nlm.nih.gov/15951574/>
   *치료가 중대한 주산기 복합 결과를 4% → 1%로 감소.* `[→ 치료효과 크기]`

3. Landon MB, et al. **A multicenter, randomized trial of treatment for mild gestational diabetes.** N Engl J Med 2009;361:1339-48. (MFMU)
   <https://pubmed.ncbi.nlm.nih.gov/19797280/>
   *경증 GDM에서 LGA 14.5→7.1%, 거대아 14.3→5.9%, 견부난산 4.0→1.5%, 전자간증 13.6→8.6%, 제왕절개 33.8→26.9%.*
   `[→ PSD0/KSD, PCS0/KCS, PPE0/KPEG/KPEB 를 이 5개 쌍에서 직접 해석적으로 풀었음]`

4. Langer O, et al. **A comparison of glyburide and insulin in women with gestational diabetes mellitus.** N Engl J Med 2000;343:1134-8.
   <https://pubmed.ncbi.nlm.nih.gov/11036118/>
   *글리부라이드를 1차 대안으로 널리 쓰이게 만든 원 연구. 후속 연구들이 이 결론을 약화시킴(§5 참조).*

5. Rowan JA, et al. **Metformin versus insulin for the treatment of gestational diabetes.** N Engl J Med 2008;358:2003-15. (MiG)
   <https://pubmed.ncbi.nlm.nih.gov/18463376/>
   *복합 신생아 결과 동등. 핵심 수치: 메트포민군의 **46%가 보조 인슐린 필요**.*
   `[→ 시나리오 S8은 이 46%를 입력이 아니라 예측으로 재현하려는 시도]`

6. Sénat MV, et al. **Effect of Glyburide vs Subcutaneous Insulin on Perinatal Complications Among Women With Gestational Diabetes: A Randomized Clinical Trial.** JAMA 2018;319:1773-1780.
   <https://pubmed.ncbi.nlm.nih.gov/29715355/>
   *비열등성 입증 실패. 신생아 저혈당이 초과 위험의 주된 성분.* `[→ 글리부라이드 태아 SUR1 항]`

7. Balsells M, et al. **Glibenclamide, metformin, and insulin for the treatment of gestational diabetes: a systematic review and meta-analysis.** BMJ 2015;350:h102.
   <https://pubmed.ncbi.nlm.nih.gov/25609400/>
   *글리부라이드는 인슐린 대비 출생체중 +109 g, 거대아 RR 2.62, 신생아 저혈당 RR 2.04.*

8. Feig DS, et al. **Metformin in women with type 2 diabetes in pregnancy (MiTy): a multicentre, international, randomised, placebo-controlled trial.** Lancet Diabetes Endocrinol 2020;8:834-844.
   <https://pubmed.ncbi.nlm.nih.gov/32946820/>
   *인슐린에 메트포민 추가 시 인슐린 요구량·체중증가·거대아 감소, 그러나 SGA 증가 신호.*

9. Feig DS, et al. **Continuous glucose monitoring in pregnant women with type 1 diabetes (CONCEPTT): a multicentre international randomised controlled trial.** Lancet 2017;390:2347-2359.
   <https://pubmed.ncbi.nlm.nih.gov/28923465/>
   `[→ TIR 목표 및 CGM 지표 정의]`

10. Hillier TA, et al. **A Pragmatic, Randomized Clinical Trial of Gestational Diabetes Screening.** N Engl J Med 2021;384:895-904. (ScreenR2GDM)
    <https://pubmed.ncbi.nlm.nih.gov/33704936/>
    *1단계(IADPSG) 접근은 진단율을 두 배로 늘렸으나 주산기 결과는 개선하지 않았음 — 진단 절단점 문제를 모델링해야 하는 이유.*

11. Stewart ZA, et al. **Closed-Loop Insulin Delivery during Pregnancy in Women with Type 1 Diabetes.** N Engl J Med 2016;375:644-54.
    <https://pubmed.ncbi.nlm.nih.gov/27532830/>

12. Donovan LE, et al. **Closed-Loop Insulin Delivery in Type 1 Diabetes in Pregnancy: The CIRCUIT Randomized Clinical Trial.** JAMA 2025;334:2176-2185.
    <https://pubmed.ncbi.nlm.nih.gov/41134589/>

---

## 2. 인슐린 저항성의 기전 (Mechanisms of gestational insulin resistance)

모델의 `CID`(대항인슐린 지수) → `SIREL` 블록 근거.

13. Barbour LA, et al. **Cellular mechanisms for insulin resistance in normal pregnancy and gestational diabetes.** Diabetes Care 2007;30 Suppl 2:S112-9.
    <https://pubmed.ncbi.nlm.nih.gov/17596458/>
    *수용체 후 결함(IRS-1 단백 감소, p85 증가, GLUT4 전위 감소)이 핵심이라는 종합 리뷰.*
    `[→ SIREL 을 수용체 후 이득으로 모델링한 근거]`

14. Kirwan JP, et al. **TNF-alpha is a predictor of insulin resistance in human pregnancy.** Diabetes 2002;51:2207-13.
    <https://pubmed.ncbi.nlm.nih.gov/12086951/>
    *TNF-α가 hPL·코르티솔·에스트라디올·프로락틴·렙틴보다 인슐린 감수성과 더 강하게 역상관.*
    `[→ STNF 가 단위변화당 최대 가중치를 갖는 이유]`

15. Catalano PM, et al. **Longitudinal changes in glucose metabolism during pregnancy in obese women with normal glucose tolerance and gestational diabetes mellitus.** Am J Obstet Gynecol 1999;180:903-16.
    <https://pubmed.ncbi.nlm.nih.gov/10203659/>
    `[→ 임신 중 SI 50-60% 감소 목표치]`

16. Catalano PM, et al. **Carbohydrate metabolism during pregnancy in control subjects and women with gestational diabetes.** Am J Physiol 1993;264:E60-7.
    <https://pubmed.ncbi.nlm.nih.gov/8430789/>
    `[→ EGP·기저 포도당 대사율의 임신 중 변화]`

17. Catalano PM, et al. **Longitudinal changes in insulin release and insulin resistance in nonobese pregnant women.** Am J Obstet Gynecol 1991;165:1667-72.
    <https://pubmed.ncbi.nlm.nih.gov/1750458/>

18. Ryan EA, et al. **Insulin action during pregnancy. Studies with the euglycemic clamp technique.** Diabetes 1985;34:380-9.
    <https://pubmed.ncbi.nlm.nih.gov/3882502/>
    *클램프로 측정한 임신 중 인슐린 저항성의 고전 연구.*

19. Barbour LA, et al. **Human placental growth hormone causes severe insulin resistance in transgenic mice.** Am J Obstet Gynecol 2002;186:512-7.
    <https://pubmed.ncbi.nlm.nih.gov/11904616/>
    `[→ GH-V 를 hPL 과 함께 태반 구동에 포함시킨 근거]`

20. Barbour LA, et al. **Human placental growth hormone increases expression of the p85 regulatory unit of phosphatidylinositol 3-kinase and triggers severe insulin resistance in skeletal muscle.** Endocrinology 2004;145:1144-50.
    <https://pubmed.ncbi.nlm.nih.gov/14633976/>

21. Masuyama H, Hiramatsu Y. **Potential role of estradiol and progesterone in insulin resistance through constitutive androstane receptor.** J Mol Endocrinol 2011;47:229-39.
    <https://pubmed.ncbi.nlm.nih.gov/21768169/>
    `[→ SPROG 가중치]`

22. González C, et al. **Regulation of insulin receptor substrate-1 in the liver, skeletal muscle and adipose tissue of rats throughout pregnancy.** Gynecol Endocrinol 2003;17:187-97.
    <https://pubmed.ncbi.nlm.nih.gov/12857426/>
    `[→ 간 vs 말초 저항성 분리(HFR = 0.8)]`

---

## 3. 베타세포 적응과 그 실패 (β-cell adaptation and its failure)

**이 섹션이 모델의 중심 가설(BCAP)의 근거입니다.**

23. Butler AE, et al. **Adaptive changes in pancreatic beta cell fractional area and beta cell turnover in human pregnancy.** Diabetologia 2010;53:2167-76.
    <https://pubmed.ncbi.nlm.nih.gov/20523966/>
    *사람에서 베타세포 질량 증가는 약 1.4배에 불과 — 설치류의 2-3배와 다름. 사람의 적응은 **질량보다 기능** 중심.*
    `[→ BCM 상한 2.6 및 KSENS(곡선 좌측 이동)로 기능 성분을 분리한 이유]`

24. Kim H, et al. **Serotonin regulates pancreatic beta cell mass during pregnancy.** Nat Med 2010;16:804-8.
    <https://pubmed.ncbi.nlm.nih.gov/20581837/>
    *락토겐 → Tph1 → 세로토닌 → HTR2B(증식)/HTR3A(분비 역치 저하).*
    `[→ ADAPT · KSENS 경로의 분자 근거]`

25. Schraenen A, et al. **Placental lactogens induce serotonin biosynthesis in a subset of mouse beta cells during pregnancy.** Diabetologia 2010;53:2589-99.
    <https://pubmed.ncbi.nlm.nih.gov/20938637/>

26. Vasavada RC, et al. **Targeted expression of placental lactogen in the beta cells of transgenic mice results in beta cell proliferation, islet mass augmentation, and hypoglycemia.** J Biol Chem 2000;275:15399-406.
    <https://pubmed.ncbi.nlm.nih.gov/10809775/>
    `[→ LACTDR = HPL/(HPL+KMH) 구동항]`

27. Billestrup N, Nielsen JH. **The stimulatory effect of growth hormone, prolactin, and placental lactogen on beta-cell proliferation is not mediated by insulin-like growth factor-I.** Endocrinology 1991;129:883-8.
    <https://pubmed.ncbi.nlm.nih.gov/1677331/>

28. Qiao L, et al. **Adiponectin Promotes Maternal β-Cell Expansion Through Placental Lactogen Expression.** Diabetes 2021;70:132-142.
    <https://pubmed.ncbi.nlm.nih.gov/33087456/>
    *아디포넥틴 → 태반 락토겐 → 베타세포 확장. 저아디포넥틴이 단순한 표지가 아니라 인과적일 수 있음.*
    `[→ ADIPO 를 상태변수로 둔 이유; 현재 모델은 SI 경로만 연결하고 이 β-세포 경로는 미구현(한계로 명시)]`

29. Buchanan TA, Xiang AH. **Gestational diabetes mellitus.** J Clin Invest 2005;115:485-91.
    <https://pubmed.ncbi.nlm.nih.gov/15765129/>
    *"GDM은 만성적 베타세포 결손이 임신이라는 부하검사에서 드러난 것"이라는 프레임을 세운 리뷰.*
    `[→ 모델의 중심 논지 전체]`

30. Powe CE, et al. **Heterogeneous Contribution of Insulin Sensitivity and Secretion Defects to Gestational Diabetes Mellitus.** Diabetes Care 2016;39:1052-5.
    <https://pubmed.ncbi.nlm.nih.gov/27208340/>
    *GDM을 인슐린저항 우세형 / 분비결핍형 / 혼합형으로 분해. 아형마다 표현형과 치료반응이 다름.*
    `[→ 이 모델에서 아형은 지정되지 않고 BCAP×BMI 로부터 창발함 — 탭 1의 아형 좌표]`

31. Xiang AH, et al. **Longitudinal changes in insulin sensitivity and beta cell function between women with and without a history of gestational diabetes mellitus.** Diabetologia 2013;56:2753-60.
    <https://pubmed.ncbi.nlm.nih.gov/24030069/>
    `[→ 산후 베타세포 감소(KDECL)와 disposition index 궤적]`

---

## 4. 태반 수송과 태아 (Placental transport & the fetal compartment)

32. Freinkel N. **Banting Lecture 1980. Of pregnancy and progeny.** Diabetes 1980;29:1023-35.
    <https://pubmed.ncbi.nlm.nih.gov/7002669/>
    *"연료 매개 기형발생(fuel-mediated teratogenesis)" — 포도당만이 아니라 아미노산·지질도 태아 프로그래밍의 연료라는 명제.*
    `[→ BFAT(모체 FFA → 태아 지방) 항의 근거]`

33. Freinkel N, et al. **Diabetic embryopathy and fuel-mediated organ teratogenesis: lessons from animal models.** Horm Metab Res 1988;20:463-75.
    <https://pubmed.ncbi.nlm.nih.gov/3053387/>

34. Hahn T, et al. **Sustained hyperglycemia in vitro down-regulates the GLUT1 glucose transport system of cultured human term placental trophoblast: a mechanism to protect fetal development?** FASEB J 1998;12:1221-31.
    <https://pubmed.ncbi.nlm.nih.gov/9737725/>
    *고혈당이 GLUT1을 하향조절하는 보호 기전 — 현재 모델은 이 적응을 **구현하지 않았고**, 따라서 극단적 고혈당에서 태아 노출을 과대추정할 수 있음.*
    `[→ 한계로 명시됨]`

35. Takata K, et al. **Immunolocalization of glucose transporter GLUT1 in the rat placental barrier.** Cell Tissue Res 1994;276:411-8.
    <https://pubmed.ncbi.nlm.nih.gov/8062336/>

36. Barta E, Drugan A. **A theoretical model of glucose transport suggests symmetric GLUT1 characteristics at placental membranes.** J Membr Biol 2014;247:685-94.
    <https://pubmed.ncbi.nlm.nih.gov/24894722/>
    `[→ KTRF 를 대칭 양방향 확산으로 모델링한 근거]`

37. Catalano PM, et al. **Increased fetal adiposity: a very sensitive marker of abnormal in utero development.** Am J Obstet Gynecol 2003;189:1698-704.
    <https://pubmed.ncbi.nlm.nih.gov/14710101/>
    *출생체중이 같아도 GDM 신생아의 체지방은 더 많음 — 과성장이 비대칭이라는 직접 증거.*
    `[→ AFAT(0.55) ≫ AIGF(0.20) 의 핵심 근거; 정상 만삭 지방 13.3%]`

38. Catalano PM, et al. **Fetuses of obese mothers develop insulin resistance in utero.** Diabetes Care 2009;32:1076-80.
    <https://pubmed.ncbi.nlm.nih.gov/19460915/>

39. Schaefer-Graf UM, et al. **Differences in the implications of maternal lipids on fetal metabolism and growth between gestational diabetes mellitus and control pregnancies.** Diabet Med 2011;28:1053-9.
    <https://pubmed.ncbi.nlm.nih.gov/21658120/>
    `[→ BFAT 계수]`

40. Kemball ML, et al. **Neonatal hypoglycaemia in infants of diabetic mothers given sulphonylurea drugs in pregnancy.** Arch Dis Child 1970;45:696-701.
    <https://pubmed.ncbi.nlm.nih.gov/5477685/>
    *설폰요소제의 태아 노출이 신생아 저혈당을 일으킨다는 50년 전 관찰 — 2018년 Sénat 결과의 예고.*

41. Haworth JC, et al. **Prognosis of infants of diabetic mothers in relation to neonatal hypoglycaemia.** Dev Med Child Neurol 1976;18:471-9.
    <https://pubmed.ncbi.nlm.nih.gov/955311/>

42. Agrawal RK, et al. **Neonatal hypoglycaemia in infants of diabetic mothers.** J Paediatr Child Health 2000;36:354-6.
    <https://pubmed.ncbi.nlm.nih.gov/10940170/>
    `[→ PNH0 / PNHMAX 범위]`

---

## 5. 약물 PK/PD 및 태반 통과 (Drug PK/PD and placental passage)

**세 약물이 실제로 다른 유일한 축.**

### 메트포민 (crosses freely)

43. Eyal S, et al. **Pharmacokinetics of metformin during pregnancy.** Drug Metab Dispos 2010;38:833-40.
    <https://pubmed.ncbi.nlm.nih.gov/20118196/>
    `[→ CLM · CLMP(임신 중 신클리어런스 +25%)]`

44. Charles B, et al. **Population pharmacokinetics of metformin in late pregnancy.** Ther Drug Monit 2006;28:67-72.
    <https://pubmed.ncbi.nlm.nih.gov/16418696/>
    `[→ VMC/VMP/QM 2구획 구조]`

45. Abduljalil K, et al. **Prediction of Maternal and Fetal Acyclovir, Emtricitabine, Lamivudine, and Metformin Concentrations during Pregnancy Using a Physiologically Based Pharmacokinetic Modeling Approach.** Clin Pharmacokinet 2022;61:725-748.
    <https://pubmed.ncbi.nlm.nih.gov/35067869/>
    *모체–태아 PBPK. 이 QSP 모델의 KPLM/CLMF 가 재현해야 하는 대상.*

46. Tiley JB, et al. **Comparison of Metformin PBPK Models Incorporating Placental Transfer to Predict Fetal and Maternal Exposure.** CPT Pharmacometrics Syst Pharmacol 2026;15:e70136.
    <https://pubmed.ncbi.nlm.nih.gov/41289433/>
    `[→ umbilical:maternal ≈ 1.0 목표; 모델은 KPLM/(KPLM+CLMF) = 0.94]`

47. Gu X, et al. **Transplacental transfer of metformin and interaction of metformin with the uptake transporters of placental trophoblast cells.** Eur J Pharm Sci 2025;212:107161.
    <https://pubmed.ncbi.nlm.nih.gov/40494429/>
    `[→ OCT3/PMAT 매개 유입]`

48. Sheng B, et al. **Short-term neonatal outcomes in women with gestational diabetes treated using metformin versus insulin: a systematic review and meta-analysis of randomized controlled trials.** Acta Diabetol 2023;60:595-608.
    <https://pubmed.ncbi.nlm.nih.gov/36593391/>

### 글리부라이드 (crosses; BCRP-effluxed)

49. Hebert MF, et al. **Are we optimizing gestational diabetes treatment with glyburide? The pharmacologic basis for better clinical practice.** Clin Pharmacol Ther 2009;85:607-14.
    <https://pubmed.ncbi.nlm.nih.gov/19295505/>
    *임신 중 글리부라이드 클리어런스가 약 2배 — 용량 부족과 태아 노출이 동시에 문제.*
    `[→ CLGP = 1.0 (2배)]`

50. Gedeon C, et al. **Breast cancer resistance protein: mediating the trans-placental transfer of glyburide across the human placenta.** Placenta 2008;29:39-43.
    <https://pubmed.ncbi.nlm.nih.gov/17923155/>
    `[→ KPLG 를 BCRP 유출을 감한 순수송으로 정의]`

51. Kraemer J, et al. **Perfusion studies of glyburide transfer across the human placenta: implications for fetal safety.** Am J Obstet Gynecol 2006;195:270-4.
    <https://pubmed.ncbi.nlm.nih.gov/16579925/>

52. Elliott BD, et al. **Insignificant transfer of glyburide occurs across the human placenta.** Am J Obstet Gynecol 1991;165:807-12.
    <https://pubmed.ncbi.nlm.nih.gov/1951536/>
    *반대 결론의 원 논문 — 이후 관류·제대혈 연구들이 뒤집었습니다. 문헌이 서로 충돌할 때 모델이 어느 쪽을 택했는지 밝히는 것이 정직합니다: 이 모델은 통과를 인정하고 cord:maternal ≈ 0.7 로 설정했습니다.*

53. Yu DQ, et al. **Glycemic control and neonatal outcomes in women with gestational diabetes mellitus treated using glyburide, metformin, or insulin: a pairwise and network meta-analysis.** BMC Endocr Disord 2021;21:199.
    <https://pubmed.ncbi.nlm.nih.gov/34641848/>

### 인슐린 (does not cross)

54. Athanasiadou KI, et al. **Safety and efficacy of insulin detemir versus NPH in the treatment of diabetes during pregnancy: Systematic review and meta-analysis of randomized controlled trials.** Diabetes Res Clin Pract 2022;190:110020.
    <https://pubmed.ncbi.nlm.nih.gov/35878788/>
    `[→ KAB (기저 유사체 흡수)]`

55. Koren R, Toledano Y, Hod M. **The use of insulin detemir during pregnancy: a safety evaluation.** Expert Opin Drug Saf 2015;14:593-9.
    <https://pubmed.ncbi.nlm.nih.gov/25731934/>

---

## 6. 장기 결과 — 모체 (Long-term maternal outcomes)

56. Bellamy L, et al. **Type 2 diabetes mellitus after gestational diabetes: a systematic review and meta-analysis.** Lancet 2009;373:1773-9.
    <https://pubmed.ncbi.nlm.nih.gov/19465232/>
    *상대위험 7.43.* `[→ KDI = 3.8 을 이 RR에서 역산; DI/DIREF = 0.45 에서 RR ≈ 8]`

57. Vounzoulaki E, et al. **Progression to type 2 diabetes in women with a known history of gestational diabetes: systematic review and meta-analysis.** BMJ 2020;369:m1361.
    <https://pubmed.ncbi.nlm.nih.gov/32404325/>
    `[→ 5년 누적발생률]`

58. Gunderson EP, et al. **Lactation and Progression to Type 2 Diabetes Mellitus After Gestational Diabetes Mellitus: A Prospective Cohort Study.** Ann Intern Med 2015;163:889-98.
    <https://pubmed.ncbi.nlm.nih.gov/26595611/>
    `[→ LACT / KLACT = 0.45]`

59. Schwartz N, et al. **The prevalence of gestational diabetes mellitus recurrence--effect of ethnicity and parity: a metaanalysis.** Am J Obstet Gynecol 2015;213:310-7.
    <https://pubmed.ncbi.nlm.nih.gov/25757637/>

60. Pan Y, et al. **Gestational diabetes mellitus recurrence rate and risk factors: a systematic review and meta-analysis.** Diabetes Res Clin Pract 2025;230:112949.
    <https://pubmed.ncbi.nlm.nih.gov/41130422/>

61. Damm P, et al. **Gestational diabetes mellitus and long-term consequences for mother and offspring: a view from Denmark.** Diabetologia 2016;59:1396-1399.
    <https://pubmed.ncbi.nlm.nih.gov/27174368/>

62. Nouhjah S, et al. **Postpartum screening practices, progression to abnormal glucose tolerance and its related risk factors in Asian women with a known history of gestational diabetes: A systematic review and meta-analysis.** Diabetes Metab Syndr 2017;11 Suppl 2:S703-S712.
    <https://pubmed.ncbi.nlm.nih.gov/28571777/>

---

## 7. 장기 결과 — 자녀 (Offspring outcomes & developmental programming)

63. Lowe WL Jr, et al. **Association of Gestational Diabetes With Maternal Disorders of Glucose Metabolism and Childhood Adiposity.** JAMA 2018;320:1005-1016. (HAPO Follow-Up Study)
    <https://pubmed.ncbi.nlm.nih.gov/30208453/>

64. Hillier TA, et al. **Childhood obesity and metabolic imprinting: the ongoing effects of maternal hyperglycemia.** Diabetes Care 2007;30:2287-92.
    <https://pubmed.ncbi.nlm.nih.gov/17519427/>

65. Dabelea D, Crume T. **Maternal environment and the transgenerational cycle of obesity and diabetes.** Diabetes 2011;60:1849-55.
    <https://pubmed.ncbi.nlm.nih.gov/21709280/>

66. Rowan JA, et al. **Metformin in gestational diabetes: the offspring follow-up (MiG TOFU): body composition at 2 years of age.** Diabetes Care 2011;34:2279-84.
    <https://pubmed.ncbi.nlm.nih.gov/21949222/>
    *메트포민 노출아는 2세에 피하지방이 더 많았음.*
    `[→ 모델은 이 결과를 예측하지 못하며, dxdt_FATF 의 태아 메트포민 항을 의도적으로 0으로 두고 비활성 상태로 남겨 그 공백을 드러냄]`

67. Rowan JA, et al. **Metformin in Gestational Diabetes The Offspring Follow Up (MiGTOFU): Associations between maternal characteristics and size and adiposity of boys and girls at nine years.** Aust N Z J Obstet Gynaecol 2023;63:825-828.
    <https://pubmed.ncbi.nlm.nih.gov/37469163/>

68. Plagemann A, Harder T. **Fuel-mediated teratogenesis and breastfeeding.** Diabetes Care 2011;34:779-81.
    <https://pubmed.ncbi.nlm.nih.gov/21357365/>

---

## 8. 진단·선별·지침 (Diagnosis, screening, guidelines)

69. International Association of Diabetes and Pregnancy Study Groups Consensus Panel. **International association of diabetes and pregnancy study groups recommendations on the diagnosis and classification of hyperglycemia in pregnancy.** Diabetes Care 2010;33:676-82.
    <https://pubmed.ncbi.nlm.nih.gov/20190296/>
    *IADPSG 기준: 공복 5.1 / 1시간 10.0 / 2시간 8.5 mmol/L (92 / 180 / 153 mg/dL).*

70. American Diabetes Association. **13. Management of Diabetes in Pregnancy: Standards of Medical Care in Diabetes-2018.** Diabetes Care 2018;41:S137-S143.
    <https://pubmed.ncbi.nlm.nih.gov/29222384/>
    `[→ 공복 <95, 1시간 <140 mg/dL 치료 목표]`

71. **ACOG Practice Bulletin No. 190: Gestational Diabetes Mellitus.** Obstet Gynecol 2018;131:e49-e64.
    <https://pubmed.ncbi.nlm.nih.gov/29370047/>

72. Brady M, et al. **One-Step Compared With Two-Step Gestational Diabetes Screening and Pregnancy Outcomes: A Systematic Review and Meta-analysis.** Obstet Gynecol 2022;140:712-723.
    <https://pubmed.ncbi.nlm.nih.gov/36201772/>

---

## 9. 바이오마커 (Biomarkers)

73. Ghosh P, et al. **Plasma Glycated CD59, a Novel Biomarker for Detection of Pregnancy-Induced Glucose Intolerance.** Diabetes Care 2017;40:981-984.
    <https://pubmed.ncbi.nlm.nih.gov/28450368/>

74. Williams MA, et al. **Plasma adiponectin concentrations in early pregnancy and subsequent risk of gestational diabetes mellitus.** J Clin Endocrinol Metab 2004;89:2306-11.
    <https://pubmed.ncbi.nlm.nih.gov/15126557/>
    `[→ ADIPO 를 1분기 예측 인자로 둔 근거]`

75. Ye Y, et al. **Adiponectin, leptin, and leptin/adiponectin ratio with risk of gestational diabetes mellitus: A prospective nested case-control study among Chinese women.** Diabetes Res Clin Pract 2022;191:110039.
    <https://pubmed.ncbi.nlm.nih.gov/35985429/>

76. Durnwald C, et al. **Continuous Glucose Monitoring Profiles in Pregnancies With and Without Gestational Diabetes Mellitus.** Diabetes Care 2024;47:1333-1341.
    <https://pubmed.ncbi.nlm.nih.gov/38701400/>
    `[→ CGM 탭의 TIR·평균혈당 기준값]`

77. García-Moreno RM, et al. **Efficacy of continuous glucose monitoring on maternal and neonatal outcomes in gestational diabetes mellitus: a systematic review and meta-analysis of randomized clinical trials.** Diabet Med 2022;39:e14703.
    <https://pubmed.ncbi.nlm.nih.gov/34564868/>

78. Benhalima K, et al. **Application of continuous glucose monitoring and automated insulin delivery technologies for pregnant women with type 1, type 2, or gestational diabetes: an international consensus statement.** Lancet Diabetes Endocrinol 2026;14:157-177.
    <https://pubmed.ncbi.nlm.nih.gov/41421368/>

---

## 10. 예방·생활습관 (Prevention: lifestyle, exercise, supplements)

79. Davenport MH, et al. **Prenatal exercise for the prevention of gestational diabetes mellitus and hypertensive disorders of pregnancy: a systematic review and meta-analysis.** Br J Sports Med 2018;52:1367-1375.
    <https://pubmed.ncbi.nlm.nih.gov/30337463/>
    `[→ EXEFF = 0.20 (주 150분 이상)]`

80. Tsironikos GI, et al. **Effectiveness of exercise intervention during pregnancy on high-risk women for gestational diabetes mellitus prevention: A meta-analysis of published RCTs.** PLoS One 2022;17:e0272711.
    <https://pubmed.ncbi.nlm.nih.gov/35930592/>

81. Mashayekh-Amiri S, et al. **Myo-inositol supplementation for prevention of gestational diabetes mellitus in overweight and obese pregnant women: a systematic review and meta-analysis.** Diabetol Metab Syndr 2022;14:93.
    <https://pubmed.ncbi.nlm.nih.gov/35794663/>

---

## 11. 유전학 (Genetics of GDM susceptibility)

82. Shan D, et al. **MTNR1B rs1387153 Polymorphism and Risk of Gestational Diabetes Mellitus: Meta-Analysis and Trial Sequential Analysis.** Public Health Genomics 2023;26:201-211.
    <https://pubmed.ncbi.nlm.nih.gov/37980891/>
    `[→ BCAP 를 유전적 요소를 포함하는 환자 수준 파라미터로 둔 근거]`

83. Zhang Y, et al. **MTNR1B gene variations and high pre-pregnancy BMI increase gestational diabetes mellitus risk in Chinese women.** Gene 2024;894:148023.
    <https://pubmed.ncbi.nlm.nih.gov/38007162/>
    *BCAP(유전) × BMI(저항성) 상호작용의 임상적 대응 — 모델 탭 1의 아형 좌표가 바로 이 2차원 평면입니다.*

84. Huang LT, et al. **Adiponectin gene polymorphisms and risk of gestational diabetes mellitus: A meta-analysis.** World J Clin Cases 2019;7:572-584.
    <https://pubmed.ncbi.nlm.nih.gov/30863757/>

---

## 12. 전자간증과의 공유 병태생리 (Shared biology with preeclampsia)

85. Elgazzaz M, Lazartigues E. **Implications of pregnancy on cardiometabolic disease risk: preeclampsia and gestational diabetes.** Am J Physiol Cell Physiol 2024;327:C646-C660.
    <https://pubmed.ncbi.nlm.nih.gov/39010840/>
    *GDM과 전자간증은 상류 생물학(내피 기능장애·염증)을 공유합니다. 이 모델의 지도에는 sFlt-1/PlGF가 그려져 있으나 **미분방정식으로는 구현되지 않았고**, 전자간증은 회귀식으로만 다룹니다 — 한계로 명시.*

86. Aziz F, et al. **Gestational diabetes mellitus, hypertension, and dyslipidemia as the risk factors of preeclampsia.** Sci Rep 2024;14:6182.
    <https://pubmed.ncbi.nlm.nih.gov/38486097/>

87. Verlohren S, Dröge LA. **Clinical interpretation and implementation of the sFlt-1/PlGF ratio in the prediction, diagnosis and management of preeclampsia.** Pregnancy Hypertens 2022;27:42-50.
    <https://pubmed.ncbi.nlm.nih.gov/34915395/>

---

## 문헌이 서로 충돌하는 지점 (Where the literature disagrees)

모델은 다음 세 지점에서 한쪽을 **선택**했습니다. 선택을 밝히지 않으면
파라미터가 사실처럼 보이기 때문에 여기 적어둡니다.

| 쟁점 | 상충하는 근거 | 이 모델의 선택 |
|---|---|---|
| 글리부라이드 태반 통과 | Elliott 1991 (#52) "무의미한 통과" vs Kraemer 2006 (#51)·Gedeon 2008 (#50) 관류/BCRP 연구 | 통과 인정, cord:maternal ≈ 0.70 |
| 1단계 vs 2단계 선별 | IADPSG (#69) 권고 vs Hillier NEJM 2021 (#10) 주산기 결과 개선 없음 | 두 기준 모두 관찰값으로 출력, 어느 쪽도 정답으로 두지 않음 |
| 메트포민의 자녀 장기 영향 | MiG 단기 동등 (#5) vs MiG-TOFU 지방 증가 (#66, #67) | 태아 노출은 계산하되 성장 효과는 **0으로 고정** — 기전이 없는데 효과를 넣으면 조작이 됨 |

---

## 인용 시 유의 (Caveat)

본 문헌 목록은 모델 파라미터의 **출처**를 추적할 수 있게 하는 것이 목적입니다.
개별 파라미터가 해당 논문에서 그대로 보고된 값인 경우와, 보고된 정상상태
관찰값으로부터 손으로 역산한 값인 경우가 섞여 있습니다. 어느 쪽인지는
`gdm_mrgsolve_model.R` 의 CALIBRATION 섹션에 항목별로 적어두었습니다.
모델은 개별 환자 데이터에 적합(fit)된 바 없습니다.
