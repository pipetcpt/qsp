# 남성 성선기능저하증 QSP 모델 — 참고문헌 (References)

**Male Hypogonadism · MHG**

이 문서는 `mhg_qsp_model.dot`(기계론적 지도), `mhg_mrgsolve_model.R`(ODE 모델),
`mhg_shiny_app.R`(대시보드)에 사용된 파라미터·구조·보정 표적의 출처를 정리한다.
모델 파일의 **CALIBRATION BLOCK**에 "ANCHORED / FITTED / CALIBRATED"로 분류한
각 항목이 아래 어느 문헌에 대응하는지 표시했다.

---

## 1. 진단 기준 및 가이드라인 (Diagnosis & Guidelines)

1. Bhasin S, et al. **Testosterone Therapy in Men With Hypogonadism: An Endocrine
   Society Clinical Practice Guideline.** J Clin Endocrinol Metab. 2018;103(5):1715-1744.
   <https://pubmed.ncbi.nlm.nih.gov/29562364/>
   — 진단 알고리즘, 아침 공복 총 T 2회 측정, 300 ng/dL 역치, 모니터링 일정.
   모델의 `MHG_diagnostic_frame()`이 검증 대상으로 삼는 규칙.

2. Travison TG, et al. **Harmonized Reference Ranges for Circulating Testosterone
   Levels in Men of Four Cohort Studies in the United States and Europe.**
   J Clin Endocrinol Metab. 2017;102(4):1161-1173.
   <https://pubmed.ncbi.nlm.nih.gov/28324103/>
   — CDC 조화화 기준치(19-39세 2.5백분위수 264 ng/dL). 300 ng/dL 역치의 근거.

3. Mulhall JP, et al. **Evaluation and Management of Testosterone Deficiency:
   AUA Guideline.** J Urol. 2018;200(2):423-432.
   <https://pubmed.ncbi.nlm.nih.gov/29601923/>

4. Salonia A, et al. **European Association of Urology Guidelines on Sexual and
   Reproductive Health — 2021 Update: Male Sexual Dysfunction.** Eur Urol.
   2021;80(3):333-357. <https://pubmed.ncbi.nlm.nih.gov/34183196/>

5. Rosner W, et al. **Position statement: Utility, limitations, and pitfalls in
   measuring testosterone: an Endocrine Society position statement.**
   J Clin Endocrinol Metab. 2007;92(2):405-413.
   <https://pubmed.ncbi.nlm.nih.gov/17090633/>
   — 면역측정법이 낮은 농도에서 부정확한 이유, LC-MS/MS 기준법.

---

## 2. 결합 평형 — 이 모델의 수학적 중심 (Binding Equilibrium)

6. Vermeulen A, Verdonck L, Kaufman JM. **A critical evaluation of simple methods
   for the estimation of free testosterone in serum.** J Clin Endocrinol Metab.
   1999;84(10):3666-3672. <https://pubmed.ncbi.nlm.nih.gov/10523012/>
   — **ANCHORED**: `KS = 1.0e9 M^-1`, `KA = 3.6e4 M^-1`, 그리고 모델 `$GLOBAL`의
   `free_T_pgmL()`가 푸는 2차방정식 그 자체.

7. Zakharov MN, et al. **A multi-step, dynamic allosteric model of testosterone's
   binding to sex hormone binding globulin.** Mol Cell Endocrinol. 2015;399:190-200.
   <https://pubmed.ncbi.nlm.nih.gov/25240469/>
   — SHBG 이량체의 알로스테릭 결합 모델. Vermeulen 단순 모델의 한계를 지적하며
   낮은 SHBG에서 유리 T를 과소평가할 수 있음을 보인다. 본 모델은 계산 가능성과
   임상 관행을 이유로 Vermeulen을 채택했으며, 이는 **의도적 단순화**다.

8. Handelsman DJ. **Free Testosterone: Pumping up the Tires or Ending the Free
   Ride?** Endocr Rev. 2017;38(4):297-301.
   <https://pubmed.ncbi.nlm.nih.gov/28486701/>
   — 유리 T 계산의 방법론적 비판. 모델이 유리 T를 "진실"이 아니라 "총 T와 다른
   프레임"으로 다루는 이유.

9. Goldman AL, et al. **A Reappraisal of Testosterone's Binding in Circulation:
   Physiological and Clinical Implications.** Endocr Rev. 2017;38(4):302-324.
   <https://pubmed.ncbi.nlm.nih.gov/28673039/>

10. Dunn JF, Nisula BC, Rodbard D. **Transport of steroid hormones: binding of 21
    endogenous steroids to both testosterone-binding globulin and
    corticosteroid-binding globulin in human plasma.** J Clin Endocrinol Metab.
    1981;53(1):58-68. <https://pubmed.ncbi.nlm.nih.gov/7195404/>
    — 알부민 결합 상수의 원전.

11. Ly LP, Handelsman DJ. **Empirical estimation of free testosterone from
    testosterone and sex hormone-binding globulin immunoassays.** Eur J Endocrinol.
    2005;152(3):471-478. <https://pubmed.ncbi.nlm.nih.gov/15757866/>

---

## 3. SHBG 조절 (SHBG Regulation)

12. Selby C. **Sex hormone binding globulin: origin, function and clinical
    significance.** Ann Clin Biochem. 1990;27(Pt 6):532-541.
    <https://pubmed.ncbi.nlm.nih.gov/2080856/>

13. Simó R, et al. **Novel insights in SHBG regulation and clinical implications.**
    Trends Endocrinol Metab. 2015;26(7):376-383.
    <https://pubmed.ncbi.nlm.nih.gov/26044465/>
    — **FITTED**: `EXP_INS`(인슐린), `EXP_E2`, `THYR`. HNF-4α를 매개로 한 조절.

14. Selva DM, et al. **Monosaccharide-induced lipogenesis regulates the human
    hepatic sex hormone-binding globulin gene.** J Clin Invest. 2007;117(12):3979-3987.
    <https://pubmed.ncbi.nlm.nih.gov/17992261/>
    — 단당류/de novo 지방생성 → HNF-4α ↓ → SHBG ↓. 지도의 `DNL_SHBG` 노드.

15. Wu FCW, et al. **Identification of late-onset hypogonadism in middle-aged and
    elderly men (EMAS).** N Engl J Med. 2010;363(2):123-135.
    <https://pubmed.ncbi.nlm.nih.gov/20554979/>
    — **FITTED**: 연령에 따른 SHBG 상승, 증상-호르몬 관계의 역치.

16. Feldman HA, et al. **Age trends in the level of serum testosterone and other
    hormones in middle-aged men: longitudinal results from the Massachusetts Male
    Aging Study.** J Clin Endocrinol Metab. 2002;87(2):589-598.
    <https://pubmed.ncbi.nlm.nih.gov/11836290/>
    — **ANCHORED**: 총 T 약 1%/년, 유리 T 약 2-3%/년 감소, SHBG 약 1%/년 상승
    (`AGE_S = 0.010`).

---

## 4. 테스토스테론 생산·청소 동태 (Production & Clearance)

17. Southren AL, et al. **Plasma production rates of testosterone in normal adult
    men and women.** J Clin Endocrinol Metab. 1965;25(11):1441-1450.
    <https://pubmed.ncbi.nlm.nih.gov/5842460/>
    — **ANCHORED**: 생산율 약 6 mg/day. 모델의 `KSPILL x ITT0 = 6000 ug/day` 항등식.

18. Vierhapper H, Nowotny P, Waldhäusl W. **Determination of testosterone
    production rates in men and women using stable isotope/dilution and mass
    spectrometry.** J Clin Endocrinol Metab. 1997;82(5):1492-1496.
    <https://pubmed.ncbi.nlm.nih.gov/9141539/>

19. Coviello AD, et al. **Effects of graded doses of testosterone on erythropoiesis
    in healthy young and older men.** J Clin Endocrinol Metab. 2008;93(3):914-919.
    <https://pubmed.ncbi.nlm.nih.gov/18160461/>
    — **FITTED**: 용량 의존적 Hb/Hct 상승, 고령에서 더 큼. 적혈구 모듈의 주 보정 표적.

---

## 5. 고환내 테스토스테론과 정자형성 (ITT & Spermatogenesis)

20. Jarow JP, Chen H, Rosner W, Trentacoste S, Zirkin BR. **Assessment of the
    androgen environment within the human testis: minimally invasive method to
    obtain intratesticular fluid.** J Androl. 2001;22(4):640-645.
    <https://pubmed.ncbi.nlm.nih.gov/11451361/>
    — **ANCHORED**: ITT가 혈청의 수십 배라는 핵심 관찰.

21. Coviello AD, et al. **Low-dose human chorionic gonadotropin maintains
    intratesticular testosterone in normal men with testosterone-induced
    gonadotropin suppression.** J Clin Endocrinol Metab. 2005;90(5):2595-2602.
    <https://pubmed.ncbi.nlm.nih.gov/15713727/>
    — **FITTED**: T 200 mg 주간 투여 시 ITT −94%; hCG 125-500 IU EOD의 용량 반응.
    `MHG_ITT_collapse()`와 `HCG_POT`, `ITT50_S`의 직접 보정 표적.

22. Roth MY, et al. **Serum LH correlates highly with intratesticular steroid
    levels in normal men.** J Androl. 2010;31(2):138-145.
    <https://pubmed.ncbi.nlm.nih.gov/19726700/>
    — **ANCHORED**: ITT 절대값 정량과 혈청 LH와의 상관. 모델의 `ITT0 = 700 nmol/L`
    및 `ITT_IN`이 LH_EQ의 포화 함수라는 구조의 근거.

23. Liu PY, et al. **Rate, extent, and modifiers of spermatogenic recovery after
    hormonal male contraception: an integrated analysis.** Lancet.
    2006;367(9520):1412-1420. <https://pubmed.ncbi.nlm.nih.gov/16650651/>
    — **FITTED**: 중단 후 회복 동태. 20 M/mL 도달까지 중앙값 3.4개월, 6개월 67%,
    12개월 90%, 24개월 100%. `MHG_recovery_curve()`의 대조 표적.

24. Heller CG, Clermont Y. **Kinetics of the germinal epithelium in man.**
    Recent Prog Horm Res. 1964;20:545-575.
    <https://pubmed.ncbi.nlm.nih.gov/14285045/>
    — **ANCHORED**: 정자형성 주기 약 74일 (`TAU_SPG`).

25. World Health Organization. **WHO laboratory manual for the examination and
    processing of human semen, 6th edition.** 2021.
    <https://www.who.int/publications/i/item/9789240030787>
    — 정자 농도 하한 기준치 16 M/mL (5판 15 M/mL에서 개정).

26. Ramasamy R, et al. **Testosterone supplementation versus clomiphene citrate
    for hypogonadism: an age matched comparison of satisfaction and efficacy.**
    J Urol. 2014;192(3):875-879. <https://pubmed.ncbi.nlm.nih.gov/24657837/>

27. Wenker EP, et al. **The Use of HCG-Based Combination Therapy for Recovery of
    Spermatogenesis after Testosterone Use.** J Sex Med. 2015;12(6):1334-1337.
    <https://pubmed.ncbi.nlm.nih.gov/25904023/>

---

## 6. 제형별 약동학 (Formulation Pharmacokinetics)

28. Behre HM, Nieschlag E. **Testosterone buciclate (20 Aet-1) in hypogonadal men:
    pharmacokinetics and pharmacodynamics.** J Clin Endocrinol Metab.
    1992;75(5):1204-1210. <https://pubmed.ncbi.nlm.nih.gov/1430080/>

29. Snyder PJ, Lawrence DA. **Treatment of male hypogonadism with testosterone
    enanthate.** J Clin Endocrinol Metab. 1980;51(6):1335-1339.
    <https://pubmed.ncbi.nlm.nih.gov/7440698/>
    — **FITTED**: 에난트산 IM의 첨두-저점 파형 (`KA_IM`, `F_IM`).

30. Swerdloff RS, et al. **Long-term pharmacokinetics of transdermal testosterone
    gel in hypogonadal men.** J Clin Endocrinol Metab. 2000;85(12):4500-4510.
    <https://pubmed.ncbi.nlm.nih.gov/11134099/>
    — **FITTED**: 경피 겔의 준정상상태 프로파일 (`KA_GEL`, `F_GEL`).

31. Wang C, et al. **Long-term testosterone gel (AndroGel) treatment maintains
    beneficial effects on sexual function and mood, lean and fat mass, and bone
    mineral density in hypogonadal men.** J Clin Endocrinol Metab.
    2004;89(5):2085-2098. <https://pubmed.ncbi.nlm.nih.gov/15126525/>

32. Morgentaler A, et al. **Long acting testosterone undecanoate therapy in men
    with hypogonadism: results of a pharmacokinetic clinical study.** J Urol.
    2008;180(6):2307-2313. <https://pubmed.ncbi.nlm.nih.gov/18930277/>
    — **FITTED**: 운데칸산 IM의 매우 느린 방출 (`KA_TU`).

33. Swerdloff RS, et al. **A New Oral Testosterone Undecanoate Formulation
    Restores Testosterone to Normal Concentrations in Hypogonadal Men.**
    J Clin Endocrinol Metab. 2020;105(8):2515-2531.
    <https://pubmed.ncbi.nlm.nih.gov/32427336/>
    — **FITTED**: 경구 TU의 림프 흡수, 식사 의존성, 높은 첨두 (`F_ORAL`, `KA_ORAL`).

34. Kaminetsky JC, et al. **A 52-Week Study of Dose Adjusted Subcutaneous
    Testosterone Enanthate in Oil Self-Administered via Disposable Auto-Injector.**
    J Urol. 2019;201(3):587-594. <https://pubmed.ncbi.nlm.nih.gov/30273615/>
    — **FITTED**: 주 1회 SC 자동주사기(Xyosted)의 완만한 파형 (`KA_SC`, `F_SC`).

35. Kaminetsky J, et al. **Efficacy and Safety of Natesto, a Testosterone Nasal
    Gel, in Hypogonadal Men.** J Sex Med. 2017;14(9):1109-1116.
    <https://pubmed.ncbi.nlm.nih.gov/28781215/>
    — 비강 겔이 LH/FSH 억제를 상대적으로 덜 일으킨다는 관찰 — 짧은 첨두의 결과.

36. Kovac JR, et al. **Testosterone pellet implants: a review.** Asian J Androl.
    2015;17(6):938-943. <https://pubmed.ncbi.nlm.nih.gov/25994641/>

---

## 7. 적혈구증가증 — 모델의 두 번째 볼록성 (Erythrocytosis)

37. Bachman E, et al. **Testosterone induces erythrocytosis via increased
    erythropoietin and suppressed hepcidin: evidence for a new erythropoietin/
    hemoglobin set point.** J Gerontol A Biol Sci Med Sci. 2014;69(6):725-735.
    <https://pubmed.ncbi.nlm.nih.gov/24158761/>
    — **FITTED**: 모델의 헵시딘-EPO 설정점 구조 그 자체. `IMAX_H`, `SEPO`.

38. Guo W, et al. **Testosterone Administration Inhibits Hepcidin Transcription
    and Is Associated With Increased Iron Incorporation Into Red Blood Cells.**
    Aging Cell. 2013;12(2):280-291. <https://pubmed.ncbi.nlm.nih.gov/23399021/>
    — 헵시딘 억제 → ferroportin 유지 → 철 이용성 증가. `IRON_F` 항.

39. Ohlander SJ, Varghese B, Pastuszak AW. **Erythrocytosis Following Testosterone
    Therapy.** Sex Med Rev. 2018;6(1):77-85.
    <https://pubmed.ncbi.nlm.nih.gov/28526632/>
    — **주요 대조 표적**: 근주 제형이 경피 제형보다 적혈구증가증을 약 3배 흔하게
    유발한다는 관찰. `MHG_convexity_decomposition()`이 설명하려는 바로 그 격차.

40. Jones SD Jr, Dukovac T, Sangkum P, Yafi FA, Hellstrom WJG. **Erythrocytosis
    and Polycythemia Secondary to Testosterone Replacement Therapy in the Aging
    Male.** Sex Med Rev. 2015;3(2):101-112.
    <https://pubmed.ncbi.nlm.nih.gov/27784544/>

41. Roy CN, et al. **Association of Testosterone Levels With Anemia in Older Men:
    A Controlled Clinical Trial (T-Trials Anemia Trial).** JAMA Intern Med.
    2017;177(4):480-490. <https://pubmed.ncbi.nlm.nih.gov/28241237/>
    — **FITTED**: 원인불명 빈혈의 54% 교정 (위약 15%). 모델에서 성선기능저하 남성이
    기저에서 Hct가 낮게 출발하는 구조의 근거.

---

## 8. 골 — 에스트로겐 의존성 (Bone)

42. Finkelstein JS, et al. **Gonadal steroids and body composition, strength, and
    sexual function in men.** N Engl J Med. 2013;369(11):1011-1022.
    <https://pubmed.ncbi.nlm.nih.gov/24024838/>
    — **핵심 문헌 · FITTED**: 고세렐린 배경 + 단계적 T 겔 ± 아나스트로졸.
    제지방·근력은 T 의존, **체지방은 E2 의존**, 성기능은 양쪽 모두.
    `MHG_finkelstein()`이 재현하려는 설계. `SFAT_E`, `SLEAN_T`의 직접 보정 표적.

43. Finkelstein JS, et al. **Gonadal Steroid-Dependent Effects on Bone Turnover
    and Bone Mineral Density in Men.** J Clin Invest. 2016;126(3):1114-1125.
    <https://pubmed.ncbi.nlm.nih.gov/26901812/>
    — 골 교체와 BMD에 대한 E2 우세 효과. `E_SCL`, `E_OC`.

44. Snyder PJ, et al. **Effect of Testosterone Treatment on Volumetric Bone
    Density and Strength in Older Men With Low Testosterone (T-Trials Bone Trial).**
    JAMA Intern Med. 2017;177(4):471-479.
    <https://pubmed.ncbi.nlm.nih.gov/28241231/>
    — **FITTED**: 요추 해면골 vBMD 1년에 약 +7.5%. `KFORM`/`KRES`의 보정 표적.

45. Khosla S, et al. **Relationship of serum sex steroid levels and bone turnover
    markers with bone mineral density in men and women: a key role for bioavailable
    estrogen.** J Clin Endocrinol Metab. 1998;83(7):2266-2274.
    <https://pubmed.ncbi.nlm.nih.gov/9661593/>

46. Smith EP, et al. **Estrogen resistance caused by a mutation in the
    estrogen-receptor gene in a man.** N Engl J Med. 1994;331(16):1056-1061.
    <https://pubmed.ncbi.nlm.nih.gov/8090165/>
    — 남성에서 에스트로겐이 골에 필수적임을 보인 자연 실험.

47. Snyder PJ, et al. **Testosterone Treatment and Fractures in Men with
    Hypogonadism (TRAVERSE fracture substudy).** N Engl J Med. 2024;390(3):203-211.
    <https://pubmed.ncbi.nlm.nih.gov/38231624/>
    — **모델이 재현하지 못하는 결과**: 임상 골절 3.50% vs 2.46% (HR 1.43).
    BMD는 올라갔는데 골절은 늘었다. `MHG_trial_ledger()`가 이 불일치를 숨기지 않고
    출력하는 이유이며, BMD를 대리지표로 삼는 가정 자체의 한계를 드러낸다.

---

## 9. 체성분·근육 (Body Composition & Muscle)

48. Bhasin S, et al. **Testosterone dose-response relationships in healthy young
    men.** Am J Physiol Endocrinol Metab. 2001;281(6):E1172-E1181.
    <https://pubmed.ncbi.nlm.nih.gov/11701431/>
    — **FITTED**: 제지방량의 용량-반응. `SLEAN_T`, `KL50`.

49. Bhasin S, et al. **The effects of supraphysiologic doses of testosterone on
    muscle size and strength in normal men.** N Engl J Med. 1996;335(1):1-7.
    <https://pubmed.ncbi.nlm.nih.gov/8637535/>

50. Storer TW, et al. **Effects of Testosterone Supplementation for 3 Years on
    Muscle Performance and Physical Function in Older Men.** J Clin Endocrinol
    Metab. 2017;102(2):583-593. <https://pubmed.ncbi.nlm.nih.gov/27754805/>

---

## 10. 임상 결과 시험 (Clinical Outcome Trials)

51. Snyder PJ, et al. **Effects of Testosterone Treatment in Older Men
    (The Testosterone Trials).** N Engl J Med. 2016;374(7):611-624.
    <https://pubmed.ncbi.nlm.nih.gov/26886521/>
    — **FITTED**: 성기능(PDQ-Q4 +0.58 vs +0.10)은 개선, **활력(FACIT-F)과
    보행거리 1차 종점은 미달성**. 모델의 PRO 가중치를 일부러 작게 잡은 근거.

52. Resnick SM, et al. **Testosterone Treatment and Cognitive Function in Older
    Men With Low Testosterone and Age-Associated Memory Impairment.** JAMA.
    2017;317(7):717-727. <https://pubmed.ncbi.nlm.nih.gov/28241356/>
    — 인지 개선 없음.

53. Budoff MJ, et al. **Testosterone Treatment and Coronary Artery Plaque Volume
    in Older Men With Low Testosterone.** JAMA. 2017;317(7):708-716.
    <https://pubmed.ncbi.nlm.nih.gov/28241355/>
    — 비석회화 관상동맥 플라크 용적 증가.

54. Lincoff AM, et al. **Cardiovascular Safety of Testosterone-Replacement Therapy
    (TRAVERSE).** N Engl J Med. 2023;389(2):107-117.
    <https://pubmed.ncbi.nlm.nih.gov/37326322/>
    — **주요 대조 표적**: n=5,246, MACE 7.0% vs 7.3% (HR 0.96, 95% CI 0.78-1.17,
    비열등성 충족). 다만 심방세동 3.5% vs 2.4%, 폐색전 0.9% vs 0.5%,
    급성 신손상 2.3% vs 1.5%로 증가. 결정론적 모델은 사건 과정이 없어 이 결과를
    예측하지 않으며, `MHG_trial_ledger()`가 그 한계를 명시한다.

55. Vigen R, et al. **Association of testosterone therapy with mortality,
    myocardial infarction, and stroke in men with low testosterone levels.**
    JAMA. 2013;310(17):1829-1836. <https://pubmed.ncbi.nlm.nih.gov/24193080/>
    — TRAVERSE 이전 논쟁의 발단이 된 관찰연구 (방법론 비판 다수).

---

## 11. 원인별 병태생리 (Aetiology)

56. Bojesen A, Juul S, Gravholt CH. **Prenatal and postnatal prevalence of
    Klinefelter syndrome: a national registry study.** J Clin Endocrinol Metab.
    2003;88(2):622-626. <https://pubmed.ncbi.nlm.nih.gov/12574191/>

57. Groth KA, et al. **Klinefelter Syndrome — A Clinical Update.** J Clin
    Endocrinol Metab. 2013;98(1):20-30.
    <https://pubmed.ncbi.nlm.nih.gov/23118429/>

58. Boehm U, et al. **Expert consensus document: European Consensus Statement on
    congenital hypogonadotropic hypogonadism — pathogenesis, diagnosis and
    treatment.** Nat Rev Endocrinol. 2015;11(9):547-564.
    <https://pubmed.ncbi.nlm.nih.gov/26194704/>
    — 칼만 증후군/nIHH의 유전학(ANOS1, FGFR1, PROKR2, CHD7). 지도의 ⑩ 클러스터.

59. Daniell HW. **Hypogonadism in men consuming sustained-action oral opioids.**
    J Pain. 2002;3(5):377-384. <https://pubmed.ncbi.nlm.nih.gov/14622741/>
    — **FITTED**: `OPI_SUP`. 만성 아편유사제 사용자에서의 높은 유병률.

60. Grossmann M, Matsumoto AM. **A Perspective on Middle-Aged and Older Men With
    Functional Hypogonadism: Focus on Holistic Management.** J Clin Endocrinol
    Metab. 2017;102(3):1067-1075. <https://pubmed.ncbi.nlm.nih.gov/28359097/>
    — 기능성(가역적) 성선기능저하증 개념. 지도의 ⑪ 클러스터.

61. Corona G, et al. **Body weight loss reverts obesity-associated
    hypogonadotropic hypogonadism: a systematic review and meta-analysis.**
    Eur J Endocrinol. 2013;168(6):829-843.
    <https://pubmed.ncbi.nlm.nih.gov/23482592/>
    — **FITTED**: `WT_LOSS`. 체중 감량 정도에 따른 총 T 회복 폭.

62. Rastrelli G, et al. **Testosterone Replacement Therapy for Sexual Symptoms.**
    Sex Med Rev. 2019;7(3):464-475. <https://pubmed.ncbi.nlm.nih.gov/30926459/>

63. Rochira V, et al. **Hypogonadism in the HIV-infected man.** Endocrinol Metab
    Clin North Am. 2014;43(3):709-730.
    <https://pubmed.ncbi.nlm.nih.gov/25169563/>

---

## 12. 비안드로겐 전략 (Non-Androgen Strategies)

64. Wheeler KM, et al. **A review of the role of human chorionic gonadotropin
    (hCG) in male reproduction.** Transl Androl Urol. 2019;8(Suppl 1):S6-S17.
    <https://pubmed.ncbi.nlm.nih.gov/31143669/>

65. Rambhatla A, et al. **The Role of Estrogen Modulators in Male Hypogonadism
    and Infertility.** Rev Urol. 2016;18(2):66-72.
    <https://pubmed.ncbi.nlm.nih.gov/27601965/>
    — 클로미펜/엔클로미펜/아나스트로졸의 기전 비교. 지도의 ⑲ 클러스터.

66. Kim ED, et al. **Oral enclomiphene citrate raises testosterone and preserves
    sperm counts in obese hypogonadal men, unlike topical testosterone: restoration
    instead of replacement.** BJU Int. 2016;117(4):677-685.
    <https://pubmed.ncbi.nlm.nih.gov/26496621/>
    — 제목 자체가 이 모델의 ⑥번 분석 함수가 말하려는 바다: 대체(replacement)와
    회복(restoration)은 다른 개입이다.

67. Burnett-Bowie SA, et al. **Effects of aromatase inhibition on bone mineral
    density and bone turnover in older men with low testosterone levels.**
    J Clin Endocrinol Metab. 2009;94(12):4785-4792.
    <https://pubmed.ncbi.nlm.nih.gov/19820017/>
    — **FITTED**: `MHG_aromatase_cost()`의 대조 표적. 아로마타제 억제가 총 T를
    올리면서 요추 BMD를 떨어뜨린다.

68. Dwyer AA, et al. **Trial of Recombinant Follicle-Stimulating Hormone Pretreatment
    for GnRH-Induced Fertility in Patients with Congenital Hypogonadotropic
    Hypogonadism.** J Clin Endocrinol Metab. 2013;98(11):E1790-E1795.
    <https://pubmed.ncbi.nlm.nih.gov/24037890/>

69. Jayasena CN, et al. **Twice-weekly administration of kisspeptin-54 for 8 weeks
    stimulates release of reproductive hormones in women with hypothalamic
    amenorrhea.** Clin Pharmacol Ther. 2010;88(6):840-847.
    <https://pubmed.ncbi.nlm.nih.gov/20980998/>
    — 키스펩틴 유사체의 축 재기동 개념 (남성 적용은 연구 단계).

---

## 13. 전립선 안전성 — 포화 모델 (Prostate Safety)

70. Morgentaler A, Traish AM. **Shifting the paradigm of testosterone and prostate
    cancer: the saturation model and the limits of androgen-dependent growth.**
    Eur Urol. 2009;55(2):310-320. <https://pubmed.ncbi.nlm.nih.gov/18838208/>
    — **FITTED**: `KPSA`를 `DHT0`보다 낮게 잡아 포화가 생리적 범위 아래에서
    일어나게 한 구조적 선택의 근거.

71. Boyle P, et al. **Endogenous and exogenous testosterone and the risk of
    prostate cancer and increased prostate-specific antigen (PSA): a
    meta-analysis.** BJU Int. 2016;118(5):731-741.
    <https://pubmed.ncbi.nlm.nih.gov/27313122/>

---

## 14. KNDy / GnRH 박동 생성기 (Pulse Generator)

72. Lehman MN, Coolen LM, Goodman RL. **Minireview: kisspeptin/neurokinin B/
    dynorphin (KNDy) cells of the arcuate nucleus: a central node in the control
    of gonadotropin-releasing hormone secretion.** Endocrinology.
    2010;151(8):3479-3489. <https://pubmed.ncbi.nlm.nih.gov/20501670/>
    — 지도 ① 클러스터의 구조.

73. Pitteloud N, et al. **Inhibition of luteinizing hormone secretion by
    testosterone in men requires aromatization for its pituitary but not its
    hypothalamic effects: evidence from the tandem study of normal and
    gonadotropin-releasing hormone-deficient men.** J Clin Endocrinol Metab.
    2008;93(3):784-791. <https://pubmed.ncbi.nlm.nih.gov/18073301/>
    — **FITTED**: `WE2 = 2.2`. 남성 되먹임에서 E2가 지배적이라는 직접 증거.

74. Veldhuis JD, et al. **Testosterone blunts feedback inhibition of gonadotropin
    secretion by exogenous estradiol in healthy older but not young men.**
    J Clin Endocrinol Metab. 2008;93(3):874-880.
    <https://pubmed.ncbi.nlm.nih.gov/18160462/>

75. Keenan DM, Veldhuis JD. **Pulsatility of hypothalamo-pituitary hormones: a
    challenge in quantification.** Physiology (Bethesda). 2016;31(1):34-50.
    <https://pubmed.ncbi.nlm.nih.gov/26674550/>
    — 본 모델이 LH 박동을 명시적으로 시뮬레이션하지 않고 평균 농도로 다루는
    선택의 배경(그리고 그 대가).

---

## 15. QSP 방법론 (QSP Methodology)

76. Baron KT, Gastonguay MR. **mrgsolve: Simulate from ODE-Based Population PK/PD
    and Systems Pharmacology Models.** <https://mrgsolve.org/>

77. Musante CJ, et al. **Quantitative Systems Pharmacology: A Case for
    Disease Models.** Clin Pharmacol Ther. 2017;101(1):24-27.
    <https://pubmed.ncbi.nlm.nih.gov/27709613/>

78. Ribba B, et al. **Model-Informed Drug Development for Quantitative Systems
    Pharmacology.** CPT Pharmacometrics Syst Pharmacol. 2017;6(8):496-498.
    <https://pubmed.ncbi.nlm.nih.gov/28571121/>

---

## 부록: 모델이 재현하지 못하는 것 (What This Model Does Not Reproduce)

정직한 QSP 모델은 자신이 틀리는 지점을 명시해야 한다. 다음은 의도적으로
재현을 시도하지 않았거나, 시도했으나 실패한 항목이다.

| 항목 | 상태 | 이유 |
|------|------|------|
| TRAVERSE MACE (HR 0.96) | 재현 시도 안 함 | 5,246명의 사건 과정이며 결정론적 단일 환자 모델에는 대응하는 구조가 없다 |
| TRAVERSE 골절 (HR 1.43) | **재현 실패 (의도적 노출)** | 모델은 BMD 상승을 재현하므로 골절 증가를 예측할 수 없다. BMD를 대리지표로 쓰는 가정의 실패이지 파라미터 조정 문제가 아니다 (문헌 47) |
| LH 박동성 | 단순화 | 평균 농도로 대체. 박동 빈도가 LHB/FSHB 전사를 차등 조절하는 기전은 지도에만 있고 ODE에는 없다 (문헌 75) |
| SHBG 알로스테릭 결합 | 단순화 | Vermeulen 단순 모델 채택. 낮은 SHBG에서 유리 T를 과소평가할 수 있다 (문헌 7) |
| 개체간 변이 (AR CAG 반복 등) | 미구현 | 모델은 결정론적이며 집단 분포는 `HCT_SD` 한 곳에서만 근사적으로 쓰인다 |
| `EC50_HEPC` | **보정값이며 측정값 아님** | 제형별 적혈구증가증 격차를 설명하는 핵심 가정. `MHG_hepcidin_sensitivity()`가 이 값에 대한 결론의 민감도를 출력한다 |

---

> **면책 조항 (Disclaimer)** — 본 모델과 문헌 정리는 교육 및 연구 목적입니다.
> 임상 의사결정, 처방, 규제 제출에 직접 사용해서는 안 됩니다.
