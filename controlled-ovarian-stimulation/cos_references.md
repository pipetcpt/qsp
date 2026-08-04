# 조절 난소 자극 (Controlled Ovarian Stimulation) — 참고문헌

이 QSP 모델(기계론적 지도 · 65-ODE mrgsolve 모델 · Shiny 대시보드)의 근거가 된
문헌 목록입니다. **모든 PMID는 NCBI E-utilities(`esearch` → `esummary`)로 개별
조회하여 제1저자·연도·저널·제목이 일치하는지 확인했습니다.** 해석되지 않거나
제목이 일치하지 않은 후보는 목록에서 제외했습니다.

모델의 수치 보정에 직접 사용된 문헌은 ★로 표시하고, 어떤 파라미터의 근거인지
적었습니다.

---

## 1. FSH 임계값·창(threshold/window) 개념과 난포 선택

이 모델의 2번 클러스터 전체(임계값 분포·선택·우성화)는 이 개념 위에 서 있습니다.
난포 수를 지정하는 파라미터가 없는 것도 이 때문입니다.

1. ★ Brown JB. Pituitary control of ovarian function—concepts derived from gonadotrophin therapy. *Aust N Z J Obstet Gynaecol* 1978. — [PMID 278588](https://pubmed.ncbi.nlm.nih.gov/278588/) — 임계값/창 개념의 원전. 모델의 `T50`·`SIGT`
2. ★ Baird DT. A model for follicular selection and ovulation: lessons from superovulation. *J Steroid Biochem* 1987. — [PMID 3121918](https://pubmed.ncbi.nlm.nih.gov/3121918/) — "외인성 FSH는 난포를 만들지 않고 폐쇄를 막는다"는 구조의 근거
3. ★ Fauser BC, Van Heusden AM. Manipulation of human ovarian function: physiological concepts and clinical consequences. *Endocr Rev* 1997. — [PMID 9034787](https://pubmed.ncbi.nlm.nih.gov/9034787/) — 모델의 자연/자극 주기 대칭성
4. ★ Schipper I, et al. The follicle-stimulating hormone (FSH) threshold/window concept examined by different interventions with exogenous FSH during the follicular phase of the normal menstrual cycle. *J Clin Endocrinol Metab* 1998. — [PMID 9543158](https://pubmed.ncbi.nlm.nih.gov/9543158/) — `HT = 4` (임계값 Hill 계수)의 근거
5. ★ Hillier SG. Current concepts of the roles of follicle stimulating hormone and luteinizing hormone in folliculogenesis. *Hum Reprod* 1994. — [PMID 8027271](https://pubmed.ncbi.nlm.nih.gov/8027271/) — 11 mm 이후 LH 의존으로의 전환 (`KLHA`, `DLH`)
6. Macklon NS, Fauser BC. Aspects of ovarian follicle development throughout life. *Horm Res* 1999. — [PMID 10725781](https://pubmed.ncbi.nlm.nih.gov/10725781/)
7. Mumusoglu S, et al. Initial and cyclic recruitment of ovarian follicles: a quarter-century update. *Reprod Biomed Online* 2025. — [PMID 41067165](https://pubmed.ncbi.nlm.nih.gov/41067165/)
8. Lyu Z, et al. Stochastic mechanism of dominant follicle selection: selection of one suppresses selection of others. *J R Soc Interface* 2026. — [PMID 42014065](https://pubmed.ncbi.nlm.nih.gov/42014065/) — 본 모델과 다른(확률적) 우성화 기전. 대안 가설로 기재

## 2. 스테로이드 생성 — 2세포 2고나도트로핀 구조

3번 클러스터. 테카(LH) → 안드로겐 → 과립세포(FSH) 아로마타제 → E2 라는 곱셈
구조가 모델에서 "LH를 깊게 억제하면 난포당 E2가 떨어진다"를 만들어 냅니다.

9. ★ Hillier SG, Whitelaw PF, Smyth CD. Follicular oestrogen synthesis: the 'two-cell, two-gonadotrophin' model revisited. *Mol Cell Endocrinol* 1994. — [PMID 8056158](https://pubmed.ncbi.nlm.nih.gov/8056158/) — `KTHECA = 0.40 IU/L`
10. Hillier SG. Gonadotropic control of ovarian follicular growth and development. *Mol Cell Endocrinol* 2001. — [PMID 11420129](https://pubmed.ncbi.nlm.nih.gov/11420129/)
11. Magoffin DA. The ovarian androgen-producing cells: a 2001 perspective. *Rev Endocr Metab Disord* 2002. — [PMID 11883104](https://pubmed.ncbi.nlm.nih.gov/11883104/)
12. Hattori N, et al. Epoxide hydrolase affects estrogen production in the human ovary. *Endocrinology* 2000. — [PMID 10965908](https://pubmed.ncbi.nlm.nih.gov/10965908/)

## 3. AMH — 예비능 표지자이자 임계값 인자

모델에서 AMH는 두 역할을 합니다: 코호트 크기의 관측 가능한 대리변수이면서,
동시에 FSH 임계값을 올리는 항(`KAMHT`)입니다.

13. ★ Weenen C, et al. Anti-Müllerian hormone expression pattern in the human ovary: potential implications for initial and cyclic follicle recruitment. *Mol Hum Reprod* 2004. — [PMID 14742691](https://pubmed.ncbi.nlm.nih.gov/14742691/) — AMH가 소난포에서만 분비된다는 사실 = 모델의 `DAMH = 8 mm`
14. ★ Pellatt L, et al. Granulosa cell production of anti-Müllerian hormone is increased in polycystic ovaries. *J Clin Endocrinol Metab* 2007. — [PMID 17062765](https://pubmed.ncbi.nlm.nih.gov/17062765/) — PCOS 표현형의 `KAMHT` 경로
15. Dewailly D, et al. The physiology and clinical utility of anti-Müllerian hormone in women. *Hum Reprod Update* 2014. — [PMID 24430863](https://pubmed.ncbi.nlm.nih.gov/24430863/)
16. La Marca A, et al. Anti-Müllerian hormone (AMH) as a predictive marker in assisted reproductive technology (ART). *Hum Reprod Update* 2010. — [PMID 19793843](https://pubmed.ncbi.nlm.nih.gov/19793843/)
17. ★ Broer SL, et al. AMH and AFC as predictors of excessive response in controlled ovarian hyperstimulation: a meta-analysis. *Hum Reprod Update* 2011. — [PMID 20667894](https://pubmed.ncbi.nlm.nih.gov/20667894/) — AFC 25–32 고반응군 정의
18. Broekmans FJ, et al. The antral follicle count: practical recommendations for better standardization. *Fertil Steril* 2010. — [PMID 19589513](https://pubmed.ncbi.nlm.nih.gov/19589513/) — `AFC` 정의 (2–9 mm)

## 4. 고나도트로핀 약동학

모델의 rFSH PK(F 0.80 · V/F 26 L · 다회투여 t½ 42 h → 150 IU/일에서 11.9 IU/L)와
corifollitropin 구획의 근거입니다.

19. ★ le Cotonnec JY, et al. Clinical pharmacology of recombinant human follicle-stimulating hormone (FSH). I. Comparative pharmacokinetics with urinary human FSH. *Fertil Steril* 1994. — [PMID 8150109](https://pubmed.ncbi.nlm.nih.gov/8150109/) — `KAF`·`KELFX`·`VF`
20. ★ Duijkers IJ, et al. Single dose pharmacokinetics and effects on follicular growth and serum hormones of a long-acting recombinant FSH preparation (FSH-CTP) in healthy pituitary-suppressed females. *Hum Reprod* 2002. — [PMID 12151425](https://pubmed.ncbi.nlm.nih.gov/12151425/) — corifollitropin t½ 69 h
21. Fauser BC, et al. Advances in recombinant DNA technology: corifollitropin alfa, a hybrid molecule with sustained follicle-stimulating activity and reduced injection frequency. *Hum Reprod Update* 2009. — [PMID 19182099](https://pubmed.ncbi.nlm.nih.gov/19182099/)
22. Corifollitropin Alfa Dose-finding Study Group. A randomized dose-response trial of a single injection of corifollitropin alfa to sustain multifollicular growth during controlled ovarian stimulation. *Hum Reprod* 2008. — [PMID 18684735](https://pubmed.ncbi.nlm.nih.gov/18684735/) — `POTCO`(FSH 등가 역가) 보정
23. Olsson H, Sandström R, Grundemar L. Different pharmacokinetic and pharmacodynamic properties of recombinant follicle-stimulating hormone (rFSH) derived from a human cell line compared with rFSH from a non-human cell line. *J Clin Pharmacol* 2014. — [PMID 24800998](https://pubmed.ncbi.nlm.nih.gov/24800998/) — follitropin delta

## 5. GnRH 길항제 — 약동학과 서지 차단

24. ★ Oberyé JJ, et al. Pharmacokinetic and pharmacodynamic characteristics of ganirelix (Antagon/Orgalutran). Part I. Absolute bioavailability of 0.25 mg of ganirelix after a single subcutaneous injection in healthy female volunteers. *Fertil Steril* 1999. — [PMID 10593371](https://pubmed.ncbi.nlm.nih.gov/10593371/) — `FBIOANT 0.91`, Cmax 11.2 ng/mL, tmax 1.1 h
25. ★ Oberyé JJ, et al. Pharmacokinetic and pharmacodynamic characteristics of ganirelix (Antagon/Orgalutran). Part II. Dose-proportionality and gonadotropin suppression after multiple doses of ganirelix in healthy female volunteers. *Fertil Steril* 1999. — [PMID 10593372](https://pubmed.ncbi.nlm.nih.gov/10593372/) — LH 억제 70–80% = 모델의 77% (`IC50ANT`, `FANTL`)
26. ★ Ganirelix Dose-finding Study Group. A double-blind, randomized, dose-finding study to assess the efficacy of the gonadotrophin-releasing hormone antagonist ganirelix (Org 37462) to prevent premature luteinizing hormone surges in women undergoing ovarian stimulation with recombinant follicle stimulating hormone (Puregon). *Hum Reprod* 1998. — [PMID 9853849](https://pubmed.ncbi.nlm.nih.gov/9853849/) — 0.25 mg 선택의 근거
27. Lambalk CB, et al. GnRH antagonist versus long agonist protocols in IVF: a systematic review and meta-analysis accounting for patient type. *Hum Reprod Update* 2017. — [PMID 28903472](https://pubmed.ncbi.nlm.nih.gov/28903472/)
28. ★ Kuang Y, et al. Medroxyprogesterone acetate is an effective oral alternative for preventing premature luteinizing hormone surges in women undergoing controlled ovarian hyperstimulation for in vitro fertilization. *Fertil Steril* 2015. — [PMID 25956370](https://pubmed.ncbi.nlm.nih.gov/25956370/) — PPOS arm (`PPOS`, `KP4S`)

## 6. 트리거 — 이 모델의 핵심 주장이 검증되는 곳

hCG와 GnRH 작용제는 같은 수용체를 쓰지만, 모델에서는 성숙(앞머리)과 VEGF(면적)를
서로 다른 커널로 읽기 때문에 두 약이 갈라집니다.

29. ★ Humaidan P, et al. GnRH agonist (buserelin) or hCG for ovulation induction in GnRH antagonist IVF/ICSI cycles: a prospective randomized study. *Hum Reprod* 2005. — [PMID 15760966](https://pubmed.ncbi.nlm.nih.gov/15760966/) — 작용제 트리거의 황체기 결핍. 모델은 방향만 재현(−41%)하고 크기는 과소예측
30. ★ Kolibianakis EM, et al. A lower ongoing pregnancy rate can be expected when GnRH agonist is used for triggering final oocyte maturation instead of HCG in patients undergoing IVF with GnRH antagonists. *Hum Reprod* 2005. — [PMID 15979994](https://pubmed.ncbi.nlm.nih.gov/15979994/)
31. ★ Fauser BC, et al. Endocrine profiles after triggering of final oocyte maturation with GnRH agonist after cotreatment with the GnRH antagonist ganirelix during ovarian hyperstimulation for in vitro fertilization. *J Clin Endocrinol Metab* 2002. — [PMID 11836309](https://pubmed.ncbi.nlm.nih.gov/11836309/) — 작용제 플레어의 LH 최고치·지속시간 (`AMPA`, `SLHMAX`)
32. ★ Youssef MA, et al. Gonadotropin-releasing hormone agonist versus HCG for oocyte triggering in antagonist-assisted reproductive technology. *Cochrane Database Syst Rev* 2014. — [PMID 25358904](https://pubmed.ncbi.nlm.nih.gov/25358904/) — OHSS 감소·신선이식 임신율 저하의 메타분석
33. ★ Griffin D, et al. Dual trigger of oocyte maturation with gonadotropin-releasing hormone agonist and low-dose human chorionic gonadotropin to optimize live birth rates in high responders. *Fertil Steril* 2012. — [PMID 22480822](https://pubmed.ncbi.nlm.nih.gov/22480822/) — 이중 트리거 arm
34. ★ Shapiro BS, et al. Comparison of "triggers" using leuprolide acetate alone or in combination with low-dose human chorionic gonadotropin. *Fertil Steril* 2011. — [PMID 21550042](https://pubmed.ncbi.nlm.nih.gov/21550042/)
35. Humaidan P, et al. GnRH agonist for triggering of final oocyte maturation: time for a change of practice? *Hum Reprod Update* 2011. — [PMID 21450755](https://pubmed.ncbi.nlm.nih.gov/21450755/)
36. Humaidan P, et al. GnRHa trigger and individualized luteal phase hCG support according to ovarian response to stimulation: two prospective randomized controlled multi-centre studies in IVF patients. *Hum Reprod* 2013. — [PMID 23753114](https://pubmed.ncbi.nlm.nih.gov/23753114/)
37. Yding Andersen C, et al. Improving the luteal phase after ovarian stimulation: reviewing new options. *Reprod Biomed Online* 2014. — [PMID 24656557](https://pubmed.ncbi.nlm.nih.gov/24656557/)
38. ★ Beckers NG, et al. Nonsupplemented luteal phase characteristics after the administration of recombinant human chorionic gonadotropin, recombinant luteinizing hormone, or gonadotropin-releasing hormone (GnRH) agonist to induce final oocyte maturation in in vitro fertilization patients after ovarian stimulation with recombinant follicle-stimulating hormone and GnRH antagonist cotreatment. *J Clin Endocrinol Metab* 2003. — [PMID 12970285](https://pubmed.ncbi.nlm.nih.gov/12970285/) — 황체기 P4 프로파일. 모델의 `KP4CL`·`KSATCL`·`KLYSM`

## 7. 난자 성숙과 채취 시각 — 두 시계 사이의 간격

39. ★ Park JY, et al. EGF-like growth factors as mediators of LH action in the ovulatory follicle. *Science* 2004. — [PMID 14726596](https://pubmed.ncbi.nlm.nih.gov/14726596/) — 성숙이 "앞머리 신호 + 자율 진행"이라는 구조의 분자적 근거
40. ★ Mansour RT, Aboulghar MA, Serour GI. Study of the optimum time for human chorionic gonadotropin–ovum pickup interval in in vitro fertilization. *J Assist Reprod Genet* 1994. — [PMID 7633170](https://pubmed.ncbi.nlm.nih.gov/7633170/) — 34–38 h 창. 모델의 `TMAT 0.85 d`·`TRUP 1.55 d`·`RUPX 0.90 d`
41. Dozortsev DI, Diamond MP. Luteinizing hormone-independent rise of progesterone as the physiological trigger of the ovulatory gonadotropins surge in the human. *Fertil Steril* 2020. — [PMID 32741458](https://pubmed.ncbi.nlm.nih.gov/32741458/) — 모델이 채택하지 않은 대안 트리거 가설
42. Jia Z, et al. Effects of C-type natriuretic peptide on meiotic arrest and developmental competence of oocytes. *Sci Rep* 2020. — [PMID 33106527](https://pubmed.ncbi.nlm.nih.gov/33106527/) — cGMP 감수분열 정지 (지도 8번 클러스터)

## 8. OHSS — VEGF 매개 혈관 투과성

모델에서 OHSS는 중증도 척도가 아니라 (과립세포 질량 × LHCGR 점유 면적)의 결과로
계산됩니다.

43. ★ McClure N, et al. Vascular endothelial growth factor as capillary permeability agent in ovarian hyperstimulation syndrome. *Lancet* 1994. — [PMID 7913160](https://pubmed.ncbi.nlm.nih.gov/7913160/) — VEGF → 투과성 축의 원전
44. ★ Gómez R, et al. Vascular endothelial growth factor receptor-2 activation induces vascular permeability in hyperstimulated rats, and this effect is prevented by dopamine receptor-2 activation. *Endocrinology* 2002. — [PMID 12399430](https://pubmed.ncbi.nlm.nih.gov/12399430/) — 카베르골린 작용점 (`EMAXCAB`)
45. ★ Soares SR, et al. Targeting the vascular endothelial growth factor system to prevent ovarian hyperstimulation syndrome. *Hum Reprod Update* 2008. — [PMID 18385260](https://pubmed.ncbi.nlm.nih.gov/18385260/)
46. ★ Golan A, et al. Ovarian hyperstimulation syndrome: an update review. *Obstet Gynecol Surv* 1989. — [PMID 2660037](https://pubmed.ncbi.nlm.nih.gov/2660037/) — 모델이 계산해 내는 등급 정의(Hct 45%/55%)
47. ★ Practice Committee of the American Society for Reproductive Medicine. Prevention of moderate and severe ovarian hyperstimulation syndrome: a guideline. *Fertil Steril* 2024. — [PMID 38099867](https://pubmed.ncbi.nlm.nih.gov/38099867/)
48. ★ Papanikolaou EG, et al. Incidence and prediction of ovarian hyperstimulation syndrome in women undergoing gonadotropin-releasing hormone antagonist in vitro fertilization cycles. *Fertil Steril* 2006. — [PMID 16412740](https://pubmed.ncbi.nlm.nih.gov/16412740/) — 난포 수 기반 위험 예측
49. ★ Steward RG, et al. Oocyte number as a predictor for ovarian hyperstimulation syndrome and live birth: an analysis of 256,381 in vitro fertilization cycles. *Fertil Steril* 2014. — [PMID 24462057](https://pubmed.ncbi.nlm.nih.gov/24462057/) — 난자 수–OHSS–출생률의 동시 관계. 모델의 D 스윕과 대응
50. ★ Devroey P, Polyzos NP, Blockeel C. An OHSS-Free Clinic by segmentation of IVF treatment. *Hum Reprod* 2011. — [PMID 21828116](https://pubmed.ncbi.nlm.nih.gov/21828116/) — 작용제 트리거 + 전동결 arm
51. ★ Alvarez C, et al. Dopamine agonist cabergoline reduces hemoconcentration and ascites in hyperstimulated women undergoing assisted reproduction. *J Clin Endocrinol Metab* 2007. — [PMID 17456571](https://pubmed.ncbi.nlm.nih.gov/17456571/) — 카베르골린 arm의 보정 표적
52. Korhonen KV, et al. C-reactive protein response is higher in early than in late ovarian hyperstimulation syndrome. *Eur J Obstet Gynecol Reprod Biol* 2016. — [PMID 27865939](https://pubmed.ncbi.nlm.nih.gov/27865939/) — early/late 구분
53. Youssef MA, et al. Volume expanders for the prevention of ovarian hyperstimulation syndrome. *Cochrane Database Syst Rev* 2016. — [PMID 27577848](https://pubmed.ncbi.nlm.nih.gov/27577848/)
54. Bassiouny YA, et al. Randomized trial of combined cabergoline and coasting in preventing ovarian hyperstimulation syndrome. *Int J Gynaecol Obstet* 2018. — [PMID 29055130](https://pubmed.ncbi.nlm.nih.gov/29055130/) — coasting arm
55. Wu D, et al. Comparison of the effectiveness of various medicines in the prevention of ovarian hyperstimulation syndrome. *Front Endocrinol (Lausanne)* 2022. — [PMID 35154015](https://pubmed.ncbi.nlm.nih.gov/35154015/)
56. Salari E, et al. Comparative study of cabergoline and hydroxychloroquine to prevent ovarian hyperstimulation syndrome. *J Ovarian Res* 2025. — [PMID 40442814](https://pubmed.ncbi.nlm.nih.gov/40442814/)
57. Timmons D, et al. Ovarian hyperstimulation syndrome: a review for emergency clinicians. *Am J Emerg Med* 2019. — [PMID 31097257](https://pubmed.ncbi.nlm.nih.gov/31097257/) — 중증 합병증(혈전·ARDS)
58. Braat DD, et al. Maternal death related to IVF in the Netherlands 1984–2008. *Hum Reprod* 2010. — [PMID 20488805](https://pubmed.ncbi.nlm.nih.gov/20488805/)

## 9. 트리거일 프로게스테론 상승과 전동결

모델에서 트리거일 P4는 별도의 사건이 아니라 과립세포 질량의 네 번째 읽기입니다.

59. ★ Venetis CA, et al. Progesterone elevation and probability of pregnancy after IVF: a systematic review and meta-analysis of over 60,000 cycles. *Hum Reprod Update* 2013. — [PMID 23827986](https://pubmed.ncbi.nlm.nih.gov/23827986/) — 1.5 ng/mL 역치
60. ★ Bosch E, et al. Circulating progesterone levels and ongoing pregnancy rates in controlled ovarian stimulation cycles for in vitro fertilization: analysis of over 4000 cycles. *Hum Reprod* 2010. — [PMID 20539042](https://pubmed.ncbi.nlm.nih.gov/20539042/)
61. Labarta E, et al. Endometrial receptivity is affected in women with high circulating progesterone levels at the end of the follicular phase: a functional genomics analysis. *Hum Reprod* 2011. — [PMID 21540246](https://pubmed.ncbi.nlm.nih.gov/21540246/)
62. Shapiro BS, et al. Evidence of impaired endometrial receptivity after ovarian stimulation for in vitro fertilization: a prospective randomized trial. *Fertil Steril* 2011. — [PMID 21737072](https://pubmed.ncbi.nlm.nih.gov/21737072/)
63. Roque M, et al. Fresh versus elective frozen embryo transfer in IVF/ICSI cycles: a systematic review and meta-analysis of reproductive outcomes. *Hum Reprod Update* 2019. — [PMID 30388233](https://pubmed.ncbi.nlm.nih.gov/30388233/)

## 10. 용량–반응과 개별화 투여

모델의 A 스윕(용량 3배 → 난자 +15%, 복수 +125%)이 대응하는 임상 문헌입니다.

64. ★ Nyboe Andersen A, et al. Individualized versus conventional ovarian stimulation for in vitro fertilization: a multicenter, randomized, controlled, assessor-blinded, phase 3 noninferiority trial (ESTHER-1). *Fertil Steril* 2017. — [PMID 27912901](https://pubmed.ncbi.nlm.nih.gov/27912901/) — 표준 arm(AFC 12 · 150 IU · 난자 ~10)의 보정 표적
65. Nelson SM, et al. Individualized versus conventional ovarian stimulation for in vitro fertilization (authors' reply). *Fertil Steril* 2024. — [PMID 38266801](https://pubmed.ncbi.nlm.nih.gov/38266801/)
66. Ngwenya O, et al. Individualised gonadotropin dose selection using markers of ovarian reserve for women undergoing IVF/ICSI. *Cochrane Database Syst Rev* 2024. — [PMID 38174816](https://pubmed.ncbi.nlm.nih.gov/38174816/)
67. van Tilborg TC, et al. Individualized versus standard FSH dosing in women starting IVF/ICSI: an RCT. Part 1: The predicted poor responder (OPTIMIST). *Hum Reprod* 2017. — [PMID 29121326](https://pubmed.ncbi.nlm.nih.gov/29121326/) — 용량 증량이 결과를 바꾸지 못한다는 결과
68. Arce JC, et al. Antimüllerian hormone in gonadotropin releasing-hormone antagonist cycles: prediction of ovarian response and cumulative treatment outcome in good-prognosis patients. *Fertil Steril* 2013. — [PMID 23394782](https://pubmed.ncbi.nlm.nih.gov/23394782/)
69. Blockeel C, et al. Follicular phase endocrine characteristics during ovarian stimulation and GnRH antagonist cotreatment for IVF. *J Clin Endocrinol Metab* 2011. — [PMID 21307142](https://pubmed.ncbi.nlm.nih.gov/21307142/) — 자극 중 내인성 FSH·LH·E2 궤적

## 11. 난자 수 · 배수성 · 누적 생존출생률

모델의 배아 사슬(2PN → 배아낭 → 정상배수체 → CLBR)의 계수는 전부 여기서 왔습니다.

70. ★ Sunkara SK, et al. Association between the number of eggs and live birth in IVF treatment: an analysis of 400,135 treatment cycles. *Hum Reprod* 2011. — [PMID 21558332](https://pubmed.ncbi.nlm.nih.gov/21558332/) — 난자 수–출생률의 포화 곡선
71. ★ Drakopoulos P, et al. Conventional ovarian stimulation and single embryo transfer for IVF/ICSI. How many oocytes do we need to maximize cumulative live birth rates after utilization of all fresh and frozen embryos? *Hum Reprod* 2016. — [PMID 26724797](https://pubmed.ncbi.nlm.nih.gov/26724797/) — 모델 D 스윕의 대응 문헌
72. ★ Franasiak JM, et al. The nature of aneuploidy with increasing age of the female partner: a review of 15,169 consecutive trophectoderm biopsies evaluated with comprehensive chromosomal screening. *Fertil Steril* 2014. — [PMID 24355045](https://pubmed.ncbi.nlm.nih.gov/24355045/) — 정상배수체 분율 곡선 `0.85/(1+exp((AGE−38.2)/4))`
73. Goldman RH, et al. Predicting the likelihood of live birth for elective oocyte cryopreservation: a counseling tool for physicians and patients. *Hum Reprod* 2017. — [PMID 28166330](https://pubmed.ncbi.nlm.nih.gov/28166330/) — 난자→출생 사슬의 독립적 구성
74. Toftager M, et al. Cumulative live birth rates after one ART cycle including all subsequent frozen-thaw cycles in 1050 women. *Hum Reprod* 2017. — [PMID 28130435](https://pubmed.ncbi.nlm.nih.gov/28130435/)
75. Totonchi M, et al. Preimplantation genetic screening and the success rate of in vitro fertilization: a three-years study on Iranian population. *Cell J* 2021. — [PMID 32347040](https://pubmed.ncbi.nlm.nih.gov/32347040/)

## 12. PCOS — 고반응군의 기전

모델은 PCOS를 "고반응군"으로 지정하지 않습니다. AMH가 임계값을 올리고 LH 긴장도가
LHCGR 획득 전 난포를 정지시키는 두 항만으로 6–9 mm 정지(다낭성 형태)와 자극 시의
과다반응이 동시에 나옵니다.

76. ★ Jonard S, Dewailly D. The follicular excess in polycystic ovaries, due to intra-ovarian hyperandrogenism, may be the main culprit for the follicular arrest. *Hum Reprod Update* 2004. — [PMID 15073141](https://pubmed.ncbi.nlm.nih.gov/15073141/) — 난포 정지 기전
77. ★ Willis DS, et al. Premature response to luteinizing hormone of granulosa cells from anovulatory women with polycystic ovary syndrome: relevance to mechanism of anovulation. *J Clin Endocrinol Metab* 1998. — [PMID 9814480](https://pubmed.ncbi.nlm.nih.gov/9814480/) — 모델의 `KLHARR`·`NLHARR`

## 13. 저반응군 · LH 보충

78. ★ Ferraretti AP, et al. ESHRE consensus on the definition of 'poor response' to ovarian stimulation for in vitro fertilization: the Bologna criteria. *Hum Reprod* 2011. — [PMID 21505041](https://pubmed.ncbi.nlm.nih.gov/21505041/)
79. ★ Humaidan P, et al. The novel POSEIDON stratification of 'Low prognosis patients in Assisted Reproductive Technology' and its proposed marker of successful outcome. *F1000Res* 2016. — [PMID 28232864](https://pubmed.ncbi.nlm.nih.gov/28232864/) — 저반응군 arm(AFC 5 · T50 12)
80. ★ Mochtar MH, et al. Recombinant luteinizing hormone (rLH) and recombinant follicle stimulating hormone (rFSH) for ovarian stimulation in IVF/ICSI cycles. *Cochrane Database Syst Rev* 2017. — [PMID 28537052](https://pubmed.ncbi.nlm.nih.gov/28537052/) — LH 보충이 난자 수를 늘리지 못한다는 모델 예측의 대응 문헌

## 14. 황체기 지지

81. Fatemi HM, et al. An update of luteal phase support in stimulated IVF cycles. *Hum Reprod Update* 2007. — [PMID 17626114](https://pubmed.ncbi.nlm.nih.gov/17626114/)
82. Fatemi HM. Simplifying luteal phase support in stimulated assisted reproduction cycles. *Fertil Steril* 2018. — [PMID 30396544](https://pubmed.ncbi.nlm.nih.gov/30396544/)
83. ★ Tavaniotou A, Devroey P. Effect of human chorionic gonadotropin on luteal luteinizing hormone concentrations in natural cycles. *Fertil Steril* 2003. — [PMID 12969719](https://pubmed.ncbi.nlm.nih.gov/12969719/) — 프로게스테론–LH 음성 되먹임(`KP4LH`)의 직접 근거

## 15. 월경주기·난포 동역학의 수학적 모델 (선행 모델)

이 모델은 아래 계열의 연장선에 있으며, 차이는 (i) 난포를 임계값 분위 슬롯으로
이산화한 점, (ii) LHCGR 신호를 세 개의 커널로 나눈 점, (iii) OHSS 체액 구획을
포함한 점입니다.

84. ★ Reinecke I, Deuflhard P. A complex mathematical model of the human menstrual cycle. *J Theor Biol* 2007. — [PMID 17448501](https://pubmed.ncbi.nlm.nih.gov/17448501/)
85. ★ Röblitz S, et al. A mathematical model of the human menstrual cycle for the administration of GnRH analogues. *J Theor Biol* 2013. — [PMID 23206386](https://pubmed.ncbi.nlm.nih.gov/23206386/) — GnRH 유사체 투여를 포함한 선행 모델
86. Röblitz S, et al. A computational systems biology view on the role of the menstrual cycle. *J Mol Endocrinol* 2026. — [PMID 42306999](https://pubmed.ncbi.nlm.nih.gov/42306999/)
87. ★ Clark LR, Schlosser PM, Selgrade JF. Multiple stable periodic solutions in a model for hormonal control of the menstrual cycle. *Bull Math Biol* 2003. — [PMID 12597121](https://pubmed.ncbi.nlm.nih.gov/12597121/) — 주기 모델의 쌍안정성

---

## 인용 요약

| 구분 | 편수 |
|------|------|
| ★ 파라미터 보정에 직접 사용 | 48 |
| 배경·기전·검증 | 39 |
| **합계 (전부 PMID 개별 확인)** | **87** |

## 모델이 문헌과 어긋나는 지점 (의도적으로 남긴 것)

정직하게 기록해 둡니다. 아래는 모델의 구조가 만들어 낸 예측이며, 반증 가능한
형태로 남겨 둔 것입니다.

1. **황체기 레트로졸**: 모델에서 E2를 79% 낮추지만 복수는 4% 미만 변화합니다.
   E2를 표지자, VEGF를 매개자로 본 구조의 직접적 결과입니다. 레트로졸이 OHSS를
   줄인다는 무작위 자료가 확립되면 E2 의존 투과성 항의 누락이 반증됩니다.
2. **작용제 트리거 후 자발 배란**: 파열 적분(`RUPX`)에 도달하지 못하므로 44시간
   에서도 배란이 거의 없다고 예측합니다. hCG보다 관대한 채취 시각이 허용된다는
   뜻이며, 임상적으로 확인되지 않았습니다.
3. **coasting**: 소난포 과립세포 질량을 죽이는 방식으로만 작동하므로 난자 손실
   (−35%) 없이 OHSS를 줄일 수 없습니다.
4. **작용제 트리거의 황체기 결핍**: 방향은 재현하나 크기(−41%)가 임상보다 작습니다
   (참고문헌 29·30). 누락된 기전은 황체화 자체의 질적 결함으로 추정합니다.
5. **카베르골린**: 메타분석의 RR 0.38보다 보수적입니다(복수 −17%).
6. **자연 주기 우성난포 직경**: 16.9 mm에서 멈춥니다(관찰 20–22 mm). 서지를
   발화시키는 E2 농도가 직경보다 먼저 도달하기 때문입니다.
