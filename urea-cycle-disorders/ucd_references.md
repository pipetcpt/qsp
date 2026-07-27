# Urea Cycle Disorders (UCD) — QSP 모델 참고문헌
### References for the UCD / OTC-deficiency Quantitative Systems Pharmacology model

이 문서는 `ucd_qsp_model.dot` (기계론적 지도), `ucd_mrgsolve_model.R` (ODE 모델),
`ucd_shiny_app.R` (대시보드)에 들어간 **구조·파라미터·임상 앵커의 근거**를 정리한
것입니다. 각 절의 첫머리에 그 문헌군이 모델의 **어느 부분을 규정했는지** 적었습니다.

---

## 1. 종설 · 진료지침 (Reviews & guidelines)

> 모델의 전체 구조, 진단 알고리즘(시트룰린 + 오로트산 판별), 암모니아 문턱값
> (100 / 200 / 360 µmol/L), 응급 처치 순서를 규정.

1. Häberle J, et al. **Suggested guidelines for the diagnosis and management of urea cycle disorders: first revision.** *J Inherit Metab Dis* 2019;42:1192-1230. — 개정 유럽 지침. 암모니아 문턱, 투석 적응증, 스캐빈저 용량의 1차 근거. <https://pubmed.ncbi.nlm.nih.gov/30982989/>
2. Häberle J, et al. **Suggested guidelines for the diagnosis and management of urea cycle disorders.** *Orphanet J Rare Dis* 2012;7:32. <https://pubmed.ncbi.nlm.nih.gov/22642880/>
3. Summar ML, et al. **The incidence of urea cycle disorders.** *Mol Genet Metab* 2013;110:179-180. <https://pubmed.ncbi.nlm.nih.gov/23972786/>
4. Brusilow SW, Maestri NE. **Urea cycle disorders: diagnosis, pathophysiology, and therapy.** *Adv Pediatr* 1996;43:127-170. — 대체 질소 배출 경로 개념의 원전. <https://pubmed.ncbi.nlm.nih.gov/8794176/>
5. Matsumoto S, et al. **Urea cycle disorders — update.** *J Hum Genet* 2019;64:833-847. <https://pubmed.ncbi.nlm.nih.gov/31110235/>
6. Ah Mew N, et al. **Urea Cycle Disorders Overview.** *GeneReviews* (updated 2017). <https://pubmed.ncbi.nlm.nih.gov/20301396/>
7. Batshaw ML, Tuchman M, Summar M, Seminara J; Members of the UCDC. **A longitudinal study of urea cycle disorders.** *Mol Genet Metab* 2014;113:127-130. — Urea Cycle Disorders Consortium 종단 코호트. <https://pubmed.ncbi.nlm.nih.gov/25135652/>
8. Stepien KM, et al. **Challenges in the management of adult patients with urea cycle disorders.** *J Clin Med* 2022;11:2029. <https://pubmed.ncbi.nlm.nih.gov/35407637/>

## 2. 요소회로 효소학 · NAG-CPS1 조절

> `VMAX_UC`(2400 µmol N/kg/h), `KMNH4`(120 µmol/L), NAG 활성화 항 `f_NAG`,
> 아르기닌의 NAGS 자극(`GARG`), 카글룸산의 CPS1 활성화(`GNCG`)를 규정.

9. Meijer AJ, Lamers WH, Chamuleau RA. **Nitrogen metabolism and ornithine cycle function.** *Physiol Rev* 1990;70:701-748. — 간 문맥주위/중심정맥주위 구획화와 요소회로 조절의 표준 종설. <https://pubmed.ncbi.nlm.nih.gov/2194222/>
10. Morris SM Jr. **Regulation of enzymes of the urea cycle and arginine metabolism.** *Annu Rev Nutr* 2002;22:87-105. <https://pubmed.ncbi.nlm.nih.gov/12055339/>
11. Caldovic L, Tuchman M. **N-acetylglutamate and its changing role through evolution.** *Biochem J* 2003;372:279-290. <https://pubmed.ncbi.nlm.nih.gov/12633501/>
12. Shi D, Allewell NM, Tuchman M. **The N-acetylglutamate synthase family: structures, function and mechanisms.** *Int J Mol Sci* 2015;16:13004-13022. <https://pubmed.ncbi.nlm.nih.gov/26068232/>
13. Nissim I, et al. **Effects of a glucokinase activator on hepatic intermediary metabolism: study with 13C-isotopomer-based metabolomics.** *Biochem J* 2012;444:537-551. — 간 요소생성 플럭스 정량. <https://pubmed.ncbi.nlm.nih.gov/22448977/>
14. Rubio V, Grisolía S. **Human carbamoylphosphate synthetase I.** *Enzyme* 1981;26:233-239. <https://pubmed.ncbi.nlm.nih.gov/7318135/>
15. Yudkoff M, et al. **In vivo measurement of ureagenesis with stable isotopes.** *J Inherit Metab Dis* 1998;21(Suppl 1):21-29. — 15N/13C 요소생성 플럭스 측정법; 유전자치료 시험의 1차 약력학 지표. <https://pubmed.ncbi.nlm.nih.gov/9686341/>

## 3. OTC 결핍 — 유전학, X-불활성화, 표현형 스펙트럼

> `OTCACT` 파라미터의 의미, 이형접합 여성의 간 모자이시즘(`Xinact` 노드),
> 신생아형 vs 후기발병형 구분(`Onset_class`)의 근거.

16. Tuchman M, et al. **Cross-sectional multicenter study of patients with urea cycle disorders in the United States.** *Mol Genet Metab* 2008;94:397-402. <https://pubmed.ncbi.nlm.nih.gov/18562231/>
17. Yamaguchi S, et al. **Mutations and polymorphisms in the human ornithine transcarbamylase (OTC) gene.** *Hum Mutat* 2006;27:626-632. <https://pubmed.ncbi.nlm.nih.gov/16786505/>
18. Maestri NE, Brusilow SW, Clissold DB, Bassett SS. **Long-term treatment of girls with ornithine transcarbamylase deficiency.** *N Engl J Med* 1996;335:855-859. <https://pubmed.ncbi.nlm.nih.gov/8778603/>
19. Batshaw ML, Msall M, Beaudet AL, Trojak J. **Risk of serious illness in heterozygotes for ornithine transcarbamylase deficiency.** *J Pediatr* 1986;108:236-241. <https://pubmed.ncbi.nlm.nih.gov/3944709/>
20. Gyato K, Wray J, Huang ZJ, Yudkoff M, Batshaw ML. **Metabolic and neuropsychological phenotype in women heterozygous for ornithine transcarbamylase deficiency.** *Ann Neurol* 2004;55:80-86. — '무증상' 보인자에서도 실행기능 결함. 모델의 ADHD-유사 표현형 노드. <https://pubmed.ncbi.nlm.nih.gov/14705115/>
21. Lichter-Konecki U, Caldovic L, Morizono H, Simpson K. **Ornithine Transcarbamylase Deficiency.** *GeneReviews* (updated 2022). <https://pubmed.ncbi.nlm.nih.gov/24006547/>

## 4. 암모니아 뇌독성 — 별아교세포 글루타민 삼투질 가설

> 모델에서 가장 중요한 CNS 모듈(`GLNB` · `MINSB` · `BRWAT` · `ICP`)의 근거.
> 만성 고암모니아혈증에서 myo-inositol이 이미 고갈되어 있어 급성 상승 시
> 완충이 없다는 `FCOMP` 상호작용이 여기서 나옵니다.

22. Brusilow SW, Koehler RC, Traystman RJ, Cooper AJL. **Astrocyte glutamine synthetase: importance in hyperammonemic syndromes and potential target for therapy.** *Neurotherapeutics* 2010;7:452-470. — 삼투질 가설의 결정판. <https://pubmed.ncbi.nlm.nih.gov/20880508/>
23. Butterworth RF. **Pathophysiology of hepatic encephalopathy: a new look at ammonia.** *Metab Brain Dis* 2002;17:221-227. <https://pubmed.ncbi.nlm.nih.gov/12602499/>
24. Norenberg MD, Rao KV, Jayakumar AR. **Mechanisms of ammonia-induced astrocyte swelling.** *Metab Brain Dis* 2005;20:303-318. — 미토콘드리아 투과성 전이(mPTP), ROS/RNS 경로. <https://pubmed.ncbi.nlm.nih.gov/16382342/>
25. Häussinger D, Laubenberger J, vom Dahl S, et al. **Proton magnetic resonance spectroscopy studies on human brain myo-inositol in hypo-osmolarity and hepatic encephalopathy.** *Gastroenterology* 1994;107:1475-1480. — 글루타민↑ / myo-inositol↓ 상호 교환의 직접 증거. <https://pubmed.ncbi.nlm.nih.gov/7926510/>
26. Connelly A, Cross JH, Gadian DG, Hunter JV, Kirkham FJ, Leonard JV. **Magnetic resonance spectroscopy shows increased brain glutamine in ornithine carbamoyl transferase deficiency.** *Pediatr Res* 1993;33:77-81. <https://pubmed.ncbi.nlm.nih.gov/8433886/>
27. Gropman AL, et al. **1H MRS identifies symptomatic and asymptomatic subjects with partial ornithine transcarbamylase deficiency.** *Mol Genet Metab* 2008;95:21-30. <https://pubmed.ncbi.nlm.nih.gov/18662894/>
28. Gropman AL, Summar M, Leonard JV. **Neurological implications of urea cycle disorders.** *J Inherit Metab Dis* 2007;30:865-879. <https://pubmed.ncbi.nlm.nih.gov/18038189/>
29. Cooper AJL, Plum F. **Biochemistry and physiology of brain ammonia.** *Physiol Rev* 1987;67:440-519. — NH3/NH4+ 분배, 이온 트래핑, 뇌 pH 의존성. `PKANH3`·`PHBR` 파라미터의 근거. <https://pubmed.ncbi.nlm.nih.gov/2882529/>
30. Bachmann C. **Mechanisms of hyperammonemia.** *Clin Chem Lab Med* 2002;40:653-662. <https://pubmed.ncbi.nlm.nih.gov/12241009/>
31. Cagnon L, Braissant O. **Hyperammonemia-induced toxicity for the developing central nervous system.** *Brain Res Rev* 2007;56:183-197. <https://pubmed.ncbi.nlm.nih.gov/17881060/>

## 5. 신경학적 예후 — 혼수 지속시간이 결정한다

> `NEURO` 누적 손상 적분기와 `IQEST` 출력의 근거. 모델은 손상을
> `(암모니아 − 150) × 시간 × (1 + 2×뇌수분)`으로 적분합니다.

32. Msall M, Batshaw ML, Suss R, Brusilow SW, Mellits ED. **Neurologic outcome in children with inborn errors of urea synthesis: outcome of urea-cycle enzymopathies.** *N Engl J Med* 1984;310:1500-1505. — 신생아 고암모니아혈증 **혼수 지속시간**이 이후 IQ의 최강 예측인자임을 처음 보인 논문. <https://pubmed.ncbi.nlm.nih.gov/6717538/>
33. Bachmann C. **Outcome and survival of 88 patients with urea cycle disorders: a retrospective evaluation.** *Eur J Pediatr* 2003;162:410-416. <https://pubmed.ncbi.nlm.nih.gov/12684900/>
34. Krivitzky L, et al. **Intellectual, adaptive, and behavioral functioning in children with urea cycle disorders.** *Pediatr Res* 2009;66:96-101. <https://pubmed.ncbi.nlm.nih.gov/19287347/>
35. Enns GM, et al. **Survival after treatment with phenylacetate and benzoate for urea-cycle disorders.** *N Engl J Med* 2007;356:2282-2292. — Ammonul 등록 데이터; 혼수 지속과 생존의 관계. <https://pubmed.ncbi.nlm.nih.gov/17538087/>
36. Waisbren SE, et al. **Neuropsychological outcomes in individuals with urea cycle disorders.** *Mol Genet Metab* 2016;119:37-45. <https://pubmed.ncbi.nlm.nih.gov/27380995/>
37. Posset R, et al. **Impact of diagnosis and therapy on cognitive function in urea cycle disorders.** *Ann Neurol* 2019;86:116-128. <https://pubmed.ncbi.nlm.nih.gov/31018026/>

## 6. 질소 스캐빈저 — 페닐부티레이트 축의 PK/PD

> `KHYDG`(GPB 가수분해), `CLPBA`·`FMPAA`(β-산화), `VMPAGN`·`KMPAA`·`KMGLNP`
> (GLYATL1 포합), `CLPGN`(OAT 분비), 그리고 NaPBA vs GPB 프로파일 차이를 규정.

38. Brusilow SW, Danney M, Waber LJ, et al. **Treatment of episodic hyperammonemia in children with inborn errors of urea synthesis.** *N Engl J Med* 1984;310:1630-1634. — 벤조산·페닐아세트산 요법의 원전. <https://pubmed.ncbi.nlm.nih.gov/6427608/>
39. Brusilow SW. **Phenylacetylglutamine may replace urea as a vehicle for waste nitrogen excretion.** *Pediatr Res* 1991;29:147-150. — PAGN이 분자당 **질소 2개**를 나른다는 핵심 화학량론. <https://pubmed.ncbi.nlm.nih.gov/1903526/>
40. Diaz GA, Krivitzky LS, Mokhtarani M, et al. **Ammonia control and neurocognitive outcome among urea cycle disorder patients treated with glycerol phenylbutyrate.** *Hepatology* 2013;57:2171-2179. — GPB vs NaPBA 교차 설계. **동일 PBA 몰수에서 GPB는 Cmax가 낮고 24시간 암모니아 프로파일이 평탄**. 모델 시나리오 3 vs 4의 앵커. <https://pubmed.ncbi.nlm.nih.gov/22961727/>
41. Smith W, Diaz GA, Lichter-Konecki U, et al. **Ammonia control in children ages 2 months through 5 years with urea cycle disorders: comparison of sodium phenylbutyrate and glycerol phenylbutyrate.** *J Pediatr* 2013;162:1228-1234. <https://pubmed.ncbi.nlm.nih.gov/23324524/>
42. Monteleone JPR, Mokhtarani M, Diaz GA, et al. **Population pharmacokinetic modeling and dosing simulations of nitrogen-scavenging compounds: disposition of glycerol phenylbutyrate and sodium phenylbutyrate in adult and pediatric patients with urea cycle disorders.** *J Clin Pharmacol* 2013;53:699-710. — 이 모델의 스캐빈저 PK 구조(전구체 → PBA → PAA → PAGN)와 U-PAGN 회수율(~60-75%)의 직접 근거. <https://pubmed.ncbi.nlm.nih.gov/23686462/>
43. Mokhtarani M, Diaz GA, Rhead W, et al. **Urinary phenylacetylglutamine as dosing biomarker for patients with urea cycle disorders.** *Mol Genet Metab* 2012;107:308-314. — U-PAGN을 용량 적정 바이오마커로 쓰는 근거. <https://pubmed.ncbi.nlm.nih.gov/22958974/>
44. Mokhtarani M, Diaz GA, Rhead W, et al. **Elevated phenylacetic acid levels do not correlate with adverse events in patients with urea cycle disorders or hepatic encephalopathy and can be predicted based on the plasma PAA to PAGN ratio.** *Mol Genet Metab* 2013;110:446-453. — **PAA:PAGN 비 > 2.5 (µg/mL 기준) = 포합 포화**. 모델의 `PARATIO` 출력과 시나리오 13의 근거. <https://pubmed.ncbi.nlm.nih.gov/24144944/>
45. Thibault A, et al. **Phase I study of phenylacetate administered twice daily to patients with cancer.** *Cancer* 1995;75:2932-2938. — PAA 신경독성(졸림·오심)이 나타나는 혈중 농도 범위(≈500 µg/mL). <https://pubmed.ncbi.nlm.nih.gov/7773945/>
46. Berry SA, et al. **Safety and efficacy of glycerol phenylbutyrate for management of urea cycle disorders in patients aged 2 months to 2 years.** *Mol Genet Metab* 2017;122:46-53. <https://pubmed.ncbi.nlm.nih.gov/28916119/>
47. Scaglia F, Carter S, O'Brien WE, Lee B. **Effect of alternative pathway therapy on branched chain amino acid metabolism in urea cycle disorder patients.** *Mol Genet Metab* 2004;81(Suppl 1):S79-S85. — 페닐부티레이트 특이적 **BCAA 고갈**. 모델의 `BCAAP` 구획과 `KBPAA` 파라미터. <https://pubmed.ncbi.nlm.nih.gov/15050979/>
48. Burrage LC, et al. **Sodium phenylbutyrate decreases plasma branched-chain amino acids in patients with urea cycle disorders.** *Mol Genet Metab* 2014;113:131-135. <https://pubmed.ncbi.nlm.nih.gov/25042691/>

## 7. 벤조산-글리신 축, 카글룸산, 회로 기질

> `VMHIP`·`KMBZ`·`KMGLY`(GLYAT 포합, **질소 1개**), 글리신 고갈,
> 카글룸산의 NAGS 대체, 시트룰린의 OTC-우회 질소 배출(`KCITN`)을 규정.

49. Batshaw ML, MacArthur RB, Tuchman M. **Alternative pathway therapy for urea cycle disorders: twenty years later.** *J Pediatr* 2001;138(1 Suppl):S46-S55. — 벤조산 vs 페닐부티레이트의 질소 화학량론 비교. <https://pubmed.ncbi.nlm.nih.gov/11148549/>
50. Kasapkara ÇS, et al. **N-carbamylglutamate treatment for acute neonatal hyperammonemia in isolated methylmalonic acidemia.** *Eur J Pediatr* 2011;170:799-801. <https://pubmed.ncbi.nlm.nih.gov/21170724/>
51. Daniotti M, la Marca G, Fiorini P, Filippi L. **New developments in the treatment of hyperammonemia: emerging use of carglumic acid.** *Int J Gen Med* 2011;4:21-28. <https://pubmed.ncbi.nlm.nih.gov/21403788/>
52. Caldovic L, et al. **N-acetylglutamate synthase deficiency: a heterogeneous disorder.** *Mol Genet Metab* 2007;91:36-41. — NAGS 결핍이 카글룸산으로 **표현형적으로 치유**되는 유일한 UCD임을 확립. 모델 시나리오 11. <https://pubmed.ncbi.nlm.nih.gov/17368065/>
53. Brusilow SW. **Arginine, an indispensable amino acid for patients with inborn errors of urea synthesis.** *J Clin Invest* 1984;74:2144-2148. <https://pubmed.ncbi.nlm.nih.gov/6511918/>
54. Nagasaka H, et al. **Effects of arginine treatment on nutrition, growth and urea cycle function in seven Japanese boys with late-onset ornithine transcarbamylase deficiency.** *Eur J Pediatr* 2006;165:618-624. <https://pubmed.ncbi.nlm.nih.gov/16691408/>
55. Erez A, Nagamani SCS, Shchelochkov OA, et al. **Requirement of argininosuccinate lyase for systemic nitric oxide production.** *Nat Med* 2011;17:1619-1626. — ASL이 NOS 복합체의 필수 구성요소라는 발견. 모델의 `NO_deficit`·`Hypertension_ASL` 노드. <https://pubmed.ncbi.nlm.nih.gov/22081021/>
56. Nagamani SCS, et al. **Nitric-oxide supplementation for treatment of long-term complications in argininosuccinic aciduria.** *Am J Hum Genet* 2012;90:836-846. <https://pubmed.ncbi.nlm.nih.gov/22541557/>

## 8. 글루타민 — 완충고이자 조기경보

> `ucd_vgs()` 함수(고친화 + 저친화 GS 성분), Gln:NH3 관계,
> 그리고 글루타민이 암모니아보다 **먼저** 오른다는 임상 관찰의 근거.

57. Darmaun D, Matthews DE, Bier DM. **Glutamine and glutamate kinetics in humans.** *Am J Physiol* 1986;251:E117-E126. — 전신 글루타민 플럭스 ≈ 300 µmol/kg/h. `KGSKG`·`VG0KG` 보정의 직접 근거. <https://pubmed.ncbi.nlm.nih.gov/2873746/>
58. Lee B, Diaz GA, Rhead W, et al. **Blood ammonia and glutamine as predictors of hyperammonemic crises in patients with urea cycle disorder.** *Genet Med* 2015;17:561-568. — **글루타민 > 1000 µmol/L가 위기를 예측**. 모델의 조기경보 논리. <https://pubmed.ncbi.nlm.nih.gov/25341094/>
59. Maestri NE, McGowan KD, Brusilow SW. **Plasma glutamine concentration: a guide in the management of urea cycle disorders.** *J Pediatr* 1992;121:259-261. <https://pubmed.ncbi.nlm.nih.gov/1640295/>
60. Häussinger D. **Nitrogen metabolism in liver: structural and functional organization and physiological relevance.** *Biochem J* 1990;267:281-290. — 문맥주위(고용량 요소생성) / 중심정맥주위(고친화 글루타민 합성) 구획화. 모델의 `Hepatocyte_periportal` / `Hepatocyte_perivenous` 노드. <https://pubmed.ncbi.nlm.nih.gov/1970241/>
61. Weiner ID, Verlander JW. **Renal ammonia metabolism and transport.** *Compr Physiol* 2013;3:201-220. — 신장 글루타미나아제 유래 암모니아 배설(모델의 `CLRGLN` 경로). <https://pubmed.ncbi.nlm.nih.gov/23720285/>

## 9. 오로트산 — OTC를 CPS1/NAGS와 가르는 신호

62. Bachmann C, Colombo JP. **Diagnostic value of orotic acid excretion in heritable disorders of the urea cycle and in hyperammonemia due to organic acidurias.** *Eur J Pediatr* 1980;134:109-113. <https://pubmed.ncbi.nlm.nih.gov/7449803/>
63. Salerno C, Crifò C. **Diagnostic value of urinary orotic acid levels: applicable separation methods.** *J Chromatogr B* 2002;781:57-71. <https://pubmed.ncbi.nlm.nih.gov/12450653/>

## 10. 응급 처치 · 체외 제거 · 반동

> `CLHD`(투석 청소율), `FHDGLN`, 그리고 심부 글루타민 저장고(`GLNM`)에서
> 오는 **투석 후 반동**의 근거. 시나리오 9 vs 10.

64. Picca S, et al. **Extracorporeal dialysis in neonatal hyperammonemia: modalities and prognostic indicators.** *Pediatr Nephrol* 2001;16:862-867. <https://pubmed.ncbi.nlm.nih.gov/11685590/>
65. Spinale JM, et al. **High-dose continuous renal replacement therapy for neonatal hyperammonemia.** *Pediatr Nephrol* 2013;28:983-986. <https://pubmed.ncbi.nlm.nih.gov/23515666/>
66. Schaefer F, et al. **Dialysis in neonates with inborn errors of metabolism.** *Nephrol Dial Transplant* 1999;14:910-918. <https://pubmed.ncbi.nlm.nih.gov/10328466/>
67. Wiegand C, Thompson T, Bock GH, Mathis RK, Kjellstrand CM, Mauer SM. **The management of life-threatening hyperammonemia: a comparison of several therapeutic modalities.** *J Pediatr* 1980;96:142-144. — 혈액투석이 복막투석·교환수혈보다 압도적으로 빠름. <https://pubmed.ncbi.nlm.nih.gov/7350297/>
68. Ah Mew N, et al. **Comparison of sodium phenylbutyrate and glycerol phenylbutyrate in the acute management of hyperammonemia.** *Mol Genet Metab* 2018;124:1-6. <https://pubmed.ncbi.nlm.nih.gov/29655841/>

## 11. 간이식 · 유전자치료 · mRNA 치료

> `ACTX`/`TTX`(이식), `TRANSG`·`OTCX`·`KEXPR`(AAV/ mRNA)의 근거.
> 시나리오 12.

69. Yu L, et al. **Liver transplantation for urea cycle disorders: analysis of the United Network for Organ Sharing database.** *Transplant Proc* 2015;47:2413-2418. <https://pubmed.ncbi.nlm.nih.gov/26518941/>
70. Kido J, et al. **Long-term outcome and intervention of urea cycle disorders in Japan.** *J Inherit Metab Dis* 2012;35:777-785. — 이식은 요소생성을 표현형적으로 치유하지만 **기존 CNS 손상은 되돌리지 못한다**. <https://pubmed.ncbi.nlm.nih.gov/22167275/>
71. Raper SE, et al. **Fatal systemic inflammatory response syndrome in a ornithine transcarbamylase deficient patient following adenoviral gene transfer.** *Mol Genet Metab* 2003;80:148-158. — 유전자치료 역사의 전환점(Gelsinger 사례). 모델의 `Immune_AAV` 노드. <https://pubmed.ncbi.nlm.nih.gov/14567964/>
72. Wang L, et al. **AAV gene therapy corrects OTC deficiency and prevents liver fibrosis in aged neonatal ornithine transcarbamylase-deficient mice.** *Mol Genet Metab* 2017;120:299-305. <https://pubmed.ncbi.nlm.nih.gov/28202336/>
73. Harding CO, et al. **Safety and efficacy of DTX301, an AAV8-mediated gene transfer for adults with late-onset ornithine transcarbamylase deficiency: interim results.** *Mol Genet Metab* 2021;132:S49. — 요소생성 플럭스 회복과 스캐빈저 감량. (프로그램 개요) <https://pubmed.ncbi.nlm.nih.gov/33642232/>
74. Prieve MG, et al. **Targeted mRNA therapy for ornithine transcarbamylase deficiency.** *Mol Ther* 2018;26:801-813. — LNP-mRNA로 반복 투여하며 간 OTC 활성을 펄스처럼 회복. 모델의 mRNA 모드. <https://pubmed.ncbi.nlm.nih.gov/29433939/>
75. Diez-Fernandez C, Häberle J. **Targeting CPS1 in the treatment of carbamoyl phosphate synthetase 1 deficiency.** *Expert Opin Ther Targets* 2017;21:391-399. <https://pubmed.ncbi.nlm.nih.gov/28281899/>

## 12. 유발인자 · 약물상호작용 (발프로산 등)

> 시나리오 14와 지도의 `Valproate`·`Corticosteroid_DDI`·`Probenecid` 노드.

76. Coulter DL, Allen RJ. **Secondary hyperammonaemia: a possible mechanism for valproate encephalopathy.** *Lancet* 1980;1:1310-1311. <https://pubmed.ncbi.nlm.nih.gov/6104119/>
77. Aires CCP, et al. **New insights on the mechanisms of valproate-induced hyperammonemia: inhibition of hepatic N-acetylglutamate synthase activity by valproyl-CoA.** *J Hepatol* 2011;55:426-434. — 발프로산이 NAGS를 직접 억제. 모델의 `VPAI` 파라미터. <https://pubmed.ncbi.nlm.nih.gov/21147182/>
78. Tuchman M, Yudkoff M. **Blood levels of ammonia and nitrogen scavenging amino acids in patients with inherited hyperammonemia.** *Mol Genet Metab* 1999;66:10-15. <https://pubmed.ncbi.nlm.nih.gov/9973543/>
79. Nott L, Price TJ, Pittman K, Patterson K, Fletcher J. **Hyperammonemia encephalopathy: an important cause of neurological deterioration following chemotherapy.** *Leuk Lymphoma* 2007;48:1702-1711. <https://pubmed.ncbi.nlm.nih.gov/17786705/>
80. Lipskind S, Loanzon S, Simi E, Ouyang DW. **Hyperammonemic coma in ornithine transcarbamylase deficiency: a case of postpartum decompensation.** *Obstet Gynecol* 2011;117:503-505. — 산후 자궁퇴축이 대표적 이화 유발인자. <https://pubmed.ncbi.nlm.nih.gov/21252804/>

## 13. 영양 · 단백 처방 · 성장

> `PROT`·`FDISP`·`NOBLKG` 파라미터와 `PROTTOL`(천연 단백 내성) 출력의 근거.

81. Adam S, et al. **Dietary management of urea cycle disorders: European practice.** *Mol Genet Metab* 2013;110:439-445. <https://pubmed.ncbi.nlm.nih.gov/24113687/>
82. Singh RH. **Nutritional management of patients with urea cycle disorders.** *J Inherit Metab Dis* 2007;30:880-887. <https://pubmed.ncbi.nlm.nih.gov/17957501/>
83. Boyle M, et al. **Growth in patients with urea cycle disorders.** *Mol Genet Metab* 2014;113:220-224. <https://pubmed.ncbi.nlm.nih.gov/25266922/>
84. WHO/FAO/UNU. **Protein and amino acid requirements in human nutrition.** *WHO Tech Rep Ser* 935, 2007. — 의무적 질소 손실(obligatory nitrogen loss) 값의 출처. <https://pubmed.ncbi.nlm.nih.gov/18330140/>

## 14. QSP · mrgsolve 방법론

85. Baron KT, et al. **mrgsolve: Simulate from ODE-Based Models.** R package. <https://mrgsolve.org/>
86. Nijsen MJMA, et al. **Preclinical QSP modeling in the pharmaceutical industry: an IQ consortium survey.** *CPT Pharmacometrics Syst Pharmacol* 2018;7:135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/>
87. Bai JPF, et al. **Quantitative systems pharmacology: landscape analysis of regulatory submissions to the US FDA.** *CPT Pharmacometrics Syst Pharmacol* 2021;10:1479-1484. <https://pubmed.ncbi.nlm.nih.gov/34617412/>
88. Gadkar K, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2016;5:235-249. <https://pubmed.ncbi.nlm.nih.gov/27299936/>

---

## 부록 A — 모델 파라미터가 어느 문헌에서 왔는가

| 파라미터 | 값 | 근거 |
|----------|-----|------|
| `UCAP` (최대 요소생성 용량) | 2400 µmol N/kg/h | 최대 요소생성률이 습관적 질소 부하의 4-5배 [9, 13, 15] |
| `KMNH4` | 120 µmol/L | 건강인이 암모니아 30 µmol/L에서 용량의 ~20%를 사용하도록 하는 내부 정합 값 [9, 29] |
| `VG0KG` / `KGSKG` / `KGSAT` | 184.6 / 3.843 / 1200 | 전신 글루타민 플럭스 ≈300 µmol/kg/h [57], 관찰된 Gln-NH3 관계 [58, 59] |
| `CLRGKG` | 0.0275 L/h/kg | 신장 글루타민 유래 질소 배설 30-50 mmol/day [61] |
| `VMPKG` / `KMPAA` / `KMGLNP` | 115 µmol/kg/h / 300 / 400 | U-PAGN 회수율 60-75% [42, 43] |
| `KHYDG` (GPB 가수분해) | 0.30 /h | GPB의 낮은 Cmax·평탄한 프로파일 [40, 42] |
| PAA 독성 문턱 | 500 µg/mL | [44, 45] |
| PAA:PAGN 문턱 | 2.5 (µg/mL) | [44] |
| `PSBBB` / `PKANH3` / `PHBR` | 400 L/h / 9.15 / 7.05 | 비이온 확산 + 뇌 이온 트래핑 [29] |
| `VMGSB` / `KOUTGB` | 892 µmol/h / 0.05 /h | 뇌 글루타민 5 → 15-20 mmol/L 상승 [25, 26, 27] |
| `KINMI` / `KEXTMI` / `FCOMP` | 0.10 / 0.02 / 0.60 | myo-inositol의 **부분적** 삼투 보상 [25, 22] |
| `KINJ` / `IQ50` | 0.10 / 50 | 혼수 지속시간-IQ 관계 [32, 33, 37] |
| 암모니아 문턱 100 / 200 / 360 | — | [1] |
| `NOBLKG` (의무적 질소 유출) | 43 µmol N/kg/h | 무단백 식이 시 의무적 질소 손실 [84] |

## 부록 B — 모델이 재현하는 임상 관찰

| 관찰 | 재현 방식 |
|------|-----------|
| 질환은 경사가 아니라 **절벽** | Michaelis-Menten 요소생성 + 글루타민 완충고 → 잔존 활성 5%에서 암모니아 484, 40%에서 61 µmol/L (탭 1의 cliff 그래프) |
| 글루타민이 암모니아보다 **먼저** 오른다 | `ucd_vgs()`가 암모니아에서 글루타민으로 질소를 옮기는 완충 구조 [58] |
| **만성** 고암모니아혈증 환자가 급성 유발에 더 취약하다 | 기저에서 이미 고갈된 `MINSB` → 삼투 완충 없음 [25] |
| GPB는 같은 몰수에서 **Cmax가 낮다** | `GPBG → PBAG` 리파아제 의존 가수분해 단계 [40] |
| 벤조산은 몰당 질소 1개, 페닐부티레이트는 2개 | `dxdt_NH4C`의 `-vhip` (1 N) vs PAGN 경로 (2 N) [39, 49] |
| 신생아 위기에서 **투석이 약물보다 결정적** | 시나리오 9 vs 10 (암모니아 1182 → 236 vs 1182 → 31) [67] |
| **투석 후 반동** | 심부 글루타민 저장고 `GLNM`의 재분포 (시나리오 10, 86-120 h) |
| 카글룸산은 NAGS 결핍을 **표현형적으로 치유** | `GNCG`가 `f_NAG`를 정상화 → 단백 내성 0.46 → 1.35 g/kg/d (시나리오 11) [52] |
| 예후를 정하는 것은 치료 방식이 아니라 **진단까지의 지연** | 시나리오 9/10에서 60시간 지연 후에는 투석을 해도 IQ 회복이 3점에 그침 [32] |

---

> **면책 (Disclaimer).** 본 모델과 문헌 정리는 교육·연구 목적입니다.
> 파라미터는 공개 문헌에서 유도한 근사치이며 독립 검증을 거치지 않았습니다.
> 임상 의사결정·처방·규제 제출에 사용해서는 안 됩니다.
