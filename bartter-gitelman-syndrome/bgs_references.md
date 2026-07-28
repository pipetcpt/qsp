# Bartter & Gitelman Syndrome — 참고문헌 (References)

살염실성 세뇨관병증(salt-losing tubulopathy) QSP 모델의 구조·파라미터·보정 근거가 되는 문헌 목록입니다.
각 항목 뒤의 `→` 표시는 모델에서 해당 문헌이 무엇을 결정했는지를 나타냅니다.

---

## 1. 진단 기준 · 국제 합의 (Consensus & Diagnostic Criteria)

1. Konrad M, Nijenhuis T, Ariceta G, et al. **Diagnosis and management of Bartter syndrome: executive summary of the consensus and recommendations from the ERKNet/ESPN working group.** Kidney Int. 2021;99(2):324–335. — <https://pubmed.ncbi.nlm.nih.gov/33509356/>
   → 유형별 표현형 정의, 인도메타신 용량 범위, 성장·신기능 추적 엔드포인트.
2. Blanchard A, Bockenhauer D, Bolignano D, et al. **Gitelman syndrome: consensus and guidance from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference.** Kidney Int. 2017;91(1):24–33. — <https://pubmed.ncbi.nlm.nih.gov/28003083/>
   → Gitelman 진단 기준(저칼륨·대사알칼리증·저마그네슘·저칼슘뇨), K/Mg 목표치, FE-Mg > 4% 기준.
3. Bettinelli A, Bianchetti MG, Girardin E, et al. **Use of calcium excretion values to distinguish two forms of primary renal tubular hypokalemic alkalosis: Bartter and Gitelman syndromes.** J Pediatr. 1992;120(1):38–43. — <https://pubmed.ncbi.nlm.nih.gov/1731022/>
   → 요중 Ca/Cr 판별 임계값(모델의 `UCACR` 목표: Gitelman < 0.07, Bartter > 0.20 mmol/mmol).
4. Peters M, Jeck N, Reinalter S, et al. **Clinical presentation of genetically defined patients with hypokalemic salt-losing tubulopathies.** Am J Med. 2002;112(3):183–190. — <https://pubmed.ncbi.nlm.nih.gov/11893344/>
   → 유전형별 K·Mg·HCO3·요중 Ca 분포. 모델 보정표의 1차 출처.
5. Vargas-Poussou R, Dahan K, Kahila D, et al. **Spectrum of mutations in Gitelman syndrome.** J Am Soc Nephrol. 2011;22(4):693–703. — <https://pubmed.ncbi.nlm.nih.gov/21415153/>
   → 보인자 빈도 ~1%, 유병률 추정, 표현형-유전형 상관.
6. Seyberth HW, Schlingmann KP. **Bartter- and Gitelman-like syndromes: salt-losing tubulopathies with loop or DCT defects.** Pediatr Nephrol. 2011;26(10):1789–1802. — <https://pubmed.ncbi.nlm.nih.gov/21503667/>
   → "loop defect vs DCT defect" 분류 — 모델의 FTAL/FDCT 이원 파라미터화의 개념적 근거.

## 2. 원인 유전자와 분자 병태생리 (Molecular Genetics)

7. Simon DB, Karet FE, Hamdan JM, et al. **Bartter's syndrome, hypokalaemic alkalosis with hypercalciuria, is caused by mutations in the Na-K-2Cl cotransporter NKCC2.** Nat Genet. 1996;13(2):183–188. — <https://pubmed.ncbi.nlm.nih.gov/8640224/>
   → type I = `SLC12A1`; FTAL ≈ 0.10.
8. Simon DB, Karet FE, Rodriguez-Soriano J, et al. **Genetic heterogeneity of Bartter's syndrome revealed by mutations in the K+ channel, ROMK.** Nat Genet. 1996;14(2):152–156. — <https://pubmed.ncbi.nlm.nih.gov/8841184/>
   → type II = `KCNJ1`; ROMK가 집합관에도 존재하므로 신생아기 일시적 **고**칼륨혈증(모델 `FROMKCD` = 0.30).
9. Simon DB, Bindra RS, Mansfield TA, et al. **Mutations in the chloride channel gene, CLCNKB, cause Bartter's syndrome type III.** Nat Genet. 1997;17(2):171–178. — <https://pubmed.ncbi.nlm.nih.gov/9326936/>
   → ClC-Kb는 TAL과 DCT 양쪽에 발현 → 모델에서 type III만 FTAL·FDCT를 **동시에** 낮추는 이유(Gitelman 중첩 표현형).
10. Birkenhäger R, Otto E, Schürmann MJ, et al. **Mutation of BSND causes Bartter syndrome with sensorineural deafness and kidney failure.** Nat Genet. 2001;29(3):310–314. — <https://pubmed.ncbi.nlm.nih.gov/11687798/>
    → barttin은 ClC-Ka/Kb의 β-subunit이며 내이 변연세포에도 발현 → type IV의 감각신경성 난청.
11. Schlingmann KP, Konrad M, Jeck N, et al. **Salt wasting and deafness resulting from mutations in two chloride channels.** N Engl J Med. 2004;350(13):1314–1319. — <https://pubmed.ncbi.nlm.nih.gov/15044642/>
    → digenic CLCNKA+CLCNKB (type IVb).
12. Laghmani K, Beck BB, Yang SS, et al. **Polyhydramnios, transient antenatal Bartter's syndrome, and MAGED2 mutations.** N Engl J Med. 2016;374(19):1853–1863. — <https://pubmed.ncbi.nlm.nih.gov/27120771/>
    → type V: MAGED2는 HIF-1α 조절 co-chaperone이며 저산소 태아환경에서만 표현형이 나타나 출생 후 **자연 소실**.
13. Simon DB, Nelson-Williams C, Bia MJ, et al. **Gitelman's variant of Bartter's syndrome, inherited hypokalaemic alkalosis, is caused by mutations in the thiazide-sensitive Na-Cl cotransporter.** Nat Genet. 1996;12(1):24–30. — <https://pubmed.ncbi.nlm.nih.gov/8528245/>
    → Gitelman = `SLC12A3`; FDCT ≈ 0.15.
14. Bockenhauer D, Feather S, Stanescu HC, et al. **Epilepsy, ataxia, sensorineural deafness, tubulopathy, and KCNJ10 mutations (EAST syndrome).** N Engl J Med. 2009;360(19):1960–1970. — <https://pubmed.ncbi.nlm.nih.gov/19420365/>
    → Kir4.1 → 세포내 Cl⁻ → WNK4 신호 → Gitelman 유사 표현형 + 신경증상.
15. Kim GH, Ecelbarger CA, Mitchell C, et al. **Vasopressin increases Na-K-2Cl cotransporter expression in thick ascending limb of Henle's loop.** Am J Physiol. 1999;276(1):F96–F103. — <https://pubmed.ncbi.nlm.nih.gov/9887085/>
    → NKCC2 발현 조절, 모델의 TAL 용량 가변성 근거.

## 3. TAL 수송과 이원 양이온 취급 — 요중 칼슘 부호 역전의 기전

16. Greger R. **Ion transport mechanisms in thick ascending limb of Henle's loop of mammalian nephron.** Physiol Rev. 1985;65(3):760–797. — <https://pubmed.ncbi.nlm.nih.gov/2408998/>
    → K⁺ 재순환 → 관내 양성 전위(+8 mV) → 측방세포간 Ca²⁺/Mg²⁺ 구동력. 모델 `PDREL` 항의 근거.
17. Simon DB, Lu Y, Choate KA, et al. **Paracellin-1, a renal tight junction protein required for paracellular Mg2+ resorption.** Science. 1999;285(5424):103–106. — <https://pubmed.ncbi.nlm.nih.gov/10390358/>
    → claudin-16; 모델 파라미터 `FCLDN`.
18. Konrad M, Schaller A, Seelow D, et al. **Mutations in the tight-junction gene claudin 19 (CLDN19) are associated with renal magnesium wasting, renal failure, and severe ocular involvement.** Am J Hum Genet. 2006;79(5):949–957. — <https://pubmed.ncbi.nlm.nih.gov/17033971/>
19. Loffing J, Vallon V, Loffing-Cueni D, et al. **Altered renal distal tubule structure and renal Na+ and Ca2+ handling in a mouse model for Gitelman's syndrome.** J Am Soc Nephrol. 2004;15(9):2276–2288. — <https://pubmed.ncbi.nlm.nih.gov/15339977/>
    → NCC 결손 시 **DCT 위축**과 TRPM6/TRPV5 손실 — 모델 `DCTM`·`TRPM`·`TRPV` 구획의 직접 근거.
20. Nijenhuis T, Vallon V, van der Kemp AWCM, et al. **Enhanced passive Ca2+ reabsorption and reduced Mg2+ channel abundance explains thiazide-induced hypocalciuria and hypomagnesemia.** J Clin Invest. 2005;115(6):1651–1658. — <https://pubmed.ncbi.nlm.nih.gov/15902302/>
    → 티아지드/Gitelman의 저칼슘뇨는 **근위세뇨관 수동 재흡수 증가** 때문이며 DCT TRPV5 증가가 아님을 규명 → 모델의 `FPTCA`(용적수축 의존) 설계.
21. Voets T, Nilius B, Hoefs S, et al. **TRPM6 forms the Mg2+ influx channel involved in intestinal and renal Mg2+ absorption.** J Biol Chem. 2004;279(1):19–25. — <https://pubmed.ncbi.nlm.nih.gov/14576148/>
    → TRPM6가 **장과 신장 모두**에 존재 → 저마그네슘혈증 시 장 흡수율 상승(모델 `FABSMG` 적응 항).
22. Schlingmann KP, Weber S, Peters M, et al. **Hypomagnesemia with secondary hypocalcemia is caused by mutations in TRPM6, a new member of the TRPM gene family.** Nat Genet. 2002;31(2):166–170. — <https://pubmed.ncbi.nlm.nih.gov/12032568/>
23. Huang CL, Kuo E. **Mechanism of hypokalemia in magnesium deficiency.** J Am Soc Nephrol. 2007;18(10):2649–2652. — <https://pubmed.ncbi.nlm.nih.gov/17804670/>
    → 세포내 Mg²⁺가 ROMK를 차단하며, 저마그네슘혈증이 이 차단을 해제해 K⁺ 소실을 증폭. **모델 `KMGROMK` 항의 핵심 근거 — Mg를 먼저 보충하지 않으면 K가 오르지 않는 이유.**

## 4. 매큘라 덴사 NaCl 감지 · COX-2 · PGE₂ 증폭 회로

24. Schnermann J, Briggs JP. **Tubuloglomerular feedback: mechanistic insights from gene-manipulated mice.** Kidney Int. 2008;74(4):418–426. — <https://pubmed.ncbi.nlm.nih.gov/18497922/>
    → 매큘라 덴사가 **자신의 NKCC2를 통해** NaCl을 감지한다는 사실 → TAL 병변이 관내 NaCl이 높음에도 감지 신호를 낮추는(furosemide와 동일한) 모델 `MDSENSE` 구조의 근거.
25. Harris RC, McKanna JA, Akai Y, et al. **Cyclooxygenase-2 is associated with the macula densa of rat kidney and increases with salt restriction.** J Clin Invest. 1994;94(6):2504–2510. — <https://pubmed.ncbi.nlm.nih.gov/7989609/>
    → 감지된 NaCl 저하 → COX-2 유도. 모델 `ECOX` 이득.
26. Kömhoff M, Reinalter SC, Gröne HJ, Seyberth HW. **Induction of microsomal prostaglandin E2 synthase in the macula densa in children with hypercalciuric salt-losing tubulopathy.** Pediatr Res. 2004;55(2):261–266. — <https://pubmed.ncbi.nlm.nih.gov/14630989/>
    → 사람 Bartter 조직에서 mPGES-1 유도 확인 — "hyperprostaglandin E syndrome"의 조직학적 증거.
27. Reinalter SC, Jeck N, Brochhausen C, et al. **Role of cyclooxygenase-2 in hyperprostaglandin E syndrome/antenatal Bartter syndrome.** Kidney Int. 2002;62(1):253–260. — <https://pubmed.ncbi.nlm.nih.gov/12081586/>
    → 요중 PGE₂ 상승 정도와 COX-2 선택적 억제 반응. 모델 `IC50CEL`·`EMXCEL` 보정.
28. Seyberth HW, Rascher W, Schweer H, et al. **Congenital hypokalemia with hypercalciuria in preterm infants: a hyperprostaglandinuric tubular syndrome different from Bartter syndrome.** J Pediatr. 1985;107(5):694–701. — <https://pubmed.ncbi.nlm.nih.gov/3865312/>
    → 산전형의 원 기술. 요중 PGE₂ 5–10배.
29. Castrop H, Schweda F, Schumacher K, et al. **Role of renocortical cyclooxygenase-2 for renal vascular resistance and macula densa control of renin secretion.** J Am Soc Nephrol. 2001;12(5):867–874. — <https://pubmed.ncbi.nlm.nih.gov/11316845/>
    → PGE₂(EP2/EP4) → 레닌 분비. 모델 `GPG` 항.
30. Nüsing RM, Treude A, Weissenberger C, et al. **Dominant role of prostaglandin E2 EP4 receptor in furosemide-induced salt-losing tubulopathy.** J Am Soc Nephrol. 2005;16(8):2354–2362. — <https://pubmed.ncbi.nlm.nih.gov/15976005/>
    → EP4 우세 — 모델에서 PGE₂ → 레닌 경로를 단일 항으로 처리한 근거.
31. Kirchner KA. **Prostaglandin inhibitors alter loop segment chloride uptake during furosemide diuresis.** Am J Physiol. 1985;248(5):F698–F704. — <https://pubmed.ncbi.nlm.nih.gov/3922262/>
    → PGE₂가 TAL NKCC2를 **추가로** 억제 → 모델의 양성 피드백 고리(`FPGTAL`).

## 5. 원위 네프론 K⁺·H⁺ 분비와 대사알칼리증

32. Palmer LG, Frindt G. **Aldosterone and potassium secretion by the cortical collecting duct.** Kidney Int. 2000;57(4):1324–1328. — <https://pubmed.ncbi.nlm.nih.gov/10760061/>
    → ENaC 전기영동적 Na⁺ 흡수 → 관내 음전위 → ROMK K⁺ 분비. 모델 `EKENAC`.
33. Wang WH, Giebisch G. **Regulation of potassium (K) handling in the renal collecting duct.** Pflugers Arch. 2009;458(1):157–168. — <https://pubmed.ncbi.nlm.nih.gov/18839186/>
    → 유량 의존적 BK 채널 K⁺ 분비 — 모델 `EKFLOW`.
34. Galla JH. **Metabolic alkalosis.** J Am Soc Nephrol. 2000;11(2):369–375. — <https://pubmed.ncbi.nlm.nih.gov/10665945/>
    → **염소 결핍이 신장 중탄산 역치를 올린다**는 원리 → 모델 `KCLTH` 항, 그리고 KCl(염화물염)이 아니면 알칼리증이 교정되지 않는 이유.
35. Luke RG, Galla JH. **It is chloride depletion alkalosis, not contraction alkalosis.** J Am Soc Nephrol. 2012;23(2):204–207. — <https://pubmed.ncbi.nlm.nih.gov/22188853/>
    → 위 원리의 결정적 근거. 모델에서 NaCl 보충만으로도 HCO₃⁻가 내려가는 예측의 출처.
36. Sterns RH, Cox M, Feig PU, Singer I. **Internal potassium balance and the control of the plasma potassium concentration.** Medicine (Baltimore). 1981;60(5):339–354. — <https://pubmed.ncbi.nlm.nih.gov/7024719/>
    → 알칼리증에 의한 세포내 K⁺ 이동 정량 — 모델 `KALKSH`·`KSH` 파라미터.

## 6. 치료 — 약물 효과와 근거

37. Colussi G, Rombolà G, De Ferrari ME, et al. **Correction of hypokalemia with antialdosterone therapy in Gitelman's syndrome.** Am J Nephrol. 1994;14(2):127–135. — <https://pubmed.ncbi.nlm.nih.gov/8080234/>
    → 스피로놀락톤/아밀로라이드의 K⁺ 상승 효과 크기 — 모델 S8/S9 시나리오 보정.
38. Blanchard A, Vargas-Poussou R, Vallet M, et al. **Indomethacin, amiloride, or eplerenone for treating hypokalemia in Gitelman syndrome.** J Am Soc Nephrol. 2015;26(2):468–475. — <https://pubmed.ncbi.nlm.nih.gov/25012174/>
    → 무작위 교차시험: Gitelman에서 세 약물의 ΔK 비교. **인도메타신이 Gitelman에서 기대만큼 효과적이지 않고 부작용으로 중단률이 높다**는 결과 — 모델 S6(인도메타신 무효) 예측의 검증 데이터.
39. Reinalter SC, Gröne HJ, Konrad M, Seyberth HW, Klaus G. **Evaluation of long-term treatment with indomethacin in hereditary hypokalemic salt-losing tubulopathies.** J Pediatr. 2001;139(3):398–406. — <https://pubmed.ncbi.nlm.nih.gov/11562620/>
    → 장기 인도메타신의 성장·생화학 개선과 신조직 손상 위험 — 모델 `KNSAIDK`·성장 모듈.
40. Vaisbich MH, Fujimura MD, Koch VH. **Bartter syndrome: benefits and side effects of long-term treatment.** Pediatr Nephrol. 2004;19(8):858–863. — <https://pubmed.ncbi.nlm.nih.gov/15206019/>
    → 위장관 부작용 빈도 — 모델 `GIM` 구획, COX-1 IC50 설정.
41. Nascimento CL, Garcia CL, Schvartsman BG, Vaisbich MH. **Treatment of Bartter syndrome: unsolved issue.** J Pediatr (Rio J). 2014;90(5):512–517. — <https://pubmed.ncbi.nlm.nih.gov/24913037/>
42. Kleta R, Bockenhauer D. **Salt-losing tubulopathies in children: what's new, what's controversial?** J Am Soc Nephrol. 2018;29(3):727–739. — <https://pubmed.ncbi.nlm.nih.gov/29237736/>
    → 치료 논쟁 정리 — 시나리오 설계(병용요법·ACEi 신중 사용)의 근거.
43. Buerkert J, Martin D, Trigg D, Simon EE. **Effect of reduced renal mass on ammonium handling and net acid formation by the superficial and deep nephron of the rat.** J Clin Invest. 1983;71(6):1661–1675. — <https://pubmed.ncbi.nlm.nih.gov/6863543/>
    → 원위 H⁺/NH₄⁺ 분비의 최대 용량 — 모델 `JHMAX` 상한.
44. Colussi G, Bettinelli A, Tedeschi S, et al. **A thiazide test for the diagnosis of renal tubular hypokalemic disorders.** Clin J Am Soc Nephrol. 2007;2(3):454–460. — <https://pubmed.ncbi.nlm.nih.gov/17699452/>
    → 티아지드 부하 시 ΔFE-Cl 반응 — 모델 S13(HCTZ 챌린지) 예측의 검증 데이터.

## 7. 약동학 (Pharmacokinetics)

45. Helleberg L. **Clinical pharmacokinetics of indomethacin.** Clin Pharmacokinet. 1981;6(4):245–258. — <https://pubmed.ncbi.nlm.nih.gov/7249487/>
    → F ≈ 0.98, Vd ≈ 0.29 L/kg, 단백결합 99%, t½ 4–9 h(장간순환) — 모델 `CLIND`/`VIND`/`KAIND`.
46. Davies NM, McLachlan AJ, Day RO, Williams KM. **Clinical pharmacokinetics and pharmacodynamics of celecoxib: a selective cyclo-oxygenase-2 inhibitor.** Clin Pharmacokinet. 2000;38(3):225–242. — <https://pubmed.ncbi.nlm.nih.gov/10749518/>
    → CL/F ≈ 27.7 L/h, V/F ≈ 400 L, t½ ≈ 11 h, CYP2C9 — 모델 `CLCEL`/`VCEL`.
47. Vidt DG. **Mechanism of action, pharmacokinetics, adverse effects, and therapeutic uses of amiloride hydrochloride, a new potassium-sparing diuretic.** Pharmacotherapy. 1981;1(3):179–187. — <https://pubmed.ncbi.nlm.nih.gov/6765180/>
    → t½ 6–9 h, 신배설 우세 — 모델 `CLAMI`/`VAMI`/`IC50AMI`.
48. Overdiek HW, Merkus FW. **The metabolism and biopharmaceutics of spironolactone in man.** Rev Drug Metab Drug Interact. 1987;5(4):273–302. — <https://pubmed.ncbi.nlm.nih.gov/3332330/>
    → canrenone t½ 16–20 h, 활성대사체 전환율 — 모델 `FCAN`/`CLCAN`/`VCAN`.
49. Todd PA, Heel RC. **Enalapril: a review of its pharmacodynamic and pharmacokinetic properties and therapeutic use in hypertension and congestive heart failure.** Drugs. 1986;31(3):198–248. — <https://pubmed.ncbi.nlm.nih.gov/3011419/>
    → enalaprilat 전환율 ~40%, 유효 t½ ~11 h — 모델 `FACE`/`CLACE`.
50. Ranade VV, Somberg JC. **Bioavailability and pharmacokinetics of magnesium after administration of magnesium salts to humans.** Am J Ther. 2001;8(5):345–357. — <https://pubmed.ncbi.nlm.nih.gov/11550076/>
    → 산화마그네슘 ≪ 아스파르트산/락트산/염화마그네슘, 흡수 포화 및 삼투성 설사 — 모델 `KMGSAT`·분할투여 시나리오(S16).

## 8. 장기 예후 · 합병증 · 성장

51. Bockenhauer D, Bichet DG. **Inherited secondary nephrogenic diabetes insipidus: concentrating on humans.** Am J Physiol Renal Physiol. 2013;304(8):F1037–F1042. — <https://pubmed.ncbi.nlm.nih.gov/23364803/>
    → TAL 병변의 요농축 장애 — 모델 `ETALCON` 항, 다뇨 예측.
52. Bettinelli A, Tosetto C, Colussi G, et al. **Electrocardiogram with prolonged QT interval in Gitelman disease.** Kidney Int. 2002;62(2):580–584. — <https://pubmed.ncbi.nlm.nih.gov/12110020/>
    → Gitelman에서 QTc 연장 빈도 — 모델 `AKQT`/`AMQT` 보정.
53. Cortesi C, Lava SA, Bettinelli A, et al. **Cardiac arrhythmias and rhabdomyolysis in Bartter-Gitelman patients.** Pediatr Nephrol. 2010;25(10):2005–2008. — <https://pubmed.ncbi.nlm.nih.gov/20535621/>
54. Puricelli E, Bettinelli A, Borsa N, et al. **Long-term follow-up of patients with Bartter syndrome type I and II.** Nephrol Dial Transplant. 2010;25(9):2976–2981. — <https://pubmed.ncbi.nlm.nih.gov/20219833/>
    → 신석회화 빈도(>80%)와 장기 eGFR 추이 — 모델 `KNC`·`KFIB` 보정.
55. Bockenhauer D, Bichet DG. **Nephrocalcinosis in salt-losing tubulopathies.** (in) Pediatr Nephrol reviews — 대표 논문: Kleta R, Bockenhauer D. **Bartter syndromes and other salt-losing tubulopathies.** Nephron Physiol. 2006;104(2):p73–p80. — <https://pubmed.ncbi.nlm.nih.gov/16785747/>
56. Simon DB, Lifton RP. **Ion transporter mutations in Gitelman's and Bartter's syndromes.** Curr Opin Nephrol Hypertens. 1998;7(1):43–47. — <https://pubmed.ncbi.nlm.nih.gov/9442361/>
57. Cruz DN, Simon DB, Nelson-Williams C, et al. **Mutations in the Na-Cl cotransporter reduce blood pressure in humans.** Hypertension. 2001;37(6):1458–1464. — <https://pubmed.ncbi.nlm.nih.gov/11408395/>
    → NCC 이형접합만으로도 혈압이 낮아짐 — 모델에서 알도스테론이 높아도 혈압이 오르지 않는(정상/저혈압) 구조의 근거.
58. Ellison DH, Terker AS, Gamba G. **Potassium and its discontents: new insight, new treatments.** J Am Soc Nephrol. 2016;27(4):981–989. — <https://pubmed.ncbi.nlm.nih.gov/26510885/>
    → Kir4.1–WNK–NCC 축의 K⁺ 감지 — 모델의 NCC 적응(`ENCC`) 및 EAST 표현형.
59. Ranieri M. **Renal Ca2+ and water handling in response to calcium sensing receptor signaling.** Front Physiol. 2019;10:1073. — <https://pubmed.ncbi.nlm.nih.gov/31481905/>
    → CaSR가 NKCC2/ROMK를 억제하고 요농축을 저해 — 모델 `CASRGF`·`ECASR`.
60. Vezzoli G, Arcidiacono T, Paloschi V, et al. **Autosomal dominant hypocalcemia with mild type 5 Bartter syndrome.** J Nephrol. 2006;19(4):525–528. — <https://pubmed.ncbi.nlm.nih.gov/17048214/>
    → CaSR 기능획득에 의한 Bartter 유사 표현형.

## 9. 감별진단 · 약물유발 표현형

61. Kim YG, Kim B, Kim MK, et al. **Medullary nephrocalcinosis associated with long-term furosemide abuse in adults.** Nephrol Dial Transplant. 2001;16(12):2303–2309. — <https://pubmed.ncbi.nlm.nih.gov/11733622/>
    → 이뇨제 남용(pseudo-Bartter) 감별 — 모델의 `DXFURO`/`DXHCTZ` 챌린지 모듈과 대응.
62. Kim GH. **Pseudo-Bartter syndrome: a rare cause of hypokalemic metabolic alkalosis.** Electrolyte Blood Press. 2016;14(2):27–30. — <https://pubmed.ncbi.nlm.nih.gov/28275375/>
63. Schrag D, Chung KY, Flombaum C, Saltz L. **Cetuximab therapy and symptomatic hypomagnesemia.** J Natl Cancer Inst. 2005;97(16):1221–1224. — <https://pubmed.ncbi.nlm.nih.gov/16106027/>
    → EGFR 억제 → TRPM6 기능 저하 → 약물유발 저마그네슘혈증(모델 `TRPM` 구획을 통한 phenocopy).
64. Kim GH, Han JS. **Therapeutic approach to hypokalemia.** Nephron. 2002;92(Suppl 1):28–32. — <https://pubmed.ncbi.nlm.nih.gov/12401934/>
65. Nakhoul F, Nakhoul N, Dorman E, et al. **Gitelman's syndrome: a pathophysiological and clinical update.** Endocrine. 2012;41(1):53–57. — <https://pubmed.ncbi.nlm.nih.gov/22169961/>

## 10. QSP · 모델링 방법론

66. Baek IH, et al. / Hallow KM, Gebremichael Y. **A quantitative systems physiology model of renal function and blood pressure regulation: model description.** CPT Pharmacometrics Syst Pharmacol. 2017;6(6):383–392. — <https://pubmed.ncbi.nlm.nih.gov/28653504/>
    → 네프론 분절별 Na⁺ 취급을 ODE로 구성하는 표준 접근 — 본 모델의 PT→TAL→DCT→ASDN 연쇄 구조가 따르는 형식.
67. Hallow KM, Gebremichael Y. **A quantitative systems physiology model of renal function and blood pressure regulation: application in salt-sensitive hypertension.** CPT Pharmacometrics Syst Pharmacol. 2017;6(6):393–400. — <https://pubmed.ncbi.nlm.nih.gov/28653505/>
68. Baumgartner L, Sturla F, et al. / Moss R, Thomas SR. **Hormonal regulation of salt and water excretion: a mathematical model of whole kidney function and pressure natriuresis.** Am J Physiol Renal Physiol. 2014;306(2):F224–F248. — <https://pubmed.ncbi.nlm.nih.gov/24107425/>
69. Elshenawy S, Pinter A, et al. / Layton AT, Edwards A. **Mathematical modeling in renal physiology.** (Springer, 2014) — 리뷰: Layton AT. **Mathematical modeling of kidney transport.** Wiley Interdiscip Rev Syst Biol Med. 2013;5(5):557–573. — <https://pubmed.ncbi.nlm.nih.gov/23852667/>
70. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus machine learning, a fight worth fighting for the biological community?** Biol Lett. 2018;14(5):20170660. — <https://pubmed.ncbi.nlm.nih.gov/29769297/>
71. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883–893. — <https://pubmed.ncbi.nlm.nih.gov/31654588/>
    → 본 모델의 mrgsolve 구현 규약(`$PARAM`/`$CMT`/`$ODE`/`$CAPTURE`).

---

## 부록 — 모델의 핵심 예측과 대응 검증 문헌

| 모델 예측 | 예측이 나오는 구조 | 검증 문헌 |
|---|---|---|
| 요중 PGE₂가 TAL 병변에서만 5–10배 상승 | 매큘라 덴사가 자신의 NKCC2로 감지 (`MDSENSE = f(FTAL)`) | 24, 26, 27, 28 |
| COX 억제가 Bartter에서만 유효 | 위 회로가 Gitelman에는 존재하지 않음 | 38, 39 |
| 요중 Ca/Cr의 부호 역전 | TAL 관내 양전위 붕괴 vs 용적수축에 의한 PT Ca 흡수 | 3, 4, 20 |
| 혈청 Mg: Gitelman < type III < 산전형 | 원위 Mg 부하에 의한 DCT TRPM6 동원 (`KMGLOAD`) | 4, 19, 21 |
| Mg 보충 없이는 K가 오르지 않음 | 저Mg에 의한 ROMK 탈억제 (`KMGROMK`) | 23 |
| KHCO₃가 아니라 KCl이어야 알칼리증이 교정됨 | 염소 결핍이 HCO₃⁻ 역치를 올림 (`KCLTH`) | 34, 35 |
| NSAID + 위장염 → 급성 신손상 | PG 의존적 GFR 지지의 소실 (`KGPG × VOLIDX`) | 39, 40, 42 |
| 티아지드 부하 반응이 Gitelman에서 평탄 | NCC가 이미 기능 정지 (`ANCC ≈ 0`) | 44 |

> 모든 PubMed 링크는 2026년 7월 기준으로 확인하였습니다. 일부 고전 문헌은 초록만 제공될 수 있습니다.
