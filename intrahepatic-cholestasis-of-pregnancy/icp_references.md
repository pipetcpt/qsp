# 임신성 간내 담즙정체 (ICP) — QSP 모델 참고문헌

> **이 목록의 작성 방식.** 아래 모든 PMID 는 NCBI E-utilities (`esearch` +
> `esummary`) 로 실제 조회하여 **반환된 레코드만** 옮긴 것입니다. 제1저자·연도·
> 저널명은 PubMed 가 돌려준 값 그대로이며, 기억이나 추정으로 적은 항목은
> 없습니다. 각 항목 뒤의 `→` 는 그 문헌이 모델의 **어느 파라미터 또는 어느
> 구조적 결정**에 쓰였는지를 가리킵니다. 근거가 없는 값은 "모델 선택"이라고
> 명시했습니다.
>
> **How this list was built.** Every PMID below was resolved against NCBI
> E-utilities and only records that came back are listed; first author, year
> and journal are as PubMed returned them. The `→` note on each entry says
> which parameter or which structural decision in the model it supports.
> Where the model has no literature anchor, that is stated instead of
> concealed.

---

## 1. 진단 기준 · 진료지침 · 역학 (Diagnosis, guidelines, epidemiology)

모델의 중증도 층화(≥10 진단 · ≥40 중증 · ≥100 최중증)와 분만 시기 비교의
기준점.

1. Dajti E, 2025, *Nat Rev Dis Primers* — Intrahepatic cholestasis of pregnancy. [PMID 40707479](https://pubmed.ncbi.nlm.nih.gov/40707479/) → 질환 전반의 최신 종합; 진단 문턱·역학·치료 현황의 기준 서술
2. European Association for the Study of the Liver, 2023, *J Hepatol* — EASL Clinical Practice Guidelines on the management of liver diseases in pregnancy. [PMID 37394016](https://pubmed.ncbi.nlm.nih.gov/37394016/) → 진단 기준과 분만 시기 권고 (모델 G 절의 비교 대상)
3. Society for Maternal-Fetal Medicine (SMFM), 2021, *Am J Obstet Gynecol* — Consult Series #56: Intrahepatic cholestasis of pregnancy. [PMID 33197417](https://pubmed.ncbi.nlm.nih.gov/33197417/) → ≥100 µmol/L 에서 36–37주 분만 권고; 모델 최적화 결과와 대조
4. Hobson SR, 2024, *J Obstet Gynaecol Can* — Guideline No. 452: Diagnosis and Management of Intrahepatic Cholestasis of Pregnancy. [PMID 39089469](https://pubmed.ncbi.nlm.nih.gov/39089469/) → 캐나다 지침, 층화별 관리
5. Hague WM, 2023, *Aust N Z J Obstet Gynaecol* — Diagnosis and management: consensus statement (SOMANZ). [PMID 37431680](https://pubmed.ncbi.nlm.nih.gov/37431680/) → 층화 기준
6. Kothari S, 2024, *Gastroenterology* — AGA Clinical Practice Update on Pregnancy-Related Gastrointestinal and Liver Disease. [PMID 39140906](https://pubmed.ncbi.nlm.nih.gov/39140906/) → 소화기 관점의 관리 권고
7. Giouleka S, 2025, *Obstet Gynecol Surv* — Intrahepatic Cholestasis of Pregnancy: A Comparative Review of Guidelines. [PMID 41194372](https://pubmed.ncbi.nlm.nih.gov/41194372/) → 지침 간 불일치의 범위 (모델이 왜 지침 하나에 맞추지 않는지의 근거)
8. Bicocca MJ, 2018, *Eur J Obstet Gynecol Reprod Biol* — Review of six national and regional guidelines. [PMID 30396107](https://pubmed.ncbi.nlm.nih.gov/30396107/) → 같은 취지
9. Mitchell AL, 2021, *BJOG* — Re-evaluating diagnostic thresholds for intrahepatic cholestasis of pregnancy. [PMID 33586324](https://pubmed.ncbi.nlm.nih.gov/33586324/) → 비공복 ≥19 µmol/L 문턱; 모델이 "검사값"과 "내인성 분율"을 분리해 보고하는 이유
10. Rodriguez N, 2026, *Am J Gastroenterol* — Evaluating the 2023 EASL Guidelines: Risk Stratification and Outcomes. [PMID 41351236](https://pubmed.ncbi.nlm.nih.gov/41351236/) → 층화의 외부 검증
11. Tran TT, 2016, *Am J Gastroenterol* — ACG Clinical Guideline: Liver Disease and Pregnancy. [PMID 26832651](https://pubmed.ncbi.nlm.nih.gov/26832651/) → 간질환 맥락의 위치
12. Contreras Vidal C, 2026, *Medwave* — ICP in Chile: epidemiological change and a microbiological hypothesis. [PMID 42008798](https://pubmed.ncbi.nlm.nih.gov/42008798/) → 지역·시대 간 발생률 변동 (모델은 유전 소인 벡터로 표현)
13. Piechota J, 2020, *J Clin Med* — Intrahepatic Cholestasis in Pregnancy: Review of the Literature. [PMID 32384779](https://pubmed.ncbi.nlm.nih.gov/32384779/) → 배경 서술

---

## 2. 주산기 결과와 사산 위험 — 모델이 적합한 세 개의 숫자

이 절의 1번 문헌이 모델에서 **유일하게 결과 역학에 적합된 세 상수**
(`HSB0`, `HSBSC`, `HN`) 의 출처입니다.

14. **Ovadia C, 2019, *Lancet*** — Association of adverse perinatal outcomes of ICP with biochemical markers: aggregate and individual patient data meta-analyses. [PMID 30773280](https://pubmed.ncbi.nlm.nih.gov/30773280/) → **사산 0.13% / 0.28% / 3.44% (TBA <40 / 40–99 / ≥100 µmol/L).** `icp_calibration.py` C 절에서 이 세 값에 `HSB0`·`HSBSC`·`HN` 을 적합했고, D 절에서 같은 세 값을 세 개의 서로 다른 구동변수에 다시 적합해 문턱의 급경사가 어디에 있는지 측정했습니다
15. Sarker M, 2022, *Am J Obstet Gynecol* — Beyond stillbirth: association of ICP severity and adverse outcomes. [PMID 36008054](https://pubmed.ncbi.nlm.nih.gov/36008054/) → 사산 외 엔드포인트(태변·조산·NICU)의 층화별 빈도 → `HM0`, `HMSC`, `HPT0` 의 기준
16. Zhou Q, 2024, *BMC Pregnancy Childbirth* — Severity of ICP increases risks of adverse outcomes beyond stillbirth: 15,826 patients. [PMID 38997626](https://pubmed.ncbi.nlm.nih.gov/38997626/) → 같은 목적, 대규모 검증
17. Huang X, 2024, *PLoS One* — Systematic review and meta-analysis of obstetric and neonatal outcomes. [PMID 38833446](https://pubmed.ncbi.nlm.nih.gov/38833446/) → NICU 입원률의 기준 범위
18. Yao L, 2025, *BMC Pregnancy Childbirth* — Total bile acid concentrations and adverse perinatal outcomes, including asymptomatic hypercholanemia. [PMID 41102717](https://pubmed.ncbi.nlm.nih.gov/41102717/) → 무증상 고담즙산혈증의 존재; 모델에서 가려움 축과 담즙산 축이 분리되어 있어야 하는 임상적 근거
19. Çalık MG, 2026, *BMC Pregnancy Childbirth* — Perinatal outcomes in relation to transaminase and bile acid levels. [PMID 41987074](https://pubmed.ncbi.nlm.nih.gov/41987074/) → ALT 와 TBA 가 서로 다른 정보를 담는다는 관찰
20. Granese R, 2023, *J Clin Med* — Maternal and Neonatal Outcomes in ICP. [PMID 37445442](https://pubmed.ncbi.nlm.nih.gov/37445442/) → 태변 착색 빈도
21. Feng F, 2024, *Sci Rep* — Clinical subtypes and bile acid levels with pregnancy outcomes. [PMID 38806569](https://pubmed.ncbi.nlm.nih.gov/38806569/) → 아형별 결과
22. Sarker MR, 2025, *Am J Perinatol* — Adverse Outcomes Associated with Progressive ICP. [PMID 39592109](https://pubmed.ncbi.nlm.nih.gov/39592109/) → 시간에 따라 악화하는 경우의 결과; 모델 G 절이 지적한 "매끄러운 궤적" 한계와 직접 관련
23. Sarker MR, 2025, *Pregnancy (Hoboken)* — Earlier diagnosis of ICP and adverse pregnancy outcomes. [PMID 40959761](https://pubmed.ncbi.nlm.nih.gov/40959761/) → 조기 발현의 예후
24. Hocaoglu M, 2025, *Fetal Pediatr Pathol* — Early- and Late-Onset ICP: maternal and neonatal outcomes. [PMID 39932778](https://pubmed.ncbi.nlm.nih.gov/39932778/) → 발현 시기와 결과
25. Niculae LE, 2025, *Biomedicines* — Neonatal Impact Through the Lens of Current Evidence. [PMID 41007629](https://pubmed.ncbi.nlm.nih.gov/41007629/) → 신생아 결과
26. Nielsen JH, 2021, *Acta Obstet Gynecol Scand* — Differentiated timing of induction for women with ICP. [PMID 32970824](https://pubmed.ncbi.nlm.nih.gov/32970824/) → 층화별 유도분만 시기의 관찰 데이터; 모델 G 절의 대조
27. Zehner L, 2023, *Arch Gynecol Obstet* — Evaluation of obstetric management in German maternity units. [PMID 36030428](https://pubmed.ncbi.nlm.nih.gov/36030428/) → 실제 진료의 분만 시기 분포
28. Xu T, 2024, *Int J Gynaecol Obstet* — Risk-stratified management strategies: tertiary centre review. [PMID 37470272](https://pubmed.ncbi.nlm.nih.gov/37470272/) → 층화 관리의 실제
29. Capatina N, 2024, *Obstet Med* — Meta-analyses in cholestatic pregnancy: the outstanding clinical questions. [PMID 39262915](https://pubmed.ncbi.nlm.nih.gov/39262915/) → 미해결 질문 목록; 이 모델이 겨냥한 두 가지 역설의 서술

---

## 3. 유전적 소인 (Genetic susceptibility)

모델의 `GBSEP`, `GMDR3`, `GFIC1`, `GSULT`, `GFXR` 다섯 개 환자 기술자의 근거.

30. Dixon PH, 2022, *Nat Commun* — GWAS meta-analysis of ICP implicates multiple hepatic genes and regulatory elements. [PMID 35977952](https://pubmed.ncbi.nlm.nih.gov/35977952/) → ICP 가 다유전자성이라는 것; 모델이 단일 변이가 아니라 연속적 수송능 배수를 쓰는 이유
31. Tyrmi JS, 2026, *Nat Commun* — Genome-wide meta-analysis identifies genetic drivers of bile acid metabolism in ICP. [PMID 42178310](https://pubmed.ncbi.nlm.nih.gov/42178310/) → 같은 취지, 담즙산 대사 축 중심
32. Dixon PH, 2009, *Gut* — Contribution of variant alleles of ABCB11 to susceptibility to ICP. [PMID 18987030](https://pubmed.ncbi.nlm.nih.gov/18987030/) → **`GBSEP` (ABCB11/BSEP 수송능 저하)** 의 근거
33. Zhang D, 2025, *Sci Rep* — Association between ABCB4 variants and intrahepatic cholestasis of pregnancy. [PMID 39865141](https://pubmed.ncbi.nlm.nih.gov/39865141/) → **`GMDR3`** 의 근거
34. Prescher M, 2019, *Biol Chem* — ABCB4/MDR3 in health and disease. [PMID 30730833](https://pubmed.ncbi.nlm.nih.gov/30730833/) → MDR3 가 인지질 floppase 라는 기능적 근거 → 모델의 담즙 담즙산:PC 비 (`BSPC`) 구조
35. Sticova E, 2020, *Ann Hepatol* — ABCB4 disease: many faces of one gene deficiency. [PMID 31759867](https://pubmed.ncbi.nlm.nih.gov/31759867/) → ABCB4 표현형 스펙트럼 (GGT 상승 여부 포함)
36. Nayagam JS, 2020, *Aliment Pharmacol Ther* — Liver disease in adults with variants in ABCB11, ABCB4 and ATP8B1. [PMID 33070363](https://pubmed.ncbi.nlm.nih.gov/33070363/) → **`GFIC1` (ATP8B1)** 포함 세 유전자의 성인 표현형
37. Yeap SP, 2019, *J Gastroenterol Hepatol* — Biliary transporter gene mutations in severe ICP. [PMID 29992621](https://pubmed.ncbi.nlm.nih.gov/29992621/) → 중증 ICP 에서 변이 빈도가 높다는 관찰; 모델의 중증도-유전형 대응
38. Aydın GA, 2020, *Taiwan J Obstet Gynecol* — The role of genetic mutations in ICP. [PMID 32917322](https://pubmed.ncbi.nlm.nih.gov/32917322/) → 같은 취지
39. Stieger B, 2011, *Expert Opin Drug Metab Toxicol* — Genetic variations of bile salt transporters as predisposing factors. [PMID 21320040](https://pubmed.ncbi.nlm.nih.gov/21320040/) → 약물유발 담즙정체와 ICP 의 공통 취약성
40. Pauli-Magnus C, 2010, *Semin Liver Dis* — Genetic determinants of drug-induced cholestasis and ICP. [PMID 20422497](https://pubmed.ncbi.nlm.nih.gov/20422497/) → 같은 취지
41. Degiorgio D, 2016, *J Gastroenterol* — ABCB4 mutations in adult patients with cholestatic liver disease. [PMID 26324191](https://pubmed.ncbi.nlm.nih.gov/26324191/) → 변이별 잔존 기능; `GMDR3` 값 범위의 근거
42. Falcão D, 2022, *Dig Liver Dis* — The wide phenotypic and genetic spectrum of ABCB4 deficiency. [PMID 34376370](https://pubmed.ncbi.nlm.nih.gov/34376370/) → 같은 취지
43. Henkel SA, 2019, *World J Hepatol* — Expanding etiology of progressive familial intrahepatic cholestasis. [PMID 31183005](https://pubmed.ncbi.nlm.nih.gov/31183005/) → PFIC 스펙트럼 (TJP2, NR1H4 포함)

---

## 4. 성호르몬이 BSEP 를 억제한다 — 병변 그 자체

모델에서 질환을 일으키는 두 개의 억제항 (`KI_E2G`, `KI_P4S`) 과, 제3삼분기
발현 및 분만 후 소실이 **파라미터 변경 없이** 나오는 이유.

44. **Vallejo M, 2006, *J Hepatol*** — Potential role of trans-inhibition of the bile salt export pump by progesterone metabolites in the etiopathogenesis of ICP. [PMID 16458994](https://pubmed.ncbi.nlm.nih.gov/16458994/) → **`KI_P4S` (프로게스테론 설페이트의 BSEP trans-억제)** 의 직접 근거
45. **Abu-Hayyeh S, 2013, *Hepatology*** — ICP levels of sulfated progesterone metabolites inhibit FXR resulting in a cholestatic phenotype. [PMID 22961653](https://pubmed.ncbi.nlm.nih.gov/22961653/) → P4-S 의 FXR 길항 (모델 지도의 `P4S → FXR` 억제 엣지)
46. Abu-Hayyeh S, 2010, *J Biol Chem* — Inhibition of NTCP-mediated bile acid transport by cholestatic sulfated progesterone metabolites. [PMID 20177056](https://pubmed.ncbi.nlm.nih.gov/20177056/) → NTCP 억제; 모델의 `updep` 항이 표현하는 흡수 저해와 정합
47. Abu-Hayyeh S, 2016, *Hepatology* — Prognostic and mechanistic potential of progesterone sulfates in ICP and pruritus gravidarum. [PMID 26426865](https://pubmed.ncbi.nlm.nih.gov/26426865/) → **P4-S 가 담즙산보다 가려움과 잘 맞는다**는 관찰 → 모델이 오토탁신을 P4S 로 구동하는 구조적 결정의 핵심 근거 (`EA_P4S`, `KA_P4S`, `HP4S`)
48. Meng LJ, 1997, *Hepatology* — Progesterone metabolites and bile acids in serum of patients with ICP: effect of UDCA therapy. [PMID 9398000](https://pubmed.ncbi.nlm.nih.gov/9398000/) → **UDCA 가 프로게스테론 설페이트를 정상화하지 못한다**는 관찰; 모델에서 UDCA 가 가려움을 거의 못 줄이는 이유
49. Meng LJ, 1997, *J Hepatol* — Effects of UDCA on conjugated bile acids and progesterone metabolites in serum and urine of ICP. [PMID 9453429](https://pubmed.ncbi.nlm.nih.gov/9453429/) → 같은 취지
50. Meng LJ, 1997, *J Hepatol* — Profiles of bile acids and progesterone metabolites in urine and serum of women with ICP. [PMID 9288610](https://pubmed.ncbi.nlm.nih.gov/9288610/) → P4-S 농도 범위 → `KFP4S`, `KOUTP4S` 규모
51. Reyes H, 2000, *Ann Med* — Bile acids and progesterone metabolites in ICP. [PMID 10766400](https://pubmed.ncbi.nlm.nih.gov/10766400/) → 두 축의 병행 상승
52. Mitchell AL, 2025, *Am J Physiol Gastrointest Liver Physiol* — Progesterone sulfates are enterohepatically recycled and stimulate GPBAR1-mediated gut hormone release. [PMID 39888313](https://pubmed.ncbi.nlm.nih.gov/39888313/) → P4-S 가 장간순환한다는 점 (모델에서는 단순화하여 혈중 상태변수로만 둠 — **명시적 단순화**)
53. Sanchon-Sanchez P, 2024, *Biochim Biophys Acta Mol Basis Dis* — Relationship between cholestasis and altered progesterone metabolism in the placenta-maternal liver tandem. [PMID 37956602](https://pubmed.ncbi.nlm.nih.gov/37956602/) → 태반-모체 간 축의 프로게스테론 대사
54. Ren L, 2019, *Biochem Biophys Res Commun* — Epiallopregnanolone sulfate induces hepatic bile acid accumulation and liver injury. [PMID 31575408](https://pubmed.ncbi.nlm.nih.gov/31575408/) → P4-S 가 인과적이라는 실험적 근거
55. Liu H, 2023, *Nat Commun* — Structural basis of bile salt extrusion and small-molecule inhibition in human BSEP. [PMID 37949847](https://pubmed.ncbi.nlm.nih.gov/37949847/) → BSEP 억제의 구조적 기반; 경쟁적 억제 형태의 타당성
56. Hosamani S, 2024, *J Phys Chem Lett* — Cholesterol allosterically modulates the structure and dynamics of ABCB11. [PMID 39058973](https://pubmed.ncbi.nlm.nih.gov/39058973/) → 막 환경이 BSEP 활성을 바꾼다는 근거 → `GFIC1` (ATP8B1) 이 BSEP 를 곱하는 구조

---

## 5. 담즙산 합성 · 되먹임 · 수송 (Synthesis, feedback, transport)

57. Chiang JY, 2009, *J Lipid Res* — Bile acids: regulation of synthesis. [PMID 19346330](https://pubmed.ncbi.nlm.nih.gov/19346330/) → `VSYN`, `FRCA` (CYP7A1/CYP8B1 분기), FXR-SHP 되먹임 구조 전체
58. Holt JA, 2003, *Genes Dev* — Definition of a novel growth factor-dependent signal cascade for the suppression of bile acid biosynthesis. [PMID 12815072](https://pubmed.ncbi.nlm.nih.gov/12815072/) → 회장 FXR–FGF19 축 → `KS_F19`, `EC50FXI`, `KF19`
59. Song KH, 2009, *Hepatology* — Bile acids activate FGF19 signaling in human hepatocytes to inhibit CYP7A1. [PMID 19085950](https://pubmed.ncbi.nlm.nih.gov/19085950/) → 사람 세포에서의 FGF19–CYP7A1 억제
60. Zollner G, 2007, *Liver Int* — Expression of bile acid synthesis and detoxification enzymes and the alternative efflux pump MRP4 in PBC. [PMID 17696930](https://pubmed.ncbi.nlm.nih.gov/17696930/) → **담즙정체에서 MRP4 가 유도된다**는 사람 데이터 → `EM4`, `KM4`, `HM4` (기저에서는 잠잠하고 담즙정체에서 크게 유도되는 Hill 형태의 근거)
61. Schaffner CA, 2015, *Liver Int* — The organic solute transporters alpha and beta are induced by hypoxia in human hepatocytes. [PMID 24703425](https://pubmed.ncbi.nlm.nih.gov/24703425/) → OSTα-β 유도성
62. Xu S, 2014, *Am J Physiol Gastrointest Liver Physiol* — A novel RARα/CAR-mediated mechanism for regulation of human OSTβ. [PMID 24264050](https://pubmed.ncbi.nlm.nih.gov/24264050/) → 기저측 배출로의 전사 조절
63. van Dijk R, 2015, *Liver Int* — Characterization and treatment of persistent hepatocellular secretory failure. [PMID 24905729](https://pubmed.ncbi.nlm.nih.gov/24905729/) → 간세포 분비 실패의 임상 표현형
64. Kastrinou Lampou V, 2023, *Toxicol In Vitro* — Novel insights into bile acid detoxification via CYP, UGT and SULT enzymes. [PMID 36473578](https://pubmed.ncbi.nlm.nih.gov/36473578/) → **`V6A` (CYP3A4 6α-수산화)** 와 **`VLCAS` (SULT2A1 LCA 설페이트화)** 두 해독 경로의 근거
65. Bansal S, 2016, *J Chromatogr B* — Fast and sensitive quantification of human liver cytosolic lithocholic acid sulfation. [PMID 26773894](https://pubmed.ncbi.nlm.nih.gov/26773894/) → 사람 간 LCA 설페이트화 용량 → `VLCAS`, `KMLCAS` 규모
66. Mao F, 2019, *J Biol Chem* — Increased sulfation of bile acids in NTCP deficiency. [PMID 31201272](https://pubmed.ncbi.nlm.nih.gov/31201272/) → 담즙산 배출이 막히면 설페이트화·신배설이 대체 경로가 된다 → `ESULT`
67. Chatterjee B, 2005, *Methods Enzymol* — Vitamin D receptor regulation of SULT2A1. [PMID 16399349](https://pubmed.ncbi.nlm.nih.gov/16399349/) → LCA–VDR–SULT2A1 되먹임 (지도에 반영, ODE 에서는 PXR 항으로 통합 — **명시적 단순화**)
68. Mok HY, 1977, *Gastroenterology* — Regulation of pool size of bile acids in man. [PMID 892372](https://pubmed.ncbi.nlm.nih.gov/892372/) → **사람 담즙산 풀 크기와 하루 순환 횟수** → `KBILE`, `KTR`, `VASBT` 의 규모 결정
69. Low-Beer TS, 1973, *Br Med J* — Regulation of bile salt pool size in man. [PMID 4704519](https://pubmed.ncbi.nlm.nih.gov/4704519/) → 같은 목적
70. Hepner GW, 1975, *Gastroenterology* — Effect of decreased gallbladder stimulation on enterohepatic cycling and kinetics of bile acids. [PMID 1132637](https://pubmed.ncbi.nlm.nih.gov/1132637/) → 담낭 저류와 순환 빈도
71. Angelin B, 1978, *J Lipid Res* — Individual serum bile acid concentrations in normo- and hyperlipoproteinemia: relation to pool size. [PMID 670831](https://pubmed.ncbi.nlm.nih.gov/670831/) → **정상 혈중 담즙산 농도 (모델 정상 임신 4–5 µmol/L 의 기준)**
72. Hofmann AF, 1987, *Gastroenterology* — Simulation of the metabolism and enterohepatic circulation of endogenous deoxycholic acid in humans using a physiologic pharmacokinetic model. [PMID 3623017](https://pubmed.ncbi.nlm.nih.gov/3623017/) → 이 모델의 장간순환 구획 구조의 직계 선행 연구; 대장 7α-탈수산화 (`KDH_CA`) 처리 방식
73. Cravetto C, 1988, *Hepatology* — Computer simulation of portal venous shunting and other isolated hepatobiliary defects of the enterohepatic circulation using a physiologic model. [PMID 3391514](https://pubmed.ncbi.nlm.nih.gov/3391514/) → 간 초회통과 추출률 결함의 모델링 선례 → `EUP*`, `JHALF`, `KHLSAT`
74. Hulzebos CV, 2001, *J Lipid Res* — Measurement of parameters of cholic acid kinetics in plasma using stable isotope dilution. [PMID 11714862](https://pubmed.ncbi.nlm.nih.gov/11714862/) → 콜산 동태 파라미터
75. Rutgeerts P, 1983, *J Lipid Res* — The enterohepatic circulation of bile acids during continuous liquid formula perfusion of the duodenum. [PMID 6875385](https://pubmed.ncbi.nlm.nih.gov/6875385/) → 장내 담즙산 농도 (mmol/L 규모) → `MIC50` 및 미셀 형성 항

---

## 6. 태반 담즙산 수송 (Placental bile acid transport)

모델의 `VP`, `KMP*`, `PS*`, `TRANSIN`. **이 절의 문헌은 태반 수송이
포화성이며 양방향임을 지지하지만, 사람에서 `VP` 의 절대값을 직접 측정한
연구는 찾지 못했습니다 — `VP` 는 정상 임신의 제대혈 농도에 맞춘 모델
선택입니다.**

76. Monte MJ, 1995, *Pediatr Res* — Relationship between bile acid transplacental gradients and transport across the fetal-facing plasma membrane of the human trophoblast. [PMID 7478809](https://pubmed.ncbi.nlm.nih.gov/7478809/) → **태아-모체 농도 기울기와 영양막 수송의 관계** → `PS*` 와 `JF2M` 의 구조
77. Marin JJ, 1995, *Am J Physiol* — ATP-dependent bile acid transport across microvillous membrane of human term trophoblast. [PMID 7733292](https://pubmed.ncbi.nlm.nih.gov/7733292/) → 능동(ATP 의존) 수송의 존재 → 포화 가능한 `VP`
78. Marín JJ, 2005, *Ann Hepatol* — Molecular bases of the excretion of fetal bile acids and pigments through the fetal liver-placenta-maternal liver pathway. [PMID 16010240](https://pubmed.ncbi.nlm.nih.gov/16010240/) → **태아 담즙산의 유일한 실질적 배출 경로가 태반이라는 구조** → 모델에서 태아 구획의 escape route 설계
79. Macias RI, 2000, *Hepatology* — Effect of maternal cholestasis on bile acid transfer across the rat placenta-maternal liver tandem. [PMID 10733555](https://pubmed.ncbi.nlm.nih.gov/10733555/) → 모체 담즙정체가 태반 이동을 손상시킨다 → `CAPUTIL` 포화 개념
80. Macias RI, 2009, *World J Gastroenterol* — Excretion of biliary compounds during intrauterine life. [PMID 19230042](https://pubmed.ncbi.nlm.nih.gov/19230042/) → 태아 배출 경로의 종합
81. Ontsouka E, 2021, *Int J Mol Sci* — Placental Expression of Bile Acid Transporters in Intrahepatic Cholestasis of Pregnancy. [PMID 34638773](https://pubmed.ncbi.nlm.nih.gov/34638773/) → ICP 태반의 수송체 발현 변화 (사람) → `GP` 를 환자 기술자로 둔 근거
82. St-Pierre MV, 2000, *Am J Physiol Regul Integr Comp Physiol* — Expression of members of the multidrug resistance protein family in human term placenta. [PMID 11004020](https://pubmed.ncbi.nlm.nih.gov/11004020/) → 태반 MRP 발현
83. Blazquez AG, 2014, *Toxicol Appl Pharmacol* — Acetaminophen impairs the placental barrier to bile acids via BCRP. [PMID 24631341](https://pubmed.ncbi.nlm.nih.gov/24631341/) → 태반 장벽이 약물로 조절 가능하다는 근거 (모델의 `GP` 조작 시나리오)
84. Serrano MA, 1998, *J Hepatol* — Beneficial effect of UDCA on alterations induced by cholestasis of pregnancy in bile acid transport across the human placenta. [PMID 9625319](https://pubmed.ncbi.nlm.nih.gov/9625319/) → **UDCA 가 태반 수송에 개입한다는 직접 근거** → `KMP5`, `TRANSIN` (UDCA 의 수송체 경합) 의 출처
85. Serrano MA, 2003, *J Pharmacol Exp Ther* — Effect of UDCA on the impairment induced by maternal cholestasis in the rat placenta-maternal liver tandem. [PMID 12606635](https://pubmed.ncbi.nlm.nih.gov/12606635/) → 같은 취지
86. Estiú MC, 2015, *Br J Clin Pharmacol* — Effect of UDCA on the altered progesterone and bile acid homeostasis in the mother-placenta-foetus trio. [PMID 25099365](https://pubmed.ncbi.nlm.nih.gov/25099365/) → 모체-태반-태아 삼자에서 UDCA 의 종별 효과
87. Basu S, 2025, *Am J Physiol Gastrointest Liver Physiol* — Unresolved alterations in bile acid composition and dyslipidemia in maternal and cord blood after UDCA treatment for ICP. [PMID 39947696](https://pubmed.ncbi.nlm.nih.gov/39947696/) → **UDCA 치료 후에도 제대혈 담즙산 조성 이상이 남는다** → 모델이 "검사값은 내려가도 태아 소수성 부하는 그만큼 내려가지 않는다"고 말하는 부분의 직접적 실험적 지지
88. Wang G, 2026, *Mol Cell Biochem* — Bile acid-induced TBG depletion promotes trophoblast apoptosis in ICP. [PMID 42371391](https://pubmed.ncbi.nlm.nih.gov/42371391/) → 영양막 세포자멸사 (`TROPH`)
89. Wu WB, 2015, *Placenta* — FXR agonist protects against bile acid induced damage and oxidative stress in mouse placenta. [PMID 25747729](https://pubmed.ncbi.nlm.nih.gov/25747729/) → 태반 산화 스트레스
90. Keitel V, 2013, *Placenta* — Effect of maternal cholestasis on TGR5 expression in human and rat placenta at term. [PMID 23849932](https://pubmed.ncbi.nlm.nih.gov/23849932/) → 태반 TGR5

---

## 7. 태아 심근 — 문턱이 실제로 있는 곳

**모델의 가장 중요한 절.** `IMAXGJ`, `IC50GJ`, `HGJ`, `ECA`, `KCA` 는 모두
이 절에서 왔고, `icp_calibration.py` I 절의 절제 실험은 **여기의 Hill 형태를
선형으로 바꾸면 100 µmol/L 문턱이 사라진다**는 것을 보입니다.

91. **Williamson C, 2001, *Clin Sci (Lond)*** — The bile acid taurocholate impairs rat cardiomyocyte function: a proposed mechanism for intra-uterine fetal death in obstetric cholestasis. [PMID 11256973](https://pubmed.ncbi.nlm.nih.gov/11256973/) → **모델의 태아 심근 모듈 전체의 출발점**
92. **Gorelik J, 2002, *Clin Sci (Lond)*** — Taurocholate induces changes in rat cardiomyocyte contraction and calcium dynamics. [PMID 12149111](https://pubmed.ncbi.nlm.nih.gov/12149111/) → **`ECA`, `KCA` (칼슘 과부하)** 의 직접 근거
93. **Gorelik J, 2004, *BJOG*** — Comparison of the arrhythmogenic effects of tauro- and glycoconjugates of cholic acid in an in vitro study of rat cardiomyocytes. [PMID 15270939](https://pubmed.ncbi.nlm.nih.gov/15270939/) → **담즙산 종별로 부정맥 유발능이 다르다** → 모델이 종별 세포독성 가중치 `WTOX` 를 쓰는 근거
94. Sheikh Abdul Kadir SH, 2010, *PLoS One* — Bile acid-induced arrhythmia is mediated by muscarinic M2 receptors in neonatal rat cardiomyocytes. [PMID 20300620](https://pubmed.ncbi.nlm.nih.gov/20300620/) → M2 수용체 경로 (지도의 `MRECEPT`)
95. Ibrahim E, 2018, *Sci Rep* — Bile acids and their respective conjugates elicit different responses in neonatal cardiomyocytes: role of Gi protein, muscarinic receptors and TGR5. [PMID 29740092](https://pubmed.ncbi.nlm.nih.gov/29740092/) → 종별·접합체별 반응 차이의 정량적 근거 → `WTOX` 순서
96. Miragoli M, 2011, *Hepatology* — A protective antiarrhythmic role of ursodeoxycholic acid in an in vitro rat model of the cholestatic fetal heart. [PMID 21809354](https://pubmed.ncbi.nlm.nih.gov/21809354/) → **UDCA 의 직접적 항부정맥 효과**; 모델은 이를 `WTOX[UDCA] = 0.02` 로만 표현하고 별도의 심근 보호항을 두지 않았습니다 — **모델이 UDCA 의 태아 심장 보호를 과소평가할 수 있는 지점으로 명시**
97. Adeyemi O, 2017, *PLoS One* — UDCA prevents ventricular conduction slowing and arrhythmia by restoring T-type calcium current in fetuses during cholestasis. [PMID 28934223](https://pubmed.ncbi.nlm.nih.gov/28934223/) → 같은 취지, in vivo 태아
98. Gorelik J, 2003, *BJOG* — Dexamethasone and UDCA protect against the arrhythmogenic effect of taurocholate. [PMID 12742331](https://pubmed.ncbi.nlm.nih.gov/12742331/) → 덱사메타손의 in vitro 보호 (임상 RCT 는 음성 — 지도에 그렇게 표기)
99. Gorelik J, 2006, *BJOG* — Genes encoding bile acid, phospholipid and anion transporters are expressed in a human fetal cardiomyocyte culture. [PMID 16637898](https://pubmed.ncbi.nlm.nih.gov/16637898/) → 태아 심근세포가 담즙산 수송체를 발현한다 (사람)
100. Abdul Kadir SH, 2009, *J Cell Mol Med* — Embryonic stem cell-derived cardiomyocytes as a model to study fetal arrhythmia related to maternal disease. [PMID 19438812](https://pubmed.ncbi.nlm.nih.gov/19438812/) → 사람 유래 모델
101. Schultz F, 2016, *Prog Biophys Mol Biol* — The protective effect of UDCA in an in vitro model of the human fetal heart occurs via targeting cardiac fibroblasts. [PMID 26777584](https://pubmed.ncbi.nlm.nih.gov/26777584/) → 섬유아세포 매개 기전 (모델에 미포함 — **명시적 누락**)
102. Williamson C, 2011, *Dig Dis* — Bile acid signaling in fetal tissues: implications for ICP. [PMID 21691106](https://pubmed.ncbi.nlm.nih.gov/21691106/) → 태아 조직 담즙산 신호의 종합
103. Guerra M, 2025, *Am J Perinatol* — Fetal TEI Index in Pregnancies with ICP: A Case-Control Study. [PMID 40112872](https://pubmed.ncbi.nlm.nih.gov/40112872/) → **사람 태아에서 심기능 지표가 실제로 변한다는 임상 관찰** → 모델의 심근 모듈이 사람에서도 검증 가능하다는 근거
104. Şık ME, 2026, *J Perinat Med* — Evaluation of Doppler parameters and obstetric outcomes in ICP. [PMID 41770037](https://pubmed.ncbi.nlm.nih.gov/41770037/) → 도플러가 급성 사건을 예측하지 못한다는 관찰 (모델의 임상적 함의와 정합)
105. Sepúlveda WH, 1991, *Eur J Obstet Gynecol Reprod Biol* — Vasoconstrictive effect of bile acids on isolated human placental chorionic veins. [PMID 1773876](https://pubmed.ncbi.nlm.nih.gov/1773876/) → **`EV`, `KV` (융모판 정맥 수축)** 의 직접 근거

---

## 8. 가려움 축 — 담즙산 축과 분리되어 있다

모델의 두 번째 축 (`ATX`, `LPA`, `ITCHC`). **이 절의 핵심 관찰은 "혈중
오토탁신은 가려움과 상관하지만 총담즙산은 그렇지 않다"이며, 그 관찰이
모델에서 `EA_CH` 를 작게 두고 `KA_CH` 를 낮게 (즉 포화되게) 둔 이유입니다.**

106. **Kremer AE, 2010, *Gastroenterology*** — Lysophosphatidic acid is a potential mediator of cholestatic pruritus. [PMID 20546739](https://pubmed.ncbi.nlm.nih.gov/20546739/) → LPA 가 매개자라는 최초의 제시 → 모델의 `LPA` 상태변수
107. **Kremer AE, 2012, *Hepatology*** — Serum autotaxin is increased in pruritus of cholestasis, but not of other origin, and responds to therapeutic interventions. [PMID 22473838](https://pubmed.ncbi.nlm.nih.gov/22473838/) → **`ERIFA` (리팜피신이 오토탁신을 낮춘다) 와, UDCA 가 그렇지 않다는 관찰의 출처.** 모델의 두 축 분리 구조 전체가 이 논문에 걸려 있습니다
108. Beuers U, 2023, *Nat Rev Gastroenterol Hepatol* — Mechanisms of pruritus in cholestasis: understanding and treating the itch. [PMID 36307649](https://pubmed.ncbi.nlm.nih.gov/36307649/) → 최신 종합; 치료제별 작용 축
109. Beuers U, 2014, *Hepatology* — Pruritus in cholestasis: facts and fiction. [PMID 24807046](https://pubmed.ncbi.nlm.nih.gov/24807046/) → 담즙산-가려움 상관의 약함을 명시
110. Oude Elferink RP, 2011, *Dig Dis* — The molecular mechanism of cholestatic pruritus. [PMID 21691108](https://pubmed.ncbi.nlm.nih.gov/21691108/) → 기전 종합
111. Oude Elferink RP, 2015, *Biochim Biophys Acta* — Lysophosphatidic acid and signaling in sensory neurons. [PMID 25218302](https://pubmed.ncbi.nlm.nih.gov/25218302/) → LPA→감각신경 → `HLP` (다단계 증폭의 겉보기 협동성) 의 근거
112. Wunsch E, 2016, *Sci Rep* — Serum autotaxin is a marker of the severity of liver injury and overall survival in cholestatic liver diseases. [PMID 27506882](https://pubmed.ncbi.nlm.nih.gov/27506882/) → 오토탁신 농도 범위 (모델의 ATX 배수 1.5–4)
113. Keune WJ, 2016, *Nat Commun* — Steroid binding to autotaxin links bile salts and lysophosphatidic acid signalling. [PMID 27075612](https://pubmed.ncbi.nlm.nih.gov/27075612/) → **담즙산이 오토탁신에 직접 결합한다** → 모델의 약한 `EA_CH` 항의 구조적 근거
114. Langedijk JAGM, 2021, *Biochim Biophys Acta Mol Basis Dis* — Inhibition of autotaxin by bile salts and bile salt-like molecules increases its expression by feedback regulation. [PMID 34389475](https://pubmed.ncbi.nlm.nih.gov/34389475/) → **담즙산이 오토탁신을 억제하면 발현이 되먹임으로 올라간다** → 담즙산-가려움 관계가 단조롭지 않은 이유; 모델의 포화 항과 정합
115. **Yu H, 2019, *Elife*** — MRGPRX4 is a bile acid receptor for human cholestatic itch. [PMID 31500698](https://pubmed.ncbi.nlm.nih.gov/31500698/) → **`WBA`, `KBA` (직접 담즙산-가려움 경로)** 의 근거
116. Meixiong J, 2019, *Proc Natl Acad Sci U S A* — MRGPRX4 is a GPCR activated by bile acids that may contribute to cholestatic pruritus. [PMID 31068464](https://pubmed.ncbi.nlm.nih.gov/31068464/) → 같은 발견의 독립 보고
117. Yu H, 2021, *Semin Liver Dis* — MRGPRX4 in Cholestatic Pruritus. [PMID 34161994](https://pubmed.ncbi.nlm.nih.gov/34161994/) → 종합
118. Rodrigo M, 2026, *Liver Int* — Defining the contribution of genetic variants in MRGPRX4 with pruritus in paediatric cholestasis. [PMID 41817014](https://pubmed.ncbi.nlm.nih.gov/41817014/) → MRGPRX4 변이가 가려움 정도를 바꾼다 (모델에 미포함 — **확장 후보로 명시**)
119. Yang J, 2026, *Nat Chem Biol* — Development of a clinically viable MRGPRX4 inverse agonist for cholestatic itch treatment. [PMID 41957282](https://pubmed.ncbi.nlm.nih.gov/41957282/) → 가려움 축을 직접 겨냥한 신약; 모델의 `WBA` 항이 표적으로 삼는 경로
120. Yang J, 2024, *Cell* — Structure-guided discovery of bile acid derivatives for treating liver diseases without causing itch. [PMID 39476841](https://pubmed.ncbi.nlm.nih.gov/39476841/) → 담즙산 축과 가려움 축이 약리적으로 분리 가능하다는 직접적 증거
121. Wolters F, 2025, *Biochim Biophys Acta Mol Basis Dis* — The role of bile salts in itch receptor activation. [PMID 40628364](https://pubmed.ncbi.nlm.nih.gov/40628364/) → 종별 수용체 활성화
122. Alemi F, 2013, *J Clin Invest* — The TGR5 receptor mediates bile acid-induced itch and analgesia. [PMID 23524965](https://pubmed.ncbi.nlm.nih.gov/23524965/) → TGR5 경로 (지도에 반영)
123. Lieu T, 2014, *Gastroenterology* — The bile acid receptor TGR5 activates the TRPA1 channel to induce itch in mice. [PMID 25194674](https://pubmed.ncbi.nlm.nih.gov/25194674/) → TRPA1 하류
124. Song MH, 2022, *Biomol Ther (Seoul)* — Lithocholic acid activates MAS-related GPCRs, contributing to itch in mice. [PMID 34263729](https://pubmed.ncbi.nlm.nih.gov/34263729/) → **LCA 가 가장 강한 가려움 유발 종**이라는 근거; `WTOX` 와 `WBA` 의 정합성
125. Akiyama T, 2013, *Neuroscience* — Neural processing of itch. [PMID 23891755](https://pubmed.ncbi.nlm.nih.gov/23891755/) → 척수 GRPR 및 중추 감작 → `KSENS`, `KSC` (긁기-가려움 되먹임)
126. Vander Does A, 2022, *Am J Clin Dermatol* — Cholestatic Itch: Our Current Understanding of Pathophysiology and Treatments. [PMID 35900649](https://pubmed.ncbi.nlm.nih.gov/35900649/) → 치료제별 축 배정의 근거
127. Kanda T, 2025, *Int J Mol Sci* — Pruritus in Chronic Cholestatic Liver Diseases. [PMID 40076514](https://pubmed.ncbi.nlm.nih.gov/40076514/) → 최신 종합
128. Trivella J, 2021, *Expert Opin Drug Saf* — Safety considerations for the management of cholestatic itch. [PMID 33836644](https://pubmed.ncbi.nlm.nih.gov/33836644/) → 안전성 (임신에서의 제약)
129. Düll MM, 2019, *Curr Gastroenterol Rep* — Treatment of Pruritus Secondary to Liver Disease. [PMID 31367993](https://pubmed.ncbi.nlm.nih.gov/31367993/) → 단계적 치료 알고리듬
130. Kremer AE, 2017, *Hautarzt* — Intrahepatic cholestasis of pregnancy: rare but important. [PMID 28074213](https://pubmed.ncbi.nlm.nih.gov/28074213/) → ICP 가려움의 임상적 특징 (야간·손발바닥)
131. Ovadia C, 2018, *J Clin Apher* — Therapeutic plasma exchange as a novel treatment for severe ICP. [PMID 30321466](https://pubmed.ncbi.nlm.nih.gov/30321466/) → 혈장교환이 효과가 있다는 관찰; 지도의 `NBD` 노드 (오토탁신 구동 자체를 제거)

---

## 9. UDCA — 시험 결과와 기전

132. **Chappell LC, 2019, *Lancet*** — Ursodeoxycholic acid versus placebo in women with ICP (PITCHES): a randomised controlled trial. [PMID 31378395](https://pubmed.ncbi.nlm.nih.gov/31378395/) → **모델 E 절이 재현하는 시험.** 복합 엔드포인트 24.8% 대 27.8%, 가려움 −0.7 cm/10 cm
133. Chappell LC, 2018, *Trials* — PITCHES protocol. [PMID 30482254](https://pubmed.ncbi.nlm.nih.gov/30482254/) → 복합 엔드포인트의 정확한 정의 (모델의 분해가 이 정의를 따름)
134. Fleminger J, 2021, *BJOG* — UDCA in ICP: a secondary analysis of the PITCHES trial. [PMID 33063439](https://pubmed.ncbi.nlm.nih.gov/33063439/) → 하위군 분석; 중증도별 효과
135. **Ovadia C, 2021, *Lancet Gastroenterol Hepatol*** — UDCA in ICP: systematic review and individual participant data meta-analysis. [PMID 33915090](https://pubmed.ncbi.nlm.nih.gov/33915090/) → **PITCHES 와 방향이 다른 결과.** 모델 E 절은 두 연구가 서로 다른 것을 측정했을 뿐 모순이 아니라고 주장하고 그 산수를 제시합니다
136. Walker KF, 2020, *Cochrane Database Syst Rev* — Pharmacological interventions for treating ICP. [PMID 32716060](https://pubmed.ncbi.nlm.nih.gov/32716060/) → 약물별 효과크기의 체계적 종합
137. Kong X, 2016, *Medicine (Baltimore)* — Effectiveness and safety of UDCA in ICP: meta-analysis. [PMID 27749550](https://pubmed.ncbi.nlm.nih.gov/27749550/) → **UDCA 의 총담즙산 감소폭 (모델 −34~−46% 의 대조값)**
138. Grand'Maison S, 2014, *J Obstet Gynaecol Can* — UDCA treatment for ICP: meta-analysis. [PMID 25184983](https://pubmed.ncbi.nlm.nih.gov/25184983/) → 같은 목적
139. Shen Y, 2019, *Clin Drug Investig* — Is it necessary to perform pharmacological interventions for ICP? Bayesian network meta-analysis. [PMID 30357607](https://pubmed.ncbi.nlm.nih.gov/30357607/) → 약물 간 순위 (모델 F 절의 대조)
140. Zhang L, 2015, *Eur Rev Med Pharmacol Sci* — UDCA and SAMe in ICP: multi-centred randomized controlled trial. [PMID 26502869](https://pubmed.ncbi.nlm.nih.gov/26502869/) → **`ESAM`, `ESAM_R` (SAMe)** 의 효과크기 근거
141. Zhang Y, 2016, *Hepat Mon* — UDCA and SAMe for ICP: meta-analysis. [PMID 27799965](https://pubmed.ncbi.nlm.nih.gov/27799965/) → 같은 목적
142. Marschall HU, 2015, *Expert Rev Gastroenterol Hepatol* — Management of ICP. [PMID 26313609](https://pubmed.ncbi.nlm.nih.gov/26313609/) → UDCA 용량 (10–20 mg/kg/d) → 모델의 투여 시나리오
143. Setoguchi T, 1984, *J Lipid Res* — Epimerization of the four 3,7-dihydroxy bile acid epimers by human fecal microorganisms. [PMID 6520544](https://pubmed.ncbi.nlm.nih.gov/6520544/) → **`KDH_UD` (UDCA → LCA 전환)** 의 근거
144. Hirano S, 1981, *J Lipid Res* — In vitro transformation of chenodeoxycholic acid and ursodeoxycholic acid by human intestinal flora. [PMID 7288282](https://pubmed.ncbi.nlm.nih.gov/7288282/) → 같은 목적, 전환 속도의 규모
145. Edenharder R, 1981, *J Lipid Res* — Epimerization of chenodeoxycholic acid to ursodeoxycholic acid by human intestinal Clostridia. [PMID 7276738](https://pubmed.ncbi.nlm.nih.gov/7276738/) → 역방향 에피머화 (모델에 미포함 — **명시적 단순화**)
146. Lee JY, 2013, *J Lipid Res* — Contribution of the 7β-hydroxysteroid dehydrogenase from *Ruminococcus gnavus* to UDCA formation in the human colon. [PMID 23729502](https://pubmed.ncbi.nlm.nih.gov/23729502/) → 관여 균종
147. Ferrari A, 1988, *Proc Soc Exp Biol Med* — In vitro transformation of cheno- and ursodeoxycholic acids by human intestinal microflora. [PMID 3368473](https://pubmed.ncbi.nlm.nih.gov/3368473/) → 같은 목적

---

## 10. 리팜피신 · 담즙산 결합제 · IBAT 억제제 · 기타

148. Hague WM, 2021, *BMC Pregnancy Childbirth* — A multi-centre, open label, randomised trial comparing UDCA with RIFampicin in ICP (TURRIFIC). [PMID 33435904](https://pubmed.ncbi.nlm.nih.gov/33435904/) → **모델 F 절이 예측하는 비교를 실제로 수행하는 시험 프로토콜.** 모델의 예측(리팜피신이 가려움에서 UDCA 를 크게 앞서고 담즙산에서는 뒤진다)은 이 시험으로 직접 검증 가능
149. Markus C, 2021, *Clin Chem Lab Med* — The BACH project: international total Bile Acid Comparison and Harmonisation project, sub-study of TURRIFIC. [PMID 34355544](https://pubmed.ncbi.nlm.nih.gov/34355544/) → **총담즙산 검사가 기관 간 조화되어 있지 않다**는 점; 모델이 "검사값"과 "내인성 분율"을 분리해 보고하는 실무적 이유
150. Kremer AE, 2014, *Dig Dis* — Advances in pathogenesis and management of pruritus in cholestasis. [PMID 25034299](https://pubmed.ncbi.nlm.nih.gov/25034299/) → 리팜피신·콜레스티라민·날트렉손·설트랄린의 단계적 위치
151. Mittal A, 2016, *Curr Probl Dermatol* — Cholestatic Itch Management. [PMID 27578083](https://pubmed.ncbi.nlm.nih.gov/27578083/) → 같은 목적
152. Peverelle M, 2026, *Aliment Pharmacol Ther* — Review: Ileal Bile Acid Transport (IBAT) inhibitors as an emerging treatment for cholestatic liver disease. [PMID 41953994](https://pubmed.ncbi.nlm.nih.gov/41953994/) → **`KIODEV` (IBAT 억제) 시나리오**의 약리적 근거. 임신 적응증은 없으며 모델에서도 가설로만 표기
153. Lacey G, 2025, *Clin Ther* — Indirect comparison of maralixibat and odevixibat for PFIC. [PMID 40544071](https://pubmed.ncbi.nlm.nih.gov/40544071/) → 효과크기의 규모
154. Muntaha HST, 2022, *J Clin Med* — IBAT blockers for cholestatic liver disease in Alagille syndrome: systematic review and meta-analysis. [PMID 36556142](https://pubmed.ncbi.nlm.nih.gov/36556142/) → 담즙산 감소폭
155. Jeyaraj R, 2024, *Lancet Child Adolesc Health* — Paediatric research sets new standards for therapy in paediatric and adult cholestasis. [PMID 38006895](https://pubmed.ncbi.nlm.nih.gov/38006895/) → IBAT 억제제의 위치
156. ElSalem SA, 2026, *Cureus* — Refractory ICP in twin gestation managed with UDCA and adjunctive cholestyramine. [PMID 41694888](https://pubmed.ncbi.nlm.nih.gov/41694888/) → 콜레스티라민 병용의 임상 사례
157. Pongcharoen P, 2016, *Eur J Pain* — An evidence-based review of systemic treatments for itch. [PMID 26416344](https://pubmed.ncbi.nlm.nih.gov/26416344/) → 날트렉손·설트랄린·항히스타민의 근거 수준 → `ENTX`, `EAH` 를 작게 둔 이유
158. Anderson SL, 2026, *Drugs Context* — New and emerging treatments for PBC-related pruritus. [PMID 42232646](https://pubmed.ncbi.nlm.nih.gov/42232646/) → 가려움 축 신약 파이프라인
159. Levy C, 2023, *Clin Gastroenterol Hepatol* — New Treatment Paradigms in Primary Biliary Cholangitis. [PMID 36809835](https://pubmed.ncbi.nlm.nih.gov/36809835/) → 담즙정체 치료의 두 축 구조 (다른 질환에서의 확인)
160. Dong XR, 2023, *World J Clin Cases* — Effect of polyene phosphatidylcholine/UDCA/ademetionine on pregnancy outcomes in ICP. [PMID 37900240](https://pubmed.ncbi.nlm.nih.gov/37900240/) → 병용요법
161. Liao E, 2025, *J Mol Histol* — Inhibition of the GPR30-PI3K pathway by 4-phenylbutyric acid in the treatment of ICP. [PMID 40063113](https://pubmed.ncbi.nlm.nih.gov/40063113/) → 새로운 표적 (모델에 미포함)
162. Zeng Z, 2025, *J Reprod Immunol* — Probiotic VSL#3 alleviates ICP by upregulating FXR-FGF15. [PMID 41066868](https://pubmed.ncbi.nlm.nih.gov/41066868/) → 장내 미생물 개입 (모델의 `KDH_*`, `EC50FXI` 를 표적으로 삼는 접근)

---

## 11. 비타민 K · 응고 · 산후 출혈

163. **Cemortan M, 2025, *BMC Pregnancy Childbirth*** — Comparative analysis of vitamin K levels in women with ICP. [PMID 40200185](https://pubmed.ncbi.nlm.nih.gov/40200185/) → **`MIC50`, `KVK`, `WINR` (비타민 K 상태와 응고)** 의 사람 데이터
164. Arslanoğlu T, 2025, *BMC Pregnancy Childbirth* — ICP and coagulation: a dual risk of hypercoagulability and bleeding. [PMID 40281473](https://pubmed.ncbi.nlm.nih.gov/40281473/) → **ICP 의 응고 이상이 양방향**이라는 관찰; 모델은 출혈 쪽만 표현 — **명시적 누락**
165. Mao J, 2025, *J Glob Health* — Association between perinatal complications and venous thromboembolism in postpartum women. [PMID 40375726](https://pubmed.ncbi.nlm.nih.gov/40375726/) → 혈전 쪽 위험 (모델 미포함)

---

## 12. 쌍태 임신 · 재발 · 장기 예후

166. **Axelsen SM, 2024, *Acta Obstet Gynecol Scand*** — The effect of twin pregnancy in ICP: a case control study. [PMID 39058263](https://pubmed.ncbi.nlm.nih.gov/39058263/) → **`TWIN = 1.55`** 의 근거
167. Gu L, 2026, *J Matern Fetal Neonatal Med* — Maternal and neonatal outcomes of ICP in twin pregnancies: systematic review and meta-analysis. [PMID 42392869](https://pubmed.ncbi.nlm.nih.gov/42392869/) → 쌍태 ICP 의 결과
168. Xu T, 2022, *BMC Pregnancy Childbirth* — Perinatal outcomes associated with ICP in twin pregnancies were worse than singletons. [PMID 36335293](https://pubmed.ncbi.nlm.nih.gov/36335293/) → 같은 목적
169. Zhao Y, 2025, *BMC Pregnancy Childbirth* — Preterm birth and stillbirth: total bile acid levels in ICP and outcomes of twin pregnancies. [PMID 40389846](https://pubmed.ncbi.nlm.nih.gov/40389846/) → 쌍태에서의 층화
170. Mitta K, 2023, *Case Rep Womens Health* — Selective feticide reverses ICP in twins discordant for growth. [PMID 37534193](https://pubmed.ncbi.nlm.nih.gov/37534193/) → **태반 부하를 줄이면 질환이 역전된다**는 관찰; 모델에서 `TWIN` 이 성호르몬 궤적에만 작용하는 구조를 지지하는 가장 직접적인 증거
171. Rosenberg HM, 2026, *Obstet Gynecol* — ICP Recurrence in a Subsequent Pregnancy. [PMID 40811826](https://pubmed.ncbi.nlm.nih.gov/40811826/) → 재발률; 모델에서 재발은 유전 소인 벡터가 고정되어 있다는 결과
172. Zloto K, 2026, *Int J Gynaecol Obstet* — Outcomes of subsequent pregnancy following ICP. [PMID 42460752](https://pubmed.ncbi.nlm.nih.gov/42460752/) → 같은 목적
173. Sarker MR, 2024, *Am J Perinatol* — History of cholestasis is not associated with worsening outcomes in subsequent pregnancy. [PMID 38423120](https://pubmed.ncbi.nlm.nih.gov/38423120/) → 재발 시 중증도가 반드시 악화하지 않음
174. **Marschall HU, 2013, *Hepatology*** — ICP and associated hepatobiliary disease: a population-based cohort study. [PMID 23564560](https://pubmed.ncbi.nlm.nih.gov/23564560/) → 장기 간담도 질환 위험 (지도의 `E_LATER`)
175. Wikström Shemer EA, 2015, *J Hepatol* — ICP and cancer, immune-mediated and cardiovascular diseases: population-based cohort study. [PMID 25772037](https://pubmed.ncbi.nlm.nih.gov/25772037/) → 장기 예후의 범위
176. Odutola PO, 2023, *ILIVER* — ICP is associated with increased risk of hepatobiliary disease and adverse fetal outcomes: systematic review and meta-analysis. [PMID 40636920](https://pubmed.ncbi.nlm.nih.gov/40636920/) → 같은 목적
177. Wei Y, 2025, *Clin Rheumatol* — Postpartum may be a risk factor for biochemical flares in PBC. [PMID 40670882](https://pubmed.ncbi.nlm.nih.gov/40670882/) → 분만 후 담즙정체 질환의 거동
178. Dai S, 2026, *Arch Gynecol Obstet* — ICP and offspring neurodevelopment: bile acid-mediated mechanisms and long-term outcomes. [PMID 41949636](https://pubmed.ncbi.nlm.nih.gov/41949636/) → 자녀의 장기 신경발달 (모델 미포함 — **확장 후보**)
179. Guner Yilmaz B, 2025, *Children (Basel)* — Early metabolic profile in neonates with maternal ICP. [PMID 41462795](https://pubmed.ncbi.nlm.nih.gov/41462795/) → 신생아 대사

---

## 13. 동반질환 · 조절 인자

180. Han Z, 2026, *Hepatol Commun* — Clinical outcomes of ICP with versus without chronic hepatitis B. [PMID 41493830](https://pubmed.ncbi.nlm.nih.gov/41493830/) → HBV 동반 (모델 미포함)
181. Gao Q, 2024, *BMC Pregnancy Childbirth* — ICP combined with different stages of HBV infection on pregnancy outcomes. [PMID 38582906](https://pubmed.ncbi.nlm.nih.gov/38582906/) → 같은 목적
182. Li X, 2024, *Diabetol Metab Syndr* — Gestational diabetes aggravates adverse perinatal outcomes in women with ICP. [PMID 38429774](https://pubmed.ncbi.nlm.nih.gov/38429774/) → GDM 동반
183. Yang C, 2026, *Medicine (Baltimore)* — Association of GDM with ICP: Mendelian randomization. [PMID 41861211](https://pubmed.ncbi.nlm.nih.gov/41861211/) → 인과 방향
184. Lv M, 2026, *J Glob Health* — Impact of hypothyroidism on the risk of ICP. [PMID 41524244](https://pubmed.ncbi.nlm.nih.gov/41524244/) → 갑상선 기능과의 관계
185. Kwon Y, 2025, *J Pharm Pharm Sci* — ICP associated with azathioprine: disproportionality analysis. [PMID 41477595](https://pubmed.ncbi.nlm.nih.gov/41477595/) → 약물 유발 ICP (`GBSEP` 를 약리적으로 낮추는 상황)
186. Misirlioglu R, 2026, *Medicine (Baltimore)* — Prevalence and clinical significance of proteinuria in ICP. [PMID 42065144](https://pubmed.ncbi.nlm.nih.gov/42065144/) → 단백뇨 동반
187. Watad H, 2024, *PLoS One* — Proteinuria is a clinical characteristic of ICP but not a marker of severity. [PMID 39259746](https://pubmed.ncbi.nlm.nih.gov/39259746/) → 같은 주제, 중증도 지표로는 부적합

---

## 14. 진단 예측 모델 (모델의 비교 대상)

188. Han S, 2026, *J Matern Fetal Neonatal Med* — Innovative nomogram integrating bile acid metabolomics for early diagnosis of ICP. [PMID 42161851](https://pubmed.ncbi.nlm.nih.gov/42161851/) → **담즙산 메타볼로믹스(종별 측정)가 진단을 개선한다** → 이 QSP 모델이 종별로 계산하는 이유의 임상적 확인
189. Ding F, 2026, *BMC Pregnancy Childbirth* — Nomogram integrating dynamic total bile acid monitoring for preterm birth prediction in twin pregnancies. [PMID 42104311](https://pubmed.ncbi.nlm.nih.gov/42104311/) → **연속 측정이 단일 측정보다 낫다** → 모델 G 절이 스스로 지적한 한계(매끄러운 궤적)의 근거
190. Xie W, 2025, *BMC Pregnancy Childbirth* — Nomogram for predicting preterm birth in ICP. [PMID 39984873](https://pubmed.ncbi.nlm.nih.gov/39984873/) → 조산 예측
191. Wen M, 2025, *Medicine (Baltimore)* — Clinical prediction models for preterm birth in ICP. [PMID 41465959](https://pubmed.ncbi.nlm.nih.gov/41465959/) → 같은 목적
192. Yılmaz EBS, 2025, *BMC Pregnancy Childbirth* — Diagnostic and prognostic value of APRI and de Ritis ratio in ICP. [PMID 40877816](https://pubmed.ncbi.nlm.nih.gov/40877816/) → 간효소 비율 지표
193. Kahraman NÇ, 2026, *J Perinat Med* — FAR, PAR, APRI and adverse neonatal outcomes in ICP. [PMID 41351215](https://pubmed.ncbi.nlm.nih.gov/41351215/) → 같은 목적
194. Gregorc P, 2025, *Diagnostics (Basel)* — Evaluating the effect of bile acid levels on maternal and perinatal outcomes in ICP. [PMID 40941672](https://pubmed.ncbi.nlm.nih.gov/40941672/) → 담즙산-결과 관계의 외부 검증

---

## 15. QSP 방법론 · mrgsolve

195. **Elmokadem A, 2019, *CPT Pharmacometrics Syst Pharmacol*** — Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/) → 이 저장소가 따르는 구현 관례
196. Lu T, 2024, *CPT Pharmacometrics Syst Pharmacol* — gPKPDviz: A flexible R shiny tool for PK/PD simulations using mrgsolve. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/) → Shiny 앱 구조의 참조
197. Nigam SK, 2020, *Clin Pharmacol Ther* — The Systems Biology of Drug Metabolizing Enzymes and Transporters: Relevance to QSP. [PMID 32119114](https://pubmed.ncbi.nlm.nih.gov/32119114/) → 수송체 중심 QSP 의 방법론적 근거
198. Generaux G, 2019, *Pharmacol Res Perspect* — QST reproduces species differences in PF-04895162 liver safety due to combined mitochondrial and bile acid toxicity. [PMID 31624633](https://pubmed.ncbi.nlm.nih.gov/31624633/) → 담즙산 독성을 QSP 로 다룬 선례
199. Woodhead JL, 2024, *Xenobiotica* — Prediction of the liver safety profile of a first-in-class myeloperoxidase inhibitor using QST modeling. [PMID 38874513](https://pubmed.ncbi.nlm.nih.gov/38874513/) → 같은 계열
200. Mayo AK, 2026, *Clin Pharmacol Ther* — QST Model Predicts Obeticholic Acid-Associated Liver Injury in MASLD. [PMID 42332345](https://pubmed.ncbi.nlm.nih.gov/42332345/) → 담즙산 수용체 작용제의 간손상 예측
201. Beaudoin JJ, 2023, *Int J Mol Sci* — Combination of a human biomimetic liver microphysiology system with a QST modeling platform. [PMID 37298645](https://pubmed.ncbi.nlm.nih.gov/37298645/) → in vitro–모델 결합
202. Nyholm I, 2025, *J Hepatol* — Accumulation of altered serum bile acids predicts liver injury after portoenterostomy in biliary atresia. [PMID 39889904](https://pubmed.ncbi.nlm.nih.gov/39889904/) → **담즙산 조성(종별)이 총량보다 예측력이 높다**는 다른 질환에서의 확인

---

## 모델이 문헌과 어긋나거나 문헌 없이 정한 것들 (정직한 목록)

이 절은 검토자가 가장 먼저 읽어야 하는 부분입니다.

| 항목 | 상태 |
|------|------|
| **분만 시기 40–99 µmol/L 구간** | 모델은 39–40주를 최적으로 계산하며, 지침의 37–38주 권고를 **재현하지 못합니다**. `icp_calibration.py` G 절에 그 이유(Ovadia 의 해당 구간 초과 위험이 0.15%p 에 불과)와, 37주가 최적이 되기 위해 필요한 신생아 이환 가중치(<0.01)를 계산해 두었습니다. 이것은 모델의 결함일 수도 있고 권고의 근거가 사산 숫자에 있지 않다는 뜻일 수도 있으며, 어느 쪽인지 모델만으로는 판정할 수 없습니다 |
| **`VP` (태반 수송 최대능)** | 사람에서 직접 측정한 값을 찾지 못했습니다. 정상 임신 제대혈 농도에 맞춘 **모델 선택** |
| **`WTOX` 세포독성 가중치** | 소수성 서열과 심근세포 실험(문헌 93, 95, 124)의 **순서**는 근거가 있으나 **간격**은 모델 선택 |
| **`HLP` = 4 (LPA→가려움 Hill)** | 단일 결합 사건의 협동성이 아니라 LPAR–TRPV1–GRPR 다단계 증폭의 겉보기 값. 기전적으로 해석하지 말 것 |
| **UDCA 의 태아 심장 직접 보호** | 문헌 96, 97 은 UDCA 가 태아 심근을 직접 보호한다고 보고합니다. 모델은 이를 `WTOX[UDCA]=0.02` 로만 표현하고 별도 항을 두지 않았으므로 **UDCA 의 태아 보호를 과소평가할 가능성**이 있습니다 |
| **프로게스테론 설페이트의 장간순환** | 문헌 52 가 보고하지만 모델은 혈중 상태변수로만 두었습니다 (단순화) |
| **CDCA↔UDCA 역방향 에피머화** | 문헌 145 가 보고하지만 모델에는 UDCA→LCA 만 있습니다 (단순화) |
| **ICP 의 혈전 위험** | 문헌 164, 165 는 출혈과 혈전 양쪽 위험을 보고합니다. 모델은 출혈(비타민 K)만 표현합니다 |
| **자녀의 장기 신경발달** | 문헌 178. 모델 범위 밖 |
| **MRGPRX4 유전형** | 문헌 118. 가려움 개인차의 유력한 요인이나 모델에 미포함 |
| **모체 담즙산의 일내·일간 변동** | 모델은 매끄러운 궤적을 돌립니다. 문헌 189 는 연속 측정이 더 낫다고 보고하며, 이것이 모델의 가장 중요한 확장 방향입니다 |
| **적합된 상수의 총수** | `HSB0`, `HSBSC`, `HN` (사산, 문헌 14 의 세 값에 적합) · `HPT0`, `HM0`, `HMSC` (배경률) · `VSCALE` (VAS 척도) · `WMORB` (신생아 이환 대 사산의 효용비, 판단값) — 이상 8개. 나머지는 모두 수송·결합·회전 상수이거나 정상 임신이 정상값에 오도록 정한 척도입니다 |

---

## 파일 안내

| 파일 | 내용 |
|------|------|
| `icp_qsp_model.dot` / `.svg` / `.png` | 20개 클러스터 기계론적 지도 |
| `icp_reference_model.py` | 의존성 없는 Python RK4 참조 구현 (78 ODE) — 모든 방정식이 먼저 여기서 실행되었습니다 |
| `icp_calibration.py` | 적합·절제·시험 재현 분석 (A–I 절) |
| `icp_calibration_output.txt` | 위 스크립트의 실행 결과 |
| `icp_mrgsolve_model.R` | mrgsolve 모델 (78 ODE, 35 시나리오) |
| `icp_shiny_app.R` | 10탭 Shiny 대시보드 |
| `README.md` | 모델의 구조와 결과 요약 |
