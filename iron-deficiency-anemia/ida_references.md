# 철결핍성 빈혈 (Iron Deficiency Anaemia) — QSP 모델 참고문헌

이 문서는 `ida_qsp_model.dot` / `ida_mrgsolve_model.R` / `ida_shiny_app.R`의
구조와 파라미터 근거가 되는 문헌을 섹션별로 정리한 것입니다. 각 항목의
`→` 표시는 그 문헌이 모델의 **어느 부분**을 지지하는지 나타냅니다.

> ⚠️ 본 모델은 교육·연구 목적입니다. 파라미터는 문헌에 기반한 근사치이며
> 개별 환자 진료에 직접 사용할 수 없습니다.

---

## 1. 역학 · 진단 · 임상 개괄 (Epidemiology, diagnosis, overview)

1. Camaschella C. **Iron deficiency.** *Blood* 2019;133(1):30-39.
   <https://pubmed.ncbi.nlm.nih.gov/30401704/>
   → 전체 병태생리 골격, 절대적/기능적 철결핍 구분.
2. Camaschella C. **Iron-deficiency anemia.** *N Engl J Med* 2015;372(19):1832-1843.
   <https://pubmed.ncbi.nlm.nih.gov/25946282/>
   → 진단 알고리즘, 페리틴 <15 ng/mL 기준, 치료 반응 기대치.
3. Pasricha SR, Tye-Din J, Muckenthaler MU, Swinkels DW. **Iron deficiency.**
   *Lancet* 2021;397(10270):233-248.
   <https://pubmed.ncbi.nlm.nih.gov/33285139/>
   → 원인별 분류, 경구/정맥 철제 선택 기준.
4. Kassebaum NJ, et al. **A systematic analysis of global anemia burden from 1990
   to 2010.** *Blood* 2014;123(5):615-624.
   <https://pubmed.ncbi.nlm.nih.gov/24297872/>
   → 전 세계 빈혈 부담의 절반 가까이가 철결핍이라는 규모 근거.
5. Peyrin-Biroulet L, Williet N, Cacoub P. **Guidelines on the diagnosis and
   treatment of iron deficiency across indications: a systematic review.**
   *Am J Clin Nutr* 2015;102(6):1585-1594.
   <https://pubmed.ncbi.nlm.nih.gov/26561626/>
   → 적응증별 진단 역치(페리틴·TSAT)의 이질성.
6. Weiss G, Ganz T, Goodnough LT. **Anemia of inflammation.**
   *Blood* 2019;133(1):40-50.
   <https://pubmed.ncbi.nlm.nih.gov/30401705/>
   → 염증 시나리오(IL6_IN)의 근거: 페리틴↑ + TSAT↓ 해리 현상.

## 2. 철 흡수 생리 — 장세포 게이트 (Intestinal absorption)

7. Frazer DM, Anderson GJ. **The regulation of iron transport.**
   *Biofactors* 2014;40(2):206-214.
   <https://pubmed.ncbi.nlm.nih.gov/24132807/>
   → DMT1–ferroportin 축, 장세포 페리틴의 역할 (A_ENT / A_EFT 구조).
8. Gunshin H, et al. **Cloning and characterization of a mammalian
   proton-coupled metal-ion transporter (DMT1).** *Nature* 1997;388(6641):482-488.
   <https://pubmed.ncbi.nlm.nih.gov/9242408/>
   → 정점막(apical) Fe²⁺ 수송체 동정 — VMAX_D / KM_D의 대상.
9. McKie AT, et al. **An iron-regulated ferric reductase associated with the
   absorption of dietary iron (DCYTB).** *Science* 2001;291(5509):1755-1759.
   <https://pubmed.ncbi.nlm.nih.gov/11230685/>
   → Fe³⁺→Fe²⁺ 환원 단계(위산·아스코르브산 의존성).
10. Donovan A, et al. **The iron exporter ferroportin/Slc40a1 is essential for
    iron homeostasis.** *Cell Metab* 2005;1(3):191-200.
    <https://pubmed.ncbi.nlm.nih.gov/16054062/>
    → 기저측(basolateral) 배출의 유일한 경로 = FPN_ENT 인자의 정당성.
11. Vulpe CD, et al. **Hephaestin, a ceruloplasmin homologue implicated in
    intestinal iron transport.** *Nat Genet* 1999;21(2):195-199.
    <https://pubmed.ncbi.nlm.nih.gov/9988272/>
    → 트랜스페린 적재를 위한 재산화 단계.
12. Brasse-Lagnel C, et al. **Intestinal DMT1 cotransporter is down-regulated
    by hepcidin via proteasome internalization.**
    *Gastroenterology* 2011;140(4):1261-1271.
    <https://pubmed.ncbi.nlm.nih.gov/21199652/>
    → 헵시딘의 정점막 직접 작용 = 파라미터 `K_APIC`.
13. Hahn PF, Bale WF, Ross JF, Balfour WM, Whipple GH. **Radioactive iron
    absorption by gastro-intestinal tract: influence of anemia, anoxia, and
    antecedent feeding.** *J Exp Med* 1943;78(3):169-188.
    <https://pubmed.ncbi.nlm.nih.gov/19871320/>
    → "점막 차단(mucosal block)"의 원전 — 본 모델 THESIS 1의 역사적 근거.
14. Hurrell R, Egli I. **Iron bioavailability and dietary reference values.**
    *Am J Clin Nutr* 2010;91(5):1461S-1467S.
    <https://pubmed.ncbi.nlm.nih.gov/20200263/>
    → 식이 철 생체이용률 (`F_DIET` = 0.30의 근거 범위).
15. Lynch S, et al. **Biomarkers of Nutrition for Development (BOND) — Iron
    review.** *J Nutr* 2018;148(suppl_1):1001S-1067S.
    <https://pubmed.ncbi.nlm.nih.gov/29878148/>
    → 흡수율·체내 철 분포·바이오마커 기준값 종합.

## 3. 헵시딘 — 조절자 (Hepcidin: the controller)

16. Ganz T. **Systemic iron homeostasis.** *Physiol Rev* 2013;93(4):1721-1741.
    <https://pubmed.ncbi.nlm.nih.gov/24137020/>
    → 전신 철 흐름의 정량적 지도: 재활용 20-25 mg/day vs 흡수 1-2 mg/day.
    본 모델 THESIS 3의 핵심 근거.
17. Nemeth E, et al. **Hepcidin regulates cellular iron efflux by binding to
    ferroportin and inducing its internalization.**
    *Science* 2004;306(5704):2090-2093.
    <https://pubmed.ncbi.nlm.nih.gov/15514116/>
    → 헵시딘→ferroportin 분해 = `KDEG_FPE` / `KDEG_FPR`의 기전.
18. Nemeth E, et al. **IL-6 mediates hypoferremia of inflammation by inducing
    the synthesis of the iron regulatory hormone hepcidin.**
    *J Clin Invest* 2004;113(9):1271-1276.
    <https://pubmed.ncbi.nlm.nih.gov/15124018/>
    → `EMAX_IL6` / `KM_IL6` (STAT3 경로).
19. Andriopoulos B Jr, et al. **BMP6 is a key endogenous regulator of hepcidin
    expression and iron metabolism.** *Nat Genet* 2009;41(4):482-487.
    <https://pubmed.ncbi.nlm.nih.gov/19252486/>
    → 저장철 감지(BMP6/SMAD) = `HEP_STORE_E` / `KM_HEP_LIV`.
20. Kautz L, et al. **Identification of erythroferrone as an erythroid regulator
    of iron metabolism.** *Nat Genet* 2014;46(7):678-684.
    <https://pubmed.ncbi.nlm.nih.gov/24880340/>
    → ERFE 상태변수와 `KI_ERFE`의 근거.
21. Arezes J, et al. **Erythroferrone inhibits the induction of hepcidin by BMP6.**
    *Blood* 2018;132(14):1473-1477.
    <https://pubmed.ncbi.nlm.nih.gov/30097509/>
    → ERFE가 BMP를 포획하는 방식 (모델에서 f_erfe 곱셈항).
22. Finberg KE, et al. **Mutations in TMPRSS6 cause iron-refractory iron
    deficiency anemia (IRIDA).** *Nat Genet* 2008;40(5):569-571.
    <https://pubmed.ncbi.nlm.nih.gov/18408718/>
    → `TMPRSS6` 파라미터와 IRIDA 시나리오.
23. Silvestri L, et al. **The serine protease matriptase-2 (TMPRSS6) inhibits
    hepcidin activation by cleaving membrane hemojuvelin.**
    *Cell Metab* 2008;8(6):502-511.
    <https://pubmed.ncbi.nlm.nih.gov/18976966/>
    → matriptase-2 → HJV 절단 기전.
24. Ganz T, et al. **Immunoassay for human serum hepcidin.**
    *Blood* 2008;112(10):4292-4297.
    <https://pubmed.ncbi.nlm.nih.gov/18689548/>
    → 혈청 헵시딘 정상치 범위 (모델 healthy 6.4 ng/mL, IDA 0.32 ng/mL 대조).
25. Girelli D, Nemeth E, Swinkels DW. **Hepcidin in the diagnosis of iron
    disorders.** *Blood* 2016;127(23):2809-2813.
    <https://pubmed.ncbi.nlm.nih.gov/27044621/>
    → 헵시딘이 경구 철 반응성을 예측한다는 임상적 활용(모델 NOTE E).

## 4. 투여 간격 · 경구 철 약동학 (Dosing interval — the central thesis)

26. Moretti D, et al. **Oral iron supplements increase hepcidin and decrease
    iron absorption from daily or twice-daily doses in iron-depleted young
    women.** *Blood* 2015;126(17):1981-1989.
    <https://pubmed.ncbi.nlm.nih.gov/26289639/>
    → 60 mg 투여 후 헵시딘 상승이 24시간 지속되며 같은 날 두 번째 용량의
    흡수를 감소시킨다는 직접 근거. 모델 검증: probe +4 h = 75.4 %, +24 h = 91.5 %.
27. Stoffel NU, et al. **Iron absorption from oral iron supplements given on
    consecutive versus alternate days and as single morning doses versus
    twice-daily split dosing in iron-depleted women: two open-label,
    randomised controlled trials.** *Lancet Haematol* 2017;4(11):e524-e533.
    <https://pubmed.ncbi.nlm.nih.gov/29032957/>
    → 격일 투여의 분획흡수율 우위. 관찰된 daily:alternate 비 ≈ 0.75,
    본 모델 0.87 (모델이 보수적 — 한계로 명시).
28. Stoffel NU, et al. **Iron absorption from supplements is greater with
    alternate day than with consecutive day dosing in iron-deficient anemic
    women.** *Haematologica* 2020;105(5):1232-1239.
    <https://pubmed.ncbi.nlm.nih.gov/31413088/>
    → 빈혈이 있는(비단순 결핍) 여성에서도 동일 결론 — 본 모델의 IDA
    기저상태에서 시뮬레이션한 근거.
29. Kaundal R, Bhatia P, Jain A, et al. **Randomized controlled trial of
    twice-daily versus alternate-day oral iron therapy in the treatment of
    iron-deficiency anemia.** *Ann Hematol* 2020;99(1):57-63.
    <https://pubmed.ncbi.nlm.nih.gov/31811360/>
    → 임상 결과 수준의 검증: 격일 투여가 1일 2회 투여와 유사한 Hb 반응을
    내면서 부작용은 적었다. 총 흡수량과 분획흡수율이 다른 순서를 가진다는
    THESIS 2의 임상적 대응.
30. Rimon E, et al. **Are we giving too much iron? Low-dose iron therapy is
    effective in octogenarians.** *Am J Med* 2005;118(10):1142-1147.
    <https://pubmed.ncbi.nlm.nih.gov/16194646/>
    → 저용량(15-50 mg)이 고용량과 유사한 Hb 반응 — 모델의 분획흡수율
    포화(30 mg에서 FIA 25.2 %) 예측과 부합.
31. Schrier SL. **So you know how to treat iron deficiency anemia.**
    *Blood* 2015;126(17):1971.
    <https://pubmed.ncbi.nlm.nih.gov/26494915/>
    → Moretti 연구에 대한 논평: 임상 처방 관행 전환의 논거.
32. Zimmermann MB, Hurrell RF. **Nutritional iron deficiency.**
    *Lancet* 2007;370(9586):511-520.
    <https://pubmed.ncbi.nlm.nih.gov/17693180/>
    → 아스코르브산·피틴산 등 식이 상호작용.
33. Tolkien Z, Stecher L, Mander AP, Pereira DI, Powell JJ. **Ferrous sulfate
    supplementation causes significant gastrointestinal side-effects in
    adults: a systematic review and meta-analysis.**
    *PLoS One* 2015;10(2):e0117383.
    <https://pubmed.ncbi.nlm.nih.gov/25700159/>
    → GI 부작용·순응도 항(`K_GI`, `EMAX_ADH`)의 근거.
34. Gasche C, et al. **Ferric maltol is effective in correcting iron deficiency
    anemia in patients with inflammatory bowel disease (AEGIS trials).**
    *Inflamm Bowel Dis* 2015;21(3):579-588.
    <https://pubmed.ncbi.nlm.nih.gov/25545376/>
    → 대체 경구 제형(ferric maltol) — 지도의 oral pharmacology 클러스터.

## 5. 정맥 철제 — 약동학과 안전성 (Intravenous iron)

35. Auerbach M, Macdougall I. **The available intravenous iron formulations:
    History, efficacy, and toxicology.** *Hemodial Int* 2017;21 Suppl 1:S83-S92.
    <https://pubmed.ncbi.nlm.nih.gov/28371203/>
    → 제형별 탄수화물 껍질·RES 처리·유리철 방출 (A_COL / F_COL_RES / K_COL).
36. Geisser P, Burckhardt S. **The pharmacokinetics and pharmacodynamics of
    iron preparations.** *Pharmaceutics* 2011;3(1):12-33.
    <https://pubmed.ncbi.nlm.nih.gov/24310424/>
    → FCM 혈중 반감기 7-12 h (`K_COL` = 0.0693/h ≈ t½ 10 h).
37. Onken JE, et al. **A multicenter, randomized, active-controlled study to
    investigate the efficacy and safety of ferric carboxymaltose in patients
    with iron deficiency anemia (REPAIR-IDA).**
    *Transfusion* 2014;54(2):306-315.
    <https://pubmed.ncbi.nlm.nih.gov/23772856/>
    → FCM 750 mg×2 요법의 Hb 반응 크기 (모델 시나리오 9).
38. Van Wyck DB, Mangione A, Morrison J, Hadley PE, Jehle JA, Goodnough LT.
    **Large-dose intravenous ferric carboxymaltose injection for iron
    deficiency anemia in heavy uterine bleeding.**
    *Transfusion* 2009;49(12):2719-2728.
    <https://pubmed.ncbi.nlm.nih.gov/19682342/>
    → 자궁 출혈에 의한 IDA에서 대용량 FCM — 본 모델의 기저 환자상과 동일.
39. Rognoni C, Venturini S, Meregaglia M, Marmifero M, Tarricone R.
    **Efficacy and safety of ferric carboxymaltose and other formulations in
    iron-deficient patients: a systematic review and network meta-analysis.**
    *Clin Drug Investig* 2016;36(3):177-194.
    <https://pubmed.ncbi.nlm.nih.gov/26692005/>
    → 제형 간 Hb 반응의 유사성 — THESIS 4("동등한 효능, 상이한 대가")의 전제.
40. Rampton D, et al. **Hypersensitivity reactions to intravenous iron:
    guidance for risk minimization and management.**
    *Haematologica* 2014;99(11):1671-1676.
    <https://pubmed.ncbi.nlm.nih.gov/25420283/>
    → CARPA·Fishbane 반응 (지도의 안전성 노드).

## 6. FGF23 — 인 축 (The carboxymaltose phosphate penalty)

41. Wolf M, et al. **Effects of iron isomaltoside vs ferric carboxymaltose on
    hypophosphatemia in iron-deficiency anemia: two randomized clinical trials
    (PHOSPHARE-IDA).** *JAMA* 2020;323(5):432-443.
    <https://pubmed.ncbi.nlm.nih.gov/32016310/>
    → 본 모델 FGF23-인 축의 주 보정 근거: FCM에서 인 <2.0 mg/dL 발생률이
    현저히 높고 derisomaltose에서는 드물다. 모델: FCM 1000 mg → 인 최저
    1.98 mg/dL(6.8일), <2.0 mg/dL 4.2일; derisomaltose → 3.02 mg/dL, 0일.
42. Wolf M, Koch TA, Bregman DB. **Effects of iron deficiency anemia and its
    treatment on fibroblast growth factor 23 and phosphate homeostasis in
    women.** *J Bone Miner Res* 2013;28(8):1793-1803.
    <https://pubmed.ncbi.nlm.nih.gov/23505057/>
    → iFGF23 상승이 FCM 특이적이며 24-48시간에 정점 (`E_CLV`, `KOUT_CLV`).
43. Schouten BJ, Hunt PJ, Livesey JH, Frampton CM, Soule SG.
    **FGF23 elevation and hypophosphatemia after intravenous iron
    polymaltose: a prospective study.**
    *J Clin Endocrinol Metab* 2009;94(7):2332-2337.
    <https://pubmed.ncbi.nlm.nih.gov/19366850/>
    → 제형별 차이의 초기 관찰.
44. Zoller H, Schaefer B, Glodny B. **Iron-induced hypophosphatemia: an
    emerging complication.** *Curr Opin Nephrol Hypertens* 2017;26(4):266-275.
    <https://pubmed.ncbi.nlm.nih.gov/28399017/>
    → 절단 억제(cleavage inhibition) 가설 — 모델의 CLV 상태변수 구조.
45. Shimada T, et al. **FGF-23 is a potent regulator of vitamin D metabolism
    and phosphate homeostasis.** *J Bone Miner Res* 2004;19(3):429-435.
    <https://pubmed.ncbi.nlm.nih.gov/15040831/>
    → NaPi-2a/2c 억제 + CYP27B1 억제 (`IMAX_FGF_P`, `IMAX_FGF_D`).
46. Klein K, et al. **Severe FGF23-based hypophosphataemic osteomalacia due to
    ferric carboxymaltose administration.** *BMJ Case Rep* 2018;2018:bcr-2017-222851.
    <https://pubmed.ncbi.nlm.nih.gov/29298794/>
    → 장기 저인산혈증의 임상 결과(골연화증) — 지도의 OSTEOMAL 노드.

## 7. 적혈구 생성 · 철 제한 적혈구조혈 (Erythropoiesis and iron restriction)

47. Koury MJ, Ponka P. **New insights into erythropoiesis: the roles of folate,
    vitamin B12, and iron.** *Annu Rev Nutr* 2004;24:105-131.
    <https://pubmed.ncbi.nlm.nih.gov/15189115/>
    → 철 제한 시 적혈구조혈의 비효율성(apoptosis) — `KAPO_MAX`.
48. Rivella S. **Iron metabolism under conditions of ineffective
    erythropoiesis in beta-thalassemia.** *Blood* 2019;133(1):51-58.
    <https://pubmed.ncbi.nlm.nih.gov/30401707/>
    → 비효율적 적혈구조혈에서 철의 재순환(모델: apo_fe → A_RES).
49. Ponka P. **Tissue-specific regulation of iron metabolism and heme
    synthesis: distinct control mechanisms in erythroid cells.**
    *Blood* 1997;89(1):1-25.
    <https://pubmed.ncbi.nlm.nih.gov/8978272/>
    → ALAS2·미토페린·ferrochelatase 경로, ZPP 생성.
50. Shayeghi M, et al. **Identification of an intestinal heme transporter.**
    *Cell* 2005;122(5):789-801.
    <https://pubmed.ncbi.nlm.nih.gov/16143108/>
    → 헴철 흡수 경로 (지도의 HCP1 노드).
51. Brugnara C. **Iron deficiency and erythropoiesis: new diagnostic
    approaches.** *Clin Chem* 2003;49(10):1573-1578.
    <https://pubmed.ncbi.nlm.nih.gov/14500582/>
    → 망상적혈구 헤모글로빈 함량(CHr) — 모델 `CHR` 출력(IDA 20.5 pg).
52. Thomas C, Thomas L. **Biochemical markers and hematologic indices in the
    diagnosis of functional iron deficiency.**
    *Clin Chem* 2002;48(7):1066-1076.
    <https://pubmed.ncbi.nlm.nih.gov/12089176/>
    → sTfR / log(ferritin) 지표 (`STFR` 상태변수).
53. Beguin Y. **Soluble transferrin receptor for the evaluation of
    erythropoiesis and iron status.** *Clin Chim Acta* 2003;329(1-2):9-22.
    <https://pubmed.ncbi.nlm.nih.gov/12589962/>
    → sTfR 정상 범위 및 철결핍 시 상승폭 (`STFR0`, `STFR_E`).

## 8. 비적혈구 철 기능 · 증상 (Non-erythroid iron and symptoms)

54. Haas JD, Brownlie T 4th. **Iron deficiency and reduced work capacity:
    a critical review of the research to determine a causal relationship.**
    *J Nutr* 2001;131(2S-2):676S-690S.
    <https://pubmed.ncbi.nlm.nih.gov/11160598/>
    → 조직 철(myoglobin·Fe-S 효소)과 운동능력 — `W_FACIT_TISS`의 근거.
55. Allen RP, et al. **Restless legs syndrome/Willis-Ekbom disease
    pathophysiology.** *Sleep Med Clin* 2015;10(3):207-214.
    <https://pubmed.ncbi.nlm.nih.gov/26329429/>
    → 뇌 철 결핍-도파민 축 → IRLS 엔드포인트.
56. Krayenbuehl PA, Battegay E, Breymann C, Furrer J, Schulthess G.
    **Intravenous iron for the treatment of fatigue in nonanemic,
    premenopausal women with low serum ferritin concentration.**
    *Blood* 2011;118(12):3222-3227.
    <https://pubmed.ncbi.nlm.nih.gov/21705493/>
    → 빈혈이 없어도 조직 철 보충이 피로를 개선 — 모델 NOTE H(조직 철이
    Hb보다 늦게/다르게 회복)의 임상적 대응.
57. Falkingham M, et al. **The effects of oral iron supplementation on
    cognition in older children and adults: a systematic review and
    meta-analysis.** *Nutr J* 2010;9:4.
    <https://pubmed.ncbi.nlm.nih.gov/20100340/>
    → 인지 엔드포인트.
58. Lozoff B, Beard J, Connor J, Barbara F, Georgieff M, Schallert T.
    **Long-lasting neural and behavioral effects of iron deficiency in
    infancy.** *Nutr Rev* 2006;64(5 Pt 2):S34-S43.
    <https://pubmed.ncbi.nlm.nih.gov/16770951/>
    → 영아 신경발달 결과.

## 9. 적응증별 임상시험 (Indication-specific trials)

59. Anker SD, et al. **Ferric carboxymaltose in patients with heart failure and
    iron deficiency (FAIR-HF).** *N Engl J Med* 2009;361(25):2436-2448.
    <https://pubmed.ncbi.nlm.nih.gov/19920054/>
    → 빈혈 유무와 무관한 기능 개선 — 조직 철 축의 임상 근거.
60. Ponikowski P, et al. **Ferric carboxymaltose for iron deficiency at
    discharge after acute heart failure (AFFIRM-AHF).**
    *Lancet* 2020;396(10266):1895-1904.
    <https://pubmed.ncbi.nlm.nih.gov/33197395/>
    → 재입원 감소 엔드포인트.
61. Muñoz M, et al. **International consensus statement on the peri-operative
    management of anaemia and iron deficiency.**
    *Anaesthesia* 2017;72(2):233-247.
    <https://pubmed.ncbi.nlm.nih.gov/27996086/>
    → 수술 전 최적화·수혈 회피 엔드포인트.
62. Pavord S, et al. **UK guidelines on the management of iron deficiency in
    pregnancy.** *Br J Haematol* 2020;188(6):819-830.
    <https://pubmed.ncbi.nlm.nih.gov/31578718/>
    → 임신 중 요구량(+1000 mg)과 치료 역치.
63. Dignass AU, Gasche C, et al. **European consensus on the diagnosis and
    management of iron deficiency and anaemia in inflammatory bowel diseases.**
    *J Crohns Colitis* 2015;9(3):211-222.
    <https://pubmed.ncbi.nlm.nih.gov/25518052/>
    → 염증성 장질환에서 경구 철 회피 권고 — 모델 NOTE E의 임상적 대응.
64. Munro MG, Mast AE, Powers JM, et al. **The relationship between heavy
    menstrual bleeding, iron deficiency, and iron deficiency anemia.**
    *Am J Obstet Gynecol* 2023;229(1):1-9.
    <https://pubmed.ncbi.nlm.nih.gov/36706856/>
    → `VBLEED` 파라미터화의 임상 근거(>80 mL/cycle 정의, 철 손실 환산).

## 10. QSP · 모델링 방법론 (QSP methodology)

65. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and
    physiologically-based pharmacokinetic modeling with mrgsolve: a
    hands-on tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019;8(12):883-893.
    <https://pubmed.ncbi.nlm.nih.gov/31652028/>
    → mrgsolve 구현 관례.
66. Parmar JH, Mendes P. **A computational model to understand mouse iron
    physiology and disease.** *PLoS Comput Biol* 2019;15(1):e1006680.
    <https://pubmed.ncbi.nlm.nih.gov/30608934/>
    → 전신 철 대사 ODE 모델의 선행 사례(구획 구성 비교 근거).
67. Sarkar J, Potdar AA, Saidel GM. **Whole-body iron transport and metabolism:
    mechanistic, multi-scale model to improve treatment of anemia in chronic
    kidney disease.** *PLoS Comput Biol* 2018;14(4):e1006060.
    <https://pubmed.ncbi.nlm.nih.gov/29659573/>
    → 다중 규모 철 수송 모델 — 본 모델의 재활용/저장 구획 설계 비교.
68. Enculescu M, et al. **Modelling systemic iron regulation during dietary
    iron overload and acute inflammation: role of hepcidin-independent
    mechanisms.** *PLoS Comput Biol* 2017;13(1):e1005322.
    <https://pubmed.ncbi.nlm.nih.gov/28068331/>
    → 헵시딘 의존/비의존 기전의 분리 — 본 모델의 두 ferroportin 시계 설계.
69. Ganz T, Nemeth E. **Hepcidin and iron homeostasis.**
    *Biochim Biophys Acta* 2012;1823(9):1434-1443.
    <https://pubmed.ncbi.nlm.nih.gov/22306005/>
    → 정량적 파라미터 종합(헵시딘 반감기, ferroportin 회복 시간).

---

## 모델 보정 요약 (Calibration anchors and where they came from)

| 모델 예측값 | 관찰/문헌값 | 근거 문헌 |
|---|---|---|
| 60 mg 단회 분획흡수율 22.0 % | 철결핍 여성 약 20-22 % | 27, 28 |
| 헵시딘 24 h 후 1.36배 잔존 | 24 h까지 상승 지속 | 26 |
| +4 h probe 흡수 75.4 % | 같은 날 2회 투여 시 흡수 저하 | 26 |
| daily : alternate 분획흡수 비 0.87 | 관찰 0.75 (모델이 보수적) | 27 |
| 혈청철 정점 5.8 h | 2-5 h (모델이 다소 늦음) | 26 |
| 망상적혈구 정점 12-14일 | 7-14일 | 2 |
| FCM 1000 mg → 인 최저 1.98 mg/dL (6.8일) | 최저 ~1.9-2.0 mg/dL, 7-14일 | 41, 42 |
| derisomaltose → iFGF23 무변화 | FCM 특이적 현상 | 41, 44 |
| 골수 철 소비 13.3 mg/day (재활용 ~85 %) | 20-25 mg/day (체격 보정 시 부합) | 16 |
| 정상 헵시딘 6.4 / IDA 0.32 ng/mL | 정상 1-20, IDA 검출한계 이하 | 24, 25 |
| 질량보존 오차 ~1e-12 % | — | (내부 감사, NOTE I) |

## 명시적 한계 (Documented limitations)

1. **격일 투여 이득이 관찰값보다 작다.** 모델의 daily:alternate 분획흡수율
   비는 0.87이지만 Stoffel 2017의 관찰값은 약 0.75이다. 즉 본 모델은
   refractory penalty를 **과소평가**한다. 헵시딘의 순환 반감기(2.5 h)를
   고정한 채 장세포 배출능의 회복 시계만으로 24시간 기억을 만들면 깊이와
   지속성이 상충하기 때문이다.
2. **혈청철 정점이 늦다** (모델 5.8 h vs 관찰 2-5 h): 흡수 창(`KTR`)과
   기저측 배출이 직렬로 연결되어 있어 흡수 지연이 누적된다.
3. **IRIDA 표현형이 경증**이다. `TMPRSS6` = 0.30으로 얻은 기저 Hb는
   12.8 g/dL로, 실제 IRIDA의 중등도 미세적혈구 빈혈보다 가볍다. 경구
   불응성(oral +0.32 vs +2.59 g/dL)이라는 정성적 결론은 유지된다.
4. **GI 부작용–순응도 항은 반정량적**이다. `K_GI`, `EMAX_ADH`는 메타분석의
   부작용 발생률(33)과 정합하는 크기로 설정했을 뿐, 개별 환자의 중단률을
   예측하도록 보정된 것은 아니다.
5. **단일 환자 체격**(60 kg 여성)으로 보정되어 있다. 소아·임신·CKD 등은
   `BV_L`, `PV_L`, `VBLEED`, `IL6_IN`을 바꾼 뒤 **재평형**이 필요하다
   (Shiny 앱의 "기저상태 재평형" 옵션).
