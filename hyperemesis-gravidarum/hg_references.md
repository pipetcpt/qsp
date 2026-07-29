# Hyperemesis Gravidarum (HG) — QSP Model References

임신 오조 (hyperemesis gravidarum) QSP 모델의 근거 문헌.
모든 PMID는 NCBI E-utilities로 저자·저널·연도·제목을 대조 검증했습니다 (2026-07-29).

정렬 원칙: 모델의 **구조적 주장**을 지지하거나 반증할 수 있는 문헌을 먼저 배치했습니다.
각 섹션 머리말은 그 문헌군이 모델의 어느 방정식·파라미터에 들어갔는지 밝힙니다.

---

## 1. 핵심 기전 — GDF15·GFRAL 축과 "배수 변화(fold-change)" 가설

> 모델의 1번 구조적 가정(`SP`·`ALPHA`·`TAU_SP`, 클러스터 3)이 여기서 나옵니다.
> Fejzo 2024가 결정적입니다: 태아 유래 GDF15 **생산량**과 산모의 **감수성**이
> 각각 독립적으로 위험에 기여하고, 감수성은 임신 전 GDF15 노출로 정해지며,
> 마우스에서 GDF15 사전 노출은 **탈감작**을 유발합니다.

| # | 문헌 | PMID |
|---|------|------|
| 1 | Fejzo M, et al. **GDF15 linked to maternal risk of nausea and vomiting during pregnancy.** *Nature* 2024;625:760-767. | [38092039](https://pubmed.ncbi.nlm.nih.gov/38092039/) |
| 2 | Fejzo MS, et al. **Placenta and appetite genes GDF15 and IGFBP7 are associated with hyperemesis gravidarum.** *Nat Commun* 2018. | [29563502](https://pubmed.ncbi.nlm.nih.gov/29563502/) |
| 3 | Fejzo M, et al. **Multi-ancestry genome-wide association study of severe pregnancy nausea and vomiting.** *Nat Genet* 2026. | [41981316](https://pubmed.ncbi.nlm.nih.gov/41981316/) |
| 4 | Mullican SE, et al. **GFRAL is the receptor for GDF15 and the ligand promotes weight loss in mice and nonhuman primates.** *Nat Med* 2017. | [28846097](https://pubmed.ncbi.nlm.nih.gov/28846097/) |
| 5 | Breit SN, et al. **The GDF15-GFRAL Pathway in Health and Metabolic Disease: Friend or Foe?** *Annu Rev Physiol* 2021. | [33228454](https://pubmed.ncbi.nlm.nih.gov/33228454/) |
| 6 | Hes C, et al. **GDNF family receptor alpha-like (GFRAL) expression is restricted to the caudal brainstem.** *Mol Metab* 2025. | [39608751](https://pubmed.ncbi.nlm.nih.gov/39608751/) |
| 7 | Coll AP, et al. **GDF15 mediates the effects of metformin on body weight and energy balance.** *Nature* 2020. | [31875646](https://pubmed.ncbi.nlm.nih.gov/31875646/) |
| 8 | Alshaikh ABA, et al. **Hyperemesis gravidarum revisited: from GDF15 biology to precision multimodal therapy.** *Naunyn Schmiedebergs Arch Pharmacol* 2026. | [41942591](https://pubmed.ncbi.nlm.nih.gov/41942591/) |

**모델에 들어간 지점.** 4번·6번은 GFRAL이 area postrema/NTS에만 발현한다는 사실을
근거로, 배수 변화 검출기의 출력이 오직 후뇌 노드로만 들어가게 제한합니다(클러스터 4).
7번은 메트포르민 → 기저 GDF15 상승 경로(`MET_EMAX_GDF`)의 유일한 근거입니다.

---

## 2. 임신 중 GDF15의 시간 경과 — "농도는 계속 오르는데 증상은 사라진다"

> 모델이 설명해야 하는 가장 껄끄러운 관찰. 산모 GDF15는 1분기를 지나서도
> 떨어지지 않는데(오히려 만삭까지 오름), 증상은 임신 9-11주에 정점을 찍고
> 16-20주에 사라집니다. 농도 기반 모델로는 불가능합니다.

| # | 문헌 | PMID |
|---|------|------|
| 9 | Marjono AB, et al. **Macrophage inhibitory cytokine-1 in gestational tissues and maternal serum in normal and pre-eclamptic pregnancy.** *Placenta* 2003. | [12495665](https://pubmed.ncbi.nlm.nih.gov/12495665/) |
| 10 | Tong S, et al. **Serum concentrations of macrophage inhibitory cytokine 1 (MIC 1) as a predictor of miscarriage.** *Lancet* 2004. | [14726168](https://pubmed.ncbi.nlm.nih.gov/14726168/) |
| 11 | Kaitu'u-Lino TJ, et al. **Plasma MIC-1 and PAPP-A levels are decreased among women presenting to an early pregnancy assessment unit...** *PLoS One* 2013. | [24069146](https://pubmed.ncbi.nlm.nih.gov/24069146/) |
| 12 | Chen Q, et al. **Serum levels of GDF15 are reduced in preeclampsia...** *Cytokine* 2016. | [27173615](https://pubmed.ncbi.nlm.nih.gov/27173615/) |
| 13 | Lyu C, et al. **Insufficient GDF15 expression predisposes women to unexplained recurrent pregnancy loss by impairing extravillous trophoblast invasion.** *Cell Prolif* 2023. | [37272232](https://pubmed.ncbi.nlm.nih.gov/37272232/) |
| 14 | Yang SH, et al. **GDF15 promotes trophoblast invasion and pregnancy success via the BMPR1A/BMPR2/p-SMAD1 pathway.** *Life Sci* 2025. | [40157640](https://pubmed.ncbi.nlm.nih.gov/40157640/) |

**주의.** 13번·14번은 GDF15가 영양막 침습에 **필요**하다는 증거이고, 이것이
모델의 anti-GDF15 항체 시나리오에 붙은 태아 안전성 경고(`FETAL_SAF` 노드)의
근거입니다. 배수 변화 논리만 보면 리간드 차단이 최선이지만, 그 리간드는
1분기에 태반이 필요로 하는 것이기도 합니다.

---

## 3. 임신 전 GDF15 노출이 보호적이라는 세 가지 독립 증거

> 모델의 가장 강한 예측이 시험되는 지점. β-thalassemia(만성 상승),
> 메트포르민(약 2배 상승), 흡연(약 1/3 상승) — 세 노출 모두 **임신 전**
> 작용하며 모두 보호적이고, 상승 폭 순서대로 보호 강도가 정렬됩니다.
> 모델은 `ALPHA` 하나로 이 순서를 재현합니다(적합이 아니라 예측).

| # | 문헌 | PMID |
|---|------|------|
| 15 | Sharma N, et al. **Prepregnancy metformin use associated with lower risk of severe nausea and vomiting of pregnancy and hyperemesis gravidarum.** *Am J Obstet Gynecol* 2025. — 메트포르민 aRR 0.29 (95% CI 0.12-0.71), 흡연 aRR 0.51 (0.30-0.86), n=5414 | [40588059](https://pubmed.ncbi.nlm.nih.gov/40588059/) |
| 16 | Ranjbaran R, et al. **GDF-15 negatively regulates excess erythropoiesis and its overexpression is involved in erythroid hyperplasia.** *Exp Cell Res* 2020. | [33164866](https://pubmed.ncbi.nlm.nih.gov/33164866/) |
| 17 | Piolatto A, et al. **GDF15 as a Marker of Ineffective Erythropoiesis and Erythroid Expansion in Thalassemia: a Clinical Perspective.** *Clin Lab* 2026. | [42159121](https://pubmed.ncbi.nlm.nih.gov/42159121/) |

---

## 4. 정의·역학·자연사 — PUQE와 Windsor 기준

> `puqe24()` 함수의 구간 경계, HG 사례 정의(`hg_case()`), 그리고 자연사 적합
> 표적(정점 9-11주, 16주 완화)이 여기서 나옵니다.

| # | 문헌 | PMID |
|---|------|------|
| 18 | Fejzo MS, et al. **Nausea and vomiting of pregnancy and hyperemesis gravidarum.** *Nat Rev Dis Primers* 2019. | [31515515](https://pubmed.ncbi.nlm.nih.gov/31515515/) |
| 19 | Jansen LAW, et al. **The Windsor definition for hyperemesis gravidarum: a multistakeholder international consensus definition.** *Eur J Obstet Gynecol Reprod Biol* 2021. | [34555550](https://pubmed.ncbi.nlm.nih.gov/34555550/) |
| 20 | Koren G, et al. **Validation studies of the Pregnancy Unique-Quantification of Emesis (PUQE) scores.** *J Obstet Gynaecol* 2005. | [16147725](https://pubmed.ncbi.nlm.nih.gov/16147725/) |
| 21 | Koren G, et al. **Measuring the severity of nausea and vomiting of pregnancy; a 20-year perspective on the use of PUQE.** *J Obstet Gynaecol* 2021. | [32811235](https://pubmed.ncbi.nlm.nih.gov/32811235/) |
| 22 | Yilmaz T, et al. **Psychometric properties of the Pregnancy-Unique Quantification of Emesis (PUQE-24) Scale.** *J Obstet Gynaecol* 2022. | [35253594](https://pubmed.ncbi.nlm.nih.gov/35253594/) |
| 23 | Jansen LAW, et al. **Diagnosis and treatment of hyperemesis gravidarum.** *CMAJ* 2024. | [38621783](https://pubmed.ncbi.nlm.nih.gov/38621783/) |
| 24 | Ioannidou P, et al. **Predictive factors of Hyperemesis Gravidarum: a systematic review.** *Eur J Obstet Gynecol Reprod Biol* 2019. | [31126753](https://pubmed.ncbi.nlm.nih.gov/31126753/) |
| 25 | Beyene GA, et al. **Prevalence and determinants of hyperemesis gravidarum among pregnant women in Ethiopia: systematic review and meta-analysis.** *PLoS One* 2024. | [39625915](https://pubmed.ncbi.nlm.nih.gov/39625915/) |
| 26 | Niemeijer MN, et al. **Diagnostic markers for hyperemesis gravidarum: a systematic review and metaanalysis.** *Am J Obstet Gynecol* 2014. | [24530975](https://pubmed.ncbi.nlm.nih.gov/24530975/) |
| 27 | Vadakekut ES, et al. **Hyperemesis Gravidarum.** *StatPearls* (updated 2026). | [30422512](https://pubmed.ncbi.nlm.nih.gov/30422512/) |
| 28 | Nurmi M, et al. **Recurrence patterns of hyperemesis gravidarum.** *Am J Obstet Gynecol* 2018. | [30121224](https://pubmed.ncbi.nlm.nih.gov/30121224/) |
| 29 | Fassett MJ, et al. **Hyperemesis Gravidarum: Risk of Recurrence in Subsequent Pregnancies.** *Reprod Sci* 2023. | [36163577](https://pubmed.ncbi.nlm.nih.gov/36163577/) |
| 30 | Fejzo MS, et al. **Recurrence risk of hyperemesis gravidarum.** *J Midwifery Womens Health* 2011. | [21429077](https://pubmed.ncbi.nlm.nih.gov/21429077/) |

---

## 5. 진료 지침

| # | 문헌 | PMID |
|---|------|------|
| 31 | Nelson-Piercy C, et al. **The Management of Nausea and Vomiting in Pregnancy and Hyperemesis Gravidarum (RCOG Green-top Guideline No. 69).** *BJOG* 2024. | [38311315](https://pubmed.ncbi.nlm.nih.gov/38311315/) |
| 32 | Committee on Practice Bulletins—Obstetrics. **ACOG Practice Bulletin No. 189: Nausea and Vomiting of Pregnancy.** *Obstet Gynecol* 2018. | [29266076](https://pubmed.ncbi.nlm.nih.gov/29266076/) |
| 33 | Clark SM, et al. **Inpatient Management of Hyperemesis Gravidarum.** *Obstet Gynecol* 2024. | [38301258](https://pubmed.ncbi.nlm.nih.gov/38301258/) |
| 34 | Spinosa D, et al. **Management Considerations for Recalcitrant Hyperemesis.** *Obstet Gynecol Surv* 2020. | [31999353](https://pubmed.ncbi.nlm.nih.gov/31999353/) |

---

## 6. 무작위 대조시험 — 모델 검증의 정량 기준점

> 이 표의 숫자가 `hg_reference_impl.py`의 VALIDATION 표에 그대로 들어갑니다.
> **35번(VOMIT 시험)이 모델 전체의 출발점**입니다: 지침 권고 약물인
> 온단세트론이 −0.51(비유의)에 그쳤고 미르타자핀이 −1.86을 냈으며, 차이가
> 4일 이후 벌어졌습니다. 두 파라미터(`W_VAG`, `E0`)만 이 시험에 적합했고
> 나머지 시험(36-39번)은 예측입니다.

| # | 문헌 | 주요 수치 | PMID |
|---|------|-----------|------|
| 35 | Ostenfeld A, et al. **Mirtazapine or ondansetron for hyperemesis gravidarum: a randomized placebo-controlled trial (VOMIT).** *Am J Obstet Gynecol* 2026. | ΔPUQE-24 day 2: 미르타자핀 −1.86 (95% CI −3.61 to −0.12); 온단세트론 −0.51 (−2.32 to 1.30); n=59, 7개 병원, 덴마크 | [41478546](https://pubmed.ncbi.nlm.nih.gov/41478546/) |
| 36 | Koren G, et al. **Effectiveness of delayed-release doxylamine and pyridoxine for nausea and vomiting of pregnancy: a randomized placebo controlled trial.** *Am J Obstet Gynecol* 2010. | ΔPUQE −4.8 vs −3.9 (P=.006), 중등도 NVP | [20843504](https://pubmed.ncbi.nlm.nih.gov/20843504/) |
| 37 | Guttuso T Jr, et al. **Effect of gabapentin on hyperemesis gravidarum: a double-blind, randomized controlled trial.** *Am J Obstet Gynecol MFM* 2021. | day 5-7 PUQE 감소 52% 우월 (95% CI 16-88), n=21 | [33451591](https://pubmed.ncbi.nlm.nih.gov/33451591/) |
| 38 | Maina A, et al. **Transdermal clonidine in the treatment of severe hyperemesis: a pilot randomised control trial (CLONEMESI).** *BJOG* 2014. | PUQE 개선 95% CI 0.43-3.24; 수축기혈압 −6 mmHg | [24684734](https://pubmed.ncbi.nlm.nih.gov/24684734/) |
| 39 | Yost NP, et al. **A randomized, placebo-controlled trial of corticosteroids for hyperemesis due to pregnancy.** *Obstet Gynecol* 2003. | 재입원 34% vs 35% (P=.89) — **음성 결과**, n=110 | [14662211](https://pubmed.ncbi.nlm.nih.gov/14662211/) |
| 40 | Smith C, et al. **A randomized controlled trial of ginger to treat nausea and vomiting in pregnancy.** *Obstet Gynecol* 2004. | | [15051552](https://pubmed.ncbi.nlm.nih.gov/15051552/) |
| 41 | McParlin C, et al. **Treatments for Hyperemesis Gravidarum and Nausea and Vomiting in Pregnancy: a Systematic Review.** *JAMA* 2016. | | [27701665](https://pubmed.ncbi.nlm.nih.gov/27701665/) |
| 42 | O'Donnell A, et al. **Treatments for hyperemesis gravidarum and nausea and vomiting in pregnancy: a systematic review and economic assessment.** *Health Technol Assess* 2016. | | [27731292](https://pubmed.ncbi.nlm.nih.gov/27731292/) |
| 43 | Abramowitz A, et al. **Treatment options for hyperemesis gravidarum.** *Arch Womens Ment Health* 2017. | | [28070660](https://pubmed.ncbi.nlm.nih.gov/28070660/) |

---

## 7. 약동학·수용체 약리 — 노드 위치 법칙의 입력값

> 모델의 2번 구조적 가정. 각 약물의 Ki·단백결합률·뇌투과비·PK는 모두 아래
> 문헌에서 가져왔고, 적합하지 않았습니다. 따라서 독시라민·가바펜틴·
> 클로니딘·프로메타진의 효과 크기는 **예측**입니다.

| # | 문헌 | PMID |
|---|------|------|
| 44 | Elkomy MH, et al. **Ondansetron pharmacokinetics in pregnant women and neonates.** *Clin Pharmacol Ther* 2015. | [25670522](https://pubmed.ncbi.nlm.nih.gov/25670522/) |
| 45 | Siu SS, et al. **Placental transfer of ondansetron during early human pregnancy.** *Clin Pharmacokinet* 2006. | [16584287](https://pubmed.ncbi.nlm.nih.gov/16584287/) |
| 46 | Lemon LS, et al. **Ondansetron Exposure Changes in a Pregnant Woman.** *Pharmacotherapy* 2016. | [27374186](https://pubmed.ncbi.nlm.nih.gov/27374186/) |
| 47 | Davis R, Whittington R. **Mirtazapine: a review of its pharmacology and therapeutic potential.** *CNS Drugs* 1996. | [26071050](https://pubmed.ncbi.nlm.nih.gov/26071050/) |
| 48 | Andrews PL, et al. **The abdominal visceral innervation and the emetic reflex: pathways, pharmacology, and plasticity.** *Can J Physiol Pharmacol* 1990. | [2178756](https://pubmed.ncbi.nlm.nih.gov/2178756/) |
| 49 | Rudd JA, et al. **The involvement of TRPV1 in emesis and anti-emesis.** *Temperature* 2015. | [27227028](https://pubmed.ncbi.nlm.nih.gov/27227028/) |
| 50 | Wang XF, et al. **A meta-analysis of olanzapine for the prevention of chemotherapy-induced nausea and vomiting.** *Sci Rep* 2014. | [24770591](https://pubmed.ncbi.nlm.nih.gov/24770591/) |
| 51 | Langley-DeGroot M, et al. **Olanzapine in the treatment of refractory nausea and vomiting.** *J Pain Palliat Care Pharmacother* 2015. | [26095486](https://pubmed.ncbi.nlm.nih.gov/26095486/) |
| 52 | Grimison P, et al. **Oral Cannabis Extract for Secondary Prevention of Chemotherapy-Induced Nausea and Vomiting.** *J Clin Oncol* 2024. | [39151115](https://pubmed.ncbi.nlm.nih.gov/39151115/) |

---

## 8. 약물 안전성

| # | 문헌 | PMID |
|---|------|------|
| 53 | Huybrechts KF, et al. **Association of Maternal First-Trimester Ondansetron Use With Cardiac Malformations and Oral Clefts in Offspring.** *JAMA* 2018. | [30561479](https://pubmed.ncbi.nlm.nih.gov/30561479/) |
| 54 | Picot C, et al. **Risk of malformation after ondansetron in pregnancy: an updated systematic review and meta-analysis.** *Birth Defects Res* 2020. | [32420702](https://pubmed.ncbi.nlm.nih.gov/32420702/) |

---

## 9. 체액·전해질 — 저염소성 대사알칼리증

> 모델 3번 구조적 가정의 절반. 염소 결핍이 신장의 중탄산 배설을 막기 때문에
> 염소를 보충하지 않으면 다른 무엇을 주어도 알칼리증이 교정되지 않습니다
> (`cl_avail` 게이팅). 중탄산뇨가 나트륨 배설을 강제한다는 점이 저나트륨혈증의
> 기전입니다(`NA_HCO3_COUP`).

| # | 문헌 | PMID |
|---|------|------|
| 55 | Signorelli GC, et al. **Dietary Chloride Deficiency Syndrome: Pathophysiology, History, and Systematic Literature Review.** *Nutrients* 2020. | [33182508](https://pubmed.ncbi.nlm.nih.gov/33182508/) |
| 56 | Adrogué HJ, Madias NE. **Diagnosis and Management of Hyponatremia: a Review.** *JAMA* 2022. | [35852524](https://pubmed.ncbi.nlm.nih.gov/35852524/) |

---

## 10. 티아민과 Wernicke 뇌증 — 가장 느린 시계

> 모델 3번 구조적 가정의 나머지 절반, 그리고 임상적으로 가장 실행 가능한
> 예측: **티아민 없는 포도당 수액은 수액을 주지 않는 것보다 Wernicke 위험이
> 높다.** 모델에서 P(WE) 15.7% vs 0.0%.

| # | 문헌 | PMID |
|---|------|------|
| 57 | Oudman E, et al. **Wernicke's encephalopathy in hyperemesis gravidarum: a systematic review.** *Eur J Obstet Gynecol Reprod Biol* 2019. | [30889425](https://pubmed.ncbi.nlm.nih.gov/30889425/) |
| 58 | Erick M. **Gestational malnutrition, hyperemesis gravidarum, and Wernicke's encephalopathy: what is missing?** *Nutr Clin Pract* 2022. | [36250744](https://pubmed.ncbi.nlm.nih.gov/36250744/) |
| 59 | Maslin K, Dean C. **Nutritional consequences and management of hyperemesis gravidarum: a narrative review.** *Nutr Res Rev* 2022. | [34526158](https://pubmed.ncbi.nlm.nih.gov/34526158/) |
| 60 | Stokke G, et al. **Hyperemesis gravidarum, nutritional treatment by nasogastric tube feeding: a 10-year retrospective cohort study.** *Acta Obstet Gynecol Scand* 2015. | [25581215](https://pubmed.ncbi.nlm.nih.gov/25581215/) |
| 61 | Hsu JJ, et al. **Nasogastric enteral feeding in the management of hyperemesis gravidarum.** *Obstet Gynecol* 1996. | [8752236](https://pubmed.ncbi.nlm.nih.gov/8752236/) |

---

## 11. 내분비 — hCG의 TSH 수용체 교차자극

> hCG와 GDF15가 같은 세포(융합영양막)에서 나오므로, 모델은 생화학적
> 갑상선중독증이 **태아 생산 축**(`TROPH_GAIN`)을 따라가고 산모 감수성
> 축(`SENS`)은 따라가지 않는다고 예측합니다. PUQE가 같은 두 여성의 TSH가
> 완전히 다를 수 있다는 뜻이고, 짝지은 GDF15/TSH 측정으로 검증 가능합니다.

| # | 문헌 | PMID |
|---|------|------|
| 62 | Rodien P, et al. **Abnormal stimulation of the thyrotrophin receptor during gestation.** *Hum Reprod Update* 2004. | [15073140](https://pubmed.ncbi.nlm.nih.gov/15073140/) |
| 63 | Albaar MT, Adam JMF. **Gestational transient thyrotoxicosis.** *Acta Med Indones* 2009. | [19390130](https://pubmed.ncbi.nlm.nih.gov/19390130/) |
| 64 | Iijima S. **Pitfalls in the assessment of gestational transient thyrotoxicosis.** *Gynecol Endocrinol* 2020. | [32301638](https://pubmed.ncbi.nlm.nih.gov/32301638/) |
| 65 | Zimmerman CF, et al. **Thyroid Storm Caused by Hyperemesis Gravidarum.** *AACE Clin Case Rep* 2022. | [35602873](https://pubmed.ncbi.nlm.nih.gov/35602873/) |

---

## 12. 위장관 인자와 감별진단

| # | 문헌 | PMID |
|---|------|------|
| 66 | Ng QX, et al. **A meta-analysis of the association between Helicobacter pylori infection and hyperemesis gravidarum.** *Helicobacter* 2018. | [29178407](https://pubmed.ncbi.nlm.nih.gov/29178407/) |
| 67 | Sorensen CJ, et al. **Cannabinoid Hyperemesis Syndrome: Diagnosis, Pathophysiology, and Treatment — a Systematic Review.** *J Med Toxicol* 2017. | [28000146](https://pubmed.ncbi.nlm.nih.gov/28000146/) |

---

## 13. 산과·태아·장기 소아 결과

| # | 문헌 | PMID |
|---|------|------|
| 68 | Vandraas KF, et al. **Hyperemesis gravidarum and birth outcomes — a population-based cohort study of 2.2 million births in the Norwegian Birth Registry.** *BJOG* 2013. | [24021026](https://pubmed.ncbi.nlm.nih.gov/24021026/) |
| 69 | Nijsten K, et al. **Long-term health outcomes of children born to mothers with hyperemesis gravidarum: a systematic review and meta-analysis.** *Am J Obstet Gynecol* 2022. | [35367190](https://pubmed.ncbi.nlm.nih.gov/35367190/) |
| 70 | Pont S, et al. **Long-term health, neurodevelopmental, and educational outcomes of children born to mothers with hyperemesis gravidarum: a population-based sibling-design study.** *Am J Obstet Gynecol* 2025. | [40064411](https://pubmed.ncbi.nlm.nih.gov/40064411/) |
| 71 | Auger N, et al. **Hyperemesis gravidarum and the risk of offspring morbidity: a longitudinal cohort study.** *Eur J Pediatr* 2024. | [38884821](https://pubmed.ncbi.nlm.nih.gov/38884821/) |
| 72 | Koren G, et al. **The protective effects of nausea and vomiting of pregnancy against adverse fetal outcome — a systematic review.** *Reprod Toxicol* 2014. | [24893173](https://pubmed.ncbi.nlm.nih.gov/24893173/) |
| 73 | Boskovic R, et al. **Is lack of morning sickness teratogenic? A prospective controlled study.** *Birth Defects Res A* 2004. | [15329830](https://pubmed.ncbi.nlm.nih.gov/15329830/) |

**모델과의 긴장.** 72번·73번은 NVP가 **있는** 임신의 결과가 더 좋다는
관찰입니다. GDF15 신호를 없애는 치료가 그 자체로 무해하다고 가정할 수 없다는
뜻이고, 모델은 이 긴장을 해소하지 않습니다 — 다만 클러스터 16에 명시해
두었습니다.

---

## 14. 기전 기반 치료제 — GDF15/GFRAL 표적 (임신에서는 전부 미검증)

> 모델의 가장 큰 예측치(anti-GDF15 항체)와 가장 큰 미지수가 여기 있습니다.
> 종양 악성 소모증후군에서 ponsegromab의 인체 데이터는 리간드 차단이
> 안전하고 효과적일 수 있음을 보여주지만, **임신 1분기 사용 근거는 전무**하며
> 13번·14번 문헌은 정반대 방향의 우려를 제기합니다.

| # | 문헌 | PMID |
|---|------|------|
| 74 | Groarke JD, et al. **Ponsegromab for the Treatment of Cancer Cachexia.** *N Engl J Med* 2024. | [39282907](https://pubmed.ncbi.nlm.nih.gov/39282907/) |
| 75 | Groarke JD, et al. **Phase 2 study of the efficacy and safety of ponsegromab in patients with cancer cachexia: PROACC-1 study design.** *J Cachexia Sarcopenia Muscle* 2024. | [38500292](https://pubmed.ncbi.nlm.nih.gov/38500292/) |
| 76 | Crawford J, et al. **A Phase Ib First-In-Patient Study Assessing the Safety, Tolerability, Pharmacokinetics, and Pharmacodynamics of Ponsegromab in Participants with Cancer Cachexia.** *Clin Cancer Res* 2024. | [37982848](https://pubmed.ncbi.nlm.nih.gov/37982848/) |
| 77 | Lee BY, et al. **GDNF family receptor alpha-like antagonist antibody alleviates chemotherapy-induced cachexia in melanoma-bearing mice.** *J Cachexia Sarcopenia Muscle* 2023. | [37017344](https://pubmed.ncbi.nlm.nih.gov/37017344/) |

---

## 요약 — 어떤 문헌이 어떤 파라미터를 결정했는가

| 파라미터 | 적합/고정 | 근거 |
|----------|-----------|------|
| `ALPHA`, `TAU_SP` | **적합** (자연사) | #9, #18 — 정점 9-11주, 16주 완화, GDF15는 만삭까지 상승 |
| `W_VAG` | **적합** (온단세트론) | #35 — ΔPUQE −0.51 |
| `E0` | **적합** (미르타자핀) | #35 — ΔPUQE −1.86 |
| `R_H1`,`R_A2`,`R_A2D`,… | 고정 (약리) | #47, #48 — 수용체 부류별 Ki와 NTS 전달 관문 |
| 각 약물 PK/Ki/Kp | 고정 (문헌) | #44-#47 |
| `MET_EMAX_GDF` | 고정 | #7 — 메트포르민 → GDF15 약 2배 |
| `TOBACCO_F` | 고정 | #15 — 흡연은 메트포르민보다 작은 상승 |
| `THAL_FOLD` | 고정 | #16, #17 |
| `THI0`, `KEL_THI`, `THI_PER_CHO` | 고정 (생리) | #57, #58 |
| `GJ_CL`, `GJ_H`, `KREN_HCO3`, `NA_HCO3_COUP` | 고정 (생리) | #55, #56 |
| `AH_HCG` | 고정 | #62, #63 |
| PUQE 구간 경계 | 고정 (도구 정의) | #20, #21 |

**예측(적합하지 않은 결과)**: β-thalassemia 보호, 메트포르민 임신 전 보호와
임신 중 무효, 흡연의 중간 보호, 독시라민 효과 크기, 가바펜틴 우월성,
클로니딘 효과 크기, 코르티코스테로이드 음성 결과, 태아 생산 축에서만 나타나는
갑상선중독증, 티아민 없는 포도당의 해악.

**명시된 실패**: 메토클로프라미드를 사실상 무효로 예측하는데, 임상시험에서는
프로메타진과 비슷합니다. 모델을 반증하기 가장 쉬운 지점입니다.

---

## 도구 문헌

- mrgsolve: <https://mrgsolve.org/>
- gPKPDviz (mrgsolve 기반 Shiny): <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · <https://github.com/Genentech/gPKPDviz/>
- Graphviz: <https://graphviz.org/>
