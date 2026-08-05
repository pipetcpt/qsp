# 독사 교상(뱀독 중독) + 항독소 QSP 모델 — 참고문헌
# Snakebite Envenoming and Antivenom QSP Model — References

이 목록은 `sbe_qsp_model.dot`(기계론적 지도 · 145 노드 · 21 클러스터),
`sbe_mrgsolve_model.R`(50-ODE 모델), `sbe_reference_model.py`(독립 검증 구현)에
들어간 구조적 선택과 파라미터의 근거입니다.

**모든 항목은 PubMed E-utilities(esearch + esummary)로 직접 조회한 실제
레코드입니다** — PMID·저자·연도·저널·제목은 조회 결과를 그대로 옮겼고, 기억에
의존해 작성한 인용은 하나도 없습니다.

Every entry below is a real record retrieved directly from PubMed via the NCBI
E-utilities API. PMIDs, first authors, years, journals and titles are
transcribed from that retrieval rather than written from memory.

---

## 이 모델이 문헌에서 가져온 것과, 가져오지 않은 것
## What is sourced, and what is not

정직하게 구분하는 것이 이 파일의 절반입니다. QSP 모델에서 가장 위험한 것은
보정 파라미터가 관측값처럼 읽히는 것입니다.

**문헌에서 직접 가져온 것 (sourced from the literature)**

- 항독소 절편별 종말 반감기 — ovine Fab 12–23 h, equine F(ab')₂ ≈ 130 h,
  equine whole IgG 90–200 h. 모델의 2-구획 파라미터는 이 값을 재현하도록
  선택되었고, `report_product_pk()`가 실제로 계산해 확인합니다.
- 인도 다가 항독소(ASV)의 **표시 역가** — 1 mL당 *Daboia russelii* 독 0.6 mg
  중화, 즉 10 mL 바이알당 6 mgNE. 모델의 항독소 단위(mgNE)는 이 규제 수치
  자체이며 적합된 값이 아닙니다.
- 독액 계열별 질량 분율 (venomics/proteomics: Echis 68% SVMP, Naja 55% 3FTx,
  Bungarus 66% PLA2 등).
- 20분 전혈응고검사(20WBCT)가 피브리노겐 약 0.5 g/L에서 양성으로 전환된다는
  진단 임계값.
- 피브리노겐 정상치 2.8 g/L, 반감기 96 h, 급성기 합성 상향 조절 범위.
- α-신경독소의 준가역적 해리(수십 시간)와 β-신경독소(PLA2)의 **비가역적**
  운동종말 파괴라는 구분, 그리고 항콜린에스테라제가 후자에 효과가 없다는
  임상 사실.
- 건성 교상(dry bite) 빈도, 항독소 급성 이상반응 및 혈청병 빈도의 제품별 차이,
  헤파린이 VICC에 효과가 없다는 무작위 시험 결과.
- 러셀살무사 교상의 임상 시간 경과 — VICC 발생·회복 시점, 신경독성 징후,
  급성 신손상의 3–5일 정점.

**모델 수준의 가설이며 관측이 아닌 것 (model-level hypotheses, NOT observations)**

- **ε (계열별 중화 효율) 벡터 = SVMP 1.20 / SVSP 1.00 / PLA2 0.60 / 3FTx 0.35.**
  항독소가 큰 효소를 잘 덮고 작고 면역원성이 낮은 독소를 못 덮는다는 것은
  antivenomics 문헌의 정성적 결론이지만, **이 네 개의 숫자는 정량적으로 보고된
  값이 아니고 모델 파라미터입니다.** 코브라 교상에 20–30 바이알이 필요하다는
  임상 관행을 산수로 재현하도록 선택되었습니다.
- **느린 저장고 (f_slow, ka_s).** 독액 항원이 수일간 검출된다는 관측은 실재하지만
  "주입량의 35%가 깊은 구획에 격리되어 t½ 45 h로 방출된다"는 구획화는 모델의
  구조적 가설입니다. §E의 민감도 분석이 이 가설의 한계를 명시적으로 검정합니다.
- 교상부 국소 중화 벌점 (k_b,loc/k_b0 = 0.038)의 정확한 크기.
- 신경근 안전계수 SF₀ = 4.0, 안검하수 임계값 SF < 2.2, 인공환기 임계값
  VC_frac < 0.15 — 순서와 방향은 견고하지만 절대값은 보정값입니다.
- 신장 손실 항의 모든 속도 상수 (k_n,fib · k_n,cast · k_n,isch · k_n,dir),
  그리고 손실의 10%가 영구 반흔이 된다는 가정.
- 사망·출혈 위험함수의 기저 위험도.
- 독액-항독소 복합체의 해리 속도 k_off = 0 (기본값). §E에서 이 값을 sweep 하여
  임상적 재발성 응고장애를 설명하려면 얼마여야 하는지 계산합니다.

보고된 **비율·순서·시간 창**은 절대값보다 훨씬 견고합니다. 이 모델의 네 가지
주장은 모두 비율에 관한 것입니다.

---

## 1. 질병부담 · 역학 · 세계 보건 전략
## Burden, epidemiology and the global health strategy

Why this is a pharmacology problem and not only a first-aid problem: the deaths
are concentrated where the drug supply and the ventilator are not.

- GBD 2019 Snakebite Envenomation Collaborators. **Global mortality of snakebite envenoming between 1990 and 2019.** *Nat Commun*. 2022. [PMID 36284094](https://pubmed.ncbi.nlm.nih.gov/36284094/)
- Afroz A, et al. **Snakebite envenoming: A systematic review and meta-analysis of global morbidity and mortality.** *PLoS Negl Trop Dis*. 2024. [PMID 38574167](https://pubmed.ncbi.nlm.nih.gov/38574167/)
- Longbottom J, et al. **Vulnerability to snakebite envenoming: a global mapping of hotspots.** *Lancet*. 2018. [PMID 30017551](https://pubmed.ncbi.nlm.nih.gov/30017551/)
- Iliyasu G, et al. **Case fatality rate and burden of snakebite envenoming in children — A systematic review and meta-analysis.** *Toxicon*. 2023. [PMID 37739273](https://pubmed.ncbi.nlm.nih.gov/37739273/)
- Gololo AA, et al. **Epidemiological models to estimate the burden of snakebite envenoming: A systematic review.** *Trop Med Int Health*. 2025. [PMID 39743841](https://pubmed.ncbi.nlm.nih.gov/39743841/)
- Menon JC, et al. **ICMR task force project — survey of the incidence, mortality, morbidity and socio-economic burden of snakebite in India: A study protocol.** *PLoS One*. 2022. [PMID 35994445](https://pubmed.ncbi.nlm.nih.gov/35994445/)
- Alcoba G, et al. **Snakebite epidemiology and health-seeking behavior in Akonolinga health district, Cameroon: Cross-sectional study.** *PLoS Negl Trop Dis*. 2020. [PMID 32584806](https://pubmed.ncbi.nlm.nih.gov/32584806/)
- Pach S, et al. **Paediatric snakebite envenoming: the world's most neglected 'Neglected Tropical Disease'?** *Arch Dis Child*. 2020. [PMID 32998874](https://pubmed.ncbi.nlm.nih.gov/32998874/)
- Minghui R, et al. **WHO's Snakebite Envenoming Strategy for prevention and control.** *Lancet Glob Health*. 2019. [PMID 31129124](https://pubmed.ncbi.nlm.nih.gov/31129124/)
- Williams DJ, et al. **Strategy for a globally coordinated response to a priority neglected tropical disease: Snakebite envenoming.** *PLoS Negl Trop Dis*. 2019. [PMID 30789906](https://pubmed.ncbi.nlm.nih.gov/30789906/)
- Chippaux JP. **The WHO strategy for prevention and control of snakebite envenoming: a sub-Saharan Africa plan.** *J Venom Anim Toxins Incl Trop Dis*. 2019. [PMID 31839803](https://pubmed.ncbi.nlm.nih.gov/31839803/)
- Harrison RA, et al. **The time is now: a call for action to translate recent momentum on tackling tropical snakebite into sustained benefit for victims.** *Trans R Soc Trop Med Hyg*. 2019. [PMID 30668842](https://pubmed.ncbi.nlm.nih.gov/30668842/)

## 2. 통합 리뷰 (모델의 전체 골격)
## Integrative reviews — the skeleton the map is built on

- Gutiérrez JM, Calvete JJ, Habib AG, Harrison RA, Williams DJ, Warrell DA. **Snakebite envenoming.** *Nat Rev Dis Primers*. 2017. [PMID 28905944](https://pubmed.ncbi.nlm.nih.gov/28905944/)
- Gutiérrez JM, et al. **Understanding and confronting snakebite envenoming: The harvest of cooperation.** *Toxicon*. 2016. [PMID 26615826](https://pubmed.ncbi.nlm.nih.gov/26615826/)
- Clare RH, et al. **Small Molecule Drug Discovery for Neglected Tropical Snakebite.** *Trends Pharmacol Sci*. 2021. [PMID 33773806](https://pubmed.ncbi.nlm.nih.gov/33773806/)

## 3. 독액 단백질체 — 4개 계열로 나눈 근거
## The venom proteome — the basis for four toxin classes

The mass fractions in `SNAKES` come from this literature. They are the reason
*Echis* is a pure coagulopathy, *Naja* a pure paralysis, and *Daboia* both.

- Calvete JJ, et al. **Exploring the venom proteome of the western diamondback rattlesnake, Crotalus atrox, via snake venomics and combinatorial peptide ligand library approaches.** *J Proteome Res*. 2009. [PMID 19371136](https://pubmed.ncbi.nlm.nih.gov/19371136/)
- Serrano SM, et al. **A multifaceted analysis of viperid snake venoms by two-dimensional gel electrophoresis: an approach to understanding venom proteomics.** *Proteomics*. 2005. [PMID 15627971](https://pubmed.ncbi.nlm.nih.gov/15627971/)
- Slagboom J, et al. **High-Throughput Venomics.** *J Proteome Res*. 2023. [PMID 37010854](https://pubmed.ncbi.nlm.nih.gov/37010854/)
- Oh AMF, et al. **Venomics of Bungarus caeruleus (Indian krait): Comparable venom profiles, variable immunoreactivities among specimens from Sri Lanka, India and Pakistan.** *J Proteomics*. 2017. [PMID 28476572](https://pubmed.ncbi.nlm.nih.gov/28476572/)
- Hia YL, et al. **Comparative venom proteomics of banded krait (Bungarus fasciatus) from five geographical locales.** *Acta Trop*. 2020. [PMID 32278639](https://pubmed.ncbi.nlm.nih.gov/32278639/)
- Bittenbinder MA, et al. **Coagulotoxic Cobras: Clinical Implications of Strong Anticoagulant Actions of African Spitting Naja Venoms That Are Not Neutralised by Antivenom but Are by LY315920 (Varespladib).** *Toxins (Basel)*. 2018. [PMID 30518149](https://pubmed.ncbi.nlm.nih.gov/30518149/)

## 4. SVMP — 금속단백분해효소: 출혈독소 · 프로트롬빈/FX 활성화 · 기저막 소화
## Snake venom metalloproteinases

- Kamiguti AS, Hay CR, Theakston RD, Zuzel M. **Insights into the mechanism of haemorrhage caused by snake venom metalloproteinases.** *Toxicon*. 1996. [PMID 8817809](https://pubmed.ncbi.nlm.nih.gov/8817809/)
- Moura-da-Silva AM, et al. **Jararhagin, a hemorrhagic snake venom metalloproteinase from Bothrops jararaca.** *Toxicon*. 2012. [PMID 22534074](https://pubmed.ncbi.nlm.nih.gov/22534074/)
- Valente RH, et al. **BJ46a, a snake venom metalloproteinase inhibitor. Isolation, characterization, cloning and insights into its mechanism of action.** *Eur J Biochem*. 2001. [PMID 11358523](https://pubmed.ncbi.nlm.nih.gov/11358523/)
- Bickler PE. **Amplification of Snake Venom Toxicity by Endogenous Signaling Pathways.** *Toxins (Basel)*. 2020. [PMID 31979014](https://pubmed.ncbi.nlm.nih.gov/31979014/)

## 5. SVSP — 세린단백분해효소: 트롬빈 유사 효소와 피브리노겐 소모
## Snake venom serine proteases and thrombin-like enzymes

The `kcat_svsp_fg` term and the Michaelis–Menten form of fibrinogen consumption
come from this pharmacology: these enzymes cleave fibrinopeptide A only, so the
fibrin they make is friable and the fibrinogen they consume is simply gone.

- Stocker K, Barlow GH. **Thrombin-like snake venom proteinases.** *Toxicon*. 1982. [PMID 7043783](https://pubmed.ncbi.nlm.nih.gov/7043783/)
- Aronson DL. **Comparison of the actions of thrombin and the thrombin-like venom enzymes ancrod and batroxobin.** *Thromb Haemost*. 1976. [PMID 1036831](https://pubmed.ncbi.nlm.nih.gov/1036831/)
- Hahn BS, et al. **Purification and molecular cloning of calobin, a thrombin-like enzyme from Agkistrodon caliginosus (Korean viper).** *J Biochem*. 1996. [PMID 8797081](https://pubmed.ncbi.nlm.nih.gov/8797081/)
- Cho SY, et al. **Purification and characterization of calobin II, a second type of thrombin-like enzyme from Agkistrodon caliginosus (Korean viper).** *Toxicon*. 2001. [PMID 11024490](https://pubmed.ncbi.nlm.nih.gov/11024490/)

## 6. svPLA2 — 근독성 · 전시냅스 신경독성 · 부종
## Secreted phospholipase A₂: myotoxicity, presynaptic neurotoxicity, oedema

- Gutiérrez JM, Rucavado A, Chaves F, Díaz C, Escalante T. **Experimental pathology of local tissue damage induced by Bothrops asper snake venom.** *Toxicon*. 2009. [PMID 19303033](https://pubmed.ncbi.nlm.nih.gov/19303033/)
- Gutiérrez JM, León G, Rojas G, Lomonte B, Rucavado A, Chaves F. **Neutralization of local tissue damage induced by Bothrops asper (terciopelo) snake venom.** *Toxicon*. 1998. [PMID 9792169](https://pubmed.ncbi.nlm.nih.gov/9792169/)
- Chaves F, et al. **Histopathological and biochemical alterations induced by intramuscular injection of Bothrops asper (terciopelo) venom in mice.** *Toxicon*. 1989. [PMID 2815106](https://pubmed.ncbi.nlm.nih.gov/2815106/)
- Rucavado A, et al. **Analysis of wound exudates reveals differences in the patterns of tissue damage and inflammation induced by the venoms of Daboia russelii and Bothrops asper.** *Toxicon*. 2020. [PMID 32781076](https://pubmed.ncbi.nlm.nih.gov/32781076/)
- Rucavado A, et al. **Increments in cytokines and matrix metalloproteinases in skeletal muscle after injection of tissue-damaging toxins from the venom of the snake Bothrops asper.** *Mediators Inflamm*. 2002. [PMID 12061424](https://pubmed.ncbi.nlm.nih.gov/12061424/)
- Silva A, et al. **Clinical and Pharmacological Investigation of Myotoxicity in Sri Lankan Russell's Viper (Daboia russelii) Envenoming.** *PLoS Negl Trop Dis*. 2016. [PMID 27911900](https://pubmed.ncbi.nlm.nih.gov/27911900/)
- Sanhajariya S, et al. **Investigating myotoxicity following Australian red-bellied black snake (Pseudechis porphyriacus) envenomation.** *PLoS One*. 2021. [PMID 34506531](https://pubmed.ncbi.nlm.nih.gov/34506531/)
- Silva MDS, et al. **NLRP3 inflammasome activation in human peripheral blood mononuclear cells induced by venoms secreted PLA₂s.** *Int J Biol Macromol*. 2022. [PMID 35074331](https://pubmed.ncbi.nlm.nih.gov/35074331/)

## 7. 3FTx — 삼지 독소: 후시냅스 α-신경독소와 세포독소
## Three-finger toxins: postsynaptic α-neurotoxins and cytotoxins

The apparent perijunctional K_d and the slow off-rate (t½ ≈ 35 h, and ~230 h for
α-bungarotoxin) come from this literature, and they are why an occupancy decays
on its own clock rather than being stripped off by antivenom.

- Nirthanan S. **Snake three-finger α-neurotoxins and nicotinic acetylcholine receptors: molecules, mechanisms and medicine.** *Biochem Pharmacol*. 2020. [PMID 32710970](https://pubmed.ncbi.nlm.nih.gov/32710970/)
- Nirthanan S, Gwee MC. **Three-finger alpha-neurotoxins and the nicotinic acetylcholine receptor, forty years on.** *J Pharmacol Sci*. 2004. [PMID 14745112](https://pubmed.ncbi.nlm.nih.gov/14745112/)
- Shenkarev ZO, et al. **Membrane-mediated interaction of non-conventional snake three-finger toxins with nicotinic acetylcholine receptors.** *Commun Biol*. 2022. [PMID 36477694](https://pubmed.ncbi.nlm.nih.gov/36477694/)
- Bekbossynova A, et al. **Venom-Derived Neurotoxins Targeting Nicotinic Acetylcholine Receptors.** *Molecules*. 2021. [PMID 34204855](https://pubmed.ncbi.nlm.nih.gov/34204855/)
- Choudhury M, et al. **Snake Venom Three-Finger Neurotoxins and Neurotoxin-Like Proteins: Insights Into Their Structural and Functional Aspects.** *Basic Clin Pharmacol Toxicol*. 2025. [PMID 41230907](https://pubmed.ncbi.nlm.nih.gov/41230907/)
- Tamiya N, Sato A. **Studies on sea snake venom.** *Proc Jpn Acad Ser B Phys Biol Sci*. 2011. [PMID 21422738](https://pubmed.ncbi.nlm.nih.gov/21422738/)

## 8. β-신경독소 — 전시냅스 파괴: 항독소가 되돌릴 수 없는 손상
## β-neurotoxins: the presynaptic destruction antivenom cannot undo

This section is the evidentiary basis of the model's claim 3. The `dTERM/dt`
equation contains a destruction term and a 5-day regeneration term, and neither
antivenom nor neostigmine appears in it.

- Galappaththige J, et al. **Neurotoxicity of Sri Lankan Krait (Bungarus ceylonicus) and Common Krait (Bungarus caeruleus) Venoms and Their Neutralisation by Commercial Antivenoms.** *Toxins (Basel)*. 2025. [PMID 41003503](https://pubmed.ncbi.nlm.nih.gov/41003503/)
- Liang Q, et al. **In Vitro Neurotoxicity of Chinese Krait (Bungarus multicinctus) Venom and Neutralization by Antivenoms.** *Toxins (Basel)*. 2021. [PMID 33440641](https://pubmed.ncbi.nlm.nih.gov/33440641/)
- Herkert M, et al. **Beta-bungarotoxin is a potent inducer of apoptosis in cultured rat neurons by receptor-mediated internalization.** *Eur J Neurosci*. 2001. [PMID 11576186](https://pubmed.ncbi.nlm.nih.gov/11576186/)
- Othman IB, Spokes JW, Dolly JO. **Preparation of neurotoxic ³H-beta-bungarotoxin: demonstration of saturable binding to brain synapses and its inhibition by toxin I.** *Eur J Biochem*. 1982. [PMID 7173209](https://pubmed.ncbi.nlm.nih.gov/7173209/)

## 9. 신경독성 교상의 임상 — 크레이트와 코브라
## Clinical neurotoxic envenoming — krait and cobra

The ventilation durations (4–10 days), the delayed krait onset, and the
observation that anticholinesterases help the cobra and not the krait, are all
read off these series.

- Kularatne SA. **Common krait (Bungarus caeruleus) bite in Anuradhapura, Sri Lanka: a prospective clinical study, 1996-98.** *Postgrad Med J*. 2002. [PMID 12151569](https://pubmed.ncbi.nlm.nih.gov/12151569/)
- Bawaskar HS, Bawaskar PH. **Envenoming by the common krait (Bungarus caeruleus) and Asian cobra (Naja naja): clinical manifestations and their management in a rural setting.** *Wilderness Environ Med*. 2004. [PMID 15636376](https://pubmed.ncbi.nlm.nih.gov/15636376/)
- Theakston RD, et al. **Envenoming by the common krait (Bungarus caeruleus) and Sri Lankan cobra (Naja naja naja): efficacy and complications of therapy with Haffkine antivenom.** *Trans R Soc Trop Med Hyg*. 1990. [PMID 2389328](https://pubmed.ncbi.nlm.nih.gov/2389328/)
- Pannu AK, et al. **Efficacy of a low dose of antivenom for severe neuroparalysis in Bungarus caeruleus (common krait) envenomation: a pilot study.** *Toxicol Res (Camb)*. 2024. [PMID 38450179](https://pubmed.ncbi.nlm.nih.gov/38450179/)
- Gupta A, et al. **Unusually prolonged neuromuscular weakness caused by krait (Bungarus caeruleus) bite: Two case reports.** *Toxicon*. 2021. [PMID 33497743](https://pubmed.ncbi.nlm.nih.gov/33497743/)
- Faiz MA, et al. **Bites by the Monocled Cobra, Naja kaouthia, in Chittagong Division, Bangladesh: Epidemiology, Clinical Features of Envenoming and Management of 70 Identified Cases.** *Am J Trop Med Hyg*. 2017. [PMID 28138054](https://pubmed.ncbi.nlm.nih.gov/28138054/)
- Sarin K, et al. **Clinical profile & complications of neurotoxic snake bite & comparison of two regimens of polyvalent anti-snake venom in its treatment.** *Indian J Med Res*. 2017. [PMID 28574015](https://pubmed.ncbi.nlm.nih.gov/28574015/)
- Silva A, et al. **Neurotoxicity in Russell's viper (Daboia russelii) envenoming in Sri Lanka: a clinical and neurophysiological study.** *Clin Toxicol (Phila)*. 2016. [PMID 26923566](https://pubmed.ncbi.nlm.nih.gov/26923566/)

## 10. 독소 유발 소모성 응고장애 (VICC)
## Venom-induced consumption coagulopathy

- Maduwage K, Isbister GK. **Current treatment for venom-induced consumption coagulopathy resulting from snakebite.** *PLoS Negl Trop Dis*. 2014. [PMID 25340841](https://pubmed.ncbi.nlm.nih.gov/25340841/)
- Maduwage K, Buckley NA, de Silva HJ, Lalloo DG, Isbister GK. **Snake antivenom for snake venom induced consumption coagulopathy.** *Cochrane Database Syst Rev*. 2015. [PMID 26058967](https://pubmed.ncbi.nlm.nih.gov/26058967/)
- Wedasingha S, et al. **Bedside Coagulation Tests in Diagnosing Venom-Induced Consumption Coagulopathy in Snakebite.** *Toxins (Basel)*. 2020. [PMID 32927702](https://pubmed.ncbi.nlm.nih.gov/32927702/)
- Valenta J, et al. **Fibrinogenolysis in Venom-Induced Consumption Coagulopathy after Viperidae Snakebites: A Pilot Study.** *Toxins (Basel)*. 2022. [PMID 36006200](https://pubmed.ncbi.nlm.nih.gov/36006200/)
- Little M, Isbister GK. **Is D-dimer the new test for venom-induced consumption coagulopathy after snakebite?** *Med J Aust*. 2022. [PMID 35843724](https://pubmed.ncbi.nlm.nih.gov/35843724/)
- Rajkumar B, et al. **Venom induced consumption coagulopathy and performance of 20-min whole blood clotting test for its detection in viperid envenomation.** *J R Coll Physicians Edinb*. 2022. [PMID 36300884](https://pubmed.ncbi.nlm.nih.gov/36300884/)
- Silva A, et al. **Indian Polyvalent Antivenom Accelerates Recovery From Venom-Induced Consumption Coagulopathy (VICC) in Sri Lankan Russell's Viper (Daboia russelii) Envenoming.** *Front Med (Lausanne)*. 2022. [PMID 35321467](https://pubmed.ncbi.nlm.nih.gov/35321467/)

**이 절이 claim 1의 근거입니다.** 특히 Brown 2009는 응고인자 보충이 회복
'속도'가 아니라 '깊이'를 바꾼다는 본 모델의 구조를 임상에서 직접 다룹니다.

- Brown SG, Caruso N, Borland ML, McCoubrie DL, Celenza A, Isbister GK. **Clotting factor replacement and recovery from snake venom-induced consumptive coagulopathy.** *Intensive Care Med*. 2009. [PMID 19547954](https://pubmed.ncbi.nlm.nih.gov/19547954/)

## 11. 20분 전혈응고검사 — 비용이 들지 않는 관문
## The 20-minute whole blood clotting test

The model's `FG_unclot = 0.5 g/L` threshold and the "dry bite gate" argument in
D6 rest on this diagnostic literature.

- Lamb T, et al. **The 20-minute whole blood clotting test (20WBCT) for snakebite coagulopathy — A systematic review and meta-analysis of diagnostic test accuracy.** *PLoS Negl Trop Dis*. 2021. [PMID 34375338](https://pubmed.ncbi.nlm.nih.gov/34375338/)
- Lamb T, et al. **Correction: The 20-minute whole blood clotting test (20WBCT) for snakebite coagulopathy.** *PLoS Negl Trop Dis*. 2023. [PMID 36693087](https://pubmed.ncbi.nlm.nih.gov/36693087/)
- Thongtonyong N, Chinthammitr Y. **Sensitivity and specificity of 20-minute whole blood clotting test, prothrombin time, activated partial thromboplastin time tests in diagnosis of defibrination syndrome.** *Toxicon*. 2020. [PMID 32712023](https://pubmed.ncbi.nlm.nih.gov/32712023/)
- Tianyi FL, et al. **Diagnostic characteristics of the 20-minute whole blood clotting test in detecting venom-induced consumptive coagulopathy following carpet viper envenoming.** *PLoS Negl Trop Dis*. 2023. [PMID 37363905](https://pubmed.ncbi.nlm.nih.gov/37363905/)
- Hamza M, et al. **Performance of the 20 minutes Whole Blood Clotting Test in detection, monitoring and antivenom therapy of West African Carpet viper (Echis romani) envenoming.** *Toxicon*. 2023. [PMID 36640811](https://pubmed.ncbi.nlm.nih.gov/36640811/)
- Suseel A, et al. **Comparing modified Lee and White method against 20-minute whole blood clotting test as bedside coagulation screening test in snake envenomation victims.** *J Venom Anim Toxins Incl Trop Dis*. 2023. [PMID 37342654](https://pubmed.ncbi.nlm.nih.gov/37342654/)

## 12. 혈소판 감소 — C형 렉틴과 GPIb
## Venom-induced thrombocytopenia

- Lu Q, Clemetson JM, Clemetson KJ. **Snake venom C-type lectins interacting with platelet receptors. Structure-function relationships and effects on haemostasis.** *Toxicon*. 2005. [PMID 15876445](https://pubmed.ncbi.nlm.nih.gov/15876445/)
- Long C, et al. **Potential Role of Platelet-Activating C-Type Lectin-Like Proteins in Viper Envenomation Induced Thrombotic Microangiopathy Symptom.** *Toxins (Basel)*. 2020. [PMID 33260875](https://pubmed.ncbi.nlm.nih.gov/33260875/)
- Shen C, et al. **Viper venoms drive the macrophages and hepatocytes to sequester and clear platelets: novel mechanism and therapeutic strategy for venom-induced thrombocytopenia.** *Arch Toxicol*. 2021. [PMID 34519865](https://pubmed.ncbi.nlm.nih.gov/34519865/)
- Thomazini CM, et al. **Involvement of von Willebrand factor and botrocetin in the thrombocytopenia induced by Bothrops jararaca snake venom.** *PLoS Negl Trop Dis*. 2021. [PMID 34478462](https://pubmed.ncbi.nlm.nih.gov/34478462/)
- Xu G, et al. **How does agkicetin-C bind on platelet glycoprotein Ibalpha and achieve its platelet effects?** *Toxicon*. 2005. [PMID 15777951](https://pubmed.ncbi.nlm.nih.gov/15777951/)

## 13. 독액 및 항독소 약동학 — 이 모델의 PK 골격
## Venom and antivenom pharmacokinetics

Isbister 2015 is the anchor for the antivenom disposition parameters; Morris 2022
is the closest prior computational model of the same system and the reason this
model was written as two clocks rather than one.

- Isbister GK, et al. **Population Pharmacokinetics of an Indian F(ab')₂ Snake Antivenom in Patients with Russell's Viper (Daboia russelii) Bites.** *PLoS Negl Trop Dis*. 2015. [PMID 26135318](https://pubmed.ncbi.nlm.nih.gov/26135318/)
- Sanhajariya S, et al. **Population pharmacokinetics of Pseudechis porphyriacus (red-bellied black snake) venom in snakebite patients.** *Clin Toxicol (Phila)*. 2021. [PMID 33832399](https://pubmed.ncbi.nlm.nih.gov/33832399/)
- Morris NM, et al. **Developing a computational pharmacokinetic model of systemic snakebite envenomation and antivenom treatment.** *Toxicon*. 2022. [PMID 35716719](https://pubmed.ncbi.nlm.nih.gov/35716719/)
- Audebert F, et al. **Viper bites in France: clinical and biological evaluation; kinetics of envenomations.** *Hum Exp Toxicol*. 1994. [PMID 7826686](https://pubmed.ncbi.nlm.nih.gov/7826686/)
- Theakston RD. **Snake venoms in science and clinical medicine. 2. Applied immunology in snake venom research.** *Trans R Soc Trop Med Hyg*. 1989. [PMID 2617643](https://pubmed.ncbi.nlm.nih.gov/2617643/)
- Chavez-Olortegui C, et al. **An enzyme linked immunosorbent assay (ELISA) that discriminates between Bothrops atrox and Lachesis muta muta venoms.** *Toxicon*. 1993. [PMID 8503131](https://pubmed.ncbi.nlm.nih.gov/8503131/)

## 14. 재발 — 짧은 절편이 긴 꼬리를 놓치는 현상
## Recurrence after a short-lived fragment

**Claim 2의 임상적 근거.** 오빈 Fab 투여 후 재발성 응고장애·혈소판감소·항원혈증이
반복적으로 보고되었고, CroFab 라벨의 유지용량 스케줄(2 바이알 q6h ×3)은 그
결과입니다. §E는 이 문헌이 보고하는 *중증* 재발이 본 모델의 저장고 방출만으로는
설명되지 않는다는 음성 결과를 명시합니다.

- Ruha AM, Curry SC, Albrecht C, Riley B, Pizon A. **Late hematologic toxicity following treatment of rattlesnake envenomation with crotalidae polyvalent immune Fab antivenom.** *Toxicon*. 2011. [PMID 20920516](https://pubmed.ncbi.nlm.nih.gov/20920516/)
- Miller AD, et al. **Recurrent coagulopathy and thrombocytopenia in children treated with crotalidae polyvalent immune fab: a case series.** *Pediatr Emerg Care*. 2010. [PMID 20693856](https://pubmed.ncbi.nlm.nih.gov/20693856/)
- Fazelat J, Teperman SH, Touger M. **Recurrent hemorrhage after western diamondback rattlesnake envenomation treated with crotalidae polyvalent immune fab (ovine).** *Clin Toxicol (Phila)*. 2008. [PMID 18608290](https://pubmed.ncbi.nlm.nih.gov/18608290/)
- Lavonas EJ, et al. **Initial experience with Crotalidae polyvalent immune Fab (ovine) antivenom in the treatment of copperhead snakebite.** *Ann Emerg Med*. 2004. [PMID 14747809](https://pubmed.ncbi.nlm.nih.gov/14747809/)
- Keating GM. **Crotalidae polyvalent immune Fab: in patients with North American crotaline envenomation.** *BioDrugs*. 2011. [PMID 21443271](https://pubmed.ncbi.nlm.nih.gov/21443271/)
- Keating GM. **Crotalidae polyvalent immune Fab: a guide to its use in North American crotaline envenomation.** *Clin Drug Investig*. 2012. [PMID 22765769](https://pubmed.ncbi.nlm.nih.gov/22765769/)
- Boyer L, et al. **Safety of intravenous equine F(ab')₂: insights following clinical trials involving 1534 recipients of scorpion antivenom.** *Toxicon*. 2013. [PMID 23916602](https://pubmed.ncbi.nlm.nih.gov/23916602/)

## 15. 항독소 용량 · 역가 · 임상시험
## Antivenom dose, potency and trials

The labelled-potency unit (mgNE) and the "10 vials is a 0.98× margin" arithmetic
in section B come from this literature.

- Abubakar IS, et al. **Randomised controlled double-blind non-inferiority trial of two antivenoms for saw-scaled or carpet viper (Echis ocellatus) envenoming in Nigeria.** *PLoS Negl Trop Dis*. 2010. [PMID 20668549](https://pubmed.ncbi.nlm.nih.gov/20668549/)
- Abubakar SB, et al. **Pre-clinical and preliminary dose-finding and safety studies to identify candidate antivenoms for treatment of envenoming by saw-scaled or carpet vipers (Echis ocellatus) in Nigeria.** *Toxicon*. 2010. [PMID 19874841](https://pubmed.ncbi.nlm.nih.gov/19874841/)
- Meyer WP, et al. **First clinical experiences with a new ovine Fab Echis ocellatus snake bite antivenom in Nigeria: randomized comparative trial with Institute Pasteur Serum.** *Am J Trop Med Hyg*. 1997. [PMID 9129531](https://pubmed.ncbi.nlm.nih.gov/9129531/)
- Sagar P, et al. **Comparison of two Anti Snake Venom protocols in hemotoxic snake bite: A randomized trial.** *J Forensic Leg Med*. 2020. [PMID 32658754](https://pubmed.ncbi.nlm.nih.gov/32658754/)
- Patra A, Mukherjee AK. **Assessment of quality and pre-clinical efficacy of a newly developed polyvalent antivenom against the medically important snakes of Sri Lanka.** *Sci Rep*. 2021. [PMID 34521877](https://pubmed.ncbi.nlm.nih.gov/34521877/)
- Lim ASS, et al. **Immunoreactivity and neutralization efficacy of Pakistani Viper Antivenom (PVAV) against venoms of Saw-scaled Vipers (Echis carinatus subspp.) and Western Russell's Viper.** *Acta Trop*. 2024. [PMID 38097152](https://pubmed.ncbi.nlm.nih.gov/38097152/)
- Chippaux JP. **[Guidelines for the production, control and regulation of snake antivenom immunoglobulins].** *Biol Aujourdhui*. 2010. [PMID 20950580](https://pubmed.ncbi.nlm.nih.gov/20950580/)
- World Health Organization. **WHO Expert Committee on Biological Standardization.** *World Health Organ Tech Rep Ser*. 2012. [PMID 22900409](https://pubmed.ncbi.nlm.nih.gov/22900409/)
- Scheske L, Ruitenberg J, Bissumbhar B. **Needs and availability of snake antivenoms: relevance and application of international guidelines.** *Int J Health Policy Manag*. 2015. [PMID 26188809](https://pubmed.ncbi.nlm.nih.gov/26188809/)

## 16. 항독소 이상반응 — 장부의 비용 쪽
## Antivenom adverse reactions — the cost side of the ledger

The rate-driven anaphylactoid term (`rho * AVRATE`) and the immune-complex →
serum-sickness cascade are parameterised against these product-specific rates.

- Habib AG. **Effect of pre-medication on early adverse reactions following antivenom use in snakebite: a systematic review and meta-analysis.** *Drug Saf*. 2011. [PMID 21879781](https://pubmed.ncbi.nlm.nih.gov/21879781/)
- Morais V. **Antivenom therapy: efficacy of premedication for the prevention of adverse reactions.** *J Venom Anim Toxins Incl Trop Dis*. 2018. [PMID 29507580](https://pubmed.ncbi.nlm.nih.gov/29507580/)
- Ryan NM, Kearney RT, Brown SG, Isbister GK. **Incidence of serum sickness after the administration of Australian snake antivenom (ASP-22).** *Clin Toxicol (Phila)*. 2016. [PMID 26490786](https://pubmed.ncbi.nlm.nih.gov/26490786/)
- Tu M, et al. **Safety profile of antivenom in a cohort of patients envenomed by Deinagkistrodon acutus in Hangzhou, Zhejiang Province, Southeast China.** *Clin Toxicol (Phila)*. 2025. [PMID 40106271](https://pubmed.ncbi.nlm.nih.gov/40106271/)
- LoVecchio F, et al. **Incidence of immediate and delayed hypersensitivity to Centruroides antivenom.** *Ann Emerg Med*. 1999. [PMID 10533009](https://pubmed.ncbi.nlm.nih.gov/10533009/)
- Sheikh A. **Glucocorticosteroids for the treatment and prevention of anaphylaxis.** *Curr Opin Allergy Clin Immunol*. 2013. [PMID 23507835](https://pubmed.ncbi.nlm.nih.gov/23507835/)

## 17. 독소 유발 급성 신손상 — 적분으로서의 신장
## Venom-induced acute kidney injury

Claim 4. The day-3-to-5 creatinine peak, the fibrin-microthrombus and pigment
cast pathways, and bilateral cortical necrosis as the irreversible end of the
same integral, come from this literature.

- Sitprija V. **Snakebite nephropathy.** *Nephrology (Carlton)*. 2006. [PMID 17014559](https://pubmed.ncbi.nlm.nih.gov/17014559/)
- Chugh KS. **Snake-bite-induced acute renal failure in India.** *Kidney Int*. 1989. [PMID 2651763](https://pubmed.ncbi.nlm.nih.gov/2651763/)
- Vikrant S, Jaryal A, Parashar A. **Clinicopathological spectrum of snake bite-induced acute kidney injury from India.** *World J Nephrol*. 2017. [PMID 28540205](https://pubmed.ncbi.nlm.nih.gov/28540205/)
- Rao PSK, et al. **Snakebite envenomation-associated acute kidney injury: a South-Asian perspective.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749490](https://pubmed.ncbi.nlm.nih.gov/39749490/)
- Alvitigala BY, et al. **Snakebite-associated acute kidney injury in South Asia: narrative review on epidemiology, pathogenesis and management.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749470](https://pubmed.ncbi.nlm.nih.gov/39749470/)
- Meena P, et al. **The kidney histopathological spectrum of patients with kidney injury following snakebite envenomation in India: scoping review of five decades.** *BMC Nephrol*. 2024. [PMID 38515042](https://pubmed.ncbi.nlm.nih.gov/38515042/)
- Prakash J, Singh TB, Ghosh B, et al. **Changing picture of renal cortical necrosis in acute kidney injury in developing country.** *World J Nephrol*. 2015. [PMID 26558184](https://pubmed.ncbi.nlm.nih.gov/26558184/)

## 18. 국소 손상 · 구획증후군 · 근막절개
## Local injury, compartment syndrome and fasciotomy

Claim 5: the compartment intravenous antibody cannot reach.

- Spyres MB, et al. **Compartment Syndrome after Crotalid Envenomation in the United States: A Review of the North American Snakebite Registry from 2013 to 2021.** *Wilderness Environ Med*. 2023. [PMID 37474357](https://pubmed.ncbi.nlm.nih.gov/37474357/)
- Kim YH, et al. **Fasciotomy in compartment syndrome from snakebite.** *Arch Plast Surg*. 2019. [PMID 30685944](https://pubmed.ncbi.nlm.nih.gov/30685944/)
- Chen C, et al. **Severe compartment syndrome following snakebite in a man: A case report.** *Asian J Surg*. 2024. [PMID 39237419](https://pubmed.ncbi.nlm.nih.gov/39237419/)
- Sassoè-Pognetto M, et al. **Acute compartment syndrome and fasciotomy after a viper bite in Italy: a case report.** *Ital J Pediatr*. 2024. [PMID 38627836](https://pubmed.ncbi.nlm.nih.gov/38627836/)
- Hou YT, et al. **Prediction of Compartment Syndrome after Protobothrops mucrosquamatus Snakebite by Diastolic Retrograde Arterial Flow: A Case Report.** *Medicina (Kaunas)*. 2022. [PMID 35893111](https://pubmed.ncbi.nlm.nih.gov/35893111/)

## 19. 러셀살무사 — 이 모델의 기준 종
## Daboia russelii — the model's reference species

- Kularatne SA. **Epidemiology and clinical picture of the Russell's viper (Daboia russelii russelii) bite in Anuradhapura, Sri Lanka: a prospective study of 336 patients.** *Southeast Asian J Trop Med Public Health*. 2003. [PMID 15115100](https://pubmed.ncbi.nlm.nih.gov/15115100/)
- Kularatne SA, et al. **Revisiting Russell's viper (Daboia russelii) bite in Sri Lanka: is abdominal pain an early feature of systemic envenoming?** *PLoS One*. 2014. [PMID 24587278](https://pubmed.ncbi.nlm.nih.gov/24587278/)

## 20. 효과가 없는 것 — 헤파린과 트라넥사믹산
## What does NOT work: heparin and tranexamic acid

**본 모델이 D1에서 계산으로 재현하는 임상 음성 결과입니다.** 헤파린은 캐스케이드
밖에서 생성되는 트롬빈을 규제할 수 없고, 트라넥사믹산은 피브리노겐 손실의 95%가
직접 효소적 절단이기 때문에 손댈 수 있는 부분이 5%뿐입니다.

- Tin Na Swe, et al. **Heparin therapy in Russell's viper bite victims with disseminated intravascular coagulation: a controlled trial.** *Southeast Asian J Trop Med Public Health*. 1992. [PMID 1345132](https://pubmed.ncbi.nlm.nih.gov/1345132/)
- Myint-Lwin, et al. **Heparin therapy in Russell's viper bite victims with impending DIC (a controlled trial).** *Southeast Asian J Trop Med Public Health*. 1989. [PMID 2532790](https://pubmed.ncbi.nlm.nih.gov/2532790/)
- Paul V, et al. **Trial of heparin in viper bites.** *J Assoc Physicians India*. 2003. [PMID 12725259](https://pubmed.ncbi.nlm.nih.gov/12725259/)
- Paul V, et al. **Trial of low molecular weight heparin in the treatment of viper bites.** *J Assoc Physicians India*. 2007. [PMID 17844693](https://pubmed.ncbi.nlm.nih.gov/17844693/)
- Huang YN, et al. **Effects of heparin on venom-induced consumption coagulopathy: a meta-analysis of randomized controlled trials.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749494](https://pubmed.ncbi.nlm.nih.gov/39749494/)
- Larréché S, et al. **[Evaluation of the Efficacy and Tolerance of Tranexamic Acid Combined with Inoserp™ PAN-AFRICA Antivenom in the Treatment of Hemorrhagic Syndrome].** *Med Trop Sante Int*. 2026. [PMID 42221451](https://pubmed.ncbi.nlm.nih.gov/42221451/)
- Yamamoto A, et al. **Therapeutic Effects of Single and Combined Anti-Disseminated Intravascular Coagulation (DIC) Drugs in a Rat Venom-Induced Consumption Coagulopathy Model.** *Toxins (Basel)*. 2026. [PMID 41893574](https://pubmed.ncbi.nlm.nih.gov/41893574/)
- Abraham SV, et al. **Hematotoxic Snakebite Victim with Trauma: The Role of Guided Transfusion, Rotational Thromboelastometry, and Tranexamic Acid.** *Wilderness Environ Med*. 2020. [PMID 33162320](https://pubmed.ncbi.nlm.nih.gov/33162320/)

## 21. 소분자 억제제 — 항독소가 못 가는 곳으로 가는 약
## Small-molecule inhibitors — the drug that goes where antibody cannot

The `VAR_*` compartments and `INH_PLA2` factor, and the D5 result, are built on
this programme. The BRAVO trial is the reason varespladib is in the ODEs at all
rather than only in the map.

- Gerardo CJ, et al. **Oral varespladib for the treatment of snakebite envenoming in India and the USA (BRAVO): a phase II randomised clinical trial.** *BMJ Glob Health*. 2024. [PMID 39442939](https://pubmed.ncbi.nlm.nih.gov/39442939/)
- Carter RW, et al. **The BRAVO Clinical Study Protocol: Oral Varespladib for Inhibition of Secretory Phospholipase A2 in the Treatment of Snakebite Envenoming.** *Toxins (Basel)*. 2022. [PMID 36668842](https://pubmed.ncbi.nlm.nih.gov/36668842/)
- Quiroz S, et al. **Inhibitory Effects of Varespladib, CP471474, and Their Potential Synergistic Activity on Bothrops asper and Crotalus durissus cumanensis Venoms.** *Molecules*. 2022. [PMID 36500682](https://pubmed.ncbi.nlm.nih.gov/36500682/)
- Woliver C, et al. **Phospholipase A2 inhibitor may shorten the duration of clinical signs in the treatment of neurotoxicity caused by eastern coral snake (Micrurus fulvius) envenomation.** *J Am Vet Med Assoc*. 2025. [PMID 40738154](https://pubmed.ncbi.nlm.nih.gov/40738154/)

## 22. 차세대 항독소 — 재조합 항체와 나노바디
## Next-generation antivenom

Where the ε vector stops being a fact of nature and becomes a design variable.

- Kini RM, et al. **Biosynthetic Oligoclonal Antivenom (BOA) for Snakebite and Next-Generation Treatments for Snakebite Victims.** *Toxins (Basel)*. 2018. [PMID 30551565](https://pubmed.ncbi.nlm.nih.gov/30551565/)
- Knudsen C, Laustsen AH. **Engineering and design considerations for next-generation snakebite antivenoms.** *Toxicon*. 2019. [PMID 31173790](https://pubmed.ncbi.nlm.nih.gov/31173790/)
- Fernandes CFC, et al. **Engineering of single-domain antibodies for next-generation snakebite antivenoms.** *Int J Biol Macromol*. 2021. [PMID 34118288](https://pubmed.ncbi.nlm.nih.gov/34118288/)
- Jenkins TP, Laustsen AH. **Cost of Manufacturing for Recombinant Snakebite Antivenoms.** *Front Bioeng Biotechnol*. 2020. [PMID 32766215](https://pubmed.ncbi.nlm.nih.gov/32766215/)
- Jenkins TP, et al. **Toxin Neutralization Using Alternative Binding Proteins.** *Toxins (Basel)*. 2019. [PMID 30658491](https://pubmed.ncbi.nlm.nih.gov/30658491/)
- Knudsen C, et al. **Novel Snakebite Therapeutics Must Be Tested in Appropriate Rescue Models to Robustly Assess Their Preclinical Efficacy.** *Toxins (Basel)*. 2020. [PMID 32824899](https://pubmed.ncbi.nlm.nih.gov/32824899/)
- Tabares Vélez S, et al. **Standard Quality Characteristics and Efficacy of a New Third-Generation Antivenom Developed in Colombia Covering Micrurus spp. Venoms.** *Toxins (Basel)*. 2024. [PMID 38668608](https://pubmed.ncbi.nlm.nih.gov/38668608/)

## 23. 응급처치와 림프 흐름 — 저장고를 조절하려는 시도
## First aid and lymphatic flow — trying to control the depot

The pressure-immobilisation node in cluster 2, and the model's ka being a
lymphatic rate constant rather than a capillary one, come from here.

- Parker-Cote J, Meggs WJ. **First Aid and Pre-Hospital Management of Venomous Snakebites.** *Trop Med Infect Dis*. 2018. [PMID 30274441](https://pubmed.ncbi.nlm.nih.gov/30274441/)
- van Helden DF, et al. **Pharmacological approaches that slow lymphatic flow as a snakebite first aid.** *PLoS Negl Trop Dis*. 2014. [PMID 24587472](https://pubmed.ncbi.nlm.nih.gov/24587472/)
- Howarth DM, Southee AE, Whyte IM. **Lymphatic flow rates and first-aid in simulated peripheral snake or spider envenomation.** *Med J Aust*. 1994. [PMID 7830641](https://pubmed.ncbi.nlm.nih.gov/7830641/)

## 24. 공급 · 비용 · 접근성 — 실제 결과를 결정하는 변수
## Supply, cost and access — the variable that actually decides outcome

A QSP model can compute the optimal vial count; it cannot compute whether the
vial is in the building. These papers are why the "no ventilator" and "no
dialysis" switches exist in the model at all.

- Hamza M, et al. **Cost-Effectiveness of Antivenoms for Snakebite Envenoming in 16 Countries in West Africa.** *PLoS Negl Trop Dis*. 2016. [PMID 27027633](https://pubmed.ncbi.nlm.nih.gov/27027633/)
- Habib AG, et al. **Cost-effectiveness of antivenoms for snakebite envenoming in Nigeria.** *PLoS Negl Trop Dis*. 2015. [PMID 25569252](https://pubmed.ncbi.nlm.nih.gov/25569252/)
- Ooms GI, et al. **Availability, affordability and stock-outs of commodities for the treatment of snakebite in Kenya.** *PLoS Negl Trop Dis*. 2021. [PMID 34398889](https://pubmed.ncbi.nlm.nih.gov/34398889/)
- Gampini S, et al. **Retrospective study on the incidence of envenomation and accessibility to antivenom in Burkina Faso.** *J Venom Anim Toxins Incl Trop Dis*. 2016. [PMID 26985188](https://pubmed.ncbi.nlm.nih.gov/26985188/)

## 25. 장기 후유증 — 사망률이 놓치는 부담
## Long-term sequelae — the burden mortality statistics miss

The model's `FIBR` (permanent renal scar) and `NEC` (19-day healing half-life)
states exist because of this literature: the patient who survives is not the same as
the patient who recovers.

- Waiddyanatha S, Silva A, Siribaddana S, Isbister GK. **Long-term Effects of Snake Envenoming.** *Toxins (Basel)*. 2019. [PMID 30935096](https://pubmed.ncbi.nlm.nih.gov/30935096/)
- Jayawardana S, Gnanathasan A, Arambepola C, Chang T. **Long-term health complications following snake envenoming.** *J Multidiscip Healthc*. 2018. [PMID 29983571](https://pubmed.ncbi.nlm.nih.gov/29983571/)

## 26. 방법론 — mrgsolve와 QSP
## Methodology — mrgsolve and QSP

- Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol*. 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)

---

## 검증 (Verification)

모든 PMID는 다음 절차로 확보되었습니다.

1. 주제별 질의 약 50건을 NCBI `esearch.fcgi`(db=pubmed, sort=relevance)로 실행
2. 반환된 PMID를 `esummary.fcgi`로 조회해 저자·연도·저널·제목을 획득
3. **조회 결과에 없는 항목은 이 파일에 넣지 않았습니다.** 기억으로 쓴 인용,
   추정한 PMID, "아마 이 논문일 것" 같은 항목은 하나도 없습니다.
4. 각 PMID 링크는 `https://pubmed.ncbi.nlm.nih.gov/<PMID>/` 형식으로,
   위 조회에서 확인된 숫자만 사용했습니다.

이 파일이 뒷받침하지 **못하는** 것은 위의 "가져오지 않은 것" 절에 명시했습니다.
QSP 모델에서 가장 흔한 부정직은 없는 근거를 있는 척하는 것이 아니라, 보정
파라미터를 관측값처럼 배치하는 것입니다.
