# 심상성 여드름 (Acne vulgaris) — QSP 모델 참고문헌
### Reference list for `acn_qsp_model.dot` / `acn_mrgsolve_model.R` / `acn_shiny_app.R`

이 목록의 모든 PMID는 NCBI E-utilities로 조회하여 제목·저널·연도를 확인한 것입니다.
각 항목 아래 회색 주석은 **이 문헌이 모델의 어느 부분을 뒷받침하는지** 를 적은 것입니다.
(모델 파라미터는 이들 문헌의 보고값을 참고한 근사치이며 형식적 적합(fitting)을 거치지 않았습니다.)


## A. 개관 · 역학 · 질병 부담 (Overview, epidemiology, burden)

1. Williams HC, Dellavalle RP, Garner S. Acne vulgaris. *Lancet*. 2012;379:361-72. [PMID 21880356](https://pubmed.ncbi.nlm.nih.gov/21880356/)  
   <sub>Lancet 세미나. 유병률·자연경과·치료 계층화의 표준 요약. 모델의 '중등증' 기준 표현형 정의에 사용.</sub>
2. Moradi Tuchayi S, Makrantonaki E, Ganceviciene R, et al. Acne vulgaris. *Nat Rev Dis Primers*. 2015;1:15033. [PMID 27227877](https://pubmed.ncbi.nlm.nih.gov/27227877/)  
   <sub>Nature Reviews Disease Primers. 네 개의 병인 기둥(과피지·과각화·C. acnes·염증) 구조 자체가 이 모델의 골격.</sub>
3. Yuan R, Long H, Li T, et al. Comparison of the Burden of Acne Vulgaris in China vs Globally: Insights from the Global Burden of Disease Study (2021) and Projections to 2050. *Acta Derm Venereol*. 2026;106. [PMID 42051025](https://pubmed.ncbi.nlm.nih.gov/42051025/)  
   <sub>GBD 기반 질병부담 비교. 청소년기 정점과 성인 지속형의 인구학적 배경.</sub>
4. Samuels DV, Rosenthal R, Lin R, et al. Acne vulgaris and risk of depression and anxiety: A meta-analytic review. *J Am Acad Dermatol*. 2020;83:532-541. [PMID 32088269](https://pubmed.ncbi.nlm.nih.gov/32088269/)  
   <sub>우울·불안 위험 메타분석. Shiny 앱의 QoL 맥락과 조기 치료 논거.</sub>
5. Veldi VDK, Metta AK, Metta S, et al. Living With Acne Vulgaris in Young Adults: A Holistic Examination of Its Impact on Quality of Life Using the Dermatology Life Quality Index (DLQI). *Cureus*. 2025;17:e77167. [PMID 39925569](https://pubmed.ncbi.nlm.nih.gov/39925569/)  
   <sub>청년층 여드름의 삶의 질 영향에 대한 질적·정량적 검토.</sub>
6. Salari Y, Latt M, Lau E, et al. The Prevalence and Burden of Truncal Acne. *Dermatol Ther (Heidelb)*. 2026;16:3245-3251. [PMID 41963698](https://pubmed.ncbi.nlm.nih.gov/41963698/)  
   <sub>체간 여드름의 유병률과 부담. 트리파로텐 체간 적응증 시나리오의 배경.</sub>
7. Friedman GD. Twin studies of disease heritability based on medical records: application to acne vulgaris. *Acta Genet Med Gemellol (Roma)*. 1984;33:487-95. [PMID 6241419](https://pubmed.ncbi.nlm.nih.gov/6241419/)  
   <sub>의무기록 기반 쌍둥이 연구 — 여드름 유전율. 모델의 체질 파라미터 SEVX/LESX 도입 근거.</sub>
8. Maxwell J, Mitchell BL, Du-Harpur X, et al. Genome-wide association meta-regression identifies stem cell lineage orchestration as a key driver of acne risk. *medRxiv*. 2025. [PMID 40666320](https://pubmed.ncbi.nlm.nih.gov/40666320/)  
   <sub>여드름 GWAS 메타회귀 — 모낭 줄기세포 계통 유전자. 기계론적 지도 cluster 1.</sub>

## B. 병태생리 총론 · 미세면포 · 각질화 (Pathogenesis, the microcomedone)

9. Oulès B, Saurat JH. Strategic Targets in Acne, Update 2025: The Microcomedone Is Not Just a Plug, It Is an Egg. *Dermatology*. 2026;242:8-14. [PMID 40924653](https://pubmed.ncbi.nlm.nih.gov/40924653/)  
   <sub>2025 업데이트. '미세면포는 단순한 마개가 아니다' — 모델에서 MC를 별도 상태변수로 둔 직접적 근거.</sub>
10. Jeremy AH, Holland DB, Roberts SG, et al. Inflammatory events are involved in acne lesion initiation. *J Invest Dermatol*. 2003;121:20-7. [PMID 12839559](https://pubmed.ncbi.nlm.nih.gov/12839559/)  
   <sub>Jeremy 등. 임상적으로 정상인 모낭에서 이미 염증 사건이 선행함을 보인 고전. de novo 염증성 구진 항(KDENOVO)의 근거.</sub>
11. Persson G, Johansson-Jänkänpää E, Ganceviciene R, et al. No evidence for follicular keratinocyte hyperproliferation in acne lesions as compared to autologous healthy hair follicles. *Exp Dermatol*. 2018;27:668-671. [PMID 29582469](https://pubmed.ncbi.nlm.nih.gov/29582469/)  
   <sub>여드름 병변에서 모낭 각질세포 과증식의 증거가 없다는 반대 연구. 모델이 KER를 '증식'이 아닌 '응집·분화 이상' 지수로 다루는 이유.</sub>
12. Perisho K, Wertz PW, Madison KC, et al. Fatty acids of acylceramides from comedones and from the skin surface of acne patients and control subjects. *J Invest Dermatol*. 1988;90:350-3. [PMID 2964492](https://pubmed.ncbi.nlm.nih.gov/2964492/)  
   <sub>면포와 여드름 피부 표면 아실세라마이드의 지방산 조성 — 리놀레산 희석 가설(상태변수 LA)의 실측 근거.</sub>
13. Das S, Reynolds RV. Recent advances in acne pathogenesis: implications for therapy. *Am J Clin Dermatol*. 2014;15:479-88. [PMID 25388823](https://pubmed.ncbi.nlm.nih.gov/25388823/)  
   <sub>여드름 발병기전의 최근 진전과 치료적 함의 — 기둥 간 상호작용 정리.</sub>

## C. 안드로겐 축과 피지선 생물학 (Androgen axis, sebaceous biology)

14. Chen WC, Zouboulis CC. Hormones and the pilosebaceous unit. *Dermatoendocrinol*. 2009;1:81-6. [PMID 20224689](https://pubmed.ncbi.nlm.nih.gov/20224689/)  
   <sub>호르몬과 모낭피지선 단위. ARS → SGM/LIP 연결의 총론.</sub>
15. Picardo M, Ottaviani M, Camera E, et al. Sebaceous gland lipids. *Dermatoendocrinol*. 2009;1:68-71. [PMID 20224686](https://pubmed.ncbi.nlm.nih.gov/20224686/)  
   <sub>피지선 지질 총설 — 스쿠알렌·왁스에스터·트리글리세라이드 조성비(지도 cluster 4).</sub>
16. Khondker L, Khan SI. Acne vulgaris related to androgens - a review. *Mymensingh Med J*. 2014;23:181-5. [PMID 24584396](https://pubmed.ncbi.nlm.nih.gov/24584396/)  
   <sub>안드로겐과 여드름 총설. 유리 안드로겐 지수(FAI) 사용 근거.</sub>
17. Samson M, Labrie F, Zouboulis CC, et al. Biosynthesis of dihydrotestosterone by a pathway that does not require testosterone as an intermediate in the SZ95 sebaceous gland cell line. *J Invest Dermatol*. 2010;130:602-4. [PMID 19812596](https://pubmed.ncbi.nlm.nih.gov/19812596/)  
   <sub>피지세포에서 테스토스테론을 거치지 않는 DHT 생합성 경로 — 국소 스테로이드 생성 노드.</sub>
18. Leyden J, Bergfeld W, Drake L, et al. A systemic type I 5 alpha-reductase inhibitor is ineffective in the treatment of acne vulgaris. *J Am Acad Dermatol*. 2004;50:443-7. [PMID 14988688](https://pubmed.ncbi.nlm.nih.gov/14988688/)  
   <sub>전신 1형 5α-환원효소 억제제가 여드름에 무효했던 시험. 모델에서 E5ARI 기본값을 0으로 두고 5AR 단독 차단을 약한 지렛대로 취급하는 이유.</sub>
19. Cotterill JA, Cunliffe WJ, Williamson B. Severity of acne and sebum excretion rate. *Br J Dermatol*. 1971;85:93-4. [PMID 4254153](https://pubmed.ncbi.nlm.nih.gov/4254153/)  
   <sub>여드름 중증도와 피지 분비율의 관계 — SER 기준값(정상 ~0.9 vs 여드름 ~2 µg/cm²/min) 보정.</sub>
20. Pan J, Wang Q, Tu P. A Topical Medication of All-Trans Retinoic Acid Reduces Sebum Excretion Rate in Patients With Forehead Acne. *Am J Ther*. 2017;24:e207-e212. [PMID 26872139](https://pubmed.ncbi.nlm.nih.gov/26872139/)  
   <sub>국소 all-trans 레티노산의 피지 분비율 감소 관찰 — 레티노이드의 피지 효과 크기 상한 설정.</sub>
21. Chiba K, Yoshizawa K, Makino I, et al. Comedogenicity of squalene monohydroperoxide in the skin after topical application. *J Toxicol Sci*. 2000;25:77-83. [PMID 10845185](https://pubmed.ncbi.nlm.nih.gov/10845185/)  
   <sub>스쿠알렌 모노하이드로퍼옥사이드의 면포유발성 — SQOX → 면포 형성/개방면포 산화 착색 경로.</sub>
22. Zouboulis CC, Oeff MK, Hiroi N, et al. Involvement of Pattern Recognition Receptors in the Direct Influence of Bacterial Components and Standard Antiacne Compounds on Human Sebaceous Gland Cells. *Skin Pharmacol Physiol*. 2021;34:19-29. [PMID 33601383](https://pubmed.ncbi.nlm.nih.gov/33601383/)  
   <sub>피지세포의 패턴인식수용체 관여 — 피지선이 단순 분비 기관이 아니라 면역 활성 조직임.</sub>

## D. 영양-대사 신호 (Insulin · IGF-1 · mTORC1 · FoxO1)

23. Cappel M, Mauger D, Thiboutot D. Correlation between serum levels of insulin-like growth factor 1, dehydroepiandrosterone sulfate, and dihydrotestosterone and acne lesion counts in adult women. *Arch Dermatol*. 2005;141:333-8. [PMID 15781674](https://pubmed.ncbi.nlm.nih.gov/15781674/)  
   <sub>성인 여성에서 IGF-1·DHEAS·DHT와 병변 수의 상관 — IGF1 → LIP 가중치(WIGF) 설정.</sub>
24. Smith R, Mann N, Mäkeläinen H, et al. A pilot study to determine the short-term effects of a low glycemic load diet on hormonal markers of acne: a nonrandomized, parallel, controlled feeding trial. *Mol Nutr Food Res*. 2008;52:718-26. [PMID 18496812](https://pubmed.ncbi.nlm.nih.gov/18496812/)  
   <sub>저혈당부하 식이의 무작위 예비연구 — GLYLOAD 파라미터의 효과 크기 상한.</sub>
25. Adebamowo CA, Spiegelman D, Berkey CS, et al. Milk consumption and acne in teenaged boys. *J Am Acad Dermatol*. 2008;58:787-93. [PMID 18194824](https://pubmed.ncbi.nlm.nih.gov/18194824/)  
   <sub>청소년 남성의 우유 섭취와 여드름 — DAIRY 파라미터.</sub>
26. Juhl CR, Bergholdt HKM, Miller IM, et al. Dairy Intake and Acne Vulgaris: A Systematic Review and Meta-Analysis of 78,529 Children, Adolescents, and Young Adults. *Nutrients*. 2018;10. [PMID 30096883](https://pubmed.ncbi.nlm.nih.gov/30096883/)  
   <sub>유제품 섭취와 여드름: 78,000명 대상 체계적 문헌고찰·메타분석.</sub>
27. Melnik BC. Lifetime Impact of Cow's Milk on Overactivation of mTORC1: From Fetal to Childhood Overgrowth, Acne, Diabetes, Cancers, and Neurodegeneration. *Biomolecules*. 2021;11. [PMID 33803410](https://pubmed.ncbi.nlm.nih.gov/33803410/)  
   <sub>우유에 의한 mTORC1 과활성화 — 류신/BCAA 노드와 NUTR 항.</sub>
28. Agamia NF, Hussein OM, Abdelmaksoud RE, et al. Effect of oral isotretinoin on the nucleo-cytoplasmic distribution of FoxO1 and FoxO3 proteins in sebaceous glands of patients with acne vulgaris. *Exp Dermatol*. 2018;27:1344-1351. [PMID 30240097](https://pubmed.ncbi.nlm.nih.gov/30240097/)  
   <sub>경구 이소트레티노인이 FoxO1/FoxO3의 핵-세포질 분포를 바꾼다는 인체 자료 — 모델에서 이소트레티노인이 LIP(mTORC1/SREBP)를 직접 억제하는 항(EISOLIP)의 근거.</sub>

## E. Cutibacterium acnes — 계통형 · 바이오필름 · 대사산물

29. Fitz-Gibbon S, Tomida S, Chiu BH, et al. Propionibacterium acnes strain populations in the human skin microbiome associated with acne. *J Invest Dermatol*. 2013;133:2152-60. [PMID 23337890](https://pubmed.ncbi.nlm.nih.gov/23337890/)  
   <sub>Fitz-Gibbon 등의 이정표 연구. 균 '양'이 아니라 '균주 구성'이 여드름과 연관 — 모델이 PHYLIA(IA1 우세도)를 별도 파라미터로 둔 이유.</sub>
30. McDowell A, Perry AL, Lambert PA, et al. A new phylogenetic group of Propionibacterium acnes. *J Med Microbiol*. 2008;57:218-224. [PMID 18201989](https://pubmed.ncbi.nlm.nih.gov/18201989/)  
   <sub>C. acnes의 새로운 계통군 기술.</sub>
31. Dagnelie MA, Corvec S, Saint-Jean M, et al. Cutibacterium acnes phylotypes diversity loss: a trigger for skin inflammatory process. *J Eur Acad Dermatol Venereol*. 2019;33:2340-2348. [PMID 31299116](https://pubmed.ncbi.nlm.nih.gov/31299116/)  
   <sub>계통형 다양성 상실이 피부 염증의 방아쇠라는 직접 증거.</sub>
32. Dreno B, Dekio I, Baldwin H, et al. Acne microbiome: From phyla to phylotypes. *J Eur Acad Dermatol Venereol*. 2024;38:657-664. [PMID 37777343](https://pubmed.ncbi.nlm.nih.gov/37777343/)  
   <sub>여드름 마이크로바이옴: 문(phyla)에서 계통형까지 — 최신 총설.</sub>
33. Coenye T, Spittaels KJ, Achermann Y. The role of biofilm formation in the pathogenesis and antimicrobial susceptibility of Cutibacterium acnes. *Biofilm*. 2022;4:100063. [PMID 34950868](https://pubmed.ncbi.nlm.nih.gov/34950868/)  
   <sub>바이오필름 형성이 병인과 항균제 감수성에 미치는 역할 — CAB 구획과 보호계수 PROTB.</sub>
34. Ruffier d'Epenoux L, Fayoux E, Veziers J, et al. Biofilm of Cutibacterium acnes: a target of different active substances. *Int J Dermatol*. 2024;63:1541-1550. [PMID 38760974](https://pubmed.ncbi.nlm.nih.gov/38760974/)  
   <sub>C. acnes 바이오필름을 표적으로 하는 활성물질 검토.</sub>
35. Johnson T, Kang D, Barnard E, et al. Strain-Level Differences in Porphyrin Production and Regulation in Propionibacterium acnes Elucidate Disease Associations. *mSphere*. 2016;1. [PMID 27303708](https://pubmed.ncbi.nlm.nih.gov/27303708/)  
   <sub>균주 수준의 포르피린 생성 차이와 조절 — PORP 상태변수와 PHYLIA 의존성.</sub>
36. Nakatsuji T, Tang DC, Zhang L, et al. Propionibacterium acnes CAMP factor and host acid sphingomyelinase contribute to bacterial virulence: potential targets for inflammatory acne treatment. *PLoS One*. 2011;6:e14797. [PMID 21533261](https://pubmed.ncbi.nlm.nih.gov/21533261/)  
   <sub>CAMP factor와 숙주 산성 스핑고미엘리나제 — 각질세포 세포독성 노드.</sub>
37. Nakase K, Momose M, Yukawa T, et al. Development of skin sebum medium and inhibition of lipase activity in Cutibacterium acnes by oleic acid. *Access Microbiol*. 2022;4:acmi000397. [PMID 36415741](https://pubmed.ncbi.nlm.nih.gov/36415741/)  
   <sub>피지 배지 개발과 C. acnes 리파아제 활성 억제 — 리파아제 → 유리지방산(FFA) 경로.</sub>
38. Wang Y, Kao MS, Yu J, et al. A Precision Microbiome Approach Using Sucrose for Selective Augmentation of Staphylococcus epidermidis Fermentation against Propionibacterium acnes. *Int J Mol Sci*. 2016;17. [PMID 27834859](https://pubmed.ncbi.nlm.nih.gov/27834859/)  
   <sub>S. epidermidis 발효를 이용한 C. acnes 억제 — 길항 균총 노드.</sub>

## F. 선천면역 · 인플라마좀 · Th17 (Innate and adaptive inflammation)

39. Kim J. Review of the innate immune response in acne vulgaris: activation of Toll-like receptor 2 in acne triggers inflammatory cytokine responses. *Dermatology*. 2005;211:193-8. [PMID 16205063](https://pubmed.ncbi.nlm.nih.gov/16205063/)  
   <sub>여드름 선천면역 반응 총설 — TLR2 활성화. 상태변수 TLR의 근거.</sub>
40. Qin M, Pirouz A, Kim MH, et al. Propionibacterium acnes Induces IL-1β secretion via the NLRP3 inflammasome in human monocytes. *J Invest Dermatol*. 2014;134:381-388. [PMID 23884315](https://pubmed.ncbi.nlm.nih.gov/23884315/)  
   <sub>C. acnes가 인간 단핵구에서 NLRP3 인플라마좀을 통해 IL-1β를 분비시킴 — IL1B 생성항과 NLRP3G 2차 신호.</sub>
41. Nagy I, Pivarcsi A, Koreck A, et al. Distinct strains of Propionibacterium acnes induce selective human beta-defensin-2 and interleukin-8 expression in human keratinocytes through toll-like receptors. *J Invest Dermatol*. 2005;124:931-8. [PMID 15854033](https://pubmed.ncbi.nlm.nih.gov/15854033/)  
   <sub>균주에 따라 hBD-2와 IL-8 발현이 선택적으로 유도됨 — 계통형 의존적 염증.</sub>
42. Agak GW, Qin M, Nobe J, et al. Propionibacterium acnes Induces an IL-17 Response in Acne Vulgaris that Is Regulated by Vitamin A and Vitamin D. *J Invest Dermatol*. 2014;134:366-373. [PMID 23924903](https://pubmed.ncbi.nlm.nih.gov/23924903/)  
   <sub>C. acnes가 IL-17 반응을 유도하며 비타민 A/D가 이를 조절 — IL17 구획.</sub>
43. Kelhälä HL, Palatsi R, Fyhrquist N, et al. IL-17/Th17 pathway is activated in acne lesions. *PLoS One*. 2014;9:e105238. [PMID 25153527](https://pubmed.ncbi.nlm.nih.gov/25153527/)  
   <sub>여드름 병변에서 IL-17/Th17 경로가 실제로 활성화되어 있음을 보인 인체 조직 연구.</sub>
44. Kang S, Cho S, Chung JH, et al. Inflammation and extracellular matrix degradation mediated by activated transcription factors nuclear factor-kappaB and activator protein-1 in inflammatory acne lesions in vivo. *Am J Pathol*. 2005;166:1691-9. [PMID 15920154](https://pubmed.ncbi.nlm.nih.gov/15920154/)  
   <sub>여드름 병변 내 NF-κB·AP-1 활성화와 기질금속단백분해효소 — MMP 상태변수와 흉터 경로.</sub>

## G. 신경내분비 · 스트레스 축 (Neuroendocrine axis)

45. Krause K, Schnitger A, Fimmel S, et al. Corticotropin-releasing hormone skin signaling is receptor-mediated and is predominant in the sebaceous glands. *Horm Metab Res*. 2007;39:166-70. [PMID 17326013](https://pubmed.ncbi.nlm.nih.gov/17326013/)  
   <sub>CRH 피부 신호전달이 수용체 매개이며 피지선에서 우세 — CRH → 피지 생성(GCRH) 항.</sub>
46. Chiu A, Chon SY, Kimball AB. The response of skin disease to stress: changes in the severity of acne vulgaris as affected by examination stress. *Arch Dermatol*. 2003;139:897-900. [PMID 12873885](https://pubmed.ncbi.nlm.nih.gov/12873885/)  
   <sub>시험 스트레스에 따른 여드름 중증도 변화 — STRESS 파라미터의 임상적 근거.</sub>

## H. 흉터 · 염증후 색소침착 (Scarring and PIH)

47. Connolly D, Vu HL, Mariwalla K, et al. Acne Scarring-Pathogenesis, Evaluation, and Treatment Options. *J Clin Aesthet Dermatol*. 2017;10:12-23. [PMID 29344322](https://pubmed.ncbi.nlm.nih.gov/29344322/)  
   <sub>여드름 흉터의 발병기전·평가·치료 총설 — SCAR 상태변수 설계.</sub>
48. Tan J, Beissert S, Cook-Bolden F, et al. Evaluation of psychological well-being and social impact of atrophic acne scarring: A multinational, mixed-methods study. *JAAD Int*. 2022;6:43-50. [PMID 35005652](https://pubmed.ncbi.nlm.nih.gov/35005652/)  
   <sub>위축성 여드름 흉터의 심리적·사회적 영향.</sub>
49. Callender VD, St Surin-Lord S, Davis EC, et al. Postinflammatory hyperpigmentation: etiologic and therapeutic considerations. *Am J Clin Dermatol*. 2011;12:87-99. [PMID 21348540](https://pubmed.ncbi.nlm.nih.gov/21348540/)  
   <sub>염증후 과색소침착의 원인과 치료 — PIH 상태변수와 SKINPIG(피츠패트릭) 파라미터.</sub>
50. Dréno B, Bissonnette R, Gagné-Henley A, et al. Prevention and Reduction of Atrophic Acne Scars with Adapalene 0.3%/Benzoyl Peroxide 2.5% Gel in Subjects with Moderate or Severe Facial Acne: Results of a 6-Month Randomized, Vehicle-Controlled Trial Using Intra-Individual Comparison. *Am J Clin Dermatol*. 2018;19:275-286. [PMID 29549588](https://pubmed.ncbi.nlm.nih.gov/29549588/)  
   <sub>아다팔렌 0.3%/BPO에 의한 위축성 흉터 예방·감소 — '흉터를 줄이는 유일한 방법은 결절 존재 시간을 줄이는 것'이라는 모델 결론의 임상 대응물.</sub>

## I. 국소 치료 (Topical therapy)

51. Kolli SS, Pecone D, Pona A, et al. Topical Retinoids in Acne Vulgaris: A Systematic Review. *Am J Clin Dermatol*. 2019;20:345-365. [PMID 30674002](https://pubmed.ncbi.nlm.nih.gov/30674002/)  
   <sub>국소 레티노이드 체계적 문헌고찰 — 12주 병변 감소율 보정 목표.</sub>
52. Thiboutot DM, Weiss J, Bucko A, et al. Adapalene-benzoyl peroxide, a fixed-dose combination for the treatment of acne vulgaris: results of a multicenter, randomized double-blind, controlled study. *J Am Acad Dermatol*. 2007;57:791-9. [PMID 17655969](https://pubmed.ncbi.nlm.nih.gov/17655969/)  
   <sub>아다팔렌–BPO 고정용량 복합제 — 시나리오 5의 근거.</sub>
53. Del Rosso JQ, Lain E, York JP, et al. Trifarotene 0.005% Cream in the Treatment of Facial and Truncal Acne Vulgaris in Patients with Skin of Color: a Case Series. *Dermatol Ther (Heidelb)*. 2022;12:2189-2200. [PMID 35994159](https://pubmed.ncbi.nlm.nih.gov/35994159/)  
   <sub>트리파로텐 0.005% 크림의 안면·체간 여드름 시험 — ACN_RETPOT의 상대역가 설정.</sub>
54. Tanghetti E, Dhawan S, Green L, et al. Randomized comparison of the safety and efficacy of tazarotene 0.1% cream and adapalene 0.3% gel in the treatment of patients with at least moderate facial acne vulgaris. *J Drugs Dermatol*. 2010;9:549-58. [PMID 20480800](https://pubmed.ncbi.nlm.nih.gov/20480800/)  
   <sub>타자로텐 0.1% 크림과 아다팔렌의 무작위 비교.</sub>
55. Fulton JE Jr, Farzad-Bakshandeh A, Bradley S. Studies on the mechanism of action to topical benzoyl peroxide and vitamin A acid in acne vulgaris. *J Cutan Pathol*. 1974;1:191-200. [PMID 4283462](https://pubmed.ncbi.nlm.nih.gov/4283462/)  
   <sub>국소 벤조일퍼옥사이드와 비타민 A 산의 작용기전 고전 연구.</sub>
56. Leyden JJ, Preston N, Osborn C, et al. In-vivo Effectiveness of Adapalene 0.1%/Benzoyl Peroxide 2.5% Gel on Antibiotic-sensitive and Resistant Propionibacterium acnes. *J Clin Aesthet Dermatol*. 2011;4:22-6. [PMID 21607190](https://pubmed.ncbi.nlm.nih.gov/21607190/)  
   <sub>아다팔렌 0.1%/BPO 2.5% 젤의 항생제 내성 P. acnes에 대한 생체 내 효과 — BPO가 내성 분획을 정화한다는 KBPORES 항의 직접 근거.</sub>
57. Webster G, Thiboutot DM, Chen DM, et al. Impact of a fixed combination of clindamycin phosphate 1.2%-benzoyl peroxide 2.5% aqueous gel on health-related quality of life in moderate to severe acne vulgaris. *Cutis*. 2010;86:263-7. [PMID 21214129](https://pubmed.ncbi.nlm.nih.gov/21214129/)  
   <sub>클린다마이신 1.2%–BPO 고정 복합제의 영향 — 시나리오 7.</sub>
58. Gold M, Lain T, Harper JC, et al. Efficacy and Safety of Clindamycin Phosphate 1.2%/Adapalene 0.15%/Benzoyl Peroxide 3.1% Gel: Post Hoc Analysis by Baseline Disease Severity. *Dermatol Ther (Heidelb)*. 2025;15:1867-1882. [PMID 40377868](https://pubmed.ncbi.nlm.nih.gov/40377868/)  
   <sub>클린다마이신 1.2%/아다팔렌 0.15%/BPO 3.1% 삼중 고정복합제의 유효성·안전성.</sub>
59. Draelos ZD, Rodriguez DA, Kempers SE, et al. Treatment Response With Once-Daily Topical Dapsone Gel, 7.5% for Acne Vulgaris: Subgroup Analysis of Pooled Data from Two Randomized, Double-Blind Stu. *J Drugs Dermatol*. 2017;16:591-598. [PMID 28686777](https://pubmed.ncbi.nlm.nih.gov/28686777/)  
   <sub>답손 겔 7.5% 1일 1회 — 호중구(MPO) 억제 항(EDAP).</sub>
60. Gowda CM, Wairkar S. Azelaic acid-based lyotropic liquid crystals gel for acne vulgaris: Formulation optimization, antimicrobial activity and dermatopharmacokinetic study. *Int J Pharm*. 2024;667:124879. [PMID 39490554](https://pubmed.ncbi.nlm.nih.gov/39490554/)  
   <sub>아젤라산 제형 연구 — 각질정상화·항균·티로시나아제 억제의 다중 작용.</sub>
61. Makrantonaki E. Topical minocycline foam: A new option for acne treatment?. *J Eur Acad Dermatol Venereol*. 2025;39:885-886. [PMID 40277219](https://pubmed.ncbi.nlm.nih.gov/40277219/)  
   <sub>국소 미노사이클린 폼 — 전신 노출을 최소화한 국소 항생제 옵션.</sub>
62. Stein Gold L, Kwong P, Draelos Z, et al. Impact of Topical Vehicles and Cutaneous Delivery Technologies on Patient Adherence and Treatment Outcomes in Acne and Rosacea. *J Clin Aesthet Dermatol*. 2023;16:26-34. [PMID 37288283](https://pubmed.ncbi.nlm.nih.gov/37288283/)  
   <sub>제형(vehicle)과 전달 기술이 순응도에 미치는 영향 — ADHERE 파라미터와 자극-순응도 되먹임.</sub>

## J. 전신 항생제 · 항생제 내성과 청지기 (Systemic antibiotics, resistance, stewardship)

63. Garner SE, Eady A, Bennett C, et al. Minocycline for acne vulgaris: efficacy and safety. *Cochrane Database Syst Rev*. 2012;2012:CD002086. [PMID 22895927](https://pubmed.ncbi.nlm.nih.gov/22895927/)  
   <sub>미노사이클린의 유효성·안전성(Cochrane) — 계열 내 우월성 근거 부재.</sub>
64. Zhanel G, Critchley I, Lin LY, et al. Microbiological Profile of Sarecycline, a Novel Targeted Spectrum Tetracycline for the Treatment of Acne Vulgaris. *Antimicrob Agents Chemother*. 2019;63. [PMID 30397052](https://pubmed.ncbi.nlm.nih.gov/30397052/)  
   <sub>사레사이클린의 미생물학적 프로파일 — 협범위(NARROW) 파라미터.</sub>
65. Zhang J, He L, Chen X, et al. Efficacy and Safety of Sarecycline in Chinese Patients with Moderate-to-Severe Acne Vulgaris: Randomized Phase 3 Clinical Trial with Open-Label Follow-Up. *Dermatol Ther (Heidelb)*. 2025;15:3285-3300. [PMID 40877730](https://pubmed.ncbi.nlm.nih.gov/40877730/)  
   <sub>중등증–중증 여드름에서 사레사이클린의 유효성·안전성.</sub>
66. Toossi P, Farshchian M, Malekzad F, et al. Subantimicrobial-dose doxycycline in the treatment of moderate facial acne. *J Drugs Dermatol*. 2008;7:1149-52. [PMID 19137768](https://pubmed.ncbi.nlm.nih.gov/19137768/)  
   <sub>아항균 용량 독시사이클린의 중등증 안면 여드름 치료 — 항균 arm 없이도 효과가 있다는 핵심 근거.</sub>
67. Moore A, Ling M, Bucko A, et al. Efficacy and Safety of Subantimicrobial Dose, Modified-Release Doxycycline 40 mg Versus Doxycycline 100 mg Versus Placebo for the treatment of Inflammatory Lesions in Moderate and Severe Acne: A Randomized, Double-Blinded, Controlled Study. *J Drugs Dermatol*. 2015;14:581-6. [PMID 26091383](https://pubmed.ncbi.nlm.nih.gov/26091383/)  
   <sub>서방형 아항균 용량 독시사이클린 40 mg의 유효성·안전성 — EC5AI ≪ EC5TET 로 두 arm을 분리한 설계의 근거.</sub>
68. Ross JI, Snelling AM, Eady EA, et al. Phenotypic and genotypic characterization of antibiotic-resistant Propionibacterium acnes isolated from acne patients attending dermatology clinics in Europe, the U.S.A., Japan and Australia. *Br J Dermatol*. 2001;144:339-46. [PMID 11251569](https://pubmed.ncbi.nlm.nih.gov/11251569/)  
   <sub>여드름 환자 유래 항생제 내성 P. acnes의 표현형·유전형 특성.</sub>
69. Ross JI, Eady EA, Cove JH, et al. 16S rRNA mutation associated with tetracycline resistance in a gram-positive bacterium. *Antimicrob Agents Chemother*. 1998;42:1702-5. [PMID 9661007](https://pubmed.ncbi.nlm.nih.gov/9661007/)  
   <sub>테트라사이클린 내성과 연관된 16S rRNA 돌연변이 — RESF의 분자적 실체.</sub>
70. Zhu C, Wei B, Li Y, et al. Antibiotic resistance rates in Cutibacterium acnes isolated from patients with acne vulgaris: a systematic review and meta-analysis. *Front Microbiol*. 2025;16:1565111. [PMID 40535003](https://pubmed.ncbi.nlm.nih.gov/40535003/)  
   <sub>C. acnes 항생제 내성률의 최신 조사.</sub>
71. Barbieri JS, Spaccarelli N, Margolis DJ, et al. Approaches to limit systemic antibiotic use in acne: Systemic alternatives, emerging topical therapies, dietary modification, and laser and light-based treatments. *J Am Acad Dermatol*. 2019;80:538-549. [PMID 30296534](https://pubmed.ncbi.nlm.nih.gov/30296534/)  
   <sub>여드름에서 전신 항생제 사용을 제한하는 접근법.</sub>
72. Rosenberg AL, Shah M, Del Rosso JQ, et al. Optimal Use Recommendations and Stewardship Principles with Oral Antibiotics in Acne Vulgaris Management: An Expert Consensus Panel. *J Clin Aesthet Dermatol*. 2025;18:21-29. [PMID 41640785](https://pubmed.ncbi.nlm.nih.gov/41640785/)  
   <sub>경구 항생제의 최적 사용 권고와 청지기 원칙 — 시나리오 6 vs 7 비교의 임상 대응물.</sub>

## K. 호르몬 요법 (Hormonal therapy)

73. Arowojolu AO, Gallo MF, Lopez LM, et al. Combined oral contraceptive pills for treatment of acne. *Cochrane Database Syst Rev*. 2012;2012:CD004425. [PMID 22786490](https://pubmed.ncbi.nlm.nih.gov/22786490/)  
   <sub>복합경구피임제의 여드름 치료 효과(Cochrane) — 시나리오 12의 효과 크기.</sub>
74. Koo EB, Petersen TD, Kimball AB. Meta-analysis comparing efficacy of antibiotics versus oral contraceptives in acne vulgaris. *J Am Acad Dermatol*. 2014;71:450-9. [PMID 24880665](https://pubmed.ncbi.nlm.nih.gov/24880665/)  
   <sub>항생제 대 경구피임제의 효능 비교 메타분석 — 6개월 시점에서 동등.</sub>
75. Santer M, Lawrence M, Pyne S, et al. Clinical and cost-effectiveness of spironolactone in treating persistent facial acne in women: SAFA double-blinded RCT. *Health Technol Assess*. 2024;28:1-86. [PMID 39268864](https://pubmed.ncbi.nlm.nih.gov/39268864/)  
   <sub>SAFA 무작위 이중맹검 시험 — 성인 여성 여드름에서 스피로놀락톤의 임상·비용 효과.</sub>
76. Kow CS, Ramachandram DS, Hasan SS, et al. Spironolactone for the Treatment of Moderate to Severe Acne in Adult Women: A Systematic Review and Meta-Analysis of Randomised Controlled Trials. *Australas J Dermatol*. 2025;66:165-168. [PMID 39912292](https://pubmed.ncbi.nlm.nih.gov/39912292/)  
   <sub>성인 여성 중등증–중증 여드름에서 스피로놀락톤.</sub>
77. Basu P, Elman SA, Abudu B, et al. High-dose spironolactone for acne in patients with polycystic ovarian syndrome: A single-institution retrospective study. *J Am Acad Dermatol*. 2021;85:740-741. [PMID 31400460](https://pubmed.ncbi.nlm.nih.gov/31400460/)  
   <sub>PCOS 환자에서 고용량 스피로놀락톤 — pcos 표현형 시나리오.</sub>
78. Patiyasikunt M, Chancheewa B, Asawanonda P, et al. Efficacy and tolerability of low-dose spironolactone and topical benzoyl peroxide in adult female acne: A randomized, double-blind, placebo-controlled trial. *J Dermatol*. 2020;47:1411-1416. [PMID 32857471](https://pubmed.ncbi.nlm.nih.gov/32857471/)  
   <sub>저용량 스피로놀락톤과 국소 BPO 병용의 유효성·내약성.</sub>
79. Plovanich M, Weng QY, Mostaghimi A. Low Usefulness of Potassium Monitoring Among Healthy Young Women Taking Spironolactone for Acne. *JAMA Dermatol*. 2015;151:941-4. [PMID 25796182](https://pubmed.ncbi.nlm.nih.gov/25796182/)  
   <sub>건강한 젊은 여성에서 칼륨 모니터링의 낮은 유용성 — KSER 추적 상태변수의 해석.</sub>
80. Rosette C, Agan FJ, Mazzetti A, et al. Cortexolone 17α-propionate (Clascoterone) Is a Novel Androgen Receptor Antagonist that Inhibits Production of Lipids and Inflammatory Cytokines from Sebocytes In Vitro. *J Drugs Dermatol*. 2019;18:412-418. [PMID 31141847](https://pubmed.ncbi.nlm.nih.gov/31141847/)  
   <sub>클라스코테론(코르텍솔론 17α-프로피오네이트)이 피지세포 AR 길항제로서 지질·염증성 사이토카인 생성을 억제 — KICLAS 경쟁적 길항 항.</sub>
81. Alkhodaidi ST, Al Hawsawi KA, Alkhudaidi IT, et al. Efficacy and safety of topical clascoterone cream for treatment of acne vulgaris: A systematic review and meta-analysis of randomized placebo-controlled trials. *Dermatol Ther*. 2021;34:e14609. [PMID 33258536](https://pubmed.ncbi.nlm.nih.gov/33258536/)  
   <sub>국소 클라스코테론 크림의 유효성·안전성 — 시나리오 11.</sub>
82. Damoulaki E, Sioutis D, Sarli V, et al. Polycystic Ovary Syndrome-Associated Acne: The Interplay of Hyperandrogenism, Insulin Resistance, and Therapeutic Strategies. *Cureus*. 2025;17:e98103. [PMID 41473651](https://pubmed.ncbi.nlm.nih.gov/41473651/)  
   <sub>PCOS 연관 여드름과 고안드로겐혈증의 상호작용.</sub>

## L. 이소트레티노인 — PK · 작용기전 · 누적용량 · 안전성

83. Colburn WA, Gibson DM, Wiens RE, et al. Food increases the bioavailability of isotretinoin. *J Clin Pharmacol*. 1983;23:534-9. [PMID 6582073](https://pubmed.ncbi.nlm.nih.gov/6582073/)  
   <sub>음식이 이소트레티노인의 생체이용률을 높인다 — FOOD/FOODEF 파라미터의 원 근거.</sub>
84. Webster GF, Leyden JJ, Gross JA. Results of a Phase III, double-blind, randomized, parallel-group, non-inferiority study evaluating the safety and efficacy of isotretinoin-Lidose in patients with severe recalcitrant nodular acne. *J Drugs Dermatol*. 2014;13:665-70. [PMID 24918555](https://pubmed.ncbi.nlm.nih.gov/24918555/)  
   <sub>Lidose 이소트레티노인의 제3상 비열등성 시험 — LIDOSE 플래그(식이 의존성 소실).</sub>
85. Almond-Roesler B, Blume-Peytavi U, Bisson S, et al. Monitoring of isotretinoin therapy by measuring the plasma levels of isotretinoin and 4-oxo-isotretinoin. A useful tool for management of severe acne. *Dermatology*. 1998;196:176-81. [PMID 9557257](https://pubmed.ncbi.nlm.nih.gov/9557257/)  
   <sub>이소트레티노인과 4-옥소-이소트레티노인 혈장 농도 모니터링 — 대사체 구획과 POTOXO.</sub>
86. Landthaler M, Kummermehr J, Wagner A, et al. Inhibitory effects of 13-cis-retinoic acid on human sebaceous glands. *Arch Dermatol Res*. 1980;269:297-309. [PMID 6453562](https://pubmed.ncbi.nlm.nih.gov/6453562/)  
   <sub>13-cis-레티노산의 인간 피지선 억제 효과 — SER 90% 감소 목표치.</sub>
87. Nelson AM, Gilliland KL, Cong Z, et al. 13-cis Retinoic acid induces apoptosis and cell cycle arrest in human SEB-1 sebocytes. *J Invest Dermatol*. 2006;126:2178-89. [PMID 16575387](https://pubmed.ncbi.nlm.nih.gov/16575387/)  
   <sub>13-cis 레티노산이 인간 SEB-1 피지세포에서 세포자멸과 세포주기 정지를 유도 — KAPO 증강항(EISOAP).</sub>
88. Nelson AM, Zhao W, Gilliland KL, et al. Neutrophil gelatinase-associated lipocalin mediates 13-cis retinoic acid-induced apoptosis of human sebaceous gland cells. *J Clin Invest*. 2008;118:1468-78. [PMID 18317594](https://pubmed.ncbi.nlm.nih.gov/18317594/)  
   <sub>NGAL(리포칼린-2)이 13-cis 레티노산 유도 피지세포 자멸을 매개 — 지도 cluster 15의 분자 경로.</sub>
89. Layton AM, Knaggs H, Taylor J, et al. Isotretinoin for acne vulgaris--10 years later: a safe and successful treatment. *Br J Dermatol*. 1993;129:292-6. [PMID 8286227](https://pubmed.ncbi.nlm.nih.gov/8286227/)  
   <sub>이소트레티노인 10년 경험 — 재발의 임상적 결정인자.</sub>
90. Blasiak RC, Stamey CR, Burkhart CN, et al. High-dose isotretinoin treatment and the rate of retrial, relapse, and adverse effects in patients with acne vulgaris. *JAMA Dermatol*. 2013;149:1392-8. [PMID 24173086](https://pubmed.ncbi.nlm.nih.gov/24173086/)  
   <sub>고용량 이소트레티노인과 재시도·재발·이상반응 비율.</sub>
91. Borghi A, Mantovani L, Minghetti S, et al. Low-cumulative dose isotretinoin treatment in mild-to-moderate acne: efficacy in achieving stable remission. *J Eur Acad Dermatol Venereol*. 2011;25:1094-8. [PMID 21198947](https://pubmed.ncbi.nlm.nih.gov/21198947/)  
   <sub>저누적용량 이소트레티노인의 안정적 관해 달성 능력 — 시나리오 15의 대조.</sub>
92. Rademaker M. Making sense of the effects of the cumulative dose of isotretinoin in acne vulgaris. *Int J Dermatol*. 2016;55:518-23. [PMID 26471145](https://pubmed.ncbi.nlm.nih.gov/26471145/)  
   <sub>여드름에서 이소트레티노인 누적용량 효과의 해석 — CUMISO → DURAB 힐 함수(CD50 ≈ 85 mg/kg)의 근거.</sub>
93. Zane LT, Leyden WA, Marqueling AL, et al. A population-based analysis of laboratory abnormalities during isotretinoin therapy for acne vulgaris. *Arch Dermatol*. 2006;142:1016-22. [PMID 16924051](https://pubmed.ncbi.nlm.nih.gov/16924051/)  
   <sub>이소트레티노인 치료 중 검사실 이상의 인구기반 분석 — TG·ALT 상태변수 보정.</sub>
94. Xia E, Han J, Faletsky A, et al. Isotretinoin Laboratory Monitoring in Acne Treatment: A Delphi Consensus Study. *JAMA Dermatol*. 2022;158:942-948. [PMID 35704293](https://pubmed.ncbi.nlm.nih.gov/35704293/)  
   <sub>이소트레티노인 검사 모니터링에 관한 델파이 합의.</sub>
95. Huang YC, Cheng YC. Isotretinoin treatment for acne and risk of depression: A systematic review and meta-analysis. *J Am Acad Dermatol*. 2017;76:1068-1076.e9. [PMID 28291553](https://pubmed.ncbi.nlm.nih.gov/28291553/)  
   <sub>이소트레티노인과 우울증 위험의 체계적 문헌고찰·메타분석.</sub>
96. Coberly S, Lammer E, Alashari M. Retinoic acid embryopathy: case report and review of literature. *Pediatr Pathol Lab Med*. 1996;16:823-36. [PMID 9025880](https://pubmed.ncbi.nlm.nih.gov/9025880/)  
   <sub>레티노산 배아병증 — 최기형성과 iPLEDGE 요구사항의 근거.</sub>
97. Fraunfelder FW, Fraunfelder FT, Corbett JJ. Isotretinoin-associated intracranial hypertension. *Ophthalmology*. 2004;111:1248-50. [PMID 15177980](https://pubmed.ncbi.nlm.nih.gov/15177980/)  
   <sub>이소트레티노인 연관 두개내압상승 — 테트라사이클린 병용 금기의 근거.</sub>
98. Scaramuzzino L, Coronella L, Lauletta G, et al. Recalcitrant isotretinoin-induced acne fulminans successfully treated with oral dapsone. *Ital J Dermatol Venerol*. 2025;160:381-382. [PMID 40292612](https://pubmed.ncbi.nlm.nih.gov/40292612/)  
   <sub>이소트레티노인 유발 전격성 여드름 — 중증 결절성에서 초기 악화와 스테로이드 선행의 필요성.</sub>

## M. 진료 지침 (Guidelines)

99. Reynolds RV, Yeung H, Cheng CE, et al. Guidelines of care for the management of acne vulgaris. *J Am Acad Dermatol*. 2024;90:1006.e1-1006.e30. [PMID 38300170](https://pubmed.ncbi.nlm.nih.gov/38300170/)  
   <sub>미국피부과학회(AAD) 2024 여드름 진료지침 — 시나리오 라이브러리의 요법 구성 기준.</sub>
100. Zaenglein AL, Pathy AL, Schlosser BJ, et al. Guidelines of care for the management of acne vulgaris. *J Am Acad Dermatol*. 2016;74:945-73.e33. [PMID 26897386](https://pubmed.ncbi.nlm.nih.gov/26897386/)  
   <sub>AAD 2016 여드름 진료지침.</sub>
101. Nast A, Al Wattar BH, Beylot Barry M, et al. Update of the EuroGuiDerm evidence-based guideline for the treatment of acne-Short version. *J Eur Acad Dermatol Venereol*. 2026;40:1162-1172. [PMID 41847993](https://pubmed.ncbi.nlm.nih.gov/41847993/)  
   <sub>EuroGuiDerm 근거기반 여드름 치료 지침 업데이트.</sub>
102. Nast A, Rosumeck S, Sammain A, et al. Methods report on the development of the European S3 guidelines for the treatment of acne. *J Eur Acad Dermatol Venereol*. 2012;26 Suppl 1:e1-41. [PMID 22356612](https://pubmed.ncbi.nlm.nih.gov/22356612/)  
   <sub>유럽 S3 지침 개발 방법론 보고.</sub>
103. Tan J, Alexis A, Baldwin H, et al. The Personalised Acne Care Pathway-Recommendations to guide longitudinal management from the Personalising Acne: Consensus of Experts. *JAAD Int*. 2021;5:101-111. [PMID 34816135](https://pubmed.ncbi.nlm.nih.gov/34816135/)  
   <sub>개인화 여드름 관리 경로 권고 — 표현형별 치료 선택.</sub>
104. Thiboutot DM, Shalita AR, Yamauchi PS, et al. Adapalene gel, 0.1%, as maintenance therapy for acne vulgaris: a randomized, controlled, investigator-blind follow-up of a recent combination study. *Arch Dermatol*. 2006;142:597-602. [PMID 16702497](https://pubmed.ncbi.nlm.nih.gov/16702497/)  
   <sub>아다팔렌 0.1% 젤의 유지요법 무작위 시험 — 시나리오 17(유도 후 유지)의 직접 근거.</sub>

## N. 변이 표현형 · 감별 (Variant phenotypes)

105. Gorji M, Joseph J, Pavlakis N, et al. Prevention and management of acneiform rash associated with EGFR inhibitor therapy: A systematic review and meta-analysis. *Asia Pac J Clin Oncol*. 2022;18:526-539. [PMID 35352492](https://pubmed.ncbi.nlm.nih.gov/35352492/)  
   <sub>EGFR 억제제 연관 여드름양 발진의 예방과 관리.</sub>
106. Wang M, Zhang P, Shen M. A de novo heterozygous PSTPIP1 variant associated with PAPA syndrome: a Chinese case report and literature review. *Front Genet*. 2026;17:1825761. [PMID 42358432](https://pubmed.ncbi.nlm.nih.gov/42358432/)  
   <sub>PSTPIP1 변이와 PAPA 증후군 — 자가염증성 중복 표현형.</sub>

## O. QSP 방법론 · 구현 (QSP methodology and implementation)

107. Kapitanov GI, Earp JC, Gadkar K, et al. Bridging the Gap: Integrating Quantitative Systems Pharmacology and Pharmacometrics in Drug Development. *Clin Pharmacol Ther*. 2026;119:830-833. [PMID 41472478](https://pubmed.ncbi.nlm.nih.gov/41472478/)  
   <sub>QSP와 파머코메트릭스의 통합.</sub>
108. Cheng Y, Straube R, Alnaif AE, et al. Virtual Populations for Quantitative Systems Pharmacology Models. *Methods Mol Biol*. 2022;2486:129-179. [PMID 35437722](https://pubmed.ncbi.nlm.nih.gov/35437722/)  
   <sub>QSP 모델의 가상 인구(virtual population) 구축 — 본 모델의 P_success(로지스틱 변환) 접근의 방법론적 배경.</sub>
109. Coto-Segura P, Segú-Vergés C, Martorell A, et al. A quantitative systems pharmacology model for certolizumab pegol treatment in moderate-to-severe psoriasis. *Front Immunol*. 2023;14:1212981. [PMID 37809085](https://pubmed.ncbi.nlm.nih.gov/37809085/)  
   <sub>피부질환(건선)에서의 QSP 모델 사례 — 피부과 QSP의 선례.</sub>
110. Lu T, Poon V, Brooks L, et al. gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve. *CPT Pharmacometrics Syst Pharmacol*. 2024;13:341-358. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)  
   <sub>gPKPDviz — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구.</sub>

---

## 문헌이 모델 구조로 들어간 지점 (mapping)

| 모델 요소 | 근거 문헌 (섹션) |
|---|---|
| 네 개의 병인 기둥과 미세면포 저장고 (`MC`) | A · B |
| 안드로겐 → `ARS` → 피지선 질량 `SGM` · 피지율 `SER` | C |
| IGF-1 / 인슐린 / mTORC1 · FoxO1 → `LIP` | D |
| `SER` 상승 → 리놀레산 희석 `LA` → 과각화 `KER` | B · C |
| 스쿠알렌 과산화물 `SQOX` (면포유발 + NLRP3 2차 신호) | C · F |
| `CAP`/`CAB` 로지스틱 성장, 니치 = 피지, 바이오필름 보호 | E |
| 내성 분획 `RESF` — 항생제 선택압 ↑ / BPO 정화 ↓ | J |
| TLR2 → NLRP3 → IL-1β → IL-8 → 호중구 → MMP | F |
| Th17/IL-17 증폭 고리 | F |
| CRH·스트레스 → 피지 생성 | G |
| 병변 전이 사슬과 흉터·PIH 누적 | B · H |
| 국소 레티노이드/BPO/클린다마이신/답손/아젤라산 PD | I |
| 테트라사이클린의 항균 arm 과 비항균 항염 arm 분리 | J |
| COC(SHBG·LH) · 스피로놀락톤 · 클라스코테론 | K |
| 이소트레티노인 PK(식이효과·4-옥소), 피지선 자멸, 누적용량–재발 | L |
| 안전성(TG·ALT·K⁺·점막피부·최기형성·IIH) | K · L |
| 시나리오 설계와 엔드포인트 정의 | I · J · K · M |
| QSP 방법론·mrgsolve 구현 | O |

---

**총 110건.** 모든 링크는 PubMed로 연결됩니다.
교육·연구 목적의 모델 문서이며, 임상 진료 지침을 대체하지 않습니다.

