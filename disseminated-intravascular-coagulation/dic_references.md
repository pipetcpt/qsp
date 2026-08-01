# 파종성 혈관내 응고 (DIC) — 참고문헌
# Disseminated Intravascular Coagulation — Reference List

이 목록은 `dic_qsp_model.dot` / `dic_mrgsolve_model.R` 의 구조와 파라미터
근거를 섹션별로 정리한 것입니다. 각 항목의 PMID 링크는 PubMed로 연결됩니다.
모델의 각 방정식이 어느 문헌에 근거하는지는 §12 의 대응표를 참조하십시오.

---

## 1. 정의 · 진단기준 · 역학 (Definition, scoring, epidemiology)

1. Taylor FB Jr, Toh CH, Hoots WK, Wada H, Levi M. **Towards definition, clinical and laboratory criteria, and a scoring system for disseminated intravascular coagulation.** *Thromb Haemost.* 2001;86(5):1327-30. — 원본 ISTH overt-DIC 점수. 모델의 `ISTH` 출력이 그대로 구현. [PMID 11816725](https://pubmed.ncbi.nlm.nih.gov/11816725/)
2. Levi M, Ten Cate H. **Disseminated intravascular coagulation.** *N Engl J Med.* 1999;341(8):586-92. [PMID 10451465](https://pubmed.ncbi.nlm.nih.gov/10451465/)
3. Levi M, van der Poll T. **Coagulation and sepsis.** *Thromb Res.* 2017;149:38-44. [PMID 27886531](https://pubmed.ncbi.nlm.nih.gov/27886531/)
4. Iba T, Nisio MD, Levy JH, Kitamura N, Thachil J. **New criteria for sepsis-induced coagulopathy (SIC) following the revised sepsis definition.** *BMJ Open.* 2017;7(9):e017046. — SIC 점수. 섬유소원을 의도적으로 제외한 이유가 모델에서 유도됨. [PMID 28963294](https://pubmed.ncbi.nlm.nih.gov/28963294/)
5. Iba T, Levy JH, Warkentin TE, Thachil J, van der Poll T, Levi M. **Diagnosis and management of sepsis-induced coagulopathy and disseminated intravascular coagulation.** *J Thromb Haemost.* 2019;17(11):1989-94. [PMID 31410983](https://pubmed.ncbi.nlm.nih.gov/31410983/)
6. Gando S, Iba T, Eguchi Y, et al. **A multicenter, prospective validation of disseminated intravascular coagulation diagnostic criteria for critically ill patients: comparing current criteria.** *Crit Care Med.* 2006;34(3):625-31. — JAAM 기준. [PMID 16521260](https://pubmed.ncbi.nlm.nih.gov/16521260/)
7. Gando S, Levi M, Toh CH. **Disseminated intravascular coagulation.** *Nat Rev Dis Primers.* 2016;2:16037. [PMID 27250996](https://pubmed.ncbi.nlm.nih.gov/27250996/)
8. Singer M, Deutschman CS, Seymour CW, et al. **The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3).** *JAMA.* 2016;315(8):801-10. — SOFA 정의. [PMID 26903338](https://pubmed.ncbi.nlm.nih.gov/26903338/)

## 2. 조직인자 · 트롬빈 생성 (Tissue factor and thrombin generation)

9. Mackman N. **Role of tissue factor in hemostasis, thrombosis, and vascular development.** *Arterioscler Thromb Vasc Biol.* 2004;24(6):1015-22. [PMID 15117736](https://pubmed.ncbi.nlm.nih.gov/15117736/)
10. Osterud B, Bjorklid E. **Sources of tissue factor.** *Semin Thromb Hemost.* 2006;32(1):11-23. [PMID 16479458](https://pubmed.ncbi.nlm.nih.gov/16479458/)
11. Hemker HC, Giesen P, Al Dieri R, et al. **Calibrated automated thrombin generation measurement in clotting plasma.** *Pathophysiol Haemost Thromb.* 2003;33(1):4-15. — 트롬빈 생성 곡선의 정량적 규모(peak 200-400 nM)를 제공. [PMID 12853707](https://pubmed.ncbi.nlm.nih.gov/12853707/)
12. Mann KG, Butenas S, Brummel K. **The dynamics of thrombin formation.** *Arterioscler Thromb Vasc Biol.* 2003;23(1):17-25. [PMID 12524220](https://pubmed.ncbi.nlm.nih.gov/12524220/)
13. Hoffman M, Monroe DM 3rd. **A cell-based model of hemostasis.** *Thromb Haemost.* 2001;85(6):958-65. — 모델이 `PSURF`(활성화 혈소판 인지질 표면)를 별도 상태변수로 둔 근거. [PMID 11434702](https://pubmed.ncbi.nlm.nih.gov/11434702/)

## 3. 내피 · 글리코칼릭스 · 트롬보모듈린 (Endothelium, glycocalyx, thrombomodulin)

14. Esmon CT. **The protein C pathway.** *Chest.* 2003;124(3 Suppl):26S-32S. [PMID 12970121](https://pubmed.ncbi.nlm.nih.gov/12970121/)
15. Faust SN, Levin M, Harrison OB, et al. **Dysfunction of endothelial protein C activation in severe meningococcal sepsis.** *N Engl J Med.* 2001;345(6):408-16. — 패혈증에서 TM/EPCR 발현이 소실되어 protein C 활성화 자체가 실패한다는 직접 조직학적 증거. 모델의 4항 곱(THR x TM x EPCR x PC)의 근거. [PMID 11496851](https://pubmed.ncbi.nlm.nih.gov/11496851/)
16. Ait-Oufella H, Maury E, Lehoux S, Guidet B, Offenstadt G. **The endothelium: physiological functions and role in microcirculatory failure during severe sepsis.** *Intensive Care Med.* 2010;36(8):1286-98. [PMID 20443110](https://pubmed.ncbi.nlm.nih.gov/20443110/)
17. Chappell D, Jacob M, Becker BF, Hofmann-Kiefer K, Conzen P, Rehm M. **Expedition glycocalyx: a newly discovered "Great Barrier Reef".** *Anaesthesist.* 2008;57(10):959-69. [PMID 18773178](https://pubmed.ncbi.nlm.nih.gov/18773178/)
18. Ostrowski SR, Johansson PI. **Endothelial glycocalyx degradation induces endogenous heparinization in patients with severe injury and early traumatic coagulopathy.** *J Trauma Acute Care Surg.* 2012;73(1):60-6. — syndecan-1 shedding. 모델의 `GLX` 및 모세혈관 누출항. [PMID 22743373](https://pubmed.ncbi.nlm.nih.gov/22743373/)
19. Ono T, Mimuro J, Madoiwa S, et al. **Severe secondary deficiency of von Willebrand factor-cleaving protease (ADAMTS13) in patients with sepsis-induced disseminated intravascular coagulation.** *Blood.* 2006;107(2):528-34. [PMID 16189276](https://pubmed.ncbi.nlm.nih.gov/16189276/)

## 4. NET · 히스톤 · 혈전염증 (NETs, histones, thromboinflammation)

20. Xu J, Zhang X, Pelayo R, et al. **Extracellular histones are major mediators of death in sepsis.** *Nat Med.* 2009;15(11):1318-21. — 순환 히스톤 농도(10-50 µg/mL)와 세포독성. 모델의 `HIST` 및 그 하류 효과. [PMID 19855397](https://pubmed.ncbi.nlm.nih.gov/19855397/)
21. Fuchs TA, Brill A, Duerschmied D, et al. **Extracellular DNA traps promote thrombosis.** *Proc Natl Acad Sci USA.* 2010;107(36):15880-5. [PMID 20798043](https://pubmed.ncbi.nlm.nih.gov/20798043/)
22. Semeraro F, Ammollo CT, Morrissey JH, et al. **Extracellular histones promote thrombin generation through platelet-dependent mechanisms: involvement of platelet TLR2 and TLR4.** *Blood.* 2011;118(7):1952-61. — 모델에서 히스톤이 혈소판 활성화의 주된 구동항인 근거. [PMID 21673343](https://pubmed.ncbi.nlm.nih.gov/21673343/)
23. Ammollo CT, Semeraro F, Xu J, Esmon NL, Esmon CT. **Extracellular histones increase plasma thrombin generation by impairing thrombomodulin-dependent protein C activation.** *J Thromb Haemost.* 2011;9(9):1795-803. [PMID 21711444](https://pubmed.ncbi.nlm.nih.gov/21711444/)
24. Engelmann B, Massberg S. **Thrombosis as an intravascular effector of innate immunity.** *Nat Rev Immunol.* 2013;13(1):34-45. [PMID 23222502](https://pubmed.ncbi.nlm.nih.gov/23222502/)

## 5. 급성기 반응 — CLOCK 1 (Acute-phase response)

25. Gabay C, Kushner I. **Acute-phase proteins and other systemic responses to inflammation.** *N Engl J Med.* 1999;340(6):448-54. — 섬유소원은 양성, 안티트롬빈·protein C는 음성 급성기 단백. 모델 CLOCK 1의 직접적 근거. [PMID 9971870](https://pubmed.ncbi.nlm.nih.gov/9971870/)
26. Fuller GM, Zhang Z. **Transcriptional control mechanism of fibrinogen gene expression.** *Ann N Y Acad Sci.* 2001;936:469-79. — IL-6에 의한 섬유소원 합성 2-5배 증가. `APOSF = 3`의 근거. [PMID 11460504](https://pubmed.ncbi.nlm.nih.gov/11460504/)
27. Dahlbäck B. **C4b-binding protein: a forgotten factor in thrombosis and hemostasis.** *Semin Thromb Hemost.* 2011;37(4):355-61. — C4BP 상승 → 유리 protein S 감소. 모델의 `ANEGPS`. [PMID 21805442](https://pubmed.ncbi.nlm.nih.gov/21805442/)
28. Collen D, Tytgat GN, Claeys H, Piessens R. **Metabolism and distribution of fibrinogen. I. Fibrinogen turnover in physiological conditions in humans.** *Br J Haematol.* 1972;22(6):681-700. — 섬유소원 반감기 약 100 h, 일일 turnover 2-3 g. `THFIB` 및 소비율 규모의 근거. [PMID 5064500](https://pubmed.ncbi.nlm.nih.gov/5064500/)
29. Collen D, Schetz J, de Cock F, Holmer E, Verstraete M. **Metabolism of antithrombin III (heparin cofactor) in man: effects of venous thrombosis and of heparin administration.** *Eur J Clin Invest.* 1977;7(1):27-35. — AT 반감기 약 2.8일, 그리고 헤파린 투여가 AT turnover를 증가시킨다는 관찰. 모델의 `THAT` 및 헤파린-AT 소비 결합항. [PMID 65284](https://pubmed.ncbi.nlm.nih.gov/65284/)

## 6. 섬유소용해 — CLOCK 2 (Fibrinolysis)

30. Bachmann F. **Plasminogen-plasmin enzyme system.** In: *Hemostasis and Thrombosis.* — 플라스미노겐 2 µM, α2-antiplasmin 1 µM 등 기준 농도.
31. Kruithof EK, Tran-Thang C, Gudinchet A, et al. **Fibrinolysis in pregnancy: a study of plasminogen activator inhibitors.** *Blood.* 1987;69(2):460-6. [PMID 3099859](https://pubmed.ncbi.nlm.nih.gov/3099859/)
32. Mesters RM, Flörke N, Ostermann H, Kienast J. **Increase of plasminogen activator inhibitor levels predicts outcome of leukocytopenic patients with sepsis.** *Thromb Haemost.* 1996;75(6):902-7. — 패혈증에서 PAI-1 200-1000 ng/mL. CLOCK 2의 정량적 근거. [PMID 8822583](https://pubmed.ncbi.nlm.nih.gov/8822583/)
33. Hermans PW, Hibberd ML, Booy R, et al. **4G/5G promoter polymorphism in the plasminogen-activator-inhibitor-1 gene and outcome of meningococcal disease.** *Lancet.* 1999;354(9178):556-60. — PAI-1 유전형이 예후를 바꾼다 = 섬유소용해 차단이 인과적임을 시사. [PMID 10470699](https://pubmed.ncbi.nlm.nih.gov/10470699/)
34. Bajzar L, Morser J, Nesheim M. **TAFI, or plasma procarboxypeptidase B, couples the coagulation and fibrinolytic cascades through the thrombin-thrombomodulin complex.** *J Biol Chem.* 1996;271(28):16603-8. — 같은 트롬빈-TM 복합체가 protein C와 TAFI를 동시에 활성화. 모델에서 rTM의 양면성이 여기서 나옴. [PMID 8663147](https://pubmed.ncbi.nlm.nih.gov/8663147/)
35. Boffa MB, Koschinsky ML. **Curiouser and curiouser: recent advances in measurement of thrombin-activatable fibrinolysis inhibitor (TAFI) and in understanding its molecular genetics, gene regulation, and biological roles.** *Clin Biochem.* 2007;40(7-8):431-42. — 활성 TAFIa 반감기 약 8-15분. 모델의 `KDTAFI`. [PMID 17331491](https://pubmed.ncbi.nlm.nih.gov/17331491/)
36. Mosnier LO, Bouma BN. **Regulation of fibrinolysis by thrombin activatable fibrinolysis inhibitor, an unstable carboxypeptidase B that unites the pathways of coagulation and fibrinolysis.** *Arterioscler Thromb Vasc Biol.* 2006;26(11):2445-53. [PMID 16960106](https://pubmed.ncbi.nlm.nih.gov/16960106/)

## 7. 급성 전골수구성 백혈병의 응고병증 (APL coagulopathy)

37. Menell JS, Cesarman GM, Jacovina AT, McLaughlin MA, Lev EA, Hajjar KA. **Annexin II and bleeding in acute promyelocytic leukemia.** *N Engl J Med.* 1999;340(13):994-1004. — APL 아세포 표면 annexin A2가 tPA 의존적 플라스민 생성을 증폭. 모델의 `KANX`. [PMID 10099140](https://pubmed.ncbi.nlm.nih.gov/10099140/)
38. Falanga A, Marchetti M. **Anticancer treatment and thrombosis.** *Thromb Res.* 2012;129(3):353-9. [PMID 22119389](https://pubmed.ncbi.nlm.nih.gov/22119389/)
39. Stein E, McMahon B, Kwaan H, Altman JK, Frankfurt O, Tallman MS. **The coagulopathy of acute promyelocytic leukaemia revisited.** *Best Pract Res Clin Haematol.* 2009;22(1):153-63. — APL이 과섬유소용해형이라는 병태생리 정리. [PMID 19285282](https://pubmed.ncbi.nlm.nih.gov/19285282/)
40. Sanz MA, Fenaux P, Tallman MS, et al. **Management of acute promyelocytic leukemia: updated recommendations from an expert panel of the European LeukemiaNet.** *Blood.* 2019;133(15):1630-43. — 항섬유소용해제를 일상적으로 권고하지 않으며 ATRA 병용 시 혈전 위험을 경고. 모델의 시나리오 Q 결과와 일치. [PMID 30803991](https://pubmed.ncbi.nlm.nih.gov/30803991/)
41. Avvisati G, ten Cate JW, Büller HR, Mandelli F. **Tranexamic acid for control of haemorrhage in acute promyelocytic leukaemia.** *Lancet.* 1989;2(8655):122-4. — APL에서 tranexamic acid의 무작위 시험: 이익 없음. [PMID 2567895](https://pubmed.ncbi.nlm.nih.gov/2567895/)
42. Lo-Coco F, Avvisati G, Vignetti M, et al. **Retinoic acid and arsenic trioxide for acute promyelocytic leukemia.** *N Engl J Med.* 2013;369(2):111-21. — APL0406. 2년 EFS 97% vs 86%. [PMID 23841729](https://pubmed.ncbi.nlm.nih.gov/23841729/)
43. Tallman MS, Andersen JW, Schiffer CA, et al. **All-trans-retinoic acid in acute promyelocytic leukemia.** *N Engl J Med.* 1997;337(15):1021-8. [PMID 9321529](https://pubmed.ncbi.nlm.nih.gov/9321529/)
44. Breccia M, Latagliata R, Cannella L, et al. **Early hemorrhagic death before starting therapy in acute promyelocytic leukemia: association with high WBC count, late diagnosis and delayed treatment initiation.** *Haematologica.* 2010;95(5):853-4. [PMID 20015874](https://pubmed.ncbi.nlm.nih.gov/20015874/)

## 8. 항응고 치료 임상시험 (Anticoagulant trials — the arms the model must reproduce)

45. Bernard GR, Vincent JL, Laterre PF, et al. **Efficacy and safety of recombinant human activated protein C for severe sepsis (PROWESS).** *N Engl J Med.* 2001;344(10):699-709. — 28일 사망률 24.7% vs 30.8%. 모델 시나리오 I. [PMID 11236773](https://pubmed.ncbi.nlm.nih.gov/11236773/)
46. Ranieri VM, Thompson BT, Barie PS, et al. **Drotrecogin alfa (activated) in adults with septic shock (PROWESS-SHOCK).** *N Engl J Med.* 2012;366(22):2055-64. — 26.4% vs 24.2%. 약물 철회로 이어짐. 모델 시나리오 J. [PMID 22616830](https://pubmed.ncbi.nlm.nih.gov/22616830/)
47. Warren BL, Eid A, Singer P, et al. **Caring for the critically ill patient. High-dose antithrombin III in severe sepsis: a randomized controlled trial (KyberSept).** *JAMA.* 2001;286(15):1869-78. — 전체 38.9% vs 38.7%; 헤파린 비병용 아군에서만 이익; 출혈 22.0% vs 12.8%. 모델 시나리오 F/G. [PMID 11597289](https://pubmed.ncbi.nlm.nih.gov/11597289/)
48. Vincent JL, Francois B, Zabolotskikh I, et al. **Effect of a recombinant human soluble thrombomodulin on mortality in patients with sepsis-associated coagulopathy: the SCARLET randomized clinical trial.** *JAMA.* 2019;321(20):1993-2002. — 26.8% vs 29.4%, p=0.32. 모델 시나리오 H. [PMID 31104069](https://pubmed.ncbi.nlm.nih.gov/31104069/)
49. Vincent JL, Ramesh MK, Ernest D, et al. **A randomized, double-blind, placebo-controlled, phase 2b study to evaluate the safety and efficacy of recombinant human soluble thrombomodulin, ART-123, in patients with sepsis and suspected disseminated intravascular coagulation.** *Crit Care Med.* 2013;41(9):2069-79. [PMID 23979365](https://pubmed.ncbi.nlm.nih.gov/23979365/)
50. Saito H, Maruyama I, Shimazaki S, et al. **Efficacy and safety of recombinant human soluble thrombomodulin (ART-123) in disseminated intravascular coagulation: results of a phase III, randomized, double-blind clinical trial.** *J Thromb Haemost.* 2007;5(1):31-41. [PMID 17059423](https://pubmed.ncbi.nlm.nih.gov/17059423/)
51. Zarychanski R, Abou-Setta AM, Kanji S, et al. **The efficacy and safety of heparin in patients with sepsis: a systematic review and metaanalysis.** *Crit Care Med.* 2015;43(3):511-8. — 헤파린의 실제 효과 크기. 모델이 이 값을 과대추정한다는 점을 판정하는 기준. [PMID 25493972](https://pubmed.ncbi.nlm.nih.gov/25493972/)
52. Umemura Y, Yamakawa K, Ogura H, Yuhara H, Fujimi S. **Efficacy and safety of anticoagulant therapy in three specific populations with sepsis: a meta-analysis of randomized controlled trials.** *J Thromb Haemost.* 2016;14(3):518-30. — DIC 아군에서만 이익이 보인다는 계층화 결과. [PMID 26670422](https://pubmed.ncbi.nlm.nih.gov/26670422/)
53. Yamakawa K, Umemura Y, Hayakawa M, et al. **Benefit profile of anticoagulant therapy in sepsis: a nationwide multicentre registry in Japan.** *Crit Care.* 2016;20(1):229. [PMID 27472991](https://pubmed.ncbi.nlm.nih.gov/27472991/)
54. Jaimes F, De La Rosa G, Morales C, et al. **Unfractioned heparin for treatment of sepsis: a randomized clinical trial (The HETRASE Study).** *Crit Care Med.* 2009;37(4):1185-96. [PMID 19242322](https://pubmed.ncbi.nlm.nih.gov/19242322/)

## 9. 헤파린 약리 · 안티트롬빈 · 헤파린 저항성 (Heparin pharmacology)

55. Hirsh J, Anand SS, Halperin JL, Fuster V. **Guide to anticoagulant therapy: heparin.** *Circulation.* 2001;103(24):2994-3018. — 헤파린이 AT-트롬빈 반응을 약 1000-4000배 가속한다는 표준 기술, 그리고 AT를 촉매적으로 이용한다는 점. 모델이 이 항을 AT에 대해 Michaelis-Menten으로 쓴 근거. [PMID 11413093](https://pubmed.ncbi.nlm.nih.gov/11413093/)
56. Olson ST, Björk I, Sheffer R, Craig PA, Shore JD, Choay J. **Role of the antithrombin-binding pentasaccharide in heparin acceleration of antithrombin-proteinase reactions.** *J Biol Chem.* 1992;267(18):12528-38. — AT-헤파린 결합의 Kd 및 촉매 화학량론. `KMAT`의 근거. [PMID 1618758](https://pubmed.ncbi.nlm.nih.gov/1618758/)
57. Levy JH, Connors JM. **Heparin resistance — clinical perspectives and management strategies.** *N Engl J Med.* 2021;385(9):826-32. — 헤파린 저항성이 AT 결핍만의 문제가 아니며 PF4·히스톤·급성기 단백 결합이 관여한다는 정리. 모델의 `KHEPN` 항. [PMID 34437782](https://pubmed.ncbi.nlm.nih.gov/34437782/)
58. Lane DA, Denton J, Flynn AM, Thunberg L, Lindahl U. **Anticoagulant activities of heparin oligosaccharides and their neutralization by platelet factor 4.** *Biochem J.* 1984;218(3):725-32. [PMID 6721826](https://pubmed.ncbi.nlm.nih.gov/6721826/)
59. Hirsh J, Bauer KA, Donati MB, Gould M, Samama MM, Weitz JI. **Parenteral anticoagulants: ACCP evidence-based clinical practice guidelines.** *Chest.* 2008;133(6 Suppl):141S-159S. — UFH·에녹사파린·폰다파리눅스·아르가트로반의 PK 파라미터 출처. [PMID 18574264](https://pubmed.ncbi.nlm.nih.gov/18574264/)

## 10. 항섬유소용해제 (Antifibrinolytics)

60. CRASH-2 trial collaborators. **Effects of tranexamic acid on death, vascular occlusive events, and blood transfusion in trauma patients with significant haemorrhage (CRASH-2).** *Lancet.* 2010;376(9734):23-32. — 조기 투여 시 이익, 3시간 이후 투여 시 해로움. 모델의 TXA 부호 역전과 같은 구조. [PMID 20554319](https://pubmed.ncbi.nlm.nih.gov/20554319/)
61. WOMAN Trial Collaborators. **Effect of early tranexamic acid administration on mortality, hysterectomy, and other morbidities in women with post-partum haemorrhage (WOMAN).** *Lancet.* 2017;389(10084):2105-16. [PMID 28456509](https://pubmed.ncbi.nlm.nih.gov/28456509/)
62. HALT-IT Trial Collaborators. **Effects of a high-dose 24-h infusion of tranexamic acid on death and thromboembolic events in patients with acute gastrointestinal bleeding (HALT-IT).** *Lancet.* 2020;395(10241):1927-36. — 이익 없음, 정맥혈전색전증 증가. [PMID 32563378](https://pubmed.ncbi.nlm.nih.gov/32563378/)
63. Pabinger I, Fries D, Schöchl H, Streif W, Toller W. **Tranexamic acid for treatment and prophylaxis of bleeding and hyperfibrinolysis.** *Wien Klin Wochenschr.* 2017;129(9-10):303-16. — TXA의 PK 및 IC50. [PMID 28251277](https://pubmed.ncbi.nlm.nih.gov/28251277/)

## 11. 수혈 · 보충요법 · 가이드라인 (Blood products and guidelines)

64. Wada H, Thachil J, Di Nisio M, et al. **Guidance for diagnosis and treatment of disseminated intravascular coagulation from harmonization of the recommendations from three guidelines.** *J Thromb Haemost.* 2013;11(4):761-7. [PMID 23379279](https://pubmed.ncbi.nlm.nih.gov/23379279/)
65. Levi M, Toh CH, Thachil J, Watson HG. **Guidelines for the diagnosis and management of disseminated intravascular coagulation.** *Br J Haematol.* 2009;145(1):24-33. [PMID 19222477](https://pubmed.ncbi.nlm.nih.gov/19222477/)
66. Squizzato A, Hunt BJ, Kinasewitz GT, et al. **Supportive management strategies for disseminated intravascular coagulation. An international consensus.** *Thromb Haemost.* 2016;115(5):896-904. [PMID 26791023](https://pubmed.ncbi.nlm.nih.gov/26791023/)
67. Kozek-Langenecker SA, Ahmed AB, Afshari A, et al. **Management of severe perioperative bleeding: guidelines from the European Society of Anaesthesiology. First update 2016.** *Eur J Anaesthesiol.* 2017;34(6):332-95. — 섬유소원 농축제 1 g 당 약 25-30 mg/dL 상승. 모델의 `FIBINF` 용량 환산. [PMID 28459785](https://pubmed.ncbi.nlm.nih.gov/28459785/)
68. Estcourt LJ, Birchall J, Allard S, et al. **Guidelines for the use of platelet transfusions.** *Br J Haematol.* 2017;176(3):365-94. — 성인 1 단위 성분채집 혈소판 당 약 30-50 x10^9/L 상승. [PMID 28009056](https://pubmed.ncbi.nlm.nih.gov/28009056/)
69. Levy JH, Sniecinski RM, Welsby IJ, Levi M. **Antithrombin: anti-inflammatory properties and clinical applications.** *Thromb Haemost.* 2016;115(4):712-28. — AT 농축제 1 IU/kg 당 약 1.4% 상승. 모델의 `ATC_LOAD` 환산. [PMID 26676886](https://pubmed.ncbi.nlm.nih.gov/26676886/)
70. Iba T, Levy JH, Raj A, Warkentin TE. **Advance in the management of sepsis-induced coagulopathy and disseminated intravascular coagulation.** *J Clin Med.* 2019;8(5):728. [PMID 31121897](https://pubmed.ncbi.nlm.nih.gov/31121897/)

## 12. 임신·외상·악성종양 DIC (Obstetric, trauma and malignancy-associated DIC)

71. Erez O, Mastrolia SA, Thachil J. **Disseminated intravascular coagulation in pregnancy: insights in pathophysiology, diagnosis and management.** *Am J Obstet Gynecol.* 2015;213(4):452-63. [PMID 25840271](https://pubmed.ncbi.nlm.nih.gov/25840271/)
72. Brohi K, Cohen MJ, Ganter MT, et al. **Acute coagulopathy of trauma: hypoperfusion induces systemic anticoagulation and hyperfibrinolysis.** *Ann Surg.* 2007;245(5):812-8. — 외상성 응고병증이 활성 protein C 경로를 통해 DIC와 구별된다는 근거. [PMID 17457176](https://pubmed.ncbi.nlm.nih.gov/17457176/)
73. Moore HB, Moore EE, Gonzalez E, et al. **Hyperfibrinolysis, physiologic fibrinolysis, and fibrinolysis shutdown: the spectrum of postinjury fibrinolysis and relevance to antifibrinolytic therapy.** *J Trauma Acute Care Surg.* 2014;77(6):811-7. — "fibrinolysis shutdown" 개념의 정량화. 모델 CLOCK 2의 임상적 대응물. [PMID 25051384](https://pubmed.ncbi.nlm.nih.gov/25051384/)
74. Levi M. **Disseminated intravascular coagulation in cancer patients.** *Best Pract Res Clin Haematol.* 2009;22(1):129-36. [PMID 19285280](https://pubmed.ncbi.nlm.nih.gov/19285280/)

## 13. QSP 모델링 방법론 (QSP and systems modelling of coagulation)

75. Hockin MF, Jones KC, Everse SJ, Mann KG. **A model for the stoichiometric regulation of blood coagulation.** *J Biol Chem.* 2002;277(21):18322-33. — 응고 캐스케이드의 고전적 ODE 모델. [PMID 11893748](https://pubmed.ncbi.nlm.nih.gov/11893748/)
76. Chelle P, Morin C, Montmartin A, Piot M, Cournil M, Tardy-Poncet B. **Evaluation and calibration of a mechanistic evolutionary model of blood coagulation.** *PLoS One.* 2018;13(3):e0193658. [PMID 29505602](https://pubmed.ncbi.nlm.nih.gov/29505602/)
77. Wajima T, Isbister GK, Duffull SB. **A comprehensive model for the humoral coagulation network in humans.** *Clin Pharmacol Ther.* 2009;86(3):290-8. — 임상 PK/PD와 결합 가능한 응고망 모델의 표준 참조. 본 모델의 인자 반감기·기저 농도 다수가 여기와 정합. [PMID 19516255](https://pubmed.ncbi.nlm.nih.gov/19516255/)
78. Nayak S, Lee D, Patel-Hett S, et al. **Using a systems pharmacology model of the blood coagulation network to predict the effects of various therapies on biomarkers.** *CPT Pharmacometrics Syst Pharmacol.* 2015;4(7):396-405. [PMID 26312163](https://pubmed.ncbi.nlm.nih.gov/26312163/)
79. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn.* 2015;42:S84-5. — 본 저장소가 사용하는 시뮬레이션 엔진.
80. Gulati A, Faed JM, Isbister GK, Duffull SB. **Application of a systems model of coagulation to the prediction of the effect of warfarin.** *CPT Pharmacometrics Syst Pharmacol.* 2014;3:e93. [PMID 24429593](https://pubmed.ncbi.nlm.nih.gov/24429593/)

---

## 14. 방정식 ↔ 문헌 대응표 (Equation-to-source map)

모델의 각 구조적 선택이 어느 문헌에 근거하는지의 요약입니다.

| 모델 요소 | 근거 문헌 |
|---|---|
| ISTH / SIC / JAAM 점수를 상태벡터의 **출력**으로 계산 | 1, 4, 6 |
| CLOCK 1 — 섬유소원은 양성, AT·PC·PS는 음성 급성기 단백 | 25, 26, 27 |
| 섬유소원 반감기 100 h, 일일 turnover 2-3 g | 28 |
| AT 반감기 65 h, 헤파린이 AT 소비를 가속 | 29 |
| Protein C 반감기 6 h (가장 먼저 바닥에 닿는 항응고 인자) | 14, 77 |
| CLOCK 2 — PAI-1 200-1000 ng/mL, t-PA 20-30 ng/mL | 32, 33 |
| PAI-1 : t-PA 몰비가 인과적으로 예후를 바꾼다 | 33 |
| annexin A2 → APL의 t-PA 의존적 플라스민 증폭 | 37 |
| α2-antiplasmin 고갈이 전신 fibrinogenolysis의 관문 | 30, 39 |
| APC 생성 = THR x TM x EPCR x PC (4항의 곱) | 14, 15, 23 |
| 같은 트롬빈-TM 복합체가 TAFI도 활성화 | 34, 36 |
| 활성 TAFIa 반감기 약 10분 | 35 |
| ART-123은 protein C는 잘, TAFI는 잘 활성화하지 못함 | 34, 49, 50 |
| 헤파린은 AT에 대해 **촉매**이므로 AT 항이 포화 | 55, 56 |
| PF4·히스톤에 의한 비-AT 헤파린 저항성 | 57, 58 |
| 히스톤이 혈소판 활성화의 주 구동항 | 20, 22 |
| 글리코칼릭스 소실 → AT/PC의 모세혈관 누출 | 17, 18 |
| 시나리오 F/G (KyberSept) | 47 |
| 시나리오 H (SCARLET) | 48, 49, 50 |
| 시나리오 I/J (PROWESS, PROWESS-SHOCK) | 45, 46 |
| 시나리오 D/E 헤파린 (모델이 **과대추정**하는 지점) | 51, 52, 54 |
| 시나리오 N/O/Q (APL, ATRA, 항섬유소용해제) | 40, 41, 42, 43 |
| TXA의 부호 역전 (조기 vs 지연, 출혈형 vs shutdown형) | 60, 62, 73 |
| 혈액제제 용량 → 농도 환산 | 67, 68, 69 |
| 응고망 ODE 모델의 선행 구조 | 75, 76, 77, 78, 80 |

---

## 15. 이 모델이 문헌과 **어긋나는** 지점 (Where this model disagrees with the evidence)

정직하게 기록해 둡니다. 모델이 만든 예측 중 아래 두 가지는 현재의 임상근거와
맞지 않으며, 반증 가능한 가설로 읽어야 합니다.

1. **미분획 헤파린의 효과 크기.** 모델은 septic DIC 28일 사망률을 43.3% →
   28.3% (−15.0%p) 로 낮춥니다. 실제 메타분석(51, 52, 54)에서는 잘해야 몇
   %p 수준이며 그마저도 DIC 아군에 한정됩니다. 모델에 빠져 있을 가능성이
   가장 큰 기전은 (i) 이미 손상된 장기로의 헤파린 유발 출혈, (ii) DIC에서
   aPTT/anti-Xa 모니터링이 어려워 실제로는 목표 농도에 도달하지 못하는 점
   입니다.
2. **AT 농축제 + 헤파린 병용.** 모델은 병용이 단독보다 낫다고 예측하지만
   KyberSept(47)에서는 병용군에서 이익이 사라지고 출혈이 늘었습니다. 모델은
   방향(병용 시 AT의 추가 이익이 작아지고 출혈 부담이 커짐)은 재현하지만
   크기가 부족합니다.

반대로 모델이 **저자의 사전 가설을 반증한** 지점도 기록합니다.
"침착된 피브린은 재고(stock)이고 항응고제는 속도(rate)만 건드리므로 늦게
주면 소용없다"는 설명은 이 모델에서 PROWESS-SHOCK를 설명하지 못했습니다.
drotrecogin alfa 투여 시작 시각을 0/6/12/24/30/48/72/96시간으로 훑으면
28일 사망률이 37.8 / 37.4 / 36.9 / 36.1 / 35.8 / 35.9 / 37.3 / 39.1 % 로
나와, 이익이 **24-48시간에서 최대**가 되고 72시간을 넘어서야 무너집니다.
타이밍은 중요하지만 "너무 늦어서"라는 방향은 아니었습니다.
