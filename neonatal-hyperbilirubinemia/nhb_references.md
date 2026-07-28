# 신생아 고빌리루빈혈증 QSP 모델 — 참고문헌
# Neonatal Hyperbilirubinaemia (NHB) QSP Model — References

이 문서는 `nhb_qsp_model.dot`, `nhb_mrgsolve_model.R`, `nhb_shiny_app.R`,
`nhb_reference_check.py`에 들어간 **모든 구조적 가정과 파라미터의 근거**를 정리한 것이다.

**총 97편** · 모든 PMID는 PubMed E-utilities(esearch + esummary)로 **조회하여**
저자·저널·연도·제목을 그대로 가져왔다. 기억에 의존해 적은 PMID는 하나도 없으며,
서지정보는 API가 반환한 문자열을 기계적으로 삽입했다.

각 항목의 주석은 "이 논문이 모델의 **어느 파라미터 또는 어느 구조적 선택**을
지지하는가"를 적은 것이다. 모델에 쓰이지 않은 문헌은 넣지 않았다.

정량적 주장의 산출값은 모두 `nhb_reference_output.txt`(독립 python 구현의 실행
출력)에서 왔고, 본문의 A1-A14는 그 파일의 절 번호를 가리킨다.

---

## 1. 임계값과 진료 지침 — 모델의 `thr_pt` / `thr_et`는 이 곡선의 해석적 근사다 (Guidelines and thresholds)

AAP 2022 곡선은 재태주수·시간·신경독성 위험인자의 함수다. 모델은 이 구조를 `plateau(GA) x (0.42 + 0.58(1-e^{-t/36h}))`로 근사하고, 위험인자는 plateau에서 2 mg/dL를 감한다. **A1b에서 알부민 3.0 g/dL만으로 -2.65 mg/dL가 계산되어 나오는데, 이는 지침을 보여주지 않은 상태에서 결합 화학량론만으로 재현된 값이다.**

1. **Kemper AR** et al. *Clinical Practice Guideline Revision: Management of Hyperbilirubinemia in the Newborn Infant 35 or More Weeks of Gestation.* Pediatrics. 2022.
   [PMID 35927462](https://pubmed.ncbi.nlm.nih.gov/35927462/)  
   → AAP 2022 개정 지침 — 광선치료·교환수혈 임계값 곡선, 신경독성 위험인자 목록(GA<38주·알부민<3.0·동종면역 용혈·G6PD·패혈증·임상적 불안정), 가정 광선치료 허용. 모델의 `GA`, `RF`, `thr_pt`, `thr_et`의 직접 근거.

2. **Slaughter JL** et al. *Technical Report: Diagnosis and Management of Hyperbilirubinemia in the Newborn Infant 35 or More Weeks of Gestation.* Pediatrics. 2022.
   [PMID 35927519](https://pubmed.ncbi.nlm.nih.gov/35927519/)  
   → 위 지침의 기술 보고서 — 임계값 도출 근거와 IVIG·교환수혈 권고의 증거 등급. `ESCAL = ET - 2`(escalation-of-care) 정의의 출처.

3. **Bhutani VK** et al. *Predictive ability of a predischarge hour-specific serum bilirubin for subsequent significant hyperbilirubinemia in healthy term and near-term newborns.* Pediatrics. 1999.
   [PMID 9917432](https://pubmed.ncbi.nlm.nih.gov/9917432/)  
   → Bhutani 노모그램 원논문 — 퇴원 전 시간별 백분위가 이후 중증 고빌리루빈혈증을 예측. 모델이 시간축을 출생 후 절대 시각(h)으로 다루는 이유.

4. **Bhutani VK** et al. *Noninvasive measurement of total serum bilirubin in a multiracial predischarge newborn population to assess the risk of severe hyperbilirubinemia.* Pediatrics. 2000.
   [PMID 10920173](https://pubmed.ncbi.nlm.nih.gov/10920173/)  
   → 다인종 코호트에서 비침습 총빌리루빈 측정의 위험 평가 성능 — 모델의 TcB 출력(`TCB = 0.92 x C_ex`)이 스크리닝 지표로서 갖는 지위.

5. **Bhutani VK** et al. *Phototherapy to prevent severe neonatal hyperbilirubinemia in the newborn infant 35 or more weeks of gestation.* Pediatrics. 2011.
   [PMID 21949150](https://pubmed.ncbi.nlm.nih.gov/21949150/)  
   → 광선치료 기술 보고서 — 조도·노출 면적·광원 스펙트럼의 정량적 권고. 모델의 `IRRSET`, `FBSASET`, 'intensive = 30 µW/cm²/nm' 정의.

6. **Wickremasinghe AC** et al. *Neonatal Hyperbilirubinemia.* Pediatr Clin North Am. 2025.
   [PMID 40619190](https://pubmed.ncbi.nlm.nih.gov/40619190/)  
   → 신생아 고빌리루빈혈증 최신 종합 리뷰 — 모델 전체 구조의 임상 프레임.

---

## 2. 유리 빌리루빈과 알부민 결합 — 모델의 제1 주제 (Free bilirubin and albumin binding)

측정되는 값(TSB)과 손상시키는 값(Bf)은 포화 결합 등온선으로 연결된다. 모델은 고친화도 1차 부위(KA1 ≈ 4.5x10⁷ M⁻¹)와 저친화도 2차 부위(KA2 ≈ 1x10⁶ M⁻¹)의 **두 부위 등온선**을 이분법으로 푼다. 단일 부위 등온선은 B/A 몰비가 1을 넘으면 물리적으로 불가능한 Bf를 예측하므로 사용할 수 없다.

7. **Ahlfors CE** et al. *Unbound (free) bilirubin: improving the paradigm for evaluating neonatal jaundice.* Clin Chem. 2009.
   [PMID 19423734](https://pubmed.ncbi.nlm.nih.gov/19423734/)  
   → 유리 빌리루빈 패러다임 — 왜 TSB가 아니라 Bf가 위험 지표인가. 모델이 `BF`를 별도 출력으로 계산하고 BBB 유입을 Bf로만 구동하는 이유.

8. **Hulzebos CV** et al. *Bilirubin-albumin binding, bilirubin/albumin ratios, and free bilirubin levels: where do we stand?* Semin Perinatol. 2014.
   [PMID 25304058](https://pubmed.ncbi.nlm.nih.gov/25304058/)  
   → 빌리루빈-알부민 결합, B/A 비, 유리 빌리루빈의 현재 위치를 정리한 리뷰. 모델 A1b의 'B/A는 알부민에 독립이지만 친화도에는 독립이 아니다'라는 유도 결과의 배경.

9. **Hulzebos CV** et al. *The bilirubin albumin ratio in the management of hyperbilirubinemia in preterm infants to improve neurodevelopmental outcome: a randomized controlled trial--BARTrial.* PLoS One. 2014.
   [PMID 24927259](https://pubmed.ncbi.nlm.nih.gov/24927259/)  
   → BARTrial — 미숙아에서 B/A 비를 치료 결정에 추가한 무작위 시험. B/A를 임상 지표로 쓰는 근거이자 그 한계의 실증.

10. **Hegyi T** et al. *Neonatal hyperbilirubinemia and the role of unbound bilirubin.* J Matern Fetal Neonatal Med. 2022.
   [PMID 34957902](https://pubmed.ncbi.nlm.nih.gov/34957902/)  
   → 비결합 빌리루빈의 역할 — Bf 측정과 임상 해석. 모델의 Bf 손상 임계 30-35 nM 설정 근거.

11. **Wennberg RP** et al. *Toward understanding kernicterus: a challenge to improve the management of jaundiced newborns.* Pediatrics. 2006.
   [PMID 16452368](https://pubmed.ncbi.nlm.nih.gov/16452368/)  
   → 핵황달 이해를 위한 관점 — TSB만으로 위험을 규정할 수 없음. 모델 제1 주제의 임상적 진술.

12. **Amin SB** et al. *Effect of ibuprofen on bilirubin-albumin binding affinity in premature infants.* J Perinat Med. 2011.
   [PMID 20954849](https://pubmed.ncbi.nlm.nih.gov/20954849/)  
   → 이부프로펜이 미숙아의 빌리루빈-알부민 결합 친화도를 낮춘다 — 모델 `FDISP` 파라미터의 직접 근거.

13. **Amin SB** et al. *Intravenous lipid and bilirubin-albumin binding variables in premature infants.* Pediatrics. 2009.
   [PMID 19564302](https://pubmed.ncbi.nlm.nih.gov/19564302/)  
   → 정맥 지질(유리지방산)이 결합 변수를 변화시킴 — `FDISP`에 유리지방산 경로가 포함되는 이유.

14. **Kleinfeld A** et al. *Displacement of Bilirubin From Albumin by Omegaven and Intralipid.* J Pediatr Surg. 2025.
   [PMID 40946844](https://pubmed.ncbi.nlm.nih.gov/40946844/)  
   → Omegaven/Intralipid에 의한 알부민으로부터의 빌리루빈 치환 — `FDISP`의 최신 정량 근거.

15. **Martin E** et al. *Ceftriaxone--bilirubin-albumin interactions in the neonate: an in vivo study.* Eur J Pediatr. 1993.
   [PMID 8335024](https://pubmed.ncbi.nlm.nih.gov/8335024/)  
   → 세프트리악손-빌리루빈-알부민 상호작용의 생체 내 연구 — `FDISP = 0.70` 시나리오.

16. **Gulian JM** et al. *Bilirubin displacement by ceftriaxone in neonates: evaluation by determination of 'free' bilirubin and erythrocyte-bound bilirubin.* J Antimicrob Chemother. 1987.
   [PMID 3610909](https://pubmed.ncbi.nlm.nih.gov/3610909/)  
   → 세프트리악손에 의한 치환을 '유리' 빌리루빈과 적혈구 결합 빌리루빈으로 측정 — 치환이 Bf를 올린다는 직접 증거.

17. **Wadsworth SJ** et al. *In vitro displacement of bilirubin by antibiotics and 2-hydroxybenzoylglycine in newborns.* Antimicrob Agents Chemother. 1988.
   [PMID 3190184](https://pubmed.ncbi.nlm.nih.gov/3190184/)  
   → 신생아 혈청에서 항생제에 의한 빌리루빈 치환의 시험관 정량 — 약물별 치환 강도의 순서.

18. **Amin SB** et al. *Unbound Bilirubin and Auditory Neuropathy Spectrum Disorder in Late Preterm and Term Infants with Severe Jaundice.* J Pediatr. 2016.
   [PMID 26952116](https://pubmed.ncbi.nlm.nih.gov/26952116/)  
   → 중증 황달 만삭·후기 조산아에서 비결합 빌리루빈이 청각신경병증(ANSD)과 연관 — 모델 `BBRTHRA`(ABR 임계)를 Bf 기반으로 둔 근거.

19. **Ahlfors CE** et al. *Unbound bilirubin predicts abnormal automated auditory brainstem response in a diverse newborn population.* J Perinatol. 2009.
   [PMID 19242487](https://pubmed.ncbi.nlm.nih.gov/19242487/)  
   → 비결합 빌리루빈이 자동 ABR 이상을 예측(TSB보다 우수) — Bf → ABR 경로의 정량 근거.

---

## 3. 빌리루빈 생성과 헴 대사, 그리고 ETCOc — 모델의 제2 주제 (Production, haem catabolism, end-tidal CO)

TSB는 두 플럭스의 차이의 적분이므로 그 자체로는 원인을 식별하지 못한다. HO-1이 빌리루빈 1분자당 CO를 정확히 1분자 방출하므로 COHb/ETCOc는 **생성 플럭스의 직접 판독값**이다. 모델은 `dCOHB/dt = KCO x PROD/W - KCOOUT x COHB`로 이를 구현한다.

20. **Stevenson DK** et al. *Bilirubin production in healthy term infants as measured by carbon monoxide in breath.* Clin Chem. 1994.
   [PMID 7923775](https://pubmed.ncbi.nlm.nih.gov/7923775/)  
   → 건강한 만삭아의 빌리루빈 생성을 호기 CO로 측정 — 모델의 생성량 8.5 mg/kg/day(성인의 약 2배) 및 `BILPERHB = 34 mg/g`, `FEARLY = 0.25`의 근거.

21. **Stevenson DK** et al. *Prediction of hyperbilirubinemia in near-term and term infants.* J Perinatol. 2001.
   [PMID 11803421](https://pubmed.ncbi.nlm.nih.gov/11803421/)  
   → 근접만삭·만삭아의 고빌리루빈혈증 예측 — 생성 지표와 TSB 궤적의 관계.

22. **Stevenson DK** et al. *Increased Carbon Monoxide Washout Rates in Newborn Infants.* Neonatology. 2020.
   [PMID 31634890](https://pubmed.ncbi.nlm.nih.gov/31634890/)  
   → 신생아의 CO 세척(washout) 속도 증가 — 모델 `KCOOUT`(t½ ≈ 5 h)의 근거.

23. **Tidmarsh GF** et al. *End-tidal carbon monoxide and hemolysis.* J Perinatol. 2014.
   [PMID 24743136](https://pubmed.ncbi.nlm.nih.gov/24743136/)  
   → 호기말 CO와 용혈 — ETCOc가 용혈 유무를 가르는 임상 지표라는 근거. 모델 A3의 핵심 주장.

24. **Bhutani VK** et al. *Identification of risk for neonatal haemolysis.* Acta Paediatr. 2018.
   [PMID 29532503](https://pubmed.ncbi.nlm.nih.gov/29532503/)  
   → 신생아 용혈 위험의 식별 — 생성 병변과 제거 병변을 임상에서 구별하는 방법.

25. **Herschel M** et al. *Isoimmunization is unlikely to be the cause of hemolysis in ABO-incompatible but direct antiglobulin test-negative neonates.* Pediatrics. 2002.
   [PMID 12093957](https://pubmed.ncbi.nlm.nih.gov/12093957/)  
   → ABO 부적합이지만 DAT 음성인 신생아에서는 동종면역이 용혈의 원인이 아니다 — 모델이 용혈 부하를 `ABMAT`(항체 부하)라는 연속 변수로 두고 혈액형 부적합 자체로 두지 않은 이유.

26. **Wong RJ** et al. *In vitro inhibition of heme oxygenase isoenzymes by metalloporphyrins.* J Perinatol. 2011.
   [PMID 21448202](https://pubmed.ncbi.nlm.nih.gov/21448202/)  
   → 메탈로포르피린의 HO 이성질체 시험관 억제 — 모델 `KISNMP`(경쟁적 억제 상수)와 스타노포르핀의 작용 지점.

27. **Fujioka K** et al. *Inhibition of heme oxygenase activity using a microparticle formulation of zinc protoporphyrin in an acute hemolytic newborn mouse model.* Pediatr Res. 2016.
   [PMID 26488552](https://pubmed.ncbi.nlm.nih.gov/26488552/)  
   → 급성 용혈 신생아 마우스에서 아연 프로토포르피린 미립자 제형에 의한 HO 활성 억제 — 생성 차단이 빌리루빈을 낮춘다는 인과 증명.

---

## 4. UGT1A1 성숙·유전형·간 흡수 — 제거 플럭스 (Conjugation, ontogeny, genotype)

모델은 `UGT_target = GENO x ont(t) x (1 + E_pb C/(EC50+C)) + TGX`로 쓴다. 이 곱셈 구조에서 **GENO = 0(CN-I)이면 유도항이 함께 소멸**하며, 이것이 페노바르비탈이 CN-II에서는 작동하고 CN-I에서는 작동하지 않는 이유를 별도 규칙 없이 만들어낸다.

28. **Nie YL** et al. *Histone Modifications Regulate the Developmental Expression of Human Hepatic UDP-Glucuronosyltransferase 1A1.* Drug Metab Dispos. 2017.
   [PMID 29025858](https://pubmed.ncbi.nlm.nih.gov/29025858/)  
   → 인간 간 UGT1A1의 발달적 발현이 히스톤 수식으로 조절됨 — `TAUUGT` 출생 후 성숙 시간상수의 분자적 근거.

29. **Nie YL** et al. *Hepatic expression of transcription factors affecting developmental regulation of UGT1A1 in the Han Chinese population.* Eur J Clin Pharmacol. 2017.
   [PMID 27704169](https://pubmed.ncbi.nlm.nih.gov/27704169/)  
   → 한족 집단에서 UGT1A1 발달 조절 전사인자의 간 발현 — 성숙 곡선의 인구집단 변이.

30. **Travan L** et al. *Severe neonatal hyperbilirubinemia and UGT1A1 promoter polymorphism.* J Pediatr. 2014.
   [PMID 24726540](https://pubmed.ncbi.nlm.nih.gov/24726540/)  
   → 중증 신생아 고빌리루빈혈증과 UGT1A1 프로모터 다형성 — `GENO`(Gilbert 0.35)의 임상 근거.

31. **Žaja O** et al. *Correlation of UGT1A1 TATA-box polymorphism and jaundice in breastfed newborns-early presentation of Gilbert's syndrome.* J Matern Fetal Neonatal Med. 2014.
   [PMID 23981182](https://pubmed.ncbi.nlm.nih.gov/23981182/)  
   → UGT1A1 TATA-box 다형성과 모유수유 신생아 황달의 상관 — 모델 A9의 유전형 x 수유 상호작용.

32. **Yamamoto A** et al. *Gly71Arg mutation of the bilirubin UDP-glucuronosyltransferase 1A1 gene is associated with neonatal hyperbilirubinemia in the Japanese population.* Kobe J Med Sci. 2002.
   [PMID 12502904](https://pubmed.ncbi.nlm.nih.gov/12502904/)  
   → 일본인 집단에서 Gly71Arg(UGT1A1*6)와 신생아 고빌리루빈혈증의 연관 — `GENO = 0.40`(동아시아 변이).

33. **Prachukthum S** et al. *Association between UGT 1A1 Gly71Arg (G71R) polymorphism and neonatal hyperbilirubinemia.* J Med Assoc Thai. 2012.
   [PMID 23964438](https://pubmed.ncbi.nlm.nih.gov/23964438/)  
   → 태국인에서 UGT1A1 G71R과 신생아 고빌리루빈혈증 — 동아시아·동남아 변이의 재현.

34. **Talebi H** et al. *Association of SLCO1B1 genetic variants with neonatal hyperbilirubinemia: a consolidated analysis of 36 studies.* BMC Pediatr. 2025.
   [PMID 40155882](https://pubmed.ncbi.nlm.nih.gov/40155882/)  
   → 36개 연구를 통합한 SLCO1B1 변이와 신생아 고빌리루빈혈증 — 모델 `FOATP`(간 흡수 유전형 계수).

35. **Fan J** et al. *Associations between UGT1A1, SLCO1B1, SLCO1B3, BLVRA and HMOX1 polymorphisms and susceptibility to neonatal severe hyperbilirubinemia in Chinese Han population.* BMC Pediatr. 2024.
   [PMID 38279097](https://pubmed.ncbi.nlm.nih.gov/38279097/)  
   → UGT1A1·SLCO1B1·SLCO1B3·BLVRA·HMOX1 다형성과 중증 고빌리루빈혈증 — 모델이 흡수·포합·생성 각 단계에 유전 계수를 둔 근거.

36. **Wang WF** et al. *Variants in UGT1A1 and SLCO1B1 increase the risk of neonatal hyperbilirubinemia: a case-control study in subtropical China.* Front Pediatr. 2026.
   [PMID 42422451](https://pubmed.ncbi.nlm.nih.gov/42422451/)  
   → 아열대 중국 코호트에서 UGT1A1과 SLCO1B1 변이의 위험 증가 — 흡수와 포합이 독립적 위험 축이라는 근거.

---

## 5. Crigler-Najjar와 효소 대체 — 모델의 제4 주제 (Crigler-Najjar, gene transfer)

루미루빈은 포합을 필요로 하지 않으므로 광선치료는 UGT1A1을 완전히 우회하는 유일한 경로다. 모델은 성인 UGT1A1 활성의 약 10 %만으로 CN-I이 CN-II 표현형으로 바뀌며 TSB가 램프 없이 ~9 mg/dL에 머무는 것을 계산한다(A11c).

37. **D'Antiga L** et al. *Gene Therapy in Patients with the Crigler-Najjar Syndrome.* N Engl J Med. 2023.
   [PMID 37585628](https://pubmed.ncbi.nlm.nih.gov/37585628/)  
   → GNT0003 최초 인체 유전자 치료 — CN-I 환자에서 AAV8-hUGT1A1 투여 후 광선치료 없이 빌리루빈이 조절됨. 모델 `TGXMAX = 0.10`과 `KTGON`(t½ 10일)의 직접 근거.

38. **Collaud F** et al. *Preclinical Development of an AAV8-hUGT1A1 Vector for the Treatment of Crigler-Najjar Syndrome.* Mol Ther Methods Clin Dev. 2019.
   [PMID 30705921](https://pubmed.ncbi.nlm.nih.gov/30705921/)  
   → AAV8-hUGT1A1 벡터의 전임상 개발 — 발현 수준과 용량 반응.

39. **Aronson SJ** et al. *Prevalence and Relevance of Pre-Existing Anti-Adeno-Associated Virus Immunity in the Context of Gene Therapy for Crigler-Najjar Syndrome.* Hum Gene Ther. 2019.
   [PMID 31502485](https://pubmed.ncbi.nlm.nih.gov/31502485/)  
   → CN 유전자 치료에서 기존 항-AAV 면역의 유병률과 의미 — 모델이 다루지 않는 실패 양식(중화항체)을 명시하기 위한 인용.

40. **Agati G** et al. *Bilirubin photoisomerization products in serum and urine from a Crigler-Najjar type I patient treated by phototherapy.* J Photochem Photobiol B. 1998.
   [PMID 10093917](https://pubmed.ncbi.nlm.nih.gov/10093917/)  
   → 광선치료를 받은 CN-I 환자의 혈청·소변 광이성질체 — **인간에서 광화학 산물이 실제로 배출된다는 직접 증거이자, 광선치료가 포합과 무관하게 작동한다는 결정적 관찰.** 모델 제3·제4 주제의 실험적 기반.

41. **Jesus Sá C** et al. *Between Crigler-Najjar Syndrome Type II and Gilbert Syndrome: Expanding the Spectrum of Uridine Diphosphate Glucuronosyltransferase 1A1 (UGT1A1)-Related Hyperbilirubinemia.* Cureus. 2026.
   [PMID 41658774](https://pubmed.ncbi.nlm.nih.gov/41658774/)  
   → CN-II와 Gilbert 사이의 스펙트럼 — `GENO`가 연속적 활성 스케일이라는 근거(CN-I 0 → CN-II 0.05 → Gilbert 0.35 → 야생형 1.0).

---

## 6. 광화학과 광선치료 용량-반응 — 모델의 제3 주제 (Photochemistry)

광선치료는 두 반응이다. 빠르고 **가역적인** (4Z,15Z) ⇄ (4Z,15E) 배위 이성질화가 광정상상태(모델에서 측정 TSB의 20.6 %)에 도달하고, E-이성질체의 **비가역적** 고리화가 루미루빈을 만들어 실제 배출을 담당한다. 두 광이성질체는 통상 분석법에서 '빌리루빈'으로 측정되므로 보고된 TSB의 26 %가 광이성질체이며, 램프를 끄면 TSB가 먼저 **하강**(루미루빈 세척, t½ 1 h)한 뒤 열적 역이성질화로 재상승한다.

42. **Tan KL** *The nature of the dose-response relationship of phototherapy for neonatal hyperbilirubinemia.* J Pediatr. 1977.
   [PMID 839340](https://pubmed.ncbi.nlm.nih.gov/839340/)  
   → 광선치료의 용량-반응 관계의 본질 — 조도를 올릴 때의 수확 체감. 모델 `I50 = 45 µW/cm²/nm`(광학적 포화)의 고전적 근거.

43. **Tan KL** *The pattern of bilirubin response to phototherapy for neonatal hyperbilirubinaemia.* Pediatr Res. 1982.
   [PMID 7110789](https://pubmed.ncbi.nlm.nih.gov/7110789/)  
   → 광선치료에 대한 빌리루빈 반응의 패턴 — 24시간 감소율의 크기(모델: 집중 광선치료에서 39 %).

44. **Pratesi R** et al. *Phototherapy for neonatal hyperbilirubinemia.* Photodermatol. 1989.
   [PMID 2700092](https://pubmed.ncbi.nlm.nih.gov/2700092/)  
   → 신생아 고빌리루빈혈증의 광선치료 — 광원 스펙트럼, 침투 깊이, 광화학 경로 정리. 모델이 광자 전달을 `조도 x 노출면적 x 광학 접근도`로 분해한 근거.

45. **Borden AR** et al. *Variation in the Phototherapy Practices and Irradiance of Devices in a Major Metropolitan Area.* Neonatology. 2018.
   [PMID 29393277](https://pubmed.ncbi.nlm.nih.gov/29393277/)  
   → 대도시권 장비의 광선치료 실무와 실제 조도의 변이 — 처방된 조도와 전달된 조도가 다르다는 점. 모델의 `IRRSET`이 전달 조도임을 명시하는 이유.

46. **Chang PW** et al. *A Clinical Prediction Rule for Rebound Hyperbilirubinemia Following Inpatient Phototherapy.* Pediatrics. 2017.
   [PMID 28196932](https://pubmed.ncbi.nlm.nih.gov/28196932/)  
   → 광선치료 중단 후 재상승의 임상 예측 규칙 — 모델의 재상승(A5: +2.2 mg/dL)이 일부는 광화학적 장부(E→Z 역전환)라는 해석의 대조 기준.

47. **Slusher TM** et al. *A Randomized Trial of Phototherapy with Filtered Sunlight in African Neonates.* N Engl J Med. 2015.
   [PMID 26376136](https://pubmed.ncbi.nlm.nih.gov/26376136/)  
   → 여과 태양광 광선치료 무작위 시험 — 광자 전달이 본질이고 전원이 아님을 보여주는 자원 제한 환경의 증거.

48. **Slusher TM** et al. *Filtered sunlight versus intensive electric powered phototherapy in moderate-to-severe neonatal hyperbilirubinaemia: a randomised controlled non-inferiority trial.* Lancet Glob Health. 2018.
   [PMID 30170894](https://pubmed.ncbi.nlm.nih.gov/30170894/)  
   → 여과 태양광 대 집중 전기 광선치료 비열등성 시험(중등도-중증) — 위와 같음.

49. **Magee S** et al. *Multi-directional phototherapy vs. unidirectional phototherapy from below for severe neonatal jaundice: a randomized pilot trial in home phototherapy.* Eur J Pediatr. 2026.
   [PMID 41975040](https://pubmed.ncbi.nlm.nih.gov/41975040/)  
   → 다방향 대 하방 단일방향 광선치료 무작위 시험 — **노출 체표면적이 조도보다 중요하다는 모델 A4의 결론에 대응하는 임상 시험.**

50. **Cordero N** et al. *Transcutaneous bilirubin measurements in preterm infants: the impact of race, age, and phototherapy.* J Perinatol. 2026.
   [PMID 41490938](https://pubmed.ncbi.nlm.nih.gov/41490938/)  
   → 미숙아에서 인종·연령·광선치료가 TcB에 미치는 영향 — 모델이 TcB를 피부 구획(`BEX`)에서 계산하고 램프 아래에서 혈청과 해리되도록 만든 이유.

---

## 7. 장간 순환 — 신생아의 두 번째 간, 그러나 반대 방향 (Enterohepatic shunt)

포합 빌리루빈은 신생아에서 출구가 아니라 **재순환**이다. β-글루쿠로니다제가 이를 탈포합하고 공장이 재흡수해 문맥을 통해 간으로 되돌린다. 따라서 유효 제거율은 `CL_conj x (1 - f_reabsorb)`이며, 수유 지원은 약효와 같은 단위를 갖는다.

51. **Kreamer BL** et al. *A novel inhibitor of beta-glucuronidase: L-aspartic acid.* Pediatr Res. 2001.
   [PMID 11568288](https://pubmed.ncbi.nlm.nih.gov/11568288/)  
   → β-글루쿠로니다제 억제제(L-아스파르트산) — 탈포합 단계가 약리학적 표적이 될 수 있다는 증명. 모델 `KBGLUC` 항의 개입 가능성.

52. **Alonso EM** et al. *Enterohepatic circulation of nonconjugated bilirubin in rats fed with human milk.* J Pediatr. 1991.
   [PMID 1999786](https://pubmed.ncbi.nlm.nih.gov/1999786/)  
   → 모유를 먹인 쥐에서 비포합 빌리루빈의 장간 순환 — `BGABM`(모유에 의한 β-글루쿠로니다제 증가) 항의 근거.

53. **Yau KI** et al. *Factors affecting the severity of neonatal jaundice of unknown etiology: the role of enterohepatic circulation.* Zhonghua Min Guo Xiao Er Ke Yi Xue Hui Za Zhi. 1992.
   [PMID 1626448](https://pubmed.ncbi.nlm.nih.gov/1626448/)  
   → 원인 불명 신생아 황달의 중증도에 장간 순환이 기여 — 재순환 분율 ~50 %의 근거.

54. **Vítek L** et al. *The impact of intestinal microflora on serum bilirubin levels.* J Hepatol. 2005.
   [PMID 15664250](https://pubmed.ncbi.nlm.nih.gov/15664250/)  
   → 장내 미생물총이 혈청 빌리루빈 수준에 미치는 영향 — 무균 장에서 유로빌리노이드 전환이 없다는 점. 모델 `TAUBGA`(flora 성숙 21일).

55. **Vítek L** et al. *Gut microbiota and bilirubin metabolism: unveiling new pathways in health and disease.* Trends Mol Med. 2025.
   [PMID 39757046](https://pubmed.ncbi.nlm.nih.gov/39757046/)  
   → 장내 미생물총과 빌리루빈 대사의 최신 리뷰 — 위 경로의 현대적 정리.

56. **Vítek L** et al. *Identification of bilirubin reduction products formed by Clostridium perfringens isolated from human neonatal fecal flora.* J Chromatogr B Analyt Technol Biomed Life Sci. 2006.
   [PMID 16504607](https://pubmed.ncbi.nlm.nih.gov/16504607/)  
   → 신생아 대변 flora의 Clostridium perfringens에 의한 빌리루빈 환원 산물 — flora 확립이 재순환을 끊는 기전.

57. **Abdel-Aziz Ali SM** et al. *Efficacy of oral agar in management of indirect hyperbilirubinemia in full-term neonates.* J Matern Fetal Neonatal Med. 2022.
   [PMID 32192396](https://pubmed.ncbi.nlm.nih.gov/32192396/)  
   → 만삭아 간접 고빌리루빈혈증에서 경구 아가의 효능 — 모델 `GBIND`/`KB50`(장내 결합제) 항의 임상 근거.

58. **Lazarus G** et al. *Role of ursodeoxycholic acid in neonatal indirect hyperbilirubinemia: a systematic review and meta-analysis of randomized controlled trials.* Ital J Pediatr. 2022.
   [PMID 36253867](https://pubmed.ncbi.nlm.nih.gov/36253867/)  
   → 신생아 간접 고빌리루빈혈증에서 UDCA의 역할: 무작위 시험 메타분석 — `EMAXU`, `EC50UDCA`.

59. **Tabrizi M** et al. *Ursodeoxycholic Acid in the Management of Prolonged Neonatal Hyperbilirubinemia: A Randomized Controlled Clinical Trial.* Iran J Med Sci. 2026.
   [PMID 42100234](https://pubmed.ncbi.nlm.nih.gov/42100234/)  
   → 지속성 신생아 고빌리루빈혈증에서 UDCA 무작위 시험 — 위와 같음, 지속성 황달 대상.

60. **de Oliveira HM** et al. *Zinc sulfate on neonatal hyperbilirubinemia: an updated systematic review and meta-analysis.* Eur J Pediatr. 2024.
   [PMID 39671002](https://pubmed.ncbi.nlm.nih.gov/39671002/)  
   → 아연 설페이트의 갱신된 체계적 고찰·메타분석 — 장내 침전을 통한 재흡수 차단(모델에서 `GBIND` 계열).

61. **Gholitabar M** et al. *Clofibrate in combination with phototherapy for unconjugated neonatal hyperbilirubinaemia.* Cochrane Database Syst Rev. 2012.
   [PMID 23235669](https://pubmed.ncbi.nlm.nih.gov/23235669/)  
   → 클로피브레이트 + 광선치료 Cochrane 고찰 — PPARα를 통한 UGT1A1 유도(모델 `CAR` 경로의 대체 작용제).

62. **Fei Q** et al. *Phototherapy Alone or Combined with Adjuvant Drugs for Neonatal Hyperbilirubinemia: A Systematic Review and Network Meta-Analysis.* Children (Basel). 2026.
   [PMID 42073150](https://pubmed.ncbi.nlm.nih.gov/42073150/)  
   → 광선치료 단독 대 보조 약물 병용의 네트워크 메타분석 — 모델의 보조요법 순위(A10)와 비교할 최신 정량 기준.

63. **Tavares LC** et al. *Oral adjuvants supplemented to phototherapy as a therapeutic strategy for neonatal hyperbilirubinemia: a systematic review.* J Perinat Med. 2026.
   [PMID 41811691](https://pubmed.ncbi.nlm.nih.gov/41811691/)  
   → 광선치료 보조 경구 약물의 체계적 고찰 — 위와 같음.

---

## 8. 교환수혈과 IVIG — 부하·결합능·생성에 동시에 작용 (Exchange transfusion, IVIG)

모델에서 교환수혈은 순간 사건이 아니라 2-4시간에 걸친 실제 시술이며, 제거량은 `∫ Q_ET x C_p dt`의 적분으로 **계산되어 나온다**. 공혈 혈장이 알부민을 재설정하고 모체 항체를 제거하므로 Bf의 하강폭(82 %)이 TSB의 하강폭(40 %)보다 크다.

64. **Zwiers C** et al. *Immunoglobulin for alloimmune hemolytic disease in neonates.* Cochrane Database Syst Rev. 2018.
   [PMID 29551014](https://pubmed.ncbi.nlm.nih.gov/29551014/)  
   → 신생아 동종면역 용혈질환에서 면역글로불린 Cochrane 고찰 — `IMAXIVIG = 0.65`, `IC50IVIG`의 근거이자 효과 크기의 불확실성 명시.

65. **Huang M** et al. *Efficacy and safety of different doses of intravenous immunoglobulin combined with phototherapy in neonatal hemolytic disease: a systematic review and meta-analysis.* Eur J Pediatr. 2026.
   [PMID 42458120](https://pubmed.ncbi.nlm.nih.gov/42458120/)  
   → IVIG 용량별 효능·안전성 메타분석 — 1 g/kg 용량 선택의 근거.

66. **van Klink JM** et al. *Immunoglobulins in Neonates with Rhesus Hemolytic Disease of the Fetus and Newborn: Long-Term Outcome in a Randomized Trial.* Fetal Diagn Ther. 2016.
   [PMID 26159803](https://pubmed.ncbi.nlm.nih.gov/26159803/)  
   → Rh 용혈질환 신생아에서 면역글로불린의 장기 결과(무작위 시험) — IVIG가 교환수혈 필요를 줄이는지에 대한 상반 증거. 모델이 IVIG 효과를 Emax 0.65로 제한한 이유.

67. **Legler TJ** *RhIg for the prevention Rh immunization and IVIg for the treatment of affected neonates.* Transfus Apher Sci. 2020.
   [PMID 33004277](https://pubmed.ncbi.nlm.nih.gov/33004277/)  
   → 예방적 항-D(RhIg)와 치료적 IVIg — 지도에서 RhIG를 '1차 예방'으로 둔 근거.

68. **Dash N** et al. *Pre exchange Albumin Administration in Neonates with Hyperbilirubinemia: A Randomized Controlled Trial.* Indian Pediatr. 2015.
   [PMID 26519710](https://pubmed.ncbi.nlm.nih.gov/26519710/)  
   → 교환수혈 전 알부민 투여 무작위 시험 — **결합능을 올리는 것이 부하를 내리는 것과 별개의 치료 축이라는 임상 증명.** 모델 `ALBINF`/`ALBDON` 항.

69. **Shahian M** et al. *Effect of albumin administration prior to exchange transfusion in term neonates with hyperbilirubinemia--a randomized controlled trial.* Indian Pediatr. 2010.
   [PMID 19578230](https://pubmed.ncbi.nlm.nih.gov/19578230/)  
   → 교환수혈 전 알부민 투여의 무작위 시험(만삭아) — 위와 같음.

70. **Hemmati F** et al. *Exchange Transfusion Trends and Risk Factors for Extreme Neonatal Hyperbilirubinemia over 10 Years in Shiraz, Iran.* Iran J Med Sci. 2024.
   [PMID 38952637](https://pubmed.ncbi.nlm.nih.gov/38952637/)  
   → 10년간 극단적 고빌리루빈혈증의 교환수혈 추이와 위험인자 — 실제 교환수혈 대상군의 특성.

71. **Sarı EE** et al. *The '72-Hour Divergence' in severe hyperbilirubinemia: sepsis is associated with altered renal and hematologic adaptation in a 10-year exchange transfusion cohort.* BMC Pediatr. 2026.
   [PMID 41776517](https://pubmed.ncbi.nlm.nih.gov/41776517/)  
   → 중증 고빌리루빈혈증의 '72시간 분기' — 10년 교환수혈 코호트에서 패혈증이 신장·혈액학적 적응을 변화시킴. 모델 `FBBB`(패혈증에 의한 장벽 투과성 증가)와 `FACID`의 임상 대응.

---

## 9. 헴 산화효소 억제 — 생성 차단이라는 별개의 치료 축 (Production blockade)

스타노포르핀은 HO-1의 경쟁적 억제제로, 모델에서 `F_SNMP = 1/(1 + C/KISNMP)`로 **생성 항에만** 작용한다. 이것이 ETCOc로 확인되는 생성 병변에서만 의미가 있는 이유이며, 제거 병변(CN, Gilbert)에서는 무의미하다.

72. **Rosenfeld WN** et al. *Stannsoporfin with phototherapy to treat hyperbilirubinemia in newborn hemolytic disease.* J Perinatol. 2022.
   [PMID 34635771](https://pubmed.ncbi.nlm.nih.gov/34635771/)  
   → 신생아 용혈질환에서 스타노포르핀 + 광선치료 임상시험 — 단회 IM 4.5 mg/kg 용량과 효과 크기(모델 `KASNMP`, `VSNMP`, `KESNMP`, `KISNMP`).

73. **Bhutani VK** et al. *Clinical trial of tin mesoporphyrin to prevent neonatal hyperbilirubinemia.* J Perinatol. 2016.
   [PMID 26938918](https://pubmed.ncbi.nlm.nih.gov/26938918/)  
   → 주석 메소포르피린의 신생아 고빌리루빈혈증 예방 임상시험 — 생성 차단 40-75 %의 근거.

74. **Poudel P** et al. *Efficacy and Safety Concerns with Sn-Mesoporphyrin as an Adjunct Therapy in Neonatal Hyperbilirubinemia: A Literature Review.* Int J Pediatr. 2022.
   [PMID 35898803](https://pubmed.ncbi.nlm.nih.gov/35898803/)  
   → Sn-메소포르피린 보조요법의 효능·안전성 문헌 고찰 — 미승인 상태와 안전성 논점의 명시.

---

## 10. 신경독성 — 무엇이 장벽을 넘고 무엇이 손상을 남기는가 (Neurotoxicity)

모델은 `dBBR/dt = KINBBB x FBBB x Bf - KOUTBBB x BBR`로 **유리 빌리루빈만** 장벽을 넘게 하고, 손상은 임계 초과분에 비례해 누적시킨다. ABR 결손은 가역, 누적 손상은 비가역으로 분리했다.

75. **Watchko JF** et al. *Bilirubin-induced neurologic damage--mechanisms and management approaches.* N Engl J Med. 2013.
   [PMID 24256380](https://pubmed.ncbi.nlm.nih.gov/24256380/)  
   → 빌리루빈에 의한 신경 손상 — 기전과 관리(권위 있는 리뷰). 모델 제7 클러스터 전체의 구조적 출처.

76. **Ostrow JD** et al. *Molecular basis of bilirubin-induced neurotoxicity.* Trends Mol Med. 2004.
   [PMID 15102359](https://pubmed.ncbi.nlm.nih.gov/15102359/)  
   → 빌리루빈 신경독성의 분자적 기초 — 미토콘드리아 이탈결합, 산화 스트레스, 세포사. 모델 `KINJ` 경로의 기전.

77. **Watchko JF** et al. *Brain bilirubin content is increased in P-glycoprotein-deficient transgenic null mutant mice.* Pediatr Res. 1998.
   [PMID 9803459](https://pubmed.ncbi.nlm.nih.gov/9803459/)  
   → P-당단백 결핍 마우스에서 뇌 빌리루빈 함량 증가 — **모델 `KOUTBBB`에 P-gp 유출이 포함되는 직접 증거.**

78. **Hansen TW** et al. *Rates of bilirubin clearance from rat brain regions.* Biol Neonate. 1995.
   [PMID 8534773](https://pubmed.ncbi.nlm.nih.gov/8534773/)  
   → 쥐 뇌 영역별 빌리루빈 제거 속도 — `KOUTBBB`(t½ ≈ 7 h)의 정량 근거이자 뇌 축적이 가역적일 수 있음의 근거.

79. **Amin SB** et al. *Bilirubin and serial auditory brainstem responses in premature infants.* Pediatrics. 2001.
   [PMID 11335741](https://pubmed.ncbi.nlm.nih.gov/11335741/)  
   → 미숙아에서 빌리루빈과 연속 ABR — Bf 상승 시 ABR 잠복기 연장, 치료 후 회복. 모델의 **가역적** `ABRD` 상태의 근거.

80. **Amin SB** *Clinical assessment of bilirubin-induced neurotoxicity in premature infants.* Semin Perinatol. 2004.
   [PMID 15686265](https://pubmed.ncbi.nlm.nih.gov/15686265/)  
   → 미숙아 빌리루빈 신경독성의 임상 평가 — BIND류 지표와 재태주수 의존성.

81. **Alexandra Brito M** et al. *Bilirubin toxicity to human erythrocytes: a review.* Clin Chim Acta. 2006.
   [PMID 16887110](https://pubmed.ncbi.nlm.nih.gov/16887110/)  
   → 빌리루빈의 적혈구 독성 리뷰 — 막·미토콘드리아 수준의 공통 기전.

82. **Lu Y** et al. *Bilirubin Oxidation End Products (BOXes) Induce Neuronal Oxidative Stress Involving the Nrf2 Pathway.* Oxid Med Cell Longev. 2021.
   [PMID 34373769](https://pubmed.ncbi.nlm.nih.gov/34373769/)  
   → 빌리루빈 산화 최종산물(BOXes)이 Nrf2 경로를 통해 신경 산화 스트레스를 유발 — 모델 `KALT`(비-UGT 산화 분해 경로)가 무해한 배출이 아닐 수 있다는 경고로 인용.

83. **Seidel RA** et al. *Impact of higher-order heme degradation products on hepatic function and hemodynamics.* J Hepatol. 2017.
   [PMID 28412296](https://pubmed.ncbi.nlm.nih.gov/28412296/)  
   → 고차 헴 분해 산물이 간 기능과 혈역학에 미치는 영향 — 위와 같음.

84. **Wisnowski JL** et al. *Magnetic resonance imaging of bilirubin encephalopathy: current limitations and future promise.* Semin Perinatol. 2014.
   [PMID 25267277](https://pubmed.ncbi.nlm.nih.gov/25267277/)  
   → 빌리루빈 뇌증의 MRI — 창백핵 신호 변화의 한계와 전망.

85. **Gburek-Augustat J** et al. *Acute and Chronic Kernicterus: MR Imaging Evolution of Globus Pallidus Signal Change during Childhood.* AJNR Am J Neuroradiol. 2023.
   [PMID 37620154](https://pubmed.ncbi.nlm.nih.gov/37620154/)  
   → 급성·만성 핵황달의 창백핵 신호 변화의 소아기 MRI 진화 — 급성 가역 성분과 만성 비가역 성분의 영상학적 분리.

86. **Yi M** et al. *Globus pallidus/putamen T(1)WI signal intensity ratio in grading and predicting prognosis of neonatal acute bilirubin encephalopathy.* Front Pediatr. 2023.
   [PMID 37842026](https://pubmed.ncbi.nlm.nih.gov/37842026/)  
   → 창백핵/선조체 T1WI 신호비의 등급화와 예후 예측 — `INJ` 지표의 임상 대응물.

87. **Gelineau-Morel R** et al. *Predictive and diagnostic measures for kernicterus spectrum disorder: a prospective cohort study.* Pediatr Res. 2024.
   [PMID 37689774](https://pubmed.ncbi.nlm.nih.gov/37689774/)  
   → 핵황달 스펙트럼 장애의 예측·진단 지표(전향 코호트) — `KERN` 로지스틱 출력의 대조 기준.

88. **Johnson L** et al. *Clinical report from the pilot USA Kernicterus Registry (1992 to 2004).* J Perinatol. 2009.
   [PMID 19177057](https://pubmed.ncbi.nlm.nih.gov/19177057/)  
   → 미국 핵황달 등록사업 임상 보고(1992-2004) — 실제 발생 사례의 TSB 분포와 위험인자 조합.

89. **Bhutani VK** et al. *Synopsis report from the pilot USA Kernicterus Registry.* J Perinatol. 2009.
   [PMID 19177058](https://pubmed.ncbi.nlm.nih.gov/19177058/)  
   → 미국 핵황달 등록사업 개요 보고 — 위와 같음.

---

## 11. G6PD 결핍 — 세계적으로 가장 흔한 핵황달 원인 (G6PD deficiency)

모델 A8의 화학량론적 요점: 헤모글로빈 1 g/dL의 파괴는 29 mg/kg의 빌리루빈, 즉 2.5 dL/kg의 분포 공간에서 **11.6 mg/dL의 TSB**를 만든다. 따라서 3 g/dL의 Hb 감소는 약 35 mg/dL의 색소를 실어오며, 이것이 헤모글로빈이 거의 정상으로 보이는데도 교환수혈 수준의 빌리루빈에 도달하는 정량적 이유다.

90. **Howes RE** et al. *Spatial distribution of G6PD deficiency variants across malaria-endemic regions.* Malar J. 2013.
   [PMID 24228846](https://pubmed.ncbi.nlm.nih.gov/24228846/)  
   → 말라리아 유행 지역의 G6PD 결핍 변이 분포 — `G6PD` 플래그가 전 세계 수억 명에 해당한다는 배경.

91. **Bancone G** et al. *Contribution of genetic factors to high rates of neonatal hyperbilirubinaemia on the Thailand-Myanmar border.* PLOS Glob Public Health. 2022.
   [PMID 36962413](https://pubmed.ncbi.nlm.nih.gov/36962413/)  
   → 태국-미얀마 국경의 높은 신생아 고빌리루빈혈증률에 대한 유전 요인의 기여 — G6PD와 UGT1A1 변이의 결합 효과.

92. **Kaplan M** et al. *Grand Rounds Hyperbilirubinemia following Phototherapy in Glucose-6-Phosphate Dehydrogenase-Deficient Neonates: Not Out of the Woods.* J Pediatr. 2023.
   [PMID 37169338](https://pubmed.ncbi.nlm.nih.gov/37169338/)  
   → G6PD 결핍 신생아에서 광선치료 후의 고빌리루빈혈증 — 광선치료를 끝낸 뒤에도 위험이 남는다는 임상 관찰. 모델의 재상승(A5)과 폐쇄 루프 광선치료(A13) 해석.

93. **Somanath S** et al. *Kernicterus in a late preterm infant with overlooked G6PD deficiency.* BMJ Case Rep. 2025.
   [PMID 40876910](https://pubmed.ncbi.nlm.nih.gov/40876910/)  
   → 간과된 G6PD 결핍이 있는 후기 조산아의 핵황달 — 모델 S6·S7의 임상 대응 사례.

94. **Zahedpasha Y** et al. *Relation between Neonatal Icter and Gilbert Syndrome in Gloucose-6-Phosphate Dehydrogenase Deficient Subjects.* J Clin Diagn Res. 2014.
   [PMID 24783083](https://pubmed.ncbi.nlm.nih.gov/24783083/)  
   → G6PD 결핍자에서 신생아 황달과 Gilbert 증후군의 관계 — **생성 병변과 제거 병변의 곱셈적 상호작용**(모델 A9의 유전형 x 부하 구조).

95. **Zakerihamidi M** et al. *Comparison of severity and prognosis of jaundice due to Rh incompatibility and G6PD deficiency.* Transfus Apher Sci. 2023.
   [PMID 37164807](https://pubmed.ncbi.nlm.nih.gov/37164807/)  
   → Rh 부적합과 G6PD 결핍에 의한 황달의 중증도·예후 비교 — 같은 용혈 부하가 다른 궤적을 만드는 이유(모델 A8: 펄스 대 지속).

96. **Xu JX** et al. *Etiology analysis for term newborns with severe hyperbilirubinemia in eastern Guangdong of China.* World J Clin Cases. 2023.
   [PMID 37123300](https://pubmed.ncbi.nlm.nih.gov/37123300/)  
   → 중국 동부 광둥의 만삭아 중증 고빌리루빈혈증 원인 분석 — 원인별 기여도의 실제 분포.

---

## 12. 체표면적과 성장 — 제4 주제의 정량적 근거 (Allometry and growth)

모델의 모든 광반응 속도는 `광자 전달 x 조사층 내 농도`이며 결코 `x 양`이 아니다. 따라서 광자 전달은 체표면적에 비례하고 kg당 광제거율은 `BSA/W ~ W^-0.30`으로 떨어진다. **그런데 A11b는 이 항이 kg당 빌리루빈 생성량의 감소(7.7 → 3.4 mg/kg/day)로 거의 정확히 상쇄됨을 보인다** — 즉 '표면적/질량비가 줄어서 램프가 듣지 않게 된다'는 통상의 설명만으로는 임상 경과를 설명할 수 없다.

97. **Haycock GB** et al. *Geometric method for measuring body surface area: a height-weight formula validated in infants, children, and adults.* J Pediatr. 1978.
   [PMID 650346](https://pubmed.ncbi.nlm.nih.gov/650346/)  
   → Haycock 체표면적 공식(신생아·소아·성인에서 검증) — 모델 `bsa_cm2()`의 출처. BSA = 0.024265 x W^0.5378 x H^0.3964.

---

## 모델이 의도적으로 다루지 않는 것 (Explicit scope limits)

아래 항목은 문헌이 존재하지만 모델에 **넣지 않았다**. 결과를 해석할 때의 경계다.

- **깊고 느린 조직 구획.** 혈관외 구획을 단일 균질 구획으로 두었기 때문에 이중용적
  교환수혈이 체내 총 부하의 40-45 %를 제거한다(고전적 교육 수치 ~25 %). 같은 이유로
  CN-I의 광선치료 한계값이 낙관적이다(4.9-8.4 대 보고된 15-25 mg/dL). 두 편향은
  **같은 원인**을 가지며 `PSXKG` 하나로 동시에 맞출 수 없다.
- **항-AAV 중화항체**에 의한 유전자 치료 실패(위 31502485).
- **청동아 증후군**의 담즙정체 기전은 지도에는 있으나 ODE에는 `FCHOL` 계수로만 존재한다.
- **모아 분리·수유 중단**의 행동적 결과, 광선치료의 망막·DNA 산화 손상은 지도에만 있다.
- **인종·피부 색소에 따른 TcB 편향**은 출력에 반영되지 않았다(위 41490938).
- **AAP 2022 임계값은 해석적 근사**이며 표값이 아니다. 임상 판단에 사용할 수 없다.

---

## ⚠️ 면책

본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이다. 공개 문헌을 근거로
구성했으나 독립적으로 검증·인증되지 않았으며, **실제 임상 의사결정, 처방, 규제**
**제출에 직접 사용해서는 안 된다.**
