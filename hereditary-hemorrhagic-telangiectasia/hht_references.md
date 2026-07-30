# 유전성 출혈성 모세혈관확장증 (HHT) — 참고문헌

**Hereditary Haemorrhagic Telangiectasia (Rendu-Osler-Weber disease) — annotated references**


이 목록은 `hht_mrgsolve_model.R`의 각 파라미터·구조·검증 축이 실제로 어디에서
왔는지를 추적할 수 있도록 정리한 것이다. 굵게 표시된 주석은 그 문헌이 모델의
어느 부분을 결정했는지를 명시한다. **모든 PMID는 PubMed E-utilities API로
직접 조회하여 서지사항을 확인했으며, 조회되지 않은 항목은 추측하지 않고
목록에서 제외했다.**

모델이 정량적으로 재현하거나 명시적으로 반박하는 핵심 수치:

| 출처 | 관찰값 | 모델에서의 역할 |
|------|--------|-----------------|
| PATH-HHT (PMID 39292928) | ESS −1.84 vs 위약 −0.90, 차이 −0.94 | 1차 보정 — 위약이 치료군 개선의 48.9%를 차지 |
| PATH-HHT | 4주당 정맥철 중앙값 333 mg(위약) vs 0 mg | 수요기반 철 정책의 **예측값**(모델 327 mg)이자 평균 출혈량 31 mL/일의 역산 근거 |
| PATH-HHT | 혈색소 24주 차이 없음, 중단 4주 후 +1.09 | AXIS 2 — 철 구조요법이 약효와 같은 축에서 반대 방향으로 작용하여 상쇄 |
| PATH-HHT 생검 | 내피증식 억제 없음, 벽세포 피복 증가 | 포말리도마이드를 혈관 성숙제로 구현 |
| Dupuis-Girod 2012 (PMID 22396517) | 심장지수 5.05→4.2, 코피 221→134→43 분/월 | AXIS 4(기각된 가설) · AXIS 6(교란) |
| Azzopardi 2015 (PMID 25751241) | 3/2/1개월 유지 → 41/45/50% 및 34/43/60% | AXIS 5 — 모델이 **재현하지 못하는** 투여간격 의존성 |
| NOSE (PMID 27599329) | 식염수군도 ESS 유의 개선 | 위약 표류의 독립적 근거 |
| Hoag 2010 (PMID 20087969) | ESS에 빈혈·수혈 항목 포함 | AXIS 3 — 점수 자체가 철에 반응 |
| Lebrin 2010 (PMID 20364125) | 탈리도마이드가 혈관 성숙 촉진 | MURALC 상태변수의 기전적 선례 |
| Lu 2008 (PMID 18205003) | 반감기 ~20일, CL ~0.2 L/일 | 몰 과량 8.5×10⁴ 계산 → AXIS 4의 가설 기각 |


---

## 1. 진단 기준 · 국제 진료지침 (Diagnosis and international guidelines)

1. Faughnan ME et al. *Second International Guidelines for the Diagnosis and Management of Hereditary Hemorrhagic Telangiectasia*. Ann Intern Med 174:1035-1036 (2021). [PMID 34280351](https://pubmed.ncbi.nlm.nih.gov/34280351/)
   - 제2차 국제 진료지침(Ann Intern Med 2021). 폐 동정맥기형은 조영증강 심초음파로 **전원 선별**하도록 권고하며 맥박산소측정으로 대체하지 않는다 — AXIS 7이 정량적으로 재현하는 권고.

2. Faughnan ME et al. *Second International Guidelines for the Diagnosis and Management of Hereditary Hemorrhagic Telangiectasia*. Ann Intern Med 173:989-1001 (2020). [PMID 32894695](https://pubmed.ncbi.nlm.nih.gov/32894695/)
   - 제2차 국제 진료지침 원문. 코피·철결핍·빈혈 관리, 장기별 선별검사 간격, 임신·시술 시 주의사항.

3. Rössler J et al. *EHA Endorsement of the Second International Guidelines for the Diagnosis and Management of Hereditary Hemorrhagic Telangiectasia*. Hemasphere 5:e647 (2021). [PMID 34646980](https://pubmed.ncbi.nlm.nih.gov/34646980/)

4. Adam MP et al. *Hereditary Hemorrhagic Telangiectasia*.  (1993). [PMID 20301525](https://pubmed.ncbi.nlm.nih.gov/20301525/)

5. Anderson E et al. *Pulmonary arteriovenous malformations may be the only clinical criterion present in genetically confirmed hereditary haemorrhagic telangiectasia*. Thorax 77:628-630 (2022). [PMID 35165143](https://pubmed.ncbi.nlm.nih.gov/35165143/)

6. Westermann CJ et al. *The prevalence and manifestations of hereditary hemorrhagic telangiectasia in the Afro-Caribbean population of the Netherlands Antilles: a family screening*. Am J Med Genet A 116A:324-8 (2003). [PMID 12522784](https://pubmed.ncbi.nlm.nih.gov/12522784/)


---

## 2. 유전학 — ENG / ACVRL1 / SMAD4 / GDF2 반접합부족 (Genetics)

7. McAllister KA et al. *Endoglin, a TGF-beta binding protein of endothelial cells, is the gene for hereditary haemorrhagic telangiectasia type 1*. Nat Genet 8:345-51 (1994). [PMID 7894484](https://pubmed.ncbi.nlm.nih.gov/7894484/)
   - 엔도글린이 HHT1의 원인 유전자임을 밝힌 원 논문.

8. McAllister KA et al. *Six novel mutations in the endoglin gene in hereditary hemorrhagic telangiectasia type 1 suggest a dominant-negative effect of receptor function*. Hum Mol Genet 4:1983-5 (1995). [PMID 8595426](https://pubmed.ncbi.nlm.nih.gov/8595426/)

9. Gallione CJ et al. *Mutation and expression analysis of the endoglin gene in hereditary hemorrhagic telangiectasia reveals null alleles*. Hum Mutat 11:286-94 (1998). [PMID 9554745](https://pubmed.ncbi.nlm.nih.gov/9554745/)

10. Abdalla SA et al. *Disease-associated mutations in conserved residues of ALK-1 kinase domain*. Eur J Hum Genet 11:279-87 (2003). [PMID 12700602](https://pubmed.ncbi.nlm.nih.gov/12700602/)

11. Gallione CJ et al. *SMAD4 mutations found in unselected HHT patients*. J Med Genet 43:793-7 (2006). [PMID 16613914](https://pubmed.ncbi.nlm.nih.gov/16613914/)
   - 선별되지 않은 HHT 환자에서 발견된 SMAD4 변이 — JP-HHT 중복 표현형.

12. Gallione C et al. *Overlapping spectra of SMAD4 mutations in juvenile polyposis (JP) and JP-HHT syndrome*. Am J Med Genet A 152A:333-9 (2010). [PMID 20101697](https://pubmed.ncbi.nlm.nih.gov/20101697/)

13. Pawlikowska L et al. *The ACVRL1 c.314-35A>G polymorphism is associated with organ vascular malformations in hereditary hemorrhagic telangiectasia patients with ENG mutations, but not in patients with ACVRL1 mutations*. Am J Med Genet A 167:1262-7 (2015). [PMID 25847705](https://pubmed.ncbi.nlm.nih.gov/25847705/)

14. Kilian A et al. *Genotype-Phenotype Correlations in Children with HHT*. J Clin Med 9 (2020). [PMID 32842615](https://pubmed.ncbi.nlm.nih.gov/32842615/)
   - 소아 HHT의 유전형-표현형 상관. ENG는 폐·뇌 동정맥기형, ACVRL1은 간 혈관기형으로 치우친다는 장기 편향의 근거 — AXIS 9.

15. Karlsson T et al. *Mutations in the ENG, ACVRL1, and SMAD4 genes and clinical manifestations of hereditary haemorrhagic telangiectasia: experience from the Center for Osler's Disease, Uppsala University Hospital*. Ups J Med Sci 123:153-157 (2018). [PMID 30251589](https://pubmed.ncbi.nlm.nih.gov/30251589/)
   - ENG·ACVRL1·SMAD4 변이별 임상 소견(오슬러 센터 경험).

16. Sturiale CL et al. *Genotype-phenotype correlations and protein domain-level predictors of cerebrovascular malformations in hereditary hemorrhagic telangiectasia*. J Neurol 273 (2026). [PMID 41915210](https://pubmed.ncbi.nlm.nih.gov/41915210/)
   - 유전형-표현형 상관 및 단백질 도메인 수준 예측인자(뇌혈관기형).


---

## 3. BMP9–ALK1–ENG–SMAD1/5/8 신호전달과 혈관기형 형성 기전 (Mechanism)

17. Arthur HM et al. *An update on preclinical models of hereditary haemorrhagic telangiectasia: Insights into disease mechanisms*. Front Med (Lausanne) 9:973964 (2022). [PMID 36250069](https://pubmed.ncbi.nlm.nih.gov/36250069/)
   - HHT 전임상 모형 최신 종설: BMP9/10–ALK1–ENG–SMAD1/5/8 축, second hit, 혈류 의존성.

18. Cirulli A et al. *Vascular endothelial growth factor serum levels are elevated in patients with hereditary hemorrhagic telangiectasia*. Acta Haematol 110:29-32 (2003). [PMID 12975554](https://pubmed.ncbi.nlm.nih.gov/12975554/)
   - HHT 환자에서 혈청 VEGF가 상승되어 있음. 본 모델 VEGF 기저값의 근거.

19. Kaitu'u-Lino TJ et al. *MT-MMPs in pre-eclamptic placenta: relationship to soluble endoglin production*. Placenta 34:168-73 (2013). [PMID 23261267](https://pubmed.ncbi.nlm.nih.gov/23261267/)
   - MT-MMP에 의한 가용성 엔도글린(sENG) 절단. 지도의 ENDOSOL 노드 근거.

20. Ruiz S et al. *Tacrolimus rescues the signaling and gene expression signature of endothelial ALK1 loss-of-function and improves HHT vascular pathology*. Hum Mol Genet 26:4786-4798 (2017). [PMID 28973643](https://pubmed.ncbi.nlm.nih.gov/28973643/)
   - 타크로리무스가 ALK1 기능소실 내피세포의 신호전달과 유전자 발현 서명을 회복시킨다. 본 모델에서 타크로리무스가 SMAD1/5/8 신호강도를 직접 올리는 유일한 약물로 구현된 근거 — 결과가 아니라 병변 자체를 겨냥하는 유일한 기전.


---

## 4. 임상 표현형 · 자연사 · 역학 (Clinical spectrum, natural history, epidemiology)

21. Al-Samkari H et al. *Clinical spectrum of hereditary hemorrhagic telangiectasia: data from the Comprehensive HHT Outcomes Registry of the US (CHORUS)*. Blood 148:417-432 (2026). [PMID 41843464](https://pubmed.ncbi.nlm.nih.gov/41843464/)
   - CHORUS(미국 HHT 통합 결과 레지스트리) 임상 스펙트럼. 장기별 동정맥기형 유병률과 출혈 부담의 최근 대규모 기술.

22. Donaldson JW et al. *Complications and mortality in hereditary hemorrhagic telangiectasia: A population-based study*. Neurology 84:1886-93 (2015). [PMID 25862798](https://pubmed.ncbi.nlm.nih.gov/25862798/)
   - 인구기반 연구: HHT의 합병증과 사망률. 본 모델이 예측하는 자연사(수혈 의존으로의 진행)의 임상적 대응.

23. Shovlin CL et al. *Reported cardiac phenotypes in hereditary hemorrhagic telangiectasia emphasize burdens from arrhythmias, anemia and its treatments, but suggest reduced rates of myocardial infarction*. Int J Cardiol 215:179-85 (2016). [PMID 27116331](https://pubmed.ncbi.nlm.nih.gov/27116331/)
   - HHT의 심장 표현형 — 부정맥·빈혈과 그 치료의 부담. 항응고 역설(AF·정맥혈전증이 흔한데 질병 자체가 출혈이라는 문제)의 근거.

24. Silva BM et al. *Lifestyle and dietary influences on nosebleed severity in hereditary hemorrhagic telangiectasia*. Laryngoscope 123:1092-9 (2013). [PMID 23404156](https://pubmed.ncbi.nlm.nih.gov/23404156/)

25. Squitín Tasende M et al. *Hereditary hemorrhagic telangiectasia in pediatrics: descriptive study in a specialized unit*. Arch Argent Pediatr 123:e202510661 (2025). [PMID 40762448](https://pubmed.ncbi.nlm.nih.gov/40762448/)

26. Schleupner MC et al. *Cobalamin and iron deficiency still presents a challenge in hereditary hemorrhagic telangiectasia*. Sci Rep 15:28463 (2025). [PMID 40759730](https://pubmed.ncbi.nlm.nih.gov/40759730/)


---

## 5. 코피 중증도 점수(ESS)와 삶의 질 지표 — 본 모델의 1차 평가변수 (Endpoints)

27. Hoag JB et al. *An epistaxis severity score for hereditary hemorrhagic telangiectasia*. Laryngoscope 120:838-43 (2010). [PMID 20087969](https://pubmed.ncbi.nlm.nih.gov/20087969/)
   - ESS(Epistaxis Severity Score) 원 개발 논문, 21개국 900명. 중증도 결정 인자로 빈도·지속시간·강도·의료개입 필요성과 함께 **빈혈**과 **수혈 필요성**을 포함한다. 가중치가 closed form으로 공표되지 않았기 때문에 본 모델은 6항목 대리척도를 쓰고 빈혈+수혈 가중치 비중(W_FE_SHARE)을 고정하지 않고 sweep한다 — AXIS 3.

28. Yin LX et al. *The minimal important difference of the epistaxis severity score in hereditary hemorrhagic telangiectasia*. Laryngoscope 126:1029-32 (2016). [PMID 26393959](https://pubmed.ncbi.nlm.nih.gov/26393959/)
   - ESS의 최소 임상적 중요 차이(MCID) 0.71점. PATH-HHT가 −0.94를 임상적으로 의미있다고 판정한 기준.

29. Hayama M et al. *Validation of epistaxis severity score for hereditary hemorrhagic telangiectasia in Japan*. Auris Nasus Larynx 49:415-420 (2022). [PMID 34857410](https://pubmed.ncbi.nlm.nih.gov/34857410/)

30. Merlo CA et al. *The effects of epistaxis on health-related quality of life in patients with hereditary hemorrhagic telangiectasia*. Int Forum Allergy Rhinol 4:921-5 (2014). [PMID 25145809](https://pubmed.ncbi.nlm.nih.gov/25145809/)

31. Al-Samkari H et al. *Validation and clinical application of the hereditary hemorrhagic telangiectasia-specific quality of life scale*. J Thromb Haemost 24:108-118 (2026). [PMID 41092987](https://pubmed.ncbi.nlm.nih.gov/41092987/)
   - HHT 특이 삶의 질 척도(HHT-QoL)의 검증 — PATH-HHT의 주요 2차 평가변수.


---

## 6. 무작위배정 임상시험 — 위약효과의 크기가 기록된 두 연구 (Randomised trials)

32. Al-Samkari H et al. *Pomalidomide for Epistaxis in Hereditary Hemorrhagic Telangiectasia*. N Engl J Med 391:1015-1027 (2024). [PMID 39292928](https://pubmed.ncbi.nlm.nih.gov/39292928/)
   - **본 모델의 1차 보정 근거.** PATH-HHT: 포말리도마이드 4 mg/일 vs 위약 2:1, n=144, 24주. ESS 변화 −1.84 (95% CI −2.25~−1.43) vs −0.90 (−1.39~−0.40), 차이 −0.94 (−1.57~−0.31), p=0.004. HHT-QoL −2.7 vs −1.2 (차이 −1.4). 기저 ESS 5.0±1.5, HHT-QoL 6.3±3.1, 연령 58.8±12.2, 69% 빈혈, 이전 6개월간 84% 정맥철·19% 수혈. 12~24주 수혈률 9% vs 18%, 4주당 정맥철 중앙값 0 mg vs 333 mg. **혈색소는 24주에 군간 차이 없었고, 투약 중단 4주 후에야 +1.09 g/dL (0.38~1.80)로 벌어졌다** — AXIS 2 전체가 이 한 쌍의 관찰에서 나온다. 코 점막 생검에서 내피세포 증식 억제는 없었고 벽세포(mural cell) 피복 증가만 관찰되어, 본 모델이 포말리도마이드를 항혈관신생제가 아니라 혈관 성숙제로 구현한 근거가 되었다.

33. Zhang E et al. *Pomalidomide for hereditary hemorrhagic telangiectasia: after trial longitudinal assessment study (PATH-HHT ATLAS)*. Blood Adv 10:1799-1808 (2026). [PMID 41512167](https://pubmed.ncbi.nlm.nih.gov/41512167/)

34. Al-Samkari H et al. *Characteristics associated with clinical response to pomalidomide in hereditary hemorrhagic telangiectasia*. Blood Adv 10:2967-2976 (2026). [PMID 41719457](https://pubmed.ncbi.nlm.nih.gov/41719457/)

35. Whitehead KJ et al. *Effect of Topical Intranasal Therapy on Epistaxis Frequency in Patients With Hereditary Hemorrhagic Telangiectasia: A Randomized Clinical Trial*. JAMA 316:943-51 (2016). [PMID 27599329](https://pubmed.ncbi.nlm.nih.gov/27599329/)
   - **위약 축의 독립적 근거.** NOSE 시험, n=121, 이중맹검 위약대조: 국소 베바시주맙·에스트리올·트라넥삼산 모두 생리식염수와 차이가 없었고, **모든 군(식염수 포함)에서 12·24주 ESS가 유의하게 개선**되었다. HHT 코피 지표의 비약물성 표류가 실재하며 크다는 것을 무작위 설계로 보여준 연구.

36. Gaillard S et al. *Tranexamic acid for epistaxis in hereditary hemorrhagic telangiectasia patients: a European cross-over controlled trial in a rare disease*. J Thromb Haemost 12:1494-502 (2014). [PMID 25040799](https://pubmed.ncbi.nlm.nih.gov/25040799/)


---

## 7. 베바시주맙 — 단일군 연구와 노출-반응 분석 (Bevacizumab)

37. Dupuis-Girod S et al. *Bevacizumab in patients with hereditary hemorrhagic telangiectasia and severe hepatic vascular malformations and high cardiac output*. JAMA 307:948-55 (2012). [PMID 22396517](https://pubmed.ncbi.nlm.nih.gov/22396517/)
   - **AXIS 4·6의 대상 데이터.** 단일군 2상, n=25, 베바시주맙 5 mg/kg 격주 6회(마지막 투여 70일). 심장지수 중앙값 5.05 → 4.2(3개월) → 4.1(6개월) L/min/m²; 코피 지속시간 평균 221 → 134(3개월) → 43(6개월) 분/월. **최저값이 마지막 투여 110일 후에 나타난다는 점**이 본 모델이 검증하려 한 구조적 주장의 출발점이었다(그리고 그 주장은 기각되었다).

38. Azzopardi N et al. *Dose - response relationship of bevacizumab in hereditary hemorrhagic telangiectasia*. MAbs 7:630-7 (2015). [PMID 25751241](https://pubmed.ncbi.nlm.nih.gov/25751241/)
   - **AXIS 5의 대상 데이터.** 같은 25명에 대한 노출-반응 분석. 심장지수는 transit compartment(지연), 코피 지속시간은 direct inhibition(무지연)으로 모형화. 3·2·1개월 유지투여 시뮬레이션에서 24개월 시점 심장지수 <4 L/min/m² 달성률 41/45/50%, 코피 <20분/월 달성률 34/43/60%. 본 모델은 이 투여간격 의존성을 재현하지 못하며, 그 불일치가 '노출'을 농도로 볼 것인가 표적점유율로 볼 것인가의 구조적 차이에서 나온다는 점을 AXIS 5에서 다룬다.

39. Al-Samkari H et al. *An international survey to evaluate systemic bevacizumab for chronic bleeding in hereditary haemorrhagic telangiectasia*. Haemophilia 26:1038-1045 (2020). [PMID 32432841](https://pubmed.ncbi.nlm.nih.gov/32432841/)
   - 국제 HHT 우수센터 설문(291명 치료 경험): 모든 센터가 유도요법에 5 mg/kg를 쓰고 대부분 격주 6회를 투여하나 유지 요법은 표준화되어 있지 않다.

40. Haahr PD et al. *Availability, use, efficacy and safety of bevacizumab in European hereditary haemorrhagic telangiectasia centres*. Br J Clin Pharmacol 92:535-544 (2026). [PMID 40985322](https://pubmed.ncbi.nlm.nih.gov/40985322/)

41. Buscarini E et al. *Safety of thalidomide and bevacizumab in patients with hereditary hemorrhagic telangiectasia*. Orphanet J Rare Dis 14:28 (2019). [PMID 30717761](https://pubmed.ncbi.nlm.nih.gov/30717761/)

42. Lu JF et al. *Clinical pharmacokinetics of bevacizumab in patients with solid tumors*. Cancer Chemother Pharmacol 62:779-86 (2008). [PMID 18205003](https://pubmed.ncbi.nlm.nih.gov/18205003/)
   - 베바시주맙 임상 약동학: 반감기 약 20일, 청소율 약 0.2 L/일, 중심구획 약 3 L. 본 모델 2구획 PK 파라미터의 출처이며, AXIS 4의 몰 과량 계산이 여기에 의존한다.


---

## 8. 그 외 전신 약물치료 (Other systemic therapies)

43. Lebrin F et al. *Thalidomide stimulates vessel maturation and reduces epistaxis in individuals with hereditary hemorrhagic telangiectasia*. Nat Med 16:420-8 (2010). [PMID 20364125](https://pubmed.ncbi.nlm.nih.gov/20364125/)
   - **포말리도마이드 기전 구현의 직접적 선례.** 탈리도마이드가 혈관 성숙(mural cell 피복)을 촉진하여 HHT 코피를 감소시킨다는 것을 동물모형과 환자에서 함께 보인 연구. 본 모델의 MURALC 상태변수와 그것이 파열빈도·출혈유속·지속시간 세 축에 동시에 작용하는 구조가 여기서 나왔다.

44. Parambil JG et al. *Pazopanib for severe bleeding and transfusion-dependent anemia in hereditary hemorrhagic telangiectasia*. Angiogenesis 25:87-97 (2022). [PMID 34292451](https://pubmed.ncbi.nlm.nih.gov/34292451/)

45. Faughnan ME et al. *Pazopanib may reduce bleeding in hereditary hemorrhagic telangiectasia*. Angiogenesis 22:145-155 (2019). [PMID 30191360](https://pubmed.ncbi.nlm.nih.gov/30191360/)

46. Al-Samkari H. *Systemic Antiangiogenic Therapies for Bleeding in Hereditary Hemorrhagic Telangiectasia: A Practical, Evidence-Based Guide for Clinicians*. Semin Thromb Hemost 48:514-528 (2022). [PMID 35226946](https://pubmed.ncbi.nlm.nih.gov/35226946/)
   - HHT 출혈에 대한 전신 항혈관신생 치료의 실무 지침 — 유도·유지 용량과 반응 평가.

47. Albiñana V et al. *Review of Pharmacological Strategies with Repurposed Drugs for Hereditary Hemorrhagic Telangiectasia Related Bleeding*. J Clin Med 9 (2020). [PMID 32517280](https://pubmed.ncbi.nlm.nih.gov/32517280/)

48. Hsu YP et al. *Medical Treatment for Epistaxis in Hereditary Hemorrhagic Telangiectasia: A Meta-analysis*. Otolaryngol Head Neck Surg 160:22-35 (2019). [PMID 30200816](https://pubmed.ncbi.nlm.nih.gov/30200816/)


---

## 9. 국소·시술적 치료 (Local and procedural treatment)

49. Lund VJ et al. *Nasal closure for severe hereditary haemorrhagic telangiectasia in 100 patients. The Lund modification of the Young's procedure: a 22-year experience*. Rhinology 55:135-141 (2017). [PMID 28064338](https://pubmed.ncbi.nlm.nih.gov/28064338/)
   - 영(Young) 술식 코폐쇄 100례, 22년 경험. AXIS 8이 정량화하는 시술: B_nose를 0으로 만들지만 B_gi는 남긴다.

50. Andersen JH et al. *Patient-recorded benefit from nasal closure in a Danish cohort of patients with hereditary haemorrhagic telangiectasia*. Eur Arch Otorhinolaryngol 277:791-800 (2020). [PMID 31845036](https://pubmed.ncbi.nlm.nih.gov/31845036/)

51. Manfredi G et al. *Gastrointestinal bleeding in patients with hereditary hemorrhagic telangiectasia: Long-term results of endoscopic treatment*. Endosc Int Open 11:E1145-E1152 (2023). [PMID 38108019](https://pubmed.ncbi.nlm.nih.gov/38108019/)
   - HHT 위장관 출혈의 내시경적 치료 장기 성적. 본 모델 GIENDOF 파라미터의 대응.


---

## 10. 폐 동정맥기형 — 역설적 색전과 산소포화도의 해리 (Pulmonary AVM)

52. Al-Saleh S et al. *Screening for pulmonary and cerebral arteriovenous malformations in children with hereditary haemorrhagic telangiectasia*. Eur Respir J 34:875-81 (2009). [PMID 19386691](https://pubmed.ncbi.nlm.nih.gov/19386691/)

53. Topiwala KK et al. *Ischemic Stroke in Patients With Pulmonary Arteriovenous Fistulas*. Stroke 52:e311-e315 (2021). [PMID 34082575](https://pubmed.ncbi.nlm.nih.gov/34082575/)
   - 폐 동정맥루 환자의 허혈성 뇌졸중. 역설적 색전 위험이 산소포화도와 무관하게 해부학으로 결정된다는 AXIS 7 주장의 임상적 근거.

54. Bodilsen J et al. *Pulmonary arteriovenous malformations in patients with previous brain abscess: a cross-sectional population-based study*. Eur J Neurol 31:e16176 (2024). [PMID 38064178](https://pubmed.ncbi.nlm.nih.gov/38064178/)
   - 뇌 농양 병력이 있는 환자에서의 폐 동정맥기형 — 인구기반 단면연구. 색전 경로의 역방향 근거.

55. Le Banner E et al. *A contemporary picture of bacterial infections in patients with hereditary hemorrhagic telangiectasia: A nationwide cohort study*. J Infect 92:106686 (2026). [PMID 41544822](https://pubmed.ncbi.nlm.nih.gov/41544822/)

56. Hering K et al. *Posttreatment Monitoring of Pulmonary Arteriovenous Malformations: Challenges and Approaches*. Chest 168:1471-1480 (2025). [PMID 40543746](https://pubmed.ncbi.nlm.nih.gov/40543746/)

57. Dupuis-Girod S et al. *The Lung in Hereditary Hemorrhagic Telangiectasia*. Respiration 94:315-330 (2017). [PMID 28850955](https://pubmed.ncbi.nlm.nih.gov/28850955/)

58. Jutant EM et al. *Pulmonary hypertension associated with hereditary haemorrhagic telangiectasia: from genetics to clinical management*. Eur Respir J 67 (2026). [PMID 41713948](https://pubmed.ncbi.nlm.nih.gov/41713948/)
   - HHT 관련 폐고혈압: 유전학에서 임상관리까지. 본 모델이 다루는 후모세혈관성 고심박출성 폐고혈압과, 범위 밖으로 명시한 ACVRL1 연관 진성 PAH의 구분.


---

## 11. 간 혈관기형과 고심박출성 심부전 (Hepatic VM and high-output failure)

59. Haahr PD et al. *Haemodynamic alterations associated with hepatic vascular malformations in hereditary haemorrhagic telangiectasia: A cohort study*. JHEP Rep 8:101867 (2026). [PMID 42019870](https://pubmed.ncbi.nlm.nih.gov/42019870/)
   - HHT 간 혈관기형의 혈역학적 변화 코호트 연구. 심장지수·간동맥 확장·문맥압의 관계.

60. Muller YD et al. *Hereditary haemorrhagic telangiectasia: to transplant or not to transplant - is there a right time for liver transplantation?*. Liver Int 36:1735-1740 (2016). [PMID 27864873](https://pubmed.ncbi.nlm.nih.gov/27864873/)
   - HHT에서 간이식의 시점 — 고심박출성 심부전의 유일한 근治적 치료가 왜 약물이 아닌지에 대한 임상적 논의. AXIS 5가 FRES_L ≈ 0.4–0.55로 정량화하는 사실.

61. Feng Q et al. *Imaging features, diagnosis, and clinical management of hepatic involvement in hereditary hemorrhagic telangiectasia: a case series report and literature review*. Front Med (Lausanne) 13:1869789 (2026). [PMID 42517041](https://pubmed.ncbi.nlm.nih.gov/42517041/)


---

## 12. 뇌·척수 혈관기형 (Cerebral and spinal vascular malformations)

62. Kofoed MS et al. *Cerebral vascular malformation screening in hereditary hemorrhagic telangiectasia: Balancing low diagnostic yield against high-risk hemorrhage*. Clin Neurol Neurosurg 264:109348 (2026). [PMID 41713125](https://pubmed.ncbi.nlm.nih.gov/41713125/)

63. Palermo M et al. *Cerebrovascular Malformations Associated With Hereditary Hemorrhagic Telangiectasia and HHT-Like Syndromes: A Comparative Overview*. Eur J Neurol 33:e70523 (2026). [PMID 41704211](https://pubmed.ncbi.nlm.nih.gov/41704211/)

64. Xie K et al. *Decoding brain arteriovenous malformations: from genetic insights to modeling the vascular maze*. Acta Neuropathol Commun 14 (2026). [PMID 42106885](https://pubmed.ncbi.nlm.nih.gov/42106885/)


---

## 13. 철 대사 — 헵시딘 축과 철 제한적 적혈구생성 (Iron metabolism)

65. Ganz T. *Systemic Iron Metabolism*. Adv Exp Med Biol 1480:33-45 (2025). [PMID 40603782](https://pubmed.ncbi.nlm.nih.gov/40603782/)
   - 전신 철 대사 종설. 헵시딘–페로포틴 축, 경구 철 흡수의 포화, 정맥철이 헵시딘 관문을 우회한다는 점 — 본 모델 A_net 구조의 생리학적 근거.

66. Camaschella C et al. *The mutual crosstalk between iron and erythropoiesis*. Int J Hematol 116:182-191 (2022). [PMID 35618957](https://pubmed.ncbi.nlm.nih.gov/35618957/)
   - 철과 적혈구생성의 상호 크로스토크(에리스로페론). 본 모델에서 HHT 빈혈이 골수 용량이 아니라 **철 공급**에 의해 제한된다는 구조의 근거.


---

총 **66편**. 모든 링크는 PubMed 영구 링크이며 서지사항은 API 조회 결과와 일치한다.


## 이 목록에 의도적으로 넣지 않은 것

- **ESS 항목 가중치의 closed form.** 검색으로 찾을 수 없었다. 추정치를 만들어
  인용하는 대신 모델에서 해당 비중을 sweep하는 파라미터로 남겼다.
- **HHT 코호트의 mL/일 단위 실측 출혈량.** 그런 측정을 보고한 연구를 찾지
  못했다. 본 모델의 B_nose·B_gi·f_GI는 모두 철 균형을 통해 역산한 값이며,
  관측값이 아니다.
- **빈혈 보상 지수 α의 HHT 내 추정치.** Fick 원리상 심장지수는 Hb^(−α)로
  변하지만, HHT 코호트에서 심장지수와 혈색소를 같은 시점에 짝지어 보고한
  연구가 없다. AXIS 6은 점추정 대신 α 범위 전체를 보고한다.
- **최저 유효 트로프 유리 VEGF 농도.** AXIS 5가 지목하는 결정적 측정값이지만
  보고된 바 없다. 3개월 트로프에서 유리 VEGF가 여전히 억제되어 있는지 여부가
  베바시주맙 유지 간격 문제를 판정한다.

## 면책 (Disclaimer)

본 참고문헌 목록과 이에 근거한 모델은 **교육·연구 목적**이다. 임상 의사결정,
처방, 규제 제출에 직접 사용해서는 안 된다.
