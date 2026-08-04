# 미숙아 동맥관 개존증 (Patent Ductus Arteriosus of Prematurity) — 참고문헌

> QSP 모델 `pda_qsp_model.dot` · `pda_mrgsolve_model.R` · `pda_reference_model.py` 의 근거 문헌.
> 모든 PMID는 NCBI E-utilities로 개별 조회하여 제목·저자·저널·연도를 확인했습니다
> (memory에서 옮겨 적은 것이 아닙니다). 본문에 인용된 수치는 해당 논문의
> abstract를 직접 받아 대조했으며, 대조 과정에서 실제로 두 개의 오류를
> 발견하여 고쳤습니다 — 아래 §1 각주 참조.
>
> All PMIDs were resolved individually through NCBI E-utilities and the numbers
> quoted below were checked against the fetched abstracts. That check caught two
> errors of recollection, both corrected — see the footnote to §1.

---

## 1. 모델을 교정(calibrate)하거나 반증(falsify)하는 데 쓴 핵심 임상시험

| # | 문헌 | PMID | 모델에서의 역할 |
|---|------|------|-----------------|
| 1 | **Semberova J, et al.** Spontaneous Closure of Patent Ductus Arteriosus in Infants ≤1500 g. *Pediatrics* 2017;140(2):e20164258 | [28701390](https://pubmed.ncbi.nlm.nih.gov/28701390/) | **적합 (FITTED).** 진정한 보존적 관리 280명 중 237명(85%)이 퇴원 전 자연 폐쇄. 폐쇄까지 중앙값 **71일(<26+0주) · 13일(26+0–27+6) · 8일(28+0–29+6) · 6일(≥30주)**. 이 4개 값이 `TAUSYN0` · `KTAUSYN` · `KINVGA` · `NET50` · `NETW` 를 결정합니다. |
| 2 | **Gupta S, et al.** (Baby-OSCAR) Trial of Selective Early Treatment of Patent Ductus Arteriosus with Ibuprofen. *N Engl J Med* 2024;390(4):314–325 | [38265644](https://pubmed.ncbi.nlm.nih.gov/38265644/) | **적합 (FITTED).** 653명, 23–28주, 72시간 내 큰 PDA. 1차 결과(사망 또는 중등도-중증 BPD) **69.2% (이부프로펜) 대 63.5% (위약)**, 보정 RR 1.09 (95% CI 0.98–1.20). 사망 13.6% 대 10.3%. 위험계수 `B_GA` · `B_BUR` · `B_VENT` 3개를 양 군에 동시 적합. |
| 3 | **Gupta S, et al.** Two-year outcomes after selective early treatment of PDA with ibuprofen. *EClinicalMedicine* 2025 | [40896454](https://pubmed.ncbi.nlm.nih.gov/40896454/) | Baby-OSCAR 2년 추적. 단기 결과 무개선이 장기에서도 유지되는지 확인. |
| 4 | **Hundscheid T, et al.** (BeNeDuctus) Expectant Management or Early Ibuprofen for Patent Ductus Arteriosus. *N Engl J Med* 2023;388(11):980–990 | [36477458](https://pubmed.ncbi.nlm.nih.gov/36477458/) | **예측 (PREDICTED).** 273명, 중앙 재태 26주·출생체중 845 g. 복합(NEC/BPD/사망) **46.3% (기대요법) 대 63.5% (조기 이부프로펜)**, 위험차 −17.2%p, 비열등성 충족. BPD 33.3% 대 50.9%. 여기에는 어떤 파라미터도 맞추지 않았습니다. |
| 5 | **Schmidt B, et al.** (TIPP) Long-term effects of indomethacin prophylaxis in extremely-low-birth-weight infants. *N Engl J Med* 2001;344(26):1966–1972 | [11430325](https://pubmed.ncbi.nlm.nih.gov/11430325/) | **예측 (PREDICTED).** 1202명. 예방적 인도메타신은 PDA를 **24% 대 50%**, 중증 IVH를 **9% 대 13%** 로 줄였으나 18개월 사망/장애는 **47% 대 46%** (OR 1.1). 모델의 IVH 항이 서로 반대 부호의 두 항(혈소판 COX-1 억제 = 해, 배아기질 혈관 성숙 = 이득)으로 되어 있는 이유. |
| 6 | **Clyman RI, et al.** (PDA-TOLERATE) An Exploratory Randomized Controlled Trial of Treatment of Moderate-to-Large Patent Ductus Arteriosus at 1 Week of Age. *J Pediatr* 2019;205:41–48 | [30340932](https://pubmed.ncbi.nlm.nih.gov/30340932/) | **예측.** <28주, 6–14일에 조기 상용 치료 대 보존. 결찰/퇴원 시 PDA 32% 대 39%, NEC 16% 대 19%, BPD 49% 대 53%, BPD/사망 58% 대 57%. 구제기준 충족 31% 대 62%. |
| 7 | **Mitra S, et al.** Association of Placebo, Indomethacin, Ibuprofen, and Acetaminophen With Closure of Hemodynamically Significant PDA: Systematic Review and Network Meta-analysis. *JAMA* 2018;319(12):1221–1238 | [29584842](https://pubmed.ncbi.nlm.nih.gov/29584842/) | **적합 (FITTED).** 약제별 폐쇄율 순위가 `KI_IBU` · `KI_IND` · `IC50_APAP` 를 결정. 고용량 경구 이부프로펜이 최상위. |
| 8 | **Buvaneswarran S, et al.** Active Treatment vs Expectant Management of Patent Ductus Arteriosus in Preterm Infants: A Meta-Analysis. *JAMA Pediatr* 2025 | [40423988](https://pubmed.ncbi.nlm.nih.gov/40423988/) | 최신 통합 근거 — "관은 닫히지만 결과는 개선되지 않는다"는 모델의 중심 주장에 대한 외부 검증. |
| 9 | **Mitra S, et al.** Early treatment versus expectant management of hemodynamically significant PDA. *Cochrane Database Syst Rev* 2025 | [40548426](https://pubmed.ncbi.nlm.nih.gov/40548426/) | 조기 치료 대 기대요법 Cochrane 통합. |
| 10 | **Mitra S, et al.** Interventions for patent ductus arteriosus (PDA) in preterm infants: an overview of Cochrane Systematic Reviews. *Cochrane Database Syst Rev* 2023 | [37039501](https://pubmed.ncbi.nlm.nih.gov/37039501/) | 전체 개입 지형 요약. |
| 11 | **Mitra S, et al.** Prophylactic cyclo-oxygenase inhibitor drugs for the prevention of morbidity and mortality in extremely preterm infants. *Arch Dis Child Fetal Neonatal Ed* 2024 | [37419686](https://pubmed.ncbi.nlm.nih.gov/37419686/) | 예방적 COX 억제제 통합 — 시나리오 S7의 근거. |
| 12 | **Hundscheid T, et al.** Conservative Management of Patent Ductus Arteriosus in Preterm Infants — A Systematic Review and Meta-Analyses. *Front Pediatr* 2021 | [33718300](https://pubmed.ncbi.nlm.nih.gov/33718300/) | 보존적 관리의 근거 기반. |

> **각주 — 대조에서 발견한 두 오류.** (i) Baby-OSCAR 1차 결과를 69.4%로 기억하고
> 있었으나 실제 published 값은 **69.2%** 였습니다. (ii) 더 중요한 것은,
> Semberova의 자연 폐쇄 중앙값을 재태연령에 대해 완만한 곡선으로 상정하고
> 26주 목표를 48일로 잡았는데, 논문의 실제 값은 26+0–27+6주에서 **13일**
> 이었습니다. <26주(71일)와 26–27주(13일) 사이의 **5.5배 불연속**은 약 1주의
> 재태연령 차이에서 일어나며, 이는 이 모델이 구조적 폐쇄 문턱(TMAXGA가
> 남기는 잔여 내경이 폐쇄 기준을 넘는 지점, 약 25.5주)을 두고 있는 곳과
> 정확히 일치합니다. 틀린 목표값을 쓰고 있었다는 사실을 발견한 것이,
> 결과적으로 모델 구조에 대한 독립적 지지 근거를 준 셈입니다.

---

## 2. 동맥관의 정상 생리 — 프로스타글란딘 의존적 개존

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 13 | **Coceani F, Olley PM.** Lamb ductus arteriosus: effect of prostaglandin synthesis inhibitors on the muscle tone and the response to prostaglandin E2. *Prostaglandins* 1975 | [1135442](https://pubmed.ncbi.nlm.nih.gov/1135442/) | 원 발견. 동맥관 개존이 국소 프로스타글란딘 합성에 의존한다는 것을 확립. |
| 14 | **Segi E, et al.** Patent ductus arteriosus and neonatal death in prostaglandin receptor EP4-deficient mice. *Biochem Biophys Res Commun* 1998 | [9600059](https://pubmed.ncbi.nlm.nih.gov/9600059/) | EP4 유전자 결손 마우스에서 PDA와 신생기 사망 — 모델의 `EP4` 상태변수가 관 개존의 필수 매개체인 근거. |
| 15 | **Rheinlaender C, et al.** Changing expression of cyclooxygenases and prostaglandin receptor EP4 during development of the human ductus arteriosus. *Pediatr Res* 2006 | [16857763](https://pubmed.ncbi.nlm.nih.gov/16857763/) | 인간 동맥관에서 COX·EP4 발현이 발달에 따라 변화 — `FSYN` 관 국소 합성 involution과 `EP4GA`·`KEP4` 의 직접 근거. |
| 16 | **Trivedi DB, et al.** Attenuated cyclooxygenase-2 expression contributes to patent ductus arteriosus in preterm mice. *Pediatr Res* 2006 | [17065565](https://pubmed.ncbi.nlm.nih.gov/17065565/) | 미숙 마우스에서 COX-2 발현 감소가 PDA에 기여 — 재태연령 의존적 프로스타노이드 항의 근거. |
| 17 | **Momma K, Takeuchi H.** Constriction of fetal ductus arteriosus by nonsteroidal antiinflammatory drugs. *Adv Prostaglandin Thromboxane Leukot Res* 1983 | [6221638](https://pubmed.ncbi.nlm.nih.gov/6221638/) | NSAID에 의한 태아 동맥관 수축 — 약력학의 고전적 근거. |
| 18 | **Printz MP, et al.** Studies of pulmonary prostaglandin biosynthetic and catabolic enzymes as factors in ductus arteriosus patency. *Pediatr Res* 1984 | [6422431](https://pubmed.ncbi.nlm.nih.gov/6422431/) | 폐의 프로스타글란딘 생성·분해 효소가 관 개존의 인자 — `PGDH` 상태변수(폐 15-PGDH)의 근거. |
| 19 | **Takizawa T, et al.** Inhibitory effect of indomethacin on neonatal lung catabolism of prostaglandin E2. *J Toxicol Sci* 1996 | [8959648](https://pubmed.ncbi.nlm.nih.gov/8959648/) | 신생기 폐의 PGE2 분해와 그 억제 — 순환 PGE2 청소가 폐혈류·15-PGDH의 곱으로 쓰인 이유. |
| 20 | **Ovalı F.** Molecular and Mechanical Mechanisms Regulating Ductus Arteriosus Closure in Preterm Infants. *Front Pediatr* 2020 | [32984222](https://pubmed.ncbi.nlm.nih.gov/32984222/) | 미숙아 관 폐쇄의 분자·기계적 기전 종합 리뷰. |

---

## 3. 산소 감지와 평활근 수축 기구

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 21 | **Archer SL, et al.** O2 sensing in the human ductus arteriosus: redox-sensitive K+ channels are regulated by mitochondria-derived hydrogen peroxide. *Biol Chem* 2004 | [15134333](https://pubmed.ncbi.nlm.nih.gov/15134333/) | 인간 동맥관의 산소 감지 — 미토콘드리아 유래 H2O2가 redox 민감성 K+ 채널을 조절. 지도의 `O2SENSE → MEMPOT → CAL` 경로. |
| 22 | **Bentley RET, et al.** The molecular mechanisms of oxygen-sensing in human ductus arteriosus smooth muscle cells: A comprehensive transcriptome profile. *Genomics* 2021 | [34245829](https://pubmed.ncbi.nlm.nih.gov/34245829/) | 인간 관 평활근 산소 감지 전사체 — `P50_O2` · `P50_GA` 의 기전적 배경. |
| 23 | **Kajimoto H, et al.** Oxygen activates the Rho/Rho-kinase pathway and induces RhoB and ROCK-1 expression in human and rabbit ductus arteriosus. *Circulation* 2007 | [17353442](https://pubmed.ncbi.nlm.nih.gov/17353442/) | 산소가 Rho/ROCK 경로를 활성화 — 지도의 `ROCK` (칼슘 감작) 노드. |
| 24 | **Hong Z, et al.** Role of store-operated calcium channels and calcium sensitization in normoxic contraction of the ductus arteriosus. *Circulation* 2006 | [16982938](https://pubmed.ncbi.nlm.nih.gov/16982938/) | 저장 조절 칼슘 채널과 칼슘 감작 — `CASMC` · `MLC20` 경로. |

---

## 4. 수축 ≠ 폐쇄 — 관벽 저산소증과 해부학적 재형성 (모델의 이중안정 스위치)

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 25 | **Kajino H, et al.** Vasa vasorum hypoperfusion is responsible for medial hypoxia and anatomic remodeling in the newborn lamb ductus arteriosus. *Pediatr Res* 2002 | [11809919](https://pubmed.ncbi.nlm.nih.gov/11809919/) | **모델 구조의 핵심 근거.** 수축 후 **vasa vasorum 저관류**가 중막 저산소증과 해부학적 재형성을 일으킴. `WALLO2` 를 "얇은 벽 = 확산으로 충분 / 두꺼운 벽 = vasa 의존, 수축으로 차단"의 두 항으로 쓴 이유. |
| 26 | **Clyman RI, et al.** VEGF regulates remodeling during permanent anatomic closure of the ductus arteriosus. *Am J Physiol Regul Integr Comp Physiol* 2002 | [11742839](https://pubmed.ncbi.nlm.nih.gov/11742839/) | 영구 폐쇄 시 VEGF가 재형성을 조절 — `VEGFD → TGFBD → REMOD` 연쇄. |
| 27 | **Chorne N, et al.** Risk factors for persistent ductus arteriosus patency during indomethacin treatment. *J Pediatr* 2007 | [18035143](https://pubmed.ncbi.nlm.nih.gov/18035143/) | 인도메타신 치료 중 지속 개존의 위험인자 — 재개통·무반응의 임상적 상관. |
| 28 | **Keller RL, Clyman RI.** Persistent Doppler flow predicts lack of response to multiple courses of indomethacin in premature infants with recurrent patent ductus arteriosus. *Pediatrics* 2003 | [12949288](https://pubmed.ncbi.nlm.nih.gov/12949288/) | 반복 치료 무반응의 예측 — 2차 코스 시나리오(S12)의 근거. |
| 29 | **Louis D, et al.** Factors associated with non-response to second course indomethacin for PDA treatment in preterm neonates. *J Matern Fetal Neonatal Med* 2018 | [28391737](https://pubmed.ncbi.nlm.nih.gov/28391737/) | 2차 코스 무반응 인자. |

---

## 5. 약동학 — 미숙아에서의 세 약물

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 30 | **Aranda JV, Thomas R.** Systematic review: intravenous ibuprofen in preterm newborns. *Semin Perinatol* 2006 | [16813969](https://pubmed.ncbi.nlm.nih.gov/16813969/) | 정맥 이부프로펜 체계적 문헌고찰. |
| 31 | **Van Overmeire B, et al.** Ibuprofen pharmacokinetics in preterm infants with patent ductus arteriosus. *Clin Pharmacol Ther* 2001 | [11673749](https://pubmed.ncbi.nlm.nih.gov/11673749/) | 미숙아 이부프로펜 PK. |
| 32 | **Hirt D, et al.** An optimized ibuprofen dosing scheme for preterm neonates with patent ductus arteriosus, based on a population pharmacokinetic and pharmacodynamic study. *Br J Clin Pharmacol* 2008 | [18307541](https://pubmed.ncbi.nlm.nih.gov/18307541/) | popPK/PD 기반 최적 용법 — **출생 후 연령에 따른 청소율 성숙**(`CLMAT_IBU`)과 후기 치료 시 고용량 필요성의 근거. |
| 33 | **Yaffe SJ, et al.** The disposition of indomethacin in preterm babies. *J Pediatr* 1980 | [7441407](https://pubmed.ncbi.nlm.nih.gov/7441407/) | 미숙아 인도메타신 처분 — `CL_IND0` · `V1_IND`. |
| 34 | **Allegaert K, et al.** Intravenous paracetamol (propacetamol) pharmacokinetics in term and preterm neonates. *Eur J Clin Pharmacol* 2004 | [15071761](https://pubmed.ncbi.nlm.nih.gov/15071761/) | 신생아 정맥 파라세타몰 PK — `CL_APAP0` · `V1_APAP`. |
| 35 | **Padavia F, et al.** Population Pharmacokinetics of Intravenous Paracetamol and Its Metabolites in Extreme Preterm Neonates. *Clin Pharmacokinet* 2024 | [39578300](https://pubmed.ncbi.nlm.nih.gov/39578300/) | 초극소 미숙아 파라세타몰·대사체 popPK — 글루쿠로나이드/설페이트 및 CYP 산화 분율(`FNAPQI`). |
| 36 | **Pacifici GM.** Clinical pharmacology of ibuprofen and indomethacin in preterm infants with patent ductus arteriosus. *Curr Pediatr Rev* 2014 | [25088343](https://pubmed.ncbi.nlm.nih.gov/25088343/) | 두 NSAID의 임상약리 비교. |
| 37 | **Thibaut C, et al.** Effect of ibuprofen on bilirubin-albumin binding during the treatment of patent ductus arteriosus in preterm infants. *J Matern Fetal Neonatal Med* 2011 | [21815876](https://pubmed.ncbi.nlm.nih.gov/21815876/) | 이부프로펜의 빌리루빈-알부민 결합 치환 — `BILIDISP_IBU` · `BFREE` 상태변수. |
| 38 | **Diot C, et al.** Effect of ibuprofen on bilirubin-albumin binding in vitro at concentrations observed during treatment of patent ductus arteriosus. *Early Hum Dev* 2010 | [20472375](https://pubmed.ncbi.nlm.nih.gov/20472375/) | 치료 농도에서의 in vitro 치환 정도. |

---

## 6. 아세트아미노펜 — 왜 "다른 부위"인가 (모델의 검증 가능한 예측)

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 39 | **Ouellet M, Percival MD.** Mechanism of acetaminophen inhibition of cyclooxygenase isoforms. *Arch Biochem Biophys* 2001 | [11370851](https://pubmed.ncbi.nlm.nih.gov/11370851/) | **모델 예측의 근거.** 아세트아미노펜은 아라키돈산 채널이 아니라 **peroxidase 부위**에서 작용하며, 따라서 **peroxide와 경쟁적**입니다. `IC50_APAP × (1 + PEROX/KPEROX)` 형태의 직접 근거. |
| 40 | **Graham GG, Scott KF.** Mechanism of action of paracetamol. *Am J Ther* 2005 | [15662292](https://pubmed.ncbi.nlm.nih.gov/15662292/) | 기전 리뷰 — peroxide 긴장도가 높은 조직에서 효력이 떨어진다는 개념. |
| 41 | **Ayoub SS.** Paracetamol (acetaminophen): A familiar drug with an unexplained mechanism of action. *Temperature (Austin)* 2021 | [34901318](https://pubmed.ncbi.nlm.nih.gov/34901318/) | 기전 논쟁의 최신 정리 — 모델의 peroxide 항이 가설임을 명시하기 위해 인용. |
| 42 | **Hammerman C, et al.** Ductal closure with paracetamol: a surprising new approach to patent ductus arteriosus treatment. *Pediatrics* 2011 | [22065264](https://pubmed.ncbi.nlm.nih.gov/22065264/) | 파라세타몰로 관을 닫은 첫 보고. |
| 43 | **Oncel MY, et al.** Oral paracetamol versus oral ibuprofen in the management of patent ductus arteriosus in preterm infants: a randomized controlled trial. *J Pediatr* 2014 | [24359938](https://pubmed.ncbi.nlm.nih.gov/24359938/) | 무작위 비교 — 폐쇄율 유사. |
| 44 | **Ohlsson A, Shah PS.** Paracetamol (acetaminophen) for patent ductus arteriosus in preterm or low-birth-weight infants. *Cochrane Database Syst Rev* 2015 | [25758061](https://pubmed.ncbi.nlm.nih.gov/25758061/) | Cochrane 통합. |
| 45 | **Pranata R, et al.** The efficacy and safety of oral paracetamol versus oral ibuprofen for patent ductus arteriosus closure in preterm neonates. *Indian Heart J* 2020 | [32768013](https://pubmed.ncbi.nlm.nih.gov/32768013/) | 메타분석 — 효능 유사·안전성 우위. |

---

## 7. 장기별 독성 — 같은 효소, 다른 장기

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 46 | **Ohlsson A, et al.** Ibuprofen for the treatment of patent ductus arteriosus in preterm or low birth weight (or both) infants. *Cochrane Database Syst Rev* 2020 | [32045960](https://pubmed.ncbi.nlm.nih.gov/32045960/) | 이부프로펜 대 인도메타신: 폐쇄율 유사, 이부프로펜에서 NEC·일시적 신기능저하 감소. `KI_*_K`·`KI_*_G` 및 `VC_IND` 대 `VC_IBU` 대비의 근거. |
| 47 | **Pacifici GM.** Differential renal adverse effects of ibuprofen and indomethacin in preterm infants: a review. *Clin Pharmacol* 2014 | [25114597](https://pubmed.ncbi.nlm.nih.gov/25114597/) | 두 NSAID의 신장 부작용 차이 — `GFR_PGE2` · `SCR` · `UO` 출력의 대조 근거. |
| 48 | **Patel J, et al.** Randomized double-blind controlled trial comparing the effects of ibuprofen with indomethacin on cerebral hemodynamics in preterm infants with patent ductus arteriosus. *Pediatr Res* 2000 | [10625080](https://pubmed.ncbi.nlm.nih.gov/10625080/) | **`VC_IND` 대 `VC_IBU` 의 직접 근거.** 인도메타신은 뇌 혈역학을 급성으로 저하시키고 이부프로펜은 그렇지 않음. |
| 49 | **Mosca F, et al.** Comparative evaluation of the effects of indomethacin and ibuprofen on cerebral perfusion and oxygenation in preterm infants with patent ductus arteriosus. *J Pediatr* 1997 | [9386657](https://pubmed.ncbi.nlm.nih.gov/9386657/) | 뇌 관류·산소화 비교. |
| 50 | **Pezzati M, et al.** Effects of indomethacin and ibuprofen on mesenteric and renal blood flow in preterm infants with patent ductus arteriosus. *J Pediatr* 1999 | [10586177](https://pubmed.ncbi.nlm.nih.gov/10586177/) | 장·신장 혈류 — `QMESREL` · `QRENREL` 의 근거. |
| 51 | **Watterberg KL, et al.** Prophylaxis of early adrenal insufficiency to prevent bronchopulmonary dysplasia: a multicenter trial. *Pediatrics* 2004 | [15574629](https://pubmed.ncbi.nlm.nih.gov/15574629/) | 조기 하이드로코르티손 + 인도메타신 병용 시 자연 장천공 증가 — `B_SIP_HC` 및 시나리오 S13. |
| 52 | **Ment LR, et al.** Randomized low-dose indomethacin trial for prevention of intraventricular hemorrhage in very low birth weight neonates. *J Pediatr* 1988 | [3373405](https://pubmed.ncbi.nlm.nih.gov/3373405/) | 저용량 인도메타신의 IVH 예방 — `IVH_IND_PROT` (배아기질 혈관 성숙) 항의 근거. |
| 53 | **Sallmon H, et al.** Mature and immature platelets during the first week after birth and incidence of patent ductus arteriosus. *Cardiol Young* 2020 | [32340633](https://pubmed.ncbi.nlm.nih.gov/32340633/) | 혈소판과 PDA — `TXA2` · `BTIME` 경로. |

---

## 8. 혈역학 — 단락의 크기와 그 결과

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 54 | **Lindner W, et al.** Stroke volume and left ventricular output in preterm infants with patent ductus arteriosus. *Pediatr Res* 1990 | [2320395](https://pubmed.ncbi.nlm.nih.gov/2320395/) | PDA에서 좌심실 박출 — `QMAX_LV` · `QP` 스케일. |
| 55 | **Evans N, Kluckow M.** Assessment of ductus arteriosus shunt in preterm infants supported by mechanical ventilation: effect of interatrial shunting. *J Pediatr* 1994 | [7965434](https://pubmed.ncbi.nlm.nih.gov/7965434/) | 단락 평가 방법론. |
| 56 | **Fajardo MF, et al.** Effect of positive end-expiratory pressure on ductal shunting and systemic blood flow in preterm infants with patent ductus arteriosus. *Neonatology* 2014 | [24193163](https://pubmed.ncbi.nlm.nih.gov/24193163/) | PEEP가 단락과 전신 혈류에 미치는 영향 — 지도의 `PEEPADJ` 노드. |
| 57 | **Paradisis M, et al.** Randomized trial of milrinone versus placebo for prevention of low systemic blood flow in very preterm infants. *J Pediatr* 2009 | [18822428](https://pubmed.ncbi.nlm.nih.gov/18822428/) | 낮은 전신 혈류 — `QSDEF` (전신 관류 결핍)의 임상적 실재. |
| 58 | **Zonnenberg I, de Waal K.** The definition of a haemodynamic significant duct in randomized controlled trials: a systematic literature review. *Acta Paediatr* 2012 | [21913976](https://pubmed.ncbi.nlm.nih.gov/21913976/) | RCT들이 "혈역학적으로 유의한 관"을 서로 다르게 정의해 왔다는 문헌고찰 — `QSIG` 문턱을 상태변수로 둔 이유. |
| 59 | **El-Khuffash A, et al.** A Patent Ductus Arteriosus Severity Score Predicts Chronic Lung Disease or Death before Discharge. *J Pediatr* 2015 | [26474706](https://pubmed.ncbi.nlm.nih.gov/26474706/) | PDA 중증도 점수가 결과를 예측 — `PDABUR` (부담)과 표적 치료 전략(S16)의 근거. |
| 60 | **Groves AM, et al.** Introduction to neonatologist-performed echocardiography. *Pediatr Res* 2018 | [30072808](https://pubmed.ncbi.nlm.nih.gov/30072808/) | 신생아 심초음파 — 지도의 `DALPA` · `LAAO` · `ABSDIAST` 측정 노드. |
| 61 | **Levy PT, et al.** Application of Neonatologist Performed Echocardiography in the Assessment and Management of Neonatal Heart Disease. *Pediatr Res* 2018 | [30072802](https://pubmed.ncbi.nlm.nih.gov/30072802/) | 임상 적용. |
| 62 | **MacLellan A, et al.** Fluid restriction for treatment of symptomatic patent ductus arteriosus in preterm infants. *Cochrane Database Syst Rev* 2024 | [39692231](https://pubmed.ncbi.nlm.nih.gov/39692231/) | 수분 제한 — 지도의 `FLUIDR` 노드. |

---

## 9. PDA 부담과 결과 — BPD·NEC·IVH

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 63 | **Clyman RI, et al.** Patent ductus arteriosus, tracheal ventilation, and the risk of bronchopulmonary dysplasia. *Pediatr Res* 2022 | [33790415](https://pubmed.ncbi.nlm.nih.gov/33790415/) | **`PDABUR` × `VENTIX` 곱 구조의 근거.** BPD 위험은 관의 존재 여부가 아니라 노출 기간과 인공환기의 상호작용. |
| 64 | **Nawaytou H, et al.** Patent ductus arteriosus and the risk of bronchopulmonary dysplasia-associated pulmonary hypertension. *Pediatr Res* 2023 | [36804505](https://pubmed.ncbi.nlm.nih.gov/36804505/) | PDA와 BPD 관련 폐고혈압. |
| 65 | **Hamrick SEG, et al.** Patent Ductus Arteriosus of the Preterm Infant. *Pediatrics* 2020 | [33093140](https://pubmed.ncbi.nlm.nih.gov/33093140/) | 종합 리뷰 — 병태생리부터 관리 논쟁까지. |
| 66 | **Benitz WE; AAP Committee on Fetus and Newborn.** Patent Ductus Arteriosus in Preterm Infants. *Pediatrics* 2016 | [26672023](https://pubmed.ncbi.nlm.nih.gov/26672023/) | 미국소아과학회 임상 보고 — "치료가 결과를 개선한다는 근거가 없다"는 입장의 표준 인용. |
| 67 | **Härkin P, et al.** Morbidities associated with patent ductus arteriosus in preterm infants. Nationwide cohort study. *J Matern Fetal Neonatal Med* 2018 | [28651469](https://pubmed.ncbi.nlm.nih.gov/28651469/) | 전국 코호트에서의 동반 질환 — `H_NEC0` · `H_IVH0` 기저 위험 규모. |

---

## 10. 침습적 치료 — 결찰과 경피적 폐쇄

| # | 문헌 | PMID | 내용 |
|---|------|------|------|
| 68 | **Weisz DE, et al.** Association of Patent Ductus Arteriosus Ligation With Death or Neurodevelopmental Impairment Among Extremely Preterm Infants. *JAMA Pediatr* 2017 | [28264088](https://pubmed.ncbi.nlm.nih.gov/28264088/) | 결찰과 사망/신경발달장애의 연관. |
| 69 | **Sathanandam S, et al.** Consensus Guidelines for the Prevention and Management of Periprocedural Complications of Transcatheter Patent Ductus Arteriosus Closure in Extremely Low Birth Weight Infants. *Pediatr Cardiol* 2021 | [34195869](https://pubmed.ncbi.nlm.nih.gov/34195869/) | 경피적 폐쇄 합의 지침. |
| 70 | **Baruteau AE, et al.** The Transcatheter Closure of Patent Ductus Arteriosus in Extremely Low-Birth-Weight Infants: Technique and Outcomes. *J Cardiovasc Dev Dis* 2023 | [38132644](https://pubmed.ncbi.nlm.nih.gov/38132644/) | 술기와 결과. |
| 71 | **Morray BH, et al.** 3-year follow-up of a prospective, multicenter study of the Amplatzer Piccolo Occluder for transcatheter patent ductus arteriosus closure. *J Perinatol* 2023 | [37587183](https://pubmed.ncbi.nlm.nih.gov/37587183/) | Piccolo 3년 추적. |
| 72 | **Duboue PM, et al.** Post-ligation cardiac syndrome after surgical versus transcatheter closure of patent ductus arteriosus. *Eur J Pediatr* 2024 | [38381375](https://pubmed.ncbi.nlm.nih.gov/38381375/) | 결찰 후 심장 증후군 — 지도의 `POSTLIG` 노드. |
| 73 | **Serrano RM, et al.** Comparison of 'post-patent ductus arteriosus ligation syndrome' in premature infants after surgical versus percutaneous closure. *J Perinatol* 2020 | [31578421](https://pubmed.ncbi.nlm.nih.gov/31578421/) | 두 방법 간 비교. |
| 74 | **Mosalli R, Alfaleh K.** Prophylactic surgical ligation of patent ductus arteriosus for prevention of mortality and morbidity in extremely low birth weight infants. *Cochrane Database Syst Rev* 2008 | [18254095](https://pubmed.ncbi.nlm.nih.gov/18254095/) | 예방적 결찰 Cochrane. |
| 75 | **Liebowitz M, Clyman RI.** Prophylactic Indomethacin Compared with Delayed Conservative Management of the Patent Ductus Arteriosus in Extremely Preterm Infants. *J Pediatr* 2017 | [28396025](https://pubmed.ncbi.nlm.nih.gov/28396025/) | 예방적 인도메타신 대 지연 보존 — 시나리오 S7 대 S1. |

---

## 11. 모델링 도구

- **mrgsolve** (R, ODE 기반 PK/PD): <https://mrgsolve.org/>
- **mrgsolve를 이용한 R 기반 QSP**: <https://vantage-research.net/qsp-in-r/>
- **gPKPDviz** — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구:
  논문 <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · 코드 <https://github.com/Genentech/gPKPDviz/>
- **Graphviz** (기계론적 지도 렌더링): <https://graphviz.org/>

---

## ⚠️ 면책 조항

본 모델과 문헌 정리는 **교육 및 연구 목적**입니다. 파라미터는 공개 문헌에서
취했거나 위 §1의 시험 결과에 적합시킨 근사치이며, 독립적으로 검증·인증되지
않았습니다. **실제 임상 의사결정, 처방, 규제 제출에 사용해서는 안 됩니다.**
특히 미숙아 동맥관 개존증은 "치료해야 하는가" 자체가 미해결 문제이며, 이
모델은 그 불확실성을 해소하는 것이 아니라 **정량적으로 서술**하기 위한
도구입니다.
