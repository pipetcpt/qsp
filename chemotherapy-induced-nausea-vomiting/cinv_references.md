# 항암화학요법 유발 구역·구토 (CINV) — 참고문헌
### Chemotherapy-Induced Nausea and Vomiting · QSP model references

모든 PMID는 NCBI E-utilities (`esearch` + `esummary`)로 **실제 조회하여 검증**했습니다.
서지사항(저자·연도·저널·제목)은 PubMed 레코드에서 그대로 가져온 것입니다.
검증 스크립트가 반환하지 않은 인용은 이 목록에 넣지 않았습니다.

> **표기 규칙.** ★ = 모델 보정(calibration)에 **직접 사용된 앵커**.
> ⚑ = 모델의 가정을 **반박하거나 제한하는** 문헌 (일부러 포함시켰습니다).
> ⊕ = 모델의 held-out 예측을 검증하는 데 사용한 문헌.

---

## 1. 보정 앵커 (Calibration anchors) — 이 9개 숫자만 적합에 사용

이 모델의 구토(emesis) 앵커는 **단 하나의 무작위배정 시험**에서 가져왔습니다.
편의에 따라 여러 코호트에서 숫자를 골라 조립하지 않기 위한 의도적 제약입니다.

1. ★ Hesketh PJ, et al. **The oral neurokinin-1 antagonist aprepitant for the prevention of chemotherapy-induced nausea and vomiting: a multinational, randomized, double-blind, placebo-controlled trial in patients receiving high-dose cisplatin.** *J Clin Oncol* 2003. — [PMID 14559886](https://pubmed.ncbi.nlm.nih.gov/14559886/)
   *두 군 × 세 구간 = 6개 앵커 (CR 급성/지연/전체: 대조군 78.1 / 55.8 / 52.3%, 아프레피탄트군 89.2 / 75.4 / 72.7%).*
2. ★ Navari RM, et al. **Olanzapine for Chemotherapy-Induced Nausea and Vomiting.** *N Engl J Med* 2016. — [PMID 27705262](https://pubmed.ncbi.nlm.nih.gov/27705262/)
   ***대조군만** 사용 (구역 없음 0–120 h = 22%)하여 구역 척도를 고정. 올란자핀군(37%)은 held-out.*
3. ★ Cubeddu LX. **Serotonin mechanisms in chemotherapy-induced emesis in cancer patients.** *Oncology* 1996. — [PMID 8692546](https://pubmed.ncbi.nlm.nih.gov/8692546/)
   *미치료 시스플라틴의 자연경과(24시간 구토 횟수 규모)에 대한 근거.*

---

## 2. Held-out 검증에 사용한 임상시험 (적합에 사용하지 않음)

4. ⊕ Poli-Bigelli S, et al. **Addition of the neurokinin 1 receptor antagonist aprepitant to standard antiemetic therapy improves control of chemotherapy-induced nausea and vomiting.** *Cancer* 2003. — [PMID 12784346](https://pubmed.ncbi.nlm.nih.gov/12784346/)
5. ⊕ Warr DG, et al. **The oral NK(1) antagonist aprepitant for the prevention of acute and delayed chemotherapy-induced nausea and vomiting.** *Eur J Cancer* 2005. — [PMID 15939263](https://pubmed.ncbi.nlm.nih.gov/15939263/)
6. ⊕ Herrstedt J, et al. **Efficacy and tolerability of aprepitant for the prevention of chemotherapy-induced nausea and emesis over multiple cycles of moderately emetogenic chemotherapy.** *Cancer* 2005. — [PMID 16104039](https://pubmed.ncbi.nlm.nih.gov/16104039/)
7. ⊕ Gralla R, et al. **Palonosetron improves prevention of chemotherapy-induced nausea and vomiting following moderately emetogenic chemotherapy.** *Ann Oncol* 2003. — [PMID 14504060](https://pubmed.ncbi.nlm.nih.gov/14504060/)
8. ⊕ Hesketh PJ, et al. **Efficacy and safety of NEPA, an oral combination of netupitant and palonosetron, for prevention of chemotherapy-induced nausea and vomiting following highly emetogenic chemotherapy.** *Ann Oncol* 2014. — [PMID 24608196](https://pubmed.ncbi.nlm.nih.gov/24608196/)
9. ⊕ Rapoport BL, et al. **Safety and efficacy of rolapitant for prevention of chemotherapy-induced nausea and vomiting after administration of cisplatin-based highly emetogenic chemotherapy.** *Lancet Oncol* 2015. — [PMID 26272769](https://pubmed.ncbi.nlm.nih.gov/26272769/)
10. ⊕ Schwartzberg LS, et al. **Safety and efficacy of rolapitant for prevention of chemotherapy-induced nausea and vomiting after administration of moderately emetogenic chemotherapy.** *Lancet Oncol* 2015. — [PMID 26272768](https://pubmed.ncbi.nlm.nih.gov/26272768/)
11. ⊕ Hashimoto H, et al. **Olanzapine 5 mg plus standard antiemetic therapy for the prevention of chemotherapy-induced nausea and vomiting (J-FORCE): a multicentre, randomised, double-blind, placebo-controlled, phase 3 trial.** *Lancet Oncol* 2020. — [PMID 31838011](https://pubmed.ncbi.nlm.nih.gov/31838011/)
12. Hashimoto H, et al. **Study protocol for J-SUPPORT 1604 (J-FORCE).** *Jpn J Clin Oncol* 2018. — [PMID 30124989](https://pubmed.ncbi.nlm.nih.gov/30124989/)
13. ⊕ Zhao Y, et al. **A multicenter, randomized, double-blind, placebo-controlled, phase 3 trial of olanzapine plus triple antiemetic therapy.** *EClinicalMedicine* 2023. — [PMID 36712888](https://pubmed.ncbi.nlm.nih.gov/36712888/)
14. ⊕ Boccia R, et al. **Randomized phase III trial of APF530 versus palonosetron in the prevention of chemotherapy-induced nausea and vomiting.** *BMC Cancer* 2016. — [PMID 26921245](https://pubmed.ncbi.nlm.nih.gov/26921245/)
15. ⊕ Schnadig ID, et al. **APF530 versus ondansetron, each in a guideline-recommended three-drug regimen.** *Cancer Manag Res* 2017. — [PMID 28579832](https://pubmed.ncbi.nlm.nih.gov/28579832/)
16. ⊕ Weinstein C, et al. **Single-dose fosaprepitant for the prevention of chemotherapy-induced nausea and vomiting.** *BMC Cancer* 2020. — [PMID 32988373](https://pubmed.ncbi.nlm.nih.gov/32988373/)
17. ⊕ Yahata H, et al. **Efficacy of aprepitant for the prevention of chemotherapy-induced nausea and vomiting with a moderately emetogenic chemotherapy regimen.** *Int J Clin Oncol* 2016. — [PMID 26662632](https://pubmed.ncbi.nlm.nih.gov/26662632/)
18. Longo F, et al. **Combination of aprepitant, palonosetron and dexamethasone as antiemetic prophylaxis in lung cancer patients receiving multiple cycles of cisplatin.** *Int J Clin Pract* 2012. — [PMID 22805267](https://pubmed.ncbi.nlm.nih.gov/22805267/)
19. Kubota K, et al. **Control of nausea with palonosetron versus granisetron, both combined with dexamethasone.** *Support Care Cancer* 2016. — [PMID 27129842](https://pubmed.ncbi.nlm.nih.gov/27129842/)
20. Chow R, et al. **Efficacy of the combination neurokinin-1 receptor antagonist, palonosetron, and dexamethasone.** *Ann Palliat Med* 2018. — [PMID 29764184](https://pubmed.ncbi.nlm.nih.gov/29764184/)
21. Iimura Y, et al. **A randomized, double-blind, placebo-controlled phase III study evaluating the preventive effect.** *Int J Clin Oncol* 2025. — [PMID 40369354](https://pubmed.ncbi.nlm.nih.gov/40369354/)
22. Zhang M, et al. **Efficacy and safety of Aprepitant-containing triple therapy for the prevention and treatment of CINV.** *Medicine (Baltimore)* 2023. — [PMID 38013306](https://pubmed.ncbi.nlm.nih.gov/38013306/)
23. Ünlü B, et al. **Aprepitant and Fosaprepitant for Preventing Nausea and Vomiting in Patients Receiving Highly Emetogenic Chemotherapy.** *Curr Oncol* 2026. — [PMID 42505183](https://pubmed.ncbi.nlm.nih.gov/42505183/)

---

## 3. 메타분석 · 계통적 문헌고찰

24. Piechotta V, et al. **Antiemetics for adults for prevention of nausea and vomiting caused by moderately or highly emetogenic chemotherapy: a network meta-analysis.** *Cochrane Database Syst Rev* 2021. — [PMID 34784425](https://pubmed.ncbi.nlm.nih.gov/34784425/)
25. Sutherland A, et al. **Olanzapine for the prevention and treatment of cancer-related nausea and vomiting in adults.** *Cochrane Database Syst Rev* 2018. — [PMID 30246876](https://pubmed.ncbi.nlm.nih.gov/30246876/)
26. Smith LA, et al. **Cannabinoids for nausea and vomiting in adults with cancer receiving chemotherapy.** *Cochrane Database Syst Rev* 2015. — [PMID 26561338](https://pubmed.ncbi.nlm.nih.gov/26561338/)
27. Di Maio M, et al. **Efficacy of neurokinin-1 receptor antagonists in the prevention of chemotherapy-induced nausea and vomiting.** *Crit Rev Oncol Hematol* 2018. — [PMID 29548482](https://pubmed.ncbi.nlm.nih.gov/29548482/)
28. Tricco AC, et al. **Comparative safety and effectiveness of serotonin receptor antagonists in patients undergoing chemotherapy.** *BMC Med* 2016. — [PMID 28007031](https://pubmed.ncbi.nlm.nih.gov/28007031/)
29. ⊕ Chow R, et al. **Olanzapine (10 mg vs 5 mg vs 2.5 mg) for the prophylaxis of chemotherapy-induced nausea and vomiting.** *Support Care Cancer* 2026. — [PMID 42461312](https://pubmed.ncbi.nlm.nih.gov/42461312/)
30. Chow R, et al. **Multi-day vs single-day dexamethasone for the prophylaxis of chemotherapy-induced nausea and vomiting.** *Support Care Cancer* 2024. — [PMID 39432169](https://pubmed.ncbi.nlm.nih.gov/39432169/)
31. Nimmons D, et al. **Clinical effectiveness of pharmacological and non-pharmacological treatments for the management of nausea.** *Neurosci Biobehav Rev* 2024. — [PMID 38097097](https://pubmed.ncbi.nlm.nih.gov/38097097/)
32. Ezzo J, et al. **Acupuncture-point stimulation for chemotherapy-induced nausea or vomiting.** *Cochrane Database Syst Rev* 2014 (withdrawn). — [PMID 25412832](https://pubmed.ncbi.nlm.nih.gov/25412832/)

---

## 4. 급성기 기전 — 장내 엔테로크로마핀 세포와 세로토닌 (모델의 ACUTE 항)

33. ★ Cubeddu LX. **Serotonin mechanisms in chemotherapy-induced emesis in cancer patients.** *Oncology* 1996. — [PMID 8692546](https://pubmed.ncbi.nlm.nih.gov/8692546/)
34. ⊕ **Higa GM, et al. 5-Hydroxyindoleacetic acid and substance P profiles in patients receiving emetogenic chemotherapy.** *J Oncol Pharm Pract* 2006. — [PMID 17156592](https://pubmed.ncbi.nlm.nih.gov/17156592/)
    *이 모델의 급성(5-HIAA burst)/지연(substance P ramp) 시간 분리를 검증하는 데 쓴 핵심 문헌.*
35. Castejon AM, et al. **Use of intravenous microdialysis to monitor changes in serotonin release and metabolism induced by cisplatin.** *J Pharmacol Exp Ther* 1999. — [PMID 10565811](https://pubmed.ncbi.nlm.nih.gov/10565811/)
36. Nemani MR, et al. **The mucosal-neural axis in chemotherapy-induced gastrointestinal toxicity.** *Oncol Rev* 2026. — [PMID 42199990](https://pubmed.ncbi.nlm.nih.gov/42199990/)
37. Andrews PL, Sanger GJ. **Signals for nausea and emesis: Implications for models of upper gastrointestinal diseases.** *Auton Neurosci* 2006. — [PMID 16556512](https://pubmed.ncbi.nlm.nih.gov/16556512/)
38. Ullah I, et al. **Phytotherapeutic Approach in the Management of Cisplatin Induced Vomiting; Neurochemical Considerations.** *Oxid Med Cell Longev* 2022. — [PMID 36148411](https://pubmed.ncbi.nlm.nih.gov/36148411/)
39. Yuan W, et al. **Prevention of cisplatin-induced nausea and vomiting by seabuckthorn.** *Iran J Basic Med Sci* 2021. — [PMID 33953865](https://pubmed.ncbi.nlm.nih.gov/33953865/)

## 5. 5-HT3 수용체 — 구조, 결합동역학, 내재화

40. Ye ZJ, et al. **The serotonin-gated 5-HT3 receptor: a tale of functions and structures of a prototypical pentameric ligand-gated ion channel.** *Acta Pharmacol Sin* 2026. — [PMID 41813973](https://pubmed.ncbi.nlm.nih.gov/41813973/)
41. Hope AG, et al. **Characterization of a human 5-hydroxytryptamine3 receptor type A subunit stably expressed in a cell line.** *Br J Pharmacol* 1996. — [PMID 8818349](https://pubmed.ncbi.nlm.nih.gov/8818349/)
42. **Rojas C, et al. Palonosetron exhibits unique molecular interactions with the 5-HT3 receptor.** *Anesth Analg* 2008. — [PMID 18633025](https://pubmed.ncbi.nlm.nih.gov/18633025/)
    *팔로노세트론의 알로스테릭 결합 — 모델에서 느린 `KOFFP`와 `KINT`(수용체 내재화)의 근거.*
43. **Rojas C, et al. The antiemetic 5-HT3 receptor antagonist Palonosetron inhibits substance P-mediated responses in vitro and in vivo.** *J Pharmacol Exp Ther* 2010. — [PMID 20724484](https://pubmed.ncbi.nlm.nih.gov/20724484/)
    *모델의 `XTALK` 항(팔로노세트론의 부분적 지연기 커버)의 근거.*
44. Rojas C, Slusher BS. **Pharmacological mechanisms of 5-HT₃ and tachykinin NK₁ receptor antagonism to prevent chemotherapy-induced nausea and vomiting.** *Eur J Pharmacol* 2012. — [PMID 22425650](https://pubmed.ncbi.nlm.nih.gov/22425650/)
45. Panuganti VK, et al. **A novel Ondansetron extended release injectable suspension (OERIS).** *Cancer Chemother Pharmacol* 2025. — [PMID 41284044](https://pubmed.ncbi.nlm.nih.gov/41284044/)

## 6. 지연기 기전 — Substance P / NK1 (모델의 DELAYED 항)

46. **Hesketh PJ, et al. Differential time course of action of 5-HT3 and NK1 receptor antagonists when used with highly emetogenic chemotherapy.** *Support Care Cancer* 2011. — [PMID 20623144](https://pubmed.ncbi.nlm.nih.gov/20623144/)
    *이 모델의 조직 원리 #2와 **독립적으로 동일한 결론**에 도달한 임상 관찰.*
47. Andrews PLR, Sanger GJ. **An assessment of the effects of neurokinin(1) receptor antagonism against nausea and vomiting.** *Br J Clin Pharmacol* 2023. — [PMID 37452618](https://pubmed.ncbi.nlm.nih.gov/37452618/)
48. Minami M, et al. **Antiemetic effects of sendide, a peptide tachykinin NK1 receptor antagonist, in the ferret.** *Eur J Pharmacol* 1998. — [PMID 9877081](https://pubmed.ncbi.nlm.nih.gov/9877081/)
49. Andrews PL, Bhandari P. **The anti-emetic action of the neurokinin(1) receptor antagonist CP-99,994 does not require the presence of the area postrema.** *Neurosci Lett* 2001. — [PMID 11698156](https://pubmed.ncbi.nlm.nih.gov/11698156/)
50. Ariumi H, et al. **The role of tachykinin NK-1 receptors in the area postrema of ferrets in emesis.** *Neurosci Lett* 2000. — [PMID 10825652](https://pubmed.ncbi.nlm.nih.gov/10825652/)
51. Rudd JA, et al. **Profile of Antiemetic Activity of Netupitant Alone or in Combination with Palonosetron and Dexamethasone in Ferrets and Suncus murinus.** *Front Pharmacol* 2016. — [PMID 27630563](https://pubmed.ncbi.nlm.nih.gov/27630563/)
52. Hargreaves R, et al. **Development of aprepitant, the first neurokinin-1 receptor antagonist for the prevention of chemotherapy-induced nausea and vomiting.** *Ann N Y Acad Sci* 2011. — [PMID 21434941](https://pubmed.ncbi.nlm.nih.gov/21434941/)
53. Di Fabio R, et al. **Discovery and biological characterization of a potent NK1 receptor antagonist.** *J Med Chem* 2011. — [PMID 21229983](https://pubmed.ncbi.nlm.nih.gov/21229983/)
54. Ruzza C, et al. **In vitro and in vivo pharmacological characterization of Pronetupitant, a prodrug of the neurokinin 1 receptor antagonist netupitant.** *Peptides* 2015. — [PMID 25843024](https://pubmed.ncbi.nlm.nih.gov/25843024/)
55. Walter T, et al. **Broadening the View: Substance P and Its Metabolism in Pruritus-Related Diseases.** *J Dermatol* 2026. — [PMID 42244108](https://pubmed.ncbi.nlm.nih.gov/42244108/)
56. Wu Z, et al. **The role of neprilysin in musculoskeletal diseases.** *Tissue Cell* 2026. — [PMID 41529369](https://pubmed.ncbi.nlm.nih.gov/41529369/)
57. Li WW, et al. **Neuropeptide regulation of adaptive immunity.** *J Neuroinflammation* 2018. — [PMID 29642930](https://pubmed.ncbi.nlm.nih.gov/29642930/)

## 7. NK1 수용체 점유율 (held-out 검증)

58. ⊕ **Bergström M, et al. Human positron emission tomography studies of brain neurokinin 1 receptor occupancy by aprepitant.** *Biol Psychiatry* 2004. — [PMID 15121485](https://pubmed.ncbi.nlm.nih.gov/15121485/)
    *모델이 예측한 24시간 지속 NK1 점유율의 held-out 검증에 사용.*
59. Hart XM, et al. **Update Lessons from Positron Emission Tomography Imaging Part I.** *Ther Drug Monit* 2024. — [PMID 38018857](https://pubmed.ncbi.nlm.nih.gov/38018857/)

## 8. 중추 회로 — 최후영역(area postrema), NTS, 섬피질

60. Borison HL. **Area postrema: chemoreceptor circumventricular organ of the medulla oblongata.** *Prog Neurobiol* 1989. — [PMID 2660187](https://pubmed.ncbi.nlm.nih.gov/2660187/)
    *혈액뇌장벽이 없는 최후영역 — 모델의 `WAP`(수용체와 무관한 항)의 해부학적 근거.*
61. Miller AD, Leslie RA. **The area postrema and vomiting.** *Front Neuroendocrinol* 1994. — [PMID 7895890](https://pubmed.ncbi.nlm.nih.gov/7895890/)
62. Hornby PJ. **Central neurocircuitry associated with emesis.** *Am J Med* 2001. — [PMID 11749934](https://pubmed.ncbi.nlm.nih.gov/11749934/)
63. Sanger GJ, Andrews PL. **Treatment of nausea and vomiting: gaps in our knowledge.** *Auton Neurosci* 2006. — [PMID 16934536](https://pubmed.ncbi.nlm.nih.gov/16934536/)
64. Parvizi J, et al. **Functional Architecture of the Human Insula Revealed by Causal Intracranial Mapping.** *Res Sq* 2026. — [PMID 41836516](https://pubmed.ncbi.nlm.nih.gov/41836516/)
65. Hagiwara K, et al. **Insular lobe epilepsy. Part 1: semiology.** *Rinsho Shinkeigaku* 2024. — [PMID 39069491](https://pubmed.ncbi.nlm.nih.gov/39069491/)
66. Bolender A, et al. **Altered parasympathetic outflow and central sensitization response to continuous pain in cyclic vomiting.** *Am J Physiol Gastrointest Liver Physiol* 2025. — [PMID 39545877](https://pubmed.ncbi.nlm.nih.gov/39545877/)
67. Belkacemi L, Darmani NA. **Evidence for Bell-Shaped Dose-Response Emetic Effects of Temsirolimus and Analogs.** *Front Pharmacol* 2022. — [PMID 35444553](https://pubmed.ncbi.nlm.nih.gov/35444553/)

## 9. 위 운동성 · 서파 · 바소프레신 (구역 축의 주변부 항)

68. van Helden DF, et al. **Mechanisms underlying gastric antral slow waves — a human perspective.** *J Physiol* 2026. — [PMID 42285761](https://pubmed.ncbi.nlm.nih.gov/42285761/)
69. Bashashati M, et al. **Gastric electrophysiology.** *Expert Rev Gastroenterol Hepatol* 2026. — [PMID 42018791](https://pubmed.ncbi.nlm.nih.gov/42018791/)
70. Feinle C, et al. **Relationship between increasing duodenal lipid doses, gastric perception, and plasma hormone levels.** *Am J Physiol* 2000. — [PMID 10801290](https://pubmed.ncbi.nlm.nih.gov/10801290/)
71. Wellington J, et al. **Effect of endoscopic pyloric therapies for patients with nausea and vomiting and functional obstructive gastroparesis.** *Auton Neurosci* 2017. — [PMID 27460691](https://pubmed.ncbi.nlm.nih.gov/27460691/)
72. Shakil A, et al. **Gastrointestinal complications of diabetes.** *Am Fam Physician* 2008. — [PMID 18619079](https://pubmed.ncbi.nlm.nih.gov/18619079/)

## 10. 예기성 구역 — 연합학습 (모델의 ANTIC 항)

73. Kamen C, et al. **Anticipatory nausea and vomiting due to chemotherapy.** *Eur J Pharmacol* 2014. — [PMID 24157982](https://pubmed.ncbi.nlm.nih.gov/24157982/)
74. Akechi T, et al. **Anticipatory nausea among ambulatory cancer patients undergoing chemotherapy: prevalence, associated factors, and impact on quality of life.** *Cancer Sci* 2010. — [PMID 20946120](https://pubmed.ncbi.nlm.nih.gov/20946120/)
75. Rivi V, et al. **A translational and multidisciplinary approach to studying the Garcia effect, a higher form of learning.** *J Exp Biol* 2024. — [PMID 38639079](https://pubmed.ncbi.nlm.nih.gov/38639079/)
76. Peng TR, et al. **Revisiting the Role of Lorazepam as an Adjunct in the Management of Chemotherapy-Induced Nausea and Vomiting.** *Biomedicines* 2026. — [PMID 42072465](https://pubmed.ncbi.nlm.nih.gov/42072465/)

## 11. 코르티코스테로이드 — 유전체(genomic) 작용의 느린 시간상수

77. Lyu D, et al. **Glucocorticoid signaling mitigates visceral hypersensitivity by suppressing inflammatory signalling.** *Int Immunopharmacol* 2026. — [PMID 41554225](https://pubmed.ncbi.nlm.nih.gov/41554225/)
78. Agar M, et al. **Validating self-report and proxy reports of the Dexamethasone Symptom Questionnaire.** *Support Care Cancer* 2016. — [PMID 26294320](https://pubmed.ncbi.nlm.nih.gov/26294320/)
79. Watanabe D, et al. **One-Day Versus Three-Day Dexamethasone with NK1RA for Patients Receiving Carboplatin and Moderately Emetogenic Chemotherapy.** *Oncologist* 2022. — [PMID 35427418](https://pubmed.ncbi.nlm.nih.gov/35427418/)
80. Celio L, et al. **Dexamethasone-sparing regimens with NEPA for the prevention of chemotherapy-induced nausea and vomiting.** *J Geriatr Oncol* 2023. — [PMID 37290207](https://pubmed.ncbi.nlm.nih.gov/37290207/)
81. Lee J, et al. **Reduced-dose dexamethasone premedication for weekly paclitaxel.** *Support Care Cancer* 2026. — [PMID 42371187](https://pubmed.ncbi.nlm.nih.gov/42371187/)
82. Johnson TN, et al. **Bioavailability of Oral Hydrocortisone Corrected for Binding Proteins.** *J Bioequivalence Bioavailab* 2018. — [PMID 29795974](https://pubmed.ncbi.nlm.nih.gov/29795974/)

## 12. 약물상호작용 (CYP3A4 / CYP2D6) — 덱사메타손 용량 보정의 근거

83. ★⊕ **McCrea JB, et al. Effects of the neurokinin1 receptor antagonist aprepitant on the pharmacokinetics of dexamethasone and methylprednisolone.** *Clin Pharmacol Ther* 2003. — [PMID 12844131](https://pubmed.ncbi.nlm.nih.gov/12844131/)
84. ⊕ Natale JJ, et al. **Drug-drug interaction profile of components of a fixed combination of netupitant and palonosetron.** *J Oncol Pharm Pract* 2016. — [PMID 25998320](https://pubmed.ncbi.nlm.nih.gov/25998320/)
85. Marbury TC, et al. **Pharmacokinetics of oral dexamethasone and midazolam when administered with single-dose intravenous fosaprepitant.** *J Clin Pharmacol* 2011. — [PMID 21209230](https://pubmed.ncbi.nlm.nih.gov/21209230/)
86. ⚑ **Nijstad AL, et al. Overestimation of the effect of (fos)aprepitant on intravenous dexamethasone pharmacokinetics.** *Support Care Cancer* 2022. — [PMID 36287279](https://pubmed.ncbi.nlm.nih.gov/36287279/)
    *이 문헌은 모델이 앵커로 쓴 2.2배 AUC 증가가 **과대추정**일 수 있음을 지적합니다. README의 "보고된 실패" 절에서 이 반박을 명시적으로 다룹니다.*
87. ⊕ Wang J, et al. **Effects of rolapitant administered orally on the pharmacokinetics of dextromethorphan (CYP2D6), tolbutamide (CYP2C9), and efavirenz.** *Support Care Cancer* 2019. — [PMID 30084103](https://pubmed.ncbi.nlm.nih.gov/30084103/)
    *롤라피탄트가 CYP2D6만 억제하고 CYP3A4는 건드리지 않는다는 "계열 분화"의 근거.*
88. Rapoport BL. **Differential pharmacology and clinical utility of rolapitant in chemotherapy-induced nausea and vomiting.** *Cancer Manag Res* 2017. — [PMID 28260945](https://pubmed.ncbi.nlm.nih.gov/28260945/)
89. Deb S, et al. **Simulation of drug-drug interactions between breast cancer chemotherapeutic agents and antiemetics.** *Daru* 2023. — [PMID 37223851](https://pubmed.ncbi.nlm.nih.gov/37223851/)
90. Imbs DC, et al. **Pharmacokinetic interaction between pazopanib and cisplatin regimen.** *Cancer Chemother Pharmacol* 2016. — [PMID 26779916](https://pubmed.ncbi.nlm.nih.gov/26779916/)
91. Shahzad I, et al. **Development and Evaluation of Physiologically Based Pharmacokinetic (PBPK) Models.** *Pharmaceuticals* 2026. — [PMID 42515786](https://pubmed.ncbi.nlm.nih.gov/42515786/)

## 13. 항구토제 약동학 (모델 PK 파라미터의 출처)

92. Jann MW, et al. **Relative bioavailability of ondansetron 8-mg oral tablets versus extemporaneous suppositories.** *Pharmacotherapy* 1998. — [PMID 9545148](https://pubmed.ncbi.nlm.nih.gov/9545148/)
93. Zhang ZY, et al. **Absorption, metabolism, and excretion of the antiemetic rolapitant, a selective neurokinin-1 receptor antagonist.** *Invest New Drugs* 2019. — [PMID 30032410](https://pubmed.ncbi.nlm.nih.gov/30032410/)
94. Kurteva G, et al. **Pharmacokinetic profile and safety of intravenous NEPA, a fixed combination of fosnetupitant and palonosetron.** *Eur J Pharm Sci* 2019. — [PMID 31404621](https://pubmed.ncbi.nlm.nih.gov/31404621/)
95. Gasanin E, et al. **Bioequivalence of netupitant and palonosetron (NEPA) oral suspension and hard capsules.** *Support Care Cancer* 2026. — [PMID 41627515](https://pubmed.ncbi.nlm.nih.gov/41627515/)
96. Callaghan JT, et al. **Olanzapine. Pharmacokinetic and pharmacodynamic profile.** *Clin Pharmacokinet* 1999. — [PMID 10511917](https://pubmed.ncbi.nlm.nih.gov/10511917/)
97. Atanasovska E, et al. **Comparative, Single-Dose, 2-Way Cross-Over Bioavailability Study of Two Olanzapine 10 mg Tablets.** *Pril (Makedon Akad Nauk)* 2022. — [PMID 35843923](https://pubmed.ncbi.nlm.nih.gov/35843923/)
98. Cherniakov I, et al. **Safety, Tolerability, and Pharmacokinetics of Subcutaneous Extended-Release Injectable Olanzapine.** *Clin Drug Investig* 2026. — [PMID 41632432](https://pubmed.ncbi.nlm.nih.gov/41632432/)
99. Greenblatt DJ, et al. **Clinical pharmacokinetics of anxiolytics and hypnotics in the elderly.** *Clin Pharmacokinet* 1991. — [PMID 1684924](https://pubmed.ncbi.nlm.nih.gov/1684924/)
100. Ghiasi N, et al. **Lorazepam.** *StatPearls*. — [PMID 30422485](https://pubmed.ncbi.nlm.nih.gov/30422485/)
101. Jang JH, et al. **Physiologically based pharmacokinetic modeling of granisetron for optimizing antiemetic therapy.** *Cancer Chemother Pharmacol* 2026. — [PMID 42159758](https://pubmed.ncbi.nlm.nih.gov/42159758/)
102. Li Q, et al. **Pharmacokinetics, safety, and population pharmacokinetic profiles of Ritanine® in Chinese subjects.** *Front Pharmacol* 2026. — [PMID 42305433](https://pubmed.ncbi.nlm.nih.gov/42305433/)
103. Eiseman JL, et al. **Plasma pharmacokinetics and tissue and brain distribution of cisplatin in musk shrews.** *Cancer Chemother Pharmacol* 2015. — [PMID 25398697](https://pubmed.ncbi.nlm.nih.gov/25398697/)
104. Tao J, et al. **Systematic Identification of Proteins Binding with Cisplatin in Blood.** *J Proteome Res* 2021. — [PMID 34427088](https://pubmed.ncbi.nlm.nih.gov/34427088/)

## 14. 수용체 결합 프로파일 (올란자핀의 다중 수용체 작용)

105. **Bymaster FP, et al. In vitro and in vivo biochemistry of olanzapine: a novel, atypical antipsychotic drug.** *J Clin Psychiatry* 1997. — [PMID 9265914](https://pubmed.ncbi.nlm.nih.gov/9265914/)
     *모델의 `KI_OLA_D2 / _2A / _H1 / _M1` 네 개 상수의 출처.*
106. Mamo D, et al. **A PET study of dopamine D2 and serotonin 5-HT2 receptor occupancy in patients with schizophrenia.** *Am J Psychiatry* 2004. — [PMID 15121646](https://pubmed.ncbi.nlm.nih.gov/15121646/)
107. Sun Y, et al. **The Cannabinoid CB(1) Receptor Inverse Agonist/Antagonist SR141716A Activates the Adenylate Cyclase pathway.** *Int J Mol Sci* 2025. — [PMID 41155180](https://pubmed.ncbi.nlm.nih.gov/41155180/)

## 15. 약물유전체학

108. Moore C, et al. **Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for CYP2D6 Genotype and Use of Ondansetron and Tropisetron.** *Clin Pharmacol Ther* 2026. — [PMID 41979467](https://pubmed.ncbi.nlm.nih.gov/41979467/)
109. Moore C, et al. **CYP2D6 genotype and associated 5-HT(3) receptor antagonist outcomes: A systematic review and meta-analysis.** *Clin Transl Sci* 2025. — [PMID 39899439](https://pubmed.ncbi.nlm.nih.gov/39899439/)
110. Perwitasari DA, et al. **Association of ABCB1, 5-HT3B receptor and CYP2D6 genetic polymorphisms with ondansetron and metoclopramide antiemetic response.** *Jpn J Clin Oncol* 2011. — [PMID 21840870](https://pubmed.ncbi.nlm.nih.gov/21840870/)
111. Farhat K, et al. **Association of anti-emetic efficacy of Ondansetron with 18792A>G polymorphism in a drug target gene.** *J Pak Med Assoc* 2018. — [PMID 29885172](https://pubmed.ncbi.nlm.nih.gov/29885172/)
112. Jin Y, et al. **A randomized trial evaluating the association between related gene polymorphism and nausea and vomiting.** *BMC Med Genomics* 2023. — [PMID 37924126](https://pubmed.ncbi.nlm.nih.gov/37924126/)

## 16. 위험인자 · 예측모형 (모델에서는 출력이지 규칙이 아님)

113. Dranitsaris G, et al. **The development of a prediction tool to identify cancer patients at high risk for chemotherapy-induced nausea and vomiting.** *Ann Oncol* 2017. — [PMID 28398530](https://pubmed.ncbi.nlm.nih.gov/28398530/)
114. Warr DG, et al. **Evaluation of risk factors predictive of nausea and vomiting with current standard-of-care antiemetic treatment.** *Support Care Cancer* 2011. — [PMID 20461438](https://pubmed.ncbi.nlm.nih.gov/20461438/)
115. Bouganim N, et al. **Prospective validation of risk prediction indexes for acute and delayed chemotherapy-induced nausea and vomiting.** *Curr Oncol* 2012. — [PMID 23300365](https://pubmed.ncbi.nlm.nih.gov/23300365/)
116. Zhao Y, et al. **Validation of different personalized risk models of chemotherapy-induced nausea and vomiting.** *Cancer Commun* 2023. — [PMID 36545810](https://pubmed.ncbi.nlm.nih.gov/36545810/)

## 17. 구토유발성 분류와 진료지침

117. Grunberg SM, et al. **Evaluation of new antiemetic agents and definition of antineoplastic agent emetogenicity — an update.** *Support Care Cancer* 2005. — [PMID 15599601](https://pubmed.ncbi.nlm.nih.gov/15599601/)
118. Grunberg SM, et al. **Evaluation of new antiemetic agents and definition of antineoplastic agent emetogenicity — state of the art.** *Support Care Cancer* 2011. — [PMID 20972805](https://pubmed.ncbi.nlm.nih.gov/20972805/)
119. Kris MG, et al. **Consensus proposals for the prevention of acute and delayed vomiting and nausea following high-emetic-risk chemotherapy.** *Support Care Cancer* 2005. — [PMID 15565277](https://pubmed.ncbi.nlm.nih.gov/15565277/)
120. Herrstedt J, et al. **2016 Updated MASCC/ESMO Consensus Recommendations: Prevention of Nausea and Vomiting Following High Emetic Risk Chemotherapy.** *Support Care Cancer* 2017. — [PMID 27443154](https://pubmed.ncbi.nlm.nih.gov/27443154/)
121. Einhorn LH, et al. **2016 updated MASCC/ESMO consensus recommendations: prevention of nausea and vomiting following multiple-day chemotherapy.** *Support Care Cancer* 2017. — [PMID 27815710](https://pubmed.ncbi.nlm.nih.gov/27815710/)
122. Hesketh PJ, et al. **Antiemetics: ASCO Guideline Update.** *J Clin Oncol* 2020. — [PMID 32658626](https://pubmed.ncbi.nlm.nih.gov/32658626/)
123. Ruhlmann CH, et al. **2016 updated MASCC/ESMO consensus recommendations: prevention of radiotherapy-induced nausea and vomiting.** *Support Care Cancer* 2017. — [PMID 27624464](https://pubmed.ncbi.nlm.nih.gov/27624464/)
124. Feyer PC, et al. **Radiotherapy-induced nausea and vomiting (RINV): MASCC/ESMO guideline.** *Support Care Cancer* 2011. — [PMID 20697746](https://pubmed.ncbi.nlm.nih.gov/20697746/)
125. McKenzie E, et al. **Radiation-induced nausea and vomiting: a comparison between MASCC/ESMO, ASCO, and NCCN antiemetic guidelines.** *Support Care Cancer* 2019. — [PMID 30607675](https://pubmed.ncbi.nlm.nih.gov/30607675/)
126. Gan TJ, et al. **Fourth Consensus Guidelines for the Management of Postoperative Nausea and Vomiting.** *Anesth Analg* 2020. — [PMID 32467512](https://pubmed.ncbi.nlm.nih.gov/32467512/)
127. Navari RM, Aapro M. **Antiemetic Prophylaxis for Chemotherapy-Induced Nausea and Vomiting.** *N Engl J Med* 2016. — [PMID 27050207](https://pubmed.ncbi.nlm.nih.gov/27050207/)
128. Navari RM. **Management of chemotherapy-induced nausea and vomiting: focus on newer agents and new uses for older agents.** *Drugs* 2013. — [PMID 23404093](https://pubmed.ncbi.nlm.nih.gov/23404093/)
129. Grunberg S. **Patient-centered management of chemotherapy-induced nausea and vomiting.** *Cancer Control* 2012. — [PMID 22488023](https://pubmed.ncbi.nlm.nih.gov/22488023/)
130. Ryan MF, et al. **Moving beyond promethazine: Advancing precision in antiemetic therapy.** *SAGE Open Med* 2026. — [PMID 42205501](https://pubmed.ncbi.nlm.nih.gov/42205501/)

## 18. 안전성 — QTc, 추체외로증상, 기타

131. ⊕ **Zuo P, et al. Integration of modeling and simulation to support changes to ondansetron dosing following a randomized, double-blind, placebo-, and active-controlled thorough QT study.** *J Clin Pharmacol* 2014. — [PMID 24782199](https://pubmed.ncbi.nlm.nih.gov/24782199/)
     *32 mg 정맥 용량 철회의 근거이며, 모델의 `SQ_OND` 기울기를 검증하는 held-out 자료.*
132. Obal D, et al. **Perioperative doses of ondansetron or dolasetron do not lengthen the QT interval.** *Mayo Clin Proc* 2014. — [PMID 24388024](https://pubmed.ncbi.nlm.nih.gov/24388024/)
133. Hoffman RJ, et al. **Effect of intravenous ondansetron on QTc interval in children with gastroenteritis.** *Am J Emerg Med* 2018. — [PMID 29029798](https://pubmed.ncbi.nlm.nih.gov/29029798/)
134. Yürük Mısırlıoğlu H, et al. **The effect of intravenous ondansetron on QT interval in the emergency department.** *Am J Emerg Med* 2024. — [PMID 39153265](https://pubmed.ncbi.nlm.nih.gov/39153265/)
135. Dias CS, et al. **Uncommon Extrapyramidal Reaction to Ondansetron.** *Cureus* 2025. — [PMID 40206890](https://pubmed.ncbi.nlm.nih.gov/40206890/)
136. Navari RM. **HTX-019: polysorbate 80- and synthetic surfactant-free neurokinin 1 receptor antagonist.** *Future Oncol* 2019. — [PMID 30304952](https://pubmed.ncbi.nlm.nih.gov/30304952/)
137. Kassar O, et al. **Efficacy and safety of ondansetron in schizophrenia.** *Gen Hosp Psychiatry* 2025. — [PMID 40482395](https://pubmed.ncbi.nlm.nih.gov/40482395/)

## 19. 전신 결과 — 신독성, 삶의 질, 용량강도

138. Crona DJ, et al. **A Systematic Review of Strategies to Prevent Cisplatin-Induced Nephrotoxicity.** *Oncologist* 2017. — [PMID 28438887](https://pubmed.ncbi.nlm.nih.gov/28438887/)
139. Mathavan A, et al. **Beneficial Effects of SGLT2 Inhibitors in Cisplatin Nephrotoxicity.** *Kidney360* 2026. — [PMID 41854743](https://pubmed.ncbi.nlm.nih.gov/41854743/)
140. Bloechl-Daum B, et al. **Delayed nausea and vomiting continue to reduce patients' quality of life after highly and moderately emetogenic chemotherapy despite antiemetic treatment.** *J Clin Oncol* 2006. — [PMID 16983116](https://pubmed.ncbi.nlm.nih.gov/16983116/)
141. Cohen L, et al. **Chemotherapy-induced nausea and vomiting: incidence and impact on patient quality of life.** *Support Care Cancer* 2007. — [PMID 17103197](https://pubmed.ncbi.nlm.nih.gov/17103197/)
142. Martin AR, et al. **Assessing the impact of chemotherapy-induced nausea and vomiting on patients' daily lives: a modified version of the Functional Living Index-Emesis (FLIE).** *Support Care Cancer* 2003. — [PMID 12827483](https://pubmed.ncbi.nlm.nih.gov/12827483/)
143. Han L, et al. **Development and validation of the Pediatrics Functional Living Index-Emesis scale.** *Front Oncol* 2025. — [PMID 40630203](https://pubmed.ncbi.nlm.nih.gov/40630203/)
144. Ihbe-Heffinger A, et al. **The impact of delayed chemotherapy-induced nausea and vomiting on patients, health resource utilization and costs.** *Ann Oncol* 2004. — [PMID 14998860](https://pubmed.ncbi.nlm.nih.gov/14998860/)
145. Haiderali A, et al. **Impact on daily functioning and indirect/direct costs associated with chemotherapy-induced nausea and vomiting.** *Support Care Cancer* 2011. — [PMID 20532923](https://pubmed.ncbi.nlm.nih.gov/20532923/)
146. Farrell C, et al. **The impact of chemotherapy-related nausea on patients' nutritional status, psychological distress and quality of life.** *Support Care Cancer* 2013. — [PMID 22610269](https://pubmed.ncbi.nlm.nih.gov/22610269/)
147. Aapro M, et al. **Oncologist perspectives on chemotherapy-induced nausea and vomiting (CINV) management and outcomes.** *Cancer Rep* 2018. — [PMID 32729252](https://pubmed.ncbi.nlm.nih.gov/32729252/)
148. Yeo W, et al. **Dataset on chemotherapy-induced nausea and vomiting (CINV) and quality of life (QOL) during multiple cycles of chemotherapy.** *Data Brief* 2020. — [PMID 32215313](https://pubmed.ncbi.nlm.nih.gov/32215313/)
149. Zhang M, et al. **Pharmacoeconomic Evaluation of Netupitant Palonosetron for Chemotherapy-Induced Nausea and Vomiting.** *Ther Innov Regul Sci* 2026. — [PMID 42207482](https://pubmed.ncbi.nlm.nih.gov/42207482/)
150. Homma S, et al. **Mallory-Weiss Syndrome Without Vomiting.** *Cureus* 2026. — [PMID 42453852](https://pubmed.ncbi.nlm.nih.gov/42453852/)

## 20. 보조요법

151. Ryan JL, et al. **Ginger (Zingiber officinale) reduces acute chemotherapy-induced nausea: a URCC CCOP study of 576 patients.** *Support Care Cancer* 2012. — [PMID 21818642](https://pubmed.ncbi.nlm.nih.gov/21818642/)
152. Bhargava R, et al. **The effect of ginger (Zingiber officinale Roscoe) in patients with advanced cancer.** *Support Care Cancer* 2020. — [PMID 31745695](https://pubmed.ncbi.nlm.nih.gov/31745695/)
153. Wenxi YU, et al. **Neiguan (PC6) acupoint stimulation for preventing chemotherapy-induced nausea and vomiting.** *J Tradit Chin Med* 2024. — [PMID 38767643](https://pubmed.ncbi.nlm.nih.gov/38767643/)
154. Zhu L, et al. **Progress in traditional Chinese medicine for chemotherapy-induced nausea and vomiting.** *Am J Transl Res* 2026. — [PMID 42491038](https://pubmed.ncbi.nlm.nih.gov/42491038/)
155. Chou R, et al. **Cannabis-Based Products for Chronic Pain: An Updated Systematic Review.** *Ann Intern Med* 2026. — [PMID 41429020](https://pubmed.ncbi.nlm.nih.gov/41429020/)
156. A Mahrous M, et al. **Evaluation of clinical outcomes and efficacy of palonosetron and granisetron in combination regimens.** *Cancer Chemother Pharmacol* 2021. — [PMID 33835230](https://pubmed.ncbi.nlm.nih.gov/33835230/)

## 21. 방법론 — QSP, mrgsolve, frailty 모형

157. **Elmokadem A, et al. Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. — [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
     *이 저장소의 ODE 구현 도구.*
158. **Takawira HT, et al. Beyond the Cox model: a log-logistic accelerated failure time frailty model for analysing first-event times.** *BMC Med Res Methodol* 2026. — [PMID 42332583](https://pubmed.ncbi.nlm.nih.gov/42332583/)
     *CR = E_S[exp(−S·I)] 의 lognormal frailty 구조에 대한 방법론적 근거.*
159. Wu H, et al. **Proportional likelihood ratio mixed model for discrete longitudinal data.** *Stat Med* 2021. — [PMID 33588517](https://pubmed.ncbi.nlm.nih.gov/33588517/)
160. Bender G, et al. **Pharmacokinetic-pharmacodynamic analysis of a drug response using a mechanistic model.** *J Pharmacol Exp Ther* 2010. — [PMID 20444880](https://pubmed.ncbi.nlm.nih.gov/20444880/)
161. Dogné JM, et al. **From animal testing to model-informed drug development: building on ICH M15 and EMA initiatives.** *Front Pharmacol* 2026. — [PMID 42499501](https://pubmed.ncbi.nlm.nih.gov/42499501/)
162. Choi Y, et al. **Integrating Real-World Data and Pharmacometrics to Bridge Evidence Gaps in Special Populations.** *Pharmaceutics* 2026. — [PMID 42514883](https://pubmed.ncbi.nlm.nih.gov/42514883/)
163. Roila F. **[Research on antiemetics: an Italian model of success].** *Tumori* 1998. — [PMID 9617377](https://pubmed.ncbi.nlm.nih.gov/9617377/)

---

## 부록 A — 앵커 숫자와 그 출처의 일대일 대응

| # | 앵커 | 값 | 출처 | PMID |
|---|------|-----|------|------|
| 1 | CR 급성 (0–24 h), ondansetron + dexamethasone | 78.1% | Hesketh 2003 대조군 | 14559886 |
| 2 | CR 지연 (24–120 h), 동일 군 | 55.8% | Hesketh 2003 대조군 | 14559886 |
| 3 | CR 전체 (0–120 h), 동일 군 | 52.3% | Hesketh 2003 대조군 | 14559886 |
| 4 | CR 급성, + aprepitant | 89.2% | Hesketh 2003 시험군 | 14559886 |
| 5 | CR 지연, + aprepitant | 75.4% | Hesketh 2003 시험군 | 14559886 |
| 6 | CR 전체, + aprepitant | 72.7% | Hesketh 2003 시험군 | 14559886 |
| 7 | 미치료 시스플라틴 24시간 구토 횟수 | ≈6회 | Cubeddu 1996 | 8692546 |
| 8 | 미치료 구역 없음 (0–120 h) | ≈2% | 자연경과 | 8692546 |
| 9 | 3제 병용 구역 없음 (0–120 h) | 22% | Navari 2016 **대조군** | 27705262 |

**적합에 사용하지 않은 것 (전부 held-out):** 팔로노세트론·NEPA·롤라피탄트·올란자핀의
모든 효능 수치, AC/카보플라틴/옥살리플라틴 요법, 덱사메타손 AUC 비, NK1 PET 점유율,
5-HIAA 시간경과, QTc 기울기, CYP2D6 초고속대사자, 위험인자 공변량, 다주기 예기성 이월.

---

## 부록 B — 검증 방법

```bash
# 이 목록의 모든 PMID는 아래와 같이 실제 조회하여 서지사항을 확인했습니다.
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=<PMID>"
```

레코드가 조회되지 않은 인용은 이 파일에서 제외했습니다. 따라서 위 163개 항목의
저자·연도·저널·제목은 모두 PubMed에서 확인된 실제 문헌입니다. 다만 **본문의 세부
수치를 그 문헌에 귀속시킬 때는** 제목 수준의 확인을 넘어서는 검토가 필요하며,
이 모델이 앵커로 사용한 9개 숫자(부록 A)만이 원문에서 직접 대조된 값입니다.
그 밖의 인용은 해당 기전·파라미터 영역의 **근거 문헌 목록**으로 이해하시기 바랍니다.
