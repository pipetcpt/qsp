# 선천성 고인슐린증 (Congenital Hyperinsulinism, CHI) — 참고문헌

각 절의 마지막 줄은 **이 모델의 어느 항이 그 문헌에 의존하는지**를 밝힌다.
보정에 실제로 쓴 숫자(N1–N8, C1)와 예측(P1–P10)의 정의는
[`README.md`](README.md) 및 [`chi_mrgsolve_model.R`](chi_mrgsolve_model.R) 머리말에 있다.

---

## 1. 총설 · 진단 · 관리 지침 (Reviews, diagnosis, guidelines)

1. Stanley CA. **Perspective on the genetics and diagnosis of congenital hyperinsulinism disorders.** J Clin Endocrinol Metab. 2016;101(3):815–26. <https://pubmed.ncbi.nlm.nih.gov/26908103/>
2. De Leon DD, Stanley CA. **Congenital hypoglycemia disorders: new aspects of etiology, diagnosis, treatment and outcomes.** Pediatr Diabetes. 2017;18(1):3–9. <https://pubmed.ncbi.nlm.nih.gov/28074627/>
3. Lord K, De León DD. **Monogenic hyperinsulinemic hypoglycemia: current insights into the pathogenesis and management.** Int J Pediatr Endocrinol. 2013;2013(1):3. <https://pubmed.ncbi.nlm.nih.gov/23448374/>
4. Thornton PS, Stanley CA, De Leon DD, et al. **Recommendations from the Pediatric Endocrine Society for evaluation and management of persistent hypoglycemia in neonates, infants, and children.** J Pediatr. 2015;167(2):238–45. <https://pubmed.ncbi.nlm.nih.gov/25957977/>
5. Banerjee I, Salomon-Estebanez M, Shah P, et al. **Therapies and outcomes of congenital hyperinsulinism-induced hypoglycaemia.** Diabet Med. 2019;36(1):9–21. <https://pubmed.ncbi.nlm.nih.gov/30246418/>
6. Demirbilek H, Hussain K. **Congenital hyperinsulinism: diagnosis and treatment update.** J Clin Res Pediatr Endocrinol. 2017;9(Suppl 2):69–87. <https://pubmed.ncbi.nlm.nih.gov/29280739/>
7. Rosenfeld E, Ganguly A, De Leon DD. **Congenital hyperinsulinism disorders: genetic and clinical characteristics.** Am J Med Genet C. 2019;181(4):682–92. <https://pubmed.ncbi.nlm.nih.gov/31414570/>
8. Galcheva S, Demirbilek H, Al-Khawaga S, Hussain K. **The genetic and molecular mechanisms of congenital hyperinsulinism.** Front Endocrinol. 2019;10:111. <https://pubmed.ncbi.nlm.nih.gov/30873120/>

> 모델 의존: GIR > 8 mg/kg/min 이라는 CHI의 작업정의(진단 노드), critical sample 항목,
> 목표 혈당 70 mg/dL 라는 관리 기준 — 모델은 이 목표를 **가정하지 않고**
> 뇌 연료 산술로 재도출한다(P5).

---

## 2. K_ATP 채널 — 이 모델의 축(g)이 되는 분자 (ABCC8 / KCNJ11)

9. Thomas PM, Cote GJ, Wohllk N, et al. **Mutations in the sulfonylurea receptor gene in familial persistent hyperinsulinemic hypoglycemia of infancy.** Science. 1995;268(5209):426–9. <https://pubmed.ncbi.nlm.nih.gov/7716548/>
10. Nichols CG, Shyng SL, Nestorowicz A, et al. **Adenosine diphosphate as an intracellular regulator of insulin secretion.** Science. 1996;272(5269):1785–7. <https://pubmed.ncbi.nlm.nih.gov/8650576/>
11. Ashcroft FM. **ATP-sensitive potassium channelopathies: focus on insulin secretion.** J Clin Invest. 2005;115(8):2047–58. <https://pubmed.ncbi.nlm.nih.gov/16075048/>
12. Ashcroft FM, Rorsman P. **K_ATP channels and islet hormone secretion: new insights and controversies.** Nat Rev Endocrinol. 2013;9(11):660–9. <https://pubmed.ncbi.nlm.nih.gov/24042324/>
13. Shyng SL, Nichols CG. **Octameric stoichiometry of the K_ATP channel complex.** J Gen Physiol. 1997;110(6):655–64. <https://pubmed.ncbi.nlm.nih.gov/9382894/>
14. Cartier EA, Conti LR, Vandenberg CA, Shyng SL. **Defective trafficking and function of K_ATP channels caused by a sulfonylurea receptor 1 mutation associated with persistent hyperinsulinemic hypoglycemia of infancy.** Proc Natl Acad Sci USA. 2001;98(5):2882–7. <https://pubmed.ncbi.nlm.nih.gov/11226335/>
15. Snider KE, Becker S, Boyajian L, et al. **Genotype and phenotype correlations in 417 children with congenital hyperinsulinism.** J Clin Endocrinol Metab. 2013;98(2):E355–63. <https://pubmed.ncbi.nlm.nih.gov/23275527/>
16. Pinney SE, MacMullen C, Becker S, et al. **Clinical characteristics and biochemical mechanisms of congenital hyperinsulinism associated with dominant K_ATP channel mutations.** J Clin Invest. 2008;118(8):2877–86. <https://pubmed.ncbi.nlm.nih.gov/18596924/>
17. Macmullen CM, Zhou Q, Snider KE, et al. **Diazoxide-unresponsive congenital hyperinsulinism in children with dominant mutations of the beta-cell sulfonylurea receptor SUR1.** Diabetes. 2011;60(6):1797–804. <https://pubmed.ncbi.nlm.nih.gov/21536951/>
18. Saint-Martin C, Arnoux JB, de Lonlay P, Bellanné-Chantelot C. **KATP channel mutations in congenital hyperinsulinism.** Semin Pediatr Surg. 2011;20(1):18–22. <https://pubmed.ncbi.nlm.nih.gov/?term=KATP+channel+mutations+congenital+hyperinsulinism+Saint-Martin>

> 모델 의존: **핵심.** 표면 채널 밀도 g가 단 하나의 병변 변수라는 구조,
> 그리고 열성 대립유전자 다수가 *트래피킹* 실패로 g→0을 만든다는 점(14번)이
> "디아족사이드는 열 채널이 없으면 아무것도 못 한다"(P2)의 분자적 근거다.
> 15·16·17번은 유전형–표현형 분리를 제공하며, 모델은 그 분리를 **재현할 뿐**
> 그것으로 보정하지 않는다.

---

## 3. 막전위 · Ca²⁺ · 분비 커플링 (전압 분배기의 실증)

19. Ashcroft FM, Harrison DE, Ashcroft SJ. **Glucose induces closure of single potassium channels in isolated rat pancreatic beta-cells.** Nature. 1984;312(5993):446–8. <https://pubmed.ncbi.nlm.nih.gov/6095103/>
20. Rorsman P, Ashcroft FM. **Pancreatic β-cell electrical activity and insulin secretion: of mice and men.** Physiol Rev. 2018;98(1):117–214. <https://pubmed.ncbi.nlm.nih.gov/29212789/>
21. Henquin JC. **Triggering and amplifying pathways of regulation of insulin secretion by glucose.** Diabetes. 2000;49(11):1751–60. <https://pubmed.ncbi.nlm.nih.gov/11078440/>
22. Henquin JC. **Regulation of insulin secretion: a matter of phase control and amplitude modulation.** Diabetologia. 2009;52(5):739–51. <https://pubmed.ncbi.nlm.nih.gov/19288070/>
23. Straub SG, Sharp GW. **Hypothesis: one rate-limiting step controls the magnitude of both phases of glucose-stimulated insulin secretion.** Am J Physiol Cell Physiol. 2004;287(3):C565–71. <https://pubmed.ncbi.nlm.nih.gov/15308464/>
24. Matschinsky FM. **Banting Lecture 1995: A lesson in metabolic regulation inspired by the glucokinase glucose sensor paradigm.** Diabetes. 1996;45(2):223–41. <https://pubmed.ncbi.nlm.nih.gov/8549869/>

> 모델 의존: `Popen → G_KATP → V_m → fCa → 분비`의 사슬,
> 그리고 **증폭 경로**(21·22번)가 저혈당에서 35 %의 바닥을 갖는다는 항.
> 이 바닥이 g=0에서 분비가 무한이 아니라 유한한 이유다.

---

## 4. 비-K_ATP 유전형 — 같은 채널을 다른 쪽에서 건드리는 방법

25. Stanley CA, Lieu YK, Hsu BY, et al. **Hyperinsulinism and hyperammonemia in infants with regulatory mutations of the glutamate dehydrogenase gene.** N Engl J Med. 1998;338(19):1352–7. <https://pubmed.ncbi.nlm.nih.gov/9571255/>
26. Stanley CA, Fang J, Kutyna K, et al. **Molecular basis and characterization of the hyperinsulinism/hyperammonemia syndrome: predominance of mutations in exons 11 and 12 of the glutamate dehydrogenase gene.** Diabetes. 2000;49(4):667–73. <https://pubmed.ncbi.nlm.nih.gov/10871207/>
27. Kelly A, Ng D, Ferry RJ Jr, et al. **Acute insulin responses to leucine in children with the hyperinsulinism/hyperammonemia syndrome.** J Clin Endocrinol Metab. 2001;86(8):3724–8. <https://pubmed.ncbi.nlm.nih.gov/11502802/>
28. Glaser B, Kesavan P, Heyman M, et al. **Familial hyperinsulinism caused by an activating glucokinase mutation.** N Engl J Med. 1998;338(4):226–30. <https://pubmed.ncbi.nlm.nih.gov/9435328/>
29. Christesen HB, Jacobsen BB, Odili S, et al. **The second activating glucokinase mutation (A456V): implications for glucose homeostasis and diabetes therapy.** Diabetes. 2002;51(4):1240–6. <https://pubmed.ncbi.nlm.nih.gov/11916951/>
30. Clayton PT, Eaton S, Aynsley-Green A, et al. **Hyperinsulinism in short-chain L-3-hydroxyacyl-CoA dehydrogenase deficiency reveals the importance of beta-oxidation in insulin secretion.** J Clin Invest. 2001;108(3):457–65. <https://pubmed.ncbi.nlm.nih.gov/11489939/>
31. Pearson ER, Boj SF, Steele AM, et al. **Macrosomia and hyperinsulinaemic hypoglycaemia in patients with heterozygous mutations in the HNF4A gene.** PLoS Med. 2007;4(4):e118. <https://pubmed.ncbi.nlm.nih.gov/17407387/>
32. Otonkoski T, Jiao H, Kaminen-Ahola N, et al. **Physical exercise-induced hypoglycemia caused by failed silencing of monocarboxylate transporter 1 in pancreatic beta cells.** Am J Hum Genet. 2007;81(3):467–74. <https://pubmed.ncbi.nlm.nih.gov/17701894/>
33. González-Barroso MM, Giurgea I, Bouillaud F, et al. **Mutations in UCP2 in congenital hyperinsulinism reveal a role for regulation of insulin secretion.** PLoS One. 2008;3(12):e3850. <https://pubmed.ncbi.nlm.nih.gov/19065272/>
34. Flanagan SE, Vairo F, Johnson MB, et al. **A CACNA1D mutation in a patient with persistent hyperinsulinaemic hypoglycaemia, heart defects, and severe hypotonia.** Pediatr Diabetes. 2017;18(4):320–3. <https://pubmed.ncbi.nlm.nih.gov/28107785/>

> 모델 의존: 이 유전형들은 **모두 g=1을 유지한 채** Rt(대사신호)나 GCK의 S₀.₅,
> 또는 GDH의 GTP 브레이크만 바꾼다. 그래서 모델은 이들이 디아족사이드에
> **잘 반응한다**는 것을 예측으로 내놓는다(P2, 시나리오 4).
> 27번은 류신 부하 예측(P8, 37 vs 16 mg/dL)의 비교 기준이다.

---

## 5. 국소형 대 미만성 · 병리 · ¹⁸F-DOPA PET

35. de Lonlay P, Fournet JC, Rahier J, et al. **Somatic deletion of the imprinted 11p15 region in sporadic persistent hyperinsulinemic hypoglycemia of infancy is specific of focal adenomatous hyperplasia.** J Clin Invest. 1997;100(4):802–7. <https://pubmed.ncbi.nlm.nih.gov/9259578/>
36. Rahier J, Guiot Y, Sempoux C. **Persistent hyperinsulinaemic hypoglycaemia of infancy: a heterogeneous syndrome unrelated to nesidioblastosis.** Arch Dis Child Fetal Neonatal Ed. 2000;82(2):F108–12. <https://pubmed.ncbi.nlm.nih.gov/10685983/>
37. Ribeiro MJ, De Lonlay P, Delzescaux T, et al. **Characterization of hyperinsulinism in infancy assessed with PET and ¹⁸F-fluoro-L-DOPA.** J Nucl Med. 2005;46(4):560–6. <https://pubmed.ncbi.nlm.nih.gov/15809475/>
38. Hardy OT, Hernandez-Pampaloni M, Saffer JR, et al. **Accuracy of ¹⁸F-fluorodopa positron emission tomography for diagnosing and localizing focal congenital hyperinsulinism.** J Clin Endocrinol Metab. 2007;92(12):4706–11. <https://pubmed.ncbi.nlm.nih.gov/17895314/>
39. Laje P, States LJ, Zhuang H, et al. **Accuracy of PET/CT scan in the diagnosis of the focal form of congenital hyperinsulinism.** J Pediatr Surg. 2013;48(2):388–93. <https://pubmed.ncbi.nlm.nih.gov/23414871/>
40. Sempoux C, Capito C, Bellanné-Chantelot C, et al. **Morphological mosaicism of the pancreatic islets: a novel anatomopathological form of persistent hyperinsulinemic hypoglycemia of infancy.** J Clin Endocrinol Metab. 2011;96(12):3785–93. <https://pubmed.ncbi.nlm.nih.gov/21956412/>

> 모델 의존: 두 개의 베타세포 집단(`w_ab`, `dens_a`)이라는 구조.
> 36·40번의 "확대된 핵을 가진 고밀도 베타세포" 기술이 P10의 정량적 요구
> (병소 분율 × 국소 밀도 ≈ 1.0)와 맞물린다. 모델은 이를 **예측이자 요구**로
> 제시한다: 작은 병소가 미만성과 같은 중증도를 내려면 밀도가 ~10배여야 한다.

---

## 6. 디아족사이드 — g에 곱해지는 약

41. Drash A, Wolff F. **Drug therapy in leucine-sensitive hypoglycemia.** Metabolism. 1964;13:487–92. <https://pubmed.ncbi.nlm.nih.gov/14169215/>
42. Trube G, Rorsman P, Ohno-Shosaku T. **Opposite effects of tolbutamide and diazoxide on the ATP-dependent K⁺ channel in mouse pancreatic beta-cells.** Pflugers Arch. 1986;407(5):493–9. <https://pubmed.ncbi.nlm.nih.gov/3786110/>
43. Shyng SL, Ferrigni T, Nichols CG. **Regulation of K_ATP channel activity by diazoxide and MgADP: distinct functions of the two nucleotide binding folds of the sulfonylurea receptor.** J Gen Physiol. 1997;110(6):643–54. <https://pubmed.ncbi.nlm.nih.gov/9382893/>
44. Herman GA, Hussain K. **Diazoxide in the management of congenital hyperinsulinism: dose, response and adverse effects.** (임상 시리즈 리뷰) Arch Dis Child. 2015;100(9):884–5. <https://pubmed.ncbi.nlm.nih.gov/25977564/>
45. Welters A, Lerch C, Kummer S, et al. **Long-term medical treatment in congenital hyperinsulinism: a descriptive analysis in a large cohort of patients from different clinical centers.** Orphanet J Rare Dis. 2015;10:150. <https://pubmed.ncbi.nlm.nih.gov/26608306/>
46. Thornton P, Truong L, Reynolds C, et al. **Rate of serious adverse events associated with diazoxide treatment of patients with hyperinsulinism.** Horm Res Paediatr. 2019;91(1):25–32. <https://pubmed.ncbi.nlm.nih.gov/30836359/>
47. Timlin MR, Black AB, Delaney HM, et al. **Development of pulmonary hypertension during treatment with diazoxide: a case series and literature review.** Pediatr Cardiol. 2017;38(6):1247–50. <https://pubmed.ncbi.nlm.nih.gov/28642988/>
48. Pruitt LG, Bloomstone J, et al. **Diazoxide pharmacokinetics in the newborn** — 신생아 t½ 9–30 h 및 단백결합 >90 %에 관한 자료. Clin Pharmacol Ther. (고전 자료) <https://pubmed.ncbi.nlm.nih.gov/?term=diazoxide+pharmacokinetics+newborn>

> 모델 의존: **43번이 결정적이다.** 디아족사이드가 SUR1의 NBD2에서 MgADP-의존적으로
> 작용한다는 사실 때문에, 모델은 이 약을 "일정 비율의 채널을 강제로 여는 것"이 아니라
> **채널의 ATP/ADP 설정점(R50)을 오른쪽으로 미는 것**으로 쓴다. 처음에 전자로 썼을 때는
> 2.5 mg/kg/day가 정상 베타세포의 분비를 없애버렸고(혈당 260 mg/dL), 이는 검증 과정에서
> 발견해 고쳤다. PK는 48번, 독성 노드는 46·47번.

---

## 7. 소마토스타틴 유사체 — 분배기에 더해지는 약

49. Yamada Y, Post SR, Wang K, et al. **Cloning and functional characterization of a family of human and mouse somatostatin receptors expressed in brain, gastrointestinal tract, and kidney.** Proc Natl Acad Sci USA. 1992;89(1):251–5. <https://pubmed.ncbi.nlm.nih.gov/1346068/>
50. Rorsman P, Braun M. **Regulation of insulin secretion in human pancreatic islets.** Annu Rev Physiol. 2013;75:155–79. <https://pubmed.ncbi.nlm.nih.gov/22974438/>
51. Kailey B, van de Bunt M, Cheley S, et al. **SSTR2 is the functionally dominant somatostatin receptor in human pancreatic β- and α-cells.** Am J Physiol Endocrinol Metab. 2012;303(9):E1107–16. <https://pubmed.ncbi.nlm.nih.gov/22932785/>
52. Le Quan Sang KH, Arnoux JB, Mamoune A, et al. **Successful treatment of congenital hyperinsulinism with long-acting release octreotide.** Eur J Endocrinol. 2012;166(2):333–9. <https://pubmed.ncbi.nlm.nih.gov/22084154/>
53. Demirbilek H, Shah P, Arya VB, et al. **Long-term follow-up of children with congenital hyperinsulinism on octreotide therapy.** J Clin Endocrinol Metab. 2014;99(10):3660–7. <https://pubmed.ncbi.nlm.nih.gov/25004250/>
54. Laje P, Halaby L, Adzick NS, Stanley CA. **Necrotizing enterocolitis in neonates receiving octreotide for the management of congenital hyperinsulinism.** Pediatr Diabetes. 2010;11(2):142–7. <https://pubmed.ncbi.nlm.nih.gov/19558634/>
55. Corda H, Kummer S, Welters A, et al. **Treatment with long-acting lanreotide autogel in early infancy in patients with severe neonatal hyperinsulinism.** Orphanet J Rare Dis. 2017;12(1):108. <https://pubmed.ncbi.nlm.nih.gov/28615019/>

> 모델 의존: **51번이 이 모델의 두 번째 축이다.** SSTR2가 Gi를 통해 cAMP를 낮추는
> 동시에 **G-단백 조절 K⁺ 전도도(GIRK)를 연다**는 사실 때문에, 옥트레오타이드는
> 돌연변이를 지니지 않은 **두 번째 채널**을 가져온다. 이것이 g=0에서도 효과가 남는
> 이유(P3: 59 %)이고, 디아족사이드와 정반대 결론이 나오는 유일한 이유다.
> `gGIRK = 1.2`는 이 모델이 CHI에서 가져온 **단 하나의 숫자**(C1)이며 53번에 맞췄다.
> 54번은 신생아 NEC 위험 노드.

---

## 8. 글루카곤 · 인슐린수용체 표적 · 그 외 약제

56. Neylon OM, Moran MM, Pellicano A, et al. **Successful subcutaneous glucagon use for persistent hypoglycaemia in congenital hyperinsulinism.** J Pediatr Endocrinol Metab. 2013;26(11-12):1157–61. <https://pubmed.ncbi.nlm.nih.gov/23751387/>
57. Hawkes CP, Lado JJ, Givler S, De Leon DD. **The effect of continuous intravenous glucagon on glucose requirements in infants with congenital hyperinsulinism.** JIMD Rep. 2019;45:45–50. <https://pubmed.ncbi.nlm.nih.gov/30569021/>
58. Ferrara CT, Boodhansingh KE, Paradies E, et al. **Novel hypoglycemia phenotype in congenital hyperinsulinism due to dominant mutations of uncoupling protein 2.** J Clin Endocrinol Metab. 2017;102(3):942–9. <https://pubmed.ncbi.nlm.nih.gov/28324055/>
59. Calabria AC, Li C, Gallagher PR, et al. **GLP-1 receptor antagonist exendin-(9-39) elevates fasting blood glucose levels in congenital hyperinsulinism owing to inactivating mutations in the ATP-sensitive K⁺ channel.** Diabetes. 2012;61(10):2585–91. <https://pubmed.ncbi.nlm.nih.gov/22855730/>
60. Ng CM, Tang F, Seeholzer SH, et al. **Population pharmacokinetics of exendin-(9-39) and clinical dose selection in patients with congenital hyperinsulinism.** Br J Clin Pharmacol. 2018;84(3):520–32. <https://pubmed.ncbi.nlm.nih.gov/29148594/>
61. Senniappan S, Alexandrescu S, Tatevian N, et al. **Sirolimus therapy in infants with severe hyperinsulinemic hypoglycemia.** N Engl J Med. 2014;370(12):1131–7. <https://pubmed.ncbi.nlm.nih.gov/24645944/>
62. Szymanowski M, Estebanez MS, Padidela R, et al. **mTOR inhibitors for the treatment of severe congenital hyperinsulinism: perspectives on limited therapeutic success.** J Clin Endocrinol Metab. 2016;101(12):4719–29. <https://pubmed.ncbi.nlm.nih.gov/27715326/>
63. Corbin JA, Bhaskar V, Goldfine ID, et al. **Improved glucose metabolism in vitro and in vivo by an allosteric monoclonal antibody that increases insulin receptor binding affinity.** PLoS One. 2014;9(2):e88684. <https://pubmed.ncbi.nlm.nih.gov/24533135/>
64. Johnson MB, De Franco E, et al. / **RZ358 (ersodetug) — insulin receptor antagonist antibody in congenital hyperinsulinism** (개발 프로그램 및 초기 임상). 검색: <https://pubmed.ncbi.nlm.nih.gov/?term=RZ358+OR+ersodetug+hyperinsulinism>
65. Müller D, Zimmering M, Roehr CC. **Should nifedipine be used to counter low blood sugar levels in children with persistent hyperinsulinaemic hypoglycaemia?** Arch Dis Child. 2004;89(1):83–5. <https://pubmed.ncbi.nlm.nih.gov/14709521/>

> 모델 의존: 57번은 글루카곤 지속주입의 포도당 요구량 감소(모델 35 %)의 비교 기준.
> 59·60번은 exendin(9-39) 노드와 그 PK. 61 대 62번의 **대비**가 시롤리무스를
> "기전 불명, 효과 modest"로 쓰게 한 근거이며 모델은 19 %를 낸다(T4에서 그 불확실성을 명시).
> 63·64번은 에르소데투그가 **수용체 수준**, 즉 g와 무관하게 작동함을 근거지운다(T3에 정량적 한계 기재).
> 65번은 니페디핀 실패의 임상적 확인이며, 모델은 그 실패를
> "기전은 옳지만 베타세포 EC₅₀가 혈장 Cmax의 약 8배"라는 **약동학적 이유**로 설명한다(P9).

---

## 9. 수술 — B를 줄이고 g는 남기는 개입

66. Adzick NS, De Leon DD, States LJ, et al. **Surgical treatment of congenital hyperinsulinism: results from 500 pancreatectomies in neonates and children.** J Pediatr Surg. 2019;54(1):27–32. <https://pubmed.ncbi.nlm.nih.gov/30343978/>
67. Lord K, Dzata E, Snider KE, et al. **Clinical presentation and management of children with diffuse and focal hyperinsulinism: a review of 223 cases.** J Clin Endocrinol Metab. 2013;98(11):E1786–9. <https://pubmed.ncbi.nlm.nih.gov/24014812/>
68. Beltrand J, Caquard M, Arnoux JB, et al. **Glucose metabolism in 105 children and adolescents after pancreatectomy for congenital hyperinsulinism.** Diabetes Care. 2012;35(2):198–203. <https://pubmed.ncbi.nlm.nih.gov/22190679/>
69. Arya VB, Senniappan S, Demirbilek H, et al. **Pancreatic endocrine and exocrine function in children following near-total pancreatectomy for diffuse congenital hyperinsulinism.** PLoS One. 2014;9(5):e98054. <https://pubmed.ncbi.nlm.nih.gov/24840042/>
70. Lord K, Radcliffe J, Gallagher PR, et al. **High risk of diabetes and neurobehavioral deficits in individuals with surgically treated hyperinsulinism.** J Clin Endocrinol Metab. 2015;100(11):4133–9. <https://pubmed.ncbi.nlm.nih.gov/26327482/>

> 모델 의존: **68·70번이 P7의 비교 기준이다.** 모델에서 수술은 `BMASS0`만 바꾸고
> `g`는 건드리지 않으므로, 하나의 방정식이 (i) 잔존량 0.5에서의 지속 저혈당,
> (ii) 0.2–0.3에서의 정상 혈당, (iii) ≤0.1에서의 즉각적 당뇨를 모두 만든다.
> 성장이 같은 잔존량을 다시 읽는다는 논증(3.5 kg의 2 %가 20 kg에서 0.35 %)이
> 68·70번의 "청소년기까지 약 절반이 당뇨"를 **새로운 가정 없이** 설명한다.

---

## 10. 뇌 연료 — 인슐린이 삭제하는 예비 연료

71. Cremer JE. **Substrate utilization and brain development.** J Cereb Blood Flow Metab. 1982;2(4):394–407. <https://pubmed.ncbi.nlm.nih.gov/6754748/>
72. Nehlig A. **Brain uptake and metabolism of ketone bodies in animal models.** Prostaglandins Leukot Essent Fatty Acids. 2004;70(3):265–75. <https://pubmed.ncbi.nlm.nih.gov/14769485/>
73. Vannucci SJ, Simpson IA. **Developmental switch in brain nutrient transporter expression in the rat.** Am J Physiol Endocrinol Metab. 2003;285(5):E1127–34. <https://pubmed.ncbi.nlm.nih.gov/14534077/>
74. Pellerin L, Pellegri G, Bittar PG, et al. **Evidence supporting the existence of an activity-dependent astrocyte-neuron lactate shuttle.** Dev Neurosci. 1998;20(4-5):291–9. <https://pubmed.ncbi.nlm.nih.gov/9778565/>
75. Bougneres PF, Lemmel C, Ferré P, Bier DM. **Ketone body transport in the human neonate and infant.** J Clin Invest. 1986;77(1):42–8. <https://pubmed.ncbi.nlm.nih.gov/2867094/>
76. Kinnala A, Rikalainen H, Lapinleimu H, et al. **Cerebral magnetic resonance imaging and ultrasonography findings after neonatal hypoglycemia.** Pediatrics. 1999;103(4 Pt 1):724–9. <https://pubmed.ncbi.nlm.nih.gov/10103293/>
77. Burns CM, Rutherford MA, Boardman JP, Cowan FM. **Patterns of cerebral injury and neurodevelopmental outcomes after symptomatic neonatal hypoglycemia.** Pediatrics. 2008;122(1):65–74. <https://pubmed.ncbi.nlm.nih.gov/18595988/>
78. Menni F, de Lonlay P, Sevin C, et al. **Neurologic outcomes of 90 neonates and infants with persistent hyperinsulinemic hypoglycemia.** Pediatrics. 2001;107(3):476–9. <https://pubmed.ncbi.nlm.nih.gov/11230585/>
79. Avatapalle HB, Banerjee I, Shah S, et al. **Abnormal neurodevelopmental outcomes are common in children with transient congenital hyperinsulinism.** Front Endocrinol. 2013;4:60. <https://pubmed.ncbi.nlm.nih.gov/23730298/>
80. Ludwig A, Enke S, Heindorf J, et al. **Formal neurocognitive testing in 60 patients with congenital hyperinsulinism.** Horm Res Paediatr. 2018;89(1):1–6. <https://pubmed.ncbi.nlm.nih.gov/29131017/>

> 모델 의존: **P5의 전부.** 71·72·75번이 신생아 뇌가 케톤을 실제로 대량 산화한다는
> 근거이고, 브레인 MCT의 Km(≈6 mM, 72·74번)이 *속도제한 단계*라는 점이 결정적이다.
> Km을 6 mM로 두면 BOHB 2 mM에서 케톤은 요구량의 19 %만 대며, 그 결과 케톤성 저혈당의
> 허용 혈당이 25 mg/dL로 계산되어 임상 관찰(30–40 mg/dL대)에 근접한다.
> Km을 1.5 mM로 잘못 두면 허용 혈당이 16 mg/dL로 계산되어 임상과 어긋났고,
> 이 불일치가 Km 선택을 강제했다 — 보정이 아니라 **구조적 제약**이다.
> 76–80번은 결과(NDD) 노드.

---

## 11. 인슐린의 용량–반응 분리 (IC₅₀ 서열) 및 정상 신생아 대사

81. Rizza RA, Mandarino LJ, Gerich JE. **Dose-response characteristics for effects of insulin on production and utilization of glucose in man.** Am J Physiol. 1981;240(6):E630–9. <https://pubmed.ncbi.nlm.nih.gov/7018254/>
82. Nurjhan N, Campbell PJ, Kennedy FP, et al. **Insulin dose-response characteristics for suppression of glycerol release and conversion to glucose in humans.** Diabetes. 1986;35(12):1326–31. <https://pubmed.ncbi.nlm.nih.gov/3536590/>
83. Keller U, Schnell H, Sonnenberg GE, et al. **Role of insulin in the regulation of ketone body metabolism in man.** Diabetes Care. 1980;3(1):137–41. <https://pubmed.ncbi.nlm.nih.gov/6997936/>
84. Bier DM, Leake RD, Haymond MW, et al. **Measurement of "true" glucose production rates in infancy and childhood with 6,6-dideuteroglucose.** Diabetes. 1977;26(11):1016–23. <https://pubmed.ncbi.nlm.nih.gov/198425/>
85. Kalhan SC, Savin SM, Adam PA. **Measurement of glucose turnover in human newborn with glucose-1-¹³C.** J Clin Endocrinol Metab. 1976;43(3):704–7. <https://pubmed.ncbi.nlm.nih.gov/956353/>
86. Shelley HJ, Neligan GA. **Neonatal hypoglycaemia.** Br Med Bull. 1966;22(1):34–9. <https://pubmed.ncbi.nlm.nih.gov/5321967/>
87. Hume R, Burchell A, Williams FL, Koh DK. **Glucose homeostasis in the newborn.** Early Hum Dev. 2005;81(1):95–101. <https://pubmed.ncbi.nlm.nih.gov/15707720/>
88. Stanley CA, Rozance PJ, Thornton PS, et al. **Re-evaluating "transitional neonatal hypoglycemia": mechanism and implications for management.** J Pediatr. 2015;166(6):1520–5. <https://pubmed.ncbi.nlm.nih.gov/25819173/>

> 모델 의존: **N8(IC₅₀ 서열)과 N6·N7의 근거.** 81–83번이 지방분해·케톤생성이
> 간 포도당 생성보다 **훨씬 낮은** 인슐린 농도에서 억제된다는 것을 보여주며,
> 이 서열이 "CHI 저혈당은 언제나 무케톤성"이라는 결론을 자동으로 만든다.
> 84·85번은 신생아 포도당 회전율 4–6 mg/kg/min, 86·87번은 만삭 간 글리코겐 저장량,
> 88번은 정상 신생아의 이행기 저혈당이 **같은 K_ATP 설정점 이동**이라는 관점을 제공한다.

---

## 12. 진단 검사 — 각 검사가 모델의 어느 항을 읽는가

89. Finegold DN, Stanley CA, Baker L. **Glycemic response to glucagon during fasting hypoglycemia: an aid in the diagnosis of hyperinsulinism.** J Pediatr. 1980;96(2):257–9. <https://pubmed.ncbi.nlm.nih.gov/6766008/>
90. Ferrara C, Patel P, Becker S, et al. **Biomarkers of insulin for the diagnosis of hyperinsulinemic hypoglycemia in infants and children.** J Pediatr. 2016;168:212–9. <https://pubmed.ncbi.nlm.nih.gov/26490124/>
91. Palladino AA, Bennett MJ, Stanley CA. **Hyperinsulinism in infancy and childhood: when an insulin level is not always enough.** Clin Chem. 2008;54(2):256–63. <https://pubmed.ncbi.nlm.nih.gov/18156285/>
92. Hussain K, Hindmarsh P, Aynsley-Green A. **Neonates with symptomatic hyperinsulinemic hypoglycemia generate inappropriately low serum cortisol counterregulatory hormonal responses.** J Clin Endocrinol Metab. 2003;88(9):4342–7. <https://pubmed.ncbi.nlm.nih.gov/12970308/>

> 모델 의존: **89번이 P6의 직접적 대응물이다.** 글루카곤 자극검사가 진단적인 이유는
> 인슐린이 글리코겐을 *쓰지 못하게 지켜두기* 때문이며, 모델에서 글리코겐이
> **상태변수**라는 선택이 이 검사의 판별력을 자동으로 만든다(CHI +62 vs 정상 20 h 금식 +16 mg/dL).
> 91번은 인슐린:C-펩타이드 비 노드, 92번은 반대조절 실패 노드 — 모델에서는
> 별도의 병변이 아니라 α세포가 국소 인슐린에 억제된 **결과**로 나온다.

---

## 13. 시상하부 K_ATP · 저혈당 인식 (같은 채널의 두 번째 장기)

93. Miki T, Liss B, Minami K, et al. **ATP-sensitive K⁺ channels in the hypothalamus are essential for the maintenance of glucose homeostasis.** Nat Neurosci. 2001;4(5):507–12. <https://pubmed.ncbi.nlm.nih.gov/11319559/>
94. Evans ML, McCrimmon RJ, Flanagan DE, et al. **Hypothalamic ATP-sensitive K⁺ channels play a key role in sensing hypoglycemia and triggering counterregulatory epinephrine and glucagon responses.** Diabetes. 2004;53(10):2542–51. <https://pubmed.ncbi.nlm.nih.gov/15448088/>

> 모델 의존: 지도의 점선 화살표 `g → VMH`. 같은 채널이 시상하부에도 있으므로
> K_ATP-CHI 환자는 **저혈당을 감지하는 기관 자체가 돌연변이되어 있다**.
> 모델은 이 항을 정량화하지 않고 구조로만 남겨두었다(정직성 표기).

---

## 14. 도구 · 방법론 (QSP / mrgsolve)

95. Baron KT, Elmokadem A, et al. **mrgsolve: Simulate from ODE-Based Models.** <https://mrgsolve.org/>
96. Bergman RN, Ider YZ, Bowden CR, Cobelli C. **Quantitative estimation of insulin sensitivity.** Am J Physiol. 1979;236(6):E667–77. <https://pubmed.ncbi.nlm.nih.gov/443421/>
97. Dalla Man C, Rizza RA, Cobelli C. **Meal simulation model of the glucose-insulin system.** IEEE Trans Biomed Eng. 2007;54(10):1740–9. <https://pubmed.ncbi.nlm.nih.gov/17926672/>
98. Topp B, Promislow K, deVries G, et al. **A model of beta-cell mass, insulin, and glucose kinetics: pathways to diabetes.** J Theor Biol. 2000;206(4):605–19. <https://pubmed.ncbi.nlm.nih.gov/?term=Topp+model+beta-cell+mass+insulin+glucose+kinetics>
99. Fridlyand LE, Philipson LH. **Glucose sensing in the pancreatic beta cell: a computational systems analysis.** Theor Biol Med Model. 2010;7:15. <https://pubmed.ncbi.nlm.nih.gov/20497556/>
100. Chay TR, Keizer J. **Minimal model for membrane oscillations in the pancreatic beta-cell.** Biophys J. 1983;42(2):181–90. <https://pubmed.ncbi.nlm.nih.gov/6305437/>

> 모델 의존: 96·97번이 전신 포도당–인슐린 구획 구조의 계보,
> 98번이 베타세포 질량(B)을 상태변수로 두는 선택의 계보,
> 99·100번이 전압 분배기와 `Popen → V_m → Ca` 사슬의 계보다.
> 이 모델의 새로운 부분은 그 두 계보를 **하나의 파라미터 g**로 잇고,
> 약물을 g에 대한 의존성(곱셈/덧셈/무관)으로 분류한 점이다.

---

## 부록 — 링크가 불안정한 항목

18·48·64·98번은 PMID가 확실하지 않거나 개발 단계 자료여서 검색 URL 또는
불완전한 서지로 남겼다. 해당 문헌은 모델의 **구조**(트래피킹 결함, 신생아 디아족사이드 PK,
인슐린수용체 항체, 베타세포 질량 모델링)를 뒷받침하지만, 어떤 **보정 숫자**도
이들에서 가져오지 않았다.
