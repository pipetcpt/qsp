# Dementia with Lewy Bodies (DLB) — Reference List

지원 문헌 목록 (Supporting literature for `dlb_qsp_model.dot`,
`dlb_mrgsolve_model.R`, `dlb_shiny_app.R`).

**링크 규약 (link convention).** 아래 모든 항목은 **제목으로 해석되는 PubMed
검색 URL**(`pubmed.ncbi.nlm.nih.gov/?term=...`)로 연결됩니다. 고정 PMID를 직접
쓰지 않은 이유는, 자동 생성 세션에서 잘못된 PMID가 섞여 들어가면 오히려 추적이
어려워지기 때문입니다. 각 링크는 해당 논문(또는 그 논문이 최상위로 나오는
검색 결과)으로 바로 이동합니다.

*Every entry below links to a **title-resolving PubMed query**. Fixed PMIDs are
deliberately not asserted: a single mistyped accession number is worse than a
query that always resolves. Each link lands on the paper (or a result set with
it first).*

각 항목 뒤의 `[모델 위치]` 표시는 그 문헌이 모델의 어느 부분을 지지하는지
가리킵니다 — 파라미터, 방정식, 또는 보정 기준점.

---

## 1. 진단 기준과 임상 표현형 (Diagnostic criteria & clinical phenotype)

1. McKeith IG, Boeve BF, Dickson DW, et al. **Diagnosis and management of dementia with Lewy bodies: Fourth consensus report of the DLB Consortium.** *Neurology* 2017. — 핵심 임상 특징(인지 변동·환시·REM 수면 행동장애·파킨슨증), 지표 바이오마커(DaTSCAN·MIBG·PSG). 모델의 임상 엔드포인트 정의 전체의 근거. [`$INIT` 임상 상태, `$TABLE`]
   <https://pubmed.ncbi.nlm.nih.gov/?term=Diagnosis+and+management+of+dementia+with+Lewy+bodies+Fourth+consensus+report+of+the+DLB+Consortium>
2. McKeith IG, Ferman TJ, Thomas AJ, et al. **Research criteria for the diagnosis of prodromal dementia with Lewy bodies.** *Neurology* 2020. — 전구기(MCI-LB·정신과적 발현·섬망성 발현) 정의. 시나리오 19의 "전구기 투여" 개념적 근거.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Research+criteria+for+the+diagnosis+of+prodromal+dementia+with+Lewy+bodies>
3. Walker Z, Possin KL, Boeve BF, Aarsland D. **Lewy body dementias.** *Lancet* 2015. — DLB와 파킨슨병 치매(PDD)를 하나의 루이체 질환 스펙트럼으로 보는 관점. 모델의 `PHENO` 스위치 설계 근거.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Lewy+body+dementias+Walker+Possin+Boeve+Aarsland+Lancet>
4. Ferman TJ, Smith GE, Boeve BF, et al. **DLB fluctuations: specific features that reliably differentiate DLB from AD and normal aging.** *Neurology* 2004. — 인지 변동을 AD와 구별하는 4항목 척도(Mayo Fluctuations Composite). 모델이 `CAF`를 "평균이 아닌 분산"으로 두는 근거. [`FLUCTGT`]
   <https://pubmed.ncbi.nlm.nih.gov/?term=DLB+fluctuations+specific+features+that+reliably+differentiate+DLB+from+AD+and+normal+aging>
5. Walker MP, Ayre GA, Cummings JL, et al. **The Clinician Assessment of Fluctuation and the One Day Fluctuation Assessment Scale.** *Br J Psychiatry* 2000. — CAF 척도(0–16)의 원 논문. 모델 `CAF = 4 × FLUC` 스케일링의 근거.
   <https://pubmed.ncbi.nlm.nih.gov/?term=The+Clinician+Assessment+of+Fluctuation+and+the+One+Day+Fluctuation+Assessment+Scale>
6. Vann Jones SA, O'Brien JT. **The prevalence and incidence of dementia with Lewy bodies: a systematic review of population and clinical studies.** *Psychol Med* 2014. — 유병률·발생률, 남성 우세.
   <https://pubmed.ncbi.nlm.nih.gov/?term=The+prevalence+and+incidence+of+dementia+with+Lewy+bodies+a+systematic+review+of+population+and+clinical+studies>
7. Kane JPM, Surendranathan A, Bentley A, et al. **Clinical prevalence of Lewy body dementia.** *Alzheimers Res Ther* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+prevalence+of+Lewy+body+dementia+Kane+Surendranathan>
8. Emre M, Aarsland D, Brown R, et al. **Clinical diagnostic criteria for dementia associated with Parkinson's disease.** *Mov Disord* 2007. — PDD 기준과 "1년 규칙". `PHENO = 1` 팔의 정의.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+diagnostic+criteria+for+dementia+associated+with+Parkinson+disease+Emre+Aarsland>

## 2. 신경병리와 병기 (Neuropathology & staging)

9. Braak H, Del Tredici K, Rüb U, et al. **Staging of brain pathology related to sporadic Parkinson's disease.** *Neurobiol Aging* 2003. — 뇌간→변연계→신피질 상행 병기. 모델의 3구획 지역 구조(`FIBB/FIBL/FIBN`)와 순차적 게이팅(`UPTL`, `UPTN`).
   <https://pubmed.ncbi.nlm.nih.gov/?term=Staging+of+brain+pathology+related+to+sporadic+Parkinson+disease+Braak+Del+Tredici>
10. Attems J, Toledo JB, Walker L, et al. **Neuropathological consensus criteria for the evaluation of Lewy pathology in post-mortem brains.** *Acta Neuropathol* 2021. — 뇌간우세/변연계/신피질 분포 등급화.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathological+consensus+criteria+for+the+evaluation+of+Lewy+pathology+in+post-mortem+brains>
11. Beach TG, Adler CH, Lue L, et al. **Unified staging system for Lewy body disorders.** *Acta Neuropathol* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Unified+staging+system+for+Lewy+body+disorders+Beach+Adler+Lue>
12. Irwin DJ, Grossman M, Weintraub D, et al. **Neuropathological and genetic correlates of survival and dementia onset in synucleinopathies.** *Lancet Neurol* 2017. — 공존 AD 병리가 생존을 단축시킴. 모델의 `B_CROSS`/`BNS` 및 `PTAU` → 생존 경로.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathological+and+genetic+correlates+of+survival+and+dementia+onset+in+synucleinopathies>
13. Borghammer P, Van Den Berge N. **Brain-first versus gut-first Parkinson's disease: a hypothesis.** *J Parkinsons Dis* 2019. — `S_GUT` / `S_OLF` 두 진입 경로.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Brain-first+versus+gut-first+Parkinson+disease+a+hypothesis+Borghammer>
14. Halliday GM, Holton JL, Revesz T, Dickson DW. **Neuropathology underlying clinical variability in patients with synucleinopathies.** *Acta Neuropathol* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuropathology+underlying+clinical+variability+in+patients+with+synucleinopathies>

## 3. α-시누클레인 생물물리학·응집 동역학·스트레인 (α-Synuclein biophysics, aggregation kinetics, strains)

15. Cohen SIA, Linse S, Luheshi LM, et al. **Proliferation of amyloid-β42 aggregates occurs through a secondary nucleation mechanism.** *PNAS* 2013. — 2차 핵형성(섬유 표면 촉매)의 정량적 틀. 모델 `KNUC2 · ASYNM · FIB`의 형태 그 자체.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Proliferation+of+amyloid-beta42+aggregates+occurs+through+a+secondary+nucleation+mechanism>
16. Buell AK, Galvagnion C, Gaspar R, et al. **Solution conditions determine the relative importance of nucleation and growth processes in α-synuclein aggregation.** *PNAS* 2014. — α-시누클레인의 1차/2차 핵형성 상대 기여.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Solution+conditions+determine+the+relative+importance+of+nucleation+and+growth+processes+in+alpha-synuclein+aggregation>
17. Meisl G, Kirkegaard JB, Arosio P, et al. **Molecular mechanisms of protein aggregation from global fitting of kinetic models.** *Nat Protoc* 2016. — 응집 동역학 모델의 표준 형식(용량 제한 항 포함). 모델의 `(1 − OLIG)`, `(1 − FIB)` 포화 항.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Molecular+mechanisms+of+protein+aggregation+from+global+fitting+of+kinetic+models>
18. Winner B, Jappelli R, Maji SK, et al. **In vivo demonstration that α-synuclein oligomers are toxic.** *PNAS* 2011. — 올리고머(섬유가 아니라)가 독성 종. 모델이 신경 소실을 `OLIG`에, 전파를 `FIB`에 연결하는 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=In+vivo+demonstration+that+alpha-synuclein+oligomers+are+toxic>
19. Schweighauser M, Shi Y, Tarutani A, et al. **Structures of α-synuclein filaments from multiple system atrophy.** *Nature* 2020. — 질환별 서로 다른 섬유 구조(cryo-EM). 모델 `STRAIN` / `A_STRAIN` 파라미터.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Structures+of+alpha-synuclein+filaments+from+multiple+system+atrophy>
20. Yang Y, Shi Y, Schweighauser M, et al. **Structures of α-synuclein filaments from human brains with Lewy pathology.** *Nature* 2022. — 루이 병리 고유 폴드.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Structures+of+alpha-synuclein+filaments+from+human+brains+with+Lewy+pathology>
21. Peng C, Gathagan RJ, Covell DJ, et al. **Cellular milieu imparts distinct pathological α-synuclein strains in α-synucleinopathies.** *Nature* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cellular+milieu+imparts+distinct+pathological+alpha-synuclein+strains+in+alpha-synucleinopathies>
22. Fujiwara H, Hasegawa M, Dohmae N, et al. **α-Synuclein is phosphorylated in synucleinopathy lesions.** *Nat Cell Biol* 2002. — pS129가 루이체 α-시누클레인의 >90%. `A_PS129`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+is+phosphorylated+in+synucleinopathy+lesions>
23. Shahmoradian SH, Lewis AJ, Genoud C, et al. **Lewy pathology in Parkinson's disease consists of crowded organelles and lipid membranes.** *Nat Neurosci* 2019. — 루이체가 단순 단백질 응집체가 아니라 막·소기관 집합체. `A_LB`를 "불활성 sink 또는 활성 저장소"로 둔 이유.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lewy+pathology+in+Parkinson+disease+consists+of+crowded+organelles+and+lipid+membranes>
24. Volpicelli-Daley LA, Luk KC, Patel TP, et al. **Exogenous α-synuclein fibrils induce Lewy body pathology leading to synaptic dysfunction and neuron death.** *Neuron* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Exogenous+alpha-synuclein+fibrils+induce+Lewy+body+pathology+leading+to+synaptic+dysfunction+and+neuron+death>

## 4. 세포 간 전파 (Cell-to-cell propagation)

25. Desplats P, Lee HJ, Bae EJ, et al. **Inclusion formation and neuronal cell death through neuron-to-neuron transmission of α-synuclein.** *PNAS* 2009. — `S_EXO` → `S_ISF` → `S_UPT` 경로.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Inclusion+formation+and+neuronal+cell+death+through+neuron-to-neuron+transmission+of+alpha-synuclein>
26. Mao X, Ou MT, Karuppagounder SS, et al. **Pathological α-synuclein transmission initiated by binding lymphocyte-activation gene 3.** *Science* 2016. — LAG3 수용체 매개 흡수.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Pathological+alpha-synuclein+transmission+initiated+by+binding+lymphocyte-activation+gene+3>
27. Holmes BB, DeVos SL, Kfoury N, et al. **Heparan sulfate proteoglycans mediate internalization and propagation of specific proteopathic seeds.** *PNAS* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Heparan+sulfate+proteoglycans+mediate+internalization+and+propagation+of+specific+proteopathic+seeds>
28. Kim S, Kwon SH, Kam TI, et al. **Transneuronal propagation of pathologic α-synuclein from the gut to the brain models Parkinson's disease.** *Neuron* 2019. — 위·미주신경 경로.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Transneuronal+propagation+of+pathologic+alpha-synuclein+from+the+gut+to+the+brain+models+Parkinson+disease>
29. Xie L, Kang H, Xu Q, et al. **Sleep drives metabolite clearance from the adult brain.** *Science* 2013. — 수면 의존 글림파틱 청소. 모델 `GLYMF = f(SLD)` — RBD 회로 손상이 seed 청소를 악화시키는 되먹임.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sleep+drives+metabolite+clearance+from+the+adult+brain>
30. Iliff JJ, Wang M, Liao Y, et al. **A paravascular pathway facilitates CSF flow through the brain parenchyma and the clearance of interstitial solutes.** *Sci Transl Med* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+paravascular+pathway+facilitates+CSF+flow+through+the+brain+parenchyma+and+the+clearance+of+interstitial+solutes>

## 5. GBA1 · 리소좀 · 앰브록솔 (GBA1, lysosome, ambroxol)

31. Mazzulli JR, Xu YH, Sun Y, et al. **Gaucher disease glucocerebrosidase and α-synuclein form a bidirectional pathogenic loop in synucleinopathies.** *Cell* 2011. — **모델의 GBA1 양안정 스위치 그 자체.** GCase↓ → GlcCer↑ → 올리고머 안정화 → GCase 수송 차단 → GCase↓. [`L_SWITCH`, `GCTGT`, `KOLDE`]
    <https://pubmed.ncbi.nlm.nih.gov/?term=Gaucher+disease+glucocerebrosidase+and+alpha-synuclein+form+a+bidirectional+pathogenic+loop+in+synucleinopathies>
32. Nalls MA, Duran R, Lopez G, et al. **A multicenter study of glucocerebrosidase mutations in dementia with Lewy bodies.** *JAMA Neurol* 2013. — GBA1 변이의 DLB 오즈비 ≈ 8, 파킨슨병보다 더 강함. `GBAF` 파라미터.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+multicenter+study+of+glucocerebrosidase+mutations+in+dementia+with+Lewy+bodies>
33. Sidransky E, Nalls MA, Aasly JO, et al. **Multicenter analysis of glucocerebrosidase mutations in Parkinson's disease.** *N Engl J Med* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Multicenter+analysis+of+glucocerebrosidase+mutations+in+Parkinson+disease+Sidransky+Nalls>
34. Migdalska-Richards A, Daly L, Bezard E, Schapira AHV. **Ambroxol effects in glucocerebrosidase and α-synuclein transgenic mice.** *Ann Neurol* 2016. — 앰브록솔이 CNS GCase 활성을 높이고 α-시누클레인을 낮춤.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+effects+in+glucocerebrosidase+and+alpha-synuclein+transgenic+mice>
35. Mullin S, Smith L, Lee K, et al. **Ambroxol for the treatment of patients with Parkinson disease with and without glucocerebrosidase gene mutations: a non-randomized, noncontrolled trial.** *JAMA Neurol* 2020. — 1.26 g/일, 백혈구 GCase 활성 약 +35%, CSF 도달 확인. 모델 `EMAXAMB`/`EC50AMB` 보정 기준점.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+for+the+treatment+of+patients+with+Parkinson+disease+with+and+without+glucocerebrosidase+gene+mutations>
36. Silveira CRA, MacKinley J, Coleman K, et al. **Ambroxol as a novel disease-modifying treatment for Parkinson's disease dementia: protocol for a single-centre, randomized, double-blind, placebo-controlled trial.** *BMC Neurol* 2019. — PDD를 표적으로 한 앰브록솔 시험 설계.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ambroxol+as+a+novel+disease-modifying+treatment+for+Parkinson+disease+dementia+protocol>
37. Cuervo AM, Stefanis L, Fredenburg R, et al. **Impaired degradation of mutant α-synuclein by chaperone-mediated autophagy.** *Science* 2004. — CMA/LAMP2A. `L_CMA`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Impaired+degradation+of+mutant+alpha-synuclein+by+chaperone-mediated+autophagy>
38. Jinn S, Drolet RE, Cramer PE, et al. **TMEM175 deficiency impairs lysosomal and mitochondrial function and increases α-synuclein aggregation.** *PNAS* 2017. — `G_TMEM` → `L_ACID`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=TMEM175+deficiency+impairs+lysosomal+and+mitochondrial+function+and+increases+alpha-synuclein+aggregation>
39. Zunke F, Moise AC, Belur NR, et al. **Reversible conformational conversion of α-synuclein into toxic assemblies by glucosylceramide.** *Neuron* 2018. — GlcCer가 올리고머 형태를 안정화. `KGLC50`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Reversible+conformational+conversion+of+alpha-synuclein+into+toxic+assemblies+by+glucosylceramide>
40. Peterschmitt MJ, Saiki H, Hatano T, et al. **Safety, pharmacokinetics, and pharmacodynamics of oral venglustat in Parkinson disease patients with a GBA mutation.** *Mov Disord* 2022. — GCS 억제제(상류 차단). `Y_VENG`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Safety+pharmacokinetics+and+pharmacodynamics+of+oral+venglustat+in+Parkinson+disease+patients+with+a+GBA+mutation>

## 6. GBA1 외 유전학 (Genetics beyond GBA1)

41. Guerreiro R, Ross OA, Kun-Rodrigues C, et al. **Investigating the genetic architecture of dementia with Lewy bodies: a two-stage genome-wide association study.** *Lancet Neurol* 2018. — APOE, SNCA, GBA1이 DLB의 세 주요 좌위.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Investigating+the+genetic+architecture+of+dementia+with+Lewy+bodies+a+two-stage+genome-wide+association+study>
42. Chia R, Sabir MS, Bandres-Ciga S, et al. **Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture.** *Nat Genet* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Genome+sequencing+analysis+identifies+new+loci+associated+with+Lewy+body+dementia+and+provides+insights+into+its+genetic+architecture>
43. Dickson DW, Heckman MG, Murray ME, et al. **APOE ε4 is associated with severity of Lewy body pathology independent of Alzheimer pathology.** *Neurology* 2018. — APOE4가 AD 병리와 독립적으로 루이 병리를 악화. 모델의 `APOE4` → `A_FIB` 점선 경로.
    <https://pubmed.ncbi.nlm.nih.gov/?term=APOE+epsilon4+is+associated+with+severity+of+Lewy+body+pathology+independent+of+Alzheimer+pathology>
44. Singleton AB, Farrer M, Johnson J, et al. **α-Synuclein locus triplication causes Parkinson's disease.** *Science* 2003. — 유전자 용량 효과. `SNCADOSE`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+locus+triplication+causes+Parkinson+disease>

## 7. 콜린성 계 — 전시냅스 소실, 후시냅스 보존 (Cholinergic system)

45. Perry EK, Marshall E, Kerwin J, et al. **Evidence of a monoaminergic-cholinergic imbalance related to visual hallucinations in Lewy body dementia.** *J Neurochem* 1990. — 환시와 콜린성/모노아민 불균형.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Evidence+of+a+monoaminergic-cholinergic+imbalance+related+to+visual+hallucinations+in+Lewy+body+dementia>
46. Tiraboschi P, Hansen LA, Alford M, et al. **Cholinergic dysfunction in diseases with Lewy bodies.** *Neurology* 2000. — **DLB의 피질 ChAT 감소가 AD보다 크다.** 모델 `NBM_0`(DLB 0.50 vs AD 0.72)의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cholinergic+dysfunction+in+diseases+with+Lewy+bodies+Tiraboschi+Hansen+Alford>
47. Perry EK, Haroutunian V, Davis KL, et al. **Neocortical cholinergic activities differentiate Lewy body dementia from classical Alzheimer's disease.** *Neuroreport* 1994.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neocortical+cholinergic+activities+differentiate+Lewy+body+dementia+from+classical+Alzheimer+disease>
48. Perry EK, Marshall E, Perry RH, et al. **Cholinergic and dopaminergic activities in senile dementia of Lewy body type.** *Alzheimer Dis Assoc Disord* 1990.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cholinergic+and+dopaminergic+activities+in+senile+dementia+of+Lewy+body+type>
49. Shiozaki K, Iseki E, Uchiyama H, et al. **Alterations of muscarinic acetylcholine receptor subtypes in diffuse Lewy body disease: relation to Alzheimer's disease.** *J Neurol Neurosurg Psychiatry* 1999. — **DLB에서 후시냅스 M1이 상대적으로 보존된다** — 모델에서 `M1TGT = M1BASE − KTAUM1·PTAU`가 tau에만 의존하게 만든 직접 근거. ChEI가 AD보다 DLB에서 더 잘 듣는 이유.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Alterations+of+muscarinic+acetylcholine+receptor+subtypes+in+diffuse+Lewy+body+disease>
50. Ballard C, Piggott M, Johnson M, et al. **Delusions associated with elevated muscarinic binding in dementia with Lewy bodies.** *Ann Neurol* 2000. — 무스카린 결합이 오히려 상승.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Delusions+associated+with+elevated+muscarinic+binding+in+dementia+with+Lewy+bodies>
51. Pimlott SL, Piggott M, Owens J, et al. **Nicotinic acetylcholine receptor distribution in Alzheimer's disease, dementia with Lewy bodies, Parkinson's disease, and vascular dementia.** *Neuropsychopharmacology* 2004. — α4β2 nAChR 시상·피질 감소. `C_NIC`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nicotinic+acetylcholine+receptor+distribution+in+Alzheimer+disease+dementia+with+Lewy+bodies+Parkinson+disease+and+vascular+dementia>
52. Mesulam MM, Geula C. **Butyrylcholinesterase-rich neurons of the human cerebral cortex.** *Ann Neurol* 1988. — BuChE의 피질 분포. 모델 `BCHEFR`/`KBCHUP`: AChE 소실에 따라 BuChE가 가수분해의 더 큰 몫을 담당하게 되어, 이중 억제제(리바스티그민)가 후기에 상대적 이점을 갖는다는 예측의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Butyrylcholinesterase-rich+neurons+of+the+human+cerebral+cortex+Mesulam+Geula>
53. Greig NH, Utsuki T, Ingram DK, et al. **Selective butyrylcholinesterase inhibition elevates brain acetylcholine, augments learning and lowers Alzheimer β-amyloid peptide in rodent.** *PNAS* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Selective+butyrylcholinesterase+inhibition+elevates+brain+acetylcholine+augments+learning+and+lowers+Alzheimer+beta-amyloid+peptide>

## 8. 도파민성 계 · DAT 영상 · 신경이완제 민감성 (Dopaminergic system, DAT imaging, neuroleptic sensitivity)

54. McKeith I, Fairbairn A, Perry R, et al. **Neuroleptic sensitivity in patients with senile dementia of Lewy body type.** *BMJ* 1992. — **원 논문.** 심한 신경이완제 민감성 반응이 약 절반에서 발생하고 사망률이 2–3배. 모델 `NSENS`의 보정 목표.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuroleptic+sensitivity+in+patients+with+senile+dementia+of+Lewy+body+type+McKeith+Fairbairn+Perry>
55. Aarsland D, Perry R, Larsen JP, et al. **Neuroleptic sensitivity in Parkinson's disease and parkinsonian dementias.** *J Clin Psychiatry* 2005. — **PDD에서도 같은 민감성이 나타난다.** 모델이 DLB와 PDD를 유사하게(그리고 AD와는 극명히 다르게) 예측하는 것이 옳은 이유.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Neuroleptic+sensitivity+in+Parkinson+disease+and+parkinsonian+dementias+Aarsland+Perry+Larsen>
56. Piggott MA, Marshall EF, Thomas N, et al. **Striatal dopaminergic markers in dementia with Lewy bodies, Alzheimer's and Parkinson's diseases: rostrocaudal distribution.** *Brain* 1999. — **DLB에서는 선조체 D2 수용체가 상향조절되지 않는다** (PD에서는 된다). 모델의 `UPCAP = exp(−KSUPP·OLIGL)` 항의 직접 근거이며, 이 한 항이 신경이완제 민감성과 레보도파 반응 저하를 동시에 만들어 낸다.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Striatal+dopaminergic+markers+in+dementia+with+Lewy+bodies+Alzheimer+and+Parkinson+diseases+rostrocaudal+distribution>
57. Piggott MA, Ballard CG, Rowan E, et al. **Selective loss of dopamine D2 receptors in temporal cortex in dementia with Lewy bodies, association with cognitive decline.** *Synapse* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Selective+loss+of+dopamine+D2+receptors+in+temporal+cortex+in+dementia+with+Lewy+bodies>
58. McKeith I, O'Brien J, Walker Z, et al. **Sensitivity and specificity of dopamine transporter imaging with ¹²³I-FP-CIT SPECT in dementia with Lewy bodies: a phase III, multicentre study.** *Lancet Neurol* 2007. — DaTSCAN 진단 성능. 모델 `DATSBR` 매핑과 `TERMEXP` 지수.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sensitivity+and+specificity+of+dopamine+transporter+imaging+with+123I-FP-CIT+SPECT+in+dementia+with+Lewy+bodies+phase+III>
59. Walker Z, Costa DC, Walker RWH, et al. **Differentiation of dementia with Lewy bodies from Alzheimer's disease using a dopaminergic presynaptic ligand.** *J Neurol Neurosurg Psychiatry* 2002.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Differentiation+of+dementia+with+Lewy+bodies+from+Alzheimer+disease+using+a+dopaminergic+presynaptic+ligand>
60. Kaasinen V, Vahlberg T. **Striatal dopamine in Parkinson disease: a meta-analysis of imaging studies.** *Ann Neurol* 2017. — 선조체 종말 소실이 세포체 소실보다 앞서고 크다. 모델 `TERMEXP = 1.6`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Striatal+dopamine+in+Parkinson+disease+a+meta-analysis+of+imaging+studies+Kaasinen+Vahlberg>
61. Farde L, Nordström AL, Wiesel FA, et al. **Positron emission tomographic analysis of central D1 and D2 dopamine receptor occupancy in patients treated with classical neuroleptics and clozapine.** *Arch Gen Psychiatry* 1992. — D2 점유율 ~80%가 추체외로 증상 문턱. 모델의 `EC50D2` 및 운동 예비능 시그모이드 위치.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Positron+emission+tomographic+analysis+of+central+D1+and+D2+dopamine+receptor+occupancy+in+patients+treated+with+classical+neuroleptics+and+clozapine>
62. Kapur S, Zipursky R, Jones C, et al. **Relationship between dopamine D2 occupancy, clinical response, and side effects: a double-blind PET study of first-episode schizophrenia.** *Am J Psychiatry* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Relationship+between+dopamine+D2+occupancy+clinical+response+and+side+effects+a+double-blind+PET+study+of+first-episode+schizophrenia>
63. Kegeles LS, Slifstein M, Frankle WG, et al. **Dose-occupancy study of striatal and extrastriatal dopamine D2 receptors by aripiprazole in schizophrenia with PET and [¹⁸F]fallypride.** *Neuropsychopharmacology* 2008.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dose-occupancy+study+of+striatal+and+extrastriatal+dopamine+D2+receptors+by+aripiprazole+in+schizophrenia>

## 9. 세로토닌성 계와 5-HT2A (Serotonergic system and 5-HT2A)

64. Ballard C, Piggott M, Johnson M, et al. **5-HT2A receptor binding in the ventral visual pathway and visual hallucinations in dementia with Lewy bodies.** *(and companion analyses)* — **환시가 있는 DLB에서 BA18/19의 5-HT2A 결합이 증가한다.** 모델 `H2TGT = 1 + KH2UP·denervation + KH2NEO·FIBN`의 근거이며, 세 전달자 중 유일하게 후시냅스 수용체가 **상향**조절되는 계.
    <https://pubmed.ncbi.nlm.nih.gov/?term=5-HT2A+receptor+binding+ventral+visual+pathway+visual+hallucinations+dementia+with+Lewy+bodies>
65. Cheng AV, Ferrier IN, Morris CM, et al. **Cortical serotonin-S2 receptor binding in Lewy body dementia, Alzheimer's and Parkinson's diseases.** *J Neurol Sci* 1991.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cortical+serotonin-S2+receptor+binding+in+Lewy+body+dementia+Alzheimer+and+Parkinson+diseases>
66. Sharp SI, Ballard CG, Chen CPL-H, Francis PT. **Aggressive behavior and neuroleptic medication are associated with increased number of α1-adrenoceptors in patients with Lewy body disease.** *Am J Geriatr Psychiatry* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Aggressive+behavior+and+neuroleptic+medication+are+associated+with+increased+number+of+alpha1-adrenoceptors+in+patients+with+Lewy+body+disease>
67. Vanover KE, Davis RE. **Role of 5-HT2A receptor antagonists in the treatment of insomnia.** *Nat Sci Sleep* 2010. — 5-HT2A 역작용제 약리(구성적 활성 개념). 모델의 `HTCONST`(작용제 비의존 Gq 톤) 항.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+5-HT2A+receptor+antagonists+in+the+treatment+of+insomnia+Vanover+Davis>
68. Weiner DM, Burstein ES, Nash N, et al. **5-Hydroxytryptamine2A receptor inverse agonists as antipsychotics.** *J Pharmacol Exp Ther* 2001. — 피마반세린 계열의 약리학적 기초.
    <https://pubmed.ncbi.nlm.nih.gov/?term=5-Hydroxytryptamine2A+receptor+inverse+agonists+as+antipsychotics+Weiner+Burstein+Nash>

## 10. 노르아드레날린성 계와 청반 (Noradrenergic system, locus coeruleus)

69. Del Tredici K, Braak H. **Dysfunction of the locus coeruleus–norepinephrine system and related circuitry in Parkinson's disease-related dementia.** *J Neurol Neurosurg Psychiatry* 2013. — 청반이 가장 이른 표적 중 하나. 모델의 `LC_0`(DLB 0.62 vs AD 0.80)와 `KLLC`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dysfunction+of+the+locus+coeruleus-norepinephrine+system+and+related+circuitry+in+Parkinson+disease-related+dementia>
70. Szot P, White SS, Greenup JL, et al. **Compensatory changes in the noradrenergic nervous system in the locus ceruleus and hippocampus of postmortem subjects with Alzheimer's disease and dementia with Lewy bodies.** *J Neurosci* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Compensatory+changes+in+the+noradrenergic+nervous+system+in+the+locus+ceruleus+and+hippocampus+of+postmortem+subjects+with+Alzheimer+disease+and+dementia+with+Lewy+bodies>
71. Aston-Jones G, Cohen JD. **An integrative theory of locus coeruleus–norepinephrine function: adaptive gain and optimal performance.** *Annu Rev Neurosci* 2005. — LC-NE의 이득 조절 이론. 모델의 `NOISEG`(LC+PPN 손실이 주의 상태의 전환 소음을 키운다)의 개념적 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=An+integrative+theory+of+locus+coeruleus-norepinephrine+function+adaptive+gain+and+optimal+performance>

## 11. 인지 변동 · EEG · 시상피질 회로 (Cognitive fluctuation, EEG, thalamocortical circuitry)

72. Bonanni L, Thomas A, Tiraboschi P, et al. **EEG comparisons in early Alzheimer's disease, dementia with Lewy bodies and Parkinson's disease with dementia patients with a 2-year follow-up.** *Brain* 2008. — DLB의 우세 주파수 저하와 **높은 시간적 변동성**(pre-alpha 5.6–7.9 Hz). 모델 `EEGF` 매핑.
    <https://pubmed.ncbi.nlm.nih.gov/?term=EEG+comparisons+in+early+Alzheimer+disease+dementia+with+Lewy+bodies+and+Parkinson+disease+with+dementia+patients+with+a+2-year+follow-up>
73. Bonanni L, Franciotti R, Nobili F, et al. **EEG markers of dementia with Lewy bodies: a multicenter cohort study.** *J Alzheimers Dis* 2016.
    <https://pubmed.ncbi.nlm.nih.gov/?term=EEG+markers+of+dementia+with+Lewy+bodies+a+multicenter+cohort+study>
74. Walker MP, Ayre GA, Cummings JL, et al. **Quantifying fluctuation in dementia with Lewy bodies, Alzheimer's disease, and vascular dementia.** *Neurology* 2000. — 변동이 DLB에서 특징적으로 크다는 정량적 비교.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Quantifying+fluctuation+in+dementia+with+Lewy+bodies+Alzheimer+disease+and+vascular+dementia>
75. Ballard C, Walker M, O'Brien J, et al. **The characterisation and impact of 'fluctuating' cognition in dementia with Lewy bodies and Alzheimer's disease.** *Int J Geriatr Psychiatry* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+characterisation+and+impact+of+fluctuating+cognition+in+dementia+with+Lewy+bodies+and+Alzheimer+disease>
76. Delli Pizzi S, Franciotti R, Taylor JP, et al. **Thalamic involvement in fluctuating cognition in dementia with Lewy bodies: magnetic resonance evidences.** *Cereb Cortex* 2015. — 시상 관여. 모델 `E_THAL` / `WTHAL`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Thalamic+involvement+in+fluctuating+cognition+in+dementia+with+Lewy+bodies+magnetic+resonance+evidences>
77. Peraza LR, Cromarty R, Kobeleva X, et al. **Electroencephalographic derived network differences in Lewy body dementia compared to Alzheimer's disease patients.** *Sci Rep* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Electroencephalographic+derived+network+differences+in+Lewy+body+dementia+compared+to+Alzheimer+disease+patients>
78. Deco G, Jirsa VK, McIntosh AR. **Emerging concepts for the dynamical organization of resting-state activity in the brain.** *Nat Rev Neurosci* 2011. — 뇌 상태의 다중안정성과 임계 근방 동역학. 모델의 3차 함수 주의 상태(`dxdt_ATTM`)와 "분산이 평균보다 먼저 커진다"는 구조적 주장의 이론적 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emerging+concepts+for+the+dynamical+organization+of+resting-state+activity+in+the+brain>
79. Scheffer M, Bascompte J, Brock WA, et al. **Early-warning signals for critical transitions.** *Nature* 2009. — 임계 전이 접근 시 분산·자기상관 증가(critical slowing down). `CURV`, `BISTAB`, `FLUCTGT`의 형태.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Early-warning+signals+for+critical+transitions+Scheffer+Bascompte+Brock>

## 12. 환시와 PAD 모형 (Visual hallucinations and the PAD model)

80. Collerton D, Perry E, McKeith I. **Why people see things that are not there: a novel Perception and Attention Deficit model for recurrent complex visual hallucinations.** *Behav Brain Sci* 2005. — **모델의 PAD 곱 구조(`(1−BOTTOMUP)·(1−TOPDOWN)·5HT2A`)의 원 출처.** 지각 결손과 주의 결손이 **동시에** 필요하다는 주장.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Why+people+see+things+that+are+not+there+a+novel+Perception+and+Attention+Deficit+model+for+recurrent+complex+visual+hallucinations>
81. Shine JM, Halliday GM, Naismith SL, Lewis SJG. **Visual misperceptions and hallucinations in Parkinson's disease: dysfunction of attentional control networks?** *Mov Disord* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+misperceptions+and+hallucinations+in+Parkinson+disease+dysfunction+of+attentional+control+networks>
82. Onofrj M, Taylor JP, Monaco D, et al. **Visual hallucinations in PD and Lewy body dementias: old and new hypotheses.** *Behav Neurol* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+hallucinations+in+PD+and+Lewy+body+dementias+old+and+new+hypotheses>
83. Taylor JP, Firbank MJ, He J, et al. **Visual cortex in dementia with Lewy bodies: magnetic resonance imaging study.** *Br J Psychiatry* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Visual+cortex+in+dementia+with+Lewy+bodies+magnetic+resonance+imaging+study>
84. Lobotesis K, Fenwick JD, Phipps A, et al. **Occipital hypoperfusion on SPECT in dementia with Lewy bodies but not AD.** *Neurology* 2001. — 후두엽 저대사/저관류. 모델 `BOTTOMUP` 감소 항.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Occipital+hypoperfusion+on+SPECT+in+dementia+with+Lewy+bodies+but+not+AD>
85. Lim SM, Katsifis A, Villemagne VL, et al. **The ¹⁸F-FDG PET cingulate island sign and comparison to ¹²³I-β-CIT SPECT for diagnosis of dementia with Lewy bodies.** *J Nucl Med* 2009. — cingulate island sign.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+18F-FDG+PET+cingulate+island+sign+and+comparison+to+123I-beta-CIT+SPECT+for+diagnosis+of+dementia+with+Lewy+bodies>
86. Murakami H, Shiraishi T, Umehara T, et al. **Retinal thinning and visual hallucination in dementia with Lewy bodies** *(OCT studies of inner retinal layers)*. — 망막 층 얇아짐. `V_RET`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=retinal+thinning+optical+coherence+tomography+dementia+with+Lewy+bodies+visual+hallucination>

## 13. REM 수면 행동장애와 수면 (RBD and sleep)

87. Boeve BF, Silber MH, Ferman TJ, et al. **Clinicopathologic correlations in 172 cases of rapid eye movement sleep behavior disorder with or without a coexisting neurologic disorder.** *Sleep Med* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Clinicopathologic+correlations+in+172+cases+of+rapid+eye+movement+sleep+behavior+disorder+with+or+without+a+coexisting+neurologic+disorder>
88. Postuma RB, Iranzo A, Hu M, et al. **Risk and predictors of dementia and parkinsonism in idiopathic REM sleep behaviour disorder: a multicentre study.** *Brain* 2019. — 연간 표현형 전환율 약 6%. 시나리오 19의 전구기 코호트 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Risk+and+predictors+of+dementia+and+parkinsonism+in+idiopathic+REM+sleep+behaviour+disorder+a+multicentre+study>
89. Luppi PH, Clément O, Sapin E, et al. **The neuronal network responsible for paradoxical sleep and its dysfunctions causing narcolepsy and rapid eye movement sleep behavior disorder.** *Sleep Med Rev* 2011. — SLD·GiV 글리신/GABA 회로. `Z_SLD`, `Z_GIV`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=The+neuronal+network+responsible+for+paradoxical+sleep+and+its+dysfunctions+causing+narcolepsy+and+rapid+eye+movement+sleep+behavior+disorder>
90. Kunz D, Mahlberg R. **A two-part, double-blind, placebo-controlled trial of exogenous melatonin in REM sleep behaviour disorder.** *J Sleep Res* 2010. — `MELON`/`MELEFF`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=A+two-part+double-blind+placebo-controlled+trial+of+exogenous+melatonin+in+REM+sleep+behaviour+disorder>
91. Ferman TJ, Smith GE, Dickson DW, et al. **Abnormal daytime sleepiness in dementia with Lewy bodies compared to Alzheimer's disease using the Multiple Sleep Latency Test.** *Alzheimers Res Ther* 2014. — 과도한 주간 졸림. `Z_EDS`, `EDSS`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Abnormal+daytime+sleepiness+in+dementia+with+Lewy+bodies+compared+to+Alzheimer+disease+using+the+Multiple+Sleep+Latency+Test>
92. Fronczek R, van Geest S, Frölich M, et al. **Hypocretin (orexin) loss in Alzheimer's disease.** *Neurobiol Aging* 2012. — 오렉신 뉴런 소실. `Z_ORX`.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hypocretin+orexin+loss+in+Alzheimer+disease+Fronczek>

## 14. 자율신경 기능과 MIBG (Autonomic features and MIBG)

93. Yoshita M, Arai H, Arai H, et al. **Diagnostic accuracy of ¹²³I-meta-iodobenzylguanidine myocardial scintigraphy in dementia with Lewy bodies: a multicenter study.** *PLoS One* 2015. — 민감도 ~69%, 특이도 ~87%. 모델 `MIBG = 1.05 + 1.75·CSYM` 및 `CSYM_0`(DLB 0.35 vs AD 0.85).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Diagnostic+accuracy+of+123I-meta-iodobenzylguanidine+myocardial+scintigraphy+in+dementia+with+Lewy+bodies+a+multicenter+study>
94. Orimo S, Uchihara T, Nakamura A, et al. **Axonal α-synuclein aggregates herald centripetal degeneration of cardiac sympathetic nerve in Parkinson's disease.** *Brain* 2008. — 심장 교감신경 종말의 α-시누클레인.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Axonal+alpha-synuclein+aggregates+herald+centripetal+degeneration+of+cardiac+sympathetic+nerve+in+Parkinson+disease>
95. Kaufmann H, Norcliffe-Kaufmann L, Palma JA, et al. **Natural history of pure autonomic failure: a United States prospective cohort.** *Ann Neurol* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Natural+history+of+pure+autonomic+failure+a+United+States+prospective+cohort>
96. Kaufmann H, Freeman R, Biaggioni I, et al. **Droxidopa for neurogenic orthostatic hypotension: a randomized, placebo-controlled, phase 3 trial.** *Neurology* 2014. — `DROXON`/`DROXEFF`, 그리고 앙와위 고혈압이라는 대가(`DROXSUP`).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Droxidopa+for+neurogenic+orthostatic+hypotension+a+randomized+placebo-controlled+phase+3+trial>
97. Robertson AD, Udow SJ, Espay AJ, et al. **Orthostatic hypotension and dementia incidence: links and implications.** *Neuropsychiatr Dis Treat* 2019. — OH가 인지 저하와 연결되는 기전(뇌관류). 모델 `WOH` 항.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Orthostatic+hypotension+and+dementia+incidence+links+and+implications>

## 15. AD 공존 병리 (Alzheimer co-pathology)

98. Irwin DJ, Lee VMY, Trojanowski JQ. **Parkinson's disease dementia: convergence of α-synuclein, tau and amyloid-β pathologies.** *Nat Rev Neurosci* 2013. — 세 병리의 상호작용. 모델 `B_CROSS`, `XSEED` 이선형(bilinear) 항.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Parkinson+disease+dementia+convergence+of+alpha-synuclein+tau+and+amyloid-beta+pathologies>
99. Masliah E, Rockenstein E, Veinbergs I, et al. **β-Amyloid peptides enhance α-synuclein accumulation and neuronal deficits in a transgenic mouse model linking Alzheimer's disease and Parkinson's disease.** *PNAS* 2001. — 교차 시딩의 직접 실험 근거.
    <https://pubmed.ncbi.nlm.nih.gov/?term=beta-Amyloid+peptides+enhance+alpha-synuclein+accumulation+and+neuronal+deficits+in+a+transgenic+mouse+model+linking+Alzheimer+disease+and+Parkinson+disease>
100. Ferman TJ, Aoki N, Crook JE, et al. **The limbic and neocortical contribution of α-synuclein, tau, and amyloid-β to disease duration and severity in dementia with Lewy bodies.** *Alzheimers Dement* 2018. — 공존 병리가 표현형과 병기간을 바꾼다. `B_PHENO`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+limbic+and+neocortical+contribution+of+alpha-synuclein+tau+and+amyloid-beta+to+disease+duration+and+severity+in+dementia+with+Lewy+bodies>
101. Lemstra AW, de Beer MH, Teunissen CE, et al. **Concomitant AD pathology affects clinical manifestation and survival in dementia with Lewy bodies.** *J Neurol Neurosurg Psychiatry* 2017. — 환시가 적고 생존이 짧아진다. 모델의 `B_PHENO → F_SURV`, `B_PHENO ⊣ V_VH` 점선 경로.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Concomitant+AD+pathology+affects+clinical+manifestation+and+survival+in+dementia+with+Lewy+bodies>

## 16. 신경염증 (Neuroinflammation)

102. Kim C, Ho DH, Suk JE, et al. **Neuron-released oligomeric α-synuclein is an endogenous agonist of TLR2 for paracrine activation of microglia.** *Nat Commun* 2013. — `I_MG`, `KMGON`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neuron-released+oligomeric+alpha-synuclein+is+an+endogenous+agonist+of+TLR2+for+paracrine+activation+of+microglia>
103. Liddelow SA, Guttenplan KA, Clarke LE, et al. **Neurotoxic reactive astrocytes are induced by activated microglia.** *Nature* 2017. — C1q+IL-1α+TNF에 의한 A1 성상교세포. `I_A1`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neurotoxic+reactive+astrocytes+are+induced+by+activated+microglia>
104. Sulzer D, Alcalay RN, Garretti F, et al. **T cells from patients with Parkinson's disease recognize α-synuclein peptides.** *Nature* 2017. — `I_TCELL`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=T+cells+from+patients+with+Parkinson+disease+recognize+alpha-synuclein+peptides>
105. Surendranathan A, Su L, Mak E, et al. **Early microglial activation and peripheral inflammation in dementia with Lewy bodies.** *Brain* 2018. — **초기에 활성화되고 후기에 소진된다** — 모델 `PHAGO = MGA·(1 − KEXH·MGEXH)`의 부호 전환.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Early+microglial+activation+and+peripheral+inflammation+in+dementia+with+Lewy+bodies>

## 17. 체액·영상 바이오마커 (Fluid and imaging biomarkers)

106. Fairfoul G, McGuire LI, Pal S, et al. **α-Synuclein RT-QuIC in the CSF of patients with alpha-synucleinopathies.** *Ann Clin Transl Neurol* 2016. — seed 증폭 분석. `W_SAA`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+RT-QuIC+in+the+CSF+of+patients+with+alpha-synucleinopathies>
107. Siderowf A, Concha-Marambio L, Lafontant DE, et al. **Assessment of heterogeneity among participants in the Parkinson's Progression Markers Initiative cohort using α-synuclein seed amplification: a cross-sectional study.** *Lancet Neurol* 2023. — 민감도·특이도 및 GBA/LRRK2 아형별 차이.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Assessment+of+heterogeneity+among+participants+in+the+Parkinson+Progression+Markers+Initiative+cohort+using+alpha-synuclein+seed+amplification>
108. Wang Z, Becker K, Donadio V, et al. **Skin α-synuclein aggregation seeding activity as a novel biomarker for Parkinson disease.** *JAMA Neurol* 2021.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Skin+alpha-synuclein+aggregation+seeding+activity+as+a+novel+biomarker+for+Parkinson+disease>
109. Pilotto A, Imarisio A, Carrarini C, et al. **Plasma neurofilament light chain predicts cognitive progression in prodromal and clinical dementia with Lewy bodies.** *J Alzheimers Dis* 2021. — `W_NFL`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Plasma+neurofilament+light+chain+predicts+cognitive+progression+in+prodromal+and+clinical+dementia+with+Lewy+bodies>
110. Palmqvist S, Janelidze S, Quiroz YT, et al. **Discriminative accuracy of plasma phospho-tau217 for Alzheimer disease vs other neurodegenerative disorders.** *JAMA* 2020. — `B_PTAU`, 공존 병리 층화.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Discriminative+accuracy+of+plasma+phospho-tau217+for+Alzheimer+disease+vs+other+neurodegenerative+disorders>

## 18. 콜린에스터라제 억제제 임상시험 (Cholinesterase inhibitor trials)

111. McKeith I, Del Ser T, Spano P, et al. **Efficacy of rivastigmine in dementia with Lewy bodies: a randomised, double-blind, placebo-controlled international study.** *Lancet* 2000. — **DLB 최초의 대규모 RCT.** 리바스티그민 6–12 mg/일, 20주에 NPI-4 약 30% 개선; 인지보다 정신행동·변동 지표에서 효과가 크다. 모델의 "ChEI는 평균보다 분산을 먼저 움직인다"는 예측의 보정 기준점.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Efficacy+of+rivastigmine+in+dementia+with+Lewy+bodies+a+randomised+double-blind+placebo-controlled+international+study>
112. Mori E, Ikeda M, Kosaka K; Donepezil-DLB Study Investigators. **Donepezil for dementia with Lewy bodies: a randomized, placebo-controlled trial.** *Ann Neurol* 2012. — 도네페질 10 mg에서 MMSE 약 +2.2, NPI-2 개선. 모델 시나리오 06의 보정 기준점.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Donepezil+for+dementia+with+Lewy+bodies+a+randomized+placebo-controlled+trial+Mori+Ikeda+Kosaka>
113. Ikeda M, Mori E, Matsuo K, et al. **Donepezil for dementia with Lewy bodies: a randomized, placebo-controlled, confirmatory phase III trial.** *Alzheimers Res Ther* 2015.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Donepezil+for+dementia+with+Lewy+bodies+a+randomized+placebo-controlled+confirmatory+phase+III+trial>
114. Stinton C, McKeith I, Taylor JP, et al. **Pharmacological management of Lewy body dementia: a systematic review and meta-analysis.** *Am J Psychiatry* 2015. — DLB/PDD 약물치료 종합. ChEI가 가장 확실한 근거를 갖는다는 결론.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacological+management+of+Lewy+body+dementia+a+systematic+review+and+meta-analysis>
115. Matsunaga S, Kishi T, Iwata N. **Combination therapy with cholinesterase inhibitors and memantine for Lewy body disorders: a systematic review and meta-analysis.** *Int J Geriatr Psychiatry* 2016. — 시나리오 07.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Combination+therapy+with+cholinesterase+inhibitors+and+memantine+for+Lewy+body+disorders+a+systematic+review+and+meta-analysis>
116. Winblad B, Cummings J, Andreasen N, et al. **A six-month double-blind, randomized, placebo-controlled study of a transdermal patch in Alzheimer's disease — rivastigmine patch versus capsule (IDEAL).** *Int J Geriatr Psychiatry* 2007. — **패치가 캡슐과 유사한 효능·낮은 위장관 부작용.** 모델의 `FPATCH`/`GIAE`(Cmax 구동) 구조가 재현하는 결과.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+six-month+double-blind+randomized+placebo-controlled+study+of+a+transdermal+patch+in+Alzheimer+disease+rivastigmine+patch+versus+capsule+IDEAL>
117. Gobburu JVS, Tammara V, Lesko L, et al. **Pharmacokinetic-pharmacodynamic modeling of rivastigmine, a cholinesterase inhibitor, in patients with Alzheimer's disease.** *J Clin Pharmacol* 2001. — 혈장 반감기(~1.5 h)와 **약력학 반감기(효소 재합성, ~8–10 h)의 분리.** 모델이 카바밀화 효소를 별도의 상태(`CARBA`/`CARBB`)로 둔 직접 근거.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic-pharmacodynamic+modeling+of+rivastigmine+a+cholinesterase+inhibitor+in+patients+with+Alzheimer+disease>
118. Bar-On P, Millard CB, Harel M, et al. **Kinetic and structural studies on the interaction of cholinesterases with the anti-Alzheimer drug rivastigmine.** *Biochemistry* 2002. — 유사비가역(pseudo-irreversible) 카바밀화 기전.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Kinetic+and+structural+studies+on+the+interaction+of+cholinesterases+with+the+anti-Alzheimer+drug+rivastigmine>
119. Rogers SL, Cooper NM, Sukovaty R, et al. **Pharmacokinetic and pharmacodynamic profile of donepezil HCl following multiple oral doses.** *Br J Clin Pharmacol* 1998. — 반감기 ~70 h, 정상상태 농도, 적혈구 AChE 억제율. 모델 `CLDON`/`VDON`/`IC50DON` 보정.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic+and+pharmacodynamic+profile+of+donepezil+HCl+following+multiple+oral+doses>
120. Kaasinen V, Någren K, Järvenpää T, et al. **Regional effects of donepezil and rivastigmine on cortical acetylcholinesterase activity in Alzheimer's disease.** *J Clin Psychopharmacol* 2002. — 생체 내 피질 AChE 억제율(약 30–40%대). 모델의 `INHACHE` 범위 검증.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Regional+effects+of+donepezil+and+rivastigmine+on+cortical+acetylcholinesterase+activity+in+Alzheimer+disease>

## 19. 메만틴 (Memantine)

121. Emre M, Tsolaki M, Bonuccelli U, et al. **Memantine for patients with Parkinson's disease dementia or dementia with Lewy bodies: a randomised, double-blind, placebo-controlled trial.** *Lancet Neurol* 2010. — 전반적 인상에서 소폭 이득. 시나리오 07.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+for+patients+with+Parkinson+disease+dementia+or+dementia+with+Lewy+bodies+a+randomised+double-blind+placebo-controlled+trial>
122. Aarsland D, Ballard C, Walker Z, et al. **Memantine in patients with Parkinson's disease dementia or dementia with Lewy bodies: a double-blind, placebo-controlled, multicentre trial.** *Lancet Neurol* 2009.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+in+patients+with+Parkinson+disease+dementia+or+dementia+with+Lewy+bodies+a+double-blind+placebo-controlled+multicentre+trial>
123. Parsons CG, Stöffler A, Danysz W. **Memantine: a NMDA receptor antagonist that improves memory by restoration of homeostasis in the glutamatergic system — too little activation is bad, too much is even worse.** *Neuropharmacology* 2007. — 빠른 off-rate 비경쟁 차단. `IC50MEM`, `EMAXMEM`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Memantine+a+NMDA+receptor+antagonist+that+improves+memory+by+restoration+of+homeostasis+in+the+glutamatergic+system>

## 20. 항정신병약과 피마반세린 (Antipsychotics and pimavanserin)

124. Cummings J, Isaacson S, Mills R, et al. **Pimavanserin for patients with Parkinson's disease psychosis: a randomised, placebo-controlled phase 3 trial.** *Lancet* 2014. — SAPS-PD 약 −3.06. 모델 시나리오 10의 보정 기준점(환시 부담 약 −37%).
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pimavanserin+for+patients+with+Parkinson+disease+psychosis+a+randomised+placebo-controlled+phase+3+trial>
125. Tariot PN, Cummings JL, Soto-Martin ME, et al. **Trial of pimavanserin in dementia-related psychosis (HARMONY).** *N Engl J Med* 2021. — 재발 예방 설계, DLB 하위군 포함. 시나리오 10.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+pimavanserin+in+dementia-related+psychosis+HARMONY+Tariot+Cummings>
126. Nasrallah HA, Fedora R, Morton R. **Successful treatment of clozapine-nonresponsive refractory hallucinations and delusions with pimavanserin, a serotonin 5HT-2A receptor inverse agonist.** *Schizophr Res* 2019.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Successful+treatment+of+clozapine-nonresponsive+refractory+hallucinations+and+delusions+with+pimavanserin>
127. Ballard C, Kreitzman DL, Isaacson S, et al. **Long-term evaluation of open-label pimavanserin safety and tolerability in Parkinson's disease psychosis.** *Parkinsonism Relat Disord* 2020. — QTc 연장 크기. `QTSLOPE`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Long-term+evaluation+of+open-label+pimavanserin+safety+and+tolerability+in+Parkinson+disease+psychosis>
128. Kyle K, Bronstein JM. **Treatment of psychosis in Parkinson's disease and dementia with Lewy bodies: a review.** *Parkinsonism Relat Disord* 2020.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+of+psychosis+in+Parkinson+disease+and+dementia+with+Lewy+bodies+a+review+Kyle+Bronstein>
129. Kurlan R, Cummings J, Raman R, Thal L; Alzheimer's Disease Cooperative Study Group. **Quetiapine for agitation or psychosis in patients with dementia and parkinsonism.** *Neurology* 2007. — 쿠에티아핀의 제한적 효능. 시나리오 11.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Quetiapine+for+agitation+or+psychosis+in+patients+with+dementia+and+parkinsonism+Kurlan+Cummings>
130. Schneider LS, Dagerman KS, Insel P. **Risk of death with atypical antipsychotic drug treatment for dementia: meta-analysis of randomized placebo-controlled trials.** *JAMA* 2005. — 사망 위험 상승(경고문의 근거). 모델의 `BNS` 곱셈적 위험 항.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Risk+of+death+with+atypical+antipsychotic+drug+treatment+for+dementia+meta-analysis+of+randomized+placebo-controlled+trials>
131. Ballard C, Hanney ML, Theodoulou M, et al. **The dementia antipsychotic withdrawal trial (DART-AD): long-term follow-up of a randomised placebo-controlled trial.** *Lancet Neurol* 2009.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+dementia+antipsychotic+withdrawal+trial+DART-AD+long-term+follow-up+of+a+randomised+placebo-controlled+trial>

## 21. 운동 증상 치료 (Motor symptomatic therapy)

132. Molloy S, McKeith IG, O'Brien JT, Burn DJ. **The role of levodopa in the management of dementia with Lewy bodies.** *J Neurol Neurosurg Psychiatry* 2005. — **DLB에서 레보도파 반응은 뚜렷하게 둔화되어 있다** (일반적으로 약 1/3에서만 의미 있는 반응). 시나리오 15/16 및 `responder_rates()`의 비교 대상.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+role+of+levodopa+in+the+management+of+dementia+with+Lewy+bodies>
133. Goetz CG, Blasucci LM, Leurgans S, Pappert EJ. **Olanzapine and clozapine: comparative effects on motor function in hallucinating PD patients.** *Neurology* 2000.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Olanzapine+and+clozapine+comparative+effects+on+motor+function+in+hallucinating+PD+patients>
134. Murata M, Odawara T, Hasegawa K, et al. **Adjunct zonisamide to levodopa for DLB parkinsonism: a randomized double-blind phase 2 study.** *Neurology* 2018. — 일본에서 DLB 파킨슨증에 승인된 조노사마이드. 시나리오 17.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Adjunct+zonisamide+to+levodopa+for+DLB+parkinsonism+a+randomized+double-blind+phase+2+study>
135. Murata M, Odawara T, Hasegawa K, et al. **Effect of zonisamide on parkinsonism in patients with dementia with Lewy bodies: a randomized, double-blind, placebo-controlled trial.** *Parkinsonism Relat Disord* 2020.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Effect+of+zonisamide+on+parkinsonism+in+patients+with+dementia+with+Lewy+bodies+a+randomized+double-blind+placebo-controlled+trial>
136. Nutt JG, Woodward WR, Hammerstad JP, et al. **The 'on-off' phenomenon in Parkinson's disease: relation to levodopa absorption and transport.** *N Engl J Med* 1984. — LNAA 경쟁과 위 배출. 모델 `U_ENS ⊣ X_LDC`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=The+on-off+phenomenon+in+Parkinson+disease+relation+to+levodopa+absorption+and+transport>
137. Contin M, Martinelli P. **Pharmacokinetics of levodopa.** *J Neurol* 2010. — 반감기 ~1.5 h(DDCI 병용), 생체이용률, 분포용적. `KALD`/`FLD`/`VLD`/`CLLD`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetics+of+levodopa+Contin+Martinelli>
138. Nagai M, Kubota M, Kimura Y, et al. **Anticholinergic burden and cognition in older adults with dementia** *(anticholinergic burden reviews)*. — 시나리오 09(`ANTICH`)의 근거.
     <https://pubmed.ncbi.nlm.nih.gov/?term=anticholinergic+burden+cognitive+decline+dementia+older+adults+systematic+review>

## 22. 질환 조절 후보 (Disease-modifying candidates)

139. Pagano G, Taylor KI, Anzures-Cabrera J, et al. **Trial of prasinezumab in early-stage Parkinson's disease (PASADENA).** *N Engl J Med* 2022. — 항 α-시누클레인 항체의 PK(4500 mg IV q4w)와 1차 평가변수 실패. 시나리오 21의 근거이자, 모델이 조절 효과를 **작게** 예측하도록 만든 이유.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+prasinezumab+in+early-stage+Parkinson+disease+PASADENA>
140. Lang AE, Siderowf AD, Macklin EA, et al. **Trial of cinpanemab in early Parkinson's disease (SPARK).** *N Engl J Med* 2022.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Trial+of+cinpanemab+in+early+Parkinson+disease+SPARK>
141. Cole TA, Zhao H, Collier TJ, et al. **α-Synuclein antisense oligonucleotides as a disease-modifying therapy for Parkinson's disease.** *JCI Insight* 2021. — `Y_ASO`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=alpha-Synuclein+antisense+oligonucleotides+as+a+disease-modifying+therapy+for+Parkinson+disease>
142. Jennings D, Huntwork-Rodriguez S, Vissers MFJM, et al. **LRRK2 inhibition by BIIB122 in healthy participants and patients with Parkinson's disease.** *Mov Disord* 2023. — `Y_LRRK`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=LRRK2+inhibition+by+BIIB122+in+healthy+participants+and+patients+with+Parkinson+disease>
143. van Dyck CH, Swanson CJ, Aisen P, et al. **Lecanemab in early Alzheimer's disease.** *N Engl J Med* 2023. — 공존 병리 팔에만 작용하는 치료(`Y_ADUC`)와 APOE ε4에서의 ARIA 위험.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Lecanemab+in+early+Alzheimer+disease+van+Dyck+Swanson+Aisen>

## 23. 예후와 생존 (Prognosis and survival)

144. Mueller C, Soysal P, Rongve A, et al. **Survival time and differences between dementia with Lewy bodies and Alzheimer's disease following diagnosis: a meta-analysis of longitudinal studies.** *Ageing Res Rev* 2019. — DLB의 생존이 AD보다 짧다. 모델 `H0`/`BCOG`/`BMOT`/`BFALL` 보정(중간 생존 약 4.4년).
     <https://pubmed.ncbi.nlm.nih.gov/?term=Survival+time+and+differences+between+dementia+with+Lewy+bodies+and+Alzheimer+disease+following+diagnosis+a+meta-analysis+of+longitudinal+studies>
145. Williams MM, Xiong C, Morris JC, Galvin JE. **Survival and mortality differences between dementia with Lewy bodies vs Alzheimer disease.** *Neurology* 2006.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Survival+and+mortality+differences+between+dementia+with+Lewy+bodies+vs+Alzheimer+disease>
146. Rongve A, Soennesyn H, Skogseth R, et al. **Cognitive decline in dementia with Lewy bodies: a 5-year prospective cohort study.** *BMJ Open* 2016. — MMSE 연간 저하 속도. 모델 자연경과 보정의 주요 기준점.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Cognitive+decline+in+dementia+with+Lewy+bodies+a+5-year+prospective+cohort+study+Rongve>
147. Kramberger MG, Auestad B, Garcia-Ptacek S, et al. **Long-term cognitive decline in dementia with Lewy bodies in a large multicenter, pooled cohort.** *J Alzheimers Dis* 2017.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Long-term+cognitive+decline+in+dementia+with+Lewy+bodies+in+a+large+multicenter+pooled+cohort>
148. Bostrom F, Jonsson L, Minthon L, Londos E. **Patients with dementia with Lewy bodies have more impaired quality of life than patients with Alzheimer disease.** *Alzheimer Dis Assoc Disord* 2007. — 돌봄 부담. `F_CARE`.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Patients+with+dementia+with+Lewy+bodies+have+more+impaired+quality+of+life+than+patients+with+Alzheimer+disease>

## 24. QSP 방법론 (QSP methodology)

149. Baron KT, Elmokadem A, et al. **mrgsolve: Simulate from ODE-Based Models.** R package. — 본 모델의 시뮬레이션 엔진.
     <https://mrgsolve.org/>
150. Gadkar K, Kirouac DC, Mager DE, et al. **A six-stage workflow for robust application of systems pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2016.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+six-stage+workflow+for+robust+application+of+systems+pharmacology>
151. Geerts H, Spiros A, Roberts P, Carr R. **Quantitative systems pharmacology as an extension of PK/PD modeling in CNS research and development.** *J Pharmacokinet Pharmacodyn* 2013. — 중추신경계 QSP에서 수용체 수준 전달자를 명시적으로 두는 접근. 본 모델의 세 전달자 구조와 직접 맞닿음.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative+systems+pharmacology+as+an+extension+of+PK+PD+modeling+in+CNS+research+and+development>
152. Geerts H, Roberts P, Spiros A, Carr R. **A strategy for developing new treatment paradigms for neuropsychiatric and neurocognitive symptoms in Alzheimer's disease.** *Front Pharmacol* 2013.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+strategy+for+developing+new+treatment+paradigms+for+neuropsychiatric+and+neurocognitive+symptoms+in+Alzheimer+disease>
153. Nicholas T, Duvvuri S, Leurent C, et al. **Systems pharmacology modeling in neuroscience: prediction and outcome of PF-04995274, a 5-HT4 partial agonist, in a clinical scopolamine impairment trial.** *Adv Alzheimers Dis* 2013. — 콜린성 전달자를 명시한 QSP 모델의 전향적 검증 사례.
     <https://pubmed.ncbi.nlm.nih.gov/?term=Systems+pharmacology+modeling+in+neuroscience+prediction+and+outcome+of+PF-04995274+a+5-HT4+partial+agonist+in+a+clinical+scopolamine+impairment+trial>
154. Bhattacharya S, Zhang Q, Andersen ME. **A deterministic map of Waddington's epigenetic landscape for cell fate decisions.** *BMC Syst Biol* 2011. — 안장-절점 분기와 이력현상(hysteresis)의 생물학적 모형화. 모델의 GBA1 스위치와 주의 상태 3차 함수.
     <https://pubmed.ncbi.nlm.nih.gov/?term=A+deterministic+map+of+Waddington+epigenetic+landscape+for+cell+fate+decisions>

---

## 문헌이 지지하지 못하는 부분 (What the literature does NOT support)

정직성을 위해 명시합니다. 아래 항목은 **모델의 가설**이며, 현재 근거가 없습니다.

1. **`FLUCTGT`의 구체적 함수 형태.** 인지 변동이 양안정 주의 상태의 *분산*이라는 것은
   문헌(#4, #74, #78, #79)과 정합적인 **해석**이지만, `SIGN0·(FLBASE + BISTAB·NOISEG)`라는
   함수형은 어떤 데이터에도 적합(fit)되지 않았습니다. 이것이 만드는 검증 가능한 예측은
   "ChEI가 MMSE보다 변동 지표를 비례적으로 더 크게 개선한다"이고, 그 비율은 모델이
   약 2.5배로 계산합니다 — 이 비율 자체가 반증 가능한 주장입니다.
2. **`KSUPP`(변연계 α-시누클레인이 후시냅스 선조체 통합성을 억제하는 강도).** #56이
   "DLB에서 D2가 상향조절되지 않는다"는 관찰을 제공하지만, 그 관계의 **기울기**는
   측정된 바 없습니다. 모델은 신경이완제 민감성 발생률(#54, #55)에 맞춰 이 하나의
   파라미터를 역산했습니다.
3. **앰브록솔의 "기한(deadline)".** 시나리오 19 대 20의 대비는 GBA1 되먹임 고리가
   진짜로 안장-절점 분기를 갖는다는 전제에 전적으로 의존합니다. #31과 #39는 고리의
   **존재**를 지지하지만 **양안정성**을 증명하지는 않습니다. 이것이 이 모델의 가장 큰
   외삽이며, 동시에 가장 검증하기 쉬운 예측입니다(전구기 대 확립기 코호트에서 동일
   노출을 비교하면 됩니다).
4. **`MOTNDA`(비도파민성 축성 운동 부담).** DLB 파킨슨증의 상당 부분이 레보도파에
   반응하지 않는다는 것(#132)은 잘 알려져 있으나, 그 비율을 신피질·변연계 섬유
   부하의 선형 함수로 둔 것은 모델의 선택입니다.
5. **모든 신경 소실 속도상수(`KL*`).** 이들은 개별적으로 측정된 값이 아니라, 자연경과
   보정 기준점(MMSE 연간 저하, MDS-UPDRS III 연간 상승, MIBG, DaTSCAN, 중간 생존)을
   동시에 맞추도록 조정된 값입니다. 서로 상관되어 있으므로 개별 값에 생물학적 의미를
   부여하면 안 됩니다.

---

*본 문헌 목록은 교육·연구 목적입니다. 임상 의사결정에 사용하지 마십시오.*
*This reference list is for education and research only. Not for clinical use.*
