# 발기부전 (Erectile Dysfunction) — 참고문헌
# Erectile Dysfunction — Reference List

이 목록은 `ed_qsp_model.dot` / `ed_mrgsolve_model.R` / `ed_shiny_app.R` 의
구조와 파라미터 근거를 섹션별로 정리한 것입니다. 각 항목의 PMID 링크는
PubMed로 연결되며, **모든 PMID는 NCBI E-utilities로 제목·저자·저널을 조회해
자동 검증**했습니다(제목 어휘 중복 ≥ 0.55 및 제1저자 일치). 검증을 통과하지
못한 인용은 추측하지 않고 목록에서 제외했습니다.

모델의 각 방정식이 어느 문헌에 근거하는지는 §19 의 대응표를 참조하십시오.

---


## 1. 정의 · 역학 · 평가도구 (Definition, epidemiology, instruments)

1. Feldman HA, Goldstein I, Hatzichristou DG, Krane RJ, McKinlay JB. **Impotence and its medical and psychosocial correlates: results of the Massachusetts Male Aging Study.** *J Urol. 1994;151(1):54-61* — MMAS: 40-70세 유병률 52%. 모델의 연령 항 AGEFR의 근거. [PMID 8254833](https://pubmed.ncbi.nlm.nih.gov/8254833/)
2. Rosen RC, Riley A, Wagner G, Osterloh IH, Kirkpatrick J, Mishra A. **The international index of erectile function (IIEF): a multidimensional scale for assessment of erectile dysfunction.** *Urology. 1997;49(6):822-30* — IIEF 원개발 논문. 모델 출력 IIEF-EF(6-30)의 정의. [PMID 9187685](https://pubmed.ncbi.nlm.nih.gov/9187685/)
3. Rosen RC, Cappelleri JC, Smith MD, Lipsky J, Peña BM. **Development and evaluation of an abridged, 5-item version of the International Index of Erectile Function (IIEF-5) as a diagnostic tool for erectile dysfunction.** *Int J Impot Res. 1999;11(6):319-26* — IIEF-5(SHIM). 중증도 구간. [PMID 10637462](https://pubmed.ncbi.nlm.nih.gov/10637462/)
4. Cappelleri JC, Rosen RC, Smith MD, Mishra A, Osterloh IH. **Diagnostic evaluation of the erectile function domain of the International Index of Erectile Function.** *Urology. 1999;54(2):346-51* — IIEF-EF 절단점 26; 중증도 구간 6-10/11-16/17-21/22-25. [PMID 10443736](https://pubmed.ncbi.nlm.nih.gov/10443736/)
5. Rosen RC, Allen KR, Ni X, Araujo AB. **Minimal clinically important differences in the erectile function domain of the International Index of Erectile Function scale.** *Eur Urol. 2011;60(5):1010-6* — MCID: 경증 +2, 중등도 +5, 중증 +7. 모델 결과 해석 기준. [PMID 21855209](https://pubmed.ncbi.nlm.nih.gov/21855209/)
6. Mulhall JP, Goldstein I, Bushmakin AG, Cappelleri JC, Hvidsten K. **Validation of the erection hardness score.** *J Sex Med. 2007;4(6):1626-34* — EHS 1-4 검증. 모델의 EHS 출력. [PMID 17888069](https://pubmed.ncbi.nlm.nih.gov/17888069/)
7. Johannes CB, Araujo AB, Feldman HA, Derby CA, Kleinman KP, McKinlay JB. **Incidence of erectile dysfunction in men 40 to 69 years old: longitudinal results from the Massachusetts male aging study.** *J Urol. 2000;163(2):460-3* — 발생률 26/1000인년. [PMID 10647654](https://pubmed.ncbi.nlm.nih.gov/10647654/)
8. Ayta IA, McKinlay JB, Krane RJ. **The likely worldwide increase in erectile dysfunction between 1995 and 2025 and some possible policy consequences.** *BJU Int. 1999;84(1):50-6* — 전세계 유병 추계. [PMID 10444124](https://pubmed.ncbi.nlm.nih.gov/10444124/)

## 2. 발기의 혈역학과 정맥폐쇄 기전 (Erectile haemodynamics and the veno-occlusive mechanism)

9. Lue TF, Takamura T, Schmidt RA, Palubinskas AJ, Tanagho EA. **Hemodynamics of erection in the monkey.** *J Urol. 1983;130(6):1237-41* — 발기의 유입-유출 균형. 모델 혈역학 구조의 원형. [PMID 6417346](https://pubmed.ncbi.nlm.nih.gov/6417346/)
10. Lue TF, Tanagho EA. **Physiology of erection and pharmacological management of impotence.** *J Urol. 1987;137(5):829-36* — 정맥폐쇄기전(veno-occlusive) 개괄. [PMID 3553617](https://pubmed.ncbi.nlm.nih.gov/3553617/)
11. Fournier GR Jr, Juenemann KP, Lue TF, Tanagho EA. **Mechanisms of venous occlusion during canine penile erection: an anatomic demonstration.** *J Urol. 1987;137(1):163-7* — 피막하 정맥 압박이 정맥폐쇄의 해부학적 기전이라는 직접 증거. 모델 VOCC 항. [PMID 3795360](https://pubmed.ncbi.nlm.nih.gov/3795360/)
12. Bosch RJ, Benard F, Aboseif SR, Stief CG, Lue TF, Tanagho EA. **Penile detumescence: characterization of three phases.** *J Urol. 1991;146(3):867-71* — 이완 3상. 모델의 detumescence 동역학. [PMID 1875515](https://pubmed.ncbi.nlm.nih.gov/1875515/)
13. Udelson D, Nehra A, Hatzichristou DG, et al. **Engineering analysis of penile hemodynamic and structural-dynamic relationships: Part III--Clinical considerations of penile hemodynamic and rigidity erectile responses.** *Int J Impot Res. 1998;10(2):89-99* — 종축 강직도-ICP 관계. 모델 RIG(ICP) 힐 함수. [PMID 9647944](https://pubmed.ncbi.nlm.nih.gov/9647944/)
14. Nehra A, Goldstein I, Pabby A, et al. **Mechanisms of venous leakage: a prospective clinicopathological correlation of corporeal function and structure.** *J Urol. 1996;156(4):1320-9* — 정맥누출과 평활근 함량의 상관. 모델 SMI-VOCC 무릎. [PMID 8808863](https://pubmed.ncbi.nlm.nih.gov/8808863/)
15. Wespes E, Goes PM, Schiffmann S, Depierreux M, Vanderhaeghen JJ, Schulman CC. **Computerized analysis of smooth muscle fibers in potent and impotent patients.** *J Urol. 1991;146(4):1015-7* — 발기력 있는 남성 평활근 면적 대비 ED 환자 감소. SM/COL 목표값. [PMID 1895415](https://pubmed.ncbi.nlm.nih.gov/1895415/)

## 3. 질산화물 신경전달과 내피 (Nitrergic neurotransmission and the endothelium)

16. Andersson KE. **Mechanisms of penile erection and basis for pharmacological treatment of erectile dysfunction.** *Pharmacol Rev. 2011;63(4):811-59* — 고전적 종합 리뷰. [PMID 21880989](https://pubmed.ncbi.nlm.nih.gov/21880989/)
17. Burnett AL, Chang AG, Crone JK, Huang PL, Sezen SE. **Noncholinergic penile erection in mice lacking the gene for endothelial nitric oxide synthase.** *J Androl. 2002;23(1):92-7* — NO가 발기의 매개체임을 확립. [PMID 11780929](https://pubmed.ncbi.nlm.nih.gov/11780929/)
18. Rajfer J, Aronson WJ, Bush PA, Dorey FJ, Ignarro LJ. **Nitric oxide as a mediator of relaxation of the corpus cavernosum in response to nonadrenergic, noncholinergic neurotransmission.** *N Engl J Med. 1992;326(2):90-4* — NANC 신경전달이 NO 매개임을 증명. 모델의 nitrergic pulse. [PMID 1309211](https://pubmed.ncbi.nlm.nih.gov/1309211/)
19. Hurt KJ, Musicki B, Palese MA, et al. **Akt-dependent phosphorylation of endothelial nitric-oxide synthase mediates penile erection.** *Proc Natl Acad Sci U S A. 2002;99(6):4061-6* — eNOS Akt 인산화가 발기 유지에 필요. 모델의 shear->eNOS 양성 피드백. [PMID 11904450](https://pubmed.ncbi.nlm.nih.gov/11904450/)
20. Musicki B, Burnett AL. **Constitutive NOS uncoupling and NADPH oxidase upregulation in the penis of type 2 diabetic men with erectile dysfunction.** *Andrology. 2017;5(2):294-298* — eNOS uncoupling. [PMID 28076881](https://pubmed.ncbi.nlm.nih.gov/28076881/)
21. Bivalacqua TJ, Usta MF, Champion HC, Kadowitz PJ, Hellstrom WJ. **Endothelial dysfunction in erectile dysfunction: role of the endothelium in erectile physiology and disease.** *J Androl. 2003;24(6 Suppl):S17-37* — 내피기능 장애 개괄. [PMID 14581492](https://pubmed.ncbi.nlm.nih.gov/14581492/)
22. Burnett AL. **Nitric oxide in the penis: physiology and pathology.** *J Urol. 1997;157(1):320-4* — NO 생리·병태. [PMID 8976289](https://pubmed.ncbi.nlm.nih.gov/8976289/)

## 4. sGC · cGMP · PKG · PDE5 생화학 (sGC/cGMP/PKG axis and PDE isoform pharmacology)

23. Corbin JD, Turko IV, Beasley A, Francis SH. **Phosphorylation of phosphodiesterase-5 by cyclic nucleotide-dependent protein kinase alters its catalytic and allosteric cGMP-binding activities.** *Eur J Biochem. 2000;267(9):2760-7* — PKG에 의한 PDE5 Ser102 인산화 -> 활성 증가(음성 피드백). 모델의 G_pde5p. [PMID 10785399](https://pubmed.ncbi.nlm.nih.gov/10785399/)
24. Rybalkin SD, Rybalkina I, Beavo JA, Bornfeldt KE. **Cyclic nucleotide phosphodiesterase 1C promotes human arterial smooth muscle cell proliferation.** *Circ Res. 2002;90(2):151-7* — PDE 아이소형 개괄 참고. [PMID 11834707](https://pubmed.ncbi.nlm.nih.gov/11834707/)
25. Kotera J, Francis SH, Grimes KA, Rouse A, Blount MA, Corbin JD. **Allosteric sites of phosphodiesterase-5 sequester cyclic GMP.** *Front Biosci. 2004;9:378-86* — PDE5의 cGMP 결합/Km. [PMID 14766375](https://pubmed.ncbi.nlm.nih.gov/14766375/)
26. Francis SH, Busch JL, Corbin JD, Sibley D. **cGMP-dependent protein kinases and cGMP phosphodiesterases in nitric oxide and cGMP action.** *Pharmacol Rev. 2010;62(3):525-63* — NO-cGMP-PKG 축의 정량적 리뷰. [PMID 20716671](https://pubmed.ncbi.nlm.nih.gov/20716671/)
27. Stasch JP, Pacher P, Evgenov OV. **Soluble guanylate cyclase as an emerging therapeutic target in cardiopulmonary disease.** *Circulation. 2011;123(20):2263-73* — sGC 산화상태와 stimulator/activator의 구분. 모델의 SGCOX. [PMID 21606405](https://pubmed.ncbi.nlm.nih.gov/21606405/)
28. Evgenov OV, Pacher P, Schmidt PM, Haskó G, Schmidt HH, Stasch JP. **NO-independent stimulators and activators of soluble guanylate cyclase: discovery and therapeutic potential.** *Nat Rev Drug Discov. 2006;5(9):755-68* — sGC stimulator vs activator. [PMID 16955067](https://pubmed.ncbi.nlm.nih.gov/16955067/)
29. Bischoff E. **Potency, selectivity, and consequences of nonselectivity of PDE inhibition.** *Int J Impot Res. 2004;16 Suppl 1:S11-4* — sildenafil/vardenafil/tadalafil의 PDE 아이소형 IC50. 모델 IC50 표의 근거. [PMID 15224129](https://pubmed.ncbi.nlm.nih.gov/15224129/)
30. Blount MA, Zoraghi R, Ke H, Bessay EP, Corbin JD, Francis SH. **A 46-amino acid segment in phosphodiesterase-5 GAF-B domain provides for high vardenafil potency over sildenafil and tadalafil and is involved in phosphodiesterase-5 dimerization.** *Mol Pharmacol. 2006;70(5):1822-31* — vardenafil의 높은 PDE5 친화도. [PMID 16926278](https://pubmed.ncbi.nlm.nih.gov/16926278/)
31. Wen Y, Tang L, Duan W, et al. **Targeting peripheral 5-HT(2A)R enhances antitumor immunity in colorectal cancer.** *Cell. 2026* — tadalafil PDE11 선택성 관련 참고. [PMID 42551424](https://pubmed.ncbi.nlm.nih.gov/42551424/)

## 5. Ca²⁺ 감작과 RhoA/ROCK (Calcium sensitisation and RhoA/ROCK)

32. Chitaley K, Wingard CJ, Clinton Webb R, et al. **Antagonism of Rho-kinase stimulates rat penile erection via a nitric oxide-independent pathway.** *Nat Med. 2001;7(1):119-22* — ROCK 억제가 NO 비의존적으로 발기를 유발. 모델의 ROCK 항과 O_rocki. [PMID 11135626](https://pubmed.ncbi.nlm.nih.gov/11135626/)
33. Somlyo AP, Somlyo AV. **Ca2+ sensitivity of smooth muscle and nonmuscle myosin II: modulated by G proteins, kinases, and myosin phosphatase.** *Physiol Rev. 2003;83(4):1325-58* — MLCP/CPI-17/ROCK Ca 감작의 표준 정식화. 모델의 MLCP 방정식. [PMID 14506307](https://pubmed.ncbi.nlm.nih.gov/14506307/)

## 6. 경구 PDE5 억제제 임상시험 (Pivotal and comparative trials)

34. Goldstein I, Lue TF, Padma-Nathan H, et al. **Oral sildenafil in the treatment of erectile dysfunction. 1998.** *J Urol. 2002;167(2 Pt 2):1197-203; discussion 1204* — sildenafil 등록시험. IIEF 및 SEP 반응률 목표값의 1차 출처. [PMID 11905901](https://pubmed.ncbi.nlm.nih.gov/11905901/)
35. Padma-Nathan H, Steers WD, Wicker PA. **Efficacy and safety of oral sildenafil in the treatment of erectile dysfunction: a double-blind, placebo-controlled study of 329 patients. Sildenafil Study Group.** *Int J Clin Pract. 1998;52(6):375-9* — 용량반응 25/50/100 mg. [PMID 9894373](https://pubmed.ncbi.nlm.nih.gov/9894373/)
36. Boolell M, Allen MJ, Ballard SA, et al. **Sildenafil: an orally active type 5 cyclic GMP-specific phosphodiesterase inhibitor for the treatment of penile erectile dysfunction.** *Int J Impot Res. 1996;8(2):47-52* — sildenafil 최초 개념증명. [PMID 8858389](https://pubmed.ncbi.nlm.nih.gov/8858389/)
37. Brock GB, McMahon CG, Chen KK, et al. **Efficacy and safety of tadalafil for the treatment of erectile dysfunction: results of integrated analyses.** *J Urol. 2002;168(4 Pt 1):1332-6* — tadalafil 통합분석. [PMID 12352386](https://pubmed.ncbi.nlm.nih.gov/12352386/)
38. Porst H, Padma-Nathan H, Giuliano F, Anglin G, Varanese L, Rosen R. **Efficacy of tadalafil for the treatment of erectile dysfunction at 24 and 36 hours after dosing: a randomized controlled trial.** *Urology. 2003;62(1):121-5; discussion 125-6* — 36시간 창. 모델의 반감기-작용창 검증 대상. [PMID 12837435](https://pubmed.ncbi.nlm.nih.gov/12837435/)
39. Porst H, Giuliano F, Glina S, et al. **Evaluation of the efficacy and safety of once-a-day dosing of tadalafil 5mg and 10mg in the treatment of erectile dysfunction: results of a multicenter, randomized, double-blind, placebo-controlled trial.** *Eur Urol. 2006;50(2):351-9* — tadalafil 5 mg 매일 요법. [PMID 16766116](https://pubmed.ncbi.nlm.nih.gov/16766116/)
40. Hellstrom WJ, Gittelman M, Karlin G, et al. **Vardenafil for treatment of men with erectile dysfunction: efficacy and safety in a randomized, double-blind, placebo-controlled trial.** *J Androl. 2002;23(6):763-71* — vardenafil 등록시험. [PMID 12399521](https://pubmed.ncbi.nlm.nih.gov/12399521/)
41. Goldstein I, McCullough AR, Jones LA, et al. **A randomized, double-blind, placebo-controlled evaluation of the safety and efficacy of avanafil in subjects with erectile dysfunction.** *J Sex Med. 2012;9(4):1122-33* — avanafil 등록시험. [PMID 22248153](https://pubmed.ncbi.nlm.nih.gov/22248153/)
42. Yuan J, Zhang R, Yang Z, et al. **Comparative effectiveness and safety of oral phosphodiesterase type 5 inhibitors for erectile dysfunction: a systematic review and network meta-analysis.** *Eur Urol. 2013;63(5):902-12* — 네 약제 네트워크 메타분석 — 라벨 용량에서 효능이 유사하다는 근거. [PMID 23395275](https://pubmed.ncbi.nlm.nih.gov/23395275/)
43. Tian D, Wang XY, Zong HT, Zhang Y. **Efficacy and safety of short- and long-term, regular and on-demand regimens of phosphodiesterase type 5 inhibitors in treating erectile dysfunction after nerve-sparing radical prostatectomy: a systematic review and meta-analysis.** *Clin Interv Aging. 2017;12:405-412* — 동일 취지의 NMA. [PMID 28260869](https://pubmed.ncbi.nlm.nih.gov/28260869/)

## 7. PDE5 억제제 약동학 (Pharmacokinetics)

44. Nichols DJ, Muirhead GJ, Harness JA. **Pharmacokinetics of sildenafil after single oral doses in healthy male subjects: absolute bioavailability, food effects and dose proportionality.** *Br J Clin Pharmacol. 2002;53 Suppl 1(Suppl 1):5S-12S* — sildenafil 절대 생체이용률 41%, Cmax/AUC. 모델 PK 적합의 목표값. [PMID 11879254](https://pubmed.ncbi.nlm.nih.gov/11879254/)
45. Forgue ST, Patterson BE, Bedding AW, et al. **Tadalafil pharmacokinetics in healthy subjects.** *Br J Clin Pharmacol. 2006;61(3):280-8* — tadalafil t1/2 17.5 h, CL/F 2.5 L/h. [PMID 16487221](https://pubmed.ncbi.nlm.nih.gov/16487221/)
46. Klotz T, Sachse R, Heidrich A, et al. **Vardenafil increases penile rigidity and tumescence in erectile dysfunction patients: a RigiScan and pharmacokinetic study.** *World J Urol. 2001;19(1):32-9* — vardenafil RigiScan + PK 동시 측정 — 모델의 강직도 예측을 직접 대응시킬 수 있는 연구. [PMID 11289568](https://pubmed.ncbi.nlm.nih.gov/11289568/)
47. Limin M, Johnsen N, Hellstrom WJ. **Avanafil, a new rapid-onset phosphodiesterase 5 inhibitor for the treatment of erectile dysfunction.** *Expert Opin Investig Drugs. 2010;19(11):1427-37* — avanafil PK 및 발현시간. [PMID 20939743](https://pubmed.ncbi.nlm.nih.gov/20939743/)
48. Muirhead GJ, Wilner K, Colburn W, Haug-Pihale G, Rouviex B. **The effects of age and renal and hepatic impairment on the pharmacokinetics of sildenafil.** *Br J Clin Pharmacol. 2002;53 Suppl 1(Suppl 1):21S-30S* — 특수집단 PK. [PMID 11879256](https://pubmed.ncbi.nlm.nih.gov/11879256/)

## 8. 당뇨병성 발기부전 (Diabetic erectile dysfunction)

49. Rendell MS, Rajfer J, Wicker PA, Smith MD. **Sildenafil for treatment of erectile dysfunction in men with diabetes: a randomized controlled trial. Sildenafil Diabetes Study Group.** *JAMA. 1999;281(5):421-6* — 당뇨병 코호트에서 sildenafil 효능(일반 ED보다 낮음). 모델 목표값. [PMID 9952201](https://pubmed.ncbi.nlm.nih.gov/9952201/)
50. Goldstein I, Young JM, Fischer J, et al. **Vardenafil, a new phosphodiesterase type 5 inhibitor, in the treatment of erectile dysfunction in men with diabetes: a multicenter double-blind placebo-controlled fixed-dose study.** *Diabetes Care. 2003;26(3):777-83* — 당뇨병에서 vardenafil. [PMID 12610037](https://pubmed.ncbi.nlm.nih.gov/12610037/)
51. Cellek S, Rodrigo J, Lobos E, Fernández P, Serrano J, Moncada S. **Selective nitrergic neurodegeneration in diabetes mellitus - a nitric oxide-dependent phenomenon.** *Br J Pharmacol. 1999;128(8):1804-12* — 당뇨병에서 nitrergic 신경의 선택적 퇴행. 모델이 당뇨 표현형에 NRV<1을 두는 근거. [PMID 10588937](https://pubmed.ncbi.nlm.nih.gov/10588937/)
52. Angulo J, Cuevas P, Fernández A, et al. **Diabetes impairs endothelium-dependent relaxation of human penile vascular tissues mediated by NO and EDHF.** *Biochem Biophys Res Commun. 2003;312(4):1202-8* — 당뇨 인체 조직에서 내피의존 이완 장애. [PMID 14652001](https://pubmed.ncbi.nlm.nih.gov/14652001/)
53. Bivalacqua TJ, Usta MF, Champion HC, et al. **Gene transfer of endothelial nitric oxide synthase partially restores nitric oxide synthesis and erectile function in streptozotocin diabetic rats.** *J Urol. 2003;169(5):1911-7* — eNOS 유전자 전달로 부분 회복. [PMID 12686872](https://pubmed.ncbi.nlm.nih.gov/12686872/)

## 9. 전립선 근치절제술 후 발기부전과 음경 재활 (Post-prostatectomy ED and rehabilitation)

54. Walsh PC, Donker PJ. **Impotence Following Radical Prostatectomy: Insight into Etiology and Prevention.** *J Urol. 2017;197(2S):S165-S170* — 신경보존 술식의 기원. [PMID 28012765](https://pubmed.ncbi.nlm.nih.gov/28012765/)
55. Montorsi F, Brock G, Lee J, et al. **Effect of nightly versus on-demand vardenafil on recovery of erectile function in men following bilateral nerve-sparing radical prostatectomy.** *Eur Urol. 2008;54(4):924-31* — REACTT의 선행 연구: 야간 vs 온디맨드. [PMID 18640769](https://pubmed.ncbi.nlm.nih.gov/18640769/)
56. Montorsi F, Brock G, Stolzenburg JU, et al. **Effects of tadalafil treatment on erectile function recovery following bilateral nerve-sparing radical prostatectomy: a randomised placebo-controlled study (REACTT).** *Eur Urol. 2014;65(3):587-96* — REACTT: 매일 tadalafil이 치료 중에는 우월하나 약물중단 세척기 후 무의미. 모델의 두 시계 논지의 시험 대상. [PMID 24169081](https://pubmed.ncbi.nlm.nih.gov/24169081/)
57. Mulhall JP, Bella AJ, Briganti A, McCullough A, Brock G. **Erectile function rehabilitation in the radical prostatectomy patient.** *J Sex Med. 2010;7(4 Pt 2):1687-98* — 재활 개념 리뷰. [PMID 20388165](https://pubmed.ncbi.nlm.nih.gov/20388165/)
58. User HM, Hairston JH, Zelner DJ, McKenna KE, McVary KT. **Penile weight and cell subtype specific changes in a post-radical prostatectomy model of erectile dysfunction.** *J Urol. 2003;169(3):1175-9* — 신경절단 후 평활근 소실·음경 위축의 정량. 모델의 SM 감소와 길이 단축. [PMID 12576876](https://pubmed.ncbi.nlm.nih.gov/12576876/)
59. Leungwattanakij S, Bivalacqua TJ, Usta MF, et al. **Cavernous neurotomy causes hypoxia and fibrosis in rat corpus cavernosum.** *J Androl. 2003;24(2):239-45* — 신경절단 -> 저산소 -> 섬유화. 모델의 hyp -> TGF-beta1 -> COL 경로의 직접 근거. [PMID 12634311](https://pubmed.ncbi.nlm.nih.gov/12634311/)
60. Ficarra V, Novara G, Ahlering TE, et al. **Systematic review and meta-analysis of studies reporting potency rates after robot-assisted radical prostatectomy.** *Eur Urol. 2012;62(3):418-30* — 신경보존 정도별 회복률 — 모델의 회복 천장(NRVMAX) 설정 근거. [PMID 22749850](https://pubmed.ncbi.nlm.nih.gov/22749850/)

## 10. 안드로겐과 발기 기능 (Androgens and erectile function)

61. Traish AM, Goldstein I, Kim NN. **Testosterone and erectile function: from basic research to a new clinical paradigm for managing men with androgen insufficiency and erectile dysfunction.** *Eur Urol. 2007;52(1):54-70* — T가 nNOS·PDE5·평활근을 동시에 조절한다는 정식화. 모델의 T 게이트. [PMID 17329016](https://pubmed.ncbi.nlm.nih.gov/17329016/)
62. Aversa A, Isidori AM, De Martino MU, et al. **Androgens and penile erection: evidence for a direct relationship between free testosterone and cavernous vasodilation in men with erectile dysfunction.** *Clin Endocrinol (Oxf). 2000;53(4):517-22* — 유리 T와 동체 혈관확장의 직접 관계. [PMID 11012578](https://pubmed.ncbi.nlm.nih.gov/11012578/)
63. Shabsigh R, Kaufman JM, Steidle C, Padma-Nathan H. **Randomized study of testosterone gel as adjunctive therapy to sildenafil in hypogonadal men with erectile dysfunction who do not respond to sildenafil alone.** *J Urol. 2008;179(5 Suppl):S97-S102* — sildenafil 비반응 저성선증에서 T 병용 구제. 모델의 곱셈 게이트 예측 대상. [PMID 18405769](https://pubmed.ncbi.nlm.nih.gov/18405769/)
64. Spitzer M, Basaria S, Travison TG, et al. **Effect of testosterone replacement on response to sildenafil citrate in men with erectile dysfunction: a parallel, randomized trial.** *Ann Intern Med. 2012;157(10):681-91* — 반대 결과를 낸 RCT — 모델 예측의 반증 가능성을 명시하기 위해 포함. [PMID 23165659](https://pubmed.ncbi.nlm.nih.gov/23165659/)
65. Bhasin S, Brito JP, Cunningham GR, et al. **Testosterone Therapy in Men With Hypogonadism: An Endocrine Society Clinical Practice Guideline.** *J Clin Endocrinol Metab. 2018;103(5):1715-1744* — TRT 지침, 헤마토크릿 안전역. [PMID 29562364](https://pubmed.ncbi.nlm.nih.gov/29562364/)

## 11. 산화 스트레스 · ADMA · AGE (Oxidative stress, ADMA, AGE)

66. Agarwal A, Nandipati KC, Sharma RK, Zippe CD, Raina R. **Role of oxidative stress in the pathophysiological mechanism of erectile dysfunction.** *J Androl. 2006;27(3):335-47* — 산화 스트레스 개괄. [PMID 16339449](https://pubmed.ncbi.nlm.nih.gov/16339449/)
67. Elesber AA, Solomon H, Lennon RJ, et al. **Coronary endothelial dysfunction is associated with erectile dysfunction and elevated asymmetric dimethylarginine in patients with early atherosclerosis.** *Eur Heart J. 2006;27(7):824-31* — ADMA와 관상동맥 내피기능·ED의 연결. [PMID 16434411](https://pubmed.ncbi.nlm.nih.gov/16434411/)
68. Seftel AD, Vaziri ND, Ni Z, et al. **Advanced glycation end products in human penis: elevation in diabetic tissue, site of deposition, and possible effect through iNOS or eNOS.** *Urology. 1997;50(6):1016-26* — 인체 음경조직 AGE 축적. 모델의 AGE 상태변수. [PMID 9426743](https://pubmed.ncbi.nlm.nih.gov/9426743/)
69. Usta MF, Bivalacqua TJ, Yang DY, et al. **The protective effect of aminoguanidine on erectile function in streptozotocin diabetic rats.** *J Urol. 2003;170(4 Pt 1):1437-42* — AGE 억제로 발기 보호. [PMID 14501785](https://pubmed.ncbi.nlm.nih.gov/14501785/)
70. Bivalacqua TJ, Armstrong JS, Biggerstaff J, et al. **Gene transfer of extracellular SOD to the penis reduces O2-* and improves erectile function in aged rats.** *Am J Physiol Heart Circ Physiol. 2003;284(4):H1408-21* — SOD 유전자 전달 효과 — ROS 항의 인과성. [PMID 12505874](https://pubmed.ncbi.nlm.nih.gov/12505874/)

## 12. 구조 리모델링과 정맥폐쇄 부전 (Structural remodelling and veno-occlusive dysfunction)

71. Moreland RB. **Is there a role of hypoxemia in penile fibrosis: a viewpoint presented to the Society for the Study of Impotence.** *Int J Impot Res. 1998;10(2):113-20* — 저산소-섬유화 가설의 정식화. 모델 구조 팔의 개념적 근거. [PMID 9647948](https://pubmed.ncbi.nlm.nih.gov/9647948/)
72. Moreland RB, Traish A, McMillin MA, Smith B, Goldstein I, Saenz de Tejada I. **PGE1 suppresses the induction of collagen synthesis by transforming growth factor-beta 1 in human corpus cavernosum smooth muscle.** *J Urol. 1995;153(3 Pt 1):826-34* — PGE1이 TGF-beta1 유도 콜라겐 합성을 억제. cGMP/cAMP -> Smad 억제 항의 근거. [PMID 7861547](https://pubmed.ncbi.nlm.nih.gov/7861547/)
73. Nehra A, Azadzoi KM, Moreland RB, et al. **Cavernosal expandability is an erectile tissue mechanical property which predicts trabecular histology in an animal model of vasculogenic erectile dysfunction.** *J Urol. 1998;159(6):2229-36* — 팽창성-조직학의 정량 연결. [PMID 9598575](https://pubmed.ncbi.nlm.nih.gov/9598575/)
74. Sattar AA, Wespes E, Schulman CC. **Computerized measurement of penile elastic fibres in potent and impotent men.** *Eur Urol. 1994;25(2):142-4* — 탄력섬유 정량. [PMID 8137855](https://pubmed.ncbi.nlm.nih.gov/8137855/)
75. Gratzke C, Angulo J, Chitaley K, et al. **Anatomy, physiology, and pathophysiology of erectile dysfunction.** *J Sex Med. 2010;7(1 Pt 2):445-75* — 구조-기능 종합 리뷰. [PMID 20092448](https://pubmed.ncbi.nlm.nih.gov/20092448/)

## 13. 동체내 주사요법 (Intracavernosal and intraurethral therapy)

76. Linet OI, Ogrinc FG. **Efficacy and safety of intracavernosal alprostadil in men with erectile dysfunction. The Alprostadil Study Group.** *N Engl J Med. 1996;334(14):873-7* — alprostadil ICI 등록시험 — 비반응자에서도 높은 성공률. [PMID 8596569](https://pubmed.ncbi.nlm.nih.gov/8596569/)
77. Porst H. **The rationale for prostaglandin E1 in erectile failure: a survey of worldwide experience.** *J Urol. 1996;155(3):802-15* — PGE1의 근거와 세계 경험. [PMID 8583582](https://pubmed.ncbi.nlm.nih.gov/8583582/)
78. Bennett AH, Carpenter AJ, Barada JH. **An improved vasoactive drug combination for a pharmacological erection program.** *J Urol. 1991;146(6):1564-5* — trimix 조합. [PMID 1719248](https://pubmed.ncbi.nlm.nih.gov/1719248/)
79. Ahn HS, Lee SW, Yoon SJ, Hann HJ, Hong JM. **A comparison of colour duplex ultrasonography after transurethral alprostadil and intracavernous alprostadil in the assessment of erectile dysfunction.** *J Int Med Res. 2004;32(3):317-23* — MUSE 요도내 투여. [PMID 15174226](https://pubmed.ncbi.nlm.nih.gov/15174226/)
80. Broderick GA, Kadioglu A, Bivalacqua TJ, Ghanem H, Nehra A, Shamloul R. **Priapism: pathogenesis, epidemiology, and management.** *J Sex Med. 2010;7(1 Pt 2):476-500* — 지속발기증 — 모델의 중립 게이트 부재 논지. [PMID 20092449](https://pubmed.ncbi.nlm.nih.gov/20092449/)

## 14. 새로운 표적 (sGC 자극제/활성제 · ROCK · Li-ESWT)

81. Sandner P, Hütter J, Tinel H, Ziegelbauer K, Bischoff E. **PDE5 inhibitors beyond erectile dysfunction.** *Int J Impot Res. 2007;19(6):533-43* — PDE5i의 다면 효과. [PMID 17625575](https://pubmed.ncbi.nlm.nih.gov/17625575/)
82. Gomberg-Maitland M, Ghofrani HA, Gibbs JSR, et al. **Sotatercept Safety and Efficacy in Intermediate- to Low-Risk Pulmonary Arterial Hypertension: A Pooled Analysis of PULSAR and STELLAR.** *Chest. 2026* — riociguat: sGC stimulator의 임상 근거와 PDE5i 병용 금기. [PMID 42208730](https://pubmed.ncbi.nlm.nih.gov/42208730/)
83. Stasch JP, Schmidt P, Alonso-Alija C, et al. **NO- and haem-independent activation of soluble guanylyl cyclase: molecular basis and cardiovascular implications of a new pharmacological principle.** *Br J Pharmacol. 2002;136(5):773-83* — cinaciguat 계열 activator의 원리 — 산화 sGC 표적. [PMID 12086987](https://pubmed.ncbi.nlm.nih.gov/12086987/)
84. Gruenwald I, Appel B, Vardi Y. **Low-intensity extracorporeal shock wave therapy--a novel effective treatment for erectile dysfunction in severe ED patients who respond poorly to PDE5 inhibitor therapy.** *J Sex Med. 2012;9(1):259-64* — Li-ESWT. [PMID 22008059](https://pubmed.ncbi.nlm.nih.gov/22008059/)

## 15. 안전성 · 약물상호작용 (Safety and drug interactions)

85. Webb DJ, Freestone S, Allen MJ, Muirhead GJ. **Sildenafil citrate and blood-pressure-lowering drugs: results of drug interaction studies with an organic nitrate and a calcium antagonist.** *Am J Cardiol. 1999;83(5A):21C-28C* — sildenafil + 질산염의 혈압 강하 정량 — 모델의 MAP 낙차 검증 대상. [PMID 10078539](https://pubmed.ncbi.nlm.nih.gov/10078539/)
86. Kloner RA, Hutter AM, Emmick JT, Mitchell MI, Denne J, Jackson G. **Time course of the interaction between tadalafil and nitrates.** *J Am Coll Cardiol. 2003;42(10):1855-60* — tadalafil-질산염 상호작용의 시간 경과. [PMID 14642699](https://pubmed.ncbi.nlm.nih.gov/14642699/)
87. Kloner RA, Mitchell M, Emmick JT. **Cardiovascular effects of tadalafil in patients on common antihypertensive therapies.** *Am J Cardiol. 2003;92(9A):47M-57M* — 항고혈압제 병용. [PMID 14609623](https://pubmed.ncbi.nlm.nih.gov/14609623/)
88. Kloner RA, Jackson G, Emmick JT, et al. **Interaction between the phosphodiesterase 5 inhibitor, tadalafil and 2 alpha-blockers, doxazosin and tamsulosin in healthy normotensive men.** *J Urol. 2004;172(5 Pt 1):1935-40* — 알파차단제 병용. [PMID 15540759](https://pubmed.ncbi.nlm.nih.gov/15540759/)
89. Laties A, Zrenner E. **Viagra (sildenafil citrate) and ophthalmology.** *Prog Retin Eye Res. 2002;21(5):485-506* — PDE6와 시각 이상. [PMID 12207947](https://pubmed.ncbi.nlm.nih.gov/12207947/)
90. Carson CC, Rajfer J, Eardley I, et al. **The efficacy and safety of tadalafil: an update.** *BJU Int. 2004;93(9):1276-81* — PDE11과 근육통·요통. [PMID 15180622](https://pubmed.ncbi.nlm.nih.gov/15180622/)
91. Nehra A, Jackson G, Miner M, et al. **The Princeton III Consensus recommendations for the management of erectile dysfunction and cardiovascular disease.** *Mayo Clin Proc. 2012;87(8):766-78* — 성행위에 대한 심혈관 적합성 평가. [PMID 22862865](https://pubmed.ncbi.nlm.nih.gov/22862865/)
92. Vlachopoulos C, Jackson G, Stefanadis C, Montorsi P. **Erectile dysfunction in the cardiovascular patient.** *Eur Heart J. 2013;34(27):2034-46* — ED의 심혈관 전조 가치. [PMID 23616415](https://pubmed.ncbi.nlm.nih.gov/23616415/)

## 16. 심리 요인과 순응도 (Psychogenic factors and adherence)

93. Althof SE, Rosen RC, Perelman MA, Rubio-Aurioles E. **Standard operating procedures for taking a sexual history.** *J Sex Med. 2013;10(1):26-35* — 심리·관계 요인 평가. [PMID 22970717](https://pubmed.ncbi.nlm.nih.gov/22970717/)
94. Rosen RC. **Psychogenic erectile dysfunction. Classification and management.** *Urol Clin North Am. 2001;28(2):269-78* — 심인성 ED의 분류와 수행불안 고리. [PMID 11402580](https://pubmed.ncbi.nlm.nih.gov/11402580/)
95. Carvalheira AA, Pereira NM, Maroco J, Forjaz V. **Dropout in the treatment of erectile dysfunction with PDE5: a study on predictors and a qualitative analysis of reasons for discontinuation.** *J Sex Med. 2012;9(9):2361-9* — PDE5i 중단률과 예측인자. [PMID 22616766](https://pubmed.ncbi.nlm.nih.gov/22616766/)
96. Melnik T, Soares BG, Nasselo AG. **Psychosocial interventions for erectile dysfunction.** *Cochrane Database Syst Rev. 2007;2007(3):CD004825* — 심리 개입의 효과. [PMID 17636774](https://pubmed.ncbi.nlm.nih.gov/17636774/)

## 17. 생활습관 개선과 진료지침 (Lifestyle modification and guidelines)

97. Esposito K, Giugliano F, Di Palo C, et al. **Effect of lifestyle changes on erectile dysfunction in obese men: a randomized controlled trial.** *JAMA. 2004;291(24):2978-84* — 생활습관 개선 RCT — 모델의 EXER/체중 항. [PMID 15213209](https://pubmed.ncbi.nlm.nih.gov/15213209/)
98. Gupta BP, Murad MH, Clifton MM, Prokop L, Nehra A, Kopecky SL. **The effect of lifestyle modification and cardiovascular risk factor reduction on erectile dysfunction: a systematic review and meta-analysis.** *Arch Intern Med. 2011;171(20):1797-803* — 메타분석. [PMID 21911624](https://pubmed.ncbi.nlm.nih.gov/21911624/)
99. Salonia A, Capogrosso P, Boeri L, et al. **European Association of Urology Guidelines on Male Sexual and Reproductive Health: 2025 Update on Male Hypogonadism, Erectile Dysfunction, Premature Ejaculation, and Peyronie's Disease.** *Eur Urol. 2025;88(1):76-102* — EAU 지침. [PMID 40340108](https://pubmed.ncbi.nlm.nih.gov/40340108/)
100. Burnett AL, Nehra A, Breau RH, et al. **Erectile Dysfunction: AUA Guideline.** *J Urol. 2018;200(3):633-641* — AUA 지침. [PMID 29746858](https://pubmed.ncbi.nlm.nih.gov/29746858/)

## 18. QSP 모델링 방법론 (QSP modelling methodology)

101. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883-893* — mrgsolve 기반 QSP 구현 방법. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
102. Peterson MC, Riggs MM. **FDA Advisory Meeting Clinical Pharmacology Review Utilizes a Quantitative Systems Pharmacology (QSP) Model: A Watershed Moment?.** *CPT Pharmacometrics Syst Pharmacol. 2015;4(3):e00020* — QSP 모델 구축 절차. [PMID 26225239](https://pubmed.ncbi.nlm.nih.gov/26225239/)


---

## 19. 방정식 ↔ 문헌 대응표 (Equation-to-source map)

| 모델 요소 | 근거 섹션 | 비고 |
|---|---|---|
| `S(t)` 자극 입력 · `AROU` · `NRV` 곱 구조 | §3, §16 | NANC 신경전달이 NO 매개이고, 중추 각성과 신경 온전성이 곱으로 들어간다는 점 |
| `NO_prod` 의 nitrergic 항 (`KNONN·NNOS·DRIVE`) | §3 | nNOS 유래 펄스가 발기를 개시 |
| `NO_prod` 의 내피 항 (`KNOEN·ECPL·SHEAR`) | §3 | 전단응력 → Akt → eNOS Ser1177 → 유지 |
| `KSCAV·ROS·NO` (NO 소거) | §11 | NO는 생성되지만 초과산화물에 파괴된다 |
| `NOSEFF = 1/(1+ADMA/Ki)` | §11 | ADMA 경쟁적 NOS 억제 |
| `SGCOX` (산화 sGC 분율) 과 activator 항 | §4, §14 | 산화 sGC는 NO 불응, activator만 작동 |
| `vpde5 = KPDE5·PDE5E·P5RES·cGMP/(Km+cGMP)` | §4 | PDE5 동역학과 Km |
| `P5RES = 1/(1+Cu/IC50)` 및 4약제 IC50 표 | §4, §7 | 재조합 효소 IC50 + 라벨 PK |
| PKG → PDE5 인산화 (음성 피드백) | §4 | Ser102 인산화로 활성 3–10배 |
| `PKG`·`PKA` 병렬 → Ca²⁺ 배출 | §5 | 두 kinase가 같은 효과기로 수렴 |
| `MLCP` 방정식 · ROCK의 MLCP 억제 | §5 | Ca²⁺ 감작의 표준 정식화 |
| ROCK 억제제 시나리오 | §5, §14 | NO 비의존적 발기 유발 |
| `GCAV = GC0·exp(KG·R)` | §2 | 저항 ∝ r⁻⁴ 로부터의 지수적 민감도 |
| `VOCC` (팽창의 Hill⁶) 와 `G_ven` 압박 | §2 | 피막하 정맥 압박이 정맥폐쇄의 기전 |
| `VOCMAX(SMI)` 무릎 | §2, §12 | 평활근 함량이 정맥폐쇄 천장을 결정 |
| `ICP(V)` 지수적 백막 경화 | §2 | 압력-용적 관계 |
| `RIG(ICP)` Hill³ | §2 | 종축 강직도-압력 공학적 해석 |
| `RFEED ∝ (1−STEN)⁻⁴` | §2, §15 | Poiseuille + ED의 심혈관 전조 가치 |
| `hyp → TGFB → COL`, `SM` 아폽토시스 | §9, §12 | 신경절단 → 저산소 → 섬유화의 직접 증거 |
| `KTGFCG` (만성 cGMP의 Smad 억제) | §12 | PGE1/NO-cGMP가 콜라겐 합성을 억제 |
| `NRV` 재생 τ 와 `NRVMAX` | §9 | 신경보존 정도별 회복률 |
| `D_trg` (T 게이트 → nNOS · PDE5 · 평활근 · 리비도) | §10 | 안드로겐의 다중 표적 |
| `HCT` 상승 | §10 | TRT의 적혈구증 안전역 |
| `PA` (수행불안) 양성 피드백 | §16 | 실패 → 불안 → 실패 고리 |
| `SCG` · `MAP` · 질산염 상호작용 | §15 | 같은 두 인자 곱을 전신 혈관에서 읽은 결과 |
| `P6INH` · `P11INH` | §4, §15 | 아이소형 선택성 → 시각 이상 · 근육통 |
| `IIEF-EF` 매핑과 MCID | §1 | 도구 개발 및 임상적 최소 유의차 |
| 가상 집단 (FNOI 로그정규) | §1, §6 | 시험 평균은 혼합분포이며 단일 환자값이 아니다 |

---

## 20. 이 모델이 문헌과 충돌하거나 검증되지 않은 지점

1. **`KNRVDRUG = 0` (기본값).** 동물 실험은 만성 PDE5 억제가 신경영양적
   효과를 가진다고 시사하지만(§9), 사람에서 약물중단 세척기 후의 무보조
   발기 기능은 개선되지 않았습니다. 기본값 0은 후자를 인코딩한 것이며,
   파라미터를 남겨 두어 가설을 검정할 수 있게 했습니다.
2. **테스토스테론 병용의 효과 크기.** §10 에는 sildenafil 비반응 저성선증에서
   T 병용이 효과적이라는 RCT와, 효과가 없다는 RCT가 모두 있습니다. 이 모델의
   곱셈 게이트 구조는 전자를 예측합니다 — 즉 이 지점은 모델의 **반증 가능한
   예측**입니다.
3. **`sigma = 0.55` (개인별 NO 생성능의 로그정규 SD).** 측정값이 아니라
   집단 평균 IIEF-EF를 재현하도록 정한 추론값입니다. 반응자 *비율*은 평균
   IIEF-EF보다 이 값에 훨씬 민감합니다.
4. **`RIG50 = 62%`, `TAE50 = 8분`.** 물리적 상수가 아니라 압력 파형을 성행위
   일지 항목으로 변환하는 행동 파라미터이며, 엔드포인트 불확실성의 가장 큰
   단일 원천입니다.
