# 내장 리슈만편모충증 QSP 모델 — 참고문헌
# Visceral Leishmaniasis QSP Model — References

이 목록은 `vl_qsp_model.dot`(기계론적 지도), `vl_mrgsolve_model.R`(73-ODE 모델),
`vl_reference_model.py`(독립 검증 구현)에 들어간 구조적 선택과 파라미터의 근거입니다.

**모든 항목은 PubMed E-utilities로 직접 조회한 실제 레코드입니다** — 제목·저자·
저널·연도·PMID는 조회 결과를 그대로 옮겼고, 기억에 의존해 작성한 인용은 없습니다.
조회 스크립트는 41개 질의(주제별)로 esearch/esummary를 호출했고 중복을 제거해
아래 목록을 만들었습니다.

Every entry below is a real record retrieved directly from PubMed via the NCBI
E-utilities API; titles, authors, journals, years and PMIDs are transcribed
from that retrieval rather than written from memory.

---

## 이 모델이 문헌에서 가져온 것과, 문헌에 없는 것 (What is sourced, and what is not)

**문헌에서 직접 가져온 것 (sourced):** 약물별 혈장 약동학 파라미터(청소율,
분포용적, 흡수속도, 반감기), 알로메트릭 지수(0.75/1.0), 밀테포신 소아 노출 부족,
비하르 안티몬 내성 빈도, 각 요법의 임상 완치율과 이상반응 빈도, IL-10이 활성
VL의 표지라는 사실, 레이슈마닌 피부반응이 치료 후 전환된다는 사실, PKDL이
치료 후에 나타난다는 사실.

**모델 수준의 가설이며 관측이 아닌 것 (model-level hypotheses, NOT observations):**
세포내 대식세포 구획의 EC50 값(자유 세포내 농도는 측정 불가하므로 임상 완치율에
맞춘 보정 파라미터), 리포솜 청소율의 MPS 분배 비율(FSP_A/FLI_A/FBM_A/FSK_A),
정지기(quiescent) 개체군의 크기와 EC50 배수, TMEM·IL-10·MACT 모듈의 모든 속도
상수, 그리고 CD4에 따른 분리선(separatrix)의 절대 위치. 이들은 구조적 주장을
검정 가능하게 만들기 위한 것이며, 보고된 비율·순서·창(window)은 절대값보다
훨씬 견고합니다.

---

## 1. 역학·질병부담·퇴치 프로그램 / Epidemiology, burden and the elimination programme

Where the disease is, how much of it there is, and why the 2030 elimination target makes the PKDL reservoir the central pharmacological problem.

- Addisu A, Bayiyana A, Reis Cunha JL, et al. **Systemic inflammatory markers of visceral leishmaniasis treatment response in East Africa.** *PLoS Negl Trop Dis*. 2026. [PMID 41758879](https://pubmed.ncbi.nlm.nih.gov/41758879/)
- Pareyn M, Alves F, Burza S, et al. **Leishmaniasis.** *Nat Rev Dis Primers*. 2025. [PMID 41266459](https://pubmed.ncbi.nlm.nih.gov/41266459/)
- Sundar S. **The story of elimination of visceral leishmaniasis (kala-azar) in India-Challenges towards sustainment.** *PLoS Negl Trop Dis*. 2025. [PMID 40828749](https://pubmed.ncbi.nlm.nih.gov/40828749/)
- Saout M, Pinheiro T, da Costa V, et al. **Fatal case of imported visceral leishmaniasis in a dog caused by Leishmania infantum in French Guiana.** *Vet Parasitol Reg Stud Reports*. 2024. [PMID 39326963](https://pubmed.ncbi.nlm.nih.gov/39326963/)
- Musa AM, Mbui J, Mohammed R, et al. **Paromomycin and Miltefosine Combination as an Alternative to Treat Patients With Visceral Leishmaniasis in Eastern Africa: A Randomized, Controlled, Multicountry Trial.** *Clin Infect Dis*. 2023. [PMID 36164254](https://pubmed.ncbi.nlm.nih.gov/36164254/)
- Kumar A, Singh VK, Tiwari R, et al. **Post kala-azar dermal leishmaniasis in the Indian sub-continent: challenges and strategies for elimination.** *Front Immunol*. 2023. [PMID 37638047](https://pubmed.ncbi.nlm.nih.gov/37638047/)
- Garlapati R, Iniguez E, Serafim TD, et al. **Towards a Sustainable Vector-Control Strategy in the Post Kala-Azar Elimination Era.** *Front Cell Infect Microbiol*. 2021. [PMID 33768013](https://pubmed.ncbi.nlm.nih.gov/33768013/)
- Ghosh P, Roy P, Chaudhuri SJ, et al. **Epidemiology of Post-Kala-azar Dermal Leishmaniasis.** *Indian J Dermatol*. 2021. [PMID 33911289](https://pubmed.ncbi.nlm.nih.gov/33911289/)
- Thakur CP, Thakur M. **Accelerating kala-azar elimination in India.** *Indian J Med Res*. 2020. [PMID 34145092](https://pubmed.ncbi.nlm.nih.gov/34145092/)
- Carnielli JBT, Monti-Rocha R, Costa DL, et al. **Natural Resistance of Leishmania infantum to Miltefosine Contributes to the Low Efficacy in the Treatment of Visceral Leishmaniasis in Brazil.** *Am J Trop Med Hyg*. 2019. [PMID 31436148](https://pubmed.ncbi.nlm.nih.gov/31436148/)
- Gonçalves AAM, Leite JC, Resende LA, et al. **An Overview of Immunotherapeutic Approaches Against Canine Visceral Leishmaniasis: What Has Been Tested on Dogs and a New Perspective on Improving Treatment Efficacy.** *Front Cell Infect Microbiol*. 2019. [PMID 31921703](https://pubmed.ncbi.nlm.nih.gov/31921703/)
- Dos Santos Nogueira F, Avino VC, Galvis-Ovallos F, et al. **Use of miltefosine to treat canine visceral leishmaniasis caused by Leishmania infantum in Brazil.** *Parasit Vectors*. 2019. [PMID 30736866](https://pubmed.ncbi.nlm.nih.gov/30736866/)
- Carnielli JBT, Crouch K, Forrester S, et al. **A Leishmania infantum genetic marker associated with miltefosine treatment failure for visceral leishmaniasis.** *EBioMedicine*. 2018. [PMID 30268832](https://pubmed.ncbi.nlm.nih.gov/30268832/)
- Burza S, Croft SL, Boelaert M. **Leishmaniasis.** *Lancet*. 2018. [PMID 30126638](https://pubmed.ncbi.nlm.nih.gov/30126638/)
- Bhattacharya SK, Dash AP. **Elimination of Kala-Azar from the Southeast Asia Region.** *Am J Trop Med Hyg*. 2017. [PMID 28115678](https://pubmed.ncbi.nlm.nih.gov/28115678/)
- Kimutai R, Musa AM, Njoroge S, et al. **Safety and Effectiveness of Sodium Stibogluconate and Paromomycin Combination for the Treatment of Visceral Leishmaniasis in Eastern Africa: Results from a Pharmacovigilance Programme.** *Clin Drug Investig*. 2017. [PMID 28066878](https://pubmed.ncbi.nlm.nih.gov/28066878/)
- Dorlo TPC, Kip AE, Younis BM, et al. **Visceral leishmaniasis relapse hazard is linked to reduced miltefosine exposure in patients from Eastern Africa: a population pharmacokinetic/pharmacodynamic study.** *J Antimicrob Chemother*. 2017. [PMID 28961737](https://pubmed.ncbi.nlm.nih.gov/28961737/)
- Gadisa E, Tsegaw T, Abera A, et al. **Eco-epidemiology of visceral leishmaniasis in Ethiopia.** *Parasit Vectors*. 2015. [PMID 26187584](https://pubmed.ncbi.nlm.nih.gov/26187584/)
- Ready PD. **Epidemiology of visceral leishmaniasis.** *Clin Epidemiol*. 2014. [PMID 24833919](https://pubmed.ncbi.nlm.nih.gov/24833919/)
- Khalil EA, Weldegebreal T, Younis BM, et al. **Safety and efficacy of single dose versus multiple doses of AmBisome for treatment of visceral leishmaniasis in eastern Africa: a randomised trial.** *PLoS Negl Trop Dis*. 2014. [PMID 24454970](https://pubmed.ncbi.nlm.nih.gov/24454970/)
- Hailu A, Musa A, Wasunna M, et al. **Geographical variation in the response of visceral leishmaniasis to paromomycin in East Africa: a multicentre, open-label, randomized trial.** *PLoS Negl Trop Dis*. 2010. [PMID 21049059](https://pubmed.ncbi.nlm.nih.gov/21049059/)
- Mondal D, Singh SP, Kumar N, et al. **Visceral leishmaniasis elimination programme in India, Bangladesh, and Nepal: reshaping the case finding/case management strategy.** *PLoS Negl Trop Dis*. 2009. [PMID 19159009](https://pubmed.ncbi.nlm.nih.gov/19159009/)
- Sundar S, Mondal D, Rijal S, et al. **Implementation research to support the initiative on the elimination of kala azar from Bangladesh, India and Nepal--the challenges for diagnosis and treatment.** *Trop Med Int Health*. 2008. [PMID 18290995](https://pubmed.ncbi.nlm.nih.gov/18290995/)

## 2. 매개체·전파·접종량 / Vector, transmission and inoculum

The sandfly end of the cycle. The inoculum size and the promastigote secretory gel set the starting condition of the model (`stage="naive"`).

- Roy L, Cloots K, Uranw S, et al. **The ongoing risk of Leishmania donovani transmission in eastern Nepal: an entomological investigation during the elimination era.** *Parasit Vectors*. 2023. [PMID 37932813](https://pubmed.ncbi.nlm.nih.gov/37932813/)
- Hendrickx S, Van Bockstal L, Aslan H, et al. **Transmission potential of paromomycin-resistant Leishmania infantum and Leishmania donovani.** *J Antimicrob Chemother*. 2020. [PMID 31886863](https://pubmed.ncbi.nlm.nih.gov/31886863/)
- Seblova V, Dujardin JC, Rijal S, et al. **ISC1, a new Leishmania donovani population emerging in the Indian sub-continent: Vector competence of Phlebotomus argentipes.** *Infect Genet Evol*. 2019. [PMID 31629887](https://pubmed.ncbi.nlm.nih.gov/31629887/)
- Sadlova J, Myskova J, Lestinova T, et al. **Leishmania donovani development in Phlebotomus argentipes: comparison of promastigote- and amastigote-initiated infections.** *Parasitology*. 2017. [PMID 27876097](https://pubmed.ncbi.nlm.nih.gov/27876097/)
- Cameron MM, Acosta-Serrano A, Bern C, et al. **Understanding the transmission dynamics of Leishmania donovani to provide robust evidence for interventions to eliminate visceral leishmaniasis in Bihar, India.** *Parasit Vectors*. 2016. [PMID 26812963](https://pubmed.ncbi.nlm.nih.gov/26812963/)
- Ghosh KN, Mukhopadhyay J. **The effect of anti-sandfly saliva antibodies on Phlebotomus argentipes and Leishmania donovani.** *Int J Parasitol*. 1998. [PMID 9512990](https://pubmed.ncbi.nlm.nih.gov/9512990/)

## 3. 세포내 생존 지위와 대식세포 생물학 / The intracellular niche and macrophage biology

Why the drug has to get inside a macrophage before it can do anything — the single structural commitment of this model.

- Hong J, Mukherjee B, Sanjoba C, et al. **Upregulation of ATP6V0D2 benefits intracellular survival of Leishmania donovani in erythrocytes-engulfing macrophages.** *Front Cell Infect Microbiol*. 2024. [PMID 38357442](https://pubmed.ncbi.nlm.nih.gov/38357442/)
- Lopes ME, Dos Santos LM, Sacks D, et al. **Resistance Against Leishmania major Infection Depends on Microbiota-Guided Macrophage Activation.** *Front Immunol*. 2021. [PMID 34745100](https://pubmed.ncbi.nlm.nih.gov/34745100/)
- Badirzadeh A, Montakhab-Yeganeh H, Miandoabi T. **Arginase/nitric oxide modifications using live non-pathogenic Leishmania tarentolae as an effective delivery system inside the mammalian macrophages.** *J Parasit Dis*. 2021. [PMID 33746388](https://pubmed.ncbi.nlm.nih.gov/33746388/)
- Giri S, Shaha C. **Leishmania donovani parasite requires Atg8 protein for infectivity and survival under stress.** *Cell Death Dis*. 2019. [PMID 31649242](https://pubmed.ncbi.nlm.nih.gov/31649242/)
- Heidari-Kharaji M, Fallah-Omrani V, Badirzadeh A, et al. **Sambucus ebulus extract stimulates cellular responses in cutaneous leishmaniasis.** *Parasite Immunol*. 2019. [PMID 30472741](https://pubmed.ncbi.nlm.nih.gov/30472741/)
- Parmar N, Chandrakar P, Vishwakarma P, et al. **Leishmania donovani Exploits Tollip, a Multitasking Protein, To Impair TLR/IL-1R Signaling for Its Survival in the Host.** *J Immunol*. 2018. [PMID 29907707](https://pubmed.ncbi.nlm.nih.gov/29907707/)
- Acuña SM, Aoki JI, Laranjeira-Silva MF, et al. **Arginase expression modulates nitric oxide production in Leishmania (Leishmania) amazonensis.** *PLoS One*. 2017. [PMID 29135983](https://pubmed.ncbi.nlm.nih.gov/29135983/)
- Yau WL, Pescher P, MacDonald A, et al. **The Leishmania donovani chaperone cyclophilin 40 is essential for intracellular infection independent of its stage-specific phosphorylation status.** *Mol Microbiol*. 2014. [PMID 24811325](https://pubmed.ncbi.nlm.nih.gov/24811325/)
- Balaña-Fouce R, Calvo-Álvarez E, Álvarez-Velilla R, et al. **Role of trypanosomatid's arginase in polyamine biosynthesis and pathogenesis.** *Mol Biochem Parasitol*. 2012. [PMID 22033378](https://pubmed.ncbi.nlm.nih.gov/22033378/)
- Mukbel RM, Patten C Jr, Gibson K, et al. **Macrophage killing of Leishmania amazonensis amastigotes requires both nitric oxide and superoxide.** *Am J Trop Med Hyg*. 2007. [PMID 17426168](https://pubmed.ncbi.nlm.nih.gov/17426168/)
- Mundodi V, Kucknoor AS, Gedamu L. **Role of Leishmania (Leishmania) chagasi amastigote cysteine protease in intracellular parasite survival: studies by gene disruption and antisense mRNA inhibition.** *BMC Mol Biol*. 2005. [PMID 15691375](https://pubmed.ncbi.nlm.nih.gov/15691375/)
- Moore KJ, Labrecque S, Matlashewski G. **Alteration of Leishmania donovani infection levels by selective impairment of macrophage signal transduction.** *J Immunol*. 1993. [PMID 8482844](https://pubmed.ncbi.nlm.nih.gov/8482844/)
- Chang KP, Fong D. **Cell biology of host-parasite membrane interactions in leishmaniasis.** *Ciba Found Symp*. 1983. [PMID 6357669](https://pubmed.ncbi.nlm.nih.gov/6357669/)
- Chang KP, Dwyer DM. **Leishmania donovani. Hamster macrophage interactions in vitro: cell entry, intracellular survival, and multiplication of amastigotes.** *J Exp Med*. 1978. [PMID 564391](https://pubmed.ncbi.nlm.nih.gov/564391/)

## 4. 면역병태생리: IL-10, Th1, 그리고 프라이밍 경쟁 / Immunopathogenesis: IL-10, Th1 and the priming race

The evidence behind the TMEM / IL-10 / MACT module and the claim that subclinical versus clinical infection is a race rather than a switch.

- Amorim Sacramento L, Gonzalez-Lombana C, Scott P. **Malnutrition disrupts adaptive immunity during visceral leishmaniasis by enhancing IL-10 production.** *PLoS Pathog*. 2024. [PMID 39527629](https://pubmed.ncbi.nlm.nih.gov/39527629/)
- Bragato JP, Rebech GT, Freitas JH, et al. **miRNA-21 regulates CD69 and IL-10 expression in canine leishmaniasis.** *PLoS One*. 2022. [PMID 35324917](https://pubmed.ncbi.nlm.nih.gov/35324917/)
- Teles LF, Viana AG, Cardoso MS, et al. **Evaluation of medullary cytokine expression and clinical and laboratory aspects in severe human visceral leishmaniasis.** *Parasite Immunol*. 2021. [PMID 34558674](https://pubmed.ncbi.nlm.nih.gov/34558674/)
- Moulik S, Karmakar J, Joshi S, et al. **Status of IL-4 and IL-10 driven markers in experimental models of Visceral Leishmaniasis.** *Parasite Immunol*. 2021. [PMID 32734677](https://pubmed.ncbi.nlm.nih.gov/32734677/)
- Murray HW. **Targeting IL-27 and/or IL-10 in Experimental Murine Visceral Leishmaniasis.** *Am J Trop Med Hyg*. 2020. [PMID 32815498](https://pubmed.ncbi.nlm.nih.gov/32815498/)
- Roy S, Mukhopadhyay D, Mukherjee S, et al. **An IL-10 dominant polarization of monocytes is a feature of Indian Visceral Leishmaniasis.** *Parasite Immunol*. 2018. [PMID 29745990](https://pubmed.ncbi.nlm.nih.gov/29745990/)
- Murray HW, Lu CM, Mauze S, et al. **Interleukin-10 (IL-10) in experimental visceral leishmaniasis and IL-10 receptor blockade as immunotherapy.** *Infect Immun*. 2002. [PMID 12379707](https://pubmed.ncbi.nlm.nih.gov/12379707/)
- Melby PC, Tryon VV, Chandrasekar B, et al. **Cloning of Syrian hamster (Mesocricetus auratus) cytokine cDNAs and analysis of cytokine mRNA expression in experimental visceral leishmaniasis.** *Infect Immun*. 1998. [PMID 9573100](https://pubmed.ncbi.nlm.nih.gov/9573100/)

## 5. 혈액학적 이상·비장비대·HLH / Haematology, splenomegaly and HLH

The two independent mechanisms — hypersplenic destruction and marrow occupancy — that the HGB / PLT / WBC equations keep separate.

- Bansal A, Gupta A, Gupta M, et al. **Visceral leishmaniasis concealed under a diagnosis of haemophagocytic lymphohistiocytosis.** *BMJ Case Rep*. 2025. [PMID 40588301](https://pubmed.ncbi.nlm.nih.gov/40588301/)
- Freitas AFA, Neto ASL, de Moura CMC, et al. **Relationship Between Hepcidin, Iron Metabolism, Inflammation and Hypersplenism in Anaemia of Kala-Azar.** *Parasite Immunol*. 2025. [PMID 40637365](https://pubmed.ncbi.nlm.nih.gov/40637365/)
- Scarpini S, Dondi A, Totaro C, et al. **Visceral Leishmaniasis: Epidemiology, Diagnosis, and Treatment Regimens in Different Geographical Areas with a Focus on Pediatrics.** *Microorganisms*. 2022. [PMID 36296164](https://pubmed.ncbi.nlm.nih.gov/36296164/)
- Chen H, Fan C, Gao H, et al. **Leishmaniasis Diagnosis via Metagenomic Next-Generation Sequencing.** *Front Cell Infect Microbiol*. 2020. [PMID 33072623](https://pubmed.ncbi.nlm.nih.gov/33072623/)
- Kalmi G, Vignon-Pennamen MD, Ram-Wolff C, et al. **Visceral leishmaniasis in patients with lymphoma: Case reports and review of the literature.** *Medicine (Baltimore)*. 2020. [PMID 33157924](https://pubmed.ncbi.nlm.nih.gov/33157924/)
- Ahmed MA, Ahmed AA, Omar SM, et al. **Epidemiology of visceral leishmaniasis among children in Gadarif hospital, eastern Sudan.** *BMC Public Health*. 2016. [PMID 27927185](https://pubmed.ncbi.nlm.nih.gov/27927185/)
- Prajapati R, Kumar A, Sharma P, et al. **A Rare Presentation of Leishmaniasis.** *J Clin Exp Hepatol*. 2016. [PMID 27493461](https://pubmed.ncbi.nlm.nih.gov/27493461/)
- Koster KL, Laws HJ, Troeger A, et al. **Visceral Leishmaniasis as a Possible Reason for Pancytopenia.** *Front Pediatr*. 2015. [PMID 26176005](https://pubmed.ncbi.nlm.nih.gov/26176005/)
- Gallina V, Binazzi R, Golemi A, et al. **Imported visceral leishmaniasis - unexpected bone marrow diagnosis in a patient with fever, pancytopenia, and splenomegaly.** *Am J Blood Res*. 2014. [PMID 25755910](https://pubmed.ncbi.nlm.nih.gov/25755910/)
- Bhauwala A, Dhir V, Gill A, et al. **Visceral leishmaniasis complicated by haemophagocytosis.** *Natl Med J India*. 2013. [PMID 25017849](https://pubmed.ncbi.nlm.nih.gov/25017849/)
- Cascio A, Colomba C. **[Childhood Mediterranean visceral leishmaniasis].** *Infez Med*. 2003. [PMID 12719664](https://pubmed.ncbi.nlm.nih.gov/12719664/)
- Al Sineidi K, Wali YA, Pathare AV, et al. **Visceral leishmaniasis and haemophagocytic syndrome in an Omani child.** *J Sci Res Med Sci*. 2002. [PMID 24019726](https://pubmed.ncbi.nlm.nih.gov/24019726/)
- al-Jurayyan NA, al-Nasser MN, al-Fawaz IM, et al. **The haematological manifestations of visceral leishmaniasis in infancy and childhood.** *J Trop Pediatr*. 1995. [PMID 7636932](https://pubmed.ncbi.nlm.nih.gov/7636932/)

## 6. 진단, 치료판정, 혈청검사의 한계 / Diagnosis, test of cure and the serology trap

Why the model reports a splenic-aspirate grade and why it predicts that rK39 cannot be used to call relapse.

- Hagos DG, Schallig HDFH, Kiros YK, et al. **Performance of rapid rk39 tests for the diagnosis of visceral leishmaniasis in Ethiopia: a systematic review and meta-analysis.** *BMC Infect Dis*. 2021. [PMID 34789175](https://pubmed.ncbi.nlm.nih.gov/34789175/)
- Verrest L, Kip AE, Musa AM, et al. **Blood Parasite Load as an Early Marker to Predict Treatment Response in Visceral Leishmaniasis in Eastern Africa.** *Clin Infect Dis*. 2021. [PMID 33580234](https://pubmed.ncbi.nlm.nih.gov/33580234/)
- Mann S, Frasca K, Scherrer S, et al. **A Review of Leishmaniasis: Current Knowledge and Future Directions.** *Curr Trop Med Rep*. 2021. [PMID 33747716](https://pubmed.ncbi.nlm.nih.gov/33747716/)
- van Griensven J, Diro E. **Visceral Leishmaniasis: Recent Advances in Diagnostics and Treatment Regimens.** *Infect Dis Clin North Am*. 2019. [PMID 30712769](https://pubmed.ncbi.nlm.nih.gov/30712769/)
- Garg S, Tripathi R, Tripathi K. **Oral mucosal involvement in visceral leishmaniasis.** *Asian Pac J Trop Med*. 2013. [PMID 23375044](https://pubmed.ncbi.nlm.nih.gov/23375044/)
- Quinnell RJ, Carson C, Reithinger R, et al. **Evaluation of rK39 rapid diagnostic tests for canine visceral leishmaniasis: longitudinal study and meta-analysis.** *PLoS Negl Trop Dis*. 2013. [PMID 23326615](https://pubmed.ncbi.nlm.nih.gov/23326615/)
- Maia Z, Lírio M, Mistro S, et al. **Comparative study of rK39 Leishmania antigen for serodiagnosis of visceral leishmaniasis: systematic review with meta-analysis.** *PLoS Negl Trop Dis*. 2012. [PMID 22303488](https://pubmed.ncbi.nlm.nih.gov/22303488/)
- Sundar S, Jha TK, Thakur CP, et al. **Oral miltefosine for Indian visceral leishmaniasis.** *N Engl J Med*. 2002. [PMID 12456849](https://pubmed.ncbi.nlm.nih.gov/12456849/)
- Mukherjee C, Mitra PK, Dasgupta A, et al. **Clinicopathological study of clinically undiagnosed cases of kala-azar with special reference to grading of parasites.** *J Indian Med Assoc*. 1995. [PMID 8522820](https://pubmed.ncbi.nlm.nih.gov/8522820/)
- Chulay JD, Bryceson AD. **Quantitation of amastigotes of Leishmania donovani in smears of splenic aspirates from patients with visceral leishmaniasis.** *Am J Trop Med Hyg*. 1983. [PMID 6859397](https://pubmed.ncbi.nlm.nih.gov/6859397/)

## 7. 리포솜 암포테리신 B: 약동학과 조직 표적화 / Liposomal amphotericin B: PK and tissue targeting

The delivery step that claim 1 is built on: the liposome is cleared by the mononuclear phagocyte system, which is the cell that harbours the parasite.

- Chu WY, Singh OP, Sundar S, et al. **Bidirectional Interaction Between Liposomal Amphotericin B Pharmacokinetics and Parasite Dynamics in Patients With Post-Kala-Azar Dermal Leishmaniasis: Potential Implications for Optimal Dosing.** *Clin Pharmacol Ther*. 2026. [PMID 41288134](https://pubmed.ncbi.nlm.nih.gov/41288134/)
- Siddiqui NA, Ansari MZ, Sinha SK, et al. **Treatment Outcomes of Single-Dose Liposomal Amphotericin B-Treated Visceral Leishmaniasis Patients and Factors Affecting Outcome in Bihar, India.** *Am J Trop Med Hyg*. 2024. [PMID 39378874](https://pubmed.ncbi.nlm.nih.gov/39378874/)
- Vaish E, Gupta KK, Ansari SA, et al. **Amphotericin B induced pancytopenia.** *J Family Med Prim Care*. 2022. [PMID 36505588](https://pubmed.ncbi.nlm.nih.gov/36505588/)
- Ekram MR, Amin MR, Hasan MJ, et al. **Efficacy and safety of single-dose liposomal amphotericin B in patients with visceral leishmaniasis in Bangladesh: a real-life experience.** *J Parasit Dis*. 2021. [PMID 34789971](https://pubmed.ncbi.nlm.nih.gov/34789971/)
- Sundar S, Singh A, Agrawal N, et al. **Effectiveness of Single-Dose Liposomal Amphotericin B in Visceral Leishmaniasis in Bihar.** *Am J Trop Med Hyg*. 2019. [PMID 31436156](https://pubmed.ncbi.nlm.nih.gov/31436156/)
- Adler-Moore J, Lewis RE, Brüggemann RJM, et al. **Preclinical Safety, Tolerability, Pharmacokinetics, Pharmacodynamics, and Antifungal Activity of Liposomal Amphotericin B.** *Clin Infect Dis*. 2019. [PMID 31222254](https://pubmed.ncbi.nlm.nih.gov/31222254/)
- Wijnant GJ, Van Bocxlaer K, Fortes Francisco A, et al. **Local Skin Inflammation in Cutaneous Leishmaniasis as a Source of Variable Pharmacokinetics and Therapeutic Efficacy of Liposomal Amphotericin B.** *Antimicrob Agents Chemother*. 2018. [PMID 30082295](https://pubmed.ncbi.nlm.nih.gov/30082295/)
- Rodrigo C, Weeratunga P, Fernando SD, et al. **Amphotericin B for treatment of visceral leishmaniasis: systematic review and meta-analysis of prospective comparative clinical studies including dose-ranging studies.** *Clin Microbiol Infect*. 2018. [PMID 29138100](https://pubmed.ncbi.nlm.nih.gov/29138100/)
- Patere SN, Pathak PO, Kumar Shukla A, et al. **Surface-Modified Liposomal Formulation of Amphotericin B: In vitro Evaluation of Potential Against Visceral Leishmaniasis.** *AAPS PharmSciTech*. 2017. [PMID 27222025](https://pubmed.ncbi.nlm.nih.gov/27222025/)
- Voak AA, Harris A, Qaiser Z, et al. **Pharmacodynamics and Biodistribution of Single-Dose Liposomal Amphotericin B at Different Stages of Experimental Visceral Leishmaniasis.** *Antimicrob Agents Chemother*. 2017. [PMID 28630200](https://pubmed.ncbi.nlm.nih.gov/28630200/)
- Kshirsagar NA. **Single-dose liposomal amphotericin B for visceral leishmaniasis.** *Lancet Glob Health*. 2014. [PMID 25103058](https://pubmed.ncbi.nlm.nih.gov/25103058/)
- Hamill RJ. **Amphotericin B formulations: a comparative review of efficacy and toxicity.** *Drugs*. 2013. [PMID 23729001](https://pubmed.ncbi.nlm.nih.gov/23729001/)
- Sundar S, Chakravarty J, Agarwal D, et al. **Single-dose liposomal amphotericin B for visceral leishmaniasis in India.** *N Engl J Med*. 2010. [PMID 20147716](https://pubmed.ncbi.nlm.nih.gov/20147716/)
- Smith PJ, Olson JA, Constable D, et al. **Effects of dosing regimen on accumulation, retention and prophylactic efficacy of liposomal amphotericin B.** *J Antimicrob Chemother*. 2007. [PMID 17400589](https://pubmed.ncbi.nlm.nih.gov/17400589/)
- Groll AH, Lyman CA, Petraitis V, et al. **Compartmentalized intrapulmonary pharmacokinetics of amphotericin B and its lipid formulations.** *Antimicrob Agents Chemother*. 2006. [PMID 17005824](https://pubmed.ncbi.nlm.nih.gov/17005824/)
- Goldsmith DR, Perry CM. **Amphotericin B lipid complex: in visceral leishmaniasis.** *Drugs*. 2004. [PMID 15329037](https://pubmed.ncbi.nlm.nih.gov/15329037/)
- Sundar S, Jha TK, Thakur CP, et al. **Single-dose liposomal amphotericin B in the treatment of visceral leishmaniasis in India: a multicenter study.** *Clin Infect Dis*. 2003. [PMID 12955641](https://pubmed.ncbi.nlm.nih.gov/12955641/)

## 8. 암포테리신 B 데옥시콜레이트와 신독성 / Amphotericin B deoxycholate and nephrotoxicity

The other integral of the same dose: free plasma drug, renal cortex, tubular injury, potassium and magnesium wasting.

- Pasqualotto AC, Lana DD, Godoy CSM, et al. **Single High Dose of Liposomal Amphotericin B in Human Immunodeficiency Virus/AIDS-Related Disseminated Histoplasmosis: A Randomized Trial.** *Clin Infect Dis*. 2023. [PMID 37232940](https://pubmed.ncbi.nlm.nih.gov/37232940/)
- Stott KE, Moyo M, Ahmadu A, et al. **Population pharmacokinetics of liposomal amphotericin B in adults with HIV-associated cryptococcal meningoencephalitis.** *J Antimicrob Chemother*. 2022. [PMID 36411251](https://pubmed.ncbi.nlm.nih.gov/36411251/)
- Kably B, Launay M, Derobertmasure A, et al. **Antifungal Drugs TDM: Trends and Update.** *Ther Drug Monit*. 2022. [PMID 34923544](https://pubmed.ncbi.nlm.nih.gov/34923544/)
- Tragiannidis A, Gkampeta A, Vousvouki M, et al. **Antifungal agents and the kidney: pharmacokinetics, clinical nephrotoxicity, and interactions.** *Expert Opin Drug Saf*. 2021. [PMID 33896310](https://pubmed.ncbi.nlm.nih.gov/33896310/)
- Alenazi SA, Elmorsy E, Al-Ghafari A, et al. **Effect of amphotericin B-deoxycholate (Fungizone) on the mitochondria of Wistar rats' renal proximal tubules cells.** *J Appl Toxicol*. 2021. [PMID 33740284](https://pubmed.ncbi.nlm.nih.gov/33740284/)
- Liamis G, Hoorn EJ, Florentin M, et al. **An overview of diagnosis and management of drug-induced hypomagnesemia.** *Pharmacol Res Perspect*. 2021. [PMID 34278747](https://pubmed.ncbi.nlm.nih.gov/34278747/)
- Le T, Ly VT, Thu NTM, et al. **Population Pharmacodynamics of Amphotericin B Deoxycholate for Disseminated Infection Caused by Talaromyces marneffei.** *Antimicrob Agents Chemother*. 2019. [PMID 30420478](https://pubmed.ncbi.nlm.nih.gov/30420478/)
- Groll AH, Rijnders BJA, Walsh TJ, et al. **Clinical Pharmacokinetics, Pharmacodynamics, Safety and Efficacy of Liposomal Amphotericin B.** *Clin Infect Dis*. 2019. [PMID 31222253](https://pubmed.ncbi.nlm.nih.gov/31222253/)
- Stone NR, Bicanic T, Salim R, et al. **Liposomal Amphotericin B (AmBisome(®)): A Review of the Pharmacokinetics, Pharmacodynamics, Clinical Experience and Future Directions.** *Drugs*. 2016. [PMID 26818726](https://pubmed.ncbi.nlm.nih.gov/26818726/)
- Trejtnar F, Mandíková J, Kočíncová J, et al. **Renal handling of amphotericin B and amphotericin B-deoxycholate and potential renal drug-drug interactions with selected antivirals.** *Antimicrob Agents Chemother*. 2014. [PMID 24957831](https://pubmed.ncbi.nlm.nih.gov/24957831/)
- Chai LY, Netea MG, Tai BC, et al. **An elevated pro-inflammatory cytokine response is linked to development of amphotericin B-induced nephrotoxicity.** *J Antimicrob Chemother*. 2013. [PMID 23557927](https://pubmed.ncbi.nlm.nih.gov/23557927/)
- Moen MD, Lyseng-Williamson KA, Scott LJ. **Liposomal amphotericin B: a review of its use as empirical therapy in febrile neutropenia and in the treatment of invasive fungal infections.** *Drugs*. 2009. [PMID 19275278](https://pubmed.ncbi.nlm.nih.gov/19275278/)
- Deray G. **Amphotericin B nephrotoxicity.** *J Antimicrob Chemother*. 2002. [PMID 11801579](https://pubmed.ncbi.nlm.nih.gov/11801579/)
- Sabra R, Zeinoun N, Sharaf LH, et al. **Role of humoral mediators in, and influence of a liposomal formulation on, acute amphotericin B nephrotoxicity.** *Pharmacol Toxicol*. 2001. [PMID 11322173](https://pubmed.ncbi.nlm.nih.gov/11322173/)
- Edoute Y, Silbershatz I. **[Hypokalemia induced by amphotericin B].** *Harefuah*. 1999. [PMID 10959350](https://pubmed.ncbi.nlm.nih.gov/10959350/)
- Sabra R, Branch RA. **Amphotericin B nephrotoxicity.** *Drug Saf*. 1990. [PMID 2182052](https://pubmed.ncbi.nlm.nih.gov/2182052/)

## 9. 밀테포신: 약동학, 소아 노출, 내성 / Miltefosine: PK, paediatric exposure and resistance

The allometry trap (claim 2), the long terminal tail, adherence, and LdMT/LdRos3 transporter loss.

- Mazariegos Herrera A, Karlsson MO, Svensson EM, et al. **Weight-band-based simplification of oral allometric miltefosine dosing in paediatric patients with visceral leishmaniasis.** *J Antimicrob Chemother*. 2026. [PMID 41574656](https://pubmed.ncbi.nlm.nih.gov/41574656/)
- Ghousepeer GD, Singh SK, Singh S, et al. **Potential efficacy of repurposed drugs for treating drug-resistant visceral leishmaniasis.** *Naunyn Schmiedebergs Arch Pharmacol*. 2026. [PMID 42363946](https://pubmed.ncbi.nlm.nih.gov/42363946/)
- Masini T, Maia-Elkhoury ANS, Mondal D, et al. **Optimization and prioritization of paediatric drugs for visceral leishmaniasis.** *Front Pediatr*. 2025. [PMID 40980123](https://pubmed.ncbi.nlm.nih.gov/40980123/)
- Chu WY, Verrest L, Younis BM, et al. **Disease-Specific Differences in Pharmacokinetics of Paromomycin and Miltefosine Between Post-Kala-Azar Dermal Leishmaniasis and Visceral Leishmaniasis Patients in Eastern Africa.** *J Infect Dis*. 2024. [PMID 39166299](https://pubmed.ncbi.nlm.nih.gov/39166299/)
- Verrest L, Roseboom IC, Wasunna M, et al. **Population pharmacokinetics of a combination of miltefosine and paromomycin in Eastern African children and adults with visceral leishmaniasis.** *J Antimicrob Chemother*. 2023. [PMID 37726401](https://pubmed.ncbi.nlm.nih.gov/37726401/)
- Saboia-Vahia L, Cuervo P, Wiśniewski JR, et al. **In-Depth Quantitative Proteomics Characterization of In Vitro Selected Miltefosine Resistance in Leishmania infantum.** *Proteomes*. 2022. [PMID 35466238](https://pubmed.ncbi.nlm.nih.gov/35466238/)
- Van Bockstal L, Bulté D, Hendrickx S, et al. **Impact of clinically acquired miltefosine resistance by Leishmania infantum on mouse and sand fly infection.** *Int J Parasitol Drugs Drug Resist*. 2020. [PMID 32388220](https://pubmed.ncbi.nlm.nih.gov/32388220/)
- Palić S, Kip AE, Beijnen JH, et al. **Characterizing the non-linear pharmacokinetics of miltefosine in paediatric visceral leishmaniasis patients from Eastern Africa.** *J Antimicrob Chemother*. 2020. [PMID 32780098](https://pubmed.ncbi.nlm.nih.gov/32780098/)
- Mbui J, Olobo J, Omollo R, et al. **Pharmacokinetics, Safety, and Efficacy of an Allometric Miltefosine Regimen for the Treatment of Visceral Leishmaniasis in Eastern African Children: An Open-label, Phase II Clinical Trial.** *Clin Infect Dis*. 2019. [PMID 30188978](https://pubmed.ncbi.nlm.nih.gov/30188978/)
- Veronica J, Chandrasekaran S, Dayakar A, et al. **Iron superoxide dismutase contributes to miltefosine resistance in Leishmania donovani.** *FEBS J*. 2019. [PMID 31087522](https://pubmed.ncbi.nlm.nih.gov/31087522/)
- Alves F, Bilbe G, Blesson S, et al. **Recent Development of Visceral Leishmaniasis Treatments: Successes, Pitfalls, and Perspectives.** *Clin Microbiol Rev*. 2018. [PMID 30158301](https://pubmed.ncbi.nlm.nih.gov/30158301/)
- Sundar S, Agarwal D. **Visceral Leishmaniasis-Optimum Treatment Options in Children.** *Pediatr Infect Dis J*. 2018. [PMID 29280784](https://pubmed.ncbi.nlm.nih.gov/29280784/)
- Villa-Pulgarín JA, Gajate C, Botet J, et al. **Mitochondria and lipid raft-located FOF1-ATP synthase as major therapeutic targets in the antileishmanial and anticancer activities of ether lipid edelfosine.** *PLoS Negl Trop Dis*. 2017. [PMID 28829771](https://pubmed.ncbi.nlm.nih.gov/28829771/)
- Monge-Maillo B, López-Vélez R. **Treatment Options for Visceral Leishmaniasis and HIV Coinfection.** *AIDS Rev*. 2016. [PMID 26936761](https://pubmed.ncbi.nlm.nih.gov/26936761/)
- Wasunna M, Njenga S, Balasegaram M, et al. **Efficacy and Safety of AmBisome in Combination with Sodium Stibogluconate or Miltefosine and Miltefosine Monotherapy for African Visceral Leishmaniasis: Phase II Randomized Trial.** *PLoS Negl Trop Dis*. 2016. [PMID 27627654](https://pubmed.ncbi.nlm.nih.gov/27627654/)
- Dorlo TP, Rijal S, Ostyn B, et al. **Failure of miltefosine in visceral leishmaniasis is associated with low drug exposure.** *J Infect Dis*. 2014. [PMID 24443541](https://pubmed.ncbi.nlm.nih.gov/24443541/)
- Ostyn B, Hasker E, Dorlo TP, et al. **Failure of miltefosine treatment for visceral leishmaniasis in children and men in South-East Asia.** *PLoS One*. 2014. [PMID 24941345](https://pubmed.ncbi.nlm.nih.gov/24941345/)
- Dorlo TP, Huitema AD, Beijnen JH, et al. **Optimal dosing of miltefosine in children and adults with visceral leishmaniasis.** *Antimicrob Agents Chemother*. 2012. [PMID 22585212](https://pubmed.ncbi.nlm.nih.gov/22585212/)
- Pérez-Victoria JM, Bavchvarov BI, Torrecillas IR, et al. **Sitamaquine overcomes ABC-mediated resistance to miltefosine and antimony in Leishmania.** *Antimicrob Agents Chemother*. 2011. [PMID 21646479](https://pubmed.ncbi.nlm.nih.gov/21646479/)
- Pérez-Victoria FJ, Sánchez-Cañete MP, Seifert K, et al. **Mechanisms of experimental resistance of Leishmania to miltefosine: Implications for clinical use.** *Drug Resist Updat*. 2006. [PMID 16814199](https://pubmed.ncbi.nlm.nih.gov/16814199/)
- Bhattacharya SK, Jha TK, Sundar S, et al. **Efficacy and tolerability of miltefosine for childhood visceral leishmaniasis in India.** *Clin Infect Dis*. 2004. [PMID 14699453](https://pubmed.ncbi.nlm.nih.gov/14699453/)
- Sundar S, Rosenkaimer F, Makharia MK, et al. **Trial of oral miltefosine for visceral leishmaniasis.** *Lancet*. 1998. [PMID 9851383](https://pubmed.ncbi.nlm.nih.gov/9851383/)

## 10. 파로모마이신: 효능, 약동학, 이독성 / Paromomycin: efficacy, PK and ototoxicity

Slow pinocytic entry is why the course is 21 days; near-irreversible cochlear accumulation is the cost.

- Wu T, Zhao Z, Wang P, et al. **Drug-induced hearing loss: a real-world pharmacovigilance study using the FDA adverse event reporting system database.** *Hear Res*. 2025. [PMID 40188564](https://pubmed.ncbi.nlm.nih.gov/40188564/)
- Van Bocxlaer K, Dixon J, Platteeuw JJ, et al. **Efficacy of oleylphosphocholine in experimental cutaneous leishmaniasis.** *J Antimicrob Chemother*. 2023. [PMID 37229566](https://pubmed.ncbi.nlm.nih.gov/37229566/)
- Verrest L, Wasunna M, Kokwaro G, et al. **Geographical Variability in Paromomycin Pharmacokinetics Does Not Explain Efficacy Differences between Eastern African and Indian Visceral Leishmaniasis Patients.** *Clin Pharmacokinet*. 2021. [PMID 34105063](https://pubmed.ncbi.nlm.nih.gov/34105063/)
- Kasabalis D, Chatzis MK, Apostolidis K, et al. **Evaluation of nephrotoxicity and ototoxicity of aminosidine (paromomycin)-allopurinol combination in dogs with leishmaniosis due to Leishmania infantum: A randomized, blinded, controlled study.** *Exp Parasitol*. 2019. [PMID 31539540](https://pubmed.ncbi.nlm.nih.gov/31539540/)
- Matsushita T, Chen W, Juskeviciene R, et al. **Influence of 4'-O-Glycoside Constitution and Configuration on Ribosomal Selectivity of Paromomycin.** *J Am Chem Soc*. 2015. [PMID 26024064](https://pubmed.ncbi.nlm.nih.gov/26024064/)
- Fortin A, Hendrickx S, Yardley V, et al. **Efficacy and tolerability of oleylphosphocholine (OlPC) in a laboratory model of visceral leishmaniasis.** *J Antimicrob Chemother*. 2012. [PMID 22782488](https://pubmed.ncbi.nlm.nih.gov/22782488/)
- Musa AM, Younis B, Fadlalla A, et al. **Paromomycin for the treatment of visceral leishmaniasis in Sudan: a randomized, open-label, dose-finding study.** *PLoS Negl Trop Dis*. 2010. [PMID 21049063](https://pubmed.ncbi.nlm.nih.gov/21049063/)
- Sundar S, Chakravarty J. **Paromomycin in the treatment of leishmaniasis.** *Expert Opin Investig Drugs*. 2008. [PMID 18447603](https://pubmed.ncbi.nlm.nih.gov/18447603/)
- Sundar S, Jha TK, Thakur CP, et al. **Injectable paromomycin for Visceral leishmaniasis in India.** *N Engl J Med*. 2007. [PMID 17582067](https://pubmed.ncbi.nlm.nih.gov/17582067/)

## 11. 5가 안티몬: 활성화, 내성, 심독성 / Pentavalent antimony: activation, resistance and cardiotoxicity

Sb(V) as a prodrug, thiol-dependent efflux as the resistance mechanism (claim 5), and the QTc/pancreatic ceiling on dose escalation.

- Kumari S, Kumar V, Tiwari RK, et al. **Amphotericin B: A drug of choice for Visceral Leishmaniasis.** *Acta Trop*. 2022. [PMID 35998680](https://pubmed.ncbi.nlm.nih.gov/35998680/)
- Singh N, Kumar M, Singh RK. **Leishmaniasis: current status of available drugs and new potential drug targets.** *Asian Pac J Trop Med*. 2012. [PMID 22575984](https://pubmed.ncbi.nlm.nih.gov/22575984/)
- Mohamed-Ahmed AH, Brocchini S, Croft SL. **Recent advances in development of amphotericin B formulations for the treatment of visceral leishmaniasis.** *Curr Opin Infect Dis*. 2012. [PMID 23147810](https://pubmed.ncbi.nlm.nih.gov/23147810/)
- Sinha S, Sundaram S, Kumar V, et al. **Antimony resistance during Visceral Leishmaniasis: A possible consequence of serial mutations in ABC transporters of Leishmania species.** *Bioinformation*. 2011. [PMID 21584185](https://pubmed.ncbi.nlm.nih.gov/21584185/)
- Mukhopadhyay R, Mukherjee S, Mukherjee B, et al. **Characterisation of antimony-resistant Leishmania donovani isolates: biochemical and biophysical studies and interaction with host cells.** *Int J Parasitol*. 2011. [PMID 21920365](https://pubmed.ncbi.nlm.nih.gov/21920365/)
- Palumbo E. **Oral miltefosine treatment in children with visceral leishmaniasis: a brief review.** *Braz J Infect Dis*. 2008. [PMID 18553005](https://pubmed.ncbi.nlm.nih.gov/18553005/)
- Samant M, Sahasrabuddhe AA, Singh N, et al. **Proteophosphoglycan is differentially expressed in sodium stibogluconate-sensitive and resistant Indian clinical isolates of Leishmania donovani.** *Parasitology*. 2007. [PMID 17362540](https://pubmed.ncbi.nlm.nih.gov/17362540/)
- Jha TK. **Drug unresponsiveness & combination therapy for kala-azar.** *Indian J Med Res*. 2006. [PMID 16778318](https://pubmed.ncbi.nlm.nih.gov/16778318/)
- Sinha PK, Ranjan A, Singh VP, et al. **Visceral leishmaniasis (kala-azar)--the Bihar (India) perspective.** *J Infect*. 2006. [PMID 16269185](https://pubmed.ncbi.nlm.nih.gov/16269185/)
- Alvarez M, Malécot CO, Gannier F, et al. **Antimony-induced cardiomyopathy in guinea-pig and protection by L-carnitine.** *Br J Pharmacol*. 2005. [PMID 15644865](https://pubmed.ncbi.nlm.nih.gov/15644865/)
- Das VN, Ranjan A, Bimal S, et al. **Magnitude of unresponsiveness to sodium stibogluconate in the treatment of visceral leishmaniasis in Bihar.** *Natl Med J India*. 2005. [PMID 16130613](https://pubmed.ncbi.nlm.nih.gov/16130613/)
- Sundar S, Rai M. **Advances in the treatment of leishmaniasis.** *Curr Opin Infect Dis*. 2002. [PMID 12821836](https://pubmed.ncbi.nlm.nih.gov/12821836/)
- Croft SL, Yardley V. **Chemotherapy of leishmaniasis.** *Curr Pharm Des*. 2002. [PMID 11860369](https://pubmed.ncbi.nlm.nih.gov/11860369/)
- Berhe N, Abraham Y, Hailu A, et al. **Electrocardiographic findings in Ethiopians on pentavalent antimony therapy for visceral leishmaniasis.** *East Afr Med J*. 2001. [PMID 12219968](https://pubmed.ncbi.nlm.nih.gov/12219968/)
- Ortega-Carnicer J, Alcázar R, De la Torre M, et al. **Pentavalent antimonial-induced torsade de pointes.** *J Electrocardiol*. 1997. [PMID 9141610](https://pubmed.ncbi.nlm.nih.gov/9141610/)
- Chulay JD, Fleckenstein L, Smith DH. **Pharmacokinetics of antimony during treatment of visceral leishmaniasis with sodium stibogluconate or meglumine antimoniate.** *Trans R Soc Trop Med Hyg*. 1988. [PMID 2845611](https://pubmed.ncbi.nlm.nih.gov/2845611/)
- Chulay JD, Spencer HC, Mugambi M. **Electrocardiographic changes during treatment of leishmaniasis with pentavalent antimony (sodium stibogluconate).** *Am J Trop Med Hyg*. 1985. [PMID 2992303](https://pubmed.ncbi.nlm.nih.gov/2992303/)

## 12. 병용요법과 단기요법 / Combination and short-course therapy

The trials behind claim 6: synergy that lives in time rather than in concentration.

- Sundar S, Singh A. **Chemotherapeutics of visceral leishmaniasis: present and future developments.** *Parasitology*. 2018. [PMID 29215329](https://pubmed.ncbi.nlm.nih.gov/29215329/)
- Sundar S, Chakravarty J. **An update on pharmacotherapy for leishmaniasis.** *Expert Opin Pharmacother*. 2015. [PMID 25346016](https://pubmed.ncbi.nlm.nih.gov/25346016/)
- Sundar S, Chakravarty J. **Leishmaniasis: an update of current pharmacotherapy.** *Expert Opin Pharmacother*. 2013. [PMID 23256501](https://pubmed.ncbi.nlm.nih.gov/23256501/)
- van Griensven J, Diro E. **Visceral leishmaniasis.** *Infect Dis Clin North Am*. 2012. [PMID 22632641](https://pubmed.ncbi.nlm.nih.gov/22632641/)
- van Griensven J, Boelaert M. **Combination therapy for visceral leishmaniasis.** *Lancet*. 2011. [PMID 21255829](https://pubmed.ncbi.nlm.nih.gov/21255829/)
- Matlashewski G, Arana B, Kroeger A, et al. **Visceral leishmaniasis: elimination with existing interventions.** *Lancet Infect Dis*. 2011. [PMID 21453873](https://pubmed.ncbi.nlm.nih.gov/21453873/)
- Olliaro PL. **Drug combinations for visceral leishmaniasis.** *Curr Opin Infect Dis*. 2010. [PMID 20871400](https://pubmed.ncbi.nlm.nih.gov/20871400/)

## 13. HIV 동반감염과 재발 / HIV co-infection and relapse

The separatrix argument of claim 3: the same residual burden, a different finish line.

- Kennedy GC, O'Brien K, Nyakundi H, et al. **Visceral leishmaniasis follow-up and treatment outcomes in Tiaty East and West sub-counties, Kenya: Cure, relapse, and Post Kala-azar Dermal Leishmaniasis.** *PLoS One*. 2024. [PMID 38917127](https://pubmed.ncbi.nlm.nih.gov/38917127/)
- Goyal V, Das VNR, Singh SN, et al. **Long-term incidence of relapse and post-kala-azar dermal leishmaniasis after three different visceral leishmaniasis treatment regimens in Bihar, India.** *PLoS Negl Trop Dis*. 2020. [PMID 32687498](https://pubmed.ncbi.nlm.nih.gov/32687498/)
- Coutinho JVSC, Santos FSD, Ribeiro RDSP, et al. **Visceral leishmaniasis and leishmaniasis-HIV coinfection: comparative study.** *Rev Soc Bras Med Trop*. 2017. [PMID 29160515](https://pubmed.ncbi.nlm.nih.gov/29160515/)
- Lindoso JA, Cunha MA, Queiroz IT, et al. **Leishmaniasis-HIV coinfection: current challenges.** *HIV AIDS (Auckl)*. 2016. [PMID 27785103](https://pubmed.ncbi.nlm.nih.gov/27785103/)
- Diro E, Lynen L, Ritmeijer K, et al. **Visceral Leishmaniasis and HIV coinfection in East Africa.** *PLoS Negl Trop Dis*. 2014. [PMID 24968313](https://pubmed.ncbi.nlm.nih.gov/24968313/)
- Lindoso JA, Cota GF, da Cruz AM, et al. **Visceral leishmaniasis and HIV coinfection in Latin America.** *PLoS Negl Trop Dis*. 2014. [PMID 25233461](https://pubmed.ncbi.nlm.nih.gov/25233461/)
- Monge-Maillo B, Norman FF, Cruz I, et al. **Visceral leishmaniasis and HIV coinfection in the Mediterranean region.** *PLoS Negl Trop Dis*. 2014. [PMID 25144380](https://pubmed.ncbi.nlm.nih.gov/25144380/)
- van Griensven J, Diro E, Lopez-Velez R, et al. **HIV-1 protease inhibitors for treatment of visceral leishmaniasis in HIV-co-infected individuals.** *Lancet Infect Dis*. 2013. [PMID 23427890](https://pubmed.ncbi.nlm.nih.gov/23427890/)
- Pandey K, Pun SB, Pandey BD. **Relapse of kala-azar after use of multiple drugs: a case report and brief review of literature.** *Indian J Med Microbiol*. 2012. [PMID 22664444](https://pubmed.ncbi.nlm.nih.gov/22664444/)
- Cota GF, de Sousa MR, Rabello A. **Predictors of visceral leishmaniasis relapse in HIV-infected patients: a systematic review.** *PLoS Negl Trop Dis*. 2011. [PMID 21666786](https://pubmed.ncbi.nlm.nih.gov/21666786/)
- Murray HW. **Kala-azar as an AIDS-related opportunistic infection.** *AIDS Patient Care STDS*. 1999. [PMID 10800524](https://pubmed.ncbi.nlm.nih.gov/10800524/)
- Grevelink SA, Lerner EA. **Leishmaniasis.** *J Am Acad Dermatol*. 1996. [PMID 8642091](https://pubmed.ncbi.nlm.nih.gov/8642091/)
- MORTON TC. **Kala-azar; relapse following splenectomy.** *Trans R Soc Trop Med Hyg*. 1949. [PMID 18139102](https://pubmed.ncbi.nlm.nih.gov/18139102/)

## 14. PKDL / Post-kala-azar dermal leishmaniasis

A lesion made by returning immunity, in the organ with the lowest drug exposure.

- Sundar S, Chakravarty J, Singh J, et al. **Amphotericin B Deoxycholate Treatment of Post-Kala-Azar Dermal Leishmaniasis in India.** *Am J Trop Med Hyg*. 2024. [PMID 38834057](https://pubmed.ncbi.nlm.nih.gov/38834057/)
- Roy S, Moulik S, Roy M, et al. **Monitoring the Long-Term Effectiveness of Miltefosine in Indian Post-Kala-Azar Dermal Leishmaniasis.** *Am J Trop Med Hyg*. 2024. [PMID 38442428](https://pubmed.ncbi.nlm.nih.gov/38442428/)
- Sundar S, Singh J, Dinkar A, et al. **Safety and Effectiveness of Miltefosine in Post-Kala-Azar Dermal Leishmaniasis: An Observational Study.** *Open Forum Infect Dis*. 2023. [PMID 37234513](https://pubmed.ncbi.nlm.nih.gov/37234513/)
- Monge-Maillo B, Norman FF, Chamorro-Tojeiro S, et al. **Post-kala-azar dermal leishmaniasis due to Leishmania infantum in an HIV-negative patient treated with miltefosine.** *J Travel Med*. 2022. [PMID 34668558](https://pubmed.ncbi.nlm.nih.gov/34668558/)
- Pijpers J, den Boer ML, Essink DR, et al. **The safety and efficacy of miltefosine in the long-term treatment of post-kala-azar dermal leishmaniasis in South Asia - A review and meta-analysis.** *PLoS Negl Trop Dis*. 2019. [PMID 30742620](https://pubmed.ncbi.nlm.nih.gov/30742620/)
- Pradhan A, Basak S, Chowdhury T, et al. **Keratitis After Post-Kala-Azar Dermal Leishmaniasis.** *Cornea*. 2018. [PMID 29040115](https://pubmed.ncbi.nlm.nih.gov/29040115/)
- Zijlstra EE, Musa AM, Khalil EA, et al. **Post-kala-azar dermal leishmaniasis.** *Lancet Infect Dis*. 2003. [PMID 12560194](https://pubmed.ncbi.nlm.nih.gov/12560194/)

## 15. 영양·숙주 인자 / Nutrition and host factors

Malnutrition as a modifier of T-cell priming (the MALNUT parameter), and treatment in pregnancy.

- Karampas G, Koulouraki S, Daikos GL, et al. **Visceral Leishmaniasis in a Twin Pregnancy: A Case Report and Review of the Literature.** *J Clin Med*. 2024. [PMID 38673673](https://pubmed.ncbi.nlm.nih.gov/38673673/)
- Dahal P, Singh-Phulgenda S, Maguire BJ, et al. **Visceral Leishmaniasis in pregnancy and vertical transmission: A systematic literature review on the therapeutic orphans.** *PLoS Negl Trop Dis*. 2021. [PMID 34375339](https://pubmed.ncbi.nlm.nih.gov/34375339/)
- Custodio E, López-Alcalde J, Herrero M, et al. **Nutritional supplements for patients being treated for active visceral leishmaniasis.** *Cochrane Database Syst Rev*. 2018. [PMID 29578237](https://pubmed.ncbi.nlm.nih.gov/29578237/)
- Atia AM, Mumina A, Tayler-Smith K, et al. **Sodium stibogluconate and paromomycin for treating visceral leishmaniasis under routine conditions in eastern Sudan.** *Trop Med Int Health*. 2015. [PMID 26427033](https://pubmed.ncbi.nlm.nih.gov/26427033/)
- Oliveira JM, Fernandes AC, Dorval ME, et al. **[Mortality due to visceral leishmaniasis: clinical and laboratory characteristics].** *Rev Soc Bras Med Trop*. 2010. [PMID 20464151](https://pubmed.ncbi.nlm.nih.gov/20464151/)
- Figueiró-Filho EA, Duarte G, El-Beitune P, et al. **Visceral leishmaniasis (kala-azar) and pregnancy.** *Infect Dis Obstet Gynecol*. 2004. [PMID 15460194](https://pubmed.ncbi.nlm.nih.gov/15460194/)
- Caldas AJ, Costa JM, Silva AA, et al. **Risk factors associated with asymptomatic infection by Leishmania chagasi in north-east Brazil.** *Trans R Soc Trop Med Hyg*. 2002. [PMID 11925984](https://pubmed.ncbi.nlm.nih.gov/11925984/)
- Rosenblatt JE. **Antiparasitic agents.** *Mayo Clin Proc*. 1999. [PMID 10560606](https://pubmed.ncbi.nlm.nih.gov/10560606/)

## 16. 약동-약력학 모델링과 QSP 방법론 / PK-PD modelling and QSP methodology

Prior exposure-response work in leishmaniasis and the animal models the kill-rate parameters are anchored to.

- Alexander N, Giraldo-Parra L, Rebellón-Sánchez DE, et al. **Modelling immune gene expression profiles as pharmacodynamic endpoints of antileishmanial treatment.** *Br J Clin Pharmacol*. 2026. [PMID 42397049](https://pubmed.ncbi.nlm.nih.gov/42397049/)
- Delorenzi JC, Melo TS, Silveira FT, et al. **Oral 18-methoxycoronaridine activity in simian and murine Leishmania amazonensis infection.** *Front Pharmacol*. 2026. [PMID 42038284](https://pubmed.ncbi.nlm.nih.gov/42038284/)
- Roy M, Sarkar D, Chatterjee M. **Quantitative monitoring of experimental and human leishmaniasis employing amastigote-specific genes.** *Parasitology*. 2022. [PMID 35535469](https://pubmed.ncbi.nlm.nih.gov/35535469/)
- Klionsky DJ, Abdel-Aziz AK, Abdelfatah S, et al. **Guidelines for the use and interpretation of assays for monitoring autophagy (4th edition)(1).** *Autophagy*. 2021. [PMID 33634751](https://pubmed.ncbi.nlm.nih.gov/33634751/)
- Das S, Banerjee A, Kamran M, et al. **A chemical inhibitor of heat shock protein 78 (HSP78) from Leishmania donovani represents a potential antileishmanial drug candidate.** *J Biol Chem*. 2020. [PMID 32471865](https://pubmed.ncbi.nlm.nih.gov/32471865/)
- Tran HTL, Morikawa K, Anggakusuma, et al. **OCIAD1 is a host mitochondrial substrate of the hepatitis C virus NS3-4A protease.** *PLoS One*. 2020. [PMID 32697788](https://pubmed.ncbi.nlm.nih.gov/32697788/)
- Kip AE, Castro MDM, Gomez MA, et al. **Simultaneous population pharmacokinetic modelling of plasma and intracellular PBMC miltefosine concentrations in New World cutaneous leishmaniasis and exploration of exposure-response relationships.** *J Antimicrob Chemother*. 2018. [PMID 29757380](https://pubmed.ncbi.nlm.nih.gov/29757380/)
- Vishwakarma P, Parmar N, Yadav PK, et al. **15d-Prostaglandin J2 induced reactive oxygen species-mediated apoptosis during experimental visceral leishmaniasis.** *J Mol Med (Berl)*. 2016. [PMID 26830627](https://pubmed.ncbi.nlm.nih.gov/26830627/)
- Bhuniya D, Mukkavilli R, Shivahare R, et al. **Aminothiazoles: Hit to lead development to identify antileishmanial agents.** *Eur J Med Chem*. 2015. [PMID 26318065](https://pubmed.ncbi.nlm.nih.gov/26318065/)
- Singh MK, Paul J, De T, et al. **Bioactivity guided fractionation of Moringa oleifera Lam. flower targeting Leishmania donovani.** *Indian J Exp Biol*. 2015. [PMID 26669018](https://pubmed.ncbi.nlm.nih.gov/26669018/)
- Tomioka H. **Editorial: Current status and perspective on drug targets in tubercle bacilli and drug design of antituberculous agents based on structure-activity relationship.** *Curr Pharm Des*. 2014. [PMID 24245755](https://pubmed.ncbi.nlm.nih.gov/24245755/)
- Datta N, Mukherjee S, Das L, et al. **Targeting of immunostimulatory DNA cures experimental visceral leishmaniasis through nitric oxide up-regulation and T cell activation.** *Eur J Immunol*. 2003. [PMID 12778468](https://pubmed.ncbi.nlm.nih.gov/12778468/)
- Mullen AB, Carter KC, Baillie AJ. **Comparison of the efficacies of various formulations of amphotericin B against murine visceral leishmaniasis.** *Antimicrob Agents Chemother*. 1997. [PMID 9333030](https://pubmed.ncbi.nlm.nih.gov/9333030/)

## 17. 신약 파이프라인과 면역치료 / Pipeline compounds and immunotherapy

What an oral, short-course, all-ages antileishmanial would have to beat, and the therapeutic-vaccine route into the same host arm.

- Aronson NE, Musa AM, Satoskar AR. **Leishmaniasis.** *N Engl J Med*. 2026. [PMID 42202321](https://pubmed.ncbi.nlm.nih.gov/42202321/)
- Gupta J, Menon Y, Kumar S, et al. **Vaccine Designing Technology against Leishmaniasis: Current Challenges and Implication.** *Curr Drug Discov Technol*. 2025. [PMID 38798212](https://pubmed.ncbi.nlm.nih.gov/38798212/)
- Sundar S, Singh VK, Agrawal N, et al. **Investigational new drugs for the treatment of leishmaniasis.** *Expert Opin Investig Drugs*. 2024. [PMID 39225742](https://pubmed.ncbi.nlm.nih.gov/39225742/)
- Kancharla P, Ortiz D, Fargo CM, et al. **Discovery and Optimization of Tambjamines as a Novel Class of Antileishmanial Agents.** *J Med Chem*. 2024. [PMID 38722757](https://pubmed.ncbi.nlm.nih.gov/38722757/)
- Kelleci K, Allahverdiyev A, Bağirova M, et al. **Particulate and non-particle adjuvants in Leishmaniasis vaccine designs: A review.** *J Vector Borne Dis*. 2023. [PMID 37417162](https://pubmed.ncbi.nlm.nih.gov/37417162/)
- de Santiago-Silva KM, da Silva Gomes GF, Perez CC, et al. **Molecular Targets for Chalcones in Antileishmanial Drug Discovery.** *Mini Rev Med Chem*. 2023. [PMID 36705240](https://pubmed.ncbi.nlm.nih.gov/36705240/)
- Kaye PM, Matlashewski G, Mohan S, et al. **Vaccine value profile for leishmaniasis.** *Vaccine*. 2023. [PMID 37951693](https://pubmed.ncbi.nlm.nih.gov/37951693/)
- Shah AP, Hura N, Kishore Babu N, et al. **A Core-Linker-Polyamine (CLP) Strategy Enables Rapid Discovery of Antileishmanial Aminoalkylquinolinecarboxamides That Target Oxidative Stress Mechanism.** *ChemMedChem*. 2022. [PMID 35638162](https://pubmed.ncbi.nlm.nih.gov/35638162/)
- Imran M, Khan SA, Abida, et al. **Small molecules as kinetoplastid specific proteasome inhibitors for leishmaniasis: a patent review from 1998 to 2021.** *Expert Opin Ther Pat*. 2022. [PMID 35220857](https://pubmed.ncbi.nlm.nih.gov/35220857/)
- Mantravadi PK, Parthasarathy A, Kalesh K. **Antileishmanial Drug Development: A Review of Modern Molecular Chemical Tools and Research Strategies.** *Curr Med Chem*. 2021. [PMID 33238841](https://pubmed.ncbi.nlm.nih.gov/33238841/)
- Nagle A, Biggart A, Be C, et al. **Discovery and Characterization of Clinical Candidate LXE408 as a Kinetoplastid-Selective Proteasome Inhibitor for the Treatment of Leishmaniases.** *J Med Chem*. 2020. [PMID 32667203](https://pubmed.ncbi.nlm.nih.gov/32667203/)
- Deeks ED. **Fexinidazole: First Global Approval.** *Drugs*. 2019. [PMID 30635838](https://pubmed.ncbi.nlm.nih.gov/30635838/)
- Ghorbani M, Farhoudi R. **Leishmaniasis in humans: drug or vaccine therapy?.** *Drug Des Devel Ther*. 2018. [PMID 29317800](https://pubmed.ncbi.nlm.nih.gov/29317800/)
- Gillespie PM, Beaumier CM, Strych U, et al. **Status of vaccine research and development of vaccines for leishmaniasis.** *Vaccine*. 2016. [PMID 26973063](https://pubmed.ncbi.nlm.nih.gov/26973063/)
- Barrett MP, Croft SL. **Management of trypanosomiasis and leishmaniasis.** *Br Med Bull*. 2012. [PMID 23137768](https://pubmed.ncbi.nlm.nih.gov/23137768/)
- Wyllie S, Patterson S, Stojanovski L, et al. **The anti-trypanosome drug fexinidazole shows potential for treating visceral leishmaniasis.** *Sci Transl Med*. 2012. [PMID 22301556](https://pubmed.ncbi.nlm.nih.gov/22301556/)
- Croft SL, Engel J. **Miltefosine--discovery of the antileishmanial activity of phospholipid derivatives.** *Trans R Soc Trop Med Hyg*. 2006. [PMID 16904717](https://pubmed.ncbi.nlm.nih.gov/16904717/)
- Jha TK, Sundar S, Thakur CP, et al. **A phase II dose-ranging study of sitamaquine for the treatment of visceral leishmaniasis in India.** *Am J Trop Med Hyg*. 2005. [PMID 16354802](https://pubmed.ncbi.nlm.nih.gov/16354802/)
- Wasunna MK, Rashid JR, Mbui J, et al. **A phase II dose-increasing study of sitamaquine for the treatment of visceral leishmaniasis in Kenya.** *Am J Trop Med Hyg*. 2005. [PMID 16282296](https://pubmed.ncbi.nlm.nih.gov/16282296/)
- Sangraula H, Sharma KK, Rijal S, et al. **Orally effective drugs for kala-azar (visceral leishmaniasis): focus on miltefosine and sitamaquine.** *J Assoc Physicians India*. 2003. [PMID 14621038](https://pubmed.ncbi.nlm.nih.gov/14621038/)

---

## 지침·1차 자료 (Guidelines and primary sources, consulted directly)

- WHO. *Control of the leishmaniases.* WHO Technical Report Series 949 — the
  source of the regimen definitions used as scenarios S02–S20.
- WHO. *Guideline for the treatment of visceral leishmaniasis in HIV
  co-infected patients in East Africa and South-East Asia* (2022) — the source
  of the L-AmB 30–40 mg/kg total plus miltefosine 28-day combination modelled
  in S18, and of the secondary-prophylaxis logic in the map.
- WHO. *Regional Strategic Framework for accelerating and sustaining
  elimination of kala-azar in the South-East Asia Region.*
- DNDi target product profiles for visceral leishmaniasis — the specification
  an oral, short-course, all-ages replacement would have to meet, which is the
  yardstick used in the pipeline cluster of the map.

## 참고문헌 총계 (Reference count)

- PubMed 레코드: **223편** (전부 직접 조회, 중복 제거)
- 분류 섹션: **17개**
