#!/usr/bin/env python3
"""
mpm_fetch_references.py -- build mpm_references.md from PubMed records that were
ACTUALLY RETURNED by NCBI E-utilities.

Nothing is written from memory: every entry is a SEARCH, and only the record
NCBI hands back (its real PMID, title, journal, year, first author) is written
into the reference file.  Searches that return nothing are reported as misses
at the end of the run rather than filled in from recollection.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                print("   ! %s" % e, file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))


def esearch(term, n=1):
    u = (EUTILS + "esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d&term=%s"
         % (n, urllib.parse.quote(term)))
    d = _get(u)
    if not d:
        return []
    return d.get("esearchresult", {}).get("idlist", [])


def esummary(pmids):
    if not pmids:
        return {}
    u = EUTILS + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(pmids)
    d = _get(u)
    if not d:
        return {}
    return d.get("result", {})


# (section key, search term, required substring in the title, "" = no check)
QUERIES = [
    # ---------------- 1. epidemiology, exposure, burden ----------------
    ("epi", "Robinson Lake advances malignant mesothelioma New England Journal Medicine", "mesothelioma"),
    ("epi", "Carbone mesothelioma scientific clues prevention diagnosis therapy CA Cancer J Clin", "mesothelioma"),
    ("epi", "Mesothelioma Nature Reviews Disease Primers Carbone Adusumilli", "esotheli"),
    ("epi", "Peto European mesothelioma epidemic projections", "mesothelioma"),
    ("epi", "Odgerel Global epidemiology mesothelioma incidence estimation", "mesothelioma"),
    ("epi", "Furuya Global asbestos disaster international journal environmental research public health", "sbestos"),
    ("epi", "Delgermaa global mesothelioma deaths reported World Health Organization", "esotheli"),
    ("epi", "Henley mesothelioma incidence United States trends", "esotheli"),
    ("epi", "Bianchi malignant mesothelioma global incidence", "esotheli"),
    ("epi", "Marinaccio pleural mesothelioma occupational exposure Italian national registry", "esotheli"),
    ("epi", "Lacourt occupational asbestos exposure risk pleural mesothelioma time since first exposure", "esotheli"),
    ("epi", "Frost latency period mesothelioma British asbestos workers survey", "atenc"),
    ("epi", "Baumann Ambrosi Carbone mesothelioma erionite Cappadocia Turkey", "erionite"),
    ("epi", "Carbone erionite exposure North Dakota mesothelioma risk", "erionite"),
    ("epi", "Goldberg Luce health impact non-occupational asbestos exposure", "sbestos"),
    ("epi", "Reid household domestic asbestos exposure mesothelioma risk women", "sbestos"),
    ("epi", "Mensi mesothelioma women non-occupational exposure", "esotheli"),
    ("epi", "Teta radiation induced mesothelioma second malignancy Hodgkin", "esotheli"),
    ("epi", "peritoneal mesothelioma incidence survival SEER analysis", "esotheli"),
    ("epi", "van Gerwen Taioli mesothelioma survival disparities United States", "esotheli"),

    # ---------------- 2. fibre biology, biopersistence, frustrated phagocytosis ----------------
    ("fibre", "Stanton relation particle dimension carcinogenicity mesothelioma pleural", "arcinogen"),
    ("fibre", "Donaldson Poland Duffin asbestos fibre length pleural inflammation paradigm", "fibre"),
    ("fibre", "Boulanger Quantification short long asbestos fibers lung mesothelioma", "asbestos"),
    ("fibre", "Bernstein biopersistence chrysotile amphibole asbestos inhalation", "biopersist"),
    ("fibre", "Mossman Lippmann asbestos fiber type mesothelioma pathogenesis mechanisms", "sbestos"),
    ("fibre", "Broaddus asbestos pleural mesothelial cell apoptosis reactive oxygen species", "sbestos"),
    ("fibre", "Miserocchi translocation particles parietal pleura lymphatic stomata", "pleur"),
    ("fibre", "Boutin black spots parietal pleura asbestos fibers thoracoscopy", "pleura"),
    ("fibre", "Donaldson pleural translocation nanoparticles fibres mesothelium", "pleur"),
    ("fibre", "Dostert Petrilli innate immune activation NALP3 inflammasome asbestos silica", "sbestos"),
    ("fibre", "Yang Rivera Carbone TNF-alpha inhibits asbestos-induced cytotoxicity mesothelial", "sbestos"),
    ("fibre", "Yang Carbone programmed necrosis HMGB1 asbestos mesothelioma", "HMGB1"),
    ("fibre", "Jube Carbone HMGB1 mesothelioma biomarker therapeutic target", "HMGB1"),
    ("fibre", "Napolitano Carbone HMGB1 secretion asbestos mesothelioma prevention", "HMGB1"),
    ("fibre", "Hillegass Mossman inflammasome mesothelioma pleural inflammation IL-1", "flammasome"),
    ("fibre", "Kadariya Testa inflammation-related IL-1 receptor antagonist anakinra mesothelioma", "esotheli"),
    ("fibre", "Chew Toyokuni iron asbestos carcinogenesis oxidative stress mesothelioma", "iron"),
    ("fibre", "Toyokuni iron overload asbestos-induced mesothelial carcinogenesis", "iron"),
    ("fibre", "Wang Xu frustrated phagocytosis macrophage asbestos fiber", "phagocytosis"),
    ("fibre", "Sayan Mossman NLRP3 inflammasome mesothelioma asbestos review", "flammasome"),

    # ---------------- 3. genomics and tumour suppressors ----------------
    ("gene", "Bueno comprehensive genomic analysis malignant pleural mesothelioma Nature Genetics", "esotheli"),
    ("gene", "Hmeljak Cherniack Integrative molecular characterization malignant pleural mesothelioma TCGA", "esotheli"),
    ("gene", "Bott Ladanyi NF2 CUL1 mesothelioma genome sequencing BAP1", "esotheli"),
    ("gene", "Testa Carbone germline BAP1 mutations predispose mesothelioma", "BAP1"),
    ("gene", "Bononi Carbone BAP1 calcium mitochondria apoptosis tumour suppression", "BAP1"),
    ("gene", "Carbone BAP1 cancer syndrome tumor predisposition Nature Reviews Cancer", "BAP1"),
    ("gene", "Nasu Carbone Testa low incidence germline BAP1 mutation mesothelioma", "BAP1"),
    ("gene", "Cigognetti BAP1 immunohistochemistry distinguishes malignant mesothelioma reactive", "BAP1"),
    ("gene", "Sekido Inactivation merlin NF2 tumor suppressor mesothelioma", "NF2"),
    ("gene", "Sekido molecular pathogenesis malignant mesothelioma review", "esotheli"),
    ("gene", "Cheng p16 CDKN2A homozygous deletion malignant mesothelioma FISH", "CDKN2A"),
    ("gene", "Illei Ladanyi Homozygous deletion CDKN2A methylthioadenosine phosphorylase mesothelioma", "esotheli"),
    ("gene", "Hida Sekido MTAP immunohistochemistry surrogate CDKN2A deletion mesothelioma", "MTAP"),
    ("gene", "Chiosea Dacic p16 FISH pleural malignant mesothelioma reactive mesothelial", "p16"),
    ("gene", "Tallet TERT promoter mutations malignant pleural mesothelioma", "TERT"),
    ("gene", "Guo whole exome sequencing malignant pleural mesothelioma", "esotheli"),
    ("gene", "Quetel Jean genomic epigenomic immune subtypes malignant pleural mesothelioma", "esotheli"),
    ("gene", "Blum Jaurand Molecular classification malignant pleural mesothelioma transcriptome", "esotheli"),
    ("gene", "Alcala Fernandez-Cuesta integrative continuous molecular profiles pleural mesothelioma", "esotheli"),
    ("gene", "LaFave Wilson BAP1 loss EZH2 dependency malignant mesothelioma", "EZH2"),
    ("gene", "Mavrakis Smolen disordered methionine metabolism MTAP PRMT5 vulnerability", "PRMT5"),
    ("gene", "Kryukov MTAP deletion confers PRMT5 dependency cancer", "PRMT5"),
    ("gene", "Szlosarek arginine deprivation ASS1 argininosuccinate synthetase mesothelioma", "rginin"),

    # ---------------- 4. Hippo-YAP, mechanotransduction, cell cycle ----------------
    ("hippo", "Sekido Hippo pathway dysfunction malignant mesothelioma YAP", "Hippo"),
    ("hippo", "Yokoyama Sekido YAP1 is involved mesothelioma development", "YAP"),
    ("hippo", "Zhang Yang Hippo YAP TEAD mesothelioma therapeutic target", "YAP"),
    ("hippo", "Dupont Piccolo role YAP TAZ mechanotransduction Nature", "mechanotransduction"),
    ("hippo", "Tranchant Poulikakos Merlin NF2 FAK inhibitor sensitivity mesothelioma", "FAK"),
    ("hippo", "Shapiro Merlin deficiency predicts FAK inhibitor response mesothelioma", "FAK"),
    ("hippo", "Miyanaga Sekido Hippo pathway gene mutations malignant mesothelioma cell lines", "esotheli"),
    ("hippo", "Felley-Bosco Stahel MicroRNA and Hippo signalling mesothelioma", "esotheli"),

    # ---------------- 5. histology, EMT, grading, pathology ----------------
    ("path", "Husain Guidelines pathologic diagnosis malignant mesothelioma 2020 update consensus", "esotheli"),
    ("path", "Nicholson EURACAN IASLC proposals pathologic diagnosis mesothelioma architectural patterns", "esotheli"),
    ("path", "Sauter WHO classification tumors thoracic 2021 mesothelioma", "esotheli"),
    ("path", "Kadota nuclear grading system epithelioid malignant pleural mesothelioma", "grad"),
    ("path", "Rosen Nuclear grade mitosis necrosis score epithelioid mesothelioma validation", "esotheli"),
    ("path", "Galateau-Salle New insights diagnosis classification malignant mesothelioma", "esotheli"),
    ("path", "Meyerhoff Klabatsa impact sarcomatoid histology survival malignant pleural mesothelioma", "arcomatoid"),
    ("path", "Schunselaar Baas sarcomatoid mesothelioma epithelial mesenchymal transition", "esenchymal"),
    ("path", "Fassina epithelial mesenchymal transition malignant mesothelioma", "esenchymal"),
    ("path", "Ordonez value calretinin immunostaining diagnosis epithelioid mesothelioma", "alretinin"),
    ("path", "Chang Pastan Molecular cloning mesothelin differentiation antigen mesothelioma ovarian", "mesothelin"),
    ("path", "Hassan Pastan mesothelin tumor differentiation antigen therapy", "esothelin"),
    ("path", "Rump Pastan binding mesothelin MUC16 CA125 cell adhesion", "MUC16"),
    ("path", "Berg Churg localized malignant mesothelioma diffuse distinction", "esotheli"),
    ("path", "Klebe desmoplastic mesothelioma clinicopathological", "desmoplastic"),

    # ---------------- 6. imaging, mRECIST, tumour volume ----------------
    ("imag", "Byrne Nowak Modified RECIST criteria assessment response malignant pleural mesothelioma", "esotheli"),
    ("imag", "Armato Nowak Revised modified response evaluation criteria solid tumors mesothelioma mRECIST 1.1", "esotheli"),
    ("imag", "Armato Imaging response assessment mesothelioma tumor thickness measurement variability", "esotheli"),
    ("imag", "Nowak measurement pleural thickness computed tomography mesothelioma reproducibility", "esotheli"),
    ("imag", "Pass Temeck preoperative tumor volume prognostic malignant pleural mesothelioma", "volume"),
    ("imag", "Gill Naidich tumor volume mesothelioma prognostic computed tomography", "volume"),
    ("imag", "Rusch initial analysis volumetric staging malignant pleural mesothelioma", "volum"),
    ("imag", "Frauenfelder Weder volumetry malignant pleural mesothelioma response assessment", "volum"),
    ("imag", "Armato Observer variability mesothelioma tumor thickness CT measurements", "esotheli"),
    ("imag", "Labby Armato Optimization response classification criteria mesothelioma tumor measurement", "esotheli"),
    ("imag", "Sensakovic Armato computerized segmentation mesothelioma tumor CT", "esotheli"),
    ("imag", "Katz Volumetric CT mesothelioma response criteria comparison", "esotheli"),
    ("imag", "Tsim Blyth Diagnostic imaging pleural malignancy review", "pleural"),

    # ---------------- 7. staging and prognosis ----------------
    ("stage", "Rusch Nowak Initial analysis IASLC mesothelioma staging database", "esotheli"),
    ("stage", "Rice Nowak IASLC mesothelioma staging project eighth edition TNM", "esotheli"),
    ("stage", "Curran Sahmoud prognostic factors malignant mesothelioma EORTC", "prognostic"),
    ("stage", "Herndon Green Factors predictive survival CALGB mesothelioma trials", "urvival"),
    ("stage", "Edwards prognostic factors malignant mesothelioma validation EORTC CALGB", "prognostic"),
    ("stage", "Kao Neutrophil-to-lymphocyte ratio prognostic malignant mesothelioma", "ymphocyte"),
    ("stage", "Meniawy Nowak Existing models prognostication mesothelioma validation", "esotheli"),
    ("stage", "Brims novel clinical prediction model mesothelioma decision tree", "esotheli"),

    # ---------------- 8. biomarkers ----------------
    ("bio", "Robinson Creaney Mesothelin-family proteins diagnosis mesothelioma Lancet", "esothelin"),
    ("bio", "Creaney Robinson serum soluble mesothelin related protein tumour volume mesothelioma", "esothelin"),
    ("bio", "Hollevoet Creaney serum mesothelin diagnosis pleural mesothelioma individual patient data meta-analysis", "esothelin"),
    ("bio", "Beyer MESOMARK potential diagnostic test mesothelioma soluble mesothelin", "esothelin"),
    ("bio", "Pass Goparaju fibulin-3 plasma effusion pleural mesothelioma New England Journal", "ibulin"),
    ("bio", "Creaney Robinson combined CA125 mesothelin fibulin-3 mesothelioma diagnosis", "esothelin"),
    ("bio", "Pass Lott asbestos exposure pleural mesothelioma serum osteopontin", "osteopontin"),
    ("bio", "Grigoriu utility osteopontin diagnosis prognosis malignant pleural mesothelioma", "osteopontin"),
    ("bio", "Creaney Nowak Serum mesothelin monitoring treatment response mesothelioma", "esothelin"),
    ("bio", "Wheatley-Price serum mesothelin response chemotherapy mesothelioma prognostic", "esothelin"),
    ("bio", "Hollevoet renal function serum mesothelin confounding", "esothelin"),
    ("bio", "Creaney hyaluronic acid pleural fluid mesothelioma diagnostic", "yaluron"),
    ("bio", "Sneddon Creaney Malignant cells pleural fluid mesothelioma diagnosis biomarker panel", "esotheli"),

    # ---------------- 9. pleural physiology, effusion, trapped lung ----------------
    ("pleura", "Lai-Fook Pleural mechanics fluid exchange physiological reviews", "leural"),
    ("pleura", "Miserocchi Physiology pathophysiology pleural fluid turnover", "leural"),
    ("pleura", "Wang Anatomy pleura clinics chest medicine stomata lymphatic", "pleura"),
    ("pleura", "Zocchi Physiology pathophysiology pleural fluid turnover European Respiratory Journal", "leural"),
    ("pleura", "Grabczak Krenke pleural manometry pleural elastance trapped lung", "manometry"),
    ("pleura", "Feller-Kopman Ernst Large volume thoracentesis pleural pressure", "horacentesis"),
    ("pleura", "Lan Lo Chen elastance pleural space predicting pleurodesis success", "lastance"),
    ("pleura", "Huggins Doelken Sahn Pleural manometry unexpandable lung", "nexpandable"),
    ("pleura", "Doelken Sahn trapped lung entrapped lung pathophysiology", "rapped lung"),
    ("pleura", "Bhatnagar Maskell Outcomes talc pleurodesis malignant pleural effusion trial", "leurodesis"),
    ("pleura", "Dresler Olak phase III intrapleural talc insufflation slurry malignant pleural effusion", "talc"),
    ("pleura", "Davies Lee Rahman indwelling pleural catheter chest tube talc TIME2 randomized trial", "atheter"),
    ("pleura", "Thomas Lee indwelling pleural catheter talc pleurodesis AMPLE randomized clinical trial", "atheter"),
    ("pleura", "Bhatnagar Maskell Indwelling pleural catheter talc pleurodesis IPC-PLUS", "atheter"),
    ("pleura", "Roberts Maskell British Thoracic Society pleural disease guideline malignant effusion", "leural"),
    ("pleura", "Feller-Kopman ATS STS STR clinical practice guideline management malignant pleural effusions", "leural"),
    ("pleura", "Clive Maskell Predicting survival malignant pleural effusion LENT score", "ffusion"),
    ("pleura", "Zahid VEGF vascular endothelial growth factor malignant pleural effusion mesothelioma", "ffusion"),
    ("pleura", "Yano Nishioka vascular endothelial growth factor pleural effusion formation", "ffusion"),

    # ---------------- 10. angiogenesis and VEGF ----------------
    ("vegf", "Strizzi Vascular endothelial growth factor autocrine growth factor malignant mesothelioma", "esotheli"),
    ("vegf", "Ohta VEGF expression microvessel density malignant pleural mesothelioma prognosis", "esotheli"),
    ("vegf", "Kumar-Singh Van Marck vascular endothelial growth factor angiogenesis mesothelioma prognostic", "esotheli"),
    ("vegf", "Zalcman Mazieres Bevacizumab pemetrexed cisplatin malignant pleural mesothelioma MAPS randomised", "evacizumab"),
    ("vegf", "Grosso Scagliotti Nintedanib pemetrexed cisplatin LUME-Meso phase III", "intedanib"),
    ("vegf", "Scagliotti Nintedanib pemetrexed cisplatin unresectable malignant pleural mesothelioma phase II LUME-Meso", "intedanib"),
    ("vegf", "Jain normalization tumor vasculature emerging concept antiangiogenic therapy Science", "ormaliz"),
    ("vegf", "Jain Normalizing tumor microenvironment treat cancer bench bedside", "ormaliz"),

    # ---------------- 11. drug delivery, penetration, IFP ----------------
    ("deliv", "Minchinton Tannock Drug penetration solid tumours Nature Reviews Cancer", "enetration"),
    ("deliv", "Tannock Lee Limited penetration anticancer drugs through tumor tissue", "enetration"),
    ("deliv", "Thurber Wittrup mechanistic compartmental model antibody targeting binding site barrier", "ntibody"),
    ("deliv", "Thurber Schmidt Wittrup Antibody tumor penetration transport opposed systemic clearance", "enetration"),
    ("deliv", "Weinstein van Osdol binding site barrier antibody tumor penetration", "arrier"),
    ("deliv", "Jain Stylianopoulos Delivering nanomedicine solid tumors interstitial pressure", "olid tumor"),
    ("deliv", "Heldin Rubin Pietras high interstitial fluid pressure barrier drug", "nterstitial"),
    ("deliv", "Netti Berk Role extracellular matrix interstitial transport tumors", "nterstitial"),
    ("deliv", "Dedrick Myers Bungay DeVita Pharmacokinetic rationale peritoneal drug administration", "eritoneal"),
    ("deliv", "Los Ruevekamp Penetration platinum peritoneal tumor after intraperitoneal cisplatin", "ntraperitoneal"),
    ("deliv", "van de Vaart intraperitoneal cisplatin DNA adduct formation tumor penetration depth", "enetration"),
    ("deliv", "Ceelen Flessner intraperitoneal chemotherapy transport barriers pharmacokinetics", "eritoneal"),
    ("deliv", "Sugarbaker Zellos intracavitary hyperthermic cisplatin lavage mesothelioma", "esotheli"),
    ("deliv", "Ried Hyperthermic intrathoracic chemotherapy HITHOC pleural malignancies", "ntrathoracic"),
    ("deliv", "Markman intracavitary therapy pharmacologic principles", "ntracavitary"),

    # ---------------- 12. chemotherapy trials and PK ----------------
    ("chemo", "Vogelzang Rusthoven Phase III study pemetrexed cisplatin versus cisplatin alone malignant pleural mesothelioma", "emetrexed"),
    ("chemo", "van Meerbeeck Gaafar Randomized phase III raltitrexed cisplatin EORTC mesothelioma", "altitrexed"),
    ("chemo", "Muers Rudd Active symptom control chemotherapy MS01 mesothelioma randomised trial", "esotheli"),
    ("chemo", "Ceresoli Zucali pemetrexed carboplatin elderly malignant pleural mesothelioma phase II", "arboplatin"),
    ("chemo", "Santoro Vogelzang pemetrexed cisplatin carboplatin expanded access program mesothelioma", "emetrexed"),
    ("chemo", "Byrne Nowak Cisplatin gemcitabine treatment malignant mesothelioma phase II", "emcitabine"),
    ("chemo", "Nowak Byrne multicentre phase II trial cisplatin gemcitabine malignant mesothelioma", "emcitabine"),
    ("chemo", "Fennell Vinorelbine relapsed mesothelioma VIM randomised phase", "inorelbine"),
    ("chemo", "Zucali Ceresoli Vinorelbine gemcitabine pretreated malignant pleural mesothelioma", "esotheli"),
    ("chemo", "Niyikiza Hanauske pemetrexed homocysteine methylmalonic acid vitamin supplementation toxicity", "emetrexed"),
    ("chemo", "Vogelzang folic acid vitamin B12 supplementation pemetrexed toxicity mesothelioma", "emetrexed"),
    ("chemo", "Latz Karlsson population pharmacokinetic pemetrexed renal function", "emetrexed"),
    ("chemo", "Ouellet Periclou pemetrexed population pharmacokinetics renal impairment", "emetrexed"),
    ("chemo", "Urien Lokiec population pharmacokinetics total unbound plasma cisplatin", "isplatin"),
    ("chemo", "Shirazi Rogers pharmacokinetics ultrafilterable platinum cisplatin", "latin"),
    ("chemo", "Adjei pemetrexed pharmacology mechanism action multitargeted antifolate", "emetrexed"),
    ("chemo", "Shih Chen LY231514 multitargeted antifolate thymidylate synthase inhibition", "ntifolate"),
    ("chemo", "Chattopadhyay Moran Goldman pemetrexed biochemical cellular pharmacology polyglutamation", "emetrexed"),
    ("chemo", "Kelland resurgence platinum-based cancer chemotherapy Nature Reviews Cancer", "latin"),
    ("chemo", "Siddik Cisplatin mode cytotoxic action molecular basis resistance", "isplatin"),
    ("chemo", "Friberg Henningsson Model myelosuppression semi-physiological chemotherapy", "yelosuppression"),
    ("chemo", "Kloft Friberg population pharmacokinetic pharmacodynamic model neutropenia", "eutropenia"),
    ("chemo", "Calvert Newell Carboplatin dosage prospective evaluation formula renal function", "arboplatin"),

    # ---------------- 13. immunotherapy ----------------
    ("io", "Baas Scherpereel Nivolumab ipilimumab versus chemotherapy malignant pleural mesothelioma CheckMate 743", "esotheli"),
    ("io", "Peters Scherpereel first-line nivolumab ipilimumab CheckMate 743 3-year update", "esotheli"),
    ("io", "Fennell Ewings Nivolumab versus placebo relapsed malignant mesothelioma CONFIRM randomised", "ivolumab"),
    ("io", "Popat Curioni-Fontecedro multicentre randomised phase III nivolumab pembrolizumab relapsed mesothelioma", "esotheli"),
    ("io", "Chu Piccirillo Pembrolizumab chemotherapy versus chemotherapy pleural mesothelioma IND227", "embrolizumab"),
    ("io", "Forde Sun Durvalumab chemotherapy first-line mesothelioma PrE0505 phase 2", "urvalumab"),
    ("io", "Nowak Lesterhuis DREAM durvalumab first-line chemotherapy mesothelioma phase 2", "urvalumab"),
    ("io", "Scherpereel Mazieres nivolumab ipilimumab second-line MAPS2 randomised phase 2", "ivolumab"),
    ("io", "Alley Molife Pembrolizumab PD-L1 positive mesothelioma KEYNOTE-028 phase 1b", "esotheli"),
    ("io", "Maio Scherpereel Tremelimumab second-line malignant mesothelioma DETERMINE randomised", "remelimumab"),
    ("io", "Mansfield Roden PD-L1 expression sarcomatoid mesothelioma", "PD-L1"),
    ("io", "Cedres Felip PD-L1 expression malignant pleural mesothelioma prognosis", "PD-L1"),
    ("io", "Brosseau Scherpereel Shall we stop calling mesothelioma cold tumor immune", "esotheli"),
    ("io", "Chee Byrne Nowak Immune microenvironment malignant pleural mesothelioma review", "esotheli"),
    ("io", "Ujiie Adusumilli solid tumor infiltrating lymphocytes prognosis mesothelioma", "esotheli"),
    ("io", "Yamada Tumor-infiltrating regulatory T cells mesothelioma prognosis", "esotheli"),
    ("io", "Burt Zhu circulating macrophage phenotype mesothelioma prognosis", "acrophage"),
    ("io", "Marcq Pauwels prognostic immune checkpoint tumour infiltrating lymphocytes mesothelioma", "esotheli"),
    ("io", "Agarwal Zhang tumour mutational burden mesothelioma immunotherapy", "esotheli"),
    ("io", "Brahmer Drake Phase I safety clinical activity anti-PD-1 antibody receptor occupancy", "PD-1"),
    ("io", "Bajaj Gupta Model-based population pharmacokinetic analysis nivolumab", "ivolumab"),
    ("io", "Feng Roy Exposure response ipilimumab population pharmacokinetics", "pilimumab"),

    # ---------------- 14. targeted and cellular therapy ----------------
    ("target", "Szlosarek Steele ADI-PEG20 pemetrexed cisplatin ATOMIC-Meso randomized clinical trial", "rgin"),
    ("target", "Szlosarek Steele arginine deprivation pegylated arginine deiminase ADAM randomised mesothelioma", "rginine"),
    ("target", "Zauderer Krug Phase 2 tazemetostat BAP1-inactivated malignant pleural mesothelioma", "azemetostat"),
    ("target", "Fennell Defactinib maintenance therapy malignant pleural mesothelioma COMMAND randomised", "efactinib"),
    ("target", "Fennell Krebs MiST multi-arm phase 2 trial mesothelioma stratified", "esotheli"),
    ("target", "Kindler Novello Anetumab ravtansine versus vinorelbine mesothelin mesothelioma randomised", "netumab"),
    ("target", "Hassan Pastan Mesothelin targeted immunotoxin SS1P mesothelioma clinical", "esothelin"),
    ("target", "Adusumilli Zauderer regional delivery mesothelin-targeted CAR T cell therapy pleural", "CAR"),
    ("target", "Adusumilli Rusch regional intrapleural CAR T cells mesothelin phase 1", "CAR"),
    ("target", "Ceresoli Grosso Tumour treating fields chemotherapy first-line mesothelioma STELLAR", "esotheli"),
    ("target", "Zauderer Mesothelioma treatment landscape emerging therapies review", "esotheli"),
    ("target", "Yap Aerts Novel insights therapeutic targets malignant pleural mesothelioma", "esotheli"),

    # ---------------- 15. surgery and radiotherapy ----------------
    ("surg", "Treasure Lang-Lazdunski Extra-pleural pneumonectomy versus no extra-pleural pneumonectomy MARS randomised feasibility", "neumonectomy"),
    ("surg", "Lim Waller MARS 2 extended pleurectomy decortication chemotherapy mesothelioma randomised", "esotheli"),
    ("surg", "Rusch Baldini role surgical cytoreduction treatment malignant pleural mesothelioma", "urgical"),
    ("surg", "Flores Rusch Extrapleural pneumonectomy versus pleurectomy decortication survival mesothelioma", "leurectomy"),
    ("surg", "Sugarbaker Flores Resection margins extrapleural nodal status trimodality therapy mesothelioma", "esotheli"),
    ("surg", "Cao Yan systematic review extrapleural pneumonectomy malignant pleural mesothelioma", "neumonectomy"),
    ("surg", "Rintoul Treasure MARS trial feasibility surgery mesothelioma commentary", "esotheli"),
    ("surg", "Rimner Zauderer Phase II hemithoracic intensity-modulated pleural radiation therapy IMPRINT", "adiation"),
    ("surg", "Cho Feld Anraku Surgery for Mesothelioma After Radiation Therapy SMART", "SMART"),
    ("surg", "Rusch Rosenzweig Feasibility hemithoracic intensity-modulated radiation therapy mesothelioma", "adiation"),
    ("surg", "Clive Maskell Prophylactic radiotherapy prevention procedure-tract metastases SMART randomised", "adiotherapy"),
    ("surg", "Bayman Lansdown prophylactic irradiation tracts mesothelioma PIT randomised", "rradiation"),

    # ---------------- 16. guidelines, systematic reviews, QoL ----------------
    ("guide", "Popat Baas ESMO Clinical Practice Guidelines malignant pleural mesothelioma diagnosis treatment", "esotheli"),
    ("guide", "Kindler Ismaila Treatment malignant pleural mesothelioma American Society Clinical Oncology guideline", "esotheli"),
    ("guide", "Woolhouse Bishop British Thoracic Society guideline investigation management malignant pleural mesothelioma", "esotheli"),
    ("guide", "Scherpereel Opitz ERS ESTS EACTS ESTRO guidelines management malignant pleural mesothelioma", "esotheli"),
    ("guide", "Hollen Gralla Measurement quality life lung cancer symptom scale LCSS mesothelioma", "esotheli"),
    ("guide", "Nowak Chansky IASLC staging project mesothelioma database update", "esotheli"),
    ("guide", "Bibby Maskell Malignant pleural mesothelioma an update pathogenesis diagnosis treatment", "esotheli"),
    ("guide", "Janes Alrifai Perspectives on mesothelioma New England Journal Medicine review", "esotheli"),

    # ---------------- 17. QSP / modelling methodology ----------------
    ("qsp", "Baker Simulating clinical trials quantitative systems pharmacology oncology", "harmacology"),
    ("qsp", "Chelliah Kimko Quantitative systems pharmacology approaches oncology drug development", "harmacology"),
    ("qsp", "Simeoni Rocchetti predictive pharmacokinetic pharmacodynamic model tumor growth inhibition", "umor growth"),
    ("qsp", "Claret Bruno model-based prediction phase III overall survival tumor size", "urvival"),
    ("qsp", "Bruno Mercier Bhatt model-based drug development oncology tumour growth inhibition metrics", "umour"),
    ("qsp", "Ribba Holford Mentre grading system tumour growth inhibition models", "umour"),
    ("qsp", "Wang Sung Elassaiss-Schaap exposure response tumor size overall survival modeling framework", "urvival"),
    ("qsp", "Kosinsky Bai quantitative systems pharmacology model immuno-oncology checkpoint", "mmuno"),
    ("qsp", "Elmokadem mrgsolve simulation pharmacometrics R", "mrgsolve"),
    ("qsp", "Gaudreau Novakovic mechanistic model tumor immune checkpoint inhibitor response", "umor"),
]


STOP = set("""a an and the of for in on to with without versus vs by from at as is are
randomized randomised phase trial study review clinical journal new england lancet nature
cancer oncology medicine research""".split())


def content_words(term):
    return [w.lower() for w in re.findall(r"[A-Za-z0-9-]{5,}", term)
            if w.lower() not in STOP]


def overlap_ok(title, term, need=2):
    t = title.lower()
    cw = content_words(term)
    return sum(1 for w in cw if w in t) >= min(need, len(cw))


# Words too generic to identify a paper on their own.  A relaxed-tier hit is
# only accepted if the returned title shares a token with the query that is NOT
# on this list -- without it, "mesothelioma incidence United States trends"
# happily matches a paper on the incidence of osteoporotic fracture in the
# United States, which is exactly the failure this guard exists to stop.
GENERIC = set("""incidence exposure exposures patients patient treatment treatments
analysis analyses outcome outcomes disease diseases therapy therapies effect effects
response responses clinical national registry population populations survival
prognostic prognosis mortality malignant management multicentre multicenter
international guidelines guideline update updated systematic overall national
american british european japanese chinese italian turkish""".split())



# A relaxed-tier hit must ALSO land inside the subject matter.  Without this,
# "mesothelioma incidence United States trends" matched a paper on the incidence
# of osteoporotic fracture in the United States: "united" and "states" are
# non-generic six-letter tokens, so the rare-token guard alone was not enough.
ANCHORS = """mesotheli mesothelial mesothelin pleura pleural pleuro intrapleural
pleurodesis pleurectomy pneumonectomy decortication thoracoscop thoracentesis
thoracic intrathoracic hemithora asbest erionite amphibole chrysotile fibre fiber
effusion empyema lung respiratory pulmonary
bap1 nf2 cdkn2a mtap prmt5 ezh2 yap taz hippo merlin calretinin telomerase tert
inflammasome nlrp3 hmgb1 macrophage lymphatic stomata
pemetrexed cisplatin carboplatin platinum antifolate gemcitabine vinorelbine
raltitrexed nivolumab ipilimumab pembrolizumab durvalumab tremelimumab
bevacizumab nintedanib arginine deiminase tazemetostat defactinib anetumab
ravtansine immunotox chimeric antigen checkpoint pd-1 pd-l1 ctla-4 interferon
tumor tumour tumour-treating oncolog carcinogen neoplas malign cancer
penetration interstitial intraperitoneal peritoneal intracavitary hyperthermic
pharmacokinet pharmacodynam pharmacolog myelosuppression neutropenia folate
folic elastance manometry trapped unexpandable catheter talc radiotherapy
radiation recist volumetr imaging computed tomography biomarker prognostic
angiogenesis vascular endothelial growth vegf transforming fibrosis collagen
model modelling modeling simulation survival""".split()


def on_topic(title):
    t = title.lower()
    return any(w in t for w in ANCHORS)


def rare_tokens(term):
    return [w for w in content_words(term) if w not in GENERIC and len(w) >= 6]


def try_query(term, need_sub, tier):
    """Tier 0: full query, title must contain need_sub.
       Tier 1: full query, no substring filter, >=2 content words must overlap.
       Tier 2: author surnames dropped.
       Tier 3: six most distinctive words only.
       Tiers 1-3 additionally require the returned title to contain at least one
       NON-GENERIC token from the query, so a relaxed search cannot drift onto a
       different disease that happens to share filler words."""
    toks = term.split()
    if tier == 0:
        q, sub, ov = term, need_sub, 0
    elif tier == 1:
        q, sub, ov = term, "", 2
    elif tier == 2:
        q, sub, ov = " ".join(toks[2:]) if len(toks) > 5 else term, "", 2
    else:
        cw = content_words(term)
        q, sub, ov = " ".join(cw[:6]), "", 3
    if not q.strip():
        return None
    ids = esearch(q, n=3)
    if not ids:
        return None
    summ = esummary(ids)
    rare = rare_tokens(term)
    for pid in ids:
        r = summ.get(pid)
        if not r or "title" not in r:
            continue
        title = re.sub(r"<[^>]+>", "", r.get("title", "")).strip()
        if sub and sub.lower() not in title.lower():
            continue
        if ov and not overlap_ok(title, term, ov):
            continue
        if tier > 0 and rare and not any(w in title.lower() for w in rare):
            continue
        if tier > 0 and not on_topic(title):
            continue
        return (pid, r, title, tier)
    return None


def main():
    seen = {}
    misses = []
    sections = {}
    for i, (sec, term, need) in enumerate(QUERIES, 1):
        hit = None
        for tier in (0, 1, 2, 3):
            hit = try_query(term, need, tier)
            time.sleep(0.34)
            if hit:
                break
        if not hit:
            misses.append((sec, term, "no PubMed record matched at any tier"))
            print("%3d/%d  MISS  %s" % (i, len(QUERIES), term[:64]), file=sys.stderr)
            continue
        pid, r, title, tier = hit
        if pid in seen:
            print("%3d/%d  DUP   %s" % (i, len(QUERIES), term[:64]), file=sys.stderr)
            continue
        seen[pid] = True
        au = r.get("sortfirstauthor") or (r.get("authors") or [{}])[0].get("name", "")
        yr = (r.get("pubdate") or "")[:4]
        jr = r.get("source", "")
        sections.setdefault(sec, []).append((au, yr, title, jr, pid))
        print("%3d/%d  ok(t%d) %s (%s)" % (i, len(QUERIES), tier, title[:54], pid),
              file=sys.stderr)

    write_markdown(sections, misses)


SECTION_TITLES = [
    ("epi", "1. Epidemiology, exposure and disease burden"),
    ("fibre", "2. Fibre biology, biopersistence and the inflammatory field"),
    ("gene", "3. Genomic landscape and tumour-suppressor loss"),
    ("hippo", "4. Hippo-YAP signalling, mechanotransduction and the cell cycle"),
    ("path", "5. Histopathology, subtype and the epithelial-mesenchymal axis"),
    ("imag", "6. Imaging, mRECIST and tumour volumetry"),
    ("stage", "7. Staging and prognostic models"),
    ("bio", "8. Biomarkers (mesothelin, fibulin-3, osteopontin)"),
    ("pleura", "9. Pleural physiology, effusion, manometry and pleurodesis"),
    ("vegf", "10. VEGF, angiogenesis and anti-angiogenic therapy"),
    ("deliv", "11. Drug delivery, tissue penetration and intracavitary pharmacology"),
    ("chemo", "12. Cytotoxic chemotherapy: trials, pharmacology and PK/PD"),
    ("io", "13. Immune contexture and checkpoint blockade"),
    ("target", "14. Targeted, cellular and device-based therapy"),
    ("surg", "15. Surgery and radiotherapy"),
    ("guide", "16. Guidelines, reviews and patient-reported outcomes"),
    ("qsp", "17. QSP and pharmacometric methodology"),
]


def write_markdown(sections, misses):
    n = sum(len(v) for v in sections.values())
    lines = []
    lines.append("# Malignant Pleural Mesothelioma (MPM) — QSP Model References")
    lines.append("")
    lines.append("악성 흉막 중피종 QSP 모델의 근거 문헌 목록입니다. "
                 "**총 %d편**이며, 모든 항목은 `mpm_fetch_references.py`가 "
                 "NCBI E-utilities로 실제 조회하여 반환받은 레코드만 기록한 것입니다 "
                 "(제목·저널·연도·PMID 모두 PubMed 응답 그대로)." % n)
    lines.append("")
    lines.append("Every entry below was returned by a live PubMed query; nothing is "
                 "written from memory. Regenerate with `python3 mpm_fetch_references.py`.")
    lines.append("")
    lines.append("---")
    lines.append("")
    for key, title in SECTION_TITLES:
        items = sections.get(key, [])
        if not items:
            continue
        lines.append("## %s" % title)
        lines.append("")
        lines.append("<sub>%d편</sub>" % len(items))
        lines.append("")
        for au, yr, ti, jr, pid in items:
            lines.append("- **%s (%s)** — %s *%s*. "
                         "[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)"
                         % (au or "—", yr or "—", ti, jr, pid, pid))
        lines.append("")
    if misses:
        lines.append("---")
        lines.append("")
        lines.append("## Queries that returned nothing usable")
        lines.append("")
        lines.append("이 검색어들은 PubMed가 조건에 맞는 레코드를 돌려주지 않았습니다. "
                     "기억으로 채우지 않고 그대로 남겨 둡니다.")
        lines.append("")
        for sec, term, why in misses:
            lines.append("- `[%s]` %s — %s" % (sec, term, why))
        lines.append("")
    open("mpm_references.md", "w").write("\n".join(lines) + "\n")
    print("\nwrote mpm_references.md with %d references (%d misses)" % (n, len(misses)),
          file=sys.stderr)


if __name__ == "__main__":
    main()
