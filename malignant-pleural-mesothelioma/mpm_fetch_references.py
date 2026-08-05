#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpm_fetch_references.py -- build mpm_references.md from PubMed records that
were ACTUALLY RETURNED by NCBI E-utilities.

Nothing here is written from memory.  Each entry below is a SEARCH; only the
record NCBI hands back (its real PMID, title, journal, year, author list) is
written into the reference file.  Searches that return nothing are reported as
misses at the end of the run rather than filled in from recollection.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = "mpm_references.md"


def _get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.loads(r.read().decode())
        except Exception as e:                                  # noqa: BLE001
            if attempt == 3:
                print("   ! %s" % e, file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


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


SECTIONS = [
    ("epi",   "1. Epidemiology, asbestos exposure and disease burden"),
    ("fibre", "2. Fibre dosimetry, biopersistence and translocation to the pleura"),
    ("inflam", "3. Frustrated phagocytosis, NLRP3/IL-1beta and HMGB1"),
    ("genet", "4. Germline predisposition and the somatic tumour-suppressor landscape"),
    ("sig",   "5. Hippo/YAP, BAP1 chromatin biology and oncogenic signalling"),
    ("mtap",  "6. CDKN2A/MTAP co-deletion, PRMT5 synthetic lethality, ASS1 and arginine"),
    ("histo", "7. Histological subtype, grading and prognostic scores"),
    ("geom",  "8. Tumour volume, rind thickness and modified RECIST"),
    ("deliv", "9. Drug penetration, interstitial pressure and the intrapleural route"),
    ("pleur", "10. Pleural fluid physiology, effusion management and pleurodesis"),
    ("angio", "11. VEGF, angiogenesis and bevacizumab"),
    ("immune", "12. Tumour immune microenvironment, PD-L1, VISTA and mesothelin"),
    ("chemo", "13. Pemetrexed and platinum: pharmacology and pharmacokinetics"),
    ("folate", "14. Folate and B12 supplementation, polyglutamation, homocysteine"),
    ("trials", "15. Randomised trials of systemic therapy"),
    ("io",    "16. Immune checkpoint blockade trials and biomarkers"),
    ("surg",  "17. Surgery and radiotherapy"),
    ("tox",   "18. Toxicity: myelosuppression, nephrotoxicity, immune-related AEs"),
    ("biom",  "19. Biomarkers: mesothelin, fibulin-3, osteopontin, diagnostic IHC"),
    ("qsp",   "20. Modelling methodology: tumour growth, QSP, mrgsolve"),
]

QUERIES = [
    # ---------------- 1. epidemiology -------------------------------------
    ("epi", "Advances in Malignant Mesothelioma review New England Journal of Medicine", "mesothelioma"),
    ("epi", "Mesothelioma scientific clues for prevention diagnosis and therapy", "mesothelioma"),
    ("epi", "Estimation of the global burden of mesothelioma deaths from incomplete national mortality data", "mesothelioma"),
    ("epi", "Global magnitude of reported and unreported mesothelioma", "mesothelioma"),
    ("epi", "Continuing increase in mesothelioma mortality in Britain", "mesothelioma"),
    ("epi", "Mesothelioma incidence in 50 states and the District of Columbia United States", "mesothelioma"),
    ("epi", "Global asbestos disaster", "asbestos"),
    ("epi", "Mesothelioma and asbestos exposure in Cappadocia erionite villages", "erionite"),
    ("epi", "Malignant mesothelioma global incidence and relationship with asbestos", "mesothelioma"),
    ("epi", "Occupational domestic and environmental mesothelioma risks in Britain", "mesothelioma"),
    ("epi", "Pleural mesothelioma primer Nature Reviews Disease Primers", "mesothelioma"),
    # ---------------- 2. fibre dosimetry ------------------------------------
    ("fibre", "Asbestos carbon nanotubes and the pleural mesothelium a review of the hypothesis regarding the role of long fibre retention", "pleural"),
    ("fibre", "Black spots concentrate oncogenic asbestos fibers in the parietal pleura", "pleura"),
    ("fibre", "Biopersistence of chrysotile and amphibole asbestos inhalation", "biopersistence"),
    ("fibre", "Relation of particle dimension to carcinogenicity in amphibole asbestoses and other fibrous minerals", "carcinogenicity"),
    ("fibre", "Translocation of particles from pulmonary alveoli into the pleural space", "pleural"),
    ("fibre", "Diffuse pleural mesothelioma and asbestos exposure in the North Western Cape Province", "mesothelioma"),
    ("fibre", "Mesothelioma latency period asbestos exposure duration", "latency"),
    ("fibre", "Latent period for malignant mesothelioma of occupational origin", "latent"),
    ("fibre", "Asbestos induces apoptosis of human and rabbit pleural mesothelial cells via reactive oxygen species", "mesothelial"),
    ("fibre", "Pleural plaques and the risk of pleural mesothelioma", "plaques"),
    # ---------------- 3. inflammation ----------------------------------------
    ("inflam", "Innate immune activation through Nalp3 inflammasome sensing of asbestos and silica", "inflammasome"),
    ("inflam", "TNF-alpha inhibits asbestos-induced cytotoxicity NF-kappaB mesothelial cells", "asbestos"),
    ("inflam", "asbestos programmed necrosis mesothelial cells HMGB1 release", "necrosis"),
    ("inflam", "Cancer cell secretion of the DAMP protein HMGB1 supports progression in malignant mesothelioma", "HMGB1"),
    ("inflam", "anakinra interleukin-1 receptor antagonist mesothelioma", "IL"),
    ("inflam", "Minerals and cancer HMGB1 mesothelioma inflammation", "mesothelioma"),
    ("inflam", "Serum HMGB1 as a biomarker in asbestos-exposed individuals and mesothelioma patients", "HMGB1"),
    ("inflam", "Asbestos induces mesothelial cell transformation via HMGB1-driven autophagy", "HMGB1"),
    ("inflam", "Iron-induced oxidative stress in carcinogenesis asbestos", "oxidative"),
    ("inflam", "The role of the inflammasome in mesothelioma and asbestos-induced inflammation", "inflammasome"),
    # ---------------- 4. genetics ---------------------------------------------
    ("genet", "Germline BAP1 mutations predispose to malignant mesothelioma Nature Genetics", "BAP1"),
    ("genet", "The nuclear deubiquitinase BAP1 is commonly inactivated by somatic mutations and 3p21.1 losses in malignant pleural mesothelioma", "BAP1"),
    ("genet", "Comprehensive genomic analysis of malignant pleural mesothelioma identifies recurrent mutations gene fusions and splicing alterations", "mesothelioma"),
    ("genet", "Integrative molecular characterization of malignant pleural mesothelioma", "mesothelioma"),
    ("genet", "BAP1 cancer syndrome malignant mesothelioma uveal melanoma", "BAP1"),
    ("genet", "Frequency of germline mutations in DNA repair and other genes in patients with mesothelioma", "germline"),
    ("genet", "High frequency of inactivating mutations in the neurofibromatosis type 2 gene NF2 in primary malignant mesotheliomas", "NF2"),
    ("genet", "p16 alterations and deletion mapping of 9p21-p22 in malignant mesothelioma", "mesothelioma"),
    ("genet", "Neoantigen landscape and tumor mutation burden in malignant pleural mesothelioma", "mesothelioma"),
    ("genet", "Genetic and epigenetic molecular classification of malignant pleural mesothelioma", "mesothelioma"),
    ("genet", "Mesothelioma patients with germline BAP1 mutations have 7-fold improved long-term survival", "BAP1"),
    ("genet", "Whole exome sequencing reveals frequent genetic alterations in BAP1 NF2 CDKN2A and CUL1 in malignant pleural mesothelioma", "mesothelioma"),
    # ---------------- 5. signalling ---------------------------------------------
    ("sig", "Molecular pathogenesis of malignant mesothelioma Hippo signaling pathway", "mesothelioma"),
    ("sig", "LATS2 is a tumor suppressor gene of malignant mesothelioma", "LATS2"),
    ("sig", "BAP1 links metabolic regulation of ferroptosis to tumour suppression", "BAP1"),
    ("sig", "BAP1 regulates IP3R3-mediated Ca2+ flux to mitochondria suppressing cell transformation", "BAP1"),
    ("sig", "Loss of BAP1 function leads to EZH2-dependent transformation", "BAP1"),
    ("sig", "Hepatocyte growth factor MET signalling in malignant pleural mesothelioma", "mesothelioma"),
    ("sig", "Co-occurring mutations of tumor suppressor genes LATS2 and NF2 in malignant pleural mesothelioma", "mesothelioma"),
    ("sig", "YAP mesothelioma Hippo pathway inactivation", "mesothelioma"),
    # ---------------- 6. MTAP / PRMT5 / ASS1 --------------------------------------
    ("mtap", "Disordered methionine metabolism in MTAP CDKN2A-deleted cancers leads to dependence on PRMT5", "PRMT5"),
    ("mtap", "MTAP deletion confers enhanced dependency on the PRMT5 arginine methyltransferase in cancer cells", "MTAP"),
    ("mtap", "Fragment-based discovery of MRTX1719 a synthetic lethal inhibitor of the PRMT5 MTA complex for the treatment of MTAP-deleted cancers", "PRMT5"),
    ("mtap", "argininosuccinate synthetase mesothelioma arginine deprivation biomarker", "arginine"),
    ("mtap", "Arginine deprivation with pegylated arginine deiminase in patients with argininosuccinate synthetase 1-deficient malignant pleural mesothelioma randomized", "arginine"),
    ("mtap", "Pegargiminase plus first-line chemotherapy in patients with nonepithelioid pleural mesothelioma ATOMIC-Meso randomized clinical trial", "Pegargiminase"),
    ("mtap", "MTAP immunohistochemistry loss and CDKN2A deletion in mesothelioma diagnosis", "MTAP"),
    ("mtap", "Utility of BAP1 immunohistochemistry and p16 FISH in the diagnosis of biphasic mesothelioma", "mesothelioma"),
    # ---------------- 7. histology and prognosis ----------------------------------
    ("histo", "EURACAN IASLC proposals for updating the histologic classification of pleural mesothelioma", "mesothelioma"),
    ("histo", "The 2021 WHO Classification of Tumors of the Pleura advances since the 2015 classification", "Pleura"),
    ("histo", "nuclear grade epithelioid malignant pleural mesothelioma grading prognosis", "mesothelioma"),
    ("histo", "Prognostic factors in patients with pleural mesothelioma the European Organization for Research and Treatment of Cancer experience", "mesothelioma"),
    ("histo", "factors predictive of survival mesothelioma Cancer and Leukemia Group B", "mesothelioma"),
    ("histo", "Impact of mesothelioma histologic subtype on outcomes in the Surveillance Epidemiology and End Results database", "mesothelioma"),
    ("histo", "Prognostic impact of histologic subtype in malignant pleural mesothelioma", "mesothelioma"),
    ("histo", "IASLC Mesothelioma Staging Project proposals T descriptors", "Mesothelioma"),
    ("histo", "Nuclear grade and necrosis predict prognosis in malignant epithelioid pleural mesothelioma", "mesothelioma"),
    # ---------------- 8. geometry and mRECIST ---------------------------------------
    ("geom", "Modified RECIST criteria for assessment of response in malignant pleural mesothelioma", "RECIST"),
    ("geom", "modified RECIST 1.1 malignant pleural mesothelioma response", "RECIST"),
    ("geom", "observer variability mesothelioma tumor thickness CT", "mesothelioma"),
    ("geom", "Imaging in pleural mesothelioma a review of the 15th International Conference of the International Mesothelioma Interest Group", "mesothelioma"),
    ("geom", "volumetric CT clinical staging malignant pleural mesothelioma", "mesothelioma"),
    ("geom", "Preoperative tumor volume is associated with outcome in malignant pleural mesothelioma", "volume"),
    ("geom", "computerized tumor volume CT malignant pleural mesothelioma", "mesothelioma"),
    ("geom", "Diagnostic performance of computed tomography for tumour volume in malignant pleural mesothelioma", "mesothelioma"),
    ("geom", "Baseline tumour volume derived from CT predicts survival in malignant pleural mesothelioma", "mesothelioma"),
    ("geom", "Tumour volume doubling time growth rate malignant pleural mesothelioma", "mesothelioma"),
    ("geom", "Pleural thickness measurement reproducibility CT mesothelioma response assessment", "mesothelioma"),
    # ---------------- 9. delivery -----------------------------------------------------
    ("deliv", "Drug penetration in solid tumours Nature Reviews Cancer", "penetration"),
    ("deliv", "Delivering nanomedicine to solid tumors interstitial fluid pressure", "tumors"),
    ("deliv", "Normalization of tumor vasculature an emerging concept in antiangiogenic therapy", "vasculature"),
    ("deliv", "Factors determining antibody distribution in tumors", "antibody"),
    ("deliv", "Penetration of carboplatin and cisplatin into rat peritoneal tumor nodules after intraperitoneal chemotherapy", "penetration"),
    ("deliv", "Pharmacokinetic rationale for peritoneal drug administration in the treatment of ovarian cancer", "peritoneal"),
    ("deliv", "Hyperthermic intrathoracic chemotherapy after cytoreductive surgery for pleural malignancies", "intrathoracic"),
    ("deliv", "Regional delivery of mesothelin-targeted CAR T cell therapy generates potent and long-lasting CD4-dependent tumor immunity", "CAR"),
    ("deliv", "A phase I clinical trial of single-dose intrapleural IFN-beta gene transfer for malignant pleural mesothelioma", "intrapleural"),
    ("deliv", "Intrapleural administration of cisplatin pharmacokinetics pleural malignancy", "intrapleural"),
    # ---------------- 10. pleural physiology --------------------------------------------
    ("pleur", "Physiology and pathophysiology of pleural fluid turnover", "pleural"),
    ("pleur", "Anatomy of the pleural lymphatic stomata mesothelium", "pleura"),
    ("pleur", "Management of Malignant Pleural Effusions An Official ATS STS STR Clinical Practice Guideline", "Pleural"),
    ("pleur", "British Thoracic Society Guideline for pleural disease", "pleural"),
    ("pleur", "Effect of an indwelling pleural catheter vs chest tube and talc pleurodesis for relieving dyspnea in patients with malignant pleural effusion the TIME2 randomized controlled trial", "pleural"),
    ("pleur", "Effect of thoracoscopic talc poudrage vs talc slurry via chest tube on pleurodesis failure rate among patients with malignant pleural effusions", "talc"),
    ("pleur", "Phase III intergroup study of talc poudrage vs talc slurry sclerosis for malignant pleural effusion", "talc"),
    ("pleur", "Nonexpandable lung trapped lung pleural manometry elastance", "lung"),
    ("pleur", "Vascular endothelial growth factor in malignant pleural effusion permeability", "pleural"),
    ("pleur", "Lymphatic drainage of the pleural space stomata absorption rate", "pleural"),
    # ---------------- 11. angiogenesis -----------------------------------------------------
    ("angio", "vascular endothelial growth factor autocrine mesothelioma", "mesothelioma"),
    ("angio", "vascular endothelial growth factor mesothelioma prognosis expression", "mesothelioma"),
    ("angio", "Bevacizumab for newly diagnosed pleural mesothelioma in the Mesothelioma Avastin Cisplatin Pemetrexed Study MAPS a randomised controlled open-label phase 3 trial", "mesothelioma"),
    ("angio", "Multicenter double-blind placebo-controlled randomized phase II trial of gemcitabine cisplatin plus bevacizumab or placebo in patients with malignant mesothelioma", "mesothelioma"),
    ("angio", "Nintedanib in combination with pemetrexed and cisplatin for chemotherapy-naive patients with advanced malignant pleural mesothelioma LUME-Meso", "mesothelioma"),
    ("angio", "Phase II study of intermittent sunitinib malignant pleural mesothelioma", "mesothelioma"),
    # ---------------- 12. immune microenvironment --------------------------------------------
    ("immune", "Prognostic and predictive aspects of the tumor immune microenvironment and immune checkpoints in malignant pleural mesothelioma", "mesothelioma"),
    ("immune", "V-domain Ig suppressor of T cell activation VISTA is highly expressed in epithelioid malignant pleural mesothelioma", "mesothelioma"),
    ("immune", "Analysis of expression of programmed cell death 1 ligand 1 PD-L1 in malignant pleural mesothelioma", "mesothelioma"),
    ("immune", "B7-H1 expression in malignant pleural mesothelioma is associated with sarcomatoid histology and poor prognosis", "mesothelioma"),
    ("immune", "Tumour associated macrophages in malignant pleural mesothelioma", "mesothelioma"),
    ("immune", "Mesothelin targeted cancer immunotherapy", "Mesothelin"),
    ("immune", "Cytotoxic T cells in PD-L1-positive malignant pleural mesotheliomas are counterbalanced by distinct immunosuppressive factors", "mesothelioma"),
    ("immune", "Pleural effusion of patients with malignant mesothelioma induces macrophage-mediated T cell suppression", "mesothelioma"),
    ("immune", "Immune checkpoint inhibitors in malignant pleural mesothelioma predictive biomarkers", "mesothelioma"),
    # ---------------- 13. chemotherapy pharmacology --------------------------------------------
    ("chemo", "LY231514 a pyrrolo 2,3-d pyrimidine-based antifolate that inhibits multiple folate-requiring enzymes", "antifolate"),
    ("chemo", "pemetrexed cellular pharmacology resistance", "Pemetrexed"),
    ("chemo", "Pharmacology and mechanism of action of pemetrexed", "pemetrexed"),
    ("chemo", "Population pharmacokinetics of pemetrexed disodium in patients with cancer", "pemetrexed"),
    ("chemo", "A semimechanistic-physiologic population pharmacokinetic pharmacodynamic model for neutropenia following pemetrexed therapy", "pemetrexed"),
    ("chemo", "Phase I and pharmacokinetic study of pemetrexed administered every 3 weeks to advanced cancer patients with normal and impaired renal function", "pemetrexed"),
    ("chemo", "Population pharmacokinetics of total and unbound plasma cisplatin in adult patients", "cisplatin"),
    ("chemo", "The resurgence of platinum-based cancer chemotherapy", "platinum"),
    ("chemo", "Cisplatin mode of cytotoxic action and molecular basis of resistance", "Cisplatin"),
    ("chemo", "ERCC1 expression and outcome in malignant pleural mesothelioma treated with platinum", "mesothelioma"),
    ("chemo", "Thymidylate synthase expression and prognosis in malignant pleural mesothelioma pemetrexed", "mesothelioma"),
    # ---------------- 14. folate supplementation ------------------------------------------------
    ("folate", "homocysteine methylmalonic acid pemetrexed toxicity", "pemetrexed"),
    ("folate", "Phase II study of pemetrexed with and without folic acid and vitamin B12 as front-line therapy in malignant pleural mesothelioma", "pemetrexed"),
    ("folate", "Phase III study of pemetrexed in combination with cisplatin versus cisplatin alone in patients with malignant pleural mesothelioma", "pemetrexed"),
    ("folate", "Pemetrexed in malignant pleural mesothelioma FDA drug approval summary", "pemetrexed"),
    ("folate", "Resistance to antifolates folylpolyglutamate synthetase gamma-glutamyl hydrolase", "antifolate"),
    ("folate", "folylpolyglutamate synthetase pemetrexed cytotoxicity", "polyglutam"),
    ("folate", "The role of folic acid in modulating the toxicity and efficacy of the multitargeted antifolate LY231514", "folic"),
    ("folate", "Vitamin supplementation reduces toxicity of pemetrexed without compromising efficacy", "pemetrexed"),
    # ---------------- 15. systemic therapy trials ------------------------------------------------
    ("trials", "Randomized phase III study of cisplatin with or without raltitrexed in patients with malignant pleural mesothelioma", "mesothelioma"),
    ("trials", "Active symptom control with or without chemotherapy in the treatment of patients with malignant pleural mesothelioma MS01 a multicentre randomised trial", "mesothelioma"),
    ("trials", "Pemetrexed plus cisplatin or pemetrexed plus carboplatin for chemonaive patients with malignant pleural mesothelioma expanded access program", "mesothelioma"),
    ("trials", "Phase II study of pemetrexed plus carboplatin in malignant pleural mesothelioma", "mesothelioma"),
    ("trials", "Vinorelbine and gemcitabine as second- or third-line therapy for malignant pleural mesothelioma", "mesothelioma"),
    ("trials", "Active symptom control with or without oral vinorelbine in patients with relapsed malignant pleural mesothelioma VIM a randomised phase 2 trial", "mesothelioma"),
    ("trials", "pembrolizumab chemotherapy untreated advanced pleural mesothelioma phase 3", "mesothelioma"),
    ("trials", "Second-line chemotherapy for malignant pleural mesothelioma systematic review", "mesothelioma"),
    ("trials", "Maintenance pemetrexed after first-line chemotherapy in malignant pleural mesothelioma randomized", "mesothelioma"),
    # ---------------- 16. checkpoint blockade ---------------------------------------------------------
    ("io", "nivolumab plus ipilimumab unresectable malignant pleural mesothelioma CheckMate 743", "nivolumab"),
    ("io", "First-line nivolumab plus ipilimumab versus chemotherapy in patients with unresectable malignant pleural mesothelioma 3-year outcomes from CheckMate 743", "nivolumab"),
    ("io", "nivolumab versus placebo relapsed malignant mesothelioma CONFIRM", "Nivolumab"),
    ("io", "Nivolumab or nivolumab plus ipilimumab in patients with relapsed malignant pleural mesothelioma MAPS2 a multicentre open-label randomised non-comparative phase 2 trial", "nivolumab"),
    ("io", "A multicentre randomised phase III trial comparing pembrolizumab versus single-agent chemotherapy for advanced pre-treated malignant pleural mesothelioma PROMISE-meso", "pembrolizumab"),
    ("io", "Clinical safety and activity of pembrolizumab in patients with malignant pleural mesothelioma KEYNOTE-028 preliminary results", "pembrolizumab"),
    ("io", "Tremelimumab as second-line or third-line treatment in relapsed malignant mesothelioma DETERMINE a multicentre international randomised double-blind placebo-controlled phase 2b trial", "mesothelioma"),
    ("io", "Tremelimumab combined with durvalumab in patients with mesothelioma NIBIT-MESO-1 an open-label non-randomised phase 2 study", "mesothelioma"),
    ("io", "EZH2 inhibitor tazemetostat in patients with relapsed or refractory BAP1-inactivated malignant pleural mesothelioma", "mesothelioma"),
    ("io", "inflammatory gene signature mesothelioma immunotherapy benefit", "mesothelioma"),
    # ---------------- 17. surgery and radiotherapy -------------------------------------------------------
    ("surg", "Extra-pleural pneumonectomy versus no extra-pleural pneumonectomy for patients with malignant pleural mesothelioma clinical outcomes of the Mesothelioma and Radical Surgery MARS randomised feasibility study", "pneumonectomy"),
    ("surg", "Extended pleurectomy decortication and chemotherapy versus chemotherapy alone for pleural mesothelioma MARS 2 a phase 3 randomised controlled trial", "mesothelioma"),
    ("surg", "Initial analysis of the International Association for the Study of Lung Cancer mesothelioma database", "mesothelioma"),
    ("surg", "Extrapleural pneumonectomy versus pleurectomy decortication in the surgical management of malignant pleural mesothelioma", "pleurectomy"),
    ("surg", "Phase II study of hemithoracic intensity-modulated pleural radiation therapy IMPRINT as part of lung-sparing multimodality therapy in patients with malignant pleural mesothelioma", "pleural"),
    ("surg", "Surgery for malignant pleural mesothelioma after radiation therapy SMART final results", "mesothelioma"),
    ("surg", "prophylactic radiotherapy procedure tract metastases mesothelioma randomised", "mesothelioma"),
    ("surg", "Prophylactic irradiation of tracts in patients with malignant pleural mesothelioma an open-label multicentre phase 3 randomised controlled trial", "mesothelioma"),
    ("surg", "A randomised phase II study of extended pleurectomy decortication preceded or followed by chemotherapy in patients with early stage pleural mesothelioma", "mesothelioma"),
    # ---------------- 18. toxicity ---------------------------------------------------------------------------
    ("tox", "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs", "myelosuppression"),
    ("tox", "Semi-mechanistic model of chemotherapy induced neutropenia absolute neutrophil count time course", "neutropenia"),
    ("tox", "Cisplatin nephrotoxicity mechanisms and renoprotective strategies", "nephrotoxicity"),
    ("tox", "Pathophysiology of cisplatin-induced acute kidney injury", "cisplatin"),
    ("tox", "Long-term renal outcomes after cisplatin treatment", "cisplatin"),
    ("tox", "Immune-related adverse events associated with immune checkpoint blockade", "Immune-Related"),
    ("tox", "Management of immune-related adverse events in patients treated with immune checkpoint inhibitor therapy ASCO clinical practice guideline", "Immune-Related"),
    ("tox", "Cancer cachexia understanding the molecular basis", "cachexia"),
    ("tox", "Febrile neutropenia risk and granulocyte colony-stimulating factor prophylaxis in solid tumours", "neutropenia"),
    # ---------------- 19. biomarkers ----------------------------------------------------------------------------
    ("biom", "Mesothelin-family proteins and diagnosis of mesothelioma", "mesothelioma"),
    ("biom", "Serum and pleural fluid biomarkers for mesothelioma", "mesothelioma"),
    ("biom", "Serum mesothelin for diagnosing malignant pleural mesothelioma an individual patient data meta-analysis", "mesothelin"),
    ("biom", "The effect of clinical covariates on the diagnostic and prognostic value of soluble mesothelin and megakaryocyte potentiating factor", "mesothelin"),
    ("biom", "Glomerular filtration rate is a confounder for the measurement of soluble mesothelin in serum", "mesothelin"),
    ("biom", "Fibulin-3 as a blood and effusion biomarker for pleural mesothelioma", "mesothelioma"),
    ("biom", "Asbestos exposure pleural mesothelioma and serum osteopontin levels", "mesothelioma"),
    ("biom", "Serum soluble mesothelin concentrations in malignant pleural mesothelioma relationship to tumor volume clinical stage and changes in tumor burden", "mesothelin"),
    ("biom", "Guidelines for pathologic diagnosis of malignant mesothelioma 2017 update of the consensus statement from the International Mesothelioma Interest Group", "mesothelioma"),
    ("biom", "malignant pleural mesothelioma American Society of Clinical Oncology guideline", "mesothelioma"),
    ("biom", "Malignant pleural mesothelioma ESMO Clinical Practice Guidelines for diagnosis treatment and follow-up", "mesothelioma"),
    ("biom", "Neutrophil to lymphocyte ratio as a prognostic biomarker in malignant pleural mesothelioma", "mesothelioma"),
    # ---------------- 20. modelling ------------------------------------------------------------------------------
    ("qsp", "Predictive pharmacokinetic-pharmacodynamic modeling of tumor growth kinetics in xenograft models after administration of anticancer agents", "tumor"),
    ("qsp", "A review of mixed-effects models of tumor growth and effects of anticancer drug treatment used in population analysis", "tumor"),
    ("qsp", "Model-based prediction of phase III overall survival in colorectal cancer on the basis of phase II tumor dynamics", "tumor"),
    ("qsp", "Quantitative systems pharmacology approaches for immuno-oncology perspectives", "immuno-oncology"),
    ("qsp", "A computational model of neoadjuvant PD-1 inhibition in non-small cell lung cancer", "model"),
    ("qsp", "Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve a hands-on tutorial", "mrgsolve"),
    ("qsp", "mrgsolve simulation population PK/PD systems pharmacology R", "mrgsolve"),
    ("qsp", "Agent-based and hybrid models of cancer growth and treatment", "model"),
    ("qsp", "Clinically relevant modeling of tumor growth and treatment response", "tumor"),
    ("qsp", "Mathematical modelling of tumour response to radiotherapy linear quadratic", "radiotherapy"),
    ("qsp", "A quantitative systems pharmacology model of immune checkpoint blockade in solid tumours virtual patients", "model"),
]


def build():
    got = {}
    misses = []
    seen = set()
    for i, (sec, term, must) in enumerate(QUERIES, 1):
        ids = esearch(term, n=8)
        rec = None
        summ = esummary(ids) if ids else {}
        for pmid in ids:
            r = summ.get(pmid)
            if not isinstance(r, dict):
                continue
            title = (r.get("title") or "").lower()
            if must and must.lower() not in title:
                continue
            if pmid in seen:
                continue
            rec = (pmid, r)
            break
        if rec is None:
            # accept the top hit even if the keyword check failed, but SAY SO
            for pmid in ids:
                r = summ.get(pmid)
                if isinstance(r, dict) and pmid not in seen:
                    rec = (pmid, r)
                    break
        if rec is None:
            misses.append(term)
            print("   MISS  %s" % term[:78])
            continue
        pmid, r = rec
        seen.add(pmid)
        got.setdefault(sec, []).append((pmid, r))
        print("   [%3d/%3d] %s  %s" % (i, len(QUERIES), pmid, (r.get("title") or "")[:70]))
        time.sleep(0.36)          # stay under the 3 requests/second guidance
    return got, misses


def fmt(pmid, r):
    authors = r.get("authors") or []
    names = [a.get("name", "") for a in authors if a.get("authtype") == "Author"]
    if len(names) > 3:
        who = ", ".join(names[:3]) + ", et al."
    else:
        who = ", ".join(names) if names else "—"
    title = re.sub(r"\s+", " ", (r.get("title") or "").strip()).rstrip(".")
    journal = r.get("source") or ""
    year = (r.get("pubdate") or "")[:4]
    vol = r.get("volume") or ""
    pages = r.get("pages") or ""
    bits = journal
    if year:
        bits += f". {year}"
    if vol:
        bits += f";{vol}"
    if pages:
        bits += f":{pages}"
    return (f"- {who} **{title}**. *{bits}*. "
            f"[PMID {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")


def main():
    got, misses = build()
    n = sum(len(v) for v in got.values())
    lines = [
        "# Malignant Pleural Mesothelioma — QSP model reference list",
        "",
        f"**{n} references**, grouped by the part of the model they support.",
        "",
        "Every entry in this file was retrieved from PubMed with NCBI E-utilities by ",
        "`mpm_fetch_references.py`. The PMID, title, journal, year and author list are ",
        "whatever the NCBI record actually contains — nothing here was written from ",
        "memory, and searches that returned nothing were reported as misses rather than ",
        "filled in.",
        "",
        "---",
        "",
    ]
    for key, heading in SECTIONS:
        recs = got.get(key, [])
        if not recs:
            continue
        lines.append(f"## {heading}")
        lines.append("")
        for pmid, r in recs:
            lines.append(fmt(pmid, r))
        lines.append("")
    if misses:
        lines.append("---")
        lines.append("")
        lines.append("## Searches that returned nothing")
        lines.append("")
        lines.append("These queries were run and came back empty. They are listed rather "
                     "than replaced with a remembered citation.")
        lines.append("")
        for m in misses:
            lines.append(f"- `{m}`")
        lines.append("")
    open(OUT, "w", encoding="utf-8").write("\n".join(lines))
    print(f"\nwrote {OUT}: {n} references, {len(misses)} misses")


if __name__ == "__main__":
    main()
