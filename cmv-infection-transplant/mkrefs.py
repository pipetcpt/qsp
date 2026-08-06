#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build cmv_references.md from PubMed.

Every title, journal, year, author and PMID in the bibliography is resolved LIVE
from NCBI esearch + esummary.  Nothing is written from memory, because a
plausible-looking PMID written from memory is a fabricated citation.

Each entry is a (section, intent, query) triple.  The top RELEVANCE-ranked hit is
taken and printed together with its INTENT, so a retrieved paper that does not
match what it was asked for is visible at a glance instead of hidden.

E-utilities defaults to DATE order, not relevance, so sort=relevance is
mandatory: without it the "top hit" is simply the newest paper matching any term.

    python3 mkrefs.py             # build (uses refs_raw.json cache if present)
    python3 mkrefs.py --refresh   # ignore the cache and re-query
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = "refs_raw.json"
META = "refs_meta.json"
OUT = "cmv_references.md"


def _get(url):
    for k in range(6):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode(
                "utf-8", "replace")
        except Exception:
            time.sleep(2.0 + 3 * k)
    return ""


def esearch(term, n=1):
    u = B + ("esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d"
             "&term=%s" % (n, urllib.parse.quote(term)))
    try:
        return json.loads(_get(u))["esearchresult"]["idlist"]
    except Exception:
        return []


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 120):
        chunk = pmids[i:i + 120]
        u = B + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(chunk)
        try:
            r = json.loads(_get(u))["result"]
        except Exception:
            continue
        for pid in chunk:
            if pid in r:
                out[pid] = r[pid]
        time.sleep(0.34)
    return out


def author_str(authors):
    if not authors:
        return "?"
    n = len(authors)
    a = authors[0].get("name", "?")
    if n == 1:
        return a
    if n == 2:
        return "%s & %s" % (a, authors[1].get("name", "?"))
    return "%s et al." % a


def year_of(rec):
    d = rec.get("pubdate", "") or rec.get("epubdate", "")
    return (d.split(" ")[0] or "?")


# ---------------------------------------------------------------------------
#  (section, intent, query)
#  The INTENT column is the point: it records what the citation is being asked
#  to support, so that a mismatch is a visible defect rather than a silent one.
# ---------------------------------------------------------------------------
S1 = "1 · 개괄 · 역학 · 진료지침 (Overview · epidemiology · guidelines)"
S2 = "2 · 바이러스 동역학: 배가시간 · 감소반감기 (Viral kinetics in vivo)"
S3 = "3 · 잠복과 재활성화 (Latency, reactivation, myeloid reservoir)"
S4 = "4 · CMV 특이 세포면역과 방어 역치 (CMV-specific T cells and protection)"
S5 = "5 · 면역억제제와 CMV 위험 (Immunosuppression, ATG, mTOR inhibitors)"
S6 = "6 · 간시클로버 · 발간시클로버 PK/PD (Ganciclovir / valganciclovir)"
S7 = "7 · 레테르모비르 (Letermovir: terminase, PK, DDI, trials)"
S8 = "8 · 마리바비르 (Maribavir: pUL97, antagonism with ganciclovir, SOLSTICE)"
S9 = "9 · 포스카넷 · 시도포비어 · 구제요법 (Foscarnet, cidofovir, salvage)"
S10 = "10 · 내성 기전과 유전자형 (Resistance: UL97, UL54, UL56, UL27)"
S11 = "11 · 예방 vs 선제치료 · 감시간격 (Prophylaxis vs pre-emptive; monitoring)"
S12 = "12 · 후기발현 CMV와 면역유도 전략 (Late-onset CMV; immune-guided stopping)"
S13 = "13 · 골수억제 · 호중구감소 모델링 (Myelosuppression; Friberg model)"
S14 = "14 · 신독성 · 용량조절 (Nephrotoxicity and renal dose adjustment)"
S15 = "15 · 간접효과 · 이식편 결과 (Indirect effects; rejection; graft outcome)"
S16 = "16 · 장기별 질환 · 조직 구획화 (End-organ disease; compartmentalisation)"
S17 = "17 · 검사 표준화 · IU/mL (Assay standardisation; WHO standard)"
S18 = "18 · 세포치료 · 면역글로불린 · 백신 (Cell therapy, CMV-Ig, vaccines)"
S19 = "19 · QSP · 모델링 방법론 (QSP and viral-dynamics modelling methodology)"

Q = [
 # ---------------------------------------------------------------- 1
 (S1, "modern review of CMV in solid-organ transplantation",
  "cytomegalovirus solid organ transplantation review management"),
 (S1, "The Transplantation Society international consensus guidelines on CMV",
  "Transplantation Society international consensus guidelines cytomegalovirus solid organ transplant"),
 (S1, "AST Infectious Diseases Community of Practice CMV guideline",
  "American Society Transplantation Infectious Diseases Community Practice cytomegalovirus guidelines"),
 (S1, "CMV after haematopoietic cell transplantation, current management",
  "cytomegalovirus hematopoietic cell transplantation management review"),
 (S1, "epidemiology and risk factors for CMV infection after kidney transplant",
  "cytomegalovirus infection kidney transplantation risk factors epidemiology"),
 (S1, "donor/recipient serostatus as the dominant risk factor (D+/R-)",
  "cytomegalovirus donor positive recipient negative seronegative risk transplant"),
 (S1, "economic burden and cost of CMV after transplantation",
  "cytomegalovirus transplantation cost economic burden"),

 # ---------------------------------------------------------------- 2
 (S2, "in-vivo CMV replication dynamics and doubling time in transplant patients",
  "dynamics human cytomegalovirus replication in vivo doubling time"),
 (S2, "viral load kinetics and the decline half-life on ganciclovir",
  "cytomegalovirus viral load decline half-life ganciclovir kinetics"),
 (S2, "viral load doubling time as a predictor of CMV disease",
  "cytomegalovirus DNA doubling time predictor disease transplant"),
 (S2, "quantitative relationship between peak viral load and end-organ disease",
  "cytomegalovirus viral load threshold end organ disease quantitative"),
 (S2, "target-cell-limited modelling of herpesvirus replication in vivo",
  "mathematical model cytomegalovirus replication kinetics target cell"),
 (S2, "non-encapsidated versus virion-associated CMV DNA in plasma",
  "cytomegalovirus plasma DNA fragmented encapsidated"),

 # ---------------------------------------------------------------- 3
 (S3, "CMV latency in CD34+ progenitors and CD14+ monocytes",
  "cytomegalovirus latency CD34 progenitor CD14 monocyte reservoir"),
 (S3, "the major immediate-early promoter as the latency switch",
  "cytomegalovirus major immediate early promoter latency reactivation chromatin"),
 (S3, "TNF-alpha, catecholamines and cAMP-PKA signalling as reactivation triggers",
  "cytomegalovirus reactivation tumor necrosis factor catecholamine cyclic AMP"),
 (S3, "myeloid differentiation and reactivation from latency",
  "cytomegalovirus reactivation monocyte macrophage differentiation latency"),
 (S3, "allogeneic stimulation, inflammation and CMV reactivation in critical illness",
  "cytomegalovirus reactivation critically ill sepsis inflammation"),
 (S3, "UL138 and LUNA latency-associated transcripts",
  "cytomegalovirus UL138 LUNA latency associated transcript"),

 # ---------------------------------------------------------------- 4
 (S4, "CMV-specific CD8 and CD4 T cells and protection from CMV disease",
  "cytomegalovirus specific CD8 CD4 T cell protection disease transplant recipients"),
 (S4, "QuantiFERON-CMV to predict late-onset disease after prophylaxis",
  "QuantiFERON-CMV predict late onset cytomegalovirus disease transplant"),
 (S4, "ELISPOT / interferon-gamma release assay thresholds for CMV protection",
  "cytomegalovirus ELISPOT interferon gamma release assay threshold protection transplant"),
 (S4, "cell-mediated immunity guided management of CMV",
  "cytomegalovirus cell mediated immunity guided prophylaxis randomized"),
 (S4, "kinetics of CMV-specific T-cell reconstitution after transplantation",
  "cytomegalovirus specific T cell reconstitution kinetics after transplantation"),
 (S4, "adaptive NKG2C+ NK cell expansion during CMV infection",
  "NKG2C adaptive natural killer cell expansion cytomegalovirus"),
 (S4, "memory inflation and T-cell exhaustion in chronic CMV infection",
  "cytomegalovirus memory inflation T cell exhaustion PD-1"),
 (S4, "neutralising antibody to the pentamer complex and CMV neutralisation",
  "cytomegalovirus pentamer gH gL UL128 neutralizing antibody"),
 (S4, "US2-US11 mediated MHC class I down-regulation and CD8 evasion",
  "cytomegalovirus US2 US11 MHC class I downregulation immune evasion"),

 # ---------------------------------------------------------------- 5
 (S5, "antithymocyte globulin and the risk of CMV infection",
  "antithymocyte globulin induction cytomegalovirus infection risk transplant"),
 (S5, "alemtuzumab lymphodepletion and opportunistic viral infection",
  "alemtuzumab lymphocyte depletion cytomegalovirus infection transplantation"),
 (S5, "mTOR inhibitors reduce CMV infection after transplantation",
  "mTOR inhibitor everolimus sirolimus reduced cytomegalovirus infection transplant"),
 (S5, "mycophenolate and the risk of CMV disease",
  "mycophenolate mofetil cytomegalovirus disease risk transplant"),
 (S5, "net state of immunosuppression as a determinant of infection risk",
  "net state of immunosuppression infection risk transplant recipient"),
 (S5, "tacrolimus exposure, calcineurin inhibition and antiviral T-cell function",
  "tacrolimus calcineurin inhibitor antiviral T cell function cytomegalovirus"),

 # ---------------------------------------------------------------- 6
 (S6, "valganciclovir pharmacokinetics and bioavailability in transplant recipients",
  "valganciclovir pharmacokinetics bioavailability ganciclovir transplant recipients"),
 (S6, "ganciclovir population pharmacokinetics and renal clearance",
  "ganciclovir population pharmacokinetics renal clearance creatinine"),
 (S6, "ganciclovir triphosphate intracellular half-life and mechanism",
  "ganciclovir triphosphate intracellular"),
 (S6, "pUL97 phosphorylation is required to activate ganciclovir",
  "UL97 phosphotransferase ganciclovir phosphorylation activation cytomegalovirus"),
 (S6, "valganciclovir 900 mg once daily for prophylaxis: PK-PD target attainment",
  "valganciclovir exposure response cytomegalovirus prophylaxis"),
 (S6, "IMPACT trial: 200 versus 100 days of valganciclovir prophylaxis",
  "valganciclovir 200 days versus 100 days prophylaxis kidney transplant IMPACT"),
 (S6, "oral valganciclovir versus intravenous ganciclovir for CMV disease treatment",
  "valganciclovir versus intravenous ganciclovir treatment cytomegalovirus disease VICTOR"),

 # ---------------------------------------------------------------- 7
 (S7, "letermovir mechanism: inhibition of the CMV terminase complex",
  "letermovir terminase complex UL56 mechanism of action cytomegalovirus"),
 (S7, "letermovir pharmacokinetics, protein binding and dose",
  "letermovir pharmacokinetics protein binding dose healthy subjects"),
 (S7, "letermovir prophylaxis after allogeneic HSCT (phase 3)",
  "letermovir prophylaxis cytomegalovirus allogeneic hematopoietic cell transplantation phase 3"),
 (S7, "letermovir versus valganciclovir prophylaxis in kidney transplant D+/R-",
  "letermovir versus valganciclovir prophylaxis kidney transplant recipients randomized"),
 (S7, "letermovir drug interaction with tacrolimus and cyclosporine",
  "letermovir tacrolimus drug interaction"),
 (S7, "letermovir in-vitro potency EC50 against cytomegalovirus",
  "letermovir in vitro antiviral activity EC50 cytomegalovirus"),
 (S7, "extended letermovir prophylaxis to 200 days after HSCT",
  "letermovir extended prophylaxis 200 days cytomegalovirus transplantation"),
 (S7, "letermovir has no activity against other herpesviruses",
  "letermovir antiviral spectrum herpesviruses"),

 # ---------------------------------------------------------------- 8
 (S8, "maribavir mechanism: inhibition of pUL97 kinase",
  "maribavir UL97 kinase inhibitor mechanism cytomegalovirus"),
 (S8, "maribavir for refractory or resistant CMV (SOLSTICE)",
  "maribavir refractory resistant cytomegalovirus transplant randomized SOLSTICE"),
 (S8, "maribavir antagonises ganciclovir in vitro",
  "maribavir ganciclovir antagonism in vitro combination cytomegalovirus"),
 (S8, "maribavir pharmacokinetics and dysgeusia",
  "maribavir pharmacokinetics safety"),
 (S8, "maribavir resistance mutations UL97 T409M H411Y C480F",
  "maribavir resistance UL97 T409M H411Y C480F mutation"),
 (S8, "maribavir failure of the earlier prophylaxis phase 3 programme",
  "maribavir prophylaxis stem cell transplant trial"),

 # ---------------------------------------------------------------- 9
 (S9, "foscarnet mechanism, pharmacokinetics and dosing",
  "foscarnet pharmacokinetics cytomegalovirus"),
 (S9, "foscarnet nephrotoxicity and electrolyte wasting",
  "foscarnet nephrotoxicity electrolyte"),
 (S9, "cidofovir pharmacology and intracellular metabolite half-life",
  "cidofovir pharmacokinetics intracellular diphosphate half life"),
 (S9, "cidofovir proximal tubular toxicity and probenecid",
  "cidofovir nephrotoxicity proximal tubule probenecid"),
 (S9, "management of ganciclovir-resistant CMV: comparative salvage options",
  "ganciclovir resistant cytomegalovirus salvage foscarnet treatment outcome transplant"),
 (S9, "combination antiviral therapy for resistant CMV",
  "combination antiviral therapy resistant cytomegalovirus ganciclovir foscarnet"),

 # ---------------------------------------------------------------- 10
 (S10, "UL97 mutations conferring ganciclovir resistance",
  "UL97 mutation ganciclovir resistance cytomegalovirus M460V A594V"),
 (S10, "UL54 DNA polymerase mutations and cross-resistance patterns",
  "UL54 DNA polymerase mutation cytomegalovirus cross resistance cidofovir foscarnet"),
 (S10, "UL56 C325Y and letermovir resistance",
  "UL56 C325Y letermovir resistance cytomegalovirus mutation"),
 (S10, "letermovir resistance emerging in clinical use",
  "letermovir resistance clinical breakthrough UL56 mutation transplant"),
 (S10, "replicative fitness cost of antiviral resistance mutations in CMV",
  "cytomegalovirus antiviral resistance fitness"),
 (S10, "genotypic resistance testing and next-generation sequencing for CMV",
  "cytomegalovirus genotypic antiviral resistance testing next generation sequencing"),
 (S10, "incidence and risk factors for ganciclovir resistance in D+/R- recipients",
  "ganciclovir resistance incidence solid organ transplant"),
 (S10, "recombinant phenotyping to quantify fold-change in EC50",
  "cytomegalovirus recombinant phenotyping antiviral susceptibility"),

 # ---------------------------------------------------------------- 11
 (S11, "prophylaxis versus pre-emptive therapy: randomised comparison",
  "prophylaxis versus preemptive therapy cytomegalovirus randomized transplant"),
 (S11, "meta-analysis of prophylaxis versus pre-emptive therapy",
  "meta-analysis prophylaxis preemptive cytomegalovirus solid organ transplantation"),
 (S11, "PCR monitoring frequency and the viral load threshold for pre-emption",
  "cytomegalovirus preemptive therapy viral load threshold monitoring frequency"),
 (S11, "how often to monitor: weekly versus less frequent CMV PCR",
  "weekly cytomegalovirus PCR monitoring frequency preemptive interval"),
 (S11, "hybrid prophylaxis followed by surveillance strategies",
  "hybrid prophylaxis surveillance cytomegalovirus strategy transplant"),
 (S11, "cost-effectiveness of prophylaxis versus pre-emptive therapy",
  "cost effectiveness prophylaxis preemptive cytomegalovirus transplantation"),

 # ---------------------------------------------------------------- 12
 (S12, "late-onset CMV disease after discontinuing prophylaxis",
  "late onset cytomegalovirus disease after prophylaxis discontinuation transplant"),
 (S12, "prophylaxis delays but does not prevent CMV-specific immune reconstitution",
  "prophylaxis delayed cytomegalovirus specific immunity reconstitution transplant"),
 (S12, "subclinical replication during prophylaxis and T-cell priming",
  "cytomegalovirus specific immunity during valganciclovir prophylaxis"),
 (S12, "immune-guided duration of prophylaxis: trials and feasibility",
  "immune guided duration prophylaxis cytomegalovirus cell mediated immunity trial"),
 (S12, "relapse and secondary prophylaxis after treating a CMV episode",
  "cytomegalovirus recurrence relapse secondary prophylaxis after treatment"),

 # ---------------------------------------------------------------- 13
 (S13, "Friberg semi-mechanistic model of chemotherapy-induced myelosuppression",
  "Friberg semimechanistic model myelosuppression neutropenia chemotherapy"),
 (S13, "valganciclovir-associated neutropenia incidence and management",
  "valganciclovir neutropenia incidence transplant recipients management"),
 (S13, "ganciclovir myelotoxicity mechanism",
  "ganciclovir bone marrow toxicity"),
 (S13, "G-CSF for antiviral-associated neutropenia in transplant recipients",
  "granulocyte colony stimulating factor neutropenia transplant valganciclovir"),
 (S13, "letermovir has less myelosuppression than valganciclovir",
  "letermovir versus valganciclovir myelosuppression leukopenia neutropenia"),

 # ---------------------------------------------------------------- 14
 (S14, "renal dose adjustment of valganciclovir and the risk of under-dosing",
  "valganciclovir renal dose adjustment underdosing overdosing transplant"),
 (S14, "ganciclovir overexposure in renal impairment and toxicity",
  "ganciclovir renal impairment dose"),
 (S14, "therapeutic drug monitoring of ganciclovir",
  "ganciclovir therapeutic drug monitoring target exposure"),
 (S14, "CMV and allograft function: eGFR trajectory after infection",
  "cytomegalovirus allograft function glomerular filtration rate kidney transplant"),

 # ---------------------------------------------------------------- 15
 (S15, "indirect effects of CMV: rejection, graft loss and mortality",
  "cytomegalovirus indirect effects acute rejection graft loss mortality transplant"),
 (S15, "CMV and chronic lung allograft dysfunction / bronchiolitis obliterans",
  "cytomegalovirus chronic lung allograft dysfunction bronchiolitis obliterans"),
 (S15, "CMV and cardiac allograft vasculopathy",
  "cytomegalovirus cardiac allograft vasculopathy heart transplant"),
 (S15, "CMV endothelial infection and pro-inflammatory cytokine induction",
  "cytomegalovirus endothelial cell infection inflammatory cytokine interleukin 6"),
 (S15, "CMV and new-onset diabetes after transplantation",
  "cytomegalovirus new onset diabetes after transplantation"),
 (S15, "CMV as an immunomodulator predisposing to other opportunistic infection",
  "cytomegalovirus immunomodulation opportunistic infection transplant recipient"),

 # ---------------------------------------------------------------- 16
 (S16, "CMV gastrointestinal disease in transplant recipients",
  "cytomegalovirus gastrointestinal colitis disease transplant recipients"),
 (S16, "compartmentalised CMV disease with negative plasma PCR",
  "cytomegalovirus gastrointestinal disease negative viremia"),
 (S16, "CMV pneumonitis after transplantation",
  "cytomegalovirus pneumonia pneumonitis transplantation outcome"),
 (S16, "CMV retinitis in transplant recipients and drug penetration",
  "cytomegalovirus retinitis ganciclovir intraocular"),
 (S16, "CMV central nervous system disease and CSF drug penetration",
  "cytomegalovirus encephalitis central nervous system cerebrospinal fluid ganciclovir"),
 (S16, "definitions of CMV infection and disease for clinical trials",
  "definitions cytomegalovirus infection disease clinical trials consensus"),

 # ---------------------------------------------------------------- 17
 (S17, "WHO international standard for CMV DNA and IU/mL harmonisation",
  "WHO international standard cytomegalovirus DNA international units harmonization"),
 (S17, "inter-laboratory variability in CMV viral load quantification",
  "interlaboratory variability cytomegalovirus viral load quantification"),
 (S17, "plasma versus whole blood for CMV DNA quantification",
  "plasma versus whole blood cytomegalovirus DNA quantification comparison"),
 (S17, "thresholds for initiating pre-emptive therapy are assay dependent",
  "cytomegalovirus viral load threshold assay dependent preemptive therapy"),

 # ---------------------------------------------------------------- 18
 (S18, "adoptive transfer of CMV-specific T cells",
  "adoptive transfer cytomegalovirus specific T cells transplant"),
 (S18, "third-party virus-specific T cells for refractory CMV",
  "third party virus specific T cells refractory cytomegalovirus"),
 (S18, "CMV hyperimmune globulin in transplantation",
  "cytomegalovirus hyperimmune globulin transplantation prophylaxis"),
 (S18, "CMV vaccine candidates including mRNA-1647",
  "cytomegalovirus vaccine candidate mRNA-1647 glycoprotein B clinical trial"),

 # ---------------------------------------------------------------- 19
 (S19, "quantitative systems pharmacology: scope and best practice",
  "quantitative systems pharmacology model best practice drug development"),
 (S19, "mrgsolve for ODE-based pharmacometric simulation",
  "mrgsolve simulation pharmacometric ordinary differential equation R"),
 (S19, "viral dynamics modelling to design antiviral trials",
  "viral dynamics model antiviral efficacy"),
 (S19, "the critical efficacy threshold for antiviral suppression",
  "antiviral efficacy threshold critical drug efficacy viral dynamics model"),
 (S19, "modelling emergence of antiviral drug resistance",
  "mathematical model emergence antiviral drug resistance mutant selection window"),
 (S19, "target-cell-limited versus immune-mediated control in viral infection models",
  "target cell limited immune mediated control viral infection model"),
 (S19, "identifiability and calibration of QSP models",
  "identifiability parameter estimation quantitative systems pharmacology calibration"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()

    cache = {} if a.refresh else (
        json.load(open(CACHE)) if os.path.exists(CACHE) else {})
    resolved = OrderedDict()
    for i, (sec, intent, q) in enumerate(Q):
        key = q
        if key not in cache:
            cache[key] = esearch(q, 1)
            time.sleep(0.40)
            sys.stderr.write("[%3d/%3d] %s -> %s\n"
                             % (i + 1, len(Q), q[:58],
                                cache[key][0] if cache[key] else "MISS"))
        resolved[(sec, intent, q)] = cache[key][0] if cache[key] else None
    json.dump(cache, open(CACHE, "w"), indent=1)

    pmids = sorted({v for v in resolved.values() if v})
    meta = json.load(open(META)) if (os.path.exists(META) and not a.refresh) else {}
    need = [p for p in pmids if p not in meta]
    if need:
        meta.update(esummary(need))
        json.dump(meta, open(META, "w"), indent=1)

    secs = OrderedDict()
    seen = set()
    dupes = []
    for (sec, intent, q), pid in resolved.items():
        secs.setdefault(sec, [])
        if pid is None:
            secs[sec].append((intent, None, None))
            continue
        if pid in seen:
            dupes.append((sec, intent, pid))
        seen.add(pid)
        secs[sec].append((intent, pid, meta.get(pid, {})))

    n = sum(1 for s in secs.values() for e in s if e[1])
    L = []
    L.append("# 이식 후 거대세포바이러스(CMV) 감염 — 참고문헌")
    L.append("### Cytomegalovirus in Transplant Recipients · References")
    L.append("")
    L.append("**총 %d편** · %d개 섹션 · 모든 PMID·제목·저널·연도는 `mkrefs.py`가 "
             "NCBI E-utilities로 **실시간 조회**한 것이며 기억으로 작성된 것이 "
             "없습니다." % (n, len(secs)))
    L.append("")
    L.append("각 항목의 `의도(intent)`는 그 인용이 **무엇을 뒷받침하도록 "
             "요청되었는지**를 기록합니다. 검색 결과가 의도와 맞지 않으면 그 "
             "불일치가 숨지 않고 한눈에 보이도록 하기 위한 것입니다.")
    L.append("")
    L.append("재생성:")
    L.append("```bash")
    L.append("python3 mkrefs.py --refresh    # PubMed 재조회")
    L.append("```")
    L.append("")
    L.append("---")
    L.append("")
    for sec, entries in secs.items():
        L.append("## %s" % sec)
        L.append("")
        for intent, pid, rec in entries:
            if not pid:
                L.append("- *(unresolved query — no PubMed hit)* — 의도: %s" % intent)
                continue
            title = (rec.get("title") or "").strip().rstrip(".")
            src = rec.get("source", "")
            L.append("- **%s** (%s, %s). %s  "
                     % (author_str(rec.get("authors")), src, year_of(rec), title))
            L.append("  [PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/) · "
                     "의도: *%s*" % (pid, pid, intent))
        L.append("")
    if dupes:
        L.append("---")
        L.append("")
        L.append("## 중복 조회 기록 (queries that resolved to an already-cited paper)")
        L.append("")
        L.append("서로 다른 의도의 검색이 같은 논문으로 수렴한 경우입니다. "
                 "숨기지 않고 남겨 둡니다.")
        L.append("")
        for sec, intent, pid in dupes:
            L.append("- [PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/) — %s / *%s*"
                     % (pid, pid, sec.split(" · ")[0], intent))
        L.append("")
    L.append("---")
    L.append("")
    L.append("## 인용 위치 (where these are used in the model)")
    L.append("")
    L.append("| 모델 구성요소 | 근거 섹션 |")
    L.append("|---|---|")
    L.append("| `DBL0` 1.2 d · `THALFTX` 2.4 d → `KPROD`, e\\* = 0.667 | 2 |")
    L.append("| 잠복 저장소 · `KREACT`·`AREACT` 재활성화 항 | 3 |")
    L.append("| `KE8`·`E8MAX`·`KMEM` → E8\\* = 4.81 /µL (역치 2) | 4 |")
    L.append("| `KAG` 항원 의존 증식 = 후기발현 CMV의 기전 | 4, 12 |")
    L.append("| `ISI` 가중치 · ATG·mTORi 항 | 5 |")
    L.append("| `FVGCV`·`CLGCV`·`FRENGCV`·`KPHOS`·`EC50GTP` | 6 |")
    L.append("| `FULTV`·`EC50LTV`·`DDITLTV` 2.4배 상호작용 | 7 |")
    L.append("| `KIMBVK` (pUL97 활성화 차단 = 길항작용) | 8 |")
    L.append("| `EC50FOS`·`AFOS`·`FMG` · 시도포비어 `KPPCDV` | 9, 14 |")
    L.append("| `FK_A` 0.125 · `RESBLTV` 3000 · `FIT_A`·`FIT_B` | 10 |")
    L.append("| 감시간격 시계 2^(Δt/1.2 d) · 역치 설정 | 11, 17 |")
    L.append("| `GAM`·`EMAXMYE`·`EMPA`·`KGCSF` (Friberg 사슬) | 13 |")
    L.append("| `H0REJ`·`CREJV`·`CREJIS` · eGFR 궤적 | 15 |")
    L.append("| `KDIS`·`KD50`·`PVIRT` · sanctuary 침투계수 | 16 |")
    L.append("| `ACT`·`CMVIG` 입력 | 18 |")
    L.append("| 모델 구조 · 임계효능 개념 · 검증 방법론 | 19 |")
    L.append("")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("wrote %s: %d references, %d sections, %d duplicate resolutions"
          % (OUT, n, len(secs), len(dupes)))
    miss = [(s, i) for s, e in secs.items() for i, p, _ in e if not p]
    if miss:
        print("UNRESOLVED:")
        for s, i in miss:
            print("  ", s, "|", i)


if __name__ == "__main__":
    main()
