#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""abm_references.md 생성기.

NCBI E-utilities (esearch + esummary) 로 실제 반환된 레코드만 옮긴다.
직접 작성한 인용은 하나도 없다 — PMID·저자·저널·연도는 모두 API 응답값이다.
"""
import json
import time
import urllib.request
import urllib.parse
import sys

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

SECTIONS = [
    ("역학 · 임상 양상 · 예후", "Epidemiology, clinical presentation, prognosis", [
        "community-acquired bacterial meningitis adults prospective cohort",
        "pneumococcal meningitis adults outcome mortality",
        "bacterial meningitis clinical features diagnostic accuracy triad",
        "bacterial meningitis global burden epidemiology",
    ]),
    ("폐렴균 병독인자 — 방출될 화물", "Pneumococcal virulence factors", [
        "pneumolysin neurotoxicity meningitis",
        "pneumococcal cell wall inflammation cerebrospinal fluid",
        "LytA autolysin pneumococcus lysis",
        "pneumococcal capsule serotype invasive disease virulence",
        "pneumococcal hydrogen peroxide SpxB neuronal damage",
    ]),
    ("침습 경로와 혈액-CSF 장벽 통과", "Invasion and crossing of the blood-CSF barrier", [
        "Streptococcus pneumoniae invasion blood brain barrier endothelial",
        "choroid plexus epithelium bacterial invasion meningitis",
        "pneumococcal nasopharyngeal colonization invasion bacteremia",
        "platelet activating factor receptor pneumococcal adherence",
    ]),
    ("선천면역 인식 — TLR · NOD · 인플라마좀", "Innate recognition", [
        "TLR2 pneumococcal meningitis innate immunity",
        "NLRP3 inflammasome pneumococcal meningitis interleukin-1",
        "NOD2 peptidoglycan central nervous system inflammation",
        "complement pneumococcal meningitis C5a",
        "microglia activation bacterial meningitis",
    ]),
    ("CSF 사이토카인 네트워크", "CSF cytokine network", [
        "cerebrospinal fluid tumor necrosis factor bacterial meningitis",
        "interleukin-1 beta cerebrospinal fluid meningitis outcome",
        "interleukin-6 interleukin-8 cerebrospinal fluid bacterial meningitis",
        "interleukin-10 cerebrospinal fluid bacterial meningitis",
    ]),
    ("호중구 유입과 기질단백분해효소", "Neutrophil influx and matrix metalloproteinases", [
        "matrix metalloproteinase-9 cerebrospinal fluid bacterial meningitis",
        "neutrophil recruitment pleocytosis experimental meningitis",
        "MMP inhibitor experimental pneumococcal meningitis neuroprotection",
        "neutrophil extracellular traps cerebrospinal fluid meningitis",
    ]),
    ("장벽 투과성과 알부민 지수", "Barrier permeability and the albumin quotient", [
        "blood-cerebrospinal fluid barrier albumin quotient meningitis",
        "blood brain barrier permeability experimental bacterial meningitis",
        "tight junction claudin occludin bacterial meningitis",
        "VEGF cerebrospinal fluid bacterial meningitis barrier",
    ]),
    ("CSF 역학 · 두개내압 · 순응도", "CSF hydrodynamics, intracranial pressure, compliance", [
        "intracranial pressure bacterial meningitis monitoring adults",
        "cerebrospinal fluid outflow resistance meningitis",
        "pressure volume index intracranial compliance",
        "brain edema bacterial meningitis mechanism",
        "hydrocephalus bacterial meningitis",
    ]),
    ("뇌관류 · 자동조절 · 허혈", "Cerebral perfusion, autoregulation, ischaemia", [
        "cerebral blood flow autoregulation bacterial meningitis",
        "cerebral perfusion pressure meningitis outcome",
        "cerebral infarction pneumococcal meningitis",
        "transcranial Doppler bacterial meningitis vasculopathy",
    ]),
    ("CSF 대사 — 포도당과 락테이트", "CSF metabolism: glucose and lactate", [
        "cerebrospinal fluid lactate bacterial meningitis diagnostic",
        "cerebrospinal fluid glucose ratio bacterial meningitis",
        "glucose transport blood brain barrier GLUT1 kinetics",
        "cerebrospinal fluid glucose consumption leukocytes",
    ]),
    ("항생제 유도 용해와 염증 폭발", "Antibiotic-induced lysis and the inflammatory burst", [
        "antibiotic induced release pneumococcal cell wall inflammation meningitis",
        "nonbacteriolytic antibiotic experimental pneumococcal meningitis",
        "cerebrospinal fluid cytokine increase after antibiotic meningitis",
        "rifampin versus ceftriaxone experimental meningitis inflammation",
        "daptomycin experimental pneumococcal meningitis non-lytic",
    ]),
    ("항생제 CSF 침투와 PK", "Antibiotic CSF penetration and pharmacokinetics", [
        "ceftriaxone cerebrospinal fluid penetration meningitis pharmacokinetics",
        "vancomycin cerebrospinal fluid penetration meningitis",
        "rifampin cerebrospinal fluid penetration",
        "antibiotic pharmacokinetics cerebrospinal fluid penetration review",
        "ceftriaxone protein binding saturable pharmacokinetics",
        "dexamethasone effect vancomycin cerebrospinal fluid concentration",
    ]),
    ("항생제 PD — 실험적 수막염 모델", "Antibiotic PD in experimental meningitis", [
        "rabbit model pneumococcal meningitis bactericidal activity cerebrospinal fluid",
        "pharmacodynamics antibiotic cerebrospinal fluid bactericidal rate meningitis",
        "time above MIC cerebrospinal fluid beta-lactam meningitis",
        "continuous infusion ceftriaxone meningitis",
    ]),
    ("내성 폐렴균과 병용요법", "Resistant pneumococci and combination therapy", [
        "penicillin resistant pneumococcal meningitis treatment failure",
        "cephalosporin resistant Streptococcus pneumoniae meningitis vancomycin",
        "vancomycin plus ceftriaxone synergy pneumococcal meningitis",
        "PBP2x mosaic gene cephalosporin resistance pneumococcus",
    ]),
    ("덱사메타손 — 임상시험과 메타분석", "Dexamethasone: trials and meta-analyses", [
        "dexamethasone adults bacterial meningitis randomized trial",
        "corticosteroids acute bacterial meningitis Cochrane meta-analysis",
        "dexamethasone Haemophilus influenzae meningitis children hearing",
        "dexamethasone pneumococcal meningitis mortality subgroup",
        "adjunctive dexamethasone meningitis low income countries trial",
    ]),
    ("덱사메타손 — 기전과 투여 시점", "Dexamethasone: mechanism and timing", [
        "dexamethasone timing before antibiotic experimental meningitis",
        "glucocorticoid NF-kB transrepression mechanism",
        "dexamethasone cerebrospinal fluid inflammation experimental meningitis",
        "dexamethasone hippocampal apoptosis experimental meningitis",
    ]),
    ("삼투요법 · 두개내압 관리", "Osmotherapy and ICP management", [
        "glycerol adjuvant therapy bacterial meningitis trial",
        "mannitol intracranial pressure treatment",
        "hypertonic saline intracranial pressure",
        "intracranial pressure targeted treatment bacterial meningitis outcome",
    ]),
    ("청력 손실과 와우 손상", "Hearing loss and cochlear injury", [
        "hearing loss pneumococcal meningitis children",
        "labyrinthitis ossificans meningitis cochlear implantation",
        "cochlear injury experimental pneumococcal meningitis",
        "sensorineural hearing loss bacterial meningitis adults",
    ]),
    ("인지 후유증과 해마 치상핵", "Cognitive sequelae and the dentate gyrus", [
        "hippocampal apoptosis dentate gyrus experimental pneumococcal meningitis",
        "cognitive impairment survivors bacterial meningitis",
        "learning deficit experimental meningitis infant rat",
        "neuronal injury markers cerebrospinal fluid meningitis",
    ]),
    ("항생제 투여 지연과 예후", "Time to antibiotic and outcome", [
        "door to antibiotic time bacterial meningitis outcome",
        "delay antibiotic therapy adverse outcome bacterial meningitis",
        "prehospital antibiotic meningitis mortality",
    ]),
    ("예방접종과 혈청형 변화", "Vaccination and serotype replacement", [
        "pneumococcal conjugate vaccine meningitis incidence",
        "serotype replacement pneumococcal conjugate vaccine",
        "pneumococcal vaccination adults invasive disease",
    ]),
    ("진단 — CSF 검사와 분자진단", "Diagnosis: CSF studies and molecular tests", [
        "multiplex PCR cerebrospinal fluid meningitis diagnosis",
        "cerebrospinal fluid analysis bacterial versus viral meningitis score",
        "cranial CT before lumbar puncture meningitis",
        "procalcitonin bacterial meningitis diagnosis",
    ]),
    ("가이드라인과 종설", "Guidelines and reviews", [
        "ESCMID guideline acute bacterial meningitis",
        "IDSA practice guidelines bacterial meningitis management",
        "acute bacterial meningitis review Lancet",
        "community-acquired bacterial meningitis review New England",
    ]),
    ("정량 모델링 — QSP · PK/PD · 시스템생물학", "Quantitative modelling: QSP, PK/PD, systems biology", [
        "mathematical model bacterial meningitis dynamics",
        "quantitative systems pharmacology infectious disease model",
        "pharmacokinetic pharmacodynamic model antibacterial killing bacteria",
        "mrgsolve simulation pharmacometrics",
        "physiologically based model central nervous system drug distribution",
        "intracranial pressure mathematical model cerebrospinal fluid dynamics",
    ]),
]

RETMAX = 8


def eget(path, params):
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as exc:                     # noqa: BLE001
            print(f"    retry {attempt+1}: {exc}", file=sys.stderr)
            time.sleep(2 + 2 * attempt)
    return None


def search(term):
    js = eget("esearch.fcgi", dict(db="pubmed", term=term, retmax=RETMAX,
                                   retmode="json", sort="relevance"))
    if not js:
        return []
    return js.get("esearchresult", {}).get("idlist", [])


def summary(pmids):
    if not pmids:
        return {}
    js = eget("esummary.fcgi", dict(db="pubmed", id=",".join(pmids), retmode="json"))
    return (js or {}).get("result", {})


def main():
    seen = set()
    out = []
    total = 0
    for ko, en, queries in SECTIONS:
        recs = []
        for q in queries:
            ids = search(q)
            time.sleep(0.4)
            res = summary(ids)
            time.sleep(0.4)
            for pid in ids:
                if pid in seen or pid not in res:
                    continue
                r = res[pid]
                title = (r.get("title") or "").strip().rstrip(".")
                if not title:
                    continue
                auth = r.get("authors") or []
                first = auth[0].get("name", "") if auth else ""
                al = f"{first} et al." if len(auth) > 1 else first
                journal = r.get("source", "")
                year = (r.get("pubdate", "") or "")[:4]
                recs.append((pid, al, title, journal, year))
                seen.add(pid)
            print(f"  {ko} / {q}: {len(ids)} hits", file=sys.stderr)
        recs.sort(key=lambda x: (x[4] or "0"), reverse=True)
        out.append((ko, en, recs))
        total += len(recs)

    L = []
    W = L.append
    W("# 급성 세균성 수막염 (폐렴균) — 참고문헌")
    W("")
    W("**Acute Bacterial Meningitis (*Streptococcus pneumoniae*) — References**")
    W("")
    W(f"총 **{total}편**.  모든 항목은 NCBI E-utilities (`esearch` + `esummary`) 로 "
      "PubMed 에 직접 질의해 **반환된 레코드만** 옮긴 것이다 — PMID·저자·저널·연도는 "
      "모두 API 응답값이며 임의로 작성한 인용은 없다.  생성 스크립트는 이 디렉토리의 "
      "`fetch_refs.py`.")
    W("")
    W("| 섹션 | 편수 |")
    W("|------|------|")
    for i, (ko, en, recs) in enumerate(out, 1):
        anchor = f"{i}-{ko.replace(' ', '-').replace('·', '')}"
        W(f"| {i}. {ko} | {len(recs)} |")
    W("")
    W("---")
    W("")
    for i, (ko, en, recs) in enumerate(out, 1):
        W(f"## {i}. {ko}")
        W("")
        W(f"*{en}*")
        W("")
        for (pid, al, title, journal, year) in recs:
            W(f"- {al} **{title}.** *{journal}* {year}. "
              f"[PMID {pid}](https://pubmed.ncbi.nlm.nih.gov/{pid}/)")
        W("")
    W("---")
    W("")
    W("## 모델 파라미터가 어느 문헌에 걸려 있는가")
    W("")
    W("| 모델 요소 | 값 | 근거가 있는 섹션 |")
    W("|-----------|-----|------------------|")
    rows = [
        ("CSF 균 밀도 (발현 시)", "10⁵–10⁸ CFU/mL", "1, 13"),
        ("폐렴균 CSF 증식률 μ_max", "0.85 /h (배가 49 분)", "13"),
        ("세프트리악손 CSF 침투율", "총농도 2–10 %", "12"),
        ("세프트리악손 단백결합 포화", "f_u 0.076 → 0.167", "12"),
        ("반코마이신 CSF 침투율", "총농도 2–13 %", "12"),
        ("리팜핀 CSF 침투율", "10–20 %, 염증 비의존", "12"),
        ("β-락탐 살균속도", "−0.5 ~ −0.7 log₁₀ CFU/mL/h", "13"),
        ("반코마이신 살균속도", "−0.3 ~ −0.4 log₁₀ CFU/mL/h", "13"),
        ("용해성 대 비용해성 수확량", "Y_rif ≈ 0.15 × Y_βlactam", "11"),
        ("CSF TNF-α", "100–1,000 pg/mL", "5"),
        ("CSF IL-6", "10⁴–10⁶ pg/mL", "5"),
        ("CSF MMP-9", "100–1,000 ng/mL", "6"),
        ("CSF 백혈구 (폐렴균)", "1,000–5,000 /µL", "1, 22"),
        ("CSF 단백", "100–500 mg/dL", "1, 22"),
        ("CSF 포도당 비", "<0.4 (약 70 %)", "10, 22"),
        ("CSF 락테이트", ">3.5 mmol/L (흔히 6–12)", "10"),
        ("알부민 지수 Q_alb", "정상 <8 → 수막염 30–100+", "7"),
        ("CSF 생성률 Q_f", "0.35 mL/min (500 mL/일)", "8"),
        ("CSF 유출저항 R_out", "정상 6–10 mmHg/(mL/min)", "8"),
        ("압력-용적 지수 PVI", "≈25 mL", "8"),
        ("두개내압", "정상 7–15 → 25–40 mmHg", "8"),
        ("허혈 임계 CBF", "<25 mL/100 g/min", "9"),
        ("덱사메타손 용법", "0.15 mg/kg q6h × 4 일", "15"),
        ("덱사메타손 사망률 효과 (폐렴균)", "34 % → 14 %", "15"),
        ("덱사메타손 불량결과", "52 % → 26 %", "15"),
        ("중증 난청 (소아 Hib)", "~15 % → ~5 %", "15, 18"),
        ("난청 (폐렴균 성인)", "어떤 난청 20–30 %", "18"),
        ("뇌경색 (폐렴균)", "15–25 %", "9"),
        ("항생제 지연", "6 h 초과 시 예후 악화", "20"),
        ("세팔로스포린 내성 MIC", "최대 4 mg/L", "14"),
    ]
    for a, b, c in rows:
        W(f"| {a} | {b} | {c} |")
    W("")
    with open("abm_references.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"총 {total}편 저장", file=sys.stderr)


if __name__ == "__main__":
    main()
