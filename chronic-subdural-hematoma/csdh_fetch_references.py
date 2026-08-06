"""Programmatically retrieve real PubMed records for the cSDH reference list."""
import json, time, urllib.request, urllib.parse, sys

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

SECTIONS = [
 ("역학·자연사 (Epidemiology & Natural History)", [
    ("chronic subdural hematoma epidemiology incidence ageing population", 4),
    ("chronic subdural hematoma natural history spontaneous resolution conservative", 4),
    ("chronic subdural hematoma recurrence risk factors predictors", 4),
 ]),
 ("병태생리: 신생막과 신생혈관 (Neomembrane & Neovascularisation)", [
    ("chronic subdural hematoma outer membrane neomembrane histology", 4),
    ("chronic subdural hematoma vascular endothelial growth factor VEGF", 5),
    ("chronic subdural hematoma angiopoietin Tie2 vascular permeability", 3),
    ("chronic subdural hematoma neovascularization immature capillary fenestration", 3),
 ]),
 ("병태생리: 염증과 헴 (Inflammation & Haem)", [
    ("chronic subdural hematoma inflammation interleukin-6 interleukin-8", 5),
    ("chronic subdural hematoma eosinophil macrophage infiltration membrane", 3),
    ("chronic subdural hematoma matrix metalloproteinase MMP-9", 4),
    ("subdural hematoma heme oxygenase iron oxidative stress", 3),
 ]),
 ("병태생리: 국소 과섬유소용해 (Local Hyperfibrinolysis)", [
    ("chronic subdural hematoma fibrinolysis tissue plasminogen activator", 4),
    ("chronic subdural hematoma fibrin degradation products D-dimer fluid", 4),
    ("chronic subdural hematoma coagulopathy haemostasis membrane bleeding", 3),
 ]),
 ("체액 역학과 청소 경로 (Fluid Dynamics & Clearance)", [
    ("chronic subdural hematoma osmotic theory protein content fluid", 4),
    ("meningeal lymphatic vessels cerebrospinal fluid drainage dura", 4),
    ("meningeal lymphatics subdural hematoma clearance", 3),
    ("chronic subdural hematoma computed tomography density hounsfield classification", 3),
 ]),
 ("두개내 역학·뇌 재팽창 (Intracranial Mechanics & Brain Re-expansion)", [
    ("chronic subdural hematoma brain re-expansion postoperative residual space", 4),
    ("chronic subdural hematoma midline shift brain atrophy volume", 4),
    ("Monro-Kellie intracranial compliance pressure volume", 3),
 ]),
 ("수술적 치료 (Surgical Management)", [
    ("burr hole craniostomy chronic subdural hematoma randomized", 4),
    ("subdural drain chronic subdural haematoma randomised controlled trial Santarius", 4),
    ("chronic subdural hematoma irrigation versus no irrigation trial", 3),
    ("chronic subdural hematoma craniotomy membranectomy septated recurrence", 3),
 ]),
 ("중경막동맥 색전술 (Middle Meningeal Artery Embolisation)", [
    ("middle meningeal artery embolization chronic subdural hematoma randomized trial", 6),
    ("EMBOLISE middle meningeal artery embolization subdural hematoma", 3),
    ("MAGIC-MT middle meningeal artery embolization trial subdural", 3),
    ("middle meningeal artery embolization technique anatomy dangerous anastomoses", 3),
 ]),
 ("스테로이드 (Corticosteroids)", [
    ("dexamethasone chronic subdural haematoma randomised trial", 5),
    ("corticosteroid chronic subdural hematoma outcome meta-analysis", 3),
    ("dexamethasone adverse effects hyperglycaemia infection elderly", 3),
 ]),
 ("아토르바스타틴 (Atorvastatin)", [
    ("atorvastatin chronic subdural hematoma randomized clinical trial", 5),
    ("atorvastatin dexamethasone chronic subdural haematoma ATOCH", 3),
    ("statin angiogenesis vascular maturation endothelial progenitor pleiotropic", 3),
 ]),
 ("트라넥삼산·항혈전제 (Tranexamic Acid & Antithrombotics)", [
    ("tranexamic acid chronic subdural hematoma trial", 5),
    ("anticoagulation resumption chronic subdural hematoma timing", 4),
    ("antiplatelet anticoagulant chronic subdural haematoma recurrence risk", 3),
 ]),
 ("약동학 (Pharmacokinetics)", [
    ("dexamethasone pharmacokinetics population adults", 3),
    ("atorvastatin pharmacokinetics bioavailability metabolites", 3),
    ("tranexamic acid pharmacokinetics population model", 3),
    ("apixaban pharmacokinetics population elderly renal", 3),
 ]),
 ("결과·예후 (Outcomes & Prognosis)", [
    ("chronic subdural hematoma modified Rankin Scale functional outcome elderly", 4),
    ("chronic subdural hematoma mortality long term survival", 3),
    ("chronic subdural hematoma cognitive outcome dementia", 3),
 ]),
 ("QSP 방법론 (QSP Methodology)", [
    ("quantitative systems pharmacology model development validation", 3),
    ("mrgsolve simulation pharmacometrics R", 2),
    ("physiologically based pharmacokinetic model verification good practice", 2),
 ]),
]


def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode()
        except Exception as e:
            if attempt == 3:
                print("   ! failed:", e, file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))


def esearch(term, n):
    q = urllib.parse.urlencode(dict(db="pubmed", term=term, retmax=n,
                                    retmode="json", sort="relevance"))
    r = get(E + "esearch.fcgi?" + q)
    if not r: return []
    try:
        return json.loads(r)["esearchresult"].get("idlist", [])
    except Exception:
        return []


def esummary(pmids):
    if not pmids: return {}
    q = urllib.parse.urlencode(dict(db="pubmed", id=",".join(pmids),
                                    retmode="json"))
    r = get(E + "esummary.fcgi?" + q)
    if not r: return {}
    try:
        return json.loads(r)["result"]
    except Exception:
        return {}


if __name__ == "__main__":
    seen, out = set(), []
    for sect, queries in SECTIONS:
        recs = []
        for term, n in queries:
            ids = [i for i in esearch(term, n + 3) if i not in seen]
            ids = ids[:n]
            seen.update(ids)
            res = esummary(ids)
            for pid in ids:
                d = res.get(pid)
                if not isinstance(d, dict) or 'title' not in d:
                    continue
                authors = d.get('authors', [])
                a1 = authors[0]['name'] if authors else '?'
                al = f"{a1} et al." if len(authors) > 1 else a1
                recs.append(dict(pmid=pid, title=d['title'].rstrip('.'),
                                 authors=al,
                                 journal=d.get('source', ''),
                                 year=(d.get('pubdate', '') or '')[:4]))
            time.sleep(0.4)
        out.append((sect, recs))
        print(f"{sect}: {len(recs)}", flush=True)
    json.dump(out, open("_refs.json", "w"), ensure_ascii=False, indent=1)
    print("total:", sum(len(r) for _, r in out))
