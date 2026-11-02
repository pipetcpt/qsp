"""Study 1 data layer: mine the QSP disease-model library.

Walks every disease directory, counts ODE compartments, mechanistic-map
nodes/clusters, and references, classifies each model into a therapeutic
area from the README gallery table, and scores each model for
safety/toxicology relevance by keyword. Also reconstructs the growth
curve from git history. Writes JSON for the plotting layer.
"""
import json, re, subprocess, os
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path("/home/user/qsp")
OUT = ROOT / "presentation/2026-tox/data"
OUT.mkdir(parents=True, exist_ok=True)

# ---- README gallery: number, category, dir slug ---------------------------
readme = (ROOT / "README.md").read_text(encoding="utf-8")
cat_by_dir = {}
row_re = re.compile(r"^\|\s*(\d+)\s*\|\s*([^|]+?)\s*\|.*?\]\(([a-z0-9\-]+)/\)", re.M)
for m in row_re.finditer(readme):
    num, cat, slug = m.group(1), m.group(2).strip(), m.group(3).strip()
    cat_by_dir[slug] = cat
print("README rows:", len(cat_by_dir))

# published category counts line (ground truth for the headline bar chart)
cat_counts = {}
mline = re.search(r"\*\*분류별 모델 수\*\*:\s*(.+)", readme)
if mline:
    for part in mline.group(1).split("·"):
        part = part.strip()
        mm = re.match(r"(.+?)\s+(\d+)$", part)
        if mm:
            cat_counts[mm.group(1).strip()] = int(mm.group(2))
print("categories:", len(cat_counts), "sum", sum(cat_counts.values()))

# ---- safety / tox keyword lexicon -----------------------------------------
SAFETY_KW = [
    "toxicity", "toxic", "adverse", "safety", "hepatotox", "nephrotox",
    "cardiotox", "QT", "torsad", "myelosupp", "neutropenia", "cytokine release",
    "immune-related", "irAE", "bleeding", "hemorrhage", "haemorrhage",
    "overdose", "poisoning", "drug-induced", "DILI", "AKI", "ototox",
    "neurotox", "bone loss", "hypoglyc", "lactic acidosis", "serotonin syndrome",
    "infection risk", "PML", "anaphyl", "teratogen", "withdrawal",
]
safe_re = re.compile("|".join(re.escape(k) for k in SAFETY_KW), re.I)

def count_odes(rtext):
    # prefer explicit dxdt_ assignments; fall back to $CMT token count
    n = len(set(re.findall(r"dxdt_([A-Za-z0-9_]+)\s*=", rtext)))
    if n:
        return n
    m = re.search(r"\$CMT\b(.*?)(\$[A-Z]|$)", rtext, re.S)
    if m:
        toks = re.findall(r"[A-Za-z][A-Za-z0-9_]*", m.group(1))
        toks = [t for t in toks if t.upper() != t or True]
        return len(toks)
    return 0

def map_stats(dot):
    # nodes: identifiers that have a [ ... ] attribute list and are not edges/graph kw
    clusters = len(re.findall(r"subgraph\s+cluster", dot))
    nodes = 0
    for line in dot.splitlines():
        s = line.strip()
        if "->" in s or s.startswith(("//", "#")):
            continue
        if "[" in s and "]" in s and "=" in s.split("[", 1)[1]:
            head = s.split("[", 1)[0].strip()
            if head and not head.startswith(("graph", "node", "edge", "subgraph", "rankdir", "label", "{", "}")):
                nodes += 1
    return nodes, clusters

models = []
for d in sorted(ROOT.iterdir()):
    if not d.is_dir() or d.name.startswith((".", "scripts", "docs", "presentation")):
        continue
    rfiles = list(d.glob("*_mrgsolve_model.R")) + list(d.glob("*mrgsolve*.R"))
    if not rfiles:
        continue
    rtext = rfiles[0].read_text(encoding="utf-8", errors="ignore")
    odes = count_odes(rtext)
    # scenarios: count ev( or scenario labels
    scen = len(re.findall(r"\bev\s*\(", rtext))
    nodes = clusters = 0
    dots = list(d.glob("*.dot"))
    if dots:
        nodes, clusters = map_stats(dots[0].read_text(encoding="utf-8", errors="ignore"))
    refs = 0
    reftxt = ""
    reff = list(d.glob("*_references.md")) + list(d.glob("*reference*.md"))
    if reff:
        reftxt = reff[0].read_text(encoding="utf-8", errors="ignore")
        refs = len(re.findall(r"PMID[:\s]*\d+", reftxt, re.I)) or len(re.findall(r"https?://pubmed", reftxt))
    hay = rtext + "\n" + reftxt
    safety_hits = len(safe_re.findall(hay))
    models.append({
        "dir": d.name,
        "category": cat_by_dir.get(d.name, "기타"),
        "odes": odes, "scenarios": scen,
        "nodes": nodes, "clusters": clusters,
        "refs": refs, "safety_hits": safety_hits,
        "safety": safety_hits >= 3,
    })

print("models mined:", len(models))

# ---- growth curve from git (first appearance of each model dir) ------------
growth = []
try:
    # list commits oldest->newest with date, then diff names; cheaper: use log --name-only --diff-filter=A on README not reliable.
    out = subprocess.run(
        ["git", "-C", str(ROOT), "log", "--reverse", "--date=short",
         "--pretty=format:%cd", "--name-only", "--diff-filter=A"],
        capture_output=True, text=True, timeout=120).stdout
    seen = set()
    cur_date = None
    by_date = defaultdict(int)
    for line in out.splitlines():
        if re.match(r"^\d{4}-\d{2}-\d{2}$", line):
            cur_date = line
        elif "/" in line and cur_date:
            slug = line.split("/", 1)[0]
            if slug in cat_by_dir and slug not in seen:
                seen.add(slug)
                by_date[cur_date] += 1
    cum = 0
    for dt in sorted(by_date):
        cum += by_date[dt]
        growth.append({"date": dt, "new": by_date[dt], "cum": cum})
    print("growth points:", len(growth), "final cum", cum)
except Exception as e:
    print("growth failed", e)

summary = {
    "n_models": len(models),
    "n_safety": sum(m["safety"] for m in models),
    "total_refs": sum(m["refs"] for m in models),
    "total_odes": sum(m["odes"] for m in models),
    "total_clusters": sum(m["clusters"] for m in models),
    "median_odes": sorted(m["odes"] for m in models)[len(models)//2],
    "median_refs": sorted(m["refs"] for m in models)[len(models)//2],
    "median_nodes": sorted(m["nodes"] for m in models)[len(models)//2],
    "cat_counts": cat_counts,
    "safety_by_cat": dict(Counter(m["category"] for m in models if m["safety"])),
    "total_by_cat": dict(Counter(m["category"] for m in models)),
}
(OUT / "library_models.json").write_text(json.dumps(models, ensure_ascii=False, indent=1))
(OUT / "library_growth.json").write_text(json.dumps(growth, ensure_ascii=False, indent=1))
(OUT / "library_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1))
print(json.dumps(summary, ensure_ascii=False, indent=1))
