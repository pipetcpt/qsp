import json, time
from refs import esearch, esummary
R3=[("9","PRIME trial sacubitril valsartan mitral regurgitation","PRIME: ARNI reduced functional MR beyond valsartan alone -- the model's combined afterload-plus-remodelling prediction."),
    ("9","Brater diuretic therapy New England Journal","Loop diuretic pharmacology, including the braking phenomenon the model implements."),
    ("2","proximal flow convergence hemispheric assumption regurgitant orifice shape","Where the hemispheric assumption behind PISA fails."),
    ("3","Yellin dynamic determinants acute mitral regurgitation left atrial","Acute MR haemodynamics and the role of atrial properties."),
    ("2","vena contracta area three-dimensional functional mitral regurgitation crescent","The crescentic functional orifice imaged directly."),
    ("7","reverse remodeling after mitral valve surgery left ventricular","Reverse remodelling after correction, and its incompleteness."),
]
rows=[]
for sec,q,note in R3:
    pid=esearch(q); rows.append((sec,q,note,pid)); time.sleep(0.34)
got=[r[3] for r in rows if r[3]]; meta=esummary(sorted(set(got)))
json.dump({"rows":rows,"meta":meta},open("refs_raw3.json","w"))
for sec,q,note,pid in rows:
    if not pid: print("  UNRESOLVED:",q); continue
    m=meta.get(pid,{}); au=m.get("authors",[])
    print(f"  [{sec:>2}] {pid:>9} {(au[0]['name'] if au else '?'):<16} {m.get('source','?'):<20} {(m.get('pubdate','?') or '?')[:4]} {m.get('title','?')[:86]}")
