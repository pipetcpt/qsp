# Idiopathic Multicentric Castleman Disease (iMCD) — QSP Model

> HHV-8-negative lymphoproliferative cytokine storm driven by IL-6 dysregulation,
> with overlapping mTORC1 hyperactivation in the TAFRO subtype.
> Models siltuximab (FDA-approved), tocilizumab, sirolimus, rituximab, anakinra,
> ruxolitinib, CHOP-like cytotoxics, and prednisone.

[![map](imcd_qsp_model.png)](imcd_qsp_model.svg)

## 📦 Deliverables in this directory

| File | Purpose |
|------|---------|
| [`imcd_qsp_model.dot`](imcd_qsp_model.dot) | Graphviz source for the 14-cluster / 150+ node mechanistic map |
| [`imcd_qsp_model.svg`](imcd_qsp_model.svg) / [`.png`](imcd_qsp_model.png) | Rendered mechanistic map |
| [`imcd_mrgsolve_model.R`](imcd_mrgsolve_model.R) | mrgsolve ODE model (19 compartments, 7 scenarios) |
| [`imcd_shiny_app.R`](imcd_shiny_app.R) | Interactive 8-tab Shiny dashboard |
| [`imcd_references.md`](imcd_references.md) | 60+ curated PubMed references |

## 🧬 Disease Summary

iMCD is the HHV-8-negative form of multicentric Castleman disease — a rare but
life-threatening polyclonal lymphoproliferative disorder. Patients present with
multistation lymphadenopathy and characteristic lymph-node histology
(hyperplastic / plasmacytic / mixed / hypervascular), together with a systemic
cytokine storm dominated by **IL-6**. Manifestations include B-symptoms,
hepatosplenomegaly, anasarca/effusions, anemia, polyclonal hypergammaglobulinemia,
elevated CRP/ESR/ferritin, and (in TAFRO subtype) thrombocytopenia and reticulin
fibrosis. Diagnosis follows the CDCN 2017 criteria (Fajgenbaum 2017 *Blood*).

The **TAFRO** subtype (Thrombocytopenia, Anasarca, Fever, Reticulin fibrosis,
Organomegaly) often is refractory to IL-6 blockade alone and shows hyperactive
PI3K/AKT/mTORC1 signaling in Tfh cells (Fajgenbaum 2019 *JCI*), motivating
sirolimus.

## 🎯 Core Pathway Map (14 clusters)

```
①  Genetic & triggers       ②  LN histopathology         ③  Cellular players
④  IL-6 / sIL-6R / gp130    ⑤  IL-1, TNF, VEGF, IL-18    ⑥  PI3K-AKT-mTORC1
⑦  NLRP3 inflammasome       ⑧  Hepatic acute phase       ⑨  Anemia / cytopenias
⑩  Multi-organ effects       ⑪  Drug PK compartments      ⑫  Drug MoA
⑬  CDCN diagnostic & resp.   ⑭  Safety / AE
```

## ⚙️ ODE Model Summary

19 compartments :
- Drug PK (siltuximab, tocilizumab, sirolimus, rituximab, anakinra,
  ruxolitinib, doxorubicin, cyclophosphamide, prednisone)
- IL-6 total / free (siltuximab traps IL-6; free IL-6 drops, total rises)
- Lymph node composite size · Plasmablast fraction · Memory B cell
- CRP · Hb (hepcidin–mediated anemia) · IgG (polyclonal)
- VEGF · Anasarca · Platelet (TAFRO) · mTORC1 activity
- Cumulative OS hazard

7 pre-defined scenarios :
1. Untreated natural history
2. **Siltuximab 11 mg/kg q3w** (CONCERT label)
3. **Tocilizumab 8 mg/kg q2w** (Nishimoto regimen)
4. **Sirolimus 2 mg QD** (TAFRO Fajgenbaum 2019)
5. **Rituximab × 4 + Prednisone**
6. **CHOP + Siltuximab**
7. **Siltuximab + Sirolimus + Anakinra** (refractory TAFRO triple)

Calibration anchors (CONCERT, Nishimoto 2005, Fajgenbaum 2019,
van Rhee 2018 consensus, Pierson 2023 real-world) are documented inline.

## 📊 Shiny Dashboard

Eight tabs let users:
1. Configure patient profile (subtype, weight, baseline labs)
2. View drug PK profiles
3. Inspect IL-6 axis and mTOR dynamics
4. Track acute-phase response (CRP, Hb, IgG)
5. Observe TAFRO/VEGF/anasarca/platelet
6. Watch lymph node and plasmablast burden
7. Compare pre-defined scenarios head-to-head
8. Read the CDCN composite response and survival curve

## ▶️ Running

```bash
# render map
dot -Tsvg idiopathic-multicentric-castleman-disease/imcd_qsp_model.dot \
    -o imcd_qsp_model.svg
dot -Tpng -Gdpi=150 idiopathic-multicentric-castleman-disease/imcd_qsp_model.dot \
    -o imcd_qsp_model.png
```

```r
# install dependencies (once)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny", "DT"))

# source the ODE model
source("idiopathic-multicentric-castleman-disease/imcd_mrgsolve_model.R")

# launch the dashboard
shiny::runApp("idiopathic-multicentric-castleman-disease/imcd_shiny_app.R")
```

## ⚠️ Disclaimer

This is an educational, research-only QSP model. It is not a clinical decision
tool. All parameters are illustrative and anchored to published literature; they
have **not** been formally validated against prospective patient-level data.
Clinical decisions about iMCD must be guided by the CDCN consensus treatment
guidelines (van Rhee 2018 *Blood*) and a qualified hematology team.

## 📚 Key references

- Fajgenbaum DC et al. CDCN diagnostic criteria. *Blood* 2017
- van Rhee F et al. CONCERT trial (siltuximab). *Lancet Oncol* 2014
- van Rhee F et al. Consensus treatment guidelines. *Blood* 2018
- Fajgenbaum DC et al. PI3K/mTOR in IL-6-refractory iMCD. *JCI* 2019
- Carbone A et al. Castleman disease review. *Nat Rev Dis Primers* 2021

See [imcd_references.md](imcd_references.md) for the full curated list.
