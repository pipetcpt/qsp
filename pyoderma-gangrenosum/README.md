# Pyoderma Gangrenosum (PG) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Pyoderma%20Gangrenosum-purple)]()
[![Category](https://img.shields.io/badge/Category-Neutrophilic%20Dermatosis-orange)]()
[![Nodes](https://img.shields.io/badge/Mechanistic%20Map-200%2B%20nodes-blue)]()
[![ODEs](https://img.shields.io/badge/ODEs-30%20compartments-green)]()
[![Drugs](https://img.shields.io/badge/Drugs-6%20modelled-red)]()
[![Scenarios](https://img.shields.io/badge/Scenarios-9-yellow)]()

## Disease Overview

**Pyoderma Gangrenosum (PG, 괴저성 농피증)** is a rare, ulcerative, neutrophilic dermatosis characterized by rapidly progressive painful skin ulcers with violaceous undermined borders, sterile neutrophilic infiltrates, and a strong pathergy response (lesions provoked or worsened by minor trauma). It is classified within the **autoinflammatory spectrum** and frequently coexists with **inflammatory bowel disease (IBD)**, **rheumatoid arthritis**, **hematologic disorders (MDS/AML, monoclonal gammopathy)**, and the syndromic forms **PAPA**, **PASH**, and **PAPASH** driven by PSTPIP1 mutations. Incidence ≈ 3–10 per million per year; mortality 16% at 5 years (driven by comorbidities and immunosuppression).

---

## Key Pathophysiological Mechanisms

| Cluster | Key Components |
|---------|---------------|
| **Genetics / trigger** | PSTPIP1 / MEFV / NLRP3 GoF · pathergy (trauma) · drug triggers (cocaine, PTU, G-CSF, TKIs) |
| **Inflammasome** | NLRP3 · Pyrin · AIM2 → Caspase-1 → IL-1β / IL-18 · gasdermin-D → pyroptosis |
| **Neutrophil dysreg.** | IL-8 / CXCL1-2 / LTB4 / C5a recruitment · NOX2 → ROS · NETosis (PAD4 → CitH3) · MPO / elastase / Cathepsin G |
| **Th17 / IL-23 axis** | DC → IL-23 → JAK1-2 / TYK2 → STAT3 → RORγt → IL-17A/F / IL-22 / GM-CSF |
| **TNF / IL-6 amplification** | TNF-α → NF-κB → IL-6 → STAT3 → CRP/SAA/Ferritin |
| **Tissue destruction** | MMP-2/3/9/10 + Elastase → dermal matrix digestion → sterile ulcer · keratinocyte apoptosis · violaceous edge |
| **Healing failure** | TGF-β / VEGF impaired by TNF-α / IL-17 → cribriform scar |
| **Comorbid bridges** | IBD gut translocation · synovitis (RA) · clonal hematopoiesis (MDS/AML) · monoclonal IgA (MGUS) |

---

## Drug PK/PD Parameters (modeled)

| Drug | Mechanism | F | t½ | EC₅₀ | Trial / Evidence |
|------|-----------|---|-----|------|------------------|
| **Prednisone** | GR → ↓NF-κB / AP-1 | 1.0 | 3–4 h | ~5 ng/mL | STOP-GAP (Ormerod 2015) |
| **Cyclosporine** | Calcineurin / NFAT → ↓Th17, ↓IL-2 | 0.30 | 8–27 h | 100 ng/mL | STOP-GAP (Ormerod 2015) |
| **Infliximab** | Anti-TNF-α | – | 9 d | 0.10 μg/mL | Brooklyn 2006 (only RCT) |
| **Adalimumab** | Anti-TNF-α | 0.64 | 14 d | 0.20 μg/mL | Yamasaki 2022 (Japanese RCT) |
| **Anakinra** | IL-1Rα | 0.95 | 4–6 h | 0.50 μg/mL | Brenner 2009 (PAPA) |
| **Ustekinumab** | Anti-IL-12/23 (p40) | 0.57 | 21 d | 0.50 μg/mL | Goldminz 2012 |

---

## Clinical Endpoint Definitions

| Endpoint | Definition |
|----------|-----------|
| **PARACELSUS** | 0–60 weighted score (Jockenhöfer 2019) — combines progression, exclusion of differentials, violaceous border, response to immunosuppression, pain, size, depth, pathergy |
| **Ulcer area** | Direct planimetry (cm²); primary endpoint in most series |
| **Time to healing** | Days/weeks to complete re-epithelialization (no open ulcer) |
| **Pain VAS** | 0–10 visual analogue |
| **DLQI** | Dermatology Life Quality Index 0–30 |
| **% BSA** | Body-surface-area ulceration |
| **Relapse rate** | Number of recurrences per year |
| **Complete remission** | No active lesions for ≥ 8 weeks off ≥ 2 immunosuppressants |

---

## Model Deliverables

| File | Description |
|------|-------------|
| `pg_qsp_model.dot` | Graphviz mechanistic map source (200+ nodes, 10 clusters) |
| `pg_qsp_model.svg` | Scalable vector map (high-resolution) |
| `pg_qsp_model.png` | Raster map at 150 dpi |
| `pg_mrgsolve_model.R` | 30-compartment mrgsolve QSP model; 9 treatment scenarios; virtual population helper |
| `pg_shiny_app.R` | 7-tab interactive Shiny dashboard |
| `pg_references.md` | 40+ PubMed-linked references organized by category |

---

## Quick Start

```bash
# Render mechanistic map
dot -Tsvg pg_qsp_model.dot -o pg_qsp_model.svg
dot -Tpng -Gdpi=150 pg_qsp_model.dot -o pg_qsp_model.png
```

```r
# Run mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr"))
source("pg_mrgsolve_model.R")

# Single scenario
df <- simulate_pg(pg_qsp_model, "Infliximab", tmax = 168)
plot(df$time, df$Ulcer, type = "l")

# Launch Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT","shinycssloaders"))
shiny::runApp("pg_shiny_app.R")
```

---

## Shiny Dashboard Tabs

| Tab | Contents |
|-----|----------|
| **1. Patient Profile** | Age, sex, weight, PG variant, location, prior therapy, severity index |
| **2. Drug PK** | Concentration–time profiles for ADA, IFX, ANA, CsA, UST, PRED |
| **3. PD — Cytokines** | TNF-α, IL-1β, IL-17A, IL-6, IL-23, IL-8; % suppression vs. untreated |
| **4. Cellular Inflammation** | Neutrophils, NET burden, Th17/Treg, M1 macrophage, MMP-9, ROS, CRP, Calprotectin |
| **5. Clinical Endpoints** | Ulcer area (cm²), % healed, PARACELSUS, Pain VAS, DLQI, complete healing |
| **6. Scenario Comparison** | 9-arm comparative ulcer trajectory + endpoint table @ Week 24 |
| **7. Virtual Population** | VPOP (n configurable) — distribution of healed fraction, PARACELSUS spaghetti |

---

## Modeled Treatment Scenarios

1. **NoTreatment** — natural history reference
2. **Prednisone_SOC** — 60 mg/d × 2 wk → taper × 12 wk
3. **Cyclosporine** — 4 mg/kg/d (BID), 6 months
4. **Infliximab** — 5 mg/kg @ 0, 2, 6 wk → q8w
5. **Adalimumab** — 80 mg loading → 40 mg q1w
6. **Anakinra** — 100 mg/d SC (PAPA-spectrum)
7. **Ustekinumab** — 90 mg q12w
8. **Combo_PRED_CSA** — STOP-GAP-style induction + maintenance
9. **Combo_IFX_low_CsA** — biologic + low-dose calcineurin inhibitor

---

## References (Key)

- Brooklyn TN et al. *Gut.* 2006;55:505 (Infliximab RCT — only RCT in PG)
- Ormerod AD et al. *BMJ.* 2015;350:h2958 (STOP-GAP — Prednisone vs. CsA)
- Maverakis E et al. *JAMA Dermatol.* 2018;154:461 (Delphi diagnostic criteria)
- Jockenhöfer F et al. *Br J Dermatol.* 2019;180:615 (PARACELSUS score)
- Brenner M et al. *Br J Dermatol.* 2009;161:1199 (Anakinra in PAPA)
- Goldminz AM et al. *J Am Acad Dermatol.* 2012;67:e237 (Ustekinumab)
- Yamasaki K et al. *J Dermatol.* 2022;49:479 (Adalimumab phase 3 — Japan)
- Full bibliography: [`pg_references.md`](pg_references.md)
