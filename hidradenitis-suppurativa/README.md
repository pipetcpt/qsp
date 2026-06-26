# Hidradenitis Suppurativa (HS) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Hidradenitis%20Suppurativa-purple)]()
[![Category](https://img.shields.io/badge/Category-Autoimmune%20%2F%20Chronic%20Skin-orange)]()
[![Nodes](https://img.shields.io/badge/Mechanistic%20Map-160%2B%20nodes-blue)]()
[![ODEs](https://img.shields.io/badge/ODEs-20%20compartments-green)]()
[![Drugs](https://img.shields.io/badge/Drugs-5%20modelled-red)]()

## Disease Overview

**Hidradenitis Suppurativa (HS, 화농성 한선염)** is a chronic, relapsing, inflammatory skin disorder affecting apocrine gland-bearing areas (axillae, groin, inframammary folds, perianal region). Characterized by recurrent painful nodules, abscesses, and sinus tracts/fistulae, it severely impairs quality of life. Prevalence ~1–4% globally; predominantly affects women of reproductive age. The pathophysiology integrates follicular occlusion, dysbiotic microbiome, innate/adaptive immune dysregulation, hormonal factors, and progressive fibrosis.

---

## Key Pathophysiological Mechanisms

| Cluster | Key Components |
|---------|---------------|
| **Follicular occlusion** | γ-Secretase complex → ↓Notch signaling → hyperkeratosis → comedone → rupture |
| **Innate immunity** | NLRP3 inflammasome → IL-1β, IL-18; TLR2/4/9 → NF-κB → TNF-α, IL-6, IL-8 |
| **Adaptive immunity** | Th17 (IL-17A/F, IL-22) + Th1 (IFN-γ); ↓Treg; IL-23 amplification loop |
| **TNF-α signaling** | TNFR1→TRAF2→IKKβ→NF-κB feedforward; ADAM17 shedding |
| **Hormonal factors** | DHT↑ → sebum↑, keratinocyte proliferation; obesity → insulin resistance → mTOR |
| **Microbiome/biofilm** | S. aureus, anaerobes → biofilm → AMR; ↓β-defensins, skin barrier dysfunction |
| **Fibrosis** | TGF-β → myofibroblast → collagen → sinus tracts; MMP-3/9 remodelling |

---

## Drug PK/PD Parameters

| Drug | Mechanism | F | t½ | EC₅₀ | Trial Reference |
|------|-----------|---|-----|------|-----------------|
| **Adalimumab** | Anti-TNF-α | 0.64 | 14 d | 0.20 μg/mL | PIONEER I/II |
| **Secukinumab** | Anti-IL-17A | 0.73 | 28 d | 0.08 μg/mL | SUNSHINE/SUNRISE |
| **Bimekizumab** | Anti-IL-17A/F | 0.74 | 23 d | 0.06 μg/mL | BE HEARD I/II |
| **Ixekizumab** | Anti-IL-17A | 0.60 | 13 d | 0.10 μg/mL | IXORA-HS |
| **Risankizumab** | Anti-IL-23 (p19) | 0.62 | 28 d | 0.12 μg/mL | Phase 2 data |

---

## Clinical Endpoint Definitions

| Endpoint | Definition |
|----------|-----------|
| **HiSCR** | ≥50% reduction in AN count (abscesses + nodules) with no increase in abscess or sinus tract count |
| **IHS4** | (Nodules × 1) + (Abscesses × 2) + (Draining fistulae × 4); mild <4, moderate 4–10, severe >10 |
| **AN count** | Total number of inflammatory nodules and abscesses |
| **DLQI** | Dermatology Life Quality Index (0–30, higher=worse) |
| **VAS pain** | Visual analogue scale 0–10 |
| **Hurley stage** | I (isolated abscess), II (recurrent + tracts), III (diffuse/extensive) |

---

## Model Deliverables

| File | Description |
|------|-------------|
| `hs_qsp_model.dot` | Graphviz mechanistic map source (160+ nodes, 10 clusters) |
| `hs_qsp_model.svg` | Scalable vector map (high-resolution) |
| `hs_qsp_model.png` | Raster map at 150 dpi |
| `hs_mrgsolve_model.R` | 20-ODE mrgsolve model; 5 treatment scenarios; virtual population |
| `hs_shiny_app.R` | 6-tab interactive Shiny dashboard |
| `hs_references.md` | 37 PubMed-linked references by category |

---

## Quick Start

```bash
# Render mechanistic map
dot -Tsvg hs_qsp_model.dot -o hs_qsp_model.svg
dot -Tpng -Gdpi=150 hs_qsp_model.dot -o hs_qsp_model.png
```

```r
# Run mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("hs_mrgsolve_model.R")

# Launch Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT","shinycssloaders"))
shiny::runApp("hs_shiny_app.R")
```

---

## Shiny Dashboard Tabs

| Tab | Contents |
|-----|----------|
| **1. Patient Profile** | Age, sex, BMI, Hurley stage, prior therapy, baseline IHS4/DLQI |
| **2. Drug PK** | Concentration–time profiles for ADA, SEC, BIM; trough monitoring |
| **3. PD — Cytokines** | TNF-α, IL-17A, IL-6, IL-1β, IL-23; Th17/M1/Neutrophil indices |
| **4. Clinical Endpoints** | AN count, IHS4, HiSCR, DLQI, VAS pain, fistula score |
| **5. Scenario Comparison** | 5-arm efficacy comparison; comparative table; HiSCR bar chart |
| **6. Biomarker & VPop** | Spider plot, Th17/Treg, virtual population HiSCR distribution, smoking impact |

---

## References (Key)

- Kimball AB et al. *N Engl J Med.* 2016;375:422 (PIONEER I/II — Adalimumab)
- Kimball AB et al. *Lancet.* 2023;401:747 (SUNSHINE/SUNRISE — Secukinumab)
- Mughal AA et al. *Lancet.* 2023;402:1881 (BE HEARD I/II — Bimekizumab)
- Zouboulis CC et al. *Br J Dermatol.* 2017;177:1401 (IHS4 validation)
- Full bibliography: [`hs_references.md`](hs_references.md)
