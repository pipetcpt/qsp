# Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN) — QSP Model

> Drug-induced severe cutaneous adverse reaction (SCAR) modeled as a coupled
> pharmacokinetic / immune / epidermal / clinical-outcome system.
> The model links culprit drug exposure (carbamazepine, allopurinol,
> sulfamethoxazole, lamotrigine, etc.) and HLA-driven antigen presentation
> to cytotoxic CD8⁺ CTL expansion, the granulysin / sFasL / TNF-α effector axis,
> full-thickness keratinocyte apoptosis, BSA detachment, SCORTEN-based
> mortality hazard, and re-epithelialization under treatment.

| Output | File |
|---|---|
| 🗺️ Mechanistic map (.dot) | [`sjsten_qsp_model.dot`](sjsten_qsp_model.dot) |
| 🖼️ Rendered map (SVG / PNG) | [`sjsten_qsp_model.svg`](sjsten_qsp_model.svg) · [`sjsten_qsp_model.png`](sjsten_qsp_model.png) |
| ⚙️ mrgsolve model | [`sjsten_mrgsolve_model.R`](sjsten_mrgsolve_model.R) |
| 📊 Shiny dashboard | [`sjsten_shiny_app.R`](sjsten_shiny_app.R) |
| 📚 References (PubMed) | [`sjsten_references.md`](sjsten_references.md) |

## 1. Disease background

SJS/TEN is a delayed-type, T-cell–mediated severe cutaneous adverse reaction
characterized by widespread epidermal apoptosis and mucosal erosions.
Severity is defined by body-surface-area (BSA) detachment:
- **SJS** (<10% BSA), **SJS/TEN overlap** (10–30%), **TEN** (>30%).
- Mortality: SJS ≈ 5 %; TEN 25–35 %; up to ~50 % in older patients.
- Mucosal involvement (>90%): ocular (symblepharon, blindness),
  oral, urogenital strictures.

### Key pharmacogenomic risk alleles

| Allele | Drug | Population | Odds ratio |
|---|---|---|---|
| HLA-B*15:02 | Carbamazepine, phenytoin | Han Chinese, SE Asian | ~113 |
| HLA-B*58:01 | Allopurinol | Han Chinese, Thai, Korean | ~580 |
| HLA-A*31:01 | Carbamazepine | European, Japanese | 8–25 |
| HLA-B*57:01 | Abacavir (AHS) | Caucasian | ~960 |
| HLA-B*13:01 | Dapsone | Han Chinese | ~20 |

## 2. Mechanistic map (.dot)

10 functional clusters and **130+ nodes**:

| # | Cluster | Highlights |
|---|---|---|
| ① | Culprit drugs & PK | CBZ, OXC, LTG, PHT, allopurinol→oxypurinol, SMX, dapsone, nevirapine, oxicam NSAIDs, abacavir; reactive metabolites (NHOH-SMX, CBZ-EPX) |
| ② | Pharmacogenomics | HLA-B*15:02 / B*58:01 / A*31:01 / B*57:01 / B*13:01, CYP2C9*3, public TCR Vβ clonotypes |
| ③ | Drug-HLA-TCR | Hapten · p-i concept · altered peptide repertoire · Langerhans cell cross-presentation |
| ④ | Cytotoxic effectors | CD8⁺ CTL, NK/NKT, granulysin, perforin/GzmB, sFasL, TNF-α, IFN-γ, IL-15, HMGB1 |
| ⑤ | Keratinocyte death | Fas/Caspase-8/3 apoptosis, MLKL/RIPK3 necroptosis, granulysin pore lysis |
| ⑥ | Clinical phenotype | BSA detachment, SCORTEN (age, BSA, HR, BUN, glucose, HCO₃, malignancy) |
| ⑦ | Systemic complications | Sepsis, ARDS, AKI, GI bleed, hypothermia, hypovolemic shock |
| ⑧ | Pharmacologic Rx | Drug withdrawal, IVIG, cyclosporine, etanercept, infliximab, methylprednisolone, TPE, JAK inhibitors, supportive care |
| ⑨ | Biomarkers | Serum granulysin, sFasL, IL-15, HMGB1, LTT, IFN-γ ELISpot |
| ⑩ | Re-epithelialization | Hair-follicle bulge stem cells, EGF/KGF, scarring, dry eye, PTSD, re-exposure |

## 3. mrgsolve ODE model

**23 compartments**, including:
- Culprit drug PK (oral depot + central)
- Antigen-HLA complex
- Activated CTL clones
- TNF-α, IFN-γ, IL-15, granulysin, sFasL, HMGB1 (cytokine/effector pool)
- Keratinocyte alive mass, BSA detachment, cumulative hazard, survival
- PK for IVIG, etanercept, infliximab, cyclosporine, methylprednisolone, JAK inhibitor

**7 treatment scenarios**:
1. Supportive only
2. IVIG 1 g/kg/d × 4
3. Cyclosporine 3 mg/kg/d (Valeyrie-Allanore 2010)
4. Etanercept 25 mg SC × 2 (Wang 2018 JCI RCT)
5. Methylprednisolone 1 mg/kg pulse
6. Cyclosporine + etanercept combination
7. JAK inhibitor (tofacitinib 5 mg BID) — investigational

The mortality hazard is modulated by SCORTEN (Bastuji-Garin 2000) and a
treatment effect multiplier benchmarked to Wang 2018 (etanercept RR ≈ 0.18)
and González-Herrada 2017 (cyclosporine SMR ≈ 0.43).

## 4. Shiny dashboard tabs

1. **Patient profile** — age, weight, HLA risk, culprit drug, dose
2. **Drug PK** — culprit drug + biologic / immunomodulator PK
3. **Immune drivers (PD)** — activated CTL, TNF-α, IFN-γ, IL-15
4. **Clinical endpoints** — BSA detachment, SCORTEN, mortality probability
5. **Scenario comparison** — table of max-BSA, max-SCORTEN, day-14 mortality, time to 50% re-epi
6. **Biomarkers** — granulysin, sFasL, IL-15, HMGB1
7. **References** — embedded literature list

## 5. Usage

```bash
# Render mechanistic map
dot -Tsvg sjsten_qsp_model.dot -o sjsten_qsp_model.svg
dot -Tpng -Gdpi=150 sjsten_qsp_model.dot -o sjsten_qsp_model.png
```

```r
# Run mrgsolve simulation
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread("sjsten_mrgsolve_model.R")
out <- mod %>% param(SCEN_CSA=1, SCEN_ETAN=1) %>%
  ev(amt=400, ii=12, addl=2, cmt="A_drug_dep") %>%
  ev(amt=25, ii=3, addl=1, cmt="ETAN_dep") %>%
  ev(amt=100, ii=0.5, addl=20, cmt="CSA_dep") %>%
  mrgsim(end=30, delta=0.1)
plot(out, BSA_loss + PredMort ~ time)

# Run interactive dashboard
shiny::runApp("sjsten_shiny_app.R")
```

## 6. Caveats

- Parameters are illustrative and intended for **teaching / hypothesis
  generation**, not bedside decision support.
- The SCORTEN hazard mapping uses a simplified exponential link; real-world
  SCORTEN bins map to specific predicted mortalities.
- Treatment effect multipliers reflect observational SMRs and a single RCT;
  RegiSCAR meta-analyses are mixed for IVIG.
- HLA effects are encoded as a binary risk multiplier; population-specific
  carrier frequencies and screening cut-points are not modeled.
