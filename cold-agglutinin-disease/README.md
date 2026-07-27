# Cold Agglutinin Disease (CAD) — QSP Model

> Integrated Quantitative Systems Pharmacology model of Cold Agglutinin
> Disease, linking clonal marrow B-cell/lymphoplasmacytic monoclonal IgM
> cold-agglutinin production to cold-dependent RBC binding, classical
> complement pathway activation (C1q-C1r-C1s → C4 → C3), C3b-mediated
> **extravascular** hepatic Kupffer-cell hemolysis, C5b-9 (MAC)-mediated
> **intravascular** hemolysis, and downstream hemolysis biomarkers (Hb,
> reticulocytes, bilirubin, LDH, haptoglobin) — with the modern pharmacology
> stack (rituximab, bendamustine, sutimlimab, pegcetacoplan) and the clinical
> endpoints used to track and manage this **steroid-refractory** disease.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`cad_qsp_model.dot`](cad_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`cad_qsp_model.svg`](cad_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`cad_qsp_model.png`](cad_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`cad_mrgsolve_model.R`](cad_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`cad_shiny_app.R`](cad_shiny_app.R) |
| 📚 References             | [`cad_references.md`](cad_references.md) |

---

## 1. Disease in one paragraph

Cold Agglutinin Disease (CAD) is a rare, clonal B-cell/lymphoplasmacytic
lymphoproliferative disorder of the bone marrow (IGHV4-34-restricted,
typically **MYD88 wild-type** — distinguishing it from Waldenström
macroglobulinemia) in which the malignant clone secretes a **monoclonal IgM**
with cold-dependent anti-I/i specificity. At peripheral/acral temperatures
below core body temperature, this IgM binds RBC surface antigens and
**agglutinates** red cells, and its pentameric structure is a potent activator
of the **classical complement pathway**: C1q binding drives C1r/C1s serine
protease activation, C4 and C3 cleavage, and **C3b deposition** on the RBC
surface. Opsonized cells are predominantly cleared by **hepatic Kupffer-cell
CR1/CRIg-mediated erythrophagocytosis (extravascular hemolysis)**, while a
variable fraction proceeds to full **C5b-9 membrane attack complex (MAC)
formation and intravascular hemolysis**. The result is chronic hemolytic
anemia with reticulocytosis, indirect hyperbilirubinemia, elevated LDH, and
haptoglobin depletion, punctuated by **acute cold-exposure hemolytic crises**
and accompanied by acrocyanosis/Raynaud-like phenomena. Critically, CAD is
**characteristically steroid-refractory** (unlike warm AIHA), so therapy
targets either the clonal B-cell compartment (**rituximab**, alone or with
**bendamustine**) or the complement cascade directly (**sutimlimab**, an
anti-C1s antibody approved for CAD; **pegcetacoplan**, an investigational C3
inhibitor with dual classical+alternative pathway blockade).

## 2. Mechanistic clusters (14 in the DOT map, 116 nodes)

1. Genetic / clonal predisposition (MYD88 wild-type, KMT2D mutation, CXCR4 variant, trisomy 3, age-related immunosenescence, Mycoplasma/EBV infection triggers for secondary CAD)
2. Clonal B-cell lymphoproliferation & IgM production (bone marrow clonal B-cell population, IL-6/BAFF microenvironment, plasmacytoid differentiation, IgM pentamer assembly, anti-I/anti-i specificity, free light chain)
3. Cold-dependent antibody-RBC binding (peripheral/acral cooling, thermal amplitude, IgM conformational change, RBC I/i antigen, agglutination, microvascular sludging, rewarming dissociation)
4. Classical complement pathway activation (C1q-C1r-C1s, C4/C2 cleavage, C3 convertase C4b2a, C3b/iC3b/C3d deposition, RBC complement regulators CD55/CD59/MCP)
5. Terminal complement / MAC (C5 convertase, C5a anaphylatoxin, C5b-6-7 complex, C5b-9 MAC, RBC membrane pore, intravascular hemolysis, free plasma hemoglobin)
6. Extravascular hemolysis — liver Kupffer cells (CR1/CRIg receptor, phagocytosis, erythrophagocytosis, partial spherocyte formation, hepatic iron recycling)
7. Hemolysis consequences / biomarkers (hemoglobin, anemia severity, reticulocytes, indirect bilirubin, LDH, haptoglobin depletion, hemoglobinuria, compensatory erythropoiesis/EPO)
8. Clinical manifestations (fatigue, acrocyanosis, Raynaud-like phenomenon, livedo reticularis, cold-induced hemolytic crisis, jaundice, exertional dyspnea, cognitive fatigue, QoL/FACIT-fatigue)
9. Associated/secondary conditions (lymphoplasmacytic lymphoma transformation, Waldenström overlap, secondary cold agglutinin syndrome, complement-mediated thrombosis risk)
10. Drug PK/PD — rituximab (anti-CD20 B-cell depletion, ADCC/CDC lysis, IgM titer decline, relapse via B-cell repopulation)
11. Drug PK/PD — bendamustine (DNA cross-linking, clonal apoptosis, marrow suppression, synergy with rituximab)
12. Drug PK/PD — sutimlimab (anti-C1s target binding/TMDD, classical-pathway blockade, rapid Hb response, infusion-reaction risk)
13. Drug PK/PD — pegcetacoplan (C3/C3b binding, C3-convertase blockade, dual classical+alternative pathway inhibition, breakthrough-hemolysis risk)
14. Clinical endpoints & management (transfusion requirement, Hb-response criteria, steroid-refractory status, cold avoidance, folate supplementation, thromboprophylaxis)

## 3. mrgsolve model (21 ODE compartments)

* **Drug PK (8 compartments)** — rituximab (central + peripheral), bendamustine
  (central + peripheral), sutimlimab (central + peripheral), pegcetacoplan (SC
  depot + central). As with other surrogates in this library, each drug
  compartment is a dose-proportional "exposure" state mapped through an Emax
  function onto its mechanistic target, not a literal multi-compartment
  plasma PK/TMDD model (though sutimlimab's saturable elimination is
  approximated).
* **Clonal / complement axis (7 compartments)** — clonal B-cell burden,
  monoclonal IgM titer, an acute cold-exposure pulse driver, active C1s
  enzymatic state, C4, C3, C3b RBC-surface opsonization, and C5b-9 MAC (8
  compartments total including the cold-pulse driver).
* **Hemolysis biomarkers (5 compartments)** — hemoglobin, reticulocyte
  percentage, indirect bilirubin, LDH, haptoglobin.

### 7 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history (cold avoidance only) | Berentsen 2020 (232-patient natural history) |
| 2 | Rituximab monotherapy (weekly x4) | Berentsen 2004 Blood 103:2925 (~54% ORR, relapse common) |
| 3 | Rituximab + Bendamustine (RB, q28d x4 cycles) | Berentsen 2017 Blood 130:537 (71% ORR, 40% CR) |
| 4 | Sutimlimab load (d0/d7) + biweekly maintenance | CARDINAL (Röth 2021 NEJM 384:1323, ~54% response, rapid onset) |
| 5 | Pegcetacoplan SC biweekly (investigational) | Gertz 2025 Blood 145:397 (dual-pathway C3 blockade) |
| 6 | Sutimlimab discontinued at day 90 (relapse) | Post-discontinuation hemolysis recurrence (CARDINAL/CADENZA extensions) |
| 7 | Acute cold-exposure hemolytic crisis (on sutimlimab) | Jäger 2020 Blood Rev 41:100648 (thermal-amplitude-driven crises) |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — severity slider, secondary-CAD flag, presenting Hb, ambient cold exposure, cold-avoidance counseling flag.
2. **Drug PK** — fractional drug-effect trajectories for all four agents.
3. **Complement pathway PD** — classical-pathway flux/C4/C3, C3b opsonization & MAC, clonal burden & IgM titer.
4. **Hemolysis biomarkers** — hemoglobin & reticulocytes, bilirubin & LDH, haptoglobin.
5. **Clinical endpoints** — endpoint summary table (Day-0, Day-30, Day-90, Day-end).
6. **Scenario comparison** — runs all 7 built-in scenarios with the chosen patient profile.
7. **Cold-exposure challenge** — dedicated view of the acute cold-pulse scenario's transient MAC/hemolysis spike and Hb dip.
8. **References** — links to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg cad_qsp_model.dot -o cad_qsp_model.svg
dot -Tpng -Gdpi=150 cad_qsp_model.dot -o cad_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("cad_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=1.3, cmt="SUT_CENT", ii=7, addl=1) +
                    ev(amt=1.1, cmt="SUT_CENT", time=14, ii=14, addl=11), end=180)
plot(out, c("Hemoglobin","MAC_level","Clonal_burden","LDH_level"))

# 3) Launch the dashboard
shiny::runApp("cad_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Composite response | Sutimlimab (CARDINAL) | ~54% at 26 weeks, Hb rise within days |
| Composite response | Sutimlimab (CADENZA, no recent transfusion) | ~73% vs ~15% placebo |
| Overall response rate | Rituximab-bendamustine vs rituximab alone | 71% (40% CR) vs ~45-54% (CR rare) |
| Median response duration | Rituximab monotherapy | ~11-15 months (relapse common) |
| Hb change at week 48 | Pegcetacoplan (CAD arm) | Median +2.4 g/dL |
| Corticosteroid response | CAD vs warm AIHA | CAD is characteristically steroid-refractory |
| Diagnostic threshold | Cold agglutinin titer | ≥64 at 4°C, strongly positive DAT for C3d |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* Drug "exposure" compartments are pragmatic dose-proportional surrogates for
  sustained mechanistic exposure, not literal multi-compartment plasma
  PK/TMDD models.
* Complement compartments (C1s, C4, C3, C3b, MAC) are simplified lumped
  states representing pathway *flux*, not a full biophysical cascade model
  with individual rate constants for every convertase step.
* The extravascular-vs-intravascular hemolysis split (`W_C3B_EV` vs
  `W_MAC_IV`) is a qualitative weighting consistent with CAD's known
  predominance of extravascular (Kupffer-cell) clearance, not a validated
  quantitative partition.
* Corticosteroids are intentionally **not** modeled as a treatment arm,
  reflecting CAD's well-documented steroid-refractory nature — this is a
  deliberate omission, not an oversight.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
