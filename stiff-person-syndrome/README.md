# Stiff Person Syndrome (SPS) — QSP Model

Quantitative Systems Pharmacology package for **Stiff Person Syndrome
(SPS)**, an autoimmune CNS disorder driven primarily by IgG
autoantibodies against the 65-kDa isoform of glutamic acid
decarboxylase (anti-GAD65), with intrathecal antibody synthesis, loss
of GABAergic inhibition, and resulting α-motor-neuron
hyper-excitability that produces axial rigidity, painful spasms, and
task-specific phobia.

## Disease in one paragraph

Anti-GAD65 IgG (often > 10,000 U/mL serum, with high CSF/serum index)
inhibits presynaptic GABA synthesis, reduces vesicular GABA loading,
and is permissive for additional surface-antigen autoantibodies
(GABA-A receptor-associated protein, gephyrin, glycine receptor in
the PERM variant, amphiphysin in paraneoplastic SPS). The resulting
GABA-deficit removes Renshaw, Ia-reciprocal and presynaptic
inhibition over spinal α-motor neurons, while cortical and brainstem
GABAergic systems become hyper-excitable. Patients present with
fluctuating stiffness (truncal > limb), stimulus-sensitive spasms,
agoraphobia, and (in severe phenotypes) PERM with brainstem
involvement and dysautonomia. Roughly one-third have type-1 diabetes
and other organ-specific autoimmune comorbidity.

## Files

| File                            | Purpose                                                              |
| ------------------------------- | -------------------------------------------------------------------- |
| `sps_qsp_model.dot`             | Graphviz mechanistic map (124 nodes, 11 clusters)                    |
| `sps_qsp_model.svg`             | Rendered vector mechanistic map                                      |
| `sps_qsp_model.png`             | Rendered raster mechanistic map (150 dpi)                            |
| `sps_mrgsolve_model.R`          | 26-compartment mrgsolve ODE model + 6 scenario presets               |
| `sps_shiny_app.R`               | 8-tab Shiny dashboard (patient, antibody, PK, GABAergic, endpoints, scenario compare, safety, VPC) |
| `sps_references.md`             | 70 curated PubMed references organised in 14 sections                 |
| `README.md`                     | This file                                                            |

## Mechanistic-map clusters

1. Genetic susceptibility & triggers (HLA-DRB1*03:01, PTPN22, AIRE, ICI exposure, paraneoplastic)
2. Adaptive immune cascade (APC, Th1/Th17/Tfh, GC, plasmablasts, LLPC, intrathecal B cells)
3. Autoantibody effectors (anti-GAD65 serum & CSF, anti-amphiphysin, anti-GlyR, anti-DPPX, GABARAP, gephyrin)
4. CNS GABAergic neurochemistry (GAD65, GAD67, VGAT, GAT-1, GABA-A/B, GlyR, NMDA, AMPA)
5. Spinal motor circuitry (Renshaw, Ia-inhibitory IN, presynaptic inhibition, α/γ-MN, H/F-reflex, co-contraction)
6. Supraspinal & autonomic (motor cortex, brainstem RF, cerebellum, amygdala, sympathetic surge)
7. Clinical phenotype (axial/limb stiffness, spasms, gait, falls, PERM, SPS-plus, comorbidities)
8. PK compartments (Diazepam + DMD, Baclofen oral / IT, Gabapentin, Tiagabine, IVIG, Rituximab, PLEX, Prednisolone, MMF, Tacrolimus, BTKi, Bortezomib)
9. Pharmacodynamics & MoA (BZD-PAM, GABA-B agonism, α2δ, GAT-1 block, GABA-T inhibition, CD20/CD19 depletion, FcRn block, IVIG anti-idiotype, GR transrepression, HSCT reset)
10. Biomarkers & endpoints (HSI, stiffness distribution & severity scores, spasm diary, anti-GAD titer, IgG index, EMG, Hmax/Mmax, MRS GABA, falls, hospitalisation)
11. Safety / adverse events (BZD tolerance, baclofen withdrawal, PML, infection, steroid AE, IVIG thrombosis, HSCT TRM)

## ODE model (mrgsolve) at a glance

26 compartments, including:

- 7 PK compartments for diazepam (gut/central/peripheral) + active
  metabolite desmethyldiazepam, plus baclofen (gut/central/CSF) and
  gabapentin (gut/central).
- IVIG, rituximab (2-cmt), and prednisolone PK.
- B-cell / plasmablast / long-lived plasma cell dynamics with
  CD20-mediated rituximab kill and prednisolone-mediated suppression
  of B-cell production.
- Anti-GAD65 antibody compartments (serum & CSF) with IVIG
  anti-idiotype neutralisation.
- GAD65 enzyme activity, central GABA pool, α-MN excitability,
  stiffness score, spasm frequency, BMD trajectory, HSI surrogate.

Drug effects combine multiplicatively on a saturable inhibition term
that suppresses MN excitability:

```
INH_total = E_GABA + E_BZD + E_BAC + E_GAB    (clipped at 0.95)
dMN/dt    = k_in*(1 - INH_total) - k_out*MN
```

Six presets in `sps_run()`:

1. `dx_bzd_only` — newly-diagnosed, diazepam mono
2. `bzd_baclofen_combo` — BZD + oral baclofen
3. `bzd_ivig_q4w` — BZD + IVIG q4w
4. `rtx_induction` — BZD + rituximab 1 g × 2 + steroid pulse
5. `plex_rescue` — refractory crisis with 5-session PLEX
6. `intrathecal_baclofen` — IT pump on chronic baseline BZD

## How to run

```r
# In R:
source("sps_mrgsolve_model.R")
res <- sps_run("rtx_induction", severity = "severe", horizon_d = 180)
head(res)

# Shiny dashboard
shiny::runApp("sps_shiny_app.R")
```

System dependencies (one-time):

```bash
# Graphviz (for re-rendering the mechanistic map)
sudo apt-get install -y graphviz

# R packages
R -e 'install.packages(c("mrgsolve","dplyr","ggplot2","shiny","bslib","plotly","tidyr"))'
```

## Re-render the mechanistic map

```bash
dot -Tsvg sps_qsp_model.dot -o sps_qsp_model.svg
dot -Tpng -Gdpi=150 sps_qsp_model.dot -o sps_qsp_model.png
```

## Calibration anchors (from references)

| Quantity                                  | Value used | Source                                                                                  |
| ----------------------------------------- | ---------- | --------------------------------------------------------------------------------------- |
| Serum anti-GAD65 in active SPS            | 10⁴ U/mL (moderate); 4×10⁴ U/mL (crisis) | Dalakas, *Neurology* 2001                                |
| Serum IgG half-life                       | ≈ 21 d     | Looney, *Curr Allergy Asthma Rep* 2007                                                  |
| Diazepam terminal t½ (parent / DMD)       | ≈ 30 / 50 h | Friedman, *CPT* 1985                                                                    |
| Baclofen oral F / renal CL                | 0.8 / ≈ 11 L/h | Wuis, *Eur J Clin Pharmacol* 1989                                                    |
| Rituximab CL / Vc                         | 0.27 L/d / 3 L | Ng, *J Clin Pharmacol* 2005                                                          |
| IVIG dose                                  | 2 g/kg per cycle | Dalakas RCT, *NEJM* 2001                                                            |
| Brain GABA reduction in SPS               | ~30% lower vs control | Levy, *Arch Neurol* 2005                                                       |
| HLA-DRB1*03:01 OR                         | ≈ 3-5       | Pugliese, *J Clin Endocrinol Metab* 1993                                                |
| α-MN hyperexcitability                    | H/M ratio ↑ ~50% | Floeter, *Neurology* 1998                                                          |

See `sps_references.md` for full citations.

## Caveats

- Parameters are *illustrative* and tuned to qualitative-match
  published SPS observations; they are not from a formal NLME fit.
- The PLEX implementation is a deterministic between-event reset of
  the antibody compartments (~40-60% removal per session); replace
  with stoichiometric pharmacometric model if exchange volumes are
  known.
- This model does not yet include the cerebellar ataxia subtype of
  GAD-antibody disease in mechanistic detail — that requires
  Purkinje-cell-specific synaptic dynamics that are out of scope
  for v1.

## License

Distributed under the repo-root `LICENSE`.
