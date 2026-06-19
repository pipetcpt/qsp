# NAFLD/NASH QSP Model — Stability & Re-calibration Design Brief

This note documents why the disease model in `nafld_mrgsolve_model.R` and
`nafld_shiny_app.R` is written the way it is, and the literature it is based on.
It was produced after an earlier version **diverged** (placebo grew ~10×/week,
liver fat reaching ~1e+85%) because the disease pools were not initialized at
their own steady state and a positive-feedback loop was unbounded.

## Root cause (original bug)

1. **Production terms were absolute, not baseline-balanced.** e.g.
   `dxdt_KUPFFER = (KOUT_KUP*KUP0 + KLIP_KUP*LIPOTOX/LIPOTOX0) − KOUT_KUP*KUPFFER`
   has algebraic steady state 4.5, but the state is initialized at 0.5 → it leaves
   baseline immediately. Same pattern for TNFA, IL6C, TGFB, HSC, COLLAGEN, ALT.
2. **Liver-fat influx 5× too high:** `KLIN_LF=0.003` vs `KOUT_LF·LF0 = 0.0006`.
3. **Unbounded positive feedback:** liver fat → lipotoxicity → Kupffer → TNF-α →
   insulin resistance → de novo lipogenesis (DNL) → liver fat, with open-loop
   gain ≥ 1, plus an ALT term that multiplied two unbounded drivers.

## Design principles (literature-backed)

- **P1 — Baseline IS the steady state.** Every disease pool is a turnover /
  indirect-response (IDR) pool `dxdt_X = kin_X − kout_X·X` with **`kin_X = kout_X·X0`**.
  Combined with `X(0)=X0`, the initial derivative is exactly 0 → placebo is flat by
  construction. *(Dayneka, Garg & Jusko 1993; Jusko & Ko 1994; Woo, Pawaskar & Jusko 2009.)*
- **P2 — Fold-change normalized drivers.** Cross-talk enters as `(driver/driver0 − 1)`,
  which is 0 at baseline, so the whole loop (not just isolated nodes) is self-consistent.
  *(Goentoro et al. 2009, fold-change detection.)*
- **P3 — Saturable / bounded activators.** The Kupffer drive uses a saturating Emax
  form `1 + g·(d−1)/(1 + γ·(d−1))` (finite ceiling) so a spike cannot run away.
  *(Krzyzanski & Jusko 2006; Goldbeter & Koshland 1981.)*
- **P4 — Loop gain < 1.** A pure positive loop is locally stable iff the product of
  per-edge sensitivities `G < 1`. Here
  `G ≈ GLIP_KUP·GKUP_TNF·KTNF_IR·KDNL_IR·WDNL + KFFA_IR·KDNL_IR·WDNL ≈ 0.03 ≪ 1`.
  *(Angeli, Ferrell & Sontag 2004; Mager, Wyska & Jusko 2003.)*
- **P5 — Timescale separation.** Cytokines (hours) ≪ fat / HOMA-IR (days–weeks) ≪
  collagen (months → `KOUT_COL` made the slowest pool). *(Decaris 2015; Thorsted 2019; Wang 2024.)*
- **P6 — Drugs are bounded multiplicative factors on kin or kout**, never additive
  sources, so placebo is untouched and direction cannot flip. *(Dayneka 1993.)*

## Per-compartment form (as implemented)

Canonical: `dxdt_X = KOUT_X·X0·(1 + Σ gᵢ·(driverᵢ/driver0ᵢ − 1))/INH − KOUT_X·X`.

| Pool | Driver(s) | Drug effects |
|------|-----------|--------------|
| LIVER_FAT | influx = `KOUT_LF·LF0·(WDNL·DNL_n + WUPT·(BODY_WT/WT0))` | RSM/OCA ↓DNL; RSM/SEM ↑efflux |
| INS_RES | set-point: liver fat, TNF-α, adiponectin | SEM/EMP/RSM ↓ |
| KUPFFER | saturable lipotoxicity (=LIVER_FAT/LF0) | OCA/SEM suppress |
| TNFA | Kupffer; adiponectin-protected | OCA/SEM ↓ |
| IL6C | Kupffer | OCA ↓ |
| TGFB | Kupffer + lipotoxicity | OCA/RSM anti-fibrotic |
| HSC | TGF-β | OCA/RSM/SEM block activation |
| COLLAGEN | HSC (slowest pool) | OCA/RSM inhibit synthesis |
| ALT | **sum** of normalized TNF-α + lipotoxicity | falls as drivers fall |
| ADIPONECTIN, BODY_WT | set-point forms (already SS-consistent) | SEM/EMP |

## Verification (see `_test_recalib.R` runs)

1. Placebo flat to machine precision (rel. drift 0) for all 11 states over 72 wk.
2. Perturbation (LIVER_FAT→0.30) recovers to 0.20; no blow-up.
3. No divergence in any arm; states stay physiological.
4. Drug directions at wk72: resmetirom ↓liver fat ~40% (cf. MAESTRO −35/−39%),
   ↓ALT; semaglutide ↓fat/↓HOMA-IR/↓weight/↑adiponectin; OCA ↓collagen/fibrosis;
   triple ≥ monotherapy.

## Key references

- Dayneka NL, Garg V, Jusko WJ (1993). Four basic models of indirect PD responses. *J Pharmacokinet Biopharm* 21:457–478.
- Jusko WJ, Ko HC (1994). Physiologic indirect response models. *Clin Pharmacol Ther* 56:406–419.
- Krzyzanski W, Jusko WJ (2006). IDR models with physiological limits. *J Pharmacokinet Pharmacodyn* 33:635–655.
- Woo S, Pawaskar D, Jusko WJ (2009). Baseline handling for IDR models. *J Pharmacokinet Pharmacodyn* 36:381–405.
- Angeli D, Ferrell JE, Sontag ED (2004). Positive-feedback bistability. *PNAS* 101:1822–1827.
- Goentoro L, et al. (2009). Fold-change detection. *Mol Cell* 36:894–899.
- Mager DE, Wyska E, Jusko WJ (2003). Mechanism-based PD models. *Drug Metab Dispos* 31:510–518.
- Rieger TR, Allen RJ, Musante CJ (2022). A quantitative systems pharmacology model of liver lipid metabolism for investigation of non-alcoholic fatty liver disease. *Front Pharmacol* 13 (PMC9343875).
- Decaris ML, et al. (2015). Turnover of hepatic collagen in humans. *PLOS ONE*.
- Harrison SA, et al. (2023). Resmetirom MAESTRO-NAFLD-1. *Nat Med* 29:2919–2928.
- Newsome PN, et al. (2021). Semaglutide in NASH (Ph2). *N Engl J Med* 384:1113–1124.
- Sanyal AJ, et al. (2023). Obeticholic acid, REGENERATE final analysis. *J Hepatol* 79:1110–1120.
- Sanyal AJ, et al. (2025). Semaglutide in MASH, ESSENCE Ph3. *N Engl J Med* 392:2089–2099.
- Donnelly KL, et al. (2005). Sources of hepatic TAG in NAFLD. *J Clin Invest* 115:1343–1351.
- Kuchay MS, et al. (2018). Empagliflozin on liver fat (E-LIFT). *Diabetes Care* 41:1801–1808.

## v4 recalibration & adversarial-review notes (2026-06-19)

Effect sizes are **calibrated to PLACEBO-CORRECTED (between-group) trial endpoints**, because the QSP
placebo arm is flat by construction (kin = KOUT·X0).

- **Resmetirom**: LF efflux gain 0.6→0.48 and DNL inhibition 0.4→0.30 → 100 mg monotherapy MRI-PDFF
  **−34.0% @wk52**, matching MAESTRO-NAFLD-1 wk52 placebo-adjusted −33.9% (Harrison, Nat Med 2023, Table 3).
  The steady-state model is flat from <wk16, so it maps to the wk52 plateau and does **not** reproduce the
  trial's wk16 −38.6% → wk52 −33.9% within-study decline.
- **Empagliflozin**: EC50_EMP 0.15→0.015 (10 mg was near-inert, E_EMP 0.15→0.60) + new phenomenological
  `WEMP_LF = 0.47` hepatic-fat efflux term → monotherapy liver fat **−24.6%**, matching the E-LIFT
  **placebo-corrected** between-group effect (~−24.7%; empa 16.2→11.3 minus control 16.4→15.5), **not** the
  raw within-arm −30%. Single 20-wk n=50 trial anchor → treat WEMP_LF as uncertain; sensitivity-analyse.

**Documented limitations (from adversarial QSP review):**
- Empagliflozin weight (0.05·E_EMP≈3.3%) and direct HOMA-IR (0.20·E_EMP≈13%) coefficients were *not*
  re-tuned after the EC50 change; the IR term lacks an empagliflozin-paper anchor. Revisit before reporting
  any empagliflozin-containing scenario.
- `WDNL = 0.4` lumps DNL (26%) + dietary (15%) per Donnelly; drug DNL-suppression also acts on the dietary
  fraction (modest over-credit). Future: split WDNL=0.25 + constant WDIET=0.15.
- NAS sub-scores for inflammation (ALT>40) and ballooning (TNFA>1.5×) are crude binary proxies (cap NAS at 5/8);
  use NAS **deltas**, not absolute values.
- OCA activity-pathway coefficients are non-trivial although REGENERATE's NASH-resolution co-primary failed
  (final 6.5% vs 3.5%); realized OCA-monotherapy NAS shows zero net change, so it does not manifest at the
  endpoint. OCA FXR safety liabilities (pruritus, LDL rise) are unrepresented.
- The per-pool IDR/turnover form is a standard Dayneka/Jusko-Ko construct, **not** Rieger's mass-action
  structure; "Rieger" is cited only as QSP precedent + steady-state-initialization rationale.
