# Annotation Scaling Decision for AlphaGenome -> SuSiNE

**Date:** 2026-04-09  
**Status:** Decision / execution handoff for HPC implementation  
**Scope:** Real-data case study across 10 loci using AlphaGenome-derived signed annotations with `susine_rss()`

## 1. Problem Statement

We need a coherent way to convert AlphaGenome variant effect predictions into a prior mean vector `mu_0` for SuSiNE. The core mismatch is:

- eQTL summary data come from variation across people and are usually consumed as `z`, `R`, and `n`, or sometimes `slope`, `slope_se`, `R`, and `n`.
- AlphaGenome produces sequence-to-function predictions such as raw predicted log fold change (LFC) and signed quantile scores.
- There is no inherent shared unit system between these sources.

The goal is not to find a perfect "true" conversion from AlphaGenome outputs into eQTL effects. The goal is to build a robust, approximate, and auditable scaling pipeline that places the annotation into the right ballpark for SuSiNE, then hedge residual uncertainty with a grid over the annotation scale.

## 2. Main Decision

The annotation vector should target the scale of the **prior mean effect conditional on inclusion**, `mu_0`, which is in the same unit system as the **posterior conditional effect mean**, `mu_1`.

It should **not** target:

- raw `z`
- unconditional posterior `b_hat`

Reason:

- In the SER update, `mu_1` is the posterior mean effect **conditional on the SNP being the selected SNP for the effect**.
- `b_hat = alpha * mu_1`, so `b_hat` is shrunk by posterior inclusion probability and is not the right calibration target.
- `beta_hat` is the marginal univariate effect estimate used in the Bayes factor and lives in the same effect-size scale family as `mu_0` and `mu_1`.

So the correct conceptual target is:

```text
mu_0  <->  mu_1
```

and the practical pre-fit calibration object is:

```text
beta_hat
```

because `beta_hat` is available directly from the input summary statistics and does not require fitting SuSiE/SuSiNE first.

## 3. Relevant Code Facts from This Repo

These are the code facts that matter for implementation:

- `mu_1` is the posterior conditional mean effect in `SER()` / `SER_ss()`.
- `b_hat = alpha * mu_1`.
- `beta_hat = Xty / d` in full data and `beta_hat = Xty / dXtX` in summary-stat SER.
- `susine_rss()` has two paths:
  - `z`, `R`, `n` only: standardized-effect path
  - `slope`, `slope_se`, `var_y`, `R`, `n`: original-scale path
- Current autoscaling code estimates annotation scale using the unconditional model

```text
beta_hat_j ~ N(c * a_j, shat2_j + tau^2)
```

and returns `c` and `tau^2`.

Files reviewed:

- `susine/R/SER.R`
- `susine/R/SER_ss.R`
- `susine/R/susine.R`
- `susine/R/susine_ss.R`
- `susine/R/susine_rss.R`
- `susine/R/initialize.R`
- `susine/R/update_priors.R`
- `susine/R/update_priors_ss.R`

## 4. High-Level Strategy

We will prepare `mu_0` in two layers:

1. Build a signed, robust, locus-level annotation template `a`
2. Multiply by a global scale `c`

```text
mu_0 = c * a
```

The annotation transform defines the **shape**. The scalar `c` defines the **units**.

We will estimate a pooled `c0` across loci, then use a grid around `c0` rather than trusting a single estimate.

## 5. Input Annotation Choice

### Decision

Use the **signed AlphaGenome quantile score** as the primary source for `a`, not raw LFC.

### Rationale

- raw LFC appears highly concentrated near 0 with a small number of extreme outliers
- raw LFC scale varies strongly across tracks / assays
- signed quantile scores are rank-based and more robust to track-specific scaling and batch effects
- this makes them a better base object for a directional prior

## 6. Recommended Transform from Quantile Score to Annotation Template

Let `q_j in [-1, 1]` be the signed quantile score for SNP `j`.

### Step 1: clip away exact boundaries

```text
q*_j = clip(q_j, -1 + eps, 1 - eps)
```

Recommended default:

```text
eps = 1e-4
```

### Step 2: inverse-normal transform

Use:

```text
a_raw_j = qnorm((q*_j + 1) / 2)
```

This gives the desired behavior:

- `q = -0.5 -> qnorm(0.25)`
- `q = 0 -> 0`
- `q = 0.5 -> qnorm(0.75)`

### Step 3: winsorize

```text
a_clip_j = clip(a_raw_j, -2.5, 2.5)
```

Start with `2.5`. A sensitivity check using `3.0` is reasonable.

### Step 4: RMS-normalize within locus

Do **not** mean-center, because `0` has semantic meaning: "no predicted directional effect."

Instead:

```text
s_a = sqrt(mean(a_clip^2))
a_j = a_clip_j / s_a
```

If `s_a = 0`, set `a_j = 0` for that locus.

This preserves:

- sign
- ordering
- `0 -> 0`

and makes the scale factor `c` interpretable across loci.

## 7. What Scale Are We Calibrating To?

We are calibrating `a` to the same effect-size scale as `beta_hat`, `mu_0`, and `mu_1`.

Important clarification:

- "standardized effect scale" does **not** mean the empirical distribution of effects should resemble `N(0,1)`
- it only means the coefficients are expressed in the units induced by the summary-stat conversion / standardization

Therefore, the fact that approximate `mu_1` values from the real-data runs do **not** look standard normal is not a problem.

The nonpathological loci in the current 10-locus case study suggest a rough conditional effect scale around:

```text
|mu_1| ~ 0.1 to 0.3
```

with one clear pathological outlier locus (`MEI1`) showing implausibly huge values that should not be used to set the prior scale.

## 8. How to Compute beta_hat and shat2 Without Running SuSiE First

This is the key implementation point.

We do **not** need to run baseline SuSiE in order to obtain `beta_hat`.

`beta_hat` is a marginal summary-statistic object derived directly from the input data.

### 8.1 If using the `z`-only RSS path

For each SNP:

```text
adj_j = (n_j - 1) / (z_j^2 + n_j - 2)
z_tilde_j = sqrt(adj_j) * z_j
beta_hat_j = z_tilde_j / sqrt(n_j - 1)
shat2_j = 1 / (n_j - 1)
```

Equivalent compact form:

```text
beta_hat_j = z_j * sqrt(adj_j / (n_j - 1))
```

This is the standardized-effect-scale version implied by the current `susine_rss()` implementation.

### 8.2 If using `slope` and `slope_se`

Let `z_j = slope_j / slope_se_j`.

Then:

```text
adj_j = (n_j - 1) / (z_j^2 + n_j - 2)
beta_hat_j = slope_j / sqrt(adj_j)
shat2_j = slope_se_j^2 / adj_j
```

These are the effect estimate and sampling variance on the same adjusted scale implied by the `susine_rss()` conversion.

## 9. Decision on Using slope / slope_se vs z

### Decision

Prefer `slope` and `slope_se` when they are available and trustworthy, but only route the full model through the original-scale branch if `var_y` is available or can be justified.

### Why

The directional prior `mu_0` is fundamentally an effect-size prior, not a z-score prior. So `slope` / `slope_se` are closer to the object of interest.

However, the current `susine_rss()` original-scale branch requires `var_y`. Without `var_y`, we cannot cleanly claim to know the original effect scale used by the model.

## 10. How to Check Whether var_y Is Effectively 1

### Main rule

If the phenotype was already standardized or inverse-normal transformed upstream, `var_y` may be close enough to 1 that using `var_y = 1` is appropriate.

### Proposed diagnostic

For each SNP where the needed inputs are available, estimate:

```text
R2_j = z_j^2 / (z_j^2 + n_j - 2)
```

If genotype variance `var(x_j)` is known or can be approximated, then:

```text
var_y_hat_j = slope_j^2 * var(x_j) / R2_j
```

Equivalent expression using `slope_se_j`:

```text
var_y_hat_j = slope_se_j^2 * (n_j - 1) * var(x_j) / (1 - R2_j)
```

where:

- `var(x_j)` can be computed from the actual genotype matrix if available
- if only MAF is available, a dosage-scale approximation is:

```text
var(x_j) ~= 2 * maf_j * (1 - maf_j)
```

### Interpretation

- If the distribution of `var_y_hat_j` is tightly centered near 1 across loci, then using `var_y = 1` is defensible.
- If it is far from 1 or highly inconsistent across SNPs, flag the issue.

### Important limitation

If we do **not** know genotype variance per SNP, then `slope`, `slope_se`, and `z` alone are not enough to recover `var_y` on the raw original phenotype scale.

So the decision rule is:

- if `var_y_hat` is estimable and close to 1 -> route `var_y = 1`
- if `var_y_hat` is estimable and not close to 1 -> flag, and do not pretend the original scale is recovered unless actual `var_y` is available
- if `var_y_hat` is not estimable because genotype variance is missing -> fall back to the standardized `z` path for the actual SuSiNE runs

## 11. Estimating the Annotation Scale c Across Loci

### Core model

After constructing `a_j`, estimate the annotation scale using:

```text
beta_hat_j ~ N(c * a_j, shat2_j + tau^2)
```

This is the same unconditional scaling logic already implemented in the package.

### Recommended procedure

1. Build transformed annotations `a_j` for all SNPs across all 10 loci.
2. Compute `beta_hat_j` and `shat2_j` directly from summary data.
3. Estimate `c_l` separately within each locus.
4. Estimate a pooled `c0` across all loci.
5. Compare pooled and per-locus results.
6. Use `c0` only as the center of a grid, not as a single truth.

### Why estimate per-locus c_l as well?

Because one pooled estimate can be distorted by:

- one pathological locus
- different locus-specific signal strengths
- weak or noisy annotation-track alignment

The per-locus estimates act as a stability diagnostic.

### Recommended summary of c

Compute:

- pooled `c0`
- vector of per-locus `c_l`
- median of `c_l`
- MAD or IQR of `c_l`

If pooled `c0` and median `c_l` are reasonably similar, that supports a shared-scale assumption.

If they disagree strongly, widen the grid and avoid strong conclusions about a single "best" scale.

## 12. Is a Shared c Across Loci Reasonable?

### Decision

Yes, as a first-pass approximation, provided that the loci all come from the same overall analysis regime:

- same tissue
- same phenotype scaling
- same summary-stat generation pipeline
- same AlphaGenome output type and aggregation rule

### Caveat

Shared `c` is not assumed to be exact. It is just a practical working model for finding the right ballpark.

That is why the final SuSiNE runs should still use a small grid around the estimated center.

## 13. Recommended c Grid for SuSiNE Runs

Let `c0` be the pooled center.

Recommended default:

```text
c_grid = {0, 0.5*c0, 1*c0, 2*c0, 4*c0}
```

If the pooled estimate is very noisy or near zero, use a manually widened exploratory grid around the rough expected effect scale instead.

If sign reversals are biologically plausible, optionally include mirrored values:

```text
{-4*c0, -2*c0, -1*c0, -0.5*c0, 0, 0.5*c0, 1*c0, 2*c0, 4*c0}
```

If sign reversals are not plausible, constrain to nonnegative `c`.

## 14. Practical Starting Ballpark from Current Real-Data Results

From the current vanilla `susie_rss` outputs:

- nonpathological loci suggest approximate conditional effects around `|mu_1| ~ 0.1 to 0.3`
- one locus (`MEI1`) shows extreme apparent values (`|mu_1| > 20`) that should be treated as pathological and excluded from scale calibration

So on the standardized-effect path, the prior should aim for **modest conditional effects**, not huge ones.

Practical interpretation:

- if RMS-normalized `a` has top `|a|` around `2 to 3`
- then a plausible `c` center may be somewhere around `0.05 to 0.12`
- because that yields top `|mu_0|` around `0.1 to 0.3`

This is a sanity-check range, not a hard constraint.

## 15. Visual Diagnostics to Produce

The HPC workflow should generate at least the following plots.

### 15.1 Distribution of beta_hat across all loci

Create a histogram of pooled `beta_hat`.

Recommended variants:

- all pooled SNPs
- pooled SNPs after trimming top and bottom 0.5%
- faceted by locus

Goal:

- understand the effect-size scale being targeted
- detect pathologies / outliers
- confirm that most loci live in a modest range while a few may be extreme

### 15.2 Distribution of transformed annotation a

Plot histograms or density plots for:

- raw signed quantile
- inverse-normal transformed `a_raw`
- clipped `a_clip`
- final RMS-normalized `a`

Goal:

- verify that the transform behaves as intended
- confirm that winsorization is not too aggressive or too weak

### 15.3 Per-locus c_l estimates

Create a visual summary of `c_l` for each locus.

Recommended display:

- point-range plot of `c_l`
- add pooled `c0` as a horizontal reference line
- optionally color loci by whether they were flagged as pathological

Goal:

- see whether the shared-scale assumption is plausible
- identify loci that strongly disagree with the pooled center

### 15.4 Estimated var_y diagnostics

If genotype variance is available and `var_y_hat` can be computed:

- histogram of `var_y_hat`
- per-locus boxplots / violin plots of `var_y_hat`
- scatter of `var_y_hat` vs MAF or `|z|`

Goal:

- assess whether `var_y = 1` is defensible
- identify whether the original-scale branch is actually usable

## 16. Locus Exclusion / Flagging Rules

For scale estimation, exclude or at least flag loci with any of the following:

- absurdly large inferred conditional effects compared with the rest of the study
- obvious summary-stat pathologies
- severe LD instability
- near-degenerate or flat annotations
- strong disagreement between sign of annotation and marginal effects at the lead signals

For the current case study, `MEI1` should be treated as a flagged locus for the initial scale-estimation pass.

## 17. Step-by-Step Execution Guide for the HPC Agent

### Phase A: assemble inputs

1. Load the 10-locus summary-stat tables.
2. Load the AlphaGenome signed quantile annotations for the same variants.
3. Harmonize variant IDs, alleles, and sign conventions.
4. Build one master table with:
   - `locus_id`
   - `variant_id`
   - `z`
   - `n`
   - `slope` if available
   - `slope_se` if available
   - genotype variance proxy if available (`var_x` or MAF)
   - signed AlphaGenome quantile score

### Phase B: transform annotation

5. Compute clipped quantile `q*`.
6. Compute inverse-normal transform `a_raw`.
7. Winsorize to `a_clip`.
8. RMS-normalize within locus to get final `a`.
9. Save intermediate columns for diagnostics.

### Phase C: compute effect-scale summaries

10. Compute `beta_hat` and `shat2` using:
    - `z` formulas for the standardized path
    - `slope` / `slope_se` formulas where applicable
11. If genotype variance is available, compute `var_y_hat`.
12. Summarize `var_y_hat` overall and by locus.
13. Decide whether `var_y = 1` is justified or whether the original-scale path should be flagged as unsupported.

### Phase D: estimate annotation scale

14. Estimate `c_l` and `tau_l^2` within each locus using:

```text
beta_hat_j ~ N(c_l * a_j, shat2_j + tau_l^2)
```

15. Estimate pooled `c0` and `tau0^2` across all loci.
16. Compute robust summaries of `c_l`:
    - median
    - IQR
    - MAD
17. Compare pooled `c0` to median `c_l`.

### Phase E: plotting

18. Plot pooled `beta_hat` histogram.
19. Plot annotation transform diagnostics.
20. Plot per-locus `c_l` with pooled `c0`.
21. Plot `var_y_hat` diagnostics if computable.

### Phase F: choose run grid

22. Define `c_grid` around `c0`.
23. Build `mu_0 = c * a` for each `c` in the grid.
24. Run SuSiNE over that grid.
25. Treat `c0` as a center, not as a single ground-truth value.

## 18. Deliverables Expected from the HPC Agent

The HPC-side implementation should return:

1. A cleaned merged table of all loci and variants used for scale estimation
2. A table of transformed annotations per variant
3. A table of `beta_hat`, `shat2`, and, if possible, `var_y_hat`
4. A per-locus table of `c_l` and `tau_l^2`
5. A pooled estimate table for `c0` and `tau0^2`
6. The requested plots
7. A recommended `c_grid` for downstream SuSiNE runs
8. A short note listing any loci that were flagged or excluded and why

## 19. Final Operational Decision

The working plan is:

- use signed quantile annotations
- inverse-normal transform them
- clip and RMS-normalize within locus
- estimate scale using `beta_hat`, not `b_hat`
- pool across loci to get a rough `c0`
- inspect per-locus `c_l` for stability
- use a grid around `c0`
- only use the original-scale `slope` / `slope_se` route if `var_y` is available or can be credibly shown to be approximately 1
- otherwise use the standardized `z` path and treat the prior as operating on the standardized-effect scale

This keeps the annotation scaling coherent, auditable, and robust to the fact that the true mapping from AlphaGenome outputs to eQTL effect sizes is only approximate.
