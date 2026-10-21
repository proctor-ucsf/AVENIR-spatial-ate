# AVENIR-spatial-ate

Estimation of fine-scale spatial heterogeneity in treatment effects in a large cluster-randomized trial

This repository contains the full analysis code for a secondary, spatial
analysis of the AVENIR cluster-randomized trial of biannual mass azithromycin
distribution in the Dosso region of Niger. The analysis fits geospatial models 
to village-level mortality counts to estimate continuous surfaces of
mortality and of treatment effect across the trial area, and asks whether the
effect of azithromycin varies systematically over space.

## Data availability

Village geographic coordinates data required to reproduce this analysis are 
protected and are not available to protect participant confidentiality.

## Repository structure

```
analysis/     all analysis code (see the pipeline table below)
data/         inputs and intermediate objects  [gitignored]
results/
  figs/       exploratory and diagnostic figures  [gitignored]
  manuscript/ manuscript and supplementary figures  [gitignored]
renv/         renv library metadata
renv.lock     pinned package versions
```

## Reproducing the analysis

Requires R 4.5.2 (the version `renv.lock` was built against), Quarto, and
access to the data described above.

```bash
git clone https://github.com/proctor-ucsf/AVENIR-spatial-ate.git
cd AVENIR-spatial-ate
```

Then, from R in the project root:

```r
source("analysis/run-all.R")    # render the full pipeline in order
```

[`analysis/run-all.R`](analysis/run-all.R) renders every pipeline script in
dependency order. Each script writes its intermediate objects to `data/` and
its figures to `results/`, and later scripts read what earlier ones wrote.

### Configuration

[`analysis/0-config.R`](analysis/0-config.R) holds the data path, the package
attachments, and every tunable parameter — the two coordinate reference systems,
the prediction grid resolution, the number of posterior simulations, and the
refit flags.

The refit flags matter quite a bit. Fitting the spatial models takes roughly 90
minutes; the posterior simulation takes about 13 minutes per model. All
flags default to `FALSE`, which loads cached `.rds` objects instead, so a
normal run takes minutes. Set a flag to `TRUE` only when you intend to refit
and overwrite the cache.

## Pipeline

Scripts run in numeric order. `0-` prefixed files are shared setup and are
sourced by the others rather than run on their own.

| Script | What it does |
|---|---|
| [`0-config.R`](analysis/0-config.R) | Data paths, package attachments, geospatial parameters, refit flags |
| [`0-functions.R`](analysis/0-functions.R) | Shared helpers: bootstrap scaffolding, prediction and rasterisation utilities |
| [`01-data-prep.Rmd`](analysis/01-data-prep.Rmd) | Builds the grappe-level analytic dataset and the spatial prediction grid |
| [`02A-spatial-models-fit.Rmd`](analysis/02A-spatial-models-fit.Rmd) | Fits every spatial model (Poisson + Matérn, via `spaMM`); predicts both arms over the grid |
| [`02B-posterior-prob-sim.Rmd`](analysis/02B-posterior-prob-sim.Rmd) | Simulates the predictive distribution; maps P(ATE < 0) |
| [`03-csi-analysis.qmd`](analysis/03-csi-analysis.qmd) | Repeats the spatial analysis at the CSI (health facility catchment) scale |
| [`04-overdispersion.qmd`](analysis/04-overdispersion.qmd) | Partitions variance across grappe, CSI and district scales; ICCs with bootstrap CIs |
| [`05-results.qmd`](analysis/05-results.qmd) | Consolidates every number the manuscript text cites |
| [`06-check-linearity.qmd`](analysis/06-check-linearity.qmd) | Sensitivity analysis: is the primary model's linear ATE plane adequate? |
| [`07x-data_prep.R`](analysis/07x-data_prep.R) | Builds the phase-level mortality and dose inputs `07` needs |
| [`07x-functions.R`](analysis/07x-functions.R) | Helpers for the targeting analysis |
| [`07-targeting-analysis.qmd`](analysis/07-targeting-analysis.qmd) | Targeting addendum: does ranking by mortality recover the benefit of ranking by effect? |
| [`08A-fig1-context_map.qmd`](analysis/08A-fig1-context_map.qmd) | Figure 1: study area context map |
| [`08B-fig2-primary_maps.qmd`](analysis/08B-fig2-primary_maps.qmd) | Figure 2: the six primary maps |
| [`08C-fig3-targeting_curves.qmd`](analysis/08C-fig3-targeting_curves.qmd) | Figure 3: targeting curves |
| [`08D-supp_mtrl.qmd`](analysis/08D-supp_mtrl.qmd) | Supplementary table and figures |

The `08` scripts are a pure presentation layer: they attach no models and
estimate nothing, they read cached objects and compose figures.

## Side analyses

Scripts prefixed `x-` are not pipeline steps and are excluded from
`run-all.R`. They answer specific questions raised in review.

| Script | What it does |
|---|---|
| [`x-mordor-question.qmd`](analysis/x-mordor-question.qmd) | Does the ATE differ in the three districts that overlapped the earlier MORDOR trial? Approaches it by *adjusting*: fuzzy-matches villages across the two trials and adds a MORDOR indicator |
| [`x-district-exclusion-check.qmd`](analysis/x-district-exclusion-check.qmd) | The complementary approach: *removes* those districts, refits the primary model, and compares the surfaces |
| [`x-assess-spatial-scale-of-correlation.qmd`](analysis/x-assess-spatial-scale-of-correlation.qmd) | Standalone assessment of the spatial scale of correlation (semivariogram, Matérn range, variance partitioning). **Does not currently run** — it sources two helper files that are not in this repo |
| [`x-slide_deck-figs.qmd`](analysis/x-slide_deck-figs.qmd) | Figures for talks, at presentation proportions |
