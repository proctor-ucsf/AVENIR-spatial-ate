# AVENIR-spatial-ate

Estimating spatially continuous under-5 mortality rates and average treatment
effects (ATE) in the AVENIR trial.

This repository contains the full analysis code for a secondary, spatial
analysis of the AVENIR cluster-randomized trial of biannual mass azithromycin
distribution in the Dosso region of Niger. Rather than reporting a single
pooled effect, the analysis fits geostatistical models to grappe-level
(community cluster) mortality counts to estimate continuous surfaces of
mortality and of treatment effect across the trial area, and asks whether the
effect of azithromycin varies systematically over space.

A companion analysis asks a programmatic question that follows from those
surfaces: if a programme cannot treat everywhere, does prioritising places by
predicted **mortality** recover most of the benefit of prioritising by
predicted **treatment effect**?

- **Analysis:** UCSF Francis I. Proctor Foundation
- **Contact:** Ariktha Srivathsan (ariktha.srivathsan@ucsf.edu), Mia Navarro (mia.navarro@ucsf.edu)
- **License:** CC0 1.0 Universal (see [LICENSE](LICENSE))

## Data availability

**The trial data are not in this repository and cannot be distributed with it.**

All scripts read the AVENIR mortality extracts from a path set as `final_data`
in [`analysis/0-config.R`](analysis/0-config.R), which points at a Box folder
available to the study team. Running the pipeline from scratch therefore
requires access to those extracts; everything downstream of them is fully
reproducible from this code.

Three other inputs live under `data/`, which is gitignored:

| Path | Source |
|---|---|
| `data/untouched/ner_admin_boundaries/` | Niger administrative boundaries (ADM1, ADM2), from the [Humanitarian Data Exchange](https://data.humdata.org/dataset/cod-ab-ner) |
| `data/untouched/from-mordor/` | MORDOR trial village GPS, used only by the MORDOR side analyses |
| `data/*.rds` | Intermediate objects written by the pipeline itself |

Requests for the underlying trial data should be directed to the AVENIR study
team.

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
renv::restore()                 # install the pinned package versions
source("analysis/run-all.R")    # render the full pipeline in order
```

[`analysis/run-all.R`](analysis/run-all.R) renders every pipeline script in
dependency order. Each script writes its intermediate objects to `data/` and
its figures to `results/`, and later scripts read what earlier ones wrote, so
the order is not optional.

### Configuration

[`analysis/0-config.R`](analysis/0-config.R) holds the data path, the package
attachments, and every tunable parameter — the two coordinate reference systems,
the prediction grid resolution, the number of posterior simulations, and the
refit flags.

The refit flags matter most. Fitting the spatial models takes roughly **90
minutes**; the posterior simulation takes about **13 minutes per model**. All
flags default to `FALSE`, which loads cached `.rds` objects instead, so a
normal run takes minutes. Set a flag to `TRUE` only when you intend to refit
and overwrite the cache.

| Flag | Controls |
|---|---|
| `refit_models_flag` | the spaMM fits in `02A` |
| `refit_csi_models_flag` | the CSI-level fits in `03` |
| `refit_overdispersion_models_flag` | the variance-partitioning fits in `04` |
| `re_simulate` | the posterior simulation in `02B` |
| `rebuild_targeting_preds_flag` | the per-village predictions in `07` |
| `rebuild_targeting_boot_flag` | the predictive bootstrap in `07` |
| `refit_exclusion_models_flag` | the subset refits in `x-district-exclusion-check` |

The two `07` caches carry an attribute stamp of the settings they were built
under, and `07` stops rather than renders if a stamp no longer matches — so
leaving those flags at `FALSE` cannot silently produce stale numbers.

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

### The primary model

Model comparison in `02A` and `06` is **not** model selection. The primary
specification is fixed in advance: generalized kriging with a linear
treatment-by-coordinate interaction (`fitmort_child_gen_krig_int`), which
constrains the log rate ratio to a plane in latitude and longitude. The more
flexible specifications — additive splines at two knot counts, a quadratic
treatment contrast, an arm-specific Matérn deviation field — exist to bound how
much curvature a plane could be missing. It stays primary regardless of what
they show.

## Side analyses

Scripts prefixed `x-` are not pipeline steps and are excluded from
`run-all.R`. They answer specific questions raised in review or in analysis
meetings and are kept for the record.

| Script | What it does |
|---|---|
| [`x-mordor-question.qmd`](analysis/x-mordor-question.qmd) | Does the ATE differ in the three districts that overlapped the earlier MORDOR trial? Approaches it by *adjusting*: fuzzy-matches villages across the two trials and adds a MORDOR indicator |
| [`x-district-exclusion-check.qmd`](analysis/x-district-exclusion-check.qmd) | The complementary approach: *removes* those districts, refits the primary model, and compares the surfaces |
| [`x-assess-spatial-scale-of-correlation.qmd`](analysis/x-assess-spatial-scale-of-correlation.qmd) | Standalone assessment of the spatial scale of correlation (semivariogram, Matérn range, variance partitioning). **Does not currently run** — it sources two helper files that are not in this repo |
| [`x-slide_deck-figs.qmd`](analysis/x-slide_deck-figs.qmd) | Figures for talks, at presentation proportions |

## Citation

Manuscript in preparation. Please contact the authors before citing this
repository.
