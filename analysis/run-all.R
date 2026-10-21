#--------------------------------------------
#
# AVENIR-spatial-ate/analysis/run-all.R
#
# Orchestrator: renders every analysis script in pipeline order.
#
# Run from the project root. Scripts write their intermediate objects to
# data/ and their figures to results/, and each step below reads what the
# steps above it wrote, so the order here is a hard dependency order, not a
# convenience.
#
# Runtime is dominated by two steps, both controlled by flags in 0-config.R:
# model fitting in 02A (~90 min) and posterior simulation in 02B (~13 min per
# model). With the flags at their defaults (FALSE) both load cached .rds files
# instead and the whole pipeline runs in minutes.
#
# Scripts prefixed `x-` are side analyses and slide-deck figures. They are not
# pipeline steps and are deliberately excluded here; see README.md.
#
#--------------------------------------------

renv::restore()

library(here)

# Helper: render .Rmd (rmarkdown) or .qmd (quarto) by file extension.
render_script <- function(filename) {
  path <- here("analysis", filename)
  message("\n--- Rendering: ", filename, " ---")
  if (grepl("\\.qmd$", filename)) {
    quarto::quarto_render(path)
  } else {
    rmarkdown::render(path)
  }
}

# 01 — Data preparation ---------------------------------------------------
# Reads the AVENIR extracts from Box (path `final_data` in 0-config.R) and
# writes the grappe-level analytic dataset and the prediction grid to data/.

render_script("01-data-prep.Rmd")

# 02 — Spatial model fitting and posterior probability simulation ----------
# 02A fits the spaMM models (M1-M5) and saves them to data/, then predicts
# both arms over the prediction grid into data/preds_ATE.rds. It can take
# ~90 min on first run; subsequent runs load the cached .rds files when
# refit_models_flag = FALSE (set in 0-config.R).
# 02B simulates the posterior probability that the ATE is beneficial at each
# location; controlled by re_simulate in 0-config.R.

if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = TRUE)
}

render_script("02A-spatial-models-fit.Rmd")
render_script("02B-posterior-prob-sim.Rmd")

# 03 — CSI-level spatial analysis -----------------------------------------
# Repeats the primary model at the CSI (health facility catchment) scale.

render_script("03-csi-analysis.qmd")

# 04 — Variance decomposition / overdispersion ----------------------------
# Partitions mortality variance across grappe, CSI and district scales.

render_script("04-overdispersion.qmd")

# 05 — Results summary ----------------------------------------------------
# Consolidates the headline estimates cited in the manuscript text.

render_script("05-results.qmd")

# 06 — Sensitivity: adequacy of the linear ATE plane ----------------------
# Loads the cached fits (M1-M5) and cached prediction surfaces; fits nothing.
# Compares four more flexible treatment surfaces to the primary linear plane.

render_script("06-check-linearity.qmd")

# 07 — Targeting-rule addendum --------------------------------------------
# 07x-data_prep.R builds the two phase-level inputs 07 reads (grappe x phase
# mortality, grappe x round doses) directly from the Box extracts; it is a
# plain script, so it is sourced rather than rendered. It is sourced into its
# own environment because, like every other script here, it opens with
# rm(list = ls()) — into the global environment that would delete
# render_script() and stop the pipeline at the next call.
# 07 then evaluates whether ranking villages by predicted mortality recovers
# the benefit of ranking by predicted treatment effect, and saves the objects
# Figure 3 is drawn from. Its two rebuild flags (rebuild_targeting_preds_flag,
# rebuild_targeting_boot_flag) are in 0-config.R; both cached objects are
# stamped with the settings they were built under and 07 stops rather than
# renders if the stamp no longer matches.

source(here("analysis", "07x-data_prep.R"), local = new.env())
render_script("07-targeting-analysis.qmd")

# 08 — Manuscript figures -------------------------------------------------
# Pure presentation layer: these attach no models and fit nothing, they read
# the cached objects written above. 08C must follow 07, which writes
# data/targeting_fig3_objects.rda.

render_script("08A-fig1-context_map.qmd")
render_script("08B-fig2-primary_maps.qmd")
render_script("08C-fig3-targeting_curves.qmd")
render_script("08D-supp_mtrl.qmd")

message("\n--- run-all.R complete ---")
