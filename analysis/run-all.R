#--------------------------------------------
#
# AVENIR-spatial-ate-code-check/analysis/run-all.R
#
# Orchestrator: renders every analysis script in pipeline order.
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

render_script("01-data-prep.Rmd")

# 02 — Spatial model fitting and posterior probability simulation ----------
# NOTE: 02A fits spaMM models and saves them to data/. It can take ~90 min
# on first run; subsequent runs load cached .rds files when
# refit_models_flag = FALSE (set in 0-config.R).

render_script("02A-spatial-models-fit.Rmd")
render_script("02B-posterior-prob-sim.Rmd")

# 03 — CSI-level spatial analysis -----------------------------------------

render_script("03-csi-analysis.qmd")

# 04 — Variance decomposition / overdispersion ----------------------------

render_script("04-overdispersion.qmd")

# 05 — Results summary ----------------------------------------------------

render_script("05-results.qmd")

# 06 — Manuscript figures -------------------------------------------------

render_script("06A-fig1-context_map.qmd")
render_script("06B-fig2-primary_maps.qmd")

# scripts that are not pipeline steps are intentionally excluded here

message("\n--- run-all.R complete ---")
