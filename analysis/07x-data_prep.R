#--------------------------------------------
#
# AVENIR-spatial-ate/analysis/07x-data_prep.R
#
# Builds the two phase-level inputs 07-targeting-analysis.qmd reads: grappe x
# phase mortality, and grappe x round dose counts. run-all.R sources this
# immediately before rendering 07; run it by hand first if rendering 07 alone.
#
# It is sourced rather than rendered, so run-all.R sources it into its own
# environment — the rm(list = ls()) below would otherwise wipe the caller's
# globals.
#
# Inputs:  AVENIR extracts under `final_data` (Box path, set in 0-config.R):
#            AVENIR_grappe_phase_2023-07-31.rds
#            AVENIR_child_census_visit_2025-08-25.rds
# Outputs: data/grappe_mortality_data_by_phase.rds
#          data/grappe_dose_counts_by_phase.rds
#
#--------------------------------------------

rm(list = ls())

# 0-config.R is sourced for `final_data` (the Box path to the AVENIR extracts)
# and the tidyverse attachments. 0-functions.R is not needed here.
source(here::here("analysis", "0-config.R"))

# Phase-level mortality data: one row per grappe × study phase
grappe_phase_mort_raw <-
  readRDS(paste0(final_data,"AVENIR_grappe_phase_2023-07-31.rds"))

child_census_raw <-
  readRDS(paste0(final_data,"AVENIR_child_census_visit_2025-08-25.rds"))

# Only the columns that survive the final select() are constructed here. In
# particular district is *not* harmonized to the shapefile spelling: this table
# is joined to grappe_mortality_data.rds by site_code inside
# build_prediction_table(), and district comes from there already harmonized by
# 01-data-prep.Rmd. isTreated_alltx and the `dummy` covariate belong to the
# primary analysis and have no role in the targeting addendum.
grappe_phase_mort <-
  grappe_phase_mort_raw %>%
  mutate(
    tx = case_when(
      tx == "1 - Azithro 1-59m" ~ "Child Azithro",
      tx == "2 - Azithro 1-11m" ~ "Infant Azithro",
      tx == "3 - Placebo 1-59m" ~ "Placebo"),

    # Primary analysis: "Child Azithro" vs. all others
    isTreated_child = case_when(
      tx == "Child Azithro" ~ "Treatment",
      .default = "Control"
    ),

    tx              = factor(tx,              levels = c("Child Azithro", "Infant Azithro", "Placebo")),
    isTreated_child = factor(isTreated_child, levels = c("Control", "Treatment")),

    # shorten adaptive randomization labels
    adaptive_rand = factor(adaptive_rand,
                           levels = c("Random Allocation 1",
                                      "Random Allocation 2",
                                      "Random Allocation 3"),
                           labels = c("Randomization 1",
                                      "Randomization 2",
                                      "Randomization 3")),
    phase = case_when(
      phase == "0-6" ~ "Round 1",
      phase == "6-12" ~ "Round 2",
      phase == "12-18" ~ "Round 3",
      phase == "18-24" ~ "Round 4",
      TRUE ~ as.character(NA)
    )
  ) %>%
  filter(adaptive_rand == "Randomization 1") %>%
  dplyr::select(site_code, grappe, phase, tx, ends_with("_1to59"), isTreated_child)

# The phase recode above maps the four AVENIR rounds and sends anything else to
# NA. The dose extract below is hard-filtered to the same four rounds, so an
# unrecognised phase would silently inflate the eligible-child denominator
# without touching the dose count. Fail loudly instead.
stopifnot(!any(is.na(grappe_phase_mort$phase)))

saveRDS(grappe_phase_mort, here("data", "grappe_mortality_data_by_phase.rds"))

dose_counts <- child_census_raw %>%
  dplyr::filter(is_administered == "Yes") %>%
  dplyr::filter(phase %in% c("0", "6", "12", "18")) %>%
  dplyr::filter(site_code %in% unique(grappe_phase_mort$site_code)) %>%
  dplyr::select(phase, site_code, grappe) %>%
  group_by(phase, site_code, grappe) %>%
  dplyr::summarise(n_doses = n(), .groups = "drop") %>%
  dplyr::mutate(phase = case_when(
    phase == "0" ~ "Round 1",
    phase == "6" ~ "Round 2",
    phase == "12" ~ "Round 3",
    phase == "18" ~ "Round 4",
    TRUE ~ as.character(NA)
  ))

saveRDS(dose_counts, here("data", "grappe_dose_counts_by_phase.rds"))
