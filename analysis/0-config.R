# Spatial analysis of mortality rates in Dosso
# Secondary data analysis of AVENIR study data

# Configure data paths & load required packages

# by Ariktha Srivathsan (ariktha.srivathsan@ucsf.edu)

# Load required packages ----

library(colorspace)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(ggspatial)
library(ggtext)
library(here)
library(lme4)
library(paletteer)
library(parallel)
library(patchwork)
library(purrr)
library(raster)
library(renv)
library(scales)
library(sf)
library(spaMM)
# splines must be attached, not called as splines::ns() inside a spaMM formula:
# spaMM's .expand_multIMRFs() compares each term's head against "multIMRF", and a
# namespaced call head is length 3, which errors under R >= 4.2.
library(splines)
library(terra)
library(tidyverse)
library(tictoc)
library(viridis)

# Define file paths ----

final_data <- "~/Library/CloudStorage/Box-Box/AVENIR/Data/Mortality/rounds_1to5/final/"

figure_path <- paste0(here::here(), "/results/figs/")
manuscript_figs_path <- paste0(here::here(), "/results/manuscript/")
shpfile_path <- paste0(here::here(), "/data/untouched/")


if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = TRUE)
}

if (!dir.exists(manuscript_figs_path)) {
  dir.create(manuscript_figs_path, recursive = TRUE)
}


# Set up ----

## Geospatial parameters
global_crs <- 4326 # WGS84, for lat/lon coordinates
proj_crs <- 32631 # UTM zone 31N, for Niger, for projected distances in meters

## Parallel processing
avail_thr <- parallel::detectCores(logical = FALSE) - 1L

## Grid dimensions for prediction layers
grid_res <- 200

## Simulations for posterior predictions
nsim_input <- 1000

## Flags to control if models are re-fit
refit_models_flag <- TRUE # Set to TRUE to refit models, which can take a long time. Set to FALSE to load previously fitted models.
refit_csi_models_flag <- TRUE
refit_overdispersion_models_flag <- TRUE
## Refits the primary model on district subsets in x-district-exclusion-check.qmd
## (~4-5 min per subset). The full-data fit is always loaded, never refitted there.
refit_exclusion_models_flag <- TRUE

## Flags to control if posterior probability simulations are re-run in script 2B
re_simulate <- TRUE # Set TRUE to re-run simulations (slow: ~13 min per model)

## Flags to control if the targeting inputs in script 07 are rebuilt.
## Both cached objects carry an attribute stamp of the settings they were built
## under, and 07 stops rather than renders if the stamp no longer matches, so
## leaving these FALSE cannot silently produce stale numbers.
rebuild_targeting_preds_flag <- TRUE # Set TRUE to re-predict both arms at every grappe (loads a ~54 MB model per specification)
rebuild_targeting_boot_flag  <- TRUE # Set TRUE to redraw the predictive bootstrap (slow: the joint covariance over villages x arms)
