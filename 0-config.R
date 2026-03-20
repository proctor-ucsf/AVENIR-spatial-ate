# Spatial analysis of mortality rates in Dosso
# Secondary data analysis of AVENIR study data

# Configure data paths & load required packages

# by Mia Navarro (mia.navarro@ucsf.edu); Ariktha Srivathsan (ariktha.srivathsan@ucsf.edu)

# Load required packages ----

library(dplyr)
library(ggpubr)
library(ggspatial)
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
library(terra)
library(tidyverse)
library(tictoc)
library(viridis)

# Define file paths ----

final_data <- "~/Library/CloudStorage/Box-Box/AVENIR/Data/Mortality/rounds_1to5/final/"

figure_path = paste0(here::here(), "/results/figs/")
table_path = paste0(here::here(), "/results/tables/")

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
