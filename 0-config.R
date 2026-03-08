##############################################
# Spatial analysis of mortality rates in Dosso
# Secondary data analysis of AVENIR study data

# Configure data paths & load required packages

# by Mia Navarro (mia.navarro@ucsf.edu)
##############################################

##############################################
# Set wd
##############################################
wd_path = here::here()
# setwd(wd_path)

##############################################
# Load required packages
##############################################
# library(assertthat)
# library(cowplot)
library(dplyr)
# library(foreign)
# library(gee)
# library(geeM) #> Warning: package 'geeM' was built under R version 3.6.2
# library(geepack) #> Warning: package 'geepack' was built under R version 3.6.2
# library(geosphere)
# library(geodist)
# library(ggplot2) # tidyverse data visualization package
# library(ggraph)
# library(ggmap)
library(ggpubr)
library(ggspatial)
# library(gratia)
# library(grid)
# library(gridExtra)
# library(gtsummary)
library(here)
# library(igraph)
# library(leaflet) # for interactive maps
library(lme4)
# library(MASS)
# library(Matrix)
# library(mgcv)
# library(mgcViz)
library(paletteer)
# library(parallel) # use mc.cores=7
# library(parallelly)
library(patchwork)
# library(psych)
library(purrr)
library(raster)
# library(readr)
# library(readxl)
# library(recipes)
library(renv)
# library(rgl)
library(scales)
library(sf)
# library(sjPlot)
library(spaMM)
# library(spData)
# library(spDataLarge)
# library(splines)
# library(tableone)
# library(terra)
# library(tidymv)
library(tidyverse)
# library(tmap)    # for static and interactive maps
library(tictoc)

##############################################
# Define data paths
##############################################
# raw_data_path = paste0(here::here(), "/data/untouched-raw/")
final_data <- "~/Library/CloudStorage/Box-Box/AVENIR/Data/Mortality/rounds_1to5/final/"
shp_data <- "~/Library/CloudStorage/Box-Box/AVENIR/Publications/Primary Mortality/avenir-primary/data/ner_adm_shapefiles/"

# generated data paths
my_temp_data_path = paste0(here::here(), "/data/temp/")
my_final_data_path = paste0(here::here(), "/data/final/")

# load base functions
# source(paste0(here::here(),"/analysis/0-functions-processing.R"))
# source(paste0(here::here(),"/analysis/0-functions-analysis.R"))

##############################################
# Define figure and table paths
##############################################
figure_path = paste0(here::here(), "/results/figs/")
table_path = paste0(here::here(), "/results/tables/")
