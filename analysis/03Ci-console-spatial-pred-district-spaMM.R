# 03C.i-console-spatial-pred-district-spaMM

results_all <- list()
results_all$func <- rep("spaMM", 15)
results_all$analysis <- rep("2ry: by district", 15)
results_all$district <- unique(df_all$district)


######################################################################
# Identify the bounding boxes around districts -- IN results_all ####
######################################################################
results_all$b_boxgrid <- list()

for (i in 1:15) {
  
  b_boxgrid_temp <- st_make_grid(
    niger_adm2 %>% filter(adm_02 == results_all$district[i]),
    n = c(100, 100),
    what = "centers"
  )
  
  results_all$b_boxgrid[[i]] <- b_boxgrid_temp
  
}

######################################################################
# Fit spatial models for mortality rate with new binary tx vars -- IN results_all ####
######################################################################
print(Sys.time())
create_fitmort_d <- function(i) {
  temp_fitmort <- create_fitme(df = df_all %>%
                                 dplyr::filter(district == districts[i]))
  
  temp_fitmort_int <- create_fitme_int(df = df_all %>%
                                         dplyr::filter(district == districts[i]))
  
  results_all$fitmort[[i]] <- temp_fitmort
  results_all$fitmort_int[[i]] <- temp_fitmort_int
  
}

mclapply(1:15, create_fitmort_d, mc.cores = 6)


######################################################################
# Generate and store matern correlation plot -- IN results_all #### 
######################################################################
results_all$materncorr <- list()

create_materncorr_d <- function(i) {
  temp_materncorr <- plot_matern_corr(
    results_all$fitmort[[i]],
    paste0("avenir-materncorr-",results_all$district[i]))
  
  results_all$materncorr[[i]] <- temp_materncorr
  
}

mclapply(1:15, create_materncorr_d, mc.cores = 6)


######################################################################
# Create model predictions per tx arm -- started 5/23 3:19 PM ####
######################################################################
results_all$b_preds_mort_tx0 <- list()
results_all$b_preds_mort_tx1 <- list()

create_preds_d <- function(i) {
  temp_pred_tx0 <- get_grid_preds(
    input_grid = results_all$b_boxgrid[[i]],
    spamm_model_fit = results_all$fitmort[[i]],
    tx_group = unique(df_all$isTreated_child)[1]
  )
  
  temp_pred_tx1 <- get_grid_preds(
    input_grid = results_all$b_boxgrid[[i]],
    spamm_model_fit = results_all$fitmort[[i]],
    tx_group = unique(df_all$isTreated_child)[2]
  )
  
  results_all$b_preds_mort_tx0[[i]] <- temp_pred_tx0
  results_all$b_preds_mort_tx1[[i]] <- temp_pred_tx1
  
}

mclapply(1:15, create_preds_d, mc.cores = 6)


######################################################################
# Convert predictions into mapping elements ####
######################################################################
create_mapping_elements_d <- function(i) {
  temp_mort_tx0 <- b_preds_mort_tx0_d[[i]]
  temp_mort_tx1 <- b_preds_mort_tx1_d[[i]]
  
  temp_coords_tx0 <- st_coordinates(temp_mort_tx0)
  temp_coords_tx1 <- st_coordinates(temp_mort_tx1)
  
  temp_raster_tx0 <- points_to_raster(
    x = temp_coords_tx0[, 1],
    y = temp_coords_tx0[, 2],
    z = temp_mort_tx0$pred,
    mask1 = niger_adm2,
    crop1 = niger_adm2
  )
  
  temp_raster_tx1 <- points_to_raster(
    x = temp_coords_tx1[, 1],
    y = temp_coords_tx1[, 2],
    z = temp_mort_tx1$pred,
    mask1 = niger_adm2,
    crop1 = niger_adm2
  )
  
  temp_raster2_tx0 <- raster_to_tibble(temp_raster_tx0)
  temp_raster2_tx1 <- raster_to_tibble(temp_raster_tx1)
  
  results_all$b_preds_mort_coords_tx0[[i]] <- temp_coords_tx0
  results_all$b_preds_mort_coords_tx1[[i]] <- temp_coords_tx1
  
  results_all$b_preds_mort_raster_tx0[[i]] <- temp_raster_tx0
  results_all$b_preds_mort_raster_tx1[[i]] <- temp_raster_tx1
  
  results_all$b_preds_mort_raster2_tx0[[i]] <- temp_raster2_tx0
  results_all$b_preds_mort_raster2_tx1[[i]] <- temp_raster2_tx1
  
}

mclapply(1:15, create_mapping_elements_d, mc.cores = 6)


######################################################################
# Generate diff (ATE) raster2 element ####
######################################################################
create_diff_raster2_d <- function(i) {
  raster2_tx0 <- results_all$b_preds_mort_raster2_tx0
  raster2_tx1 <- results_all$b_preds_mort_raster2_tx1
  
  results_all$b_preds_mort_raster2_diff[[i]] <-
    raster2_tx0[[i]] %>%
    left_join(raster2_tx1[[i]], by = join_by(x, y)) %>%
    mutate(value = value.x - value.y) %>%
    rename(variable = variable.y) %>%
    dplyr::select(x, y, value, variable)
  
}

mclapply(1:15, create_diff_raster2_d, mc.cores = 6)


######################################################################
# Calculate uncertainty from simulations (`simulate()`). ####
######################################################################
# Steps:
#   Simulate for placebo + tx arm (foreach::doParallel)
#   Subtract from each other for ATE
#   Use sims_diar()

relevel(df_all$isTreated_child, ref = "Control")

create_newds_d <- function(i) {
  temp_newd2 <- b_boxgrid_d[[i]] %>%
    st_centroid() %>%
    st_as_sf(crs = 4326) %>%
    mutate(
      lon = sf::st_coordinates(.)[, 1],
      lat = sf::st_coordinates(.)[, 2],
      isTreated_child = "Treatment",
      site_code = 0,
      pt_year_1to59 = 1000) %>%
    st_drop_geometry()
  
  #site_code and pt_year_1to59 vals are dummy vals
  
  temp_newd3 <- temp_newd2 %>% mutate(isTreated_child = "Control")
  
  results_all$newd2[[i]] <- temp_newd2
  results_all$newd3[[i]] <- temp_newd3
  
}

nsim_input <- 1000

create_simulated_response_d <- function(i) {
  temp_sim_resp_tx0 <- simulate(
    fitmort_d[[i]],
    type = "predVar",
    variances = list(linPred = TRUE, disp = TRUE),
    nsim = nsim_input,
    newdata = newd3_d[[i]],
    seed = 12345
  )
  
  temp_sim_resp_tx1 <- simulate(
    fitmort_d[[i]],
    type = "predVar",
    variances = list(linPred = TRUE, disp = TRUE),
    nsim = nsim_input,
    newdata = newd2_d[[i]],
    seed = 12345
  )
  
  temp_sim_resp_ate <- temp_sim_resp_tx1 - temp_sim_resp_tx0
  
  results_all$simulated_responses_tx0[[i]] <- temp_sim_resp_tx0
  results_all$simulated_responses_tx1[[i]] <- temp_sim_resp_tx1
  results_all$simulated_responses_ate[[i]] <- temp_sim_resp_ate
  
}

mclapply(1:15, create_newds_d, mc.cores = 6)
mclapply(1:15, create_simulated_response_d, mc.cores = 6)


######################################################################
# Compute predicted ATE means, standard errors, and 95% CIs. ####
######################################################################

create_ate_stats_d <- function(i) {
  temp_sim_resp_ate <- results_all$simulated_responses_ate[[i]]
  
  temp_sim_means <- rowMeans(temp_sim_resp_ate) # Mean of simulations
  temp_sim_sds <- apply(temp_sim_resp_ate, 1, sd) # Standard deviation of simulations
  temp_sim_lower <- apply(temp_sim_resp_ate, 1, quantile, probs = 0.025)
  temp_sim_upper <- apply(temp_sim_resp_ate, 1, quantile, probs = 0.975)
  
  # Combine results in a data frame for easier inspection
  temp_newd2 <- results_all$newd2[[i]]
  
  results <- data.frame(
    lon = temp_newd2$lon,
    lat = temp_newd2$lat,
    Simulated_Mean = temp_sim_means,
    SE = temp_sim_sds,
    CI_Lower = temp_sim_lower,
    CI_Upper = temp_sim_upper
  )
  
  results_all$sim_means[[i]] <- temp_sim_means
  results_all$sim_sds[[i]] <- temp_sim_sds
  results_all$sim_lower[[i]] <- temp_sim_lower
  results_all$sim_upper[[i]] <- temp_sim_upper
  results_all$results[[i]] <- results
  
}

mclapply(1:15, create_ate_stats_d, mc.cores = 6)


######################################################################
# Compute posterior prob that the diff in Mortality (ATE) <0
######################################################################
create_probeff_d <- function(i) {
  temp_sims_probeff <- apply(results_all$simulated_responses_ate[[i]],
                           1, 
                           FUN = function(x) mean(ifelse(x<0,1,0),na.rm=T))
  
  results_all$sims_mort_probeff[[i]] <- temp_sims_probeff
}

mclapply(1:15, create_probeff_d, mc.cores = 6)


######################################################################
# Convert probeff surfaces to raster then to tibble for ggplot2 ####
######################################################################
create_probeff_rasters_d <- function(i) {
  temp_coords_tx1 <- results_all$b_preds_mort_coords_tx1[[i]]
  temp_sims_probeff <- results_all$sims_mort_probeff[[i]]
  
  temp_probeff_raster <- points_to_raster(
    x = temp_coords_tx1[, 1],
    y = temp_coords_tx1[, 2],
    z = temp_sims_probeff,
    mask1 = niger_adm2 %>% filter(adm_02==districts[i]),
    crop1 = niger_adm2 %>% filter(adm_02==districts[i])
  )
  
  temp_probeff_raster2 <- raster_to_tibble(temp_probeff_raster)
  
  results_all$probeff_raster[[i]] <- temp_probeff_raster
  results_all$probeff_raster2[[i]] <- temp_probeff_raster2
  
}

mclapply(1:15, create_probeff_rasters_d, mc.cores = 6)


######################################################################
# Done! Save workspace image.
######################################################################
save.image(paste0(
  "~/Documents/GitHub/AVENIR-spatial-ate/data-",
  Sys.Date(),
  ".RData"
))
