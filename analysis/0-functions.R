#--------------------------------------------
#
# AVENIR-spatial-ate-code-check/analysis/0-functions.R
#
# Custom helper functions for the AVENIR spatial ATE analysis.
#
#--------------------------------------------

#--------------------------------------
# run_boot()
#
# Stratified bootstrap scaffold. Resamples data within each district stratum
# with replacement, applies a caller-supplied summary function to each bootstrap
# dataset, and returns the list of per-replicate results. Execution is
# parallelised with parallel::mclapply().
#
# Arguments:
# - data:      data.frame; the grappe-level dataset to resample. Must contain
#              a "district" column used to define strata.
# - districts: character vector; unique district names that define the strata
# - n_boot:    integer; number of bootstrap replicates
# - mc_cores:  integer; number of cores passed to parallel::mclapply()
# - iter_fn:   function(boot_data) -> any; applied to each resampled dataset.
#              Must return a consistently-shaped object across all replicates.
#
# Details:
# - Resampling is stratified: within each district, rows are sampled with
#   replacement at the same size as the observed stratum, preserving
#   district-level sample sizes across replicates.
# - iter_fn receives the full resampled data.frame, not row indices.
#
# Requirements:
# - Packages: parallel, purrr, dplyr
#
# Returns:
# - A list of length n_boot; each element is the return value of iter_fn()
#--------------------------------------

run_boot <- function(data, districts, n_boot, mc_cores, iter_fn) {
  parallel::mclapply(seq_len(n_boot), function(i) {
    boot_data <- purrr::map_dfr(districts, function(d) {
      d_dat <- dplyr::filter(data, district == d)
      d_dat[sample(nrow(d_dat), size = nrow(d_dat), replace = TRUE), ]
    })
    iter_fn(boot_data)
  }, mc.cores = mc_cores)
}


#--------------------------------------
# boot_mort_rates()
#
# Bootstrap district-level and pooled mortality rates (deaths per 1000
# person-years), stratified by a caller-specified grouping variable. A thin
# wrapper around run_boot() that fixes the mortality-rate summary function
# used throughout the AVENIR analysis.
#
# Arguments:
# - data:      data.frame; grappe-level data containing n_deaths_1to59,
#              pt_year_1to59, district, and the column named by group_var
# - districts: character vector; unique district names (strata for resampling)
# - n_boot:    integer; number of bootstrap replicates
# - mc_cores:  integer; number of cores passed to run_boot()
# - group_var: character; name of the grouping column — either
#              "isTreated_child" (treatment vs. control) or "tx"
#              (three-arm randomization label)
#
# Requirements:
# - Depends on: run_boot()
# - Packages: dplyr, tidyr, purrr
#
# Returns:
# - A list of length n_boot; each element is a data.frame with columns
#   district, MortRate, and <group_var>, with one row per district x
#   group combination plus an "Overall" pooled row
#--------------------------------------

boot_mort_rates <- function(data, districts, n_boot, mc_cores, group_var) {
  run_boot(data, districts, n_boot, mc_cores, function(boot_data) {
    district_rows <- boot_data %>%
      group_by(district, .data[[group_var]]) %>%
      summarise(
        n_deaths_1to59 = sum(n_deaths_1to59, na.rm = TRUE),
        pt_year_1to59  = sum(pt_year_1to59,  na.rm = TRUE),
        MortRate       = n_deaths_1to59 * 1e3 / pt_year_1to59,
        .groups        = "drop"
      ) %>%
      dplyr::select(district, MortRate, all_of(group_var))
    
    overall_row <- boot_data %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        n_deaths_1to59 = sum(n_deaths_1to59, na.rm = TRUE),
        pt_year_1to59  = sum(pt_year_1to59,  na.rm = TRUE),
        MortRate       = n_deaths_1to59 * 1e3 / pt_year_1to59,
        .groups        = "drop"
      ) %>%
      mutate(district = "Overall") %>%
      dplyr::select(district, MortRate, all_of(group_var))
    
    bind_rows(district_rows, overall_row)
  })
}


#--------------------------------------
# boot_ate()
#
# Bootstrap district-level and pooled average treatment effects (ATE) as the
# difference in pooled mortality rates (Treatment - Control, deaths per 1000
# person-years). Stratified resampling is performed within districts; the
# treatment indicator is fixed to isTreated_child.
#
# Arguments:
# - data:      data.frame; grappe-level data containing n_deaths_1to59,
#              pt_year_1to59, district, and isTreated_child
# - districts: character vector; unique district names (strata for resampling)
# - n_boot:    integer; number of bootstrap replicates
# - mc_cores:  integer; number of cores passed to run_boot()
#
# Requirements:
# - Depends on: run_boot()
# - Packages: dplyr, tidyr
#
# Returns:
# - A list of length n_boot; each element is a data.frame with columns
#   district and ATE, with one row per district plus an "Overall" pooled row
#--------------------------------------

boot_ate <- function(data, districts, n_boot, mc_cores) {
  run_boot(data, districts, n_boot, mc_cores, function(boot_data) {
    district_rows <- boot_data %>%
      group_by(district, isTreated_child) %>%
      summarise(
        n_deaths_1to59 = sum(n_deaths_1to59, na.rm = TRUE),
        pt_year_1to59  = sum(pt_year_1to59,  na.rm = TRUE),
        MortRate       = n_deaths_1to59 * 1e3 / pt_year_1to59,
        .groups        = "drop"
      ) %>%
      pivot_wider(
        names_from  = isTreated_child,
        values_from = c(n_deaths_1to59, pt_year_1to59, MortRate)
      ) %>%
      mutate(ATE = MortRate_Treatment - MortRate_Control) %>%
      dplyr::select(district, ATE)
    
    overall_row <- boot_data %>%
      group_by(isTreated_child) %>%
      summarise(
        n_deaths_1to59 = sum(n_deaths_1to59, na.rm = TRUE),
        pt_year_1to59  = sum(pt_year_1to59,  na.rm = TRUE),
        MortRate       = n_deaths_1to59 * 1e3 / pt_year_1to59,
        .groups        = "drop"
      ) %>%
      pivot_wider(
        names_from  = isTreated_child,
        values_from = c(n_deaths_1to59, pt_year_1to59, MortRate)
      ) %>%
      mutate(district = "Overall", ATE = MortRate_Treatment - MortRate_Control) %>%
      dplyr::select(district, ATE)
    
    bind_rows(district_rows, overall_row)
  })
}


#--------------------------------------
# plot_matern_corr_proj()
#
# Plot and save the Matérn spatial correlation as a function of distance,
# using the estimated nu and rho parameters from a fitted geostatistical
# model. Bins pairwise distances at 25 m resolution before plotting to
# reduce overplotting; the plot is also saved to disk as a .png.
#
# Arguments:
# - spamm_model_fit: fitted spaMM HLfit object; must contain "lon" and
#                    "lat" columns in $data (UTM-projected metres)
# - filenameval:     output filename stem (no extension); saved under the
#                    global `figure_path` directory as a .png
# - modelname:       character label printed in the top-left of the plot
#
# Requirements:
# - Depends on global: figure_path (defined in 0-config.R)
# - Packages: ggplot2, spaMM, dplyr, tibble
#
# Returns:
# - A ggplot object of the Matérn correlation vs. distance (invisible side
#   effect: .png saved to figure_path)
#--------------------------------------

plot_matern_corr_proj <- function(
  spamm_model_fit,
  filenameval = "avenir-materncorr",
  modelname = ""
) {
  # create a distance matrix between all points
  coords <- as.matrix(spamm_model_fit$data[, c("lon", "lat")])
  d_vec_m <- as.numeric(dist(coords))

  # get estimated parameters from the Matern function
  est_nu <- spamm_model_fit$CorrEst_and_RanFix$corrPars$`1`$nu
  est_rho <- spamm_model_fit$CorrEst_and_RanFix$corrPars$`1`$rho

  rho_sci <- format(est_rho, scientific = TRUE, digits = 3)
  range_km <- 1 / (1000 * est_rho)

  # estimate the correlation between each location
  matern_vec <- spaMM::MaternCorr(d = d_vec_m, nu = est_nu, rho = est_rho)

  # NOTE: label_text is unused — annotations below use bquote() for math formatting instead
  # label_text <- sprintf("nu = %0.3f\nrho = %0.3f\nrange = %0.2f km", est_nu, est_rho, range_km)

  # bin pairwise distances at 25 m resolution to reduce overplotting
  d_matern_cor <- tibble(
    distm = d_vec_m,
    distkm = distm / 1000,
    materncor = matern_vec
  ) %>%
    mutate(
      bin_km = cut(
        distkm,
        breaks = seq(0, max(distkm), by = 0.025),
        include.lowest = TRUE
      )
    ) %>%
    group_by(bin_km) %>%
    summarise(
      distkm = mean(distkm),
      materncor = mean(materncor),
      .groups = "drop"
    ) %>%
    arrange(distkm)

  # plot the correlation function
  plot_materncor <- ggplot(data = d_matern_cor %>% dplyr::filter(distkm <= 50)) +
    geom_point(
      aes(x = distkm, y = materncor),
      color = "black",
      size = 0.1
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(0, 50, by = 10)) +
    coord_cartesian(xlim = c(0, 50), ylim = c(0, 1)) +
    labs(x = "distance (km)", y = "Matern correlation") +
    annotate(
      "text",
      x = 15,
      y = 1,
      label = modelname,
      hjust = 0,
      size = 6
    ) +
    annotate(
      "text",
      x = 15,
      y = 0.90,
      label = as.expression(bquote(
        nu ==
          .(
            format(est_nu, digits = 3)
          )
      )),
      hjust = 0,
      size = 6
    ) +
    annotate(
      "text",
      x = 15,
      y = 0.80,
      label = as.expression(bquote(
        rho == .(rho_sci) ~ m^{
          -1
        }
      )),
      hjust = 0,
      size = 6
    ) +
    annotate(
      "text",
      x = 15,
      y = 0.65,
      label = as.expression(bquote(
        range ==
          .(
            format(range_km, digits = 4)
          ) ~ km
      )),
      hjust = 0,
      size = 6
    ) +
    theme_classic() +
    theme(
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18)
    )

  ggsave(
    filename = paste0(figure_path, filenameval, ".png"),
    plot = plot_materncor,
    device = "png",
    width = 7,
    height = 5
  )

  return(plot_materncor)
}


#--------------------------------------
# get_grid_preds()
#
# Get predicted mean mortality rate from a spaMM::fitme() geostatistical
# model at new locations defined by a spatial grid. Predictions are on the
# response scale (deaths per 1000 person-years) using a fixed offset of 1000.
# NOTE: variance prediction (predVar) is available but currently commented
#       out — see the get_predVar() call in the function body.
#
# Arguments:
# - input_grid:      sf points object; grid of prediction locations
# - spamm_model_fit: spaMM HLfit object; fitted geostatistical model
# - tx_group:        character; treatment arm label assigned to all grid cells,
#                    one of "Treatment" (default) or "Control". site_code is
#                    set to "0" (dummy) and pt_year_1to59 to 1000 for prediction.
#
# Requirements:
# - Depends on global: proj_crs (CRS for output sf object, defined in 0-config.R)
# - Packages: sf, dplyr, spaMM
#
# Returns:
# - An sf points dataframe with the same geometry as input_grid, plus columns
#   `pred` (numeric, predicted mortality rate) and `isTreated_child` (tx_group)
#--------------------------------------

get_grid_preds <- function(
  input_grid,
  spamm_model_fit,
  tx_group = "Treatment"
) {
  # site_code is a dummy value (not in model); offset is per 1000 person-years
  newd <- input_grid %>%
    st_centroid() %>%
    mutate(
      lon = sf::st_coordinates(.)[, 1],
      lat = sf::st_coordinates(.)[, 2],
      isTreated_child = as.factor(tx_group),
      site_code = "0",
      pt_year_1to59 = 1000
    ) %>%
    st_drop_geometry()

  temp_preds <- predict(
    object = spamm_model_fit,
    type = "response",
    newdata = newd
  )

  # get estimated variance at each location
  # temp_predVar <- spaMM::get_predVar(object = spamm_model_fit,
  #                                    newdata = newd,
  #                                    which = "predVar")

  # attach predictions to the original grid; mutate() preserves sf geometry and CRS
  ret <- input_grid %>%
    mutate(
      pred = as.numeric(temp_preds),
      # pred_var = as.numeric(temp_predVar),
      isTreated_child = tx_group
    )

  return(ret)
}

#--------------------------------------
# get_grid_preds_both_tx()
#
# Wrapper around get_grid_preds() that runs predictions for both treatment
# arms ("Control" and "Treatment") and returns a single long-format sf
# dataframe with an isTreated_child column distinguishing the two groups.
#
# Arguments:
# - input_grid:      sf points object; grid of prediction locations
# - spamm_model_fit: spaMM HLfit object; fitted geostatistical model
#
# Requirements:
# - Depends on: get_grid_preds()
#
# Returns:
# - Row-bound sf dataframe with 2 * nrow(input_grid) rows; one set per arm
#--------------------------------------

get_grid_preds_both_tx <- function(input_grid, spamm_model_fit) {
  preds_tx0 <- get_grid_preds(
    input_grid = input_grid,
    spamm_model_fit = spamm_model_fit,
    tx_group = "Control"
  )

  preds_tx1 <- get_grid_preds(
    input_grid = input_grid,
    spamm_model_fit = spamm_model_fit,
    tx_group = "Treatment"
  )

  preds_long <- bind_rows(preds_tx0, preds_tx1)

  return(preds_long)
}

# Geospatial helper functions from Ben Arnold's repo - https://osf.io/cxb5e
# Geographic pair matching in large-scale cluster randomized trials
# Functions from file https://osf.io/eqcvh

#--------------------------------------
# points_to_raster()
#
# Convert a predicted surface of point observations into a regular raster
# grid. Adapted from Ben Arnold's geographic pair-matching repo (osf.io/cxb5e).
# NOTE: uses the legacy `raster` package; see make_raster_map() for the
#       terra equivalent.
#
# Arguments:
# - x, y, z: numeric vectors; longitude, latitude, and the value to rasterize
# - mask1:   optional sf/Spatial object; mask the raster (if NULL, skipped)
# - mask2:   optional; apply a second mask (if NULL, skipped)
# - crop1:   optional; crop the raster to this object's extent (if NULL, skipped)
# - crs_inp: CRS for the output raster; defaults to global `global_crs`
#
# Requirements:
# - Depends on global: global_crs (defined in 0-config.R)
# - Packages: raster
#
# Returns:
# - A raster object masked and cropped as specified
#--------------------------------------

points_to_raster <- function(
  x,
  y,
  z,
  mask1 = NULL,
  mask2 = NULL,
  crop1 = NULL,
  crs_inp = global_crs
) {
  pred_raster <- raster::rasterFromXYZ(
    xyz = data.frame(x = x, y = y, z = z),
    crs = crs_inp
  )

  if (!is.null(mask1)) {
    pred_raster <- mask(x = pred_raster, mask = mask1)
  }
  if (!is.null(mask2)) {
    pred_raster <- mask(x = pred_raster, mask = mask2)
  }
  if (!is.null(crop1)) {
    pred_raster <- crop(x = pred_raster, y = extent(crop1))
  }

  return(pred_raster)
}

#--------------------------------------
# make_raster_map()
#
# Convert an sf points dataframe into a terra SpatRaster projected to WGS84.
# Terra-based alternative to points_to_raster() for use with get_grid_preds()
# output.
#
# Arguments:
# - df:         sf points dataframe containing the column named by `colname`
# - proj_crs:   integer or CRS; CRS of the input data (default: global `proj_crs`).
#               NOTE: parameter name shadows the global — intentional.
# - global_crs: integer or CRS; target CRS for output raster (default: global
#               `global_crs`). NOTE: parameter name shadows the global — intentional.
# - colname:    character; name of the column in df to rasterize (default: "pred")
#
# Requirements:
# - Packages: terra, sf, dplyr
#
# Returns:
# - A terra SpatRaster reprojected to global_crs
#--------------------------------------

make_raster_map <- function(
  df,
  proj_crs = proj_crs,
  global_crs = global_crs,
  colname = "pred"
) {
  df2 <- df %>%
    mutate(
      x = sf::st_coordinates(.)[, 1],
      y = sf::st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry()

  r_proj <- terra::rast(df2[, c("x", "y", colname)], type = "xyz")
  terra::crs(r_proj) <- sf::st_crs(proj_crs)$wkt

  r_global <- terra::project(
    r_proj,
    sf::st_crs(global_crs)$wkt,
    method = "bilinear"
  )

  return(r_global)
}

#--------------------------------------
# sim2tib()
#
# Convert the matrix output of spaMM::simulate() into a tidy long-format
# tibble with one row per grid location × treatment arm × simulation draw.
# Used in 02B to post-process posterior predictive simulations before
# computing the spatial ATE and its posterior probability.
#
# Arguments:
# - sim_responses: matrix; output of spaMM::simulate() with rows = newdat
#                  rows and columns = simulation draws
# - newdat:        data.frame; the prediction dataset passed to simulate()
#                  (must be in the same row order as sim_responses)
# - modelname:     character; label added to the output for faceting/grouping
#                  (default: NA_character_)
#
# Requirements:
# - Packages: dplyr, tidyr, tibble
#
# Returns:
# - A tibble with columns: lon, lat, isTreated_child, sim_n (integer),
#   pred (simulated count), model (modelname)
#--------------------------------------

sim2tib <- function(sim_responses, newdat, modelname = NA_character_) {
  sim_tibs <- as_tibble(
    sim_responses,
    .name_repair = ~ paste0("sim_", seq_along(.))
  )

  combined <- bind_cols(
    newdat %>% dplyr::select(lon, lat, isTreated_child),
    sim_tibs
  )

  long_preds <- combined %>%
    pivot_longer(
      cols = starts_with("sim_"),
      names_to = "sim",
      values_to = "pred"
    ) %>%
    mutate(
      sim_n = as.integer(sub("^sim_", "", sim))
    ) %>%
    dplyr::select(lon, lat, isTreated_child, sim_n, pred) %>%
    mutate(model = modelname)

  return(long_preds)
}
