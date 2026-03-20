## 0-functions-analysis.R

#--------------------------------------
# plot_matern_corr_proj()
#
# estimate the effect correlation by
# distance relationship given fitted
# Matern parameters in a geostatistical model
# of class HLfit from spaMM
#
# the model fit must have variables
# "lon" and "lat" as the x,y coordinates
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

  # store in a data frame for plotting
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

  label_text <- sprintf(
    "nu = %0.3f\rho = %0.3f\nrange = %0.2f km",
    est_nu,
    est_rho,
    range_km
  )

  # plot the correlation function
  plot_materncor <- ggplot(data = d_matern_cor %>% filter(distkm <= 50)) +
    # geom_point(aes(x=distms,y=distkm)) +
    geom_point(
      aes(x = distkm, y = materncor),
      color = "black",
      alpha = 1,
      size = 0.1
    ) +
    # geom_vline(xintercept= 3 ,color="gray60", lwd=0.1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = c(0, seq(10, 50, by = 10))) +
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
# updated 3/10/2026
#
# get predictions (mean, var) from
# a spaMM::fitme() geospatial model
# along a grid
# of new locations
#
# input_grid       : grid to get values over, class = sf (points)
# spamm_model_fit  : spaMM geospatial model fit object, class = HLfit
#--------------------------------------

get_grid_preds <- function(
  input_grid,
  spamm_model_fit,
  tx_group = "Treatment"
) {
  # convert grid into a data.frame with lon/lat variables
  newd <- input_grid %>%
    st_centroid() %>%
    # st_transform(crs = proj_crs) %>%
    mutate(
      lon = sf::st_coordinates(.)[, 1],
      lat = sf::st_coordinates(.)[, 2],
      isTreated_child = as.factor(tx_group),
      site_code = "0", #dummy
      pt_year_1to59 = 1000
    ) %>% #per 1000 person-years
    st_drop_geometry()

  # str(newd)

  # get predictions (response)

  temp_preds <- predict(
    object = spamm_model_fit,
    type = "response",
    newdata = newd
  )

  # get estimated variance at each location
  # temp_predVar <- spaMM::get_predVar(object = spamm_model_fit,
  #                                    newdata = newd,
  #                                    which = "predVar")

  # create a dataframe
  # for predictions using
  # centroid coordinates from the grid
  # add predictions and variance to the grid
  ret <- input_grid %>%
    as.data.frame() %>%
    st_as_sf(crs = proj_crs) %>%
    mutate(
      pred = as.numeric(temp_preds),
      # pred_var = as.numeric(temp_predVar),
      isTreated_child = tx_group
    )

  return(ret)
}

# Wrapper for get_grid_preds to get predictions for both treatment groups and combine into one dataframe

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

  preds_long <- rbind(preds_tx0, preds_tx1)

  return(preds_long)
}

# Geospatial helper functions from Ben Arnold's repo - https://osf.io/cxb5e
# Geographic pair matching in large-scale cluster randomized trials
# Functions from file https://osf.io/eqcvh

# points_to_raster()
#
# convert predicted surface of points
# to a raster to make it a regular
# set of tiles (rather than points),
# this function takes as input
# an object of class sf dataframe
# in particular, the results of the function
# above, get_grid_preds()
#
# @x,y,z  : longitude, latitude, and some value to plot in raster
# @mask1  : mask the raster by a variable? if NULL, skips.
# @mask2  : mask the raster by a second variable? if NULL skips.
# @crop1  : crop the raster by the first variable?  if NULL, skips

points_to_raster <- function(
  x,
  y,
  z,
  mask1 = NULL,
  mask2 = NULL,
  crop1 = NULL,
  crs_inp = global_crs
) {
  # make raster from points
  pred_raster <- raster::rasterFromXYZ(
    xyz = data.frame(x = x, y = y, z = z),
    crs = crs_inp
  )

  # mask and crop raster if specified
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

make_raster_map <- function(
  df,
  proj_crs = proj_crs,
  global_crs = global_crs,
  colname = "pred"
) {
  # extract coordinates
  df2 <- df %>%
    mutate(
      x = sf::st_coordinates(.)[, 1],
      y = sf::st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry()

  # rasterize
  r_proj <- terra::rast(df2[, c("x", "y", colname)], type = "xyz")
  terra::crs(r_proj) <- sf::st_crs(proj_crs)$wkt

  # project
  r_global <- terra::project(
    r_proj,
    sf::st_crs(global_crs)$wkt,
    method = "bilinear"
  )

  return(r_global)
}
