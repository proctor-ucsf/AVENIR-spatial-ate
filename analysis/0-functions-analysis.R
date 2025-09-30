## 0-functions-analysis.R

#--------------------------------------
# create_fitme()
# updated 11/7/2024
#--------------------------------------
create_fitme <- function (df  = df_all,
                          age = "1to59") {
  
    b_fitmort_temp <- spaMM::fitme(n_deaths_1to59 ~ 
                                     1 + lat + lon + 
                                     Matern(1|lat + lon) +
                                     isTreated_child +
                                     (1|site_code) +
                                     offset(log(pt_year_1to59)),
                                   data = df,
                                   family = poisson(link = "log"))
  
  return(b_fitmort_temp)
  
}

create_fitme_int <- function (df  = df_all,
                          age = "1to59") {
  
  b_fitmort_temp <- spaMM::fitme(n_deaths_1to59 ~ 
                                   1 + lat + lon + 
                                   Matern(1|lat + lon) +
                                   isTreated_child +
                                   isTreated_child*lat + isTreated_child*lon +
                                   (1|site_code) +
                                   offset(log(pt_year_1to59)),
                                 data = df,
                                 family = poisson(link = "log"))
  
  return(b_fitmort_temp)
  
}

#--------------------------------------
# get_grid_preds()
# updated 11/7/2024
#
# get predictions (mean, var) from
# a spaMM::fitme() geospatial model 
# along a grid
# of new locations
#
# input_grid       : grid to get values over, class = sf (points)
# spamm_model_fit  : spaMM geospatial model fit object, class = HLfit
#--------------------------------------
get_grid_preds <- function(input_grid, 
                           spamm_model_fit, 
                           tx_group="Treatment"){
  
  # convert grid into a data.frame with lon/lat variables
  newd <- input_grid %>%
    st_centroid() %>% 
    st_as_sf(crs = 4326) %>% 
    mutate(lon = sf::st_coordinates(.)[,1],
           lat = sf::st_coordinates(.)[,2],
           isTreated_child = tx_group,
           site_code = 0, #dummy
           pt_year_1to59 = 1000) %>% #dummy
    st_drop_geometry()
  
  print(newd)
  
  # get predictions (response)
  
  temp_preds <- predict(object = spamm_model_fit,
                        type = "response",
                        newdata = newd)
  
  # get estimated variance at each location
  temp_predVar <- spaMM::get_predVar(object = spamm_model_fit, 
                                     newdata = newd, 
                                     which = "predVar")
  
  # create a dataframe 
  # for predictions using 
  # centroid coordinates from the grid
  # add predictions and variance to the grid
  ret <- input_grid %>% 
    as.data.frame() %>% 
    st_as_sf(crs = 4326) %>% 
    mutate(pred = as.numeric(temp_preds),
           pred_var = as.numeric(temp_predVar))
  
  return(ret)
  
}


#--------------------------------------
# plot_matern_corr()
#
# estimate the effect correlation by
# distance relationship given fitted
# Matern parameters in a geostatistical model
# of class HLfit from spaMM
#
# the model fit must have variables 
# "lon" and "lat" as the x,y coordinates
#--------------------------------------
plot_matern_corr <- function(spamm_model_fit, 
                             filenameval="avenir-materncorr") {
  
  # create a distance matrix between all points
  # compute great circle distances in meters
  dist_mat_m <- spamm_model_fit$data %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>%
    st_distance(which = "Great Circle") %>%
    as.numeric()
  
  model_dist_mat <- as.matrix(dist(spamm_model_fit$data [,c("lon","lat")]))
  
  # get estimated parameters from the Matern function
  est_nu <- spamm_model_fit$CorrEst_and_RanFix$corrPars$`1`$nu
  est_rho <- spamm_model_fit$CorrEst_and_RanFix$corrPars$`1`$rho
  # estimate the correlation between each location
  matern_cor <- MaternCorr(d = model_dist_mat, nu=est_nu, rho=est_rho)
  # store in a data frame for plotting
  d_matern_cor <- data.frame(
    distms = as.numeric(model_dist_mat),
    distkm = as.numeric(dist_mat_m)/1000,
    materncor = as.numeric(matern_cor)) %>%
    arrange(distms) %>%
    group_by(distkm) %>%
    slice(1) %>%
    ungroup()
  
  # there is a very, very slightly imperfect translation
  # between units in lat/lon used in the geostatistical
  # model and the actual great circle distances on earth
  # regress the two to get a very accurate estimate for plotting
  mskmfit <- lm(distkm~distms,data=d_matern_cor)
  d_matern_cor$pdistkm <- predict(mskmfit)
  
  #--------------------------------------
  # plot the correlation function
  #--------------------------------------
  plot_materncor <- ggplot(data = d_matern_cor %>% filter(pdistkm<=40)) +
    # geom_point(aes(x=distms,y=distkm)) +
    geom_point(aes(x=pdistkm, y = materncor),color="black",alpha=1,size=0.1) + 
    # geom_vline(xintercept= 3 ,color="gray60", lwd=0.1) +
    scale_y_continuous(breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(breaks=c(0,seq(10,40,by=10))) +
    coord_cartesian(xlim=c(0,40),ylim=c(0,1)) +
    labs(x="distance (km)",y = "Matern correlation") +
    annotate("text",x=15,y=0.9,label=paste("nu ==", deparse(sprintf("%1.1f",est_nu))),parse=TRUE, hjust =0, size=8)+
    annotate("text",x=15,y=0.7,label=paste("rho ==", deparse(sprintf("%1.1f",est_rho))),parse=TRUE, hjust =0, size=8)+
    theme_classic() +
    theme(
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size=16),
      axis.title = element_text(size=18)
    )
  
  ggsave(filename=paste0(figure_path,filenameval,".png"),
         plot=plot_materncor,
         device="png",
         width=3.3,
         height=3.5)
  
  return(plot_materncor)
  
}


#--------------------------------------
# get_grid_uncert()
# 
# get uncertainty from bootstrapping
#--------------------------------------
get_grid_uncert <- function(input_grid, 
                            spamm_model_fit, 
                            nsim_input=10,
                            tx_group="Treatment"){
  
  # convert grid into a data.frame with lon/lat variables
  newd <- input_grid %>%
    st_centroid() %>% 
    st_as_sf(crs = 4326) %>% 
    mutate(lon = sf::st_coordinates(.)[,1],
           lat = sf::st_coordinates(.)[,2],
           isTreated_child = tx_group,
           site_code = 0, #dummy
           pt_year_1to59 = 1000) %>% #dummy
    st_drop_geometry()
  
  # get bootstrap simulations
  temp_sims <- simulate(object = spamm_model_fit,
                        nsim = nsim_input,
                        newdata = newd)
  
  # get estimated variance at each location
  temp_predVar <- spaMM::get_predVar(object = spamm_model_fit, 
                                     newdata = newd, 
                                     which = "predVar")
  
  # create a dataframe 
  # for predictions using 
  # centroid coordinates from the grid
  # add predictions and variance to the grid
  ret <- input_grid %>% 
    as.data.frame() %>% 
    st_as_sf(crs = 4326) %>% 
    mutate(pred = as.numeric(temp_preds),
           pred_var = as.numeric(temp_predVar))
  
  return(ret)
  
}

