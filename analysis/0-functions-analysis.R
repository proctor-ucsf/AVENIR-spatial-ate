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
#--------------------------------------
points_to_raster <- function(x,y,z,mask1=NULL,mask2=NULL,crop1=NULL) {
  
  # make raster from points
  pred_raster <- rasterFromXYZ(xyz=data.frame(x = x, y = y, z = z),crs = 4326)
  
  # mask and crop raster if specified
  if(!is.null(mask1)) {
    pred_raster <- mask(x = pred_raster,mask = mask1)
  }
  if(!is.null(mask2)) {
    pred_raster <- mask(x = pred_raster,mask = mask2)
  }
  if(!is.null(crop1)) {
    pred_raster <- crop(x = pred_raster,y = extent(crop1))
  }
  
  return(pred_raster)
}

#------------------------------------
# raster_to_tibble()
#------------------------------------
#' Transform raster as data.frame to be later used with ggplot
#' Modified from rasterVis::gplot
#' From: https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r
#'
#' @param x A Raster* object
#' @param maxpixels Maximum number of pixels to use
#'
#' @details rasterVis::gplot is nice to plot a raster in a ggplot but
#' if you want to plot different rasters on the same plot, you are stuck.
#' If you want to add other information or transform your raster as a
#' category raster, you can not do it. With `SDMSelect::gplot_data`, you retrieve your
#' raster as a tibble that can be modified as wanted using `dplyr` and
#' then plot in `ggplot` using `geom_tile`.
#' If Raster has levels, they will be joined to the final tibble.
#'
#' @export

raster_to_tibble <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- tibble::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  return(dat)
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





#------------------------------------
# create_mort_map()
#------------------------------------
create_mort_map <- function (raster2     = b_preds_mort_raster2_s0,
                             agetxt      = "(1 to 59)",
                             tagtxt      = "Overall",
                             legtxt      = "Mortality rate (per 1000 person-years)",
                             paletteopt  = 1,
                             min = 0, max = 35) {
  
  map_mort_temp <- ggplot(data=raster2) + 
    geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
    geom_sf(data=d_sf, alpha=0.9, size=0.5, shape=1) + #grappes
    geom_sf(data=niger_adm2,fill=NA) +      #boundaries
    coord_sf(crs = 4326) +
    labs(x="Longitude",y="Latitude",title=tagtxt) +
    theme_minimal() +
    theme(plot.tag = element_text(face="bold",size=16),
          legend.position = "bottom")
  
  if(paletteopt==1){ 
    map_mort_temp <- 
      map_mort_temp + 
      scale_fill_viridis(option = "D", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         breaks = seq(min,max,5),
                         guide=guide_colorbar(title=legtxt, 
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white") )
  } else if (paletteopt==2) {
    map_mort_temp <- 
      map_mort_temp + 
      scale_fill_gradient2(name=legtxt,
                           midpoint=0,
                           low = muted("dodgerblue"),
                           mid = "white",
                           high = muted("red"),
                           limits = c(min,max),
                           breaks = seq(min,max,0.5),
                           labels = as.character(seq(min,max,0.5)),
                           na.value = NA,
                           guide=guide_colorbar(title=legtxt, 
                                                title.position = "top",
                                                title.hjust = 0.5,
                                                direction="horizontal",
                                                label.position = "bottom",
                                                label.hjust = 0.5,
                                                barheight=unit(10,"pt"),
                                                barwidth=unit(200,"pt"),
                                                ticks.colour = "white",
                                                ticks.linewidth = 1,
                                                frame.colour = "white",
                                                draw.llim = T,
                                                draw.ulim = T))
    
  } else if (paletteopt==2.1) {
    map_mort_temp <- 
      map_mort_temp + 
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         end = max,
                         labels = as.character(seq(min,max,0.5)),
                         guide=guide_colorbar(title=legtxt, 
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white",
                                              draw.llim = T,
                                              draw.ulim = T) )
    
    
  } else if (paletteopt==3) {
    map_mort_temp <- 
      map_mort_temp + 
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         begin = min,
                         # breaks = seq(0,10,by=2),
                         guide=guide_colorbar(title=legtxt, 
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white") )
    
    
  }
  
  ggsave(filename=paste0(figure_path,"/spatial-map-panels/avenir-spatial-ATE-map-",substr(tagtxt,1,1),".png"),
         plot=map_mort_temp,
         device="png",
         width=3.3,
         height=3.5)
  
  return(map_mort_temp)
}

create_mort_map_districts <- function (raster2,
                                       agetxt      = "(1 to 59)",
                                       tagtxt      = "Overall",
                                       legtxt      = "Mortality rate (per 1000 person-years)",
                                       distr       = districts[1],
                                       paletteopt  = 1,
                                       min = 0, max = 35, interval = 1) {
  
  map_mort_temp <- ggplot(data=raster2) + 
    geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
    geom_sf(data=d_sf %>% filter(district==distr), alpha=0.9, size=0.5, shape=1) + #grappes
    geom_sf(data=niger_adm2 %>% filter(adm_02==distr),fill=NA) +      #boundaries
    coord_sf(crs = 4326) +
    labs(x="Longitude",y="Latitude",title=tagtxt) +
    theme_minimal() +
    theme(plot.tag = element_text(face="bold",size=16),
          legend.position = "bottom")
  
  if(paletteopt==1){
    map_mort_temp <-
      map_mort_temp +
      scale_fill_viridis(option = "D", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         breaks = seq(min,max,interval),
                         guide=guide_colorbar(title=legtxt,
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white") )
  } else if (paletteopt==2) {
    map_mort_temp <-
      map_mort_temp +
      scale_fill_gradient2(name=legtxt,
                           midpoint=0,
                           low = muted("dodgerblue"),
                           mid = "white",
                           high = muted("red"),
                           limits = c(min,max),
                           breaks = seq(min,max,interval),
                           labels = as.character(seq(min,max,interval)),
                           na.value = NA,
                           guide=guide_colorbar(title=legtxt,
                                                title.position = "top",
                                                title.hjust = 0.5,
                                                direction="horizontal",
                                                label.position = "bottom",
                                                label.hjust = 0.5,
                                                barheight=unit(10,"pt"),
                                                barwidth=unit(200,"pt"),
                                                ticks.colour = "white",
                                                ticks.linewidth = 1,
                                                frame.colour = "white",
                                                draw.llim = T,
                                                draw.ulim = T))

  } else if (paletteopt==2.1) {
    map_mort_temp <-
      map_mort_temp +
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         end = max,
                         labels = as.character(seq(min,max,interval)),
                         guide=guide_colorbar(title=legtxt,
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white",
                                              draw.llim = T,
                                              draw.ulim = T) )


  } else if (paletteopt==3) {
    map_mort_temp <-
      map_mort_temp +
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         begin = min,
                         breaks = seq(min,max,by=interval),
                         guide=guide_colorbar(title=legtxt,
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              direction="horizontal",
                                              label.position = "bottom",
                                              label.hjust = 0.5,
                                              barheight=unit(10,"pt"),
                                              barwidth=unit(200,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white") )


  }
  
  ggsave(filename=paste0(figure_path,"/spatial-map-panels/avenir-spatial-ATE-map-",substr(tagtxt,1,1),".png"),
         plot=map_mort_temp,
         device="png",
         width=3.3,
         height=3.5)
  
  return(map_mort_temp)
}
