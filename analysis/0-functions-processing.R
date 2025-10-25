## 0-functions-processing.R

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
    geom_sf(data=niger_adm2 %>% filter(adm_01=="Dosso"),fill=NA) +      #boundaries
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
                           breaks = seq(min,max,1),
                           labels = as.character(seq(min,max,1)),
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


create_mort_map_2 <- function (raster2     = b_preds_mort_raster2_s0,
                             agetxt      = "(1 to 59)",
                             tagtxt      = "Overall",
                             legtxt      = "Mortality rate (per 1000 person-years)",
                             paletteopt  = 1,
                             min = 0, max = 35,
                             add_grappes = T) {
  
  map_mort_temp <- ggplot(data=raster2) + 
    geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
    geom_sf(data=niger_adm2 %>% filter(adm_01=="Dosso"),fill=NA,color="white") +      #boundaries
    coord_sf(crs = 4326) +
    labs(x="Longitude",y="Latitude",title=tagtxt) +
    theme_minimal() +
    theme(plot.tag = element_text(face="bold",size=16),
          legend.position = "bottom")
  
  
  if(add_grappes==T) {
    map_mort_temp <- 
      map_mort_temp +
      geom_sf(data=d_sf, alpha=0.9, size=0.5, shape=1) #grappes
  }
  
  
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
                           breaks = seq(min,max,1),
                           labels = as.character(seq(min,max,1)),
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
      scale_fill_viridis(option = "B", direction = 1,
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
