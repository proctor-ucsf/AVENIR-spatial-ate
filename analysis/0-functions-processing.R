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
                             grappes = T,
                             min = 0, max = 35) {
  
  if(!grappes) {
    
    map_mort_temp <- ggplot(data=raster2) + 
      geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
      geom_sf(data=d_sf, alpha=0, size=0, shape=1) + #grappes
      geom_sf(data=niger_adm2,fill=NA) + #boundaries
      coord_sf(crs = 4326) +
      labs(x="Longitude",y="Latitude",title=tagtxt) +
      theme_minimal() +
      theme(plot.tag = element_text(face="bold",size=16),
            legend.position = "bottom")
    
  } else if(grappes){
    
    map_mort_temp <- ggplot(data=raster2) + 
      geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
      geom_sf(data=d_sf, alpha=0.9, size=0.5, shape=1) + #grappes
      geom_sf(data=niger_adm2,fill=NA) + #boundaries
      coord_sf(crs = 4326) +
      labs(x="Longitude",y="Latitude",title=tagtxt) +
      theme_minimal() +
      theme(plot.tag = element_text(face="bold",size=16),
            legend.position = "bottom")
    
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


#------------------------------------
# create_mort_map_districts()
#------------------------------------
create_mort_map_districts <- function (raster2,
                                       agetxt      = "(1 to 59)",
                                       tagtxt      = "Overall",
                                       legtxt      = "Mortality rate (per 1000 person-years)",
                                       distr       = districts[1],
                                       paletteopt  = 1,
                                       min = -17, max = 5, interval = 1,
                                       include_axes_labs = T,
                                       include_legend = T) {
  
  
  map_mort_temp <- ggplot(data=raster2) + 
    geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
    geom_sf(data=d_sf %>% filter(district==distr), alpha=0.9, size=0.3, shape=1) + #grappes
    geom_sf(data=niger_adm2 %>% filter(adm_02==distr),fill=NA) +      #boundaries
    coord_sf(crs = 4326) +
    labs(x="Longitude",y="Latitude",title=tagtxt) +
    ggspatial::annotation_scale(location = "bl",
                                line_width = unit(0.8, "cm"),
                                height = unit(0.1, "cm")) +
    # ggspatial::annotation_north_arrow(location = "tl",,
    #                                   height = unit(0.5, "cm"),
    #                                   width = unit(0.5, "cm"),
    #                                   pad_x = unit(0.4, "cm"),
    #                                   pad_y = unit(0.8, "cm"),
    #                                   style = north_arrow_fancy_orienteering
    #                                   ) +
    theme_minimal()
  
  
  if (include_axes_labs && include_legend) {
    
    map_mort_temp <-
      map_mort_temp +
      theme(axis.text = element_text(size=5),
            axis.title = element_text(size=8),
            plot.title = element_text(size=10),
            plot.tag = element_text(face="bold",size=8),
            legend.position = "bottom",
            legend.title = element_text(size=8))
    
  } else if (!include_axes_labs && include_legend) {
    
    map_mort_temp <-
      map_mort_temp +
      theme(axis.text=element_text(size=5),
            axis.title=element_blank(),
            plot.title = element_text(size=10),
            plot.tag = element_text(face="bold",size=8),
            legend.position = "bottom",
            legend.title = element_text(size=8))
    
  } else if (include_axes_labs && !include_legend) {
    
    map_mort_temp <-
      map_mort_temp +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=8),
            plot.title = element_text(size=10),
            plot.tag = element_text(face="bold",size=8),
            legend.position = "none")
    
  } else if (!include_axes_labs && !include_legend) {
    
    map_mort_temp <-
      map_mort_temp +
      theme(axis.text=element_text(size=6),
            axis.title=element_blank(),
            plot.title = element_text(size=10),
            plot.tag = element_text(face="bold",size=8),
            legend.position = "none")
    
  }
  
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
                                              barheight=unit(5,"pt"),
                                              barwidth=unit(100,"pt"),
                                              ticks.colour = "white",
                                              ticks.linewidth = 1,
                                              frame.colour = "white")
      )
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
                                                barheight=unit(5,"pt"),
                                                barwidth=unit(100,"pt"),
                                                ticks.colour = "white",
                                                ticks.linewidth = 1,
                                                frame.colour = "white",
                                                draw.llim = T,
                                                draw.ulim = T)
      )
    
  } else if (paletteopt==2.1) {
    map_mort_temp <-
      map_mort_temp +
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         end = max,
                         labels = as.character(seq(min,max,interval)),
                         guide=ifelse(include_legend,guide_colorbar(title=legtxt,
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
                                                                    draw.ulim = T), "none" ))
    
    
  } else if (paletteopt==3) {
    map_mort_temp <-
      map_mort_temp +
      scale_fill_viridis(option = "E", direction = 1,
                         limits = c(min,max),
                         na.value = NA,
                         begin = min,
                         breaks = seq(min,max,by=interval),
                         guide=ifelse(include_legend,guide_colorbar(title=legtxt,
                                                                    title.position = "top",
                                                                    title.hjust = 0.5,
                                                                    direction="horizontal",
                                                                    label.position = "bottom",
                                                                    label.hjust = 0.5,
                                                                    barheight=unit(10,"pt"),
                                                                    barwidth=unit(100,"pt"),
                                                                    ticks.colour = "white",
                                                                    ticks.linewidth = 1,
                                                                    frame.colour = "white"), "none" ))
    
    
  }
  
  # ggsave(filename=paste0(figure_path,"/spatial-map-panels/avenir-spatial-ATE-map-",distr,".png"),
  #        plot=map_mort_temp,
  #        device="png",
  #        width=3.3,
  #        height=3.5)
  
  return(map_mort_temp)
}