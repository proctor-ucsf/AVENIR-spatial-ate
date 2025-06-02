i=1

print(paste(districts[i],Sys.time()))

simulated_responses_tx0_districts <- list()
simulated_responses_tx1_districts <- list()
simulated_responses_ate_districts <- list()

simulated_responses_tx0_districts[[i]] <- simulate(results_spaMM_districts$fitmort[[i]],
                                                   type = "predVar",
                                                   variances = list(linPred=TRUE, disp=TRUE),
                                                   nsim = nsim_input,
                                                   newdata = results_spaMM_districts$newd3[[i]],
                                                   seed = 12345
)

simulated_responses_tx1_districts[[i]] <- simulate(results_spaMM_districts$fitmort[[i]],
                                                   type = "predVar",
                                                   variances = list(linPred=TRUE, disp=TRUE),
                                                   nsim = nsim_input,
                                                   newdata = results_spaMM_districts$newd2[[i]],
                                                   seed = 12345
)

simulated_responses_ate_districts[[i]] <- 
  simulated_responses_tx1_districts[[i]] - 
  simulated_responses_tx0_districts[[i]]



print(paste("Completed",Sys.time()))
  

simulated_responses_ate_districts[[1]]
simulated_responses_tx1_districts[[1]]
simulated_responses_tx0_districts[[1]]



i=1

#1 and 14 have ATE > 0; this tracks based on the spaMM estimates

d_sf2 <- st_transform(d_sf,
                      crs    = 32631) #UTM ZONE 31N (specific to that region, includes southern Niger -- for sure Dosso, checking Tahoua)

colnames(d_sf2)
# floor(min(results_spaMM_districts$b_preds_mort_raster2_tx0[[i]]$value, na.rm=T))
# ceiling(max(results_spaMM_districts$b_preds_mort_raster2_tx0[[i]]$value, na.rm=T))
# 
create_mort_map_districts(raster2 = results_spaMM_districts$b_preds_mort_raster2_diff[[i]],
                agetxt  = "(1 to 59)",
                tagtxt  = "C. ATE",
                min = floor(min(results_spaMM_districts$b_preds_mort_raster2_diff[[i]]$value, na.rm=T)),
                max = ceiling(max(results_spaMM_districts$b_preds_mort_raster2_diff[[i]]$value, na.rm=T)),
                interval = 2,
                paletteopt = 2,
                distr = districts[i],
                legtxt = "ATE: difference in mortality per 1000 person-years")


results_spaMM_districts$b_preds_mort_tx0[[2]]

ggplot(data=results_spaMM_districts$b_preds_mort_raster2_tx0[[2]]) + 
  geom_tile(aes(x=x,y=y,fill=value)) +    #gradient
  geom_sf(data=d_sf %>% filter(district==districts[2]), alpha=0.9, size=0.5, shape=1) + #grappes
  geom_sf(data=niger_adm2 %>% filter(adm_02==districts[2]),fill=NA) +      #boundaries
  coord_sf(crs = 4326) +
  labs(x="Longitude",y="Latitude",title="C. ATE") +
  theme_minimal() +
  theme(plot.tag = element_text(face="bold",size=16),
        legend.position = "bottom")
