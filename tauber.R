# Project: Pollen production and exposure in NYC
# Class: PLSCI 4450/6450 Urban Plants and Public Health
# Contact: Dan Katz, Cornell University
# 
# This script is for calculating dispersal kernels for pollen by combining pollen production estimates
# with Tauber samples of airborne pollen from 2013 from Kate Weinberger


#set up work environment
library(tidycensus)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(purrr)
library(cowplot)
library(readr)

focal_genus_list <- c("Acer", "Betula", "Gleditsia", "Morus", "Platanus", "Quercus", "Ulmus", "Populus", "Juglans")


### load in pollen sample data ######################################
# # #original version sent by Kate
# taub_raw <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/2013_NYCPSdataentry_fordankatz_data_tab.csv")
# 
# #use the same samples in analysis as Kate did
# taub_f_orig <- taub_raw %>%
#   #use the same samples in analysis as Kate did
#   #   filter(Main_analysis == 1) %>%
# 
#   # take an average of the replicates when available
#   group_by(Longitude, Latitude) %>%
#   summarize(Influx_trees = mean(Influx_trees),
#             Influx_plat = mean(Influx_plat),
#             Influx_que	= mean(Influx_que),
#             Influx_acer = mean(Influx_acer),
#             Influx_bet = mean(Influx_bet)) %>%
#   filter(!is.na(Influx_trees)) %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# new version including data from 2013 and 2014 sent by Guy
taub_locations <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/Kate sampling locations.csv") %>% 
  rename(nyccas_id = SITEID)

taub2013 <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/Kate_Guy_tauber_2013_400grains_rawcounts.csv")%>% 
  mutate(Influx_acer = ((acer_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_bet = ((bet_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_gle = ((gledit_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_mor = ((morus_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_plat = ((platanus_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_que = ((quercus_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_ulm = ((ulmus_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_pop = ((pop_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_jug = ((jug_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_trees = ((total_tree_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_grass = ((poa_count_2013 * mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013,
         Influx_herb =  (((amb_count_2013 + art_count_2013 + otherasther_count_2013 + chenoam_count_2013 +
                             gallium_count_2013 + plantago_count_2013 + poa_count_2013 + rumex_count_2013 +
                             taraxacum_count_2013 + typha_m_count_2013 + typha_t_count_2013 + urtica_count_2013) * 
                            mean_lycopodium_tablet_2013)/lycopodium_count_2013)/tauber_area_2013) %>% 
    left_join(., taub_locations) %>%  
    select(nyccas_id, Longitude, Latitude, year, contains("Influx_")) %>% 
    filter(!nyccas_id %in% c("Calder", "Law School")) %>%  #law school is at a different height so shouldn't be included 
    #take an average of the replicates when available
    group_by(Longitude, Latitude, year) %>%
    summarize( Influx_acer = mean(Influx_acer),
               Influx_bet = mean(Influx_bet),
               Influx_gle = mean(Influx_gle),
               Influx_mor = mean(Influx_mor),
               Influx_plat = mean(Influx_plat),
               Influx_que = mean(Influx_que),
               Influx_ulm = mean(Influx_ulm),
               Influx_pop = mean(Influx_pop),
               Influx_jug = mean(Influx_jug),
               Influx_trees	= mean(Influx_trees),
               n = n()) %>%
    filter(!is.na(Influx_trees))
  
  
taub2014 <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/Kate_Guy_tauber_2014_rawcounts.csv") %>% 
  janitor::remove_empty(which = "rows") %>% 
  mutate(Influx_acer = ((maple_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_bet = ((birch_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_gle = ((honeylocust_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_mor = ((morus_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_plat = ((sycamore_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_que = ((oak_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_ulm = ((elm_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_pop = ((populus_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_jug = ((walnut_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_trees = ((total_tree_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_grass = ((grass_count_2014 * mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014,
         Influx_herb =  (((aster_count_2014 + cheno_count_2014 + plantago_count_2014 + grass_count_2014 +
                             mugwort_count_2014 + ragweed_count_2014 + vibarnum_count_2014) * 
                            mean_lycopodium_tablet_2014)/lycopodium_count_2014)/tauber_area_2014) %>% 
  left_join(., taub_locations) %>%  #missing coordinates for the site that is listed just as "law school" - no meta-data for this site
  select(nyccas_id, Longitude, Latitude, year, contains("Influx_")) %>% 
  #take an average of the replicates when available
  group_by(Longitude, Latitude, year) %>%
  summarize( Influx_acer = mean(Influx_acer),
             Influx_bet = mean(Influx_bet),
             Influx_gle = mean(Influx_gle),
             Influx_mor = mean(Influx_mor),
             Influx_plat = mean(Influx_plat),
             Influx_que = mean(Influx_que),
             Influx_ulm = mean(Influx_ulm),
             Influx_pop = mean(Influx_pop),
             Influx_jug = mean(Influx_jug),
             Influx_trees	= mean(Influx_trees),
             n = n()) %>%
  filter(!is.na(Influx_trees)) 
 
bind_rows(taub2013, taub2014) %>% 
ggplot(aes(x = Influx_gle)) + geom_histogram() + facet_wrap(~year)

taub_f <- bind_rows(taub2013, taub2014) %>% 
  #ggplot(aes(x = Longitude, y = Latitude, color = Influx_plat)) + geom_point() + scale_color_viridis_c(trans = "log10") + facet_wrap(~year)
  
  #average across years
  group_by(Longitude, Latitude) %>%
  summarize( Influx_acer = mean(Influx_acer),
             Influx_bet = mean(Influx_bet),
             Influx_gle = mean(Influx_gle),
             Influx_mor = mean(Influx_mor),
             Influx_plat = mean(Influx_plat),
             Influx_que = mean(Influx_que),
             Influx_ulm = mean(Influx_ulm),
             Influx_pop = mean(Influx_pop),
             Influx_jug = mean(Influx_jug),
             Influx_trees	= mean(Influx_trees),
             n = sum(n)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) #%>%   filter(n > 1)



### load in pollen production rasters ###############################

  #load in pollen production rasters
  prod_rast <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/production_1ha_Quercus.tif")
  plot(prod_rast)
  
### threshold function to compare airborne pollen to predicted pollen within X meters ######################

fun_thresh <- function(focal_genus, param_dist){ #focal_genus <- "Juglans"
  
  #load in the pollen production raster for that genus 
  prod_rast_focal_raw <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/production_1ha_",
  focal_genus, ".tif"))
  prod_rast_focal_raw <- subst(prod_rast_focal_raw, NA, 0)
    
  #choose airborne pollen genus
  taub_f_focal <- taub_f %>% 
    mutate(Influx_focal = case_when(focal_genus == "Quercus" ~ Influx_que,
                                    focal_genus == "Acer" ~ Influx_acer,
                                    focal_genus == "Platanus" ~ Influx_plat,
                                    focal_genus == "Betula" ~ Influx_bet,
                                    focal_genus == "Gleditsia" ~ Influx_gle,
                                    focal_genus == "Morus" ~ Influx_mor,
                                    focal_genus == "Ulmus" ~ Influx_ulm,
                                    focal_genus == "Populus" ~ Influx_pop,
                                    focal_genus == "Juglans" ~ Influx_jug,
                                    focal_genus == "trees" ~ Influx_trees)) #Influx_plat Influx_acer Influx_bet Influx_que  Influx_trees

    
  # Create a circular focal window
  focal_matrix <- focalMat(prod_rast_focal_raw, d = param_dist, type = "circle", fillNA = TRUE)
  #focal_matrix <- focalMat(prod_rast_focal_raw, d = 400, type = "circle", fillNA = TRUE)
  focal_matrix_no_weights <- focal_matrix
  focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window
  
  #calculate pollen production within distance
  prod_rast_focal <-  focal( prod_rast_focal_raw, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
  names(prod_rast_focal) <- "prod_within_dist"
  #return(prod_rast_focal)
  
  #extract data from production surface 
  taub_f_comp <- taub_f_focal %>% 
    cbind(., poll_prod_genus = terra::extract(prod_rast_focal, taub_f, ID = FALSE))
  
  #linear model
  if(sum(taub_f_comp$prod_within_dist > 0, na.rm = TRUE) > 5){
  fit <- lm(Influx_focal ~ 
              prod_within_dist, data = taub_f_comp)
  lm_r2 <- summary(fit)$r.squared
  lm_r2_adj <- summary(fit)$adj.r.squared
  aic_value <- AIC(fit)} else { lm_r2 <- NA} #protect against empty dataframes for sparse taxa
  
  ## return aic or r2
  return(lm_r2)
  
   }




### apply function to compare airborne pollen to predicted pollen production at X distance ##########################################

#set up dataframe to hold results  
dist_df <- expand_grid(focal_genus =  c("Quercus", "Platanus", "Acer", "Betula", "Gleditsia", "Morus", "Ulmus", "Populus", "Juglans"), #   
                       param_dist = seq.int(100, 5000, by = 100)) 
    
#apply the comparison threshold function to each row of the dataframe
  dist_result <-  dist_df %>% 
    mutate(r2 = map2_dbl(focal_genus, param_dist, fun_thresh))
  
  # for loop version for error checks 
  # for(k in 2:2){
  #   dist_df$r2[k] <- fun_thresh(focal_genus = dist_df$focal_genus[k], param_dist = dist_df$param_dist[k])  #test <- fun_thresh(500)
  # #output_r2[i] <- fun_compare(test)
  # }
  dist_result %>% group_by(focal_genus) %>% slice_max(r2)

### Fig SI 4: Relationship between airborne pollen and pollen measurements across threshold distances
ggplot(dist_result, aes(x = param_dist, y = r2, color = focal_genus)) + geom_point() + geom_line() + theme_bw() + 
  scale_color_discrete(name = "genus")+ xlab("threshold distance (m)") + ylab(bquote("linear model fit (R"^2~")"))+
  theme(legend.text = element_text(face = "italic"))

write_csv(dist_result, "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/dist_r2_by_genus.csv")
dist_result <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/dist_r2_by_genus.csv")



### scatter plot of pollen production to airborne pollen ###################################

air_vs_prod_fun <- function(focal_genus){

  #focal_genus <- "Juglans"
  #load in the pollen production raster for that genus 
  prod_rast_focal_raw <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/production_1ha_",
                                     focal_genus, ".tif"))
  # prod_rast_focal_raw <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
  #                                   "production_within_400m_Quercus.tif"))
  #plot(prod_rast_focal_raw)
  prod_rast_focal_raw <- subst(prod_rast_focal_raw, NA, 0)
  
  
  #choose airborne pollen genus
  taub_f_focal <- taub_f %>% 
    mutate(Influx_focal = case_when(focal_genus == "Quercus" ~ Influx_que,
                                    focal_genus == "Acer" ~ Influx_acer,
                                    focal_genus == "Platanus" ~ Influx_plat,
                                    focal_genus == "Betula" ~ Influx_bet,
                                    focal_genus == "Gleditsia" ~ Influx_gle,
                                    focal_genus == "Morus" ~ Influx_mor,
                                    focal_genus == "Ulmus" ~ Influx_ulm,
                                    focal_genus == "Populus" ~ Influx_pop,
                                    focal_genus == "Juglans" ~ Influx_jug,
                                    focal_genus == "trees" ~ Influx_trees)) #Influx_plat Influx_acer Influx_bet Influx_que  Influx_trees
  
  #get pollen production within the best distance
  #dist_result <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/dist_r2_by_genus.csv")
  focal_genus2 <- focal_genus #work around to to prevent data masking error with filter
  focal_genus_dist <- dist_result %>% 
    filter(focal_genus == focal_genus2) %>% 
    slice_max(r2) %>% 
    pull(param_dist)
    
  #Create a circular focal window
  focal_matrix <- focalMat(prod_rast_focal_raw, d = focal_genus_dist, type = "circle", fillNA = TRUE)
  #focal_matrix <- focalMat(prod_rast_focal_raw, d = 400, type = "circle", fillNA = TRUE)
  focal_matrix_no_weights <- focal_matrix
  focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window
  
  #calculate pollen production within distance
  prod_rast_focal <-  focal( prod_rast_focal_raw, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
  names(prod_rast_focal) <- "prod_within_dist"
  #return(prod_rast_focal)
  
  #extract data from production surface 
  taub_f_plot <- taub_f_focal %>% 
    cbind(., poll_prod_focal = terra::extract(prod_rast_focal, taub_f, ID = FALSE)) %>% 
    rename(poll_prod_focal = prod_within_dist)
  
  #compare measurements to production surface
  fit <- lm(Influx_focal ~ poll_prod_focal, data = taub_f_plot)
  fit_r2 <- round(summary(fit)$r.squared,2)
  # fit_pval <- summary(fit)$coefficients[, "Pr(>|t|)"][2]
  # fit_pval <- ifelse(fit_pval < 0.001, "0.001", round(p_value, 3))
  
  xlab_text <- paste("pollen production within", focal_genus_dist, "m (trillions of grains)")
  air_prod_scatter <-
    ggplot(taub_f_plot, aes(x = poll_prod_focal/1000 , y = Influx_focal )) + geom_point() + theme_bw() + geom_smooth(method = "lm") +
    xlab(xlab_text) + ylab(bquote(airborne~pollen~(grains/cm^2)))+
    annotate("text", hjust = -0.5, vjust = 1.5,x = -Inf, y = Inf, label = focal_genus, fontface = "italic") + 
    annotate("text", hjust = -0.5, vjust = 2.5,x = -Inf, y = Inf, 
             label = paste0("R^2 == ", fit_r2), parse = TRUE) 
    # annotate("text", hjust = -0.5, vjust = 5.5,x = -Inf, y = Inf, 
    #          label = paste0("p < ", fit_pval), parse = FALSE) 
  
  #a more detailed map for comparison
  # air_prod_map <-
  #   ggplot() + ggthemes::theme_few() + ggtitle("Quercus") + 
  #   geom_spatraster_rgb(data = nyc_topo_spatrast) +
  #   geom_spatraster(data = prod_400m_focal/1000, alpha = 0.6) +
  #   geom_sf(data = taub_f_plot, aes(color = Influx_focal), size = 3) + 
  #   scale_fill_viridis_c(na.value = "transparent", 
  #                        #option = "plasma",
  #                        name = "pollen produced within 1 km \n (trillions of grains)",
  #                        labels = scales::label_comma())+
  #   scale_color_viridis_c(option = "magma", name = bquote(airborne~pollen~(grains/cm^2)))

return(air_prod_scatter)
}


scatter_panels <- map( c("Quercus", "Juglans","Ulmus", "Platanus", "Populus", "Acer", "Morus", "Gleditsia", "Betula"), air_vs_prod_fun)
plot_grid(plotlist = scatter_panels)


#plot(prod_1km_focal_sum)
# writeRaster(prod_1km_focal_sum, 
#             paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
#                    "production_within_1km_", focal_genus, ".tif"), overwrite = TRUE)




# ### function to compare map of pollen production to airborne pollen ###################################
# 
# fun_compare <- function(prod_rast_focal){
#   #extract data from production surface 
#   taub_f_comp <- taub_f %>% 
#     cbind(., poll_prod_genus = terra::extract(prod_rast_focal, taub_f, ID = FALSE))
#   
#   #linear model
#   fit <- lm(poll_prod_genus ~ prod_within_dist, data = taub_f_comp)
#   lm_r2 <- summary(fit)$r.squared
#   aic_value <- AIC(fit)
#   
#   ## return aic or r2
#   return(lm_r2)
# }

# ### exponential function #################### 
#   fun_exp <- function(params){ #param_dist is in 100 m 
#     
#     # Create a circular exponential focal matrix
#     # Parameters
#       # radius_orig <- params[1]  # radius of the circular focal matrix
#       # radius_int <- round(radius)
#       # radius <- if(radius_int %% 2 == 0){radius_int + 1} else {radius_int}
#       
#       radius <- 15 #temporarily hard coding
#       decay <- params[1]  # decay rate for exponential distribution
#     
#       
#     # Calculate matrix size
#       size <- 2 * radius + 1
#     
#     # Calculate the center of the matrix
#      center <- radius + 1
#     
#     # Create distance matrix from center
#       focal_mat <- matrix(NA, nrow = size, ncol = size)
#     
#     for (i in 1:size) {
#       for (j in 1:size) {
#         # Calculate Euclidean distance from center
#         distance <- sqrt((i - center)^2 + (j - center)^2)
#         
#         # Only assign weights within the circular radius
#         if (distance <= radius) {
#           # Apply exponential decay
#           #focal_mat[i, j] <- exp(-decay * distance)
#           
#           #sub in Rayleigh
#           focal_mat[i, j] <- (distance / decay^2) * exp(-distance^2 / (2 * decay^2))
#           
#         }
#       }
#     }
#     #focal_mat
#     
#     #calculate pollen production within distance
#     prod_rast_focal <-  focal( prod_rast, w = focal_mat, fun = "sum", na.rm = TRUE)
#     names(prod_rast_focal) <- "prod_within_dist"
#     
#     #extract data from production surface 
#     taub_f_comp <- taub_f %>% 
#       cbind(., poll_prod_que = terra::extract(prod_rast_focal, taub_f, ID = FALSE))
#     
#     #linear model
#     fit <- lm(Influx_que ~ prod_within_dist, data = taub_f_comp)
#     lm_r2 <- summary(fit)$r.squared
#     
#     sum_resids <- sum(residuals(fit)^2) # Return the sum of squared residuals
#     ## return r2
#     return(-lm_r2)
#     
#     #return(prod_rast_focal)
#     #plot(prod_rast_focal) #plot(prod_rast)
#     # writeRaster(prod_1km_focal_sum, 
#     #             paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
#     #                    "production_within_1km_", focal_genus, ".tif"), overwrite = TRUE)
#   }
#   
#   
#   # #optimizing the search through parameter space for an exponential function  
#   # initial_params <- c(param_dist = 20, param_1 = 0.5)
#   # 
#   test <- optim(par = c(0.5), fn = fun_exp)
# 
#   # #exponential
#   exp_param_space <- seq(0.1, 1, by = 0.2)
#   thresh_dist_list <- seq.int(1, 50, by = 5) #in hundreds of meters
#   model_run_results <- expand_grid(thresh_dist_list, exp_param_space) %>% 
#     mutate(r2 = NA)
#   
#   for(k in 1:nrow(model_run_results)){
#     test <- fun_exp(param_dist = model_run_results$thresh_dist_list[k],
#                     param_1 = model_run_results$exp_param_space[k])
#     model_run_results$r2[k]  <- fun_compare(test)
#   }
#   
#   model_run_results %>% arrange(-r2)
#   
#   ggplot(model_run_results, aes(x = thresh_dist_list, y = exp_param_space, color = r2)) + geom_point(size = 2) +
#     scale_color_viridis_c() + theme_bw()
#   