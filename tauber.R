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

focal_genus_list <- c("Acer", "Betula", "Gleditsia", "Morus", "Platanus", "Quercus", "Ulmus", "Populus", "Juglans")


### load in pollen sample data ######################################
taub_raw <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/2013_NYCPSdataentry_forDanKatz_data_tab.csv")

#use the same samples in analysis as Kate did
taub_f <- taub_raw %>% 
  #use the same samples in analysis as Kate did
  #   filter(Main_analysis == 1) %>% 
  
  # take an average of the replicates when available
  group_by(Longitude, Latitude) %>% 
  summarize(Influx_trees = mean(Influx_trees),
            Influx_plat = mean(Influx_plat),
            Influx_que	= mean(Influx_que), 
            Influx_acer = mean(Influx_acer),
            Influx_bet = mean(Influx_bet)) %>% 
  filter(!is.na(Influx_trees)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)



### load in pollen production rasters ###############################

  #load in pollen production rasters
  prod_rast <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/production_1ha_Quercus.tif")
  plot(prod_rast)
  
### threshold function ######################
fun_thresh <- function(param_dist){
  
  # Create a circular focal window
  focal_matrix <- focalMat(prod_rast, d = param_dist, type = "circle", fillNA = TRUE)
  #focal_matrix <- focalMat(prod_rast, d = 400, type = "circle", fillNA = TRUE)
  focal_matrix_no_weights <- focal_matrix
  focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window
  
  #calculate pollen production within distance
  prod_rast_focal <-  focal( prod_rast, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
  names(prod_rast_focal) <- "prod_within_dist"
  #return(prod_rast_focal)
  
  #extract data from production surface 
  taub_f_comp <- taub_f %>% 
    cbind(., poll_prod_genus = terra::extract(prod_rast_focal, taub_f, ID = FALSE))
  
  #linear model
  fit <- lm(Influx_que ~ #NEED TO CHANGE THIS TO MATCH INPUT GENUS
              prod_within_dist, data = taub_f_comp)
  lm_r2 <- summary(fit)$r.squared
  lm_r2_adj <- summary(fit)$adj.r.squared
  aic_value <- AIC(fit)
  
  ## return aic or r2
  return(lm_r2)
  
   }




### apply to list of values ##########################################

#threshold  
dist_df <- data.frame(dist = seq.int(100, 3000, by = 100),
                      r2 = NA)
  
for(k in 1:nrow(dist_df)){
  dist_df$r2[k] <- fun_thresh(dist_df$dist[k])  #test <- fun_thresh(500)
#output_r2[i] <- fun_compare(test)
}
  

dist_df  %>% arrange(-r2) %>% head()
  
ggplot(dist_df, aes(x = dist, y = r2)) + geom_point() + theme_bw()



  


### detailed investigation of a particular scenario to compare map of pollen production to airborne pollen ###################################
prod_400m_focal_sum <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                                  "production_within_400m_Quercus.tif"))
#plot(prod_400m_focal_sum)

#extract data from production surface 
taub_f_plot <- taub_f %>% 
  cbind(., poll_prod_que = terra::extract(prod_1km_focal_sum, taub_f, ID = FALSE))


#compare measurements to production surface
ggplot(taub_f_plot, aes(x = prod_within_400m , y = Influx_que )) + geom_point() + theme_bw() + geom_smooth(method = "lm")

fit <- lm(Influx_que ~ prod_within_1km, data = taub_f)
summary(fit)

## save outputs


#a more detailed map for comparison
ggplot() + ggthemes::theme_few() + ggtitle("Quercus") + 
  geom_spatraster_rgb(data = nyc_topo_spatrast) +
  geom_spatraster(data = prod_rast/1000, alpha = 0.6) +
  geom_sf(data = taub_f, aes(color = Influx_que), size = 3) + 
  scale_fill_viridis_c(na.value = "transparent", 
                       #option = "plasma",
                       name = "pollen produced within 1 km \n (trillions of grains)",
                       labels = scales::label_comma())+scale_color_viridis_c(option = "magma")

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