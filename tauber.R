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


### load in pollen production rasters ###############################
#load in pollen production rasters
prod_rast_Quercus <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/production_1ha_Quercus.tif")


### load in pollen sample data ######################################
taub_raw <- read_csv("C:/Users/dsk273/Box/NYC projects/pollen data/2013_NYCPSdataentry_forDanKatz_data_tab.csv")
taub <- st_as_sf(taub_raw,  coords = c("Longitude", "Latitude"), crs = 4326)

#use the same samples in analysis as Kate did
taub_f <- taub %>% 
  filter(Main_analysis == 1) %>% 
  filter(!is.na(Influx_trees))






### compare map of pollen production to airborne pollen ###################################
prod_1km_focal_sum <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                                  "production_within_1km_Quercus.tif"))

#extract data from production surface 
taub_f <- taub_f %>% 
  cbind(., poll_prod_que = terra::extract(prod_1km_focal_sum, taub_f, ID = FALSE))

test <- cbind(taub_f, poll_prod_que)
ggplot(test, aes(x = prod_within_1km , y = Influx_que )) + geom_point() + theme_bw() + geom_smooth(method = "lm")

fit <- lm(test$Influx_que ~ test$prod_within_1km)
summary(fit)


#a more detailed map
ggplot() + ggthemes::theme_few() + ggtitle("Quercus") + 
  geom_spatraster_rgb(data = nyc_topo_spatrast) +
  geom_spatraster(data = prod_1km_focal_sum/1000, alpha = 0.6) +
  geom_sf(data = taub_f, aes(color = Influx_que), size = 3) + 
  scale_fill_viridis_c(na.value = "transparent", 
                       #option = "plasma",
                       name = "pollen produced within 1 km \n (trillions of grains)",
                       labels = scales::label_comma())+scale_color_viridis_c(option = "magma")
