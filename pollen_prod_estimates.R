#this script is for calculating pollen production from individual trees of known or predicted identity
# the original version of this script is available here:
# "C:\Users\dsk273\Box\classes\plants and public health fall 2025\class project analysis\calculating pollen production from NYC tree classification.R"

library(sf)
library(dplyr)
library(readr)
library(basemaps)
library(terra)
library(tidyterra)
library(ggplot2)

#read in prediction polygons from Dave
trees_raw <- st_read("C:/Users/dsk273/Box/Katz lab/NYC/classifications/practical/v1/cls_poly_practical_v1_monthcomp.gpkg")
# head(trees)
# plot(trees[1,1])
# sf::plot_sf(trees[1:2, 1:2])

#which genera will be included in the analysis
genera_with_equations_and_classifications <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus",
                                               "Populus", "Juglans", "Morus")

#process the input dataset
trees <- trees_raw %>% 
  mutate(tree_area  = SHAPE_Area * 0.092903) %>% #convert area to square meters from square feet
  rename(Genus = Genus_Merged, #rename some columns for convenience
         Species = Species_Ref) %>% 
  select(Poly_ID, Genus, Species, tree_area)  %>%  #only retain the relevant columns
  filter(Genus %in% genera_with_equations_and_classifications) #only retain the relevant genera

tr <- trees %>% st_drop_geometry() #remove the geometry column for faster processing times in the analysis


### calculate pollen production for each individual tree, using a loop to propagate uncertainty in pollen production equations ################
for(i in 1:100){ #for the final version of the manuscript, let's bump this up to 1000
  
  ### parameters for canopy area calculations from Katz et al. 2020
  # note that some parameters have been updated from the original publication to
  # allow for more robust results.
  # the script for that is available here: "C:/Users/dsk273/Box/MIpostdoc/trees/pollen per tree/pollen per tree analysis and figs 251031.R"
  
  Acne_param_a <- rnorm(n = 1, mean = 0.49, sd = 0.09)
  Acne_param_b <- rnorm(n = 1, mean = -3.16, sd = 3.72)
  Acpl_param_a <- rnorm(n = 1, mean = 0.04, sd = 0.01)
  Acpl_param_b <- rnorm(n = 1, mean = 0.53, sd = 0.42)
  Acru_param_a <- rnorm(n = 1, mean = 0.10, sd = 0.02)
  Acru_param_b <- rnorm(n = 1, mean = -0.10, sd = 0.64)
  Acsa_param_a <- rnorm(n = 1, mean = 0.05, sd =0.02)
  Acsa_param_b <- rnorm(n = 1, mean = 2.59, sd =2.33)
  Bepa_param_a <- rnorm(n = 1, mean = 1.27, sd = 0.55)
  Bepa_param_b <- rnorm(n = 1, mean = -4.73, sd = 8.92)
  Gltr_param_a <- rnorm(n = 1, mean = 0.89, sd = 0.14)
  Gltr_param_b <- rnorm(n = 1, mean = -6.13, sd = 2.58)
  Juni_param_a <- rnorm(n = 1, mean = 0.66, sd = 0.15)
  Juni_param_b <- rnorm(n = 1, mean = 1.42, sd = 8.30)
  Mosp_param_a <- rnorm(n = 1, mean = 0.035375, sd = 0.005136) 
  Mosp_param_b <- rnorm(n = 1, mean = 24.260318, sd = 0.424734)
  Plac_param_a <- rnorm(n = 1, mean = 1.87, sd = 0.48)
  Plac_param_b <- rnorm(n = 1, mean = -26.43, sd = 16.48)
  Posp_param_a <- rnorm(n = 1, mean = 1.2995, sd = 0.3054) 
  Posp_param_b <- rnorm(n = 1, mean = -34.7360, sd = 54.8874)
  Qusp_param_a <- rnorm(n = 1, mean = 0.97, sd = 0.16)
  Qusp_param_b <- rnorm(n = 1, mean = 17.02, sd = 9.54)
  Qupa_param_a <- rnorm(n = 1, mean = 0.65, sd =0.19)
  Qupa_param_b <- rnorm(n = 1, mean = 8.30, sd = 6.88)
  Ulsp_param_a <- rnorm(n = 1, mean = 1.263, sd = 0.420) 
  Ulsp_param_b <- rnorm(n = 1, mean = -38.904, sd = 98.023)  #rnorm(n = 1, mean = 23.11, sd = 0.15)
  
  it_dbh_genus_np_i <-  
    tr %>%  
    
    #protect against unrealistic large trees having overly large numbers by setting canopy area equal to the max recorded in the underlying dataset
    # note: this is only an issue for non-linear relationships, which have been substituted for linear relationships for everything besides Morus
    mutate(tree_area_c = case_when(
      Genus == "Morus" & tree_area > 162 ~ 162,
      .default = tree_area
    )) %>% 
    
    #calculate per tree pollen production as a function of canopy area
    mutate(per_tree_pollen_prod = case_when(
      Genus == "Acer" & Species == "Acer negundo"  ~ ( Acne_param_a * tree_area_c + Acne_param_b) *0.558, #.558 is the sex ratio,
      Genus == "Acer" & Species == "Acer platanoides" ~ Acpl_param_a * tree_area_c + Acpl_param_b,
      Genus == "Acer" & Species == "Acer rubrum"  ~ ( Acru_param_a * tree_area_c + Acru_param_b) * 0.106, #.106 is the sex ratio
      Genus == "Acer" & Species == "Acer saccharinum"~ Acsa_param_a * tree_area_c + Acsa_param_b, 
      # for Acer trees with unknown species from classification, weight by average basal area of Acer species from 2013 i-Tree Eco sampling 
      Genus == "Acer" & is.na(Species)  ~   0.6311 * (Acpl_param_a * tree_area_c + Acpl_param_b) +  #Acer platanoides
                                            0.1713 * (( Acru_param_a * tree_area_c + Acru_param_b) * 0.106) + #Acer rubrum
                                            0.1668 * (Acsa_param_a * tree_area_c + Acsa_param_b) + #Acer saccharinum
                                            0.0011 * (( Acne_param_a * tree_area_c + Acne_param_b) *0.558), #Acer negundo
      # assuming other species produce no pollen
      # relative basal area of other Acer spp (unknown individuals, palmatum, psuedoplatanus, buergerianum, saccharum, campestre = 0.0377)
      Genus == "Betula"  ~ Bepa_param_a* tree_area_c + Bepa_param_b,
      Genus == "Gleditsia"  ~ Gltr_param_a * tree_area_c + Gltr_param_b,
      Genus == "Juglans"  ~ Juni_param_a * tree_area_c + Juni_param_b,
      Genus == "Morus"  ~ ((exp( Mosp_param_a * tree_area_c + Mosp_param_b) )/1000000000 ) * 0.578, #convert to billions and adjust for sex ratio
      Genus == "Platanus"  ~ Plac_param_a * tree_area_c + Plac_param_b,
      Genus == "Populus"  ~ ( Posp_param_a * tree_area_c + Posp_param_b) * 0.482, #adjust for sex ratio
      Genus == "Quercus" & !is.na(Species)  ~ Qusp_param_a * tree_area_c + Qusp_param_b, #red oaks and other non-pin oaks
      Genus == "Quercus" & Species == "Quercus palustris"  ~ Qupa_param_a * tree_area_c + Qupa_param_b, #pin oaks
      #for Quercus species from classification, weight by average basal area of Quercus species from 2013 i-Tree Eco sampling 
      Genus == "Quercus" & is.na(Species)  ~ 0.6546 * (Qusp_param_a * tree_area_c + Qusp_param_b) + #red oaks and other oaks with species
                                             0.3454 * (Qupa_param_a * tree_area_c + Qupa_param_b), 
      Genus == "Ulmus"  ~ ( Ulsp_param_a * tree_area_c + Ulsp_param_b)
    )) %>% 
    
    #protect against small trees having negative numbers
    mutate(per_tree_pollen_prod = case_when(per_tree_pollen_prod < 0 ~ 0,
                                            per_tree_pollen_prod > 0 ~ per_tree_pollen_prod)) %>% 
    mutate(iter = i ) 
  
  
  ifelse(i == 1,
         it_dbh_genus_np_all <- it_dbh_genus_np_i,
         it_dbh_genus_np_all <- bind_rows(it_dbh_genus_np_all, it_dbh_genus_np_i))
  print(i)
}

#summarize across each iteration to calculate both the mean and the standard deviation in pollen production for each tree
indiv_tree_pol_pred <- it_dbh_genus_np_all %>% 
  group_by(Poly_ID) %>% 
  dplyr::summarize(
    pol_mean = mean(per_tree_pollen_prod), 
    pol_sd = sd(per_tree_pollen_prod)
  ) #head(indiv_tree_pol_pred)

#join the pollen production results back to the version that retains geometry. 
tr_export <- left_join(trees, indiv_tree_pol_pred)  #head(tr_export)

#export files for later stages of the analysis
# st_write(tr_export,  "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251115.gpkg", 
#          driver = "GPKG")

#export a csv that has polygon centroid
#get centroids of polygons
tr_export_centroids <- st_centroid(tr_export) 


#convert tr_export to 
tr_export_centroids_proj <- st_transform(tr_export_centroids, crs = 32618) #convert to EPSG 32618 for UTM 18N

#create columns for the x and y coordinates in the 32618
tr_export_centroids_proj_full <- tr_export_centroids_proj %>% 
  mutate(x_EPSG_32618 = sf::st_coordinates(.)[,1],
         y_EPSG_32618 = sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry(.)

# #create columns for the x and y coordinates in the native projected CRS (EPSG 4269)
# tr_export_centroids <- tr_export_centroids %>% 
#   mutate(x_EPSG_4269 = sf::st_coordinates(.)[,1],
#          y_EPSG_4269 = sf::st_coordinates(.)[,2])

# #Also add unprojected lat and long
# tr_export_centroids_4326 <- st_transform(tr_export_centroids, crs = 4326) %>% 
#   mutate(lat = sf::st_coordinates(.)[,1],
#          lon = sf::st_coordinates(.)[,2]) %>% 
#   st_drop_geometry(.)

#export as a csv
write_csv(tr_export_centroids_proj_full, 
          "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251117.csv")



### calculate pollen production within 1 km for each genus #####################################

    # Get extent of NYC
    bbox <- st_bbox(tr_export_centroids_proj)
    
    # download basemaps
      nyc_topo_rast <- basemap_raster(bbox, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
      nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
      

    # Create empty raster template with 100m resolution
    raster_template <- rast(
      xmin = bbox["xmin"],
      xmax = bbox["xmax"],
      ymin = bbox["ymin"],
      ymax = bbox["ymax"],
      resolution = 100,  # 100 meters
      crs = st_crs(tr_export_centroids_proj)$wkt
    )
    
    focal_genus_list <- unique(tr_export_centroids_proj$Genus)
    
    for(i in 1:length(focal_genus_list)){
    focal_genus <-  focal_genus_list[i] 
    #focal_genus <- "Platanus" #Morus Acer Gleditsia Platanus
    tr_export_centroids_proj_quercus <- filter(tr_export_centroids_proj, Genus == focal_genus)
    
    # Convert sf to SpatVector for terra
    tr_vect <- vect(tr_export_centroids_proj_quercus)
    
    # Rasterize using population density (per sq mile)
    prod_raster <- rasterize(
      tr_vect, 
      raster_template, 
      field = "pol_mean",
      fun = "sum"  
    ) #plot(prod_raster)
    
    
    # Create a circular focal window
    focal_matrix <- focalMat(prod_raster, d = 1000, type = "circle", fillNA = TRUE)
    focal_matrix_no_weights <- focal_matrix
    focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window

    #calculate pollen production within 1 km
    prod_1km_focal_sum <-  focal( prod_raster, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
    names(prod_1km_focal_sum) <- "prod_within_1km"
    
    plot(prod_1km_focal_sum)

    writeRaster(prod_1km_focal_sum, 
                paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                "production_within_1km_", focal_genus, ".tif"), overwrite = TRUE)
    
    #a more detailed map
    ggplot() + ggthemes::theme_few() + ggtitle(focal_genus) + 
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = prod_1km_focal_sum/1000, alpha = 0.6) +
      scale_fill_viridis_c(na.value = "transparent", 
                           #option = "plasma",
                           name = "pollen produced within 1 km \n (trillions of grains)",
                           labels = scales::label_comma())
    
    }
