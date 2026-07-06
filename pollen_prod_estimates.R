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
library(data.table)
library(fs)

#read in prediction polygons from Dave; doing the unzipping directly in the command with /vsizip/
#trees_raw <- st_read("/vsizip/C:/Users/dsk273/Box/Katz lab/NYC/classifications/practical/v2/nyc_class_tree_genus_polygons_v2.zip/nyc_class_tree_genus_polygons_v2.gpkg")
trees_raw <- st_read("C:/Users/dsk273/Desktop/nyc_class_tree_genus_polygons_v2/nyc_class_tree_genus_polygons_v2.gpkg")



# head(trees_raw)
# plot(trees[1,1])
# sf::plot_sf(trees[1:2, 1:2])


#which genera will be included in the analysis
genera_with_equations <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus", #genera identified in classification
                                               "Populus", "Juglans", "Morus") #genera not identified

#process the input dataset
trees <- trees_raw %>% 
  mutate(tree_area  = round(TNC_SHAPE_Area * 0.092903, 2)) %>% #convert area to square meters from square feet
  rename(#Genus = Genus_Merged, #rename some columns for convenience
         Species = Species_Ref) %>% 
  select(Poly_ID, Genus_Predicted, Genus_Ref, Species, tree_area,  #only retain the relevant columns
         Acer, Betula, Gleditsia, Platanus, Quercus, Ulmus)  %>% 
  filter(Genus_Ref %in% genera_with_equations | is.na(Genus_Ref)) #only retain trees that are potentially relevant genera (removing trees that are known to be otherwise)

tr <- trees %>% st_drop_geometry() #remove the geometry column for faster processing times in the analysis



  
### calculate pollen production per tree, using a loop to propagate uncertainty in pollen production equations, classification ################
for(k in 1:6){ #including uncertainty from classification
 
  #focal_genus_id <- "Gleditsia" #focal_genus_id <- genus_list[k]
  focal_genus_id <- genera_with_equations[k]
    
  for(j in 1:100){ #dumping out the results locally on a hard drive so I don't run out of ram
    for(i in 1:10){ #
    
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
    
    

    # df <- tibble(
    #   id = 1:10,
    #   probability = c(0.1, 0.9, 0.5, NA, 0.25, 0.6, 0.4, 0.95, 0.05, 0.8)
    # )
    # 
    # # Simulate binary outcome
    # df <- df %>%
    #   mutate(outcome = rbinom(n(), size = 1, prob = probability))
    # 
    # print(df)
    
    
    it_dbh_genus_np_i <-  
      tr %>%  
      #sample_n(100) %>% 
      
      #assign each classified tree to a genus based on probability
      mutate(focal_genus_prob_NAs = as.numeric(.data[[focal_genus_id]]),
             focal_genus_prob = replace_na(focal_genus_prob_NAs, 0),
             Genus_class_outcome = rbinom(n(), size = 1, prob = focal_genus_prob),
             Genus_it = case_when(Genus_class_outcome == 1 ~ focal_genus_id,
                                  Genus_Ref == focal_genus_id ~ focal_genus_id,
                                  .default = "not focal genus")) %>% 
      
      #filter out trees that are not assigned to the genus on this iteration
      filter(Genus_it == focal_genus_id) %>%
      
      #filter out trees that are known to be another genus
      filter(Genus_Ref == focal_genus_id | is.na(Genus_Ref)) %>% 
      
      #protect against unrealistic large trees having overly large numbers by setting canopy area equal to the max recorded in the underlying dataset
      # note: this is only an issue for non-linear relationships, which have been substituted for linear relationships for everything besides Morus
      mutate(tree_area_c = case_when(
        Genus_it == "Morus" & tree_area > 162 ~ 162,
        .default = tree_area
      )) %>% 
      
      #assign species probabilistically when there are different equations available within the same genus
      mutate(species_prob = runif(n(), 0, 1),
        species_it = case_when(!is.na(Species) ~ Species, 
                               focal_genus_id == "Acer" & Genus_it == "Acer" & species_prob < 0.6311 ~ "Acer platanoides",
                               focal_genus_id == "Acer" & Genus_it == "Acer" & species_prob > 0.6311 & species_prob < 0.8024  ~ "Acer rubrum",
                               focal_genus_id == "Acer" & Genus_it == "Acer" & species_prob > 0.8024 & species_prob < 0.9692  ~ "Acer saccharinum",
                               focal_genus_id == "Acer" & Genus_it == "Acer" & species_prob > 0.9692  ~ "Acer saccharinum",
                               focal_genus_id == "Quercus" & Genus_it == "Quercus" & species_prob < 0.3454  ~ "Quercus palustris",
                               .default = NA
                               )) %>% 
      
      #calculate per tree pollen production as a function of canopy area
      mutate(per_tree_pollen_prod = case_when(
        Genus_it == "Acer" & Species == "Acer negundo"  ~ ( Acne_param_a * tree_area_c + Acne_param_b) *0.558, #.558 is the sex ratio,
        Genus_it == "Acer" & Species == "Acer platanoides" ~ Acpl_param_a * tree_area_c + Acpl_param_b,
        Genus_it == "Acer" & Species == "Acer rubrum"  ~ ( Acru_param_a * tree_area_c + Acru_param_b) * 0.106, #.106 is the sex ratio
        Genus_it == "Acer" & Species == "Acer saccharinum"~ Acsa_param_a * tree_area_c + Acsa_param_b, 
        # for Acer trees with unknown species from classification, weight by average basal area of Acer species from 2013 i-Tree Eco sampling 
        Genus_it == "Acer" & is.na(Species)  ~   0.6311 * (Acpl_param_a * tree_area_c + Acpl_param_b) +  #Acer platanoides
                                              0.1713 * (( Acru_param_a * tree_area_c + Acru_param_b) * 0.106) + #Acer rubrum
                                              0.1668 * (Acsa_param_a * tree_area_c + Acsa_param_b) + #Acer saccharinum
                                              0.0011 * (( Acne_param_a * tree_area_c + Acne_param_b) *0.558), #Acer negundo
        # assuming other species produce no pollen
        # relative basal area of other Acer spp (unknown individuals, palmatum, psuedoplatanus, buergerianum, saccharum, campestre = 0.0377)
        Genus_it == "Betula"  ~ Bepa_param_a* tree_area_c + Bepa_param_b,
        Genus_it == "Gleditsia"  ~ Gltr_param_a * tree_area_c + Gltr_param_b,
        Genus_it == "Juglans"  ~ Juni_param_a * tree_area_c + Juni_param_b,
        Genus_it == "Morus"  ~ ((exp( Mosp_param_a * tree_area_c + Mosp_param_b) )/1000000000 ) * 0.578, #convert to billions and adjust for sex ratio
        Genus_it == "Platanus"  ~ Plac_param_a * tree_area_c + Plac_param_b,
        Genus_it == "Populus"  ~ ( Posp_param_a * tree_area_c + Posp_param_b) * 0.482, #adjust for sex ratio
        Genus_it == "Quercus" & !is.na(Species)  ~ Qusp_param_a * tree_area_c + Qusp_param_b, #red oaks and other non-pin oaks
        Genus_it == "Quercus" & Species == "Quercus palustris"  ~ Qupa_param_a * tree_area_c + Qupa_param_b, #pin oaks
        #for Quercus species from classification, weight by average basal area of Quercus species from 2013 i-Tree Eco sampling 
        Genus_it == "Quercus" & is.na(Species)  ~ 0.6546 * (Qusp_param_a * tree_area_c + Qusp_param_b) + #red oaks and other oaks with species
                                               0.3454 * (Qupa_param_a * tree_area_c + Qupa_param_b), 
        Genus_it == "Ulmus"  ~ ( Ulsp_param_a * tree_area_c + Ulsp_param_b),
        .default = 0 #otherwise assigning the tree to 0 pollen production
      )) %>% 
      
      #round pollen production
      mutate(per_tree_pollen_prod = round(per_tree_pollen_prod, 2)) %>% 
      
      #protect against small trees having negative numbers
      mutate(per_tree_pollen_prod = case_when(per_tree_pollen_prod < 0 ~ 0,
                                              per_tree_pollen_prod > 0 ~ per_tree_pollen_prod)) %>% 
      mutate(iter = i, 
             chunk_j = j,
             focal_genus = focal_genus_id) %>%  
      
      #remove some columns to make the files smaller
      select(Poly_ID, Genus_Ref, Species, tree_area_c, species_it, per_tree_pollen_prod, iter, chunk_j, focal_genus)
    
    
    ifelse(i == 1,
           it_dbh_genus_np_all <- it_dbh_genus_np_i,
           it_dbh_genus_np_all <- bind_rows(it_dbh_genus_np_all, it_dbh_genus_np_i))
    #print(i)
  }
    
    
    write_csv(it_dbh_genus_np_all, paste0("C:/Users/dsk273/Desktop/prod_chunk/indiv_tree_pol_pred_",focal_genus_id,"_", j, ".csv"))
    #write_parquet(x = it_dbh_genus_np_all, sink = paste0("C:/Users/dsk273/Desktop/prod_chunk/indiv_tree_pol_pred_",focal_genus_id,"_", j, ".parquet"))
    print(paste(focal_genus_id,": chunk", j, "out of 100. Completed:", Sys.time()))
  }
} #end genus loop


### process simulated pollen production per grid cell for each genus===============================
  
  ### create a lookup table for tree location and grid cell
      #get coordinates for each tree in EPSG 32618
      tr_export_centroids <- st_centroid(trees) #simplify trees to points
      tr_export_centroids_proj <- tr_export_centroids %>% 
                                  select(Poly_ID) %>% 
                                  st_transform(., crs = 32618) %>% #convert to EPSG 32618 for UTM 18N
                                  mutate(x_EPSG_32618 = sf::st_coordinates(.)[,1],
                                          y_EPSG_32618 = sf::st_coordinates(.)[,2]) #%>%  st_drop_geometry(.)
    
      # Get extent of NYC 
      bbox <- st_bbox(tr_export_centroids_proj)
      
      # Create empty raster template with 100m resolution
      raster_template <- rast(
        xmin = bbox["xmin"] , #
        xmax = bbox["xmax"] ,
        ymin = bbox["ymin"] ,
        ymax = bbox["ymax"] ,
        resolution = 100,  # 100 meters
        crs = st_crs(tr_export_centroids_proj)$wkt
      )
      
      # Get the cell number for each id 
      id_lookup <- tr_export_centroids_proj %>%
        mutate(grid_cell = cellFromXY(raster_template, cbind(x_EPSG_32618, y_EPSG_32618))) %>%
        select(Poly_ID, grid_cell)
   
      #check to make sure that there aren't NAs (i.e., points outside of the area)
        sum(is.na(id_lookup$grid_cell))  
      
      # # Pull their original coordinates to inspect
      # na_coords <- id_coords %>% filter(id %in% na_ids$id)
      # na_coords
  
      all_cells <- id_lookup %>% distinct(grid_cell)
      
      
      # ---- Step 1: function to process one file ----
      process_file <- function(path, out_dir) {
        dt <- fread(path, select = c("Poly_ID", "per_tree_pollen_prod"))  # only read needed columns
       # dt <- fread(files_to_read_genus[1], select = c("Poly_ID", "per_tree_pollen_prod"))  # only read needed columns
        
        summary_tbl <- as_tibble(dt) %>%
          inner_join(id_lookup, by = "Poly_ID") %>%        # inner_join drops unmatched ids; use left_join + check if you need to know about them
          group_by(grid_cell) %>%
          summarize(
            sum_pol = sum(per_tree_pollen_prod, na.rm = TRUE),
            n = n(),
            had_obs = TRUE,   # flag real (non-filled) rows
            .groups = "drop"
          ) %>% 
          complete( #add in zero values at this stage
            grid_cell = all_cells$grid_cell,
            fill = list(sum_pol = 0, n = 0, had_obs = FALSE)
          )
        
        out_path <- path(out_dir, paste0(path_ext_remove(path_file(path)), "_summary.rds"))
        write_rds(summary_tbl, out_path)
        rm(dt, summary_tbl); gc()
        out_path
      }
      
      #process_file(path = files_to_read_genus[1], out_dir = "C:/Users/dsk273/Desktop/prod_chunk_pixel_tibble/")
      
  
      # ---- Run over all files ----
   
      for(j in 1:6){   #start genus loop
        
        #focal_genus_id <- "Acer"
        focal_genus_id <- genera_with_equations[i]
        
        files_to_read <- dir("C:/Users/dsk273/Desktop/prod_chunk/", full.names = TRUE)
        files_to_read_genus <- stringr::str_subset(files_to_read, focal_genus_id)
      
      
      # csv_files <- dir_ls("raw_data/", glob = "*.csv")
        #dir_create("intermediate_summaries/")
        
        # Simple sequential loop with progress
        summary_paths <- character(length(files_to_read_genus))
        #for (i in seq_along(csv_files)) {
        for (i in 1:2) { #this will be up to 100 for all files
          summary_paths[i] <- process_file(path = files_to_read_genus[i], out_dir = "C:/Users/dsk273/Desktop/prod_chunk_pixel_tibble/intermediate_summaries/")
          message(sprintf("Done %d/%d: %s", i, length(files_to_read_genus), files_to_read_genus[i]))
        }
      
        # ---- Step 2: combine and compute final stats ----
        final_summary <- summary_paths[1:2] %>%
          map(read_rds) %>%
          list_rbind() %>%
          group_by(grid_cell) %>%
          summarize(
            mean_sum      = mean(sum_pol, na.rm = TRUE),
            sd_sum        = sd(sum_pol, na.rm = TRUE),
            mean_n        = mean(n, na.rm = TRUE),
            sd_n          = sd(n, na.rm = TRUE),
            files_with_an_obs  = sum(had_obs),   # how many of the 1000 files had any obs in this cell
            .groups = "drop"
          ) %>% filter(!is.na(grid_cell))
        
        # ---- Step 3: save rasters of pollen mean, sd, and n for each genus 
        r_mean <- raster_template  # 1 ha raster created earlier
        values(r_mean) <- NA
        r_mean[final_summary$grid_cell] <- final_summary$mean_sum #plot(r_mean)
        
        r_sd <- raster_template
        values(r_sd) <- NA
        r_sd[final_summary$grid_cell] <- final_summary$sd_sum #plot(r_sd)
        
        r_n <- raster_template
        values(r_n) <- NA
        r_n[final_summary$grid_cell] <- final_summary$mean_n #plot(r_n)
        
        genus_raster_name <- paste0( "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis/",
                                     focal_genus_id, "_1ha_summary_raster.tif")
        writeRaster(c(r_mean, r_sd, r_n), genus_raster_name, names = c("mean", "sd", "mean_n"), overwrite = TRUE)
        print(Sys.time())
      } #end genus loop
      
      
      
    #leaving off here   
      
      
      
      
      
      
      
      
      
      
      

      
      
  # download basemaps
      nyc_topo_rast <- basemap_raster(bbox, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
      nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
      
      
  





#load in the results from temporary local harddrive storage
files_to_read <- dir("C:/Users/dsk273/Desktop/prod_chunk/", full.names = TRUE)

files_to_read_genus <- stringr::str_subset(files_to_read, "Quercus")
it_dbh_genus_np_all <- purrr::map_dfr(files_to_read_genus, read_csv) #test <- read_csv(files_to_read[1]) head(test)

# test <- read_csv(files_to_read_genus[1])
ggplot(it_dbh_genus_np_all, aes(x = tree_area_c, y = per_tree_pollen_prod, color = species_it)) + geom_point() + facet_wrap(~focal_genus)

#summarize across each iteration to calculate both the mean and the standard deviation in pollen production for each tree
indiv_tree_pol_pred <- it_dbh_genus_np_all %>% 
  group_by(Poly_ID) %>% 
  dplyr::summarize(
    pol_mean = mean(per_tree_pollen_prod), 
    pol_sd = sd(per_tree_pollen_prod),
    n = n()
  ) #head(indiv_tree_pol_pred)


test <- filter(it_dbh_genus_np_all, Poly_ID == 7)
test2 <- filter(tr, Poly_ID == 7)
#write_csv(indiv_tree_pol_pred, "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/indiv_tree_pol_pred.csv")

#join the pollen production results back to the version that retains geometry. 
tr_export <- left_join(trees, indiv_tree_pol_pred)  #head(tr_export) #rm(it_dbh_genus_np_all)

#export files for later stages of the analysis
# st_write(tr_export,  "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.gpkg",
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
          "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.csv")






### calculate pollen production on a 1 ha raster for each genus #####################################

#reload data from previous step if necessary
# tr_export_centroids_proj_full <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.csv")
# tr_export_centroids_proj <- st_as_sf(tr_export_centroids_proj_full, coords = c("x_EPSG_32618", "y_EPSG_32618"), crs = 32618)

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
    #focal_genus <- "Quercus" #Morus Acer Gleditsia Platanus
    tr_export_centroids_proj_quercus <- filter(tr_export_centroids_proj, Genus == focal_genus)
    
    # Convert sf to SpatVector for terra
    tr_vect <- vect(tr_export_centroids_proj_quercus)
    
    # Rasterize to pollen production per ha
    prod_raster <- rasterize(
      tr_vect, 
      raster_template, 
      field = "pol_mean",
      fun = "sum"  
    ) #plot(prod_raster)
    
    #export just the pollen production raster
    writeRaster(prod_raster, 
                paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                       "production_1ha_", focal_genus, ".tif"), overwrite = TRUE)
    
    # #a more detailed map of just pollen production
    ggplot() + ggthemes::theme_few() + ggtitle(focal_genus) +
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = prod_raster, alpha = 0.6) +
      scale_fill_viridis_c(na.value = "transparent",
                           #option = "plasma",
                           name = "pollen produced  \n (billions of grains)",
                           labels = scales::label_comma())

    }

    
### calculate pollen production within 400 m for each genus #####################################
  
  #reload data from previous step if necessary
  # tr_export_centroids_proj_full <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.csv")
  # tr_export_centroids_proj <- st_as_sf(tr_export_centroids_proj_full, coords = c("x_EPSG_32618", "y_EPSG_32618"), crs = 32618)
  
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
    #focal_genus <- "Quercus" #Morus Acer Gleditsia Platanus
    tr_export_centroids_proj_quercus <- filter(tr_export_centroids_proj, Genus == focal_genus)
    
    # Convert sf to SpatVector for terra
    tr_vect <- vect(tr_export_centroids_proj_quercus)
    
    # Rasterize to pollen production per ha
    prod_raster <- rasterize(
      tr_vect, 
      raster_template, 
      field = "pol_mean",
      fun = "sum"  
    ) #plot(prod_raster)
    
    # Create a circular focal window
    focal_matrix <- focalMat(prod_raster, d = 400, type = "circle", fillNA = TRUE)
    focal_matrix_no_weights <- focal_matrix
    focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window
    
    #calculate pollen production within 400 m
    prod_1km_focal_sum <-  focal( prod_raster, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
    names(prod_1km_focal_sum) <- "prod_within_400m"
    
    plot(prod_1km_focal_sum)
    
    writeRaster(prod_1km_focal_sum, 
                paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                       "production_within_400m_", focal_genus, ".tif"), overwrite = TRUE)
    
    #a more detailed map
    ggplot() + ggthemes::theme_few() + ggtitle(focal_genus) + 
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = prod_1km_focal_sum/1000, alpha = 0.6) +
      scale_fill_viridis_c(na.value = "transparent", 
                           #option = "plasma",
                           name = "pollen produced within 400m \n (trillions of grains)",
                           labels = scales::label_comma())
    
  }
  
  
  
### compare tree pollen production with and without tree classification data ###########################
  #reload data from previous step if necessary
  # tr_export_centroids_proj_full <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.csv")
  # tr_export_centroids_proj <- st_as_sf(tr_export_centroids_proj_full, coords = c("x_EPSG_32618", "y_EPSG_32618"), crs = 32618)
  
tr_export_centroids_proj %>% 
    st_drop_geometry() %>% 
    mutate(is_st_tree = case_when(is.na(Species) ~ "classified",
                                  !is.na(Species) ~ "street tree")) %>% 
    group_by(is_st_tree, Genus) %>% 
    summarize(total_pollen = sum(pol_mean, na.rm = TRUE)) %>%
    mutate(rel_pollen = total_pollen/sum(tr_export_centroids_proj$pol_mean, na.rm = TRUE) * 100) %>% 
    ggplot(aes(x = Genus, y = total_pollen/1000000, fill = is_st_tree)) + geom_col(position = position_dodge2(preserve = "single")) + theme_bw() + 
    ylab("pollen production (quadrillions of grains)") + 
    scale_fill_discrete(name = "data source") + theme(axis.text.x = element_text(face = "italic"))
  
  
### maximum pollen production within a 1 ha pixel for each genus #########################################
  focal_genus <- "Acer"
  # focal_raster <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
  #                    "production_within_400m_", focal_genus, ".tif"))
  focal_raster <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                                  "production_1ha_", focal_genus, ".tif"))
  max(values(focal_raster), na.rm = TRUE) * 1000000000
  