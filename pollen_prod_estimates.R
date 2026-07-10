#this script is for calculating pollen production from individual trees of known or predicted identity
# the original version of this script is available here:
# "C:\Users\dsk273\Box\classes\plants and public health fall 2025\class project analysis\calculating pollen production from NYC tree classification.R"

library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(basemaps)
library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)
library(fs)
library(arrow)
#rm(list = ls())


### prepare work environment: load files and prepare some geographic stuff ########################################################

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

  #get a version of the nyc boundary and convert it to a terra SpatVector
  nyc_boundary2 <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    # st_sf(geometry = st_union(.)) %>% #combine the different boroughs
    st_transform(., crs = 32618) %>% 
    st_union(.) %>% 
    st_sf(geometry = .) %>% 
    vect() #plot(nyc_boundary2)
  
  #load in focal distance for all genera from the tauber.R script
  dist_result <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/dist_r2_by_genus.csv") %>% 
    group_by(focal_genus) %>% 
    slice_max(r2)
  
  
### calculate pollen production per tree, using a loop to propagate uncertainty in pollen production equations, classification ################
for(k in 1:6){ #including uncertainty from classification
 
  #focal_genus_id <- "Quercus" #focal_genus_id <- genus_list[k]
  focal_genus_id <- genera_with_equations[k]
    
  for(j in 1:100){ #100 #dumping out the results locally on a hard drive so I don't run out of ram
    
    #create a list for the production raster lists for each chunk
    r_sum_raster_list_chunk_j <- vector("list", length = 10)
    prod_dist_raster_list_chunk_j <- vector("list", length = 10)
    
    
    for(i in 1:10){ #10
    
      ### parameters for canopy area calculations from Katz et al. 2020
      # note that some parameters have been updated from the original publication to
      # allow for more robust results.
      # the script for that is available here: "C:/Users/dsk273/Box/MIpostdoc/trees/pollen per tree/pollen per tree analysis and figs 251031.R"
      
      #draw from coefficient distributions for each iteration
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
      
    
      #calculate pollen production per tree for the members of the focal genus =================================
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
      
    ### calculate pollen production per 1 ha grid cell for each iteration ========================================
      #calculate per grid cell as a dataframe
     cell_sum_prod_i <- it_dbh_genus_np_i %>%
        inner_join(id_lookup, by = "Poly_ID") %>%       
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
        ) %>% 
       filter(!is.na(grid_cell))
      
     #convert to a raster of total pollen production per ha (mean of each iteration)
     r_sum <- raster_template  # 1 ha raster created earlier
     values(r_sum) <- NA
     r_sum[cell_sum_prod_i$grid_cell] <- cell_sum_prod_i$sum_pol #plot(r_mean)
     
       #fill in NAs with zeros
       #get the cell numbers of the raster that fall within the nyc boundary
       cell_nums <- cells(r_sum, nyc_boundary2)[, "cell"]
       
       # among those cells, find which ones are NA
       na_in_poly <- cell_nums[is.na(r_sum[cell_nums][,1])]
       
       # replace those NA cells with 0
       r_sum[na_in_poly] <- 0
       #plot(r_sum) #check that NA values were filled in
     
       #save the resulting raster into a list
       r_sum_raster_list_chunk_j[[i]] <- wrap(r_sum) #use wrap to make it survive saving
       
       
     ### calculate average pollen production/ha within focal distance for each 1 ha cell ===================================
       # Create a circular focal window
       focal_genus_dist <- filter(dist_result, focal_genus == focal_genus_id) %>% 
         pull(param_dist)
       focal_matrix <- focalMat(r_sum, d = focal_genus_dist, type = "circle", fillNA = TRUE)
       focal_matrix_no_weights <- focal_matrix
       focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1    # Replace all values > 0 with 1 to create an unweighted window
       
       #calculate pollen production within focal_distance
       prod_dist_focal_r_sum <-  focal( r_sum, w = focal_matrix_no_weights, fun = "sum", na.rm = TRUE)
       names(prod_dist_focal_r_sum) <- "mean_prod_ha_within_dist"
       focal_window_area <- (focal_genus_dist^2 * pi)/10000 #divide by 10000 to convert to ha
       prod_dist_focal_rate_ha <- prod_dist_focal_r_sum/focal_window_area  #plot(prod_dist_focal_rate_ha)
       
       #save the resulting raster into a list
       prod_dist_raster_list_chunk_j[[i]] <- wrap(prod_dist_focal_rate_ha) #use wrap to make it survive saving
     
    } #end the i loop 
    
    
    #save the sum(prod per ha)/cell raster list
    saveRDS(r_sum_raster_list_chunk_j, paste0("C:/Users/dsk273/Desktop/prod_per_ha_raster_list/", "prod_raster_list_",focal_genus_id,"_", j, ".rds"))
    
    #save the mean nearby production raster list
    saveRDS(prod_dist_raster_list_chunk_j, 
            paste0("C:/Users/dsk273/Desktop/prod_within_dist_raster_list/", "prod_dist_raster_list_",focal_genus_id,"_", j, ".rds"))
    
    #save the intermediate results of per-tree production per iteration
    #write_csv(it_dbh_genus_np_all, paste0("C:/Users/dsk273/Desktop/prod_chunk/indiv_tree_pol_pred_",focal_genus_id,"_", j, ".csv"))
    write_parquet(x = it_dbh_genus_np_all, 
                  sink = paste0("C:/Users/dsk273/Desktop/prod_chunk/indiv_tree_pol_pred_",focal_genus_id,"_", j, ".parquet"))
    
    print(paste(focal_genus_id,": chunk", j, "out of 100. Completed:", Sys.time()))
    
  } #end the j loop (chunks of 10 due to RAM constraints)
  
  
  
  
  ### rasters for each genus of prod/ha, sd(prod/ha) ###########################
  rds_files <- dir_ls("C:/Users/dsk273/Desktop/prod_per_ha_raster_list/", glob = "*.rds") %>% 
                stringr::str_subset(focal_genus_id)
  #length(rds_files)  # should be 100
  
  # ---- Read all files, unwrap, and collect into one flat list of 1000 SpatRasters ----
  all_rasters <- vector("list", length = length(rds_files) * 10)
  idx <- 1
  
  for (f in rds_files) {
    wrapped_list <- readRDS(f)          # list of 10 PackedSpatRasters
    unwrapped <- lapply(wrapped_list, rast)  # unwrap each back to SpatRaster
    
    for (r in unwrapped) {
      all_rasters[[idx]] <- r
      idx <- idx + 1
    }
  }
  
  #length(all_rasters)  # should be 1000
  
  # ---- Combine into a single multi-layer SpatRaster ----
  r_stack <- rast(all_rasters)
  #nlyr(r_stack)  # should be 1000
  
  # ---- Compute mean and sd per cell across all 1000 layers ----
  r_stats <- app(r_stack, fun = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd   = sd(x, na.rm = TRUE)))
  names(r_stats) <- c("mean", "sd") #plot(r_stats)
  
  # save result
  writeRaster(r_stats, 
              paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis/", 
                     focal_genus_id,"_final_mean_sd.tif"), overwrite = TRUE)
  
  

  ### mean prod/ha within dist, sd prod/ha within dist ###########################
  rds_files_dist <- dir_ls("C:/Users/dsk273/Desktop/prod_within_dist_raster_list/", glob = "*.rds") %>% 
    stringr::str_subset(focal_genus_id)
  #length(rds_files)  # should be 100
  
  # ---- Read all files, unwrap, and collect into one flat list of 1000 SpatRasters ----
  all_rasters_dist <- vector("list", length = length(rds_files_dist) * 10)
  idx <- 1
  
  for (f in rds_files_dist) {
    wrapped_list <- readRDS(f)          # list of 10 PackedSpatRasters
    unwrapped <- lapply(wrapped_list, rast)  # unwrap each back to SpatRaster
    
    for (r in unwrapped) {
      all_rasters[[idx]] <- r
      idx <- idx + 1
    }
  }
  
  #length(all_rasters)  # should be 1000
  
  # ---- Combine into a single multi-layer SpatRaster ----
  r_stack_dist <- rast(all_rasters)
  #nlyr(r_stack)  # should be 1000
  
  # ---- Compute mean and sd per cell across all 1000 layers ----
  r_stats_dist <- app(r_stack_dist, fun = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd   = sd(x, na.rm = TRUE)))
  names(r_stats_dist) <- c("prod_dist_mean", "prod_dist_sd") #plot(r_stats_dist)
  
  # save result
  writeRaster(r_stats_dist, 
              paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis/", 
                     focal_genus_id,"_final_dist_mean_sd.tif"), overwrite = TRUE)
  
} #end genus loop


  
### summarize for each polygon the mean and SD of pollen production ###################
  genera_with_equations <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus", #genera identified in classification
                             "Populus", "Juglans", "Morus") #genera not identified
  
      tr2 <- trees_raw %>% 
      select(Poly_ID, Genus_Merged) %>% 
      filter(Genus_Merged %in% genera_with_equations) #only retain trees that are potentially relevant genera (removing trees that are known to be otherwise)
    
    #doing this with a genus loop
     tr_prod_focal_genus_lst <- vector("list", length = 6)
    
    for(k in 1:6){
      focal_genus_id  <- genera_with_equations[k]
      
      parquet_files_genus <- dir_ls("C:/Users/dsk273/Desktop/prod_chunk/", glob = "*.parquet") %>% 
        stringr::str_subset(focal_genus_id)
      
      ds <- open_dataset(parquet_files_genus, format = "parquet")
      
      per_tree_prod <- ds |>
        select(Poly_ID, per_tree_pollen_prod) |>
        group_by(Poly_ID) |>
        summarise(
          mean_pollen = mean(per_tree_pollen_prod, na.rm = TRUE),
          sd_pollen   = sd(per_tree_pollen_prod, na.rm = TRUE),
          .groups = "drop"
        ) |>
        collect()
      
      tr_prod_focal_genus_lst[[k]] <- 
        tr2 %>% 
        filter(Genus_Merged == focal_genus_id) %>% #only retain the most likely 
        left_join(., per_tree_prod)
    }
    
    #combine all of the different genera dataframes into a single long dataframe
    trees3 <- bind_rows(tr_prod_focal_genus_lst)
    
    #
    trees3_centroids <- st_centroid(trees3) #simplify trees to points
    trees3_centroids_proj <- trees3_centroids %>% 
      select(Poly_ID, Genus_Merged, mean_pollen, sd_pollen, geom) %>% 
      st_transform(., crs = 32618) 
    
    
    #save the long file (includes point geometry)
    st_write(trees3_centroids_proj, 
             "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis/tree_poly_prod_260710b.gpkg", driver = "GPKG")
  
   
  
  

      
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
  