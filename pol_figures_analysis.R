# Project: Pollen production and exposure in NYC
# Class: PLSCI 4450/6450 Urban Plants and Public Health
# Contact: Dan Katz, Cornell University
# 
# this script is for analyzing and visualizing pollen production in NYC and per tree exposure


#set up work environment
library(tidycensus)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(readr)
library(tidyterra)
library(ggspatial)
library(basemaps)
library(purrr)
library(cowplot) 
library(patchwork)
library(grid)

### load in tree pollen production estimates
# these were created in the 'pollen_prod_estimates.R' script
tr_export_centroids_proj_full_csv <- read_csv( "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251125.csv")
tr_export_centroids_proj <- st_as_sf(tr_export_centroids_proj_full_csv, coords = c("x_EPSG_32618" , "y_EPSG_32618"),
                                          crs = 32618)

## load in census population density
# these were created in the 'census_map.R' script
density_focal_sum <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_400m_focal_sum.tif")

## extract population density for each tree 
tr_export_centroids_proj

polpop <- tr_export_centroids_proj %>% 
          mutate(pol_mean_orig = pol_mean,
                 pop_within_1_km = terra::extract(density_focal_sum, tr_export_centroids_proj)[,2],
                 pop_pol = pol_mean * pop_within_1_km * 1000000000) # calculate pollen production * population to get potential impact

polpop_df <- polpop %>% st_drop_geometry(.)

### statistics for paper ###############################################################
#read in prediction polygons from Dave, takes a minute to load
trees_raw <- st_read("C:/Users/dsk273/Box/Katz lab/NYC/classifications/practical/v1_1/cls_poly_practical_v1_1_monthcomp_filt_w_crown_metrics.gpkg")

genera_with_equations <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus", "Populus", "Juglans", "Morus")

#process the input dataset
trees_stats <- trees_raw %>% 
  mutate(tree_area  = TNC_SHAPE_Area * 0.092903) %>% #convert area to square meters from square feet
  rename(Genus = Genus_Merged, #rename some columns for convenience
         Species = Species_Ref) %>% 
  select(Poly_ID, Genus, Species, tree_area)  %>%  #only retain the relevant columns
  mutate(pollen_calc = case_when(Genus %in% genera_with_equations ~ "pollen_calculated",
                                 !(Genus %in% genera_with_equations) ~ "pollen_not_calculated")) %>% 
 st_drop_geometry() #remove the geometry column for faster processing times in the analysis

#what portion of canopy area do we have 
trees_stats %>% 
  group_by(pollen_calc) %>% 
  summarize(total_area = sum(tree_area),
            n_trees = n())

#total pollen production by genus
tr_export_centroids_proj_full_csv %>% 
  group_by(Genus) %>% 
  summarize(total_pol = sum(pol_mean, na.rm = TRUE) * 1000000000,
            n_trees = n()) %>% 
  arrange(-total_pol)



### visualize results ###################################################

# #Histogram of pollen exposure for each tree by genus
# ggplot(polpop_df, aes(x = pop_pol)) + geom_histogram() + facet_wrap(~Genus, scales = "free") + theme_bw() + scale_x_log10()

#cdf plot
ggplot(polpop_df, aes(x = pop_pol + 1, color = Genus)) + stat_ecdf(geom = "step", pad = FALSE) + theme_bw() + scale_x_log10() +
  xlab(expression("exposure per tree (pollen grains" %*% "population within 400 m)")) + ylab("cumulative distribution (proportion)") +
  theme(legend.text = element_text(face = "italic"))

# #map of a particular genus
# polpop_sub <- polpop %>% 
#   filter(Genus == "Quercus") %>% 
#   slice_sample(n = 20000) %>% #to reduce plotting time
#   filter(tree_area > 125 & tree_area < 150)

# #map of pollen production by a particular genus
# polpop_sub <- polpop %>% 
#   filter(Genus == "Quercus") %>% 
#   slice_sample(n = 10000) # %>%  filter(tree_area > 125 & tree_area < 150)
# 
# ggplot() + 
#   geom_sf(data = polpop_sub, fill = "lightgray") + # Draw the base map
#   geom_sf(data = polpop_sub, aes(color = pol_mean), size = 1, alpha = 0.1) + # Add points
#   theme_minimal() + 
#   scale_color_viridis_c(name = "pollen produced by tree x population within 400 m")# Apply a minimal theme


#comparison of tree size to pop_pol
ggplot(polpop, aes(x = tree_area, y = pop_pol + 1)) + geom_hex(name = "density") + facet_wrap(~Genus) + theme_bw() +
  scale_fill_viridis_c(trans = "log10", name = "number of trees") + xlab(bquote("tree area ("~m^2~")")) + ylab("per-tree exposure (people within 400 m x pollen)") +
  scale_y_log10() + scale_x_continuous(limits = c(0, 1000)) +  theme(strip.text = element_text(face = "italic"))
       

# # map of tree pollen within 400 m for a single genus
#     prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
#                        "production_within_400m_", "Quercus", ".tif"))
#     plot(prod_m_focal_mean)
#     hist(prod_m_focal_mean)
#     
#     #load in nyc boundary polygon
#     nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
#       st_union() %>% #combine the different boroughs
#       st_transform(., crs = 32618)
#         #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
#         # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
#     
#     nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
#     nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
# 
#     # create map
#       ggplot() + ggthemes::theme_few() +   
#       geom_spatraster_rgb(data = nyc_topo_spatrast) +
#       geom_spatraster(data = prod_m_focal_mean/1000, alpha = 0.6) +
#       scale_fill_viridis_c(na.value = "transparent", 
#                            #option = "magma",
#                            name = paste0("mean pollen production \n(trillions of grains/ha within ", focal_genus_dist, "m)"),
#                            labels = scales::label_comma()) +
#         annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
#                          bar_cols = c("black", "white"), # Colors of the scale bar segments
#                          style = "ticks",
#                          text_cex = 0.8) +  # Text size for the scale bar label
#         annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
#                                pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
#         theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
#                 legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)
# 

        
        
        
### Fig 2: panel figure of 1 ha, sd, and mean pollen production within focal distance for all tree taxa #############################################################
 #pre load map elements
    #load in nyc boundary polygon
      nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
        st_union() %>% #combine the different boroughs
        st_transform(., crs = 32618)
      #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
      # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
      
      nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
      nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
      
    #load in the distance for each genus
      dist_result <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/dist_r2_by_genus.csv") %>% 
        group_by(focal_genus) %>% 
        slice_max(r2)
      
  #function to create a map panel for: sum production per ha 
    fun_sum_1ha_prod_map_genus <- function(focal_genus){    
      focal_genus_id <- focal_genus #focal_genus_id <- "Quercus"
      
      prod_ha_rast_name <- dir("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis", full.names = TRUE) %>% 
        stringr::str_subset("1ha_summary_raster_v2.tif") %>% 
        stringr::str_subset(focal_genus_id)
      
      prod_ha_rast <- rast(prod_ha_rast_name, lyr = "mean")
      
      # #shrink everything above the 99th percentile of values to 99% for more effective visualization
      # p99 <- unlist(global(prod_ha_rast, fun = quantile, probs = 0.99, na.rm = TRUE))
      # prod_ha_rast[prod_ha_rast[] > p99] <- p99
      
      # create map
      focal_map_panel <-   
        ggplot() + #ggthemes::theme_few() +   
        geom_spatraster_rgb(data = nyc_topo_spatrast) +
        geom_spatraster(data = prod_ha_rast) +
        scale_fill_viridis_c(na.value = "transparent", 
                             option = "turbo",
                             name = expression(atop(atop(textstyle("mean pollen"),
                                                         textstyle("production")),
                                                    "("~10^9~"grains/ha)")),
                             labels = scales::label_comma()) +
        theme(  
          legend.position = c(0.05, 0.95),  # Places the legend at the top-left corner
          legend.justification = c(0.05, 0.95),# Aligns the legend box to its top-left corner)
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
          #plot.background = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
      
      return(focal_map_panel)
    } #end map making function
    
  #function to create a map panel for: sd production per ha 
  fun_sd_1ha_prod_map_genus <- function(focal_genus){    
      focal_genus_id <- focal_genus #focal_genus_id <- "Quercus"
      
      sd_ha_rast_name <- dir("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis", full.names = TRUE) %>% 
        stringr::str_subset("1ha_summary_raster_v2.tif") %>% 
        stringr::str_subset(focal_genus_id)
      
      sd_ha_rast <- rast(sd_ha_rast_name, lyr = "sd")
      
      # #shrink everything above the 99th percentile of values to 99% for more effective visualization
      # p99 <- unlist(global(sd_ha_rast, fun = quantile, probs = 0.99, na.rm = TRUE))
      # sd_ha_rast[sd_ha_rast[] > p99] <- p99
      
      cv_ha_rast <- (sd_ha_rast/prod_ha_rast) * 100
      ggplot() +  geom_spatraster(data = prod_ha_rast) + scale_fill_viridis_c( option = "turbo")
      
      # create map
      focal_map_panel <-   
        ggplot() + #ggthemes::theme_few() +   
        geom_spatraster_rgb(data = nyc_topo_spatrast) +
        geom_spatraster(data = sd_ha_rast) +
        scale_fill_viridis_c(na.value = "transparent", 
                             option = "turbo",
                             name = expression(atop(atop(textstyle("SD pollen"),
                                                         textstyle("production")),
                                                    "("~10^9~"grains/ha)")),
                             labels = scales::label_comma()) +
        theme(  
          legend.position = c(0.05, 0.95),  # Places the legend at the top-left corner
          legend.justification = c(0.05, 0.95),# Aligns the legend box to its top-left corner)
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
          #plot.background = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
      
      return(focal_map_panel)
    } #end map making function
    
  #create a function that produces a map panel of mean pollen production within focal distance
    fun_mean_prod_focal_dist_map_genus <- function(focal_genus){    
        focal_genus_id <- focal_genus #focal_genus_id <- "Quercus"
      
    
      # load raster of mean pollen production within x m for a genus
        focal_genus_dist <- filter(dist_result, focal_genus == focal_genus_id) %>% 
          pull(param_dist)
      prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/July26_reanalysis/",
                                       "mean_production_within_",focal_genus_dist, "m_", focal_genus_id, ".tif"))
 
      #shrink everything above the 99th percentile of values to 99% for more effective visualization
      p99 <- unlist(global(prod_m_focal_mean, fun = quantile, probs = 0.99, na.rm = TRUE))
      prod_m_focal_mean[prod_m_focal_mean[] > p99] <- p99
      
      # create map
    focal_map_panel <-   
      ggplot() + #ggthemes::theme_few() +   
        geom_spatraster_rgb(data = nyc_topo_spatrast) +
        geom_spatraster(data = prod_m_focal_mean) +
        scale_fill_viridis_c(na.value = "transparent", 
                             option = "turbo",
                             name = expression(atop(atop(textstyle("neighborhood"),
                                                    textstyle("pollen prod.")),
                                               "("~10^9~"grains/ha)")),
                             labels = scales::label_comma()) +
              theme(  
                legend.position = c(0.05, 0.95),  # Places the legend at the top-left corner
                legend.justification = c(0.05, 0.95),# Aligns the legend box to its top-left corner)
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 8),
                axis.title=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(), 
                plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
                #plot.background = element_blank(),
                plot.background = element_rect(fill = "white", color = "white"),
                panel.background = element_rect(fill = "white"),
                plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
    
    return(focal_map_panel)
  } #end map making function
  
  

  
#create list of maps
  focal_genus_panels <- c("Quercus", "Platanus", "Ulmus", "Acer")
  
  map_list_1ha_sum <- map(focal_genus_panels, fun_sum_1ha_prod_map_genus)
  map_list_1ha_sd <- map(focal_genus_panels, fun_sd_1ha_prod_map_genus)
  map_list_1ha_mean_dist <- map(focal_genus_panels, fun_mean_prod_focal_dist_map_genus)
  
  map_list <- c(map_list_1ha_sum, map_list_1ha_sd, map_list_1ha_mean_dist)
  
  # replacing old cowplot approach with patchwork for greater control
  # plot_grid(plotlist = map_list, ncol = 3,
  #           # labels = c( bquote(paste("A) "~italic("Quercus"))),
  #           #             bquote(paste("B) "~italic("Platanus"))),
  #           #             bquote(paste("C) "~italic("Ulmus"))),
  #           #             bquote(paste("D) "~italic("Acer")))),
  #                       # bquote(paste("E) "~italic("Quercus"))),
  #                       # bquote(paste("F) "~italic("Ulmus")))), 
  #           label_fontface = "italic",
  #           label_x = 0.3, label_y = 1.01, byrow = FALSE)
      

  p1 <- map_list_1ha_sum[[1]] 
  p2 <- map_list_1ha_sd[[1]]
  p3 <- map_list_1ha_mean_dist[[1]]
  p4 <- map_list_1ha_sum[[2]] 
  p5 <- map_list_1ha_sd[[2]] 
  p6 <- map_list_1ha_mean_dist[[2]]
  p7 <- map_list_1ha_sum[[3]]
  p8 <- map_list_1ha_sd[[3]] 
  p9 <- map_list_1ha_mean_dist[[3]]
  p10 <- map_list_1ha_sum[[4]] 
  p11 <- map_list_1ha_sd[[4]] 
  p12 <- map_list_1ha_mean_dist[[4]]
  

  
  # helper to make a rotated, centered row label
  make_row_label <- function(label, size = 12) {
    wrap_elements(grid::textGrob(
      label,
      rot = 90,
      gp = gpar(fontsize = size, fontface = "italic"),
      hjust = 0.5, vjust = 0.5   # centers text within the grob's allotted space
    )) +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  # helper to make a column label (no rotation needed)
  make_col_label <- function(label, size = 12) {
    wrap_elements(grid::textGrob(
      label,
      gp = gpar(fontsize = size, fontface = "bold"),
      hjust = 0.5, vjust = 0.5
    )) +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  row_label_1 <- make_row_label("Quercus")
  row_label_2 <- make_row_label("Platanus")
  row_label_3 <- make_row_label("Ulmus")
  row_label_4 <- make_row_label("Acer")
  
  col_label_1 <- make_col_label("Mean")
  col_label_2 <- make_col_label("SD")
  col_label_3 <- make_col_label("Mean within focal distance")
  
  # build each row: label + 3 maps, with the label column narrow
  row1 <- row_label_1 + p1  + p2  + p3  + plot_layout(widths = c(0.06, 1, 1, 1))
  row2 <- row_label_2 + p4  + p5  + p6  + plot_layout(widths = c(0.06, 1, 1, 1))
  row3 <- row_label_3 + p7  + p8  + p9  + plot_layout(widths = c(0.06, 1, 1, 1))
  row4 <- row_label_4 + p10 + p11 + p12 + plot_layout(widths = c(0.06, 1, 1, 1))
  
  # top strip of column labels — spacer aligns with row-label column
  col_header <- wrap_elements(grid::nullGrob()) + col_label_1 + col_label_2 + col_label_3 +
    plot_layout(widths = c(0.02, 1, 1, 1))
  
  # stack it all together
  final_plot <- col_header / row1 / row2 / row3 / row4 +
    plot_layout(heights = c(0.02, 1, 1, 1, 1))
  
  final_plot <- final_plot &
    theme(plot.margin = margin(2, 2, 2, 2))  # top, right, bottom, left, in pt
  
  final_plot
  
  ggsave(final_plot, filename = "C:/Users/dsk273/Box/writing/UPPH 25 NYC pollen production/July 2026 submission/fig2_260708.tif",
         units = "in", height = 15, width = 12, dpi = 300, compression = "lzw")
        
### visualizing a single map of pollen exposure (pollen within 400 m x population density) #######################################
 #loop through species
 focal_genus_list <- c("Acer", "Betula", "Gleditsia", "Morus", "Platanus", "Quercus", "Ulmus", "Populus", "Juglans")
 for(i in 1:length(focal_genus_list)){
   focal_genus <-  focal_genus_list[i] 
   #focal_genus <- "Quercus" #Morus Acer Gleditsia Platanus
     
 #load population density raster
  density_raster <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_nyc.tif")
  #plot(density_raster)    
  
 #load tree pollen produced within 400 m
   prod_m_focal_mean <-  rast(
                  paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                         "production_within_400m_", focal_genus, ".tif")) %>% 
     extend(., density_raster) %>% 
     resample(., density_raster, method="bilinear") #resample to get two rasters to match up precisely
  #plot(prod_m_focal_mean)    
  
  #multiply the rasters to get exposure in people x pollen within 400 m
  pop_prod_rast <- density_raster * prod_m_focal_mean
  
  #export exposure raster
  writeRaster(pop_prod_rast, 
              paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                     "exposure_raster_", focal_genus, ".tif"), overwrite = TRUE)
  
  #shrink everything above the 99th percentile of values to 99% for more effective visualization
  p99 <- unlist(global(pop_prod_rast, fun = quantile, probs = 0.99, na.rm = TRUE))
  pop_prod_rast[pop_prod_rast[] > p99] <- p99
  
  # create map
  ggplot() + ggthemes::theme_few() + ggtitle(focal_genus) + 
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    geom_spatraster(data = pop_prod_rast) +
    scale_fill_viridis_c(na.value = "transparent", 
                         option = "turbo",
                         limits = c(0, p99),
                         name = "pollen exposure \n(pollen produced within 400 m x \n population density)",
                         labels = scales::label_comma())
  
 } # end genus loop
      
    
### panel of potential human exposure within 400 m for all tree taxa #############################################################
  #pre load map elements
  #load in nyc boundary polygon
  nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    st_union() %>% #combine the different boroughs
    st_transform(., crs = 32618)
  #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
  # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
  
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  #load population density raster
  density_raster <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_nyc.tif")
  #plot(density_raster)    
  
  
    #create a function that produces a map panel
    fun_expo_map_genus <- function(focal_genus){    
      
      # load raster of pollen production within 400 m for a genus
      prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                                         "production_within_400m_", focal_genus, ".tif")) %>% 
        # prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
        #                                    "production_within_400m_", "Quercus", ".tif")) %>% 
        extend(., density_raster) %>% 
        resample(., density_raster, method="bilinear")
      
      #multiply the rasters to get exposure in people x pollen within 400 m
      pop_prod_rast <- density_raster * prod_m_focal_mean
      
      #convert from billions of pollen grains to trillions
      pop_prod_rast <- (pop_prod_rast/ 1000) + 1
      
      #shrink everything above the 99th percentile of values to 99% for more effective visualization
      p99 <- unlist(global(pop_prod_rast, fun = quantile, probs = 0.99, na.rm = TRUE))
      pop_prod_rast[pop_prod_rast[] > p99] <- p99
      
      # create map
      focal_map_panel <-   
        ggplot() + #ggthemes::theme_few() +   
        geom_spatraster_rgb(data = nyc_topo_spatrast) +
        geom_spatraster(data = pop_prod_rast) +
        scale_fill_viridis_c(na.value = "transparent", 
                             option = "rocket",
                             name = expression(     atop(
                                                    #atop(textstyle("potential"),
                                                         textstyle("exposure"),
                                                    "("~10^18~"grains x people)")),
                             labels = scales::label_comma(),
                             limits = c(0, p99)  #trans = "log10"
                             ) +
        theme(  
          legend.position = c(0.05, 0.95),  # Places the legend at the top-left corner
          legend.justification = c(0.05, 0.98),# Aligns the legend box to its top-left corner)
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
          #plot.background = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
      
      return(focal_map_panel)
    } #end map making function
    
    
    #create list of maps
    focal_genus_panels <- c("Acer", "Gleditsia", "Platanus", "Quercus", "Ulmus") #"Betula", 
    map_panel_list <- map(focal_genus_panels, fun_expo_map_genus)
    
    #create a custom map for density of people 
    #shrink everything above the 99th percentile of values to 99% for more effective visualization
    p99 <- unlist(global(density_raster, fun = quantile, probs = 0.99, na.rm = TRUE))
    density_raster_p99 <- density_raster
    density_raster_p99[density_raster_p99[] > p99] <- p99
    
    # create map
    pop_den_map_panel <- ggplot() +  
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = density_raster_p99) +
      scale_fill_viridis_c(na.value = "transparent", 
                           option = "plasma",
                           name = "people \n(people per ha)",
                           labels = scales::label_comma()) +
      theme(  
        legend.position = c(0.00, 0.95),  # Places the legend at the top-left corner
        legend.justification = c(0.00, 0.98),# Aligns the legend box to its top-left corner)
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white")) #+ ggtitle(focal_genus)
    
    
    #with the population density map as the first panel
    plot_grid(plotlist = c(pop_den_map_panel, map_panel_list), ncol = 2,
              labels = c( "A: population density",
                          bquote(paste("B) "~italic("Acer"))),
                          bquote(paste("C) "~italic("Gleditsia"))),
                          bquote(paste("D) "~italic("Platanus"))),
                          bquote(paste("E) "~italic("Quercus"))),
                          bquote(paste("F) "~italic("Ulmus")))),
              label_x = 0.1, label_y = 1.01)
     
  #version without the first panel being a genus
  # plot_grid(plotlist = map_panel_list, ncol = 2,
  #           labels = c( bquote(paste("A) "~italic("Acer"))),
  #                       bquote(paste("B) "~italic("Betula"))),
  #                       bquote(paste("C) "~italic("Gleditsia"))),
  #                       bquote(paste("D) "~italic("Platanus"))),
  #                       bquote(paste("E) "~italic("Quercus"))),
  #                       bquote(paste("F) "~italic("Ulmus")))), 
  #           label_fontface = "italic",
  #           label_x = 0.3, label_y = 1.01)
  

  
  
### panel of potential asthma exposure within 400 m for all tree taxa #############################################################
  #pre load map elements
  #load in nyc boundary polygon
  nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    st_union() %>% #combine the different boroughs
    st_transform(., crs = 32618)
  #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
  # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
  
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  #load population density raster
  asthma_raster <- rast("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_people_with_asthma.tif")
  #plot(asthma_raster)    
  
  
  
  #create a function that produces a map panel
  fun_asthma_expo_map_genus <- function(focal_genus){    
    
    # load raster of pollen production within 400 m for a genus
    prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
                                       "production_within_400m_", focal_genus, ".tif")) %>% 
      # prod_m_focal_mean <- rast(paste0("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/",
      #                                    "production_within_400m_", "Quercus", ".tif")) %>% 
      extend(., asthma_raster) %>% 
      resample(., asthma_raster, method="bilinear")
    
    #multiply the rasters to get exposure in people x pollen within 400 m
    pop_prod_rast <- asthma_raster * prod_m_focal_mean
    
    #convert from billions of pollen grains to trillions
    pop_prod_rast <- (pop_prod_rast/ 1000) + 1
    
    #shrink everything above the 99th percentile of values to 99% for more effective visualization
    pop_prod_rast_p99 <- pop_prod_rast
    p99 <- unlist(global(pop_prod_rast, fun = quantile, probs = 0.99, na.rm = TRUE))
    pop_prod_rast_p99[pop_prod_rast_p99[] > p99] <- p99
    
    
    # create map
    focal_map_panel <-   
      ggplot() + #ggthemes::theme_few() +   
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = pop_prod_rast_p99) +
      scale_fill_viridis_c(na.value = "transparent", 
                           option = "rocket",
                           name = expression(     atop(
                             #atop(textstyle("potential"),
                             textstyle("exposure"),
                             "("~10^18~"grains x people)")),
                           labels = scales::label_comma()
                           #trans = "log10"
                           #limits = c(10000, max(pop_prod_rast))
      ) +
      theme(  
        legend.position = c(0.00, 0.95),  # Places the legend at the top-left corner
        legend.justification = c(0.00, 0.98),# Aligns the legend box to its top-left corner)
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
        #plot.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
    
    return(focal_map_panel)
  } #end map making function
  
  
  #create list of maps
  focal_genus_panels <- c("Acer", "Gleditsia", "Platanus", "Quercus", "Ulmus")
  map_panel_list <- map(focal_genus_panels, fun_asthma_expo_map_genus)
  
  
  #create a custom map for density of people with asthma 
    #shrink everything above the 99th percentile of values to 99% for more effective visualization
      p99 <- unlist(global(asthma_raster, fun = quantile, probs = 0.99, na.rm = TRUE))
      asthma_raster_p99 <- asthma_raster
      asthma_raster_p99[asthma_raster_p99[] > p99] <- p99
      
    
    # create map
    asthma_map_panel <- ggplot() +  
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = asthma_raster_p99) +
      scale_fill_viridis_c(na.value = "transparent", 
                           option = "inferno",
                           name = "people with asthma \n(people per ha)",
                           labels = scales::label_comma()) +
      theme(  
        legend.position = c(0.00, 0.95),  # Places the legend at the top-left corner
        legend.justification = c(0.00, 0.98),# Aligns the legend box to its top-left corner)
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white")) #+ ggtitle(focal_genus)
  
  
  #with the asthma map as the first panel
  plot_grid(plotlist = c(asthma_map_panel, map_panel_list), ncol = 2,
            labels = c( "A: people with asthma",
                        bquote(paste("B) "~italic("Acer"))),
                        bquote(paste("C) "~italic("Gleditsia"))),
                        bquote(paste("D) "~italic("Platanus"))),
                        bquote(paste("E) "~italic("Quercus"))),
                        bquote(paste("F) "~italic("Ulmus")))),
            label_x = 0, label_y = 1.01)
  
  
  # #without the asthma map as a panel #will need to rerun to include Betula here
  # plot_grid(plotlist = map_panel_list, ncol = 2,
  #           labels = c( bquote(paste("A) "~italic("Acer"))),
  #                       bquote(paste("B) "~italic("Betula"))),
  #                       bquote(paste("C) "~italic("Gleditsia"))),
  #                       bquote(paste("D) "~italic("Platanus"))),
  #                       bquote(paste("E) "~italic("Quercus"))),
  #                       bquote(paste("F) "~italic("Ulmus")))), 
  #           label_fontface = "italic",
  #           label_x = 0.3, label_y = 1.01)
  