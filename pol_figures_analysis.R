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

### load in tree pollen production estimates
# these were created in the 'pollen_prod_estimates.R' script
tr_export_centroids_proj_full_csv <- read_csv( "C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/NYC_pollen_prod_estimates_251117.csv")
tr_export_centroids_proj <- st_as_sf(tr_export_centroids_proj_full_csv, coords = c("x_EPSG_32618" , "y_EPSG_32618"),
                                          crs = 32618)

## load in census population density
# these were created in the 'census_map.R' script
density_focal_sum <- rast("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_400m_focal_sum.tif")

## extract population density for each tree 
tr_export_centroids_proj

polpop <- tr_export_centroids_proj %>% 
          mutate(pop_within_1_km = terra::extract(density_focal_sum, tr_export_centroids_proj)[,2],
                 pop_pol = pol_mean * pop_within_1_km) # calculate pollen production * population to get potential impact

polpop_df <- polpop %>% st_drop_geometry(.)


### visualize results ###################################################

#Histogram of pollen exposure for each tree by genus
ggplot(polpop_df, aes(x = pop_pol)) + geom_histogram() + facet_wrap(~Genus, scales = "free") + theme_bw() + scale_x_log10()

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
       

# map of tree pollen within 400 m for a single genus
    prod_400m_focal_sum <- rast(paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
                       "production_within_400m_", "Quercus", ".tif"))
    plot(prod_400m_focal_sum)
    hist(prod_400m_focal_sum)
    
    #load in nyc boundary polygon
    nyc_boundary <- st_read( "C:/Users/danka/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
      st_union() %>% #combine the different boroughs
      st_transform(., crs = 32618)
        #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
        # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
    
    nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
    nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 

    # create map
      ggplot() + ggthemes::theme_few() +   
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = prod_400m_focal_sum/1000, alpha = 0.6) +
      scale_fill_viridis_c(na.value = "transparent", 
                           #option = "magma",
                           name = "pollen production \n(trillions of grains within 400 m)",
                           labels = scales::label_comma()) +
        annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
                         bar_cols = c("black", "white"), # Colors of the scale bar segments
                         style = "ticks",
                         text_cex = 0.8) +  # Text size for the scale bar label
        annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
                               pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
        theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
                legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)


        
        
        
### panel of pollen production within 400 m for all tree taxa #############################################################
 #pre load map elements
    #load in nyc boundary polygon
      nyc_boundary <- st_read( "C:/Users/danka/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
        st_union() %>% #combine the different boroughs
        st_transform(., crs = 32618)
      #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
      # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
      
      nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
      nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
      
      
#create a function that produces a map panel
  fun_prod_map_genus <- function(focal_genus){    
      
      # load raster of pollen production within 400 m for a genus
      prod_400m_focal_sum <- rast(paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
                                         "production_within_400m_", focal_genus, ".tif"))
 
      # create map
    focal_map_panel <-   
      ggplot() + #ggthemes::theme_few() +   
        geom_spatraster_rgb(data = nyc_topo_spatrast) +
        geom_spatraster(data = prod_400m_focal_sum/1000, alpha = 0.6) +
        scale_fill_viridis_c(na.value = "transparent", 
                             #option = "magma",
                             name = expression(atop(atop(textstyle("pollen"),
                                                    textstyle("production")),
                                               "("~10^12~"grains)")),
                             labels = scales::label_comma()) +
              theme(  
                legend.position = c(0.05, 0.85),  # Places the legend at the top-left corner
                legend.justification = c(0.05, 0.95),# Aligns the legend box to its top-left corner)
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 8),
                axis.title=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(), 
                plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
                #plot.background = element_blank(),
                plot.background = element_rect(fill = "white", color = "white"),
                panel.background = element_rect(fill = "white"),
                plot.title = element_text(face = "italic")) #+ ggtitle(focal_genus)
    
    return(focal_map_panel)
  } #end map making function
  
  
#create list of maps
  focal_genus_panels <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus")
  test <- map(focal_genus_panels, fun_prod_map_genus)
   
  plot_grid(plotlist = test, ncol = 2,
            labels = c( bquote(paste("A) "~italic("Acer"))),
                        bquote(paste("B) "~italic("Betula"))),
                        bquote(paste("C) "~italic("Gleditsia"))),
                        bquote(paste("D) "~italic("Platanus"))),
                        bquote(paste("E) "~italic("Quercus"))),
                        bquote(paste("F) "~italic("Ulmus")))), 
            label_fontface = "italic",
            label_x = 0.3, label_y = 1.01)
      
      
### visualizing a map of pollen exposure (pollen within 400 m x population density) #######################################
 #loop through species
 focal_genus_list <- c("Acer", "Betula", "Gleditsia", "Morus", "Platanus", "Quercus", "Ulmus", "Populus", "Juglans")
 for(i in 1:length(focal_genus_list)){
   focal_genus <-  focal_genus_list[i] 
   #focal_genus <- "Juglans" #Morus Acer Gleditsia Platanus
     
 #load population density raster
  density_raster <- rast("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_nyc.tif")
  #plot(density_raster)    
  
 #load tree pollen produced within 400 m
   prod_400m_focal_sum <-  rast(
                  paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
                         "production_within_400m_", focal_genus, ".tif")) %>% 
     extend(., density_raster) %>% 
     resample(., density_raster, method="bilinear") #resample to get two rasters to match up precisely
  #plot(prod_400m_focal_sum)    
  
  #multiply the rasters to get exposure in people x pollen within 400 m
  pop_prod_rast <- density_raster * prod_400m_focal_sum
  
  #export exposure raster
  writeRaster(pop_prod_rast, 
              paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
                     "exposure_raster_", focal_genus, ".tif"), overwrite = TRUE)
  
  # create map
  ggplot() + ggthemes::theme_few() + ggtitle(focal_genus) + 
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    geom_spatraster(data = pop_prod_rast, alpha = 0.6) +
    scale_fill_viridis_c(na.value = "transparent", 
                         option = "plasma",
                         name = "pollen exposure \n(pollen produced within 400 m x \n population density)",
                         labels = scales::label_comma())
  
 } # end genus loop
      
    
### panel of potential human exposure within 400 m for all tree taxa #############################################################
  #pre load map elements
  #load in nyc boundary polygon
  nyc_boundary <- st_read( "C:/Users/danka/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    st_union() %>% #combine the different boroughs
    st_transform(., crs = 32618)
  #nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
  # nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)
  
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  #load population density raster
  density_raster <- rast("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_nyc.tif")
  #plot(density_raster)    
  

  
  #create a function that produces a map panel
  fun_expo_map_genus <- function(focal_genus){    
    
    # load raster of pollen production within 400 m for a genus
    prod_400m_focal_sum <- rast(paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
                                       "production_within_400m_", focal_genus, ".tif")) %>% 
      # prod_400m_focal_sum <- rast(paste0("C:/Users/danka/Box/classes/plants and public health fall 2025/class project analysis/",
      #                                    "production_within_400m_", "Quercus", ".tif")) %>% 
      extend(., density_raster) %>% 
      resample(., density_raster, method="bilinear")
    
    #multiply the rasters to get exposure in people x pollen within 400 m
    pop_prod_rast <- density_raster * prod_400m_focal_sum
    
    #convert from billions of pollen grains to trillions
    pop_prod_rast <- pop_prod_rast/ 1000
    
    # create map
    focal_map_panel <-   
      ggplot() + #ggthemes::theme_few() +   
      geom_spatraster_rgb(data = nyc_topo_spatrast) +
      geom_spatraster(data = pop_prod_rast, alpha = 0.8) +
      scale_fill_viridis_c(na.value = "transparent", 
                           option = "turbo",
                           name = expression(     atop(
                                                  #atop(textstyle("potential"),
                                                       textstyle("exposure"),
                                                  "("~10^18~"grains x people)")),
                           labels = scales::label_comma()) +
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
  focal_genus_panels <- c("Acer", "Betula", "Gleditsia", "Platanus", "Quercus", "Ulmus")
  map_panel_list <- map(focal_genus_panels, fun_expo_map_genus)
  
  plot_grid(plotlist = map_panel_list, ncol = 2,
            labels = c( bquote(paste("A) "~italic("Acer"))),
                        bquote(paste("B) "~italic("Betula"))),
                        bquote(paste("C) "~italic("Gleditsia"))),
                        bquote(paste("D) "~italic("Platanus"))),
                        bquote(paste("E) "~italic("Quercus"))),
                        bquote(paste("F) "~italic("Ulmus")))), 
            label_fontface = "italic",
            label_x = 0.3, label_y = 1.01)
  
