# Project: Pollen production and exposure in NYC
# Class: PLSCI 4450/6450 Urban Plants and Public Health
# Contact: Dan Katz, Cornell University
# 
# This script is for downloading population density from the US census 
# and for calculating the number of people that live within a certain distance of each point


#set up work environment
library(tidycensus)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(ggspatial) 

#census_api_key("INSERT YOUR CENSUS KEY HERE", install = TRUE, overwrite = TRUE)

### download census population data and rasterize ################################################

# Get population data at tract level for NYC
nyc_blocks <- get_decennial(
      geography = "block",
      variables = "P1_001N",  # Total population
      state = "NY",
      county = c("New York", "Kings", "Queens", "Bronx", "Richmond"),
      year = 2020,
      geometry = TRUE
    ) 
#head(nyc_blocks)

# Calculate area in square km and convert to population density and transform to UTM 18N
nyc_block_density <- nyc_blocks %>%
  rename(population = value) %>%
  mutate(
     area_sqkm = as.numeric(st_area(geometry)) / 1e6,         # Convert m² to sq km
     pop_density_sqkm = ifelse(area_sqkm > 0, population / area_sqkm, 0),
     pop_density_ha = pop_density_sqkm / 100 # Convert to population density per ha
  ) %>% 
  st_transform(., crs = 32618) #convert to EPSG 32618 for UTM 18N
  
#head(nyc_block_density)
#gut check: does this match the size of NYC? 
#sum(nyc_block_density$area_sqkm) #yes, compared to ~1,220 km2
#plot(nyc_block_density)
#sum(nyc_block_density$population)



### convert population density from polygons to raster ################################################

# Get extent of NYC
bbox <- st_bbox(nyc_block_density)

# Create empty raster template with 100m resolution
raster_template <- rast(
  xmin = bbox["xmin"],
  xmax = bbox["xmax"],
  ymin = bbox["ymin"],
  ymax = bbox["ymax"],
  resolution = 100,  # 100 meters
  crs = st_crs(nyc_block_density)$wkt
)

# Convert sf to SpatVector for terra
nyc_vect <- vect(nyc_block_density)

# Rasterize using population density (per ha)
density_raster <- rasterize(
  nyc_vect, 
  raster_template, 
  field = "pop_density_ha",
  fun = "mean"  # Use mean if multiple blocks overlap a cell
)

  writeRaster(density_raster, "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_nyc.tif", overwrite = TRUE)

#plot(density_raster)

# Create a more detailed ggplot visualization
# Convert raster to dataframe for ggplot
density_df <- as.data.frame(density_raster, xy = TRUE) %>%
  filter(!is.na(pop_density_ha)) %>% 
  filter(pop_density_ha > 0)

ggplot() +
  geom_raster(data = density_df, aes(x = x, y = y, fill = pop_density_ha + 1)) +
  scale_fill_viridis_c(
    name = "Population\nDensity\n(per ha)",
    trans = "log10",
    labels = scales::comma,
    #option = "plasma",
    na.value = "transparent"  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank() ) +
  labs( title = "NYC Population Density - 100m Resolution",
    subtitle = "2020 Decennial Census",
    x = NULL,
    y = NULL)





#load in nyc boundary polygon

nyc_topo_rast <- basemap_raster(bbox, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 

density_raster_no_zeros <- density_raster
density_raster_no_zeros[density_raster_no_zeros == 0] <- NA

# create map
ggplot() + ggthemes::theme_few() +   
  geom_spatraster_rgb(data = nyc_topo_spatrast) +
  geom_spatraster(data = density_raster, alpha = 0.7) +
  scale_fill_viridis_c(na.value = "transparent", 
                       #option = "magma",
                       name = "population density \n(residents per ha)",
                       labels = scales::label_comma())


### asthma prevalence based on CDC Places dataset ################################
  #download from the CDC Places website after filtering to New York at the Census Tract level
  # https://data.cdc.gov/500-Cities-Places/PLACES-Census-Tract-Data-GIS-Friendly-Format-2024-/yjkw-uj5s/about_data
  tract_asthma_places <- read_csv("C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/PLACES__Census_Tract_Data_asthma.csv") %>% 
    mutate(TractFIPS = as.character(TractFIPS))
  
  # Get population data at tract level for NYC
      nyc_tracts <- get_decennial(
        geography = "tract",
        variables = "P1_001N",  # Total population
        state = "NY",
        county = c("New York", "Kings", "Queens", "Bronx", "Richmond"),
        year = 2020,
        geometry = TRUE
      ) 
      nyc_tracts <- nyc_tracts %>%
        mutate(TractFIPS = as.character(GEOID))
  
      plot(nyc_tracts)

  #join asthma and census tract data and calculate area
   nyc_mean_asthma <- left_join(nyc_tracts, tract_asthma_places) %>% st_drop_geometry() %>% 
     summarize(CASTHMA_CrudePrev_mean_nyc = mean(CASTHMA_CrudePrev, na.rm=TRUE))
      
      
  tract_asthma <- left_join(nyc_tracts, tract_asthma_places) %>%
    st_transform(., crs = 32618) %>% #convert to EPSG 32618 for UTM 18N
    mutate(
      population_census = value, #use the directly downloaded population data; the Places data assumes 0 people when it's low
      area_ha = as.numeric(st_area(geometry)) / 10000,         # Convert m² to sq km
      pop_density_ha = ifelse(area_ha > 0, population_census / area_ha, 0),
      CASTHMA_CrudePrev_interp = case_when(!is.na(CASTHMA_CrudePrev) ~ CASTHMA_CrudePrev,
                                           is.na(CASTHMA_CrudePrev) ~ nyc_mean_asthma$CASTHMA_CrudePrev_mean_nyc[1]), #assume missing tracts have same crude rate
                                                                                        #as NYC in general
      people_with_asthma_ha = pop_density_ha * (CASTHMA_CrudePrev_interp/100)
    )  #plot(tract_asthma)
       #tract_asthma %>% dplyr::select(people_with_asthma_ha) %>% plot()
  
  #how many people was this value imputed for
  sum(tract_asthma$value[is.na(tract_asthma$CASTHMA_CrudePrev) & tract_asthma$value > 0])
  
  summary(tract_asthma)
  # Get extent of NYC
  bbox <- st_bbox(nyc_block_density)
  
  # Create empty raster template with 100m resolution
  raster_template <- rast(
    xmin = bbox["xmin"],
    xmax = bbox["xmax"],
    ymin = bbox["ymin"],
    ymax = bbox["ymax"],
    resolution = 100,  # 100 meters
    crs = st_crs(nyc_block_density)$wkt
  )
  
  # Convert sf to SpatVector for terra
  asthma_vect <- vect(tract_asthma)
  
  # Rasterize using population density (per ha)
  asthma_raster <- rasterize(
    asthma_vect, 
    raster_template, 
    field = "people_with_asthma_ha",
    fun = "mean"  # Use mean if multiple blocks overlap a cell
  )
  #plot(asthma_raster)
  
  writeRaster(asthma_raster, "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_people_with_asthma.tif", overwrite = TRUE)
  
  # Create a more detailed ggplot visualization
  # Convert raster to dataframe for ggplot
 asthma_df <- as.data.frame(asthma_raster, xy = TRUE) %>%
    filter(!is.na(people_with_asthma_ha))
  
  ggplot() +
    geom_raster(data = asthma_df, aes(x = x, y = y, fill = people_with_asthma_ha + 1)) +
    scale_fill_viridis_c(
      name = "Population\nDensity\n(per ha)",
      trans = "log10",
      labels = scales::comma,
      #option = "plasma",
      na.value = "transparent"  ) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank() ) +
    labs( title = "people with asthma - 100m Resolution",
          subtitle = "2020 Decennial Census",
          x = NULL,
          y = NULL)

  #load in nyc boundary polygon
  nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
    st_union() %>% #combine the different boroughs
    st_transform(., crs = 32618)
  
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  # create map
  ggplot() + ggthemes::theme_few() +   
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    geom_spatraster(data = asthma_raster, alpha = 0.6) +
    scale_fill_viridis_c(na.value = "transparent", 
                         option = "magma",
                         name = "people with asthma \n(people per ha)",
                         labels = scales::label_comma()) +
    annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
                     bar_cols = c("black", "white"), # Colors of the scale bar segments
                     style = "ticks",
                     text_cex = 0.8) +  # Text size for the scale bar label
    annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
                           pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
    theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
            legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)
  
  

### sum population within 400 m of a cell #######################################

  # Create a circular focal window
  focal_matrix <- focalMat(density_raster, d = 400, type = "circle", fillNA = TRUE)
  focal_matrix_no_weights <- focal_matrix

  # Replace all values > 0 with 1 to create an unweighted window
  focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1

  density_focal_sum <-
    focal(
      density_raster,
      w = focal_matrix_no_weights,
      fun = "sum",
      na.rm = TRUE)

  names(density_focal_sum) <- "pop_within_400_m"

  plot(density_focal_sum)


  # visualize the data
  focal_df <- as.data.frame(density_focal_sum, xy = TRUE) %>%
     filter(!is.na(pop_within_400_m))

  ggplot() +
    geom_raster(data = focal_df, aes(x = x, y = y, fill = pop_within_400_m)) +
    scale_fill_viridis_c(
      name = "residents within 1 km",
      #trans = "log10",
      labels = scales::comma,
      option = "plasma"
    ) +
    coord_equal() +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())

  # Save the focal sum raster
  writeRaster(density_focal_sum, "C:/Users/dsk273/Box/classes/plants and public health fall 2025/class project analysis/nyc_pop_density_400m_focal_sum.tif", overwrite = TRUE)

