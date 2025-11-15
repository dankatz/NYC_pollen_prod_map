# Project: UPPH '25 map of pollen production in NYC
# For any questions or comments, contact Dan Katz, Cornell University


#set up work environment
library(tidycensus)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)

census_api_key("INSERT YOUR CENSUS KEY HERE", install = TRUE, overwrite = TRUE)

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
  
head(nyc_block_density)
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

# Rasterize using population density (per sq mile)
density_raster <- rasterize(
  nyc_vect, 
  raster_template, 
  field = "pop_density_ha",
  fun = "mean"  # Use mean if multiple blocks overlap a cell
)

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


### sum population within 1 km of a cell #######################################

# Create a circular focal window
focal_matrix <- focalMat(density_raster, d = 1000, type = "circle", fillNA = TRUE)
focal_matrix_no_weights <- focal_matrix

# Replace all values > 0 with 1 to create an unweighted window
focal_matrix_no_weights[focal_matrix_no_weights > 0] <- 1

density_focal_sum <- 
  focal(
    density_raster,
    w = focal_matrix_no_weights,
    fun = "sum",
    na.rm = TRUE)

names(density_focal_sum) <- "pop_within_1_km"

plot(density_focal_sum)


# visualize the data
focal_df <- as.data.frame(density_focal_sum, xy = TRUE) %>%
   filter(!is.na(pop_within_1_km))

ggplot() +
  geom_raster(data = focal_df, aes(x = x, y = y, fill = pop_within_1_km)) +
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
writeRaster(density_focal_sum, "nyc_pop_density_1km_focal_sum.tif", overwrite = TRUE)

