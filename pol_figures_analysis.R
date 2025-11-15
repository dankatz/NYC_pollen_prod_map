## this script is for analyzing and visualizing pollen production in NYC and for detecting the influence of individual trees

## set up work environment

## load in tree pollen production estimates

## load in census population density
density_focal_sum


## extract population density for each tree 
tr_export_centroids_proj

polpop <- tr_export_centroids_proj %>% 
          mutate(pop_within_1_km = terra::extract(density_focal_sum, tr_export_centroids_proj)[,2],
                 pop_pol = pol_mean * pop_within_1_km) # calculate pollen production * population to get potential impact

polpop_df <- polpop %>% st_drop_geometry(.)

## visualize results
ggplot(polpop_df, aes(x = pop_pol)) + geom_histogram() + facet_wrap(~Genus, scales = "free") + theme_bw() + scale_x_log10()

#cdf plot
ggplot(polpop_df, aes(x = pop_pol + 1, color = Genus)) + stat_ecdf(geom = "step", pad = FALSE) + theme_bw() + scale_x_log10() +
  xlab("exposure from tree (pollen grains x population within 1 km)") + ylab("cumulative distribution (proportion)")

#map of a particular genus
polpop_sub <- polpop %>% 
  filter(Genus == "Quercus") %>% 
  slice_sample(n = 20000) %>% 
  filter(tree_area > 125 & tree_area < 150)

#map of pollen production by a particular genus
polpop_sub <- polpop %>% 
  filter(Genus == "Quercus") %>% 
  slice_sample(n = 10000) # %>%  filter(tree_area > 125 & tree_area < 150)

ggplot() + 
  geom_sf(data = polpop_sub, fill = "lightgray") + # Draw the base map
  geom_sf(data = polpop_sub, aes(color = pol_mean), size = 1, alpha = 0.1) + # Add points
  theme_minimal() + 
  scale_color_viridis_c(name = "pollen produced by tree x population within 1 km")# Apply a minimal theme


#comparison of size to pop_pol
ggplot(polpop, aes(x = tree_area, y = pop_pol)) + geom_hex(name = "density") + facet_wrap(~Genus, scales = "free") + theme_bw() +
  scale_fill_viridis_c(trans = "log10") + xlab("tree area (m2)") + ylab("per-tree exposure (people within 1 km x pollen)")
       