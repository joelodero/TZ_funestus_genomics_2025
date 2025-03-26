#-------------------------------------------------------------------------------#
# This script contains R and Rstan code required to reproduce figure 2          #
# and the associated supplementary material for the AsGARD MS                   #
# Modelling spatial and environmental genomic data for An. stephensi            #
#-------------------------------------------------------------------------------#
# Tristan Dennis, 12.2024

# Import libs
libs <- c('tidyverse','data.table','ggthemes','ggquiver','bayesplot','ggeffects','cowplot','maps','terra','sf','geosphere','rstanarm','loo','projpred',
  'tidybayes','corrplot','posterior','lme4','MuMIn','DHARMa','akima','scales','glmmTMB','dggridR')
lapply(libs, library, character.only = TRUE)

# Read metadata
df_samples = read.csv('/Users/dennistpw/Projects/funestus_tz/feems/df_samples.csv')

# Load border data
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")
# Load shapefiles gor plotting
admin_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")
river_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
lake = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster//ne_50m_lakes/ne_50m_lakes.shp")
ocean_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_ocean/ne_50m_ocean.shp")

df_samples_filtered <- df_samples %>% 
  group_by(location_grouped) %>% 
  filter(max(row_number()) > 5) %>%
  ungroup()

#save sample coords


# Region extent
region_bbox <- c(xmin=26, ymin=-25, ymax=4.8, xmax=45)

write.csv(x = df_samples_filtered, file ='/Users/dennistpw/Projects/funestus_tz/feems/df_samples_filtered.csv')


ggplot()+
  geom_sf(data = border_shape, fill = '#f7f7f7',col =gray(0.1))+
  geom_sf(data=st_geometry(lake),colour = '#4a80f5', fill='#9bbff4')+
  coord_sf(xlim = c(26, 45), ylim = c(-25, 4.8), expand = FALSE)+
  geom_point(data=df_samples, aes(x=longitude, y=latitude))
  

# Define some helper variables

# Define the coordinate reference system. Getting this wrong really screws everything up so I've been probably a bit overzealous in specifying the CRS
# Pretty much everywhere where there is the option to - even when many objects in spatial analyses inherit it as a property
crs.geo <- 4326



# Sampling locations
points <- df_samples_filtered %>% dplyr::select(country, location_grouped, latitude, longitude) %>% unique()


# Generate a dggs specifying an intercell spacing of ~150 miles
dggs <- dgconstruct(spacing=75, metric=TRUE, resround='nearest', aperture = 4, topology = "TRIANGLE",  projection = "ISEA")

# Make square box covering sampled region in KSA/HoA/Yemen
region_bbox = sf::st_bbox(region_bbox, crs = crs.geo)
regionshape <- st_as_sfc(region_bbox)

# WO
write.csv(regionshape[[1]][[1]], '/Users/dennistpw/Projects/funestus_tz/feems/riftregionoutline.csv', quote = FALSE, row.names = FALSE, col.names = FALSE)
st_write(regionshape,  '/Users/dennistpw/Projects/funestus_tz/feems/riftregionoutline.shp', append=FALSE)

# Make cropped grid
su_grid <- dgshptogrid(dggs, '/Users/dennistpw/Projects/funestus_tz/feems/riftregionoutline.shp')
st_write(su_grid, '/Users/dennistpw/Projects/funestus_tz/feems/riftregionoutline.TRI.75K.shp', append=FALSE)

plot(su_grid)




# These grids can be used for feems

#######################
# Plot fEEMS output on a map
#######################

# Function for reading feems output
prepare_data <- function(edge_file, node_file){
  edges <- fread(edge_file, col.names = c("from_id", "to_id", "edge_weight"))
  nodes <- fread(node_file, col.names = c("Longitude", "Latitude","N")) %>% mutate(V1 = row_number() - 1)
  
  # Convert necessary columns to integer
  edges$from_id <- as.integer(edges$from_id)
  edges$to_id <- as.integer(edges$to_id)
  nodes$V1 <- as.integer(nodes$V1)
  
  # Join edges and nodes data to get the start and end points of each edge
  edges <- edges %>%
    left_join(nodes, by = c("from_id" = "V1")) %>%
    left_join(nodes, by = c("to_id" = "V1"), suffix = c(".from", ".to")) %>%
    mutate(weight = log10(edge_weight)-mean(log10(edge_weight)))
  
  # Create a list of linestrings, each defined by a pair of points
  edges$geometry <- mapply(function(lon_from, lat_from, lon_to, lat_to) {
    st_linestring(rbind(c(lon_from, lat_from), c(lon_to, lat_to)))
  }, edges$Longitude.from, edges$Latitude.from, edges$Longitude.to, edges$Latitude.to, SIMPLIFY = FALSE)
  
  # Convert edges to an sf object
  edges_sf <- st_as_sf(edges)
  
  # Set the CRS
  st_crs(edges_sf) <- crs.geo
  
  # Convert nodes data.table to an sf object
  nodes_sf <- st_as_sf(nodes, coords = c("Longitude", "Latitude"), crs = crs.geo)
  
  list(edges_sf = edges_sf, nodes_sf = nodes_sf)
}

# Load feems output
edge_file = '/Users/dennistpw/Projects/funestus_tz/feems/rift_incksaedgew.csv'
node_file = '/Users/dennistpw/Projects/funestus_tz/feems/rift_incksanodepos.csv'

edges <- fread(edge_file, col.names = c("from_id", "to_id", "edge_weight"))
nodes <- fread(node_file, col.names = c("Longitude", "Latitude","N")) %>% mutate(V1 = row_number() - 1)

# Convert necessary columns to integer
edges$from_id <- as.integer(edges$from_id)
edges$to_id <- as.integer(edges$to_id)
nodes$V1 <- as.integer(nodes$V1)

# Join edges and nodes data to get the start and end points of each edge
edges <- edges %>%
  left_join(nodes, by = c("from_id" = "V1")) %>%
  left_join(nodes, by = c("to_id" = "V1"), suffix = c(".from", ".to")) %>%
  mutate(weight = log10(edge_weight)-mean(log10(edge_weight)))

# Create a list of linestrings, each defined by a pair of points
edges$geometry <- mapply(function(lon_from, lat_from, lon_to, lat_to) {
  st_linestring(rbind(c(lon_from, lat_from), c(lon_to, lat_to)))
}, edges$Longitude.from, edges$Latitude.from, edges$Longitude.to, edges$Latitude.to, SIMPLIFY = FALSE)

# Convert edges to an sf object
edges_sf <- st_as_sf(edges)

# Set the CRS
st_crs(edges_sf) <- 4326

# Convert nodes data.table to an sf object
nodes_sf <- st_as_sf(nodes, coords = c("Longitude", "Latitude"), crs = 4326)



# Make new lat and long for more concise plotting
df_samples$newlong <- trunc(df_samples$longitude * 100) / 100
df_samples$newlat <- trunc(df_samples$latitude * 100) / 100

# Load collection site status data from WP1a
#site_status <- fread('~/Projects/AsGARD/data/feems_20240920/status_table.csv')

# Plot with feems weighted graph grid on top
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

# Plotting below

# Load ggrastr for good plotting
# Set dpi
options("ggrastr.default.dpi" = 500) 
library(ggrastr)

#set midpoint for plotting
midpoint = median(log10(edges_sf$edge_weight))  # 75th percentile


sample_map <- ggplot()+
  rasterize(geom_sf(data = border_shape, fill = '#f7f7f7',col =gray(0.1))) + 
  geom_sf(data=st_geometry(lake),colour = '#4a80f5', fill='#9bbff4')+
  geom_sf(data=st_geometry(river_shape),colour = '#4a80f5', fill='#9bbff4')+
  geom_sf(data=edges_sf, aes(colour=log10(edges_sf$edge_weight)), alpha=0.7)+
  scale_color_gradient2(low = "#ad3c07",mid='white', high = "#4dc1ff", midpoint=midpoint) +       # Color gradient for weight
  scale_fill_manual(values = c('#808080','#2e8fff','#fcba03'))+
  # coord_sf(xlim=c(28, 52),ylim=c(4,28), expand=FALSE)+
  labs(colour="log10(w)", x='Long.', y='Lat.') +
  #theme_void() +
  coord_sf(xlim = c(26, 45), ylim = c(-25, 4.8), expand = FALSE)+
  geom_point(data=df_samples, aes(x=longitude, y=latitude))+
theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        #legend.position = 'none',
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "top", 
        axis.line = element_blank(),
        axis.text = element_text(size = 10),  # Increase legend text size
        #axis.title = element_blank(),
        legend.text = element_text(size = 10),  # Increase legend text size
        legend.title = element_text(size = 12)) +  
  theme(legend.key=element_blank())

sample_map
ggsave(filename = '~/Projects/AsGARD/figures/feems_map.svg', sample_map, width=7, height=7)


log10(edges_sf$edge_weight)
# Align right-hand plots properly
modelplots <- cowplot::plot_grid(fst_modelplot, pipred, rohpred, 
                                 ncol = 1, 
                                 align = "v", 
                                 axis = "tb",
                                 rel_heights = c(1, 1.15, 1)) # Ensures equal height
ggsave(filename = '~/Projects/AsGARD/figures/modelplots.svg', modelplots, width=4, height=7)

# Fix the overall layout
final_plot <- cowplot::plot_grid(sample_map, modelplots, 
                                 ncol = 2, 
                                 rel_widths = c(2, 1),  # Adjust widths if needed
                                 align = "hv",  # Ensures vertical & horizontal alignment
                                 axis = "tblr")


print(final_plot)



r.squaredGLMM(mod_pi_int_gam)

r.squaredGLMM(mod_roh_int_gaus)


