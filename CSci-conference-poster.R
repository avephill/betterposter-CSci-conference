library(sf)
library(terra)
library(tidyverse)


results.dir <- "../results/sdm/run6_addHI/"



# Methods -----------------------------------------------------------------


inat_cleaned.sf <- st_read(paste0(results.dir, "inat_obs.gpkg"))
species <- inat_cleaned.sf$species %>% unique()

length(species)
nrow(inat_cleaned.sf)

inat_cleaned.sf %>% slice_min(year, n = 10) %>% pull(year)
hist(inat_cleaned.sf$year)

predictor.stack <- rast(paste0(results.dir, "predictors.tif"))
predictor.stack %>% names() %>%  length()



# Main Map ----------------------------------------------------------------

library(basemaps)

study.sf <- st_read("../data/study_area.gpkg") %>% st_transform(crs= 3857)
urban.sf <- st_read("../data/urban_area.gpkg") %>% st_transform(crs= 3857)
inat_cleaned.sf <- st_read(paste0(results.dir, "inat_obs.gpkg"))
terrestrial.sf <- st_read("~/Data/Boundaries/Political/CLN2.0/CLN_Database_2_0_1/data/Base/vector/terrestrial.shp")
# Bay counties
counties.sf <- st_read("../../../Data/Boundaries/Political/CA_Counties/CA_Counties_TIGER2016.shp") %>% 
    filter(NAME %in% c("Solano", "Contra Costa", "Alameda", "Santa Clara", "San Mateo",
                       "San Francisco", "Marin", "Sonoma", "Napa", "Santa Cruz")) %>% 
    st_transform(crs = 4326)
target_counties.sf <- 
    counties.sf %>% 
    st_union() 

inat24.sf <- inat_cleaned.sf


us.sf <- st_read("../../../Data/Boundaries/Political/USA/USA.shp") %>%
    st_transform(crs=3857)

background_landmask.sf <- study.sf %>%
    # Remove the ocean
    st_intersection(terrestrial.sf %>% st_transform(3857)) %>%
    dplyr::select(geom) %>% 
    # limit to just counties
    st_intersection(target_counties.sf %>% 
                        st_transform(crs = 3857) %>% 
                        st_intersection(study.sf))

cln_habitat.sf <- st_read("../../../Data/Boundaries/Political/CLN2.0/CLN_Database_2_0_1/data/CLN/final_cln2_2019_dissolved.shp")
cln_habitat.sf %>% pull(Symbol) %>% unique()
# cln_cons.sf <- cln_habitat.sf %>% 
#     filter(Symbol %in% c("Essential","Contributor","Important","Connector")) %>% 
#     st_union() %>% 
#     st_transform(crs= 3857)
# 
# write_sf(cln_cons.sf, "../data/cln_habitat.gpkg")
cln_cons.sf <- st_read("../data/cln_habitat.gpkg")

area_ext <- study.sf %>% st_bbox()
get_maptypes()

# basemap_magick(area_ext, map_service = "osm_stamen", map_type = "toner")

bay_base.g <- 
    ggplot() + 
    # basemap_gglayer(study.sf, map_service = "osm_stamen", map_type = "toner", map_res = 1) +
    geom_sf(data = background_landmask.sf, fill = NA) +
    theme_void() +
    coord_sf() +
    scale_fill_identity()

# ggsave("../results/figures/baybase.png",bay_base.g, dpi = 300, width = 8, height = 8)

# bay_base.g + 
#     geom_sf(data = cln_cons.sf, fill="darkgoldenrod1", alpha = .5)

# ggsave("../results/figures/baybase2.png",dpi = 300, width = 8, height = 8)

bay_base.g + 
    geom_sf(data = cln_cons.sf, aes(fill="darkgoldenrod1"), color = NA, alpha = .7) +
    geom_sf(data = urban.sf, aes(fill = "steelblue"), color = "steelblue", alpha = .3) +
    scale_fill_identity(breaks = c("darkgoldenrod1", "steelblue"), 
                        labels = c("Regions of \nConservation Value", "Urban Lands"), 
                        guide = "legend")+
    theme(legend.text = element_text(face = "bold", size = 32),
          legend.position = c(.25,.25),
          panel.border = element_rect(colour = "white", fill=NA, linewidth=5),
          legend.title = element_blank())
    

ggsave("img/urbancln.png", dpi = 300, width = 8, height = 8)
# ggsave("img/urbancln.eps", dpi = 300, width = 8, height = 8, device = cairo_ps)


full.plot <- bay_base.g + 
    geom_sf(data = cln_cons.sf, fill="darkgoldenrod1") +
    geom_sf(data = urban.sf, fill = "steelblue", color = NA) #+
    # geom_sf(data = inat24.sf %>% filter(coordinateuncertaintyinmeters < 500))

ggsave("../results/figures/baybase4.png",
       full.plot,
       dpi = 300, width = 8, height = 8)
