#' Produce materials and figures for the CUB CSci conference.
#' Also included is calculation of MaxTSS thresholds and binary thresholding
#' of SDMs

library(sf)
library(terra)
library(tidyverse)
library(gt)
library(stars)


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

# library(basemaps)

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

text_color <- "white"

cln_color <- "#F9CB40" #"#B1D8B7"
urb_color <- "#6C91C2"

urbancln.g <- 
    bay_base.g + 
    geom_sf(data = cln_cons.sf, aes(fill=cln_color), color = NA, alpha = 1) +
    geom_sf(data = urban.sf, aes(fill = urb_color), color = urb_color, alpha = 1) +
    scale_fill_identity(breaks = c(cln_color, urb_color), 
                        labels = c("Recognized\nConservation\nValue", "Urban Lands"), 
                        guide = "legend")+
    theme(legend.text = element_text(face = "bold", size = 28, color = text_color),
          legend.position = c(.2,.25),
          legend.spacing.y = unit(.25, 'in'),
          panel.border = element_rect(colour = "white", fill=NA, linewidth=5),
          legend.title = element_blank()#
          # plot.margin = unit(-c(2,2,2,2), "in"),
          # panel.spacing = unit(c(4,4,4,4), "in"),
          # panel.spacing.x = unit(4, "in")
    ) +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE))


ggsave("img/urbancln.png", urbancln.g, dpi = 300, width = 10, height = 10)
# ggsave("img/urbancln.eps", dpi = 300, width = 8, height = 8, device = cairo_ps)


# full.plot <- bay_base.g + 
#     geom_sf(data = cln_cons.sf, fill="darkgoldenrod1") +
#     geom_sf(data = urban.sf, fill = "steelblue", color = NA) #+
# # geom_sf(data = inat24.sf %>% filter(coordinateuncertaintyinmeters < 500))
# 
# ggsave("../results/figures/baybase4.png",
#        full.plot,
#        dpi = 300, width = 8, height = 8)




# Main Plot ---------------------------------------------------------------
library(SDMtune)
library(maxnet)
library(furrr)

# Get SDM fps
results.dir <- "../results/sdm/run6_addHI/"
spec.dirs <- list.dirs(results.dir, recursive = T, full.names = T)[-1]

# Just get dirs with all the sdm stuff in it
full_dirs <- 
    map(spec.dirs, 
        function(x){
            fnames <- list.files(x)
            targetfiles <- fnames %>% str_subset("prediction.tif|sdm_model.rds|var_imp.csv")
            if(length(targetfiles) == 3){
                return(x)
            }
        }) %>% 
    unlist()



# Get thresholds for each model
plan(multisession, workers = 5)
threshes.df <- 
    future_map(full_dirs, 
               function(fp){
                   library(maxnet)
                   model <- readRDS(paste0(fp, "/sdm_model.rds"))
                   cv.folds <- model@models %>% length()
                   spec_name <- model@models[[1]]@data@species
                   
                   # Get maxTSS thresh from each fold and just average
                   avg_thresh <- map(1:cv.folds, function(x){
                       threshes <- SDMtune::thresholds(model@models[[x]], type = "cloglog")
                       threshes %>% 
                           filter(Threshold == "Maximum training sensitivity plus specificity") %>% 
                           select(`Cloglog value`)
                   }) %>% unlist() %>% mean()
                   
                   df <- tibble(
                       Species = spec_name,
                       `MaxTSS Cloglog Thresh` = avg_thresh)
                   
                   return(df)
               }) %>% bind_rows()



# Make all previous predictions binary and then put them in a list
binary_predictions <- 
    future_map(full_dirs, 
               function(fp){
                   prediction.sr <- rast(paste0(fp, "/prediction.tif"))
                   
                   # Get spec name
                   fp_broken <- ( fp %>% str_split_1("/") )
                   spec_name <- fp_broken[length(fp_broken)] %>% str_replace("_", " ")
                   
                   # Get thresh
                   thrsh <- threshes.df %>% 
                       filter(Species == spec_name) %>% 
                       pull(`MaxTSS Cloglog Thresh`)
                   
                   reclass.m <- 
                       c(-Inf, thrsh, 0,
                         thrsh, Inf, 1) %>% 
                       matrix(ncol = 3, byrow = T)
                   
                   prediction_bin.sr <- classify(prediction.sr, reclass.m)
                   names(prediction_bin.sr) <- spec_name
                   
                   return(wrap(prediction_bin.sr))
               })

# Concatenate species predictions into a raster
binary_pred.sr <- rast()
for(pred in binary_predictions){
    binary_pred.sr <- c(binary_pred.sr, unwrap(pred))
}

# Add up observations
spec_rich.sr <- app(binary_pred.sr, sum, cores = 5)
names(spec_rich.sr) <- "Species Richness"

writeRaster(spec_rich.sr, 
            "../results/sdm/run6_addHI/species_richness.tif",
            overwrite = T)



# Now actually plot it
specrich.strs <- 
    read_stars("../results/sdm/run6_addHI/species_richness.tif") %>% 
    st_downsample(10) %>% 
    st_transform(3857)



specrich.sf <- specrich.strs %>% st_as_sf()

specrich.g <- ggplot() +
    # geom_stars(data = specrich.strs, aes(fill = species_richness.tif)) +
    geom_sf(data = specrich.sf, aes(fill = species_richness.tif), 
            color = NA, lwd = 0) +
    # geom_sf(data = urban.sf, aes(color = urb_color), fill = NA, lwd = 2) +
    scale_fill_gradient2(low = "#FFEAEE",
                         mid = "#828F7E",
                         high = "#2F5233",
                         midpoint = 10,
                         aesthetics = "fill",
                         name = "Estimated Count\nof Species with\nConservation Value") +
    # scale_color_identity(breaks = c(urb_color),
    #                      labels = c("Urban Lands"),
    #                      guide = "legend",
    #                      name = "") +
    theme_void() +
    theme(legend.position=c(.2, .3),
          legend.text = element_text(face = "bold", size = 28, color = text_color),
          legend.title = element_text(size = 32, face = "bold", color = text_color),
          legend.key.size = unit(3, "line")#,
          # legend.spacing.y = unit(2, "cm"),
          # legend.margin = margin(1,1,1,1, "cm")
          ) +
    guides(fill = guide_legend(byrow = TRUE))


# wid <- 7
# hei <- 10
# dev.new(width=wid, height=hei, unit="in")
# print(specrich.g)
# dev.off()

asp_rat <- .947
height <- 20

ggsave("img/species_richness.png", specrich.g,
       width = height * asp_rat, height = height, dpi = 300, units = "in")

ggsave("img/species_richness.eps", #width = 7.2,
       specrich.g,
       height = height,
       width = height * asp_rat,
       units = "in",
       device = cairo_ps)


# system('convert /home/ahill/Projects/together-bay-area/CUB-poster/img/species_richness.png -trim /home/ahill/Projects/together-bay-area/CUB-poster/img/species_richness.png')
# knitr::plot_crop("img/species_richness.png")

# Tables ------------------------------------------------------------------
# library(taxize)
# Species Table

sum.df <- read_csv("../together-bay-area/SDM-Report-1/data/sdm/sdm_sum.csv") %>% 
    mutate(AUC = round(AUC, 2),
           TSS = round(TSS, 2),
           `Most Imp. Variables` = `Most Imp. Variables` %>% str_replace_all("_", " "))

# uids <- get_uid(sum.df$Species)
# common <- sci2comm(uids, db = "ncbi")

# make_hyperlink <- function(spec) {
#     paste0('<a href="https://www.inaturalist.org/taxa/',
#            str_replace(spec, " ", "-"),'" target="_blank">', spec, "</a>")
# }

spec_tab.gt <- 
    sum.df %>% 
    select(-`Most Imp. Variables`) %>%
    gt(rowname_col = "Species") %>% 
    # tab_header(title = "Quick SDM metrics by species") %>%
    tab_style(
        style = cell_text(style = "italic"),
        locations = cells_stub()
    ) %>% 
    tab_footnote(footnote = "Each species name is hyperlinked to an iNaturalist page (some links may be broken)",
                 locations = cells_column_labels(
                     columns = Species
                 ))


gtsave(spec_tab.gt,
       "img/species_table.tex")
# gtsave(spec_tab.gt,
#        "img/species_table.png", vwidth = 1500, vheight = 1200)


# Predictors --------------------------------------------------------------

spec.fp <- species.fp %>%
    str_subset(input$species_name_map)

varimp_tot.df <- 
    map(full_dirs, function(fp){
        
        # Get spec name
        fp_broken <- ( fp %>% str_split_1("/") )
        spec_name <- fp_broken[length(fp_broken)] %>% str_replace("_", " ")
        
        # Read csv with Variable importances
        varimp.df <- read_csv(paste0(fp, "/var_imp.csv"))
        varimp_piece.df <- 
            varimp.df %>% 
            mutate(Species = spec_name)
        return(varimp_piece.df)
    }) %>% bind_rows()

varimp_processed.df <- 
    varimp_tot.df %>% 
    select(-sd) %>% 
    mutate(Variable = Variable %>% str_replace_all("_", " ")) %>% 
    group_by(Variable) %>% 
    summarise(mean = mean(Permutation_importance),
              sd = sd(Permutation_importance)) %>% 
    mutate(sd_hi = mean + sd,
           sd_low = mean - sd) %>% 
    mutate(sd_low = case_when(sd_low < 0 ~ 0, 
                              T ~ sd_low)) %>% 
    arrange(desc(mean))


varimp.plot <- ggplot(varimp_processed.df, aes(fct_reorder(Variable, mean, .desc = F), 
                                               mean)) +
    geom_col(fill = "#6C91C2") +
    # geom_errorbar(aes(ymin = sd_low, 
    #                   ymax = sd_hi),
    #               width = .2) +
    labs(x = "", y = "") +
    # labs(x = "", y = "Relative Variable Importance") +
    scale_y_continuous(position = "right") +
    
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 15, face = "bold"),
          # axis.title.x = element_text(size = 20, face = "bold")
          axis.text.x = element_text(size = 18, face = "bold")
          )

ggsave("img/predictor.png", varimp.plot, 
       height = 10, width = 7.5, dpi = 300)


