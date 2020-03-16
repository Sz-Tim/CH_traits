
library(readxl); library(googlesheets)
library(viridis); library(sf); library(tidyverse); theme_set(theme_bw())
walk(paste0("../1_opfo/code/", c("lc_cols", "00_opfo_fn"), ".R"), source)

dem <- raster::raster("../2_gis/data/VD_21781/dem_VD_21781.tif")
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp") %>%
  filter(!grepl("Lac ", NAME)) 
VD <- st_union(VD_raw)
dem_VD <- raster::crop(dem, VD_raw) %>% raster::mask(., st_zm(VD_raw))

## site and plot
sites.sf <- st_read("../2_gis/data/opfo/MilieuxBDMdiss.shp") %>%
  st_transform(3857)
site.box.sf <- st_read("../2_gis/data/VD_21781/BDMplus_VD_21781.shp") %>%
  st_transform(st_crs(sites.sf))

plot_i <- read_csv("../1_opfo/data/opfo_envDataProcessed.csv")
plot.sf <- st_read("../2_gis/data/opfo/opfo_soil_25per.shp") %>% 
  st_transform(st_crs(sites.sf)) %>%
  mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
  select(Plot_id, propArea) %>% left_join(., plot_i, by="Plot_id") %>% 
  group_by(BDM_id) %>% mutate(nPlots=n()) %>% ungroup %>%
  mutate(el_plot=raster::extract(dem, .))

max_per_sp <- 20

ant <- load_ant_data(clean_spp=T)
ant.all <- ant$all %>%
  mutate(mnt25=raster::extract(dem, .),
         elBin=mnt25 %/% 100 * 100)

measure.df <- read_csv("data/tubes_trait_measurements.csv")




ggplot(filter(measure.df, Genus=="Myrm"), aes(x=mnt25, y=SPECIESID)) + 
  geom_line() + geom_point(alpha=0.5)

p <- ant.all %>% filter(TubeNo %in% filter(measure.df, Genus=="Myrm")$TubeNo) %>% 
  ggplot() + geom_sf(data=VD, fill="gray90", colour="gray70") + 
  geom_sf(aes(colour=mnt25)) + geom_sf(shape=1) + facet_wrap(~SPECIESID, ncol=5) +
  scale_colour_viridis("Elevation (m)") + theme(legend.position=c(0.9, 0.15))
ggsave("~/Desktop/Myrmica_map.pdf", p, height=8, width=12)



pin.df <- ant.all %>% 
  filter(TubeNo %in% filter(measure.df, Genus=="Myrm")$TubeNo) %>%
  left_join(., select(measure.df, TubeNo, BAGNUMBER, Genus), by="TubeNo") %>%
  mutate(lon=st_coordinates(.)[,1], lat=st_coordinates(.)[,2])
  





nora.dir <- "~/Documents/unil/teaching/2020_Khelidj/data/"
write_csv(pin.df %>% st_set_geometry(NULL), 
          paste0(nora.dir, "Myrmica_to_pin.csv"))
write_csv(ant$str %>% st_set_geometry(NULL), 
          paste0(nora.dir, "2019_opfo_strSamp_all.csv"))
write_csv(ant$str %>% st_set_geometry(NULL) %>% filter(GENUSID=="Myrmica"), 
          paste0(nora.dir, "2019_opfo_strSamp_Myrmica.csv"))
write_csv(ant$all %>% 
            mutate(lon=st_coordinates(.)[,1], 
                   lat=st_coordinates(.)[,2],
                   Genus=str_split_fixed(SPECIESID, "_", 2)[,1],
                   mnt25=raster::extract(dem, .),
                   elBin=mnt25 %/% 100 * 100) %>% 
            st_set_geometry(NULL) %>% filter(Genus=="Myrm"), 
          paste0(nora.dir, "2019_opfo_strpub_Myrmica.csv"))
meta.df <- rbind(data.frame(df="Myrmica_to_pin", 
                            names=names(st_set_geometry(pin.df, NULL))),
                 data.frame(df="strSamp_all", 
                            names=names(st_set_geometry(ant$str, NULL))),
                 data.frame(df="strSamp_Myrmica", 
                            names=names(st_set_geometry(ant$str, NULL))),
                 data.frame(df="strpub_Myrmica", 
                            names=names(st_set_geometry(ant$all %>% 
                                                          mutate(lon=st_coordinates(.)[,1], 
                                                                 lat=st_coordinates(.)[,2],
                                                                 Genus=str_split_fixed(SPECIESID, "_", 2)[,1],
                                                                 mnt25=raster::extract(dem, .),
                                                                 elBin=mnt25 %/% 100 * 100), NULL))))

write_csv(meta.df, paste0(nora.dir, "2019_metadata.csv"))
