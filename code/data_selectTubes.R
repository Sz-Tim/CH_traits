# Script for selecting tube numbers

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

# ant.all.sum <- ant.all %>% group_by(SPECIESID) %>%
#   summarise(nTubes=n()) %>%
#   mutate(nToPin=nTubes)
# ant.all.sum$nToPin[ant.all.sum$nToPin > max_per_sp] <- max_per_sp
# sum(ant.all.sum$nToPin)
# 
# 
# sp_counts <- ant.all %>% filter(!is.na(SPECIESID)) %>%
#   group_by(SPECIESID) %>% summarise(nTubes=n()) 
# sp_leMax <- sp_counts$SPECIESID[sp_counts$nTubes <= max_per_sp]
# sp_gMax <- sp_counts$SPECIESID[sp_counts$nTubes > max_per_sp]
# 
# measure.df <- ant.all %>% st_set_geometry(NULL) %>% filter(SPECIESID %in% sp_leMax)
# 
# for(i in seq_along(sp_gMax)) {
#   sp_df <- filter(ant.all, SPECIESID==sp_gMax[i]) %>% st_set_geometry(NULL)
#   el_bins <- table(sp_df$elBin)
#   max_tubes_per_bin <- trunc(max_per_sp/length(el_bins))
#   underrep_bins <- which(el_bins <= max_tubes_per_bin)
#   
#   tubes_per_bin <- el_bins; tubes_per_bin[] <- 0
#   tubes_per_bin[underrep_bins] <- el_bins[underrep_bins]
#   tubes_per_bin[-underrep_bins] <- max_tubes_per_bin
#   tubes_remaining <- el_bins - tubes_per_bin
#   while(sum(tubes_per_bin) < max_per_sp) {
#     rand_bin <- sample(names(tubes_remaining)[tubes_remaining>0], 1)
#     tubes_per_bin[rand_bin] <- tubes_per_bin[rand_bin] + 1
#     tubes_remaining[rand_bin] <- tubes_remaining[rand_bin] - 1
#   }
#   
#   measure.df <- rbind(measure.df, 
#                       map2_df(names(tubes_per_bin), tubes_per_bin, 
#                       ~sp_df %>% filter(elBin==.x) %>% 
#                         sample_n(.y, replace=F,
#                                  weight=(.$source=="s")+.000000001)))
# }
# 
# measure.df <- measure.df %>% 
#   left_join(., ant$str %>% st_set_geometry(NULL) %>% 
#               select(TubeNo, BAGNUMBER), by="TubeNo") %>%
#   mutate(Genus=str_split_fixed(SPECIESID, "_", 2)[,1])
# write_csv(measure.df, "data/pinned_tubes_20_FINAL.csv")



ggplot(measure.df, aes(mnt25)) + geom_histogram(binwidth=100)
ggplot(ant.all, aes(mnt25)) + geom_histogram(binwidth=100)


ggplot() + 
  geom_density(data=measure.df, aes(mnt25), colour=NA, fill="red", alpha=0.5) + 
  geom_density(data=ant.all %>% filter(!is.na(SPECIESID)), aes(mnt25)) + 
  facet_wrap(~SPECIESID, scales="free_y")








ggplot(measure.df, aes(mnt25, SPECIESID)) + 
  geom_point(aes(colour=source), alpha=0.5) + geom_line()




ant.all %>% group_by(SPECIESID) %>%
  summarise(nTubes=n()) %>%
  ggplot(aes(SPECIESID, nTubes)) + geom_point() +
  coord_flip()



p <- ggplot() + geom_sf(data=VD_raw, fill="white", colour="gray") + 
  geom_sf(data=ant.all, alpha=0.3) +  facet_wrap(~SPECIESID)
ggsave("~/Desktop/sp_VD.pdf", p, width=10, height=10)



ant.rng.src <- ant.all %>% group_by(SPECIESID, source) %>%
  summarise(minEl=min(elBin, na.rm=T), 
            medEl=median(elBin, na.rm=T),
            maxEl=max(elBin, na.rm=T),
            minEl.raw=min(mnt25, na.rm=T), 
            medEl.raw=median(mnt25, na.rm=T),
            maxEl.raw=max(mnt25, na.rm=T)) 
ant.rng.W <- filter(ant.rng.src, source=="public")
ant.rng.Y <- filter(ant.rng.src, source=="structured")
ant.rng <- ant.all %>% group_by(SPECIESID) %>%
  summarise(minEl=min(elBin, na.rm=T), 
            medEl=median(elBin, na.rm=T),
            maxEl=max(elBin, na.rm=T)) %>%
  arrange(minEl, maxEl) %>%
  mutate(elOrder=factor(SPECIESID, levels=SPECIESID[row_number()])) %>%
  arrange(medEl) %>%
  mutate(elOrderMed=factor(SPECIESID, levels=SPECIESID[row_number()]))

ant.S <- data.frame(elBin=sort(unique(ant.all$elBin)),
                    nPlots=NA,
                    S_u=NA,
                    S_i=NA,
                    S_u.W=NA,
                    S_i.W=NA,
                    S_u.Y=NA,
                    S_i.Y=NA,
                    area=NA)
for(i in seq_along(ant.S$elBin)) {
  ant.all.i <- filter(ant.all, elBin==ant.S$elBin[i])
  ant.W.i <- filter(ant.all, elBin==ant.S$elBin[i] & source=="public")
  ant.Y.i <- filter(ant.all, elBin==ant.S$elBin[i] & source=="structured")
  
  ant.S$nPlots[i] <- sum(plot.sf$el_plot < (ant.S$elBin[i]+100) &
                           plot.sf$el_plot >= ant.S$elBin[i])
  ant.S$S_u[i] <- n_distinct(ant.all.i$SPECIESID)
  ant.S$S_i[i] <- sum(ant.rng$minEl < (ant.S$elBin[i]+100) &
                             ant.rng$maxEl >= ant.S$elBin[i])
  ant.S$S_u.W[i] <- n_distinct(ant.W.i$SPECIESID)
  ant.S$S_i.W[i] <- sum(ant.rng.W$minEl < (ant.S$elBin[i]+100) &
                          ant.rng.W$maxEl >= ant.S$elBin[i])
  ant.S$S_u.Y[i] <- n_distinct(ant.Y.i$SPECIESID)
  ant.S$S_i.Y[i] <- sum(ant.rng.Y$minEl < (ant.S$elBin[i]+100) &
                          ant.rng.Y$maxEl >= ant.S$elBin[i])
  ant.S$area[i] <- sum(dem_VD@data@values %/% 100 * 100 == ant.S$elBin[i], na.rm=T)
}

ggplot(ant.rng.src, aes(x=SPECIESID, y=medEl, ymin=minEl, 
                        ymax=maxEl, colour=source)) +
  geom_pointrange(fatten=0.5, position=position_dodge(width=0.4)) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
ggplot(full_join(ant.rng.W, ant.rng.Y, by="SPECIESID", suffix=c(".W", ".Y")),
       aes(minEl.raw.W, minEl.raw.Y)) + 
  geom_abline() + geom_point(alpha=0.5) +
  labs(x="Minimum elevation (public)", y="Minimum elevation (structured)")
ggplot(full_join(ant.rng.W, ant.rng.Y, by="SPECIESID", suffix=c(".W", ".Y")),
       aes(medEl.raw.W, medEl.raw.Y)) + 
  geom_abline() + geom_point(alpha=0.5) +
  labs(x="Median elevation (public)", y="Median elevation (structured)")
ggplot(full_join(ant.rng.W, ant.rng.Y, by="SPECIESID", suffix=c(".W", ".Y")),
       aes(maxEl.raw.W, maxEl.raw.Y)) + 
  geom_abline() + geom_point(alpha=0.5) +
  labs(x="Maximum elevation (public)", y="Maximum elevation (structured)")

ggplot(ant.rng.src, aes(x=maxEl-minEl, colour=source)) + geom_density()

gather(ant.S, source, S, 3:8) %>%
  ggplot(aes(x=elBin, y=S)) + stat_smooth(method="loess", se=F, span=1.5) + 
  geom_point() + facet_wrap(~source)
gather(ant.S, source, S, 3:8) %>%
  ggplot(aes(x=log(area), y=S)) + stat_smooth(method="lm", se=F) + 
  geom_point() + facet_wrap(~source)
gather(ant.S, source, S, 3:8) %>%
  ggplot(aes(x=nPlots, y=S)) + stat_smooth(method="loess", se=F, span=1.5) + 
  geom_point() + facet_wrap(~source)
ggplot(ant.S, aes(area, nPlots)) + geom_point()

ggplot(ant.S, aes(elBin)) + 
  geom_point(aes(y=S_i.W-S_i.Y), colour="black", size=3, shape=1) +
  geom_point(aes(y=S_u.W-S_u.Y), colour="red", size=3, shape=1) 


ggplot(ant.rng, aes(x=elOrder, y=medEl, ymin=minEl, ymax=maxEl)) + 
  geom_pointrange(fatten=0.5) + coord_flip()
ggplot(ant.rng, aes(x=elOrderMed, y=medEl, ymin=minEl, ymax=maxEl)) + 
  geom_pointrange(fatten=0.5) + coord_flip()

plot(ant.S$elBin, ant.S$S_i.W, ylim=c(0,55), xlab="Elevation (m)", 
     ylab="Richness", type="l", lwd=2)
lines(ant.S$elBin, ant.S$S_u.W, type="l", lwd=2, col="red")

ggplot(ant.S, aes(log(area))) + 
  stat_smooth(aes(y=S_u.W), method="lm", se=F, colour="red") +
  geom_point(aes(y=S_u.W), colour="red") +
  stat_smooth(aes(y=S_i.W), method="lm", se=F, colour="black") +
  geom_point(aes(y=S_i.W), colour="black")
summary(lm(S_i.W ~ log(area), data=ant.S))
summary(lm(S_u.W ~ log(area), data=ant.S))


# Fully identified genera (no NA if genus is defined):
# - Aphaenogaster
# - Colobopsis
# - Formicoxenus
# - Lasius
# - Leptothorax
# - Manica
# - Myrmecina
# - Solenopsis






sp.prop <- ant.all %>% #st_set_geometry(NULL) %>% 
  filter(!is.na(SPECIESID)) %>%
  group_by(SPECIESID, source) %>% summarise(nTubes=n()) %>% 
  spread(source, nTubes)
sp.prop[is.na(sp.prop)] <- 0
sp.prop <- sp.prop %>% ungroup %>%
  mutate(public_prop=public/sum(public),
         structured_prop=structured/sum(structured)) %>%
  arrange(public_prop-structured_prop) %>%
  mutate(diffOrder=factor(SPECIESID, levels=SPECIESID[row_number()]),
         genus=str_sub(SPECIESID, 1, 5))  %>%
  arrange(public_prop) %>% 
  mutate(publicOrder=factor(SPECIESID, levels=SPECIESID[row_number()])) %>%
  arrange(structured_prop) %>% 
  mutate(structuredOrder=factor(SPECIESID, levels=SPECIESID[row_number()])) 

ggplot(sp.prop, aes(x=public_prop, y=structured_prop)) + geom_point(alpha=0.5)
ggplot(sp.prop, aes(x=public_prop/structured_prop)) + geom_density()
ggplot(sp.prop, aes(x=public_prop/(public_prop+structured_prop))) + geom_density()
ggplot(sp.prop, aes(x=diffOrder, y=public_prop-structured_prop)) + 
  geom_hline(yintercept=0, linetype=2, colour="gray") + 
  geom_text(aes(label=diffOrder), size=3) + coord_flip() + 
  theme(axis.text.y=element_blank(), panel.grid=element_blank()) +
  facet_wrap(~genus, scales="free_y")

ggplot(sp.prop, aes(y=public_prop, x=publicOrder)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=(public_prop-structured_prop)>0), stat="identity") + 
  scale_fill_brewer("Greater\nprevalence", type="qual", palette=2, direction=-1,
                    labels=c("Structured", "Public"), 
                    guide=guide_legend(reverse=T)) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position=c(0.1, 0.87)) +
  # coord_flip() + facet_wrap(~genus, scales="free_y") +
  labs(x="Species", y="Proportion of tubes (public)")
ggsave("../../opfo_str_sampling/eda/sp_relAbund_cs.pdf", width=9, height=6, units="in")
ggplot(sp.prop, aes(y=structured_prop, x=structuredOrder)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=(public_prop-structured_prop)>0), stat="identity") + 
  scale_fill_brewer("Greater\nprevalence", type="qual", palette=2, direction=-1,
                    labels=c("Structured", "Public"), 
                    guide=guide_legend(reverse=T)) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position=c(0.1, 0.87)) +
  # coord_flip() + facet_wrap(~genus, scales="free_y") +
  labs(x="Species", y="Proportion of tubes (structured)")
ggsave("../../opfo_str_sampling/eda/sp_relAbund_ss.pdf", width=9, height=6, units="in")
ggplot(sp.prop, aes(y=public_prop-structured_prop, x=diffOrder)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=(public_prop-structured_prop)>0), stat="identity") + 
  scale_fill_brewer("Greater\nprevalence", type="qual", palette=2, direction=-1,
                    labels=c("Structured", "Public"), 
                    guide=guide_legend(reverse=T)) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position=c(0.1, 0.87)) +
  # coord_flip() + facet_wrap(~genus, scales="free_y") +
  labs(x="Species", y="Difference in relative abundance (public - structured)")
ggsave("../../opfo_str_sampling/eda/sp_relAbund.pdf", width=9, height=6, units="in")
sp.prop %>% arrange(public_prop) %>% 
  select(SPECIESID, public_prop, structured_prop) %>%
  mutate(pubOrder=factor(SPECIESID, levels=SPECIESID[row_number()])) %>%
  gather(source, proportion, 2:3) %>%
  ggplot(aes(y=proportion, x=pubOrder)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=source), stat="identity", position="dodge") + 
  scale_fill_manual("Survey", values=c("red", "gray10"),
                    labels=c("Public", "Structured")) +
  theme(panel.grid=element_blank()) + coord_flip() + 
  labs(x="Species", y="Relative abundance")
sp.prop %>% arrange(structured_prop) %>% 
  select(SPECIESID, public_prop, structured_prop) %>%
  mutate(strOrder=factor(SPECIESID, levels=SPECIESID[row_number()])) %>%
  gather(source, proportion, 2:3) %>%
  ggplot(aes(y=proportion, x=strOrder)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=source), stat="identity", position="dodge") + 
  scale_fill_manual("Survey", values=c("red", "gray10"),
                    labels=c("Public", "Structured")) +
  theme(panel.grid=element_blank()) + coord_flip() + 
  labs(x="Species", y="Relative abundance")



sp.prop.bin <- ant.all %>% #st_set_geometry(NULL) %>% 
  filter(!is.na(SPECIESID)) %>%
  group_by(SPECIESID, elBin, source) %>% summarise(nTubes=n()) %>% 
  spread(source, nTubes)
sp.prop.bin[is.na(sp.prop.bin)] <- 0
sp.prop.bin <- sp.prop.bin %>% group_by(elBin) %>%
  mutate(public_prop=public/sum(public),
         structured_prop=structured/sum(structured)) %>% ungroup %>%
  filter(public > 0 | structured > 0) %>%
  arrange(elBin, public_prop-structured_prop) %>%
  mutate(diffOrder=factor(paste(elBin,SPECIESID), 
                          levels=paste(elBin, SPECIESID)[row_number()]))
sp.prop.bin$structured_prop[is.nan(sp.prop.bin$structured_prop)] <- 0

ggplot(sp.prop.bin, aes(x=diffOrder, y=public_prop-structured_prop)) + 
  geom_hline(yintercept=0, linetype=2, colour="gray") + 
  geom_text(aes(label=SPECIESID), size=2) + coord_flip() + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid=element_blank()) + 
  facet_wrap(~elBin, scales="free_y")
ggplot(sp.prop.bin, aes(x=elBin, y=public_prop-structured_prop)) + 
  geom_hline(yintercept=0, linetype=2, colour="gray") + 
  geom_boxplot(aes(group=elBin))# + geom_point(alpha=0.5)
ggplot(sp.prop.bin, aes(x=elBin, y=public_prop-structured_prop)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=(public_prop-structured_prop)>0), stat="identity") + 
  scale_fill_manual("Greater\nPrevalence", values=c("gray10", "red"), 
                    labels=c("Structured", "Public")) +
  guides(fill=guide_legend(reverse=TRUE)) + 
  theme(panel.grid=element_blank()) + facet_wrap(~SPECIESID) +
  labs(x="Elevation (m)", 
       y="Difference in relative abundance (public - structured)")
ggplot(sp.prop.bin, aes(x=diffOrder, y=public_prop-structured_prop)) + 
  geom_hline(yintercept=0, colour="gray90", size=0.25) +
  geom_bar(aes(fill=(public_prop-structured_prop)>0), stat="identity") + 
  scale_fill_manual("Greater\nPrevalence", values=c("gray10", "red"), 
                    labels=c("Structured", "Public")) +
  theme(panel.grid=element_blank()) + facet_wrap(~elBin, scales="free_y") + 
  coord_flip() +
  labs(x="Species", y="Difference in relative abundance (public - structured)")


