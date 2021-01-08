
library(tidyverse); library(readxl); library(ape); library(nlme); library(lmerTest)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")
walk(paste0("../1_opfo/code/", c("lc_cols", "00_fn"), ".R"), source)
lc_i <- readxl::read_xlsx("../1_opfo/data/landcover_id.xlsx", 1) %>%
  mutate(lcNum=as.numeric(LC_ID))

gis_dir <- "../2_gis/data/VD_21781/"
rast_end <- "_VD_21781.tif"
msr_dir <- "data/img/"
col_dir <- "data/img/"

talk_fonts <- theme(panel.grid=element_blank(),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=10),
                    strip.text=element_text(size=16),
                    plot.title=element_text(size=16))

trait_names <- list(lat=c("WebersLength", "HindTibia", "MidTibia"),
                    fro=c("HeadLength", "HeadWidth", 
                          "InterocularDistance", "ScapeLength"),
                    dor=c("MesosomaWidth", "MesosomaLength", 
                          "HindFemur", "MidFemur"))


ant.ls <- load_ant_data(clean_spp=T)
ant.ls$all$TubeNo <- str_remove(ant.ls$all$TubeNo, "\\.")

ant.ls$all <- ant.ls$all %>% 
  mutate(SampleDate=lubridate::yday(SampleDate),
         GENUSID=str_split_fixed(SPECIESID, "_", 2)[,1]) %>%
  add_covariates(list(mnt25=paste0(gis_dir, "dem", rast_end),
                      GDD0=paste0(gis_dir, "growingDegDays0_envirem", rast_end),
                      AP=paste0(gis_dir, "AP_chelsa", rast_end),
                      npp=paste0(gis_dir, "MODIS_2010-2019", rast_end),
                      TAR=paste0(gis_dir, "TAR_chelsa", rast_end),
                      aspect=paste0(gis_dir, "aspect", rast_end),
                      lc=paste0(gis_dir, "lc_21781.tif"))) %>%
  mutate(aspectN=cos(aspect*pi/180),
         CnpyClosed=(lc_i$Canopy[match(lc, lc_i$lcNum)]=="Closed")*1,
         CnpyMixed=(lc_i$Canopy[match(lc, lc_i$lcNum)]=="Mixed")*1,
         CnpyOpen=(lc_i$Canopy[match(lc, lc_i$lcNum)]=="Open")*1)
trts <- load_traits(ant_i=ant.ls$all, msr_dir=msr_dir, col_dir=col_dir, 
                    na.thresh=0.05, lat_names=trait_names$lat,
                    fro_names=trait_names$fro, dor_names=trait_names$dor)


trait_match <- data.frame(orig=c("v", "HeadLength", "HeadWidth", 
                                 "InterocularDistance", "DVE",
                                 "ScapeLength", "WebersLength",
                                 "HindTibia", "MidTibia", "MesosomaWidth",
                                 "MesosomaLength", "HindFemur", "MidFemur",
                                 "MidLen", "HindLen", "MesoSA", 
                                 "HeadShape", "PronotExp", "ScapeProp"),
                          full=c("Color (lightness)", "Head Length", "Head Width",
                                 "Interocular Dist. (rel)", "Dorsoventral Eye Pos.",
                                 "Scape Length", "Webers Length",
                                 "Hind Tibia", "Mid Tibia", "Pronotum Width",
                                 "Mesosoma Length", "Hind Femur", "Mid Femur",
                                 "Rel. Leg Length (mid)", "Rel. Leg Length (hind)", 
                                 "Mesosoma SA", "Head Shape",
                                 "Pronotal Expansion", "Scape Proportion"))


focus_traits <- trait_match$full
focus_traits <- c("Color (lightness)", "Webers Length", 
                  #"Pronotum Width", 
                  "Dorsoventral Eye Pos.", "Scape Length",
                  "Head Width", #"Head Length", #"Interocular Dist. (rel)",
                  "Rel. Leg Length (hind)", "Head Shape",
                  "Pronotal Expansion", "Scape Proportion")


wkr.df_3 <- filter(trts$wkr.wide, n_clny > 3) %>% select(-`NA`) %>%
  # filter(GENUSID %in% c("Myrm")) %>%
  mutate(MidLen=MidLen/WebersLength, HindLen=HindLen/WebersLength,
         DVE=InterocularDistance/HeadWidth,
         InterocularDistance=(HeadWidth-InterocularDistance)/HeadLength,
         HeadShape=HeadWidth/HeadLength,
         PronotExp=MesosomaWidth/WebersLength,
         ScapeProp=ScapeLength/HeadLength) %>%
  select(one_of("Worker", "TubeNo", "SPECIESID", 
                trait_match$orig[trait_match$full %in% focus_traits])) %>%
  pivot_longer(cols=3+(seq_along(focus_traits)), names_to="Trait", values_to="Value") %>%
  group_by(Trait) %>% 
  mutate(Value_std=scale(Value)[,1],
         GENUSID=str_sub(SPECIESID, 1, 4),
         Trait_orig=Trait,
         Trait=trait_match$full[match(Trait, trait_match$orig)]) %>%
  filter(!is.na(Value))
clny.df_3 <- wkr.df_3 %>% group_by(TubeNo, SPECIESID, Trait) %>%
  summarise(mnValue=mean(Value, na.rm=T),
            sdValue=sd(Value, na.rm=T),
            mnValue_std=mean(Value_std, na.rm=T),
            sdValue_std=sd(Value_std, na.rm=T))

varPart_cols <- c("Within colonies"="#a6611a", 
                  "Among colonies"="#dfc27d",
                  "Among species"="#80cdc1",
                  "Among genera"="#018571")
Trait_varcomp <- map_dfr(unique(wkr.df_3$Trait), 
    ~varcomp(lme(Value ~ 1, random=~1|GENUSID/SPECIESID/TubeNo, 
                 data=filter(wkr.df_3, Trait==.x & !is.na(Value))), TRUE)) %>%
  mutate(Trait=unique(wkr.df_3$Trait)) %>%
  rename(`Within colonies`=Within, 
         `Among colonies`=TubeNo, 
         `Among species`=SPECIESID,
         `Among genera`=GENUSID) %>%
  pivot_longer(-5, names_to="Effect", values_to="Variance") %>%
  mutate(Effect=factor(Effect, levels=rev(unique(Effect))))
Trait_varcomp %>%
  filter(Trait %in% focus_traits) %>% 
  # filter(grepl("Colonies", Effect)) %>%
  mutate(Trait=factor(Trait, levels=focus_traits)) %>%
  ggplot(aes(Trait, Variance, fill=Effect)) + 
  geom_bar(position="fill", stat="identity", colour="gray30") + coord_flip() +
  scale_fill_manual("", values=varPart_cols) + 
  labs(x="", y="Proportion of variance") + talk_fonts
ggsave("eda/varpart_all.jpg", width=10, height=6)
Trait_varcomp %>%
  filter(Trait %in% focus_traits) %>% 
  filter(!grepl("genera", Effect)) %>%
  mutate(Trait=factor(Trait, levels=focus_traits)) %>%
  ggplot(aes(Trait, Variance, fill=Effect)) + 
  geom_bar(position="fill", stat="identity", colour="gray30") + coord_flip() +
  scale_fill_manual("", values=varPart_cols) + 
  labs(x="", y="Proportion of variance (within genera)") + talk_fonts
ggsave("eda/varpart_WithinGen.jpg", width=10, height=6)
Trait_varcomp %>%
  filter(Trait %in% focus_traits) %>% 
  filter(grepl("colonies", Effect)) %>%
  mutate(Trait=factor(Trait, levels=focus_traits)) %>%
  ggplot(aes(Trait, Variance, fill=Effect)) + 
  geom_bar(position="fill", stat="identity", colour="gray30") + coord_flip() +
  scale_fill_manual("", values=varPart_cols) + 
  labs(x="", y="Proportion of variance (intraspecific)") + talk_fonts
ggsave("eda/varpart_WithinSpp.jpg", width=10, height=6)

  
# var_part.df <- wkr.df_3 %>% group_by(SPECIESID, Trait) %>%
#   summarise(#Fval=rand(lmer(Value ~ (1|TubeNo)))$`F value`[1],
#             # MS_a=rand(lmer(Value ~ (1|TubeNo)))$`Mean Sq`[1],
#             # MS_w=rand(lmer(Value ~ (1|TubeNo)))$`Mean Sq`[2],
#             Pval=rand(lmer(Value ~ (1|TubeNo)))$`Pr(>Chisq)`[2],
#             #Fval_std=rand(lmer(Value_std ~ (1|TubeNo)))$`F value`[1],
#             #MS_a_std=rand(lmer(Value_std ~ (1|TubeNo)))$`Mean Sq`[1],
#             #MS_w_std=rand(lmer(Value_std ~ (1|TubeNo)))$`Mean Sq`[2],
#             Pval_std=rand(lmer(Value_std ~ (1|TubeNo)))$`Pr(>Chisq)`[2]) %>%
#   mutate(GENUSID=str_sub(SPECIESID, 1, 4))
# 
# ggplot(var_part.df, aes(Trait, Fval)) + geom_boxplot() + ylim(0,20) + 
#   labs(y="MS[among colonies] / MS[within colonies]") + coord_flip() + 
#   facet_wrap(~GENUSID)
# 
# ggplot(var_part.df, aes(SPECIESID, Fval)) + geom_boxplot() + ylim(0,20) + 
#   labs(y="MS[among colonies] / MS[within colonies]") + coord_flip()
# 
# ggplot(var_part.df, aes(GENUSID, Fval)) + geom_boxplot() + ylim(0,20) + 
#   labs(y="MS[among colonies] / MS[within colonies]") + coord_flip() + 
#   facet_wrap(~Trait)
# 
# var_part.df %>%
#   filter(GENUSID != "Form") %>%
#   filter(Trait %in% c("v", "WebersLength", 
#                       "HeadWidth",# "HeadLength", 
#                       "InterocularDistance", #"ScapeLength", 
#                       "MidLen", "HindLen")) %>%
# ggplot(aes(Trait, Fval, colour=GENUSID)) + geom_boxplot() + 
#   facet_grid(Trait~., scales="free_y") +
#   labs(y="MS[among colonies] / MS[within colonies]") + coord_flip() 
# 
# var_part.df %>%
#   filter(GENUSID != "Form") %>%
#   filter(Trait %in% c("v", "WebersLength", 
#                       "HeadWidth",# "HeadLength", 
#                       "InterocularDistance", #"ScapeLength", 
#                       "MidLen", "HindLen")) %>%
#   ggplot(aes(Trait, MS_a_std, colour=GENUSID)) + geom_boxplot() + 
#   facet_grid(Trait~., scales="free_y") +
#   labs(y="MS[among colonies]") + coord_flip() 
# 
# var_part.df %>%
#   filter(GENUSID != "Form") %>%
#   filter(Trait %in% c("v", "WebersLength", 
#                       "HeadWidth",# "HeadLength", 
#                       "InterocularDistance", #"ScapeLength", 
#                       "MidLen", "HindLen")) %>%
#   ggplot(aes(Trait, MS_w_std, colour=GENUSID)) + geom_boxplot() + 
#   facet_grid(Trait~., scales="free_y") +
#   labs(y="MS[within colonies]") + coord_flip() 
# 
# var_part.df %>%
#   filter(GENUSID != "Form") %>%
#   filter(Trait %in% c("v", "WebersLength", 
#                       "HeadWidth",# "HeadLength", 
#                       "InterocularDistance", #"ScapeLength", 
#                       "MidLen", "HindLen")) %>%
#   ggplot(aes(MS_w_std, MS_a_std, colour=GENUSID)) + geom_point() + 
#   geom_abline() + facet_wrap(~Trait) + xlim(0,2) + ylim(0,2)





wkr.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  ggplot(aes(Value, SPECIESID, fill=GENUSID)) + 
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001) +
  scale_fill_brewer("", type="qual", palette=2) + 
  facet_wrap(~Trait, scales="free_x") +
  labs(x="Value", y="", title="Worker traits") + 
  theme(legend.position="none") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))
ggsave("eda/var_violin_wkrs.jpg", width=12, height=12)

clny.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  ggplot(aes(mnValue, SPECIESID, fill=GENUSID)) + 
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001) +
  scale_fill_brewer("", type="qual", palette=2) + 
  facet_wrap(~Trait, scales="free_x") +
  labs(x="Within-colony mean", y="", title="Means within colonies") +
  theme(legend.position="none") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))
ggsave("eda/var_violin_clny.jpg", width=12, height=12)

wkr.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(mn=mean(Value, na.rm=T), sd=sd(Value, na.rm=T), CV=sd/mn) %>%
  ggplot(aes(CV, SPECIESID, colour=GENUSID)) + 
  geom_vline(xintercept=0) + geom_point() + 
  scale_colour_brewer("", type="qual", palette=2) + 
  facet_wrap(~Trait) + 
  scale_x_continuous(labels=scales::percent) +
  labs(x="CV", y="", title="Within-species CV") + 
  theme(legend.position="none") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))

wkr.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(mn=mean(Value, na.rm=T), sd=sd(Value, na.rm=T), CV=sd/mn) %>%
  ggplot(aes(CV, Trait, fill=GENUSID)) + 
  geom_vline(xintercept=0) + 
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001, alpha=0.5) +
  scale_fill_brewer("", type="qual", palette=2) + 
  scale_x_continuous(labels=scales::percent) +
  labs(x="CV", y="", title="Within-species CV") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))

clny.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(SD=sd(mnValue_std, na.rm=T)) %>%
  ggplot(aes(SPECIESID, SD, group=SPECIESID, colour=GENUSID)) + 
  geom_point(size=4) + ylim(0, NA) + 
  scale_colour_brewer(type="qual", palette=2) + facet_wrap(~Trait, scales="free_x") +
  theme(legend.position="none") + #talk_fonts + 
  labs(x="", y="Among-colony sd", title="Variation among colonies") +
  coord_flip() + theme(panel.grid=element_blank())

clny.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(SD=sd(mnValue_std, na.rm=T)) %>%
  ggplot(aes(GENUSID, SD, group=GENUSID, colour=GENUSID)) + 
  geom_point(size=2) + ylim(0, NA) + 
  theme(legend.position="none") + #talk_fonts + 
  scale_colour_brewer(type="qual", palette=2) + facet_wrap(~Trait) +
  labs(x="", y="Among-colony sd (z-transformed)", title="Variation among colonies") +
  coord_flip() + theme(panel.grid=element_blank())

clny.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(SD=sd(mnValue_std, na.rm=T)) %>%
  ggplot(aes(SD, Trait)) + ggridges::geom_density_ridges() +
  labs(x="Among-colony sd (z-transformed)", y="", title="Variation among colonies") +
  theme(panel.grid=element_blank())

clny.df_3 %>%# filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(CV=sd(mnValue, na.rm=T)/mean(mnValue, na.rm=T)) %>%
  ggplot(aes(SPECIESID, CV, group=SPECIESID, colour=GENUSID)) + 
  # geom_boxplot(outlier.size=0.25) + 
  geom_point(size=4) + ylim(0, NA) + 
  scale_colour_brewer(type="qual", palette=2) + facet_wrap(~Trait) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="", y="Among-colony sd / among-colony mean", title="CV among colonies") +
  coord_flip() + theme(panel.grid=element_blank())

clny.df_3 %>%# filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(CV=sd(mnValue, na.rm=T)/mean(mnValue, na.rm=T)) %>%
  ggplot(aes(CV, Trait)) + ggridges::geom_density_ridges() +
  labs(x="Among-colony sd / among-colony mean", y="", title="CV among colonies") +
  theme(panel.grid=element_blank())

clny.df_3 %>%# filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  group_by(GENUSID, SPECIESID, Trait) %>%
  summarise(CV=sd(mnValue, na.rm=T)/mean(mnValue, na.rm=T)) %>%
  ggplot(aes(CV, Trait, fill=GENUSID)) + 
  geom_vline(xintercept=0) + 
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001, alpha=0.5) +
  scale_fill_brewer("", type="qual", palette=2) + 
  scale_x_continuous(labels=scales::percent) +
  labs(y="", x="Among-colony sd / among-colony mean", title="CV among colonies") +
  theme(panel.grid=element_blank())
ggsave("eda/var_violin_clny_CV_bySpp.jpg", width=12, height=12)

clny.df_3 %>% #filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  # filter((sdValue/mnValue) < 0.5) %>%
  ggplot(aes(log(sdValue), SPECIESID, fill=GENUSID)) + 
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001) +
  scale_fill_brewer("", type="qual", palette=2) + facet_wrap(~Trait, scales="free_x") +
  labs(x="Within-colony sd", y="", title="standard deviation within colonies") +
  theme(legend.position="none") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))
ggsave("eda/var_violin_clny_sd.jpg", width=12, height=12)

clny.df_3 %>%# filter(SPECIESID != "Form_cuni") %>%
  mutate(GENUSID=str_sub(SPECIESID, 1, 4)) %>%
  filter(Trait %in% focus_traits) %>%
  filter((sdValue/mnValue) < 0.5) %>%
  ggplot(aes(sdValue/mnValue, SPECIESID, fill=GENUSID)) + 
  geom_vline(xintercept=0, colour="gray80", size=0.2) +
  ggridges::geom_density_ridges(size=0.25, rel_min_height=0.0001) +
  scale_fill_brewer("", type="qual", palette=2) + facet_wrap(~Trait) +
  scale_x_continuous(labels=scales::percent_format(accuracy=1)) + 
  labs(x="Within-colony CV", y="", title="CV within colonies") +
  theme(legend.position="none") + #talk_fonts + 
  theme(panel.grid.major.y=element_line(colour="gray80", size=0.2))
ggsave("eda/var_violin_clny_CV.jpg", width=12, height=12)








###--- Questions: 
# Which traits are more variable within species?

# Do colonies differ in their trait distributions? Or do they represent random draws from the species-level distribution?








### PCA
library(vegan); library(viridis)
spp_10 <- unique(filter(trts$wkr.wide, n_clny>3)$SPECIESID)

# for(i in seq_along(spp_10)) {
  wkr.df_3_wide <- wkr.df_3 %>% 
    # filter(SPECIESID %in% spp_10) %>%
    # filter(GENUSID == "Myrm") %>%
    # filter(SPECIESID==spp_10[i]) %>%
    select(-Value_std, -Trait_orig) %>%
    pivot_wider(names_from="Trait", values_from=c("Value")) %>%
    mutate(el=trts$clny.wide$mnt25[match(TubeNo, trts$clny.wide$TubeNo)])
  complete.df <- wkr.df_3_wide[complete.cases(wkr.df_3_wide),]
  # complete.df <- complete.df %>% mutate(across(focus_traits[-1], ~./`Webers Length`))
  #pairs(complete.df[,5:19], lower.panel=panel.smooth, cex=0.5)
  
  pca <- princomp(select(complete.df, all_of(focus_traits))) 
  pc <- prcomp(select(complete.df, all_of(focus_traits)))
  km_scree <- map_dbl(1:20, ~sum(kmeans(pc$x[,1:2], ., iter.max=1e3, nstart=25)$withinss))
  
  km <- kmeans(pc$x[,1:2], 4, iter.max=1e4, nstart=50)
  complete.df$cluster <- as.character(km$cluster)
  complete.df$PC1 <- pc$x[,1]
  complete.df$PC2 <- pc$x[,2]
  complete.df$PC3 <- pc$x[,3]
  
  p <- ggplot(complete.df, aes(x=PC1, y=PC2, shape=SPECIESID, colour=GENUSID)) + 
    geom_point(size=3) + 
    scale_shape_manual(values=LETTERS[1:n_distinct(complete.df$SPECIESID)]) +
    scale_colour_brewer(type="qual", palette=2) #+ ggtitle(spp_10[i])
  print(p)

# }


pca <- rda(select(complete.df, all_of(focus_traits)))


{biplot(pca, type=c("text", "points"), cex=1.5)
  ordihull(pca, group=complete.df$SPECIESID, lwd=2,
           col=viridis::viridis_pal()(n_distinct(complete.df$SPECIESID)))
  legend("topright", col=viridis::viridis_pal()(n_distinct(complete.df$SPECIESID)),
         lty=1, lwd=2, levels(factor(complete.df$SPECIESID)))
}











ant.ls$str %>% st_set_geometry(NULL) %>% #filter(TypeOfSample=="soil") %>%
  mutate(montane=mnt25>1000) %>% 
  group_by(GENUSID, montane) %>% summarise(N=n()) %>% 
  group_by(montane) %>% mutate(P=N/sum(N)) %>% 
  ggplot(aes(montane, P, fill=GENUSID)) + 
  geom_hline(yintercept=0, colour="gray30", size=0.5) + 
  geom_bar(stat="identity", colour="gray30", size=0.5) + 
  scale_fill_viridis_d("") +
  scale_x_discrete(labels=c("<1000", ">1000")) + 
  labs(x="Elevation (m)", y="Proportion of collections", title="Structured") + 
  talk_fonts
ggsave("eda/prop_tubes_genus.jpg", width=6, height=8)









library(sf)
wkr.cor <- filter(trts$wkr.wide, n_clny > 3) %>% 
  select(-`NA`, -unit, -contains("mn_"), -contains("med_"), -contains("sd_"), 
         -rgbCol, -h, -s, -source, -geometry) %>%
  mutate(MidLen=MidLen/WebersLength, HindLen=HindLen/WebersLength,
         InterocularDistance=(HeadWidth-InterocularDistance)/HeadLength,
         HeadShape=HeadWidth/HeadLength,
         PronotExp=MesosomaWidth/WebersLength,
         ScapeProp=ScapeLength/HeadLength)
round(cor(select(wkr.cor, -Worker, -TubeNo, -SPECIESID, -GENUSID), use="pairwise"), 3)

pairs(select(wkr.cor, -Worker, -TubeNo, -SPECIESID, -GENUSID), 
      lower.panel=panel.smooth, cex=0.5)





rel.traits.df <- wkr.df_3 %>% 
  filter(Trait %in% c("Color (lightness)", 
                      "Dorsoventral Eye Pos.",
                      "Head Shape", 
                      "Webers Length",
                      # "Interocular Dist. (rel)", 
                      "Pronotal Expansion", 
                      "Rel. Leg Length (hind)",
                      # "Rel. Leg Length (mid)", 
                      "Scape Proportion")) %>% 
  select(-Value, -Trait_orig) %>%
  pivot_wider(names_from="Trait", values_from="Value_std") #%>%
  # left_join(., select(trts$clny.df, TubeNo, GDD0, AP, npp, aspectN))
complete.rel <- rel.traits.df[complete.cases(rel.traits.df),]
rel.traits.clny <- rel.traits.df %>% group_by(GENUSID, SPECIESID, TubeNo) %>%
  summarise(across(where(is.numeric), list(mn=mean, sd=sd)))

library(GGally)
rel.traits.df %>% #filter(GENUSID %in% c("Lasi")) %>%
  ggpairs(columns=5:ncol(.), 
        aes(colour=GENUSID, group=SPECIESID, fill=GENUSID),
        lower=list(continuous=wrap('smooth', shape=1, size=0.5, se=F)),
        diag=list(continuous=wrap('densityDiag', alpha=0.1, size=0.25)),
        upper=list(continuous=wrap('cor', size=3))) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) 

rel.traits.clny %>% #filter(GENUSID %in% c("Lasi")) %>%
  select(GENUSID, SPECIESID, TubeNo, contains("_mn")) %>% 
  ggpairs(columns=4:ncol(.), 
          aes(colour=GENUSID, group=SPECIESID, fill=GENUSID),
          lower=list(continuous=wrap('smooth', shape=1, size=0.5, se=F)),
          diag=list(continuous=wrap('densityDiag', alpha=0.1, size=0.25)),
          upper=list(continuous=wrap('cor', size=3))) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) 

rel.traits.clny %>% #filter(GENUSID %in% c("Lasi")) %>%
  select(GENUSID, SPECIESID, TubeNo, contains("_sd")) %>% 
  mutate(across(contains("_sd"), log)) %>%
  ggpairs(columns=4:ncol(.), 
          aes(colour=GENUSID, group=SPECIESID, fill=GENUSID),
          lower=list(continuous=wrap('smooth', shape=1, size=0.5, se=F)),
          diag=list(continuous=wrap('densityDiag', alpha=0.1, size=0.25)),
          upper=list(continuous=wrap('cor', size=3))) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) 

rel.traits.clny %>% #filter(GENUSID %in% c("Lasi")) %>%
  select(TubeNo, contains("_sd")) %>% 
  mutate(across(contains("_sd"), log)) %>%
  right_join(select(trts$clny.wide, TubeNo, GENUSID, SPECIESID, 
                    GDD0, AP, aspectN) %>%
               mutate(Cnpy=case_when(trts$clny.wide$CnpyClosed==1~"Closed",
                                     trts$clny.wide$CnpyMixed==1~"Mixed",
                                     trts$clny.wide$CnpyOpen==1~"Open")), .) %>%
  ggpairs(columns=4:ncol(.), 
          aes(colour=GENUSID, group=SPECIESID, fill=GENUSID),
          lower=list(continuous=wrap('smooth', shape=1, size=0.5, se=F)),
          diag=list(continuous=wrap('densityDiag', alpha=0.1, size=0.25)),
          upper=list(continuous=wrap('cor', size=3))) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) 

rel.traits.clny %>% #filter(GENUSID %in% c("Lasi")) %>%
  select(TubeNo, contains("_mn")) %>% 
  right_join(select(trts$clny.wide, TubeNo, GENUSID, SPECIESID, 
                    GDD0, AP, aspectN), .) %>%
  ggpairs(columns=4:ncol(.), 
          aes(colour=GENUSID, group=SPECIESID, fill=GENUSID),
          lower=list(continuous=wrap('smooth', shape=1, size=0.5, se=F)),
          diag=list(continuous=wrap('densityDiag', alpha=0.1, size=0.25)),
          upper=list(continuous=wrap('cor', size=3))) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) 




library(vegan)
pca.data <- filter(complete.rel, !GENUSID %in% c("Form"))
pca.data <- complete.rel %>% filter(GENUSID=="Myrm")
pca <- rda(pca.data[,-(1:4)])


{biplot(pca, type=c("text", "points"), cex=1.5)
  ordihull(pca, group=pca.data$SPECIESID, lwd=2,
           col=viridis::viridis_pal()(n_distinct(pca.data$SPECIESID)))
  legend("topleft", col=viridis::viridis_pal()(n_distinct(pca.data$SPECIESID)),
         lty=1, lwd=2, levels(factor(pca.data$SPECIESID)), cex=0.5)
}
