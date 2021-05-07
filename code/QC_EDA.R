# CH ant trait biogeography
# Trait measurements and color processing


library(tidyverse); library(readxl); library(googlesheets); library(viridis); 
library(GGally); library(sf)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")
walk(paste0("../1_opfo/code/", c("lc_cols", "00_fn"), ".R"), source)

gis_dir <- "../2_gis/data/VD_21781/"
rast_end <- "_VD_21781.tif"
msr_dir <- "data/img/"
col_dir <- "data/img/"


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
                      annPET=paste0(gis_dir, "annualPET_envirem", rast_end),
                      AP=paste0(gis_dir, "AP_chelsa", rast_end),
                      npp=paste0(gis_dir, "MODIS_2010-2019", rast_end),
                      DTR=paste0(gis_dir, "DTR_chelsa", rast_end),
                      TAR=paste0(gis_dir, "TAR_chelsa", rast_end),
                      aspect=paste0(gis_dir, "aspect", rast_end))) %>%
  mutate(aspectN=cos(aspect*pi/180))
trts <- load_traits(ant_i=ant.ls$all, msr_dir=msr_dir, col_dir=col_dir, 
                    na.thresh=0.05, lat_names=trait_names$lat,
                    fro_names=trait_names$fro, dor_names=trait_names$dor)




# Data QC
qc.done <- read_csv("data/data_QC/msr_QC_checked.csv")
trait_name.df <- enframe(trait_names) %>% unnest(cols=c(value))

map(unique(trts$wkr.df$SPECIESID),
    ~ggplot(filter(trts$wkr.std, !is.na(Value) & SPECIESID==.), 
            aes(x=mnt25, y=Value)) + facet_wrap(~Trait) + 
      ggtitle(.) + geom_point(aes(colour=abs(Value)>qnorm(0.995) & 
                                    !TubeNo %in% qc.done$TubeNo)) + 
      geom_text(data=filter(trts$wkr.std, 
                            !is.na(Value) & 
                              SPECIESID==. &
                              abs(Value)>qnorm(0.995) &
                              !TubeNo %in% qc.done$TubeNo),
                aes(label=paste(TubeNo, Worker, sep="-")), size=2.5) +
      scale_colour_manual(values=c(`FALSE`="gray30", `TRUE`="red")) +
      theme(legend.position="none"))

qc.todo <- trts$wkr.std %>% 
  filter(abs(Value)>qnorm(0.995) & !TubeNo %in% qc.done$TubeNo) %>%
  mutate(img=paste0(trait_name.df$name[match(Trait, trait_name.df$value)],
                    Worker)) %>% 
  select(SPECIESID, TubeNo, img, Trait, Value) %>%
  arrange(SPECIESID, TubeNo, img, Trait)
qc.todo %>% print.AsIs
write_csv(qc.todo, "data/data_QC/msr_QC_todo.csv")





# Quantiles
ggplot(filter(trts$wkr.wide, !is.na(v) & n_clny>2), aes(x=mnt25, y=v)) + 
  geom_point(alpha=0.25) + 
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="gray", linetype=2) +
  stat_smooth(method="lm", se=F, aes(group=SPECIESID)) +
  theme(legend.position="none") + facet_wrap(~SPECIESID)
ggplot(filter(trts$clny.wide, !is.na(v) & n_clny>2), aes(x=mnt25, y=v)) + 
  geom_point(alpha=0.25) + facet_wrap(~GENUSID) +
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="gray", linetype=2) +
  stat_smooth(method="lm", se=F, aes(group=SPECIESID)) +
  theme(legend.position="none") 
ggplot(filter(trts$wkr.df, !is.na(Value) & n_clny>2), aes(x=mnt25, y=Value)) + 
  # geom_point(alpha=0.15, shape=1) + 
  stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x, se=F, size=0.2) +
  stat_smooth(method="lm", formula=y~x, se=F, size=0.5, colour="black") +
  facet_wrap(~Trait)

ggplot(filter(trts$wkr.std, !is.na(Value)), aes(x=mnt25, y=Value)) + 
  geom_point(alpha=0.25, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", se=F, colour="black", linetype=2, size=0.5) +
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="black", linetype=3) +
  facet_wrap(~Trait)
ggplot(filter(trts$wkr.std, !is.na(Value)), aes(x=mnt25, y=Value)) + 
  geom_point(alpha=0.25, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", se=F, colour="black", size=0.5) +
  stat_smooth(method="lm", se=F, colour="black", linetype=2, size=0.3, 
              aes(group=SPECIESID)) +
  facet_grid(Trait~GENUSID, scales="free")
ggplot(filter(trts$wkr.std, !is.na(Value)), aes(x=mnt25, y=Value)) + 
  geom_point(alpha=0.25, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", se=F, colour="black", linetype=2, size=0.5) +
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="black", linetype=3) +
  facet_grid(Trait~SPECIESID, scales="free")
ggplot(filter(trts$wkr.df, !is.na(Value)), aes(x=mnt25, y=Value)) + 
  geom_point(aes(colour=GENUSID), alpha=0.25) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", se=F, colour="black", linetype=2, size=0.5) +
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="black", linetype=3) +
  facet_wrap(~Trait)
ggplot(filter(trts$clny.std, !is.na(mnValue)), 
       aes(x=mnt25, y=mnValue, colour=GENUSID)) + 
  geom_point(alpha=0.25) + 
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="red", linetype=2) +
  facet_wrap(~Trait)
ggplot(filter(trts$clny.df, !is.na(mnValue)), aes(x=mnt25, y=mnValue)) + 
  geom_point(aes(colour=GENUSID), alpha=0.25) + 
  stat_quantile(quantiles=c(.1,.9), formula=y~x, colour="red", linetype=2) +
  facet_wrap(~Trait)

ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=WebersLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=HeadWidth)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=(HindLen+HindFemur)/WebersLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=ScapeLength/HeadLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=InterocularDistance/HeadWidth)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$wkr.wide, n_clny>2), 
       aes(x=mnt25, y=MesosomaWidth/WebersLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")

ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=mnValue_WebersLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=log(sdValue_WebersLength))) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=mnValue_HeadWidth)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=log(sdValue_HeadWidth))) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=mnValue_HeadLength)) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")
ggplot(filter(trts$clny.wide, n_clny>2), 
       aes(x=mnt25, y=log(sdValue_HeadLength))) + 
  geom_point(shape=1, aes(colour=GENUSID)) + 
  scale_colour_brewer(palette=2, type="qual") +
  stat_smooth(method="lm", colour="black", linetype=2, size=0.5) +
  facet_wrap(~SPECIESID, scales="free")






# colors
trts$wkr.df %>% filter(!is.na(med_R)) %>%
  # filter(GENUSID %in% c("Lept", "Temn")) %>%
  plot_colors_by_v("mnt25", "Elevation (m)")


# traits across environmental variables
ggplot(trts$clny.std %>% filter(n_clny > 2), aes(GDD5, mnValue)) + 
  facet_wrap(~Trait, scales="free", drop=T) +
  stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x, se=F, size=0.3) +
  stat_smooth(method="lm", formula=y~x, se=F, size=0.7, colour="red") 
ggplot(trts$clny.std %>% filter(n_clny > 2), aes(AP, log(sdValue))) + 
  facet_wrap(~Trait, scales="free", drop=T) +
  stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x, se=F, size=0.3) +
  stat_smooth(method="lm", formula=y~x, se=F, size=0.7, colour="red") 




# trait variability within colonies
ggplot(trts$clny.df %>% filter(n_clny > 2), aes(mnt25, log(sdValue))) + 
  facet_wrap(~Trait, scales="free") +
  # geom_point(alpha=0.1) +
  stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x, se=F, size=0.3) +
  stat_smooth(method="lm", formula=y~x, se=F, size=0.7, colour="red")
ggplot(trts$clny.df, aes(mnt25, sdValue/mnValue, group=SPECIESID)) + 
  facet_wrap(~Trait, scales="free") +
  geom_point(alpha=0.4, shape=1) +
  stat_smooth(method="lm", se=F, size=0.5, colour="gray30") 
ggplot(trts$clny.wide, aes(mnt25, v_CV)) + geom_point() + 
  stat_smooth(method="lm", se=F, colour="gray30") + #facet_wrap(~SPECIESID, scales="free") + 
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray", size=.5) 
ggplot(trts$clny.std, aes(mnt25, log(v_var))) + geom_point()  + 
  stat_smooth(method="lm", se=F, colour="gray30") +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray", size=.5) 
ggplot(trts$clny.std, aes(v, log(v_var))) + geom_point(alpha=0.8)  + 
  stat_smooth(method="lm", se=F, colour="gray30") +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray", size=.5) 





trts$wkr.std %>% group_by(TubeNo, Worker) %>% 
  summarise(across(everything(), first)) %>%
  ggplot(aes(SampleDate, v)) + geom_point(alpha=0.5, shape=1) +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") +
  facet_wrap(~SPECIESID)
ggplot(trts$clny.std, aes(SampleDate, log(v_var))) + geom_point()  + 
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") 

ggplot(trts$wkr.wide, aes(SampleDate, WebersLength)) + geom_point(alpha=0.5) +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") +
  facet_wrap(~SPECIESID, scales="free_y")

ggplot(trts$wkr.wide, aes(SampleDate, HindLen/WebersLength)) + geom_point(alpha=0.5) +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") +
  facet_wrap(~SPECIESID, scales="free_y")

ggplot(trts$wkr.wide, aes(aspectN, v)) + geom_point(alpha=0.5) +
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") +
  facet_wrap(~SPECIESID, scales="free_y")
ggplot(trts$clny.df, aes(aspectN, log(v_var))) + geom_point()  + 
  stat_smooth(aes(group=SPECIESID), method="lm", se=F, colour="gray") 


ggplot(trts$wkr.wide, aes(WebersLength, v)) + geom_point() + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, colour="gray30") + 
  facet_wrap(~SPECIESID, scales="free")

ggplot(trts$wkr.std %>% filter(Trait=="WebersLength"), aes(Value, v)) +
  geom_point() + 
  stat_smooth(aes(group=SPECIESID, colour=GENUSID), method="lm", 
              formula=y~x, se=F, size=1) +
  stat_smooth(method="lm", formula=y~x, se=F, size=2, colour="red") +
  scale_colour_brewer("", type="qual", palette=2)

ggplot(trts$wkr.df %>% filter(Trait=="WebersLength"), aes(Value, v)) +
  geom_point(shape=1) + 
  stat_smooth(aes(group=SPECIESID), method="lm", 
              formula=y~x, se=F, size=0.5) +
  stat_smooth(method="lm", formula=y~x, se=F, size=1, colour="red") +
  facet_wrap(~GENUSID, scales="free")

ggplot(trts$wkr.std, aes(Value, v)) + geom_point(alpha=0.1) + 
  stat_smooth(aes(group=SPECIESID, colour=GENUSID), method="lm", 
              formula=y~x, se=F, size=0.3) +
  stat_smooth(method="lm", formula=y~x, se=F, size=0.7, colour="red") +
  scale_colour_brewer("", type="qual", palette=2) + 
  facet_wrap(~Trait) + xlim(-4, 4) + ylim(-4, 4)

ggplot(trts$wkr.wide, aes(mnt25, InterocularDistance/HeadWidth)) + geom_point() + 
  stat_smooth(method="lm") + facet_wrap(~SPECIESID, scales="free")




library(quantreg)
rq(v ~ I(mnt25/100), tau=c(0.05, 0.5, 0.95), data=trts$wkr.wide)
library(lme4)
summary(lmer(v ~ I(mnt25/100) + (1|SPECIESID), data=trts$clny.wide))
summary(lmer(mnValue_WebersLength ~ I(mnt25/100) + (1|SPECIESID), 
             data=trts$clny.wide))












# Some PCA
complete.df <- trts$wkr.wide[complete.cases(trts$wkr.wide[,c(16,34:40,42:47)]),]
complete.df[,c(34:40,42:47)] <- complete.df[,c(34:40,42:47)]/complete.df$WebersLength
pairs(complete.df[,c(16,34:40,42:47)], lower.panel=panel.smooth, cex=0.5)

library(vegan)
pca <- rda(complete.df[,c(16,34:40,42:47)])

pdf("eda/00_PCA_Myrmica.pdf", width=8, height=8)
{biplot(pca, type=c("text", "points"), cex=1.5)
ordihull(pca, group=complete.df$SPECIESID, lwd=2,
         col=viridis::viridis_pal()(n_distinct(complete.df$GENUSID)))
legend("topleft", col=viridis::viridis_pal()(n_distinct(complete.df$GENUSID)),
       lty=1, lwd=2, levels(factor(complete.df$GENUSID)))
}
dev.off()

pca <- princomp(complete.df[,c(16,34:40,42:47)]) # 16 = v
pc <- prcomp(complete.df[,c(16,34:40,42:47)])
km_scree <- map(1:20, ~sum(kmeans(pc$x[,1:2], ., iter.max=1e3, nstart=25)$withinss))
plot(1:length(km_scree), km_scree, type="b")

km <- kmeans(pc$x[,1:2], 20, iter.max=1e4, nstart=50)
complete.df$cluster <- km$cluster
complete.df$PC1 <- pc$x[,1]
complete.df$PC2 <- pc$x[,2]
complete.df$PC3 <- pc$x[,3]

plot(pc$x[,1], pc$x[,2], col=km$cluster)

ggplot(complete.df, aes(x=cluster, y=SPECIESID)) + geom_point(alpha=0.2, size=4)
ggplot(complete.df, aes(PC1, PC2, colour=GENUSID)) + geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + facet_wrap(~cluster)
ggplot(complete.df, aes(PC1, PC2, colour=SPECIESID)) + geom_point() + 
  scale_colour_viridis_d() + facet_wrap(~cluster)
ggplot(complete.df, aes(PC1, PC2, colour=as.factor(cluster))) + geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + facet_wrap(~SPECIESID)
ggplot(complete.df, aes(PC1, PC2, colour=as.factor(cluster))) + geom_point() + 
  scale_colour_viridis_d() + facet_wrap(~SPECIESID)

round(cor(complete.df[c(22:26,28,16,30:43)])[,1:6],3)

ggplot(complete.df, aes(mnt25, PC1, colour=as.factor(cluster))) + 
  stat_smooth(method="loess", se=F) +
  geom_point() + scale_colour_viridis_d()# + facet_wrap(~SPECIESID)
ggplot(complete.df, aes(AP, PC1, group=SPECIESID)) + 
  stat_smooth(method="lm", se=F) +
  geom_point() #+ facet_wrap(~SPECIESID)
ggplot(complete.df, aes(GDD0, PC2, group=SPECIESID)) + 
  stat_smooth(method="lm", se=F) +
  geom_point() #+ scale_colour_viridis_d() #+ facet_wrap(~SPECIESID)






# Variance across levels
ggplot(trts$clny.wide, aes(y=SPECIESID, x=v_var)) + geom_boxplot()
ggplot(trts$clny.df, aes(y=Trait, x=sdValue/mnValue)) + geom_boxplot() + 
  facet_wrap(~SPECIESID)

pairs(trts$clny.wide[,c(5,86:96)], lower.panel=panel.smooth, cex=0.5)

# Are some traits more variable than others?
ggplot(trts$clny.df, aes(y=Trait, x=100*sdValue/mnValue)) + geom_boxplot()
anova(lm((100*sdValue/mnValue) ~ Trait, data=trts$clny.df))
cv_tukey <- TukeyHSD(aov(lm((100*sdValue/mnValue) ~ Trait, data=trts$clny.df)))
cv_tukey$Trait[cv_tukey$Trait[,4]<0.05,]

# How much variation occurs within colonies vs. within species?
full_join(trts$clny.df %>% group_by(SPECIESID, Trait) %>%
            summarise(cv_within_clny=mean(sdValue/mnValue, na.rm=T)*100,
                      var_within_clny=mean(sdValue^2, na.rm=T),
                      n_col=n()),
          trts$wkr.df %>% group_by(SPECIESID, Trait) %>%
            summarise(cv_within_spp=sd(Value, na.rm=T)/mean(Value, na.rm=T)*100,
                      var_within_spp=var(Value, na.rm=T)),
          by=c("SPECIESID", "Trait")) %>%
  filter(n_col > 3) %>%
  mutate(cv_clny_div_spp=cv_within_clny/cv_within_spp,
         var_clny_div_spp=var_within_clny/var_within_spp) %>%
  group_by(Trait) %>% summarise(mn=mean(var_clny_div_spp, na.rm=T)) %>%
  arrange(mn)
  # ggplot(aes(x=var_clny_div_spp, y=SPECIESID)) + geom_point() + facet_wrap(~Trait)


full_join(trts$clny.wide %>% group_by(SPECIESID) %>%
            summarise(cv_within_clny=mean(sqrt(v_var)/v, na.rm=T)*100,
                      var_within_clny=mean(v_var, na.rm=T),
                      n_col=n()),
          trts$wkr.wide %>% group_by(SPECIESID) %>%
            summarise(cv_within_spp=sd(v, na.rm=T)/mean(v, na.rm=T)*100,
                      var_within_spp=var(v, na.rm=T)),
          by=c("SPECIESID")) %>%
  filter(n_col > 3) %>%
  mutate(cv_clny_div_spp=cv_within_clny/cv_within_spp,
         var_clny_div_spp=var_within_clny/var_within_spp) %>%
  summarise(mn=mean(var_clny_div_spp, na.rm=T))
# ggplot(aes(x=var_clny_div_spp, y=SPECIESID)) + geom_point() + facet_wrap(~Trait)




nora.spp <- read_csv("~/Documents/unil/teaching/2020_Khelidj/code/RE__DPCOA/Species_fil.csv")
gen_tree <- ape::read.tree("~/Desktop/Formicidae_genus.nwk")
gen_tree.nora <- ape::drop.tip(gen_tree, gen_tree$tip.label[-match(unique(nora.spp$Genus), gen_tree$tip.label, nomatch=0)])
ape::cophenetic.phylo(gen_tree.nora)
