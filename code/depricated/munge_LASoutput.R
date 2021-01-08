

library(tidyverse); library(readxl); library(googlesheets); 
theme_set(theme_bw() + theme(panel.grid=element_blank()))
library(viridis); library(GGally)
source("code/fn_aux.R")
source("~/Documents/unil/opfo_str_sampling/code/opfo_fn.R")
source("~/Documents/unil/ms_unil/CH_diversity/code/00_fn.R")
source("~/Documents/unil/opfo_str_sampling/code/lc_cols.R")

msr_dir <- "/Volumes/BOOTCAMP/Users/tsz/Desktop/opfo_images/"
col_dir <- "data/trait_test/"

trait_names <- list(lat=c("WebersLength", "HindTibia", "MidTibia"),
                    fro=c("HeadLength", "HeadWidth", 
                          "InterocularDistance", "ScapeLength"),
                    dor=c("MesosomaWidth", "MesosomaLength", 
                          "HindFemur", "MidFemur"))

ant.ls <- load_ant_data(pub.dir="data/", clean_spp=T)
traits_all <- load_traits(ant_i=ant.ls$str, msr_dir=msr_dir, col_dir=col_dir, 
                          na.thresh=0.05, lat_names=trait_names$lat,
                          fro_names=trait_names$fro, dor_names=trait_names$dor)
traits_all <- map(traits_all, ~filter(., SPECIESID!="Myrm_rugi"))
traits_all$wkr.std <- traits_all$wkr.df %>% 
  group_by(SPECIESID, Trait) %>%
  mutate(Value=(Value-mean(Value, na.rm=T))/sd(Value, na.rm=T), 
         v=(v-mean(v, na.rm=T))/sd(v, na.rm=T))
traits_all$clny.std <- traits_all$clny.df %>% 
  group_by(SPECIESID, Trait) %>%
  mutate(mnValue=(mnValue-mean(mnValue, na.rm=T))/sd(mnValue, na.rm=T), 
         v=(v-mean(v, na.rm=T))/sd(v, na.rm=T))
wide.wkr <- traits_all$wkr.df %>% 
  pivot_wider(names_from=Trait, values_from=Value) %>%
  mutate(MidLen=MidFemur+MidTibia, HindLen=HindFemur+HindTibia, 
         MesoSA=MesosomaLength*MesosomaWidth)
wide.clny <- traits_all$clny.df %>% 
  pivot_wider(names_from=Trait, values_from=c(mnValue, sdValue)) %>%
  mutate(MidLen=mnValue_MidFemur+mnValue_MidTibia, 
         HindLen=mnValue_HindFemur+mnValue_HindTibia, 
         MesoSA=mnValue_MesosomaLength*mnValue_MesosomaWidth)


plot_colors_by_v(filter(traits_all$wkr.df, !is.na(med_R)), 
                 "mnt25", "Elevation (m)")


ggplot(trts$clny.std, aes(mnt25, mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(data=trts$wkr.std, aes(y=Value), shape=1, alpha=0.7, size=0.7) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 
ggplot(traits_all$clny.std, aes(bio1_tmean_8110, mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(data=traits_all$wkr.std, aes(y=Value), shape=1, alpha=0.7, size=0.7) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 
ggplot(traits_all$clny.std, aes(SoilTemp, mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(data=traits_all$wkr.std, aes(y=Value), shape=1, alpha=0.7, size=0.7) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 


ggplot(traits_all$clny.df, aes(mnt25, sdValue/mnValue, colour=SPECIESID)) + 
  facet_wrap(~Trait, scales="free") +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 
ggplot(traits_all$clny.df, aes(bio1_tmean_8110, sdValue/mnValue, colour=SPECIESID)) + 
  facet_wrap(~Trait, scales="free") +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 
ggplot(traits_all$clny.df, aes(SoilTemp, sdValue/mnValue, colour=SPECIESID)) + 
  facet_wrap(~Trait, scales="free") +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 




ggplot(wide.wkr, aes(WebersLength, v, colour=SPECIESID)) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + 
  geom_point(data=wide.clny, aes(x=mnValue_WebersLength)) + 
  facet_wrap(~SPECIESID, scales="free")
ggplot(wide.clny, aes(Grass, HindLen/mnValue_WebersLength, colour=SPECIESID)) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1)
ggplot(wide.clny, aes(mnt25, mnValue_InterocularDistance/mnValue_HeadLength, 
                      colour=SPECIESID)) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 



ggplot(traits_all$wkr.df, aes(Value, colour=SPECIESID)) + 
  geom_density() + facet_wrap(~Trait, scales="free_y")




ggplot(traits_all$wkr.df, aes(mnt25, v)) + geom_point() + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  facet_wrap(~SPECIESID, scales="free_y")


ggplot(traits_all$clny.std, aes(mnValue, v)) + 
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  facet_grid(SPECIESID~Trait, scales="free")


ggplot(traits_all$wkr.df, aes(mnt25, Value, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")

ggplot(traits_all$clny.df, aes(mnt25, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")
ggplot(traits_all$clny.df, aes(mnt25, sdValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")

ggplot(traits_all$clny.df, aes(bio1_tmean_8110, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")
ggplot(traits_all$clny.df, aes(SoilTemp, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")
ggplot(traits_all$clny.df, aes(SoilTempAnomaly, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")

ggplot(traits_all$clny.df, aes(Grass, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")
ggplot(traits_all$clny.df, aes(Bare, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")

ggplot(traits_all$clny.df, aes(SampleDate, mnValue, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")



ggplot(traits_all$wkr.std, aes(SampleDate, Value, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")



ggplot(traits_all$wkr.std, aes(bio1_tmean_8110/100 + SoilTempAnomaly, 
                              Value, colour=SPECIESID)) + 
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1) + 
  facet_wrap(~Trait, scales="free_y")


ggplot(wide.clny, aes(mnt25, HindLen/mnValue_WebersLength, colour=SPECIESID)) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, shape=1)
ggplot(wide.wkr, aes(mnt25, WebersLength, colour=SPECIESID)) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) + 
  geom_point(size=2, alpha=0.5) + facet_wrap(~SPECIESID)
  



ggplot(traits_all$wkr.std, aes(Categorie, Value, fill=Categorie)) + 
  geom_boxplot() + facet_grid(SPECIESID~Trait) +
  scale_fill_manual(values=lc_cols) + coord_flip()
ggplot(traits_all$wkr.df, aes(bio1_tmean_8110, Value)) + 
  stat_smooth(method="lm") +
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")
ggplot(traits_all$wkr.df, aes(bio7_tar_8110, Value)) + 
  stat_smooth(method="lm") +
  geom_point(size=2, shape=1) + facet_wrap(~Trait, scales="free_y")



ggplot(traits_all$wkr.df, aes(SoilTempAnomaly, Value)) + 
  geom_point(aes(colour=bio1_tmean_8110), size=2) + scale_colour_viridis () +
  facet_wrap(~Trait, scales="free_y")




ggpairs(filter(wide.wkr, !is.na(HeadWidth)), 
        aes(colour=SPECIESID), columns=c(16, 80:90, 92:94),
        diag=list(continuous=wrap("densityDiag", alpha=0.5)),
        lower=list(continuous=wrap("smooth", se=F, alpha=0.5)))
ggcorr(wide.wkr[,c(4:12, 16, 44:66, 71:78, 80:94)])

ggpairs(filter(wide.clny, !is.na(mnValue_HeadWidth)), 
        aes(colour=SPECIESID), columns=c(6, 10, 75:85),
        lower=list(continuous=wrap("smooth", se=F, alpha=0.5)))
ggcorr(wide.clny[,c(6, 10, 38:60, 65:72, 75:85, 87:100)], size=3)



summary(lmerTest::lmer(mnValue_WebersLength ~ mnt25 + (1|SPECIESID), 
                       data=wide.clny))
summary(lmerTest::lmer(mnValue_WebersLength ~ I(bio1_tmean_8110/100) + (1|SPECIESID), 
                       data=wide.clny))
summary(lmerTest::lmer(mnValue_WebersLength ~ SoilTemp + (1|SPECIESID), 
                       data=wide.clny))

summary(lmerTest::lmer(mnValue_WebersLength ~ SampleDate + mnt25 + (1|SPECIESID), 
                       data=wide.clny))















# Data QC
qc.done <- read_csv("msr_QC_checked.csv")
trait_name.df <- enframe(trait_names) %>% unnest(cols=c(value))

map(unique(traits_all$wkr.df$SPECIESID)[1:6],
     ~ggplot(filter(traits_all$wkr.std, !is.na(Value) & SPECIESID==.), 
             aes(x=mnt25, y=Value)) + facet_wrap(~Trait) + 
      ggtitle(.) + geom_point(aes(colour=abs(Value)>qnorm(0.995) & 
                                    !TubeNo %in% qc.done$TubeNo)) + 
      geom_text(data=filter(traits_all$wkr.std, 
                            !is.na(Value) & 
                              SPECIESID==. &
                              abs(Value)>qnorm(0.995) &
                              !TubeNo %in% qc.done$TubeNo),
                aes(label=paste(TubeNo, Worker, sep="-")), size=2.5) +
      scale_colour_manual(values=c(`FALSE`="gray30", `TRUE`="red")) +
       theme(panel.grid.major=element_line(colour="gray80"),
             panel.grid.minor=element_line(colour="gray90"),
             legend.position="none"))


qc.todo <- traits_all$wkr.std %>% 
  filter(abs(Value)>qnorm(0.995) & !TubeNo %in% qc.done$TubeNo) %>%
  mutate(img=paste0(trait_name.df$name[match(Trait, trait_name.df$value)],
                    Worker)) %>% 
  select(SPECIESID, TubeNo, img, Trait, Value) %>%
  arrange(SPECIESID, TubeNo, img, Trait)
qc.todo %>% print.AsIs
write_csv(qc.todo, "msr_QC_todo.csv")
