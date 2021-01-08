


# photoshop export exploration
library(tidyverse); library(readxl); library(googlesheets); theme_set(theme_bw())
library(viridis); library(foreach); library(doSNOW)
source("~/Documents/unil/opfo_str_sampling/code/lc_cols.R")
source("code/fn_aux.R")
pinned.df <- gs_read(ss=gs_title("Operation Fourmis pinned specimens"), ws="Sheet1") %>%
  mutate(SPECIESID=paste(str_sub(Genus, 1, 4), str_sub(Species, 1, 4), sep="_"))
tube_inc.df <- read_csv("data/tubes_trait_measurements.csv") %>%
  mutate(TubePinned=TubeNo %in% pinned.df$TubeNo,
         SppPinned=SPECIESID %in% pinned.df$SPECIESID)



ch_tanja <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer") %>%
  rename(TubeNo=CATALOGUENUMBER) %>%
  mutate(SampleDate=lubridate::dmy(SampleDate))

img_dir <- "data/trait_test/"
tube.dirs <- dir(img_dir, "^999", include.dirs=T, recursive=T)
# tube.dirs <- tube.dirs[str_split_fixed(tube.dirs, "/", 3)[,3] %in%
#                          tube_inc.df$TubeNo]

prl_c <- makeCluster(4); registerDoSNOW(prl_c)
tube.hist <- foreach(tube=seq_along(tube.dirs), .errorhandling="remove",
                     .packages=c("tidyverse")) %dopar% {
  tubeNo <- tube.dirs[tube]
  nWorkers <- max(length(dir(paste0(img_dir, tubeNo), "dor[1-9].txt")),
                  length(dir(paste0(img_dir, tubeNo), "meso[1-9].txt")))
  dor.hist <- meso.hist <- vector("list", length=nWorkers)
  for(i in 1:nWorkers) {
    meso_i.f <- dir(paste0(img_dir, tubeNo, "/meso", i, " Data/"), full.names=T)
    if(length(meso_i.f) > 0) {
      meso.hist[[i]] <- map(meso_i.f, ~t(as.matrix(read.csv(., header=F)))) %>%
        Reduce('+', .) %>%
        as.data.frame(.) %>% setNames(., c("R", "G", "B")) %>%
        mutate(value=1:256, 
               TubeNo=str_split_fixed(tube.dirs[tube], "/", 3)[3], 
               worker=i) 
    }
  }
  if(is.null(meso.hist[[1]])) {
    rbind(do.call('rbind', dor.hist) %>% mutate(img="dor"))
  } else if (is.null(dor.hist[[1]])) {
    rbind(do.call('rbind', meso.hist) %>% mutate(img="meso"))
  } else {
    rbind(do.call('rbind', dor.hist) %>% mutate(img="dor"), 
          do.call('rbind', meso.hist) %>% mutate(img="meso")) 
  }
}
stopCluster(prl_c)


all.hist <- do.call('rbind', tube.hist) %>%
  mutate(tubeWorker=paste0(TubeNo, "_", worker)) %>%
  gather(channel, count, 1:3) %>%
  group_by(tubeWorker, img) %>% mutate(prop=count/sum(count)) %>%
  ungroup() %>% left_join(., ch_tanja, by="TubeNo")
all.sum <- all.hist %>% group_by(TubeNo, worker, tubeWorker, channel, img) %>%
  summarise(mnValue=sum(value*count)/sum(count),
            medValue=median(rep(value, count)), 
            sdValue=sd(rep(value, count))) %>%
  left_join(., ch_tanja, by="TubeNo")
all.sum.wide <- all.sum %>% ungroup %>% select(-c(6,8)) %>% 
  spread(channel, medValue) %>%
  mutate(rgbCol=paste0("rgb(", R/256, ",", G/256, ",", B/256, ")"),
         h=rgb2hsv(R, G, B, maxColorValue=255)[1,],
         s=rgb2hsv(R, G, B, maxColorValue=255)[2,],
         v=rgb2hsv(R, G, B, maxColorValue=255)[3,],
         tubeCt=as.numeric(factor(TubeNo)))
colony.mean <- all.sum.wide %>% group_by(TubeNo, img) %>%
  summarise(R=mean(R), G=mean(G), B=mean(B), V_var=var(v), V_CV=sd(v)/mean(v)) %>%
  mutate(rgbCol=paste0("rgb(", R/256, ",", G/256, ",", B/256, ")"),
         h=rgb2hsv(R, G, B, maxColorValue=255)[1,],
         s=rgb2hsv(R, G, B, maxColorValue=255)[2,],
         v=rgb2hsv(R, G, B, maxColorValue=255)[3,]) %>%
  left_join(., ch_tanja, by="TubeNo")
  
img_focus <- c("meso", "dor")[1]
par(mfrow=c(1,1)) 
with(filter(all.sum.wide, img==img_focus), 
     {plot(mnt25, v, pch=19, cex=1, main="All species",
           xlab="Elevation (m)", ylab="Brightness (V)",
           col=rgb(R, G, B, maxColorValue=255))
       abline(lm(v ~ mnt25))}
)
with(filter(colony.mean, img==img_focus), 
     points(mnt25, v, pch=1, cex=2, col=rgb(R, G, B, maxColorValue=255)))

par(mfrow=c(3,3))
for(sp in unique(all.sum.wide$SPECIESID)) {
  # full color
  with(filter(all.sum.wide, SPECIESID==sp & img==img_focus), 
       {plot(mnt25, v, pch=19, cex=2, main=sp,
             xlab="Elevation (m)", ylab="Brightness (V)",
             col=rgb(R, G, B, maxColorValue=255))
         if(n_distinct(mnt25) > 1) abline(lm(v ~ mnt25))}
       )
  with(filter(colony.mean, SPECIESID==sp & img==img_focus), 
       points(mnt25, v, pch=1, cex=3, col=rgb(R, G, B, maxColorValue=255)))
}





ggplot(all.sum.wide, aes(x=bio1_tmean_8110/100, y=v)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="bio1: Mean Annual Temperature (ºC)", y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(all.sum.wide, aes(x=mnt25, y=v)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="Elevation (m)", y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2) 
ggplot(all.sum.wide, aes(x=SoilTemp, y=v)) + 
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="Soil Temperature (ºC)", y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(all.sum.wide, aes(x=bio1_tmean_8110/100+SoilTempAnomaly, y=v)) + 
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="MAT + Site-level Soil Temperature Anomaly (ºC)", 
       y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(all.sum.wide, aes(x=bio7_tar_8110/100, y=v)) + 
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="bio7: Temperature Annual Range (ºC)", 
       y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(all.sum.wide, aes(x=SampleDate, y=v)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", linetype=2, fill="gray85") +
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  geom_point(data=colony.mean, size=3, alpha=0.5, aes(colour=SPECIESID)) + 
  labs(x="Sample date", y="Brightness (mesosoma dorsum)") +
  facet_wrap(~SPECIESID, scales="free_y") +
  scale_colour_brewer(type="qual", palette=2)


ggplot(colony.mean, aes(x=bio1_tmean_8110, y=V_var))  + 
  stat_smooth(aes(colour=SPECIESID), method="lm", formula=y~x+I(x^2), 
              linetype=2, fill="gray85") +
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  labs(x="bio1: Mean Annual Temperature (ºC)",
       y="Intra-colony brightness variance") +
  facet_wrap(~SPECIESID, scales="free") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(colony.mean, aes(x=bio7_tar_8110, y=V_var))  + 
  stat_smooth(method="lm", colour="black") +
  geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", se=F, linetype=2) +
  labs(x="bio7: Temperature Annual Range (ºC)", 
       y="Intra-colony brightness variance") +
  facet_grid(.~img) + scale_colour_brewer(type="qual", palette=2)
ggplot(colony.mean, aes(x=bio1_tmean_8110, y=V_var))  + 
  stat_smooth(method="lm", colour="black") +
  # geom_point(alpha=0.7, size=2, shape=1, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", se=F, linetype=2) +
  labs(x="bio1: Mean Annual Temperature (ºC)", 
       y="Intra-colony brightness variance") +
  facet_grid(.~img) + scale_colour_brewer(type="qual", palette=2)
ggplot(colony.mean, aes(x=bio1_tmean_8110, y=v))  + 
  stat_smooth(method="lm", colour="black") +
  geom_point(alpha=0.7, size=2, aes(colour=SPECIESID)) + 
  stat_smooth(aes(colour=SPECIESID), method="lm", se=F, linetype=2) +
  labs(x="bio1: Mean Annual Temperature (ºC)", 
       y="Mean colony brightness") +
  scale_colour_brewer(type="qual", palette=2)

ggplot(colony.mean, aes(x=bio1_tmean_8110, y=V_CV))  + 
  stat_smooth(method="lm", colour="black", fill="gray85") +
  geom_point(alpha=0.7, size=2, aes(colour=SPECIESID)) +
  stat_smooth(aes(colour=SPECIESID), method="lm", se=F, linetype=2) +
  labs(x="bio1: Mean Annual Temperature (ºC)", 
       y="Intra-colony brightness CV") +
  scale_colour_brewer(type="qual", palette=2)
ggplot(colony.mean, aes(x=mnt25, y=V_CV))  + 
  stat_smooth(method="lm", colour="black", fill="gray85") +
  geom_point(alpha=0.7, size=2, aes(colour=SPECIESID)) +
  stat_smooth(aes(colour=SPECIESID), method="lm", se=F, linetype=2, size=0.5) +
  labs(x="Elevation (m)", 
       y="Intra-colony brightness CV") +
  scale_colour_brewer(type="qual", palette=2)




ggplot(colony.mean, aes(x=v, y=V_var, colour=mnt25)) + 
  geom_point() + facet_wrap(~SPECIESID, scales="free") + 
  scale_colour_viridis()



all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(h_std=(h-mean(h))/sd(h)) %>%
  ggplot(aes(x=mnt25, y=h_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)
all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(s_std=(s-mean(s))/sd(s)) %>%
  ggplot(aes(x=mnt25, y=s_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)
all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(v_std=(v-mean(v))/sd(v)) %>%
  ggplot(aes(x=mnt25, y=v_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)

all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(R_std=(R-mean(R))/sd(R)) %>%
  ggplot(aes(x=mnt25, y=R_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)
all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(G_std=(G-mean(G))/sd(G)) %>%
  ggplot(aes(x=mnt25, y=G_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)
all.sum.wide %>% group_by(SPECIESID) %>%
  mutate(B_std=(B-mean(B))/sd(B)) %>%
  ggplot(aes(x=mnt25, y=B_std)) + geom_point(alpha=0.5) + 
  stat_smooth(method="lm", col="gray30") + 
  stat_smooth(aes(group=SPECIESID), method="lm", size=0.5, linetype=2, se=F)


par(mar=c(6,6,4,4))
all.sum.wide %>% group_by(SPECIESID) %>% 
  mutate(mnt25_std=(mnt25-min(mnt25))/(max(mnt25)-min(mnt25))) %>%
  plot_colors_by_v(., "mnt25_std", "Location in elevational range")
plot_colors_by_v(all.sum.wide, "mnt25", "Elevation (m)")
plot_colors_by_v(all.sum.wide, "bio1_tmean_8110", "Mean Annual Temperature (ºC)")
plot_colors_by_v(all.sum.wide, "SoilTemp", "Soil Temperature (ºC)")





library(lme4)
mod <- lmer(I(v*100) ~ I(bio1_tmean_8110/100) + (1|SPECIESID) + (1|TubeNo), 
            data=all.sum.wide)
summary(mod)
confint(mod)
plot(mod)

summary(lm(I(v*100) ~ I(bio1_tmean_8110/100) + SPECIESID, data=all.sum.wide))
mod <- lmer(I(v*100) ~ mnt25 + (1|SPECIESID) + (1|TubeNo), data=all.sum.wide)
summary(mod)
confint(mod)
plot(mod)

mod <- lmer(V_CV*100 ~ mnt25 + (1|SPECIESID), data=colony.mean)
summary(mod)
confint(mod)




tube_inc.df %>% filter(grepl("Lasi_", SPECIESID)) %>%
  ggplot(aes(mnt25, forcats::fct_rev(factor(SPECIESID)), 
             colour=source, shape=TubePinned)) + 
  geom_point(alpha=0.5, size=3) + scale_shape_manual(values=c(1,19)) +
  scale_colour_manual(values=c("gray30","red")) + labs(x="Elevation (m)", y="")
tube_inc.df %>% group_by(SPECIESID, source) %>% 
  summarise(nTube=n()) %>% 
  spread(source, nTube) %>%
  filter(!is.na(structured)) %>%
  print.AsIs()
tube_inc.df %>% filter(SPECIESID=="Form_cuni" & source=="structured") %>%
  arrange(TubeNo) %>% print.AsIs()





