# CH ant trait biogeography
# Trait measurements and color processing


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
trts <- load_traits(ant_i=ant.ls$str, msr_dir=msr_dir, col_dir=col_dir, 
                    na.thresh=0.05, lat_names=trait_names$lat,
                    fro_names=trait_names$fro, dor_names=trait_names$dor)



# Data QC
qc.done <- read_csv("msr_QC_checked.csv")
trait_name.df <- enframe(trait_names) %>% unnest(cols=c(value))

map(unique(trts$wkr.df$SPECIESID)[1:7],
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
      theme(panel.grid.major=element_line(colour="gray80"),
            panel.grid.minor=element_line(colour="gray90"),
            legend.position="none"))

qc.todo <- trts$wkr.std %>% 
  filter(abs(Value)>qnorm(0.995) & !TubeNo %in% qc.done$TubeNo) %>%
  mutate(img=paste0(trait_name.df$name[match(Trait, trait_name.df$value)],
                    Worker)) %>% 
  select(SPECIESID, TubeNo, img, Trait, Value) %>%
  arrange(SPECIESID, TubeNo, img, Trait)
qc.todo %>% print.AsIs
write_csv(qc.todo, "msr_QC_todo.csv")






# colors
plot_colors_by_v(filter(trts$wkr.df, !is.na(med_R)), 
                 "mnt25", "Elevation (m)")


# traits across environmental variables
ggplot(trts$clny.std, aes(mnt25, mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(data=trts$wkr.std, aes(y=Value), shape=1, alpha=0.7, size=0.7) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 
ggplot(trts$clny.std, aes(bio1_tmean_8110, mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(data=trts$wkr.std, aes(y=Value), shape=1, alpha=0.7, size=0.7) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) +
  stat_smooth(method="lm", formula=y~x+I(x^2), se=F, linetype=2, size=0.5)



# trait variability within colonies
ggplot(trts$clny.df, aes(mnt25, sdValue/mnValue)) + 
  facet_grid(SPECIESID~Trait, scales="free") +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm", se=F, linetype=2, size=0.5) 

