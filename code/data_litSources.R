

# Fill out by measuring images on AntWeb? Why hasn't anyone done that?  


library(tidyverse); library(googlesheets); theme_set(theme_bw())


# temp for merging species lists

ch_tanja <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Species") %>%
  filter(!is.na(DESCRIPTEUR)) %>%
  rename(Subfamily=SOUSFAMILLE, Genus=GENRE, 
         Species=ESPECE, Binomial=GENREESPECE) %>%
  select(Subfamily, Genus, Species, Binomial) %>%
  mutate(Tanja_CH=1)
ch_aline <- read.csv("data/CH_spp_aline.csv") %>% 
  mutate(Aline_CH_List=1)
ch_aw <- read.csv("../../prep/CH_sp.csv") %>% 
  mutate(Binomial=str_trim(paste(Genus, Species, Subspecies)))
ch_both <- full_join(ch_aw, ch_aline, 
                     by=c("Subfamily", "Genus", "Species", "Binomial")) %>%
  full_join(., ch_tanja, by=c("Subfamily", "Genus", "Species", "Binomial"))
ch_both$AW_CH_List[is.na(ch_both$AW_CH_List)] <- 0
ch_both$AW_CH_specimen[is.na(ch_both$AW_CH_specimen)] <- 0
ch_both$AW_Photo[is.na(ch_both$AW_Photo)] <- 0
ch_both$Aline_CH_List[is.na(ch_both$Aline_CH_List)] <- 0
write_csv(ch_both, "data/CH_spp.csv")


ch_sp <- read.csv("data/CH_spp.csv")



########
## Identify overlapping species
########

### Arnan et al 2014
arnan2014 <- readxl::read_xlsx("data/lit/arnan2014/arnan2014_S2.xlsx", 2) %>%
  rename(Binomial=Species)
arnan2014_ch <- ch_sp$Binomial %in% unique(arnan2014$Binomial)
sum(arnan2014_ch)  # 47
ch_sp$Arnan2014 <- as.numeric(arnan2014_ch)

### Arnan et al 2017
arnan2017 <- readxl::read_xlsx("data/lit/arnan2017/arnan2017.xlsx", 1)
arnan2017_ch <- ch_sp$Binomial %in% unique(arnan2017$Binomial)
sum(arnan2017_ch)  # 67
ch_sp$Arnan2017 <- as.numeric(arnan2017_ch)

### GlobalAnts
ga.df <- read_csv("data/lit/GlobalAnts/GlobalAnts.csv") %>%
  mutate(Binomial=paste(Genus, Species, sep=" "))
ga_ch <- ch_sp$Binomial %in% unique(ga.df$Binomial)
sum(ga_ch)  # 25
ch_sp$GlobalAnts <- as.numeric(ga_ch)

## Burchill & Moreau 2016
bm2016 <- readxl::read_xls("data/lit/BurchillMoreau2016/BM_S1.xls", 1)
bm_ch <- ch_sp$Binomial %in% unique(bm2016$Binomial)
sum(bm_ch)  # 28
ch_sp$BurchillMoreau <- as.numeric(bm_ch)


## Seifert 2018
seifert <- readxl::read_xlsx("data/lit/seifert2018/seifert.xlsx", 1) %>%
  mutate(Binomial=paste(Genus, Species))






########
## Combine traits into one dataset
########

trait.ls <- list(Arnan2014=filter(arnan2014, Binomial %in% ch_sp$Binomial), 
                 Arnan2017=filter(arnan2017, Binomial %in% ch_sp$Binomial),  
                 GlobalAnts=filter(ga.df, Binomial %in% ch_sp$Binomial), 
                 BM2016=filter(bm2016, Binomial %in% ch_sp$Binomial))
size_col.df <- map_dfr(trait.ls, ~select(., one_of("Binomial", "WS", "CS")), 
                       .id="Source") %>%
  filter(!(is.na(WS) & is.na(CS)))
size_col.sum <- size_col.df %>% group_by(Binomial) %>%
  summarise(WS=mean(WS, na.rm=T),
            CS=mean(CS, na.rm=T)) %>%
  mutate(Genus=ch_sp$Genus[match(.$Binomial, ch_sp$Binomial)],
         Subfamily=ch_sp$Subfamily[match(.$Binomial, ch_sp$Binomial)])
  









sum(ch_sp$Arnan2014 | ch_sp$Arnan2017 | ch_sp$GlobalAnts | ch_sp$BurchillMoreau)
