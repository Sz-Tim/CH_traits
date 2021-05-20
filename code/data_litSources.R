

# Fill out by measuring images on AntWeb? Why hasn't anyone done that?  


library(tidyverse); library(googlesheets); theme_set(theme_bw())
source("code/00_fn.R")
walk(paste0("../1_opfo/code/", c("lc_cols", "00_fn"), ".R"), source)
tax_i <- read_csv("../3_diversity/data/tax_i.csv") %>%
  mutate(Binomial=paste(FullGen, FullSpp))


# temp for merging species lists

ch_sp <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Species") %>%
  filter(!is.na(DESCRIPTEUR)) %>%
  rename(Subfamily=SOUSFAMILLE, Genus=GENRE, 
         Species=ESPECE, Binomial=GENREESPECE) %>%
  select(Subfamily, Genus, Species, Binomial) %>%
  mutate(Tanja_CH=1)
ch_obs <- load_ant_data(clean_spp=T)$all %>% st_set_geometry(NULL) %>%
  group_by(SPECIESID) %>% summarise(N=n()) %>%
  ungroup %>% mutate(Binomial=str_replace(SPECIESID, "_", " "))


ch_sp <- tax_i %>% select(Binomial)

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


## Ant Profiler
profiler <- read_tsv("data/lit/ANTPROFILER.txt") %>%
  mutate(Binomial=str_replace(Species, "_", " "))
prof_ch <- tax_i$Binomial %in% unique(profiler$Binomial)
sum(prof_ch)  # 28
ch_sp$prof <- as.numeric(prof_ch)




profiler %>% filter(Binomial %in% tax_i$Binomial) %>%
  filter(!is.na(WorkerPolymorphism)) %>% arrange(Binomial) %>%
  select(Binomial, WorkerPolymorphism) %>% print.AsIs





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
