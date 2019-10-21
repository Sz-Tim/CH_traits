

# Fill out by measuring images on AntWeb? Why hasn't anyone done that?  


library(tidyverse); theme_set(theme_bw())

ch_sp <- read.csv("data/CH_spp_aline.csv")

########
## Regional studies with published databases
########

### Arnan et al 2014
arnan2014 <- readxl::read_xlsx("data/lit/arnan2014/arnan2014_S2.xlsx", 2) %>%
  rename(Binomial=Species)
arnan2014_ch <- ch_sp$Binomial %in% unique(arnan2014$Binomial)
sum(arnan2014_ch)  # 46
ch_sp$Arnan2014 <- as.numeric(arnan2014_ch)

### Arnan et al 2017
arnan2017 <- readxl::read_xlsx("data/lit/arnan2017/arnan2017.xlsx", 1)
arnan2017_ch <- ch_sp$Binomial %in% unique(arnan2017$Binomial)
sum(arnan2017_ch)  # 66
ch_sp$Arnan2017 <- as.numeric(arnan2017_ch)

### GlobalAnts
ga.df <- read_csv("data/lit/GlobalAnts/GlobalAnts.csv") %>%
  mutate(Binomial=paste(Genus, Species, sep=" "))
ga_ch <- ch_sp$Binomial %in% unique(ga.df$Binomial)
sum(ga_ch)  # 24
ch_sp$GlobalAnts <- as.numeric(ga_ch)
