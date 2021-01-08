# Simulating nested trait data


# Heavily modified, but originally copied from:
# https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

# devtools::install_github("rmcelreath/rethinking")
library(rethinking); library(tidyverse); library(rstan); library(mvtnorm)
library(readxl); library(googlesheets); library(bayesplot)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")
walk(paste0("../1_opfo/code/", c("lc_cols", "00_fn"), ".R"), source)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

lc_i <- readxl::read_xlsx("../1_opfo/data/landcover_id.xlsx", 1) %>%
  mutate(lcNum=as.numeric(LC_ID))

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
                      AP=paste0(gis_dir, "AP_chelsa", rast_end),
                      TwarmQ=paste0(gis_dir, "TwarmQ_chelsa", rast_end),
                      PwarmQ=paste0(gis_dir, "PwarmQ_chelsa", rast_end),
                      PcoldQ=paste0(gis_dir, "PcoldQ_chelsa", rast_end),
                      minTwarmest=paste0(gis_dir, "minTempWarmest_envirem", rast_end),
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


# model details
std <- F  # standardize (scale) within each species
clny_min <- 3
response_vars <- c("v", 
                   #"HeadWidth", "HeadLength", 
                   "HeadShape",
                   "DVE",
                   # "ScapeLength",
                   "WebersLength", 
                   "PronotExp",
                   "ScapeProp",
                   #"MesosomaWidth", #"MesosomaLength", "MesoSA",
                   # "MidTibia", "MidFemur",
                   # "MidLen",
                   # "HindTibia", "HindFemur", 
                   "RelLegHind")
X_vars_mn <- c(GDD="GDD0",
               AP="AP",
               # NPP="npp",
               # TAR="TAR",
               North="aspectN",
               # CnpyM="CnpyMixed",
               # CnpyO="CnpyOpen",
               # Elev="mnt25",
               Day="SampleDate"
               ) 
X_vars_sd <- c(#GDD="GDD0",
               AP="AP",
               # NPP="npp",
               TAR="TAR",
               North="aspectN",
               CnpyM="CnpyMixed",
               CnpyO="CnpyOpen",
               Elev="mnt25",
               Day="SampleDate"
               ) 
genera_incl <- c("Myrm", 
                 "Mani",
                 "Lasi",
                 "Temn",
                 "Tetr",
                 "Lept",
                 "Apha"
)

wkr.df <- trts$wkr.wide %>% 
  mutate(RelIntOc=(HeadWidth-InterocularDistance)/HeadLength,
         DVE=InterocularDistance/HeadWidth,
         RelLegHind=HindLen/WebersLength,
         HeadShape=HeadWidth/HeadLength,
         PronotExp=MesosomaWidth/WebersLength,
         ScapeProp=ScapeLength/HeadLength) %>%
  select(Worker, TubeNo, GENUSID, SPECIESID, SampleDate, mnt25,
         one_of(response_vars), one_of(X_vars_mn, X_vars_sd)) #%>%
# mutate(across(one_of(X_vars_mn, X_vars_sd), ~c(scale(.))))
if(std) {
  wkr.df <- wkr.df %>% group_by(SPECIESID) %>%
    mutate(across(one_of(response_vars, X_vars_mn, X_vars_sd), ~c(scale(.))))
  # wkr.df <- wkr.df %>% group_by(SPECIESID) %>%
  # mutate_if(is.numeric, ~c(scale(.)))
} else {
  wkr.df <- wkr.df %>% 
    mutate(across(one_of(response_vars, X_vars_mn, X_vars_sd), ~c(scale(.))))
  # wkr.df <- wkr.df %>% 
  #   mutate_if(is.numeric, ~c(scale(.)))
}
if("CnpyMixed" %in% c(X_vars_mn, X_vars_sd)) wkr.df$CnpyMixed <- trts$wkr.wide$CnpyMixed
if("CnpyOpen" %in% c(X_vars_mn, X_vars_sd)) wkr.df$CnpyOpen <- trts$wkr.wide$CnpyOpen
wkr.df <- wkr.df %>% 
  group_by(SPECIESID) %>% mutate(nColony=n_distinct(TubeNo)) %>% ungroup %>%
  filter(nColony >= clny_min) %>%
  filter(GENUSID %in% genera_incl) %>%
  filter(!is.na(CnpyOpen)) %>%
  arrange(SPECIESID, TubeNo) 


nSims <- 1e3
sim.ls <- priors.ls <- vector("list", nSims)

for(k in 1:length(response_vars)) {
  
  Y_var <- response_vars[k]; cat("----Starting", Y_var, "\n")
  
  # remove NAs and make list of data inputs
  if(any(is.na(wkr.df[[Y_var]]))) {
    wkr.df_i <- wkr.df[-which(is.na(wkr.df[[Y_var]])),]
  } else {
    wkr.df_i <- wkr.df
  }
  wkr.df_i <- wkr.df_i %>% 
    mutate(clny_id=as.numeric(factor(TubeNo, levels=unique(TubeNo))))
  clny.df_i <- wkr.df_i %>% 
    group_by(TubeNo, GENUSID, SPECIESID) %>% 
    summarise(across(where(is.numeric), list("mn"=mean, "sd"=sd)),
              nWkr=n()) %>%
    rename_with(~str_sub(., 1L, -4L), contains("_mn")) %>%
    ungroup %>% arrange(SPECIESID, TubeNo) %>%
    mutate(spp_id=as.numeric(factor(SPECIESID, levels=unique(SPECIESID))))
  
  d.ls <- list(N_clny=nrow(clny.df_i),
               N_wkr=nrow(wkr.df_i),
               S=n_distinct(wkr.df_i$SPECIESID),
               G=n_distinct(wkr.df_i$GENUSID),
               tax_i=cbind(as.numeric(factor(unique(wkr.df_i$SPECIESID))),
                           as.numeric(factor(str_sub(unique(wkr.df_i$SPECIESID),1,4)))),
               P_mn=length(X_vars_mn)+1,
               P_sd=length(X_vars_sd)+1,
               x_mn=cbind(int=1, as.matrix(clny.df_i[X_vars_mn])),
               x_sd=cbind(int=1, as.matrix(clny.df_i[X_vars_sd])),
               y=wkr.df_i[[Y_var]],
               clny_id=wkr.df_i$clny_id,
               spp_id=clny.df_i$spp_id)
  
  # sensitivity to wkr number
  clny_id_alt <- unlist(map(1:d.ls$N_clny, ~rep(.x, sample(1:3, 1))))
  d.ls$N_wkr <- length(clny_id_alt)
  d.ls$clny_id <- clny_id_alt
  
  for(s in 1:nSims) {
    # hyperpriors <- list(sigma_clny_1_global=truncnorm::rtruncnorm(1, 0, Inf, 3, 2),
    #                     sigma_clny_2_global=truncnorm::rtruncnorm(1, 0, Inf, 5, 2),
    #                     alpha=rnorm(d.ls$P_sd, 0, 0.5),
    #                     beta=rnorm(d.ls$P_mn, 0, 0.5),
    #                     sigma_A=abs(rnorm(d.ls$P_sd, 0, 0.1)),
    #                     sigma_B=abs(rnorm(d.ls$P_mn, 0, 0.2)),
    #                     sigma_a=abs(rnorm(d.ls$P_sd, 0, 0.1)),
    #                     sigma_b=abs(rnorm(d.ls$P_mn, 0, 0.2)))
    ppcheck <- prior_predictive_check(d.ls, hyperpriors)
    priors.ls[[s]] <- ppcheck$priors
    sim.ls[[s]] <- ppcheck$sims
  }
  
  sims <- list(y=unlist(map(sim.ls, ~.x$y)),
               y_bar=unlist(map(sim.ls, ~.x$y_bar)),
               mu=unlist(map(sim.ls, ~.x$mu)),
               d=unlist(map(sim.ls, ~.x$d)),
               delta=unlist(map(sim.ls, ~.x$delta)))
  priors <- list(sigma_clny_1=unlist(map(priors.ls, ~.x$sigma_clny_1)),
                 sigma_clny_2=unlist(map(priors.ls, ~.x$sigma_clny_2)),
                 sigma_clny_1_global=unlist(map(priors.ls, ~.x$sigma_clny_1_global)),
                 sigma_clny_2_global=unlist(map(priors.ls, ~.x$sigma_clny_2_global)))
  sims.rnk <- map(sims, ~percent_rank(.x))
  priors.rnk <- map(priors, ~percent_rank(.x))
  
  
  lo <- 0.025
  hi <- 0.975
  ybar_sample <- sample(which(between(sims.rnk$y_bar, lo, hi)), 500, F)
  d_sample <- sample(which(between(sims.rnk$d, lo, hi)), 500, F)
  wkr_sample <- sample(which(between(sims.rnk$y, lo, hi)), 500, F)
  map(sims, summary)
  par(mfrow=c(3,3))
  pwalk(list(a=sims, b=sims.rnk, c=names(sims)), 
        function(a, b, c) plot(density(a[which(between(b, lo, hi))]), main=c))
  plot(exp(sims$delta[d_sample]), sims$d[d_sample], col=rgb(0,0,0,0.2),
       xlab=expression(e^delta), ylab=expression(d[sim]))
  abline(a=0,b=1)
  plot(sims$mu[ybar_sample], sims$y_bar[ybar_sample], col=rgb(0,0,0,0.2),
       xlab=expression(mu), ylab=expression(bar(y)[sim]))
  abline(a=0,b=1)
  plot(density(priors$sigma_clny_1[between(priors.rnk$sigma_clny_1, lo, hi)]), 
       main=paste("sigma_1: ", round(mean(priors$sigma_clny_1_global, 2))))
  plot(density(priors$sigma_clny_2[between(priors.rnk$sigma_clny_2, lo, hi)]), 
       main=paste("sigma_2: ", round(mean(priors$sigma_clny_2_global, 2))))
  
  
  # circular simulation
  stan_d <- d.ls %>% list_modify(y=sim.ls[[s]]$y)
  out <- stan("code/mods/7_aRE_bRE.stan", data=stan_d, chains=3,
              warmup=1500, iter=2000, control=list(adapt_delta=0.95))
  
  wkr_out <- tibble(clny_id=d.ls$clny_id,
                    spp_id=d.ls$spp_id[clny_id],
                    y_obs=stan_d$y)
  clny_out <- wkr_out %>% group_by(spp_id, clny_id) %>%
    summarise(y_bar_obs=mean(y_obs), 
              d_obs=sd(y_obs)) %>% 
    bind_cols(map(c("y_bar", "d", "mu", "delta_exp"), 
                  ~summary_hdi(out, .x, 0.95) %>% 
                    set_names(paste0(.x, "_", names(.)))))
    
  # clny_out.8_12 <- clny_out
  # clny_out.1_4 <- clny_out
  # clny_out.1_3 <- clny_out
  
  # clny_out <- clny_out.8_12
  # clny_out <- clny_out.1_4
  
  ggplot(clny_out, aes(y_bar_obs, mu_mean, ymin=mu_L025, ymax=mu_L975)) + 
    geom_abline() + geom_point() + geom_linerange()
  ggplot(clny_out, aes(y_bar_obs, y_bar_mean, ymin=y_bar_L025, ymax=y_bar_L975)) + 
    geom_abline() + geom_point() + geom_linerange()
  ggplot(clny_out, aes(d_obs, delta_exp_mean, ymin=delta_exp_L025, ymax=delta_exp_L975)) + 
    geom_abline() + geom_point() + geom_linerange() + facet_wrap(~spp_id)
  ggplot(clny_out, aes(d_obs, d_mean, ymin=d_L025, ymax=d_L975)) + 
    geom_abline() + geom_point() + geom_linerange() + facet_wrap(~spp_id)
  
  out_alpha <- summary_hdi(out, "alpha") %>%
    mutate(true=hyperpriors$alpha)
  out_beta <- summary_hdi(out, "beta") %>%
    mutate(true=hyperpriors$beta)
  
  ggplot(out_alpha, aes(true, mean, ymin=L025, ymax=L975)) + 
    geom_abline() + geom_point() + geom_linerange()
  ggplot(out_beta, aes(true, mean, ymin=L025, ymax=L975)) + 
    geom_abline() + geom_point() + geom_linerange()
  
  
  
  hyperpriors <- list(sigma_clny_1_global=summary(out, "sigma_clny_1_global")$summary[,1],
                      sigma_clny_2_global=summary(out, "sigma_clny_2_global")$summary[,1],
                      alpha=summary(out, "alpha")$summary[,1],
                      beta=summary(out, "beta")$summary[,1],
                      sigma_A=summary(out, "sigma_A")$summary[,1],
                      sigma_B=summary(out, "sigma_B")$summary[,1],
                      sigma_a=summary(out, "sigma_a")$summary[,1],
                      sigma_b=summary(out, "sigma_b")$summary[,1])
  ppcheck <- prior_predictive_check(d.ls, hyperpriors)
  
  sim.df <- tibble(clny_id=d.ls$clny_id,
                   spp_id=d.ls$spp_id[clny_id],
                   y_obs=stan_d$y,
                   y=ppcheck$sims$y, 
                   y_bar=ppcheck$sims$y_bar[clny_id],
                   d=ppcheck$sims$d[clny_id],
                   mu=ppcheck$sims$mu[clny_id],
                   delta=ppcheck$sims$delta[clny_id])
  
  ggplot(sim.df, aes(x=y, y=y_obs)) + facet_wrap(~spp_id) + geom_abline() + geom_point()
  ggplot(sim.df, aes(x=exp(delta), y=d)) + facet_wrap(~spp_id) + geom_abline() + geom_point()
  
  
  for(s in 1:nSims) {
    
    stan_d <- d.ls %>% list_modify(y=sim.ls[[s]]$y)
    out <- stan("code/mods/7_aRE_bRE.stan", data=stan_d, chains=3,
                warmup=1500, iter=2000, control=list(adapt_delta=0.99))
    
    summary(out, pars="alpha")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=priors.ls[[s]]$alpha) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    summary(out, pars="beta")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=priors.ls[[s]]$beta) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    
    summary(out, pars="A")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=c(t(priors.ls[[s]]$A))) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    summary(out, pars="B")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=c(t(priors.ls[[s]]$B))) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    
    summary(out, pars="a")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=c(t(priors.ls[[s]]$a))) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    summary(out, pars="b")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=c(t(priors.ls[[s]]$b))) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    
    summary(out, pars="sigma_clny_1")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=priors.ls[[s]]$sigma_clny_1) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    summary(out, pars="sigma_clny_2")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=priors.ls[[s]]$sigma_clny_2) %>%
      ggplot(aes(x=rowname)) + geom_point(aes(y=mean)) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5) + 
      geom_point(aes(y=true), col="red", shape=1, size=3) + 
      labs(x="", y="value") + coord_flip()
    
    summary(out, pars="mu")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=sim.ls[[s]]$mu) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="fitted", title="mu")
    summary(out, pars="delta")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=sim.ls[[s]]$delta) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="fitted", title="delta")
    
    summary(out, pars="y_bar")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=sim.ls[[s]]$y_bar) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="fitted", title="y_bar")
    summary(out, pars="d_log")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=log(sim.ls[[s]]$d)) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="fitted", title="log(d)")
    
    summary(out, pars="y_bar_pred")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=sim.ls[[s]]$y_bar) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="predicted", title="y_bar")
    summary(out, pars="d_pred")$summary %>% 
      as.data.frame() %>% rownames_to_column() %>%
      mutate(true=sim.ls[[s]]$d) %>%
      ggplot(aes(true, mean)) + geom_abline(colour="gray") + 
      geom_point(alpha=0.5) +
      geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1, alpha=0.5) + 
      geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), size=0.5, alpha=0.5) + 
      labs(x="true", y="predicted", title="d") + scale_y_log10() + scale_x_log10()
  }
  
  
  
  
  
  
}


cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
plot(rep(1:d.ls$P_sd, each=d.ls$S), c(par.trans$a), 
     col=cols[rep(d.ls$tax_i[,2], times=d.ls$P_sd)])
plot(rep(1:d.ls$P_mn, each=d.ls$S), c(par.trans$b), 
     col=cols[rep(d.ls$tax_i[,2], times=d.ls$P_mn)])

out <- stan("code/mods/7_aRE_bRE.stan", data=stan_d, chains=3, 
            warmup=1500, iter=2000, control=list(adapt_delta=0.99))

mcmc_areas(out, pars=c(paste0("alpha[", 1:d.ls$P_sd, "]"),
                       paste0("beta[", 1:d.ls$P_mn, "]")))

post.beta <- as_tibble(rstan::extract(out, 
                                      pars=paste0("beta[", 1:d.ls$P_mn, "]"))) %>%
  rename_with(~paste0("beta_", 1:d.ls$P_mn))
for(i in 2:(d.ls$P_mn)) {
  p <- ggplot(post.beta) + xlim(-2, 2) + ylim(-2, 2) +
    geom_abline(aes_string(intercept="beta_1", slope=paste0("beta_", i)), 
                alpha=0.05, colour="cadetblue") +
    geom_abline(intercept=priors.ls[[s]]$beta[1], slope=priors.ls[[s]]$beta[i]) + 
    ggtitle(paste("beta", i))
  print(p)
}

post.alpha <- as_tibble(rstan::extract(out, 
                                      pars=paste0("alpha[", 1:d.ls$P_sd, "]"))) %>%
  rename_with(~paste0("alpha_", 1:d.ls$P_sd))
for(i in 2:(d.ls$P_sd)) {
  p <- ggplot(post.alpha) + xlim(-2, 2) + ylim(-2, 2) + 
    geom_abline(aes_string(intercept="alpha_1", slope=paste0("alpha_", i)), 
                alpha=0.05, colour="cadetblue") +
    geom_abline(intercept=priors.ls[[s]]$alpha[1], slope=priors.ls[[s]]$alpha[i]) + 
    ggtitle(paste("alpha", i))
  print(p)
}


post.y <- as_tibble(rstan::extract(out, pars="y_pred")[[1]], .name_repair="unique") %>%
  rename_with(~paste0("wkr_", 1:length(est.ls$y))) %>%
  pivot_longer(starts_with("wkr"), names_to="wkr", values_to="value") %>%
  mutate(wkr_id=as.numeric(str_sub(wkr, 5, -1L)),
         wkr_obs=est.ls$y[wkr_id],
         SPECIESID=wkr.df_i$SPECIESID[wkr_id],
         GENUSID=wkr.df_i$GENUSID[wkr_id]) %>%
  mutate(resid=value-wkr_obs) 
ggplot(post.y, aes(resid)) + geom_density() + facet_wrap(~SPECIESID, scales="free")
ggplot(post.y, aes(wkr_obs, resid, colour=factor(GENUSID))) + 
  geom_point(alpha=0.25) + facet_wrap(~SPECIESID, scales="free") + 
  scale_colour_brewer(palette=2, type="qual")

post.mu <- as_tibble(rstan::extract(out, pars="mu")[[1]], .name_repair="unique") %>%
  rename_with(~paste0("clny_", 1:nrow(sim_clny))) %>%
  pivot_longer(starts_with("clny"), names_to="clny", values_to="value") %>%
  mutate(clny_id=as.numeric(str_sub(clny, 6, -1L))) %>%
  full_join(., sim_clny, by="clny_id") %>% 
  mutate(resid=value-clny_mn,
         genus=tax_i[spp_id,2]) 
ggplot(post.mu[1:1e4,], aes(resid, group=clny_id)) +
  geom_vline(xintercept=0, linetype=3) +
  geom_line(stat="density", alpha=0.5) + facet_wrap(~spp_id)
  


out_b <- summary(out, pars="b")$summary[,1]
plot(c(b), out_b, col=cols[rep(tax_i[,2], times=P_mn)]); abline(0, 1, lty=3)

out_a <- summary(out, pars="a")$summary[,1]
plot(c(a), out_a, col=cols[rep(tax_i[,2], times=P_sd)]); abline(0, 1, lty=3)

out_alpha <- summary(out, pars="alpha")$summary[,1]
plot(alpha, out_alpha); abline(0, 1, lty=3)

out_beta <- summary(out, pars="beta")$summary[,1]
plot(beta, out_beta); abline(0, 1, lty=3)

out_sigma <- summary(out, pars="delta_exp")$summary[,1]
plot(sim_clny$delta, out_sigma); abline(0, 1, lty=3)

out_mu <- summary(out, pars="mu")$summary[,1]
plot(sim_clny$mu, out_mu); abline(0, 1, lty=3)


out.clny <- sim_clny %>%
  mutate(pred_mu=out_mu,
         pred_clny=summary(out, pars="clny_mn")$summary[,1],
         pred_sigma=out_sigma)

ggplot(out.clny, aes(pred_mu, mu)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(pred_clny, clny_mn)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(pred_sigma, delta)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)

ggplot(out.clny, aes(X2, pred_clny)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(X3, pred_clny)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(X4, pred_clny)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)

ggplot(out.clny, aes(X2, pred_sigma)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(X3, pred_sigma)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)
ggplot(out.clny, aes(X4, pred_sigma)) + geom_point() + 
  stat_smooth(method="lm", colour="black", size=0.3) + facet_wrap(~spp_id)

