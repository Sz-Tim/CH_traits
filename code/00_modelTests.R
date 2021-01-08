
library(tidyverse); library(readxl); library(googlesheets); library(rstan)
options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank(),
                             panel.grid.major=element_line(colour="gray90",
                                                           size=0.25)))

source("code/00_fn.R")
walk(paste0("../1_opfo/code/", c("lc_cols", "00_fn"), ".R"), source)
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
trts$spp_rng <- trts$clny.wide %>% group_by(SPECIESID) %>% 
  summarise(lo=min(mnt25), hi=max(mnt25), rng=hi-lo)



# model details
std <- F  # standardize (scale) within each species
clny_min <- 3
rng_thresh <- 500
response_vars <- c("v", 
                   "HeadWidth", "HeadLength",
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
X_vars_mn <- c(#GDD="GDD0",
               # Twarm="TwarmQ",
               minTwarm="minTwarmest",
               # AP="AP",
               # Pwarm="PwarmQ",
               # Pcold="PcoldQ",
               # NPP="npp",
               # TAR="TAR",
               # North="aspectN",
               # CnpyM="CnpyMixed",
               # CnpyO="CnpyOpen",
               # Elev="mnt25",
               Day="SampleDate") 
X_vars_sd <- c(#GDD="GDD0",
               # Twarm="TwarmQ",
               minTwarm="minTwarmest",
               # AP="AP",
               # Pwarm="PwarmQ",
               # Pcold="PcoldQ",
               # NPP="npp",
               # TAR="TAR",
               # North="aspectN",
               # CnpyM="CnpyMixed",
               # CnpyO="CnpyOpen",
               # Elev="mnt25",
               Day="SampleDate") 
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
  wkr.df <- wkr.df %>% #group_by(SPECIESID) %>%
    # mutate(across(one_of(response_vars, X_vars_mn, X_vars_sd), ~c(scale(.))))
    mutate(across(one_of(X_vars_mn, X_vars_sd), ~c(scale(.)))) %>%
    group_by(SPECIESID) %>%
    mutate(across(one_of(response_vars), ~c(scale(.)))) 
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
  filter(SPECIESID %in% filter(trts$spp_rng, rng >= rng_thresh)$SPECIESID) %>%
  filter(GENUSID %in% genera_incl) %>%
  # filter(!is.na(CnpyOpen)) %>%
  arrange(SPECIESID, TubeNo) 
write_csv(wkr.df, "data/wkr_df_temp.csv")




for(k in 1:length(response_vars)) {
  
  Y_var <- response_vars[k]; cat("----Starting", Y_var, "\n")
  
  # remove NAs and make list of data inputs
  if(any(is.na(wkr.df[[Y_var]]))) {
    wkr.df_i <- wkr.df[-which(is.na(wkr.df[[Y_var]])),]
  } else {
    wkr.df_i <- wkr.df
  }
  if(any(is.na(wkr.df_i[c(X_vars_mn, X_vars_sd)]))) {
    wkr.df_i <- wkr.df_i[-which(is.na(rowSums(wkr.df_i[c(X_vars_mn, X_vars_sd)]))),]
  }
  wkr.df_i <- wkr.df_i %>% 
    mutate(clny_id=as.numeric(factor(TubeNo, levels=unique(TubeNo))))
  clny.df_i <- wkr.df_i %>% 
    group_by(TubeNo, SPECIESID) %>% 
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
  
  
  # run models
  mod_dir <- "code/mods/convergence_trials"
  mods <- dir(mod_dir, ".stan", full.names=T, recursive=F)
  mods <- setNames(mods, str_sub(dir(mod_dir, ".stan", recursive=F), 3L, -6L))
  out <- map(mods, ~stan(., data=d.ls, chains=3, iter=3000, warmup=2000,
                         control=list(adapt_delta=0.9)))
  
  # compare models
  loo_comp <- loo::loo_compare(map(out, loo::loo)); loo_comp
  best <- rownames(loo_comp)[1]
  write.csv(as.data.frame(loo_comp),
            paste0("eda/loo_", Y_var, ifelse(std, "_std", ""), ".csv"))
  saveRDS(out[[best]], paste0("eda/0_out_", Y_var, ifelse(std, "_std", ""), ".rds"))
  f.i <- paste0(Y_var, ifelse(std, "_std", ""), ".pdf")

  
  # summarize output
  pars <- c("mu", "y_bar", "delta", "delta_exp", "d", "d_log")
  # pars <- c("y_bar", "d", "d_log")
  
  out.clny <- clny.df_i %>% ungroup %>%
    mutate(El_m=trts$clny.wide$mnt25[match(TubeNo, trts$clny.wide$TubeNo)]) %>%
    bind_cols(map(pars, ~summary_hdi(out[[best]], .x) %>%
                    set_names(paste0(.x, "_", names(.)))))
  out.wkr <- wkr.df_i %>% ungroup %>%
    full_join(select(out.clny, clny_id, all_of(paste0(pars, "_mean"))), 
              by="clny_id") %>%
    bind_cols(map(c("y_pred", "y_pred_"), ~summary_hdi(out[[best]], .x) %>%
                    set_names(paste0(.x, "_", names(.)))))
  out_alpha <- summary_hdi(out[[best]], "alpha", varNames=c("int", names(X_vars_sd)))
  out_beta <- summary_hdi(out[[best]], "beta", varNames=c("int", names(X_vars_mn)))


  
  
  
  
  # visualize
  ## parameter estimates
  ggplot(filter(out_alpha, !X_var %in% c("int")), 
         aes(x=X_var, y=median)) + 
    geom_hline(yintercept=0, colour="gray70") + geom_point(size=3) + 
    geom_linerange(aes(ymin=L25, ymax=L75), size=1.1) +
    geom_linerange(aes(ymin=L05, ymax=L95), size=0.5) +
    geom_errorbar(aes(ymin=L025, ymax=L975), size=0.15, width=0.1) +
    theme(panel.grid=element_blank(),
          axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
    labs(title=paste(Y_var, "sd"), x="", y="") #+ talk_fonts
  ggsave(paste0("eda/alpha_", f.i), width=3, height=6)
  ggplot(filter(out_beta, !X_var %in% c("int")), 
         aes(x=X_var, y=median)) + 
    geom_hline(yintercept=0, colour="gray70") + geom_point(size=3) + 
    geom_linerange(aes(ymin=L25, ymax=L75), size=1.1) +
    geom_linerange(aes(ymin=L05, ymax=L95), size=0.5) +
    geom_errorbar(aes(ymin=L025, ymax=L975), size=0.15, width=0.1) +
    theme(panel.grid=element_blank(),
          axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
    labs(title=paste(Y_var, "mean"), x="", y="") #+ talk_fonts
  ggsave(paste0("eda/beta_", f.i), width=3, height=6)
  
  if(any(grepl("^a_mn", names(out[[best]])))) {
    out_a <- summary_hdi(out[[best]], "a_mn")
    n_a_pars <- nrow(out_a)/d.ls$S
    out_a %>%
      mutate(spp=rep(unique(clny.df_i$SPECIESID), n_a_pars),
             par=rep(1:n_a_pars, each=d.ls$S), 
             parName=rep(colnames(d.ls$x_sd)[1:n_a_pars], each=d.ls$S)) %>%
      ggplot(aes(x=mean, y=spp)) + 
      geom_vline(xintercept=0, colour="gray90") + 
      geom_point(alpha=0.9, size=2) + 
      geom_linerange(aes(xmin=L25, xmax=L75), size=1.1) +
      geom_linerange(aes(xmin=L05, xmax=L95), size=0.5) +
      geom_errorbar(aes(xmin=L025, xmax=L975), size=0.15, width=0.1) +
      theme(panel.grid.major.x=element_blank()) + 
      facet_wrap(~parName, scales="free_x") +
      labs(x="Mean species slopes: Trait sd", y="", title=Y_var)
    ggsave(paste0("eda/a_", f.i), width=5.25, height=5)
  }
  if(any(grepl("^A_mn", names(out[[best]])))) {
    out_A <- summary_hdi(out[[best]], "A_mn")
    n_A_pars <- nrow(out_A)/d.ls$G
    out_A %>%
      mutate(gen=rep(unique(str_split_fixed(clny.df_i$SPECIESID, "_", 2)[,1]), n_A_pars),
             par=rep(1:n_A_pars, each=d.ls$G), 
             parName=rep(colnames(d.ls$x_sd)[1:n_A_pars], each=d.ls$G)) %>%
      ggplot(aes(x=mean, y=gen)) + 
      geom_vline(xintercept=0, colour="gray90") + 
      geom_point(alpha=0.9, size=2) + 
      geom_linerange(aes(xmin=L25, xmax=L75), size=1.1) +
      geom_linerange(aes(xmin=L05, xmax=L95), size=0.5) +
      geom_errorbar(aes(xmin=L025, xmax=L975), size=0.15, width=0.1) +
      theme(panel.grid.major.x=element_blank()) + 
      facet_wrap(~parName, scales="free_x") +
      labs(x="Mean genus slopes: Trait sd", y="", title=Y_var)
    ggsave(paste0("eda/Ag_", f.i), width=5.25, height=5)
  }
  if(any(grepl("^b_mn", names(out[[best]])))) {
    out_b <- summary_hdi(out[[best]], "b_mn")
    n_b_pars <- nrow(out_b)/d.ls$S
    out_b %>%
      mutate(spp=rep(unique(clny.df_i$SPECIESID), n_b_pars),
             par=rep(1:n_b_pars, each=d.ls$S), 
             parName=rep(colnames(d.ls$x_mn)[1:n_b_pars], each=d.ls$S)) %>%
      ggplot(aes(x=mean, y=spp)) + 
      geom_vline(xintercept=0, colour="gray90") + 
      geom_point(alpha=0.9, size=2) + 
      geom_linerange(aes(xmin=L25, xmax=L75), size=1.1) +
      geom_linerange(aes(xmin=L05, xmax=L95), size=0.5) +
      geom_errorbar(aes(xmin=L025, xmax=L975), size=0.15, width=0.1) +
      theme(panel.grid.major.x=element_blank()) + 
      facet_wrap(~parName, scales="free_x") +
      labs(x="Mean species slopes: Trait mean", y="", title=Y_var)
    ggsave(paste0("eda/b_", f.i), width=5.25, height=5)
  }
  if(any(grepl("^B_mn", names(out[[best]])))) {
    out_B <- summary_hdi(out[[best]], "B_mn")
    n_B_pars <- nrow(out_B)/d.ls$G
    out_B %>%
      mutate(gen=rep(unique(str_split_fixed(clny.df_i$SPECIESID, "_", 2)[,1]), n_B_pars),
             par=rep(1:n_B_pars, each=d.ls$G), 
             parName=rep(colnames(d.ls$x_mn)[1:n_B_pars], each=d.ls$G)) %>%
      ggplot(aes(x=mean, y=gen)) + 
      geom_vline(xintercept=0, colour="gray90") + 
      geom_point(alpha=0.9, size=2) + 
      geom_linerange(aes(xmin=L25, xmax=L75), size=1.1) +
      geom_linerange(aes(xmin=L05, xmax=L95), size=0.5) +
      geom_errorbar(aes(xmin=L025, xmax=L975), size=0.15, width=0.1) +
      theme(panel.grid.major.x=element_blank()) + 
      facet_wrap(~parName, scales="free_x") +
      labs(x="Mean genus slopes: Trait mn", y="", title=Y_var)
    ggsave(paste0("eda/Bg_", f.i), width=5.25, height=5)
  }
  
  ## predictions vs. observed values
  ggplot(out.clny, aes_string("mu_mean", Y_var)) + 
    geom_abline(linetype=3, size=0.5) + geom_point() + 
    stat_smooth(method="lm", formula=y~x, colour="black", size=0.3) +
    labs(x=expression(mu), y="Observed colony mean", title=Y_var)
  ggsave(paste0("eda/mod_obs_mu_", f.i), width=5.25, height=5)
  ggplot(out.clny, aes_string("y_bar_mean", Y_var)) + 
    geom_abline(linetype=3, size=0.5) + geom_point() + 
    stat_smooth(method="lm", formula=y~x, colour="black", size=0.3) +
    labs(x=expression(bar(y)), y="Observed colony mean", title=Y_var)
  ggsave(paste0("eda/mod_obs_y_bar_", f.i), width=5.25, height=5)
  ggplot(out.clny, aes_string("delta_mean", paste0("log(", Y_var, "_sd)"))) +
    geom_abline(linetype=3, size=0.5) + geom_point() + 
    stat_smooth(method="lm", formula=y~x, colour="black", size=0.3) +
    labs(x=expression(delta), y="Observed within-colony log(sd)", title=Y_var)
  ggsave(paste0("eda/mod_obs_delta_", f.i), width=5.25, height=5)
  ggplot(out.clny, aes_string("d_log_mean", paste0("log(", Y_var, "_sd)"))) +
    geom_abline(linetype=3, size=0.5) + geom_point() + 
    stat_smooth(method="lm", formula=y~x, colour="black", size=0.3) +
    labs(x=expression(log(d)), y="Observed within-colony log(sd)", title=Y_var)
  ggsave(paste0("eda/mod_obs_d_", f.i), width=5.25, height=5)
  
  
  
  ## across elevation
  ggplot(out.clny, aes(El_m, mu_mean)) + geom_point() + 
    geom_point(aes(y=y_bar_mean), shape=1, alpha=0.5) +
    {if(std) stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
                         colour="gray40", size=0.2, se=F) } +
    {if(std) stat_smooth(method="lm", formula=y~x, colour="gray40", se=F, size=1.5) } +
    {if(!std) stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
                          colour="gray40", size=0.3, se=F) } +
    labs(x="Elevation (m)", y="Mean", title=Y_var)
  ggsave(paste0("eda/el_mn_all_", f.i), width=5.25, height=5)
  ggplot(out.clny, aes(El_m, delta_mean)) + geom_point() + 
    geom_point(aes(y=d_log_mean), shape=1, alpha=0.5) +
    {if(std) stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
                         colour="gray40", size=0.2, se=F) } +
    {if(std) stat_smooth(method="lm", formula=y~x, colour="gray40", se=F, size=1.5) } +
    {if(!std) stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
                          colour="gray40", size=0.3, se=F) } +
    labs(x="Elevation (m)", y="Within-colony log sd", title=Y_var)
  ggsave(paste0("eda/el_sd_all_", f.i), width=5.25, height=5)
  
# 
#   mn_plots <- setNames(vector("list", length(X_vars_mn)), names(X_vars_mn))
#   sd_plots <- setNames(vector("list", length(X_vars_sd)), names(X_vars_sd))
#   for(j in 1:length(X_vars_mn)) {
#     i <- X_vars_mn[j]
#     mn_plots[[j]] <- ggplot(out.wkr, aes_string(i, "mu_mean")) +
#       geom_point(data=out.clny, colour="#2171b5") +
#       geom_point(aes_string(y=Y_var), shape=1, size=0.3, alpha=0.3) +
#       stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
#                   colour="gray40", size=0.3, se=F) +
#       labs(x=i, y="")
#   }
#   for(j in 1:length(X_vars_sd)) {
#     i <- X_vars_sd[j]
#     sd_plots[[j]] <- ggplot(out.clny, aes_string(i, "delta_mean")) +
#       geom_point(colour="#2171b5") +
#       geom_point(aes_string(y=paste0("log(", Y_var, "_sd)")), 
#                  shape=1, size=0.3, alpha=0.5) +
#       stat_smooth(aes(group=SPECIESID), method="lm", formula=y~x,
#                   colour="gray40", size=0.3, se=F) +
#       labs(x=i, y="") 
#   }
#   mn_p <- ggpubr::ggarrange(plotlist=mn_plots) %>%
#     ggpubr::annotate_figure(top=Y_var, left="Mean")
#   ggsave(paste0("eda/margEff_mn_", Y_var, ifelse(std, "_std", ""), ".pdf"),
#          mn_p, width=10, height=7)
#   sd_p <- ggpubr::ggarrange(plotlist=sd_plots) %>%
#     ggpubr::annotate_figure(top=Y_var, left="Within-colony sd")
#   ggsave(paste0("eda/margEff_sd_", Y_var, ifelse(std, "_std", ""), ".pdf"),
#          sd_p, width=10, height=7)
#   
#   wkr.sample <- sample(1:nrow(out.wkr), min(1000, nrow(out.wkr)), replace=F)
#   wkr.post <- map(wkr.sample, 
#                ~rstan::extract(out[[best]], pars=paste0("y_pred[", .x, "]")) %>%
#                  as_tibble %>% rename_at(1, ~c("y_post")))
#   wkr.post.df <- out.wkr[wkr.sample,] %>%
#     mutate(wkr_id=wkr.sample, post=wkr.post) %>% unnest(post)
#   wkr.post.sum <- wkr.post.df %>% group_by(wkr_id) %>%
#     rename_at(Y_var, ~c("obs")) %>%
#     summarise(obs=first(obs), SPECIESID=first(SPECIESID), GENUSID=first(GENUSID),
#               mnt25=first(mnt25),
#               y_pred=first(y_pred_mean), 
#               y_q025=quantile(y_post, probs=0.025),
#               y_q05=quantile(y_post, probs=0.05),
#               y_q25=quantile(y_post, probs=0.25),
#               y_q50=quantile(y_post, probs=0.5),
#               y_q75=quantile(y_post, probs=0.75),
#               y_q95=quantile(y_post, probs=0.95),
#               y_q975=quantile(y_post, probs=0.975))
#   
#   ggplot(wkr.post.sum, aes(obs, y=y_q50, colour=GENUSID)) + 
#     scale_colour_brewer("", type="qual", palette=2) +
#     geom_abline() + 
#     geom_point(alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q25, ymax=y_q75), size=1, alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q05, ymax=y_q95), size=0.5, alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q025, ymax=y_q975), size=0.2, alpha=0.5) + 
#     labs(x=expression(y[obs]), y=expression(y[mod]))
#   ggsave(paste0("eda/post_y_line_", f.i), width=6, height=5)
#   
#   
#   # clny.sample <- sample(1:nrow(out.clny), 25, replace=F)
#   clny.post.mu <- map(1:nrow(out.clny), 
#                   ~rstan::extract(out[[best]], pars=paste0("mu[", .x, "]")) %>%
#                     as_tibble %>% rename_at(1, ~c("mu_post")))
#   clny.post.y <- map(1:nrow(out.clny), 
#                      ~rstan::extract(out[[best]], pars=paste0("y_bar[", .x, "]")) %>%
#                        as_tibble %>% rename_at(1, ~c("y_bar_post")))
#   clny.post.delta <- map(1:nrow(out.clny), 
#                    ~rstan::extract(out[[best]], pars=paste0("delta[", .x, "]")) %>%
#                      as_tibble %>% rename_at(1, ~c("delta_post")))
#   clny.post.d <- map(1:nrow(out.clny), 
#                      ~rstan::extract(out[[best]], pars=paste0("d_log[", .x, "]")) %>%
#                        as_tibble %>% rename_at(1, ~c("d_post")))
#   clny.post.df <- out.clny[1:nrow(out.clny),] %>%
#     mutate(clny_id=1:nrow(out.clny), 
#            post.mu=clny.post.mu, post.y=clny.post.y, 
#            post.delta=clny.post.delta, post.d=clny.post.d) %>% 
#     unnest(c(post.mu, post.y, post.delta, post.d)) %>% 
#     mutate(GENUSID=str_split_fixed(SPECIESID, "_", 2)[,1]) 
#   clny.post.sum <- clny.post.df %>% 
#     group_by(clny_id) %>%
#     rename_at(Y_var, ~c("obs")) %>% rename_at(paste0(Y_var, "_sd"), ~c("obs_sd")) %>%
#     summarise(obs_y=first(obs), obs_d=log(first(obs_sd)),
#               SPECIESID=first(SPECIESID), GENUSID=first(GENUSID),
#               mnt25=first(mnt25),
#               y_pred=first(y_bar_mean), 
#               mu_q025=quantile(mu_post, probs=0.025),
#               mu_q05=quantile(mu_post, probs=0.05),
#               mu_q25=quantile(mu_post, probs=0.25),
#               mu_q50=quantile(mu_post, probs=0.5),
#               mu_q75=quantile(mu_post, probs=0.75),
#               mu_q95=quantile(mu_post, probs=0.95),
#               mu_q975=quantile(mu_post, probs=0.975),
#               y_q025=quantile(y_bar_post, probs=0.025),
#               y_q05=quantile(y_bar_post, probs=0.05),
#               y_q25=quantile(y_bar_post, probs=0.25),
#               y_q50=quantile(y_bar_post, probs=0.5),
#               y_q75=quantile(y_bar_post, probs=0.75),
#               y_q95=quantile(y_bar_post, probs=0.95),
#               y_q975=quantile(y_bar_post, probs=0.975),
#               d_pred=first(d_mean), 
#               delta_q025=quantile(delta_post, probs=0.025),
#               delta_q05=quantile(delta_post, probs=0.05),
#               delta_q25=quantile(delta_post, probs=0.25),
#               delta_q50=quantile(delta_post, probs=0.5),
#               delta_q75=quantile(delta_post, probs=0.75),
#               delta_q95=quantile(delta_post, probs=0.95),
#               delta_q975=quantile(delta_post, probs=0.975),
#               d_q025=quantile(d_post, probs=0.025),
#               d_q05=quantile(d_post, probs=0.05),
#               d_q25=quantile(d_post, probs=0.25),
#               d_q50=quantile(d_post, probs=0.5),
#               d_q75=quantile(d_post, probs=0.75),
#               d_q95=quantile(d_post, probs=0.95),
#               d_q975=quantile(d_post, probs=0.975))
#   
#   ggplot(clny.post.sum, aes(obs_y, y=mu_q50, colour=GENUSID)) + 
#     scale_colour_brewer("", type="qual", palette=2) +
#     geom_abline() + 
#     geom_point(alpha=0.5) + 
#     geom_linerange(aes(ymin=mu_q25, ymax=mu_q75), size=1, alpha=0.5) + 
#     geom_linerange(aes(ymin=mu_q05, ymax=mu_q95), size=0.5, alpha=0.5) + 
#     geom_linerange(aes(ymin=mu_q025, ymax=mu_q975), size=0.2, alpha=0.5) +
#     labs(x=expression(bar(y)[obs]), y=expression(mu[post]))
#   ggsave(paste0("eda/post_mu_line_", f.i), width=6, height=5)
#   ggplot(clny.post.sum, aes(obs_y, y=y_q50, colour=GENUSID)) + 
#     scale_colour_brewer("", type="qual", palette=2) +
#     geom_abline() + 
#     geom_point(alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q25, ymax=y_q75), size=1, alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q05, ymax=y_q95), size=0.5, alpha=0.5) + 
#     geom_linerange(aes(ymin=y_q025, ymax=y_q975), size=0.2, alpha=0.5) +
#     labs(x=expression(bar(y)[obs]), y=expression(bar(y)[post]))
#   ggsave(paste0("eda/post_y_bar_line_", f.i), width=6, height=5)
#   ggplot(clny.post.sum, aes(x=obs_d, y=delta_q50, colour=GENUSID)) + 
#     scale_colour_brewer("", type="qual", palette=2) +
#     geom_abline() + 
#     geom_point(alpha=0.5) +
#     geom_linerange(aes(ymin=delta_q25, ymax=delta_q75), size=1, alpha=0.5) +
#     geom_linerange(aes(ymin=delta_q05, ymax=delta_q95), size=0.5, alpha=0.5) +
#     geom_linerange(aes(ymin=delta_q025, ymax=delta_q975), size=0.2, alpha=0.5) +
#     labs(x=expression(log(d[obs])), y=expression(delta[post]))
#   ggsave(paste0("eda/post_delta_line_", f.i), width=6, height=5)
#   ggplot(clny.post.sum, aes(x=obs_d, y=d_q50, colour=GENUSID)) + 
#     scale_colour_brewer("", type="qual", palette=2) +
#     geom_abline() + 
#     geom_point(alpha=0.5) +
#     geom_linerange(aes(ymin=d_q25, ymax=d_q75), size=1, alpha=0.5) +
#     geom_linerange(aes(ymin=d_q05, ymax=d_q95), size=0.5, alpha=0.5) +
#     geom_linerange(aes(ymin=d_q025, ymax=d_q975), size=0.2, alpha=0.5) +
#     labs(x=expression(log(d[obs])), y=expression(log(d[post])))
#   ggsave(paste0("eda/post_d_line_", f.i), width=6, height=5)
#   
#   
  
  
}





p_y <- p_d <- resid_y <- resid_d <- vector("list", length(response_vars))
R2.df <- data.frame(Trait=response_vars, R2_mn=NA, R2_sd=NA)

for(k in 1:length(response_vars)) {
  Y_var <- response_vars[k]
  
  out <- readRDS(paste0("eda/0_out_", Y_var, ifelse(std, "_std", ""), ".rds"))
  
  if(any(is.na(wkr.df[[Y_var]]))) {
    wkr.df_i <- wkr.df[-which(is.na(wkr.df[[Y_var]])),]
  } else {
    wkr.df_i <- wkr.df
  }
  clny.df_i <- wkr.df_i %>% 
    group_by(TubeNo, SPECIESID) %>% 
    summarise(across(where(is.numeric), list("mn"=mean, "sd"=sd))) %>%
    rename_with(~str_sub(., 1L, -4L), contains("_mn")) %>%
    ungroup %>% arrange(SPECIESID, TubeNo)
  
  # summarize output
  probs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  prob_names <- paste0("q", str_sub(as.character(probs), 3, -1))
  pars <- c("mu", "y_bar", "delta", "delta_exp", "d", "d_log")
  out.clny <- clny.df_i %>% ungroup %>%
    mutate(El_m=trts$clny.wide$mnt25[match(TubeNo, trts$clny.wide$TubeNo)]) %>%
    bind_cols(map_dfc(setNames(pars, pars), ~summary(out, pars=.x)$summary[,1]))
  
  SS_tot_mn <- sum((out.clny[[Y_var]] - mean(out.clny[[Y_var]]))^2)
  SS_resid_mn <- sum((out.clny$mu - out.clny[[Y_var]])^2)
  R2_mn <- 1 - SS_resid_mn/SS_tot_mn
  
  SS_tot_sd <- sum((out.clny[[paste0(Y_var, "_sd")]] - 
                      mean(out.clny[[paste0(Y_var, "_sd")]], na.rm=T))^2, na.rm=T)
  SS_resid_sd <- sum((out.clny$delta_exp - 
                        out.clny[[paste0(Y_var, "_sd")]])^2, na.rm=T)
  R2_sd <- 1 - SS_resid_sd/SS_tot_sd
  
  R2.df$R2_mn[k] <- R2_mn
  R2.df$R2_sd[k] <- R2_sd
  
  p_y[[k]] <- ggplot(out.clny, aes(mu, group=SPECIESID)) + 
    geom_abline(colour="gray") + 
    geom_point(aes_string(y=response_vars[k])) + 
    stat_smooth(aes_string(y=response_vars[k]), method="lm", 
                formula=y~x, se=F, colour="cadetblue", size=0.75) +
    theme(panel.grid=element_blank()) + facet_wrap(~SPECIESID) +
    labs(x=expression(mu), y="Observed mean", title=response_vars[k])
  p_d[[k]] <- ggplot(out.clny, aes(delta_exp, group=SPECIESID)) + 
    geom_abline(colour="gray") + 
    geom_point(aes_string(y=paste0(response_vars[k], "_sd"))) + 
    stat_smooth(aes_string(y=paste0(response_vars[k], "_sd")), method="lm", 
                formula=y~x, se=F, colour="cadetblue", size=0.75) +
    theme(panel.grid=element_blank()) + facet_wrap(~SPECIESID) + 
    labs(x=expression(e^delta), y="Observed sd", title=response_vars[k])
  resid_y[[k]] <- ggplot(out.clny, aes(mu)) + 
    geom_hline(yintercept=0, colour="gray") + 
    geom_point(aes_string(y=paste0(response_vars[k], "- mu"))) + 
    theme(panel.grid=element_blank()) + #facet_wrap(~SPECIESID) +
    labs(x=expression(mu), y="Residual error", title=response_vars[k])
  resid_d[[k]] <- ggplot(out.clny, aes(delta)) + 
    geom_hline(yintercept=0, colour="gray") + 
    geom_point(aes_string(y=paste0("log(", response_vars[k], "_sd) - delta"))) + 
    theme(panel.grid=element_blank()) + #facet_wrap(~SPECIESID) + 
    labs(x=expression(delta), y="Residual error", title=response_vars[k])
}


R2.df
walk(p_y, print)
walk(p_d, print)
walk(resid_y, print)
walk(resid_d, print)



p1.sp <- p2.sp <- p3.sp <- vector("list", n_distinct(wkr.df_i$SPECIESID))
a12.sp <- a13.sp <- a23.sp <- vector("list", n_distinct(wkr.df_i$SPECIESID))
b12.sp <- b13.sp <- b23.sp <- vector("list", n_distinct(wkr.df_i$SPECIESID))
for(i in 1:n_distinct(wkr.df_i$SPECIESID)) {
  p1.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[1,", "b[1,"), i, "]"))
  p2.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[2,", "b[2,"), i, "]"))
  p3.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[3,", "b[3,"), i, "]"))
  
  a12.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[1,", "a[2,"), i, "]"))
  a13.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[1,", "a[3,"), i, "]"))
  a23.sp[[i]] <- stan_scat(out[[3]], paste0(c("a[2,", "a[3,"), i, "]"))
  
  b12.sp[[i]] <- stan_scat(out[[3]], paste0(c("b[1,", "b[2,"), i, "]"))
  b13.sp[[i]] <- stan_scat(out[[3]], paste0(c("b[1,", "b[3,"), i, "]"))
  b23.sp[[i]] <- stan_scat(out[[3]], paste0(c("b[2,", "b[3,"), i, "]"))
}

ggpubr::ggarrange(plotlist=p1.sp)
ggpubr::ggarrange(plotlist=p2.sp)
ggpubr::ggarrange(plotlist=p3.sp)

ggpubr::ggarrange(plotlist=a12.sp)
ggpubr::ggarrange(plotlist=a13.sp)
ggpubr::ggarrange(plotlist=a23.sp)

ggpubr::ggarrange(plotlist=b12.sp)
ggpubr::ggarrange(plotlist=b13.sp)
ggpubr::ggarrange(plotlist=b23.sp)



wkr.df_i %>% mutate(Cnpy=case_when(CnpyOpen==1 ~ "Open",
                                   CnpyMixed==1 ~ "Mixed",
                                   CnpyOpen==0 & CnpyMixed==0 ~ "Closed")) %>%
  ggplot(aes(v, colour=Cnpy)) + geom_density() + facet_wrap(~SPECIESID)





spp.tb <- bind_rows(
  as_tibble(summary_hdi(out[["clny_est"]], "delta_exp")) %>% 
    bind_cols(out.clny %>% group_by(SPECIESID) %>% 
                summarise(obs=mean(RelLegHind_sd, na.rm=T))) %>%
    mutate(par="delta_exp"),
  as_tibble(summary_hdi(out[["clny_est"]], "mu")) %>% 
    bind_cols(out.clny %>% group_by(SPECIESID) %>% 
                summarise(obs=mean(RelLegHind, na.rm=T))) %>%
    mutate(par="mu")
  ) 

spp.tb %>% filter(par=="mu") %>%
  ggplot(aes(x=obs, y=median, ymin=L05, ymax=L95)) + 
  geom_abline() +
  geom_pointrange()
spp.tb %>% filter(par=="delta_exp") %>%
  ggplot(aes(x=obs, y=median, ymin=L05, ymax=L95)) + 
  geom_abline() +
  geom_pointrange()

ggplot(out.clny, aes(RelLegHind, y_bar_median, ymin=y_bar_L05, ymax=y_bar_L95)) + 
  geom_abline() +
  geom_linerange(size=0.2) + geom_point(shape=1, size=0.5)

ggplot(out.clny, aes(log(RelLegHind_sd), d_log_median, ymin=d_log_L05, ymax=d_log_L95)) + 
  geom_abline() +
  geom_linerange(size=0.2) + geom_point(shape=1, size=0.5)



