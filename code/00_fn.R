# Auxilliary functions for keeping things cleaner




#' LAS measurements: Aggregate measurements for a single worker
#' 
#' @param f_LAS file name of LAS .xslx report, including directory
#' @param na.thresh minimum trait value; values smaller are marked NA
#' @return Dataframe with a row for each measurement and columns for Trait
#'   (numeric), Value (mm), and Worker

aggregate_measurements <- function(f_LAS, na.thresh=0.05) {
  map2_dfr(f_LAS, grep("Profile", excel_sheets(f_LAS), value=T),
           ~read_xlsx(.x, .y, skip=1, progress=F,
                      col_types=rep(c("numeric", "skip"), each=2)) %>% 
             filter(row_number() %in% c(1,n())) %>% 
             mutate(unit=str_sub(names(.)[1], -3, -2)) %>% 
             setNames(c("x", "y", "unit")) %>%
             summarise(Value=sqrt(diff(x)^2 + diff(y)^2) %>%
                         units::as_units(., first(unit)) %>%
                         units::set_units(., "mm") %>%
                         as.numeric %>%
                         replace(., . < na.thresh, NA)) %>%
             mutate(Worker=str_sub(.x, -6, -6)), 
           .id="Trait") 
}







#' Add values from raster covariates
#'
#' @param ant_i dataframe with ant species, covariates, etc for each tube number
#' @param cov_f named list of raster file paths
#' @return Updated ant_i with the values for each point from each raster in
#'   cov_f, where the column names are the names of the items in cov_f

add_covariates <- function(ant_i, cov_f) {
  
  for(i in seq_along(cov_f)) {
    ant_i[[names(cov_f)[i]]] <- raster::extract(raster::raster(cov_f[[i]]),
                                                ant_i, fun=mean)
  }
  
  return(ant_i)
}






#' Load and aggregate all measurements (LAS, color) within specified directories
#'
#' @param ant_i dataframe with ant species, covariates, etc for each tube number
#' @param msr_dir directory containing LAS xlsx reports; recursive search
#' @param col_dir directory containing color histograms; recursive search
#' @param bg_dir directory containing gray background values
#' @param na.thresh minimum trait value; values smaller are marked NA
#' @param lat_names vector of trait names for lateral view
#' @param fro_names vector of trait names for frontal view
#' @param dor_names vector of trait names for dorsal view
#' @return List of two data frames: wkr.df = worker-level traits, clny.df =
#'   colony-level mean and sd values for each trait. Each is joined by all
#'   columns from ant_i, but only tubes with trait measurements are included

load_traits <- function(ant_i, msr_dir, col_dir, bg_dir, na.thresh=0.05,
                        lat_names, fro_names, dor_names) {
  # setup
  msr.tubes <- dir(msr_dir, "^[1-9]", include.dirs=T, recursive=T)
  col.tubes <- dir(col_dir, "^[1-9]", include.dirs=T, recursive=T)
  fro_all <- lat_all <- dor_all <- vector("list", length(msr.tubes))
  
  pb <- txtProgressBar(min=1, max=length(msr.tubes)+length(col.tubes))
  for(tube in seq_along(msr.tubes)) {
    setTxtProgressBar(pb, tube)
    
    # files to read
    tubeSubDir <- msr.tubes[tube]
    tubeNo <- str_split(tubeSubDir, "/")[[1]] %>% last()
    f_fro <- dir(paste0(msr_dir, tubeSubDir), "fro[0-9]\\.xlsx", full.names=T)
    f_lat <- dir(paste0(msr_dir, tubeSubDir), "lat[0-9]\\.xlsx", full.names=T)
    f_dor <- dir(paste0(msr_dir, tubeSubDir), "dor[0-9]\\.xlsx", full.names=T)
    
  
    if(length(f_fro) > 0) {
      fro_all[[tube]] <- map_dfr(f_fro, aggregate_measurements) %>%
        mutate(TubeNo=tubeNo, unit="mm", Trait=fro_names[as.numeric(Trait)])
    }
    
    if(length(f_lat) > 0) {
      lat_all[[tube]] <- map_dfr(f_lat, aggregate_measurements) %>%
        mutate(TubeNo=tubeNo, unit="mm", Trait=lat_names[as.numeric(Trait)])
    }
    
    if(length(f_dor) > 0) {
      dor_all[[tube]] <- map_dfr(f_dor, aggregate_measurements) %>%
        mutate(TubeNo=tubeNo, unit="mm", Trait=dor_names[as.numeric(Trait)])
    }
  }
  col.ls <- grey.ls <- vector("list", length(col.tubes))
  for(tube in seq_along(col.tubes)) {
    setTxtProgressBar(pb, tube+length(msr.tubes))
    tubeNo <- col.tubes[tube]
    wkr.f <- dir(paste0(col_dir, tubeNo), "meso[1-9].txt", full.names=T)
    grey.ls[[tube]] <- imap(wkr.f, ~read.table(.x, fileEncoding="UTF-16LE", 
                                               header=T, nrows=1) %>%
                              select(Gray.Value..Median., Gray.Value..Mean.) %>%
                              mutate(TubeNo=str_split_fixed(col.tubes[tube], "/", 3)[3], 
                                     Worker=as.character(.y)) %>%
                              rename(grey_md=Gray.Value..Median.,
                                     grey_mn=Gray.Value..Mean.)) %>%
      do.call('rbind', .)
    nWorkers <- length(wkr.f)
    meso.hist <- vector("list", length=nWorkers)
    for(i in 1:nWorkers) {
      meso_i.f <- dir(str_replace(wkr.f[i], ".txt", " Data/"), full.names=T)
      if(length(meso_i.f) > 0) {
        meso.hist[[i]] <- map(meso_i.f, ~t(as.matrix(read.csv(., header=F)))) %>%
          Reduce('+', .) %>%
          as.data.frame(.) %>% setNames(., c("R", "G", "B")) %>%
          mutate(value=1:256,
                 TubeNo=str_split_fixed(col.tubes[tube], "/", 3)[3],
                 Worker=as.character(i))
      }
    }
    col.ls[[tube]] <- rbind(do.call('rbind', meso.hist))
  }
  close(pb)
  
  # gray backgrounds
  grey_bg <- dir(bg_dir, full.names=T) %>%
    map_dfr(., ~suppressMessages(read_csv(.x))) %>%
    mutate(TubeNo=str_split_fixed(Label, "-", 2)[,1], 
           Worker=str_sub(str_split_fixed(Label, "-", 2)[,2], 1, 1)) %>%
    select(TubeNo, Worker, Mode) %>%
    rename(bg_v=Mode)
  
  # color dataframes
  col.wkr <- do.call('rbind', col.ls) %>%
    gather(channel, count, 1:3) %>%
    group_by(TubeNo, Worker, channel) %>%
    summarise(mn=sum(value*count)/sum(count),
              med=median(rep(value, count)), 
              sd=sd(rep(value, count))) %>% 
    ungroup %>% pivot_wider(names_from=channel, values_from=c(mn, med, sd)) %>%
    mutate(rgbCol=paste0("rgb(", med_R/256, ",", med_G/256, ",", med_B/256, ")"),
           h=rgb2hsv(med_R, med_G, med_B, maxColorValue=255)[1,],
           s=rgb2hsv(med_R, med_G, med_B, maxColorValue=255)[2,],
           v=rgb2hsv(med_R, med_G, med_B, maxColorValue=255)[3,]) %>%
    full_join(., do.call('rbind', grey.ls), by=c("TubeNo", "Worker")) %>%
    full_join(., grey_bg, by=c("TubeNo", "Worker")) %>%
    mutate(median_bg=median(bg_v, na.rm=T),
           v=v/(bg_v/median_bg),
           grey_md=grey_md/(bg_v/median_bg),
           grey_mn=grey_mn/(bg_v/median_bg))
  col.clny <- col.wkr %>% group_by(TubeNo) %>%
    summarise(R=mean(med_R), G=mean(med_G), B=mean(med_B), 
              v_var=var(v), v_CV=sd(v)/mean(v), v=mean(v), 
              grey_var=var(grey_md), grey_CV=sd(grey_md)/mean(grey_mn),
              grey_mn=mean(grey_mn), grey_md=mean(grey_md)) %>%
    mutate(rgbCol=paste0("rgb(", R/256, ",", G/256, ",", B/256, ")")) 
  
  # aggregated dataframes
  wkr.df <- rbind(do.call("rbind", fro_all), 
                     do.call("rbind", lat_all), 
                     do.call("rbind", dor_all)) %>%
    filter(!is.na(Value)) %>%
    full_join(., col.wkr, by=c("TubeNo", "Worker")) %>%
    left_join(., ant_i, by="TubeNo") %>%
    group_by(SPECIESID) %>% mutate(n_clny=n_distinct(TubeNo)) %>% ungroup
  
  clny.df <- wkr.df %>% group_by(TubeNo, Trait) %>% 
    summarise(mnValue=mean(Value, na.rm=T),
              sdValue=sd(Value, na.rm=T)) %>%
    full_join(., col.clny, by=c("TubeNo")) %>%
    left_join(., ant_i, by="TubeNo") %>%
    group_by(SPECIESID) %>% mutate(n_clny=n_distinct(TubeNo)) %>% ungroup
  
  # standardized by species
  wkr.std <- wkr.df %>% 
    group_by(SPECIESID, Trait) %>%
    mutate(Value=(Value-mean(Value, na.rm=T))/sd(Value, na.rm=T), 
           v=(v-mean(v, na.rm=T))/sd(v, na.rm=T),
           grey_mn=(grey_mn-mean(grey_mn, na.rm=T))/sd(grey_mn, na.rm=T),
           grey_md=(grey_md-mean(grey_md, na.rm=T))/sd(grey_md, na.rm=T))
  clny.std <- clny.df %>% 
    group_by(SPECIESID, Trait) %>%
    mutate(mnValue=(mnValue-mean(mnValue, na.rm=T))/sd(mnValue, na.rm=T), 
           v=(v-mean(v, na.rm=T))/sd(v, na.rm=T),
           grey_mn=(grey_mn-mean(grey_mn, na.rm=T))/sd(grey_mn, na.rm=T),
           grey_md=(grey_md-mean(grey_md, na.rm=T))/sd(grey_md, na.rm=T))
  
  # wide format
  wkr.wide <- wkr.df %>% 
    pivot_wider(names_from=Trait, values_from=Value) %>%
    mutate(MidLen=MidFemur+MidTibia, HindLen=HindFemur+HindTibia, 
           MesoSA=MesosomaLength*MesosomaWidth)
  clny.wide <- clny.df %>% 
    pivot_wider(names_from=Trait, values_from=c(mnValue, sdValue)) %>%
    mutate(MidLen=mnValue_MidFemur+mnValue_MidTibia, 
           HindLen=mnValue_HindFemur+mnValue_HindTibia, 
           MesoSA=mnValue_MesosomaLength*mnValue_MesosomaWidth)
  
  return(list(wkr.df=wkr.df, clny.df=clny.df, 
              wkr.std=wkr.std, clny.std=clny.std, 
              wkr.wide=wkr.wide, clny.wide=clny.wide))
}






#' Prior predictive check
#'
#' Generate simulated trait distributions for colonies and workers based on
#' specified priors. Useful for confirming that priors are constrained to
#' reasonable values, particularly for the heteroskedastic term, which can
#' generate absurdly large values for standard deviations if priors are overly
#' broad.
#'
#' @param d.ls list to use as basis for generated traits. Must include elements
#'   \code{.$S}, \code{.$G}, \code{.$P_mn}, \code{.$P_sd}, \code{.$N_clny},
#'   \code{.$spp_id}, \code{.$tax_i}, \code{.$x_mn}, \code{.$x_sd}, and
#'   \code{.$clny_id}
#' @param hyperpriors list of hyperprior values. Must include
#'   \code{.$sigma_clny_1_global}, \code{.$sigma_clny_2_global}, \code{.$alpha},
#'   \code{.$beta}, \code{.$sigma_A}, \code{.$sigma_B}, \code{.$sigma_a},
#'   \code{.$sigma_b}, \code{.$L_A}, and \code{.$L_B}
#' @return list with \code{.$priors} and \code{.$sims}, with all prior
#'   distributions and predictions for colony-level \code{.$mu}, \code{.$delta},
#'   \code{.$y_bar}, \code{.$d}, and worker-level \code{.$y}, where the values
#'   correspond with those in \code{d.ls}

prior_predictive_check <- function(d.ls, hyperpriors) {
  library(tidyverse)
  
  # Generate priors for distributions dependent on hyperpriors
  priors <- hyperpriors %>%
    list_modify(sigma_clny_1=rexp(d.ls$S, .$sigma_clny_1_global),
                sigma_clny_2=rexp(d.ls$S, .$sigma_clny_2_global),
                z_A=matrix(rnorm(d.ls$P_sd*d.ls$G, 0, 1), d.ls$P_sd, d.ls$G),
                z_B=matrix(rnorm(d.ls$P_mn*d.ls$G, 0, 1), d.ls$P_mn, d.ls$G),
                z_a=matrix(rnorm(d.ls$P_sd*d.ls$S, 0, 1), d.ls$P_sd, d.ls$S),
                z_b=matrix(rnorm(d.ls$P_mn*d.ls$S, 0, 1), d.ls$P_mn, d.ls$S),
                L_A=rethinking::rlkjcorr(n=1, K=d.ls$P_sd, eta=1),
                L_B=rethinking::rlkjcorr(n=1, K=d.ls$P_mn, eta=1)) %>%
    list_modify(#resid_clny_1=rnorm(d.ls$N_clny, 0, .$sigma_clny_1[d.ls$spp_id]),
                #resid_clny_2=rnorm(d.ls$N_clny, 0, .$sigma_clny_2[d.ls$spp_id]),
                A=(diag(.$sigma_A) %*% .$L_A) %*% .$z_A,
                B=(diag(.$sigma_B) %*% .$L_B) %*% .$z_B) %>%
    list_modify(a=t(sapply(1:d.ls$P_sd, function(x) .$A[x,d.ls$tax_i[,2]] + 
                             .$z_a[x,] * .$sigma_a[x])),
                b=t(sapply(1:d.ls$P_mn, function(x) .$B[x,d.ls$tax_i[,2]] + 
                             .$z_b[x,] * .$sigma_b[x]))) 
  
  sims <- list(mu=map_dbl(1:d.ls$N_clny, 
                          ~d.ls$x_mn[.x,] %*% 
                            (priors$beta + priors$b[,d.ls$spp_id[.x]])),
               delta=map_dbl(1:d.ls$N_clny, 
                             ~d.ls$x_sd[.x,] %*% 
                               (priors$alpha + priors$a[,d.ls$spp_id[.x]]))) %>%
    # list_merge(y_bar=.$mu + priors$resid_clny_1,
    #            d=exp(.$delta + priors$resid_clny_2)) %>%
    list_merge(y_bar=rnorm(d.ls$N_clny, .$mu, priors$sigma_clny_1[d.ls$spp_id]),
               d=truncnorm::rtruncnorm(d.ls$N_clny, 0, Inf,
                                       exp(.$delta), priors$sigma_clny_2[d.ls$spp_id])) %>%
    list_merge(y=rnorm(d.ls$N_wkr, .$y_bar[d.ls$clny_id], .$d[d.ls$clny_id]))
  
  return(list(priors=priors, sims=sims))
}








#' Summarise Stan output with mean, HDIs, n_eff, and Rhat
#'
#' @param out.stan stan object with posterior samples
#' @param param parameter of interest
#' @param creds vector of credible intervals
#' @param varNames optional vector of variable names (e.g., for beta covariates)
#' @return Tibble with mean, median, n_eff, Rhat, X_var (if provided), and HDI limits
#' 
summary_hdi <- function(out.stan, param, creds=c(0.5, 0.9, 0.95), varNames=NULL) {
  library(rstan); library(tidyverse)
  
  out.df <- summary(out.stan, param)$summary %>%
    as.data.frame %>% rename(median=`50%`) %>%
    select(mean, median, n_eff, Rhat) 
  if(!is.null(varNames)) {
    out.df <- out.df %>% mutate(X_var=varNames)
  }
  
  # par.post <- rstan::extract(out.stan, param)
  par.post <- rstan::As.mcmc.list(out.stan, param)
  out.df <- bind_cols(
    out.df, 
    map(creds, ~t(HDInterval::hdi(par.post[[1]], .x)) %>%
          magrittr::set_colnames(
            paste0("L", str_sub( c(0.5-.x/2, 0.5+.x/2), 3, -1) )) %>%
          as_tibble) %>%
      do.call('cbind', .)
  ) 
  
  return(as_tibble(out.df))
}







#' Summarise Stan output to the colony level
#'
#' @param out.stan stan object with posterior samples
#' @param params parameter of interest
#' @param clny_df dataframe of worker-level traits
#' @param trts list of trait dataframes produced by \code{load_traits()}
#' @return Updated clny_df with HDI summaries for params
#' 
summarise_clny <- function(out.stan, params, clny_df, trts) {
  library(tidyverse)
  clny_df %>% ungroup %>%
    mutate(El_m=trts$clny.wide$mnt25[match(TubeNo, trts$clny.wide$TubeNo)]) %>%
    bind_cols(map(params,
                  ~summary_hdi(out.stan, .x) %>%
                    set_names(paste0(.x, "_", names(.)))))
}






#' Summarise Stan output to the worker level
#'
#' @param out.stan stan object with posterior samples
#' @param params parameter of interest
#' @param wkr_df dataframe of worker-level traits
#' @param clny_out Updated clny_df from \code{summarise_clny()}
#' @return Updated wkr_df with HDI summaries for params
#' 
summarise_wkr <- function(out.stan, params, wkr_df, clny_out) {
  library(tidyverse)
  wkr_df %>% ungroup %>%
    full_join(clny_out %>% select(clny_id, all_of(paste0(params, "_mean"))), 
              by="clny_id") %>%
    bind_cols(map(c("y_pred", "y_pred_"), 
                  ~summary_hdi(out.stan, .x) %>%
                    set_names(paste0(.x, "_", names(.)))))
}









#' Save output from summarised stanfit object
#'
#' @param out_dir directory to save output
#' @param Y_var response variable of interest
#' @param mod_i name of model from `mods`
#' @param mods.ls lookup table with model names and file extension abbreviations
#' @param out.ls named list of summarised output (summarised with
#'   `summarise_wkr()`, `summarise_clny()`, or `summary_hdi()`)
#' @return Updated wkr_df with HDI summaries for params
#' 
save_summaries <- function(out_dir, Y_var, mod_i, mods.lu, out.ls) {
  library(tidyverse)
  imap(out.ls, ~write_csv(.x, paste0(out_dir, 
                                     filter(mods.lu, name==mod_i)$abbr, "_", 
                                     .y, "_", 
                                     Y_var, ".csv")))
}








#' Plot points and linear smoother in GGally
ggally_pts_smooth <- function(data, mapping, method="lm", se=F,
                              size_pt=0.4, size_smooth=0.3,
                              alpha_pt=0.4, alpha_smooth=1, ...) {
  p <- ggplot(data=data, mapping=mapping) + 
    geom_point(shape=1, size=size_pt, alpha=alpha_pt) + 
    geom_smooth(method=method, formula=y~x, se=se, 
                size=size_smooth, alpha=alpha_smooth)
  p
}






#' Plot of species-standardized color against a chosen variable
#'
#' @param tr.df dataframe with color values; currently 'all.sum.wide' in munge_color.R
#' @param x_var character for which column to plot on the x-axis 
#' @param xlab optional x-axis label
#' @return Plot of x_var on x-axis and brightness standardized within each
#'   species as (v - mean(v))/sd(v), with solid points for colonies and hollow
#'   points for individual workers, colour as the actual mean colour, simple
#'   linear regression best-fit lines (dotted = species, solid = total)

plot_colors_by_v <- function(tr.df, x_var, xlab=NULL) {
  library(tidyverse)
  if(is.null(xlab)) xlab <- x_var
  sp.df <- tr.df %>% group_by(SPECIESID) %>% mutate(v_std=(v-mean(v))/sd(v))
  tube.df <- sp.df %>%  group_by(TubeNo) %>% 
    summarise_if(is.numeric, mean)
  reg_sp <- sp.df %>% summarise(n=n_distinct(TubeNo)) %>% filter(n > 2)
  plot(sp.df[[x_var]], sp.df$v_std, 
       pch=1, cex=1, col=rgb(sp.df$med_R, sp.df$med_G, sp.df$med_B, maxColorValue=255),
       xlab=xlab, 
       ylab=expression(Standardized~brightness:~frac(v-bar(v[s]),sigma[v[s]])))
  points(tube.df[[x_var]], tube.df$v_std, pch=19, cex=1.5, 
         col=rgb(tube.df$med_R, tube.df$med_G, tube.df$med_B, maxColorValue=255))
  abline(lm(sp.df$v_std ~ sp.df[[x_var]]), lwd=2)
  for(i in 1:n_distinct(reg_sp$SPECIESID)) {
    sp_i.df <- filter(sp.df, SPECIESID==reg_sp$SPECIESID[i])
    clip(range(sp_i.df[[x_var]], na.rm=T)[1], 
         range(sp_i.df[[x_var]], na.rm=T)[2], 
         range(sp_i.df$v_std)[1], range(sp_i.df$v_std)[2])
    abline(lm(sp_i.df$v_std ~ sp_i.df[[x_var]]), lty=2)
  }
}





