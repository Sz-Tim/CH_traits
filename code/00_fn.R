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
#' @param na.thresh minimum trait value; values smaller are marked NA
#' @param lat_names vector of trait names for lateral view
#' @param fro_names vector of trait names for frontal view
#' @param dor_names vector of trait names for dorsal view
#' @return List of two data frames: wkr.df = worker-level traits, clny.df =
#'   colony-level mean and sd values for each trait. Each is joined by all
#'   columns from ant_i, but only tubes with trait measurements are included

load_traits <- function(ant_i, msr_dir, col_dir, na.thresh=0.05,
                        lat_names, fro_names, dor_names) {
  # setup
  msr.tubes <- dir(msr_dir, "^999", include.dirs=T, recursive=T)
  col.tubes <- dir(col_dir, "^999", include.dirs=T, recursive=T)
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
  col.ls <- vector("list", length(col.tubes))
  for(tube in seq_along(col.tubes)) {
    setTxtProgressBar(pb, tube+length(msr.tubes))
    tubeNo <- col.tubes[tube]
    nWorkers <- length(dir(paste0(col_dir, tubeNo), "meso[1-9].txt"))
    meso.hist <- vector("list", length=nWorkers)
    for(i in 1:nWorkers) {
      meso_i.f <- dir(paste0(col_dir, tubeNo, "/meso", i, " Data/"), full.names=T)
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
           v=rgb2hsv(med_R, med_G, med_B, maxColorValue=255)[3,])
  col.clny <- col.wkr %>% group_by(TubeNo) %>%
    summarise(R=mean(med_R), G=mean(med_G), B=mean(med_B), 
              v_var=var(v), v_CV=sd(v)/mean(v)) %>%
    mutate(rgbCol=paste0("rgb(", R/256, ",", G/256, ",", B/256, ")"),
           h=rgb2hsv(R, G, B, maxColorValue=255)[1,],
           s=rgb2hsv(R, G, B, maxColorValue=255)[2,],
           v=rgb2hsv(R, G, B, maxColorValue=255)[3,]) 
  
  # aggregated dataframes
  wkr.df <- rbind(do.call("rbind", fro_all), 
                     do.call("rbind", lat_all), 
                     do.call("rbind", dor_all)) %>%
    full_join(., col.wkr, by=c("TubeNo", "Worker")) %>%
    left_join(., ant_i, by="TubeNo")
  
  clny.df <- wkr.df %>% group_by(TubeNo, Trait) %>% 
    summarise(mnValue=mean(Value),
              sdValue=sd(Value)) %>%
    full_join(., col.clny, by=c("TubeNo")) %>%
    left_join(., ant_i, by="TubeNo")
  
  # standardized by species
  wkr.std <- wkr.df %>% 
    group_by(SPECIESID, Trait) %>%
    mutate(Value=(Value-mean(Value, na.rm=T))/sd(Value, na.rm=T), 
           v=(v-mean(v, na.rm=T))/sd(v, na.rm=T))
  clny.std <- clny.df %>% 
    group_by(SPECIESID, Trait) %>%
    mutate(mnValue=(mnValue-mean(mnValue, na.rm=T))/sd(mnValue, na.rm=T), 
           v=(v-mean(v, na.rm=T))/sd(v, na.rm=T))
  
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


