# ETP by season 
# By: Cesar A
# CIAT
# 2023
# =-----------------------------------------------=
# R options
g <- gc(reset = T); rm(list = ls())    # Empty garbage collector
# .rs.restartR()                       # Restart R session
options(warn = -1, scipen = 999)       # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

OSys <<- Sys.info()[1]
root <<- '//CATALOGUE/BanRep_papa' 
# =-----------------------------------------------=
# inputs and outputs
rgn = 'Antioquia' #'Antioquia', 'Boyaca', 'Cauca', 'Cundinamarca', 'N_Santander', 'Narino', 'Santander', 'Tolima'
seasons <- list(s1 = c(3:8), s2 = c(9:12,1:2)) # s2 = c(10:12,1:2)
shp_fl <- paste0(root, '/1.Data/Shapefile/',rgn,'/Municipios_Paperos/','Municipios_Paperos_',rgn,'.shp')
outfile <- paste0(root, '/7.Results/agro_indices/',rgn, '/'); if(!dir.exists(outfile)){dir.create(outfile, F, T)}
# =-----------------------------------------------=
# Function to compute basic Agro-climatic indices
calc_AgrClm <- function(season = season, shp_fl = shp_fl){
  
  ## ROI: regions of interest
  shp <- terra::vect(shp_fl)
  
  ## Daily files
  # Precipitation
  chr_pth <- paste0(root, '/1.Data/Chirps')
  chr_fls <- gtools::mixedsort(list.files(chr_pth, pattern = '*.tif$', full.names = T))
  chr_dts <- strsplit(x = chr_fls, split = 'chirps-', fixed = T) %>% purrr::map(2) %>% unlist()
  chr_dts <- strsplit(x = chr_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  chr_dts <- as.Date(gsub('_', '-', chr_dts, fixed = T))
  
  # Tmax
  era5Dir <- paste0(root, '/1.Data/agERA5')
  tmx_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
  tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '_maximum_.*.tif$', full.names = T))
  tmx_dts <- strsplit(x = tmx_fls, split = '2m_temperature-24_hour_maximum_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmx_dts <- strsplit(x = tmx_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  tmx_dts <- as.Date(gsub('_', '-', tmx_dts, fixed = T))
  
  # Tmin
  tmn_pth <- paste0(era5Dir,'/2m_temperature-24_hour_minimum')
  tmn_fls <- gtools::mixedsort(list.files(tmn_pth, pattern = '_minimum_.*.tif$', full.names = T))
  tmn_dts <- strsplit(x = tmn_fls, split = '2m_temperature-24_hour_minimum_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmn_dts <- strsplit(x = tmn_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  tmn_dts <- as.Date(gsub('_', '-', tmn_dts, fixed = T))
  
  # Solar radiation
  srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
  srd_fls <- gtools::mixedsort(list.files(srd_pth, pattern = '*.tif$', full.names = T))
  srd_dts <- strsplit(x = srd_fls, split = 'solar_radiation_flux_', fixed = T) %>% purrr::map(2) %>% unlist()
  srd_dts <- strsplit(x = srd_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  srd_dts <- as.Date(gsub('_', '-', srd_dts, fixed = T))
  
  # Filtering days within the season
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  
  tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
  tmn_fls <- tmn_fls[lubridate::year(tmn_dts) %in% yrs]
  srd_fls <- srd_fls[lubridate::year(srd_dts) %in% yrs]
  
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  tmn_dts <- tmn_dts[lubridate::year(tmn_dts) %in% yrs]
  srd_dts <- srd_dts[lubridate::year(srd_dts) %in% yrs]
  
  if(length(season) < 12){
    cnd <- lubridate::month(tmx_dts) %in% season # Days within the season
    yrs_dts <<- split(tmx_dts[cnd],cumsum(c(1,diff(tmx_dts[cnd])!=1)))
  } else {
    yrs <- lubridate::year(tmx_dts)
    grp <- with(rle(yrs), rep(seq_along(values), lengths))
    yrs_dts <<- split(tmx_dts, grp)
  }
  
  # Precipitation indices 
  cat('..... Computing: ATP. Priestley-Taylor Model Evapotranspiration. \n')
  ETP <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      ref <<- terra::rast(paste0(root,"/1.Data/chirps-v2.0.2020.01.01.tif"))
      #
      cat(paste0('..... Computing: ATP ', lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-'), '\n'))
      cat('..... Read vars: \n')
      cat('..... Tmax \n')
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- terra::resample(tmx,ref) %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      #
      cat('..... Tmin \n')
      tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
      tmn <- terra::resample(tmn,ref) %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmn <- tmn - 273.15
      #
      cat('..... Srad \n')
      sr <- terra::rast(srd_fls[srd_dts %in% yrs_dts[[i]]])
      sr <- terra::resample(sr,ref) %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      sr <- sr/1000000 # solar radiation from j to kj
      sr <- (1/2.45)*(sr) # Solar radiation from kj to mm day-1
      # 
      cat('..... calc Tmean \n')
      tmean <<- (tmx+tmn)/2
      #
      cat('..... read elevation \n')
      elev <- terra::rast(paste0(root,'/1.Data/Elevation/COL_elev.tif'))
      elev <- terra::resample(elev,ref) %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      lst <- list()
      for(i in 1:length(yrs_dts[[i]])){
        lst[[i]] <- elev
      }
      elv <- terra::rast(lst)
      #
      cat('..... Computing: calc constants \n')
      # Priestley and Taylor evaporative coefficient a = 1.26
      a <- tmean
      terra::values(a) <- 1.26  
      # Latent heat of vaporization at 20°C  b = 2.45
      b <- tmean
      terra::values(b) <- 2.45
      # Relationship between saturation vapor pressure and air temperature
      d <- (4098*(0.6108*(exp((17.27*tmean)/(tmean+237.3)))))/((tmean+237.3)*(tmean+237.3)) 
      # Psychrometric constant
      e <- (-0.000007*elv)+0.0666
      #
      ## ETP. evapotranspiracion potencial
      calc_etp <- function(a, b, d, e, sr){
        ETP = a*(d/(d+e))*(sr/b)  #Priestley-Taylor Model
        return(ETP)
      }
      #
      cat('..... Computing: calc ETP \n')
      ETP <- terra::lapp(x = terra::sds(a, b, d, e, sr), 
                         fun = calc_etp) %>% sum()
      # names(ETP) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(ETP)
    }) %>% terra::rast()
  ETP <- ETP %>% terra::mask(shp)
  #
  cat('..... End.\n')
  return(list(ETP = ETP))
}
# =-----------------------------------------------=
# Loop through seasons
1:length(seasons) %>%
  purrr::map(.f = function(s){
    cat(paste0('Processing season ',names(seasons)[s],':\n'))
    # Indices calculation
    indices <- calc_AgrClm(seasons[[s]], shp_fl)
    # Load 5 km raster template
    tmp <- terra::rast(paste0(root, '/1.Data/chirps-v2.0.2020.01.01.tif'))
    shp <- terra::vect(shp_fl)
    tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
    tmp[!is.na(tmp)] <- 1
    # Indices resampling
    indices <- indices %>% purrr::map(.f = function(r){r <- r %>% terra::resample(x = ., y = tmp) %>% terra::mask(shp); return(r)})
    # Saving results
    out <- paste0(outfile, names(seasons)[s]); if(!dir.exists(out)){dir.create(out,F,T)}
    1:length(names(indices)) %>%
      purrr::map(.f = function(j){
        terra::writeRaster(x = indices[[j]], filename = paste0(out,'/',names(indices)[j],'.tif'), overwrite = T)
      })
    return(cat('Process finished successfully!\n'))
  })
# =-----------------------------------------------=
