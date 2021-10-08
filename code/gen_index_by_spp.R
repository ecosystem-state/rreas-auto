library(mgcv)
library(rerddap)
library(dplyr)
library(lubridate)
library(sf)
library(tidyr)
library(purrr)

# settings
source("code/set_control_params.R")

# load species tot cpue
tot_cpue <- readRDS("indices/tot_cpue_species.rds")
#tot_cpue <- dplyr::arrange(tot_cpue,-tot_cpue)
#tot_cpue <- tot_cpue[1:top_species,]

# load prediction grid
pred_grid <- readRDS("indices/pred_grid.rds")

# model function - edit to change form of model
gam_pres_fit <- function(df) {
  gam(pres ~ jday + s(latitude, longitude) + year +I(jday^2),
      data = df,
      family = "binomial"
  )
}
gam_pos_fit <- function(df) {
  gam(count ~ jday + s(latitude, longitude) + year +I(jday^2),
      data = df[which(df$count > 0),],
      family = "poisson"
  )
}

# grab data for all species
unique_files <- unique(tot_cpue$erddap)
for (i in 1:length(unique_files)) {
  out <- info(as.character(unique_files[i]))
  # station_dat <- tabledap(out, fields = c(
  #   "station", "line", "latitude",
  #   "longitude", "time", "scientific_name", "larvae_10m2"
  # ))
  dat <- tabledap(
    out,
    fields = c(
      "common_name", "latitude",
      "longitude", "maturity",
      "sci_name", "species_group",
      "station_bottom_depth","time",
      "catch","cruise", "haul_no",
      "strata","station"
    )
  )
  # filter out species in question
  #species <- dplyr::filter(tot_cpue, erddap == unique_files[i])
  #dat <- dplyr::filter(station_dat,
  #                     sci_name %in% tot_cpue$sci_name)

  # convert date
  dat <- as.data.frame(dat)
  dat$date <- lubridate::as_date(dat$time)
  dat$year <- lubridate::year(dat$date)
  dat$month <- lubridate::month(dat$date)
  dat$jday <- lubridate::yday(dat$date)

  # convert to UTM - kms
  # make the UTM cols spatial (X/Easting/lon, Y/Northing/lat)
  dat$latitude <- as.numeric(dat$latitude)
  dat$longitude <- as.numeric(dat$longitude)
  dat <-
    st_as_sf(dat,
             coords = c("longitude", "latitude"),
             crs = 4326)
  dat <- st_transform(x = dat, crs = 32610)
  dat$longitude = st_coordinates(dat)[, 1]
  dat$latitude = st_coordinates(dat)[, 2]
  dat <- as.data.frame(dat)
  dat$longitude <- dat$longitude / 1000 # to kms
  dat$latitude <- dat$latitude / 1000 # to kms

  dat$uniqueID<-do.call(paste, c(dat[c("cruise", "haul_no")], sep = "_"))

  # sum total counts for spp we don't care about maturity for
  dat$uniqueID<-do.call(paste, c(dat[c("cruise", "haul_no")], sep = "_"))
  spp_all<-c('Sardinops sagax','Doryteuthis opalescens','Euphausiacea')
  dat_allStages<-dat[dat$sci_name %in% spp_all,]
  dat_all<-dat_allStages %>%
    group_by(year,uniqueID,latitude,longitude,jday,sci_name,common_name) %>%
    summarise(count=sum(catch))%>%
    data.frame

  ## For those species that we want only juvenile stages (or only juv stages identified)
  #keep those data where maturity == "Y"
  spp_juvenile<-c('Sardinops sagax','Engraulis mordax','Citharichthys stigmaeus',
                  'Merluccius productus','Sebastes goodei','Sebastes entomelas',
                  'Sebastes jordani','Sebastes mystinus','Sebastes paucispinis',
                  'Sebastes semicinctus','Citharichthys sordidus')
  spp_extra <- c("Sebastes melanostomus", "Sebastes wilsoni","Sebastes pinniger",
                 "Sebastes auriculatus","Sebastes flavidus","Sebastes spp. caurinus complex",
                 "Sebastes hopkinsi","Sebastes saxicola")
  yoy <- dplyr::filter(dat, sci_name%in% c(spp_extra, spp_juvenile),
                       maturity=="Y") %>%
        dplyr::rename(count = catch)

  yoy_juv <- yoy[,names(yoy) %in% names(dat_all)]
  yoy_juv$common_name = paste0(yoy_juv$common_name,"-juv")
  ### Bind the two data frames
  dat<-rbind(yoy_juv,dat_all)

  # format response
  dat$count <- as.numeric(dat$count)
  dat$count[which(is.na(dat$count))] <- 0
  dat$pres <- ifelse(dat$count > 0, 1, 0)

  # remove species with 0 records
  dat = dplyr::group_by(dat, sci_name) %>%
    dplyr::mutate(tot = sum(count)) %>%
    dplyr::filter(tot > 0) %>%
    dplyr::select(-tot)

  if (nrow(dat) > 0) {
    # expand predicted grid to have separate rows for each spp
    new_grid <-
      expand.grid(
        "common_name" = unique(dat$common_name),
        "station" = unique(pred_grid$station)
      )
    new_grid <- dplyr::left_join(new_grid, pred_grid) %>%
      dplyr::rename(species = common_name) %>%
      dplyr::filter(year >= min(dat$year), year <= max(dat$year))
    new_grid$jday <- median(dat$jday)

    # this is a pain, but with factors we can't make predictions to year-species
    # combinations with no data, so remove them from prediction grid
    zeros <- dplyr::group_by(dat, common_name, year) %>%
      dplyr::summarise(n = length(which(count>0)),
                       drop = ifelse(n>0,0,1)) %>%
      dplyr::rename(species = common_name)
    new_grid <- dplyr::left_join(zeros,new_grid) %>%
      dplyr::filter(drop==0) %>%
      dplyr::select(-n, -drop)

    # make sure year is a factor
    dat$year <- as.factor(dat$year)
    new_grid$year <- as.factor(new_grid$year)

    # nest fitted and predicted data
    dat_nested <-
      dat %>%
      dplyr::rename(species = common_name) %>%
      nest(-species) %>%
      rename(myorigdata = data)

    # create second dataset - predictions
    pred_nested <-
      new_grid %>%
      nest(-species) %>%
      rename(mynewdata = data)

    # fit presence-absence GAMs
    predictions_pres <-
      dat_nested %>%
      mutate(my_model = map(myorigdata, gam_pres_fit)) %>%
      full_join(pred_nested, by = "species") %>%
      mutate(my_new_pred = map2(my_model, mynewdata, predict)) %>%
      select(species, mynewdata, my_new_pred) %>%
      unnest(mynewdata, my_new_pred)# %>%
      #dplyr::group_by(species, year) %>%
      #dplyr::summarise(index = log(sum(exp(my_new_pred)))) %>%
      #as.data.frame()

    predictions_pos <-
      dat_nested %>%
      mutate(my_model = map(myorigdata, gam_pos_fit)) %>%
      full_join(pred_nested, by = "species") %>%
      mutate(my_new_pred = map2(my_model, mynewdata, predict)) %>%
      select(species, mynewdata, my_new_pred) %>%
      unnest(mynewdata, my_new_pred) #%>%
      #dplyr::group_by(species, year) %>%
      #dplyr::summarise(index = log(sum(exp(my_new_pred)))) %>%
      #as.data.frame()
    # also calculate summaries from data

    pred_all <- predictions_pos
    pred_all$my_new_pred <- exp(pred_all$my_new_pred) * plogis(predictions_pres$my_new_pred)
    pred_all <- dplyr::group_by(pred_all, species, year)%>%
      dplyr::summarise(index = log(sum(my_new_pred))) %>%
      as.data.frame()

    summaries <- dplyr::rename(dat, species = common_name) %>%
      dplyr::group_by(species, year) %>%
      dplyr::summarize(mean_catch = mean(count),
                       n_pos_catch = length(which(count > 0)))
    predictions_all <- dplyr::left_join(pred_all, summaries)

    if (i == 1) {
      all_pred <- predictions_all
    } else {
      all_pred <- rbind(all_pred, predictions_all)
    }
  }
}
all_pred$year = as.numeric(as.character(all_pred$year))
saveRDS(all_pred, "indices/predicted_indices.rds")

# Filter out experimental stations
# https://calcofi.org/field-work/station-positions.html
# stations <- read.csv("data/CalCOFIStationOrder.csv")
# stations <- dplyr::rename(stations, station = Station)
# dat <- dplyr::left_join(dat, stations[, c("station", "StaType")])
# dat <- dplyr::filter(dat, StaType == "ROS")

