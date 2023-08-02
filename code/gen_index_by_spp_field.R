library(mgcv)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(usethis)

# settings
source("code/set_control_params.R")

# url_str <- "https://github.com/ecosystem-state/ecodata/blob/main/inst/rreas_index_data.rds"
# usethis::use_github_file(url_str,
#                          save_as = "data/raw_data.rds")
# url_str <- "https://github.com/ecosystem-state/ecodata/blob/main/inst/tot_cpue_rreas.rds"
# usethis::use_github_file(url_str,
#                          save_as = "indices/tot_cpue_species.rds")
#
# # load species tot cpue
# tot_cpue <- readRDS("indices/tot_cpue_species.rds")
#tot_cpue <- dplyr::arrange(tot_cpue,-tot_cpue)
#tot_cpue <- tot_cpue[1:top_species,]

# load prediction grid
pred_grid <- readRDS("indices/pred_grid.rds")

# model function - edit to change form of model
gam_pres_fit <- function(df) {
  gam(pres ~ jday + s(latitude, longitude) + as.factor(year) +I(jday^2),
      data = df,
      family = "binomial"
  )
}

gam_pos_fit <- function(df) {
  if(length(unique(df$year[which(df$pres == 1)])) > 1) {
  gam(count ~ jday + s(latitude, longitude) + as.factor(year) +I(jday^2),
      data = df[which(df$count > 0),],
      family = "poisson"
  )
  } else {
    gam(count ~ 1,
        data = df[which(df$count > 0),],
        family = "poisson"
    )
  }
}

# grab data for all species
dat_old <- readRDS("data/raw_data.rds")

dat <- read.csv("data/rreas_mary_eric_20230308.csv")
names(dat) <- tolower(names(dat))

name <- names(dat)[18:36]
t <- matrix(0, length(name), length(unique(dat$year)))
years <- sort(unique(dat$year))

for(i in 1:length(name)) {
  for(j in 1:length(years)) {
    indx <- dat[which(dat$year==years[j]), name[i]]
    indx <- sum(ceiling(indx/1.0e10))
    t[i,j] <- indx
  }
}
write.csv(t,"data_summary.csv")


dat <- sdmTMB::add_utm_columns(dat, ll_names = c("londd","latdd"))
dat <- dplyr::rename(dat, latitude=Y, longitude=X)

rreas_stations <- dplyr::group_by(dat, station) %>%
    dplyr::summarise(latitude = latitude[1],
                     longitude = longitude[1])
write.csv(rreas_stations, "data/rreas_stations_field.csv", row.names = FALSE)

# generate pred grid
dat$floor_lon <- floor(dat$longitude / resolution)
dat$floor_lat <- floor(dat$latitude / resolution)
dat$station <- paste(dat$floor_lon, dat$floor_lat)

pred_grid <- expand.grid(
  station = unique(dat$station),
  year = unique(dat$year)
)

  # convert date
  # dat <- as.data.frame(dat)
  # dat$date <- lubridate::as_date(dat$time)
  # dat$year <- lubridate::year(dat$date)
  # dat$month <- lubridate::month(dat$date)
  # dat$jday <- lubridate::yday(dat$date)

  dat$uniqueID<-do.call(paste, c(dat[c("cruise", "haul_no")], sep = "_"))

  # sum total counts for spp we don't care about maturity for
  # dat$uniqueID<-do.call(paste, c(dat[c("cruise", "haul_no")], sep = "_"))
  # spp_all<-c('Sardinops sagax','Doryteuthis opalescens','Euphausiacea')
  # dat_allStages<-dat[dat$sci_name %in% spp_all,]
  # dat_all<-dat_allStages %>%
  #   group_by(year,uniqueID,latitude,longitude,jday,sci_name,common_name) %>%
  #   summarise(count=sum(catch))%>%
  #   data.frame

  # ## For those species that we want only juvenile stages (or only juv stages identified)
  # #keep those data where maturity == "Y"
  # spp_juvenile<-c('Sardinops sagax','Engraulis mordax','Citharichthys stigmaeus',
  #                 'Merluccius productus','Sebastes goodei','Sebastes entomelas',
  #                 'Sebastes jordani','Sebastes mystinus','Sebastes paucispinis',
  #                 'Sebastes semicinctus','Citharichthys sordidus')
  # spp_extra <- c("Sebastes melanostomus", "Sebastes wilsoni","Sebastes pinniger",
  #                "Sebastes auriculatus","Sebastes flavidus","Sebastes spp. caurinus complex",
  #                "Sebastes hopkinsi","Sebastes saxicola")
  # yoy <- dplyr::filter(dat, sci_name%in% c(spp_extra, spp_juvenile),
  #                      maturity=="Y") %>%
  #       dplyr::rename(count = catch)
  #
  # yoy_juv <- yoy[,names(yoy) %in% names(dat_all)]
  # yoy_juv$common_name = paste0(yoy_juv$common_name,"-juv")
  # ### Bind the two data frames
  # dat<-rbind(yoy_juv,dat_all)

  # wide to long conversion
  dat_long <- dat %>%
    pivot_longer(
      cols = c(adult_anchovy, adult_sardine, market_squid,
               octopus, total_krill, total_myctophids,
               yoy_anchovy, yoy_bocaccio, yoy_chilipepper,
               yoy_halfbanded, yoy_lingcod, yoy_pac_sanddab,
               yoy_pacific_hake, yoy_rockfish, yoy_sardine,
               yoy_shortbelly, yoy_speckled_sanddab, yoy_widow,
               yoy_yellowtail), # columns you want to gather into a single column
      names_to = "species", # new column for gathered column names
      values_to = "values" # new column for gathered column values
    )
  dat <- dat_long

  # format response
  dat$count <- as.numeric(dat$values)
  dat$count[which(is.na(dat$count))] <- 0
  dat$pres <- ifelse(dat$count > 0, 1, 0)

  # remove species with 0 records
  dat = dplyr::group_by(dat, species) %>%
    dplyr::mutate(tot = sum(count)) %>%
    dplyr::filter(tot > 0) %>%
    dplyr::select(-tot)

if (nrow(dat) > 0) {
    # expand predicted grid to have separate rows for each spp
    new_grid <-
      expand.grid(
        "common_name" = unique(dat$species),
        "station" = unique(pred_grid$station)
      )
    new_grid <- dplyr::left_join(new_grid, pred_grid) %>%
      dplyr::rename(species = common_name) %>%
      dplyr::filter(year >= min(dat$year), year <= max(dat$year))
    new_grid$jday <- median(dat$jday)
    new_grid$longitude <- as.numeric(unlist(lapply(strsplit(as.character(new_grid$station), " "), getElement, 1)))
    new_grid$latitude <- as.numeric(unlist(lapply(strsplit(as.character(new_grid$station), " "), getElement, 2)))

    # this is a pain, but with factors we can't make predictions to year-species
    # combinations with no data, so remove them from prediction grid
    zeros <- dplyr::group_by(dat, species, year) %>%
      dplyr::summarise(n = length(which(count>0)),
                       drop = ifelse(n>0,0,1))
    new_grid <- dplyr::left_join(zeros,new_grid) %>%
      dplyr::filter(drop==0) %>%
      dplyr::select(-n, -drop)

    # make sure year is a factor
    dat$year <- as.factor(dat$year)
    new_grid$year <- as.factor(new_grid$year)

    # nest fitted and predicted data
    dat_nested <-
      dat %>%
      dplyr::group_by(species) %>%
      nest()
    dat_nested <- dat_nested %>% dplyr::rename(myorigdata = data)

    # create second dataset - predictions
    pred_nested <-
      group_by(new_grid, species) %>%
      tidyr::nest() %>%
      dplyr::rename(mynewdata = data)

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

    summaries <- dat %>%
      dplyr::group_by(species, year) %>%
      dplyr::summarize(mean_catch = mean(count),
                       n_pos_catch = length(which(count > 0)))
    predictions_all <- dplyr::left_join(pred_all, summaries)


    all_pred <- predictions_all
  }

all_pred$year = as.numeric(as.character(all_pred$year))
saveRDS(all_pred, "indices/predicted_indices_field.rds")

# Filter out experimental stations
# https://calcofi.org/field-work/station-positions.html
# stations <- read.csv("data/CalCOFIStationOrder.csv")
# stations <- dplyr::rename(stations, station = Station)
# dat <- dplyr::left_join(dat, stations[, c("station", "StaType")])
# dat <- dplyr::filter(dat, StaType == "ROS")
