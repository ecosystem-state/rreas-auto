library(dplyr)
library(lubridate)

# settings
source("code/set_control_params.R")
#
rreas_erddap <- "erddap	species
FED_Rockfish_Catch	Diaphus theta
FED_Rockfish_Catch	Doryteuthis opalescens
FED_Rockfish_Catch	Engraulis mordax
FED_Rockfish_Catch	Euphausiacea
FED_Rockfish_Catch	Lampadena urophaos
FED_Rockfish_Catch	Merluccius productus
FED_Rockfish_Catch	Myctophidae
FED_Rockfish_Catch	Sebastes auriculatus
FED_Rockfish_Catch	Sebastes aurora
FED_Rockfish_Catch	Sebastes babcocki
FED_Rockfish_Catch	Symbolophorus californiensis
FED_Rockfish_Catch	Sebastes crameri
FED_Rockfish_Catch	Tarletonbeania crenularis
FED_Rockfish_Catch	Sebastes dallii
FED_Rockfish_Catch	Sebastes diploproa
FED_Rockfish_Catch	Sebastes elongatus
FED_Rockfish_Catch	Sebastes emphaeus
FED_Rockfish_Catch	Sebastes entomelas
FED_Rockfish_Catch	Sebastes flavidus
FED_Rockfish_Catch	Sebastes goodei
FED_Rockfish_Catch	Sebastes hopkinsi
FED_Rockfish_Catch	Sebastes jordani
FED_Rockfish_Catch	Sebastes levis
FED_Rockfish_Catch	Sebastes melanops
FED_Rockfish_Catch	Sebastes melanostomus
FED_Rockfish_Catch	Sebastes miniatus
FED_Rockfish_Catch	Sebastes mystinus
FED_Rockfish_Catch	Sebastes ovalis
FED_Rockfish_Catch	Sebastes paucispinis
FED_Rockfish_Catch	Sebastes pinniger
FED_Rockfish_Catch	Sebastes rastrelliger
FED_Rockfish_Catch	Sebastes ruberrimus
FED_Rockfish_Catch	Sebastes rufus
FED_Rockfish_Catch	Sardinops sagax
FED_Rockfish_Catch	Sebastes saxicola
FED_Rockfish_Catch	Sebastes semicinctus
FED_Rockfish_Catch	Sebastes serranoides
FED_Rockfish_Catch	Sebastes serriceps
FED_Rockfish_Catch	Sebastes spp.
FED_Rockfish_Catch	Sebastomus spp
FED_Rockfish_Catch	Sebastes wilsoni
FED_Rockfish_Catch	Sebastes zacentrus
FED_Rockfish_Catch	Sebastes spp. mel-flav complex
FED_Rockfish_Catch	Sebastes spp. caurinus complex"
rreas = read.table(textConnection(rreas_erddap), header=TRUE, sep="\t")

station_dat <- readRDS("data/raw_data.rds")

  station_dat <- as.data.frame(station_dat)
  station_dat$date <- lubridate::as_date(station_dat$time)
  station_dat$year <- lubridate::year(station_dat$date)
  station_dat$month <- lubridate::month(station_dat$date)
  station_dat$yday <- lubridate::yday(station_dat$date)

  # filter out recent years with consistent sampling
  #station_dat <- dplyr::filter(station_dat, year >= min_year)
  # format response
station_dat$catch <- as.numeric(station_dat$catch)
station_dat$catch[which(is.na(station_dat$larvae_10m2))] <- 0

station_dat$name <- paste0(station_dat$common_name, "-",station_dat$maturity)
station_dat$file = rreas$erddap[1]
station_dat <- dplyr::group_by(station_dat, name) %>%
  dplyr::summarize(tot_cpue = sum(catch),
                   sci_name = sci_name[1],
                     erddap = file[1])

dat <- station_dat
dat <- dplyr::filter(dat, sci_name %in% rreas$species)

saveRDS(dat, "indices/tot_cpue_species.rds")
