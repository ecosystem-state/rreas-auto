library(sdmTMB)
library(rerddap)
library(dplyr)
library(lubridate)
library(sf)
library(tidyr)
library(purrr)

# settings
source("code/set_control_params.R")

# load prediction grid
pred_grid <- readRDS("indices/pred_grid.rds")

# load species tot cpue
dat <- readRDS("data/raw_data.rds")

rreas_stations <- dplyr::group_by(dat, station) %>%
  dplyr::summarise(latitude = latitude[1],
                   longitude = longitude[1])
write.csv(rreas_stations, "data/rreas_stations.csv", row.names = FALSE)

# convert date
dat <- as.data.frame(dat)
dat$date <- lubridate::as_date(dat$time)
dat$year <- lubridate::year(dat$date)
dat$month <- lubridate::month(dat$date)
dat$jday <- lubridate::yday(dat$date)

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
    new_grid$fyear <- as.factor(new_grid$year)

    for(spp in 1:length(unique(dat$common_name))) {
      sub <- dplyr::filter(dat, common_name == unique(dat$common_name)[spp])
      if(spp==1) {
        mesh = make_mesh(sub, xy_cols = c("longitude","latitude"),
                       cutoff = 50)
      }
      sub$fyear = as.factor(sub$year)
      # Tweedie converges fine, delta_ models do not
      sub$year <- as.numeric(as.character(sub$year))
      fit_spde_tw <- sdmTMB(count ~ -1 + fyear + s(jday,k=3),
                   spatiotemporal = "off",
                   time="year",
                   spatial="on",
                   family = tweedie(),
                   mesh=mesh,
                   data=sub)

      fit_spde_delta <- sdmTMB(count ~ -1 + fyear + s(jday,k=3),
                              spatiotemporal = "off",
                              time="year",
                              spatial="on",
                              family = delta_truncated_nbinom2(),
                              mesh=mesh,
                              data=sub)

      fit_spde_nb2 <- sdmTMB(count ~ -1 + fyear + s(jday,k=3),
                                 spatiotemporal = "off",
                                 time="year",
                                 spatial="on",
                                 family = nbinom2(),
                                 mesh=mesh,
                                 data=sub)

      fit_smooth_tw <- sdmTMB(count ~ -1 + s(longitude, latitude) + fyear + s(jday,k=3),
                              spatiotemporal = "off",
                              time="year",
                              spatial="off",
                              family = tweedie(),
                              mesh=mesh,
                              data=sub)

      fit_smooth_nb2 <- sdmTMB(count ~ -1 + s(longitude, latitude) + fyear + s(jday,k=3),
                              spatiotemporal = "off",
                              time="year",
                              spatial="off",
                              family = nbinom2(),
                              mesh=mesh,
                              data=sub)

      fit_smooth_delta <- sdmTMB(count ~ -1 + s(longitude, latitude) + fyear + s(jday,k=3),
                               spatiotemporal = "off",
                               time="year",
                               spatial="off",
                               family = delta_truncated_nbinom2(),
                               mesh=mesh,
                               data=sub)
      df <- tidy(fit_smooth_nb2)
      df$model <- "smooth_nb2"
      df2 <- tidy(fit_smooth_tw)
      df2$model <- "smooth_tw"
      df3 <- tidy(fit_spde_nb2)
      df3$model <- "spde_nb2"
      df4 <- tidy(fit_spde_tw)
      df4$model <- "spde_tw"
      df_all <- rbind(df, df2, df3, df4)
      df_all$species <- sub$common_name[1]
      df_all$term <- substr(df_all$term,6,9)

      if(spp == 1) {
        coefs <- df_all
      } else {
        coefs <- rbind(coefs, df_all)
      }
    }
    write.csv(coefs,"coefs_for_swfsc.csv",row.names = FALSE)

}

# Filter out experimental stations
# https://calcofi.org/field-work/station-positions.html
# stations <- read.csv("data/CalCOFIStationOrder.csv")
# stations <- dplyr::rename(stations, station = Station)
# dat <- dplyr::left_join(dat, stations[, c("station", "StaType")])
# dat <- dplyr::filter(dat, StaType == "ROS")

