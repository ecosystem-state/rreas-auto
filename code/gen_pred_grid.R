library(dplyr)
library(lubridate)
library(sf)

source("code/set_control_params.R")

out <- info(as.character(erddap[i]))

station_dat <- readRDS("data/raw_data.rds")

dat <- station_dat

dat <- as.data.frame(dat)
dat$date <- lubridate::as_date(dat$time)
dat$year <- lubridate::year(dat$date)
dat$month <- lubridate::month(dat$date)
dat$yday <- lubridate::yday(dat$date)

spp <- unique(dat$sci_name)[1]
dat = dplyr::filter(dat,
                    sci_name == spp)
dat$latitude <- as.numeric(dat$latitude)
dat$longitude <- as.numeric(dat$longitude)
dat = dplyr::filter(dat, latitude < lat_max, latitude > lat_min)

# convert to UTM - kms
# make the UTM cols spatial (X/Easting/lon, Y/Northing/lat)
dat <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326)
# transform to UTM
dat<- st_transform(x = dat, crs = 32610)
dat$longitude = st_coordinates(dat)[,1]
dat$latitude = st_coordinates(dat)[,2]

dat <- as.data.frame(dat)
dat$longitude <- dat$longitude / 1000
dat$latitude <- dat$latitude / 1000

# come up with prediction grid
resolution <- pred_resolution
dat$floor_lon <- floor(dat$longitude / resolution)
dat$floor_lat <- floor(dat$latitude / resolution)
dat$station <- paste(dat$floor_lon, dat$floor_lat)

pred_grid <- expand.grid(
  station = unique(dat$station),
  year = unique(dat$year)
)

station_df <- data.frame(station = unique(dat$station))
station_df$longitude <- as.numeric(unlist(lapply(strsplit(as.character(station_df$station), " "), getElement, 1))) * resolution
station_df$latitude <- as.numeric(unlist(lapply(strsplit(as.character(station_df$station), " "), getElement, 2))) * resolution
pred_grid <- dplyr::left_join(pred_grid, station_df)

saveRDS(pred_grid, "indices/pred_grid.rds")