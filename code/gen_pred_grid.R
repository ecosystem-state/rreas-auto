library(dplyr)
library(lubridate)
library(usethis)

source("code/set_control_params.R")
url_str <- "https://github.com/ecosystem-state/ecodata/blob/main/inst/rreas_index_data.rds"
usethis::use_github_file(url_str,
                         save_as = "data/raw_data.rds")

dat <- readRDS("data/raw_data.rds")
dat <- as.data.frame(dat)
dat$date <- lubridate::as_date(dat$time)
dat$year <- lubridate::year(dat$date)
dat$month <- lubridate::month(dat$date)
dat$yday <- lubridate::yday(dat$date)

spp <- unique(dat$sci_name)[1]
dat = dplyr::filter(dat,
                    sci_name == spp)
dat = dplyr::filter(dat, lat_dd < lat_max, lat_dd > lat_min)

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
