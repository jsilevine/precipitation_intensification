library(ggplot2)
library(raster)
library(cowplot)
library(measurements)
library(mgcv)
library(sf)
library(evgam)
library(lqmm)
library(ncdf4)
library(foreach)
library(brms)
library(metR)
library(ggdist)

xylem <- read.csv("data/xylem_traits.csv")
##---------------------------------------------------------------
## xylem
##---------------------------------------------------------------

xylem$P50..Mpa. <- as.numeric(xylem$P50..MPa.)

for (i in 1:nrow(xylem)) {

  nlat <- xylem[i, "Latitude"]
  nlat <- gsub("'", "", nlat)
  nlat <- gsub("N", "", nlat)
  if (!grepl("\"", nlat)) {
    nlat <- paste0(nlat, "00")
  } else {
    nlat <- gsub("\"", "", nlat)
  }
  xylem[i, "Latitude"] <- as.numeric(conv_unit(nlat, from = "deg_min_sec", to = "dec_deg"))

  nlon <- xylem[i, "Longitude"]
  nlon <- gsub("'", "", nlon)
  nlon <- gsub("N", "", nlon)
  if (!grepl("\"", nlon)) {
    nlon <- paste0(nlon, "00")
  } else {
    nlon <- gsub("\"", "", nlon)
  }
  xylem[i, "Longitude"] <- as.numeric(conv_unit(nlon, from = "deg_min_sec", to = "dec_deg"))

}

sites <- unique(xylem$Site)
site_data <- data.frame(site = sites)
site_data$site_name <- NA
site_data$suf_data <- 0
site_data$map <- NA
site_data$lat <- NA
site_data$lon <- NA
site_data$range_P50 <- NA
site_data$min_P50 <- NA
site_data$max_P50 <- NA
site_data$nspp <- NA
site_data$var_PPT <- NA
for (i in 1:nrow(site_data)) {
  if (length(unique(xylem[xylem$Site == site_data[i, "site"], "Cleaned.binomial"])) > 3) {
    site_data[i, "suf_data"] <- 1
    site_data[i, "nspp"] <- length(unique(xylem[xylem$Site == site_data[i, "site"], "Cleaned.binomial"]))
    site_data[i, "map"] <- xylem[xylem$Site == site_data[i, "site"], "PPTbest"][1]
    site_data[i, "site_name"] <- xylem[xylem$Site == site_data[i, "site"], "site.name"][1]
    site_data[i, "lat"] <- xylem[xylem$Site == site_data[i, "site"], "LATnew"][1]
    site_data[i, "lon"] <- xylem[xylem$Site == site_data[i, "site"], "LONnew"][1]
    site_data[i, "range_P50"] <- range(as.numeric(xylem[xylem$Site == site_data[i, "site"], "P50..MPa."]), na.rm = TRUE)[2] -
      range(as.numeric(xylem[xylem$Site == site_data[i, "site"], "P50..MPa."]), na.rm = TRUE)[1]
    site_data[i, "min_P50"] <- min(as.numeric(xylem[xylem$Site == site_data[i, "site"], "P50..MPa."]), na.rm = TRUE)
    site_data[i, "max_P50"] <- max(as.numeric(xylem[xylem$Site == site_data[i, "site"], "P50..MPa."]), na.rm = TRUE)
    site_data[i, "var_PPT"] <- xylem[xylem$Site == site_data[i, "site"], 115:126] 

  }
}

site_data <- site_data[site_data$suf_data == 1, ]
site_data <- site_data[site_data$range_P50 != -Inf,]
site_data <- site_data[!is.na(site_data$map),]


ggplot(site_data, aes(x = map, y = range_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 4)) +
  theme_bw()

ggplot(site_data, aes(x = map, y = range_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  geom_smooth() +
  theme_bw()

ggplot(site_data, aes(x = map, y = min_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  geom_smooth() +
  theme_bw()

ggplot(site_data, aes(x = map, y = range_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  theme_bw()

ggplot(site_data, aes(x = map, y = max_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  #geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  theme_bw()

ggplot(site_data, aes(x = var_PPT / map, y = range_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  #geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  theme_bw()

ggplot(site_data, aes(x = var_PPT / map, y = min_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  #geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  theme_bw()

ggplot(site_data, aes(x = map, y = min_P50, weight = nspp)) +
  geom_point(aes(size = nspp)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  theme_bw()

ggplot(site_data, aes(x = map, y = nspp)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

xylem$PPTbest <- as.numeric(xylem$PPTbest)
xylem[,"P50..MPa."] <- as.numeric(xylem[,"P50..MPa."])
xylem[,"Ks..kg.m.1.MPa.1.s.1."] <- as.numeric(xylem[,"Ks..kg.m.1.MPa.1.s.1."])

ggplot(xylem, aes(x = PPTbest, y = log(Ks..kg.m.1.MPa.1.s.1.))) +
  geom_point() +
  theme_bw()

ggplot(xylem, aes(x = PPTbest, y = -1*log(abs(P50..MPa.)))) +
  geom_point() +
  theme_bw()

ggplot(xylem, aes(x = PPTbest, y = P50..MPa.)) +
  geom_point() +
  theme_bw()

summary(lm(range_P50 ~ poly(map, 4), data = site_data))
gam_data <- site_data
gm <- gam(range_P50 ~ s(map, k = 7), data = gam_data, weights = gam_data$nspp)
summary(gm)
AIC(gm)
lm <- lm(range_P50 ~ map, data = gam_data, weights = gam_data$nspp)
AIC(lm)
lm2 <- lm(range_P50 ~ poly(map, 2), data = gam_data, weights = gam_data$nspp)
AIC(lm2)
lm3 <- lm(range_P50 ~ poly(map, 3), data = gam_data, weights = gam_data$nspp)
AIC(lm3)
lm4 <- lm(range_P50 ~ poly(map, 4), data = gam_data, weights = gam_data$nspp)
AIC(lm4)

gm <- gam(min_P50 ~ s(map, k = 7), data = gam_data, weights = gam_data$nspp)
summary(gm)
plot(gm)
AIC(gm)
lm <- lm(min_P50 ~ map, data = gam_data, weights = gam_data$nspp)
AIC(lm)
lm2 <- lm(min_P50 ~ poly(map, 2), data = gam_data, weights = gam_data$nspp)
AIC(lm2)
lm3 <- lm(min_P50 ~ poly(map, 3), data = gam_data, weights = gam_data$nspp)
AIC(lm3)
lm4 <- lm(min_P50 ~ poly(map, 4), data = gam_data, weights = gam_data$nspp)
AIC(lm4)

gm <- gam(max_P50 ~ s(map, k = 7), data = gam_data, weights = gam_data$nspp)
summary(gm)
plot(gm)
AIC(gm)
lm <- lm(max_P50 ~ map, data = gam_data, weights = gam_data$nspp)
AIC(lm)
lm2 <- lm(max_P50 ~ poly(map, 2), data = gam_data, weights = gam_data$nspp)
AIC(lm2)
lm3 <- lm(range_P50 ~ poly(map, 3), data = gam_data, weights = gam_data$nspp)
AIC(lm3)
lm4 <- lm(range_P50 ~ poly(map, 4), data = gam_data, weights = gam_data$nspp)
AIC(lm4)

pdata <- data.frame(map = seq(min(site_data$map), max(site_data$map), length.out = 100))
p <- predict(gm, pdata, "response", TRUE)
pdata$range_P50 <- p$fit
pdata$upper <- pdata$range_P50 + 1.97 * p$se.fit
pdata$lower <- pdata$range_P50 - 1.97 * p$se.fit
pdata[pdata$lower < 0.0, "lower"] <- 0.0

ggplot(site_data, aes(x = map, y = range_P50)) +
  geom_point(aes(size = nspp)) +
  geom_ribbon(data = pdata, aes(ymin = lower, ymax = upper), alpha = 0.1)+
  geom_line(data = pdata, color = "#3182bd", size = 2.5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0.0, 0.05), limits = c(0.0, NA)) +
  theme_bw() +
  ggtitle("Hydraulic trait diversity vs. mean precipitation") +
  labs(size = "# species") +
  ylab("Range P50 (Mpa)") +
  xlab("Mean annual precipitation (mm)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("hydraulic_traits/figures/trait_diversity.pdf")

world <- st_read("hydraulic_traits/data/World_Countries_Generalized.shp")
world <- st_transform(world, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
loc <- st_as_sf(site_data, coords = c("lon", "lat"), crs = 4326)

ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = map), size = 5) +
  theme_bw() +
  labs(color = "Mean annual \n precipitation") +
  scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Xylem trait data origins")
ggsave("hydraulic_traits/figures/trait_locations.pdf")


##---------------------------------------------------------------
## Extreme value analysis of raw trait data
##---------------------------------------------------------------

xylem <- read.csv("data/xylem_traits.csv")

traits_full <- xylem[xylem[,"Natural...Greenhouse"] == "N" &
                    xylem[,"Plant.organ"] %in% c("S", "P") &
                    xylem[,"Developmental.stage"] == "A", ]

traits_full[!(traits_full$P50.method %in% c("DH", "CE", "CA", "AD", "AS", "Pn")), "P50"] <- NA

## clean latitude and longitude data
for (i in 1:nrow(traits_full)) {

  if (is.na(traits_full[i, "LATnew"]) & !is.na(traits_full[i, "Latitude"])) {

    nlat <- traits_full[i, "Latitude"]
    nlat <- gsub("'", "", nlat)
    if (grepl("S", nlat)) nlat <- paste0("-", nlat)
    nlat <- gsub("N", "", nlat)
    nlat <- gsub("S", "", nlat)
    if (nchar(nlat) == 5) {
      nlat <- paste0(nlat, " 00")
    } else if (!grepl("\"", nlat)) {
      nlat <- paste0(nlat, "00")
    } else {
      nlat <- gsub("\"", "", nlat)
    }
    traits_full[i, "LATnew"] <- as.numeric(conv_unit(nlat, from = "deg_min_sec", to = "dec_deg"))

    nlon <- traits_full[i, "Longitude"]
    nlon <- gsub("'", "", nlon)
    if (grepl("W", nlon)) nlon <- paste0("-", nlon)
    nlon <- gsub("W", "", nlon)
    nlon <- gsub("E", "", nlon)
    if (nchar(nlon) == 5) {
      nlon <- paste0(nlon, " 00")
    } else if (!grepl("\"", nlon)) {
      nlon <- paste0(nlon, "00")
    } else {
      nlon <- gsub("\"", "", nlon)
    }
    traits_full[i, "LONnew"] <- as.numeric(conv_unit(nlon, from = "deg_min_sec", to = "dec_deg"))

  }
}

tokeep <- c("Cleaned.family", "Cleaned.genus", "Cleaned.species", "Cleaned.binomial", "P50..MPa.", "P12..MPa.", "P88..MPa.",
            "Ks..kg.m.1.MPa.1.s.1.", "Location", "Country", "LATnew", "LONnew", "MAP..mean.annual.precipitation..mm.",
            "MAT..mean.annual.temperature...C..", "Biome", "site.name", "PPTbest", "PPTjan", "PPTfeb", "PPTmar", "PPTapr",
            "PPTmay", "PPTjun", "PPTjul", "PPTaug", "PPTsep", "PPToct", "PPTnov", "PPTdec", "PPTcru", "NewID", "Amax..micromol.m2.s.")

traits_full <- traits_full[, tokeep]
newnames <- c("family", "genus", "species", "binomial", "P50", "P12", "P88",
            "Ks", "location", "country", "lat", "lon", "map_old",
            "mat", "biome", "site_name", "map", "map_jan", "map_feb", "map_mar", "map_apr",
            "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru", "id", "a_m")
colnames(traits_full) <- newnames

for (i in 1:nrow(traits_full)) {
  if (is.na(traits_full[i, "map"]) && !is.na(traits_full[i, "map_old"])) {
    traits_full[i, "map"] <- traits_full[i, "map_old"]
  }
}

traits_full <- traits_full[, colnames(traits_full) != "map_old"]

tonum <- c( "P50", "P12", "P88",
            "Ks", "lat", "lon",
            "mat", "map", "map_jan", "map_feb", "map_mar", "map_apr",
            "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru",
           "a_m")

for (c in tonum) {
  traits_full[,c] <- as.numeric(traits_full[,c])
}

traits_full$sd_map <- NA
for (i in 1:nrow(traits_full)) {
  traits_full[i, "sd_map"] <- sd(traits_full[i, 17:28], na.rm = TRUE)
}
traits_full$sd_map <- traits_full$sd_map / traits_full$map

traits_full$log_P50 <- -log(abs(traits_full$P50))
traits_full$log_Ks <- log(traits_full$Ks)

##---------------------------------------------------------------
## Reanalysis processing
##---------------------------------------------------------------

traits_full <- traits_full[!is.na(traits_full$P50 | !is.na(traits_full$Ks)),]


traits_full$latlon <- NA
for (i in 1:nrow(traits_full)) {
  print(i)
  if (!is.na(traits_full[i, "lat"])) {
    traits_full[i, "latlon"] <- paste0(traits_full[i, "lat"], traits_full[i, "lon"])
  }
}


site_data <- data.frame(latlon = unique(traits_full[!is.na(traits_full$latlon), "latlon"]))
for (i in 1:nrow(site_data)) {
  cordata <- traits_full[traits_full$latlon == site_data[i, "latlon"], ]
  cordata <- cordata[!is.na(cordata$lat),]
  site_data[i, "lat"] <- cordata$lat[1]
  site_data[i, "lon"] <- cordata$lon[1]
  site_data[i, "site_name"] <- cordata$site_name[1]
  site_data[i, "id"] <- cordata$id[1]
  site_data[i, "map"] <- cordata$map[1]
}

world <- st_read("hydraulic_traits/data/World_Countries_Generalized.shp")
world <- st_transform(world, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
loc <- st_as_sf(site_data[!is.na(site_data$lon),], coords = c("lon", "lat"), crs = 4326)

ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = map), size = 5) +
  theme_bw() +
  labs(color = "Mean annual \n precipitation") +
  scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Xylem trait data origins")

ggsave("hydraulic_traits/figures/trait_locations.pdf")



mean_storm_size <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
median_storm_size <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
sd_storm_size <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
storm_freq <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
map_reanalysis <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
mean_temp <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
median_temp <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)
sd_temp <- matrix(rep(0, 110*nrow(site_data)), ncol = 110)


leap_years <- seq(1904, 2008, by = 4)

clim_files <- list.files("hydraulic_traits/data/reanalysis/")
pr_files <- clim_files[grepl("pr", clim_files)]
tas_files <- clim_files[grepl("tas", clim_files)]

## growing season files
gs_first <- raster("hydraulic_traits/data/gs_first.tif")
gs_last <- raster("hydraulic_traits/data/gs_last.tif")

site_data$gs_first <- extract(gs_first, site_data[,c("lon", "lat")])
site_data$gs_last <- extract(gs_last, site_data[,c("lon", "lat")])

site_data$gs_hem <- "N"
site_data[site_data$gs_last < site_data$gs_first & !is.na(site_data$gs_first), "gs_hem"] <- "S"
site_data[is.na(site_data$gs_first), "gs_hem"] <- "O"

cfun <- function(x, y) {
  for (i in 1:length(x)) {
    x[[i]] <- cbind(x[[i]], y[[i]])
  }
  return(x)
}

library(doParallel)

cl <- makeCluster(6, outfile = "")
registerDoParallel(cl)

clim_data <- foreach(f = 1:length(pr_files),
                     .combine = 'cfun',
                     .packages = c("ncdf4")) %dopar% {

  mean_storm_size <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  median_storm_size <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  sd_storm_size <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  storm_freq <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  map_reanalysis <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  mean_storm_size_fy <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  median_storm_size_fy <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  sd_storm_size_fy <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  storm_freq_fy <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  map_reanalysis_fy <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  mean_temp <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  median_temp <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)
  sd_temp <- matrix(rep(0, 10*nrow(site_data)), ncol = 10)

  pr_data <- nc_open(paste0("hydraulic_traits/data/reanalysis/", pr_files[f]))
  tas_data <- nc_open(paste0("hydraulic_traits/data/reanalysis/", tas_files[f]))

  lon <- ncvar_get(pr_data, "lon")
  lat <- ncvar_get(pr_data, "lat")
  time <- ncvar_get(pr_data, "time")

  for (i in 1:nrow(site_data)) {

    print(paste0("starting row: ", i, " of ", nrow(site_data), " \n for file: ", f))

    x <- as.numeric(site_data[i, c("lon", "lat")])

    latindex = which(abs(lat - x[2]) == min(abs(lat - x[2])))[1]
    lonindex = which(abs(lon - x[1]) == min(abs(lon - x[1])))[1]
    start <- c(lonindex, latindex, 1)
    count <- c(1, 1, -1)

    pr <- as.numeric(ncvar_get(pr_data, varid = "pr", start = start, count))
    pr <- pr * 24 * 60 * 60 ## convert to mm
    tas <- as.numeric(ncvar_get(tas_data, varid = "tas", start = start, count))
    tas <- tas - 273.15

    year <- as.numeric(strsplit(pr_files[f], "_")[[1]][2])
    rstart <- 1

    for (y in 1:10) {

      if (year %in% leap_years) {
        range <- rstart:(rstart+365)
        ed <- 366
      } else {
        range <- rstart:(rstart+364)
        ed <- 365
      }

      if (site_data[i, "gs_hem"] == "O") {
        subrange <- range
      } else if (site_data[i, "gs_hem"] == "N") {
        subrange <- range[site_data[i, "gs_first"]:site_data[i, "gs_last"]]
      } else if (site_data[i, "gs_hem"] == "S") {
        subrange <- range[c(site_data[i, "gs_first"]:ed, 1:site_data[i, "gs_last"])]
      }

      tas_sub <- tas[subrange]

      mean_temp[i,y] <- mean(tas_sub)
      median_temp[i,y] <- median(tas_sub)
      sd_temp[i,y] <- sd(tas_sub)

      pr_sub <- pr[subrange]
      storms <- pr_sub[pr_sub > 2.5]

      print(length(storms) / (length(subrange) / 30))

      mean_storm_size[i, y] <- mean(storms)
      median_storm_size[i, y] <- median(storms)
      sd_storm_size[i, y] <- sd(storms)
      storm_freq[i, y] <- length(storms) / (length(subrange) / 30)
      map_reanalysis[i, y] <- sum(pr[subrange])

      pr_sub <- pr[range]
      storms <- pr_sub[pr_sub > 2.5]

      mean_storm_size_fy[i, y] <- mean(storms)
      median_storm_size_fy[i, y] <- median(storms)
      sd_storm_size_fy[i, y] <- sd(storms)
      storm_freq_fy[i, y] <- length(storms) / (length(range) / 30)
      map_reanalysis_fy[i, y] <- sum(pr[range])

      year <- year+1
      rstart <- range[length(range)]+1

    }
  }

  list(mean_storm_size, median_storm_size, sd_storm_size, storm_freq, map_reanalysis, mean_temp, median_temp, sd_temp,
       mean_storm_size_fy, median_storm_size_fy, sd_storm_size_fy, storm_freq_fy, map_reanalysis_fy)

}

saveRDS(clim_data, "hydraulic_traits/data/reanalysis/clim_data_24.rds")
clim_data <- readRDS("hydraulic_traits/data/reanalysis/clim_data_24.rds")

site_data$mean_storm_size <- rowMeans(clim_data[[1]], na.rm = TRUE)
site_data$median_storm_size <- rowMeans(clim_data[[2]], na.rm = TRUE)
site_data$sd_storm_size <- rowMeans(clim_data[[3]], na.rm = TRUE)
site_data$storm_freq <- rowMeans(clim_data[[4]], na.rm = TRUE)
site_data$map_reanalysis <- rowMeans(clim_data[[5]], na.rm = TRUE)
site_data$mean_temp <- rowMeans(clim_data[[6]], na.rm = TRUE)
site_data$median_temp <- rowMeans(clim_data[[7]], na.rm = TRUE)
site_data$sd_temp <- rowMeans(clim_data[[8]], na.rm = TRUE)
site_data$mean_storm_size_fy <- rowMeans(clim_data[[9]], na.rm = TRUE)
site_data$median_storm_size_fy <- rowMeans(clim_data[[10]], na.rm = TRUE)
site_data$sd_storm_size_fy <- rowMeans(clim_data[[11]], na.rm = TRUE)
site_data$storm_freq_fy <- rowMeans(clim_data[[12]], na.rm = TRUE)
site_data$map_reanalysis_fy <- rowMeans(clim_data[[13]], na.rm = TRUE)


world <- st_read("hydraulic_traits/data/World_Countries_Generalized.shp")
world <- st_transform(world, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
loc <- st_as_sf(site_data[!is.na(site_data$lon),], coords = c("lon", "lat"), crs = 4326)

p1 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = map_reanalysis), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Mean annual \n precipitation") +
  scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Mean annual precipitation (mm)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")

p2 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = median_storm_size), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Median storm size") +
  scale_color_gradient2(low = "#ccece6", mid = "#66c2a4", high = "#006d2c", midpoint = 5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Median storm size (mm)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")

p3 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = storm_freq), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Storm frequency \n (n. per year)") +
  scale_color_gradient2(low = "#ccece6", mid = "#66c2a4", high = "#006d2c", midpoint = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Storm frequency (n. per year)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")

p4 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = loc, aes(color = log(median_storm_size / storm_freq)), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "log(storm size / \n storm freq.)") +
  scale_color_gradient2(low = "#fcc5c0", mid = "#f768a1", high = "#7a0177", midpoint = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Storm size / storm frequency, or the 'few large vs many small' index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")
p4

plot_grid(p1, p3, p4, ncol = 1)
ggsave("hydraulic_traits/figures/figure1.pdf", width = 12, height = 20)




traits_full <- merge(traits_full, site_data, by = c("lat", "lon"))
traits_full <- traits_full[,c("family", "genus", "species", "binomial", "P50", "P12", "P88",
                              "Ks", "location", "country", "lat", "lon",
                              "mat", "biome", "site_name.x", "map.x", "map_jan", "map_feb", "map_mar", "map_apr",
                              "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru",
                              "sd_map", "log_P50", "log_Ks", "mean_storm_size", "median_storm_size", "sd_storm_size", "storm_freq",
                              "map_reanalysis", "a_m", "mean_temp")]

colnames(traits_full) <- c("family", "genus", "species", "binomial", "P50", "P12", "P88",
                              "Ks", "location", "country", "lat", "lon",
                              "mat", "biome", "site_name", "map", "map_jan", "map_feb", "map_mar", "map_apr",
                              "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru",
                              "sd_map", "log_P50", "log_Ks", "mean_storm_size", "median_storm_size", "sd_storm_size", "storm_freq",
                              "map_reanalysis", "a_m", "mean_temp")
write.csv(traits_full, "hydraulic_traits/data/traits_cleaned_3.csv")


##---------------------------------------------------------------
## Analysis
##---------------------------------------------------------------
traits_full <- read.csv("data/traits_cleaned_3.csv")

## take mean for each species within each site's observation -- cant get RE models to converge or make sense.
bal_data <- data.frame(binomial = character(0),
                       P50 = numeric(0),
                       Ks = numeric(0),
                       lat = numeric(0),
                       lon = numeric(0),
                       log_P50 = numeric(0),
                       log_Ks = numeric(0),
                       mean_storm_size = numeric(0),
                       median_storm_size = numeric(0),
                       sd_storm_size = numeric(0),
                       storm_freq = numeric(0),
                       map_reanalysis = numeric(0),
                       nspp = numeric(0),
                       site_name = character(0),
                       mat = numeric(0),
                       a_m = numeric(0))



for (i in 1:nrow(traits_full)) {
  traits_full[i, "latlon"] <- paste0(traits_full[i, "lat"], traits_full[i, "lon"])
}

for (i in unique(traits_full$latlon)) {

  if (length(unique(traits_full[traits_full$latlon == i, "binomial"])) > 3) {
    for (s in unique(traits_full[traits_full$latlon == i, "binomial"])) {

      sdata <- traits_full[traits_full$latlon == i & traits_full$binomial == s,]
      ndata <- data.frame(binomial = s,
                         P50 = mean(sdata$P50),
                         Ks = mean(sdata$Ks),
                         lat = sdata[1,"lat"],
                         lon = sdata[1,"lon"],
                         log_P50 = -log(abs(mean(sdata$P50))),
                         log_Ks = log(mean(sdata$Ks)),
                         mean_storm_size = sdata[1,"mean_storm_size"],
                         median_storm_size = sdata[1,"median_storm_size"],
                         sd_storm_size = sdata[1,"sd_storm_size"],
                         storm_freq = sdata[1,"storm_freq"],
                         map_reanalysis = sdata[1,"map_reanalysis"],
                         mean_temp = sdata[1,"mean_temp"],
                         nspp = length(unique(traits_full[traits_full$latlon == i, "binomial"])),
                         site_name = sdata[1, "site_name"],
                         mat = sdata[1, "mat"],
                         a_m = mean(sdata$a_m))
      bal_data <- rbind(bal_data, ndata)
    }
  }
}

bal_data$map_reanalysis_scaled <- scale(bal_data$map_reanalysis)
bal_data$storm_freq_scaled <- scale(bal_data$storm_freq)
bal_data$mean_storm_size_scaled <- scale(bal_data$mean_storm_size)
bal_data$mean_temp_scaled <- scale(bal_data$mean_temp)

##---------------------------------------------------------------
## Linear models
##---------------------------------------------------------------

fit_min_P50 <- brm(bf(P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
fit_max_P50 <- brm(bf(P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_P50)
summary(fit_max_P50)

## log P50

fit_mean_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50)


fit_min_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50)
plot(fit_min_log_P50)
bal_data$mean_storm_size_scaled <- scale(bal_data$mean_storm_size)
pp_check(fit_min_log_P50)
pp_check(fit_min_log_P50, type = "ecdf_overlay")

fit_min_log_P50_draws <- as_draws_df(fit_min_log_P50)
colnames(fit_min_log_P50_draws)[4] <- "freq_map"
ggplot(data = fit_min_log_P50_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("hydraulic_traits/figures/min_P50_post.pdf")

fit_max_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50)
plot(fit_max_log_P50)
pp_check(fit_max_log_P50)
pp_check(fit_max_log_P50, type = "ecdf_overlay")

fit_max_log_P50_draws <- as_draws_df(fit_max_log_P50)
colnames(fit_max_log_P50_draws)[4] <- "freq_map"
ggplot(data = fit_max_log_P50_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("hydraulic_traits/figures/max_P50_post.pdf")


fit_mean_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50)


fit_min_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50_temp)

fit_max_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50_temp)


## Amax
fit_min_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am)
plot(fit_min_am)
pp_check(fit_min_am)
pp_check(fit_min_am, type = "ecdf_overlay")


fit_min_am_draws <- as_draws_df(fit_min_am)
colnames(fit_min_am_draws)[4] <- "freq_map"
ggplot(data = fit_min_am_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("hydraulic_traits/figures/min_Ks_post.pdf")


fit_max_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am)
plot(fit_max_am)
pp_check(fit_max_am)
pp_check(fit_max_am, type = "ecdf_overlay")

fit_max_am_draws <- as_draws_df(fit_max_am)
colnames(fit_max_am_draws)[4] <- "freq_map"
ggplot(data = fit_max_am_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("hydraulic_traits/figures/max_Ks_post.pdf")

fit_mean_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_am)


fit_min_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am_temp)

fit_max_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am_temp)


##---------------------------------------------------------------
## Log P50 prediction plots
##---------------------------------------------------------------
lev <- quantile(unique(bal_data$map_reanalysis_scaled), c(0.165, 0.495, 0.9))

## Log P50
pdata <- data.frame(storm_freq_scaled = rep(seq(-1.5, 2.7, length.out = 800), 3),
                    map_reanalysis_scaled = rep(c(lev[1], lev[2], lev[3]), each = 800),
                    mean_temp_scaled = rep(0.0, times = 800 * 3))

pred <- posterior_epred(fit_min_log_P50, pdata, dpar = "mu")
pdata$min_mean <- colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "min_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "min_mean_lower"] <- quantile(pred[,i], 0.025)
}

t <- (c(5.0, 12.0) - mean(bal_data$storm_freq)) / sd(bal_data$storm_freq)

a1 <- colMeans(posterior_epred(fit_min_log_P50, data.frame(storm_freq_scaled = t, map_reanalysis_scaled = c(0.0, 0.0)), dpar = "mu"))
a2 <- colMeans(posterior_epred(fit_max_log_P50, data.frame(storm_freq_scaled = t, map_reanalysis_scaled = c(0.0, 0.0)), dpar = "mu"))
a3 <- a2 - a1
a3[1] - a3[2]

a3[1] / a3[2]


pred <- posterior_epred(fit_max_log_P50, pdata, dpar = "mu")
pdata$max_mean <-colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "max_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "max_mean_lower"] <- quantile(pred[,i], 0.025)
}

pred <- posterior_epred(fit_mean_log_P50, pdata, dpar = "mu")
pdata$mean <-colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "mean_lower"] <- quantile(pred[,i], 0.025)
}

pdata$storm_freq_scaled <- (pdata$storm_freq_scaled * sd(bal_data$storm_freq)) + mean(bal_data$storm_freq)

pdata$map_reanalysis_scaled <- as.factor(pdata$map_reanalysis_scaled)






bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.4), "inplot"] <- "Y"

p1 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#253494", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#253494", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#253494", alpha = 0.4) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, y = mean), color = "#253494", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = log_P50, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 15.0)) +
  coord_cartesian(ylim = c(-3, 2.5)) +
  ggtitle(paste0("MAP < ",
                 round(quantile(bal_data$map_reanalysis, 0.33), -2))) +
  ylab("-log(|P50|)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p1

bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.4) &
         bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.75), "inplot"] <- "Y"

p2 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#41b6c4", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#41b6c4", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#41b6c4", alpha = 0.4) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, y = mean), color = "#41b6c4", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = log_P50, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 15.0)) +
  coord_cartesian(ylim = c(-3, 2.5)) +
  ggtitle(paste0(round(quantile(bal_data$map_reanalysis, 0.33), -2),
                 " < MAP < ",
                 round(quantile(bal_data$map_reanalysis, 0.66), -2))) +
  ylab("-log(|P50|)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p2


bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.75), "inplot"] <- "Y"

p3 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#a1dab4", alpha = 0.3) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#a1dab4", alpha = 0.3) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#a1dab4", alpha = 0.3) +
   geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, y = mean), color = "#a1dab4", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = log_P50, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high = "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 20.0)) +
  coord_cartesian(ylim = c(-3, 2.5)) +
  ggtitle(paste0("MAP > ",
                 round(quantile(bal_data$map_reanalysis, 0.66), 0))) +
  ylab("-log(|P50|)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p3

legend <- get_legend(p3)


plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          p3 + theme(legend.position = "none"),
          legend,
          ncol = 4, rel_widths = c(0.3, 0.3, 0.3, 0.1))

ggsave("hydraulic_traits/figures/log_P50_lm.pdf", width = 15, height = 5, units = "in")

##---------------------------------------------------------------
## Log Ks prediction plots
##---------------------------------------------------------------

lev <- quantile(bal_data$map_reanalysis_scaled, c(0.165, 0.495, 0.9))

## Log Ks
pdata <- data.frame(storm_freq_scaled = rep(seq(-2, 2.5, length.out = 200), 3),
                    map_reanalysis_scaled = rep(c(lev[1], lev[2], lev[3]), each = 200))

pred <- posterior_epred(fit_min_log_Ks, pdata, dpar = "mu")
pdata$min_mean <- colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "min_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "min_mean_lower"] <- quantile(pred[,i], 0.025)
}


pred <- posterior_epred(fit_max_log_Ks, pdata, dpar = "mu")
pdata$max_mean <-colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "max_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "max_mean_lower"] <- quantile(pred[,i], 0.025)
}

pdata$storm_freq_scaled <- (pdata$storm_freq_scaled * sd(bal_data$storm_freq)) + mean(bal_data$storm_freq)

pdata$map_reanalysis_scaled <- as.factor(pdata$map_reanalysis_scaled)

bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.33), "inplot"] <- "Y"

p1 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean_upper), fill = "red", alpha = 0.2) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = max_mean_lower, ymax = max_mean_upper), fill = "red", alpha = 0.2) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
            aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
            aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = log_Ks, color = map_reanalysis, alpha = inplot),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", breaks = c(500, 1500, 3000)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(-3.5, 5.5)) +
  ggtitle(paste0("MAP < ", round(quantile(bal_data$map_reanalysis, 0.33), 0))) +
  ylab("log(Ks)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")


bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.33) &
         bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.66), "inplot"] <- "Y"

p2 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean_upper), fill = "red", alpha = 0.2) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = max_mean_lower, ymax = max_mean_upper), fill = "red", alpha = 0.2) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
            aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
            aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = log_Ks, color = map_reanalysis, alpha = inplot),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", breaks = c(500, 1500, 3000)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(-3.5, 5.5)) +
  ggtitle(paste0(round(quantile(bal_data$map_reanalysis, 0.33), 0),
                 " < MAP < ",
                 round(quantile(bal_data$map_reanalysis, 0.66), 0))) +
  ylab("log(Ks)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")


plot(site_data$storm_freq, site_data$storm_freq_fy)
bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.66), "inplot"] <- "Y"

p3 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean_upper), fill = "red", alpha = 0.2) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = max_mean_lower, ymax = max_mean_upper), fill = "red", alpha = 0.2) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
            aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
            aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = log_Ks, color = map_reanalysis, alpha = inplot),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", breaks = c(500, 1500, 3000)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(-3.5, 5.5)) +
  ggtitle(paste0(round(quantile(bal_data$map_reanalysis, 0.66), 0),
                 " < MAP")) +
  ylab("log(Ks)") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(3, "cm"),
        legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))

legend <- get_legend(p3)


plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          p3 + theme(legend.position = "none"),
          legend,
          ncol = 2)

ggsave("hydraulic_traits/figures/log_Ks_lm.pdf")



##---------------------------------------------------------------
## a max prediction plots
##---------------------------------------------------------------


lev <- quantile(bal_data$map_reanalysis_scaled, c(0.2, 0.495, 0.9))

## Log Ks
pdata <- data.frame(storm_freq_scaled = rep(seq(-2, 2.5, length.out = 200), 3),
                    map_reanalysis_scaled = rep(c(lev[1], lev[2], lev[3]), each = 200))

pred <- posterior_epred(fit_min_am, pdata, dpar = "mu")
pdata$min_mean <- colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "min_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "min_mean_lower"] <- quantile(pred[,i], 0.025)
}


pred <- posterior_epred(fit_max_am, pdata, dpar = "mu")
pdata$max_mean <-colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "max_mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "max_mean_lower"] <- quantile(pred[,i], 0.025)
}


pred <- posterior_epred(fit_mean_am, pdata, dpar = "mu")
pdata$mean <-colMeans(pred)
for (i in 1:nrow(pdata)) {
  pdata[i, "mean_upper"] <- quantile(pred[,i], 0.975)
  pdata[i, "mean_lower"] <- quantile(pred[,i], 0.025)
}

t <- (c(5.0, 12.0) - mean(bal_data$storm_freq)) / sd(bal_data$storm_freq)

a1 <- colMeans(posterior_epred(fit_min_am, data.frame(storm_freq_scaled = t, map_reanalysis_scaled = c(0.0, 0.0)), dpar = "mu"))
a2 <- colMeans(posterior_epred(fit_max_am, data.frame(storm_freq_scaled = t, map_reanalysis_scaled = c(0.0, 0.0)), dpar = "mu"))
a3 <- a2 - a1
(a3[1] - a3[2]) / a3[1]
(a3[1] / a3[2])



a1
a2

pdata$storm_freq_scaled <- (pdata$storm_freq_scaled * sd(bal_data$storm_freq)) + mean(bal_data$storm_freq)

pdata$map_reanalysis_scaled <- as.factor(pdata$map_reanalysis_scaled)

bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.4), "inplot"] <- "Y"

p1 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#253494", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#253494", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#253494", alpha = 0.4) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[1],],
              aes(x = storm_freq_scaled, y = mean), color = "#253494", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = a_m, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 13.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle(paste0("MAP < ",
                 round(quantile(bal_data$map_reanalysis, 0.33), -2))) +
  ylab("a max") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p1

bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.4) &
         bal_data$map_reanalysis_scaled < quantile(bal_data$map_reanalysis_scaled, 0.75), "inplot"] <- "Y"

p2 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#41b6c4", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#41b6c4", alpha = 0.4) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#41b6c4", alpha = 0.4) +
  geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, y = mean), color = "#41b6c4", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = a_m, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 13.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle(paste0(round(quantile(bal_data$map_reanalysis, 0.33), -2),
                 " < MAP < ",
                 round(quantile(bal_data$map_reanalysis, 0.66), -2))) +
  ylab("a max") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p2


bal_data$inplot <- "N"
bal_data[bal_data$map_reanalysis_scaled > quantile(bal_data$map_reanalysis_scaled, 0.75), "inplot"] <- "Y"

p3 <- ggplot() +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = min_mean, ymax = max_mean), fill = "#a1dab4", alpha = 0.3) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = min_mean_lower, ymax = min_mean), fill = "#a1dab4", alpha = 0.3) +
  geom_ribbon(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[3],],
              aes(x = storm_freq_scaled, ymin = max_mean, ymax = max_mean_upper), fill = "#a1dab4", alpha = 0.3) +
   geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
              aes(x = storm_freq_scaled, y = mean), color = "#a1dab4", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = max_mean), color = "red", size = 3) +
  #geom_line(data = pdata[pdata$map_reanalysis_scaled == unique(pdata$map_reanalysis_scaled)[2],],
  #          aes(x = storm_freq_scaled, y = min_mean), color = "red", size = 3) +
  geom_point(data = bal_data[bal_data$inplot == "Y",],
             aes(x = storm_freq, y = a_m, color = map_reanalysis),
             size = 6) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high = "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 13.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle(paste0("MAP > ",
                 round(quantile(bal_data$map_reanalysis, 0.66), 0))) +
  ylab("a max") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p3

legend <- get_legend(p3)


plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          p3 + theme(legend.position = "none"),
          legend,
          ncol = 4, rel_widths = c(0.3, 0.3, 0.3, 0.1))

ggsave("hydraulic_traits/figures/am_lm.pdf", width = 15, height = 5, units = "in")




##---------------------------------------------------------------
## Smoothing splines
##---------------------------------------------------------------

fit_min_log_P50_smooth <- brm(bf(log_P50 ~ s(storm_freq_scaled, map_reanalysis_scaled, k = 7), quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50_smooth)
plot(fit_min_log_P50_smooth)
ms1 <- conditional_smooths(fit_min_log_P50_smooth)
plot(ms1, stype = "raster")
pp_check(fit_min_log_P50_smooth)

fit_max_log_P50_smooth <- brm(bf(log_P50 ~ s(storm_freq_scaled, map_reanalysis_scaled, k = 7), quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50_smooth)
plot(fit_max_log_P50_smooth)
ms2 <- conditional_smooths(fit_max_log_P50_smooth)
plot(ms2, stype = "raster")
pp_check(fit_max_log_P50_smooth)

##---------------------------------------------------------------
##Log P50 smooth prediction plots
##---------------------------------------------------------------

storm_freq_scaled <- seq(min(bal_data$storm_freq_scaled)-1, max(bal_data$storm_freq_scaled)+1, length.out = 300)
map_reanalysis_scaled <- seq(min(bal_data$map_reanalysis_scaled)-1, max(bal_data$map_reanalysis_scaled)+1, length.out = 300)

pdata <- data.frame(expand.grid(storm_freq_scaled, map_reanalysis_scaled))

colnames(pdata) <- c("storm_freq_scaled", "map_reanalysis_scaled")
pdata$val <- 1

env <- st_as_sf(bal_data, coords = c("storm_freq_scaled", "map_reanalysis_scaled"))
hull <- st_concave_hull(st_union(st_buffer(env, 0.3)), 0.8)

pdata <- rasterFromXYZ(pdata)
env <- mask(pdata, as_Spatial(hull))
pdata <- as.data.frame(env, xy = TRUE, na.rm = TRUE)[,1:2]
colnames(pdata) <- c("storm_freq_scaled", "map_reanalysis_scaled")

pred <- posterior_epred(fit_min_log_P50_smooth, pdata, dpar = "mu")
pdata$min_mean <- colMeans(pred)

pred <- posterior_epred(fit_max_log_Ks, pdata, dpar = "mu")
pdata$max_mean <-colMeans(pred)

pdata$storm_freq_scaled <- (pdata$storm_freq_scaled * sd(bal_data$storm_freq)) + mean(bal_data$storm_freq)
pdata$map_reanalysis_scaled <- (pdata$map_reanalysis_scaled * sd(bal_data$map_reanalysis)) + mean(bal_data$map_reanalysis)


p1 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = min_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(-2.5, 2, length.out = 10), 2)) +
  labs(fill = "Minimum -log(|P50|)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted minimum -log(|P50|)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")

p2 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = max_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(1, 4.5, length.out = 10), 2)) +
  labs(fill = "Maximum -log(|P50|)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted maximum -log(|P50|)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")

p3 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = max_mean - min_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(1, 6.5, length.out = 10), 2)) +
  labs(fill = "Range -log(|P50|)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted range in -log(|P50|)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")


plot_grid(p1, p2, p3, ncol = 2)
ggsave("hydraulic_traits/figures/spline_P50.pdf")


##---------------------------------------------------------------
## Ks smoothing splines
##---------------------------------------------------------------


fit_min_log_Ks_smooth <- brm(bf(log_Ks ~ s(storm_freq_scaled, map_reanalysis_scaled, k = 7), quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_Ks_smooth)
plot(fit_min_log_Ks_smooth)
ms1 <- conditional_smooths(fit_min_log_Ks_smooth)
plot(ms1, stype = "raster")
pp_check(fit_min_log_Ks_smooth)

fit_max_log_Ks_smooth <- brm(bf(log_Ks ~ s(storm_freq_scaled, map_reanalysis_scaled, k = 7), quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_Ks_smooth)
plot(fit_max_log_Ks_smooth)
ms2 <- conditional_smooths(fit_max_log_Ks_smooth)
plot(ms2, stype = "raster")
pp_check(fit_max_log_Ks_smooth)

##---------------------------------------------------------------
##Log Ks smooth prediction plots
##---------------------------------------------------------------





storm_freq_scaled <- seq(min(bal_data$storm_freq_scaled)-1, max(bal_data$storm_freq_scaled)+1, length.out = 300)
map_reanalysis_scaled <- seq(min(bal_data$map_reanalysis_scaled)-1, max(bal_data$map_reanalysis_scaled)+1, length.out = 300)

pdata <- data.frame(expand.grid(storm_freq_scaled, map_reanalysis_scaled))

colnames(pdata) <- c("storm_freq_scaled", "map_reanalysis_scaled")
pdata$val <- 1

env <- st_as_sf(bal_data, coords = c("storm_freq_scaled", "map_reanalysis_scaled"))
hull <- st_concave_hull(st_union(st_buffer(env, 0.3)), 0.8)

pdata <- rasterFromXYZ(pdata)
env <- mask(pdata, as_Spatial(hull))
pdata <- as.data.frame(env, xy = TRUE, na.rm = TRUE)[,1:2]
colnames(pdata) <- c("storm_freq_scaled", "map_reanalysis_scaled")

pred <- posterior_epred(fit_min_log_Ks_smooth, pdata, dpar = "mu")
pdata$min_mean <- colMeans(pred)

pred <- posterior_epred(fit_max_log_Ks, pdata, dpar = "mu")
pdata$max_mean <-colMeans(pred)

pdata$storm_freq_scaled <- (pdata$storm_freq_scaled * sd(bal_data$storm_freq)) + mean(bal_data$storm_freq)
pdata$map_reanalysis_scaled <- (pdata$map_reanalysis_scaled * sd(bal_data$map_reanalysis)) + mean(bal_data$map_reanalysis)


p1 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = min_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(-3.5, 1, length.out = 10), 2)) +
  labs(fill = "Minimum log(Ks)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted minimum log(Ks)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")

p2 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = max_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(1, 4.5, length.out = 10), 2)) +
  labs(fill = "Maximum log(Ks)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted maximum log(Ks)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")

p3 <- ggplot() +
  geom_contour_filled(data = pdata, aes(x = storm_freq_scaled, y = map_reanalysis_scaled, z = max_mean - min_mean),
                      bins = 20) +
  geom_point(data = bal_data,
             aes(x = storm_freq, y = map_reanalysis),
             size = 6) +
  scale_fill_discretised(breaks = round(seq(1, 6.5, length.out = 10), 2)) +
  labs(fill = "Range log(Ks)") +
  ylab("Mean annual precip. (mm)") +
  xlab("Storm day frequency (n. per year)") +
  ggtitle("Predicted range in log(Ks)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")


plot_grid(p1, p2, p3, ncol = 2)
ggsave("hydraulic_traits/figures/spline_Ks.pdf")











##---------------------------------------------------------------
## Scratch
##---------------------------------------------------------------


P50mod1 <- lqm(log_P50 ~ mean_storm_size + storm_freq_scaled, tau = c(0.05, 0.95),
                data = bal_data[!is.na(bal_data$P50),])
summary(P50mod1)

P50mod1 <- lqm(log_P50 ~ map_reanalysis_scaled * storm_freq_scaled, tau = c(0.05, 0.95),
                data = bal_data[!is.na(bal_data$P50),])
summary(P50mod1)


log_P50_max_form <- list(log_P50 ~ storm_freq_scaled * map_reanalysis_scaled)
log_P50_max_evgam <- evgam(log_P50_max_form, bal_data,
                       family = "ald", ald.args = list(tau = 1 - zeta))
summary(log_P50_max_evgam)

log_P50_max_form <- list(log_P50 ~ storm_freq_scaled * map_reanalysis_scaled)
log_P50_max_evgam <- evgam(log_P50_max_form, bal_data,
                       family = "ald", ald.args = list(tau = 0.05))
summary(log_P50_max_evgam)

## prediction plots
range(bal_data$map_reanalysis_scaled)
range(bal_data$map_reanalysis_scaled)
pdata <- data.frame(storm_freq_scaled = rep(seq(0.1, 1.15, length.out = 200), 3),
                    map_reanalysis_scaled = rep(c(0.5, 1.5, 2.5), each = 200))

pred <- predict(P50mod1, pdata)
pdata$min_mean <- pred[,1]
pdata$max_mean <- pred[,2]

ggplot() +
  geom_line(data = pdata,
            aes(x = storm_freq_scaled, y = max_mean, linetype = as.factor(map_reanalysis_scaled)), color = "red",
            size = 3) +
  geom_line(data = pdata,
            aes(x = storm_freq_scaled, y = min_mean, linetype = as.factor(map_reanalysis_scaled)),
            size = 3) +
  geom_point(data = bal_data, aes(x = storm_freq_scaled, y = log_P50, color = log(map_reanalysis)), size = 4) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Ksmod1 <- lqm(log_Ks ~ mean_storm_size + storm_freq_scaled, tau = c(0.05, 0.95),
                data = bal_data[!is.na(bal_data$log_Ks),])
summary(P50mod1)

Ksmod1 <- lqm(log_Ks ~ map_reanalysis_scaled * storm_freq_scaled, tau = c(0.05, 0.95),
                data = bal_data[!is.na(bal_data$log_Ks),])
summary(P50mod1)

## prediction plots
range(bal_data$map_reanalysis_scaled)
range(bal_data$map_reanalysis_scaled)
pdata <- data.frame(storm_freq_scaled = rep(seq(0.1, 1.15, length.out = 200), 3),
                    map_reanalysis_scaled = rep(c(0.75, 1.5, 2.5), each = 200))

pred <- predict(Ksmod1, pdata)
pdata$min_mean <- pred[,1]
pdata$max_mean <- pred[,2]
pdata$storm_freq_scaled <- pdata$storm_freq_scaled * 100
pdata$map_reanalysis_scaled <- pdata$map_reanalysis_scaled * 1000

pdata$map_reanalysis_scaled <- as.factor(pdata$map_reanalysis_scaled)
ggplot() +
  geom_line(data = pdata,
            aes(x = storm_freq_scaled, y = max_mean, linetype = map_reanalysis_scaled), size = 3, color = "red") +
  geom_line(data = pdata,
            aes(x = storm_freq_scaled, y = min_mean, linetype = map_reanalysis_scaled), size = 3) +
  geom_point(data = bal_data, aes(x = storm_freq, y = log_Ks, color = log(map_reanalysis)), size = 4) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
