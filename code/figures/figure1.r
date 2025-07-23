
##---------------------------------------------------------------
## Figure 1 - Map of functional trait data collection sites
##
## Author: Jacob I. Levine; email: jacob.levine@utah.edu
##---------------------------------------------------------------

library(sf)
library(ggplot2)
library(cowplot)

site_data <- read.csv("data/site_data.csv")
traits_full <- read.csv("data/traits_and_climate.csv")

unique(traits_full$latlon)

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

  if (length(unique(traits_full[traits_full$latlon == i, "binomial"])) > 1) {
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

for (i in 1:nrow(bal_data)) {
  bal_data[i, "latlon"] <- paste0(bal_data[i, "lat"], bal_data[i, "lon"])
}

bal_data <- bal_data[!is.na(bal_data$P50) | !is.na(bal_data$a_m),]

mapdata <- unique(bal_data[,c("lat", "lon")])


for (i in 1:nrow(mapdata)) {
  sd <- bal_data[bal_data$lat == mapdata[i,"lat"] & bal_data$lon == mapdata[i, "lon"],]
  mapdata[i,"P50_range"] <- max(sd$log_P50) - min(sd$log_P50)
  mapdata[i,"Amax_range"] <- max(sd$a_m) - min(sd$a_m)
  mapdata[i,"map_reanalysis"] <- sd[1, "map_reanalysis"]
  mapdata[i,"storm_freq"] <- sd[1, "storm_freq"]
}
nrow(mapdata)

sf_mapdata <- st_as_sf(mapdata, coords = c("lon", "lat"), crs = 4326)
#sf_mapdata <- sf_mapdata[!is.na(sf_mapdata$P50_range) | !is.na(sf_mapdata$Amax_range),]

world <- st_read("data/World_Countries_Generalized.shp")
world <- st_transform(world, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

p1 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_mapdata, aes(color = map_reanalysis), size = 5, alpha = 0.9) +
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
p1

p2 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_mapdata, aes(color = storm_freq), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Storm frequency \n (n. per year)") +
  scale_color_gradient2(low = "#bfd3e6", mid = "#8c96c6", high = "#810f7c", midpoint = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Storm frequency (n. per year)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")
p2

p3 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_mapdata[!is.na(sf_mapdata$P50_range),], aes(color = P50_range), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Range in \n log P50") +
  scale_color_gradient2(low = "#fcc5c0", mid = "#f768a1", high = "#7a0177", midpoint = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Range in log P50 (max - min)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")
p3


p4 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_mapdata[!is.na(sf_mapdata$Amax_range),], aes(color = Amax_range), size = 5, alpha = 0.9) +
  theme_bw() +
  labs(color = "Range in \n Amax") +
  scale_color_gradient2(low = "#ccece6", mid = "#66c2a4", high = "#006d2c", midpoint = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Range in Amax (max - min)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "right")
p4


plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv")
ggsave("figures/figure1/figure1.pdf", width = 25, height = 15)
