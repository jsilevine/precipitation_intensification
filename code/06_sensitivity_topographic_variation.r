
##---------------------------------------------------------------
## Understanding whether topographic/edaphic complexity varies
## with precipitation frequency and magnitude
##
## Author: jacob levine; email: jacob.levine@utah.edu
##---------------------------------------------------------------

library(ggplot2)
library(elevatr)
library(raster)
library(terra)
library(sf)
library(tidyterra)
library(cowplot)

traits_full <- read.csv("data/traits_and_climate.csv")

unq_loc <- unique(traits_full[, c("lat", "lon")])

spatial <- st_as_sf(unq_loc, coords = c("lon", "lat"), crs = 4326)

## for 3km buffer
circles <- st_buffer(spatial, 3000)

unq_loc$elev_var <- NA

for (i in 1:nrow(circles)) {
  print(paste("running", i, "of", nrow(circles)))
  dem <- rast(get_elev_raster(circles[i,], z = 14, clip = "bbox"))
  dem <- mask(dem, circles[i,])
  dem.df <- as.data.frame(dem)
  unq_loc[i, "elev_var"] <- sd(dem.df[,1])
}

unq_loc <- merge(unq_loc,
                     traits_full[,c("lat", "lon", "map_reanalysis", "storm_freq", "mean_temp")],
                     by = c("lat", "lon"), all.x = FALSE)
unq_loc <- unq_loc[!duplicated(unq_loc),]

write.csv(unq_loc, "data/location_topo_3km.csv")


mod1 <- glm(elev_var ~ storm_freq * map_reanalysis, data = unq_loc, family = Gamma(link = "log"))
summary(mod1)

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(250, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

a <- ggplot(data = unq_loc[unq_loc$map_reanalysis < 400,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 400") +
  theme_bw()
a


pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(700, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

b <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 400 & unq_loc$map_reanalysis < 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("400 < MAP < 1000") +
  theme_bw()

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(1500, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

c <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 22), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 1000") +
  theme_bw()

plot_grid(a, b, c, nrow = 1)

ggsave("figures/figureS1a.png", width = 12, height = 4)


##---------------------------------------------------------------
## 1km
##---------------------------------------------------------------


unq_loc <- unique(traits_full[, c("lat", "lon")])

spatial <- st_as_sf(unq_loc, coords = c("lon", "lat"), crs = 4326)

circles <- st_buffer(spatial, 1000)

for (i in 1:nrow(circles)) {
  print(paste("running", i, "of", nrow(circles)))
  dem <- rast(get_elev_raster(circles[i,], z = 14, clip = "bbox"))
  dem <- mask(dem, circles[i,])
  dem.df <- as.data.frame(dem)
  unq_loc[i, "elev_var"] <- sd(dem.df[,1])
}

unq_loc <- merge(unq_loc,
                     traits_full[,c("lat", "lon", "map_reanalysis", "storm_freq", "mean_temp")],
                     by = c("lat", "lon"))
unq_loc <- unq_loc[!duplicated(unq_loc),]

write.csv(unq_loc, "data/location_topo_1km.csv")

mod1 <- glm(elev_var ~ storm_freq * map_reanalysis, data = unq_loc, family = Gamma(link = "log"))
summary(mod1)

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(250, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

a <- ggplot(data = unq_loc[unq_loc$map_reanalysis < 400,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 400") +
  theme_bw()
a


pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(700, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

b <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 400 & unq_loc$map_reanalysis < 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("400 < MAP < 1000") +
  theme_bw()

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(1500, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

c <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 22), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 1000") +
  theme_bw()

plot_grid(a, b, c, nrow = 1)
ggsave("figures/figureS1b.png", width = 12, height = 4)

##---------------------------------------------------------------
## 5km
##---------------------------------------------------------------

unq_loc <- unique(traits_full[, c("lat", "lon")])

spatial <- st_as_sf(unq_loc, coords = c("lon", "lat"), crs = 4326)

circles <- st_buffer(spatial, 5000)

for (i in 1:nrow(circles)) {
  print(paste("running", i, "of", nrow(circles)))
  dem <- rast(get_elev_raster(circles[i,], z = 14, clip = "bbox"))
  dem <- mask(dem, circles[i,])
  dem.df <- as.data.frame(dem)
  unq_loc[i, "elev_var"] <- sd(dem.df[,1])
}

unq_loc <- merge(unq_loc,
                     traits_full[,c("lat", "lon", "map_reanalysis", "storm_freq", "mean_temp")],
                     by = c("lat", "lon"))
unq_loc <- unq_loc[!duplicated(unq_loc),]

write.csv(unq_loc, "data/location_topo_5km.csv")
unq_loc <- read.csv("data/location_topo_5km.csv", row.names = 1)
head(unq_loc)

mod1 <- glm(elev_var ~ storm_freq * map_reanalysis, data = unq_loc, family = Gamma(link = "log"))
summary(mod1)

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(250, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

a <- ggplot(data = unq_loc[unq_loc$map_reanalysis < 400,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 400") +
  theme_bw()
a


pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(700, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

b <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 400 & unq_loc$map_reanalysis < 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("400 < MAP < 1000") +
  theme_bw()

pdata <- data.frame(storm_freq = seq(0, 25, length.out = 26),
                    map_reanalysis = rep(1500, times = 26))

pred <- predict(mod1, pdata, se.fit = TRUE, type = "response")
pdata$mean_prediction <- pred$fit
pdata$lower_prediction <- pred$fit - 1.97 * pred$se.fit
pdata$upper_prediction <- pred$fit + 1.97 * pred$se.fit

c <- ggplot(data = unq_loc[unq_loc$map_reanalysis > 1000,]) +
  geom_point(aes(x = storm_freq, y = elev_var), size = 2.5) +
  geom_ribbon(data = pdata, aes(x = storm_freq, ymin = lower_prediction, ymax = upper_prediction),
              fill = "lightblue", alpha = 0.6) +
  geom_line(data = pdata, aes(x = storm_freq, y = mean_prediction), size = 1.5, color = "black") +
  scale_x_continuous(limits = c(0, 22), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  ylab("std. dev. elevation (m)") +
  xlab("storm frequency") +
  ggtitle("MAP < 1000") +
  theme_bw()

plot_grid(a, b, c, nrow = 1)
ggsave("figures/figureS1c.png", width = 12, height = 4)
