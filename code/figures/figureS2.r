library(ggplot2)
library(raster)
library(cowplot)
library(brms)
library(ggdist)

traits_full <- read.csv("data/traits_and_climate.csv")

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

bal_data$map_reanalysis_scaled <- scale(bal_data$map_reanalysis)
bal_data$storm_freq_scaled <- scale(bal_data$storm_freq)
bal_data$mean_storm_size_scaled <- scale(bal_data$mean_storm_size)
bal_data$mean_temp_scaled <- scale(bal_data$mean_temp)


##---------------------------------------------------------------
## Figure 2A-C
##---------------------------------------------------------------

fit_min_log_P50 <- readRDS("data/model_objects/fit_min_log_P50_temp.rds")
fit_max_log_P50 <- readRDS("data/model_objects/fit_max_log_P50_temp.rds")
fit_mean_log_P50 <- readRDS("data/model_objects/fit_mean_log_P50_temp.rds")

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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 15.0)) +
  coord_cartesian(ylim = c(-2.5, 2)) +
  ggtitle("MAP < 500") +
  ylab("-log(|P50|)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 15.0)) +
  coord_cartesian(ylim = c(-2.5, 2)) +
  ggtitle("500 < MAP < 1500") +
  xlab("Storm day frequency (n. per year)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high = "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 20.0)) +
  coord_cartesian(ylim = c(-2.5, 2)) +
  ggtitle("MAP > 1500") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p3

legend <- get_legend(p3)

pg1 <- plot_grid(p1 + theme(legend.position = "none"),
                 p2 + theme(legend.position = "none"),
                 p3 + theme(legend.position = "none"),
                 nrow = 1,
                 align = "hv")
plot_grid(pg1, legend, nrow = 1, rel_widths = c(0.9, 0.1))

ggsave("figures/figureS2/figureS2A.pdf", width = 15, height = 4.5, units = "in")

##---------------------------------------------------------------
## a max prediction plots
##---------------------------------------------------------------

fit_min_am <- readRDS("data/model_objects/fit_min_am_temp.rds")
fit_max_am <- readRDS("data/model_objects/fit_max_am_temp.rds")
fit_mean_am <- readRDS("data/model_objects/fit_mean_am_temp.rds")

lev <- quantile(bal_data$map_reanalysis_scaled, c(0.2, 0.495, 0.9))

## Log Ks
pdata <- data.frame(storm_freq_scaled = rep(seq(-2, 2.5, length.out = 200), 3),
                    map_reanalysis_scaled = rep(c(lev[1], lev[2], lev[3]), each = 200),
                    mean_temp_scaled = rep(0.0, times = 600))

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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 10.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle("MAP < 500") +
  ylab("Amax") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high =  "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 10.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle("500 < MAP < 1500") +
  xlab("Storm freq. (n. per month)") +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
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
             size = 4) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_color_gradient(trans = "log", limits = c(100,3500), breaks = c(500, 1500, 3000), low = "#252525", high = "#bdbdbd") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(1.0, 20.0)) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle(paste0("MAP > 1500")) +
  labs(color = "MAP (mm)") +
  guides(color = "colorbar", alpha = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.position = "right")
p3

legend <- get_legend(p3)


pg1 <- plot_grid(p1 + theme(legend.position = "none"),
                 p2 + theme(legend.position = "none"),
                 p3 + theme(legend.position = "none"),
                 nrow = 1,
                 align = "hv")
plot_grid(pg1, legend, nrow = 1, rel_widths = c(0.9, 0.1))

ggsave("figures/figureS2/figureS2B.pdf", width = 15, height = 4.5, units = "in")
