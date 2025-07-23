library(brms)
library(ggplot2)
library(ggdist)

##---------------------------------------------------------------
## More stringent diversity requirements
##---------------------------------------------------------------

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
## P50
##---------------------------------------------------------------

fit_mean_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50)
saveRDS(fit_mean_log_P50, "data/model_objects/fit_mean_log_P50_sens.rds")

fit_min_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50)
saveRDS(fit_min_log_P50, "data/model_objects/fit_min_log_P50_sens.rds")

fit_max_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50)
saveRDS(fit_max_log_P50, "data/model_objects/fit_max_log_P50_sens.rds")



fit_mean_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50_temp)
saveRDS(fit_mean_log_P50_temp, "data/model_objects/fit_mean_log_P50_temp_sens.rds")

fit_min_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50_temp)
saveRDS(fit_min_log_P50_temp, "data/model_objects/fit_min_log_P50_temp_sens.rds")

fit_max_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50_temp)
saveRDS(fit_max_log_P50_temp, "data/model_objects/fit_max_log_P50_temp_sens.rds")



##---------------------------------------------------------------
## Amax
##---------------------------------------------------------------

fit_min_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am)
saveRDS(fit_min_am, "data/model_objects/fit_min_am_sens.rds")


fit_max_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am)
saveRDS(fit_max_am, "data/model_objects/fit_max_am_sens.rds")

fit_mean_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_am)
saveRDS(fit_mean_am, "data/model_objects/fit_mean_am_sens.rds")

fit_mean_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_am_temp)
saveRDS(fit_mean_am_temp, "data/model_objects/fit_mean_am_temp_sens.rds")

fit_min_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am_temp)
saveRDS(fit_min_am_temp, "data/model_objects/fit_min_am_temp_sens.rds")

fit_max_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am_temp)
saveRDS(fit_max_am_temp, "data/model_objects/fit_max_am_temp_sens.rds")
