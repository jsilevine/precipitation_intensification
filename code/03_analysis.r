
library(brms)
library(ggplot2)
library(ggdist)

##---------------------------------------------------------------
## Analysis
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
## P50
##---------------------------------------------------------------

fit_mean_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50)
saveRDS(fit_mean_log_P50, "data/model_objects/fit_mean_log_P50.rds")

fit_min_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50)
saveRDS(fit_min_log_P50, "data/model_objects/fit_min_log_P50.rds")

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
ggsave("figures/model_results/min_P50_post.pdf")

fit_max_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50)
saveRDS(fit_max_log_P50, "data/model_objects/fit_max_log_P50.rds")

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
ggsave("figures/model_results/max_P50_post.pdf")


fit_mean_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_log_P50_temp)
saveRDS(fit_mean_log_P50_temp, "data/model_objects/fit_mean_log_P50_temp.rds")

fit_min_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50_temp)
saveRDS(fit_min_log_P50_temp, "data/model_objects/fit_min_log_P50_temp.rds")

fit_max_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50_temp)
saveRDS(fit_max_log_P50_temp, "data/model_objects/fit_max_log_P50_temp.rds")



##---------------------------------------------------------------
## Amax
##---------------------------------------------------------------

fit_min_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am)
saveRDS(fit_min_am, "data/model_objects/fit_min_am.rds")


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
ggsave("figures/model_resultsmin_am_post.pdf")


fit_max_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am)
saveRDS(fit_max_am, "data/model_objects/fit_max_am.rds")

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
ggsave("figures/model_resultsmax_am_post.pdf")

fit_mean_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_am)
saveRDS(fit_mean_am, "data/model_objects/fit_mean_am.rds")

fit_mean_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_mean_am_temp)
saveRDS(fit_mean_am_temp, "data/model_objects/fit_mean_am_temp.rds")

fit_min_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am_temp)
saveRDS(fit_min_am_temp, "data/model_objects/fit_min_am_temp.rds")

fit_max_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am_temp)
saveRDS(fit_max_am_temp, "data/model_objects/fit_max_am_temp.rds")


