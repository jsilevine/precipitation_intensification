##---------------------------------------------------------------
## 03_data_analysis.r -- Global Hydraulic Trait Analysis
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

##---------------------------------------------------------------
## 00. Load required libraries and data
##---------------------------------------------------------------

invisible(lapply(c("brms", "ggplot2", "ggdist"), library, character.only = TRUE))

traits_full <- read.csv("data/traits_and_climate.csv")

##---------------------------------------------------------------
## 01. Prepare data for analysis
##---------------------------------------------------------------

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

## unique latlon identifier for each site
traits_full$latlon <- paste0(traits_full$lat, traits_full$lon)

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

## Scale selected variables for modeling
scaled_vars <- c("map_reanalysis", "storm_freq", "mean_storm_size", "mean_temp")
for (var in scaled_vars) {
  bal_data[[paste0(var, "_scaled")]] <- scale(bal_data[[var]])
}
write.csv(bal_data, "data/bal_data.csv", row.names = FALSE)

##---------------------------------------------------------------
## 02. Fit models for log P50
##---------------------------------------------------------------

## fit model for median log P50
fit_median_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_median_log_P50)
saveRDS(fit_median_log_P50, "data/model_objects/fit_median_log_P50.rds")

## Fit model for minimum log P50
fit_min_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
                       data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50)
saveRDS(fit_min_log_P50, "data/model_objects/fit_min_log_P50.rds")

## Plot and perform posterior predictive checks for minimum log P50
plot(fit_min_log_P50)
pp_check(fit_min_log_P50)
pp_check(fit_min_log_P50, type = "ecdf_overlay")

## Visualize posterior distributions for minimum log P50
fit_min_log_P50_draws <- as_draws_df(fit_min_log_P50)
colnames(fit_min_log_P50_draws)[4] <- "freq_map"
ggplot(data = fit_min_log_P50_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("figures/model_results/min_P50_post.pdf")

## Fit model for maximum log P50
fit_max_log_P50 <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
                       data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50)
saveRDS(fit_max_log_P50, "data/model_objects/fit_max_log_P50.rds")

## Plot and perform posterior predictive checks for maximum log P50
plot(fit_max_log_P50)
pp_check(fit_max_log_P50)
pp_check(fit_max_log_P50, type = "ecdf_overlay")

## Visualize posterior distributions for maximum log P50
fit_max_log_P50_draws <- as_draws_df(fit_max_log_P50)
colnames(fit_max_log_P50_draws)[4] <- "freq_map"
ggplot(data = fit_max_log_P50_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("figures/model_results/max_P50_post.pdf")

##---------------------------------------------------------------
## 02.5. Fit models for log P50 with temperature
##---------------------------------------------------------------

fit_median_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_median_log_P50_temp)
saveRDS(fit_median_log_P50_temp, "data/model_objects/fit_median_log_P50_temp.rds")

fit_min_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_log_P50_temp)
saveRDS(fit_min_log_P50_temp, "data/model_objects/fit_min_log_P50_temp.rds")

fit_max_log_P50_temp <- brm(bf(log_P50 ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_log_P50_temp)
saveRDS(fit_max_log_P50_temp, "data/model_objects/fit_max_log_P50_temp.rds")

##---------------------------------------------------------------
## 03. Fit models for Amax
##---------------------------------------------------------------

## fit model for median Amax
fit_median_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_median_am)
saveRDS(fit_median_am, "data/model_objects/fit_median_am.rds")

## Fit model for minimum Amax
fit_min_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.05),
                  data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am)
saveRDS(fit_min_am, "data/model_objects/fit_min_am.rds")

## Plot and perform posterior predictive checks for minimum Amax
plot(fit_min_am)
pp_check(fit_min_am)
pp_check(fit_min_am, type = "ecdf_overlay")

## Visualize posterior distributions for minimum Amax
fit_min_am_draws <- as_draws_df(fit_min_am)
colnames(fit_min_am_draws)[4] <- "freq_map"
ggplot(data = fit_min_am_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("figures/model_resultsmin_am_post.pdf")

## Fit model for maximum Amax
fit_max_am <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled, quantile = 0.95),
                  data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am)
saveRDS(fit_max_am, "data/model_objects/fit_max_am.rds")

## Plot and perform posterior predictive checks for maximum Amax
plot(fit_max_am)
pp_check(fit_max_am)
pp_check(fit_max_am, type = "ecdf_overlay")

## Visualize posterior distributions for maximum Amax
fit_max_am_draws <- as_draws_df(fit_max_am)
colnames(fit_max_am_draws)[4] <- "freq_map"
ggplot(data = fit_max_am_draws) +
  stat_halfeye(aes(x = b_Intercept, y = "Intercept")) +
  stat_halfeye(aes(x = b_storm_freq_scaled, y = "Storm frequency")) +
  stat_halfeye(aes(x = b_map_reanalysis_scaled, y = "MAP")) +
  stat_halfeye(aes(x = freq_map, y = "Storm frequency : MAP")) +
  theme_bw()
ggsave("figures/model_resultsmax_am_post.pdf")

##---------------------------------------------------------------
## 03.5. Fit models for Amax with temperature
##---------------------------------------------------------------

fit_median_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.5),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_median_am_temp)
saveRDS(fit_median_am_temp, "data/model_objects/fit_median_am_temp.rds")

fit_min_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.05),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_min_am_temp)
saveRDS(fit_min_am_temp, "data/model_objects/fit_min_am_temp.rds")

fit_max_am_temp <- brm(bf(a_m ~ storm_freq_scaled*map_reanalysis_scaled + mean_temp_scaled, quantile = 0.95),
           data = bal_data, family = asym_laplace(), cores = 4)
summary(fit_max_am_temp)
saveRDS(fit_max_am_temp, "data/model_objects/fit_max_am_temp.rds")


