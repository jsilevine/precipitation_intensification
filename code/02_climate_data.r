##---------------------------------------------------------------
## Process re-analysis climate data and extract relevant metrics
##
## Author: Jacob Levine; email: jacob.levine@utah.edu
##---------------------------------------------------------------

##---------------------------------------------------------------
## 00. load required libraries and data
##---------------------------------------------------------------

lapply(c("ggplot2", "terra", "sf", "ncdf4", "foreach", "doParallel"), library, character.only = TRUE)

traits_full <- read.csv("data/traits_full.csv")

##---------------------------------------------------------------
## 01. Process reanalysis climate data
##---------------------------------------------------------------

# Create a dataframe containing unique site locations
site_data <- data.frame(latlon = unique(traits_full[!is.na(traits_full$latlon), "latlon"]))

# Loop through each site to extract latitude, longitude, and associated metadata
for (i in 1:nrow(site_data)) {
  cordata <- traits_full[traits_full$latlon == site_data[i, "latlon"], ]  # Filter data for the current site
  cordata <- cordata[!is.na(cordata$lat), ]  # Remove rows with missing latitude values

  # Assign latitude, longitude, site name, ID, and mean annual precipitation (MAP) to site_data
  site_data[i, "lat"] <- cordata$lat[1]
  site_data[i, "lon"] <- cordata$lon[1]
  site_data[i, "site_name"] <- cordata$site_name[1]
  site_data[i, "id"] <- cordata$id[1]
  site_data[i, "map"] <- cordata$map[1]
}

# Load growing season raster files
gs_first <- rast("data/growing_season/gs_first.tif")
gs_last <- rast("data/growing_season/gs_last.tif")

# Extract the first and last growing season days for each site
site_data$gs_first <- terra::extract(gs_first, site_data[, c("lon", "lat")])[,2]
site_data$gs_last <- terra::extract(gs_last, site_data[, c("lon", "lat")])[,2]

# Determine the hemisphere for the growing season
site_data$gs_hem <- "N"  # Default to Northern Hemisphere
site_data[site_data$gs_last < site_data$gs_first & !is.na(site_data$gs_first), "gs_hem"] <- "S"  # Southern Hemisphere
site_data[is.na(site_data$gs_first), "gs_hem"] <- "O"  #  undefined growing season

# Define a function to combine climate data matrices
cfun <- function(x, y) {
  for (i in 1:length(x)) {
    x[[i]] <- cbind(x[[i]], y[[i]])
  }
  return(x)
}


cl <- makeCluster(8, outfile = "")
registerDoParallel(cl)

clim_data <- foreach(f = 1:length(pr_files),
                     .combine = 'cfun',
                     .packages = c("ncdf4")) %dopar% {

  # Initialize matrices to store climate statistics for each site
  mean_storm_size <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  median_storm_size <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  sd_storm_size <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  storm_freq <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  map_reanalysis <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)

  mean_storm_size_fy <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  median_storm_size_fy <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  sd_storm_size_fy <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  storm_freq_fy <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  map_reanalysis_fy <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)

  mean_temp <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  median_temp <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)
  sd_temp <- matrix(rep(0, 10 * nrow(site_data)), ncol = 10)

  # Open reanalysis data for precipitation (pr) and temperature (tas) data
  pr_data <- nc_open(paste0("data/reanalysis/", pr_files[f]))
  tas_data <- nc_open(paste0("data/reanalysis/", tas_files[f]))

  lon <- ncvar_get(pr_data, "lon")
  lat <- ncvar_get(pr_data, "lat")
  time <- ncvar_get(pr_data, "time")

  # Loop through each site
  for (i in 1:nrow(site_data)) {

    print(paste0("starting row: ", i, " of ", nrow(site_data), " \n for file: ", f))

    # Extract site coordinates
    x <- as.numeric(site_data[i, c("lon", "lat")])

    # Find the closest grid point in the dataset
    latindex <- which(abs(lat - x[2]) == min(abs(lat - x[2])))[1]
    lonindex <- which(abs(lon - x[1]) == min(abs(lon - x[1])))[1]

    start <- c(lonindex, latindex, 1)
    count <- c(1, 1, -1) # Extract all time points at this location

    # Extract precipitation and temperature data for the selected location
    pr <- as.numeric(ncvar_get(pr_data, varid = "pr", start = start, count))
    pr <- pr * 24 * 60 * 60  # Convert precipitation from kg/m²/s to mm/day
    tas <- as.numeric(ncvar_get(tas_data, varid = "tas", start = start, count))
    tas <- tas - 273.15  # Convert temperature from Kelvin to Celsius

    # Extract year from filename
    year <- as.numeric(strsplit(pr_files[f], "_")[[1]][2])
    rstart <- 1  # Index tracking the start of each year

    # Looping over a 10-year period
    for (y in 1:10) {

      # Handling leap years
      if (year %in% leap_years) {
        range <- rstart:(rstart + 365) # 366 days
        ed <- 366
      } else {
        range <- rstart:(rstart + 364) # 365 days
        ed <- 365
      }

      # Subsetting based on the growing season in the Northern or Southern Hemisphere
      if (site_data[i, "gs_hem"] == "O") { # No specific growing season
        subrange <- range
      } else if (site_data[i, "gs_hem"] == "N") { # Northern Hemisphere growing season
        subrange <- range[site_data[i, "gs_first"]:site_data[i, "gs_last"]]
      } else if (site_data[i, "gs_hem"] == "S") { # Southern Hemisphere growing season
        subrange <- range[c(site_data[i, "gs_first"]:ed, 1:site_data[i, "gs_last"])]
      }

      # Extract temperature for the growing season and compute statistics
      tas_sub <- tas[subrange]
      mean_temp[i, y] <- mean(tas_sub)
      median_temp[i, y] <- median(tas_sub)
      sd_temp[i, y] <- sd(tas_sub)

      # Extract precipitation for the growing season and identify storms (>2.5 mm)
      pr_sub <- pr[subrange]
      storms <- pr_sub[pr_sub > 2.5]

      # Compute storm statistics for the growing season
      mean_storm_size[i, y] <- mean(storms)
      median_storm_size[i, y] <- median(storms)
      sd_storm_size[i, y] <- sd(storms)
      storm_freq[i, y] <- length(storms) / (length(subrange) / 30)
      map_reanalysis[i, y] <- sum(pr[subrange])

      # Compute storm statistics for the full year
      pr_sub <- pr[range]
      storms <- pr_sub[pr_sub > 2.5]

      mean_storm_size_fy[i, y] <- mean(storms)
      median_storm_size_fy[i, y] <- median(storms)
      sd_storm_size_fy[i, y] <- sd(storms)
      storm_freq_fy[i, y] <- length(storms) / (length(range) / 30)
      map_reanalysis_fy[i, y] <- sum(pr[range])

      # Move to the next year in the dataset
      year <- year + 1
      rstart <- range[length(range)] + 1

    }
  }

  # Return the computed matrices for all sites in the current file
  list(mean_storm_size, median_storm_size, sd_storm_size, storm_freq, map_reanalysis, mean_temp, median_temp, sd_temp,
       mean_storm_size_fy, median_storm_size_fy, sd_storm_size_fy, storm_freq_fy, map_reanalysis_fy)
}

## save processed climate data
saveRDS(clim_data, "data/reanalysis/clim_data_24.rds")

clim_columns <- c("mean_storm_size", "median_storm_size", "sd_storm_size", "storm_freq",
                  "map_reanalysis", "mean_temp", "median_temp", "sd_temp",
                  "mean_storm_size_fy", "median_storm_size_fy", "sd_storm_size_fy",
                  "storm_freq_fy", "map_reanalysis_fy")

# Append climate data to site_data
for (i in seq_along(clim_columns)) {
  site_data[[clim_columns[i]]] <- rowMeans(clim_data[[i]], na.rm = TRUE)
}

write.csv(site_data, "data/site_data.csv")

## merge trait and climate data
traits_full <- merge(traits_full, site_data, by = c("lat", "lon"))

## remove unnecessary columns
traits_full <- traits_full[,c("family", "genus", "species", "binomial", "P50", "P12", "P88",
                              "Ks", "location", "country", "lat", "lon",
                              "mat", "biome", "site_name.x", "map.x", "map_jan", "map_feb", "map_mar", "map_apr",
                              "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru",
                              "sd_map", "log_P50", "log_Ks", "mean_storm_size", "median_storm_size", "sd_storm_size", "storm_freq",
                              "map_reanalysis", "a_m", "mean_temp")]

## clean up column names
colnames(traits_full) <- c("family", "genus", "species", "binomial", "P50", "P12", "P88",
                              "Ks", "location", "country", "lat", "lon",
                              "mat", "biome", "site_name", "map", "map_jan", "map_feb", "map_mar", "map_apr",
                              "map_may", "map_jun", "map_jul", "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru",
                              "sd_map", "log_P50", "log_Ks", "mean_storm_size", "median_storm_size", "sd_storm_size", "storm_freq",
                              "map_reanalysis", "a_m", "mean_temp")

write.csv(traits_full, "data/traits_and_climate.csv")
