##---------------------------------------------------------------
## Process re-analysis climate data and extract relevant metrics
##
## Author: Jacob Levine; email: jacob.levine@utah.edu
##---------------------------------------------------------------


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
library(doParallel)

traits_full <- read.csv("data/traits_full.csv")

##---------------------------------------------------------------
## Reanalysis processing
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

# Load the world map shapefile for visualization
world <- st_read("data/World_Countries_Generalized.shp")

# Transform the coordinate reference system (CRS) to a Robinson projection
world <- st_transform(world, crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

# Convert site data to a spatial object (Simple Features format) using latitude and longitude
loc <- st_as_sf(site_data[!is.na(site_data$lon), ], coords = c("lon", "lat"), crs = 4326)

# Plot the world map and overlay site locations with color representing mean annual precipitation (MAP)
ggplot(world) +
  geom_sf(fill = "gray", color = "white") +  # Base map
  geom_sf(data = loc, aes(color = map), size = 5) +  # Add site points
  theme_bw() +  # Apply a clean theme
  labs(color = "Mean annual \n precipitation") +  # Legend label
  scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +  # Color gradient for MAP
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on x-axis
  scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on y-axis
  ggtitle("Xylem trait data origins")  # Title of the plot

# Save the plot as a PDF
ggsave("figures/trait_locations.pdf")

# Initialize matrices to store climate statistics for 110 years across all sites
mean_storm_size <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
median_storm_size <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
sd_storm_size <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
storm_freq <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
map_reanalysis <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
mean_temp <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
median_temp <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)
sd_temp <- matrix(rep(0, 110 * nrow(site_data)), ncol = 110)

# Define leap years from 1904 to 2008 for time series calculations
leap_years <- seq(1904, 2008, by = 4)

# Retrieve a list of climate data files
clim_files <- list.files("data/reanalysis/")

# Filter the list to separate precipitation (pr) and temperature (tas) files
pr_files <- clim_files[grepl("pr", clim_files)]
tas_files <- clim_files[grepl("tas", clim_files)]

# Load growing season raster files
gs_first <- raster("data/gs_first.tif")
gs_last <- raster("data/gs_last.tif")

# Extract the first and last growing season days for each site
site_data$gs_first <- extract(gs_first, site_data[, c("lon", "lat")])
site_data$gs_last <- extract(gs_last, site_data[, c("lon", "lat")])

# Determine the hemisphere for the growing season
site_data$gs_hem <- "N"  # Default to Northern Hemisphere
site_data[site_data$gs_last < site_data$gs_first & !is.na(site_data$gs_first), "gs_hem"] <- "S"  # Southern Hemisphere
site_data[is.na(site_data$gs_first), "gs_hem"] <- "O"  # Outlier or undefined growing season

# Define a function to combine climate data matrices
cfun <- function(x, y) {
  for (i in 1:length(x)) {
    x[[i]] <- cbind(x[[i]], y[[i]])  # Combine corresponding elements in the lists
  }
  return(x)
}

# Set up parallel computing with 8 cores
cl <- makeCluster(8, outfile = "")
registerDoParallel(cl)


# Parallel processing loop using foreach and dopar for climate data processing
clim_data <- foreach(f = 1:length(pr_files),
                     .combine = 'cfun', # Combining results using a custom function
                     .packages = c("ncdf4")) %dopar% { # Running in parallel with ncdf4 package

  # Initializing matrices to store various climate statistics for each site
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

  # Opening NetCDF files for precipitation (pr) and temperature (tas) data
  pr_data <- nc_open(paste0("data/reanalysis/", pr_files[f]))
  tas_data <- nc_open(paste0("data/reanalysis/", tas_files[f]))

  # Extracting spatial and temporal variables from NetCDF files
  lon <- ncvar_get(pr_data, "lon")
  lat <- ncvar_get(pr_data, "lat")
  time <- ncvar_get(pr_data, "time")

  # Loop through each site in site_data
  for (i in 1:nrow(site_data)) {

    print(paste0("starting row: ", i, " of ", nrow(site_data), " \n for file: ", f))

    # Extracting site coordinates
    x <- as.numeric(site_data[i, c("lon", "lat")])

    # Finding the closest grid point in the dataset
    latindex <- which(abs(lat - x[2]) == min(abs(lat - x[2])))[1]
    lonindex <- which(abs(lon - x[1]) == min(abs(lon - x[1])))[1]

    # Defining extraction parameters for NetCDF variables
    start <- c(lonindex, latindex, 1) # Start at the selected grid point
    count <- c(1, 1, -1) # Extract all time points at this location

    # Extract precipitation and temperature data for the selected location
    pr <- as.numeric(ncvar_get(pr_data, varid = "pr", start = start, count))
    pr <- pr * 24 * 60 * 60  # Convert precipitation from kg/m²/s to mm/day
    tas <- as.numeric(ncvar_get(tas_data, varid = "tas", start = start, count))
    tas <- tas - 273.15  # Convert temperature from Kelvin to Celsius

    # Extract year from the filename
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

      print(length(storms) / (length(subrange) / 30)) # Print storm frequency

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

saveRDS(clim_data, "data/reanalysis/clim_data_24.rds")
clim_data <- readRDS("data/reanalysis/clim_data_24.rds")


## append climate data to site data
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

write.csv(site_data, "data/site_data.csv")

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

write.csv(traits_full, "data/traits_and_climate.csv")
