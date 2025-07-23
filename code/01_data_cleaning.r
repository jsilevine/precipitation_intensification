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

xylem <- read.csv("data/xylem_traits.csv")

# Filter the dataset to only include adult plant samples of stems (S) or petioles (P) from natural populations
traits_full <- xylem[xylem[,"Natural...Greenhouse"] == "N" &
                    xylem[,"Plant.organ"] %in% c("S", "P") &
                    xylem[,"Developmental.stage"] == "A", ]

# Remove untrusted methods for determining P50
traits_full[!(traits_full$P50.method %in% c("DH", "CE", "CA", "AD", "AS", "Pn")), "P50..MPa."] <- NA

## Clean and standardize latitude and longitude data
for (i in 1:nrow(traits_full)) {

  # Convert Latitude values if LATnew is missing
  if (is.na(traits_full[i, "LATnew"]) & !is.na(traits_full[i, "Latitude"])) {

    # Extract and clean latitude values
    nlat <- traits_full[i, "Latitude"]
    nlat <- gsub("'", "", nlat)  # Remove apostrophes
    if (grepl("S", nlat)) nlat <- paste0("-", nlat)  # Mark southern latitudes as negative
    nlat <- gsub("N|S", "", nlat)  # Remove N and S indicators

    # Standardize format for conversion
    if (nchar(nlat) == 5) {
      nlat <- paste0(nlat, " 00")
    } else if (!grepl("\"", nlat)) {
      nlat <- paste0(nlat, "00")
    } else {
      nlat <- gsub("\"", "", nlat)  # Remove quotation marks if present
    }

    # Convert to decimal degrees
    traits_full[i, "LATnew"] <- as.numeric(conv_unit(nlat, from = "deg_min_sec", to = "dec_deg"))

    # Extract and clean longitude values
    nlon <- traits_full[i, "Longitude"]
    nlon <- gsub("'", "", nlon)
    if (grepl("W", nlon)) nlon <- paste0("-", nlon)  # Mark western longitudes as negative
    nlon <- gsub("W|E", "", nlon)

    # Standardize format for conversion
    if (nchar(nlon) == 5) {
      nlon <- paste0(nlon, " 00")
    } else if (!grepl("\"", nlon)) {
      nlon <- paste0(nlon, "00")
    } else {
      nlon <- gsub("\"", "", nlon)
    }

    # Convert to decimal degrees
    traits_full[i, "LONnew"] <- as.numeric(conv_unit(nlon, from = "deg_min_sec", to = "dec_deg"))
  }
}

tokeep <- c("Cleaned.family", "Cleaned.genus", "Cleaned.species", "Cleaned.binomial",
            "P50..MPa.", "P12..MPa.", "P88..MPa.", "Ks..kg.m.1.MPa.1.s.1.",
            "Location", "Country", "LATnew", "LONnew", "MAP..mean.annual.precipitation..mm.",
            "MAT..mean.annual.temperature...C..", "Biome", "site.name", "PPTbest", "PPTjan",
            "PPTfeb", "PPTmar", "PPTapr", "PPTmay", "PPTjun", "PPTjul", "PPTaug", "PPTsep",
            "PPToct", "PPTnov", "PPTdec", "PPTcru", "NewID", "Amax..micromol.m2.s.")
traits_full <- traits_full[, tokeep]

# Rename columns for consistency and readability
newnames <- c("family", "genus", "species", "binomial", "P50", "P12", "P88",
              "Ks", "location", "country", "lat", "lon", "map_old",
              "mat", "biome", "site_name", "map", "map_jan", "map_feb", "map_mar",
              "map_apr", "map_may", "map_jun", "map_jul", "map_aug", "map_sep",
              "map_oct", "map_nov", "map_dec", "map_cru", "id", "a_m")
colnames(traits_full) <- newnames

# Fill missing "MAP" values with "MAP_old" values if available
for (i in 1:nrow(traits_full)) {
  if (is.na(traits_full[i, "map"]) && !is.na(traits_full[i, "map_old"])) {
    traits_full[i, "map"] <- traits_full[i, "map_old"]
  }
}

# Remove the redundant "map_old" column
traits_full <- traits_full[, colnames(traits_full) != "map_old"]

# Convert specific columns to numeric format
tonum <- c("P50", "P12", "P88", "Ks", "lat", "lon", "mat", "map", "map_jan",
           "map_feb", "map_mar", "map_apr", "map_may", "map_jun", "map_jul",
           "map_aug", "map_sep", "map_oct", "map_nov", "map_dec", "map_cru", "a_m")

for (c in tonum) {
  traits_full[,c] <- as.numeric(traits_full[,c])
}

# Compute standard deviation of monthly precipitation
traits_full$sd_map <- NA
for (i in 1:nrow(traits_full)) {
  traits_full[i, "sd_map"] <- sd(traits_full[i, 17:28], na.rm = TRUE)
}
# Normalize standard deviation by mean annual precipitation
traits_full$sd_map <- traits_full$sd_map / traits_full$map

# Compute logarithmic transformations for key traits
traits_full$log_P50 <- -log(abs(traits_full$P50))
traits_full$log_Ks <- log(traits_full$Ks)

# Filter out rows where all key traits (P50, Ks, a_m) are missing
traits_full <- traits_full[!is.na(traits_full$P50) | !is.na(traits_full$Ks) | !is.na(traits_full$a_m),]

# Generate a combined lat-lon identifier for each entry
traits_full$latlon <- NA
for (i in 1:nrow(traits_full)) {
  print(i)  # Print iteration number for tracking progress
  if (!is.na(traits_full[i, "lat"])) {
    traits_full[i, "latlon"] <- paste0(traits_full[i, "lat"], traits_full[i, "lon"])
  }
}

# Save the cleaned and processed dataset
write.csv(traits_full, "data/traits_full.csv")


nrow(traits_full[is.na(traits_full$P50) & is.na(traits_full$a_m),])
nrow(traits_full)
