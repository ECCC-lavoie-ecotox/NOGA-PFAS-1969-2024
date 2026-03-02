###############################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 2 - Suess Effect Correction using the SuessR R Package
### Script originally written by Rose Lacombe
### Adapted and modified by Anaïs Fournier on 05-14-2025
###############################

# Load required packages
library(readxl)    
library(SuessR)    
library(ggplot2)  
library(dplyr)     
library(writexl)    

# Import dataset
dat <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 2 - SuessR.csv")
dat[, 3:7] <- lapply(dat[, 3:7], as.numeric)

# Select and rename relevant columns for consistency
all <- dat[, c("FID", "Location", "Year", "d13C")]
colnames(all) <- c("id", "location", "year", "d13c")

# Remove rows containing missing values
all <- na.omit(all)

# Prepare dataset for SuessR correction
y <- all[, c("id", "year", "d13c")]
y$region <- "Subpolar North Atlantic"  # Closest marine region available for δ13C correction
y <- y[, c("id", "year", "region", "d13c")]  # Reorder columns as required by SuessR

# Apply Suess correction to δ13C values (standardized to 1950 atmospheric CO2 baseline)
x <- SuessR(y, correct.to = 1950)
x <- as.data.frame(x)  # Convert the result to a data frame

# Plot the magnitude of Suess correction over time
ggplot(x, aes(x = year)) +
  geom_smooth(aes(y = Suess.cor, color = "SuessR Correction"), method = "loess", se = FALSE, size = 1) +
  labs(title = "Suess effect correction over time",
       x = "Year",
       y = "Correction Value (‰)",
       color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot δ13C values before and after correction
ggplot(x, aes(x = year)) +
  geom_smooth(aes(y = d13c.uncor, color = "Before Correction"), method = "loess", se = FALSE, size = 1) +
  geom_smooth(aes(y = d13c.cor, color = "After SuessR Correction"), method = "loess", se = FALSE, size = 1) +
  labs(title = "δ13C values before and after Suess effect correction",
       x = "Year",
       y = "δ13C (‰)",
       color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Merge the corrected values back into the original dataset
# Keep only 'id' and 'd13c.cor' from the corrected dataset 'x'
corrected <- x[, c("id", "d13c.cor")]
colnames(corrected) <- c("FID", "d13CS")  # Rename columns to match 'dat' and name new corrected column

# Join the corrected d13C values (d13CS) into the original dataset 'dat'
dat <- left_join(dat, corrected, by = "FID")

# Remove the original d13C column
dat <- dat %>% select(-d13C) 

# Now 'dat' contains the corrected d13C values under the name 'd13CS'
# (If you prefer keeping both, just skip the select/rename part)

# Reorder columns to match original order
dat <- dat[, c(
  "FID", "Location", "Year", "Moisture", "Lipid", "d15N", "d13CS",
  "PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA",
  "PFHxS", "PFOS", "PFDS", "PFCA", "PFSA", "PFAS" 
)]

write.csv2(dat, "/Users/anaisfournier/Desktop/Stats/Step 3&4 - Correlation of Moisture, Lipid, d15N, d13CS and PFAS Concentrations.csv", row.names = FALSE)

