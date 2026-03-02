#############################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 1 - Dealing with non-detects
### Script originally created by Rose Lacombe
### Modified by Anaïs Fournier on 05-14-2025
#############################

# Load required packages
library(dplyr)
library(stringr)
library(NADA)
library(NADA2)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)

# Function to censor concentration values based on detection limits.
# Flags non-detect values (denoted by '<' or NA) with a Boolean indicator variable.
censor <- function(data, t) {
  # Store column names for dynamic naming later
  cols <- colnames(data)
  
  # Iterate through compound columns starting from column 't'
  for (i in t:ncol(data)) {
    
    # Create a censoring indicator column: TRUE if value is below detection limit ('<') or missing (NA), FALSE otherwise
    data$new <- if_else(grepl("<", data[, i]), TRUE,
                        if_else(is.na(data[, i]), TRUE, FALSE))
    
    # Rename the censoring indicator column using the compound name followed by "cen"
    names(data)[names(data) == 'new'] <- paste0(cols[[i]], "cen")
    
    # Remove the '<' character to facilitate numeric conversion
    data[, i] <- str_replace(data[, i], "<", "")
    
    # Convert values to numeric for further statistical processing
    data[, i] <- as.numeric(as.character(data[, i]))
  }
  
  return(data)
}


# Function to impute censored (non-detect) values using the Regression on Order Statistics (ROS) method.
# Performs imputation separately for each group defined by the 'Year' variable.
impute <- function(data, h) {
  
  # Extract covariates for later reintegration (e.g., for statistical modeling)
  SI <- subset(data, select = c(FID, Location, Year, Moisture, Lipid, d15N, d13C))  # Covariates to retain
  
  # Remove covariate columns to focus on compound and censoring data
  dat2 <- subset(data, select = -c(Moisture, Lipid, d15N, d13C))
  
  # Store column names for dynamic processing
  cols <- colnames(dat2)
  
  # Split data into groups by 'Year' for separate ROS modeling
  colony <- dat2 %>% group_by(Year) %>% group_split()
  
  # Initialize an empty long-format data frame to store results
  stats <- data.frame(FID = as.character(),
                      Location = as.character(),
                      Year = as.numeric(),
                      Compound = as.character(),
                      Concentration = as.numeric())
  
  # Loop over each 'Year' group
  for (i in 1:length(colony)) {
    
    loc <- colony[[i]]  # Isolate the group
    
    # Determine the split point between concentration values and censor indicators
    k <- ifelse(ncol(loc) %% 2 == 0, round(ncol(loc)/2), round(ncol(loc)/2) + 1)
    
    # Loop over compound columns only (excluding censor flags)
    for (f in h:k) {
      
      # Identify the corresponding censor indicator column for the compound
      x <- grep(cols[f], colnames(loc), value = FALSE)
      
      # Extract FID, compound values, and censor indicators
      comp <- as.data.frame(loc[, c(1, x)])
      colnames(comp) <- c("FID", "Compound", "Cen")
      
      # Order for compatibility with cenros function
      comp2 <- comp[order(comp$Compound, decreasing = FALSE, na.last = FALSE),]
      
      # Skip groups with only one observation (ROS requires at least two)
      d <- nrow(loc)
      e <- d - 1
      if (d < 2) { print(loc$FID) }
      
      # Skip compounds with only non-detects or missing values
      if (is.na(sum(comp[, 2])) == TRUE) next
      if (length(grep("TRUE", comp2[, 3])) >= e) next
      
      # Apply ROS method to model censored data
      compros <- cenros(comp2[, 2], comp2[, 3])
      
      # Create a tidy-format data frame to store imputed concentrations
      add <- data.frame(FID = comp2$FID,
                        Location = loc$Location,
                        Year = loc$Year,
                        Compound = rep(cols[f], times = d),
                        Concentration = dput(as.numeric(compros$modeled)))
      
      # Append results to the main statistics data frame
      stats <- rbind(stats, add)
    }
  }
  
  # Reshape the long-format imputed dataset to a wide format
  interp <- stats %>%
    group_by(Location, Year, Compound) %>%
    pivot_wider(id_cols = c(FID, Location, Year), 
                names_from = Compound, values_from = Concentration)
  
  # Merge imputed data with original covariates
  final <- merge(SI, interp, by = c("FID", "Location", "Year"))
  
  return(final)
}

# Load dataset
data <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 1 - Dealing with non-detects.csv", stringsAsFactors = FALSE)

# Apply censoring function from column 8 onwards
PFAScen <- censor(data, 8)

# Apply imputation function starting from column 4
PFASint <- impute(PFAScen, 4)

# --- Special case handling for individual compounds such as PFDS ---
# Note: The main imputation loop may fail to process certain compounds even when enough detected values are available. 
# Therefore, it is important to verify manually whether all compounds have been correctly processed. 
# If not, reprocess them separately as below.

# Handle PFDS separately if missing
data2 <- data.frame(
  FID = data$FID,
  Location = data$Location,
  Year = data$Year,
  Moisture = data$Moisture,
  Lipid = data$Lipid,
  d15N = data$d15N,
  d13C = data$d13C,
  PFOS = data$PFOS,
  PFDS = data$PFDS
)

data2cen <- censor(data2, 8)
data2int <- impute(data2cen, 4)

# Merge the separately imputed PFDS values into the main dataset
PFASint <- merge(PFASint, data2int[, c("FID", "PFDS")], by = "FID", all.x = TRUE)

# Reorder columns to match original order (compound columns might be incomplete if fully censored)
PFASint <- PFASint[, c(
  "FID", "Location", "Year", "Moisture", "Lipid", "d15N", "d13C",
  "PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA",
  "PFHxS", "PFOS", "PFDS"
)]

# Calculate summed variables per row
PFASint$PFCA <- rowSums(PFASint[, c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA")], na.rm = TRUE)
PFASint$PFSA <- rowSums(PFASint[, c("PFHxS", "PFOS", "PFDS")], na.rm = TRUE)
PFASint$PFAS <- rowSums(PFASint[, c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA",
                                    "PFHxS", "PFOS", "PFDS")], na.rm = TRUE)

write.csv2(PFASint, file = "/Users/anaisfournier/Desktop/Stats/Step 2 - SuessR.csv", row.names = FALSE)
