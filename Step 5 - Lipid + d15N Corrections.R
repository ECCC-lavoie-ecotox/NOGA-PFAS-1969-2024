###############################
### SI - Half a century of per- and polyfluoroalkyl substances in Northern Gannet eggs: Impact of regulations
### Step 5 - Lipid + d15N corrections
### Created by Anaïs Fournier on 05-14-2025
###############################

# Import dataset from Excel
data <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 5 - Lipid-Corrected Concentrations for d15N Corection.csv")
data[, 3:7] <- lapply(data[, 3:7], as.numeric)

# Non Lipid-corrected concentrations
data2 <- read.csv2("/Users/anaisfournier/Desktop/Stats/Step 3&4 - Correlation of Moisture, Lipid, d15N, d13CS and PFAS Concentrations.csv")
data2[, 3:7] <- lapply(data2[, 3:7], as.numeric)

# Remove all rows from the dataset where d15N is missing (NA)
data <- data[!is.na(data$d15N), ]

## The lipid + d15N-corrected concentration will be calculated using the following formula:
# More details are provided in the Supporting Information file (page S6)
# C_corrected = 10^((C_lipidcorrected + TMS_δ15N × (δ15N_mean - δ15N_sample)))
# where:
#   - C_lipidcorrected is the log10-transformed lipid corrected concentration
#   - δ15N_sample is the nitrogen isotope ratio of the sample
#   - δ15N_mean is the average δ15N value across all samples
#   - TMS_δ15N (Trophic Magnification Slope based on δ15N) must be estimated


## To associate the TMF-TL and TMS-d15N values only with the specified compounds
# Define the selected compounds
compounds <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFTrDA", "PFTeDA", 
               "PFHxS", "PFOS", "PFDS", "PFCA", "PFSA", "PFAS")

# Define the TMF-TL values for each selected compound (Miranda et al., 2021)
# For the sums, we calculated by summing the TMF-TL values of the individual compounds
tmf_tl_values <- c(3.2, 2.6, 4.8, 2.1, 1.7, 1.7, 1.6, 1.1, 5.2, 4.3, 
                   1.7, 2.8, 1.9)

# Calculate the TMS-d15N values using the formula: ABS(LOG10(TMF-TL) / 3.4)
tms_d15n_values <- abs(log10(tmf_tl_values) / 3.4)

# Create a data frame associating the compounds with the TMF-TL and TMS-d15N values
results <- data.frame(Compound = compounds, TMF_TL = tmf_tl_values, TMS_d15N = tms_d15n_values)

# Display the results
print(results)


## Lipid + d15N Corrections
# Extract the concentrations from columns 8 to 20 (compounds)
measured_concentrations <- data[, 8:20]

# Extract the δ15N values for the sample (column 'd15N')
d15N_sample <- data$d15N
d15N_sample <- as.numeric(data$d15N)

# Define the average value of δ15N (‰)
d15N_mean <- 14.88

# Repeat the TMS-d15N values for each compound for each row of the dataset
# This will ensure that each compound gets its corresponding TMS-d15N value
tms_d15n_expanded <- rep(tms_d15n_values, times = nrow(data))

# Reshape measured_concentrations to have one compound per column (for easier calculation)
measured_concentrations_matrix <- as.matrix(measured_concentrations)

# Apply log10 to the concentrations before calculation
log_measured_concentrations <- log10(measured_concentrations_matrix)

# Calculate the corrected concentrations (Ccorrected) using the formula for each compound
corrected_concentrations <- 10 ^ (log_measured_concentrations + 
                                    matrix(tms_d15n_expanded, ncol = ncol(measured_concentrations), byrow = TRUE) * 
                                    (d15N_mean - d15N_sample))

# Replace the original concentrations (columns 8 to 20) with the corrected concentrations
data[, 8:20] <- corrected_concentrations

# View the updated data frame with corrected concentrations
view(data)

write.csv2(data, "/Users/anaisfournier/Desktop/Stats/Step 6 - Lipid+d15N-Corrected Concentrations for GAMs.csv", row.names = FALSE)
