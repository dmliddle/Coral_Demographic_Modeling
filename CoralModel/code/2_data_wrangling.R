# Wrangling Coral and Temperature Datasets

# Extraction of demographic data
# install.packages('reshape2')
library(readr)
library(dplyr)
library(zoo)
library(ggplot2)
library(tidyr)

# ---- Step 1: Read in and Combine Data Files ----

# Read in ocean temperature data
dhw_data <- read_csv("~/Documents/Coral_Model/data/raw_temp_data/northern_line_islands_full_dataset.csv")

# Set the folder path for coral data
folder_path <- "~/Documents/Coral_Model/data/processed_data/species_data"

# Get list of CSV files
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and add Class column for each file
data_list <- lapply(csv_files, function(file) {
  
  # Read the CSV
  df <- read.csv(file)
  
  # Extract the species name from the filename
  # Example: final_merge_Pocillopora.csv -> Pocillopora
  species_name <- sub("final_merge_(.*)\\.csv", "\\1", basename(file))
  
  # Add it as a new column
  df$Class <- species_name
  
  # Move Class to 4th column (makes looking at data easier later)
  df <- df[, c(1:3, ncol(df), 4:(ncol(df)-1)) ]
  
  return(df)
})

# Combine into one dataframe
combined_data <- do.call(rbind, data_list)

#write.csv(combined_data, "~/Documents/Coral_Model/data/processed_data/models/combined_data.csv", row.names = FALSE)


# ---- Step 2: Pivot to Long Format ----

# Pivot data into long format
coral_long <- combined_data %>%
  pivot_longer(
    cols = starts_with("size"),    # Select size columns
    names_to = "Year",             # Create a Year column
    names_prefix = "size_",        # Remove prefix from column names
    values_to = "CoralSize"        # Values go into CoralSize column
  ) %>%
  # Add log-transformed CoralSize for modeling
  mutate(
    # Adding a small offset (1) to avoid log(0) = -Inf if any colonies have size 0
    log_CoralSize = log(CoralSize)
  )

coral_long$log_CoralSize[is.infinite(coral_long$log_CoralSize)] <- log(0.01)

# Remove the 1st, 3rd, and 5th through 10th columns
coral_long <- coral_long %>%
  dplyr::select(-c(3,5:10))

# ---- Step 3: Add Growth, Recruitment, and Mortality Factors ----

coral_long <- coral_long %>%
  arrange(Year, GenSiteID) %>%
  group_by(GenSiteID) %>%
  dplyr::mutate(
    Recruitment = ifelse(lag(CoralSize) == 0 & CoralSize > 0, 1, 0), # Newly appeared coral
    Mortality = ifelse(CoralSize == 0 & lag(CoralSize) > 0, 1, 0),    # Coral that disappeared
    GrowthRate = ifelse(lag(CoralSize) > 0, log(CoralSize / lag(CoralSize)), NA)
  )

# Handle NAs and create next year's mortality column
coral_long <- coral_long %>%
  group_by(GenSiteID) %>%  # Ensure lead() works per colony
  mutate(
    # Handle missing growth rates by setting them to 0
    GrowthRate = ifelse(is.na(GrowthRate), 0, GrowthRate),
    
    # For the first year (2013), set mortality and recruitment to 0
    Mortality = ifelse(Year == 2013, 0, Mortality),
    Recruitment = ifelse(Year == 2013, 0, Recruitment),
    
    # Create a column for next year's mortality
    # This is useful for predicting survival in the following year
    Next_Year_Mortality = lead(Mortality, default = 0)
  ) %>%
  ungroup()  # Shift mortality status to next year, this is so that we can capture the size of the coral before mortality, otherwise the size will just be 0, every time.

# Change the Year column from a Character to Numeric

coral_long <- coral_long %>%
  mutate(Year = as.numeric(Year))

# Save the merged dataset (if needed)

#write.csv(coral_long, "~/Documents/Coral_Model/data/processed_data/coral_long.csv", row.names = FALSE)

# ---- Step 4: Merge in DHW ----

dhw_yearly <- dhw_data %>%
  group_by(Year = as.character(YYYY)) %>%
  summarise(Max_DHW = max(DHW, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year = as.numeric(Year))

merged_data <- coral_long %>%
  left_join(dhw_yearly, by = "Year")

# ---- Step 5: Create Lagged Variables ----

merged_data <- merged_data %>%
  arrange(Year, GenSiteID) %>%
  group_by(GenSiteID) %>%  # Group by GenSiteID to calculate lags separately for each coral genet
  mutate(
    # Lag predictors
    Max_DHW_Lag1 = lag(Max_DHW, 1),   # Previous year's heat stress
    Max_DHW_Lag2 = lag(Max_DHW, 2),   # Two years prior
    Mortality_Lag1 = lag(Mortality, 1) # Last year's mortality
  ) %>%
  ungroup() %>%
  mutate(
    # Handle NA values in early years for DHW lags
    Max_DHW_Lag1 = ifelse(Year %in% c(2013, 2014) & is.na(Max_DHW_Lag1), 0, Max_DHW_Lag1),
    Max_DHW_Lag2 = ifelse(Year %in% c(2013, 2014) & is.na(Max_DHW_Lag2), 0, Max_DHW_Lag2),
    
    # New: Add FirstYear indicator and clean Mortality_Lag1
    FirstYear = ifelse(Year == 2013, 1, 0),
    Mortality_Lag1 = ifelse(is.na(Mortality_Lag1), 0, Mortality_Lag1)
  )

# Add a 
merged_data <- merged_data %>%
  group_by(GenSiteID) %>%
  mutate(
    # First year the colony recruited
    first_recruitment_year = ifelse(
      any(Recruitment > 0 & CoralSize > 0), 
      min(Year[Recruitment > 0 & CoralSize > 0], na.rm = TRUE), 
      NA
    ),
    
    # First year the colony died
    first_mortality_year = ifelse(
      any(Mortality == 1), 
      min(Year[Mortality == 1], na.rm = TRUE), 
      NA
    )
  ) %>%
  ungroup()

# ---- Step 6: Add Site-Level Lagged Mortality ----

merged_data <- merged_data %>%
  # First, calculate site-year mortality proportion
  group_by(Site, Year) %>%
  mutate(
    SiteMortality_Prop = mean(Mortality, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  # Next, calculate lagged site-level mortality
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    SiteMortality_Prop_Lag1 = lag(SiteMortality_Prop, 1),
    SiteMortality_Prop_Lag1 = ifelse(is.na(SiteMortality_Prop_Lag1), 0, SiteMortality_Prop_Lag1)
  ) %>%
  ungroup()

#write.csv(merged_data, "~/Documents/Coral_Model/data/processed_data/models/merged_data.csv", row.names = FALSE)


# ---- Step 7: Growth Model Dataset ----

# We create the dataframe for the growth model by subsetting values in which the coral size is non-zero (not dead or pre-recruited),
# Excludes mortality events, and initial recruitment growth (0 → first size), keeps growth modeling focused on established colonies.

coral_growth_df <- merged_data %>%
  filter(CoralSize > 0, Mortality == 0, Recruitment == 0, Year != 2013)

#write.csv(merged_data, "~/Documents/Coral_Model/data/processed_data/models/growth_data.csv", row.names = FALSE)

# ---- Step 8: Recruitment Dataset ----

# - Keep 0's before first recruitment and the first 1 for recruiting colonies
# - Keep only the first 0 for colonies that never recruited (avoids oversampling zeros)
# - Ensures each GenSiteID contributes a single trial without redundant post-recruitment zeros

recruitment_df <- merged_data %>%
  arrange(GenSiteID, Year) %>%
  group_by(GenSiteID) %>%
  filter(
    # Include only the first year this genet appears
    Year == min(Year[CoralSize > 0], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # Recruitment is 1 if CoralSize > 0 and CoralSize in previous year was 0 or missing
    Recruitment = ifelse(CoralSize > 0 & (lag(CoralSize) == 0 | is.na(lag(CoralSize))), 1, 0)
  ) %>%
  filter(!is.na(Recruitment))

#write.csv(merged_data, "~/Documents/Coral_Model/data/processed_data/models/recruitment_data.csv", row.names = FALSE)

# ---- Step 9: Mortality Dataset ----

mortality_df <- merged_data %>%
  group_by(GenSiteID) %>%
  filter(
    (Year >= first_recruitment_year | is.na(first_recruitment_year)) |
      (Year == 2013 & CoralSize > 0)
  ) %>%
  # Keep rows *up to the year before the first mortality event*, because that's when Next_Year_Mortality == 1
  filter(
    is.na(first_mortality_year) | Year < first_mortality_year
  ) %>%
  ungroup()


#write.csv(merged_data, "~/Documents/Coral_Model/data/processed_data/models/mortality_data.csv", row.names = FALSE)

# ---- Step 10: Subset for 2019 Alive Colonies ----

subset_2019 <- merged_data %>%
  filter(Year == 2019, CoralSize > 0)

# Write CSV (If Needed)
#write.csv(subset_2019, "~/Documents/Coral_Model/data/processed_data/forward_simulation/subset_2019.csv", row.names = FALSE)

# ---- Step 10: Obtain Summary Stats to Inform Priors ----

# Stats for Growth Model
summary_stats <- coral_growth_df %>%
  summarise(
    mean_growth = mean(GrowthRate, na.rm = TRUE),
    sd_growth = sd(GrowthRate, na.rm = TRUE),
    median_growth = median(GrowthRate, na.rm = TRUE),
    iqr_growth = IQR(GrowthRate, na.rm = TRUE)
  )
summary_stats

# Recruitment Statistics
# Proportion of colonies that ever recruit (1 at least once)
baseline_recruitment <- recruitment_df %>%
  group_by(GenSiteID) %>%
  summarise(ever_recruited = max(Recruitment, na.rm = TRUE)) %>%
  summarise(
    proportion_recruiting = mean(ever_recruited)
  )

print(baseline_recruitment)

# Predictor Ranges (Min, Max, Median)
predictor_summary <- recruitment_df %>%
  summarise(
    # Max Degree Heating Week (DHW)
    Max_DHW_min = min(Max_DHW, na.rm = TRUE),
    Max_DHW_median = median(Max_DHW, na.rm = TRUE),
    Max_DHW_max = max(Max_DHW, na.rm = TRUE),
    
    # Max DHW with 2-year lag
    Max_DHW_Lag2_min = min(Max_DHW_Lag2, na.rm = TRUE),
    Max_DHW_Lag2_median = median(Max_DHW_Lag2, na.rm = TRUE),
    Max_DHW_Lag2_max = max(Max_DHW_Lag2, na.rm = TRUE),
    
    # Mortality Lag 1 (binary, so proportion is useful)
    Mortality_Lag1_min = min(Mortality_Lag1, na.rm = TRUE),
    Mortality_Lag1_prop = mean(Mortality_Lag1, na.rm = TRUE),
    Mortality_Lag1_max = max(Mortality_Lag1, na.rm = TRUE)
  )

print(predictor_summary)

# Recruitment Probability by Mortality
# Helps visualize effect size for Mortality_Lag1
recruitment_by_mortality <- recruitment_df %>%
  group_by(Mortality_Lag1) %>%
  summarise(
    n = n(),
    mean_recruitment = mean(Recruitment, na.rm = TRUE)
  )

print(recruitment_by_mortality)

# ---- 1. Check predictor ranges ----
mortality_df %>%
  summarise(
    min_log_CoralSize = min(log_CoralSize, na.rm = TRUE),
    max_log_CoralSize = max(log_CoralSize, na.rm = TRUE),
    mean_log_CoralSize = mean(log_CoralSize, na.rm = TRUE),
    
    min_DHW = min(Max_DHW, na.rm = TRUE),
    max_DHW = max(Max_DHW, na.rm = TRUE),
    mean_DHW = mean(Max_DHW, na.rm = TRUE),
    
    min_DHW_L1 = min(Max_DHW_Lag1, na.rm = TRUE),
    max_DHW_L1 = max(Max_DHW_Lag1, na.rm = TRUE),
    mean_DHW_L1 = mean(Max_DHW_Lag1, na.rm = TRUE),
    
    min_DHW_L2 = min(Max_DHW_Lag2, na.rm = TRUE),
    max_DHW_L2 = max(Max_DHW_Lag2, na.rm = TRUE),
    mean_DHW_L2 = mean(Max_DHW_Lag2, na.rm = TRUE)
  )

# ---- 2. Check baseline mortality rate ----
mortality_df %>%
  summarise(
    overall_mortality_rate = mean(Next_Year_Mortality, na.rm = TRUE)
  )

# ---- 3. Mortality rate by stress level and size ----
# This gives a quick feel for effect sizes
mortality_df %>%
  mutate(size_class = cut(log_CoralSize, breaks = quantile(log_CoralSize, probs = seq(0, 1, 0.25), na.rm = TRUE),
                          include.lowest = TRUE),
         stress_class = cut(Max_DHW, breaks = quantile(Max_DHW, probs = seq(0, 1, 0.25), na.rm = TRUE),
                            include.lowest = TRUE)) %>%
  group_by(size_class, stress_class) %>%
  summarise(
    mean_mortality = mean(Next_Year_Mortality, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

## WRITE CSV IF NEEDED

