
## FORWARD SIMULATION FRAMEWORK

library(readr)
library(dplyr)
library(brms)
library(ggplot2)
library(patchwork)
library(cowplot)
library(tidyr)
library(stringr)

# Load processed coral dataset and models
merged_df <- read_csv("~/Documents/Coral_Model/data/processed_data/forward_simulation/subset_2019.csv")
load("~/Documents/Coral_Model/models/coral_growth_model.RData")
load("~/Documents/Coral_Model/models/coral_mortality_model.RData")
load("~/Documents/Coral_Model/models/coral_recruitment_model.RData")

set.seed(123)

# Avg recruit size by class
avg_recruit_size_by_class <- merged_df %>%
  filter(Recruitment == 1) %>%
  group_by(Class) %>%
  summarize(Avg_Recruit_Size = mean(CoralSize, na.rm = TRUE), .groups = "drop")

simulate_coral_cover <- function(initial_df,
                                 recruitment_model,
                                 mortality_model,
                                 growth_rate_model,
                                 future_dhw_values,
                                 do_outplanting = FALSE,
                                 outplant_years = NULL,
                                 outplant_species = NULL,
                                 outplant_n = 0) {
  
  results_list <- list()
  
  # Ensure Class is character
  current_df <- initial_df %>%
    mutate(Source = "Survivor", Class = as.character(Class))
  
  for (year_idx in seq_along(future_dhw_values)) {
    year <- 2020 + (year_idx - 1)
    
    # Assign current and lagged DHW values
    current_df <- current_df %>%
      mutate(
        Year = year,
        Max_DHW = future_dhw_values[year_idx],
        Max_DHW_Lag1 = ifelse(year_idx >= 2, future_dhw_values[year_idx - 1], 0),
        Max_DHW_Lag2 = ifelse(year_idx >= 3, future_dhw_values[year_idx - 2], 0)
      )
    
    # Predict recruitment, mortality, and growth rate
    current_df$Predicted_Recruitment <- colMeans(predict(recruitment_model, newdata = current_df, allow_new_levels = TRUE, summary = FALSE))
    current_df$Predicted_Mortality <- colMeans(predict(mortality_model, newdata = current_df, allow_new_levels = TRUE, summary = FALSE))
    current_df$Predicted_GrowthRate <- colMeans(predict(growth_rate_model, newdata = current_df, allow_new_levels = TRUE, summary = FALSE))
    
    # Apply mortality and update CoralSize
    current_df <- current_df %>%
      mutate(
        Survived = runif(n()) > Predicted_Mortality,
        GrowthRate = Predicted_GrowthRate
      ) %>%
      filter(Survived) %>%
      mutate(
        CoralSize = CoralSize * exp(GrowthRate),
        log_CoralSize = log(CoralSize)
      )
    
    # Generate new recruits
    new_recruits <- current_df %>%
      filter(runif(n()) < Predicted_Recruitment) %>%
      left_join(avg_recruit_size_by_class, by = "Class") %>%  # ✅ safer than match()
      mutate(
        GenSiteID = paste0("Recruit_", year, "_", row_number()),
        CoralSize = Avg_Recruit_Size,
        log_CoralSize = log(CoralSize),
        Year = year,
        Source = "Recruit"
      ) %>%
      select(-Avg_Recruit_Size)
    
    # Predict GrowthRate for recruits
    if (nrow(new_recruits) > 0) {
      new_recruits$GrowthRate <- colMeans(
        predict(growth_rate_model, newdata = new_recruits, allow_new_levels = TRUE, summary = FALSE)
      )
    }
    
    # Generate outplants
    if (do_outplanting && year %in% outplant_years) {
      acropora_size <- avg_recruit_size_by_class %>%
        filter(Class == outplant_species) %>%
        pull(Avg_Recruit_Size)
      
      outplant_recruits <- tibble(
        Class = outplant_species,
        GenSiteID = paste0("Outplant_", year, "_", 1:outplant_n),
        CoralSize = rep(acropora_size, outplant_n),
        log_CoralSize = log(acropora_size),
        Max_DHW = future_dhw_values[year_idx],
        Max_DHW_Lag1 = ifelse(year_idx >= 2, future_dhw_values[year_idx - 1], 0),
        Max_DHW_Lag2 = ifelse(year_idx >= 3, future_dhw_values[year_idx - 2], 0),
        Year = year,
        Source = "Outplant"
      )
      
      outplant_recruits$GrowthRate <- colMeans(
        predict(growth_rate_model, newdata = outplant_recruits, allow_new_levels = TRUE, summary = FALSE)
      )
      
    } else {
      outplant_recruits <- tibble()  # ✅ empty tibble prevents NA coercion
    }
    
    # Combine survivors, recruits, and outplants safely
    current_df <- bind_rows(current_df, new_recruits, outplant_recruits)
    
    results_list[[as.character(year)]] <- current_df
  }
  
  final_df <- bind_rows(results_list)
  return(final_df)
}



# OBSERVED STRESS + Predicted Bleaching
observed_dhw <- sapply(2020:2030, function(y) {
  if (y == 2020) return(0.3)
  if (y == 2023) return(17.2)
  if (y == 2024) return(17.1)
  if (y == 2025) return(3.94)
  if (y == 2026) return(1.02)
  if (y == 2027) return(3.05)
  if (y == 2028) return(12.50)
  if (y == 2029) return(9.46)
  if (y == 2030) return(7.03)
  return(0)
})


# OBSERVED STRESS + Worst Case Bleaching
observed_dhw_worst <- sapply(2020:2030, function(y) {
  if (y == 2020) return(0.3)
  if (y == 2023) return(17.2)
  if (y == 2024) return(17.1)
  if (y == 2025) return(9.82)
  if (y == 2026) return(6.9)
  if (y == 2027) return(8.95)
  if (y == 2028) return(18.37)
  if (y == 2029) return(15.36)
  if (y == 2030) return(12.92)
  return(0)
})

# OBSERVED STRESS + Best Case Bleaching
observed_dhw_best <- sapply(2020:2030, function(y) {
  if (y == 2020) return(0.3)
  if (y == 2023) return(17.2)
  if (y == 2024) return(17.1)
  if (y == 2025) return(0)
  if (y == 2026) return(0)
  if (y == 2027) return(0)
  if (y == 2028) return(6.62)
  if (y == 2029) return(3.55)
  if (y == 2030) return(1.13)
  return(0)
})


# Create a DHW vector of 0s for 2020–2030
baseline_dhw <- rep(0, length(2020:2030))

future_no_stress <- simulate_coral_cover(
  initial_df = merged_df,
  recruitment_model = coral_recruitment_model,
  mortality_model = coral_mortality_model,
  growth_rate_model = coral_growth_model,
  future_dhw_values = baseline_dhw,
  do_outplanting = FALSE
)

future_with_stress <- simulate_coral_cover(
  initial_df = merged_df,
  recruitment_model = coral_recruitment_model,
  mortality_model = coral_mortality_model,
  growth_rate_model = coral_growth_model,
  future_dhw_values = observed_dhw,
  do_outplanting = FALSE
)

future_outplanted_stress <- simulate_coral_cover(
  initial_df = merged_df,
  recruitment_model = coral_recruitment_model,
  mortality_model = coral_mortality_model,
  growth_rate_model = coral_growth_model,
  future_dhw_values = observed_dhw,
  do_outplanting = TRUE,
  outplant_years = 2025:2030,
  outplant_species = "coral_A",
  outplant_n = 200
)

future_outplanted_stress_best <- simulate_coral_cover(
  initial_df = merged_df,
  recruitment_model = coral_recruitment_model,
  mortality_model = coral_mortality_model,
  growth_rate_model = coral_growth_model,
  future_dhw_values = observed_dhw_best,
  do_outplanting = TRUE,
  outplant_years = 2025:2030,
  outplant_species = "coral_A",
  outplant_n = 200
)

future_outplanted_stress_worst <- simulate_coral_cover(
  initial_df = merged_df,
  recruitment_model = coral_recruitment_model,
  mortality_model = coral_mortality_model,
  growth_rate_model = coral_growth_model,
  future_dhw_values = observed_dhw_worst,
  do_outplanting = TRUE,
  outplant_years = 2025:2030,
  outplant_species = "coral_A",
  outplant_n = 200
)

# Total survey area
total_survey_area <- 7 * 10000  # 70,000 m²

compute_percent_cover_by_class <- function(df) {
  df %>%
    group_by(Year, Class) %>%
    summarize(
      Total_Coral_Area = sum(CoralSize, na.rm = TRUE) / 100,  # cm² to m²
      Percent_Cover = (Total_Coral_Area / total_survey_area) * 100,
      .groups = "drop"
    )
}

forecast_data1 <- compute_percent_cover_by_class(future_no_stress)
forecast_data2 <- compute_percent_cover_by_class(future_with_stress)
forecast_data3 <- compute_percent_cover_by_class(future_outplanted_stress)
forecast_data4 <- compute_percent_cover_by_class(future_outplanted_stress_best)
forecast_data5 <- compute_percent_cover_by_class(future_outplanted_stress_worst)

#write.csv(forecast_data1, file = "~/Documents/Masters_Project/Demographic_Model/no_outplant_no_stress.csv")
#write.csv(forecast_data2, file = "~/Documents/Masters_Project/Demographic_Model/no_outplant_stress.csv")
#write.csv(forecast_data3, file = "~/Documents/Masters_Project/Demographic_Model/outplant_stress.csv")

forecast_data1 <- forecast_data1 %>%
  mutate(Scenario = "No Stress, No Outplanting")
forecast_data2 <- forecast_data2 %>%
  mutate(Scenario = "Stress, No Outplanting")
forecast_data3 <- forecast_data3 %>%
  mutate(Scenario = "Stress + Outplanting")
forecast_data4 <- forecast_data4 %>%
  mutate(Scenario = "Best Case Stress + Outplanting")
forecast_data5 <- forecast_data5 %>%
  mutate(Scenario = "Worst Case Stress + Outplanting")

# Combine all into one dataset
all_data <- bind_rows(forecast_data1, forecast_data2, forecast_data3, forecast_data4, forecast_data5)

# Summarize species-level cover for 2020 and 2030
species_cover_summary <- all_data %>%
  filter(Year %in% c(2020, 2030)) %>%
  group_by(Scenario, Class, Year) %>%
  summarize(Percent_Cover = sum(Percent_Cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Percent_Cover, names_prefix = "Cover_") %>%
  rename(Beginning = Cover_2020, End = Cover_2030)

# Summarize total cover per scenario (summing across all classes)
total_cover_summary <- all_data %>%
  filter(Year %in% c(2020, 2030)) %>%
  group_by(Scenario, Year) %>%
  summarize(Percent_Cover = sum(Percent_Cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Percent_Cover, names_prefix = "Cover_") %>%
  mutate(Class = "Total") %>%
  rename(Beginning = Cover_2020, End = Cover_2030)


## VISUALIZE RESULTS


custom_colors <- c(
  "coral_A" = "#f35d59",
  "coral_B" = "#c87e02",
  "coral_C" = "#819c03",
  "coral_D" = "#1eb231",
  "coral_E" = "#2dbc96",
  "coral_F" = "#13aadc",
  "coral_G" = "#4d85ff",
  "coral_H" = "#d052fa",
  "coral_I" = "#fb45b6"
)

max_cover <- max(all_data %>% group_by(Year, Scenario) %>%
                   summarize(total = sum(Percent_Cover), .groups = "drop") %>%
                   pull(total), na.rm = TRUE)

shared_y <- scale_y_continuous(
  limits = c(0, ceiling(max_cover / 5) * 5),
  labels = scales::percent_format(scale = 1)
)


# Plot!
p1 <- ggplot(forecast_data1, aes(x = Year, y = Percent_Cover, fill = Class)) +
  geom_area(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  labs(
    title = "No Stress, No Outplanting",
    x = "Year", y = "Percent Cover (%)", fill = "Coral Species"
  ) +
  theme(legend.position = "right", text = element_text(size = 14))

p2 <- ggplot(forecast_data2, aes(x = Year, y = Percent_Cover, fill = Class)) +
  geom_area(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  labs(
    title = "With Stress, No Outplanting",
    x = "Year", y = "Percent Cover (%)", fill = "Coral Species"
  ) +
  theme(legend.position = "right", text = element_text(size = 14))

p3 <- ggplot(forecast_data3, aes(x = Year, y = Percent_Cover, fill = Class)) +
  geom_area(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  labs(
    title = "Acropora Outplanting with Predicted Stress",
    x = "Year", y = "Percent Cover (%)", fill = "Coral Species"
  ) +
  theme(legend.position = "right", text = element_text(size = 14))

p4 <- ggplot(forecast_data4, aes(x = Year, y = Percent_Cover, fill = Class)) +
  geom_area(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  labs(
    title = "Acropora Outplanting with Low Stress",
    x = "Year", y = "Percent Cover (%)", fill = "Coral Species"
  ) +
  theme(legend.position = "right", text = element_text(size = 14))

p5 <- ggplot(forecast_data5, aes(x = Year, y = Percent_Cover, fill = Class)) +
  geom_area(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  labs(
    title = "Acropora Outplanting with High Stress",
    x = "Year", y = "Percent Cover (%)", fill = "Coral Species"
  ) +
  theme(legend.position = "right", text = element_text(size = 14))


# Shared theme and axis settings
p1 <- p1 + shared_y + theme(legend.position = "none")
p2 <- p2 + shared_y + theme(axis.title.y = element_blank(), legend.position = "none")
p3 <- p3 + shared_y + theme(axis.title.y = element_blank(), legend.position = "none")

# Stack p1, p2, p3
combined_plot_1 <- (p1 / p2 / p3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.margin = margin(t = 5),
    legend.box.margin = margin(t = -5)
  )

# Add y-axis label
final_plot_1 <- ggdraw() +
  draw_plot(combined_plot_1, x = 0.05, y = 0, width = 0.95, height = 1) +
  draw_label("Percent Cover (%)",
             x = 0.02, y = 0.5,
             angle = 90,
             fontface = "bold",
             size = 14,
             hjust = 0.5)

# Show or save
final_plot_1

# Apply consistent axis scale and formatting
p4 <- p4 + shared_y + theme(legend.position = "none")
p5 <- p5 + shared_y + theme(axis.title.y = element_blank())

# Stack p4 and p5 with shared legend
combined_plot_2 <- (p4 / p5) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.margin = margin(t = 5),
    legend.box.margin = margin(t = -5)
  )

# Add shared y-axis
final_plot_2 <- ggdraw() +
  draw_plot(combined_plot_2, x = 0.05, y = 0, width = 0.95, height = 1) +
  draw_label("Percent Cover (%)",
             x = 0.02, y = 0.5,
             angle = 90,
             fontface = "bold",
             size = 14,
             hjust = 0.5)

# Show or save
final_plot_2

## IN THIS SECTION WE SHOW WHAT THE PREDICTED SURVIVAL OF OUTPLANTS ARE

future_outplants <- future_outplanted_stress %>%
  filter(Source == "Outplant") %>%
  count(Year)

# Step 1: Count cumulative number of outplanted corals by year
# Assuming 100 corals are planted every year starting in 2025
outplant_summary <- future_outplanted_stress %>%
  filter(Source == "Outplant") %>%
  group_by(Year) %>%
  summarize(Surviving_Outplants = n(), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(
    Total_Outplanted = (Year - 2024) * 200,  # 200 outplanted starting in 2025
    Survival_Rate = Surviving_Outplants / Total_Outplanted
  )


ggplot(outplant_summary, aes(x = Year, y = Survival_Rate)) +
  geom_line(size = 1.2, color = "#4d85ff") +
  geom_point(size = 2, color = "#4d85ff") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 2025:2030) +
  labs(
    title = "Cumulative Survival Rate of Outplanted Corals",
    x = "Year",
    y = "Survival Rate (%)"
  ) +
  theme(
    text = element_text(size = 14)
  )

cohort_survival <- future_outplanted_stress %>%
  # Step 1: Extract cohort year from ID
  mutate(
    Outplant_Year = as.numeric(str_extract(GenSiteID, "(?<=Outplant_)\\d{4}")),
    Age = Year - Outplant_Year
  )

# Step 2: Count survivors per cohort per age
survival_summary <- cohort_survival %>%
  group_by(Outplant_Year, Age) %>%
  summarise(
    Survivors = n(), .groups = "drop"
  ) %>%
  arrange(Outplant_Year, Age)

initial_cohorts <- data.frame(
  Outplant_Year = 2025:2030,
  Initial_N = 200
)

# Step 3: Merge and calculate survival rate
survival_summary <- survival_summary %>%
  left_join(initial_cohorts, by = "Outplant_Year") %>%
  group_by(Outplant_Year) %>%
  mutate(Survival_Rate = Survivors / Initial_N) %>%
  ungroup()

