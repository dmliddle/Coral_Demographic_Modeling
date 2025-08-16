library(dplyr)
library(purrr)
library(ggplot2)

# Run the Forward Prediction Script 100 times to assess variation between runs

# WARNING: Running this can take several hours

n_simulations <- 100
set.seed(1234)

scenario3_sims <- vector("list", n_simulations)

for (i in 1:n_simulations) {
  scenario3_sims[[i]] <- simulate_coral_cover(
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
}

# Compute percent cover for each simulation
cover_results <- map(scenario3_sims, compute_percent_cover_by_class)

# Add simulation ID to each
cover_results <- map2(cover_results, 1:n_simulations, ~mutate(.x, SimID = .y))

# Combine all into one dataframe
cover_all <- bind_rows(cover_results)

cover_all <- cover_all %>%
  filter(Percent_Cover <= 100)

cover_summary <- cover_all %>%
  group_by(Year, Class) %>%
  summarize(
    Mean_Cover = mean(Percent_Cover, na.rm = TRUE),
    CI_Lower = quantile(Percent_Cover, 0.025, na.rm = TRUE),
    CI_Upper = quantile(Percent_Cover, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

cover_summary_smoothed <- cover_summary %>%
  group_by(Class) %>%
  arrange(Year) %>%
  mutate(
    Mean_Cover_smooth = loess(Mean_Cover ~ Year)$fitted,
    CI_Lower_smooth = loess(CI_Lower ~ Year)$fitted,
    CI_Upper_smooth = loess(CI_Upper ~ Year)$fitted
  ) %>%
  ungroup()


ggplot(cover_summary_smoothed, aes(x = Year, y = Mean_Cover_smooth, fill = Class, color = Class)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_Lower_smooth, ymax = CI_Upper_smooth), alpha = 0.2, color = NA) +
  scale_fill_manual(values = custom_colors) + 
  scale_x_continuous(breaks = seq(2020, 2030, by = 1), labels = as.character(seq(2020, 2030, by = 1))) +
  facet_wrap(~Class, scales = "free_y") +
  labs(
    title = "Smoothed Projected Coral Cover (Stress + Outplanting)",
    y = "Percent Cover",
    x = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

