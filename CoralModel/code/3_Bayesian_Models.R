#### RUN BAYESIAN MODELS

# Load libraries
library(brms)

# Read in demographic dataframes

coral_growth_df <- read_csv("data/processed_data/models/growth_data.csv")
recruitment_df <- read_csv("data/processed_data/models/recruitment_data.csv")
mortality_df <- read_csv("data/processed_data/models/mortality_data.csv")

## GROWTH MODEL

coral_growth_model <- brm(
  formula = GrowthRate ~ Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + log_CoralSize + (1 | Class),
  data = coral_growth_df,
  family = student(),
  prior = c(
    # Intercept ~ mean log-growth rate (≈0) with SD ≈ 0.75
    prior(normal(0, 0.75), class = "Intercept"),
    
    # Slopes ~ weakly informative
    # - A 1-unit change in predictors rarely changes log-growth > ~0.5
    prior(normal(0, 0.3), class = "b"),
    
    # Random species (Class) intercept SD ~ allows ~95% range ±2 log units
    prior(exponential(1), class = "sd"),
    
    # Student-t ν (controls tail heaviness)
    prior(gamma(5, 1), class = "nu")
  ),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

# Generate a summary table rounded out to 6 decimal places (if necessary)
#posterior_summary(coral_growth_model, probs = c(0.025, 0.975)) %>% round(6)

# Check Posterior Distribution if Necessary
pp_check(coral_growth_model, ndraws = 100) +
  coord_cartesian(xlim = c(-5, 5))

# Save Model (if necessary)
#save(coral_growth_model, file =  "~/Documents/Coral_Model/models/coral_growth_model.RData")

### RECRUITMENT MODEL

coral_recruitment_model <- brm(
  Recruitment ~ Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + SiteMortality_Prop_Lag1 + (1 | Class) + (1 | GenSiteID),
  family = bernoulli(link = "logit"),
  data = recruitment_df,
  prior = c(
    prior(normal(-4, 1), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

# Generate a summary table rounded out to 6 decimal places (if necessary)
#posterior_summary(coral_recruitment_model, probs = c(0.025, 0.975)) %>% round(6)

# Save Model (if necessary)
#save(coral_recruitment_model, file =  "~/Documents/Coral_Model/models/coral_recruitment_model.RData")

### MORTALITY MODEL

coral_mortality_model <- brm(
  Next_Year_Mortality ~ log_CoralSize + Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + (1 | Class),
  family = bernoulli(),
  data = mortality_df,
  prior = c(
    prior(normal(-0.85, 0.5), class = "Intercept"),
    prior(normal(-0.4, 0.3), class = "b", coef = "log_CoralSize"),
    prior(normal(0.02, 0.02), class = "b", coef = "Max_DHW"),
    prior(normal(0.02, 0.02), class = "b", coef = "Max_DHW_Lag1"),
    prior(normal(0.02, 0.02), class = "b", coef = "Max_DHW_Lag2"),
    prior(cauchy(0, 1), class = "sd")
  ),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  control = list(adapt_delta = 0.90, max_treedepth = 12)
)

# Generate a summary table rounded out to 6 decimal places (if necessary)
#posterior_summary(coral_mortality_model, probs = c(0.025, 0.975)) %>% round(6)

# Save Model (if necessary)
#save(coral_mortality_model, file =  "~/Documents/Coral_Model/models/coral_mortality_model.RData")

bayes_R2(coral_growth_model)
bayes_R2(coral_recruitment_model)
bayes_R2(coral_mortality_model)
