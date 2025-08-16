
Coral Model Metadata
David Liddle; david.liddle@duke.edu

Raw-ish Data Metadata

| Supplied Name | Supplied Description                                | Units          |
| ------------- | --------------------------------------------------- | -------------- |
| GenSiteID     | Unique identifier for each colony at a given site   | unitless       |
| Site          | Field site where coral was observed                 | unitless       |
| Class         | Taxonomic grouping of coral (e.g., Pocillopora)     | unitless       |
| size_YEAR     | Colony planar area for each survey year             | cm²            |
| Fate_YEAR     | Colony outcome (alive, dead, missing) for each year | unitless       |
| Notes         | Field comments about colony                         | text           |



Processed Metadata

| Supplied Name             | Supplied Description                                | Units.         |
| ------------------------- | --------------------------------------------------- | -------------- |
| Site                      | Site name of coral colony                           | unitless       |
| GenSiteID                 | Unique colony identifier across sites               | unitless       |
| Class                     | Coral taxonomic class/group                         | unitless       |
| Year                      | Year of survey (2013–2019)                          | year           |
| CoralSize                 | Colony planar size                                  | cm²            |
| log_CoralSize             | Natural log of CoralSize (with offset for 0)        | log(cm²)       |
| Recruitment               | Indicator: 1 = colony recruited this year           | binary         |
| Mortality                 | Indicator: 1 = colony died this year                | binary         |
| GrowthRate                | Log ratio of CoralSize between years                | log ratio      |
| Next_Year_Mortality       | Mortality in the following year                     | binary         |
| Max_DHW                   | Maximum annual Degree Heating Weeks                 | °C-weeks       |
| Max_DHW_Lag1              | Max DHW one year prior                              | °C-weeks       |
| Max_DHW_Lag2              | Max DHW two years prior                             | °C-weeks       |
| Mortality_Lag1            | Colony mortality in previous year                   | binary         |
| FirstYear                 | Indicator: 1 = first survey year (2013)             | binary         |
| first_recruitment_year    | Year colony first recruited                         | year           |
| first_mortality_year      | Year colony first died                              | year           |
| SiteMortality_Prop        | Proportion of colonies that died at site-year level | proportion     |
| SiteMortality_Prop_Lag1   | Lagged site mortality proportion                    | proportion     |


Model Prior Justifications:

1. Growth Model

Formula:
GrowthRate ~ Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + log_CoralSize + (1 | Class)
Family: Student-t

Intercept: normal(0, 0.75)
Mean log-growth rate = 0.0652, SD = 0.735.
Center  prior near zero (close to 0.0652) and set the SD equal to the observed spread (0.75)

Slopes: normal(0, 0.3)
Predictors (DHW, lagged DHW, log size) are in natural units with observed ranges far smaller than the variation in growth rate.
A prior SD of 0.3 is less than half the observed SD of growth, constraining effect sizes so that a 1-unit change in predictors rarely shifts log-growth by more than ±0.5.

Random intercept SD: exponential(1)
Allows for ~95% of species-level intercepts to vary by ±2 log units while favoring smaller differences. Based on observed interspecies spread in growth rates.

ν parameter: gamma(5, 1)
Distribution shows occasional extreme growth/shrinkage values. gamma(5,1) gives the model flexibility to expect a few extreme cases without letting them dominate the results.



2. Recruitment Model

Formula:
Recruitment ~ Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + SiteMortality_Prop_Lag1 + (1 | Class) + (1 | GenSiteID)
Family: Bernoulli with logit link

Intercept: normal(-4, 1)
Lifetime proportion of colonies that ever recruit = 0.00272 (~0.27%).
On the logit scale, -4 corresponds to ~1.8% annual recruitment, which is slightly higher than the lifetime rate to reflect that the lifetime measure dilutes annual probability.

Slopes: normal(0, 0.5)
Weakly informative to allow for moderate predictor effects without permitting extreme recruitment probabilities without strong data support.

Random intercept SD: exponential(1)
Recruitment rates can differ between species and even between individual colonies or sites.
This prior favors smaller differences but allows the model to recognize if some groups are more prone to recruitment.



3. Mortality Model

Formula:
Next_Year_Mortality ~ log_CoralSize + Max_DHW + Max_DHW_Lag1 + Max_DHW_Lag2 + (1 | Class)
Family: Bernoulli

Intercept: normal(-0.85, 0.5)
Observed overall annual mortality rate = 0.304 (~30%).
Logit(-0.85) ≈ 30% probability, 0.5 SD adds modest uncertainty, giving the model flexibility to shift.

log_CoralSize slope: normal(-0.4, 0.3)
Mortality breakdown shows smallest colonies have ~48–55% mortality while largest colonies drop to ~8–16%.
This prior reflects our expectation that mortality risk goes down as colony size goes up, without locking in the exact relationship.

Max_DHW & lag slopes: normal(0.02, 0.02)
Small positive association between heat stress and mortality.
Narrow SD reflects strong prior belief (from literature + data) in a consistently small positive effect per unit DHW.

Random intercept SD: cauchy(0, 1)
Some coral groups are much more fragile than others. This prior gives the model the freedom to capture large differences between groups while still assuming most differences are moderate.
