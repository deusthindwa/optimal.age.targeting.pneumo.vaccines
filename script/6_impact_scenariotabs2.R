# written by Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021

#===============================================================================

# Calculation of preventable cases by vaccination age, vaccine type, age dependency, waning, and country.
impact_scene <- 
  ipd_mc %>%
  mutate(cases = map(.x = mc, .f = ~group_by(.x, sim) %>%
                       crossing(Vac.age = seq(55, 85, by = 5)) %>%
                       filter(agey >= Vac.age) %>%
                       group_by(Vac.age, sim) %>%
                       summarise(cases = sum(fit)))) %>%
  select(-data, -model, -mc) %>%
  unnest(cases) %>%
  inner_join(VE_impact_by_age) %>%
  mutate(rel_impact = Impact/cases) %>%
  group_by_at(.vars = vars(-c(sim, cases, Impact, rel_impact))) %>%
  nest %>%
  mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

#===============================================================================

# (Table S1) A scenario of PPV23 use under fast waning efficacy/effectiveness. 
# Comparing preventable cases between different vaccination age cohorts (e.g., 55 vs 70 years cohorts) relative to unvaccinated.
Table_S1 <-
impact_scene %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1, path = here("output", "Table_S1_scenario.csv"))

#===============================================================================

# (Table S2) A scenario of vaccinating 65 years old cohort under fast waning efficacy/effectiveness. 
# Comparing preventable cases between use of different vaccines (e.g., PCV20 vs PPV23) relative to unvaccinated.
Table_S2 <-
  impact_scene %>%
  filter(Vac.age == 65 & Waning == "Fast waning" & age_dep == FALSE)
readr::write_csv(x = Table_S2, path = here("output", "Table_S2_scenario.csv"))

#===============================================================================

# (Table S3) A scenario of PPV23 use in 65 years old cohort. 
# Comparing preventable cases between fast and slow waning efficacy/effectiveness relative to unvaccinated.
Table_S3 <-
  impact_scene %>%
  filter(Vac.age == 65 & age_dep == FALSE)
readr::write_csv(x = Table_S3, path = here("output", "Table_S3_scenario.csv"))

#===============================================================================

# (Table S4) A scenario of PPV23 use under fast waning efficacy/effectiveness. 
# Comparing preventable cases per vaccinee (e.g., individual aged 55 vs 85 years old) relative to unvaccinated.



readr::write_csv(x = Table_S4, path = here("output", "Table_S4_scenario.csv"))

#===============================================================================

# (Table S5) A scenario of PPV23 use under fast waning efficacy/effectiveness. (VE) and age-dependent initial VE. 
# Comparing preventable cases between vaccination age cohorts (e.g., 55 vs 65 vs 75 years old) relative to unvaccinated.
Table_S5 <-
  impact_scene %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == TRUE)
readr::write_csv(x = Table_S5, path = here("output", "Table_S5_scenario.csv"))

#===============================================================================

# (Table S6) A scenario of PPV23 use under fast waning efficacy/effectiveness (VE) and age-dependent initial VE. 
# Comparing preventable cases per vaccinee (e.g., individual aged 55 vs 65 vs 75 years old) relative to unvaccinated.


readr::write_csv(x = Table_S6, path = here("output", "Table_S6_scenario.csv"))
