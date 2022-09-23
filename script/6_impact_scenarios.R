# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

#===============================================================================

# Calculation of the proportion of prevented VT IPD of all IPD in 55+y
impact_total <-
VE_impact_by_age %>%
  inner_join(filter(pop_cases, serogroup == "All") %>% 
               select(country, Vac.age, cases)) %>%
  dplyr::group_by(age_dep,
                  country,
                  Waning,
                  sim) %>%
  mutate(tot_cases = sum(cases),
         reImpact = Impact/tot_cases*100) %>%
select(sim, country, serogroup, age_dep, Waning, Vac.age, reImpact) %>%
  nest(data = c(sim, reImpact)) %>%
  mutate(Q = map(data, ~quantile(x     = .x$reImpact,
                                 probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

#===============================================================================

# Calculation of the number of individuals needed to vaccinate to prevent a case
efficiency <- 
  dplyr::select(pop_country_df, country, agey, ntotal, N) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  dplyr::group_by(serogroup,
                  Waning,
                  age_dep,
                  sim,
                  Vac.age,
                  country) %>%
  summarise(eff = ntotal/Impact) %>%
  ungroup() %>%
  nest(data = c(sim, eff)) %>%
  mutate(Q = map(data, ~quantile(x     = .x$eff,
                                 probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

#===============================================================================

# (Table S1) A scenario of PPV23 use under fast waning efficacy/effectiveness. 
# Comparing preventable cases between different vaccination age cohorts (e.g., 55 vs 70 years cohorts) relative to unvaccinated.
Table_S1a <-
  impact_total %>%
  filter(serogroup == "PCV13" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1a, file = here("output", "Table_S1_scenario_pcv13.csv"))

Table_S1b <-
  impact_total %>%
  filter(serogroup == "PCV15" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1b, file = here("output", "Table_S1_scenario_pcv15.csv"))

Table_S1c <-
  impact_total %>%
  filter(serogroup == "PCV20" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1c, file = here("output", "Table_S1_scenario_pcv20.csv"))

Table_S1d <-
  impact_total %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1d, file = here("output", "Table_S1_scenario_pp23.csv"))

#===============================================================================

# (Table S2) A scenario of vaccinating 65 years old cohort under fast waning efficacy/effectiveness. 
# Comparing preventable cases between use of different vaccines (e.g., PCV20 vs PPV23) relative to unvaccinated.
Table_S2 <-
  impact_total %>%
  filter((Vac.age == 65 | Vac.age == 55) & Waning == "Fast waning" & age_dep == FALSE)
readr::write_csv(x = Table_S2, file = here("output", "Table_S2_scenario.csv"))

#===============================================================================

# (Table S3) A scenario of PPV23 use in 65 years old cohort. 
# Comparing preventable cases between fast and slow waning efficacy/effectiveness relative to unvaccinated.
Table_S3 <-
  impact_total %>%
  filter(Vac.age == 65 & age_dep == FALSE)
readr::write_csv(x = Table_S3, file = here("output", "Table_S3_scenario.csv"))

#===============================================================================

# (Table S4) A scenario of PPV23 use under fast waning efficacy/effectiveness. 
# Comparing number of individuals needed to vaccinate to prevent a case (e.g., in 55 vs 85 years old)
Table_S4 <-
  efficiency %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == FALSE)
readr::write_csv(x = Table_S4, file = here("output", "Table_S4_scenario.csv"))
