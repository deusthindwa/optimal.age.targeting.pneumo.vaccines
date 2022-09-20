# written by Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021

#===============================================================================

# Calculation of the proportion of prevented cases relative to unvaccinated
impact_total <- 
  ipd_mc %>%
  mutate(cases = map(.x = mc, .f = ~group_by(.x, sim) %>%
                       crossing(Vac.age = seq(-55, 85, by = 5)) %>%
                       filter(agey >= Vac.age) %>%
                       group_by(Vac.age, sim) %>%
                       summarise(cases = sum(fit)))) %>%
  select(-data, -model, -mc) %>%
  unnest(cases) %>%
  inner_join(VE_impact_by_age) %>%
  mutate(rel_impact = Impact/cases*100) %>%
  group_by_at(.vars = vars(-c(sim, cases, Impact, rel_impact))) %>%
  nest %>%
  mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

#===============================================================================

# Calculation of the number of individuals needed to vaccinate to prevent a case
impact_case <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = ntotal/Impact) %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

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
  filter(Vac.age == 65 & Waning == "Fast waning" & age_dep == FALSE)
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
impact_case %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == FALSE)
readr::write_csv(x = Table_S4, file = here("output", "Table_S4_scenario.csv"))
