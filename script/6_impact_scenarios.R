# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

#===============================================================================

#get full incidence uncertainty dataset
ipd_mc3 <- 
  unnest(ipd_mc, mc) %>%
  select(serogroup, country, agey, sim, fit) %>%
  arrange(serogroup, country, agey, sim) %>%
  filter(serogroup == "All")

#compute full case uncertainty dataset
pop_cases3 <- 
  dplyr::left_join(ipd_mc3, pop_country_df, by = c("agey", "country")) %>% 
  dplyr::mutate(cases  = fit/scale*ntotal, Vac.age = agey) %>%
  rename("age" = "agey")

# Calculation of the proportion of prevented VT IPD of all IPD in 55+y
impact_total <-
VE_impact_by_age %>%
  inner_join(filter(pop_cases3) %>% 
               select(country, Vac.age, sim, cases)) %>% 
  dplyr::group_by(Waning,
                  age_dep,
                  sim,
                  country) %>%
  mutate(tot_cases = sum(cases), 
         reImpact = Impact/tot_cases*100) %>%
  select(sim, country, serogroup, age_dep, Waning, Vac.age, reImpact) %>%
  nest(data = c(sim, reImpact)) %>%
  mutate(Q = map(data, ~quantile(x = .x$reImpact, probs = c(0.025, 0.5, 0.975)))) %>%
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
Table_S1 <-
  impact_total %>%
  filter(serogroup == "PPV23" & Waning == "Fast waning" & age_dep == FALSE) 
readr::write_csv(x = Table_S1, file = here("output", "Table_S1_scenario.csv"))

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
