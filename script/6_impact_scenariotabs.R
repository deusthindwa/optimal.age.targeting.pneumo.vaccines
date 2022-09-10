# written by Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021


#===============================================================================

# (Table S1) comparing proportion of cases prevented between 55y vs 70y among all cases averted

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

Table_S1 <- VE_impact_by_age %>%
  dplyr::filter(!is.na(country) & age_dep == FALSE & Waning == "Fast waning", serogroup == "PPV23") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

Table_S1 <- 
  Table_S1 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact, -cases) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Vac.age, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  filter(`Age dep.` == FALSE)

readr::write_csv(x = Table_S1, path = here("output", "Table_S1_scenario.csv"))

#===============================================================================

# (Table S2) comparing proportion of cases prevented by PCV20 vs PPV23 among all cases averted

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

Table_S2 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65 & !is.na(country) & age_dep == FALSE & Waning == "Fast waning") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

Table_S2 <- 
  Table_S2 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact, -cases) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  filter(`Age dep.` == FALSE)

readr::write_csv(x = Table_S2, path = here("output", "Table_S2_scenario.csv"))

#===============================================================================

# (Table S3) comparing proportion of cases prevented under fast vs slow waning VE among all cases averted

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

Table_S3 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 80 & !is.na(country) & age_dep == FALSE, serogroup == "PPV23") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

Table_S3 <- 
  Table_S3 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact, -cases) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  filter(`Age dep.` == FALSE)

readr::write_csv(x = Table_S3, path = here("output", "Table_S3_scenario.csv"))

#===============================================================================

# (Table S4) comparing proportion of cases prevented by 55y vs 85y waning VE among all cases averted

Table_S4 <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = Impact*scale/ntotal) %>%
  dplyr::filter(!is.na(country) & age_dep == FALSE & Waning == "Fast waning", serogroup == "PPV23")

Table_S4 <- 
  Table_S4 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Vac.age, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  filter(`Age dep.` == FALSE)

readr::write_csv(x = Table_S4, path = here("output", "Table_S4_scenario.csv"))

#===============================================================================

# (Table S5) comparing proportion of cases prevented by 55y vs 65y vs 75y waning VE among all averted cases (age dependent)

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

Table_S5 <- VE_impact_by_age %>%
  dplyr::filter(!is.na(country) & age_dep == TRUE & Waning == "Fast waning", serogroup == "PPV23") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

Table_S5 <- 
  Table_S5 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact, -cases) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Vac.age, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) 

readr::write_csv(x = Table_S5, path = here("output", "Table_S5_scenario.csv"))

#===============================================================================

# (Table S6) comparing proportion of cases prevented per vaccinee by 55y vs 65y vs 75y waning VE among all averted cases (age dependent)

Table_S6 <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = Impact*scale/ntotal) %>%
  dplyr::filter(!is.na(country) & age_dep == TRUE & Waning == "Fast waning", serogroup == "PPV23")

Table_S6 <- 
  Table_S6 %>%
  dplyr::ungroup() %>%
  dplyr::group_by(country, sim) %>%
  dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
  dplyr::select(-Impact) %>%
  tidyr::nest(data = c(sim, rel_impact)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country, serogroup, Vac.age, Waning, `Age dep.` = age_dep, contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country))

readr::write_csv(x = Table_S6, path = here("output", "Table_S6_scenario.csv"))
