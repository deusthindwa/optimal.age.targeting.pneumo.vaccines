# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021


#===============================================================================

# calculate the share of preventable cases by vaccinating with different vaccine products
# a scenario of 65 year old in England assuming coverage target of 100%

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65 & !is.na(country)) %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prop_averted_cases_65y_vax <- 
  A65 %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Waning = if_else(Waning == "Fast waning \n(5 years' delay)", "Fast waning",
                                 if_else(Waning == "No waning", "Slow waning", Waning))) %>%
  dplyr::group_by(country, sim, age_dep, Waning) %>%
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
  filter(`Age dep.` == TRUE)

readr::write_csv(x    = prop_averted_cases_65y_vax, 
                 path = here("output", "Table_S2_prop_averted_cases_65y_vax.csv"))

#===============================================================================

# Optimal age for vaccination (cohort vaccination) 55 vs 70 years old

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A5570 <- VE_impact_by_age %>%
  dplyr::filter(!is.na(country) & age_dep == FALSE & Waning == "Fast waning", serogroup == "PPV23") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prop_averted_cases_5570y_vax <- 
  A5570 %>%
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

readr::write_csv(x    = prop_averted_cases_5570y_vax, 
                 path = here("output", "Table_S3_prop_averted_cases_55_70y.csv"))

#===============================================================================

# Optimal age for vaccination (cohort vaccination) PCV20 vs PPV23

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65 & !is.na(country) & age_dep == FALSE & Waning == "Fast waning") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prop_averted_cases_65y_vax <- 
  A65 %>%
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

readr::write_csv(x    = prop_averted_cases_65y_vax, 
                 path = here("output", "Table_S4_prop_averted_cases_65y.csv"))

#===============================================================================

# Optimal age for vaccination (cohort vaccination) fast vs slow waning

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65wane <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65 & !is.na(country) & age_dep == FALSE, serogroup == "PPV23") %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prop_averted_cases_65y_wane <- 
  A65wane %>%
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

readr::write_csv(x    = prop_averted_cases_65y_wane, 
                 path = here("output", "Table_S5_prop_averted_cases_wane.csv"))
