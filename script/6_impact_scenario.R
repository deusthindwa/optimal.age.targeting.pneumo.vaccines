# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# Calculate the impact of an intervention targeting 65 year-old cohort for life time
# Assume coverage target of 70%

coverage <- 0.7

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age >= 65) %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

impact1 <- A65 %>%
  dplyr::filter(Vac.age == 65) %>%
  dplyr::group_by(country, Waning, sim) %>%
  dplyr::mutate(value  = coverage*Impact/sum(cases),
                Waning = sub(pattern     = "\\swaning", 
                             replacement = "", 
                             x           = Waning)) %>% 
  dplyr::select(-delay, -Impact, -cases) %>%
  tidyr::nest(data = c(sim, value)) %>%
  dplyr::mutate(Q = purrr::map(data, 
                               ~quantile(.x$value,
                                         probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")),
                   .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country,
                serogroup,
                Waning,
                `Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  mutate(Waning = if_else(Waning == "Fast \n(5 years' delay)", "Fast (5 years' delay)", Waning)) %>%
  arrange(-desc(country))


readr::write_csv(x    = impact1, 
                 path = here("output", "impact_65y_70pc.csv"))


#===============================================================================

# calculate the impact on preventable burden comparing vaccine products
# a scenario of 65 year-olds in England assuming coverage target of 70%

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65) %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prevantable_burden_65y_vax <- A65 %>%
  mutate(Waning = if_else(Waning == "Fast waning \n(5 years' delay)", "Fast waning",
                          if_else(Waning == "No waning", "Slow waning", Waning))) %>%
  dplyr::group_by(country, sim, age_dep, Waning) %>%
  dplyr::mutate(value  = Impact/sum(Impact)) %>%
  dplyr::select(-delay, -Impact, -cases) %>%
  tidyr::nest(data = c(sim, value)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$value, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")),
                   .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country,
                serogroup,
                Waning,
                #`Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  distinct(Impact)

readr::write_csv(x    = prevantable_burden_65y_vax, 
                 path = here("output", "prevantable_burden_65y_vax.csv"))



#===============================================================================

# calculate preventable burden comparing vaccine products
# a scenario of 65 year-olds in England assuming coverage target of 70%

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(Vac.age == 65) %>% 
  dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prevantable_cases_65y_vax <- A65 %>%
  mutate(Waning = if_else(Waning == "Fast waning \n(5 years' delay)", "Fast waning",
                          if_else(Waning == "No waning", "Slow waning", Waning))) %>%
  dplyr::group_by(country, sim, age_dep, Waning) %>%
  dplyr::mutate(value  = cases/sum(cases)) %>%
  dplyr::select(-delay, -Impact, -cases) %>%
  tidyr::nest(data = c(sim, value)) %>%
  dplyr::mutate(Q = purrr::map(data, ~quantile(.x$value, probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")),
                   .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(country,
                serogroup,
                Waning,
                #`Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
  arrange(-desc(country)) %>%
  distinct(Impact)

readr::write_csv(x    = prevantable_cases_65y_vax, 
                 path = here("output", "prevantable_cases_65y_vax.csv"))


