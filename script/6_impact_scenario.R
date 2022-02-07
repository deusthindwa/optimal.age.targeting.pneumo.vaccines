# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# Calculate the impact of an intervention targeting 65 year olds
# Assume coverage target of 70%

coverage <- 0.7

Cases <- dplyr::inner_join(unnest(ipd_mc, mc), 
                           pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  #dplyr::filter(serogroup == "PPV23",
  #              country   == "England",
  #              Vac.age >= 65) %>% 
  
  dplyr::filter(country   == "England",
                Vac.age >= 65) %>% 
  
  dplyr::left_join(
    dplyr::select(Cases,
                  serogroup,
                  Vac.age, 
                  country, 
                  sim,
                  cases))

impact_65y_70pc <- A65 %>%
  dplyr::filter(Vac.age == 65) %>%
  dplyr::group_by(Waning, sim) %>%
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
  dplyr::select(serogroup,
                Waning,
                `Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`))


readr::write_csv(x    = impact_65y_70pc, 
                 path = here("output", "impact_65y_70pc.csv"))


#===============================================================================

# Calculate the impact of an intervention targeting 65 year olds
# Assume coverage target of 70%

coverage <- 0.7

Cases <- dplyr::inner_join(unnest(ipd_mc, mc), 
                           pop_country_df) %>%
  dplyr::filter(serogroup != "All") %>%
  dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)

A65 <- VE_impact_by_age %>%
  dplyr::filter(country   == "England",
                Vac.age >= 55) %>% 
  
  dplyr::left_join(
    dplyr::select(Cases,
                  serogroup,
                  Vac.age, 
                  country, 
                  sim,
                  cases))

A65 %>%
  dplyr::filter(Vac.age == 65) %>%
  dplyr::group_by(Waning, sim) %>%
  dplyr::mutate(value  = cases/sum(cases),
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
  dplyr::select(serogroup,
                Waning,
                `Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(PropVT = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`))


readr::write_csv(x    = impact_65ye_70pc, 
                 path = here("output", "impact_65y_70pc.csv"))
