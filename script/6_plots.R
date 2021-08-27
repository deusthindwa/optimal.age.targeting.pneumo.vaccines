# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021


# compute maximum vaccine impact
 VE_impact_max <- VE_impact_by_age %>%
   dplyr::group_by_at(.vars = dplyr::vars(-c(Vac.age, Impact))) %>%
   dplyr::filter(Impact == max(Impact))

# 65y old programme impact (%) at 70% coverage
coverage <- 0.7

A65 <- impact_by_age_to_plot %>%
  dplyr::filter(serogroup == "PPV23",
                Country   == "England/Wales",
                Vac.age >= 65) %>% 
  dplyr::inner_join(
    dplyr::select(Cases,
                  serogroup,
                  Vac.age, 
                  Country, 
                  sim,
                  cases))

impact_65y_70pc <-
  dplyr::group_by(A65, Waning, sim) %>%
  dplyr::mutate(value  = coverage*Impact/sum(cases),
                Waning = sub(pattern     = "\\swaning", 
                             replacement = "", 
                             x           = Waning)) %>% 
  dplyr::filter(Vac.age == 65) %>%
  dplyr::select(-delay, -Impact, -cases) %>%
  tidyr::nest(data = c(sim, value)) %>%
  dplyr::mutate(Q = purrr::map(data, 
                               ~quantile(.x$value,
                                         probs = c(0.025, 0.5, 0.975)))) %>%
  tidyr::unnest_wider(Q) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate_at(.vars =  dplyr::vars(contains("%")),
                   .funs = ~scales::percent(., 0.1)) %>%
  dplyr::select(Waning,
                `Age dep.` = age_dep,
                contains("%")) %>%
  dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
  dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`))


readr::write_csv(x    = impact_65y_70pc, 
                 path = here("output", "impact_65y_70pc.csv"))

