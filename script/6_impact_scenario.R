# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021

#===============================================================================

# Optimal age for vaccination (cohort vaccination) 55 vs 70 years old

cases <- dplyr::inner_join(unnest(ipd_mc, mc), pop_country_df) %>%
    select(-model, -data) %>%
    dplyr::filter(serogroup != "All") %>%
    dplyr::mutate(cases = fit*ntotal/scale, Vac.age = agey)


scenarios <- list(
    `55-70`      = list(
        filter        = data.frame(age_dep   = FALSE,
                                   Waning    = "Fast waning",
                                   serogroup = "PPV23"),
        grouping_vars = c("country", "age_dep", "sim",
                          "Study.waning",
                          "serogroup",
                          "Waning")),
    `65`         = list(
        filter        = data.frame(Vac.age   = 65,
                                   age_dep   = FALSE,
                                   Waning    = "Fast waning",
                                   serogroup = c("PCV20", "PPV23")),
        grouping_vars = c("country",
                          "sim",
                          "age_dep")),
    `65 waning`  = list(
        filter        = data.frame(Vac.age   = 65,
                                   age_dep   = FALSE, 
                                   serogroup = "PPV23"),
        grouping_vars = c("country", "age_dep", "sim",
                          #"Study.waning",
                          "serogroup"#,
                          #"Waning"
        )),
    `80`         = list(
        filter        = data.frame(age_dep   = FALSE,
                                   Waning    = "Fast waning", 
                                   serogroup = "PPV23",
                                   Vac.age   = c(55, 85)),
        grouping_vars = c("country", "age_dep", "sim",
                          "Study.waning",
                          "serogroup",
                          "Waning")),
    `65, 75, 85` = list(
        filter        = data.frame(age_dep   = TRUE,
                                   Waning    = "Fast waning",
                                   serogroup = "PPV23"),
        grouping_vars = c("country", "age_dep", "sim",
                          "Study.waning",
                          "serogroup",
                          "Waning"))
)


# turn into a function

make_averted_cases <- function(x, grouping_vars = c("country", "age_dep", "sim",
                                                    "Study.waning",
                                                    "serogroup",
                                                    "Waning")){
    
    dplyr::ungroup(x) %>%
        dplyr::group_by_at(.vars = vars(any_of(grouping_vars))) %>%
        dplyr::mutate(Impact_total = sum(Impact)) %>%
        dplyr::mutate(rel_impact  = Impact/sum(Impact)) %>%
        dplyr::select(-any_of(c("Impact", "cases", "Impact_total"))) %>%
        ungroup %>%
        tidyr::nest(data = c(sim, rel_impact)) %>%
        dplyr::mutate(Q = purrr::map(data, ~quantile(x     = .x$rel_impact,
                                                     probs = c(0.025, 0.5, 0.975)))) %>%
        tidyr::unnest_wider(Q) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate_at(.vars =  dplyr::vars(contains("%")), 
                         .funs = ~scales::percent(., 0.1)) %>%
        dplyr::select(any_of(grouping_vars), contains("%")) %>%
        dplyr::rename(`Age dep.` = age_dep) %>%
        dplyr::group_by_at(.vars = dplyr::vars(-contains("%"))) %>%
        dplyr::transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`)) %>%
        arrange(-desc(country))
}

impact_scenarios <- scenarios %>% 
    {map(.x = ., ~inner_join(VE_impact_by_age, .x$filter) %>%
             na.omit %>%
             dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases)) %>%
             make_averted_cases(grouping_vars = .x$grouping_vars))} 

# this probably needs a better naming scheme
impact_scenarios %>%
    {map2(.x = ., .y = names(.),
          ~write_csv(x = .x,
                     path = sprintf("output/scenario %s.csv", .y)))}

# 55-70
# readr::write_csv(x    = prop_averted_cases_5570y_vax, 
#                  path = "output/Table_S2_prop_averted_cases_55_70y.csv")

#===============================================================================

# Optimal age for vaccination (cohort vaccination) PCV20 vs PPV23
# 65
# readr::write_csv(x    = prop_averted_cases_65y_vax, 
#                  path = "output/Table_S3_prop_averted_cases_65y.csv")

#===============================================================================

# Optimal age for vaccination (cohort vaccination) fast vs slow waning
# 65 waning
# readr::write_csv(x    = prop_averted_cases_65y_wane, 
#                  path = "output/Table_S4_prop_averted_cases_wane.csv")

#===============================================================================

# Optimal age for vaccination (individual vaccination) 55 vs 80 years old
A80 <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
    dplyr::rename(Vac.age = agey) %>% 
    dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
    mutate(Impact = Impact*scale/ntotal) %>%
    dplyr::filter(!is.na(country) ,
                  age_dep   == FALSE,
                  Waning    == "Fast waning",
                  serogroup == "PPV23",
                  Vac.age   %in% c(55, 80)) %>%
    select(-ntotal)

prop_averted_cases_80y_vax <- make_averted_cases(A80,
                                                 grouping_vars = c("country", "age_dep", "sim",
                                                                   "Study.waning",
                                                                   "serogroup",
                                                                   "Waning") )

readr::write_csv(x    = prop_averted_cases_80y_vax, 
                 path = "output/Table_S5_prop_averted_cases_80y.csv")

#===============================================================================

# Optimal age for vaccination (cohort vaccination + age dependent) 55 vs 65 vs 75 years old

A556575 <- VE_impact_by_age %>%
    dplyr::filter(!is.na(country) & age_dep == TRUE & Waning == "Fast waning", serogroup == "PPV23") %>% 
    dplyr::left_join(dplyr::select(cases, serogroup, Vac.age, country, sim, cases))

prop_averted_cases_556575y_vax <- 
    A556575 %>%
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

readr::write_csv(x    = prop_averted_cases_556575y_vax, 
                 path = "output/Table_S6_prop_averted_cases_55_65_75y.csv")

#===============================================================================

# Optimal age for vaccination (individual vaccination + age depenndent) 55 vs 65 vs 75 years old
A556575 <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
    dplyr::rename(Vac.age = agey) %>% 
    dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
    mutate(Impact = Impact*scale/ntotal) %>%
    dplyr::filter(!is.na(country) & age_dep == TRUE & Waning == "Fast waning", serogroup == "PPV23")

prop_averted_cases_556575y_vax <- 
    A556575 %>%
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

readr::write_csv(x    = prop_averted_cases_556575y_vax, 
                 path = "output/Table_S7_prop_averted_cases_55_65_57y_vax.csv")
