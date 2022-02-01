# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

#add empty row in VE table
VE_table <- tibble::add_row(VE_table, Study = "None", rate = 0, `Half-life` = Inf)

# function to assign initial VE values. make a list check it twice
initial_VE <- function(age, serogroup, age_dep = FALSE){
  
  # age-dependent vaccine efficacy at time of vaccination may be superseded
  dplyr::case_when(
    serogroup == "PCV13" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE    ~ NA_real_),
    
    # age-dependent vaccine efficacy at time of vaccination may be superseded
    serogroup == "PCV20" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE    ~ NA_real_),
    
    # if no age dependency, we need to just use the value from relevant study, which we can handle outside this
    serogroup == "PPV23" ~
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE    ~ NA_real_),
    TRUE    ~ NA_real_)
}

df_from_study_ <- distinct(df_from_study, Study, VE, rate, sim)

# create scenarios table based on initial VE values, assumptions, vaccine type and age
scenarios <- list(`PPV23` = crossing(delay = 0,
                                     Study.waning = c("Andrews et al. (2012)",
                                                      "Djennad et al. (2018)"),
                                     age_dep = c(FALSE, TRUE),
                                     serogroup = "PPV23"),
                  `PCV`   = crossing(delay = c(0,5),
                                     age_dep = c(FALSE, TRUE),
                                     serogroup = c("PCV13",
                                                   "PCV20")) %>%
                    mutate(Study.waning = ifelse(delay == 0,
                                                 "None",
                                                 "Andrews et al. (2012)"))) %>%
  bind_rows %>%
  mutate(Study.VE = ifelse(Study.waning == "Djennad et al. (2018)",
                           "Djennad et al. (2018)",
                           "Andrews et al. (2012)")) %>%
  crossing(Vac.age = seq(55, 85, by = 5),
           age     = seq(55, 85, by = 1)) %>%
  mutate(t     = age - Vac.age) %>%
  mutate(t_eff = pmax(0, t - delay)) %>%
  mutate(scale_initial = initial_VE(Vac.age, serogroup, age_dep),
         scale_initial = ifelse(t < 0, 0, scale_initial)) 

# we want the curve to be at VE if age >= vac.age + delay. when delay > 0, we want to subtract delay off
VE_by_Vac.age <- 
  scenarios %>%
  inner_join(select(df_from_study, Study.VE = Study, fit, sim)) %>% # get initial waning
  mutate(VE = fit * scale_initial/100) %>%
  dplyr::mutate(Impact = VE*cases) 

# need to deal with agey_since in terms of dealing with 5 year delay
