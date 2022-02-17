# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# function to assign initial VE values. make a list check it twice
initial_VE <- function(age, serogroup, age_dep = FALSE){
  
  # age-dependent vaccine efficacy at time of vaccination may be superseded
  dplyr::case_when(
    serogroup == "PCV13" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE      ~ 0),
    
    # age-dependent vaccine efficacy at time of vaccination may be superseded
    serogroup == "PCV15" ~ 
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
                       TRUE      ~ 0),
    TRUE    ~ NA_real_)
}

# df_from_study_ <- distinct(df_from_study, Study, VE, rate, sim)

# create scenarios table based on initial VE values, assumptions, vaccine type and age
scenarios <- crossing(Study.waning = c("Fast", "Slow"),
                      age_dep = c(FALSE, TRUE),
                      serogroup = c("PPV23", "PCV13", "PCV15", "PCV20")) %>%
  mutate(Study.waning = case_when(Study.waning == "Slow" ~ "Djennad et al. (2018)",
                                  Study.waning == "Fast" ~ "Andrews et al. (2012)",
                                  TRUE                   ~ NA_character_)) %>%
  
  crossing(Vac.age = seq(55, 85, by = 5),
           age     = seq(55, 85, by = 1)) %>%
  mutate(t     = age - Vac.age) %>%
  mutate(scale_initial = initial_VE(Vac.age, serogroup, age_dep),
         scale_initial = ifelse(t < 0, 0, scale_initial))

# ensure we scale the PCV ones properly. may need to merge in the info from df 
# to get the waning right
# basically we need to make sure the PCV VE at t=0 is 0.75 and then decay according
# to either djennad or andrews' scalings.

scenarios %<>%
  mutate(scale_initial = case_when(
    Study.waning == "Andrews et al. (2012)"   &
      grepl(pattern = 'PCV', x = serogroup) &
      t >= 5 ~ scale_initial * 15/31.8,
    Study.waning == "Djennad et al. (2018)"   & 
      grepl(pattern = 'PCV', x = serogroup) &
      t >= 5 ~ scale_initial * 23/36.8,
    TRUE                                      ~ scale_initial))


VE_by_Vac.age <- 
  scenarios %>%
  inner_join(select(df_from_study, t,
                    Study.waning = Study, fit, sim)) %>% # get initial waning
  mutate(VE = fit/100 * scale_initial) 

VE_by_Vac.age <- pop_cases %>%
  select(serogroup, country, age = agey, cases, Vac.age)  %>%
  right_join(VE_by_Vac.age) %>%
  dplyr::mutate(Impact = VE*cases)

VE_time <- 
  VE_by_Vac.age %>%
  filter(Vac.age %in% c(55, 65, 75)) %>%
  select(-scale_initial, -fit) %>%
  nest(data = c(sim, VE, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$VE, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  ggplot(data = ., aes(x = t, y = `50%`)) +
  xlab("Years since vaccination") +
  ylab("Vaccine efficacy") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  RcmdrPlugin.KMggplot2::geom_stepribbon(aes(fill = age_dep,
                                             ymin = `2.5%`,
                                             ymax = `97.5%`),
                                         color = NA, alpha = 0.2) +
  geom_step(aes(color = age_dep)) + 
  facet_nested(Vac.age ~ Study.waning  + serogroup, 
               labeller = labeller(Vac.age = function(x){paste("Vacc. age:", x)})) +
  scale_fill_brewer(palette = 'Set1', name = 'Age dependent vaccine effectiveness') +
  scale_color_brewer(palette = 'Set1', name = 'Age dependent vaccine effectiveness') +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_continuous(limits = c(0,10), breaks = ~pretty.default(., n=3)) +
  theme(legend.position = 'bottom', panel.grid.minor = element_blank()) 

ggsave(filename = "output/Fig_vaccine_effectiveness_time.png", 
       plot = VE_time,
       width = 14, height = 8, units = "in", dpi = 300)


