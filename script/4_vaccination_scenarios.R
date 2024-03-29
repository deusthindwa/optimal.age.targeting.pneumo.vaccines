# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

# function to assign initial VE values
initial_VE <- function(age, serogroup, age_dep = FALSE){
  
  # age-dependent initial vaccine efficacy at the time of vaccination may be superseded
  dplyr::case_when(
    serogroup == "PCV13" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE      ~ 0),
  
    serogroup == "PCV15" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE      ~ 0),
    
    serogroup == "PCV20" ~ 
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE      ~ 0),
      
    serogroup == "PPV23" ~
      dplyr::case_when(age_dep == FALSE ~ 1,
                       age >= 75 ~ 0.30/0.54,
                       age >= 65 ~ 0.36/0.54,
                       age >= 55 ~ 0.54/0.54,
                       TRUE      ~ 0),
      TRUE    ~ NA_real_)
}

# create scenarios table based on initial VE values, assumptions, vaccine type and age
  scenarios <- crossing(Study.waning = c("Fast", "Slow"),
                        age_dep = c(FALSE, TRUE),
                        serogroup = c("PPV23", "PCV13", "PCV15", "PCV20")) %>%
  mutate(Study.waning = case_when(Study.waning == "Slow" ~ "Djennad et al. (2018)",
                                  Study.waning == "Fast" ~ "Andrews et al. (2012)",
                                  TRUE                   ~ NA_character_)) %>%

crossing(Vac.age = seq(55, 85, by = 5),
         age     = seq(55, 90, by = 1)) %>%
  mutate(t     = age - Vac.age) %>%
  mutate(scale_initial = initial_VE(Vac.age, serogroup, age_dep),
         scale_initial = ifelse(t < 0, 0, scale_initial))

# ensure we scale the PCV ones properly. may need to merge with df_from_study_mc to get the waning right
# basically we need to make sure the PCV VE at t = 0-4 is 0.75 and then decays accordingly

    scenarios %<>%
    mutate(scale_initial = case_when(
    Study.waning == "Andrews et al. (2012)" & 
      grepl(pattern = 'PCV', x = serogroup) &
      t > 4 & t <= 6 ~ scale_initial * 65/100, # fast waning after 5 years for the next 2 years
    
    Study.waning == "Andrews et al. (2012)" & 
      grepl(pattern = 'PCV', x = serogroup) &
      t > 6 ~ scale_initial * 48/100, # fast waning from 7 years onward
    
    Study.waning == "Djennad et al. (2018)" & 
      grepl(pattern = 'PCV', x = serogroup) &
      t > 4 & t <= 6 ~ scale_initial * 83/100, # slow waning after 5 years for the next 2 years
    
    Study.waning == "Djennad et al. (2018)" & 
      grepl(pattern = 'PCV', x = serogroup) &
      t > 6 ~ scale_initial * 68/100, # slow waning from 7 years onward

     TRUE                                      ~ scale_initial))

# merge vaccination scenarios and waning VE (df_from_study_mc) separately for PPV and PCVs to avoid errors
VE_by_Vac.age <-
  rbind( 
    filter(scenarios, serogroup == "PPV23") %>%
      inner_join(select(df_from_study_mc, Study.waning, serogroup, sim, t, fit)), # get initial VE for PPV
      
      filter(scenarios, serogroup != "PPV23") %>%
      inner_join(select(filter(df_from_study_mc, serogroup != "PPV23"), serogroup, sim, t, fit))) %>% # get initial VE for PCVs
  
  mutate(VE = fit/100 * scale_initial)

# plot waning VE against time since vaccination
VE_time <- 
  VE_by_Vac.age %>%
  filter(serogroup == "PCV13" | serogroup == "PPV23") %>%
  mutate(Study.waning = if_else(Study.waning == "Andrews et al. (2012)", "Fast waning", "Slow waning"),
         serogroup = if_else(serogroup == "PCV13", "PCVs", "PPV23")) %>%
  filter(Vac.age %in% c(55, 65, 75)) %>%
  filter(!is.na(VE)) %>%
  select(-scale_initial, -fit) %>%
  nest(data = c(sim, VE)) %>%
  mutate(Q = map(data, ~quantile(.x$VE, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  
  ggplot(data = ., aes(x = t, y = `50%`)) +
  xlab("Years since vaccination") +
  ylab("Vaccine efficacy/effectiveness") +
  scale_y_continuous(limits = c(NA, 1.05), labels = scales::percent_format(accuracy = 1)) +
  RcmdrPlugin.KMggplot2::geom_stepribbon(aes(fill = age_dep, ymin = `2.5%`, ymax = `97.5%`), color = NA, alpha = 0.2) +
  geom_step(aes(color = age_dep)) + 
  facet_nested(Vac.age ~ Study.waning  + serogroup, labeller = labeller(Vac.age = function(x){paste("Vacc. age:", x)})) +
  scale_fill_brewer(palette = 'Set1', name = 'Age dependent vaccine efficacy/effectiveness') +
  scale_color_brewer(palette = 'Set1', name = 'Age dependent vaccine efficacy/effectiveness') +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_continuous(limit = c(0, 10), breaks = ~pretty.default(., n=3)) +
  theme(legend.position = 'bottom', panel.grid.minor = element_blank()) 

ggsave(filename = "output/S5_Fig_vaccine_efficacy_time.png", 
       plot = VE_time,
       width = 8, height = 7, units = "in", dpi = 300)


# combine estimated VE and popn cases and demography (based on fitted incidence and smooth population)
VE_by_Vac.ageX <-
  VE_by_Vac.age %>%
  inner_join(pop_cases2 %>% 
               select(serogroup, sim, Vac.age, country, cases)) %>%
  dplyr::mutate(Impact = VE*cases)
