# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021

# function for plotting vaccine impact
make_grid_plot <- function(x, ylab = NULL, percent = FALSE, ylim = c(0,NA)){
  
  
  p <- ggplot(x,
              aes(x = Vac.age, y = `50%`,
                  color = factor(age_dep),
                  group = interaction(Waning, age_dep, serogroup, country))) +
    geom_line() + 
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                    fill = factor(age_dep)), color = NA, alpha = 0.2) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                    fill = factor(age_dep)), color = NA, alpha = 0.2) +
    facet_nested(country ~ serogroup + Waning, scales = "free_y") +
    theme_bw() +
    scale_x_continuous(breaks = seq(60, 90, 10)) +
    theme(axis.text=element_text(size=10, color="black")) +
    xlab("Vaccination Age (years)") +
    ylab(ylab) +
    theme_bw(base_size = 14, base_family = "Lato") +
    theme(axis.text        = element_text(face = "bold"),
          strip.background = element_rect(fill = "white"),
          panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
    theme(legend.position = "none") 
    #scale_color_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1") + 
    #scale_fill_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1")
  
  if (percent){
    p <- p + scale_y_continuous(labels = function(x){sprintf("%g%%",x*100)},
                                limits = ylim)
  } else {
    p <- p + scale_y_continuous(limits = ylim)
  }
  
  p
}

# function to identify the age for maximum impact for each vaccine type and country
q <- function(z){
VE_impact_max <- 
  left_join(
    z %>%
      mutate(maxI = `50%`) %>%
      select(serogroup, Waning, country, maxI) %>%
      group_by(serogroup, Waning, country) %>%
      summarise(Impactmax = max(maxI, na.rm=TRUE)),
    
    z %>%
      mutate(maxI = `50%`) %>%
      select(serogroup, Waning, country, age_dep, Vac.age, maxI) %>%
      group_by(serogroup, Waning, country, age_dep, Vac.age) %>%
      summarise(Impactmax = max(maxI, na.rm=TRUE)) 
  )
VE_impact_max
}

# function to identify the age for vaccination efficiency for each vaccine type and country
o <- function(r){
  VE_impact_min <- 
    left_join(
      r %>%
        mutate(minI = `50%`) %>%
        select(serogroup, Waning, country, minI) %>%
        group_by(serogroup, Waning, country) %>%
        summarise(Impactmin = min(minI, na.rm=TRUE)),
      
      r %>%
        mutate(minI = `50%`) %>%
        select(serogroup, Waning, country, age_dep, Vac.age, minI) %>%
        group_by(serogroup, Waning, country, age_dep, Vac.age) %>%
        summarise(Impactmin = min(minI, na.rm=TRUE)) 
    )
  VE_impact_min
}
#===============================================================================================

# define efficacy waning using studies
VE_impact_by_age <- VE_by_Vac.age %>%
    dplyr::group_by(serogroup,
                    Study.waning,
                    age_dep,
                    sim,
                    Vac.age,
                    country) %>%
    dplyr::summarise(Impact = sum(Impact)) %>%
    ungroup %>%
    dplyr::mutate(Waning = sub(pattern     = "None", 
                               replacement = "No waning",
                               x           = Study.waning),
                  Waning = sub(pattern     = "Andrews et al. (2012)",
                               replacement = "Fast waning",
                               x           = Waning),
                  Waning = sub(pattern     = "Djennad et al. (2018)",
                               replacement = "Slow waning", 
                               x           = Waning)) %>%
  mutate(Waning = if_else(Waning == "Andrews et al. (2012)", "Fast waning", "Slow waning"))

# add uncertainty to VE impact (age-independent)
VE_impact_by_ageI_ <- VE_impact_by_age %>%
    filter(!is.na(Impact)) %>%
    nest(data = c(sim, Impact)) %>%
    mutate(Q = map(data, ~quantile(x     = .x$Impact,
                                   probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q) %>%
    select(-data) %>%
  filter(age_dep == FALSE)

# plot age-independent vaccine impact (expected total number of cases averted)
VE_A1 <- make_grid_plot(x = VE_impact_by_ageI_, ylab = "Vaccine impact (expected total cases averted)") +
  geom_point(data = q, aes(x = Vac.age, y = Impactmax), shape = 4, stroke = 1, size = 1)
  
ggsave(filename = "output/Fig2_vaccine_impact_cohort_age_indep.png", 
       plot = VE_A1,
       width = 14, height = 8, units = "in", dpi = 300)

# add uncertainty to VE impact (age-dependent)
VE_impact_by_ageD_ <- VE_impact_by_age %>%
  filter(!is.na(Impact)) %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(x     = .x$Impact,
                                 probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data) %>%
  filter(age_dep == TRUE)

# plot age-dependent vaccine impact (expected total number of cases averted)
VE_A2 <- make_grid_plot(x = VE_impact_by_ageD_, ylab = "Vaccine impact (expected total cases averted)") +
  geom_point(data = q, aes(x = Vac.age, y = Impactmax), shape = 4, stroke = 1, size = 1)

ggsave(filename = "output/S6_Fig_vaccine_impact_cohort_age_dep.png", 
       plot = VE_A2,
       width = 14, height = 8, units = "in", dpi = 300)

#===============================================================================================

# vaccine impact per 100,000 age cohort vaccinees  (age-independent)
VE_impact_validated_ageI_ <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
    dplyr::rename(Vac.age = agey) %>% 
    dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = ntotal/Impact) %>%
    nest(data = c(sim, Impact)) %>%
    mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q) %>%
  filter(age_dep == FALSE)

# impact per 100,000 older adults vaccinated in specific age cohort (ntotal)
VE_B1 <- make_grid_plot(x = VE_impact_validated_ageI_, ylab = "Impact (Number of individuals needed to vaccinate to prevent a case)") +
  geom_point(data = o, aes(x = Vac.age, y = Impactmin), shape = 4, stroke = 1, size = 1)

ggsave(filename = "output/Fig3_vaccine_impact_vaccinee_age_indep.png", 
       plot = VE_B1,
       width = 14, height = 8, units = "in", dpi = 300)

# vaccine impact per 100,000 age cohort vaccinees (age-dependent)
VE_impact_validated_ageD_ <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = ntotal/Impact) %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  filter(age_dep == TRUE)

# impact per 100,000 older adults vaccinated in specific age cohort (ntotal)
VE_B2 <- make_grid_plot(x = VE_impact_validated_ageD_, ylab = "Impact (Number of individuals needed to vaccinate to prevent a case)") +
  geom_point(data = o, aes(x = Vac.age, y = Impactmin), shape = 4, stroke = 1, size = 1)

ggsave(filename = "output/S7_Fig_vaccine_impact_vaccinee_age_dep.png", 
       plot = VE_B2,
       width = 14, height = 8, units = "in", dpi = 300)

#===============================================================================================

#vaccine impact (proportion of preventable IPD cases) e.g., cases that could have been prevented due to vaccination 
impact_per_case <- ipd_mc %>%
    mutate(cases = map(.x = mc, .f = ~group_by(.x, sim) %>%
                           # make it per vaccinee
                           crossing(Vac.age = seq(55, 85, by = 5)) %>%
                           filter(agey >= Vac.age) %>%
                           # end per vaccinee
                           group_by(Vac.age, sim) %>%
                           summarise(cases = sum(fit)))) %>%
    select(-data, -model, -mc) %>%
    unnest(cases) %>%
    inner_join(VE_impact_by_age) %>%
    mutate(rel_impact = Impact/cases) %>%
    group_by_at(.vars = vars(-c(sim, cases, Impact, rel_impact))) %>%
    nest %>%
    mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q)

# plot_impact_per_case <- 
VE_C1 <- make_grid_plot(x = impact_per_case,  percent = TRUE, ylab = "Vaccine impact (proportion of preventable IPD cases)") + 
  geom_point(data = q, aes(x = Vac.age, y = Impactmax), shape = 4, stroke = 1, size = 1) + 
  theme(legend.position = "bottom") + 
  scale_color_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1") + 
  scale_fill_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1")

ggsave(filename = "output/S8_Fig_vaccine_impact_per_vaccinee.png", 
       plot = VE_C1,
       width = 14, height = 8, units = "in", dpi = 300)

#===============================================================================================

# vaccine impact per 100,000 total population (55+y pop)
impact_per_vaccinee <- 
    
    # get the population over 55, 60, etc. as potential vaccinees
    pop_country_df %>%
    crossing(Vac.age = seq(55, 85, by = 5)) %>%
    filter(agey >= Vac.age) %>%
    group_by(country, Vac.age) %>%
    summarise(pop = sum(ntotal)) %>%
    
    # merge with Impact data (averted cases, absolute)
    inner_join(VE_impact_by_age) %>%
    
    # relative impact is per 100,000 total population 55+y
    mutate(rel_impact = Impact/pop*scale) %>%
    group_by_at(.vars = vars(-c(sim, pop, Impact, rel_impact))) %>%
    nest %>%
    mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q)

# plot_impact_per_vaccinee <- 
VE_C2 <- make_grid_plot(x = impact_per_vaccinee, ylab = "Vaccine impact (Cases averted per 100,000 population)") +
  geom_point(data = q, aes(x = Vac.age, y = Impactmax), shape = 4, stroke = 1, size = 1) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1") + 
  scale_fill_brewer(name = "Age dependent vaccine efficacy/effectiveness", palette = "Set1")

ggsave(filename = "output/S9_Fig_vaccine_impact_per_100k_pop.png", 
       plot = VE_C2,
       width = 14, height = 8, units = "in", dpi = 300)
