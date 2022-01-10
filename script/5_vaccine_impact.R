# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# define efficacy waning using studies
VE_impact_by_age <- VE_by_Vac.age %>%
  dplyr::group_by(serogroup,
                  Study.waning,
                  Study.VE,
                  age_dep,
                  delay,
                  sim,
                  Vac.age,
                  country) %>%
  
  dplyr::summarise(Impact = sum(Impact)) %>%
  
  dplyr::mutate(Waning = dplyr::case_when(
        Study.waning == "None"           ~ "No waning",
        Study.waning == "Andrews et al. (2012)" ~ "Fast waning",
        Study.waning == "Djennad et al. (2018)" ~ "Slow waning"),
        Waning = ifelse(delay > 0, paste(Waning, sprintf("\n(%i years' delay)", delay)), Waning)) 

# add uncertainty to VE impact
VE_impact_by_age_ <- VE_impact_by_age %>% 
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

# plot vaccine impact as in expected number of cases averted
VE_A <- ggplot(VE_impact_by_age_, aes(x = Vac.age, y = `50%`, color = factor(age_dep), group = interaction(Waning, age_dep, serogroup, delay, country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(age_dep)), color = NA, alpha = 0.2) +
  facet_grid(country ~ serogroup + Waning, scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  theme(axis.text=element_text(size=10, color="black")) +
  xlab("Vaccination Age (years)") +
  ylab("Impact (expected total cases averted)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent vaccine efficacy", palette = "Set1") + 
  scale_fill_brewer(name = "Age dependent vaccine efficacy", palette = "Set1")

ggsave(filename = "output/Fig3_vaccine_impact.png", 
       plot = VE_A,
       width = 10, height = 7, units = "in", dpi = 300)

#===============================================================================================

# vaccine impact per 10000 older adults
VE_impact_validated <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
  dplyr::rename(Vac.age = agey) %>% 
  dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
  mutate(Impact = Impact*100000/ntotal) %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

#impact per 10000 older adults vaccinated
VE_B <- ggplot(VE_impact_validated, aes(x = Vac.age, y= `50%`, color = factor(age_dep), group = interaction(Waning, age_dep, serogroup, delay, country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(age_dep)), color = NA, alpha = 0.2) +
  facet_grid(country ~ serogroup + Waning, scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  theme(axis.text=element_text(size=10, color="black")) +
  xlab("Vaccination Age") +
  ylab("Impact (cases averted per 100,000 older adults vaccinated)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent vaccine efficacy", palette = "Set1") +
  scale_fill_brewer(name = "Age dependent vaccine efficacy", palette = "Set1") +
  theme(panel.grid.minor.y = element_blank())

ggsave(filename = "output/Fig4_vaccine_impact_herd.png", 
       plot = VE_B,
       width = 10, height = 7, units = "in", dpi = 300)

