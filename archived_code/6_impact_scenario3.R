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