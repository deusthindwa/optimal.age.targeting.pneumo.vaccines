
impact_by_age_to_plot_ <-
  impact_by_age_to_plot %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

#vaccine impact plot
impact_by_age_plot <- 
  ggplot(data = impact_by_age_to_plot_,
         aes(x = Vac.age, y= `50%`,
             color  = factor(age_dep),
             group = interaction(Waning, age_dep, serogroup, delay,
                                 Country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = factor(age_dep)),
              color = NA,
              alpha = 0.2) +
  facet_grid(Country ~ serogroup + Waning,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (expected total cases averted)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                     palette = "Set1") + 
  scale_fill_brewer(name = "Age dependent\nvaccine efficacy",
                    palette = "Set1")

ggsave(filename = "output/impact_by_vac_age.png", 
       plot = impact_by_age_plot,
       width = 7, height = 4, units = "in", dpi = 300)


#impact per 10000 older adults vaccinated
impact_validated_plot <- 
  ggplot(data = impact_validated,
         aes(x = Vac.age, y= `50%`,
             color  = factor(age_dep),
             group = interaction(Waning, age_dep, serogroup, delay,
                                 Country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = factor(age_dep)),
              color = NA,
              alpha = 0.2) +
  facet_grid(Country ~ serogroup + Waning,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (cases averted per 10K older adults vaccinated)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                     palette = "Set1") +
  scale_fill_brewer(name = "Age dependent\nvaccine efficacy",
                    palette = "Set1") +
  theme(panel.grid.minor.y = element_blank())



ggsave(filename = "output/impact_validated_by_vac_age.png", 
       plot = impact_validated_plot,
       width = 7, height = 4, units = "in", dpi = 300)
