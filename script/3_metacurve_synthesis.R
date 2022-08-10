waning_rate <- dat_ %>%
  filter(Study != "Patterson et al. (2016)") %>%
  mutate(time = map2(.x = xmin, .y = xmax, ~data.frame(t = seq(.x, .y-1)))) %>%
  unnest(time) %>%
  select(-xmin, -xmax) %>%
  group_by(serogroup, t) %>%
  summarise_at(.vars = vars(Mean, Min, Max), .funs = mean) %>%
  mutate(t = cut(t, c(0, 2, 5, Inf),
                 include.lowest = T, right = FALSE)) %>%
  distinct %>%
  mutate(t = gsub(x= t, pattern = "(\\[|\\]|\\)|\\()", 
                  replacement = "")) %>%
  separate(t, into = c('xmin', 'xmax'), sep = ",") %>%
  mutate_at(.vars = vars(xmin, xmax), 
            .funs = ~ifelse(. == "Inf", Inf, parse_number(.))) %>%
  distinct 

waning_rate <- waning_rate %>% 
  ungroup %>%
  select(-serogroup) %>%
  mutate_at(
    .vars = vars(Mean, Min, Max),
    .funs = ~multiply_by(e1 = .,
                         e2 = 1/waning_rate$Mean[1])) 

VE_meta <- dat_ %>%
  filter(Study == "Patterson et al. (2016)") %>%
  select(serogroup_VE = serogroup, Mean) %>%
  rename(y = Mean) %>%
  crossing(waning_rate) %>%
  mutate_at(.vars = vars(Mean, Min, Max),
            .funs = ~multiply_by(., y)) %>%
  select(-y) 

VE_plot_meta <- 
  ggplot(data = VE_meta) +
  geom_segment(aes(x=xmin, xend = xmax, y = Mean, yend = Mean)) +
  
  geom_rect(color = NA, alpha = 0.2, aes(xmin = xmin, xmax = xmax,
                                         ymin = Min,  ymax = Max)) +
  
  labs(x = "Years since vaccination", y = "Vaccine efficacy (VE, %)") +
  facet_wrap( ~ serogroup_VE ) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.minor.x = element_blank()) +
  ylim(c(NA,100)) +
  xlim(c(0, 20))


ggsave("output/S4_Fig_vaccine_efficacy_meta.png",
       plot = VE_plot_meta,
       width = 9, height = 4.5, unit="in", dpi = 300)
