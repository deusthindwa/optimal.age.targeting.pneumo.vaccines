
# proportion of people <65y by country
pop_country_df %>% 
  filter(agey<65) %>%
  group_by(country) %>%
  tally(p)


# total number IPD cases and proportion of cases <65y by country
cbind (
  ipd %>% 
    group_by(country) %>%
    tally(cases) %>%
    rename("total" = n),
  
  ipd %>% 
    filter(agey<65) %>%
    group_by(country) %>%
    tally(cases) %>%
    select(n)
  ) %>%
  mutate(p = n/total)


# proportion of serotypes targeted by vaccines by country
ipd %>% 
  filter(serogroup != "All") %>%
  group_by(country, serogroup) %>%
  tally(cases) %>%
  mutate(p = n/sum(n))


# plot population, raw IPD cases, and incidence 
# requires that incidence.R and pops.R are run
A <- ggplot(data = pop_country_df, aes(x = agey, y = p)) +
  geom_col(width = 0.8) +
  scale_x_continuous(breaks  = seq(60, 90, by = 10),
                     labels  = function(x){gsub(pattern = "90", replacement = "90+", x = x)}) + 
  scale_y_continuous(labels  = ~scales::percent(x = ., accuracy = 1),
                     limits  = c(0, NA)) + 
  labs(x = "Age (years)", 
       y = paste0("Share of ",
                  ifelse(pop_use_totals, "total national", "55y+"),
                  " population")) +
  #coord_flip() +
  facet_grid(country ~ .) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) 


B <-ipd %>% 
  #filter(serogroup != "All") %>%
  mutate(ageg = if_else(agegroup == "55-59", "55-",
                        if_else(agegroup == "60-64", "60-",
                                if_else(agegroup == "65-69", "65-",
                                        if_else(agegroup == "70-74", "70-",
                                                if_else(agegroup == "75-79", "75-",
                                                        if_else(agegroup == "80-84", "80-", "85-"))))))) %>%
  select(country, ageg, serogroup, cases) %>%
  ggplot() +
  geom_bar(aes(x=ageg, y=cases), stat = "identity", width = 0.6, position = position_dodge()) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_x_discrete(breaks = c("60-", "70-", "80-")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age group (years)", y = "Number of IPD cases") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text = element_text(face = "bold"),strip.background = element_rect(fill = "white"), panel.border = element_rect(colour = "black", fill = NA, size = 1)) 


C <- 
  ggplot(data = ipd_curves) +
  geom_line(aes(x = agey, y = `50%`), size = 1) +
  geom_ribbon(aes(x = agey, y = `50%`,
                  ymin = `2.5%`, ymax = `97.5%`),
              alpha = 0.2, color = NA) +
  geom_point(data = ipd, aes(x = agey, y = obs), size = 2) +
  geom_linerange(data = ipd,
                 aes(x = agey, 
                     ymin = obs_lci, 
                     ymax = obs_uci)) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), labels = ~sprintf("%g", .)) +
  labs(x = "Age (years)", 
       y = "IPD incidence per 100,000 population") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) 

ggsave(here("output", "Figure 1.png"),
       plot = (A + B + C + plot_layout(ncol = 3, widths = c(1,3,3))), #G,H,I are plots from the spatial distance - Fig2b_spatially_contacts
       width = 20, height = 10, unit="in", dpi = 300)
