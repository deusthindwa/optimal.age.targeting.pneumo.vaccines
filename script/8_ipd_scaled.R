# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 31/02/2022

# calculate and plot scaled incidence
ipd_scaled <- ipd %>% 
  dplyr::group_by(country, serogroup) %>% 
  dplyr::mutate(p = incidence/sum(incidence))

B <- 
  ggplot(data = ipd_scaled) + 
  geom_col(aes(x = agey, y = p),
           color = NA) +
  theme_bw(base_size = 14, base_family = 'Lato') +
  labs(x = "Age (years)", y = "Scaled Incidence") +
  scale_y_continuous(limits = c(0, NA), labels = ~sprintf("%g", .)) +
  scale_x_continuous(breaks = seq(60, 90, 10), limits = c(NA, 90)) +
  theme(axis.text        = element_text(face = "bold"),
        legend.position  = "none", 
        legend.text      = element_text(size=12), 
        legend.title     = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  facet_grid(country ~ serogroup) 

# combined incidence plot
ggsave(here("output", "S3_Fig_scaled_incidence.png"),
       plot = B,
       width = 10, height = 8, unit="in", dpi = 300)
