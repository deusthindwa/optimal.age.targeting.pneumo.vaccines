# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases and estimate uncertainty of observed IPD cases
yrcases <- readr::read_csv(here("data", "yearly_cases.csv")) 

A <- yrcases %>% filter(!is.na(ipd) & serogroup != "All" & agey<90) %>%
  ggplot() +
  geom_line(aes(x = agey, y = ipd, color = survyr), size = 0.8) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age (years)", y = "Reported IPD cases") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Year of IPD surveillance"))
  
B <- yrcases %>% filter(!is.na(ipd) & serogroup == "All" & agey<90) %>%
  mutate(serogroup = if_else(serogroup == "All", "All serotype groups", serogroup)) %>%
  ggplot() +
  geom_line(aes(x = agey, y = ipd, color = survyr), size = 0.8) +
  scale_color_hue(l = 40, c = 35) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age (years)", y = "") +
  theme_bw() + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = ""))
  

  
# combined incidence plot
ggsave(here("output", "S3_yearly_cases.png"),
       plot = (A | B | plot_layout(ncol = 2, width = c(4,1))),
       width = 14, height = 8, unit="in", dpi = 300)

