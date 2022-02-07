# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases and estimate uncertainty of observed IPD cases
yrcases <- readr::read_csv(here("data", "yearly_cases.csv")) 

A <- yrcases %>% filter(!is.na(ipd) & serogroup != "All" & agey<90) %>%
  ggplot() +
  geom_line(aes(x = agey, y = ipd, color = survyr), size = 1.3) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age (years)", y = "Reported IPD cases") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background=element_rect(fill="white")) +
  theme(axis.text.x = element_text(face = "bold", size = 18), axis.text.y = element_text(face = "bold", size = 18)) +
  theme(legend.position = "bottom", legend.text=element_text(size = 18), legend.title = element_text(size = 18)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = ""))
  
B <- yrcases %>% filter(!is.na(ipd) & serogroup == "All" & agey<90) %>%
  mutate(serogroup = if_else(serogroup == "All", "All serotype groups", serogroup)) %>%
  ggplot() +
  geom_line(aes(x = agey, y = ipd, color = survyr), size = 1.3) +
  scale_color_hue(l = 40, c = 35) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age (years)", y = "") +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background=element_rect(fill="white")) +
  theme(axis.text.x = element_text(face = "bold", size = 18), axis.text.y = element_text(face = "bold", size = 18)) +
  theme(legend.position = "bottom", legend.text=element_text(size = 18), legend.title = element_text(size = 18)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = ""))
  

  
# combined incidence plot
ggsave(here("output", "S3_yearly_cases.png"),
       plot = (A | B | plot_layout(ncol = 2, width = c(4,1))),
       width = 20, height = 9, unit="in", dpi = 300)

