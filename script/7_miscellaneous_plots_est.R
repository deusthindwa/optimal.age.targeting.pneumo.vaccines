# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/12/2021

# load the IPD cases and estimate uncertainty of observed IPD cases
yrcases <- readr::read_csv("data/yearly_cases.csv")

A <- yrcases %>% filter(!is.na(ipd) & serogroup != "All" & agey<90) %>%
  ggplot() +
  geom_line(aes(x = agey, y = ipd, color = survyr), size = 1.3) +
  facet_grid(country  ~ serogroup, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Age (years)", y = "Reported IPD cases") +
  theme_bw() +
  theme(axis.title       = element_text(size = 20),
        strip.text       = element_text(size = 20), 
        strip.background = element_rect(fill="white"),
        axis.text        = element_text(face = "bold", size = 18),
        legend.position  = "bottom",
        legend.text      = element_text(size = 18),
        legend.title     = element_text(size = 18)) +
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
ggsave("output/S1_Fig_yearly_cases.png",
       plot = (A | B | plot_layout(ncol = 2, width = c(4,1))),
       width = 20, height = 9, unit="in", dpi = 300)

#================================================================

#plot absolute number of cases observed from each country
pop_burden_plot <- 
  pop_cases %>% 
  filter(agey <= 85) %>%
  ggplot(aes(x = agey, y = cases, color = serogroup, fill  = serogroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcases, ymax = ucases), alpha = 0.2, color = NA) +
  theme_bw() +
  facet_wrap(~country, scales = "free_y") + 
  ylim(c(0, NA)) + 
  labs(x = "Age (years)", y = "Expected number of cases") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette  = "Dark2") 

ggsave(filename = "output/S2_Fig_ipd_burden.png", 
       plot = pop_burden_plot,
       width = 8, height = 4, units = "in", dpi = 300)

#================================================================

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
ggsave("output/S3_Fig_scaled_incidence.png",
       plot = B,
       width = 10, height = 8, unit="in", dpi = 300)

#================================================================

# proportion of people <65y by country
pop_country_df %>% 
  filter(agey<65) %>%
  group_by(country) %>%
  tally(p)

#================================================================

# total number IPD cases and proportion of cases <65y by country
cbind (
  ipd %>% 
    filter(serogroup == "All") %>%
    group_by(country) %>%
    tally(cases/survyr) %>%
    rename("total" = n),
  
  ipd %>% 
    filter(agey<65 & serogroup == "All") %>%
    group_by(country) %>%
    tally(cases/survyr) %>%
    select(n)
  ) %>%
  mutate(p = n/total)

#================================================================

# proportion of serotypes targeted by vaccines by country
ipd %>% 
  #filter(serogroup != "All") %>%
  group_by(country, serogroup) %>%
  tally(cases) %>%
  mutate(p = n/sum(n))

#================================================================

# plot population, raw IPD cases, and incidence (requires that incidence.R and pops.R are run)
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
  geom_point(data = filter(ipd, !is.na(encases)), aes(x = agey, y = incidencex), size = 1, color = "red") +
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

ggsave("output/Fig1_pop_cases_incid.png",
       plot = (A + B + C + plot_layout(ncol = 3, widths = c(1,3,3))),
       width = 20, height = 10, unit="in", dpi = 300)

#================================================================

# incidence rate ratios between 55y and 85y
IncidR <- 
  
cbind(  
  left_join(
    ipd_curves %>% 
      filter((agey == 55) & serogroup == "All") %>% 
      group_by(country) %>%
      dplyr::select(serogroup, country, agey, `2.5%`) %>%
      rename("incid1" = `2.5%`),
    
    ipd %>% 
      filter((agey == 55) & serogroup == "All") %>% 
      group_by(country) %>%
      dplyr::select(country, agey, cases, npop) %>%
      rename("cases1" = cases, "npop1" = npop)
    ),

  left_join(
    ipd_curves %>% 
      filter((agey == 85) & serogroup == "All") %>% 
      group_by(country) %>%
      dplyr::select(serogroup, country, agey, `2.5%`) %>%
      rename("incid2" = `2.5%`),
    
    ipd %>% 
      filter((agey == 85) & serogroup == "All") %>% 
      group_by(country) %>%
      dplyr::select(country, agey, cases, npop)
    ) %>%
    rename("cases2" = cases, "npop2" = npop) %>%
    ungroup() %>%
    dplyr::select(incid2, cases2, npop2)
)

IncidR %>% 
  mutate(irrsd = (1/cases1 + 1/cases2)^0.5,
         irr = incid2/incid1,
         irrl = irr - 1.96*irrsd,
         irru = irr + 1.96*irrsd)

#================================================================

# HIV propensity to determine immunocompetency
propensity <- readr::read_csv("data/hiv_propensity.csv") %>% filter(hiv != "Unknown")

# perform year stratified random sampling
set.seed(1988) #reproducibility
propensity_samp <- stratified(propensity, c("year"), 15) #sampling
propensity_samp %>% group_by(hiv) %>% tally(totalST)

#================================================================
