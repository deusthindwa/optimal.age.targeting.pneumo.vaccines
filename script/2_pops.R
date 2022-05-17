# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 31/02/2022

#loaad datasets into memory
pop_en <- read_csv(here("data", "total_pop_en.csv"))
pop_mw <- read_csv(here("data", "total_pop_mw.csv"))
pop_sa <- read_csv(here("data", "total_pop_sa.csv"))
pop_br <- read_csv(here("data", "total_pop_br.csv"))

#ggplot comparing % populations in England/Wales versus UN SDG regions
pop_country <- c("Brazil" = "blue",
                 "England" = "red", 
                 "Malawi" = "Black",
                 "South Africa" = "orange")

pop_totals <- list(`England` = 16492465, # mid-2017 in England
                   `Malawi`        = 66589, # mid-2018 pop in Blantyre
                   `South Africa`  = 6907974, # projected from 2011 to mid-2016 using Thembisa HIV model 
                   `Brazil`        = 27734412) %>% # projected from 2011 to mid-2016 per world popn growth stats
  map_df(.id = "country", ~data.frame(N = .x))

#set smooth or unsmooth conditions
pop_use_totals <- FALSE
pop_smooth    <- TRUE

#unsmoothed population values
pop_country_df <- list(`England` = pop_en,
                       `Malawi`        = pop_mw,
                       `South Africa`  = pop_sa,
                       `Brazil`        = pop_br) %>%
  bind_rows(.id = "country") %>%
  group_by(country) %>%
  inner_join(pop_totals) %>%
  mutate(p = ntotal/(pop_use_totals*N + (1-pop_use_totals)*sum(ntotal)))

#smoothed population values
if (pop_smooth){
  
  pop_country_df %<>%
    split(.$country) %>%
    map(~mutate(.x, N_ = sum(ntotal),
                ntotal = c(head(.x$ntotal, 2),
                           zoo::rollmean(.x$ntotal, 5),
                           tail(.x$ntotal, 2)))) %>%
    map_df(~mutate(.x,
                   p = ntotal/sum(ntotal)))
  
}

pop_cases <- dplyr::left_join(ipd_curves, pop_country_df, by = c("agey", "country")) %>% 
  dplyr::mutate(cases  = `50%`/scale*ntotal,
                lcases = `2.5%`/scale*ntotal, 
                ucases = `97.5%`/scale*ntotal, Vac.age = agey) %>% 
  filter(serogroup != "All") 

#plot smoothed population values
pop_country_plot <- ggplot(data = pop_country_df, aes(x = agey, y = p)) +
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
  facet_wrap(. ~ country) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) 


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
  #scale_x_continuous(breaks = seq(60, 90, 10), limits = c(55, 85)) + 
  labs(x = "Age (years)", y = "Expected number of cases") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette  = "Dark2") 
  
#ggsave(filename = here("output","Fig1_popn_burden.png"), 
#       plot = pop_country_plot + pop_burden_plot,
#       width = 12, height = 6, units = "in", dpi = 300)

ggsave(filename = here("output","S1_Fig_ipd_burden.png"), 
       plot = pop_burden_plot,
       width = 8, height = 4, units = "in", dpi = 300)
