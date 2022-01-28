# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

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

pop_totals <- list(`England` = 16492465, # mid-2017 in 
                   `Malawi`        = 66589, # mid-2018 pop in Blantyre
                   `South Africa`  = 6907974, # mid-2016
                   `Brazil`        = 27734412) %>% # projected from 2011 to 2016
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
  dplyr::mutate(cases = `50%`/1e5*ntotal, lcases = `2.5%`/1e5*ntotal, ucases = `97.5%`/1e5*ntotal, Vac.age = agey) %>% 
  filter(serogroup != "All serotypes") 

#plot smoothed population values
pop_country_plot <- ggplot(data = pop_country_df, aes(x = agey, y = p)) +
  geom_col(aes(fill = country), width = 0.8, position = position_dodge(width = 0.8)) +
  scale_x_continuous(breaks  = seq(55, 90, by = 5), labels  = function(x){gsub(pattern = "90", replacement = "90+", x = x)}) + 
  scale_y_continuous(labels  = scales::percent, limits = c(0, NA)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.01)) +
  labs(title="A", x="Age in years (y)", y=paste0("Share of ", ifelse(pop_use_totals, "total national","55y+"), " population"), color = "Countries") +
  scale_fill_manual(values = pop_country) +
  theme(axis.text=element_text(size=10, color="black"), legend.position = "right") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(fill=guide_legend(title="Serotype group"))


#plot absolute number of cases observed from each country
pop_burden_plot <- pop_cases %>% 
  mutate(countryx = if_else(country == "Brazil", "B, Brazil",
                            if_else(country == "England", "C, England",
                                    if_else(country == "Malawi", "D, Malawi", "E, South Africa")))) %>%
  ggplot() +
  geom_line(aes(x = agey, y = cases, color = serogroup, fill  = serogroup)) +
  geom_ribbon(aes(x = agey, y = cases, color = serogroup, fill = serogroup, ymin = lcases, ymax = ucases), alpha = 0.2, color = NA) +
  theme_bw() +
  facet_wrap(countryx~., scales = "free_y", nrow = 1) + 
  ylim(c(0, NA)) + 
  scale_x_continuous(breaks = seq(55, 85, 5), limits = c(55, 85)) + 
  labs(x = "Age (years)", y = "Absolute number of cases\n(adjusted for smoothed population)") +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "right") +
  theme(strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(strip.text.x = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(fill = guide_legend(title="Serotype group"), color = guide_legend(title="Serotype group"))

ggsave(filename = here("output","Fig1_popn_burden.png"), 
       plot = pop_country_plot/pop_burden_plot,
       width = 12, height = 6, units = "in", dpi = 300)
