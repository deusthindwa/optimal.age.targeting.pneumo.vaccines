# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

#loaad datasets into memory
pop_ew <- read_csv(here("data", "total_pop_EW.csv"))
pop_mw <- read_csv(here("data", "total_pop_MW.csv"))
pop_sa <- read_csv(here("data", "total_pop_SA.csv"))
pop_br <- read_csv(here("data", "total_pop_BR.csv"))

#ggplot comparing % populations in England/Wales versus UN SDG regions
pop_country <- c("England/Wales" = "red", 
                 "Malawi" = "Black",
                 "South Africa" = "orange", 
                 "Brazil" = "blue")

pop_totals <- list(`England/Wales` = 56286961 + 3152879, # mid-2019
                   `Malawi`        = 17210000,
                   `South Africa`  = 56200000,
                   `Brazil`        = 206200000) %>%
  map_df(.id = "country", ~data.frame(N = .x))

#set smooth or unsmooth conditions
pop_use_totals <- FALSE
pop_smooth    <- TRUE

#unsmoothed population values
pop_country_df <- list(`England/Wales` = pop_ew,
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

#pop_cases <- dplyr::inner_join(bind_rows(ipd_mc, .id = "serogroup"), 
#                               pop_country_df, by = c("agey")) %>% 
#  dplyr::mutate(cases = fit/1e5*ntotal, Vac.age = agey)

#pop_cases <- rbind(
  
#  filter(pop_cases, serogroup == "All serotypes.Brazil" | serogroup == "PCV13.Brazil" | serogroup == "PPV23.Brazil") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Brazil"),
  
#  filter(pop_cases, serogroup == "All serotypes.England/Wales" | serogroup == "PCV13.England/Wales" | serogroup == "PPV23.England/Wales") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "England/Wales"),
  
# filter(pop_cases, serogroup == "All serotypes.Malawi" | serogroup == "PCV13.Malawi" | serogroup == "PPV23.Malawi") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Malawi"),
  
#  filter(pop_cases, serogroup == "All serotypes.South Africa" | serogroup == "PCV13.South Africa" | serogroup == "PPV23.South Africa") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "South Africa")
#) %>% 
#  dplyr::filter(serogroup != "All serotypes")

#plot smoothed population values and absolute number of cases
pop_country_plot <- ggplot(data = pop_country_df, aes(x = agey, y = p)) +
  geom_col(aes(fill = country), width = 0.8, position = position_dodge(width = 0.8)) +
  scale_x_continuous(breaks  = seq(55, 90, by = 5), labels  = function(x){gsub(pattern = "90", replacement = "90+", x = x)}) + 
  scale_y_continuous(labels  = scales::percent, limits = c(0, NA)) + 
  theme_bw() + 
  labs(title="", x="Age in years (y)", y=paste0("Share of ", ifelse(pop_use_totals, "total national","55y+"), " population"), color = "Countries") +
  scale_fill_manual(values=pop_country) +
  theme(axis.text=element_text(size=10, color="black"), legend.position = "right")

pop_burden_plot <- ggplot(data = pop_cases, aes(x = agey, y = cases, color = serogroup, fill  = serogroup)) +
  geom_line() +
  theme_bw() +
  #geom_ribbon(aes(ymin = lcases, ymax = ucases), alpha = 0.2, color = NA) +
  facet_grid(. ~ country) +
  ylim(c(0, NA)) +
  #coord_cartesian(ylim = c(0, 40)) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  labs(x = "Age (years)", y = "Absolute number of cases") +
  theme(axis.text=element_text(size=10, color="black"), legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(strip.text.x = element_text(size = 14))
  
ggsave(filename = here("output","Fig2_popn_burden.png"), 
       plot = pop_country_plot/pop_burden_plot,
       width = 10, height = 5, units = "in", dpi = 300)