# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

#loaad datasets into memory
pop.ew <- read_csv(here("data", "total_pop_EW.csv"))
pop.mw <- read_csv(here("data", "total_pop_MW.csv"))
pop.sa <- read_csv(here("data", "total_pop_SA.csv"))
pop.br <- read_csv(here("data", "total_pop_BR.csv"))

#ggplot comparing % populations in England/Wales versus UN SDG regions
countries <- c("England/Wales" = "red", "Malawi" = "black", "South Africa" = "orange", "Brazil" = "blue")

pop.totals <- list(`England/Wales` = 56286961 + 3152879, # mid-2019
                   `Malawi`        = 18628747,
                   `South Africa`        = 18628747,
                   `Brazil`        = 18628747) %>%
  map_df(.id = "Country", ~data.frame(N = .x))

#set smooth or unsmooth conditions
use.pop.totals <- FALSE
smooth.pops    <- TRUE

#unsmoothed population values
countries_df <- list(`England/Wales` = pop.ew,
                     `Malawi`        = pop.mw,
                     `South Africa`  = pop.mw,
                     `Brazil`        = pop.mw) %>%
  bind_rows(.id = "Country") %>%
  group_by(Country) %>%
  inner_join(pop.totals) %>%
  mutate(p = ntotal/(use.pop.totals*N + (1-use.pop.totals)*sum(ntotal)))

#smoothed population values
if (smooth.pops){
  
  countries_df %<>%
    split(.$Country) %>%
    map(~mutate(.x, N_ = sum(ntotal),
                ntotal = c(head(.x$ntotal, 2),
                           zoo::rollmean(.x$ntotal, 5),
                           tail(.x$ntotal, 2)))) %>%
    map_df(~mutate(.x,
                   p = ntotal/sum(ntotal)))
  
}

#plot smoothed population values
countries_plot <- ggplot(data = countries_df, aes(x = agey, y = p)) +
  geom_col(aes(fill = Country), width = 0.8, position = position_dodge(width = 0.8)) +
  scale_x_continuous(breaks  = seq(55, 90, by = 5), labels  = function(x){gsub(pattern = "90", replacement = "90+", x = x)}) + 
  scale_y_continuous(labels  = scales::percent, limits = c(0, NA)) + 
  theme_bw() + 
  labs(title="", x="Age in years (y)", y=paste0("Age as share of ", ifelse(use.pop.totals, "total national","55y+"), " population"), color = "Countries") +
  scale_fill_manual(values=countries) +
  theme(axis.text=element_text(face="bold", size=10, color="black"), legend.position = "bottom")

ggsave(filename = here("output","Fig1_countries_popn.png"), 
       plot = countries_plot,
       width = 7, height = 3.5, units = "in", dpi = 300)
