# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

# vaccine efficacy from published papers
dat <- list(
  `Andrews et al. (2012)` = 
    list(
      PPV23 = list(`0-2` = c(48, 32, 60),
                   `2-5` = c(31, 03, 36),
                   `5-Inf` = c(15, -3, 30))
    ),
  
  `Wright et al. (2013)` = 
    list(
      PPV23 = list(`0-5` = c(-9, -119, 43),
                   `5-10` = c(38.3, -6, 64),
                   `10-Inf` = c(-21, -137, 35))
    ),
  
  `Rudnick et al. (2013)` = 
    list(
      PPV23 = list(`0-5` = c(41.3, 20, 57),
                   `5-Inf` = c(34.1, 6, 54))
    ),
  
  `Guttierez et al. (2014)` = 
    list(
      PPV23 = list(`0-5` = c(44.5, 19, 62),
                   `5-Inf` = c(32.5, -6, 57))
    ),
  
  `Djennad et al. (2018)` = 
    list(
      PPV23 = list(`0-2` = c(41, 23, 54),
                   `2-5` = c(34, 16, 48),
                   `5-Inf` = c(23, 12, 32))
    ),
  
  `Patterson et al. (2016)` = 
    list(
      PCV13 = list(`0-Inf` = c(75, 41, 91))
    )
)

# function to create dataset for published VE excluding Wright's paper
dat_ <- lapply(X = dat, 
               FUN = function(x){
                 map(x,~map_df(.x, ~set_names(.x, c("Mean", "Min", "Max")), .id = "Ages"))
               }) %>%
  map_df(~bind_rows(.x, .id = "serogroup"), .id = "Study") %>%
  separate(Ages, into = c("xmin", "xmax")) %>%
  mutate_at(.vars = vars(xmin, xmax), .funs = parse_number) %>%
  mutate(xmax = ifelse(is.na(xmax), 30, xmax))

# plot of VE and waning rate
VE_plot <- 
  dat_ %>% filter(Study != "Wright et al. (2013)") %>% rename(y = "Mean") %>%
  mutate(xmax = if_else(Study == "Patterson et al. (2016)", 5, xmax)) %>%
  ggplot() +
  geom_segment(aes(x=xmin, xend = xmax, y = y, yend = y)) +
  geom_rect(color = NA, alpha = 0.2, aes(xmin = xmin, xmax = xmax,
                                         ymin = Min,  ymax = Max)) +
  labs(x = "Years since vaccination", y = "Vaccine efficacy/effectiveness (%)") +
  facet_grid( ~ serogroup + Study) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.minor.x = element_blank()) +
  ylim(c(NA,100))

ggsave("output/S4_Fig_vaccine_efficacy.png",
       plot = VE_plot,
       width = 14, height = 4, unit="in", dpi = 300)

# a function to simulate VE and waning VE dataset
df_from_study <- 
  dat_ %>%
  filter(!grepl('Wright', Study), !grepl('Rudnick', Study), !grepl('Guttierez', Study)) %>%
  mutate(sd = (Max - Min)/2/1.96,
         xmax = xmax - 0.5)  %>%
  nest(data = -c(Study, serogroup, xmin, xmax)) %>%
  mutate(newdata = map2(.x = xmin, .y = xmax, ~data.frame(t = seq(.x, min(.y, 50))))) %>%
  mutate(newdata = map(.x = newdata, ~crossing(sim = 1:nsims, .x))) %>%
  mutate(newdata = map2(.x = newdata, .y = data,
                        .f = ~mutate(.x, fit = rnorm(n    = nrow(.x),
                                                     mean = .y$Mean,
                                                     sd   = .y$sd)))) %>%
  select(-data) %>%
  unnest(newdata) %>%
  arrange(Study, sim, t, fit)

# simulated VE dataset with mean VE and 95%CI
# used to check if reflects observed data
df_from_studyCk <- 
  df_from_study %>%
  nest(data = -c(Study, serogroup, t)) %>%
  mutate(Q = map(data, ~quantile(.x$fit, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

# manipulate metacurve so PCV13 is duplicated to make PCV15 and PCV20 for use in vaccination scenarios
df_from_study_mc <-
  rbind(
    df_from_study, 
    
    df_from_study %>% 
      filter(Study == "Patterson et al. (2016)") %>%
      mutate(serogroup = "PCV15"),
    
    df_from_study %>% 
      filter(Study == "Patterson et al. (2016)") %>%
      mutate(serogroup = "PCV20")) %>%
  rename("Study.waning" = "Study")
