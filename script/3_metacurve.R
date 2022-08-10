# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# piecewise constant models of VE.
# 1/08/2021-30/09/2021

# vaccine efficacy from published papers

dat <- list(
  `Andrews et al. (2012)` = 
    list(
      PPV23 = list(`0-2` = c(48, 32, 60),
                   `2-5` = c(21, 03, 36),
                   `5-Inf` = c(15, -3, 30))
    ),
  
  # `Wright et al. (2013)` = 
  #   list(
  #     PPV23 = list(`0-5` = c(-9, -119, 43),
  #                  `5-10` = c(38.3, -6, 64),
  #                  `10-Inf` = c(-21, -137, 35))
  #   ),
  
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
      All   = list(`0-Inf` = c(52, 22, 77)),
      PCV13 = list(`0-Inf` = c(75, 41, 91))
    )
)

# function to create dataset for published VE excluding Wright's paper
dat_ <- lapply(X = dat, 
               FUN = function(x){
                 map(x,~map_df(.x, ~set_names(.x, c("Mean", "Min", "Max")), 
                               .id = "Ages"))
               }) %>%
  map_df(~bind_rows(.x, .id = "serogroup"), .id = "Study") %>%
  separate(Ages, into = c("xmin", "xmax")) %>%
  mutate_at(.vars = vars(xmin, xmax), .funs = parse_integer) %>%
  mutate(xmax = ifelse(is.na(xmax), 60L, xmax))

source('script/3_metacurve_synthesis.R')

# df <- dat_ %>% 
#   # filter(Study != "Wright et al. (2013)") %>%
#   rename(y = "Mean")


df_from_study <- 
  VE_meta %>%
  #filter(!grepl('Wright', Study)) %>%
  mutate(sd = (Max - Min)/2/1.96,
         xmax = xmax - 0.5)  %>%
  nest(data = -c(serogroup_VE, xmin, xmax)) %>%
  mutate(newdata = map2(.x = xmin, .y = xmax, ~data.frame(t = seq(.x, min(.y, 50))))) %>%
  mutate(newdata = map(.x = newdata, ~crossing(sim = 1:nsims, .x))) %>%
  mutate(newdata = map2(.x = newdata, .y = data,
                        .f = ~mutate(.x, fit = rnorm(n    = nrow(.x),
                                                     mean = .y$Mean,
                                                     sd   = .y$sd)))) %>%
  select(-data) %>%
  unnest(newdata) %>%
  arrange(sim, t, fit)



df_from_study_meta <-
  df_from_study %>%
  select(serogroup_VE, sim, t, fit) %>%
  mutate(t = cut(t, c(0, 2, 5, Inf),
                 include.lowest = T)) %>%
  group_by(serogroup_VE, t, sim) %>%
  summarise(fit = mean(fit)) %>%
  group_by(serogroup_VE, t) %>%
  summarise(y = mean(fit),
            Min = quantile(fit, 0.025),
            Max = quantile(fit, 0.975)) %>%
  ungroup %>%
  mutate(t = gsub(x= t, pattern = "(\\[|\\]|\\)|\\()", 
                  replacement = "")) %>%
  separate(t, into = c('xmin', 'xmax'), sep = ",") %>%
  mutate_at(.vars = vars(xmin, xmax), .funs = parse_number) %>%
  mutate(xmax = ifelse(is.na(xmax), 60, xmax)) 

ggplot(df_from_study_meta) +
  geom_segment(aes(x=xmin, xend = xmax, y = y, yend = y)) +
  geom_rect(color = NA, alpha = 0.2, aes(xmin = xmin, xmax = xmax,
                                         ymin = Min,  ymax = Max)) +
  labs(x = "Years since vaccination", y = "Vaccine efficacy (VE, %)") +
  facet_wrap( ~ serogroup_VE) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.minor.x = element_blank()) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 100))


# 
# # plot of VE and waning rate
# VE_plot <- ggplot(data=df) +
#   geom_segment(aes(x=xmin, xend = xmax, y = y, yend = y)) +
#   geom_segment(data = df_from_study_meta,
#                lty = 2,
#                aes(x=xmin, xend = xmax, y = y, yend = y),
#                color = 'forestgreen') +
#   geom_rect(color = NA, alpha = 0.2, aes(xmin = xmin, xmax = xmax,
#                                          ymin = Min,  ymax = Max)) +
#   geom_rect(data = df_from_study_meta,
#             color = NA, alpha = 0.1, aes(xmin = xmin, xmax = xmax,
#                                          ymin = Min,  ymax = Max),
#             fill = "forestgreen") +
#   labs(x = "Years since vaccination", y = "Vaccine efficacy (VE, %)") +
#   facet_wrap( ~ serogroup_VE) +
#   theme_bw(base_size = 14, base_family = "Lato") +
#   theme(axis.text        = element_text(face = "bold"),
#         strip.background = element_rect(fill = "white"),
#         panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
#   theme(panel.grid.minor.x = element_blank()) +
#   ylim(c(NA,100))
# 
# ggsave("output/S4_Fig_vaccine_efficacy.png",
#        plot = VE_plot,
#        width = 9, height = 6, unit="in", dpi = 300)

