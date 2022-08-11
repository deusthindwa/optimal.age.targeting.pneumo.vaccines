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
  mutate_at(.vars = vars(xmin, xmax), .funs = parse_number) %>%
  mutate(xmax = ifelse(is.na(xmax), 60, xmax))

source('script/3_metacurve_synthesis.R')

# df <- VE_meta %>% 
#   # filter(Study != "Wright et al. (2013)") %>%
#   rename(y = "Mean")


df_from_study <- 
  VE_meta %>%
  #filter(!grepl('Wright', Study)) %>%
  mutate(sd = (Max - Min)/2/1.96,
         xmax = xmax - 0.5)  %>%
  nest(data = -c(serogroup, xmin, xmax)) %>%
  mutate(newdata = map2(.x = xmin, .y = xmax, ~data.frame(t = seq(.x, min(.y, 60))))) %>%
  mutate(newdata = map(.x = newdata, ~crossing(sim = 1:nsims, .x))) %>%
  mutate(newdata = map2(.x = newdata, .y = data,
                        .f = ~mutate(.x, fit = rnorm(n    = nrow(.x),
                                                     mean = .y$Mean,
                                                     sd   = .y$sd)))) %>%
  select(-data) %>%
  unnest(newdata) %>%
  arrange(sim, t, fit)

