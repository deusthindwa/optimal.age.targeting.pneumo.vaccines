# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# vaccine efficacy from published papers
dat <- list(
  `Andrews et al. (2012)` = 
    list(
      PPV23 = list(`0-2` = c(48, 32, 60),
                   `2-5` = c(21, 03, 36),
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
      All   = list(`0-Inf` = c(52, 22, 77)),
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
  mutate(xmax = ifelse(is.na(xmax), 20, xmax))

df <- dat_ %>% filter(Study != "Wright et al. (2013)") %>% rename(y = "Mean")

# function to compute mean and variance given range of VE data 
f <- function(parms, df){
  g <- function(parms, xmin, xmax){
    mean(parms[1]*exp(parms[2]*seq(xmin, xmax, by = 1)))
  }
  df <- df %>% 
    rowwise %>% 
    mutate(ynew = g(parms, xmin, xmax)) %>% 
    ungroup
  return( sum( (df$y - df$ynew)^2 ))
}

# simulating exponential decay model to compute mean VE and 95%CIs
ans_by_study <- df %>% 
  split(.$Study) %>% 
  map(~optim(par = c(50, -0.5), fn = f, df = .x, hessian = TRUE)) 

simulate_from_ans <- function(x, nsim = 1e4, newdata){
  f <- function(x){
    mutate(x, fit = X1*exp(X2*t))
  }
  dat <- data.frame(rmvnorm(n = nsim, mean = x$M, sigma = x$V)) %>% 
    mutate(sim = 1:n()) %>% 
    crossing(newdata) %>% 
    f 
  return(dat)
}


df_from_study <- 
  dat_ %>%
  filter(!grepl('Wright', Study)) %>%
  mutate(sd = (Max - Min)/2/1.96,
         xmax = xmax - 0.5)  %>%
  nest(data = -c(Study, xmin, xmax)) %>%
  mutate(newdata = map2(.x = xmin, .y = xmax, ~data.frame(t = seq(.x, min(.y, 50))))) %>%
  mutate(newdata = map(.x = newdata, ~crossing(sim = 1:nsims, .x))) %>%
  mutate(newdata = map2(.x = newdata, .y = data,
                        .f = ~mutate(.x, fit = rnorm(n    = nrow(.x),
                                                     mean = .y$Mean,
                                                     sd   = .y$sd)))) %>%
  select(-data) %>%
  unnest(newdata) %>%
  arrange(Study, sim, t, fit)

# simulated VE dataset
# df_from_study <- ans_by_study %>%
#   map(~list(M = .x$par, V = solve(.x$hessian))) %>%
#   map_df(~simulate_from_ans(.x, nsim = nsims, 
#                             newdata = data.frame(t = seq(0,20))),
#          .id = "Study") %>%
#   rename(VE = X1, rate = X2) %>%
#   mutate(VE = VE/100)

# simulated VE dataset with mean VE and 95%CI
df_by_study_q <- df_from_study %>%
  nest(data = -c(Study, t)) %>%
  mutate(Q = map(data, ~quantile(.x$fit, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

# plot of VE and waning rate
VE_plot <- ggplot(data=df) +
  geom_segment(aes(x=xmin, xend = xmax, y = y, yend = y)) +
  geom_rect(color = NA, alpha = 0.2, aes(xmin = xmin, xmax = xmax,
                                         ymin = Min,  ymax = Max)) +
  # ylim(c(0, NA)) +
  # geom_line(data= df_by_study_q, aes(x=t, y = `50%`)) +
  # geom_ribbon(data = df_by_study_q, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  labs(x = "Years since vaccination", y = "Vaccine efficacy (VE, %)") +
  facet_wrap(~ Study + serogroup) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(axis.text        = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.minor.x = element_blank())

ggsave("output/S2Fig_vaccine_efficacy.png",
       plot = VE_plot,
       width = 6, height = 4.5, unit="in", dpi = 300)

