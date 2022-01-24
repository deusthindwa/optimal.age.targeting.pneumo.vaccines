# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# vaccine efficacy from published papers
dat <- list(`Andrews et al. (2012)` = list(`0-2` = c(48, 32, 60),
                                    `2-5` = c(21, 03, 06),
                                    `5-Inf` = c(15, -3, 30)),
            
            `Wright et al. (2013)` = list(`0-5` = c(-9, -119, 43),
                                   `5-10` = c(38.3, -6, 64),
                                   `10-Inf` = c(-21, -137, 35)),
            
            `Rudnick et al. (2013)` = list(`0-5` = c(41.3, 20, 57),
                                    `5-Inf` = c(34.1, 6, 54)),
            
            `Guttierez et al. (2014)` = list(`0-5` = c(44.5, 19, 62),
                                      `5-Inf` = c(32.5, -6, 57)),
            
            `Djennad et al. (2018)` = list(`0-2` = c(41, 23, 54),
                                    `2-5` = c(34, 16, 48),
                                    `5-Inf` = c(23, 12, 32)))

# function to create dataset for published VE excluding Wright's paper
dat_ <- lapply(X = dat, FUN = function(x){
  map_df(x, ~set_names(.x, c("Mean", "Min", "Max")), .id = "Ages")
  }) %>%
  bind_rows(.id = "Study") %>%
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

# simulated VE dataset
df_from_study <- ans_by_study %>%
  map(~list(M = .x$par, V = solve(.x$hessian))) %>%
  map_df(~simulate_from_ans(.x, nsim = nsims, newdata = data.frame(t = seq(0,20))), .id = "Study") %>%
  rename(VE = X1, rate = X2) %>%
  mutate(VE = VE/100)

# simulated VE dataset with mean VE and 95%CI
df_by_study_q <- df_from_study %>%
  nest(data = -c(Study, t)) %>%
  mutate(Q = map(data, ~quantile(.x$fit, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

# summary of initial VE and waning rate and 95%CI from model
ans_by_study_parms_from_model <-
  df_from_study %>%
  distinct(Study, VE, sim, rate) %>%
  mutate(rate = -rate) %>%
  gather(key, value, VE, rate) %>%
  nest(data = c(sim, value)) %>%
  mutate(Q = map(data, ~quantile(.x$value, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data) %>%
  group_by(Study, key) %>%
  transmute(value = sprintf("%0.2f (%0.2f, %0.2f)", `50%`, `2.5%`, `97.5%`)) %>%
  spread(key, value) %>%
  select(Study, `Initial efficacy` = VE, `Waning rate` = rate) %>%
  write_csv(here("output", "initial_VE_and_waning_rates.csv"))

# summary of initial VE and waning rate from optim
ans_by_study_parms <-  ans_by_study %>%
  map("par") %>%
  map_df(.id = "Study", .f = ~data.frame(A = .x[[1]], B = .x[[2]]))

# save to CSV initial VE and waning rate from optim
VE_table <- ans_by_study_parms %>%
  mutate(`Half-life` = -log(2)/B) %>%
  rename(VE = A, rate = B) %>%
  mutate(VE = VE/100)

# plot of VE and waning rate
VE_plot <- ggplot(data=df) +
  geom_segment(aes(x=xmin, xend = xmax, y = y, yend = y)) +
  ylim(c(0, NA)) +
  geom_line(data= df_by_study_q, aes(x=t, y = `50%`)) +
  geom_ribbon(data = df_by_study_q, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.25) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = "Years since vaccination time", y = "Vaccine efficacy (VE, %)") +
  scale_color_brewer(palette = "Set1") +
  guides(col = guide_legend(ncol = 2)) +
  facet_grid(.~Study) +
  theme(strip.text.x = element_text(size = 16)) +
  geom_text(data = ans_by_study_parms, x = 10, y = 50, parse = T, aes(label = paste("VE == ", round(A,1), "*e^{", round(B,3), "*t}"))) 

ggsave("output/Fig S2_vaccine_efficacy.png",
       plot = VE_plot,
       width = 10, height = 4, unit="in", dpi = 300)

