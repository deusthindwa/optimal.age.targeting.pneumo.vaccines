# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# load the IPD cases
ipd <- readr::read_csv(here("data", "ipd_incid_EW.csv"))
ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup, 1, 2)))
Nsims <- 1e3 # number of simulations to use for all uncertainty analysis

# estimate the rest of parameters using a simple linear model
theta0 <- min(ipd$incidence, na.rm = TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data = ipd)  
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

# fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
# we parameterise in terms of log-rates to ensure the estimates are positive
ipd_modelEW <- filter(ipd, country == "Englad/Wales") %>% 
    split(.$serogroup) %>%
    purrr::map(~nls(data = .x, 
                    incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                    nls.control(maxiter = 200),
                    start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0))))

ipd_modelMW <- filter(ipd, country == "Malawi") %>% 
  split(.$serogroup) %>%
  purrr::map(~nls(data = .x, 
                  incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                  nls.control(maxiter = 200),
                  start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))

# function to simulate model to generate uncertainty 
simulate_from_model <- function(x, nsim = 1e4, newdata){
    V <- vcov(x)
    M <- coef(x)
    
    dat <- data.frame(rmvnorm(n = nsim, mean = M, sigma = V)) %>%
        mutate(sim = 1:n()) %>%
        crossing(newdata) %>%
        mutate(fit = predict(x, .)) %>%
        dplyr::select(fit, sim, one_of(names(newdata)))
    
    return(dat)
}

# function to wrap up 95% uncertainty levels
summarise_from_model <- function(x, probs = c(0.025, 0.5, 0.975)){
  nest(x, data = -agey) %>%
    mutate(Q = map(data, ~quantile(.x$fit, probs = probs))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    return
}

# run the simulations
ipd_mcEW <- ipd_modelEW %>% map(~simulate_from_model(.x, newdata = ipd_x, nsim = Nsims))
ipd_mcMW <- ipd_modelMW %>% map(~simulate_from_model(.x, newdata = ipd_x, nsim = Nsims))

# summarise uncertainty
ipd_curvesEW <- ipd_mcEW %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "Englad/Wales")
ipd_curvesMW <- ipd_mcMW %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "Malawi") 

# calculate scaled incidence
ipd_scaled <- ipd %>% dplyr::group_by(serogroup) %>% dplyr::mutate(p = incidence/sum(incidence))

# generate IPD cases from total pop and IPD incidence annually
CasesEW <- dplyr::inner_join(bind_rows(ipd_mcEW, .id="serogroup"), countries_df, by = "agey") %>% dplyr::filter(serogroup != "All serotypes") %>% dplyr::mutate(cases = fit/1e5*ntotal, Vac.age = agey)
CasesMW <- dplyr::inner_join(bind_rows(ipd_mcMW, .id="serogroup"), countries_df, by = "agey") %>% dplyr::filter(serogroup != "All serotypes") %>% dplyr::mutate(cases = fit/1e5*ntotal, Vac.age = agey)

# plot backward or forward extrapolation incidence
incidence_plot <-  
  ggplot(data = rbind(ipd_curvesEW, ipd_curvesMW), aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  facet_grid(serogroup ~ country) +
  ylim(c(0, NA)) + 
  theme_bw() +
  xlab("Age (years)") +
  ylab("Incidence (cases per 100,000)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_point(data = ipd, aes(y = incidence))

ggsave(here("output","incidence_plot.png"),
       plot = incidence_plot,
       width = 7, height = 5, unit="in", dpi = 300)


#scaled incidence plot
scaled_incidence_plot <-
  ggplot(data = ipd_scaled, aes(x = agey, y = p, color = serogroup)) + 
  geom_line() +
  theme_bw() +
  xlab("Age (years)") +
  ylab("Observed incidence") +
  ylim(c(0, NA)) +
  scale_color_brewer(palette = "Dark2", guide = FALSE)


ggsave(here("output","scaled_plot.png"),
       plot = scaled_incidence_plot,
       width = 7, height = 5, unit="in", dpi = 300)


incidence_plots <- 
  scaled_incidence_plot +
  incidence_plot +
  patchwork::plot_layout(nrow = 1,
                         guides = "collect") &
  theme(legend.position='bottom')

ggsave(here("output", "incidence_plots.png"),
       plot = incidence_plots,
       width = 7, height = 4, unit="in", dpi = 300)

