# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# load the IPD cases
ipd <- readr::read_csv(here("data", "ipd_cases_incid.csv"))
ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup, 1, 2)), cases = cases/survyr, survyr = NULL)
Nsims <- 1e3 # number of simulations to use for all uncertainty analysis

# estimate the rest of parameters using a simple linear model
#theta0 <- min(ipd$incidence, na.rm = TRUE)*0.5
#model0 <- lm(log(incidence-theta0) ~ agey, data = ipd)
#alpha0 <- exp(coef(model0)[1])
#beta0  <- coef(model0)[2]

# fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
# we parameterise in terms of log-rates to ensure the estimates are positive

EW <- filter(ipd, country == "England/Wales")
theta0 <- min(EW$incidence, na.rm = TRUE)*0.5
model0 <- lm(log(incidence-theta0) ~ agey, data = EW)
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

ipd_modelEW <- filter(ipd, country == "England/Wales") %>% 
    split(.$serogroup) %>%
    purrr::map(~nls(data = .x, 
                    incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                    nls.control(maxiter = 2000),
                    start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0))))

#-------------------------------------------------------

MW <- filter(ipd, country == "Malawi")
theta0 <- min(MW$incidence, na.rm = TRUE)*0.5
model0 <- lm(log(incidence-theta0) ~ agey, data = MW)
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

ipd_modelMW <- filter(ipd, country == "Malawi") %>% 
  split(.$serogroup) %>%
  purrr::map(~nls(data = .x, 
                  incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                  nls.control(maxiter = 2000),
                  start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0))))

#-------------------------------------------------------

SA <- filter(ipd, country == "South Africa")
theta0 <- min(SA$incidence, na.rm = TRUE)*0.5
model0 <- lm(log(incidence-theta0) ~ agey, data = SA)
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

ipd_modelSA <- filter(ipd, country == "South Africa") %>% 
  split(.$serogroup) %>%
  purrr::map(~nls(data = .x, 
                  incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                  nls.control(maxiter = 2000),
                  start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0))))

#-------------------------------------------------------

BR <- filter(ipd, country == "Brazil")
theta0 <- min(BR$incidence, na.rm = TRUE)*0.5
model0 <- lm(log(incidence-theta0) ~ agey, data = BR)
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

ipd_modelBR <- filter(ipd, country == "Brazil") %>% 
  split(.$serogroup) %>%
  purrr::map(~nls(data = .x, 
                  incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                  nls.control(maxiter = 2000),
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
ipd_mcSA <- ipd_modelSA %>% map(~simulate_from_model(.x, newdata = ipd_x, nsim = Nsims))
ipd_mcBR <- ipd_modelBR %>% map(~simulate_from_model(.x, newdata = ipd_x, nsim = Nsims))

# summarise uncertainty
ipd_curvesEW <- ipd_mcEW %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "England/Wales")
ipd_curvesMW <- ipd_mcMW %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "Malawi") 
ipd_curvesSA <- ipd_mcSA %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "South Africa") 
ipd_curvesBR <- ipd_mcBR %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)), country = "Brazil") 

# calculate scaled incidence
ipd_scaled <- ipd %>% dplyr::group_by(country, serogroup) %>% dplyr::mutate(p = incidence/sum(incidence))

# plot scaled incidence
ipd_A <- ggplot(data = ipd_scaled, aes(x = agey, y = p, color = serogroup)) + 
  geom_line() +
  theme_bw() +
  facet_grid(. ~ country) +
  labs(x = "", y = "Scaled incidence") +
  ylim(c(0, NA)) +
  xlim(55, 90) +
  scale_color_brewer(palette = "Dark2", guide = FALSE) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# plot backward or forward extrapolation incidence
ipd_B <- ggplot(data = rbind(ipd_curvesEW, ipd_curvesMW, ipd_curvesSA, ipd_curvesBR), aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = ipd, aes(y = incidence)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  facet_grid(. ~ country) +
  #ylim(c(0, NA)) + 
  coord_cartesian(ylim = c(0, 50)) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  labs(x = "Age (years)", y = "Incident cases per \n100,000 population") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

# combined incidence plot
ggsave(here("output", "Fig1_ipd_incidence.png"),
       plot = ipd_A/ipd_B,
       width = 10, height = 5, unit="in", dpi = 300)
