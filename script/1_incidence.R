# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# generalised additive models, exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases
scale = 100000
ipd <- readr::read_csv(here("data", "total_incidence.csv")) %>%
  mutate(agey = readr::parse_number(substr(agegroup, 1, 2)),
         obs = (cases/npop)*scale,
         obs_lci = (exactci(cases, npop, 0.95)$conf.int[1:70])*scale,
         obs_uci = (exactci(cases, npop, 0.95)$conf.int[71:140])*scale) %>%
  filter(country == "England")
  #filter(serogroup == "All serotypes") 

#------------------------------------------------------- OPTION 1 WITH NLS (but fits badly for South Africa)

# estimate the rest of parameters using a simple linear model
# log(y-theta0) = log(alpha0) + beta0*age
fit_model <- function(x){
  
  #set initial parameter values for the model
  theta0 <- min(x$incidence, na.rm = TRUE)*0.05
  model0 <- lm(log(incidence-theta0) ~ agey, data = x)
  start = list(alpha = exp(coef(model0)[1]), beta  = coef(model0)[2], theta = theta0)
  
  #fit and NLS model
  nls(data = x,
  incidence ~ exp(alpha) * exp(beta*agey) + theta,
  nls.control(maxiter = 2000),
  start = start
  )
  
}


ipd_model <- ipd %>% 
  split(list(.$serogroup, .$country)) %>%
  purrr::map(~fit_model(.x))

#------------------------------------------------------- OPTION 2 WITH GAM

#fit a GAM to the incidence data for interpolation- and extrapolation to yearly ages
# y = f(x) + e
fit_model <- function(x){ 
  #gam(incidence ~ te(agey, bs = "ps"), 
  #    family = gaussian(link = "identity"),
  #    data = x)
  gam(incidence ~ s(agey, k = 3), data = x, method = "REML")
  }

ipd_model <- ipd %>% 
  split(list(.$country)) %>%
  purrr::map(~fit_model(.x))




#--------------------------------------------------

#fit GAM
m <- gam(incidence ~ s(agey, k = 3), data = ipd, method = "REML")

#create multivariate random normal function
rmvn <- function(n, mu, sig){
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

#extract from fitted GAM
Vb <- vcov(m, unconditional = TRUE)
newd <- data.frame(agey = seq(55, 90, by = 1))
pred <- predict(m, newd, se.fit = TRUE)
se.fit <- pred$se.fit

#generate simulations of maximum absolute standardized deviation of the fitted model
set.seed(42)
N <- 10000
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
Cg <- predict(m, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit <- quantile(masd, prob = 0.95, type = 8)

predDS <- as_tibble(pred) %>% mutate(fitL = fit - (2 * se.fit), 
                                     fitU = fit + (2 * se.fit),
                                     simL = fit - (crit * se.fit),
                                     simU = fit + (crit * se.fit)
                                     ) 
predDS <- cbind(newd, predDS)


# plot backward or forward extrapolation incidence
ggplot(data = predDS) +
  geom_point(data = filter(ipd, country == "Malawi"), aes(x = agey, y = obs, size = cases), color = "black",shape = 1, stroke = 1.5) +
  geom_errorbar(data = filter(ipd, country == "Malawi"), aes(agey, ymin = obs_lci, ymax = obs_uci), color = "black", width = 0, size = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(aes(x = agey, y = fit), color = "red", size = 1) +
  theme_bw() +
  geom_ribbon(aes(x = agey, y = fit, ymin = fitL, ymax = fitU), alpha = 0.2, fill = "green") +
  geom_ribbon(aes(x = agey, y = fit, ymin = simL, ymax = simU), alpha = 0.1, fill = "blue") +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 1)) +
  labs(title = "", subtitle = "", x = "Age (years)", y = "IPD incidence per 100,000 population") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


#-------------------------------------------------------

# function to simulate model to generate uncertainty 
simulate_from_model <- function(x, newdata, nsim){
  V <- vcov(x)
  M <- coef(x)
  
  dat <- data.frame(rmvnorm(n = nsim, mean = M, sigma = V)) %>%
    mutate(sim = 1:n()) %>%
    crossing(newdata) %>%
    mutate(fit = predict(x, .)) %>%
    dplyr::select(fit, sim, one_of(names(newdata)))
  
  return(dat)
  
  
  
}

#-------------------------------------------------------

# run the simulations
nsims <- 1e3 # number of simulations to use for all uncertainty analysis
ipd_x <- data.frame(agey = seq(55, 90, by = 1))
ipd_mc <- ipd_model %>% map(~simulate_from_model(.x, newdata = ipd_x, nsim = nsims))

#-------------------------------------------------------

# function to wrap up 95% uncertainty levels
summarise_from_model <- function(x, probs = c(0.025, 0.5, 0.975)){
  nest(x, data = -agey) %>%
    mutate(Q = map(data, ~quantile(.x$fit, probs = probs))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    return
}

# summarise uncertainty
#ipd_curves <- ipd_mc %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)))
ipd_curves <- ipd_mc %>% map_df(summarise_from_model, .id = "country") %>% mutate(serogroup = fct_inorder(factor(country)))

#-------------------------------------------------------

# plot scaled incidence
#ipd_curves <- rbind(
  
#  filter(ipd_curves, serogroup == "All serotypes.Brazil" | serogroup == "PCV13.Brazil" | serogroup == "PCV20.Brazil" | serogroup == "PPV23.Brazil") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Brazil"),
  
  #filter(ipd_curves, serogroup == "All serotypes.England and Wales" | serogroup == "PCV13.England and Wales" | serogroup == "PCV20.England and Wales" | serogroup == "PPV23.England and Wales") %>% 
   # mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "England and Wales"),
  
  #filter(ipd_curves, serogroup == "All serotypes.Malawi" | serogroup == "PCV13.Malawi" | serogroup == "PCV20.Malawi" | serogroup == "PPV23.Malawi") %>% 
    #mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Malawi"),
  
#  filter(ipd_curves, serogroup == "All serotypes.South Africa" | serogroup == "PCV13.South Africa" | serogroup == "PCV20.South Africa" | serogroup == "PPV23.South Africa") %>% 
#    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "South Africa")
#) %>%
#  mutate(`50%` = if_else(`50%` <0, 0, `50%`), `2.5%` = if_else(`2.5%` <0, 0, `2.5%`), `97.5%` = if_else(`97.5%` <0, 0, `97.5%`))

#============================================================================


# plot backward or forward extrapolation incidence
ggplot() +
  geom_line(data = filter(ipd_curves), aes(x = agey, y = `50%`, color = country, fill  = serogroup), size = 1) +
  geom_ribbon(data = filter(ipd_curves), aes(x = agey, y = `50%`, color = country, fill  = serogroup, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  geom_point(data = filter(ipd), aes(x = agey, y = obs, color = country, size = cases), shape = 1, stroke = 1.5) +
  geom_errorbar(data = filter(ipd), aes(agey, ymin = obs_lci, ymax = obs_uci, color = country), width = 0, size = 0.3, position = position_dodge(width = 0.05)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  facet_wrap(.~country, scales = "free") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 1)) +
  labs(title = "", subtitle = "", x = "Age (years)", y = "IPD incidence per 100,000 population") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# combined incidence plot
ggsave(here("output", "Fig2_ipd_incidence.png"),
       plot = (A),
       width = 12, height = 10, unit="in", dpi = 300)


#============================================================================


# calculate scaled incidence
ipd_scaled <- ipd %>% dplyr::group_by(country, serogroup) %>% dplyr::mutate(p = incidence/sum(incidence))


A <- filter(ipd_scaled, country == "Brazil") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "Brazil", subtitle = "A", x = "", y = "Scaled Incidence") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(55, 90, 5), limits = c(55, 90)) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = c(0.3, 0.8), legend.text=element_text(size=12), legend.title = element_text(size = 16))


B <- filter(ipd_scaled, country == "England and Wales") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "England and Wales", subtitle = "C", x = "", y = "") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(55, 90, 5), limits = c(55, 90)) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 

C <- filter(ipd_scaled, country == "Malawi") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "Malawi", subtitle = "E", x = "", y = "") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(55, 90, 5), limits = c(55, 90)) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = "none")

D <- filter(ipd_scaled, country == "South Africa") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "South Africa", subtitle = "G", x = "", y = "") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(55, 90, 5), limits = c(55, 90)) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = "none")

# combined incidence plot
ggsave(here("output", "Fig S1_ipd_incidence.png"),
       plot = (A | B | C | D),
       width = 17, height = 10, unit="in", dpi = 300)

