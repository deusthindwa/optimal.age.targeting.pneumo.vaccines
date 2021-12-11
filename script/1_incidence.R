# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# generalised additive models, exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases
ipd <- readr::read_csv(here("data", "ipd_cases_incid.csv"))

ipd <- dplyr::mutate(ipd, 
                     agey = readr::parse_number(substr(agegroup, 1, 2)), 
                     cases = cases/survyr, 
                     incidence = incidence/survyr, 
                     survyr = NULL)

ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup, 1, 2))) 

#------------------------------------------------------- OPTION 1 WITH NLS

# estimate the rest of parameters using a simple linear model
# log(y-theta0) = log(alpha0) + beta0*age
#fit_model <- function(x){
  #set initial parameter values for the model
#  theta0 <- min(x$incidence, na.rm = TRUE)*0.5
#  model0 <- lm(log(incidence-theta0) ~ agey, data = x)
#  alpha0 <- exp(coef(model0)[1])
#  beta0  <- coef(model0)[2]
  
  #fit and NLS model
#  nls(data = x,
#  incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
#  nls.control(maxiter = 2000),
#  start = list(log_alpha = (alpha0), beta  = (beta0), theta = (theta0)))
#}

#ipd_model <- ipd %>% 
#  split(list(.$serogroup, .$country)) %>%
#  purrr::map(~fit_model(.x))

#------------------------------------------------------- OPTION 2 WITH GAM

#fit a GAM to the incidence data for interpolation- and extrapolation to yearly ages
# y = f(x) + e
fit_model <- function(x){ 
  gam(incidence ~ te(agey, bs = "tp"), 
      family = gaussian(link = "identity"),
      data = x)
  }

ipd_model <- ipd %>% 
  split(list(.$serogroup, .$country)) %>%
  purrr::map(~fit_model(.x))

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
ipd_curves <- ipd_mc %>% map_df(summarise_from_model, .id = "serogroup") %>% mutate(serogroup = fct_inorder(factor(serogroup)))

#-------------------------------------------------------

# plot scaled incidence
ipd_curves <- rbind(
  
  filter(ipd_curves, serogroup == "All serotypes.Brazil" | serogroup == "PCV13.Brazil" | serogroup == "PPV23.Brazil") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Brazil"),
  
  filter(ipd_curves, serogroup == "All serotypes.England/Wales" | serogroup == "PCV13.England/Wales" | serogroup == "PPV23.England/Wales") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "England/Wales"),
  
  filter(ipd_curves, serogroup == "All serotypes.Malawi" | serogroup == "PCV13.Malawi" | serogroup == "PPV23.Malawi") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Malawi"),
  
  filter(ipd_curves, serogroup == "All serotypes.South Africa" | serogroup == "PCV13.South Africa" | serogroup == "PPV23.South Africa") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "South Africa")
) %>%
  mutate(`50%` = if_else(`50%` <0, 0, `50%`), `2.5%` = if_else(`2.5%` <0, 0, `2.5%`), `97.5%` = if_else(`97.5%` <0, 0, `97.5%`))

# calculate scaled incidence
ipd_scaled <- ipd %>% dplyr::group_by(country, serogroup) %>% dplyr::mutate(p = incidence/sum(incidence))


#============================================================================

A <- filter(ipd_scaled, country == "England/Wales") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "Scaled Incidence", subtitle = "A, England/Wales", x = "", y = "") +
  ylim(c(0, NA)) +
  xlim(55, 90) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = c(0.2, 0.6))

# plot backward or forward extrapolation incidence
B <- filter(ipd_curves, country == "England/Wales") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "England/Wales"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 60, 10), position = "right") +
  labs(title = "Incidence per 100,000 population per year", x = "Age (years)", subtitle = "B", y = "") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02)) + 
  theme(plot.title = element_text(size = 20)) + 
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

C <- filter(ipd_scaled, country == "Brazil") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "C, Brazil", x = "", y = "") +
  ylim(c(0, NA)) +
  xlim(55, 90) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# plot backward or forward extrapolation incidence
D <- filter(ipd_curves, country == "Brazil") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "Brazil"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 2.5, 0.5), position = "right") +
  labs(title = "D", x = "Age (years)", y = "") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

E <- filter(ipd_scaled, country == "South Africa") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "E, South Africa", x = "", y = "") +
  ylim(c(0, NA)) +
  xlim(55, 90) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# plot backward or forward extrapolation incidence
F <- filter(ipd_curves, country == "South Africa") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "South Africa"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 6.5, 1), position = "right") +
  labs(title = "F", x = "Age (years)", y = "") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

G <- filter(ipd_scaled, country == "Malawi") %>% 
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "G, Malawi", x = "Age (years)", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# plot backward or forward extrapolation incidence
H <- filter(ipd_curves, country == "Malawi") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "Malawi"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 14, 2), position = "right") +
  labs(title = "H", x = "Age (years)", y = "") +
  theme(plot.title = element_text(size = 18, margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


# combined incidence plot
ggsave(here("output", "Fig1_ipd_incidence.png"),
       plot = ((A | B)/(C | D)/(E | F)/(G | H)),
       width = 12, height = 16, unit="in", dpi = 300)

