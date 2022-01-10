# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# generalised additive models, exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases
ipd <- readr::read_csv(here("data", "ipd_incid_all.csv")) 
#%>% filter(country != "Malawi")

ipd <- dplyr::mutate(ipd, 
                     agey = readr::parse_number(substr(agegroup, 1, 2)), 
                     cases = cases, 
                     incidence = incidence, 
                     survyr = NULL) 

ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup, 1, 2))) 

#------------------------------------------------------- OPTION 1 WITH NLS (but fits badly for South Africa)

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
  
  filter(ipd_curves, serogroup == "All serotypes.England and Wales" | serogroup == "PCV13.England and Wales" | serogroup == "PPV23.England and Wales") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "England and Wales"),
  
  filter(ipd_curves, serogroup == "All serotypes.Malawi" | serogroup == "PCV13.Malawi" | serogroup == "PPV23.Malawi") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Malawi"),
  
  filter(ipd_curves, serogroup == "All serotypes.South Africa" | serogroup == "PCV13.South Africa" | serogroup == "PPV23.South Africa") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "South Africa")
) %>%
  mutate(`50%` = if_else(`50%` <0, 0, `50%`), `2.5%` = if_else(`2.5%` <0, 0, `2.5%`), `97.5%` = if_else(`97.5%` <0, 0, `97.5%`))

# calculate scaled incidence
ipd_scaled <- ipd %>% dplyr::group_by(country, serogroup) %>% dplyr::mutate(p = incidence/sum(incidence))


#============================================================================

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


# plot backward or forward extrapolation incidence
B <- filter(ipd_curves, country == "Brazil") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "Brazil"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  labs(title = "", subtitle = "B", x = "Age (years)", y = "Incidence per 100,000 population") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

C <- filter(ipd_scaled, country == "England and Wales") %>% 
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

# plot backward or forward extrapolation incidence
D <- filter(ipd_curves, country == "England and Wales") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "England and Wales"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.1)) +
  labs(title = "", x = "Age (years)", subtitle = "D", y = "") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) + 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(plot.title = element_text(size = 20)) + 
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

E <- filter(ipd_scaled, country == "Malawi") %>% 
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

# plot backward or forward extrapolation incidence
F <- filter(ipd_curves, country == "Malawi") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "Malawi"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.1)) +
  labs(title = "", subtitle = "F", x = "Age (years)", y = "") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

G <- filter(ipd_scaled, country == "South Africa") %>% 
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

# plot backward or forward extrapolation incidence
H <- filter(ipd_curves, country == "South Africa") %>% 
  ggplot(aes(x = agey, y = `50%`, color = serogroup, fill  = serogroup)) +
  geom_point(data = filter(ipd, country == "South Africa"), aes(y = incidence), size = 4) +
  geom_line(size = 1) +
  theme_bw() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  labs(title = "", subtitle = "H", x = "Age (years)", y = "") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#============================================================================

# combined incidence plot
ggsave(here("output", "Fig2_ipd_incidence.png"),
       plot = ((A | C | E | G)/(B | D | F | H)),
       width = 17, height = 10, unit="in", dpi = 300)

