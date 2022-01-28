# written by Deus Thindwa & Samuel Clifford
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# generalised additive models, exponential decay and growth models.
# 1/08/2021-30/12/2021

# load the IPD cases and estimate uncertainty of observed IPD cases
ipd <- readr::read_csv(here("data", "total_incidence.csv")) %>%
  mutate(agey = readr::parse_number(substr(agegroup, 1, 2)),
         obs = (cases/npop)*1e5,
         obs_lci = (exactci(cases, npop, 0.95)$conf.int[1:112])*1e5,
         obs_uci = (exactci(cases, npop, 0.95)$conf.int[113:224])*1e5) 

#---------- FIT USING NLS

# estimate the rest of parameters using a simple linear model
# log(y) = log(alpha0) + beta0*age
fit_model <- function(x){
  
  #set initial parameter values for the model
  model0 <- lm(log(incidence) ~ agey, data = x)
  start = list(alpha = exp(coef(model0)[1]), beta  = coef(model0)[2])
  
  #fit and NLS model
  nls(data = x,
  incidence ~ exp(alpha) * exp(beta*agey),
  nls.control(maxiter = 2000),
  start = start
  )
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

# generate relation table for ggplotting
ipd_curves <- rbind(
  
  filter(ipd_curves, serogroup == "All serotypes.Brazil" | serogroup == "PCV13.Brazil" | serogroup == "PCV20.Brazil" | serogroup == "PPV23.Brazil") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Brazil"),
  
  filter(ipd_curves, serogroup == "All serotypes.England" | serogroup == "PCV13.England" | serogroup == "PCV20.England" | serogroup == "PPV23.England") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "England"),
  
  filter(ipd_curves, serogroup == "All serotypes.Malawi" | serogroup == "PCV13.Malawi" | serogroup == "PCV20.Malawi" | serogroup == "PPV23.Malawi") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "Malawi"),
  
  filter(ipd_curves, serogroup == "All serotypes.South Africa" | serogroup == "PCV13.South Africa" | serogroup == "PCV20.South Africa" | serogroup == "PPV23.South Africa") %>% 
    mutate(serogroup = substr(serogroup,1,5)) %>% mutate(serogroup = if_else(serogroup == "All s", "All serotypes", serogroup), country = "South Africa")
) %>% mutate(`50%` = if_else(`50%` <0, 0, `50%`), `2.5%` = if_else(`2.5%` <0, 0, `2.5%`), `97.5%` = if_else(`97.5%` <0, 0, `97.5%`))

#============================================================================

# plot fitted IPD incidence along with observed IPD cases with uncertainty
  #arrange serogroups so position of line plots correspond to position in the legend
ipd_curves <- ipd_curves %>% mutate(serogroup = factor(serogroup, levels(factor(serogroup))[c(1,4,3,2)]))
ipd <- ipd %>% mutate(serogroup = factor(serogroup, levels(factor(serogroup))[c(1,4,3,2)]))

A <- ggplot() +
  geom_line(data = filter(ipd_curves), aes(x = agey, y = `50%`, color = serogroup), size = 1) +
  geom_ribbon(data = filter(ipd_curves), aes(x = agey, y = `50%`, color = serogroup, fill = serogroup, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, color = NA) +
  geom_point(data = filter(ipd), aes(x = agey, y = obs, color = serogroup, size = cases), shape = 1, stroke = 1.5, position = position_dodge(width = 1)) +
  geom_errorbar(data = filter(ipd), aes(agey, ymin = obs_lci, ymax = obs_uci, color = serogroup), width = 0, size = 0.3, position = position_dodge(width = 1)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  facet_wrap(~country, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 1)) +
  labs(title = "", subtitle = "", x = "Age (years)", y = "IPD incidence per 100,000 population") +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(legend.position = "right") +
  guides(color=guide_legend(title="Serogroup"), fill=guide_legend(title="Serogroup"), size=guide_legend(title="Observed cases")) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

# combined incidence plot
ggsave(here("output", "Fig2_ipd_incidence.png"),
       plot = A,
       width = 10, height = 8, unit="in", dpi = 300)

#============================================================================

# calculate and plot scaled incidence
B <- ipd %>% dplyr::group_by(country, serogroup) %>% dplyr::mutate(p = incidence/sum(incidence)) %>%
  ggplot() + 
  geom_line(aes(x = agey, y = p, color = serogroup), size = 1) +
  theme_bw() +
  labs(title = "", subtitle = "", x = "Age (years)", y = "Scaled IPD incidence") +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(55, 90, 5), limits = c(55, 90)) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(.~country, scales = "free") +
  theme(axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(plot.subtitle = element_text(size = 18, face = "bold", margin = margin(t = 10, b = -25), hjust = 0.02)) +
  theme(plot.title = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(legend.position = "right", legend.text=element_text(size=12), legend.title = element_text(size = 16)) +
  guides(color=guide_legend(title="Serogroup"))

# combined incidence plot
ggsave(here("output", "S1Fig_scaled_incidence.png"),
       plot = (B),
       width = 10, height = 8, unit="in", dpi = 300)
