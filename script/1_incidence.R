# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

# back-inflation of IPD cases in England due to current PPV23 program based on Djennad et al.
cov = 0.70 # vaccination coverage in England
VEa = 0.41 # VE of 41% in 65-69y (though 41% in 65-66y and 34% in 67-69y)
VEb = 0.23 # VE of 23% in 70+y 

ipd <- readr::read_csv("data/total_incidence.csv") %>%
  rename("casesx" = "cases", "incidencex" = "incidence") %>%
  mutate(cases = if_else(country == "England" & agegroup == "65-69", casesx/(cov*(1-VEa) + (1-cov)),
                          if_else( country == "England" & (agegroup == "70-74" | agegroup == "75-79" | agegroup == "80-84" | agegroup == "85+"), casesx/(cov*(1-VEb) + (1-cov)), casesx)),
         incidence = cases/npop*100000,
         encases = if_else(casesx != cases, casesx, NA_real_))

# load the IPD cases and estimate uncertainty of observed IPD cases
scale = 100000
ipd <- 
  ipd %>%
  mutate(agey = readr::parse_number(substr(agegroup, 1, 2)),
         obs = (cases/npop)*scale) %>%
  mutate(serogroup = ifelse(serogroup == "All serotypes", "All", serogroup)) %>%
  mutate(serogroup = factor(serogroup,
                            levels = c("All",
                                       "PPV23",
                                       "PCV20",
                                       "PCV15",
                                       "PCV13"))) %>% 
  dplyr::filter(!is.na(cases))

ipd %<>% nest(data = c(cases, npop)) %>%
  mutate(CI = map(.x = data, ~exactci(.x$cases, .x$npop, conf.level = 0.95)) %>%
           map('conf.int') %>%
           map(~data.frame(obs_lci = .x[1]*scale,
                           obs_uci = .x[2]*scale))) %>%
  unnest_wider(CI) %>%
  unnest_wider(data)


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

ipd_model <- ipd %>% nest(data = -c(serogroup, country)) %>%
  mutate(model = map(.x = data, ~fit_model(.x)))

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
ipd_mc <- mutate(ipd_model, mc = map(.x = model, ~simulate_from_model(.x, ipd_x, nsims)))

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
ipd_curves <- 
  mutate(ipd_mc, curves = map(.x = mc, summarise_from_model)) %>%
  select(-data, -model, -mc) %>%
  unnest(curves)

#-------------------------------------------------------

# generate relation table for ggplotting
# can't have less than 0
ipd_curves %<>% 
  mutate_at(.vars = vars(`2.5%`, `50%`, `97.5%`),
            .funs = ~pmax(0, .))

