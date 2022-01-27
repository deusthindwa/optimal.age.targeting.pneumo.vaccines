
M <- c(-0.09431694, 0.05897346, 0.37410736, 0.68986941 )
V <- rbind(c(0.006599189, 0.005693222, 0.001009099, -0.004859266),
           c(0.005693222, 0.012744563, 0.010991425, 0.001501227),
           c(0.001009099, 0.010991425, 0.020098142, 0.016832205),
           c(-0.004859266, 0.001501227, 0.016832205, 0.036053482)
                  )
rownames(V) <- c("a","b", "c", "d")
colnames(V) <- c("a","b", "c", "d")


dat <- data.frame(rmvnorm(n = 1000, mean = M, sigma = V)) %>%
  mutate(sim = 1:n()) %>%
  crossing(ipd_x) %>%
  mutate(fit = predict(ipd_model$`All serotypes.Brazil`,.)) %>%
  dplyr::select(fit, sim, one_of(names(ipd_x)))

X <- ipd_mc$`All serotypes.Brazil`
X %>% group_by(agey) %>% summarise(quantile(fit, probs = 0.025), quantile(fit, probs = 0.5), quantile(fit, probs = 0.975))




slope = function(data, indices){
  data = data[indices,]
  incidence = data[,1]
  agey = data[,2]
  model_crude = gam(incidence ~ te(agey, bs = "tp"), family = gaussian(link = "identity"))
  predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
}

data = select(filter(ipd, country == "England/Wales" & serogroup == "PCV13"), incidence, agey)
set.seed(1988)
slopeb = as_tibble(t(boot(data, slope, R = 1000, parallel = "multicore", ncpus = 3)$t))

slopec <- slopeb %>% 
  rowwise() %>% 
  mutate(fit = mean(c_across(c(V1:V1000)), na.rm = TRUE), 
         se = (sd(c_across(c(V1:V1000)), na.rm = TRUE))/length(c_across(c(V1:V1000))),
         fit_lci = fit - qt(1-(0.05/2), length(c_across(c(V1:V1000)))-1)*se,
         fit_uci = fit + qt(1-(0.05/2), length(c_across(c(V1:V1000)))-1)*se) %>%
  select(fit, fit_lci, fit_uci)

crude <- cbind(crude, slopec)




# Fit models to 100 bootstrap replicates of the data

predictions = replicate(
  1000,{
    boot = data[sample.int(nrow(data), replace = TRUE), ]
    model = gam(incidence ~ te(agey, bs = "tp"), family = nb(link = "log"), data = data)
    
    # Output predictions at each point that we'll want to plot later
    predict(model, data.frame(x = ipd_x))
  }
)
plot(predictions)

