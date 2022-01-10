# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

#==================================================================

confint <- function (x, ci = 0.95)
{
  `%>%` <- magrittr::`%>%`
  standard_deviation <- sd(x)
  sample_size <- length(x)
  Margin_Error <- abs(qnorm((1-ci)/2))* standard_deviation/sqrt(sample_size)
  df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x),
                        Margin_Error=Margin_Error,
                        'CI lower limit'=(mean(x) - Margin_Error),
                        'CI Upper limit'=(mean(x) + Margin_Error)) %>%
    tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 )
  return(df_out)
}

#==================================================================

#table 1 scenario 1
scene55 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 2
scene55 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 3
scene55 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "None", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "None", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == TRUE & Study.waning == "None", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 4
scene55 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "None", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "None", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PCV13" & age_dep == FALSE & Study.waning == "None", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 5
scene55 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Andrews et al. (2012)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 6
scene55 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Andrews et al. (2012)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 7
scene55 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Djennad et al. (2018)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Djennad et al. (2018)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == TRUE & Study.waning == "Djennad et al. (2018)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

#table 1 scenario 8
scene55 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Djennad et al. (2018)", Vac.age==55)
scene65 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Djennad et al. (2018)", Vac.age==65)
scene75 <- scenariosp %>% filter(serogroup == "PPV23" & age_dep == FALSE & Study.waning == "Djennad et al. (2018)", Vac.age==75)
confint(scene55$VE); confint(scene55$rate)
confint(scene65$VE); confint(scene65$rate)
confint(scene75$VE); confint(scene75$rate)

