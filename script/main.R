# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 1/08/2021-30/09/2021

# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", "here", "scales", "magrittr",  "mvtnorm", "zoo", "patchwork"))

options(stringsAsFactors = FALSE)
setwd(here::here())


# plot smoothed populations for all countries
source(here::here("script", "pops.R"))


# model incidences in each age group using 
source(here::here("script", "incidence.R"))


# estimate vaccine impact against all IPD serotypes
source(here::here("script", "metacurve.R"))


# generate scenarios by serogroup, country, VE, age and waning
source(here::here("script", "vaccination_scenarios.R"))


# compute vaccine impact
source(here::here("script", "vaccine_impact.R"))


# make plots
source(here::here("script", "plots.R"))
