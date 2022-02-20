# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 31/02/2022

# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
}

pacman::p_load(char = c("tidyverse", "here","readstata13", 
                        "scales", "boot", "magrittr",  "mvtnorm", 
                        "zoo", "patchwork", "mgcv", "PropCIs", "ggh4x"))

options(stringsAsFactors = FALSE)
setwd(here::here())


# model incidences in each age group using 
source(here::here("script", "1_incidence.R"))

# plot smoothed populations for all countries
source(here::here("script", "2_pops.R"))

# estimate vaccine impact against all IPD serotypes
source(here::here("script", "3_metacurve.R"))

# generate scenarios by serogroup, country, VE, age and waning
source(here::here("script", "4_vaccination_scenarios.R"))

# compute vaccine impact
source(here::here("script", "5_vaccine_impact.R"))

# make plots
source(here::here("script", "6_plots.R"))
