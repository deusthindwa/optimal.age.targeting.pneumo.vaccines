# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# exponential decay and growth models.
# 31/02/2022

# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
}

pacman::p_load(char = c("tidyverse", "here","readstata13", 
                        "scales", "boot", "magrittr",  "mvtnorm", 
                        "zoo", "patchwork", "mgcv", "PropCIs", "showtext", "ggh4x"))

options(stringsAsFactors = FALSE)
setwd(here::here())


# model incidences in each age group using (main)
source(here::here("script", "1_incidence.R"))

# plot smoothed populations for all countries (main)
source(here::here("script", "2_pops.R"))

# estimate initial efficacy and waning rates (supplementary)
source(here::here("script", "3_metacurve.R"))

# generate scenarios by serogroup, country, VE, age and waning (supplementary)
source(here::here("script", "4_vaccination_scenarios.R"))

# compute vaccine impact, cases averted (main)
source(here::here("script", "5_vaccine_impact.R"))

# compute isolated impact scenarios (supplementary)
source(here::here("script", "6_impact scenario.R"))

# make plots on raw IPD data from various countries (supplementary)
source(here::here("script", "7_yearly_cases.R"))

# estimate initial efficacy and waning rates (supplementary)
source(here::here("script", "8_ipd_scaled.R"))

