# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 22/09/2022

# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
}

pacman::p_load(char = c("tidyverse", "here","readstata13", 
                        "scales", "boot", "magrittr",  "mvtnorm", 
                        "zoo", "patchwork", "mgcv", "PropCIs", "showtext", "ggh4x", "splitstackshape"))

options(stringsAsFactors = FALSE)
set.seed = 1988

# model incidences in each age group using (main)
source("script/1_incidence.R")

# plot smoothed populations for all countries (main)
source("script/2_pops.R")

# estimate initial efficacy and waning rates (supplementary)
source("script/3_metacurve.R")

# generate scenarios by serogroup, country, VE, age and waning (supplementary)
source("script/4_vaccination_scenarios.R")

# compute vaccine impact, cases averted (main)
source("script/5_vaccine_impact.R")

# compute isolated impact scenarios (supplementary)
source("script/6_impact_scenarios.R")

# plots and estimates of various stats for the results (main+supplementary)
source("script/7_miscellaneous_plots_est.R")
