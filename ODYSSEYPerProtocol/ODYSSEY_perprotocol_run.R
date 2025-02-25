################################################################################################
# NAME: ODYSSEY_run.r
# AUTHORS: Catherine Li
# DATE STARTED: 05/02/2024
# PURPOSE: Sample code to run per-protocol analyses models using ODYSSEY Per-Protocol Function
# UPDATES:  
################################################################################################

# Clear workspace
rm(list=ls())

# Load Packages
packages <- c("tidyverse", "broom", "data.table", "tidyr", "tidyselect", #For data cleaning
              "splines", #for splines
              "survival", "flexsurv", "ggfortify", #for survival data
              "lmtp") #for ICE g-comp and TMLE
for (package in packages) {
  require(package, character.only=T)
}


# Set Working Directory
setwd("[SET WORKING DIRECTORY HERE]")

# Set time of run
st=format(Sys.time(), "%Y%m%d")

#Load Functions
load("ODYSSEY_ppfxn_20250214_0221.R")
add_spline <- function(data, var, names, df, id = TRIALNO) {
  data <- data %>% ungroup()
  spline <- data %>% select(var) %>% data.frame()
  spline_dat <- data.frame(ns(spline[,1], df = df))
  colnames(spline_dat) <- names
  dat <- cbind(data, spline_dat)
  dat <- dat %>% relocate(names, .after = var) %>% select(-var)
  return(dat)
}

# Import Data
load("Data/dat_96_prot1_v2.2.Rdata") 

#Specification of the IPCW models
#denominator
d_mod <- nltfu_t ~ tot_ab*(as.factor(t_start) + sex + ns(bage_v2, knots = knots_age, Boundary.knots = knots_age_boundary) + mode_3cat + country_5cat + 
  ns(bweight, knots = knots_weight, Boundary.knots = knots_weight_boundary) + ns(bheight, knots = knots_height, Boundary.knots = knots_height_boundary)+
  ns(cd4, knots = knots_cd4, Boundary.knots = knots_cd4_boundary) + ns(lag1_cd4, knots = knots_cd4_1, Boundary.knots = knots_cd4_1_boundary) + ns(lag2_cd4, knots = knots_cd4_2, Boundary.knots = knots_cd4_2_boundary) + 
  ns(lag1_log10vl, knots = knots_vl_1, Boundary.knots = knots_vl_1_boundary) + ns(lag2_log10vl, knots = knots_vl_2, Boundary.knots = knots_vl_2_boundary))

#numerator
n_mod <- nltfu_t ~ as.factor(t_start) + tot_ab

#truncation
tmax <- 15

#Bootstraps
b <- 2000

###############
#Primary analyses (toxicities NOT censored)
###############

#IPCW
#Protocol 1+ missed doses
prot1_ipw <- ODYSSEY_ppfxn(boot=0, #Number of bootstraps, set to 0 if want to use original data 
                       seed=252, #To set seed for the bootstraps replication
                       st, #time of run for saving data
                       data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                       g_mod = d_mod, #IPCW denominator mod
                       n_mod = n_mod, #IPCW numerator mod
                       truncate = F, #if TRUE, truncate IPCW to tmax
                       truncate_max = tmax,
                       estimator="IPW")
prot1_ipw_2000bs <- ODYSSEY_ppfxn(boot=b, #Number of bootstraps
                   seed=252, #To set seed for the bootstraps replication
                   st, #time of run for saving data
                   data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                   g_mod = d_mod, #IPCW denominator mod
                   n_mod = n_mod, #IPCW numerator mod
                   truncate = T, #if TRUE, truncate IPCW to tma
                   truncate_max = tmax,
                   estimator="IPW")

#ICE G-computation
#Protocol 1+ missed doses
prot1_gcomp <- ODYSSEY_ppfxn(boot=0, #Number of bootstraps, set to 0 if want to use original data 
                           seed=252, #To set seed for the bootstraps replication
                           st, #time of run for saving data
                           data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                           estimator="gcomp")
prot1_gcomp_2000bs <- ODYSSEY_ppfxn(boot=b, #Number of bootstraps
                                 seed=252, #To set seed for the bootstraps replication
                                 st, #time of run for saving data
                                 data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                                 estimator="gcomp")

#TMLE
#Protocol 1+ missed doses
prot1_tmle <- ODYSSEY_ppfxn(boot=0, #Number of bootstraps, set to 0 if want to use original data 
                           seed=252, #To set seed for the bootstraps replication
                           st, #time of run for saving data
                           data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                           estimator="TMLE")
prot1_tmle_2000bs <- ODYSSEY_ppfxn(boot=b, #Number of bootstraps
                                 seed=252, #To set seed for the bootstraps replication
                                 st, #time of run for saving data
                                 data = dat_96_prot1_v2, #Data set for Primary Analysis, Protocol 1+
                                 estimator="TMLE")


#END
