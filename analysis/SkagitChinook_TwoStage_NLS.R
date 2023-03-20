######################################################################
######################################################################
# Non-least squares final models for freshwater and marine life stages
# From: Integrated Population Models with application to Skagit River 
#       Chinook Recovery Evaluation - Near-Term Action (NTA) 2018-0697
# Authors: Sarah L. Schooler, Michael LeMoine
# Collaborators: Casey Ruff, Eric Beamer, Catherine Austin
# Code Prepared By: Sarah L. Schooler
# Contact: sschooler@skagitcoop.org
# 17 March 2023
######################################################################
######################################################################


# load necessary packages
if(!require("here")) {
  install.packages("here")
  library("here")
}
if(!require("nlstools")) {
  install.packages("nlstools")
  library("nlstools")
}
if(!require("nlraa")) {
  install.packages("nlraa")
  library("nlraa")
}
if(!require("FSA")) {
  install.packages("FSA")
  library("FSA")
}
if(!require("FSAmisc")) {
  install.packages("FSAmisc")
  library("FSAmisc")
}
if(!require("stringr")) {
  install.packages("stringr")
  library("stringr")
}


## set directory locations
datadir <- here("data")
jagsdir <- here("jags")
analdir <- here("analysis")
savedir <- here("analysis/cache")

##########################
#### initialize data #####
##########################

# read in smolts, spawners, recruitment estimates, and covariates
smolts <- read.csv(file.path(datadir, "skagit_chinook_smolt.csv"))
spawners <- read.csv(file.path(datadir, "skagit_chinook_esc.csv"))
recruits <- read.csv(file.path(datadir, "skagit_chinook_rec_est.csv"))
dat_cvrs <- read.csv("./data/skagit_chinook_covars.csv")
dat_cvrs <- dat_cvrs[dat_cvrs$year>1992&dat_cvrs$year<2015, ]
dat_cvrs <- dat_cvrs[,-1] 
scl_cvrs <- as.matrix(scale(dat_cvrs)) 
## seperate adult-smolt productivity and smolt-adult productivity covariates
## I've labeled covariates adult-smolt prod. as "as_" and smolt-adult as "sa_"
scl_cvrs_sa <- scl_cvrs[,grep("sa", colnames(scl_cvrs))]
scl_cvrs_as <- scl_cvrs[,grep("as", colnames(scl_cvrs))]

# initialize data frames
# Note: identification as "stock" or "recruit" varies by life stage 
# (marine [smolt to adult, sa], or freshwater [adult to smolt, as])
psa.df <- data.frame("stock" = smolts[smolts$year<2015,"smolt"], 
                     "recruits" = recruits[recruits$year>1992,"recruitment"],
                     "logR" = log(recruits[recruits$year>1992,"recruitment"]),
                     scl_cvrs[,colnames(scl_cvrs_sa)])
pas.df <- data.frame("recruits" = smolts[smolts$year<2015,"smolt"], 
                     "stock" = spawners[spawners$year>1992&spawners$year<2015,
                                        "escapement"],
                     "logR" = log(smolts[smolts$year<2015,"smolt"]),
                     scl_cvrs[,colnames(scl_cvrs_as)])
colnames(psa.df) <- str_remove(names(psa.df), "sa_")
colnames(pas.df) <- str_remove(names(pas.df), "as_")

##########################
####### fit models #######
##########################

# marine productivity (smolt to adult)
# parameterized as a linearized Ricker model
# note: for correct plotting in Skagit_Chinook_Two_Stage_plots,
# covariates in the initial list and equation must be in the same
# order as in the Bayesian covariates file

psa_nls <- nls(
  logR ~ log(a) + log(stock) + gNPGO*NPGO + gCUTI*CUTI + gMEI*MEI +  
    gSST*SST - b*stock, 
  data = psa.df, 
  start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gCUTI = 0, gMEI = 0, gSST = 0))

# freshwater productivity (adult to smolt)
# parameterized as a non-density dependent model
pas_nls <- nls(logR ~ log(a) + log(stock) + gfloodRI*floodRI + gFlowPRI*propRI1 +
                     gFlowSp*flowFebJun + gAirSp*airTJanApr, 
                   data = pas.df,
                   start = list(a = 150, gfloodRI = 0, gFlowPRI = 0, 
                                gFlowSp= 0, gAirSp = 0))

# check results
summary(pas_nls)
summary(psa_nls)

# save model and data files for plotting
save(psa_nls, pas_nls, psa.df, pas.df, file = paste0(savedir, "/", "nls_fit_", 
                            format(Sys.Date(), "%m%d%y")))

##########################
#### model diagnostics ###
##########################
## Examine model fit metrics
residPlot(psa_nls)
residPlot(pas_nls)
