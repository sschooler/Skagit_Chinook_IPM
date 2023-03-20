######################################################################
######################################################################
# Non-least squares model selection for testing correlated covariates
# From: Integrated Population Models with application to Skagit River 
#       Chinook Recovery Evaluation - Near-Term Action (NTA) 2018-0697
# Authors: Sarah L. Schooler, Michael LeMoine
# Collaborators: Casey Ruff, Eric Beamer, Catherine Austin
# Code Prepared By: Sarah L. Schooler
# Contact: sschooler@skagitcoop.org
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
# if(!require("FSA")) {
#   install.packages("FSA")
#   library("FSA")
# }
# if(!require("FSAmisc")) {
#   install.packages("FSAmisc")
#   library("FSAmisc")
# }
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
dat_cvrs <- read.csv("./data/skagit_chinook_covars_all.csv")
dat_cvrs <- dat_cvrs[dat_cvrs$year>1992&dat_cvrs$year<2015, ]
# create marine (smolt-adult [sa] data frame)
psa.df <- data.frame("smolts" = smolts[smolts$year<2015,"smolt"], 
                     "recruits" = recruits[recruits$year>1992,"recruitment"],
                     "logR" = log(recruits[recruits$year>1992,"recruitment"]),
                     scale(dat_cvrs[,grep("sa", names(dat_psa_cvrs))]))
colnames(psa.df) <- str_remove(names(psa.df), "sa_")

# create freshwater (adult-smolt [as] data frame)
pas.df <-  data.frame("smolts" = smolts[smolts$year<2015,"smolt"], 
                      "escapement" = spawners[spawners$year>1992&spawners$year<2015,
                                         "escapement"],
                      "logS" = log(smolts[smolts$year<2015,"smolt"]),
                      scale(dat_cvrs[,grep("as", names(dat_cvrs))]))
colnames(pas.df) <- str_remove(names(pas.df), "as_")

##########################
# freshwater (as) models #
##########################
# iteratively run all combinations of uncorrelated parameters
# keeping covariates that improve the model
# initial values for covariates = 0; initial values for parameters
# are estimates from null model

# check correlations
round(cor(pas.df[,-c(1:3)]), 2) 
# freshwater (as) correlations (> 0.7): flowPeakW-floodRI, 
#                                       airTJanApr-airTann (derivative)
# run models
null_pas <-  nls(logS ~ log(a*escapement), data = pas.df, start = list(a = 150))
mod1_pas <- nls(logS ~ log(a) + log(escapement) + gFlowW*flowPeakW +  
              gFlowPRI*propRI1 + gFlowSp*flowFebJun + gAirSp*airTJanApr, 
              data = pas.df, start = list(a = 275, gFlowW = 0, gFlowSp = 0, 
                                          gAirSp = 0, gFlowPRI= 0 ))
mod2_pas <- nls(logS ~ log(a) + log(escapement) + gFloodRI*floodRI +  
                gFlowPRI*propRI1 + gFlowSp*flowFebJun + gAirSp*airTJanApr, 
              data = pas.df, start = list(a = 275, gFloodRI = 0, gFlowSp = 0, 
                                          gAirSp = 0, gFlowPRI= 0 ))
mod3_pas <- nls(logS ~ log(a) + log(escapement) + gFlowW*flowPeakW +  
                  gFlowPRI*propRI1 + gFlowSp*flowFebJun + gAirAnn*airTann, 
                data = pas.df, start = list(a = 275, gFlowW = 0, gFlowSp = 0, 
                                            gAirAnn = 0, gFlowPRI= 0 ))
mod4_pas <- nls(logS ~ log(a) + log(escapement) + gFloodRI*floodRI +  
                  gFlowPRI*propRI1 + gFlowSp*flowFebJun + gAirAnn*airTann, 
                data = pas.df, start = list(a = 275, gFloodRI = 0, gFlowSp = 0, 
                                            gAirAnn = 0, gFlowPRI= 0 ))
# examine models
AIC(mod1_pas, mod2_pas, mod3_pas, mod4_pas)
View(data.frame(AIC(mod1_pas, mod2_pas, mod3_pas, mod4_pas), 
                formula = as.character(lapply(list(
                  mod1_pas, mod2_pas, mod3_pas, mod4_pas), formula))))

##########################
# marine (sa) models #
##########################
# iteratively run all combinations of uncorrelated parameters
# keeping covariates that improve the model
# initial values for covariates = 0; initial values for parameters
# are estimates from null model

# examine correlation between predictors
round(cor(psa.df[,-c(1:3)]), 2) 
# marine (sa) correlations (> 0.7): CUTI-BEUTI, BEUTI-WASST, MEI-ONI, 
#                                   ONI-PDO, ONI-deltaT, ONI-SST, PDO-SST, 
#                                   PDO-WASST, deltaT-SST, SST-WASST
# run models
null_psa <- nls(logR ~ log(a) + log(smolts) - b*smolts, data = psa.df, 
            start = list(a = .01, b = 1.5e-7))
mod1_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gCUTI*CUTI + 
    gONI*ONI + gWASST*WASST - b*smolts, data = psa.df, 
  start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gCUTI = 0, gONI = 0, 
               gWASST = 0))
mod2_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gBEUTI*BEUTI + 
                  gONI*ONI - b*smolts, data = psa.df, 
                start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gBEUTI = 0, gONI = 0))
mod3_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gCUTI*CUTI + 
                  gMEI*MEI + gPDO*PDO + gdeltaT*deltaT - b*smolts, data = psa.df, 
                start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gCUTI = 0, gMEI = 0, 
                             gPDO = 0, gdeltaT = 0))
mod4_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gBEUTI*BEUTI + 
                  gMEI*MEI + gPDO*PDO + gdeltaT*deltaT - b*smolts, data = psa.df, 
                start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gBEUTI = 0, gMEI = 0, 
                             gPDO = 0, gdeltaT = 0))
mod5_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gCUTI*CUTI + 
                  gMEI*MEI + gWASST*WASST + gdeltaT*deltaT - b*smolts, data = psa.df, 
                start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gCUTI = 0, gMEI = 0, 
                             gWASST = 0, gdeltaT = 0))
mod6_psa <- nls(logR ~ log(a) + log(smolts) + gNPGO*NPGO + gCUTI*CUTI + 
                  gMEI*MEI + gSST*SST - b*smolts, data = psa.df, 
                start = list(a = 0.01, b = 1.5e-7, gNPGO = 0, gCUTI = 0, gMEI = 0, 
                             gSST = 0))

# examine models
AIC(mod1_psa, mod2_psa, mod3_psa, mod4_psa, mod5_psa, mod6_psa)
View(data.frame(AIC(mod1_psa, mod2_psa, mod3_psa, mod4_psa, mod5_psa, mod6_psa), 
                formula = as.character(lapply(list(
                  mod1_psa, mod2_psa, mod3_psa, mod4_psa, mod5_psa, mod6_psa),
                  formula))))

