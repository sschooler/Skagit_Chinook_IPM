######################################################################
######################################################################
# Non-least squares model selection for testing density dependence
# From: Integrated Population Models with application to Skagit River 
#       Chinook Recovery Evaluation - Near-Term Action (NTA) 2018-0697
# Authors: Sarah L. Schooler, Michael LeMoine
# Collaborators: Casey Ruff, Eric Beamer, Catherine Austin
# Code Prepared By: Sarah L. Schooler
# Contact: sschooler@skagitcoop.org
# 17 March 2023
# Code partially from: 
# http://derekogle.com/fishR/examples/oldFishRVignettes/StockRecruit.pdf
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
if(!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}

## set directory locations
datadir <- here("data")
jagsdir <- here("jags")
analdir <- here("analysis")
savedir <- here("analysis/cache")

##########################
#### initialize data #####
##########################

# load smolts, spawners, recruitment estimates
smolts <- read.csv(file.path(datadir, "skagit_chinook_smolt.csv"))
spawners <- read.csv(file.path(datadir, "skagit_chinook_esc.csv"))
recruits <- read.csv(file.path(datadir, "skagit_chinook_rec_est.csv"))

# initialize data frames
# Note: identification as "stock" or "recruit" varies by life stage 
psa.df <- data.frame("stock" = smolts[smolts$year<2015,"smolt"], 
                     "recruits" = recruits[recruits$year>1992,"recruitment"],
                     "logR" = log(recruits[recruits$year>1992,"recruitment"]))
pas.df <- data.frame("recruits" = smolts[smolts$year<2015,"smolt"], 
                     "stock" = spawners[spawners$year>1992&spawners$year<2015,
                                        "escapement"],
                     "logR" = log(smolts[smolts$year<2015,"smolt"]))

# equations
rkr_eq <- logR~log(a*stock*exp(-b*stock)) #
bh_eq <- logR~log((a*stock)/(1+b*stock))
null_eq <- logR~log(a*stock)

##########################
####### freshwater #######
## adults to smolts [as] #
##########################

## calculate initial values ##
pas_rkr_s <- srStarts(recruits~stock, data = pas.df, 
                    type = "Ricker", param = 1)
pas_bh_s <- srStarts(recruits~stock, data = pas.df, 
                    type = "BevertonHolt", param = 1)

## fit models ##
# linear - no density dependence - uses bev-holt initial alpha
pas_l.mod <- nls(null_eq, data = pas.df, start = list(a = pas_bh_s$a))
# Beverton-Holt
pas_bh.mod <- nls(bh_eq, data = pas.df, start = pas_bh_s)
# Ricker
pas_rkr.mod <- nls(rkr_eq, data = pas.df, start = pas_rkr_s)

## compare and check models ##
# check residuals
residPlot(pas_l.mod); residPlot(pas_bh.mod); residPlot(pas_rkr.mod)
# compare models
anova(pas_l.mod, pas_bh.mod, pas_rkr.mod)
AIC(pas_l.mod, pas_bh.mod, pas_rkr.mod)
# bootstrapped estimates
pas_l.boot <- boot_nls(pas_l.mod); plot(pas_l.boot) # need to use function from nlraa
pas_bh.boot <- nlsBoot(pas_bh.mod); confint(pas_bh.boot, plot = T)
pas_rkr.boot <- nlsBoot(pas_rkr.mod); confint(pas_rkr.boot, plot = T)
# examine correlation
plot(pas_bh.boot); plot(pas_rkr.boot)


## examine plot ##
ggplot(data = pas.df, aes(x = stock, y = recruits/1e6)) + geom_point() + 
  theme_classic() +
  labs(x = "Escapement", y = "Smolts (millions)") + 
  geom_line(aes(y = coef(pas_bh.mod)[1]*stock/(1+coef(pas_bh.mod)[2]*stock)/1e6, 
                col = "Beverton Holt"), linewidth = 1) +
  geom_line(aes(y = coef(pas_rkr.mod)[1]*stock*exp(-coef(pas_rkr.mod)[2]*stock)/1e6, 
                col = "Ricker"), linewidth = 1) +
  geom_line(aes(y = (coef(pas_l.mod)[1]*stock)/1e6, col = "Density ind."), linewidth = 1) +
  scale_color_manual(name = "",
                     breaks = c("Density ind.","Beverton Holt", "Ricker"),
                     values=c("Density ind." = "blue","Beverton Holt" = "red", 
                              "Ricker" = "purple3")) +
  theme(legend.position = "top") 


##########################
######### marine #########
## smolts to adults [sa] #
##########################

## calculate initial values ##
psa_rkr_s <- srStarts(recruits~stock, data = psa.df, 
                      type = "Ricker", param = 1)
psa_bh_s <- srStarts(recruits~stock, data = psa.df, 
                     type = "BevertonHolt", param = 1)

## fit models ##
# linear - no density dependence - uses bev-holt initial alpha
psa_l.mod <- nls(null_eq, data = psa.df, start = list(a = psa_bh_s$a))
# Beverton-Holt
psa_bh.mod <- nls(bh_eq, data = psa.df, start = psa_bh_s)
# Ricker
psa_rkr.mod <- nls(rkr_eq, data = psa.df, start = psa_rkr_s)

## compare and check models ##
# check residuals
residPlot(psa_l.mod); residPlot(psa_bh.mod); residPlot(psa_rkr.mod)
# compare models
anova(psa_l.mod, psa_bh.mod, psa_rkr.mod)
AIC(psa_l.mod, psa_bh.mod, psa_rkr.mod)
# bootstrapped estimates
psa_l.boot <- boot_nls(psa_l.mod); plot(psa_l.boot) # need to use function from nlraa
psa_bh.boot <- nlsBoot(psa_bh.mod); confint(psa_bh.boot, plot = T)
psa_rkr.boot <- nlsBoot(psa_rkr.mod); confint(psa_rkr.boot, plot = T)
# examine correlation
plot(psa_bh.boot); plot(psa_rkr.boot)

## examine plot
ggplot(data = psa.df, aes(x = stock, y = recruits/1e6)) + geom_point() + 
  theme_classic() +
  labs(x = "Escapement", y = "Smolts (millions)") + 
  geom_line(aes(y = coef(psa_bh.mod)[1]*stock/(1+coef(psa_bh.mod)[2]*stock)/1e6, 
                col = "Beverton Holt"), linewidth = 1) +
  geom_line(aes(y = coef(psa_rkr.mod)[1]*stock*exp(-coef(psa_rkr.mod)[2]*stock)/1e6, 
                col = "Ricker"), linewidth = 1) +
  geom_line(aes(y = (coef(psa_l.mod)[1]*stock)/1e6, col = "Density ind."), linewidth = 1) +
  scale_color_manual(name = "",
                     breaks = c("Density ind.","Beverton Holt", "Ricker"),
                     values=c("Density ind." = "blue","Beverton Holt" = "red", 
                              "Ricker" = "purple3")) +
  theme(legend.position = "top") 
