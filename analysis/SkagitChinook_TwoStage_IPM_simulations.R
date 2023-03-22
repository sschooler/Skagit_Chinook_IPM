######################################################################
######################################################################
# State-space life-cycle Bayesian integrated population model (IPM):
# application to simulated data sets
# From: Integrated Population Models with application to Skagit River 
#       Chinook Recovery Evaluation - Near-Term Action (NTA) 2018-0697
# Authors: Sarah L. Schooler, Michael LeMoine
# Collaborators: Casey Ruff, Eric Beamer, Catherine Austin
# Code Prepared By: Sarah L. Schooler, Casey Ruff
# Contact: sschooler@skagitcoop.org
# 17 March 2023
######################################################################
######################################################################

# load necessary packages
if(!require("here")) {
  install.packages("here")
  library("here")
}
if(!require("rjags")) {
  install.packages("rjags")
  library("rjags")
}
if(!require("jagsUI")) {
  install.packages("jagsUI")
  library("jagsUI")
}
if(!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
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
## first & last years of fish data
# change these for simulations 15-19 for numbers of years
yr_frst <- 1992
yr_last <- 2021

## min & max adult age classes
age_min <- 2
age_max <- 5
# range of ages
A <- age_max - age_min + 1

## load data frame
sim.num <- 19 # change to apply each simulation (note: 15-20 use same file, change yrs)
dat <- read.csv(file.path(datadir, paste0("simulations/Sim", sim.num, "_Data.csv")))
dat <- dat[which(dat$year %in% yr_frst:yr_last),]
dat_yrs <- dat$year; n_yrs <- length(dat_yrs)
# define escapement
dat_esc <- dat$escapement
ln_dat_esc <- log(dat_esc) ## log of escapement
# define subyearling smolt abundance
dat_smolt <- dat$smolts
ln_dat_smolt <- log(dat_smolt) ## log of subyearling smolt abundance
# define age comp data & drop first age_min rows
dat_age <- dat[, grepl("age", colnames(dat))][-(1:(age_min)),]
# get total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age, 1, sum)
# convert class
dat_age <- as.matrix(dat_age)
# define harvest data
dat_harv <- dat$harvest

##########################
### define jags model ####
##########################

cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## Note: the term E means expected value of the true state ###
  
  ## freshwater (adult-smolt) productivity (pas)
  pas <- exp(mu_pas)
  mu_pas ~ dunif(0, 8)

  ## marine (smolt-adult) productivity (psa)
  psa <- exp(mu_R_psa)
  mu_R_psa ~ dunif(-10, 0)

  ## marine density dependence (beta)
  beta ~ dunif(0, .1)
  
  ## observation variance for spawners
  sigma_sp ~ dunif(0, 10)
  tau_sp <- 1/sigma_sp;

  ## process variance for smolts
  sigma_pas ~ dunif(0, 10)
  tau_pas <- 1/sigma_pas;
  
  ## process variance for spawners
  sigma_psa ~ dunif(0, 10)
  tau_psa <- 1/sigma_psa
  
  ## observation variance for smolts
  sigma_sm ~ dunif(0, 10) 
  tau_sm <- 1/sigma_sm;

  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0, 0.001) 
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
	ttl_run_mu ~ dunif(1,5);
	ttl_run_tau ~ dunif(1,20);

  ## get total cal yr returns for first age_min yrs
  for(i in 1:age_min) {
		ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
		tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## maturity schedule
  ## uniform vector for Dirichlet prior
  theta <- c(.2,.2,.2,.2)
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dnorm(0, 0.01) T(0,);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ## estimated harvest rate
  for(t in 1:(n_yrs)) { h_rate[t] ~ dunif(0,1) }
  av_hrate <- mean(h_rate)

  ##------------
  ## LIKELIHOOD
  ##------------
  ## Marine productivity (smolt-adult) model as ricker density dependence
  
  ## brood years 1:(n_yrs-age_min)
  for(t in 1:(n_yrs-age_min)) {
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
    
    ln_R_psa[t] ~ dnorm(mu_R_psa, tau_psa) ;
    E_ln_Rec[t] <- E_ln_smolt[t] + ln_R_psa[t] - beta*smolt[t];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t], tau_psa)
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
  
    # for posterior check
    tot_ln_Rec_new[t] ~  dnorm(E_ln_Rec[t], tau_psa)
    res_ln_Rec_new[t] <-  tot_ln_Rec_new[t] - E_ln_Rec[t];

  } ## end t loop over year
  
  # sums for Bayesian p-value
  fit_lnRec <- sum(res_ln_Rec[])
  fit.new_lnRec <- sum(res_ln_Rec_new[])
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min
  for(i in 1:(age_max-age_min)) {
    ## projected recruits
    for(a in 1:(i)) {
      Run[i,a] <- Rec[(i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinomial
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A):(n_yrs-age_min)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## freshwater productivity
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs)) {
    # estimating spawners
    est_harv[t] = h_rate[t] * tot_Run[t];
    dat_harv[t] ~ dlnorm(log(est_harv[t]), 20);
    Sp[t] = tot_Run[t] - est_harv[t];
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_sp);

    ## process model for smolts with normally distributed error
    ln_pas[t] ~ dnorm(mu_pas, tau_pas); # process error
    E_ln_smolt[t] <- ln_pas[t] + ln_Sp[t] 
    smolt[t] <- exp(E_ln_smolt[t])
    
    ## observation  model for smolt production
    ln_dat_smolt[t] ~ dnorm(E_ln_smolt[t], tau_sm);
    
    # model predictive check for Bayesian p-value
    res_ln_smolt[t] <- ln_dat_smolt[t] - E_ln_smolt[t];
    ln_dat_smolt_new[t] ~  dnorm(E_ln_smolt[t], tau_sm)
    res_ln_smolt_new[t] <-ln_dat_smolt_new[t] - E_ln_smolt[t];
    
  }
  
  # Summary for fit statistics
  fit_lnSmolt <- sum(res_ln_smolt[])
  fit.new_lnSmolt <- sum(res_ln_smolt_new[])
    
  
} ## end model description

", file=file.path(jagsdir, "Skagit_Chinook_TwoStage_IPM_sims.txt"))

##########################
####### JAGS setup #######
##########################
## data to pass to JAGS:
dat_jags <- list(dat_age = dat_age, ln_dat_esc = ln_dat_esc, 
                 ln_dat_smolt = ln_dat_smolt, dat_harv = dat_harv,
                 A = A, age_min = age_min, n_yrs = n_yrs, 
                 age_max = age_max) 

## define initial values 
init_vals <- function() {
  list(pi_tau = 10, pi_eta = rep(1,A),
       pi_vec = matrix(c(0.02,0.20,0.50,0.18), n_yrs-age_min, A, byrow = TRUE),
       tot_ln_Rec = rep(log(1000), n_yrs - age_min))
}

## model parameters/states for JAGS to return:
par_jags <- c("pas","mu_pas","psa","mu_R_psa", "beta",
              "av_hrate", "sigma_sp","sigma_sm","sigma_pas",
              "Sp","smolt","Rec","tot_ln_Rec",
              "pi_eta","pi_tau", "res_ln_Rec",
              "fit_lnRec", "fit.new_lnRec",
              "fit_lnSmolt", "fit.new_lnSmolt")

# 4 chains, 50000 adaptations (similar to burn-in),
# thin rate of 400, and 400000 iterations
# note: for simulations with > 25 years, number of iterations
# and adaptations must be increased for convergence
# may need to add burn-in
nc <- 4
ni <- 400000
na <- 50000
nt <- 400
nb <- 0

##########################
####### fit model ########
##########################
# fit model through jagsUI using parallel processing
# where n.cores = n.chains
mod_fit <- jags(data = dat_jags, inits = init_vals,
                model.file = file.path(jagsdir, "Skagit_Chinook_TwoStage_IPM_sims.txt"),
                parameters.to.save = par_jags, parallel = T, 
                n.chains = nc, n.adapt = na, n.thin = nt, n.iter = ni,
                n.burnin = nb)
beep()
# examine results
summary(mod_fit)
head(formatC(mod_fit$summary, format = "e", digits = 2), n = 10)
jags.View(mod_fit)

##########################
#### model diagnostics ###
##########################

# examine traceplots and density plots for convergence
## params of interest
par_conv <- c("pas","mu_pas","psa","E_R_psa","mu_R_psa",
              "beta", "sigma_r","sigma_sp","sigma_sm", "sigma_pas")
traceplot(mod_fit, parameters = par_conv)
densityplot(mod_fit, parameters = par_conv)

# examine Bayesian p-values
pp.check(mod_fit, observed = "fit_lnRec", simulated = "fit.new_lnRec")
pp.check(mod_fit, observed = "fit_lnSmolt", simulated = "fit.new_lnSmolt")

## examine autocorrelation
t(round(autocorr.diag(
  mod_fit$samples[,c("pas","mu_pas","psa","mu_R_psa",
                     "beta", "sigma_sp","sigma_sm", "sigma_pas")],
  lags = seq(mod_fit$mcmc.info$n.thin,  4*mod_fit$mcmc.info$n.thin, 
             mod_fit$mcmc.info$n.thin),  relative = FALSE), 2))


##########################
####### make plots #######
##########################

model_results <- mod_fit$summary

# Escapement estimates
as.data.frame(model_results[grep("Sp", rownames(model_results)),]) %>%
  mutate("median" = `50%`, "lower" = `2.5%`, 
         "upper" = `97.5%`, "data" = dat_esc,
         year = seq(yr_frst, length.out = n_yrs), .keep = "none") %>%
  
  ggplot(aes(x = year, y = median/1000)) +
    geom_ribbon(aes(ymin = lower/1000, ymax = upper/1000), alpha = .2, fill = "blue") +
    geom_line(color = "blue", linewidth = 1) + theme_classic() + 
    labs(x = "Year", y = "Escapement (thousands)") + 
    scale_x_continuous(breaks = seq(1980, 2015, 5), expand = c(0, 0)) +
    geom_point(aes(y = data/1000))

as.data.frame(model_results[grep("smolt", rownames(model_results)),]) %>%
  mutate("median" = `50%`, "lower" = `2.5%`, 
         "upper" = `97.5%`, "data" = dat_smolt,
         year = seq(yr_frst, length.out = n_yrs), .keep = "none") %>%
  
  ggplot(aes(x = year, y = median/1e6)) +
    geom_ribbon(aes(ymin = lower/1e6, ymax = upper/1e6), alpha = .2, fill = "blue") +
    geom_line(color = "blue", linewidth = 1) + theme_classic() + 
    labs(x = "Year", y = "Smolt (millions)") + 
    scale_x_continuous(breaks = seq(1980, 2015, 5), expand = c(0, 0)) +
    geom_point(aes(y = data/1e6))


