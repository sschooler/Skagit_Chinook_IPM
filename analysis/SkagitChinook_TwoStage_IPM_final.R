######################################################################
######################################################################
# State-space life-cycle Bayesian integrated population model (IPM)
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


## set directory locations
datadir <- here("data")
jagsdir <- here("jags")
analdir <- here("analysis")
savedir <- here("analysis/cache")

##########################
#### initialize data #####
##########################
## first & last years of fish data
yr_frst <- 1993
yr_last <- 2016

## min & max adult age classes
age_min <- 2
age_max <- 5
# range of ages
A <- age_max - age_min + 1

## escapement data
dat_esc <- read.csv(file.path(datadir,  "skagit_chinook_esc.csv"))
dat_esc <- dat_esc[which(dat_esc$year %in% seq(yr_frst,yr_last,1)),]
# years of data
dat_yrs <- dat_esc$year
# number of years of data
n_yrs <- length(dat_yrs)
# log of escapement
ln_dat_esc <- log(dat_esc$escapement)


## subyearling smolt abundance
dat_smolt <- read.csv(file.path(datadir, "skagit_chinook_smolt.csv"))
dat_smolt <- dat_smolt[which(dat_smolt$year %in% seq(yr_frst,yr_last,1)),]
ln_dat_smolt <- log(dat_smolt$smolt)


## age composition
dat_age <- read.csv(file.path(datadir,  "skagit_chinook_age.csv"))
dat_age <- dat_age[which(dat_age$year %in% seq(yr_frst,yr_last,1)),]
# drop year col & first (age_min ) rows
dat_age <- dat_age[-(1:(age_min)),-1]
# num of age classes
A <- age_max - age_min + 1
# total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age, 1, sum)
# row indices for any years with no obs age comp
idx_NA_yrs <- which(dat_age$sum<A, TRUE)
# replace 0's in yrs w/o any obs with NA's
dat_age[idx_NA_yrs,(1:A)] <- NA
# change total in yrs w/o any obs from 0 to A to help dmulti()
dat_age[idx_NA_yrs,"sum"] <- A
# convert class
dat_age <- as.matrix(dat_age)

## harvest
dat_harv <- read.csv(file.path(datadir,"skagit_chinook_aeq_harvest.csv"))
dat_harv <- dat_harv[which(dat_harv$year %in% seq(yr_frst,yr_last,1)),2]

## covariate(s)
dat_cvrs <- read.csv(file.path(datadir, "skagit_chinook_covars.csv"))
dat_cvrs <- dat_cvrs[dat_cvrs$year>=yr_frst & dat_cvrs$year<=yr_last,]
## drop year col
dat_cvrs <- dat_cvrs[,-1] 
scl_cvrs <- as.matrix(scale(dat_cvrs)) 
# separate adult-smolt productivity and smolt-adult productivity covariates
# adult-smolt prod. = "as_"; smolt-adult = "sa_"
scl_cvrs_sa <- scl_cvrs[,grep("sa", colnames(scl_cvrs))]
scl_cvrs_as <- scl_cvrs[,grep("as", colnames(scl_cvrs))]
# extract covariates from nls model selection
covs_names_sa <- c("sa_NPGO","sa_CUTI", "sa_MEI", "sa_SST")
covs_names_as <- c("as_floodRI", "as_propRI1", "as_flowFebJun", "as_airTJanApr")
# get covariate names
scl_cvrs_sa <- scl_cvrs_sa[, covs_names_sa]
scl_cvrs_as <- scl_cvrs_as[, covs_names_as]
## total number of covariates
n_cov_sa <- dim(scl_cvrs_sa)[2] 
n_cov_as <- dim(scl_cvrs_as)[2]

##########################
### define jags model ####
##########################

cat("

model {
  
  ##--------
  ## PRIORS
  ##--------

  ## freshwater (adult-smolt) productivity (pas)
  pas <- exp(mu_pas)
  mu_pas ~ dunif(0, 8)

  ## marine (smolt-adult) productivity (psa)
  psa <- exp(mu_R_psa)
  mu_R_psa ~ dunif(-10, 0)
  E_R_psa <- mu_R_psa + sigma_r/(2 - 2*phi^2)
  
  ## marine density dependence (beta)
  beta ~ dunif(0, .1)
  
  # covariates for each stage (gamma)
  for(i in 1:n_cov_sa) { 
    gamma_sa[i] ~ dnorm(0,0.01) 
  }
  for(i in 1:n_cov_as){
    gamma_as[i] ~ dnorm(0, 0.01)
  }
  
  ## auto-regressive coef for proc errors
  phi ~ dnorm(.1,1e-6)I(-1,1)
  
  ## process variance for adult recruits model (feeds into autogregression)
  sigma_r ~ dunif(0, 10) 
  tau_r <- 1/sigma_r;

  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## observation variance for spawners
  sigma_sp ~ dunif(0, 10)
  tau_sp <- 1/sigma_sp;

  ## process variance for smolts
  sigma_pas ~ dunif(0, 10)
  tau_pas <- 1/sigma_pas;
  
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
  for(i in 1:(age_min)) {
		ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
		tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## maturity schedule
  ## unif vec for Dirch prior
  theta <- c(.2,.2,.2,.2)
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dnorm(0, 0.01) T(0,);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ## estimated harvest rate
  for(t in 1:(n_yrs)) { h_rate[t] ~ dunif(0,1) }

  ##------------
  ## LIKELIHOOD
  ##------------
  ## Marine productivity (smolt-adult) model as ricker density dependence
  ## with autoregressive error structure
  
  ## 1st brood yr requires different innovation
  covar_sa[1] <- inprod(gamma_sa,mod_cvrs_sa[1,]);
  ln_R_psa[1] <- mu_R_psa + covar_sa[1]; 
  E_ln_Rec[1] <- E_ln_smolt[1] + ln_R_psa[1] - beta*smolt[1] + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  
  # for posterior check
  tot_ln_Rec_new[1] ~  dnorm(E_ln_Rec[1],tau_r)
  res_ln_Rec_new[1] <-  tot_ln_Rec_new[1] - E_ln_Rec[1];
  
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## brood-yr recruits by age
  for(a in 1:A) {
    Rec[1,a] <- tot_Rec[1] * pi_vec[1,a];
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min)) {
    covar_sa[t] <- inprod(gamma_sa, mod_cvrs_sa[t,])
    ln_R_psa[t] <- mu_R_psa + covar_sa[t]
    E_ln_Rec[t] <- E_ln_smolt[t] + ln_R_psa[t] - beta*smolt[t] + phi*res_ln_Rec[t-1]
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r)
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t]
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
  
    # model predictive check for Bayesian p-value
    tot_ln_Rec_new[t] ~  dnorm(E_ln_Rec[t],tau_r)
    res_ln_Rec_new[t] <-  tot_ln_Rec_new[t] - E_ln_Rec[t];
  
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- tot_Rec[t] * pi_vec[t,a];
    }
  } ## end t loop over year
  
  # sums for Bayesian p-value
  fit_lnRec <- sum(res_ln_Rec[])
  fit.new_lnRec <- sum(res_ln_Rec_new[])
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
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
  
  ## step 2: info from complete broods
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
    covar_as[t] <- inprod(gamma_as, mod_cvrs_as[t,]);
    ln_pas[t] ~ dnorm(mu_pas, tau_pas); # process error
    E_ln_smolt[t] <- ln_pas[t] + covar_as[t] + ln_Sp[t] 
    smolt[t] <- exp(E_ln_smolt[t])
    
    ## observation  model for smolt production
    ln_dat_smolt[t] ~ dnorm(E_ln_smolt[t], tau_sm);
    
    # model predictive check for Bayesian p-value
    res_ln_smolt[t] <- ln_dat_smolt[t] - E_ln_smolt[t];
    ln_dat_smolt_new[t] ~  dnorm(E_ln_smolt[t],tau_sm)
    res_ln_smolt_new[t] <-ln_dat_smolt_new[t] - E_ln_smolt[t];
    
  }
  
  # Summary for fit statistics
  fit_lnSmolt <- sum(res_ln_smolt[])
  fit.new_lnSmolt <- sum(res_ln_smolt_new[])
    
  
} ## end model description

", file=file.path(jagsdir, "Skagit_Chinook_TwoStage_IPM.txt"))


##########################
####### JAGS setup #######
##########################
## data to pass to JAGS:
dat_Rcov_jags <- list(dat_age = dat_age,
                      ln_dat_esc = ln_dat_esc, ln_dat_smolt = ln_dat_smolt,
                      dat_harv = dat_harv, A = A, age_min = age_min,
                      age_max = age_max, n_yrs = n_yrs, 
                      mod_cvrs_sa = scl_cvrs_sa, n_cov_sa = n_cov_sa,
                      mod_cvrs_as = scl_cvrs_as, n_cov_as = n_cov_as) 

## model parameters/states for JAGS to return:
par_Rcov_jags <- c("pas","mu_pas","psa","E_R_psa","mu_R_psa", "beta",
                   "sigma_r","sigma_sp","sigma_sm","sigma_pas",
                   "gamma_sa", "gamma_as", 
                   "Sp","smolt","Rec","tot_ln_Rec",
                   "pi_eta","pi_tau", "res_ln_Rec",
                   "fit_lnRec", "fit.new_lnRec",
                   "fit_lnSmolt", "fit.new_lnSmolt")

## define initial values
init_vals_AR <- function() {
  list(mu_pas = log(250), mu_psa=log(.01), beta = 1.5e-7, 
       pi_tau = 10, pi_eta = rep(1,A),
       pi_vec = matrix(c(0.02,0.20,0.50,0.18), n_yrs-age_min, A, 
                       byrow = TRUE),
       Rec_mu = log(1000), Rec_sig = 0.1, 
       tot_ln_Rec = rep(log(1000), n_yrs - age_min),
       innov_1 = 0, phi = 0.5)
}

# 4 chains, a burn-in (n.adapt) of 50000,
# thin rate of 400, and 400000 iterations
nc <- 4
ni <- 400000
na <- 50000
nt <- 400

##########################
####### fit model ########
##########################

# fit model through jagsUI using parallel processing
# with 4 cores; 2.30 GHz processor takes ~ 7 minutes
mod_fit <- jags(data = dat_Rcov_jags, inits = init_vals_AR,
                model.file = file.path(jagsdir, "Skagit_Chinook_TwoStage_IPM.txt"),
                parameters.to.save = par_Rcov_jags, parallel = T, 
                n.chains = nc, n.adapt = na,
                n.thin = nt, n.iter = ni)
# examine model results and check convergence
print(mod_fit)
model_results <- mod_fit$summary
# examine results
head(formatC(mod_fit$summary, format = "e", digits = 2), n = 20)
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
  mod_fit$samples[,c("pas","mu_pas","psa","E_R_psa","mu_R_psa",
                     "beta", "sigma_r","sigma_sp","sigma_sm", "sigma_pas")],
  lags = seq(mod_fit$mcmc.info$n.thin,  4*mod_fit$mcmc.info$n.thin, 
             mod_fit$mcmc.info$n.thin),  relative = FALSE), 2))
