######################################################################
######################################################################
# Plots comparing IPM and NLS results
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
if(!require("nlstools")) {
  install.packages("nlstools")
  library("nlstools")
}
if(!require("tidyr")) {
  install.packages("tidyr")
  library("tidyr")
}
if(!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
if(!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}
if(!require("ggh4x")) {
  install.packages("ggh4x")
  library("ggh4x")
}
if(!require("stringr")) {
  install.packages("stringr")
  library("stringr")
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

###############################
#### load models and data #####
###############################
lapply(file.path(savedir, list.files(savedir)), load, .GlobalEnv)
# read in smolts, spawners, recruitment estimates, and covariates
smolts <- read.csv(file.path(datadir, "skagit_chinook_smolt.csv"))
spawners <- read.csv(file.path(datadir, "skagit_chinook_esc.csv"))
recruits <- read.csv(file.path(datadir, "skagit_chinook_rec_est.csv"))
dat_cvrs <- read.csv("./data/skagit_chinook_covars.csv")
dat_cvrs <- dat_cvrs[dat_cvrs$year>1992&dat_cvrs$year<2015, ]
yr_frst <- min(dat_cvrs$year)
n_yrs <- nrow(dat_cvrs)

###############################
## posterior plot data frames #
###############################

## Bayesian IPM parameter posteriors
all_samples.df <- data.frame(do.call("rbind", mod_fit$samples))
# subset columns to plot
plot_samples.df <- all_samples.df[, c("pas", "psa", "beta", 
                                      grep("gamma", colnames(all_samples.df), 
                                           value = T))]
# rename columns
names(plot_samples.df)[grep("gamma", colnames(plot_samples.df))] <- 
  gsub("\\.", "", names(plot_samples.df)[grep("gamma", colnames(plot_samples.df))])

## Bootstrapped NLS posteriors
# smolt to adult (marine) 
boot_psa <- nlsBoot(psa_nls, niter = 1000)
nls_sa_results <- data.frame(psa = boot_psa$coefboot[,"a"], 
                             beta = boot_psa$coefboot[,"b"],
                             boot_psa$coefboot[,grep("g", colnames(boot_psa$coefboot))])
colnames(nls_sa_results)[grep("g", colnames(nls_sa_results))] <-
  paste0(rep("gamma_sa", length(grep("g", colnames(nls_sa_results)))), 
         c(1:length(grep("g", colnames(nls_sa_results)))))
# adult to smolt (freshwater)
boot_pas <- nlsBoot(pas_nls, niter = 1000); #confint(boot_pas, plot = T)
nls_as_results <- data.frame(pas = boot_pas$coefboot[,"a"],
                             boot_pas$coefboot[,grep("g", colnames(boot_pas$coefboot))])
colnames(nls_as_results)[grep("g", colnames(nls_as_results))] <-
  paste0(rep("gamma_as", length(grep("g", colnames(nls_as_results)))), 
         c(1:length(grep("g", colnames(nls_as_results)))))

## combine into a single data frame of posteiors
all_results.df <- rbind(nls_as_results %>% pivot_longer(everything()) %>% 
                          mutate(model = "NLS", stage = "AS"),
                        nls_sa_results %>% pivot_longer(everything()) %>% 
                          mutate(model = "NLS", stage = "SA"),
                        plot_samples.df %>% pivot_longer(everything()) %>% 
                          mutate(model = "IPM", stage = "NA"))
all_results.df[all_results.df$model == "IPM" &
                 grepl("as", pull(all_results.df[,"name"])),"stage"] <- "AS"
all_results.df[all_results.df$model == "IPM" &
                 (grepl("sa", pull(all_results.df[,"name"])) |
                    all_results.df$name == "beta"),"stage"] <- "SA"
all_results.df$name <- factor(all_results.df$name, 
                              levels = c("pas", "psa", "beta", 
                                         grep("gamma", levels(factor(all_results.df$name)), 
                                              value = T)))
# add names of covariates for plot labels
# note: ensure the order is correct and each gamma matches
# with the correct covariate
all_results.df$plot_names <- recode(factor(all_results.df$name), 
                                    beta = "Density dependence",
                                    pas = "Freshwater productivity", 
                                    psa = "Marine productivity",
                                    gamma_as1 = "Flood RI",
                                    gamma_as2 = "Prop. above 1-yr RI",
                                    gamma_as3 = "Mean spring flow",
                                    gamma_as4 = "Mean spring air temp.",
                                    gamma_sa1 = "NPGO",
                                    gamma_sa2 = "CUTI",
                                    gamma_sa3 = "MEI",
                                    gamma_sa4 = "Sea surface temp.")


###############################
##### calculate quantiles #####
###############################
# productivity posterior quantiles
prodest_quantiles.df <-
  all_results.df[!grepl("gamma", all_results.df$name)
                 ,] %>%
  group_by(model, stage, plot_names) %>%
  summarise(
    "lower" = quantile(value, 0.025, na.rm = T),
    "upper" = quantile(value, 0.975, na.rm = T),
    "median" = quantile(value, .5, na.rm = T)
  ) 
# remove extreme values
all_results.df[all_results.df$name == "psa"&all_results.df$value>0.025, "value"] <- NA

# coefficient posterior quantiles
coef_quantiles.df <-
  all_results.df[grep("gamma", all_results.df$name),] %>%
  group_by(model, stage, plot_names) %>%
  summarise(
    "lower" = quantile(value, 0.025),
    "upper" = quantile(value, 0.975),
    "median" = quantile(value, .5)
  ) 

###############################
####### posterior plots #######
###############################
# plot productivity parameter estimates with 97.5% confidence/credible
# intervals and a line at the median
parameterests.gg <- 
  ggplot(data = all_results.df[!grepl("gamma", all_results.df$name), ], 
         aes(x = value, group = plot_names)) + 
  geom_rect(data = prodest_quantiles.df,
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf), fill = "lightblue1",
            inherit.aes = F) +
  geom_histogram(aes(y = after_stat(density)) ,color = "blue3", fill = "blue3", alpha = .3) +
  geom_vline(data = prodest_quantiles.df, aes(xintercept = median), color = "darkblue",
             linewidth = 1) +
  theme_classic() + facet_grid2(model~plot_names, scales = "free", switch = "both",
                                independent = "y") +
  labs(x = "Estimates") + #xlim(c(-.8, .7)) + 
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank(), plot.title = element_text(hjust = .5),
        strip.placement = "outside", strip.background.x = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"), 
        strip.text = element_text(size = 9), panel.spacing.x = unit(8, "pt"))
# examine plot
parameterests.gg


# marine (smolt-adult [SA]) covariate parameter estimates with 97.5% 
# confidence/credible intervals and a line at the median
SmoltAdultcovs.gg <- ggplot(data = all_results.df[all_results.df$stage == "SA" &
                                                    grepl("gamma", all_results.df$name), ], 
                            aes(x = value, group = plot_names)) + 
  geom_rect(data = coef_quantiles.df[coef_quantiles.df$stage == "SA",],
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf), fill = "lightblue1",
            inherit.aes = F) +  
  geom_vline(data = coef_quantiles.df[coef_quantiles.df$stage == "SA",], 
             aes(xintercept = median), color = "darkblue", linewidth = 1) +
  #geom_vline(xintercept = 0) +
  geom_histogram(aes(y = after_stat(density)) ,color = "blue3", fill = "blue3", alpha = .3) +
  theme_classic() + 
  facet_grid2(model~plot_names, scales = "free", switch = "y", independent = "y") +
  labs(x = "Coefficient estimate") + #xlim(c(-.8, .7)) + 
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9), plot.title = element_text(hjust = .5),
        axis.text.x = element_text(size = 8, color = "black"), strip.text = element_text(size = 9),
        panel.spacing.x = unit(9, "pt"))
# examine plot
SmoltAdultcovs.gg

# freshwater (adult-smolt [AS]) covariate parameter estimates with 97.5% 
# confidence/credible intervals and a line at the median
AdultSmoltcovs.gg <- ggplot(data = all_results.df[all_results.df$stage == "AS"&
                                                    grepl("gamma", all_results.df$name), ], 
                            aes(x = value, group = plot_names)) + 
  geom_rect(data = coef_quantiles.df[coef_quantiles.df$stage=="AS",],
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf), fill = "lightblue1",
            inherit.aes = F) +  
  geom_vline(data = coef_quantiles.df[coef_quantiles.df$stage=="AS",], 
             aes(xintercept = median), color = "darkblue", linewidth = 1) +
  geom_histogram(aes(y = after_stat(density)) ,color = "blue3", fill = "blue3", alpha = .3) +
  theme_classic() + facet_grid2(model~plot_names, scales = "free", switch = "y", independent = "y") +
  labs(x = "Coefficient estimate") + #xlim(c(-.6, .6)) + 
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), strip.text = element_text(size = 9))

# examine plot
AdultSmoltcovs.gg


###############################
# population est. data frames #
###############################
# Bayesian estimates
fullEst_Bayes <- as.data.frame(mod_fit$summary)
# NLS bootstraps
PSAest_NLS <- as.data.frame(exp(nlsBootPredict(boot_psa, interval = "confidence")))
PASest_NLS <- as.data.frame(exp(nlsBootPredict(boot_pas, interval = "confidence")))

# combines into a single data frame with some modifications
# exponentiates Bayesian recruitment estimates
# translates to thousands (recruitment) and millions (smolts)

popests.df <- 
  rbind(fullEst_Bayes[grep("tot_ln_Rec", rownames(fullEst_Bayes)),] %>% 
          mutate(lower = exp(`2.5%`)/1000, upper = exp(`97.5%`)/1000, median = exp(`50%`)/1000,
                 Year = seq(yr_frst, length.out = n_yrs), data = psa.df$recruits/1000, model = "IPM",
                 type = "Recruitment (thousands)",
                 .keep = "none"),
        fullEst_Bayes[grep("smolt", rownames(fullEst_Bayes)),] %>% 
          mutate(lower = `2.5%`/1E6, upper = `97.5%`/1E6, median = `50%`/1E6,
                 Year = seq(yr_frst, length.out = n_yrs+2), data = c(pas.df$recruits/1E6, NA, NA), model = "IPM",
                 type = "Smolt (millions)",
                 .keep = "none"),
        PSAest_NLS %>% 
          mutate(lower = `2.5%`/1000, upper = `97.5%`/1000, median = Median/1000,
                 Year = seq(yr_frst, length.out = n_yrs), data = psa.df$recruits/1000, model = "NLS",
                 type = "Recruitment (thousands)",
                 .keep = "none"),
        PASest_NLS %>% 
          mutate(lower = `2.5%`/1E6, upper = `97.5%`/1E6, median = Median/1E6,
                 Year = seq(yr_frst, length.out = n_yrs), data = pas.df$recruits/1E6, model = "NLS",
                 type = "Smolt (millions)",
                 .keep = "none"))

# creates plot
popestimates.gg <- 
  ggplot(data = popests.df, aes(x = Year, y = median)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_line(color = "blue", linewidth = 1) + geom_point(aes(y = data)) +
  facet_grid2(type~model, switch = "y", scales = "free_y", axes = "all", remove_labels = "all") + 
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(0, 45), expand = c(0,0)),
                               scale_y_continuous(limits = c(0, 9.9), expand = c(0,0)))) +
  theme_classic() + 
  theme(axis.title.y = element_blank(), strip.placement = "outside", 
        strip.background.y = element_blank(), 
        axis.text = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(size = 9),
        strip.text = element_text(size = 9)) + 
  scale_x_continuous(limits = c(1992, 2014), expand = c(0,0))
popestimates.gg
