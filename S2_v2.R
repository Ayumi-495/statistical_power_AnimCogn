pacman::p_load(
  # data wrangling
  tidyverse, dplyr, tidyr, purrr, stringr,
  # visualisation
  ggplot2, gghighlight, ggrepel,patchwork,
  # meta-analysis
  metafor, esc, orchaRd,
  # other
  here, janitor, kableExtra, readr, gt, lme4
  )

# Custom functions from Yang et al(2023) ----

# Function to calculate power (two-tail) for meta-analysis
power.ma_Shinichi <- function(mu, SE, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/SE)-pnorm(qnorm(1-alpha/2)+abs(mu)/SE)
} # or power.ma_Shinichi1 <- function(mu,SE){1 - pnorm(qnorm(1-0.05/2)-abs(mu)/SE) + pnorm(-qnorm(1-0.05/2)-abs(mu)/SE)}


# Function for power analysis for empirical data point
power.individual_Shinichi <- function(mu, se, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/se)-pnorm(qnorm(1-alpha/2)+abs(mu)/se)} # two-tailed power


# Function for Type S error for empirical data point
error_S <- function(mu, se, alpha = 0.05){
  #z <- qnorm(1 - alpha/2) # Z-score or quantile
  p.u <- 1 - pnorm(qnorm(1 - alpha/2) - abs(mu)/se) # upper-tail probability
  p.l <- pnorm(-qnorm(1 - alpha/2) - abs(mu)/se) # lower-tail probability
  power <- p.u + p.l # upper + lower
  errorS <- p.l/power # percentage of the opposite direction
  return(errorS)
} 

# Function for Type M error for empirical data point
error_M <- function(mu, se, alpha = 0.05, N = 10000) {
  est.random <- rnorm(n=N, mean = mu, sd = se)
  # est.random <- mu + se*rnorm(n=N, mean=0, sd=1)
  sig.index <- abs(est.random) > se*qnorm(1 - alpha/2)
  overestimate <- mean(abs(est.random)[sig.index])/abs(mu) # ratio is regardless of sign, so we need absolute value
  absolute_error <- overestimate*abs(mu) - abs(mu)
  relative_error <- absolute_error/(overestimate*abs(mu))
  return(abs(overestimate) |> round(3))
}


error_M2 <- function(mu, se, alpha = 0.05, N = 10000) {
  est.random <- rnorm(n=N, mean = mu, sd = se)
  # est.random <- mu + se*rnorm(n=N, mean=0, sd=1)
  sig.index <- abs(est.random) > se*qnorm(1 - alpha/2)
  overestimate <- mean(abs(est.random)[sig.index])/abs(mu) # ratio is regardless of sign, so we need absolute value
  absolute_error <- overestimate*abs(mu) - abs(mu)
  relative_error <- absolute_error/(overestimate*abs(mu))
  return(abs(relative_error) |> round(3))
} # relative error: (M - 1) / M


# meta-analysis of magnitude
## folded effect size
folded_es <-function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
}
## folded error
folded_error <- function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <- fold_se^2
  fold_v
}


# custom function for extracting mean and CI from each metafor model
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}

###
# 1. Load the data and calculate effect sizes ----
###

## read lnRR, SMD and Zr datasets, which contain the calculated effect size (es) and sampling variance (var)
## need to use full name of each dataset, otherwise read_csv is not able to read it

### lnRR
lnRR_csv <- list.files(path = "lnRR", pattern = "*.csv", full.names = TRUE) |> lapply(read_csv)
### SMD
SMD_csv <- list.files(path = "SMD", pattern = "*.csv", full.names = TRUE) |> lapply(read_csv)
### Zr
Zr_csv <- list.files(path = "Zr", pattern = "*.csv", full.names = TRUE) |> lapply(read_csv)

## get names of each .csv file
lnRR_filenames <- list.files(path = "lnRR", pattern = "*.csv", full.names = FALSE)
SMD_filenames <- list.files(path = "SMD", pattern = "*.csv", full.names = FALSE)
Zr_filenames <- list.files(path = "Zr", pattern = "*.csv", full.names = FALSE)

### rename the elements of the list
names(lnRR_csv) <- lnRR_filenames
names(SMD_csv) <- SMD_filenames
names(Zr_csv) <- Zr_filenames

## calculate the effect size (es) and sampling variance (var) from raw data
### get file names and paths
SMD_des_csv <- list.files(path = "./SMD/des_stat", pattern = "*.csv", full.names = TRUE) |> 
  lapply(read_csv)
SMD_des_filenames <- list.files(path = "./SMD/des_stat", pattern = "*.csv", full.names = FALSE) # extract file names, which will be used later

## for datasets with descriptive statistics, we also need to create 'effective sample size'-related variables for each of them 
## these variables make the examination of the small-study effect more statistically sound

### function to calculate effective sample size
ess_cal <- function(dat){(4*dat$C_n*dat$T_n) / (dat$C_n + dat$T_n)}

### calculate effective sample size for SMD with descriptive statistics
ess <- NA
for (i in 1:length(SMD_des_csv)) {
  ess[i] <- ess_cal(SMD_des_csv[[i]]) |> list()}

### allocate each set of effective sample size into corresponding dataset
for (i in 1:length(SMD_des_csv)) {
  SMD_des_csv[[i]]$ess <- ess[[i]]
}

#### create inverse of effective sample size - "effective sample size" based "sampling variance" 
#### function to calculate inverse of effective sample size 
ess.var_cal <- function(dat){1/dat$C_n + 1/dat$T_n}

#### calculations for SMD
ess.var <- NA
for (i in 1:length(SMD_des_csv)) {
  ess.var[i] <- ess.var_cal(SMD_des_csv[[i]]) |> list()}
#### allocate each set of effective sample size into corresponding dataset
for (i in 1:length(SMD_des_csv)) {
  SMD_des_csv[[i]]$ess.var <- ess.var[[i]]
}
#### calculate inverse sqrt of effective sample size - "effective sample size" based "sampling error" 
for (i in 1:length(SMD_des_csv)) {
  SMD_des_csv[[i]]$ess.sei <- sqrt(SMD_des_csv[[i]]$ess.var)
}

SMD_es <- NA
for (i in 1:length(SMD_des_csv)) {
  SMD_es[i] <- escalc(measure = "SMD",
                      m1i = T_mean,
                      m2i = C_mean,
                      sd1i = T_sd,
                      sd2i = C_sd,
                      n1i = T_n,
                      n2i = C_n,
                      data = SMD_des_csv[[i]]) |> list()
}
### rename the elements of the list
names(SMD_es) <- SMD_des_filenames

### SMD
for (i in 1:length(SMD_es)) {
  names(SMD_es[[i]])[names(SMD_es[[i]]) == "yi"] <- "es"
  names(SMD_es[[i]])[names(SMD_es[[i]]) == "vi"] <- "var"
}

## combine two sets of dataset for SMD 
SMD <- append(SMD_csv, SMD_es)
lnRR <- lnRR_csv # for consistence, create lnRR to instead of lnRR_csv
Zr <- Zr_csv # for consistence, create Zr to instead of Zr_csv

## remove NAs, zero variance…
### lnRR
#### delete NAs, zero variance
for (i in 1:length(lnRR)) {
  lnRR[[i]] <- lnRR[[i]][!is.na(lnRR[[i]]$es) & !is.na(lnRR[[i]]$var) & lnRR[[i]]$var != 0 & !is.na(lnRR[[i]]$year), ]
}

### SMD
#### delete NAs, zero variance
for (i in 1:length(SMD)) {
  SMD[[i]] <- SMD[[i]][!is.na(SMD[[i]]$es) & !is.na(SMD[[i]]$var) & SMD[[i]]$var != 0 & !is.na(SMD[[i]]$year), ]
}

### Zr
#### delete NAs, zero variance
for (i in 1:length(Zr)) {
  Zr[[i]] <- Zr[[i]][!is.na(Zr[[i]]$es) & !is.na(Zr[[i]]$var) & Zr[[i]]$var != 0 & !is.na(Zr[[i]]$year), ]
}

# create the variable of latest-year-centring publication year ---- 
# this was used as a predictor to test time-lag bias (decline effect)
# The reason why creating this variable is to set the intercept conditional on the latest year rather than zero year

## lnRR
for (i in 1:length(lnRR)) {
  lnRR[[i]]$year_pub.l <- as.vector(lnRR[[i]]$year - max(lnRR[[i]]$year))
}
## SMD
for (i in 1:length(SMD)) {
  SMD[[i]]$year_pub.l <- as.vector(SMD[[i]]$year - max(SMD[[i]]$year))
}

## Zr
for (i in 1:length(Zr)) {
  Zr[[i]]$year_pub.l <- as.vector(Zr[[i]]$year - max(Zr[[i]]$year))
}


# create the variable of sampling error ----
# this was used as a predictor to test small-study effect

## lnRR
for (i in 1:length(lnRR)) {
  lnRR[[i]]$sei <- sqrt(lnRR[[i]]$var)
}

## SMD
for (i in 1:length(SMD)) {
  SMD[[i]]$sei <- sqrt(SMD[[i]]$var)
}

## Zr
for (i in 1:length(Zr)) {
  Zr[[i]]$sei <- sqrt(Zr[[i]]$var)
}

# transform effect size, sei, and year_pub.l prior to model fitting ----
# This is to eliminate scale-dependency and to allow for aggregations of model coefficients over different effect size metrics in subsequent second-order meta-analysis

## scale data using respective standard deviation ----
# var_zscore and sei_zscore are Equation 7 in Yang et al. 2023 for explanations

### lnRR
for (i in 1:length(lnRR)) {
  lnRR[[i]]$es_zscore <- scale(lnRR[[i]]$es, center = F, scale = TRUE) # without centring, which is used to estimate intercept
  lnRR[[i]]$var_zscore <- scale(lnRR[[i]]$var, scale = TRUE) - ( (0-mean(lnRR[[i]]$var))/sd(lnRR[[i]]$var) ) 
  lnRR[[i]]$sei_zscore <- scale(lnRR[[i]]$sei, scale = TRUE) - ( (0-mean(lnRR[[i]]$sei))/sd(lnRR[[i]]$sei) ) 
  lnRR[[i]]$year_pub.l_zscore <- scale(lnRR[[i]]$year_pub.l, scale = TRUE)
}

### SMD
for (i in 1:length(SMD)) {
  SMD[[i]]$es_zscore <- scale(SMD[[i]]$es, center = F, scale = TRUE) # without centring, which is used to estimate intercept
  SMD[[i]]$var_zscore <- scale(SMD[[i]]$var, scale = TRUE) - ( (0-mean(SMD[[i]]$var))/sd(SMD[[i]]$var) ) 
  SMD[[i]]$sei_zscore <- scale(SMD[[i]]$sei, scale = TRUE) - ( (0-mean(SMD[[i]]$sei))/sd(SMD[[i]]$sei) ) 
  SMD[[i]]$year_pub.l_zscore <- scale(SMD[[i]]$year_pub.l, scale = TRUE)
}

## also need to scale effective sample size related variables
for (i in 1:length(SMD[SMD_des_filenames])) {
  SMD[SMD_des_filenames][[i]]$ess.var_zscore <- scale(SMD[SMD_des_filenames][[i]]$ess.var, scale = TRUE) - ( (0-mean(SMD[SMD_des_filenames][[i]]$ess.var))/sd(SMD[SMD_des_filenames][[i]]$ess.var) )
  SMD[SMD_des_filenames][[i]]$ess.sei_zscore <- scale(SMD[SMD_des_filenames][[i]]$ess.sei, scale = TRUE) - ( (0-mean(SMD[SMD_des_filenames][[i]]$ess.sei))/sd(SMD[SMD_des_filenames][[i]]$ess.sei) )
}

SMD[SMD_des_filenames][[5]]

### Zr
for (i in 1:length(Zr)) {
  Zr[[i]]$es_zscore <- scale(Zr[[i]]$es, center = F, scale = TRUE) # without centering, which is used to estimate intercept
  Zr[[i]]$var_zscore <- scale(Zr[[i]]$var, scale = TRUE) - ( (0-mean(Zr[[i]]$var))/sd(Zr[[i]]$var) ) # see Equation 7 for explanations
  Zr[[i]]$sei_zscore <- scale(Zr[[i]]$sei, scale = TRUE) - ( (0-mean(Zr[[i]]$sei))/sd(Zr[[i]]$sei) ) # see Equation 7 for explanations
  Zr[[i]]$year_pub.l_zscore <- scale(Zr[[i]]$year_pub.l, scale = TRUE)
}


###
# 2. Multilevel meta-analytic modelling  ----
###

# we will use multilevel meta-analytic model to fit two types data:
# - original scale data
# (and scaled data?)
# 
# For each type of data, we:
# (i) estimate the meta-analytic overall mean, model intercept (beta0)
# (ii) detect potential publication bias - test for small-study (sign & significance of beta1) and decline effects (sign & significance of beta1)
# (iii) correct for publication bias and estimate bias-corrected overall mean (beta0_c)
# 
## original scale data ----
# 
#### (1) estimates of beta0 ----
# fit intercept-only multilevel model to each dataset

####
# meta-analytic overall mean #
####
# lnRR 
model_lnRR <- NA
for (i in 1:length(lnRR)) {
  model_lnRR[i] <- rma.mv(yi = es, V = var, 
                          random = list(~1|study_ID/obs_ID), 
                          method = "REML", test = "t", 
                          data = lnRR[[i]], sparse = TRUE, 
                          control = list(optimizer = "optim")) |> list()
  }

# SMD
model_SMD <- NA
for (i in 1:length(SMD)) {
  model_SMD[i] <- rma.mv(yi = es, V = var, 
                         random = list(~1|study_ID/obs_ID), 
                         method = "REML", test = "t", 
                         data = SMD[[i]], sparse = TRUE, 
                         control = list(optimizer = "optim")) |> list()
  }

# Zr
## Model fitting failed using `optimizer = "optim"` - I use `nlminb`
model_Zr <- NA
for (i in 1:length(Zr)) {
  model_Zr[i] <- rma.mv(yi = es, V = var, 
                        random = list(~1|study_ID/obs_ID), 
                        method = "REML", test = "t", 
                        data = Zr[[i]], sparse=TRUE, 
                        control=list(optimizer="nlminb")) |> list()
  }


#### (2) detect publication bias ----
# we aim for detecting two forms of publication bias: small-study effect  and decline effect.
# we use a full model with sampling error (sei) and publication year (year_pub.l) as moderators to detect publication bias.
# of relevance, sei's slope (beta1) and year_pub.l's slope (beta2) can be used to indicate the occurrence of small-study effect and decline effect, respectively.
# note that for SMD we need to model 'effective sample size' based sampling error analogue (ees.sei), where possible

##### (a) fit full models with sampling error and year as predictors ----
# the point estimate of SMD and lnRR are inherently correlated with their sampling variances. 
# To avoid such ‘artefactual’ correlation between effect size and sampling error, 
# we need to use "effective sample size" based sampling error to let its estimate get rid of point estimate

####
# Full model with error and latest year as predictors #
####
## use sampling error (sei) and latest year (year_pub.l) as predictors for those with calculated effect sizes and sampling variance

# lnRR
model_lnRR_sei.year <- NA
for (i in 1:length(lnRR[lnRR_filenames])) {
  model_lnRR_sei.year[i] <- rma.mv(yi = es, V = var, 
                                   random = list(~1|study_ID/obs_ID), 
                                   mods = ~ sei + year_pub.l, 
                                   method = "REML", test = "t", 
                                   data = lnRR[lnRR_filenames][[i]], sparse = TRUE, 
                                   control=list(optimizer = "optim")) |> list()
  } 

# SMD
model_SMD_sei.year <- NA
for (i in 1:length(SMD[SMD_filenames])) {
  model_SMD_sei.year[i] <- rma.mv(yi = es, V = var, 
                                  random = list(~1|study_ID/obs_ID), 
                                  mods = ~ sei + year_pub.l, 
                                  method = "REML", test = "t", 
                                  data = SMD[SMD_filenames][[i]], sparse = TRUE, 
                                  control=list(optimizer = "nlminb")) |> list()
  }

model_SMD_ess.sei.year <- NA
for (i in 1:length(SMD[SMD_des_filenames])) {
  model_SMD_ess.sei.year[i] <- rma.mv(yi = es, V = var, 
                                      random = list(~1|study_ID/obs_ID), 
                                      mods = ~ ess.sei + year_pub.l, 
                                      method = "REML", test = "t",
                                      data = SMD[SMD_des_filenames][[i]], sparse = TRUE, 
                                      control=list(optimizer = "optim")) |> list()
  } 

# Zr
## Zr does not have the concern of ‘artefactual’ correlation between effect size and sampling error (because the formula to Zr's estimate sampling error has no component of point estimate: 1/(n-3)). So we only need to fit the regression model with sampling error (sei) as a predictor

model_Zr_sei.year <- NA
for (i in 1:length(Zr)) {
  model_Zr_sei.year[i] <- rma.mv(yi = es, V = var, random = list(~1|study_ID/obs_ID), mods = ~ sei + year_pub.l, method = "REML", test = "t", data = Zr[[i]], sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

##### (b) identify the presence of publication bias ----
#'[please refer to Yang et al. (2023) for details]
# we next aim for identify the presence of the small-study effect and decline effect for each meta-analysis. 

# our rational is: for an effect that is expected to be positive, a small study effect and decline effect would be expressed in a positive value of beta1 and negative value of beta2, respectively. In this respect, a slope (beta1 or beta2)) with opposing direction (unexpected sign) indicates no detectable publication bias and subsequently does not require correction for such a bias.

# we use the product of beta0 and beta1 (i.e., beta×beta1) as the signal, that is, if beta0×beta1 is positive, it indicates the examined meta-analysis has a small-study effect (beta1 is in a correct direction).

# similarly, when the product of beta0×beta2 is negative, the examined meta-analysis has a decline effect (beta 2 is in a correct direction)

# check the significance and direction of model regressions from each meta-analysis to identify whether it presents a small-study effect or decline effect

# first to create a data frame containing full-model's parameter estimates
# combine two types full-models (with sei and ess.sei as a predictor, respectively)

###### lnRR ----
model_lnRR_pb <- model_lnRR_sei.year

## extract model model coefficients and their significance test results
model_est_lnRR <- data.frame(case = names(lnRR),
                             es_type = rep("lnRR", length(lnRR)),
                             beta0 = sapply(model_lnRR, function(x) x$beta), # meta-analytic overall mean
                             se_beta0 = sapply(model_lnRR, function(x) x$se), # standard error of beta0
                             pval_beta0 = sapply(model_lnRR, function(x) x$pval), # p value of beta0
                             beta0_c = sapply(model_lnRR_pb, function(x) x$beta[1]), # bias corrected overall mean
                             se_beta0_c = sapply(model_lnRR_pb, function(x) x$se[1]), # standard error of beta0_c
                             pval_beta0_c = sapply(model_lnRR_pb, function(x) x$pval[1]), # p valuer of beta0_c
                             beta1 = sapply(model_lnRR_pb, function(x) x$beta[2]), # sampling error 
                             se_beta1 = sapply(model_lnRR_pb, function(x) x$se[2]), # standard error of beta1
                             pval_beta1 = sapply(model_lnRR_pb, function(x) x$pval[2]), # p value of beta1
                             beta2 = sapply(model_lnRR_pb, function(x) x$beta[3]), # slope of year 
                             se_beta2 = sapply(model_lnRR_pb, function(x) x$se[3]), # standard error of beta2
                             pval_beta2 = sapply(model_lnRR_pb, function(x) x$pval[3]) # p value of beta2
                             )

# then, identify the presence of the small-study effect and decline effect for each meta-analysis. We also figure out the wrong directions of slope, which will be used to inform the parameterization of reduced models

# our rational is: for an effect that is expected to be positive, a small study effect and decline effect would be expressed in a positive value of beta1 and negative value of beta2, respectively. In such a case, a slope (beta1] or beta2)) with opposing direction (unexpected sign) indicates no detectable publication bias and subsequently does not require correction for such a bias

# we use the product of beta0 and beta1 (i.e., beta*beta1) as the signal, that is, if beta0*beta1 is positive, it indicates the examined meta-analysis has a small-study effect (beta1 is in a correct direction)

# of relevance, when the value of beta0*beta2 is negative, the examined meta-analysis has a decline effect (beta 2 is in a correct direction)

# so we first to create two new columns to contain the two products:  beta0*beta1 and beta0*beta2

# ncol(model_est_lnRR) -> 14
model_est_lnRR[15:16] <- data.frame(beta0Tbeta1 = model_est_lnRR$beta0 * model_est_lnRR$beta1, beta0Tbeta2 = model_est_lnRR$beta0 * model_est_lnRR$beta2)

# visual check
model_est_lnRR

# identify the small-study effect - significant beta1 with correct sign
sse_lnRR <- model_est_lnRR |> subset(model_est_lnRR$pval_beta1 < 0.05 & model_est_lnRR$beta0Tbeta1 > 0)
# check which meta-analyses have small-study effects
sse_lnRR$case 


## identify the decline effect - significant beta2 with correct sign 
de_lnRR <- model_est_lnRR |> subset(model_est_lnRR$pval_beta2 < 0.05 & model_est_lnRR$beta0Tbeta2 < 0)
## check which meta-analyses have decline effects
de_lnRR$case 


## identify the concurrence of the small-study effect and decline effect
sse_de_lnRR <- model_est_lnRR |> subset(model_est_lnRR$pval_beta1 < 0.05 & model_est_lnRR$beta0Tbeta1 > 0 & model_est_lnRR$pval_beta2 < 0.05 & model_est_lnRR$beta0Tbeta2 < 0) 
sse_de_lnRR$case


###### SMD ----

model_SMD_pb <- append(model_SMD_sei.year, model_SMD_ess.sei.year)

model_est_SMD <- data.frame(case = names(SMD),
                            es_type = rep("SMD", length(SMD)),
                            beta0 = sapply(model_SMD, function(x) x$beta), # meta-analytic overall mean
                            se_beta0 = sapply(model_SMD, function(x) x$se), # standard error of beta0
                            pval_beta0 = sapply(model_SMD, function(x) x$pval), # p value of beta0
                            beta0_c = sapply(model_SMD_pb, function(x) x$beta[1]), # bias corrected overall mean
                            se_beta0_c = sapply(model_SMD_pb, function(x) x$se[1]), # standard error of beta0_c
                            pval_beta0_c = sapply(model_SMD_pb, function(x) x$pval[1]), # p valuer of beta0_c
                            beta1 = sapply(model_SMD_pb, function(x) x$beta[2]), # slope of sampling error 
                            se_beta1 = sapply(model_SMD_pb, function(x) x$se[2]), # standard error of beta1
                            pval_beta1 = sapply(model_SMD_pb, function(x) x$pval[2]), # p value of beta1
                            beta2 = sapply(model_SMD_pb, function(x) x$beta[3]), #  slope of year 
                            se_beta2 = sapply(model_SMD_pb, function(x) x$se[3]), # standard error of beta2
                            pval_beta2 = sapply(model_SMD_pb, function(x) x$pval[3]) # p value of beta2
                            )

# ncol(model_est_SMD) ->14
model_est_SMD[15:16] <- data.frame(beta0Tbeta1 = model_est_SMD$beta0 * model_est_SMD$beta1, beta0Tbeta2 = model_est_SMD$beta0 * model_est_SMD$beta2) 

# identify the small-study effect - significant beta1 with correct sign
sse_SMD <- model_est_SMD |> subset(model_est_SMD$pval_beta1 < 0.05 & model_est_SMD$beta0Tbeta1 > 0)
# check which meta-analyses have small-study effects
sse_SMD$case 


# identify the decline effect - significant beta2 with correct sign 
de_SMD <- model_est_SMD |> subset(model_est_SMD$pval_beta2 < 0.05 & model_est_SMD$beta0Tbeta2 < 0)
# check which meta-analyses have decline effects
de_SMD$case # character(0)


# identify the concurrence of the small-study effect and decline effect
sse_de_SMD <- model_est_SMD |> subset(model_est_SMD$pval_beta1 < 0.05 & model_est_SMD$beta0Tbeta1 > 0 & model_est_SMD$pval_beta2 < 0.05 & model_est_SMD$beta0Tbeta2 < 0) 
sse_de_SMD$case # character(0)


###### Zr -----
model_Zr_pb <- model_Zr_sei.year
model_est_Zr <- data.frame(case = names(Zr),
                          es_type = rep("Zr", length(Zr)),
                          beta0 = sapply(model_Zr, function(x) x$beta),
                          se_beta0 = sapply(model_Zr, function(x) x$se),
                          pval_beta0 = sapply(model_Zr, function(x) x$pval),
                          beta0_c = sapply(model_Zr_pb, function(x) x$beta[1]),
                          se_beta0_c = sapply(model_Zr_pb, function(x) x$se[1]),
                          pval_beta0_c = sapply(model_Zr_pb, function(x) x$pval[1]),
                          beta1 = sapply(model_Zr_pb, function(x) x$beta[2]), 
                          se_beta1 = sapply(model_Zr_pb, function(x) x$se[2]), 
                          pval_beta1 = sapply(model_Zr_pb, function(x) x$pval[2]), 
                          beta2 = sapply(model_Zr_pb, function(x) x$beta[3]), 
                          se_beta2 = sapply(model_Zr_pb, function(x) x$se[3]), 
                          pval_beta2 = sapply(model_Zr_pb, function(x) x$pval[3])
                          )

# ncol(model_est_Zr)
model_est_Zr[15:16] <- data.frame(beta0Tbeta1 = model_est_Zr$beta0 * model_est_Zr$beta1, beta0Tbeta2 = model_est_Zr$beta0 * model_est_Zr$beta2) 

# identify the small-study effect - significant beta1 with correct sign
sse_Zr <- model_est_Zr |> subset(model_est_Zr$pval_beta1 < 0.05 & model_est_Zr$beta0Tbeta1 > 0)
sse_Zr$case

# identify the decline effect - significant beta2 with correct sign
de_Zr <- model_est_Zr |> subset(model_est_Zr$pval_beta2 < 0.05 & model_est_Zr$beta0Tbeta2 <0)
de_Zr$case

# both effects
sse_de_Zr <- model_est_Zr |> subset(model_est_Zr$pval_beta1 < 0.05 & model_est_Zr$beta0Tbeta1 > 0&
                                      model_est_Zr$pval_beta2 < 0.05 & model_est_Zr$beta0Tbeta2 <0)

#### (3) correct for publication bias ---- 
# we aim for estimating bias-corrected beta0 according to 4 scenarios (you can check them in Yafeng et al. 2023)
# if a model slope (beta1/beta2) has wrong direction, we need to take out it when fitting model to estimate the bias-corrected mean

##### Multilevel models to estimate bias-corrected effect ----

###### lnRR ----
#'[scenario 1 - both beta1 and beta2 has correct direction]
# in scenario 1, two slopes have correct direction; we can use full model directly, no need to take out any predictor

beta1c_beta2c_lnRR <- model_est_lnRR |> subset(model_est_lnRR$beta0Tbeta1 > 0 &
                                               model_est_lnRR$beta0Tbeta2 < 0)
beta1c_beta2c_lnRR$case

#make a data list only contains scenario1's data
s1_file <- beta1c_beta2c_lnRR$case

## subset of sei
s1_sei_file <- lnRR_filenames[lnRR_filenames %in% s1_file] # this subset should use sei as a predictor and belong to scenario 1

# model fitting - fit full model
model_lnRR_sei_s1 <- NA
for (i in 1: length(s1_sei_file)){
  model_lnRR_sei_s1[i] <- rma.mv(yi = es, V = var,
                                 random = list(~1|study_ID/obs_ID),
                                 mod = sei + year_pub.l, method = "REML",
                                 test = "t", data = lnRR[s1_sei_file][[i]],
                                 sparse = TRUE, control = list(optimizer = "optim")) |>  list()
  }

## replace sei by var
model_lnRR_var_s1 <- NA
for (i in 1: length(s1_sei_file)){
  model_lnRR_var_s1[i] <- rma.mv(yi = es, V = var,
                                 random = list(~1|study_ID/obs_ID),
                                 mod = var + year_pub.l, method = "REML",
                                 test = "t", data = lnRR[s1_sei_file][[i]],
                                 sparse = TRUE, control = list(optimizer = "optim")) |>  list()
  }

# extract model coefficients and their significance test results for 'sei' in scenario 1
model_est_lnRR_sei_s1 <- data.frame(case = s1_sei_file,
                                    es_type = rep("lnRR", length(s1_sei_file)),
                                    beta0 = model_est_lnRR[model_est_lnRR$case %in% s1_sei_file, ]$beta0,
                                    se_beta0 = model_est_lnRR[model_est_lnRR$case %in% s1_sei_file, ]$se_beta0,
                                    pval_beta0 = model_est_lnRR[model_est_lnRR$case %in% s1_sei_file, ]$pval_beta0, 
                                    beta0_c = sapply(model_lnRR_sei_s1, function(x) x$beta[1]), 
                                    se_beta0_c = sapply(model_lnRR_sei_s1, function(x) x$se[1]), 
                                    pval_beta0_c = sapply(model_lnRR_sei_s1, function(x) x$pval[1]),
                                    beta1 = sapply(model_lnRR_sei_s1, function(x) x$beta[2]),
                                    se_beta1 = sapply(model_lnRR_sei_s1, function(x) x$se[2]),
                                    pval_beta1 = sapply(model_lnRR_sei_s1, function(x) x$pval[2]), 
                                    beta2 = sapply(model_lnRR_sei_s1, function(x) x$beta[3]), 
                                    se_beta2 = sapply(model_lnRR_sei_s1, function(x) x$se[3]), 
                                    pval_beta2 = sapply(model_lnRR_sei_s1, function(x) x$pval[3]), 
                                    beta0_c2 = sapply(model_lnRR_var_s1, function(x) x$beta[1]), 
                                    se_beta0_c2 = sapply(model_lnRR_var_s1, function(x) x$se[1]),
                                    pval_beta0_c2 = sapply(model_lnRR_var_s1, function(x) x$pval[1]) 
                                    )


#'[scenario 2 - beta1 has wrong direction / beta2 has correct direction]
beta1w_beta2c_lnRR <- model_est_lnRR |> subset(model_est_lnRR$beta0Tbeta1 < 0 &
                                                 model_est_lnRR$beta0Tbeta2 < 0) 
beta1w_beta2c_lnRR$case 
# character(0)

#'[scenario 3 - beta1 has correct direction / beta2 has wrong direction]
beta1c_beta2w_lnRR <- model_est_lnRR |> subset(model_est_lnRR$beta0Tbeta1 > 0 &
                                                model_est_lnRR$beta0Tbeta2 > 0) 
beta1c_beta2w_lnRR$case 

## data list which only contains scenario3's data
s3_file <- beta1c_beta2w_lnRR$case

## subset of sei
# this subset should fit sei as a predictor and belong to scenario 3
s3_sei_file <- lnRR_filenames[lnRR_filenames %in% s3_file] 

## model fitting - keep beta1-related predictor (sei) and take out beta2-related predictor (year_pub.l)
model_lnRR_sei_s3 <- NA
for (i in 1:length(s3_sei_file)) {
  model_lnRR_sei_s3[i] <- rma.mv(yi = es, V = var, 
                                 random = list(~1|study_ID/obs_ID), 
                                 mods = ~ sei, method = "REML",
                                 test = "t", data = lnRR[s3_sei_file][[i]], 
                                 sparse = TRUE, control = list(optimizer = "optim")) |> list()
  }

# replace var by sei
model_lnRR_var_s3 <- NA
for (i in 1:length(s3_sei_file)) {
  model_lnRR_var_s3[i] <- rma.mv(yi = es, V = var, 
                                 random = list(~1|study_ID/obs_ID), 
                                 mods = ~ var, method = "REML", 
                                 test = "t", data = lnRR[s3_sei_file][[i]], 
                                 sparse = TRUE, control = list(optimizer = "optim")) |> list()
  }

model_est_lnRR_sei_s3 <- data.frame(case = s3_sei_file,
                                    es_type = rep("lnRR", length(s3_sei_file)),
                                    beta0 = model_est_lnRR[model_est_lnRR$case %in% s3_sei_file, ]$beta0,
                                    se_beta0 = model_est_lnRR[model_est_lnRR$case %in% s3_sei_file, ]$se_beta0,
                                    pval_beta0 = model_est_lnRR[model_est_lnRR$case %in% s3_sei_file, ]$pval_beta0,
                                    beta0_c = sapply(model_lnRR_sei_s3, function(x) x$beta[1]), 
                                    se_beta0_c = sapply(model_lnRR_sei_s3, function(x) x$se[1]), 
                                    pval_beta0_c = sapply(model_lnRR_sei_s3, function(x) x$pval[1]),
                                    beta1 = sapply(model_lnRR_sei_s3, function(x) x$beta[2]), 
                                    se_beta1 = sapply(model_lnRR_sei_s3, function(x) x$se[2]), 
                                    pval_beta1 = sapply(model_lnRR_sei_s3, function(x) x$pval[2]),
                                    beta2 = 0, # slope of year: no year term 
                                    se_beta2 = 0,
                                    pval_beta2 = 0, 
                                    beta0_c2 = sapply(model_lnRR_var_s3, function(x) x$beta[1]),
                                    se_beta0_c2 = sapply(model_lnRR_var_s3, function(x) x$se[1]), 
                                    pval_beta0_c2 = sapply(model_lnRR_var_s3, function(x) x$pval[1]) 
                                    )


#'[scenario 4 - beta1 and beta2 have wrong directions]
# needs to take out both of the two predictors (sei and pub_year.l). 
# this is equivalent to a null model (intercept-only model), which is used to estimate (uncorrected) meta-analytic overall mean

beta1w_beta2w_lnRR <- model_est_lnRR |>  subset(model_est_lnRR$beta0Tbeta1 < 0 &
                                                model_est_lnRR$beta0Tbeta2 > 0) 
beta1w_beta2w_lnRR$case 

s4_file <-  beta1w_beta2w_lnRR$case
model_lnRR_s4 <- NA
for (i in 1:length(s4_file)) {
  model_lnRR_s4[i] <- rma.mv(yi = es, V = var, 
                             random = list(~1|study_ID/obs_ID), 
                             method = "REML", test = "t", 
                             data = lnRR[s4_file][[i]], sparse = TRUE, 
                             control=list(optimizer = "optim")) |> list()
 }

model_est_lnRR_s4 <- data.frame(case = s4_file,
                                es_type = rep("lnRR", length(s4_file)),
                                beta0 = model_est_lnRR[model_est_lnRR$case %in% s4_file, ]$beta0, 
                                se_beta0 = model_est_lnRR[model_est_lnRR$case %in% s4_file, ]$se_beta0,
                                pval_beta0 = model_est_lnRR[model_est_lnRR$case %in% s4_file, ]$pval_beta0, 
                                beta0_c = sapply(model_lnRR_s4, function(x) x$beta[1]), 
                                se_beta0_c = sapply(model_lnRR_s4, function(x) x$se[1]), 
                                pval_beta0_c = sapply(model_lnRR_s4, function(x) x$pval[1]), 
                                beta1 = 0, # slope of sampling error: no error term
                                se_beta1 = 0, # standard error of beta1
                                pval_beta1 = 0, # p value of beta1
                                beta2 = 0, # slope of year 
                                se_beta2 = 0, # standard error of beta2
                                pval_beta2 = 0, # p value of beta2
                                beta0_c2 = sapply(model_lnRR_s4, function(x) x$beta[1]), # this is equal to uncorrected overall mean (i.e., beta0)
                                se_beta0_c2 = sapply(model_lnRR_s4, function(x) x$se[1]), # standard error of beta0_c
                                pval_beta0_c2 = sapply(model_lnRR_s4, function(x) x$pval[1])
                                )

###### SMD ----
#'[scenario 1 - both beta1 and beta2 has correct direction]

beta1c_beta2c_SMD <- model_est_SMD |> subset(model_est_SMD$beta0Tbeta1 > 0 & 
                                             model_est_SMD$beta0Tbeta2 < 0) 
beta1c_beta2c_SMD$case  

s1_file <- beta1c_beta2c_SMD$case
## fit ess.sei where possible
## subset of sei
s1_sei_file <- SMD_filenames[SMD_filenames %in% s1_file] 

model_SMD_sei_s1 <- NA
for (i in 1:length(s1_sei_file)) {
  model_SMD_sei_s1[i] <- rma.mv(yi = es, V = var, 
                                random = list(~1|study_ID/obs_ID), 
                                mods = ~ sei + year_pub.l, method = "REML", 
                                test = "t", data = SMD[s1_sei_file][[i]], 
                                sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

## replace sei by var
model_SMD_var_s1 <- NA
for (i in 1:length(s1_sei_file)) {
  model_SMD_var_s1[i] <- rma.mv(yi = es, V = var, 
                                random = list(~1|study_ID/obs_ID),
                                mods = ~ var + year_pub.l, method = "REML", 
                                test = "t", data = SMD[s1_sei_file][[i]], 
                                sparse = TRUE, control=list(optimizer = "optim")) |> list()
 }

## extract model model coefficients and their significance test results for 'sei' in scenario 1
model_est_SMD_sei_s1 <- data.frame(case = s1_sei_file,
                                   es_type = rep("SMD", length(s1_sei_file)),
                                   beta0 = model_est_SMD[model_est_SMD$case %in% s1_sei_file, ]$beta0, 
                                   se_beta0 = model_est_SMD[model_est_SMD$case %in% s1_sei_file, ]$se_beta0,
                                   pval_beta0 = model_est_SMD[model_est_SMD$case %in% s1_sei_file, ]$pval_beta0,
                                   beta0_c = sapply(model_SMD_sei_s1, function(x) x$beta[1]), 
                                   se_beta0_c = sapply(model_SMD_sei_s1, function(x) x$se[1]), 
                                   pval_beta0_c = sapply(model_SMD_sei_s1, function(x) x$pval[1]), 
                                   beta1 = sapply(model_SMD_sei_s1, function(x) x$beta[2]), 
                                   se_beta1 = sapply(model_SMD_sei_s1, function(x) x$se[2]), 
                                   pval_beta1 = sapply(model_SMD_sei_s1, function(x) x$pval[2]), 
                                   beta2 = sapply(model_SMD_sei_s1, function(x) x$beta[3]), 
                                   se_beta2 = sapply(model_SMD_sei_s1, function(x) x$se[3]), 
                                   pval_beta2 = sapply(model_SMD_sei_s1, function(x) x$pval[3]), 
                                   beta0_c2 = sapply(model_SMD_var_s1, function(x) x$beta[1]), 
                                   se_beta0_c2 = sapply(model_SMD_var_s1, function(x) x$se[1]),
                                   pval_beta0_c2 = sapply(model_SMD_var_s1, function(x) x$pval[1]) 
                                   )

# subset of ess.sei
# this subset should fit ess.sei as a predictor and belong to scenario 1
s1_ess.sei_file <- SMD_des_filenames[SMD_des_filenames %in% s1_file] 

## model fitting - full model, which does not need to take out any predictor
model_SMD_ess.sei_s1 <- NA
for (i in 1:length(s1_ess.sei_file)) {
  model_SMD_ess.sei_s1[i] <- rma.mv(yi = es, V = var, 
                                    random = list(~1|study_ID/obs_ID), 
                                    mods = ~ ess.sei + year_pub.l, method = "REML",
                                    test = "t", data = SMD[s1_ess.sei_file][[i]], 
                                    sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

## replace sei by var
model_SMD_ess.var_s1 <- NA
for (i in 1:length(s1_ess.sei_file)) {
  model_SMD_ess.var_s1[i] <- rma.mv(yi = es, V = var, 
                                    random = list(~1|study_ID/obs_ID), 
                                    mods = ~ ess.var + year_pub.l, method = "REML", 
                                    test = "t", data = SMD[s1_ess.sei_file][[i]], 
                                    sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

## extract model model coefficients and their significance test results 'see.sei' for scenario 2 
model_est_SMD_ess.sei_s1 <- data.frame(case = s1_ess.sei_file,
                                       es_type = rep("SMD", length(s1_ess.sei_file)),
                                       beta0 = model_est_SMD[model_est_SMD$case %in% s1_ess.sei_file, ]$beta0,
                                       se_beta0 = model_est_SMD[model_est_SMD$case %in% s1_ess.sei_file, ]$se_beta0,
                                       pval_beta0 = model_est_SMD[model_est_SMD$case %in% s1_ess.sei_file, ]$pval_beta0,
                                       beta0_c = sapply(model_SMD_ess.sei_s1, function(x) x$beta[1]), 
                                       se_beta0_c = sapply(model_SMD_ess.sei_s1, function(x) x$se[1]), 
                                       pval_beta0_c = sapply(model_SMD_ess.sei_s1, function(x) x$pval[1]), 
                                       beta1 = sapply(model_SMD_ess.sei_s1, function(x) x$beta[2]),
                                       se_beta1 = sapply(model_SMD_ess.sei_s1, function(x) x$se[2]), 
                                       pval_beta1 = sapply(model_SMD_ess.sei_s1, function(x) x$pval[2]),
                                       beta2 = sapply(model_SMD_ess.sei_s1, function(x) x$beta[3]), 
                                       se_beta2 = sapply(model_SMD_ess.sei_s1, function(x) x$se[3]), 
                                       pval_beta2 = sapply(model_SMD_ess.sei_s1, function(x) x$pval[3]),
                                       beta0_c2 = sapply(model_SMD_ess.var_s1, function(x) x$beta[1]), 
                                       se_beta0_c2 = sapply(model_SMD_ess.var_s1, function(x) x$se[1]), 
                                       pval_beta0_c2 = sapply(model_SMD_ess.var_s1, function(x) x$pval[1])
                                       )
#'[scenario 2 - beta1 has wrong direction / beta2 has correct direction]

beta1w_beta2c_SMD <- model_est_SMD |> subset(model_est_SMD$beta0Tbeta1 < 0 & 
                                               model_est_SMD$beta0Tbeta2 < 0)
beta1w_beta2c_SMD$case 

# make a data list which only contains scenario2's data
s2_file <-  beta1w_beta2c_SMD$case

## model fitting -  take out beta1-related predictor (sei or esss.sei) and keep beta2-related predictor (year_pub.l)
## beta1-related predictor is removed, so no need to fit ess.sei where possible
model_SMD_s2 <- NA
for (i in 1:length(s2_file)) {
  model_SMD_s2[i] <- rma.mv(yi = es, V = var, 
                            random = list(~1|study_ID/obs_ID), 
                            mods = ~ year_pub.l, method = "REML", 
                            test = "t", data = SMD[s2_file][[i]], 
                            sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

model_est_SMD_s2 <- data.frame(case = s2_file,
                              es_type = rep("SMD", length(s2_file)),
                              beta0 = model_est_SMD[model_est_SMD$case %in% s2_file, ]$beta0, 
                              se_beta0 = model_est_SMD[model_est_SMD$case %in% s2_file, ]$se_beta0, 
                              pval_beta0 = model_est_SMD[model_est_SMD$case %in% s2_file, ]$pval_beta0, 
                              beta0_c = sapply(model_SMD_s2, function(x) x$beta[1]), 
                              se_beta0_c = sapply(model_SMD_s2, function(x) x$se[1]),
                              pval_beta0_c = sapply(model_SMD_s2, function(x) x$pval[1]), 
                              beta1 = 0,
                              se_beta1 = 0, 
                              pval_beta1 = 0,
                              beta2 = sapply(model_SMD_s2, function(x) x$beta[2]), 
                              se_beta2 = sapply(model_SMD_s2, function(x) x$se[2]), 
                              pval_beta2 = sapply(model_SMD_s2, function(x) x$pval[2]),
                              beta0_c2 = sapply(model_SMD_s2, function(x) x$beta[1]), 
                              se_beta0_c2 = sapply(model_SMD_s2, function(x) x$se[1]),
                              pval_beta0_c2 = sapply(model_SMD_s2, function(x) x$pval[1])
                              )

#'[scenario 3 - beta1 has correct direction / beta2 has wrong direction]

beta1c_beta2w_SMD <- model_est_SMD |> subset(model_est_SMD$beta0Tbeta1 > 0 & 
                                             model_est_SMD$beta0Tbeta2 > 0)

beta1c_beta2w_SMD$case

# make a data list which only contains scenario3's data
s3_file <- beta1c_beta2w_SMD$case
## subset of sei
s3_sei_file <- SMD_filenames[SMD_filenames %in% s3_file]

## model fitting - keep beta1-related predictor (sei) and take out beta2-related predictor (year_pub.l)
model_SMD_sei_s3 <- NA
for (i in 1:length(s3_sei_file)) {
  model_SMD_sei_s3[i] <- rma.mv(yi = es, V = var, 
                                random = list(~1|study_ID/obs_ID), 
                                mods = ~ sei, method = "REML", 
                                test = "t", data = SMD[s3_sei_file][[i]], 
                                sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

## replace sei by var
model_SMD_var_s3 <- NA
for (i in 1:length(s3_sei_file)) {
  model_SMD_var_s3[i] <- rma.mv(yi = es, V = var, 
                                random = list(~1|study_ID/obs_ID), 
                                mods = ~ var, method = "REML", 
                                test = "t", data = SMD[s3_sei_file][[i]], 
                                sparse = TRUE, control=list(optimizer = "optim")) |> list()
}

# subset of ess.sei
## this subset should fit ess.sei as a predictor and belong to scenario 3
s3_ess.sei_file <- SMD_des_filenames[SMD_des_filenames %in% s3_file] 

## model fitting - keep beta1-related predictor (sei) and take out beta2-related predictor (year_pub.l)
model_SMD_ess.sei_s3 <- NA
for (i in 1:length(s3_ess.sei_file)) {
  model_SMD_ess.sei_s3[i] <- rma.mv(yi = es, V = var, 
                                    random = list(~1|study_ID/obs_ID), 
                                    mods = ~ ess.sei, method = "REML", 
                                    test = "t", data = SMD[s3_ess.sei_file][[i]], 
                                    sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

# replace sei by var
model_SMD_ess.var_s3 <- NA
for (i in 1:length(s3_ess.sei_file)) {
  model_SMD_ess.var_s3[i] <- rma.mv(yi = es, V = var, 
                                    random = list(~1|study_ID/obs_ID), 
                                    mods = ~ ess.var, method = "REML", 
                                    test = "t", data = SMD[s3_ess.sei_file][[i]],
                                    sparse = TRUE, control=list(optimizer = "optim")) |> list()
  }

## extract model model coefficients and their significance test results for 'sei' in scenario 3
model_est_SMD_sei_s3 <- data.frame(case = s3_sei_file,
                                   es_type = rep("SMD", length(s3_sei_file)),
                                   beta0 = model_est_SMD[model_est_SMD$case %in% s3_sei_file, ]$beta0, 
                                   se_beta0 = model_est_SMD[model_est_SMD$case %in% s3_sei_file, ]$se_beta0, 
                                   pval_beta0 = model_est_SMD[model_est_SMD$case %in% s3_sei_file, ]$pval_beta0, 
                                   beta0_c = sapply(model_SMD_sei_s3, function(x) x$beta[1]), 
                                   se_beta0_c = sapply(model_SMD_sei_s3, function(x) x$se[1]), 
                                   pval_beta0_c = sapply(model_SMD_sei_s3, function(x) x$pval[1]), 
                                   beta1 = sapply(model_SMD_sei_s3, function(x) x$beta[2]), 
                                   se_beta1 = sapply(model_SMD_sei_s3, function(x) x$se[2]), 
                                   pval_beta1 = sapply(model_SMD_sei_s3, function(x) x$pval[2]), 
                                   beta2 = 0, 
                                   se_beta2 = 0,
                                   pval_beta2 = 0, 
                                   beta0_c2 = sapply(model_SMD_var_s3, function(x) x$beta[1]), 
                                   se_beta0_c2 = sapply(model_SMD_var_s3, function(x) x$se[1]), 
                                   pval_beta0_c2 = sapply(model_SMD_var_s3, function(x) x$pval[1]) 
                                   )

## extract model model coefficients and their significance test results for 'ess.sei' in scenario 3
model_est_SMD_ess.sei_s3 <- data.frame(case = s3_ess.sei_file,
                                       es_type = rep("SMD", length(s3_ess.sei_file)),
                                       beta0 = model_est_SMD[model_est_SMD$case %in% s3_ess.sei_file, ]$beta0, 
                                       se_beta0 = model_est_SMD[model_est_SMD$case %in% s3_ess.sei_file, ]$se_beta0, 
                                       pval_beta0 = model_est_SMD[model_est_SMD$case %in% s3_ess.sei_file, ]$pval_beta0, 
                                       beta0_c = sapply(model_SMD_ess.sei_s3, function(x) x$beta[1]), 
                                       se_beta0_c = sapply(model_SMD_ess.sei_s3, function(x) x$se[1]), 
                                       pval_beta0_c = sapply(model_SMD_ess.sei_s3, function(x) x$pval[1]), 
                                       beta1 = sapply(model_SMD_ess.sei_s3, function(x) x$beta[2]), 
                                       se_beta1 = sapply(model_SMD_ess.sei_s3, function(x) x$se[2]), 
                                       pval_beta1 = sapply(model_SMD_ess.sei_s3, function(x) x$pval[2]), 
                                       beta2 = 0, 
                                       se_beta2 = 0, 
                                       pval_beta2 = 0, 
                                       beta0_c2 = sapply(model_SMD_ess.var_s3, function(x) x$beta[1]), 
                                       se_beta0_c2 = sapply(model_SMD_ess.var_s3, function(x) x$se[1]),
                                       pval_beta0_c2 = sapply(model_SMD_ess.var_s3, function(x) x$pval[1])
                                       )

#'[scenario 4 - beta1 and beta2 have wrong directions]
# this reduced model needs to take out both of the two predictors (sei and pub_year.l). This is equivalent to a null model (intercept-only model), which is used to estimate (uncorrected) meta-analytic overall mean
beta1w_beta2w_SMD <- model_est_SMD |> subset(model_est_SMD$beta0Tbeta1 < 0 & 
                                             model_est_SMD$beta0Tbeta2 > 0)
beta1w_beta2w_SMD$case
s4_file <- beta1w_beta2w_SMD$case

# no need to subset sei and ess.sei because scenario4 fits a null model (without any predictor)
# model fitting
## take out both beta1-related predictor (sei or ess.sei) and  beta2-related predictor (year_pub.l)
model_SMD_s4 <- NA
for (i in 1:length(s4_file)) {
  model_SMD_s4[i] <- rma.mv(yi = es, V = var, 
                            random = list(~1|study_ID/obs_ID), 
                            method = "REML", test = "t", 
                            data = SMD[s4_file][[i]], sparse = TRUE, 
                            control=list(optimizer = "optim")) |>  list()
}

## extract model model coefficients and their significance test results scenario 4
model_est_SMD_s4 <- data.frame(case = s4_file,
                               es_type = rep("SMD", length(s4_file)),
                               beta0 = model_est_SMD[model_est_SMD$case %in% s4_file, ]$beta0, 
                               se_beta0 = model_est_SMD[model_est_SMD$case %in% s4_file, ]$se_beta0, 
                               pval_beta0 = model_est_SMD[model_est_SMD$case %in% s4_file, ]$pval_beta0, 
                               beta0_c = sapply(model_SMD_s4, function(x) x$beta[1]), 
                               se_beta0_c = sapply(model_SMD_s4, function(x) x$se[1]), 
                               pval_beta0_c = sapply(model_SMD_s4, function(x) x$pval[1]), 
                               beta1 = 0, 
                               se_beta1 = 0, 
                               pval_beta1 = 0, 
                               beta2 = 0, 
                               se_beta2 = 0, 
                               pval_beta2 = 0, 
                               beta0_c2 = sapply(model_SMD_s4, function(x) x$beta[1]), 
                               se_beta0_c2 = sapply(model_SMD_s4, function(x) x$se[1]), 
                               pval_beta0_c2 = sapply(model_SMD_s4, function(x) x$pval[1]) 
                               )

###### Zr ----
#'[scenario 1 - both beta1 and beta2 has correct direction]
# if a model slope (beta1 and beta2) has a wrong direction, we need to take out it when fitting model to estimate the bias-corrected mean
# in scenario 1, both of the two slopes have a correct direction, we can use full model directly, no need to use take out any predictor
beta1c_beta2c_Zr <- model_est_Zr |> subset(model_est_Zr$beta0Tbeta1 > 0 & 
                                            model_est_Zr$beta0Tbeta2 < 0) 
beta1c_beta2c_Zr$case 
s1_file <- beta1c_beta2c_Zr$case
s1_file <- Zr_filenames[Zr_filenames %in% s1_file] # for Zr, do not need to subset sei and ess.sei

## model fitting - fit a full model
model_Zr_s1 <- NA
for (i in 1:length(s1_file)) {
  model_Zr_s1[i] <- rma.mv(yi = es, V = var, 
                           random = list(~1|study_ID/obs_ID), 
                           mods = ~ sei + year_pub.l, method = "REML", 
                           test = "t", data = Zr[s1_file][[i]], 
                           sparse = TRUE, control=list(optimizer = "optim")) |>  list()
  }

# replace sei by var
model_Zr_var.s1 <- NA
for (i in 1:length(s1_file)) {
  model_Zr_var.s1[i] <- rma.mv(yi = es, V = var, 
                               random = list(~1|study_ID/obs_ID), 
                               mods = ~ var + year_pub.l, method = "REML", 
                               test = "t", data = Zr[s1_file][[i]], sparse = TRUE, 
                               control=list(optimizer = "optim")) |> list()
  }

## extract model coefficients and their significance test results for 'sei' in scenario 1
model_est_Zr_s1 <- data.frame(case = s1_file,
                              es_type = rep("Zr", length(s1_file)),
                              beta0 = model_est_Zr[model_est_Zr$case %in% s1_file, ]$beta0, 
                              se_beta0 = model_est_Zr[model_est_Zr$case %in% s1_file, ]$se_beta0, 
                              pval_beta0 = model_est_Zr[model_est_Zr$case %in% s1_file, ]$pval_beta0, 
                              beta0_c = sapply(model_Zr_s1, function(x) x$beta[1]), 
                              se_beta0_c = sapply(model_Zr_s1, function(x) x$se[1]), 
                              pval_beta0_c = sapply(model_Zr_s1, function(x) x$pval[1]), 
                              beta1 = sapply(model_Zr_s1, function(x) x$beta[2]), 
                              se_beta1 = sapply(model_Zr_s1, function(x) x$se[2]), 
                              pval_beta1 = sapply(model_Zr_s1, function(x) x$pval[2]), 
                              beta2 = sapply(model_Zr_s1, function(x) x$beta[3]), 
                              se_beta2 = sapply(model_Zr_s1, function(x) x$se[3]), 
                              pval_beta2 = sapply(model_Zr_s1, function(x) x$pval[3]), 
                              beta0_c2 = sapply(model_Zr_var.s1, function(x) x$beta[1]), 
                              se_beta0_c2 = sapply(model_Zr_var.s1, function(x) x$se[1]),
                              pval_beta0_c2 = sapply(model_Zr_var.s1, function(x) x$pval[1])
                              )

#'[scenario 2 - beta1 has wrong direction / beta2 has correct direction]

beta1w_beta2c_Zr <- model_est_Zr |> subset(model_est_Zr$beta0Tbeta1 < 0 & 
                                            model_est_Zr$beta0Tbeta2 < 0)
beta1w_beta2c_Zr$case
s2_file <-  beta1w_beta2c_Zr$case

# model fitting - take out beta1-related predictor (sei or esss.sei) and keep beta2-related predictor (year_pub.l)
# beta1-related predictor is removed, so no need to fit ess.sei where possible

model_Zr_s2 <- NA
for (i in 1:length(s2_file)) {
  model_Zr_s2[i] <- rma.mv(yi = es, V = var, 
                           random = list(~1|study_ID/obs_ID), mods = ~ year_pub.l, 
                           method = "REML", test = "t", data = Zr[s2_file][[i]], 
                           sparse = TRUE, control=list(optimizer = "optim")) |> list()
}

model_est_Zr_s2 <- data.frame(case = s2_file,
                              es_type = rep("Zr", length(s2_file)),
                              beta0 = model_est_Zr[model_est_Zr$case %in% s2_file, ]$beta0, 
                              se_beta0 = model_est_Zr[model_est_Zr$case %in% s2_file, ]$se_beta0, 
                              pval_beta0 = model_est_Zr[model_est_Zr$case %in% s2_file, ]$pval_beta0, 
                              beta0_c = sapply(model_Zr_s2, function(x) x$beta[1]), 
                              se_beta0_c = sapply(model_Zr_s2, function(x) x$se[1]), 
                              pval_beta0_c = sapply(model_Zr_s2, function(x) x$pval[1]),
                              beta1 = 0, 
                              se_beta1 = 0, 
                              pval_beta1 = 0, 
                              beta2 = sapply(model_Zr_s2, function(x) x$beta[2]), 
                              se_beta2 = sapply(model_Zr_s2, function(x) x$se[2]),
                              pval_beta2 = sapply(model_Zr_s2, function(x) x$pval[2]), 
                              beta0_c2 = sapply(model_Zr_s2, function(x) x$beta[1]), 
                              se_beta0_c2 = sapply(model_Zr_s2, function(x) x$se[1]), 
                              pval_beta0_c2 = sapply(model_Zr_s2, function(x) x$pval[1])
                              )

#'[scenario 3 - beta1 has correct direction / beta2 has wrong direction]
beta1c_beta2w_Zr <- model_est_Zr |> subset(model_est_Zr$beta0Tbeta1 > 0 & 
                                           model_est_Zr$beta0Tbeta2 > 0)

beta1c_beta2w_Zr$case
s3_file <- beta1c_beta2w_Zr$case

## Zr does not need to subset sei and ess.sei
s3_file <- Zr_filenames[Zr_filenames %in% s3_file] 

# model fitting - keep beta1-related predictor (sei) and take out beta2-related predictor (year_pub.l)
model_Zr_s3 <- NA
for (i in 1:length(s3_file)) {
  model_Zr_s3[i] <- rma.mv(yi = es, V = var, 
                           random = list(~1|study_ID/obs_ID), 
                           mods = ~ sei, method = "REML", 
                           test = "t", data = Zr[s3_file][[i]], 
                           sparse = TRUE, control=list(optimizer = "optim")) |> list()
}

# replace sei by var
model_Zr_var.s3 <- NA
for (i in 1:length(s3_file)) {
  model_Zr_var.s3[i] <- rma.mv(yi = es, V = var, 
                               random = list(~1|study_ID/obs_ID), 
                               mods = ~ var, method = "REML", test = "t", 
                               data = Zr[s3_file][[i]], 
                               sparse = TRUE, control=list(optimizer = "optim")) |> list()
}

## extract model coefficients and their significance test results for 'sei' in scenario 3
model_est_Zr_s3 <- data.frame(case = s3_file,
                              es_type = rep("Zr", length(s3_file)),
                              beta0 = model_est_Zr[model_est_Zr$case %in% s3_file, ]$beta0, 
                              se_beta0 = model_est_Zr[model_est_Zr$case %in% s3_file, ]$se_beta0, 
                              pval_beta0 = model_est_Zr[model_est_Zr$case %in% s3_file, ]$pval_beta0, 
                              beta0_c = sapply(model_Zr_s3, function(x) x$beta[1]), 
                              se_beta0_c = sapply(model_Zr_s3, function(x) x$se[1]), 
                              pval_beta0_c = sapply(model_Zr_s3, function(x) x$pval[1]), 
                              beta1 = sapply(model_Zr_s3, function(x) x$beta[2]), 
                              se_beta1 = sapply(model_Zr_s3, function(x) x$se[2]),
                              pval_beta1 = sapply(model_Zr_s3, function(x) x$pval[2]), 
                              beta2 = 0, 
                              se_beta2 = 0, 
                              pval_beta2 = 0, 
                              beta0_c2 = sapply(model_Zr_var.s3, function(x) x$beta[1]), 
                              se_beta0_c2 = sapply(model_Zr_var.s3, function(x) x$se[1]),
                              pval_beta0_c2 = sapply(model_Zr_var.s3, function(x) x$pval[1])
                              )

#' [scenario 4 - beta1 and beta2 have wrong directions]
beta1w_beta2w_Zr <- model_est_Zr |> subset(model_est_Zr$beta0Tbeta1 < 0 & 
                                             model_est_Zr$beta0Tbeta2 > 0) 
beta1w_beta2w_Zr$case 
# character(0)

###### Contrust bias-corrected dataframe ----

# compare the absolute values of beta0 and beta0_c2 to decide whether the publication bias is correctly adjusted. The principle is that when the bias-corrected meta-analytic overall mean (beta0_c2) is smaller than the original meta-analytic overall mean (beta0), the publication bias is correctly adjusted. Otherwise, the publication bias is not correctly adjusted because the slope of sampling variance (beta1) and publication year (beta2) are not in a expected directionality (which is caused by the un-accounted heterogeneity).

# combine reduced model estimates from lnRR, SMD, and Zr
model_est_all_corrected_original <- rbind(model_est_lnRR_sei_s1,
                                          model_est_lnRR_sei_s3,
                                          model_est_lnRR_s4,
                                          model_est_SMD_sei_s1,
                                          model_est_SMD_ess.sei_s1,
                                          model_est_SMD_s2,
                                          model_est_SMD_sei_s3,
                                          model_est_SMD_ess.sei_s3,
                                          model_est_SMD_s4,
                                          model_est_Zr_s1,
                                          model_est_Zr_s2,
                                          model_est_Zr_s3
                                          )

# create a variable to decide whether beta0 is larger than beta0_c2
model_est_all_corrected_original$beta0Mbeta0_c2 <- abs(model_est_all_corrected_original$beta0) - abs(model_est_all_corrected_original$beta0_c2)

# create a variable to contain the 'real' bias-corrected meta-analytic mean
model_est_all_corrected_original$beta0_c3 <- ifelse(model_est_all_corrected_original$beta0Mbeta0_c2 > 0,
                                                    model_est_all_corrected_original$beta0_c2,
                                                    model_est_all_corrected_original$beta0)

##### Meta-meta-analysis of (original) beta0 ----

# we use random-effect meta-analytic models to aggregate the (original) beta0 to see the overall magnitude of these meta-analyses.

# combine model estimates from lnRR, SMD, and Zr
model_est_all_original <- rbind(model_est_lnRR, model_est_SMD, model_est_Zr) 

# add unique identifier to each meta-analysis paper
model_est_all_original$MA <- model_est_all_original$case |>
  str_extract("^[A-Za-z]+[0-9]+") 

# get folded mean and variance
model_est_all_original$beta0_folded <- folded_es(mean = model_est_all_original$beta0, 
                                                 variance = model_est_all_original$se_beta0^2)
model_est_all_original$beta0_var_folded <- folded_error(mean = model_est_all_original$beta0, 
                                                        variance = model_est_all_original$se_beta0^2)

# overall magnitude of beta0
MMA_beta0_all_original <- rma.mv(yi = beta0_folded, V = beta0_var_folded, 
                                 random = list(~1 | MA, ~1 | case), 
                                 method = "REML", 
                                 data = model_est_all_original) 

# overall decline in different effect size measures
MMA_beta0_all_original_es_type <- rma.mv(yi = beta0_folded, V = beta0_var_folded, 
                                         random = list(~1 | MA, ~1 | case), 
                                         mods = ~ es_type - 1, method = "REML", 
                                         data = model_est_all_original
                                         )
# orchard plot
## for overall mean
p1_MMA_beta0_all_original <- orchard_plot(
  MMA_beta0_all_original, mod = "1", group = "case",
  k = TRUE, g = TRUE, 
  xlab = "Overall mean",
  angle = 90
  )
p1_MMA_beta0_all_original 

## for different effect size measures
p2_MMA_beta0_all_original_es_type <- orchard_plot(
  MMA_beta0_all_original_es_type, mod = "es_type", group = "case",
  k = TRUE, g = TRUE, 
  xlab = "Overall mean",
  angle = 90
  )
p2_MMA_beta0_all_original_es_type

#### (4) Meta-science ----
##### lnRR ----
###### meta-analysis level ----

#' *two-tailed power for meta-analyses*
## power to detect meta-analytic overall mean
model_est_lnRR$MA.power <- power.ma_Shinichi(mu = model_est_lnRR$beta0, SE = model_est_lnRR$se_beta0)

## power to detect bias-corrected meta-analytic overall mean
### add bias-corrected mean to model_est_lnRR
beta0_c3_lnRR <- (model_est_all_corrected_original |> 
                    filter(es_type == "lnRR")) |> 
                    select(case, beta0_c3)
model_est_lnRR <- left_join(model_est_lnRR, beta0_c3_lnRR, by = "case")                                

## calculate bias-corrected power
### still use the unconditional se rather than the conditional se (se of beta0c_3)
model_est_lnRR$MA.power_c <- power.ma_Shinichi(mu = model_est_lnRR$beta0_c3, SE = model_est_lnRR$se_beta0)

## power to detect small-study effect
model_est_lnRR$sse.power <- power.ma_Shinichi(mu = model_est_lnRR$beta1, SE = model_est_lnRR$se_beta1)

## power to detect decline effect
model_est_lnRR$de.power <- power.ma_Shinichi(mu = model_est_lnRR$beta2, SE = model_est_lnRR$se_beta2)


#' *type M error for meta-analysis*

## meta-analytic overall mean
MA.power.M <- NA
for (i in 1:length(model_est_lnRR$case)) {
  MA.power.M[i] <- error_M(mu = model_est_lnRR$beta0[i],
                           se = model_est_lnRR$se_beta0[i],
                           alpha = 0.05, N = 10000) |>  unlist()
  }
model_est_lnRR$MA.power.M <- MA.power.M

## bias-corrected version
MA.power.M_c <- NA
for (i in 1:length(model_est_lnRR$case)) {
  MA.power.M_c[i] <- error_M(mu = model_est_lnRR$beta0_c3[i],
                             se = model_est_lnRR$se_beta0[i],
                             alpha = 0.05, N=10000) |> unlist()
  }
model_est_lnRR$MA.power.M_c <- MA.power.M_c

#' *type S error for meta-analysis*

## meta-analytic overall mean
MA.power.S <- NA
for (i in 1:length(model_est_lnRR$case)) {
  MA.power.S[i] <- error_S(mu = model_est_lnRR$beta0[i],
                           se = model_est_lnRR$se_beta0[i],
                           alpha = 0.05) |> unlist()
  }
model_est_lnRR$MA.power.S <- MA.power.S

## bias-corrected version
MA.power.S_c <- NA
for (i in 1:length(model_est_lnRR$case)) {
  MA.power.S_c[i] <- error_S(mu = model_est_lnRR$beta0_c3[i],
                             se = model_est_lnRR$se_beta0[i],
                             alpha = 0.05) |> unlist()
  }
model_est_lnRR$MA.power.S_c <- MA.power.S_c

###### primary study level ----
# power for primary studies within meta-analysis  
#'*two-tailed power for primary studies*

# for every primary-study data set stored in the list lnRR, calculate the two-tailed statistical power of each experiment and save the results back into the same data set.

## meta-analytic overall mean
power_lnRR <- NA
for (i in 1:length(lnRR)) {
  power_lnRR[i] <- power.individual_Shinichi(mu = model_est_lnRR$beta0[i], se = lnRR[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_lnRR)) {
  lnRR[[i]]$power_lnRR <- power_lnRR[[i]]
  }

## bias-corrected version
power_c_lnRR <- NA
for (i in 1:length(lnRR)) {
  power_c_lnRR[i] <- power.individual_Shinichi(mu = model_est_lnRR$beta0_c3[i], se = lnRR[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_lnRR)) {
  lnRR[[i]]$power_c_lnRR <- power_c_lnRR[[i]]
}


#' *type M error for primary studies*
# meta-analytic overall mean
power.M_lnRR <- NA
for (i in 1:length(lnRR)) {
  power.M_lnRR[i] <- mapply(error_M, mu = model_est_lnRR$beta0[i], se = lnRR[[i]]$sei) |> list()}

for (i in 1:length(power.M_lnRR)) {
  lnRR[[i]]$power.M_lnRR <- power.M_lnRR[[i]]
}

# bias-corrected version
power.M_c_lnRR <- NA
for (i in 1:length(lnRR)) {
  power.M_c_lnRR[i] <- mapply(error_M, mu = model_est_lnRR$beta0_c3[i], se = lnRR[[i]]$sei) |> list()}

for (i in 1:length(power.M_lnRR)) {
  lnRR[[i]]$power.M_c_lnRR <- power.M_c_lnRR[[i]]
}


#' *type S error for primary studies*
# meta-analytic overall mean
power.S_lnRR <- NA
for (i in 1:length(lnRR)) {
  power.S_lnRR[i] <- mapply(error_S,mu=model_est_lnRR$beta0[i], se=lnRR[[i]]$sei) |> list()}

for (i in 1:length(power.S_lnRR)) {
  lnRR[[i]]$power.S_lnRR <- power.S_lnRR[[i]]
}

# bias-corrected version
power.S_c_lnRR <- NA
for (i in 1:length(lnRR)) {
  power.S_c_lnRR[i] <- mapply(error_S,mu=model_est_lnRR$beta0_c3[i], se=lnRR[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.S_lnRR)) {
  lnRR[[i]]$power.S_c_lnRR <- power.S_c_lnRR[[i]]
}

# summary of primary study power #

## two tailed power
### meta-analytic overall mean
power_summary_lnRR <- data.frame(case = names(lnRR),
                                 Minimum = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[1,],      
                                 `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[2,],
                                 Median = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[3,],
                                 Mean = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[4,], 
                                 `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[5,], 
                                 Maximum = mapply(summary, sapply(lnRR, function(x) x$power_lnRR))[6,]) 

### bias-corrected version
power_c_summary_lnRR <- data.frame(case = names(lnRR),
                                   Minimum = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[1,],      
                                   `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[2,],
                                   Median = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[3,],
                                   Mean = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[4,], 
                                   `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[5,], 
                                   Maximum = mapply(summary, sapply(lnRR, function(x) x$power_c_lnRR))[6,]) 

## type M error
### meta-analytic overall mean
power.M_summary_lnRR <-  data.frame(case = names(lnRR),
                                    Minimum = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[1,],      
                                    `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[2,],
                                    Median = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[3,],
                                    Mean = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[4,], 
                                    `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[5,], 
                                    Maximum = mapply(summary, sapply(lnRR, function(x) x$power.M_lnRR))[6,]) 

### bias-corrected version
power.M_c_summary_lnRR <-  data.frame(case = names(lnRR),
                                      Minimum = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[1,],    
                                      `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[2,],
                                      Median = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[3,],
                                      Mean = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[4,], 
                                      `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[5,], 
                                      Maximum = mapply(summary, sapply(lnRR, function(x) x$power.M_c_lnRR))[6,]) 

## type S error 
### meta-analytic overall mean
power.S_summary_lnRR <-  data.frame(case = names(lnRR),
                                    Minimum = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[1,],      
                                    `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[2,],
                                    Median = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[3,],
                                    Mean = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[4,], 
                                    `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[5,], 
                                    Maximum = mapply(summary, sapply(lnRR, function(x) x$power.S_lnRR))[6,]) 

### bias-corrected version
power.S_c_summary_lnRR <-  data.frame(case = names(lnRR),
                                      Minimum = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[1,],      
                                      `First quarter` = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[2,],
                                      Median = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[3,],
                                      Mean = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[4,], 
                                      `Third quarter` = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[5,], 
                                      Maximum = mapply(summary, sapply(lnRR, function(x) x$power.S_c_lnRR))[6,]) 


###### aggregation ----
# This section is used to obtain overall estimates of the three parameters across different meta-analyses (which provided us with comparable summaries of the three parameters).
# We used weighted regression to statistically aggregate over the three parameters obtained at the within-meta-analysis level whereas we used mixed effects models to aggregate these parameters at the primary study level. Both procedures involved aggregating the parameters across meta-analyses (i.e., between-meta-analysis modelling). 

#' [estimate overall power for meta-analysis level power]

## add study_ID
model_est_lnRR$base <- sub("\\.csv$", "", model_est_lnRR$case)
model_est_lnRR$study_ID <- sub("_.*", "", model_est_lnRR$base)

## add N and k
N_lnRR <- NA
for (i in 1:length(lnRR)) {
  N_lnRR[i] <- lnRR[[i]]$study_ID |> unique() |> length()
  }
model_est_lnRR$N <- N_lnRR

k_lnRR <- NA
for (i in 1:length(lnRR)) {
  k_lnRR[i] <- lnRR[[i]]$obs_ID |> length()
  }
model_est_lnRR$k <- k_lnRR

##' *two tailed power* ##

### log
###' *lmer* ###
#' TODO -[ASK: which model should I use? When I run model with lmer(), I cannot get the 95%CI using method = 'profile' (need to use 'Wald'), but when I use lm(),I can get more accurate 95%CI (i.e., method = 'profile') - the result did not change the dramatically - for now, I am using lm()]

# log scale 
#MMA_MA.power_lnRR <- lmer(log(MA.power) ~ 1 + (1 | study_ID), weights = k, data = model_est_lnRR)
# 
# ### original scale
# MMA_MA.power_lnRR2 <- lmer(MA.power ~ 1 + (1 | study_ID), weights = k, data = model_est_lnRR)
# 
# exp(fixef(MMA_MA.power_lnRR)) # median
# 
# exp(fixef(MMA_MA.power_lnRR)[1] + 0.5 * var(log(model_est_lnRR$MA.power))) # mean
# #confidence interval of median
# confint(MMA_MA.power_lnRR, method = "Wald") |> exp()

###' *lm* ###
# log scale
MMA_MA.power_lnRR <- lm(log(MA.power) ~ 1, weights = k, data = model_est_lnRR)

# original scale
MMA_MA.power_lnRR2 <- lm(MA.power ~ 1, weights = k, data = model_est_lnRR)
MMA_MA.power_lnRR2$coefficients

# this is median
MMA_MA.power_lnRR$coefficients  |> exp() 
# this is mean
(MMA_MA.power_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power))) |> exp() 
#confidence interval of median
confint(MMA_MA.power_lnRR) |> exp()

# compare residuals to check normality - which one should we use? 
## we should use log-transformed version
par(mfrow = c(1, 2))
residuals(MMA_MA.power_lnRR) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_MA.power_lnRR2) |> hist(main = paste("original power"), xlab = "Residual") 


# bias-corrected version
# log
MMA_MA.power_c_lnRR <- lm(log(MA.power_c) ~ 1, weights = k, data = model_est_lnRR)

# this is median
MMA_MA.power_c_lnRR$coefficients  |> exp() 

# this is mean
(MMA_MA.power_c_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power))) |> exp() 
#confidence interval of median
confint(MMA_MA.power_c_lnRR) |> exp()

##' *type M error * ##
## meta-analytic overall mean
# log
MMA_MA.power.M_lnRR <- lm(log(MA.power.M) ~ 1, weights = k, data = model_est_lnRR)

# this is median
MMA_MA.power.M_lnRR$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power.M))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_lnRR) |> exp()

## bias-corrected version
# log
MMA_MA.power.M_c_lnRR <- lm(log(MA.power.M_c) ~ 1, weights = k, data = model_est_lnRR)

# this is median
MMA_MA.power.M_c_lnRR$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_c_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power.M_c))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_c_lnRR) |> exp()

##' *type S error * ##
## meta-analytic overall mean
# log
MMA_MA.power.S_lnRR <- lm(log(MA.power.S+0.025) ~ 1, weights = k, data = model_est_lnRR) # add an offset of 0.025(25%) to avoid ln(0) = infinity 

# this is median
MMA_MA.power.S_lnRR$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power.S + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_lnRR) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it

## bias-corrected version
# log
MMA_MA.power.S_c_lnRR <- lm(log(MA.power.S_c + 0.025) ~ 1, weights = k, data = model_est_lnRR) # add an offset of 0.025(25%) to avoid ln(0) = infinity 

# this is median
MMA_MA.power.S_c_lnRR$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_c_lnRR$coefficients + 0.5*var(log(model_est_lnRR$MA.power.S_c + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_c_lnRR) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it

#' [estimate overall power for primary study level power]

# first, we extract and compiles statistical metrics (power, Type S error, Type M error) for individual primary studies included in multiple meta-analyses.
# it works with a list of data frames (lnRR), where each element corresponds to one meta-analysis.
# for each primary study within these meta-analyses, it extracts:
#    raw power (probability of detecting the effect),
#    type S error (probability of getting the sign wrong),
#    type M error (magnitude exaggeration),
#                               and their bias-corrected versions.

study_ID_lnRR <- sapply(lnRR, function(x) x$study_ID) |> unlist()
power_lnRR <- sapply(lnRR, function(x) x$power_lnRR) |> unlist()
power.S_lnRR <- sapply(lnRR, function(x) x$power.S_lnRR) |> unlist()
power.M_lnRR <- sapply(lnRR, function(x) x$power.M_lnRR) |> unlist()
power_c_lnRR <- sapply(lnRR, function(x) x$power_c_lnRR) |> unlist()
power.S_c_lnRR <- sapply(lnRR, function(x) x$power.S_c_lnRR) |> unlist()
power.M_c_lnRR <- sapply(lnRR, function(x) x$power.M_c_lnRR) |> unlist()

# then combines these extracted values into a single data frame called individual_est_lnRR, where each row represents one primary study.
individual_est_lnRR <- data.frame("study_ID_lnRR" = study_ID_lnRR,
                                  "power_lnRR" = power_lnRR,
                                  "power.S_lnRR" = power.S_lnRR,
                                  "power.M_lnRR" = power.M_lnRR,
                                  "power_c_lnRR" = power_c_lnRR,
                                  "power.S_c_lnRR" = power.S_c_lnRR,
                                  "power.M_c_lnRR" = power.M_c_lnRR
                                  )

##' *two tailed power* ##
# log
MMA_EXP.power_lnRR <- lmer(log(power_lnRR) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR)
# original scale
MMA_EXP.power_lnRR2 <- lmer(power_lnRR ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR)
summary(MMA_EXP.power_lnRR2)$coefficients[1]

# this is median 
summary(MMA_EXP.power_lnRR)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_lnRR)$coefficients[1] + 0.5*var(log(individual_est_lnRR$power_lnRR))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_lnRR) |> exp()

# compare residual to check normarity - log-transformed is better
par(mfrow = c(1, 2))
residuals(MMA_EXP.power_lnRR) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_EXP.power_lnRR2) |> hist(main = paste("orignal power"), xlab = "Residual")

## bias-corrected version
# log
MMA_EXP.power_c_lnRR <- lmer(log(power_c_lnRR) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR)

# this is median 
summary(MMA_EXP.power_c_lnRR)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_c_lnRR)$coefficients[1] + 0.5*var(log(individual_est_lnRR$power_lnRR))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_c_lnRR) |> exp()

#' *type M error* ##

## meta-analytic overall mean
# log
MMA_EXP.power.M_lnRR <- lmer(log(power.M_lnRR) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR)

# this is median 
summary(MMA_EXP.power.M_lnRR)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_lnRR)$coefficients[1] + 0.5*var(log(individual_est_lnRR$power.M_lnRR))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_lnRR) |> exp()

## bias-corrected version
# log
MMA_EXP.power.M_c_lnRR <- lmer(log(power.M_c_lnRR) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR)

# this is median 
summary(MMA_EXP.power.M_c_lnRR)$coefficients[1] |> exp()
# this is mean
(summary(MMA_EXP.power.M_c_lnRR)$coefficients[1] + 0.5*var(log(individual_est_lnRR$power.M_lnRR))) |> exp()
# confidence interval of median
confint(MMA_EXP.power.M_c_lnRR) |> exp()


#' *type S error* ##
## meta-analytic overall mean
# log
MMA_EXP.power.S_lnRR <- lmer( log(power.S_lnRR + 0.025) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR) # add an offset of 0.025 to avoid log(0) = inf

# this is median 
summary(MMA_EXP.power.S_lnRR)$coefficients[1] |> exp() - 0.025 # - 0.025 is important

# this is mean
(summary(MMA_EXP.power.S_lnRR)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_lnRR)$varcor[[1]][[1]] + # sigma^2 for study level
           summary(MMA_EXP.power.S_lnRR)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
) |> exp() - 0.025

#confidence interval of median
# if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it
confint(MMA_EXP.power.S_lnRR) |> exp() - 0.025 

## bias-corrected version
# log
MMA_EXP.power.S_c_lnRR <- lmer( log(power.S_c_lnRR + 0.025) ~ 1 + (1 | study_ID_lnRR), data = individual_est_lnRR) 

# this is median 
summary(MMA_EXP.power.S_c_lnRR)$coefficients[1] |> exp() - 0.025

# this is mean
(summary(MMA_EXP.power.S_c_lnRR)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_lnRR)$varcor[[1]][[1]] + 
           summary(MMA_EXP.power.S_lnRR)$sigma^2) 
) |> exp() - 0.025

#confidence interval of median
confint(MMA_EXP.power.S_c_lnRR) |> exp() - 0.025 

##### SMD ----
###### meta-analysis level ----

#' *two-tailed power for meta-analyses*
# power to detect meta-analytic overall mean
model_est_SMD$MA.power <- power.ma_Shinichi(mu = model_est_SMD$beta0,SE = model_est_SMD$se_beta0)

# power to detect bias-corrected meta-analytic overall mean
## add bias-corrected mean to model_est_SMD
beta0_c3_SMD <- (model_est_all_corrected_original |> filter(es_type == "SMD")) |> select(case, beta0_c3)

model_est_SMD <- left_join(model_est_SMD, beta0_c3_SMD, by="case") # reorder rows according to vector of case

## calculate bias-corrected power
model_est_SMD$MA.power_c <- power.ma_Shinichi(mu = model_est_SMD$beta0_c3,SE = model_est_SMD$se_beta0) # note to still use the unconditional se rather than the conditional se (se of beta0c_3)

# power to detect small-study effect
model_est_SMD$sse.power <- power.ma_Shinichi(mu = model_est_SMD$beta1, SE = model_est_SMD$se_beta1)

# power to detect decline effect
model_est_SMD$de.power <- power.ma_Shinichi(mu = model_est_SMD$beta2, SE = model_est_SMD$se_beta2)


#' *type M error for meta-analysis*
# meta-analytic overall mean
MA.power.M <- NA
for (i in 1:length(model_est_SMD$case)) {
  MA.power.M[i] <- error_M(mu = model_est_SMD$beta0[i], se = model_est_SMD$se_beta0[i], 
                           alpha = 0.05, N = 10000) |> unlist()
}

model_est_SMD$MA.power.M <- MA.power.M

# bias-corrected version
MA.power.M_c <- NA
for (i in 1:length(model_est_SMD$case)) {
  MA.power.M_c[i] <- error_M(mu = model_est_SMD$beta0_c3[i], se = model_est_SMD$se_beta0[i],
                             alpha = 0.05,N = 10000) |> unlist()
}

model_est_SMD$MA.power.M_c <- MA.power.M_c


#' *type S error for meta-analysis*

# meta-analytic overall mean
MA.power.S <- NA
for (i in 1:length(model_est_SMD$case)) {
  MA.power.S[i] <- error_S(mu = model_est_SMD$beta0[i], se = model_est_SMD$se_beta0[i], 
                           alpha=0.05) |> unlist()
  }

model_est_SMD$MA.power.S <- MA.power.S

# bias-corrected version
MA.power.S_c <- NA
for (i in 1:length(model_est_SMD$case)) {
  MA.power.S_c[i] <- error_S(mu = model_est_SMD$beta0_c3[i],se = model_est_SMD$se_beta0[i],
                             alpha=0.05) |> unlist()
}

model_est_SMD$MA.power.S_c <- MA.power.S_c

###### primary study level ----
# power for peimary studies within meta-analysis  
#'*two-tailed power for primary studies*

# meta-analytic overall mean
# SMD
power_SMD <- NA
for (i in 1:length(SMD)) {
  power_SMD[i] <- power.individual_Shinichi(mu = model_est_SMD$beta0[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_SMD)) {
  SMD[[i]]$power_SMD <- power_SMD[[i]]
}

# bias-corrected version
# SMD
power_c_SMD <- NA
for (i in 1:length(SMD)) {
  power_c_SMD[i] <- power.individual_Shinichi(mu = model_est_SMD$beta0_c3[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_SMD)) {
  SMD[[i]]$power_c_SMD <- power_c_SMD[[i]]
}

#' *type M error*
# meta-analytic overall mean
power.M_SMD <- NA
for (i in 1:length(SMD)) {
  power.M_SMD[i] <- mapply(error_M, mu = model_est_SMD$beta0[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.M_SMD)) {
  SMD[[i]]$power.M_SMD <- power.M_SMD[[i]]
}

# bias-corrected version
power.M_c_SMD <- NA
for (i in 1:length(SMD)) {
  power.M_c_SMD[i] <- mapply(error_M, mu = model_est_SMD$beta0_c3[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.M_SMD)) {
  SMD[[i]]$power.M_c_SMD <- power.M_c_SMD[[i]]
}


#' *type S error*
# meta-analytic overall mean
power.S_SMD <- NA
for (i in 1:length(SMD)) {
  power.S_SMD[i] <- mapply(error_S, mu = model_est_SMD$beta0[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.S_SMD)) {
  SMD[[i]]$power.S_SMD <- power.S_SMD[[i]]
}

# bias-corrected version
power.S_c_SMD <- NA
for (i in 1:length(SMD)) {
  power.S_c_SMD[i] <- mapply(error_S, mu = model_est_SMD$beta0_c3[i], se = SMD[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.S_SMD)) {
  SMD[[i]]$power.S_c_SMD <- power.S_c_SMD[[i]]
}


# summary of primary study power #
#'*two-tailed power*
# meta-analytic overall mean
power_summary_SMD <- data.frame(case = names(SMD),
                                Minimum = mapply(summary, sapply(SMD, function(x) x$power_SMD))[1,],      
                                `First quarter` = mapply(summary, sapply(SMD, function(x) x$power_SMD))[2,],
                                Median = mapply(summary, sapply(SMD, function(x) x$power_SMD))[3,],
                                Mean = mapply(summary, sapply(SMD, function(x) x$power_SMD))[4,], 
                                `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power_SMD))[5,], 
                                Maximum = mapply(summary, sapply(SMD, function(x) x$power_SMD))[6,]) 

# bias-corrected version
power_c_summary_SMD <- data.frame(case = names(SMD),
                                  Minimum = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[1,],      
                                  `First quarter` = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[2,],
                                  Median = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[3,],
                                  Mean = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[4,], 
                                  `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[5,], 
                                  Maximum = mapply(summary, sapply(SMD, function(x) x$power_c_SMD))[6,]) 

#' *type M error*
# meta-analytic overall mean
power.M_summary_SMD <-  data.frame(case = names(SMD),
                                   Minimum = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[1,],      
                                   `First quarter` = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[2,],
                                   Median = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[3,],
                                   Mean = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[4,], 
                                   `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[5,], 
                                   Maximum = mapply(summary, sapply(SMD, function(x) x$power.M_SMD))[6,]) 

# bias-corrected version
power.M_c_summary_SMD <-  data.frame(case = names(SMD),
                                     Minimum = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[1,],      
                                     `First quarter` = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[2,],
                                     Median = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[3,],
                                     Mean = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[4,], 
                                     `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[5,], 
                                     Maximum = mapply(summary, sapply(SMD, function(x) x$power.M_c_SMD))[6,]) 


#' *type S error*
# meta-analytic overall mean
power.S_summary_SMD <-  data.frame(case = names(SMD),
                                   Minimum = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[1,],      
                                   `First quarter` = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[2,],
                                   Median = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[3,],
                                   Mean = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[4,], 
                                   `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[5,], 
                                   Maximum = mapply(summary, sapply(SMD, function(x) x$power.S_SMD))[6,]) 

# bias-corrected version
power.S_c_summary_SMD <-  data.frame(case = names(SMD),
                                     Minimum = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[1,],      
                                     `First quarter` = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[2,],
                                     Median = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[3,],
                                     Mean = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[4,], 
                                     `Third quarter` = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[5,], 
                                     Maximum = mapply(summary, sapply(SMD, function(x) x$power.S_c_SMD))[6,]) 


###### aggregation ----
#' [estimate overall power for meta-analysis level power]

# add N and k
N_SMD <- NA
for (i in 1:length(SMD)) {
  N_SMD[i] <- SMD[[i]]$study_ID |> unique() |> length()
}
model_est_SMD$N <- N_SMD

k_SMD <- NA
for (i in 1:length(SMD)) {
  k_SMD[i] <- SMD[[i]]$obs_ID |> length()
}
model_est_SMD$k <- k_SMD


#' *two tailed power*
## meta-analytic overall mean
# log
MMA_MA.power_SMD <- lm(log(MA.power) ~ 1, weights = k, data = model_est_SMD)
# original scale
MMA_MA.power_SMD2 <- lm(MA.power ~ 1, weights = k, data = model_est_SMD)
MMA_MA.power_SMD2$coefficients

# this is median (log)
MMA_MA.power_SMD$coefficients  |> exp() 

# this is mean (log)
(MMA_MA.power_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power))) |> exp() 

#confidence interval of median (log)
confint(MMA_MA.power_SMD) |> exp()

# compare residuals (log)
par(mfrow = c(1, 2))
residuals(MMA_MA.power_SMD) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_MA.power_SMD2) |> hist(main = paste("original power"), xlab = "Residual") 

## bias-corrected version
# log
MMA_MA.power_c_SMD <- lm(log(MA.power_c) ~ 1, weights = k, data = model_est_SMD)

# this is median
MMA_MA.power_c_SMD$coefficients  |> exp() 

# this is mean
(MMA_MA.power_c_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power))) |> exp() 

#confidence interval of median
confint(MMA_MA.power_c_SMD) |> exp()


#' *type M error*
## meta-analytic overall mean
# log
MMA_MA.power.M_SMD <- lm(log(MA.power.M) ~ 1, weights = k, data = model_est_SMD)

# this is median
MMA_MA.power.M_SMD$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power.M))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_SMD) |> exp()


## bias-corrected version
# log
MMA_MA.power.M_c_SMD <- lm(log(MA.power.M_c) ~ 1, weights = k, data = model_est_SMD)
# this is median
MMA_MA.power.M_c_SMD$coefficients |> exp() 
# this is mean
(MMA_MA.power.M_c_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power.M_c))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_c_SMD) |> exp()


#' *type S error*
## meta-analytic overall mean
# log
MMA_MA.power.S_SMD <- lm(log(MA.power.S+0.025) ~ 1, weights = k, data = model_est_SMD) # add an offset of 0.025(25%) to avoid ln(0) = infinity 

# this is median
MMA_MA.power.S_SMD$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power.S + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_SMD) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it


# bias-corrected version
# log
MMA_MA.power.S_c_SMD <- lm(log(MA.power.S_c + 0.025) ~ 1, weights = k, data = model_est_SMD) # add an offset of 0.025(25%) to avoid ln(0) = infinity 

# this is median
MMA_MA.power.S_c_SMD$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_c_SMD$coefficients + 0.5*var(log(model_est_SMD$MA.power.S_c + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_c_SMD) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it


#' [estimate overall power for primary study level power]

study_ID_SMD <- sapply(SMD, function(x) x$study_ID) |> unlist()
power_SMD <- sapply(SMD, function(x) x$power_SMD) |> unlist()
power.S_SMD <- sapply(SMD, function(x) x$power.S_SMD) |> unlist()
power.M_SMD <- sapply(SMD, function(x) x$power.M_SMD) |> unlist()
power_c_SMD <- sapply(SMD, function(x) x$power_c_SMD) |> unlist()
power.S_c_SMD <- sapply(SMD, function(x) x$power.S_c_SMD) |> unlist()
power.M_c_SMD <- sapply(SMD, function(x) x$power.M_c_SMD) |> unlist()

individual_est_SMD <- data.frame("study_ID_SMD" = study_ID_SMD,
                                 "power_SMD" = power_SMD,
                                 "power.S_SMD" = power.S_SMD,
                                 "power.M_SMD" = power.M_SMD,
                                 "power_c_SMD" = power_c_SMD,
                                 "power.S_c_SMD" = power.S_c_SMD,
                                 "power.M_c_SMD" = power.M_c_SMD
                                 )


#' *two tailed power*
## meta-analytic overall mean
# log
MMA_EXP.power_SMD <- lmer(log(power_SMD) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD)

# original scale
MMA_EXP.power_SMD2 <- lmer(power_SMD ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD)
summary(MMA_EXP.power_SMD2)$coefficients[1]

# this is median 
summary(MMA_EXP.power_SMD)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_SMD)$coefficients[1] + 0.5*var(log(individual_est_SMD$power_SMD))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_SMD) |> exp()

# compare residual
par(mfrow = c(1, 2))
residuals(MMA_EXP.power_SMD) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_EXP.power_SMD2) |> hist(main = paste("orignal power"), xlab = "Residual")

## bias-corrected version
# use log
MMA_EXP.power_c_SMD <- lmer(log(power_c_SMD) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD)

# this is median 
summary(MMA_EXP.power_c_SMD)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_c_SMD)$coefficients[1] + 0.5*var(log(individual_est_SMD$power_SMD))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_c_SMD) |> exp()


#' *type M error* 
## meta-analytic overall mean
# log
MMA_EXP.power.M_SMD <- lmer(log(power.M_SMD) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD)

# this is median 
summary(MMA_EXP.power.M_SMD)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_SMD)$coefficients[1] + 0.5*var(log(individual_est_SMD$power.M_SMD))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_SMD) |> exp()

# bias-corrected version
# log
MMA_EXP.power.M_c_SMD <- lmer(log(power.M_c_SMD) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD)

# this is median 
summary(MMA_EXP.power.M_c_SMD)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_c_SMD)$coefficients[1] + 0.5*var(log(individual_est_SMD$power.M_SMD))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_c_SMD) |> exp()


#' *type S error* 
## meta-analytic overall mean
# log
MMA_EXP.power.S_SMD <- lmer( log(power.S_SMD + 0.025) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD) # add an offset of 0.025 to avoid log(0) = inf

# this is median 
summary(MMA_EXP.power.S_SMD)$coefficients[1] |> exp() - 0.025 # - 0.025 is important

# this is mean
(summary(MMA_EXP.power.S_SMD)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_SMD)$varcor[[1]][[1]] + # sigma^2 for study level
           summary(MMA_EXP.power.S_SMD)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
) |> exp() - 0.025

#confidence interval of median
# if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it
confint(MMA_EXP.power.S_SMD) |> exp() - 0.025 

## bias-corrected version
# log
MMA_EXP.power.S_c_SMD <- lmer( log(power.S_c_SMD + 0.025) ~ 1 + (1 | study_ID_SMD), data = individual_est_SMD) 

# this is median 
summary(MMA_EXP.power.S_c_SMD)$coefficients[1] |> exp() - 0.025 

# this is mean
(summary(MMA_EXP.power.S_c_SMD)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_SMD)$varcor[[1]][[1]] +
           summary(MMA_EXP.power.S_SMD)$sigma^2) 
) |> exp() - 0.025

#confidence interval of median
confint(MMA_EXP.power.S_c_SMD) |> exp() - 0.025 

##### Zr ----
###### meta-analysis level ----

#' *two-tailed power for meta-analyses*
# power to detect meta-analytic overall mean
model_est_Zr$MA.power <- power.ma_Shinichi(mu = model_est_Zr$beta0, SE = model_est_Zr$se_beta0)

# power to detect bias-corrected meta-analytic overall mean
## add bias-corrected mean to model_est_Zr
beta0_c3_Zr <- (model_est_all_corrected_original |> 
                filter(es_type == "Zr")) |> select(case, beta0_c3)

model_est_Zr <- left_join(model_est_Zr, beta0_c3_Zr, by = "case") # reorder rows according to vector of case

## calculate bias-corrected power
model_est_Zr$MA.power_c <- power.ma_Shinichi(mu = model_est_Zr$beta0_c3, SE = model_est_Zr$se_beta0) # note to still use the unconditional se rather than the conditional se (se of beta0c_3)

# power to detect small-study effect
model_est_Zr$sse.power <- power.ma_Shinichi(mu = model_est_Zr$beta1,SE = model_est_Zr$se_beta1)

# power to detect decline effect
model_est_Zr$de.power <- power.ma_Shinichi(mu = model_est_Zr$beta2, SE = model_est_Zr$se_beta2)


#' *type M error for meta-analyses*
# meta-analytic overall mean
MA.power.M <- NA
for (i in 1:length(model_est_Zr$case)) {
  MA.power.M[i] <- error_M(mu = model_est_Zr$beta0[i], se = model_est_Zr$se_beta0[i],
                           alpha = 0.05,N = 10000) |> unlist()
  }

model_est_Zr$MA.power.M <- MA.power.M

# bias-corrected version
MA.power.M_c <- NA
for (i in 1:length(model_est_Zr$case)) {
  MA.power.M_c[i] <- error_M(mu = model_est_Zr$beta0_c3[i],se = model_est_Zr$se_beta0[i],
                             alpha = 0.05,N = 10000) |> unlist()
  }

model_est_Zr$MA.power.M_c <- MA.power.M_c


#' *type S error for meta-analyses*
# meta-analytic overall mean
MA.power.S <- NA
for (i in 1:length(model_est_Zr$case)) {
  MA.power.S[i] <- error_S(mu = model_est_Zr$beta0[i], se = model_est_Zr$se_beta0[i],
                           alpha = 0.05) |> unlist()
  }

model_est_Zr$MA.power.S <- MA.power.S

# bias-corrected version
MA.power.S_c <- NA
for (i in 1:length(model_est_Zr$case)) {
  MA.power.S_c[i] <- error_S(mu = model_est_Zr$beta0_c3[i], se = model_est_Zr$se_beta0[i],
                             alpha = 0.05) |> unlist()
  }

model_est_Zr$MA.power.S_c <- MA.power.S_c

###### primary study level ----
# power for primary studies within meta-analysis  

#'*two-tailed power for primary studies*
# meta-analytic overall mean
power_Zr <- NA
for (i in 1:length(Zr)) {
  power_Zr[i] <- power.individual_Shinichi(mu = model_est_Zr$beta0[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_Zr)) {
  Zr[[i]]$power_Zr <- power_Zr[[i]]
  }

# bias-corrected version
power_c_Zr <- NA
for (i in 1:length(Zr)) {
  power_c_Zr[i] <- power.individual_Shinichi(mu = model_est_Zr$beta0_c3[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of power into corresponding dataset
for (i in 1:length(power_Zr)) {
  Zr[[i]]$power_c_Zr <- power_c_Zr[[i]]
  }


#'*type M error*
# meta-analytic overall mean
power.M_Zr <- NA
for (i in 1:length(Zr)) {
  power.M_Zr[i] <- mapply(error_M,mu = model_est_Zr$beta0[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.M_Zr)) {
  Zr[[i]]$power.M_Zr <- power.M_Zr[[i]]
  }

# bias-corrected version
power.M_c_Zr <- NA
for (i in 1:length(Zr)) {
  power.M_c_Zr[i] <- mapply(error_M,mu = model_est_Zr$beta0_c3[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.M_Zr)) {
  Zr[[i]]$power.M_c_Zr <- power.M_c_Zr[[i]]
  }


#'*type S error*
# meta-analytic overall mean

power.S_Zr <- NA
for (i in 1:length(Zr)) {
  power.S_Zr[i] <- mapply(error_S, mu = model_est_Zr$beta0[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.S_Zr)) {
  Zr[[i]]$power.S_Zr <- power.S_Zr[[i]]
  }

# bias-corrected version
power.S_c_Zr <- NA
for (i in 1:length(Zr)) {
  power.S_c_Zr[i] <- mapply(error_S, mu = model_est_Zr$beta0_c3[i], se = Zr[[i]]$sei) |> list()}

# allocate each set of type S error into corresponding dataset
for (i in 1:length(power.S_Zr)) {
  Zr[[i]]$power.S_c_Zr <- power.S_c_Zr[[i]]
  }


## summary of primary study power ##
#'*two tailed power*
# meta-analytic overall mean
power_summary_Zr <- data.frame(case = names(Zr),
                               Minimum = mapply(summary, sapply(Zr, function(x) x$power_Zr))[1,],      
                               `First quarter` = mapply(summary, sapply(Zr, function(x) x$power_Zr))[2,],
                               Median = mapply(summary, sapply(Zr, function(x) x$power_Zr))[3,],
                               Mean = mapply(summary, sapply(Zr, function(x) x$power_Zr))[4,], 
                               `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power_Zr))[5,], 
                               Maximum = mapply(summary, sapply(Zr, function(x) x$power_Zr))[6,]) 

# bias-corrected version
power_c_summary_Zr <- data.frame(case = names(Zr),
                                 Minimum = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[1,],      
                                 `First quarter` = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[2,],
                                 Median = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[3,],
                                 Mean = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[4,], 
                                 `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[5,], 
                                 Maximum = mapply(summary, sapply(Zr, function(x) x$power_c_Zr))[6,]) 

#'*type M error*
# meta-analytic overall mean
power.M_summary_Zr <-  data.frame(case = names(Zr),
                                  Minimum = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[1,],      
                                  `First quarter` = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[2,],
                                  Median = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[3,],
                                  Mean = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[4,], 
                                  `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[5,], 
                                  Maximum = mapply(summary, sapply(Zr, function(x) x$power.M_Zr))[6,]) 

# bias-corrected version
power.M_c_summary_Zr <-  data.frame(case = names(Zr),
                                    Minimum = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[1,],      
                                    `First quarter` = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[2,],
                                    Median = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[3,],
                                    Mean = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[4,], 
                                    `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[5,], 
                                    Maximum = mapply(summary, sapply(Zr, function(x) x$power.M_c_Zr))[6,]) 


#'*type S error*
# meta-analytic overall mean
power.S_summary_Zr <-  data.frame(case = names(Zr),
                                  Minimum = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[1,],      
                                  `First quarter` = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[2,],
                                  Median = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[3,],
                                  Mean = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[4,], 
                                  `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[5,], 
                                  Maximum = mapply(summary, sapply(Zr, function(x) x$power.S_Zr))[6,]) 

# bias-corrected version
power.S_c_summary_Zr <-  data.frame(case = names(Zr),
                                    Minimum = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[1,],      
                                    `First quarter` = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[2,],
                                    Median = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[3,],
                                    Mean = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[4,], 
                                    `Third quarter` = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[5,], 
                                    Maximum = mapply(summary, sapply(Zr, function(x) x$power.S_c_Zr))[6,]) 


###### aggregation ----
#' [estimate overall power for meta-analysis level power]

# add N and k
N_Zr <- NA
for (i in 1:length(Zr)) {
  N_Zr[i] <- Zr[[i]]$study_ID |> unique() |> length()
  }
model_est_Zr$N <- N_Zr

k_Zr <- NA
for (i in 1:length(Zr)) {
  k_Zr[i] <- Zr[[i]]$obs_ID |> length()
  }
model_est_Zr$k <- k_Zr

#' *two tailed power*
## meta-analytic overall mean
# log
MMA_MA.power_Zr <- lm(log(MA.power) ~ 1, weights = k, data = model_est_Zr)

# original scale
MMA_MA.power_Zr2 <- lm(MA.power ~ 1, weights = k, data = model_est_Zr)
MMA_MA.power_Zr2$coefficients

# this is median
MMA_MA.power_Zr$coefficients  |> exp() 

# this is mean
(MMA_MA.power_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power))) |> exp() 

#confidence interval of median
confint(MMA_MA.power_Zr) |> exp()

# compare residuals
par(mfrow = c(1, 2))
residuals(MMA_MA.power_Zr) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_MA.power_Zr2) |> hist(main = paste("original power"), xlab = "Residual") 


## bias-corrected version
# log
MMA_MA.power_c_Zr <- lm(log(MA.power_c) ~ 1, weights = k, data = model_est_Zr)

# this is median
MMA_MA.power_c_Zr$coefficients  |> exp() 

# this is mean
(MMA_MA.power_c_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power))) |> exp() 
#confidence interval of median
confint(MMA_MA.power_c_Zr) |> exp()


#' *type M error*
## meta-analytic overall mean
# log
MMA_MA.power.M_Zr <- lm(log(MA.power.M) ~ 1, weights = k, data = model_est_Zr)
# this is median
MMA_MA.power.M_Zr$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power.M))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_Zr) |> exp()


## bias-corrected version
# log
MMA_MA.power.M_c_Zr <- lm(log(MA.power.M_c) ~ 1, weights = k, data = model_est_Zr)
# this is median
MMA_MA.power.M_c_Zr$coefficients |> exp() 
# this is mean
(MMA_MA.power.M_c_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power.M_c))) |> exp() 

#confidence interval of median
confint(MMA_MA.power.M_c_Zr) |> exp()


#' *type S roor*
## meta-analytic overall mean
# log
MMA_MA.power.S_Zr <- lm(log(MA.power.S+0.025) ~ 1, weights = k, data = model_est_Zr) # add an offset of 0.025(25%) to avoid ln(0) = infinity 

# this is median
MMA_MA.power.S_Zr$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power.S + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_Zr) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it

## bias-corrected version
# log
MMA_MA.power.S_c_Zr <- lm(log(MA.power.S_c + 0.025) ~ 1, weights = k, data = model_est_Zr) # add an offset of 0.025(25%) to avoid ln(0) = infinity

# this is median
MMA_MA.power.S_c_Zr$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_c_Zr$coefficients + 0.5*var(log(model_est_Zr$MA.power.S_c + 0.025))) |> exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_c_Zr) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it


#' [estimate overall power for primary study level power]

study_ID_Zr <- sapply(Zr, function(x) x$study_ID) |> unlist()
power_Zr <- sapply(Zr, function(x) x$power_Zr) |> unlist()
power.S_Zr <- sapply(Zr, function(x) x$power.S_Zr) |> unlist()
power.M_Zr <- sapply(Zr, function(x) x$power.M_Zr) |> unlist()
power_c_Zr <- sapply(Zr, function(x) x$power_c_Zr) |> unlist()
power.S_c_Zr <- sapply(Zr, function(x) x$power.S_c_Zr) |> unlist()
power.M_c_Zr <- sapply(Zr, function(x) x$power.M_c_Zr) |> unlist()

individual_est_Zr <- data.frame("study_ID_Zr" = study_ID_Zr,
                                "power_Zr" = power_Zr,
                                "power.S_Zr" = power.S_Zr,
                                "power.M_Zr" = power.M_Zr,
                                "power_c_Zr" = power_c_Zr,
                                "power.S_c_Zr" = power.S_c_Zr,
                                "power.M_c_Zr" = power.M_c_Zr
                                )

#' *two tailed power*
## meta-analytic overall mean
# log
MMA_EXP.power_Zr <- lmer(log(power_Zr) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr)

# original scale
MMA_EXP.power_Zr2 <- lmer(power_Zr ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr)
summary(MMA_EXP.power_Zr2)$coefficients[1]

# this is median 
summary(MMA_EXP.power_Zr)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_Zr)$coefficients[1] + 0.5*var(log(individual_est_Zr$power_Zr))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_Zr) |> exp()

# compare residual
par(mfrow = c(1, 2))
residuals(MMA_EXP.power_Zr) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_EXP.power_Zr2) |> hist(main = paste("orignal power"), xlab = "Residual")

## bias-corrected version
# log
MMA_EXP.power_c_Zr <- lmer(log(power_c_Zr) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr)

# this is median 
summary(MMA_EXP.power_c_Zr)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power_c_Zr)$coefficients[1] + 0.5*var(log(individual_est_Zr$power_Zr))) |> exp()

# confidence interval of median
confint(MMA_EXP.power_c_Zr) |> exp()

#' *type M error*
## meta-analytic overall mean
# log
MMA_EXP.power.M_Zr <- lmer(log(power.M_Zr) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr)

# this is median 
summary(MMA_EXP.power.M_Zr)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_Zr)$coefficients[1] + 0.5*var(log(individual_est_Zr$power.M_Zr))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_Zr) |> exp()

## bias-corrected version
# log
MMA_EXP.power.M_c_Zr <- lmer(log(power.M_c_Zr) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr)

# this is median 
summary(MMA_EXP.power.M_c_Zr)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_c_Zr)$coefficients[1] + 0.5*var(log(individual_est_Zr$power.M_Zr))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_c_Zr) |> exp()


#' *type S error*
## meta-analytic overall mean
# log
MMA_EXP.power.S_Zr <- lmer( log(power.S_Zr + 0.025) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr) # add an offset of 0.025 to avoid log(0) = inf

# this is median 
summary(MMA_EXP.power.S_Zr)$coefficients[1] |> exp() - 0.025 # - 0.025 is important
# this is mean
(summary(MMA_EXP.power.S_Zr)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_Zr)$varcor[[1]][[1]] + # sigma^2 for study level
           summary(MMA_EXP.power.S_Zr)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
) |> exp() - 0.025

#confidence interval of median
# if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it
confint(MMA_EXP.power.S_Zr) |> exp() - 0.025 

## bias-corrected version
# log
MMA_EXP.power.S_c_Zr <- lmer( log(power.S_c_Zr + 0.025) ~ 1 + (1 | study_ID_Zr), data = individual_est_Zr) 

# this is median 
summary(MMA_EXP.power.S_c_Zr)$coefficients[1] |> exp() - 0.025 
# this is mean
(summary(MMA_EXP.power.S_c_Zr)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_Zr)$varcor[[1]][[1]] + 
           summary(MMA_EXP.power.S_Zr)$sigma^2) 
  ) |> exp() - 0.025

#confidence interval of median
confint(MMA_EXP.power.S_c_Zr) |> exp() - 0.025 


##### Across 3 effect size measures (lnRR, SMD, Zr) ----
#' [estimate overall power for meta-analysis level power]

#' *two tailed power*
## power to detect meta-analytic overall mean
model_est_all_corrected_original$MA.power <- power.ma_Shinichi(mu = model_est_all_original$beta0, 
                                                               SE =model_est_all_corrected_original$se_beta0)

## power to detect bias-corrected meta-analytic overall mean
model_est_all_corrected_original$MA.power_c <- power.ma_Shinichi(mu = model_est_all_corrected_original$beta0_c3, 
                                                                 SE = model_est_all_corrected_original$se_beta0) # note to still use the unconditional se rather than the conditional se (se of beta0c_3)

## power to detect small-study effect
model_est_all_corrected_original$sse.power <- power.ma_Shinichi(mu = model_est_all_corrected_original$beta1,
                                                                SE = model_est_all_corrected_original$se_beta1)

## power to detect decline effect
model_est_all_corrected_original$de.power <- power.ma_Shinichi(mu = model_est_all_corrected_original$beta2,
                                                               SE = model_est_all_corrected_original$se_beta2)


#' *type M error*
## meta-analytic overall mean
MA.power.M <- NA
for (i in 1:length(model_est_all_corrected_original$case)) {
  MA.power.M[i] <- error_M(mu = model_est_all_corrected_original$beta0[i],
                           se = model_est_all_corrected_original$se_beta0[i],
                           alpha = 0.05, N = 10000) |> unlist()
}

model_est_all_corrected_original$MA.power.M <- MA.power.M

## bias-corrected version
MA.power.M_c <- NA
for (i in 1:length(model_est_all_corrected_original$case)) {
  MA.power.M_c[i] <- error_M(mu = model_est_all_corrected_original$beta0_c3[i],
                             se = model_est_all_corrected_original$se_beta0[i],
                             alpha = 0.05, N = 10000) |> unlist()
}

model_est_all_corrected_original$MA.power.M_c <- MA.power.M_c

#' *type S error*
## meta-analytic overall mean
MA.power.S <- NA
for (i in 1:length(model_est_all_corrected_original$case)) {
  MA.power.S[i] <- error_S(mu = model_est_all_corrected_original$beta0[i],
                           se = model_est_all_corrected_original$se_beta0[i],
                           alpha = 0.05) |> unlist()
}

model_est_all_corrected_original$MA.power.S <- MA.power.S

## bias-corrected version
MA.power.S_c <- NA
for (i in 1:length(model_est_all_corrected_original$case)) {
  MA.power.S_c[i] <- error_S(mu = model_est_all_corrected_original$beta0_c3[i],
                             se = model_est_all_corrected_original$se_beta0[i],
                             alpha=0.05) |> unlist()
}

model_est_all_corrected_original$MA.power.S_c <- MA.power.S_c


#' [estimate overall power for primary study level power]

individual_est_all <- data.frame(study_ID_all = c(individual_est_lnRR$study_ID_lnRR, individual_est_SMD$study_ID_SMD, individual_est_Zr$study_ID_Zr),
                                 power_all = c(individual_est_lnRR$power_lnRR,individual_est_SMD$power_SMD,individual_est_Zr$power_Zr),
                                 power.S_all = c(individual_est_lnRR$power.S_lnRR,individual_est_SMD$power.S_SMD,individual_est_Zr$power.S_Zr),
                                 power.M_all = c(individual_est_lnRR$power.M_lnRR,individual_est_SMD$power.M_SMD,individual_est_Zr$power.M_Zr),
                                 power_c_all = c(individual_est_lnRR$power_c_lnRR,individual_est_SMD$power_c_SMD,individual_est_Zr$power_c_Zr),
                                 power.S_c_all = c(individual_est_lnRR$power.S_c_lnRR,individual_est_SMD$power.S_c_SMD,individual_est_Zr$power.S_c_Zr),
                                 power.M_c_all = c(individual_est_lnRR$power.M_c_lnRR,individual_est_SMD$power.M_c_SMD,individual_est_Zr$power.M_c_Zr))

###### aggregation ----
#' [estimate overall power for meta-analysis level power] 
# combine lnRR, SMD, and Zr
cols_to_drop <- c("base", "study_ID")
model_est_lnRR <- select(model_est_lnRR, -all_of(cols_to_drop))
model_est_all_corrected_original2 <- rbind(model_est_lnRR,model_est_SMD,model_est_Zr)

#' *two tailed power*

## meta-analytic overall mean
# log
MMA_MA.power_all <- lm(log(MA.power) ~ 1, weights = k, data = model_est_all_corrected_original2)

# original scale
MMA_MA.power_all2 <- lm(MA.power ~ 1, weights = k, data = model_est_all_corrected_original2)
MMA_MA.power_all2$coefficients

# this is median
MMA_MA.power_all$coefficients |> exp() 

# this is mean
(MMA_MA.power_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power))) |> exp() 

# confidence interval of median
confint(MMA_MA.power_all) |> exp()

# compare residuals
par(mfrow = c(1, 2))
residuals(MMA_MA.power_all) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_MA.power_all2) |> hist(main = paste("original power"), xlab = "Residual") 

## bias-corrected version
# log
MMA_MA.power_c_all <- lm(log(MA.power_c) ~ 1, weights = k, data = model_est_all_corrected_original2)

# this is median
MMA_MA.power_c_all$coefficients |> exp() 

# this is mean
(MMA_MA.power_c_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power))) |> exp() 

# confidence interval of median
confint(MMA_MA.power_c_all) |> exp()


#' *type M error*

## meta-analytic overall mean
# log
MMA_MA.power.M_all <- lm(log(MA.power.M) ~ 1, weights = k, data = model_est_all_corrected_original2)

# this is median
MMA_MA.power.M_all$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power.M))) |> exp() 

# confidence interval of median
confint(MMA_MA.power.M_all) |> exp()


## bias-corrected version
# log
MMA_MA.power.M_c_all <- lm(log(MA.power.M_c) ~ 1, weights = k, data = model_est_all_corrected_original2)

# this is median
MMA_MA.power.M_c_all$coefficients |> exp() 

# this is mean
(MMA_MA.power.M_c_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power.M_c))) |> exp() 

# confidence interval of median
confint(MMA_MA.power.M_c_all) |> exp()


#' *typ S error*

## meta-analytic overall mean
# log
# add an offset of 0.025(25%) to avoid ln(0) = infinity 
MMA_MA.power.S_all <- lm(log(MA.power.S + 0.025) ~ 1, weights = k, data = model_est_all_corrected_original2) 

# this is median
MMA_MA.power.S_all$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power.S + 0.025))) |> exp() - 0.025

#confidence interval of median
# if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it
confint(MMA_MA.power.S_all) |> exp() - 0.025 

# bias-corrected version
# log
MMA_MA.power.S_c_all <- lm(log(MA.power.S_c + 0.025) ~ 1, weights = k, 
                           data = model_est_all_corrected_original2) 

# this is median
MMA_MA.power.S_c_all$coefficients |> exp() - 0.025

# this is mean
(MMA_MA.power.S_c_all$coefficients + 0.5*var(log(model_est_all_corrected_original2$MA.power.S_c + 0.025))) |>  exp() - 0.025

#confidence interval of median
confint(MMA_MA.power.S_c_all) |> exp() - 0.025 


#' [estimate overall power for primary study level power] 

#' *two tailed power*

## meta-analytic overall mean
# log
MMA_EXP.power_all <- lmer(log(power_all) ~ 1 + (1 | study_ID_all), data = individual_est_all)
# original scale
MMA_EXP.power_all2 <- lmer(power_all ~ 1 + (1 | study_ID_all), data = individual_est_all)
summary(MMA_EXP.power_all2)$coefficients[1]

# this is median 
summary(MMA_EXP.power_all)$coefficients[1] |> exp()
# this is mean
(summary(MMA_EXP.power_all)$coefficients[1] + 0.5*var(log(individual_est_all$power_all))) |> exp()
# confidence interval of median
confint(MMA_EXP.power_all) |> exp()

# compare residual
par(mfrow = c(1, 2))
residuals(MMA_EXP.power_all) |> hist(main = paste("log power"), xlab = "Residual")
residuals(MMA_EXP.power_all2) |> hist(main = paste("orignal power"), xlab = "Residual")

# bias-corrected version
# log
MMA_EXP.power_c_all <- lmer(log(power_c_all) ~ 1 + (1 | study_ID_all), data = individual_est_all)

# this is median 
summary(MMA_EXP.power_c_all)$coefficients[1] |> exp()
# this is mean
(summary(MMA_EXP.power_c_all)$coefficients[1] + 0.5*var(log(individual_est_all$power_all))) |> exp()
# confidence interval of median
confint(MMA_EXP.power_c_all) |> exp()


#' *type M error*

## meta-analytic overall mean
# log
MMA_EXP.power.M_all <- lmer(log(power.M_all) ~ 1 + (1 | study_ID_all), data = individual_est_all)

# this is median 
summary(MMA_EXP.power.M_all)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_all)$coefficients[1] + 0.5*var(log(individual_est_all$power.M_all))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_all) |> exp()

## bias-corrected version
# log
MMA_EXP.power.M_c_all <- lmer(log(power.M_c_all) ~ 1 + (1 | study_ID_all), data = individual_est_all)

# this is median 
summary(MMA_EXP.power.M_c_all)$coefficients[1] |> exp()

# this is mean
(summary(MMA_EXP.power.M_c_all)$coefficients[1] + 0.5*var(log(individual_est_all$power.M_all))) |> exp()

# confidence interval of median
confint(MMA_EXP.power.M_c_all) |> exp()


#' *type S error*

## meta-analytic overall mean
# log
MMA_EXP.power.S_all <- lmer( log(power.S_all + 0.025) ~ 1 + (1 | study_ID_all), data = individual_est_all) # add an offset of 0.025 to avoid log(0) = inf

# this is median 
summary(MMA_EXP.power.S_all)$coefficients[1] |> exp() - 0.025 # - 0.025 is important

# this is mean
(summary(MMA_EXP.power.S_all)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_all)$varcor[[1]][[1]] + # sigma^2 for study level
           summary(MMA_EXP.power.S_all)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
) |> exp() - 0.025

#confidence interval of median
confint(MMA_EXP.power.S_all) |> exp() - 0.025 # if the lower boundary was negative (probably caused by the negative variance), we used 0 to replace it

## bias-corrected version
# log
MMA_EXP.power.S_c_all <- lmer( log(power.S_c_all + 0.025) ~ 1 + (1 | study_ID_all), data = individual_est_all) # add an offset of 0.025 to avoid log(0) = inf

# this is median 
summary(MMA_EXP.power.S_c_all)$coefficients[1] |> exp() - 0.025 # - 0.025 is important

# this is mean
(summary(MMA_EXP.power.S_c_all)$coefficients[1] + 
    0.5*(summary(MMA_EXP.power.S_all)$varcor[[1]][[1]] + # sigma^2 for study level
           summary(MMA_EXP.power.S_all)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
) |> exp() - 0.025

#confidence interval of median
confint(MMA_EXP.power.S_c_all) |> exp() - 0.025

###
# 3. Plots ----
###
## firepower plot - power ----

### lnRR
nrow(model_est_lnRR)
MA_case_lnRR <- paste("lnRR", sep = " ", 1:5)

intercept_lnRR <- data.frame(MA_case = MA_case_lnRR,
                             es = model_est_lnRR$beta0,
                             power = power_summary_lnRR$Median,
                             effect = rep(c("Sampling"), 5),
                             es_cat = rep(c("a.Intercept_es"), 5))


intercept_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                      es = model_est_lnRR$beta0_c3,
                                      power = power_c_summary_lnRR$Median,
                                      effect = rep(c("cSampling"), 5),
                                      es_cat = rep(c("a.Intercept_es"), 5))

MA.power_lnRR <- data.frame(MA_case = MA_case_lnRR,
                            es = model_est_lnRR$beta0,
                            power = model_est_lnRR$MA.power,
                            effect = rep(c("Meta-analysis"), 5),
                            es_cat = rep(c("b.MA.power"), 5))

MA.power_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                     es = model_est_lnRR$beta0_c3,
                                     power = model_est_lnRR$MA.power_c,
                                     effect = rep(c("cMeta-analysis"), 5),
                                     es_cat = rep(c("b.MA.power"), 5))

power_firepower_lnRR <- rbind(intercept_adjusted_lnRR,
                              intercept_lnRR,
                              MA.power_adjusted_lnRR,
                              MA.power_lnRR)

# convert MA_case to factor with desired ordering of levels
power_firepower_lnRR$MA_case <- factor(power_firepower_lnRR$MA_case, levels=rev(MA_case_lnRR))


### SMD
nrow(model_est_SMD)
MA_case_SMD <- paste("SMD", sep = " ", 1:32)

intercept_SMD <- data.frame(MA_case = MA_case_SMD,
                            es = model_est_SMD$beta0,
                            power = power_summary_SMD$Median,
                            effect = rep(c("Sampling"), 32),
                            es_cat = rep(c("a.Intercept_es"), 32))

intercept_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                     es = model_est_SMD$beta0_c3,
                                     power = power_c_summary_SMD$Median,
                                     effect = rep(c("cSampling"), 32),
                                     es_cat = rep(c("a.Intercept_es"), 32))

MA.power_SMD <- data.frame(MA_case = MA_case_SMD,
                           es = model_est_SMD$beta0,
                           power = model_est_SMD$MA.power,
                           effect = rep(c("Meta-analysis"), 32),
                           es_cat = rep(c("b.MA.power"), 32))

MA.power_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                    es = model_est_SMD$beta0_c3,
                                    power = model_est_SMD$MA.power_c,
                                    effect = rep(c("cMeta-analysis"), 32),
                                    es_cat = rep(c("b.MA.power"), 32))

power_firepower_SMD <- rbind(intercept_adjusted_SMD,
                             intercept_SMD,
                             MA.power_adjusted_SMD,
                             MA.power_SMD)

# convert MA_case to factor with desired ordering of levels
power_firepower_SMD$MA_case <- factor(power_firepower_SMD$MA_case, levels=rev(MA_case_SMD))

### Zr
nrow(model_est_Zr)
MA_case_Zr <- paste("Zr", sep = " ", 1:11)

intercept_Zr <- data.frame(MA_case = MA_case_Zr,
                           es = model_est_Zr$beta0,
                           power = power_summary_Zr$Median,
                           effect = rep(c("Sampling"), 11),
                           es_cat = rep(c("a.Intercept_es"), 11))


intercept_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                    es = model_est_Zr$beta0_c3,
                                    power = power_c_summary_Zr$Median,
                                    effect = rep(c("cSampling"), 11),
                                    es_cat = rep(c("a.Intercept_es"), 11))

MA.power_Zr <- data.frame(MA_case = MA_case_Zr,
                          es = model_est_Zr$beta0,
                          power = model_est_Zr$MA.power,
                          effect = rep(c("Meta-analysis"), 11),
                          es_cat = rep(c("b.MA.power"), 11))

MA.power_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                   es = model_est_Zr$beta0_c3,
                                   power = model_est_Zr$MA.power_c,
                                   effect = rep(c("cMeta-analysis"), 11),
                                   es_cat = rep(c("b.MA.power"), 11))


power_firepower_Zr <- rbind(intercept_adjusted_Zr,
                            intercept_Zr,
                            MA.power_adjusted_Zr,
                            MA.power_Zr)

# convert MA_case to factor with desired ordering of levels
power_firepower_Zr$MA_case <- factor(power_firepower_Zr$MA_case, levels=rev(MA_case_Zr))

## all
power_firepower_all <- rbind(power_firepower_lnRR,power_firepower_SMD,power_firepower_Zr)

MA_case_all <- c(
  paste("lnRR", 1:5),
  paste("SMD", 1:32),
  paste("Zr",  1:11)
)

power_firepower_all <- rbind(
  power_firepower_lnRR,
  power_firepower_SMD,
  power_firepower_Zr
)

power_firepower_all$MA_case <- factor(power_firepower_all$MA_case,
                                      levels = rev(MA_case_all)) 

power_firepower_all$effect <- factor(power_firepower_all$effect,
                                 levels = c("Sampling", "cSampling", "Meta-analysis", "cMeta-analysis"))

firepower_plot_all <- ggplot(power_firepower_all) +
  geom_tile(aes(x = effect, y = MA_case, fill = power),
            width = 0.95, height = 0.5) +
  scale_fill_gradient(
    name = "Power",
    low  = "#FF7F50",
    high = "#FFFFFF",
    limits   = c(0, 1),
    na.value = "#FFFFFF"
  ) +
  facet_grid(~ es_cat, scale = "free_x", space = "free_x") +
  # scale_y_discrete(limits = rev) +       
  scale_x_discrete(position = "top") +   
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 12, colour = "#000000",
                               margin = margin(t = 0, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(size = 12, colour = "#000000"),
    axis.title.x = element_text(size = 12, colour = "#000000", face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width  = unit(0.5, "in"),
    legend.key.height = unit(0.4, "in"),
    legend.position   = "bottom",
    panel.background = element_rect(fill = "#FFFFFF", colour = NA),
    plot.background  = element_rect(fill = "#FFFFFF", colour = NA),
    panel.grid       = element_blank()
  ) +
  labs(x = "", y = "") +              
  labs(title = "Statistical power")   

firepower_plot_all


## firepower plot - type M error ----
### lnRR
intercept_lnRR <- data.frame(MA_case = MA_case_lnRR,
                             es = model_est_lnRR$beta0,
                             power = power.M_summary_lnRR$Median,
                             effect = rep(c("Sampling"), 5),
                             es_cat = rep(c("a.Intercept_es"), 5))

intercept_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                      es = model_est_lnRR$beta0_c3,
                                      power = power.M_c_summary_lnRR$Median,
                                      effect = rep(c("cSampling"), 5),
                                      es_cat = rep(c("a.Intercept_es"), 5))

MA.power_lnRR <- data.frame(MA_case = MA_case_lnRR,
                            es = model_est_lnRR$beta0,
                            power = model_est_lnRR$MA.power.M,
                            effect = rep(c("Meta-analysis"), 5),
                            es_cat = rep(c("b.MA.power"), 5))

MA.power_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                     es = model_est_lnRR$beta0_c3,
                                     power = model_est_lnRR$MA.power.M_c,
                                     effect = rep(c("cMeta-analysis"), 5),
                                     es_cat = rep(c("b.MA.power"), 5))

power.M_firepower_lnRR <- rbind(intercept_adjusted_lnRR,
                                intercept_lnRR,
                                MA.power_adjusted_lnRR,
                                MA.power_lnRR)

# convert MA_case to factor with desired ordering of levels
power.M_firepower_lnRR$MA_case <- factor(power_firepower_lnRR$MA_case, levels = rev(MA_case_lnRR))

### SMD

MA_case_SMD <- paste("SMD", sep = " ", 1:32) 

intercept_SMD <- data.frame(MA_case = MA_case_SMD,
                            es = model_est_SMD$beta0,
                            power = power.M_summary_SMD$Median,
                            effect = rep(c("Sampling"), 32),
                            es_cat = rep(c("a.Intercept_es"), 32))

intercept_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                     es = model_est_SMD$beta0_c3,
                                     power = power.M_c_summary_SMD$Median,
                                     effect = rep(c("cSampling"), 32),
                                     es_cat = rep(c("a.Intercept_es"), 32))

MA.power_SMD <- data.frame(MA_case = MA_case_SMD,
                           es = model_est_SMD$beta0,
                           power = model_est_SMD$MA.power.M,
                           effect = rep(c("Meta-analysis"), 32),
                           es_cat = rep(c("b.MA.power"), 32))

MA.power_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                    es = model_est_SMD$beta0_c3,
                                    power = model_est_SMD$MA.power.M_c,
                                    effect = rep(c("cMeta-analysis"), 32),
                                    es_cat = rep(c("b.MA.power"), 32))

power.M_firepower_SMD <- rbind(intercept_adjusted_SMD,
                               intercept_SMD,
                               MA.power_adjusted_SMD,
                               MA.power_SMD)

# convert MA_case to factor with desired ordering of levels
power.M_firepower_SMD$MA_case <- factor(power.M_firepower_SMD$MA_case, levels=rev(MA_case_SMD))

### Zr

MA_case_Zr <- paste("Zr", sep = " ", 1:11) # nrow(model_est_Zr)

intercept_Zr <- data.frame(MA_case = MA_case_Zr,
                           es = model_est_Zr$beta0,
                           power = power.M_summary_Zr$Median,
                           effect = rep(c("Sampling"), 11),
                           es_cat = rep(c("a.Intercept_es"), 11))

intercept_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                    es = model_est_Zr$beta0_c3,
                                    power = power.M_c_summary_Zr$Median,
                                    effect = rep(c("cSampling"), 11),
                                    es_cat = rep(c("a.Intercept_es"), 11))

MA.power_Zr <- data.frame(MA_case = MA_case_Zr,
                          es = model_est_Zr$beta0,
                          power = model_est_Zr$MA.power.M,
                          effect = rep(c("Meta-analysis"), 11),
                          es_cat = rep(c("b.MA.power"), 11))

MA.power_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                   es = model_est_Zr$beta0_c3,
                                   power = model_est_Zr$MA.power.M_c,
                                   effect = rep(c("cMeta-analysis"), 11),
                                   es_cat = rep(c("b.MA.power"), 11))


power.M_firepower_Zr <- rbind(intercept_adjusted_Zr,
                              intercept_Zr,
                              MA.power_adjusted_Zr,
                              MA.power_Zr)

# convert MA_case to factor with desired ordering of levels
power.M_firepower_Zr$MA_case <- factor(power.M_firepower_Zr$MA_case, levels=rev(MA_case_Zr))


# all
firepower.M_all <- rbind(power.M_firepower_lnRR,power.M_firepower_SMD,power.M_firepower_Zr)

firepower.M_all$MA_case <- factor(firepower.M_all$MA_case, levels = rev(MA_case_all))

firepower.M_all$effect <- factor(firepower.M_all$effect,
                                 levels = c("Sampling", "cSampling", "Meta-analysis", "cMeta-analysis"))

firepower.M_plot_all <- ggplot(data = firepower.M_all) +
  geom_tile(aes(x = effect, y = MA_case, fill = power, width = 0.95, height = 0.5)) +
  scale_fill_gradient(
    name = "Type M error",
    low = "#FFFFFF", 
    #mid = "white",
    high = "#00868B", 
    #midpoint = 0.5, # if using scale_fill_gradient2, it needs mid color and midpoint
    limits = c(0,20),
    guide = "colourbar",
  ) + 
  facet_grid( ~ es_cat, scale = 'free_x', space = "free_x") +
  # scale_y_discrete(limits = rev) +
  theme(strip.text.x = element_blank()) + 
  labs(x ="", y = "") + 
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 12, colour = "#000000", 
                               margin = margin(t = 0, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(size = 12, colour = "#000000"),
    axis.title.x = element_text(size = 12, colour = "#000000", face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.5, "in"),
    legend.key.height = unit(0.4, "in"),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#FFFFFF", colour = NA),
    plot.background = element_rect(fill = "#FFFFFF", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_x_discrete(position = "top") + labs(title = "Type M error")

firepower.M_plot_all


## firepower plot - S error ----
### lnRR

intercept_lnRR <- data.frame(MA_case = MA_case_lnRR,
                             es = model_est_lnRR$beta0,
                             power = power.S_summary_lnRR$Median,
                             effect = rep(c("Sampling"), 5),
                             es_cat = rep(c("a.Intercept_es"),5))


intercept_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                      es = model_est_lnRR$beta0_c3,
                                      power = power.S_c_summary_lnRR$Median,
                                      effect = rep(c("cSampling"), 5),
                                      es_cat = rep(c("a.Intercept_es"), 5))

MA.power_lnRR <- data.frame(MA_case = MA_case_lnRR,
                            es = model_est_lnRR$beta0,
                            power = model_est_lnRR$MA.power.S,
                            effect = rep(c("Meta-analysis"), 5),
                            es_cat = rep(c("b.MA.power"), 5))

MA.power_adjusted_lnRR <- data.frame(MA_case = MA_case_lnRR,
                                     es = model_est_lnRR$beta0_c3,
                                     power = model_est_lnRR$MA.power.S_c,
                                     effect = rep(c("cMeta-analysis"), 5),
                                     es_cat = rep(c("b.MA.power"), 5))


power.S_firepower_lnRR <- rbind(intercept_adjusted_lnRR,
                                intercept_lnRR,
                                MA.power_adjusted_lnRR,
                                MA.power_lnRR)

# convert MA_case to factor with desired ordering of levels
power.S_firepower_lnRR$MA_case <- factor(power_firepower_lnRR$MA_case, levels=rev(MA_case_lnRR))


### SMD
MA_case_SMD <- paste("SMD", sep = " ", 1:32) # nrow(model_est_SMD)

intercept_SMD <- data.frame(MA_case = MA_case_SMD,
                            es = model_est_SMD$beta0,
                            power = power.S_summary_SMD$Median,
                            effect = rep(c("Sampling"), 32),
                            es_cat = rep(c("a.Intercept_es"), 32))


intercept_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                     es = model_est_SMD$beta0_c3,
                                     power = power.S_c_summary_SMD$Median,
                                     effect = rep(c("cSampling"), 32),
                                     es_cat = rep(c("a.Intercept_es"), 32))

MA.power_SMD <- data.frame(MA_case = MA_case_SMD,
                           es = model_est_SMD$beta0,
                           power = model_est_SMD$MA.power.S,
                           effect = rep(c("Meta-analysis"), 32),
                           es_cat = rep(c("b.MA.power"), 32))

MA.power_adjusted_SMD <- data.frame(MA_case = MA_case_SMD,
                                    es = model_est_SMD$beta0_c3,
                                    power = model_est_SMD$MA.power.S_c,
                                    effect = rep(c("cMeta-analysis"), 32),
                                    es_cat = rep(c("b.MA.power"), 32))


power.S_firepower_SMD <- rbind(intercept_adjusted_SMD,
                               intercept_SMD,
                               MA.power_adjusted_SMD,
                               MA.power_SMD)

# convert MA_case to factor with desired ordering of levels
power.S_firepower_SMD$MA_case <- factor(power_firepower_SMD$MA_case, levels = rev(MA_case_SMD))

### Zr

MA_case_Zr <- paste("Zr", sep = " ", 1:11) # nrow(model_est_Zr)

intercept_Zr <- data.frame(MA_case = MA_case_Zr,
                           es = model_est_Zr$beta0,
                           power = power.S_summary_Zr$Median,
                           effect = rep(c("Sampling"),11),
                           es_cat = rep(c("a.Intercept_es"),11))


intercept_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                    es = model_est_Zr$beta0_c3,
                                    power = power.S_c_summary_Zr$Median,
                                    effect = rep(c("cSampling"),11),
                                    es_cat = rep(c("a.Intercept_es"),11))

MA.power_Zr <- data.frame(MA_case = MA_case_Zr,
                          es = model_est_Zr$beta0,
                          power = model_est_Zr$MA.power.S,
                          effect = rep(c("Meta-analysis"),11),
                          es_cat = rep(c("b.MA.power"),11))

MA.power_adjusted_Zr <- data.frame(MA_case = MA_case_Zr,
                                   es = model_est_Zr$beta0_c3,
                                   power = model_est_Zr$MA.power.S_c,
                                   effect = rep(c("cMeta-analysis"),11),
                                   es_cat = rep(c("b.MA.power"),11))


power.S_firepower_Zr <- rbind(intercept_adjusted_Zr,
                              intercept_Zr,
                              MA.power_adjusted_Zr,
                              MA.power_Zr)

# convert MA_case to factor with desired ordering of levels
power.S_firepower_Zr$MA_case <- factor(power_firepower_Zr$MA_case, levels = rev(MA_case_Zr))

### all
firepower.S_all <- rbind(power.S_firepower_lnRR,power.S_firepower_SMD,power.S_firepower_Zr)

firepower.S_all$MA_case <- factor(firepower.S_all$MA_case, levels = rev(MA_case_all))


firepower.S_all$effect <- factor(firepower.S_all$effect,
                                 levels = c("Sampling", "cSampling", "Meta-analysis", "cMeta-analysis"))

firepower.S_plot_all <- ggplot(data = firepower.S_all) +
  geom_tile(aes(x = effect, y = MA_case, fill = power, width = 0.95, height = 0.5)) +
  scale_fill_gradient(
    name = "Type S error",
    low = "#FFFFFF",
    #mid = "white",
    high = "#CDCD00", 
    #midpoint = 0.5, # if using scale_fill_gradient2, it needs mid color and midpoint
    limits = c(0,1),
    guide = "colourbar",
  ) + 
  facet_grid(~ es_cat, scale = "free_x", space = "free_x") +
  # scale_y_discrete(limits = rev) +
  theme(strip.text.x = element_blank()) + 
  labs(x ="", y = "") + 
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 12, colour = "#000000", 
                               margin = margin(t = 0, r = 0, b = 5, l = 0)),
    axis.text.y = element_text(size = 12, colour = "#000000"),
    axis.title.x = element_text(size = 12, colour = "#000000", face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.5, "in"),
    legend.key.height = unit(0.4, "in"),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#FFFFFF", colour = NA),
    plot.background = element_rect(fill = "#FFFFFF", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_x_discrete(position = "top") + labs(title = "Type S error")

firepower.S_plot_all

## check all figures ----
firepower_plot_all
firepower.M_plot_all
firepower.S_plot_all

nrow(firepower.S_all$MA_case)


# significant -> non-significant ----
## originally significant
initially_significant <- subset(model_est_all_original, pval_beta0 < 0.05)

## changed
lost_significance <- subset(initially_significant, pval_beta0_c >= 0.05)

num_initially_significant <- nrow(initially_significant)
num_lost_significance <- nrow(lost_significance)

percentage_lost <- (num_lost_significance / num_initially_significant) * 100

print(paste("Initially significant:", num_initially_significant))
print(paste("Lost significance:", num_lost_significance))
print(paste("Percentage:", percentage_lost))

