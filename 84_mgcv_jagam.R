
library(mgcv)
library(purrr)
library(dplyr)
library(ggplot2)


?jagam

## the following illustrates a typical workflow. To run the 
## 'Not run' code you need rjags (and JAGS) to be installed.
require(mgcv)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Simulate data ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

set.seed(2) ## simulate some data... 
n <- 400
dat <- gamSim(1,n=n,dist="normal",scale=2)

par(mfrow = c(2,2), mar = c(4,5,2,1))
plot(y~x0, data = dat)
plot(y~x1, data = dat)
plot(y~x2, data = dat)
plot(y~x3, data = dat)

## regular gam fit for comparison...
b0 <- gam(y ~ s(x0)+s(x1)+ s(x2)+s(x3),data=dat,method="REML")

plot(b0)


## Set directory and file name for file containing jags code.
## In real use you would *never* use tempdir() for this. It is
## only done here to keep CRAN happy, and avoid any chance of
## an accidental overwrite. Instead you would use
## setwd() to set an appropriate working directory in which
## to write the file, and just set the file name to what you
## want to call it (e.g. "test.jags" here). 

# jags.file <- paste(tempdir(),"/test.jags",sep="") 

## Set up JAGS code and data. In this one might want to diagonalize
## to use conjugate samplers. Usually call 'setwd' first, to set
## directory in which model file ("test.jags") will be written.

# All variables 
# jd <- jagam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat,file=jags.file,
#             sp.prior="gamma",diagonalize=TRUE)
# edit(file = jags.file)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Testing jagam ----
#
# Analyse effect of just one variabl, x2 (for simplicity)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Just once: create directory for files 
# dir.create("C:/Data/R_test/84_mgcv_jagam")


jags.file_cr <- paste("C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_cr_orig.txt") 
jd <- jagam(y ~ s(x2, bs="cr"), data=dat,file=jags.file_cr,
            sp.prior="gamma", diagonalize=TRUE)


jags.file <- paste("C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp_orig.txt") 
jd <- jagam(y ~ s(x2, bs="tp"), data=dat,file=jags.file,
            sp.prior="gamma", diagonalize=TRUE)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Set up models/data ----
#
# Thin plate splines, fixed max knots  
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jags.file <- paste("C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tpX_orig.txt") 

k <- 3
jags.file_tp3 <- sub("X", k, jags.file)
jd_tp3 <- jagam(y ~ s(x2, bs="tp", k=k), data=dat,file=sub("X", k, jags.file),
                sp.prior="gamma", diagonalize=TRUE)

# If you want to specify placemets of the knots:
# knots = list(x2 = seq(0, 1, length = k)))

# Three places to change 3 to 4: k, jags.file_tp4, jd_tp4
k <- 4
jags.file_tp4 <- sub("X", k, jags.file)
jd_tp4 <- jagam(y ~ s(x2, bs="tp", k=k), data=dat,file=sub("X", k, jags.file),
                sp.prior="gamma", diagonalize=TRUE)

k <- 5
jags.file_tp5 <- sub("X", k, jags.file)
# debugonce(jagam)
jd_tp5 <- jagam(y ~ s(x2, bs="tp", k=k), data=dat,file=sub("X", k, jags.file),
                sp.prior="gamma", diagonalize=TRUE)


# jd <- jagam(y ~ s(x2, bs="cr", fx = TRUE), data=dat,file=jags.file,
#             sp.prior="gamma", diagonalize=TRUE)
# subscript out of bounds


## In normal use the model in "test.jags" would now be edited to add 
## the non-standard stochastic elements that require use of JAGS....


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Estimate one of the models above ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jm <- jags.model(jags.file_tp4, data=jd_tp4$jags.data, inits=jd_tp4$jags.ini, n.chains=2)  # changed from n.chains=1 
list.samplers(jm)
update(jm, n.iter=10000)
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Gelman test 
sam_mcmc <- as.mcmc(sam)
# gelman.diag(sam)  - doesn't work
gelman.diag(sam$b)
gelman.diag(sam$rho)

# Raftery test
sam_mcmc <- mcmcUpgrade(sam_mcmc)
raftery.diag(sam_mcmc)
X <- raftery.diag(sam_mcmc)
str(X)
X$resmatrix[2]

# Make into GAM object and plot normalized effect
jam <- sim2jam(sam, jd_tp4$pregam)
plot(jam, pages=1)
# plot(jam, pages=1, residuals = TRUE) - not supported  

# Plot effect on original scale
get_plotdata <- function(jam, n = 40, x_name = "x2"){
  df_fit <- data.frame(x = seq(0,1,length = n))
  names(df_fit) <- x_name
  pred <- predict(jam, newdata = df_fit, se.fit = TRUE)
  df_fit$y <- pred$fit
  df_fit$y_lo <- pred$fit - 1.96*pred$se.fit
  df_fit$y_hi <- pred$fit + 1.96*pred$se.fit
  df_fit
}
df_fit <- get_plotdata(jam)
plot(dat$x2, dat$y)
lines(y ~ x2, data = df_fit)
lines(y_lo ~ x2, data = df_fit, lty = "dashed")
lines(y_hi ~ x2, data = df_fit, lty = "dashed")

# DIC
get_dic <- function(jagsmodel){
  dic.pd <- rjags::dic.samples(model = jagsmodel, n.iter=1000, type="pD")
  # Select the observations for which we got penalties
  dic.sel.pd <- !is.nan(dic.pd$penalty )
  # Get penalties and deviances for those
  pd <- dic.pd$penalty[dic.sel.pd]
  deviance <- dic.pd$deviance[dic.sel.pd]
  # Calculate DIC
  dic <- sum(deviance) + sum(pd)
  dic
}
get_dic(jm)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# For all non-linear models ---- 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

library(purrr)

input_list <- list(
  file = list(jags.file_tp3, jags.file_tp4, jags.file_tp5),
  jd = list(jd_tp3, jd_tp4, jd_tp5))

jm_list <- pmap(input_list, 
                function(file, jd) jags.model(file, data=jd$jags.data, inits=jd$jags.ini, n.chains=2)  # changed from n.chains=1 
)

# Update
walk(jm_list, update, n.iter = 10000)

# Get samples
sam_list <- map(jm_list, jags.samples, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Make into GAM objects
jam_list <- map2(sam_list, input_list$jd, ~sim2jam(sam = .x, pregam = .y$pregam))

# Plot models
par(mfrow = c(2,2), mar = c(4,5,2,1))
for (i in 1:3) plot(jam_list[[i]])

# DIC
map_dbl(jm_list, get_dic)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Linear model ---- 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jags.file_linear <- "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_linear.txt"  
data_linear <- list(X = dat$x2, y = dat$y, n = nrow(dat))

jm_linear <- jags.model(jags.file_linear, data=data_linear, n.chains=2)  # changed from n.chains=1 
update(jm_linear, n.iter = 10000)

sam_linear <- jags.samples(jm_linear, c("intercept", "slope", "scale"), n.iter=10000, thin=10)

# Coefficients
intercept_q <- quantile(as.numeric(sam_linear$intercept), c(0.025, 0.5, 0.975))
slope_q <- quantile(as.numeric(sam_linear$slope), c(0.025, 0.5, 0.975))
sigma_q <- quantile(as.numeric(sam_linear$scale), c(0.025, 0.5, 0.975))

# Get predicted lines
dims <- dim(sam_linear$intercept)
input <- list(i = 1:dims[2], j = 1:dims[3])
idx <- expand.grid(input)
linearfit_n_samples <- 1000
idx <- input_comb[sample(nrow(input_comb), linearfit_n_samples),]  # pick 1000 random 
idx[1,]

get_linearfit <- function(i, xvalues, samples){
  data.frame(x = xvalues, y = samples$intercept[1,idx[i,1],idx[i,2]] + samples$slope[1,idx[i,1],idx[i,2]]*xvalues)
}
# test
# get_linearfit(1, seq(0,1,length=10), sam_linear)

linearfit_samples <- seq_len(linearfit_n_samples) %>% 
  map_dfr(
    get_linearfit, 
    xvalues = seq(0, 1, length = 20), samples = sam_linear)

linearfit_quantiles <- linearfit_samples %>%
  dplyr::group_by(x) %>%
  dplyr::summarise(
    y_lo = quantile(y, 0.025),
    y_med = quantile(y, 0.5),
    y_hi = quantile(y, 0.975)
  )

# Plot fit and 95% confidence interval
ggplot(linearfit_quantiles, aes(x)) +
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi), fill = "lightblue") +
  geom_abline(intercept = intercept_q["50%"], slope = slope_q["50%"])

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# For linear plus non-linear models ---- 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

library(purrr)

# GAM models
input_list_gam <- list(
  file = list(jags.file_tp3, jags.file_tp4, jags.file_tp5),
  jd = list(jd_tp3, jd_tp4, jd_tp5))

jm_list_gam <- pmap(input_list_jagam, 
                    function(file, jd) jags.model(file, data=jd$jags.data, inits=jd$jags.ini, n.chains=2)  # changed from n.chains=1 
)

# Update
walk(jm_list, update, n.iter = 10000)

# Name
names(jm_list_gam) <- c("gam, df = 2", "gam, df = 3", "gam, df = 4")

# Make jm_list for linear + gam models    
jm_list = append(list(Linear = jm_linear), jm_list_gam)

# Get samples (gam only)
sam_list <- map(jm_list[-1], jags.samples, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Make into GAM objects (linear + gam)
jam_list <- map2(sam_list, input_list$jd, ~sim2jam(sam = .x, pregam = .y$pregam))

# Plot GAM models  
par(mfrow = c(2,2), mar = c(4,5,2,1))
for (i in 1:3) plot(jam_list[[i]])

# DIC (linear + gam)
map_dbl(jm_list, get_dic)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Left-censored data, try 1 ---- 
#
# Tried to use sim2jam() to 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

ggplot(dat, aes(x2, y)) +
  geom_point()

# Create censored data
thresh <- 4
dat_cens <- dat[c("x2","y")]
dat_cens$y_orig <- dat_cens$y
sel <- dat_cens$y > thresh
dat_cens$y[!sel] <- NA
dat_cens$cut <- thresh
dat_cens$cut[sel] <- NA
dat_cens$uncensored <- 0
dat_cens$uncensored[sel] <- 1

# Order file with uncensored data first
dat_cens_ordered <- bind_rows(
  dat_cens[sel,],
  dat_cens[!sel,]
)

jags.file_tp5_lc_try1 <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_try1.txt"
jags.file_tp5_orig <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_orig.txt"
jags.file_tp5_forfit <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_forfit.txt"
# edit(file=jags.file_tp5_lc)
# debugonce(jagam)
jd_tp5_lc <- jagam(y ~ s(x2, bs="tp", k=k), 
                   data = dat_cens_ordered, 
                   file = jags.file_tp5_orig,   # this file will be overwritten (was used as basis for "_leftcens" file)
                   sp.prior = "gamma", 
                   diagonalize = TRUE)

# Modify jags.data object
str(jd_tp5_lc$jags.data, 1)
n <- sum(sel)
m <- sum(!sel)
jd_tp5_lc$jags.data$n <- n
jd_tp5_lc$jags.data$m <- m
jd_tp5_lc$jags.data$Z <- c(rep(0,n), rep(1, m))
jd_tp5_lc$jags.data$cut <- dat_cens_ordered$cut[!sel]

jm <- jags.model(jags.file_tp5_lc_try1, 
                 data=jd_tp5_lc$jags.data, 
                 inits=jd_tp5_lc$jags.ini, 
                 n.chains=2)  # changed from n.chains=1 
list.samplers(jm)
update(jm, n.iter=10000)
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Make into GAM object and plot normalized effect
jam <- sim2jam(sam, jd_tp5_lc$pregam)   # doesn't work - "coefficient simulation data is missing"
plot(jam, pages=1)
# plot(jam, pages=1, residuals = TRUE) - not supported  

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Left-censored data, try 2 ---- 
#
# Trying to get the fitted line by sending the x's for the fit to
#   a second jagam call, and using the X matrix form resulting jags.data
#   for the 'X_fit' matrix in the 
# This approach didn't work
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jags.file_tp5_lc <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens.txt"
jags.file_tp5_lc_try2 <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_try2.txt"
jags.file_tp5_orig <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_orig.txt"
jags.file_tp5_forfit <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_forfit.txt"
# edit(file=jags.file_tp5_lc)
# debugonce(jagam)

# Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
# Make (2) jd_tp5_lc$jags.data (will be manpulated below)
jd_tp5_lc <- jagam(y ~ s(x2, bs="tp", k=k), 
                   data = dat_cens_ordered, 
                   file = jags.file_tp5_orig,   # this file will be overwritten (was used as basis for "_leftcens" file)
                   sp.prior = "gamma", 
                   diagonalize = TRUE)

# Modify jags.data object
str(jd_tp5_lc$jags.data, 1)
n <- sum(dat_cens_ordered$uncensored == 1)
m <- sum(dat_cens_ordered$uncensored == 0)
jd_tp5_lc$jags.data$n <- n
jd_tp5_lc$jags.data$m <- m
jd_tp5_lc$jags.data$Z <- c(rep(0,n), rep(1, m))
jd_tp5_lc$jags.data$cut <- dat_cens_ordered$cut[dat_cens_ordered$uncensored == 0]

# Make 'X matrix" for the fitted data  
dat_for_fit <- data.frame(
  x2 = seq(0, 1, length = 30),
)
# Use ordinary gam (on uncensored data) to add the y
# (Reason: perhaps the y value affects some sort of normalization)  
dat_for_fit$y <- predict(
  mgcv::gam(y~x2, data = subset(dat_cens_ordered, !is.na(y))),
  newdata = dat_for_fit)

jd_tp5_lc_forfit <- jagam(y ~ s(x2, bs="tp", k=k), 
                          data = dat_for_fit, 
                          file = jags.file_tp5_forfit,   # this file will be overwritten (and not used)
                          sp.prior = "gamma", 
                          diagonalize = TRUE)
jd_tp5_lc_forfit$jags.data$X

# Use this in the first 'jags.data' object
jd_tp5_lc$jags.data$X_fit <- jd_tp5_lc_forfit$jags.data$X

# Check
str(jd_tp5_lc$jags.data, 1)

# Run model
# edit(file = jags.file_tp5_lc)
jm <- jags.model(jags.file_tp5_lc, 
                 data=jd_tp5_lc$jags.data, 
                 inits=jd_tp5_lc$jags.ini, 
                 n.chains=2)  # changed from n.chains=1 
list.samplers(jm)
update(jm, n.iter=10000)

# Sample main varables  
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Sample mu_fit  
sam_fit <- jags.samples(jm, c("mu_fit"), n.iter=5000, thin=5)

sam_fit$mu_fit %>% str()

fitted_n <- dim(sam_fit$mu_fit)[1]
qs <- c(0.025, 0.5, 0.975)
fitted <- seq_len(fitted_n) %>% map_dfr(~quantile(sam_fit$mu_fit[.x,,], qs))
names(fitted) <- paste0("Q_", qs)
fitted <- data.frame(x2 = dat_for_fit$x2, fitted)

ggplot(fitted, aes(x2))+
  geom_ribbon(aes(ymin = Q_0.025, ymax = Q_0.975), fill = "lightblue") +
  geom_line(aes(y = Q_0.5)) +
  geom_point(data = dat_cens_ordered, aes(x2, y)) +
  geom_point(data = dat_cens_ordered, aes(x2, cut), shape = 4)
# That didn't look good....


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Left-censored, working version ---- 
#
# Strategy: add the x's to be used for the fitted line
#   to the end of the data
# Set n and m (number of uncensored and censored data) so that these
#   last lines are never touched, but they are included in the first line
#   which calculates mu (expected values)
#
# This worked :-)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jags.file_tp5_lc <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens.txt"
jags.file_tp5_orig <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_orig.txt"
jags.file_tp5_forfit <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_forfit.txt"
# edit(file=jags.file_tp5_lc)
# debugonce(jagam)

# Make 'X matrix" for the fitted data  
dat_for_fit <- data.frame(
  x2 = seq(0, 1, length = 30),
)
dat_cens_ordered[sel,] %>% View

# dat_cens_ordered1 = data file
# dat_cens_ordered2 = with addition for daa for fitted line (3o points)
dat_cens_ordered1 <- dat_cens_ordered   
dat_cens_ordered2 <- bind_rows(
  dat_cens_ordered,
  dat_for_fit
)

nrow(dat_cens_ordered1)
nrow(dat_cens_ordered2)

# Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
# Make (2) jd_tp5_lc$jags.data (will be manpulated below)
jd_tp5_lc <- jagam(y ~ s(x2, bs="tp", k=k), 
                   data = dat_cens_ordered2,    # file no. 2 here
                   file = jags.file_tp5_orig,   # this file will be overwritten (was used as basis for "_leftcens" file)
                   sp.prior = "gamma", 
                   diagonalize = TRUE)

# Modify jags.data object
str(jd_tp5_lc$jags.data, 1)
n <- sum(dat_cens_ordered1$uncensored %in% 1)   # file no. 1 here 
m <- sum(dat_cens_ordered1$uncensored %in% 0)   # file no. 1 here
jd_tp5_lc$jags.data$n <- n   # - makes sure only these data are use for the likelihood
jd_tp5_lc$jags.data$m <- m   #    - " -
jd_tp5_lc$jags.data$Z <- c(rep(0,n), rep(1, m))
jd_tp5_lc$jags.data$cut <- dat_cens_ordered1$cut[dat_cens_ordered1$uncensored %in% 0]

# Check
str(jd_tp5_lc$jags.data, 1)

# View(jd_tp5_lc$jags.data$X)

# Run model
# edit(file = jags.file_tp5_lc)
jm <- jags.model(jags.file_tp5_lc, 
                 data=jd_tp5_lc$jags.data, 
                 inits=jd_tp5_lc$jags.ini, 
                 n.chains=2)  # changed from n.chains=1 
list.samplers(jm)
update(jm, n.iter=10000)

# Sample main varables  
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Sample varaibles that have been inserted only to get the fitted line
sam_fit <- jags.samples(jm, c("mu[401:430]"), n.iter=2000, thin=2)

# Make 'fitted'  
fitted_n <- dim(sam_fit[[1]])[1]
qs <- c(0.025, 0.5, 0.975)
fitted <- seq_len(fitted_n) %>% map_dfr(~quantile(sam_fit[[1]][.x,,], qs))
names(fitted) <- paste0("Q_", qs)
fitted <- data.frame(x2 = dat_for_fit$x2, fitted)

ggplot(fitted, aes(x2))+
  geom_ribbon(aes(ymin = Q_0.025, ymax = Q_0.975), fill = "lightblue") +
  geom_line(aes(y = Q_0.5)) +
  geom_point(data = dat_cens_ordered, aes(x2, y)) +
  geom_point(data = dat_cens_ordered, aes(x2, cut), shape = 4)
# Worked :-)

get_dic(jm)

effectiveSize(as.mcmc.list(sam_fit$`mu[401:430]`))

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Left-censored, testing by fitting only to the uncensored data ---- 
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

jags.file_tp5_orig <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_orig.txt"

dat_cens_fortest <- dat_cens_ordered1 %>%
  filter(uncensored %in% 1)

nrow(dat_cens_ordered1)
nrow(dat_cens_fortest)

# Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
# Make (2) jd_tp5_lc$jags.data (will be manpulated below)
jd_tp5_lc <- jagam(y ~ s(x2, bs="tp", k=k), 
                   data = dat_cens_fortest,    # file no. 2 here
                   file = jags.file_tp5_orig,   # this file will be overwritten (was used as basis for "_leftcens" file)
                   sp.prior = "gamma", 
                   diagonalize = TRUE)

# Run model
# edit(file = jags.file_tp5_lc)
jm <- jags.model(jags.file_tp5_orig, 
                 data = jd_tp5_lc$jags.data, 
                 inits = jd_tp5_lc$jags.ini, 
                 n.chains = 2)  # changed from n.chains=1 
list.samplers(jm)
update(jm, n.iter=10000)

# Sample main varables  
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Make into GAM object and plot normalized effect
jam <- sim2jam(sam, jd_tp5_lc$pregam)
plot(jam, pages=1)

df_fit_test <- get_plotdata(jam)


ggplot(fitted, aes(x2))+
  geom_ribbon(aes(ymin = Q_0.025, ymax = Q_0.975), fill = "lightblue", color = "blue") +
  geom_ribbon(data = df_fit_test, aes(ymin = y_lo, ymax = y_hi), fill = "red", color = "red3", alpha = 0.5) +
  geom_line(aes(y = Q_0.5)) +
  geom_point(data = dat_cens_ordered, aes(x2, y)) +
  geom_point(data = dat_cens_ordered, aes(x2, cut), shape = 4)
# It does differ


jags.file_tp5_orig_tp5_ver1 <- jd_tp5
jd_tp5_ver1$jags.data$X_pred = seq(0)

model_runjags <- autorun.jags(
  model = "../../R_test/84_jagsmodel_tp5_ver1.txt", 
  data = jd_tp5$jags.data,
  monitor = c("b","rho","lambda","scale"), n.chains = 4)
model_mcmc <- coda::as.mcmc(model_runjags)
q_tp5 <- summary(model_mcmc)$quantiles



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Left-censored, final version ---- 
#
# Start with n objects in memory  
# Strategy, see "working version"  
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Simulate data
set.seed(2) ## simulate some data... 
n <- 400
dat <- gamSim(1,n=n,dist="normal",scale=2)  

# we will use only x2 and y, and x2 is renamed 'x'
dat <- dat[c("x2", "y")]
names(dat)[1] <- "x"

ggplot(dat, aes(x, y)) +
  geom_point()

#
# Create censored data
#

# Here: fixed threshold (but that is not necessary)
thresh <- 4
dat_cens <- dat[c("x","y")]
dat_cens$y_orig <- dat_cens$y       # original (will not be used)
sel_uncens <- dat_cens$y > thresh
dat_cens$y[!sel_uncens] <- NA
dat_cens$cut <- thresh
dat_cens$cut[sel_uncens] <- NA
dat_cens$uncensored <- 0
dat_cens$uncensored[sel_uncens] <- 1

# Order file with uncensored data first
dat_cens_ordered1 <- bind_rows(
  dat_cens[sel_uncens,],
  dat_cens[!sel_uncens,]
)
# y_comb is the combination of y and cut
# - will be used as the response in the analysis
# - will only affect the uncensored values
dat_cens_ordered1$y_comb <- c(dat_cens[sel_uncens, "y"], 
                              dat_cens[!sel_uncens, "cut"])

ggplot() +
  geom_point(data = dat_cens[sel_uncens,], aes(x = x, y = y)) +
  geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)

#
# Model files (specific for k=5)
#

# Model file that will be used for left-censored analysis
# (made based on the 'orig' model)
jags.file_tp5_lc <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens.txt"
# Model file that will be made by 'jagam'  
jags.file_tp5_orig <-  "C:/Data/R_test/84_mgcv_jagam/84_jagsmodel_tp5_leftcens_orig.txt"

#
# Make x data that will be used for making the fittd line + conf.int.
# - to be added to the rest of the data
# - we also make a y (could perhaps be NA?)
#
# Make x
dat_for_fit <- data.frame(
  x = seq(0, 1, length = 30)
)
# Make y
# Use ordinary gam (on uncensored data) to add the y
# (Reason: perhaps the y value affects some sort of normalization)
dat_for_fit$y_comb <- predict(
  mgcv::gam(y_comb~x, data = dat_cens_ordered1),
  newdata = dat_for_fit)

#
# Add "x data for fitted line" to the "actual" data
#
# dat_cens_ordered1 = data file
# dat_cens_ordered2 = with addition for daa for fitted line (3o points)
dat_cens_ordered2 <- bind_rows(
  dat_cens_ordered1,
  dat_for_fit
)

nrow(dat_cens_ordered1)
nrow(dat_cens_ordered2)

# Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
# Make (2) jd_tp5_lc$jags.data (will be manpulated below)
jd_tp5_lc <- jagam(y_comb ~ s(x, bs="tp", k=5), 
                   data = dat_cens_ordered2,    # file no. 2 here
                   file = jags.file_tp5_orig,   # this file will be overwritten (was used as basis for "_leftcens" file)
                   sp.prior = "gamma", 
                   diagonalize = TRUE)

# Modify jags.data object
n <- sum(dat_cens_ordered1$uncensored %in% 1)   # file no. 1 here 
m <- sum(dat_cens_ordered1$uncensored %in% 0)   # file no. 1 here
jd_tp5_lc$jags.data$n <- n   # - makes sure only these data are use for the likelihood
jd_tp5_lc$jags.data$m <- m   #    - " -
jd_tp5_lc$jags.data$Z <- c(rep(0,n), rep(1, m))
jd_tp5_lc$jags.data$cut <- dat_cens_ordered1$cut[dat_cens_ordered1$uncensored %in% 0]

# Check
str(jd_tp5_lc$jags.data, 1)

#
# . run model - alt. 1 ---- 
# Using rjags::jags.model 
#
# edit(file = jags.file_tp5_lc)
jm <- jags.model(jags.file_tp5_lc, 
                 data=jd_tp5_lc$jags.data, 
                 inits=jd_tp5_lc$jags.ini, 
                 n.chains=2)  # changed from n.chains=1 
# list.samplers(jm)
update(jm, n.iter=10000)

# Sample main varables  
sam <- jags.samples(jm, c("b","rho","lambda","scale"), n.iter=10000, thin=10)

# Sample varaibles that have been inserted only to get the fitted line
mu_fitted_names <- paste0("mu[", nrow(dat_cens_ordered1)+1, ":", nrow(dat_cens_ordered2), "]")

sam_fit <- jags.samples(jm, mu_fitted_names, n.iter=2000, thin=2)

# Make 'plot_data'  
fitted_n <- dim(sam_fit[[1]])[1]
qs <- c(0.025, 0.5, 0.975)
plot_data <- seq_len(fitted_n) %>% map_dfr(~quantile(sam_fit[[1]][.x,,], qs))
names(plot_data) <- paste0("Q_", qs)
plot_data <- data.frame(x = dat_for_fit$x, plot_data)

ggplot(plot_data, aes(x))+
  geom_ribbon(aes(ymin = Q_0.025, ymax = Q_0.975), fill = "lightblue") +
  geom_line(aes(y = Q_0.5)) +
  geom_point(data = dat_cens_ordered1, aes(x, y)) +
  geom_point(data = dat_cens_ordered1, aes(x, cut), shape = 6)
# Worked :-)

#
# . run model - alt. 2 ---- 
# Using runjags::autorun.jags
#

# Choose the parameters to watch
model_parameters_for_convergence <- c("b","rho","lambda","scale")
# Sample varaibles that have been inserted only to get the fitted line
mu_fitted_names1 <- paste0("mu[", nrow(dat_cens_ordered1)+1, ":", nrow(dat_cens_ordered2), "]")
mu_fitted_names2 <- paste0("mu[", seq(nrow(dat_cens_ordered1)+1, nrow(dat_cens_ordered2)), "]")

n.burnin <- 4000
n.thin <- 2
n.iter <- 4000

### Run model
# Initial run  
model_converged <- runjags::autorun.jags(
  data = jd_tp5_lc$jags.data,
  monitor = model_parameters_for_convergence,     
  inits = jd_tp5_lc$jags.ini,
  model = jags.file_tp5_lc,
  n.chains = 2,    # Number of different starting positions
  startsample = 4000,     # Number of iterations
  startburnin = n.burnin, # Number of iterations to remove at start
  thin = n.thin)          # Amount of thinning

# Add all model parameters and get samples for them
model_result <- runjags::extend.jags(model_converged, 
                                     add.monitor = mu_fitted_names1,
                                     sample = n.iter)

# model_result
model_mcmc <- coda::as.mcmc(model_result)

summary <- summary(model_mcmc)
quants <- summary$quantiles
pick_rownames <- rownames(quants) %in% mu_fitted_names2
# Make 'plot_data'  
# y and lower and upper CI  values are back-transformed (un-normalized) using unnorm:
plot_data <- data.frame(
  x = dat_for_fit$x, 
  y = quants[pick_rownames,"50%"],
  y_lo = quants[pick_rownames,"2.5%"],
  y_hi = quants[pick_rownames,"97.5%"]  )

ggplot(plot_data, aes(x))+
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi), fill = "lightblue") +
  geom_line(aes(y = y)) +
  geom_point(data = dat_cens_ordered1, aes(x, y)) +
  geom_point(data = dat_cens_ordered1, aes(x, cut), shape = 6)
# Worked :-)


get_dic(jm)

effectiveSize(as.mcmc.list(sam_fit[[mu_fitted_names1]]))


# 
# Appendix ----
#

effectiveSize(as.mcmc.list(sam$b))
