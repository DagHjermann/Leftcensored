#
# From 
# https://onofriandreapg.github.io/agriCensData/RandomEffects.html
#
# Just to test that pnrom in fact works
#

# install.packages("agriCensData")
# Package not on CRAN anymore, but I downloaded the data from
#   https://github.com/OnofriAndreaPG/agriCensData/blob/master/data/starchGrainU.rda


# Save BUGS description of the model as a string
modelSpec <- "
data{
for (i in 1:N) {
  zeros[i] <- 0
}}

model{
for (i in 1:N) {
  exp[i] <- mu[Group[i]] + gamma[Photo[i]]
}

for (i in 1:N1) {
  #Likelihood for left-censored
  S2[i] <- pnorm(high[i], exp[i], tau.e)
  L[i] <- S2[i]     #(Equation 3)       
  phi[i] <- -(log(L[i]))
  zeros[i] ~ dpois(phi[i])    
}

for (i in (N1+1):N2) {
  #Likelihood for interval-censored
  S[i] <- pnorm(low[i], exp[i], tau.e)
  S2[i] <- pnorm(high[i], exp[i], tau.e)
  L[i] <- S2[i] - S[i]  #(Equation 4)      
  phi[i] <- -(log(L[i]))
  zeros[i] ~ dpois(phi[i])    
}

for (i in (N2+1):N) {
  #Likelihood for right-censored
  S[i] <- pnorm(low[i], exp[i], tau.e)
  L[i] <- 1 - S[i]  #(Equation 5)      
  phi[i] <- -(log(L[i]))
  zeros[i] ~ dpois(phi[i])    
}

#Priors
  sigma.e ~ dunif(0, 100)
  sigma.P ~ dunif(0, 100)
for(i in 1:2){
  mu[i] ~ dnorm(0, 0.000001)
} for(i in 1:24){
   gamma[i] ~ dnorm(0, tau.P)
}

#Derived quantities
  sigma2p <- sigma.P*sigma.P
  sigma2e <- sigma.e*sigma.e
  tau.P <- 1 / sigma2p
  tau.e <- 1 / sigma2e
  diff <- mu[1] - mu[2]
}
"


# Data

# data(starchGrainU)
load("starchGrainU.rda")

dataset_jags <- starchGrainU[order(starchGrainU$Class), ]  # sorting data set
N1 <- 639
N2 <- 2130
N <- 2441


# JAGS

library(rjags)
win.data <- list(low = dataset_jags$sizeLow, 
                 high = dataset_jags$sizeUp, 
                 N1 = N1, N2 = N2, N = N, 
                 Group = factor(dataset_jags$Group),
                 Photo = factor(dataset_jags$Photo)
                 )

init <- list(mu = c(7.3, 8.7), sigma.e = 1.7, sigma.P = 0.5)
mcmc <- jags.model(textConnection(modelSpec),
                   data = win.data, 
                   inits = init, 
                   n.chains = 4, 
                   n.adapt = 100)
params <- c("mu", "sigma.e", "sigma.P", 
            "sigma2p", "sigma2e", "diff")
res3 <- coda.samples(mcmc, params, n.iter = 10000)
burn.in <- 1000

