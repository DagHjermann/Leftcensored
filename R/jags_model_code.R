
get_jags_model_code <- function(bs = "tp",
                                k = 5,
                                type = "leftcensored"){
  
  #
  # Models without measurement error ----
  #
  
  if (bs == "tp" & k == 3 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1:2,1:2] * lambda[1]  + S1[1:2,3:4] * lambda[2]
  b[2:3] ~ dmnorm(zero[2:3],K1) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  
  } else if (bs == "tp" & k == 4 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:3)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(4)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  } else if (bs == "tp" & k == 5 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:4)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(5)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

#
# Models with measurement error ----
#

} else if (bs == "tp" & k == 3 & type == "leftcensored_measerror"){
  
  code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1:2,1:2] * lambda[1]  + S1[1:2,3:4] * lambda[2]
  b[2:3] ~ dmnorm(zero[2:3],K1) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

} else if (bs == "tp" & k == 4 & type == "leftcensored_measerror"){
  
  code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:3)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(4)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else if (bs == "tp" & k == 5 & type == "leftcensored_measerror"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:4)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(5)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else {
    
    stop("The given combination of bs = ", sQuote(bs), ", k = ", k, ", and type = ", sQuote(type),
         " has not been implemented.")
    
  }

code

}


