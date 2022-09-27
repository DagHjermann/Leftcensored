

jagscode_txt <- '
model {
  mu <- X %*% b ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:(k-1))) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(k)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
' 
