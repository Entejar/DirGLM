## Required Packages ----

require(mvtnorm)
require(gldrm)
require(extraDistr)
require(tidyverse)

## Loading ----
source("mcmc_helpers.R")
source("dir-glm.R")


## Simulation scenario 1 ----
sim_dat_1 <- function(n, p) {
  X <- matrix(c(rep(1, n), rnorm(n, mean = 0, sd = 1)),
              byrow = F,
              nrow = n,
              ncol = p)
  beta <- numeric(p)
  beta[1:2] <- c(-0.7, 0.2)     # true beta
  mu <- exp(X %*% beta)         # true mu, g(mu) = log(mu) = X*beta
  spt <- 0:5                    # support of y
  f0 <- dtpois(spt, lambda = 1, b = 5)     # true f0
  theta <- gldrm:::getTheta(
    spt = spt,
    f0 = f0,
    mu = mu,
    sampprobs = NULL,
    ySptIndex = NULL
  )$theta                       # true theta
  prob <- t(sapply(1:n, function(i)
    exp(theta[i] * spt) * f0))
  y <- sapply(1:n, function(i)
    sample(spt, 1, prob = prob[i, ])) # simulated y
  sim.data <- data.frame(X, y)
  return(sim.data)
}


# Simulation scenario 2 ----
sim_dat_2 <- function(n, p) {
  X <- matrix(c(rep(1, n), rnorm(n, mean = 0, sd = 1)),
              byrow = F,
              nrow = n,
              ncol = p)
  beta <- numeric(p)
  beta[1:2] <- c(-0.7, 0.2)     # true beta
  mu <- exp(X %*% beta)         # true mu, g(mu) = log(mu) = X*beta
  spt <- 0:5                    # support of y
  temp <- dpois(spt, lambda = 1)
  f0 <- numeric(6)
  f0[1] <- 3 * temp[1] / (3 * temp[1] + sum(temp[-1]))
  f0[-1] <- temp[-1] / (3 * temp[1] + sum(temp[-1]))  # true f0
  theta <- gldrm:::getTheta(
    spt = spt,
    f0 = f0,
    mu = mu,
    sampprobs = NULL,
    ySptIndex = NULL
  )$theta                        # true theta
  prob <- t(sapply(1:n, function(i)
    exp(theta[i] * spt) * f0))
  y <- sapply(1:n, function(i)
    sample(spt, 1, prob = prob[i, ])) # simulated y
  sim.data <- data.frame(X, y)
  return(sim.data)
}


## Simulation ----
## Example 1

spt <- 0:5
mu0 <- 1
data <- sim_dat_1(n = 500, p = 2)
X <- data[, -3]
y <- data[, 3]
init <- gldrm(y ~ X2,
              link = "log",
              mu0 = mean(y),
              data = data)
burnin <- 1000
thin <- 1
save <- 1000

mcmc_samples <- dir_glm(
  X,
  y,
  spt = spt,
  init = init,
  rho = 1,
  burnin = burnin,
  thin = thin,
  save = save
)

beta_samples <- mcmc_samples$beta_samples
f0_samples <- mcmc_samples$f0_samples

f0star_samples <- matrix(0, nrow = nrow(f0_samples), ncol = length(spt))
for (iter in 1:nrow(f0_samples)) {
  wh <- f0_samples[iter, ]
  theta0 <- gldrm:::getTheta(
    spt = spt,
    f0 = wh,
    mu = mu0,
    sampprobs = NULL,
    ySptIndex = NULL
  )$theta
  wh <- wh * exp(theta0 * spt)
  wh <- wh / sum(wh)
  f0star_samples[iter, ] <- wh
}

f0_samples <- f0star_samples  # projected f0 samples

colMeans(beta_samples)        # true beta is c(-0.7, 0.2)
colMeans(f0_samples)          # true f0 is truncated Poisson = (0.368, 0.368, 0.184, 0.061, 0.015, 0.003) [rounded]


