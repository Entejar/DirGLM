# Clear environment
rm(list = ls())

# Required packages
require(foreach)
require(doParallel)

require(mvtnorm)
require(gldrm)
require(extraDistr)
require(tidyr)

# Parallelization setup
num_cores <- 2

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load functions for Dir - SPGLM
source("MCMC functions.R")
source("Main function.R")


# Simulation design: DGM 1
sim_dat_1 <- function(n, p){
  X <- matrix(c(rep(1, n), rnorm(n, mean = 0, sd = 1)), byrow = F,
              nrow = n, ncol = p)
  beta <- c(- 0.7, 0.2)     ## true beta
  mu <- exp(X %*% beta)     ## true mu, g(mu) = log(mu) = X*beta
  spt <- 0:5                ## support of y
  f0 <- dtpois(spt, lambda = 1, b = 5)     ## true f0
  theta <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu, sampprobs = NULL,
                            ySptIndex = NULL)$theta    ## true theta
  prob <- t(sapply(1:n, function(i) exp(theta[i]*spt)*f0))
  y <- sapply(1:n, function(i) sample(spt, 1, prob = prob[i, ])) ## simulated y
  sim.data <- data.frame(X, y)
  return(sim.data)
}


# Simulation design: DGM 2
sim_dat_2 <- function(n, p){
  X <- matrix(c(rep(1, n), rnorm(n, mean = 0, sd = 1)), byrow = F,
              nrow = n, ncol = p)
  beta <- c(- 0.7, 0.2)     ## true beta
  mu <- exp(X %*% beta)     ## true mu, g(mu) = log(mu) = X*beta
  spt <- 0:5                ## support of y
  temp <- dpois(spt, lambda = 1)
  f0 <- numeric(6)
  f0[1] <- 3*temp[1]/(3*temp[1] + sum(temp[-1]))
  f0[-1] <- temp[-1]/(3*temp[1] + sum(temp[-1]))
  theta <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu, sampprobs = NULL,
                            ySptIndex = NULL)$theta    ## true theta
  prob <- t(sapply(1:n, function(i) exp(theta[i]*spt)*f0))
  y <- sapply(1:n, function(i) sample(spt, 1, prob = prob[i, ])) ## simulated y
  sim.data <- data.frame(X, y)
  return(sim.data)
}


# Simulation study
num_datasets <- 1000
sim_study <- list()
settings <- c(1, 2)
seed_init <-  sample.int(.Machine$integer.max, 1)

sim_study <- foreach(setting = settings,
                     .packages = c("gldrm", "extraDistr", "mvtnorm", "tidyr")) %:%
  foreach(dat_index = 1:num_datasets) %dopar% {
    seed_val <- seed_init + dat_index
    set.seed(seed_val)
    n <- 250
    dat <- if (setting %in% c(1, 2)) sim_dat_1(n = n, p = 2) else sim_dat_2(n = n, p = 2)
    X <- dat[, -3]
    y <- dat[, 3]
    fit <- gldrm(y ~ X2, link = "log", mu0 = mean(y), data = dat)
    mcmc_samples <- dir_spglm(X, y, spt = 0:5, init = fit, rho = 1, iter = 5000)
    out <- list(seed = seed_val,
                set = setting,
                data = dat,
                mcmc_samples = mcmc_samples
                )
    return(out)
  }

# Stop parallel backend
stopCluster(cl)

# Save result
saveRDS(sim_study, file = "dir-spglm_sim_study.rds")