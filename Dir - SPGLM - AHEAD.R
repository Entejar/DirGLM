# Clear the workspace
rm(list = ls())

# Load required packages
library(foreach)
library(doParallel)

library(readstata13)
library(mvtnorm)
library(gldrm)
library(tidyverse)

# Parallelization setup
num_cores <- 2

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load functions for Dir - SPGLM
source("MCMC functions.R")
source("Main function.R")

dat <- read.dta13("~/Documents/Projects/Dir - SPGLM - R/Main/ah1993.dta")
dat$sex <- dat$sex-1
dat1 <- model.matrix(~numiadl + age + factor(sex) + iwr + factor(netwc), data = dat)
X <- dat1[, -2]
y <- dat1[, 2]
n <- length(y)

# case 1: small sample size (n = 100) and case 2: full data (n = 6441)
set.seed(500)
s_trn_ind <- sample(1:n, 100)
f_trn_ind <- 1:n
trn_ls <- list(s_trn_ind, f_trn_ind)
data <- list(X[trn_ls[[1]], ], y[trn_ls[[1]]], X[trn_ls[[2]], ], y[trn_ls[[2]]])

dataa <- model.matrix(~ numiadl + age + sex + iwr + netwc, data = dat)
iter_v <- c(10000, 10000)

ahead_study <- list()
trn_ind <- list()

ahead_study <- foreach(
  index = 1:2,
  .packages = c("gldrm", "extraDistr", "mvtnorm", "tidyr")  
) %dopar% {
  X <- as.matrix(data[[(2 * index) - 1]])
  X[, c(2,4)] <- scale(X[, c(2,4)])
  y <- as.matrix(data[[2 * index]])
  mu0 <- mean(y)
  
  trn_ind <- trn_ls[[index]]
  dat_spglm <- dataa[trn_ind, ]
  dat_spglm[, c(3, 5)] <- scale(dat_spglm[, c(3, 5)])
  dat_spglm <- dat_spglm %>% as.data.frame()
  dat_spglm$sex <- factor(dat_spglm$sex)
  dat_spglm$netwc <- factor(dat_spglm$netwc)
  fit <- gldrm(numiadl ~ age + sex + iwr + netwc, link = "log", mu0 = mu0, 
               data = dat_spglm)
  
  iter <- iter_v[index]
  mcmc_samples <- dir_spglm(X, y, spt = 0:5, init = fit, rho = 1, iter = iter)
  
  out <- list(trn_ind = trn_ind,
              mcmc_samples = mcmc_samples
              )
  return(out)
}

# Stop the parallel backend
stopCluster(cl)

# Save the results
saveRDS(ahead_study, file = "dir-spglm_ahead_study.rds")


