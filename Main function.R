
#' Title: Main function for Dir - SPGLM
#'
#' @param X design matrix
#' @param y response variable
#' @param spt support of y 
#' @param init MCMC initialization for beta 
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#' @param iter no of MCMC iterations
#'
#' @return
#' @export
#'
#' @examples
dir_spglm <- function(X, y, spt, init, rho, iter){
  X <- X %>% as.matrix()
  y <- y
  
  n <- length(y)
  p <- dim(X)[2]
  l <- length(spt)
  max_spt <- spt[l]
  mu0 <- mean(y)
  
  # Initialization 
  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  f0_samples <- matrix(NA, nrow = iter, ncol = l)
  
  beta_samples[1, ] <- init$beta
  
  f0 <- rep(1/l, l)
  tht0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL,
                           ySptIndex = NULL)$theta
  f0star <- (f0 * exp(tht0 * spt)) %>% `/` (sum(.))
  f0_samples[1, ] <- f0star
  
  dir_pr_parm <- dtpois(spt, lambda = mu0, b = max_spt)
  ind_mt <- outer(y, spt, `==`) * 1
  
  beta <- beta_samples[1, ]
  f0 <- f0_samples[1, ]
  mu <- (X %*% beta) %>% exp()
  out <- tht_sol(spt, f0, mu, NULL)
  tht <- out$tht
  btht <- out$btht
  bpr2 <- out$bpr2
  f0_y <- f0y(y, spt, f0)
  
  # MH 
  for(j in 2:iter){
    
    # beta update
    Sig <- Sigma_beta(X, mu, bpr2, rho)
    out <- beta_update_separate(X, y, spt, beta, Sig, f0, tht, bpr2, btht, rho)
    beta <- out$cr_bt
    tht <- out$cr_tht
    btht <- out$cr_btht
    bpr2 <- out$cr_bpr2
    mu <- (X %*% beta) %>% exp()
    
    
    
    
    # f0 update
    propsl_dir_parm <- dir_parm(y, tht, btht, dir_pr_parm, ind_mt)
    # 1. without tilting to achieve mean mu0 (each time)
    out <- f0_update(y, spt, f0, f0_y, propsl_dir_parm, mu, tht, bpr2, btht,
                     dir_pr_parm, ind_mt)
    f0 <- out$cr_f0
    f0_y <- out$cr_f0y
    tht <- out$cr_tht
    btht <- out$cr_btht
    bpr2 <- out$cr_bpr2
    
    # # 2. tilting step to achieve mean mu0 (each time)
    # out <- f0_update(y, spt, f0, f0y, propsl_dir_parm, mu, tht, bpr2, btht,
    #                  dir_pr_parm, ind_mt)
    # f0 <- out$cr_f0
    # tht0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL,
    #                          ySptIndex = NULL)$theta
    # f0star <- (f0 * exp(tht0 * spt)) %>% `/` (sum(.))
    # f0 <- f0star
    # out <- tht_sol(spt, f0, mu, tht)
    # tht <- out$tht
    # btht <- out$btht
    # bpr2 <- out$bpr2
    
    
    # storing beta & f0
    beta_samples[j, ] <- beta
    f0_samples[j, ] <- f0
  }
  return(list(beta_samples = beta_samples, f0_samples = f0_samples))
}


