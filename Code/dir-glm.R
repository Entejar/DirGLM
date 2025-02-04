## Note that: The functions below are written w.r.t. the log link as in our data
## applications.For using a different link function for your analysis, a very few
## places need modification --- check for `# due to log link ...` comment.

#' Title: Main function for Dir-GLM
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
#'

dir_glm <- function(X, y, spt, init, rho, burnin, thin, save) {
  X <- X %>% as.matrix()
  y <- y

  n <- length(y)
  p <- dim(X)[2]
  l <- length(spt)
  max_spt <- spt[l]
  mu0 <- mean(y)

  ## Initialization ----
  iter <- burnin + thin * save
  beta_samples <- matrix(NA, nrow = save, ncol = p)
  f0_samples <- matrix(NA, nrow = save, ncol = l)

  beta_samples[1, ] <- init$beta

  f0 <- rep(1 / l, l)
  tht0 <- gldrm:::getTheta(
    spt = spt,
    f0 = f0,
    mu = mu0,
    sampprobs = NULL,
    ySptIndex = NULL
  )$theta
  f0star <- (f0 * exp(tht0 * spt)) %>% `/` (sum(.))
  f0_samples[1, ] <- f0star

  ind_mt <- outer(y, spt, `==`) * 1
  alpha <- 1
  dir_pr_parm <- alpha * colMeans(ind_mt)
  eps <- 1e-6
  dir_pr_parm <- dir_pr_parm + eps

  beta <- beta_samples[1, ]
  f0 <- f0_samples[1, ]
  mu <- (X %*% beta) %>% exp()    # due to log link: g^{-1}(y) = exp(y)
  out <- tht_sol(spt, f0, mu, NULL)
  tht <- out$tht
  btht <- out$btht
  bpr2 <- out$bpr2
  f0_y <- f0y(y, spt, f0)

  ## MH loop ----
  for (r in 2:iter) {
    ## beta update ----
    Sig <- Sigma_beta(X, mu, bpr2, rho)
    out <- beta_update_joint(X, y, spt, beta, Sig, f0, tht, bpr2, btht, rho) # for one-by-one updating, change to beta_update_separate()
    beta <- out$cr_bt
    tht <- out$cr_tht
    btht <- out$cr_btht
    bpr2 <- out$cr_bpr2
    mu <- (X %*% beta) %>% exp()     # due to log link: g^{-1}(y) = exp(y)

    ## f0 update ----
    propsl_dir_parm <- dir_parm(y, tht, btht, dir_pr_parm, ind_mt)
    out <- f0_update(y,
                     spt,
                     f0,
                     f0_y,
                     propsl_dir_parm,
                     mu,
                     tht,
                     bpr2,
                     btht,
                     dir_pr_parm,
                     ind_mt)
    f0 <- out$cr_f0
    f0_y <- out$cr_f0y
    tht <- out$cr_tht
    btht <- out$cr_btht
    bpr2 <- out$cr_bpr2


    ## storing beta & f0 ----
    if (r > burnin & r %% thin == 0) {
      j <- (r - burnin) / thin
      beta_samples[j, ] <- beta
      f0_samples[j, ] <- f0
    }
  }
  return(list(beta_samples = beta_samples, f0_samples = f0_samples))
}


