rm(list = ls())

#...............................Packages........................................
require(mvtnorm)
require(latex2exp)
require(gldrm)
require(extraDistr)
require(tidyverse)
require(ggpubr)
#.............................Data Generation...................................
## trying this -- Peter
set.seed(442)
n <- 250; p <- 2
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
data <- data.frame("x1" = X[, 2], y)
data$x1 <- scale(data$x1)

#.........................MH Algorithm..........................................
itr <- 2000  
unique.y <- sort(unique(y))
freq <- as.vector(table(y))
K <- length(unique.y)
S <- 100    
H <- S*K
M <- 1
mu0 <- mean(y)

### Initialization
beta.mat <- matrix(0, nrow = itr, ncol = 2)
beta.mat[1, ] <- coef(glm(formula = y ~ x1, family = poisson, data = data))

f0.mat <- matrix(list(), nrow = itr, ncol = 2)
f0.mat[1, ] <- list(rep(unique.y, S), rep(rep(1/K, K), S))
f0.mat[1, ][[2]] <- f0.mat[1, ][[2]]/sum(f0.mat[1, ][[2]])
mh.tilde <- rep(0, H)

## indicator and indicator.tilde: mixture component indicator (for m_h) 
indicator <- rep(0, H)   
indicator.tilde <- rep(0, H) 

for(j in 2:itr){
  ### solve theta [based on x, and current beta & f_0]
  mu <- exp(X %*% beta.mat[j-1, ])
  mh <- f0.mat[j-1, ][[1]]
  wh <- f0.mat[j-1, ][[2]]
  A <- gldrm:::getTheta(spt = mh, f0 = wh, mu = mu, sampprobs = NULL, 
                        ySptIndex = NULL)
  theta <- A$theta
  
  ### calculate sigma matrix for generating beta (tilde) [proposal]
  g.prime <- 1/A$bPrime
  bprime2 <- A$bPrime2  
  g2b <- ((g.prime)^2)*bprime2
  Xtemp <- X/sqrt(g2b)
  info_mat <- t(Xtemp) %*% Xtemp
  info_mat.inv <- solve(info_mat)
  
  ### generate beta (tilde) [proposal]
  beta.tilde <- rmvnorm(1, mean = beta.mat[j-1, ], sigma = info_mat.inv)
  
  ### generate m_h (tilde) and indicator (tilde) [proposal]
  btheta <- t(sapply(1:n, function(i) log(sum(exp(theta[i]*mh)*wh))))
  c.prob <- 1/exp(theta*y - btheta)
  a.prob <- ((1/n)*c.prob)/sum(c.prob)
  result <- sapply(1:H, function(h) {
    if (runif(1) < M/(M+n)) {
      mh.tilde[h] <<- rtpois(n = 1, lambda = mu0, b = 5)
      indicator.tilde[h] <<- 1
    } else {
      mh.tilde[h] <<- sample(y, 1, prob = a.prob)
      indicator.tilde[h] <<- 0
    }
  }
  )
  
  ### generate omega_h (tilde) [proposal]
  temp <- rbeta(H-1, 1, M+n)
  wh.tilde <- c(temp[1], c(temp[-1], 1)*cumprod(1 - temp))
  
  ### calculate theta (tilde) [while loop to ensure existence of theta (tilde)]
  mu.tilde <- exp(X %*% as.numeric(beta.tilde))
  while(sum(mu.tilde < min(mh.tilde)) >= 1 || sum(mu.tilde > max(mh.tilde)) >= 1
  ){
    result <- sapply(1:H, function(h) {
      if (runif(1) < M/(M+n)) {
        mh.tilde[h] <<- rtpois(n = 1, lambda = mu0, b = 5)
        indicator.tilde[h] <<- 1
      } else {
        mh.tilde[h] <<- sample(y, 1, prob = a.prob)
      }
    }
    )
    beta.tilde <- rmvnorm(1, mean = beta.mat[j-1, ], sigma = info_mat.inv) 
    mu.tilde <- exp(X %*% as.numeric(beta.tilde))
  }
  
  Atilde <- gldrm:::getTheta(spt = mh.tilde, f0 = wh.tilde, mu = mu.tilde, 
                             sampprobs = NULL, ySptIndex = NULL)
  thetatilde <- Atilde$theta
  btildetheta <- t(sapply(1:n, function(i) log(sum(exp(thetatilde[i]*mh.tilde)*
                                                     wh.tilde))))
  c.probtilde <- 1/exp(thetatilde*y - btildetheta)
  a.probtilde <- ((1/n)*c.probtilde)/sum(c.probtilde)
  
  ### calculate sigma matrix for q(beta) [based on beta (tilde) & f_0 (tilde)] 
  gtilde.prime <- 1/Atilde$bPrime
  btildeprime2 <- Atilde$bPrime2  
  g2btilde <- ((gtilde.prime)^2)*btildeprime2
  Xtildetemp <- X/sqrt(g2btilde)
  info_mattilde <- t(Xtildetemp) %*% Xtildetemp
  info_mat.invtilde <- solve(info_mattilde)
  
  ### calculate f_0(y) and q(m_h)  
  f0.y <- f0tilde.y <- numeric(n)
  q.mh <- qtilde.mh <- numeric(H)
  indices <- expand.grid(i = 1:n, h = 1:H)
  result <- sapply(1:nrow(indices), function(k) {
    i <- indices$i[k]
    h <- indices$h[k]
    if (y[i] == mh.tilde[h]) {
      f0tilde.y[i] <<- f0tilde.y[i] + wh.tilde[h]
      qtilde.mh[h] <<- qtilde.mh[h] + a.prob[i]
    }
    if (y[i] == mh[h]) {
      f0.y[i] <<- f0.y[i] + wh[h]
      q.mh[h] <<- q.mh[h] + a.probtilde[i]
    }
  }
  )
  
  ### calculate acceptance probability [alpha, in log scale]
  
  if(sum(qtilde.mh == 0) == 0 && sum(q.mh == 0) == 0){
    llik.tilde <- sum(thetatilde*y - btildetheta + log(f0tilde.y))
    llik <- sum(theta*y - btheta + log(f0.y))
    pbeta.tilde <- dmvnorm(beta.tilde, log = T) 
    pbeta <- dmvnorm(beta.mat[j-1, ], log = T)
    pmh.tilde <- sum(dtpois(mh.tilde, lambda = mu0, b = 5, log = T))  
    pmh <- sum(dtpois(mh, lambda = mu0, b = 5, log = T))
    pwh.tilde <- sum(dbeta(wh.tilde[-H], shape1 = 1, shape2 = M, log = T))
    pwh <- sum(dbeta(wh[-H], shape1 = 1, shape2 = M, log = T))  
    qbeta <- dmvnorm(beta.mat[j-1, ], mean = beta.tilde, 
                     sigma = info_mat.invtilde, log = T)  
    qbeta.tilde <- dmvnorm(beta.tilde, mean = beta.mat[j-1, ], sigma = 
                             info_mat.inv, log = T)
    qwh <- sum(dbeta(wh[-H], shape1 = M, shape2 = M+n, log = T)) 
    qwh.tilde <- sum(dbeta(wh.tilde[-H], shape1 = M, shape2 = M+n, log = T))  
    qmh <- sum((log(M/(M+n))*indicator + 
                  dtpois(mh, lambda = mu0, b = 5, log = T))*indicator) +  
      sum((log(n/(M+n))*(1 - indicator)) + log(q.mh)*(1 - indicator))
    qmh.tilde <- sum((log(M/(M+n))*indicator.tilde + 
                        dtpois(mh.tilde, lambda = mu0, b = 5, log = T))*
                       indicator.tilde) +
      sum((log(n/(M+n))*(1 - indicator.tilde)) + log(qtilde.mh)*
            (1 - indicator.tilde)) 
    alpha <- min(0, (llik.tilde - llik + pbeta.tilde - pbeta + pmh.tilde - pmh + 
                       pwh.tilde - pwh + qbeta - qbeta.tilde + qwh - qwh.tilde + 
                       qmh - qmh.tilde))
  } else if(sum(q.mh == 0) >= 1 && sum(qtilde.mh == 0) == 0){
    llik <- sum(theta*y - btheta + log(f0.y))
    llik.tilde <- sum(thetatilde*y - btildetheta + log(f0tilde.y))
    alpha <- ifelse(llik - llik.tilde < 0, 0, -Inf)
    ## 0/(non-zero) case in acceptance probability
  } else if(sum(q.mh == 0) >= 1 && sum(qtilde.mh == 0) >= 1){
    llik <- sum(theta*y - btheta + log(f0.y))
    llik.tilde <- sum(thetatilde*y - btildetheta + log(f0tilde.y))
    alpha <- ifelse(llik - llik.tilde < 0, 0, -Inf)
    ## 0/0 case in acceptance probability. accept if log-likelihood of proposed 
    ## is greater than that of current
  } else{
    llik <- sum(theta*y - btheta + log(f0.y))
    llik.tilde <- sum(thetatilde*y - btildetheta + log(f0tilde.y))
    alpha <- ifelse(llik - llik.tilde < 0, 0, -Inf)
    ## (non-zero)/0 case in acceptance probability
  }

  ### MH update step
  
  if(log(runif(1)) < alpha){
    beta.mat[j, ] <- beta.tilde
    f0.mat[j, ] <- list(mh.tilde, wh.tilde)
    indicator <- indicator.tilde
    } else{
    beta.mat[j, ] <- beta.mat[j-1, ]
    f0.mat[j, ] <- f0.mat[j-1, ]
    }
}


#.................................Results.......................................
### tilt to achieve pre-specified mean of f_0 [simulation truth = 1]
f0mat.tilt <- matrix(list(), nrow = itr, ncol = 2)
for(r in 1:itr){
  mh <- f0.mat[r, ][[1]]
  wh <- f0.mat[r, ][[2]]
  theta0 <- gldrm:::getTheta(spt = mh, f0 = wh, mu = 1, 
                             sampprobs = NULL, ySptIndex = NULL)$theta
  wh <- wh*exp(theta0*mh)
  wh <- wh/sum(wh)
  f0mat.tilt[r, ] <- list(mh, wh)
}

est_f0 <- matrix(0, nrow = itr, ncol = 6)

for(j in 1:itr){
  mh <- f0mat.tilt[j,][[1]]
  wh <- f0mat.tilt[j,][[2]]/sum(f0mat.tilt[j,][[2]])
  for(k in 1:6){
    for(h in 1:H){
      if(spt[k] == mh[h]){est_f0[j, k] <- est_f0[j, k] + wh[h]}
    }
  }
}
colMeans(est_f0) ## f_0 estimates
f0               ## f_0 true values with lambda = 1

### Trace plots of beta
index <- 1:itr
df <- cbind(index, as.data.frame(beta.mat))

p1 <- ggplot(df, aes(x = index, y = V1)) +
  geom_line(col = "blue") + labs(title = "Trace plots", y = TeX("$\\beta_0$")) +
  geom_hline(aes(yintercept = -0.7), col = "red") +
  theme_bw() + theme(axis.title.x = element_blank()) + ylim(-1, -0.4)

p2 <- ggplot(df, aes(x = index, y = V2)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\beta_1$")) + 
  geom_hline(aes(yintercept = 0.2), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) + ylim(0, 0.6)

ggarrange(p1, p2, ncol = 1, nrow = 2)

### Trace plots of f_0
index <- 1:itr
df <- cbind(index, as.data.frame(est_f0))

p1 <- ggplot(df, aes(x = index, y = V1)) +
  geom_line(col = "blue") + labs(title = "Trace plots", y = TeX("$\\f_0(0)$")) +
  geom_hline(aes(yintercept = 0.368), col = "red") +
  theme_bw() + theme(axis.title.x = element_blank()) 

p2 <- ggplot(df, aes(x = index, y = V2)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\f_0(1)$")) + 
  geom_hline(aes(yintercept = 0.368), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) 

p3 <- ggplot(df, aes(x = index, y = V3)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\f_0(2)$")) + 
  geom_hline(aes(yintercept = 0.184), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) 

p4 <- ggplot(df, aes(x = index, y = V4)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\f_0(3)$")) + 
  geom_hline(aes(yintercept = 0.061), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) 

p5 <- ggplot(df, aes(x = index, y = V5)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\f_0(4)$")) + 
  geom_hline(aes(yintercept = 0.015), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) 

p6 <- ggplot(df, aes(x = index, y = V6)) + 
  geom_line(col = "blue") + labs(y = TeX("$\\f_0(5)$")) + 
  geom_hline(aes(yintercept = 0.003), col = "red") + theme_bw() +
  theme(axis.title.x = element_blank()) 

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
