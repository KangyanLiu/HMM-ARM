ArHardAssign <- function(dat, nphi, nvar,
                           control=list(epsilon=1e-3, maxit=200),
                           screen=TRUE){
  # Initialization
  IJ <- 3 * nvar
  s <- 3 * nvar * nphi
  # Row 1 for mu, row 2 for s2r, row 3 for phi.
  comb <- t(expand.grid(1:3, 1:nvar, 1:nphi)) 
  m1 <- mean(dat[which(dat > 0)])
  m2 <- mean(dat[which(dat < 0)])
  s2r <- c(0.01, 0.05, 0.1, 0.15)[1:nvar]
  para <- c(runif(nphi, 0.4, 0.6), m1, m2, s2r)
  N <- length(dat)
  n <-  N - 1
  x <- dat[-N]
  y <- dat[-1]
  Y <- rep(y, s)
  vec <- rep(comb[3, ], each=n) - 1
  Xphi <- table(1:length(vec), as.factor(vec))
  vec <- rep(comb[1, ], each=n) - 1
  Xmu <- table(1:length(vec), as.factor(vec))[,2:3]
  X <- cbind(Xphi * x, Xmu)
  loglikelihood.old <- -1.797693e+308
  for (it in 1:control$maxit) {
    ######### E step ##########
    mu <- c(0, para[1 + nphi], para[2 + nphi])
    curtmu <- matrix(rep(mu[comb[1, ]], each=n), ncol=s)
    std <- sqrt(para[(3 + nphi):(2 + nphi + nvar)])
    curtstd <- matrix(rep(std[comb[2, ]], each=n), ncol=s)
    phi <- para[1:nphi]
    curtphi <- matrix(rep(phi[comb[3, ]], each=n), ncol=s)
    z <- matrix(rep(y, s) - curtphi * rep(x, s), nrow=n)
    lemission <- t(dnorm(z, curtmu, curtstd, log=TRUE))
    loglikelihood <- sum(apply(lemission, 2, max))
    ######### M step ##########
    states <- apply(lemission, 2, which.max)
    mu.states <- comb[1, states]
    s2r.states <- comb[2, states]
    phi.states <- comb[3, states]
    xphi <- table(1:length(phi.states), as.factor(phi.states))
    xmu <- table(1:length(mu.states), as.factor(mu.states))[, 2:3]
    xx <- cbind(xphi * x, xmu)
    s2r <- para[(3 + nphi):(2 + nphi + nvar)][s2r.states] 
    W <- diag(s2r)
    b <- solve(t(xx) %*% W %*% xx) %*% t(xx) %*% W %*% y
    para[1:(2 + nphi)] <- b
    se <- (y - xx %*% b)^2
    for (i in 1:nvar){
      sei <- se[s2r.states == i]
      para[2 + nphi + i] <- sum(sei) / length(sei)
    }
    if (screen) {
      cat(para, loglikelihood, "\n")
    }
    if (loglikelihood - loglikelihood.old > control$epsilon) {
      loglikelihood.old <- loglikelihood
    }else{
      return(list(parhat=para, logll=loglikelihood,
                  dat=dat, states=states, it=it))
    }
  }
}
