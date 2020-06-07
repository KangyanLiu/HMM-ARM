#library(Rcpp)
#sourceCpp("Hmm.cpp")

ECMinit <-function(dat, nphi, nvar) {
  m1 <- mean(dat[which(dat > 0)]) * 0.5
  m2 <- mean(dat[which(dat < 0)]) * 0.5
  s <- 3 * nvar * nphi  # number of states
  A <- matrix(1 / s, s, s) # Uniform initial probability
  # s2r <- rep(0.01,nvar)+abs(rnorm(nvar,0,0.01))
  s2r <- c(0.01, 0.05, 0.1, 0.15)[1:nvar]
  para <- c(runif(nphi, 0.4, 0.6), m1, m2, s2r)
  Pi <- solve(t(diag(s) - A + 1), rep(1, s))
  return(list(A=A, Pi=Pi, para=para))
}
# Only a part of the likelihood depends on phi. We just need to maxmize this part instead of the whole likelihood.
# This part can be vectorized.
llphi <- function(phi, y, x, D, curtmu, curtstd, p, IJ) {
  n <- length(y)
  curtmu <- curtmu * (1 - phi^D) / (1 - phi)
  curtvar <- curtstd^2 * (1 - phi^(2 * D)) / (1 - phi^2)
  ll <-sum((-log(curtvar) / 2 - (rep(y, IJ) -
           phi^D * rep(x, IJ) - curtmu)^2 / (2 * curtvar)) * p)
  return(-ll)
}

ECM <- function(dat, init=NULL, nphi, nvar,
                control=list(epsilon=1e-6, maxit=200),
                screen=TRUE, underflow='Scale') {
  ind <- which(!is.na(dat))
  D <- ind[-1] - ind[-length(ind)]
  dat <- dat[!is.na(dat)]
  if (is.null(init)){
    init <- ECMinit(dat, nphi, nvar)
  }
  IJ <- 3 * nvar
  s <- 3 * nvar * nphi
  comb <- t(expand.grid(1:3, 1:nvar, 1:nphi))
  para <- init$par
  ind <- which(!is.na(dat))
  N <- length(dat)
  n <-  N - 1
  x <- dat[-N]
  y <- dat[-1]
  Pi <- init$Pi
  A <- init$A
  Y <- rep(y, s)
  vec <- rep(comb[3, ], each=n) - 1
  Xphi <- table(1:length(vec), as.factor(vec))
  vec <- rep(comb[1,], each=n) - 1
  Xmu <- table(1:length(vec), as.factor(vec))[, 2:3]
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
    curtmu <- curtmu * (1 - curtphi^D) / (1 - curtphi)
    curtstd <- curtstd * sqrt((1 - curtphi^(2 * D)) / (1 - curtphi^2))
    z <- matrix(rep(y, s) - curtphi^D * rep(x, s), nrow=n)
    if(underflow == 'Scale'){
      emission <- t(dnorm(z, curtmu, curtstd))
      abg <- AlphaBetaGamma_scaled(emission, A, Pi)
      Gamma <- abg$Gamma
      Pi.new <- Gamma[, 1]
      A.new <- abg$A
      loglikelihood <- -sum(log(abg$Scale))
    }else if (underflow == 'logsum'){
      lemission <- t(dnorm(z, curtmu, curtstd, log=TRUE))
      labg <- lAlphaBetaGamma(lemission, log(A), log(Pi))
      Gamma <- exp(labg$lGamma)
      Pi.new <- Gamma[, 1]
      A.new <- exp(labg$lA)
      loglikelihood <- labg$logll
    }else{
      stop("Wrong method for underflow control")
    }
    ######### M step ##########
    p <- as.numeric(t(Gamma))
    s2r <- rep(para[(3 + nphi):(2 + nphi + nvar)][comb[2, ]], each=n)
    W <- 1 / (s2r / p)
    #Solve directly?
    Xc <- as.numeric((1 - curtphi^D) / (1 - curtphi)) * Xmu
    zz <- as.numeric(z)
    b <- solve(t(Xc * W) %*% Xc, t(Xc * W) %*% zz)
    para[(1 + nphi):(2 + nphi)] <- b
    squaredloss <- (zz - Xc %*% b)^2
    et <- as.numeric((1 - curtphi^(2 * D)) / (1 - curtphi^2))
    for (i in 1:nvar) {
      ind <- which(rep(comb[2, ], each=n) == i)
      para[2 + nphi + i] <- sum(squaredloss[ind] * p[ind] / et[ind]) / sum(p[ind])
    }
    mu <- c(0, para[1 + nphi], para[2 + nphi])
    curtmu <- matrix(rep(mu[comb[1, ]], each=n), ncol=s)
    std <- sqrt(para[(3 + nphi):(2 + nphi + nvar)])
    curtstd <- matrix(rep(std[comb[2, ]], each=n), ncol=s)
    for (i in 1:nphi){
      ind <- which(rep(comb[3, ], each=n) == i)
      op <- optim(par=para[i], fn=llphi, y=y, x=x, D=D, curtmu=curtmu[ind],
                  curtstd=curtstd[ind], p=p[ind], IJ=IJ,
                  method="Brent",lower=0,upper=1)
      para[i] <- op$par
    }
    if (screen) {
      cat(para, loglikelihood, "\n")
    }
    if (abs(loglikelihood - loglikelihood.old) > control$epsilon){
      loglikelihood.old <- loglikelihood
      A <- A.new
      Pi <- Pi.new
    }else{
      states <- apply(Gamma, 2, which.max)
      statehat <- Viterbi(A, Pi, t(dnorm(z, curtmu, curtstd, log=TRUE)))
      return(list(parhat=para, Trans=A.new, Initprob=Pi.new,
                  statehat=statehat + 1, logll=loglikelihood,
                  dat=dat, states=states, it=it))
    }
  }
  states <- apply(Gamma, 2, which.max)
  statehat <- Viterbi(A, Pi, t(dnorm(z, curtmu , curtstd, log=TRUE)))
  return(list(parhat=para, Trans=A.new, Initprob=Pi.new,
              statehat=statehat + 1, logll=loglikelihood,
              dat=dat, states=states, it=it))
}
