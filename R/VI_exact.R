#library(Rcpp)
# If a step can not be vectorized, a function written in Cpp is used to reduce computation time.
# UpdateTheta, UpdateMu, UpdateS2r are written in Cpp.
# ELBO consists of 3 parts. l1 and l2 computes two parts. (A.14)
#sourceCpp("VICpp_meanfield.cpp")

# Main variational inference function for ARMHMM
VI <- function(dat, A=NULL, Pi=NULL,
               control=list(epsilon=1e-6, maxit=200),
               screen=T) {
  # Initialization
  phi <- runif(1, 0.25, 0.35) # Start from a small autocorrelation
  ELBO <- -1.797693e+308 # Smallest number in R
  ind <- which(!is.na(dat)) # Indicator for non-missing position
  y <- dat[ind[-1]]
  x <- dat[ind[-length(ind)]]
  D <- ind[-1] - ind[-length(ind)] # Distance between two non-missing postion
  CD <- c(0, cumsum(D)) # Cumulative distance used to calcualte likelihood
  CD <- CD[-length(CD)]
  yt <- y - phi^D * x
  m1 <- mean(yt[which(yt > 0)])
  m2 <- mean(yt[which(yt < 0)])
  mu <- matrix(c(0, m1, m2), 3, 1)
  s2r <- var(yt) / 10
  Init.states <- rep(1, length(dat) - 1)
  Init.states[which(yt > m1)] <- 2
  Init.states[which(yt < m2)] <- 3
  theta <- t(table(1:length(Init.states), as.factor(Init.states))) # States matrix
  # If initial probability or transition matrix is not provided,
  # assign values uniformly.
  if (is.null(A)) {
    A <- matrix(rep(1 / 3, 9), 3, 3)
  }
  if (is.null(Pi)){
    Pi <- matrix(rep(1 / 3, 3), 3, 1)
  }
  logA <- log(A)
  logPi <- log(Pi)
  for (it in 1:control$maxit) {
    yt <- y - phi^D * x
    C <- (1 - phi^(2 * D)) / (1 - phi^2)
    cnt1 <- 0 # A counter for iterations of E step.
    theta.new <- theta
    ELBO.new <- -1.797693e+308
    # E Step by VI
    # Convergence is measured by ELBO without normalization term for E step.
    # We may not need to wait until convergence since the ELBO is non-decreasing.
    # However, only one step does not work. More tests for convergence are needed.
    while ((cnt1 == 0) || ((ELBO.new - ELBO>control$epsilon) &&
           (cnt1 < control$maxit))){
      ELBO <- ELBO.new
      cnt1 <- cnt1 + 1
      theta <- theta.new
      theta.new <- UpdateTheta(yt, logPi, logA, D, CD, C * s2r , phi, mu, theta)
      ELBO.new <- -l1(phi, y=y, x=x, D=D, CD=CD, s2r=s2r, mu=mu, theta=theta.new) / 2 +
                   l2(logPi, logA, matrix(unlist(theta.new), nrow=3)) -
                   sum(theta.new * log(theta.new + (theta.new == 0) * 4.940656e-324))
    }
    # When we change the theta to maximize the ELBO,
    # it's possible that new theta gives a smaller ELBO.
    # Use theta instead of theta.new. This gives the highest current ELBO.
    Pi <- theta[, 1]
    A <- matrix(NA, 3, 3)
    for (i in 1:3)  {
      for (j in 1:3) {
        A[i, j] <- sum(theta[i, -1] * theta[j, -ncol(theta)]) /
                   sum(theta[j, -ncol(theta)])
      }
    }
    Pi[which(Pi == 0)] <- 4.940656e-324
    logA <- log(A)
    logPi <- log(Pi)
    # M step
    # Iteration is not needed.
    mu <- c(0, UpdateMu(yt, D, C, CD, phi, theta))
    s2r <- UpdateS2r(yt, D, CD, C, phi, mu, theta)
    op <- optim(par=phi, fn=l1, y=y, x=x, D=D, CD=CD, s2r=s2r, mu=mu, theta=theta,
                method="Brent", lower=0, upper=1)
    phi <- op$par
    ELBO.new <- -l1(phi, y=y, x=x, D=D, CD=CD, s2r=s2r, mu=mu, theta=theta) / 2 +
                 l2(logPi, logA, matrix(unlist(theta), nrow=3)) -
                 sum(theta.new * log(theta.new + (theta.new ==0 ) * 4.940656e-324))
    if (screen) {
      cat(phi, mu, s2r, ELBO, cnt1, ELBO.new - ELBO, "\n")
    }
    if (abs(ELBO.new - ELBO) > control$epsilon) {
      ELBO <- ELBO.new
    }else{
      states <- matrix(theta, 3, length(dat) - 1)
      states <- apply(states, 2, which.max)
      return(list(phi=phi, mu=mu, s2r=s2r, states=states, A=A, it=it, ELBO=ELBO))
    }
  }
  states <- matrix(theta, 3, length(dat) - 1)
  states <- apply(states, 2, which.max)
  return(list(phi=phi, mu=mu, s2r=s2r, states=states, A=A, it=it, ELBO=ELBO))
}
