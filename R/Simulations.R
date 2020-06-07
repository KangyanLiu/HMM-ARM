# Simulations by ARMHMM setting
# 1.Phi and sigma2 do not change
SimHomo <- function(l, mu1, mu2, sigma2, phi, missrate) {
  stay.prob <- 0.9  # probability of staying at current state
  trans <- matrix((1 - stay.prob) / 2, 3, 3) +
  		   diag(stay.prob - (1 - stay.prob) / 2, 3, 3)
  mu <- matrix(c(0, mu1, mu2), 1, 3)
  state <- diag(rep(1, 3))
  current.state <- 1
  x <- 0
  res <- c(x)
  state.true <- c(current.state)
  for (i in 2:l) {
    current.state <- t(rmultinom(1, 1, prob=trans[current.state, ])) %*%
                     matrix(c(1:3), 3, 1)
    state.true <- c(state.true, current.state)
    x <- phi * x + as.numeric(mu %*% state[current.state, ])+
         rnorm(1, 0, sqrt(sigma2))
    res <- c(res, x)
  }
  x_full <- res
  res[sample(1:length(res), length(res) * missrate)] <- NA
  while (is.na(res[1])) {
    res <- res[-1]
    x_full <- x_full[-1]
    state.true <- state.true[-1]
  }
  while(is.na(res[length(res)])) {
    res <- res[-length(res)]
    x_full <- x_full[-length(x_full)]
    state.true <- state.true[-length(state.true)]
  }
  result <- list(res, state.true, x_full, trans)
  names(result) <- c("x", "state.true", 'x_full', 'trans')
  return(result)
}

# 2. Phi and sigma2 change.
SimHetero <- function(l, mu1, mu2, sigma2, phi, missrate){
  stay.prob <- 0.9  # probability of staying at current state
  nvar <- length(sigma2)
  nphi <- length(phi)
  IJ <- 3 * nvar
  s <- 3 * nvar * nphi
  trans <- matrix((1 - stay.prob) / (s - 1), s, s) +
           diag(stay.prob - (1 - stay.prob) / (s - 1), s, s)
  comb <- t(expand.grid(1:3, 1:nvar, 1:nphi))
  mu <- matrix(c(0, mu1, mu2), 1, 3)
  state <- diag(rep(1, 3))
  current.state <- 1
  x <- 0
  res <- c(x)
  state.true <- c(current.state)
  for (i in 2:l){
    current.state <- t(rmultinom(1, 1, prob=trans[current.state, ])) %*%
    			     matrix(c(1:s), s, 1)
    state.true <- c(state.true, current.state)
    x <- phi[comb[2, current.state]] * x +
         as.numeric(mu[comb[1, current.state]]) +
         rnorm(1, 0, sqrt(sigma2[comb[3, current.state]]))
    res <- c(res, x)
  }
  x_full <- res
  res[sample(1:length(res), length(res) * missrate)] <- NA
  while (is.na(res[1])){
    res <- res[-1]
    x_full <- x_full[-1]
    state.true <- state.true[-1]
  }
  while(is.na(res[length(res)])){
    res <- res[-length(res)]
    x_full <- x_full[-length(x_full)]
    state.true <- state.true[-length(state.true)]
  }
  state.true <- state.true %% 3
  state.true[which(state.true == 0)] <- 3
  result <- list(res, state.true, x_full, trans)
  names(result) <- c("x", "state.true", 'x_full', 'trans')
  return(result)
}

# True likelihood
Truell <- function(phi, dat, s2r, mu, theta, trans) {
  ind <- which(!is.na(dat))
  y <- dat[ind[-1]]
  x <- dat[ind[-length(ind)]]
  D <- ind[-1]- ind[-length(ind)]
  CD <- c(0, cumsum(D))
  CD <- CD[-length(CD)]
  logA <- log(trans)
  pi <- theta[, 1]
  pi[which(pi == 0)] <- 4.940656e-324
  logpi <- matrix(log(pi), 3, 1)
  theta <- theta[, -1]
  ll <- -l1(phi, y=y, x=x, D=D, CD=CD, s2r=s2r, mu=mu, theta=theta) / 2 +
        l2(logpi, logA, matrix(unlist(theta), nrow=3)) -
        sum(theta * log(theta + (theta == 0) * 4.940656e-324))
}

# Simulations for AR without missing
# 3.Constant phi
SimHomoAR <- function(l, mu1, mu2, sigma2, phi){
  prob_vec <- rep(1 / 3, 3)   # Equal probability for each state
  mu <- matrix(c(0, mu1, mu2), 1, 3)
  state <- diag(rep(1, 3))
  current.state <- 1
  x <- 0
  res <- c(x)
  state.true <- c(current.state)
  for (i in 2:l) {
    current.state <- t(rmultinom(1, 1, prob=prob_vec)) %*% matrix(c(1:3), 3, 1)
    state.true <- c(state.true, current.state)
    x <- phi * x +
         as.numeric(mu %*% state[current.state, ]) +
         rnorm(1, 0, sqrt(sigma2))
    res <- c(res, x)
  }
  result <- list(res, state.true)
  names(result) <- c("x", "state.true")
  return(result)
}

# 4 Multiple phi and sigma2.
SimHeteroAR <- function(l, mu1, mu2, sigma2, phi){
  nvar <- length(sigma2)
  nphi <- length(phi)
  IJ <- 3 * nvar
  s <- 3 * nvar * nphi
  prob_vec <- rep(1 / s, s)
  comb <- t(expand.grid(1:3, 1:nvar, 1:nphi))
  mu <- matrix(c(0, mu1, mu2), 1, 3)
  state <- diag(rep(1, 3))
  current.state <- 1
  x <- 0
  res <- c(x)
  state.true <- c(current.state)
  for (i in 2:l) {
    current.state <- t(rmultinom(1, 1, prob=prob_vec)) %*% matrix(c(1:s), s, 1)
    state.true <- c(state.true, current.state)
    x <- phi[comb[2, current.state]] * x +
    	 as.numeric(mu[comb[1, current.state]]) +
         rnorm(1, 0, sqrt(sigma2[comb[3, current.state]]))
    res <- c(res, x)
  }
  state.true <- state.true %% 3
  state.true[which(state.true == 0)] <- 3
  result <- list(res, state.true)
  names(result) <- c("x", "state.true")
  return(result)
}


