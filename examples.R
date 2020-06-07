rm(list=ls())
library(HMMARM)
# inArgs <- commandArgs(trailingOnly = TRUE)
# missrate <- as.numeric(inArgs[1L])
# ID <- as.integer(inArgs[2L])

missrate <- 0
ID <-1

phi <- 0
mu <- 2
v <- 2

set.seed(ID)
l <- 1000
muset <- matrix(c(2, -1, 1, -1), 2, 2)
vset <- c(0.01, 0.1)
ll <- round(l / (1 - missrate))
mu1 <- muset[1, mu]
mu2 <- muset[2, mu]
sigma2 <- vset[v]

sim <- SimHomo(ll, mu1, mu2, sigma2, phi, missrate)
dat <- sim$x

# 1. ECM
ind <- which(!is.na(dat[-1]))
c1 <- sim$state.true[-1]
start.time <- Sys.time()
res <- try(ECM(dat, nphi=1, nvar=1, screen=FALSE, underflow='Scale'))
end.time <- Sys.time()
if(class(res) == "try-error") {
  # cat(c(cnt,rep(NA,11)), "\n",file=f,append = TRUE)
} else {
  c2 <- res$states %% 3
  c2[which(c2 == 0)] <- 3
  acc <- sum(c1[ind] == c2) / length(c1[ind])
  c3 <- c1
  c3[ind] <- c2
  c3[-ind] <- NA
  for (i in (length(c3) - 1):1) {
    if (is.na(c3[i])) {
      c3[i] <- c3[i + 1]}
  }
  acc2 <- sum(c1 == c3) / length(c3)
  cat(ID, 'ECM', missrate, mu1, mu2, res$parhat[2:3], acc, acc2,
      difftime(end.time, start.time, units="secs"), res$it, "\n")
}
cat(res$parhat)

# 2. ECM2
start.time <- Sys.time()
res <- try(ECM(dat, nphi=2, nvar=2, screen=FALSE, underflow='Scale'))
end.time <- Sys.time()
if(class(res) == "try-error") {
  # cat(c(cnt,rep(NA,11)), "\n",file=f,append = TRUE)
} else {
  c2 <- res$states %% 3
  c2[which(c2 == 0)] <- 3
  acc <- sum(c1[ind] == c2) / length(c1[ind])
  c3 <- c1
  c3[ind] <- c2
  c3[-ind] <- NA
  for (i in (length(c3) - 1):1){
    if (is.na(c3[i])) {
      c3[i] <- c3[i + 1]}
  }
  acc2 <- sum(c1 == c3) / length(c3)
  cat(ID, 'ECM2', missrate, mu1, mu2, res$parhat[3:4], acc, acc2,
      difftime(end.time, start.time, units ="secs"), res$it, "\n")
}
cat(res$parhat)


# 3. VI
theta.ture <- t(table(1:length(sim$state.true), as.factor(sim$state.true)))
Tll <- Truell(phi, dat, sigma2, c(0, mu1, mu2), theta.ture, sim$trans)
start.time <- Sys.time()
control <- list(epsilon=1e-3, maxit=200)
res <- try(VI(dat, screen=FALSE, control=control))
end.time <- Sys.time()
if(class(res) == "try-error") {
  # cat(ID,missrate, phi,mu1, mu2, sigma2,NA,NA,NA,NA,NA,NA,NA,NA,"\n")
} else {
  c1 <- res$states %% 3
  c1[which(c1 == 0)] <- 3
  ind <- which(!is.na(sim$x[-1]))
  acc <- sum(c1[ind] == sim$state.true[-1][ind]) / length(c1[ind])
  acc2 <- sum(c1 == sim$state.true[-1]) / length(c1)
  cat(ID, 'VI', missrate, mu1, mu2, res$mu[2:3], acc, acc2,
      difftime(end.time, start.time, units ="secs"), res$it, "\n")
}
cat(res$phi, res$s2r)
# 4. quartile threshold
start.time <- Sys.time()
lowerth <- summary(dat)[2]
higherth <- summary(dat)[5]
end.time <- Sys.time()
dat.nona <- dat[!is.na(dat)][-1]
c1 <- rep(1, length(dat.nona))
c1[dat.nona > higherth] <- 2
c1[dat.nona < lowerth] <- 3
acc <- sum(c1 == sim$state_true[-1][ind]) / length(c1[ind])
c3 <- sim$state_true[-1]
c3[ind] <- c1
c3[-ind] = NA
for (i in (length(c3) - 1):1) {
  if (is.na(c3[i])) {
    c3[i] <- c3[i + 1]}
}
acc2 <- sum(sim$state.true[-1] == c3) / length(c3)
cat(ID, 'threshold', missrate, mu1, mu2, higherth, lowerth, acc, acc2,
    difftime(end.time, start.time, units ="secs"), res$it, "\n")
