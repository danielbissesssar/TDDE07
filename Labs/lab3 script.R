#lab 3
#assignment 1
#a)

library(LaplacesDemon)
set.seed(12345)
data <- read.table("rainfall.dat")
lny <- log(data$V1)
n <- length(lny)
nDraws <- 10000
mu0 <- mean(lny)
xMean <- mean(lny)
taosqr0 <- 1
v0 <- 1
vN <- n+v0
sigmasqr0 <- 1

mu <- mu0
sigmasqr <- sigmasqr0

muVector <- rep(0,nDraws)
sigmaVector <- rep(0,nDraws)
for (i in 1:nDraws) {
  w <- (n/sigmasqr)/((n/sigmasqr)+(1/taosqr0))
  muN <- w*xMean + (1-w)*mu0
  taosqrN <- 1/((n/sigmasqr)+(1/taosqr0))
  
  mu <- rnorm(1,muN,taosqrN)
  muVector[i] <- mu
  
  
  sigma2N <- (v0*sigmasqr0 + sum((lny-mu)^2))/vN
  
  sigmasqr <- rinvchisq(1, vN, sigma2N)
  sigmaVector[i] <- sigmasqr
}

plot(1:nDraws, muVector, type = "l",col="red")
plot(1:nDraws, sigmaVector, type = "l",col="red")

mu_acf <- acf(muVector)
sigma_acf <- acf(sigmaVector)
barplot(height = mu_acf$acf[-1],col="blue")
barplot(height = sigma_acf$acf[-1],col="blue")
hist(muVector)
hist(sigmaVector)
mu_if <- 1+2*sum(mu_acf$acf)
sigma_if <- 1+2*sum(sigma_acf$acf)

#b)
hist(lny, freq = FALSE)
vec <- seq(0, 6, 0.01)
yTilde <- dnorm(vec, mu, sigmasqr)
lines(vec, yTilde)
#OK, not great

#assignment 2
#a)
library(glmnet)
OGdata <- read.table("eBayNumberOfBidderData.dat", header=TRUE)
data <- OGdata[,-2]
glmmodel <- glm(nBids~.,family = poisson,data=data)
abs_coeff <- abs(glmmodel$coefficients)
#VerifyID, Sealed, MinBidShare seem to have a larger impact and on the result and thus more significant

#b)
library(mvtnorm)
X <- as.matrix(OGdata[,-1])
Y <- as.matrix(OGdata[,1])
sd <- 100*solve((t(X)%*%X))
mu <- c(0,0,0,0,0,0,0,0,0)
B <- c(0,0,0,0,0,0,0,0,0)

logPost <- function(betas, y, x, mu, Sigma){
  pred <- x%*%betas
  lik <- sum(-log(factorial(y)) + y * pred - exp(pred))
  #lik <- sum(pred*y - log(1 + exp(pred)))
  prior <- dmvnorm(betas, mu, Sigma, log=TRUE)
  return(lik + prior)
}

res <- optim(B, logPost, y = Y, x = X, mu = mu, Sigma = sd, control=list(fnscale=-1), hessian = TRUE)
res$par
inv.hess <- -solve(res$hessian)
inv.hess

#c)
set.seed(12345)
av <- rep(1,5000)
RWMsampler <- function(logPostFunc, c, theta, BigSigma, iterations, ...){
  mat <- matrix(nrow=iterations, ncol=length(theta))
  mat[1,] <- theta
  
  for (i in 2:iterations){
    sample <- rmvnorm(1,mat[i-1,],c*BigSigma)
    prob_i <- logPostFunc(as.vector(sample), ...)
    prob_im1 <- logPostFunc(as.vector(mat[i-1,]), ...)
    q <- exp(prob_i-prob_im1)
    alpha <- min(1, q)
    av[i] <<- alpha
    if ( alpha >= runif(1) ){
      mat[i,] <- sample
    } else {
      mat[i,] <- mat[i-1,]
    }
  }
  return(mat)
}

rwm <- RWMsampler(logPost, 0.5, B, inv.hess, 5000, Y, X, mu, sd)
mean(av)
hist(rwm[,9])

#must show with graphs

#d)

newx <- c(1,1,1,1,0,1,0,1,0.7)

rwmB <- apply(rwm, 2, median)

draws <- rpois(1000, exp(t(newx)%*%rwmB))
hist(draws)

p <- sum(draws==0)/length(draws)
#p is the probability of having no bidders

# assignment 3
# a)

arp <- function(mu, phi, sigma2, t){
  x <- rep(0,t)
  e1 <- rnorm(1, mean=0, sigma2)
  x[1] <- mu + e1
  for (i in 2:t){
    et <- rnorm(1, mean=0, sigma2)
    x[i] <- mu + phi*(x[i-1]-mu) + et
  }
  return(x)
}
test <- arp(20, -0.1, 4, 200)
plot(test, type='l')

# b)
library(rstan)
library(BH)

phi9 <- arp(20, 0.9, 4, 200)
phi3 <- arp(20, 0.3, 4, 200)

StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  for (n in 2:N){
    y[n] ~ normal(mu + phi * y[n-1], sigma);
  }
}'

fit3 <- stan(model_code=StanModel, data=list(N=200,y=phi3), warmup=1000, iter=2000, chains=4)

postdraws3 <- extract(fit3)

fit9 <- stan(model_code=StanModel, data=list(N=200,y=phi9), warmup=1000, iter=2000, chains=4)

postdraws9 <- extract(fit9)

print(fit3, digits_summary = 3)
#for fit3 it was alright. Mu was a little bit low. With a Rhat lower than 1.05 it has converged, which is the case
print(fit9, digits_summary = 3)
#for fit9 sigma and phi was ok. Mu was very off. With a Rhat lower than 1.05 it has converged, which i sthe case

plot(postdraws3$mu, type='l')
#hard to tell
plot(postdraws9$mu, type='l')
#
plot(postdraws3$phi, type='l')
#
plot(postdraws9$phi, type='l')
#

plot(postdraws3$phi, postdraws3$mu)
plot(postdraws9$phi, postdraws9$mu)
