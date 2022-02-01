# Bayesian Learning Exam 2020-06-04
# Run this file once during the exam to get all the required data and functions for the exam in working memory 
# Author: Per Siden

###############################
########## Problem 1 ########## 
############################### 
#b)
x <- 33
gridstep <- 0.001
thetaGrid <- seq(0,1,gridstep)
unnormpost <- dbinom(x,50,thetaGrid)*(thetaGrid>0.3)*(thetaGrid<0.7)
postB <- 1/gridstep*unnormpost/sum(unnormpost)
plot(thetaGrid,postB, type = 'l')
lines(thetaGrid, dbeta(thetaGrid,x+1,51-x), col = 2)
#c)
probA <- pbeta(0.5,x+1,51-x)
probB <- sum(postB[thetaGrid<=0.5]*gridstep)
print(c(probA, probB))
###############################
########## Problem 2 ########## 
############################### 
x <- c(log(20),log(20),log(50),log(40))
y <- c(5,3,17,8)
logPost <- function(x, y, beta) {
  linpred <- x*beta
  logLik <- sum(dpois(y, exp(linpred), log = TRUE))
  logPrior <-dnorm(beta, 1, 0.1, log = TRUE)
  return(logLik+logPrior)
}

optim <- optim(1, logPost, x = x, y = y, control = list(fnscale = -1), hessian = TRUE, method = 'BFGS')
mean <- optim$par
sd <- sqrt(-solve(optim$hessian))

#b)
n <- 10000
x5 <- log(20)
betasim <- rnorm(n,mean,sd)
y <- rpois(n,exp(x5*betasim))
l1 <- 4 + exp(x5)/50 -sqrt(y)
l1 <- mean(l1)

x5 <- log(40)
y <- rpois(n,exp(x5*betasim))
l2 <- 4 + exp(x5)/50 -sqrt(y)
l2 <- mean(l2)
###############################
########## Problem 3 ########## 
############################### 

# Reading the data from file
load(file = 'C:/Users/Daniel Bissessar/Desktop/TDDE07/Exams/20200604/Code/titanic.RData')

if("mvtnorm" %in% rownames(installed.packages()) == FALSE) {install.packages("mvtnorm")}
if("msm" %in% rownames(installed.packages()) == FALSE) {install.packages("msm")}
library(mvtnorm) # For mulitvariate normal
library(msm) # For truncated normal

BayesProbReg <- function(y, X, mu_0, tau, nIter){
  # Gibbs sampling in probit regression using data augmentation:
  #
  # beta | tau ~ N(mu_0, tau^2*I)
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   tau - prior standard deviation for beta
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   betaSample    - Posterior samples of beta.          nIter-by-nCovs matrix
  
  # Prior
  nPara <- dim(X)[2] # this line was missing in the exam
  priorCov <- tau^2*diag(nPara)
  priorPrec <- solve(priorCov)
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  n1 = sum(y)
  n0 = n - n1
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  
  # The actual sampling
  betaSample = matrix(NA, nIter, nCovs)
  u <- matrix(NA, n, 1)
  beta <- solve(XX,crossprod(X,y)) # OLS estimate as initial value
  for (i in 1:nIter){
    
    xBeta <- X%*%beta
    
    # Draw u | beta
    u[y == 0] <- rtnorm(n = n0, mean = xBeta[y==0], sd = 1, lower = -Inf, upper = 0)
    u[y == 1] <- rtnorm(n = n1, mean = xBeta[y==1], sd = 1, lower = 0, upper = Inf)
    
    # Draw beta | u
    betaHat <- solve(XX,t(X)%*%u)
    postPrec <- XX + priorPrec
    postCov <- solve(postPrec)
    betaMean <- solve(postPrec,XX%*%betaHat + priorPrec%*%mu_0)
    beta <- t(rmvnorm(n = 1, mean = betaMean, sigma = postCov))
    betaSample[i,] <- t(beta)
    
  }
  
  return(betaSample=betaSample)
}
iterations <- 1000
y <- as.matrix(titanic[,1])
x <- as.matrix(titanic[,-1])
tao <- 50
mu_0 <- c(mean(x[1,]),mean(x[,2]),mean(x[,3]),mean(x[,4]),mean(x[,5]))
sampling <- BayesProbReg(y,x,mu_0,tao,iterations)

hist(sampling[,1],50,freq = FALSE)
hist(sampling[,2],50,freq = FALSE)
hist(sampling[,3],50,freq = FALSE)
hist(sampling[,4],50,freq = FALSE)
hist(sampling[,5],50,freq = FALSE)

#b)

#c)
prob <- mean(sampling[,2]+sampling[,5]>0)
###############################
########## Problem 4 ########## 
###############################
xbar_FL <- 14
xbar_FW <- 300
xbar_ML <- 12
xbar_MW <- 280
sigma_L <- 2
sigma_W <- 50

prob_unnorm_F <- dnorm(10,14,sigma_L*sqrt(1+1/16))*dnorm(250,300,sigma_W*sqrt(1+1/16))*0.75
prob_unnorm_M <- dnorm(10,12,sigma_L*sqrt(1+1/4))*dnorm(250,280,sigma_W*sqrt(1+1/4))*0.25
prob_F <- prob_unnorm_F/(prob_unnorm_F+prob_unnorm_M)
print(prob_F)
