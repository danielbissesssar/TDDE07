# Bayesian Learning Exam 2020-08-19
# Run this file once during the exam to get all the required data and functions for the exam in working memory 
# Author: Per Siden
set.seed(1)
###############################
########## Problem 1 ########## 
###############################
samples <- 10000
ndraw <- rnorm(samples,10000,500)
interval <- quantile(ndraw, probs=c(0.05,0.95))
print(interval)
#b)
patientdraws <- rpois(samples,ndraw)
patientinterval <- quantile(patientdraws, probs=c(0.05,0.95))
print(patientinterval)

#c)
seq <- seq(20,30,1)
lossfunction <- function(p,d){
  return((d + 0.04*pmax(p-500*d,0))**2)
}
losses <- 0*seq
count <- 1
for (i in seq) {
  losses[count] <- mean(lossfunction(patientdraws,i))
  count = count + 1
}
plot(seq,losses)
print(which.min(losses))
###############################
########## Problem 2 ########## 
############################### 

if("MASS" %in% rownames(installed.packages()) == FALSE) {install.packages("MASS")}
if("mvtnorm" %in% rownames(installed.packages()) == FALSE) {install.packages("mvtnorm")}
library(MASS)
library(mvtnorm)
muscle <- muscle

# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
  # Direct sampling from a Gaussian linear regression with conjugate prior:
  #
  # beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))
  # sigma2 ~ Inv-Chi2(v_0,sigma2_0)
  # 
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   Omega_0  - prior precision matrix for beta
  #   v_0      - degrees of freedom in the prior for sigma2
  #   sigma2_0 - location ("best guess") in the prior for sigma2
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   results$betaSample     - Posterior sample of beta.     nIter-by-nCovs matrix
  #   results$sigma2Sample   - Posterior sample of sigma2.   nIter-by-1 vector
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)
  
  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){
    
    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, df = v_n, scale = sigma2_n)
    sigma2Sample[i] = sigma2
    
    # Simulate from p(beta | sigma2, y, X)
    beta_ = rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_
    
  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}
y <- as.matrix(muscle[,3])
x <- as.matrix(muscle[,2])
x2 <- cbind(rep(1,60),x,x**2)
omega0 <- 0.001
v0 <- 1
sigma02 <- 10
iterations <- 5000
mu0m1 <- 0
mu0m2 <- c(0,0,0)
m1 <- BayesLinReg(y,x,mu0m1,omega0,v0,sigma02,iterations)
m2 <- BayesLinReg(y,x2,mu0m2,omega0*diag(3),v0,sigma02,iterations)
m1interval <- quantile(m1$betaSample, probs=c(0.025,0.975))
m2interval1 <- quantile(m2$betaSample[,1], probs=c(0.025,0.975))
m2interval2 <- quantile(m2$betaSample[,2], probs=c(0.025,0.975))
m2interval3 <- quantile(m2$betaSample[,3], probs=c(0.025,0.975))

#b)
seq <- seq(0,5,0.01)
m1beta <- mean(m1$betaSample)
count <- 1
lengthres <- rep(0,length(seq))
lengthbands <- matrix(0,length(seq),2)
for (i in seq) {
  lengthres[count] <- m1beta*i
  lengthbands[count,] <- quantile(lengthres[count], probs=c(0.025,0.975))
  count <- count + 1
}
plot(seq,lengthres,ylim = c(0,40))
points(x,y)
points(seq,lengthbands[,1])
points(seq,lengthbands[,2])


m2beta0 <- mean(m2$betaSample[,1])
m2beta1 <- mean(m2$betaSample[,2])
m2beta2 <- mean(m2$betaSample[,3])
count <- 1
lengthres <- rep(0,length(seq))
lengthbands <- matrix(0,length(seq),2)
for (i in seq) {
  lengthres[count] <- m2beta0+m2beta1*i+m2beta2*i**2
  lengthbands[count,] <- quantile(lengthres[count], probs=c(0.025,0.975))
  count <- count + 1
}
plot(seq,lengthres, ylim = c(0,40))
points(x,y)
points(seq,lengthbands[,1])

###############################
########## Problem 3 ########## 
############################### 

# Reading the data from file
load(file = 'C:/Users/Daniel Bissessar/Desktop/TDDE07/Exams/20200819/code/delivery.RData')
n <- length(delivery)
k <- c(0.5,1.5,2.5)
samp1 <- rgamma(10000,2+n,2+sum(delivery**k[1]))
invsamp1 <- 1/samp1
dens1 <- density(invsamp1)
pmode1 <- dens1$x[which.max(dens1$y)]
samp2 <- rgamma(10000,2+n,2+sum(delivery**k[2]))
invsamp2 <- 1/samp2
dens2 <- density(invsamp2)
pmode2 <- dens2$x[which.max(dens2$y)]
samp3 <- rgamma(10000,2+n,2+sum(delivery**k[3]))
invsamp3 <- 1/samp3
dens3 <- density(invsamp3)
pmode3 <- dens3$x[which.max(dens3$y)]

xgrid <- seq(0,10, 0.01)
wdens1 <- k[1]/pmode1*xgrid^(k[1]-1)*exp(-xgrid**k[1]/pmode1)
wdens2 <- k[2]/pmode2*xgrid^(k[2]-1)*exp(-xgrid**k[2]/pmode2)
wdens3 <- k[3]/pmode3*xgrid^(k[3]-1)*exp(-xgrid**k[3]/pmode3)


hist(delivery,20,  freq = FALSE)
lines(xgrid,wdens1, col = 2)
lines(xgrid,wdens2, col = 3)
lines(xgrid,wdens3, col = 4)

###############################
########## Problem 4 ########## 
###############################
logunorm <- function(x,k) {
  logprior <- dexp(k,1,log = TRUE)
  logposterior <- lgamma(102) - lgamma(2) + 2*log(2) + n*log(k) +
    sum((k-1)*log(x)) - (102)*log(2+sum(x^k))
  return(logprior+logposterior)
}
xgrid <- seq(0,3,0.01)
count <- 1
for (i in xgrid) {
sugkuk[count] <- logunorm(delivery,i)
count <- count + 1
}
sugkuknorm <- 1/0.01*exp(sugkuk)/sum(exp(sugkuk))
plot(xgrid,sugkuknorm, type = 'l')

kPostCDF <- cumsum(sugkuknorm)*0.01
lowerBound <- xgrid[which.min(abs(kPostCDF-0.025))]
upperBound <- xgrid[which.min(abs(kPostCDF-0.975))]
print(c(lowerBound,upperBound))
