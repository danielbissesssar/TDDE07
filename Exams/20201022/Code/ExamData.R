# Bayesian Learning Exam 2020-08-19
# Run this file once during the exam to get all the required data and functions for the exam in working memory 
# Author: Per Siden

###############################
########## Problem 1 ########## 
###############################
draw1 <- pbeta(0.9,96,6)
draw2 <- pbeta(0.9,88,14)
print(c(1-draw1,1-draw2))
#b)
n <- 10000
draw1 <- rbeta(n,96,6)
draw2 <- rbeta(n,88,14)
prob <- sum((draw1>draw2))/n
print(prob)
#c)
diff <- draw1-draw2
dens <- density(diff)
kPostCDF <- cumsum(dens$x)
lowerBound <- dens$x[which.min(abs(kPostCDF-0.025))]
upperBound <- dens$x[which.min(abs(kPostCDF-0.975))]
print(c(lowerBound,upperBound))

#d)
print(c(96/102,88/102))

###############################
########## Problem 2 ########## 
############################### 
n <- 3
x1 <- 9
x2 <- 1
x3 <- 9
alpha <- 1
beta <- 1
S = x1+x2+x3
mlike <- function(K){
  num <- choose(K,x1)*choose(K,x2)*choose(K,x3)*gamma(alpha+beta)*gamma(S+alpha)*gamma(n*K-S+beta)
  den <- gamma(alpha)*gamma(beta)*gamma(alpha+beta+n*K)
  return(num/den)
}
mlike1 <- mlike(10)
mlike2 <- mlike(20)

postprob <- c(mlike1,mlike2)/(mlike1+mlike2)
print(postprob)
###############################
########## Problem 3 ########## 
############################### 

# Reading the data from file
load(file = 'C:/Users/Daniel Bissessar/Desktop/TDDE07/Exams/20201022/Code/listeners.RData')

GibbsMixNormal <- function(x, nComp, nIter){
  
  # Gibbs sampling for a mixture of normals
  # Author: Mattias Villani and Per Siden, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   x - vector with data observations (counts)
  #   nComp - Number of mixture components to be fitted
  #   nIter - Number of Gibbs iterations
  #
  # OUTPUTS:
  #   results$muSample - Gibbs sample of mixture component means.   nIter-by-nComp matrix
  #   results$sigma2Sample - Gibbs sample of mixture component variances.   nIter-by-nComp matrix
  #   results$piSample     - Gibbs sample of mixture component weights. nIter-by-nComp matrix
  #   results$mixDensGrid - Grid of values over which the density can be plotted
  #   results$mixDensMean - Posterior mean of the estimated mixture density over xGrid.
  
  ###### Defining a function that simulates from the 
  rScaledInvChi2 <- function(n, df, scale){
    return((df*scale)/rchisq(n,df=df))
  }
  
  ####### Defining a function that simulates from a Dirichlet distribution
  rDirichlet <- function(param){
    nCat <- length(param)
    muDraws <- matrix(NA,nCat,1)
    for (j in 1:nCat){
      muDraws[j] <- rgamma(1,param[j],1)
    }
    muDraws = muDraws/sum(muDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
    return(muDraws)
  }
  
  # Simple function that converts between two different representations of the mixture allocation
  S2alloc <- function(S){
    n <- dim(S)[1]
    alloc <- rep(0,n)
    for (i in 1:n){
      alloc[i] <- which(S[i,] == 1)
    }
    return(alloc)
  }
  
  x <- as.matrix(x)
  
  # Prior parameters
  alpha <- 1*rep(1,nComp) # Dirichlet(alpha)
  muPrior <- rep(mean(x),nComp) # Prior mean of mu
  tau2Prior <- rep(10*sd(x),nComp) # Prior std mu
  sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
  nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2
  
  # Initial value for the MCMC
  nObs <- length(x)
  S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
  mu <- quantile(x, probs = seq(0,1,length = nComp))
  sigma2 <- rep(var(x),nComp)
  probObsInComp <- rep(NA, nComp)
  alloc <- S2alloc(S)
  nAlloc <- colSums(S)
  
  # Setting up the plot
  xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
  xGridMin <- min(xGrid)
  xGridMax <- max(xGrid)
  mixDensMean <- rep(0,length(xGrid))
  effIterCount <- 0
  
  # Save-matrices
  muSample = matrix(0,nIter,nComp)
  sigma2Sample = matrix(0,nIter,nComp)
  piSample = matrix(0,nIter,nComp)
  
  for (k in 1:nIter){
    # Update components probabilities
    pi_ <- rDirichlet(alpha + nAlloc)
    
    # Update mu's
    for (j in 1:nComp){
      precPrior <- 1/tau2Prior[j]
      precData <- nAlloc[j]/sigma2[j]
      precPost <- precPrior + precData
      piPrior <- precPrior/precPost
      muPost <- piPrior*muPrior + (1-piPrior)*mean(x[alloc == j])
      tau2Post <- 1/precPost
      mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
    }
    
    # Update sigma2's
    for (j in 1:nComp){
      sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
    }
    
    # Update allocation
    for (i in 1:nObs){
      for (j in 1:nComp){
        probObsInComp[j] <- pi_[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
      }
      S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    alloc <- S2alloc(S)
    nAlloc <- colSums(S)
    
    # Printing the fitted density against data histogram
    if (k%%1 == 0){
      effIterCount <- effIterCount + 1
      mixDens <- rep(0,length(xGrid))
      components <- c()
      for (j in 1:nComp){
        compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
        mixDens <- mixDens + pi_[j]*compDens
        components[j] <- paste("Component ",j)
      }
      mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    }
    
    
    # Save
    muSample[k,] = mu
    sigma2Sample[k,] = sigma2
    piSample[k,] = pi_
    
  }
  list(muSample=muSample,sigma2Sample=sigma2Sample,piSample=piSample,mixDensGrid=xGrid,mixDensMean=mixDensMean)
}
x <- listeners
nComp <- 2
nIter <- 1000
set.seed(100)
gibbs <- GibbsMixNormal(x, nComp, nIter)
par(mfrow=c(2,1))
plot(sqrt(gibbs$sigma2Sample[,1]),type = 'l')
plot(gibbs$muSample[,2], type = 'l')
piMean <- mean(gibbs$piSample[51:nIter])
mu1Mean <- mean(gibbs$muSample[51:nIter,1])
mu2Mean <- mean(gibbs$muSample[51:nIter,2])
sigma21Mean <- mean(gibbs$sigma2Sample[51:nIter,1])
sigma22Mean <- mean(gibbs$sigma2Sample[51:nIter,2])

hist(listeners, 20, freq = F)
lines(gibbs$mixDensGrid,gibbs$mixDensMean, col=2)
lines(gibbs$mixDensGrid, piMean*dnorm(gibbs$mixDensGrid, mu1Mean, sqrt(sigma21Mean)), col = 3)
lines(gibbs$mixDensGrid, (1-piMean)*dnorm(gibbs$mixDensGrid, mu2Mean, sqrt(sigma22Mean)), col = 4)

#3b
print(c(mu1Mean, mu2Mean, piMean))

###############################
########## Problem 4 ########## 
###############################

# Reading the data from file
load(file = 'zinc.RData')