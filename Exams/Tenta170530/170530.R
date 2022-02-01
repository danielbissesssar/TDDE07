library(Bessel)

#compute log posterior 
# Since we don't have any information about the prior we assuem it to be constant. Ex 1 => log(constant) = constant
# when normalazing this all vlues will be affected equally hence it does not make any difference
log.post = function(theta, x, v){ 
  log_like = sum(log(x/v)*(-(x^2+theta^2)) + log(besselI(x*theta/v, 0))) + 2
return(log_like)
}

riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)
sequence = seq(0.1, 5, 0.1)
theta = 1
v = 1
prior = 1

# Plot the posterior over a grid of theta values. Gives us the probability depending on different theta values 
gridWidth <- 0.01
thetaGrid <- seq(0.01, 3, by = gridWidth)
logRicePostGrid <- rep(0,length(thetaGrid))
count <- 0
for (theta in thetaGrid){
  count <- count + 1
  logRicePostGrid[count] <-log.post(theta, riceData, v)
}

posterior = exp(logRicePostGrid)
posterior_normalized = (1/gridWidth)*posterior/(sum(posterior))

plot(posterior_normalized)


#b) Normal approximation of the posterior distr. of theta 

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- c(1); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode
OptParams<-optim(initVal,log.post,gr=NULL,riceData, v ,method=c("L-BFGS-B"),lower = 0, control=list(fnscale=-1),hessian=TRUE)

my_posterior = OptParams$par
# takes negative so that the posterior can be approx. as normal
# J = -second derivate evaluated in theta_hat
hessian_posterior = -OptParams$hessian

# take inverse for using it in the formula
hessian_posterior = (solve(hessian_posterior))
sigma = sqrt(hessian_posterior[1,1])

#Draw samples from Betas posterior distribution. 
set.seed(12345)

#We now have a posterior distributin which is ~N(my_posterior, sigma_square). We want to do the same thing as before
# plot the pdf over a grid of values => Since we have it in a good form we can do it immedeatly and directly get
# the normalized values 

approximate_density_distribution = dnorm(thetaGrid, mean = my_posterior, sd = sigma)
plot(thetaGrid, posterior_normalized)
lines(thetaGrid, approximate_density_distribution, col = "red", type ="b")


#C) 
#1. Explain on paper how the predictive distribution for a new observation is ocmputed by integration. Only
#The general dormula
#2. Compute the predictive distribution for a new observation by simulation. USe the approximate posterior from b
# simulator for r_rice is given by exam file

rRice <-function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}
#Draw from the posterior distribution and use these draws to plug in and make draw from test data
outcome =c()
set.seed(12345)
for (i in 1:1000){
  #nake draw from posterior distr.  
  theta_draw = rnorm(n = 1, mean = my_posterior, sd = sigma)
  #plug in draw and make fraw from predictions distr. 
  rice = rRice(1,theta_draw,1)
  outcome = append(outcome, rice)
}

#Histrogram of the predictive 
hist(outcome, main="Predictive distribution x_tilde", xlab="rice", ylab="Acumulation of rice out of 1000 draws")


#Assignemnt 2

# eBay bids data
load(file = 'bids.RData')    # Loading the vector 'bids' into workspace
bids = bids
bidsCounts <- table(bids)  # data2Counts is a frequency table of counts.
xGrid <- seq(min(bids),max(bids))  # A grid used as input to GibbsMixPois.R over which the mixture density is evaluated.

#a
#Copute the posterior distribution for theta and plot it => the posterior pdf
n = auctions = length(bids)
alpha = 1
beta = 1
amount_bids = sum(bids)

#Posterior parameters for gamma distriution
alpha_posterior = amount_bids + alpha
beta_posterior = n + beta

#Plot the posterior distribution for theta (gamma funktion) 
# => draws from gamma density distribution over a seequence of different values
gamma_grid = seq(3.4,3.9, length = 1000)
plot(gamma_grid, dgamma(gamma_grid, alpha_posterior, beta_posterior), xlab = "Theta", ylab = "Density")

#b
#Does the poisson model describe the distribution of the data well ?
#=> Model evaluation => Make draws from the posterior of theta. Use these to plug in
# to the distribution of the data (in this casepossion) and make new replica draws
data_density = density(bids)
plot(data_density)

xGrid = seq(min(bids),max(bids))

#Manually create the density of the data
#Normalize
data_dens = bidsCounts/sum(bidsCounts)
plot(xGrid, data_dens, type = 'l')

#Make repliica draws using draws from posterior
nDraws = 1000
set.seed(12345)
theta_draw = rgamma(nDraws, shape = alpha_posterior, rate = beta_posterior)
#For each draw theta we want to compute the possion distribution (distribution for
# x values)

density_mean = rep(0, length(xGrid))
for (i in 1:n){
#Compute poission distribution over Xgrid = över de x värden vi vill distribuera 
# proba
density_mean = density_mean + dpois(xGrid, theta_draw[i])
}
#Make a T function (condition) that we can compare with the original data. in this 
#case we compute the mean density for each grid value
density_mean = density_mean/n

plot(xGrid, density_mean, type = 'o', ylim = c(0,0.2), pch = 'o', xlab = "bids", ylab = "density")
lines(xGrid, data_dens, col = "red", type = 'b', pch = 'o')




#C

######################## Give Code######################
###############################
########## Problem 2 ########## 
############################### 

# eBay bids data
load(file = 'bids.RData')    # Loading the vector 'bids' into workspace
bidsCounts <- table(bids)  # data2Counts is a frequency table of counts.
xGrid <- seq(min(bids),max(bids))  # A grid used as input to GibbsMixPois.R over which the mixture density is evaluated.

# Code for Problem 3 - Exam in Bayesian Learning 2017-05-30
GibbsMixPois <- function(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter){
  
  # Gibbs sampling for a mixture of Poissons
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   x - vector with data observations (counts)
  #   nComp - Number of mixture components to be fitted
  #   alpha - The prior on the mixture component weights is w ~ Dir(alpha, alpha,..., alpha) 
  #   alphaGamma and betaGamma - 
  #              The prior on the mean (theta) of the Poisson mixture components is 
  #              theta ~ Gamma(alphaGamma, betaGamma) [rate parametrization of the Gamma dist]
  #   xGrid - the grid of data values over which the mixture is evaluated and plotted
  #   nIter - Number of Gibbs iterations
  #
  # OUTPUTS:
  #   results$wSample     - Gibbs sample of mixture component weights. nIter-by-nComp matrix
  #   results$thetaSample - Gibbs sample of mixture component means.   nIter-by-nComp matrix
  #   results$mixDensMean - Posterior mean of the estimated mixture density over xGrid.
  
  
  ####### Defining a function that simulates from a Dirichlet distribution
  rDirichlet <- function(param){
    nCat <- length(param)
    thetaDraws <- matrix(NA,nCat,1)
    for (j in 1:nCat){
      thetaDraws[j] <- rgamma(1,param[j],1)
    }
    thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
    return(thetaDraws)
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
  
  # Initial values for the Gibbs sampling
  nObs <- length(x)
  S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
  theta <- rep(mean(x), nComp) # Each component is initialized at the mean of the data
  
  # Setting up the grid where the mixture density is evaluated.
  mixDensMean <- rep(0,length(xGrid))
  effIterCount <- 0
  
  # Setting up matrices to store the draws
  wSample <- matrix(0, nIter, nComp)
  thetaSample <- matrix(0, nIter, nComp)
  probObsInComp <- rep(NA, nComp)
  
  # Setting up the priors - the same prior for all components
  alpha <- rep(alpha, nComp) 
  alphaGamma <- rep(alphaGamma, nComp) 
  betaGamma <- rep(betaGamma, nComp) 
  
  # HERE STARTS THE ACTUAL GIBBS SAMPLING
  
  for (k in 1:nIter){
    message(paste('Iteration number:',k))
    alloc <- S2alloc(S) # Function that converts between different representations of the group allocations
    nAlloc <- colSums(S)
    
    # Step 1 - Update components probabilities
    w <- rDirichlet(alpha + nAlloc)
    wSample[k,] <- w
    
    # Step 2 - Update theta's in Poisson components
    for (j in 1:nComp){
      theta[j] <- rgamma(1, shape = alphaGamma + sum(x[alloc == j]), rate = betaGamma + nAlloc[j])
    }
    thetaSample[k,] <- theta
    
    # Step 3 - Update allocation
    for (i in 1:nObs){
      for (j in 1:nComp){
        probObsInComp[j] <- w[j]*dpois(x[i], lambda = theta[j])
      }
      S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    
    # Computing the mixture density at the current parameters, and averaging that over draws.
    effIterCount <- effIterCount + 1
    mixDens <- rep(0,length(xGrid))
    for (j in 1:nComp){
      compDens <- dpois(xGrid, lambda = theta[j])
      mixDens <- mixDens + w[j]*compDens
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
  }
  return(results = list(wSample = wSample, thetaSample = thetaSample, mixDensMean = mixDensMean))
}
##################### Given COde#####################################


p = GibbsMixPois(bids, 2, alpha, alpha, beta, xGrid, 500)
p$mixDensMean
lines(xGrid, p$mixDensMean)

p = GibbsMixPois(bids, 3, alpha, alpha, beta, xGrid, 500)
p$mixDensMean
lines(xGrid, p$mixDensMean)




#Assignment 3
load("cars.RData")

install.packages("mvtnorm")
library(mvtnorm)

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

data = cars
data = as.matrix(cars)
Y = data[,1]
X_matrix = data[,-1]

my_0 = c(0,0,0,0)
v_0 = 4
omega_0 = diag(1,4)
omega_0 = omega_0*0.01
omega_invers = solve(omega_0)

sigma2_0 = 36
nIter = 1000

Joint_posterior = BayesLinReg(Y, X_matrix, my_0, omega_0, v_0, sigma2_0, nIter)
betas = Joint_posterior$betaSample
beta1 = betas[,1]
beta2 = betas[,2]
beta3 = betas[,3]
beta4 = betas[,4]


betas_density1 = density(beta1)
betas_density2 = density(beta2)
betas_density3 = density(beta3)
betas_density4 = density(beta4)

plot(betas_density1)
plot(betas_density2)
plot(betas_density3)
plot(betas_density4)

sigma = Joint_posterior$sigma2Sample
sigma_density = density(sigma)
plot(sigma_density)

#ii
#Compute a point estimate for each regression coefficient b0, b1, b2... assuming the linear loss fuunction L(Bk, alpha) = |Bk - alpha|
#Where Bk is the kth regression koefficient

beta1_estimate = median(beta1)
beta2_estimate = median(beta2)
beta3_estimate = median(beta3)
beta4_estimate = median(beta4)

quantile(Joint_posterior$sigma2Sample, c(0.025, 0.975))

#Construct a 95% equal tail probability interval for each parameter
quantile(beta1, probs = c(0.25,0.9725))
quantile(beta2, probs = c(0.25,0.9725))
quantile(beta3, probs = c(0.25,0.9725))
quantile(beta4, probs = c(0.25,0.9725))


#B
#Investigate if the effect on mpg is different in cars with 6 cylinders compared with 8 cylindes
# 6 => (1, 0)
# 8 => (0, 1)

difference = beta3-beta4
plot(density(difference))
hist(difference)
#Since the difference has a fairy high density aroun 0 it means that we cannot exclude the fact that they have the same distribution =>
# That they might have the same effect on response variable mpg

#C
#Compute by simulation the predictive dsitribution for a new car (1, 3.5, 0, 0)
Joint_posterior = BayesLinReg(Y, X_matrix, my_0, omega_0, v_0, sigma2_0, nIter)


# expression for the response variable y (logistic regression)
mpg_calc = function(prediction_data, betas, sigma){
  mpg = sum(prediction_data*betas) + rnorm(1, 0, sqrt(sigma))
  #mpg = betas[,1] + betas[,2]*weight + betas[,3]*sixcyl + betas[,4]*eightcyl
  return(mpg)
}

prediction_data = c(1, 3.5, 0, 0)
mpg = c()
set.seed(12345)
for (i in 1:nIter){
  Joint_posterior = BayesLinReg(Y, X_matrix, my_0, omega_0, v_0, sigma2_0, 1)
  betas = Joint_posterior$betaSample
  sigma = Joint_posterior$sigma2Sample
  mpg = append(mpg,mpg_calc(prediction_data, betas, sigma))
}

hist(mpg, 40, freq = FALSE)


#4
#let x1...xn be geomerically distributed (theta) for x >= 0 otherwise 0
#a - derive posterior distribution for thera using conjugate prior beta(alpha, beta)

#outcome = W L L W WW L L L W W 
#Geometric conversion of this data -> showing how many games it goes before we win is like
x = c(0,2,3,0)
alpha = 1
beta = 1
n = length(x)
xSum = sum(x)

#Predictive dsitribution wil be proportional to (according to calucaltion)
#gamma(k, + sum(x) + b)/gamma(k + alpha + n + sum(x) + 1 )

utility = 0
prob = c()
win = c()
for (i in 0:10){
prob = append(prob, (gamma(i + xSum + beta)/gamma(i + alpha + n + xSum + 1 )))
win = append(win, ((2^i -1) -2))
}
#Normalize so that sum of probabilities sum to 1
prob = prob/sum(prob)

Expected_utility = sum(prob*win)
Expected_utility





