#
# Reading the data vector yVect from file
load(file = 'CauchyData.RData')

dCauchy <- function(x, theta = 0, gamma = 1){
  return(dens = (1/(pi*gamma))*(1/(1+((x-theta)/gamma)^2)))
}

theta_grid = seq(3,8, length = 1000)

dlognormal <- function(x, mu, sigma2){
  return(dens = (1/(sqrt(2*pi*sigma2)*x))*exp((-1/(2*sigma2))*(log(x)-mu)^2))
}

log.posterior = function(yVect, theta){
  likelihood= sum(log(dCauchy(yVect, theta, 1)))
  log.prior = dnorm(theta, 0, 10, log = TRUE)
  log.post = likelihood + log.prior
  return(log.post)
}

posterior = c()
for (theta in theta_grid){
  log_posterior = log.posterior(yVect, theta)
  posterior = append(posterior, exp(log_posterior))
}

gridWidth = theta_grid[2] - theta_grid[1]
post_norm = (1/gridWidth)*posterior/sum(posterior)
plot(theta_grid, post_norm)




### b

dlognormal <- function(x, mu, sigma2){
  return(dens = (1/(sqrt(2*pi*sigma2)*x))*exp((-1/(2*sigma2))*(log(x)-mu)^2))
}

#Är det inte märkligt att köra X en dragning fråån log normal ? som i lab1 vad händer ?
log.posterior = function(params, yVect){
  theta = params[1]
  sigma = params[2]
  log.likelihood= sum(log(dCauchy(yVect, theta, sigma)))
  log.theta.prior = dnorm(theta, 0, 10, log = TRUE)
  log.sigma.prior = log(dlognormal(sigma, 0, 1))
  log.post = log.likelihood + log.theta.prior + log.sigma.prior
  return(log.post)
}


# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,2)); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode
OptParams<-optim(initVal,log.posterior,gr=NULL,yVect,method=c("L-BFGS-B"),lower = c(-Inf,0.0001), upper = c(Inf,Inf), control=list(fnscale=-1),hessian=TRUE)

betas_posterior = OptParams$par
# takes negative so that the posterior can be approx. as normal
# J = -second derivate evaluated in theta_hat
hessian_posterior = -OptParams$hessian

# take inverse for using it in the formula
hessian_posterior = solve(hessian_posterior)



## use the normal approximation above to obtain the marginal posterior for the 99% percentile of teh candy distribution 
# Theta + gamma*tan(pi(0.99 - 0.5))
library(mvtnorm)

#Fattar ej
draws = rmvnorm(5000, mean = betas_posterior, sigma =hessian_posterior )
cauchy = draws[,1] + draws[,2]*tan(pi*(0.99 - 0.5))

hist(cauchy, 50)



######## 2



# Reading the data from file
#library(MASS)
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)

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
my_0 = rep(0, 14)
Omega_0 = diag(0, 14)
Omega_0 = Omega_0*0.01
sigma = 6

joint_posterior = BayesLinReg(y, X, mu_0 = my_0,Omega_0 = Omega_0,  v_0 = 4, sigma2_0 = sigma^2 , nIter = 5000 )
betas = joint_posterior$betaSample
sigmas = joint_posterior$sigma2Sample

betas_estimate = c()
betas_quantile = matrix(0, ncol(betas), 2)
#Point estimator - Quadratiq loss => posterior mean
for (i in 1:ncol(betas)){
  betas_estimate = append(betas_estimate, mean(betas[,i]))
  quantiles_full = quantile(betas[,i], probs = seq(0,1, 0.025))
  betas_quantile[i,] = c(quantiles_full[2], quantiles_full[40])
}

sigma_quantile =  quantile(sigma, probs = seq(0,1, 0.025))
sigma_quantile = c(sigma_quantile[2], sigma_quantile[40])

#Nrooms is the 8th covariate
betas_quantile[7,]
#says that one unit increase of number of rooms will increase the median value with between 1.4611 and 5.76 X1000 $



# find the probability of a new observation for a new obervation holding X data
prediction_data = c(1.0000, 10, 0.0000, 18.1000, 0.0000, 0.6710, 6.9680, 91.9000, 1.4165, 24.0000, 666.0000, 20.2000, 396.9000, 17.2100 )

set.seed(12345)
mdv_pred = c()
for (i in 1:nrow(betas)){
pred = sum(betas[i,]*prediction_data) + rnorm(1, 0, sqrt(sigmas[i])) 
mdv_pred = append(mdv_pred, pred)
}

hist(mdv_pred)
mean(mdv_pred)

sum(mdv_pred>30)/length(mdv_pred)
#mdv_pred = ifelse(mdv_pred >= 30, 1, 0)
probability = sum(mdv_pred)/length(mdv_pred)
probability

##4
#let x1..x2 | theta ~ exp(theta)
nSim = 1000
Y =c(220, 323, 174, 229)
n = 4
sumY = sum(Y)
alpha = 25
Beta = 0.1
alpha_posterior = sumY + alpha
beta_posterior = n + Beta
gamma_draws = rgamma(nSim, sumY + alpha, n + Beta)

poisson_draws = c()
for (i in 1:nSim){
poisson_draws = append(poisson_draws, rpois(1, gamma_draws[i]))

}
hist(poisson_draws, 50, freq = FALSE)
sum(poisson_draws<=200)/nSim

#when n goes toward large the expected value will go towards the mean
expected_X = round(mean(poisson_draws))

aGrid = seq(expected_X-100, expected_X+100)
utility_vec = c()
for (i in aGrid){
  utility_vec = append(utility_vec, mean(utility(i, poisson_draws)))
}
plot(aGrid,utility_vec)

aGrid[which.max(utility_vec)]


