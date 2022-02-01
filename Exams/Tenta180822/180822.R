x = lions
sigma2 = 0.04

log.posterior = function(my, x, sigma2){
  likelihood = sum(dlnorm(x, my, sqrt(sigma2), log = TRUE))
  #likelihood = sum(log(1/(sqrt(2*pi*sigma)*x))-1/(2*sigma2)*(log(x)-my)^2)
  log.prior = dnorm(my, 5, 1, log = TRUE)
  log.posterior = likelihood + log.prior
  return(exp(log.posterior))
}

myGrid = seq(5.15,5.4,length = 1000)

posterior = sapply(myGrid, log.posterior, x = lions, sigma2 = sigma2 )
posterior = c()
for(my in myGrid){
  posterior = append(posterior, log.posterior(my, lions, sigma2))
}

Gridwidth = myGrid[2] -myGrid[1]
posterior_normalized = (1/Gridwidth)*posterior/sum(posterior)
plot(myGrid, posterior_normalized, type = 'l')


hist(x)

#assume now 

#library(rstan)

StanModel= '
data {
  int<lower=0> N;
  real x[N];
}
parameters {
  real mu;
  real<lower=0> sigma2;
}
model {
  mu ~ normal(5, sqrt(0.04));
  sigma2 ~ inv_gamma(5,1);
  for (n in 1:N){
  x[n] ~ lognormal(mu, sqrt(sigma2));
  }
}
'


DataRStan<-
  list(N = length(x),
       x = x) 

Stan_Model<-stan(model_code=StanModel,
                  data=DataRStan,
                  warmup=500,
                  iter=2000,
                  chains=1)

print(Stan_Model,digits=4, pars =c("mu", "sigma2"), probs = c(0.025, 0.975))
pairs(fit, pars = c("mu", "sigma2"))
res<-extract(Stan_Model)
res

hist(res$mu, freq = FALSE)
plot(res$mu, res$sigma2)
hist(res$sigma, freq = FALSE)


#c
#give an average weigth of of male lions. give an 95% credible interval of the average weight of the male
#lions based on the posterior calculated in previous step
# use formula for expected value for log normal

mean_weights = exp(res$mu + (1/2)*res$sigma2)
mean(mean_weights)
quantile(mean, probs = c(0.025, 0.975))


#2
titanic
covNames <- names(titanic)[2:length(names(titanic))]
y = titanic$survived
X = titanic[,2:6] # Adding a column of ones for the intercept
y <- as.numeric(y)
X <- as.matrix(X)

T2 = 50
n = length(y) # Number of observations
nCovs = dim(X)[2] # Number of covariates
my = rep(0, nCovs)
sigma_prior = diag(T2, nrow = nCovs)
library(mvtnorm)

log.posterior = function(betas, X,Y,my,sigma_prior){
  
  # is simply log of the product of density function
  log_likelihood = sum((X%*%betas)*Y-log(1+exp(X%*%betas)))
  
  #log of deensity for nultivariate normal => dmvnorm(density multivariate)
  log_prior = dmvnorm(x=betas,mean=my,sigma=sqrt(sigma_prior), log=TRUE)
  
  return(log_likelihood + log_prior)
}


# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,nCovs)); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode
OptimResults<-optim(initVal,log.posterior,gr=NULL,X,y,my,sigma_prior,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)


#Plot the marginal (approximated) posteriors of beta. Since they are all normal they can be written separately


postMode <- OptimResults$par
postCov <- -solve(OptimResults$hessian) # Posterior covariance matrix is -inv(Hessian)
names(postMode) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(postCov)) # Computing approximate standard deviations.
names(approxPostStd) <- covNames # Naming the coefficient by covariates

par(mfrow=c(2,2))
for(i in 2:5){
  grid = seq(postMode[i] - 3*approxPostStd[i],postMode[i] + 3*approxPostStd[i],length=1000)
  plot(grid,dnorm(grid,postMode[i],approxPostStd[i]),type="l",main=c(names(postMode)[i],"posterior"),ylab="density")
}

man_data = c(1,1,1,0,0)
woman_data = c(1,1,0,1,0)

betaDraws = rmvnorm(1000, mean = betas_posterior, sigma = postCov)

#compute the probabilty 
#Function for computing p(y = 1)
logistic_regression = function(data, betaDraws){
  exp(betaDraws%*%data) / (1+exp((betaDraws%*%data)))
}

prob_man = logistic_regression(man_data, betaDraws)
outcome_man = rbinom(1000,1, prob = prob_man)

prob_woman = logistic_regression(woman_data, betaDraws)
outcome_woman = rbinom(1000,1, prob = prob_woman)

outcome = cbind(outcome_man, outcome_woman)
sum(outcome[,1]==0 & outcome[,2] ==1)/1000

#OR 
ySim <- as.numeric(outcome_woman & !outcome_man)
mean(ySim)



#3c
alpha = 0.5
beta = 0.5
sumx = 15
n = 3

log.marginalLikelihoodP1 = lbeta(alpha + n, sumx+ beta) - lbeta(alpha,beta)
theta = 0.5
log.marginalLikelihoodP2 = sumx*log(1-theta) + n*log(theta)


unnorm_modelprob = c(exp(log.marginalLikelihoodP1)/10,exp(log.marginalLikelihoodP2)*(9/10))
modelprob = unnorm_modelprob / sum(unnorm_modelprob)
modelprob


#4a
# 4b
# Setting up data and prior
y <- c(184,67,149)
alpha <- c(1,1,1) # Dirichlet prior hyperparameters
nIter <- 10000 # Number of posterior draws
# Defining a function that simulates from a Dirichlet distribution
SimDirichlet <- function(nIter, param){
  nCat <- length(param)
  thetaDraws <- as.data.frame(matrix(NA, nIter, nCat)) # Storage.
  for (j in 1:nCat){
    thetaDraws[,j] <- rgamma(nIter,param[j],1)
  }
  for (i in 1:nIter){
    thetaDraws[i,] = thetaDraws[i,]/sum(thetaDraws[i,])
  }
  return(thetaDraws)
}
# Posterior sampling from Dirichlet posterior

#get theta draws from its posterior distribution
theta_draws = SimDirichlet(nIter, alpha + y)

#B
#Compute the probability that party A gets a majority of the votes (more than 50%)
#in the election. Assume that everyone in the population is voting.
sum(theta_draws[,1]>0.5)/nIter

#C
# Compute the probability that party A becomes the largest party. Assume that
#everyone in the population is voting.
sum(theta_draws[,1] > theta_draws[,2] & theta_draws[,1]> theta_draws[,3])/nIter

mean(theta_draws[,1])


#Calculate the expected value and standards deviation
# p = 0.971
p = 0.9
stdev = sqrt(p*(1-p)/10000)
c(p-1.96*stdev,p+1.96*stdev)


#e how many draws to make the previous confindence band in half

# according to monte carlo estimate the variance of an estimate will be sigma^2/N where sigma^2 is the varance for each obsevation
"The width from previous sample can be calculated with 2*sqrt(sigma^2/1000)
Solve the equation sqrt(sigma^2/1000) = 2sqrt(sigma^2/N)