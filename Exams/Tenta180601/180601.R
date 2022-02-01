## 1

#Plot the prior and posterior distribution 
# Plotted its true distirbution over a gammaGrid

alpha = 20
beta = 2
sumx = 500
n = 50

plot(GammaGrid, dgamma(GammaGrid, 20, 2))
lines(GammaGrid, dgamma(GammaGrid, 10*50 +20, 50 +2 ))

priordraws = rgamma(1000, alpha, beta)
posteriordraws = rgamma(1000, sumx + alpha, n + beta)

#Plot the prior and posterior dsitribution using the samples
hist(priordraws, freq = FALSE)
hist(posteriordraws, freq = FALSE)

#Plot the prior and posterior using their mathematical expression 
#(The true distribution)
#=> Do it over a grid and plot the sensity on this grid
GammaGrid = seq(0, 20, length = 1000)
lines(GammaGrid, dgamma(GammaGrid, alpha, beta))
lines(GammaGrid, dgamma(GammaGrid, sumx + alpha, n + beta))

#B simulate 1000 draws from the predictive distribution of a new observation
#and plot the distribution using the samples 
predictiveDraws = rpois(1000, posteriordraws)
hist(predictiveDraws, freq = FALSE, 50)

#Get the probabilty for that x == 10 by simulation
sum(predictiveDraws==10)/1000



#2
fish
library(mvtnorm)
mu = c(0, 0, 0)
omega_0 = 0.01*diag(1,3)
v_o = 1
sigma2_0 = 10000
nIter = 5000


data = fish
data = as.matrix(fish)
Y = data[,1]
X_matrix = data[,-1]


joint_posterior = BayesLinReg(Y, X_matrix, mu, omega_0, v_0, sigma2_0, nIter)
betas = joint_posterior$betaSample
sigmas = joint_posterior$sigma2Sample

hist(sigmas)
for (i in 1:ncol(betas)){
  hist(betas[,i])
}

cat("equal tail interval for B1:", quantile(betas[,2], probs = c(0.05, 0.95)))
# 2.285509 2.964409
# an increase in the age (amount of days) of the fish will with 95% probability increase the length of the fish with between
#2.66 and 2.96 mm

#c
# Evidence says that b2 might not be useful. How to handle this ? 
#=> similar like ridge regression but more simple. Adjust the prior for Beta2|sigma
# atm it is ~N(0, sigma^2/0.01) ( precision matrix ) => it will have a big varance resulting in a difference in the model
# by changing the precision matrix to be (3,3) = 1000 it will minimizethe variance applied for beta2 and since it is centered around 0 it will not have any impact


#d)
prediction_data_30 = c(1, 30, 30)
prediction_data_100 = c(1, 100, 30)

ytilde <- rep(0,nIter)
xtilde <- c(1,3.5,0,0)
for (i in 1:nIter){
  if (runif(1) > .5) {
    xtilde = c(1,30,30)
  } else {
    xtilde = c(1,100,30)
  }
  ytilde[i] = sum(xtilde*Results$betaSample[i,]) + rnorm(n=1, mean = 0, sd = sqrt(Results$sigma2Sample[i]))
}
hist(ytilde,50, main ="Predictive distribution of y", freq = FALSE)

## 3

#M1
n = 10
x = 3
a = 1 
b = 1
p_m1 = p_m2 = p_m3 = 1/3

log.marginalLikelihoodP1 = lchoose(n,x) + lbeta(x+a,n-x+b) - lbeta(a,b)

a = 4
b = 4
log.marginalLikelihoodP2 = lchoose(n,x) + lbeta(x+a,n-x+b) - lbeta(a,b)
p = 0.5
log.marginalLikelihoodP3 = lchoose(n,x) + x*log(0.5) + (n-x)*log(0.5)

unnorm_modelprob = c(exp(log.marginalLikelihoodP1)/3,exp(log.marginalLikelihoodP2)/3,exp(log.marginalLikelihoodP3)/3)
modelprob = unnorm_modelprob / sum(unnorm_modelprob)
modelprob



#4

data = sulfur
data = subset(data, data >200)
x = sulfur[sulfur > 200]
myGrid = seq(100,400,1)
sigma =100
L = 200

log.posterior = function(my, x, L, sigma){
  sum(dnorm((x-my)/sigma, log = TRUE) - log(sigma) - log(1 - pnorm((L-my)/sigma)))
  #dnorm((x-mu)/sigma,log=TRUE) - log(sigma) - log(1-pnorm((L-mu)/sigma))
}

posterior = c()
for(my in myGrid){
posterior = append(posterior, log.posterior(my, data, L, sigma))  
}

posterior_normalized = posterior/sum(posterior)
hist(posterior_normalized)
plot(myGrid, posterior_normalized)



#4b

x = sulfur
#library(rstan)
T = length(sulfur)
T_cens = sum(sulfur <= 200)

censData <- list(T=T, T_cens = T_cens, x=sulfur, L=200)

# Model
censModel <- '
data {
  int<lower=0> T;       // Total number of time points
  int<lower=0> T_cens;  // Number of censored time points
  real x[T];            // Partly censored data
  real<upper=max(x)> L; // Lower truncation point
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<upper=L> x_cens[T_cens]; // Censored values
}
model {
  int t_cens = 0;
  for (t in 1:T){
    if (x[t] > L) 
      x[t] ~ normal(mu,sigma);
    else {
      t_cens += 1;
      x_cens[t_cens] ~ normal(mu,sigma);
    }
  }
}'



fit_cens =stan(model_code=censModel, 
               data=censData,
               warmup  = 500,
               iter = 2000,
               chains = 1)


print(fit_cens,digits=4)
res<-extract(fit_cens)

plot(res$sigma)
res$x_cens

hist(res$sigma, freq = FALSE)
hist(res$mu, freq = FALSE)



normal_Model = '

data {

  int<lower=0> T;       // Total number of time points
  int<lower=0> T_cens;  // Number of censored time points
  real x[T];            // Partly censored data
  real<upper=max(x)> L; // Lower truncation point
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<upper=L> x_cens[T_cens]; // Censored values
  real z[T];
  real phi;
}

model {
  int t_cens = 0;
  phi ~ uniform(-1,1);
  mu ~ normal(300, 10);
  for (t in 2:T){
    z[t] ~ normal(mu + phi * (z[t-1] - mu), sigma);
  }
  for (t in 1:T){
     if (x[t] > L) 
      x[t] ~ normal(z[t], 20);
    else {
      t_cens += 1;
      x_cens[t_cens] ~ normal(mu,sigma);
    }
    
  }

}'


fit_cens =stan(model_code=normal_Model, 
               data=censData,
               warmup  = 500,
               iter = 2000,
               chains = 1)

res = extract(fit_cens)
hist(res$phi)
z = res$z
plot(sulfur, type ="p", ylim = c(0, 500), xlab = "time")


postThetaMean = colMeans(res$z)
postThetaBands = apply(res$z, 2, quantile, probs = c(0.025, 0.975))

lines(postThetaMean,type="l",col=2)
lines(postThetaBands[1,],type="l",col=c(3))
lines(postThetaBands[2,],type="l",col=c(3))
