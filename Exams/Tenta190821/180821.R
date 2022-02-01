nCovs = dim(X)[2]
mu_0 = rep(0,nCovs)
Omega_0 = (1/100)*diag(nCovs)
v_0 = 1
sigma2_0 = 5^2  

joint_posterior = BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter = 5000)
betas = joint_posterior$betaSample
sigmas = joint_posterior$sigma2Sample


Bmean = colMeans(joint_posterior$betaSample)
Bq025 = apply(joint_posterior$betaSample,2,quantile,.025)
Bq975 = apply(joint_posterior$betaSample,2,quantile,.975)
print(data.frame(round(cbind(Bmean,Bq025,Bq975),3)),row.names=covNames)#B compute density 
sigma_density = density(joint_posterior$sigma2Sample)
plot(sigma_density)
summary(sigma_density)


sorted_normalized_y = sort(sigma_density$y, decreasing = TRUE)/sum(sigma_density$y)
sorted_x = sigma_density$x[order(-sigma_density$y)]


count = 0
summa = 0
while(summa <= 0.95){
  count = count + 1
  summa = sorted_normalized_y[count] + summa
}
a = min(sorted_x[1:count-1])
b = max(sorted_x[1:count-1])

mode = sorted_x[which.max(sorted_normalized_y)]

XNewHouse
nSim <- dim(joint_posterior$betaSample)[1]

ypred = c()

for(i in 1:5000){
  ypred = append(ypred, XNewHouse%*%betas[i,]) + rnorm(1, 0, sqrt(sigmas[i]))
}

sum(ypred>20)/nSim

#2
alpha = 1
beta = 1
n = 5
sumx = 65
nsim = 5000

theta = rbeta(5000, n + alpha, sumx + beta)
predicted_observations = rgeom(5000, theta)

#Using quantile extraciting 0.95 we get the value for x where we with 95% probability can say that the next earthwake will have occured
quantile(predicted_observations, 0.95)


# 3
yData
xData
m = length(xData)
sumx = sum(xData)
log.posterior = function(xData, yData, theta){
  m = length(xData)
  sumx = sum(xData)
  log.theta_given_x = dgamma(theta, m + 3, sumx + 2, log = TRUE)
  #log.likelihood = sum(-3*log(1+(1/5)*(yData-log(theta))^2))
  log.likelihood = sum(dt(yData-log(theta), 5, log = TRUE))
  log.post = log.theta_given_x + log.likelihood
  return(exp(log.post))
}

posterior = c()
thetaGrid = seq(0,2, 0.01)
for(theta in thetaGrid){
  posterior = append(posterior, log.posterior(xData, yData, theta))
}


posterior_normalized = (1/0.01)*posterior/sum(posterior)
plot(thetaGrid, posterior_normalized)

logdNegBin <- function(x,mu,phi){
  lchoose(x+phi-1,x) + x*log(mu/(mu+phi)) + phi*log(phi/(mu+phi))
}

log.posterior = function(param, x){
  #initval[1] = my
  #initval[2] = phi
  #temp = factorial(x+phi-1)/(factorial(phi-1)*factorial(x))*(x*log(my/(my+phi))*phi*log(phi/(my+phi)))
  theta1 = param[1]
  theta2 = param[2]
  logPost = sum(logdNegBin(x, theta1, theta2))  -  2*log(theta2)
  return(logPost)
  #return(temp)
}


#

library(mvtnorm)
metropolis = function(n, c, initval, hessian, posterior_density, x){
  
  # this step depends on previous position. Previous position becomes this turns mean. 
  proposal_draws_previous = initval;
  
  acceptedDraws = matrix(0, ncol=2,nrow=n)
  accprobvec <- rep(0,n)
  
  set.seed(12345)
  for(i in 1:n){
    # draws (theta_p) from the proposal distribution ~ N(theta_p-1, c*hessian)
    proposal_draws = rmvnorm(1, proposal_draws_previous, c*hessian)
    proposal_draws[proposal_draws <= 0] = 1e-6
    # create a ratio depending on if this draw is better than previous, take exp to remove logarithm (logposterior)
    # posterior_density = log.posterior => exp of the division => logA -logB 
    acceptance_ratio = min(1,exp(posterior_density(proposal_draws, x)-posterior_density(proposal_draws_previous, x)))
    # draw a random uniformed variable to compare wiht acceptance ratio
    random_acceptance = runif(1,0,1)
    # if acceptance ratio is bigger than random variable than we move to the new position, otherwise we stay
    accprobvec[i] <- min(acceptance_ratio,1)
    if(acceptance_ratio >= random_acceptance){
      proposal_draws_previous = proposal_draws
      params = proposal_draws
    }
    acceptedDraws[i,] = params
    
  }
  return(list(draws = acceptedDraws, prob = accprobvec))
}

c = 0.1
initval =c(200,20)
hessian = diag(100,2)
hessian = c*diag(c(100,5))
x = incidents$incidents

testDraws = metropolis(10000, c=0.1, initval, hessian, log.posterior, x)

for(i in 1:2){
  plot(testDraws[,i], type='s')
  #a = c(rep(betas_posterior[i],nrow(testDraws)))
  #lines(a, col='red')
}


c = .8
MetropolisHastings <- function(c,niter,warmup,initVal,logPostFunc,x) {
  
  theta <- initVal
  thetamat <- matrix(0,length(theta),warmup+niter)
  thetamat[,1] <- theta
  accprobvec <- rep(0,warmup+niter)
  
  for(i in 2:(warmup+niter)) {
    thetaProp <- rmvgamma(theta,c)
    accprob <- exp(logPostFunc(thetaProp,x) - logPostFunc(theta,x) + ldmvgamma(theta,thetaProp,c))
    accprobvec[i] <- min(accprob,1)
    if(runif(1) < accprob) {
      theta <- thetaProp
    } 
    thetamat[,i] <- theta
  }
  
  return(list(thetamat=thetamat,accprobvec=accprobvec))
}

mp <- MetropolisHastings(c,niter,warmup,initVal,logPostNegBin,x)

# Optionally compute posterior mean and variance
# rowMeans(mp$thetamat[,warmup+1:niter])
# apply(mp$thetamat,1,var)

par(mfrow=c(1,2))
plot(mp$thetamat[1,warmup+1:niter], type="l",main="Traceplot",xlab="Iteration",ylab="mu")
plot(mp$thetamat[2,warmup+1:niter], type="l",main="Traceplot",xlab="Iteration",ylab="phi")
mean(mp$accprobvec)


