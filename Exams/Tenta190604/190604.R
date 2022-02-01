
c1 = 4
c2 = 16

par(mfrow = c(1,2))
hist(rbeta(1000, sqrt(c1), 20), 30, main = "c = 4", freq = FALSE)
hist(rbeta(1000, sqrt(c2), 20), 30, main = "c = 16", freq = FALSE)

#lowe tail = false => we get the probability that theta (generated from ex rbeta will be larger than value x (0.1)) 
pbeta(0.1, sqrt(c1), 20, lower.tail = FALSE)
pbeta(0.1, sqrt(c2), 20, lower.tail = FALSE)

cGrid = seq(4,20,0.5)


#Using the utility function (quadratic loss) loss find the optimal value for c. For each value c we compute a vector of posssible theta values draws from rbeta
#We use these to compute a range of possible different utility values (for a given c). In order to choose 1 we have to make a point estimate of the utilty and since
#we are using quadratic loss => mean. OBS do not take mean of the theta values
utility_function = function(c){
  theta_draws = rbeta(10000, sqrt(c), 20)
  util = 0
  util = 100 + 20*log(theta_draws) - c
  return (mean(util))
}

par(mfrow=c(1,1))
utility_c1 =c()
for(c in cGrid){
  utility_c1 =  append(utility_c1, utility_function(c))
}

optC = cGrid[which.max(utility_c1)]

plot(cGrid, utility_c1, type = 'l', main = "utility curve for differenc values of c", xlab = "Cost c", ylab = "utlity", lwd = 2)
abline(v = optC, col = "blue", lwd = 2, lty = 2)
lines(cGrid[which.max(utility_c1)], utility_c1[which.max(utility_c1)], col = "red", type = 'p', pch = 21, cex = (1.2))
legend("topright", c("utility", "optimalc", "maxUtility"), col =c("black", "blue", "red"), cex = (1), lty = c(1,2), lwd = c(2, 2))



## 2


data = ebay
x = data
bidsCounts <- table(ebay)  # data2Counts is a frequency table of counts.
xGrid <- seq(min(ebay),max(ebay))  # A grid used as input to GibbsMixPois.R over which the mixture density is evaluated.
n = length(data)
N = 50

thetaGrid = seq(0,1,length = 1000)
Gridwidth = thetaGrid[2]-thetaGrid[1]


#log.likelihood2 = sum(sum(x)*log(0.5) + (n*N-sum(x))*log(1-0.5))
#log.likelihood3 = sum(x*log(0.5) + (N-x)*log(1-0.5)) # should be same as above
#log.likelihood = sum(dbinom(x, N, 0.5, log = TRUE)) #nromalized



log.posterior = function(theta, x){
  #log.likelihood = sum(x)*log(theta) + (n*N-sum(x))*log(1-theta)
  print(log.likelihood_p)
  #log.likelihood_test = sum(x*log(theta) + (N-x)*log(1-theta)) # should be same as above
  log.likelihood = sum(dbinom(x, N, theta, log = TRUE)) #nromalized
  log.prior = 2*log(1-theta)
  log.posterior = log.likelihood_p + log.prior
  return(exp(log.posterior))
}

test = sapply(thetaGrid, log.posterior, x = data)

posterior = c()
for(theta in thetaGrid){
  posterior = append(posterior, log.posterior(theta, x))
}

posterior_normalized = (1/Gridwidth)*posterior/sum(posterior)
test_normalized = (1/Gridwidth)*test/sum(test)

plot(thetaGrid, posterior_normalized)
plot(thetaGrid, test_normalized)
theta_mode = thetaGrid[which.max(posterior_normalized)]
max = max(posterior_normalized)
lines(theta_mode, max, col = "red", type = 'p')


#B

x

nComp = 2
nIter = 500
alpha = 1
betaGamma = 1
alphaGamma = 1

set.seed(100)
GibbsMixPois <- GibbsMixPois(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter)

thetas = GibbsMixPois$thetaSample

cum_mean1 = rep(0,nIter)
cum_mean2 = rep(0,nIter)
for (i in 1:nrow(thetas)){
  cum_mean1[i] = mean(thetas[,1][1:i])
  cum_mean2[i] = mean(thetas[,2][1:i])
}

plot(thetas[,1], ylim = c(0,10), type = 'l')
lines(thetas[,2])
lines(cum_mean1)
lines(cum_mean2)

#
par
hist(ebay)
hist(ebay_normalized)
plot(exp(GibbsMixPois$mixDensMean))
ebay_normalized = bidsCounts/sum(bidsCounts)
hist(ebay)
hist(ebay_normalized, 40, freq = FALSE)
plot(xGrid, ebay)
plot(dens)
lines(GibbsMixPois$mixDensMean)

plot(xGrid, ebay_normalized)
lines(dbinom(xGrid, 50, 0.103))



#
x = cellphones
n = length(x)
alpha1 = 2
beta1 = 1
alpha2 = 10
beta2 = 10
thetaGrid = seq(0,4,.01)
plot(thetaGrid,dgamma(thetaGrid,alpha1,beta1),type="l",ylab="",xlab="theta",main="Gamma prior densities")
lines(thetaGrid,dgamma(thetaGrid,alpha2,beta2),col="red")
legend(x=3,y=.3,c("M1","M2"),col = c("black","red"), lty = c(1,1))


### 4b
logprior <- function(theta,alpha,beta){
  return(dgamma(theta,shape=alpha,rate=beta,log=T))
}
loglik <- function(x,theta){
  return(sum(dexp(x,theta,log=T)))
}
logposterior <- function(theta,alpha,beta,x){
  return(dgamma(theta,shape=alpha+length(x),rate=beta+sum(x),log=T))
}
logmarglik <- function(theta,alpha,beta,x){
  return(loglik(x,theta) + logprior(theta,alpha,beta) - logposterior(theta,alpha,beta,x))
}

#
lms = c(logmarglik(0.5,alpha1,beta1,x),logmarglik(0.5,alpha2,beta2,x))
unnormProbs = .5*exp(lms)
probs = unnormProbs/sum(unnormProbs)
print(probs)

x = cellphones
n = length(x)
alpha1 = 2
beta1 = 1



p1 = ((beta1^alpha1)*gamma(n+alpha1))/(gamma(alpha1)*(beta1 + sum(x))^(alpha1+n))
p2 = ((beta2^alpha2)*gamma(n+alpha2))/(gamma(alpha2)*(beta2 + sum(x))^(alpha2+n))

#Multiplly my model prior probabilty = 0.5 for both
p1 = p1*0.5
p2 = p2*0.5
summa = p1 + p2

#normalize so they sum to 1
p1 = p1/summa
p2 = p2/summa

nSim = 100000
predictive_draw = c()

# OBS we get the model prior probability by rbinom. p2 is the lower probability so it will favor zeros. 
#therfor we compensate by + 1 to chose model 1. would also have worked with p = 1 which will favor
# 1. if M = 1 then we choose model 1. 
for (i in 1:nSim){
  M = rbinom(1,1, p2) + 1 # Simulate which model to use
  if(M==1){
  draws = rgamma(1, alpha1 + n, beta1 + sum(x))
} else{
  draws = rgamma(1, alpha2 + n, beta2 + sum(x))
}
  predictive_draw = append(predictive_draw, rexp(1, draws))
}

print(quantile(predictive_draw,probs = c(.05,.95)))


                                     