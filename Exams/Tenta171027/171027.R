  yProp

  thetaGrid = seq(0.01,15, length =1000)  
  
  lambda = 1

 posterior.calc = function(y, theta){
    
 #  for(i in length(y)){
  #   likelihood = likelihood + (theta-1)*log(y[i])+(theta-1)*log(1-y[i])
  # }
  # log.post = likelihood + log(lambda*exp(-lambda*theta))
   likelihood = sum(log(dgamma(y, theta, theta)))
   log.prior = dexp(theta, rate = 1, log = TRUE)
   log.post = likelihood + log.prior
   return(exp(log.post))
 }
  
  posterior = c()
  for(theta in thetaGrid){
    posterior = append(posterior, posterior.calc(yProp, theta))
  }


  gridwidth = thetaGrid[2] - thetaGrid[1]
  posterior_nromalized = (1/gridwidth)*posterior/sum(posterior)
  
  
  plot(thetaGrid, posterior_nromalized)
  
  #zero 1 loss => posterior mode
  posterior_mode = thetaGrid[which.max(posterior_nromalized)]
  
  
## b
  log.posterior = function(params, yVect){
    theta1 = params[1]
    theta2 = params[2]
    likelihood = sum(log(dgamma(y, theta1, theta2)))
    log.prior1 = dexp(theta1, rate = 1, log = TRUE)
    log.prior2 = dexp(theta2, rate = 1, log = TRUE)
    log.post = likelihood + log.prior1 + log.prior2
    return(log.post)
  }
  
  initVal <- as.vector(rep(0,2)); 
  OptParams<-optim(initVal,log.posterior,gr=NULL,yVect,method=c("L-BFGS-B"),lower = c(0.0001,0.0001), control=list(fnscale=-1),hessian=TRUE)
  my_posterior = OptParams$par
  hessian_posterior = -OptParams$hessian

  
  #2a)
  mu_0 = rep(0, ncol(X))
  Omega_0 = diag(1, ncol(X))
  Omega_0 = Omega_0*0.01
  sigma = 6^2
  v_0 = 1
  nSim = 5000
  
joint_posterior <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0 = sigma, nIter = nSim)
sort(joint_posterior$betaSample[,14])
betas = joint_posterior$betaSample
sigmas = joint_posterior$sigma2Sample

#ibrary(HDInterval)
hdiRange = hdi(joint_posterior$betaSample[,14], credMass=0.95)
quantile(joint_posterior$betaSample[,14], probs = seq(0,1, 0.025))
hdiRange

#2b

prediction_data_original = X[9,] 
prediction_data_new = X[9,]
prediction_data_new[14] =  prediction_data_new[14]*0.7


set.seed(12345)
mdv_pred_original = c()
mdv_pred_new = c()
for (i in 1:nrow(betas)){
  pred = sum(betas[i,]*prediction_data_original) + rnorm(1, 0, sqrt(sigmas[i])) 
  mdv_pred_original = append(mdv_pred_original, pred)
  pred = sum(betas[i,]*prediction_data_new) + rnorm(1, 0, sqrt(sigmas[i])) 
  mdv_pred_new = append(mdv_pred_new, pred)
}

difference = mdv_pred_new - mdv_pred_original
sum(difference>0)/nSim
hist(difference, 50)

diffPrice = joint_posterior$betaSample[,14]*(29.93*0.7-29.93)
hist(diffPrice,50)


# Hur verifiera att priset ökar ?????
# 95% HPD for the price difference
quantile(diffPrice,c(0.02,0.975))

# which does no contain zero, and the expected price increase is 
mean(diffPrice)

diffPrice = joint_posterior$betaSample[,14]*(29.93*0.7-29.93)
hist(diffPrice,50)


##4

Y = c(195, 191, 196, 197, 189)
Y_mean = mean(Y)
n = length(Y)
sigma2_posterior = (10^2)/n
#simulate 1000 draws from from the predicitive distribution
#Can calculate the parameters for the predicitive distributin and make draws immedeatly

set.seed(12345)
posterior_draws = rnorm(5000, Y_mean, sqrt(sigma2_posterior/n + sigma2_posterior))
hist(posterior_draws, 50)
sum(posterior_draws>230)

set.seed(12345)
posterior_draws2 = rnorm(5000, Y_mean, sqrt(sigma2_posterior/n))
weight_prediciton = c()
for(mean in posterior_draws2){
 weight_prediciton =  append(weight_prediciton, rnorm(1, mean = mean, sd = 10))
}
hist(weight_prediciton, 50)

#What is the prob that a weight in the next 365 days will have be larger than 230
max_weight = c()
for (i in 1:nSim){
prediction_365_days = max(rnorm(365, Y_mean, sqrt(sigma2_posterior/n + sigma2_posterior)))
max_weight = append(max_weight, prediction_365_days)

}
hist(max_weight, 50)

expected_a = Y_mean/10


predicted_y_wights = rnorm(365, Y_mean, sqrt(sigma2_posterior/n))
aGrid = seq(expected_a-50, expected_a + 50, length = 365)



ExpectedLoss<-function(a, maxWeightYear){
  ProbCollapse = sum(maxWeightYear>10*a)/nSim
  EL = (1-ProbCollapse)*(a) + ProbCollapse*(100+a)
  return(EL)
}


aGrid = seq(20,30,by = 0.01)
EL = rep(NA,length(aGrid),1)
count = 0
for (a in aGrid){
  count = count + 1
  EL[count] = ExpectedLoss(a, maxWeightYear)
}

plot(aGrid, EL, type = "l")
aOpt = aGrid[which.min(EL)] # This is the optimal a
points(x= aOpt, y = ExpectedLoss(a=aOpt, maxWeightYear), col = "red")


