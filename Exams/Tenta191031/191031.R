
# Compute the bayecian decision 
p_up = 0.6
p_down = 1 - p_up
EU_buy = p_up*30 + p_down*(-10)
EU_notBuy = p_up*90 + p_down*(-120)

#EU_byy > nnot byt => should buy

p_posterior_up = 65/105
p_posterior_down = 1-p_posterior_up

EU_buy2 = p_posterior_up*30 + p_posterior_down*(-10)
EU_notBuy2 = p_posterior_up*90 + p_posterior_down*(-120)


#2
n = nrow(Traffic)
sumY = sum(Traffic$y)


thetaGrid = seq(18,24,.01)
alpha = 20
beta = 1
alpha_posterior = sumY + alpha
beta_posterior = n + beta
plot(thetaGrid, dgamma(thetaGrid, alpha_posterior, beta_posterior))

prob = pgamma(21, alpha_posterior, beta_posterior)

#b
data = Traffic
yes_subset <- subset(data, limit == 'yes')
no_subset = subset(data, limit == 'no')

draws_yes = rgamma(5000, sum(yes_subset$y) + alpha, nrow(yes_subset) + beta)
draws_no = rgamma(5000, sum(no_subset$y) + alpha, nrow(no_subset) + beta)
difference = draws_no - draws_yes 
hist(difference)

mean(draws_no>draws_yes)
#andelen av fall där no limit kommer att resultera att number off accidents (y ) är större än antalet accidents när speed limit fanns

#tar antalet accidents när inget limit fanns * 0.85 => ex 12 accidents * 0.85. Kollar nu andelen där de fortfarande skulle vara större n där det fanns limit
# Ser då att andelen minskar till 0.869. Det är är därmed rimligt att 
mean(draws_no*0.9>draws_yes)

#c
set.seed(12345)
x = 20
lambda = 30
alpha = 2
beta = 2
nDraws = 2000


v_draws = c(30)
pi_draws = c(0.5)
pi_draw = 1
for (i in 1:nDraws){
  #Compute Full conditional posterior (Normal model with conditionally conjugate prior)
  #my <- rnorm(1, mean = myn, sd = sqrt(taonSquared))
  z <- rpois(1, lambda*(1-pi_draw))
  v = z + x
  v_draws = append(v_draws, v)
  
  pi_draw = rbeta(1, alpha + x, beta + v -x)
  pi_draws = append(pi_draws, pi_draw)
}

v_draws = v_draws[500:nDraws]
pi_draws = pi_draws[500:nDraws]
par(mfrow=c(2,1))
#plot(v_draws, type = 'l')
#plot(pi_draws, type = 'l')

hist(pi_draws[500:nDraws], 30)
hist(v_draws[500:nDraws], 30)

set.seed(1235)
# start values
burnin = 500
niter = 2000
nu <- 30
pi <- .5
nu_vec <- rep(0,burnin+niter)
pi_vec <- rep(0,burnin+niter)
nu_vec[1] <- nu
pi_vec[1] <- pi
for(i in 2:(burnin+niter)){
  z <- rpois(1,lambda*(1-pi))
  nu = z + x
  nu_vec[i] <- nu
  pi <- rbeta(1,alpha+x,beta+nu-x)
  pi_vec[i] <- pi
}


par(mfrow=c(2,1))
hist(nu_vec[(burnin+1):(burnin+niter)],30,prob=TRUE,main="Posterior",xlab = "nu")
hist(pi_vec[(burnin+1):(burnin+niter)],30,prob=TRUE,main="Posterior",xlab = "pi")
par(mfrow=c(2,1))
plot(nu_vec[(burnin+1):(burnin+niter)],type="l",main="Traceplot",ylab = "nu",xlab="Iteration")
plot(pi_vec[(burnin+1):(burnin+niter)],type="l",main="Traceplot",ylab = "pi",xlab="Iteration")

#The Markov chain seems to have good mixing, since it rapidly explores the posterior, so the convergence is good.


#4

cars

library(rstan)

LinRegModel <- '
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma2;
}
model {
  sigma2 ~ scaled_inv_chi_square(5,10);
  for (n in 1:N)
    y[n] ~ normal(alpha + beta * x[n], sqrt(sigma2));
}'

DataRStan<-
  list(N = nrow(cars),
       x = cars$speed,
       y = cars$dist) 

fit_Model<-stan(model_code=LinRegModel,
                  data=DataRStan,
                  warmup=500,
                  iter=2000,
                  chains = 1)

print(fit_Model,digits=4)
res<-extract(fit_Model)
res


plot(DataRStan)
plot(cars)

xGrid = seq(0,25)
y_prediction_mean = c()
lower = c()
upper = c()
for (x in xGrid){
  ypred = res$alpha + res$beta*x  + rnorm(1500, 0, sqrt(res$sigma2))
  y_prediction_mean = append(y_prediction_mean, mean(ypred))
  qc = quantile(ypred, probs = seq(0,1, 0.05))
  lower = append(lower, qc[2])
  upper = append(upper, qc[20])
}

qc = quantile(ypred, probs = seq(0,1, 0.05))
qc[20]
lines(y_prediction_mean)
lines(upper)
lines(lower)

quantile(res$alpha, probs = c(0.05, 0.95))
res$sigma2
lines(res$sigma2)



#C

LinRegModel2 <- '
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real phi;
  real<lower=0> sigma2[N];
}
model {
  
  for (n in 1:N){
    sigma2[n] ~ scaled_inv_chi_square(5, exp(gamma + phi * x[n]));
    y[n] ~ normal(alpha + beta * x[n], sqrt(sigma2[n]));
  }
}'

fit_Model2<-stan(model_code=LinRegModel2,
                data=DataRStan,
                warmup=500,
                iter=2000,
                chains = 1)

print(fit_Model2,digits=4)
res<-extract(fit_Model2)
res$sigma2
res
res

plot(DataRStan)
plot(cars)

xGrid = seq(0,25)
y_prediction_mean = c()
lower = c()
upper = c()
for (x in xGrid){
  ypred = res$alpha + res$beta*x  + rnorm(1500, 0, sqrt(res$sigma2))
  y_prediction_mean = append(y_prediction_mean, mean(ypred))
  qc = quantile(ypred, probs = seq(0,1, 0.05))
  lower = append(lower, qc[2])
  upper = append(upper, qc[20])
}

qc = quantile(ypred, probs = seq(0,1, 0.05))
qc[20]
lines(y_prediction_mean)
lines(upper)
lines(lower)

quantile(res$alpha, probs = c(0.05, 0.95))
res$sigma2
lines(res$sigma2)
